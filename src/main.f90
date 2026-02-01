program main
  use libMath, only: isInverse, cross_product, matmulAX, zAxis
  use libCommon
  use libPostprocess
  implicit none

  ! Ensure all necessary files exist
  inquire(file='config.nml', exist=fileExists)
  if (.not. fileExists) error stop 'ERROR: config.nml does not exist'

  ! Read config.in file
  print*
  call print_status('Reading file '//'config.nml')
  call readConfig('config.nml', ResultsDir//'config.nml')
  call print_status()    ! SUCCESS

  ! Allocate rotor objects
  allocate (rotor(nr))

  ! Read rotor??.in files
  do ir = 1, nr
    write (rotorChar, '(I0.2)') ir
    rotorFile = 'geom'//rotorChar//'.nml'
    inquire(file=rotorFile, exist=fileExists)
    call print_status('Reading file '//rotorFile)
    if (.not. fileExists) error stop 'ERROR: A geomXX.nml file does not exist'
    call rotor(ir)%readGeom(rotorFile, ResultsDir//rotorFile)
    call print_status()    ! SUCCESS
  enddo

  ! Rotor and wake initialization
  do ir = 1, nr
    if (rotor(ir)%surfaceType .ge. 0) then
      call rotor(ir)%init(ir, density, dt, nt, switches)
    else
      call rotor(ir)%init(ir, density, dt, nt, switches, &
        & rotor(rotor(ir)%imageRotorNum))
    endif
    call params2file(rotor(ir), nt, dt, nr, density, velSound, switches)
  enddo

  ! Rotate wing pc, vr, cp and nCap by initial pitch angle
  do ir = 1, nr
    do ib = 1, rotor(ir)%nb
      if (rotor(ir)%surfaceType == -1) then
        rotor(ir)%blade(ib)%theta = &
          & rotor(rotor(ir)%imageRotorNum)%blade(ib)%theta
        if (rotor(ir)%imagePlane == 3) then
          rotor(ir)%blade(ib)%theta = -1._dp*rotor(ir)%blade(ib)%theta
        endif
      else
        rotor(ir)%blade(ib)%theta = rotor(ir)%gettheta(rotor(ir)%psiStart, ib)
      endif

      call rotor(ir)%blade(ib)%rot_pitch( &
        sign(1._dp, rotor(ir)%Omega)*rotor(ir)%blade(ib)%theta)
    enddo
  enddo

  ! Plot wing surface geometry
  do ir = 1, nr
    call geomSurface2file(rotor(ir))
  enddo

  ! Compute AIC and AIC_inv matrices for rotors
  do ir = 1, nr
    if (rotor(ir)%surfaceType > 0) then
      call print_status('Computing AIC matrix')
      call rotor(ir)%calcAIC()

      ! Check if inverse was correctly computed
      if (isInverse(rotor(ir)%AIC, rotor(ir)%AIC_inv)) then
        call print_status()    ! SUCCESS
      else
        ! if (abs(rotor(ir)%surfaceType) == 1) then
        stop 'Warning: Computed AIC_inv does not seem &
          & to be correct within given tolerance'
        ! endif
      endif
    endif
  enddo

  ! Initialize vel probes
  if (switches%probe .eq. 1) then
    open(unit=10, file='probes.in', status='old', action='read')
    read(10, *) switches%nProbes
    allocate(probe(3, switches%nProbes))
    allocate(probeVel(3, switches%nProbes))
    do i = 1, switches%nProbes
      read(10, *) probe(:, i), probeVel(:, i)
    enddo
    close(10)
  endif

  ! Obtain initial solution without wake
  call print_status('Computing initial solution')
  if (switches%slowStart .gt. 0) then
    do ir = 1, nr
      rotor(ir)%omegaSlow = 0._dp
    enddo
  else
    do ir = 1, nr
      rotor(ir)%omegaSlow = rotor(ir)%Omega
    enddo
  endif

  t = 0._dp
  iter = 0
  write (timestamp, '(I0.5)') iter

  ! Custom trajectory
  do ir = 1, nr
    if (rotor(ir)%customTrajectorySwitch .eq. 1) then
      rotor(ir)%velBody = rotor(ir)%velBodyHistory(:, 1)
      rotor(ir)%omegaBody = rotor(ir)%omegaBodyHistory(:, 1)
    endif
  enddo

  ! Compute RHS for initial solution without wake
  ntSubInitLoop: do i = 0, switches%ntSubInit
    do ir = 1, nr
      if (rotor(ir)%surfaceType .gt. 0) then
        ! Compute velCP and RHS for lifting and non-lifting surfaces
        !$omp parallel do collapse(3) private(row, jr) schedule(runtime)
        do ib = 1, rotor(ir)%nbConvect
          do is = 1, rotor(ir)%ns
            do ic = 1, rotor(ir)%nc
              row = ic + rotor(ir)%nc*(is - 1) &
                + rotor(ir)%ns*rotor(ir)%nc*(ib - 1)

              ! Translational, rotational, omega, flap vel
              rotor(ir)%blade(ib)%wiP(ic, is)%velCP = &
                & -1._dp*rotor(ir)%velBody &
                & -cross_product(rotor(ir)%omegaBody, &
                & rotor(ir)%blade(ib)%wiP(ic, is)%CP - rotor(ir)%cgCoords) &
                & -cross_product(rotor(ir)%omegaSlow*rotor(ir)%shaftAxis, &
                & rotor(ir)%blade(ib)%wiP(ic, is)%CP - rotor(ir)%hubCoords) &
                & -rotor(ir)%blade(ib)%secMFlapArm(is)* &
                & rotor(ir)%blade(ib)%dflap

              ! Record velocities due to motion for 
              ! computing lift and drag directions
              rotor(ir)%blade(ib)%wiP(ic, is)%velCPm = &
                rotor(ir)%blade(ib)%wiP(ic, is)%velCP

              ! Velocity due to wing vortices of other rotors
              do jr = 1, nr
                if (ir .ne. jr) then
                  rotor(ir)%blade(ib)%wiP(ic, is)%velCP = &
                    rotor(ir)%blade(ib)%wiP(ic, is)%velCP + &
                    rotor(jr)%vind_bywing(rotor(ir)%blade(ib)%wiP(ic, is)%CP)
                endif
              enddo

              rotor(ir)%RHS(row) = &
                dot_product(rotor(ir)%blade(ib)%wiP(ic, is)%velCP, &
                rotor(ir)%blade(ib)%wiP(ic, is)%nCap)

              ! Pitch vel
              !rotor(ir)%blade%(ib)%wing(ic,is)%velPitch= &
              !  rotor(ir)%thetadot_pitch(0._dp,ib)* &
              !  rotor(ir)%blade(ib)%wiP(ic,is)%rHinge
              !rotor(ir)%RHS(row)= RHS(row)+wing(ib,ic,is)%velPitch
            enddo
          enddo
        enddo
        !$omp end parallel do

        axisymRHS0: if (rotor(ir)%axisymmetrySwitch .eq. 1) then
          do ib = 2, rotor(ir)%nb
            rotor(ir)%RHS( &
              & rotor(ir)%nc*rotor(ir)%ns*(ib-1)+1: &
              & rotor(ir)%nc*rotor(ir)%ns*ib) = &
              & rotor(ir)%RHS(1:rotor(ir)%nc*rotor(ir)%ns)
          enddo
        endif axisymRHS0

      else
        ! For image lifting and non-lifting surfaces,
        ! copy velCP and velCPm
        call rotor(ir)%mirrorVelCP(rotor(rotor(ir)%imageRotorNum))
      endif

      rotor(ir)%RHS = -1._dp*rotor(ir)%RHS
    enddo

    do ir = 1, nr
      rotor(ir)%gamVecPrev = rotor(ir)%gamVec
      if (rotor(ir)%surfaceType .gt. 0) then
        rotor(ir)%gamVec = matmulAX(rotor(ir)%AIC_inv, rotor(ir)%RHS)
      else  ! Image lifting or non-lifting surface
        call rotor(ir)%mirrorGamma(rotor(rotor(ir)%imageRotorNum))
      endif

      ! Map gamVec to wing gam for each blade in rotor
      call rotor(ir)%map_gam()
    enddo

    ! Check if initialization using converged soultion is requested
    if (switches%ntSubInit .ne. 0) then
      subIterResidual = 0._dp
      do ir = 1, nr
        subIterResidual = max(subIterResidual, &
          & norm2(rotor(ir)%gamVec - rotor(ir)%gamVecPrev))
      enddo
      do ir = 1, nr
        if (subIterResidual .le. eps) &
          exit ntSubInitLoop
      enddo
    endif
  enddo ntSubInitLoop

  if (switches%ntSubInit .ne. 0) then
    if (i .gt. switches%ntSubInit) then
      print*
      print*, "Initial solution did not converge. & 
        & Try increasing sub-iterations."
      print*, 'Sub-iterations ', i, subIterResidual
    else
      call print_status()    ! SUCCESS
      print*, 'Sub-iterations ', i, subIterResidual
    endif
  else
    call print_status()    ! SUCCESS
  endif

  do ir = 1, nr
    rotor(ir)%rowFar = rotor(ir)%nFwake + 1
    ! Since assignshed TE assigns to rowNear-1 panel
    rotor(ir)%rowNear = rotor(ir)%nNwake + 1
    if (rotor(ir)%nNwake > 0) then
      call rotor(ir)%assignshed('TE')  ! Store shed vortex as TE
    endif
  enddo

  ! Compute forces

  if (switches%rotorForcePlot .ne. 0) then
    call init_plots(nr)    ! Create headers for plot files
    do ir = 1, nr

      select case (rotor(ir)%forceCalcSwitch)

      case (0)  ! Compute using wing circulation
        ! Compute and plot alpha if requested
        ! Compute alpha
        !$omp parallel do collapse(3) private(jr) schedule(runtime)
        do ib = 1, rotor(ir)%nbConvect
          do is = 1, rotor(ir)%ns
            do ic = 1, rotor(ir)%nc
              ! Compute local velocity vector
              ! (excluding induced velocities from wing bound vortices)
              rotor(ir)%blade(ib)%wiP(ic, is)%velCPTotal = &
                rotor(ir)%blade(ib)%wiP(ic, is)%velCP

              ! Neglect velocity due to spanwise vortices for all wings
              do jr = 1, nr
                rotor(ir)%blade(ib)%wip(ic, is)%velCPTotal = &
                  rotor(ir)%blade(ib)%wip(ic, is)%velCPTotal - &
                  rotor(jr)%vind_bywing_boundVortices( &
                  rotor(ir)%blade(ib)%wiP(ic, is)%CP)
              enddo

              ! Add self induced velocity due to wing vortices
              rotor(ir)%blade(ib)%wiP(ic, is)%velCPTotal = &
                rotor(ir)%blade(ib)%wiP(ic, is)%velCPTotal + &
                rotor(ir)%vind_bywing(rotor(ir)%blade(ib)%wiP(ic, is)%CP)
            enddo
          enddo
        enddo
        !$omp end parallel do

        axisym00: if (rotor(ir)%axisymmetrySwitch .eq. 1) then
          do ib = 2, rotor(ir)%nb
            do is = 1, rotor(ir)%ns
              do ic = 1, rotor(ir)%nc
                rotor(ir)%blade(ib)%wiP(ic, is)%velCPTotal = &
                  & rotor(ir)%blade(1)%wiP(ic, is)%velCPTotal
              enddo
            enddo
          enddo
        endif axisym00

        call rotor(ir)%calc_secAlpha()
        call rotor(ir)%calc_force(density, dt)

        ! For the first iteration, assign the first flap moment to 
        ! prev flap moment for use in flap dynamics equation
        do ib = 1, rotor(ir)%nb
          rotor(ir)%blade(ib)%MflapLiftPrev = rotor(ir)%blade(ib)%MflapLift
        enddo

      case (1)  ! Compute using alpha
        ! Compute alpha
        !$omp parallel do collapse(3) private(jr) schedule(runtime)
        do ib = 1, rotor(ir)%nbConvect
          do is = 1, rotor(ir)%ns
            do ic = 1, rotor(ir)%nc
              ! Compute local velocity vector
              ! (excluding induced velocities from wing bound vortices)
              rotor(ir)%blade(ib)%wiP(ic, is)%velCPTotal = &
                rotor(ir)%blade(ib)%wiP(ic, is)%velCP

              ! Neglect velocity due to spanwise vortices for all wings
              do jr = 1, nr
                rotor(ir)%blade(ib)%wip(ic, is)%velCPTotal = &
                  rotor(ir)%blade(ib)%wip(ic, is)%velCPTotal - &
                  rotor(jr)%vind_bywing_boundVortices( &
                  rotor(ir)%blade(ib)%wiP(ic, is)%CP)
              enddo

              ! Add self induced velocity due to wing vortices
              rotor(ir)%blade(ib)%wiP(ic, is)%velCPTotal = &
                rotor(ir)%blade(ib)%wiP(ic, is)%velCPTotal + &
                rotor(ir)%vind_bywing(rotor(ir)%blade(ib)%wiP(ic, is)%CP)
            enddo
          enddo
        enddo
        !$omp end parallel do

        axisym01: if (rotor(ir)%axisymmetrySwitch .eq. 1) then
          do ib = 2, rotor(ir)%nb
            do is = 1, rotor(ir)%ns
              do ic = 1, rotor(ir)%nc
                rotor(ir)%blade(ib)%wiP(ic, is)%velCPTotal = &
                  & rotor(ir)%blade(1)%wiP(ic, is)%velCPTotal
              enddo
            enddo
          enddo
        endif axisym01

        call rotor(ir)%calc_secAlpha()
        call rotor(ir)%calc_force_alpha(density, velSound)

      case (2)  ! Compute lift using alpha approximated from sec circulation
        call rotor(ir)%calc_force_alphaGamma(density, velSound, dt)

      end select

      ! Initial force value
      call force2file(timestamp, rotor(ir))

      ! Flap dynamics
      if (rotor(ir)%bladeDynamicsSwitch .ne. 0) then
        call dynamics2file(timestamp, rotor(ir))
        if (rotor(ir)%surfaceType == -1) then
          do ib = 1, rotor(ir)%nb
            if (rotor(ir)%imagePlane == 3) then
              rotor(ir)%blade(ib)%dflap = &
                & -1._dp*rotor(rotor(ir)%imageRotorNum)%blade(ib)%dflap
              rotor(ir)%blade(ib)%flap = &
                & -1._dp*rotor(rotor(ir)%imageRotorNum)%blade(ib)%flap
              rotor(ir)%blade(ib)%flapPrev = &
                & -1._dp*rotor(rotor(ir)%imageRotorNum)%blade(ib)%flapPrev
            else
              rotor(ir)%blade(ib)%dflap = &
                & rotor(rotor(ir)%imageRotorNum)%blade(ib)%dflap
              rotor(ir)%blade(ib)%flap = &
                & rotor(rotor(ir)%imageRotorNum)%blade(ib)%flap
              rotor(ir)%blade(ib)%flapPrev = &
                & rotor(rotor(ir)%imageRotorNum)%blade(ib)%flapPrev
            endif
          enddo
        else
          call rotor(ir)%computeBladeDynamics(dt)
        endif
      endif

      ! Body dynamics
      if (rotor(ir)%bodyDynamicsSwitch .ne. 0) then
        call dynamics2file(timestamp, rotor(ir))
        if (rotor(ir)%surfaceType == -1) then
          rotor(ir)%velBody =  rotor(rotor(ir)%imageRotorNum)%velBody
          rotor(ir)%velBody(rotor(ir)%imagePlane) = &
            & -1._dp*rotor(ir)%velBody(rotor(ir)%imagePlane)
        else
          call rotor(ir)%computeBodyDynamics(dt)
        endif
      endif

    enddo
  endif

  currentTime = ''

  ! ------- MAIN LOOP START -------

  iterStart = 1

  if (switches%restartFromNt .gt. 0) then
    write (timestamp, '(I0.5)') switches%restartFromNt
    open(unit=23, file='Restart/restart'//timestamp//'.dat', &
      & status='old', action='read', form='unformatted')
    read(23) t
    read(23) rotor
    close(23)
    iterStart = switches%restartFromNt + 1
  endif

  do iter = iterStart, nt
    t = t + dt

    ! In case current time is required uncomment below line
    ! call date_and_time(time=currentTime)

    print *, currentTime, iter, nt

    write (timestamp, '(I0.5)') iter

    ! rowNear and rowFar keep track of what row
    ! the current iteration is in for near wake and far wake
    do ir = 1, nr
      rotor(ir)%rowNear = max(rotor(ir)%rowNear - 1, 1)
      if (iter > rotor(ir)%nNwake) then
        rotor(ir)%rowFar = max(rotor(ir)%rowFar - 1, 1)
      endif
    enddo

    ! Use custom trajectory if specified
    do ir = 1, nr
      if (rotor(ir)%customTrajectorySwitch .eq. 1) then
        rotor(ir)%velBody = rotor(ir)%velBodyHistory(:, iter)
        rotor(ir)%omegaBody = rotor(ir)%omegaBodyHistory(:, iter)
      endif
    enddo

    ! In case of slow start, determine RPM
    select case (switches%slowStart)
    case (-1)   ! Body dynamics governed
      continue
    case (0)    ! No slow start
      do ir = 1, nr
        rotor(ir)%omegaSlow = rotor(ir)%Omega
      enddo
    case (1)    ! Linear slope
      do ir = 1, nr
        rotor(ir)%omegaSlow = &
          min(real(switches%slowStartNt), real(iter + 1))*rotor(ir)%Omega/switches%slowStartNt
      enddo
    case (2)    ! tanh slope
      do ir = 1, nr
        rotor(ir)%omegaSlow = tanh(5._dp*iter/switches%slowStartNt)*rotor(ir)%Omega
      enddo
    case (3)    ! tanh slope
      do ir = 1, nr
        rotor(ir)%omegaSlow = &
          (tanh(6._dp*real(iter)/real(switches%slowStartNt) - 3._dp) + 1._dp)* &
          0.5_dp*rotor(ir)%Omega
      enddo
    case default
      error stop "Assign correct slowStart"
    end select

    ! Move wing to next position
    do ir = 1, nr
      call rotor(ir)%move(rotor(ir)%velBody*dt)
      call rotor(ir)%rot_pts(rotor(ir)%omegaBody*dt, rotor(ir)%cgCoords, 1)
      call rotor(ir)%rot_advance(rotor(ir)%omegaSlow*dt)

      if (rotor(ir)%bladeDynamicsSwitch == 1) then
        call rotor(ir)%rot_flap()
      endif
    enddo

    ! Assign LE of near wake
    if (switches%wakeSuppress == 0) then
      do ir = 1, nr
        if (rotor(ir)%nNwake > 0) then
          call rotor(ir)%assignshed('LE')  ! Store shed vortex as LE
        endif
      enddo

      ! Update wake age
      do ir = 1, nr
        if (rotor(ir)%nNwake > 0) then
          call rotor(ir)%age_wake(dt)
        endif
      enddo

      ! Dissipate wake
      if (switches%wakeDissipation .eq. 1) then
        do ir = 1, nr
          if (rotor(ir)%nNwake > 0) then
            ! Wake tip dissipation
            call rotor(ir)%dissipate_wake(dt, kinematicVisc)
          endif
        enddo
      endif

      ! Burst wake
      do ir = 1, nr
        if (rotor(ir)%nNwake > 0) then
          if (switches%wakeBurst .ne. 0) then
            if (mod(iter, switches%wakeBurst) .eq. 0) &
              call rotor(ir)%burst_wake()
          endif
          ! Plot wake skew parameter
          if (rotor(ir)%skewPlotSwitch .ne. 0) then
            if (mod(iter, rotor(ir)%skewPlotSwitch) .eq. 0) then
              call rotor(ir)%calc_skew()
              call skew2file(timestamp, rotor(ir))
            endif
          endif
        endif
      enddo
    endif

    ! Write out wing n' wake
    do ir = 1, nr
      if (switches%wakePlot .ne. 0) then
        if (mod(iter, switches%wakePlot) .eq. 0) &
          call geom2file(timestamp, rotor(ir), switches%wakeSuppress)
      endif

      if (switches%wakeTipPlot .ne. 0) then
        if (mod(iter, switches%wakeTipPlot) .eq. 0) &
          call tip2file(timestamp, rotor(ir))
      endif
    enddo

    ! Compute RHS
    ntSubLoop: do i = 0, switches%ntSub
      do ir = 1, nr

        rotor(ir)%RHS = 0._dp
        if (rotor(ir)%surfaceType .gt. 0) then
          ! Compute velCP and RHS for lifting and non-lifting surfaces
          !$omp parallel do collapse(3) private(row, jr) schedule(runtime)
          do ib = 1, rotor(ir)%nbConvect
            do is = 1, rotor(ir)%ns
              do ic = 1, rotor(ir)%nc
                row = ic + rotor(ir)%nc*(is - 1) &
                  + rotor(ir)%ns*rotor(ir)%nc*(ib - 1)

                ! Translational, rotationa, omega, flap vel
                rotor(ir)%blade(ib)%wiP(ic, is)%velCP = &
                  & -1._dp*rotor(ir)%velBody &
                  & -cross_product(rotor(ir)%omegaBody, &
                  & rotor(ir)%blade(ib)%wiP(ic, is)%CP - rotor(ir)%cgCoords) &
                  & -cross_product(rotor(ir)%omegaSlow*rotor(ir)%shaftAxis, &
                  & rotor(ir)%blade(ib)%wiP(ic, is)%CP - rotor(ir)%hubCoords) &
                  & -rotor(ir)%blade(ib)%secMFlapArm(is)* &
                  & rotor(ir)%blade(ib)%dflap

                ! Record velocities due to motion for induced drag computation
                rotor(ir)%blade(ib)%wiP(ic, is)%velCPm = &
                  rotor(ir)%blade(ib)%wiP(ic, is)%velCP

                do jr = 1, nr
                  ! Wake induced vel due to all rotors
                  rotor(ir)%blade(ib)%wiP(ic, is)%velCP = &
                    rotor(ir)%blade(ib)%wiP(ic, is)%velCP + &
                    rotor(jr)%vind_bywake(rotor(ir)%blade(ib)%wiP(ic, is)%CP)

                  ! Wing induced vel due to all rotors except self
                  if (ir .ne. jr) then
                    rotor(ir)%blade(ib)%wiP(ic, is)%velCP = &
                      rotor(ir)%blade(ib)%wiP(ic, is)%velCP + &
                      rotor(jr)%vind_bywing(rotor(ir)%blade(ib)%wiP(ic, is)%CP)
                  endif
                enddo

                rotor(ir)%RHS(row) = &
                  dot_product(rotor(ir)%blade(ib)%wiP(ic, is)%velCP, &
                  rotor(ir)%blade(ib)%wiP(ic, is)%nCap)

                ! Pitch vel
                !wing(ib,ic,is)%velPitch=thetadot*wing(ib,ic,is)%rHinge
                !RHS(row)=RHS(row)+wing(ib,ic,is)%velPitch
              enddo
            enddo
          enddo
          !$omp end parallel do

          axisymRHS: if (rotor(ir)%axisymmetrySwitch .eq. 1) then
            do ib = 2, rotor(ir)%nb
              rotor(ir)%RHS( &
                & rotor(ir)%nc*rotor(ir)%ns*(ib-1)+1: &
                & rotor(ir)%nc*rotor(ir)%ns*ib) = &
                & rotor(ir)%RHS(1:rotor(ir)%nc*rotor(ir)%ns)
            enddo
          endif axisymRHS

        else
          ! For image lifting and non-lifting surfaces,
          ! copy velCP and velCPm
          call rotor(ir)%mirrorVelCP(rotor(rotor(ir)%imageRotorNum))
        endif

        rotor(ir)%RHS = -1._dp*rotor(ir)%RHS
      enddo

      do ir = 1, nr
        rotor(ir)%gamvecPrev = rotor(ir)%gamVec
        if (rotor(ir)%surfaceType .gt. 0) then
          rotor(ir)%gamVec = matmulAX(rotor(ir)%AIC_inv, rotor(ir)%RHS)
        else
          call rotor(ir)%mirrorGamma(rotor(rotor(ir)%imageRotorNum))
        endif

        ! Map gamVec to wing gam for each blade in rotor
        call rotor(ir)%map_gam()
      enddo

      ! Check if sub-iterations are requested
      if (switches%ntSub .ne. 0) then
        subIterResidual = 0._dp
        do ir = 1, nr
          subIterResidual = max(subIterResidual, &
            & norm2(rotor(ir)%gamVec - rotor(ir)%gamVecPrev))
          if (subIterResidual .le. eps) &
            & exit ntSubLoop
        enddo
      endif
    enddo ntSubLoop

    if (switches%ntSub .ne. 0) then
      print *, 'Sub-iterations ', i, subIterResidual
      if (i .gt. switches%ntSub) &
        & print*, 'Solution did not converge. Try increasing sub-iterations.'
    endif

    ! Compute forces
    if (switches%rotorForcePlot .ne. 0) then
      if (mod(iter, switches%rotorForcePlot) .eq. 0) then
        do ir = 1, nr

          select case (rotor(ir)%forceCalcSwitch)

          case (0)  ! Compute using wing circulation
            ! Compute alpha
            !$omp parallel do collapse(3) private(jr) schedule(runtime)
            do ib = 1, rotor(ir)%nbConvect
              do is = 1, rotor(ir)%ns
                do ic = 1, rotor(ir)%nc
                  ! Compute local velocity vector
                  ! (excluding induced velocities from wing bound vortices)
                  rotor(ir)%blade(ib)%wiP(ic, is)%velCPTotal = &
                    & rotor(ir)%blade(ib)%wiP(ic, is)%velCP

                  ! Neglect velocity due to spanwise vortices for all wings
                  do jr = 1, nr
                    rotor(ir)%blade(ib)%wip(ic, is)%velCPTotal = &
                      & rotor(ir)%blade(ib)%wip(ic, is)%velCPTotal - &
                      & rotor(jr)%vind_bywing_boundVortices( &
                      & rotor(ir)%blade(ib)%wiP(ic, is)%CP)
                  enddo

                  ! Add self induced velocity due to wing vortices
                  rotor(ir)%blade(ib)%wiP(ic, is)%velCPTotal = &
                    & rotor(ir)%blade(ib)%wiP(ic, is)%velCPTotal + &
                    & rotor(ir)%vind_bywing(rotor(ir)%blade(ib)%wiP(ic, is)%CP)
                enddo
              enddo
            enddo
            !$omp end parallel do

            axisymx0: if (rotor(ir)%axisymmetrySwitch .eq. 1) then
              do ib = 2, rotor(ir)%nb
                do is = 1, rotor(ir)%ns
                  do ic = 1, rotor(ir)%nc
                    rotor(ir)%blade(ib)%wiP(ic, is)%velCPTotal = &
                      & rotor(ir)%blade(1)%wiP(ic, is)%velCPTotal
                  enddo
                enddo
              enddo
            endif axisymx0

            call rotor(ir)%calc_secAlpha()
            call rotor(ir)%calc_force(density, dt)

          case (1)  ! Compute using alpha
            ! Compute alpha
            !$omp parallel do collapse(3) private(jr) schedule(runtime)
            do ib = 1, rotor(ir)%nbConvect
              do is = 1, rotor(ir)%ns
                do ic = 1, rotor(ir)%nc
                  ! Compute local velocity vector
                  ! (excluding induced velocities from wing bound vortices)
                  rotor(ir)%blade(ib)%wiP(ic, is)%velCPTotal = &
                    rotor(ir)%blade(ib)%wiP(ic, is)%velCP

                  ! Neglect velocity due to spanwise vortices for all wings
                  do jr = 1, nr
                    rotor(ir)%blade(ib)%wip(ic, is)%velCPTotal = &
                      rotor(ir)%blade(ib)%wip(ic, is)%velCPTotal - &
                      rotor(jr)%vind_bywing_boundVortices( &
                      rotor(ir)%blade(ib)%wiP(ic, is)%CP)
                  enddo

                  ! Add self induced velocity due to wing vortices
                  rotor(ir)%blade(ib)%wiP(ic, is)%velCPTotal = &
                    rotor(ir)%blade(ib)%wiP(ic, is)%velCPTotal + &
                    rotor(ir)%vind_bywing(rotor(ir)%blade(ib)%wiP(ic, is)%CP)
                enddo
              enddo
            enddo
            !$omp end parallel do

            axisymx1: if (rotor(ir)%axisymmetrySwitch .eq. 1) then
              do ib = 2, rotor(ir)%nb
                do is = 1, rotor(ir)%ns
                  do ic = 1, rotor(ir)%nc
                    rotor(ir)%blade(ib)%wiP(ic, is)%velCPTotal = &
                      & rotor(ir)%blade(1)%wiP(ic, is)%velCPTotal
                  enddo
                enddo
              enddo
            endif axisymx1

            call rotor(ir)%calc_secAlpha()
            call rotor(ir)%calc_force_alpha(density, velSound)

          case (2)  ! Compute lift using alpha approximated from sec circulation
            call rotor(ir)%calc_force_alphaGamma(density, velSound, dt)

          end select

          call force2file(timestamp, rotor(ir))

        enddo
      endif

      ! Flap dynamics
      do ir = 1, nr
        if (rotor(ir)%bladeDynamicsSwitch .ne. 0) then
          call dynamics2file(timestamp, rotor(ir))
          if (rotor(ir)%surfaceType == -1) then
            do ib = 1, rotor(ir)%nb
              if (rotor(ir)%imagePlane == 3) then
                rotor(ir)%blade(ib)%dflap = &
                  & -1._dp*rotor(rotor(ir)%imageRotorNum)%blade(ib)%dflap
                rotor(ir)%blade(ib)%flap = &
                  & -1._dp*rotor(rotor(ir)%imageRotorNum)%blade(ib)%flap
                rotor(ir)%blade(ib)%flapPrev = &
                  & -1._dp*rotor(rotor(ir)%imageRotorNum)%blade(ib)%flapPrev
              else
                rotor(ir)%blade(ib)%dflap = &
                  & rotor(rotor(ir)%imageRotorNum)%blade(ib)%dflap
                rotor(ir)%blade(ib)%flap = &
                  & rotor(rotor(ir)%imageRotorNum)%blade(ib)%flap
                rotor(ir)%blade(ib)%flapPrev = &
                  & rotor(rotor(ir)%imageRotorNum)%blade(ib)%flapPrev
              endif
            enddo
          else
            call rotor(ir)%computeBladeDynamics(dt)
          endif
        endif

        ! Body dynamics
        if (rotor(ir)%bodyDynamicsSwitch .ne. 0) then
          call dynamics2file(timestamp, rotor(ir))
          if (rotor(ir)%surfaceType == -1) then
            rotor(ir)%velBody =  rotor(rotor(ir)%imageRotorNum)%velBody
            rotor(ir)%velBody(rotor(ir)%imagePlane) = &
              & -1._dp*rotor(ir)%velBody(rotor(ir)%imagePlane)
          else
            call rotor(ir)%computeBodyDynamics(dt)
          endif
        endif

      enddo
    endif


    ! Plot inflow
    do ir = 1, nr
      if (rotor(ir)%inflowPlotSwitch .ne. 0) then
        if (mod(iter, rotor(ir)%inflowPlotSwitch) .eq. 0) then
          call inflow2file(timestamp, rotor, ir, -zAxis)
        endif
      endif
    enddo

    ! Plot sec wing bound circulation
    do ir = 1, nr
      if (rotor(ir)%gammaPlotSwitch .ne. 0) then
        if (mod(iter, rotor(ir)%gammaPlotSwitch) .eq. 0) then
          call gamma2file(timestamp, rotor(1))
        endif
      endif
    enddo

    ! Record filaments for grid plot computation
    if (switches%gridPlot .ne. 0) then
      if (mod(iter, switches%gridPlot) .eq. 0) then
        call filaments2file(timestamp, rotor)
      endif
    endif

    ! Compute probe velocities
    if (switches%probe .ne. 0) then
      if (mod(iter, switches%probe) .eq. 0) then
        call probes2file(timestamp, probe, probeVel, rotor, t)
      endif
    endif

    ! Wake convection
    if (switches%wakeSuppress == 0) then
      ! Initialise wake velocity matrices
      do ir = 1, nr
        if (rotor(ir)%nNwake > 0) then
          !$omp parallel do schedule(runtime)
          do ib = 1, rotor(ir)%nbConvect
            rotor(ir)%blade(ib)%velNwake(:, rotor(ir)%rowNear:rotor(ir)%nNwake, :) = 0._dp
            rotor(ir)%blade(ib)%velFwake(:, rotor(ir)%rowFar:rotor(ir)%nFwake) = 0._dp
          enddo
          !$omp end parallel do
        endif
      enddo

      ! Compute induced velocity due to rotors in domain
      do ir = 1, nr
        if (rotor(ir)%nNwake > 0) then
          do ib = 1, rotor(ir)%nbConvect
            do jr = 1, nr
              rotor(ir)%blade(ib)%velNwake(:, rotor(ir)%rowNear:rotor(ir)%nNwakeEnd, :) = &
                rotor(ir)%blade(ib)%velNwake(:, rotor(ir)%rowNear:rotor(ir)%nNwakeEnd, :) + &
                vind_onNwake_byRotor(rotor(jr), &
                rotor(ir)%blade(ib)%waN(rotor(ir)%rowNear:rotor(ir)%nNwakeEnd, :))
              rotor(ir)%blade(ib)%velFwake(:, rotor(ir)%rowFar:rotor(ir)%nFwakeEnd) = &
                rotor(ir)%blade(ib)%velFwake(:, rotor(ir)%rowFar:rotor(ir)%nFwakeEnd) + &
                vind_onFwake_byRotor(rotor(jr), &
                rotor(ir)%blade(ib)%waF(rotor(ir)%rowFar:rotor(ir)%nFwakeEnd))
            enddo

            ! Add initial wake velocity if provided
            if (iter < switches%initWakeVelNt) then
              do i = 1, 3
                rotor(ir)%blade(ib)%velNwake(i, rotor(ir)%rowNear:rotor(ir)%nNwakeEnd, :) = &
                  rotor(ir)%blade(ib)%velNwake(i, rotor(ir)%rowNear:rotor(ir)%nNwakeEnd, :) + &
                  rotor(ir)%initWakeVel*rotor(ir)%shaftAxis(i)
                rotor(ir)%blade(ib)%velFwake(i, rotor(ir)%rowFar:rotor(ir)%nFwakeEnd) = &
                  rotor(ir)%blade(ib)%velFwake(i, rotor(ir)%rowFar:rotor(ir)%nFwakeEnd) &
                  - rotor(ir)%initWakeVel*rotor(ir)%shaftAxis(i)
              enddo
            endif
          enddo
        endif
      enddo

      ! Update wake vortex locations if lifting surface
      select case (switches%fdScheme)

      case (0)    ! Explicit Forward Diff (1st order)
        do ir = 1, nr
          if (rotor(ir)%nNwake > 0) then
            if (abs(rotor(ir)%surfaceType) == 1) then
              if (rotor(ir)%surfaceType .gt. 0) then
                ! Lifting surface
                call rotor(ir)%convectwake(iter, dt, 'C')
              else
                ! Image lifting surface
                call rotor(ir)%mirrorWake(rotor(rotor(ir)%imageRotorNum), 'C')
              endif
            endif
          endif
        enddo

      case (1)    ! Predictor-Corrector (2nd order)
        ! Compute predicted wake
        do ir = 1, nr
          if (rotor(ir)%nNwake > 0) then
            if (abs(rotor(ir)%surfaceType) == 1) then
              if (rotor(ir)%surfaceType .gt. 0) then
                ! Lifting surface
                do ib = 1, rotor(ir)%nbConvect
                  rotor(ir)%blade(ib)%waNPredicted(rotor(ir)%rowNear:rotor(ir)%nNwakeEnd, :) = &
                    rotor(ir)%blade(ib)%waN(rotor(ir)%rowNear:rotor(ir)%nNwakeEnd, :)
                  rotor(ir)%blade(ib)%waFPredicted(rotor(ir)%rowFar:rotor(ir)%nFwakeEnd) = &
                    rotor(ir)%blade(ib)%waF(rotor(ir)%rowFar:rotor(ir)%nFwakeEnd)
                enddo
                call rotor(ir)%convectwake(iter, dt, 'P')
              else
                ! Image lifting surface
                call rotor(ir)%mirrorWake(rotor(rotor(ir)%imageRotorNum), 'P')
              endif
            endif
          endif
        enddo

        ! Compute velocity on predicted wake
        do ir = 1, nr
          if (rotor(ir)%nNwake > 0) then
            if (abs(rotor(ir)%surfaceType) == 1) then
              if (rotor(ir)%surfaceType .gt. 0) then
                ! Lifting surface
                do ib = 1, rotor(ir)%nbConvect
                  rotor(ir)%blade(ib)%velNwakePredicted(:, rotor(ir)%rowNear:rotor(ir)%nNwakeEnd, :) = 0._dp
                  rotor(ir)%blade(ib)%velFwakePredicted(:, rotor(ir)%rowFar:rotor(ir)%nFwakeEnd) = 0._dp
                  do jr = 1, nr
                    rotor(ir)%blade(ib)%velNwakePredicted(:, rotor(ir)%rowNear:rotor(ir)%nNwakeEnd, :) = &
                      rotor(ir)%blade(ib)%velNwakePredicted(:, rotor(ir)%rowNear:rotor(ir)%nNwakeEnd, :) + &
                      vind_onNwake_byRotor(rotor(jr), &
                      rotor(ir)%blade(ib)%waNPredicted(rotor(ir)%rowNear:rotor(ir)%nNwakeEnd, :), 'P')
                    rotor(ir)%blade(ib)%velFwakePredicted(:, rotor(ir)%rowFar:rotor(ir)%nFwakeEnd) = &
                      rotor(ir)%blade(ib)%velFwakePredicted(:, rotor(ir)%rowFar:rotor(ir)%nFwakeEnd) + &
                      vind_onFwake_byRotor(rotor(jr), &
                      rotor(ir)%blade(ib)%waFPredicted(rotor(ir)%rowFar:rotor(ir)%nFwakeEnd), 'P')
                  enddo
                  if (iter < switches%initWakeVelNt) then
                    do i = 1, 3
                      rotor(ir)%blade(ib)%velNwakePredicted(i, rotor(ir)%rowNear:rotor(ir)%nNwakeEnd, :) = &
                        rotor(ir)%blade(ib)%velNwakePredicted(i, rotor(ir)%rowNear:rotor(ir)%nNwakeEnd, :) - &
                        rotor(ir)%initWakeVel*rotor(ir)%shaftAxis(i)
                    enddo
                    do i = 1, 3
                      rotor(ir)%blade(ib)%velFwakePredicted(i, rotor(ir)%rowFar:rotor(ir)%nFwakeEnd) = &
                        rotor(ir)%blade(ib)%velFwakePredicted(i, rotor(ir)%rowFar:rotor(ir)%nFwakeEnd) - &
                        rotor(ir)%initWakeVel*rotor(ir)%shaftAxis(i)
                    enddo
                  endif
                enddo
              endif
            endif
          endif
        enddo

        ! Compute averaged velocity and convect wake
        do ir = 1, nr
          if (rotor(ir)%nNwake > 0) then
            if (abs(rotor(ir)%surfaceType) == 1) then
              if (rotor(ir)%surfaceType .gt. 0) then
                ! Lifting surface
                do ib = 1, rotor(ir)%nbConvect
                  rotor(ir)%blade(ib)%velNwake(:, rotor(ir)%rowNear:rotor(ir)%nNwakeEnd, :) = &
                    vel_order2_Nwake(rotor(ir)%blade(ib)%velNwake(:, &
                    rotor(ir)%rowNear:rotor(ir)%nNwakeEnd, :), &
                    rotor(ir)%blade(ib)%velNwakePredicted(:, rotor(ir)%rowNear:rotor(ir)%nNwakeEnd, :))
                  !rotor(ir)%blade(ib)%velNwake(:,rotor(ir)%rowNear:rotor(ir)%nNwakeEnd,:)= &
                  !0.5_dp*(rotor(ir)%blade(ib)%velNwake(:,rotor(ir)%rowNear:rotor(ir)%nNwakeEnd,:)+ &
                  !rotor(ir)%blade(ib)%velNwakePredicted(:,rotor(ir)%rowNear:rotor(ir)%nNwakeEnd,:))
                  rotor(ir)%blade(ib)%velFwake(:, rotor(ir)%rowFar:rotor(ir)%nFwakeEnd) = &
                    vel_order2_Fwake(rotor(ir)%blade(ib)%velFwake(:, &
                    rotor(ir)%rowFar:rotor(ir)%nFwakeEnd), &
                    rotor(ir)%blade(ib)%velFwakePredicted(:, rotor(ir)%rowFar:rotor(ir)%nFwakeEnd))
                  !rotor(ir)%blade(ib)%velFwake(:,rotor(ir)%rowFar:rotor(ir)%nFwakeEnd)= &
                  !0.5_dp*(rotor(ir)%blade(ib)%velFwake(:,rotor(ir)%rowFar:rotor(ir)%nFwakeEnd)+ &
                  !rotor(ir)%blade(ib)%velFwakePredicted(:,rotor(ir)%rowFar:rotor(ir)%nFwakeEnd))
                enddo
                call rotor(ir)%convectwake(iter, dt, 'C')
              else
                ! Image lifting surface
                call rotor(ir)%mirrorWake(rotor(rotor(ir)%imageRotorNum), 'C')
              endif
            endif
          endif
        enddo

      case (2)    ! Explicit Adams-Bashforth (2nd order)
        if (iter == 1) then
          do ir = 1, nr
            if (rotor(ir)%nNwake > 0) then
              if (abs(rotor(ir)%surfaceType) == 1) then
                if (rotor(ir)%surfaceType .gt. 0) then
                  ! Lifting surface
                  call rotor(ir)%convectwake(iter, dt, 'C')
                  do ib = 1, rotor(ir)%nbConvect
                    rotor(ir)%blade(ib)%velNwake1 = rotor(ir)%blade(ib)%velNwake
                    rotor(ir)%blade(ib)%velFwake1 = rotor(ir)%blade(ib)%velFwake
                  enddo
                else
                  ! Image lifting surface
                  call rotor(ir)%mirrorWake(rotor(rotor(ir)%imageRotorNum), 'C')
                endif
              endif
            endif
          enddo
        else
          do ir = 1, nr
            if (rotor(ir)%nNwake > 0) then
              if (abs(rotor(ir)%surfaceType) == 1) then
                if (rotor(ir)%surfaceType .gt. 0) then
                  ! Lifting surface
                  do ib = 1, rotor(ir)%nbConvect
                    rotor(ir)%blade(ib)%velNwakeStep = &
                      0.5_dp*(3._dp*rotor(ir)%blade(ib)%velNwake - &
                      rotor(ir)%blade(ib)%velNwake1)
                    rotor(ir)%blade(ib)%velFwakeStep = &
                      0.5_dp*(3._dp*rotor(ir)%blade(ib)%velFwake - &
                      rotor(ir)%blade(ib)%velFwake1)

                    ! For next step
                    rotor(ir)%blade(ib)%velNwake1 = rotor(ir)%blade(ib)%velNwakeStep
                    rotor(ir)%blade(ib)%velFwake1 = rotor(ir)%blade(ib)%velFwakeStep

                    ! For convection
                    rotor(ir)%blade(ib)%velNwake = rotor(ir)%blade(ib)%velNwakeStep
                    rotor(ir)%blade(ib)%velFwake = rotor(ir)%blade(ib)%velFwakeStep
                  enddo
                  call rotor(ir)%convectwake(iter, dt, 'C')
                else
                  ! Image lifting surface
                  call rotor(ir)%mirrorWake(rotor(rotor(ir)%imageRotorNum), 'C')
                endif
              endif
            endif
          enddo
        endif

      case (3)    ! Predictor-Corrector Adams-Moulton (2nd order)
        if (iter == 1) then
          do ir = 1, nr
            if (rotor(ir)%nNwake > 0) then
              if (abs(rotor(ir)%surfaceType) == 1) then
                if (rotor(ir)%surfaceType .gt. 0) then
                  ! Lifting surface
                  call rotor(ir)%convectwake(iter, dt, 'C')
                  do ib = 1, rotor(ir)%nbConvect
                    rotor(ir)%blade(ib)%velNwake1 = rotor(ir)%blade(ib)%velNwake
                    rotor(ir)%blade(ib)%velFwake1 = rotor(ir)%blade(ib)%velFwake
                  enddo
                else
                  ! Image lifting surface
                  call rotor(ir)%mirrorWake(rotor(rotor(ir)%imageRotorNum), 'C')
                endif
              endif
            endif
          enddo
        else
          do ir = 1, nr
            if (rotor(ir)%nNwake > 0) then
              if (abs(rotor(ir)%surfaceType) == 1) then
                if (rotor(ir)%surfaceType .gt. 0) then
                  ! Lifting surface
                  do ib = 1, rotor(ir)%nbConvect
                    rotor(ir)%blade(ib)%waNPredicted(rotor(ir)%rowNear:rotor(ir)%nNwakeEnd, :) = &
                      rotor(ir)%blade(ib)%waN(rotor(ir)%rowNear:rotor(ir)%nNwakeEnd, :)
                    ! Store Nwake to Nwake_step for later use
                    rotor(ir)%blade(ib)%velNwakeStep = rotor(ir)%blade(ib)%velNwake
                    rotor(ir)%blade(ib)%velNwake = &
                      0.5_dp*(3._dp*rotor(ir)%blade(ib)%velNwake - &
                      rotor(ir)%blade(ib)%velNwake1)
                    rotor(ir)%blade(ib)%waFPredicted(rotor(ir)%rowFar:rotor(ir)%nFwakeEnd) = &
                      rotor(ir)%blade(ib)%waF(rotor(ir)%rowFar:rotor(ir)%nFwakeEnd)
                    ! Store Fwake to Fwake_step for later use
                    rotor(ir)%blade(ib)%velFwakeStep = rotor(ir)%blade(ib)%velFwake
                    rotor(ir)%blade(ib)%velFwake = &
                      0.5_dp*(3._dp*rotor(ir)%blade(ib)%velFwake - &
                      rotor(ir)%blade(ib)%velFwake1)
                  enddo
                  call rotor(ir)%convectwake(iter, dt, 'P')
                else
                  ! Image lifting surface
                  call rotor(ir)%mirrorWake(rotor(rotor(ir)%imageRotorNum), 'P')
                endif
              endif
            endif
          enddo

          do ir = 1, nr
            if (rotor(ir)%nNwake > 0) then
              if (abs(rotor(ir)%surfaceType) == 1) then
                if (rotor(ir)%surfaceType .gt. 0) then
                  ! Lifting surface
                  do ib = 1, rotor(ir)%nbConvect
                    rotor(ir)%blade(ib)%velNwakePredicted(:, rotor(ir)%rowNear:rotor(ir)%nNwakeEnd, :) = 0._dp
                    rotor(ir)%blade(ib)%velFwakePredicted(:, rotor(ir)%rowFar:rotor(ir)%nFwakeEnd) = 0._dp
                    do jr = 1, nr
                      rotor(ir)%blade(ib)%velNwakePredicted(:, rotor(ir)%rowNear:rotor(ir)%nNwakeEnd, :) = &
                        rotor(ir)%blade(ib)%velNwakePredicted(:, rotor(ir)%rowNear:rotor(ir)%nNwakeEnd, :) + &
                        vind_onNwake_byRotor(rotor(jr), &
                        rotor(ir)%blade(ib)%waNPredicted(rotor(ir)%rowNear:rotor(ir)%nNwakeEnd, :), 'P')
                      rotor(ir)%blade(ib)%velFwakePredicted(:, rotor(ir)%rowFar:rotor(ir)%nFwakeEnd) = &
                        rotor(ir)%blade(ib)%velFwakePredicted(:, rotor(ir)%rowFar:rotor(ir)%nFwakeEnd) + &
                        vind_onFwake_byRotor(rotor(jr), &
                        rotor(ir)%blade(ib)%waFPredicted(rotor(ir)%rowFar:rotor(ir)%nFwakeEnd), 'P')
                    enddo
                    if (iter < switches%initWakeVelNt) then
                      do i = 1, 3
                        rotor(ir)%blade(ib)%velNwakePredicted(i, rotor(ir)%rowNear:rotor(ir)%nNwakeEnd, :) = &
                          rotor(ir)%blade(ib)%velNwakePredicted(i, rotor(ir)%rowNear:rotor(ir)%nNwakeEnd, :) - &
                          rotor(ir)%initWakeVel*rotor(ir)%shaftAxis(i)
                      enddo
                      do i = 1, 3
                        rotor(ir)%blade(ib)%velFwakePredicted(i, rotor(ir)%rowFar:rotor(ir)%nFwakeEnd) = &
                          rotor(ir)%blade(ib)%velFwakePredicted(i, rotor(ir)%rowFar:rotor(ir)%nFwakeEnd) - &
                          rotor(ir)%initWakeVel*rotor(ir)%shaftAxis(i)
                      enddo
                    endif
                  enddo
                endif
              endif
            endif
          enddo

          do ir = 1, nr
            if (rotor(ir)%nNwake > 0) then
              if (abs(rotor(ir)%surfaceType) == 1) then
                if (rotor(ir)%surfaceType .gt. 0) then
                  ! Lifting surface
                  do ib = 1, rotor(ir)%nbConvect
                    rotor(ir)%blade(ib)%velNwake = &
                      (rotor(ir)%blade(ib)%velNwakePredicted &
                      + rotor(ir)%blade(ib)%velNwakeStep)*0.5_dp
                    rotor(ir)%blade(ib)%velFwake = &
                      (rotor(ir)%blade(ib)%velFwakePredicted &
                      + rotor(ir)%blade(ib)%velFwakeStep)*0.5_dp
                  enddo
                  call rotor(ir)%convectwake(iter, dt, 'C')

                  do ib = 1, rotor(ir)%nbConvect
                    ! For next step
                    rotor(ir)%blade(ib)%velNwake1 = rotor(ir)%blade(ib)%velNwakeStep
                    rotor(ir)%blade(ib)%velFwake1 = rotor(ir)%blade(ib)%velFwakeStep
                  enddo
                else
                  ! Image lifting surface
                  call rotor(ir)%mirrorWake(rotor(rotor(ir)%imageRotorNum), 'C')
                endif
              endif
            endif
          enddo
        endif

      case (4)    ! Predictor-Corrector Adams-Moulton (3rd order)
        if (iter == 0) then
          do ir = 1, nr
            if (rotor(ir)%nNwake > 0) then
              if (abs(rotor(ir)%surfaceType) == 1) then
                if (rotor(ir)%surfaceType .gt. 0) then
                  call rotor(ir)%convectwake(iter, dt, 'C')
                  do ib = 1, rotor(ir)%nbConvect
                    rotor(ir)%blade(ib)%velNwake1 = rotor(ir)%blade(ib)%velNwake
                    rotor(ir)%blade(ib)%velFwake1 = rotor(ir)%blade(ib)%velFwake
                  enddo
                else
                  ! Image lifting surface
                  call rotor(ir)%mirrorWake(rotor(rotor(ir)%imageRotorNum), 'C')
                endif
              endif
            endif
          enddo
        elseif (iter == 2) then
          do ir = 1, nr
            if (rotor(ir)%nNwake > 0) then
              if (abs(rotor(ir)%surfaceType) == 1) then
                if (rotor(ir)%surfaceType .gt. 0) then
                  call rotor(ir)%convectwake(iter, dt, 'C')
                  do ib = 1, rotor(ir)%nbConvect
                    rotor(ir)%blade(ib)%velNwake2 = rotor(ir)%blade(ib)%velNwake
                    rotor(ir)%blade(ib)%velFwake2 = rotor(ir)%blade(ib)%velFwake
                  enddo
                else
                  ! Image lifting surface
                  call rotor(ir)%mirrorWake(rotor(rotor(ir)%imageRotorNum), 'C')
                endif
              endif
            endif
          enddo
        else
          do ir = 1, nr
            if (rotor(ir)%nNwake > 0) then
              if (abs(rotor(ir)%surfaceType) == 1) then
                if (rotor(ir)%surfaceType .gt. 0) then
                  do ib = 1, rotor(ir)%nbConvect
                    rotor(ir)%blade(ib)%waNPredicted(rotor(ir)%rowNear:rotor(ir)%nNwakeEnd, :) = &
                      rotor(ir)%blade(ib)%waN(rotor(ir)%rowNear:rotor(ir)%nNwakeEnd, :)
                    ! Store Nwake to Nwake_step for later use
                    rotor(ir)%blade(ib)%velNwakeStep = rotor(ir)%blade(ib)%velNwake
                    rotor(ir)%blade(ib)%velNwake = &
                      (23._dp*rotor(ir)%blade(ib)%velNwake &
                      - 16._dp*rotor(ir)%blade(ib)%velNwake2 &
                      + 05._dp*rotor(ir)%blade(ib)%velNwake1)/12._dp
                    rotor(ir)%blade(ib)%waFPredicted(rotor(ir)%rowFar:rotor(ir)%nFwakeEnd) = &
                      rotor(ir)%blade(ib)%waF(rotor(ir)%rowFar:rotor(ir)%nFwakeEnd)
                    ! Store Fwake to Fwake_step for later use
                    rotor(ir)%blade(ib)%velFwakeStep = rotor(ir)%blade(ib)%velFwake
                    rotor(ir)%blade(ib)%velFwake = &
                      (23._dp*rotor(ir)%blade(ib)%velFwake &
                      - 16._dp*rotor(ir)%blade(ib)%velFwake2 &
                      + 05._dp*rotor(ir)%blade(ib)%velFwake1)/12._dp
                  enddo
                  call rotor(ir)%convectwake(iter, dt, 'P')
                else
                  ! Image lifting surface
                  call rotor(ir)%mirrorWake(rotor(rotor(ir)%imageRotorNum), 'P')
                endif
              endif
            endif
          enddo

          do ir = 1, nr
            if (rotor(ir)%nNwake > 0) then
              if (abs(rotor(ir)%surfaceType) == 1) then
                if (rotor(ir)%surfaceType .gt. 0) then
                  do ib = 1, rotor(ir)%nbConvect
                    rotor(ir)%blade(ib)%velNwakePredicted(:, rotor(ir)%rowNear:rotor(ir)%nNwakeEnd, :) = 0._dp
                    rotor(ir)%blade(ib)%velFwakePredicted(:, rotor(ir)%rowFar:rotor(ir)%nFwakeEnd) = 0._dp
                    do jr = 1, nr
                      rotor(ir)%blade(ib)%velNwakePredicted(:, rotor(ir)%rowNear:rotor(ir)%nNwakeEnd, :) = &
                        rotor(ir)%blade(ib)%velNwakePredicted(:, rotor(ir)%rowNear:rotor(ir)%nNwakeEnd, :) + &
                        vind_onNwake_byRotor(rotor(jr), &
                        rotor(ir)%blade(ib)%waNPredicted(rotor(ir)%rowNear:rotor(ir)%nNwakeEnd, :), 'P')
                      rotor(ir)%blade(ib)%velFwakePredicted(:, rotor(ir)%rowFar:rotor(ir)%nFwakeEnd) = &
                        rotor(ir)%blade(ib)%velFwakePredicted(:, rotor(ir)%rowFar:rotor(ir)%nFwakeEnd) + &
                        vind_onFwake_byRotor(rotor(jr), &
                        rotor(ir)%blade(ib)%waFPredicted(rotor(ir)%rowFar:rotor(ir)%nFwakeEnd), 'P')
                    enddo
                    if (iter < switches%initWakeVelNt) then
                      do i = 1, 3
                        rotor(ir)%blade(ib)%velNwakePredicted(i, rotor(ir)%rowNear:rotor(ir)%nNwakeEnd, :) = &
                          rotor(ir)%blade(ib)%velNwakePredicted(i, rotor(ir)%rowNear:rotor(ir)%nNwakeEnd, :) - &
                          rotor(ir)%initWakeVel*rotor(ir)%shaftAxis(i)
                      enddo
                      do i = 1, 3
                        rotor(ir)%blade(ib)%velFwakePredicted(i, rotor(ir)%rowFar:rotor(ir)%nFwakeEnd) = &
                          rotor(ir)%blade(ib)%velFwakePredicted(i, rotor(ir)%rowFar:rotor(ir)%nFwakeEnd) - &
                          rotor(ir)%initWakeVel*rotor(ir)%shaftAxis(i)
                      enddo
                    endif
                  enddo
                endif
              endif
            endif
          enddo

          do ir = 1, nr
            if (rotor(ir)%nNwake > 0) then
              if (abs(rotor(ir)%surfaceType) == 1) then
                if (rotor(ir)%surfaceType .gt. 0) then
                  do ib = 1, rotor(ir)%nbConvect
                    rotor(ir)%blade(ib)%velNwake = &
                      (05._dp*rotor(ir)%blade(ib)%velNwakePredicted &
                      + 08._dp*rotor(ir)%blade(ib)%velNwakeStep &
                      - 01._dp*rotor(ir)%blade(ib)%velNwake2)/12._dp
                    rotor(ir)%blade(ib)%velFwake = &
                      (05._dp*rotor(ir)%blade(ib)%velFwakePredicted &
                      + 08._dp*rotor(ir)%blade(ib)%velFwakeStep &
                      - 01._dp*rotor(ir)%blade(ib)%velFwake2)/12._dp
                  enddo
                  call rotor(ir)%convectwake(iter, dt, 'C')
                  do ib = 1, rotor(ir)%nbConvect
                    rotor(ir)%blade(ib)%velNwake1 = rotor(ir)%blade(ib)%velNwake2
                    rotor(ir)%blade(ib)%velNwake2 = rotor(ir)%blade(ib)%velNwakeStep

                    rotor(ir)%blade(ib)%velFwake1 = rotor(ir)%blade(ib)%velFwake2
                    rotor(ir)%blade(ib)%velFwake2 = rotor(ir)%blade(ib)%velFwakeStep
                  enddo
                else
                  ! Image lifting surface
                  call rotor(ir)%mirrorWake(rotor(rotor(ir)%imageRotorNum), 'C')
                endif
              endif
            endif
          enddo
        endif

      case (5)    ! Predictor-Corrector Adams-Moulton (4th order)
        if (iter == 1) then
          do ir = 1, nr
            if (rotor(ir)%nNwake > 0) then
              if (abs(rotor(ir)%surfaceType) == 1) then
                if (rotor(ir)%surfaceType .gt. 0) then
                  call rotor(ir)%convectwake(iter, dt, 'C')
                  do ib = 1, rotor(ir)%nbConvect
                    rotor(ir)%blade(ib)%velNwake1 = rotor(ir)%blade(ib)%velNwake
                    rotor(ir)%blade(ib)%velFwake1 = rotor(ir)%blade(ib)%velFwake
                  enddo
                else
                  ! Image lifting surface
                  call rotor(ir)%mirrorWake(rotor(rotor(ir)%imageRotorNum), 'C')
                endif
              endif
            endif
          enddo
        elseif (iter == 2) then
          do ir = 1, nr
            if (rotor(ir)%nNwake > 0) then
              if (abs(rotor(ir)%surfaceType) == 1) then
                if (rotor(ir)%surfaceType .gt. 0) then
                  call rotor(ir)%convectwake(iter, dt, 'C')
                  do ib = 1, rotor(ir)%nbConvect
                    rotor(ir)%blade(ib)%velNwake2 = rotor(ir)%blade(ib)%velNwake
                    rotor(ir)%blade(ib)%velFwake2 = rotor(ir)%blade(ib)%velFwake
                  enddo
                else
                  ! Image lifting surface
                  call rotor(ir)%mirrorWake(rotor(rotor(ir)%imageRotorNum), 'C')
                endif
              endif
            endif
          enddo
        elseif (iter == 3) then
          do ir = 1, nr
            if (rotor(ir)%nNwake > 0) then
              if (abs(rotor(ir)%surfaceType) == 1) then
                if (rotor(ir)%surfaceType .gt. 0) then
                  call rotor(ir)%convectwake(iter, dt, 'C')
                  do ib = 1, rotor(ir)%nbConvect
                    rotor(ir)%blade(ib)%velNwake3 = rotor(ir)%blade(ib)%velNwake
                    rotor(ir)%blade(ib)%velFwake3 = rotor(ir)%blade(ib)%velFwake
                  enddo
                else
                  ! Image lifting surface
                  call rotor(ir)%mirrorWake(rotor(rotor(ir)%imageRotorNum), 'C')
                endif
              endif
            endif
          enddo
        else
          do ir = 1, nr
            if (rotor(ir)%nNwake > 0) then
              if (abs(rotor(ir)%surfaceType) == 1) then
                if (rotor(ir)%surfaceType .gt. 0) then
                  do ib = 1, rotor(ir)%nbConvect
                    rotor(ir)%blade(ib)%waNPredicted(rotor(ir)%rowNear:rotor(ir)%nNwakeEnd, :) = &
                      rotor(ir)%blade(ib)%waN(rotor(ir)%rowNear:rotor(ir)%nNwakeEnd, :)
                    ! Store Nwake to Nwake_step for later use
                    rotor(ir)%blade(ib)%velNwakeStep = rotor(ir)%blade(ib)%velNwake
                    rotor(ir)%blade(ib)%velNwake = &
                      & (55._dp*rotor(ir)%blade(ib)%velNwake &  ! Overwrite Nwake
                      & - 59._dp*rotor(ir)%blade(ib)%velNwake3 &
                      & + 37._dp*rotor(ir)%blade(ib)%velNwake2 &
                      & - 09._dp*rotor(ir)%blade(ib)%velNwake1)/24._dp
                    rotor(ir)%blade(ib)%waFPredicted(rotor(ir)%rowFar:rotor(ir)%nFwakeEnd) = &
                      rotor(ir)%blade(ib)%waF(rotor(ir)%rowFar:rotor(ir)%nFwakeEnd)
                    ! Store Fwake to Fwake_step for later use
                    rotor(ir)%blade(ib)%velFwakeStep = rotor(ir)%blade(ib)%velFwake
                    rotor(ir)%blade(ib)%velFwake = &
                      & (55._dp*rotor(ir)%blade(ib)%velFwake &  ! Overwrite Fwake
                      & - 59._dp*rotor(ir)%blade(ib)%velFwake3 &
                      & + 37._dp*rotor(ir)%blade(ib)%velFwake2 &
                      & - 09._dp*rotor(ir)%blade(ib)%velFwake1)/24._dp
                  enddo
                  call rotor(ir)%convectwake(iter, dt, 'P')
                else
                  ! Image lifting surface
                  call rotor(ir)%mirrorWake(rotor(rotor(ir)%imageRotorNum), 'P')
                endif
              endif
            endif
          enddo

          do ir = 1, nr
            if (rotor(ir)%nNwake > 0) then
              if (abs(rotor(ir)%surfaceType) == 1) then
                if (rotor(ir)%surfaceType .gt. 0) then
                  do ib = 1, rotor(ir)%nbConvect
                    rotor(ir)%blade(ib)%velNwakePredicted(:, rotor(ir)%rowNear:rotor(ir)%nNwakeEnd, :) = 0._dp
                    rotor(ir)%blade(ib)%velFwakePredicted(:, rotor(ir)%rowFar:rotor(ir)%nFwakeEnd) = 0._dp
                    do jr = 1, nr
                      rotor(ir)%blade(ib)%velNwakePredicted(:, rotor(ir)%rowNear:rotor(ir)%nNwakeEnd, :) = &
                        rotor(ir)%blade(ib)%velNwakePredicted(:, rotor(ir)%rowNear:rotor(ir)%nNwakeEnd, :) + &
                        vind_onNwake_byRotor(rotor(jr), &
                        rotor(ir)%blade(ib)%waNPredicted(rotor(ir)%rowNear:rotor(ir)%nNwakeEnd, :), 'P')
                      rotor(ir)%blade(ib)%velFwakePredicted(:, rotor(ir)%rowFar:rotor(ir)%nFwakeEnd) = &
                        rotor(ir)%blade(ib)%velFwakePredicted(:, rotor(ir)%rowFar:rotor(ir)%nFwakeEnd) + &
                        vind_onFwake_byRotor(rotor(jr), &
                        rotor(ir)%blade(ib)%waFPredicted(rotor(ir)%rowFar:rotor(ir)%nFwakeEnd), 'P')
                    enddo
                    if (iter < switches%initWakeVelNt) then
                      do i = 1, 3
                        rotor(ir)%blade(ib)%velNwakePredicted(i, rotor(ir)%rowNear:rotor(ir)%nNwakeEnd, :) = &
                          rotor(ir)%blade(ib)%velNwakePredicted(i, rotor(ir)%rowNear:rotor(ir)%nNwakeEnd, :) - &
                          rotor(ir)%initWakeVel*rotor(ir)%shaftAxis(i)
                      enddo
                      do i = 1, 3
                        rotor(ir)%blade(ib)%velFwakePredicted(i, rotor(ir)%rowFar:rotor(ir)%nFwakeEnd) = &
                          rotor(ir)%blade(ib)%velFwakePredicted(i, rotor(ir)%rowFar:rotor(ir)%nFwakeEnd) - &
                          rotor(ir)%initWakeVel*rotor(ir)%shaftAxis(i)
                      enddo
                    endif
                  enddo
                endif
              endif
            endif
          enddo

          do ir = 1, nr
            if (rotor(ir)%nNwake > 0) then
              if (abs(rotor(ir)%surfaceType) == 1) then
                if (rotor(ir)%surfaceType .gt. 0) then
                  do ib = 1, rotor(ir)%nbConvect
                    rotor(ir)%blade(ib)%velNwake = &
                      & (09._dp*rotor(ir)%blade(ib)%velNwakePredicted &
                      & + 19._dp*rotor(ir)%blade(ib)%velNwakeStep &
                      & - 05._dp*rotor(ir)%blade(ib)%velNwake3 &
                      & + 01._dp*rotor(ir)%blade(ib)%velNwake2)/24._dp
                    rotor(ir)%blade(ib)%velFwake = &
                      & (09._dp*rotor(ir)%blade(ib)%velFwakePredicted &
                      & + 19._dp*rotor(ir)%blade(ib)%velFwakeStep &
                      & - 05._dp*rotor(ir)%blade(ib)%velFwake3 &
                      & + 01._dp*rotor(ir)%blade(ib)%velFwake2)/24._dp
                  enddo
                  call rotor(ir)%convectwake(iter, dt, 'C')
                  do ib = 1, rotor(ir)%nbConvect
                    rotor(ir)%blade(ib)%velNwake1 = rotor(ir)%blade(ib)%velNwake2
                    rotor(ir)%blade(ib)%velNwake2 = rotor(ir)%blade(ib)%velNwake3
                    rotor(ir)%blade(ib)%velNwake3 = rotor(ir)%blade(ib)%velNwakeStep

                    rotor(ir)%blade(ib)%velFwake1 = rotor(ir)%blade(ib)%velFwake2
                    rotor(ir)%blade(ib)%velFwake2 = rotor(ir)%blade(ib)%velFwake3
                    rotor(ir)%blade(ib)%velFwake3 = rotor(ir)%blade(ib)%velFwakeStep
                  enddo
                else
                  ! Image lifting surface
                  call rotor(ir)%mirrorWake(rotor(rotor(ir)%imageRotorNum), 'C')
                endif
              endif
            endif
          enddo
        endif

      end select

      ! Strain wake
      if (switches%wakeStrain .eq. 1) then
        do ir = 1, nr
          if (rotor(ir)%nNwake > 0) then
            if (abs(rotor(ir)%surfaceType) == 1) &
              & call rotor(ir)%strain_wake()
          endif
        enddo
      endif

      ! Assign TE of wake and compute rollup
      do ir = 1, nr
        if (rotor(ir)%nNwake > 0) then
          if (abs(rotor(ir)%surfaceType) == 1) then
            if (rotor(ir)%surfaceType .gt. 0) then
              if (rotor(ir)%rowNear .eq. 1) then
                ! Last step of near wake or later steps
                call rotor(ir)%rollup()       ! Rollup wake for next far wake panel
                ! Store shed vortex as TE for next near wake panel
                call rotor(ir)%assignshed('TE')
              else
                if (rotor(ir)%surfaceType .gt. 0) then
                  call rotor(ir)%assignshed('TE')
                endif
              endif
            else
              ! Image lifting surface
              call rotor(ir)%mirrorWake(rotor(rotor(ir)%imageRotorNum), 'A')
            endif
          endif
        endif
      enddo
    endif

    ! Write out restart file in binary format
    if (switches%restartWriteNt .ne. 0) then
      if (mod(iter, switches%restartWriteNt) .eq. 0) then
        open(unit=24, file='Restart/restart'//timestamp//'.dat', &
          & status='replace', action='write', form='unformatted')
        write(24) t
        write(24) rotor
        close(24)
      endif
    endif
  enddo

  ! Deinitialize all variables
  do ir = 1, nr
    call rotor(ir)%deinit(switches)
  enddo

end program main
