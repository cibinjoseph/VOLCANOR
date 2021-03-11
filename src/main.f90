program main
  use libCommon
  use libPostprocess

  ! Ensure all necessary files exist
  inquire(file='config.in', exist=fileExists)
  if (.not. fileExists) error stop 'ERROR: config.in does not exist'

  ! Read config.in file
  print*
  call print_status('Reading file '//'config.in')
  call read_config('config.in')
  call print_status()    ! SUCCESS

  ! Allocate rotor objects
  allocate (rotor(nr))

  ! Read rotor??.in files
  do ir = 1, nr
    write (rotorChar, '(I0.2)') ir
    rotorFile = 'geom'//rotorChar//'.in'
    inquire(file=rotorFile, exist=fileExists)
    if (.not. fileExists) error stop 'ERROR: A rotorXX.in file does not exist'
    call print_status('Reading file '//rotorFile)
    call rotor(ir)%getdata(rotorFile)
    call print_status()    ! SUCCESS
  enddo

  ! Rotor and wake initialization
  do ir = 1, nr
    call rotor(ir)%init(ir, density, dt, nt, switches)
    call params2file(rotor(ir), ir, nt, dt, nr, density, velSound, switches)
  enddo

  ! Rotate wing pc, vr, cp and nCap by initial pitch angle
  do ir = 1, nr
    do ib = 1, rotor(ir)%nb
      rotor(ir)%blade(ib)%theta = rotor(ir)%gettheta(rotor(ir)%psiStart, ib)
      call rotor(ir)%blade(ib)%rot_pitch( &
        sign(1._dp, rotor(ir)%Omega)*rotor(ir)%blade(ib)%theta)
    enddo
  enddo

  ! Wing geometry plot for checks
  ! call wingverify(rotor(1)%blade(1)%wiP)

  ! Compute AIC and AIC_inv matrices for rotors
  do ir = 1, nr
    call print_status('Computing AIC matrix')
    call rotor(ir)%calcAIC()
    if (isInverse(rotor(ir)%AIC, rotor(ir)%AIC_inv)) then
      call print_status()    ! SUCCESS
    else
      ! if (rotor(ir)%surfaceType == 0) then
      stop 'Warning: Computed AIC_inv does not seem &
        & to be correct within given tolerance'
      ! endif
    endif
  enddo

  ! Initialize plot switches
  if (switches%wakePlot .lt. 0) &
    & switches%wakePlot = int(nt/abs(switches%wakePlot))
  if (switches%wakeTipPlot.lt. 0) &
    & switches%wakeTipPlot = int(nt/abs(switches%wakeTipPlot))
  if (switches%rotorForcePlot .lt. 0) &
    & switches%rotorForcePlot = int(nt/abs(switches%rotorForcePlot))
  if (switches%gridPlot .lt. 0) &
    & switches%gridPlot = int(nt/abs(switches%gridPlot))

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
  if (switches%slowStart .ne. 0) then
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
      do ib = 1, rotor(ir)%nb
        do is = 1, rotor(ir)%ns
          do ic = 1, rotor(ir)%nc
            row = ic + rotor(ir)%nc*(is - 1) &
              + rotor(ir)%ns*rotor(ir)%nc*(ib - 1)

            ! Translational vel
            rotor(ir)%blade(ib)%wiP(ic, is)%velCP = &
              -1._dp*rotor(ir)%velBody &
              ! Rotational vel
            - cross_product(rotor(ir)%omegaBody, &
              rotor(ir)%blade(ib)%wiP(ic, is)%CP - rotor(ir)%cgCoords) + &
              ! Omega vel
            cross_product(-rotor(ir)%omegaSlow*rotor(ir)%shaftAxis, &
              rotor(ir)%blade(ib)%wiP(ic, is)%CP - rotor(ir)%hubCoords)

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
      call rotor(ir)%dirLiftDrag()
      rotor(ir)%RHS = -1._dp*rotor(ir)%RHS
    enddo

    do ir = 1, nr
      rotor(ir)%gamVecPrev = rotor(ir)%gamVec
      rotor(ir)%gamVec = matmul(rotor(ir)%AIC_inv, rotor(ir)%RHS)

      ! Map gamVec to wing gam for each blade in rotor
      call rotor(ir)%map_gam()

      ! Check if initialization using converged soultion is requested
      if (switches%ntSubInit .ne. 0) then
        subIterResidual = norm2(rotor(ir)%gamVec - rotor(ir)%gamVecPrev)
        if ( subIterResidual .le. eps) &
          exit ntSubInitLoop
      endif
    enddo
  enddo ntSubInitLoop

  if ((i .gt. switches%ntSubInit) .and. (switches%ntSubInit .ne. 0)) then
    print *, 'Sub-iterations ', i, subIterResidual
    print*, "Initial solution did not converge. Try increasing sub-iterations."
  else
    call print_status()    ! SUCCESS
  endif

  do ir = 1, nr
    rotor(ir)%rowFar = rotor(ir)%nFwake + 1
    ! Since assignshed TE assigns to rowNear-1 panel
    rotor(ir)%rowNear = rotor(ir)%nNwake + 1
    call rotor(ir)%assignshed('TE')  ! Store shed vortex as TE
  enddo

  ! Compute forces
  if (switches%rotorForcePlot .ne. 0) then
    call init_plots(nr)    ! Create headers for plot files
    do ir = 1, nr
      ! Compute sec freestream velocity for secCL
      do ib = 1, rotor(ir)%nb
        do is = 1, rotor(ir)%ns
          rotor(ir)%blade(ib)%secVelFreestream(:, is) = &
            -1._dp*rotor(ir)%velBody &
            - cross_product(rotor(ir)%omegaBody, &
            rotor(ir)%blade(ib)%secCP(:, is) - rotor(ir)%cgCoords) &
            + cross_product(-rotor(ir)%omegaSlow*rotor(ir)%shaftAxis, &
            rotor(ir)%blade(ib)%secCP(:, is) - rotor(ir)%hubCoords)
        enddo
      enddo

      select case (rotor(ir)%forceCalcSwitch)

      case (0)  ! Compute using wing circulation
        call rotor(ir)%calc_force_gamma(density, dt)

        ! Compute and plot alpha if requested
        if (rotor(ir)%alphaPlotSwitch .ne. 0) then
          if (mod(iter, rotor(ir)%bladeforcePlotSwitch) .eq. 0) then
            ! Compute alpha
            do ib = 1, rotor(ir)%nb
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

            call rotor(ir)%calc_secAlpha()
          endif
        endif

      case (1)  ! Compute using alpha
        ! Compute alpha
        do ib = 1, rotor(ir)%nb
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

        call rotor(ir)%calc_secAlpha()
        call rotor(ir)%calc_force_alpha(density, velSound)

      case (2)  ! Compute lift using alpha approximated from sec circulation
        call rotor(ir)%calc_force_alphaGamma(density, velSound, dt)

      end select

      ! Initial force value
      call force2file(timestamp, rotor(ir), ir)

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
      rotor(ir)%rowNear = max(rotor(ir)%nNwake - (iter - 1), 1)
      rotor(ir)%rowFar = nt - (iter - 1)
      ! 0 => no roll up
      if (iter <= rotor(ir)%nNwake) rotor(ir)%rowFar = rotor(ir)%nFwake + 1
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
    enddo

    ! Assign LE of near wake
    do ir = 1, nr
      call rotor(ir)%assignshed('LE')  ! Store shed vortex as LE
    enddo

    ! Dissipate wake
    if (switches%wakeDissipation .eq. 1) then
      do ir = 1, nr
        ! Wake tip dissipation
        call rotor(ir)%dissipate_wake(dt)
      enddo
    endif

    ! Burst wake
    do ir = 1, nr
      if (switches%wakeBurst .ne. 0) then
        if (mod(iter, switches%wakeBurst) .eq. 0) &
          call rotor(ir)%burst_wake()
      endif
      ! Plot wake skew parameter
      call rotor(ir)%calc_skew()
      call skew2file(timestamp, rotor(ir), ir)
    enddo

    ! Write out wing n' wake
    do ir = 1, nr
      if (switches%wakePlot .ne. 0) then
        if (mod(iter, switches%wakePlot) .eq. 0) &
          call rotor2file(timestamp, rotor(ir), ir)
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
        do ib = 1, rotor(ir)%nb
          do is = 1, rotor(ir)%ns
            do ic = 1, rotor(ir)%nc
              row = ic + rotor(ir)%nc*(is - 1) &
                + rotor(ir)%ns*rotor(ir)%nc*(ib - 1)

              ! Translational vel
              rotor(ir)%blade(ib)%wiP(ic, is)%velCP = &
                -1._dp*rotor(ir)%velBody &
                ! Rotational vel
              - cross_product(rotor(ir)%omegaBody, &
                rotor(ir)%blade(ib)%wiP(ic, is)%CP - rotor(ir)%cgCoords) + &
                ! Omega vel
              cross_product(-rotor(ir)%omegaSlow*rotor(ir)%shaftAxis, &
                rotor(ir)%blade(ib)%wiP(ic, is)%CP - rotor(ir)%hubCoords)

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
        call rotor(ir)%dirLiftDrag()
        rotor(ir)%RHS = -1._dp*rotor(ir)%RHS
      enddo

      do ir = 1, nr
        rotor(ir)%gamvecPrev = rotor(ir)%gamVec
        rotor(ir)%gamVec = matmul(rotor(ir)%AIC_inv, rotor(ir)%RHS)

        ! Map gamVec to wing gam for each blade in rotor
        call rotor(ir)%map_gam()

        ! Check if sub-iterations are requested
        if (switches%ntSub .ne. 0) then
          subIterResidual = norm2(rotor(ir)%gamVec - rotor(ir)%gamVecPrev)
          if (subIterResidual .le. eps) &
            exit ntSubLoop
        endif
      enddo
    enddo ntSubLoop

    if (switches%ntSub .ne. 0) then
      print *, 'Sub-iterations ', i, subIterResidual
      if (i .gt. switches%ntSub) &
        print*, 'Solution did not converge. Try increasing sub-iterations.'
    endif

    ! Compute forces
    if (switches%rotorForcePlot .ne. 0) then
      if (mod(iter, switches%rotorForcePlot) .eq. 0) then
        do ir = 1, nr
          ! Compute sec freestream velocity for secCL
          do ib = 1, rotor(ir)%nb
            do is = 1, rotor(ir)%ns
              rotor(ir)%blade(ib)%secVelFreestream(:, is) = &
                -1._dp*rotor(ir)%velBody &
                - cross_product(rotor(ir)%omegaBody, &
                rotor(ir)%blade(ib)%secCP(:, is) - rotor(ir)%cgCoords) &
                + cross_product(-rotor(ir)%omegaSlow*rotor(ir)%shaftAxis, &
                rotor(ir)%blade(ib)%secCP(:, is) - rotor(ir)%hubCoords)
            enddo
          enddo

          select case (rotor(ir)%forceCalcSwitch)

          case (0)  ! Compute using wing circulation
            call rotor(ir)%calc_force_gamma(density, dt)

            ! Compute and plot alpha if requested
            if (rotor(ir)%alphaPlotSwitch .ne. 0) then
              if (mod(iter, rotor(ir)%bladeforcePlotSwitch) .eq. 0) then
                ! Compute alpha
                do ib = 1, rotor(ir)%nb
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

                call rotor(ir)%calc_secAlpha()
                ! call alpha2file(timestamp, rotor(ir), ir)
              endif
            endif

          case (1)  ! Compute using alpha
            ! Compute alpha
            do ib = 1, rotor(ir)%nb
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

            call rotor(ir)%calc_secAlpha()
            call rotor(ir)%calc_force_alpha(density, velSound)

          case (2)  ! Compute lift using alpha approximated from sec circulation
            call rotor(ir)%calc_force_alphaGamma(density, velSound, dt)

          end select

          call force2file(timestamp, rotor(ir), ir)

        enddo
      endif
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
          call gamma2file(timestamp, rotor(1), ir)
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
    ! Initialise wake velocity matrices
    do ir = 1, nr
      do ib = 1, rotor(ir)%nb
        rotor(ir)%blade(ib)%velNwake(:, rotor(ir)%rowNear:rotor(ir)%nNwake, :) = &
          0._dp
        rotor(ir)%blade(ib)%velFwake(:, rotor(ir)%rowFar:rotor(ir)%nFwake) = &
          0._dp
      enddo
    enddo

    ! Compute induced velocity due to rotors in domain
    do ir = 1, nr
      do ib = 1, rotor(ir)%nb
        do jr = 1, nr
          rotor(ir)%blade(ib)%velNwake(:, rotor(ir)%rowNear:rotor(ir)%nNwake, :) = &
            rotor(ir)%blade(ib)%velNwake(:, rotor(ir)%rowNear:rotor(ir)%nNwake, :) + &
            vind_onNwake_byRotor(rotor(jr), &
            rotor(ir)%blade(ib)%waP(rotor(ir)%rowNear:rotor(ir)%nNwake, :))
          rotor(ir)%blade(ib)%velFwake(:, rotor(ir)%rowFar:rotor(ir)%nFwake) = &
            rotor(ir)%blade(ib)%velFwake(:, rotor(ir)%rowFar:rotor(ir)%nFwake) + &
            vind_onFwake_byRotor(rotor(jr), &
            rotor(ir)%blade(ib)%waF(rotor(ir)%rowFar:rotor(ir)%nFwake))
        enddo

        ! Add initial wake velocity if provided
        if (iter < switches%initWakeVelNt) then
          do i = 1, 3
            rotor(ir)%blade(ib)%velNwake(i, rotor(ir)%rowNear:rotor(ir)%nNwake, :) = &
              rotor(ir)%blade(ib)%velNwake(i, rotor(ir)%rowNear:rotor(ir)%nNwake, :) + &
              rotor(ir)%initWakeVel*rotor(ir)%shaftAxis(i)
          enddo
          do i = 1, 3
            rotor(ir)%blade(ib)%velFwake(i, rotor(ir)%rowFar:rotor(ir)%nFwake) = &
              rotor(ir)%blade(ib)%velFwake(i, rotor(ir)%rowFar:rotor(ir)%nFwake) &
              - rotor(ir)%initWakeVel*rotor(ir)%shaftAxis(i)
          enddo
        endif
      enddo
    enddo

    ! Update wake vortex locations if lifting surface
    select case (switches%fdScheme)

    case (0)    ! Explicit Forward Diff (1st order)
      do ir = 1, nr
        if (rotor(ir)%surfaceType == 0) then
          do ib = 1, rotor(ir)%nb
            call rotor(ir)%blade(ib)%convectwake( &
              rotor(ir)%rowNear, rotor(ir)%rowFar, dt, 'C')
          enddo
        endif
      enddo

    case (1)    ! Predictor-Corrector (2nd order)
      ! Compute predicted wake
      do ir = 1, nr
        if (rotor(ir)%surfaceType == 0) then
        do ib = 1, rotor(ir)%nb
          rotor(ir)%blade(ib)%waPPredicted(rotor(ir)%rowNear:rotor(ir)%nNwake, :) = &
            rotor(ir)%blade(ib)%waP(rotor(ir)%rowNear:rotor(ir)%nNwake, :)
          rotor(ir)%blade(ib)%waFPredicted(rotor(ir)%rowFar:rotor(ir)%nFwake) = &
            rotor(ir)%blade(ib)%waF(rotor(ir)%rowFar:rotor(ir)%nFwake)
          call rotor(ir)%blade(ib)%convectwake(rotor(ir)%rowNear, &
            rotor(ir)%rowFar, dt, 'P')
        enddo
      endif
      enddo

      ! Compute velocity on predicted wake
      do ir = 1, nr
        if (rotor(ir)%surfaceType == 0) then
        do ib = 1, rotor(ir)%nb
          rotor(ir)%blade(ib)%velNwakePredicted(:, rotor(ir)%rowNear:rotor(ir)%nNwake, :) = 0._dp
          rotor(ir)%blade(ib)%velFwakePredicted(:, rotor(ir)%rowFar:rotor(ir)%nFwake) = 0._dp
          do jr = 1, nr
            rotor(ir)%blade(ib)%velNwakePredicted(:, rotor(ir)%rowNear:rotor(ir)%nNwake, :) = &
              rotor(ir)%blade(ib)%velNwakePredicted(:, rotor(ir)%rowNear:rotor(ir)%nNwake, :) + &
              vind_onNwake_byRotor(rotor(jr), &
              rotor(ir)%blade(ib)%waPPredicted(rotor(ir)%rowNear:rotor(ir)%nNwake, :), 'P')
            rotor(ir)%blade(ib)%velFwakePredicted(:, rotor(ir)%rowFar:rotor(ir)%nFwake) = &
              rotor(ir)%blade(ib)%velFwakePredicted(:, rotor(ir)%rowFar:rotor(ir)%nFwake) + &
              vind_onFwake_byRotor(rotor(jr), &
              rotor(ir)%blade(ib)%waFPredicted(rotor(ir)%rowFar:rotor(ir)%nFwake), 'P')
          enddo
          if (iter < switches%initWakeVelNt) then
            do i = 1, 3
              rotor(ir)%blade(ib)%velNwakePredicted(i, rotor(ir)%rowNear:rotor(ir)%nNwake, :) = &
                rotor(ir)%blade(ib)%velNwakePredicted(i, rotor(ir)%rowNear:rotor(ir)%nNwake, :) - &
                rotor(ir)%initWakeVel*rotor(ir)%shaftAxis(i)
            enddo
            do i = 1, 3
              rotor(ir)%blade(ib)%velFwakePredicted(i, rotor(ir)%rowFar:rotor(ir)%nFwake) = &
                rotor(ir)%blade(ib)%velFwakePredicted(i, rotor(ir)%rowFar:rotor(ir)%nFwake) - &
                rotor(ir)%initWakeVel*rotor(ir)%shaftAxis(i)
            enddo
          endif
        enddo
      endif
      enddo

      ! Compute averaged velocity and convect wake
      do ir = 1, nr
        if (rotor(ir)%surfaceType == 0) then
        do ib = 1, rotor(ir)%nb
          rotor(ir)%blade(ib)%velNwake(:, rotor(ir)%rowNear:rotor(ir)%nNwake, :) = &
            vel_order2_Nwake(rotor(ir)%blade(ib)%velNwake(:, &
            rotor(ir)%rowNear:rotor(ir)%nNwake, :), &
            rotor(ir)%blade(ib)%velNwakePredicted(:, rotor(ir)%rowNear:rotor(ir)%nNwake, :))
          !rotor(ir)%blade(ib)%velNwake(:,rotor(ir)%rowNear:rotor(ir)%nNwake,:)= &
          !0.5_dp*(rotor(ir)%blade(ib)%velNwake(:,rotor(ir)%rowNear:rotor(ir)%nNwake,:)+ &
          !rotor(ir)%blade(ib)%velNwakePredicted(:,rotor(ir)%rowNear:rotor(ir)%nNwake,:))
          rotor(ir)%blade(ib)%velFwake(:, rotor(ir)%rowFar:rotor(ir)%nFwake) = &
            vel_order2_Fwake(rotor(ir)%blade(ib)%velFwake(:, &
            rotor(ir)%rowFar:rotor(ir)%nFwake), &
            rotor(ir)%blade(ib)%velFwakePredicted(:, rotor(ir)%rowFar:rotor(ir)%nFwake))
          !rotor(ir)%blade(ib)%velFwake(:,rotor(ir)%rowFar:rotor(ir)%nFwake)= &
          !0.5_dp*(rotor(ir)%blade(ib)%velFwake(:,rotor(ir)%rowFar:rotor(ir)%nFwake)+ &
          !rotor(ir)%blade(ib)%velFwakePredicted(:,rotor(ir)%rowFar:rotor(ir)%nFwake))
          call rotor(ir)%blade(ib)%convectwake(rotor(ir)%rowNear, &
            rotor(ir)%rowFar, dt, 'C')
        enddo
      endif
    enddo

  case (2)    ! Explicit Adams-Bashforth (2nd order)
    if (iter == 1) then
      do ir = 1, nr
        if (rotor(ir)%surfaceType == 0) then
          do ib = 1, rotor(ir)%nb
            call rotor(ir)%blade(ib)%convectwake(rotor(ir)%rowNear, &
              rotor(ir)%rowFar, dt, 'C')
            rotor(ir)%blade(ib)%velNwake1 = rotor(ir)%blade(ib)%velNwake
            rotor(ir)%blade(ib)%velFwake1 = rotor(ir)%blade(ib)%velFwake
          enddo
        endif
      enddo
    else
      do ir = 1, nr
        if (rotor(ir)%surfaceType == 0) then
          do ib = 1, rotor(ir)%nb
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

            call rotor(ir)%blade(ib)%convectwake(rotor(ir)%rowNear, &
              rotor(ir)%rowFar, dt, 'C')
          enddo
        endif
      enddo
    endif

  case (3)    ! Predictor-Corrector Adams-Moulton (2nd order)
    if (iter == 1) then
        do ir = 1, nr
        if (rotor(ir)%surfaceType == 0) then
          do ib = 1, rotor(ir)%nb
            call rotor(ir)%blade(ib)%convectwake(rotor(ir)%rowNear, &
              rotor(ir)%rowFar, dt, 'C')
            rotor(ir)%blade(ib)%velNwake1 = rotor(ir)%blade(ib)%velNwake
            rotor(ir)%blade(ib)%velFwake1 = rotor(ir)%blade(ib)%velFwake
          enddo
        endif
        enddo
      else
        do ir = 1, nr
        if (rotor(ir)%surfaceType == 0) then
          do ib = 1, rotor(ir)%nb
            rotor(ir)%blade(ib)%waPPredicted(rotor(ir)%rowNear:rotor(ir)%nNwake, :) = &
              rotor(ir)%blade(ib)%waP(rotor(ir)%rowNear:rotor(ir)%nNwake, :)
            ! Store Nwake to Nwake_step for later use
            rotor(ir)%blade(ib)%velNwakeStep = rotor(ir)%blade(ib)%velNwake
            rotor(ir)%blade(ib)%velNwake = &
              0.5_dp*(3._dp*rotor(ir)%blade(ib)%velNwake - &
              rotor(ir)%blade(ib)%velNwake1)
            rotor(ir)%blade(ib)%waFPredicted(rotor(ir)%rowFar:rotor(ir)%nFwake) = &
              rotor(ir)%blade(ib)%waF(rotor(ir)%rowFar:rotor(ir)%nFwake)
            ! Store Fwake to Fwake_step for later use
            rotor(ir)%blade(ib)%velFwakeStep = rotor(ir)%blade(ib)%velFwake
            rotor(ir)%blade(ib)%velFwake = &
              0.5_dp*(3._dp*rotor(ir)%blade(ib)%velFwake - &
              rotor(ir)%blade(ib)%velFwake1)

            call rotor(ir)%blade(ib)%convectwake(rotor(ir)%rowNear, &
              rotor(ir)%rowFar, dt, 'P')
          enddo
        endif
        enddo

        do ir = 1, nr
          if (rotor(ir)%surfaceType == 0) then
            do ib = 1, rotor(ir)%nb
              rotor(ir)%blade(ib)%velNwakePredicted(:, rotor(ir)%rowNear:rotor(ir)%nNwake, :) = 0._dp
              rotor(ir)%blade(ib)%velFwakePredicted(:, rotor(ir)%rowFar:rotor(ir)%nFwake) = 0._dp
              do jr = 1, nr
                rotor(ir)%blade(ib)%velNwakePredicted(:, rotor(ir)%rowNear:rotor(ir)%nNwake, :) = &
                  rotor(ir)%blade(ib)%velNwakePredicted(:, rotor(ir)%rowNear:rotor(ir)%nNwake, :) + &
                  vind_onNwake_byRotor(rotor(jr), &
                  rotor(ir)%blade(ib)%waPPredicted(rotor(ir)%rowNear:rotor(ir)%nNwake, :), 'P')
                rotor(ir)%blade(ib)%velFwakePredicted(:, rotor(ir)%rowFar:rotor(ir)%nFwake) = &
                  rotor(ir)%blade(ib)%velFwakePredicted(:, rotor(ir)%rowFar:rotor(ir)%nFwake) + &
                  vind_onFwake_byRotor(rotor(jr), &
                  rotor(ir)%blade(ib)%waFPredicted(rotor(ir)%rowFar:rotor(ir)%nFwake), 'P')
              enddo
              if (iter < switches%initWakeVelNt) then
                do i = 1, 3
                  rotor(ir)%blade(ib)%velNwakePredicted(i, rotor(ir)%rowNear:rotor(ir)%nNwake, :) = &
                    rotor(ir)%blade(ib)%velNwakePredicted(i, rotor(ir)%rowNear:rotor(ir)%nNwake, :) - &
                    rotor(ir)%initWakeVel*rotor(ir)%shaftAxis(i)
                enddo
                do i = 1, 3
                  rotor(ir)%blade(ib)%velFwakePredicted(i, rotor(ir)%rowFar:rotor(ir)%nFwake) = &
                    rotor(ir)%blade(ib)%velFwakePredicted(i, rotor(ir)%rowFar:rotor(ir)%nFwake) - &
                    rotor(ir)%initWakeVel*rotor(ir)%shaftAxis(i)
                enddo
              endif
            enddo
          endif
        enddo

        do ir = 1, nr
          if (rotor(ir)%surfaceType == 0) then
            do ib = 1, rotor(ir)%nb
              rotor(ir)%blade(ib)%velNwake = &
                (rotor(ir)%blade(ib)%velNwakePredicted &
                + rotor(ir)%blade(ib)%velNwakeStep)*0.5_dp
              rotor(ir)%blade(ib)%velFwake = &
                (rotor(ir)%blade(ib)%velFwakePredicted &
                + rotor(ir)%blade(ib)%velFwakeStep)*0.5_dp
              call rotor(ir)%blade(ib)%convectwake(rotor(ir)%rowNear, &
                rotor(ir)%rowFar, dt, 'C')

              ! For next step
              rotor(ir)%blade(ib)%velNwake1 = rotor(ir)%blade(ib)%velNwakeStep
              rotor(ir)%blade(ib)%velFwake1 = rotor(ir)%blade(ib)%velFwakeStep
            enddo
          endif
        enddo
      endif

    case (4)    ! Predictor-Corrector Adams-Moulton (3rd order)
      if (iter == 1) then
        do ir = 1, nr
          if (rotor(ir)%surfaceType == 0) then
            do ib = 1, rotor(ir)%nb
              call rotor(ir)%blade(ib)%convectwake(rotor(ir)%rowNear, &
                rotor(ir)%rowFar, dt, 'C')
              rotor(ir)%blade(ib)%velNwake1 = rotor(ir)%blade(ib)%velNwake
              rotor(ir)%blade(ib)%velFwake1 = rotor(ir)%blade(ib)%velFwake
            enddo
          endif
        enddo
      elseif (iter == 2) then
        do ir = 1, nr
          if (rotor(ir)%surfaceType == 0) then
            do ib = 1, rotor(ir)%nb
              call rotor(ir)%blade(ib)%convectwake(rotor(ir)%rowNear, &
                rotor(ir)%rowFar, dt, 'C')
              rotor(ir)%blade(ib)%velNwake2 = rotor(ir)%blade(ib)%velNwake
              rotor(ir)%blade(ib)%velFwake2 = rotor(ir)%blade(ib)%velFwake
            enddo
          endif
        enddo
      else
        do ir = 1, nr
          if (rotor(ir)%surfaceType == 0) then
            do ib = 1, rotor(ir)%nb
              rotor(ir)%blade(ib)%waPPredicted(rotor(ir)%rowNear:rotor(ir)%nNwake, :) = &
                rotor(ir)%blade(ib)%waP(rotor(ir)%rowNear:rotor(ir)%nNwake, :)
              ! Store Nwake to Nwake_step for later use
              rotor(ir)%blade(ib)%velNwakeStep = rotor(ir)%blade(ib)%velNwake
              rotor(ir)%blade(ib)%velNwake = &
                (23._dp*rotor(ir)%blade(ib)%velNwake &
                - 16._dp*rotor(ir)%blade(ib)%velNwake2 &
                + 05._dp*rotor(ir)%blade(ib)%velNwake1)/12._dp
              rotor(ir)%blade(ib)%waFPredicted(rotor(ir)%rowFar:rotor(ir)%nFwake) = &
                rotor(ir)%blade(ib)%waF(rotor(ir)%rowFar:rotor(ir)%nFwake)
              ! Store Fwake to Fwake_step for later use
              rotor(ir)%blade(ib)%velFwakeStep = rotor(ir)%blade(ib)%velFwake
              rotor(ir)%blade(ib)%velFwake = &
                (23._dp*rotor(ir)%blade(ib)%velFwake &
                - 16._dp*rotor(ir)%blade(ib)%velFwake2 &
                + 05._dp*rotor(ir)%blade(ib)%velFwake1)/12._dp
              call rotor(ir)%blade(ib)%convectwake(rotor(ir)%rowNear, &
                rotor(ir)%rowFar, dt, 'P')
            enddo
          endif
        enddo

        do ir = 1, nr
          if (rotor(ir)%surfaceType == 0) then
            do ib = 1, rotor(ir)%nb
              rotor(ir)%blade(ib)%velNwakePredicted(:, rotor(ir)%rowNear:rotor(ir)%nNwake, :) = 0._dp
              rotor(ir)%blade(ib)%velFwakePredicted(:, rotor(ir)%rowFar:rotor(ir)%nFwake) = 0._dp
              do jr = 1, nr
                rotor(ir)%blade(ib)%velNwakePredicted(:, rotor(ir)%rowNear:rotor(ir)%nNwake, :) = &
                  rotor(ir)%blade(ib)%velNwakePredicted(:, rotor(ir)%rowNear:rotor(ir)%nNwake, :) + &
                  vind_onNwake_byRotor(rotor(jr), &
                  rotor(ir)%blade(ib)%waPPredicted(rotor(ir)%rowNear:rotor(ir)%nNwake, :), 'P')
                rotor(ir)%blade(ib)%velFwakePredicted(:, rotor(ir)%rowFar:rotor(ir)%nFwake) = &
                  rotor(ir)%blade(ib)%velFwakePredicted(:, rotor(ir)%rowFar:rotor(ir)%nFwake) + &
                  vind_onFwake_byRotor(rotor(jr), &
                  rotor(ir)%blade(ib)%waFPredicted(rotor(ir)%rowFar:rotor(ir)%nFwake), 'P')
              enddo
              if (iter < switches%initWakeVelNt) then
                do i = 1, 3
                  rotor(ir)%blade(ib)%velNwakePredicted(i, rotor(ir)%rowNear:rotor(ir)%nNwake, :) = &
                    rotor(ir)%blade(ib)%velNwakePredicted(i, rotor(ir)%rowNear:rotor(ir)%nNwake, :) - &
                    rotor(ir)%initWakeVel*rotor(ir)%shaftAxis(i)
                enddo
                do i = 1, 3
                  rotor(ir)%blade(ib)%velFwakePredicted(i, rotor(ir)%rowFar:rotor(ir)%nFwake) = &
                    rotor(ir)%blade(ib)%velFwakePredicted(i, rotor(ir)%rowFar:rotor(ir)%nFwake) - &
                    rotor(ir)%initWakeVel*rotor(ir)%shaftAxis(i)
                enddo
              endif
            enddo
          endif
        enddo

        do ir = 1, nr
          if (rotor(ir)%surfaceType == 0) then
            do ib = 1, rotor(ir)%nb
              rotor(ir)%blade(ib)%velNwake = &
                (05._dp*rotor(ir)%blade(ib)%velNwakePredicted &
                + 08._dp*rotor(ir)%blade(ib)%velNwakeStep &
                - 01._dp*rotor(ir)%blade(ib)%velNwake2)/12._dp
              rotor(ir)%blade(ib)%velFwake = &
                (05._dp*rotor(ir)%blade(ib)%velFwakePredicted &
                + 08._dp*rotor(ir)%blade(ib)%velFwakeStep &
                - 01._dp*rotor(ir)%blade(ib)%velFwake2)/12._dp
              call rotor(ir)%blade(ib)%convectwake(rotor(ir)%rowNear, &
                rotor(ir)%rowFar, dt, 'C')

              rotor(ir)%blade(ib)%velNwake1 = rotor(ir)%blade(ib)%velNwake2
              rotor(ir)%blade(ib)%velNwake2 = rotor(ir)%blade(ib)%velNwakeStep

              rotor(ir)%blade(ib)%velFwake1 = rotor(ir)%blade(ib)%velFwake2
              rotor(ir)%blade(ib)%velFwake2 = rotor(ir)%blade(ib)%velFwakeStep
            enddo
          endif
        enddo
      endif

    case (5)    ! Predictor-Corrector Adams-Moulton (4th order)
      if (iter == 1) then
        do ir = 1, nr
          if (rotor(ir)%surfaceType == 0) then
            do ib = 1, rotor(ir)%nb
              call rotor(ir)%blade(ib)%convectwake(rotor(ir)%rowNear, &
                rotor(ir)%rowFar, dt, 'C')
              rotor(ir)%blade(ib)%velNwake1 = rotor(ir)%blade(ib)%velNwake
              rotor(ir)%blade(ib)%velFwake1 = rotor(ir)%blade(ib)%velFwake
            enddo
          endif
        enddo
      elseif (iter == 2) then
        do ir = 1, nr
          if (rotor(ir)%surfaceType == 0) then
            do ib = 1, rotor(ir)%nb
              call rotor(ir)%blade(ib)%convectwake(rotor(ir)%rowNear, &
                rotor(ir)%rowFar, dt, 'C')
              rotor(ir)%blade(ib)%velNwake2 = rotor(ir)%blade(ib)%velNwake
              rotor(ir)%blade(ib)%velFwake2 = rotor(ir)%blade(ib)%velFwake
            enddo
          endif
        enddo
      elseif (iter == 3) then
        do ir = 1, nr
          if (rotor(ir)%surfaceType == 0) then
            do ib = 1, rotor(ir)%nb
              call rotor(ir)%blade(ib)%convectwake(rotor(ir)%rowNear, &
                rotor(ir)%rowFar, dt, 'C')
              rotor(ir)%blade(ib)%velNwake3 = rotor(ir)%blade(ib)%velNwake
              rotor(ir)%blade(ib)%velFwake3 = rotor(ir)%blade(ib)%velFwake
            enddo
          endif
        enddo
      else
        do ir = 1, nr
          if (rotor(ir)%surfaceType == 0) then
            do ib = 1, rotor(ir)%nb
              rotor(ir)%blade(ib)%waPPredicted(rotor(ir)%rowNear:rotor(ir)%nNwake, :) = &
                rotor(ir)%blade(ib)%waP(rotor(ir)%rowNear:rotor(ir)%nNwake, :)
              ! Store Nwake to Nwake_step for later use
              rotor(ir)%blade(ib)%velNwakeStep = rotor(ir)%blade(ib)%velNwake
              rotor(ir)%blade(ib)%velNwake = &
                55._dp/24._dp*rotor(ir)%blade(ib)%velNwake &  ! Overwrite Nwake
                - 59._dp/24._dp*rotor(ir)%blade(ib)%velNwake3 &
                + 37._dp/24._dp*rotor(ir)%blade(ib)%velNwake2 &
                - 09._dp/24._dp*rotor(ir)%blade(ib)%velNwake1
              rotor(ir)%blade(ib)%waFPredicted(rotor(ir)%rowFar:rotor(ir)%nFwake) = &
                rotor(ir)%blade(ib)%waF(rotor(ir)%rowFar:rotor(ir)%nFwake)
              ! Store Fwake to Fwake_step for later use
              rotor(ir)%blade(ib)%velFwakeStep = rotor(ir)%blade(ib)%velFwake
              rotor(ir)%blade(ib)%velFwake = &
                55._dp/24._dp*rotor(ir)%blade(ib)%velFwake &  ! Overwrite Fwake
                - 59._dp/24._dp*rotor(ir)%blade(ib)%velFwake3 &
                + 37._dp/24._dp*rotor(ir)%blade(ib)%velFwake2 &
                - 09._dp/24._dp*rotor(ir)%blade(ib)%velFwake1
              call rotor(ir)%blade(ib)%convectwake(rotor(ir)%rowNear, &
                rotor(ir)%rowFar, dt, 'P')
            enddo
          endif
        enddo

        do ir = 1, nr
          if (rotor(ir)%surfaceType == 0) then
            do ib = 1, rotor(ir)%nb
              rotor(ir)%blade(ib)%velNwakePredicted(:, rotor(ir)%rowNear:rotor(ir)%nNwake, :) = 0._dp
              rotor(ir)%blade(ib)%velFwakePredicted(:, rotor(ir)%rowFar:rotor(ir)%nFwake) = 0._dp
              do jr = 1, nr
                rotor(ir)%blade(ib)%velNwakePredicted(:, rotor(ir)%rowNear:rotor(ir)%nNwake, :) = &
                  rotor(ir)%blade(ib)%velNwakePredicted(:, rotor(ir)%rowNear:rotor(ir)%nNwake, :) + &
                  vind_onNwake_byRotor(rotor(jr), &
                  rotor(ir)%blade(ib)%waPPredicted(rotor(ir)%rowNear:rotor(ir)%nNwake, :), 'P')
                rotor(ir)%blade(ib)%velFwakePredicted(:, rotor(ir)%rowFar:rotor(ir)%nFwake) = &
                  rotor(ir)%blade(ib)%velFwakePredicted(:, rotor(ir)%rowFar:rotor(ir)%nFwake) + &
                  vind_onFwake_byRotor(rotor(jr), &
                  rotor(ir)%blade(ib)%waFPredicted(rotor(ir)%rowFar:rotor(ir)%nFwake), 'P')
              enddo
              if (iter < switches%initWakeVelNt) then
                do i = 1, 3
                  rotor(ir)%blade(ib)%velNwakePredicted(i, rotor(ir)%rowNear:rotor(ir)%nNwake, :) = &
                    rotor(ir)%blade(ib)%velNwakePredicted(i, rotor(ir)%rowNear:rotor(ir)%nNwake, :) - &
                    rotor(ir)%initWakeVel*rotor(ir)%shaftAxis(i)
                enddo
                do i = 1, 3
                  rotor(ir)%blade(ib)%velFwakePredicted(i, rotor(ir)%rowFar:rotor(ir)%nFwake) = &
                    rotor(ir)%blade(ib)%velFwakePredicted(i, rotor(ir)%rowFar:rotor(ir)%nFwake) - &
                    rotor(ir)%initWakeVel*rotor(ir)%shaftAxis(i)
                enddo
              endif
            enddo
          endif
        enddo

        do ir = 1, nr
          if (rotor(ir)%surfaceType == 0) then
            do ib = 1, rotor(ir)%nb
              rotor(ir)%blade(ib)%velNwake = &
                09._dp/24._dp*rotor(ir)%blade(ib)%velNwakePredicted &
                + 19._dp/24._dp*rotor(ir)%blade(ib)%velNwakeStep &
                - 05._dp/24._dp*rotor(ir)%blade(ib)%velNwake3 &
                + 01._dp/24._dp*rotor(ir)%blade(ib)%velNwake2
              rotor(ir)%blade(ib)%velFwake = &
                09._dp/24._dp*rotor(ir)%blade(ib)%velFwakePredicted &
                + 19._dp/24._dp*rotor(ir)%blade(ib)%velFwakeStep &
                - 05._dp/24._dp*rotor(ir)%blade(ib)%velFwake3 &
                + 01._dp/24._dp*rotor(ir)%blade(ib)%velFwake2
              call rotor(ir)%blade(ib)%convectwake(rotor(ir)%rowNear, &
                rotor(ir)%rowFar, dt, 'C')

              rotor(ir)%blade(ib)%velNwake1 = rotor(ir)%blade(ib)%velNwake2
              rotor(ir)%blade(ib)%velNwake2 = rotor(ir)%blade(ib)%velNwake3
              rotor(ir)%blade(ib)%velNwake3 = rotor(ir)%blade(ib)%velNwakeStep

              rotor(ir)%blade(ib)%velFwake1 = rotor(ir)%blade(ib)%velFwake2
              rotor(ir)%blade(ib)%velFwake2 = rotor(ir)%blade(ib)%velFwake3
              rotor(ir)%blade(ib)%velFwake3 = rotor(ir)%blade(ib)%velFwakeStep
            enddo
          endif
        enddo
      endif

    end select

    ! Strain wake
    if (switches%wakeStrain .eq. 1) then
      do ir = 1, nr
        if (rotor(ir)%surfaceType == 0) &
          & call rotor(ir)%strain_wake()
      enddo
    endif

    ! Assign TE of wake and compute rollup
    do ir = 1, nr
      if (rotor(ir)%surfaceType == 0) then
        if ((rotor(ir)%rowNear .eq. 1) .and. (rotor(ir)%rowFar /= 1)) then
          ! Last step of near wake or later steps
          call rotor(ir)%rollup()       ! Rollup wake for next far wake panel
          call rotor(ir)%shiftwake()    ! Shift wake
          ! Store shed vortex as TE for next near wake panel
          call rotor(ir)%assignshed('TE')
        else
          call rotor(ir)%assignshed('TE')
        endif
      endif
    enddo

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
