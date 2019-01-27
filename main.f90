program main
  use library
  use postproc

  ! Initialize variables
  include "init_file.f90"
  print*

  ! Read config.in file
  call print_status('Reading file '//'config.in')
  open(unit=11,file='config.in')
  call skiplines(11,2)
  read(11,*) nt,dt,nr
  call skiplines(11,4)
  read(11,*) density, turbulentViscosity
  call skiplines(11,4)
  read(11,*) spanSpacingSwitch
  call skiplines(11,4)
  read(11,*) tipDissipationSwitch, wakeStrainSwitch
  call skiplines(11,4)
  read(11,*) slowStartSwitch, slowStartNt
  call skiplines(11,4)
  read(11,*) wakePlotSwitch, gridPlotSwitch, forcePlotSwitch
  call skiplines(11,4)
  read(11,*) fdSchemeSwitch
  call skiplines(11,4)
  read(11,*) wakeIgnoreNt
  call skiplines(11,4)
  read(11,*) initWakeVelNt
  close(11)
  call print_status()    ! SUCCESS

  ! Allocate rotor objects
  allocate(rotor(nr))

  ! Read rotor??.in files
  do ir=1,nr
    write(rotorChar,'(I0.2)') ir
    rotorFile='rotor'//rotorChar//'.in'
    call print_status('Reading file '//rotorFile)
    call rotor(ir)%getdata(rotorFile,nt)
    call print_status()    ! SUCCESS
  enddo

  ! Rotor and wake initialization
  do ir=1,nr
    call rotor(ir)%init(density,dt,spanSpacingSwitch,fdSchemeSwitch)
  enddo

  ! Rotate wing pc, vr, cp and nCap by initial pitch angle 
  do ir=1,nr
    do ib=1,rotor(ir)%nb
      rotor(ir)%blade(ib)%theta=rotor(ir)%gettheta(rotor(ir)%psiStart,ib)
      call rotor(ir)%blade(ib)%rot_pitch(sign(1._dp,rotor(ir)%Omega)*rotor(ir)%blade(ib)%theta)
    enddo
  enddo

  ! Compute AIC and AIC_inv matrices for rotors
  call print_status('Computing AIC matrices')
  do ir=1,nr
    call rotor(ir)%calcAIC()
    if (isInverse(rotor(ir)%AIC,rotor(ir)%AIC_inv)) then
      call print_status()    ! SUCCESS
    else
      stop 'Warning: Computed AIC_inv does not seem to be correct within given tolerance'
      read*
    endif
  enddo

  ! Obtain initial solution
  call print_status('Computing initial solution')
  if (slowStartSwitch .ne. 0) then
    do ir=1,nr
      rotor(ir)%omegaSlow=0._dp
    enddo
  else 
    do ir=1,nr
      rotor(ir)%omegaSlow=rotor(ir)%Omega
    enddo
  endif

  t=0._dp
  iter=0
  write(timestamp,'(I0.5)') iter

  ! Compute RHS
  do ir=1,nr
    do ib=1,rotor(ir)%nb
      do is=1,rotor(ir)%ns
        do ic=1,rotor(ir)%nc
          row=ic+rotor(ir)%nc*(is-1)+rotor(ir)%ns*rotor(ir)%nc*(ib-1)

          ! Rotational vel
          rotor(ir)%RHS(row) = dot_product(rotor(ir)%velWind,rotor(ir)%blade(ib)%wiP(ic,is)%nCap)

          ! Rotational vel
          rotor(ir)%RHS(row)=rotor(ir)%RHS(row)+dot_product(cross3(rotor(ir)%omegaWind  &
            ,rotor(ir)%blade(ib)%wiP(ic,is)%CP-rotor(ir)%cgCoords),rotor(ir)%blade(ib)%wiP(ic,is)%nCap)

          ! Omega vel
          rotor(ir)%RHS(row)=rotor(ir)%RHS(row)+dot_product(cross3(-rotor(ir)%omegaSlow*rotor(ir)%shaftAxis  &
            ,rotor(ir)%blade(ib)%wiP(ic,is)%CP-rotor(ir)%hubCoords),rotor(ir)%blade(ib)%wiP(ic,is)%nCap)

          ! Pitch vel
          !rotor(ir)%blade%(ib)%wing(ic,is)%velPitch=rotor(ir)%thetadot_pitch(0._dp,ib)*rotor(ir)%blade(ib)%wiP(ic,is)%rHinge
          !rotor(ir)%RHS(row)= RHS(row)+wing(ib,ic,is)%velPitch
        enddo
      enddo
    enddo
    rotor(ir)%RHS=-1._dp*rotor(ir)%RHS
    rotor(ir)%gamVec=matmul(rotor(ir)%AIC_inv,rotor(ir)%RHS)
  enddo
  call print_status()    ! SUCCESS

  ! Map gamVec to wing gam for each blade in rotor
  do ir=1,nr
    call rotor(ir)%map_gam()
  enddo

  do ir=1,nr
    rotor(ir)%rowFar=0
    rotor(ir)%rowNear=rotor(ir)%nNwake+1  ! Since assignshed TE assigns to rowNear-1 panel
    call rotor(ir)%assignshed('TE')  ! Store shed vortex as TE
  enddo

  ! Compute forces
  if (forcePlotSwitch .ne. 0) then
    call init_plots(nr)    ! Create headers for plot files
    do ir=1,nr
      call rotor(ir)%calc_force(density,dt)
      call force2file(timestamp,rotor(ir),ir,-zAxis)  ! Negative sign due to negative inflow or gamma

      ! Compute alpha
      do ib=1,rotor(ir)%nb
        do is=1,rotor(ir)%ns
          do ic=1,rotor(ir)%nc
            ! Compute local velocity vector (excluding induced velocities from wing bound vortices)
            do jr=1,nr
              rotor(ir)%blade(ib)%wiP(ic,is)%velLocal=rotor(ir)%blade(ib)%wiP(ic,is)%velLocal+  &
                rotor(ir)%blade(ib)%wiP(ic,is)%velCP+  &
                rotor(jr)%vind_bywing(rotor(ir)%blade(ib)%wiP(ic,is)%CP)-  &
                rotor(jr)%vind_bywing_boundVortices(rotor(ir)%blade(ib)%wiP(ic,is)%CP)
            enddo

            rotor(ir)%blade(ib)%wiP(ic,is)%alpha=  &
              dot_product(rotor(ir)%blade(ib)%wiP(ic,is)%velLocal,  &
              rotor(ir)%blade(ib)%wiP(ic,is)%tauCapChord)/  &
              norm2(rotor(ir)%blade(ib)%wiP(ic,is)%velLocal)
          enddo
        enddo
      enddo
    enddo
  endif

  ! ------- MAIN LOOP START -------
  do iter=1,nt
    t=t+dt
    print*,iter,nt
    write(timestamp,'(I0.5)') iter
    do ir=1,nr
      rotor(ir)%rowNear=max(rotor(ir)%nNwake-(iter-1),1)
      rotor(ir)%rowFar=nt-(iter-1)
      if (iter<=rotor(ir)%nNwake) rotor(ir)%rowFar=0    ! 0 implies roll up has not occured
    enddo

    select case (slowStartSwitch)
    case (0)    ! No slow start
      do ir=1,nr
        rotor(ir)%omegaSlow=rotor(ir)%Omega
      enddo
    case (1)    ! Linear slope
      do ir=1,nr
        rotor(ir)%omegaSlow=min(real(slowStartNt),real(iter+1))*rotor(ir)%Omega/slowStartNt
      enddo
    case (2)    ! tanh slope
      do ir=1,nr
        rotor(ir)%omegaSlow=tanh(5._dp*iter/slowStartNt)*rotor(ir)%Omega
      enddo
    case (3)    ! tanh slope
      do ir=1,nr
        rotor(ir)%omegaSlow=(tanh(6._dp*real(iter)/real(slowStartNt)-3._dp)+1._dp)*0.5_dp*rotor(ir)%Omega
      enddo
    case default
      error stop "Assign correct slowStartSwitch"
    end select

    if (tipDissipationSwitch .eq. 1) then
      do ir=1,nr
        if (rotor(ir)%rowFar .ne. 0) then
          ! Age vortex filaments
          call rotor(ir)%age_wake(dt)

          ! Wake tip dissipation
          call rotor(ir)%dissipate_tip(turbulentViscosity)
        endif
      enddo
    endif

    ! Move wing
    do ir=1,nr
      call rotor(ir)%move(rotor(ir)%velBody*dt)
      call rotor(ir)%rot_pts(rotor(ir)%omegaBody*dt,rotor(ir)%cgCoords,1)
      call rotor(ir)%rot_advance(rotor(ir)%omegaSlow*dt)
    enddo

    ! Assign LE of near wake
    do ir=1,nr
      call rotor(ir)%assignshed('LE')  ! Store shed vortex as LE
    enddo

    ! Write out wing n' wake
    do ir=1,nr
      if (wakePlotSwitch .ne. 0) then
        if ((mod(iter,wakePlotSwitch) .eq. 0) .and.(rotor(ir)%rowFar .ne. 0) ) call rotor2file(timestamp,rotor(ir))
      endif
    enddo

    ! Compute RHS
    do ir=1,nr
      rotor(ir)%RHS=0._dp
      do ib=1,rotor(ir)%nb
        do is=1,rotor(ir)%ns
          do ic=1,rotor(ir)%nc
            row=ic+rotor(ir)%nc*(is-1)+rotor(ir)%ns*rotor(ir)%nc*(ib-1)

            ! Translational vel
            rotor(ir)%blade(ib)%wiP(ic,is)%velCP=rotor(ir)%velWind

            ! Rotational vel
            rotor(ir)%blade(ib)%wiP(ic,is)%velCP=rotor(ir)%blade(ib)%wiP(ic,is)%velCP  &
              +cross3(rotor(ir)%omegaWind,rotor(ir)%blade(ib)%wiP(ic,is)%CP-rotor(ir)%cgCoords)

            ! Omega vel
            rotor(ir)%blade(ib)%wiP(ic,is)%velCP=rotor(ir)%blade(ib)%wiP(ic,is)%velCP  &
              +cross3(-rotor(ir)%omegaSlow*rotor(ir)%shaftAxis,rotor(ir)%blade(ib)%wiP(ic,is)%CP-rotor(ir)%hubCoords)

            ! Wake vel
            do jr=1,nr
              rotor(ir)%blade(ib)%wiP(ic,is)%velCP=rotor(ir)%blade(ib)%wiP(ic,is)%velCP  &
                +rotor(jr)%vind_bywake(rotor(ir)%blade(ib)%wiP(ic,is)%CP)
            enddo

            rotor(ir)%RHS(row)=dot_product(rotor(ir)%blade(ib)%wiP(ic,is)%velCP,rotor(ir)%blade(ib)%wiP(ic,is)%nCap)

            ! Pitch vel
            !wing(ib,ic,is)%velPitch=thetadot*wing(ib,ic,is)%rHinge
            !RHS(row)=RHS(row)+wing(ib,ic,is)%velPitch
          enddo
        enddo
      enddo
      rotor(ir)%RHS=-1._dp*rotor(ir)%RHS
    enddo

    do ir=1,nr
      call rotor(ir)%record_gamPrev()    ! For calculating dGam/dT
      rotor(ir)%gamVec=matmul(rotor(ir)%AIC_inv,rotor(ir)%RHS)
    enddo

    ! Map gamVec to wing gam for each blade in rotor
    do ir=1,nr
      call rotor(ir)%map_gam()
    enddo

    ! Compute forces
    if (forcePlotSwitch .ne. 0) then
      if (mod(iter,forcePlotSwitch) .eq. 0) then 
        do ir=1,nr
          call rotor(ir)%calc_force(density,dt)
          call force2file(timestamp,rotor(ir),ir,-zAxis)  ! -ve sign due to negative inflow or gamma

          ! Compute alpha
          do ib=1,rotor(ir)%nb
            do is=1,rotor(ir)%ns
              do ic=1,rotor(ir)%nc
                ! Compute local velocity vector (excluding induced velocities from wing bound vortices)
                do jr=1,nr
                  rotor(ir)%blade(ib)%wiP(ic,is)%velLocal=rotor(ir)%blade(ib)%wiP(ic,is)%velLocal+  &
                    rotor(ir)%blade(ib)%wiP(ic,is)%velCP+  &
                    rotor(jr)%vind_bywing(rotor(ir)%blade(ib)%wiP(ic,is)%CP)-  &
                    rotor(jr)%vind_bywing_boundVortices(rotor(ir)%blade(ib)%wiP(ic,is)%CP)
                enddo

                rotor(ir)%blade(ib)%wiP(ic,is)%alpha=  &
                  dot_product(rotor(ir)%blade(ib)%wiP(ic,is)%velLocal,  &
                  rotor(ir)%blade(ib)%wiP(ic,is)%tauCapChord)/  &
                  norm2(rotor(ir)%blade(ib)%wiP(ic,is)%velLocal)
              enddo
            enddo
          enddo
        enddo
      endif
    endif

    ! Inflow plot
    do ir=1,nr
      if (rotor(ir)%inflowPlotSwitch .ne. 0) then
        if (mod(iter,rotor(ir)%inflowPlotSwitch) .eq. 0) then 
          call inflow2file(timestamp,rotor,ir,-zAxis)
        endif
      endif
    enddo

    ! Gamma plot
    do ir=1,nr
      if (rotor(ir)%gammaPlotSwitch .ne. 0) then
        if (mod(iter,rotor(ir)%gammaPlotSwitch) .eq. 0) then 
          call gamma2file(timestamp,rotor(1),ir)
        endif
      endif
    enddo

    ! Grid plot computation
    if (gridPlotSwitch .ne. 0) then
      if (mod(iter,gridPlotSwitch) .eq. 0) then 
        call filaments2file(timestamp,rotor)
      endif
    endif

    ! Induced vel on wake vortices
    do ir=1,nr
      do ib=1,rotor(ir)%nb
        rotor(ir)%blade(ib)%velNwake(:,rotor(ir)%rowNear:rotor(ir)%nNwake,:)=0._dp
        if (rotor(ir)%rowFar .ne. 0)  rotor(ir)%blade(ib)%velFwake(:,rotor(ir)%rowFar:rotor(ir)%nFwake)=0._dp
      enddo
    enddo

    do ir=1,nr
      do ib=1,rotor(ir)%nb
        do jr=1,nr
          rotor(ir)%blade(ib)%velNwake(:,rotor(ir)%rowNear:rotor(ir)%nNwake,:)=rotor(ir)%blade(ib)%velNwake(:,rotor(ir)%rowNear:rotor(ir)%nNwake,:)  &
            +  vind_onNwake_byRotor(rotor(jr),rotor(ir)%blade(ib)%waP(rotor(ir)%rowNear:rotor(ir)%nNwake,:))
          if (rotor(ir)%rowFar .ne. 0) then
            rotor(ir)%blade(ib)%velFwake(:,rotor(ir)%rowFar:rotor(ir)%nFwake)=rotor(ir)%blade(ib)%velFwake(:,rotor(ir)%rowFar:rotor(ir)%nFwake)  &
              +  vind_onFwake_byRotor(rotor(jr),rotor(ir)%blade(ib)%waF(rotor(ir)%rowFar:rotor(ir)%nFwake))
          endif
        enddo
        if (iter < initWakeVelNt) then
          do i=1,3
            rotor(ir)%blade(ib)%velNwake(i,rotor(ir)%rowNear:rotor(ir)%nNwake,:)=rotor(ir)%blade(ib)%velNwake(i,rotor(ir)%rowNear:rotor(ir)%nNwake,:)  &
              -  rotor(ir)%initWakeVel*rotor(ir)%shaftAxis(i)
          enddo
          if (rotor(ir)%rowFar .ne. 0) then
            do i=1,3
              rotor(ir)%blade(ib)%velFwake(i,rotor(ir)%rowFar:rotor(ir)%nFwake)=rotor(ir)%blade(ib)%velFwake(i,rotor(ir)%rowFar:rotor(ir)%nFwake)  &
                -  rotor(ir)%initWakeVel*rotor(ir)%shaftAxis(i)
            enddo
          endif
        endif
      enddo
    enddo

    ! Update wake vortex locations
    select case (fdSchemeSwitch)

    case (0)    ! Explicit forward diff (1st order)
      do ir=1,nr
        do ib=1,rotor(ir)%nb
          call rotor(ir)%blade(ib)%convectwake(rotor(ir)%rowNear,rotor(ir)%rowFar,dt,'C')
        enddo
      enddo


    case (1)    ! Predictor-Corrector (2nd order)
      ! Compute predicted wake
      do ir=1,nr
        do ib=1,rotor(ir)%nb
          rotor(ir)%blade(ib)%waPPredicted(rotor(ir)%rowNear:rotor(ir)%nNwake,:)=  &
            rotor(ir)%blade(ib)%waP(rotor(ir)%rowNear:rotor(ir)%nNwake,:)
          if (rotor(ir)%rowFar .ne. 0) then
            rotor(ir)%blade(ib)%waFPredicted(rotor(ir)%rowFar:rotor(ir)%nFwake)=  &
              rotor(ir)%blade(ib)%waF(rotor(ir)%rowFar:rotor(ir)%nFwake)
          endif
          call rotor(ir)%blade(ib)%convectwake(rotor(ir)%rowNear,rotor(ir)%rowFar,dt,'P')
        enddo
      enddo

      ! Compute velocity on predicted wake
      do ir=1,nr
        do ib=1,rotor(ir)%nb
          rotor(ir)%blade(ib)%velNwakePredicted(:,rotor(ir)%rowNear:rotor(ir)%nNwake,:)=0._dp
          if (rotor(ir)%rowFar .ne. 0)  rotor(ir)%blade(ib)%velFwakePredicted(:,rotor(ir)%rowFar:rotor(ir)%nFwake)=0._dp
          do jr=1,nr
            rotor(ir)%blade(ib)%velNwakePredicted(:,rotor(ir)%rowNear:rotor(ir)%nNwake,:)=rotor(ir)%blade(ib)%velNwakePredicted(:,rotor(ir)%rowNear:rotor(ir)%nNwake,:)  &
              +  vind_onNwake_byRotor(rotor(jr),rotor(ir)%blade(ib)%waPPredicted(rotor(ir)%rowNear:rotor(ir)%nNwake,:),'P')
            if (rotor(ir)%rowFar .ne. 0) then
              rotor(ir)%blade(ib)%velFwakePredicted(:,rotor(ir)%rowFar:rotor(ir)%nFwake)=rotor(ir)%blade(ib)%velFwakePredicted(:,rotor(ir)%rowFar:rotor(ir)%nFwake)  &
                +  vind_onFwake_byRotor(rotor(jr),rotor(ir)%blade(ib)%waFPredicted(rotor(ir)%rowFar:rotor(ir)%nFwake),'P')
            endif
          enddo
          if (iter < initWakeVelNt) then
            do i=1,3
              rotor(ir)%blade(ib)%velNwakePredicted(i,rotor(ir)%rowNear:rotor(ir)%nNwake,:)=rotor(ir)%blade(ib)%velNwakePredicted(i,rotor(ir)%rowNear:rotor(ir)%nNwake,:)  &
                -  rotor(ir)%initWakeVel*rotor(ir)%shaftAxis(i)
            enddo
            if (rotor(ir)%rowFar .ne. 0) then
              do i=1,3
                rotor(ir)%blade(ib)%velFwakePredicted(i,rotor(ir)%rowFar:rotor(ir)%nFwake)=rotor(ir)%blade(ib)%velFwakePredicted(i,rotor(ir)%rowFar:rotor(ir)%nFwake)  &
                  -  rotor(ir)%initWakeVel*rotor(ir)%shaftAxis(i)
              enddo
            endif
          endif
        enddo
      enddo

      ! Compute averaged velocity and convect wake
      do ir=1,nr
        do ib=1,rotor(ir)%nb
          rotor(ir)%blade(ib)%velNwake(:,rotor(ir)%rowNear:rotor(ir)%nNwake,:)=vel_order2_Nwake(  &
            rotor(ir)%blade(ib)%velNwake(:,rotor(ir)%rowNear:rotor(ir)%nNwake,:),rotor(ir)%blade(ib)%velNwakePredicted(:,rotor(ir)%rowNear:rotor(ir)%nNwake,:))
          !rotor(ir)%blade(ib)%velNwake(:,rotor(ir)%rowNear:rotor(ir)%nNwake,:)=0.5_dp*(  &
          !  rotor(ir)%blade(ib)%velNwake(:,rotor(ir)%rowNear:rotor(ir)%nNwake,:)+rotor(ir)%blade(ib)%velNwakePredicted(:,rotor(ir)%rowNear:rotor(ir)%nNwake,:))
          if (rotor(ir)%rowFar .ne. 0) then
            rotor(ir)%blade(ib)%velFwake(:,rotor(ir)%rowFar :rotor(ir)%nFwake  )=vel_order2_Fwake(  &
              rotor(ir)%blade(ib)%velFwake(:,rotor(ir)%rowFar :rotor(ir)%nFwake  ),rotor(ir)%blade(ib)%velFwakePredicted(:,rotor(ir)%rowFar :rotor(ir)%nFwake  ))
            !rotor(ir)%blade(ib)%velFwake(:,rotor(ir)%rowFar :rotor(ir)%nFwake  )=0.5_dp*(  &
            !  rotor(ir)%blade(ib)%velFwake(:,rotor(ir)%rowFar :rotor(ir)%nFwake  )+rotor(ir)%blade(ib)%velFwakePredicted(:,rotor(ir)%rowFar :rotor(ir)%nFwake  ))
          endif
          call rotor(ir)%blade(ib)%convectwake(rotor(ir)%rowNear,rotor(ir)%rowFar,dt,'C')
        enddo
      enddo


    case (2)    ! Adam-Bashforth (2nd order)
      if (iter == 1) then
        do ir=1,nr
          do ib=1,rotor(ir)%nb
            call rotor(ir)%blade(ib)%convectwake(rotor(ir)%rowNear,rotor(ir)%rowFar,dt,'C')
            rotor(ir)%blade(ib)%velNwake1=rotor(ir)%blade(ib)%velNwake
            if (rotor(ir)%rowFar .ne. 0)  rotor(ir)%blade(ib)%velFwake1=rotor(ir)%blade(ib)%velFwake
          enddo
        enddo
      else
        do ir=1,nr
          do ib=1,rotor(ir)%nb
            rotor(ir)%blade(ib)%velNwakeStep=0.5_dp*(3._dp*rotor(ir)%blade(ib)%velNwake-rotor(ir)%blade(ib)%velNwake1)
            if (rotor(ir)%rowFar .ne. 0) then
              rotor(ir)%blade(ib)%velFwakeStep=0.5_dp*(3._dp*rotor(ir)%blade(ib)%velFwake-rotor(ir)%blade(ib)%velFwake1)
            endif
            rotor(ir)%blade(ib)%velNwake1=rotor(ir)%blade(ib)%velNwake
            rotor(ir)%blade(ib)%velNwake=rotor(ir)%blade(ib)%velNwakeStep
            rotor(ir)%blade(ib)%velFwake1=rotor(ir)%blade(ib)%velFwake
            rotor(ir)%blade(ib)%velFwake=rotor(ir)%blade(ib)%velFwakeStep
            call rotor(ir)%blade(ib)%convectwake(rotor(ir)%rowNear,rotor(ir)%rowFar,dt,'C')
          enddo
        enddo
      endif


    case (3)    ! Predictor-Corrector Adam-Bashforth (4th order)
      if (iter == 1) then
        do ir=1,nr
          do ib=1,rotor(ir)%nb
            call rotor(ir)%blade(ib)%convectwake(rotor(ir)%rowNear,rotor(ir)%rowFar,dt,'C')
            rotor(ir)%blade(ib)%velNwake1=rotor(ir)%blade(ib)%velNwake
            if (rotor(ir)%rowFar .ne. 0)  rotor(ir)%blade(ib)%velFwake1=rotor(ir)%blade(ib)%velFwake
          enddo
        enddo
      elseif (iter == 2) then
        do ir=1,nr
          do ib=1,rotor(ir)%nb
            call rotor(ir)%blade(ib)%convectwake(rotor(ir)%rowNear,rotor(ir)%rowFar,dt,'C')
            rotor(ir)%blade(ib)%velNwake2=rotor(ir)%blade(ib)%velNwake
            if (rotor(ir)%rowFar .ne. 0)  rotor(ir)%blade(ib)%velFwake2=rotor(ir)%blade(ib)%velFwake
          enddo
        enddo
      elseif (iter == 3) then
        do ir=1,nr
          do ib=1,rotor(ir)%nb
            call rotor(ir)%blade(ib)%convectwake(rotor(ir)%rowNear,rotor(ir)%rowFar,dt,'C')
            rotor(ir)%blade(ib)%velNwake3=rotor(ir)%blade(ib)%velNwake
            if (rotor(ir)%rowFar .ne. 0)  rotor(ir)%blade(ib)%velFwake3=rotor(ir)%blade(ib)%velFwake
          enddo
        enddo
      else
        do ir=1,nr
          do ib=1,rotor(ir)%nb
            rotor(ir)%blade(ib)%waPPredicted(rotor(ir)%rowNear:rotor(ir)%nNwake,:)=rotor(ir)%blade(ib)%waP(rotor(ir)%rowNear:rotor(ir)%nNwake,:)
            rotor(ir)%blade(ib)%velNwakeStep=rotor(ir)%blade(ib)%velNwake  ! Store Nwake to Nwake_step for later use
            rotor(ir)%blade(ib)%velNwake = 55._dp/24._dp*rotor(ir)%blade(ib)%velNwake  &  ! Overwrite Nwake
              -59._dp/24._dp*rotor(ir)%blade(ib)%velNwake3 & 
              +37._dp/24._dp*rotor(ir)%blade(ib)%velNwake2 & 
              -09._dp/24._dp*rotor(ir)%blade(ib)%velNwake1  
            if (rotor(ir)%rowFar .ne. 0) then
              rotor(ir)%blade(ib)%waFPredicted(rotor(ir)%rowFar:rotor(ir)%nFwake)=rotor(ir)%blade(ib)%waF(rotor(ir)%rowFar:rotor(ir)%nFwake)
              rotor(ir)%blade(ib)%velFwakeStep=rotor(ir)%blade(ib)%velFwake  ! Store Fwake to Fwake_step for later use
              rotor(ir)%blade(ib)%velFwake = 55._dp/24._dp*rotor(ir)%blade(ib)%velFwake  &  ! Overwrite Fwake
                -59._dp/24._dp*rotor(ir)%blade(ib)%velFwake3 & 
                +37._dp/24._dp*rotor(ir)%blade(ib)%velFwake2 & 
                -09._dp/24._dp*rotor(ir)%blade(ib)%velFwake1  
            endif
            call rotor(ir)%blade(ib)%convectwake(rotor(ir)%rowNear,rotor(ir)%rowFar,dt,'P')
          enddo
        enddo

        do ir=1,nr
          do ib=1,rotor(ir)%nb
            rotor(ir)%blade(ib)%velNwakePredicted(:,rotor(ir)%rowNear:rotor(ir)%nNwake,:)=0._dp
            if (rotor(ir)%rowFar .ne. 0)  rotor(ir)%blade(ib)%velFwakePredicted(:,rotor(ir)%rowFar:rotor(ir)%nFwake)=0._dp
            do jr=1,nr
              rotor(ir)%blade(ib)%velNwakePredicted(:,rotor(ir)%rowNear:rotor(ir)%nNwake,:)=rotor(ir)%blade(ib)%velNwakePredicted(:,rotor(ir)%rowNear:rotor(ir)%nNwake,:)  &
                +  vind_onNwake_byRotor(rotor(jr),rotor(ir)%blade(ib)%waPPredicted(rotor(ir)%rowNear:rotor(ir)%nNwake,:),'P')
              if (rotor(ir)%rowFar .ne. 0) then
                rotor(ir)%blade(ib)%velFwakePredicted(:,rotor(ir)%rowFar:rotor(ir)%nFwake)=rotor(ir)%blade(ib)%velFwakePredicted(:,rotor(ir)%rowFar:rotor(ir)%nFwake)  &
                  +  vind_onFwake_byRotor(rotor(jr),rotor(ir)%blade(ib)%waFPredicted(rotor(ir)%rowFar:rotor(ir)%nFwake),'P')
              endif
            enddo
            if (iter < initWakeVelNt) then
              do i=1,3
                rotor(ir)%blade(ib)%velNwakePredicted(i,rotor(ir)%rowNear:rotor(ir)%nNwake,:)=rotor(ir)%blade(ib)%velNwakePredicted(i,rotor(ir)%rowNear:rotor(ir)%nNwake,:)  &
                  -  rotor(ir)%initWakeVel*rotor(ir)%shaftAxis(i)
              enddo
              if (rotor(ir)%rowFar .ne. 0) then
                do i=1,3
                  rotor(ir)%blade(ib)%velFwakePredicted(i,rotor(ir)%rowFar:rotor(ir)%nFwake)=rotor(ir)%blade(ib)%velFwakePredicted(i,rotor(ir)%rowFar:rotor(ir)%nFwake)  &
                    -  rotor(ir)%initWakeVel*rotor(ir)%shaftAxis(i)
                enddo
              endif
            endif
          enddo
        enddo

        do ir=1,nr
          do ib=1,rotor(ir)%nb
            rotor(ir)%blade(ib)%velNwake = 09._dp/24._dp*rotor(ir)%blade(ib)%velNwakePredicted  & 
              +19._dp/24._dp*rotor(ir)%blade(ib)%velNwakeStep  & 
              -05._dp/24._dp*rotor(ir)%blade(ib)%velNwake3 &  
              +01._dp/24._dp*rotor(ir)%blade(ib)%velNwake2  
            if (rotor(ir)%rowFar .ne. 0) then
              rotor(ir)%blade(ib)%velFwake = 09._dp/24._dp*rotor(ir)%blade(ib)%velFwakePredicted  & 
                +19._dp/24._dp*rotor(ir)%blade(ib)%velFwakeStep  & 
                -05._dp/24._dp*rotor(ir)%blade(ib)%velFwake3 &  
                +01._dp/24._dp*rotor(ir)%blade(ib)%velFwake2  
            endif
            call rotor(ir)%blade(ib)%convectwake(rotor(ir)%rowNear,rotor(ir)%rowFar,dt,'C') 

            rotor(ir)%blade(ib)%velNwake1=rotor(ir)%blade(ib)%velNwake2
            rotor(ir)%blade(ib)%velNwake2=rotor(ir)%blade(ib)%velNwake3
            rotor(ir)%blade(ib)%velNwake3=rotor(ir)%blade(ib)%velNwakeStep

            rotor(ir)%blade(ib)%velFwake1=rotor(ir)%blade(ib)%velFwake2
            rotor(ir)%blade(ib)%velFwake2=rotor(ir)%blade(ib)%velFwake3
            rotor(ir)%blade(ib)%velFwake3=rotor(ir)%blade(ib)%velFwakeStep
          enddo
        enddo
      endif

    end select

    ! Strain wake
    if (wakeStrainSwitch .eq. 1) then
      do ir=1,nr
        if (rotor(ir)%rowFar .ne. 0)  call rotor(ir)%strain_wake()
      enddo
    endif

    do ir=1,nr
      if ((rotor(ir)%rowNear .eq. 1) .and. (rotor(ir)%rowFar/=1)) then  ! Last step of near wake or later steps
        call rotor(ir)%rollup()    ! Rollup wake for next far wake panel
        call rotor(ir)%shiftwake()    ! Shift wake 
        call rotor(ir)%assignshed('TE')  ! Store shed vortex as TE for next near wake panel
      else
        call rotor(ir)%assignshed('TE')  
      endif
    enddo

  enddo

  if (wakePlotSwitch .eq. 1) then
    do ir=1,nr
      call rotor2file(timestamp,rotor(ir))
    enddo
  endif

  do ir=1,nr
    call rotor(ir)%deinit(fdSchemeSwitch)
  enddo

end program main
