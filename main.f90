program main
  use library
  use postproc

  ! Variable Initialization
  include "init_file.f90"
  print*

  ! Read config.in file
  call print_status('Reading file '//'config.in')
  open(unit=11,file='config.in')
  call skiplines(11,2)
  read(11,*) nt,dt,nr
  call skiplines(11,4)
  read(11,*) density, turb_visc
  call skiplines(11,4)
  read(11,*) span_spacing_switch
  call skiplines(11,4)
  read(11,*) tip_diss_switch, wakestrain_switch
  call skiplines(11,4)
  read(11,*) slowstart_switch, slowstart_nt
  call skiplines(11,4)
  read(11,*) wakeplot_switch, gridplot_switch, forceplot_switch
  call skiplines(11,4)
  read(11,*) FDscheme_switch
  call skiplines(11,4)
  read(11,*) wake_ignore_nt
  call skiplines(11,4)
  read(11,*) init_wake_vel_nt
  close(11)
  call print_status()    ! SUCCESS

  ! Allocate rotor objects
  allocate(rotor(nr))

  ! Read rotor??.in files
  do ir=1,nr
    write(rotor_char,'(I0.2)') ir
    rotorfile='rotor'//rotor_char//'.in'
    call print_status('Reading file '//rotorfile)
    call rotor(ir)%getdata(rotorfile,nt)
    call print_status()    ! SUCCESS
  enddo

  ! Rotor and wake initialization
  do ir=1,nr
    call rotor(ir)%init(dt,span_spacing_switch,FDscheme_switch)
  enddo

  ! Rotate wing pc, vr, cp and ncap by initial pitch angle 
  do ir=1,nr
    do ib=1,rotor(ir)%nb
      rotor(ir)%blade(ib)%theta=rotor(ir)%gettheta(rotor(ir)%psi_start,ib)
      call rotor(ir)%blade(ib)%rot_pitch(sign(1._dp,rotor(ir)%Omega)*rotor(ir)%blade(ib)%theta)
    enddo
  enddo

  ! Compute AIC and AIC_inv matrices for rotors
  call print_status('Computing AIC matrices')
  do ir=1,nr
    call rotor(ir)%calcAIC()
  enddo
  call print_status()    ! SUCCESS

  ! Initial Solution
  call print_status('Computing initial solution')
  if (slowstart_switch .ne. 0) then
    do ir=1,nr
      rotor(ir)%Omega_slow=0._dp
    enddo
  else 
    do ir=1,nr
      rotor(ir)%Omega_slow=rotor(ir)%Omega
    enddo
  endif
  t=0._dp

  ! Compute RHS
  do ir=1,nr
    do ib=1,rotor(ir)%nb
      do is=1,rotor(ir)%ns
        do ic=1,rotor(ir)%nc
          row=ic+rotor(ir)%nc*(is-1)+rotor(ir)%ns*rotor(ir)%nc*(ib-1)

          ! Rotational vel
          rotor(ir)%RHS(row) = dot_product(rotor(ir)%v_wind,rotor(ir)%blade(ib)%wiP(ic,is)%ncap)

          ! Rotational vel
          rotor(ir)%RHS(row)=rotor(ir)%RHS(row)+dot_product(cross3(rotor(ir)%om_wind  &
            ,rotor(ir)%blade(ib)%wiP(ic,is)%CP-rotor(ir)%CG_coords),rotor(ir)%blade(ib)%wiP(ic,is)%ncap)

          ! Omega vel
          rotor(ir)%RHS(row)=rotor(ir)%RHS(row)+dot_product(cross3(-rotor(ir)%Omega_slow*rotor(ir)%shaft_axis  &
            ,rotor(ir)%blade(ib)%wiP(ic,is)%CP-rotor(ir)%hub_coords),rotor(ir)%blade(ib)%wiP(ic,is)%ncap)

          ! Pitch vel
          !rotor(ir)%blade%(ib)%wing(ic,is)%vel_pitch=rotor(ir)%thetadot_pitch(0._dp,ib)*rotor(ir)%blade(ib)%wiP(ic,is)%r_hinge
          !rotor(ir)%RHS(row)= RHS(row)+wing(ib,ic,is)%vel_pitch
        enddo
      enddo
    enddo
    rotor(ir)%RHS=-1._dp*rotor(ir)%RHS
    rotor(ir)%gamvec=matmul(rotor(ir)%AIC_inv,rotor(ir)%RHS)
  enddo
  call print_status()    ! SUCCESS

  do ir=1,nr
    rotor(ir)%gamvec_prev=rotor(ir)%gamvec
  enddo

  ! Map gamvec to wing gam for each blade in rotor
  do ir=1,nr
    call rotor(ir)%map_gam()
  enddo

  do ir=1,nr
    rotor(ir)%row_far=0
    rotor(ir)%row_near=rotor(ir)%nNwake+1  ! Since assignshed TE assigns to row_near-1 panel
    call rotor(ir)%assignshed('TE')  ! Store shed vortex as TE
  enddo

  ! ------- MAIN LOOP START -------
  do iter=1,nt
    t=t+dt
    print*,iter,nt
    write(timestamp,'(I0.5)') iter
    do ir=1,nr
      rotor(ir)%row_near=max(rotor(ir)%nNwake-(iter-1),1)
      rotor(ir)%row_far=nt-(iter-1)
      if (iter<=rotor(ir)%nNwake) rotor(ir)%row_far=0    ! 0 implies roll up has not occured
    enddo

    select case (slowstart_switch)
    case (0)    ! No slow start
      do ir=1,nr
        rotor(ir)%Omega_slow=rotor(ir)%Omega
      enddo
    case (1)    ! Linear slope
      do ir=1,nr
        rotor(ir)%Omega_slow=min(real(slowstart_nt),real(iter+1))*rotor(ir)%Omega/slowstart_nt
      enddo
    case (2)    ! tanh slope
      do ir=1,nr
        rotor(ir)%Omega_slow=tanh(5._dp*iter/slowstart_nt)*rotor(ir)%Omega
      enddo
    case (3)    ! tanh slope
      do ir=1,nr
        rotor(ir)%Omega_slow=(tanh(6._dp*real(iter)/real(slowstart_nt)-3._dp)+1._dp)*0.5_dp*rotor(ir)%Omega
      enddo
    case default
      error stop "Assign correct slowstart_switch"
    end select

    if (tip_diss_switch .eq. 1) then
      do ir=1,nr
        if (rotor(ir)%row_far .ne. 0) then
          ! Age vortex filaments
          call rotor(ir)%age_wake(dt)

          ! Wake tip dissipation
          call rotor(ir)%dissipate_tip(turb_visc)
        endif
      enddo
    endif

    ! Wing motion 
    do ir=1,nr
      call rotor(ir)%move(rotor(ir)%v_body*dt)
      call rotor(ir)%rot_pts(rotor(ir)%om_body*dt,rotor(ir)%CG_coords,1)
      call rotor(ir)%rot_advance(rotor(ir)%Omega_slow*dt)
    enddo

    ! Assign LE of near wake
    do ir=1,nr
      call rotor(ir)%assignshed('LE')  ! Store shed vortex as LE
    enddo

    ! Write out wing n' wake
    do ir=1,nr
      if (wakeplot_switch .ne. 0) then
        if ((mod(iter,wakeplot_switch) .eq. 0) .and.(rotor(ir)%row_far .ne. 0) ) call rotor2file(rotor(ir),timestamp)
      endif
    enddo
    if (gridplot_switch .ne. 0) then
      if ((mod(iter,gridplot_switch) .eq. 0) .and. (minval(rotor%row_far) .ne. 0)) call filaments2file(rotor,timestamp)
    endif

    !    call tip2file(wing,wake(row_near:nt,:),'Results/tip'//timestamp//'.plt')
    !    gam_sectional=calcgam(wing)
    !    call gam2file(yvec,gam_sectional,'Results/gam'//timestamp//'.curve')

    ! Compute RHS
    do ir=1,nr
      rotor(ir)%RHS=0._dp
      do ib=1,rotor(ir)%nb
        do is=1,rotor(ir)%ns
          do ic=1,rotor(ir)%nc
            row=ic+rotor(ir)%nc*(is-1)+rotor(ir)%ns*rotor(ir)%nc*(ib-1)

            ! Translational vel
            rotor(ir)%blade(ib)%wiP(ic,is)%velCP=rotor(ir)%v_wind

            ! Rotational vel
            rotor(ir)%blade(ib)%wiP(ic,is)%velCP=rotor(ir)%blade(ib)%wiP(ic,is)%velCP  &
              +cross3(rotor(ir)%om_wind,rotor(ir)%blade(ib)%wiP(ic,is)%CP-rotor(ir)%CG_coords)

            ! Omega vel
            rotor(ir)%blade(ib)%wiP(ic,is)%velCP=rotor(ir)%blade(ib)%wiP(ic,is)%velCP  &
              +cross3(-rotor(ir)%Omega_slow*rotor(ir)%shaft_axis,rotor(ir)%blade(ib)%wiP(ic,is)%CP-rotor(ir)%hub_coords)

            ! Wake vel
            do jr=1,nr
              rotor(ir)%blade(ib)%wiP(ic,is)%velCP=rotor(ir)%blade(ib)%wiP(ic,is)%velCP  &
                +rotor(jr)%vind_bywake(rotor(ir)%blade(ib)%wiP(ic,is)%CP)
            enddo

            rotor(ir)%RHS(row)=dot_product(rotor(ir)%blade(ib)%wiP(ic,is)%velCP,rotor(ir)%blade(ib)%wiP(ic,is)%ncap)

            ! Pitch vel
            !wing(ib,ic,is)%vel_pitch=thetadot*wing(ib,ic,is)%r_hinge
            !RHS(row)=RHS(row)+wing(ib,ic,is)%vel_pitch
          enddo
        enddo
      enddo
      rotor(ir)%RHS=-1._dp*rotor(ir)%RHS
    enddo

    do ir=1,nr
      rotor(ir)%gamvec_prev=rotor(ir)%gamvec    ! For calculating dGam/dT
      rotor(ir)%gamvec=matmul(rotor(ir)%AIC_inv,rotor(ir)%RHS)
    enddo

    ! Map gamvec to wing gam for each blade in rotor
    do ir=1,nr
      call rotor(ir)%map_gam()
    enddo

    ! Forces computation
    do ir=1,nr
      call rotor(ir)%calc_thrust(density)
    enddo
    !    drag(iter)=calcdrag(wing,gamvec_prev,dt)

    if (forceplot_switch .ne. 0) then
      do ir=1,nr
        if (mod(iter,forceplot_switch) .eq. 0) call thrust2file(rotor(ir),ir,timestamp)
      enddo
    endif

    ! Induced vel on wake vortices
    do ir=1,nr
      do ib=1,rotor(ir)%nb
        rotor(ir)%blade(ib)%vind_Nwake(:,rotor(ir)%row_near:rotor(ir)%nNwake,:)=0._dp
        if (rotor(ir)%row_far .ne. 0)  rotor(ir)%blade(ib)%vind_Fwake(:,rotor(ir)%row_far:rotor(ir)%nFwake)=0._dp
      enddo
    enddo

    do ir=1,nr
      do ib=1,rotor(ir)%nb
        do jr=1,nr
          rotor(ir)%blade(ib)%vind_Nwake(:,rotor(ir)%row_near:rotor(ir)%nNwake,:)=rotor(ir)%blade(ib)%vind_Nwake(:,rotor(ir)%row_near:rotor(ir)%nNwake,:)  &
            +  vind_onNwake_byrotor(rotor(jr),rotor(ir)%blade(ib)%waP(rotor(ir)%row_near:rotor(ir)%nNwake,:))
          if (rotor(ir)%row_far .ne. 0) then
            rotor(ir)%blade(ib)%vind_Fwake(:,rotor(ir)%row_far:rotor(ir)%nFwake)=rotor(ir)%blade(ib)%vind_Fwake(:,rotor(ir)%row_far:rotor(ir)%nFwake)  &
              +  vind_onFwake_byrotor(rotor(jr),rotor(ir)%blade(ib)%waF(rotor(ir)%row_far:rotor(ir)%nFwake))
          endif
        enddo
        if (iter < init_wake_vel_nt) then
          do i=1,3
            rotor(ir)%blade(ib)%vind_Nwake(i,rotor(ir)%row_near:rotor(ir)%nNwake,:)=rotor(ir)%blade(ib)%vind_Nwake(i,rotor(ir)%row_near:rotor(ir)%nNwake,:)  &
              -  rotor(ir)%init_wake_vel*rotor(ir)%shaft_axis(i)
          enddo
          if (rotor(ir)%row_far .ne. 0) then
            do i=1,3
              rotor(ir)%blade(ib)%vind_Fwake(i,rotor(ir)%row_far:rotor(ir)%nFwake)=rotor(ir)%blade(ib)%vind_Fwake(i,rotor(ir)%row_far:rotor(ir)%nFwake)  &
                -  rotor(ir)%init_wake_vel*rotor(ir)%shaft_axis(i)
            enddo
          endif
        endif
      enddo
    enddo

    ! Update wake vortex locations
    select case (FDscheme_switch)

    case (0)    ! Explicit forward diff (1st order)
      do ir=1,nr
        do ib=1,rotor(ir)%nb
          call rotor(ir)%blade(ib)%convectwake(rotor(ir)%row_near,rotor(ir)%row_far,dt,'C')
        enddo
      enddo


    case (1)    ! Predictor-Corrector (2nd order)
      ! Compute predicted wake
      do ir=1,nr
        do ib=1,rotor(ir)%nb
          rotor(ir)%blade(ib)%waP_predicted(rotor(ir)%row_near:rotor(ir)%nNwake,:)=  &
            rotor(ir)%blade(ib)%waP(rotor(ir)%row_near:rotor(ir)%nNwake,:)
          if (rotor(ir)%row_far .ne. 0) then
            rotor(ir)%blade(ib)%waF_predicted(rotor(ir)%row_far:rotor(ir)%nFwake)=  &
              rotor(ir)%blade(ib)%waF(rotor(ir)%row_far:rotor(ir)%nFwake)
          endif
          call rotor(ir)%blade(ib)%convectwake(rotor(ir)%row_near,rotor(ir)%row_far,dt,'P')
        enddo
      enddo

      ! Compute velocity on predicted wake
      do ir=1,nr
        do ib=1,rotor(ir)%nb
          rotor(ir)%blade(ib)%vind_Nwake_predicted(:,rotor(ir)%row_near:rotor(ir)%nNwake,:)=0._dp
          if (rotor(ir)%row_far .ne. 0)  rotor(ir)%blade(ib)%vind_Fwake_predicted(:,rotor(ir)%row_far:rotor(ir)%nFwake)=0._dp
          do jr=1,nr
            rotor(ir)%blade(ib)%vind_Nwake_predicted(:,rotor(ir)%row_near:rotor(ir)%nNwake,:)=rotor(ir)%blade(ib)%vind_Nwake_predicted(:,rotor(ir)%row_near:rotor(ir)%nNwake,:)  &
              +  vind_onNwake_byrotor(rotor(jr),rotor(ir)%blade(ib)%waP_predicted(rotor(ir)%row_near:rotor(ir)%nNwake,:),'P')
            if (rotor(ir)%row_far .ne. 0) then
              rotor(ir)%blade(ib)%vind_Fwake_predicted(:,rotor(ir)%row_far:rotor(ir)%nFwake)=rotor(ir)%blade(ib)%vind_Fwake_predicted(:,rotor(ir)%row_far:rotor(ir)%nFwake)  &
                +  vind_onFwake_byrotor(rotor(jr),rotor(ir)%blade(ib)%waF_predicted(rotor(ir)%row_far:rotor(ir)%nFwake),'P')
            endif
          enddo
          if (iter < init_wake_vel_nt) then
            do i=1,3
              rotor(ir)%blade(ib)%vind_Nwake_predicted(i,rotor(ir)%row_near:rotor(ir)%nNwake,:)=rotor(ir)%blade(ib)%vind_Nwake_predicted(i,rotor(ir)%row_near:rotor(ir)%nNwake,:)  &
                -  rotor(ir)%init_wake_vel*rotor(ir)%shaft_axis(i)
            enddo
            if (rotor(ir)%row_far .ne. 0) then
              do i=1,3
                rotor(ir)%blade(ib)%vind_Fwake_predicted(i,rotor(ir)%row_far:rotor(ir)%nFwake)=rotor(ir)%blade(ib)%vind_Fwake_predicted(i,rotor(ir)%row_far:rotor(ir)%nFwake)  &
                  -  rotor(ir)%init_wake_vel*rotor(ir)%shaft_axis(i)
              enddo
            endif
          endif
        enddo
      enddo

      ! Compute averaged velocity and convect wake
      do ir=1,nr
        do ib=1,rotor(ir)%nb
          rotor(ir)%blade(ib)%vind_Nwake(:,rotor(ir)%row_near:rotor(ir)%nNwake,:)=vel_order2_Nwake(  &
            rotor(ir)%blade(ib)%vind_Nwake(:,rotor(ir)%row_near:rotor(ir)%nNwake,:),rotor(ir)%blade(ib)%vind_Nwake_predicted(:,rotor(ir)%row_near:rotor(ir)%nNwake,:))
          !rotor(ir)%blade(ib)%vind_Nwake(:,rotor(ir)%row_near:rotor(ir)%nNwake,:)=0.5_dp*(  &
          !  rotor(ir)%blade(ib)%vind_Nwake(:,rotor(ir)%row_near:rotor(ir)%nNwake,:)+rotor(ir)%blade(ib)%vind_Nwake_predicted(:,rotor(ir)%row_near:rotor(ir)%nNwake,:))
          if (rotor(ir)%row_far .ne. 0) then
            rotor(ir)%blade(ib)%vind_Fwake(:,rotor(ir)%row_far :rotor(ir)%nFwake  )=vel_order2_Fwake(  &
              rotor(ir)%blade(ib)%vind_Fwake(:,rotor(ir)%row_far :rotor(ir)%nFwake  ),rotor(ir)%blade(ib)%vind_Fwake_predicted(:,rotor(ir)%row_far :rotor(ir)%nFwake  ))
            !rotor(ir)%blade(ib)%vind_Fwake(:,rotor(ir)%row_far :rotor(ir)%nFwake  )=0.5_dp*(  &
            !  rotor(ir)%blade(ib)%vind_Fwake(:,rotor(ir)%row_far :rotor(ir)%nFwake  )+rotor(ir)%blade(ib)%vind_Fwake_predicted(:,rotor(ir)%row_far :rotor(ir)%nFwake  ))
          endif
          call rotor(ir)%blade(ib)%convectwake(rotor(ir)%row_near,rotor(ir)%row_far,dt,'C')
        enddo
      enddo


    case (2)    ! Adam-Bashforth (2nd order)
      if (iter == 1) then
        do ir=1,nr
          do ib=1,rotor(ir)%nb
            call rotor(ir)%blade(ib)%convectwake(rotor(ir)%row_near,rotor(ir)%row_far,dt,'C')
            rotor(ir)%blade(ib)%vind_Nwake1=rotor(ir)%blade(ib)%vind_Nwake
            if (rotor(ir)%row_far .ne. 0)  rotor(ir)%blade(ib)%vind_Fwake1=rotor(ir)%blade(ib)%vind_Fwake
          enddo
        enddo
      else
        do ir=1,nr
          do ib=1,rotor(ir)%nb
            rotor(ir)%blade(ib)%vind_Nwake_step=0.5_dp*(3._dp*rotor(ir)%blade(ib)%vind_Nwake-rotor(ir)%blade(ib)%vind_Nwake1)
            if (rotor(ir)%row_far .ne. 0) then
              rotor(ir)%blade(ib)%vind_Fwake_step=0.5_dp*(3._dp*rotor(ir)%blade(ib)%vind_Fwake-rotor(ir)%blade(ib)%vind_Fwake1)
            endif
            rotor(ir)%blade(ib)%vind_Nwake1=rotor(ir)%blade(ib)%vind_Nwake
            rotor(ir)%blade(ib)%vind_Nwake=rotor(ir)%blade(ib)%vind_Nwake_step
            rotor(ir)%blade(ib)%vind_Fwake1=rotor(ir)%blade(ib)%vind_Fwake
            rotor(ir)%blade(ib)%vind_Fwake=rotor(ir)%blade(ib)%vind_Fwake_step
            call rotor(ir)%blade(ib)%convectwake(rotor(ir)%row_near,rotor(ir)%row_far,dt,'C')
          enddo
        enddo
      endif


    case (3)    ! Predictor-Corrector Adam-Bashforth (4th order)
      if (iter == 1) then
        do ir=1,nr
          do ib=1,rotor(ir)%nb
            call rotor(ir)%blade(ib)%convectwake(rotor(ir)%row_near,rotor(ir)%row_far,dt,'C')
            rotor(ir)%blade(ib)%vind_Nwake1=rotor(ir)%blade(ib)%vind_Nwake
            if (rotor(ir)%row_far .ne. 0)  rotor(ir)%blade(ib)%vind_Fwake1=rotor(ir)%blade(ib)%vind_Fwake
          enddo
        enddo
      elseif (iter == 2) then
        do ir=1,nr
          do ib=1,rotor(ir)%nb
            call rotor(ir)%blade(ib)%convectwake(rotor(ir)%row_near,rotor(ir)%row_far,dt,'C')
            rotor(ir)%blade(ib)%vind_Nwake2=rotor(ir)%blade(ib)%vind_Nwake
            if (rotor(ir)%row_far .ne. 0)  rotor(ir)%blade(ib)%vind_Fwake2=rotor(ir)%blade(ib)%vind_Fwake
          enddo
        enddo
      elseif (iter == 3) then
        do ir=1,nr
          do ib=1,rotor(ir)%nb
            call rotor(ir)%blade(ib)%convectwake(rotor(ir)%row_near,rotor(ir)%row_far,dt,'C')
            rotor(ir)%blade(ib)%vind_Nwake3=rotor(ir)%blade(ib)%vind_Nwake
            if (rotor(ir)%row_far .ne. 0)  rotor(ir)%blade(ib)%vind_Fwake3=rotor(ir)%blade(ib)%vind_Fwake
          enddo
        enddo
      else
        do ir=1,nr
          do ib=1,rotor(ir)%nb
            rotor(ir)%blade(ib)%waP_predicted(rotor(ir)%row_near:rotor(ir)%nNwake,:)=rotor(ir)%blade(ib)%waP(rotor(ir)%row_near:rotor(ir)%nNwake,:)
            rotor(ir)%blade(ib)%vind_Nwake_step=rotor(ir)%blade(ib)%vind_Nwake  ! Store Nwake to Nwake_step for later use
            rotor(ir)%blade(ib)%vind_Nwake = 55._dp/24._dp*rotor(ir)%blade(ib)%vind_Nwake  &  ! Overwrite Nwake
              -59._dp/24._dp*rotor(ir)%blade(ib)%vind_Nwake3 & 
              +37._dp/24._dp*rotor(ir)%blade(ib)%vind_Nwake2 & 
              -09._dp/24._dp*rotor(ir)%blade(ib)%vind_Nwake1  
            if (rotor(ir)%row_far .ne. 0) then
              rotor(ir)%blade(ib)%waF_predicted(rotor(ir)%row_far:rotor(ir)%nFwake)=rotor(ir)%blade(ib)%waF(rotor(ir)%row_far:rotor(ir)%nFwake)
              rotor(ir)%blade(ib)%vind_Fwake_step=rotor(ir)%blade(ib)%vind_Fwake  ! Store Fwake to Fwake_step for later use
              rotor(ir)%blade(ib)%vind_Fwake = 55._dp/24._dp*rotor(ir)%blade(ib)%vind_Fwake  &  ! Overwrite Fwake
                -59._dp/24._dp*rotor(ir)%blade(ib)%vind_Fwake3 & 
                +37._dp/24._dp*rotor(ir)%blade(ib)%vind_Fwake2 & 
                -09._dp/24._dp*rotor(ir)%blade(ib)%vind_Fwake1  
            endif
            call rotor(ir)%blade(ib)%convectwake(rotor(ir)%row_near,rotor(ir)%row_far,dt,'P')
          enddo
        enddo

        do ir=1,nr
          do ib=1,rotor(ir)%nb
            rotor(ir)%blade(ib)%vind_Nwake_predicted(:,rotor(ir)%row_near:rotor(ir)%nNwake,:)=0._dp
            if (rotor(ir)%row_far .ne. 0)  rotor(ir)%blade(ib)%vind_Fwake_predicted(:,rotor(ir)%row_far:rotor(ir)%nFwake)=0._dp
            do jr=1,nr
              rotor(ir)%blade(ib)%vind_Nwake_predicted(:,rotor(ir)%row_near:rotor(ir)%nNwake,:)=rotor(ir)%blade(ib)%vind_Nwake_predicted(:,rotor(ir)%row_near:rotor(ir)%nNwake,:)  &
                +  vind_onNwake_byrotor(rotor(jr),rotor(ir)%blade(ib)%waP_predicted(rotor(ir)%row_near:rotor(ir)%nNwake,:),'P')
              if (rotor(ir)%row_far .ne. 0) then
                rotor(ir)%blade(ib)%vind_Fwake_predicted(:,rotor(ir)%row_far:rotor(ir)%nFwake)=rotor(ir)%blade(ib)%vind_Fwake_predicted(:,rotor(ir)%row_far:rotor(ir)%nFwake)  &
                  +  vind_onFwake_byrotor(rotor(jr),rotor(ir)%blade(ib)%waF_predicted(rotor(ir)%row_far:rotor(ir)%nFwake),'P')
              endif
            enddo
            if (iter < init_wake_vel_nt) then
              do i=1,3
                rotor(ir)%blade(ib)%vind_Nwake_predicted(i,rotor(ir)%row_near:rotor(ir)%nNwake,:)=rotor(ir)%blade(ib)%vind_Nwake_predicted(i,rotor(ir)%row_near:rotor(ir)%nNwake,:)  &
                  -  rotor(ir)%init_wake_vel*rotor(ir)%shaft_axis(i)
              enddo
              if (rotor(ir)%row_far .ne. 0) then
                do i=1,3
                  rotor(ir)%blade(ib)%vind_Fwake_predicted(i,rotor(ir)%row_far:rotor(ir)%nFwake)=rotor(ir)%blade(ib)%vind_Fwake_predicted(i,rotor(ir)%row_far:rotor(ir)%nFwake)  &
                    -  rotor(ir)%init_wake_vel*rotor(ir)%shaft_axis(i)
                enddo
              endif
            endif
          enddo
        enddo

        do ir=1,nr
          do ib=1,rotor(ir)%nb
            rotor(ir)%blade(ib)%vind_Nwake = 09._dp/24._dp*rotor(ir)%blade(ib)%vind_Nwake_predicted  & 
              +19._dp/24._dp*rotor(ir)%blade(ib)%vind_Nwake_step  & 
              -05._dp/24._dp*rotor(ir)%blade(ib)%vind_Nwake3 &  
              +01._dp/24._dp*rotor(ir)%blade(ib)%vind_Nwake2  
            if (rotor(ir)%row_far .ne. 0) then
              rotor(ir)%blade(ib)%vind_Fwake = 09._dp/24._dp*rotor(ir)%blade(ib)%vind_Fwake_predicted  & 
                +19._dp/24._dp*rotor(ir)%blade(ib)%vind_Fwake_step  & 
                -05._dp/24._dp*rotor(ir)%blade(ib)%vind_Fwake3 &  
                +01._dp/24._dp*rotor(ir)%blade(ib)%vind_Fwake2  
            endif
            call rotor(ir)%blade(ib)%convectwake(rotor(ir)%row_near,rotor(ir)%row_far,dt,'C') 

            rotor(ir)%blade(ib)%vind_Nwake1=rotor(ir)%blade(ib)%vind_Nwake2
            rotor(ir)%blade(ib)%vind_Nwake2=rotor(ir)%blade(ib)%vind_Nwake3
            rotor(ir)%blade(ib)%vind_Nwake3=rotor(ir)%blade(ib)%vind_Nwake_step

            rotor(ir)%blade(ib)%vind_Fwake1=rotor(ir)%blade(ib)%vind_Fwake2
            rotor(ir)%blade(ib)%vind_Fwake2=rotor(ir)%blade(ib)%vind_Fwake3
            rotor(ir)%blade(ib)%vind_Fwake3=rotor(ir)%blade(ib)%vind_Fwake_step
          enddo
        enddo
      endif

    end select

    ! Strain wake
    if (wakestrain_switch .eq. 1) then
      do ir=1,nr
        if (rotor(ir)%row_far .ne. 0)  call rotor(ir)%strain_wake()
      enddo
    endif

    do ir=1,nr
      if ((rotor(ir)%row_near .eq. 1) .and. (rotor(ir)%row_far/=1)) then  ! Last step of near wake or later steps
        call rotor(ir)%rollup()    ! Rollup wake for next far wake panel
        call rotor(ir)%shiftwake()    ! Shift wake 
        call rotor(ir)%assignshed('TE')  ! Store shed vortex as TE for next near wake panel
      else
        call rotor(ir)%assignshed('TE')  
      endif
    enddo

  enddo

  ! Postprocesing
  !do ir=1,nr
  !  call thrust2file(rotor(ir),ir,timestamp)
  !enddo

  !  call drag2file(drag,'Results/drag.curve',(/dt,om_body(3),span,vwind(1)/))

  if (wakeplot_switch .eq. 1) then
    do ir=1,nr
      call rotor2file(rotor(ir),timestamp)
    enddo
  endif

  do ir=1,nr
    call rotor(ir)%deinit(FDscheme_switch)
  enddo

end program main
