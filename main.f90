program main
  use library
  use postproc

  ! Variable Initialization
  include "init_file.f90"

  ! Read config.in file
  open(unit=11,file='config.in')
  call skiplines(11,2)
  read(11,*) nt,dt,nr
  call skiplines(11,4)
  read(11,*) density
  call skiplines(11,4)
  read(11,*) span_spacing_switch
  call skiplines(11,4)
  read(11,*) tip_diss_switch, wakestrain_switch
  call skiplines(11,4)
  read(11,*) slowstart_switch, slowstart_nt
  call skiplines(11,4)
  read(11,*) wakeplot_switch 
  call skiplines(11,4)
  read(11,*) FDscheme_switch
  call skiplines(11,4)
  read(11,*) wake_ignore_nt
  call skiplines(11,4)
  read(11,*) init_wake_vel_nt
  close(11)

  ! Allocate rotor objects
  allocate(rotor(nr))
  allocate(row_now(nr))

  ! Read rotor??.in files
  do ir=1,nr
    write(rotor_char,'(I0.2)') ir
    rotorfile='rotor'//rotor_char//'.in'
    call rotor(ir)%getdata(rotorfile,nt)
  enddo

  ! Rotor and wake initialization
  do ir=1,nr
    call rotor(ir)%init(nt,dt,span_spacing_switch,FDscheme_switch)
  enddo

  ! Rotate wing pc, vr, cp and ncap by initial pitch angle 
  do ir=1,nr
    do ib=1,rotor(ir)%nb
      rotor(ir)%blade(ib)%theta=rotor(ir)%gettheta(rotor(ir)%psi_start,ib)
      call rotor(ir)%blade(ib)%rot_pitch(rotor(ir)%blade(ib)%theta)
    enddo
  enddo

  ! Compute AIC and AIC_inv matrices for rotors
  do ir=1,nr
    call rotor(ir)%calcAIC()
  enddo

  ! Initial Solution
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

  do ir=1,nr
    do ib=1,rotor(ir)%nb
      do is=1,rotor(ir)%ns
        do ic=1,rotor(ir)%nc
          row=ic+rotor(ir)%nc*(is-1)+rotor(ir)%ns*rotor(ir)%nc*(ib-1)
          rotor(ir)%RHS(row) = dot_product(rotor(ir)%v_wind,rotor(ir)%blade(ib)%wiP(ic,is)%ncap)

          ! Pitch vel
          !rotor(ir)%blade%(ib)%wing(ic,is)%vel_pitch=rotor(ir)%thetadot_pitch(0._dp,ib)*rotor(ir)%blade(ib)%wiP(ic,is)%r_hinge
          !rotor(ir)%RHS(row)= RHS(row)+wing(ib,ic,is)%vel_pitch

          ! pqr vel
          rotor(ir)%RHS(row)=rotor(ir)%RHS(row)+dot_product(cross3(rotor(ir)%om_wind-rotor(ir)%Omega_slow*rotor(ir)%shaft_axis  &
            ,rotor(ir)%blade(ib)%wiP(ic,is)%cp),rotor(ir)%blade(ib)%wiP(ic,is)%ncap)
        enddo
      enddo
    enddo
    rotor(ir)%RHS=-1._dp*rotor(ir)%RHS
    rotor(ir)%gamvec=matmul(rotor(ir)%AIC_inv,rotor(ir)%RHS)
  enddo

  do ir=1,nr
    rotor(ir)%gamvec_prev=rotor(ir)%gamvec
  enddo

  ! Map gamvec to wing gam for each blade in rotor
  do ir=1,nr
    call rotor(ir)%map_gam()
  enddo

  do ir=1,nr
    call rotor(ir)%assignshed(nt,'TE')  ! Store shed vortex as TE
  enddo

  ! ------- MAIN LOOP START -------
  do iter=1,nt
    t=t+dt
    print*,iter,nt
    write(timestamp,'(I0.5)') iter
    do ir=1,nr
      row_now(ir)=rotor(ir)%nNwake-(iter-1)
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
        ! Age vortex filaments
        call rotor(ir)%age_wake(row_now(ir),dt)

        ! Wake tip dissipation
        call rotor(ir)%dissipate_tip(row_now(ir))
      enddo
    endif

    ! Wing motion 
    do ir=1,nr
      call rotor(ir)%move(rotor(ir)%v_body*dt)
      call rotor(ir)%rot_pts(rotor(ir)%om_body*dt,rotor(ir)%CG_coords,1)
      call rotor(ir)%rot_advance(rotor(ir)%Omega_slow*dt)
    enddo

    do ir=1,nr
      call rotor(ir)%assignshed(row_now(ir),'LE')  ! Store shed vortex as TE
    enddo

    !    ! Write out wing n' wake
    do ir=1,nr
      if (wakeplot_switch .eq. 2) call rotor2file(rotor(ir),row_now(ir),'Results/wNw'//timestamp//'.tec')
    enddo
    !    call tip2file(wing,wake(row_now:nt,:),'Results/tip'//timestamp//'.tec')
    !    gam_sectional=calcgam(wing)
    !    call gam2file(yvec,gam_sectional,'Results/gam'//timestamp//'.curve')

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
              +cross3(rotor(ir)%om_wind-rotor(ir)%Omega_slow*rotor(ir)%shaft_axis,rotor(ir)%blade(ib)%wiP(ic,is)%cp)

            ! Wake vel
            do jr=1,nr
              rotor(ir)%blade(ib)%wiP(ic,is)%velCP=rotor(ir)%blade(ib)%wiP(ic,is)%velCP  &
                +rotor(jr)%vind_bywake(row_now(ir),rotor(ir)%blade(ib)%wiP(ic,is)%cp)
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

    !    ! Forces computation
    !    call calc_wingalpha(wing)
    !    lift(iter)=calclift(wing,gamvec_prev,dt)
    !    drag(iter)=calcdrag(wing,gamvec_prev,dt)

    ! Induced vel on wake vortices
    do ir=1,nr
      do ib=1,rotor(ir)%nb
        rotor(ir)%blade(ib)%vind_wake(:,row_now(ir):nt,:)=0._dp
      enddo
    enddo

    do ir=1,nr
      do ib=1,rotor(ir)%nb
        do jr=1,nr
          rotor(ir)%blade(ib)%vind_wake(:,row_now(ir):nt,:)=rotor(ir)%blade(ib)%vind_wake(:,row_now(ir):nt,:)  &
            +  vind_onwake_byrotor(rotor(jr),rotor(ir)%blade(ib)%waP(row_now(ir):nt,:))
        enddo
        if (iter < init_wake_vel_nt) then
          do i=1,3
            rotor(ir)%blade(ib)%vind_wake(i,row_now(ir):nt,:)=rotor(ir)%blade(ib)%vind_wake(i,row_now(ir):nt,:)  &
              -  rotor(ir)%init_wake_vel*rotor(ir)%shaft_axis(i)
          enddo
        endif
      enddo
    enddo

    ! Update wake vortex locations
    select case (FDscheme_switch)

    case (0)    ! Explicit forward diff (1st order)
      do ir=1,nr
        do ib=1,rotor(ir)%nb
          call rotor(ir)%blade(ib)%convectwake(rotor(ir)%blade(ib)%vind_wake(:,row_now(ir):nt,:)*dt)
        enddo
      enddo


    case (1)    ! Predictor-Corrector (2nd order)
      do ir=1,nr
        do ib=1,rotor(ir)%nb
          rotor(ir)%blade(ib)%Pwake(row_now(ir):nt,:)=rotor(ir)%blade(ib)%waP(row_now(ir):nt,:)
          call rotor(ir)%blade(ib)%convectwake(rotor(ir)%blade(ib)%vind_wake(:,row_now(ir):nt,:)*dt,'P')
        enddo
      enddo

      do ir=1,nr
        do ib=1,rotor(ir)%nb
          rotor(ir)%blade(ib)%Pvind_wake(:,row_now(ir):nt,:)=0._dp
          do jr=1,nr
            rotor(ir)%blade(ib)%Pvind_wake(:,row_now(ir):nt,:)=rotor(ir)%blade(ib)%Pvind_wake(:,row_now(ir):nt,:)  &
              +  vind_onwake_byrotor(rotor(jr),rotor(ir)%blade(ib)%Pwake(row_now(ir):nt,:),'P')
          enddo
          if (iter < init_wake_vel_nt) then
            do i=1,3
              rotor(ir)%blade(ib)%Pvind_wake(i,row_now(ir):nt,:)=rotor(ir)%blade(ib)%Pvind_wake(i,row_now(ir):nt,:)  &
                -  rotor(ir)%init_wake_vel*rotor(ir)%shaft_axis(i)
            enddo
          endif
        enddo
      enddo

      do ir=1,nr
        do ib=1,rotor(ir)%nb
          call rotor(ir)%blade(ib)%convectwake(vel_order2(rotor(ir)%blade(ib)%vind_wake(:,row_now(ir):nt,:)  &
            ,  rotor(ir)%blade(ib)%Pvind_wake(:,row_now(ir):nt,:))*dt)
        enddo
      enddo

      !do ir=1,nr
      !  do ib=1,rotor(ir)%nb
      !    call rotor(ir)%blade(ib)%convectwake((rotor(ir)%blade(ib)%vind_wake(:,row_now(ir):nt,:)  &
      !      +  rotor(ir)%blade(ib)%Pvind_wake(:,row_now(ir):nt,:))*dt*0.5_dp)
      !  enddo
      !enddo

    case (2)    ! Adam-Bashforth (2nd order)
      if (iter == 1) then
        do ir=1,nr
          do ib=1,rotor(ir)%nb
            call rotor(ir)%blade(ib)%convectwake(rotor(ir)%blade(ib)%vind_wake(:,row_now(ir):nt,:)*dt)
            rotor(ir)%blade(ib)%vind_wake1=rotor(ir)%blade(ib)%vind_wake
          enddo
        enddo
      else
        do ir=1,nr
          do ib=1,rotor(ir)%nb
            rotor(ir)%blade(ib)%vind_wake_step=0.5_dp*(3._dp*rotor(ir)%blade(ib)%vind_wake-rotor(ir)%blade(ib)%vind_wake1)
            call rotor(ir)%blade(ib)%convectwake(rotor(ir)%blade(ib)%vind_wake_step(:,row_now(ir):nt,:)*dt)
            rotor(ir)%blade(ib)%vind_wake1=rotor(ir)%blade(ib)%vind_wake
          enddo
        enddo
      endif


    case (3)    ! Predictor-Corrector Adam-Bashforth (4th order)
      if (iter == 1) then
        do ir=1,nr
          do ib=1,rotor(ir)%nb
            call rotor(ir)%blade(ib)%convectwake(rotor(ir)%blade(ib)%vind_wake(:,row_now(ir):nt,:)*dt)
            rotor(ir)%blade(ib)%vind_wake1=rotor(ir)%blade(ib)%vind_wake
          enddo
        enddo
      elseif (iter == 2) then
        do ir=1,nr
          do ib=1,rotor(ir)%nb
            call rotor(ir)%blade(ib)%convectwake(rotor(ir)%blade(ib)%vind_wake(:,row_now(ir):nt,:)*dt)
            rotor(ir)%blade(ib)%vind_wake2=rotor(ir)%blade(ib)%vind_wake
          enddo
        enddo
      elseif (iter == 3) then
        do ir=1,nr
          do ib=1,rotor(ir)%nb
            call rotor(ir)%blade(ib)%convectwake(rotor(ir)%blade(ib)%vind_wake(:,row_now(ir):nt,:)*dt)
            rotor(ir)%blade(ib)%vind_wake3=rotor(ir)%blade(ib)%vind_wake
          enddo
        enddo
      else
        do ir=1,nr
          do ib=1,rotor(ir)%nb
            rotor(ir)%blade(ib)%Pwake(row_now(ir):nt,:)=rotor(ir)%blade(ib)%waP(row_now(ir):nt,:)
            rotor(ir)%blade(ib)%vind_wake_step =55._dp/24._dp*rotor(ir)%blade(ib)%vind_wake  &
              -59._dp/24._dp*rotor(ir)%blade(ib)%vind_wake3 & 
              +37._dp/24._dp*rotor(ir)%blade(ib)%vind_wake2 & 
              -09._dp/24._dp*rotor(ir)%blade(ib)%vind_wake1  
            call rotor(ir)%blade(ib)%convectwake(rotor(ir)%blade(ib)%vind_wake_step(:,row_now(ir):nt,:)*dt,'P')
          enddo
        enddo

        do ir=1,nr
          do ib=1,rotor(ir)%nb
            rotor(ir)%blade(ib)%Pvind_wake(:,row_now(ir):nt,:)=0._dp
            do jr=1,nr
              rotor(ir)%blade(ib)%Pvind_wake(:,row_now(ir):nt,:)=rotor(ir)%blade(ib)%Pvind_wake(:,row_now(ir):nt,:)  &
                +  vind_onwake_byrotor(rotor(jr),rotor(ir)%blade(ib)%Pwake(row_now(ir):nt,:),'P')
            enddo
            if (iter < init_wake_vel_nt) then
              do i=1,3
                rotor(ir)%blade(ib)%Pvind_wake(i,row_now(ir):nt,:)=rotor(ir)%blade(ib)%Pvind_wake(i,row_now(ir):nt,:)  &
                  -  rotor(ir)%init_wake_vel*rotor(ir)%shaft_axis(i)
              enddo
            endif
          enddo
        enddo

        do ir=1,nr
          do ib=1,rotor(ir)%nb
            rotor(ir)%blade(ib)%vind_wake_step =09._dp/24._dp*rotor(ir)%blade(ib)%Pvind_wake  & 
              +19._dp/24._dp*rotor(ir)%blade(ib)%vind_wake  & 
              -05._dp/24._dp*rotor(ir)%blade(ib)%vind_wake3 &  
              +01._dp/24._dp*rotor(ir)%blade(ib)%vind_wake2  
            call rotor(ir)%blade(ib)%convectwake(rotor(ir)%blade(ib)%vind_wake_step(:,row_now(ir):nt,:)*dt)

            rotor(ir)%blade(ib)%vind_wake1=rotor(ir)%blade(ib)%vind_wake2
            rotor(ir)%blade(ib)%vind_wake2=rotor(ir)%blade(ib)%vind_wake3
            rotor(ir)%blade(ib)%vind_wake3=rotor(ir)%blade(ib)%vind_wake
          enddo
        enddo
      endif

    end select

    do ir=1,nr
      ! Strain wake
      if (wakestrain_switch .eq. 1) call rotor(ir)%strain_wake(row_now(ir))
    enddo

    do ir=1,nr
      ! Store shed vortex as TE for next wake panel
      if (row_now(ir)>1) call rotor(ir)%assignshed(row_now(ir)-1,'TE')  
    enddo

  enddo

  !  ! Postprocesing
  !  call lift2file(lift,'Results/lift.curve',(/dt,om_body(3),span,vwind(1)/))
  !  call drag2file(drag,'Results/drag.curve',(/dt,om_body(3),span,vwind(1)/))
  !
  !  if (wakeplot_switch .eq. 1) call mesh2file(wing,wake(row_now(ir):nt,:),'Results/wNw'//timestamp//'.tec')

  do ir=1,nr
    call rotor(ir)%deinit(FDscheme_switch)
  enddo

end program main
