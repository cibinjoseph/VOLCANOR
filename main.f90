program main
  use rotor_classdef
  !use postproc

  ! Variable Initialization
  include "init_file.f90"

  ! Read config.in file
  open(unit=11,file='config.in')
  call skiplines(11,2)
  read(11,*) nt,dt,nr
  call skiplines(11,5)
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

  ! Read rotor??.in files
  do ir=1,nr
    write(rotor_char,'(I0.2)') ir
    rotorfile='rotor'//rotor_char//'.in'
    call rotor(ir)%getdata(rotorfile,nt)
  enddo

  ! Conversions
  !vbody=-1._dp*vwind
  !pqr=-1._dp*om_body

  ! Rotor and wake initialization
  do ir=1,nr
    call rotor(ir)%init(span_spacing_switch,nt,dt)
  enddo

  ! Rotate wing pc, vr, cp and ncap by initial pitch angle 
  do ir=1,nr
    call rotor(ir)%pitch(wing,theta_pitch,pivotLE)
  enddo

  ! Compute AIC and AIC_inv matrices for rotors
  do ir=1,nr
    call rotor(ir)%calcAIC()
  enddo

  ! Initial Solution
  if (slowstart_switch .eq. 0) then
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
            rotor(ir)%RHS(row)=rotor(ir)%RHS(row)+dot_product(cross3(rotor(ir)%om_wind-rotor(ir)%Omega*this%shaft_axis  &
              ,rotor(ir)%blade(ib)%wiP(ic,is)%cp),rotor(ir)%blade(ib)%wiP(ic,is)%ncap)
          enddo
        enddo
      enddo
      rotor(ir)%RHS=-1._dp*rotor(ir)%RHS
      rotor(ir)%gamvec=matmul(rotor(ir)%AIC_inv,rotor(ir)%RHS)
    enddo

  else
    do ir=1,nr
      rotor(ir)%gamvec=0._dp
    enddo
  endif

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
    row_now=nt-(iter-1)

    select case (slowstart_switch)
    case (0)    ! No slow start
      rotor(ir)%Omega_slow=rotor(ir)%Omega
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
        call rotor(ir)%age_wake(row_now,dt)

        ! Wake tip dissipation
        call rotor(ir)%dissipate_tip(row_now)
      enddo
    endif

    ! Wing motion 
    do ir=1,nr
      call rotor(ir)%move(rotor(ir)%v_body*dt)
      call rotor(ir)%rot_pts(rotor(ir)%om_body*dt)
      call rotor(ir)%rot_advance(rotor(ir)%Omega_slow)
    enddo

    do ir=1,nr
      call rotor(ir)%assignshed(row_now,'LE')  ! Store shed vortex as TE
    enddo

    !    ! Write out wing n' wake
    !    if (wakeplot_switch .eq. 2) call mesh2file(wing,wake(row_now:nt,:),'Results/wNw'//timestamp//'.tec')
    !    call tip2file(wing,wake(row_now:nt,:),'Results/tip'//timestamp//'.tec')
    !    gam_sectional=calcgam(wing)
    !    call gam2file(yvec,gam_sectional,'Results/gam'//timestamp//'.curve')
    !
    ! Induced vel at coll. point (excluding pitch and wing induced velocities)
    !call vind_CP(wing,vwind-vel_plunge,pqr,wake(row_now:nt,:))
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
            +cross3(rotor(ir)%om_wind-rotor(ir)%Omega_slow*this%shaft_axis,rotor(ir)%blade(ib)%wiP(ic,is)%cp)

            ! Wake vel
            do jb=1,nb
              rotor(ir)%blade(ib)%wiP(ic,is)%velCP=rotor(ir)%blade(ib)%wiP(ic,is)%velCP  &
              +vind_bywake(rotor(ir)%blade(jb)%waP(row_now:nt,:),rotor(ir)%blade(ib)%wiP(ic,is)%cp)
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

    !
    !    gamvec_prev=gamvec    ! For calculating dGam/dT
    !    gamvec=matmul(Amat_inv,RHS)
    !
    !    ! Map gamvec to wing gam
    !    wing%vr%gam=reshape(gamvec,(/nc,ns/))    ! ns,nc due to transpose
    !
    !    ! Forces computation
    !    call calc_wingalpha(wing)
    !    lift(iter)=calclift(wing,gamvec_prev,dt)
    !    drag(iter)=calcdrag(wing,gamvec_prev,dt)
    !
    !    ! Induced vel on wake vortices
    !    vind_wake(:,row_now:nt,:)=vind_onwake(wing,wake(row_now:nt,:))
    !    if (iter > wake_ignore_nt .or. wake_ignore_nt .eq. 0) then 
    !      vind_wake(:,row_now:nt,:)=vind_wake(:,row_now:nt,:)+vind_onwake(wake(row_now:nt,:),wake(row_now:nt,:))
    !    endif
    !    if (iter < init_wake_vel_nt .or. init_wake_vel .ne. 0)  vind_wake(3,row_now:nt,:)=vind_wake(3,row_now:nt,:)+init_wake_vel
    !
    !    ! Update wake vortex locations
    !    select case (FDscheme_switch)
    !
    !    case (0)    ! Explicit forward diff (1st order)
    !      call convectwake(wake(row_now:nt,:),vind_wake(:,row_now:nt,:)*dt)
    !
    !
    !    case (1)    ! Predictor-Corrector (2nd order)
    !      Pwake(row_now:nt,:)=wake(row_now:nt,:)
    !      call convectwake(Pwake(row_now:nt,:),vind_wake(:,row_now:nt,:)*dt)
    !
    !      Pvind_wake(:,row_now:nt,:)=vind_onwake(wing,Pwake(row_now:nt,:))
    !      if (iter > wake_ignore_nt .or. wake_ignore_nt .eq. 0) then 
    !        Pvind_wake(:,row_now:nt,:)=Pvind_wake(:,row_now:nt,:)+vind_onwake(Pwake(row_now:nt,:),Pwake(row_now:nt,:))
    !      endif 
    !      if (iter < init_wake_vel_nt .or. init_wake_vel_nt .ne. 0) Pvind_wake(3,row_now:nt,:)=Pvind_wake(3,row_now:nt,:)+init_wake_vel
    !
    !      call convectwake(wake(row_now:nt,:),(vind_wake(:,row_now:nt,:)+Pvind_wake(:,row_now:nt,:))*dt*0.5_dp)
    !
    !
    !    case (2)    ! Adam-Bashforth (2nd order)
    !      if (iter == 1) then
    !        call convectwake(wake(row_now:nt,:),vind_wake(:,row_now:nt,:)*dt)
    !        vind_wake1=vind_wake
    !      else
    !        vind_wake_step=0.5_dp*(3._dp*vind_wake-vind_wake1)
    !        call convectwake(wake(row_now:nt,:),vind_wake_step(:,row_now:nt,:)*dt)
    !        vind_wake1=vind_wake
    !      endif
    !
    !
    !    case (3)    ! Predictor-Cprrector Adam-Bashforth (4th order)
    !      if (iter == 1) then
    !        call convectwake(wake(row_now:nt,:),vind_wake(:,row_now:nt,:)*dt)
    !        vind_wake1=vind_wake
    !      elseif (iter == 2) then
    !        call convectwake(wake(row_now:nt,:),vind_wake(:,row_now:nt,:)*dt)
    !        vind_wake2=vind_wake
    !      elseif (iter == 3) then
    !        call convectwake(wake(row_now:nt,:),vind_wake(:,row_now:nt,:)*dt)
    !        vind_wake3=vind_wake
    !      else
    !        Pwake(row_now:nt,:)=wake(row_now:nt,:)
    !        vind_wake_step =55._dp/24._dp*vind_wake  &
    !          -59._dp/24._dp*vind_wake3 & 
    !          +37._dp/24._dp*vind_wake2 & 
    !          -09._dp/24._dp*vind_wake1  
    !        call convectwake(Pwake(row_now:nt,:),vind_wake_step(:,row_now:nt,:)*dt)
    !
    !        Pvind_wake(:,row_now:nt,:)=vind_onwake(wing,Pwake(row_now:nt,:))
    !        if (iter > wake_ignore_nt .or. wake_ignore_nt .eq. 0) then 
    !          Pvind_wake(:,row_now:nt,:)=Pvind_wake(:,row_now:nt,:)+vind_onwake(Pwake(row_now:nt,:),Pwake(row_now:nt,:))
    !        endif 
    !
    !        if (iter < init_wake_vel_nt .or. init_wake_vel_nt .ne. 0)  &
    !          Pvind_wake(3,row_now:nt,:)=Pvind_wake(3,row_now:nt,:)+init_wake_vel
    !
    !        vind_wake_step =09._dp/24._dp*Pvind_wake  & 
    !          +19._dp/24._dp* vind_wake  & 
    !          -05._dp/24._dp* vind_wake3 &  
    !          +01._dp/24._dp* vind_wake2  
    !        call convectwake(wake(row_now:nt,:),vind_wake_step(:,row_now:nt,:)*dt)
    !
    !        vind_wake1=vind_wake2
    !        vind_wake2=vind_wake3
    !        vind_wake3=vind_wake
    !      endif
    !
    !    end select
    !
    !    ! Strain wake
    !    if (wakestrain_switch .eq. 1) call strain_wake(wake(row_now:nt,:))
    !
    !    ! Store shed vortex as TE for next wake panel
    !    if (row_now>1) call assignshed(wake(row_now-1,:),wing(nc,:),'TE')  
    !
  enddo
  !
  !  ! Postprocesing
  !  call lift2file(lift,'Results/lift.curve',(/dt,om_body(3),span,vwind(1)/))
  !  call drag2file(drag,'Results/drag.curve',(/dt,om_body(3),span,vwind(1)/))
  !
  !  if (wakeplot_switch .eq. 1) call mesh2file(wing,wake(row_now:nt,:),'Results/wNw'//timestamp//'.tec')
  !
  !  ! Deallocate variables
  !  select case (FDscheme_switch)
  !  case (1)
  !    deallocate(Pwake)
  !    deallocate(vind_wake1)
  !    deallocate(Pvind_wake)
  !  case (2)
  !    deallocate(vind_wake1)
  !    deallocate(vind_wake_step)
  !  case (3)
  !    deallocate(Pwake)
  !    deallocate(vind_wake1)
  !    deallocate(vind_wake2)
  !    deallocate(vind_wake3)
  !    deallocate(Pvind_wake)
  !    deallocate(vind_wake_step)
  !  end select

end program main
