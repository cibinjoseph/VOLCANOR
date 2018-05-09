program main
  use library
  use postproc

  ! Variables
  include "init_file.f90"

  ! Read config.in file
  open(unit=11,file='config.in')
  call skiplines(11,2)
  read(11,*) nt
  call skiplines(11,3)
  read(11,*) dt
  call skiplines(11,3)
  read(11,*) nr
  close(11)

  ! Allocate rotor objects
  allocate(rotor(nr))

  ! Read rotor??.in files
  do irotor=1,nr
    write(rotor_char,'(I0.2)') irotor
    rotorfile='rotor'//rotor_char//'.in'
    call rotor(irotor)%getdata(rotorfile,nt)
  enddo

  ! *** CHECK THIS ***
  ! Conversions
  vbody=-1._dp*vwind
  pqr=-1._dp*om_body

  ! Rotor and wake initialization
  do irotor=1,nr
    call rotor(irotor)%init(span_spacing_switch,nt,dt)
  enddo

  gamvec_prev=0._dp

  ! Rotate wing pc, vr, cp and ncap by initial pitch angle 
  do irotor=1,nr
    call rotor(irotor)%pitch(wing,theta_pitch,pivotLE)
  enddo

  ! Influence Coefficient Matrix 
  do irotor=1,nr
    do iblade=1,rotor(irotor)%nb
      do ispan=1,rotor(irotor)%ns      ! Collocation point loop
        do ichord=1,rotor(irotor)%nc
          row=ichord+rotor(irotor)%nc*(ispan-1)+rotor(irotor)%ns*rotor(irotor)%nc*(iblade-1)

          do jblade=1,rotor(irotor)%nb
            do j=1,rotor(irotor)%ns       ! Vortex ring loop
              do i=1,rotor(irotor)%nc   
                col=i+nc*(j-1)+rotor(irotor)%ns*rotor(irotor)%nc*(jblade-1)
                vec_dummy=wing(i,j)%vr%vind(wing(ichord,ispan)%cp)
                Amat(row,col,irotor)=dot_product(vec_dummy,wing(ichord,ispan)%ncap)
              enddo
            enddo
          enddo

        enddo
      enddo
    enddo
  enddo

  Amat_invi(:,:,irotor)=inv(Amat(:,:,irotor))

  ! Initial Solution
  if (slowstart_switch .eq. 0) then
    t=0._dp

    thetadot=thetas*om_theta*cos(om_theta*t)
    hdot=om_h*h0*cos(om_h*t)
    vel_plunge=(/0._dp,0._dp,hdot/)  

    indx=1
    do ispan=1,ns
      do ichord=1,nc
        RHS(indx) = dot_product(vwind-vel_plunge,wing(ichord,ispan)%ncap)

        ! Pitch vel
        wing(ichord,ispan)%vel_pitch=thetadot*wing(ichord,ispan)%r_hinge
        RHS(indx)= RHS(indx)+wing(ichord,ispan)%vel_pitch

        ! pqr vel
        RHS(indx)= RHS(indx)+dot_product(cross3(pqr,wing(ichord,ispan)%cp),wing(ichord,ispan)%ncap)
        indx=indx+1
      enddo
    enddo
    RHS=-1._dp*RHS

    gamvec=matmul(Amat_inv,RHS)
  else
    gamvec=0._dp
  endif

  gamvec_prev=gamvec

  ! Map gamvec to wing gam
  wing%vr%gam=reshape(gamvec,(/nc,ns/))    ! ns,nc due to transpose

  call assignshed(wake(nt,:),wing(nc,:),'TE')  ! Store shed vortex as TE

  dpts=om_body*dt     ! dphi dtheta dpsi
  pts=0._dp        ! phi theta psi

  ! ------- MAIN LOOP START -------
  do iter=1,nt
    print*,iter,nt
    write(timestamp,'(I0.5)') iter
    row_now=nt-(iter-1)

    select case (slowstart_switch)
    case (0)
      ! No slow start
    case (1)    ! Linear slope
      om_body_slow=min(real(slowstart_nt),real(iter+1))*om_body/slowstart_nt
      dpts=om_body_slow*dt     ! dphi dtheta dpsi
    case (2)    ! tanh slope
      om_body_slow=tanh(5._dp*iter/slowstart_nt)*om_body
      dpts=om_body_slow*dt     ! dphi dtheta dpsi
    case (3)    ! tanh slope
      om_body_slow=(tanh(6._dp*real(iter)/real(slowstart_nt)-3._dp)+1._dp)*0.5_dp*om_body
      dpts=om_body_slow*dt     ! dphi dtheta dpsi
    case default
      error stop "Assign correct slowstart_switch"
    end select

    t=t+dt
    pts=pts+dpts
    dtheta_pitch=theta_pitch
    theta_pitch=theta0+thetas*sin(om_theta*t)
    dtheta_pitch=dtheta_pitch-theta_pitch

    if (tip_diss_switch .eq. 1) then
      ! Age vortex filaments
      call age_wake(wake(row_now:nt,:),dt)

      ! Wake tip dissipation
      call dissipate_tip(wake(row_now:nt,:))
    endif

    ! Wing velocities
    thetadot=thetas*om_theta*cos(om_theta*t)
    hdot=om_h*h0*cos(om_h*t)
    vel_plunge=(/0._dp,0._dp,hdot/)

    ! Wing motion 
    call mov_wing(wing,(vbody+vel_plunge)*dt)    ! Translate wing
    hub_coords=hub_coords+(vbody+vel_plunge)*dt
    call rot_wing(wing,dpts,hub_coords,1)    ! Wing Global rotation
    call rot_about_axis(wing,dtheta_pitch,pivotLE)    ! Wing Pitch rotation

    call assignshed(wake(row_now,:),wing(nc,:),'LE')    ! Store shed vortex as LE

    ! Write out wing n' wake
    if (wakeplot_switch .eq. 2) call mesh2file(wing,wake(row_now:nt,:),'Results/wNw'//timestamp//'.tec')
    call tip2file(wing,wake(row_now:nt,:),'Results/tip'//timestamp//'.tec')
    gam_sectional=calcgam(wing)
    call gam2file(yvec,gam_sectional,'Results/gam'//timestamp//'.curve')

    ! Induced vel at coll. point (excluding pitch and wing induced velocities)
    call vind_CP(wing,vwind-vel_plunge,pqr,wake(row_now:nt,:))
    RHS=0._dp
    indx=1
    do ispan=1,ns
      do ichord=1,nc
        ! Normal component
        RHS(indx)=dot_product(wing(ichord,ispan)%velCP,wing(ichord,ispan)%ncap) 

        ! Pitch vel
        wing(ichord,ispan)%vel_pitch=thetadot*wing(ichord,ispan)%r_hinge
        RHS(indx)=RHS(indx)+wing(ichord,ispan)%vel_pitch

        indx=indx+1
      enddo
    enddo
    RHS=-1._dp*RHS

    gamvec_prev=gamvec    ! For calculating dGam/dT
    gamvec=matmul(Amat_inv,RHS)

    ! Map gamvec to wing gam
    wing%vr%gam=reshape(gamvec,(/nc,ns/))    ! ns,nc due to transpose

    ! Forces computation
    call calc_wingalpha(wing)
    lift(iter)=calclift(wing,gamvec_prev,dt)
    drag(iter)=calcdrag(wing,gamvec_prev,dt)

    ! Induced vel on wake vortices
    vind_wake(:,row_now:nt,:)=vind_onwake(wing,wake(row_now:nt,:))
    if (iter > wake_ignore_nt .or. wake_ignore_nt .eq. 0) then 
      vind_wake(:,row_now:nt,:)=vind_wake(:,row_now:nt,:)+vind_onwake(wake(row_now:nt,:),wake(row_now:nt,:))
    endif
    if (iter < init_wake_vel_nt .or. init_wake_vel .ne. 0)  vind_wake(3,row_now:nt,:)=vind_wake(3,row_now:nt,:)+init_wake_vel

    ! Update wake vortex locations
    select case (FDscheme_switch)

    case (0)    ! Explicit forward diff (1st order)
      call convectwake(wake(row_now:nt,:),vind_wake(:,row_now:nt,:)*dt)


    case (1)    ! Predictor-Corrector (2nd order)
      Pwake(row_now:nt,:)=wake(row_now:nt,:)
      call convectwake(Pwake(row_now:nt,:),vind_wake(:,row_now:nt,:)*dt)

      Pvind_wake(:,row_now:nt,:)=vind_onwake(wing,Pwake(row_now:nt,:))
      if (iter > wake_ignore_nt .or. wake_ignore_nt .eq. 0) then 
        Pvind_wake(:,row_now:nt,:)=Pvind_wake(:,row_now:nt,:)+vind_onwake(Pwake(row_now:nt,:),Pwake(row_now:nt,:))
      endif 
      if (iter < init_wake_vel_nt .or. init_wake_vel_nt .ne. 0) Pvind_wake(3,row_now:nt,:)=Pvind_wake(3,row_now:nt,:)+init_wake_vel

      call convectwake(wake(row_now:nt,:),(vind_wake(:,row_now:nt,:)+Pvind_wake(:,row_now:nt,:))*dt*0.5_dp)


    case (2)    ! Adam-Bashforth (2nd order)
      if (iter == 1) then
        call convectwake(wake(row_now:nt,:),vind_wake(:,row_now:nt,:)*dt)
        vind_wake1=vind_wake
      else
        vind_wake_step=0.5_dp*(3._dp*vind_wake-vind_wake1)
        call convectwake(wake(row_now:nt,:),vind_wake_step(:,row_now:nt,:)*dt)
        vind_wake1=vind_wake
      endif


    case (3)    ! Predictor-Cprrector Adam-Bashforth (4th order)
      if (iter == 1) then
        call convectwake(wake(row_now:nt,:),vind_wake(:,row_now:nt,:)*dt)
        vind_wake1=vind_wake
      elseif (iter == 2) then
        call convectwake(wake(row_now:nt,:),vind_wake(:,row_now:nt,:)*dt)
        vind_wake2=vind_wake
      elseif (iter == 3) then
        call convectwake(wake(row_now:nt,:),vind_wake(:,row_now:nt,:)*dt)
        vind_wake3=vind_wake
      else
        Pwake(row_now:nt,:)=wake(row_now:nt,:)
        vind_wake_step =55._dp/24._dp*vind_wake  &
          -59._dp/24._dp*vind_wake3 & 
          +37._dp/24._dp*vind_wake2 & 
          -09._dp/24._dp*vind_wake1  
        call convectwake(Pwake(row_now:nt,:),vind_wake_step(:,row_now:nt,:)*dt)

        Pvind_wake(:,row_now:nt,:)=vind_onwake(wing,Pwake(row_now:nt,:))
        if (iter > wake_ignore_nt .or. wake_ignore_nt .eq. 0) then 
          Pvind_wake(:,row_now:nt,:)=Pvind_wake(:,row_now:nt,:)+vind_onwake(Pwake(row_now:nt,:),Pwake(row_now:nt,:))
        endif 

        if (iter < init_wake_vel_nt .or. init_wake_vel_nt .ne. 0)  &
          Pvind_wake(3,row_now:nt,:)=Pvind_wake(3,row_now:nt,:)+init_wake_vel

        vind_wake_step =09._dp/24._dp*Pvind_wake  & 
          +19._dp/24._dp* vind_wake  & 
          -05._dp/24._dp* vind_wake3 &  
          +01._dp/24._dp* vind_wake2  
        call convectwake(wake(row_now:nt,:),vind_wake_step(:,row_now:nt,:)*dt)

        vind_wake1=vind_wake2
        vind_wake2=vind_wake3
        vind_wake3=vind_wake
      endif

    end select

    ! Strain wake
    if (wakestrain_switch .eq. 1) call strain_wake(wake(row_now:nt,:))

    ! Store shed vortex as TE for next wake panel
    if (row_now>1) call assignshed(wake(row_now-1,:),wing(nc,:),'TE')  

  enddo

  ! Postprocesing
  call lift2file(lift,'Results/lift.curve',(/dt,om_body(3),span,vwind(1)/))
  call drag2file(drag,'Results/drag.curve',(/dt,om_body(3),span,vwind(1)/))

  if (wakeplot_switch .eq. 1) call mesh2file(wing,wake(row_now:nt,:),'Results/wNw'//timestamp//'.tec')

  ! Deallocate variables
  select case (FDscheme_switch)
  case (1)
    deallocate(Pwake)
    deallocate(vind_wake1)
    deallocate(Pvind_wake)
  case (2)
    deallocate(vind_wake1)
    deallocate(vind_wake_step)
  case (3)
    deallocate(Pwake)
    deallocate(vind_wake1)
    deallocate(vind_wake2)
    deallocate(vind_wake3)
    deallocate(Pvind_wake)
    deallocate(vind_wake_step)
  end select

end program main
