program main
  use library
  use postproc

  ! Variables
  include "init_file.f90"

  ! Read wing data from file
  open(unit=12,file='inputfile')
  read(12,*)
  read(12,*)
  read(12,*) chord,span,dt,vwind(1),vwind(2),vwind(3)
  read(12,*)
  read(12,*)
  read(12,*)
  read(12,*) theta0,thetac,thetas,om_theta
  read(12,*)
  read(12,*)
  read(12,*)
  read(12,*) h0,om_h
  read(12,*)
  read(12,*)
  read(12,*)
  read(12,*) om_body(1),om_body(2),om_body(3)
  close(12)

  ! Conversions
  call degtorad(theta0)
  call degtorad(thetac)
  call degtorad(thetas)
  om_theta=2._dp*pi*om_theta
  om_h    =2._dp*pi*om_h
  vbody=-1._dp*vwind
  pqr=-1._dp*om_body

  ! Geometry Definition
  xvec=linspace(0._dp,chord,nc+1)
  select case (span_spacing_switch)
  case (1)
    yvec=linspace(0.05_dp*span,span+0.05_dp*span,ns+1)
  case (2)
    yvec=cosspace(0.05_dp*span,span+0.05_dp*span,ns+1)
  end select

  ! Initialize wake geometry and core radius
  call init_wake(wake,0.075*span)
  gamvec_prev=0._dp

  ! Initialize wing geometry, vr, cp, ncap coordinates and core radius
  call init_wing(wing,xvec,yvec,0.001_dp*span)

  ! Rotate wing pc, vr, cp and ncap by initial pitch angle 
  theta_pitch=theta0
  call rot_wing(wing,(/0._dp,theta_pitch,0._dp/),1)

  !  TE vortex position
  v_shed=0.20*vwind
  do ispan=1,ns
    call wing(nc,ispan)%vr%shiftdP(2,v_shed*dt)
    call wing(nc,ispan)%vr%shiftdP(3,v_shed*dt)
  enddo

  ! Influence Coefficient Matrix 
  row=1
  col=1
  do ispan=1,ns      ! Collocation point loop
    do ichord=1,nc
      row=ichord+nc*(ispan-1)
      do j=1,ns       ! Vortex ring loop
        do i=1,nc   
          col=i+nc*(j-1)
          vec_dummy=wing(i,j)%vr%vind(wing(ichord,ispan)%cp)
          Amat(row,col)=dot_product(vec_dummy,wing(ichord,ispan)%ncap)
        enddo
      enddo
    enddo
  enddo

  Amat_inv=inv(Amat)

  ! Initial Solution
  t=0._dp

  thetadot=thetas*om_theta*cos(om_theta*t)
  hdot=om_h*h0*cos(om_h*t)
  vel_plunge=(/0._dp,0._dp,-hdot/)  

  indx=1
  do ispan=1,ns
    do ichord=1,nc
      RHS(indx) = dot_product(vwind+vel_plunge,wing(ichord,ispan)%ncap)

      ! Pitch vel
      vel_pitch=thetadot*wing(ichord,ispan)%r_hinge
      RHS(indx)= RHS(indx)+vel_pitch

      ! pqr vel
      RHS(indx)= RHS(indx)+dot_product(cross3(pqr,wing(ichord,ispan)%cp),wing(ichord,ispan)%ncap)
      indx=indx+1
    enddo
  enddo
  RHS=-1._dp*RHS

  gamvec=matmul(Amat_inv,RHS)
  gamvec_prev=gamvec

  ! Map gamvec to wing gam
  wing%vr%gam=reshape(gamvec,(/nc,ns/))    ! ns,nc due to transpose

  call assignshed(wake(nt,:),wing(nc,:),'TE')  ! Store shed vortex as TE

  dpts=om_body*dt     ! dphi dtheta dpsi
  pts=0._dp        ! phi theta psi

  open(unit=12,file='Results/lift.tec',position='append')

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
    vel_plunge=(/0._dp,0._dp,-hdot/)

    ! Wing motion 
    call mov_wing(wing,(vbody+vel_plunge)*dt) ! Translate wing
    call rot_wing(wing,dpts,1)                ! Wing Global rotation
    call pitch_wing(wing,dtheta_pitch,pts)    ! Wing Pitch rotation

    call assignshed(wake(row_now,:),wing(nc,:),'LE')  ! Store shed vortex as LE

    ! Write out wing n' wake
    call mesh2file(wing,wake(row_now:nt,:),'Results/wNw'//timestamp//'.tec')
    call tip2file(wing,wake(row_now:nt,:),'Results/tip'//timestamp//'.tec')

    ! Induced vel at coll. point (excluding pitch and wing induced velocities)
    call vind_CP(wing,vwind+vel_plunge,pqr,wake(row_now:nt,:))
    RHS=0._dp
    indx=1
    do ispan=1,ns
      do ichord=1,nc
        ! Normal component
        RHS(indx)=dot_product(wing(ichord,ispan)%velCP,wing(ichord,ispan)%ncap) 

        ! Pitch vel
        vel_pitch=thetadot*wing(ichord,ispan)%r_hinge
        RHS(indx)=RHS(indx)+vel_pitch

        indx=indx+1
      enddo
    enddo
    RHS=-1._dp*RHS

    gamvec_prev=gamvec    ! For calculating dGam/dT
    gamvec=matmul(Amat_inv,RHS)

    ! Map gamvec to wing gam
    wing%vr%gam=reshape(gamvec,(/nc,ns/))    ! ns,nc due to transpose

    ! Induced vel on wake vortices
    vind_wake(:,row_now:nt,:)=vind_onwake(wing,wake(row_now:nt,:))
    vind_wake(:,row_now:nt,:)=vind_wake(:,row_now:nt,:)+vind_onwake(wake(row_now:nt,:),wake(row_now:nt,:))

    ! Update wake vortex locations
    call convectwake(wake(row_now:nt,:),vind_wake(:,row_now:nt,:)*dt)
    call wake_continuity(wake(row_now:nt,:))   

    ! Strain wake
    ! call strain_wake(wake(row_now,:,:))

    if (row_now>1) then
      call assignshed(wake(row_now-1,:),wing(nc,:),'TE')  ! Store shed vortex as TE for next wake panel
    endif

    ! Lift computation
    call calclift(wing,gamvec_prev,dt)
    call lift2file(wing,'Results/lift.tec',(/t,chord,span,vwind(1)/))

  enddo


end program main
