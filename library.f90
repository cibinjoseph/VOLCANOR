module library

  use wingpanel_classdef
  use wakepanel_classdef
  implicit none

  ! Input parameters
  !integer, parameter :: nt = 500
  !integer, parameter :: ns = 13
  !integer, parameter :: nc = 1

  ! Global env parameters
  real(dp), parameter :: density = 1.2_dp

  ! Overloaded functions
  interface vind_panelgeo
    module procedure vind_panelgeo_wing, vind_panelgeo_wake
  end interface 
  interface vind_onwake
    module procedure vind_onwake_bywing, vind_onwake_bywake
  end interface 

contains

  !--------------------------------------------------------!
  !                Initialization Functions                !
  !--------------------------------------------------------!

  ! Assigns coordinates to all corners of pc and vr
  subroutine init_wing(wing_array,xvec,yvec,mid_core_radius,tip_core_radius)
    type(wingpanel_class), intent(out), dimension(:,:) :: wing_array
    real(dp), intent(in), dimension(:) :: xvec 
    real(dp), intent(in), dimension(:) :: yvec
    real(dp), intent(in) :: mid_core_radius, tip_core_radius
    real(dp) :: xshiftLE, xshiftTE
    integer :: i,j,rows,cols

    rows=size(wing_array,1)
    cols=size(wing_array,2)

    if ((size(xvec) .ne. rows+1) .and. (size(yvec) .ne. cols+1)) then
      error stop 'Size mismatch between xvec, yvec and panel_array'
    endif

    ! Initialize panel coordinates
    do j=1,cols
      do i=1,rows
        call wing_array(i,j)%assignP(1,(/xvec(i  ),yvec(j  ),0._dp/))
        call wing_array(i,j)%assignP(2,(/xvec(i+1),yvec(j  ),0._dp/))
        call wing_array(i,j)%assignP(3,(/xvec(i+1),yvec(j+1),0._dp/))
        call wing_array(i,j)%assignP(4,(/xvec(i  ),yvec(j+1),0._dp/))
      enddo
    enddo

    ! Initialize vr coords of all panels except last row (to accomodate mismatch of vr coords when using unequal spacing)
    do i=1,rows-1
      xshiftLE=(xvec(i+1)-xvec(i))*0.25_dp  ! Shift x coord by dx/4
      xshiftTE=(xvec(i+2)-xvec(i+1))*0.25_dp  ! Shift x coord by dx/4
      do j=1,cols
        call wing_array(i,j)%vr%assignP(1,(/xvec(i  )+xshiftLE,yvec(j  ),0._dp/))
        call wing_array(i,j)%vr%assignP(2,(/xvec(i+1)+xshiftTE,yvec(j  ),0._dp/))
        call wing_array(i,j)%vr%assignP(3,(/xvec(i+1)+xshiftTE,yvec(j+1),0._dp/))
        call wing_array(i,j)%vr%assignP(4,(/xvec(i  )+xshiftLE,yvec(j+1),0._dp/))
      enddo
    enddo

    ! Initializing vr coords of last row
    xshiftLE=(xvec(rows+1)-xvec(rows))*0.25_dp  ! Shift x coord by dx/4
    xshiftTE=0._dp
    do j=1,cols
      call wing_array(rows,j)%vr%assignP(1,(/xvec(rows  )+xshiftLE,yvec(j  ),0._dp/))
      call wing_array(rows,j)%vr%assignP(2,(/xvec(rows+1)+xshiftTE,yvec(j  ),0._dp/))
      call wing_array(rows,j)%vr%assignP(3,(/xvec(rows+1)+xshiftTE,yvec(j+1),0._dp/))
      call wing_array(rows,j)%vr%assignP(4,(/xvec(rows  )+xshiftLE,yvec(j+1),0._dp/))
    enddo

    ! Initialize CP coords, ncap, panel_area and r_hinge
    do j=1,cols
      do i=1,rows
        call wing_array(i,j)%calcCP()
        call wing_array(i,j)%calcN()
        wing_array(i,j)%r_hinge=length3d((wing_array(1,j)%pc(:,1)+wing_array(1,j)%pc(:,4))*0.5_dp,wing_array(i,j)%cp)
        call wing_array(i,j)%calc_area()
      enddo
    enddo

    ! Initialize gamma
    wing_array%vr%gam=0._dp

    ! Initialize tag
    wing_array%tag=2.

    ! Initialize mid vortex core radius
    do i=1,4
      wing_array%vr%vf(i)%r_vc0=mid_core_radius
      wing_array%vr%vf(i)%r_vc =mid_core_radius
      wing_array%vr%vf(i)%age=0._dp
    enddo

    ! Initialize tip vortex core radius
    do i=1,rows
      wing_array(i,1)%vr%vf(1)%r_vc0    = tip_core_radius
      wing_array(i,1)%vr%vf(1)%r_vc     = tip_core_radius
      wing_array(i,cols)%vr%vf(3)%r_vc0 = tip_core_radius
      wing_array(i,cols)%vr%vf(3)%r_vc  = tip_core_radius
    enddo

    ! Verify CP is outside vortex core for boundary panels
    if (isCPinsidecore(wing_array(1,1))) then
      print*,'Warning: CP inside vortex core at panel LU'
      print*,'Any key to continue. Ctrl-C to exit'
      read(*,*)
    endif
    if (isCPinsidecore(wing_array(rows,1))) then
      print*,'Warning: CP inside vortex core at panel LB'
      print*,'Any key to continue. Ctrl-C to exit'
      read(*,*)
    endif
    if (isCPinsidecore(wing_array(1,cols))) then
      print*,'Warning: CP inside vortex core at panel RU'
      print*,'Any key to continue. Ctrl-C to exit'
      read(*,*)
    endif
    if (isCPinsidecore(wing_array(rows,cols))) then
      print*,'Warning: CP inside vortex core at panel RB'
      print*,'Any key to continue. Ctrl-C to exit'
      read(*,*)
    endif

  end subroutine init_wing

  ! Assigns vortex code radii to all filaments
  subroutine init_wake(wake_array,mid_core_radius,tip_core_radius,starting_vortex_core)
    type(wakepanel_class), intent(out), dimension(:,:) :: wake_array
    real(dp), intent(in) :: mid_core_radius, tip_core_radius, starting_vortex_core
    integer :: i,j,cols,rows

    cols=size(wake_array,2)
    rows=size(wake_array,1)

    ! Assign core_radius to mid vortices
    do i=1,4
      wake_array%vr%vf(i)%r_vc0=mid_core_radius
      wake_array%vr%vf(i)%r_vc =mid_core_radius
      wake_array%vr%vf(i)%age=0._dp
    enddo

    wake_array%tag=-1
    wake_array%vr%gam=0._dp

    ! Assign core_radius to tip vortices
    do i=1,rows
      ! Root vortex 
      wake_array(i,1)%vr%vf(1)%r_vc0      = tip_core_radius 
      wake_array(i,1)%vr%vf(1)%r_vc       = tip_core_radius 
      wake_array(i,1)%vr%vf(3)%r_vc0      = tip_core_radius 
      wake_array(i,1)%vr%vf(3)%r_vc       = tip_core_radius 

      wake_array(i,2)%vr%vf(1)%r_vc0      = tip_core_radius 
      wake_array(i,2)%vr%vf(1)%r_vc       = tip_core_radius 
      wake_array(i,2)%vr%vf(3)%r_vc0      = tip_core_radius 
      wake_array(i,2)%vr%vf(3)%r_vc       = tip_core_radius 

      ! Tip vortex 
      wake_array(i,cols)%vr%vf(1)%r_vc0   = tip_core_radius 
      wake_array(i,cols)%vr%vf(1)%r_vc    = tip_core_radius 
      wake_array(i,cols)%vr%vf(3)%r_vc0   = tip_core_radius 
      wake_array(i,cols)%vr%vf(3)%r_vc    = tip_core_radius 

      wake_array(i,cols-1)%vr%vf(3)%r_vc0 = tip_core_radius 
      wake_array(i,cols-1)%vr%vf(3)%r_vc  = tip_core_radius 
      wake_array(i,cols-1)%vr%vf(3)%r_vc0 = tip_core_radius 
      wake_array(i,cols-1)%vr%vf(3)%r_vc  = tip_core_radius 
    enddo

    if (starting_vortex_core > eps) then
      ! Assign core_radius to starting vortices
      do i=1,cols
        do j=2,4,2
          wake_array(rows,i)%vr%vf(j)%r_vc0 = starting_vortex_core
          wake_array(rows,i)%vr%vf(j)%r_vc  = starting_vortex_core
          wake_array(rows-1,i)%vr%vf(j)%r_vc0 = starting_vortex_core
          wake_array(rows-1,i)%vr%vf(j)%r_vc  = starting_vortex_core
        enddo
      enddo
    endif

  end subroutine init_wake

  ! Checks whether CP lies inside viscous core region of vortex ring
  function isCPinsidecore(wing_panel)
    type(wingpanel_class), intent(in) :: wing_panel
    logical :: isCPinsidecore
    real(dp) :: deltaxby4, deltayby2

    deltaxby4=0.25_dp*abs(wing_panel%vr%vf(1)%fc(1,1)-wing_panel%vr%vf(2)%fc(1,1))
    deltayby2=0.5_dp *abs(wing_panel%vr%vf(1)%fc(2,1)-wing_panel%vr%vf(4)%fc(2,1))

    isCPinsidecore = .false.
    if (deltayby2 .lt. wing_panel%vr%vf(1)%r_vc) then
      isCPinsidecore = .true.    ! Left edge
    elseif (deltayby2 .lt. wing_panel%vr%vf(3)%r_vc) then 
      isCPinsidecore = .true.  ! Right edge
    elseif (deltaxby4 .lt. wing_panel%vr%vf(2)%r_vc) then
      isCPinsidecore = .true.  ! Upper edge
    elseif (3._dp*deltaxby4 .lt. wing_panel%vr%vf(4)%r_vc) then
      isCPinsidecore = .true.  ! Bottom edge
    endif
  end function isCPinsidecore
  !--------------------------------------------------------!
  !                 Wing Motion Functions                  !
  !--------------------------------------------------------!

  ! Transformation matrix bg
  function Tbg(cs_phi,cs_theta,cs_psi)
    real(dp), dimension(2), intent(in) :: cs_phi, cs_theta, cs_psi  ! cos and sin
    real(dp), dimension(3,3) :: Tbg
    Tbg(1,:)=(/cs_psi(1)*cs_theta(1),cs_theta(1)*cs_psi(2),-1._dp*cs_theta(2)/)
    Tbg(2,1)=cs_psi(1)*cs_phi(2)*cs_theta(2)-cs_phi(1)*cs_psi(2) 
    Tbg(2,2)=cs_phi(1)*cs_psi(1)+cs_phi(2)*cs_psi(2)*cs_theta(2)
    Tbg(2,3)=cs_theta(1)*cs_phi(2)
    Tbg(3,1)=cs_phi(1)*cs_psi(1)*cs_theta(2)+cs_phi(2)*cs_psi(2)
    Tbg(3,2)=cs_phi(1)*cs_psi(2)*cs_theta(2)-cs_psi(1)*cs_phi(2)
    Tbg(3,3)=cs_phi(1)*cs_theta(1)
  end function Tbg

  function Tgb(cs_phi,cs_theta,cs_psi)
    real(dp), dimension(2), intent(in) :: cs_phi, cs_theta, cs_psi  ! cos and sin
    real(dp), dimension(3,3) :: Tgb
    Tgb(1,1)=cs_psi(1)*cs_theta(1)
    Tgb(1,2)=cs_phi(2)*cs_theta(2)*cs_psi(1)-cs_psi(2)*cs_phi(1)
    Tgb(1,3)=cs_phi(2)*cs_psi(2)+cs_theta(2)*cs_phi(1)*cs_psi(1)
    Tgb(2,1)=cs_psi(2)*cs_theta(1)
    Tgb(2,2)=cs_phi(2)*cs_psi(2)*cs_theta(2)+cs_phi(1)*cs_psi(1)
    Tgb(2,3)=cs_psi(2)*cs_theta(2)*cs_phi(1)-cs_phi(2)*cs_psi(1)
    Tgb(3,1)=-cs_theta(2)
    Tgb(3,2)=cs_phi(2)*cs_theta(1)
    Tgb(3,3)=cs_phi(1)*cs_theta(1)
  end function Tgb

  !ALTERNATE blade_class: rot_pts
  !subroutine rot_wing(wing_array,pts,origin,order)  
  !  type(wingpanel_class), intent(inout), dimension(:,:) :: wing_array
  !  real(dp), dimension(3), intent(in) :: pts    ! pts => phi,theta,psi
  !  real(dp), dimension(3), intent(in) :: origin ! rotation about
  !  integer :: i, j
  !  integer :: order    ! [1]gb & +ve theta , [2]bg & -ve theta       
  !  real(dp), dimension(3,3) :: TMat

  !  select case (order)
  !  case (2)
  !    TMat=Tbg((/cos(pts(1)),sin(pts(1))/),&
  !      (/cos(pts(2)),sin(pts(2))/),&
  !      (/cos(pts(3)),sin(pts(3))/))
  !  case (1)
  !    TMat=Tgb((/cos(pts(1)),sin(pts(1))/),&
  !      (/cos(pts(2)),sin(pts(2))/),&
  !      (/cos(pts(3)),sin(pts(3))/))
  !  case default
  !    error stop 'Error: wrong option for order'
  !  end select

  !  do j=1,size(wing_array,2)
  !    do i=1,size(wing_array,1)
  !      call wing_array(i,j)%shiftdP(-origin)
  !      call wing_array(i,j)%rot(TMat)
  !      call wing_array(i,j)%shiftdP(origin)
  !    enddo
  !  enddo

  !end subroutine rot_wing

  !ALTERNATE blade_class: rot_axis
  !subroutine rot_about_axis(wing_array,theta,pivotLE)  !pitch about pivotLE from LE
  !  ! pivot point calculated using straight line joining LE and TE of root panels
  !  type(wingpanel_class), intent(inout), dimension(:,:) :: wing_array
  !  real(dp), intent(in) :: theta, pivotLE
  !  real(dp), dimension(3,3) :: Tmat
  !  integer :: i,j,rows
  !  real(dp), dimension(3) :: u  ! Axis vector
  !  real(dp), dimension(3) :: dshift  
  !  real(dp) :: ct,st,omct

  !  if (abs(theta)>eps) then
  !    ! Translate to origin
  !    rows=size(wing_array,1)
  !    dshift=wing_array(1,1)%pc(:,1)*(1._dp-pivotLE)+wing_array(rows,1)%pc(:,2)*pivotLE
  !    call mov_wing(wing_array,-dshift)

  !    ! Construct axes of rotation from LE of first panel
  !    u=wing_array(1,1)%pc(:,4)-wing_array(1,1)%pc(:,1)

  !    ! Calculate TMat
  !    ct=cos(theta)
  !    st=sin(theta)
  !    omct=1-ct

  !    Tmat(:,1)=(/      ct+u(1)*u(1)*omct,  u(3)*st+u(2)*u(1)*omct, -u(2)*st+u(3)*u(1)*omct/)
  !    Tmat(:,2)=(/-u(3)*st+u(1)*u(2)*omct,       ct+u(2)*u(2)*omct,  u(1)*st+u(3)*u(2)*omct/)
  !    Tmat(:,3)=(/ u(2)*st+u(1)*u(3)*omct, -u(1)*st+u(2)*u(3)*omct,       ct+u(3)*u(3)*omct/)

  !    ! Rotate about axis
  !    do j=1,size(wing_array,2)
  !      do i=1,size(wing_array,1)
  !        call wing_array(i,j)%rot(TMat)
  !      enddo
  !    enddo

  !    ! Untranslate from origin
  !    call mov_wing(wing_array,dshift)
  !  endif
  !end subroutine rot_about_axis

  ! OBSOLETE SUBROUTINE - use rot_about_axis
  !subroutine pitch_wing(wing_array,theta_pitch,pts)  !pitch about a fixed point fp
  !  ! assuming motion is 2dimensional and body has not undergone 3d rotations
  !  type(wingpanel_class), intent(inout), dimension(:,:) :: wing_array
  !  real(dp), intent(in) :: theta_pitch
  !  real(dp), dimension(3), intent(in) :: pts
  !  real(dp), dimension(3) :: dshift, origin

  !  origin=0._dp

  !  if (abs(theta_pitch)>eps) then
  !    ! Translate to origin
  !    dshift=(/wing_array(1,1)%pc(:,1)/)
  !    call mov_wing(wing_array,-dshift)

  !    ! Rotate global angles
  !    call rot_wing(wing_array,pts,origin,1)

  !    ! Rotate pitch angle
  !    call rot_wing(wing_array,(/0._dp,-theta_pitch,0._dp/),origin,1)

  !    ! Unrotate global angles
  !    call rot_wing(wing_array,-1._dp*pts,origin,1)

  !    ! Untranslate from origin
  !    call mov_wing(wing_array,dshift)
  !  endif
  !end subroutine pitch_wing

  ! ALTERNATE blade_class : mov_blade()
  ! ALTERNATE rotor_class : mov_rotor()
  !subroutine mov_wing(wing_array,dshift)
  !  type(wingpanel_class), intent(inout), dimension(:,:) :: wing_array
  !  real(dp), intent(in), dimension(3) :: dshift
  !  integer :: i,j

  !  do j=1,size(wing_array,2)
  !    do i=1,size(wing_array,1)
  !      call wing_array(i,j)%shiftdP(dshift)
  !    enddo
  !  enddo
  !end subroutine mov_wing

  !--------------------------------------------------------!
  !                    Wake Functions                      !
  !--------------------------------------------------------!

  subroutine assignshed(wake_row,wing_row,edge)      
    ! assigns coordinates to wake_row from wing_row
    type(wakepanel_class), intent(inout), dimension(:) :: wake_row
    type(wingpanel_class), intent(in), dimension(:) :: wing_row
    character(len=2), intent(in) :: edge
    integer :: i

    wake_row%vr%gam=wing_row%vr%gam

    select case (edge)
    case ('LE')    ! assign to LE 
      do i=1,size(wing_row)
        call wake_row(i)%vr%assignP(1,wing_row(i)%vr%vf(2)%fc(:,1))
        call wake_row(i)%vr%assignP(4,wing_row(i)%vr%vf(3)%fc(:,1))
        call wake_row(i)%vr%calclength(.TRUE.)    ! TRUE => record original length
      enddo
      wake_row%tag=1
    case ('TE')    ! assign to TE
      do i=1,size(wing_row)
        call wake_row(i)%vr%assignP(2,wing_row(i)%vr%vf(2)%fc(:,1))
        call wake_row(i)%vr%assignP(3,wing_row(i)%vr%vf(3)%fc(:,1))
      enddo
    case default
      error stop 'Error: Wrong option for edge'
    end select

  end subroutine assignshed

  ! Convect wake using dP_array=vind_array*dt
  subroutine convectwake(wake_array,dP_array)
    type(wakepanel_class), intent(inout), dimension(:,:) :: wake_array
    real(dp), intent(in), dimension(:,:,:) :: dP_array
    integer :: i,j,rows,cols

    rows=size(wake_array,1)
    cols=size(wake_array,2)

    !$omp parallel do collapse(2)
    do j=1,cols
      do i=1,rows
        call wake_array(i,j)%vr%shiftdP(2,dP_array(:,i,j))
      enddo
    enddo
    !$omp end parallel do

    !$omp parallel do
    do i=1,rows
      call wake_array(i,cols)%vr%shiftdP(3,dP_array(:,i,cols+1))
    enddo
    !$omp end parallel do
    call wake_continuity(wake_array)
  end subroutine convectwake

  ! Maintain continuity between vortex ring elements after convection
  ! of vortex ring corners
  subroutine wake_continuity(wake_array)
    type(wakepanel_class), intent(inout), dimension(:,:) :: wake_array
    integer :: i,j,rows,cols

    rows=size(wake_array,1)
    cols=size(wake_array,2)

    !$omp parallel do collapse(2)
    do j=1,cols-1
      do i=2,rows
        call wake_array(i,j)%vr%assignP(1,wake_array(i-1,j)%vr%vf(2)%fc(:,1))
        call wake_array(i,j)%vr%assignP(3,wake_array(i,j+1)%vr%vf(2)%fc(:,1))
        call wake_array(i,j)%vr%assignP(4,wake_array(i-1,j+1)%vr%vf(2)%fc(:,1))
      enddo
    enddo
    !$omp end parallel do

    !$omp parallel do
    do j=1,cols-1
      call wake_array(1,j)%vr%assignP(3,wake_array(1,j+1)%vr%vf(2)%fc(:,1))
    enddo
    !$omp end parallel do

    !$omp parallel do
    do i=2,rows
      call wake_array(i,cols)%vr%assignP(1,wake_array(i-1,cols)%vr%vf(2)%fc(:,1))
      call wake_array(i,cols)%vr%assignP(4,wake_array(i-1,cols)%vr%vf(3)%fc(:,1))
    enddo
    !$omp end parallel do
  end subroutine wake_continuity

  subroutine age_wake(wake_array,dt)
    type(wakepanel_class), intent(inout), dimension(:,:) :: wake_array
    real(dp),intent(in) :: dt
    integer :: i
    !$omp parallel do
    do i=1,4
      wake_array%vr%vf(i)%age=wake_array%vr%vf(i)%age+dt
    enddo
    !$omp end parallel do
  end subroutine age_wake

  subroutine dissipate_tip(wake_array)
    type(wakepanel_class), intent(inout), dimension(:,:) :: wake_array
    real(dp) :: oseen_param, turb_visc, kin_visc, new_radius
    integer :: ii,tip
    oseen_param= 1.2564_dp
    kin_visc   = 0.0000181_dp
    turb_visc  = 500._dp
    tip=size(wake_array,2)

    do ii=1,size(wake_array,1)
      ! Root vortex core
      new_radius=sqrt(wake_array(ii,1)%vr%vf(1)%r_vc**2._dp &
        +4._dp*oseen_param*turb_visc*kin_visc*wake_array(ii,1)%vr%vf(1)%age)
      wake_array(ii,1)%vr%vf(1)%r_vc=new_radius

      ! Tip vortex core
      new_radius=sqrt(wake_array(ii,tip)%vr%vf(3)%r_vc**2._dp &
        +4._dp*oseen_param*turb_visc*kin_visc*wake_array(ii,tip)%vr%vf(3)%age)
      wake_array(ii,tip)%vr%vf(3)%r_vc=new_radius
    enddo
  end subroutine dissipate_tip

  subroutine strain_wake(wake_array)
    type(wakepanel_class), intent(inout), dimension(:,:) :: wake_array
    integer :: i,j
    !$omp parallel do collapse(2)
    do j=1,size(wake_array,2)
      do i=1,size(wake_array,1)
        call wake_array(i,j)%vr%calclength(.FALSE.)    ! Update current length
        call wake_array(i,j)%vr%strain() 
      enddo
    enddo
    !$omp end parallel do

  end subroutine strain_wake

  !--------------------------------------------------------!
  !                Induced Velocity Functions              !
  !--------------------------------------------------------!

  ! Calculates local velocity at CP velCP and velCPm on wing
  ! Includes uvw, pqr, wake induced velocity
  ! Excludes pitch velocity, wing self-induced velocity
  subroutine vind_CP(wing_array,uvw,pqr,wake_array)
    type(wingpanel_class), intent(inout), dimension(:,:) :: wing_array
    type(wakepanel_class), intent(inout), dimension(:,:) :: wake_array
    real(dp), intent(in), dimension(3) :: uvw, pqr
    integer :: i,j

    do j=1,size(wing_array,2)
      do i=1,size(wing_array,1)
        wing_array(i,j)%velCPm=uvw+cross3(pqr,wing_array(i,j)%cp)
        wing_array(i,j)%velCP=wing_array(i,j)%velCPm+vind_panelgeo(wake_array,wing_array(i,j)%cp)
      enddo
    enddo
  end subroutine vind_CP

  ! Calculates induced vel at P by chordwise vortices of wing_array
  function vind_chordvortex(wing_array,P) result(velind)
    type(wingpanel_class), intent(in), dimension(:,:) :: wing_array
    real(dp), intent(in), dimension(3) :: P
    real(dp), dimension(3) :: velind
    integer :: i,j
    velind=0._dp
    do j=1,size(wing_array,2)
      do i=1,size(wing_array,1)
        velind=velind+wing_array(i,j)%vr%vf(1)%vind(P)*wing_array(i,j)%vr%gam
        velind=velind+wing_array(i,j)%vr%vf(3)%vind(P)*wing_array(i,j)%vr%gam
      enddo
    enddo
  end function vind_chordvortex

  ! Induced velocity by a wing array on point P
  function vind_panelgeo_wing(wing_array,P) result(velind)
    type(wingpanel_class), intent(in), dimension(:,:) :: wing_array
    real(dp), intent(in), dimension(3) :: P
    real(dp), dimension(3,size(wing_array,1),size(wing_array,2)) :: velind_mat
    real(dp), dimension(3) :: velind
    integer :: i,j

    velind_mat=0._dp
    !$omp parallel do collapse(2) shared(wing_array)
    do j=1,size(wing_array,2)
      do i=1,size(wing_array,1)
        velind_mat(:,i,j)=wing_array(i,j)%vr%vind(P)*wing_array(i,j)%vr%gam
      enddo
    enddo
    !$omp end parallel do

    do i=1,3
      velind(i)=sum(velind_mat(i,:,:))
    enddo
  end function vind_panelgeo_wing

  ! Induced velocity by a wake array on point P
  function vind_panelgeo_wake(wake_array,P) result(velind)
    type(wakepanel_class), intent(in), dimension(:,:) :: wake_array
    real(dp), intent(in), dimension(3) :: P
    real(dp), dimension(3,size(wake_array,1),size(wake_array,2)) :: velind_mat
    real(dp), dimension(3) :: velind
    integer :: i,j

    velind_mat=0._dp
    !$omp parallel do collapse(2) shared(wake_array,velind_mat)
    do j=1,size(wake_array,2)
      do i=1,size(wake_array,1)
        velind_mat(:,i,j)=wake_array(i,j)%vr%vind(P)*wake_array(i,j)%vr%gam
      enddo
    enddo
    !$omp end parallel do

    !$omp parallel do
    do i=1,3
      velind(i)=sum(velind_mat(i,:,:))
    enddo
    !$omp end parallel do
  end function vind_panelgeo_wake

  ! Induced velocity by wing_array on wake_array corner points
  function vind_onwake_bywing(wing_array,wake_array) result(vind_array)
    type(wingpanel_class), intent(in), dimension(:,:) :: wing_array
    type(wakepanel_class), intent(in), dimension(:,:) :: wake_array
    real(dp), dimension(3,size(wake_array,1),size(wake_array,2)+1) :: vind_array
    integer :: i,j,rows,cols

    rows=size(wake_array,1)
    cols=size(wake_array,2)

    !$omp parallel do collapse(2) shared(wake_array,wing_array,vind_array)
    do j=1,cols
      do i=1,rows
        vind_array(:,i,j)=vind_panelgeo(wing_array,wake_array(i,j)%vr%vf(2)%fc(:,1))
      enddo
    enddo
    !$omp end parallel do

    !$omp parallel do shared(wake_array,wing_array,vind_array)
    do i=1,rows
      vind_array(:,i,cols+1)=vind_panelgeo(wing_array,wake_array(i,cols)%vr%vf(3)%fc(:,1))
    enddo
    !$omp end parallel do
  end function vind_onwake_bywing

  ! Induced velocity by bywake_array on wake_array corner points
  function vind_onwake_bywake(bywake_array,wake_array) result(vind_array)
    type(wakepanel_class), intent(in), dimension(:,:) :: bywake_array
    type(wakepanel_class), intent(in), dimension(:,:) :: wake_array
    real(dp), dimension(3,size(wake_array,1),size(wake_array,2)+1) :: vind_array
    integer :: i,j,rows,cols

    rows=size(wake_array,1)
    cols=size(wake_array,2)

    !$omp parallel do collapse(2) shared(wake_array,vind_array)
    do j=1,cols
      do i=1,rows
        vind_array(:,i,j)=vind_panelgeo(bywake_array,wake_array(i,j)%vr%vf(2)%fc(:,1))
      enddo
    enddo
    !$omp end parallel do

    !$omp parallel do shared(wake_array,vind_array)
    do i=1,rows
      vind_array(:,i,cols+1)=vind_panelgeo(bywake_array,wake_array(i,cols)%vr%vf(3)%fc(:,1))
    enddo
    !$omp end parallel do
  end function vind_onwake_bywake

  !--------------------------------------------------------!
  !               Force Computation Functions              !
  !--------------------------------------------------------!

  subroutine calc_wingalpha(wing_array)
    type(wingpanel_class), intent(in), dimension(:,:) :: wing_array
    integer :: i,j
    do j=1,size(wing_array,2)
      do i=1,size(wing_array,1)
        call wing_array(i,j)%calc_alpha()
      enddo
    enddo
  end subroutine calc_wingalpha

  function calcgam(wg)
    type(wingpanel_class), intent(inout), dimension(:,:) :: wg  !short form for wing_array
    real(dp), dimension(size(wg,2)) :: calcgam
    integer :: j,rows,cols

    rows=size(wg,1)
    cols=size(wg,2)

    ! Check if this is correct way of calculating sectional circulation
    do j=2,cols
      calcgam(j)=wg(rows,j)%vr%gam
    enddo

  end function calcgam

  function calclift(wg,gamvec_prev,dt)
    type(wingpanel_class), intent(inout), dimension(:,:) :: wg  !short form for wing_array
    real(dp), intent(in), dimension(:) :: gamvec_prev
    real(dp), intent(in) :: dt
    real(dp) :: calclift
    real(dp), dimension(size(wg,1),size(wg,2)) :: gam_prev
    real(dp), dimension(3) :: tau_c, tau_s
    integer :: i,j,rows,cols
    ! Inherent assumption that panels have subdivisions along chord and not inclined to it
    ! while calculating tangent vector
    ! LE and left sides used for calculating tangent vectors

    rows=size(wg,1)
    cols=size(wg,2)

    gam_prev=reshape(gamvec_prev,(/rows,cols/))
    do j=2,cols
      do i=2,rows
        tau_c=wg(i,j)%pc(:,2)-wg(i,j)%pc(:,1)
        tau_s=wg(i,j)%pc(:,4)-wg(i,j)%pc(:,1)
        wg(i,j)%delP=dot_product(wg(i,j)%velCP,tau_c)*(wg(i,j)%vr%gam-wg(i-1,j)%vr%gam)/dot_product(tau_c,tau_c) &
          +          dot_product(wg(i,j)%velCP,tau_s)*(wg(i,j)%vr%gam-wg(i,j-1)%vr%gam)/dot_product(tau_s,tau_s) &
          +          (wg(i,j)%vr%gam-gam_prev(i,j))/dt
      enddo
    enddo

    do j=2,cols
      tau_c=wg(1,j)%pc(:,2)-wg(1,j)%pc(:,1)
      tau_s=wg(1,j)%pc(:,4)-wg(1,j)%pc(:,1)
      wg(1,j)%delP=dot_product(wg(1,j)%velCP,tau_c)*(wg(1,j)%vr%gam)/dot_product(tau_c,tau_c) &
        +          dot_product(wg(1,j)%velCP,tau_s)*(wg(1,j)%vr%gam-wg(1,j-1)%vr%gam)/dot_product(tau_s,tau_s) &
        +          (wg(1,j)%vr%gam-gam_prev(1,j))/dt
    enddo

    tau_c=wg(1,1)%pc(:,2)-wg(1,1)%pc(:,1)
    tau_s=wg(1,1)%pc(:,4)-wg(1,1)%pc(:,1)
    wg(1,1)%delP=dot_product(wg(1,1)%velCP,tau_c)*(wg(1,1)%vr%gam)/dot_product(tau_c,tau_c) &
      +          dot_product(wg(1,1)%velCP,tau_s)*(wg(1,1)%vr%gam)/dot_product(tau_s,tau_s) &
      +          (wg(1,1)%vr%gam-gam_prev(1,1))/dt

    do i=2,rows
      tau_c=wg(i,1)%pc(:,2)-wg(i,1)%pc(:,1)
      tau_s=wg(i,1)%pc(:,4)-wg(i,1)%pc(:,1)
      wg(i,1)%delP=dot_product(wg(i,1)%velCP,tau_c)*(wg(i,1)%vr%gam-wg(i-1,1)%vr%gam)/dot_product(tau_c,tau_c) &
        +          dot_product(wg(i,1)%velCP,tau_s)*(wg(i,1)%vr%gam)/dot_product(tau_s,tau_s) &
        +          (wg(i,1)%vr%gam-gam_prev(i,1))/dt
    enddo
    wg%delP=density*wg%delP

    do j=1,cols
      do i=1,rows
        wg(i,j)%dLift=-(wg(i,j)%delP*wg(i,j)%panel_area)*cos(wg(i,j)%alpha)
      enddo
    enddo

    calclift=0._dp
    do j=1,cols
      do i=1,rows
        calclift=calclift+wg(i,j)%dlift
      enddo
    enddo
  end function calclift

  function calcdrag(wg,gamvec_prev,dt)
    type(wingpanel_class), intent(inout), dimension(:,:) :: wg !short form for wing_array
    real(dp), intent(in), dimension(:) :: gamvec_prev
    real(dp) :: calcdrag
    real(dp), intent(in) :: dt
    real(dp) :: vel_drag
    real(dp) :: drag1, drag2
    real(dp), dimension(size(wg,1),size(wg,2)) :: gam_prev
    integer :: i,j,rows,cols
    ! Inherent assumption that panels have subdivisions along chord and not inclined to it
    ! while calculating tangent vector
    ! LE and left sides used for calculating tangent vectors

    ! *** !! PREDICTS DRAG1 INCORRECTLY !! ***
    ! *** !! PREDICTS DRAG1 INCORRECTLY !! ***
    ! *** !! PREDICTS DRAG1 INCORRECTLY !! ***

    rows=size(wg,1)
    cols=size(wg,2)

    gam_prev=reshape(gamvec_prev,(/rows,cols/))
    do j=1,cols
      do i=2,rows
        !vel_drag=dot_product((vind_panelgeo(wake_array,wg(i,j)%cp))+vind_chordvortex(wg,wg(i,j)%cp),&
        !  (/0._dp,0._dp,1._dp/))
        vel_drag=dot_product((wg(i,j)%velCP-wg(i,j)%velCPm)+vind_chordvortex(wg,wg(i,j)%CP),&
          matmul(wg(i,j)%orthproj(),wg(i,j)%ncap))
        drag2=(wg(i,j)%vr%gam-gam_prev(i,j))*wg(i,j)%panel_area*sin(wg(i,j)%alpha)/dt
        drag1=-vel_drag*(wg(i,j)%vr%gam-wg(i-1,j)%vr%gam)*norm2(wg(i,j)%pc(:,4)-wg(i,j)%pc(:,1))
        wg(i,j)%dDrag=drag1-drag2
      enddo
    enddo

    ! i=1
    do j=2,cols
      !vel_drag=dot_product((vind_panelgeo(wake_array,wg(1,j)%cp))+vind_chordvortex(wg,wg(1,j)%cp),&
      !  (/0._dp,0._dp,1._dp/))
      vel_drag=dot_product((wg(1,j)%velCP-wg(1,j)%velCPm)+vind_chordvortex(wg,wg(1,j)%CP),&
        matmul(wg(1,j)%orthproj(),wg(1,j)%ncap))
      drag2=(wg(1,j)%vr%gam-gam_prev(1,j))*wg(1,j)%panel_area*sin(wg(1,j)%alpha)/dt
      drag1=-vel_drag*(wg(1,j)%vr%gam)*norm2(wg(1,j)%pc(:,4)-wg(1,j)%pc(:,1))
      wg(1,j)%dDrag=drag1-drag2
    enddo

    wg%dDrag=density*wg%dDrag

    calcdrag=0._dp
    do j=1,cols
      do i=1,rows
        calcdrag=calcdrag+wg(i,j)%dDrag
      enddo
    enddo
  end function calcdrag

  !|------+----------------------+------|
  !| ++++ | Bookeeping functions | ++++ |
  !|------+----------------------+------|

  subroutine skiplines(fileunit,nlines)
    integer, intent(in) :: fileunit,nlines
    integer :: i
    do i=1,nlines
      read(fileunit,*)
    enddo
  end subroutine skiplines

end module library
