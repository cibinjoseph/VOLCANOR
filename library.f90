module library
  use wingpanel_classdef
  use wakepanel_classdef
  implicit none

  ! Input parameters
  integer, parameter  :: nt = 100
  integer, parameter  :: ns = 8
  integer, parameter  :: nc = 3

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
  subroutine init_wing(wing_array,xvec,yvec,core_radius)
    type(wingpanel_class), intent(out), dimension(:,:) :: wing_array
    real(dp), intent(in), dimension(:) :: xvec 
    real(dp), intent(in), dimension(:) :: yvec
    real(dp), intent(in) :: core_radius
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

    ! Initialize core radius
    do i=1,4
      wing_array%vr%vf(i)%r_vc0=core_radius
      wing_array%vr%vf(i)%r_vc =core_radius
      wing_array%vr%vf(i)%age=0._dp
    enddo

  end subroutine init_wing

  ! Assigns vortex code radii to all filaments
  subroutine init_wake(wake_array,core_radius)
    type(wakepanel_class), intent(out), dimension(:,:) :: wake_array
    real(dp), intent(in) :: core_radius
    integer :: i

    do i=1,4
      wake_array%vr%vf(i)%r_vc0=core_radius
      wake_array%vr%vf(i)%r_vc =core_radius
      wake_array%vr%vf(i)%age=0._dp
    enddo
    wake_array%tag=-1
    wake_array%vr%gam=0._dp

  end subroutine init_wake

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

  subroutine rot_wing(wing_array,pts,order)  
    type(wingpanel_class), intent(inout), dimension(:,:) :: wing_array
    real(dp), dimension(3), intent(in) :: pts  ! pts => phi,theta,psi
    integer :: i, j
    integer :: order    ! [1]gb & +ve theta , [2]bg & -ve theta       
    real(dp), dimension(3,3) :: TMat

    select case (order)
    case (2)
      TMat=Tbg((/cos(pts(1)),sin(pts(1))/),&
        (/cos(pts(2)),sin(pts(2))/),&
        (/cos(pts(3)),sin(pts(3))/))
    case (1)
      TMat=Tgb((/cos(pts(1)),sin(pts(1))/),&
        (/cos(pts(2)),sin(pts(2))/),&
        (/cos(pts(3)),sin(pts(3))/))
    case default
      error stop 'Error: wrong option for order'
    end select

    do j=1,ns
      do i=1,nc
        call wing_array(i,j)%rot(TMat)
      enddo
    enddo

  end subroutine rot_wing

  subroutine pitch_wing(wing_array,theta_pitch,pts)  !pitch about a fixed point fp
    ! assuming motion is 2dimensional and body has not undergone 3d rotations
    type(wingpanel_class), intent(inout), dimension(:,:) :: wing_array
    real(dp), intent(in) :: theta_pitch
    real(dp), dimension(3), intent(in) :: pts
    real(dp), dimension(3) :: dshift

    if (theta_pitch>eps) then
      ! Translate to origin
      dshift=(/wing_array(1,1)%pc(:,1)/)
      call mov_wing(wing_array,-dshift)

      ! Rotate global angles
      call rot_wing(wing_array,pts,1)

      ! Rotate pitch angle
      call rot_wing(wing_array,(/0._dp,theta_pitch,0._dp/),1)

      ! Unrotate global angles
      call rot_wing(wing_array,-1._dp*pts,1)

      ! Untranslate from origin
      call mov_wing(wing_array,dshift)
    endif

  end subroutine pitch_wing

  subroutine mov_wing(wing_array,dshift)
    type(wingpanel_class), intent(inout), dimension(:,:) :: wing_array
    real(dp), intent(in), dimension(3) :: dshift
    integer :: i,j

    do j=1,size(wing_array,2)
      do i=1,size(wing_array,1)
        call wing_array(i,j)%shiftdP(dshift)
      enddo
    enddo
  end subroutine mov_wing

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

    do j=1,cols
      do i=1,rows
        call wake_array(i,j)%vr%shiftdP(2,dP_array(:,i,j))
      enddo
    enddo

    do i=1,rows
      call wake_array(i,cols)%vr%shiftdP(3,dP_array(:,i,cols+1))
    enddo
  end subroutine convectwake

  subroutine wake_continuity(wake_array)
    type(wakepanel_class), intent(inout), dimension(:,:) :: wake_array
    integer :: i,j,rows,cols

    rows=size(wake_array,1)
    cols=size(wake_array,2)

    do j=1,cols-1
      do i=2,rows
        call wake_array(i,j)%vr%assignP(1,wake_array(i-1,j)%vr%vf(2)%fc(:,1))
        call wake_array(i,j)%vr%assignP(3,wake_array(i,j+1)%vr%vf(2)%fc(:,1))
        call wake_array(i,j)%vr%assignP(4,wake_array(i-1,j+1)%vr%vf(2)%fc(:,1))
      enddo
    enddo

    do j=1,cols-1
      call wake_array(1,j)%vr%assignP(3,wake_array(1,j+1)%vr%vf(2)%fc(:,1))
    enddo

    do i=2,rows
      call wake_array(i,cols)%vr%assignP(1,wake_array(i-1,cols)%vr%vf(2)%fc(:,1))
      call wake_array(i,cols)%vr%assignP(4,wake_array(i-1,cols)%vr%vf(3)%fc(:,1))
    enddo
  end subroutine wake_continuity

  subroutine age_wake(wake_array,dt)
    type(wakepanel_class), intent(inout), dimension(:,:) :: wake_array
    real(dp),intent(in) :: dt
    integer :: i
    do i=1,4
      wake_array%vr%vf(i)%age=wake_array%vr%vf(i)%age+dt
    enddo
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


  !--------------------------------------------------------!
  !                Induced Velocity Functions              !
  !--------------------------------------------------------!

  ! Calculates local velocity at CP on wing
  ! Excluding Pitch velocity and wing induced velocity
  subroutine vind_CP(wing_array,uvw,pqr,wake_array)
    type(wingpanel_class), intent(inout), dimension(:,:) :: wing_array
    type(wakepanel_class), intent(inout), dimension(:,:) :: wake_array
    real(dp), intent(in), dimension(3) :: uvw, pqr
    integer :: i,j

    do j=1,size(wing_array,2)
      do i=1,size(wing_array,1)
        wing_array(i,j)%velCP=uvw
        wing_array(i,j)%velCP=wing_array(i,j)%velCP+vind_panelgeo(wake_array,wing_array(i,j)%cp)
        wing_array(i,j)%velCP=wing_array(i,j)%velCP+cross3(pqr,wing_array(i,j)%cp)
      enddo
    enddo
  end subroutine vind_CP

  ! Induced velocity by a wing array on point P
  function vind_panelgeo_wing(wing_array,P) result(velind)
    type(wingpanel_class), intent(in), dimension(:,:) :: wing_array
    real(dp), intent(in), dimension(3) :: P
    real(dp), dimension(3,size(wing_array,1),size(wing_array,2)) :: velind_mat
    real(dp), dimension(3) :: velind
    integer :: i,j

    velind_mat=0._dp
    !$omp parallel do collapse(2)
    do j=1,size(wing_array,2)
      do i=1,size(wing_array,1)
        velind_mat(:,i,j)=wing_array(i,j)%vr%vind(P)*wing_array(i,j)%vr%gam
      enddo
    enddo
    !$omp end parallel do
    velind(1)=sum(velind_mat(1,:,:))
    velind(2)=sum(velind_mat(2,:,:))
    velind(3)=sum(velind_mat(3,:,:))
  end function vind_panelgeo_wing

  ! Induced velocity by a wake array on point P
  function vind_panelgeo_wake(wake_array,P) result(velind)
    type(wakepanel_class), intent(in), dimension(:,:) :: wake_array
    real(dp), intent(in), dimension(3) :: P
    real(dp), dimension(3,size(wake_array,1),size(wake_array,2)) :: velind_mat
    real(dp), dimension(3) :: velind
    integer :: i,j

    velind_mat=0._dp
    !$omp parallel do collapse(2)
    do j=1,size(wake_array,2)
      do i=1,size(wake_array,1)
        velind_mat(:,i,j)=wake_array(i,j)%vr%vind(P)*wake_array(i,j)%vr%gam
      enddo
    enddo
    !$omp end parallel do
    velind(1)=sum(velind_mat(1,:,:))
    velind(2)=sum(velind_mat(2,:,:))
    velind(3)=sum(velind_mat(3,:,:))
  end function vind_panelgeo_wake

  ! Induced velocity by wing_array on wake_array corner points
  function vind_onwake_bywing(wing_array,wake_array) result(vind_array)
    type(wingpanel_class), intent(in), dimension(:,:) :: wing_array
    type(wakepanel_class), intent(in), dimension(:,:) :: wake_array
    real(dp), dimension(3,size(wake_array,1),size(wake_array,2)+1) :: vind_array
    integer :: i,j,rows,cols

    rows=size(wake_array,1)
    cols=size(wake_array,2)

    do j=1,cols
      do i=1,rows
        vind_array(:,i,j)=vind_panelgeo(wing_array,wake_array(i,j)%vr%vf(2)%fc(:,1))
      enddo
    enddo

    do i=1,rows
      vind_array(:,i,cols+1)=vind_panelgeo(wing_array,wake_array(i,cols)%vr%vf(3)%fc(:,1))
    enddo
  end function vind_onwake_bywing

  ! Induced velocity by bywake_array on wake_array corner points
  function vind_onwake_bywake(bywake_array,wake_array) result(vind_array)
    type(wakepanel_class), intent(in), dimension(:,:) :: bywake_array
    type(wakepanel_class), intent(in), dimension(:,:) :: wake_array
    real(dp), dimension(3,size(wake_array,1),size(wake_array,2)+1) :: vind_array
    integer :: i,j,rows,cols

    rows=size(wake_array,1)
    cols=size(wake_array,2)

    do j=1,cols
      do i=1,rows
        vind_array(:,i,j)=vind_panelgeo(bywake_array,wake_array(i,j)%vr%vf(2)%fc(:,1))
      enddo
    enddo

    do i=1,rows
      vind_array(:,i,cols+1)=vind_panelgeo(bywake_array,wake_array(i,cols)%vr%vf(3)%fc(:,1))
    enddo
  end function vind_onwake_bywake

  !--------------------------------------------------------!
  !                Lift Computation Functions              !
  !--------------------------------------------------------!

  subroutine calclift(wg,gamvec_prev,dt)
    type(wingpanel_class), intent(inout), dimension(:,:) :: wg !short form for wing_array
    real(dp), intent(in), dimension(:) :: gamvec_prev
    real(dp), intent(in) :: dt
    real(dp) :: density
    real(dp), dimension(size(wg,1),size(wg,2)) :: gam_prev
    real(dp), dimension(3) :: tau_i, tau_j
    integer :: i,j,rows,cols
    ! Inherent assumption that panels have subdivisions along chord and not inclined to it
    ! while calculating tangent vector
    ! LE and left sides used for calculating tangent vectors
    density=1.2_dp

    rows=size(wg,1)
    cols=size(wg,2)

    gam_prev=reshape(gamvec_prev,(/rows,cols/))
    do j=2,cols
      do i=2,rows
        tau_i=wg(i,j)%pc(:,2)-wg(i,j)%pc(:,1)
        tau_j=wg(i,j)%pc(:,4)-wg(i,j)%pc(:,1)
        ! Adding self induced velocity at collocation point
        wg(i,j)%velCP=wg(i,j)%velCP+vind_panelgeo_wing(wg,wg(i,j)%cp)
        wg(i,j)%delP=dot_product(wg(i,j)%velCP,tau_i)*(wg(i,j)%vr%gam-wg(i-1,j)%vr%gam)/dot_product(tau_i,tau_i) &
          +          dot_product(wg(i,j)%velCP,tau_j)*(wg(i,j)%vr%gam-wg(i,j-1)%vr%gam)/dot_product(tau_j,tau_j) &
          +          (wg(i,j)%vr%gam-gam_prev(i,j))/dt
      enddo
    enddo

    do j=1,cols
      tau_i=wg(1,j)%pc(:,2)-wg(1,j)%pc(:,1)
      tau_j=wg(1,j)%pc(:,4)-wg(1,j)%pc(:,1)
      ! Adding self induced velocity at collocation point
      wg(1,j)%velCP=wg(1,j)%velCP+vind_panelgeo_wing(wg,wg(1,j)%cp)
      wg(1,j)%delP=dot_product(wg(1,j)%velCP,tau_i)*(wg(1,j)%vr%gam)/dot_product(tau_i,tau_i) &
        +          dot_product(wg(1,j)%velCP,tau_j)*(wg(1,j)%vr%gam)/dot_product(tau_j,tau_j) &
        +          (wg(1,j)%vr%gam-gam_prev(1,j))/dt
    enddo

    do i=2,rows
      tau_i=wg(i,1)%pc(:,2)-wg(i,1)%pc(:,1)
      tau_j=wg(i,1)%pc(:,4)-wg(i,1)%pc(:,1)
      ! Adding self induced velocity at collocation point
      wg(i,1)%velCP=wg(i,1)%velCP+vind_panelgeo_wing(wg,wg(i,1)%cp)
      wg(i,1)%delP=dot_product(wg(i,1)%velCP,tau_i)*(wg(i,1)%vr%gam)/dot_product(tau_i,tau_i) &
        +          dot_product(wg(i,1)%velCP,tau_j)*(wg(i,1)%vr%gam)/dot_product(tau_j,tau_j) &
        +          (wg(i,1)%vr%gam-gam_prev(i,1))/dt
    enddo
    wg%delP=density*wg%delP

    do j=1,cols
      do i=1,rows
        wg(i,j)%dLift=-(wg(i,j)%delP*wg(i,j)%panel_area)*wg(i,j)%ncap
      enddo
    enddo
  end subroutine calclift

end module library
