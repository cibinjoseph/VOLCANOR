! Contains variable declarations
type(rotor_class), allocatable, dimension(:) :: rotor
!type(wakepanel_class), allocatable, dimension(:,:) :: Pwake
!real(dp), dimension(ns) :: gam_sectional

! Kinematics
integer :: nt,nr
real(dp) :: dt
real(dp) :: theta_pitch, dtheta_pitch
!real(dp), dimension(3) :: om_body_slow

! Other Variables
real(dp) :: t
real(dp) :: density
!real(dp), dimension(nt) :: lift, drag
!real(dp), dimension(3,nt,ns+1) :: vind_wake
!real(dp), allocatable, dimension(:,:,:) :: vind_wake1, vind_wake2, vind_wake3
!real(dp), allocatable, dimension(:,:,:) :: Pvind_wake, vind_wake_step
character(len=5) :: timestamp
character(len=2) :: rotor_char
character(len=10) :: rotorfile

! Iterators
integer :: is,ic,row,col,i,j,iter,row_now,ir,jr,ib,jb

! Switches
integer :: span_spacing_switch
integer :: tip_diss_switch, wakestrain_switch
integer :: slowstart_switch, slowstart_nt
integer :: wakeplot_switch
integer :: FDscheme_switch
integer :: wake_ignore_nt
integer :: init_wake_vel_nt

! Allocate vars required for wake convection
! on the basis of finite diff scheme
!select case (FDscheme_switch)
!case (1)
!  allocate(Pwake(nt,ns))
!  allocate(vind_wake1(3,nt,ns+1))
!  allocate(Pvind_wake(3,nt,ns+1))
!case (2)
!  allocate(vind_wake1(3,nt,ns+1))
!  allocate(vind_wake_step(3,nt,ns+1))
!case (3)
!  allocate(Pwake(nt,ns))
!  allocate(vind_wake1(3,nt,ns+1))
!  allocate(vind_wake2(3,nt,ns+1))
!  allocate(vind_wake3(3,nt,ns+1))
!  allocate(Pvind_wake(3,nt,ns+1))
!  allocate(vind_wake_step(3,nt,ns+1))
!end select

