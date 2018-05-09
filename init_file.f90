! Contains variable declarations
type(rotor_class), allocatable, dimension(:) :: rotor
!real(dp), dimension(nc*ns,nc*ns) :: Amat, Amat_inv
type(wakepanel_class), allocatable, dimension(:,:) :: Pwake
real(dp), dimension(nc*ns) :: RHS, gamvec, gamvec_prev
real(dp), dimension(ns) :: gam_sectional

! Wing variables
real(dp) :: chord, span, root_cut

! Kinematics
integer :: nt,nr
real(dp) :: dt, om_theta, om_h
real(dp) :: theta_pitch, dtheta_pitch
real(dp) :: theta0, thetac, thetas, h0, hdot, thetadot
real(dp), dimension(3) :: vwind, om_body, pqr, vbody, om_body_slow

! Other Variables
real(dp), dimension(3) :: v_shed, vec_dummy, vel_plunge
real(dp) :: t, init_wake_vel, starting_vortex_core
real(dp), dimension(nt) :: lift, drag
real(dp), dimension(3) :: pts, dpts   ! phi, theta, psi
real(dp), dimension(3,nt,ns+1) :: vind_wake
real(dp), allocatable, dimension(:,:,:) :: vind_wake1, vind_wake2, vind_wake3
real(dp), allocatable, dimension(:,:,:) :: Pvind_wake, vind_wake_step
character(len=5) :: timestamp
character(len=2) :: rotor_char
character(len=10) :: rotorfile
real(dp) :: wing_mid_core,wake_mid_core,wing_tip_core,wake_tip_core
real(dp) :: pivotLE

! Iterators
integer :: ispan, ichord, row, col, i, j, indx, iter, row_now, irotor


! Allocate vars required for wake convection
! on the basis of finite diff scheme
select case (FDscheme_switch)
case (1)
  allocate(Pwake(nt,ns))
  allocate(vind_wake1(3,nt,ns+1))
  allocate(Pvind_wake(3,nt,ns+1))
case (2)
  allocate(vind_wake1(3,nt,ns+1))
  allocate(vind_wake_step(3,nt,ns+1))
case (3)
  allocate(Pwake(nt,ns))
  allocate(vind_wake1(3,nt,ns+1))
  allocate(vind_wake2(3,nt,ns+1))
  allocate(vind_wake3(3,nt,ns+1))
  allocate(Pvind_wake(3,nt,ns+1))
  allocate(vind_wake_step(3,nt,ns+1))
end select

