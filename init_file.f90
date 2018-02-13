! Contains variable declarations
real(dp), dimension(nc+1)   :: xvec
real(dp), dimension(ns+1)     :: yvec
type(wingpanel_class), dimension(nc,ns) :: wing
type(wakepanel_class), dimension(nt,ns) :: wake
type(wakepanel_class), dimension(nt,ns) :: Pwake
real(dp), dimension(nc*ns,nc*ns) :: Amat, Amat_inv
real(dp), dimension(nc*ns) :: RHS, gamvec, gamvec_prev
real(dp), dimension(3) :: hub_coords

! Wing variables
real(dp) :: chord, span

! Kinematics
real(dp) :: dt, om_theta, om_h
real(dp) :: theta_pitch, dtheta_pitch
real(dp) :: theta0, thetac, thetas, h0, hdot, thetadot
real(dp), dimension(3) :: vwind, om_body, pqr, vbody, om_body_slow

! Other Variables
real(dp), dimension(3) :: v_shed, vec_dummy, vel_plunge
real(dp) :: t, init_wake_vel
real(dp), dimension(nt) :: lift, drag
real(dp), dimension(3) :: pts, dpts   ! phi, theta, psi
real(dp), dimension(3,nt,ns+1) :: vind_wake, vind_wake1, vind_wake2, vind_wake3
real(dp), dimension(3,nt,ns+1) :: Pvind_wake, vind_wake_step
character(len=5) :: timestamp
real(dp) :: wing_mid_core,wake_mid_core,wing_tip_core,wake_tip_core

! Iterators
integer :: ispan, ichord, row, col, i, j, indx, iter, row_now
