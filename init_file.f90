! Contains variable declarations
real(dp), dimension(nc+1)   :: xvec
real(dp), dimension(ns+1)     :: yvec
type(wingpanel_class), dimension(nc,ns) :: wing
type(wakepanel_class), dimension(nt,ns) :: wake
real(dp), dimension(nc*ns,nc*ns) :: Amat, Amat_inv
real(dp), dimension(nc*ns) :: RHS, gamvec, gamvec_prev

! Wing variables
real(dp) :: chord, span

! Kinematics
real(dp) :: dt, om_theta, om_h
real(dp) :: theta_pitch, dtheta_pitch
real(dp) :: theta0, thetac, thetas, h0, hdot, thetadot
real(dp), dimension(3) :: vwind, om_body, pqr, vbody, om_body_slow

! Other Variables
real(dp), dimension(3) :: v_shed, vec_dummy, vel_plunge
real(dp) :: vel_pitch
real(dp) :: t, liftQS
real(dp), dimension(3) :: lift
real(dp), dimension(3) :: pts, dpts   ! phi, theta, psi
real(dp), dimension(3,nt,ns+1) :: vind_wake
character(len=5) :: timestamp

! Iterators
integer :: ispan, ichord, row, col, i, j, indx, iter, row_now
