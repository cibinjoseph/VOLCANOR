! Contains variable declarations
type(rotor_class), allocatable, dimension(:) :: rotor
!real(dp), dimension(ns) :: gam_sectional

! Kinematics
integer :: nt,nr
real(dp) :: dt

! Other Variables
real(dp) :: t
real(dp) :: density
!real(dp), dimension(nt) :: lift, drag
character(len=5) :: timestamp
character(len=2) :: rotor_char
character(len=10) :: rotorfile

! Iterators
integer :: i,is,ic,row,iter,ir,jr,ib

! Switches
integer :: span_spacing_switch
integer :: tip_diss_switch, wakestrain_switch
integer :: slowstart_switch, slowstart_nt
integer :: wakeplot_switch
integer :: FDscheme_switch
integer :: wake_ignore_nt
integer :: init_wake_vel_nt
