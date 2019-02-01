! Contains variable declarations
type(rotor_class), allocatable, dimension(:) :: rotor

! Kinematics
integer :: nt,nr
real(dp) :: dt

! Other Variables
real(dp) :: t
real(dp) :: density,turbulentViscosity
character(len=5) :: timestamp
character(len=2) :: rotorChar
character(len=10) :: rotorFile

! Iterators
integer :: i,is,ic,row,iter,ir,jr,ib

! Switches
integer :: spanSpacingSwitch
integer :: tipDissipationSwitch, wakeStrainSwitch
integer :: slowStartSwitch, slowStartNt
integer :: wakePlotSwitch, gridPlotSwitch
integer :: forcePlotSwitch, forceCalcSwitch
integer :: fdSchemeSwitch
integer :: wakeIgnoreNt
integer :: initWakeVelNt
