! Contains variable declarations
type(rotor_class), allocatable, dimension(:) :: rotor

! Kinematics
integer :: nt, nr
real(dp) :: dt

! Other Variables
real(dp) :: t
real(dp) :: density, velSound
character(len=5) :: timestamp
character(len=2) :: rotorChar
character(len=10) :: rotorFile
character(len=10) :: currentTime
logical :: fileExists

! Iterators
integer :: i, is, ic, row, iter, ir, jr, ib

! Switches
integer :: ntSub, ntSubInit
integer :: spanSpacingSwitch
integer :: wakeDissipationSwitch, wakeStrainSwitch, wakeBurstSwitch
integer :: slowStartSwitch, slowStartNt
integer :: wakeTipPlotSwitch, wakePlotSwitch, gridPlotSwitch
integer :: rotorForcePlotSwitch
integer :: fdSchemeSwitch
integer :: wakeIgnoreNt
integer :: initWakeVelNt
