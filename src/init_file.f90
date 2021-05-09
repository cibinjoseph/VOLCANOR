! Contains variable declarations
type(rotor_class), allocatable, dimension(:) :: rotor

! Kinematics
integer :: nt, nr
real(dp) :: dt

! Other Variables
integer :: iterStart
real(dp) :: t
real(dp) :: density, velSound
character(len=5) :: timestamp
character(len=2) :: rotorChar
character(len=10) :: rotorFile
character(len=10) :: currentTime
logical :: fileExists
real(dp) :: subIterResidual

! Iterators
integer :: i, is, ic, row, iter, ir, jr, ib, rowErase

! Switches
type(switches_class) :: switches

! Probes for velocity
real(dp), allocatable, dimension(:, :) :: probe, probeVel
