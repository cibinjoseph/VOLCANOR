!! Module definition for libCommon

module libCommon
  !! Utility procedures and those to compute induced velocity of rotors

  use classdef, only : switches_class, rotor_class
  implicit none

  integer, parameter, private :: dp = kind(1.d0)
  !! Double precision setting

  ! Contains variable declarations
  type(rotor_class), allocatable, dimension(:) :: rotor

  ! Kinematics
  integer :: nt
  !! No. of timesteps
  integer :: nr
  !! No. of rotors
  real(dp) :: dt
  !! Timestep size

  ! Other Variables
  integer :: iterStart
  real(dp) :: t
  !! Time
  real(dp) :: density
  !! Density
  real(dp) :: velSound
  !! Velocity of sound to compute Mach
  real(dp) :: kinematicVisc
  !! Kinematic viscosity
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

contains

  subroutine readConfig(filename, outputFilename)
    !! Read config.nml input file (in namelist format)

    use classdef, only : switches_class
    character(len=*), intent(in) :: filename
    character(len=*), optional, intent(in) :: outputFilename
    character(len=10) :: fileFormatVersion, currentVersion
    integer :: restartWriteNt, restartFromNt, ntSub, ntSubInit, &
      & wakePlot, wakeTipPlot, rotorForcePlot, &
      & gridPlot, wakeDissipation, wakeStrain, wakeBurst, wakeSuppress, &
      & slowStart, slowStartNt, fdScheme, initWakeVelNt, probe

    ! Namelists
    namelist /VERSION/ fileFormatVersion
    namelist /PARAMS/ nt, dt, nr, density, velSound, kinematicVisc
    namelist /OPTIONS/ restartWriteNt, restartFromNt, ntSub, ntSubInit, &
      & wakePlot, wakeTipPlot, rotorForcePlot, &
      & gridPlot, wakeDissipation, wakeStrain, wakeBurst, wakeSuppress, &
      & slowStart, slowStartNt, fdScheme, initWakeVelNt, probe

    currentVersion = '0.5'

    open(unit=11, file=filename, status='old', action='read')
    read(unit=11, nml=VERSION)
    if (adjustl(fileFormatVersion) /= currentVersion) then
      error stop "ERROR: config.nml template version does not match"
    endif
    read(unit=11, nml=PARAMS)
    read(unit=11, nml=OPTIONS)
    close(11)

    ! Write a copy to results if requested
    if (present(outputFilename)) then
      open(unit=12, file=outputFilename, status='replace', action='write')
      write(unit=12, nml=VERSION)
      write(unit=12, nml=PARAMS)
      write(unit=12, nml=OPTIONS)
      close(12)
    endif

    switches%restartWriteNt = restartWriteNt
    switches%restartFromNt = restartFromNt
    switches%ntSub = ntSub
    switches%ntSubInit = ntSubInit
    switches%wakePlot = wakePlot
    switches%wakeTipPlot = wakeTipPlot
    switches%rotorForcePlot = rotorForcePlot
    switches%gridPlot = gridPlot
    switches%wakeDissipation = wakeDissipation
    switches%wakeStrain = wakeStrain
    switches%wakeBurst = wakeBurst
    switches%wakeSuppress = wakeSuppress
    switches%slowStart = slowStart
    switches%slowStartNt = slowStartNt
    switches%fdScheme = fdScheme
    switches%initWakeVelNt = initWakeVelNt
    switches%probe = probe
  end subroutine readConfig

  !--------------------------------------------------------!
  !                Induced Velocity Functions              !
  !--------------------------------------------------------!

  function vind_onNwake_byRotor(rotor, Nwake, optionalChar) result(vindArray)
    !! Compute induced velocity by rotor (wing + wake) on Nwake corner points

    use classdef, only : rotor_class, Nwake_class
    type(rotor_class), intent(inout) :: rotor
    !! Rotor
    type(Nwake_class), intent(in), dimension(:, :) :: Nwake
    !! Near wake
    character(len=1), optional :: optionalChar
    !! If 'P' is specified, predicted wake of rotor is used
    real(dp), dimension(3, size(Nwake, 1), size(Nwake, 2) + 1) :: vindArray
    integer :: i, j, rows, cols

    rows = size(Nwake, 1)
    cols = size(Nwake, 2)

    if (.not. present(optionalChar)) then
      ! Induced velocity due to all blades and wake
      !$omp parallel do collapse(2) schedule(runtime)
      do j = 1, cols
        do i = 1, rows
          vindArray(:, i, j) = rotor%vind_bywing(Nwake(i, j)%vr%vf(2)%fc(:, 1)) &
            + rotor%vind_bywake(Nwake(i, j)%vr%vf(2)%fc(:, 1))
        enddo
      enddo
      !$omp end parallel do

      !$omp parallel do schedule(runtime)
      do i = 1, rows
        vindArray(:, i, cols + 1) = rotor%vind_bywing(Nwake(i, cols)%vr%vf(3)%fc(:, 1)) &
          + rotor%vind_bywake(Nwake(i, cols)%vr%vf(3)%fc(:, 1))
      enddo
      !$omp end parallel do

    elseif ((optionalChar .eq. 'P') .or. (optionalChar .eq. 'p')) then

      ! Induced velocity due to all blades and Pwake
      !$omp parallel do collapse(2) schedule(runtime)
      do j = 1, cols
        do i = 1, rows
          vindArray(:, i, j) = rotor%vind_bywing(Nwake(i, j)%vr%vf(2)%fc(:, 1)) &
                               + rotor%vind_bywake(Nwake(i, j)%vr%vf(2)%fc(:, 1), 'P')
        enddo
      enddo
      !$omp end parallel do

      !$omp parallel do schedule(runtime)
      do i = 1, rows
        vindArray(:, i, cols + 1) = rotor%vind_bywing(Nwake(i, cols)%vr%vf(3)%fc(:, 1)) &
          + rotor%vind_bywake(Nwake(i, cols)%vr%vf(3)%fc(:, 1), 'P')
      enddo
      !$omp end parallel do

    else
      error stop 'ERROR: Wrong character flag for vind_onNwake_byRotor()'
    endif

  end function vind_onNwake_byRotor

  function vind_onFwake_byRotor(rotor, Fwake, optionalChar) result(vindArray)
    !! Compute induced velocity by rotor (wing + wake) on Fwake corner points

    use classdef, only : rotor_class, Fwake_class
    type(rotor_class), intent(inout) :: rotor
    !! Rotor
    type(Fwake_class), intent(in), dimension(:) :: Fwake
    !! Far wake
    character(len=1), optional :: optionalChar
    !! If 'P' is specified, predicted wake of rotor is used
    real(dp), dimension(3, size(Fwake, 1)) :: vindArray
    integer :: i, rows

    rows = size(Fwake, 1)

    if (.not. present(optionalChar)) then
      ! Induced velocity due to all blades and wake
      !$omp parallel do schedule(runtime)
      do i = 1, rows
        vindArray(:, i) = rotor%vind_bywing(Fwake(i)%vf%fc(:, 1)) &
          + rotor%vind_bywake(Fwake(i)%vf%fc(:, 1))
      enddo
      !$omp end parallel do

    elseif ((optionalChar .eq. 'P') .or. (optionalChar .eq. 'p')) then

      ! Induced velocity due to all blades and Pwake
      !$omp parallel do schedule(runtime)
      do i = 1, rows
        vindArray(:, i) = rotor%vind_bywing(Fwake(i)%vf%fc(:, 1)) &
          + rotor%vind_bywake(Fwake(i)%vf%fc(:, 1), 'P')
      enddo
      !$omp end parallel do

    else
      error stop 'ERROR: Wrong character flag for vind_onFwake_byRotor()'
    endif

  end function vind_onFwake_byRotor

  function vel_order2_Nwake(v_wake_n, v_wake_np1)
    !! Calculate 2nd order accurate induced velocity on near wake

    real(dp), intent(in), dimension(:, :, :) :: v_wake_n
    !! Induced velocity on wake at timestep n
    real(dp), intent(in), dimension(:, :, :) :: v_wake_np1
    !! Induced velocity on wake at timestep n+1
    real(dp), dimension(3, size(v_wake_n, 2), size(v_wake_n, 3)) :: vel_order2_Nwake
    integer :: i, j
    do j = 1, size(v_wake_n, 3)
      vel_order2_Nwake(:, 1, j) = (v_wake_np1(:, 1, j) &
        + v_wake_n(:, 1, j))*0.5_dp
      do i = 2, size(v_wake_n, 2) - 1
        vel_order2_Nwake(:, i, j) = (v_wake_np1(:, i, j) &
          + v_wake_np1(:, i - 1, j) &
          + v_wake_n(:, i + 1, j) &
          + v_wake_n(:, i, j))*0.25_dp
      enddo
      vel_order2_Nwake(:, size(v_wake_n, 2), j) = &
        (v_wake_np1(:, size(v_wake_n, 2), j) &
        + v_wake_n(:, size(v_wake_n, 2), j))*0.5_dp
    enddo
  end function vel_order2_Nwake

  function vel_order2_Fwake(v_wake_n, v_wake_np1)
    !! Calculate 2nd order accurate induced velocity on far wake

    real(dp), intent(in), dimension(:, :) :: v_wake_n
    !! Induced velocity on wake at timestep n
    real(dp), intent(in), dimension(:, :) :: v_wake_np1
    !! Induced velocity on wake at timestep n+1
    real(dp), dimension(3, size(v_wake_n, 2)) :: vel_order2_Fwake
    integer :: i
    vel_order2_Fwake(:, 1) = (v_wake_np1(:, 1) + v_wake_n(:, 1))*0.5_dp
    !$omp parallel do schedule(runtime)
    do i = 2, size(v_wake_n, 2) - 1
      vel_order2_Fwake(:, i) = (v_wake_np1(:, i) &
        + v_wake_np1(:, i - 1) &
        + v_wake_n(:, i + 1) &
        + v_wake_n(:, i))*0.25_dp
    enddo
    !$omp end parallel do
    vel_order2_Fwake(:, size(v_wake_n, 2)) = &
      (v_wake_np1(:, size(v_wake_n, 2)) &
      + v_wake_n(:, size(v_wake_n, 2)))*0.5_dp
  end function vel_order2_Fwake

  subroutine print_status(statusMessage)
    !! Prints status message (or SUCCESS if left blank)

    character(len=*), optional :: statusMessage
    character(len=34) :: statusPrint    ! Adjust for spacing

    if (.not. present(statusMessage)) then
      write (*, '(A)') '...   SUCCESS'
    else
      statusPrint = statusMessage
      write (*, '(A)', advance='no') statusPrint
    endif
  end subroutine print_status

end module libCommon
