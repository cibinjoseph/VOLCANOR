module libCommon
  use rotor_classdef
  use switches_classdef
  implicit none
  include "init_file.f90"

contains

  subroutine read_config(filename)
    character(len=*), intent(in) :: filename
    real :: fileFormatVersion, currentVersion

    currentVersion = 0.1

    open (unit=11, file=filename, status='old', action='read')
    call skip_comments(11)
    read (11, *)  fileFormatVersion
    if (abs(fileFormatVersion-currentVersion) < eps) then
      error stop "ERROR: config.in template version does not match"
    endif

    call skip_comments(11)
    read (11, *) nt, dt, nr
    call skip_comments(11)
    read (11, *) switches%restartWriteNt, switches%restartFromNt
    call skip_comments(11)
    read (11, *) switches%ntSub, switches%ntSubInit
    call skip_comments(11)
    read (11, *) switches%spanSpacing
    call skip_comments(11)
    read (11, *) density, velSound
    call skip_comments(11)
    read (11, *) switches%wakePlot, switches%wakeTipPlot, &
      & switches%rotorForcePlot, switches%gridPlot
    call skip_comments(11)
    read (11, *) switches%wakeDissipation, switches%wakeStrain, switches%wakeBurst
    call skip_comments(11)
    read (11, *) switches%slowStart, switches%slowStartNt
    call skip_comments(11)
    read (11, *) switches%fdScheme
    call skip_comments(11)
    read (11, *) switches%initWakeVelNt
    call skip_comments(11)
    read (11, *) switches%probe
    close (11)

  end subroutine read_config

  !--------------------------------------------------------!
  !                Induced Velocity Functions              !
  !--------------------------------------------------------!

  ! Induced velocity by rotor (wing n wake) on Nwake corner points
  function vind_onNwake_byRotor(rotor, Nwake, optionalChar) result(vindArray)
    type(rotor_class), intent(inout) :: rotor
    type(Nwake_class), intent(in), dimension(:, :) :: Nwake
    character(len=1), optional :: optionalChar
    real(dp), dimension(3, size(Nwake, 1), size(Nwake, 2) + 1) :: vindArray
    integer :: i, j, rows, cols

    rows = size(Nwake, 1)
    cols = size(Nwake, 2)

    if (.not. present(optionalChar)) then
      ! Induced velocity due to all blades and wake
      !$omp parallel do collapse(2)
      do j = 1, cols
        do i = 1, rows
          vindArray(:, i, j) = rotor%vind_bywing(Nwake(i, j)%vr%vf(2)%fc(:, 1)) &
            + rotor%vind_bywake(Nwake(i, j)%vr%vf(2)%fc(:, 1))
        enddo
      enddo
      !$omp end parallel do

      !$omp parallel do
      do i = 1, rows
        vindArray(:, i, cols + 1) = rotor%vind_bywing(Nwake(i, cols)%vr%vf(3)%fc(:, 1)) &
          + rotor%vind_bywake(Nwake(i, cols)%vr%vf(3)%fc(:, 1))
      enddo
      !$omp end parallel do

    elseif ((optionalChar .eq. 'P') .or. (optionalChar .eq. 'p')) then

      ! Induced velocity due to all blades and Pwake
      !$omp parallel do collapse(2)
      do j = 1, cols
        do i = 1, rows
          vindArray(:, i, j) = rotor%vind_bywing(Nwake(i, j)%vr%vf(2)%fc(:, 1)) &
                               + rotor%vind_bywake(Nwake(i, j)%vr%vf(2)%fc(:, 1), 'P')
        enddo
      enddo
      !$omp end parallel do

      !$omp parallel do
      do i = 1, rows
        vindArray(:, i, cols + 1) = rotor%vind_bywing(Nwake(i, cols)%vr%vf(3)%fc(:, 1)) &
                                    + rotor%vind_bywake(Nwake(i, cols)%vr%vf(3)%fc(:, 1), 'P')
      enddo
      !$omp end parallel do

    else
      error stop 'ERROR: Wrong character flag for vind_onNwake_byRotor()'
    endif

  end function vind_onNwake_byRotor

  ! Induced velocity by rotor (wing n wake) on Fwake corner points
  function vind_onFwake_byRotor(rotor, Fwake, optionalChar) result(vindArray)
    type(rotor_class), intent(inout) :: rotor
    type(Fwake_class), intent(in), dimension(:) :: Fwake
    character(len=1), optional :: optionalChar
    real(dp), dimension(3, size(Fwake, 1)) :: vindArray
    integer :: i, rows

    rows = size(Fwake, 1)

    if (.not. present(optionalChar)) then
      ! Induced velocity due to all blades and wake
      !$omp parallel do
      do i = 1, rows
        vindArray(:, i) = rotor%vind_bywing(Fwake(i)%vf%fc(:, 1)) &
                          + rotor%vind_bywake(Fwake(i)%vf%fc(:, 1))
      enddo
      !$omp end parallel do

    elseif ((optionalChar .eq. 'P') .or. (optionalChar .eq. 'p')) then

      ! Induced velocity due to all blades and Pwake
      !$omp parallel do
      do i = 1, rows
        vindArray(:, i) = rotor%vind_bywing(Fwake(i)%vf%fc(:, 1)) &
                          + rotor%vind_bywake(Fwake(i)%vf%fc(:, 1), 'P')
      enddo
      !$omp end parallel do

    else
      error stop 'ERROR: Wrong character flag for vind_onFwake_byRotor()'
    endif

  end function vind_onFwake_byRotor

  ! Calculates 2nd order accurate induced velocity on near wake
  function vel_order2_Nwake(v_wake_n, v_wake_np1)   ! np1 => n+1
    real(dp), intent(in), dimension(:, :, :) :: v_wake_n, v_wake_np1
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

  ! Calculates 2nd order accurate induced velocity on far wake
  function vel_order2_Fwake(v_wake_n, v_wake_np1)   ! np1 => n+1
    real(dp), intent(in), dimension(:, :) :: v_wake_n, v_wake_np1
    real(dp), dimension(3, size(v_wake_n, 2)) :: vel_order2_Fwake
    integer :: i
    vel_order2_Fwake(:, 1) = (v_wake_np1(:, 1) + v_wake_n(:, 1))*0.5_dp
    !$omp parallel do
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

  ! Prints status message (or SUCCESS if left blank)
  subroutine print_status(statusMessage)
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
