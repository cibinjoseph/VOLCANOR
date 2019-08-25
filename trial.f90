program trial
  use vr_classdef
  implicit none 
  type(vr_class) :: vr
  type(vf_class) :: vfx, vfy

  ! Rectangular
  call vr%assignP(1,(/0.0000_dp,0.0000_dp,0.0_dp/))
  call vr%assignP(2,(/0.3048_dp,0.0000_dp,0.0_dp/))
  call vr%assignP(3,(/0.3048_dp,0.3048_dp,0.0_dp/))
  call vr%assignP(4,(/0.3048_dp,0.0000_dp,0.0_dp/))
  P = (/0.2286_dp,0.1524_dp,0.0_dp/)

  ! Swept
  call vr%assignP(1,(/0.0000_dp,0.0000_dp,0.0_dp/))
  call vr%assignP(2,(/0.3048_dp,0.0000_dp,0.0_dp/))
  call vr%assignP(3,(/0.3048_dp,0.6096_dp,0.0_dp/))
  call vr%assignP(4,(/0.3048_dp,0.3048_dp,0.0_dp/))
  P = (/0.3810_dp,0.1524_dp,0.0_dp/)

  vfx%fc(:,1) = (/100._dp,0._dp,0._dp/)
  vfx%fc(:,2) = vr%vf(2)%fc(:,1)

  vfy%fc(:,1) = vr%vf(3)%fc(:,1)
  vfy%fc(:,2) = (/100._dp,0.3048_dp,0._dp/)

  print*, vr%vind(P)
end program trial
