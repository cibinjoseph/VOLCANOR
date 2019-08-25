program trial
  use vr_classdef
  implicit none 
  type(vr_class) :: vr
  type(vf_class) :: vfx, vfy, vfwake
  real(dp), dimension(3) :: velCP, velRel, P
  real(dp) :: gam, alpha

  ! Rectangular
  !call vr%assignP(1,(/0.0000_dp,0.0000_dp,0.0_dp/))
  !call vr%assignP(2,(/0.3048_dp,0.0000_dp,0.0_dp/))
  !call vr%assignP(3,(/0.3048_dp,0.3048_dp,0.0_dp/))
  !call vr%assignP(4,(/0.3048_dp,0.0000_dp,0.0_dp/))
  !P = (/0.2286_dp,0.1524_dp,0.0_dp/)

  ! Swept
  call vr%assignP(1,(/0.0000_dp,0.0000_dp,0.0_dp/))
  call vr%assignP(2,(/0.3048_dp,0.0000_dp,0.0_dp/))
  call vr%assignP(3,(/0.6096_dp,0.3048_dp,0.0_dp/))
  call vr%assignP(4,(/0.3048_dp,0.3048_dp,0.0_dp/))
  P = (/0.3810_dp,0.1524_dp,0.0_dp/)

  vfx%fc(:,1) = vr%vf(2)%fc(:,1)
  vfx%fc(:,2) = (/100._dp,0._dp,0._dp/)

  vfy%fc(:,1) = (/100._dp,0.3048_dp,0._dp/)
  vfy%fc(:,2) = vr%vf(3)%fc(:,1)

  vfwake%fc(:,1) = vfy%fc(:,2)
  vfwake%fc(:,2) = vfx%fc(:,1)

  velCP = vr%vind(P)+vfx%vind(P)+vfy%vind(P)+vfwake%vind(P)
  gam = -0.8716/velCP(3)
  print*, 'gam = ', gam

  velRel = gam*(vfx%vind(P)+vfy%vind(P)+vr%vf(1)%vind(P)+vr%vf(3)%vind(P))+(/9.9619_dp,0._dp,0.8716_dp/)
  alpha = atan(velRel(3)/velRel(1))

  print*, 'velRel = ', velRel
  print*, 'alpha = ', alpha*180._dp/pi
  print*, 'CL = ', alpha*2._dp*pi

end program trial
