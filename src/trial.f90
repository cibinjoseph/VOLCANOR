program trial
  use classdef
  implicit none 
  type(vr_class) :: vr
  type(vf_class) :: vfx, vfy, vfwake
  real(dp), dimension(3) :: velCP, velRel, velInd, P
  real(dp), dimension(3,4) :: C
  real(dp) :: gam, theta, alpha, sweep_rad, chord, vinf

  sweep_rad = 35._dp*pi/180._dp
  chord = 0.3048_dp
  vinf = 10._dp
  theta = 5._dp*pi/180._dp

  C(:,1) = (/0.0000_dp,0.0000_dp,0.0_dp/)
  C(:,2) = (/chord,0.0000_dp,0.0_dp/)
  C(:,3) = (/chord*tan(sweep_rad)+chord,chord,0.0_dp/)
  C(:,4) = (/chord*tan(sweep_rad),chord,0.0_dp/)

  ! Swept
  call vr%assignP(1,(/C(1,1)+chord*0.25_dp,C(2,1),C(3,1)/))
  call vr%assignP(2,C(:,2))
  call vr%assignP(3,C(:,3))
  call vr%assignP(4,(/C(1,4)+chord*0.25_dp,C(2,4),C(3,4)/))
  P = (/(C(1,1)+C(1,4))*0.5_dp+chord*0.75_dp,chord*0.5_dp,0.0_dp/)

  vfx%fc(:,1) = vr%vf(2)%fc(:,1)
  vfx%fc(:,2) = (/100._dp,0._dp,0._dp/)

  vfy%fc(:,1) = (/100._dp,chord,0._dp/)
  vfy%fc(:,2) = vr%vf(3)%fc(:,1)

  vfwake%fc(:,1) = vfy%fc(:,2)
  vfwake%fc(:,2) = vfx%fc(:,1)

  velCP = vr%vind(P)+vfx%vind(P)+vfy%vind(P)+vfwake%vind(P)
  gam = -vinf*sin(theta)/velCP(3)
  print*, 'gam = ', gam

  velInd = gam*(vfx%vind(P)+vfy%vind(P)+vr%vf(1)%vind(P)+vr%vf(3)%vind(P))
  velRel = velInd+vinf*(/cos(theta),0._dp,sin(theta)/)
  alpha = atan(velRel(3)/velRel(1))

  print*, 'velRel = ', velRel
  print*, 'alpha = ', alpha*180._dp/pi
  print*, 'CL_alpha = ', alpha*2._dp*pi
  print*, 'CL_gamma = ', 2._dp*abs(gam)/10._dp/0.3048_dp*cos(theta)

end program trial
