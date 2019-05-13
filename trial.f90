program trial
  use vf_classdef
  implicit none 
  integer, parameter :: n = 100
  type(vf_class) :: vf
  real(dp), dimension(3,n) :: r,v
  real(dp) :: a
  integer :: i

  r=0._dp
  r(1,:)=linspace(0._dp,1._dp,n)

  vf%fc(:,1)=(/0._dp, 20._dp,0._dp/)
  vf%fc(:,2)=(/0._dp,-20._dp,0._dp/)
  vf%rVc0=0.14_dp
  vf%rVc=vf%rVc0

  do i=1,n
    v(:,i)=vf%vind(r(:,i))
  enddo

  do i=1,n
    print*,r(1,i),v(3,i)
  enddo

end program trial
