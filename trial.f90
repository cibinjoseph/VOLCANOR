program trial
  use mymathlib
  implicit none 
  integer, parameter :: nx=200
  real(dp), dimension(3,nx) :: A
  real(dp) :: t
  real(dp) :: pivotLE
  integer :: i

  t=0._dp
  do i=1,nx
    t=t+0.1_dp
    A(1,i)=0.2_dp*cos(t)
    A(2,i)=0.2_dp*sin(t)
    A(3,i)=0.03_dp*t
  enddo

  open(unit=10,file='ahelix.tec')
  write(10,*) 'Title = Helix'
  write(10,*) 'Variables = "X" "Y" "Z"'
  write(10,*) 'Zone I=200 J=1 K=1'
  do i=1,nx
    write(10,*)A(1,i),A(2,i),A(3,i)
  enddo
  close(10)
end program trial
