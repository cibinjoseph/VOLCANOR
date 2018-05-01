program trial
  use mymathlib
  implicit none 
  integer, parameter :: nx=200
  real(dp), dimension(3,nx) :: A
  real(dp), dimension(3) :: u
  real(dp) :: t, theta,ct,st,onect
  real(dp), dimension(3,3) :: Tmat
  integer :: i

  ! Construct helix
  t=0._dp
  do i=1,nx
    t=t+0.1_dp
    A(1,i)=0.2_dp*cos(t)
    A(2,i)=0.2_dp*sin(t)
    A(3,i)=0.03_dp*t
  enddo

  ! Write out helix
  open(unit=10,file='ahelix1.tec')
  write(10,*) 'Title = Helix'
  write(10,*) 'Variables = "X" "Y" "Z"'
  write(10,*) 'Zone I=200 J=1 K=1'
  do i=1,nx
    write(10,*)A(1,i),A(2,i),A(3,i)
  enddo
  close(10)


  ! Specify axis of rotation
  u=(/1._dp,0._dp,1._dp/)
  u=u/norm2(u)
  print*,u

  theta=60._dp*pi/180._dp
  ct=cos(theta)
  st=sin(theta)
  onect=1-ct

  ! Rotate
  Tmat(:,1)=(/      ct+u(1)*u(1)*onect,  u(3)*st+u(2)*u(1)*onect, -u(2)*st+u(3)*u(1)*onect/)
  Tmat(:,2)=(/-u(3)*st+u(1)*u(2)*onect,       ct+u(2)*u(2)*onect,  u(1)*st+u(3)*u(2)*onect/)
  Tmat(:,3)=(/ u(2)*st+u(1)*u(3)*onect, -u(1)*st+u(2)*u(3)*onect,       ct+u(3)*u(3)*onect/)
   
  do i=1,nx
    A(:,i)=matmul(Tmat,A(:,i))
  enddo

  ! Write out rotated helix
  open(unit=10,file='ahelix2.tec')
  write(10,*) 'Title = Helix'
  write(10,*) 'Variables = "X" "Y" "Z"'
  write(10,*) 'Zone I=200 J=1 K=1'
  do i=1,nx
    write(10,*)A(1,i),A(2,i),A(3,i)
  enddo
  close(10)


end program trial
