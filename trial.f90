program trial
  use mymathlib
  implicit none 
  real(dp), dimension(3,2,2) :: A
  real(dp), dimension(2,2) :: a1,a2,a3

  a1(1,:)=(/1._dp,1._dp/)
  a1(2,:)=(/1._dp,1._dp/)
  
  a2(1,:)=(/2._dp,2._dp/)
  a2(2,:)=(/2._dp,2._dp/)
  
  a3(1,:)=(/3._dp,3._dp/)
  a3(2,:)=(/3._dp,3._dp/)

  A(1,:,:)=a1*0._dp
  A(2,:,:)=a2*0._dp
  A(3,:,:)=a3

  print*,sum(A,1)
  print*
  print*,sum(A,2)
  print*
  print*,sum(A,3)


end program trial
