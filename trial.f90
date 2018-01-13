program trial
  use mymathlib
  implicit none 
  real(dp), dimension(3,3) :: A

  A(1,:)=(/1._dp,2._dp,1._dp/)
  A(2,:)=(/3._dp,4._dp,6._dp/)
  A(3,:)=(/5._dp,1._dp,3._dp/)

  call print_mat(inv(A))


end program trial
