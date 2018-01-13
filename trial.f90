program trial
  use mymathlib
  implicit none 
  real(dp), dimension(2,2) :: A, A_inv

  A(1,:)=(/1._dp,2._dp/)
  A(2,:)=(/3._dp,4._dp/)

  A_inv=inv(A)
  call print_mat(inv(A))


end program trial
