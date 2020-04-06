module vf_test
  use fruit
  implicit none
contains

  subroutine test_vind()
    use vf_classdef, only: calclength

    call assert_equals(4, 4)

  end subroutine test_vind

end module vf_test

