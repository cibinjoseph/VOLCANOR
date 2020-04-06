module vf_test
  use fruit
  implicit none
contains

  subroutine test_vind()
    use vf_classdef

    call assert_equals(4, 4)
    call assert_equals(2, 2)

  end subroutine test_vind

end module vf_test

