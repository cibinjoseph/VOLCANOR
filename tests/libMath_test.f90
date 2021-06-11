module test_libMath
  use naturalfruit
  implicit none
  integer, parameter :: dp = kind(1.d0)

contains

  subroutine test_linspace()
    use libMath, only: linspace
    use naturalfruit
    integer, parameter :: n = 5
    real(dp), dimension(n) :: res, correct

    call testcase_initialize('test_linspace')
    res = linspace(0._dp, 1._dp, 5)
    correct = (/0._dp, 0.25_dp, 0.5_dp, 0.75_dp, 1._dp/)
    call assert_equal(res, correct)
    call testcase_finalize()
  end subroutine test_linspace

  subroutine test_trapz()
    use libMath, only: linspace, trapz
    use naturalfruit
    integer, parameter :: n = 25
    real(dp), dimension(n) :: x, y
    real(dp) :: correct

    call testcase_initialize('test_trapz')
    x = linspace(0._dp, 1._dp, n)
    y = x*x
    ! Compared against numpy.trapz()
    correct = 0.3336226851851851_dp
    call assert_equal(trapz(y, x), correct)
    call testcase_finalize()
  end subroutine test_trapz

  subroutine test_interp1()
    use libMath, only: interp1
    use naturalfruit
    real(dp), dimension(5) :: x
    real(dp), dimension(5) :: y, ya

    x = (/1._dp, 2._dp, 3._dp, 4._dp, 5._dp/)
    y = x*10._dp
    ya = y
    ya(3) = 35._dp

    call assert_equal(interp1(0.5_dp, &
      & (/0._dp, 1._dp/), (/0._dp, 0._dp/), 1), 0.0_dp, delta=1.0D-6)
    call assert_equal(interp1(1.0_dp, x, y, 1), 10._dp, delta=1.0D-6)
    call assert_equal(interp1(2.9_dp, x, y, 1), 29._dp, delta=1.0D-6)
    call assert_equal(interp1(5.0_dp, x, y, 1), 50._dp, delta=1.0D-6)

    call assert_equal(interp1(0.5_dp, &
      & (/0._dp, 1._dp/), (/0._dp, 0._dp/), 2), 0.0_dp, delta=1.0D-6)
    call assert_equal(interp1(1.0_dp, x, y, 2), 10._dp, delta=1.0D-6)
    call assert_equal(interp1(1.2_dp, x, ya, 2), 11.6_dp, delta=1.0D-6)
    call assert_equal(interp1(2.9_dp, x, y, 2), 29._dp, delta=1.0D-6)
    call assert_equal(interp1(2.9_dp, x, ya, 2), 33.275_dp, delta=1.0D-6)
    call assert_equal(interp1(4.8_dp, x, ya, 2), 47.6_dp, delta=1.0D-6)
    call assert_equal(interp1(5.0_dp, x, y, 2), 50._dp, delta=1.0D-6)
  end subroutine test_interp1

  subroutine test_inv2()
    use libMath, only: inv2
    use naturalfruit
    real(dp), dimension(3, 3) :: mat, trueval

    mat(1, :) = [1._dp, -1._dp,  2._dp]
    mat(2, :) = [4._dp,  0._dp,  6._dp]
    mat(3, :) = [0._dp,  1._dp, -1._dp]

    trueval(:, 1) = [ 3.0_dp, -2.0_dp, -2._dp]
    trueval(:, 2) = [-0.5_dp,  0.5_dp,  0.5_dp]
    trueval(:, 3) = [ 3.0_dp, -1.0_dp, -2._dp]

    call assert_equal(trueval, inv2(mat), delta=1.0D-6)
  end subroutine test_inv2

  subroutine test_matmulAX()
    use libMath, only: matmulAX
    use naturalfruit
    real(dp), dimension(3, 3) :: mat
    real(dp), dimension(3) :: x, trueval

    mat(1, :) = [1._dp, -1._dp,  2._dp]
    mat(2, :) = [4._dp,  0._dp,  6._dp]
    mat(3, :) = [0._dp,  1._dp, -1._dp]

    x = [1._dp, 2._dp, 3._dp]

    trueval = [ 5.0_dp, 22.0_dp, -1._dp]

    call assert_equal(trueval, matmulAX(mat, x), delta=1.0D-6)
  end subroutine test_matmulAX
end module test_libMath
