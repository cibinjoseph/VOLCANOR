program testrunner_libMath
  use test_libMath
  integer :: exit_code

  call testsuite_initialize()

  call test_linspace()
  call test_trapz()
  call test_interp1d()
  call test_inv2()
  call test_matmulAX()

  call testsuite_summary()
  call testsuite_finalize(exit_code)
  call exit(exit_code)
end program testrunner_libMath
