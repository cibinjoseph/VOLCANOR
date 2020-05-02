program panel2_testrunner
  use panel2_test
  integer :: exit_code

  call testsuite_initialize()

  call setup()
  call test_aic()

  call testsuite_summary()
  call testsuite_finalize(exit_code)
  call exit(exit_code)
end program panel2_testrunner