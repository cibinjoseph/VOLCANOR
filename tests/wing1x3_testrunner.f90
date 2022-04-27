program panel1x3_testrunner
  use panel1x3_test
  integer :: exit_code

  call testsuite_initialize()

  call setup()
  call test_aic()
  call test_force_gamma()

  call testsuite_summary()
  call testsuite_finalize(exit_code)
  call exit(exit_code)
end program panel1x3_testrunner
