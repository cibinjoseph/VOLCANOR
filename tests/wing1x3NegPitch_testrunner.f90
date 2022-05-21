program wing1x3NegPitch_testrunner
  use wing1x3NegPitch_test
  integer :: exit_code

  call testsuite_initialize()

  call setup()
  call test_aic()
  call test_force()

  call testsuite_summary()
  call testsuite_finalize(exit_code)
  call exit(exit_code)
end program wing1x3NegPitch_testrunner
