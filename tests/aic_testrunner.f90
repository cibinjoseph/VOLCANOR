program aic_testrunner
  use aic_test
  integer :: exit_code

  call testsuite_initialize()

  call setUp()
  call test_1panel()

  call testsuite_summary()
  call testsuite_finalize(exit_code)
  call exit(exit_code)
end program aic_testrunner
