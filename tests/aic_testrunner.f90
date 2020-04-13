program aic_testrunner
  use naturalfruit
  use aic_test

  call testsuite_initialize()


  call testsuite_summary()
  call testsuite_finalize()
  
end program aic_testrunner
