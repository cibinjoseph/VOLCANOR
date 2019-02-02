program trial
  use mymathlib
  real(dp), dimension(6) :: xData, yData 
  real(dp), dimension(3) :: out
  xData=(/0._dp, 0.5_dp, 1._dp, 1.5_dp, 2._dp, 2.5_dp/)
  yData=(/0.0674_dp, -0.9156_dp, 1.6253_dp, 3.0377_dp, 3.3535_dp, 7.9409_dp/)

  out=lsq2((/0._dp,1._dp,2._dp/),xData,yData)
  print*,out
end program trial
