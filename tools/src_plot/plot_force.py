# Plot using gnuplot
filelist=system('ls Results/r??ForceNonDim.csv')
plot for [filename in filelist] filename u 1:2 w l
pause 1
reread
