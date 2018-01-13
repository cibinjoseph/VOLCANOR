! This file contains the switches that control certain features in the code
! Vortex model: [1]Ideal Vortex [2]Rankine Vortex
integer, parameter :: model_switch=2

! Wake tip dissipation: [0]Off [1]On
integer, parameter :: tip_diss_switch=1

! Slow start to avoid large starting vortex [0]Off [1]linear [2]tanh [3]extended tanh
integer, parameter :: slowstart_switch=0
integer, parameter :: slowstart_nt=20
