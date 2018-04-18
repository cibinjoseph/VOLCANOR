! This file contains the switches that control certain features in the code

! Span discretization : [1]linear [2]cosine [3]halfsine
integer, parameter :: span_spacing_switch=1

! Wake tip dissipation: [0]Off [1]On
integer, parameter :: tip_diss_switch=0

! Slow start to avoid large starting vortex: [0]Off [1]linear [2]tanh [3]extended tanh
integer, parameter :: slowstart_switch=0
integer, parameter :: slowstart_nt=50

! Wake strain: [0]Off [1]On
integer, parameter :: wakestrain_switch=0

! Plot wake: [0]Off [1]Last timestep only [2]All timestep wakes
integer, parameter :: wakeplot_switch=2

! Wake convection finite-difference scheme [0]FT1 [1]PC2 [2]AB2 [3]AB4
integer, parameter :: FDscheme_switch=0

! Ignore Wake wake interaction [0]Off [1..]timesteps
integer, parameter :: wake_ignore_nt=0

! Initial Wake velocity [0]Off [1..]timesteps
integer, parameter :: init_wake_vel_nt=0
