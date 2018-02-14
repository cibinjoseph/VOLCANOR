! This file contains the switches that control certain features in the code
! Vortex model: [1]Ideal Vortex [2]Rankine Vortex
integer, parameter :: model_switch=1

! Span discretization : [1]linear [2]cosine [3]halfsine
integer, parameter :: span_spacing_switch=2

! Wake tip dissipation: [0]Off [1]On
integer, parameter :: tip_diss_switch=0

! Slow start to avoid large starting vortex: [0]Off [1]linear [2]tanh [3]extended tanh
integer, parameter :: slowstart_switch=1
integer, parameter :: slowstart_nt=150

! Wake strain: [0]Off [1]On
integer, parameter :: wakestrain_switch=0

! Plot wake: [0]Off [1]Last timestep only [2]All timestep wakes
integer, parameter :: wakeplot_switch=2

! Wake convection using Predictor-Corrector approach [0]Off1 [1]exPC2 [2]AB4
integer, parameter :: PCwake_switch=2

! Ignore Wake wake interaction [0]Off [1..]timesteps
integer, parameter :: wake_ignore_nt=144
