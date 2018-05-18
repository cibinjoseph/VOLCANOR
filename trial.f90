program trial
  use library
  use postproc

  ! Variable Initialization
  include "init_file.f90"

  ! Read config.in file
  open(unit=11,file='config.in_dummy')
  call skiplines(11,2)
  read(11,*) nt,dt,nr
  call skiplines(11,4)
  read(11,*) density
  call skiplines(11,4)
  read(11,*) span_spacing_switch
  call skiplines(11,4)
  read(11,*) tip_diss_switch, wakestrain_switch
  call skiplines(11,4)
  read(11,*) slowstart_switch, slowstart_nt
  call skiplines(11,4)
  read(11,*) wakeplot_switch 
  call skiplines(11,4)
  read(11,*) FDscheme_switch
  call skiplines(11,4)
  read(11,*) wake_ignore_nt
  call skiplines(11,4)
  read(11,*) init_wake_vel_nt
  close(11)

end program trial
