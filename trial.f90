program trial
  use vf_classdef
  use library
  implicit none 
  type(vf_class) :: fil
  real(dp) :: l_fil
  real(dp), dimension(100) :: yvec
  real(dp), dimension(3) :: P
  integer :: i
  real(dp), dimension(3) :: ind_vel

  l_fil=8._dp
  fil%fc(:,1)=(/0._dp,0._dp,0._dp/)
  fil%fc(:,2)=(/l_fil,0._dp,0._dp/)
  fil%r_vc0=l_fil*0.03_dp
  fil%r_vc =l_fil*0.03_dp

  yvec=linspace(l_fil*0.01_dp,l_fil,100)

  if (model_switch == 1) then
    open(unit=10,file='dummy_ideal')
  else
    open(unit=10,file='dummy_rankine')
  endif 

  do i=1,size(yvec,1)
    P=(/l_fil*0.5_dp,yvec(i),0._dp/)
    ind_vel=fil%vind(P)

    write(10,'(4F10.7)') yvec(i),ind_vel(1),ind_vel(2),ind_vel(3)
  enddo

  close(10)

  print*,'PROGRAM END'
end program trial
