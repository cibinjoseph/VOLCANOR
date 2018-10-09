module library
  use rotor_classdef
  implicit none

contains

  !--------------------------------------------------------!
  !                Induced Velocity Functions              !
  !--------------------------------------------------------!

  ! Induced velocity by rotor (wing n wake) on Nwake corner points
  function vind_onNwake_byrotor(rotor,Nwake,opt_char) result(vind_array)
    type(rotor_class), intent(inout) :: rotor
    type(Nwake_class), intent(in), dimension(:,:) :: Nwake
    character(len=1), optional :: opt_char
    real(dp), dimension(3,size(Nwake,1),size(Nwake,2)+1) :: vind_array
    integer :: i,j,rows

    rows=size(Nwake,1)

    if (.not. present(opt_char)) then
      ! Induced velocity due to all blades and wake
      !$omp parallel do collapse(2) 
      do j=1,rotor%ns
        do i=1,rows
          vind_array(:,i,j)=rotor%vind_bywing(Nwake(i,j)%vr%vf(2)%fc(:,1))  &
            +               rotor%vind_bywake(Nwake(i,j)%vr%vf(2)%fc(:,1))
        enddo
      enddo
      !$omp end parallel do

      !$omp parallel do 
      do i=1,rows
        vind_array(:,i,rotor%ns+1)=rotor%vind_bywing(Nwake(i,rotor%ns)%vr%vf(3)%fc(:,1))  &
          +                        rotor%vind_bywake(Nwake(i,rotor%ns)%vr%vf(3)%fc(:,1))
      enddo
      !$omp end parallel do

    elseif ((opt_char .eq. 'P') .or. (opt_char .eq. 'p')) then

      ! Induced velocity due to all blades and Pwake
      !$omp parallel do collapse(2) 
      do j=1,rotor%ns
        do i=1,rows
          vind_array(:,i,j)=rotor%vind_bywing(Nwake(i,j)%vr%vf(2)%fc(:,1))  &
            +               rotor%vind_bywake(Nwake(i,j)%vr%vf(2)%fc(:,1),'P')
        enddo
      enddo
      !$omp end parallel do

      !$omp parallel do 
      do i=1,rows
        vind_array(:,i,rotor%ns+1)=rotor%vind_bywing(Nwake(i,rotor%ns)%vr%vf(3)%fc(:,1))  &
          +                        rotor%vind_bywake(Nwake(i,rotor%ns)%vr%vf(3)%fc(:,1),'P')
      enddo
      !$omp end parallel do

    else
      error stop 'ERROR: Wrong character flag for vind_onNwake_byrotor()'
    endif

  end function vind_onNwake_byrotor

  ! Induced velocity by rotor (wing n wake) on Fwake corner points
  function vind_onFwake_byrotor(rotor,Fwake,opt_char) result(vind_array)
    type(rotor_class), intent(inout) :: rotor
    type(Fwake_class), intent(in), dimension(:) :: Fwake
    character(len=1), optional :: opt_char
    real(dp), dimension(3,size(Fwake,1)) :: vind_array
    integer :: i,rows

    rows=size(Fwake,1)

    if (.not. present(opt_char)) then
      ! Induced velocity due to all blades and wake
      !$omp parallel do
      do i=1,rows
        vind_array(:,i)=rotor%vind_bywing(Fwake(i)%vf%fc(:,1))  &
          +             rotor%vind_bywake(Fwake(i)%vf%fc(:,1))
      enddo
      !$omp end parallel do

    elseif ((opt_char .eq. 'P') .or. (opt_char .eq. 'p')) then

      ! Induced velocity due to all blades and Pwake
      !$omp parallel do
      do i=1,rows
        vind_array(:,i)=rotor%vind_bywing(Fwake(i)%vf%fc(:,1))  &
          +             rotor%vind_bywake(Fwake(i)%vf%fc(:,1),'P')
      enddo
      !$omp end parallel do

    else
      error stop 'ERROR: Wrong character flag for vind_onFwake_byrotor()'
    endif

  end function vind_onFwake_byrotor

  ! Calculates 2nd order accurate induced velocity on near wake
  function vel_order2_Nwake(v_wake_n,v_wake_np1)   ! np1 => n+1 
    real(dp), intent(in), dimension(:,:,:) :: v_wake_n, v_wake_np1
    real(dp), dimension(3,size(v_wake_n,2),size(v_wake_n,3)) :: vel_order2_Nwake
    integer :: i,j
    do j=1,size(v_wake_n,3)
      vel_order2_Nwake(:,1,j)=(v_wake_np1(:,1,j)+v_wake_n(:,1,j))*0.5_dp
      do i=2,size(v_wake_n,2)-1
        vel_order2_Nwake(:,i,j)=(v_wake_np1(:,i,j)+v_wake_np1(:,i-1,j)+v_wake_n(:,i+1,j)+v_wake_n(:,i,j))*0.25_dp
      enddo
      vel_order2_Nwake(:,size(v_wake_n,2),j)=(v_wake_np1(:,size(v_wake_n,2),j)+v_wake_n(:,size(v_wake_n,2),j))*0.5_dp
    enddo
  end function vel_order2_Nwake

  ! Calculates 2nd order accurate induced velocity on far wake
  function vel_order2_Fwake(v_wake_n,v_wake_np1)   ! np1 => n+1 
    real(dp), intent(in), dimension(:,:) :: v_wake_n, v_wake_np1
    real(dp), dimension(3,size(v_wake_n,2)) :: vel_order2_Fwake
    integer :: i
    vel_order2_Fwake(:,1)=(v_wake_np1(:,1)+v_wake_n(:,1))*0.5_dp
    !$omp parallel do
    do i=2,size(v_wake_n,2)-1
      vel_order2_Fwake(:,i)=(v_wake_np1(:,i)+v_wake_np1(:,i-1)+v_wake_n(:,i+1)+v_wake_n(:,i))*0.25_dp
    enddo
    !$omp end parallel do
    vel_order2_Fwake(:,size(v_wake_n,2))=(v_wake_np1(:,size(v_wake_n,2))+v_wake_n(:,size(v_wake_n,2)))*0.5_dp
  end function vel_order2_Fwake

  !--------------------------------------------------------!
  !               Force Computation Functions              !
  !--------------------------------------------------------!

  !subroutine calc_wingalpha(wing_array)
  !  type(wingpanel_class), intent(in), dimension(:,:) :: wing_array
  !  integer :: i,j
  !  do j=1,size(wing_array,2)
  !    do i=1,size(wing_array,1)
  !      call wing_array(i,j)%calc_alpha()
  !    enddo
  !  enddo
  !end subroutine calc_wingalpha

  !function calcgam(wg)
  !  type(wingpanel_class), intent(inout), dimension(:,:) :: wg  !short form for wing_array
  !  real(dp), dimension(size(wg,2)) :: calcgam
  !  integer :: j,rows,cols

  !  rows=size(wg,1)
  !  cols=size(wg,2)

  !  ! Check if this is correct way of calculating sectional circulation
  !  do j=2,cols
  !    calcgam(j)=wg(rows,j)%vr%gam
  !  enddo

  !end function calcgam

  !function calclift(wg,gamvec_prev,dt)
  !  type(wingpanel_class), intent(inout), dimension(:,:) :: wg  !short form for wing_array
  !  real(dp), intent(in), dimension(:) :: gamvec_prev
  !  real(dp), intent(in) :: dt
  !  real(dp) :: calclift
  !  real(dp), dimension(size(wg,1),size(wg,2)) :: gam_prev
  !  real(dp), dimension(3) :: tau_c, tau_s
  !  integer :: i,j,rows,cols
  !  ! Inherent assumption that panels have subdivisions along chord and not inclined to it
  !  ! while calculating tangent vector
  !  ! LE and left sides used for calculating tangent vectors

  !  rows=size(wg,1)
  !  cols=size(wg,2)

  !  gam_prev=reshape(gamvec_prev,(/rows,cols/))
  !  do j=2,cols
  !    do i=2,rows
  !      tau_c=wg(i,j)%pc(:,2)-wg(i,j)%pc(:,1)
  !      tau_s=wg(i,j)%pc(:,4)-wg(i,j)%pc(:,1)
  !      wg(i,j)%delP=dot_product(wg(i,j)%velCP,tau_c)*(wg(i,j)%vr%gam-wg(i-1,j)%vr%gam)/dot_product(tau_c,tau_c) &
  !        +          dot_product(wg(i,j)%velCP,tau_s)*(wg(i,j)%vr%gam-wg(i,j-1)%vr%gam)/dot_product(tau_s,tau_s) &
  !        +          (wg(i,j)%vr%gam-gam_prev(i,j))/dt
  !    enddo
  !  enddo

  !  do j=2,cols
  !    tau_c=wg(1,j)%pc(:,2)-wg(1,j)%pc(:,1)
  !    tau_s=wg(1,j)%pc(:,4)-wg(1,j)%pc(:,1)
  !    wg(1,j)%delP=dot_product(wg(1,j)%velCP,tau_c)*(wg(1,j)%vr%gam)/dot_product(tau_c,tau_c) &
  !      +          dot_product(wg(1,j)%velCP,tau_s)*(wg(1,j)%vr%gam-wg(1,j-1)%vr%gam)/dot_product(tau_s,tau_s) &
  !      +          (wg(1,j)%vr%gam-gam_prev(1,j))/dt
  !  enddo

  !  tau_c=wg(1,1)%pc(:,2)-wg(1,1)%pc(:,1)
  !  tau_s=wg(1,1)%pc(:,4)-wg(1,1)%pc(:,1)
  !  wg(1,1)%delP=dot_product(wg(1,1)%velCP,tau_c)*(wg(1,1)%vr%gam)/dot_product(tau_c,tau_c) &
  !    +          dot_product(wg(1,1)%velCP,tau_s)*(wg(1,1)%vr%gam)/dot_product(tau_s,tau_s) &
  !    +          (wg(1,1)%vr%gam-gam_prev(1,1))/dt

  !  do i=2,rows
  !    tau_c=wg(i,1)%pc(:,2)-wg(i,1)%pc(:,1)
  !    tau_s=wg(i,1)%pc(:,4)-wg(i,1)%pc(:,1)
  !    wg(i,1)%delP=dot_product(wg(i,1)%velCP,tau_c)*(wg(i,1)%vr%gam-wg(i-1,1)%vr%gam)/dot_product(tau_c,tau_c) &
  !      +          dot_product(wg(i,1)%velCP,tau_s)*(wg(i,1)%vr%gam)/dot_product(tau_s,tau_s) &
  !      +          (wg(i,1)%vr%gam-gam_prev(i,1))/dt
  !  enddo
  !  wg%delP=density*wg%delP

  !  do j=1,cols
  !    do i=1,rows
  !      wg(i,j)%dLift=-(wg(i,j)%delP*wg(i,j)%panel_area)*cos(wg(i,j)%alpha)
  !    enddo
  !  enddo

  !  calclift=0._dp
  !  do j=1,cols
  !    do i=1,rows
  !      calclift=calclift+wg(i,j)%dlift
  !    enddo
  !  enddo
  !end function calclift

  !function calcdrag(wg,gamvec_prev,dt)
  !  type(wingpanel_class), intent(inout), dimension(:,:) :: wg !short form for wing_array
  !  real(dp), intent(in), dimension(:) :: gamvec_prev
  !  real(dp) :: calcdrag
  !  real(dp), intent(in) :: dt
  !  real(dp) :: vel_drag
  !  real(dp) :: drag1, drag2
  !  real(dp), dimension(size(wg,1),size(wg,2)) :: gam_prev
  !  integer :: i,j,rows,cols
  !  ! Inherent assumption that panels have subdivisions along chord and not inclined to it
  !  ! while calculating tangent vector
  !  ! LE and left sides used for calculating tangent vectors

  !  ! *** !! PREDICTS DRAG1 INCORRECTLY !! ***
  !  ! *** !! PREDICTS DRAG1 INCORRECTLY !! ***
  !  ! *** !! PREDICTS DRAG1 INCORRECTLY !! ***

  !  rows=size(wg,1)
  !  cols=size(wg,2)

  !  gam_prev=reshape(gamvec_prev,(/rows,cols/))
  !  do j=1,cols
  !    do i=2,rows
  !      !vel_drag=dot_product((vind_panelgeo(wake_array,wg(i,j)%cp))+vind_chordvortex(wg,wg(i,j)%cp),&
  !      !  (/0._dp,0._dp,1._dp/))
  !      vel_drag=dot_product((wg(i,j)%velCP-wg(i,j)%velCPm)+vind_chordvortex(wg,wg(i,j)%CP),&
  !        matmul(wg(i,j)%orthproj(),wg(i,j)%ncap))
  !      drag2=(wg(i,j)%vr%gam-gam_prev(i,j))*wg(i,j)%panel_area*sin(wg(i,j)%alpha)/dt
  !      drag1=-vel_drag*(wg(i,j)%vr%gam-wg(i-1,j)%vr%gam)*norm2(wg(i,j)%pc(:,4)-wg(i,j)%pc(:,1))
  !      wg(i,j)%dDrag=drag1-drag2
  !    enddo
  !  enddo

  !  ! i=1
  !  do j=2,cols
  !    !vel_drag=dot_product((vind_panelgeo(wake_array,wg(1,j)%cp))+vind_chordvortex(wg,wg(1,j)%cp),&
  !    !  (/0._dp,0._dp,1._dp/))
  !    vel_drag=dot_product((wg(1,j)%velCP-wg(1,j)%velCPm)+vind_chordvortex(wg,wg(1,j)%CP),&
  !      matmul(wg(1,j)%orthproj(),wg(1,j)%ncap))
  !    drag2=(wg(1,j)%vr%gam-gam_prev(1,j))*wg(1,j)%panel_area*sin(wg(1,j)%alpha)/dt
  !    drag1=-vel_drag*(wg(1,j)%vr%gam)*norm2(wg(1,j)%pc(:,4)-wg(1,j)%pc(:,1))
  !    wg(1,j)%dDrag=drag1-drag2
  !  enddo

  !  wg%dDrag=density*wg%dDrag

  !  calcdrag=0._dp
  !  do j=1,cols
  !    do i=1,rows
  !      calcdrag=calcdrag+wg(i,j)%dDrag
  !    enddo
  !  enddo
  !end function calcdrag

end module library
