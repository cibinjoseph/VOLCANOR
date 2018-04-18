module postproc
  use wingpanel_classdef
  use wakepanel_classdef

contains

  subroutine mesh2file(wing_array,wake_array,filename)
    type(wingpanel_class), intent(in), dimension(:,:) :: wing_array
    type(wakepanel_class), intent(in), dimension(:,:) :: wake_array
    character(len=*), intent(in) :: filename
    character(len=5) :: nx_char, ny_char
    real(dp), dimension(3,size(wing_array,1)+1,size(wing_array,2)+1) :: wing_mesh  
    real(dp), dimension(3,size(wake_array,1)+1,size(wake_array,2)+1) :: wake_mesh  
    integer :: i,j,nx,ny

    nx=size(wing_array,1)
    ny=size(wing_array,2)
    write(nx_char,'(I5)') nx+1
    write(ny_char,'(I5)') ny+1

    open(unit=10,file=filename,position='append')
    do j=1,ny
      do i=1,nx
        wing_mesh(:,i,j)=wing_array(i,j)%pc(:,1)
      enddo
    enddo
    do i=1,nx
      wing_mesh(:,i,ny+1)=wing_array(i,ny)%pc(:,4)
    enddo
    do j=1,ny
      wing_mesh(:,nx+1,j)=wing_array(nx,j)%pc(:,2)
    enddo
    wing_mesh(:,nx+1,ny+1)=wing_array(nx,ny)%pc(:,3)

    write(10,*) 'Title = "Panel array"'
    write(10,*) 'VARIABLES = "X" "Y" "Z" "GAM" "TAG"'
    write(10,*) 'Zone I='//trim(nx_char)//' J='//trim(ny_char)//' K=1  T="Wing"'
    write(10,*) 'DATAPACKING=BLOCK'
    write(10,*) 'VARLOCATION=([4]=CELLCENTERED,[5]=CELLCENTERED)'
    write(10,*) ((wing_mesh(1,i,j),i=1,nx+1),j=1,ny+1)
    write(10,*) ((wing_mesh(2,i,j),i=1,nx+1),j=1,ny+1)
    write(10,*) ((wing_mesh(3,i,j),i=1,nx+1),j=1,ny+1)
    write(10,*) ((-1._dp*wing_array(i,j)%vr%gam,i=1,nx),j=1,ny)
    write(10,*) ((wing_array(i,j)%tag,i=1,nx),j=1,ny)

    nx=size(wake_array,1)
    ny=size(wake_array,2)
    write(nx_char,'(I5)') nx+1
    write(ny_char,'(I5)') ny+1

    !Check if necessary - $omp parallel do collapse(2)
    do j=1,ny
      do i=1,nx
        wake_mesh(:,i,j)=wake_array(i,j)%vr%vf(1)%fc(:,1)
      enddo
    enddo
    !Check if necessary -$omp end parallel do
    do i=1,nx
      wake_mesh(:,i,ny+1)=wake_array(i,ny)%vr%vf(4)%fc(:,1)
    enddo
    do j=1,ny
      wake_mesh(:,nx+1,j)=wake_array(nx,j)%vr%vf(2)%fc(:,1)
    enddo
    wake_mesh(:,nx+1,ny+1)=wake_array(nx,ny)%vr%vf(3)%fc(:,1)

    write(10,*) 'Zone I='//trim(nx_char)//' J='//trim(ny_char)//' K=1  T="Wake"'
    write(10,*) 'DATAPACKING=BLOCK'
    write(10,*) 'VARLOCATION=([4]=CELLCENTERED,[5]=CELLCENTERED)'
    write(10,*) ((wake_mesh(1,i,j),i=1,nx+1),j=1,ny+1)
    write(10,*) ((wake_mesh(2,i,j),i=1,nx+1),j=1,ny+1)
    write(10,*) ((wake_mesh(3,i,j),i=1,nx+1),j=1,ny+1)
    write(10,*) ((-1._dp*wake_array(i,j)%vr%gam,i=1,nx),j=1,ny)
    write(10,*) ((wake_array(i,j)%tag,i=1,nx),j=1,ny)

    close(10)
  end subroutine mesh2file

  subroutine wingverify(wing_array)
    ! For verifying orientation of wing panels, bound vortex rings and CPs
    type(wingpanel_class), intent(in), dimension(:,:) :: wing_array
    character(len=5) :: nx_char, ny_char
    real(dp), dimension(3,size(wing_array,1)+1,size(wing_array,2)+1) :: wing_mesh  
    integer :: i,j,nx,ny

    nx=size(wing_array,1)
    ny=size(wing_array,2)
    write(nx_char,'(I5)') nx+1
    write(ny_char,'(I5)') ny+1

    do j=1,ny
      do i=1,nx
        wing_mesh(:,i,j)=wing_array(i,j)%pc(:,1)
      enddo
    enddo
    do i=1,nx
      wing_mesh(:,i,ny+1)=wing_array(i,ny)%pc(:,4)
    enddo
    do j=1,ny
      wing_mesh(:,nx+1,j)=wing_array(nx,j)%pc(:,2)
    enddo
    wing_mesh(:,nx+1,ny+1)=wing_array(nx,ny)%pc(:,3)

    open(unit=10,file='Results/wingPC.tec',position='append')
    write(10,*) 'Title = "Panel Vertices"'
    write(10,*) 'VARIABLES = "X" "Y" "Z"'
    write(10,*) 'Zone I='//trim(nx_char)//' J='//trim(ny_char)//' K=1  T="Panel Vertices"'
    write(10,*) 'DATAPACKING=BLOCK'
    write(10,*) ((wing_mesh(1,i,j),i=1,nx+1),j=1,ny+1)
    write(10,*) ((wing_mesh(2,i,j),i=1,nx+1),j=1,ny+1)
    write(10,*) ((wing_mesh(3,i,j),i=1,nx+1),j=1,ny+1)
    close(10)

    write(nx_char,'(I5)') nx
    write(ny_char,'(I5)') ny
    do j=1,ny
      do i=1,nx
        wing_mesh(:,i,j)=wing_array(i,j)%cp
      enddo
    enddo
    open(unit=11,file='Results/wingCP.tec',position='append')
    write(11,*) 'Title = "Coll. points"'
    write(11,*) 'VARIABLES = "X" "Y" "Z"'
    write(11,*) 'Zone I='//trim(nx_char)//' J='//trim(ny_char)//' K=1  T="Coll. points"'
    write(11,*) 'DATAPACKING=BLOCK'
    write(11,*) ((wing_mesh(1,i,j),i=1,nx),j=1,ny)
    write(11,*) ((wing_mesh(2,i,j),i=1,nx),j=1,ny)
    write(11,*) ((wing_mesh(3,i,j),i=1,nx),j=1,ny)
    close(11)

    write(nx_char,'(I5)') nx+1
    write(ny_char,'(I5)') ny+1
    do j=1,ny
      do i=1,nx
        wing_mesh(:,i,j)=wing_array(i,j)%vr%vf(1)%fc(:,1)
      enddo
    enddo
    do i=1,nx
      wing_mesh(:,i,ny+1)=wing_array(i,ny)%vr%vf(4)%fc(:,1)
    enddo
    do j=1,ny
      wing_mesh(:,nx+1,j)=wing_array(nx,j)%vr%vf(2)%fc(:,1)
    enddo
    wing_mesh(:,nx+1,ny+1)=wing_array(nx,ny)%vr%vf(3)%fc(:,1)

    open(unit=12,file='Results/wingVR.tec',position='append')
    write(12,*) 'Title = "Vortex Rings"'
    write(12,*) 'VARIABLES = "X" "Y" "Z"'
    write(12,*) 'Zone I='//trim(nx_char)//' J='//trim(ny_char)//' K=1  T="Vortex Rings"'
    write(12,*) 'DATAPACKING=BLOCK'
    write(12,*) ((wing_mesh(1,i,j),i=1,nx+1),j=1,ny+1)
    write(12,*) ((wing_mesh(2,i,j),i=1,nx+1),j=1,ny+1)
    write(12,*) ((wing_mesh(3,i,j),i=1,nx+1),j=1,ny+1)
    close(12)
  end subroutine wingverify

  subroutine tip2file(wing_array,wake_array,filename)
    type(wingpanel_class), intent(in), dimension(:,:) :: wing_array
    type(wakepanel_class), intent(in), dimension(:,:) :: wake_array
    character(len=*), intent(in) :: filename
    character(len=5) :: nx_char, ny_char
    real(dp), dimension(3,size(wing_array,1)+1,size(wing_array,2)+1) :: wing_mesh  
    real(dp), dimension(3,size(wake_array,1)+1) :: wake_tip  
    integer :: i,j,nx,ny

    nx=size(wing_array,1)
    ny=size(wing_array,2)
    write(nx_char,'(I5)') nx+1
    write(ny_char,'(I5)') ny+1

    open(unit=10,file=filename,position='append')
    do j=1,ny
      do i=1,nx
        wing_mesh(:,i,j)=wing_array(i,j)%pc(:,1)
      enddo
    enddo
    do i=1,nx
      wing_mesh(:,i,ny+1)=wing_array(i,ny)%pc(:,4)
    enddo
    do j=1,ny
      wing_mesh(:,nx+1,j)=wing_array(nx,j)%pc(:,2)
    enddo
    wing_mesh(:,nx+1,ny+1)=wing_array(nx,ny)%pc(:,3)

    write(10,*) 'Title = "Panel array"'
    write(10,*) 'VARIABLES = "X" "Y" "Z"'! "gam" "tag"'
    write(10,*) 'Zone I='//trim(nx_char)//' J='//trim(ny_char)//' K=1  T="Wing"'
    write(10,*) 'DATAPACKING=BLOCK'
    write(10,*) ((wing_mesh(1,i,j),i=1,nx+1),j=1,ny+1)
    write(10,*) ((wing_mesh(2,i,j),i=1,nx+1),j=1,ny+1)
    write(10,*) ((wing_mesh(3,i,j),i=1,nx+1),j=1,ny+1)

    ! Wake root
    nx=size(wake_array,1)
    ny=size(wake_array,2)
    write(nx_char,'(I5)') nx+1

    do i=1,nx
      wake_tip(:,i)=wake_array(i,1)%vr%vf(1)%fc(:,1)
    enddo
    wake_tip(:,nx+1)=wake_array(nx,1)%vr%vf(2)%fc(:,1)

    write(10,*) 'Zone I='//trim(nx_char)//' J=1   K=1  T="wake_root"'
    write(10,*) 'DATAPACKING=BLOCK'
    write(10,*) (wake_tip(1,i),i=1,nx+1)
    write(10,*) (wake_tip(2,i),i=1,nx+1)
    write(10,*) (wake_tip(3,i),i=1,nx+1)

    ! Wake tip
    do i=1,nx
      wake_tip(:,i)=wake_array(i,ny)%vr%vf(4)%fc(:,1)
    enddo
    wake_tip(:,nx+1)=wake_array(nx,ny)%vr%vf(3)%fc(:,1)

    write(10,*) 'Zone I='//trim(nx_char)//' J=1   K=1  T="wake_tip"'
    write(10,*) 'DATAPACKING=BLOCK'
    write(10,*) (wake_tip(1,i),i=1,nx+1)
    write(10,*) (wake_tip(2,i),i=1,nx+1)
    write(10,*) (wake_tip(3,i),i=1,nx+1)
    close(10)
  end subroutine tip2file

  subroutine gam2file(yvec,gam_sectional,filename)
    real(dp), intent(in), dimension(:) :: yvec
    real(dp), intent(in), dimension(:) :: gam_sectional
    character(len=*), intent(in) :: filename
    real(dp), dimension(size(yvec)-1) :: yvec_file
    integer :: i,ny

    ny=size(gam_sectional)

    open(unit=10,file=filename,position='append')
    do i=1,ny
      yvec_file(i)=(yvec(i)+yvec(i+1))*0.5_dp
    enddo

    write(10,*) '# Gamma'
    do i=1,ny
      write(10,*) yvec_file(i), -1._dp*gam_sectional(i)
    enddo
    close(10)

  end subroutine gam2file

  subroutine lift2file(liftvec,filename,extra_params)
    real(dp), intent(in), dimension(:) :: liftvec
    character(*), intent(in) :: filename
    real(dp), intent(in), dimension(:) :: extra_params ![1]dt [2]omega(rad/s) [3]span(m) [4]speed(m/s)
    integer :: i

    open(unit=10,file=filename)
    write(10,*) '# Lift'
    do i=1,size(liftvec,1)
      write(10,*) extra_params(1)*i,liftvec(i)/(1.2_dp*(extra_params(3)*extra_params(2))**2._dp*pi*extra_params(3)**2._dp)
    enddo
    close(10)
  end subroutine lift2file

  subroutine drag2file(dragvec,filename,extra_params)
    real(dp), intent(in), dimension(:) :: dragvec
    character(*), intent(in) :: filename
    real(dp), intent(in), dimension(:) :: extra_params ![1]dt [2]chord [3]span [4]speed
    integer :: i

    open(unit=10,file=filename)
    write(10,*) '# Drag'
    do i=1,size(dragvec,1)
      write(10,*) extra_params(1)*i,dragvec(i)!/(0.5_dp*1.2_dp*extra_params(4)**2._dp*extra_params(2)*extra_params(3))
    enddo
    close(10)
  end subroutine drag2file
end module postproc

