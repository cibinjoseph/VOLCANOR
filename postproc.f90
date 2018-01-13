module postproc
  use wingpanel_classdef
  use wakepanel_classdef

contains

  subroutine meshwrite2(wing_array,wake_array,filename)
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
    write(10,*) ((wing_array(i,j)%vr%gam,i=1,nx),j=1,ny)
    write(10,*) ((wing_array(i,j)%tag,i=1,nx),j=1,ny)

    nx=size(wake_array,1)
    ny=size(wake_array,2)
    write(nx_char,'(I5)') nx+1
    write(ny_char,'(I5)') ny+1

    do j=1,ny
      do i=1,nx
        wake_mesh(:,i,j)=wake_array(i,j)%vr%vf(1)%fc(:,1)
      enddo
    enddo
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
    write(10,*) ((wake_array(i,j)%vr%gam,i=1,nx),j=1,ny)
    write(10,*) ((wake_array(i,j)%tag,i=1,nx),j=1,ny)

    close(10)
  end subroutine meshwrite2

  subroutine tipwrite(wing_array,wake_array,filename)
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
  end subroutine tipwrite

end module postproc

