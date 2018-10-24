module postproc
  use rotor_classdef

contains

  subroutine rotor2file(rotor,timestamp)
    type(rotor_class), intent(in) :: rotor
    character(len=*), intent(in) :: timestamp
    character(len=5) :: nx_char, ny_char
    real(dp), dimension(3,rotor%nc+1,rotor%ns+1) :: wing_mesh  
    real(dp), dimension(3,rotor%nNwake+1,rotor%ns+1) :: wake_mesh  
    real(dp), dimension(3,rotor%nFwake+1) :: wake_tip   ! Optimise this by only initialising reqd size
    integer :: i,j,nx,ny,ib

    if (rotor%row_far .eq. 0) error stop "ERROR: plot only after far wake is created"

    open(unit=10,file='Results/Nwake'//timestamp//'.tec',position='append')
    open(unit=11,file='Results/Fwake'//timestamp//'.tec',position='append')

    write(10,*) 'Title = "Wing and Near wake"'
    write(10,*) 'VARIABLES = "X" "Y" "Z" "GAM"'

    do ib=1,rotor%nb
      nx=rotor%nc
      ny=rotor%ns
      write(nx_char,'(I5)') nx+1
      write(ny_char,'(I5)') ny+1

      do j=1,ny
        do i=1,nx
          wing_mesh(:,i,j)=rotor%blade(ib)%wiP(i,j)%pc(:,1)
        enddo
      enddo
      do i=1,nx
        wing_mesh(:,i,ny+1)=rotor%blade(ib)%wiP(i,ny)%pc(:,4)
      enddo
      do j=1,ny
        wing_mesh(:,nx+1,j)=rotor%blade(ib)%wiP(nx,j)%pc(:,2)
      enddo
      wing_mesh(:,nx+1,ny+1)=rotor%blade(ib)%wiP(nx,ny)%pc(:,3)

      write(10,*) 'Zone I='//trim(nx_char)//' J='//trim(ny_char)//' K=1  T="Blade"'
      write(10,*) 'DATAPACKING=BLOCK'
      write(10,*) 'VARLOCATION=([4]=CELLCENTERED)'!,[5]=CELLCENTERED)'
      write(10,*) ((wing_mesh(1,i,j),i=1,nx+1),j=1,ny+1)
      write(10,*) ((wing_mesh(2,i,j),i=1,nx+1),j=1,ny+1)
      write(10,*) ((wing_mesh(3,i,j),i=1,nx+1),j=1,ny+1)
      write(10,*) ((-1._dp*rotor%blade(ib)%wiP(i,j)%vr%gam,i=1,nx),j=1,ny)

      ! Near wake 
      nx=rotor%nNwake
      ny=rotor%ns
      write(nx_char,'(I5)') nx-(rotor%row_near-1)+1
      write(ny_char,'(I5)') ny+1

      !Check if necessary - $omp parallel do collapse(2)
      do j=1,ny
        do i=rotor%row_near,nx
          wake_mesh(:,i,j)=rotor%blade(ib)%waP(i,j)%vr%vf(1)%fc(:,1)
        enddo
      enddo
      !Check if necessary -$omp end parallel do
      do i=rotor%row_near,nx
        wake_mesh(:,i,ny+1)=rotor%blade(ib)%waP(i,ny)%vr%vf(4)%fc(:,1)
      enddo
      do j=1,ny
        wake_mesh(:,nx+1,j)=rotor%blade(ib)%waP(nx,j)%vr%vf(2)%fc(:,1)
      enddo
      wake_mesh(:,nx+1,ny+1)=rotor%blade(ib)%waP(nx,ny)%vr%vf(3)%fc(:,1)

      write(10,*) 'Zone I='//trim(nx_char)//' J='//trim(ny_char)//' K=1  T="NearWake"'
      write(10,*) 'DATAPACKING=BLOCK'
      write(10,*) 'VARLOCATION=([4]=CELLCENTERED)'!,[5]=CELLCENTERED)'
      write(10,*) ((wake_mesh(1,i,j),i=rotor%row_near,nx+1),j=1,ny+1)
      write(10,*) ((wake_mesh(2,i,j),i=rotor%row_near,nx+1),j=1,ny+1)
      write(10,*) ((wake_mesh(3,i,j),i=rotor%row_near,nx+1),j=1,ny+1)
      write(10,*) ((-1._dp*rotor%blade(ib)%waP(i,j)%vr%gam,i=rotor%row_near,nx),j=1,ny)

      ! Far wake 
      nx=rotor%nFwake
      write(nx_char,'(I5)') nx-(rotor%row_far-1)+1

      !Check if necessary - $omp parallel do collapse(2)
      do i=rotor%row_far,nx
        wake_tip(:,i)=rotor%blade(ib)%waF(i)%vf%fc(:,2)
      enddo
      wake_tip(:,nx+1)=rotor%blade(ib)%waF(rotor%nFwake)%vf%fc(:,1)
      !Check if necessary -$omp end parallel do

      write(11,*) 'Title = "Far wake"'
      write(11,*) 'VARIABLES = "X" "Y" "Z"'
      write(11,*) 'Zone I='//trim(nx_char)//' J=1   K=1   T="FarWake"'
      write(11,*) 'DATAPACKING=BLOCK'
      write(11,*) (wake_tip(1,i),i=rotor%row_far,nx+1)
      write(11,*) (wake_tip(2,i),i=rotor%row_far,nx+1)
      write(11,*) (wake_tip(3,i),i=rotor%row_far,nx+1)

    enddo

    close(10)
    close(11)
  end subroutine rotor2file

  subroutine filaments2file(rotor,timestamp)
    type(rotor_class), intent(in), dimension(:) :: rotor
    character(len=*), intent(in) :: timestamp

    integer :: nvr_wing, nvr_Nwake, nvf_Fwake
    integer :: ir, ib, irow, icol
    type(vr_class), allocatable, dimension(:) :: vr_wing, vr_Nwake
    type(vf_class), allocatable, dimension(:) :: vf_Fwake
    real(dp), allocatable, dimension(:) :: gam_Fwake

    nvr_wing=0
    nvr_Nwake=0
    nvf_Fwake=0

    nr=size(rotor)

    do ir=1,nr
      if (rotor(ir) .eq. 0) error stop 'ERROR: Use filaments2file() only after development of far wake'
    enddo

    ! Compute number of each filaments
    do ir=1,nr
      nvr_wing=nvr_wing+rotor(ir)%nb*(rotor(ir)%nc*rotor(ir)%ns)
      nvr_Nwake=nvr_Nwake+rotor(ir)%nb*(rotor(ir)%nNwake*rotor(ir)%ns)
      nvf_Fwake=nvf_Fwake+(rotor(ir)%nFwake-rotor(ir)%row_far+1)*rotor(ir)%nb
    enddo

    ! Allocate filaments
    allocate(vr_wing(nvr_wing))
    allocate(vr_Nwake(nvr_Nwake))
    allocate(vf_Fwake(nvf_Fwake))
    allocate(gam_Fwake(nvf_Fwake))

    ! Extract filament properties
    ! from wing
    indx=1;
    do ir=1,nr
      do ib=1,rotor(ir)%nb
        do icol=1,rotor(ir)%ns
          do irow=1,rotor(ir)%nc
            vr_wing(indx)=rotor(ir)%blade(ib)%wiP%(irow,icol)%vr
          enddo
        enddo
        indx=indx+1
      enddo
    enddo

    ! from Nwake
    indx=1;
    do ir=1,nr
      do ib=1,rotor(ir)%nb
        do icol=1,rotor(ir)%ns
          do irow=1,rotor(ir)%nNwake
            vr_Nwake(indx)=rotor(ir)%blade(ib)%waP%(irow,icol)%vr
          enddo
        enddo
        indx=indx+1
      enddo
    enddo


    ! from Fwake
    indx=1;
    do ir=1,nr
      do ib=1,rotor(ir)%nb
        do irow=rotor(ir)%row_far,rotor(ir)%nFwake
          vf_Fwake(indx)=rotor(ir)%blade(ib)%waF(irow)%vf
          gam_Fwake(indx)=rotor(ir)%blade(ib)%waF(irow)%gam
        enddo
        indx=indx+1
      enddo
    enddo

    ! Write to filamentsXXXXX.dat binary file
    open(unit=10,file='Results/filaments'//timestamp//'.dat',form='unformatted')
    write(10) nvr_wing
    write(10) nvr_Nwake
    write(10) nvf_Fwake
    write(10) vr_wing, vr_Nwake
    write(10) vf_Fwake, gam_Fwake
    close(10)

    ! Deallocate filaments
    deallocate(vr_wing)
    deallocate(vr_Nwake)
    deallocate(vf_Fwake)
    deallocate(gam_Fwake)

  end subroutine filaments2file

  subroutine mesh2file(wing_array,wake_array,filename)
    type(wingpanel_class), intent(in), dimension(:,:) :: wing_array
    type(Nwake_class), intent(in), dimension(:,:) :: wake_array
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
    write(10,*) 'VARIABLES = "X" "Y" "Z" "GAM"'
    write(10,*) 'Zone I='//trim(nx_char)//' J='//trim(ny_char)//' K=1  T="Wing"'
    write(10,*) 'DATAPACKING=BLOCK'
    write(10,*) 'VARLOCATION=([4]=CELLCENTERED)'!,[5]=CELLCENTERED)'
    write(10,*) ((wing_mesh(1,i,j),i=1,nx+1),j=1,ny+1)
    write(10,*) ((wing_mesh(2,i,j),i=1,nx+1),j=1,ny+1)
    write(10,*) ((wing_mesh(3,i,j),i=1,nx+1),j=1,ny+1)
    write(10,*) ((-1._dp*wing_array(i,j)%vr%gam,i=1,nx),j=1,ny)

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
    write(10,*) 'VARLOCATION=([4]=CELLCENTERED)'!,[5]=CELLCENTERED)'
    write(10,*) ((wake_mesh(1,i,j),i=1,nx+1),j=1,ny+1)
    write(10,*) ((wake_mesh(2,i,j),i=1,nx+1),j=1,ny+1)
    write(10,*) ((wake_mesh(3,i,j),i=1,nx+1),j=1,ny+1)
    write(10,*) ((-1._dp*wake_array(i,j)%vr%gam,i=1,nx),j=1,ny)

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
        wing_mesh(:,i,j)=wing_array(i,j)%CP
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
    type(Nwake_class), intent(in), dimension(:,:) :: wake_array
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
    write(10,*) 'VARIABLES = "X" "Y" "Z"'
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

