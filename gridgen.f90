program gridgen
  use library

  integer :: nx,ny,nz
  real(dp), dimension(3) :: Cmin, Cmax    ! Coordinates of corners
  integer :: filerange_start, filerange_step, filerange_end

  integer :: ix,iy,iz,icrd
  real(dp), allocatable, dimension(:) :: xvec,yvec,zvec
  real(dp), allocatable, dimension(:,:,:,:) :: grid, grid_centre
  character(len=5) :: nx_char,ny_char,nz_char
  character(len=5) :: timestamp

  ! Read gridconfig.in file
  call print_status('Reading file '//'gridconfig.in')
  open(unit=11,file='gridconfig.in')
  call skiplines(11,3)
  read(11,*) nx,ny,nz
  call skiplines(11,4)
  read(11,*) Cmin(1),Cmin(2),Cmin(3)
  call skiplines(11,3)
  read(11,*) Cmax(1),Cmax(2),Cmax(3)
  call skiplines(11,5)
  read(11,*) filerange_start, filerange_step, filerange_end
  close(11)

  ! Sanity check for Cmin and Cmax values
  if (Cmin(1)>Cmax(1) .or. Cmin(2)>Cmax(2) .or. Cmin(3)>Cmax(3)) then
    error stop 'ERROR: All XYZmin values should be greater than XYZmax values'
  endif

  call print_status()    ! SUCCESS

  ! Read from filaments file
  write(timestamp,'(I0.5)') filerange_start
  open(unit=12,file='Results/filaments'//timestamp//'.dat')

  close(12)

  ! Allocate grid coordinates
  allocate(grid(3,nx,ny,nz))
  allocate(grid_centre(3,nx-1,ny-1,nz-1))
  allocate(xvec(nx))
  allocate(yvec(ny))
  allocate(zvec(nz))

  write(nx_char,'(I5)') nx
  write(ny_char,'(I5)') ny
  write(nz_char,'(I5)') nz

  xvec=linspace(Cmin(1),Cmax(1),nx)
  yvec=linspace(Cmin(2),Cmax(2),ny)
  zvec=linspace(Cmin(3),Cmax(3),nz)

  ! Create grid
  do iz=1,nz
    do iy=1,ny
      do ix=1,nx
        grid(:,ix,iy,iz)=(/xvec(ix),yvec(iy),zvec(iz)/)
      enddo
    enddo
  enddo

  ! Compute grid centres
  do iz=1,nz-1
    do iy=1,ny-1
      do ix=1,nx-1
        grid_centre(:,ix,iy,iz)=grid(:,ix,iy,iz)+grid(:,ix+1,iy,iz)+grid(:,ix+1,iy+1,iz)+grid(:,ix+1,iy+1,iz+1)  &
          + grid(:,ix,iy+1,iz)+grid(:,ix,iy+1,iz+1)+grid(:,ix,iy,iz+1)+grid(:,ix+1,iy,iz+1)
      enddo
    enddo
  enddo
  grid_centre=grid_centre*0.125_dp

  ! Write to file
  call print_status('Writing to grid file')
  write(timestamp,'(I0.5)') filerange_start
  open(unit=13,file='Results/grid'//timestamp//'.dat')

  write(13,*) 'TITLE = "Grid"'
  write(13,*) 'VARIABLES = "X" "Y" "Z" "U-vel" "V-vel" "W-vel"'
  write(13,*) 'ZONE I='//trim(nx_char)//' J='//trim(ny_char)//' K='//trim(nz_char)//' T="Data"'
  write(13,*) 'DATAPACKING=BLOCK'
  write(13,*) 'VARLOCATION=([4]=CELLCENTERED,[5]=CELLCENTERED,[6]=CELLCENTERED)'
  write(13,*) (((grid(1,ix,iy,iz),ix=1,nx),iy=1,ny),iz=1,nz) 
  write(13,*) (((grid(2,ix,iy,iz),ix=1,nx),iy=1,ny),iz=1,nz) 
  write(13,*) (((grid(3,ix,iy,iz),ix=1,nx),iy=1,ny),iz=1,nz) 
  write(13,*) (((grid_centre(1,ix,iy,iz),ix=1,nx-1),iy=1,ny-1),iz=1,nz-1)
  write(13,*) (((grid_centre(2,ix,iy,iz),ix=1,nx-1),iy=1,ny-1),iz=1,nz-1)
  write(13,*) (((grid_centre(3,ix,iy,iz),ix=1,nx-1),iy=1,ny-1),iz=1,nz-1)
  close(13)
  call print_status()

end program gridgen
