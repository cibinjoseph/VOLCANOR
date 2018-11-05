program gridgen
  use library

  integer :: nx,ny,nz
  real(dp), dimension(3) :: Cmin, Cmax    ! Coordinates of corners
  integer :: filerange_start, filerange_step, filerange_end

  integer :: ix,iy,iz,ifil
  real(dp), allocatable, dimension(:) :: xvec,yvec,zvec
  real(dp), allocatable, dimension(:,:,:,:) :: grid, grid_centre, vel
  character(len=5) :: nx_char,ny_char,nz_char
  character(len=5) :: timestamp

  integer :: nvr_wing, nvr_Nwake, nvf_Fwake
  type(vr_class), allocatable, dimension(:) :: vr_wing, vr_Nwake
  type(vf_class), allocatable, dimension(:) :: vf_Fwake
  real(dp), allocatable, dimension(:) :: gam_Fwake

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
  call print_status('Reading file '//'filaments'//timestamp//'.dat')
  open(unit=12,file='Results/filaments'//timestamp//'.dat',form='unformatted')
  read(12) nvr_wing
  read(12) nvr_Nwake
  read(12) nvf_Fwake

  allocate(vr_wing(nvr_wing))
  allocate(vr_Nwake(nvr_Nwake))
  allocate(vf_Fwake(nvf_Fwake))
  allocate(gam_Fwake(nvf_Fwake))

  read(12) vr_wing, vr_Nwake
  read(12) vf_Fwake, gam_Fwake
  close(12)
  call print_status()    ! SUCCESS

  ! Allocate grid coordinates
  allocate(grid(3,nx,ny,nz))
  allocate(grid_centre(3,nx-1,ny-1,nz-1))
  allocate(vel(3,nx-1,ny-1,nz-1))
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
  call print_status('Creating cartesian grid')
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
  call print_status()    ! SUCCESS

  ! Find induced velocities
  call print_status('Computing velocities')
  ! at cell centre
  vel=0._dp
  !$omp parallel do collapse(3)
  do iz=1,nz-1
    do iy=1,ny-1
      do ix=1,nx-1
        ! from wing
        do ifil=1,nvr_wing
          vel(:,ix,iy,iz)=vel(:,ix,iy,iz)+vr_wing(ifil)%vind(grid_centre(:,ix,iy,iz))*vr_wing(ifil)%gam
        enddo
        ! from Nwake
        do ifil=1,nvr_Nwake
          vel(:,ix,iy,iz)=vel(:,ix,iy,iz)+vr_Nwake(ifil)%vind(grid_centre(:,ix,iy,iz))*vr_Nwake(ifil)%gam
        enddo
        ! from Fwake
        do ifil=1,nvf_Fwake
          vel(:,ix,iy,iz)=vel(:,ix,iy,iz)+vf_Fwake(ifil)%vind(grid_centre(:,ix,iy,iz))*gam_Fwake(ifil)
        enddo
      enddo
    enddo
  enddo
  !$omp end parallel do
  call print_status()

  ! Write to file
  call print_status('Writing to grid file')
  write(timestamp,'(I0.5)') filerange_start
  open(unit=13,file='Results/grid'//timestamp//'.plt')

  write(13,*) 'TITLE = "Grid"'
  write(13,*) 'VARIABLES = "X" "Y" "Z" "U-vel" "V-vel" "W-vel"'
  write(13,*) 'ZONE I='//trim(nx_char)//' J='//trim(ny_char)//' K='//trim(nz_char)//' T="Data"'
  write(13,*) 'DATAPACKING=BLOCK'
  write(13,*) 'VARLOCATION=([4]=CELLCENTERED,[5]=CELLCENTERED,[6]=CELLCENTERED)'
  write(13,*) (((grid(1,ix,iy,iz),ix=1,nx),iy=1,ny),iz=1,nz) 
  write(13,*) (((grid(2,ix,iy,iz),ix=1,nx),iy=1,ny),iz=1,nz) 
  write(13,*) (((grid(3,ix,iy,iz),ix=1,nx),iy=1,ny),iz=1,nz) 
  write(13,*) (((vel(1,ix,iy,iz),ix=1,nx-1),iy=1,ny-1),iz=1,nz-1)
  write(13,*) (((vel(2,ix,iy,iz),ix=1,nx-1),iy=1,ny-1),iz=1,nz-1)
  write(13,*) (((vel(3,ix,iy,iz),ix=1,nx-1),iy=1,ny-1),iz=1,nz-1)
  close(13)
  call print_status()

  deallocate(vr_wing)
  deallocate(vr_Nwake)
  deallocate(vf_Fwake)
  deallocate(gam_Fwake)

  deallocate(grid)
  deallocate(grid_centre)
  deallocate(vel)
  deallocate(xvec)
  deallocate(yvec)
  deallocate(zvec)

end program gridgen
