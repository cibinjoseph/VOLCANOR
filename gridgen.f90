program gridgen
  use library

  integer :: nx,ny,nz
  real(dp), dimension(3) :: cMin, cMax    ! Coordinates of corners
  real(dp), dimension(3) :: vel           ! x,y,z velocities
  integer :: fileRange, fileRangeStart, fileRangeStep, fileRangeEnd

  integer :: ix,iy,iz,ifil
  real(dp), allocatable, dimension(:) :: xVec,yVec,zVec
  real(dp), allocatable, dimension(:,:,:,:) :: grid, gridCentre, velCentre
  character(len=5) :: nx_char,ny_char,nz_char
  character(len=5) :: timestamp

  integer :: nVrWing, nVrNwake, nVfFwake
  type(vr_class), allocatable, dimension(:) :: vrWing, vrNwake
  type(vf_class), allocatable, dimension(:) :: vfFwake
  real(dp), allocatable, dimension(:) :: gamFwake

  ! Read gridconfig.in file
  call print_status('Reading file '//'gridconfig.in')
  open(unit=11,file='gridconfig.in')
  call skiplines(11,3)
  read(11,*) nx,ny,nz
  call skiplines(11,4)
  read(11,*) cMin(1),cMin(2),cMin(3)
  call skiplines(11,3)
  read(11,*) cMax(1),cMax(2),cMax(3)
  call skiplines(11,4)
  read(11,*) vel(1),vel(2),vel(3)
  call skiplines(11,5)
  read(11,*) fileRangeStart, fileRangeStep, fileRangeEnd
  close(11)

  ! Sanity check for cMin and cMax values
  if (cMin(1)>cMax(1) .or. cMin(2)>cMax(2) .or. cMin(3)>cMax(3)) then
    error stop 'ERROR: All XYZmin values should be greater than XYZmax values'
  endif

  call print_status()    ! SUCCESS

  ! Allocate grid coordinates
  allocate(grid(3,nx,ny,nz))
  allocate(gridCentre(3,nx-1,ny-1,nz-1))
  allocate(velCentre(3,nx-1,ny-1,nz-1))
  allocate(xVec(nx))
  allocate(yVec(ny))
  allocate(zVec(nz))

  write(nx_char,'(I5)') nx
  write(ny_char,'(I5)') ny
  write(nz_char,'(I5)') nz

  xVec=linspace(cMin(1),cMax(1),nx)
  yVec=linspace(cMin(2),cMax(2),ny)
  zVec=linspace(cMin(3),cMax(3),nz)

  ! Create grid
  call print_status('Creating cartesian grid')
  do iz=1,nz
    do iy=1,ny
      do ix=1,nx
        grid(:,ix,iy,iz)=(/xVec(ix),yVec(iy),zVec(iz)/)
      enddo
    enddo
  enddo

  ! Compute grid centres
  do iz=1,nz-1
    do iy=1,ny-1
      do ix=1,nx-1
        gridCentre(:,ix,iy,iz)=grid(:,ix,iy,iz)+grid(:,ix+1,iy,iz)+grid(:,ix+1,iy+1,iz)+grid(:,ix+1,iy+1,iz+1)  &
          + grid(:,ix,iy+1,iz)+grid(:,ix,iy+1,iz+1)+grid(:,ix,iy,iz+1)+grid(:,ix+1,iy,iz+1)
      enddo
    enddo
  enddo
  gridCentre=gridCentre*0.125_dp
  call print_status()    ! SUCCESS

  ! Iterate through filaments files
  do fileRange = fileRangeStart, fileRangeEnd, fileRangeStep
    ! Read from filaments file
    write(timestamp,'(I0.5)') fileRange
    call print_status('Reading file '//'filaments'//timestamp//'.dat')
    open(unit=12,file='Results/filaments'//timestamp//'.dat',form='unformatted')
    read(12) nVrWing
    read(12) nVrNwake
    read(12) nVfFwake

    allocate(vrWing(nVrWing))
    allocate(vrNwake(nVrNwake))
    allocate(vfFwake(nVfFwake))
    allocate(gamFwake(nVfFwake))

    read(12) vrWing, vrNwake
    read(12) vfFwake, gamFwake
    close(12)
    call print_status()    ! SUCCESS

    ! Find induced velocities at cell centre
    call print_status('Computing velocities')
    !$omp parallel do collapse(3)
    do iz=1,nz-1
      do iy=1,ny-1
        do ix=1,nx-1
          ! from wing
          do ifil=1,nVrWing
            velCentre(:,ix,iy,iz)=vrWing(ifil)%vind(gridCentre(:,ix,iy,iz))*vrWing(ifil)%gam
          enddo
          ! from Nwake
          do ifil=1,nVrNwake
            velCentre(:,ix,iy,iz)=velCentre(:,ix,iy,iz)+vrNwake(ifil)%vind(gridCentre(:,ix,iy,iz))*vrNwake(ifil)%gam
          enddo
          ! from Fwake
          do ifil=1,nVfFwake
            velCentre(:,ix,iy,iz)=velCentre(:,ix,iy,iz)+vfFwake(ifil)%vind(gridCentre(:,ix,iy,iz))*gamFwake(ifil)
          enddo
        enddo
      enddo
    enddo
    !$omp end parallel do
    call print_status()

    ! Add freestream velocities
    velCentre(1,:,:,:)=velCentre(1,:,:,:)+vel(1)
    velCentre(2,:,:,:)=velCentre(2,:,:,:)+vel(2)
    velCentre(3,:,:,:)=velCentre(3,:,:,:)+vel(3)

    ! Write to file
    call print_status('Writing file '//'grid'//timestamp//'.plt')
    write(timestamp,'(I0.5)') fileRange
    open(unit=13,file='Results/grid'//timestamp//'.plt')

    write(13,*) 'TITLE = "Grid"'
    write(13,*) 'VARIABLES = "X" "Y" "Z" "U" "V" "W"'
    write(13,*) 'ZONE I='//trim(nx_char)//' J='//trim(ny_char)//' K='//trim(nz_char)//' T="Data"'
    write(13,*) 'DATAPACKING=BLOCK'
    write(13,*) 'VARLOCATION=([4]=CELLCENTERED,[5]=CELLCENTERED,[6]=CELLCENTERED)'
    write(13,*) (((grid(1,ix,iy,iz),ix=1,nx),iy=1,ny),iz=1,nz) 
    write(13,*) (((grid(2,ix,iy,iz),ix=1,nx),iy=1,ny),iz=1,nz) 
    write(13,*) (((grid(3,ix,iy,iz),ix=1,nx),iy=1,ny),iz=1,nz) 
    write(13,*) (((velCentre(1,ix,iy,iz),ix=1,nx-1),iy=1,ny-1),iz=1,nz-1)
    write(13,*) (((velCentre(2,ix,iy,iz),ix=1,nx-1),iy=1,ny-1),iz=1,nz-1)
    write(13,*) (((velCentre(3,ix,iy,iz),ix=1,nx-1),iy=1,ny-1),iz=1,nz-1)
    close(13)
    call print_status()

    deallocate(vrWing)
    deallocate(vrNwake)
    deallocate(vfFwake)
    deallocate(gamFwake)
  enddo

  deallocate(grid)
  deallocate(gridCentre)
  deallocate(velCentre)
  deallocate(xVec)
  deallocate(yVec)
  deallocate(zVec)


end program gridgen
