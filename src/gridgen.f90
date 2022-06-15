program gridgen
  use libCommon

  integer :: nx, ny, nz
  real(dp), dimension(3) :: xyzMin, xyzMax    ! Coordinates of corners
  real(dp), dimension(3) :: vel           ! x,y,z velocities
  integer :: fileRange, fileRangeStart, fileRangeStep, fileRangeEnd

  integer :: ix, iy, iz, ifil
  real(dp), allocatable, dimension(:) :: xVec, yVec, zVec
  real(dp), allocatable, dimension(:, :, :, :) :: grid, gridCentre, velCentre
  character(len=5) :: nx_char, ny_char, nz_char
  character(len=5) :: filetimestamp

  integer :: nVrWing, nVrNwake, nVfNwakeTE, nVfFwake
  type(vr_class), allocatable, dimension(:) :: vrWing, vrNwake
  type(vf_class), allocatable, dimension(:) :: vfFwake, vfNwakeTE
  real(dp), allocatable, dimension(:) :: gamFwake, gamNwakeTE
  character(len=10) :: fileFormatVersion, currentTemplateVersion

  currentTemplateVersion = '0.2'

  ! Read gridconfig.in file
  call print_status('Reading file '//'gridconfig.in')
  open(unit=11, file='gridconfig.in', status='old', action='read')

  namelist /VERSION/ fileFormatVersion
  read(unit=11, nml=VERSION)
  if (adjustl(fileFormatVersion) /= currentTemplateVersion) then
    error stop 'ERROR: gridconfig.nml template version does not match'
  endif

  namelist /INPUTS/ nx, ny, nz, xyzMin, xyzMax, vel, &
    & fileRangeStart, fileRangeStep, fileRangeEnd
  read(unit=11, nml=INPUTS)
  close(11)

  ! Sanity check for xyzMin and xyzMax values
  if (xyzMin(1) > xyzMax(1) .or. xyzMin(2) > xyzMax(2) .or. xyzMin(3) > xyzMax(3)) then
    error stop 'ERROR: All XYZmin values should be greater than XYZmax values'
  endif

  call print_status()    ! SUCCESS

  ! Allocate grid coordinates
  allocate (grid(3, nx, ny, nz))
  allocate (gridCentre(3, nx - 1, ny - 1, nz - 1))
  allocate (velCentre(3, nx - 1, ny - 1, nz - 1))
  allocate (xVec(nx))
  allocate (yVec(ny))
  allocate (zVec(nz))

  write (nx_char, '(I5)') nx
  write (ny_char, '(I5)') ny
  write (nz_char, '(I5)') nz

  xVec = linspace(xyzMin(1), xyzMax(1), nx)
  yVec = linspace(xyzMin(2), xyzMax(2), ny)
  zVec = linspace(xyzMin(3), xyzMax(3), nz)

  ! Create grid
  call print_status('Creating cartesian grid')
  do iz = 1, nz
    do iy = 1, ny
      do ix = 1, nx
        grid(:, ix, iy, iz) = [xVec(ix), yVec(iy), zVec(iz)]
      enddo
    enddo
  enddo

  ! Compute grid centres
  do iz = 1, nz - 1
    do iy = 1, ny - 1
      do ix = 1, nx - 1
        gridCentre(:, ix, iy, iz) = grid(:, ix, iy, iz) + grid(:, ix + 1, iy, iz) &
          + grid(:, ix + 1, iy + 1, iz) + grid(:, ix + 1, iy + 1, iz + 1) &
          + grid(:, ix, iy + 1, iz) + grid(:, ix, iy + 1, iz + 1) &
          + grid(:, ix, iy, iz + 1) + grid(:, ix + 1, iy, iz + 1)
      enddo
    enddo
  enddo
  gridCentre = gridCentre*0.125_dp
  call print_status()    ! SUCCESS

  ! Iterate through filaments files
  do fileRange = fileRangeStart, fileRangeEnd, fileRangeStep
    ! Read from filaments file
    write (filetimestamp, '(I0.5)') fileRange
    call print_status('Reading file '//'filaments'//filetimestamp//'.dat')
    open (unit=12, file='Results/filaments'//filetimestamp//'.dat', &
      & status='old', action='read', form='unformatted')
    read (12) nVrWing
    read (12) nVrNwake
    read (12) nVfNwakeTE
    read (12) nVfFwake

    allocate (vrWing(nVrWing))
    allocate (vrNwake(nVrNwake))
    allocate (vfNwakeTE(nVfNwakeTE))
    allocate (vfFwake(nVfFwake))
    allocate (gamFwake(nVfFwake))
    allocate (gamNwakeTE(nVfNwakeTE))

    read (12) vrWing, vrNwake
    read (12) vfNwakeTE, gamNwakeTE
    read (12) vfFwake, gamFwake
    close (12)
    call print_status()    ! SUCCESS

    ! Find induced velocities at cell centre
    call print_status('Computing velocities')
    !$omp parallel do collapse(3)
    do iz = 1, nz - 1
      do iy = 1, ny - 1
        do ix = 1, nx - 1
          ! from wing
          do ifil = 1, nVrWing
            velCentre(:, ix, iy, iz) = vrWing(ifil)%vind(gridCentre(:, ix, iy, iz))*vrWing(ifil)%gam
          enddo
          ! from Nwake
          do ifil = 1, nVrNwake
            velCentre(:, ix, iy, iz) = velCentre(:, ix, iy, iz) + vrNwake(ifil)%vind(gridCentre(:, ix, iy, iz))*vrNwake(ifil)%gam
          enddo
          ! from NwakeTE
          do ifil = 1, nVfNwakeTE
            velCentre(:, ix, iy, iz) = velCentre(:, ix, iy, iz) + vfNwakeTE(ifil)%vind(gridCentre(:, ix, iy, iz))*gamNwakeTE(ifil)
          enddo
          ! from Fwake
          do ifil = 1, nVfFwake
            velCentre(:, ix, iy, iz) = velCentre(:, ix, iy, iz) + vfFwake(ifil)%vind(gridCentre(:, ix, iy, iz))*gamFwake(ifil)
          enddo
        enddo
      enddo
    enddo
    !$omp end parallel do
    call print_status()

    ! Add freestream velocities
    velCentre(1, :, :, :) = velCentre(1, :, :, :) + vel(1)
    velCentre(2, :, :, :) = velCentre(2, :, :, :) + vel(2)
    velCentre(3, :, :, :) = velCentre(3, :, :, :) + vel(3)

    ! Write to file
    call print_status('Writing file '//'grid'//filetimestamp//'.plt')
    write (filetimestamp, '(I0.5)') fileRange
    open (unit=13, file='Results/grid'//filetimestamp//'.plt', &
      & status='old', action='read')

    write (13, *) 'TITLE = "Grid"'
    write (13, *) 'VARIABLES = "X" "Y" "Z" "U" "V" "W"'
    write (13, *) 'ZONE I='//trim(nx_char)//' J='//trim(ny_char)//' K='//trim(nz_char)//' T="Data"'
    write (13, *) 'DATAPACKING=BLOCK'
    write (13, *) 'VARLOCATION=([4]=CELLCENTERED,[5]=CELLCENTERED,[6]=CELLCENTERED)'
    write (13, *) (((grid(1, ix, iy, iz), ix=1, nx), iy=1, ny), iz=1, nz)
    write (13, *) (((grid(2, ix, iy, iz), ix=1, nx), iy=1, ny), iz=1, nz)
    write (13, *) (((grid(3, ix, iy, iz), ix=1, nx), iy=1, ny), iz=1, nz)
    write (13, *) (((velCentre(1, ix, iy, iz), ix=1, nx - 1), iy=1, ny - 1), iz=1, nz - 1)
    write (13, *) (((velCentre(2, ix, iy, iz), ix=1, nx - 1), iy=1, ny - 1), iz=1, nz - 1)
    write (13, *) (((velCentre(3, ix, iy, iz), ix=1, nx - 1), iy=1, ny - 1), iz=1, nz - 1)
    close (13)
    call print_status()

    deallocate (vrWing)
    deallocate (vrNwake)
    deallocate (vfNwakeTE)
    deallocate (vfFwake)
    deallocate (gamFwake)
    deallocate (gamNwakeTE)
  enddo

  deallocate (grid)
  deallocate (gridCentre)
  deallocate (velCentre)
  deallocate (xVec)
  deallocate (yVec)
  deallocate (zVec)

end program gridgen
