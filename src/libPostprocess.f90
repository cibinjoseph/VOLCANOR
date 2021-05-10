module libPostprocess
  use rotor_classdef
  use switches_classdef
  character(len=8) :: ResultsDir = 'Results/'

contains

  subroutine init_plots(numOfRotors)
    ! Initialise headers for plot files
    integer, intent(in) :: numOfRotors
    character(len=24) :: forceDimFilename
    character(len=27) :: forceNonDimFilename
    integer :: rotorNumber
    character(len=2) :: rotorNumberChar

    do rotorNumber = 1, numOfRotors
      write (rotorNumberChar, '(I0.2)') rotorNumber
      forceDimFilename = ResultsDir//'r'//rotorNumberChar//'ForceDim.csv'
      forceNonDimFilename = ResultsDir//'r'//rotorNumberChar//'ForceNonDim.csv'

      ! Add data headers
      open (unit=11, file=forceDimFilename, &
        & status='replace', action='write')
      write (11, 100) 'iter','LiftMag','DragMag', &
        & 'Lx','Ly','Lz', &
        & 'Dx','Dy','Dz', &
        & 'FInertx','FInerty','FInertz'
      close (11)

      open (unit=12, file=forceNonDimFilename, &
        & status='replace', action='write')
      write (12, 101) 'iter','CL/CT','CD/CQ','CDi','CD0','CDu', &
        & 'CFx','CFy','CFz'
      close (12)
    enddo
    100 format (A5,11(A15))
    101 format (A5,8(A15))
  end subroutine init_plots

  subroutine params2file(rotor, rotorNumber, nt, dt, nr, &
      & density, velSound, switches)
    ! Write rotor parameters to json file
    type(rotor_class), intent(in) :: rotor
    type(switches_class), intent(in) :: switches
    integer, intent(in) :: rotorNumber, nt, nr
    real(dp), intent(in) :: dt, density, velSound
    character(len=2) :: rotorNumberChar
    character(len=22) :: paramsFilename

    write (rotorNumberChar, '(I0.2)') rotorNumber
    paramsFilename = ResultsDir//'r'//rotorNumberChar//'Params.json'
    open(unit=10, file=paramsFilename, status='replace', action='write')
    write(10, *) '{'
    ! Config file 
    write(10, *) '"nt": ', nt, ','
    write(10, *) '"dt": ', dt, ','
    write(10, *) '"nr": ', nr, ','
    write(10, *) '"restartWriteNt": ', switches%restartWriteNt, ','
    write(10, *) '"restartFromNt": ', switches%restartFromNt, ','
    write(10, *) '"ntSub": ', switches%ntSub, ','
    write(10, *) '"ntSubInit": ', switches%ntSubInit, ','
    write(10, *) '"spanSpacing": ', switches%spanSpacing, ','
    write(10, *) '"density": ', density, ','
    write(10, *) '"velSound": ', velSound, ','
    write(10, *) '"wakePlot": ', switches%wakePlot, ','
    write(10, *) '"wakeTipPlot": ', switches%wakeTipPlot, ','
    write(10, *) '"rotorForcePlot": ', switches%rotorForcePlot, ','
    write(10, *) '"gridPlot": ', switches%gridPlot, ','
    write(10, *) '"wakeDissipation": ', switches%wakeDissipation, ','
    write(10, *) '"wakeStrain": ', switches%wakeStrain, ','
    write(10, *) '"wakeBurst": ', switches%wakeBurst, ','
    write(10, *) '"slowStart": ', switches%slowStart, ','
    write(10, *) '"slowStartNt": ', switches%slowStartNt, ','
    write(10, *) '"fdScheme": ', switches%fdScheme, ','
    write(10, *) '"initWakeVelNt": ', switches%initWakeVelNt, ','
    write(10, *) '"probe": ', switches%probe, ','

    ! Geom file 
    write(10, *) '"nb": ', rotor%nb, ','
    write(10, *) '"propConvention": ', rotor%propConvention, ','
    write(10, *) '"geometryFile": "', rotor%geometryFile, '",'
    write(10, *) '"nCamberFiles": "', rotor%nCamberFiles, '",'
    write(10, *) '"surfaceType": "', rotor%surfaceType, '",'
    write(10, *) '"nc": ', rotor%nc, ','
    write(10, *) '"ns": ', rotor%ns, ','
    write(10, *) '"nNwake": ', rotor%nNwake, ','
    write(10, *) '"radius": ', rotor%radius, ','
    write(10, *) '"root_cut": ', rotor%root_cut, ','
    write(10, *) '"chord": ', rotor%chord, ','
    write(10, *) '"coningAngle": ', rotor%coningAngle*180._dp/pi, ','
    write(10, *) '"Omega": ', rotor%Omega, ','
    write(10, *) '"phi": ', rotor%pts(1), ','
    write(10, *) '"theta": ', rotor%pts(2), ','
    write(10, *) '"psi": ', rotor%pts(3), ','
    write(10, *) '"theta0": ', rotor%controlPitch(1)*180._dp/pi, ','
    write(10, *) '"thetaC": ', rotor%controlPitch(2)*180._dp/pi, ','
    write(10, *) '"thetaS": ', rotor%controlPitch(3)*180._dp/pi, ','
    write(10, *) '"thetaTwist": ', rotor%thetaTwist*180._dp/pi, ','
    write(10, *) '"u": ', rotor%velBody(1), ','
    write(10, *) '"v": ', rotor%velBody(2), ','
    write(10, *) '"w": ', rotor%velBody(3), ','
    write(10, *) '"nonDimForceDenom": ', rotor%nonDimforceDenominator, ','
    if (rotor%nAirfoils .gt. 0) then
      write(10, *) '"alpha0": ', rotor%alpha0(1), ','
    else
      write(10, *) '"alpha0": ', 0._dp, ','
    endif
    write(10, *) '"wakeTruncateNt": ', rotor%wakeTruncateNt, ','
    write(10, *) '"initWakeVel": ', rotor%initWakeVel, ','
    write(10, *) '"psiStart": ', rotor%psiStart*180._dp/pi, ','
    write(10, *) '"forceCalcSwitch": ', rotor%forceCalcSwitch*180._dp/pi, ','
    if (switches%wakeDissipation .gt. 0) then
      write(10, *) '"turbulentViscosity": ', rotor%turbulentViscosity
    else
      write(10, *) '"turbulentViscosity": ', 0.0
    endif
    write(10, *) '}'
    close(10)
  end subroutine params2file

  subroutine rotor2file(timestamp, rotor, rotorNumber)
    ! Plot rotor geometry and wake to file
    type(rotor_class), intent(in) :: rotor
    character(len=*), intent(in) :: timestamp
    integer, intent(in) :: rotorNumber
    character(len=5) :: nxChar, nyChar
    character(len=2) :: rotorNumberChar, bladeNumberChar
    real(dp), dimension(3, rotor%nc + 1, rotor%ns + 1) :: wingMesh
    real(dp), dimension(3, rotor%nNwake + 1, rotor%ns + 1) :: wakeMesh
    real(dp), dimension(3, rotor%nFwake + 1) :: wakeTip   ! Optimise this by only initialising reqd size
    integer :: i, j, nx, ny, ib

    write (rotorNumberChar, '(I0.2)') rotorNumber

    open (unit=10, file=ResultsDir// &
      & 'r'//rotorNumberChar//'wingNwake'//timestamp//'.plt', position='append')

    if (abs(rotor%surfaceType) == 1) then
      ! Lifting surface
      write (10, *) 'Title = "Wing and Wake"'
      write (10, *) 'VARIABLES = "X" "Y" "Z" "GAM" "skew"' ! "Var6"'

      do ib = 1, rotor%nb
        write (bladeNumberChar, '(I0.2)') ib
        ! Wing
        nx = rotor%nc
        ny = rotor%ns
        write (nxChar, '(I5)') nx + 1
        write (nyChar, '(I5)') ny + 1

        do j = 1, ny
          do i = 1, nx
            wingMesh(:, i, j) = rotor%blade(ib)%wiP(i, j)%pc(:, 1)
          enddo
        enddo
        do i = 1, nx
          wingMesh(:, i, ny + 1) = rotor%blade(ib)%wiP(i, ny)%pc(:, 4)
        enddo
        do j = 1, ny
          wingMesh(:, nx + 1, j) = rotor%blade(ib)%wiP(nx, j)%pc(:, 2)
        enddo
        wingMesh(:, nx + 1, ny + 1) = rotor%blade(ib)%wiP(nx, ny)%pc(:, 3)

        write (10, *) 'Zone I='//trim(nxChar)//' J='//trim(nyChar)//' K=1  T="Blade '//trim(bladeNumberChar)//'"'
        write (10, *) 'DATAPACKING=BLOCK'
        write (10, *) 'VARLOCATION=([4]=CELLCENTERED,[5]=CELLCENTERED)' !,[6]=CELLCENTERED)'
        write (10, *) ((wingMesh(1, i, j), i=1, nx + 1), j=1, ny + 1)
        write (10, *) ((wingMesh(2, i, j), i=1, nx + 1), j=1, ny + 1)
        write (10, *) ((wingMesh(3, i, j), i=1, nx + 1), j=1, ny + 1)
        write (10, *) ((-1._dp*rotor%blade(ib)%wiP(i, j)%vr%gam, i=1, nx), j=1, ny)
        write(10, *) ((rotor%blade(ib)%wiP(i, j)%vr%gam*0._dp, i=1,nx), j=1,ny)
        !write(10,*) ((rotor%blade(ib)%wiP(i,j)%vr%skew,i=1,nx),j=1,ny)

        ! Near wake
        nx = rotor%nNwake
        ny = rotor%ns
        write (nxChar, '(I5)') nx - (rotor%rowNear - 1) + 1
        write (nyChar, '(I5)') ny + 1

        !Check if necessary - $omp parallel do collapse(2)
        do j = 1, ny
          do i = rotor%rowNear, nx
            wakeMesh(:, i, j) = rotor%blade(ib)%waP(i, j)%vr%vf(1)%fc(:, 1)
          enddo
        enddo
        !Check if necessary -$omp end parallel do
        do i = rotor%rowNear, nx
          wakeMesh(:, i, ny + 1) = rotor%blade(ib)%waP(i, ny)%vr%vf(4)%fc(:, 1)
        enddo
        do j = 1, ny
          wakeMesh(:, nx + 1, j) = rotor%blade(ib)%waP(nx, j)%vr%vf(2)%fc(:, 1)
        enddo
        wakeMesh(:, nx + 1, ny + 1) = rotor%blade(ib)%waP(nx, ny)%vr%vf(3)%fc(:, 1)

        write (10, *) 'Zone I='//trim(nxChar)//' J='//trim(nyChar)//' K=1  T="NearWake"'
        write (10, *) 'DATAPACKING=BLOCK'
        write (10, *) 'VARLOCATION=([4]=CELLCENTERED, [5]=CELLCENTERED)' !,[6]=CELLCENTERED)'
        write (10, *) ((wakeMesh(1, i, j), i=rotor%rowNear, nx + 1), j=1, ny + 1)
        write (10, *) ((wakeMesh(2, i, j), i=rotor%rowNear, nx + 1), j=1, ny + 1)
        write (10, *) ((wakeMesh(3, i, j), i=rotor%rowNear, nx + 1), j=1, ny + 1)
        write (10, *) ((-1._dp*rotor%blade(ib)%waP(i, j)%vr%gam, i=rotor%rowNear, nx), j=1, ny)
        write(10,*) ((rotor%blade(ib)%waP(i,j)%vr%skew, i=rotor%rowNear, nx), j=1, ny)
        !write(10,*) ((rotor%blade(ib)%waP(i,j)%vr%skew,i=rotor%rowNear,nx),j=1,ny)

        ! Far wake
        nx = rotor%nFwake
        if (rotor%rowFar .le. rotor%nFwake) then
          write (nxChar, '(I5)') nx - (rotor%rowFar - 1) + 1

          !Check if necessary - $omp parallel do collapse(2)
          do i = rotor%rowFar, nx
            wakeTip(:, i) = rotor%blade(ib)%waF(i)%vf%fc(:, 2)
          enddo
          wakeTip(:, nx + 1) = rotor%blade(ib)%waF(rotor%nFwake)%vf%fc(:, 1)
          !Check if necessary -$omp end parallel do

          write (10, *) 'Zone I='//trim(nxChar)//' J=1   K=1   T="FarWake"'
          write (10, *) 'DATAPACKING=BLOCK'
          write (10, *) 'VARLOCATION=([4]=CELLCENTERED ,[5]=CELLCENTERED)' !,[6]=CELLCENTERED)'
          write (10, *) (wakeTip(1, i), i=rotor%rowFar, nx + 1)
          write (10, *) (wakeTip(2, i), i=rotor%rowFar, nx + 1)
          write (10, *) (wakeTip(3, i), i=rotor%rowFar, nx + 1)
          write (10, *) (-1._dp*rotor%blade(ib)%waF(i)%gam, i=rotor%rowFar, nx)
          write(10,*) (rotor%blade(ib)%waF(i)%vf%rVc*0._dp, i=rotor%rowFar, nx)
          !write(10,*) (rotor%blade(ib)%waF(i)%vf%age,i=rotor%rowFar,nx)

        else  ! No far wake present

          write (nxChar, '(I5)') 2  ! Plot mesh as single redundant point
          wakeTip(:, 1) = rotor%blade(ib)%waP(rotor%nNwake, rotor%ns)%vr%vf(3)%fc(:, 1)

          write (10, *) 'Zone I='//trim(nxChar)//' J=1   K=1   T="FarWake"'
          write (10, *) 'DATAPACKING=BLOCK'
          write (10, *) 'VARLOCATION=([4]=CELLCENTERED, [5]=CELLCENTERED)' !,[6]=CELLCENTERED)'
          write (10, *) wakeTip(1, 1), wakeTip(1, 1)
          write (10, *) wakeTip(2, 1), wakeTip(2, 1)
          write (10, *) wakeTip(3, 1), wakeTip(3, 1)
          write (10, *) 0._dp
          write(10,*) 0._dp
          !write(10,*) 0._dp
        endif

      enddo
    elseif (abs(rotor%surfaceType) == 2) then
      ! Non-lifting surface
      write (10, *) 'Title = "Non-lifting surface"'
      write (10, *) 'VARIABLES = "X" "Y" "Z" "GAM"' ! "skew" "Var6"'

      do ib = 1, rotor%nb
        write (bladeNumberChar, '(I0.2)') ib

        ! Compute common nodes of non-lifting surface
        nx = rotor%blade(ib)%stlNodesCols
        ny = rotor%nc ! No. of triangular elements
        write (nxChar, '(I0.5)') nx
        write (nyChar, '(I0.5)') ny

        write(10, *) 'Zone NODES='//trim(nxChar)// &
          & ' ELEMENTS='//nyChar// &
          & ' T="Blade '//trim(bladeNumberChar)//'"'
        write(10, *) 'ZONETYPE=FETRIANGLE, DATAPACKING=BLOCK'
        write(10, *) 'VARLOCATION=(4=CELLCENTERED)'
        write(10, *) rotor%blade(ib)%stlNodes(1, 1:nx)
        write(10, *) rotor%blade(ib)%stlNodes(2, 1:nx)
        write(10, *) rotor%blade(ib)%stlNodes(3, 1:nx)
        write(10, *) rotor%blade(ib)%wiP(:, 1)%vr%gam
        do i = 1, ny
          write(10, '(3I7)') rotor%blade(ib)%stlElementNodes(:, i)
        enddo
      enddo
    endif
    close (10)
  end subroutine rotor2file

  subroutine probes2file(timestamp, probe, probeVel, rotor, t)
    ! Write velocities at probe locations
    type(rotor_class), intent(in), dimension(:) :: rotor
    character(len=*), intent(in) :: timestamp
    real(dp), intent(in) :: t
    real(dp), intent(in) , dimension(:, :) :: probe, probeVel

    integer :: ir, i
    real(dp), dimension(3) :: probeLocation, vel

    open(unit=10, file=ResultsDir//'probes'//timestamp//'.csv', &
      & status='replace',  action='write')
    write(10, 110) 'u', 'v', 'w', 'x', 'y', 'z'
    do i = 1, size(probe, 2)
      probeLocation = probe(:, i) + probeVel(:, i) * t
      vel = probeVel(:, i)
      do ir = 1, size(rotor)
        vel = vel &
          & + rotor(ir)%vind_bywing(probeLocation) &
          & + rotor(ir)%vind_bywake(probeLocation)
      enddo
      write(10, 120) vel, probeLocation
    enddo
    close(10)

    110 format(6(A15))
    120 format(6(E15.7))
  end subroutine probes2file

  subroutine filaments2file(timestamp, rotor)
    ! Write filaments to file for using with grid-based plots
    type(rotor_class), intent(in), dimension(:) :: rotor
    character(len=*), intent(in) :: timestamp

    integer :: nr, nvrWing, nvrNwake, nvfNwakeTE, nvfFwake
    integer :: ir, ib, irow, icol, indx
    type(vr_class), allocatable, dimension(:) :: vrWing, vrNwake
    type(vf_class), allocatable, dimension(:) :: vfFwake, vfNwakeTE
    real(dp), allocatable, dimension(:) :: gamFwake, gamNwakeTE

    nvrWing = 0
    nvrNwake = 0
    nvfNwakeTE = 0
    nvfFwake = 0

    nr = size(rotor)

    do ir = 1, nr
      if (rotor(ir)%rowFar .gt. rotor(ir)%nFwake) error stop 'ERROR: Use filaments2file() only after development of far wake'
    enddo

    ! Compute number of each filaments
    do ir = 1, nr
      nvrWing = nvrWing + rotor(ir)%nb*(rotor(ir)%nc*rotor(ir)%ns)
      nvrNwake = nvrNwake + rotor(ir)%nb*(rotor(ir)%nNwake*rotor(ir)%ns)
      nvfNwakeTE = nvfNwakeTE + rotor(ir)%nb*rotor(ir)%ns
      nvfFwake = nvfFwake + (rotor(ir)%nFwake - rotor(ir)%rowFar + 1)*rotor(ir)%nb
    enddo

    ! Allocate filaments
    allocate (vrWing(nvrWing))
    allocate (vrNwake(nvrNwake))
    allocate (vfNwakeTE(nvfNwakeTE))
    allocate (vfFwake(nvfFwake))
    allocate (gamFwake(nvfFwake))
    allocate (gamNwakeTE(nvfNwakeTE))

    ! Extract filament properties
    ! from wing
    indx = 1
    do ir = 1, nr
      do ib = 1, rotor(ir)%nb
        do icol = 1, rotor(ir)%ns
          do irow = 1, rotor(ir)%nc
            vrWing(indx) = rotor(ir)%blade(ib)%wiP(irow, icol)%vr
            indx = indx + 1
          enddo
        enddo
      enddo
    enddo

    ! from Nwake
    indx = 1
    do ir = 1, nr
      do ib = 1, rotor(ir)%nb
        do icol = 1, rotor(ir)%ns
          do irow = 1, rotor(ir)%nNwake
            vrNwake(indx) = rotor(ir)%blade(ib)%waP(irow, icol)%vr
            indx = indx + 1
          enddo
        enddo
      enddo
    enddo

    ! from NwakeTE
    indx = 1
    do ir = 1, nr
      irow = rotor(ir)%nNwake
      do ib = 1, rotor(ir)%nb
        do icol = 1, rotor(ir)%ns
          vfNwakeTE(indx) = rotor(ir)%blade(ib)%waP(irow, icol)%vr%vf(2)
          gamNwakeTE(indx) = rotor(ir)%blade(ib)%waP(irow, icol)%vr%gam*(-1._dp)
          indx = indx + 1
        enddo
      enddo
    enddo

    ! from Fwake
    indx = 1
    do ir = 1, nr
      do ib = 1, rotor(ir)%nb
        do irow = rotor(ir)%rowFar, rotor(ir)%nFwake
          vfFwake(indx) = rotor(ir)%blade(ib)%waF(irow)%vf
          gamFwake(indx) = rotor(ir)%blade(ib)%waF(irow)%gam
          indx = indx + 1
        enddo
      enddo
    enddo

    ! Write to filamentsXXXXX.dat binary file
    open (unit=10, file=ResultsDir//'filaments'//timestamp//'.dat', &
      & status='replace', action='write', form='unformatted')
    write (10) nvrWing
    write (10) nvrNwake
    write (10) nvfNwakeTE
    write (10) nvfFwake
    write (10) vrWing, vrNwake
    write (10) vfNwakeTE, gamNwakeTE
    write (10) vfFwake, gamFwake
    close (10)

    ! Deallocate filaments
    deallocate (vrWing)
    deallocate (vrNwake)
    deallocate (vfNwakeTE)
    deallocate (gamNwakeTE)
    deallocate (vfFwake)
    deallocate (gamFwake)
  end subroutine filaments2file

  subroutine mesh2file(wing_array, wake_array, filename)
    ! Obsolete
    type(wingpanel_class), intent(in), dimension(:, :) :: wing_array
    type(Nwake_class), intent(in), dimension(:, :) :: wake_array
    character(len=*), intent(in) :: filename
    character(len=5) :: nxChar, nyChar
    real(dp), dimension(3, size(wing_array, 1) + 1, size(wing_array, 2) + 1) :: wingMesh
    real(dp), dimension(3, size(wake_array, 1) + 1, size(wake_array, 2) + 1) :: wakeMesh
    integer :: i, j, nx, ny

    nx = size(wing_array, 1)
    ny = size(wing_array, 2)
    write (nxChar, '(I5)') nx + 1
    write (nyChar, '(I5)') ny + 1

    open (unit=10, file=filename, action='write', position='append')
    do j = 1, ny
      do i = 1, nx
        wingMesh(:, i, j) = wing_array(i, j)%pc(:, 1)
      enddo
    enddo
    do i = 1, nx
      wingMesh(:, i, ny + 1) = wing_array(i, ny)%pc(:, 4)
    enddo
    do j = 1, ny
      wingMesh(:, nx + 1, j) = wing_array(nx, j)%pc(:, 2)
    enddo
    wingMesh(:, nx + 1, ny + 1) = wing_array(nx, ny)%pc(:, 3)

    write (10, *) 'Title = "Panel array"'
    write (10, *) 'VARIABLES = "X" "Y" "Z" "GAM"'
    write (10, *) 'Zone I='//trim(nxChar)//' J='//trim(nyChar)//' K=1  T="Wing"'
    write (10, *) 'DATAPACKING=BLOCK'
    write (10, *) 'VARLOCATION=([4]=CELLCENTERED)'!,[5]=CELLCENTERED)'
    write (10, *) ((wingMesh(1, i, j), i=1, nx + 1), j=1, ny + 1)
    write (10, *) ((wingMesh(2, i, j), i=1, nx + 1), j=1, ny + 1)
    write (10, *) ((wingMesh(3, i, j), i=1, nx + 1), j=1, ny + 1)
    write (10, *) ((-1._dp*wing_array(i, j)%vr%gam, i=1, nx), j=1, ny)

    nx = size(wake_array, 1)
    ny = size(wake_array, 2)
    write (nxChar, '(I5)') nx + 1
    write (nyChar, '(I5)') ny + 1

    !Check if necessary - $omp parallel do collapse(2)
    do j = 1, ny
      do i = 1, nx
        wakeMesh(:, i, j) = wake_array(i, j)%vr%vf(1)%fc(:, 1)
      enddo
    enddo
    !Check if necessary -$omp end parallel do
    do i = 1, nx
      wakeMesh(:, i, ny + 1) = wake_array(i, ny)%vr%vf(4)%fc(:, 1)
    enddo
    do j = 1, ny
      wakeMesh(:, nx + 1, j) = wake_array(nx, j)%vr%vf(2)%fc(:, 1)
    enddo
    wakeMesh(:, nx + 1, ny + 1) = wake_array(nx, ny)%vr%vf(3)%fc(:, 1)

    write (10, *) 'Zone I='//trim(nxChar)//' J='//trim(nyChar)//' K=1  T="Wake"'
    write (10, *) 'DATAPACKING=BLOCK'
    write (10, *) 'VARLOCATION=([4]=CELLCENTERED)'!,[5]=CELLCENTERED)'
    write (10, *) ((wakeMesh(1, i, j), i=1, nx + 1), j=1, ny + 1)
    write (10, *) ((wakeMesh(2, i, j), i=1, nx + 1), j=1, ny + 1)
    write (10, *) ((wakeMesh(3, i, j), i=1, nx + 1), j=1, ny + 1)
    write (10, *) ((-1._dp*wake_array(i, j)%vr%gam, i=1, nx), j=1, ny)

    close (10)
  end subroutine mesh2file

  subroutine wingverify(wing_array)
    ! Write wing geometry to file in detail
    ! For verifying orientation of wing panels, bound vortex rings and CPs
    type(wingpanel_class), intent(in), dimension(:, :) :: wing_array
    character(len=5) :: nxChar, nyChar
    real(dp), dimension(3, size(wing_array, 1) + 1, size(wing_array, 2) + 1) :: wingMesh
    integer :: i, j, nx, ny

    nx = size(wing_array, 1)
    ny = size(wing_array, 2)
    write (nxChar, '(I5)') nx + 1
    write (nyChar, '(I5)') ny + 1

    do j = 1, ny
      do i = 1, nx
        wingMesh(:, i, j) = wing_array(i, j)%pc(:, 1)
      enddo
    enddo
    do i = 1, nx
      wingMesh(:, i, ny + 1) = wing_array(i, ny)%pc(:, 4)
    enddo
    do j = 1, ny
      wingMesh(:, nx + 1, j) = wing_array(nx, j)%pc(:, 2)
    enddo
    wingMesh(:, nx + 1, ny + 1) = wing_array(nx, ny)%pc(:, 3)

    open (unit=10, file=ResultsDir//'wingPC.plt', &
      & action='write', position='append')
    write (10, *) 'Title = "Panel Vertices"'
    write (10, *) 'VARIABLES = "X" "Y" "Z"'
    write (10, *) 'Zone I='//trim(nxChar)//' J='//trim(nyChar)//' K=1  T="Panel Vertices"'
    write (10, *) 'DATAPACKING=BLOCK'
    write (10, *) ((wingMesh(1, i, j), i=1, nx + 1), j=1, ny + 1)
    write (10, *) ((wingMesh(2, i, j), i=1, nx + 1), j=1, ny + 1)
    write (10, *) ((wingMesh(3, i, j), i=1, nx + 1), j=1, ny + 1)
    close (10)

    write (nxChar, '(I5)') nx
    write (nyChar, '(I5)') ny
    do j = 1, ny
      do i = 1, nx
        wingMesh(:, i, j) = wing_array(i, j)%CP
      enddo
    enddo
    open (unit=11, file=ResultsDir//'wingCP.plt', &
      & action='write', position='append')
    write (11, *) 'Title = "Coll. points"'
    write (11, *) 'VARIABLES = "X" "Y" "Z" "nx" "ny" "nz"' 
    write (11, *) 'Zone I='//trim(nxChar)//' J='//trim(nyChar)//' K=1  T="Coll. points"'
    write (11, *) 'DATAPACKING=BLOCK'
    write (11, *) ((wingMesh(1, i, j), i=1, nx), j=1, ny)
    write (11, *) ((wingMesh(2, i, j), i=1, nx), j=1, ny)
    write (11, *) ((wingMesh(3, i, j), i=1, nx), j=1, ny)
    write (11, *) ((wing_array(i, j)%nCap(1), i=1, nx), j=1, ny)
    write (11, *) ((wing_array(i, j)%nCap(2), i=1, nx), j=1, ny)
    write (11, *) ((wing_array(i, j)%nCap(3), i=1, nx), j=1, ny)
    close (11)

    write (nxChar, '(I5)') nx + 1
    write (nyChar, '(I5)') ny + 1
    do j = 1, ny
      do i = 1, nx
        wingMesh(:, i, j) = wing_array(i, j)%vr%vf(1)%fc(:, 1)
      enddo
    enddo
    do i = 1, nx
      wingMesh(:, i, ny + 1) = wing_array(i, ny)%vr%vf(4)%fc(:, 1)
    enddo
    do j = 1, ny
      wingMesh(:, nx + 1, j) = wing_array(nx, j)%vr%vf(2)%fc(:, 1)
    enddo
    wingMesh(:, nx + 1, ny + 1) = wing_array(nx, ny)%vr%vf(3)%fc(:, 1)

    open (unit=12, file=ResultsDir//'wingVR.plt', &
      & action='write', position='append')
    write (12, *) 'Title = "Vortex Rings"'
    write (12, *) 'VARIABLES = "X" "Y" "Z"'
    write (12, *) 'Zone I='//trim(nxChar)//' J='//trim(nyChar)//' K=1  T="Vortex Rings"'
    write (12, *) 'DATAPACKING=BLOCK'
    write (12, *) ((wingMesh(1, i, j), i=1, nx + 1), j=1, ny + 1)
    write (12, *) ((wingMesh(2, i, j), i=1, nx + 1), j=1, ny + 1)
    write (12, *) ((wingMesh(3, i, j), i=1, nx + 1), j=1, ny + 1)
    close (12)
  end subroutine wingverify

  subroutine tip2file(timestamp, rotor, rotorNumber)
    type(rotor_class), intent(in) :: rotor
    character(len=*), intent(in) :: timestamp
    integer, intent(in) :: rotorNumber
    character(len=5) :: nxChar, nyChar
    character(len=2) :: rotorNumberChar
    real(dp), dimension(3, rotor%nc + 1, rotor%ns + 1) :: wingMesh
    real(dp), dimension(3, rotor%nNwake + 1) :: nWakeTip
    real(dp), dimension(rotor%nNwake) :: gamRollup
    real(dp), dimension(3, rotor%nFwake + 1) :: fWakeTip
    real(dp) :: gamSum
    integer :: ib, i, j, nx, ny

    write (rotorNumberChar, '(I0.2)') rotorNumber

    open (unit=10, file=ResultsDir// &
      & 'r'//rotorNumberChar//'tip'//timestamp//'.plt', &
      & action='write', position='append')

    write (10, *) 'Title = "Wing and Tip"'
    write (10, *) 'VARIABLES = "X" "Y" "Z" "GAM" "ageAzimuthal"'! "Var6"'

    do ib = 1, rotor%nb
      ! Wing
      nx = rotor%nc
      ny = rotor%ns
      write (nxChar, '(I5)') nx + 1
      write (nyChar, '(I5)') ny + 1

      do j = 1, ny
        do i = 1, nx
          wingMesh(:, i, j) = rotor%blade(ib)%wiP(i, j)%pc(:, 1)
        enddo
      enddo
      do i = 1, nx
        wingMesh(:, i, ny + 1) = rotor%blade(ib)%wiP(i, ny)%pc(:, 4)
      enddo
      do j = 1, ny
        wingMesh(:, nx + 1, j) = rotor%blade(ib)%wiP(nx, j)%pc(:, 2)
      enddo
      wingMesh(:, nx + 1, ny + 1) = rotor%blade(ib)%wiP(nx, ny)%pc(:, 3)

      write (10, *) 'Zone I='//trim(nxChar)//' J='//trim(nyChar)//' K=1  T="Blade"'
      write (10, *) 'DATAPACKING=BLOCK'
      write (10, *) 'VARLOCATION=([4,5]=CELLCENTERED)'![6]=CELLCENTERED)'
      write (10, *) ((wingMesh(1, i, j), i=1, nx + 1), j=1, ny + 1)
      write (10, *) ((wingMesh(2, i, j), i=1, nx + 1), j=1, ny + 1)
      write (10, *) ((wingMesh(3, i, j), i=1, nx + 1), j=1, ny + 1)
      write (10, *) ((-1._dp*rotor%blade(ib)%wiP(i, j)%vr%gam, i=1, nx), j=1, ny)
      write(10,*) ((rotor%blade(ib)%wiP(i,j)% &
        & vr%vf(1)%ageAzimuthal*180._dp/pi, i=1,nx),j=1,ny)
      ! write(10,*) ((rotor%blade(ib)%wiP(i,j)%vr%skew,i=1,nx),j=1,ny)

      ! Near wake
      nx = rotor%nNwake
      write (nxChar, '(I5)') nx - (rotor%rowNear - 1) + 1

      ! Compute tip location for near wake
      do i = rotor%rowNear, nx
        gamRollup(i) = rotor%blade(ib)%waP(i, rotor%ns)%vr%gam
        gamSum = 0._dp
        nWakeTip(:, i) = 0._dp
        do j = rotor%rollupStart, rotor%rollupEnd
          nWakeTip(:, i) = nWakeTip(:, i) + rotor%blade(ib)%waP(i, j)%vr%vf(4)%fc(:, 1)* &
            rotor%blade(ib)%waP(i, j)%vr%gam
          gamSum = gamSum + rotor%blade(ib)%waP(i, j)%vr%gam

          ! Compute max gam for near wake filaments
          if (sign(1._dp, rotor%Omega*rotor%controlPitch(1)) > eps) then
            if (rotor%blade(ib)%waP(i, j)%vr%gam < gamRollup(i)) then    ! '<' because of negative gamma
              gamRollup(i) = rotor%blade(ib)%waP(i, j)%vr%gam
            endif
          else    ! one of Omega or pitch is negative
            if (rotor%blade(ib)%waP(i, j)%vr%gam > gamRollup(i)) then    ! '>' because of positive gamma
              gamRollup(i) = rotor%blade(ib)%waP(i, j)%vr%gam
            endif
          endif
        enddo

        if (abs(gamSum) > eps) then
          nWakeTip(:, i) = nWakeTip(:, i)/gamSum
        else
          nWakeTip(:, i) = rotor%blade(ib)%waP(i, rotor%rollupEnd)%vr%vf(4)%fc(:, 1)
        endif

      enddo

      ! For last row
      gamSum = 0._dp
      nWakeTip(:, nx + 1) = 0._dp
      do j = rotor%rollupStart, rotor%rollupEnd
        nWakeTip(:, nx + 1) = nWakeTip(:, nx + 1) + rotor%blade(ib)%waP(nx, j)%vr%vf(3)%fc(:, 1)* &
          rotor%blade(ib)%waP(nx, j)%vr%gam
        gamSum = gamSum + rotor%blade(ib)%waP(nx, j)%vr%gam
      enddo

      if (abs(gamSum) > eps) then
        nWakeTip(:, nx + 1) = nWakeTip(:, nx + 1)/gamSum
      else
        nWakeTip(:, nx + 1) = rotor%blade(ib)%waP(nx, rotor%rollupEnd)%vr%vf(3)%fc(:, 1)
      endif

      write (10, *) 'Zone I='//trim(nxChar)//' J=1    K=1  T="NearWake"'
      write (10, *) 'DATAPACKING=BLOCK'
      write (10, *) 'VARLOCATION=([4,5]=CELLCENTERED)'
      write (10, *) (nWakeTip(1, i), i=rotor%rowNear, nx + 1)
      write (10, *) (nWakeTip(2, i), i=rotor%rowNear, nx + 1)
      write (10, *) (nWakeTip(3, i), i=rotor%rowNear, nx + 1)
      write (10, *) (-1._dp*gamRollup(i), i=rotor%rowNear, nx)
      write(10,*) (rotor%blade(ib)%waP(i, rotor%ns)% &
        & vr%vf(4)%ageAzimuthal*180._dp/pi, i=rotor%rowNear, nx)
      !write(10,*) ((rotor%blade(ib)%waP(i,j)%vr%skew,i=rotor%rowNear,nx),j=1,ny)

      ! Far wake
      nx = rotor%nFwake
      if (rotor%rowFar .le. rotor%nFwake) then
        write (nxChar, '(I5)') nx - (rotor%rowFar - 1) + 1

        !Check if necessary - $omp parallel do collapse(2)
        do i = rotor%rowFar, nx
          fWakeTip(:, i) = rotor%blade(ib)%waF(i)%vf%fc(:, 2)
        enddo
        fWakeTip(:, nx + 1) = rotor%blade(ib)%waF(rotor%nFwake)%vf%fc(:, 1)
        !Check if necessary -$omp end parallel do

        write (10, *) 'Zone I='//trim(nxChar)//' J=1   K=1   T="FarWake"'
        write (10, *) 'DATAPACKING=BLOCK'
        write (10, *) 'VARLOCATION=([4,5]=CELLCENTERED)'!,[6]=CELLCENTERED)'
        write (10, *) (fWakeTip(1, i), i=rotor%rowFar, nx + 1)
        write (10, *) (fWakeTip(2, i), i=rotor%rowFar, nx + 1)
        write (10, *) (fWakeTip(3, i), i=rotor%rowFar, nx + 1)
        write (10, *) (-1._dp*rotor%blade(ib)%waF(i)%gam, i=rotor%rowFar, nx)
        write(10,*) (rotor%blade(ib)%waF(i)%vf%ageAzimuthal*180._dp/pi, &
          & i=rotor%rowFar, nx)
        !write(10,*) (rotor%blade(ib)%waF(i)%vf%age,i=rotor%rowFar,nx)

      else  ! No far wake present

        write (nxChar, '(I5)') 2  ! Plot mesh as single redundant point
        fWakeTip(:, 1) = rotor%blade(ib)%waP(rotor%nNwake, rotor%ns)%vr%vf(3)%fc(:, 1)

        write (10, *) 'Zone I='//trim(nxChar)//' J=1   K=1   T="FarWake"'
        write (10, *) 'DATAPACKING=BLOCK'
        write (10, *) 'VARLOCATION=([4,5]=CELLCENTERED)'!,[6]=CELLCENTERED)'
        write (10, *) fWakeTip(1, 1), fWakeTip(1, 1)
        write (10, *) fWakeTip(2, 1), fWakeTip(2, 1)
        write (10, *) fWakeTip(3, 1), fWakeTip(3, 1)
        write (10, *) 0._dp
        write (10, *) 0._dp
        !write(10,*) 0._dp
      endif
    enddo

    close (10)

  end subroutine tip2file

  subroutine skew2file(timestamp, rotor, rotorNumber)
    type(rotor_class), intent(in) :: rotor
    character(len=*), intent(in) :: timestamp
    integer, intent(in) :: rotorNumber
    character(len=2) :: rotorNumberChar, bladeNumberChar
    integer :: ib, irow, nrow, ncol

    write (rotorNumberChar, '(I0.2)') rotorNumber

    nrow = size(rotor%blade(1)%waP, 1)
    ncol = size(rotor%blade(1)%waP, 2)

    do ib = 1, rotor%nb
      write (bladeNumberChar, '(I0.2)') ib
      open (unit=12, file=ResultsDir// &
        & 'r'//rotorNumberChar// 'b'//bladeNumberChar// &
        & 'skew'//timestamp//'.csv', & 
        & action='write')

      write(12, 110) 'max', 'avg'
      do irow = nrow, rotor%rowNear, -1
        write(12, 120) maxval(rotor%blade(ib)%waP(irow, :)%vr%skew), &
          & sum(rotor%blade(ib)%waP(irow, :)%vr%skew)/ncol
      enddo

      110 format(2(A15))
      120 format(2(E15.7))
      close(12)
    enddo
  end subroutine skew2file

  subroutine force2file(timestamp, rotor, rotorNumber)
    ! Write sec and net force to file
    type(rotor_class), intent(in) :: rotor
    character(len=*), intent(in) :: timestamp
    integer, intent(in) :: rotorNumber
    character(len=2) :: rotorNumberChar, bladeNumberChar
    integer :: ib, ispan, iter
    character(len=24) :: forceDimFilename
    character(len=27) :: forceNonDimFilename

    write (rotorNumberChar, '(I0.2)') rotorNumber

    forceNonDimFilename = ResultsDir//'r'//rotorNumberChar//'ForceNonDim.csv'
    open (unit=11, file=forceNonDimFilename, action='write', position='append')
    write (11, 100) timestamp, &
      norm2(rotor%lift) / rotor%nonDimforceDenominator, &          ! CL
      norm2(rotor%drag) / rotor%nonDimforceDenominator, &          ! CD
      norm2(rotor%dragInduced) / rotor%nonDimforceDenominator, &   ! CDi
      norm2(rotor%dragProfile) / rotor%nonDimforceDenominator, &   ! CDo
      norm2(rotor%dragUnsteady) / rotor%nonDimforceDenominator, &  ! CDu
      rotor%forceInertial(1) / rotor%nonDimforceDenominator, &     ! CFx
      rotor%forceInertial(2) / rotor%nonDimforceDenominator, &     ! CFy
      rotor%forceInertial(3) / rotor%nonDimforceDenominator        ! CFz
    close (11)
    100 format(A, 8(E15.7))

    forceDimFilename = ResultsDir//'r'//rotorNumberChar//'ForceDim.csv'
    open (unit=12, file=forceDimFilename, action='write', position='append')
    write (12, 101) timestamp, norm2(rotor%lift), norm2(rotor%drag), &
      rotor%lift(1), rotor%lift(2), rotor%lift(3), &     ! Lift
      rotor%drag(1), rotor%drag(2), rotor%drag(3), &     ! Drag
      rotor%forceInertial(1), rotor%forceInertial(2), rotor%forceInertial(3)  ! forceInertial 
    close (12)
    101 format(A, 11(E15.7))

    if (rotor%bladeforcePlotSwitch .ne. 0) then
      read(timestamp, *) iter
      if (mod(iter, rotor%bladeforcePlotSwitch) .eq. 0) then
        do ib = 1, rotor%nb
          write (bladeNumberChar, '(I0.2)') ib
          open (unit=12, file=ResultsDir// &
            & 'r'//rotorNumberChar// 'b'//bladeNumberChar// &
            & 'ForceDist'//timestamp//'.csv', & 
            & action='write', position='append')
          write (12, 202) 'secSpan', 'secCL', 'secCD', 'secLift', 'secDrag', &
            & 'secArea', 'secVel', 'secChord', 'secAlpha'
          do ispan = 1, rotor%ns
            write (12, 102) dot_product(rotor%blade(ib)%secCP(:, ispan) - &
              & rotor%hubCoords, rotor%blade(ib)%yAxis), &
              & rotor%blade(ib)%secCL(ispan), &
              & rotor%blade(ib)%secCD(ispan), &
              & norm2(rotor%blade(ib)%secLift(:, ispan)), &
              & norm2(rotor%blade(ib)%secDrag(:, ispan)), &
              & rotor%blade(ib)%secArea(ispan), &
              & norm2(rotor%blade(ib)%secChordwiseResVel(:, ispan)), &
              & rotor%blade(ib)%secChord(ispan), &
              & rotor%blade(ib)%secAlpha(ispan)*180._dp/pi
          enddo
          close (12)
        enddo
        202 format(9(A15))
        102 format(9(E15.7))
      endif
    endif
  end subroutine force2file

  subroutine inflow2file(timestamp, rotorArray, rotorNumber, directionVector)
    ! Calculate inflow velocity along directionVector on the blades of rotor(rotorNumber)
    ! at rotor(rotorNumber)%secCP
    character(len=*), intent(in) :: timestamp
    type(rotor_class), intent(inout), dimension(:) :: rotorArray
    integer, intent(in) :: rotorNumber
    real(dp), intent(in), dimension(3) :: directionVector
    character(len=2) :: rotorNumberChar, bladeNumberChar
    integer :: is, ir, ib

    real(dp), dimension(3) :: P
    real(dp), dimension(rotorArray(rotorNumber)%ns, rotorArray(rotorNumber)%nb) :: inflowVel

    inflowVel = 0._dp
    do ir = 1, size(rotorArray)
      do ib = 1, rotorArray(rotorNumber)%nb
        do is = 1, rotorArray(rotorNumber)%ns
          P = rotorArray(rotorNumber)%blade(ib)%secCP(:, is)
          inflowVel(is, ib) = inflowVel(is, ib) + dot_product(rotorArray(ir)%vind_bywing(P), directionVector)
          inflowVel(is, ib) = inflowVel(is, ib) - dot_product(rotorArray(ir)%vind_bywing_boundVortices(P), directionVector)
          inflowVel(is, ib) = inflowVel(is, ib) + dot_product(rotorArray(ir)%vind_bywake(P), directionVector)
        enddo
      enddo
    enddo

    ! Write to file
    write (rotorNumberChar, '(I0.2)') rotorNumber
    do ib = 1, rotorArray(rotorNumber)%nb
      write (bladeNumberChar, '(I0.2)') ib
      open (unit=12, file=ResultsDir// &
        & 'r'//rotorNumberChar//'b'//bladeNumberChar// &
        & 'inflowDist'//timestamp//'.curve', & 
        & action='write', position='append')
      write(12, 100) 'secSpan', 'inflowVel'
      do is = 1, rotorArray(rotorNumber)%ns
        write (12, 101) dot_product(rotorArray(rotorNumber)%blade(ib)%secCP(:, is), &
          & rotorArray(rotorNumber)%hubCoords - &
          & rotorArray(rotorNumber)%blade(ib)%yAxis), &
          & inflowVel(is, ib)
      enddo
      close (12)
    enddo
    100 format (2(A15))
    101 format (2(F15.7))

  end subroutine inflow2file

  subroutine gamma2file(timestamp, rotor, rotorNumber)
    ! Calculate inflow velocity along directionVector on the blades of rotor(rotorNumber)
    ! at rotor(rotorNumber)%secCP
    character(len=*), intent(in) :: timestamp
    type(rotor_class), intent(in) :: rotor
    integer, intent(in) :: rotorNumber
    integer :: ib, ic, is
    character(len=2) :: rotorNumberChar, bladeNumberChar, rowNumberChar

    ! Write to file
    write (rotorNumberChar, '(I0.2)') rotorNumber
    do ib = 1, rotor%nb
      write (bladeNumberChar, '(I0.2)') ib
      open (unit=12, &
        & file=ResultsDir// &
        & 'r'//rotorNumberChar//'b'//bladeNumberChar// &
        & 'gammaDist'//timestamp//'.curve', & 
        & action='write', position='append')
      ic = 1
      write (rowNumberChar, '(I0.2)') ic
      write (12, *) '# Row'//rowNumberChar
      do is = 1, rotor%ns
        write (12, *) dot_product(rotor%blade(ib)%wiP(1, is)%CP - rotor%hubCoords, rotor%blade(ib)%yAxis), &
          -1._dp*sum(rotor%blade(ib)%wiP(:, is)%vr%gam)
      enddo
      close (12)
    enddo

  end subroutine gamma2file

  ! subroutine alpha2file(timestamp, rotor, rotorNumber)
  !   ! Write sec alpha to file
  !   character(len=*), intent(in) :: timestamp
  !   type(rotor_class), intent(inout) :: rotor
  !   integer, intent(in) :: rotorNumber
  !   integer :: ib, is
  !   character(len=2) :: rotorNumberChar, bladeNumberChar

  !   ! Write to file
  !   write (rotorNumberChar, '(I0.2)') rotorNumber
  !   do ib = 1, rotor%nb
  !     write (bladeNumberChar, '(I0.2)') ib
  !     open (unit=12, & 
  !       & file=ResultsDir// &
  !       & 'r'//rotorNumberChar//'b'//bladeNumberChar// &
  !       & 'alphaDist'//timestamp//'.curve', &
  !       & action='write', position='append')
  !     write(12, 100) 'secSpan', 'alpha'
  !     do is = 1, rotor%ns
  !       write (12, 101) dot_product(rotor%blade(ib)%secCP(:, is) - rotor%hubCoords, &
  !         & rotor%blade(ib)%yAxis), rotor%blade(ib)%secAlpha(is)*180._dp/pi
  !     enddo
  !     close (12)
  !   enddo
  !   100 format (2(A15))
  !   101 format (2(F15.7))

  ! end subroutine alpha2file

end module libPostprocess

