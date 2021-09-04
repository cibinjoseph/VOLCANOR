module libPostprocess
  use libMath, only: dp, pi, eps, radToDeg
  character(len=8) :: ResultsDir = 'Results/'

contains

  subroutine init_plots(numOfRotors)
    ! Initialise headers for plot files
    integer, intent(in) :: numOfRotors
    character(len=30) :: forceDimFilename
    character(len=30) :: forceNonDimFilename
    character(len=30) :: dynamicsFilename
    integer :: rotorNumber
    character(len=2) :: rotorNumberChar

    do rotorNumber = 1, numOfRotors
      write (rotorNumberChar, '(I0.2)') rotorNumber
      forceDimFilename = ResultsDir//'r'//rotorNumberChar//'ForceDim.csv'
      forceNonDimFilename = ResultsDir//'r'//rotorNumberChar//'ForceNonDim.csv'
      dynamicsFilename = ResultsDir//'r'//rotorNumberChar//'bladedynamics.csv'

      ! Add data headers
      open (unit=10, file=forceDimFilename, &
        & status='replace', action='write')
      write (10, 100) 'iter','LiftMag','DragMag', &
        & 'Lx','Ly','Lz', &
        & 'Dx','Dy','Dz', &
        & 'FInertx','FInerty','FInertz'
      close (10)

      open (unit=11, file=forceNonDimFilename, &
        & status='replace', action='write')
      write (11, 101) 'iter','CL/CT','CD/CQ', 'CLu','CDi','CD0','CDu', &
        & 'CFx','CFy','CFz'
      close (11)

      open(unit=12, file=dynamicsFilename, &
        & status='replace', action='write')
      write(12, 102) 'iter', 'flap', 'dflap'
      close(12)
    enddo

    100 format (A5, 11(A15))
    101 format (A5, 9(A15))
    102 format (A5, 2(A15))
  end subroutine init_plots

  subroutine params2file(rotor, nt, dt, nr, &
      & density, velSound, switches)
    ! Write rotor parameters to json file
    use classdef, only: rotor_class, switches_class
    type(rotor_class), intent(in) :: rotor
    type(switches_class), intent(in) :: switches
    integer, intent(in) :: nt, nr
    real(dp), intent(in) :: dt, density, velSound
    character(len=22) :: paramsFilename

    paramsFilename = ResultsDir//'r'//rotor%id//'Params.json'
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
    write(10, *) '"coningAngle": ', rotor%flapInitial*radToDeg, ','
    write(10, *) '"Omega": ', rotor%Omega, ','
    write(10, *) '"phi": ', rotor%pts(1), ','
    write(10, *) '"theta": ', rotor%pts(2), ','
    write(10, *) '"psi": ', rotor%pts(3), ','
    write(10, *) '"theta0": ', rotor%controlPitch(1)*radToDeg, ','
    write(10, *) '"thetaC": ', rotor%controlPitch(2)*radToDeg, ','
    write(10, *) '"thetaS": ', rotor%controlPitch(3)*radToDeg, ','
    write(10, *) '"thetaTwist": ', rotor%thetaTwist*radToDeg, ','
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
    write(10, *) '"prescWakeNt": ', rotor%prescWakeNt, ','
    write(10, *) '"initWakeVel": ', rotor%initWakeVel, ','
    write(10, *) '"psiStart": ', rotor%psiStart*radToDeg, ','
    write(10, *) '"forceCalcSwitch": ', rotor%forceCalcSwitch*radToDeg, ','
    write(10, *) '"apparentViscCoeff": ', rotor%apparentViscCoeff, ','
    write(10, *) '"decayCoeff": ', rotor%decayCoeff
    write(10, *) '}'
    close(10)
  end subroutine params2file

  subroutine geom2file(timestamp, rotor)
    ! Plot rotor geometry and wake to file
    use classdef, only: rotor_class
    type(rotor_class), intent(in) :: rotor
    character(len=*), intent(in) :: timestamp
    character(len=5) :: nxChar, nyChar
    real(dp), dimension(3, rotor%nc + 1, rotor%ns + 1) :: wingMesh
    real(dp), dimension(3, rotor%nNwake + 1, rotor%ns + 1) :: wakeMesh
    real(dp), dimension(3, rotor%nFwake + 1) :: wakeTip   ! Optimise this by only initialising reqd size
    real(dp), dimension(3, size(rotor%blade(1)%wapF%waF) + 1) :: wakeTipPresc
    integer :: i, j, nx, ny, ib

    open (unit=10, file=ResultsDir// &
      & 'r'//rotor%id//'wingNwake'//timestamp//'.plt', position='append')

    if (abs(rotor%surfaceType) == 1) then
      ! Lifting surface
      write (10, *) 'Title = "Wing and Wake"'
      write (10, *) 'VARIABLES = "X" "Y" "Z" "GAM" "skew"' ! "Var6"'

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

        write (10, *) 'Zone I='//trim(nxChar)//' J='//trim(nyChar)// &
          & ' K=1  T="Blade '//rotor%blade(ib)%id//'"'
        write (10, *) 'DATAPACKING=BLOCK'
        write (10, *) 'VARLOCATION=([4]=CELLCENTERED,[5]=CELLCENTERED)' !,[6]=CELLCENTERED)'
        write (10, *) ((wingMesh(1, i, j), i=1, nx + 1), j=1, ny + 1)
        write (10, *) ((wingMesh(2, i, j), i=1, nx + 1), j=1, ny + 1)
        write (10, *) ((wingMesh(3, i, j), i=1, nx + 1), j=1, ny + 1)
        write (10, *) ((-1._dp*rotor%blade(ib)%wiP(i, j)%vr%gam, i=1, nx), j=1, ny)
        write(10, *) ((rotor%blade(ib)%wiP(i, j)%vr%gam*0._dp, i=1,nx), j=1,ny)
        !write(10,*) ((rotor%blade(ib)%wiP(i,j)%vr%skew,i=1,nx),j=1,ny)

        ! Near wake
        nx = rotor%nNwakeEnd
        ny = rotor%ns
        write (nxChar, '(I5)') nx - (rotor%rowNear - 1) + 1
        write (nyChar, '(I5)') ny + 1

        !Check if necessary - $omp parallel do collapse(2)
        do j = 1, ny
          do i = rotor%rowNear, nx
            wakeMesh(:, i, j) = rotor%blade(ib)%waN(i, j)%vr%vf(1)%fc(:, 1)
          enddo
        enddo
        !Check if necessary -$omp end parallel do
        do i = rotor%rowNear, nx
          wakeMesh(:, i, ny + 1) = rotor%blade(ib)%waN(i, ny)%vr%vf(4)%fc(:, 1)
        enddo
        do j = 1, ny
          wakeMesh(:, nx + 1, j) = rotor%blade(ib)%waN(nx, j)%vr%vf(2)%fc(:, 1)
        enddo
        wakeMesh(:, nx + 1, ny + 1) = rotor%blade(ib)%waN(nx, ny)%vr%vf(3)%fc(:, 1)

        write (10, *) 'Zone I='//trim(nxChar)//' J='//trim(nyChar)//' K=1  T="NearWake"'
        write (10, *) 'DATAPACKING=BLOCK'
        write (10, *) 'VARLOCATION=([4]=CELLCENTERED, [5]=CELLCENTERED)' !,[6]=CELLCENTERED)'
        write (10, *) ((wakeMesh(1, i, j), i=rotor%rowNear, nx + 1), j=1, ny + 1)
        write (10, *) ((wakeMesh(2, i, j), i=rotor%rowNear, nx + 1), j=1, ny + 1)
        write (10, *) ((wakeMesh(3, i, j), i=rotor%rowNear, nx + 1), j=1, ny + 1)
        write (10, *) ((-1._dp*rotor%blade(ib)%waN(i, j)%vr%gam, i=rotor%rowNear, nx), j=1, ny)
        write(10,*) ((rotor%blade(ib)%waN(i,j)%vr%skew, i=rotor%rowNear, nx), j=1, ny)
        !write(10,*) ((rotor%blade(ib)%waN(i,j)%vr%skew,i=rotor%rowNear,nx),j=1,ny)

        ! Far wake
        nx = rotor%nFwakeEnd
        if (rotor%rowFar .le. rotor%nFwakeEnd) then
          write (nxChar, '(I5)') nx - (rotor%rowFar - 1) + 1

          !Check if necessary - $omp parallel do collapse(2)
          do i = rotor%rowFar, nx
            wakeTip(:, i) = rotor%blade(ib)%waF(i)%vf%fc(:, 2)
          enddo
          wakeTip(:, nx + 1) = rotor%blade(ib)%waF(rotor%nFwakeEnd)%vf%fc(:, 1)
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
          wakeTip(:, 1) = rotor%blade(ib)%waN(rotor%nNwake, rotor%ns)%vr%vf(3)%fc(:, 1)

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

        ! Prescribed far wake
        nx = size(rotor%blade(ib)%wapF%waF)
        if (rotor%blade(ib)%wapF%isPresent) then
          write (nxChar, '(I5)') nx + 1

          !Check if necessary - $omp parallel do collapse(2)
          do i = 1, nx
            wakeTipPresc(:, i) = rotor%blade(ib)%wapF%waF(i)%vf%fc(:, 2)
          enddo
          wakeTipPresc(:, nx + 1) = rotor%blade(ib)%wapF%waF(size(rotor%blade(1)%wapF%waF))%vf%fc(:, 1)
          !Check if necessary -$omp end parallel do

          write (10, *) 'Zone I='//trim(nxChar)//' J=1   K=1   T="PrescFarWake"'
          write (10, *) 'DATAPACKING=BLOCK'
          write (10, *) 'VARLOCATION=([4]=CELLCENTERED ,[5]=CELLCENTERED)' !,[6]=CELLCENTERED)'
          write (10, *) (wakeTipPresc(1, i), i=1, nx + 1)
          write (10, *) (wakeTipPresc(2, i), i=1, nx + 1)
          write (10, *) (wakeTipPresc(3, i), i=1, nx + 1)
          write (10, *) (-1._dp*rotor%blade(ib)%wapF%waF(i)%gam, i=1, nx)
          write(10,*) (rotor%blade(ib)%wapF%waF(i)%vf%rVc, i=1, nx)
          !write(10,*) (rotor%blade(ib)%waF(i)%vf%age,i=rotor%rowFar,nx)

        else  ! No prescribed far wake present

          write (nxChar, '(I5)') 2  ! Plot mesh as single redundant point
          wakeTip(:, 1) = rotor%blade(ib)%waN(rotor%nNwake, rotor%ns)%vr%vf(3)%fc(:, 1)

          write (10, *) 'Zone I='//trim(nxChar)//' J=1   K=1   T="PrescFarWake"'
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
        ! Compute common nodes of non-lifting surface
        nx = rotor%blade(ib)%stlNodesCols
        ny = rotor%nc ! No. of triangular elements
        write (nxChar, '(I0.5)') nx
        write (nyChar, '(I0.5)') ny

        write(10, *) 'Zone NODES='//trim(nxChar)// &
          & ' ELEMENTS='//nyChar// &
          & ' T="Blade '//rotor%blade(ib)%id//'"'
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
  end subroutine geom2file

  subroutine probes2file(timestamp, probe, probeVel, rotor, t)
    ! Write velocities at probe locations
    use classdef, only: rotor_class
    type(rotor_class), intent(inout), dimension(:) :: rotor
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
    use classdef, only: vf_class, vr_class, rotor_class
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
            vrNwake(indx) = rotor(ir)%blade(ib)%waN(irow, icol)%vr
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
          vfNwakeTE(indx) = rotor(ir)%blade(ib)%waN(irow, icol)%vr%vf(2)
          gamNwakeTE(indx) = rotor(ir)%blade(ib)%waN(irow, icol)%vr%gam*(-1._dp)
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
    use classdef, only: wingpanel_class, Nwake_class
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

  subroutine geomSurface2file(rotor)
    ! Plot surface geometry to file
    use classdef, only: rotor_class
    type(rotor_class), intent(in) :: rotor

    character(len=5) :: nxChar, nyChar
    real(dp), dimension(3, rotor%nc+1, rotor%ns+1) :: mesh
    integer :: i, j, nx, ny, ib

    do ib = 1, rotor%nb
      open (unit=10, file=ResultsDir// &
        & 'r'//rotor%id//'b'//rotor%blade(ib)%id//'Surface.plt', &
        & action='write', position='append')

      write (10, *) 'Title="r'//rotor%id//'b'//rotor%blade(ib)%id//'"'
      write (10, *) 'VARIABLES = "X" "Y" "Z" "nx" "ny" "nz"' 
      ! nx , ny, nz can be used for vectors at nodes

      ! Zone 1: Panel coordinates
      nx = rotor%nc
      ny = rotor%ns
      write (nxChar, '(I5)') nx + 1
      write (nyChar, '(I5)') ny + 1

      do j = 1, ny
        do i = 1, nx
          mesh(:, i, j) = rotor%blade(ib)%wiP(i, j)%PC(:, 1)
        enddo
      enddo
      do i = 1, nx
        mesh(:, i, ny + 1) = rotor%blade(ib)%wiP(i, ny)%PC(:, 4)
      enddo
      do j = 1, ny
        mesh(:, nx + 1, j) = rotor%blade(ib)%wiP(nx, j)%PC(:, 2)
      enddo
      mesh(:, nx + 1, ny + 1) = rotor%blade(ib)%wiP(nx, ny)%PC(:, 3)

      write (10, *) 'Zone I='//trim(nxChar)//' J='//trim(nyChar)//' K=1 T="PC"'
      write (10, *) 'DATAPACKING=BLOCK'
      write (10, *) ((mesh(1, i, j), i=1, nx + 1), j=1, ny+1)
      write (10, *) ((mesh(2, i, j), i=1, nx + 1), j=1, ny+1)
      write (10, *) ((mesh(3, i, j), i=1, nx + 1), j=1, ny+1)
      write (10, *) ((0, i=1, nx + 1), j=1, ny+1)
      write (10, *) ((0, i=1, nx + 1), j=1, ny+1)
      write (10, *) ((0, i=1, nx + 1), j=1, ny+1)

      ! Zone 2: Vortex filament coordinates
      write (nxChar, '(I5)') nx + 1
      write (nyChar, '(I5)') ny + 1
      do j = 1, ny
        do i = 1, nx
          mesh(:, i, j) = rotor%blade(ib)%wiP(i, j)%vr%vf(1)%fc(:, 1)
        enddo
      enddo
      do i = 1, nx
        mesh(:, i, ny + 1) = rotor%blade(ib)%wiP(i, ny)%vr%vf(4)%fc(:, 1)
      enddo
      do j = 1, ny
        mesh(:, nx + 1, j) = rotor%blade(ib)%wiP(nx, j)%vr%vf(2)%fc(:, 1)
      enddo
      mesh(:, nx + 1, ny + 1) = rotor%blade(ib)%wiP(nx, ny)%vr%vf(3)%fc(:, 1)

      write (10, *) 'Zone I='//trim(nxChar)//' J='//trim(nyChar)//' K=1 T="VR"'
      write (10, *) 'DATAPACKING=BLOCK'
      write (10, *) 'VARSHARELIST=([4-6]=1)'  ! Empty variables
      write (10, *) ((mesh(1, i, j), i=1, nx+1), j=1, ny+1)
      write (10, *) ((mesh(2, i, j), i=1, nx+1), j=1, ny+1)
      write (10, *) ((mesh(3, i, j), i=1, nx+1), j=1, ny+1)

      ! Zone 3: CP and nCap
      nx = rotor%nc
      ny = rotor%ns
      write (nxChar, '(I5)') nx
      write (nyChar, '(I5)') ny

      write (10, *) 'Zone I='//trim(nxChar)//' J='//trim(nyChar)//' K=1 &
        & T="CP, nCap"'
      write (10, *) 'DATAPACKING=BLOCK'
      write (10, *) ((rotor%blade(ib)%wiP(i, j)%CP(1), i=1, nx), j=1, ny)
      write (10, *) ((rotor%blade(ib)%wiP(i, j)%CP(2), i=1, nx), j=1, ny)
      write (10, *) ((rotor%blade(ib)%wiP(i, j)%CP(3), i=1, nx), j=1, ny)
      write (10, *) ((rotor%blade(ib)%wiP(i, j)%nCap(1), i=1, nx), j=1, ny)
      write (10, *) ((rotor%blade(ib)%wiP(i, j)%nCap(2), i=1, nx), j=1, ny)
      write (10, *) ((rotor%blade(ib)%wiP(i, j)%nCap(3), i=1, nx), j=1, ny)

      ! Zone4: CP and tauCapChord
      write (10, *) 'Zone I='//trim(nxChar)//' J='//trim(nyChar)//' K=1 &
        & T="CP, tauCapChord"'
      write (10, *) 'DATAPACKING=BLOCK'
      write (10, *) 'VARSHARELIST=([1-3]=3)'  ! Share CP coordinates
      write (10, *) ((rotor%blade(ib)%wiP(i, j)%tauCapChord(1), i=1, nx), j=1, ny)
      write (10, *) ((rotor%blade(ib)%wiP(i, j)%tauCapChord(2), i=1, nx), j=1, ny)
      write (10, *) ((rotor%blade(ib)%wiP(i, j)%tauCapChord(3), i=1, nx), j=1, ny)

      ! Zone 5: CP and tauCapSpan
      write (10, *) 'Zone I='//trim(nxChar)//' J='//trim(nyChar)//' K=1 &
        & T="CP, tauCapChord"'
      write (10, *) 'DATAPACKING=BLOCK'
      write (10, *) 'VARSHARELIST=([1-3]=3)'  ! Share CP coordinates
      write (10, *) ((rotor%blade(ib)%wiP(i, j)%tauCapSpan(1), i=1, nx), j=1, ny)
      write (10, *) ((rotor%blade(ib)%wiP(i, j)%tauCapSpan(2), i=1, nx), j=1, ny)
      write (10, *) ((rotor%blade(ib)%wiP(i, j)%tauCapSpan(3), i=1, nx), j=1, ny)

      ! Zone 6: Sectional CP and secNormalVec
      nx = rotor%ns
      write (nxChar, '(I5)') nx
      write (10, *) 'Zone I='//trim(nxChar)//' J=1 K=1 T="secCP, secNcap"'
      write (10, *) 'DATAPACKING=BLOCK'
      write (10, *) (rotor%blade(ib)%secCP(1, i), i=1, nx)
      write (10, *) (rotor%blade(ib)%secCP(2, i), i=1, nx)
      write (10, *) (rotor%blade(ib)%secCP(3, i), i=1, nx)
      write (10, *) (rotor%blade(ib)%secNormalVec(1, i), i=1, nx)
      write (10, *) (rotor%blade(ib)%secNormalVec(2, i), i=1, nx)
      write (10, *) (rotor%blade(ib)%secNormalVec(3, i), i=1, nx)

      ! Zone 7: Sectional CP and secNormalVec
      nx = rotor%ns
      write (nxChar, '(I5)') nx
      write (10, *) 'Zone I='//trim(nxChar)//' J=1 K=1 &
        & T="secCP, secTauCapChord"'
      write (10, *) 'DATAPACKING=BLOCK'
      write (10, *) 'VARSHARELIST=([1-3]=6)'  ! Share secCP coordinates
      write (10, *) (rotor%blade(ib)%secTauCapChord(1, i), i=1, nx)
      write (10, *) (rotor%blade(ib)%secTauCapChord(2, i), i=1, nx)
      write (10, *) (rotor%blade(ib)%secTauCapChord(3, i), i=1, nx)

      ! Zone 8: Sectional CP and secNormalVec
      nx = rotor%ns
      write (nxChar, '(I5)') nx
      write (10, *) 'Zone I='//trim(nxChar)//' J=1 K=1 &
        & T="secCP, secTauCapSpan"'
      write (10, *) 'DATAPACKING=BLOCK'
      write (10, *) 'VARSHARELIST=([1-3]=6)'  ! Share secCP coordinates
      write (10, *) (rotor%blade(ib)%secTauCapSpan(1, i), i=1, nx)
      write (10, *) (rotor%blade(ib)%secTauCapSpan(2, i), i=1, nx)
      write (10, *) (rotor%blade(ib)%secTauCapSpan(3, i), i=1, nx)

      close(10)
    enddo

  end subroutine geomSurface2file

  subroutine tip2file(timestamp, rotor)
    use classdef, only: rotor_class
    type(rotor_class), intent(in) :: rotor
    character(len=*), intent(in) :: timestamp
    character(len=5) :: nxChar, nyChar
    real(dp), dimension(3, rotor%nc + 1, rotor%ns + 1) :: wingMesh
    real(dp), dimension(3, rotor%nNwake + 1) :: nWakeTip
    real(dp), dimension(rotor%nNwake) :: gamRollup
    real(dp), dimension(3, rotor%nFwake + 1) :: fWakeTip
    real(dp) :: gamSum
    integer :: ib, i, j, nx, ny

    open (unit=10, file=ResultsDir// &
      & 'r'//rotor%id//'tip'//timestamp//'.plt', &
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
        & vr%vf(1)%ageAzimuthal*radToDeg, i=1,nx),j=1,ny)
      ! write(10,*) ((rotor%blade(ib)%wiP(i,j)%vr%skew,i=1,nx),j=1,ny)

      ! Near wake
      nx = rotor%nNwake
      write (nxChar, '(I5)') nx - (rotor%rowNear - 1) + 1

      ! Compute tip location for near wake
      do i = rotor%rowNear, nx
        gamRollup(i) = rotor%blade(ib)%waN(i, rotor%ns)%vr%gam
        gamSum = 0._dp
        nWakeTip(:, i) = 0._dp
        do j = rotor%rollupStart, rotor%rollupEnd
          nWakeTip(:, i) = nWakeTip(:, i) + rotor%blade(ib)%waN(i, j)%vr%vf(4)%fc(:, 1)* &
            rotor%blade(ib)%waN(i, j)%vr%gam
          gamSum = gamSum + rotor%blade(ib)%waN(i, j)%vr%gam

          ! Compute max gam for near wake filaments
          if (sign(1._dp, rotor%Omega*rotor%controlPitch(1)) > eps) then
            if (rotor%blade(ib)%waN(i, j)%vr%gam < gamRollup(i)) then    ! '<' because of negative gamma
              gamRollup(i) = rotor%blade(ib)%waN(i, j)%vr%gam
            endif
          else    ! one of Omega or pitch is negative
            if (rotor%blade(ib)%waN(i, j)%vr%gam > gamRollup(i)) then    ! '>' because of positive gamma
              gamRollup(i) = rotor%blade(ib)%waN(i, j)%vr%gam
            endif
          endif
        enddo

        if (abs(gamSum) > eps) then
          nWakeTip(:, i) = nWakeTip(:, i)/gamSum
        else
          nWakeTip(:, i) = rotor%blade(ib)%waN(i, rotor%rollupEnd)%vr%vf(4)%fc(:, 1)
        endif

      enddo

      ! For last row
      gamSum = 0._dp
      nWakeTip(:, nx + 1) = 0._dp
      do j = rotor%rollupStart, rotor%rollupEnd
        nWakeTip(:, nx + 1) = nWakeTip(:, nx + 1) + rotor%blade(ib)%waN(nx, j)%vr%vf(3)%fc(:, 1)* &
          rotor%blade(ib)%waN(nx, j)%vr%gam
        gamSum = gamSum + rotor%blade(ib)%waN(nx, j)%vr%gam
      enddo

      if (abs(gamSum) > eps) then
        nWakeTip(:, nx + 1) = nWakeTip(:, nx + 1)/gamSum
      else
        nWakeTip(:, nx + 1) = rotor%blade(ib)%waN(nx, rotor%rollupEnd)%vr%vf(3)%fc(:, 1)
      endif

      write (10, *) 'Zone I='//trim(nxChar)//' J=1    K=1  T="NearWake"'
      write (10, *) 'DATAPACKING=BLOCK'
      write (10, *) 'VARLOCATION=([4,5]=CELLCENTERED)'
      write (10, *) (nWakeTip(1, i), i=rotor%rowNear, nx + 1)
      write (10, *) (nWakeTip(2, i), i=rotor%rowNear, nx + 1)
      write (10, *) (nWakeTip(3, i), i=rotor%rowNear, nx + 1)
      write (10, *) (-1._dp*gamRollup(i), i=rotor%rowNear, nx)
      write(10,*) (rotor%blade(ib)%waN(i, rotor%ns)% &
        & vr%vf(4)%ageAzimuthal*radToDeg, i=rotor%rowNear, nx)
      !write(10,*) ((rotor%blade(ib)%waN(i,j)%vr%skew,i=rotor%rowNear,nx),j=1,ny)

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
        write(10,*) (rotor%blade(ib)%waF(i)%vf%ageAzimuthal*radToDeg, &
          & i=rotor%rowFar, nx)
        !write(10,*) (rotor%blade(ib)%waF(i)%vf%age,i=rotor%rowFar,nx)

      else  ! No far wake present

        write (nxChar, '(I5)') 2  ! Plot mesh as single redundant point
        fWakeTip(:, 1) = rotor%blade(ib)%waN(rotor%nNwake, rotor%ns)%vr%vf(3)%fc(:, 1)

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

  subroutine skew2file(timestamp, rotor)
    use classdef, only: rotor_class
    type(rotor_class), intent(in) :: rotor
    character(len=*), intent(in) :: timestamp
    integer :: ib, irow, nrow, ncol

    nrow = size(rotor%blade(1)%waN, 1)
    ncol = size(rotor%blade(1)%waN, 2)

    do ib = 1, rotor%nb
      open (unit=12, file=ResultsDir// &
        & 'r'//rotor%id// 'b'//rotor%blade(ib)%id// &
        & 'skew'//timestamp//'.csv', & 
        & action='write')

      write(12, 110) 'max', 'avg'
      do irow = nrow, rotor%rowNear, -1
        write(12, 120) maxval(rotor%blade(ib)%waN(irow, :)%vr%skew), &
          & sum(rotor%blade(ib)%waN(irow, :)%vr%skew)/ncol
      enddo

      110 format(2(A15))
      120 format(2(E15.7))
      close(12)
    enddo
  end subroutine skew2file

  subroutine force2file(timestamp, rotor)
    ! Write sec and net force to file
    use classdef, only: rotor_class
    type(rotor_class), intent(in) :: rotor
    character(len=*), intent(in) :: timestamp
    integer :: ib, ispan
    character(len=24) :: forceDimFilename
    character(len=27) :: forceNonDimFilename

    forceNonDimFilename = ResultsDir//'r'//rotor%id//'ForceNonDim.csv'
    open (unit=11, file=forceNonDimFilename, action='write', position='append')
    write (11, 100) timestamp, &
      norm2(rotor%lift) / rotor%nonDimforceDenominator, &          ! CL
      norm2(rotor%drag) / rotor%nonDimforceDenominator, &          ! CD
      norm2(rotor%liftUnsteady) / rotor%nonDimforceDenominator, &  ! CLu
      norm2(rotor%dragInduced) / rotor%nonDimforceDenominator, &   ! CDi
      norm2(rotor%dragProfile) / rotor%nonDimforceDenominator, &   ! CDo
      norm2(rotor%dragUnsteady) / rotor%nonDimforceDenominator, &  ! CDu
      rotor%forceInertial(1) / rotor%nonDimforceDenominator, &     ! CFx
      rotor%forceInertial(2) / rotor%nonDimforceDenominator, &     ! CFy
      rotor%forceInertial(3) / rotor%nonDimforceDenominator        ! CFz
    close (11)
    100 format(A, 9(E15.7))

    forceDimFilename = ResultsDir//'r'//rotor%id//'ForceDim.csv'
    open (unit=12, file=forceDimFilename, action='write', position='append')
    write (12, 101) timestamp, norm2(rotor%lift), norm2(rotor%drag), &
      rotor%lift(1), rotor%lift(2), rotor%lift(3), &     ! Lift
      rotor%drag(1), rotor%drag(2), rotor%drag(3), &     ! Drag
      rotor%forceInertial(1), rotor%forceInertial(2), rotor%forceInertial(3)  
    close (12)
    101 format(A, 11(E15.7))

    do ib = 1, rotor%nb
      open (unit=12, file=ResultsDir// &
        & 'r'//rotor%id// 'b'//rotor%blade(ib)%id// &
        & 'ForceDist'//timestamp//'.csv', & 
        & action='write', position='append')
      write (12, 202) 'secSpan', 'secCL', 'secCD', 'secCLu', &
        & 'secLift', 'secDrag', &
        & 'secArea', 'secVel', 'secChord', 'secAlpha', 'secPhi'
      ! secVel is resultant vel that the airfoil sees
      do ispan = 1, rotor%ns
        write (12, 102) dot_product(rotor%blade(ib)%secCP(:, ispan) - &
          & rotor%hubCoords, rotor%blade(ib)%yAxis), &
          & rotor%blade(ib)%secCL(ispan), &
          & rotor%blade(ib)%secCD(ispan), &
          & rotor%blade(ib)%secCLu(ispan), &
          & norm2(rotor%blade(ib)%secLift(:, ispan)), &
          & norm2(rotor%blade(ib)%secDrag(:, ispan)), &
          & rotor%blade(ib)%secArea(ispan), &
          & norm2(rotor%blade(ib)%secChordwiseResVel(:, ispan)), &
          & rotor%blade(ib)%secChord(ispan), &
          & rotor%blade(ib)%secAlpha(ispan)*radToDeg, &
          & rotor%blade(ib)%secPhi(ispan)*radToDeg
      enddo
      close (12)
    enddo
    202 format(11(A15))
    102 format(11(E15.7))
  end subroutine force2file

  subroutine dynamics2file(timestamp, rotor)
    use classdef, only: rotor_class
    character(len=*), intent(in) :: timestamp
    type(rotor_class), intent(in) :: rotor
    character(len=30) :: dynamicsFilename

    dynamicsFilename = ResultsDir//'r'//rotor%id//'bladedynamics.csv'
    open(unit=10, file=dynamicsFilename, action='write', position='append')
    write(10, 100) timestamp, rotor%blade(1)%flap*radToDeg, &
      & rotor%blade(1)%dflap*radToDeg
    close(10)

    100 format (A5, 2(E15.7))
  end subroutine dynamics2file

  subroutine inflow2file(timestamp, rotorArray, rotorNumber, directionVector)
    ! Calculate inflow velocity along directionVector 
    ! on the blades of rotor(rotorNumber) at rotor(rotorNumber)%secCP
    use classdef, only: rotor_class
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

  subroutine gamma2file(timestamp, rotor)
    ! Calculate inflow velocity along directionVector on the blades of rotor(rotorNumber)
    ! at rotor(rotorNumber)%secCP
    use classdef, only: rotor_class
    character(len=*), intent(in) :: timestamp
    type(rotor_class), intent(in) :: rotor
    integer :: ib, ic, is
    character(len=2) :: rowNumberChar

    ! Write to file
    do ib = 1, rotor%nb
      open (unit=12, &
        & file=ResultsDir// &
        & 'r'//rotor%id//'b'//rotor%blade(ib)%id// &
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
  !         & rotor%blade(ib)%yAxis), rotor%blade(ib)%secAlpha(is)*radToDeg
  !     enddo
  !     close (12)
  !   enddo
  !   100 format (2(A15))
  !   101 format (2(F15.7))

  ! end subroutine alpha2file

end module libPostprocess

