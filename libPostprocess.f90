module libPostprocess
  use rotor_classdef

contains

  subroutine init_plots(numOfRotors)
    ! Initialise headers for plot files
    integer, intent(in) :: numOfRotors
    character(len=24) :: forceFilename
    integer :: rotorNumber
    character(len=2) :: rotorNumberChar

    do rotorNumber=1,numOfRotors
      write(rotorNumberChar,'(I0.2)') rotorNumber
      forceFilename='Results/r'//rotorNumberChar//'forceHist.txt'

      ! Add data headers
      open(unit=11,file=forceFilename,action='write')
      write(11,'(A)') 'timestamp(iters)    CT    rotorThrust    bladeThrust1    bladeThrust2...'
    enddo
    close(11)
  end subroutine init_plots

  subroutine rotor2file(timestamp,rotor)
    ! Plot rotor geometry and wake to file
    type(rotor_class), intent(in) :: rotor
    character(len=*), intent(in) :: timestamp
    character(len=5) :: nxChar, nyChar
    real(dp), dimension(3,rotor%nc+1,rotor%ns+1) :: wingMesh  
    real(dp), dimension(3,rotor%nNwake+1,rotor%ns+1) :: wakeMesh  
    real(dp), dimension(3,rotor%nFwake+1) :: wakeTip   ! Optimise this by only initialising reqd size
    integer :: i,j,nx,ny,ib

    open(unit=10,file='Results/wingNwake'//timestamp//'.plt',position='append')

    write(10,*) 'Title = "Wing and Wake"'
    write(10,*) 'VARIABLES = "X" "Y" "Z" "GAM"'! "Var5" "Var6"'

    do ib=1,rotor%nb
      ! Wing
      nx=rotor%nc
      ny=rotor%ns
      write(nxChar,'(I5)') nx+1
      write(nyChar,'(I5)') ny+1

      do j=1,ny
        do i=1,nx
          wingMesh(:,i,j)=rotor%blade(ib)%wiP(i,j)%pc(:,1)
        enddo
      enddo
      do i=1,nx
        wingMesh(:,i,ny+1)=rotor%blade(ib)%wiP(i,ny)%pc(:,4)
      enddo
      do j=1,ny
        wingMesh(:,nx+1,j)=rotor%blade(ib)%wiP(nx,j)%pc(:,2)
      enddo
      wingMesh(:,nx+1,ny+1)=rotor%blade(ib)%wiP(nx,ny)%pc(:,3)

      write(10,*) 'Zone I='//trim(nxChar)//' J='//trim(nyChar)//' K=1  T="Blade"'
      write(10,*) 'DATAPACKING=BLOCK'
      write(10,*) 'VARLOCATION=([4]=CELLCENTERED)'!,[5]=CELLCENTERED,[6]=CELLCENTERED)'
      write(10,*) ((wingMesh(1,i,j),i=1,nx+1),j=1,ny+1)
      write(10,*) ((wingMesh(2,i,j),i=1,nx+1),j=1,ny+1)
      write(10,*) ((wingMesh(3,i,j),i=1,nx+1),j=1,ny+1)
      write(10,*) ((-1._dp*rotor%blade(ib)%wiP(i,j)%vr%gam,i=1,nx),j=1,ny)
      !write(10,*) ((rotor%blade(ib)%wiP(i,j)%vr%skew,i=1,nx),j=1,ny)
      !write(10,*) ((rotor%blade(ib)%wiP(i,j)%vr%skew,i=1,nx),j=1,ny)

      ! Near wake 
      nx=rotor%nNwake
      ny=rotor%ns
      write(nxChar,'(I5)') nx-(rotor%rowNear-1)+1
      write(nyChar,'(I5)') ny+1

      !Check if necessary - $omp parallel do collapse(2)
      do j=1,ny
        do i=rotor%rowNear,nx
          wakeMesh(:,i,j)=rotor%blade(ib)%waP(i,j)%vr%vf(1)%fc(:,1)
        enddo
      enddo
      !Check if necessary -$omp end parallel do
      do i=rotor%rowNear,nx
        wakeMesh(:,i,ny+1)=rotor%blade(ib)%waP(i,ny)%vr%vf(4)%fc(:,1)
      enddo
      do j=1,ny
        wakeMesh(:,nx+1,j)=rotor%blade(ib)%waP(nx,j)%vr%vf(2)%fc(:,1)
      enddo
      wakeMesh(:,nx+1,ny+1)=rotor%blade(ib)%waP(nx,ny)%vr%vf(3)%fc(:,1)

      write(10,*) 'Zone I='//trim(nxChar)//' J='//trim(nyChar)//' K=1  T="NearWake"'
      write(10,*) 'DATAPACKING=BLOCK'
      write(10,*) 'VARLOCATION=([4]=CELLCENTERED)'!,[5]=CELLCENTERED,[6]=CELLCENTERED)'
      write(10,*) ((wakeMesh(1,i,j),i=rotor%rowNear,nx+1),j=1,ny+1)
      write(10,*) ((wakeMesh(2,i,j),i=rotor%rowNear,nx+1),j=1,ny+1)
      write(10,*) ((wakeMesh(3,i,j),i=rotor%rowNear,nx+1),j=1,ny+1)
      write(10,*) ((-1._dp*rotor%blade(ib)%waP(i,j)%vr%gam,i=rotor%rowNear,nx),j=1,ny)
      !write(10,*) ((rotor%blade(ib)%waP(i,j)%vr%skew,i=rotor%rowNear,nx),j=1,ny)
      !write(10,*) ((rotor%blade(ib)%waP(i,j)%vr%skew,i=rotor%rowNear,nx),j=1,ny)

      ! Far wake 
      nx=rotor%nFwake
      if (rotor%rowFar .ne. 0) then
        write(nxChar,'(I5)') nx-(rotor%rowFar-1)+1

        !Check if necessary - $omp parallel do collapse(2)
        do i=rotor%rowFar,nx
          wakeTip(:,i)=rotor%blade(ib)%waF(i)%vf%fc(:,2)
        enddo
        wakeTip(:,nx+1)=rotor%blade(ib)%waF(rotor%nFwake)%vf%fc(:,1)
        !Check if necessary -$omp end parallel do

        write(10,*) 'Zone I='//trim(nxChar)//' J=1   K=1   T="FarWake"'
        write(10,*) 'DATAPACKING=BLOCK'
        write(10,*) 'VARLOCATION=([4]=CELLCENTERED)'!,[5]=CELLCENTERED,[6]=CELLCENTERED)'
        write(10,*) (wakeTip(1,i),i=rotor%rowFar,nx+1)
        write(10,*) (wakeTip(2,i),i=rotor%rowFar,nx+1)
        write(10,*) (wakeTip(3,i),i=rotor%rowFar,nx+1)
        write(10,*) (-1._dp*rotor%blade(ib)%waF(i)%gam,i=rotor%rowFar,nx)
        !write(10,*) (rotor%blade(ib)%waF(i)%vf%rVc,i=rotor%rowFar,nx)
        !write(10,*) (rotor%blade(ib)%waF(i)%vf%age,i=rotor%rowFar,nx)

      else  ! No far wake present

        write(nxChar,'(I5)') 2  ! Plot mesh as single redundant point
        wakeTip(:,1) = rotor%blade(ib)%waP(rotor%nNwake,rotor%ns)%vr%vf(3)%fc(:,1)

        write(10,*) 'Zone I='//trim(nxChar)//' J=1   K=1   T="FarWake"'
        write(10,*) 'DATAPACKING=BLOCK'
        write(10,*) 'VARLOCATION=([4]=CELLCENTERED)'!,[5]=CELLCENTERED,[6]=CELLCENTERED)'
        write(10,*) wakeTip(1,1),wakeTip(1,1)
        write(10,*) wakeTip(2,1),wakeTip(2,1)
        write(10,*) wakeTip(3,1),wakeTip(3,1)
        write(10,*) 0._dp
        !write(10,*) 0._dp
        !write(10,*) 0._dp
      endif

    enddo

    close(10)
  end subroutine rotor2file

  subroutine filaments2file(timestamp,rotor)
    ! Write filaments to file for using with grid-based plots
    type(rotor_class), intent(in), dimension(:) :: rotor
    character(len=*), intent(in) :: timestamp

    integer :: nr, nvrWing, nvrNwake, nvfFwake
    integer :: ir, ib, irow, icol, indx
    type(vr_class), allocatable, dimension(:) :: vrWing, vrNwake
    type(vf_class), allocatable, dimension(:) :: vfFwake
    real(dp), allocatable, dimension(:) :: gamFwake

    nvrWing=0
    nvrNwake=0
    nvfFwake=0

    nr=size(rotor)

    do ir=1,nr
      if (rotor(ir)%rowFar .eq. 0) error stop 'ERROR: Use filaments2file() only after development of far wake'
    enddo

    ! Compute number of each filaments
    do ir=1,nr
      nvrWing=nvrWing+rotor(ir)%nb*(rotor(ir)%nc*rotor(ir)%ns)
      nvrNwake=nvrNwake+rotor(ir)%nb*(rotor(ir)%nNwake*rotor(ir)%ns)
      nvfFwake=nvfFwake+(rotor(ir)%nFwake-rotor(ir)%rowFar+1)*rotor(ir)%nb
    enddo

    ! Allocate filaments
    allocate(vrWing(nvrWing))
    allocate(vrNwake(nvrNwake))
    allocate(vfFwake(nvfFwake))
    allocate(gamFwake(nvfFwake))

    ! Extract filament properties
    ! from wing
    indx=1
    do ir=1,nr
      do ib=1,rotor(ir)%nb
        do icol=1,rotor(ir)%ns
          do irow=1,rotor(ir)%nc
            vrWing(indx)=rotor(ir)%blade(ib)%wiP(irow,icol)%vr
            indx=indx+1
          enddo
        enddo
      enddo
    enddo

    ! from Nwake
    indx=1
    do ir=1,nr
      do ib=1,rotor(ir)%nb
        do icol=1,rotor(ir)%ns
          do irow=1,rotor(ir)%nNwake
            vrNwake(indx)=rotor(ir)%blade(ib)%waP(irow,icol)%vr
            indx=indx+1
          enddo
        enddo
      enddo
    enddo


    ! from Fwake
    indx=1
    do ir=1,nr
      do ib=1,rotor(ir)%nb
        do irow=rotor(ir)%rowFar,rotor(ir)%nFwake
          vfFwake(indx)=rotor(ir)%blade(ib)%waF(irow)%vf
          gamFwake(indx)=rotor(ir)%blade(ib)%waF(irow)%gam
          indx=indx+1
        enddo
      enddo
    enddo

    ! Write to filamentsXXXXX.dat binary file
    open(unit=10,file='Results/filaments'//timestamp//'.dat',form='unformatted')
    write(10) nvrWing
    write(10) nvrNwake
    write(10) nvfFwake
    write(10) vrWing, vrNwake
    write(10) vfFwake, gamFwake
    close(10)

    ! Deallocate filaments
    deallocate(vrWing)
    deallocate(vrNwake)
    deallocate(vfFwake)
    deallocate(gamFwake)

  end subroutine filaments2file

  subroutine mesh2file(wing_array,wake_array,filename)
    ! Obsolete 
    type(wingpanel_class), intent(in), dimension(:,:) :: wing_array
    type(Nwake_class), intent(in), dimension(:,:) :: wake_array
    character(len=*), intent(in) :: filename
    character(len=5) :: nxChar, nyChar
    real(dp), dimension(3,size(wing_array,1)+1,size(wing_array,2)+1) :: wingMesh  
    real(dp), dimension(3,size(wake_array,1)+1,size(wake_array,2)+1) :: wakeMesh  
    integer :: i,j,nx,ny

    nx=size(wing_array,1)
    ny=size(wing_array,2)
    write(nxChar,'(I5)') nx+1
    write(nyChar,'(I5)') ny+1

    open(unit=10,file=filename,position='append')
    do j=1,ny
      do i=1,nx
        wingMesh(:,i,j)=wing_array(i,j)%pc(:,1)
      enddo
    enddo
    do i=1,nx
      wingMesh(:,i,ny+1)=wing_array(i,ny)%pc(:,4)
    enddo
    do j=1,ny
      wingMesh(:,nx+1,j)=wing_array(nx,j)%pc(:,2)
    enddo
    wingMesh(:,nx+1,ny+1)=wing_array(nx,ny)%pc(:,3)

    write(10,*) 'Title = "Panel array"'
    write(10,*) 'VARIABLES = "X" "Y" "Z" "GAM"'
    write(10,*) 'Zone I='//trim(nxChar)//' J='//trim(nyChar)//' K=1  T="Wing"'
    write(10,*) 'DATAPACKING=BLOCK'
    write(10,*) 'VARLOCATION=([4]=CELLCENTERED)'!,[5]=CELLCENTERED)'
    write(10,*) ((wingMesh(1,i,j),i=1,nx+1),j=1,ny+1)
    write(10,*) ((wingMesh(2,i,j),i=1,nx+1),j=1,ny+1)
    write(10,*) ((wingMesh(3,i,j),i=1,nx+1),j=1,ny+1)
    write(10,*) ((-1._dp*wing_array(i,j)%vr%gam,i=1,nx),j=1,ny)

    nx=size(wake_array,1)
    ny=size(wake_array,2)
    write(nxChar,'(I5)') nx+1
    write(nyChar,'(I5)') ny+1

    !Check if necessary - $omp parallel do collapse(2)
    do j=1,ny
      do i=1,nx
        wakeMesh(:,i,j)=wake_array(i,j)%vr%vf(1)%fc(:,1)
      enddo
    enddo
    !Check if necessary -$omp end parallel do
    do i=1,nx
      wakeMesh(:,i,ny+1)=wake_array(i,ny)%vr%vf(4)%fc(:,1)
    enddo
    do j=1,ny
      wakeMesh(:,nx+1,j)=wake_array(nx,j)%vr%vf(2)%fc(:,1)
    enddo
    wakeMesh(:,nx+1,ny+1)=wake_array(nx,ny)%vr%vf(3)%fc(:,1)

    write(10,*) 'Zone I='//trim(nxChar)//' J='//trim(nyChar)//' K=1  T="Wake"'
    write(10,*) 'DATAPACKING=BLOCK'
    write(10,*) 'VARLOCATION=([4]=CELLCENTERED)'!,[5]=CELLCENTERED)'
    write(10,*) ((wakeMesh(1,i,j),i=1,nx+1),j=1,ny+1)
    write(10,*) ((wakeMesh(2,i,j),i=1,nx+1),j=1,ny+1)
    write(10,*) ((wakeMesh(3,i,j),i=1,nx+1),j=1,ny+1)
    write(10,*) ((-1._dp*wake_array(i,j)%vr%gam,i=1,nx),j=1,ny)

    close(10)
  end subroutine mesh2file

  subroutine wingverify(wing_array)
    ! Write wing geometry to file in detail
    ! For verifying orientation of wing panels, bound vortex rings and CPs
    type(wingpanel_class), intent(in), dimension(:,:) :: wing_array
    character(len=5) :: nxChar, nyChar
    real(dp), dimension(3,size(wing_array,1)+1,size(wing_array,2)+1) :: wingMesh  
    integer :: i,j,nx,ny

    nx=size(wing_array,1)
    ny=size(wing_array,2)
    write(nxChar,'(I5)') nx+1
    write(nyChar,'(I5)') ny+1

    do j=1,ny
      do i=1,nx
        wingMesh(:,i,j)=wing_array(i,j)%pc(:,1)
      enddo
    enddo
    do i=1,nx
      wingMesh(:,i,ny+1)=wing_array(i,ny)%pc(:,4)
    enddo
    do j=1,ny
      wingMesh(:,nx+1,j)=wing_array(nx,j)%pc(:,2)
    enddo
    wingMesh(:,nx+1,ny+1)=wing_array(nx,ny)%pc(:,3)

    open(unit=10,file='Results/wingPC.plt',position='append')
    write(10,*) 'Title = "Panel Vertices"'
    write(10,*) 'VARIABLES = "X" "Y" "Z"'
    write(10,*) 'Zone I='//trim(nxChar)//' J='//trim(nyChar)//' K=1  T="Panel Vertices"'
    write(10,*) 'DATAPACKING=BLOCK'
    write(10,*) ((wingMesh(1,i,j),i=1,nx+1),j=1,ny+1)
    write(10,*) ((wingMesh(2,i,j),i=1,nx+1),j=1,ny+1)
    write(10,*) ((wingMesh(3,i,j),i=1,nx+1),j=1,ny+1)
    close(10)

    write(nxChar,'(I5)') nx
    write(nyChar,'(I5)') ny
    do j=1,ny
      do i=1,nx
        wingMesh(:,i,j)=wing_array(i,j)%CP
      enddo
    enddo
    open(unit=11,file='Results/wingCP.plt',position='append')
    write(11,*) 'Title = "Coll. points"'
    write(11,*) 'VARIABLES = "X" "Y" "Z"'
    write(11,*) 'Zone I='//trim(nxChar)//' J='//trim(nyChar)//' K=1  T="Coll. points"'
    write(11,*) 'DATAPACKING=BLOCK'
    write(11,*) ((wingMesh(1,i,j),i=1,nx),j=1,ny)
    write(11,*) ((wingMesh(2,i,j),i=1,nx),j=1,ny)
    write(11,*) ((wingMesh(3,i,j),i=1,nx),j=1,ny)
    close(11)

    write(nxChar,'(I5)') nx+1
    write(nyChar,'(I5)') ny+1
    do j=1,ny
      do i=1,nx
        wingMesh(:,i,j)=wing_array(i,j)%vr%vf(1)%fc(:,1)
      enddo
    enddo
    do i=1,nx
      wingMesh(:,i,ny+1)=wing_array(i,ny)%vr%vf(4)%fc(:,1)
    enddo
    do j=1,ny
      wingMesh(:,nx+1,j)=wing_array(nx,j)%vr%vf(2)%fc(:,1)
    enddo
    wingMesh(:,nx+1,ny+1)=wing_array(nx,ny)%vr%vf(3)%fc(:,1)

    open(unit=12,file='Results/wingVR.plt',position='append')
    write(12,*) 'Title = "Vortex Rings"'
    write(12,*) 'VARIABLES = "X" "Y" "Z"'
    write(12,*) 'Zone I='//trim(nxChar)//' J='//trim(nyChar)//' K=1  T="Vortex Rings"'
    write(12,*) 'DATAPACKING=BLOCK'
    write(12,*) ((wingMesh(1,i,j),i=1,nx+1),j=1,ny+1)
    write(12,*) ((wingMesh(2,i,j),i=1,nx+1),j=1,ny+1)
    write(12,*) ((wingMesh(3,i,j),i=1,nx+1),j=1,ny+1)
    close(12)
  end subroutine wingverify

  ! subroutine tip2file(wing_array,wake_array,filename)
  !   type(wingpanel_class), intent(in), dimension(:,:) :: wing_array
  !   type(Nwake_class), intent(in), dimension(:,:) :: wake_array
  !   character(len=*), intent(in) :: filename
  !   character(len=5) :: nxChar, nyChar
  !   real(dp), dimension(3,size(wing_array,1)+1,size(wing_array,2)+1) :: wingMesh  
  !   real(dp), dimension(3,size(wake_array,1)+1) :: wakeTip  
  !   integer :: i,j,nx,ny

  !   nx=size(wing_array,1)
  !   ny=size(wing_array,2)
  !   write(nxChar,'(I5)') nx+1
  !   write(nyChar,'(I5)') ny+1

  !   open(unit=10,file=filename,position='append')
  !   do j=1,ny
  !     do i=1,nx
  !       wingMesh(:,i,j)=wing_array(i,j)%pc(:,1)
  !     enddo
  !   enddo
  !   do i=1,nx
  !     wingMesh(:,i,ny+1)=wing_array(i,ny)%pc(:,4)
  !   enddo
  !   do j=1,ny
  !     wingMesh(:,nx+1,j)=wing_array(nx,j)%pc(:,2)
  !   enddo
  !   wingMesh(:,nx+1,ny+1)=wing_array(nx,ny)%pc(:,3)

  !   write(10,*) 'Title = "Panel array"'
  !   write(10,*) 'VARIABLES = "X" "Y" "Z"'
  !   write(10,*) 'Zone I='//trim(nxChar)//' J='//trim(nyChar)//' K=1  T="Wing"'
  !   write(10,*) 'DATAPACKING=BLOCK'
  !   write(10,*) ((wingMesh(1,i,j),i=1,nx+1),j=1,ny+1)
  !   write(10,*) ((wingMesh(2,i,j),i=1,nx+1),j=1,ny+1)
  !   write(10,*) ((wingMesh(3,i,j),i=1,nx+1),j=1,ny+1)

  !   ! Wake root
  !   nx=size(wake_array,1)
  !   ny=size(wake_array,2)
  !   write(nxChar,'(I5)') nx+1

  !   do i=1,nx
  !     wakeTip(:,i)=wake_array(i,1)%vr%vf(1)%fc(:,1)
  !   enddo
  !   wakeTip(:,nx+1)=wake_array(nx,1)%vr%vf(2)%fc(:,1)

  !   write(10,*) 'Zone I='//trim(nxChar)//' J=1   K=1  T="wake_root"'
  !   write(10,*) 'DATAPACKING=BLOCK'
  !   write(10,*) (wakeTip(1,i),i=1,nx+1)
  !   write(10,*) (wakeTip(2,i),i=1,nx+1)
  !   write(10,*) (wakeTip(3,i),i=1,nx+1)

  !   ! Wake tip
  !   do i=1,nx
  !     wakeTip(:,i)=wake_array(i,ny)%vr%vf(4)%fc(:,1)
  !   enddo
  !   wakeTip(:,nx+1)=wake_array(nx,ny)%vr%vf(3)%fc(:,1)

  !   write(10,*) 'Zone I='//trim(nxChar)//' J=1   K=1  T="wakeTip"'
  !   write(10,*) 'DATAPACKING=BLOCK'
  !   write(10,*) (wakeTip(1,i),i=1,nx+1)
  !   write(10,*) (wakeTip(2,i),i=1,nx+1)
  !   write(10,*) (wakeTip(3,i),i=1,nx+1)
  !   close(10)
  ! end subroutine tip2file

  subroutine force2file(timestamp,rotor,rotorNumber,directionVector)
    ! Write sectional and net force to file
    type(rotor_class), intent(in) :: rotor
    character(len=*), intent(in) :: timestamp
    integer, intent(in) :: rotorNumber
    real(dp),   intent(in), dimension(3) :: directionVector
    character(len=2) :: rotorNumberChar, bladeNumberChar
    real(dp) :: rotorForce
    real(dp), dimension(rotor%nb) :: bladeForce
    integer :: ib, ispan
    character(len=24) :: forceFilename

    rotorForce = dot_product(rotor%Force,directionVector)
    do ib=1,rotor%nb
      bladeForce(ib) = dot_product(rotor%blade(ib)%Force,directionVector)
    enddo

    write(rotorNumberChar,'(I0.2)') rotorNumber
    forceFilename='Results/r'//rotorNumberChar//'forceHist.txt'

    open(unit=11,file=forceFilename,action='write',position='append')
    write(11,100) timestamp,rotorForce/rotor%nonDimForceDenominator, rotorForce, (bladeForce(ib),ib=1,rotor%nb)
    close(11)
    100 format(A,15(E15.7))

    open(unit=12,file='Results/r'//rotorNumberChar//'forceDist'//timestamp//'.curve',action='write')
    do ib=1,rotor%nb
      write(bladeNumberChar,'(I0.2)') ib
      write(12,*) '# Blade'//bladeNumberChar
      do ispan=1,rotor%ns
        write(12,*) norm2(rotor%hubCoords-rotor%blade(ib)%wiP(1,ispan)%CP), &
          dot_product(rotor%blade(ib)%sectionalForce(:,ispan),directionVector)
      enddo
    enddo
    close(12)
  end subroutine force2file

  subroutine inflow2file(timestamp,rotorArray,rotorNumber,directionVector)
    ! Calculate inflow velocity along directionVector on the blades of rotor(rotorNumber)
    ! at rotor(rotorNumber)%inflowLocations
    character(len=*), intent(in) :: timestamp
    type(rotor_class), intent(inout), dimension(:) :: rotorArray
    integer, intent(in) :: rotorNumber
    real(dp), intent(in), dimension(3) :: directionVector
    character(len=2) :: rotorNumberChar, bladeNumberChar
    integer :: il,ir,ib

    real(dp), dimension(3) :: P
    real(dp), dimension(rotorArray(rotorNumber)%ns,rotorArray(rotorNumber)%nb) :: inflowVel

    inflowVel=0._dp
    do ir=1,size(rotorArray)
      do ib=1,rotorArray(rotorNumber)%nb
        do il=1,rotorArray(rotorNumber)%ns
          P=rotorArray(rotorNumber)%blade(ib)%inflowLocations(:,il)
          inflowVel(il,ib)=inflowVel(il,ib)+dot_product(rotorArray(ir)%vind_bywing(P),directionVector) 
          inflowVel(il,ib)=inflowVel(il,ib)-dot_product(rotorArray(ir)%vind_bywing_boundVortices(P),directionVector) 
          inflowVel(il,ib)=inflowVel(il,ib)+dot_product(rotorArray(ir)%vind_bywake(P),directionVector) 
        enddo
      enddo
    enddo

    ! Write to file
    write(rotorNumberChar,'(I0.2)') rotorNumber
    open(unit=12,file='Results/r'//rotorNumberChar//'inflowDist'//timestamp//'.curve',action='write')
    do ib=1,rotorArray(rotorNumber)%nb
      write(bladeNumberChar,'(I0.2)') ib
      write(12,*) '# Blade'//bladeNumberChar
      do il=1,rotorArray(rotorNumber)%ns
        write(12,*) norm2(rotorArray(rotorNumber)%hubCoords-rotorArray(rotorNumber)%blade(ib)%inflowLocations(:,il)), &
          inflowVel(il,ib)
      enddo
    enddo
    close(12)

  end subroutine inflow2file

  subroutine gamma2file(timestamp,rotor,rotorNumber)
    ! Calculate inflow velocity along directionVector on the blades of rotor(rotorNumber)
    ! at rotor(rotorNumber)%inflowLocations
    character(len=*), intent(in) :: timestamp
    type(rotor_class), intent(in) :: rotor
    integer, intent(in) :: rotorNumber
    integer :: ib, ic, is
    character(len=2) :: rotorNumberChar, bladeNumberChar, rowNumberChar

    ! Write to file
    write(rotorNumberChar,'(I0.2)') rotorNumber
    open(unit=12,file='Results/r'//rotorNumberChar//'gammaDist'//timestamp//'.curve',action='write')
    do ib=1,rotor%nb
      write(bladeNumberChar,'(I0.2)') ib
      ic=1
      write(rowNumberChar,'(I0.2)') ic
      write(12,*) '# Blade'//bladeNumberChar//'Row'//rowNumberChar
      do is=1,rotor%ns
        write(12,*) norm2(rotor%hubCoords-rotor%blade(ib)%wiP(ic,is)%cp), &
          -1._dp*rotor%blade(ib)%wiP(ic,is)%vr%gam
      enddo
      do ic=2,rotor%nc
        write(rowNumberChar,'(I0.2)') ic
        write(12,*) '# Blade'//bladeNumberChar//'Row'//rowNumberChar
        do is=1,rotor%ns
          write(12,*) norm2(rotor%hubCoords-rotor%blade(ib)%wiP(ic,is)%cp), &
            -1._dp*(rotor%blade(ib)%wiP(ic,is)%vr%gam-rotor%blade(ib)%wiP(ic-1,is)%vr%gam)
        enddo
      enddo
    enddo
    close(12)

  end subroutine gamma2file

  subroutine alpha2file(timestamp,rotor,rotorNumber)
    ! Write sectional alpha to file
    character(len=*), intent(in) :: timestamp
    type(rotor_class), intent(inout) :: rotor
    integer, intent(in) :: rotorNumber
    integer :: ib, is
    character(len=2) :: rotorNumberChar, bladeNumberChar

    ! Write to file
    write(rotorNumberChar,'(I0.2)') rotorNumber
    open(unit=12,file='Results/r'//rotorNumberChar//'alphaDist'//timestamp//'.curve',action='write')
    do ib=1,rotor%nb
      write(bladeNumberChar,'(I0.2)') ib
      write(12,*) '# Blade'//bladeNumberChar
      do is=1,rotor%ns
        write(12,*) norm2(rotor%hubCoords-rotor%blade(ib)%wiP(1,is)%cp), &
          rotor%blade(ib)%sectionalAlpha(is)*180._dp/pi
      enddo
    enddo
    close(12)

  end subroutine alpha2file

end module libPostprocess
