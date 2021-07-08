!------+--------------------+------|
! ++++ | MODULE DEFINITIONS | ++++ |
!------+--------------------+------|
module classdef
  use libMath, only: dp, pi, eps
  use libC81, only: C81_class
  implicit none

  real(dp), parameter :: tol = 1.E-6
  real(dp), parameter :: invTol2 = 1.E06
  real(dp), parameter :: inv4pi = 0.25_dp/pi
  real(dp), parameter :: twoPi = 2.0_dp*pi

  type switches_class
    integer :: ntSub, ntSubInit
    integer :: spanSpacing, chordSpacing
    integer :: wakeDissipation, wakeStrain, wakeBurst
    integer :: slowStart, slowStartNt
    integer :: wakeTipPlot, wakePlot, gridPlot
    integer :: rotorForcePlot
    integer :: fdScheme, probe, nProbes
    integer :: initWakeVelNt
    integer :: restartFromNt, restartWriteNt
  end type switches_class

  type vf_class
    real(dp), dimension(3, 2) :: fc  ! filament coords (xyz,1:2)
    real(dp) :: l0 = 0._dp    ! original length
    real(dp) :: lc = 0._dp    ! current length
    real(dp) :: rVc0 = 0._dp ! initial vortex core radius
    real(dp) :: rVc = 0._dp  ! current vortex core radius
    real(dp) :: age = 0._dp   ! vortex age (in s)
    real(dp) :: ageAzimuthal = 0._dp   ! vortex age (in radians)
  contains
    procedure :: vind => vf_vind
    procedure :: calclength => vf_calclength
    procedure :: strain => vf_strain
  end type vf_class

  type vr_class
    type(vf_class), dimension(4) :: vf
    real(dp) :: gam
    real(dp) :: skew
  contains
    procedure :: vind => vr_vind
    procedure :: vindSource => vr_vindSource
    procedure :: assignP => vr_assignP
    procedure :: shiftdP => vr_shiftdP
    procedure :: rot => vr_rot
    procedure :: calclength => vr_calclength
    procedure :: strain => vr_strain
    procedure :: decay => vr_decay
    procedure :: calc_skew => calc_skew
    procedure :: burst => vr_burst
    procedure :: getInteriorAngles => vr_getInteriorAngles
    procedure :: getMedianAngle => vr_getMedianAngle
    procedure :: getBimedianCos => vr_getBimedianCos
    procedure :: mirror => vr_mirror
  end type vr_class

  type wingpanel_class
    ! Panel coordinates
    ! o-------> Y along span
    ! |
    ! |   1-------4  - Leading edge
    ! |   |   4   |
    ! |   |1     3|
    ! |   |   2   |
    ! |   2-------3  - Trailing edge
    ! |
    ! V X along chord
    type(vr_class) :: vr
    real(dp) :: gamPrev
    real(dp) :: gamTrapz
    real(dp), dimension(3, 4) :: pc    ! panel coords
    real(dp), dimension(3) :: cp      ! coll point coords
    real(dp), dimension(3) :: nCap    ! unit normal vector
    real(dp), dimension(3) :: tauCapChord ! unit tangential vector along chord
    real(dp), dimension(3) :: tauCapSpan  ! unit tangential vector along span
    ! local vel at CP 
    real(dp), dimension(3) :: velCP    ! excluding wing vortices
    real(dp), dimension(3) :: velCPTotal  ! including wing vortices
    real(dp), dimension(3) :: velCPm  ! rel. inertial vel at CP (due to motion)
    real(dp), dimension(3) :: normalForce ! normalForce vector (inertial frame)
    real(dp), dimension(3) :: normalForceUnsteady 
    real(dp), dimension(3) :: chordwiseResVel
    real(dp) :: velPitch             ! pitch velocity
    real(dp) :: delP                  ! Pressure difference at panel
    real(dp) :: delPUnsteady          ! Unsteady pressure at panel
    real(dp) :: delDiConstant         ! Induced drag (constant part) at panel
    real(dp) :: delDiUnsteady         ! Induced drag (unsteady part) at panel
    real(dp) :: meanChord, meanSpan ! Panel mean dimensions
    real(dp) :: panelArea            ! Panel area for computing lift
    real(dp) :: rHinge ! dist to point about which pitching occurs (LE of wing)
    real(dp) :: alpha                 ! local angle of attack
  contains
    procedure :: assignP => wingpanel_assignP
    procedure :: calcCP => wingpanel_calcCP
    procedure :: calcN => wingpanel_calcN
    procedure :: calcTau => wingpanel_calcTau
    procedure :: rot => wingpanel_rot
    procedure :: shiftdP => wingpanel_shiftdP
    procedure :: calc_chordwiseResVel => wingpanel_calc_chordwiseResVel
    procedure :: calc_area => wingpanel_calc_area
    procedure :: calc_mean_dimensions => wingpanel_calc_mean_dimensions
    procedure :: isCPinsidecore => wingpanel_isCPinsidecore
  end type wingpanel_class

  type Nwake_class
    ! VR coordinates
    ! o-------> Y along span
    ! |
    ! |   1-------4
    ! |   |   4   |
    ! |   |1     3|
    ! |   |   2   |
    ! |   2-------3
    ! |
    ! V X along chord

    type(vr_class) :: vr
  end type Nwake_class

  type Fwake_class
    ! VF coordinates
    ! o-------> Y along span
    ! |
    ! |   2 - Leading edge
    ! |   |
    ! |   |1
    ! |   |
    ! |   1 - Trailing edge
    ! |
    ! V X along chord
    type(vf_class) :: vf
    real(dp) :: gam = 0._dp
  contains
    procedure :: shiftdP => Fwake_shiftdP
    procedure :: assignP => Fwake_assignP
    procedure :: rot => Fwake_rot
    procedure :: decay => Fwake_decay
    procedure :: mirror => Fwake_mirror
  end type Fwake_class

  type pFwake_class
    !! Prescribed far wake
    ! 10 revs with 15 deg vortex filaments
    type(Fwake_class), dimension(240) :: waF
    real(dp), dimension(3, 241) :: coords
    real(dp) :: nRevs = 10.0
    real(dp) :: helixPitch = 0._dp
    real(dp) :: helixRadius = 0._dp
    real(dp) :: relaxFactor = 0.5_dp
    logical :: isClockwiseRotor = .True.
    logical :: isPresent = .false.
  contains
    procedure :: update => pFwake_update
    procedure :: rot_wake_axis => pFwake_rot_wake_axis
  end type pFwake_class

  type blade_class
    type(wingpanel_class), allocatable, dimension(:, :) :: wiP  ! Wing panel
    type(Nwake_class), allocatable, dimension(:, :) :: waN  ! Near wake
    type(Fwake_class), allocatable, dimension(:) :: waF  ! Far wake
    type(pFwake_class) :: wapF  ! Prescribed far wake
    type(pFwake_class) :: wapFPredicted  ! Prescribed far wake
    type(Nwake_class), allocatable, dimension(:, :) :: waNPredicted
    type(Fwake_class), allocatable, dimension(:) :: waFPredicted
    type(C81_class), allocatable, dimension(:) :: C81
    integer :: nc, ns
    real(dp) :: theta, psi, pivotLE
    real(dp), dimension(3) :: forceInertial
    real(dp), dimension(3) :: lift, drag
    real(dp), dimension(3) :: dragInduced, dragProfile
    real(dp), dimension(3) :: liftUnsteady, dragUnsteady
    integer, allocatable, dimension(:) :: airfoilNo
    character(len=30), allocatable, dimension(:) :: airfoilFile
    real(dp), allocatable, dimension(:) :: airfoilSectionLimit
    real(dp), allocatable, dimension(:, :, :) :: velNwake, velNwake1
    real(dp), allocatable, dimension(:, :, :) :: velNwake2, velNwake3
    real(dp), allocatable, dimension(:, :, :) :: velNwakePredicted
    real(dp), allocatable, dimension(:, :, :) :: velNwakeStep
    real(dp), allocatable, dimension(:, :) :: velFwake, velFwake1
    real(dp), allocatable, dimension(:, :) :: velFwake2, velFwake3
    real(dp), allocatable, dimension(:, :) :: velFwakePredicted
    real(dp), allocatable, dimension(:, :) :: velFwakeStep
    integer :: stlNodesCols  ! Cols of stlNodes array
    real(dp), allocatable, dimension(:, :) :: stlNodes ! (3, stlNodesCols)
    ! 3 Vertices of each element in stlElementNodes
    integer, allocatable, dimension(:, :) :: stlElementNodes
    real(dp), dimension(3) :: xAxis, yAxis, zAxis
    ! Sectional quantities
    real(dp), allocatable, dimension(:) :: secChord, secArea
    real(dp), allocatable, dimension(:, :) :: secForceInertial
    real(dp), allocatable, dimension(:, :) :: secLift, secDrag
    real(dp), allocatable, dimension(:, :) :: secLiftDir, secDragDir
    real(dp), allocatable, dimension(:, :) :: secDragInduced, secDragProfile
    real(dp), allocatable, dimension(:, :) :: secLiftUnsteady, secDragUnsteady
    real(dp), allocatable, dimension(:, :) :: secTauCapChord, secTauCapSpan
    real(dp), allocatable, dimension(:, :) :: secNormalVec, secVelFreestream
    real(dp), allocatable, dimension(:, :) :: secChordwiseResVel, secCP
    real(dp), allocatable, dimension(:) :: secAlpha, secPhi, secCL, secCLu
    real(dp), allocatable, dimension(:) :: secCD, secCM
    real(dp), allocatable, dimension(:) :: alpha0
    integer :: spanwiseLiftSwitch
  contains
    procedure :: move => blade_move
    procedure :: rot_pitch => blade_rot_pitch
    procedure :: rot_axis => blade_rot_axis
    procedure :: rot_wake_axis => blade_rot_wake_axis
    procedure :: rot_pts => blade_rot_pts
    procedure :: vind_bywing => blade_vind_bywing
    procedure :: vindSource_bywing => blade_vindSource_bywing
    procedure :: vind_bywing_boundVortices => blade_vind_bywing_boundVortices
    procedure :: vind_bywing_chordwiseVortices => &
      & blade_vind_bywing_chordwiseVortices
    procedure :: vind_boundVortex => blade_vind_boundVortex
    procedure :: vind_bywake => blade_vind_bywake
    procedure :: convectwake => blade_convectwake
    procedure :: wake_continuity => blade_wake_continuity
    procedure :: getSecDynamicPressure => blade_getSecDynamicPressure
    procedure :: calc_secArea, calc_secChord
    procedure :: calc_force_gamma => blade_calc_force_gamma
    procedure :: calc_force_alpha => blade_calc_force_alpha
    procedure :: calc_force_alphaGamma => blade_calc_force_alphaGamma
    procedure :: calc_secAlpha => blade_calc_secAlpha
    procedure :: calc_secChordwiseResVel => blade_calc_secChordwiseResVel
    procedure :: burst_wake => blade_burst_wake
    procedure :: calc_skew => blade_calc_skew
    procedure :: getSecChordwiseLocations => blade_getSecChordwiseLocations
    procedure :: lookup_secCoeffs => blade_lookup_secCoeffs
    procedure :: secCoeffsToSecForces => blade_secCoeffsTosecForces
    procedure :: dirLiftDrag => blade_dirLiftDrag
    procedure :: sumSecToNetForces => blade_sumSecToNetForces
    procedure :: calc_stlStats => blade_calc_stlStats
    ! I/O subroutines
    procedure :: blade_write
    generic :: write(unformatted) => blade_write
    procedure :: blade_read
    generic :: read(unformatted) => blade_read
  end type blade_class

  type rotor_class
    integer :: nb, ns, nc, nNwake, nFwake, nbConvect, nNwakeEnd, nFwakeEnd
    type(blade_class), allocatable, dimension(:) :: blade
    real(dp) :: Omega, omegaSlow
    real(dp), dimension(3) :: shaftAxis
    real(dp), dimension(3) :: hubCoords, cgCoords, fromCoords
    real(dp) :: radius, chord, root_cut, coningAngle
    real(dp), dimension(3) :: forceInertial, lift, drag
    real(dp), dimension(3) :: dragInduced, dragProfile
    real(dp), dimension(3) :: liftUnsteady, dragUnsteady
    real(dp), dimension(3) :: liftUnitVec, dragUnitVec, sideUnitVec
    real(dp), dimension(3) :: controlPitch  ! theta0,thetaC,thetaS
    real(dp) :: thetaTwist
    real(dp) :: pivotLE  ! pivot location from LE [x/c]
    real(dp) :: flapHinge  ! hinge location from centre [x/R]
    real(dp), dimension(3) :: velBody, omegaBody
    real(dp), allocatable, dimension(:, :) :: velBodyHistory
    real(dp), allocatable, dimension(:, :) :: omegaBodyHistory
    real(dp) :: psi
    real(dp), dimension(3) :: pts  ! phi,theta,psi about cgCoords
    character(len=1) :: streamwiseCoreSwitch
    real(dp) :: spanwiseCore
    real(dp), allocatable, dimension(:) :: streamwiseCoreVec
    real(dp), allocatable, dimension(:, :) :: AIC, AIC_inv
    real(dp), allocatable, dimension(:) :: gamVec, gamVecPrev, RHS
    real(dp), allocatable, dimension(:) :: camberSectionLimit
    real(dp), allocatable, dimension(:) :: airfoilSectionLimit
    real(dp), allocatable, dimension(:) :: alpha0
    real(dp) :: initWakeVel, psiStart, skewLimit
    real(dp) :: apparentViscCoeff, decayCoeff
    real(dp) :: rollupStartRadius, rollupEndRadius
    integer :: propConvention
    integer :: symmetricTau
    integer :: wakeTruncateNt
    integer :: prescWakeNt, prescWakeAfterTruncNt, prescWakeGenNt
    integer :: rollupStart, rollupEnd
    integer :: suppressFwakeSwitch
    integer :: forceCalcSwitch, skewPlotSwitch
    integer :: inflowPlotSwitch
    integer :: spanwiseLiftSwitch, customTrajectorySwitch
    integer :: gammaPlotSwitch
    integer :: rowNear, rowFar
    integer :: nCamberFiles, nAirfoils
    integer :: imagePlane, imageRotorNum
    integer :: surfaceType  
    integer :: imposeAxisymmetry
    character(len=30), allocatable, dimension(:) :: camberFile, airfoilFile
    character(len=30) :: geometryFile
    real(dp) :: nonDimforceDenominator
  contains
    procedure :: read_geom => rotor_read_geom
    procedure :: init => rotor_init
    procedure :: deinit => rotor_deinit
    procedure :: plot3dtoblade => rotor_plot3dtoblade
    procedure :: stltoblade => rotor_stltoblade
    procedure :: getCamber => rotor_getCamber
    procedure :: gettheta => rotor_gettheta
    procedure :: getthetadot => rotor_getthetadot
    procedure :: move => rotor_move
    procedure :: rot_pts => rotor_rot_pts
    procedure :: rot_advance => rotor_rot_advance
    procedure :: assignshed => rotor_assignshed
    procedure :: map_gam => rotor_map_gam
    procedure :: age_wake => rotor_age_wake
    procedure :: dissipate_wake => rotor_dissipate_wake
    procedure :: strain_wake => rotor_strain_wake
    procedure :: calcAIC => rotor_calcAIC
    procedure :: vind_bywing => rotor_vind_bywing
    procedure :: vind_bywing_boundVortices => rotor_vind_bywing_boundVortices
    procedure :: vind_bywake => rotor_vind_bywake
    procedure :: shiftwake => rotor_shiftwake
    procedure :: rollup => rotor_rollup
    procedure :: calc_force_gamma => rotor_calc_force_gamma
    procedure :: calc_force_alpha => rotor_calc_force_alpha
    procedure :: calc_force_alphaGamma => rotor_calc_force_alphaGamma
    procedure :: calc_secAlpha => rotor_calc_secAlpha
    procedure :: convectwake => rotor_convectwake
    procedure :: burst_wake => rotor_burst_wake
    procedure :: calc_skew => rotor_calc_skew
    procedure :: dirLiftDrag => rotor_dirLiftDrag
    procedure :: sumBladeToNetForces => rotor_sumBladeToNetForces
    procedure :: mirrorGeometry => rotor_mirrorGeometry
    procedure :: mirrorGamma => rotor_mirrorGamma
    procedure :: mirrorVelCP => rotor_mirrorVelCP
    procedure :: mirrorWake => rotor_mirrorWake
    procedure :: toChordsRevs => rotor_toChordsRevs
    procedure :: eraseNwake => rotor_eraseNwake
    procedure :: eraseFwake => rotor_eraseFwake
    procedure :: updatePrescribedWake => rotor_updatePrescribedWake
    ! I/O subroutines
    procedure :: rotor_write
    generic :: write(unformatted) => rotor_write
    procedure :: rotor_read
    generic :: read(unformatted) => rotor_read
  end type rotor_class

contains

  !------+--
  ! ++++ | vf_class Methods
  !------+--

  ! Efficient implementation to vind calculation
  function vf_vind(this, P) result(vind)
    ! Compute induced velocity by unit strength vortex filament
    use libMath, only: unitVec
  class(vf_class) :: this
    real(dp), intent(in), dimension(3) :: P
    real(dp), dimension(3) :: vind
    real(dp), dimension(3) :: r1, r2, r0, r1Xr2
    real(dp) :: r1Xr2Abs2

    r1 = P - this%fc(:, 1)
    r2 = P - this%fc(:, 2)
    r0 = r1 - r2

    ! Cross product (inlined to avoid function call)
    r1Xr2(1) = r1(2)*r2(3) - r1(3)*r2(2)
    r1Xr2(2) = r1(3)*r2(1) - r1(1)*r2(3)
    r1Xr2(3) = r1(1)*r2(2) - r1(2)*r2(1)
    r1Xr2Abs2 = dot_product(r1Xr2, r1Xr2)

    vind = 0.

    if (r1Xr2Abs2 > eps*eps) then
      ! Vatistas core model
      vind = (r1Xr2*inv4pi*dot_product(r0, unitVec(r1) - unitVec(r2))) &
        /sqrt((this%rVc*norm2(r0))**4._dp + r1Xr2Abs2**2._dp)
    endif
  end function vf_vind

  subroutine vf_calclength(this, isOriginal)
    ! Compute length of vortex filament
  class(vf_class) :: this
    logical, intent(in) :: isOriginal
    real(dp), dimension(3) :: delta

    delta = this%fc(:, 1) - this%fc(:, 2)
    this%lc = norm2(delta)
    if (isOriginal .eqv. .TRUE.) this%l0 = norm2(delta)
  end subroutine vf_calclength

  subroutine vf_strain(this)
    ! Changes core radius according to change in vortex length
  class(vf_class) :: this
    this%rVc = this%rVc0*sqrt(this%l0/this%lc)
  end subroutine vf_strain

  !------+--
  ! ++++ | vr_class Methods
  !------+--

  function vr_vind(this, P) result(vind)
    ! Compute induced velocity by
    ! unit strength 4-element vortex ring
  class(vr_class) :: this
    real(dp), dimension(3) :: P, vind
    real(dp), dimension(4, 3) :: vindMat
    integer :: i

    vind = 0._dp

    do i = 1, 4
      vindMat(i, :) = this%vf(i)%vind(P)
    enddo
    vind = sum(vindMat, 1)

  end function vr_vind

  function vr_vindSource(this, P, nCap) result(vind)
    ! Compute induced velocity by
    ! unit strength 3-element source ring
  class(vr_class) :: this
    real(dp), intent(in), dimension(3) :: P, nCap
    real(dp), dimension(3) :: vind
    integer :: i

    vind = 0._dp

    ! DEBUG
    ! Add velocity induced by 3d source triangle here

  end function vr_vindSource

  ! Panel coordinates
  ! o---------> Y along span
  ! |
  ! |   1-----------4
  ! |   |     4     |
  ! |   |           |
  ! |   |1         3|
  ! |   |           |
  ! |   |     2     |
  ! |   2-----------3
  ! |
  ! V X along chord

  subroutine vr_assignP(this, n, P)
    ! Assign coordinates to nth corner
  class(vr_class) :: this
    integer, intent(in) :: n
    real(dp), dimension(3) :: P

    select case (n)
    case (1)
      this%vf(4)%fc(:, 2) = P
      this%vf(1)%fc(:, 1) = P
    case (2)
      this%vf(1)%fc(:, 2) = P
      this%vf(2)%fc(:, 1) = P
    case (3)
      this%vf(2)%fc(:, 2) = P
      this%vf(3)%fc(:, 1) = P
    case (4)
      this%vf(3)%fc(:, 2) = P
      this%vf(4)%fc(:, 1) = P
    case default
      error stop 'n may only take values 1,2,3 or 4'
    end select

  end subroutine vr_assignP

  subroutine vr_shiftdP(this, n, dshift)
    ! Shift coordinates of nth corner by dshift distance
    ! (usually for U*dt convection)
  class(vr_class) :: this
    integer, intent(in) :: n
    real(dp), intent(in), dimension(3) :: dshift
    integer :: ifil

    select case (n)
    case (0)
      do ifil = 1, 4
        this%vf(ifil)%fc(:, 1) = this%vf(ifil)%fc(:, 1) + dshift
        this%vf(ifil)%fc(:, 2) = this%vf(ifil)%fc(:, 2) + dshift
      enddo
    case (1)
      this%vf(4)%fc(:, 2) = this%vf(4)%fc(:, 2) + dshift
      this%vf(1)%fc(:, 1) = this%vf(1)%fc(:, 1) + dshift
    case (2)
      this%vf(1)%fc(:, 2) = this%vf(1)%fc(:, 2) + dshift
      this%vf(2)%fc(:, 1) = this%vf(2)%fc(:, 1) + dshift
    case (3)
      this%vf(2)%fc(:, 2) = this%vf(2)%fc(:, 2) + dshift
      this%vf(3)%fc(:, 1) = this%vf(3)%fc(:, 1) + dshift
    case (4)
      this%vf(3)%fc(:, 2) = this%vf(3)%fc(:, 2) + dshift
      this%vf(4)%fc(:, 1) = this%vf(4)%fc(:, 1) + dshift
    case default
      error stop 'n may only take values 1,2,3 or 4'
    end select

  end subroutine vr_shiftdP

  subroutine vr_rot(this, Tmat, originVec)
    !! Rotate vortex ring using Tmat about origin
  class(vr_class) :: this
    real(dp), intent(in), dimension(3, 3) :: Tmat
    real(dp), dimension(3), optional :: originVec
    real(dp), dimension(3) :: origin
    integer :: i

    origin = [0._dp, 0._dp, 0._dp]
    if (present(originVec)) origin = originVec

    do i = 1, 4
      this%vf(i)%fc(:, 1) = matmul(Tmat, this%vf(i)%fc(:, 1)-origin)+origin
      this%vf(i)%fc(:, 2) = matmul(Tmat, this%vf(i)%fc(:, 2)-origin)+origin
    enddo

  end subroutine vr_rot

  subroutine vr_calclength(this, isOriginal)
    ! Calculate length of filaments in vortex ring
  class(vr_class) :: this
    logical, intent(in) :: isOriginal
    integer :: i
    do i = 1, 4
      call this%vf(i)%calclength(isOriginal)
    enddo
  end subroutine vr_calclength

  subroutine vr_strain(this)
  class(vr_class) :: this
    integer :: i
    do i = 1, 4
      call this%vf(i)%strain()
    enddo
  end subroutine vr_strain

  subroutine vr_decay(this, dt, decayCoeff)
  class(vr_class), intent(inout) :: this
    real(dp), intent(in) :: dt, decayCoeff

  this%gam = this%gam*exp(-decayCoeff*dt)

  end subroutine vr_decay

  function vr_getInteriorAngles(this)
    ! Obtain interior angles of vortex ring
    use libMath, only: getAngleCos
  class(vr_class) :: this
    real(dp), dimension(4) :: vr_getInteriorAngles
    real(dp), dimension(3) :: p1, p2, p3, p4

    p1 = this%vf(1)%fc(:, 1)
    p2 = this%vf(2)%fc(:, 1)
    p3 = this%vf(3)%fc(:, 1)
    p4 = this%vf(4)%fc(:, 1)

    vr_getInteriorAngles = 0._dp
    vr_getInteriorAngles(1) = getAngleCos(p2 - p1, p4 - p1)
    vr_getInteriorAngles(2) = getAngleCos(p3 - p2, p1 - p2)
    vr_getInteriorAngles(3) = getAngleCos(p4 - p3, p2 - p3)
    vr_getInteriorAngles(4) = getAngleCos(p3 - p4, p1 - p4)
  end function vr_getInteriorAngles

  function vr_getMedianAngle(this)
    ! Obtain median angle of vortex ring
    use libMath, only: getAngleCos
  class(vr_class) :: this
    real(dp) :: vr_getMedianAngle
    real(dp), dimension(3) :: p1, p2, p3, p4

    p1 = this%vf(1)%fc(:, 1)
    p2 = this%vf(2)%fc(:, 1)
    p3 = this%vf(3)%fc(:, 1)
    p4 = this%vf(4)%fc(:, 1)

    vr_getMedianAngle = getAngleCos(p3 + p4 - p1 - p2, p4 + p1 - p2 - p3)
  end function vr_getMedianAngle

  function vr_getBimedianCos(this)
  class(vr_class) :: this
    real(dp) :: vr_getBimedianCos
    real(dp), dimension(3) :: p1, p2, p3, p4
    real(dp), dimension(3) :: x1Vec, x2Vec

    p1 = this%vf(1)%fc(:, 1)
    p2 = this%vf(2)%fc(:, 1)
    p3 = this%vf(3)%fc(:, 1)
    p4 = this%vf(4)%fc(:, 1)

    x1Vec = p3 + p4 - p1 - p2
    x2Vec = p4 + p1 - p2 - p3
    vr_getBimedianCos = abs(dot_product(x1Vec, x2Vec)/ &
      sqrt(dot_product(x1Vec, x1Vec)*dot_product(x2Vec, x2Vec)))
  end function vr_getBimedianCos

  subroutine vr_mirror(this, coordNum)
    !! Mirror gamma and coordinates about a specified plane
  class(vr_class) :: this
    integer, intent(in) :: coordNum
    integer :: ifil, i

    do ifil = 1, 4
      do i = 1, 2
        this%vf(ifil)%fc(coordNum, i) = &
          & -1._dp * this%vf(ifil)%fc(coordNum, i)
      enddo
    enddo
    this%gam = -1._dp * this%gam
  end subroutine vr_mirror

  subroutine calc_skew(this)
    ! Compute skew
  class(vr_class) :: this

    if (abs(this%gam) > eps) then
      ! skew:  0-good, 1-bad
      this%skew = this%getBimedianCos()
    else
      this%skew = 0._dp
    endif
  end subroutine calc_skew

  subroutine vr_burst(this, skewLimit)
    ! Burst vortex filaments if skewLimit is hit
  class(vr_class) :: this
    real(dp), intent(in) :: skewLimit
    real(dp) :: skewVal

    if ((abs(this%gam) > eps) .and. (skewLimit > eps)) then
      ! skew:  0-good, 1-bad
      skewVal = this%getBimedianCos()
      if (skewVal .ge. skewLimit) this%gam = 0._dp
    endif
    this%skew = skewVal

  end subroutine vr_burst

  !------+--
  ! ++++ | wingpanel_class Methods
  !------+--

  subroutine wingpanel_assignP(this, n, P)
    ! Assign coordinates to nth corner
  class(wingpanel_class) :: this
    integer, intent(in) :: n
    real(dp), dimension(3) :: P

    if (n > 0 .and. n < 5) then
      this%pc(:, n) = P
    else
      error stop 'n may only take values 1,2,3 or 4'
    endif

  end subroutine wingpanel_assignP

  subroutine wingpanel_calcCP(this, isTriangle)
    ! Compute collocation point location
  class(wingpanel_class) :: this
    logical, optional :: isTriangle
    logical :: triPanel

    triPanel = .false.  ! Assumes quad panel by default
    if (present(isTriangle)) triPanel = isTriangle

    if (triPanel) then
      this%CP = (this%PC(:, 1) + this%PC(:, 2) + this%PC(:, 3))/3.0_dp
    else
      this%CP = ((this%PC(:, 1) + this%PC(:, 4))*0.25_dp &
        + (this%PC(:, 2) + this%PC(:, 3))*0.75_dp)*0.5_dp
    endif
  end subroutine wingpanel_calcCP

  subroutine wingpanel_calcN(this, isTriangle)
    ! Compute normal vector
    use libMath, only: unitVec, cross_product
  class(wingpanel_class) :: this
    logical, optional :: isTriangle
    logical :: triPanel

    triPanel = .false.  ! Assumes quad panel by default
    if (present(isTriangle)) triPanel = isTriangle

    if (triPanel) then
      this%nCap = unitvec(this%nCap)
    else
      this%nCap = unitVec(cross_product(this%pc(:, 3) - this%pc(:, 1), &
        this%pc(:, 4) - this%pc(:, 2)))
    endif
  end subroutine wingpanel_calcN

  subroutine wingpanel_calcTau(this, isTriangle)
    ! Compute chordwise and spanwise tangential vectors
    use libMath, only: unitVec
  class(wingpanel_class) :: this
    logical, optional :: isTriangle
    logical :: triPanel

    triPanel = .false.  ! Assumes quad panel by default
    if (present(isTriangle)) triPanel = isTriangle

    if (triPanel) then
      this%tauCapChord = unitvec(this%PC(:, 3)-this%PC(:, 1))
      this%tauCapSpan = unitvec(this%PC(:, 4)-this%PC(:, 2))
    else
      this%tauCapChord = unitVec(0.5_dp*((this%pc(:, 2) + this%pc(:, 3)) &
        - (this%pc(:, 1) + this%pc(:, 4))))
      this%tauCapSpan = unitVec(0.5_dp*((this%pc(:, 3) + this%pc(:, 4)) &
        - (this%pc(:, 2) + this%pc(:, 1))))
    endif
  end subroutine wingpanel_calcTau

  subroutine wingpanel_rot(this, Tmat, originVec)
    ! Rotate panel using transformation matrix
  class(wingpanel_class) :: this
    real(dp), dimension(3, 3) :: Tmat
    real(dp), dimension(3), optional :: originVec
    real(dp), dimension(3) :: origin
    integer :: i

    origin = [0._dp, 0._dp, 0._dp]
    if (present(originVec)) origin = originVec

    do i = 1, 4
      this%pc(:, i) = matmul(Tmat, this%pc(:, i)-origin)+origin
    enddo
    call this%vr%rot(Tmat, origin)
    this%CP = matmul(Tmat, this%CP-origin)+origin
    this%nCap = matmul(Tmat, this%nCap)
    this%tauCapChord = matmul(Tmat, this%tauCapChord)
    this%tauCapSpan = matmul(Tmat, this%tauCapSpan)
  end subroutine wingpanel_rot

  subroutine wingpanel_shiftdP(this, dshift)
    ! Shift corners of vortex ring by dshift
  class(wingpanel_class) :: this
    real(dp), intent(in), dimension(3) :: dshift
    integer :: i

    this%CP = this%CP + dshift
    do i = 1, 4
      this%pc(:, i) = this%pc(:, i) + dshift
      call this%vr%shiftdP(i, dshift)
    enddo

  end subroutine wingpanel_shiftdP

  subroutine wingpanel_calc_area(this)
    use libMath, only: cross_product
  class(wingpanel_class) :: this
    this%panelArea = 0.5_dp*norm2(cross_product(this%pc(:, 3) &
      - this%pc(:, 1), this%pc(:, 4) - this%pc(:, 2)))
  end subroutine wingpanel_calc_area

  subroutine wingpanel_calc_mean_dimensions(this)
    ! Calculate mean chord and mean span
  class(wingpanel_class) :: this
    this%meanSpan = 0.5_dp*(norm2(this%pc(:, 4) - this%pc(:, 1)) &
      + norm2(this%pc(:, 3) - this%pc(:, 2)))
    this%meanChord = 0.5_dp*(norm2(this%pc(:, 2) - this%pc(:, 1)) &
      + norm2(this%pc(:, 3) - this%pc(:, 4)))
  end subroutine wingpanel_calc_mean_dimensions

  function wingpanel_isCPinsidecore(this)
    ! Check whether collocation point lies
    ! inside viscous core region of vortex ring
  class(wingpanel_class), intent(in) :: this
    logical :: wingpanel_isCPinsidecore
    real(dp) :: deltaxby4, deltayby2

    deltaxby4 = 0.25_dp*abs(this%vr%vf(1)%fc(1, 1) - this%vr%vf(2)%fc(1, 1))
    deltayby2 = 0.5_dp*abs(this%vr%vf(1)%fc(2, 1) - this%vr%vf(4)%fc(2, 1))

    wingpanel_isCPinsidecore = .false.
    if (deltayby2 .lt. this%vr%vf(1)%rVc) then
      wingpanel_isCPinsidecore = .true.  ! Left edge
    elseif (deltayby2 .lt. this%vr%vf(3)%rVc) then
      wingpanel_isCPinsidecore = .true.  ! Right edge
    elseif (deltaxby4 .lt. this%vr%vf(2)%rVc) then
      wingpanel_isCPinsidecore = .true.  ! Upper edge
    elseif (3._dp*deltaxby4 .lt. this%vr%vf(4)%rVc) then
      wingpanel_isCPinsidecore = .true.  ! Bottom edge
    endif
  end function wingpanel_isCPinsidecore

  subroutine wingpanel_calc_chordwiseResVel(this)
    ! Compute panel resultant velocities using local velocities
  class(wingpanel_class), intent(inout) :: this

    this%chordwiseResVel = this%velCPTotal - &
      dot_product(this%velCPTotal, this%tauCapSpan)*this%tauCapSpan
  end subroutine wingpanel_calc_chordwiseResVel

  !subroutine wingpanel_calc_alpha(this)
  !  ! Compute panel alpha using local velocities
  !class(wingpanel_class), intent(inout) :: this
  !  real(dp) :: velCPTotalChordwiseProjectedMagnitude
  !  real(dp), dimension(3) :: velCPTotalChordwiseProjected

  !  velCPTotalChordwiseProjected=this%velCPTotal- &
  !    dot_product(this%velCPTotal,this%tauCapSpan)*this%tauCapSpan

  !  velCPTotalChordwiseProjectedMagnitude=norm2(velCPTotalChordwiseProjected)

  !  if (velCPTotalChordwiseProjectedMagnitude .gt. eps) then
  !    this%alpha=acos(dot_product(velCPTotalChordwiseProjected, &
  !      & this%tauCapChord) &
  !      & /velCPTotalChordwiseProjectedMagnitude)
  !  else
  !    this%alpha=0._dp
  !  endif
  !end subroutine wingpanel_calc_alpha

  !------+--
  ! ++++ | Nwake_class Methods
  !------+--

  !------+--
  ! ++++ | Fwake_class Methods
  !------+--

  subroutine Fwake_shiftdP(this, n, dshift)
    ! Shift coordinates of nth corner by dshift distance
    ! (usually for Udt convection)
  class(Fwake_class) :: this
    integer, intent(in) :: n
    real(dp), intent(in), dimension(3) :: dshift

    if (n < 0 .or. n > 2) &
      error stop 'n may only take values 0, 1 or 2 in Fwake_shiftdP()'
    if (n == 0) then
      this%vf%fc(:, 1) = this%vf%fc(:, 1) + dshift
      this%vf%fc(:, 2) = this%vf%fc(:, 2) + dshift
    else
      this%vf%fc(:, n) = this%vf%fc(:, n) + dshift
    endif
  end subroutine Fwake_shiftdP

  subroutine Fwake_assignP(this, n, P)
    ! Assign point to nth endpoint of filament
  class(Fwake_class) :: this
    integer, intent(in) :: n
    real(dp), intent(in), dimension(3) :: P

    if (n /= 1 .and. n /= 2) &
      error stop 'n may only take values 1 or 2 in Fwake_assignP()'
    this%vf%fc(:, n) = P
  end subroutine Fwake_assignP

  subroutine Fwake_rot(this, TMat, originVec)
    !! Rotate using TMat about originVec
  class(Fwake_class) :: this
    real(dp), intent(in), dimension(3, 3) :: TMat
    real(dp), dimension(3), optional :: originVec
    real(dp), dimension(3) :: origin

    origin = [0._dp, 0._dp, 0._dp]
    if (present(originVec)) origin = originVec

    this%vf%fc(:, 1) = matmul(TMat, this%vf%fc(:, 1)-origin)+origin
    this%vf%fc(:, 2) = matmul(TMat, this%vf%fc(:, 2)-origin)+origin
  end subroutine Fwake_rot

  subroutine Fwake_decay(this, dt, decayCoeff)
  class(Fwake_class), intent(inout) :: this
    real(dp), intent(in) :: dt, decayCoeff

    this%gam = this%gam*exp(-decayCoeff*dt)
    end subroutine Fwake_decay

  subroutine Fwake_mirror(this, coordNum)
    !! Mirror gamma and coordinates about a specified plane
  class(Fwake_class) :: this
    integer, intent(in) :: coordNum
    integer :: i

    do i = 1, 2
      this%vf%fc(coordNum, i) = &
        & -1._dp * this%vf%fc(coordNum, i)
    enddo
    this%gam = -1._dp * this%gam
  end subroutine Fwake_mirror

  !------+--
  ! ++++ | pFwake_class Methods
  !------+--
  subroutine pFwake_update(this, waF, hubCoords, shaftAxis, deltaPsi)
    use libMath, only: linspace
  class(pFwake_class), intent(inout) :: this
    type(Fwake_class), intent(in), dimension(:) :: waF
    real(dp), intent(in), dimension(3) :: hubCoords, shaftAxis
    real(dp), intent(in) :: deltaPsi
    real(dp), dimension(size(this%coords, 2)) :: theta
    real(dp), dimension(3) :: anchor
    real(dp) :: helixRadiusCurrent, helixPitchCurrent, deltaZ, dTheta
    integer :: i, npFwake, nFwake

    if (abs(shaftAxis(1)) > eps .or. abs(shaftAxis(2)) > eps) then
      error stop "Prescribed far wake only implemented for shaft along Z-axis"
    endif

    this%isPresent = .true.

    ! Find helix parameters
    nFwake = size(waF)
    anchor = waF(nFwake)%vf%fc(:, 1)

    ! Radius and pitch of helix computed using 
    ! average radius and slope of all far wake filaments
    helixPitchCurrent = 0._dp
    helixRadiusCurrent = 0._dp
    do i = 1, nFwake
      helixRadiusCurrent = helixRadiusCurrent + &
        & norm2((/waF(i)%vf%fc(2, 1), waF(i)%vf%fc(1, 1)/))
      if (i < size(waF)) then
        helixPitchCurrent = helixPitchCurrent + &
          & waF(i)%vf%fc(3, 1) - waF(i+1)%vf%fc(3, 1)
      endif
    enddo
    helixPitchCurrent = abs(helixPitchCurrent) * (-twoPi/deltaPsi)/(nFwake-1)
    helixRadiusCurrent = helixRadiusCurrent / nFwake

    ! Update pitch and radius using a relaxation factor
    ! to avoid sudden variations
    this%helixPitch = this%relaxFactor*helixPitchCurrent + &
      & (1-this%relaxFactor)*this%helixPitch
    this%helixRadius = this%relaxFactor*helixRadiusCurrent + &
      & (1-this%relaxFactor)*this%helixRadius

    ! Angle by which unit helix has to be rotated
    dTheta = atan2(anchor(2), anchor(1))

    ! delta z by which unit helix has to be translated
    deltaZ = anchor(3)-hubCoords(3)

    theta = linspace(0._dp, twoPi*this%nRevs, size(theta, 1))
    if (this%isClockwiseRotor) theta = -1._dp*theta

    this%coords(1, :) = this%helixRadius * cos(theta+dTheta)
    this%coords(2, :) = this%helixRadius * sin(theta+dTheta)
    this%coords(3, :) = this%helixPitch * abs(theta)/twoPi + deltaZ

    ! Assign to prescribed wake
    npFwake = size(this%waF)
    do i = 1, npFwake
      call this%waF(i)%assignP(2, hubCoords + this%coords(:, i))
      call this%waF(i)%assignP(1, hubCoords + this%coords(:, i+1))
    enddo

    ! To maintain continuity
    call this%waF(1)%assignP(2, anchor)

    this%waF%gam = waF(nFwake)%gam
    this%waF%vf%rVc = waF(nFwake)%vf%rVc
  end subroutine pFwake_update

  subroutine pFwake_rot_wake_axis(this, theta, axisVec, origin)
    use libMath, only: getTransformAxis
  class(pFwake_class), intent(inout) :: this
    real(dp), intent(in) :: theta
    real(dp), intent(in), dimension(3) :: axisVec
    real(dp), intent(in), dimension(3) :: origin
    real(dp), dimension(3, 3) :: Tmat
    integer :: i

    if (abs(theta) > eps) then
      TMat = getTransformAxis(theta, axisVec)

      !$omp parallel do
      do i = 1, size(this%waF)
        call this%waF(i)%rot(TMat, origin)
      enddo
      !$omp end parallel do
    endif
  end subroutine pFwake_rot_wake_axis

  !------+--
  ! ++++ | blade_class Methods
  !------+--

  subroutine blade_move(this, dshift)
    ! Move blade by dshift
  class(blade_class) :: this
    real(dp), intent(in), dimension(3) :: dshift
    integer :: i, j

    !$omp parallel do collapse(2)
    do j = 1, this%ns
      do i = 1, this%nc
        call this%wiP(i, j)%shiftdP(dshift)
      enddo
    enddo
    !$omp end parallel do

    do i = 1, size(this%secCP, 2)
      this%secCP(:, i) = this%secCP(:, i) + dshift
    enddo

  end subroutine blade_move

  subroutine blade_rot_pts(this, pts, origin, order)
    ! Rotate blade using pts => phi theta psi
    use libMath, only: Tbg, Tgb
  class(blade_class), intent(inout) :: this
    real(dp), dimension(3), intent(in) :: pts    ! pts => phi,theta,psi
    real(dp), dimension(3), intent(in) :: origin ! rotation about
    integer, intent(in) :: order    ! [1]gb & +ve theta , [2]bg & -ve theta
    integer :: i, j
    real(dp), dimension(3, 3) :: TMat

    select case (order)
    case (2)
      TMat = Tbg((/cos(pts(1)), sin(pts(1))/), &
        (/cos(pts(2)), sin(pts(2))/), &
        (/cos(pts(3)), sin(pts(3))/))
    case (1)
      TMat = Tgb((/cos(pts(1)), sin(pts(1))/), &
        (/cos(pts(2)), sin(pts(2))/), &
        (/cos(pts(3)), sin(pts(3))/))
    case default
      error stop 'ERROR: Wrong option for order'
    end select

    !$omp parallel do
    do j = 1, this%ns
      do i = 1, this%nc
        call this%wiP(i, j)%rot(TMat, origin)
      enddo

      this%secCP(:, j) = matmul(TMat, this%secCP(:, j)-origin)+origin

      ! Rotate sec vectors
      this%secTauCapChord(:, j) = matmul(TMat, this%secTauCapChord(:, j))
      this%secNormalVec(:, j) = matmul(TMat, this%secNormalVec(:, j))
    enddo
    !$omp end parallel do

    this%xAxis = matmul(Tmat, this%xAxis)
    this%yAxis = matmul(Tmat, this%yAxis)
    this%zAxis = matmul(Tmat, this%zAxis)

  end subroutine blade_rot_pts

  subroutine blade_rot_pitch(this, theta)
    ! Rotate blade by pitch angle
    ! pivot point calculated using straight line joining
    ! LE of first panel and TE of last panel
  class(blade_class), intent(inout) :: this
    real(dp), intent(in) :: theta
    real(dp), dimension(3) :: axis
    real(dp), dimension(3) :: axisOrigin!, axisEnd

    if (abs(theta) > eps) then
      axisOrigin = this%wiP(1, 1)%PC(:, 1)*(1._dp - this%pivotLE) &
        + this%wiP(this%nc, 1)%PC(:, 2)*this%pivotLE
      !axisEnd=this%wiP(1, this%ns)%PC(:,4)*(1._dp-this%pivotLE) &
      !  +this%wiP(this%nc, this%ns)%PC(:,3)*this%pivotLE
      !
      !! Construct axes of rotation from LE of first panel
      !axis=axisEnd-axisOrigin
      !axis=axis/norm2(axis)

      ! Use blade X axis for rotation
      axis = this%yAxis

      call this%rot_axis(theta, axis, axisOrigin)
    endif
  end subroutine blade_rot_pitch

  subroutine blade_rot_axis(this, theta, axisVec, origin)
    ! Rotate about axis at specified origin
    use libMath, only: getTransformAxis
  class(blade_class), intent(inout) :: this
    real(dp), intent(in) :: theta
    real(dp), intent(in), dimension(3) :: axisVec
    real(dp), intent(in), dimension(3) :: origin
    real(dp), dimension(3, 3) :: Tmat
    integer :: i, j

    if (abs(theta) > eps) then
      ! Translate to origin
      call this%move(-origin)

      ! Rotate about axisVec
      Tmat = getTransformAxis(theta, axisVec)
      do j = 1, this%ns
        do i = 1, this%nc
          call this%wiP(i, j)%rot(TMat)
        enddo
      enddo

      ! Untranslate from origin
      call this%move(origin)

      ! Rotate secCP also
      !$omp parallel do
      do i = 1, this%ns
        this%secCP(:, i) = matmul(TMat, this%secCP(:, i)-origin)+origin

      ! Rotate sec vectors also along with blade
        this%secTauCapChord(:, i) = matmul(TMat, this%secTauCapChord(:, i))
        this%secNormalVec(:, i) = matmul(TMat, this%secNormalVec(:, i))
      enddo
      !$omp end parallel do

      ! Rotate blade axes
      this%xAxis = matmul(TMat, this%xAxis)
      this%yAxis = matmul(TMat, this%yAxis)
      this%zAxis = matmul(TMat, this%zAxis)
    endif
  end subroutine blade_rot_axis

  subroutine blade_rot_wake_axis(this, theta, axisVec, origin, &
      & rowNear, rowFar, wakeType)
    ! Rotate wake about axis at specified origin
    use libMath, only: getTransformAxis
  class(blade_class), intent(inout) :: this
    real(dp), intent(in) :: theta
    real(dp), intent(in), dimension(3) :: axisVec
    real(dp), intent(in), dimension(3) :: origin
    integer, intent(in) :: rowNear, rowFar
    character(len=1), intent(in) :: wakeType  ! For predicted wake
    integer :: nNwake, nFwake
    real(dp), dimension(3, 3) :: Tmat
    integer :: i, j

    if (abs(theta) > eps) then
      TMat = getTransformAxis(theta, axisVec)

      select case (wakeType)
      case ('C')
        ! Start procedure for wake rotation
        nNwake = size(this%waN, 1)
        nFwake = size(this%waF, 1)

        !$omp parallel do collapse(2)
        do j = 1, this%ns
          do i = rowNear, nNwake
            call this%waN(i, j)%vr%rot(TMat, origin)
          enddo
        enddo
        !$omp end parallel do

        !$omp parallel do
        do i = rowFar, nFwake
          call this%waF(i)%rot(TMat, origin)
        enddo
        !$omp end parallel do

      case ('P')
        ! Start procedure for wake rotation
        nNwake = size(this%waNPredicted, 1)
        nFwake = size(this%waFPredicted, 1)

        !$omp parallel do collapse(2)
        do j = 1, this%ns
          do i = rowNear, nNwake
            call this%waNPredicted(i, j)%vr%rot(TMat, origin)
          enddo
        enddo
        !$omp end parallel do

        !$omp parallel do
        do i = rowFar, nFwake
          call this%waFPredicted(i)%rot(TMat, origin)
        enddo
        !$omp end parallel do
      end select
    endif
  end subroutine blade_rot_wake_axis

  function blade_vind_bywing(this, P)
    ! Compute induced velocity by blade bound vorticity
  class(blade_class), intent(inout) :: this
    real(dp), intent(in), dimension(3) :: P
    real(dp), dimension(3) :: blade_vind_bywing
    integer :: i, j

    blade_vind_bywing = 0._dp
    do j = 1, this%ns
      do i = 1, this%nc
        blade_vind_bywing = blade_vind_bywing + &
          & this%wiP(i, j)%vr%vind(P)*this%wiP(i, j)%vr%gam
      enddo
    enddo

  end function blade_vind_bywing

  function blade_vindSource_bywing(this, P)
    ! Compute induced velocity by blade bound vorticity
  class(blade_class), intent(in) :: this
    real(dp), intent(in), dimension(3) :: P
    real(dp), dimension(3) :: blade_vindSource_bywing
    integer :: i, j

    blade_vindSource_bywing = 0._dp
    do j = 1, this%ns
      do i = 1, this%nc
        blade_vindSource_bywing = blade_vindSource_bywing + &
          this%wiP(i, j)%vr%vindSource(P, this%wiP(i,j)%nCap)* &
          & this%wiP(i, j)%vr%gam
      enddo
    enddo
  end function blade_vindSource_bywing

  function blade_vind_bywing_boundVortices(this, P)
    ! Compute induced velocity by bound vortices alone
  class(blade_class), intent(inout) :: this
    real(dp), intent(in), dimension(3) :: P
    real(dp), dimension(3) :: blade_vind_bywing_boundVortices
    integer :: i, j

    blade_vind_bywing_boundVortices = 0._dp
    do j = 1, this%ns
      do i = 1, this%nc
        blade_vind_bywing_boundVortices = blade_vind_bywing_boundVortices + &
          & (this%wiP(i, j)%vr%vf(2)%vind(P) + &
          & this%wiP(i, j)%vr%vf(4)%vind(P))*this%wiP(i, j)%vr%gam
      enddo
    enddo

    do j = 1, this%ns
      blade_vind_bywing_boundVortices = blade_vind_bywing_boundVortices - &
        this%wiP(this%nc, j)%vr%vf(2)%vind(P)*this%wiP(this%nc, j)%vr%gam
    enddo
  end function blade_vind_bywing_boundVortices

  function blade_vind_bywing_chordwiseVortices(this, P)
    ! Compute induced velocity by bound vortices alone
  class(blade_class), intent(in) :: this
    real(dp), intent(in), dimension(3) :: P
    real(dp), dimension(3) :: blade_vind_bywing_chordwiseVortices
    integer :: i, j

    blade_vind_bywing_chordwiseVortices = 0._dp
    do j = 1, this%ns
      do i = 1, this%nc
        blade_vind_bywing_chordwiseVortices = blade_vind_bywing_chordwiseVortices + &
          (this%wiP(i, j)%vr%vf(1)%vind(P) + this%wiP(i, j)%vr%vf(3)%vind(P))* &
          this%wiP(i, j)%vr%gam
      enddo
    enddo
    do j = 1, this%ns
      blade_vind_bywing_chordwiseVortices = &
        & blade_vind_bywing_chordwiseVortices + &
        & this%wiP(this%nc, j)%vr%vf(2)%vind(P)*this%wiP(this%nc, j)%vr%gam
    enddo
  end function blade_vind_bywing_chordwiseVortices

  function blade_vind_boundVortex(this, ic, is, P)
    ! Compute induced velocity by bound vortices alone
  class(blade_class), intent(inout) :: this
    integer, intent(in) :: ic, is
    real(dp), intent(in), dimension(3) :: P
    real(dp), dimension(3) :: blade_vind_boundVortex

    if (ic > 1) then
      blade_vind_boundVortex = this%wiP(ic, is)%vr%vf(4)%vind(P) &
        *this%wiP(ic, is)%vr%gam &
        + this%wiP(ic - 1, is)%vr%vf(2)%vind(P)*this%wiP(ic - 1, is)%vr%gam
    else
      blade_vind_boundVortex = this%wiP(ic, is)%vr%vf(4)%vind(P) &
        *this%wiP(ic, is)%vr%gam
    endif
  end function blade_vind_boundVortex

  function blade_vind_bywake(this, rowNear, rowFar, P, optionalChar)
    ! Compute induced velocity by wake vortex rings
  class(blade_class), intent(in) :: this
    integer, intent(in) :: rowNear, rowFar
    real(dp), intent(in), dimension(3) :: P
    character(len=1), optional :: optionalChar
    real(dp), dimension(3) :: blade_vind_bywake
    integer :: i, j, nNwake, nFwake

    nNwake = size(this%waN, 1)
    nFwake = size(this%waF, 1)
    blade_vind_bywake = 0._dp
    if (.not. present(optionalChar)) then
      do j = 1, size(this%waN, 2)
        do i = rowNear, nNwake
          if (abs(this%waN(i, j)%vr%gam) .gt. eps) &
            blade_vind_bywake = blade_vind_bywake + &
            this%waN(i, j)%vr%vind(P)*this%waN(i, j)%vr%gam
        enddo
      enddo

      if (rowFar .le. nFwake) then
        ! Last row of Nwake is made of horseshoe vortices, if Fwake is generated
        do j = 1, size(this%waN, 2)
          blade_vind_bywake = blade_vind_bywake - &
            this%waN(nNwake, j)%vr%vf(2)%vind(P)*this%waN(nNwake, j)%vr%gam
        enddo

        do i = rowFar, nFwake
          if (abs(this%waF(i)%gam) .gt. eps) &
            blade_vind_bywake = blade_vind_bywake &
            + this%waF(i)%vf%vind(P)*this%waF(i)%gam
        enddo

        do i = 1, size(this%wapF%waF)
          if (abs(this%wapF%waF(i)%gam) .gt. eps) then
            blade_vind_bywake = blade_vind_bywake &
              & + this%wapF%waF(i)%vf%vind(P)*this%wapF%waF(i)%gam
          endif
        enddo
      endif

    elseif ((optionalChar .eq. 'P') .or. (optionalChar .eq. 'p')) then
      do j = 1, size(this%waN, 2)
        do i = rowNear, nNwake
          if (abs(this%waNPredicted(i, j)%vr%gam) .gt. eps) &
            blade_vind_bywake = blade_vind_bywake &
            + this%waNPredicted(i, j)%vr%vind(P)*this%waNPredicted(i, j)%vr%gam
        enddo
      enddo

      if (rowFar .le. nFwake) then
        ! Last row of Nwake is made of horseshoe vortices, if Fwake is generated
        do j = 1, size(this%waN, 2)
          blade_vind_bywake = blade_vind_bywake &
            - this%waNPredicted(nNwake, j)%vr%vf(2)%vind(P)*this%waNPredicted(nNwake, j)%vr%gam
        enddo

        do i = rowFar, nFwake
          if (abs(this%waFPredicted(i)%gam) .gt. eps) &
            blade_vind_bywake = blade_vind_bywake &
            + this%waFPredicted(i)%vf%vind(P)*this%waFPredicted(i)%gam
        enddo

        do i = 1, size(this%wapFPredicted%waF)
          if (abs(this%wapFPredicted%waF(i)%gam) .gt. eps) then
            blade_vind_bywake = blade_vind_bywake &
              & + this%wapFPredicted%waF(i)%vf%vind(P) &
              & * this%wapFPredicted%waF(i)%gam
          endif
        enddo
      endif
    else
      error stop 'ERROR: Wrong character flag for blade_vind_bywake()'
    endif

  end function blade_vind_bywake

subroutine blade_convectwake(this, rowNear, rowFar, dt, wakeType)
  ! Convect wake collocation points using velNwake matrix
class(blade_class), intent(inout) :: this
  integer, intent(in) :: rowNear, rowFar
  real(dp), intent(in) :: dt
  character(len=1), intent(in) :: wakeType  ! For predicted wake
  integer :: i, j, nNwake, nFwake

    nNwake = size(this%waN, 1)

    select case (wakeType)
    case ('C')    ! [C]urrent wake
      !$omp parallel do collapse(2)
      do j = 1, this%ns
        do i = rowNear, nNwake
          call this%waN(i, j)%vr%shiftdP(2, this%velNwake(:, i, j)*dt)
        enddo
      enddo
      !$omp end parallel do

      !$omp parallel do
      do i = rowNear, nNwake
        call this%waN(i, this%ns)%vr%shiftdP(3, this%velNwake(:, i, this%ns + 1)*dt)
      enddo
      !$omp end parallel do

      nFwake = size(this%waF, 1)
      !$omp parallel do
      do i = rowFar, nFwake
        call this%waF(i)%shiftdP(1, this%velFwake(:, i)*dt)  ! Shift only TE
      enddo
      !$omp end parallel do

    case ('P')    ! [P]redicted wake
      !$omp parallel do collapse(2)
      do j = 1, this%ns
        do i = 1, rowNear, nNwake
          call this%waNPredicted(i, j)%vr%shiftdP(2, this%velNwake(:, i, j)*dt)
        enddo
      enddo
      !$omp end parallel do

      !$omp parallel do
      do i = 1, rowNear, nNwake
        call this%waNPredicted(i, this%ns)%vr%shiftdP(3, this%velNwake(:, i, this%ns + 1)*dt)
      enddo
      !$omp end parallel do

      nFwake = size(this%waF, 1)
      !$omp parallel do
      do i = rowFar, nFwake
        call this%waFPredicted(i)%shiftdP(1, this%velFwake(:, i)*dt)  ! Shift only TE
      enddo
      !$omp end parallel do

    end select

    call this%wake_continuity(rowNear, rowFar, wakeType)

  end subroutine blade_convectwake

  subroutine blade_wake_continuity(this, rowNear, rowFar, wakeType)
    ! Maintain continuity between vortex ring elements after convection
    ! of wake collocation points
  class(blade_class), intent(inout) :: this
    integer, intent(in) :: rowNear, rowFar
    character(len=1), intent(in) :: wakeType  ! For predicted wake
    integer :: i, j, nNwake, nFwake

    nNwake = size(this%waN, 1)

    select case (wakeType)
    case ('C')
      !$omp parallel do collapse(2)
      do j = 1, this%ns - 1
        do i = rowNear + 1, nNwake
          call this%waN(i, j)%vr%assignP(1, this%waN(i - 1, j)%vr%vf(2)%fc(:, 1))
          call this%waN(i, j)%vr%assignP(3, this%waN(i, j + 1)%vr%vf(2)%fc(:, 1))
          call this%waN(i, j)%vr%assignP(4, this%waN(i - 1, j + 1)%vr%vf(2)%fc(:, 1))
        enddo
      enddo
      !$omp end parallel do

      !$omp parallel do
      do j = 1, this%ns - 1
        call this%waN(rowNear, j)%vr%assignP(3, this%waN(rowNear, j + 1)%vr%vf(2)%fc(:, 1))
      enddo
      !$omp end parallel do

      !$omp parallel do
      do i = rowNear + 1, nNwake
        call this%waN(i, this%ns)%vr%assignP(1, this%waN(i - 1, this%ns)%vr%vf(2)%fc(:, 1))
        call this%waN(i, this%ns)%vr%assignP(4, this%waN(i - 1, this%ns)%vr%vf(3)%fc(:, 1))
      enddo
      !$omp end parallel do

      nFwake = size(this%waF, 1)
      !$omp parallel do
      do i = rowFar + 1, nFwake
        call this%waF(i)%assignP(2, this%waF(i - 1)%vf%fc(:, 1))
      enddo
      !$omp end parallel do

    case ('P')
      ! For predicted wake

      !$omp parallel do collapse(2)
      do j = 1, this%ns - 1
        do i = rowNear + 1, nNwake
          call this%waNPredicted(i, j)%vr%assignP(1, this%waNPredicted(i - 1, j)%vr%vf(2)%fc(:, 1))
          call this%waNPredicted(i, j)%vr%assignP(3, this%waNPredicted(i, j + 1)%vr%vf(2)%fc(:, 1))
          call this%waNPredicted(i, j)%vr%assignP(4, this%waNPredicted(i - 1, j + 1)%vr%vf(2)%fc(:, 1))
        enddo
      enddo
      !$omp end parallel do

      !$omp parallel do
      do j = 1, this%ns - 1
        call this%waNPredicted(rowNear, j)%vr%assignP(3, this%waNPredicted(rowNear, j + 1)%vr%vf(2)%fc(:, 1))
      enddo
      !$omp end parallel do

      !$omp parallel do
      do i = rowNear + 1, nNwake
        call this%waNPredicted(i, this%ns)%vr%assignP(1, this%waNPredicted(i - 1, this%ns)%vr%vf(2)%fc(:, 1))
        call this%waNPredicted(i, this%ns)%vr%assignP(4, this%waNPredicted(i - 1, this%ns)%vr%vf(3)%fc(:, 1))
      enddo
      !$omp end parallel do

      nFwake = size(this%waF, 1)
      !$omp parallel do
      do i = rowFar + 1, nFwake
        call this%waFPredicted(i)%assignP(2, this%waFPredicted(i - 1)%vf%fc(:, 1))
      enddo
      !$omp end parallel do

    case default
      error stop 'ERROR: Wrong character flag for convectwake()'
    end select

  end subroutine blade_wake_continuity

  subroutine blade_calc_force_gamma(this, density, invertGammaSign, dt)
    ! Compute force using blade circulation
    use libMath, only: unitVec, cross_product, projVec
  class(blade_class), intent(inout) :: this
    real(dp), intent(in) :: density, invertGammaSign, dt
    integer :: is, ic
    real(dp), dimension(this%nc, this%ns) :: velTangentialChord
    real(dp), dimension(this%nc, this%ns) :: velTangentialSpan
    real(dp), dimension(this%nc, this%ns) :: velInduced
    real(dp), dimension(this%nc, this%ns) :: gamElementChord, gamElementSpan
    real(dp), dimension(this%ns) :: secDynamicPressure
    real(dp) :: signSecCL

    this%forceInertial = 0._dp
    this%secForceInertial = 0._dp
    this%secLift = 0._dp
    this%secDrag = 0._dp
    this%secLiftUnsteady = 0._dp

    ! Compute tangential velocity
    do is = 1, this%ns
      do ic = 1, this%nc
        velTangentialChord(ic, is) = dot_product(this%wiP(ic, is)%velCP, this%wiP(ic, is)%tauCapChord)
        velTangentialSpan(ic, is) = dot_product(this%wiP(ic, is)%velCP, this%wiP(ic, is)%tauCapSpan)
        velInduced(ic, is) = dot_product(this%wiP(ic, is)%velCP + &
          this%vind_bywing_chordwiseVortices(this%wiP(ic, is)%CP), &
          unitVec(cross_product(this%wiP(ic, is)%velCPm,this%yAxis)))
      enddo
    enddo

    ! Compute chordwise elemental circulation of edge panels
    do is = 1, this%ns
      gamElementChord(1, is) = this%wiP(1, is)%vr%gam
    enddo
    do ic = 2, this%nc
      gamElementChord(ic, 1) = this%wiP(ic, 1)%vr%gam &
        & -this%wiP(ic - 1, 1)%vr%gam
    enddo

    ! Compute spanwise elemental circulation of edge panels
    do ic = 1, this%nc
      gamElementSpan(ic, 1) = this%wiP(ic, 1)%vr%gam
    enddo
    do is = 2, this%ns
      gamElementSpan(1, is) = this%wiP(1, is)%vr%gam &
        & -this%wiP(1, is - 1)%vr%gam
    enddo

    ! Compute chordwise and spanwise elemental circulations of inner panels
    do is = 2, this%ns
      do ic = 2, this%nc
        gamElementChord(ic, is) = this%wiP(ic, is)%vr%gam - this%wiP(ic - 1, is)%vr%gam
        gamElementSpan(ic, is) = this%wiP(ic, is)%vr%gam - this%wiP(ic, is - 1)%vr%gam
      enddo
    enddo

    ! Invert gamma sign for correct computation
    gamElementSpan = -1._dp*gamElementSpan
    gamElementChord = -1._dp*gamElementChord

    ! Compute delP
    do is = 1, this%ns
      do ic = 1, this%nc
        ! Use trapezoidal rule on two points to get current gam
        ! for computing unsteady lift part
        if (ic > 1) then
          this%wiP(ic, is)%gamTrapz = -0.5_dp*(this%wiP(ic, is)%vr%gam + this%wiP(ic - 1, is)%vr%gam)
        else
          this%wiP(1, is)%gamTrapz = -0.5_dp*this%wiP(1, is)%vr%gam
        endif

        ! For checking against Katz's fixed wing code
        ! DEBUG
        ! velTangentialChord(ic,is)=10._dp*cos(5._dp*pi/180._dp)
        ! velTangentialSpan(ic,is)=0._dp

        ! -1.0 multiplied to invert sign of gamma
        this%wiP(ic, is)%delPUnsteady = density * &
          & (this%wiP(ic, is)%gamTrapz - this%wiP(ic, is)%gamPrev)/dt

        this%wiP(ic, is)%delP = this%wiP(ic, is)%delPUnsteady + &
          & density*velTangentialChord(ic, is)* &
          & gamElementChord(ic, is) / this%wiP(ic, is)%meanChord
          
        if (this%spanwiseLiftSwitch .ne. 0) then
          this%wiP(ic, is)%delP = this%wiP(ic, is)%delP + density * &
            & velTangentialSpan(ic, is)* &
            & gamElementSpan(ic, is) / this%wiP(ic, is)%meanSpan
        endif

        ! -1.0 multiplied to invert sign of gamma
        this%wiP(ic, is)%gamPrev = this%wiP(ic, is)%gamTrapz

        ! Compute induced drag
        ! velInduced in either direction doesnt change drag direction
        this%wiP(ic, is)%delDiConstant = density*abs(velInduced(ic, is))* &
          & abs(gamElementChord(ic, is))*this%wiP(ic, is)%meanSpan
        this%wiP(ic, is)%delDiUnsteady = this%wiP(ic, is)%delPUnsteady* &
          & this%wiP(ic, is)%panelArea * &
          & dot_product(this%wiP(ic, is)%nCap, unitVec(this%wiP(ic, is)%velCPm))

        ! Invert direction of normalForce according to sign of omega and collective pitch
        this%wiP(ic, is)%normalForce = this%wiP(ic, is)%delP* &
          & this%wiP(ic, is)%panelArea*this%wiP(ic, is)%nCap*invertGammaSign

        this%wiP(ic, is)%normalForceUnsteady = this%wiP(ic, is)%delPUnsteady* &
          & this%wiP(ic, is)%panelArea*this%wiP(ic, is)%nCap*invertGammaSign

        this%secForceInertial(:, is) = this%secForceInertial(:, is) + this%wiP(ic, is)%normalForce

        this%secLift(:, is) = this%secLift(:, is) + projVec(this%wiP(ic, is)%normalForce, &
          this%secLiftDir(:, is))
        this%secLiftUnsteady(:, is) = this%secLiftUnsteady(:, is) + &
          & projVec(this%wiP(ic, is)%normalForceUnsteady, this%secLiftDir(:, is))
      enddo
      this%secDragInduced(:, is) = this%secDragDir(:, is)* &
        sum(this%wiP(:, is)%delDiConstant + this%wiP(:, is)%delDiUnsteady)
      ! Drag unsteady is purely for monitoring purposes if required
      ! and is not used for computations anywhere
      this%secDragUnsteady(:, is) = this%secDragDir(:, is)* &
        sum(this%wiP(:, is)%delDiUnsteady)
    enddo

    ! To overwrite unit vectors previously assigned in main.f90
    this%secDragProfile = 0._dp  

    this%secDrag = this%secDragInduced + this%secDragProfile

    ! Compute sectional coefficients
    ! Compute secChordwiseResVel for calculating secDynamicPressure
    ! including induced velocities
    call this%calc_secChordwiseResVel()

    ! Alternately assuming sec resultant velocity is same as sec freestream vel
    ! do is = 1, this%ns
    !   this%secChordwiseResVel(:, is) = this%secVelFreestream(:, is) - &
    !     & dot_product(this%secVelFreestream(:, is), this%yAxis)*this%yAxis
    ! enddo
    secDynamicPressure = this%getSecDynamicPressure(density)

    do is = 1, this%ns
      if (abs(secDynamicPressure(is)) > eps) then
        ! Use sign of delP to obtain sign of CL
        signSecCL = sign(1._dp, sum(this%wiP(:, is)%delP))
        this%secCL(is) = norm2(this%secLift(:, is))*signSecCL/ &
          & (secDynamicPressure(is)*this%secArea(is))
        this%secCD(is) = norm2(this%secDrag(:, is))/ &
          & (secDynamicPressure(is)*this%secArea(is))
        this%secCLu(is) = norm2(this%secLiftUnsteady(:, is))*signSecCL/ &
          & (secDynamicPressure(is)*this%secArea(is))
      else
        this%secCL(is) = 0._dp
        this%secCD(is) = 0._dp
        this%secCLu(is) = 0._dp
      endif
    enddo

    call this%sumSecToNetForces()

  end subroutine blade_calc_force_gamma

  function blade_getSecDynamicPressure(this, density)
  class(blade_class), intent(in) :: this
    real(dp), intent(in) :: density
    real(dp), dimension(this%ns) :: magsecVelCPTotal
    real(dp), dimension(this%ns) :: blade_getSecDynamicPressure
    integer :: is

    do is = 1, this%ns
      magsecVelCPTotal(is) = norm2(this%secChordwiseResVel(:, is))
    enddo
    blade_getSecDynamicPressure = 0.5_dp*density*magsecVelCPTotal**2._dp
  end function blade_getSecDynamicPressure

  subroutine calc_secArea(this)
  class(blade_class), intent(inout) :: this
    integer :: is

    do is = 1, this%ns
      this%secArea(is) = sum(this%wiP(:, is)%panelArea)
    enddo
  end subroutine calc_secArea

  subroutine calc_secChord(this)
  class(blade_class), intent(inout) :: this
    integer :: is

    do is = 1, this%ns
      this%secChord(is) = norm2(0.5_dp*( &
        & (this%wiP(1, is)%PC(:,1) + this%wiP(1, is)%PC(:,4)) - &
        & (this%wiP(this%nc, is)%PC(:,2) + this%wiP(this%nc, is)%PC(:,3))))
    enddo
  end subroutine calc_secChord

  subroutine blade_calc_force_alpha(this, density, velSound)
    ! Compute force using sec alpha
  class(blade_class), intent(inout) :: this
    real(dp), intent(in) :: density, velSound

    call this%lookup_secCoeffs(velSound)

    call this%secCoeffsTosecForces(density)
    ! To overwrite unit vectors previously assigned in main.f90
    this%secDragInduced = 0._dp
    this%secDrag = this%secDragProfile + this%secDragInduced

    call this%sumSecToNetForces()

  end subroutine blade_calc_force_alpha

  subroutine blade_calc_force_alphaGamma(this, density, &
      & invertGammaSign, velSound, dt)
    ! Compute force using alpha approximated from sec circulation
    use libMath, only: unitVec, cross_product
  class(blade_class), intent(inout) :: this
    real(dp), intent(in) :: density, invertGammaSign, velSound, dt
    real(dp), dimension(3) :: secChordwiseVelFreestream, liftDir
    real(dp), dimension(this%ns) :: secDynamicPressure
    integer :: is

    ! Compute unsteady sec lift from gamma distribution
    call this%calc_force_gamma(density, invertGammaSign, dt)

    do is = 1, this%ns
      ! Compute sec freestream velocity
      secChordwiseVelFreestream = this%secVelFreestream(:, is) - &
        dot_product(this%secVelFreestream(:, is), this%yAxis)*this%yAxis

      ! Assuming sec resultant velocity is same as sec freestream vel
      ! for computing corrected alpha later
      this%secChordwiseResVel(:, is) = secChordwiseVelFreestream
    enddo

    secDynamicPressure = this%getSecDynamicPressure(density)

    do is = 1, this%ns
      ! Extract sectional lift and CL
      liftDir = cross_product(this%secChordwiseResVel(:, is), this%yAxis)
      ! Assuming gamma method only gives lift
      this%secCL(is) = dot_product(this%secForceInertial(:, is), & 
        unitVec(liftDir)) / (secDynamicPressure(is)*this%secArea(is))

      ! Compute angle of attack from linear CL
      this%secArea(is) = this%secCL(is)/twoPi + &
        & this%alpha0(this%airfoilNo(is))
    enddo
    this%secCLu = 0._dp

    ! Compute non-linear CL
    call this%lookup_secCoeffs(velSound)

    call this%secCoeffsToSecForces(density)

    this%secDrag = this%secDragInduced + this%secDragProfile

    call this%sumSecToNetForces()

  end subroutine blade_calc_force_alphaGamma

  subroutine blade_secCoeffsToSecForces(this, density)
    use libMath, only: unitVec, cross_product, lsq2
  class(blade_class), intent(inout) :: this
    real(dp), intent(in) :: density
    integer :: is
    real(dp), dimension(this%ns) :: leadingTerm

    leadingTerm = this%getSecDynamicPressure(density)*this%secArea

    do is = 1, this%ns
      ! Lift and Drag vectors
      this%secLift(:, is) = &
        & this%secLiftDir(:, is)*leadingTerm(is)*this%secCL(is) 
      this%secDragProfile(:, is) = &
        & this%secDragDir(:, is)*leadingTerm(is)*this%secCD(is)
      this%secLiftUnsteady(:, is) = &
        & this%secLiftDir(:, is)*leadingTerm(is)*this%secCLu(is) 
      ! Lift in inertial frame
      ! Warning: This would give a wrong answer if a considerable dihedral
      ! is present for the wing since the blade Y-axis is not flapped
      this%secForceInertial(:, is) = cross_product(this%yAxis, &
        this%secChordwiseResVel(:, is))
      this%secForceInertial(:, is) = sign(1._dp, sum(this%wiP(:, is)%vr%gam)) &
        * unitVec(this%secForceInertial(:, is))
      ! abs() used since direction is already captured in vector
      this%secForceInertial(:, is) = norm2(this%secLift(:,is)) * &
        this%secForceInertial(:, is)
      ! Drag in inertial frame
      this%secForceInertial(:,is) = this%secForceInertial(:,is) + &
        norm2(this%secDrag(:,is)) * unitVec(this%secChordwiseResVel(:,is))
    enddo
  end subroutine blade_secCoeffsToSecForces

  subroutine blade_lookup_secCoeffs(this, velSound)
    use libC81, only: getCL, getCD
    ! Compute sec CL, CD, CM from C81 tables and sec resultant velocity
    ! Assumes only one airfoil section present
  class(blade_class), intent(inout) :: this
    real(dp), intent(in) :: velSound
    real(dp) :: secMach, alphaDeg
    integer :: is

    do is = 1, size(this%secAlpha, 1)
      secMach = norm2(this%secChordwiseResVel(:, is))/velSound
      alphaDeg = this%secAlpha(is) * 180._dp / pi
      this%secCL(is) = this%C81(this%airfoilNo(is))%getCL(alphaDeg, secMach)
      this%secCD(is) = this%C81(this%airfoilNo(is))%getCD(alphaDeg, secMach)
      this%secCM(is) = this%C81(this%airfoilNo(is))%getCM(alphaDeg, secMach)
    enddo
    this%secCLu = 0._dp
  end subroutine blade_lookup_secCoeffs

  subroutine blade_calc_secChordwiseResVel(this)
    ! Compute sec resultant velocity by interpolating local panel velocity
    use libMath, only: lsq2
  class(blade_class), intent(inout) :: this
    integer :: i, is, ic
    real(dp), dimension(this%nc) :: xDist

    if (this%nc .ge. 3) then  
      ! Use least squares fit to get sec resultant velocity
      do is = 1, size(this%secChordwiseResVel, 2)
        do ic = 1, this%nc
          call this%wiP(ic, is)%calc_chordwiseResVel()
          xDist(ic) = dot_product(this%wiP(ic, is)%CP &
            & -this%wiP(1, is)%PC(:, 1), &
            & this%secTauCapChord(:, is))
          ! DEBUG
          if (is ==1 ) print*, xDist(ic)
        enddo
        do i = 1, 3
          this%secChordwiseResVel(i, is) = lsq2(dot_product(this%secCP(:, is) &
            & -this%wiP(1, is)%PC(:, 1), this%secTauCapChord(:, is)), xDist, &
            this%wiP(:, is)%chordwiseResVel(i))
        enddo
      enddo
    else  ! Use average of resultant velocities
      do is = 1, size(this%secChordwiseResVel, 2)
        do ic = 1, this%nc
          call this%wiP(ic, is)%calc_chordwiseResVel()
        enddo
        do i = 1, 3
          this%secChordwiseResVel(i, is) = &
            & sum(this%wiP(:, is)%chordwiseResVel(i))/this%nc
        enddo
      enddo
    endif

    ! DEBUG
    print*, 'v', this%secChordwiseResVel(:, 1)
    do ic = 1, this%nc
      print*, this%wiP(ic, 1)%chordwiseResVel
    enddo
  end subroutine blade_calc_secChordwiseResVel

  subroutine blade_calc_secAlpha(this, updateSecVel)
    ! Compute sec alpha using sec resultant velocity
    use libMath, only: getAngleTan
  class(blade_class), intent(inout) :: this
    logical, intent(in) :: updateSecVel
    integer :: is

    if (updateSecVel) call this%calc_secChordwiseResVel()

    !! Use acos() to find angle
    !do is=1,size(this%secAlpha)
    !  this%secAlpha(is)=acos(dot_product(this%secChordwiseResVel(:,is), &
    !    & this%secTauCapChord(:,is)) &
    !    & /norm2(this%secChordwiseResVel(:,is)))
    !enddo

    ! Use atan2() to find angle
    do is = 1, size(this%secAlpha)
      this%secAlpha(is) = &
        & atan2(dot_product(this%secChordwiseResVel(:, is), &
        & this%secNormalVec(:, is)), &
        & dot_product(this%secChordwiseResVel(:, is), &
        & this%secTauCapChord(:, is)))

      this%secPhi(is) = getAngleTan(this%secChordwiseResVel(:, is), &
        & this%wiP(1, is)%velCPm)
    enddo
  end subroutine blade_calc_secAlpha

  !subroutine blade_calc_secAlpha(this)
  !  ! Compute sec alpha by interpolating local panel alpha
  !class(blade_class), intent(inout) :: this
  !  integer :: is, ic, this%nc
  !  real(dp), dimension(size(this%wiP,1)) :: xDist

  !  if (this%nc .ge. 3) then  ! Use least squares fit to get alpha
  !    do is=1,size(this%secAlpha)
  !      do ic=1,this%nc
  !        xDist(ic)=dot_product(this%wiP(ic,is)%CP-this%wiP(1,is)%PC(:,1),  &
  !          this%secTauCapChord(:,is))
  !      enddo
  !      this%secAlpha(is)=lsq2(dot_product(this%secCP(:,is)-  &
  !        this%wiP(1,is)%PC(:,1),this%secTauCapChord(:,is)),xDist,this%wiP(:,is)%alpha)
  !    enddo
  !  else  ! Use average of alpha values
  !    do is=1,size(this%secAlpha)
  !      this%secAlpha(is)=sum(this%wiP(:,is)%alpha)/this%nc
  !    enddo
  !  endif
  !end subroutine blade_calc_secAlpha

  function blade_getSecChordwiseLocations(this, chordwiseFraction)
    ! Get coordinates of a point located at a fraction of chord on each section
  class(blade_class), intent(inout) :: this
    real(dp), intent(in) :: chordwiseFraction
    real(dp), dimension(3, this%ns) :: blade_getSecChordwiseLocations
    integer :: is

    do is = 1, this%ns
      blade_getSecChordwiseLocations(:, is) = &
        & (1._dp-chordwiseFraction)*(this%wiP(1, is)%PC(:, 4) &
        & + this%wiP(1, is)%PC(:, 1))*0.5_dp &
        & + chordwiseFraction*(this%wiP(this%nc, is)%PC(:, 3) &
        & + this%wiP(this%nc, is)%PC(:, 2))*0.5_dp
    enddo
  end function blade_getSecChordwiseLocations

  subroutine blade_burst_wake(this, rowFar, skewLimit, largeCoreRadius)
    use libMath, only: getAngleCos
  class(blade_class), intent(inout) :: this
    real(dp), intent(in) :: skewLimit
    integer, intent(in) :: rowFar !, rowNear
    real(dp), intent(in) :: largeCoreRadius
    integer :: irow !, icol
    real(dp) :: skewVal

    !! Burst near wake
    !do icol=1,size(this%waN,2)
    !  do irow=rowNear,size(this%waN,1)
    !    call this%waN(irow,icol)%vr%burst(skewLimit)
    !  enddo
    !enddo

    ! Burst far wake
    if (rowFar .le. size(this%waF, 1)) then
      do irow = rowFar, size(this%waF, 1) - 1
        ! NEGLECT ALREADY BURST FILAMENTS IF NECCESSARY
        !if (abs(this%waF(irow+1)%gam) > eps .and. abs(this%waF(irow)%gam) > eps) then
        skewVal = abs(getAngleCos(this%waF(irow)%vf%fc(:, 2) - this%waF(irow)%vf%fc(:, 1) &
          , this%waF(irow + 1)%vf%fc(:, 1) - this%waF(irow + 1)%vf%fc(:, 2) &
          ) - pi)/pi
        if (skewVal .ge. skewLimit) then
          !this%waF(irow+1)%gam = 0._dp
          !this%waF(irow)%gam = 0._dp
          this%waF(irow + 1)%vf%rVc = largeCoreRadius
          this%waF(irow)%vf%rVc = largeCoreRadius
        endif
        !endif
      enddo
    endif
  end subroutine blade_burst_wake

  subroutine blade_calc_skew(this, rowNear)
  class(blade_class), intent(inout) :: this
    integer, intent(in) :: rowNear
    integer :: irow, icol

    !$omp parallel do collapse(2)
    do icol = 1, size(this%waN, 2)
      do irow = rowNear, size(this%waN, 1)
        call this%waN(irow, icol)%vr%calc_skew()
      enddo
    enddo
    !$omp end parallel do
  end subroutine blade_calc_skew

  subroutine blade_dirLiftDrag(this)
    use libMath, only: unitVec, cross_product
  class(blade_class), intent(inout) :: this
    integer :: is
    do is = 1, this%ns
      this%secDragDir(:, is) = unitVec(this%wiP(1, is)%velCPm)
      this%secLiftDir(:, is) = cross_product(this%secDragDir(:, is), &
        & this%yAxis)
    enddo
  end subroutine blade_dirLiftDrag

  subroutine blade_sumSecToNetForces(this)
  class(blade_class), intent(inout) :: this

    this%forceInertial = sum(this%secForceInertial, 2)
    this%lift = sum(this%secLift, 2)
    this%drag = sum(this%secDrag, 2)
    this%liftUnsteady = sum(this%secLiftUnsteady, 2)
    this%dragProfile = sum(this%secDragProfile, 2)
    this%dragInduced = sum(this%secDragInduced, 2)
    this%dragUnsteady = sum(this%secDragUnsteady, 2)
  end subroutine blade_sumSecToNetForces

  subroutine blade_calc_stlStats(this)
  class(blade_class), intent(inout) :: this
    real(dp), dimension(3) :: vertex
    integer :: icell, ivertex, irow
    logical :: isNewVertex
    ! Dummy arrays
    real(dp), allocatable, dimension(:, :) :: stlNodes

    ! Naive allocations
    allocate(stlNodes(3, this%nc*3))
    allocate(this%stlElementNodes(3,this%nc))

    this%stlNodesCols = 0
    do icell = 1, this%nc
      do ivertex = 1, 3
        vertex = this%wiP(icell, 1)%PC(:, ivertex)

        ! Check if vertex is present in vertices array
        isNewVertex = .true.
        do irow = 1, this%stlNodesCols
          if (all(abs(stlNodes(:, irow) - vertex) < eps)) then
            isNewVertex = .false.
            exit
          endif
        enddo

        if (isNewVertex) then
          ! Record to vertices array
          this%stlNodesCols = this%stlNodesCols + 1
          stlNodes(:, this%stlNodesCols) = vertex
        endif
      enddo
    enddo

    allocate(this%stlNodes(3, this%stlNodesCols))
    this%stlNodes = stlNodes(:, 1:this%stlNodesCols)

    ! Deallocate dummy arrays
    deallocate(stlNodes)

    ! For each element find the 3 node indices in stlNodes array 
    do icell = 1, this%nc
      do ivertex = 1,3
        vertex = this%wiP(icell, 1)%PC(:, ivertex)

        do irow = 1, this%stlNodesCols
          if (all(abs(this%stlNodes(:, irow) - vertex) < eps)) then
            ! Record row number
            this%stlElementNodes(ivertex, icell) = irow
            exit
          endif
        enddo

      enddo
    enddo
  end subroutine blade_calc_stlStats

  subroutine blade_write(this, unit, iostat, iomsg)
  class(blade_class), intent(in) :: this
    integer, intent(in) :: unit
    integer, intent(out) :: iostat
    character(len=*), intent(inout) :: iomsg

    write(unit, iostat=iostat, iomsg=iomsg) this%wiP, this%waN, this%waF, &
      & this%waNPredicted, this%waFPredicted, &
      & this%theta, this%psi, &
      & this%forceInertial, this%lift, this%liftUnsteady, & 
      & this%drag, this%dragInduced, this%dragProfile, this%dragUnsteady, & 
      & this%velNwake, this%velNwake1, this%velNwake2, this%velNwake3, &
      & this%velNwakePredicted, this%velNwakeStep, &
      & this%velFwake, this%velFwake1, this%velFwake2, this%velFwake3, &
      & this%velFwakePredicted, this%velFwakeStep, &
      & this%xAxis, this%yAxis, this%zAxis, &
      & this%secChord, this%secArea, &
      & this%secForceInertial, this%secLift, this%secDrag, &
      & this%secLiftDir, this%secDragDir, &
      & this%secDragInduced, this%secDragProfile, &
      & this%secLiftUnsteady, this%secDragUnsteady, &
      & this%secTauCapChord, this%secTauCapSpan, &
      & this%secNormalVec, this%secVelFreestream, &
      & this%secChordwiseResVel, this%secCP, &
      & this%secAlpha, this%secPhi, this%secCL, this%secCD, this%secCM, &
      & this%secCLu, this%alpha0
  end subroutine blade_write

  subroutine blade_read(this, unit, iostat, iomsg)
  class(blade_class), intent(inout) :: this
    integer, intent(in) :: unit
    integer, intent(out) :: iostat
    character(len=*), intent(inout) :: iomsg

    read(unit, iostat=iostat, iomsg=iomsg) this%wiP, this%waN, this%waF, &
      & this%waNPredicted, this%waFPredicted, &
      & this%theta, this%psi, &
      & this%forceInertial, this%lift, this%liftUnsteady, & 
      & this%drag, this%dragInduced, this%dragProfile, this%dragUnsteady, & 
      & this%velNwake, this%velNwake1, this%velNwake2, this%velNwake3, &
      & this%velNwakePredicted, this%velNwakeStep, &
      & this%velFwake, this%velFwake1, this%velFwake2, this%velFwake3, &
      & this%velFwakePredicted, this%velFwakeStep, &
      & this%xAxis, this%yAxis, this%zAxis, &
      & this%secChord, this%secArea, &
      & this%secForceInertial, this%secLift, this%secDrag, &
      & this%secLiftDir, this%secDragDir, &
      & this%secDragInduced, this%secDragProfile, &
      & this%secLiftUnsteady, this%secDragUnsteady, &
      & this%secTauCapChord, this%secTauCapSpan, &
      & this%secNormalVec, this%secVelFreestream, &
      & this%secChordwiseResVel, this%secCP, &
      & this%secAlpha, this%secPhi, this%secCL, this%secCD, this%secCM, &
      & this%secCLu, this%alpha0
  end subroutine blade_read

  !------+--
  ! ++++ | rotor_class Methods
  !------+--

  subroutine rotor_read_geom(this, filename)
    use libMath, only: skip_comments
  class(rotor_class) :: this
    character(len=*), intent(in) :: filename
    integer :: i
    character(len=10) :: fileFormatVersion, currentTemplateVersion

    currentTemplateVersion = '0.10'

    open (unit=12, file=filename, status='old', action='read')
    call skip_comments(12)
    read(12, *) fileFormatVersion
    if (adjustl(fileFormatVersion) /= currentTemplateVersion) then
      error stop "ERROR: geomXX.in template version does not match"
    endif

    call skip_comments(12)
    read (12, *) this%nb, this%propConvention, &
      & this%nCamberFiles, this%geometryFile
    call skip_comments(12)

    if (this%nCamberFiles > 0) then
      allocate(this%camberSectionLimit(this%nCamberFiles))
      allocate(this%camberFile(this%nCamberFiles))
      do i = 1, this%nCamberFiles
        read(12, *) this%camberSectionLimit(i), this%camberFile(i)
        call skip_comments(12)
      enddo
    else
      allocate(this%camberSectionLimit(1))
      allocate(this%camberFile(1))
      ! Default uncambered section
      this%camberSectionLimit = 1.0
      this%camberFile = '0'
      read(12, *)
      call skip_comments(12)
    endif

    ! [0/1]Lifting [2]Non-lifting [-1]Lifting Image [-2]Non-lifting Image
    read (12, *) this%surfaceType, this%imagePlane, this%imageRotorNum
    call skip_comments(12)
    read (12, *) this%nc, this%ns, this%nNwake
    call skip_comments(12)
    read (12, *) this%hubCoords(1), this%hubCoords(2), this%hubCoords(3)
    call skip_comments(12)
    read (12, *) this%cgCoords(1), this%cgCoords(2), this%cgCoords(3)
    call skip_comments(12)
    read (12, *) this%fromCoords(1), this%fromCoords(2), this%fromCoords(3)
    call skip_comments(12)
    read (12, *) this%pts(1), this%pts(2), this%pts(3)
    call skip_comments(12)
    read (12, *) this%radius, this%root_cut, this%chord, this%coningAngle
    call skip_comments(12)
    read (12, *) this%Omega, this%shaftAxis(1), this%shaftAxis(2), this%shaftAxis(3)
    call skip_comments(12)
    read (12, *) this%controlPitch(1), this%controlPitch(2), this%controlPitch(3), this%thetaTwist
    call skip_comments(12)
    read(12, *) this%customTrajectorySwitch, this%imposeAxisymmetry
    call skip_comments(12)
    read (12, *) this%velBody(1), this%velBody(2), this%velBody(3) &
      , this%omegaBody(1), this%omegaBody(2), this%omegaBody(3)
    call skip_comments(12)
    read (12, *) this%pivotLE, this%flapHinge, this%spanwiseLiftSwitch, &
      & this%symmetricTau
    call skip_comments(12)
    read (12, *) this%apparentViscCoeff, this%decayCoeff
    call skip_comments(12)
    read(12, *) this%wakeTruncateNt, &
      & this%prescWakeAfterTruncNt, this%prescWakeGenNt
    call skip_comments(12)
    read (12, *) this%spanwiseCore, this%streamwiseCoreSwitch
    call skip_comments(12)
    allocate (this%streamwiseCoreVec(this%ns + 1))
    if (this%streamwiseCoreSwitch .eq. 'i') then  ! [i]dentical
      read (12, *) this%streamwiseCoreVec(1)
      do i = 2, this%ns + 1
        this%streamwiseCoreVec(i) = this%streamwiseCoreVec(1)
      enddo
    elseif (this%streamwiseCoreSwitch .eq. 's') then  ! [s]ectional
      read (12, *) (this%streamwiseCoreVec(i), i=1, this%ns + 1)
    else
      error stop 'ERROR: Wrong input for streamwiseCoreSwitch in geomXX.in'
    endif
    call skip_comments(12)
    ! Dimensional quantities
    read (12, *) this%rollupStartRadius, this%rollupEndRadius
    call skip_comments(12)
    read (12, *) this%initWakeVel, &
      & this%psiStart, this%skewLimit
    call skip_comments(12)
    read (12, *) this%dragUnitVec(1), this%dragUnitVec(2), this%dragUnitVec(3)
    call skip_comments(12)
    read (12, *) this%sideUnitVec(1), this%sideUnitVec(2), this%sideUnitVec(3)
    call skip_comments(12)
    read (12, *) this%liftUnitVec(1), this%liftUnitVec(2), this%liftUnitVec(3)
    call skip_comments(12)
    read (12, *) this%inflowPlotSwitch, this%gammaPlotSwitch
    call skip_comments(12)
    read (12, *) this%skewPlotSwitch
    call skip_comments(12)
    read (12, *) this%forceCalcSwitch, this%nAirfoils
    ! Ensure airfoil tables are provided when force calculation requires them
    if (this%forceCalcSwitch .gt. 0 .and. this%nAirfoils .eq. 0) then 
      error stop 'ERROR: No. of airfoil tables set to 0 in geomXX.in'
    endif
    call skip_comments(12)
    if (this%nAirfoils .gt. 0) then
      allocate (this%airfoilSectionLimit(this%nAirfoils))
      allocate (this%airfoilFile(this%nAirfoils))
      allocate (this%alpha0(this%nAirfoils))
      do i = 1, this%nAirfoils
        read (12, *) this%airfoilSectionLimit(i), this%alpha0(i), this%airfoilFile(i)
      enddo
    endif
    close (12)
  end subroutine rotor_read_geom

  subroutine rotor_init(this, rotorNumber, density, dt, nt, switches)
    ! Initialize variables of rotor geometry and wake
    use libMath
  class(rotor_class) :: this
    integer, intent(in) :: rotorNumber
    real(dp), intent(in) :: density
    real(dp) , intent(inout) :: dt
    integer, intent(inout) :: nt
    type(switches_class), intent(inout) :: switches

    real(dp), dimension(this%nc+1) :: xVec
    real(dp), dimension(this%ns+1) :: yVec
    real(dp), dimension(this%nc+1, this%ns+1) :: zVec
    real(dp), dimension(this%nc, this%ns) :: dx, dy
    real(dp), dimension(3) :: leftTipCP
    real(dp) :: dxdymin, secCPLoc, rbyR, dxMAC
    integer :: i, j, ib, is, ic
    real(dp) :: bladeOffset
    real(dp) :: velShed
    real(dp), dimension(4) :: xshift
    character(len=2) :: rotorChar
    logical :: warnUser

    ! Set dt automatically if not prescribed
    ! if dt is negative, assume no. of chords or revs
    if (dsign(1._dp, dt) < 0._dp) then
      if (abs(this%Omega) < eps) then  ! Fixed wing
        dt = abs(dt)*this%chord/norm2(this%velBody)
      else  ! Rotor
        dt = twoPi*abs(dt)/abs(this%Omega)
      endif
      print*, 'dt set to ', dt
    endif

    ! if dt is 0, assume default timesteps
    if (abs(dt) <= eps) then
      if (abs(this%Omega) < eps) then  ! Fixed wing
        dxMAC = this%chord/this%nc
        dt = dxMAC/norm2(this%velBody)
      else  ! Rotor
        ! Time for 5 deg
        dt = 5._dp*pi/(180._dp*abs(this%Omega))
      endif
      print*, 'dt set to ', dt
    endif

    ! Check if nt requires preprocessing
    if (nt .le. 0) then
      if (nt .eq. 0) nt = -10
      call this%toChordsRevs(nt, dt)
      print*, 'nt set to ', nt
    endif

    ! Preprocessing for negative value inputs
    if (switches%slowStart /= 0) then
      call this%toChordsRevs(switches%slowStartNt, dt)
    endif

    call this%toChordsRevs(switches%wakeTipPlot, dt)
    call this%toChordsRevs(switches%wakePlot, dt)
    call this%toChordsRevs(switches%gridPlot, dt)
    call this%toChordsRevs(this%wakeTruncateNt, dt)
    call this%toChordsRevs(this%prescWakeAfterTruncNt, dt)
    call this%toChordsRevs(this%prescWakeGenNt, dt)

    call this%toChordsRevs(this%nNwake, dt)
    call this%toChordsRevs(this%inflowPlotSwitch, dt)
    call this%toChordsRevs(this%gammaPlotSwitch, dt)
    call this%toChordsRevs(this%skewPlotSwitch, dt)

    if (this%wakeTruncateNt > 0 .and. this%prescWakeAfterTruncNt > 0) then
      this%prescWakeNt = this%wakeTruncateNt+this%prescWakeAfterTruncNt
    else
      this%prescWakeNt = 0
    endif

    if (abs(this%Omega) < eps .and. this%wakeTruncateNt > 0) then
      this%wakeTruncateNt = 0
    endif

    if (this%surfaceType == 0) this%surfaceType = 1
    if (abs(this%surfaceType) == 2) this%suppressFwakeSwitch = 0

    if (this%nNwake < 2) error stop 'ERROR: Atleast 2 near wake rows mandatory'

    ! Override ns to 1 if non-lifting surface
    if (abs(this%surfaceType) .eq. 2) this%ns = 1

    ! Initialize variables for use in allocating
    if (abs(this%surfaceType) == 1) then
      this%nNwake = min(this%nNwake, nt)
      this%nFwake = nt - this%nNwake
    elseif (abs(this%surfaceType) .eq. 2) then
      this%nNwake = 0
      this%nFwake = 0
    endif

    ! Define limits for use in wake convection
    this%nbConvect = this%nb
    if (this%imposeAxisymmetry == 1) this%nbConvect = 1
    this%nNwakeEnd = this%nNwake
    this%nFwakeEnd = this%nFwake

    ! Allocate rotor object variables
    allocate (this%blade(this%nb))
    allocate (this%AIC(this%nc*this%ns*this%nb, this%nc*this%ns*this%nb))
    allocate (this%AIC_inv(this%nc*this%ns*this%nb, this%nc*this%ns*this%nb))
    allocate (this%gamVec(this%nc*this%ns*this%nb))
    allocate (this%gamVecPrev(this%nc*this%ns*this%nb))
    allocate (this%RHS(this%nc*this%ns*this%nb))

    ! Read custom trajectory file if specified
    if (this%customTrajectorySwitch .eq. 1) then
      write(rotorChar, '(I0.2)') rotorNumber
      allocate(this%velBodyHistory(3, nt))
      allocate(this%omegaBodyHistory(3, nt))

      open (unit=13, file='trajectory'//rotorChar//'.in', &
        & status='old', action='read')
      call skip_comments(13)
      do i = 1, nt
        read(13, *) this%velBodyHistory(1, i), &
          & this%velBodyHistory(2, i), &
          & this%velBodyHistory(3, i), &
          & this%omegaBodyHistory(1, i), &
          & this%omegaBodyHistory(2, i), &
          & this%omegaBodyHistory(3, i)
      enddo
      close(13)
    endif

    ! Allocate blade object variables
    do ib = 1, this%nb
      allocate (this%blade(ib)%wiP(this%nc, this%ns))
      allocate (this%blade(ib)%waN(this%nNwake, this%ns))
      allocate (this%blade(ib)%waF(this%nFwake))
      allocate (this%blade(ib)%secTauCapChord(3, this%ns))
      allocate (this%blade(ib)%secTauCapSpan(3, this%ns))
      allocate (this%blade(ib)%secNormalVec(3, this%ns))
      allocate (this%blade(ib)%secVelFreestream(3, this%ns))
      allocate (this%blade(ib)%secForceInertial(3, this%ns))
      allocate (this%blade(ib)%secChord(this%ns))
      allocate (this%blade(ib)%secArea(this%ns))
      allocate (this%blade(ib)%secLift(3, this%ns))
      allocate (this%blade(ib)%secLiftUnsteady(3, this%ns))
      allocate (this%blade(ib)%secDrag(3, this%ns))
      allocate (this%blade(ib)%secLiftDir(3, this%ns))
      allocate (this%blade(ib)%secDragDir(3, this%ns))
      allocate (this%blade(ib)%secDragInduced(3, this%ns))
      allocate (this%blade(ib)%secDragProfile(3, this%ns))
      allocate (this%blade(ib)%secDragUnsteady(3, this%ns))
      allocate (this%blade(ib)%secAlpha(this%ns))
      allocate (this%blade(ib)%secPhi(this%ns))
      allocate (this%blade(ib)%airfoilNo(this%ns))
      allocate (this%blade(ib)%secCL(this%ns))
      allocate (this%blade(ib)%secCD(this%ns))
      allocate (this%blade(ib)%secCM(this%ns))
      allocate (this%blade(ib)%secCLu(this%ns))
      allocate (this%blade(ib)%secChordwiseResVel(3, this%ns))
      allocate (this%blade(ib)%secCP(3, this%ns))
    enddo

    ! Conversions
    do i = 1, 3
      call degtorad(this%controlPitch(i))
      call degtorad(this%pts(i))
    enddo
    call degtorad(this%thetaTwist)
    call degtorad(this%coningAngle)
    call degtorad(this%psiStart)

    this%spanwiseCore = this%spanwiseCore*this%chord
    this%streamwiseCoreVec = this%streamwiseCoreVec*this%chord
    this%rollupStart = ceiling(this%rollupStartRadius*this%ns)
    this%rollupEnd = floor(this%rollupEndRadius*this%ns)

    ! Rotor initialization
    this%gamVec = 0._dp
    this%gamVecPrev = 0._dp

    ! Blade initialization
    if (this%geometryFile(1:1) .eq. '0') then
      select case (switches%chordSpacing)
      case (1)
        if (this%Omega .ge. 0) then
          xVec = linspace(-this%chord, 0._dp, this%nc + 1)
        else
          xVec = linspace(this%chord, 0._dp, this%nc + 1)
        endif
      case (2)
        if (this%Omega .ge. 0) then
          xVec = cosspace(-this%chord, 0._dp, this%nc + 1)
        else
          xVec = cosspace(this%chord, 0._dp, this%nc + 1)
        endif
      case (3)
        if (this%Omega .ge. 0) then
          xVec = halfsinspace(-this%chord, 0._dp, this%nc + 1)
        else
          xVec = halfsinspace(this%chord, 0._dp, this%nc + 1)
        endif
      end select

      select case (switches%spanSpacing)
      case (1)
        yVec = linspace(this%root_cut*this%radius, this%radius, this%ns + 1)
      case (2)
        yVec = cosspace(this%root_cut*this%radius, this%radius, this%ns + 1)
      case (3)
        yVec = halfsinspace(this%root_cut*this%radius, this%radius, this%ns + 1)
      end select

      ! Initialize panel coordinates
      do ib = 1, this%nb
        ! Compute camber
        if (this%nCamberFiles > 0) then
          zVec = this%getCamber(xVec, yVec)
        else
          zVec = 0._dp
        endif
        if (this%surfaceType < 0) then
          zVec = -1._dp*zVec
        endif

        ! Assign coordinates to panels
        do j = 1, this%ns
          do i = 1, this%nc
            call this%blade(ib)%wiP(i, j)%assignP(1, &
              & (/xVec(i), yVec(j), zVec(i, j)/))
            call this%blade(ib)%wiP(i, j)%assignP(2, &
              & (/xVec(i+1), yVec(j), zVec(i+1, j)/))
            call this%blade(ib)%wiP(i, j)%assignP(3, &
              & (/xVec(i+1), yVec(j+1), zVec(i+1, j+1)/))
            call this%blade(ib)%wiP(i, j)%assignP(4, &
              & (/xVec(i), yVec(j+1), zVec(i, j+1)/))
          enddo
        enddo
      enddo

    else
      if (abs(this%surfaceType) == 1) then
        call this%plot3dtoblade(trim(this%geometryFile))
      elseif (abs(this%surfaceType) == 2) then
        call this%stltoblade(trim(this%geometryFile))
        do ib = 1, this%nb
          this%blade(ib)%nc = this%nc
          this%blade(ib)%ns = this%ns
          call this%blade(ib)%calc_stlStats()
        enddo
      endif
    endif

    ! Copy ns and nc to blade variables to avoid recomputing
    do ib = 1, this%nb
      this%blade(ib)%nc = this%nc
      this%blade(ib)%ns = this%ns
    enddo

    do ib = 1, this%nb
      ! Initialize blade axes
      ! Can go wrong if importing geometry from OpenVSP
      ! and then using these for orientation of lift vectors
      this%blade(ib)%xAxis = xAxis
      this%blade(ib)%yAxis = yAxis
      this%blade(ib)%zAxis = zAxis

      ! Initialize sec vectors
      if (abs(this%surfaceType) == 1) then
        do j = 1, this%ns
          this%blade(ib)%secTauCapChord(:, j) = &
            (this%blade(ib)%wiP(this%nc, j)%PC(:, 3) + this%blade(ib)%wiP(this%nc, j)%PC(:, 2) - &
            (this%blade(ib)%wiP(1, j)%PC(:, 4) + this%blade(ib)%wiP(1, j)%PC(:, 1)))*0.5_dp

          this%blade(ib)%secTauCapSpan(:, j) = yAxis

          this%blade(ib)%secNormalVec(:, j) = &
            cross_product(this%blade(ib)%wiP(this%nc, j)%PC(:, 2) - this%blade(ib)%wiP(1, j)%PC(:, 4), &
            this%blade(ib)%wiP(this%nc, j)%PC(:, 3) - this%blade(ib)%wiP(1, j)%PC(:, 1))

          ! Normalize
          this%blade(ib)%secTauCapChord(:, j) = unitVec(this%blade(ib)%secTauCapChord(:, j))

          this%blade(ib)%secNormalVec(:, j) = unitVec(this%blade(ib)%secNormalVec(:, j))
        enddo

        ! Initialize vr coords of all panels except last row (to accomodate mismatch of vr coords when using unequal spacing)
        do j = 1, this%ns
          do i = 1, this%nc - 1
            xshift(1) = (this%blade(ib)%wiP(i, j)%PC(1, 2) &
              - this%blade(ib)%wiP(i, j)%PC(1, 1))*0.25_dp
            xshift(2) = (this%blade(ib)%wiP(i + 1, j)%PC(1, 2) &
              - this%blade(ib)%wiP(i, j)%PC(1, 2))*0.25_dp
            xshift(3) = (this%blade(ib)%wiP(i + 1, j)%PC(1, 3) &
              - this%blade(ib)%wiP(i, j)%PC(1, 3))*0.25_dp
            xshift(4) = (this%blade(ib)%wiP(i, j)%PC(1, 3) &
              - this%blade(ib)%wiP(i, j)%PC(1, 4))*0.25_dp

            call this%blade(ib)%wiP(i, j)%vr%assignP(1, &
              (/this%blade(ib)%wiP(i, j)%PC(1, 1) &
              + xshift(1), this%blade(ib)%wiP(i, j)%PC(2, 1), &
              this%blade(ib)%wiP(i, j)%PC(3, 1)/))
            call this%blade(ib)%wiP(i, j)%vr%assignP(2, &
              (/this%blade(ib)%wiP(i, j)%PC(1, 2) &
              + xshift(2), this%blade(ib)%wiP(i, j)%PC(2, 2), &
              this%blade(ib)%wiP(i, j)%PC(3, 2)/))
            call this%blade(ib)%wiP(i, j)%vr%assignP(3, &
              (/this%blade(ib)%wiP(i, j)%PC(1, 3) &
              + xshift(3), this%blade(ib)%wiP(i, j)%PC(2, 3), &
              this%blade(ib)%wiP(i, j)%PC(3, 3)/))
            call this%blade(ib)%wiP(i, j)%vr%assignP(4, &
              (/this%blade(ib)%wiP(i, j)%PC(1, 4) &
              + xshift(4), this%blade(ib)%wiP(i, j)%PC(2, 4), &
              this%blade(ib)%wiP(i, j)%PC(3, 4)/))
          enddo
        enddo

        ! Initialize vr coords of last row
        do j = 1, this%ns
          xshift(1) = (this%blade(ib)%wiP(this%nc, j)%PC(1, 2) &
            - this%blade(ib)%wiP(this%nc, j)%PC(1, 1))*0.25_dp  ! Shift x coord by dx/4
          xshift(2) = 0._dp
          xshift(3) = 0._dp
          xshift(4) = (this%blade(ib)%wiP(this%nc, j)%PC(1, 3) &
            - this%blade(ib)%wiP(this%nc, j)%PC(1, 4))*0.25_dp  ! Shift x coord by dx/4

          call this%blade(ib)%wiP(this%nc, j)%vr%assignP(1, &
            (/this%blade(ib)%wiP(i, j)%PC(1, 1) &
            + xshift(1), this%blade(ib)%wiP(i, j)%PC(2, 1), &
            this%blade(ib)%wiP(i, j)%PC(3, 1)/))
          call this%blade(ib)%wiP(this%nc, j)%vr%assignP(2, &
            (/this%blade(ib)%wiP(i, j)%PC(1, 2) &
            + xshift(2), this%blade(ib)%wiP(i, j)%PC(2, 2), &
            this%blade(ib)%wiP(i, j)%PC(3, 2)/))
          call this%blade(ib)%wiP(this%nc, j)%vr%assignP(3, &
            (/this%blade(ib)%wiP(i, j)%PC(1, 3) &
            + xshift(3), this%blade(ib)%wiP(i, j)%PC(2, 3), &
            this%blade(ib)%wiP(i, j)%PC(3, 3)/))
          call this%blade(ib)%wiP(this%nc, j)%vr%assignP(4, &
            (/this%blade(ib)%wiP(i, j)%PC(1, 4) &
            + xshift(4), this%blade(ib)%wiP(i, j)%PC(2, 4), &
            this%blade(ib)%wiP(i, j)%PC(3, 4)/))
        enddo
      endif

      ! Find dx and dy vectors
      do is = 1, this%ns
        do ic = 1, this%nc
          dx(ic, is) = norm2(this%blade(ib)%wiP(ic, is)%PC(:, 2) - this%blade(ib)%wiP(ic, is)%PC(:, 1))
          dy(ic, is) = norm2(this%blade(ib)%wiP(ic, is)%PC(:, 3) - this%blade(ib)%wiP(ic, is)%PC(:, 2))
        enddo
      enddo
      dx = abs(dx)
      dy = abs(dy)
      dxdymin = min(minval(dx), minval(dy))

      ! Shed last row of vortices
      if (abs(this%surfaceType) == 1) then
        if (abs(this%Omega) > eps) then
          velShed = min(0.05*abs(this%Omega)*norm2(this%blade(ib)% &
            wiP(this%nc, this%ns)%vr%vf(2)%fc(:, 1) - this%hubCoords), 0.125_dp*this%chord/dt)
        else
          velShed = 0.3_dp*norm2(-1._dp*this%velBody)
        endif
        do j = 1, this%ns
          call this%blade(ib)%wiP(this%nc, j)%vr%shiftdP(2, (/sign(1._dp,this%Omega)*velShed*dt, 0._dp, 0._dp/))
          call this%blade(ib)%wiP(this%nc, j)%vr%shiftdP(3, (/sign(1._dp,this%Omega)*velShed*dt, 0._dp, 0._dp/))
        enddo
      endif

      ! Initialize CP coords, nCap, panelArea and pivotLE
      if (abs(this%surfaceType) == 1) then
        do j = 1, this%ns
          do i = 1, this%nc
            call this%blade(ib)%wiP(i, j)%calcCP()
            call this%blade(ib)%wiP(i, j)%calcN()
            call this%blade(ib)%wiP(i, j)%calcTau()
            this%blade(ib)%wiP(i, j)%rHinge = length3d((this%blade(ib)%wiP(1, j)%pc(:, 1) &
              + this%blade(ib)%wiP(1, j)%pc(:, 4))*0.5_dp, this%blade(ib)%wiP(i, j)%CP)
            call this%blade(ib)%wiP(i, j)%calc_area()
            call this%blade(ib)%wiP(i, j)%calc_mean_dimensions()
          enddo
        enddo
      elseif(abs(this%surfaceType) == 2) then
        do j = 1, this%ns
          do i = 1, this%nc
            call this%blade(ib)%wiP(i, j)%calcCP(isTriangle=.true.)
            call this%blade(ib)%wiP(i, j)%calcN(isTriangle=.true.)
            call this%blade(ib)%wiP(i, j)%calcTau(isTriangle=.true.)
            ! this%blade(ib)%wiP(i, j)%rHinge = &
            ! & length3d((this%blade(ib)%wiP(1, j)%pc(:, 1) &
            ! & + this%blade(ib)%wiP(1, j)%pc(:, 4))*0.5_dp, &
            ! & this%blade(ib)%wiP(i, j)%CP)
            call this%blade(ib)%wiP(i, j)%calc_area()
            call this%blade(ib)%wiP(i, j)%calc_mean_dimensions()
          enddo
        enddo
      endif

      ! Compute sectional areas and chords
      call this%blade(ib)%calc_secArea()
      call this%blade(ib)%calc_secChord()

      ! Assign spanwise lift term switch to blades
      this%blade(ib)%spanwiseLiftSwitch = this%spanwiseLiftSwitch

      ! Invert half of tau vectors for symmetric or swept wings
      if (this%symmetricTau .eq. 1) then
        do j = 1, (this%ns/2)
          this%blade(ib)%secTauCapSpan(:, j) = -1._dp*yAxis
          do i = 1, this%nc
            this%blade(ib)%wiP(i, j)%tauCapSpan = -1._dp*this%blade(ib)%yAxis
            this%blade(ib)%wiP(i, j)%tauCapChord = this%blade(ib)%xAxis
          enddo
        enddo
        do j = (this%ns/2) + 1, this%ns
          do i = 1, this%nc
            this%blade(ib)%wiP(i, j)%tauCapSpan = this%blade(ib)%yAxis
            this%blade(ib)%wiP(i, j)%tauCapChord = this%blade(ib)%xAxis
          enddo
        enddo
      endif

      ! Inflow calculated at mid-chord
      secCPLoc = 0.5_dp
      this%blade(ib)%secCP = this%blade(ib)%getSecChordwiseLocations(secCPLoc)

      ! Initialize gamma
      this%blade(ib)%wiP%vr%gam = 0._dp
      this%blade(ib)%wiP%vr%skew = 0._dp
      this%blade(ib)%pivotLE = this%pivotLE

      ! Initialize wake age
      do i = 1, 4
        this%blade(ib)%wiP%vr%vf(i)%age = 0._dp
        this%blade(ib)%wiP%vr%vf(i)%ageAzimuthal = 0._dp
        this%blade(ib)%waN%vr%vf(i)%age = 0._dp
        this%blade(ib)%waN%vr%vf(i)%ageAzimuthal = 0._dp
      enddo

      this%blade(ib)%waF%vf%age = 0._dp
      this%blade(ib)%waF%vf%ageAzimuthal = 0._dp

      ! Initialize all core radius of wing vortices to zero
      do i = 1, 4
        this%blade(ib)%wiP%vr%vf(i)%rVc0 = &
          & min(this%spanwiseCore, dxdymin*0.1_dp)
      enddo

      if (min(this%spanwiseCore, dxdymin*0.1_dp) < eps) then
        print*, 'Warning: Wing vortex core radius set to zero'
      endif

      ! print*,'Wing vortex core radius set to ', &
      ! & min(this%spanwiseCore,dxdymin*0.1_dp)/this%chord,'times chord'

      ! Initialize spanwise vortex core radius for last row of wing
      ! to that of wake
      this%blade(ib)%wiP(this%nc, :)%vr%vf(2)%rVc0 = this%spanwiseCore

      ! Initialize all current core radius of wing vortices 
      ! to initial core radius
      do i = 1, 4
        this%blade(ib)%wiP%vr%vf(i)%rVc = this%blade(ib)%wiP%vr%vf(i)%rVc0
      enddo

      ! Verify CP is outside vortex core for boundary panels
      warnUser = .FALSE.
      if (this%blade(ib)%wiP(1, 1)%isCPinsidecore()) then
        print *, 'Warning: CP inside vortex core at panel LU'
        warnUser = .TRUE.
      endif
      if (this%blade(ib)%wiP(this%nc, 1)%isCPinsidecore()) then
        print *, 'Warning: CP inside vortex core at panel LB'
        warnUser = .TRUE.
      endif
      if (this%blade(ib)%wiP(1, this%ns)%isCPinsidecore()) then
        print *, 'Warning: CP inside vortex core at panel RU'
        warnUser = .TRUE.
      endif
      if (this%blade(ib)%wiP(this%nc, this%ns)%isCPinsidecore()) then
        print *, 'Warning: CP inside vortex core at panel RB'
        warnUser = .TRUE.
      endif
    enddo

    ! Removing sleep() for preventing ifort warnings
    ! if (warnUser .eqv. .TRUE.) then
    !   print *
    !   call sleep(2)
    ! endif

    ! Move rotor to hub coordinates
    do ib = 1, this%nb
      call this%blade(ib)%move(this%hubCoords-this%fromCoords)
    enddo

    ! Set Coning angle
    do ib = 1, this%nb
      call this%blade(ib)%rot_axis(this%coningAngle, &
        & xAxis, (/0._dp, 0._dp, 0._dp/))
    enddo

    ! Rotate remaining blades to their positions
    ! Rotate blades for multi-bladed rotors
    do ib = 2, this%nb
      bladeOffset = twoPi/this%nb*(ib - 1)
      call this%blade(ib)%rot_axis(bladeOffset, this%shaftAxis, this%hubCoords)
    enddo

    ! Rotate rotor by phi,theta,psi about CG
    call this%rot_pts(this%pts, this%cgCoords, 1)

    ! Rotate rotor by psiStart
    call this%rot_advance(this%psiStart, nopitch=.true.)

    ! Compute denominators for non-dimensionalisation
    if (abs(this%Omega) .gt. eps) then
      if (this%propConvention .eq. 0) then
        ! Rotary-wing
        this%nonDimforceDenominator = density*(pi*this%radius**2._dp)* &
          (this%radius*this%Omega)**2._dp
      else
        ! Propeller
        this%nonDimforceDenominator = density*(this%Omega/twoPi)**2._dp* &
          (2._dp*this%radius)**4._dp
      endif
    else
      ! Fixed-wing
      this%nonDimforceDenominator = 0.5_dp*density &
        & *(this%radius*(1._dp - this%root_cut)* &
        this%chord)*(dot_product(this%velBody, this%velBody))
    endif

    ! Allocate and assign section airfoils
    if (this%nAirfoils .gt. 0) then
      do ib = 1, this%nb
        allocate (this%blade(ib)%C81(this%nAirfoils))
        do i = 1, this%nAirfoils
          if (this%airfoilFile(i) (1:1) .ne. '0') &
            call this%blade(ib)%C81(i)%readfile(trim(this%airfoilFile(i)))
        enddo
        allocate (this%blade(ib)%airfoilSectionLimit(this%nAirfoils))
        allocate (this%blade(ib)%alpha0(this%nAirfoils))
        this%blade(ib)%airfoilSectionLimit = this%airfoilSectionLimit
        do i = 1, this%nAirfoils
          this%blade(ib)%alpha0(i) = this%alpha0(i)
        enddo

        ! Assign airfoil numbers for each section
        if (this%nAirfoils .eq. 1) then
          this%blade(ib)%airfoilNo(:) = 1
        else
          leftTipCP = this%blade(ib)%wiP(1, 1)%PC(:, 4)*(1._dp - secCPLoc) &
            + this%blade(ib)%wiP(this%nc, 1)%PC(:, 3)*secCPLoc
          do is = 1, this%ns
            ! Warning: This will break if wing is centered about X-Z plane and
            ! full span length is used for reference length
            rbyR = abs(dot_product(this%blade(ib)%secCP(:, is) - leftTipCP, this%blade(ib)%xAxis)) &
              + this%root_cut
            if (rbyR .gt. 1._dp) error stop 'ERROR: r/R value is greater than 1 in airfoil selection'
            do i = 1, this%nAirfoils
              if (this%airfoilSectionLimit(i) .ge. rbyR) then
                this%blade(ib)%airfoilNo(is) = i
                exit
              endif
            enddo
          enddo
        endif
      enddo
    endif

    ! Allocate vars required for wake convection
    ! on the basis of finite diff scheme
    do ib = 1, this%nb
      allocate(this%blade(ib)%velNwake(3, this%nNwake, this%ns + 1))
      allocate(this%blade(ib)%velFwake(3, this%nFwake))
      this%blade(ib)%velNwake = 0._dp
      this%blade(ib)%velFwake = 0._dp

      select case (switches%fdScheme)
      case (0)
        ! Do nothing
      case (1)
        allocate (this%blade(ib)%waNPredicted(this%nNwake, this%ns))
        allocate (this%blade(ib)%velNwakePredicted(3, this%nNwake, this%ns + 1))
        this%blade(ib)%velNwakePredicted = 0._dp

        allocate (this%blade(ib)%waFPredicted(this%nFwake))
        allocate (this%blade(ib)%velFwakePredicted(3, this%nFwake))
        this%blade(ib)%velFwakePredicted = 0._dp
      case (2)
        allocate (this%blade(ib)%velNwake1(3, this%nNwake, this%ns + 1))
        allocate (this%blade(ib)%velNwakeStep(3, this%nNwake, this%ns + 1))
        this%blade(ib)%velNwake1 = 0._dp
        this%blade(ib)%velNwakeStep = 0._dp

        allocate (this%blade(ib)%velFwake1(3, this%nFwake))
        allocate (this%blade(ib)%velFwakeStep(3, this%nFwake))
        this%blade(ib)%velFwake1 = 0._dp
        this%blade(ib)%velFwakeStep = 0._dp

      case (3)
        allocate (this%blade(ib)%waNPredicted(this%nNwake, this%ns))
        allocate (this%blade(ib)%velNwake1(3, this%nNwake, this%ns + 1))
        allocate (this%blade(ib)%velNwakePredicted(3, this%nNwake, this%ns + 1))
        allocate (this%blade(ib)%velNwakeStep(3, this%nNwake, this%ns + 1))
        this%blade(ib)%velNwake1 = 0._dp
        this%blade(ib)%velNwakePredicted = 0._dp
        this%blade(ib)%velNwakeStep = 0._dp

        allocate (this%blade(ib)%waFPredicted(this%nFwake))
        allocate (this%blade(ib)%velFwake1(3, this%nFwake))
        allocate (this%blade(ib)%velFwakePredicted(3, this%nFwake))
        allocate (this%blade(ib)%velFwakeStep(3, this%nFwake))
        this%blade(ib)%velFwake1 = 0._dp
        this%blade(ib)%velFwakePredicted = 0._dp
        this%blade(ib)%velFwakeStep = 0._dp
      case (4)
        allocate (this%blade(ib)%waNPredicted(this%nNwake, this%ns))
        allocate (this%blade(ib)%velNwake1(3, this%nNwake, this%ns + 1))
        allocate (this%blade(ib)%velNwake2(3, this%nNwake, this%ns + 1))
        allocate (this%blade(ib)%velNwakePredicted(3, this%nNwake, this%ns + 1))
        allocate (this%blade(ib)%velNwakeStep(3, this%nNwake, this%ns + 1))
        this%blade(ib)%velNwake1 = 0._dp
        this%blade(ib)%velNwake2 = 0._dp
        this%blade(ib)%velNwakePredicted = 0._dp
        this%blade(ib)%velNwakeStep = 0._dp

        allocate (this%blade(ib)%waFPredicted(this%nFwake))
        allocate (this%blade(ib)%velFwake1(3, this%nFwake))
        allocate (this%blade(ib)%velFwake2(3, this%nFwake))
        allocate (this%blade(ib)%velFwakePredicted(3, this%nFwake))
        allocate (this%blade(ib)%velFwakeStep(3, this%nFwake))
        this%blade(ib)%velFwake1 = 0._dp
        this%blade(ib)%velFwake2 = 0._dp
        this%blade(ib)%velFwakePredicted = 0._dp
        this%blade(ib)%velFwakeStep = 0._dp
      case (5)
        allocate (this%blade(ib)%waNPredicted(this%nNwake, this%ns))
        allocate (this%blade(ib)%velNwake1(3, this%nNwake, this%ns + 1))
        allocate (this%blade(ib)%velNwake2(3, this%nNwake, this%ns + 1))
        allocate (this%blade(ib)%velNwake3(3, this%nNwake, this%ns + 1))
        allocate (this%blade(ib)%velNwakePredicted(3, this%nNwake, this%ns + 1))
        allocate (this%blade(ib)%velNwakeStep(3, this%nNwake, this%ns + 1))
        this%blade(ib)%velNwake1 = 0._dp
        this%blade(ib)%velNwake2 = 0._dp
        this%blade(ib)%velNwake3 = 0._dp
        this%blade(ib)%velNwakePredicted = 0._dp
        this%blade(ib)%velNwakeStep = 0._dp

        allocate (this%blade(ib)%waFPredicted(this%nFwake))
        allocate (this%blade(ib)%velFwake1(3, this%nFwake))
        allocate (this%blade(ib)%velFwake2(3, this%nFwake))
        allocate (this%blade(ib)%velFwake3(3, this%nFwake))
        allocate (this%blade(ib)%velFwakePredicted(3, this%nFwake))
        allocate (this%blade(ib)%velFwakeStep(3, this%nFwake))
        this%blade(ib)%velFwake1 = 0._dp
        this%blade(ib)%velFwake2 = 0._dp
        this%blade(ib)%velFwake3 = 0._dp
        this%blade(ib)%velFwakePredicted = 0._dp
        this%blade(ib)%velFwakeStep = 0._dp
      end select
    enddo

    ! Wake initialization
    ! Assign core_radius to mid vortices
    do ib = 1, this%nb
      do i = 2, 4, 2
        this%blade(ib)%waN%vr%vf(i)%rVc0 = this%spanwiseCore
        this%blade(ib)%waN%vr%vf(i)%rVc = this%spanwiseCore
      enddo

      this%blade(ib)%waN%vr%gam = 0._dp
      this%blade(ib)%waF%gam = 0._dp

      ! Assign core_radius to tip vortices
      do j = 1, this%ns
        do i = 1, this%nNwake
          this%blade(ib)%waN(i, j)%vr%vf(1)%rVc0 = this%streamwiseCoreVec(j)
          this%blade(ib)%waN(i, j)%vr%vf(1)%rVc = this%streamwiseCoreVec(j)
        enddo
      enddo

      do j = 1, this%ns
        do i = 1, this%nNwake
          this%blade(ib)%waN(i, j)%vr%vf(3)%rVc0 = this%streamwiseCoreVec(j + 1)
          this%blade(ib)%waN(i, j)%vr%vf(3)%rVc = this%streamwiseCoreVec(j + 1)
        enddo
      enddo

      !do i=1,this%nFwake
      this%blade(ib)%waF%vf%rVc0 = this%streamwiseCoreVec(this%ns + 1)
      this%blade(ib)%waF%vf%rVc = this%streamwiseCoreVec(this%ns + 1)
      !enddo

    enddo

    ! Compute direction of wind frame forces (-body)
    ! Assuming sideslip is not present
    if ((norm2(this%dragUnitVec) .le. eps) &
      .and. (norm2(this%sideUnitVec) .le. eps) &
      .and. (norm2(this%liftUnitVec) .le. eps)) then
      if (abs(this%Omega) .le. eps) then
        if (abs(this%velBody(1)) .gt. eps) then  ! v is assumed zero
          this%dragUnitVec = -1._dp*unitVec((/this%velBody(1), 0._dp, this%velBody(3)/))
          this%sideUnitVec = yAxis
        else  ! u is assumed zero
          this%dragUnitVec = -1._dp*unitVec((/0._dp, this%velBody(2), this%velBody(3)/))
          this%sideUnitVec = xAxis
        endif
        this%liftUnitVec = cross_product(this%dragUnitVec, this%sideUnitVec)
      else
        ! Drag along forward velocity direction
        if (abs(this%velBody(1)) .gt. eps) then
          this%dragUnitVec = -1._dp*unitVec((/this%velBody(1), 0._dp, this%velBody(3)/))
          this%sideUnitVec = yAxis
        elseif (abs(this%velBody(2)) .gt. eps) then
          this%dragUnitVec = -1._dp*unitVec((/this%velBody(1), 0._dp, this%velBody(3)/))
          this%sideUnitVec = xAxis
        else
          this%sideUnitVec = yAxis
          this%dragUnitVec = cross_product(this%sideUnitVec, this%shaftAxis)
        endif
        this%liftUnitVec = this%shaftAxis
      endif
    endif

  end subroutine rotor_init

  subroutine rotor_deinit(this, switches)
    ! Deinitialise rotor variables
  class(rotor_class) :: this
    type(switches_class), intent(in) :: switches
    integer :: ib
    ! Deallocate variables
    do ib = 1, this%nb
      deallocate (this%blade(ib)%velNwake)
      deallocate (this%blade(ib)%velFwake)

      select case (switches%fdScheme)
      case (0)
        ! Nothing to deallocate
      case (1)
        deallocate (this%blade(ib)%waNPredicted)
        deallocate (this%blade(ib)%velNwakePredicted)

        deallocate (this%blade(ib)%waFPredicted)
        deallocate (this%blade(ib)%velFwakePredicted)
      case (2)
        deallocate (this%blade(ib)%velNwake1)
        deallocate (this%blade(ib)%velNwakeStep)

        deallocate (this%blade(ib)%velFwake1)
        deallocate (this%blade(ib)%velFwakeStep)
      case (3)
        deallocate (this%blade(ib)%waNPredicted)
        deallocate (this%blade(ib)%velNwake1)
        deallocate (this%blade(ib)%velNwakeStep)

        deallocate (this%blade(ib)%waFPredicted)
        deallocate (this%blade(ib)%velFwake1)
        deallocate (this%blade(ib)%velFwakeStep)
      case (4)
        deallocate (this%blade(ib)%waNPredicted)
        deallocate (this%blade(ib)%velNwake1)
        deallocate (this%blade(ib)%velNwake2)
        deallocate (this%blade(ib)%velNwakeStep)

        deallocate (this%blade(ib)%waFPredicted)
        deallocate (this%blade(ib)%velFwake1)
        deallocate (this%blade(ib)%velFwake2)
        deallocate (this%blade(ib)%velFwakeStep)
      case (5)
        deallocate (this%blade(ib)%waNPredicted)
        deallocate (this%blade(ib)%velNwake1)
        deallocate (this%blade(ib)%velNwake2)
        deallocate (this%blade(ib)%velNwake3)
        deallocate (this%blade(ib)%velNwakePredicted)
        deallocate (this%blade(ib)%velNwakeStep)

        deallocate (this%blade(ib)%waFPredicted)
        deallocate (this%blade(ib)%velFwake1)
        deallocate (this%blade(ib)%velFwake2)
        deallocate (this%blade(ib)%velFwake3)
        deallocate (this%blade(ib)%velFwakePredicted)
        deallocate (this%blade(ib)%velFwakeStep)
      end select
    enddo

  end subroutine rotor_deinit

  subroutine rotor_plot3dtoblade(this, PLOT3Dfilename)
    ! Read blade geometry from PLOT3D formatted file
  class(rotor_class) :: this
    character(len=*), intent(in) :: PLOT3Dfilename
    integer :: nx, ny, nz
    real(dp), allocatable, dimension(:, :, :) :: grid
    integer :: i, j, ic, is, ib
    logical :: dataMismatch

    open (unit=10, file=PLOT3Dfilename, action='read', status='old')
    read (10, *) nx, ny, nz

    ! Verify with rotor parameters
    dataMismatch = .FALSE.
    if (nz .gt. 1) dataMismatch = .TRUE.
    if (nx .ne. (this%nc + 1)) dataMismatch = .TRUE.
    if (ny .ne. (this%ns + 1)) dataMismatch = .TRUE.

    if (dataMismatch) then
      close (10)
      print *, 'PLOT3D file   (nx,ny)    =', nx, ny
      print *, 'ROTOR file    (nc+1,ns+1)=', this%nc + 1, this%ns + 1
      error stop 'ERROR: Wrong or conflicting data in PLOT3D file '
    else
      allocate (grid(3, nx, ny))
      read (10, *) &
        ((grid(1, i, j), i=1, nx), j=1, ny), &
        ((grid(2, i, j), i=1, nx), j=1, ny), &
        ((grid(3, i, j), i=1, nx), j=1, ny)
      close (10)
    endif

    ! Assign to blades
    do ib = 1, this%nb
      do is = 1, this%ns
        do ic = 1, this%nc
          call this%blade(ib)%wiP(ic, is)%assignP(1, grid(:, ic, is))
          call this%blade(ib)%wiP(ic, is)%assignP(2, grid(:, ic + 1, is))
          call this%blade(ib)%wiP(ic, is)%assignP(3, grid(:, ic + 1, is + 1))
          call this%blade(ib)%wiP(ic, is)%assignP(4, grid(:, ic, is + 1))
        enddo
      enddo
    enddo

  end subroutine rotor_plot3dtoblade

  subroutine rotor_stltoblade(this, stlfilename)
    !! Read ASCII stl file for non-lifting surface geometry
  class(rotor_class) :: this
    character(len=*), intent(in) :: stlfilename
    character(len=5) :: facet
    character(len=6) :: normal, vertex
    integer :: i

    open(unit=10, file=stlfilename, action='read', status='old')
    read(10, *) ! solid
    do i = 1, this%nc
      read(10, *) facet, normal, this%blade(1)%wiP(i, 1)%nCap
      read(10, *) ! outer loop
      read(10, *) vertex, this%blade(1)%wiP(i, 1)%PC(:, 1)
      read(10, *) vertex, this%blade(1)%wiP(i, 1)%PC(:, 2)
      read(10, *) vertex, this%blade(1)%wiP(i, 1)%PC(:, 3)
      read(10, *) ! endloop
      read(10, *) ! endfacet
    enddo
    read(10, *) ! endsolid
    close(10)

    ! Compute vortex ring coordinates
    !$omp parallel do
    do i = 1, this%nc
      this%blade(1)%wiP(i, 1)%PC(:, 4) = this%blade(1)%wiP(i, 1)%PC(:, 3)
      call this%blade(1)%wiP(i, 1)%vr%assignP(1, &
        & this%blade(1)%wiP(i, 1)%PC(:, 1))
      call this%blade(1)%wiP(i, 1)%vr%assignP(2, &
        & this%blade(1)%wiP(i, 1)%PC(:, 2))
      call this%blade(1)%wiP(i, 1)%vr%assignP(3, &
        & this%blade(1)%wiP(i, 1)%PC(:, 3))
      call this%blade(1)%wiP(i, 1)%vr%assignP(4, &
        & this%blade(1)%wiP(i, 1)%PC(:, 4))
    enddo
    !$omp end parallel do
  end subroutine rotor_stltoblade

  function rotor_getCamber(this, x, y)
    !! Get z coordinate on wing from x, y values
    ! File format:
    ! <num of lines>
    ! <x1>  <z1>
    ! <x2>  <z2>
    ! ...
    use libMath, only: interp1d
  class(rotor_class) :: this
    real(dp), intent(in), dimension(:) :: x, y
    real(dp), dimension(size(x), size(y)) :: rotor_getCamber
    real(dp), dimension(5000, this%nCamberFiles) :: xCamber, zCamber
    real(dp) :: chord, span
    integer, dimension(this%nCamberFiles) :: nPts
    integer :: iSect, fNum, i, j

    chord = x(size(x)) - x(1)
    span = y(size(y)) - y(1)

    ! Read camber data
    do j = 1, this%nCamberFiles
      if (this%camberFile(j)(1:1) == '0') then
        nPts(j) = 2
        xCamber(1:2, j) = (/0._dp, 1._dp/)
        zCamber(1:2, j) = 0._dp
      else
        open(unit=15, file=this%camberFile(j))
        read(15, *) nPts(j)
        do i = 1, nPts(j)
          read(15, *) xCamber(i, j), zCamber(i, j)
        enddo
        close(15)
        zCamber(:, j) = zCamber(:, j)/xCamber(nPts(j), j)
        xCamber(:, j) = xCamber(:, j)/xCamber(nPts(j), j)
      endif
    enddo

    do j = 1, size(y)
      ! Check which sectionLimit y value comes under
      ! to determine airfoil file to use
      do iSect = 1, this%nCamberFiles
        if (y(j)/span <= this%camberSectionLimit(iSect)) then
          fNum = iSect
          exit
        endif
      enddo

      do i = 1, size(x)
        rotor_getCamber(i, j) = interp1d((x(i)-x(1))/chord, &
          & xCamber(1:nPts(fNum), fNum), zCamber(1:nPts(fNum), fNum), 2)
      enddo
    enddo
  end function rotor_getCamber

  function rotor_gettheta(this, psi, ib)
    ! Get pitch angle corresponding to blade azimuthal location
  class(rotor_class) :: this
    real(dp), intent(in) :: psi
    integer, intent(in) :: ib
    real(dp) :: rotor_gettheta
    real(dp) :: bladeOffset

    bladeOffset = twoPi/this%nb*(ib - 1)
    rotor_gettheta = this%controlPitch(1) &
      + this%controlPitch(2)*cos(psi + bladeOffset) &
      + this%controlPitch(3)*sin(psi + bladeOffset)

  end function rotor_gettheta

  function rotor_getthetadot(this, psi, ib)
  class(rotor_class) :: this
    real(dp), intent(in) :: psi
    integer, intent(in) :: ib
    real(dp) :: rotor_getthetadot
    real(dp) :: bladeOffset

    bladeOffset = twoPi/this%nb*(ib - 1)
    rotor_getthetadot = -this%controlPitch(2)*sin(psi + bladeOffset) &
      + this%controlPitch(3)*cos(psi + bladeOffset)

  end function rotor_getthetadot

  subroutine rotor_calcAIC(this)
    ! Compute AIC matrix for rotor
    use libMath, only: inv2
  class(rotor_class), intent(inout) :: this
    integer :: ib, jblade, is, ic, i, j, row, col
    real(dp), dimension(3) :: vec_dummy

    ! Influence Coefficient Matrix
    !$omp parallel do private(row, col) collapse(3)
    do ib = 1, this%nb
      do is = 1, this%ns      ! Collocation point loop
        do ic = 1, this%nc
          row = ic + this%nc*(is - 1) + this%ns*this%nc*(ib - 1)

          do jblade = 1, this%nb
            do j = 1, this%ns       ! Vortex ring loop
              do i = 1, this%nc
                col = i + this%nc*(j - 1) + this%ns*this%nc*(jblade - 1)
                vec_dummy = this%blade(jblade)%wiP(i, j)%vr%vind(this%blade(ib)%wiP(ic, is)%CP)
                this%AIC(row, col) = dot_product(vec_dummy, this%blade(ib)%wiP(ic, is)%nCap)
              enddo
            enddo
          enddo

        enddo
      enddo
    enddo
    !$omp end parallel do
    this%AIC_inv = inv2(this%AIC)
  end subroutine rotor_calcAIC

  subroutine rotor_map_gam(this)
    ! Map gam from vector to matrix format
  class(rotor_class), intent(inout) :: this
    integer :: ib
    if (this%imposeAxisymmetry == 0) then
      do ib = 1, this%nb
        this%blade(ib)%wiP%vr%gam &
          = reshape(this%gamVec(1+this%nc*this%ns*(ib-1):this%nc*this%ns*ib), &
          & (/this%nc, this%ns/))
      enddo
    else
      this%blade(1)%wiP%vr%gam &
        = reshape(this%gamVec(1:this%nc*this%ns), &
        & (/this%nc, this%ns/))
      do ib = 2, this%nb
        this%blade(ib)%wiP%vr%gam = this%blade(1)%wiP%vr%gam
      enddo
    endif
  end subroutine rotor_map_gam

  !-----+------------------+-----|
  ! -+- | Motion Functions | -+- |
  !-----+------------------+-----|

  subroutine rotor_move(this, dshift)
  class(rotor_class) :: this
    real(dp), intent(in), dimension(3) :: dshift

    integer :: ib

    do ib = 1, this%nb
      call this%blade(ib)%move(dshift)
    enddo
    this%hubCoords = this%hubCoords + dshift
    this%cgCoords = this%cgCoords + dshift

  end subroutine rotor_move

  subroutine rotor_rot_pts(this, pts, origin, order)
    ! Rotate using pts => phi theta psi
    use libMath, only: Tbg, Tgb
  class(rotor_class), intent(inout) :: this
    real(dp), dimension(3), intent(in) :: pts    ! pts => phi,theta,psi
    real(dp), dimension(3), intent(in) :: origin ! rotation about
    integer, intent(in) :: order    ! [1]gb & +ve theta , [2]bg & -ve theta
    integer :: ib
    real(dp), dimension(3, 3) :: TMat

    select case (order)
    case (2)
      TMat = Tbg((/cos(pts(1)), sin(pts(1))/), &
        (/cos(pts(2)), sin(pts(2))/), &
        (/cos(pts(3)), sin(pts(3))/))
    case (1)
      TMat = Tgb((/cos(pts(1)), sin(pts(1))/), &
        (/cos(pts(2)), sin(pts(2))/), &
        (/cos(pts(3)), sin(pts(3))/))
    case default
      error stop 'ERROR: wrong option for order'
    end select

    do ib = 1, this%nb
      call this%blade(ib)%rot_pts(pts, origin, order)
    enddo

    this%shaftAxis = matmul(TMat, this%shaftAxis)

    this%hubCoords = matmul(TMat, this%hubCoords-origin)+origin
    this%cgCoords = matmul(TMat, this%cgCoords-origin)+origin

  end subroutine rotor_rot_pts

  subroutine rotor_rot_advance(this, dpsi, nopitch)
    ! Rotate rotos by dpsi angle about axis
  class(rotor_class), intent(inout) :: this
    real(dp), intent(in) :: dpsi
    logical, optional ::  nopitch
    integer :: ib
    real(dp) :: dtheta

    this%psi = this%psi + dpsi
    do ib = 1, this%nb
      call this%blade(ib)%rot_axis(dpsi, this%shaftAxis, this%hubCoords)
      this%blade(ib)%psi = this%blade(ib)%psi + dpsi
      if (.not. present(nopitch)) then
        dtheta = this%gettheta(this%psi, ib) - this%blade(ib)%theta
        call this%blade(ib)%rot_pitch(dtheta)
        this%blade(ib)%theta = this%gettheta(this%psi, ib)
      elseif (nopitch .eqv. .false.) then
        dtheta = this%gettheta(this%psi, ib) - this%blade(ib)%theta
        call this%blade(ib)%rot_pitch(dtheta)
        this%blade(ib)%theta = this%gettheta(this%psi, ib)
      endif
    enddo

  end subroutine rotor_rot_advance

  !-----+---------------------------+-----|
  ! -+- | Wake Convection Functions | -+- |
  !-----+---------------------------+-----|

  subroutine rotor_assignshed(this, edge)
    ! Assign coordinates to first rowNear of wake from last row of blade
  class(rotor_class), intent(inout) :: this
    character(len=2), intent(in) :: edge
    integer :: i, ib

    select case (edge)
    case ('LE')    ! assign to LE
      do ib = 1, this%nb
        do i = 1, this%ns
          call this%blade(ib)%waN(this%rowNear, i)%vr%assignP(1, this%blade(ib)%wiP(this%nc, i)%vr%vf(2)%fc(:, 1))
          call this%blade(ib)%waN(this%rowNear, i)%vr%assignP(4, this%blade(ib)%wiP(this%nc, i)%vr%vf(3)%fc(:, 1))
          call this%blade(ib)%waN(this%rowNear, i)%vr%calclength(.TRUE.)    ! TRUE => record original length
        enddo
        this%blade(ib)%waN(this%rowNear, :)%vr%gam = this%blade(ib)%wiP(this%nc, :)%vr%gam

      enddo
    case ('TE')    ! assign to next row's TE
      do ib = 1, this%nb
        do i = 1, this%ns
          call this%blade(ib)%waN(max(this%rowNear - 1, 1), i)%vr%assignP(2, this%blade(ib)%wiP(this%nc, i)%vr%vf(2)%fc(:, 1))
          call this%blade(ib)%waN(max(this%rowNear - 1, 1), i)%vr%assignP(3, this%blade(ib)%wiP(this%nc, i)%vr%vf(3)%fc(:, 1))
        enddo
      enddo
    case default
      error stop 'ERROR: Wrong option for edge'
    end select

  end subroutine rotor_assignshed

  !-----+----------------------------+-----|
  ! -+- | Wake Dissipation Functions | -+- |
  !-----+----------------------------+-----|

  subroutine rotor_age_wake(this, dt)
    ! Update age of wake filaments
  class(rotor_class), intent(inout) :: this
    real(dp), intent(in) :: dt
    integer :: ib, ifil
    do ib = 1, this%nb
      do ifil = 1, 4
        this%blade(ib)%waN(this%rowNear:this%nNwake, :)% &
          & vr%vf(ifil)%age = &
          & this%blade(ib)%waN(this%rowNear:this%nNwake, :)% &
          & vr%vf(ifil)%age + dt

        this%blade(ib)%waN(this%rowNear:this%nNwake, :)% &
          & vr%vf(ifil)%ageAzimuthal = &
          & this%blade(ib)%waN(this%rowNear:this%nNwake, :)% &
          & vr%vf(ifil)%ageAzimuthal + dt*this%omegaSlow
      enddo
      this%blade(ib)%waF(this%rowFar:this%nFwake)%vf%age = &
        this%blade(ib)%waF(this%rowFar:this%nFwake)%vf%age + dt
      this%blade(ib)%waF(this%rowFar:this%nFwake)%vf%ageAzimuthal = &
        this%blade(ib)%waF(this%rowFar:this%nFwake)%vf%ageAzimuthal + &
        & dt*this%omegaSlow
    enddo
  end subroutine rotor_age_wake

  subroutine rotor_dissipate_wake(this, dt, kinematicViscosity)
  class(rotor_class), intent(inout) :: this
    real(dp), intent(in) :: dt, kinematicViscosity
    real(dp) :: oseenParameter
    integer :: ib, ic, is
    oseenParameter = 1.2564_dp

    ! Dissipate near wake
    do ib = 1, this%nb
      do is = 1, this%ns
        !$omp parallel do
        do ic = this%rowNear, this%nNwake
          this%blade(ib)%waN(ic, is)%vr%vf(1)%rVc = &
            & sqrt(this%blade(ib)%waN(ic, is)%vr%vf(1)%rVc**2._dp &
            & + 4._dp*oseenParameter*this%apparentViscCoeff*kinematicViscosity*dt)
          this%blade(ib)%waN(ic, is)%vr%vf(3)%rVc = &
            & this%blade(ib)%waN(ic, is)%vr%vf(1)%rVc
          ! Decay wake
          call this%blade(ib)%waN(ic, is)%vr%decay(dt, this%decayCoeff)
        enddo
        !$omp end parallel do

        ! To maintain consistency of rVc in overlapping filaments
        !$omp parallel do
        do ic = this%rowNear, this%nNwake
          this%blade(ib)%waN(ic, is)%vr%vf(2)%rVc = sqrt(this%blade(ib)%waN(ic, is)%vr%vf(2)%rVc**2._dp &
            + 4._dp*oseenParameter*this%apparentViscCoeff*kinematicViscosity*dt)
        enddo
        !$omp end parallel do

        if (this%rowNear .ne. this%nNwake) then
          !$omp parallel do
          do ic = this%rowNear + 1, this%nNwake
            this%blade(ib)%waN(ic, is)%vr%vf(4)%rVc = this%blade(ib)%waN(ic - 1, is)%vr%vf(2)%rVc
          enddo
          !$omp end parallel do
        endif
      enddo

      ! Dissipate far wake if present
      !$omp parallel do
      do ic = this%rowFar, this%nFwake
        this%blade(ib)%waF(ic)%vf%rVc = &
          & sqrt(this%blade(ib)%waF(ic)%vf%rVc**2._dp &
          & + 4._dp*oseenParameter*this%apparentViscCoeff* &
          & kinematicViscosity*dt)
          ! Decay wake
          call this%blade(ib)%waF(ic)%decay(dt, this%decayCoeff)
      enddo
      !$omp end parallel do
    enddo

  end subroutine rotor_dissipate_wake

  subroutine rotor_strain_wake(this)
  class(rotor_class), intent(inout) :: this
    integer :: i, ib

    do ib = 1, this%nb
      !$omp parallel do
      do i = this%rowFar, this%nFwake
        call this%blade(ib)%waF(i)%vf%calclength(.FALSE.)    ! Update current length
        call this%blade(ib)%waF(i)%vf%strain()
      enddo
      !$omp end parallel do
    enddo
  end subroutine rotor_strain_wake

  function rotor_vind_bywing(this, P)
    ! Compute induced velocity by all wing vortices at P
  class(rotor_class), intent(inout) :: this
    real(dp), intent(in), dimension(3) :: P
    real(dp), dimension(3) :: rotor_vind_bywing
    integer :: ib

    rotor_vind_bywing = 0._dp
    if (abs(this%surfaceType) == 1) then
      do ib = 1, this%nb
        rotor_vind_bywing = rotor_vind_bywing &
          & + this%blade(ib)%vind_bywing(P)
      enddo
    elseif (abs(this%surfaceType) == 2) then
      do ib = 1, this%nb
        rotor_vind_bywing = rotor_vind_bywing &
          & + this%blade(ib)%vindSource_bywing(P)
      enddo
    endif
  end function rotor_vind_bywing

  function rotor_vind_bywing_boundVortices(this, P)
    ! Compute induced velocity by bound vortices at P
  class(rotor_class), intent(inout) :: this
    real(dp), intent(in), dimension(3) :: P
    real(dp), dimension(3) :: rotor_vind_bywing_boundVortices
    integer :: ib

    rotor_vind_bywing_boundVortices = 0._dp
    do ib = 1, this%nb
      rotor_vind_bywing_boundVortices = rotor_vind_bywing_boundVortices &
        & + this%blade(ib)%vind_bywing_boundVortices(P)
    enddo
  end function rotor_vind_bywing_boundVortices

  function rotor_vind_bywake(this, P, optionalChar)
    ! Compute induced velocity by wake vortices at P
  class(rotor_class), intent(in) :: this
    real(dp), intent(in), dimension(3) :: P
    character(len=1), optional :: optionalChar
    real(dp), dimension(3) :: rotor_vind_bywake
    integer :: ib

    rotor_vind_bywake = 0._dp
    if (.not. present(optionalChar)) then
      do ib = 1, this%nb
        rotor_vind_bywake = rotor_vind_bywake + this%blade(ib)%vind_bywake(this%rowNear, this%rowFar, P)
      enddo
    elseif ((optionalChar .eq. 'P') .or. (optionalChar .eq. 'p')) then
      do ib = 1, this%nb
        rotor_vind_bywake = rotor_vind_bywake + this%blade(ib)%vind_bywake(this%rowNear, this%rowFar, P, 'P')
      enddo
    else
      error stop 'ERROR: Wrong character flag for rotor_vind_bywake()'
    endif
  end function rotor_vind_bywake

  subroutine rotor_shiftwake(this)
    ! Shift wake locations on rollup

  class(rotor_class), intent(inout) :: this
    integer :: ib, i

    do ib = 1, this%nb
      do i = this%nNwake, 2, -1
        this%blade(ib)%waN(i, :) = this%blade(ib)%waN(i - 1, :)
      enddo

      ! Wake age of first row has to be set to zero
      do i = 1, 4
        this%blade(ib)%waN(1, :)%vr%vf(i)%age = 0._dp
      enddo
    enddo

  end subroutine rotor_shiftwake

  subroutine rotor_rollup(this)
    !    2
    !    |    ^ Upstream
    !    |    |
    !    |
    !    1

  class(rotor_class), intent(inout) :: this
    integer :: ib, ispan, rowFarNext
    real(dp), dimension(3) :: centroidLE, centroidTE
    real(dp) :: gamRollup, ageRollup, radiusRollup, gamSum

    rowFarNext = this%rowFar - 1    ! Rollup the vortex filament of 'next' row

    do ib = 1, this%nb
      gamRollup = this%blade(ib)%waN(this%nNwake, this%ns)%vr%gam
      centroidLE = 0._dp
      centroidTE = 0._dp
      radiusRollup = 0._dp
      gamSum = 0._dp

      do ispan = this%rollupStart, this%rollupEnd
        ! Find centroid LE and TE
        centroidLE = centroidLE + this%blade(ib)%waN(this%nNwake, ispan)%vr%vf(4)%fc(:, 1)* &
          this%blade(ib)%waN(this%nNwake, ispan)%vr%gam
        centroidTE = centroidTE + this%blade(ib)%waN(this%nNwake, ispan)%vr%vf(3)%fc(:, 1)* &
          this%blade(ib)%waN(this%nNwake, ispan)%vr%gam
        gamSum = gamSum + this%blade(ib)%waN(this%nNwake, ispan)%vr%gam

        ! Assign gamRollup and radiusRollup from last row to wake filament gamma
        ! Compute gamRollup
        if (sign(1._dp, this%Omega*this%controlPitch(1)) > eps) then    ! +ve Omega or zero Omega with +ve pitch
          if (this%blade(ib)%waN(this%nNwake, ispan)%vr%gam < gamRollup) then    ! '<' because of negative gamma
            gamRollup = this%blade(ib)%waN(this%nNwake, ispan)%vr%gam
          endif
        else    ! one of Omega or pitch is negative
          if (this%blade(ib)%waN(this%nNwake, ispan)%vr%gam > gamRollup) then    ! '>' because of positive gamma
            gamRollup = this%blade(ib)%waN(this%nNwake, ispan)%vr%gam
          endif
        endif

        ! Compute radiusRollup
        radiusRollup = radiusRollup + this%blade(ib)%waN(this%nNwake, ispan)%vr%vf(3)%rVc* &
          this%blade(ib)%waN(this%nNwake, ispan)%vr%gam
      enddo

      ageRollup = this%blade(ib)%waN(this%nNwake, this%ns)%vr%vf(3)%age
      if (abs(gamSum) > eps) then
        centroidLE = centroidLE/gamSum
        centroidTE = centroidTE/gamSum
        radiusRollup = radiusRollup/gamSum
      else
        centroidLE = this%blade(ib)%waN(this%nNwake, this%rollupEnd)%vr%vf(2)%fc(:, 1)
        centroidTE = this%blade(ib)%waN(this%nNwake, this%rollupEnd)%vr%vf(3)%fc(:, 1)
        radiusRollup = this%blade(ib)%waN(this%nNwake, this%rollupEnd)%vr%vf(3)%rVc
      endif

      ! Suppress Fwake gam if required
      if (this%suppressFwakeSwitch .eq. 1) gamRollup = 0._dp

      ! Initialize far wake tip
      this%blade(ib)%waF(rowFarNext)%vf%fc(:, 2) = centroidLE
      this%blade(ib)%waF(rowFarNext)%vf%fc(:, 1) = centroidTE
      this%blade(ib)%waF(rowFarNext)%gam = gamRollup
      this%blade(ib)%waF(rowFarNext)%vf%age = ageRollup
      this%blade(ib)%waF(rowFarNext)%vf%rVc0 = radiusRollup
      this%blade(ib)%waF(rowFarNext)%vf%rVc = radiusRollup
      ! TRUE => record original length
      call this%blade(ib)%waF(rowFarNext)%vf%calclength(.TRUE.)

      ! Ensure continuity in far wake by assigning
      ! current centroidTE to LE of previous far wake filament
      ! The discontinuity would occur due to convection of
      ! last row of waN in convectwake()
      if (rowFarNext < this%nFwake) then
        this%blade(ib)%waF(rowFarNext + 1)%vf%fc(:, 2) = centroidTE
      endif
    enddo
  end subroutine rotor_rollup

  subroutine rotor_calc_force_gamma(this, density, dt)
    ! Compute force from circulation
  class(rotor_class), intent(inout) :: this
    real(dp), intent(in) :: density, dt
    integer :: ib

    this%forceInertial = 0._dp
    do ib = 1, this%nb
      call this%blade(ib)%calc_force_gamma(density, &
        & sign(1._dp, this%Omega) * sign(1._dp, this%controlPitch(1)) * &
        & sign(1._dp, this%shaftAxis(1)) * &
        & sign(1._dp, this%shaftAxis(2)) * &
        & sign(1._dp, this%shaftAxis(3)), dt)
    enddo
    call this%sumBladeToNetForces()
  end subroutine rotor_calc_force_gamma

  subroutine rotor_calc_force_alpha(this, density, velSound)
    ! Compute force from sec alpha
  class(rotor_class), intent(inout) :: this
    real(dp), intent(in) :: density, velSound
    integer :: ib

    this%forceInertial = 0._dp
    do ib = 1, this%nb
      call this%blade(ib)%calc_force_alpha(density, velSound)
    enddo
    call this%sumBladeToNetForces()
  end subroutine rotor_calc_force_alpha

  subroutine rotor_calc_force_alphaGamma(this, density, velSound, dt)
    ! Compute force from sec alpha
  class(rotor_class), intent(inout) :: this
    real(dp), intent(in) :: density, velSound, dt
    integer :: ib

    this%forceInertial = 0._dp
    do ib = 1, this%nb
      call this%blade(ib)%calc_force_alphaGamma(density, &
        & sign(1._dp, this%Omega*this%controlPitch(1)), &
        & velSound, dt)
    enddo
    call this%sumBladeToNetForces()
  end subroutine rotor_calc_force_alphaGamma

  subroutine rotor_calc_secAlpha(this, updateSecVel)
  class(rotor_class), intent(inout) :: this
    logical, intent(in) :: updateSecVel
    integer :: ib

    do ib = 1, this%nb
      call this%blade(ib)%calc_secAlpha(updateSecVel)
    enddo
  end subroutine rotor_calc_secAlpha

  subroutine rotor_convectwake(this, iter, dt, wakeType)
  class(rotor_class), intent(inout) :: this
    integer, intent(in) :: iter
    real(dp), intent(in) :: dt
    character(len=1), intent(in) :: wakeType  ! For predicted wake
    integer :: ib
    real(dp) :: bladeOffset

    if (this%imposeAxisymmetry == 0) then
      do ib = 1, this%nb
        call this%blade(ib)%convectwake(this%rowNear, this%rowFar, dt, wakeType)
      enddo
    else
      call this%blade(1)%convectwake(this%rowNear, this%rowFar, dt, wakeType)
      do ib = 2, this%nb
        bladeOffset = twoPi/this%nb*(ib - 1)
        ! Copy wakes from blade1
        select case (wakeType)
        case ('C')
          this%blade(ib)%waN(this%rowNear:, :) = &
            & this%blade(1)%waN(this%rowNear:, :)
          this%blade(ib)%waF(this%rowFar:) = &
            & this%blade(1)%waF(this%rowFar:)
        case ('P')
          this%blade(ib)%waNPredicted(this%rowNear:, :) = &
            & this%blade(1)%waNPredicted(this%rowNear:, :)
          this%blade(ib)%waFPredicted(this%rowFar:) = &
            & this%blade(1)%waFPredicted(this%rowFar:)
        case default
          error stop 'ERROR: Wrong character flag for convectwake()'
        end select
        ! Rotate wakes to correct positions
        call this%blade(ib)%rot_wake_axis(bladeOffset, &
          & this%shaftAxis, this%hubCoords, this%rowNear, this%rowFar, wakeType)
      enddo
    endif

    ! Add prescribed wake
    if (this%prescWakeNt > 0 .and. iter > this%prescWakeNt) then
      call this%updatePrescribedWake(dt, wakeType)
    endif

  end subroutine rotor_convectwake

  subroutine rotor_burst_wake(this)
  class(rotor_class), intent(inout) :: this
    integer :: ib
    do ib = 1, this%nb
      call this%blade(ib)%burst_wake(this%rowFar, this%skewLimit, this%chord)
    enddo
  end subroutine rotor_burst_wake

  subroutine rotor_calc_skew(this)
  class(rotor_class), intent(inout) :: this
    integer :: ib
    do ib = 1, this%nb
      call this%blade(ib)%calc_skew(this%rowNear)
    enddo
  end subroutine rotor_calc_skew

  subroutine rotor_dirLiftDrag(this)
  class(rotor_class), intent(inout) :: this
    integer :: ib

    do ib = 1, this%nb
      call this%blade(ib)%dirLiftDrag()
    enddo
  end subroutine rotor_dirLiftDrag

  subroutine rotor_sumBladeToNetForces(this)
  class(rotor_class), intent(inout) :: this
    integer :: ib

    this%forceInertial = 0._dp
    this%lift = 0._dp
    this%drag = 0._dp
    this%liftUnsteady = 0._dp
    this%dragInduced = 0._dp
    this%dragProfile = 0._dp
    this%dragUnsteady = 0._dp

    do ib = 1, this%nb
      this%forceInertial = this%forceInertial + this%blade(ib)%forceInertial
      this%lift = this%lift + this%blade(ib)%lift
      this%drag = this%drag + this%blade(ib)%drag
      this%liftUnsteady = this%liftUnsteady + this%blade(ib)%liftUnsteady
      this%dragInduced = this%dragInduced + this%blade(ib)%dragInduced
      this%dragProfile = this%dragProfile + this%blade(ib)%dragProfile
      this%dragUnsteady = this%dragUnsteady + this%blade(ib)%dragUnsteady
    enddo

  end subroutine rotor_sumBladeToNetForces

  subroutine rotor_mirrorGeometry(this, fromRotor)
    !! Mirrors geometry from another rotor
  class(rotor_class), intent(inout) :: this
  class(rotor_class), intent(in) :: fromRotor
    integer :: ib, ic, is

    do ib = 1, this%nb
      do is = 1, this%ns
        do ic = 1, this%nc
          this%blade(ib)%wiP(ic, is)%nCap = &
            & fromRotor%blade(ib)%wiP(ic, is)%nCap
          this%blade(ib)%wiP(ic, is)%nCap(this%imagePlane) = &
            & -1._dp*this%blade(ib)%wiP(ic, is)%nCap(this%imagePlane)
        enddo
      enddo
    enddo
  end subroutine rotor_mirrorGeometry

  subroutine rotor_mirrorGamma(this, fromRotor)
    !! Mirrors gamma from another rotor
  class(rotor_class), intent(inout) :: this
  class(rotor_class), intent(in) :: fromRotor
    this%gamVec = -1.0_dp * fromRotor%gamVec
  end subroutine rotor_mirrorGamma

  subroutine rotor_mirrorVelCP(this, fromRotor)
    !! Mirrors velCP, velCPm from another rotor
  class(rotor_class), intent(inout) :: this
  class(rotor_class), intent(in) :: fromRotor
    integer :: ic, is, ib

    do ib = 1, this%nb
      !$omp parallel do collapse (2)
      do is = 1, this%ns
        do ic = 1, this%nc
          this%blade(ib)%wiP(ic, is)%velCP = &
            & fromRotor%blade(ib)%wip(ic, is)%velCP
          this%blade(ib)%wiP(ic, is)%velCPm = &
            & fromRotor%blade(ib)%wip(ic, is)%velCPm
        enddo
      enddo
      !$omp end parallel do
    enddo

    do ib = 1, this%nb
      this%blade(ib)%wiP%velCP(this%imagePlane) = &
        & -1._dp * this%blade(ib)%wiP%velCP(this%imagePlane)
      this%blade(ib)%wiP%velCPm(this%imagePlane) = &
        & -1._dp * this%blade(ib)%wiP%velCPm(this%imagePlane)
    enddo
  end subroutine rotor_mirrorVelCP

  subroutine rotor_mirrorWake(this, fromRotor, wakeType)
    !! Mirrors wake positions from another rotor
  class(rotor_class), intent(inout) :: this
  class(rotor_class), intent(in) :: fromRotor
    character(len=1), intent(in) :: wakeType  ! For predicted wake
    integer :: ic, is, ib

    select case (wakeType)
    case ('A')  ! [A]ll current wake
      do ib = 1, this%nb
        !$omp parallel do collapse (2)
        do is = 1, this%ns
          do ic = 1, this%nNwake
            this%blade(ib)%waN(ic, is) = fromRotor%blade(ib)%waN(ic, is)
          enddo
        enddo
        !$omp end parallel do
        !$omp parallel do
        do ic = 1, this%nFwake
          this%blade(ib)%waF(ic) = fromRotor%blade(ib)%waF(ic)
        enddo
        !$omp end parallel do
      enddo

      do ib = 1, this%nb
        !$omp parallel do collapse (2)
        do is = 1, this%ns
          do ic = 1, this%nNwake
            call this%blade(ib)%waN(ic, is)%vr%mirror( &
              & this%imagePlane)
          enddo
        enddo
        !$omp end parallel do
        !$omp parallel do
        do ic = 1, this%nFwake
          call this%blade(ib)%waF(ic)%mirror(this%imagePlane)
        enddo
        !$omp end parallel do
      enddo

    case ('C')  ! [C]urrent wake
      do ib = 1, this%nb
        !$omp parallel do collapse (2)
        do is = 1, this%ns
          do ic = this%rowNear, this%nNwake
            this%blade(ib)%waN(ic, is) = fromRotor%blade(ib)%waN(ic, is)
          enddo
        enddo
        !$omp end parallel do
        !$omp parallel do 
        do ic = this%rowFar, this%nFwake
          this%blade(ib)%waF(ic) = fromRotor%blade(ib)%waF(ic)
        enddo
        !$omp end parallel do
      enddo

      do ib = 1, this%nb
        !$omp parallel do collapse (2)
        do is = 1, this%ns
          do ic = this%rowNear, this%nNwake
            call this%blade(ib)%waN(ic, is)%vr%mirror( &
              & this%imagePlane)
          enddo
        enddo
        !$omp end parallel do
        !$omp parallel do 
        do ic = this%rowFar, this%nFwake
          call this%blade(ib)%waF(ic)%mirror(this%imagePlane)
        enddo
        !$omp end parallel do
      enddo

    case ('P')  ! [P]redicted wake
      do ib = 1, this%nb
        !$omp parallel do collapse (2)
        do is = 1, this%ns
          do ic = this%rowNear, this%nNwake
            this%blade(ib)%waNPredicted(ic, is) = &
              & fromRotor%blade(ib)%waNPredicted(ic, is)
          enddo
        enddo
        !$omp end parallel do
        !$omp parallel do 
        do ic = this%rowFar, this%nFwake
          this%blade(ib)%waFPredicted(ic) = &
            & fromRotor%blade(ib)%waFPredicted(ic)
        enddo
        !$omp end parallel do
      enddo

      do ib = 1, this%nb
        !$omp parallel do collapse (2)
        do is = 1, this%ns
          do ic = this%rowNear, this%nNwake
            call this%blade(ib)%waNPredicted(ic, is)%vr%mirror( &
              & this%imagePlane)
          enddo
        enddo
        !$omp end parallel do
        !$omp parallel do 
        do ic = this%rowFar, this%nFwake
          call this%blade(ib)%waFPredicted(ic)%mirror( &
            & this%imagePlane)
        enddo
        !$omp end parallel do
      enddo
    end select
  end subroutine rotor_mirrorWake

  subroutine rotor_toChordsRevs(this, nsteps, dt)
    !! Converts -ve nsteps to nsteps for corresponding no. of chords or revs
  class(rotor_class), intent(inout) :: this
    integer, intent(inout) :: nsteps
    real(dp), intent(in) :: dt

    if (nsteps < 0) then
      if (abs(this%Omega) < eps) then ! Fixed wing
        ! nt chord distance
        nsteps = ceiling(abs(nsteps)*this%chord/(dt*norm2(this%velBody)))
      else  ! Rotor
        ! nt revs
        nsteps = ceiling(twoPi*abs(nsteps)/(abs(this%Omega)*dt))
      endif
    endif
  end subroutine rotor_toChordsRevs

  subroutine rotor_eraseNwake(this, rowErase)
    !! Erase a near wake row by setting gamma to zero
  class(rotor_class), intent(inout) :: this
    integer, intent(in) :: rowErase
    integer :: ib

    do ib = 1, this%nb
      this%blade(ib)%waN(rowErase, :)%vr%gam = 0._dp
    enddo
  end subroutine rotor_eraseNwake

  subroutine rotor_eraseFwake(this, rowErase)
    !! Erase a far wake row by setting gamma to zero
  class(rotor_class), intent(inout) :: this
    integer, intent(in) :: rowErase
    integer :: ib

    do ib = 1, this%nb
      this%blade(ib)%waF(rowErase)%gam = 0._dp
    enddo
  end subroutine rotor_eraseFwake

  subroutine rotor_updatePrescribedWake(this, dt, wakeType)
    !! Attaches prescribed far wake
  class(rotor_class), intent(inout) :: this
    real(dp), intent(in) :: dt
    character(len=1), intent(in) :: wakeType
    real(dp) :: bladeOffset
    integer :: ib, rowStart

    if (this%prescWakeGenNt == 0) then
      rowStart = this%rowFar
    else
      rowStart = this%nFwakeEnd - this%prescWakeGenNt
    endif

    ! This should ideally be handled by the blade_class
    select case (wakeType)
    case ('C')
      do ib = 1, this%nbConvect
        call this%blade(ib)%wapF%update( &
          & this%blade(ib)%waF(rowStart:this%nFwakeEnd), &
          & this%hubCoords, this%shaftAxis, this%omegaSlow*dt)
      enddo
    case ('P')
      do ib = 1, this%nbConvect
        call this%blade(ib)%wapFPredicted%update( &
          & this%blade(ib)%waFPredicted(rowStart:this%nFwakeEnd), &
          & this%hubCoords, this%shaftAxis, this%omegaSlow*dt)
      enddo
    end select

    if (this%imposeAxisymmetry == 1) then
      do ib = 2, this%nb
        bladeOffset = twoPi/this%nb*(ib - 1)
        select case (wakeType)
        case ('C')
          this%blade(ib)%wapF = this%blade(1)%wapF
          call this%blade(ib)%wapF%rot_wake_axis(bladeOffset, &
            & this%shaftAxis, this%hubCoords)
        case ('P')
          this%blade(ib)%wapFPredicted = this%blade(1)%wapFPredicted
          call this%blade(ib)%wapFPredicted%rot_wake_axis(bladeOffset, &
            & this%shaftAxis, this%hubCoords)
        end select
      enddo
    endif
  end subroutine rotor_updatePrescribedWake

  subroutine rotor_read(this, unit, iostat, iomsg)
  class(rotor_class), intent(inout) :: this
    integer, intent(in) :: unit
    integer, intent(out) :: iostat
    character(len=*), intent(inout) :: iomsg

    read(unit, iostat=iostat, iomsg=iomsg) this%blade, &
      & this%nNwake, this%nFwake, this%omegaSlow, this%shaftAxis, &
      & this%hubCoords, this%cgCoords, &
      & this%forceInertial, this%lift, this%drag, &
      & this%liftUnsteady, &
      & this%dragInduced, this%dragProfile, this%dragUnsteady, &
      & this%liftUnitVec, this%dragUnitVec, this%sideUnitVec, &
      & this%psi, this%AIC, this%AIC_inv, &
      & this%gamVec, this%gamVecPrev, this%RHS, &
      & this%rowNear, this%rowFar
  end subroutine rotor_read

  subroutine rotor_write(this, unit, iostat, iomsg)
  class(rotor_class), intent(in) :: this
    integer, intent(in) :: unit
    integer, intent(out) :: iostat
    character(len=*), intent(inout) :: iomsg

    write(unit, iostat=iostat, iomsg=iomsg) this%blade, &
      & this%nNwake, this%nFwake, this%omegaSlow, this%shaftAxis, &
      & this%hubCoords, this%cgCoords, &
      & this%forceInertial, this%lift, this%drag, &
      & this%liftUnsteady, &
      & this%dragInduced, this%dragProfile, this%dragUnsteady, &
      & this%liftUnitVec, this%dragUnitVec, this%sideUnitVec, &
      & this%psi, this%AIC, this%AIC_inv, &
      & this%gamVec, this%gamVecPrev, this%RHS, &
      & this%rowNear, this%rowFar
  end subroutine rotor_write
end module classdef
