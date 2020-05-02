!------+-------------------+------|
! ++++ | MODULE DEFINITION | ++++ |
!------+-------------------+------|
module vf_classdef
  use libMath
  implicit none

  type vf_class
    real(dp), dimension(3, 2) :: fc  ! filament coords (xyz,1:2)
    real(dp) :: l0 = 0._dp    ! original length
    real(dp) :: lc = 0._dp    ! current length
    real(dp) :: rVc0 = 0._dp ! initial vortex core radius
    real(dp) :: rVc = 0._dp  ! current vortex core radius
    real(dp) :: age = 0._dp   ! vortex age (in s)
  contains
    procedure :: vind => vfclass_vind
    procedure :: calclength => vfclass_calclength
    procedure :: strain => vfclass_strain
  end type vf_class

  real(dp), parameter :: tol = 1.E-6
  real(dp), parameter :: invTol2 = 1.E06
  real(dp), parameter :: inv4pi = 0.25_dp/pi

contains

  ! Efficient implementation to vind calculation
  function vfclass_vind(this, P) result(vind)
    ! Compute induced velocity by unit strength vortex filament
  class(vf_class) :: this
    real(dp), dimension(3) :: vind, P
    real(dp) :: r1Xr2Abs2
    real(dp), dimension(3) :: r1, r2, r0, r1Xr2

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
  end function vfclass_vind

  subroutine vfclass_calclength(this, isOriginal)
    ! Compute length of vortex filament
  class(vf_class) :: this
    logical, intent(in) :: isOriginal
    real(dp), dimension(3) :: delta

    delta = this%fc(:, 1) - this%fc(:, 2)
    this%lc = norm2(delta)
    if (isOriginal .eqv. .TRUE.) this%l0 = norm2(delta)
  end subroutine vfclass_calclength

  subroutine vfclass_strain(this)
    ! Changes core radius according to change in vortex length
  class(vf_class) :: this
    this%rVc = this%rVc0*sqrt(this%l0/this%lc)
  end subroutine vfclass_strain

end module vf_classdef

!------+-------------------+------|
! ++++ | MODULE DEFINITION | ++++ |
!------+-------------------+------|
module vr_classdef

  use vf_classdef
  implicit none
  type vr_class
    type(vf_class), dimension(4) :: vf
    real(dp) :: gam
    real(dp) :: skew
  contains
    procedure :: vind => vrclass_vind
    procedure :: assignP => vrclass_assignP
    procedure :: shiftdP => vrclass_shiftdP
    procedure :: rot => vrclass_rot
    procedure :: calclength => vrclass_calclength
    procedure :: strain => vrclass_strain
    procedure :: burst
    procedure :: getInteriorAngles
    procedure :: getMedianAngle
    procedure :: getMedianCos
  end type vr_class

contains

  function vrclass_vind(this, P) result(vind)
    ! Compute induced velocity by unit strength vortex ring
  class(vr_class) :: this
    real(dp), dimension(3) :: P, vind
    real(dp), dimension(4, 3) :: vindMat
    integer :: i

    vind = 0._dp

    do i = 1, 4
      vindMat(i, :) = this%vf(i)%vind(P)
    enddo
    vind = sum(vindMat, 1)

  end function vrclass_vind

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

  subroutine vrclass_assignP(this, n, P)
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

  end subroutine vrclass_assignP

  subroutine vrclass_shiftdP(this, n, dshift)
    ! Shift coordinates of nth corner by dshift distance
    ! (usually for Udt convection)
  class(vr_class) :: this
    integer, intent(in) :: n
    real(dp), intent(in), dimension(3) :: dshift

    select case (n)
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

  end subroutine vrclass_shiftdP

  subroutine vrclass_rot(this, Tmat)
    ! Rotate vortex ring using Tmat
  class(vr_class) :: this
    real(dp), intent(in), dimension(3, 3) :: Tmat
    integer :: i

    do i = 1, 4
      this%vf(i)%fc(:, 1) = matmul(Tmat, this%vf(i)%fc(:, 1))
      this%vf(i)%fc(:, 2) = matmul(Tmat, this%vf(i)%fc(:, 2))
    enddo

  end subroutine vrclass_rot

  subroutine vrclass_calclength(this, isOriginal)
    ! Calculate length of filaments in vortex ring
  class(vr_class) :: this
    logical, intent(in) :: isOriginal
    integer :: i
    do i = 1, 4
      call this%vf(i)%calclength(isOriginal)
    enddo
  end subroutine vrclass_calclength

  subroutine vrclass_strain(this)
  class(vr_class) :: this
    integer :: i
    do i = 1, 4
      call this%vf(i)%strain()
    enddo
  end subroutine vrclass_strain

  function getInteriorAngles(this)
    ! Obtain interior angles of vortex ring
  class(vr_class) :: this
    real(dp), dimension(4) :: getInteriorAngles
    real(dp), dimension(3) :: p1, p2, p3, p4

    p1 = this%vf(1)%fc(:, 1)
    p2 = this%vf(2)%fc(:, 1)
    p3 = this%vf(3)%fc(:, 1)
    p4 = this%vf(4)%fc(:, 1)

    getInteriorAngles = 0._dp
    getInteriorAngles(1) = getAngleCos(p2 - p1, p4 - p1)
    getInteriorAngles(2) = getAngleCos(p3 - p2, p1 - p2)
    getInteriorAngles(3) = getAngleCos(p4 - p3, p2 - p3)
    getInteriorAngles(4) = getAngleCos(p3 - p4, p1 - p4)
  end function getInteriorAngles

  function getMedianAngle(this)
    ! Obtain median angle of vortex ring
  class(vr_class) :: this
    real(dp) :: getMedianAngle
    real(dp), dimension(3) :: p1, p2, p3, p4

    p1 = this%vf(1)%fc(:, 1)
    p2 = this%vf(2)%fc(:, 1)
    p3 = this%vf(3)%fc(:, 1)
    p4 = this%vf(4)%fc(:, 1)

    getMedianAngle = getAngleCos(p3 + p4 - p1 - p2, p4 + p1 - p2 - p3)
  end function getMedianAngle

  function getMedianCos(this)
  class(vr_class) :: this
    real(dp) :: getMedianCos
    real(dp), dimension(3) :: p1, p2, p3, p4
    real(dp), dimension(3) :: x1Vec, x2Vec

    p1 = this%vf(1)%fc(:, 1)
    p2 = this%vf(2)%fc(:, 1)
    p3 = this%vf(3)%fc(:, 1)
    p4 = this%vf(4)%fc(:, 1)

    x1Vec = p3 + p4 - p1 - p2
    x2Vec = p4 + p1 - p2 - p3
    getMedianCos = abs(dot_product(x1Vec, x2Vec)/ &
      sqrt(dot_product(x1Vec, x1Vec)*dot_product(x2Vec, x2Vec)))
  end function getMedianCos

  subroutine burst(this, skewLimit)
    ! Burst vortex filaments if skewLimit is hit
  class(vr_class) :: this
    real(dp), intent(in) :: skewLimit
    real(dp) :: skewVal

    if ((abs(this%gam) > eps) .and. (skewLimit > eps)) then
      ! skew:  0-good, 1-bad
      skewVal = this%getMedianCos()
      if (skewVal .ge. skewLimit) this%gam = 0._dp
    endif
    this%skew = skewVal

  end subroutine burst

end module vr_classdef

!------+-------------------+------|
! ++++ | MODULE DEFINITION | ++++ |
!------+-------------------+------|
module wingpanel_classdef
  use vr_classdef
  implicit none
  type wingpanel_class
    type(vr_class) :: vr
    real(dp) :: gamPrev
    real(dp) :: gamTrapz
    real(dp), dimension(3, 4) :: pc    ! panel coords
    real(dp), dimension(3) :: cp      ! coll point coords
    real(dp), dimension(3) :: nCap    ! unit normal vector
    real(dp), dimension(3) :: tauCapChord ! unit tangential vector along chord
    real(dp), dimension(3) :: tauCapSpan  ! unit tangential vector along span
    real(dp), dimension(3) :: velCP    ! local vel at CP excluding wing vortices
    real(dp), dimension(3) :: velCPTotal  ! local vel at CP including wing vortices
    real(dp), dimension(3) :: velCPm  ! rel. inertial vel at CP (due to motion)
    real(dp), dimension(3) :: normalForce  ! normalForce vector (inertial frame)
    real(dp), dimension(3) :: chordwiseResVel
    real(dp) :: velPitch             ! pitch velocity
    !real(dp) :: dLift, dDrag          ! magnitudes of panel lift and drag
    real(dp) :: delP                  ! Pressure difference at panel
    real(dp) :: delDiConstant         ! Induced drag (constant part) at panel
    real(dp) :: delDiUnsteady         ! Induced drag (unsteady part) at panel
    real(dp) :: meanChord, meanSpan ! Panel mean dimensions
    real(dp) :: panelArea            ! Panel area for computing lift
    real(dp) :: rHinge  ! dist to point about which pitching occurs (LE of wing)
    real(dp) :: alpha                 ! local angle of attack
  contains
    procedure :: assignP => wingpanel_class_assignP
    procedure :: calcCP => wingpanel_class_calcCP
    procedure :: calcN => wingpanel_class_calcN
    procedure :: calcTau => wingpanel_class_calcTau
    procedure :: rot => wingpanel_class_rot
    procedure :: shiftdP => wingpanel_class_shiftdP
    procedure :: calc_chordwiseResVel => wingpanel_calc_chordwiseResVel
    procedure :: calc_area
    procedure :: calc_mean_dimensions
    procedure :: isCPinsidecore
  end type wingpanel_class

contains

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

  subroutine wingpanel_class_assignP(this, n, P)
    ! Assign coordinates to nth corner
  class(wingpanel_class) :: this
    integer, intent(in) :: n
    real(dp), dimension(3) :: P

    if (n > 0 .and. n < 5) then
      this%pc(:, n) = P
    else
      error stop 'n may only take values 1,2,3 or 4'
    endif

  end subroutine wingpanel_class_assignP

  subroutine wingpanel_class_calcCP(this)
    ! Compute collocation point location
  class(wingpanel_class) :: this
    this%CP = ((this%PC(:, 1) + this%PC(:, 4))*0.25_dp &
      + (this%PC(:, 2) + this%PC(:, 3))*0.75_dp)*0.5_dp
  end subroutine wingpanel_class_calcCP

  subroutine wingpanel_class_calcN(this)
    ! Compute normal vector
  class(wingpanel_class) :: this
    this%nCap = unitVec(cross_product(this%pc(:, 3) - this%pc(:, 1), &
      this%pc(:, 4) - this%pc(:, 2)))
  end subroutine wingpanel_class_calcN

  subroutine wingpanel_class_calcTau(this)
    ! Compute chordwise and spanwise tangential vectors
  class(wingpanel_class) :: this
    this%tauCapChord = unitVec(0.5_dp*((this%pc(:, 2) + this%pc(:, 3)) &
      - (this%pc(:, 1) + this%pc(:, 4))))
    this%tauCapSpan = unitVec(0.5_dp*((this%pc(:, 3) + this%pc(:, 4)) &
      - (this%pc(:, 2) + this%pc(:, 1))))
  end subroutine wingpanel_class_calcTau

  subroutine wingpanel_class_rot(this, Tmat)
    ! Rotate panel using transformation matrix
  class(wingpanel_class) :: this
    real(dp), dimension(3, 3) :: Tmat
    integer :: i

    do i = 1, 4
      this%pc(:, i) = matmul(Tmat, this%pc(:, i))
    enddo
    call this%vr%rot(Tmat)
    this%CP = matmul(Tmat, this%CP)
    this%nCap = matmul(Tmat, this%nCap)
    this%tauCapChord = matmul(Tmat, this%tauCapChord)
    this%tauCapSpan = matmul(Tmat, this%tauCapSpan)
  end subroutine wingpanel_class_rot

  subroutine wingpanel_class_shiftdP(this, dshift)
    ! Shift corners of vortex ring by dshift
  class(wingpanel_class) :: this
    real(dp), intent(in), dimension(3) :: dshift
    integer :: i

    this%CP = this%CP + dshift
    do i = 1, 4
      this%pc(:, i) = this%pc(:, i) + dshift
      call this%vr%shiftdP(i, dshift)
    enddo

  end subroutine wingpanel_class_shiftdP

  subroutine calc_area(this)
  class(wingpanel_class) :: this
    this%panelArea = 0.5_dp*norm2(cross_product(this%pc(:, 3) &
      - this%pc(:, 1), this%pc(:, 4) - this%pc(:, 2)))
  end subroutine calc_area

  subroutine calc_mean_dimensions(this)
    ! Calculate mean chord and mean span
  class(wingpanel_class) :: this
    this%meanSpan = 0.5_dp*(norm2(this%pc(:, 4) - this%pc(:, 1)) &
      + norm2(this%pc(:, 3) - this%pc(:, 2)))
    this%meanChord = 0.5_dp*(norm2(this%pc(:, 2) - this%pc(:, 1)) &
      + norm2(this%pc(:, 3) - this%pc(:, 4)))
  end subroutine calc_mean_dimensions

  function isCPinsidecore(this)
    ! Check whether collocation point lies
    ! inside viscous core region of vortex ring
  class(wingpanel_class), intent(in) :: this
    logical :: isCPinsidecore
    real(dp) :: deltaxby4, deltayby2

    deltaxby4 = 0.25_dp*abs(this%vr%vf(1)%fc(1, 1) - this%vr%vf(2)%fc(1, 1))
    deltayby2 = 0.5_dp*abs(this%vr%vf(1)%fc(2, 1) - this%vr%vf(4)%fc(2, 1))

    isCPinsidecore = .false.
    if (deltayby2 .lt. this%vr%vf(1)%rVc) then
      isCPinsidecore = .true.  ! Left edge
    elseif (deltayby2 .lt. this%vr%vf(3)%rVc) then
      isCPinsidecore = .true.  ! Right edge
    elseif (deltaxby4 .lt. this%vr%vf(2)%rVc) then
      isCPinsidecore = .true.  ! Upper edge
    elseif (3._dp*deltaxby4 .lt. this%vr%vf(4)%rVc) then
      isCPinsidecore = .true.  ! Bottom edge
    endif
  end function isCPinsidecore

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
  !    this%alpha=acos(dot_product(velCPTotalChordwiseProjected,this%tauCapChord) &
  !      /velCPTotalChordwiseProjectedMagnitude)
  !  else
  !    this%alpha=0._dp
  !  endif
  !end subroutine wingpanel_calc_alpha

end module wingpanel_classdef

!------+-------------------+------|
! ++++ | MODULE DEFINITION | ++++ |
!------+-------------------+------|
module Nwake_classdef
  use vr_classdef
  implicit none
  type Nwake_class
    type(vr_class) :: vr
  end type Nwake_class

contains

  ! VR coordinates
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

end module Nwake_classdef

!------+-------------------+------|
! ++++ | MODULE DEFINITION | ++++ |
!------+-------------------+------|
module Fwake_classdef
  use vf_classdef
  implicit none
  type Fwake_class
    type(vf_class) :: vf
    real(dp) :: gam

  contains
    procedure :: shiftdP => Fwake_shiftdP
    procedure :: assignP => Fwake_assignP
  end type Fwake_class

contains

  ! VF coordinates
  ! o---------> Y along span
  ! |
  ! |   1
  ! |   |
  ! |   |
  ! |   |1
  ! |   |
  ! |   |
  ! |   2
  ! |
  ! V X along chord

  subroutine Fwake_shiftdP(this, n, dshift)
    ! Shift coordinates of nth corner by dshift distance
    ! (usually for Udt convection)
  class(Fwake_class) :: this
    integer, intent(in) :: n
    real(dp), intent(in), dimension(3) :: dshift

    if (n /= 1 .and. n /= 2) &
      error stop 'n may only take values 1 or 2 in Fwake_shiftdP()'
    this%vf%fc(:, n) = this%vf%fc(:, n) + dshift
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
end module Fwake_classdef

!------+-------------------+------|
! ++++ | MODULE DEFINITION | ++++ |
!------+-------------------+------|
module blade_classdef
  use wingpanel_classdef
  use Nwake_classdef
  use Fwake_classdef
  use libC81
  implicit none
  type blade_class
    type(wingpanel_class), allocatable, dimension(:, :) :: wiP
    type(Nwake_class), allocatable, dimension(:, :) :: waP
    type(Fwake_class), allocatable, dimension(:) :: waF
    type(Nwake_class), allocatable, dimension(:, :) :: waPPredicted
    type(Fwake_class), allocatable, dimension(:) :: waFPredicted
    type(C81_class), allocatable, dimension(:) :: C81
    real(dp) :: theta, psi, pivotLE
    real(dp), dimension(3) :: forceInertial
    real(dp), dimension(3) :: lift, drag, dragInduced, dragProfile, dragUnsteady
    integer, allocatable, dimension(:) :: airfoilNo
    character(len=30), allocatable, dimension(:) :: airfoilFile
    real(dp), allocatable, dimension(:) :: airfoilSectionLimit
    real(dp), allocatable, dimension(:, :, :) :: velNwake, velNwake1, velNwake2, velNwake3
    real(dp), allocatable, dimension(:, :, :) :: velNwakePredicted, velNwakeStep
    real(dp), allocatable, dimension(:, :) :: velFwake, velFwake1, velFwake2, velFwake3
    real(dp), allocatable, dimension(:, :) :: velFwakePredicted, velFwakeStep
    real(dp), dimension(3) :: xAxis, yAxis, zAxis
    ! Sectional quantities
    real(dp), allocatable, dimension(:, :) :: secForceInertial, secLift, secDrag
    real(dp), allocatable, dimension(:, :) :: secLiftDir, secDragDir
    real(dp), allocatable, dimension(:, :) :: secDragInduced, secDragProfile, secDragUnsteady
    real(dp), allocatable, dimension(:, :) :: secTauCapChord, secTauCapSpan
    real(dp), allocatable, dimension(:, :) :: secNormalVec, secVelFreestream
    real(dp), allocatable, dimension(:, :) :: secChordwiseResVel, secCP
    real(dp), allocatable, dimension(:) :: secAlpha, secCL, secCD, secCM, CL0, CLa

  contains
    procedure :: move => blade_move
    procedure :: rot_pitch
    procedure :: rot_axis
    procedure :: rot_pts => blade_rot_pts
    procedure :: vind_bywing => blade_vind_bywing
    procedure :: vind_bywing_boundVortices => blade_vind_bywing_boundVortices
    procedure :: vind_bywing_chordwiseVortices => blade_vind_bywing_chordwiseVortices
    procedure :: vind_boundVortex => blade_vind_boundVortex
    procedure :: vind_bywake => blade_vind_bywake
    procedure :: convectwake
    procedure :: wake_continuity
    procedure :: getSecDynamicPressure
    procedure :: getSecArea
    procedure :: calc_force_gamma => blade_calc_force_gamma
    procedure :: calc_force_alpha => blade_calc_force_alpha
    procedure :: calc_force_alphaGamma => blade_calc_force_alphaGamma
    procedure :: calc_secAlpha => blade_calc_secAlpha
    procedure :: calc_secChordwiseResVel => blade_calc_secChordwiseResVel
    procedure :: burst_wake => blade_burst_wake
    procedure :: getSecChordwiseLocations
    procedure :: calc_secCoeffs
    procedure :: secCoeffsToSecForces
    procedure :: dirLiftDrag => blade_dirLiftDrag
    procedure :: sumSecToNetForces
  end type blade_class
contains

  subroutine blade_move(this, dshift)
    ! Move blade by dshift
  class(blade_class) :: this
    real(dp), intent(in), dimension(3) :: dshift
    integer :: i, j

    do j = 1, size(this%wiP, 2)
      do i = 1, size(this%wiP, 1)
        call this%wiP(i, j)%shiftdP(dshift)
      enddo
    enddo

    do i = 1, size(this%secCP, 2)
      this%secCP(:, i) = this%secCP(:, i) + dshift
    enddo

  end subroutine blade_move

  subroutine blade_rot_pts(this, pts, origin, order)
    ! Rotate blade using pts => phi theta psi
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

    do j = 1, size(this%wiP, 2)
      do i = 1, size(this%wiP, 1)
        call this%wiP(i, j)%shiftdP(-origin)
        call this%wiP(i, j)%rot(TMat)
        call this%wiP(i, j)%shiftdP(origin)
      enddo
    enddo

    this%xAxis = matmul(Tmat, this%xAxis)
    this%yAxis = matmul(Tmat, this%yAxis)
    this%zAxis = matmul(Tmat, this%zAxis)

    do i = 1, size(this%secCP, 2)
      this%secCP(:, i) = this%secCP(:, i) - origin
      this%secCP(:, i) = matmul(TMat, this%secCP(:, i))
      this%secCP(:, i) = this%secCP(:, i) + origin
    enddo

    ! Rotate sec vectors
    do j = 1, size(this%wiP, 2)
      this%secTauCapChord(:, j) = matmul(TMat, this%secTauCapChord(:, j))
      this%secNormalVec(:, j) = matmul(TMat, this%secNormalVec(:, j))
    enddo

  end subroutine blade_rot_pts

  subroutine rot_pitch(this, theta)
    ! Rotate blade by pitch angle
    ! pivot point calculated using straight line joining
    ! LE of first panel and TE of last panel
  class(blade_class), intent(inout) :: this
    real(dp), intent(in) :: theta
    real(dp), dimension(3) :: axis
    real(dp), dimension(3) :: axisOrigin!, axisEnd
    integer :: rows, cols

    if (abs(theta) > eps) then
      rows = size(this%wiP, 1)
      cols = size(this%wiP, 2)
      axisOrigin = this%wiP(1, 1)%PC(:, 1)*(1._dp - this%pivotLE) &
        + this%wiP(rows, 1)%PC(:, 2)*this%pivotLE
      !axisEnd=this%wiP(1,cols)%PC(:,4)*(1._dp-this%pivotLE) &
      !  +this%wiP(rows,cols)%PC(:,3)*this%pivotLE
      !
      !! Construct axes of rotation from LE of first panel
      !axis=axisEnd-axisOrigin
      !axis=axis/norm2(axis)

      ! Use blade X axis for rotation
      axis = this%yAxis

      call this%rot_axis(theta, axis, axisOrigin)
    endif
  end subroutine rot_pitch

  subroutine rot_axis(this, theta, axisVec, origin)
    ! Rotate about axis at specified origin
  class(blade_class), intent(inout) :: this
    real(dp), intent(in) :: theta
    real(dp), intent(in), dimension(3) :: axisVec
    real(dp), intent(in), dimension(3) :: origin
    real(dp), dimension(3, 3) :: Tmat
    real(dp), dimension(3) :: axis
    integer :: i, j, rows
    real(dp) :: ct, st, omct

    if (abs(theta) > eps) then
      ! Translate to origin
      rows = size(this%wiP, 1)
      call this%move(-origin)

      ! Ensure axis is normalized
      axis = axisVec/norm2(axisVec)

      ! Calculate TMat
      ct = cos(theta)
      st = sin(theta)
      omct = 1 - ct

      Tmat(:, 1) = (/ct + axis(1)*axis(1)*omct, &
        axis(3)*st + axis(2)*axis(1)*omct, &
        -axis(2)*st + axis(3)*axis(1)*omct/)
      Tmat(:, 2) = (/-axis(3)*st + axis(1)*axis(2)*omct, &
        ct + axis(2)*axis(2)*omct, &
        axis(1)*st + axis(3)*axis(2)*omct/)
      Tmat(:, 3) = (/axis(2)*st + axis(1)*axis(3)*omct, &
        -axis(1)*st + axis(2)*axis(3)*omct, &
        ct + axis(3)*axis(3)*omct/)

      ! Rotate about axis
      do j = 1, size(this%wiP, 2)
        do i = 1, size(this%wiP, 1)
          call this%wiP(i, j)%rot(TMat)
        enddo
      enddo

      ! Untranslate from origin
      call this%move(origin)

      ! Rotate secCP also
      do i = 1, size(this%secCP, 2)
        this%secCP(:, i) = this%secCP(:, i) - origin
        this%secCP(:, i) = matmul(TMat, this%secCP(:, i))
        this%secCP(:, i) = this%secCP(:, i) + origin
      enddo

      ! Rotate sec vectors also along with blade
      do j = 1, size(this%wiP, 2)
        this%secTauCapChord(:, j) = matmul(TMat, this%secTauCapChord(:, j))
        this%secNormalVec(:, j) = matmul(TMat, this%secNormalVec(:, j))
      enddo

      ! Rotate blade axes
      this%xAxis = matmul(TMat, this%xAxis)
      this%yAxis = matmul(TMat, this%yAxis)
      this%zAxis = matmul(TMat, this%zAxis)
    endif
  end subroutine rot_axis

  function blade_vind_bywing(this, P)
    ! Compute induced velocity by blade bound vorticity
  class(blade_class), intent(inout) :: this
    real(dp), intent(in), dimension(3) :: P
    real(dp), dimension(3) :: blade_vind_bywing
    integer :: i, j

    blade_vind_bywing = 0._dp
    do j = 1, size(this%wiP, 2)
      do i = 1, size(this%wiP, 1)
        blade_vind_bywing = blade_vind_bywing + &
          this%wiP(i, j)%vr%vind(P)*this%wiP(i, j)%vr%gam
      enddo
    enddo

  end function blade_vind_bywing

  function blade_vind_bywing_boundVortices(this, P)
    ! Compute induced velocity by bound vortices alone
  class(blade_class), intent(inout) :: this
    real(dp), intent(in), dimension(3) :: P
    real(dp), dimension(3) :: blade_vind_bywing_boundVortices
    integer :: i, j, rows, cols

    rows = size(this%wiP, 1)
    cols = size(this%wiP, 2)

    blade_vind_bywing_boundVortices = 0._dp
    do j = 1, cols
      do i = 1, rows
        blade_vind_bywing_boundVortices = blade_vind_bywing_boundVortices + &
          (this%wiP(i, j)%vr%vf(2)%vind(P) + this%wiP(i, j)%vr%vf(4)%vind(P))* &
          this%wiP(i, j)%vr%gam
      enddo
    enddo
    do j = 1, cols
      blade_vind_bywing_boundVortices = blade_vind_bywing_boundVortices - &
        this%wiP(rows, j)%vr%vf(2)%vind(P)*this%wiP(rows, j)%vr%gam
    enddo
  end function blade_vind_bywing_boundVortices

  function blade_vind_bywing_chordwiseVortices(this, P)
    ! Compute induced velocity by bound vortices alone
  class(blade_class), intent(inout) :: this
    real(dp), intent(in), dimension(3) :: P
    real(dp), dimension(3) :: blade_vind_bywing_chordwiseVortices
    integer :: i, j, rows, cols

    rows = size(this%wiP, 1)
    cols = size(this%wiP, 2)

    blade_vind_bywing_chordwiseVortices = 0._dp
    do j = 1, cols
      do i = 1, rows
        blade_vind_bywing_chordwiseVortices = blade_vind_bywing_chordwiseVortices + &
          (this%wiP(i, j)%vr%vf(1)%vind(P) + this%wiP(i, j)%vr%vf(3)%vind(P))* &
          this%wiP(i, j)%vr%gam
      enddo
    enddo
    do j = 1, cols
      blade_vind_bywing_chordwiseVortices = blade_vind_bywing_chordwiseVortices + &
        this%wiP(rows, j)%vr%vf(2)%vind(P)*this%wiP(rows, j)%vr%gam
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
  class(blade_class), intent(inout) :: this
    integer, intent(in) :: rowNear, rowFar
    real(dp), intent(in), dimension(3) :: P
    character(len=1), optional :: optionalChar
    real(dp), dimension(3) :: blade_vind_bywake
    integer :: i, j, nNwake, nFwake

    nNwake = size(this%waP, 1)
    nFwake = size(this%waF, 1)
    blade_vind_bywake = 0._dp
    if (.not. present(optionalChar)) then
      do j = 1, size(this%waP, 2)
        do i = rowNear, nNwake
          if (abs(this%waP(i, j)%vr%gam) .gt. eps) &
            blade_vind_bywake = blade_vind_bywake + &
            this%waP(i, j)%vr%vind(P)*this%waP(i, j)%vr%gam
        enddo
      enddo

      if (rowFar .le. nFwake) then
        ! Last row of Nwake is made of horseshoe vortices, if Fwake is generated
        do j = 1, size(this%waP, 2)
          blade_vind_bywake = blade_vind_bywake - &
            this%waP(nNwake, j)%vr%vf(2)%vind(P)*this%waP(nNwake, j)%vr%gam
        enddo

        do i = rowFar, nFwake
          if (abs(this%waF(i)%gam) .gt. eps) &
            blade_vind_bywake = blade_vind_bywake &
            + this%waF(i)%vf%vind(P)*this%waF(i)%gam
        enddo
      endif
    elseif ((optionalChar .eq. 'P') .or. (optionalChar .eq. 'p')) then
      do j = 1, size(this%waP, 2)
        do i = rowNear, nNwake
          if (abs(this%waPPredicted(i, j)%vr%gam) .gt. eps) &
            blade_vind_bywake = blade_vind_bywake &
            + this%waPPredicted(i, j)%vr%vind(P)*this%waPPredicted(i, j)%vr%gam
        enddo
      enddo

      if (rowFar .le. nFwake) then
        ! Last row of Nwake is made of horseshoe vortices, if Fwake is generated
        do j = 1, size(this%waP, 2)
          blade_vind_bywake = blade_vind_bywake &
            - this%waPPredicted(nNwake, j)%vr%vf(2)%vind(P)*this%waPPredicted(nNwake, j)%vr%gam
        enddo

        do i = rowFar, nFwake
          if (abs(this%waFPredicted(i)%gam) .gt. eps) &
            blade_vind_bywake = blade_vind_bywake &
            + this%waFPredicted(i)%vf%vind(P)*this%waFPredicted(i)%gam
        enddo
      endif
    else
      error stop 'ERROR: Wrong character flag for blade_vind_bywake()'
    endif

  end function blade_vind_bywake

  subroutine convectwake(this, rowNear, rowFar, dt, wakeType)
    ! Convect wake collocation points using velNwake matrix
  class(blade_class), intent(inout) :: this
    integer, intent(in) :: rowNear, rowFar
    real(dp), intent(in) :: dt
    character(len=1), intent(in) :: wakeType  ! For predicted wake
    integer :: i, j, cols, nNwake, nFwake

    cols = size(this%waP, 2)
    nNwake = size(this%waP, 1)

    select case (wakeType)
    case ('C')    ! [C]urrent wake
      !$omp parallel do collapse(2)
      do j = 1, cols
        do i = rowNear, nNwake
          call this%waP(i, j)%vr%shiftdP(2, this%velNwake(:, i, j)*dt)
        enddo
      enddo
      !$omp end parallel do

      !$omp parallel do
      do i = rowNear, nNwake
        call this%waP(i, cols)%vr%shiftdP(3, this%velNwake(:, i, cols + 1)*dt)
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
      do j = 1, cols
        do i = 1, rowNear, nNwake
          call this%waPPredicted(i, j)%vr%shiftdP(2, this%velNwake(:, i, j)*dt)
        enddo
      enddo
      !$omp end parallel do

      !$omp parallel do
      do i = 1, rowNear, nNwake
        call this%waPPredicted(i, cols)%vr%shiftdP(3, this%velNwake(:, i, cols + 1)*dt)
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

  end subroutine convectwake

  subroutine wake_continuity(this, rowNear, rowFar, wakeType)
    ! Maintain continuity between vortex ring elements after convection
    ! of wake collocation points
  class(blade_class), intent(inout) :: this
    integer, intent(in) :: rowNear, rowFar
    character(len=1), intent(in) :: wakeType  ! For predicted wake
    integer :: i, j, nNwake, nFwake, cols

    nNwake = size(this%waP, 1)
    cols = size(this%waP, 2)

    select case (wakeType)
    case ('C')
      !$omp parallel do collapse(2)
      do j = 1, cols - 1
        do i = rowNear + 1, nNwake
          call this%waP(i, j)%vr%assignP(1, this%waP(i - 1, j)%vr%vf(2)%fc(:, 1))
          call this%waP(i, j)%vr%assignP(3, this%waP(i, j + 1)%vr%vf(2)%fc(:, 1))
          call this%waP(i, j)%vr%assignP(4, this%waP(i - 1, j + 1)%vr%vf(2)%fc(:, 1))
        enddo
      enddo
      !$omp end parallel do

      !$omp parallel do
      do j = 1, cols - 1
        call this%waP(rowNear, j)%vr%assignP(3, this%waP(rowNear, j + 1)%vr%vf(2)%fc(:, 1))
      enddo
      !$omp end parallel do

      !$omp parallel do
      do i = rowNear + 1, nNwake
        call this%waP(i, cols)%vr%assignP(1, this%waP(i - 1, cols)%vr%vf(2)%fc(:, 1))
        call this%waP(i, cols)%vr%assignP(4, this%waP(i - 1, cols)%vr%vf(3)%fc(:, 1))
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
      do j = 1, cols - 1
        do i = rowNear + 1, nNwake
          call this%waPPredicted(i, j)%vr%assignP(1, this%waPPredicted(i - 1, j)%vr%vf(2)%fc(:, 1))
          call this%waPPredicted(i, j)%vr%assignP(3, this%waPPredicted(i, j + 1)%vr%vf(2)%fc(:, 1))
          call this%waPPredicted(i, j)%vr%assignP(4, this%waPPredicted(i - 1, j + 1)%vr%vf(2)%fc(:, 1))
        enddo
      enddo
      !$omp end parallel do

      !$omp parallel do
      do j = 1, cols - 1
        call this%waPPredicted(rowNear, j)%vr%assignP(3, this%waPPredicted(rowNear, j + 1)%vr%vf(2)%fc(:, 1))
      enddo
      !$omp end parallel do

      !$omp parallel do
      do i = rowNear + 1, nNwake
        call this%waPPredicted(i, cols)%vr%assignP(1, this%waPPredicted(i - 1, cols)%vr%vf(2)%fc(:, 1))
        call this%waPPredicted(i, cols)%vr%assignP(4, this%waPPredicted(i - 1, cols)%vr%vf(3)%fc(:, 1))
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

  end subroutine wake_continuity

  subroutine blade_calc_force_gamma(this, density, invertGammaSign, dt)
    ! Compute force using blade circulation
  class(blade_class), intent(inout) :: this
    real(dp), intent(in) :: density, invertGammaSign, dt
    integer :: is, ic, rows, cols
    real(dp) :: unsteadyTerm
    real(dp), dimension(size(this%wiP, 1), size(this%wiP, 2)) :: velTangentialChord, velTangentialSpan
    real(dp), dimension(size(this%wiP, 1), size(this%wiP, 2)) :: velInduced
    real(dp), dimension(size(this%wiP, 1), size(this%wiP, 2)) :: gamElementChord, gamElementSpan
    rows = size(this%wiP, 1)
    cols = size(this%wiP, 2)

    this%forceInertial = 0._dp
    this%secForceInertial = 0._dp
    this%secLift = 0._dp
    this%secDrag = 0._dp

    ! Compute tangential velocity
    do is = 1, cols
      do ic = 1, rows
        velTangentialChord(ic, is) = dot_product(this%wiP(ic, is)%velCP, this%wiP(ic, is)%tauCapChord)
        velTangentialSpan(ic, is) = dot_product(this%wiP(ic, is)%velCP, this%wiP(ic, is)%tauCapSpan)
        velInduced(ic, is) = dot_product(this%wiP(ic, is)%velCP + &
          this%vind_bywing_chordwiseVortices(this%wiP(ic, is)%CP), &
          unitVec(cross_product(this%wiP(ic, is)%velCPm,this%yAxis)))
      enddo
    enddo

    ! Compute chordwise elemental circulation of edge panels
    do is = 1, cols
      gamElementChord(1, is) = this%wiP(1, is)%vr%gam
    enddo
    do ic = 2, rows
      gamElementChord(ic, 1) = this%wiP(ic, 1)%vr%gam - this%wiP(ic - 1, 1)%vr%gam
    enddo

    ! Compute spanwise elemental circulation of edge panels
    do ic = 1, rows
      gamElementSpan(ic, 1) = this%wiP(ic, 1)%vr%gam
    enddo
    do is = 2, cols
      gamElementSpan(1, is) = this%wiP(1, is)%vr%gam - this%wiP(1, is - 1)%vr%gam
    enddo

    ! Compute chordwise and spanwise elemental circulations of inner panels
    do is = 2, cols
      do ic = 2, rows
        gamElementChord(ic, is) = this%wiP(ic, is)%vr%gam - this%wiP(ic - 1, is)%vr%gam
        gamElementSpan(ic, is) = this%wiP(ic, is)%vr%gam - this%wiP(ic, is - 1)%vr%gam
      enddo
    enddo

    ! Invert gamma sign for correct computation
    gamElementSpan = -1._dp*gamElementSpan
    gamElementChord = -1._dp*gamElementChord

    ! Compute delP
    do is = 1, cols
      do ic = 1, rows
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
        unsteadyTerm = (this%wiP(ic, is)%gamTrapz - this%wiP(ic, is)%gamPrev)/dt
        this%wiP(ic, is)%delP = density * &
          (velTangentialChord(ic, is) * gamElementChord(ic, is) / this%wiP(ic, is)%meanChord + &
          velTangentialSpan(ic, is) * gamElementSpan(ic, is) / this%wiP(ic, is)%meanSpan + &
          unsteadyTerm)
        ! -1.0 multiplied to invert sign of gamma
        this%wiP(ic, is)%gamPrev = this%wiP(ic, is)%gamTrapz

        ! Compute induced drag
        ! velInduced in either direction doesnt change drag direction
        this%wiP(ic, is)%delDiConstant = density*abs(velInduced(ic, is))* &
          gamElementChord(ic, is)*this%wiP(ic, is)%meanSpan
        this%wiP(ic, is)%delDiUnsteady = density*unsteadyTerm*this%wiP(ic, is)%panelArea * &
          dot_product(this%wiP(ic, is)%nCap, unitVec(this%wiP(ic, is)%velCPm))

        ! Invert direction of forceInertial according to sign of omega and collective pitch
        this%wiP(ic, is)%normalForce = this%wiP(ic, is)%delP* &
          this%wiP(ic, is)%panelArea*this%wiP(ic, is)%nCap*(-1._dp)*invertGammaSign

        this%secForceInertial(:, is) = this%secForceInertial(:, is) + this%wiP(ic, is)%normalForce
        this%forceInertial = this%forceInertial + this%wiP(ic, is)%normalForce

        this%secLift(:, is) = this%secLift(:, is) + projVec(this%wiP(ic, is)%normalForce, &
          this%secLiftDir(:, is))
      enddo
      this%secDragInduced(:, is) = this%secDragDir(:, is)* &
        sum(this%wiP(:, is)%delDiConstant + this%wiP(:, is)%delDiUnsteady)
      ! Drag unsteady is purely for monitoring purposes if required
      ! and is not used for computations anywhere
      this%secDragUnsteady(:, is) = this%secDragDir(:, is)* &
        sum(this%wiP(:, is)%delDiUnsteady)
    enddo
    this%secDragProfile = 0._dp  ! To overwrite unit vectors previously assigned in main.f90

    this%secDrag = this%secDragInduced + this%secDragProfile

    call this%sumSecToNetForces()

  end subroutine blade_calc_force_gamma

  function getSecDynamicPressure(this, density)
  class(blade_class), intent(inout) :: this
    real(dp), intent(in) :: density
    real(dp), dimension(size(this%wiP, 2)) :: magsecVelCPTotal
    real(dp), dimension(size(this%wiP, 2)) :: getSecDynamicPressure
    integer :: is, rows

    rows = size(this%wiP, 1)
    do is = 1, size(this%wiP, 2)
      magsecVelCPTotal(is) = norm2(this%secChordwiseResVel(:, is))
    enddo
    getSecDynamicPressure = 0.5_dp*density*magsecVelCPTotal**2._dp
  end function getSecDynamicPressure

  function getSecArea(this)
  class(blade_class), intent(inout) :: this
    real(dp), dimension(size(this%wiP, 2)) :: getSecArea
    integer :: is

    do is = 1, size(this%wiP, 2)
      getSecArea(is) = sum(this%wiP(:, is)%panelArea)
    enddo
  end function getSecArea

  subroutine blade_calc_force_alpha(this, density, velSound)
    ! Compute force using sec alpha
  class(blade_class), intent(inout) :: this
    real(dp), intent(in) :: density, velSound

    call this%calc_secCoeffs(velSound)

    call this%secCoeffsTosecForces(density)
    this%secDragInduced = 0._dp  ! To overwrite unit vectors previously assigned in main.f90
    this%secDrag = this%secDragProfile + this%secDragInduced

    call this%sumSecToNetForces()

  end subroutine blade_calc_force_alpha

  subroutine blade_calc_force_alphaGamma(this, density, invertGammaSign, velSound, dt)
    ! Compute force using alpha approximated from sec circulation
  class(blade_class), intent(inout) :: this
    real(dp), intent(in) :: density, invertGammaSign, velSound, dt
    real(dp), dimension(3) :: secChordwiseVelFreestream, liftDir
    real(dp), dimension(size(this%wiP, 2)) :: secDynamicPressure, secArea
    integer :: is, ns

    ns = size(this%wiP, 2)

    ! Compute unsteady sec lift from gamma distribution
    call this%calc_force_gamma(density, invertGammaSign, dt)

    do is = 1, ns
      ! Compute sec freestream velocity
      secChordwiseVelFreestream = this%secVelFreestream(:, is) - &
        dot_product(this%secVelFreestream(:, is), this%yAxis)*this%yAxis

      ! Assuming sec resultant velocity is same as sec freestream vel
      ! for computing corrected alpha later
      this%secChordwiseResVel(:, is) = secChordwiseVelFreestream
    enddo

    secDynamicPressure = this%getSecDynamicPressure(density)
    secArea = this%getSecArea()

    do is = 1, ns
      ! Extract sectional lift and CL
      liftDir = cross_product(this%secChordwiseResVel(:, is), this%yAxis)
      ! Assuming gamma method only gives lift
      this%secCL(is) = dot_product(this%secForceInertial(:, is), unitVec(liftDir)) &
        /(secDynamicPressure(is)*secArea(is))

      ! Compute angle of attack from linear CL
      this%secAlpha(is) = (this%secCL(is) - this%CL0(this%airfoilNo(is)))/this%CLa(this%airfoilNo(is))
    enddo

    ! Compute non-linear CL
    call this%calc_secCoeffs(velSound)

    call this%secCoeffsToSecForces(density)

    this%secDrag = this%secDragInduced + this%secDragProfile

    call this%sumSecToNetForces()

  end subroutine blade_calc_force_alphaGamma

  subroutine secCoeffsToSecForces(this, density)
  class(blade_class), intent(inout) :: this
    real(dp), intent(in) :: density
    integer :: is
    real(dp), dimension(size(this%wiP, 2)) :: leadingTerm

    leadingTerm = this%getSecDynamicPressure(density)*this%getSecArea()

    do is = 1, size(this%wiP, 2)
      ! Lift and Drag vectors
      this%secLift(:, is) = this%secLiftDir(:, is)*leadingTerm(is)*this%secCL(is) 
      this%secDragProfile(:, is) = this%secDragDir(:, is)*leadingTerm(is)*this%secCD(is)
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
  end subroutine secCoeffsToSecForces

  subroutine calc_secCoeffs(this, velSound)
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
  end subroutine calc_secCoeffs

  subroutine blade_calc_secChordwiseResVel(this)
    ! Compute sec resultant velocity by interpolating local panel velocity
  class(blade_class), intent(inout) :: this
    integer :: i, is, ic, rows
    real(dp), dimension(size(this%wiP, 1)) :: xDist

    rows = size(this%wiP, 1)
    if (rows .ge. 3) then  ! Use least squares fit to get sec resultant velocity
      do is = 1, size(this%secChordwiseResVel, 2)
        do ic = 1, rows
          call this%wiP(ic, is)%calc_chordwiseResVel()
          xDist(ic) = dot_product(this%wiP(ic, is)%CP - this%wiP(1, is)%PC(:, 1), &
            this%secTauCapChord(:, is))
        enddo
        do i = 1, 3
          this%secChordwiseResVel(i, is) = lsq2(dot_product(this%secCP(:, is) - &
            this%wiP(1, is)%PC(:, 1), this%secTauCapChord(:, is)), xDist, &
            this%wiP(:, is)%chordwiseResVel(i))
        enddo
      enddo
    else  ! Use average of resultant velocities
      do is = 1, size(this%secChordwiseResVel, 2)
        do ic = 1, rows
          call this%wiP(ic, is)%calc_chordwiseResVel()
        enddo
        do i = 1, 3
          this%secChordwiseResVel(i, is) = sum(this%wiP(:, is)%chordwiseResVel(i))/rows
        enddo
      enddo
    endif
  end subroutine blade_calc_secChordwiseResVel

  subroutine blade_calc_secAlpha(this)
    ! Compute sec alpha using sec resultant velocity
  class(blade_class), intent(inout) :: this
    integer :: is

    call this%calc_secChordwiseResVel()

    !! Use acos() to find angle
    !do is=1,size(this%secAlpha)
    !  this%secAlpha(is)=acos(dot_product(this%secChordwiseResVel(:,is),this%secTauCapChord(:,is)) &
    !    /norm2(this%secChordwiseResVel(:,is)))
    !enddo

    ! Use atan2() to find angle
    do is = 1, size(this%secAlpha)
      this%secAlpha(is) = atan2(dot_product(this%secChordwiseResVel(:, is), this%secNormalVec(:, is)), &
        dot_product(this%secChordwiseResVel(:, is), this%secTauCapChord(:, is)))
    enddo
  end subroutine blade_calc_secAlpha

  !subroutine blade_calc_secAlpha(this)
  !  ! Compute sec alpha by interpolating local panel alpha
  !class(blade_class), intent(inout) :: this
  !  integer :: is, ic, rows
  !  real(dp), dimension(size(this%wiP,1)) :: xDist

  !  rows=size(this%wiP,1)
  !  if (rows .ge. 3) then  ! Use least squares fit to get alpha
  !    do is=1,size(this%secAlpha)
  !      do ic=1,rows
  !        xDist(ic)=dot_product(this%wiP(ic,is)%CP-this%wiP(1,is)%PC(:,1),  &
  !          this%secTauCapChord(:,is))
  !      enddo
  !      this%secAlpha(is)=lsq2(dot_product(this%secCP(:,is)-  &
  !        this%wiP(1,is)%PC(:,1),this%secTauCapChord(:,is)),xDist,this%wiP(:,is)%alpha)
  !    enddo
  !  else  ! Use average of alpha values
  !    do is=1,size(this%secAlpha)
  !      this%secAlpha(is)=sum(this%wiP(:,is)%alpha)/rows
  !    enddo
  !  endif
  !end subroutine blade_calc_secAlpha

  function getSecChordwiseLocations(this, chordwiseFraction)
    ! Get coordinates of a point located at a fraction of chord on each section
  class(blade_class), intent(inout) :: this
    real(dp), intent(in) :: chordwiseFraction
    real(dp), dimension(3, size(this%wiP, 2)) :: getSecChordwiseLocations
    integer :: is, rows

    rows = size(this%wiP, 1)

    do is = 1, size(this%wiP, 2)
      getSecChordwiseLocations(:, is) = (1._dp - chordwiseFraction)*(this%wiP(1, is)%PC(:, 4) + this%wiP(1, is)%PC(:, 1))*0.5_dp + &
        chordwiseFraction*(this%wiP(rows, is)%PC(:, 3) + this%wiP(rows, is)%PC(:, 2))*0.5_dp
    enddo
  end function getSecChordwiseLocations

  subroutine blade_burst_wake(this, rowFar, skewLimit, largeCoreRadius)
  class(blade_class), intent(inout) :: this
    real(dp), intent(in) :: skewLimit
    integer, intent(in) :: rowFar !, rowNear
    real(dp), intent(in) :: largeCoreRadius
    integer :: irow !, icol
    real(dp) :: skewVal

    !! Burst near wake
    !do icol=1,size(this%waP,2)
    !  do irow=rowNear,size(this%waP,1)
    !    call this%waP(irow,icol)%vr%burst(skewLimit)
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

  subroutine blade_dirLiftDrag(this)
  class(blade_class), intent(inout) :: this
    integer :: is
    do is = 1, size(this%wiP, 2)
      this%secDragDir(:, is) = unitVec(this%wiP(1, is)%velCPm)
      this%secLiftDir(:, is) = cross_product(this%secDragDir(:, is), this%yAxis)
    enddo
  end subroutine blade_dirLiftDrag

  subroutine sumSecToNetForces(this)
  class(blade_class), intent(inout) :: this

    this%forceInertial = sum(this%secForceInertial, 2)
    this%lift = sum(this%secLift, 2)
    this%drag = sum(this%secDrag, 2)
    this%dragProfile = sum(this%secDragProfile, 2)
    this%dragUnsteady = sum(this%secDragUnsteady, 2)
    this%dragInduced = sum(this%secDragInduced, 2)
  end subroutine sumSecToNetForces

end module blade_classdef

!------+-------------------+------|
! ++++ | MODULE DEFINITION | ++++ |
!------+-------------------+------|
module rotor_classdef
  use blade_classdef
  implicit none
  type rotor_class
    integer :: nb, ns, nc, nNwake, nFwake
    type(blade_class), allocatable, dimension(:) :: blade
    real(dp) :: Omega, omegaSlow
    real(dp), dimension(3) :: shaftAxis
    real(dp), dimension(3) :: hubCoords, cgCoords, fromCoords
    real(dp) :: radius, chord, root_cut, coningAngle
    real(dp) :: CT
    real(dp), dimension(3) :: forceInertial, lift, drag
    real(dp), dimension(3) :: dragInduced, dragProfile, dragUnsteady
    real(dp), dimension(3) :: liftUnitVec, dragUnitVec, sideUnitVec
    real(dp), dimension(3) :: controlPitch  ! theta0,thetaC,thetaS
    real(dp) :: thetaTwist
    real(dp) :: pivotLE  ! pivot location from LE [x/c]
    real(dp) :: flapHinge  ! hinge location from centre [x/R]
    real(dp), dimension(3) :: velBody, omegaBody
    real(dp), dimension(3) :: velWind, omegaWind
    real(dp) :: psi
    real(dp), dimension(3) :: pts  ! phi,theta,psi about cgCoords
    character(len=1) :: streamwiseCoreSwitch
    real(dp) :: spanwiseCore
    real(dp), allocatable, dimension(:) :: streamwiseCoreVec
    real(dp), allocatable, dimension(:, :) :: AIC, AIC_inv  ! Influence coefficient matrix
    real(dp), allocatable, dimension(:) :: gamVec, gamVecPrev, RHS
    real(dp), allocatable, dimension(:) :: airfoilSectionLimit
    real(dp), allocatable, dimension(:) :: CL0, CLa
    real(dp) :: initWakeVel, psiStart, skewLimit
    real(dp) :: turbulentViscosity
    real(dp) :: rollupStartRadius, rollupEndRadius
    integer :: symmetricTau
    integer :: rollupStart, rollupEnd
    integer :: suppressFwakeSwitch
    integer :: inflowPlotSwitch, bladeforcePlotSwitch
    integer :: gammaPlotSwitch, alphaPlotSwitch
    integer :: rowNear, rowFar
    integer :: nAirfoils
    character(len=30), allocatable, dimension(:) :: airfoilFile
    character(len=30) :: geometryFile
    real(dp) :: nonDimforceDenominator
  contains
    procedure :: getdata
    procedure :: init => rotor_init
    procedure :: deinit => rotor_deinit
    procedure :: plot3dtoblade
    procedure :: gettheta
    procedure :: getthetadot
    procedure :: move => rotor_move
    procedure :: rot_pts => rotor_rot_pts
    procedure :: rot_advance => rotor_rot_advance
    procedure :: assignshed
    procedure :: map_gam
    procedure :: age_wake
    procedure :: dissipate_wake
    procedure :: strain_wake
    procedure :: calcAIC
    procedure :: vind_bywing => rotor_vind_bywing
    procedure :: vind_bywing_boundVortices => rotor_vind_bywing_boundVortices
    procedure :: vind_bywake => rotor_vind_bywake
    procedure :: shiftwake => rotor_shiftwake
    procedure :: rollup => rotor_rollup
    procedure :: calc_force_gamma => rotor_calc_force_gamma
    procedure :: calc_force_alpha => rotor_calc_force_alpha
    procedure :: calc_force_alphaGamma => rotor_calc_force_alphaGamma
    procedure :: calc_secAlpha => rotor_calc_secAlpha
    procedure :: burst_wake => rotor_burst_wake
    procedure :: dirLiftDrag => rotor_dirLiftDrag
    procedure :: sumBladeToNetForces
  end type rotor_class

contains

  !-----+--------------------------+-----|
  ! -+- | Initialization Functions | -+- |
  !-----+--------------------------+-----|

  subroutine getdata(this, filename)
  class(rotor_class) :: this
    character(len=*), intent(in) :: filename
    integer :: i

    open (unit=12, file=filename)
    call skiplines(12, 4)
    read (12, *) this%nb, this%geometryFile
    call skiplines(12, 3)
    read (12, *) this%nc, this%ns, this%nNwake
    ! If nNwake is -ve, far wake is suppressed
    if (this%nNwake < 0) then
      this%suppressFwakeSwitch = 1
      this%nNwake = abs(this%nNwake)
    else
      this%suppressFwakeSwitch = 0
    endif

    if (this%nNwake < 2) error stop 'ERROR: Atleast 2 near wake rows mandatory'
    call skiplines(12, 4)
    read (12, *) this%hubCoords(1), this%hubCoords(2), this%hubCoords(3)
    call skiplines(12, 3)
    read (12, *) this%cgCoords(1), this%cgCoords(2), this%cgCoords(3)
    call skiplines(12, 3)
    read (12, *) this%fromCoords(1), this%fromCoords(2), this%fromCoords(3)
    call skiplines(12, 3)
    read (12, *) this%pts(1), this%pts(2), this%pts(3)
    call skiplines(12, 4)
    read (12, *) this%radius, this%root_cut, this%chord, this%coningAngle
    call skiplines(12, 4)
    read (12, *) this%Omega, this%shaftAxis(1), this%shaftAxis(2), this%shaftAxis(3)
    call skiplines(12, 3)
    read (12, *) this%controlPitch(1), this%controlPitch(2), this%controlPitch(3), this%thetaTwist
    call skiplines(12, 4)
    read (12, *) this%velBody(1), this%velBody(2), this%velBody(3) &
      , this%omegaBody(1), this%omegaBody(2), this%omegaBody(3)
    call skiplines(12, 5)
    read (12, *) this%pivotLE, this%flapHinge, this%symmetricTau
    call skiplines(12, 5)
    read (12, *) this%turbulentViscosity
    call skiplines(12, 4)
    read (12, *) this%spanwiseCore, this%streamwiseCoreSwitch
    call skiplines(12, 3)
    allocate (this%streamwiseCoreVec(this%ns + 1))
    if (this%streamwiseCoreSwitch .eq. 'i') then  ! [i]dentical
      read (12, *) this%streamwiseCoreVec(1)
      do i = 2, this%ns + 1
        this%streamwiseCoreVec(i) = this%streamwiseCoreVec(1)
      enddo
    elseif (this%streamwiseCoreSwitch .eq. 's') then  ! [s]ectional
      read (12, *) (this%streamwiseCoreVec(i), i=1, this%ns + 1)
    else
      error stop 'ERROR: Wrong input for streamwiseCoreSwitch in rotorXX.in'
    endif
    call skiplines(12, 4)
    ! Dimensional quantities
    read (12, *) this%rollupStartRadius, this%rollupEndRadius
    call skiplines(12, 3)
    read (12, *) this%initWakeVel, this%psiStart, this%skewLimit
    call skiplines(12, 5)
    read (12, *) this%dragUnitVec(1), this%dragUnitVec(2), this%dragUnitVec(3)
    call skiplines(12, 2)
    read (12, *) this%sideUnitVec(1), this%sideUnitVec(2), this%sideUnitVec(3)
    call skiplines(12, 2)
    read (12, *) this%liftUnitVec(1), this%liftUnitVec(2), this%liftUnitVec(3)
    call skiplines(12, 7)
    read (12, *) this%inflowPlotSwitch, this%bladeforcePlotSwitch
    call skiplines(12, 5)
    read (12, *) this%gammaPlotSwitch, this%alphaPlotSwitch
    call skiplines(12, 4)
    read (12, *) this%nAirfoils
    call skiplines(12, 4)
    if (this%nAirfoils .gt. 0) then
      allocate (this%airfoilSectionLimit(this%nAirfoils))
      allocate (this%airfoilFile(this%nAirfoils))
      allocate (this%CL0(this%nAirfoils))
      allocate (this%CLa(this%nAirfoils))
      do i = 1, this%nAirfoils
        read (12, *) this%airfoilSectionLimit(i), this%CL0(i), this%CLa(i), this%airfoilFile(i)
      enddo
    endif
    close (12)
  end subroutine getdata

  subroutine rotor_init(this, density, dt, nt, spanSpacingSwitch, fdSchemeSwitch)
    ! Initialize variables of rotor geometry and wake
  class(rotor_class) :: this
    real(dp), intent(in) :: density, dt
    integer, intent(in) :: nt, spanSpacingSwitch, fdSchemeSwitch

    real(dp), dimension(this%nc + 1) :: xVec
    real(dp), dimension(this%ns + 1) :: yVec
    real(dp), dimension(this%nc, this%ns) :: dx, dy
    real(dp), dimension(3) :: leftTipCP
    real(dp) :: dxdymin, secCPLoc, rbyR
    integer :: i, j, ib, is, ic
    real(dp) :: bladeOffset
    real(dp) :: velShed
    real(dp), dimension(4) :: xshift
    logical :: warnUser

    ! Initialize variables for use in allocating
    this%nNwake = min(this%nNwake, nt)
    this%nFwake = nt - this%nNwake

    ! Allocate rotor object variables
    allocate (this%blade(this%nb))
    allocate (this%AIC(this%nc*this%ns*this%nb, this%nc*this%ns*this%nb))
    allocate (this%AIC_inv(this%nc*this%ns*this%nb, this%nc*this%ns*this%nb))
    allocate (this%gamVec(this%nc*this%ns*this%nb))
    allocate (this%gamVecPrev(this%nc*this%ns*this%nb))
    allocate (this%RHS(this%nc*this%ns*this%nb))

    ! Allocate blade object variables
    do ib = 1, this%nb
      allocate (this%blade(ib)%wiP(this%nc, this%ns))
      allocate (this%blade(ib)%waP(this%nNwake, this%ns))
      allocate (this%blade(ib)%waF(this%nFwake))
      allocate (this%blade(ib)%secTauCapChord(3, this%ns))
      allocate (this%blade(ib)%secTauCapSpan(3, this%ns))
      allocate (this%blade(ib)%secNormalVec(3, this%ns))
      allocate (this%blade(ib)%secVelFreestream(3, this%ns))
      allocate (this%blade(ib)%secForceInertial(3, this%ns))
      allocate (this%blade(ib)%secLift(3, this%ns))
      allocate (this%blade(ib)%secDrag(3, this%ns))
      allocate (this%blade(ib)%secLiftDir(3, this%ns))
      allocate (this%blade(ib)%secDragDir(3, this%ns))
      allocate (this%blade(ib)%secDragInduced(3, this%ns))
      allocate (this%blade(ib)%secDragProfile(3, this%ns))
      allocate (this%blade(ib)%secDragUnsteady(3, this%ns))
      allocate (this%blade(ib)%secAlpha(this%ns))
      allocate (this%blade(ib)%airfoilNo(this%ns))
      allocate (this%blade(ib)%secCL(this%ns))
      allocate (this%blade(ib)%secCD(this%ns))
      allocate (this%blade(ib)%secCM(this%ns))
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
      if (this%Omega .ge. 0) then
        xVec = linspace(-this%chord, 0._dp, this%nc + 1)
      else
        xVec = linspace(this%chord, 0._dp, this%nc + 1)
      endif
      select case (spanSpacingSwitch)
      case (1)
        yVec = linspace(this%root_cut*this%radius, this%radius, this%ns + 1)
      case (2)
        yVec = cosspace(this%root_cut*this%radius, this%radius, this%ns + 1)
      case (3)
        yVec = halfsinspace(this%root_cut*this%radius, this%radius, this%ns + 1)
      end select

      do ib = 1, this%nb
        ! Initialize panel coordinates
        do j = 1, this%ns
          do i = 1, this%nc
            call this%blade(ib)%wiP(i, j)%assignP(1, (/xVec(i), yVec(j), 0._dp/))
            call this%blade(ib)%wiP(i, j)%assignP(2, (/xVec(i + 1), yVec(j), 0._dp/))
            call this%blade(ib)%wiP(i, j)%assignP(3, (/xVec(i + 1), yVec(j + 1), 0._dp/))
            call this%blade(ib)%wiP(i, j)%assignP(4, (/xVec(i), yVec(j + 1), 0._dp/))
          enddo
        enddo
      enddo
    else
      call this%plot3dtoblade(trim(this%geometryFile))
    endif

    do ib = 1, this%nb
      ! Initialize blade axes
      ! Can go wrong if importing geometry from OpenVSP
      ! and then using these for orientation of lift vectors
      this%blade(ib)%xAxis = xAxis
      this%blade(ib)%yAxis = yAxis
      this%blade(ib)%zAxis = zAxis

      ! Initialize sec vectors
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

      ! Assign wind velocities
      this%velWind = -1._dp*this%velBody
      this%omegaWind = -1._dp*this%omegaBody

      ! Shed last row of vortices
      if (abs(norm2(this%velWind)) < eps) then
        velShed = min(0.05*this%Omega*norm2(this%blade(ib)% &
          wiP(this%nc, this%ns)%vr%vf(2)%fc(:, 1) - this%hubCoords), 0.125_dp*this%chord/dt)
      else
        velShed = 0.3_dp*norm2(this%velWind)
      endif
      do j = 1, this%ns
        call this%blade(ib)%wiP(this%nc, j)%vr%shiftdP(2, (/velShed*dt, 0._dp, 0._dp/))
        call this%blade(ib)%wiP(this%nc, j)%vr%shiftdP(3, (/velShed*dt, 0._dp, 0._dp/))
      enddo

      ! Initialize CP coords, nCap, panelArea and pivotLE
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
        this%blade(ib)%waP%vr%vf(i)%age = 0._dp
      enddo

      this%blade(ib)%waF%vf%age = 0._dp

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

      ! Initialize all core radius of wing vortices to zero
      do i = 1, 4
        this%blade(ib)%wiP%vr%vf(i)%rVc0 = min(this%spanwiseCore, dxdymin*0.1_dp)
      enddo

      !print*,'Wing vortex core radius set to ',min(this%spanwiseCore,dxdymin*0.1_dp)/this%chord,'times chord'

      ! Initialize spanwise vortex core radius for last row of wing to that of wake
      this%blade(ib)%wiP(this%nc, :)%vr%vf(2)%rVc0 = this%spanwiseCore

      ! Initialize all current core radius of wing vortices to initial core radius
      do i = 1, 4
        this%blade(ib)%wiP%vr%vf(i)%rVc = this%blade(ib)%wiP%vr%vf(i)%rVc0
      enddo

      ! Verify CP is outside vortex core for boundary panels
      warnUser = .FALSE.
      if (isCPinsidecore(this%blade(ib)%wiP(1, 1))) then
        print *, 'Warning: CP inside vortex core at panel LU'
        warnUser = .TRUE.
      endif
      if (isCPinsidecore(this%blade(ib)%wiP(this%nc, 1))) then
        print *, 'Warning: CP inside vortex core at panel LB'
        warnUser = .TRUE.
      endif
      if (isCPinsidecore(this%blade(ib)%wiP(1, this%ns))) then
        print *, 'Warning: CP inside vortex core at panel RU'
        warnUser = .TRUE.
      endif
      if (isCPinsidecore(this%blade(ib)%wiP(this%nc, this%ns))) then
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
      call this%blade(ib)%move(-1._dp*this%fromCoords)
      call this%blade(ib)%move(this%hubCoords)
    enddo

    ! Set Coning angle
    do ib = 1, this%nb
      call this%blade(ib)%rot_axis(this%coningAngle, xAxis, (/0._dp, 0._dp, 0._dp/))
    enddo

    ! Rotate remaining blades to their positions
    ! Rotate blades for multi-bladed rotors
    do ib = 2, this%nb
      bladeOffset = 2._dp*pi/this%nb*(ib - 1)
      call this%blade(ib)%rot_axis(bladeOffset, this%shaftAxis, this%hubCoords)
    enddo

    ! Rotate rotor by phi,theta,psi about CG
    call this%rot_pts(this%pts, this%cgCoords, 1)

    ! Compute denominators for non-dimensionalisation
    if (abs(this%Omega) .gt. eps) then
      ! Rotary-wing
      this%nonDimforceDenominator = density*(pi*this%radius**2._dp)* &
        (this%radius*this%Omega)**2._dp
      ! Propeller
      !this%nonDimforceDenominator = density*(this%Omega/(2._dp*pi))**2._dp* &
      !(2._dp*this%radius)**4._dp
    else
      ! Fixed-wing
      this%nonDimforceDenominator = 0.5_dp*density*(this%radius*(1._dp - this%root_cut)* &
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
        allocate (this%blade(ib)%CL0(this%nAirfoils))
        allocate (this%blade(ib)%CLa(this%nAirfoils))
        this%blade(ib)%airfoilSectionLimit = this%airfoilSectionLimit
        do i = 1, this%nAirfoils
          this%blade(ib)%CL0(i) = this%CL0(i)
          this%blade(ib)%CLa(i) = this%CLa(i)
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
      allocate (this%blade(ib)%velNwake(3, this%nNwake, this%ns + 1))
      allocate (this%blade(ib)%velFwake(3, this%nFwake))
      this%blade(ib)%velNwake = 0._dp
      this%blade(ib)%velFwake = 0._dp

      select case (fdSchemeSwitch)
      case (0)
        ! Do nothing
      case (1)
        allocate (this%blade(ib)%waPPredicted(this%nNwake, this%ns))
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
        allocate (this%blade(ib)%waPPredicted(this%nNwake, this%ns))
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
        allocate (this%blade(ib)%waPPredicted(this%nNwake, this%ns))
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
        allocate (this%blade(ib)%waPPredicted(this%nNwake, this%ns))
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
        this%blade(ib)%waP%vr%vf(i)%rVc0 = this%spanwiseCore
        this%blade(ib)%waP%vr%vf(i)%rVc = this%spanwiseCore
      enddo

      this%blade(ib)%waP%vr%gam = 0._dp
      this%blade(ib)%waF%gam = 0._dp

      ! Assign core_radius to tip vortices
      do j = 1, this%ns
        do i = 1, this%nNwake
          this%blade(ib)%waP(i, j)%vr%vf(1)%rVc0 = this%streamwiseCoreVec(j)
          this%blade(ib)%waP(i, j)%vr%vf(1)%rVc = this%streamwiseCoreVec(j)
        enddo
      enddo

      do j = 1, this%ns
        do i = 1, this%nNwake
          this%blade(ib)%waP(i, j)%vr%vf(3)%rVc0 = this%streamwiseCoreVec(j + 1)
          this%blade(ib)%waP(i, j)%vr%vf(3)%rVc = this%streamwiseCoreVec(j + 1)
        enddo
      enddo

      !do i=1,this%nFwake
      this%blade(ib)%waF%vf%rVc0 = this%streamwiseCoreVec(this%ns + 1)
      this%blade(ib)%waF%vf%rVc = this%streamwiseCoreVec(this%ns + 1)
      !enddo

    enddo

    ! Compute direction of wind frame forces
    ! Assuming sideslip is not present
    if ((norm2(this%dragUnitVec) .le. eps) &
      .and. (norm2(this%sideUnitVec) .le. eps) &
      .and. (norm2(this%liftUnitVec) .le. eps)) then
      if (abs(this%Omega) .le. eps) then
        if (abs(this%velWind(1)) .gt. eps) then  ! v is assumed zero
          this%dragUnitVec = unitVec((/this%velWind(1), 0._dp, this%velWind(3)/))
          this%sideUnitVec = yAxis
        else  ! u is assumed zero
          this%dragUnitVec = unitVec((/0._dp, this%velWind(2), this%velWind(3)/))
          this%sideUnitVec = xAxis
        endif
        this%liftUnitVec = cross_product(this%dragUnitVec, this%sideUnitVec)
      else
        ! Drag along forward velocity direction
        if (abs(this%velWind(1)) .gt. eps) then
          this%dragUnitVec = unitVec((/this%velWind(1), 0._dp, this%velWind(3)/))
          this%sideUnitVec = yAxis
        elseif (abs(this%velWind(2)) .gt. eps) then
          this%dragUnitVec = unitVec((/this%velWind(1), 0._dp, this%velWind(3)/))
          this%sideUnitVec = xAxis
        else
          this%sideUnitVec = yAxis
          this%dragUnitVec = cross_product(this%sideUnitVec, this%shaftAxis)
        endif
        this%liftUnitVec = this%shaftAxis
      endif
    endif

  end subroutine rotor_init

  subroutine rotor_deinit(this, fdSchemeSwitch)
    ! Deinitialise rotor variables
  class(rotor_class) :: this
    integer, intent(in) :: fdSchemeSwitch
    integer :: ib
    ! Deallocate variables
    do ib = 1, this%nb
      deallocate (this%blade(ib)%velNwake)
      deallocate (this%blade(ib)%velFwake)

      select case (fdSchemeSwitch)
      case (0)
        ! Nothing to deallocate
      case (1)
        deallocate (this%blade(ib)%waPPredicted)
        deallocate (this%blade(ib)%velNwakePredicted)

        deallocate (this%blade(ib)%waFPredicted)
        deallocate (this%blade(ib)%velFwakePredicted)
      case (2)
        deallocate (this%blade(ib)%velNwake1)
        deallocate (this%blade(ib)%velNwakeStep)

        deallocate (this%blade(ib)%velFwake1)
        deallocate (this%blade(ib)%velFwakeStep)
      case (3)
        deallocate (this%blade(ib)%waPPredicted)
        deallocate (this%blade(ib)%velNwake1)
        deallocate (this%blade(ib)%velNwakeStep)

        deallocate (this%blade(ib)%waFPredicted)
        deallocate (this%blade(ib)%velFwake1)
        deallocate (this%blade(ib)%velFwakeStep)
      case (4)
        deallocate (this%blade(ib)%waPPredicted)
        deallocate (this%blade(ib)%velNwake1)
        deallocate (this%blade(ib)%velNwake2)
        deallocate (this%blade(ib)%velNwakeStep)

        deallocate (this%blade(ib)%waFPredicted)
        deallocate (this%blade(ib)%velFwake1)
        deallocate (this%blade(ib)%velFwake2)
        deallocate (this%blade(ib)%velFwakeStep)
      case (5)
        deallocate (this%blade(ib)%waPPredicted)
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

  subroutine plot3dtoblade(this, PLOT3Dfilename)
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

  end subroutine plot3dtoblade

  function gettheta(this, psi, ib)
    ! Get pitch angle corresponding to blade azimuthal location
  class(rotor_class) :: this
    real(dp), intent(in) :: psi
    integer, intent(in) :: ib
    real(dp) :: gettheta
    real(dp) :: bladeOffset

    bladeOffset = 2._dp*pi/this%nb*(ib - 1)
    gettheta = this%controlPitch(1) &
      + this%controlPitch(2)*cos(psi + bladeOffset) &
      + this%controlPitch(3)*sin(psi + bladeOffset)

  end function gettheta

  function getthetadot(this, psi, ib)
  class(rotor_class) :: this
    real(dp), intent(in) :: psi
    integer, intent(in) :: ib
    real(dp) :: getthetadot
    real(dp) :: bladeOffset

    bladeOffset = 2._dp*pi/this%nb*(ib - 1)
    getthetadot = -this%controlPitch(2)*sin(psi + bladeOffset) &
      + this%controlPitch(3)*cos(psi + bladeOffset)

  end function getthetadot

  subroutine calcAIC(this)
    ! Compute AIC matrix for rotor
  class(rotor_class), intent(inout) :: this
    integer :: ib, jblade, is, ic, i, j, row, col
    real(dp), dimension(3) :: vec_dummy

    ! Influence Coefficient Matrix
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
    this%AIC_inv = inv(this%AIC)
  end subroutine calcAIC

  subroutine map_gam(this)
    ! Map gam from vector to matrix format
  class(rotor_class), intent(inout) :: this
    integer :: ib
    do ib = 1, this%nb
      this%blade(ib)%wiP%vr%gam &
        = reshape(this%gamVec(1 + this%nc*this%ns*(ib - 1):this%nc*this%ns*ib), (/this%nc, this%ns/))
    enddo
  end subroutine map_gam

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

    this%hubCoords = this%hubCoords - origin
    this%hubCoords = matmul(TMat, this%hubCoords)
    this%hubCoords = this%hubCoords + origin

    this%cgCoords = this%cgCoords - origin
    this%cgCoords = matmul(TMat, this%cgCoords)
    this%cgCoords = this%cgCoords + origin

  end subroutine rotor_rot_pts

  subroutine rotor_rot_advance(this, dpsi)
    ! Rotate rotos by dpsi angle about axis
  class(rotor_class), intent(inout) :: this
    real(dp), intent(in) :: dpsi
    integer :: ib
    real(dp) :: dtheta

    this%psi = this%psi + dpsi
    do ib = 1, this%nb
      call this%blade(ib)%rot_axis(dpsi, this%shaftAxis, this%hubCoords)
      this%blade(ib)%psi = this%blade(ib)%psi + dpsi
      dtheta = this%gettheta(this%psi, ib) - this%blade(ib)%theta
      call this%blade(ib)%rot_pitch(dtheta)
      this%blade(ib)%theta = this%gettheta(this%psi, ib)
    enddo

  end subroutine rotor_rot_advance

  !-----+---------------------------+-----|
  ! -+- | Wake Convection Functions | -+- |
  !-----+---------------------------+-----|

  subroutine assignshed(this, edge)
    ! Assign coordinates to first rowNear of wake from last row of blade
  class(rotor_class), intent(inout) :: this
    character(len=2), intent(in) :: edge
    integer :: i, ib

    select case (edge)
    case ('LE')    ! assign to LE
      do ib = 1, this%nb
        do i = 1, this%ns
          call this%blade(ib)%waP(this%rowNear, i)%vr%assignP(1, this%blade(ib)%wiP(this%nc, i)%vr%vf(2)%fc(:, 1))
          call this%blade(ib)%waP(this%rowNear, i)%vr%assignP(4, this%blade(ib)%wiP(this%nc, i)%vr%vf(3)%fc(:, 1))
          call this%blade(ib)%waP(this%rowNear, i)%vr%calclength(.TRUE.)    ! TRUE => record original length
        enddo
        this%blade(ib)%waP(this%rowNear, :)%vr%gam = this%blade(ib)%wiP(this%nc, :)%vr%gam

      enddo
    case ('TE')    ! assign to next row's TE
      do ib = 1, this%nb
        do i = 1, this%ns
          call this%blade(ib)%waP(max(this%rowNear - 1, 1), i)%vr%assignP(2, this%blade(ib)%wiP(this%nc, i)%vr%vf(2)%fc(:, 1))
          call this%blade(ib)%waP(max(this%rowNear - 1, 1), i)%vr%assignP(3, this%blade(ib)%wiP(this%nc, i)%vr%vf(3)%fc(:, 1))
        enddo
      enddo
    case default
      error stop 'ERROR: Wrong option for edge'
    end select

  end subroutine assignshed

  !-----+----------------------------+-----|
  ! -+- | Wake Dissipation Functions | -+- |
  !-----+----------------------------+-----|

  subroutine age_wake(this, dt)
    ! Update age of wake filaments
  class(rotor_class), intent(inout) :: this
    real(dp), intent(in) :: dt
    integer :: ib, ifil
    do ib = 1, this%nb
      do ifil = 1, 4
        this%blade(ib)%waP(this%rowNear:this%nNwake, :)%vr%vf(ifil)%age = &
          this%blade(ib)%waP(this%rowNear:this%nNwake, :)%vr%vf(ifil)%age + dt
      enddo
      this%blade(ib)%waF(this%rowFar:this%nFwake)%vf%age = &
        this%blade(ib)%waF(this%rowFar:this%nFwake)%vf%age + dt
    enddo
  end subroutine age_wake

  subroutine dissipate_wake(this, dt)
  class(rotor_class), intent(inout) :: this
    real(dp), intent(in) :: dt
    real(dp) :: oseenParameter, kinematicViscosity
    integer :: ib, ic, is
    oseenParameter = 1.2564_dp
    kinematicViscosity = 0.0000181_dp

    ! Update wake age
    call this%age_wake(dt)

    ! Dissipate near wake
    do ib = 1, this%nb
      do is = 1, this%ns
        !$omp parallel do
        do ic = this%rowNear, this%nNwake
          this%blade(ib)%waP(ic, is)%vr%vf(1)%rVc = sqrt(this%blade(ib)%waP(ic, is)%vr%vf(1)%rVc**2._dp &
            + 4._dp*oseenParameter*this%turbulentViscosity*kinematicViscosity*dt)
          this%blade(ib)%waP(ic, is)%vr%vf(3)%rVc = this%blade(ib)%waP(ic, is)%vr%vf(1)%rVc
        enddo
        !$omp end parallel do

        ! To maintain consistency of rVc in overlapping filaments
        !$omp parallel do
        do ic = this%rowNear, this%nNwake
          this%blade(ib)%waP(ic, is)%vr%vf(2)%rVc = sqrt(this%blade(ib)%waP(ic, is)%vr%vf(2)%rVc**2._dp &
            + 4._dp*oseenParameter*this%turbulentViscosity*kinematicViscosity*dt)
        enddo
        !$omp end parallel do

        if (this%rowNear .ne. this%nNwake) then
          !$omp parallel do
          do ic = this%rowNear + 1, this%nNwake
            this%blade(ib)%waP(ic, is)%vr%vf(4)%rVc = this%blade(ib)%waP(ic - 1, is)%vr%vf(2)%rVc
          enddo
          !$omp end parallel do
        endif
      enddo

      ! Dissipate far wake if present
      !$omp parallel do
      do ic = this%rowFar, this%nFwake
        this%blade(ib)%waF(ic)%vf%rVc = sqrt(this%blade(ib)%waF(ic)%vf%rVc**2._dp &
          + 4._dp*oseenParameter*this%turbulentViscosity*kinematicViscosity*dt)
      enddo
      !$omp end parallel do
    enddo

  end subroutine dissipate_wake

  subroutine strain_wake(this)
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
  end subroutine strain_wake

  function rotor_vind_bywing(this, P)
    ! Compute induced velocity by all wing vortices at P
  class(rotor_class), intent(inout) :: this
    real(dp), intent(in), dimension(3) :: P
    real(dp), dimension(3) :: rotor_vind_bywing
    integer :: ib

    rotor_vind_bywing = 0._dp
    do ib = 1, this%nb
      rotor_vind_bywing = rotor_vind_bywing + this%blade(ib)%vind_bywing(P)
    enddo
  end function rotor_vind_bywing

  function rotor_vind_bywing_boundVortices(this, P)
    ! Compute induced velocity by bound vortices at P
  class(rotor_class), intent(inout) :: this
    real(dp), intent(in), dimension(3) :: P
    real(dp), dimension(3) :: rotor_vind_bywing_boundVortices
    integer :: ib

    rotor_vind_bywing_boundVortices = 0._dp
    do ib = 1, this%nb
      rotor_vind_bywing_boundVortices = rotor_vind_bywing_boundVortices + this%blade(ib)%vind_bywing_boundVortices(P)
    enddo
  end function rotor_vind_bywing_boundVortices

  function rotor_vind_bywake(this, P, optionalChar)
    ! Compute induced velocity by wake vortices at P
  class(rotor_class), intent(inout) :: this
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
        this%blade(ib)%waP(i, :) = this%blade(ib)%waP(i - 1, :)
      enddo

      ! Wake age of first row has to be set to zero
      do i = 1, 4
        this%blade(ib)%waP(1, :)%vr%vf(i)%age = 0._dp
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
      gamRollup = this%blade(ib)%waP(this%nNwake, this%ns)%vr%gam
      centroidLE = 0._dp
      centroidTE = 0._dp
      radiusRollup = 0._dp
      gamSum = 0._dp

      do ispan = this%rollupStart, this%rollupEnd
        ! Find centroid LE and TE
        centroidLE = centroidLE + this%blade(ib)%waP(this%nNwake, ispan)%vr%vf(4)%fc(:, 1)* &
          this%blade(ib)%waP(this%nNwake, ispan)%vr%gam
        centroidTE = centroidTE + this%blade(ib)%waP(this%nNwake, ispan)%vr%vf(3)%fc(:, 1)* &
          this%blade(ib)%waP(this%nNwake, ispan)%vr%gam
        gamSum = gamSum + this%blade(ib)%waP(this%nNwake, ispan)%vr%gam

        ! Assign gamRollup and radiusRollup from last row to wake filament gamma
        ! Compute gamRollup
        if (sign(1._dp, this%Omega*this%controlPitch(1)) > eps) then    ! +ve Omega or zero Omega with +ve pitch
          if (this%blade(ib)%waP(this%nNwake, ispan)%vr%gam < gamRollup) then    ! '<' because of negative gamma
            gamRollup = this%blade(ib)%waP(this%nNwake, ispan)%vr%gam
          endif
        else    ! one of Omega or pitch is negative
          if (this%blade(ib)%waP(this%nNwake, ispan)%vr%gam > gamRollup) then    ! '>' because of positive gamma
            gamRollup = this%blade(ib)%waP(this%nNwake, ispan)%vr%gam
          endif
        endif

        ! Compute radiusRollup
        radiusRollup = radiusRollup + this%blade(ib)%waP(this%nNwake, ispan)%vr%vf(3)%rVc* &
          this%blade(ib)%waP(this%nNwake, ispan)%vr%gam
      enddo

      ageRollup = this%blade(ib)%waP(this%nNwake, this%ns)%vr%vf(3)%age
      if (abs(gamSum) > eps) then
        centroidLE = centroidLE/gamSum
        centroidTE = centroidTE/gamSum
        radiusRollup = radiusRollup/gamSum
      else
        centroidLE = this%blade(ib)%waP(this%nNwake, this%rollupEnd)%vr%vf(2)%fc(:, 1)
        centroidTE = this%blade(ib)%waP(this%nNwake, this%rollupEnd)%vr%vf(3)%fc(:, 1)
        radiusRollup = this%blade(ib)%waP(this%nNwake, this%rollupEnd)%vr%vf(3)%rVc
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
      call this%blade(ib)%waF(rowFarNext)%vf%calclength(.TRUE.)    ! TRUE => record original length

      ! Ensure continuity in far wake by assigning
      ! current centroidTE to LE of previous far wake filament
      ! The discontinuity would occur due to convection of
      ! last row of waP in convectwake()
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
      call this%blade(ib)%calc_force_gamma(density, sign(1._dp, this%Omega*this%controlPitch(1)), dt)
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
      call this%blade(ib)%calc_force_alphaGamma(density, sign(1._dp, this%Omega*this%controlPitch(1)), &
        velSound, dt)
    enddo
    call this%sumBladeToNetForces()
  end subroutine rotor_calc_force_alphaGamma

  subroutine rotor_calc_secAlpha(this)
  class(rotor_class), intent(inout) :: this
    integer :: ib

    do ib = 1, this%nb
      call this%blade(ib)%calc_secAlpha()
    enddo
  end subroutine rotor_calc_secAlpha

  subroutine rotor_burst_wake(this)
  class(rotor_class), intent(inout) :: this
    integer :: ib
    do ib = 1, this%nb
      call this%blade(ib)%burst_wake(this%rowFar, this%skewLimit, this%chord)
    enddo
  end subroutine rotor_burst_wake

  subroutine rotor_dirLiftDrag(this)
  class(rotor_class), intent(inout) :: this
    integer :: ib

    do ib = 1, this%nb
      call this%blade(ib)%dirLiftDrag()
    enddo
  end subroutine rotor_dirLiftDrag

  subroutine sumBladeToNetForces(this)
  class(rotor_class), intent(inout) :: this
    integer :: ib

    this%forceInertial = 0._dp
    this%lift = 0._dp
    this%drag = 0._dp
    this%dragInduced = 0._dp
    this%dragProfile = 0._dp
    this%dragUnsteady = 0._dp

    do ib = 1, this%nb
      this%forceInertial = this%forceInertial + this%blade(ib)%forceInertial
      this%lift = this%lift + this%blade(ib)%lift
      this%drag = this%drag + this%blade(ib)%drag
      this%dragInduced = this%dragInduced + this%blade(ib)%dragInduced
      this%dragProfile = this%dragProfile + this%blade(ib)%dragProfile
      this%dragUnsteady = this%dragUnsteady + this%blade(ib)%dragUnsteady
    enddo

  end subroutine sumBladeToNetForces
end module rotor_classdef