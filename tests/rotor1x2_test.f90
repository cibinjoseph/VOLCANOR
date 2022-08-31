module rotor1x2_test
  use naturalfruit
  use classdef
  implicit none

  type(rotor_class) :: rotor
  type(switches_class) :: switches
  real(dp) :: dt, density
  integer :: nt

contains

  subroutine setup()
    ! Global variables
    nt = 1
    dt = -0.014
    density = 1.2_dp

    ! Switches
    switches%fdScheme = 3
    switches%wakeDissipation = 1
    switches%slowStart = 0
    switches%ntSub = 0
    switches%ntSubInit = 0
    switches%restartWriteNt = 0

    rotor%nb = 1
    rotor%geometryFile = '0'
    rotor%spanSpacing = 2
    rotor%chordSpacing = 1

    rotor%nc = 1
    rotor%ns = 2
    rotor%nNwake = 2
    rotor%surfaceType = 1
    rotor%axisymmetrySwitch = 0

    rotor%hubCoords = [0._dp, 0._dp, 0._dp]
    rotor%cgCoords = [0._dp, 0._dp, 0._dp]
    rotor%fromCoords = [0._dp, 0._dp, 0._dp]
    rotor%pts = [0._dp, 0._dp, 0._dp]

    rotor%radius = 2._dp
    rotor%root_cut = 0.5_dp
    rotor%chord = 1._dp
    rotor%preconeAngle = 0._dp
    rotor%bladeDynamicsSwitch = 0
    rotor%pitchDynamicsSwitch = 0

    rotor%Omega = 100._dp
    rotor%shaftAxis = [0._dp, 0._dp, 1._dp]
    rotor%controlPitch = [5._dp, 0._dp, 0._dp]
    rotor%thetaTwist = 0._dp
    rotor%velBody = [0._dp, 0._dp, 0._dp]
    rotor%omegaBody = [0._dp, 0._dp, 0._dp]
    rotor%pivotLE = 0.25_dp
    rotor%flapHinge = 0._dp
    rotor%symmetricTau = 0._dp
    rotor%apparentViscCoeff = 1._dp
    rotor%decayCoeff = 0._dp

    rotor%spanwiseCore = 0.04_dp
    allocate(rotor%streamwiseCoreVec(rotor%ns+1))
    rotor%streamwiseCoreVec = 0.04_dp

    rotor%rollupStartRadius = 0.75_dp
    rotor%rollupEndRadius = 1._dp
    rotor%initWakeVel = 0._dp
    rotor%psiStart = 0._dp
    rotor%skewLimit = 0.5_dp
    rotor%dragUnitVec = 0._dp
    rotor%sideUnitVec = 0._dp
    rotor%liftUnitVec = 0._dp

    rotor%nAirfoils = 0

    call rotor%init(1, density, dt, nt, switches)
    rotor%blade(1)%theta = rotor%gettheta(rotor%psiStart, 1)
    call rotor%blade(1)%rot_pitch(rotor%blade(1)%theta)
  end subroutine setup

  subroutine test_coords()
    real(dp), dimension(3, 4) :: coords
    real(dp) :: ct, st, dtVal

    ct = cos(rotor%controlPitch(1))
    st = sin(rotor%controlPitch(1))

    call testcase_initialize('test_coords')

    dtVal = 8.7964597D-004
    call assert_equal(dtVal, dt, tol, message = 'dt does not match')

    ! Panel (1, 1) PC
    coords(:, 1) = [-(0.75+0.25*ct), 1._dp, 0.25*st]
    coords(:, 2) = [-(0.75-0.75*ct), 1.0_dp, -0.75*st]
    coords(:, 3) = [-(0.75-0.75*ct), 1.5_dp, -0.75*st]
    coords(:, 4) = [-(0.75+0.25*ct), 1.5_dp, 0.25*st]

    call assert_equal(coords, rotor%blade(1)%wiP(1, 1)%pc, &
      & message = 'Panel 1, 1 PC do not match')

    ! Panel (1, 2) PC
    coords(:, 1) = [-(0.75+0.25*ct), 1.5_dp, 0.25*st]
    coords(:, 2) = [-(0.75-0.75*ct), 1.5_dp, -0.75*st]
    coords(:, 3) = [-(0.75-0.75*ct), 2._dp, -0.75*st]
    coords(:, 4) = [-(0.75+0.25*ct), 2._dp, 0.25*st]

    call assert_equal(coords, rotor%blade(1)%wiP(1, 2)%pc, &
      & message = 'Panel 1, 2 PC do not match')

    ! Panel (1, 1) vf
    coords(:, 1) = [-0.75_dp, 1._dp, 0._dp]
    coords(:, 2) = [3.71826D-003, 1._dp, -6.59418D-002]
    coords(:, 3) = [3.71826D-003, 1.5_dp, -6.59418D-002]
    coords(:, 4) = [-0.75_dp, 1.5_dp, 0._dp]

    call assert_equal(coords(:, 1), &
      & rotor%blade(1)%wiP(1, 1)%vr%vf(1)%fc(:, 1), tol, &
      & message = 'Panel 1, 1 vf1 do not match')

    call assert_equal(coords(:, 2), &
      & rotor%blade(1)%wiP(1, 1)%vr%vf(2)%fc(:, 1), tol, &
      & message = 'Panel 1, 1 vf2 do not match')

    call assert_equal(coords(:, 3), &
      & rotor%blade(1)%wiP(1, 1)%vr%vf(3)%fc(:, 1), tol, &
      & message = 'Panel 1, 1 vf3 do not match')

    call assert_equal(coords(:, 4), &
      & rotor%blade(1)%wiP(1, 1)%vr%vf(4)%fc(:, 1), tol, &
      & message = 'Panel 1, 1 vf4 do not match')

    ! Panel (1, 2) vf
    coords(:, 1) = [-0.75_dp, 1.5_dp, 0._dp]
    coords(:, 2) = [3.71826D-003, 1.5_dp, -6.59418D-002]
    coords(:, 3) = [3.71826D-003, 2.0_dp, -6.59418D-002]
    coords(:, 4) = [-0.75_dp, 2.0_dp, 0._dp]

    call assert_equal(coords(:, 1), &
      & rotor%blade(1)%wiP(1, 2)%vr%vf(1)%fc(:, 1), tol, &
      & message = 'Panel 1, 2 vf1 do not match')

    call assert_equal(coords(:, 2), &
      & rotor%blade(1)%wiP(1, 2)%vr%vf(2)%fc(:, 1), tol, &
      & message = 'Panel 1, 2 vf2 do not match')

    call assert_equal(coords(:, 3), &
      & rotor%blade(1)%wiP(1, 2)%vr%vf(3)%fc(:, 1), tol, &
      & message = 'Panel 1, 2 vf3 do not match')

    call assert_equal(coords(:, 4), &
      & rotor%blade(1)%wiP(1, 2)%vr%vf(4)%fc(:, 1), tol, &
      & message = 'Panel 1, 2 vf4 do not match')
    
    ! Panel (1, 1) CP
    call assert_equal([0.5*ct-0.75, 1.25_dp, -0.5*st], &
      & rotor%blade(1)%wiP(1, 1)%CP, &
      & message = 'Panel 1, 1 CP do not match')

    ! Panel (1, 2) CP
    call assert_equal([0.5*ct-0.75, 1.75_dp, -0.5*st], &
      & rotor%blade(1)%wiP(1, 2)%CP, &
      & message = 'Panel 1, 2 CP do not match')

    ! Ncap
    call assert_equal([st, 0._dp, ct], rotor%blade(1)%wiP(1, 1)%nCap, &
      & message = 'Panel 1, 1 nCap do not match')

    call assert_equal(rotor%blade(1)%wiP(1, 1)%nCap, &
      & rotor%blade(1)%wiP(1, 2)%nCap, &
      & message = 'Panel 1, 2 nCap do not match')

    call testcase_finalize()
  end subroutine test_coords

  subroutine test_aic()
    real(dp), dimension(2, 2) :: AIC

    call testcase_initialize('test_aic')

    AIC(1, :) = [ 1.600113_dp, -0.281091_dp]
    AIC(2, :) = [-0.281091_dp,  1.600113_dp]

    call rotor%calcAIC()
    call assert_equal(AIC, rotor%AIC, tol, 'AICs do not match')

    call testcase_finalize()
  end subroutine test_aic

  subroutine test_force()
    use libMath, only: cross_product
    integer :: is
    real(dp), dimension(3) :: forceInertial
    real(dp), dimension(3) :: lift

    call testcase_initialize('test_force')

    do is = 1, rotor%ns
      rotor%blade(1)%wiP(1, is)%velCP = -1._dp * rotor%velBody &
        & -cross_product(rotor%Omega*rotor%shaftAxis, &
        & rotor%blade(1)%wiP(1, is)%CP-rotor%hubCoords)
      rotor%blade(1)%wiP(1, is)%velCPm = rotor%blade(1)%wiP(1, is)%velCP
      rotor%blade(1)%wiP(1, is)%velCPTotal = rotor%blade(1)%wiP(1, is)%velCP
    enddo

    rotor%RHS(1) = dot_product(rotor%blade(1)%wiP(1, 1)%velCP, &
      & rotor%blade(1)%wiP(1, 1)%nCap)
    rotor%RHS(2) = dot_product(rotor%blade(1)%wiP(1, 2)%velCP, &
      & rotor%blade(1)%wiP(1, 2)%nCap)

    rotor%RHS = -1._dp*rotor%RHS

    rotor%gamVec = matmul(rotor%AIC_inv, rotor%RHS)
    rotor%gamVecPrev = 0._dp
    call rotor%map_gam()

    call assert_equal([-8.753162_dp, -11.069649_dp], rotor%gamVec, &
      & tol, 'gamVec does not match')

    call rotor%dirLiftDrag()

    call assert_equal([0._dp, 0._dp, 1._dp], &
      & rotor%blade(1)%secLiftDir(:, 1), tol, 'secLiftDir does not match')
    call assert_equal([0._dp, 0._dp, 1._dp], &
      & rotor%blade(1)%secLiftDir(:, 2), tol, &
      & 'secLiftDir does not match')

    call assert_equal([1._dp, 0._dp, 0._dp], &
      & rotor%blade(1)%secDragDir(:, 1), tol, 'secDragDir does not match')

    call assert_equal([1._dp, 0._dp, 0._dp], &
      & rotor%blade(1)%secDragDir(:, 2), tol, &
      & 'secDragDir does not match')

    call rotor%calc_secAlpha()

    call assert_equal(5._dp*pi/180._dp*[1._dp, 1._dp], &
      & rotor%blade(1)%secTheta, tol, 'secTheta does not match')

    call assert_equal(5._dp*pi/180._dp*[1._dp, 1._dp], &
      & rotor%blade(1)%secAlpha, tol, 'secTheta does not match')

    call rotor%calc_force(density, dt)

    call assert_equal([7278.445742_dp, 9866.306155_dp], &
      & rotor%blade(1)%wiP(1, :)%delP, tol, &
      & 'delP does not match')

    call assert_equal([317.179172_dp, 0._dp ,  3625.374529_dp], &
      & rotor%blade(1)%wiP(1, 1)%normalForce, &
      & tol, 'normalForce does not match')

    call assert_equal(rotor%blade(1)%wiP(1, 1)%normalForce, &
      & rotor%blade(1)%secForceInertial(:, 1), &
      & tol, 'secForceInertial does not match')

    call assert_equal([0._dp, 0._dp, &
      & rotor%blade(1)%wiP(1, 1)%normalForce(3)], &
      & rotor%blade(1)%secLift(:, 1), tol, 'secLift does not match')

    call assert_equal([429.952620_dp, 0._dp, 4914.380941_dp], &
      & rotor%blade(1)%wiP(1, 2)%normalForce, &
      & tol, 'normalForce does not match')
    call assert_equal(rotor%blade(1)%wiP(1, 2)%normalForce, &
      & rotor%blade(1)%secForceInertial(:, 2), &
      & tol, 'secForceInertial does not match')

    call assert_equal([0._dp, 0._dp, &
      & rotor%blade(1)%wiP(1, 2)%normalForce(3)], &
      & rotor%blade(1)%secLift(:, 2), tol, 'secLift does not match')

    call assert_equal([0.773413_dp, 0.534898_dp], rotor%blade(1)%secCL, &
      & tol, 'secCL does not match')

    ! Net forces
    forceInertial = rotor%blade(1)%wiP(1, 1)%normalForce + &
    & rotor%blade(1)%wiP(1, 2)%normalForce
    call assert_equal(rotor%blade(1)%forceInertial, forceInertial, &
      & tol, 'forceInertial does not match')

    lift = rotor%blade(1)%secLift(:, 1) + rotor%blade(1)%secLift(:, 2)
    call assert_equal(rotor%blade(1)%lift, lift, &
      & tol, 'lift does not match')

    call testcase_finalize()
  end subroutine test_force


end module rotor1x2_test
