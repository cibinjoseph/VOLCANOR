module wing1x3NegPitch_test
  use naturalfruit
  use classdef
  implicit none

  type(rotor_class) :: rotor
  type(switches_class) :: switches
  real(dp) :: dt, density
  integer :: nt

  integer :: spanSpacingSwitch = 2
  integer :: chordSpacingSwitch = 1
  integer :: fdSchemeSwitch = 3

contains

  subroutine setup()
    ! Global variables
    nt = 1
    dt = 0.00625_dp
    density = 1.2_dp

    rotor%nb = 1
    rotor%geometryFile = '0'
    rotor%surfaceType = 1
    rotor%axisymmetrySwitch = 0

    rotor%nc = 1
    rotor%ns = 3
    rotor%nNwake = 2

    rotor%hubCoords = [0._dp, 0._dp, 0._dp]
    rotor%cgCoords = [0._dp, 0._dp, 0._dp]
    rotor%fromCoords = [0._dp, 0._dp, 0._dp]
    rotor%pts = [0._dp, 0._dp, 0._dp]

    rotor%radius = 2._dp
    rotor%root_cut = 0._dp
    rotor%chord = 0.3_dp
    rotor%preconeAngle = 0._dp
    rotor%bladeDynamicsSwitch = 0
    rotor%pitchDynamicsSwitch = 0

    rotor%Omega = 0._dp
    rotor%shaftAxis = [0._dp, 0._dp, 0._dp]
    rotor%controlPitch = [-7._dp, 0._dp, 0._dp]
    rotor%thetaTwist = 0._dp
    rotor%velBody = [-6._dp, 0._dp, 0._dp]
    rotor%omegaBody = [0._dp, 0._dp, 0._dp]
    rotor%pivotLE = 0.25_dp
    rotor%flapHinge = 0._dp
    rotor%symmetricTau = 1._dp
    rotor%apparentViscCoeff = 5._dp
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
    rotor%spanSpacing = spanSpacingSwitch
    rotor%chordSpacing = chordSpacingSwitch
    switches%fdScheme = fdSchemeSwitch

    call rotor%init(1, density, dt, nt, switches)
  end subroutine setup

  subroutine test_aic()
    real(dp), dimension(3, 3) :: AIC

    call testcase_initialize('test_aic')

    AIC(1, :) = [ 3.18902463326559_dp    , -0.135541879282409_dp, -0.002934863377599825_dp]
    AIC(2, :) = [-0.02762714554230721_dp ,  2.97991717943654_dp , -0.02762714554230732_dp]
    AIC(3, :) = [-0.002934863377599823_dp, -0.135541879282409_dp,  3.18902463326559_dp]

    call rotor%calcAIC()
    call assert_equal(AIC, rotor%AIC, tol, 'AICs do not match')

    call testcase_finalize()
  end subroutine test_aic

  subroutine test_force()
    integer :: is
    real(dp), dimension(3) :: delP
    real(dp), dimension(3, 3) :: normalForce
    real(dp), dimension(3, 3) :: liftDir, dragDir, secLift
    real(dp), dimension(3) :: forceInertial
    real(dp), dimension(3) :: secCL
    real(dp), dimension(3) :: lift

    call testcase_initialize('test_force')

    call rotor%blade(1)%rot_pitch(rotor%controlPitch(1))

    ! Ncap
    call assert_equal([-sin(7._dp*pi/180._dp), 0._dp, cos(7._dp*pi/180._dp)], &
      & rotor%blade(1)%wiP(1, 1)%nCap, &
      & message = 'Panel 1, 1 nCap do not match')

    call assert_equal(rotor%blade(1)%wiP(1, 1)%nCap, &
      & rotor%blade(1)%wiP(1, 2)%nCap, &
      & message = 'Panel 1, 2 nCap do not match')

    call assert_equal(rotor%blade(1)%wiP(1, 1)%nCap, &
      & rotor%blade(1)%wiP(1, 3)%nCap, &
      & message = 'Panel 1, 3 nCap do not match')

    do is = 1, rotor%ns
      rotor%blade(1)%wiP(1, is)%velCP = -1._dp * rotor%velBody
      rotor%blade(1)%wiP(1, is)%velCPm = rotor%blade(1)%wiP(1, is)%velCP
      rotor%blade(1)%wiP(1, is)%velCPTotal = rotor%blade(1)%wiP(1, is)%velCP
      rotor%RHS(is) = dot_product(rotor%blade(1)%wiP(1, is)%velCP, &
        & rotor%blade(1)%wiP(1, is)%nCap)
    enddo

    rotor%RHS = -1._dp * rotor%RHS
    rotor%gamVec = matmul(rotor%AIC_inv, rotor%RHS)
    rotor%gamVecPrev = 0._dp
    call rotor%map_gam()

    call assert_equal([0.240131_dp, 0.249833_dp, 0.240131_dp], &
      & rotor%gamVec, tol, 'gamVec mismatch')

    call rotor%dirLiftDrag()

    liftDir(:, 1) = [0._dp, 0._dp, 1._dp]
    liftDir(:, 2) = liftDir(:, 1)
    liftDir(:, 3) = liftDir(:, 1)
    call assert_equal(rotor%blade(1)%secLiftDir, liftDir, tol, &
      & 'secLiftDir does not match')

    dragDir(:, 1) = [1._dp, 0._dp, 0._dp]
    dragDir(:, 2) = dragDir(:, 1)
    dragDir(:, 3) = dragDir(:, 1)
    call assert_equal(rotor%blade(1)%secDragDir, dragDir, tol, &
      & 'secDragDir does not match')

    call rotor%calc_force(density, dt)

    delP = [-28.7727659410054_dp, -29.9353746086400_dp, -28.7727659410054_dp]
    call assert_equal(rotor%blade(1)%wiP(1, :)%delP, delP, tol, &
      & 'delP does not match')

    normalForce(:, 1) = [0.525977_dp, 0.0_dp, -4.283744_dp]
    normalForce(:, 2) = [1.094461_dp, 0.0_dp, -8.913672_dp]
    normalForce(:, 3) = [0.525977_dp, 0.0_dp, -4.283744_dp]
    call assert_equal(normalForce(:, 1), &
      & rotor%blade(1)%wiP(1, 1)%normalForce, tol, 'normalForce mismatch')
    call assert_equal(normalForce(:, 2), &
      & rotor%blade(1)%wiP(1, 2)%normalForce, tol, 'normalForce mismatch')
    call assert_equal(normalForce(:, 3), &
      & rotor%blade(1)%wiP(1, 3)%normalForce, tol, 'normalForce mismatch')

    call assert_equal(normalForce(:, 1), &
      & rotor%blade(1)%secForceInertial(:, 1), tol, 'secForceInertial mismatch')
    call assert_equal(normalForce(:, 2), &
      & rotor%blade(1)%secForceInertial(:, 2), tol, 'secForceInertial mismatch')
    call assert_equal(normalForce(:, 3), &
      & rotor%blade(1)%secForceInertial(:, 3), tol, 'secForceInertial mismatch')


    secLift(:, 1) = [0._dp, 0._dp, normalForce(3, 1)]
    secLift(:, 2) = [0._dp, 0._dp, normalForce(3, 2)]
    secLift(:, 3) = [0._dp, 0._dp, normalForce(3, 1)]
    call assert_equal(rotor%blade(1)%secLift, secLift, &
      & tol, 'secLift does not match')

    secCL = [-1.322143_dp, -1.375566_dp, -1.322143_dp]

    call assert_equal(rotor%blade(1)%secCL, secCL, tol, 'secCL does not match')

    ! Net forces
    forceInertial = [2.146416_dp, 0.0_dp, -17.481161_dp]
    call assert_equal(rotor%blade(1)%forceInertial, forceInertial, &
      & tol, 'forceInertial does not match')

    lift = [0._dp, 0._dp, -17.481161_dp]
    call assert_equal(rotor%blade(1)%lift, lift, &
      & tol, 'lift does not match')

    call testcase_finalize()
  end subroutine test_force


end module wing1x3NegPitch_test
