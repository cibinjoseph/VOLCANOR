module panel1x3NegPitch_test
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

    rotor%hubCoords = (/0._dp, 0._dp, 0._dp/)
    rotor%cgCoords = (/0._dp, 0._dp, 0._dp/)
    rotor%fromCoords = (/0._dp, 0._dp, 0._dp/)
    rotor%pts = (/0._dp, 0._dp, 0._dp/)

    rotor%radius = 2._dp
    rotor%root_cut = 0._dp
    rotor%chord = 0.3_dp
    rotor%preconeAngle = 0._dp
    rotor%bladeDynamicsSwitch = 0
    rotor%pitchDynamicsSwitch = 0

    rotor%Omega = 0._dp
    rotor%shaftAxis = (/0._dp, 0._dp, 0._dp/)
    rotor%controlPitch = (/-7._dp, 0._dp, 0._dp/)
    rotor%thetaTwist = 0._dp
    rotor%velBody = (/-6._dp, 0._dp, 0._dp/)
    rotor%omegaBody = (/0._dp, 0._dp, 0._dp/)
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
    switches%spanSpacing = spanSpacingSwitch
    switches%chordSpacing = chordSpacingSwitch
    switches%fdScheme = fdSchemeSwitch

    call rotor%init(1, density, dt, nt, switches)
  end subroutine setup

  subroutine test_aic()
    real(dp), dimension(3, 3) :: AIC

    call testcase_initialize('test_aic')

    AIC(1, :) = (/ 3.18902463326559_dp    , -0.135541879282409_dp, -0.002934863377599825_dp/)
    AIC(2, :) = (/-0.02762714554230721_dp ,  2.97991717943654_dp , -0.02762714554230732_dp/)
    AIC(3, :) = (/-0.002934863377599823_dp, -0.135541879282409_dp,  3.18902463326559_dp/)

    call rotor%calcAIC()
    call assert_equal(AIC, rotor%AIC, tol, 'AICs do not match')

    call testcase_finalize()
  end subroutine test_aic

  subroutine test_force_gamma()
    integer :: is
    real(dp), dimension(3) :: delP, delDiConstant, delDiUnsteady, forceInertial
    real(dp), dimension(3, 3) :: normalForce
    real(dp), dimension(3, 3) :: secDragInduced, secDragUnsteady
    real(dp), dimension(3, 3) :: liftDir, dragDir, secLift
    real(dp), dimension(3) :: secCL, secCD
    real(dp), dimension(3) :: lift, drag, dragUnsteady

    call testcase_initialize('test_force_gamma')

    call rotor%blade(1)%rot_pitch(rotor%controlPitch(1))

    do is = 1, rotor%ns
      rotor%blade(1)%wiP(1, is)%velCP = -1._dp * rotor%velBody
      rotor%blade(1)%wiP(1, is)%velCPm = rotor%blade(1)%wiP(1, is)%velCP
      rotor%blade(1)%wiP(1, is)%velCPTotal = rotor%blade(1)%wiP(1, is)%velCP
    enddo
    rotor%RHS = rotor%velBody(1) * &
      & sin(rotor%controlPitch(1))*(/1._dp, 1._dp, 1._dp/)

    rotor%gamVec = matmul(rotor%AIC_inv, rotor%RHS)
    rotor%gamVecPrev = 0._dp
    call rotor%map_gam()

    call rotor%dirLiftDrag()

    liftDir(:, 1) = (/0._dp, 0._dp, 1._dp/)
    liftDir(:, 2) = liftDir(:, 1)
    liftDir(:, 3) = liftDir(:, 1)
    call assert_equal(rotor%blade(1)%secLiftDir, liftDir, tol, &
      & 'secLiftDir does not match')

    dragDir(:, 1) = (/1._dp, 0._dp, 0._dp/)
    dragDir(:, 2) = dragDir(:, 1)
    dragDir(:, 3) = dragDir(:, 1)
    call assert_equal(rotor%blade(1)%secDragDir, dragDir, tol, &
      & 'secDragDir does not match')

    call rotor%calc_force_gamma(density, dt)

    delP = (/-28.7727659410054_dp, -29.9353746086400_dp, -28.7727659410054_dp/)
    call assert_equal(rotor%blade(1)%wiP(1, :)%delP, delP, tol, &
      & 'delP does not match')

    delDiConstant = (/0.07069535852052516_dp, &
      & 0.139675246588070_dp, 0.07069535852052521_dp/)
    call assert_equal(rotor%blade(1)%wiP(1, :)%delDiConstant, delDiConstant, &
      & tol, 'delDiConstant does not match')

    delDiUnsteady = (/0.421410397020871_dp, &
      & 0.876876288130308_dp, 0.421410397020871_dp/)
    call assert_equal(rotor%blade(1)%wiP(1, :)%delDiUnsteady, delDiUnsteady, &
      & tol, 'delDiUnsteady does not match')

    normalForce(:, 1) = (/-0.525977713977048_dp, 0.0_dp, 4.28374471602321_dp/)
    normalForce(:, 2) = (/-1.09446133444262_dp, 0.0_dp, 8.91367225972409_dp/)
    normalForce(:, 3) = (/-0.525977713977048_dp, 0.0_dp, 4.28374471602321_dp/)
    call assert_equal(rotor%blade(1)%wiP(1, 1)%normalForce, normalForce(:, 1), &
      & tol, 'normalForce does not match')
    call assert_equal(rotor%blade(1)%wiP(1, 2)%normalForce, normalForce(:, 2), &
      & tol, 'normalForce does not match')
    call assert_equal(rotor%blade(1)%wiP(1, 3)%normalForce, normalForce(:, 3), &
      & tol, 'normalForce does not match')

    call assert_equal(rotor%blade(1)%secForceInertial(:, 1), normalForce(:, 1), &
      & tol, 'secForceInertial does not match')
    call assert_equal(rotor%blade(1)%secForceInertial(:, 2), normalForce(:, 2), &
      & tol, 'secForceInertial does not match')
    call assert_equal(rotor%blade(1)%secForceInertial(:, 3), normalForce(:, 3), &
      & tol, 'secForceInertial does not match')


    secLift(:, 1) = (/0._dp, 0._dp, normalForce(3, 1)/)
    secLift(:, 2) = (/0._dp, 0._dp, normalForce(3, 2)/)
    secLift(:, 3) = (/0._dp, 0._dp, normalForce(3, 1)/)
    call assert_equal(rotor%blade(1)%secLift, secLift, &
      & tol, 'secLift does not match')

    secDragInduced(:, 1) = (/0.492105755541396_dp, 0._dp, 0._dp/)
    secDragInduced(:, 2) = (/1.01655153471838_dp , 0._dp, 0._dp/)
    secDragInduced(:, 3) = (/0.492105755541396_dp, 0._dp, 0._dp/)

    call assert_equal(rotor%blade(1)%secDragInduced, secDragInduced, &
      & tol, 'secDragInduced does not match')

    secDragUnsteady(:, 1) = (/delDiUnsteady(1), 0._dp, 0._dp/)
    secDragUnsteady(:, 2) = (/delDiUnsteady(2), 0._dp, 0._dp/)
    secDragUnsteady(:, 3) = (/delDiUnsteady(3), 0._dp, 0._dp/)

    call assert_equal(rotor%blade(1)%secDragUnsteady, secDragUnsteady, &
      & tol, 'secDragUnsteady does not match')

    call assert_equal(rotor%blade(1)%secDrag, secDragInduced, &
      & tol, 'secDragInduced does not match')

    secCL = (/-1.32214343087136_dp, -1.37556670674755_dp, -1.32214343087136_dp/)
    secCD = (/0.151884492451048_dp, 0.156875236839256_dp, 0.151884492451048_dp/)

    call assert_equal(rotor%blade(1)%secCL, secCL, tol, 'secCL does not match')
    call assert_equal(rotor%blade(1)%secCD, secCD, tol, 'secCD does not match')

    ! Net forces
    forceInertial = (/-2.14641676239672_dp, 0.0_dp, 17.4811616917705_dp/)
    call assert_equal(rotor%blade(1)%forceInertial, forceInertial, &
      & tol, 'forceInertial does not match')

    lift = (/0._dp, 0._dp, 17.4811616917705_dp/)
    call assert_equal(rotor%blade(1)%lift, lift, &
      & tol, 'lift does not match')

    drag = (/2.00076304580117_dp, 0._dp, 0._dp/)
    call assert_equal(rotor%blade(1)%drag, drag, &
      & tol, 'drag does not match')
    call assert_equal(rotor%blade(1)%dragInduced, drag, &
      & tol, 'dragInduced does not match')
    call assert_equal(rotor%blade(1)%dragProfile, (/0._dp, 0._dp, 0._dp/), &
      & tol, 'dragProfile does not match')

    dragUnsteady= (/1.71969708217205_dp, 0._dp, 0._dp/)
    call assert_equal(rotor%blade(1)%dragUnsteady, dragUnsteady, &
      & tol, 'dragUnsteady does not match')

    call testcase_finalize()
  end subroutine test_force_gamma


end module panel1x3NegPitch_test
