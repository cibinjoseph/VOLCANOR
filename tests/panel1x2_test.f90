module panel1x2_test
  use naturalfruit
  use rotor_classdef
  implicit none

  type(rotor_class) :: rotor
  real(dp) :: dt, density
  integer :: nt

  integer :: spanSpacingSwitch = 2
  integer :: fdSchemeSwitch = 3

contains

  subroutine setup()
    ! Global vaariables
    nt = 1
    dt = 0.00625_dp
    density = 1.2_dp

    rotor%nb = 1
    rotor%geometryFile = '0'

    rotor%nc = 1
    rotor%ns = 2
    rotor%nNwake = 2

    rotor%hubCoords = (/0._dp, 0._dp, 0._dp/)
    rotor%cgCoords = (/0._dp, 0._dp, 0._dp/)
    rotor%fromCoords = (/0._dp, 0._dp, 0._dp/)
    rotor%pts = (/0._dp, 0._dp, 0._dp/)

    rotor%radius = 2._dp
    rotor%root_cut = 0._dp
    rotor%chord = 1._dp
    rotor%coningAngle = 0._dp

    rotor%Omega = 0._dp
    rotor%shaftAxis = (/0._dp, 0._dp, 0._dp/)
    rotor%controlPitch = (/5._dp, 0._dp, 0._dp/)
    rotor%thetaTwist = 0._dp
    rotor%velBody = (/-10._dp, 0._dp, 0._dp/)
    rotor%omegaBody = (/0._dp, 0._dp, 0._dp/)
    rotor%pivotLE = 0.25_dp
    rotor%flapHinge = 0._dp
    rotor%symmetricTau = 1._dp
    rotor%turbulentViscosity = 5._dp

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

    call rotor%init(1, density, dt, nt, spanSpacingSwitch, fdSchemeSwitch)
  end subroutine setup

  subroutine test_aic()
    real(dp), dimension(2, 2) :: AIC

    call testcase_initialize('test_aic')

    AIC(1, :) = (/1.1223476_dp, -9.2667149E-002_dp/)
    AIC(2, :) = (/-9.2667149E-002_dp, 1.1223476_dp/)

    call rotor%calcAIC()
    call assert_equal(AIC, rotor%AIC, tol, 'AICs do not match')

    call testcase_finalize()
  end subroutine test_aic

  subroutine test_force_gamma()
    integer :: is
    real(dp), dimension(2) :: delP, delDiConstant, delDiUnsteady
    real(dp), dimension(3) :: normalForce, forceInertial, secLift
    real(dp), dimension(3) :: secDragInduced, secDragUnsteady
    real(dp), dimension(3) :: liftDir, dragDir
    real(dp), dimension(3) :: lift, drag, dragUnsteady

    call testcase_initialize('test_force_gamma')

    call rotor%blade(1)%rot_pitch(rotor%controlPitch(1))

    do is = 1, rotor%ns
      rotor%blade(1)%wiP(1, is)%velCP = -1._dp * rotor%velBody
      rotor%blade(1)%wiP(1, is)%velCPm = rotor%blade(1)%wiP(1, is)%velCP
    enddo
    rotor%RHS = rotor%velBody(1) * &
      & sin(rotor%controlPitch(1))*(/1._dp, 1._dp/)

    rotor%gamVec = matmul(rotor%AIC_inv, rotor%RHS)
    rotor%gamVecPrev = 0._dp
    call rotor%map_gam()

    ! Assign freestream velocity
    rotor%blade(1)%secVelFreestream = 0._dp
    do is = 1, rotor%ns
      rotor%blade(1)%secVelFreestream(:, is) = -1._dp * rotor%velBody
    enddo

    call rotor%dirLiftDrag()

    liftDir = (/0._dp, 0._dp, 1._dp/)
    call assert_equal(rotor%blade(1)%secLiftDir(:, 1), liftDir, tol, &
      & 'secliftDir does not match')
    call assert_equal(rotor%blade(1)%secLiftDir(:, 2), liftDir, tol, &
      & 'secliftDir does not match')

    dragDir = (/1._dp, 0._dp, 0._dp/)
    call assert_equal(rotor%blade(1)%secDragDir(:, 1), dragDir, tol, &
      & 'secliftDir does not match')
    call assert_equal(rotor%blade(1)%secDragDir(:, 2), dragDir, tol, &
      & 'secliftDir does not match')

    call rotor%calc_force_gamma(density, dt)

    delP = 91.3763089754279_dp * (/1._dp, 1._dp/)
    call assert_equal(rotor%blade(1)%wiP(1, :)%delP, delP, tol, &
      & 'delP does not match')

    delDiConstant = 0.656192498415084 * (/1._dp, 1._dp/)
    call assert_equal(rotor%blade(1)%wiP(1, :)%delDiConstant, delDiConstant, &
      & tol, 'delDiConstant does not match')

    delDiUnsteady = 7.08207889718725 * (/1._dp, 1._dp/)
    call assert_equal(rotor%blade(1)%wiP(1, :)%delDiUnsteady, delDiUnsteady, &
      & tol, 'delDiUnsteady does not match')

    normalForce = (/7.96397007829292_dp, 0._dp, 91.0285945325144_dp/)
    call assert_equal(rotor%blade(1)%wiP(1, 1)%normalForce, normalForce, &
      & tol, 'normalForce does not match')
    call assert_equal(rotor%blade(1)%wiP(1, 2)%normalForce, normalForce, &
      & tol, 'normalForce does not match')
    call assert_equal(rotor%blade(1)%secForceInertial(:, 1), normalForce, &
      & tol, 'secForceInertial does not match')
    call assert_equal(rotor%blade(1)%secForceInertial(:, 2), normalForce, &
      & tol, 'secForceInertial does not match')

    secLift = (/0._dp, 0._dp, normalForce(3)/)
    call assert_equal(rotor%blade(1)%secLift(:, 1), secLift, &
      & tol, 'secLift does not match')
    call assert_equal(rotor%blade(1)%secLift(:, 2), secLift, &
      & tol, 'secLift does not match')

    secDragInduced = (/7.73827139560233_dp, 0._dp, 0._dp/)
    call assert_equal(rotor%blade(1)%secDragInduced(:, 1), secDragInduced, &
      & tol, 'secDragInduced does not match')
    call assert_equal(rotor%blade(1)%secDragInduced(:, 2), secDragInduced, &
      & tol, 'secDragInduced does not match')

    secDragUnsteady = (/7.08207889718725_dp, 0._dp, 0._dp/)
    call assert_equal(rotor%blade(1)%secDragUnsteady(:, 1), secDragUnsteady, &
      & tol, 'secDragUnsteady does not match')
    call assert_equal(rotor%blade(1)%secDragUnsteady(:, 2), secDragUnsteady, &
      & tol, 'secDragUnsteady does not match')

    call assert_equal(rotor%blade(1)%secDrag(:, 1), secDragInduced, &
      & tol, 'secDragInduced does not match')
    call assert_equal(rotor%blade(1)%secDrag(:, 2), secDragInduced, &
      & tol, 'secDragInduced does not match')

    call assert_equal(rotor%blade(1)%secCL, &
      & 1.51714324220857_dp * (/1._dp, 1._dp/), tol, 'secCL does not match')
    call assert_equal(rotor%blade(1)%secCD, &
      & 0.128971189926706_dp * (/1._dp, 1._dp/), tol, 'secCD does not match')

    ! Net forces
    forceInertial = (/15.9279401565858_dp, 0._dp, 182.057189065029_dp/)
    call assert_equal(rotor%blade(1)%forceInertial, forceInertial, &
      & tol, 'forceInertial does not match')

    lift = (/0._dp, 0._dp, 182.057189065029_dp/)
    call assert_equal(rotor%blade(1)%lift, lift, &
      & tol, 'lift does not match')

    drag = (/15.4765427912047_dp, 0._dp, 0._dp/)
    call assert_equal(rotor%blade(1)%drag, drag, &
      & tol, 'drag does not match')
    call assert_equal(rotor%blade(1)%dragInduced, drag, &
      & tol, 'dragInduced does not match')
    call assert_equal(rotor%blade(1)%dragProfile, (/0._dp, 0._dp, 0._dp/), &
      & tol, 'dragProfile does not match')

    dragUnsteady= (/14.1641577943745_dp, 0._dp, 0._dp/)
    call assert_equal(rotor%blade(1)%dragUnsteady, dragUnsteady, &
      & tol, 'dragUnsteady does not match')

    call testcase_finalize()
  end subroutine test_force_gamma


end module panel1x2_test
