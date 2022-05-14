module wing1x2_test
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
    ! Global vaariables
    nt = 1
    dt = 0.00625_dp
    density = 1.2_dp

    rotor%nb = 1
    rotor%geometryFile = '0'

    rotor%nc = 1
    rotor%ns = 2
    rotor%nNwake = 2
    rotor%surfaceType = 1
    rotor%axisymmetrySwitch = 0

    rotor%hubCoords = (/0._dp, 0._dp, 0._dp/)
    rotor%cgCoords = (/0._dp, 0._dp, 0._dp/)
    rotor%fromCoords = (/0._dp, 0._dp, 0._dp/)
    rotor%pts = (/0._dp, 0._dp, 0._dp/)

    rotor%radius = 2._dp
    rotor%root_cut = 0._dp
    rotor%chord = 1._dp
    rotor%preconeAngle = 0._dp
    rotor%bladeDynamicsSwitch = 0
    rotor%pitchDynamicsSwitch = 0

    rotor%Omega = 0._dp
    rotor%shaftAxis = (/0._dp, 0._dp, 0._dp/)
    rotor%controlPitch = (/5._dp, 0._dp, 0._dp/)
    rotor%thetaTwist = 0._dp
    rotor%velBody = (/-10._dp, 0._dp, 0._dp/)
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

  subroutine test_coords()
    real(dp), dimension(3, 4) :: coords

    call testcase_initialize('test_coords')

    ! wing (1, 1) PC
    coords(:, 1) = (/-1._dp, 0._dp, 0._dp/)
    coords(:, 2) = (/ 0._dp, 0._dp, 0._dp/)
    coords(:, 3) = (/ 0._dp, 1._dp, 0._dp/)
    coords(:, 4) = (/-1._dp, 1._dp, 0._dp/)

    call assert_equal(coords, rotor%blade(1)%wiP(1, 1)%pc, &
      & message = 'wing 1, 1 PC do not match')

    ! wing (1, 2) PC
    coords(:, 1) = (/-1._dp, 1._dp, 0._dp/)
    coords(:, 2) = (/ 0._dp, 1._dp, 0._dp/)
    coords(:, 3) = (/ 0._dp, 2._dp, 0._dp/)
    coords(:, 4) = (/-1._dp, 2._dp, 0._dp/)

    call assert_equal(coords, rotor%blade(1)%wiP(1, 2)%pc, &
      & message = 'wing 1, 2 PC do not match')

    ! wing (1, 1) vf
    coords(:, 1) = (/-0.75_dp, 0._dp, 0._dp/)
    coords(:, 2) = (/ 0.01875_dp, 0._dp, 0._dp/)
    coords(:, 3) = (/ 0.01875_dp, 1._dp, 0._dp/)
    coords(:, 4) = (/-0.75_dp, 1._dp, 0._dp/)

    call assert_equal(coords(:, 1), &
      & rotor%blade(1)%wiP(1, 1)%vr%vf(1)%fc(:, 1), &
      & message = 'wing 1, 1 vf1 do not match')

    call assert_equal(coords(:, 2), &
      & rotor%blade(1)%wiP(1, 1)%vr%vf(2)%fc(:, 1), &
      & message = 'wing 1, 1 vf2 do not match')

    call assert_equal(coords(:, 3), &
      & rotor%blade(1)%wiP(1, 1)%vr%vf(3)%fc(:, 1), &
      & message = 'wing 1, 1 vf3 do not match')

    call assert_equal(coords(:, 4), &
      & rotor%blade(1)%wiP(1, 1)%vr%vf(4)%fc(:, 1), &
      & message = 'wing 1, 1 vf4 do not match')

    ! wing (1, 2) vf
    coords(:, 1) = (/-0.75_dp, 1._dp, 0._dp/)
    coords(:, 2) = (/ 0.01875_dp, 1._dp, 0._dp/)
    coords(:, 3) = (/ 0.01875_dp, 2._dp, 0._dp/)
    coords(:, 4) = (/-0.75_dp, 2._dp, 0._dp/)

    call assert_equal(coords(:, 1), &
      & rotor%blade(1)%wiP(1, 2)%vr%vf(1)%fc(:, 1), &
      & message = 'wing 1, 2 vf1 do not match')

    call assert_equal(coords(:, 2), &
      & rotor%blade(1)%wiP(1, 2)%vr%vf(2)%fc(:, 1), &
      & message = 'wing 1, 2 vf2 do not match')

    call assert_equal(coords(:, 3), &
      & rotor%blade(1)%wiP(1, 2)%vr%vf(3)%fc(:, 1), &
      & message = 'wing 1, 2 vf3 do not match')

    call assert_equal(coords(:, 4), &
      & rotor%blade(1)%wiP(1, 2)%vr%vf(4)%fc(:, 1), &
      & message = 'wing 1, 2 vf4 do not match')

    ! wing (1, 1) CP
    call assert_equal((/-0.25_dp, 0.5_dp, 0._dp/), &
      & rotor%blade(1)%wiP(1, 1)%CP, &
      & message = 'wing 1, 1 CP do not match')

    ! wing (1, 2) CP
    call assert_equal((/-0.25_dp, 1.5_dp, 0._dp/), &
      & rotor%blade(1)%wiP(1, 2)%CP, &
      & message = 'wing 1, 2 vf4 do not match')

    call testcase_finalize()
  end subroutine test_coords

  subroutine test_aic()
    real(dp), dimension(2, 2) :: AIC

    call testcase_initialize('test_aic')

    AIC(1, :) = (/1.1223476_dp, -9.2667149E-002_dp/)
    AIC(2, :) = (/-9.2667149E-002_dp, 1.1223476_dp/)

    call rotor%calcAIC()
    call assert_equal(AIC, rotor%AIC, tol, 'AICs do not match')

    call testcase_finalize()
  end subroutine test_aic

  subroutine test_force()
    integer :: is
    real(dp), dimension(2) :: delP
    real(dp), dimension(3) :: normalForce, forceInertial, secLift
    real(dp), dimension(3) :: liftDir, dragDir
    real(dp), dimension(3) :: lift

    call testcase_initialize('test_force')

    call rotor%blade(1)%rot_pitch(rotor%controlPitch(1))

    do is = 1, rotor%ns
      rotor%blade(1)%wiP(1, is)%velCP = -1._dp * rotor%velBody
      rotor%blade(1)%wiP(1, is)%velCPm = rotor%blade(1)%wiP(1, is)%velCP
      rotor%blade(1)%wiP(1, is)%velCPTotal = rotor%blade(1)%wiP(1, is)%velCP
    enddo
    rotor%RHS = rotor%velBody(1) * &
      & sin(rotor%controlPitch(1))*(/1._dp, 1._dp/)

    rotor%gamVec = matmul(rotor%AIC_inv, rotor%RHS)
    rotor%gamVecPrev = 0._dp
    call rotor%map_gam()

    call rotor%dirLiftDrag()

    liftDir = (/0._dp, 0._dp, 1._dp/)
    call assert_equal(rotor%blade(1)%secLiftDir(:, 1), liftDir, tol, &
      & 'secLiftDir does not match')
    call assert_equal(rotor%blade(1)%secLiftDir(:, 2), liftDir, tol, &
      & 'secLiftDir does not match')

    dragDir = (/1._dp, 0._dp, 0._dp/)
    call assert_equal(rotor%blade(1)%secDragDir(:, 1), dragDir, tol, &
      & 'secliftDir does not match')
    call assert_equal(rotor%blade(1)%secDragDir(:, 2), dragDir, tol, &
      & 'secliftDir does not match')

    call rotor%calc_secAlpha()

    call assert_equal(5._dp*pi/180._dp*[1._dp, 1._dp], &
      & rotor%blade(1)%secTheta, tol, 'secTheta does not match')

    call assert_equal(5._dp*pi/180._dp*[1._dp, 1._dp], &
      & rotor%blade(1)%secAlpha, tol, 'secTheta does not match')

    call rotor%calc_force(density, dt)

    delP = 91.3763089754279_dp * (/1._dp, 1._dp/)
    call assert_equal(rotor%blade(1)%wiP(1, :)%delP, delP, tol, &
      & 'delP does not match')
    
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

    call assert_equal(secLift, rotor%blade(1)%secLiftOutPlane(:, 1), &
      & tol, 'secLiftOutPlane does not match')
    call assert_equal(secLift, rotor%blade(1)%secLiftOutPlane(:, 2), &
      & tol, 'secLiftOutPlane does not match')

    call assert_equal([0._dp, 0._dp, 0._dp], &
      & rotor%blade(1)%secLiftInPlane(:, 1), &
      & tol, 'secLiftInPlane does not match')
    call assert_equal([0._dp, 0._dp, 0._dp], &
      & rotor%blade(1)%secLiftInPlane(:, 2), &
      & tol, 'secLiftInPlane does not match')

    call assert_equal(rotor%blade(1)%secCL, &
      & 1.51714324220857_dp * (/1._dp, 1._dp/), tol, 'secCL does not match')

    ! Net forces
    forceInertial = (/15.9279401565858_dp, 0._dp, 182.057189065029_dp/)
    call assert_equal(rotor%blade(1)%forceInertial, forceInertial, &
      & tol, 'forceInertial does not match')

    lift = (/0._dp, 0._dp, 182.057189065029_dp/)
    call assert_equal(rotor%blade(1)%lift, lift, &
      & tol, 'lift does not match')

    call testcase_finalize()
  end subroutine test_force


end module wing1x2_test
