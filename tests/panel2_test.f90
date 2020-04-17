module panel2_test
  use naturalfruit
  use rotor_classdef
  implicit none

  type(rotor_class) :: rotor
  real(dp) :: dt, density
  integer :: nt, spanSpacingSwitch, fdSchemeSwitch

contains

  subroutine setup()
    ! Global vaariables
    nt = 5
    dt = 0.00625_dp
    density = 1._dp

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
  end subroutine setup

  subroutine test_aic()
    integer :: spanSpacingSwitch = 2
    integer :: fdSchemeSwitch = 3
    real(dp), dimension(2, 2) :: AIC

    call testcase_initialize('test_aic')

    call rotor%init(density, dt, nt, spanSpacingSwitch, fdSchemeSwitch)
    AIC(1, :) = (/1.1223476_dp, -9.2667149E-002_dp/)
    AIC(2, :) = (/-9.2667149E-002_dp, 1.1223476_dp/)

    call rotor%calcAIC()
    call assert_equal(AIC, rotor%AIC, tol, 'AICs do not match')

    call testcase_finalize()
  end subroutine test_aic


end module panel2_test
