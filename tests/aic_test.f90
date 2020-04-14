module aic_test
  use naturalfruit
  use rotor_classdef
  implicit none

  type(rotor_class) :: rotor
  real(dp) :: dt, density
  integer :: nt, spanSpacingSwitch, fdSchemeSwitch

contains

  subroutine setUp()

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

  end subroutine setUp

  subroutine test_1panel()

    call rotor%init(density, dt, nt, spanSpacingSwitch, fdSchemeSwitch)

    call testcase_initialize()
    print*,rotor%blade(1)%wiP(1,1)%CP
    call assert_true(.true.)
    call testcase_finalize()

  end subroutine test_1panel


end module aic_test
