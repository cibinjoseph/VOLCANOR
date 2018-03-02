module vf_classdef
  use mymathlib
  implicit none

  type vf_class
    real(dp), dimension(3,2) :: fc  ! filament coords (xyz,1:2)
    real(dp) :: l0=0._dp    ! original length
    real(dp) :: lc=0._dp    ! current length
    real(dp) :: r_vc0=0._dp ! initial vortex core radius
    real(dp) :: r_vc=0._dp  ! current vortex core radius
    real(dp) :: age=0._dp   ! vortex age (in s)
  contains
    procedure :: vind => vfclass_vind   
    procedure :: calclength => vfclass_calclength
    procedure :: strain => vfclass_strain
  end type vf_class

  ! Contains switches for multiple features
  include "switches.f90"

  real(dp), parameter :: tol=1.E-6

contains

  function vfclass_vind(this,P) result(vind)
    ! Calculates induced velocity by unit gam vortex filament
  class(vf_class) :: this
    real(dp), dimension(3) :: vind, P
    real(dp) :: r1_r2_abs2, r1_abs, r2_abs, h, Kv
    real(dp), dimension(3) :: r1, r2, r0, r1_r2


    r1=P-this%fc(:,1)
    r2=P-this%fc(:,2)
    r0=r1-r2

    r1_r2=cross3(r1,r2)
    r1_r2_abs2=dot_product(r1_r2,r1_r2)

    r1_abs=norm2(r1)
    r2_abs=norm2(r2)

    vind=0.

    if (r1_r2_abs2>tol) then
      !vind = r1_r2/(4._dp*pi*r1_r2_abs2)*(dot_product(r0,r1)/r1_abs-dot_product(r0,r2)/r2_abs)
      vind=r1_r2/(4._dp*pi*r1_r2_abs2)*dot_product(r0,r1/r1_abs-r2/r2_abs)
    endif

    h=norm2(r1_r2)/norm2(r0)
    Kv=(h*h)/sqrt((this%r_vc**4._dp)+(h**4._dp))    
    if (this%r_vc < eps) Kv=1._dp

    vind=Kv*vind

  end function vfclass_vind

  subroutine vfclass_calclength(this,isoriginal) 
  class(vf_class) :: this
    logical, intent(in) :: isoriginal
    real(dp), dimension(3) :: delta

    delta = this%fc(:,1)-this%fc(:,2)
    this%lc=norm2(delta)
    if (isoriginal .eqv. .TRUE.) this%l0=norm2(delta)
  end subroutine vfclass_calclength

  subroutine vfclass_strain(this)
  class(vf_class) :: this
    this%r_vc=this%r_vc0*sqrt(this%l0/this%lc)
  end subroutine vfclass_strain

end module vf_classdef
