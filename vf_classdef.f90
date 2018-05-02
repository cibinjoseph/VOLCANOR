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

  real(dp), parameter :: tol=1.E-6
  real(dp), parameter :: inv_tol2=1.E06
  real(dp), parameter :: inv4pi=0.25_dp/pi

contains

  ! Efficient implementation to vind calculation
  function vfclass_vind(this,P) result(vind)
    ! Calculates induced velocity by unit gam vortex filament
  class(vf_class) :: this
    real(dp), dimension(3) :: vind, P
    real(dp) :: inv_r1_r2_abs2, r1_abs, r2_abs, inv_h2, Kv
    real(dp), dimension(3) :: r1, r2, r0, r1_r2

    r1=P-this%fc(:,1)
    r2=P-this%fc(:,2)
    r0=r1-r2

    ! Cross product (inlined to avoid function call)
    r1_r2(1) = r1(2)*r2(3)-r1(3)*r2(2)
    r1_r2(2) = r1(3)*r2(1)-r1(1)*r2(3)
    r1_r2(3) = r1(1)*r2(2)-r1(2)*r2(1)

    r1_abs=norm2(r1)
    r2_abs=norm2(r2)

    vind=0.

    ! Ideal vortex model (Common part)
    if (abs(sum(r1_r2))>eps) then
      inv_r1_r2_abs2=1._dp/(r1_r2(1)**2._dp+r1_r2(2)**2._dp+r1_r2(3)**2._dp)
      vind=r1_r2*inv4pi*inv_r1_r2_abs2*dot_product(r0,r1/r1_abs-r2/r2_abs)

      ! Rankine vortex model
      inv_h2=dot_product(r0,r0)*inv_r1_r2_abs2
      Kv=1._dp/sqrt(1._dp+this%r_vc**4._dp*inv_h2*inv_h2)

      vind=min(Kv,1._dp)*vind
    endif

  end function vfclass_vind

  !  function vfclass_vind(this,P) result(vind)
  !    ! Calculates induced velocity by unit gam vortex filament
  !  class(vf_class) :: this
  !    real(dp), dimension(3) :: vind, P
  !    real(dp) :: r1_r2_abs2, r1_abs, r2_abs, h, Kv
  !    real(dp), dimension(3) :: r1, r2, r0, r1_r2
  !
  !
  !    r1=P-this%fc(:,1)
  !    r2=P-this%fc(:,2)
  !    r0=r1-r2
  !
  !    r1_r2=cross3(r1,r2)
  !    r1_r2_abs2=dot_product(r1_r2,r1_r2)
  !
  !    r1_abs=norm2(r1)
  !    r2_abs=norm2(r2)
  !
  !    vind=0.
  !
  !    ! Ideal vortex model (Common part)
  !    if (r1_r2_abs2>tol) then
  !      !vind = r1_r2/(4._dp*pi*r1_r2_abs2)*(dot_product(r0,r1)/r1_abs-dot_product(r0,r2)/r2_abs)
  !      vind=r1_r2*inv4pi/r1_r2_abs2*dot_product(r0,r1/r1_abs-r2/r2_abs)
  !    endif
  !
  !    ! Rankine vortex model
  !    h=norm2(r1_r2)/norm2(r0)
  !    Kv=(h*h)/sqrt((this%r_vc**4._dp)+(h**4._dp))    
  !    !Using min function instead of if condition
  !    ! if (this%r_vc < eps) Kv=1._dp
  !    Kv=min(Kv,1._dp)
  !
  !    vind=Kv*vind
  !
  !  end function vfclass_vind

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
