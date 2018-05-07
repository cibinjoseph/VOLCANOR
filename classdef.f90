!------+-------------------+------|
! ++++ | MODULE DEFINITION | ++++ |
!------+-------------------+------|
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


!------+-------------------+------|
! ++++ | MODULE DEFINITION | ++++ |
!------+-------------------+------|
module vr_classdef
  use vf_classdef
  implicit none
  type vr_class
    type(vf_class), dimension(4) :: vf
    real(dp) :: gam
  contains
    procedure :: vind => vrclass_vind
    procedure :: assignP => vrclass_assignP
    procedure :: shiftdP  => vrclass_shiftdP
    procedure :: rot  => vrclass_rot
    procedure :: calclength => vrclass_calclength
    procedure :: strain => vrclass_strain
  end type vr_class

contains

  function vrclass_vind(this,P) result(vind)
    ! Calculates induced velocity by unit gam vortex ring
  class(vr_class) :: this
    real(dp), dimension(3) :: P, vind
    real(dp), dimension(4,3) :: vind_mat
    integer :: i

    vind=0._dp

    do i=1,4
      vind_mat(i,:)=this%vf(i)%vind(P) 
    enddo
    vind=sum(vind_mat,1)

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

  subroutine vrclass_assignP(this,n,P)   ! for assigning coordinates to nth corner
  class(vr_class) :: this
    integer, intent(in) :: n
    real(dp), dimension(3) :: P

    select case (n)
    case (1)
      this%vf(4)%fc(:,2)=P
      this%vf(1)%fc(:,1)=P
    case (2)
      this%vf(1)%fc(:,2)=P
      this%vf(2)%fc(:,1)=P
    case (3)
      this%vf(2)%fc(:,2)=P
      this%vf(3)%fc(:,1)=P
    case (4)
      this%vf(3)%fc(:,2)=P
      this%vf(4)%fc(:,1)=P
    case default
      error stop 'n may only take values 1,2,3 or 4'
    end select

  end subroutine vrclass_assignP

  subroutine vrclass_shiftdP(this,n,dshift)   ! for shifting coordinates of nth corner by dshift distance (usually for Udt convection)
  class(vr_class) :: this
    integer, intent(in) :: n
    real(dp), intent(in), dimension(3) :: dshift

    select case (n)
    case (1)
      this%vf(4)%fc(:,2)=this%vf(4)%fc(:,2)+dshift
      this%vf(1)%fc(:,1)=this%vf(1)%fc(:,1)+dshift
    case (2)             
      this%vf(1)%fc(:,2)=this%vf(1)%fc(:,2)+dshift
      this%vf(2)%fc(:,1)=this%vf(2)%fc(:,1)+dshift
    case (3)            
      this%vf(2)%fc(:,2)=this%vf(2)%fc(:,2)+dshift
      this%vf(3)%fc(:,1)=this%vf(3)%fc(:,1)+dshift
    case (4)           
      this%vf(3)%fc(:,2)=this%vf(3)%fc(:,2)+dshift
      this%vf(4)%fc(:,1)=this%vf(4)%fc(:,1)+dshift
    case default
      error stop 'n may only take values 1,2,3 or 4'
    end select

  end subroutine vrclass_shiftdP

  subroutine vrclass_rot(this,Tmat)
    ! Rotate vortex ring using Tmat
  class(vr_class) :: this
    real(dp), intent(in), dimension(3,3) :: Tmat
    integer :: i

    do i=1,4
      this%vf(i)%fc(:,1)=matmul(Tmat,this%vf(i)%fc(:,1))
      this%vf(i)%fc(:,2)=matmul(Tmat,this%vf(i)%fc(:,2))
    enddo

  end subroutine vrclass_rot

  subroutine vrclass_calclength(this,isoriginal)
  class(vr_class) :: this
    logical, intent(in) :: isoriginal
    integer :: i
    do i=1,4
      call this%vf(i)%calclength(isoriginal)
    enddo
  end subroutine vrclass_calclength

  subroutine vrclass_strain(this)
  class(vr_class) :: this
    integer :: i
    do i=1,4
      call this%vf(i)%strain()
    enddo
  end subroutine vrclass_strain

end module vr_classdef


!------+-------------------+------|
! ++++ | MODULE DEFINITION | ++++ |
!------+-------------------+------|
module wingpanel_classdef
  use vr_classdef
  implicit none
  type wingpanel_class
    type(vr_class) :: vr
    real(dp), dimension(3,4) :: pc    ! panel coords
    real(dp), dimension(3) :: cp      ! coll point coords
    real(dp), dimension(3) :: ncap    ! unit normal vector
    real(dp), dimension(3) :: velCP   ! local velocity at CP
    real(dp), dimension(3) :: velCPm  ! rel. inertial velocity at CP (due to motion)
    !real(dp), dimension(3) :: dForce  ! panel Force vector in inertial frame
    real(dp) :: vel_pitch             ! pitch velocity
    real(dp) :: dLift, dDrag          ! magnitudes of panel lift and drag
    real(dp) :: delP                  ! Pressure difference at panel
    real(dp) :: panel_area            ! Panel area for computing lift
    real(dp) :: r_hinge               ! dist to point about which pitching occurs (LE of wing)
    real(dp) :: alpha                 ! local angle of attack
    integer :: tag                    ! for identifying panel to be wing or wake
  contains
    procedure :: assignP => wingpanel_class_assignP
    procedure :: calcCP => wingpanel_class_calcCP
    procedure :: calcN => wingpanel_class_calcN
    procedure :: rot => wingpanel_class_rot
    procedure :: shiftdP => wingpanel_class_shiftdP
    procedure :: calc_alpha
    procedure :: calc_area
    procedure :: orthproj
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

  subroutine wingpanel_class_assignP(this,n,P)   ! for assigning coordinates to nth corner
  class(wingpanel_class) :: this
    integer, intent(in) :: n
    real(dp), dimension(3) :: P

    if (n>0 .and. n<5) then
      this%pc(:,n)=P
    else 
      error stop 'n may only take values 1,2,3 or 4'
    endif

  end subroutine wingpanel_class_assignP

  subroutine wingpanel_class_calcCP(this)
  class(wingpanel_class) :: this

    this%cp(1)=this%pc(1,1)+(this%pc(1,2)-this%pc(1,1))*0.75_dp
    this%cp(2)=this%pc(2,1)+(this%pc(2,4)-this%pc(2,1))*0.50_dp
    this%cp(3)=0._dp
  end subroutine wingpanel_class_calcCP

  subroutine wingpanel_class_calcN(this)
  class(wingpanel_class) :: this
    this%ncap=cross3(this%pc(:,3)-this%pc(:,1),this%pc(:,4)-this%pc(:,2))
    this%ncap=this%ncap/norm2(this%ncap)
  end subroutine wingpanel_class_calcN

  subroutine wingpanel_class_rot(this,Tmat)
  class(wingpanel_class) :: this
    real(dp), dimension(3,3) :: Tmat
    integer :: i

    do i=1,4
      this%pc(:,i)=matmul(Tmat,this%pc(:,i))
    enddo
    call this%vr%rot(Tmat)
    this%cp=matmul(Tmat,this%cp)
    this%ncap=matmul(Tmat,this%ncap)
  end subroutine wingpanel_class_rot

  subroutine wingpanel_class_shiftdP(this,dshift)
  class(wingpanel_class) :: this
    real(dp), intent(in), dimension(3) :: dshift
    integer :: i

    this%cp=this%cp+dshift
    do i=1,4
      this%pc(:,i)=this%pc(:,i)+dshift
      call this%vr%shiftdP(i,dshift)
    enddo

  end subroutine wingpanel_class_shiftdP

  subroutine calc_area(this)
  class(wingpanel_class) :: this
    this%panel_area=0.5_dp*norm2(cross3(this%pc(:,3)-this%pc(:,1),this%pc(:,4)-this%pc(:,2)))
  end subroutine calc_area

  subroutine calc_alpha(this)
  class(wingpanel_class) :: this
    real(dp), dimension(3) :: tau_c
    tau_c=this%pc(:,2)-this%pc(:,1)
    tau_c=tau_c/norm2(tau_c)
    this%alpha=0.5_dp*pi
    if (dot_product(this%velCPm,tau_c)>eps) then
      this%alpha=atan((dot_product(this%velCPm,this%ncap)+this%vel_pitch)/dot_product(this%velCPm,tau_c))
    endif
  end subroutine calc_alpha

  ! Calculates the orthogonal projection operator
  function orthproj(this)
  class(wingpanel_class) :: this
    real(dp), dimension(3,3) :: orthproj
    real(dp), dimension(3,3) :: idenmat
    real(dp), dimension(3) :: velCPm_cap
    idenmat(:,1)=(/1._dp,0._dp,0._dp/)
    idenmat(:,2)=(/0._dp,1._dp,0._dp/)
    idenmat(:,3)=(/0._dp,0._dp,1._dp/)
    velCPm_cap=this%velCPm/norm2(this%velCPm)
    orthproj=idenmat-outer_product(velCPm_cap,velCPm_cap)
  end function orthproj

end module wingpanel_classdef


!------+-------------------+------|
! ++++ | MODULE DEFINITION | ++++ |
!------+-------------------+------|
module wakepanel_classdef
  use vr_classdef
  implicit none
  type wakepanel_class
    type(vr_class) :: vr
    integer :: tag                   ! for identifying panel to be wing or wake
  end type wakepanel_class

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

end module wakepanel_classdef


!------+-------------------+------|
! ++++ | MODULE DEFINITION | ++++ |
!------+-------------------+------|
module blade_classdef
  use wingpanel_classdef
  use wakepanel_classdef
  implicit none
  type blade_class
    type(wingpanel_class), allocatable, dimension(:,:) :: wiP
    type(wakepanel_class), allocatable, dimension(:,:) :: waP
    real(dp) :: theta
    real(dp) :: psi

  end type blade_class

end module blade_classdef


!------+-------------------+------|
! ++++ | MODULE DEFINITION | ++++ |
!------+-------------------+------|
module rotor_classdef
  use blade_classdef
  use library
  implicit none
  type rotor_class
    integer :: nb,ns,nc
    type(blade_class), allocatable, dimension(:) :: blade
    real(dp), dimension :: Omega
    real(dp), dimension(3) :: shaft_axis
    real(dp), dimension(3) :: hub_coords, CG_coords
    real(dp) :: radius, chord, root_cut
    real(dp), dimension(4) :: control_pitch  ! theta0,thetaC,thetaS
    real(dp) :: theta_twist
    real(dp) :: pivotLE  ! pivot location from LE [x/c]
    real(dp) :: flap_hinge  ! hinge location from centre [x/R]
    real(dp), dimension(3) :: uvw_body, pqr_body
    real(dp) :: spanwise_core, chordwise_core
  contains
    procedure :: getdata
  end type rotor_class

contains

  subroutine getdata(this,filename,nt)
  class(rotor_class) :: this
    character(len=*), intent(in) :: filename
    integer, intent(in) :: nt  ! nt passed for allocting wake panels
    integer :: i
    real(dp) allocatable, dimension(:) :: chord_sectional

    open(unit=12,file=filename)
    call skiplines(12,2)
    read(12,*) this%nb
    call skiplines(12,3)
    read(12,*) this%ns,this%nc
    call skiplines(12,4)
    read(12,*) this%hub_coords(1),this%hub_coords(2),this%hub_coords(3)
    call skiplines(12,3)
    read(12,*) this%CG_coords(1),this%CG_coords(2),this%CG_coords(3)
    call skiplines(12,4)
    read(12,*) this%radius, this%root_cut, this%chord
    call skiplines(12,4)
    read(12,*) this%Omega, this%shaft_axis(1), this%shaft_axis(2), this%shaft_axis(3)
    call skiplines(12,3)
    read(12,*) this%control_pitch(1), this%control_pitch(2),this%control_pitch(3), this%theta_twist
    call skiplines(12,4)
    read(12,*) this%uvw_body(1), this%uvw_body(2), this%uvw_body(3) &
      ,        this%pqr_body(1), this%pqr_body(2), this%pqr_body(3)
    call skiplines(12,4)
    read(12,*) this%pivotLE, this%flap_hinge
    call skiplines(12,4)
    read(12,*) this%spanwise_core, this%chordwise_core
    call skiplines(12,4)
    read(12,*) this%init_wake_vel, this%psi_start
    close(12)

    ! Conversions
    do i=1,3
      call degtorad(this%control_pitch(i))
    enddo
    call degtorad(this%theta_twist)
    call degtorad(this%psi_start)
    this%spanwise_core=this%spanwise_core*chord
    this%chordwise_core=this%chordwise_core*chord

    ! Allocate rotor object variables
    allocate(this%blade(nb))
    ! Allocate blade object variables
    do iblade=1,this%nb
      allocate(this%blade(iblade)%wiP(this%nc,this%ns))
      allocate(this%blade(iblade)%waP(nt,this%ns))
    enddo

  end subroutine getdata

end module rotor_classdef
