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
    real(dp), dimension(3) :: dForce  ! panel Force vector in inertial frame
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
    tau_c=this%pc(:,4)-this%pc(:,1)
    tau_c=tau_c/norm2(tau_c)
    this%alpha=atan((dot_product(this%velCPm,this%ncap)+this%vel_pitch)/dot_product(this%velCPm,tau_c))
  end subroutine calc_alpha
end module wingpanel_classdef
