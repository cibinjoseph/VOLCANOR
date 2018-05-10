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
    procedure :: isCPinsidecore
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

  ! Checks whether CP lies inside viscous core region of vortex ring
  function isCPinsidecore(this)
  class(wingpanel_class), intent(in) :: this
    logical :: isCPinsidecore
    real(dp) :: deltaxby4, deltayby2

    deltaxby4=0.25_dp*abs(this%vr%vf(1)%fc(1,1)-this%vr%vf(2)%fc(1,1))
    deltayby2=0.5_dp *abs(this%vr%vf(1)%fc(2,1)-this%vr%vf(4)%fc(2,1))

    isCPinsidecore = .false.
    if (deltayby2 .lt. this%vr%vf(1)%r_vc) then
      isCPinsidecore = .true.    ! Left edge
    elseif (deltayby2 .lt. this%vr%vf(3)%r_vc) then
      isCPinsidecore = .true.  ! Right edge
    elseif (deltaxby4 .lt. this%vr%vf(2)%r_vc) then
      isCPinsidecore = .true.  ! Upper edge
    elseif (3._dp*deltaxby4 .lt. this%vr%vf(4)%r_vc) then
      isCPinsidecore = .true.  ! Bottom edge
    endif
  end function isCPinsidecore

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
  use library
  implicit none
  type blade_class
    type(wingpanel_class), allocatable, dimension(:,:) :: wiP
    type(wakepanel_class), allocatable, dimension(:,:) :: waP
    real(dp) :: theta
    real(dp) :: psi
    real(dp) :: pivotLE
  contains
    procedure :: blade_move => move
    procedure :: blade_rot_pitch => rot_pitch
    procedure :: blade_rot_axis => rot_axis
    procedure :: blade_rot_pts => rot_pts
  end type blade_class
contains

  subroutine blade_move(this,dshift)
  class(blade_class) :: this
    real(dp), intent(in), dimension(3) :: dshift
    integer :: i,j

    do j=1,size(this%wiP,2)
      do i=1,size(this%wiP,1)
        call this%wiP(i,j)%shiftdP(dshift)
      enddo
    enddo
  end subroutine blade_move

  subroutine blade_rot_pts(this,pts,origin,order)
  class(blade_class), intent(inout) :: this
    real(dp), dimension(3), intent(in) :: pts    ! pts => phi,theta,psi
    real(dp), dimension(3), intent(in) :: origin ! rotation about
    integer, intent(in) :: order    ! [1]gb & +ve theta , [2]bg & -ve theta
    integer :: i, j
    real(dp), dimension(3,3) :: TMat

    select case (order)
    case (2)
      TMat=Tbg((/cos(pts(1)),sin(pts(1))/),&
        (/cos(pts(2)),sin(pts(2))/),&
        (/cos(pts(3)),sin(pts(3))/))
    case (1)
      TMat=Tgb((/cos(pts(1)),sin(pts(1))/),&
        (/cos(pts(2)),sin(pts(2))/),&
        (/cos(pts(3)),sin(pts(3))/))
    case default
      error stop 'Error: wrong option for order'
    end select

    do j=1,size(this%wiP,2)
      do i=1,size(this%wiP,1)
        call this%wiP(i,j)%shiftdP(-origin)
        call this%wiP(i,j)%rot(TMat)
        call this%wiP(i,j)%shiftdP(origin)
      enddo
    enddo

  end subroutine blade_rot_pts

  subroutine blade_rot_pitch(this,theta)  !pitch about pivotLE from LE
    ! pivot point calculated using straight line joining LE and TE of root panels
  class(blade_class), intent(inout) :: this
    real(dp), intent(in) :: theta
    real(dp), dimension(3) :: axis  
    real(dp), dimension(3) :: origin

    if (abs(theta)>eps) then
      origin=this%wiP(1,1)%pc(:,1)*(1._dp-this%pivotLE)+this%wiP(rows,1)%pc(:,2)*this%pivotLE

      ! Construct axes of rotation from LE of first panel
      axis=this%wiP(1,1)%pc(:,4)-this%wiP(1,1)%pc(:,1)

      call this%rot_axis(theta,axis,origin)
    endif
  end subroutine blade_rot_pitch

  subroutine blade_rot_axis(this,theta,axis,origin)  !rotate about axis at origin
  class(blade_class), intent(inout) :: this
    real(dp), intent(in) :: theta
    real(dp), intent(in), dimension(3) :: axis  
    real(dp), intent(inout) :: origin
    real(dp), dimension(3,3) :: Tmat
    integer :: i,j,rows
    real(dp), dimension(3) :: dshift
    real(dp) :: ct,st,omct

    if (abs(theta)>eps) then
      ! Translate to origin
      rows=size(this%wiP,1)
      call this%move(-origin)

      ! Calculate TMat
      ct=cos(theta)
      st=sin(theta)
      omct=1-ct

      Tmat(:,1)=(/      ct+axis(1)*axis(1)*omct,  axis(3)*st+axis(2)*axis(1)*omct, -axis(2)*st+axis(3)*axis(1)*omct/)
      Tmat(:,2)=(/-axis(3)*st+axis(1)*axis(2)*omct,       ct+axis(2)*axis(2)*omct,  axis(1)*st+axis(3)*axis(2)*omct/)
      Tmat(:,3)=(/ axis(2)*st+axis(1)*axis(3)*omct, -axis(1)*st+axis(2)*axis(3)*omct,       ct+axis(3)*axis(3)*omct/)

      ! Rotate about axis
      do j=1,size(this%wiP,2)
        do i=1,size(this%wiP,1)
          call this%wiP(i,j)%rot(TMat)
        enddo
      enddo

      ! Untranslate from origin
      call this%move(origin)
    endif
  end subroutine blade_rot_axis

end module blade_classdef


!------+-------------------+------|
! ++++ | MODULE DEFINITION | ++++ |
!------+-------------------+------|
module rotor_classdef
  use blade_classdef
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
    real(dp) :: spanwise_core, streamwise_core
    real(dp), allocatable, dimension(:,:) :: AIC,AIC_inv  ! Influence coefficient matrix
  contains
    procedure :: rotor_getdata => getdata
    procedure :: rotor_init => init
    procedure :: rotor_move => move
    procedure :: rotor_rot_pts => rot_pts
    procedure :: rotor_pitch => pitch
    procedure :: assignshed
    procedure :: convectwake
  end type rotor_class

contains

  !-----+--------------------------+-----|
  ! -+- | Initialization Functions | -+- |
  !-----+--------------------------+-----|

  subroutine rotor_getdata(this,filename,nt)
  class(rotor_class) :: this
    character(len=*), intent(in) :: filename
    integer, intent(in) :: nt  ! nt passed for allocting wake panels
    integer :: i

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
    read(12,*) this%spanwise_core, this%streamwise_core
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
    this%streamwise_core=this%streamwise_core*chord

    ! Allocate rotor object variables
    allocate(this%blade(nb))
    ! Allocate blade object variables
    do iblade=1,this%nb
      allocate(this%blade(iblade)%wiP(this%nc,this%ns))
      allocate(this%blade(iblade)%waP(nt,this%ns))
      allocate(this%blade(iblade)%AIC(this%nc*this%ns,this%nc*this%ns))
      allocate(this%blade(iblade)%AIC_inv(this%nc*this%ns,this%nc*this%ns))
    enddo

  end subroutine rotor_getdata

  subroutine rotor_init(this,span_spacing_switch,nt,dt)
  class(rotor_class) :: this
    integer, intent(in) :: span_spacing_switch
    integer, intent(in) :: nt

    real(dp), dimension(this%nc+1) :: xvec
    real(dp), dimension(this%ns+1) :: yvec
    integer :: i,j,iblade
    real(dp) :: blade_offset

    ! Blade initialization
    xvec=linspace(-this%chord,0._dp,this%nc+1)
    select case (span_spacing_switch)
    case (1)
      yvec=linspace(this%root_cut*this%span,this%span,this%ns+1)
    case (2)
      yvec=cosspace(this%root_cut*this%span,this%span,this%ns+1)
    case (3)
      yvec=halfsinspace(this%root_cut*this%span,this%span,this%ns+1)
    end select

    do iblade=1,this%nb
      ! Initialize panel coordinates
      do j=1,this%ns
        do i=1,this%nc
          call this%blade(iblade)%wiP(i,j)%assignP(1,(/xvec(i  ),yvec(j  ),0._dp/))
          call this%blade(iblade)%wiP(i,j)%assignP(2,(/xvec(i+1),yvec(j  ),0._dp/))
          call this%blade(iblade)%wiP(i,j)%assignP(3,(/xvec(i+1),yvec(j+1),0._dp/))
          call this%blade(iblade)%wiP(i,j)%assignP(4,(/xvec(i  ),yvec(j+1),0._dp/))
        enddo
      enddo

      ! Initialize vr coords of all panels except last row (to accomodate mismatch of vr coords when usi    ng unequal spacing)
      do i=1,this%nc-1
        xshiftLE=(xvec(i+1)-xvec(i))*0.25_dp  ! Shift x coord by dx/4
        xshiftTE=(xvec(i+2)-xvec(i+1))*0.25_dp  ! Shift x coord by dx/4
        do j=1,this%ns
          call this%blade(iblade)%wiP(i,j)%vr%assignP(1,(/xvec(i  )+xshiftLE,yvec(j  ),0._dp/))
          call this%blade(iblade)%wiP(i,j)%vr%assignP(2,(/xvec(i+1)+xshiftTE,yvec(j  ),0._dp/))
          call this%blade(iblade)%wiP(i,j)%vr%assignP(3,(/xvec(i+1)+xshiftTE,yvec(j+1),0._dp/))
          call this%blade(iblade)%wiP(i,j)%vr%assignP(4,(/xvec(i  )+xshiftLE,yvec(j+1),0._dp/))
        enddo
      enddo

      ! Initializing vr coords of last row
      xshiftLE=(xvec(this%nc+1)-xvec(this%nc))*0.25_dp  ! Shift x coord by dx/4
      xshiftTE=0._dp
      do j=1,this%ns
        call this%blade(iblade)%wiP(this%nc,j)%vr%assignP(1,(/xvec(this%nc  )+xshiftLE,yvec(j  ),0._dp/))
        call this%blade(iblade)%wiP(this%nc,j)%vr%assignP(2,(/xvec(this%nc+1)+xshiftTE,yvec(j  ),0._dp/))
        call this%blade(iblade)%wiP(this%nc,j)%vr%assignP(3,(/xvec(this%nc+1)+xshiftTE,yvec(j+1),0._dp/))
        call this%blade(iblade)%wiP(this%nc,j)%vr%assignP(4,(/xvec(this%nc  )+xshiftLE,yvec(j+1),0._dp/))
      enddo

      ! Shed them!!
      v_shed=0.2_dp*this%radius*this%Omega
      do j=1,this%ns
        call this%blade(iblade)%wiP(this%nc,j)%vr%shiftdP(2,v_shed*dt)
        call this%blade(iblade)%wiP(this%nc,j)%vr%shiftdP(3,v_shed*dt)
      enddo

      ! Initialize CP coords, ncap, panel_area and pivotLE
      do j=1,this%ns
        do i=1,this%nc
          call this%blade(iblade)%wiP(i,j)%calcCP()
          call this%blade(iblade)%wiP(i,j)%calcN()
          this%blade(iblade)%wiP(i,j)%r_hinge=length3d((this%blade(iblade)%wiP(1,j)%pc(:,1)+this%blade(iblade)%wiP(1,j)%pc(:,4))*0.5_dp,this%blade(iblade)%(i,j)%cp)
          call this%blade(iblade)%wiP(i,j)%calc_area()
        enddo
      enddo

      ! Initialize gamma
      this%blade(iblade)%wiP%vr%gam=0._dp
      this%blade(iblade)%pivotLE=this%pivotLE

      ! Initialize mid vortex core radius
      do i=1,4
        this%blade(iblade)%wiP%vr%vf(i)%r_vc0= this%spanwise_core
        this%blade(iblade)%wiP%vr%vf(i)%r_vc = this%spanwise_core
        this%blade(iblade)%wiP%vr%vf(i)%age  = 0._dp
      enddo

      ! Initialize tip vortex core radius
      do i=1,this%nc
        this%blade(iblade)%wiP(i,1)%vr%vf(1)%r_vc0       = this%streamwise_core
        this%blade(iblade)%wiP(i,1)%vr%vf(1)%r_vc        = this%streamwise_core
        this%blade(iblade)%wiP(i,this%ns)%vr%vf(3)%r_vc0 = this%streamwise_core
        this%blade(iblade)%wiP(i,this%ns)%vr%vf(3)%r_vc  = this%streamwise_core
      enddo

      ! Verify CP is outside vortex core for boundary panels
      if (isCPinsidecore(this%blade(iblade)%wiP(1,1))) then
        print*,'Warning: CP inside vortex core at panel LU'
        print*,'Any key to continue. Ctrl-C to exit'
        read(*,*)
      endif
      if (isCPinsidecore(this%blade(iblade)%wiP(this%nc,1))) then
        print*,'Warning: CP inside vortex core at panel LB'
        print*,'Any key to continue. Ctrl-C to exit'
        read(*,*)
      endif
      if (isCPinsidecore(this%blade(iblade)%wiP(1,this%ns))) then
        print*,'Warning: CP inside vortex core at panel RU'
        print*,'Any key to continue. Ctrl-C to exit'
        read(*,*)
      endif
      if (isCPinsidecore(this%blade(iblade)%wiP(this%nc,this%ns))) then
        print*,'Warning: CP inside vortex core at panel RB'
        print*,'Any key to continue. Ctrl-C to exit'
        read(*,*)
      endif
    enddo

    ! Move rotor to hub coordinates
    do iblade=1,this%nb
      call this%blade(iblade)%move(this%hub_coords)
    enddo

    ! Rotate remaining blades to their positions
    ! Rotate blades for multi-bladed rotors
    do iblade=2,this%nb
      blade_offset=2._dp*pi/this%nb*(iblade-1)
      call this%blade(iblade)%rot_axis(blade_offset,this%shaft_axis,this%hub_coords,1)
    enddo


    ! Wake initialization
    ! Assign core_radius to mid vortices
    do iblade=1,this%nb
      do i=1,4
        this%blade(iblade)%waP%vr%vf(i)%r_vc0 = this%spanwise_core
        this%blade(iblade)%waP%vr%vf(i)%r_vc  = this%spanwise_core
        this%blade(iblade)%waP%vr%vf(i)%age=0._dp
      enddo

      this%blade(iblade)%waP%tag=-1
      this%blade(iblade)%waP%vr%gam=0._dp

      ! Assign core_radius to tip vortices
      do i=1,nt
        ! Root vortex
        this%blade(iblade)%waP(i,1)%vr%vf(1)%r_vc0      = this%streamwise_core
        this%blade(iblade)%waP(i,1)%vr%vf(1)%r_vc       = this%streamwise_core
        this%blade(iblade)%waP(i,1)%vr%vf(3)%r_vc0      = this%streamwise_core
        this%blade(iblade)%waP(i,1)%vr%vf(3)%r_vc       = this%streamwise_core

        this%blade(iblade)%waP(i,2)%vr%vf(1)%r_vc0      = this%streamwise_core
        this%blade(iblade)%waP(i,2)%vr%vf(1)%r_vc       = this%streamwise_core
        this%blade(iblade)%waP(i,2)%vr%vf(3)%r_vc0      = this%streamwise_core
        this%blade(iblade)%waP(i,2)%vr%vf(3)%r_vc       = this%streamwise_core

        ! Tip vortex
        this%blade(iblade)%waP(i,this%ns)%vr%vf(1)%r_vc0   = this%streamwise_core
        this%blade(iblade)%waP(i,this%ns)%vr%vf(1)%r_vc    = this%streamwise_core
        this%blade(iblade)%waP(i,this%ns)%vr%vf(3)%r_vc0   = this%streamwise_core
        this%blade(iblade)%waP(i,this%ns)%vr%vf(3)%r_vc    = this%streamwise_core

        this%blade(iblade)%waP(i,this%ns-1)%vr%vf(3)%r_vc0 = this%streamwise_core
        this%blade(iblade)%waP(i,this%ns-1)%vr%vf(3)%r_vc  = this%streamwise_core
        !this%blade(iblade)%waP(i,this%ns-1)%vr%vf(3)%r_vc0 = this%streamwise_core
        !this%blade(iblade)%waP(i,this%ns-1)%vr%vf(3)%r_vc  = this%streamwise_core
      enddo

      !if (starting_vortex_core > eps) then
      !  ! Assign core_radius to starting vortices
      !  do i=1,this%ns
      !    do j=2,4,2
      !      this%blade(iblade)%waP(rows,i)%vr%vf(j)%r_vc0 = starting_vortex_core
      !      this%blade(iblade)%waP(rows,i)%vr%vf(j)%r_vc  = starting_vortex_core
      !      this%blade(iblade)%waP(rows-1,i)%vr%vf(j)%r_vc0 = starting_vortex_core
      !      this%blade(iblade)%waP(rows-1,i)%vr%vf(j)%r_vc  = starting_vortex_core
      !    enddo
      !  enddo
      !endif
    enddo

  end subroutine rotor_init

  !-----+------------------+-----|
  ! -+- | Motion Functions | -+- |
  !-----+------------------+-----|

  subroutine rotor_move(this,dshift)
  class(rotor_class) :: this
    real(dp), intent(in), dimension(3) :: dshift

    integer :: iblade

    do iblade=1,this%nb
      call this%blade(iblade)%move(dshift)
    enddo
    this%hub_coords=this%hub_coords+dshift
    this%CG_coords=this%CG_coords+dshift

  end subroutine rotor_move

  subroutine rotor_rot_pts(this,pts,origin,order)
  class(rotor_class), intent(inout) :: this
    real(dp), dimension(3), intent(in) :: pts    ! pts => phi,theta,psi
    real(dp), dimension(3), intent(in) :: origin ! rotation about
    integer, intent(in) :: order    ! [1]gb & +ve theta , [2]bg & -ve theta
    integer :: iblade
    real(dp), dimension(3,3) :: TMat

    select case (order)
    case (2)
      TMat=Tbg((/cos(pts(1)),sin(pts(1))/),&
        (/cos(pts(2)),sin(pts(2))/),&
        (/cos(pts(3)),sin(pts(3))/))
    case (1)
      TMat=Tgb((/cos(pts(1)),sin(pts(1))/),&
        (/cos(pts(2)),sin(pts(2))/),&
        (/cos(pts(3)),sin(pts(3))/))
    case default
      error stop 'Error: wrong option for order'
    end select

    do iblade=1,this%nb
      call this%blade(iblade)%rot_pts(pts,origin,order)
    enddo

    this%shaft_axis=matmul(TMat,this%shaft_axis)

    this%hub_coords=this%hub_coords-origin
    this%hub_coords=matmul(TMat,this%hub_coords)
    this%hub_coords=this%hub_coords+origin

    this%CG_coords=this%CG_coords-origin
    this%CG_coords=matmul(TMat,this%CG_coords)
    this%CG_coords=this%CG_coords+origin

  end subroutine rotor_rot_pts

  subroutine rotor_pitch(this,theta_pitch)
  class(rotor_class), intent(inout) :: this
    real(dp), intent(in) :: theta_pitch
    integer :: iblade

    do iblade=1,nb
      call this%blade(iblade)%rot_pitch(theta_pitch)
    enddo
  end subroutine rotor_pitch

  !-----+----------------+-----|
  ! -+- | Wake Functions | -+- |
  !-----+----------------+-----|

  ! Assigns coordinates to first row of wake from last row of blade
  subroutine assignshed(this,row_now,edge)
  class(rotor_class), intent(inout) :: this
    integer, intent(in) :: row_now
    character(len=2), intent(in) :: edge
    integer :: i, iblade

    do iblade=1,this%nb
      this%blade(iblade)%waP%vr%gam=this%blade(iblade)%wiP%vr%gam
    enddo

    select case (edge)
    case ('LE')    ! assign to LE
      do iblade=1,this%nb
        do i=1,this%ns
          call this%blade(iblade)%waP(row_now,i)%vr%assignP(1,this%blade(iblade)%wiP(this%nc,i)%vr%vf(2)%fc(:,1))
          call this%blade(iblade)%waP(row_now,i)%vr%assignP(4,this%blade(iblade)%wiP(this%nc,i)%vr%vf(3)%fc(:,1))
          call this%blade(iblade)%waP(row_now,i)%vr%calclength(.TRUE.)    ! TRUE => record original length
        enddo
        wake_row%tag=1
      enddo
    case ('TE')    ! assign to TE
      do iblade=1,this%nb
        do i=1,this%ns
          call this%blade(iblade)%waP(row_now,i)%vr%assignP(2,this%blade(iblade)%wiP(this%nc,i)%vr%vf(2)%fc(:,1))
          call this%blade(iblade)%waP(row_now,i)%vr%assignP(3,this%blade(iblade)%wiP(this%nc,i)%vr%vf(3)%fc(:,1))
        enddo
      enddo
    case default
      error stop 'Error: Wrong option for edge'
    end select

  end subroutine assignshed

  ! Convect wake using dP_array=vind_array*dt
  subroutine convectwake(this,row_now,dP_array)
  class(rotor_class), intent(inout) :: this
    integer, intent(in) :: row_now
    real(dp), intent(in), dimension(:,:,:) :: dP_array
    integer :: i,j,rows,cols,iblade

    rows=size(this%blade(iblade)%waP,1)-row_now+1
    cols=this%ns

    do iblade=1,this%nb
      !$omp parallel do collapse(2)
      do j=1,cols
        do i=row_now,rows
          call this%blade(iblade)%waP(i,j)%vr%shiftdP(2,dP_array(:,i,j))
        enddo
      enddo
      !$omp end parallel do

      !$omp parallel do
      do i=row_now,rows
        call this%blade(iblade)%waP(i,cols)%vr%shiftdP(3,dP_array(:,i,cols+1))
      enddo
      !$omp end parallel do
    enddo

    call wake_continuity

  end subroutine convectwake

  subroutine calcAIC(this)
  class(rotor_class), intent(inout) :: this
    integer :: iblade,jblade,ispan,ichord,i,j,row,col
    real(dp), dimension(3) :: vec_dummy

    ! Influence Coefficient Matrix
    do iblade=1,this%nb
      do ispan=1,this%ns      ! Collocation point loop
        do ichord=1,this%nc
          row=ichord+this%nc*(ispan-1)+this%ns*this%nc*(iblade-1)

          do jblade=1,this%nb
            do j=1,this%ns       ! Vortex ring loop
              do i=1,this%nc
                col=i+nc*(j-1)+this%ns*this%nc*(jblade-1)
                vec_dummy=this%blade(jblade)%wiP(i,j)%vr%vind(thisi%blade(iblade)%wiP(ichord,ispan)%CP)
                this%blade(iblade)%AIC(row,col)=dot_product(vec_dummy,this%blade(iblade)%wiP(ichord,ispan)%ncap)
              enddo
            enddo
          enddo

        enddo
      enddo
      this%blade(iblade)%AIC_inv=inv(this%blade(iblade)%AIC)
    enddo
  end subroutine caclAIC

end module rotor_classdef
