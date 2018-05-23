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
    real(dp) :: pivotLE
    real(dp), allocatable, dimension(:,:,:) :: vind_wake
    real(dp), allocatable, dimension(:,:,:) :: vind_wake1, vind_wake2, vind_wake3
    real(dp), allocatable, dimension(:,:,:) :: Pvind_wake, vind_wake_step
    real(dp), allocatable, dimension(:,:) :: Pwake

  contains
    procedure :: move => blade_move
    procedure :: rot_pitch 
    procedure :: rot_axis
    procedure :: rot_pts => blade_rot_pts
    procedure :: vind_bywing => blade_vind_bywing
    procedure :: vind_bywake => blade_vind_bywake
    procedure :: convectwake
    procedure :: wake_continuity
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

  subroutine rot_pitch(this,theta)  !pitch about pivotLE from LE
    ! pivot point calculated using straight line joining LE and TE of root panels
  class(blade_class), intent(inout) :: this
    real(dp), intent(in) :: theta
    real(dp), dimension(3) :: axis  
    real(dp), dimension(3) :: origin
    integer :: rows

    if (abs(theta)>eps) then
      rows=size(this%wiP,1)
      origin=this%wiP(1,1)%pc(:,1)*(1._dp-this%pivotLE)+this%wiP(rows,1)%pc(:,2)*this%pivotLE

      ! Construct axes of rotation from LE of first panel
      axis=this%wiP(1,1)%pc(:,4)-this%wiP(1,1)%pc(:,1)

      call this%rot_axis(theta,axis,origin)
    endif
  end subroutine rot_pitch

  subroutine rot_axis(this,theta,axis,origin)  !rotate about axis at origin
  class(blade_class), intent(inout) :: this
    real(dp), intent(in) :: theta
    real(dp), intent(inout), dimension(3) :: axis  
    real(dp), intent(in), dimension(3) :: origin
    real(dp), dimension(3,3) :: Tmat
    integer :: i,j,rows
    real(dp) :: ct,st,omct

    if (abs(theta)>eps) then
      ! Translate to origin
      rows=size(this%wiP,1)
      call this%move(-origin)

      ! Ensure axis is normalized
      axis=axis/norm2(axis)

      ! Calculate TMat
      ct=cos(theta)
      st=sin(theta)
      omct=1-ct

      Tmat(:,1)=(/         ct+axis(1)*axis(1)*omct,  axis(3)*st+axis(2)*axis(1)*omct, -axis(2)*st+axis(3)*axis(1)*omct/)
      Tmat(:,2)=(/-axis(3)*st+axis(1)*axis(2)*omct,          ct+axis(2)*axis(2)*omct,  axis(1)*st+axis(3)*axis(2)*omct/)
      Tmat(:,3)=(/ axis(2)*st+axis(1)*axis(3)*omct, -axis(1)*st+axis(2)*axis(3)*omct,          ct+axis(3)*axis(3)*omct/)

      ! Rotate about axis
      do j=1,size(this%wiP,2)
        do i=1,size(this%wiP,1)
          call this%wiP(i,j)%rot(TMat)
        enddo
      enddo

      ! Untranslate from origin
      call this%move(origin)
    endif
  end subroutine rot_axis

  function blade_vind_bywing(this,P)  ! Induced velocity at a point P
    ! pivot point calculated using straight line joining LE and TE of root panels
  class(blade_class), intent(inout) :: this
    real(dp), intent(in), dimension(3) :: P
    real(dp), dimension(3) :: blade_vind_bywing
    integer :: i,j

    blade_vind_bywing=0._dp
    do j=1,size(this%wiP,2)
      do i=1,size(this%wiP,1)
        blade_vind_bywing=blade_vind_bywing+this%wiP(i,j)%vr%vind(P)*this%wiP(i,j)%vr%gam
      enddo
    enddo

  end function blade_vind_bywing

  function blade_vind_bywake(this,row_now,P)  ! Induced velocity at a point P
    ! pivot point calculated using straight line joining LE and TE of root panels
  class(blade_class), intent(inout) :: this
    integer, intent(in) :: row_now
    real(dp), intent(in), dimension(3) :: P
    real(dp), dimension(3) :: blade_vind_bywake
    integer :: i,j

    blade_vind_bywake=0._dp
    do j=1,size(this%waP,2)
      do i=row_now,size(this%waP,1)
        blade_vind_bywake=blade_vind_bywake+this%waP(i,j)%vr%vind(P)*this%waP(i,j)%vr%gam
      enddo
    enddo

  end function blade_vind_bywake

  ! Convect wake using dP_array=vind_array*dt
  subroutine convectwake(this,dP_array)
  class(blade_class), intent(inout) :: this
    real(dp), intent(in), dimension(:,:,:) :: dP_array
    integer :: i,j,rows,cols,nt,index_offset

    rows=size(dP_array,2)
    cols=size(this%waP,2)
    nt=size(this%waP,1)
    index_offset=nt-rows

    !$omp parallel do collapse(2)
    do j=1,cols
      do i=1,rows
        call this%waP(i+index_offset,j)%vr%shiftdP(2,dP_array(:,i,j))
      enddo
    enddo
    !$omp end parallel do

    !$omp parallel do
    do i=1,rows
      call this%waP(i+index_offset,cols)%vr%shiftdP(3,dP_array(:,i,cols+1))
    enddo
    !$omp end parallel do

    call this%wake_continuity(index_offset+1)

  end subroutine convectwake

  ! Maintain continuity between vortex ring elements after convection
  ! of vortex ring corners
  subroutine wake_continuity(this,row_now)
  class(blade_class), intent(inout) :: this
    integer, intent(in) :: row_now
    integer :: i,j,rows,cols

    rows=size(this%waP,1)
    cols=size(this%waP,2)

    !$omp parallel do collapse(2)
    do j=1,cols-1
      do i=row_now+1,rows
        call this%waP(i,j)%vr%assignP(1,this%waP(i-1,j)%vr%vf(2)%fc(:,1))
        call this%waP(i,j)%vr%assignP(3,this%waP(i,j+1)%vr%vf(2)%fc(:,1))
        call this%waP(i,j)%vr%assignP(4,this%waP(i-1,j+1)%vr%vf(2)%fc(:,1))
      enddo
    enddo
    !$omp end parallel do

    !$omp parallel do
    do j=1,cols-1
      call this%waP(row_now,j)%vr%assignP(3,this%waP(row_now,j+1)%vr%vf(2)%fc(:,1))
    enddo
    !$omp end parallel do

    !$omp parallel do
    do i=row_now+1,rows
      call this%waP(i,cols)%vr%assignP(1,this%waP(i-1,cols)%vr%vf(2)%fc(:,1))
      call this%waP(i,cols)%vr%assignP(4,this%waP(i-1,cols)%vr%vf(3)%fc(:,1))
    enddo
    !$omp end parallel do
  end subroutine wake_continuity
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
    real(dp) :: Omega, Omega_slow
    real(dp), dimension(3) :: shaft_axis
    real(dp), dimension(3) :: hub_coords, CG_coords
    real(dp) :: radius, chord, root_cut
    real(dp), dimension(3) :: control_pitch  ! theta0,thetaC,thetaS
    real(dp) :: theta_twist
    real(dp) :: pivotLE  ! pivot location from LE [x/c]
    real(dp) :: flap_hinge  ! hinge location from centre [x/R]
    real(dp), dimension(3) :: v_body, om_body
    real(dp), dimension(3) :: v_wind, om_wind
    real(dp) :: psi
    real(dp), dimension(3) :: pts  ! phi,theta,psi about CG_coords
    real(dp) :: spanwise_core, streamwise_core
    real(dp), allocatable, dimension(:,:) :: AIC,AIC_inv  ! Influence coefficient matrix
    real(dp), allocatable, dimension(:) :: gamvec,gamvec_prev,RHS
    real(dp) :: init_wake_vel, psi_start
  contains
    procedure :: getdata
    procedure :: init
    procedure :: deinit
    procedure :: gettheta
    procedure :: getthetadot
    procedure :: move => rotor_move
    procedure :: rot_pts => rotor_rot_pts
    procedure :: rot_advance => rotor_rot_advance
    procedure :: assignshed
    procedure :: map_gam
    procedure :: age_wake
    procedure :: dissipate_tip
    procedure :: strain_wake
    procedure :: calcAIC
    procedure :: vind_bywing => rotor_vind_bywing
    procedure :: vind_bywake => rotor_vind_bywake
  end type rotor_class

contains

  !-----+--------------------------+-----|
  ! -+- | Initialization Functions | -+- |
  !-----+--------------------------+-----|

  subroutine getdata(this,filename,nt)
  class(rotor_class) :: this
    character(len=*), intent(in) :: filename
    integer, intent(in) :: nt  ! nt passed for allocting wake panels
    integer :: i,ib

    open(unit=12,file=filename)
    call skiplines(12,2)
    read(12,*) this%nb
    call skiplines(12,3)
    read(12,*) this%ns,this%nc
    call skiplines(12,4)
    read(12,*) this%hub_coords(1),this%hub_coords(2),this%hub_coords(3)
    call skiplines(12,3)
    read(12,*) this%CG_coords(1),this%CG_coords(2),this%CG_coords(3)
    call skiplines(12,3)
    read(12,*) this%pts(1),this%pts(2),this%pts(3)
    call skiplines(12,4)
    read(12,*) this%radius, this%root_cut, this%chord
    call skiplines(12,4)
    read(12,*) this%Omega, this%shaft_axis(1), this%shaft_axis(2), this%shaft_axis(3)
    call skiplines(12,3)
    read(12,*) this%control_pitch(1), this%control_pitch(2),this%control_pitch(3), this%theta_twist
    call skiplines(12,4)
    read(12,*) this%v_body(1), this%v_body(2), this%v_body(3) &
      ,        this%om_body(1), this%om_body(2), this%om_body(3)
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
      call degtorad(this%pts(i))
    enddo
    call degtorad(this%theta_twist)
    call degtorad(this%psi_start)
    this%spanwise_core=this%spanwise_core*this%chord
    this%streamwise_core=this%streamwise_core*this%chord

    ! Allocate rotor object variables
    allocate(this%blade(this%nb))
    allocate(this%AIC(this%nc*this%ns*this%nb,this%nc*this%ns*this%nb))
    allocate(this%AIC_inv(this%nc*this%ns*this%nb,this%nc*this%ns*this%nb))
    allocate(this%gamvec(this%nc*this%ns*this%nb))
    allocate(this%gamvec_prev(this%nc*this%ns*this%nb))
    allocate(this%RHS(this%nc*this%ns*this%nb))
    ! Allocate blade object variables
    do ib=1,this%nb
      allocate(this%blade(ib)%wiP(this%nc,this%ns))
      allocate(this%blade(ib)%waP(nt,this%ns))
    enddo

  end subroutine getdata

  subroutine init(this,nt,dt,span_spacing_switch,FDscheme_switch)
  class(rotor_class) :: this
    integer, intent(in) :: nt
    real(dp), intent(in) :: dt
    integer, intent(in) :: span_spacing_switch, FDscheme_switch

    real(dp), dimension(this%nc+1) :: xvec
    real(dp), dimension(this%ns+1) :: yvec
    integer :: i,j,ib
    real(dp) :: blade_offset
    real(dp) :: xshiftLE,xshiftTE,v_shed

    ! Blade initialization
    xvec=linspace(-this%chord,0._dp,this%nc+1)
    select case (span_spacing_switch)
    case (1)
      yvec=linspace(this%root_cut*this%radius,this%radius,this%ns+1)
    case (2)
      yvec=cosspace(this%root_cut*this%radius,this%radius,this%ns+1)
    case (3)
      yvec=halfsinspace(this%root_cut*this%radius,this%radius,this%ns+1)
    end select

    do ib=1,this%nb
      ! Initialize panel coordinates
      do j=1,this%ns
        do i=1,this%nc
          call this%blade(ib)%wiP(i,j)%assignP(1,(/xvec(i  ),yvec(j  ),0._dp/))
          call this%blade(ib)%wiP(i,j)%assignP(2,(/xvec(i+1),yvec(j  ),0._dp/))
          call this%blade(ib)%wiP(i,j)%assignP(3,(/xvec(i+1),yvec(j+1),0._dp/))
          call this%blade(ib)%wiP(i,j)%assignP(4,(/xvec(i  ),yvec(j+1),0._dp/))
        enddo
      enddo

      ! Initialize vr coords of all panels except last row (to accomodate mismatch of vr coords when usi    ng unequal spacing)
      do i=1,this%nc-1
        xshiftLE=(xvec(i+1)-xvec(i))*0.25_dp  ! Shift x coord by dx/4
        xshiftTE=(xvec(i+2)-xvec(i+1))*0.25_dp  ! Shift x coord by dx/4
        do j=1,this%ns
          call this%blade(ib)%wiP(i,j)%vr%assignP(1,(/xvec(i  )+xshiftLE,yvec(j  ),0._dp/))
          call this%blade(ib)%wiP(i,j)%vr%assignP(2,(/xvec(i+1)+xshiftTE,yvec(j  ),0._dp/))
          call this%blade(ib)%wiP(i,j)%vr%assignP(3,(/xvec(i+1)+xshiftTE,yvec(j+1),0._dp/))
          call this%blade(ib)%wiP(i,j)%vr%assignP(4,(/xvec(i  )+xshiftLE,yvec(j+1),0._dp/))
        enddo
      enddo

      ! Initializing vr coords of last row
      xshiftLE=(xvec(this%nc+1)-xvec(this%nc))*0.25_dp  ! Shift x coord by dx/4
      xshiftTE=0._dp
      do j=1,this%ns
        call this%blade(ib)%wiP(this%nc,j)%vr%assignP(1,(/xvec(this%nc  )+xshiftLE,yvec(j  ),0._dp/))
        call this%blade(ib)%wiP(this%nc,j)%vr%assignP(2,(/xvec(this%nc+1)+xshiftTE,yvec(j  ),0._dp/))
        call this%blade(ib)%wiP(this%nc,j)%vr%assignP(3,(/xvec(this%nc+1)+xshiftTE,yvec(j+1),0._dp/))
        call this%blade(ib)%wiP(this%nc,j)%vr%assignP(4,(/xvec(this%nc  )+xshiftLE,yvec(j+1),0._dp/))
      enddo

      ! Shed last row of vortices
      if (abs(norm2(this%v_wind)) < eps) then
        v_shed=0.02_dp*this%chord/(dt*this%nc)
      else
        v_shed=0.2_dp*norm2(this%v_wind)
      endif
      do j=1,this%ns
        call this%blade(ib)%wiP(this%nc,j)%vr%shiftdP(2,(/v_shed*dt,0._dp,0._dp/))
        call this%blade(ib)%wiP(this%nc,j)%vr%shiftdP(3,(/v_shed*dt,0._dp,0._dp/))
      enddo

      ! Initialize CP coords, ncap, panel_area and pivotLE
      do j=1,this%ns
        do i=1,this%nc
          call this%blade(ib)%wiP(i,j)%calcCP()
          call this%blade(ib)%wiP(i,j)%calcN()
          this%blade(ib)%wiP(i,j)%r_hinge=length3d((this%blade(ib)%wiP(1,j)%pc(:,1)  &
            +                                           this%blade(ib)%wiP(1,j)%pc(:,4))*0.5_dp,this%blade(ib)%wiP(i,j)%cp)
          call this%blade(ib)%wiP(i,j)%calc_area()
        enddo
      enddo

      ! Initialize gamma
      this%blade(ib)%wiP%vr%gam=0._dp
      this%blade(ib)%pivotLE=this%pivotLE

      ! Initialize mid vortex core radius
      do i=1,4
        this%blade(ib)%wiP%vr%vf(i)%r_vc0= this%spanwise_core
        this%blade(ib)%wiP%vr%vf(i)%r_vc = this%spanwise_core
        this%blade(ib)%wiP%vr%vf(i)%age  = 0._dp
      enddo

      ! Initialize tip vortex core radius
      do i=1,this%nc
        this%blade(ib)%wiP(i,1)%vr%vf(1)%r_vc0       = this%streamwise_core
        this%blade(ib)%wiP(i,1)%vr%vf(1)%r_vc        = this%streamwise_core
        this%blade(ib)%wiP(i,this%ns)%vr%vf(3)%r_vc0 = this%streamwise_core
        this%blade(ib)%wiP(i,this%ns)%vr%vf(3)%r_vc  = this%streamwise_core
      enddo

      ! Verify CP is outside vortex core for boundary panels
      if (isCPinsidecore(this%blade(ib)%wiP(1,1))) then
        print*,'Warning: CP inside vortex core at panel LU'
        print*,'Any key to continue. Ctrl-C to exit'
        read(*,*)
      endif
      if (isCPinsidecore(this%blade(ib)%wiP(this%nc,1))) then
        print*,'Warning: CP inside vortex core at panel LB'
        print*,'Any key to continue. Ctrl-C to exit'
        read(*,*)
      endif
      if (isCPinsidecore(this%blade(ib)%wiP(1,this%ns))) then
        print*,'Warning: CP inside vortex core at panel RU'
        print*,'Any key to continue. Ctrl-C to exit'
        read(*,*)
      endif
      if (isCPinsidecore(this%blade(ib)%wiP(this%nc,this%ns))) then
        print*,'Warning: CP inside vortex core at panel RB'
        print*,'Any key to continue. Ctrl-C to exit'
        read(*,*)
      endif
    enddo

    ! Move rotor to hub coordinates
    do ib=1,this%nb
      call this%blade(ib)%move(this%hub_coords)
    enddo

    ! Rotate remaining blades to their positions
    ! Rotate blades for multi-bladed rotors
    do ib=2,this%nb
      blade_offset=2._dp*pi/this%nb*(ib-1)
      call this%blade(ib)%rot_axis(blade_offset,this%shaft_axis,this%hub_coords)
    enddo

    ! Assign wind velocities
    this%v_wind=-1._dp*this%v_body
    this%om_wind=-1._dp*this%om_body

    ! Assign pts and dpts
    !this%dpts=this%om_body*dt

    ! Allocate vars required for wake convection
    ! on the basis of finite diff scheme
    do ib=1,this%nb
      select case (FDscheme_switch)
      case (0)
        allocate(this%blade(ib)%vind_wake(3,nt,this%ns+1))
      case (1)
        allocate(this%blade(ib)%Pwake(nt,this%ns))
        allocate(this%blade(ib)%vind_wake(3,nt,this%ns+1))
        allocate(this%blade(ib)%vind_wake1(3,nt,this%ns+1))
        allocate(this%blade(ib)%Pvind_wake(3,nt,this%ns+1))
      case (2)
        allocate(this%blade(ib)%vind_wake(3,nt,this%ns+1))
        allocate(this%blade(ib)%vind_wake1(3,nt,this%ns+1))
        allocate(this%blade(ib)%vind_wake_step(3,nt,this%ns+1))
      case (3)
        allocate(this%blade(ib)%vind_wake(3,nt,this%ns+1))
        allocate(this%blade(ib)%Pwake(nt,this%ns))
        allocate(this%blade(ib)%vind_wake1(3,nt,this%ns+1))
        allocate(this%blade(ib)%vind_wake2(3,nt,this%ns+1))
        allocate(this%blade(ib)%vind_wake3(3,nt,this%ns+1))
        allocate(this%blade(ib)%Pvind_wake(3,nt,this%ns+1))
        allocate(this%blade(ib)%vind_wake_step(3,nt,this%ns+1))
      end select
    enddo

    ! Wake initialization
    ! Assign core_radius to mid vortices
    do ib=1,this%nb
      do i=1,4
        this%blade(ib)%waP%vr%vf(i)%r_vc0 = this%spanwise_core
        this%blade(ib)%waP%vr%vf(i)%r_vc  = this%spanwise_core
        this%blade(ib)%waP%vr%vf(i)%age=0._dp
      enddo

      this%blade(ib)%waP%vr%gam=0._dp

      ! Assign core_radius to tip vortices
      do i=1,nt
        ! Root vortex
        this%blade(ib)%waP(i,1)%vr%vf(1)%r_vc0      = this%streamwise_core
        this%blade(ib)%waP(i,1)%vr%vf(1)%r_vc       = this%streamwise_core
        this%blade(ib)%waP(i,1)%vr%vf(3)%r_vc0      = this%streamwise_core
        this%blade(ib)%waP(i,1)%vr%vf(3)%r_vc       = this%streamwise_core

        this%blade(ib)%waP(i,2)%vr%vf(1)%r_vc0      = this%streamwise_core
        this%blade(ib)%waP(i,2)%vr%vf(1)%r_vc       = this%streamwise_core
        this%blade(ib)%waP(i,2)%vr%vf(3)%r_vc0      = this%streamwise_core
        this%blade(ib)%waP(i,2)%vr%vf(3)%r_vc       = this%streamwise_core

        ! Tip vortex
        this%blade(ib)%waP(i,this%ns)%vr%vf(1)%r_vc0   = this%streamwise_core
        this%blade(ib)%waP(i,this%ns)%vr%vf(1)%r_vc    = this%streamwise_core
        this%blade(ib)%waP(i,this%ns)%vr%vf(3)%r_vc0   = this%streamwise_core
        this%blade(ib)%waP(i,this%ns)%vr%vf(3)%r_vc    = this%streamwise_core

        this%blade(ib)%waP(i,this%ns-1)%vr%vf(3)%r_vc0 = this%streamwise_core
        this%blade(ib)%waP(i,this%ns-1)%vr%vf(3)%r_vc  = this%streamwise_core
        !this%blade(ib)%waP(i,this%ns-1)%vr%vf(3)%r_vc0 = this%streamwise_core
        !this%blade(ib)%waP(i,this%ns-1)%vr%vf(3)%r_vc  = this%streamwise_core
      enddo

      !if (starting_vortex_core > eps) then
      !  ! Assign core_radius to starting vortices
      !  do i=1,this%ns
      !    do j=2,4,2
      !      this%blade(ib)%waP(rows,i)%vr%vf(j)%r_vc0 = starting_vortex_core
      !      this%blade(ib)%waP(rows,i)%vr%vf(j)%r_vc  = starting_vortex_core
      !      this%blade(ib)%waP(rows-1,i)%vr%vf(j)%r_vc0 = starting_vortex_core
      !      this%blade(ib)%waP(rows-1,i)%vr%vf(j)%r_vc  = starting_vortex_core
      !    enddo
      !  enddo
      !endif
    enddo

  end subroutine init

  subroutine deinit(this,FDscheme_switch)
  class(rotor_class) :: this
    integer, intent(in) :: FDscheme_switch
    integer :: ib
    ! Deallocate variables
    do ib=1,this%nb
      select case (FDscheme_switch)
      case (0)
        deallocate(this%blade(ib)%vind_wake)
      case (1)
        deallocate(this%blade(ib)%Pwake)
        deallocate(this%blade(ib)%vind_wake)
        deallocate(this%blade(ib)%vind_wake1)
        deallocate(this%blade(ib)%Pvind_wake)
      case (2)
        deallocate(this%blade(ib)%vind_wake)
        deallocate(this%blade(ib)%vind_wake1)
        deallocate(this%blade(ib)%vind_wake_step)
      case (3)
        deallocate(this%blade(ib)%Pwake)
        deallocate(this%blade(ib)%vind_wake)
        deallocate(this%blade(ib)%vind_wake1)
        deallocate(this%blade(ib)%vind_wake2)
        deallocate(this%blade(ib)%vind_wake3)
        deallocate(this%blade(ib)%Pvind_wake)
        deallocate(this%blade(ib)%vind_wake_step)
      end select
    enddo

  end subroutine deinit

  function gettheta(this,psi,ib)
  class(rotor_class) :: this
    real(dp), intent(in) :: psi
    integer, intent(in) :: ib
    real(dp) :: gettheta
    real(dp) :: blade_offset

    blade_offset=2._dp*pi/this%nb*(ib-1)
    gettheta=this%control_pitch(1)  &
      +         this%control_pitch(2)*cos(psi+blade_offset)  &
      +         this%control_pitch(3)*sin(psi+blade_offset)  

  end function gettheta

  function getthetadot(this,psi,ib)
  class(rotor_class) :: this
    real(dp), intent(in) :: psi
    integer, intent(in) :: ib
    real(dp) :: getthetadot
    real(dp) :: blade_offset

    blade_offset=2._dp*pi/this%nb*(ib-1)
    getthetadot=-this%control_pitch(2)*sin(psi+blade_offset)  &
      +          this%control_pitch(3)*cos(psi+blade_offset)  

  end function getthetadot

  subroutine calcAIC(this)
  class(rotor_class), intent(inout) :: this
    integer :: ib,jblade,is,ic,i,j,row,col
    real(dp), dimension(3) :: vec_dummy

    ! Influence Coefficient Matrix
    do ib=1,this%nb
      do is=1,this%ns      ! Collocation point loop
        do ic=1,this%nc
          row=ic+this%nc*(is-1)+this%ns*this%nc*(ib-1)

          do jblade=1,this%nb
            do j=1,this%ns       ! Vortex ring loop
              do i=1,this%nc
                col=i+this%nc*(j-1)+this%ns*this%nc*(jblade-1)
                vec_dummy=this%blade(jblade)%wiP(i,j)%vr%vind(this%blade(ib)%wiP(ic,is)%CP)
                this%AIC(row,col)=dot_product(vec_dummy,this%blade(ib)%wiP(ic,is)%ncap)
              enddo
            enddo
          enddo

        enddo
      enddo
    enddo
    this%AIC_inv=inv(this%AIC)
  end subroutine calcAIC

  subroutine map_gam(this)
  class(rotor_class), intent(inout) :: this
    integer :: ib
    do ib=1,this%nb
      this%blade(ib)%wiP%vr%gam  &
        =reshape(this%gamvec(1+this%nc*this%ns*(ib-1):this%nc*this%ns*ib),(/this%nc,this%ns/))
    enddo
  end subroutine map_gam

  !-----+------------------+-----|
  ! -+- | Motion Functions | -+- |
  !-----+------------------+-----|

  subroutine rotor_move(this,dshift)
  class(rotor_class) :: this
    real(dp), intent(in), dimension(3) :: dshift

    integer :: ib

    do ib=1,this%nb
      call this%blade(ib)%move(dshift)
    enddo
    this%hub_coords=this%hub_coords+dshift
    this%CG_coords=this%CG_coords+dshift

  end subroutine rotor_move

  subroutine rotor_rot_pts(this,pts,origin,order)
  class(rotor_class), intent(inout) :: this
    real(dp), dimension(3), intent(in) :: pts    ! pts => phi,theta,psi
    real(dp), dimension(3), intent(in) :: origin ! rotation about
    integer, intent(in) :: order    ! [1]gb & +ve theta , [2]bg & -ve theta
    integer :: ib
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

    do ib=1,this%nb
      call this%blade(ib)%rot_pts(pts,origin,order)
    enddo

    this%shaft_axis=matmul(TMat,this%shaft_axis)

    this%hub_coords=this%hub_coords-origin
    this%hub_coords=matmul(TMat,this%hub_coords)
    this%hub_coords=this%hub_coords+origin

    this%CG_coords=this%CG_coords-origin
    this%CG_coords=matmul(TMat,this%CG_coords)
    this%CG_coords=this%CG_coords+origin

  end subroutine rotor_rot_pts

  subroutine rotor_rot_advance(this,dpsi)
  class(rotor_class), intent(inout) :: this
    real(dp), intent(in) :: dpsi
    integer :: ib
    real(dp) :: dtheta

    this%psi=this%psi+dpsi
    do ib=1,this%nb
      call this%blade(ib)%rot_axis(dpsi,this%shaft_axis,this%hub_coords)
      this%blade(ib)%psi=this%blade(ib)%psi+dpsi
      dtheta=this%gettheta(this%psi,ib)-this%blade(ib)%theta
      call this%blade(ib)%rot_pitch(dtheta)
      this%blade(ib)%theta=this%gettheta(this%psi,ib)
    enddo

  end subroutine rotor_rot_advance


  !-----+---------------------------+-----|
  ! -+- | Wake Convection Functions | -+- |
  !-----+---------------------------+-----|

  ! Assigns coordinates to first row of wake from last row of blade
  subroutine assignshed(this,row_now,edge)
  class(rotor_class), intent(inout) :: this
    integer, intent(in) :: row_now
    character(len=2), intent(in) :: edge
    integer :: i, ib

    do ib=1,this%nb
      this%blade(ib)%waP(row_now,:)%vr%gam=this%blade(ib)%wiP(this%nc,:)%vr%gam
    enddo

    select case (edge)
    case ('LE')    ! assign to LE
      do ib=1,this%nb
        do i=1,this%ns
          call this%blade(ib)%waP(row_now,i)%vr%assignP(1,this%blade(ib)%wiP(this%nc,i)%vr%vf(2)%fc(:,1))
          call this%blade(ib)%waP(row_now,i)%vr%assignP(4,this%blade(ib)%wiP(this%nc,i)%vr%vf(3)%fc(:,1))
          call this%blade(ib)%waP(row_now,i)%vr%calclength(.TRUE.)    ! TRUE => record original length
        enddo

      enddo
    case ('TE')    ! assign to TE
      do ib=1,this%nb
        do i=1,this%ns
          call this%blade(ib)%waP(row_now,i)%vr%assignP(2,this%blade(ib)%wiP(this%nc,i)%vr%vf(2)%fc(:,1))
          call this%blade(ib)%waP(row_now,i)%vr%assignP(3,this%blade(ib)%wiP(this%nc,i)%vr%vf(3)%fc(:,1))
        enddo
      enddo
    case default
      error stop 'Error: Wrong option for edge'
    end select

  end subroutine assignshed


  !-----+----------------------------+-----|
  ! -+- | Wake Dissipation Functions | -+- |
  !-----+----------------------------+-----|

  subroutine age_wake(this,row_now,dt)
  class(rotor_class), intent(inout) :: this
    integer, intent(in) :: row_now
    real(dp),intent(in) :: dt
    integer :: i, ib, row_last
    row_last=size(this%blade(1)%waP,1)
    do ib=1,this%nb
      !$omp parallel do
      do i=1,4
        this%blade(ib)%waP(row_now:row_last,:)%vr%vf(i)%age=this%blade(ib)%waP(row_now:row_last,:)%vr%vf(i)%age+dt
      enddo
      !$omp end parallel do
    enddo
  end subroutine age_wake

  subroutine dissipate_tip(this,row_now)
  class(rotor_class), intent(inout) :: this
    real(dp) :: oseen_param, turb_visc, kin_visc, new_radius
    integer, intent(in) :: row_now
    integer :: row_last
    integer :: ii,tip,ib
    oseen_param= 1.2564_dp
    kin_visc   = 0.0000181_dp
    turb_visc  = 500._dp

    row_last=size(this%blade(1)%waP,1)
    tip=this%ns
    do ib=1,this%nb
      do ii=row_now,row_last
        ! Root vortex core
        new_radius=sqrt(this%blade(ib)%waP(ii,1)%vr%vf(1)%r_vc**2._dp &
          +4._dp*oseen_param*turb_visc*kin_visc*this%blade(ib)%waP(ii,1)%vr%vf(1)%age)
        this%blade(ib)%waP(ii,1)%vr%vf(1)%r_vc=new_radius

        ! Tip vortex core
        new_radius=sqrt(this%blade(ib)%waP(ii,tip)%vr%vf(3)%r_vc**2._dp &
          +4._dp*oseen_param*turb_visc*kin_visc*this%blade(ib)%waP(ii,tip)%vr%vf(3)%age)
        this%blade(ib)%waP(ii,tip)%vr%vf(3)%r_vc=new_radius
      enddo
    enddo
  end subroutine dissipate_tip

  subroutine strain_wake(this,row_now)
  class(rotor_class), intent(inout) :: this
    integer, intent(in) :: row_now
    integer :: i,j,ib

    do ib=1,this%nb
      !$omp parallel do collapse(2)
      do j=1,this%ns
        do i=row_now,size(this%blade(1)%waP)
          call this%blade(ib)%waP(i,j)%vr%calclength(.FALSE.)    ! Update current length
          call this%blade(ib)%waP(i,j)%vr%strain()
        enddo
      enddo
      !$omp end parallel do
    enddo
  end subroutine strain_wake

  function rotor_vind_bywing(this,P)
  class(rotor_class), intent(inout) :: this
    real(dp), intent(in), dimension(3) :: P
    real(dp), dimension(3) :: rotor_vind_bywing
    integer :: ib

    rotor_vind_bywing=0._dp
    do ib=1,this%nb
      rotor_vind_bywing=rotor_vind_bywing+this%blade(ib)%vind_bywing(P)
    enddo
  end function rotor_vind_bywing

  function rotor_vind_bywake(this,row_now,P)
  class(rotor_class), intent(inout) :: this
    real(dp), intent(in), dimension(3) :: P
    real(dp), dimension(3) :: rotor_vind_bywake
    integer, intent(in) :: row_now
    integer :: ib

    rotor_vind_bywake=0._dp
    do ib=1,this%nb
      rotor_vind_bywake=rotor_vind_bywake+this%blade(ib)%vind_bywake(row_now,P)
    enddo
  end function rotor_vind_bywake
end module rotor_classdef
