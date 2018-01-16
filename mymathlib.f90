module mymathlib

  implicit none
  integer, parameter :: dp = kind(1.d0)
  real(dp), parameter :: pi = atan(1._dp)*4._dp
  real(dp), parameter :: eps = epsilon(1._dp)

contains

  ! -------------------------------------------------
  !                length3d
  ! -------------------------------------------------
  function length3d(P1,P2) result(length)
    real(dp), intent(in), dimension(3) :: P1, P2
    real(dp) :: length
    real(dp), dimension(3) :: delta
    delta=P1-P2
    length=norm2(delta)
  end function length3d

  ! -------------------------------------------------
  !                linspace
  ! -------------------------------------------------
  function linspace(xstart,xend,nx) result(xout)
    real(dp), intent(in) :: xstart, xend
    integer , intent(in) :: nx
    real(dp), dimension(nx) :: xout
    integer :: i
    real(dp) :: dx

    dx = (xend-xstart)/(nx-1)

    xout = (/((i*dx),i=0,nx-1)/)
    xout = xout+xstart

  end function linspace

  ! -------------------------------------------------
  !                cosspace
  ! -------------------------------------------------
  function cosspace(xstart,xend,nx) result(xout)
    real(dp), intent(in) :: xstart, xend
    integer , intent(in) :: nx
    real(dp), dimension(nx) :: xout
    real(dp), dimension(nx) :: theta_spacing

    theta_spacing=linspace(0._dp,pi,nx)
    xout=xstart+(xend-xstart)*0.5_dp*(1._dp-cos(theta_spacing))

  end function cosspace

  ! -------------------------------------------------
  !                cross3
  ! -------------------------------------------------
  function cross3(avec,bvec)
    real(dp), intent(in), dimension(3) :: avec,bvec
    real(dp), dimension(3) :: cross3

    cross3(1) = avec(2)*bvec(3)-avec(3)*bvec(2)
    cross3(2) = avec(3)*bvec(1)-avec(1)*bvec(3)
    cross3(3) = avec(1)*bvec(2)-avec(2)*bvec(1)

  end function cross3

  ! -------------------------------------------------
  !                degtorad
  ! -------------------------------------------------
  subroutine degtorad(deg)
    real(dp), intent(inout) :: deg
    deg=deg*pi/180._dp
  end subroutine degtorad

  ! -------------------------------------------------
  !                radtodeg
  ! -------------------------------------------------
  subroutine radtodeg(rad) 
    real(dp), intent(inout) :: rad
    rad=rad*180._dp/pi
  end subroutine radtodeg

  ! -------------------------------------------------
  !                outer_product
  ! -------------------------------------------------
  function outer_product(vecA,vecB)
    real(dp), intent(in), dimension(3) :: vecA, vecB
    real(dp), dimension(3,3) :: outer_product
    outer_product(:,1)=(/vecA(1)*vecB(1),vecA(2)*vecB(1),vecA(3)*vecB(1)/)
    outer_product(:,2)=(/vecA(1)*vecB(2),vecA(2)*vecB(2),vecA(3)*vecB(2)/)
    outer_product(:,3)=(/vecA(1)*vecB(2),vecA(2)*vecB(3),vecA(3)*vecB(3)/)
  end function outer_product

  ! -------------------------------------------------
  !                inv
  ! -------------------------------------------------
  ! Matrix Inversion
  ! Ax=b
  ! PA = LU or A=P'LU
  ! P'LUx=b
  ! LUx=Pb
  ! Solve Ld=Pb using Forward  sub where d=Ux
  ! Solve Ux=d using Backward sub 
  function inv(A) result(Ainv)
    real(dp), intent(in), dimension(:,:) :: A
    real(dp), dimension(size(A,1),size(A,1)) :: Ainv

    integer :: n
    integer :: i,j,k,bb ! Running variables

    ! Variables for calculating Permutation matrix
    real(dp) :: max_elem
    real(dp), dimension(size(A,1),size(A,1)) :: A_dummy
    real(dp), dimension(size(A,1)) :: P_swap, A_swap

    ! Variables for LU Decomposition
    real(dp), dimension(size(A,1),size(A,1)) :: L,U,P
    real(dp), dimension(size(A,1)) :: Pb,d,x,bvec
    real(dp) :: sumu, suml, det

    n=size(A,1)
    A_dummy=A

    ! Find Permutation Matrix for partial pivoting
    ! Creating P as Identity Matrix
    P=0._dp
    do i=1,n
      P(i,i)=1.
    enddo

    do j=1,n
      max_elem=maxval(A_dummy(j:n,j))
      do i=j,n
        if (A(i,j) - max_elem .lt. eps) then
          P_swap=P(i,:)
          P(i,:)=P(j,:)
          P(j,:)=P_swap

          A_swap=A_dummy(i,:)
          A_dummy(i,:)=A_dummy(j,:)
          A_dummy(j,:)=A_swap
          exit
        endif
      enddo
    enddo

    ! LU decomposition using Doolittle algorithm on PA
    ! A_dummy is now P*A
    U(1,:) = A_dummy(1,:)
    L(:,1) = A_dummy(:,1)/A_dummy(1,1)

    sumu=0._dp
    suml=0._dp

    do i=2,n
      do j=2,n
        sumu=0._dp
        suml=0._dp
        do k=1,i-1
          sumu=sumu+L(i,k)*U(k,j)
          suml=suml+L(j,k)*U(k,i)
        enddo
        U(i,j)= A_dummy(i,j)-sumu
        if (abs(U(i,i)) .gt. eps) L(j,i)=(A_dummy(j,i)-suml)/U(i,i)
      enddo
    enddo

    ! Assigning all zero elements in triangular matrices
    det=1._dp
    do i=1,n
      det=det*U(i,i)
      do j=1,n
        if (i>j) then
          U(i,j)=0._dp
        elseif(j>i) then
          L(i,j)=0._dp
        endif
      enddo
    enddo

    ! Checking Determinant for singularity
    if (abs(det)<eps) then
      print*
      print*,'ERROR: Matrix is Singular or Ill-conditioned!!'
      call print_mat(A)
      print*,'Determinant was found to be:'
      print*,det
      call print_mat(U)
      stop 404
    endif

    ! Changing RHS loop
    do bb=1,n
      bvec=0._dp
      bvec(bb)=1._dp

      Pb = matmul(P,bvec)
      d=0._dp
      x=0._dp

      ! Forward Substitution
      d(1) = Pb(1)
      do i=2,n
        suml=0._dp
        do k=1,i-1
          suml=suml+L(i,k)*d(k)
        enddo
        d(i)=Pb(i)-suml
      enddo

      ! Backward Substitution
      x(n)=d(n)/U(n,n)
      do i=n-1,1,-1
        sumu=0._dp
        do k=i+1,n
          sumu=sumu+U(i,k)*x(k)
        enddo
        x(i)=(d(i)-sumu)/U(i,i)
      enddo

      Ainv(:,bb) = x
    enddo

  end function inv

  ! -------------------------------------------------
  !                print_mat
  ! -------------------------------------------------
  ! Display in matrix format
  subroutine print_mat(M)
    real(dp), intent(in), dimension(:,:)     :: M
    real(dp), dimension(size(M,1),size(M,2)) :: M_dummy
    integer :: i,j
    M_dummy=M 
    do i=1,size(M,1)
      do j=1,size(M,1)
        if (abs(M(i,j))<eps) M_dummy(i,j)=0.0
        write(*,100,advance='no') M_dummy(i,j)
      enddo
      write(*,*)
    enddo
    100 format(ES14.3)
  end subroutine print_mat

  ! -------------------------------------------------
  !                norm
  ! -------------------------------------------------
  function norm(abcvec)
    real(dp), intent(in), dimension(:) :: abcvec
    real(dp) :: norm
    integer :: is
    norm =0._dp
    do is=1,size(abcvec)
      norm = norm + abcvec(is)*abcvec(is)
    enddo
    norm=sqrt(norm)
  end function norm



end module mymathlib
