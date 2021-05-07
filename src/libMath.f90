module libMath

  implicit none
  integer, parameter :: dp = kind(1.d0)
  real(dp), parameter :: pi = atan(1._dp)*4._dp
  real(dp), parameter :: eps = epsilon(1._dp)

  real(dp), parameter, dimension(3) :: xAxis = (/1._dp, 0._dp, 0._dp/)
  real(dp), parameter, dimension(3) :: yAxis = (/0._dp, 1._dp, 0._dp/)
  real(dp), parameter, dimension(3) :: zAxis = (/0._dp, 0._dp, 1._dp/)

  character(len=1), parameter :: commentChar = '#'

  interface lsq2
    module procedure lsq2_scalar, lsq2_array
  end interface

contains

  ! -------------------------------------------------
  !                length3d
  ! -------------------------------------------------
  function length3d(P1, P2) result(length)
    real(dp), intent(in), dimension(3) :: P1, P2
    real(dp) :: length
    real(dp), dimension(3) :: delta
    delta = P1 - P2
    length = norm2(delta)
  end function length3d

  ! -------------------------------------------------
  !                linspace
  ! -------------------------------------------------
  function linspace(xstart, xend, nx) result(xout)
    real(dp), intent(in) :: xstart, xend
    integer, intent(in) :: nx
    real(dp), dimension(nx) :: xout
    integer :: i
    real(dp) :: dx

    dx = (xend - xstart)/(nx - 1)

    xout = (/((i*dx), i=0, nx - 1)/)
    xout = xout + xstart

  end function linspace

  ! -------------------------------------------------
  !                cosspace
  ! -------------------------------------------------
  function cosspace(xstart, xend, nx) result(xout)
    real(dp), intent(in) :: xstart, xend
    integer, intent(in) :: nx
    real(dp), dimension(nx) :: xout
    real(dp), dimension(nx) :: theta_spacing

    theta_spacing = linspace(0._dp, pi, nx)
    xout = xstart + (xend - xstart) * 0.5_dp * &
      (1._dp - cos(theta_spacing))

  end function cosspace

  ! -------------------------------------------------
  !                halfsinspace
  ! -------------------------------------------------
  function halfsinspace(xstart, xend, nx) result(xout)
    real(dp), intent(in) :: xstart, xend
    integer, intent(in) :: nx
    real(dp), dimension(nx) :: xout
    real(dp), dimension(nx) :: theta_spacing

    theta_spacing = linspace(0._dp, pi*0.5_dp, nx)
    xout = xstart + (xend - xstart)*sin(theta_spacing)

  end function halfsinspace

  ! -------------------------------------------------
  !                cross_product
  ! -------------------------------------------------
  function cross_product(aVec, bVec)
    real(dp), intent(in), dimension(3) :: aVec, bVec
    real(dp), dimension(3) :: cross_product

    cross_product(1) = aVec(2)*bVec(3) - aVec(3)*bVec(2)
    cross_product(2) = aVec(3)*bVec(1) - aVec(1)*bVec(3)
    cross_product(3) = aVec(1)*bVec(2) - aVec(2)*bVec(1)

  end function cross_product

  ! -------------------------------------------------
  !                degtorad
  ! -------------------------------------------------
  subroutine degtorad(deg)
    real(dp), intent(inout) :: deg
    deg = deg*pi/180._dp
  end subroutine degtorad

  ! -------------------------------------------------
  !                radtodeg
  ! -------------------------------------------------
  subroutine radtodeg(rad)
    real(dp), intent(inout) :: rad
    rad = rad*180._dp/pi
  end subroutine radtodeg

  ! -------------------------------------------------
  !                outer_product
  ! -------------------------------------------------
  function outer_product(aVec, bVec)
    real(dp), intent(in), dimension(3) :: aVec, bVec
    real(dp), dimension(3, 3) :: outer_product
    outer_product(:, 1) = (/aVec(1)*bVec(1), aVec(2)*bVec(1), aVec(3)*bVec(1)/)
    outer_product(:, 2) = (/aVec(1)*bVec(2), aVec(2)*bVec(2), aVec(3)*bVec(2)/)
    outer_product(:, 3) = (/aVec(1)*bVec(2), aVec(2)*bVec(3), aVec(3)*bVec(3)/)
  end function outer_product

  ! -------------------------------------------------
  !                getAngleTan
  ! -------------------------------------------------
  function getAngleTan(aVec, bVec)
    real(dp), intent(in), dimension(3) :: aVec, bVec
    real(dp) :: getAngleTan
    real(dp) :: mag2A, mag2B, dotAB

    ! Result will be -pi to pi
    mag2A = sum(aVec**2._dp)
    mag2B = sum(bVec**2._dp)
    dotAB = dot_product(aVec, bVec)
    getAngleTan = atan2(sqrt(mag2A*mag2B - dotAB*dotAB), dotAB)
  end function getAngleTan

  ! -------------------------------------------------
  !                getAngleCos
  ! -------------------------------------------------
  function getAngleCos(aVec, bVec)
    real(dp), intent(in), dimension(3) :: aVec, bVec
    real(dp) :: getAngleCos

    ! Assumes no angle > 180 deg exists
    getAngleCos = acos(dot_product(aVec, bVec)/ &
                       (sqrt(sum(aVec**2._dp)*sum(bVec**2._dp))))
  end function getAngleCos

  ! -------------------------------------------------
  !                unitVec
  ! -------------------------------------------------
  function unitVec(aVec)
    real(dp), intent(in), dimension(:) :: aVec
    real(dp), dimension(size(aVec)) :: unitVec
    real(dp) :: normVal

    normVal = norm2(aVec)
    if (normVal > eps) then
      unitVec = aVec/normVal
    else
      unitVec = 0.0
    endif
  end function unitVec

  ! -------------------------------------------------
  !                projVec
  ! -------------------------------------------------
  function projVec(aVec, dirVec)
    real(dp), intent(in), dimension(3) :: aVec, dirVec
    real(dp), dimension(3) :: projVec
    real(dp) :: normSq
    normSq = sum(dirVec**2._dp)
    if (normSq > eps) then
      projVec = dot_product(aVec, dirVec)*dirVec/normSq
    else
      projVec = 0._dp
    endif
  end function projVec

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
    real(dp), intent(in), dimension(:, :) :: A
    real(dp), dimension(size(A, 1), size(A, 1)) :: Ainv

    integer :: n
    integer :: i, j, k, bb ! Running variables

    ! Variables for calculating Permutation matrix
    real(dp) :: max_elem
    real(dp), dimension(size(A, 1), size(A, 1)) :: A_dummy
    real(dp), dimension(size(A, 1)) :: P_swap, A_swap

    ! Variables for LU Decomposition
    real(dp), dimension(size(A, 1), size(A, 1)) :: L, U, P
    real(dp), dimension(size(A, 1)) :: Pb, d, x, bvec
    real(dp) :: sumu, suml
    real(dp), dimension(size(A, 1)) :: diagonalTerms

    n = size(A, 1)
    A_dummy = A

    ! Find Permutation Matrix for partial pivoting
    ! Creating P as Identity Matrix
    P = 0._dp
    do i = 1, n
      P(i, i) = 1.
    enddo

    do j = 1, n
      max_elem = maxval(A_dummy(j:n, j))
      do i = j, n
        if (A(i, j) - max_elem .lt. eps) then
          P_swap = P(i, :)
          P(i, :) = P(j, :)
          P(j, :) = P_swap

          A_swap = A_dummy(i, :)
          A_dummy(i, :) = A_dummy(j, :)
          A_dummy(j, :) = A_swap
          exit
        endif
      enddo
    enddo

    ! LU decomposition using Doolittle algorithm on PA
    ! A_dummy is now P*A
    U(1, :) = A_dummy(1, :)
    L(:, 1) = A_dummy(:, 1)/A_dummy(1, 1)

    sumu = 0._dp
    suml = 0._dp

    do i = 2, n
      do j = 2, n
        sumu = 0._dp
        suml = 0._dp
        do k = 1, i - 1
          sumu = sumu + L(i, k)*U(k, j)
          suml = suml + L(j, k)*U(k, i)
        enddo
        U(i, j) = A_dummy(i, j) - sumu
        if (abs(U(i, i)) .gt. eps) L(j, i) = (A_dummy(j, i) - suml)/U(i, i)
      enddo
    enddo

    ! Assigning all zero elements in triangular matrices
    do i = 1, n
      diagonalTerms(i) = U(i, i)
      do j = 1, n
        if (i > j) then
          U(i, j) = 0._dp
        elseif (j > i) then
          L(i, j) = 0._dp
        endif
      enddo
    enddo

    ! Checking diagonal elements for zero
    ! If determinant is computed here by multiplication,
    ! for large matrices it may produce floating point overflow
    do i = 1, n
      if (abs(diagonalTerms(i)) < eps) then
        print *
        print *, 'ERROR: Matrix is Singular or Ill-conditioned!!'
        print *, 'A-matrix:'
        call print_mat(A)
        print *, 'U-matrix:'
        call print_mat(U)
        stop 404
      endif
    enddo

    ! Changing RHS loop
    do bb = 1, n
      bvec = 0._dp
      bvec(bb) = 1._dp

      Pb = matmul(P, bvec)
      d = 0._dp
      x = 0._dp

      ! Forward Substitution
      d(1) = Pb(1)
      do i = 2, n
        suml = 0._dp
        do k = 1, i - 1
          suml = suml + L(i, k)*d(k)
        enddo
        d(i) = Pb(i) - suml
      enddo

      ! Backward Substitution
      x(n) = d(n)/U(n, n)
      do i = n - 1, 1, -1
        sumu = 0._dp
        do k = i + 1, n
          sumu = sumu + U(i, k)*x(k)
        enddo
        x(i) = (d(i) - sumu)/U(i, i)
      enddo

      Ainv(:, bb) = x
    enddo

  end function inv

  function isInverse(A, Ainv)
    logical :: isInverse
    real(dp), intent(in), dimension(:, :) :: A, Ainv
    real(dp), dimension(size(A, 1), size(A, 2)) :: productMat
    integer :: i, j
    real(dp) :: tol

    productMat = matmul(A, Ainv)
    isInverse = .TRUE.
    tol = 1E-04

    do j = 1, size(A, 2)
      do i = 1, size(A, 1)
        if (i .ne. j) then
          ! Check if off-diagonal values are 0._dp
          if (productMat(i, j) > tol) then
            isInverse = .FALSE.
            !print*,i,j,'Off-diagonal values non-zero',productMat(i,j)
          endif
        else
          ! Check if on-diagonal values are 1._dp
          if (productMat(i, j) - 1._dp > tol) then
            isInverse = .FALSE.
            !print*,i,'Diagonal values non-unity',productMat(i,j)
          endif
        endif
      enddo
    enddo
  end function isInverse

  ! -------------------------------------------------
  !                print_mat
  ! -------------------------------------------------
  ! Display in matrix format
  subroutine print_mat(M)
    real(dp), intent(in), dimension(:, :)     :: M
    real(dp), dimension(size(M, 1), size(M, 2)) :: M_dummy
    integer :: i, j
    M_dummy = M
    do i = 1, size(M, 1)
      do j = 1, size(M, 1)
        if (abs(M(i, j)) < eps) M_dummy(i, j) = 0.0
        write (*, 100, advance='no') M_dummy(i, j)
      enddo
      write (*, *)
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
    norm = 0._dp
    do is = 1, size(abcvec)
      norm = norm + abcvec(is)*abcvec(is)
    enddo
    norm = sqrt(norm)
  end function norm

  ! -------------------------------------------------
  !                interp1
  ! -------------------------------------------------
  function interp1(xq, x, y, order)
    real(dp), intent(in) :: xq
    real(dp), intent(in), dimension(:) :: x, y
    integer, intent(in) :: order
    real(dp) :: interp1
    logical, dimension(size(x)) :: TFvec
    integer :: i, ix
    real(dp) :: L0, L1, L2  !! Lagrange's basis functions

    TFVec = (xq <= x)
    do i = 1, size(x)
      if (TFvec(i)) then
        ix = i
        exit
      endif
    enddo

    if (abs(xq-x(ix)) <= eps) then
      interp1 = y(ix)
    else
      select case (order)
      case (1)
        if (norm2(y(ix-1:ix)) <= eps) then
          interp1 = 0._dp
        else
          interp1 = y(ix-1) + (xq-x(ix-1))*(y(ix)-y(ix-1))/(x(ix)-x(ix-1))
        endif

      case (2)
          if (norm2(y(ix-1:ix)) <= eps) then
            interp1 = 0._dp
          else
            if (size(x) == 2) then
              interp1 = y(ix-1) + (xq-x(ix-1))*(y(ix)-y(ix-1))/(x(ix)-x(ix-1))
            elseif (ix == 2) then
              L0 = (xq-x(ix))*(xq-x(ix+1))/((x(ix-1)-x(ix))*(x(ix-1)-x(ix+1)))
              L1 = (xq-x(ix-1))*(xq-x(ix+1))/((x(ix)-x(ix-1))*(x(ix)-x(ix+1)))
              L2 = (xq-x(ix-1))*(xq-x(ix))/((x(ix+1)-x(ix-1))*(x(ix+1)-x(ix)))
              interp1 = y(ix-1)*L0 + y(ix)*L1 + y(ix+1)*L2
            else
              L0 = (xq-x(ix-1))*(xq-x(ix))/((x(ix-2)-x(ix-1))*(x(ix-2)-x(ix)))
              L1 = (xq-x(ix-2))*(xq-x(ix))/((x(ix-1)-x(ix-2))*(x(ix-1)-x(ix)))
              L2 = (xq-x(ix-2))*(xq-x(ix-1))/((x(ix)-x(ix-2))*(x(ix)-x(ix-1)))
              interp1 = y(ix-2)*L0 + y(ix-1)*L1 + y(ix)*L2
            endif
          endif

        case default 
          error stop "Specified order not implemented"

      end select
    endif

  end function interp1

  !--------------------------------------------------------!
  !        Linear Least Squares fitting (2nd order)        !
  !--------------------------------------------------------!
  function lsq2_scalar(xQuery, xData, yData)
    real(dp), intent(in) :: xQuery
    real(dp), intent(in), dimension(:) :: xData, yData
    real(dp), dimension(3) :: coeff, RHS
    real(dp), dimension(3, 3) :: Amat
    real(dp) :: lsq2_scalar

    if (size(xData) .ne. size(yData)) error stop 'ERROR: size of xData and yData have to be equal'

    Amat(1, 1) = size(xData)
    Amat(1, 2) = sum(xData)
    Amat(1, 3) = sum(xData**2._dp)
    Amat(2, 1) = Amat(1, 2)
    Amat(2, 2) = Amat(1, 3)
    Amat(2, 3) = sum(xData**3._dp)
    Amat(3, 1) = Amat(1, 3)
    Amat(3, 2) = Amat(2, 3)
    Amat(3, 3) = sum(xData**4._dp)

    RHS(1) = sum(yData)
    RHS(2) = sum(yData*xData)
    RHS(3) = sum(yData*xData**2._dp)

    coeff = matmul(inv(Amat), RHS)

    lsq2_scalar = coeff(1) + coeff(2)*xQuery + coeff(3)*xQuery*xQuery
  end function lsq2_scalar

  function lsq2_array(xQuery, xData, yData)
    real(dp), intent(in), dimension(:) :: xQuery
    real(dp), intent(in), dimension(:) :: xData, yData
    real(dp), dimension(3) :: coeff, RHS
    real(dp), dimension(3, 3) :: Amat
    real(dp), dimension(size(xQuery)) :: lsq2_array
    integer :: i

    if (size(xData) .ne. size(yData)) error stop 'ERROR: size of xData and yData have to be equal'

    Amat(1, 1) = size(xData)
    Amat(1, 2) = sum(xData)
    Amat(1, 3) = sum(xData**2._dp)
    Amat(2, 1) = Amat(1, 2)
    Amat(2, 2) = Amat(1, 3)
    Amat(2, 3) = sum(xData**3._dp)
    Amat(3, 1) = Amat(1, 3)
    Amat(3, 2) = Amat(2, 3)
    Amat(3, 3) = sum(xData**4._dp)

    RHS(1) = sum(yData)
    RHS(2) = sum(yData*xData)
    RHS(3) = sum(yData*xData**2._dp)

    coeff = matmul(inv(Amat), RHS)

    do i = 1, size(xQuery)
      lsq2_array(i) = coeff(1) + coeff(2)*xQuery(i) + coeff(3)*xQuery(i)*xQuery(i)
    enddo
  end function lsq2_array

  !--------------------------------------------------------!
  !              Transformation Functions                  !
  !--------------------------------------------------------!
  ! Code to generate Transformation matrices in Octave
  !
  ! clc; clear;
  ! pkg load symbolic;
  ! syms p t s
  ! Rp=[[1,0,0];[0,cos(p),sin(p)];[0,-sin(p),cos(p)]];
  ! Rt=[[cos(t),0,-sin(t)];[0,1,0];[sin(t),0,cos(t)]];
  ! Rs=[[cos(s),sin(s),0];[-sin(s),cos(s),0];[0,0,1]];
  ! Tbg=Rp*Rt*Rs
  ! Tgb=Rs'*Rt'*Rp'

  ! Transformation matrix bg
  function Tbg(cs_phi, cs_theta, cs_psi)
    real(dp), dimension(2), intent(in) :: cs_phi, cs_theta, cs_psi  ! cos and sin
    real(dp), dimension(3, 3) :: Tbg
    Tbg(1, :) = (/cs_psi(1)*cs_theta(1), cs_theta(1)*cs_psi(2), -1._dp*cs_theta(2)/)
    Tbg(2, 1) = cs_psi(1)*cs_phi(2)*cs_theta(2) - cs_phi(1)*cs_psi(2)
    Tbg(2, 2) = cs_phi(1)*cs_psi(1) + cs_phi(2)*cs_psi(2)*cs_theta(2)
    Tbg(2, 3) = cs_theta(1)*cs_phi(2)
    Tbg(3, 1) = cs_phi(1)*cs_psi(1)*cs_theta(2) + cs_phi(2)*cs_psi(2)
    Tbg(3, 2) = cs_phi(1)*cs_psi(2)*cs_theta(2) - cs_psi(1)*cs_phi(2)
    Tbg(3, 3) = cs_phi(1)*cs_theta(1)
  end function Tbg

  function Tgb(cs_phi, cs_theta, cs_psi)
    real(dp), dimension(2), intent(in) :: cs_phi, cs_theta, cs_psi  
    ! cs_ represent cos and sin
    real(dp), dimension(3, 3) :: Tgb
    Tgb(1, 1) = cs_psi(1)*cs_theta(1)
    Tgb(1, 2) = cs_phi(2)*cs_theta(2)*cs_psi(1) - cs_psi(2)*cs_phi(1)
    Tgb(1, 3) = cs_phi(2)*cs_psi(2) + cs_theta(2)*cs_phi(1)*cs_psi(1)
    Tgb(2, 1) = cs_psi(2)*cs_theta(1)
    Tgb(2, 2) = cs_phi(2)*cs_psi(2)*cs_theta(2) + cs_phi(1)*cs_psi(1)
    Tgb(2, 3) = cs_psi(2)*cs_theta(2)*cs_phi(1) - cs_phi(2)*cs_psi(1)
    Tgb(3, 1) = -cs_theta(2)
    Tgb(3, 2) = cs_phi(2)*cs_theta(1)
    Tgb(3, 3) = cs_phi(1)*cs_theta(1)
  end function Tgb

  function trapz(y, x)
    !! Trapezoid integration unequal intervals
    real(dp) :: trapz
    real(dp), intent(in), dimension(:) :: x, y
    integer :: i

    if (size(x) .ne. size(y)) then
      error stop 'ERROR: Sizes of vectors do not match'
    endif

    trapz = 0._dp
    do i = 2, size(x)
      trapz = trapz + (x(i)-x(i-1))*(y(i)+y(i-1))
    enddo
    trapz = 0.5_dp*trapz
  end function trapz

  ! Skips comments in input file
  ! Comment character is set using the global var commentChar
  subroutine skip_comments(fileUnit)
    integer, intent(in) :: fileUnit
    character(len=1) :: firstChar

    firstChar = commentChar
    do while (firstChar .eq. commentChar)
      read(fileUnit, '(A)') firstChar
    enddo
    backspace(fileUnit)
  end subroutine skip_comments

end module libMath
