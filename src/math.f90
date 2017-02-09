!   This file is part of EmDee.
!
!    EmDee is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    EmDee is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with EmDee. If not, see <http://www.gnu.org/licenses/>.
!
!    Author: Charlles R. A. Abreu (abreu@eq.ufrj.br)
!            Applied Thermodynamics and Molecular Simulation
!            Federal University of Rio de Janeiro, Brazil

module math

use global
use, intrinsic :: ieee_arithmetic

implicit none

real(rb), parameter, private :: a1 =  0.254829592_rb, &
                                a2 = -0.284496736_rb, &
                                a3 =  1.421413741_rb, &
                                a4 = -1.453152027_rb, &
                                a5 =  1.061405429_rb, &
                                p  =  0.327591100_rb

type, abstract :: i32rng
  logical :: seeding_required = .true.
  integer(4), private :: kn(0:127)
  real(8),    private :: wn(0:127), fn(0:127)
  contains
    procedure(i32rng_init), deferred, private :: init
    procedure(i32rng_i32),  deferred          :: i32
    procedure :: setup     => i32rng_setup
    procedure :: normal    => i32rng_normal
end type i32rng

abstract interface

  subroutine i32rng_init( a, seed )
    import :: i32rng
    class(i32rng), intent(inout) :: a
    integer,       intent(in)    :: seed
  end subroutine i32rng_init

  function i32rng_i32( a ) result( i32 )
    import :: i32rng
    class(i32rng), intent(inout) :: a
    integer                      :: i32
  end function i32rng_i32

end interface

type, extends(i32rng) :: kiss
  integer, private :: x, y, z, w
  contains
    procedure :: init => kiss_init
    procedure :: i32  => kiss_i32
end type kiss

interface operator (.dot.)
  module procedure complex_dot_product
end interface

interface operator (.cross.)
  module procedure complex_cross_product
end interface

contains

!---------------------------------------------------------------------------------------------------

  subroutine i32rng_setup( a, seed )
    class(i32rng), intent(inout) :: a
    integer(4),    intent(in)    :: seed
    integer :: i
    real(8), parameter :: m1 = 2147483648.0_8
    real(8) :: q, dn, tn, vn
    call a%init( seed )
    a%seeding_required = .false.
    dn = 3.442619855899_8
    tn = 3.442619855899_8
    vn = 0.00991256303526217_8
    q = vn*exp(0.5_8*dn*dn)
    a%kn(0) = int((dn/q)*m1)
    a%kn(1) = 0_4
    a%wn(0) = q/m1
    a%wn(127) = dn/m1
    a%fn(0) = 1.0_8
    a%fn(127) = exp( -0.5_8*dn*dn )
    do i = 126, 1, -1
      dn = sqrt( -2.0_8 * log( vn/dn + exp( -0.5_8*dn*dn ) ) )
      a%kn(i+1) = int((dn/tn)*m1)
      tn = dn
      a%fn(i) = exp(-0.5_8*dn*dn)
      a%wn(i) = dn/m1
    end do
  end subroutine i32rng_setup

!---------------------------------------------------------------------------------------------------

  function i32rng_normal( a ) result( rnor )
    class(i32rng), intent(inout) :: a
    real(8)                      :: rnor

    real(8), parameter :: r = 3.442620_8, s = 0.2328306e-9_8
    integer(4) :: hz, iz
    real(8) :: x, y

    hz = a%i32()
    iz = iand(hz, 127_4)
    if (abs(hz) < a%kn(iz)) then
      rnor = hz*a%wn(iz)
    else
      do
        if (iz == 0_4) then
          do
            x = -0.2904764_8*log(s*a%i32() + 0.5_8)
            y = -log(s*a%i32() + 0.5_8)
            if (y+y >= x*x) exit
          end do
          rnor = r + x
          if (hz <= 0_4) rnor = -rnor
          return
        end if
        x = hz * a%wn(iz)
        if ( a%fn(iz) + (s*a%i32() + 0.5_8)*(a%fn(iz-1) - a%fn(iz)) < exp(-0.5_8*x*x) ) then
          rnor = x
          return
        end if
        hz = a%i32()
        iz = iand(hz, 127_4)
        if (abs(hz) < a%kn(iz)) then
          rnor = hz*a%wn(iz)
          return
        end if
      end do
   end if
  end function i32rng_normal

!---------------------------------------------------------------------------------------------------

  subroutine kiss_init( a, seed )
    class(kiss), intent(inout) :: a
    integer(4),  intent(in)    :: seed
    a%x = m(m(m(seed,13),-17),5)
    a%y = m(m(m(a%x, 13),-17),5)
    a%z = m(m(m(a%y, 13),-17),5)
    a%w = m(m(m(a%z, 13),-17),5)
    contains
      integer(4) function m( k, n )
        integer(4), intent(in) :: k
        integer,    intent(in) :: n
        m = ieor(k,ishft(k,n))
      end function m
  end subroutine kiss_init

!---------------------------------------------------------------------------------------------------

  function kiss_i32( a ) result( i32 )
    class(kiss), intent(inout) :: a
    integer                    :: i32
    a%x = 69069_4*a%x + 1327217885_4
    a%y = ieor(a%y,ishft(a%y, 13))
    a%y = ieor(a%y,ishft(a%y,-17))
    a%y = ieor(a%y,ishft(a%y,  5))
    a%z = 18000_4*iand(a%z,65535_4) + ishft(a%z,-16)
    a%w = 30903_4*iand(a%w,65535_4) + ishft(a%w,-16)
    i32 = a%x + a%y + ishft(a%z,16) + a%w
  end function kiss_i32

!---------------------------------------------------------------------------------------------------

  pure function cross_product(a, b) result( c )
    real(rb), intent(in) :: a(3), b(3)
    real(rb)             :: c(3)
    c = [ a(2)*b(3) - a(3)*b(2), a(3)*b(1) - a(1)*b(3), a(1)*b(2) - a(2)*b(1) ]
  end function cross_product

!---------------------------------------------------------------------------------------------------

  pure function quaternion( A ) result( Q )
    real(rb), intent(in) :: A(3,3)
    real(rb)             :: Q(4)

    ! Reference: S. W. Shepperd, Journal of Guidance and Control 1, 223 (1978)

    integer :: imax
    real(rb) :: Q2(4), a11, a22, a33, Q2max

    a11 = A(1,1)
    a22 = A(2,2)
    a33 = A(3,3)
    Q2 = [ one + a11 + a22 + a33, &
           one + a11 - a22 - a33, &
           one - a11 + a22 - a33, &
           one - a11 - a22 + a33  ]
    imax = maxloc(Q2,1)
    Q2max = Q2(imax)
    select case(imax)
      case (1)
        Q = 0.5_rb*[Q2max, A(2,3) - A(3,2), A(3,1) - A(1,3), A(1,2) - A(2,1)]*sqrt(one/Q2max)
      case (2)
        Q = 0.5_rb*[A(2,3) - A(3,2), Q2max, A(1,2) + A(2,1), A(1,3) + A(3,1)]*sqrt(one/Q2max)
      case (3)
        Q = 0.5_rb*[A(3,1) - A(1,3), A(1,2) + A(2,1), Q2max, A(2,3) + A(3,2)]*sqrt(one/Q2max)
      case (4)
        Q = 0.5_rb*[A(1,2) - A(2,1), A(1,3) + A(3,1), A(2,3) + A(3,2), Q2max]*sqrt(one/Q2max)
    end select

  end function quaternion

!---------------------------------------------------------------------------------------------------

  pure function normalize( v ) result( vn )
    real(rb), intent(in) :: v(:)
    real(rb)             :: vn(size(v))
    vn = v/sqrt(sum(v*v))
  end function normalize

!---------------------------------------------------------------------------------------------------

  elemental real(rb) function phi( x )
    real(rb), intent(in) :: x
    if (abs(x) > 1E-4_rb ) then
      phi = (one - exp(-x))/x
    else
      phi = one + half*x*(third*x*(one - fourth*x) - one)
    end if
  end function phi

!---------------------------------------------------------------------------------------------------
! Numerical diagonalization of 3x3 matrcies
! Copyright (C) 2006  Joachim Kopp
!---------------------------------------------------------------------------------------------------

  pure subroutine diagonalization( matrix, q, w )
    real(rb), intent(in)  :: matrix(3,3)
    real(rb), intent(out) :: q(3,3), w(3)

    ! ----------------------------------------------------------------------------
    ! Calculates the eigenvalues and normalized eigenvectors of a symmetric 3x3
    ! matrix A using Cardano's method for the eigenvalues and an analytical
    ! method based on vector cross products for the eigenvectors.
    ! Only the diagonal and upper triangular parts of A need to contain
    ! meaningful values. However, all of A may be used as temporary storage
    ! and may hence be destroyed.
    ! ----------------------------------------------------------------------------
    ! Parameters:
    !   matrix: The symmetric input matrix
    !   q: Storage buffer for eigenvectors
    !   w: Storage buffer for eigenvalues
    ! ----------------------------------------------------------------------------

    real(rb), parameter :: eps = epsilon(one)
    real(rb), parameter :: sqrt3 = sqrt(3.0_rb)

    integer  :: i, j
    real(rb) :: a(3,3), norm, n1, n2, w1, w2, w3, thresh, t, wmax8eps
    real(rb) :: m, c1, c0, de, dd, ee, ff, p, sqrtp, r, c, s, phi
    logical  :: success

    a = matrix

    ! Calculate the eigenvalues of a symmetric 3x3 matrix a using Cardano's
    ! analytical algorithm. Only the diagonal and upper triangular parts of A are
    ! accessed. The access is read-only.
    de = a(1,2)*a(2,3)
    dd = a(1,2)**2
    ee = a(2,3)**2
    ff = a(1,3)**2
    m  = a(1,1) + a(2,2) + a(3,3)
    c1 = (a(1,1)*a(2,2) + a(1,1)*a(3,3) + a(2,2)*a(3,3)) - (dd + ee + ff)
    c0 = 27.0_rb*(a(3,3)*dd + a(1,1)*ee + a(2,2)*ff - a(1,1)*a(2,2)*a(3,3) - two*a(1,3)*de)

    p = m*m - 3.0_rb*c1
    r = m*(p - 1.5_rb*c1) - half*c0
    sqrtp = sqrt(abs(p))
    phi = third*atan2(sqrt(abs(6.75_rb*c1*c1*(p - c1) + c0*(r + 0.25_rb*c0))),r)

    c = sqrtp*cos(phi)
    s = (one/sqrt3)*sqrtp*sin(phi)

    p = third*(m - c)
    w1 = p + c
    w2 = p - s
    w3 = p + s

    ! Sort eigenvalues in decreasing order:
    if (abs(w1) < abs(w3)) call swap( w1, w3 )
    if (abs(w1) < abs(w2)) call swap( w1, w2 )
    if (abs(w2) < abs(w3)) call swap( w2, w3 )
    w = [w1, w2, w3]

    wmax8eps = 8.0_rb*eps*abs(w1)
    thresh = wmax8eps**2

    ! Prepare calculation of eigenvectors
    n1 = a(1,2)**2 + a(1,3)**2
    n2 = a(1,2)**2 + a(2,3)**2
    q(1,1) = a(1,2)*a(2,3) - a(1,3)*a(2,2)
    q(1,2) = q(1,1)
    q(2,1) = a(1,3)*a(1,2) - a(2,3)*a(1,1)
    q(2,2) = q(2,1)
    q(3,2) = a(1,2)**2

    ! Calculate first eigenvector by the formula v(1) = (A - lambda(1)).e1 x (A - lambda(1)).e2
    a(1,1) = a(1,1) - w1
    a(2,2) = a(2,2) - w1
    q(:,1) = [q(1,2) + a(1,3)*w1, q(2,2) + a(2,3)*w1, a(1,1)*a(2,2) - q(3,2)]
    call compute_eigenvector( q(:,1), a, n1, n2 )

    ! Prepare calculation of second eigenvector     
    t = w1 - w2

    ! Is this eigenvalue degenerate?
    if (abs(t) > wmax8eps) then

      ! For non-degenerate eigenvalue, calculate second eigenvector by the formula
      !         v[1] = (A - lambda[1]).e1 x (A - lambda[1]).e2
      a(1,1) = a(1,1) + t
      a(2,2) = a(2,2) + t
      q(:,2) = [q(1,2) + a(1,3)*w2, q(2,2) + a(2,3)*w2, a(1,1)*a(2,2) - q(3,2)]
      call compute_eigenvector( q(:,2), a, n1, n2 )

    else

      ! For degenerate eigenvalue, calculate second eigenvector according to
      !         v[1] = v(1) x (A - lambda[1]).e[i]

      ! This would really get too complicated if we could not assume all of A to
      !       contain meaningful values.
      a(2,1) = a(1,2)
      a(3,1) = a(1,3)
      a(3,2) = a(2,3)
      a(1,1) = a(1,1) + w1
      a(2,2) = a(2,2) + w1
      i = 0
      success = .false.
      do while ((i < 3).and.(.not.success))
        i = i + 1
        a(i,i) = a(i,i) - w2
        n1 = sum(a(:,i)**2)
        success = n1 > thresh
        if (success) then
          q(:,2) = cross_product( q(:,1), a(:,i) )
          norm = sum(q(:,2)**2)
          success = norm > (256.0_rb*eps)**2*n1
          if (success) q(:,2) = q(:,2)*sqrt(one/norm)
        end if
      end do

      ! This means that any vector orthogonal to v(1) is an EV.
      if (.not.success) then
        i = 1
        do while (q(i,1) == zero)
          i = i + 1
        end do
        j = 1 + mod(i,3)
        norm = one/sqrt(q(i,1)**2 + q(j,1)**2)
        q(i,2) =  q(j,1)*norm
        q(j,2) = -q(i,1)*norm
        q(1+mod(i+1,3),2) = zero
      end if
    end if

    ! Calculate third eigenvector according to v[2] = v(1) x v[1]
    q(:,3) = cross_product( q(:,1), q(:,2) )

    contains
      !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      pure subroutine compute_eigenvector( q, a, n1tmp, n2tmp )
        real(rb), intent(inout) :: q(3)
        real(rb), intent(in)    :: a(3,3), n1tmp, n2tmp

        real(rb) :: norm, n1, n2, error, t, f

        norm = sum(q**2)
        n1 = n1tmp + a(1,1)**2
        n2 = n2tmp + a(2,2)**2
        error = n1*n2

        ! If the first column is zero, then (1, 0, 0) is an eigenvector
        if (n1 <= thresh) then
          q = [one, zero, zero]

        ! If the second column is zero, then (0, 1, 0) is an eigenvector
        else if (n2 <= thresh) then
          q = [zero, one, zero]

        ! If angle between A(*,1) and A(*,2) is too small, don't use
        !  cross product, but calculate v ~ (1, -A0/A1, 0)
        else if (norm < (64.0_rb*eps)**2*error) then
          t = abs(a(1,2))
          f = -a(1,1)/a(1,2)
          if (abs(a(2,2)) > t) then
            t = abs(a(2,2))
            f = -a(1,2)/a(2,2)
          end if
          if (abs(a(2,3)) > t) f = -a(1,3)/a(2,3)
          norm = one/sqrt(one + f**2)
          q = [norm, f*norm, zero]

        ! This is the standard branch
        else
          q = q*sqrt(one/norm)
        end if

      end subroutine compute_eigenvector
      !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      pure subroutine swap( a, b )
        real(rb), intent(inout) :: a, b
        real(rb) :: c
        c = a; a = b; b = c
      end subroutine swap
      !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  end subroutine diagonalization

!---------------------------------------------------------------------------------------------------

  pure function staircase( x )
    real(rb), intent(in) :: x
    integer              :: staircase
    if (x > zero) then
      staircase = ceiling(x - half)
    else
      staircase = floor(x + half)
    end if
  end function staircase

!---------------------------------------------------------------------------------------------------

  pure function jacobi( u, m ) result( jac )
    real(rb), intent(in)  :: u, m
    real(rb)              :: jac(3)

    integer, parameter :: NN = 16

    integer :: n
    real(rb) :: mu(0:NN-1), nu(0:NN-1), c(0:NN-1), d(0:NN-1)
    real(rb) :: sin_umu, cos_umu, t, r, sn, cn, dn

    if (abs(m) > one) then
      jac = ieee_value(one,ieee_quiet_NaN)
    else if (abs(m) < two*epsilon(one)) then
      jac = [sin(u), cos(u), one]
    else if (abs(m - one) < two*epsilon(one)) then
      jac(1) = tanh(u)
      jac(2:3) = one/cosh(u)
    else
      n = 0
      mu(0) = one
      nu(0) = sqrt(one - m)
      do while (abs(mu(n) - nu(n)) > 4.0_rb*epsilon(one)*abs(mu(n)+nu(n)))
        mu(n+1) = half*(mu(n) + nu(n))
        nu(n+1) = sqrt(mu(n)*nu(n))
        n = n + 1
        if (n >= NN - 1) then
          jac = ieee_value(one,ieee_quiet_NaN)
          return
        end if
      end do
      sin_umu = sin(u * mu(n))
      cos_umu = cos(u * mu(n))
      if (abs(sin_umu) < abs(cos_umu)) then
        t = sin_umu/cos_umu
        c(n) = mu(n) * t
        d(n) = one
        do while (n > 0)
          c(n-1) = d(n)*c(n);
          r = (c(n)*c(n))/mu(n)
          n = n - 1
          d(n) = (r + nu(n))/(r + mu(n))
        end do
        dn = sqrt(one-m)/d(n)
        cn = dn*sign(one,cos_umu)/hypot(one, c(n))
        sn = cn*c(n)/sqrt(one-m)
      else
        t = cos_umu/sin_umu
        c(n) = mu(n) * t
        d(n) = one
        do while (n > 0)
          c(n-1) = d(n)*c(n)
          r = (c(n)*c(n))/mu(n)
          n = n - 1
          d(n) = (r + nu(n))/(r + mu(n))
        end do
        dn = d(n)
        sn = sign(one,sin_umu)/hypot(one, c(n))
        cn = c(n)*sn
      end if
      jac = [sn, cn, dn]
    end if
  end function jacobi

!---------------------------------------------------------------------------------------------------

  pure real(rb) function Carlson_RF( x, y, z ) result( RF )
    real(rb), intent(in) :: x, y, z

    real(rb), parameter :: errtol = 0.001_rb
    real(rb), parameter :: lolim  = (5.0_rb*tiny(1.0_rb))**(1.0_rb/3.0_rb)
    real(rb), parameter :: uplim  = 0.3_rb*(0.2_rb*huge(1.0_rb))**(1.0_rb/3.0_rb)

    real(rb), parameter :: c1 = 1.0_rb/24.0_rb
    real(rb), parameter :: c2 = 3.0_rb/44.0_rb
    real(rb), parameter :: c3 = 1.0_rb/14.0_rb

    real(rb) :: epslon, e2, e3, lamda
    real(rb) :: mu, s, xn, xndev
    real(rb) :: xnroot, yn, yndev, ynroot, zn, zndev, znroot

    if ((min(x,y,z) < zero).or.(max(x,y,z) > uplim).or.(min(x+y,x+z,y+z) < lolim)) then
      RF = ieee_value(one,ieee_quiet_NaN)
    else
      xn = x
      yn = y
      zn = z
      do
        mu = (xn + yn + zn)/3.0_rb
        xndev = 2.0_rb - (mu + xn)/mu
        yndev = 2.0_rb - (mu + yn)/mu
        zndev = 2.0_rb - (mu + zn)/mu
        epslon = max(abs(xndev),abs(yndev),abs(zndev))
        if (epslon < errtol) exit
        xnroot = sqrt(xn)
        ynroot = sqrt(yn)
        znroot = sqrt(zn)
        lamda = xnroot*(ynroot + znroot) + ynroot*znroot
        xn = (xn + lamda)*0.250_rb
        yn = (yn + lamda)*0.250_rb
        zn = (zn + lamda)*0.250_rb
      end do
      e2 = xndev*yndev - zndev*zndev
      e3 = xndev*yndev*zndev
      s = 1.0_rb + (c1*e2 - 0.10_rb - c2*e3)*e2 + c3*e3
      RF = s/sqrt(mu)
    end if
  end function Carlson_RF

!---------------------------------------------------------------------------------------------------

  pure real(rb) function Carlson_RC( x, y ) result( RC )
    real(rb), intent(in) :: x, y

    real(rb), parameter :: errtol = 0.001_rb
    real(rb), parameter :: lolim = 5.0_rb*tiny(one)
    real(rb), parameter :: uplim = 0.2_rb*huge(one)
    real(rb), parameter :: c1 = 1.0_rb/7.0_rb
    real(rb), parameter :: c2 = 9.0_rb/22.0_rb

    real(rb) :: lamda, mu, s, sn, xn, yn

    if ((x < zero).or.(y <= zero).or.(max(x,y) > uplim).or.(x + y < lolim)) then
      RC = ieee_value(one,ieee_quiet_NaN)
    else
      xn = x
      yn = y
      do
        mu = (xn + yn + yn)/3.0_rb
        sn = (yn + mu)/mu - 2.0_rb
        if (abs(sn) < errtol) exit
        lamda = 2.0_rb*sqrt(xn)*sqrt(yn) + yn
        xn = (xn + lamda)*0.250_rb
        yn = (yn + lamda)*0.250_rb
      end do
      s = sn*sn*(0.30_rb + sn*(c1 + sn*(0.3750_rb + sn*c2)))
      RC = (1.0_rb + s)/sqrt(mu)
    end if
  end function Carlson_RC

!---------------------------------------------------------------------------------------------------

  pure real(rb) function Carlson_RJ( x, y, z, p ) result( RJ )
    real(rb), intent(in) :: x, y, z, p

    real(rb), parameter :: errtol = 0.001_rb
    real(rb), parameter :: lolim = (5.0_rb*tiny(one))**(1.0_rb/3.0_rb)
    real(rb), parameter :: uplim = 0.30_rb*(0.2_rb*huge(one))**(1.0_rb/3.0_rb)
    real(rb), parameter :: c1 = 3.0_rb/14.0_rb
    real(rb), parameter :: c2 = 1.0_rb/3.0_rb
    real(rb), parameter :: c3 = 3.0_rb/22.0_rb
    real(rb), parameter :: c4 = 3.0_rb/26.0_rb

    real(rb) :: alfa, beta, ea, eb, ec, e2, e3, epslon
    real(rb) :: lamda, mu, pn, pndev
    real(rb) :: power4, sigma, s1, s2, s3, xn, xndev
    real(rb) :: xnroot, yn, yndev, ynroot, zn, zndev, znroot

    if ((min(x,y,z) < zero).or.(max(x,y,z,p) > uplim).or.(min(x+y,x+z,y+z,p) < lolim)) then
      RJ = ieee_value(one,ieee_quiet_NaN)
    else
      xn = x
      yn = y
      zn = z
      pn = p
      sigma  = 0.0_rb
      power4 = 1.0_rb
      do
        mu = (xn + yn + zn + pn + pn)*0.20_rb
        xndev = (mu - xn)/mu
        yndev = (mu - yn)/mu
        zndev = (mu - zn)/mu
        pndev = (mu - pn)/mu
        epslon = max(abs(xndev),abs(yndev),abs(zndev),abs(pndev))
        if (epslon < errtol) exit
        xnroot = sqrt(xn)
        ynroot = sqrt(yn)
        znroot = sqrt(zn)
        lamda = xnroot*(ynroot + znroot) + ynroot*znroot
        alfa = pn*(xnroot + ynroot + znroot) + xnroot*ynroot*znroot
        alfa = alfa*alfa
        beta = pn*(pn + lamda)*(pn + lamda)
        sigma = sigma + power4*Carlson_RC(alfa,beta)
        power4 = power4*0.250_rb
        xn = (xn + lamda)*0.250_rb
        yn = (yn + lamda)*0.250_rb
        zn = (zn + lamda)*0.250_rb
        pn = (pn + lamda)*0.250_rb
      end do
      ea = xndev*(yndev + zndev) + yndev*zndev
      eb = xndev*yndev*zndev
      ec = pndev*pndev
      e2 = ea - 3.0_rb*ec
      e3 = eb + 2.0_rb*pndev*(ea-ec)
      s1 = 1.0_rb + e2*(-c1 + 0.750_rb*c3*e2 - 1.50_rb*c4*e3)
      s2 = eb*(0.50_rb*c2 + pndev*(-c3-c3 + pndev*c4))
      s3 = pndev*ea*(c2 - pndev*c3) - c2*pndev*ec
      RJ = 3.0_rb*sigma + power4*(s1 + s2 + s3)/(mu*sqrt(mu))
    end if
  end function Carlson_RJ

!---------------------------------------------------------------------------------------------------

  pure function inverse_of_x_plus_ln_x( y ) result( x )
    real(rb), intent(in) :: y
    real(rb)             :: x
    real(rb) :: x0
    real(rb), parameter :: tol = 1.0e-12_rb
    if (y > 0.5671432904097839_rb) then
      x = y - log(y)
    else
      x = exp(y)
    end if
    x0 = x + one
    do while (abs(x - x0) > tol*x0)
      x0 = x
      x = x*(y + one - log(x))/(x + one)
    end do
  end function inverse_of_x_plus_ln_x

!---------------------------------------------------------------------------------------------------

  pure function inv_erfc( y ) result( x )
    real(rb), intent(in) :: y
    real(rb)             :: x
    real(rb), parameter :: tol = 1.0e-8_rb
    real(rb) :: x0, expmx2
    x = sqrt(-log(y))
    x0 = x + one
    do while (abs(x-x0) > tol)
      x0 = x
      expmx2 = exp(-x**2)
      x = x + half*sqrt(Pi)*(uerfc(x,expmx2) - y)/expmx2
    end do
  end function inv_erfc

!---------------------------------------------------------------------------------------------------

  elemental function normSq( x ) result( y )
    complex(rb), intent(in) :: x
    real(rb)                :: y
    y = realpart(x)**2 + imagpart(x)**2
  end function normSq

!---------------------------------------------------------------------------------------------------

  pure function uerf( x, expmx2 ) result( erfx )
    real(rb), intent(in) :: x, expmx2
    real(rb)             :: erfx
    real(rb) :: t
    t = one/(one + p*x)
    erfx = one - t*(a1 + t*(a2 + t*(a3 + t*(a4 + t*a5))))*expmx2
  end function uerf

!---------------------------------------------------------------------------------------------------

  pure function uerfc( x, expmx2 ) result( erfcx )
    real(rb), intent(in) :: x, expmx2
    real(rb)             :: erfcx
    real(rb) :: t
    t = one/(one + p*x)
    erfcx = t*(a1 + t*(a2 + t*(a3 + t*(a4 + t*a5))))*expmx2
  end function uerfc

!---------------------------------------------------------------------------------------------------

  elemental function complex_dot_product( a, b ) result( c )
    complex(rb), intent(in) :: a, b
    real(rb)                :: c
    c = realpart(a)*realpart(b) + imagpart(a)*imagpart(b)
  end function complex_dot_product

!---------------------------------------------------------------------------------------------------

  elemental function complex_cross_product( a, b ) result( c )
    complex(rb), intent(in) :: a, b
    real(rb)                :: c
    c = realpart(a)*imagpart(b) - realpart(b)*imagpart(a)
  end function complex_cross_product

!---------------------------------------------------------------------------------------------------

  elemental function symm1D( i, j ) result( k )
    integer, intent(in) :: i, j
    integer             :: k
    integer :: x, y
    x = min(i,j) - 1
    y = max(i,j) - 1
    k = x + (y + 1)*y/2 + 1
  end function symm1D

!---------------------------------------------------------------------------------------------------

end module math
