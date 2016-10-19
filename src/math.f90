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

implicit none

integer(c_int), public, bind(C) :: GSL_PREC_DOUBLE = 0

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

interface

  function gsl_sf_elljac_e( u, m, sn, cn, dn ) bind(C,name="gsl_sf_elljac_e")
    import :: c_double, c_ptr, c_int
    real(c_double), value :: u, m
    type(c_ptr),    value :: sn, cn, dn
    integer(c_int)        :: jacobi_functions
  end function gsl_sf_elljac_e

  function gsl_sf_ellint_RF( x, y, z, mode ) bind(C,name="gsl_sf_ellint_RF")
    import :: c_double, c_int
    real(c_double), value :: x, y, z
    integer(c_int), value :: mode
    real(c_double)        :: gsl_sf_ellint_RF
  end function gsl_sf_ellint_RF

  function gsl_sf_ellint_RJ( x, y, z, p, mode ) bind(C,name="gsl_sf_ellint_RJ")
    import :: c_double, c_int
    real(c_double), value :: x, y, z, p
    integer(c_int), value :: mode
    real(c_double)        :: gsl_sf_ellint_RJ
  end function gsl_sf_ellint_RJ

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
    a%kn(0) = (dn/q)*m1
    a%kn(1) = 0_4
    a%wn(0) = q/m1
    a%wn(127) = dn/m1
    a%fn(0) = 1.0_8
    a%fn(127) = exp( -0.5_8*dn*dn )
    do i = 126, 1, -1
      dn = sqrt( -2.0_8 * log( vn/dn + exp( -0.5_8*dn*dn ) ) )
      a%kn(i+1) = (dn/tn)*m1
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
    real(rb) :: a(3,3), norm, n1, n2, w1, w2, w3, thresh, wmax, t, wmax8eps
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

    wmax8eps = 8.0_rb*eps*wmax
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

end module math
