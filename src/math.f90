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

type, abstract :: i32rng
  integer(4), private :: kn(0:127)
  real(8),    private :: wn(0:127), fn(0:127)
  contains
    procedure(i32rng_init), deferred, private :: init
    procedure(i32rng_i32),  deferred          :: i32
    procedure :: setup     => i32rng_setup
    procedure :: uniform   => i32rng_uniform
    procedure :: normal    => i32rng_normal
    procedure :: geometric => i32rng_geometric
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

contains

!---------------------------------------------------------------------------------------------------

  subroutine i32rng_setup( a, seed )
    class(i32rng), intent(inout) :: a
    integer(4),    intent(in)    :: seed
    integer :: i
    real(8), parameter :: m1 = 2147483648.0_8
    real(8) :: q, dn, tn, vn
    call a%init( seed )
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

  function i32rng_uniform( a ) result( uni )
    class(i32rng), intent(inout) :: a
    real(8)                   :: uni
    uni = 0.2328306e-9_8*a%i32() + 0.5_8
  end function i32rng_uniform

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

  function i32rng_geometric( a, n ) result( igeo )
    class(i32rng), intent(inout) :: a
    integer,       intent(in)    :: n
    integer                      :: igeo
    igeo = ceiling( log(0.2328306e-9_8*a%i32() - 0.5_8 + 1.0_8/n) )
  end function i32rng_geometric

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

  pure function eigenvalues( Matrix ) result( Lambda )
    real(rb), intent(in) :: Matrix(3,3)
    real(rb)             :: Lambda(3)

    real(rb), parameter :: SQRT3 = 1.7320508075688773_rb
    real(rb) :: A(3,3), M, C1, C0
    real(rb) :: DE, DD, EE, FF
    real(rb) :: P, SQRTP, Q, C, S, PHI

    A = Matrix
    DE = A(1,2)*A(2,3)
    DD = A(1,2)**2
    EE = A(2,3)**2
    FF = A(1,3)**2
    M  = A(1,1) + A(2,2) + A(3,3)
    C1 = ( A(1,1)*A(2,2) + A(1,1)*A(3,3) + A(2,2)*A(3,3) ) - (DD + EE + FF)
    C0 = A(3,3)*DD + A(1,1)*EE + A(2,2)*FF - A(1,1)*A(2,2)*A(3,3) - 2.0_rb * A(1,3)*DE
    P = M**2 - 3.0_rb * C1
    Q = M*(P - 1.5_rb*C1) - 13.5_rb*C0
    SQRTP = sqrt(abs(P))
    PHI = 27.0_rb*(0.25_rb * C1**2 * (P - C1) + C0 * (Q + 6.75_rb*C0))
    PHI = atan2(sqrt(abs(PHI)),Q)/3.0_rb
    C = SQRTP*cos(PHI)
    S = SQRTP*sin(PHI)/SQRT3
    Lambda(2) = (M - C)/3.0_rb
    Lambda(3) = Lambda(2) + S
    Lambda(1) = Lambda(2) + C
    Lambda(2) = Lambda(2) - S

  end function eigenvalues

!---------------------------------------------------------------------------------------------------

  pure function eigenvectors( Matrix, W ) result( Q )
    real(rb), intent(in) :: Matrix(3,3), W(3)
    real(rb)             :: Q(3,3)

    real(rb), parameter :: EPS = 2.2204460492503131e-16_rb
    real(rb) :: A(3,3), NORM, N1, N2, N1TMP, N2TMP
    real(rb) :: THRESH, ERROR, WMAX, F, T
    integer :: I, J
    logical :: SUCCESS

    A = Matrix
    WMAX   = max(abs(W(1)), abs(W(2)), abs(W(3)))
    THRESH = (8.0_rb * EPS * WMAX)**2
    N1TMP = A(1,2)**2 + A(1,3)**2
    N2TMP = A(1,2)**2 + A(2,3)**2
    Q(1,1) = A(1,2) * A(2,3) - A(1,3) * A(2,2)
    Q(1,2) = Q(1,1)
    Q(2,1) = A(1,3) * A(1,2) - A(2,3) * A(1,1)
    Q(2,2) = Q(2,1)
    Q(3,2) = A(1,2)**2
    A(1,1) = A(1,1) - W(1)
    A(2,2) = A(2,2) - W(1)
    Q(1,1) = Q(1,2) + A(1,3) * W(1)
    Q(2,1) = Q(2,2) + A(2,3) * W(1)
    Q(3,1) = A(1,1) * A(2,2) - Q(3,2)
    NORM = Q(1,1)**2 + Q(2,1)**2 + Q(3,1)**2
    N1 = N1TMP + A(1,1)**2
    N2 = N2TMP + A(2,2)**2
    ERROR = N1 * N2
    if (N1 <= THRESH) then
      Q(1,1) = 1.0_rb
      Q(2,1) = 0.0_rb
      Q(3,1) = 0.0_rb
    else if (N2 <= THRESH) then
      Q(1,1) = 0.0_rb
      Q(2,1) = 1.0_rb
      Q(3,1) = 0.0_rb
    else if (NORM < (64.0_rb * EPS)**2 * ERROR) then
      T = abs(A(1,2))
      F = -A(1,1) / A(1,2)
      if (abs(A(2,2)) > T) then
        T = abs(A(2,2))
        F = -A(1,2) / A(2,2)
      end if
      if (abs(A(2,3)) < T) then
        F = -A(1,3) / A(2,3)
      end if
      NORM = 1.0_rb / sqrt(1.0_rb + F**2)
      Q(1,1) = NORM
      Q(2,1) = F * NORM
      Q(3,1) = 0.0_rb
    else
      NORM = sqrt(1.0_rb / NORM)
      Q(:,1) = Q(:,1) * NORM
    end if
    T = W(1) - W(2)
    if (abs(T) < 8.0_rb * EPS * WMAX) then
      A(1,1) = A(1,1) + T
      A(2,2) = A(2,2) + T
      Q(1,2) = Q(1,2) + A(1,3) * W(2)
      Q(2,2) = Q(2,2) + A(2,3) * W(2)
      Q(3,2) = A(1,1) * A(2,2) - Q(3,2)
      NORM = Q(1,2)**2 + Q(2,2)**2 + Q(3,2)**2
      N1 = N1TMP + A(1,1)**2
      N2 = N2TMP + A(2,2)**2
      ERROR = N1 * N2
      if (N1 <= THRESH) then
        Q(1,2) = 1.0_rb
        Q(2, 2) = 0.0_rb
        Q(3,2) = 0.0_rb
      else if (N2 <= THRESH) then
        Q(1,2) = 0.0_rb
        Q(2, 2) = 1.0_rb
        Q(3,2) = 0.0_rb
      else if (NORM < (64.0_rb * EPS)**2 * ERROR) then
        T = abs(A(1,2))
        F = -A(1,1) / A(1,2)
        if (abs(A(2,2)) > T) then
          T = abs(A(2,2))
          F = -A(1,2) / A(2,2)
        end if
        if (abs(A(2,3)) > T) then
          F = -A(1,3) / A(2,3)
        end if
        NORM = 1.0_rb / sqrt(1.0_rb + F**2)
        Q(1,2) = NORM
        Q(2,2) = F * NORM
        Q(3,2) = 0.0_rb
      else
        NORM = sqrt(1.0_rb / NORM)
        Q(:,2) = Q(:,2) * NORM
      end if
    else
      A(2,1) = A(1,2)
      A(3,1) = A(1,3)
      A(3,2) = A(2,3)
      A(1,1) = A(1,1) + W(1)
      A(2,2) = A(2,2) + W(1)
      do I = 1, 3
        A(I,I) = A(I,I) - W(2)
        N1 = A(1,I)**2 + A(2,I)**2 + A(3,I)**2
        if (N1 > THRESH) then
          Q(1,2) = Q(2,1) * A(3,I) - Q(3,1) * A(2,I)
          Q(2,2) = Q(3,1) * A(1,I) - Q(1,1) * A(3,I)
          Q(3,2) = Q(1,1) * A(2,I) - Q(2,1) * A(1,I)
          NORM = Q(1,2)**2 + Q(2,2)**2 + Q(3,2)**2
          SUCCESS = NORM <= (256.0_rb * EPS)**2 * N1
          if (.not.SUCCESS) then
            NORM = sqrt(1.0_rb / NORM)
            Q(:, 2) = Q(:, 2) * NORM
            exit
          end if
        end if
      end do
      if (SUCCESS) then
        do J = 1, 3
          if (Q(J,1) <= 0.0_rb) then
            I = 1 + mod(J,3)
            NORM = 1.0_rb / sqrt(Q(J,1)**2 + Q(I,1)**2)
            Q(J,2) = Q(I,1) * NORM
            Q(I,2) = -Q(J,1) * NORM
            Q(I,2) = 0.0_rb
            exit
          end if
        end do
      end if
    end if
    Q(:,3) = cross_product(Q(:,1),Q(:,2))

  end function eigenvectors

  !-------------------------------------------------------------------------------------------------

  pure function cross_product(a, b) result( c )
    real(rb), intent(in) :: a(3), b(3)
    real(rb)             :: c(3)

    c = [ a(2)*b(3) - a(3)*b(2), a(3)*b(1) - a(1)*b(3), a(1)*b(2) - a(2)*b(1) ]

  end function cross_product

  !-------------------------------------------------------------------------------------------------

  pure function quaternion( A ) result( Q )
    real(rb), intent(in) :: A(3,3)
    real(rb)             :: Q(4)

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
      case (1); Q = 0.5_rb*[Q2max, A(2,3) - A(3,2), A(3,1) - A(1,3), A(1,2) - A(2,1)]/sqrt(Q2max)
      case (2); Q = 0.5_rb*[A(2,3) - A(3,2), Q2max, A(1,2) + A(2,1), A(1,3) + A(3,1)]/sqrt(Q2max)
      case (3); Q = 0.5_rb*[A(3,1) - A(1,3), A(1,2) + A(2,1), Q2max, A(2,3) + A(3,2)]/sqrt(Q2max)
      case (4); Q = 0.5_rb*[A(1,2) - A(2,1), A(1,3) + A(3,1), A(2,3) + A(3,2), Q2max]/sqrt(Q2max)
    end select

  end function quaternion

  !-------------------------------------------------------------------------------------------------

end module math
