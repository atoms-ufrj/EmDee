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
!            Federal University of Rio de Janeiro, Brazilmodule lists

module kspace_ewald_module

use global
use math
use omp_lib
use kspaceModelClass

implicit none

!> Abstract class for kspace model ewald
!!
!! NOTES: 1) model parameters must be declared individually and tagged with a comment mark "!<>"
!!        2) recognizable parameter types are real(rb) and integer(ib)
!!        3) allocatable one-dimensional arrays (i.e. vectors) are permitted as parameters
!!        4) an integer(ib) scalar parameter - a size - must necessarily succeed every allocatable
!!           parameter or series of equally-sized allocatable parameters.

type, extends(cKspaceModel) :: kspace_ewald
  real(rb) :: accuracy !<> Expected accuracy in energy calculations

  real(rb) :: kmax
  integer  :: nmax(3)
  integer  :: nvecs
  integer  :: threadVecs
  real(rb) :: Volume
  integer,     allocatable :: n(:,:)
  real(rb),    allocatable :: prefac(:)
  real(rb),    allocatable :: vec(:,:)
  complex(rb), allocatable :: qeikr(:,:)
  complex(rb), allocatable :: S(:,:)

  contains
    procedure :: setup => kspace_ewald_setup
    procedure :: set_parameters => kspace_ewald_set_parameters
    procedure :: update  => kspace_ewald_update
    procedure :: prepare => kspace_ewald_prepare
    procedure :: compute => kspace_ewald_compute
end type

contains

!===================================================================================================

  subroutine kspace_ewald_setup( model, params, iparams )
    class(kspace_ewald),  intent(inout) :: model
    real(rb),  optional,  intent(in)    :: params(:)
    integer,   optional,  intent(in)    :: iparams(:)

    ! Model name:
    model%name = "ewald"

    ! Model parameters:
    model%accuracy = params(1)

    ! Fake initialization:
    allocate( model%n(0,0), model%prefac(0), model%vec(0,0), model%qeikr(0,0), model%S(0,0) )
    model%nmax = -1

  end subroutine kspace_ewald_setup

!===================================================================================================

  subroutine kspace_ewald_set_parameters( me, Rc, L, alpha )
    class(kspace_ewald), intent(inout) :: me
    real(rb),            intent(in)    :: Rc, L(3)
    real(rb),            intent(out)   :: alpha

    real(rb) :: s

    ! For simplicity, parameters alpha and kmax are calculated from the desired accuracy
    ! independently of the box geometry, such as explained in Frenkel and Smit. That is:
    !     accuracy = exp(-s^2)/s^2
    !     alpha = s/Rc
    !     kmax = 2*alpha*s
    s = sqrt(inverse_of_x_plus_ln_x(-log(me%accuracy)))
    alpha = s/Rc
    me%kmax = two*alpha*s

  end subroutine kspace_ewald_set_parameters

!===================================================================================================

  subroutine kspace_ewald_update( me, L )
    class(kspace_ewald), intent(inout) :: me
    real(rb),            intent(in)    :: L(3)

    integer  :: i, j, nvecs, n1, n2, n3, M2, M2M3, maxnvecs
    real(rb) :: kmaxSq, B, fourPiByV, factor
    integer  :: nmax(3), M(3), kvec(3)
    real(rb) :: unit(3), unitSq(3)
    integer, allocatable :: n(:,:)

    unit = twoPi/L
    unitSq = unit**2
    nmax = ceiling(me%kmax/unit)

    if (any(nmax /= me%nmax)) then
      me%nmax = nmax

      M = 2*nmax + 1
      maxnvecs = product(M)/2  ! Product of odds is odd => division by 2 is rounded down
      allocate( n(3,maxnvecs) )

      M2 = M(2)
      M2M3 = M2*M(3)
      nvecs = 0
      kmaxSq = maxval(unitSq*nmax**2)
      do i = maxnvecs+1, 2*maxnvecs
        n1 = i/M2M3
        j = i - n1*M2M3
        n2 = j/M2
        n3 = j - n2*M2
        kvec = [n1,n2,n3] - nmax
        if (sum(unitSq*kvec**2) <= kmaxSq) then
          nvecs = nvecs + 1
          n(:,nvecs) = kvec
        end if
      end do
      me%nvecs = nvecs
      me%threadVecs = (nvecs + me%nthreads - 1)/me%nthreads

      if (size(me%prefac) /= nvecs) then
        deallocate( me%prefac, me%n, me%vec, me%qeikr, me%S )
        allocate( me%prefac(nvecs), me%n(nvecs,3), me%vec(3,nvecs) )
        allocate( me%qeikr(nvecs,me%charge%nitems), me%S(nvecs,me%ntypes) )
      end if

      B = -0.25_rb/me%alpha**2
      me%Volume = product(L)
      fourPiByV = 4.0_rb*Pi/me%Volume
      !$omp parallel num_threads(me%nthreads)
      block
        integer  :: thread, v1, vN, i
        real(rb) :: k(3), ksq

        thread = omp_get_thread_num() + 1
        v1 = (thread - 1)*me%threadVecs + 1
        vN = min(thread*me%threadVecs, nvecs)
        do i = v1, vN
          me%n(i,:) = n(:,i)
          k = unit*n(:,i)
          ksq = sum(k*k)
          me%prefac(i) = fourPiByV*exp(ksq*B)/ksq
          me%vec(:,i) = two*me%prefac(i)*k
        end do
      end block
      !$omp end parallel

    else

      factor = me%Volume/product(L)
      me%prefac = factor*me%prefac
      me%vec = factor*me%vec
      me%Volume = product(L)

    end if

  end subroutine kspace_ewald_update

!===================================================================================================

  subroutine kspace_ewald_prepare( me, Rs )
    class(kspace_ewald), intent(inout) :: me
    real(rb),            intent(in)    :: Rs(:,:)

    !$omp parallel num_threads(me%nthreads)
    block
      integer  :: thread, v1, vN, i, j, itype
      real(rb) :: theta(3)
      complex(rb) :: Z(3), G1(0:me%nmax(1)), G2(-me%nmax(2):me%nmax(2)), G3(-me%nmax(3):me%nmax(3))

      thread = omp_get_thread_num() + 1

      ! Compute q*exp(i*k*R) for current thread's charges:
      G1(0) = (one,zero)
      G2(0) = (one,zero)
      G3(0) = (one,zero)
      do j = (thread - 1)*me%threadCharges + 1, min(thread*me%threadCharges, me%charge%nitems)
        i = me%charge%item(j)
        theta = twoPi*Rs(:,i)
        Z = cmplx(cos(theta),sin(theta),rb)
        G1(1:me%nmax(1)) = recursion( Z(1), me%nmax(1) )
        G2(1:me%nmax(2)) = recursion( Z(2), me%nmax(2) )
        G3(1:me%nmax(3)) = recursion( Z(3), me%nmax(3) )
        G2(-me%nmax(2):-1) = conjg(G2(me%nmax(2):1:-1))
        G3(-me%nmax(2):-1) = conjg(G3(me%nmax(2):1:-1))
        me%qeikr(:,j) = me%charge%value(j)*G1(me%n(:,1))*G2(me%n(:,2))*G3(me%n(:,3))
      end do
      !$omp barrier

      ! Compute type-specific structure factors for current thread's vectors:
      v1 = (thread - 1)*me%threadVecs + 1
      vN = min(thread*me%threadVecs, me%nvecs)
      forall (itype=1:me%ntypes)
        me%S(v1:vN,itype) = sum(me%qeikr(v1:vN,me%charge%first(itype):me%charge%last(itype)),2)
      end forall
    end block
    !$omp end parallel

    contains
      !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      pure function recursion( Z, n ) result( G )
        complex(rb), intent(in) :: Z
        integer,     intent(in) :: n
        complex(rb)             :: G(n)
        integer :: j
        complex(rb) :: B
        B = Z
        j = 1
        do while (j < n)
          G(j) = B
          B = B*Z
          j = j + 1
        end do
        G(j) = B
      end function recursion
      !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  end subroutine kspace_ewald_prepare

!===================================================================================================

  subroutine kspace_ewald_compute( me, layer, E, F )
    class(kspace_ewald), intent(inout) :: me
    integer(ib),         intent(in)    :: layer
    real(rb),            intent(inout) :: E
    real(rb), optional,  intent(inout) :: F(:,:)

    real(rb) :: Energy(me%nthreads)
    complex(rb) :: sigma(me%nvecs,me%ntypes)

    !$omp parallel num_threads(me%nthreads)
    block
      integer :: thread
      thread = omp_get_thread_num() + 1
      call perform_computation( thread, Energy(thread) )
    end block
    !$omp end parallel

    E = E + sum(Energy)

  contains
    !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    subroutine perform_computation( thread, E )
      integer,  intent(in) :: thread
      real(rb), intent(out) :: E

      integer :: v1, vN, a1, aN, i, j, itype

      ! Split wave-vectors among threads:
      v1 = (thread - 1)*me%threadVecs + 1
      vN = min(thread*me%threadVecs, me%nvecs)

      ! Compute type-specific structure factors and sigmas for current thread's vectors:
      sigma(v1:vN,:) = matmul(me%S(v1:vN,:), me%lambda(:,:,layer))

      ! Compute energy related to current thread's vectors:
      E = sum(me%prefac(v1:vN)*sum(dot(me%S(v1:vn,:), sigma(v1:vN,:)),2))
      !$omp barrier

      ! Compute forces on current thread's atoms:
      a1 = (thread - 1)*me%threadCharges + 1
      aN = min(thread*me%threadCharges, me%charge%nitems)
      itype = me%charge%object( a1 )
      do j = a1, aN
        i = me%charge%item(j)
        if (j > me%charge%last(itype)) itype = itype + 1
        F(:,i) = F(:,i) + matmul(me%vec,cross(sigma(:,itype), me%qeikr(:,j)))
      end do

    end subroutine perform_computation
    !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    elemental function dot( a, b ) result( c )
      complex(rb), intent(in) :: a, b
      real(rb)                :: c
      c = realpart(a)*realpart(b) + imagpart(a)*imagpart(b)
    end function dot
    !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    elemental function cross( a, b ) result( c )
      complex(rb), intent(in) :: a, b
      real(rb)                :: c
      c = realpart(a)*imagpart(b) - realpart(b)*imagpart(a)
    end function cross
    !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  end subroutine kspace_ewald_compute

!===================================================================================================

end module kspace_ewald_module
