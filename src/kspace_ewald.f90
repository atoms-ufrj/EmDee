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
  real(rb),    allocatable :: prefac(:)
  integer,     allocatable :: n(:,:)
  real(rb),    allocatable :: vec(:,:)
  complex(rb), allocatable :: qeikr(:,:)
  complex(rb), allocatable :: S(:,:)
  complex(rb), allocatable :: sigma(:,:)

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

    ! Initialize empty arrays:
    allocate( model%n(0,0), model%prefac(0), model%vec(0,0), &
              model%qeikr(0,0), model%S(0,0), model%sigma(0,0) )

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

!me%alpha = 5.6_rb/30.0_rb ! DELETE THIS LINE
!me%kmax = twoPi/Lsqrt(27.0_rb) ! DELETE THIS LINE

  end subroutine kspace_ewald_set_parameters

!===================================================================================================

  subroutine kspace_ewald_update( me, L )
    class(kspace_ewald), intent(inout) :: me
    real(rb),            intent(in)    :: L(3)

    integer  :: nvecs, n1, n2, n3
    real(rb) :: kmaxSq, unit(3), unitSq(3), B, k1sq, limit, twoPiByV, fourPiByV
    logical  :: nz1, nz12
    integer, allocatable :: n(:,:)

    unit = twoPi/L
    unitSq = unit**2
    me%nmax = ceiling(me%kmax/unit)
    kmaxSq = 1.0001_rb*maxval(unitSq*me%nmax**2)

    allocate( n(3,(me%nmax(1) + 1)*product(2*me%nmax(2:3) + 1)) )

    nvecs = 0
    do n1 = 0, me%nmax(1)
      nz1 = (n1 /= 0)
      k1sq = unitSq(1)*n1*n1
      do n2 = -me%nmax(2), me%nmax(2)
        nz12 = nz1.or.(n2 /= 0)
        limit = kmaxSq - (k1sq + unitSq(2)*n2*n2)
        do n3 = -me%nmax(3), me%nmax(3)
          if ((nz12.or.(n3 /= 0)).and.(unitSq(3)*n3*n3 < limit)) then
            nvecs = nvecs + 1
            n(:,nvecs) = [n1, n2, n3]
          end if
        end do
      end do
    end do

    me%nvecs = nvecs
    me%threadVecs = (nvecs + me%nthreads - 1)/me%nthreads

    if (size(me%prefac) /= nvecs) then
      deallocate( me%prefac, me%n, me%vec, me%qeikr, me%S, me%sigma )
      allocate( me%prefac(nvecs) )
      allocate( me%n(nvecs,3) )
      allocate( me%vec(3,nvecs) )
      allocate( me%qeikr(nvecs,me%natoms) )
      allocate( me%S(nvecs,me%ntypes) )
      allocate( me%sigma(nvecs,me%ntypes) )
    end if

    B = -0.25_rb/me%alpha**2
    twoPiByV = twoPi/product(L)
    fourPiByV = two*twoPiByV
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
        if (n(1,i) == 0) then
          me%prefac(i) = twoPiByV*exp(ksq*B)/ksq
        else
          me%prefac(i) = fourPiByV*exp(ksq*B)/ksq
        end if
        me%vec(:,i) = two*me%prefac(i)*k
      end do
    end block
    !$omp end parallel

  end subroutine kspace_ewald_update

!===================================================================================================

  subroutine kspace_ewald_prepare( me, thread, Rs )
    class(kspace_ewald), intent(inout) :: me
    integer,             intent(in)    :: thread
    real(rb),            intent(in)    :: Rs(:,:)

    integer  :: v1, vN, i, j, itype
    real(rb) :: theta(3)
    complex(rb) :: Z(3), G1(0:me%nmax(1)), G2(-me%nmax(2):me%nmax(2)), &
                   G3(-me%nmax(3):me%nmax(3))

    ! Compute q*exp(i*k*R) for current thread's atoms:
    G1(0) = (one,zero)
    G2(0) = (one,zero)
    G3(0) = (one,zero)
    do j = (thread - 1)*me%threadAtoms + 1, min(thread*me%threadAtoms, me%natoms)
      i = me%index(j)
      theta = twoPi*Rs(:,i)
      Z = cmplx(cos(theta),sin(theta),rb)
      G1(1:me%nmax(1)) = G( Z(1), me%nmax(1) )
      G2(1:me%nmax(2)) = G( Z(2), me%nmax(2) )
      G3(1:me%nmax(3)) = G( Z(3), me%nmax(3) )
      G2(-me%nmax(2):-1) = conjg(G2(me%nmax(2):1:-1))
      G3(-me%nmax(2):-1) = conjg(G3(me%nmax(2):1:-1))
      me%qeikr(:,j) = me%Q(j)*G1(me%n(:,1))*G2(me%n(:,2))*G3(me%n(:,3))
    end do
    !$omp barrier

    ! Compute type-specific structure factors for current thread's vectors:
    v1 = (thread - 1)*me%threadVecs + 1
    vN = min(thread*me%threadVecs, me%nvecs)
    me%S(v1:vN,:) = (zero,zero)
    do j = 1, me%natoms
      itype = me%type(j)
      me%S(v1:vN,itype) = me%S(v1:vN,itype) + me%qeikr(v1:vN,j)
    end do

    contains
      !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      pure function G( Z, n )
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
      end function G
      !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  end subroutine kspace_ewald_prepare

!===================================================================================================

  subroutine kspace_ewald_compute( me, thread, lambda, Potential, Virial, F )
    class(kspace_ewald), intent(inout) :: me
    integer,             intent(in)    :: thread
    real(rb),            intent(in)    :: lambda(:,:)
    real(rb),            intent(inout) :: Potential, Virial
    real(rb), optional,  intent(inout) :: F(:,:)

    integer  :: v1, vN, i, j, itype

    ! Split wave-vectors among threads:
    v1 = (thread - 1)*me%threadVecs + 1
    vN = min(thread*me%threadVecs, me%nvecs)

    ! Compute type-specific structure factors and sigmas for current thread's vectors:
    me%sigma(v1:vN,:) = matmul(me%S(v1:vN,:),lambda)

    ! Compute energy related to current thread's vectors:
    Potential = Potential + sum(me%prefac(v1:vN)*sum(me%S(v1:vn,:).dot.me%sigma(v1:vN,:),2))
    !$omp barrier

    ! Compute forces on current thread's atoms:
    if (present(F)) then
      do j = (thread - 1)*me%threadAtoms + 1, min(thread*me%threadAtoms, me%natoms)
        i = me%index(j)
        itype = me%type(j)
        F(:,i) = F(:,i) + matmul(me%vec, me%sigma(:,itype).cross.me%qeikr(:,j))
      end do
    end if

  end subroutine kspace_ewald_compute

!===================================================================================================

end module kspace_ewald_module
