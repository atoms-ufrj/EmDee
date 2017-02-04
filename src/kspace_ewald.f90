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

  real(rb) :: unit(3)
  integer  :: nmax(3)
  real(rb) :: kmaxSq

  integer  :: nvecs
  real(rb), allocatable :: prefac(:)
  integer,  allocatable :: n(:,:)
  real(rb), allocatable :: vec(:,:)

  integer :: vecs_per_thread

  contains
    procedure :: setup => kspace_ewald_setup
    procedure :: set_parameters => kspace_ewald_set_parameters
    procedure :: update => kspace_ewald_update
    procedure :: compute => kspace_ewald_compute
end type

contains

!===================================================================================================

  subroutine kspace_ewald_setup( model, params, iparams )
    class(kspace_ewald),  intent(inout) :: model
    real(rb), optional,   intent(in)    :: params(:)
    integer,  optional,   intent(in)    :: iparams(:)

    ! Model name:
    model%name = "ewald"

    ! Model parameters:
    model%accuracy = params(1)

    ! Initialize empty arrays:
    allocate( model%n(0,0), model%prefac(0), model%vec(0,0) )

  end subroutine kspace_ewald_setup

!===================================================================================================

  subroutine kspace_ewald_set_parameters( model, Rc, L, Q )
    class(kspace_ewald), intent(inout) :: model
    real(rb),            intent(in)    :: Rc, L(3), Q(:)

    real(rb) :: s

    ! For simplicity, parameters alpha and kmax are calculated from the desired accuracy
    ! independently of the box geometry, such as explained in Frenkel and Smit. That is:
    !     accuracy = exp(-s^2)/s^2
    !     alpha = s/Rc
    !     kmax = 2*alpha*s
    s = sqrt(inverse_of_x_plus_ln_x(-log(model%accuracy)))
    model%alpha = s/Rc
    model%kmax = two*model%alpha*s

!model%alpha = 5.6_rb/30.0_rb ! DELETE THIS LINE
!model%kmax = sqrt(27.0_rb) ! DELETE THIS LINE

  end subroutine kspace_ewald_set_parameters

!===================================================================================================

  subroutine kspace_ewald_update( model, L )
    class(kspace_ewald), intent(inout) :: model
    real(rb),            intent(in)    :: L(3)

    integer  :: maxnvecs, nvecs, n1, n2, n3, n(3)
    real(rb) :: unitSq(3), B, factor, k1sq, k12sq, ksq, twoPiByV, prefac

    model%unit = twoPi/L
    unitSq = model%unit**2
    model%nmax = ceiling(model%kmax/model%unit)
    model%kmaxSq = 1.0001_rb*maxval(unitSq*model%nmax**2)

    maxnvecs = (model%nmax(1) + 1)*product(2*model%nmax(2:3) + 1)
    if (size(model%prefac) < maxnvecs) then
      deallocate( model%prefac, model%n, model%vec )
      allocate( model%prefac(maxnvecs), model%n(maxnvecs,3), model%vec(3,maxnvecs) )
    end if

    B = -0.25_rb/model%alpha**2
    twoPiByV = twoPi/product(L)
    nvecs = 0
    do n1 = 0, model%nmax(1)
      factor = merge(1,2,n1 == 0)*twoPiByV
      k1sq = unitSq(1)*n1*n1
      do n2 = -model%nmax(2), model%nmax(2)
        k12sq = k1sq + unitSq(2)*n2*n2
        do n3 = -model%nmax(3), model%nmax(3)
          if ((n1 /= 0).or.(n2 /= 0).or.(n3 /= 0)) then
            ksq = k12sq + unitSq(3)*n3*n3
            if (ksq < model%kmaxSq) then
              nvecs = nvecs + 1
              n = [n1, n2, n3]
              prefac = factor*exp(ksq*B)/ksq
              model%n(nvecs,:) = n
              model%prefac(nvecs) = prefac
              model%vec(:,nvecs) = two*prefac*model%unit*n
            end if
          end if
        end do
      end do
    end do

    model%nvecs = nvecs
    model%vecs_per_thread = (nvecs + model%nthreads - 1)/model%nthreads

  end subroutine kspace_ewald_update

!===================================================================================================

  subroutine kspace_ewald_compute( model, R, F, Potential, Virial )
    class(kspace_ewald), intent(in)    :: model
    real(rb),            intent(in)    :: R(:,:)
    real(rb),            intent(inout) :: F(:,:), Potential, Virial

    real(rb) :: E = zero
    complex(rb) :: S(model%nvecs)
    complex(rb), allocatable :: qeikr(:,:)

    allocate( qeikr(model%nvecs,model%natoms) )

    associate( nmax => model%nmax, n => model%n(1:model%nvecs,:) )
      !$omp parallel num_threads(model%nthreads) reduction(+:E)
      block
        integer :: thread, first, last, v1, vN, i, j, k
        real(rb) :: theta(3)
        complex(rb) :: Z(3)
        complex(rb) :: G1(0:nmax(1)), G2(-nmax(2):nmax(2)), G3(-nmax(3):nmax(3))

        ! Split atoms and wave-vectors among threads:
        thread = omp_get_thread_num() + 1
        first = (thread - 1)*model%atoms_per_thread + 1
        last = min(thread*model%atoms_per_thread, model%natoms)
        v1 = (thread - 1)*model%vecs_per_thread + 1
        vN = min(thread*model%vecs_per_thread, model%nvecs)

        ! Compute qj*exp(i*k*Rj) for local atoms:
        G1(0) = (one,zero)
        G2(0) = (one,zero)
        G3(0) = (one,zero)
        do i = first, last
          j = model%index(i)
          theta = model%unit*R(:,j)
          Z = cmplx(cos(theta),sin(theta),rb)
          G1(1:nmax(1)) = recursion( Z(1), nmax(1) )
          G2(1:nmax(2)) = recursion( Z(2), nmax(2) )
          G3(1:nmax(3)) = recursion( Z(3), nmax(3) )
          G2(-nmax(2):-1) = conjg(G2(nmax(2):1:-1))
          G3(-nmax(2):-1) = conjg(G3(nmax(2):1:-1))
          qeikr(:,i) = model%Q(j)*G1(n(:,1))*G2(n(:,2))*G3(n(:,3))
        end do
        !$omp barrier

        ! Compute structure factors for local vectors:
        S(v1:vN) = sum(qeikr(v1:vN,:),2)
        !$omp barrier

        ! Compute forces of local atoms:
        do j = first, last
          i = model%index(j)
          do k = 1, model%nvecs
            F(:,i) = F(:,i) + ( S(k) .op. qeikr(k,j) )*model%vec(:,k)
          end do
        end do

        ! Compute part of energy related to local vectors:
        E = sum(model%prefac(v1:vN)*normSq(S(v1:vn)))
      end block
      !$omp end parallel
    end associate

    Potential = Potential + E
    Virial = Virial + E

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
  end subroutine kspace_ewald_compute

!===================================================================================================

end module kspace_ewald_module
