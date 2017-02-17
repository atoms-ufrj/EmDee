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

module kspace_spme_module

use global
use math
use omp_lib
use nfft
use kspaceModelClass

implicit none

integer(ib), parameter :: nfft_flags = PRE_PHI_HUT + PRE_PSI + MALLOC_X + & !MALLOC_F + MALLOC_F_HAT + &
                                       FFTW_INIT + FFT_OUT_OF_PLACE + NFFT_SORT_NODES ! + &
!                                       NFFT_OMP_BLOCKWISE_ADJOINT

integer(ib), parameter :: fftw_flags = FFTW_ESTIMATE + FFTW_DESTROY_INPUT

!> Abstract class for kspace model spme
!!
!! NOTES: 1) model parameters must be declared individually and tagged with a comment mark "!<>"
!!        2) recognizable parameter types are real(rb) and integer(ib)
!!        3) allocatable one-dimensional arrays (i.e. vectors) are permitted as parameters
!!        4) an integer(ib) scalar parameter - a size - must necessarily succeed every allocatable
!!           parameter or series of equally-sized allocatable parameters.

type, extends(cKspaceModel) :: kspace_spme
  real(rb) :: accuracy !<> Expected accuracy in energy calculations

  integer  :: nmax(3)
  real(rb) :: kmax

  type(nfft_plan), allocatable :: nfft(:)

  complex(rb), pointer :: f(:)
  complex(rb), pointer :: f_hat(:,:)
  complex(rb), allocatable :: S(:,:)

  real(rb),    allocatable :: prefac(:)
  real(rb),    allocatable :: kvec(:,:)

  integer(ib) :: nNodes
  integer(ib) :: threadNodes
  real(rb)    :: Volume

  contains
    procedure :: setup => kspace_spme_setup
    procedure :: set_parameters => kspace_spme_set_parameters
    procedure :: update  => kspace_spme_update
    procedure :: prepare => kspace_spme_prepare
    procedure :: compute => kspace_spme_compute
end type

contains

!===================================================================================================

  subroutine kspace_spme_setup( model, params, iparams )
    class(kspace_spme),  intent(inout) :: model
    real(rb),  optional,  intent(in)    :: params(:)
    integer,   optional,  intent(in)    :: iparams(:)

    ! Model name:
    model%name = "spme"

    ! Model parameters:
    model%accuracy = params(1)

  end subroutine kspace_spme_setup

!===================================================================================================

  subroutine kspace_spme_set_parameters( me, Rc, L, alpha )
    class(kspace_spme), intent(inout) :: me
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

    ! Allocate array to store the charges in complex form:
    allocate( me%f(me%charge%nitems) )
 
    ! Allocate type-wise NFFT plans:
    allocate( me%nfft(me%ntypes) )

    ! Dummy initialization:
    me%nmax = -1
    allocate( me%f_hat(0,0), me%S(0,0), me%prefac(0), me%kvec(0,0) )

  end subroutine kspace_spme_set_parameters

!===================================================================================================

  subroutine kspace_spme_update( me, L )
    class(kspace_spme), intent(inout) :: me
    real(rb),           intent(in)    :: L(3)

    integer  :: nvecs, n1, n2, n3
    real(rb) :: unit(3), unitSq(3), B, k1sq, k12sq, twoPiByV, nv(3), factor
    logical  :: nz1, nz12

    integer(ib) :: i, M_total
    real(rb)    :: ksq
    integer(ib), target :: N(3), nn(3)

    unit = twoPi/L
    N = ceiling(me%kmax/unit)

    ! Update NFFT plans if the number of nodes in any direction changed since last time:
    if (any(N /= me%nmax)) then

      deallocate( me%f_hat, me%S, me%prefac, me%kvec )
      me%nmax = N
      N = 2*N
      nn = 2*[(nfft_next_power_of_2(N(i)),i=1,3)]
      do i = 1, me%ntypes
        associate ( p => me%nfft(i) )
          M_total = me%charge%last(i) - me%charge%first(i) + 1
          call nfft_finalize( p )
          call nfft_init_guru( p, 3, c_loc(N), M_total, c_loc(nn), 6, nfft_flags, fftw_flags )
          if (i == 1) then
            me%nNodes = int(p%N_total,ib)
            allocate( me%f_hat(me%nNodes,me%ntypes) )
          end if
          p%f = c_loc(me%f(me%charge%first(i)))
          p%f_hat = c_loc(me%f_hat(1,i))
        end associate
      end do
      me%threadNodes = (me%nNodes + me%nthreads - 1)/me%nthreads
      allocate( me%S(me%nNodes,me%ntypes), me%prefac(me%nNodes), me%kvec(me%nNodes,3) )

      ! Recompute prefactors and scaled reciprocal-space vectors:
      B = -0.25_rb/me%alpha**2
      twoPiByV = twoPi/product(L)
      unitSq = unit**2
      nvecs = 0
      do n1 = -me%nmax(1), me%nmax(1) - 1
        nz1 = (n1 /= 0)
        k1sq = unitSq(1)*n1*n1
        do n2 = -me%nmax(2), me%nmax(2) - 1
          nz12 = nz1.or.(n2 /= 0)
          k12sq = k1sq + unitSq(2)*n2*n2
          do n3 = -me%nmax(3), me%nmax(3) - 1
            nvecs = nvecs + 1
            nv = [n1,n2,n3]
            ksq = k12sq + unitSq(3)*n3*n3
            if (nz12.or.(n3 /= 0)) then
              me%prefac(nvecs) = twoPiByV*exp(ksq*B)/ksq
            else
              me%prefac(nvecs) = zero
            end if
            me%kvec(nvecs,:) = two*me%prefac(nvecs)*unit*nv
          end do
        end do
      end do

    else

      ! Otherwise, only rescale prefactors and reciprocal-space vectors:
      factor = me%Volume/product(L)
      me%prefac = factor*me%prefac
      me%kvec = factor*me%kvec

    end if

    me%Volume = product(L)

  end subroutine kspace_spme_update

!===================================================================================================

  subroutine kspace_spme_prepare( me, Rs )
    class(kspace_spme), intent(inout) :: me
    real(rb),           intent(in)    :: Rs(:,:)

    integer :: i
    real(rb),    pointer :: x(:,:)

    me%f = me%charge%value
    do i = 1, me%ntypes
      associate ( p => me%nfft(i) )
        call c_f_pointer( p%x, x, [3_nfft_int,p%M_total] )
        x = Rs(:,me%charge%item(me%charge%first(i):me%charge%last(i)))
        x = x - anint(x)
        call nfft_precompute_one_psi( p )
        call nfft_adjoint( p )
      end associate
    end do

    ! Save type-wise structure factors:
    me%S = me%f_hat

  end subroutine kspace_spme_prepare

!===================================================================================================

  subroutine kspace_spme_compute( me, layer, E, F )
    class(kspace_spme), intent(inout) :: me
    integer(ib),        intent(in)    :: layer
    real(rb),           intent(inout) :: E
    real(rb), optional, intent(inout) :: F(:,:)

    integer :: i, j
    real(rb) :: Energy(me%nthreads)
    complex(rb) :: sigma(me%nNodes,me%ntypes)

    !$omp parallel num_threads(me%nthreads)
    call compute_sigma_partial( omp_get_thread_num() + 1 )
    !$omp end parallel

    E = E + sum(Energy)

    if (present(F)) then
      do i = 1, 3
        do j = 1, me%ntypes
          me%f_hat(:,j) = me%kvec(:,i)*sigma(:,j)
          call nfft_trafo( me%nfft(j) )
        end do
        F(i,me%charge%item) = F(i,me%charge%item) - me%charge%value*imagpart(me%f)
      end do
    end if

  contains
    !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    subroutine compute_sigma_partial( thread )
      integer, intent(in) :: thread

      integer :: v1, vN

      v1 = (thread - 1)*me%threadNodes + 1
      vN = min(thread*me%threadNodes, me%nNodes)
      associate( S => me%S(v1:vN,:), sig => sigma(v1:vN,:) )
        sig = matmul(S, me%lambda(:,:,layer))
        Energy(thread) = sum(me%prefac(v1:vN)*sum(dot(S,sig),2))
      end associate

    end subroutine compute_sigma_partial
    !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    elemental function dot( a, b ) result( c )
      complex(rb), intent(in) :: a, b
      real(rb)                :: c
      c = realpart(a)*realpart(b) + imagpart(a)*imagpart(b)
    end function dot
    !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  end subroutine kspace_spme_compute

!===================================================================================================

end module kspace_spme_module
