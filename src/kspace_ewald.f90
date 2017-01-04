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

! TODO: Enter a wave-vector length cutoff and use it to obtain kxmax, kymax, and kzmax. Then
!       split G in q_exp_ik_dot_r into Gx, Gy, and Gz.

module kspace_ewald_module

use global
use omp_lib
use kspaceModelClass

implicit none

!> Abstract class for kspace model ewald
!! NOTE: model parameters must be declared individually and tagged with comment mark "!<>"
type, extends(cKspaceModel) :: kspace_ewald
  real(rb) :: alpha  !<> Damping parameter
  integer(ib) :: kmax   !<> Wave-vector length cutoff

  integer :: nvectors
  integer,  allocatable :: k(:,:)
  real(rb), allocatable :: ksq(:)
  real(rb), allocatable :: prefac(:)
  contains
    procedure :: setup => kspace_ewald_setup
    procedure :: initialize => kspace_ewald_initialize
    procedure :: compute => kspace_ewald_compute
end type

contains

!---------------------------------------------------------------------------------------------------

  subroutine kspace_ewald_setup( model, params, iparams )
    class(kspace_ewald),  intent(inout) :: model
    real(rb), optional,   intent(in)    :: params(:)
    integer,  optional,   intent(in)    :: iparams(:)

    ! Model name:
    model%name = "ewald"

    ! Model parameters:
    model%alpha = params(1)
    model%kmax = iparams(1)

  end subroutine kspace_ewald_setup

!---------------------------------------------------------------------------------------------------

  subroutine kspace_ewald_initialize( model, nthreads, Rc, invL, q )
    class(kspace_ewald), intent(inout) :: model
    integer,             intent(in)    :: nthreads
    real(rb),            intent(in)    :: Rc, invL(3), q(:)

    integer  :: kmax, ksqmax, n, ksq, kx, ky, kz
    real(rb) :: B, factor

    kmax = model%kmax
    ksqmax = kmax*kmax + 2
    n = 0
    do kx = 0, kmax
      do ky = -kmax, kmax
        do kz = -kmax, kmax
          ksq = kx*kx + ky*ky + kz*kz
          if ((ksq < ksqmax).and.(ksq /= 0)) n = n + 1
        end do
      end do
    end do

    model%nvectors = n
    model%nthreads = nthreads
    allocate( model%k(n,3), model%ksq(n), model%prefac(n) )

    B = -0.25_rb/model%alpha
    n = 0
    do kx = 0, kmax
      factor = merge(4.0_rb,8.0_rb,kx == 0)*pi
      do ky = -kmax, kmax
        do kz = -kmax, kmax
          ksq = kx*kx + ky*ky + kz*kz
          if ((ksq < ksqmax).and.(ksq /= 0)) then
            n = n + 1
            model%k(n,:) = [kx, ky, kz]
            model%ksq(n) = twoPi**2*kx*kx + ky*ky + kz*kz

            model%prefac(n) = factor*exp(B*ksq)/ksq
print*,             model%prefac(n)

          end if
        end do
      end do
    end do
stop
  end subroutine kspace_ewald_initialize

!---------------------------------------------------------------------------------------------------

  subroutine kspace_ewald_compute( me, nthreads, q, Rs, V, E )
    class(kspace_ewald), intent(in)  :: me
    integer,       intent(in)  :: nthreads
    real(rb),      intent(in)  :: q(:), Rs(3,size(q)), V
    real(rb),      intent(out) :: E

    integer :: natoms
    complex(rb) :: rho(me%nvectors)

    natoms = size(q)
print*, natoms
print*, q
print*, Rs
stop
    !$omp parallel num_threads(nthreads) reduction(+:rho)
    block
      integer :: thread, per_thread, first, last
      thread = omp_get_thread_num() + 1
      per_thread = natoms/nthreads
      first = (thread - 1)*per_thread
      last = min(thread*per_thread, natoms)
      rho = kspace_charge_density( me, Q(first:last), Rs(:,first:last) )
      print*, rho
    end block
    !$omp end parallel
    E = sum(me%prefac*(realpart(rho)**2 + imagpart(rho)**2))

    contains
      !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      pure function q_exp_ik_dot_r( me, Q, Rs ) result( F )
        class(kspace_ewald), intent(in) :: me
        real(rb),      intent(in) :: Q, Rs(3)
        complex(rb)               :: F(me%nvectors)
    
        integer  :: j
        real(rb) :: twoPiRs(3)
        complex(rb) :: G(3,-me%kmax:me%kmax)

        twoPiRs = twoPi*Rs
        G(:,0) = (one,zero)
        G(:,1) = cmplx( cos(twoPiRs), sin(twoPiRs), kind = rb )
        do j = 2, me%kmax
          G(:,j) = G(:,j-1)*G(:,1)
        end do
        forall (j=1:me%kmax) G(2:3,-j) = conjg(G(2:3,j))
        F = Q*G(1,me%k(:,1))*G(2,me%k(:,2))*G(3,me%k(:,3))

      end function q_exp_ik_dot_r
      !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      pure function kspace_charge_density( me, Q, Rs ) result( rho )
        class(kspace_ewald), intent(in) :: me
        real(rb),      intent(in) :: Q(:)
        real(rb),      intent(in) :: Rs(3,size(Q))
        complex(rb)               :: rho(me%nvectors)

        integer :: i

        rho = (zero,zero)
        do i = 1, size(Q)
          rho = rho + q_exp_ik_dot_r( me, Q(i), Rs(:,i) )
        end do

      end function kspace_charge_density
      !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  end subroutine kspace_ewald_compute

!---------------------------------------------------------------------------------------------------

end module kspace_ewald_module
