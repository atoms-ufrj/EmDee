!   This file is part of EmDee.
!
!    EmDee is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Piblic License as piblished by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    EmDee is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Piblic License for more details.
!
!    You should have received a copy of the GNU General Piblic License
!    along with EmDee. If not, see <http://www.gnu.org/licenses/>.
!
!    Author: Charlles R. A. Abreu (abreu@eq.ufrj.br)
!            Applied Thermodynamics and Molecular Simulation
!            Federal University of Rio de Janeiro, Brazil

module nfft

use global

implicit none

integer, parameter :: nfft_int = c_long_long

integer(ib), parameter ::                PRE_PHI_HUT = ishft(1_ib, 0)
integer(ib), parameter ::                     FG_PSI = ishft(1_ib, 1)
integer(ib), parameter ::                PRE_LIN_PSI = ishft(1_ib, 2)
integer(ib), parameter ::                 PRE_FG_PSI = ishft(1_ib, 3)
integer(ib), parameter ::                    PRE_PSI = ishft(1_ib, 4)
integer(ib), parameter ::               PRE_FULL_PSI = ishft(1_ib, 5)
integer(ib), parameter ::                   MALLOC_X = ishft(1_ib, 6)
integer(ib), parameter ::               MALLOC_F_HAT = ishft(1_ib, 7)
integer(ib), parameter ::                   MALLOC_F = ishft(1_ib, 8)
integer(ib), parameter ::           FFT_OUT_OF_PLACE = ishft(1_ib, 9)
integer(ib), parameter ::                  FFTW_INIT = ishft(1_ib,10)
integer(ib), parameter ::            NFFT_SORT_NODES = ishft(1_ib,11)
integer(ib), parameter :: NFFT_OMP_BLOCKWISE_ADJOINT = ishft(1_ib,12)
integer(ib), parameter ::                PRE_ONE_PSI = PRE_LIN_PSI + &
                                                       PRE_FG_PSI  + &
                                                       PRE_PSI     + &
                                                       PRE_FULL_PSI
integer(ib), parameter ::               FFTW_MEASURE = ishft(0_ib, 0)
integer(ib), parameter ::         FFTW_DESTROY_INPUT = ishft(1_ib, 0)
integer(ib), parameter ::             FFTW_UNALIGNED = ishft(1_ib, 1)
integer(ib), parameter ::       FFTW_CONSERVE_MEMORY = ishft(1_ib, 2)
integer(ib), parameter ::            FFTW_EXHAUSTIVE = ishft(1_ib, 3)
integer(ib), parameter ::        FFTW_PRESERVE_INPUT = ishft(1_ib, 4)
integer(ib), parameter ::               FFTW_PATIENT = ishft(1_ib, 5)
integer(ib), parameter ::              FFTW_ESTIMATE = ishft(1_ib, 6)
integer(ib), parameter ::           FFTW_WISDOM_ONLY = ishft(1_ib,21)

type, bind(C) :: nfft_plan
  integer(nfft_int) :: N_total           ! Total number of Fourier coefficients.
  integer(nfft_int) :: M_total           ! Total number of samples.
  type(c_ptr)       :: f_hat             ! Vector of Fourier coefficients, size is N_total float_types.
  type(c_ptr)       :: f                 ! Vector of samples, size is M_total float types.
  type(c_funptr)    :: mv_trafo          ! Adjoint transform
  type(c_funptr)    :: mv_adjoint        ! Transform
  integer(nfft_int) :: d                 ! Dimension (rank)
  type(c_ptr)       :: N                 ! Multi-bandwidth
  type(c_ptr)       :: sigma             ! Oversampling factor
  type(c_ptr)       :: nn                ! Length of FFTW transforms. This is equal to sigma*N. The default is to use a power of two that satifies 2 <= sigma < 4
  integer(nfft_int) :: nn_total          ! Combined total length of FFTW transforms
  integer(nfft_int) :: m                 ! Cut-off parameter for window function. Default values for the different window functions are 6 (KAISER_BESSEL), 9 (SINC_POWER), 11 (B_SPLINE), 12 (GAUSSIAN)
  type(c_ptr)       :: b                 ! Shape parameter of the window function.
  integer(nfft_int) :: K                 ! Number of equispaced samples of window function. Used for flag PRE_LIN_PSI
  integer(ib)       :: nfft_flags        ! Flags for precomputation, (de)allocation, and FFTW usage, default setting is PRE_PHI_HUT | PRE_PSI | MALLOC_X | MALLOC_F_HAT | MALLOC_F | FFTW_INIT | FFT_OUT_OF_PLACE
  integer(ib)       :: fftw_flags        ! Flags for the FFTW, default is FFTW_ESTIMATE | FFTW_DESTROY_INPUT
  type(c_ptr)       :: x                 ! Nodes in time/spatial domain, size is d*M
  real(rb)          :: MEASURE_TIME_t(3) ! Measured time for each step if MEASURE_TIME is set
  type(c_ptr)       :: my_fftw_plan1     ! Forward FFTW plan
  type(c_ptr)       :: my_fftw_plan2     ! Backward FFTW plan
  type(c_ptr)       :: c_phi_inv         ! Precomputed data for the diagonal matrix D, size is N_0 ... N_{d-1} doubles
  type(c_ptr)       :: psi               ! Precomputed data for the sparse matrix B, size depends on precomputation scheme
  type(c_ptr)       :: psi_index_g       ! Indices in source/target vector for PRE_FULL_PSI
  type(c_ptr)       :: psi_index_f       ! Indices in source/target vector for PRE_FULL_PSI
  type(c_ptr)       :: g                 ! Oversampled vector of samples, size is n_total double complex
  type(c_ptr)       :: g_hat             ! Zero-padded vector of Fourier coefficients, size is n_total fftw_complex
  type(c_ptr)       :: g1                ! Input of fftw
  type(c_ptr)       :: g2                ! Output of fftw
  type(c_ptr)       :: spline_coeffs     ! Input for de Boor algorithm if B_SPLINE or SINC_POWER is defined
  type(c_ptr)       :: index_x           ! Index array for nodes x used when flag NFFT_SORT_NODES is set
end type

interface

  integer(ib) function nfft_next_power_of_2( x ) bind(C,name="nfft_next_power_of_2")
    import :: ib
    integer(ib), value :: x
  end function nfft_next_power_of_2

  subroutine nfft_init_3d( ths, N1, N2, N3, M ) bind(C,name="nfft_init_3d")
    import :: nfft_plan, ib
    type(nfft_plan), intent(inout) :: ths
    integer(ib),     value         :: N1, N2, N3, M
  end subroutine nfft_init_3d

  subroutine nfft_init_guru( ths, d, N, M_total, nn, m, flags, fftw_flags ) bind(C,name="nfft_init_guru")
    import :: nfft_plan, c_ptr, ib
    type(nfft_plan), intent(inout) :: ths
    integer(ib),     value         :: d, M_total, m, flags, fftw_flags
    type(c_ptr),     value         :: N, nn
  end subroutine nfft_init_guru

  subroutine nfft_precompute_one_psi( ths_plan ) bind(C,name="nfft_precompute_one_psi")
    import :: nfft_plan
    type(nfft_plan), intent(inout) :: ths_plan
  end subroutine nfft_precompute_one_psi

  subroutine nfft_adjoint( ths_plan ) bind(C,name="nfft_adjoint")
    import :: nfft_plan
    type(nfft_plan), intent(inout) :: ths_plan
  end subroutine nfft_adjoint

  subroutine nfft_trafo( ths_plan ) bind(C,name="nfft_trafo")
    import :: nfft_plan
    type(nfft_plan), intent(inout) :: ths_plan
  end subroutine nfft_trafo

  subroutine nfft_finalize( ths_plan ) bind(C,name="nfft_finalize")
    import :: nfft_plan
    type(nfft_plan), intent(inout) :: ths_plan
  end subroutine nfft_finalize

end interface

contains

  pure function match( a, b ) result( c )
    integer(ib), intent(in) :: a, b
    logical                 :: c
    c = iand(a,b) /= 0_ib
  end function match

end module nfft
