!   This file is part of EmDee.
!
!    EmDee is free software: you can redistrc_intute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    EmDee is distrc_intuted in the hope that it will be useful,
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

module EmDee

use iso_c_binding

implicit none

integer, parameter, private :: ib = c_int, rb = c_double, lb = c_bool

type, bind(C) :: tOpts
  logical(lb) :: Translate            ! Flag to activate/deactivate translations
  logical(lb) :: Rotate               ! Flag to activate/deactivate rotations
  integer(ib) :: RotationMode         ! Algorithm used for free rotation of rigid bodies
  logical(lb) :: AutoForceCompute     ! Flag to activate/deactivate automatic force computations
  logical(lb) :: AutoBodyUpdate       ! Flag to activate/deactivate automatic rigid body update
end type tOpts

type, bind(C) :: tEnergy
  real(rb)    :: Potential            ! Total potential energy of the system
  real(rb)    :: Dispersion           ! Dispersion (vdW) part of the potential energy
  real(rb)    :: Coulomb              ! Electrostatic part of the potential energy
  real(rb)    :: Bond
  real(rb)    :: Angle
  real(rb)    :: Dihedral
  real(rb)    :: Kinetic              ! Total kinetic energy of the system
  real(rb)    :: TransPart(3)         ! Translational kinetic energy at each dimension
  real(rb)    :: Rotational           ! Total rotational kinetic energy of the system
  real(rb)    :: RotPart(3)           ! Rotational kinetic energy around each principal axis
  real(rb)    :: ShadowPotential
  real(rb)    :: ShadowKinetic
  real(rb)    :: ShadowRotational
  logical(lb) :: Compute              ! Flag to activate/deactivate energy computations
  logical(lb) :: UpToDate             ! Flag to attest whether energies have been computed
end type tEnergy

type, bind(C), public :: tTime
  real(rb) :: Pair                    ! Time taken in force calculations
  real(rb) :: Motion
  real(rb) :: Neighbor
  real(rb) :: Total                   ! Total time since initialization
end type tTime

type, bind(C) :: tEmDee
  integer(ib)   :: Builds             ! Number of neighbor list builds
  type(tTime)   :: Time
  type(tEnergy) :: Energy             ! All energy terms
  real(rb)      :: Virial             ! Total internal virial of the system
  real(rb)      :: BodyVirial         ! Rigid body contribution to the internal virial
  integer(ib)   :: DoF                ! Total number of degrees of freedom
  integer(ib)   :: RotDoF             ! Number of rotational degrees of freedom
  type(c_ptr)   :: Data               ! Pointer to system data
  type(tOpts)   :: Options            ! List of options to change EmDee's behavior
end type tEmDee

interface

  function EmDee_system( threads, layers, rc, skin, N, types, masses, bodies ) &
    bind(C,name="EmDee_system")
    import :: c_int, c_double, c_ptr, tEmDee
    integer(c_int), value :: threads, layers, N
    real(c_double), value :: rc, skin
    type(c_ptr),    value :: types, masses, bodies
    type(tEmDee)          :: EmDee_system
  end function EmDee_system

  subroutine EmDee_share_phase_space( mdkeep, mdlose ) &
    bind(C,name="EmDee_share_phase_space")
    import :: tEmDee
    type(tEmDee), value         :: mdkeep
    type(tEmDee), intent(inout) :: mdlose
  end subroutine EmDee_share_phase_space

  subroutine EmDee_layer_based_parameters( md, Rc, Bonded ) &
    bind(C,name="EmDee_layer_based_parameters")
    import :: c_double, c_int, tEmDee
    type(tEmDee),   value      :: md
    real(c_double), intent(in) :: Rc(*)
    integer(c_int), intent(in) :: Bonded(*)
  end subroutine EmDee_layer_based_parameters

  subroutine EmDee_switch_model_layer( md, layer ) &
    bind(C,name="EmDee_switch_model_layer")
    import :: c_int, tEmDee
    type(tEmDee),   value :: md
    integer(c_int), value :: layer
  end subroutine EmDee_switch_model_layer

  subroutine EmDee_set_pair_model( md, itype, jtype, model, kCoul ) &
    bind(C,name="EmDee_set_pair_model")
    import :: c_int, c_ptr, c_double, tEmDee
    type(tEmDee),   value :: md
    integer(c_int), value :: itype, jtype
    type(c_ptr),    value :: model
    real(c_double), value :: kCoul
  end subroutine EmDee_set_pair_model

  subroutine EmDee_set_pair_multimodel( md, itype, jtype, model, kCoul ) &
    bind(C,name="EmDee_set_pair_multimodel")
    import :: c_int, c_ptr, c_double, tEmDee
    type(tEmDee),   value :: md
    integer(c_int), value :: itype, jtype
    type(c_ptr),    intent(in) :: model(*)
    real(c_double), intent(in) :: kCoul(*)
  end subroutine EmDee_set_pair_multimodel

  subroutine EmDee_set_kspace_model( md, model ) &
    bind(C,name="EmDee_set_kspace_model")
    import :: c_ptr, tEmDee
    type(tEmDee), value :: md
    type(c_ptr),  value :: model
  end subroutine EmDee_set_kspace_model

  subroutine EmDee_set_coul_model( md, model ) &
    bind(C,name="EmDee_set_coul_model")
    import :: c_ptr, tEmDee
    type(tEmDee), value :: md
    type(c_ptr),  value :: model
  end subroutine EmDee_set_coul_model

  subroutine EmDee_set_coul_multimodel( md, model ) &
    bind(C,name="EmDee_set_coul_multimodel")
    import :: c_ptr, tEmDee
    type(tEmDee), value      :: md
    type(c_ptr),  intent(in) :: model(*)
  end subroutine EmDee_set_coul_multimodel

  subroutine EmDee_ignore_pair( md, i, j ) &
    bind(C,name="EmDee_ignore_pair")
    import :: c_int, c_ptr, tEmDee
    type(tEmDee),   value :: md
    integer(c_int), value :: i, j
  end subroutine EmDee_ignore_pair

  subroutine EmDee_add_bond( md, i, j, model ) &
    bind(C,name="EmDee_add_bond")
    import :: c_int, c_ptr, tEmDee
    type(tEmDee),   value :: md
    integer(c_int), value :: i, j
    type(c_ptr),    value :: model
  end subroutine EmDee_add_bond

  subroutine EmDee_add_angle( md, i, j, k, model ) &
    bind(C,name="EmDee_add_angle")
    import :: c_int, c_ptr, tEmDee
    type(tEmDee),   value :: md
    integer(c_int), value :: i, j, k
    type(c_ptr),    value :: model
  end subroutine EmDee_add_angle

  subroutine EmDee_add_dihedral( md, i, j, k, l, model ) &
    bind(C,name="EmDee_add_dihedral")
    import :: c_int, c_ptr, tEmDee
    type(tEmDee),   value :: md
    integer(c_int), value :: i, j, k, l
    type(c_ptr),    value :: model
  end subroutine EmDee_add_dihedral

  subroutine EmDee_upload( md, item, address ) &
    bind(C,name="EmDee_upload")
    import :: c_ptr, c_char, tEmDee
    type(tEmDee),      intent(inout) :: md
    character(c_char), intent(in)    :: item(*)
    type(c_ptr),       value         :: address
  end subroutine EmDee_upload

  subroutine EmDee_download( md, item, address ) &
    bind(C,name="EmDee_download")
    import :: c_ptr, c_char, tEmDee
    type(tEmDee),      value      :: md
    character(c_char), intent(in) :: item(*)
    type(c_ptr),       value      :: address
  end subroutine EmDee_download

  subroutine EmDee_random_momenta( md, kT, adjust, seed ) &
    bind(C,name="EmDee_random_momenta")
    import :: c_int, c_double, c_ptr, c_bool, tEmDee
    type(tEmDee),    intent(inout) :: md
    real(c_double),  value         :: kT
    logical(c_bool), value         :: adjust
    integer(c_int),  value         :: seed
  end subroutine EmDee_random_momenta

  subroutine EmDee_boost( md, lambda, alpha, dt ) &
    bind(C,name="EmDee_boost")
    import :: c_double, c_ptr, tEmDee
    type(tEmDee),   intent(inout) :: md
    real(c_double), value         :: lambda, alpha, dt
  end subroutine EmDee_boost

  subroutine EmDee_displace( md, lambda, alpha, dt ) &
    bind(C,name="EmDee_displace")
    import :: c_double, c_ptr, tEmDee
    type(tEmDee),   intent(inout) :: md
    real(c_double), value         :: lambda, alpha, dt
  end subroutine EmDee_displace

  subroutine EmDee_verlet_step( md, dt ) &
    bind(C,name="EmDee_verlet_step")
    import :: tEmDee, rb
    type(tEmDee), intent(inout) :: md
    real(rb),     value         :: dt
  end subroutine EmDee_verlet_step

  subroutine EmDee_compute_forces( md ) &
    bind(C,name="EmDee_compute_forces")
    import :: tEmDee
    type(tEmDee), intent(inout) :: md
  end subroutine EmDee_compute_forces

  subroutine EmDee_rdf( md, bins, pairs, itype, jtype, g ) &
    bind(C,name="EmDee_rdf")
    import :: tEmDee, c_int, c_double
    type(tEmDee),   value       :: md
    integer(c_int), value       :: pairs, bins
    integer(c_int), intent(in)  :: itype(pairs), jtype(pairs)
    real(c_double), intent(out) :: g(bins,pairs)
  end subroutine EmDee_rdf

  ! MODEL MODIFIERS:
  type(c_ptr) function EmDee_shifted( model ) &
    bind(C,name="EmDee_shifted")
    import :: c_ptr
    type(c_ptr), value :: model
  end function EmDee_shifted

  type(c_ptr) function EmDee_shifted_force( model ) &
    bind(C,name="EmDee_shifted_force")
    import :: c_ptr
    type(c_ptr), value :: model
  end function EmDee_shifted_force

  type(c_ptr) function EmDee_smoothed( model, Rm ) &
    bind(C,name="EmDee_smoothed")
    import :: c_ptr, c_double
    type(c_ptr),    value :: model
    real(c_double), value :: Rm
  end function EmDee_smoothed

  type(c_ptr) function EmDee_shifted_smoothed( model, Rm ) &
    bind(C,name="EmDee_shifted_smoothed")
    import :: c_ptr, c_double
    type(c_ptr),    value :: model
    real(c_double), value :: Rm
  end function EmDee_shifted_smoothed

  type(c_ptr) function EmDee_openmm_smoothed( model, Rm ) &
    bind(C,name="EmDee_openmm_smoothed")
    import :: c_ptr, c_double
    type(c_ptr),    value :: model
    real(c_double), value :: Rm
  end function EmDee_openmm_smoothed

  ! MODELS:
  type(c_ptr) function EmDee_pair_none( ) &
    bind(C,name="EmDee_pair_none")
    import :: c_ptr
  end function EmDee_pair_none

  type(c_ptr) function EmDee_coul_none( ) &
    bind(C,name="EmDee_coul_none")
    import :: c_ptr
  end function EmDee_coul_none

  type(c_ptr) function EmDee_bond_none( ) &
    bind(C,name="EmDee_bond_none")
    import :: c_ptr
  end function EmDee_bond_none

  type(c_ptr) function EmDee_angle_none( ) &
    bind(C,name="EmDee_angle_none")
    import :: c_ptr
  end function EmDee_angle_none

  type(c_ptr) function EmDee_dihedral_none( ) &
    bind(C,name="EmDee_dihedral_none")
    import :: c_ptr
  end function EmDee_dihedral_none

end interface
end module
