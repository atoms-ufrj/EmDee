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

type, bind(C) :: tOptions
  integer(c_int) :: translate      ! Flag to activate/deactivate translations
  integer(c_int) :: rotate         ! Flag to activate/deactivate rotations
  integer(c_int) :: rotationMode   ! Algorithm used for free rotation of rigid bodies
end type tOptions

type, bind(C) :: tEmDee
  integer(c_int) :: builds         ! Number of neighbor-list builds
  real(c_double) :: pairTime       ! Time taken in force calculations
  real(c_double) :: totalTime      ! Total time since initialization
  real(c_double) :: Potential      ! Total potential energy of the system
  real(c_double) :: Kinetic        ! Total kinetic energy of the system
  real(c_double) :: Rotational     ! Rotational kinetic energy of the system
  real(c_double) :: Virial         ! Total internal virial of the system
  integer(c_int) :: DOF            ! Total number of degrees of freedom
  integer(c_int) :: RDOF           ! Number of rotational degrees of freedom
  type(c_ptr)    :: Data           ! Pointer to EmDee system data
  type(tOptions) :: Options        ! List of options to change EmDee's behavior
end type tEmDee

interface

  function EmDee_system( threads, layers, rc, skin, N, types, masses ) bind(C,name="EmDee_system")
    import
    integer(c_int), value :: threads, layers, N
    real(c_double), value :: rc, skin
    type(c_ptr),    value :: types, masses
    type(tEmDee)          :: EmDee_system
  end function EmDee_system

  subroutine EmDee_switch_model_layer( md, layer ) bind(C,name="EmDee_set_layer")
    import
    type(tEmDee),   intent(inout) :: md
    integer(c_int), value         :: layer
  end subroutine EmDee_switch_model_layer

  subroutine EmDee_set_charges( md, charges ) bind(C,name="EmDee_set_charges")
    import
    type(tEmDee), value :: md
    type(c_ptr),  value :: charges
  end subroutine EmDee_set_charges

  subroutine EmDee_set_pair_type( md, itype, jtype, model ) bind(C,name="EmDee_set_pair_type")
    import
    type(tEmDee),   value :: md
    integer(c_int), value :: itype, jtype
    type(c_ptr),    value :: model
  end subroutine EmDee_set_pair_type

  subroutine EmDee_ignore_pair( md, i, j ) bind(C,name="EmDee_ignore_pair")
    import
    type(tEmDee),   value :: md
    integer(c_int), value :: i, j
  end subroutine EmDee_ignore_pair

  subroutine EmDee_add_bond( md, i, j, model ) bind(C,name="EmDee_add_bond")
    import
    type(tEmDee),   value :: md
    integer(c_int), value :: i, j
    type(c_ptr),    value :: model
  end subroutine EmDee_add_bond

  subroutine EmDee_add_angle( md, i, j, k, model ) bind(C,name="EmDee_add_angle")
    import
    type(tEmDee),   value :: md
    integer(c_int), value :: i, j, k
    type(c_ptr),    value :: model
  end subroutine EmDee_add_angle

  subroutine EmDee_add_dihedral( md, i, j, k, l, model ) bind(C,name="EmDee_add_dihedral")
    import
    type(tEmDee),   value :: md
    integer(c_int), value :: i, j, k, l
    type(c_ptr),    value :: model
  end subroutine EmDee_add_dihedral

  subroutine EmDee_add_rigid_body( md, N, indexes ) bind(C,name="EmDee_add_rigid_body")
    import
    type(tEmDee),   value :: md
    type(c_ptr),    value :: indexes
    integer(c_int), value :: N
  end subroutine EmDee_add_rigid_body

  subroutine EmDee_upload( md, Lbox, coords, momenta, forces ) bind(C,name="EmDee_upload")
    import
    type(tEmDee), intent(inout) :: md
    type(c_ptr),  value         :: Lbox, coords, momenta, forces
  end subroutine EmDee_upload

  subroutine EmDee_download( md, Lbox, coords, momenta, forces ) bind(C,name="EmDee_download")
    import
    type(tEmDee), value :: md
    type(c_ptr),  value :: Lbox, coords, momenta, forces
  end subroutine EmDee_download

  subroutine EmDee_random_momenta( md, kT, adjust, seed ) bind(C,name="EmDee_random_momenta")
    import
    type(tEmDee),   intent(inout) :: md
    real(c_double), value         :: kT
    integer(c_int), value         :: adjust, seed
  end subroutine EmDee_random_momenta

!  subroutine EmDee_save_state( md, rigid )
!    import
!    type(tEmDee),   intent(inout) :: md
!    integer(c_int), intent(in)    :: rigid
!  end subroutine EmDee_save_state

!  subroutine EmDee_restore_state( md )
!    import
!    type(tEmDee), intent(inout) :: md
!  end subroutine EmDee_restore_state

  subroutine EmDee_boost( md, lambda, alpha, dt ) bind(C,name="EmDee_boost")
    import
    type(tEmDee),   intent(inout) :: md
    real(c_double), value         :: lambda, alpha, dt
  end subroutine EmDee_boost

  subroutine EmDee_move( md, lambda, alpha, dt ) bind(C,name="EmDee_move")
    import
    type(tEmDee),   intent(inout) :: md
    real(c_double), value         :: lambda, alpha, dt
  end subroutine EmDee_move

  subroutine EmDee_group_energy( md, na, atoms, ne, energies ) bind(C,name="EmDee_group_energy")
    import
    type(tEmDee),   value :: md
    integer(c_int), value :: na, ne
    type(c_ptr),    value :: atoms, energies
  end subroutine EmDee_group_energy

  function EmDee_model_none( ) bind(C,name="EmDee_model_none")
    import
    type(c_ptr) :: EmDee_model_none
  end function EmDee_model_none

  function EmDee_pair_none( ) bind(C,name="EmDee_pair_none")
    import
    type(c_ptr) :: EmDee_pair_none
  end function EmDee_pair_none

  function EmDee_pair_lj_cut( epsilon, sigma ) bind(C,name="EmDee_pair_lj_cut")
    import
    real(c_double), value :: epsilon, sigma
    type(c_ptr)           :: EmDee_pair_lj_cut
  end function EmDee_pair_lj_cut

  type(c_ptr) function EmDee_pair_lj_cut_coul_cut( epsilon, sigma ) bind(C,name="EmDee_pair_lj_cut_coul_cut")
    import
    real(c_double), value :: epsilon, sigma
  end function EmDee_pair_lj_cut_coul_cut

  function EmDee_pair_lj_sf( epsilon, sigma ) bind(C,name="EmDee_pair_lj_sf")
    import
    real(c_double), value :: epsilon, sigma
    type(c_ptr)           :: EmDee_pair_lj_sf
  end function EmDee_pair_lj_sf

  function EmDee_pair_lj_sf_coul_sf( epsilon, sigma ) bind(C,name="EmDee_pair_lj_sf_coul_sf")
    import
    real(c_double), value :: epsilon, sigma
    type(c_ptr)           :: EmDee_pair_lj_sf_coul_sf
  end function EmDee_pair_lj_sf_coul_sf

!  function EmDee_pair_softcore( epsilon, sigma, lambda ) bind(C,name="EmDee_pair_softcore")
!    import
!    real(c_double), value :: epsilon, sigma, lambda
!    type(c_ptr)           :: EmDee_pair_softcore
!  end function EmDee_pair_softcore

  function EmDee_bond_harmonic( k, r0 ) bind(C,name="EmDee_bond_harmonic")
    import
    real(c_double), value :: k, r0
    type(c_ptr)           :: EmDee_bond_harmonic
  end function EmDee_bond_harmonic

!  function EmDee_bond_morse( D, alpha, r0 ) bind(C,name="EmDee_bond_morse")
!    import
!    real(c_double), value :: D, alpha, r0
!    type(c_ptr)           :: EmDee_bond_morse
!  end function EmDee_bond_morse

!  function EmDee_angle_harmonic( k, theta0 ) bind(C,name="EmDee_angle_harmonic")
!    import
!    real(c_double), value :: k, theta0
!    type(c_ptr)           :: EmDee_angle_harmonic
!  end function EmDee_angle_harmonic

!  function EmDee_dihedral_harmonic( k, phi0 ) bind(C,name="EmDee_dihedral_harmonic")
!    import
!    real(c_double), value :: k, phi0
!    type(c_ptr)           :: EmDee_dihedral_harmonic
!  end function EmDee_dihedral_harmonic

  subroutine EmDee_Rotational_Energies( md, Kr ) bind(C,name="EmDee_Rotational_Energies")
    import
    type(tEmDee),   value       :: md
    real(c_double), intent(out) :: Kr(3) 
  end subroutine EmDee_Rotational_Energies

end interface

