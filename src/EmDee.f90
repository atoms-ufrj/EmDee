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

module EmDee

use iso_c_binding

implicit none

integer, parameter, private :: ib = c_int, rb = c_double

type, bind(C) :: EmDee_Model
  integer(ib) :: id = 0
  type(c_ptr) :: data = c_null_ptr
  real(rb)    :: p1 = 0.0_rb
  real(rb)    :: p2 = 0.0_rb
  real(rb)    :: p3 = 0.0_rb
  real(rb)    :: p4 = 0.0_rb
  integer(ib) :: external = 1
end type EmDee_Model

type, bind(C) :: tEmDee

  integer(ib) :: builds         ! Number of neighbor-list builds
  real(rb)    :: pairTime       ! Time taken in force calculations
  real(rb)    :: totalTime      ! Total time since initialization

  real(rb)    :: Potential      ! Total potential energy of the system
  real(rb)    :: Kinetic        ! Total kinetic energy of the system
  real(rb)    :: Rotational     ! Rotational kinetic energy of the system
  real(rb)    :: Virial         ! Total internal virial of the system

  real(rb)    :: Lbox           ! Length of the simulation box
  type(c_ptr) :: coords         ! Pointer to the coordinates of all atoms
  type(c_ptr) :: momenta        ! Pointer to the momenta of all atoms
  type(c_ptr) :: forces         ! Pointer to the resultant forces on all atoms
  type(c_ptr) :: charge         ! Pointer to the electric charges of all atoms

  real(rb)    :: Rc             ! Cut-off distance
  real(rb)    :: RcSq           ! Cut-off distance squared
  real(rb)    :: xRc            ! Extended cutoff distance (including skin)
  real(rb)    :: xRcSq          ! Extended cutoff distance squared
  real(rb)    :: skinSq         ! Square of the neighbor list skin width
  real(rb)    :: invL           ! Inverse length of the simulation box
  real(rb)    :: invL2          ! Squared inverse length of the simulation box
  real(rb)    :: totalMass      ! Sum of the masses of all atoms
  real(rb)    :: startTime      ! Time recorded at initialization
  real(rb)    :: eshift         ! Potential shifting factor for Coulombic interactions
  real(rb)    :: fshift         ! Force shifting factor for Coulombic interactions
  integer(ib) :: coulomb        ! Flag for coulombic interactions

  integer(ib) :: mcells         ! Number of cells at each dimension
  integer(ib) :: ncells         ! Total number of cells
  integer(ib) :: maxcells       ! Maximum number of cells
  integer(ib) :: maxatoms       ! Maximum number of atoms in a cell
  integer(ib) :: maxpairs       ! Maximum number of pairs formed by all atoms of a cell
  type(c_ptr) :: cell           ! Array containing all neighbor cells of each cell
  type(c_ptr) :: atomCell       ! Array containing the current cell of each atom

  integer(ib) :: natoms         ! Number of atoms in the system
  type(c_ptr) :: atomType       ! Pointer to the type indexes of all atoms
  type(c_ptr) :: atomMass       ! Pointer to the masses of all atoms
  type(c_ptr) :: invMass        ! Pointer to the inverses of atoms masses
  type(c_ptr) :: R0             ! Position of each atom at the latest neighbor list building

  integer(ib) :: ntypes         ! Number of atom types
  type(c_ptr) :: pairModel      ! Model of each type of atom pair

  type(c_ptr) :: bonds          ! List of bonds
  type(c_ptr) :: angles         ! List of angles
  type(c_ptr) :: dihedrals      ! List of dihedrals

  integer(ib) :: nbodies        ! Number of rigid bodies
  integer(ib) :: maxbodies      ! Maximum number of rigid bodies
  type(c_ptr) :: body           ! Pointer to the rigid bodies present in the system

  integer(ib) :: nfree          ! Number of independent atoms
  type(c_ptr) :: freeAtom       ! Pointer to the list of independent atoms

  integer(ib) :: nthreads       ! Number of parallel openmp threads
  integer(ib) :: threadAtoms    ! Number of atoms per parallel thread
  integer(ib) :: threadBodies   ! Number of rigid bodies per parallel thread
  type(c_ptr) :: cellAtom       ! List of atoms belonging to each cell
  type(c_ptr) :: threadCell     ! List of cells to be dealt with in each parallel thread
  type(c_ptr) :: neighbor       ! Pointer to neighbor lists
  type(c_ptr) :: excluded       ! List of pairs excluded from the neighbor lists
  type(c_ptr) :: random         ! Pointer for random number generators

end type tEmDee

interface

  function EmDee_system( threads, rc, skin, N, types, masses, seed ) result( me ) &
                                                                     bind(C,name="EmDee_system")
    import :: ib, rb, c_ptr, tEmDee
    integer(ib), value :: threads, N, seed
    real(rb),    value :: rc, skin
    type(c_ptr), value :: types, masses
    type(tEmDee)       :: me
  end function EmDee_system

  subroutine EmDee_set_charges( md, charges ) bind(C,name="EmDee_set_charges")
    import :: c_ptr
    type(c_ptr), value :: md, charges
  end subroutine EmDee_set_charges

  subroutine EmDee_set_pair_type( md, itype, jtype, model ) bind(C,name="EmDee_set_pair_type")
    import :: c_ptr, ib
    type(c_ptr), value :: md
    integer(ib), value :: itype, jtype
    type(c_ptr), value :: model
  end subroutine EmDee_set_pair_type

  subroutine EmDee_ignore_pair( md, i, j ) bind(C,name="EmDee_ignore_pair")
    import :: c_ptr, ib
    type(c_ptr), value :: md
    integer(ib), value :: i, j
  end subroutine EmDee_ignore_pair

  subroutine EmDee_add_bond( md, i, j, model ) bind(C,name="EmDee_add_bond")
    import :: c_ptr, ib
    type(c_ptr), value :: md
    integer(ib), value :: i, j
    type(c_ptr), value :: model
  end subroutine EmDee_add_bond

  subroutine EmDee_add_angle( md, i, j, k, model ) bind(C,name="EmDee_add_angle")
    import :: c_ptr, ib
    type(c_ptr), value :: md
    integer(ib), value :: i, j, k
    type(c_ptr), value :: model
  end subroutine EmDee_add_angle

  subroutine EmDee_add_dihedral( md, i, j, k, l, model ) bind(C,name="EmDee_add_dihedral")
    import :: c_ptr, ib
    type(c_ptr), value :: md
    integer(ib), value :: i, j, k, l
    type(c_ptr), value :: model
  end subroutine EmDee_add_dihedral

  subroutine EmDee_add_rigid_body( md, indexes, N ) bind(C,name="EmDee_add_rigid_body")
    import :: c_ptr, ib
    type(c_ptr), value :: md, indexes
    integer(ib), value :: N
  end subroutine EmDee_add_rigid_body

  subroutine EmDee_upload( md, L, coords, momenta, forces ) bind(C,name="EmDee_upload")
    import :: c_ptr
    type(c_ptr), value :: md, L, coords, momenta, forces
  end subroutine EmDee_upload

  subroutine EmDee_download( md, L, coords, momenta, forces ) bind(C,name="EmDee_download")
    import :: c_ptr
    type(c_ptr), value :: md, L, coords, momenta, forces
  end subroutine EmDee_download

  subroutine EmDee_random_momenta( md, kT, adjust ) bind(C,name="EmDee_random_momenta")
    import :: c_ptr, rb, ib
    type(c_ptr), value :: md
    real(rb),    value :: kT
    integer(ib), value :: adjust
  end subroutine EmDee_random_momenta

  subroutine EmDee_boost( md, lambda, alpha, dt ) bind(C,name="EmDee_boost")
    import :: c_ptr, rb
    type(c_ptr), value :: md
    real(rb),    value :: lambda, alpha, dt
  end subroutine EmDee_boost

  subroutine EmDee_move( md, lambda, alpha, dt ) bind(C,name="EmDee_move")
    import :: c_ptr, rb
    type(c_ptr), value :: md
    real(rb),    value :: lambda, alpha, dt
  end subroutine EmDee_move

  subroutine EmDee_compute( md ) bind(C,name="EmDee_compute")
    import :: c_ptr
    type(c_ptr), value :: md
  end subroutine EmDee_compute

  function EmDee_subset_energy( md, N, atoms, check ) result( Energy ) &
                                                      bind(C,name="EmDee_subset_energy")
    import :: c_ptr, ib, rb
    type(c_ptr), value :: md
    integer(ib), value :: N
    type(c_ptr), value :: atoms
    real(rb)           :: Energy
  end function EmDee_subset_energy

  function EmDee_pair_none( ) result(model) bind(C,name="EmDee_pair_none")
    import :: EmDee_Model
    type(EmDee_Model) :: model
  end function EmDee_pair_none

  function EmDee_pair_lj( epsilon, sigma ) result(model) bind(C,name="EmDee_pair_lj")
    import :: rb, EmDee_Model
    real(rb), value   :: epsilon, sigma
    type(EmDee_Model) :: model
  end function EmDee_pair_lj

  function EmDee_pair_lj_sf( epsilon, sigma, cutoff ) result(model) bind(C,name="EmDee_pair_lj_sf")
    import :: rb, EmDee_Model
    real(rb), value   :: epsilon, sigma, cutoff
    type(EmDee_Model) :: model
  end function EmDee_pair_lj_sf

  function EmDee_bond_harmonic( k, r0 ) result(model) bind(C,name="EmDee_bond_harmonic")
    import :: rb, EmDee_Model
    real(rb), value :: k, r0
    type(EmDee_Model) :: model
  end function EmDee_bond_harmonic

  function EmDee_bond_morse( D, alpha, r0 ) result(model) bind(C,name="EmDee_bond_morse")
    import :: rb, EmDee_Model
    real(rb), value :: D, alpha, r0
    type(EmDee_Model) :: model
  end function EmDee_bond_morse

  function EmDee_angle_harmonic( k, theta0 ) result(model) bind(C,name="EmDee_angle_harmonic")
    import :: rb, EmDee_Model
    real(rb), value :: k, theta0
    type(EmDee_Model) :: model
  end function EmDee_angle_harmonic

  function EmDee_dihedral_harmonic( k, phi0 ) result(model) bind(C,name="EmDee_dihedral_harmonic")
    import :: rb, EmDee_Model
    real(rb), value :: k, phi0
    type(EmDee_Model) :: model
  end function EmDee_dihedral_harmonic

end interface

end module EmDee


