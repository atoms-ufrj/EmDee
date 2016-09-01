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
  integer(ib) :: id
  type(c_ptr) :: data
  real(rb)    :: p1, p2, p3, p4
  integer(ib) :: external
  type(c_ptr) :: next
end type EmDee_Model

type, bind(C) :: tEmDee

  integer(ib) :: builds         ! Number of neighbor-list builds
  real(rb)    :: pairTime       ! Time taken in force calculations
  real(rb)    :: totalTime      ! Total time since initialization

  real(rb)    :: Potential      ! Total potential energy of the system
  real(rb)    :: Kinetic        ! Total kinetic energy of the system
  real(rb)    :: Rotational     ! Rotational kinetic energy of the system
  real(rb)    :: Virial         ! Total internal virial of the system

  type(c_ptr) :: Data           ! Pointer to EmDee system data

end type tEmDee

interface

  function EmDee_system( threads, rc, skin, N, types, masses, seed ) result( md ) &
                                                                     bind(C,name="EmDee_system")
    import :: ib, rb, c_ptr, tEmDee
    integer(ib), value :: threads, N, seed
    real(rb),    value :: rc, skin
    type(c_ptr), value :: types, masses
    type(tEmDee)       :: md
  end function EmDee_system

  subroutine EmDee_set_charges( md, charges ) bind(C,name="EmDee_set_charges")
    import :: tEmDee, c_ptr
    type(tEmDee), value :: md
    type(c_ptr),  value :: charges
  end subroutine EmDee_set_charges

  subroutine EmDee_set_pair_type( md, itype, jtype, model ) bind(C,name="EmDee_set_pair_type")
    import :: tEmDee, c_ptr, ib
    type(tEmDee), value :: md
    integer(ib),  value :: itype, jtype
    type(c_ptr),  value :: model
  end subroutine EmDee_set_pair_type

  subroutine EmDee_ignore_pair( md, i, j ) bind(C,name="EmDee_ignore_pair")
    import :: tEmDee, c_ptr, ib
    type(tEmDee), value :: md
    integer(ib),  value :: i, j
  end subroutine EmDee_ignore_pair

  subroutine EmDee_add_bond( md, i, j, model ) bind(C,name="EmDee_add_bond")
    import :: tEmDee, c_ptr, ib
    type(tEmDee), value :: md
    integer(ib),  value :: i, j
    type(c_ptr),  value :: model
  end subroutine EmDee_add_bond

  subroutine EmDee_add_angle( md, i, j, k, model ) bind(C,name="EmDee_add_angle")
    import :: tEmDee, c_ptr, ib
    type(tEmDee), value :: md
    integer(ib),  value :: i, j, k
    type(c_ptr),  value :: model
  end subroutine EmDee_add_angle

  subroutine EmDee_add_dihedral( md, i, j, k, l, model ) bind(C,name="EmDee_add_dihedral")
    import :: tEmDee, c_ptr, ib
    type(tEmDee), value :: md
    integer(ib),  value :: i, j, k, l
    type(c_ptr),  value :: model
  end subroutine EmDee_add_dihedral

  subroutine EmDee_add_rigid_body( md, N, indexes ) bind(C,name="EmDee_add_rigid_body")
    import :: tEmDee, c_ptr, ib
    type(tEmDee), value :: md
    type(c_ptr),  value :: indexes
    integer(ib),  value :: N
  end subroutine EmDee_add_rigid_body

  subroutine EmDee_upload( md, Lbox, coords, momenta, forces ) bind(C,name="EmDee_upload")
    import :: tEmDee, c_ptr
    type(tEmDee), intent(inout) :: md
    type(c_ptr),  value         :: Lbox, coords, momenta, forces
  end subroutine EmDee_upload

  subroutine EmDee_download( md, Lbox, coords, momenta, forces ) bind(C,name="EmDee_download")
    import :: tEmDee, c_ptr
    type(tEmDee), value :: md
    type(c_ptr),  value :: Lbox, coords, momenta, forces
  end subroutine EmDee_download

  subroutine EmDee_random_momenta( md, kT, adjust ) bind(C,name="EmDee_random_momenta")
    import :: tEmDee, c_ptr, rb, ib
    type(tEmDee), intent(inout) :: md
    real(rb),     value         :: kT
    integer(ib),  value         :: adjust
  end subroutine EmDee_random_momenta

  subroutine EmDee_boost( md, lambda, alpha, dt, t_flag, r_flag ) bind(C,name="EmDee_boost")
    import :: tEmDee, c_ptr, rb, ib
    type(tEmDee), intent(inout) :: md
    real(rb),     value         :: lambda, alpha, dt
    integer(ib),  value         :: t_flag, r_flag
  end subroutine EmDee_boost

  subroutine EmDee_move( md, lambda, alpha, dt ) bind(C,name="EmDee_move")
    import :: tEmDee, c_ptr, rb
    type(tEmDee), intent(inout) :: md
    real(rb),     value         :: lambda, alpha, dt
  end subroutine EmDee_move

  subroutine EmDee_compute( md ) bind(C,name="EmDee_compute")
    import :: tEmDee
    type(tEmDee), intent(inout) :: md
  end subroutine EmDee_compute

  subroutine EmDee_group_energy( md, na, atoms, ne, energies ) bind(C,name="EmDee_group_energy")
    import :: tEmDee, ib, c_ptr
    type(tEmDee), value :: md
    integer(ib),  value :: na, ne
    type(c_ptr),  value :: atoms, energies
  end subroutine EmDee_group_energy

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


