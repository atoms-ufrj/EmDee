module EmDee

type Model
  id::Int32
  data::Ptr{Void}
  p1::Float64
  p2::Float64
  p3::Float64
  p4::Float64
  external::Int32
end

Model() = Model(Int32(0),C_NULL,Float64(0),Float64(0),Float64(0),Float64(0),Int32(0))

type tEmDee

  builds::Int32           # Number of neighbor-list builds
  time::Float64           # Total time taken in force calculations
  Rc::Float64             # Cut-off distance
  RcSq::Float64           # Cut-off distance squared
  xRc::Float64            # Extended cutoff distance (including skin)
  xRcSq::Float64          # Extended cutoff distance squared
  skinSq::Float64         # Square of the neighbor list skin width
  eshift::Float64         # Energy shifting factor for Coulombic interactions
  fshift::Float64         # Force shifting factor for Coulombic interactions

  mcells::Int32           # Number of cells at each dimension
  ncells::Int32           # Total number of cells
  maxcells::Int32         # Maximum number of cells
  maxatoms::Int32         # Maximum number of atoms in a cell
  maxpairs::Int32         # Maximum number of pairs formed by all atoms of a cell
  cell::Ptr{Void}         # Array containing all neighbor cells of each cell

  natoms::Int32           # Number of atoms in the system
  atomType::Ptr{Int32}    # The type of each atom
  R0::Ptr{Float64}        # The position of each atom at the latest neighbor list building

  coulomb::Int32          # Flag for coulombic interactions
  charge::Ptr{Float64}    # Pointer to the electric charge of each atom

  ntypes::Int32           # Number of atom types
  pairData::Ptr{Void}     # Model data of each type of atom pair

  bond::Ptr{Void}         # List of bonds
  angle::Ptr{Void}        # List of angles
  dihedral::Ptr{Void}     # List of dihedrals

  Energy::Float64         # Total potential energy of the system
  Virial::Float64         # Total internal virial of the system

  nthreads::Int32         # Number of parallel openmp threads
  cellAtom::Ptr{Void}     # List of atoms belonging to each cell
  threadCell::Ptr{Void}   # List of cells to be dealt with in each parallel thread
  neighbor::Ptr{Void}     # Pointer to neighbor lists
  excluded::Ptr{Void}     # List of pairs excluded from the neighbor lists

end

#---------------------------------------------------------------------------------------------------

function system( threads::Int, rc::Real, skin::Real, atoms::Int, types::Array{Int,1} )
  return ccall( (:EmDee_system, "libemdee"), tEmDee,
                (Int32, Float64, Float64, Int32, Ptr{Int32}),
                threads, rc, skin, atoms, Array{Int32,1}(types) )
end

#---------------------------------------------------------------------------------------------------

function set_pair( md::tEmDee, itype::Int, jtype::Int, model::Model )
  ccall( (:EmDee_set_pair, "libemdee"), Void, (Ptr{Void}, Int32, Int32, Ptr{Void}),
         pointer_from_objref(md), itype, jtype, pointer_from_objref(model) )
end

#---------------------------------------------------------------------------------------------------

function set_charges( md::tEmDee, charges::Array{Float64,1} )
  ccall( (:EmDee_set_charges, "libemdee"), Void, (Ptr{Void}, Ptr{Float64}),
         pointer_from_objref(md), pointer_from_objref(charges) )
end

#---------------------------------------------------------------------------------------------------

function add_bond( md::tEmDee, i::Int, j::Int, model::Model )
  ccall( (:EmDee_add_bond, "libemdee"), Void, (Ptr{Void}, Int32, Int32, Ptr{Void}),
         pointer_from_objref(md), i, j, pointer_from_objref(model) )
end

#---------------------------------------------------------------------------------------------------

function add_angle( md::tEmDee, i::Int, j::Int, k::Int, model::Model )
  ccall( (:EmDee_add_angle, "libemdee"), Void, (Ptr{Void}, Int32, Int32, Int32, Ptr{Void}),
         pointer_from_objref(md), i, j, k, pointer_from_objref(model) )
end

#---------------------------------------------------------------------------------------------------

function add_dihedral( md::tEmDee, i::Int, j::Int, k::Int, l::Int32, model::Model )
  ccall( (:EmDee_add_dihedral, "libemdee"), Void, (Ptr{Void}, Int32, Int32, Int32, Int32, Ptr{Void}),
         pointer_from_objref(md), i, j, k, l, pointer_from_objref(model) )
end

#---------------------------------------------------------------------------------------------------

function ignore_pair( md::tEmDee, i::Int, j::Int )
  ccall( (:EmDee_ignore_pair, "libemdee"), Void, (Ptr{Void}, Int32, Int32),
         pointer_from_objref(md), i, j )
end

#---------------------------------------------------------------------------------------------------

function compute( md::tEmDee, forces::Array{Float64,2}, coords::Array{Float64,2}, L::Real )
  ccall( (:EmDee_compute, "libemdee"), Void, (Ptr{Void}, Ptr{Float64}, Ptr{Float64}, Float64),
         pointer_from_objref(md), forces, coords, L )
end

#---------------------------------------------------------------------------------------------------
#                                            M O D E L S
#---------------------------------------------------------------------------------------------------

function pair_lj( ɛ::Number, σ::Number )
  return ccall( (:EmDee_pair_lj, "libemdee"), Model, (Float64, Float64), ɛ, σ )
end

#---------------------------------------------------------------------------------------------------

function pair_lj_sf( ɛ::Number, σ::Number, rc::Number )
  return ccall( (:EmDee_pair_lj_sf, "libemdee"), Model, (Float64, Float64, Float64), ɛ, σ, rc )
end

#---------------------------------------------------------------------------------------------------

function bond_harmonic( k::Number, r₀::Number )
  return ccall( (:EmDee_bond_harmonic, "libemdee"), Model, (Float64, Float64), k, r₀ )
end

#---------------------------------------------------------------------------------------------------

function bond_morse( D::Number, α::Number, r₀::Number )
  return ccall( (:EmDee_bond_morse, "libemdee"), Model, (Float64, Float64, Float64), D, α, r₀ )
end

#---------------------------------------------------------------------------------------------------

function angle_harmonic( k::Number, θ₀::Number )
  return ccall( (:EmDee_angle_harmonic, "libemdee"), Model, (Float64, Float64), k, θ₀ )
end

#---------------------------------------------------------------------------------------------------

function dihedral_harmonic( k::Number, ϕ₀::Number )
  return ccall( (:EmDee_dihedral_harmonic, "libemdee"), Model, (Float64, Float64), k, ϕ₀ )
end

#---------------------------------------------------------------------------------------------------

end
