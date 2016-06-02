module EmDee

type Model
  data::Ptr{Void}
  params::Ptr{Void}
end

type tEmDee

  builds::Int32           # Number of neighbor-list builds
  time::Float64           # Total time taken in force calculations
  Rc::Float64             # Cut-off distance
  RcSq::Float64           # Cut-off distance squared
  xRc::Float64            # Extended cutoff distance (including skin)
  xRcSq::Float64          # Extended cutoff distance squared
  skinSq::Float64         # Square of the neighbor list skin width

  mcells::Int32           # Number of cells at each dimension
  ncells::Int32           # Total number of cells
  maxcells::Int32         # Maximum number of cells
  maxatoms::Int32         # Maximum number of atoms in a cell
  maxpairs::Int32         # Maximum number of pairs formed by all atoms of a cell
  cell::Ptr{Void}         # Array containing all neighbor cells of each cell

  natoms::Int32           # Number of atoms in the system
  atomType::Ptr{Int32}    # The type of each atom
  R0::Ptr{Float64}        # The position of each atom at the latest neighbor list building
  charge::Ptr{Float64}    # Pointer to the electric charge of each atom

  ntypes::Int32           # Number of atom types
  pairParams::Ptr{Void}   # Model parameters of each type of atom pair
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

function compute_forces( md::tEmDee, forces::Array{Float64,2}, coords::Array{Float64,2}, L::Real )
  ccall( (:EmDee_compute_forces, "libemdee"), Void, (Ptr{Void}, Ptr{Float64}, Ptr{Float64}, Float64),
         pointer_from_objref(md), forces, coords, L )
end

#---------------------------------------------------------------------------------------------------

function pair_lj( sigma::Real, epsilon::Real )
  return ccall( (:EmDee_pair_lj, "libemdee"), Model, (Float64, Float64), sigma, epsilon )
end

#---------------------------------------------------------------------------------------------------

end
