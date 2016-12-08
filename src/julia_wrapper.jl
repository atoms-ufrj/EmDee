module EmDee

type tOptions
  translate::Int32        # Flag to activate/deactivate translations
  rotate::Int32           # Flag to activate/deactivate rotations
  rotationMode::Int32     # Algorithm used for free rotation of rigid bodies
end

type tEmDee
  builds::Int32           # Number of neighbor-list builds
  pairTime::Float64       # Time taken in force calculations
  totalTime::Float64      # Total time since initialization
  Potential::Float64      # Total potential energy of the system
  Kinetic::Float64        # Total kinetic energy of the system
  Rotational::Float64     # Rotational kinetic energy of the system
  Virial::Float64         # Total internal virial of the system
  DOF::Int32              # Total number of degrees of freedom
  RDOF::Int32             # Number of rotational degrees of freedom
  Data::Ptr{Void}         # Pointer to EmDee system data
  Options::tOptions       # List of options to change EmDee's behavior
end

typealias VecOrPtr Union{Vector{Float64},Ptr{Void}}
typealias MatOrPtr Union{Matrix{Float64},Ptr{Void}}
typealias IVecOrPtr Union{Vector{Int32},Ptr{Void}}

#---------------------------------------------------------------------------------------------------
function system( threads::Int, layers::Int, rc::Float64, skin::Float64, N::Int,
                 types::Vector{Int}, masses::VecOrPtr )
  return ccall( (:EmDee_system,"libemdee"), tEmDee,
                (Int32, Int32, Float64, Float64, Int32, Ptr{Int32}, Ptr{Float64}),
                threads, layers, rc, skin, N, Vector{Int32}(types), masses )
end

function system( threads::Int, layers::Int, rc::Float64, skin::Float64, N::Int,
                 types::Ptr{Void}, masses::VecOrPtr )
  return ccall( (:EmDee_system,"libemdee"), tEmDee,
                (Int32, Int32, Float64, Float64, Int32, Ptr{Int32}, Ptr{Float64}),
                threads, layers, rc, skin, N, types, masses )
end
#---------------------------------------------------------------------------------------------------
function set_pair_model( md::tEmDee, itype::Int, jtype::Int, model::Ptr{Void} )
  ccall( (:EmDee_set_pair_model,"libemdee"), Void,
         (tEmDee, Int32, Int32, Ptr{Void}),
         md, itype, jtype, model )
end
#---------------------------------------------------------------------------------------------------
function set_pair_multimodel( md::tEmDee, itype::Int, jtype::Int, model::Vector{Ptr{Void}} )
  ccall( (:EmDee_set_pair_multimodel,"libemdee"), Void,
         (tEmDee, Int32, Int32, Ptr{Void}),
         md, itype, jtype, model )
end
#---------------------------------------------------------------------------------------------------
function switch_model_layer( md::tEmDee, layer::Int )
  ccall( (:EmDee_switch_model_layer,"libemdee"), Void, (tEmDee, Int32), md, layer )
end
#---------------------------------------------------------------------------------------------------
function set_charges( md::tEmDee, charges::Vector{Float64} )
  ccall( (:EmDee_set_charges,"libemdee"), Void, (tEmDee, Ptr{Float64}), md, charges )
end
#---------------------------------------------------------------------------------------------------
function add_bond( md::tEmDee, i::Int, j::Int, model::Ptr{Void} )
  ccall( (:EmDee_add_bond,"libemdee"), Void,
         (tEmDee, Int32, Int32, Ptr{Void}), md, i, j, model )
end
#---------------------------------------------------------------------------------------------------
function add_angle( md::tEmDee, i::Int, j::Int, k::Int, model::Ptr{Void} )
  ccall( (:EmDee_add_angle,"libemdee"), Void,
         (tEmDee, Int32, Int32, Int32, Ptr{Void}), md, i, j, k, model )
end
#---------------------------------------------------------------------------------------------------
function add_dihedral( md::tEmDee, i::Int, j::Int, k::Int, l::Int, model::Ptr{Void} )
  ccall( (:EmDee_add_dihedral,"libemdee"), Void,
         (tEmDee, Int32, Int32, Int32, Int32, Ptr{Void}), md, i, j, k, l, model )
end
#---------------------------------------------------------------------------------------------------
function add_rigid_body( md::tEmDee, N::Int, indexes::Vector{Int} )
  ccall( (:EmDee_add_rigid_body,"libemdee"), Void,
         (tEmDee, Int32, Ptr{Int32}), md, N, indexes )
end
#---------------------------------------------------------------------------------------------------
function ignore_pair( md::tEmDee, i::Int, j::Int )
  ccall( (:EmDee_ignore_pair,"libemdee"), Void, (tEmDee, Int32, Int32), md, i, j )
end
#---------------------------------------------------------------------------------------------------
function upload( md::tEmDee, Lbox::VecOrPtr, coords::MatOrPtr,
                 momenta::MatOrPtr, forces::MatOrPtr )
  ccall( (:EmDee_upload,"libemdee"), Void,
         (Ptr{Void}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}),
         pointer_from_objref(md), Lbox, coords, momenta, forces )
end
#---------------------------------------------------------------------------------------------------
function download( md::tEmDee, Lbox::VecOrPtr, coords::MatOrPtr,
                   momenta::MatOrPtr, forces::MatOrPtr )
  ccall( (:EmDee_download,"libemdee"), Void,
         (tEmDee, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}),
         md, Lbox, coords, momenta, forces )
end
#---------------------------------------------------------------------------------------------------
function random_momenta( md::tEmDee, kT::Float64, adjust::Int, seed::Int )
  ccall( (:EmDee_random_momenta,"libemdee"), Void,
         (Ptr{Void}, Float64, Int32, Int32),
         pointer_from_objref(md), kT, adjust, seed )
end
#---------------------------------------------------------------------------------------------------
function boost( md::tEmDee, lambda::Float64, alpha::Float64, dt::Float64 )
  ccall( (:EmDee_boost,"libemdee"), Void,
         (Ptr{Void}, Float64, Float64, Float64),
         pointer_from_objref(md), lambda, alpha, dt )
end
#---------------------------------------------------------------------------------------------------
function move( md::tEmDee, lambda::Float64, alpha::Float64, dt::Float64 )
  ccall( (:EmDee_move,"libemdee"), Void,
         (Ptr{Void}, Float64, Float64, Float64),
         pointer_from_objref(md), lambda, alpha, dt )
end
#---------------------------------------------------------------------------------------------------
#                                            M O D E L S
#---------------------------------------------------------------------------------------------------
function pair_none()
  return ccall( (:EmDee_pair_none,"libemdee"), Ptr{Void}, () )
end
#---------------------------------------------------------------------------------------------------
function bond_none()
  return ccall( (:EmDee_bond_none,"libemdee"), Ptr{Void}, () )
end
#---------------------------------------------------------------------------------------------------
function angle_none()
  return ccall( (:EmDee_angle_none,"libemdee"), Ptr{Void}, () )
end
#---------------------------------------------------------------------------------------------------
function dihedral_none()
  return ccall( (:EmDee_dihedral_none,"libemdee"), Ptr{Void}, () )
end
#---------------------------------------------------------------------------------------------------
end
