module EmDee

immutable tOptions
  translate::Int32           # Flag to activate/deactivate translations
  rotate::Int32              # Flag to activate/deactivate rotations
  computeProps::Int32        # Flag to activate/deactivate energy and virial computations
  rotationMode::Int32        # Algorithm used for free rotation of rigid bodies
end

immutable tEmDee
  builds::Int32              # Number of neighbor-list builds
  pairTime::Float64          # Time taken in force calculations
  totalTime::Float64         # Total time since initialization
  Potential::Float64         # Total potential energy of the system
  Kinetic::Float64           # Total kinetic energy of the system
  Rotational::Float64        # Rotational kinetic energy of the system
  Virial::Float64            # Total internal virial of the system
  DOF::Int32                 # Total number of degrees of freedom
  RDOF::Int32                # Number of rotational degrees of freedom
  UpToDate::Int32            # Flag to attest whether energy and virial have been computated
  Data::Ref{Void}            # Pointer to EmDee system data
  Options::tOptions          # List of options to change EmDee's behavior
end

typealias tModel Ref{Void}
typealias IntegerArray{T<:Integer} Array{T}
typealias RealArray{T<:Real} Array{T}

#---------------------------------------------------------------------------------------------------
function system( threads::Integer, layers::Integer, rc::Real, skin::Real, N::Integer,
                 types::IntegerArray, masses::RealArray )
  return ccall( (:EmDee_system,"libemdee"), tEmDee,
                (Int32, Int32, Float64, Float64, Int32, Ptr{Int32}, Ptr{Float64}),
                threads, layers, rc, skin, N, Vector{Int32}(types), Vector{Float64}(masses) )
end
#---------------------------------------------------------------------------------------------------
function set_pair_model( md::tEmDee, itype::Integer, jtype::Integer, model::tModel )
  ccall( (:EmDee_set_pair_model,"libemdee"), Void,
         (tEmDee, Int32, Int32, tModel),
         md, itype, jtype, model )
end
#---------------------------------------------------------------------------------------------------
function set_pair_multimodel( md::tEmDee, itype::Integer, jtype::Integer, model::Array{tModel} )
  ccall( (:EmDee_set_pair_multimodel,"libemdee"), Void,
         (tEmDee, Int32, Int32, tModel),
         md, itype, jtype, model )
end
#---------------------------------------------------------------------------------------------------
function switch_model_layer( md::tEmDee, layer::Int )
  ccall( (:EmDee_switch_model_layer,"libemdee"), Void,
         (tEmDee, Int32),
         md, layer )
end
#---------------------------------------------------------------------------------------------------
function add_bond( md::tEmDee, i::Integer, j::Integer, model::tModel )
  ccall( (:EmDee_add_bond,"libemdee"), Void,
         (tEmDee, Int32, Int32, tModel),
         md, i, j, model )
end
#---------------------------------------------------------------------------------------------------
function add_angle( md::tEmDee, i::Integer, j::Integer, k::Integer, model::tModel )
  ccall( (:EmDee_add_angle,"libemdee"), Void,
         (tEmDee, Int32, Int32, Int32, tModel),
         md, i, j, k, model )
end
#---------------------------------------------------------------------------------------------------
function add_dihedral( md::tEmDee, i::Integer, j::Integer, k::Integer, l::Integer, model::tModel )
  ccall( (:EmDee_add_dihedral,"libemdee"), Void,
         (tEmDee, Int32, Int32, Int32, Int32, tModel),
         md, i, j, k, l, model )
end
#---------------------------------------------------------------------------------------------------
function add_rigid_body( md::tEmDee, N::Integer, indexes::IntegerArray )
  ccall( (:EmDee_add_rigid_body,"libemdee"), Void,
         (tEmDee, Int32, Ptr{Int32}),
         md, N, Vector{Int32}(indexes) )
end
#---------------------------------------------------------------------------------------------------
function ignore_pair( md::tEmDee, i::Integer, j::Integer )
  ccall( (:EmDee_ignore_pair,"libemdee"), Void,
         (tEmDee, Int32, Int32),
         md, i, j )
end
#---------------------------------------------------------------------------------------------------
function upload( md::tEmDee, option::String, address::RealArray )
  ccall( (:EmDee_upload,"libemdee"), Void,
         (Ptr{tEmDee}, Cstring, Ptr{Float64}),
         Ref(md), option, Array{Float64}(address) )
end
#---------------------------------------------------------------------------------------------------
function download( md::tEmDee, option::String, address::RealArray )
  ccall( (:EmDee_download,"libemdee"), Void,
         (tEmDee, Cstring, Ptr{Float64}),
         md, option, Array{Float64}(address) )
end
#---------------------------------------------------------------------------------------------------
function random_momenta( md::tEmDee, kT::Real, adjust::Integer, seed::Integer )
  ccall( (:EmDee_random_momenta,"libemdee"), Void,
         (Ptr{tEmDee}, Float64, Int32, Int32),
         Ref(md), kT, adjust, seed )
end
#---------------------------------------------------------------------------------------------------
function boost( md::tEmDee, lambda::Real, alpha::Real, dt::Real )
  ccall( (:EmDee_boost,"libemdee"), Void,
         (Ptr{tEmDee}, Float64, Float64, Float64),
         Ref(md), lambda, alpha, dt )
end
#---------------------------------------------------------------------------------------------------
function move( md::tEmDee, lambda::Real, alpha::Real, dt::Real )
  ccall( (:EmDee_move,"libemdee"), Void,
         (Ptr{tEmDee}, Float64, Float64, Float64),
         Ref(md), lambda, alpha, dt )
end
#---------------------------------------------------------------------------------------------------
#                                            M O D E L S
#---------------------------------------------------------------------------------------------------
function pair_none()
  return ccall( (:EmDee_pair_none,"libemdee"), tModel, () )
end
#---------------------------------------------------------------------------------------------------
function bond_none()
  return ccall( (:EmDee_bond_none,"libemdee"), tModel, () )
end
#---------------------------------------------------------------------------------------------------
function angle_none()
  return ccall( (:EmDee_angle_none,"libemdee"), tModel, () )
end
#---------------------------------------------------------------------------------------------------
function dihedral_none()
  return ccall( (:EmDee_dihedral_none,"libemdee"), tModel, () )
end
#---------------------------------------------------------------------------------------------------
end
