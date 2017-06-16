module EmDee

immutable tOptions
  Translate::Int32           # Flag to activate/deactivate translations
  Rotate::Int32              # Flag to activate/deactivate rotations
  RotationMode::Int32        # Algorithm used for free rotation of rigid bodies
  AutoForceCompute::Int32    # Flag to activate/deactivate automatic force computations
  AutoBodyUpdate::Int32      # Flag to activate/deactivate automatic rigid body update
end

immutable tVec3D
  x::Float64
  y::Float64
  z::Float64
end

immutable tEnergy
  Potential::Float64            # Total potential energy of the system
  Dispersion::Float64           # Dispersion (vdW) part of the potential energy
  Coulomb::Float64              # Electrostatic part of the potential energy
  Fourier::Float64              # Reciprocal part of the electrostatic potential
  Kinetic::Float64              # Total kinetic energy of the system
  TransPart::tVec3D             # Translational kinetic energy at each dimension
  Rotational::Float64           # Rotational kinetic energy of the system
  RotPart::tVec3D               # Rotational kinetic energy around each principal axis
  LayerPotential::Ref{Float64}  # Vector with multilayer potential energy components
  LayerDispersion::Ref{Float64} # Vector with multilayer dispersion energy components
  LayerCoulomb::Ref{Float64}    # Vector with multilayer coulombic energy components
  LayerFourier::Ref{Float64}    # Vector with multilayer reciprocal energy components
  Compute::Int32                # Flag to activate/deactivate energy computations
  UpToDate::Int32               # Flag to attest whether energies have been computed
end

immutable tTime
  Pair::Float64              # Time taken in force calculations
  FastPair::Float64
  Neighbor::Float64
  Total::Float64             # Total time since initialization
end

immutable tEmDee
  Builds::Int32              # Number of neighbor-list builds
  Time::tTime
  Energy::tEnergy            # All energy terms
  Virial::Float64            # Total internal virial of the system
  BodyVirial::Float64        # Rigid body contribution to the internal virial
  DoF::Int32                 # Total number of degrees of freedom
  RotDoF::Int32              # Number of rotational degrees of freedom
  Data::Ref{Void}            # Pointer to EmDee system data
  Options::tOptions          # List of options to change EmDee's behavior
end

typealias tModel Ref{Void}
typealias IntegerArray{T<:Integer} Array{T}
typealias RealArray{T<:Real} Array{T}

#---------------------------------------------------------------------------------------------------
function system( threads::Integer, layers::Integer, rc::Real, skin::Real, N::Integer,
                 types::IntegerArray, masses::RealArray, bodies::IntegerArray )
  return ccall( (:EmDee_system,"libemdee"), tEmDee,
                (Int32, Int32, Float64, Float64, Int32, Ptr{Int32}, Ptr{Float64}),
                threads, layers, rc, skin, N, Vector{Int32}(types), Vector{Float64}(masses),
                Vector{Int32}(bodies) )
end
#---------------------------------------------------------------------------------------------------
function EmDee_share_phase_space( mdkeep::tEmDee, mdlose::tEmDee )
  ccall( (:EmDee_share_phase_space,"libemdee"), Void,
         (tEmDee, Ptr{tEmDee}),
         mdkeep, Ref(mdlose) )
end
#---------------------------------------------------------------------------------------------------
function switch_model_layer( md::tEmDee, layer::Int )
  ccall( (:EmDee_switch_model_layer,"libemdee"), Void,
         (tEmDee, Int32),
         md, layer )
end
#---------------------------------------------------------------------------------------------------
function set_pair_model( md::tEmDee, itype::Integer, jtype::Integer, model::tModel, kCoul::Real )
  ccall( (:EmDee_set_pair_model,"libemdee"), Void,
         (tEmDee, Int32, Int32, tModel, Float64),
         md, itype, jtype, model, kCoul )
end
#---------------------------------------------------------------------------------------------------
function set_pair_multimodel( md::tEmDee, itype::Integer, jtype::Integer, model::Array{tModel},
                              kCoul::RealArray )
  ccall( (:EmDee_set_pair_multimodel,"libemdee"), Void,
         (tEmDee, Int32, Int32, Ptr{tModel}, Ptr{Float64}),
         md, itype, jtype, model, kCoul )
end
#---------------------------------------------------------------------------------------------------
function set_kspace_model( md::tEmDee, model::tModel )
  ccall( (:EmDee_set_kspace_model,"libemdee"), Void,
         (tEmDee, tModel),
         md, model )
end
#---------------------------------------------------------------------------------------------------
function set_coul_model( md::tEmDee, model::tModel )
  ccall( (:EmDee_set_coul_model,"libemdee"), Void,
         (tEmDee, tModel),
         md, model )
end
#---------------------------------------------------------------------------------------------------
function set_pair_multimodel( md::tEmDee, model::Array{tModel} )
  ccall( (:EmDee_set_coul_multimodel,"libemdee"), Void,
         (tEmDee, Ptr{tModel}),
         md, model )
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
function set_respa( md::tEmDee, InRc::Real, ExRc::Real, Npair::Integer, Nbond::Integer )
  ccall( (:EmDee_set_respa,"libemdee"), Void,
         (tEmDee, Float64, Float64, Int32, Int32),
         md, InRc, ExRc, Npair, Nbond )
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
function displace( md::tEmDee, lambda::Real, alpha::Real, dt::Real )
  ccall( (:EmDee_displace,"libemdee"), Void,
         (Ptr{tEmDee}, Float64, Float64, Float64),
         Ref(md), lambda, alpha, dt )
end
#---------------------------------------------------------------------------------------------------
function compute_forces( md::tEmDee )
  ccall( (:EmDee_compute_forces,"libemdee"), Void,
         (Ptr{tEmDee}),
         Ref(md) )
end
#---------------------------------------------------------------------------------------------------
function advance( md::tEmDee, alpha_R::Real, alpha_P::Real, dt::Real )
  ccall( (:EmDee_advance,"libemdee"), Void,
         (Ptr{tEmDee}, Float64, Float64, Float64),
         Ref(md), alpha_R, alpha_P, dt )
end
#---------------------------------------------------------------------------------------------------
function rdf( md::tEmDee, bins::Int, pairs::Int, itype::IntegerArray, jtype::IntegerArray, 
              g::RealArray )
  ccall( (:EmDee_rdf,"libemdee"), Void,
         (tEmDee, Int32, Int32, Vector{Int32}, Vector{Int32}, Matrix{Float64}),
          md, bins, pairs, itype, jtype, g )
end
#---------------------------------------------------------------------------------------------------
#                                            M O D E L S
#---------------------------------------------------------------------------------------------------
function shifted_force( model::tModel )
  return ccall( (:EmDee_shifted_force,"libemdee"), tModel, (tModel), model )
end
#---------------------------------------------------------------------------------------------------
function smoothed( model::tModel, Rm::Real )
  return ccall( (:EmDee_smoothed,"libemdee"), tModel, (tModel, Float64), model, Rm )
end
#---------------------------------------------------------------------------------------------------
function pair_none()
  return ccall( (:EmDee_pair_none,"libemdee"), tModel, () )
end
#---------------------------------------------------------------------------------------------------
function coul_none()
  return ccall( (:EmDee_coul_none,"libemdee"), tModel, () )
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
