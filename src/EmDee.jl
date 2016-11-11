module EmDee

#TODO: REMOVE THIS DEBUG MESSAGE
  print("\n","running test julia from", pwd(), "\n","\n")

#TODO: 1) should this be here or elsewhere? note that the print runs twice in the terminal when running the test.
#TODO: 2) how to choose the path i.e. relative or absolute, hardcoded or read from file or passed from some makefile DEFINE.
  EmDee_libvar=Libdl.dlopen("./lib/libemdee.so")

  macro dlsym(func, lib) #REF: http://docs.julialang.org/en/release-0.5/manual/calling-c-and-fortran-code/
      z, zlocal = gensym(string(func)), gensym()
      eval(current_module(),:(global $z = C_NULL))
      z = esc(z)
      quote
          let $zlocal::Ptr{Void} = $z::Ptr{Void}
              if $zlocal == C_NULL
                 $zlocal = Libdl.dlsym($(esc(lib))::Ptr{Void}, $(esc(func)))
                 global $z = $zlocal
              end
              $zlocal
          end
      end
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
    rotationMode::Int32     # Algorithm used for free rotation of rigid bodies
    Data::Ptr{Void}         # Pointer to EmDee system data
  end

  typealias VecOrPtr Union{Vector{Float64},Ptr{Void}}
  typealias MatOrPtr Union{Matrix{Float64},Ptr{Void}}

#---------------------------------------------------------------------------------------------------
function system( threads::Int, layers::Int, rc::Float64, skin::Float64, N::Int,
                 types::Vector{Int}, masses::VecOrPtr )
  return ccall( @dlsym("EmDee_system", EmDee_libvar) tEmDee,
                (Int32, Int32, Float64, Float64, Int32, Ptr{Int32}, Ptr{Float64}),
                threads, layers, rc, skin, N, Vector{Int32}(types), masses )
end

function system( threads::Int, layers::Int, rc::Float64, skin::Float64, N::Int,
                 types::Ptr{Void}, masses::VecOrPtr )
  return ccall( @dlsym("EmDee_system", EmDee_libvar) tEmDee,
                (Int32, Int32, Float64, Float64, Int32, Ptr{Int32}, Ptr{Float64}),
                threads, layers, rc, skin, N, types, masses )
end
#---------------------------------------------------------------------------------------------------
  function set_pair_type( md::tEmDee, itype::Int, jtype::Int, model::Ptr{Void} )
    ccall( @dlsym("EmDee_set_pair_type", EmDee_libvar), Void,
           (tEmDee, Int32, Int32, Ptr{Void}),
           md, itype, jtype, model )
  end
#---------------------------------------------------------------------------------------------------
  function switch_model_layer( md::tEmDee, layer::Int )
    ccall( @dlsym("EmDee_switch_model_layer", EmDee_libvar), Void, (tEmDee, Int32), md, layer )
  end
#---------------------------------------------------------------------------------------------------
  function set_charges( md::tEmDee, charges::Vector{Float64} )
    ccall( @dlsym("EmDee_set_charges", EmDee_libvar), Void, (tEmDee, Ptr{Float64}), md, charges )
  end
#---------------------------------------------------------------------------------------------------
  function add_bond( md::tEmDee, i::Int, j::Int, model::Ptr{Void} )
    ccall( @dlsym("EmDee_add_bond", EmDee_libvar), Void,
           (tEmDee, Int32, Int32, Ptr{Void}), md, i, j, model )
  end
#---------------------------------------------------------------------------------------------------
  function add_angle( md::tEmDee, i::Int, j::Int, k::Int, model::Ptr{Void} )
    ccall( @dlsym("EmDee_add_angle", EmDee_libvar), Void,
           (tEmDee, Int32, Int32, Int32, Ptr{Void}), md, i, j, k, model )
  end
#---------------------------------------------------------------------------------------------------
  function add_dihedral( md::tEmDee, i::Int, j::Int, k::Int, l::Int, model::Ptr{Void} )
    ccall( @dlsym("EmDee_add_dihedral", EmDee_libvar), Void,
           (tEmDee, Int32, Int32, Int32, Int32, Ptr{Void}), md, i, j, k, l, model )
  end
#---------------------------------------------------------------------------------------------------
  function add_rigid_body( md::tEmDee, N::Int, indexes::Vector{Int} )
    ccall( @dlsym("EmDee_add_rigid_body", EmDee_libvar), Void,
           (tEmDee, Int32, Ptr{Int32}), md, N, indexes )
  end
#---------------------------------------------------------------------------------------------------
  function ignore_pair( md::tEmDee, i::Int, j::Int )
    ccall( @dlsym("EmDee_ignore_pair", EmDee_libvar), Void, (tEmDee, Int32, Int32), md, i, j )
  end
#---------------------------------------------------------------------------------------------------
  function upload( md::tEmDee, Lbox::VecOrPtr, coords::MatOrPtr,
                   momenta::MatOrPtr, forces::MatOrPtr )
    ccall( @dlsym("EmDee_upload", EmDee_libvar), Void,
           (Ptr{Void}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}),
           pointer_from_objref(md), Lbox, coords, momenta, forces )
  end
#---------------------------------------------------------------------------------------------------
  function download( md::tEmDee, Lbox::VecOrPtr, coords::MatOrPtr,
                     momenta::MatOrPtr, forces::MatOrPtr )
    ccall( @dlsym("EmDee_download", EmDee_libvar), Void,
           (tEmDee, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}),
           md, Lbox, coords, momenta, forces )
  end
#---------------------------------------------------------------------------------------------------
  function random_momenta( md::tEmDee, kT::Float64, adjust::Int, seed::Int )
    ccall( @dlsym("EmDee_random_momenta", EmDee_libvar), Void,
           (Ptr{Void}, Float64, Int32, Int32),
           pointer_from_objref(md), kT, adjust, seed )
  end
#---------------------------------------------------------------------------------------------------
  function boost( md::tEmDee, lambda::Float64, alpha::Float64, dt::Float64,
                        translation::Int, rotation::Int )
    ccall( @dlsym("EmDee_boost", EmDee_libvar), Void,
           (Ptr{Void}, Float64, Float64, Float64, Int32, Int32),
           pointer_from_objref(md), lambda, alpha, dt, translation, rotation )
  end
#---------------------------------------------------------------------------------------------------
  function move( md::tEmDee, lambda::Float64, alpha::Float64, dt::Float64 )
    ccall( @dlsym("EmDee_move", EmDee_libvar), Void,
           (Ptr{Void}, Float64, Float64, Float64),
           pointer_from_objref(md), lambda, alpha, dt )
  end
#---------------------------------------------------------------------------------------------------
#                                            M O D E L S
#---------------------------------------------------------------------------------------------------
  function pair_lj( ɛ::Float64, σ::Float64 )
    return ccall( @dlsym("EmDee_pair_lj", EmDee_libvar), Ptr{Void}, (Float64, Float64), ɛ, σ )
  end
#---------------------------------------------------------------------------------------------------
  function pair_lj_sf( ɛ::Number, σ::Number, rc::Number )
    return ccall( @dlsym("EmDee_pair_lj_sf", EmDee_libvar), Ptr{Void},
                  (Float64, Float64, Float64), ɛ, σ, rc )
  end
#---------------------------------------------------------------------------------------------------
  function bond_harmonic( k::Number, r₀::Number )
    return ccall( @dlsym("EmDee_bond_harmonic", EmDee_libvar), Ptr{Void}, (Float64, Float64), k, r₀ )
  end
#---------------------------------------------------------------------------------------------------
  function bond_morse( D::Number, α::Number, r₀::Number )
    return ccall( @dlsym("EmDee_bond_morse", EmDee_libvar), Ptr{Void},
                  (Float64, Float64, Float64), D, α, r₀ )
  end
#---------------------------------------------------------------------------------------------------
  function angle_harmonic( k::Number, θ₀::Number )
    return ccall( @dlsym("EmDee_angle_harmonic", EmDee_libvar), Ptr{Void}, (Float64, Float64), k, θ₀ )
  end
#---------------------------------------------------------------------------------------------------
  function dihedral_harmonic( k::Number, ϕ₀::Number )
    return ccall( @dlsym("EmDee_dihedral_harmonic", EmDee_libvar), Ptr{Void}, (Float64, Float64), k, ϕ₀ )
  end
#---------------------------------------------------------------------------------------------------
end
