
#---------------------------------------------------------------------------------------------------

function run(nthreads, file)

  N, Rc, Rs, seed, Dt, Nsteps, Nprop, rho, Temp = read_data(file)

  L = fill((N/rho)^(1.0/3.0),3)
  Dt_2 = 0.5*Dt

  md = EmDee.system( nthreads, 1, Rc, Rs, N, fill(1,N), fill(1.0,N), fill(0,N) )
  lj = EmDee.pair_lj_sf( 1.0, 1.0 )
  EmDee.set_pair_model( md, 1, 1, lj, 1.0 )

  R = generate_configuration( N, L[1] )

  EmDee.upload( md, "box", L )
  EmDee.upload( md, "coordinates", R )
  EmDee.random_momenta( md, Temp, 1, seed )

  println(0, " ", md.Potential, " ", md.Virial)
  for step = 1:Nsteps
    EmDee.boost( md, Dt_2 )
    EmDee.move( md, Dt )
    EmDee.boost( md, Dt_2 )
    mod(step,50) == 0 && println(step, " ", md.Potential, " ", md.Virial)
  end
  println("Pair time  = ",md.pairTime," s.")
  println("Total time = ",md.totalTime," s.")

end

#---------------------------------------------------------------------------------------------------

function read_data(file)
  f = open(file)
  readline(f); N = parse(Int,chomp(readline(f)))
  readline(f); Rc = parse(Float64,chomp(readline(f)))
  readline(f); Rs = parse(Float64,chomp(readline(f)))
  readline(f); seed = parse(Int,chomp(readline(f)))
  readline(f); Dt = parse(Float64,chomp(readline(f)))
  readline(f); Nsteps = parse(Int,chomp(readline(f)))
  readline(f); Nprop = parse(Int,chomp(readline(f)))
  readline(f); rho = parse(Float64,chomp(readline(f)))
  readline(f); Temp = parse(Float64,chomp(readline(f)))
  close(f)
  return N, Rc, Rs, seed, Dt, Nsteps, Nprop, rho, Temp
end

#---------------------------------------------------------------------------------------------------

function generate_configuration( N, L )
  R = Array(Float64,3,N)
  Nd = ceil(Int,N^(1.0/3.0))
  for ind = 1:N
    m = ind - 1
    k = div(m,Nd*Nd)
    j = div(m - k*Nd*Nd,Nd)
    i = m - j*Nd - k*Nd*Nd
    R[:,ind] = (L/Nd)*([i j k]' + 0.5)
  end
  return R
end

#---------------------------------------------------------------------------------------------------

push!(Base.DL_LOAD_PATH,string(DIR,"/lib"))
include(string(DIR,"/include/libemdee.jl"))
if length(ARGS) == 1
  run(1,ARGS[1])
elseif length(ARGS) == 2
  run(parse(Int,ARGS[1]),ARGS[2])
else
  println("Usage: testjulia [number-of-threads] input-file")
end
pop!(Base.DL_LOAD_PATH)

