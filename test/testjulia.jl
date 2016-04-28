function main()

N, Rc, Rs, seed, Dt, Nsteps, Nprop, rho, Temp = read_data()

L = (N/rho)^(1.0/3.0)
Dt_2 = 0.5*Dt

R, V = generate_configuration( seed, N, L, Temp )



end

function read_data()
  readline(); N = parse(Int,chomp(readline()))
  readline(); Rc = parse(Float64,chomp(readline()))
  readline(); Rs = parse(Float64,chomp(readline()))
  readline(); seed = parse(Int,chomp(readline()))
  readline(); Dt = parse(Float64,chomp(readline()))
  readline(); Nsteps = parse(Int,chomp(readline()))
  readline(); Nprop = parse(Int,chomp(readline()))
  readline(); rho = parse(Float64,chomp(readline()))
  readline(); Temp = parse(Float64,chomp(readline()))
  return N, Rc, Rs, seed, Dt, Nsteps, Nprop, rho, Temp
end

function generate_configuration( seed, N, L, Temp )
  R = Array(Float64,3,N)
  V = Array(Float64,3,N)
  RNG = MersenneTwister(seed)
  Nd = ceil(Int,N^(1.0/3.0))
  for ind = 1:N
    m = ind - 1
    k = div(m,Nd*Nd)
    j = div(m - k*Nd*Nd,Nd)
    i = m - j*Nd - k*Nd*Nd
    R[:,ind] = (L/Nd)*([i j k]' + 0.5)
    V[:,ind] = [randn(RNG) randn(RNG) randn(RNG)]'
  end
  Vcm = mean(V,2)
  for i = 1:N
    V[:,i] -= Vcm
  end
  V = sqrt(Temp*(3*N-3)/sum(V.*V))*V
  return R, V
end

main()

