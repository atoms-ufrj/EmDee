integer(ib) :: Nsteps, Nprop, i, seed
real(rb)    :: Rc, Rs, Rc2, Temp, Dt, Dt_2

real(rb)    :: mvv2e         ! Energy conversion factor
real(rb)    :: Pconv         ! Pressure conversion factor
real(rb)    :: kB            ! Boltzmann's constant
real(rb)    :: kCoul         ! Coulomb's constant

type(tEmDee), target :: md
type(c_ptr) :: pair

integer :: threads
character(256) :: filename, configFile

real(rb) :: obtained(3), expected(3)
real(rb), parameter :: tol = 1.0e-8_rb
