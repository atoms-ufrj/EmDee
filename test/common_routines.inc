!---------------------------------------------------------------------------------------------------
  subroutine command_line_arguments( filename, threads )
    character(*), intent(out) :: filename
    integer,      intent(out) :: threads

    integer :: argcount
    character(256) :: line

    argcount = command_argument_count()
    if (argcount == 1) then
      threads = 1
      call get_command_argument( 1, line )
    else if (argcount == 2) then
      call get_command_argument( 1, line )
      read(line,*) threads
      call get_command_argument( 2, line )
    else
      write(0,'("Usage: testfortran [number-of-threads] input-file")')
      stop
    end if
    filename = trim(line)
  end subroutine command_line_arguments
!---------------------------------------------------------------------------------------------------
  subroutine read_data( file )
    character(*), intent(in) :: file
    integer  :: inp, i, nseeds, seed
    open( newunit = inp, file = file, status = "old" )
    read(inp,*); read(inp,*) N
    read(inp,*); read(inp,*) Rc
    read(inp,*); read(inp,*) Rs
    read(inp,*); read(inp,*) seed
    read(inp,*); read(inp,*) Dt
    read(inp,*); read(inp,*) Nsteps
    read(inp,*); read(inp,*) Nprop
    read(inp,*); read(inp,*) rho
    read(inp,*); read(inp,*) Temp
    close(inp)
    call random_seed( size = nseeds )
    call random_seed( put = seed + 37*[(i-1,i=1,nseeds)] )
    Rc2 = Rc**2
    L = (N/rho)**(1.0_8/3.0_8)
    Dt_2 = 0.5_8*Dt
    allocate( R(3,N), V(3,N) )
  end subroutine read_data
!---------------------------------------------------------------------------------------------------
  real(rb) function random_normal()
    real(rb) :: uni(2)
    call random_number( uni )
    random_normal = sqrt(-2.0_rb*log(uni(1)))*cos(6.283185307179586_rb*uni(2))
  end function random_normal
!---------------------------------------------------------------------------------------------------
  subroutine create_configuration
    integer :: Nd, m, ind, i, j, k
    real(8) :: Vcm(3)
    Nd = ceiling(real(N,8)**(1.0_8/3.0_8))
    do ind = 1, N
      m = ind - 1
      k = m/(Nd*Nd)
      j = (m - k*Nd*Nd)/Nd
      i = m - j*Nd - k*Nd*Nd
      R(:,ind) = (L/real(Nd,8))*(real([i,j,k],8) + 0.5_8)
      V(:,ind) = [random_normal(), random_normal(), random_normal()]
    end do
    Vcm = sum(V,2)/N
    forall (i=1:N) V(:,i) = V(:,i) - Vcm
    V = sqrt(Temp*(3*N-3)/sum(V*V))*V
    Vcm = sum(V,2)/N
  end subroutine create_configuration
!---------------------------------------------------------------------------------------------------
  subroutine set_charges( types, Q )
    integer,  pointer, intent(out) :: types(:)
    real(rb), pointer, intent(out) :: Q(:)

    if (mod(N,2) /= 0) stop "PLEASE ENTER AN EVEN NUMBER OF ATOMS"
    allocate( types(N), Q(N) )
    types(1:N/2) = 1
    types(N/2+1:N) = 2
    where (types == 1)
      Q = 1.0_rb
    elsewhere
      Q = -1.0_rb
    end where
    if (mod(N,2) /= 0) Q(N) = 0.0_rb
  end subroutine set_charges
!---------------------------------------------------------------------------------------------------
  subroutine run( Nsteps, Nprop )
    integer, intent(in) :: Nsteps, Nprop
    integer :: step
    print*, 0, md%Potential, md%Virial, md%Potential + md%Kinetic
    do step = 1, Nsteps
      md%options%computeProps = mod(step,Nprop) == 0
      call EmDee_boost( md, 1.0_rb, 0.0_rb, Dt_2 )
      call EmDee_move( md, 1.0_rb, 0.0_rb, Dt )
      call EmDee_boost( md, 1.0_rb, 0.0_rb, Dt_2 )
      if (mod(step,Nprop) == 0) print*, step, md%Potential, md%Virial, md%Potential + md%Kinetic
    end do
    if (Nsteps > 0) then
      print*, "neighbor list builds = ", md%builds
      print*, "pair time = ", md%pairTime, " s."
      print*, "execution time = ", md%totalTime, " s."
    end if
  end subroutine run
!---------------------------------------------------------------------------------------------------

