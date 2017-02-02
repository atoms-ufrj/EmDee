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
    read(inp,*); read(inp,*) configFile
    read(inp,*); read(inp,*) Rc
    read(inp,*); read(inp,*) Rs
    read(inp,*); read(inp,*) seed
    read(inp,*); read(inp,*) Dt
    read(inp,*); read(inp,*) Nsteps
    read(inp,*); read(inp,*) Nprop
    read(inp,*); read(inp,*) Temp
    close(inp)
    call random_seed( size = nseeds )
    call random_seed( put = seed + 37*[(i-1,i=1,nseeds)] )
    Rc2 = Rc**2
    Dt_2 = 0.5_8*Dt
  end subroutine read_data
!---------------------------------------------------------------------------------------------------
  subroutine initialize_system( nlayers, pair )
    integer,     intent(in) :: nlayers
    type(c_ptr), intent(in) :: pair(ntypes)

    integer :: i

    md = EmDee_system( threads, nlayers, Rc, Rs, N, c_loc(atomType), c_loc(mass) )
    do i = 1, ntypes
      call EmDee_set_pair_model( md, i, i, pair(i), kCoul )
    end do
    call EmDee_upload( md, "charges"//c_null_char, c_loc(Q) )
    call EmDee_upload( md, "box"//c_null_char, c_loc(L) )
    call EmDee_upload( md, "coordinates"//c_null_char, c_loc(R(1,1)) )

  end subroutine initialize_system
!---------------------------------------------------------------------------------------------------
  real(rb) function random_normal()
    real(rb) :: uni(2)

    call random_number( uni )
    random_normal = sqrt(-2.0_rb*log(uni(1)))*cos(6.283185307179586_rb*uni(2))

  end function random_normal
!---------------------------------------------------------------------------------------------------
  subroutine set_charges( types, Q )
    integer,  intent(in)  :: types(:)
    real(rb), intent(out) :: Q(:)

    if (mod(N,2) /= 0) stop "PLEASE ENTER AN EVEN NUMBER OF ATOMS"
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
