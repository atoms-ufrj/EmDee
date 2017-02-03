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

    integer  :: inp, i, nseeds

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
  subroutine run( Nsteps, Nprop )
    integer, intent(in) :: Nsteps, Nprop

    integer :: step

    print*, 0, mvv2e*md%Potential, mvv2e*md%Virial, mvv2e*(md%Potential + md%Kinetic)
    do step = 1, Nsteps
      md%options%computeProps = mod(step,Nprop) == 0
      call EmDee_boost( md, 1.0_rb, 0.0_rb, Dt_2 )
      call EmDee_move( md, 1.0_rb, 0.0_rb, Dt )
      call EmDee_boost( md, 1.0_rb, 0.0_rb, Dt_2 )
      if (mod(step,Nprop) == 0) then
        print*, step, mvv2e*md%Potential, mvv2e*md%Virial, mvv2e*(md%Potential + md%Kinetic)
      end if
    end do
    if (Nsteps > 0) then
      print*, "neighbor list builds = ", md%builds
      print*, "pair time = ", md%pairTime, " s."
      print*, "execution time = ", md%totalTime, " s."
    end if

  end subroutine run
!---------------------------------------------------------------------------------------------------
