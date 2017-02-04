module mConfig

use iso_c_binding

implicit none

integer, parameter :: rb = c_double
integer, parameter :: ib = c_int
integer, parameter :: lb = c_bool

integer, parameter :: sl = 256
character(*), parameter :: cform = '(A256)'

integer(ib) :: ntypes
integer(ib) :: N

real(rb), allocatable :: epsilon(:)
real(rb), allocatable :: sigma(:)

real(rb), target :: L

integer(ib), allocatable, target :: atomType(:)
integer(ib), allocatable, target :: molecule(:)
real(rb),    allocatable, target :: mass(:)
real(rb),    allocatable, target :: R(:,:)
real(rb),    allocatable, target :: P(:,:)
real(rb),    allocatable, target :: Q(:)

contains

!---------------------------------------------------------------------------------------------------

  subroutine read_configuration( configFile )
    character(*), intent(in) :: configFile

    integer(ib) :: inp, i, id
    character(sl) :: Description

    open( newunit = inp, file = configFile, status = "old" )
    read(inp,*); read(inp,cform) Description
    write(*,'(A)') trim(Description)
    read(inp,*); read(inp,*) ntypes
    allocate( mass(ntypes), epsilon(ntypes), sigma(ntypes) )
    read(inp,*)
    do i = 1, ntypes
      read(inp,*) id, mass(id), epsilon(id), sigma(id)
    end do
    read(inp,*); read(inp,*) L
    read(inp,*); read(inp,*) N
    allocate( R(3,N), P(3,N), atomType(N), molecule(N), Q(N) )
    read(inp,*)
    do i = 1, N
      read(inp,*) id, molecule(id), atomType(id), Q(id), R(:,id)
    end do
    close(inp)
    P = 0.0_rb

  end subroutine read_configuration

!---------------------------------------------------------------------------------------------------

end module mConfig
