module mConfig

use iso_c_binding

implicit none

integer, parameter :: rb = c_double
integer, parameter :: ib = c_int
integer, parameter :: lb = c_bool

integer, parameter :: sl = 256
character(*), parameter :: cform = '(A256)'

integer(ib) :: ntypes
integer(ib) :: NperMol
integer(ib) :: N
real(rb)    :: mvv2e         ! Energy conversion factor
real(rb)    :: Pconv         ! Pressure conversion factor
real(rb)    :: kB            ! Boltzmann's constant
real(rb)    :: kCoul         ! Coulomb's constant

real(rb), allocatable :: epsilon(:)
real(rb), allocatable :: sigma(:)
real(rb), allocatable :: typeCharge(:)

real(rb), target :: L

integer(ib), allocatable, target :: atomType(:)
real(rb),    allocatable, target :: mass(:)
real(rb),    allocatable, target :: R(:,:)
real(rb),    allocatable, target :: P(:,:)
real(rb),    allocatable, target :: Q(:)

contains

!---------------------------------------------------------------------------------------------------

  subroutine read_configuration( configFile )
    character(*), intent(in) :: configFile

    integer(ib) :: inp, i, id, Nmol
    character(sl) :: Description

    open( newunit = inp, file = configFile, status = "old" )
    read(inp,*); read(inp,cform) Description
    write(*,'(A)') trim(Description)
    read(inp,*); read(inp,*) mvv2e
    read(inp,*); read(inp,*) Pconv
    read(inp,*); read(inp,*) kB
    read(inp,*); read(inp,*) kCoul
    read(inp,*); read(inp,*) ntypes
    allocate( mass(ntypes), epsilon(ntypes), sigma(ntypes), typeCharge(ntypes) )
    read(inp,*)
    do i = 1, ntypes
      read(inp,*) id, mass(id), epsilon(id), sigma(id), typeCharge(id)
    end do
    epsilon = kB*epsilon
    read(inp,*); read(inp,*) L
    read(inp,*); read(inp,*) NperMol
    read(inp,*); read(inp,*) Nmol
    N = NperMol*Nmol
    allocate( R(3,N), P(3,N), atomType(N), Q(N) )
    read(inp,*)
    do i = 1, N
      read(inp,*) id, R(:,id), atomType(id)
    end do
    close(inp)
    P = 0.0_rb
    Q = typeCharge(atomType)

  end subroutine read_configuration

!---------------------------------------------------------------------------------------------------

end module mConfig
