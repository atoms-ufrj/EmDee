module c_binding_extra

use, intrinsic :: iso_c_binding

contains

!---------------------------------------------------------------------------------------------------

  type(c_ptr) function malloc_int( n, value, array )
    integer(c_int), intent(in)           :: n
    integer(c_int), intent(in), optional :: value, array(:)
    integer(c_int), pointer :: ptr(:)
    allocate( ptr(n) )
    malloc_int = c_loc(ptr)
    if (present(array)) then
      ptr = array
    else if (present(value)) then
      ptr = value
    end if
    nullify( ptr )
  end function malloc_int

!---------------------------------------------------------------------------------------------------

  type(c_ptr) function malloc_real( n, value, array )
    integer(c_int), intent(in)           :: n
    real(c_double),    intent(in), optional :: value, array(:)
    real(c_double), pointer :: ptr(:)
    allocate( ptr(n) )
    malloc_real = c_loc(ptr)
    if (present(array)) then
      ptr = array
    else if (present(value)) then
      ptr = value
    end if
    nullify( ptr )
  end function malloc_real

!---------------------------------------------------------------------------------------------------

  subroutine realloc_int( ptr, size, new_size )
    type(c_ptr), intent(inout) :: ptr
    integer(c_int), intent(inout) :: size
    integer(c_int), intent(in)    :: new_size

    integer(c_int) :: n
    integer(c_int), pointer :: old(:), new(:)

    call c_f_pointer( ptr, old, [size] )
    allocate( new(new_size) )
    n = min(size,new_size)
    new(1:n) = old(1:n)
    deallocate( old )
    ptr = c_loc(new)
    size = new_size

  end subroutine realloc_int

!---------------------------------------------------------------------------------------------------

  subroutine copy_real( from, to, n )
    type(c_ptr), intent(in) :: from, to
    integer,     intent(in) :: n
    real(c_double), pointer :: F(:), T(:)
    call c_f_pointer( from, F, [n] )
    call c_f_pointer( to, T, [n] )
    T = F
    nullify( T, F )
  end subroutine copy_real

!---------------------------------------------------------------------------------------------------

end module c_binding_extra
