module global

use, intrinsic :: iso_c_binding

implicit none

integer, parameter :: ib = c_int
integer, parameter :: rb = c_double
integer, parameter :: lb = c_bool
integer, parameter :: sl = 256

real(rb), parameter :: zero   = 0.0_rb,                 &
                       one    = 1.0_rb,                 &
                       two    = 2.0_rb,                 &
                       half   = 0.5_rb,                 &
                       third  = 0.33333333333333333_rb, &
                       fourth = 0.25_rb,                &
                       pi     = 3.14159265358979324_rb, &
                       piBy2  = 0.5_rb*pi,              &
                       invSqrt3 = one/sqrt(3.0_rb)

end module global
