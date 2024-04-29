module constants

   use iso_fortran_env

   implicit none

   integer, parameter :: r8 = real64
   integer, parameter :: i8 = int64

   real(r8), parameter :: MASS_PROTON = 1.67262192369d-27
   real(r8), parameter :: MU_0 = 1.25663706212d-6
   real(r8), parameter :: PI = 4._r8 * atan(1._r8)
   real(r8), parameter :: TWOPI = 8._r8 * atan(1._r8)

end module constants
