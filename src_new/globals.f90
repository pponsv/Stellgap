module globals
!  Module for global variables

   use kind_spec

   implicit none

   integer :: irads, ir_fine_scl, irad3

   logical :: lrfp ! If true: use reversed field pinch settings

   real(r8), allocatable :: rho(:)



end module globals
