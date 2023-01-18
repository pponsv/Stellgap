module kind_spec
!
!     Double precision
!
   implicit none
   integer, parameter :: rprec = selected_real_kind(12, 100)
   integer, parameter :: r8 = selected_real_kind(12, 100)
   integer, parameter :: iprec = selected_int_kind(8)
   integer, parameter :: i8 = selected_int_kind(8)
end module kind_spec
