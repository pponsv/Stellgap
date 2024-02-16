module output

   use constants

   implicit none

contains

   subroutine write_surf_area_elements(surf_compute)
      ! use globals, only: nfp, izeta, itheta
      logical, intent(in) :: surf_compute
      integer :: io_surf_area_elements

      if (.not. surf_compute) return

      open (unit = io_surf_area_elements, &
         file = "surf_area_elements", &
         status = "unknown")

      ! write (io_surf_area_elements, *) nfp, izeta, itheta

   end subroutine write_surf_area_elements
end module output
