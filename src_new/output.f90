module output

   use kind_spec
   implicit none

contains

   subroutine write_modes
      use globals, only: rm, rn, mnmx, mn_col, im_col, in_col
      ! integer, intent(in) :: mnmx, mn_col, im_col(:), in_col(:)
      ! real(r8) :: rm(:), rn(:)
      integer :: i

      open (unit=22, file="modes", status="unknown")

      write (22, '("Equilibrium modes:",/)')
      write (22, '("meq    neq")')
      do i = 1, mnmx
         write (22, '(i3,3x,i3)') int(rm(i)), int(rn(i))
      end do
      write (22, '(///,"Eigenvector modes:",/)')
      write (22, '("m    n")')
      do i = 1, mn_col
         write (22, '(i3,3x,i3)') im_col(i), in_col(i)
      end do

      close (unit=22)

   end subroutine write_modes

   subroutine write_ion_profile
      use globals, only: bfavg, scale_khz, mu0_rho_ion, ir_fine_scl, ion_density_0, ion_density, &
         iota_r, rho_fine

      character(len=*), parameter :: fmt_ion_profile = "(e12.5, 3(3x, e12.5))"
      real(r8) :: va
      integer :: irr

      open (unit = 9, file = "ion_profile", status = "unknown")

      do irr=1, ir_fine_scl
         va = sqrt(bfavg**2 / (mu0_rho_ion(irr) / scale_khz))
         write (9, fmt_ion_profile) rho_fine(irr), ion_density_0 * ion_density(irr), iota_r(irr), va
      end do

      close (unit = 9)

   end subroutine write_ion_profile

   subroutine write_data_post(iopt, mn_col, ir_fine_scl, isym_opt)
      ! use globals, only: iopt, mn_col, ir_fine_scl, isym_opt
      integer, intent(in) :: iopt, mn_col, ir_fine_scl, isym_opt

      open (unit = 7, file = "data_post", status = "unknown")
      write (7, '(i5,3(2x,i5))') iopt, mn_col, ir_fine_scl, isym_opt
      close (unit = 7)

   end subroutine write_data_post

end module output
