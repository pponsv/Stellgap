module output

   use kind_spec
   implicit none

contains
   subroutine write_modes(mnmx, mn_col, rm, rn, im_col, in_col)
      integer, intent(in) :: mnmx, mn_col, im_col(:), in_col(:)
      real(r8) :: rm(:), rn(:)
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
end module output
