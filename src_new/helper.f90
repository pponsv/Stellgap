module helper
   use kind_spec
   implicit none

contains

   function poly_eval(x, coefs) result(out)
      real(r8), intent(in) :: x(:), coefs(:)
      real(r8) :: out(size(x))
      integer :: lenx, ncoefs, i

      lenx = size(x)
      ncoefs = size(coefs)

      out = 0.
      do i=ncoefs, 1, -1
         out = out*x + coefs(i)
      end do

   end function poly_eval

end module helper
