module helper
   use kind_spec, only: r8
   implicit none

contains

   function outer(a, b) result(c)
      real(r8) :: a(:), b(:)
      real(r8) :: c(size(a), size(b))

      integer :: i, j
      do i = 1, size(b)
         do j = 1, size(a)
            c(j, i) = a(j)*b(i)
         end do
      end do
   end function outer

   function poly_eval(x, coefs) result(out)
      real(r8), intent(in) :: x(:), coefs(:)
      real(r8) :: out(size(x))
      integer :: i

      out = 0.
      do i=size(coefs), 1, -1
         out = out*x + coefs(i)
      end do

   end function poly_eval

   function interp_3d_s(in, s, s_fine) result(out)
      use fitpack, only: curv1, curv2
      real(r8), intent(in) :: in(:,:,:), s(:), s_fine(:)

      real(r8) :: out(size(in, dim=1), size(in, dim=2), size(s_fine))

      real(r8) :: sp_fit(size(in, dim=3))
      real(r8) :: sp1, sp2, yp(3*size(s)), temp(3*size(s)), sigma_spl
      integer :: i, j, irr, ispl_opt, ierr_spl, irads

      irads = size(s)

      do i = 1, size(in, dim=1)
         do j = 1, size(in, dim=2)
            sp_fit = in(i, j, :)
            ispl_opt = 3; ierr_spl = 0; sigma_spl = 0
            call curv1(irads, s, sp_fit, sp1, sp2, ispl_opt, yp, temp, sigma_spl, ierr_spl)
            if (ierr_spl .ne. 0) write (*, '("spline error 2",i3)') ierr_spl
            do irr = 1, size(s_fine)
               out(i, j, irr) = curv2(s_fine(irr), irads, s, sp_fit, yp, sigma_spl)
            end do
         end do
      end do
   end function interp_3d_s

   subroutine write_tmp(a, b)
      real(r8) :: a(:), b(:)
      integer :: i

      open(unit=99, file='tmp', status='unknown')
      do i=1, size(a)
         write(99, *) a(i), b(i)
      end do
      close(99)
   end subroutine write_tmp

   function real_linspace(xi, xf, np) result(out)

      real(r8), intent(in)    :: xi, xf
      integer, intent(in) :: np
      real(r8)    :: out(np)

      integer :: i

      do i = 1, np
         out(i) = xi + (xf - xi)/(np - 1)*(i - 1)
      end do
   end function real_linspace

   function real_linspace_nolast(xi, xf, np) result(out)
      !       Makes a linearly spaced vector between xi, xf
      !       (xf excluded) with np elements

      real(r8), intent(in)    :: xi, xf
      integer, intent(in) :: np
      real(r8)    :: out(np)

      integer :: i

      do i = 1, np
         out(i) = xi + (xf - xi)/(np)*(i - 1)
      end do
   end function real_linspace_nolast

   subroutine lingrid(in1, in2, out1, out2)
      real(r8), intent(in) :: in1(:), in2(:)
      real(r8), intent(out) :: out1(size(in1)*size(in2)), out2(size(in1)*size(in2))

      integer :: i, j, lg

      lg = 0
      do i=1, size(in1)
         do j=1, size(in2)
            lg = lg + 1
            out1(lg) = in1(i)
            out2(lg) = in2(j)
         end do
      end do

   end subroutine lingrid


end module helper
