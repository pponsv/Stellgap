module helper

   use constants, only: r8


   implicit none
   type :: timer_
      real(r8) :: t0, t1, dt, t0_cum, t1_cum, dt_cum
   contains
      procedure :: tic => tic
      procedure :: toc  => toc
   end type timer_

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


   function interpolate_3d_s(in, s, s_fine) result(out)
      ! use fitpack, only: curv1, curv2
      real(r8), intent(in) :: in(:,:,:), s(:), s_fine(:)

      real(r8) :: out(size(in, dim=1), size(in, dim=2), size(s_fine))

      integer :: i, j

      do i = 1, size(in, dim=1)
         do j = 1, size(in, dim=2)
            out(i, j, :) = interpolate_using_spline(s, in(i, j, :), s_fine)
         end do
      end do
   end function interpolate_3d_s


   function interpolate_using_spline(x, y, x_int) result(out)
      use fitpack, only: curv1, curv2

      real(r8), intent(in) :: x(:), y(size(x)), x_int(:)
      real(r8) :: out(size(x_int))

      real(r8) :: ypi(size(x)), temp(size(x))
      integer :: ierr, i

      call curv1(n=size(x), x=x, y=y, slp1=0._r8, slpn=0._r8, islpsw=3, &
         yp=ypi, temp=temp, sigma=0._r8, ierr=ierr)
      do i = 1, size(x_int)
         out(i) = curv2(x_int(i), size(x), x, y, ypi, sigma=0._r8)
      end do
      if (ierr .ne. 0) write (*, '("spline error 1",i3)') ierr

   end function interpolate_using_spline


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
      !  Makes a flattened meshgrid from in1 and in2

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

   subroutine tic(t, text)
      class(timer_) :: t
      character(len=*), optional :: text
      integer   :: count, count_rate

      if (present(text)) print *, text
      call cpu_time(t%t0_cum)
      call system_clock(count, count_rate)
      t%t0 = real(count, kind=r8) / count_rate
   end subroutine tic

   subroutine toc(t)
      class(timer_) :: t
      integer   :: count, count_rate
      call cpu_time(t%t1_cum)
      t%dt_cum = t%t1_cum - t%t0_cum
      call system_clock(count, count_rate)
      t%t1 = real(count, kind=r8) / count_rate
      t%dt = t%t1 - t%t0
      print '(A, g10.5, A, g10.5, A)', &
         "Time elapsed: ", t%dt, "s (wall)  ", t%dt_cum, "s (cum)"
   end subroutine toc


end module helper
