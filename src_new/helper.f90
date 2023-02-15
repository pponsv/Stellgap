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

   function interp_3d_s(in, rho, rho_fine) result(out)
      use fitpack, only: curv1, curv2
      real(r8), intent(in) :: in(:,:,:), rho(:), rho_fine(:)

      real(r8) :: out(size(in, dim=1), size(in, dim=2), size(rho_fine))
      
      real(r8) :: sp_fit(size(in, dim=3))
      real(r8) :: sp1, sp2, yp(3*size(rho)), temp(3*size(rho)), sigma_spl
      integer :: i, j, irr, ispl_opt, ierr_spl, irads

      irads = size(rho)

      do i = 1, size(in, dim=1)
         do j = 1, size(in, dim=2)
            sp_fit = in(i, j, :)
            ispl_opt = 3; ierr_spl = 0; sigma_spl = 0
            call curv1(irads, rho, sp_fit, sp1, sp2, ispl_opt, yp, &
            &temp, sigma_spl, ierr_spl)
            if (ierr_spl .ne. 0) write (*, '("spline error 2",i3)') ierr_spl
            do irr = 1, size(rho_fine)
               out(i, j, irr) = curv2(rho_fine(irr), irads, rho, sp_fit, yp, sigma_spl)
            end do     !irr = 1,ir_fine_scl
         end do
      end do
   end function interp_3d_s

end module helper
