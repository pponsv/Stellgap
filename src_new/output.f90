module output

   use kind_spec
   use netcdf
   implicit none

contains

   subroutine write_modes
      use globals, only: m_fourier, n_fourier, mnmx, mn_col, im_col, in_col

      integer :: i

      open (unit=22, file="modes", status="unknown")

      write (22, '("Equilibrium modes:",/)')
      write (22, '("meq    neq")')
      do i = 1, mnmx
         write (22, '(i3,3x,i3)') int(m_fourier(i)), int(n_fourier(i))
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


   subroutine write_data_post
      use globals, only: iopt, mn_col, ir_fine_scl, ipos_def_sym
      integer :: isym_opt

      if (ipos_def_sym) isym_opt = 1
      if (.NOT. ipos_def_sym) isym_opt = 0

      open (unit = 7, file = "data_post", status = "unknown")
      write (7, '(i5,3(2x,i5))') iopt, mn_col, ir_fine_scl, isym_opt
      close (unit = 7)

   end subroutine write_data_post


   subroutine write_coef_arrays(f1_nm, f3a_nm, f3b_nm, f3c_nm)
      use globals, only: mnmx, m_fourier, n_fourier
      real(r8), intent(in), dimension(mnmx) :: f1_nm, f3a_nm, f3b_nm, f3c_nm
      integer :: mn

      open (unit = 8, file = "coef_arrays", status = "unknown")
      do mn = 1, mnmx
         write (8, '(f6.1,2x,f6.1,4(2x,e15.7))') m_fourier(mn), n_fourier(mn), &
            f1_nm(mn), f3a_nm(mn), f3b_nm(mn), f3c_nm(mn)
      end do
      close (unit = 8)

   end subroutine write_coef_arrays


   ! subroutine write_alfven_post(ipos_def_sym, iopt)
   !    use globals, only: rho_fine
   !    use fourier_lib, only: omega
   !    integer, intent(in) :: iopt
   !    logical, intent(in) :: ipos_def_sym

   !    if (iopt .eq. 1) then
   !       if (ipos_def_sym) then
   !          write (21, '(2(e15.7,2x),i4,2x,i4)') rho_fine(ir), sqrt(abs(omega(i))), m_emax, n_emax
   !       else if (.NOT. ipos_def_sym) then
   !          write (21, '(4(e15.7,2x),i4,2x,i4)') rho_fine(ir), alfr(i), alfi(i), beta(i), m_emax, n_emax
   !       end if
   !    else if (iopt .eq. 0) then
   !       if (ipos_def_sym) then
   !          write (21, '(e15.7,2x,e15.7)') rho_fine(ir), sqrt(abs(omega(i)))
   !       else if (.NOT. ipos_def_sym) then
   !          write (21, '(e15.7,3(2x,e15.7))') rho_fine(ir), alfr(i), alfi(i), beta(i)
   !       end if
   !    end if

   ! end subroutine write_alfven_post


   subroutine write_nc_all
      use globals, only: m_fourier, n_fourier, mnmx, mn_col, im_col, in_col, rho_fine

      integer :: ierr, file_id, id_meq, id_neq, id_meig, id_neig, id_rho

      ierr = nf90_create(path='stellgap_out.nc', cmode=NF90_CLOBBER, ncid=file_id)

      id_meq = def_1d_var(m_fourier, file_id, name='m_eq', long_name='Equilibrium modes: m', &
         units='None', xlabel='m', data_type=NF90_INT)
      id_neq = def_1d_var(n_fourier, file_id, name='n_eq', long_name='Equilibrium modes: n', &
         units='None', xlabel='n', data_type=NF90_INT)
      id_meig = def_1d_var(real(im_col, r8), file_id, name='m_eig', long_name='Eigenvector modes: m', &
         units='None', xlabel='m_eig', data_type=NF90_INT)
      id_neig = def_1d_var(real(in_col, r8), file_id, name='n_eig', long_name='Eigenvector modes: n', &
         units='None', xlabel='n_eig', data_type=NF90_INT)
      id_rho = def_1d_var(rho_fine, file_id, name='rho', long_name='Radial coordinate', &
         units='None', xlabel='rho', data_type=NF90_INT)

      ierr = nf90_enddef(file_id)

      ierr = nf90_put_var(file_id, id_meq, m_fourier)
      ierr = nf90_put_var(file_id, id_neq, n_fourier)
      ierr = nf90_put_var(file_id, id_meig, im_col)
      ierr = nf90_put_var(file_id, id_neig, in_col)
      ierr = nf90_put_var(file_id, id_rho, rho_fine)

      ierr = nf90_close(file_id)

   end subroutine write_nc_all

   function def_2d_var(a, file_id, name, long_name, units, xlabel, ylabel, data_type) result(array_id)
      real(r8), intent(in) :: a(:,:)
      integer, intent(in) :: file_id, data_type
      character(len=*), intent(in) :: name, units, xlabel, ylabel, long_name

      integer :: ierr, xdim_id, ydim_id, array_id

      ierr = nf90_def_dim(file_id, xlabel, size(a, dim=1), xdim_id)
      ierr = nf90_def_dim(file_id, ylabel, size(a, dim=2), ydim_id)
      ierr = nf90_def_var(file_id, name, data_type, [xdim_id, ydim_id], array_id)
      ierr = nf90_put_att(file_id, array_id, "units", units)
      ierr = nf90_put_att(file_id, array_id, "long_name", long_name)

   end function def_2d_var


   function def_1d_var(a, file_id, name, long_name, units, xlabel, data_type) result(array_id)
      real(r8), intent(in) :: a(:)
      integer, intent(in) :: file_id, data_type
      character(len=*), intent(in) :: name, units, xlabel, long_name

      integer :: ierr, xdim_id, array_id

      ierr = nf90_def_dim(file_id, xlabel, size(a), xdim_id)
      ierr = nf90_def_var(file_id, name, data_type, xdim_id, array_id)
      ierr = nf90_put_att(file_id, array_id, "units", units)
      ierr = nf90_put_att(file_id, array_id, "long_name", long_name)

   end function def_1d_var

end module output
