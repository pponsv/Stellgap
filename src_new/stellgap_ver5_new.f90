program tae_continua

   !  Serial version

   use kind_spec
   use globals
   use fourier_lib
   use fourier_lib_convolve
   use postprocess
   use fitpack
   use input
   use output
   use helper
   use eig_solver

   implicit none




   !!!!!!!  PRUEBAS

   real(r8), allocatable :: tmp(:), f2(:)
   real(r8) :: f2_avg
   integer :: tmp1
   type(timer_) :: timer

   real(r8), allocatable, dimension(:,:,:) :: f1_big, f2_big

   !!!!!!!  END PRUEBAS

   !  Dummy integers for loops

   integer :: ir, istat, i, j
   integer :: ni, nj, mi, mj, ieq, meq, neq


   real(r8) :: ccci, scsi, f1_avg, f3_avg

   real(r8), allocatable, dimension(:) :: f1_nm, f3a_nm, f3b_nm, f3c_nm
   ! real(r8), allocatable, dimension(:) :: eig_vect
   real(r8), allocatable, dimension(:,:) :: f1, f3a, f3b, f3c


   !  START PROGRAM

   ! cyl = .false.
   ! lrfp = .false.
   ! jobz = 'V'

   !  Read command-line arguments
   call read_args

   !  Read plasma.dat
   call read_plasma_dat

   !  Read fourier.dat
   call read_fourier_dat

   allocate (ion_density(ir_fine_scl), stat = istat)
   allocate (iota_r(ir_fine_scl), stat = istat)
   allocate (iota_r_inv(ir_fine_scl), stat = istat)
   allocate (iotac(irads), stat = istat)

   !    Generate Fourier arrays - TODO CHANGE
   call trig_array

   call convolution_array


   !
   !   Open files for output
   !
   !      write(*,*) trim(adjustl(outfile))

   write (*, fmt='(10(A10,3x,i4,/))') "izt:", izt, "ith:", ith, "irads:", irads, &
      "irads3:", 3*irads,"mnmx:", mnmx, "ith*izt:", ith*izt, "mn_col:", mn_col, &
      "ir_fine:", ir_fine_scl, "mpol:", mpol, "ntor:", ntor


   !    Boozer coordinates input - new ae-mode-structure input
   call read_tae_data_boozer

   !    Record equilibrium and eigenfunction Fourier mode list
   call write_modes

   !
   !   Make spline fits and fill in fine_radius_scale arrays - TODO - MAKE FUNCTION
   !
   iota_r = interp_wrap_rename(rho, iotac, rho_fine)
   iota_r_inv = interp_wrap_rename(rho, 1._r8/iotac, rho_fine) ! for lrfp .eq. .true.

   select case (ion_profile)
    case (0)
      ion_density = (iota_r / iotac(1))**2
    case (1)
      ion_density = poly_eval(rho_fine, nion)
    case (2)
      ion_density = 1._r8
    case (3)
      ion_density = (1. - aion * (rho_fine**bion))**cion
   end select

   mu0_rho_ion = MU_0 * mass_ion * ion_density_0 * scale_khz * ion_density

   call write_ion_profile

   !  Interpolate field, jacobian and gss
   bfield_lrg = interp_3d_s(bfield, rho, rho_fine)
   rjacob_lrg = interp_3d_s(rjacob, rho, rho_fine)
   gsssup_lrg = interp_3d_s(gsssup, rho, rho_fine)

   !  Allocate arrays for fourier expansions etc

   ! allocate (eig_vect(mn_col), stat = istat)
   allocate (f1_nm(mnmx), stat = istat)
   allocate (f3a_nm(mnmx), stat = istat)
   allocate (f3b_nm(mnmx), stat = istat)
   allocate (f3c_nm(mnmx), stat = istat)
   allocate (f1(izt, ith), stat = istat)
   allocate (f3a(izt, ith), stat = istat)
   allocate (f3b(izt, ith), stat = istat)
   allocate (f3c(izt, ith), stat = istat)
   allocate (f2(mnmx), stat = istat)   !PRUEBA


   open (unit = 21, file = "alfven_spec", status = "unknown")

   call initialize_solver

   f1_big = gsssup_lrg * rjacob_lrg / (bfield_lrg**2)
   f2_big = gsssup_lrg / (rjacob_lrg * (bfield_lrg**2))

   !  Main loop
   do ir = 1, ir_fine_scl

      !  Print percentage completed
      if (modulo(ir, ir_fine_scl/10) .eq. 0) write(*,fmt='(i4,"%")', advance="no") 100*ir/ir_fine_scl


      ! Make arrays to be expanded (eqs. 6, 8 of the paper)
      ! f1 = gsssup_lrg(:, :, ir) * rjacob_lrg(:, :, ir) / (bfield_lrg(:, :, ir)**2)
      ! print *, maxval(abs(f1 - f1_big(:,:,ir)))
      ! f1 = f1_big(:,:,ir)
      ! f1_avg = sum(f1_big(:,:,ir)) / real(izt*ith)
      ! if (.not. lrfp) then
      select case (lrfp)
       case (.true.)
         ! f3a = gsssup_lrg(:, :, ir) / (rjacob_lrg(:, :, ir) * (bfield_lrg(:, :, ir)**2))
         f3a = f2_big(:,:,ir)
         f3b = iota_r_inv(ir) * f3a
         f3c = iota_r_inv(ir) * f3b
         ! f3_avg = sum(f3a) / real(izt * ith) !  Unused here
       case (.false.)
         ! f3c = gsssup_lrg(:, :, ir) / (rjacob_lrg(:, :, ir) * (bfield_lrg(:, :, ir)**2))
         f3c = f2_big(:,:,ir)
         f3b = iota_r(ir) * f3c
         f3a = iota_r(ir) * f3b
         ! f3_avg = sum(f3c) / real(izt * ith) !  Unused here
         ! else if (lrfp) then
      end select

      !  Inútil?
      if (cyl) then
         if (.not. lrfp) then
            f1 = f1_avg
         else if (lrfp) then
            f1 = f1_avg
         end if
      end if

      !     Generate Fourier spectra of above gridded coefficients
      !
      f1_nm  = toFourier_new(f1_big(:,:,ir), trigtype='c') !  Si le quitas el transpose desaparecen los gaps??
      f3a_nm = toFourier_new(f3a, trigtype='c')
      f3b_nm = toFourier_new(f3b, trigtype='c')
      f3c_nm = toFourier_new(f3c, trigtype='c')
      f2 = toFourier_new(f2_big(:,:,ir), trigtype='c')

      ! print *, maxval(abs(f2 - f3c_nm))

      !     Write out coefficient spectra if half way out in flux
      !
      if (ir .eq. ir_fine_scl / 2) call write_coef_arrays(f1_nm, f3a_nm, f3b_nm, f3c_nm)  ! Why?

      !
      !     Build A and B matrices
      !
      do i = 1, mn_col
         do j = 1, mn_col
            bmat(i, j) = 0
            amat(i, j) = 0
            if ((i .lt. j) .and. ipos_def_sym)  cycle     !i.e. symmetric storage mode: only need bottom half
            ni = in_col(i)
            nj = in_col(j)
            mi = im_col(i)
            mj = im_col(j)
            do ieq = 1, mnmx
               meq = nint(m_fourier(ieq))
               neq = nint(n_fourier(ieq))
               call ccc_convolve(ccci, mi, ni, mj, nj, meq, neq)
               call scs_convolve(scsi, mi, ni, mj, nj, meq, neq)
               bmat(i, j) = bmat(i, j) - ccci * f1_nm(ieq) * mu0_rho_ion(ir)
               amat(i, j) = amat(i, j)&
                  - (scsi * (f3a_nm(ieq) * rm_col(i) * rm_col(j)&
                  - f3b_nm(ieq) * rm_col(j) * rn_col(i)&
                  - f3b_nm(ieq) * rn_col(j) * rm_col(i)&
                  + f3c_nm(ieq) * rn_col(j) * rn_col(i)))
            end do              !ieq = 1,mnmx
         end do                 !j=1,mn_col
      end do                  !i=1,mn_col

      call eigenvalue_solver

      call write_output(ir)

   end do       !ir=1,ir_fine_scl
   !
   !
   close (unit = 21)


   !  Write and deallocate
   write (*, '(/,"modes = ",i5,2x,"no. of radial points = ",i5)') mn_col, ir_fine_scl

   call write_data_post

   call post_process

   call write_nc_all

contains



!    function make_amatrix() result(amat)

!    end function make_amatrix

!    function make_bmatrix() result(bmat)

!    end function make_bmatrix
end program tae_continua
