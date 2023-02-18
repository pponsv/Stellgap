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

   real(r8), allocatable :: tmp(:), f2(:,:)
   real(r8) :: f2_avg

   !!!!!!!  END PRUEBAS

   !  Dummy integers for loops

   integer :: ir, istat, i, j, mn
   integer :: ni, nj, mi, mj, ieq, meq, neq
   integer :: j_max_index, m_emax, n_emax


   real(r8) :: eig_max, ccci, scsi, f1_avg, f3_avg

   real(r8), allocatable, dimension(:) :: f1_nm, f3a_nm, f3b_nm, f3c_nm
   real(r8), allocatable, dimension(:) :: eig_vect
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

   allocate (eig_vect(mn_col), stat = istat)

   allocate (f1_nm(mnmx), stat = istat)
   allocate (f3a_nm(mnmx), stat = istat)
   allocate (f3b_nm(mnmx), stat = istat)
   allocate (f3c_nm(mnmx), stat = istat)
   allocate (f1(izt, ith), stat = istat)
   allocate (f3a(izt, ith), stat = istat)
   allocate (f3b(izt, ith), stat = istat)
   allocate (f3c(izt, ith), stat = istat)


   allocate (f2(izt, ith), stat = istat)

   !
   !   Open files for output
   !
   !      write(*,*) trim(adjustl(outfile))

   write (*, *) izt, ith, irads, mnmx, ith*izt, mn_col


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

   open (unit = 21, file = "alfven_spec", status = "unknown")

   call initialize_solver

   !  Main loop
   do ir = (1), (ir_fine_scl)

      if (modulo(ir, ir_fine_scl/10) .eq. 0) print *, 100*ir/ir_fine_scl, '%, ', rho_fine(ir)

      ! f1_avg = 0.; f3_avg = 0.
      f1 = gsssup_lrg(:, :, ir) * rjacob_lrg(:, :, ir) / (bfield_lrg(:, :, ir)**2)
      f1_avg = sum(f1) / real(izt*ith)
      if (.not. lrfp) then
         f3c = gsssup_lrg(:, :, ir) / (rjacob_lrg(:, :, ir) * (bfield_lrg(:, :, ir)**2))
         f3b = iota_r(ir) * f3c
         f3a = iota_r(ir) * f3b
         f3_avg = sum(f3c) / real(izt * ith)
      else if (lrfp) then
         f3a = gsssup_lrg(:, :, ir) / (rjacob_lrg(:, :, ir) * (bfield_lrg(:, :, ir)**2))
         f3b = iota_r_inv(ir) * f3a
         f3c = iota_r_inv(ir) * f3b
         f3_avg = sum(f3a) / real(izt * ith)
      end if

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
      f1_nm  = toFourier_new(f1, trigtype='c') !  Si le quitas el transpose desaparecen los gaps??
      f3a_nm = toFourier_new(f3a, trigtype='c')
      f3b_nm = toFourier_new(f3b, trigtype='c')
      f3c_nm = toFourier_new(f3c, trigtype='c')

      !     Write out coefficient spectra if half way out in flux
      !
      if (ir .eq. ir_fine_scl / 2) call write_coef_arrays(f1_nm, f3a_nm, f3b_nm, f3c_nm)
         ! open (unit = 8, file = "coef_arrays", status = "unknown")
         ! do mn = 1, mnmx
         !    write (8, '(f6.1,2x,f6.1,4(2x,e15.7))') rm(mn), rn(mn), &
         !       f1_nm(mn), f3a_nm(mn), f3b_nm(mn), f3c_nm(mn)
         ! end do
         ! close (unit = 8)
      ! end if
      !
      !     Build A and B matrices
      !
      do i = 1, mn_col
         do j = 1, mn_col
            bmat(i, j) = 0
            amat(i, j) = 0
            if (i .lt. j .and. ipos_def_sym) cycle     !i.e. symmetric storage mode: only need bottom half
            ni = in_col(i)
            nj = in_col(j)
            mi = im_col(i)
            mj = im_col(j)
            do ieq = 1, mnmx
               meq = rm(ieq)
               neq = rn(ieq)
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


      !     Call matrix eigenvalue solver
      call eigenvalue_solver

      do i = 1, mn_col
         if (iopt .eq. 1) then
            do j = 1, mn_col
               if (ipos_def_sym .and. jobz .eq. 'V' .and. .not. subset_eq) then
                  eig_vect(j) = abs(amat(j, i))
               end if
               if (.not. ipos_def_sym) eig_vect(j) = abs(vr(j, i))
            end do

            eig_max = -1.d+30
            do j = 1, mn_col
               if (eig_vect(j) .gt. eig_max) then
                  eig_max = eig_vect(j)
                  j_max_index = j
               end if
            end do
            m_emax = im_col(j_max_index)
            n_emax = in_col(j_max_index)

            !       if(ir .eq. 5 .and. i .eq. mn_col/2)
            !     >   write(*,'(e15.7,2x,i4,2(2x,i4))') eig_max,m_emax,n_emax,
            !     >    j_max_index

            if (ipos_def_sym) then
               write (21, '(2(e15.7,2x),i4,2x,i4)') rho_fine(ir), sqrt(abs(omega(i))), m_emax, n_emax
            else if (.NOT. ipos_def_sym) then
               write (21, '(4(e15.7,2x),i4,2x,i4)') rho_fine(ir), alfr(i), alfi(i), beta(i), m_emax, n_emax
            end if
         else if (iopt .eq. 0) then
            if (ipos_def_sym) then
               write (21, '(e15.7,2x,e15.7)') rho_fine(ir), sqrt(abs(omega(i)))
            else if (.NOT. ipos_def_sym) then
               write (21, '(e15.7,3(2x,e15.7))') rho_fine(ir), alfr(i), alfi(i), beta(i)
            end if
         end if
      end do          !do i=1,mn_col

   end do       !ir=1,ir_fine_scl
   !
   !
   close (unit = 21)


   !  Write and deallocate
   write (*, '("modes = ",i5,2x,"no. of radial points = ",i5)') mn_col, ir_fine_scl

   ! call trg_deallocate

   call write_data_post

   call post_process

   stop
end program tae_continua
