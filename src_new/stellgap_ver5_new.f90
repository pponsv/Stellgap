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

   implicit none


   external dggev, dsygv, dsygvx ! LAPACK subroutines

   !!!!!!!  PRUEBAS

   real(r8), allocatable :: tmp(:)

   !!!!!!!  END PRUEBAS

   integer, parameter :: iopt = 1

   !  Dummy integers for loops

   integer :: ir, irr, il, iu, lg, istat, i, j, mn
   integer :: ispl_opt, ierr_spl
   integer :: ni, nj, mi, mj, ieq, meq, neq
   integer :: j_max_index
   integer :: m_emax, n_emax, isym_opt
   integer :: ldvl, ldvr, lwork, info, m_red
   integer, allocatable, dimension(:) :: ifail, iwork

   real(r8) :: sp1, sp2, sigma_spl
   real(r8) :: eig_max, ccci, scsi, egl, egu, abstol, f1_avg, f3_avg

   real(r8), allocatable, dimension(:) :: f1_nm, f3a_nm, f3b_nm, f3c_nm
   real(r8), allocatable, dimension(:) :: beta, eig_vect, omega, work, alfr, alfi
   real(r8), allocatable, dimension(:) :: ypi, tempi

   real(r8), allocatable, dimension(:,:) :: f1, f3a, f3b, f3c
   real(r8), allocatable, dimension(:,:) :: amat, bmat, vr, vl

   complex*16, allocatable :: alpha(:)

   character*1 jobz

   !  START PROGRAM

   cyl = .false.
   lrfp = .false.

   !  Read command-line arguments
   call read_args

   !  Read plasma.dat
   call read_plasma_dat

   !  Read fourier.dat
   call read_fourier_dat

   allocate (ion_density(ir_fine_scl),iota_r(ir_fine_scl), &
      iota_r_inv(ir_fine_scl), stat = istat)
   allocate (iotac(irads), stat = istat)
   ! allocate (yp(3*irads), stat = istat)
   allocate (ypi(3*irads), stat = istat)
   allocate (tempi(3*irads), stat = istat)

   jobz = 'V'

   if (ipos_def_sym) isym_opt = 1
   if (.NOT. ipos_def_sym) isym_opt = 0

   !    Generate Fourier arrays - TODO CHANGE
   call trig_array

   call convolution_array

   ! naux = 10 * mn_col
   scale_khz = (1.d+3 * 2 * PI)**2

   ldvl = mn_col
   ldvr = mn_col
   lwork = 20 * mn_col

   ! allocate (aux(naux), stat = istat)
   allocate (alpha(mn_col), stat = istat)
   allocate (beta(mn_col), stat = istat)
   allocate (eig_vect(mn_col), stat = istat)
   allocate (omega(mn_col), stat = istat)
   allocate (amat(mn_col, mn_col), stat = istat)
   allocate (bmat(mn_col, mn_col), stat = istat)
   allocate (ifail(mn_col), stat = istat)
   allocate (work(lwork), stat = istat)
   allocate (iwork(lwork), stat = istat)
   allocate (alfr(mn_col), stat = istat)
   allocate (alfi(mn_col), stat = istat)
   allocate (vl(mn_col, mn_col), stat = istat)
   allocate (vr(mn_col, mn_col), stat = istat)

   allocate (f1_nm(mnmx), stat = istat)
   allocate (f3a_nm(mnmx), stat = istat)
   allocate (f3b_nm(mnmx), stat = istat)
   allocate (f3c_nm(mnmx), stat = istat)
   allocate (f1(izt, ith), stat = istat)
   allocate (f3a(izt, ith), stat = istat)
   allocate (f3b(izt, ith), stat = istat)
   allocate (f3c(izt, ith), stat = istat)
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
   open (unit = 8, file = "coef_arrays", status = "unknown")

   !  Main loop
   do ir = (1), (ir_fine_scl)
      !       write(*,*) ir
      !      do ir=1,ir_fine_scl
      ! r_pt = rho(1) + real(ir - 1) * (rho(irads) - rho(1))&
      ! & / real(ir_fine_scl - 1)
      print *, ir, rho_fine(ir)
      f1_avg = 0.; f3_avg = 0.
      do i = 1, izt
         do j = 1, ith
            f1(i, j) = gsssup_lrg(i, j, ir) * rjacob_lrg(i, j, ir) / &
               (bfield_lrg(i, j, ir)**2)
            f1_avg = f1_avg + f1(i, j) / real(izt * ith)
            if (.not. lrfp) then
               f3c(i, j) = gsssup_lrg(i, j, ir) / &
                  (rjacob_lrg(i, j, ir) * (bfield_lrg(i, j, ir)**2))
               f3_avg = f3_avg + f3c(i, j) / real(izt * ith)
               f3b(i, j) = iota_r(ir) * f3c(i, j)
               f3a(i, j) = iota_r(ir) * f3b(i, j)
            else if (lrfp) then
               f3a(i, j) = gsssup_lrg(i, j, ir) / &
                  (rjacob_lrg(i, j, ir) * (bfield_lrg(i, j, ir)**2))
               f3b(i, j) = f3a(i, j) * iota_r_inv(ir)
               f3c(i, j) = f3b(i, j) * iota_r_inv(ir)
               f3_avg = f3_avg + f3a(i, j) / real(izt * ith)
            end if
         end do
      end do

      !  Inútil?
      if (cyl) then
         if (.not. lrfp) then
            f1 = f1_avg
         else if (lrfp) then
            f1 = f1_avg
         end if
      end if
      !
      !
      !     Generate Fourier spectra of above gridded coefficients
      !
      lg = 0
      do i = 1, izt
         do j = 1, ith
            lg = lg + 1
            f(lg) = f1(i, j)
         end do
      end do
      sin_type = 0; cos_type = 1
      call toFourier
      f1_nm(:) = fnm(:)
      !
      lg = 0
      do i = 1, izt
         do j = 1, ith
            lg = lg + 1
            f(lg) = f3a(i, j)
         end do
      end do
      sin_type = 0; cos_type = 1
      call toFourier
      f3a_nm(:) = fnm(:)
      !
      lg = 0
      do i = 1, izt
         do j = 1, ith
            lg = lg + 1
            f(lg) = f3b(i, j)
         end do
      end do
      sin_type = 0; cos_type = 1
      call toFourier
      f3b_nm(:) = fnm(:)
      !
      lg = 0
      do i = 1, izt
         do j = 1, ith
            lg = lg + 1
            f(lg) = f3c(i, j)
         end do
      end do
      sin_type = 0; cos_type = 1
      call toFourier
      f3c_nm(:) = fnm(:)
      !
      !
      !     Write out coefficient spectra if half way out in flux
      !
      if (ir .eq. ir_fine_scl / 2) then
         do mn = 1, mnmx
            write (8, '(f6.1,2x,f6.1,4(2x,e15.7))') rm(mn), rn(mn), &
            &f1_nm(mn), f3a_nm(mn), &
            &f3b_nm(mn), f3c_nm(mn)
         end do
      end if
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
               & - (scsi * (f3a_nm(ieq) * rm_col(i) * rm_col(j)&
               & - f3b_nm(ieq) * rm_col(j) * rn_col(i)&
               & - f3b_nm(ieq) * rn_col(j) * rm_col(i)&
               & + f3c_nm(ieq) * rn_col(j) * rn_col(i)))
            end do              !ieq = 1,mnmx
         end do                 !j=1,mn_col
      end do                  !i=1,mn_col
      !
      !     Call matrix eigenvalue solver
      !

      egl = 1.d-2; egu = 0.6d0; abstol = 1.d-8
      il = 0; iu = 0

      if (.not. ipos_def_sym) then
         call dggev('N', 'V', mn_col, amat, mn_col, bmat, mn_col, alfr, alfi, &
         &beta, vl, ldvl, vr, ldvr, work, lwork, info)
      else if (ipos_def_sym .and. .not. subset_eq) then
         call dsygv(iopt, jobz, 'L', mn_col, amat, mn_col, bmat, mn_col, &
         &omega, work, lwork, info)
         if (info .ne. 0) write (*, '("info = ",i8)') info

      else if (ipos_def_sym .and. subset_eq) then
         call dsygvx(iopt, 'V', 'V', 'L', mn_col, amat, mn_col, bmat, mn_col, &
         &egl, egu, il, iu, abstol, m_red, omega, vr, mn_col, &
         &work, lwork, iwork, ifail, info)
         if (info .ne. 0) write (*, '("info = ",i8)') info

      end if
      !
      do i = 1, mn_col
         if (iopt .eq. 1) then
            do j = 1, mn_col
               if (ipos_def_sym .and. jobz .eq. 'V' .and. .not. subset_eq)&
               &eig_vect(j) = abs(amat(j, i))
               if (.not. ipos_def_sym) eig_vect(j) = abs(vr(j, i))
               !        if(ir .eq. 5 .and. i .eq. mn_col/2)
               !     >   write(*,'(e15.7,2x,i4,2x,i4)') eig_vect(j),
               !     >   im_col(j), in_col(j)
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
               write (21, '(2(e15.7,2x),i4,2x,i4)') rho_fine(ir), sqrt(abs(omega(i))), &
               &m_emax, n_emax
            else if (.NOT. ipos_def_sym) then
               write (21, '(4(e15.7,2x),i4,2x,i4)') rho_fine(ir), alfr(i), &
               &alfi(i), beta(i), m_emax, n_emax
            end if
         else if (iopt .eq. 0) then
            if (ipos_def_sym) then
               write (21, '(e15.7,2x,e15.7)') rho_fine(ir), sqrt(abs(omega(i)))
            else if (.NOT. ipos_def_sym) then
               write (21, '(e15.7,3(2x,e15.7))') rho_fine(ir), alfr(i), &
               &alfi(i), beta(i)
            end if
         end if
      end do          !do i=1,mn_col

   end do       !ir=1,ir_fine_scl
   !
   !
   close (unit = 21)
   close (unit = 8)

   !  Write and deallocate
   write (*, '("modes = ",i5,2x,"no. of radial points = ",i5)') mn_col, ir_fine_scl

   call trg_deallocate
   ! deallocate (alpha, beta, aux)
   deallocate (alpha, beta)

   call write_data_post(iopt, mn_col, ir_fine_scl, isym_opt)

   call post_process

   stop
end program tae_continua
