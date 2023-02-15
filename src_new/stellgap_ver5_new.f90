
program tae_continua

   use fourier_lib
   use fourier_lib_convolve
   use kind_spec
   use postprocess
   use fitpack
   use input
   use output
   use globals
   use helper

   implicit none
   !
   !      This code is a test of a Fourier xform of a 2D function
   !      (f) vs. theta, zeta.  f has cos/sin stellarator-like symmetry.
   !      The Fourier xform is done by subroutine toFourier using
   !      cos/sin expansions. (use with Fourier_lib.f)
   !

   external dggev, dsygv, dsygvx ! LAPACK subroutines

   integer, parameter :: iopt = 1


   integer :: ir, irr, il, iu
   real(r8), dimension(:, :, :), allocatable :: bsupth, bsupzt, bfield_lrg, &
   &gsssup_lrg, bsupth_lrg, bsupzt_lrg, rjacob_lrg
   real(r8), dimension(:), allocatable :: f1_nm, f3a_nm, f3b_nm, f3c_nm, bsupth_tmp, bsupzt_tmp
   real(r8), dimension(:, :), allocatable :: f1, f3a, &
   &f3b, f3c



   real(r8), allocatable :: beta(:), aux(:), eig_vect(:)
   real(r8), allocatable :: amat(:, :), &
   &bmat(:, :), z_r(:, :), vr(:, :), vl(:, :)
   real(r8), allocatable :: omega(:), work(:), alfr(:), &
   &alfi(:)
   complex*16, allocatable :: alpha(:)

   real(r8), dimension(:), allocatable :: sp_fit
   ! real(r8), dimension(:), allocatable :: nsurf, iotac,&
   ! &phipc, rho, sp_fit, b2, iotac_inv
   real(r8) :: sp1, sp2, sigma_spl, r_pt

   real(r8), dimension(:), allocatable :: yp, temp, &
   &ypi, tempi

   real(r8) :: va, eig_max, ccci, scsi, &
      egl, egu, abstol, f1_avg, f3_avg
   integer :: ispl_opt, ierr_spl, naux
   integer :: npes, numrads
   integer :: ni, nj, mi, mj, ieq, meq, neq
   integer :: j_max_index
   integer :: m_emax, n_emax, isym_opt
   integer :: ldvl, ldvr, lwork, info, m_red
   integer, allocatable :: ifail(:), iwork(:)
   logical :: cyl
   character*20 outfile
   character*3 procnum
   character*1 jobz
   character*10 date, time, zone
   integer values(8)

   cyl = .false.
   lrfp = .false.

   call read_args

   call read_plasma_dat

   call read_fourier_dat

   allocate (bavg(ir_fine_scl), mu0_rho_ion(ir_fine_scl), &
      ion_density(ir_fine_scl),iota_r(ir_fine_scl), &
      iota_r_inv(ir_fine_scl), stat = istat)
   allocate (iotac(irads), phipc(irads), sp_fit(irads), stat = istat)
   allocate (yp(3*irads), temp(3*irads), ypi(3*irads), &
      tempi(3*irads), stat = istat)

   mype = 0
   npes = 1
   jobz = 'V'
   outfile = "alfven_spec"

   if (mype .eq. 0) then
      write (*, *) mype, npes
      write (*, *) procnum
      write (*, *) outfile
      call date_and_time(date, time, zone, values)
      write (*, *) time
   end if
   if (ipos_def_sym) isym_opt = 1
   if (.NOT. ipos_def_sym) isym_opt = 0

   !
   !    Generate Fourier arrays
   !
   ! call readin

   call trig_array

   call convolution_array

   naux = 10 * mn_col
   ! mu0 = 2.d-7*TWOPI
   scale_khz = (1.d+3 * TWOPI)**2
   ldvl = mn_col; ldvr = mn_col; lwork = 20 * mn_col
   allocate (aux(naux), stat = istat)
   allocate (alpha(mn_col), stat = istat)
   allocate (beta(mn_col), stat = istat)
   allocate (eig_vect(mn_col), stat = istat)
   allocate (z_r(mn_col, mn_col), stat = istat)
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
   allocate (bfield(izt, ith, irads), stat = istat)
   allocate (gsssup(izt, ith, irads), stat = istat)
   allocate (rjacob(izt, ith, irads), stat = istat)
   allocate (bsupth(izt, ith, irads), stat = istat)
   allocate (bsupzt(izt, ith, irads), stat = istat)
   allocate (bfield_lrg(izt, ith, ir_fine_scl), stat = istat)
   allocate (gsssup_lrg(izt, ith, ir_fine_scl), stat = istat)
   allocate (bsupth_lrg(izt, ith, ir_fine_scl), stat = istat)
   allocate (bsupzt_lrg(izt, ith, ir_fine_scl), stat = istat)
   allocate (rjacob_lrg(izt, ith, ir_fine_scl), stat = istat)
   allocate (theta_tae(ith), stat = istat)
   allocate (zeta_tae(izt), stat = istat)
   allocate (f1_nm(mnmx), stat = istat)
   allocate (f3a_nm(mnmx), stat = istat)
   allocate (f3b_nm(mnmx), stat = istat)
   allocate (f3c_nm(mnmx), stat = istat)
   allocate (bsupth_tmp(nznt), stat = istat)
   allocate (bsupzt_tmp(nznt), stat = istat)
   allocate (f1(izt, ith), stat = istat)
   allocate (f3a(izt, ith), stat = istat)
   allocate (f3b(izt, ith), stat = istat)
   allocate (f3c(izt, ith), stat = istat)
   !
   !   Open files for output
   !
   !      write(*,*) trim(adjustl(outfile))

   if (mype .eq. 0) then
      write (*, *) izt, ith, irads, mnmx, nznt, mn_col
   end if

   open (unit = 21, file = trim(adjustl(outfile)), status = "unknown")
   open (unit = 8, file = "coef_arrays", status = "unknown")
   if (mype .eq. 0) then
      open (unit = 7, file = "data_post", status = "unknown")
      ! open (unit = 9, file = "ion_profile", status = "unknown")
   end if
   !
   !    Boozer coordinates input - new ae-mode-structure input
   !
   call read_tae_data_boozer

   !
   !    Record equilibrium and eigenfunction Fourier mode list
   !
   if (mype .eq. 0) then
      call write_modes(mnmx = mnmx, mn_col = mn_col, rm = rm, rn = rn, im_col = im_col, in_col = in_col)
   end if
   !
   !   Make spline fits and fill in fine_radius_scale arrays
   !
   ! do ir = 1, irads
   !    rho(ir) = nsurf(ir) / nsurf(irads)
   ! end do
   ispl_opt = 3; ierr_spl = 0; sigma_spl = 0.
   call curv1(irads, rho, iotac, sp1, sp2, ispl_opt, ypi, tempi, sigma_spl, ierr_spl)

      ! Interpola rho a ir_fine_scl
   rho_fine = rho(1) + (rho(irads) - rho(1)) * [(irr, irr=0,ir_fine_scl-1)] / (ir_fine_scl - 1)
   do irr = 1, ir_fine_scl
      iota_r(irr) = curv2(rho_fine(irr), irads, rho, iotac, ypi, sigma_spl)
   end do

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

   if (mype .eq. 0) call write_ion_profile

   ! do irr=1, ir_fine_scl
   !    va = sqrt(bfavg**2 / (mu0_rho_ion(irr) / scale_khz))
   !    if (mype .eq. 0) then
   !       write (9, fmt_ion_profile) rho_fine(irr), ion_density_0 * ion_density(irr), iota_r(irr), va
   !    end if
   ! end do

   !   Interpolate 1/iota for lrfp = true option
   ispl_opt = 3; ierr_spl = 0; sigma_spl = 0.
   call curv1(irads, rho, 1._r8/iotac, sp1, sp2, ispl_opt, ypi, tempi, sigma_spl, ierr_spl)
   do irr = 1, ir_fine_scl
      ! r_pt = rho(1) + real(irr - 1) * (rho(irads) - rho(1))&
      ! & / real(ir_fine_scl - 1)
      iota_r_inv(irr) = curv2(rho_fine(irr), irads, rho, 1._r8/iotac, ypi, sigma_spl)
   end do

   if (ierr_spl .ne. 0) then
      write (*, '("spline error 1",i3)') ierr_spl
   end if

! 67 format(e12.5, 3(3x, e12.5))
   do i = 1, izt
      do j = 1, ith
         !
         !       Bfield array
         !
         do ir = 1, irads
            sp_fit(ir) = bfield(i, j, ir)
         end do      !ir=1,irads
         ispl_opt = 3; ierr_spl = 0; sigma_spl = 0
         call curv1(irads, rho, sp_fit, sp1, sp2, ispl_opt, yp, &
         &temp, sigma_spl, ierr_spl)
         if (ierr_spl .ne. 0) write (*, '("spline error 2",i3)') ierr_spl
         do irr = 1, ir_fine_scl
            r_pt = rho(1) + real(irr - 1) * (rho(irads) - rho(1))&
            & / real(ir_fine_scl - 1)
            bfield_lrg(i, j, irr) = curv2(r_pt, irads, rho, sp_fit, yp, sigma_spl)
         end do     !irr = 1,ir_fine_scl
         !
         !
         !       Jacobian array
         !
         do ir = 1, irads
            sp_fit(ir) = rjacob(i, j, ir)
         end do      !ir=1,irads
         ispl_opt = 3; ierr_spl = 0; sigma_spl = 0
         call curv1(irads, rho, sp_fit, sp1, sp2, ispl_opt, yp, &
         &temp, sigma_spl, ierr_spl)
         if (ierr_spl .ne. 0) write (*, '("spline error 2",i3)') ierr_spl
         do irr = 1, ir_fine_scl
            r_pt = rho(1) + real(irr - 1) * (rho(irads) - rho(1))&
            & / real(ir_fine_scl - 1)
            rjacob_lrg(i, j, irr) = curv2(r_pt, irads, rho, sp_fit, yp, sigma_spl)
         end do     !irr = 1,ir_fine_scl
         !
         !       Gsssup array
         !
         do ir = 1, irads
            sp_fit(ir) = gsssup(i, j, ir)
         end do
         ispl_opt = 3; ierr_spl = 0; sigma_spl = 0
         call curv1(irads, rho, sp_fit, sp1, sp2, ispl_opt, yp, &
         &temp, sigma_spl, ierr_spl)
         if (ierr_spl .ne. 0) write (*, '("spline error 3",i3)') ierr_spl
         do irr = 1, ir_fine_scl
            r_pt = rho(1) + real(irr - 1) * (rho(irads) - rho(1))&
            & / real(ir_fine_scl - 1)
            gsssup_lrg(i, j, irr) = curv2(r_pt, irads, rho, sp_fit, yp, sigma_spl)
         end do     !irr = 1,ir_fine_scl
      end do     !j=1,ith
   end do     !i=1,izt
   ! close (unit = 9)
   !
   !
   !     Fine_scale arrays and spline fits finished
   !
   numrads = ir_fine_scl / npes
   !       write(*,*) numrads
   do ir = (mype * numrads + 1), (mype * numrads + numrads)
      !       write(*,*) ir
      !      do ir=1,ir_fine_scl
      r_pt = rho(1) + real(ir - 1) * (rho(irads) - rho(1))&
      & / real(ir_fine_scl - 1)
      f1_avg = 0.; f3_avg = 0.
      do i = 1, izt
         do j = 1, ith
            f1(i, j) = gsssup_lrg(i, j, ir) * rjacob_lrg(i, j, ir) / &
            &(bfield_lrg(i, j, ir)**2)
            f1_avg = f1_avg + f1(i, j) / real(izt * ith)
            if (.not. lrfp) then
               f3c(i, j) = gsssup_lrg(i, j, ir) / &
               &(rjacob_lrg(i, j, ir) * (bfield_lrg(i, j, ir)**2))
               f3_avg = f3_avg + f3c(i, j) / real(izt * ith)
               f3b(i, j) = iota_r(ir) * f3c(i, j)
               f3a(i, j) = iota_r(ir) * f3b(i, j)
            else if (lrfp) then
               f3a(i, j) = gsssup_lrg(i, j, ir) / &
               &(rjacob_lrg(i, j, ir) * (bfield_lrg(i, j, ir)**2))
               f3b(i, j) = f3a(i, j) * iota_r_inv(ir)
               f3c(i, j) = f3b(i, j) * iota_r_inv(ir)
               f3_avg = f3_avg + f3a(i, j) / real(izt * ith)
            end if
         end do
      end do

      if (cyl) then
         if (.not. lrfp) then
            do i = 1, izt
               do j = 1, ith
                  f1(i, j) = f1_avg
                  !         f3a(i,j) = f3_avg*iota_r(ir)*iota_r(ir)
                  !         f3b(i,j) = f3_avg*iota_r(ir)
                  !         f3c(i,j) = f3_avg
               end do
            end do
         else if (lrfp) then
            do i = 1, izt
               do j = 1, ith
                  f1(i, j) = f1_avg
                  !         f3a(i,j) = f3_avg
                  !         f3b(i,j) = f3_avg*iota_r_inv(ir)
                  !         f3c(i,j) = f3_avg*iota_r_inv(ir)**iota_r_inv(ir)
               end do
            end do
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

      if (.NOT. ipos_def_sym) then
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
               if (.NOT. ipos_def_sym) eig_vect(j) = abs(vr(j, i))
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
               write (21, '(2(e15.7,2x),i4,2x,i4)') r_pt, sqrt(abs(omega(i))), &
               &m_emax, n_emax
            else if (.NOT. ipos_def_sym) then
               write (21, '(4(e15.7,2x),i4,2x,i4)') r_pt, alfr(i), &
               &alfi(i), beta(i), m_emax, n_emax
            end if
         else if (iopt .eq. 0) then
            if (ipos_def_sym) then
               write (21, '(e15.7,2x,e15.7)') r_pt, sqrt(abs(omega(i)))
            else if (.NOT. ipos_def_sym) then
               write (21, '(e15.7,3(2x,e15.7))') r_pt, alfr(i), &
               &alfi(i), beta(i)
            end if
         end if
      end do          !do i=1,mn_col

   end do       !ir=1,ir_fine_scl
   !
   !
   close (unit = 20)
   close (unit = 21)
   close (unit = 8)
   if (mype .eq. 0)  then
      write (*, '("modes = ",i5,2x,"no. of radial points = ",i5)') &
         mn_col, ir_fine_scl
      write (7, '(i5,3(2x,i5))') iopt, mn_col, ir_fine_scl, isym_opt
      close (unit = 7)
      call trg_deallocate
      deallocate (alpha, beta, aux)
      call date_and_time(date, time, zone, values)
      write (*, *) time
   end if

   call post_process

   stop
end program tae_continua
