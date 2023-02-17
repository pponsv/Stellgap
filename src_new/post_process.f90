module postprocess

   use kind_spec

   implicit none

contains

   subroutine post_process
      integer :: nang2, irads, iopt, m_emax, n_emax, itot, isym_pos, istat, ir, i
      real(r8) :: r_pt, lam_min, lambda_imag, lam_max, cond_no, lambda_real, fnorm
      real(r8), allocatable, dimension(:) :: temp
      real(r8), allocatable, dimension(:,:) :: omega, alpha_r, alpha_i, beta, &
         omega_r, omega_i
      complex*16, allocatable, dimension(:,:) :: omega2

      open (unit=4, file="alfven_spec", status="old")
      open (unit=7, file="alfven_post", status="unknown")
      open (unit=9, file="cond_no", status="unknown")
      open (unit=8, file="data_post", status="old")
      read (8, *) iopt, nang2, irads, isym_pos
!           if positive-definite matrix option is used in stellgap
!            (subroutine DSYGV) then isym_pos = 1
!           if general matrix option is used in stellgap
!            (subroutine DGEGV) then isym_pos = 0
!            isym_pos = 0
      write (*, *) iopt, nang2, irads, isym_pos
      allocate (omega(irads, nang2), stat=istat)
      allocate (omega2(irads, nang2), stat=istat)
      allocate (omega_r(irads, nang2), stat=istat)
      allocate (omega_i(irads, nang2), stat=istat)
      allocate (beta(irads, nang2), stat=istat)
      allocate (alpha_r(irads, nang2), stat=istat)
      allocate (alpha_i(irads, nang2), stat=istat)
      allocate (temp(nang2), stat=istat)
      itot = irads*nang2
!            fnorm = 1./(2000.*pi)    !converts omega from radians/sec to kHz
      fnorm = 1.      !conversion to kHz is now done in stellgap
!            fnorm = 1./224.
      write (7, '(i10)') itot
      lam_min = 1.d+10
      lam_max = -1.d+10
      do ir = 1, irads
         do i = 1, nang2
            if (isym_pos .eq. 1) then
               read (4, '(2(e15.7,2x),i4,2x,i4)') r_pt, omega(ir, i),&
               &m_emax, n_emax
               temp(i) = omega(ir, i)**2
               write (7, '(2(e15.7,2x),i4,2x,i4)') r_pt,&
               &omega(ir, i)*fnorm, m_emax, n_emax

            else if (isym_pos .eq. 0) then
               read (4, '(4(e15.7,2x),i4,2x,i4)') r_pt,&
               &alpha_r(ir, i),&
               &alpha_i(ir, i), beta(ir, i), m_emax, n_emax
               !   c               alpha_r(ir,i) = abs(alpha_r(ir,i))
               if (beta(ir, i) .gt. 0.d0) then
                  lambda_real = alpha_r(ir, i)/beta(ir, i)
                  if (lambda_real .gt. 0.) then
                     lam_min = min(lam_min, abs(lambda_real))
                     lam_max = max(lam_max, abs(lambda_real))
                  end if
                  lambda_imag = alpha_i(ir, i)/beta(ir, i)
                  omega2(ir, i) = dcmplx(lambda_real, lambda_imag)
               else if (beta(ir, i) .le. 0.d0) then
                  omega2(ir, i) = 1.e+10
               end if
               omega_r(ir, i) = real(cdsqrt(omega2(ir, i)))
               omega_i(ir, i) = imag(cdsqrt(omega2(ir, i)))
               if (omega_i(ir, i) .eq. 0.d0) then
                  if (r_pt .lt. 0.99) then
                     write (7, '(2(e15.7,2x),i4,2x,i4)') r_pt,&
                     &omega_r(ir, i)*fnorm,&
                     &m_emax, n_emax
                  end if
               end if

            end if

         end do
         if (isym_pos .eq. 1) then
            lam_min = minval(temp)
            lam_max = maxval(temp)
         end if
         if (lam_min .ne. 0.) cond_no = lam_max/lam_min
         write (9, '(e15.7,2x,e15.7)') r_pt, cond_no
         !   c          write(*,*) r_pt, cond_no
      end do

   end subroutine post_process

   subroutine post_process_snd

      real(r8), allocatable :: omega(:, :)
      real(r8), allocatable :: alpha_r(:, :), alpha_i(:, :),&
      &beta(:, :), omega_r(:, :), omega_i(:, :)
      complex*16, allocatable :: omega2(:, :)
      real(r8), allocatable :: temp(:)
      real(r8) :: r_pt, lam_min,&
      &lam_max, cond_no, lambda_real, xnrm,&
      &phi_norm, psi_norm, fnorm, lambda_imag
      integer :: nang2, irads, iopt, m_emax, n_emax,&
      &itot, m_emax_psi, n_emax_psi, isym_pos, istat,&
      &ir, i
      call system("rm all_alfven_spec")
      call system("cat alfven_spec* > all_alfven_spec")
      open (unit=4, file="all_alfven_spec", status="old")
      open (unit=7, file="alfven_post", status="unknown")
      open (unit=9, file="cond_no", status="unknown")
      open (unit=8, file="data_post", status="old")
      read (8, *) iopt, nang2, irads, isym_pos
!        if positive-definite matrix option is used in stellgap
!         (subroutine DSYGV) then isym_pos = 1
!        if general matrix option is used in stellgap
!         (subroutine DGEGV) then isym_pos = 0
!         isym_pos = 0
      write (*, *) iopt, nang2, irads, isym_pos
      allocate (omega(irads, nang2), stat=istat)
      allocate (omega2(irads, nang2), stat=istat)
      allocate (omega_r(irads, nang2), stat=istat)
      allocate (omega_i(irads, nang2), stat=istat)
      allocate (beta(irads, nang2), stat=istat)
      allocate (alpha_r(irads, nang2), stat=istat)
      allocate (alpha_i(irads, nang2), stat=istat)
      allocate (temp(nang2), stat=istat)
      itot = irads*nang2
!         fnorm = 1./(2000.*pi)    !converts omega from radians/sec to kHz
      fnorm = 1.      !conversion to kHz is now done in stellgap
      write (7, '(i10)') itot
      lam_min = 1.d+10
      lam_max = -1.d+10
      do ir = 1, irads
         do i = 1, nang2
            if (isym_pos .eq. 1) then
               read (4, '(2(e15.7,2x),i4,2x,i4)') r_pt, omega(ir, i),&
               &m_emax, n_emax
               temp(i) = omega(ir, i)**2
               write (7, '(2(e15.7,2x),i4,2x,i4)') r_pt,&
               &omega(ir, i)*fnorm, m_emax, n_emax

            else if (isym_pos .eq. 0) then
               read (4, '(4(e15.7,2x),i4,3(2x,i4),2(2x,e15.7))')&
               &r_pt, alpha_r(ir, i),&
               &alpha_i(ir, i), beta(ir, i), m_emax, n_emax,&
               &m_emax_psi, n_emax_psi, phi_norm, psi_norm
               xnrm = sqrt(phi_norm**2 + psi_norm**2)
!               alpha_r(ir,i) = abs(alpha_r(ir,i))
               if (beta(ir, i) .gt. 0.d0) then
                  lambda_real = alpha_r(ir, i)/beta(ir, i)
                  if (lambda_real .gt. 0.) then
                     lam_min = min(lam_min, abs(lambda_real))
                     lam_max = max(lam_max, abs(lambda_real))
                  end if
                  lambda_imag = alpha_i(ir, i)/beta(ir, i)
                  omega2(ir, i) = dcmplx(lambda_real, lambda_imag)
               else if (beta(ir, i) .le. 0.d0) then
                  omega2(ir, i) = 1.e+10
               end if
               omega_r(ir, i) = real(cdsqrt(omega2(ir, i)))
               omega_i(ir, i) = imag(cdsqrt(omega2(ir, i)))
               if (omega_i(ir, i) .eq. 0.d0) then
                  if (r_pt .lt. 0.99) then
                     write (7, '(2(e15.7,2x),i4,2x,i4,2(2x,e15.7))')&
                     &r_pt, omega_r(ir, i)*fnorm,&
                     &m_emax, n_emax, phi_norm/xnrm, psi_norm/xnrm
                  end if
               end if

            end if

         end do
         if (isym_pos .eq. 1) then
            lam_min = minval(temp)
            lam_max = maxval(temp)
         end if
         if (lam_min .ne. 0.) cond_no = lam_max/lam_min
         write (9, '(e15.7,2x,e15.7)') r_pt, cond_no
!            write(*,*) r_pt, cond_no
      end do

   end subroutine post_process_snd

end module postprocess
