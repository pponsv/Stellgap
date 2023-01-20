!
! Variables that need to be set for each time the configuration is changed:
!   ith, izt, nfp, mpol, ntor, mp_col, nt_col, mode_family, tokamak
!
module fourier_lib
	
	use kind_spec
	use input, only : ith, izt, mpol, ntor, nznt, mnmx, ntors, nfp, &
	 mode_family, nw, mwl, mwu
	
	implicit none
	
	! integer :: ith, izt, mpol, ntor, nznt, mnmx, ntors, nfp, mode_family
	real(r8), parameter :: parity_gss = 1., parity_bsupth = 1., parity_bsupzt = 1.
	integer :: i, j, lg, nu, nl, m, n, mn, istat, sin_type, cos_type
	!   Equilibrium coefficient arrays
	real(r8), allocatable :: thtgrd(:), ztgrd(:), rn(:), rm(:)
	real(r8), allocatable :: fnm(:), f(:), anm(:)
	real(r8), allocatable :: cos_ar(:, :), sin_ar(:, :)
	real(r8), allocatable :: cos_toF(:, :), sin_toF(:, :)
	
	!    The subset_eq flag allows one to just request a bracketed subset of
	!     eigenvalues, rather than all eigenvalues. This option has not been
	!     developed beyond the call. Initial tests have not indicated it
	!     speeds things up.
	logical, parameter :: subset_eq = .false.
	!   Eigenfunction arrays,variables
	integer :: mn_col, ith_col, izt_col
	real(r8), allocatable :: rn_col(:), rm_col(:)
	real(r8), allocatable :: rn2_col(:), rm2_col(:), rnm_col(:)
	integer, allocatable :: in_col(:), im_col(:)
	integer, allocatable :: in_col_aug(:), im_col_aug(:)
	! integer, allocatable :: nw(:), mwl(:), mwu(:)

contains
	!
	subroutine readin
		integer :: i, istat
		open (unit = 20, file = "fourier.dat", status = "old")
		read (20, *) nfp, ith, izt, mode_family
		mpol = ith * 2 / 5; ntor = izt * 2 / 5 ! Done to avoid aliasing issues
		nznt = ith * izt; mnmx = (2 * ntor + 1) * mpol - ntor
		read (20, *) ntors
		allocate (nw(ntors), stat = istat)
		allocate (mwl(ntors), stat = istat)
		allocate (mwu(ntors), stat = istat)
		do i = 1, ntors
			read (20, *) nw(i), mwl(i), mwu(i)
		end do
		close (unit = 20)
	end subroutine readin
	!
	!
	subroutine trig_array
		integer :: istat
		real(r8) :: dum, dnorm, arg
		allocate (thtgrd(nznt), stat = istat)
		allocate (ztgrd(nznt), stat = istat)
		allocate (rm(mnmx), stat = istat)
		allocate (rn(mnmx), stat = istat)
		allocate (fnm(mnmx), stat = istat)
		allocate (f(nznt), stat = istat)
		allocate (anm(mnmx), stat = istat)
		allocate (cos_ar(nznt, mnmx), stat = istat)
		allocate (sin_ar(nznt, mnmx), stat = istat)
		allocate (cos_toF(nznt, mnmx), stat = istat)
		allocate (sin_toF(nznt, mnmx), stat = istat)
		!    Generate theta, zeta grid
		lg = 0
		do i = 1, izt
			do j = 1, ith
				lg = lg + 1
				ztgrd(lg) = twopi * real(i - 1) / (real(nfp * izt))
				thtgrd(lg) = twopi * real(j - 1) / real(ith)
			end do
		end do
		!    Generate Fourier mode distribution
		mn = 0
		nu = ntor
		do m = 0, mpol - 1
			nl = -ntor
			if (m .eq. 0) nl = 0
			do n = nl, nu
				mn = mn + 1
				rm(mn) = real(m)
				rn(mn) = real(n * nfp)
			end do
		end do
		do i = 1, nznt ! =ith*izt
			do mn = 1, mnmx ! =(2*ntor + 1)*mpol - ntor
				arg = -rn(mn) * ztgrd(i) + rm(mn) * thtgrd(i)
				cos_ar(i, mn) = cos(arg)
				sin_ar(i, mn) = sin(arg)
				dnorm = 2. / real(nznt)
				dum = abs(rn(mn)) + abs(rm(mn))
				if (nint(dum) .eq. 0) dnorm = .5 * dnorm
				cos_toF(i, mn) = cos_ar(i, mn) * dnorm
				sin_toF(i, mn) = sin_ar(i, mn) * dnorm
			end do
		end do
	end subroutine trig_array
	
	subroutine convolution_array
		integer :: istat, count
		!     First, count the number of modes to be used
		count = 0
		do i = 1, ntors
			count = count + mwu(i) - mwl(i) + 1
		end do
		mn_col = count
		!     Allocate memory for arrays needed in eigenfunction calculations
		allocate (rm_col(mn_col), stat = istat)
		allocate (rn_col(mn_col), stat = istat)
		allocate (im_col(mn_col), stat = istat)
		allocate (in_col(mn_col), stat = istat)
		allocate (rm2_col(mn_col), stat = istat)
		allocate (rn2_col(mn_col), stat = istat)
		allocate (rnm_col(mn_col), stat = istat)
		!     Create mode number arrays needed for eigenfunction
		count = 0
		do n = 1, ntors
			do m = mwl(n), mwu(n)
				count = count + 1
				rm_col(count) = real(m)
				rn_col(count) = real(nw(n))
				im_col(count) = m
				in_col(count) = nw(n)
				rm2_col(count) = rm_col(count) * rm_col(count)
				rn2_col(count) = rn_col(count) * rn_col(count)
				rnm_col(count) = rm_col(count) * rn_col(count)
			end do
		end do
		mn_col = count
		!      write(*,*) mn_col, count
	end subroutine convolution_array
	!
	!
	subroutine scs_convolve(ans, m1, n1, m2, n2, meq, neq)
		integer :: istat
		real(r8) :: tht_int1, tht_int2, tht_int3, tht_int4, &
		 &zeta_int1, zeta_int2, zeta_int3, zeta_int4, ans
		integer :: m1, m2, n1, n2, meq, neq
		integer :: sm1, sm2, sn1, sn2, smeq, sneq
		!
		!     This subroutine calculates the 2D integral for tht and zeta running
		!     from 0 to 2*PI of:
		!     sin(m1*tht - n1*zeta)*cos(meq*tht - neq*zeta)*sin(m2*tht - n2*zeta)
		!     The value returned is actually (2/PI) times this integral.
		!
		sm1 = 1; sm2 = 1; sn1 = 1; sn2 = 1; smeq = 1; sneq = 1
		if (m1 .ne. 0) sm1 = m1 / abs(m1)
		if (m2 .ne. 0) sm2 = m2 / abs(m2)
		if (n1 .ne. 0) sn1 = n1 / abs(n1)
		if (n2 .ne. 0) sn2 = n2 / abs(n2)
		if (meq .ne. 0) smeq = meq / abs(meq)
		if (neq .ne. 0) sneq = neq / abs(neq)
		m1 = sm1 * m1; n1 = sn1 * n1
		m2 = sm2 * m2; n2 = sn2 * n2
		meq = smeq * meq; neq = sneq * neq
		call css(tht_int1, meq, m1, m2)
		call ccc(zeta_int1, n1, n2, neq)
		call ccc(tht_int2, m1, m2, meq)
		call css(zeta_int2, neq, n1, n2)
		call css(tht_int3, m2, m1, meq)
		call css(zeta_int3, n1, n2, neq)
		call css(tht_int4, m1, m2, meq)
		call css(zeta_int4, n2, n1, neq)
		ans = tht_int1 * zeta_int1 * sm1 * sm2&
		 & + tht_int2 * zeta_int2 * sn1 * sn2&
		 & - tht_int3 * zeta_int3 * sm1 * smeq * sn2 * sneq&
		 & - tht_int4 * zeta_int4 * sm2 * smeq * sn1 * sneq
		m1 = sm1 * m1; n1 = sn1 * n1
		m2 = sm2 * m2; n2 = sn2 * n2
		meq = smeq * meq; neq = sneq * neq
	end subroutine scs_convolve
	!
	!
	!
	subroutine ccc_convolve(ans, m1, n1, m2, n2, meq, neq)
		integer :: istat
		real(r8) :: tht_int1, tht_int2, tht_int3, tht_int4, &
		 &zeta_int1, zeta_int2, zeta_int3, zeta_int4, ans
		integer :: m1, m2, n1, n2, meq, neq
		integer :: sm1, sm2, sn1, sn2, smeq, sneq
		!
		!     This subroutine calculates the 2D integral for tht and zeta running
		!     from 0 to 2*PI of:
		!     cos(m1*tht - n1*zeta)*cos(meq*tht - neq*zeta)*cos(m2*tht - n2*zeta)
		!     The value returned is actually (2/PI) times this integral.
		!
		sm1 = 1; sm2 = 1; sn1 = 1; sn2 = 1; smeq = 1; sneq = 1
		if (m1 .ne. 0) sm1 = m1 / abs(m1)
		if (m2 .ne. 0) sm2 = m2 / abs(m2)
		if (n1 .ne. 0) sn1 = n1 / abs(n1)
		if (n2 .ne. 0) sn2 = n2 / abs(n2)
		if (meq .ne. 0) smeq = meq / abs(meq)
		if (neq .ne. 0) sneq = neq / abs(neq)
		m1 = sm1 * m1; n1 = sn1 * n1
		m2 = sm2 * m2; n2 = sn2 * n2
		meq = smeq * meq; neq = sneq * neq
		call ccc(tht_int1, m1, m2, meq)
		call ccc(zeta_int1, n1, n2, neq)
		call css(tht_int2, meq, m1, m2)
		call css(zeta_int2, neq, n1, n2)
		call css(tht_int3, m1, m2, meq)
		call css(zeta_int3, n1, n2, neq)
		call css(tht_int4, m2, m1, meq)
		call css(zeta_int4, n2, n1, neq)
		ans = tht_int1 * zeta_int1&
		 & + tht_int2 * zeta_int2 * sm1 * sm2 * sn1 * sn2&
		 & + tht_int3 * zeta_int3 * sm2 * smeq * sn2 * sneq&
		 & + tht_int4 * zeta_int4 * sm1 * smeq * sn1 * sneq
		m1 = sm1 * m1; n1 = sn1 * n1
		m2 = sm2 * m2; n2 = sn2 * n2
		meq = smeq * meq; neq = sneq * neq
	end subroutine ccc_convolve
	!
	!
	!
	subroutine ccc(result, i, j, k)
		integer :: istat
		real(r8), parameter :: zero = 0, one = 1, &
		 &two = 2, four = 4
		real(r8) :: result
		integer :: i, j, k, izeros
		!
		!     This subroutine calculates the 1D integral of
		!     cos(i*x)*cos(j*x)*cos(k*x)
		!     for x running from 0 to 2*PI.
		!     The value returned is actually (2/PI) times this integral.
		!
		result = zero
		izeros = 0
		if (i .eq. 0) izeros = izeros + 1
		if (j .eq. 0) izeros = izeros + 1
		if (k .eq. 0) izeros = izeros + 1
		if (izeros .ne. 0) then
			if (izeros .eq. 3) then
				result = four
				return
			else if (izeros .eq. 2) then
				result = zero
				return
			else if (izeros .eq. 1) then
				if (i .eq. 0) then
					result = zero
					if (j * j .eq. k * k) then
						result = two
					end if
					return
				else if (j .eq. 0) then
					result = zero
					if (i * i .eq. k * k) then
						result = two
					end if
					return
				else if (k .eq. 0) then
					result = zero
					if (i * i .eq. j * j) then
						result = two
					end if
					return
				end if
			end if
		else if (izeros .eq. 0) then
			if (k .eq. (i - j)) result = one
			if (k .eq. (i + j)) result = one
			if (k .eq. (j - i)) result = one
			if (k .eq. -(i + j)) result = one
		end if
	end subroutine ccc
	!
	!
	!
	subroutine css(result, k, i, j)
		integer :: istat
		real(r8), parameter :: one = 1, neg_one = -1, zero = 0, &
		 &two = 2, neg_two = -2
		integer :: i, j, k
		real(r8) :: result
		!
		!     This subroutine calculates the 1D integral of
		!     cos(k*x)*sin(i*x)*sin(j*x)
		!     for x running from 0 to 2*PI.
		!     The value returned is actually (2/PI) times this integral.
		!
		result = zero
		if (i .eq. 0 .or. j .eq. 0) return
		if (k .eq. 0) then
			if (i * i .ne. j * j) then
				return
			else if (i * i .eq. j * j) then
				result = 2
				return
			end if
		end if
		if (k .eq. (i - j)) result = one
		if (k .eq. (i + j)) result = neg_one
		if (k .eq. (j - i)) result = one
		if (k .eq. -(i + j)) result = neg_one
	end subroutine css
	!
	!
	subroutine trg_deallocate
		deallocate (thtgrd, ztgrd, rm, rn, cos_ar, sin_ar, f, fnm, anm, &
		 &rm_col, rn_col, im_col, in_col, &
		 &rm2_col, rn2_col, rnm_col)
	end subroutine trg_deallocate

end module fourier_lib
