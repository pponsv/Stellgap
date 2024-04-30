!
! Variables that need to be set for each time the configuration is changed:
!   ith, izt, nfp, mpol, ntor, mp_col, nt_col, mode_family, tokamak
!
module fourier_lib

   use constants
   ! use globals
   use helper

   implicit none


contains

   subroutine trig_array
      use globals, only: izt, ith, nfp, n_fourier, m_fourier, mnmx, mpol, ntor, fnm, &
         f, anm, cos_toF, sin_toF!, cos_ar, sin_ar

      integer :: istat, nl, i, n, m, mn
      real(r8) :: dum, dnorm, arg
      real(r8) :: zetas(izt), thetas(ith)
      real(r8) :: thtgrd(izt*ith), ztgrd(izt*ith), tmp(izt*ith)


      allocate (m_fourier(mnmx), stat = istat)
      allocate (n_fourier(mnmx), stat = istat)
      allocate (fnm(mnmx), stat = istat)
      allocate (f(ith*izt), stat = istat)
      allocate (anm(mnmx), stat = istat)
      allocate (cos_toF(ith*izt, mnmx), stat = istat)
      ! allocate (cos_ar(ith*izt, mnmx), stat = istat)
      allocate (sin_toF(ith*izt, mnmx), stat = istat)
      ! allocate (sin_ar(ith*izt, mnmx), stat = istat)

      !    Generate theta, zeta grid
      zetas  = real_linspace_nolast(0._r8, 2*PI/nfp, izt)
      thetas = real_linspace_nolast(0._r8, 2*PI, ith)
      ztgrd  = pack(spread(zetas, dim=1, ncopies=ith), .true.)
      thtgrd = pack(spread(zetas, dim=2, ncopies=ith), .true.)

      !    Generate Fourier mode distribution
      mn = 0
      do m = 0, mpol - 1
         nl = -ntor
         if (m .eq. 0) nl = 0
         do n = nl, ntor
            mn = mn + 1
            m_fourier(mn) = real(m)
            n_fourier(mn) = real(n * nfp)
         end do
      end do

      !  Make Fourier arrays
      do i = 1, ith*izt
         do mn = 1, mnmx ! =(2*ntor + 1)*mpol - ntor
            arg = -n_fourier(mn) * ztgrd(i) + m_fourier(mn) * thtgrd(i)
            ! cos_ar(i, mn) = cos(arg)
            ! sin_ar(i, mn) = sin(arg)
            dnorm = 2. / real(ith*izt)
            dum = abs(n_fourier(mn)) + abs(m_fourier(mn))
            if (nint(dum) .eq. 0) dnorm = .5 * dnorm
            ! cos_toF(i, mn) = cos_ar(i, mn) * dnorm
            ! sin_toF(i, mn) = sin_ar(i, mn) * dnorm
            cos_toF(i, mn) = cos(arg) * dnorm
            sin_toF(i, mn) = sin(arg) * dnorm
         end do
      end do
   end subroutine trig_array


   subroutine convolution_array
      use globals, only: mn_col, ntors, mwl, mwu, rm_col, rn_col, im_col, in_col, nw
      integer :: istat, count, i, m, n
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
      !     Create mode number arrays needed for eigenfunction
      count = 0
      do n = 1, ntors
         do m = mwl(n), mwu(n)
            count = count + 1
            rm_col(count) = real(m)
            rn_col(count) = real(nw(n))
            im_col(count) = m
            in_col(count) = nw(n)
         end do
      end do
      mn_col = count
   end subroutine convolution_array


   subroutine scs_convolve(ans, m1, n1, m2, n2, meq, neq)
      real(r8), intent(out) :: ans
      integer, intent(in) :: m1, m2, n1, n2, meq, neq
      real(r8) :: tht_int1, tht_int2, tht_int3, tht_int4, &
         zeta_int1, zeta_int2, zeta_int3, zeta_int4
      !
      !     This subroutine calculates the 2D integral for tht and zeta running
      !     from 0 to 2*PI of:
      !     sin(m1*tht - n1*zeta)*cos(meq*tht - neq*zeta)*sin(m2*tht - n2*zeta)
      !     The value returned is actually (2/PI) times this integral.
      !
      call css(tht_int1, abs(meq), abs(m1), abs(m2))
      call ccc(zeta_int1, abs(n1), abs(n2), abs(neq))
      call ccc(tht_int2, abs(m1), abs(m2), abs(meq))
      call css(zeta_int2, abs(neq), abs(n1), abs(n2))
      call css(tht_int3, abs(m2), abs(m1), abs(meq))
      call css(zeta_int3, abs(n1), abs(n2), abs(neq))
      call css(tht_int4, abs(m1), abs(m2), abs(meq))
      call css(zeta_int4, abs(n2), abs(n1), abs(neq))
      ans = tht_int1 * zeta_int1 * sign(1, m1) * sign(1, m2) &
         + tht_int2 * zeta_int2 * sign(1, n1) * sign(1, n2) &
         - tht_int3 * zeta_int3 * sign(1, m1) * sign(1, meq) * sign(1, n2) * sign(1, neq) &
         - tht_int4 * zeta_int4 * sign(1, m2) * sign(1, meq) * sign(1, n1) * sign(1, neq)
   end subroutine scs_convolve


   subroutine ccc_convolve(ans, m1, n1, m2, n2, meq, neq)
      integer, intent(in) :: m1, m2, n1, n2, meq, neq
      real(r8), intent(out) :: ans
      real(r8) :: tht_int1, tht_int2, tht_int3, tht_int4, &
         zeta_int1, zeta_int2, zeta_int3, zeta_int4
      !
      !     This subroutine calculates the 2D integral for tht and zeta running
      !     from 0 to 2*PI of:
      !     cos(m1*tht - n1*zeta)*cos(meq*tht - neq*zeta)*cos(m2*tht - n2*zeta)
      !     The value returned is actually (2/PI) times this integral.
      !
      call ccc(tht_int1, abs(m1), abs(m2), abs(meq))
      call ccc(zeta_int1, abs(n1), abs(n2), abs(neq))
      call css(tht_int2, abs(meq), abs(m1), abs(m2))
      call css(zeta_int2, abs(neq), abs(n1), abs(n2))
      call css(tht_int3, abs(m1), abs(m2), abs(meq))
      call css(zeta_int3, abs(n1), abs(n2), abs(neq))
      call css(tht_int4, abs(m2), abs(m1), abs(meq))
      call css(zeta_int4, abs(n2), abs(n1), abs(neq))
      ans = tht_int1 * zeta_int1 &
         + tht_int2 * zeta_int2 * sign(1, m1) * sign(1, m2) * sign(1, n1) * sign(1, n2) &
         + tht_int3 * zeta_int3 * sign(1, m2) * sign(1, meq) * sign(1, n2) * sign(1, neq) &
         + tht_int4 * zeta_int4 * sign(1, m1) * sign(1, meq) * sign(1, n1) * sign(1, neq)
   end subroutine ccc_convolve


   subroutine ccc(result, i, j, k)
      integer, intent(in) :: i, j, k
      real(r8), intent(out) :: result
      integer :: izeros
      !
      !     This subroutine calculates the 1D integral of
      !     cos(i*x)*cos(j*x)*cos(k*x)
      !     for x running from 0 to 2*PI.
      !     The value returned is actually (2/PI) times this integral.
      !
      izeros = sum(merge(1, 0, (/ i, j, k/) .eq. 0))

      result = 0
      if (izeros .eq. 3) then
         result = 4
      else if (izeros .eq. 1) then
         if ((abs(i) .eq. abs(j)) .or. (abs(j) .eq. abs(k)) .or. (abs(k) .eq. abs(i))) then
            result = 2
         end if
      else if (izeros .eq. 0) then
         if ((abs(k) .eq. abs(i - j)) .or. (abs(k) .eq. abs(i+j))) result = 1
      end if
   end subroutine ccc


   subroutine css(result, k, i, j)
      integer, intent(in) :: i, j, k
      real(r8) :: result
      !
      !     This subroutine calculates the 1D integral of
      !     cos(k*x)*sin(i*x)*sin(j*x)
      !     for x running from 0 to 2*PI.
      !     The value returned is actually (2/PI) times this integral.
      !
      result = 0
      if (i .eq. 0 .or. j .eq. 0) return
      if (k .eq. 0) then
         if (abs(i) .eq. abs(j)) result=2
         return
      end if
      if (abs(k) .eq. abs(i - j)) result = 1
      if (abs(k) .eq. abs(i + j)) result = -1
   end subroutine css


end module fourier_lib
