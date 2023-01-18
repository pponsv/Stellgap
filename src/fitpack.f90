module fitpack
   use kind_spec

   implicit none

contains

   subroutine curv1(n, x, y, slp1, slpn, islpsw, yp, temp, sigma, ierr)

!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer n, islpsw, ierr
      real(kind=rprec) slp1, slpn, sigma
      real(kind=rprec), dimension(n) :: x, y, yp, temp
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: nm1, np1, i, ibak
      real(kind=rprec) :: sigmap, slpp1, delx1, delx2,&
      &c1, c2, c3, slppn, delxn,&
      &delxnm, dx1, diag1, sdiag1, dx2, diag2, sdiag2, diag
!-----------------------------------------------
!
!
!                                 coded by alan kaylor cline
!                           from fitpack -- january 26, 1987
!                        a curve and surface fitting package
!                      a product of pleasant valley software
!                  8603 altus cove, austin, texas 78759, usa
!
! this subroutine determines the parameters necessary to
! compute an interpolatory spline under tension through
! a sequence of functional values. the slopes at the two
! ends of the curve may be specified or omitted.  for actual
! computation of points on the curve it is necessary to call
! the function curv2.
!
! on input--
!
!   n is the number of values to be interpolated (n.ge.2).
!
!   x is an array of the n increasing abscissae of the
!   functional values.
!
!   y is an array of the n ordinates of the values, (i. e.
!   y(k) is the functional value corresponding to x(k) ).
!
!   slp1 and slpn contain the desired values for the first
!   derivative of the curve at x(1) and x(n), respectively.
!   the user may omit values for either or both of these
!   parameters and signal this with islpsw.
!
!   islpsw contains a switch indicating which slope data
!   should be used and which should be estimated by this
!   subroutine,
!          = 0 if slp1 and slpn are to be used,
!          = 1 if slp1 is to be used but not slpn,
!          = 2 if slpn is to be used but not slp1,
!          = 3 if both slp1 and slpn are to be estimated
!              internally.
!
!   yp is an array of length at least n.
!
!   temp is an array of length at least n which is used for
!   scratch storage.
!
! and
!
!   sigma contains the tension factor. this value indicates
!   the curviness desired. if abs(sigma) is nearly zero
!   (e.g. .001) the resulting curve is approximately a
!   cubic spline. if abs(sigma) is large (e.g. 50.) the
!   resulting curve is nearly a polygonal line. if sigma
!   equals zero a cubic spline results.  a standard value
!   for sigma is approximately 1. in absolute value.
!
! on output--
!
!   yp contains the values of the second derivative of the
!   curve at the given nodes.
!
!   ierr contains an error flag,
!        = 0 for normal return,
!        = 1 if n is less than 2,
!        = 2 if x-values are not strictly increasing.
!
! and
!
!   n, x, y, slp1, slpn, islpsw and sigma are unaltered.
!
! this subroutine references package modules ceez, terms_1,
! and snhcsh_1.
!
!-----------------------------------------------------------
!
      nm1 = n - 1
      np1 = n + 1
      ierr = 0
      if (n > 1) then
         if (x(n) <= x(1)) go to 9
!
! denormalize tension factor
!
         sigmap = abs(sigma)*real(n - 1)/(x(n) - x(1))
!
! approximate end slopes
!
         if (islpsw < 2) then
            slpp1 = slp1
         else
            delx1 = x(2) - x(1)
            delx2 = delx1 + delx1
            if (n > 2) delx2 = x(3) - x(1)
            if (delx1 <= 0. .or. delx2 <= delx1) go to 9
            call ceez(delx1, delx2, sigmap, c1, c2, c3, n)
            slpp1 = c1*y(1) + c2*y(2)
            if (n > 2) slpp1 = slpp1 + c3*y(3)
         end if
         if (islpsw /= 1 .and. islpsw /= 3) then
            slppn = slpn
         else
            delxn = x(n) - x(nm1)
            delxnm = delxn + delxn
            if (n > 2) delxnm = x(n) - x(n - 2)
            if (delxn <= 0. .or. delxnm <= delxn) go to 9
            call ceez((-delxn), (-delxnm), sigmap, c1, c2, c3, n)
            slppn = c1*y(n) + c2*y(nm1)
            if (n > 2) slppn = slppn + c3*y(n - 2)
!
! set up right hand side and tridiagonal system for yp and
! perform forward elimination
!
         end if
         delx1 = x(2) - x(1)
         if (delx1 <= 0.) go to 9
         dx1 = (y(2) - y(1))/delx1
         call terms_1(diag1, sdiag1, sigmap, delx1)
         yp(1) = (dx1 - slpp1)/diag1
         temp(1) = sdiag1/diag1
         if (n /= 2) then
            do i = 2, nm1
               delx2 = x(i + 1) - x(i)
               if (delx2 <= 0.) go to 9
               dx2 = (y(i + 1) - y(i))/delx2
               call terms_1(diag2, sdiag2, sigmap, delx2)
               diag = diag1 + diag2 - sdiag1*temp(i - 1)
               yp(i) = (dx2 - dx1 - sdiag1*yp(i - 1))/diag
               temp(i) = sdiag2/diag
               dx1 = dx2
               diag1 = diag2
               sdiag1 = sdiag2
            end do
         end if
         diag = diag1 - sdiag1*temp(nm1)
         yp(n) = (slppn - dx1 - sdiag1*yp(nm1))/diag
!
! perform back substitution
!
         do i = 2, n
            ibak = np1 - i
            yp(ibak) = yp(ibak) - temp(ibak)*yp(ibak + 1)
         end do
         return
!
! too few points
!
      end if
      ierr = 1
      return
!
! x-values not strictly increasing
!
9     continue
      ierr = 2
      return
   end subroutine curv1
   subroutine curvs(n, x, y, d, isw, s, eps, ys, ysp, sigma, temp, ierr)

!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer n, isw, ierr
      real(kind=rprec) s, eps, sigma
      real(kind=rprec), dimension(n) :: x, y, d, ys, ysp
      real(kind=rprec), dimension(n, 9) :: temp
!-----------------------------------------------
!
!
!                                 coded by alan kaylor cline
!                           from fitpack -- january 26, 1987
!                        a curve and surface fitting package
!                      a product of pleasant valley software
!                  8603 altus cove, austin, texas 78759, usa
!
! this subroutine determines the parameters necessary to
! compute a smoothing spline under tension. for a given
! increasing sequence of abscissae (x(i)), i = 1,..., n and
! associated ordinates (y(i)), i = 1,..., n, the function
! determined minimizes the summation from i = 1 to n-1 of
! the square of the second derivative of f plus sigma
! squared times the difference of the first derivative of f
! and (f(x(i+1))-f(x(i)))/(x(i+1)-x(i)) squared, over all
! functions f with two continuous derivatives such that the
! summation of the square of (f(x(i))-y(i))/d(i) is less
! than or equal to a given constant s, where (d(i)), i = 1,
! ..., n are a given set of observation weights. the
! function determined is a spline under tension with third
! derivative discontinuities at (x(i)), i = 2,..., n-1. for
! actual computation of points on the curve it is necessary
! to call the function curv2. the determination of the curve
! is performed by subroutine curvss, the subroutine curvs
! only decomposes the workspace for curvss.
!
! on input--
!
!   n is the number of values to be smoothed (n.ge.2).
!
!   x is an array of the n increasing abscissae of the
!   values to be smoothed.
!
!   y is an array of the n ordinates of the values to be
!   smoothed, (i. e. y(k) is the functional value
!   corresponding to x(k) ).
!
!   d is a parameter containing the observation weights.
!   this may either be an array of length n or a scalar
!   (interpreted as a constant). the value of d
!   corresponding to the observation (x(k),y(k)) should
!   be an approximation to the standard deviation of error.
!
!   isw contains a switch indicating whether the parameter
!   d is to be considered a vector or a scalar,
!          = 0 if d is an array of length n,
!          = 1 if d is a scalar.
!
!   s contains the value controlling the smoothing. this
!   must be non-negative. for s equal to zero, the
!   subroutine does interpolation, larger values lead to
!   smoother funtions. if parameter d contains standard
!   deviation estimates, a reasonable value for s is
!   real(n).
!
!   eps contains a tolerance on the relative precision to
!   which s is to be interpreted. this must be greater than
!   or equal to zero and less than or equal to one. a
!   reasonable value for eps is sqrt(2./real(n)).
!
!   ys is an array of length at least n.
!
!   ysp is an array of length at least n.
!
!   sigma contains the tension factor. this value indicates
!   the degree to which the first derivative part of the
!   smoothing functional is emphasized. if sigma is nearly
!   zero (e. g. .001) the resulting curve is approximately a
!   cubic spline. if sigma is large (e. g. 50.) the
!   resulting curve is nearly a polygonal line. if sigma
!   equals zero a cubic spline results. a standard value for
!   sigma is approximately 1.
!
! and
!
!   temp is an array of length at least 9*n which is used
!   for scratch storage.
!
! on output--
!
!   ys contains the smoothed ordinate values.
!
!   ysp contains the values of the second derivative of the
!   smoothed curve at the given nodes.
!
!   ierr contains an error flag,
!        = 0 for normal return,
!        = 1 if n is less than 2,
!        = 2 if s is negative,
!        = 3 if eps is negative or greater than one,
!        = 4 if x-values are not strictly increasing,
!        = 5 if a d-value is non-positive.
!
! and
!
!   n, x, y, d, isw, s, eps, and sigma are unaltered.
!
! this subroutine references package modules curvss, terms_1,
! and snhcsh_1.
!
!-----------------------------------------------------------
!
! decompose temp into nine arrays and call curvss
!
      call curvss(n, x, y, d, isw, s, eps, ys, ysp, sigma, temp(1, 1),&
      &temp(1, 2), temp(1, 3), temp(1, 4), temp(1, 5), temp(1, 6), temp(1, 7&
      &), temp(1, 8), temp(1, 9), ierr)
      return
   end subroutine curvs
   function curv2(t, n, x, y, yp, sigma)

!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer n
      real(kind=rprec) t, sigma, curv2
      real(kind=rprec), dimension(n) :: x, y, yp
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: im1, i
      real(kind=rprec) :: sigmap, del1, del2, dels, sum, sigdel,&
      &ss, dummy, s1, s2
!-----------------------------------------------
!   E x t e r n a l   F u n c t i o n s
!-----------------------------------------------
!-----------------------------------------------
!
!
!                                 coded by alan kaylor cline
!                           from fitpack -- january 26, 1987
!                        a curve and surface fitting package
!                      a product of pleasant valley software
!                  8603 altus cove, austin, texas 78759, usa
!
! this function interpolates a curve at a given point
! using a spline under tension. the subroutine curv1 should
! be called earlier to determine certain necessary
! parameters.
!
! on input--
!
!   t contains a real value to be mapped onto the interpo-
!   lating curve.
!
!   n contains the number of points which were specified to
!   determine the curve.
!
!   x and y are arrays containing the abscissae and
!   ordinates, respectively, of the specified points.
!
!   yp is an array of second derivative values of the curve
!   at the nodes.
!
! and
!
!   sigma contains the tension factor (its sign is ignored).
!
! the parameters n, x, y, yp, and sigma should be input
! unaltered from the output of curv1.
!
! on output--
!
!   curv2 contains the interpolated value.
!
! none of the input parameters are altered.
!
! this function references package modules intrvl_1 and
! snhcsh_1.
!
!-----------------------------------------------------------
!
! determine interval
!
      im1 = intrvl_1(t, x, n)
      i = im1 + 1
!
! denormalize tension factor
!
      sigmap = abs(sigma)*real(n - 1)/(x(n) - x(1))
!
! set up and perform interpolation
!
      del1 = t - x(im1)
      del2 = x(i) - t
      dels = x(i) - x(im1)
      sum = (y(i)*del1 + y(im1)*del2)/dels
      if (sigmap == 0.) then
         curv2 = sum - del1*del2*(yp(i)*(del1 + dels) + yp(im1)*(del2 + dels))&
         &/(6.*dels)
         return
      end if
      sigdel = sigmap*dels
      call snhcsh_1(ss, dummy, sigdel, -1)
      call snhcsh_1(s1, dummy, sigmap*del1, -1)
      call snhcsh_1(s2, dummy, sigmap*del2, -1)
      curv2 = sum + (yp(i)*del1*(s1 - ss) + yp(im1)*del2*(s2 - ss))/(sigdel*&
      &sigmap*(1.+ss))
      return
   end function curv2
   function curvd(t, n, x, y, yp, sigma)

!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer n
      real(kind=rprec) t, sigma, curvd
      real(kind=rprec), dimension(n) :: x, y, yp
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: im1, i
      real(kind=rprec) :: sigmap, del1, del2, dels, sum, sigdel,&
      &ss, dummy, c1, c2
!-----------------------------------------------
!   E x t e r n a l   F u n c t i o n s
!-----------------------------------------------
!-----------------------------------------------
!
!
!                                 coded by alan kaylor cline
!                           from fitpack -- january 26, 1987
!                        a curve and surface fitting package
!                      a product of pleasant valley software
!                  8603 altus cove, austin, texas 78759, usa
!
! this function differentiates a curve at a given point
! using a spline under tension. the subroutine curv1 should
! be called earlier to determine certain necessary
! parameters.
!
! on input--
!
!   t contains a real value at which the derivative is to be
!   determined.
!
!   n contains the number of points which were specified to
!   determine the curve.
!
!   x and y are arrays containing the abscissae and
!   ordinates, respectively, of the specified points.
!
!   yp is an array of second derivative values of the curve
!   at the nodes.
!
! and
!
!   sigma contains the tension factor (its sign is ignored).
!
! the parameters n, x, y, yp, and sigma should be input
! unaltered from the output of curv1.
!
! on output--
!
!   curvd contains the derivative value.
!
! none of the input parameters are altered.
!
! this function references package modules intrvl_1 and
! snhcsh_1.
!
!-----------------------------------------------------------
!
! determine interval
!
      im1 = intrvl_1(t, x, n)
      i = im1 + 1
!
! denormalize tension factor
!
      sigmap = abs(sigma)*real(n - 1)/(x(n) - x(1))
!
! set up and perform differentiation
!
      del1 = t - x(im1)
      del2 = x(i) - t
      dels = x(i) - x(im1)
      sum = (y(i) - y(im1))/dels
      if (sigmap == 0.) then
         curvd = sum + (yp(i)*(2.*del1*del1 - del2*(del1 + dels)) - yp(im1)*(&
         &2.*del2*del2 - del1*(del2 + dels)))/(6.*dels)
         return
      end if
      sigdel = sigmap*dels
      call snhcsh_1(ss, dummy, sigdel, -1)
      call snhcsh_1(dummy, c1, sigmap*del1, 1)
      call snhcsh_1(dummy, c2, sigmap*del2, 1)
      curvd = sum + (yp(i)*(c1 - ss) - yp(im1)*(c2 - ss))/(sigdel*sigmap*(1.&
      &+ss))
      return
   end function curvd
   function curvi(xl, xu, n, x, y, yp, sigma)

!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer n
      real(kind=rprec) xl, xu, sigma, curvi
      real(kind=rprec), dimension(n) :: x, y, yp
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: ilm1, il, ium1, iu, ilp1, i
      real(kind=rprec) :: sigmap, xxl, xxu, ssign, sum, del1,&
      &del2, dels, t1, t2,&
      &dummy, c1, c2, ss, cs, delu1, delu2, dell1, dell2, deli, cu1,&
      &cu2, cl1, cl2
!-----------------------------------------------
!   E x t e r n a l   F u n c t i o n s
!-----------------------------------------------
!-----------------------------------------------
!
!
!                                 coded by alan kaylor cline
!                           from fitpack -- january 26, 1987
!                        a curve and surface fitting package
!                      a product of pleasant valley software
!                  8603 altus cove, austin, texas 78759, usa
!
! this function integrates a curve specified by a spline
! under tension between two given limits. the subroutine
! curv1 should be called earlier to determine necessary
! parameters.
!
! on input--
!
!   xl and xu contain the upper and lower limits of inte-
!   gration, respectively. (sl need not be less than or
!   equal to xu, curvi (xl,xu,...) .eq. -curvi (xu,xl,...) ).
!
!   n contains the number of points which were specified to
!   determine the curve.
!
!   x and y are arrays containing the abscissae and
!   ordinates, respectively, of the specified points.
!
!   yp is an array from subroutine curv1 containing
!   the values of the second derivatives at the nodes.
!
! and
!
!   sigma contains the tension factor (its sign is ignored).
!
! the parameters n, x, y, yp, and sigma should be input
! unaltered from the output of curv1.
!
! on output--
!
!   curvi contains the integral value.
!
! none of the input parameters are altered.
!
! this function references package modules intrvl_1 and
! snhcsh_1.
!
!-----------------------------------------------------------
!
! denormalize tension factor
!
      sigmap = abs(sigma)*real(n - 1)/(x(n) - x(1))
!
! determine actual upper and lower bounds
!
      xxl = xl
      xxu = xu
      ssign = 1.
      if (xl >= xu) then
         xxl = xu
         xxu = xl
         ssign = -1.
         if (xl <= xu) then
!
! return zero if xl .eq. xu
!
            curvi = 0.
            return
!
! search for proper intervals
!
         end if
      end if
      ilm1 = intrvl_1(xxl, x, n)
      il = ilm1 + 1
      ium1 = intrvl_1(xxu, x, n)
      iu = ium1 + 1
      if (il /= iu) then
!
! integrate from xxl to x(il)
!
         sum = 0.
         if (xxl /= x(il)) then
            del1 = xxl - x(ilm1)
            del2 = x(il) - xxl
            dels = x(il) - x(ilm1)
            t1 = (del1 + dels)*del2/(2.*dels)
            t2 = del2*del2/(2.*dels)
            sum = t1*y(il) + t2*y(ilm1)
            if (sigma /= 0.) then
               call snhcsh_1(dummy, c1, sigmap*del1, 2)
               call snhcsh_1(dummy, c2, sigmap*del2, 2)
               call snhcsh_1(ss, cs, sigmap*dels, 3)
               sum = sum + ((dels*dels*(cs - ss/2.) - del1*del1*(c1 - ss&
               &/2.))*yp(il) + del2*del2*(c2 - ss/2.)*yp(ilm1))/(sigmap*&
               &sigmap*dels*(1.+ss))
            else
               sum = sum - t1*t1*dels*yp(il)/6.-t2*(del1*(del2 + dels)&
               &+ dels*dels)*yp(ilm1)/12.
!
! integrate over interior intervals
!
            end if
         end if
         if (iu - il /= 1) then
            ilp1 = il + 1
            do i = ilp1, ium1
               dels = x(i) - x(i - 1)
               sum = sum + (y(i) + y(i - 1))*dels/2.
               if (sigma /= 0.) then
                  call snhcsh_1(ss, cs, sigmap*dels, 3)
                  sum = sum + (yp(i) + yp(i - 1))*dels*(cs - ss/2.)/(sigmap*&
                  &sigmap*(1.+ss))
               else
                  sum = sum - (yp(i) + yp(i - 1))*dels*dels*dels/24.
               end if
            end do
!
! integrate from x(iu-1) to xxu
!
         end if
         if (xxu == x(ium1)) go to 10
         del1 = xxu - x(ium1)
         del2 = x(iu) - xxu
         dels = x(iu) - x(ium1)
         t1 = del1*del1/(2.*dels)
         t2 = (del2 + dels)*del1/(2.*dels)
         sum = sum + t1*y(iu) + t2*y(ium1)
         if (sigma /= 0.) then
            call snhcsh_1(dummy, c1, sigmap*del1, 2)
            call snhcsh_1(dummy, c2, sigmap*del2, 2)
            call snhcsh_1(ss, cs, sigmap*dels, 3)
            sum = sum + (yp(iu)*del1*del1*(c1 - ss/2.) + yp(ium1)*(dels*dels&
            &*(cs - ss/2.) - del2*del2*(c2 - ss/2.)))/(sigmap*sigmap*dels*(&
            &1.+ss))
            go to 10
         end if
         sum = sum - t1*(del2*(del1 + dels) + dels*dels)*yp(iu)/12.-t2&
         &*t2*dels*yp(ium1)/6.
         go to 10
!
! integrate from xxl to xxu
!
      end if
      delu1 = xxu - x(ium1)
      delu2 = x(iu) - xxu
      dell1 = xxl - x(ium1)
      dell2 = x(iu) - xxl
      dels = x(iu) - x(ium1)
      deli = xxu - xxl
      t1 = (delu1 + dell1)*deli/(2.*dels)
      t2 = (delu2 + dell2)*deli/(2.*dels)
      sum = t1*y(iu) + t2*y(ium1)
      if (sigma /= 0.) then
         call snhcsh_1(dummy, cu1, sigmap*delu1, 2)
         call snhcsh_1(dummy, cu2, sigmap*delu2, 2)
         call snhcsh_1(dummy, cl1, sigmap*dell1, 2)
         call snhcsh_1(dummy, cl2, sigmap*dell2, 2)
         call snhcsh_1(ss, dummy, sigmap*dels, -1)
         sum = sum + (yp(iu)*(delu1*delu1*(cu1 - ss/2.) - dell1*dell1*(cl1 -&
         &ss/2.)) + yp(ium1)*(dell2*dell2*(cl2 - ss/2.) - delu2*delu2*(cu2 -&
         &ss/2.)))/(sigmap*sigmap*dels*(1.+ss))
      else
         sum = sum - t1*(delu2*(dels + delu1) + dell2*(dels + dell1))*yp&
         &(iu)/12.-t2*(dell1*(dels + dell2) + delu1*(dels + delu2))*&
         &yp(ium1)/12.
!
! correct sign and return
!
      end if
10    continue
      curvi = ssign*sum
      return
   end function curvi
   subroutine curvp1(n, x, y, p, yp, temp, sigma, ierr)

!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer n, ierr
      real(kind=rprec) p, sigma
      real(kind=rprec), dimension(n) :: x, y, yp
      real(kind=rprec), dimension(1) :: temp
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: nm1, np1, i, npi, ibak, npibak
      real(kind=rprec) :: sigmap, delx1, dx1, diag1, sdiag1,&
      &delx2, dx2, diag2,&
      &sdiag2, diag, ypn
!-----------------------------------------------
!
!
!                                 coded by alan kaylor cline
!                           from fitpack -- january 26, 1987
!                        a curve and surface fitting package
!                      a product of pleasant valley software
!                  8603 altus cove, austin, texas 78759, usa
!
! this subroutine determines the parameters necessary to
! compute a periodic interpolatory spline under tension
! through a sequence of functional values. for actual ends
! of the curve may be specified or omitted.  for actual
! computation of points on the curve it is necessary to call
! the function curvp2.
!
! on input--
!
!   n is the number of values to be interpolated (n.ge.2).
!
!   x is an array of the n increasing abscissae of the
!   functional values.
!
!   y is an array of the n ordinates of the values, (i. e.
!   y(k) is the functional value corresponding to x(k) ).
!
!   p is the period (p .gt. x(n)-x(1)).
!
!   yp is an array of length at least n.
!
!   temp is an array of length at least 2*n which is used
!   for scratch storage.
!
! and
!
!   sigma contains the tension factor.  this value indicates
!   the curviness desired. if abs(sigma) is nearly zero
!   (e.g. .001) the resulting curve is approximately a
!   cubic spline. if abs(sigma) is large (e.g. 50.) the
!   resulting curve is nearly a polygonal line. if sigma
!   equals zero a cubic spline results.  a standard value
!   for sigma is approximately 1. in absolute value.
!
! on output--
!
!   yp contains the values of the second derivative of the
!   curve at the given nodes.
!
!   ierr contains an error flag,
!        = 0 for normal return,
!        = 1 if n is less than 2,
!        = 2 if p is less than or equal to x(n)-x(1),
!        = 3 if x-values are not strictly increasing.
!
! and
!
!  n, x, y, and sigma are unaltered.
!
! this subroutine references package modules terms_1 and
! snhcsh_1.
!
!-----------------------------------------------------------
!
      nm1 = n - 1
      np1 = n + 1
      ierr = 0
      if (n > 1) then
         if (p <= x(n) - x(1) .or. p <= 0.) go to 7
!
! denormalize tension factor
!
         sigmap = abs(sigma)*real(n)/p
!
! set up right hand side and tridiagonal system for yp and
! perform forward elimination
!
         delx1 = p - (x(n) - x(1))
         dx1 = (y(1) - y(n))/delx1
         call terms_1(diag1, sdiag1, sigmap, delx1)
         delx2 = x(2) - x(1)
         if (delx2 <= 0.) go to 8
         dx2 = (y(2) - y(1))/delx2
         call terms_1(diag2, sdiag2, sigmap, delx2)
         diag = diag1 + diag2
         yp(1) = (dx2 - dx1)/diag
         temp(np1) = -sdiag1/diag
         temp(1) = sdiag2/diag
         dx1 = dx2
         diag1 = diag2
         sdiag1 = sdiag2
         if (n /= 2) then
            do i = 2, nm1
               npi = n + i
               delx2 = x(i + 1) - x(i)
               if (delx2 <= 0.) go to 8
               dx2 = (y(i + 1) - y(i))/delx2
               call terms_1(diag2, sdiag2, sigmap, delx2)
               diag = diag1 + diag2 - sdiag1*temp(i - 1)
               yp(i) = (dx2 - dx1 - sdiag1*yp(i - 1))/diag
               temp(npi) = -temp(npi - 1)*sdiag1/diag
               temp(i) = sdiag2/diag
               dx1 = dx2
               diag1 = diag2
               sdiag1 = sdiag2
            end do
         end if
         delx2 = p - (x(n) - x(1))
         dx2 = (y(1) - y(n))/delx2
         call terms_1(diag2, sdiag2, sigmap, delx2)
         yp(n) = dx2 - dx1
         temp(nm1) = temp(2*n - 1) - temp(nm1)
         if (n /= 2) then
!
! perform first step of back substitution
!
            do i = 3, n
               ibak = np1 - i
               npibak = n + ibak
               yp(ibak) = yp(ibak) - temp(ibak)*yp(ibak + 1)
               temp(ibak) = temp(npibak) - temp(ibak)*temp(ibak + 1)
            end do
         end if
         yp(n) = (yp(n) - sdiag2*yp(1) - sdiag1*yp(nm1))/(diag1 + diag2 +&
         &sdiag2*temp(1) + sdiag1*temp(nm1))
!
! perform second step of back substitution
!
         ypn = yp(n)
         yp(:nm1) = yp(:nm1) + temp(:nm1)*ypn
         return
!
! too few points
!
      end if
      ierr = 1
      return
!
! period too small
!
7     continue
      ierr = 2
      return
!
! x-values not strictly increasing
!
8     continue
      ierr = 3
      return
   end subroutine curvp1
   subroutine curvps(n, x, y, p, d, isw, s, eps, ys, ysp, sigma, temp, ierr)

!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer n, isw, ierr
      real(kind=rprec) p, s, eps, sigma
      real(kind=rprec), dimension(n) :: x, y, d, ys, ysp
      real(kind=rprec), dimension(n, 11) :: temp
!-----------------------------------------------
!
!
!                                 coded by alan kaylor cline
!                           from fitpack -- january 26, 1987
!                        a curve and surface fitting package
!                      a product of pleasant valley software
!                  8603 altus cove, austin, texas 78759, usa
!
! this subroutine determines the parameters necessary to
! compute a periodic smoothing spline under tension. for a
! given increasing sequence of abscissae (x(i)), i = 1,...,n
! and associated ordinates (y(i)), i = 1,...,n, letting p be
! the period, x(n+1) = x(1)+p, and y(n+1) = y(1), the
! function determined minimizes the summation from i = 1 to
! n of the square of the second derivative of f plus sigma
! squared times the difference of the first derivative of f
! and (f(x(i+1))-f(x(i)))/(x(i+1)-x(i)) squared, over all
! functions f with period p and two continuous derivatives
! such that the summation of the square of
! (f(x(i))-y(i))/d(i) is less than or equal to a given
! constant s, where (d(i)), i = 1,...,n are a given set of
! observation weights. the function determined is a periodic
! spline under tension with third derivative discontinuities
! at (x(i)) i = 1,...,n (and all periodic translations of
! these values). for actual computation of points on the
! curve it is necessary to call the function curvp2. the
! determination of the curve is performed by subroutine
! curvpp, the subroutin curvps only decomposes the workspace
! for curvpp.
!
! on input--
!
!   n is the number of values to be smoothed (n.ge.2).
!
!   x is an array of the n increasing abscissae of the
!   values to be smoothed.
!
!   y is an array of the n ordinates of the values to be
!   smoothed, (i. e. y(k) is the functional value
!   corresponding to x(k) ).
!
!   p is the period (p .gt. x(n)-x(1)).
!
!   d is a parameter containing the observation weights.
!   this may either be an array of length n or a scalar
!   (interpreted as a constant). the value of d
!   corresponding to the observation (x(k),y(k)) should
!   be an approximation to the standard deviation of error.
!
!   isw contains a switch indicating whether the parameter
!   d is to be considered a vector or a scalar,
!          = 0 if d is an array of length n,
!          = 1 if d is a scalar.
!
!   s contains the value controlling the smoothing. this
!   must be non-negative. for s equal to zero, the
!   subroutine does interpolation, larger values lead to
!   smoother funtions. if parameter d contains standard
!   deviation estimates, a reasonable value for s is
!   real(n).
!
!   eps contains a tolerance on the relative precision to
!   which s is to be interpreted. this must be greater than
!   or equal to zero and less than or equal to one. a
!   reasonable value for eps is sqrt(2./real(n)).
!
!   ys is an array of length at least n.
!
!   ysp is an array of length at least n.
!
!   sigma contains the tension factor. this value indicates
!   the degree to which the first derivative part of the
!   smoothing functional is emphasized. if sigma is nearly
!   zero (e. g. .001) the resulting curve is approximately a
!   cubic spline. if sigma is large (e. g. 50.) the
!   resulting curve is nearly a polygonal line. if sigma
!   equals zero a cubic spline results. a standard value for
!   sigma is approximately 1.
!
! and
!
!   temp is an array of length at least 11*n which is used
!   for scratch storage.
!
! on output--
!
!   ys contains the smoothed ordinate values.
!
!   ysp contains the values of the second derivative of the
!   smoothed curve at the given nodes.
!
!   ierr contains an error flag,
!        = 0 for normal return,
!        = 1 if n is less than 2,
!        = 2 if s is negative,
!        = 3 if eps is negative or greater than one,
!        = 4 if x-values are not strictly increasing,
!        = 5 if a d-value is non-positive,
!        = 6 if p is less than or equal to x(n)-x(1).
!
! and
!
!   n, x, y, p, d, isw, s, eps, and sigma are unaltered.
!
! this subroutine references package modules curvpp, terms_1,
! and snhcsh_1.
!
!-----------------------------------------------------------
!
! decompose temp into eleven arrays and call curvpp
!
      call curvpp(n, x, y, p, d, isw, s, eps, ys, ysp, sigma, temp(1, 1)&
      &, temp(1, 2), temp(1, 3), temp(1, 4), temp(1, 5), temp(1, 6), temp(1&
      &, 7), temp(1, 8), temp(1, 9), temp(1, 10), temp(1, 11), ierr)
      return
   end subroutine curvps
   function curvp2(t, n, x, y, p, yp, sigma)

!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer n
      real(kind=rprec) t, p, sigma, curvp2
      real(kind=rprec), dimension(n) :: x, y, yp
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: im1, i
      real(kind=rprec)::tp, sigmap, del1, del2, dels, sum,&
      &sigdel, ss, dummy, s1, s2
!-----------------------------------------------
!   E x t e r n a l   F u n c t i o n s
!-----------------------------------------------
!-----------------------------------------------
!
!
!                                 coded by alan kaylor cline
!                           from fitpack -- january 26, 1987
!                        a curve and surface fitting package
!                      a product of pleasant valley software
!                  8603 altus cove, austin, texas 78759, usa
!
! this function interpolates a curve at a given point using
! a periodic spline under tension. the subroutine curvp1
! should be called earlier to determine certain necessary
! parameters.
!
! on input--
!
!   t contains a real value to be mapped onto the interpo-
!   lating curve.
!
!   n contains the number of points which were specified to
!   determine the curve.
!
!   x and y are arrays containing the abscissae and
!   ordinates, respectively, of the specified points.
!
!   p contains the period.
!
!   yp is an array of second derivative values of the curve
!   at the nodes.
!
! and
!
!   sigma contains the tension factor (its sign is ignored).
!
! the parameters n, x, y, p, yp, and sigma should be input
! unaltered from the output of curvp1.
!
! on output--
!
!   curvp2 contains the interpolated value.
!
! none of the input parameters are altered.
!
! this function references package modules intrvp and
! snhcsh_1.
!
!-----------------------------------------------------------
!
! determine interval
!
      im1 = intrvp(t, x, n, p, tp)
      i = im1 + 1
!
! denormalize tension factor
!
      sigmap = abs(sigma)*real(n)/p
!
! set up and perform interpolation
!
      del1 = tp - x(im1)
      if (im1 /= n) then
         del2 = x(i) - tp
         dels = x(i) - x(im1)
      else
         i = 1
         del2 = x(1) + p - tp
         dels = p - (x(n) - x(1))
      end if
      sum = (y(i)*del1 + y(im1)*del2)/dels
      if (sigmap == 0.) then
         curvp2 = sum - del1*del2*(yp(i)*(del1 + dels) + yp(im1)*(del2 + dels)&
         &)/(6.*dels)
         return
      end if
      sigdel = sigmap*dels
      call snhcsh_1(ss, dummy, sigdel, -1)
      call snhcsh_1(s1, dummy, sigmap*del1, -1)
      call snhcsh_1(s2, dummy, sigmap*del2, -1)
      curvp2 = sum + (yp(i)*del1*(s1 - ss) + yp(im1)*del2*(s2 - ss))/(sigdel*&
      &sigmap*(1.+ss))
      return
   end function curvp2
   function curvpi(xl, xu, n, x, y, p, yp, sigma)

!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer n
      real(kind=rprec) xl, xu, p, sigma, curvpi
      real(kind=rprec), dimension(n) :: x, y, yp
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: uper, isign, ilm1, lper, ium1, ideltp, isave, il, iu, i&
      &, ilp1, np1, iup1, ii, im1
      real(kind=rprec) :: sigmap, x1pp, xxl, xxu, xsave, xil,&
      &xiu, s1, dels, ss, cs&
      &, s2, del1, del2, t1, t2, dummy, c1, c2, s3, s4, s5, s6, s7, s8&
      &, delu1, delu2, dell1, dell2, deli, cu1, cu2, cl1, cl2, so, si
      logical :: bdy
!-----------------------------------------------
!   E x t e r n a l   F u n c t i o n s
!-----------------------------------------------
!-----------------------------------------------
!
!
!                                 coded by alan kaylor cline
!                           from fitpack -- january 26, 1987
!                        a curve and surface fitting package
!                      a product of pleasant valley software
!                  8603 altus cove, austin, texas 78759, usa
!
! this function integrates a curve specified by a periodic
! spline under tension between two given limits. the
! subroutine curvp1 should be called earlier to determine
! necessary parameters.
!
! on input--
!
!   xl and xu contain the upper and lower limits of inte-
!   gration, respectively. (sl need not be less than or
!   equal to xu, curvpi (xl,xu,...) .eq. -curvpi (xu,xl,...) ).
!
!   n contains the number of points which were specified to
!   determine the curve.
!
!   x and y are arrays containing the abscissae and
!   ordinates, respectively, of the specified points.
!
!   p contains the period.
!
!   yp is an array from subroutine curvp1 containing
!   the values of the second derivatives at the nodes.
!
! and
!
!   sigma contains the tension factor (its sign is ignored).
!
! the parameters n, x, y, p, yp, and sigma should be input
! unaltered from the output of curvp1.
!
! on output--
!
!
!   curvpi contains the integral value.
!
! none of the input parameters are altered.
!
! this function references package modules intrvp and
! snhcsh_1.
!
!--------------------------------------------------------------
!
!
! denormalize tension factor
!
      sigmap = abs(sigma)*real(n)/p
!
! determine actual upper and lower bounds
!
      x1pp = x(1) + p
      isign = 1
      ilm1 = intrvp(xl, x, n, p, xxl)
      lper = int((xl - x(1))/p)
      if (xl < x(1)) lper = lper - 1
      ium1 = intrvp(xu, x, n, p, xxu)
      uper = int((xu - x(1))/p)
      if (xu < x(1)) uper = uper - 1
      ideltp = uper - lper
      bdy = real(ideltp)*(xxu - xxl) < 0.
      if (ideltp == 0 .and. xxu < xxl .or. ideltp < 0) isign = -1
      if (bdy) ideltp = ideltp - isign
      if (xxu < xxl) then
         xsave = xxl
         xxl = xxu
         xxu = xsave
         isave = ilm1
         ilm1 = ium1
         ium1 = isave
      end if
      il = ilm1 + 1
      if (ilm1 == n) il = 1
      xil = x(il)
      if (ilm1 == n) xil = x1pp
      iu = ium1 + 1
      if (ium1 == n) iu = 1
      xiu = x(iu)
      if (ium1 == n) xiu = x1pp
      s1 = 0.
      if (.not. (ilm1 == 1 .or. ideltp == 0 .and. .not. bdy)) then
!
! integrate from x(1) to x(ilm1), store in s1
!
         do i = 2, ilm1
            dels = x(i) - x(i - 1)
            s1 = s1 + (y(i) + y(i - 1))*dels/2.
            if (sigma /= 0.) then
               call snhcsh_1(ss, cs, sigmap*dels, 3)
               s1 = s1 + (yp(i) + yp(i - 1))*dels*(cs - ss/2.)/(sigmap*&
               &sigmap*(1.+ss))
            else
               s1 = s1 - (yp(i) + yp(i - 1))*dels*dels*dels/24.
            end if
         end do
      end if
      s2 = 0.
      if (.not. (x(ilm1) >= xxl .or. ideltp == 0 .and. .not. bdy)) then
!
! integrate from x(ilm1) to xxl, store in s2
!
         del1 = xxl - x(ilm1)
         del2 = xil - xxl
         dels = xil - x(ilm1)
         t1 = del1*del1/(2.*dels)
         t2 = (del2 + dels)*del1/(2.*dels)
         s2 = t1*y(il) + t2*y(ilm1)
         if (sigma /= 0.) then
            call snhcsh_1(dummy, c1, sigmap*del1, 2)
            call snhcsh_1(dummy, c2, sigmap*del2, 2)
            call snhcsh_1(ss, cs, sigmap*dels, 3)
            s2 = s2 + (yp(il)*del1*del1*(c1 - ss/2.) + yp(ilm1)*(dels*dels*(&
            &cs - ss/2.) - del2*del2*(c2 - ss/2.)))/(sigmap*sigmap*dels*(1.&
            &+ss))
         else
            s2 = s2 - t1*(del2*(del1 + dels) + dels*dels)*yp(il)/12.-&
            &t2*t2*dels*yp(ilm1)/6.
         end if
      end if
      s3 = 0.
      if (.not. (xxl >= xil .or. ideltp == 0 .and. bdy .or. ilm1 == ium1)) then
!
! integrate from xxl to xil, store in s3
!
         del1 = xxl - x(ilm1)
         del2 = xil - xxl
         dels = xil - x(ilm1)
         t1 = (del1 + dels)*del2/(2.*dels)
         t2 = del2*del2/(2.*dels)
         s3 = t1*y(il) + t2*y(ilm1)
         if (sigma /= 0.) then
            call snhcsh_1(dummy, c1, sigmap*del1, 2)
            call snhcsh_1(dummy, c2, sigmap*del2, 2)
            call snhcsh_1(ss, cs, sigmap*dels, 3)
            s3 = s3 + ((dels*dels*(cs - ss/2.) - del1*del1*(c1 - ss/2.))&
            &*yp(il) + del2*del2*(c2 - ss/2.)*yp(ilm1))/(sigmap*sigmap*&
            &dels*(1.+ss))
         else
            s3 = s3 - t1*t1*dels*yp(il)/6.-t2*(del1*(del2 + dels) +&
            &dels*dels)*yp(ilm1)/12.
         end if
      end if
      s4 = 0.
      if (.not. (ilm1 >= ium1 - 1 .or. ideltp == 0 .and. bdy)) then
!
! integrate from xil to x(ium1), store in s4
!
         ilp1 = il + 1
         do i = ilp1, ium1
            dels = x(i) - x(i - 1)
            s4 = s4 + (y(i) + y(i - 1))*dels/2.
            if (sigma /= 0.) then
               call snhcsh_1(ss, cs, sigmap*dels, 3)
               s4 = s4 + (yp(i) + yp(i - 1))*dels*(cs - ss/2.)/(sigmap*&
               &sigmap*(1.+ss))
            else
               s4 = s4 - (yp(i) + yp(i - 1))*dels*dels*dels/24.
            end if
         end do
      end if
      s5 = 0.
      if (.not. (x(ium1) >= xxu .or. ideltp == 0 .and. bdy .or. ilm1 == ium1)) then
!
! integrate from x(ium1) to xxu, store in s5
!
         del1 = xxu - x(ium1)
         del2 = xiu - xxu
         dels = xiu - x(ium1)
         t1 = del1*del1/(2.*dels)
         t2 = (del2 + dels)*del1/(2.*dels)
         s5 = t1*y(iu) + t2*y(ium1)
         if (sigma /= 0.) then
            call snhcsh_1(dummy, c1, sigmap*del1, 2)
            call snhcsh_1(dummy, c2, sigmap*del2, 2)
            call snhcsh_1(ss, cs, sigmap*dels, 3)
            s5 = s5 + (yp(iu)*del1*del1*(c1 - ss/2.) + yp(ium1)*(dels*dels*(&
            &cs - ss/2.) - del2*del2*(c2 - ss/2.)))/(sigmap*sigmap*dels*(1.&
            &+ss))
         else
            s5 = s5 - t1*(del2*(del1 + dels) + dels*dels)*yp(iu)/12.-&
            &t2*t2*dels*yp(ium1)/6.
         end if
      end if
      s6 = 0.
      if (.not. (xxu >= xiu .or. ideltp == 0 .and. .not. bdy)) then
!
! integrate from xxu to xiu, store in s6
!
         del1 = xxu - x(ium1)
         del2 = xiu - xxu
         dels = xiu - x(ium1)
         t1 = (del1 + dels)*del2/(2.*dels)
         t2 = del2*del2/(2.*dels)
         s6 = t1*y(iu) + t2*y(ium1)
         if (sigma /= 0.) then
            call snhcsh_1(dummy, c1, sigmap*del1, 2)
            call snhcsh_1(dummy, c2, sigmap*del2, 2)
            call snhcsh_1(ss, cs, sigmap*dels, 3)
            s6 = s6 + ((dels*dels*(cs - ss/2.) - del1*del1*(c1 - ss/2.))&
            &*yp(iu) + del2*del2*(c2 - ss/2.)*yp(ium1))/(sigmap*sigmap*&
            &dels*(1.+ss))
         else
            s6 = s6 - t1*t1*dels*yp(iu)/6.-t2*(del1*(del2 + dels) +&
            &dels*dels)*yp(ium1)/12.
         end if
      end if
      s7 = 0.
      if (.not. (iu == 1 .or. ideltp == 0 .and. .not. bdy)) then
!
! integrate from xiu to x1pp, store in s7
!
         np1 = n + 1
         iup1 = iu + 1
         do ii = iup1, np1
            im1 = ii - 1
            i = ii
            if (i == np1) i = 1
            dels = x(i) - x(im1)
            if (dels <= 0.) dels = dels + p
            s7 = s7 + (y(i) + y(im1))*dels/2.
            if (sigma /= 0.) then
               call snhcsh_1(ss, cs, sigmap*dels, 3)
               s7 = s7 + (yp(i) + yp(im1))*dels*(cs - ss/2.)/(sigmap*&
               &sigmap*(1.+ss))
            else
               s7 = s7 - (yp(i) + yp(im1))*dels*dels*dels/24.
            end if
         end do
      end if
      s8 = 0.
      if (.not. (ilm1 < ium1 .or. ideltp == 0 .and. bdy)) then
!
! integrate from xxl to xxu, store in s8
!
         delu1 = xxu - x(ium1)
         delu2 = xiu - xxu
         dell1 = xxl - x(ium1)
         dell2 = xiu - xxl
         dels = xiu - x(ium1)
         deli = xxu - xxl
         t1 = (delu1 + dell1)*deli/(2.*dels)
         t2 = (delu2 + dell2)*deli/(2.*dels)
         s8 = t1*y(iu) + t2*y(ium1)
         if (sigma /= 0.) then
            call snhcsh_1(dummy, cu1, sigmap*delu1, 2)
            call snhcsh_1(dummy, cu2, sigmap*delu2, 2)
            call snhcsh_1(dummy, cl1, sigmap*dell1, 2)
            call snhcsh_1(dummy, cl2, sigmap*dell2, 2)
            call snhcsh_1(ss, dummy, sigmap*dels, -1)
            s8 = s8 + (yp(iu)*(delu1*delu1*(cu1 - ss/2.) - dell1*dell1*(cl1 -&
            &ss/2.)) + yp(ium1)*(dell2*dell2*(cl2 - ss/2.) - delu2*delu2*(&
            &cu2 - ss/2.)))/(sigmap*sigmap*dels*(1.+ss))
         else
            s8 = s8 - t1*(delu2*(dels + delu1) + dell2*(dels + dell1))*&
            &yp(iu)/12.-t2*(dell1*(dels + dell2) + delu1*(dels +&
            &delu2))*yp(ium1)/12.
         end if
      end if
      so = s1 + s2 + s6 + s7
      si = s3 + s4 + s5 + s8
      if (.not. bdy) then
         curvpi = real(ideltp)*(so + si) + real(isign)*si
         return
      end if
      curvpi = real(ideltp)*(so + si) + real(isign)*so
      return
   end function curvpi
   subroutine kurv1(n, x, y, slp1, slpn, islpsw, xp, yp, temp, s,&
   &sigma, ierr)

!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer n, islpsw, ierr
      real(kind=rprec) slp1, slpn, sigma
      real(kind=rprec), dimension(n) :: x, y, xp, yp, temp, s
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: nm1, np1, i, im1, ibak
      real(kind=rprec) :: sigmap, slpp1x, slpp1y, dels1, dels2,&
      &c1, c2, c3, sx, sy,&
      &delt, slppnx, slppny, delsn, delsnm, dx1, dy1, diag1, sdiag1,&
      &dx2, dy2, diag2, sdiag2, diag, diagin
!-----------------------------------------------
!
!
!                                 coded by alan kaylor cline
!                           from fitpack -- january 26, 1987
!                        a curve and surface fitting package
!                      a product of pleasant valley software
!                  8603 altus cove, austin, texas 78759, usa
!
! this subroutine determines the parameters necessary to
! compute a spline under tension forming a curve in the
! plane and passing through a sequence of pairs (x(1),y(1)),
! ...,(x(n),y(n)). for actual computation of points on the
! curve it is necessary to call the subroutine kurv2.
!
! on input--
!
!   n is the number of points to be interpolated (n.ge.2).
!
!   x is an array containing the n x-coordinates of the
!   points.
!
!   y is an array containing the n y-coordinates of the
!   points. (adjacent x-y pairs must be distinct, i. e.
!   either x(i) .ne. x(i+1) or y(i) .ne. y(i+1), for
!   i = 1,...,n-1.)
!
!   slp1 and slpn contain the desired values for the angles
!   (in radians) of the slope at (x(1),y(1)) and (x(n),y(n))
!   respectively. the angles are measured counter-clock-
!   wise from the x-axis and the positive sense of the curve
!   is assumed to be that moving from point 1 to point n.
!   the user may omit values for either or both of these
!   parameters and signal this with islpsw.
!
!   islpsw contains a switch indicating which slope data
!   should be used and which should be estimated by this
!   subroutine,
!          = 0 if slp1 and slpn are to be used,
!          = 1 if slp1 is to be used but not slpn,
!          = 2 if slpn is to be used but not slp1,
!          = 3 if both slp1 and slpn are to be estimated
!              internally.
!
!   xp and yp are arrays of length at least n.
!
!   temp is an array of length at least n which is used
!   for scratch storage.
!
!   s is an array of length at least n.
!
! and
!
!   sigma contains the tension factor. this value indicates
!   the curviness desired. if abs(sigma) is nearly zero
!   (e.g. .001) the resulting curve is approximately a cubic
!   spline. if abs(sigma) is large (e. g. 50.) the resulting
!   curve is nearly a polygonal line. if sigma equals zero a
!   cubic spline results. a standard value for sigma is
!   approximately 1. in absolute value.
!
! on output--
!
!   xp and yp contain information about the curvature of the
!   curve at the given nodes.
!
!   s contains the polygonal arclengths of the curve.
!
!   ierr contains an error flag,
!        = 0 for normal return,
!        = 1 if n is less than 2,
!        = 2 if adjacent coordinate pairs coincide.
!
! and
!
!   n, x, y, slp1, slpn, islpsw, and sigma are unaltered.
!
! this subroutine references package modules ceez, terms_1,
! and snhcsh_1.
!
!-----------------------------------------------------------
!
      nm1 = n - 1
      np1 = n + 1
      ierr = 0
      if (n > 1) then
!
! determine polygonal arclengths
!
         s(1) = 0.
         do i = 2, n
            im1 = i - 1
            s(i) = s(im1) + sqrt((x(i) - x(im1))**2 + (y(i) - y(im1))**2)
         end do
!
! denormalize tension factor
!
         sigmap = abs(sigma)*real(n - 1)/s(n)
!
! approximate end slopes
!
         if (islpsw < 2) then
            slpp1x = cos(slp1)
            slpp1y = sin(slp1)
         else
            dels1 = s(2) - s(1)
            dels2 = dels1 + dels1
            if (n > 2) dels2 = s(3) - s(1)
            if (dels1 == 0. .or. dels2 == 0.) go to 12
            call ceez(dels1, dels2, sigmap, c1, c2, c3, n)
            sx = c1*x(1) + c2*x(2)
            sy = c1*y(1) + c2*y(2)
            if (n /= 2) then
               sx = sx + c3*x(3)
               sy = sy + c3*y(3)
            end if
            delt = sqrt(sx*sx + sy*sy)
            slpp1x = sx/delt
            slpp1y = sy/delt
         end if
         if (islpsw /= 1 .and. islpsw /= 3) then
            slppnx = cos(slpn)
            slppny = sin(slpn)
         else
            delsn = s(n) - s(nm1)
            delsnm = delsn + delsn
            if (n > 2) delsnm = s(n) - s(n - 2)
            if (delsn == 0. .or. delsnm == 0.) go to 12
            call ceez((-delsn), (-delsnm), sigmap, c1, c2, c3, n)
            sx = c1*x(n) + c2*x(nm1)
            sy = c1*y(n) + c2*y(nm1)
            if (n /= 2) then
               sx = sx + c3*x(n - 2)
               sy = sy + c3*y(n - 2)
            end if
            delt = sqrt(sx*sx + sy*sy)
            slppnx = sx/delt
            slppny = sy/delt
!
! set up right hand sides and tridiagonal system for xp and
! yp and perform forward elimination
!
         end if
         dx1 = (x(2) - x(1))/s(2)
         dy1 = (y(2) - y(1))/s(2)
         call terms_1(diag1, sdiag1, sigmap, s(2))
         xp(1) = (dx1 - slpp1x)/diag1
         yp(1) = (dy1 - slpp1y)/diag1
         temp(1) = sdiag1/diag1
         if (n /= 2) then
            do i = 2, nm1
               dels2 = s(i + 1) - s(i)
               if (dels2 == 0.) go to 12
               dx2 = (x(i + 1) - x(i))/dels2
               dy2 = (y(i + 1) - y(i))/dels2
               call terms_1(diag2, sdiag2, sigmap, dels2)
               diag = diag1 + diag2 - sdiag1*temp(i - 1)
               diagin = 1./diag
               xp(i) = (dx2 - dx1 - sdiag1*xp(i - 1))*diagin
               yp(i) = (dy2 - dy1 - sdiag1*yp(i - 1))*diagin
               temp(i) = sdiag2*diagin
               dx1 = dx2
               dy1 = dy2
               diag1 = diag2
               sdiag1 = sdiag2
            end do
         end if
         diag = diag1 - sdiag1*temp(nm1)
         xp(n) = (slppnx - dx1 - sdiag1*xp(nm1))/diag
         yp(n) = (slppny - dy1 - sdiag1*yp(nm1))/diag
!
! perform back substitution
!
         do i = 2, n
            ibak = np1 - i
            xp(ibak) = xp(ibak) - temp(ibak)*xp(ibak + 1)
            yp(ibak) = yp(ibak) - temp(ibak)*yp(ibak + 1)
         end do
         return
!
! too few points
!
      end if
      ierr = 1
      return
!
! coincident adjacent points
!
12    continue
      ierr = 2
      return
   end subroutine kurv1
   subroutine kurv2(t, xs, ys, n, x, y, xp, yp, s, sigma)

!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer n
      real(kind=rprec) t, xs, ys, sigma
      real(kind=rprec), dimension(n) :: x, y, xp, yp, s
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: im1, i
      real(kind=rprec) :: tn, sigmap, del1, del2, dels,&
      &sumx, sumy, d, c1, c2,&
      &sigdel, ss, dummy, s1, s2
!-----------------------------------------------
!   E x t e r n a l   F u n c t i o n s
!-----------------------------------------------
!-----------------------------------------------
!
!
!                                 coded by alan kaylor cline
!                           from fitpack -- january 26, 1987
!                        a curve and surface fitting package
!                      a product of pleasant valley software
!                  8603 altus cove, austin, texas 78759, usa
!
! this subroutine performs the mapping of points in the
! interval (0.,1.) onto a curve in the plane. the subroutine
! kurv1 should be called earlier to determine certain
! necessary parameters. the resulting curve has a parametric
! representation both of whose components are splines under
! tension and functions of the polygonal arclength
! parameter.
!
! on input--
!
!   t contains a real value to be mapped to a point on the
!   curve. the interval (0.,1.) is mapped onto the entire
!   curve, with 0. mapping to (x(1),y(1)) and 1. mapping
!   to (x(n),y(n)). values outside this interval result in
!   extrapolation.
!
!   n contains the number of points which were specified
!   to determine the curve.
!
!   x and y are arrays containing the x- and y-coordinates
!   of the specified points.
!
!   xp and yp are the arrays output from kurv1 containing
!   curvature information.
!
!   s is an array containing the polygonal arclengths of
!   the curve.
!
! and
!
!   sigma contains the tension factor (its sign is ignored).
!
! the parameters n, x, y, xp, yp, s, and sigma should be
! input unaltered from the output of kurv1.
!
! on output--
!
!   xs and ys contain the x- and y-coordinates of the image
!   point on the curve.
!
! none of the input parameters are altered.
!
! this subroutine references package modules intrvl_1 and
! snhcsh_1.
!
!-----------------------------------------------------------
!
! determine interval
!
      tn = s(n)*t
      im1 = intrvl_1(tn, s, n)
      i = im1 + 1
!
! denormalize tension factor
!
      sigmap = abs(sigma)*real(n - 1)/s(n)
!
! set up and perform interpolation
!
      del1 = tn - s(im1)
      del2 = s(i) - tn
      dels = s(i) - s(im1)
      sumx = (x(i)*del1 + x(im1)*del2)/dels
      sumy = (y(i)*del1 + y(im1)*del2)/dels
      if (sigmap == 0.) then
         d = del1*del2/(6.*dels)
         c1 = (del1 + dels)*d
         c2 = (del2 + dels)*d
         xs = sumx - xp(i)*c1 - xp(im1)*c2
         ys = sumy - yp(i)*c1 - yp(im1)*c2
         return
      end if
      sigdel = sigmap*dels
      call snhcsh_1(ss, dummy, sigdel, -1)
      call snhcsh_1(s1, dummy, sigmap*del1, -1)
      call snhcsh_1(s2, dummy, sigmap*del2, -1)
      d = sigdel*sigmap*(1.+ss)
      c1 = del1*(s1 - ss)/d
      c2 = del2*(s2 - ss)/d
      xs = sumx + xp(i)*c1 + xp(im1)*c2
      ys = sumy + yp(i)*c1 + yp(im1)*c2
      return
   end subroutine kurv2
   subroutine kurvd(t, xs, ys, xst, yst, xstt, ystt, n, x, y, xp, yp, s, sigma)

!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer n
      real(kind=rprec) t, xs, ys, xst, yst, xstt, ystt, sigma
      real(kind=rprec), dimension(n) :: x, y, xp, yp, s
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: im1, i
      real(kind=rprec) :: tn, sigmap, del1, del2, dels, sumx,&
      &sumy, sumxt, sumyt,&
      &dels6, d, c1, c2, ct1, ct2, ctt1, ctt2, sigdel, ss, dummy, s1,&
      &co1, s2, co2
!-----------------------------------------------
!   E x t e r n a l   F u n c t i o n s
!-----------------------------------------------
!-----------------------------------------------
!
!
!                                 coded by alan kaylor cline
!                           from fitpack -- january 26, 1987
!                        a curve and surface fitting package
!                      a product of pleasant valley software
!                  8603 altus cove, austin, texas 78759, usa
!
! this subroutine performs the mapping of points in the
! interval (0.,1.) onto a curve in the plane. it also
! returns the first and second derivatives of the component
! functions. the subroutine kurv1 should be called earlier
! to determine certain necessary parameters. the resulting
! curve has a parametric representation both of whose
! components are splines under tension and functions of the
! polygonal arclength parameter.
!
! on input--
!
!   t contains a real value to be mapped to a point on the
!   curve. the interval (0.,1.) is mapped onto the entire
!   curve, with 0. mapping to (x(1),y(1)) and 1. mapping
!   to (x(n),y(n)). values outside this interval result in
!   extrapolation.
!
!   n contains the number of points which were specified
!   to determine the curve.
!
!   x and y are arrays containing the x- and y-coordinates
!   of the specified points.
!
!   xp and yp are the arrays output from kurv1 containing
!   curvature information.
!
!   s is an array containing the polygonal arclengths of
!   the curve.
!
! and
!
!   sigma contains the tension factor (its sign is ignored).
!
! the parameters n, x, y, xp, yp, s, and sigma should be
! input unaltered from the output of kurv1.
!
! on output--
!
!   xs and ys contain the x- and y-coordinates of the image
!   point on the curve. xst and yst contain the first
!   derivatives of the x- and y-components of the mapping
!   with respect to t. xstt and ystt contain the second
!   derivatives of the x- and y-components of the mapping
!   with respect to t.
!
! none of the input parameters are altered.
!
! this subroutine references package modules intrvl_1 and
! snhcsh_1.
!
!-----------------------------------------------------------
!
! determine interval
!
      tn = s(n)*t
      im1 = intrvl_1(tn, s, n)
      i = im1 + 1
!
! denormalize tension factor
!
      sigmap = abs(sigma)*real(n - 1)/s(n)
!
! set up and perform interpolation
!
      del1 = tn - s(im1)
      del2 = s(i) - tn
      dels = s(i) - s(im1)
      sumx = (x(i)*del1 + x(im1)*del2)/dels
      sumy = (y(i)*del1 + y(im1)*del2)/dels
      sumxt = s(n)*(x(i) - x(im1))/dels
      sumyt = s(n)*(y(i) - y(im1))/dels
      if (sigmap == 0.) then
         dels6 = 6.*dels
         d = del1*del2/dels6
         c1 = -(del1 + dels)*d
         c2 = -(del2 + dels)*d
         dels6 = dels6/s(n)
         ct1 = (2.*del1*del1 - del2*(del1 + dels))/dels6
         ct2 = -(2.*del2*del2 - del1*(del2 + dels))/dels6
         dels = dels/(s(n)*s(n))
         ctt1 = del1/dels
         ctt2 = del2/dels
      else
         sigdel = sigmap*dels
         call snhcsh_1(ss, dummy, sigdel, -1)
         call snhcsh_1(s1, co1, sigmap*del1, 0)
         call snhcsh_1(s2, co2, sigmap*del2, 0)
         d = sigdel*sigmap*(1.+ss)
         c1 = del1*(s1 - ss)/d
         c2 = del2*(s2 - ss)/d
         ct1 = (co1 - ss)*s(n)/d
         ct2 = -(co2 - ss)*s(n)/d
         ctt1 = del1*(1.+s1)*s(n)*s(n)/(dels*(1.+ss))
         ctt2 = del2*(1.+s2)*s(n)*s(n)/(dels*(1.+ss))
      end if
      xs = sumx + c1*xp(i) + c2*xp(im1)
      ys = sumy + c1*yp(i) + c2*yp(im1)
      xst = sumxt + ct1*xp(i) + ct2*xp(im1)
      yst = sumyt + ct1*yp(i) + ct2*yp(im1)
      xstt = ctt1*xp(i) + ctt2*xp(im1)
      ystt = ctt1*yp(i) + ctt2*yp(im1)
      return
   end subroutine kurvd
   subroutine kurvp1(n, x, y, xp, yp, temp, s, sigma, ierr)

!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer n, ierr
      real(kind=rprec) sigma
      real(kind=rprec), dimension(n) :: x, y, xp, yp
      real(kind=rprec), dimension(1) :: temp
      real(kind=rprec), dimension(n) :: s
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: nm1, np1, i, im1, npi, ibak, npibak
      real(kind=rprec) :: sigmap, dels1, dx1, dy1, diag1,&
      &sdiag1, dels2, dx2, dy2,&
      &diag2, sdiag2, diag, diagin, xpn, ypn
!-----------------------------------------------
!
!
!                                 coded by alan kaylor cline
!                           from fitpack -- january 26, 1987
!                        a curve and surface fitting package
!                      a product of pleasant valley software
!                  8603 altus cove, austin, texas 78759, usa
!
! this subroutine determines the parameters necessary to
! compute a spline under tension forming a closed curve in
! the plane and passing through a sequence of pairs
! (x(1),y(1)),...,(x(n),y(n)). for actual computation of
! points on the curve it is necessary to call the subroutine
! kurvp2.
!
! on input--
!
!   n is the number of points to be interpolated (n.ge.2).
!
!   x is an array containing the n x-coordinates of the
!   points.
!
!   y is an array containing the n y-coordinates of the
!   points. (adjacent x-y pairs must be distinct, i. e.
!   either x(i) .ne. x(i+1) or y(i) .ne. y(i+1), for
!   i = 1,...,n-1 and either x(1) .ne. x(n) or y(1) .ne. y(n).)
!
!   xp and yp are arrays of length at least n.
!
!   temp is an array of length at least 2*n which is used
!   for scratch storage.
!
!   s is an array of length at least n.
!
! and
!
!   sigma contains the tension factor. this value indicates
!   the curviness desired. if abs(sigma) is nearly zero
!   (e.g. .001) the resulting curve is approximately a cubic
!   spline. if abs(sigma) is large (e. g. 50.) the resulting
!   curve is nearly a polygonal line. if sigma equals zero a
!   cubic spline results. a standard value for sigma is
!   approximately 1. in absolute value.
!
! on output--
!
!   xp and yp contain information about the curvature of the
!   curve at the given nodes.
!
!   s contains the polygonal arclengths of the curve.
!
!   ierr contains an error flag,
!        = 0 for normal return,
!        = 1 if n is less than 2,
!        = 2 if adjacent coordinate pairs coincide.
!
! and
!
!   n, x, y, and sigma are unaltered,
!
! this subroutine references package modules terms_1 and
! snhcsh_1.
!
!-----------------------------------------------------------
!
      nm1 = n - 1
      np1 = n + 1
      ierr = 0
      if (n > 1) then
!
! determine polygonal arclengths
!
         s(1) = sqrt((x(n) - x(1))**2 + (y(n) - y(1))**2)
         do i = 2, n
            im1 = i - 1
            s(i) = s(im1) + sqrt((x(i) - x(im1))**2 + (y(i) - y(im1))**2)
         end do
!
! denormalize tension factor
!
         sigmap = abs(sigma)*real(n)/s(n)
!
! set up right hand sides of tridiagonal (with corner
! elements) linear system for xp and yp
!
         dels1 = s(1)
         if (dels1 == 0.) go to 8
         dx1 = (x(1) - x(n))/dels1
         dy1 = (y(1) - y(n))/dels1
         call terms_1(diag1, sdiag1, sigmap, dels1)
         dels2 = s(2) - s(1)
         if (dels2 == 0.) go to 8
         dx2 = (x(2) - x(1))/dels2
         dy2 = (y(2) - y(1))/dels2
         call terms_1(diag2, sdiag2, sigmap, dels2)
         diag = diag1 + diag2
         diagin = 1./diag
         xp(1) = (dx2 - dx1)*diagin
         yp(1) = (dy2 - dy1)*diagin
         temp(np1) = -sdiag1*diagin
         temp(1) = sdiag2*diagin
         dx1 = dx2
         dy1 = dy2
         diag1 = diag2
         sdiag1 = sdiag2
         if (n /= 2) then
            do i = 2, nm1
               npi = n + i
               dels2 = s(i + 1) - s(i)
               if (dels2 == 0.) go to 8
               dx2 = (x(i + 1) - x(i))/dels2
               dy2 = (y(i + 1) - y(i))/dels2
               call terms_1(diag2, sdiag2, sigmap, dels2)
               diag = diag1 + diag2 - sdiag1*temp(i - 1)
               diagin = 1./diag
               xp(i) = (dx2 - dx1 - sdiag1*xp(i - 1))*diagin
               yp(i) = (dy2 - dy1 - sdiag1*yp(i - 1))*diagin
               temp(npi) = -temp(npi - 1)*sdiag1*diagin
               temp(i) = sdiag2*diagin
               dx1 = dx2
               dy1 = dy2
               diag1 = diag2
               sdiag1 = sdiag2
            end do
         end if
         dels2 = s(1)
         dx2 = (x(1) - x(n))/dels2
         dy2 = (y(1) - y(n))/dels2
         call terms_1(diag2, sdiag2, sigmap, dels2)
         xp(n) = dx2 - dx1
         yp(n) = dy2 - dy1
         temp(nm1) = temp(2*n - 1) - temp(nm1)
         if (n /= 2) then
!
! perform first step of back substitution
!
            do i = 3, n
               ibak = np1 - i
               npibak = n + ibak
               xp(ibak) = xp(ibak) - temp(ibak)*xp(ibak + 1)
               yp(ibak) = yp(ibak) - temp(ibak)*yp(ibak + 1)
               temp(ibak) = temp(npibak) - temp(ibak)*temp(ibak + 1)
            end do
         end if
         xp(n) = (xp(n) - sdiag2*xp(1) - sdiag1*xp(nm1))/(diag1 + diag2 +&
         &sdiag2*temp(1) + sdiag1*temp(nm1))
         yp(n) = (yp(n) - sdiag2*yp(1) - sdiag1*yp(nm1))/(diag1 + diag2 +&
         &sdiag2*temp(1) + sdiag1*temp(nm1))
!
! perform second step of back substitution
!
         xpn = xp(n)
         ypn = yp(n)
         xp(:nm1) = xp(:nm1) + temp(:nm1)*xpn
         yp(:nm1) = yp(:nm1) + temp(:nm1)*ypn
         return
!
! too few points
!
      end if
      ierr = 1
      return
!
! coincident adjacent points
!
8     continue
      ierr = 2
      return
   end subroutine kurvp1
   subroutine kurvp2(t, xs, ys, n, x, y, xp, yp, s, sigma)

!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer n
      real(kind=rprec) t, xs, ys, sigma
      real(kind=rprec), dimension(n) :: x, y, xp, yp, s
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: im1, i
      real(kind=rprec) :: tn, sigmap, si, del1, del2, dels,&
      &sumx, sumy, d, c1, c2,&
      &sigdel, ss, dummy, s1, s2, ci, cim1
!-----------------------------------------------
!   E x t e r n a l   F u n c t i o n s
!-----------------------------------------------
!-----------------------------------------------
!
!
!                                 coded by alan kaylor cline
!                           from fitpack -- january 26, 1987
!                        a curve and surface fitting package
!                      a product of pleasant valley software
!                  8603 altus cove, austin, texas 78759, usa
!
! this subroutine performs the mapping of points in the
! interval (0.,1.) onto a closed curve in the plane. the
! subroutine kurvp1 should be called earlier to determine
! certain necessary parameters. the resulting curve has a
! parametric representation both of whose components are
! periodic splines under tension and functions of the poly-
! gonal arclength parameter.
!
! on input--
!
!   t contains a value to be mapped onto the curve. the
!   interval (0.,1.) is mapped onto the entire closed curve
!   with both 0. and 1. mapping to (x(1),y(1)). the mapping
!   is periodic with period one thus any interval of the
!   form (tt,tt+1.) maps onto the entire curve.
!
!   n contains the number of points which were specified
!   to determine the curve.
!
!   x and y are arrays containing the x- and y-coordinates
!   of the specified points.
!
!   xp and yp are the arrays output from kurvp1 containing
!   curvature information.
!
!   s is an array containing the polygonal arclengths of
!   the curve.
!
! and
!
!   sigma contains the tension factor (its sign is ignored).
!
! the parameters n, x, y, xp, yp, s and sigma should
! be input unaltered from the output of kurvp1.
!
! on output--
!
!   xs and ys contain the x- and y-coordinates of the image
!   point on the curve.
!
! none of the input parameters are altered.
!
! this subroutine references package modules intrvl_1 and
! snhcsh_1.
!
!-----------------------------------------------------------
!
! determine interval
!
      tn = t - real(int(t))
      if (tn < 0.) tn = tn + 1.
      tn = s(n)*tn + s(1)
      im1 = n
      if (tn < s(n)) im1 = intrvl_1(tn, s, n)
      i = im1 + 1
      if (i > n) i = 1
!
! denormalize tension factor
!
      sigmap = abs(sigma)*real(n)/s(n)
!
! set up and perform interpolation
!
      si = s(i)
      if (im1 == n) si = s(n) + s(1)
      del1 = tn - s(im1)
      del2 = si - tn
      dels = si - s(im1)
      sumx = (x(i)*del1 + x(im1)*del2)/dels
      sumy = (y(i)*del1 + y(im1)*del2)/dels
      if (sigmap == 0.) then
         d = del1*del2/(6.*dels)
         c1 = (del1 + dels)*d
         c2 = (del2 + dels)*d
         xs = sumx - xp(i)*c1 - xp(im1)*c2
         ys = sumy - yp(i)*c1 - yp(im1)*c2
         return
      end if
      sigdel = sigmap*dels
      call snhcsh_1(ss, dummy, sigdel, -1)
      call snhcsh_1(s1, dummy, sigmap*del1, -1)
      call snhcsh_1(s2, dummy, sigmap*del2, -1)
      d = sigdel*sigmap*(1.+ss)
      ci = del1*(s1 - ss)/d
      cim1 = del2*(s2 - ss)/d
      xs = sumx + xp(i)*ci + xp(im1)*cim1
      ys = sumy + yp(i)*ci + yp(im1)*cim1
      return
   end subroutine kurvp2
   subroutine kurvpd(t, xs, ys, xst, yst, xstt, ystt, n, x, y, xp, yp, s, sigma)

!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer n
      real(kind=rprec) t, xs, ys, xst, yst, xstt, ystt, sigma
      real(kind=rprec), dimension(n) :: x, y, xp, yp, s
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: im1, i
      real(kind=rprec) :: tn, sigmap, si, del1, del2, dels, sumx,&
      &sumy, sumxt, sumyt&
      &, dels6, d, c1, c2, ct1, ct2, ctt1, ctt2, sigdel, ss, dummy, s1&
      &, co1, s2, co2
!-----------------------------------------------
!   E x t e r n a l   F u n c t i o n s
!-----------------------------------------------
!-----------------------------------------------
!
!
!                                 coded by alan kaylor cline
!                           from fitpack -- january 26, 1987
!                        a curve and surface fitting package
!                      a product of pleasant valley software
!                  8603 altus cove, austin, texas 78759, usa
!
! this subroutine performs the mapping of points in the
! interval (0.,1.) onto a closed curve in the plane. it also
! returns the first and second derivatives of the component
! functions. the subroutine kurvp1 should be called earlier
! to determine certain necessary parameters. the resulting
! curve has a parametric representation both of whose
! components are periodic splines under tension and
! functions of the polygonal arclength parameter.
!
! on input--
!
!   t contains a value to be mapped onto the curve. the
!   interval (0.,1.) is mapped onto the entire closed curve
!   with both 0. and 1. mapping to (x(1),y(1)). the mapping
!   is periodic with period one thus any interval of the
!   form (tt,tt+1.) maps onto the entire curve.
!
!   n contains the number of points which were specified
!   to determine the curve.
!
!   x and y are arrays containing the x- and y-coordinates
!   of the specified points.
!
!   xp and yp are the arrays output from kurvp1 containing
!   curvature information.
!
!   s is an array containing the polygonal arclengths of
!   the curve.
!
! and
!
!   sigma contains the tension factor (its sign is ignored).
!
! the parameters n, x, y, xp, yp, s and sigma should
! be input unaltered from the output of kurvp1.
!
! on output--
!
!   xs and ys contain the x- and y-coordinates of the image
!   point on the curve. xst and yst contain the first
!   derivatives of the x- and y-components of the mapping
!   with respect to t. xstt and ystt contain the second
!   derivatives of the x- and y-components of the mapping
!   with respect to t.
!
! none of the input parameters are altered.
!
! this subroutine references package modules intrvl_1 and
! snhcsh_1.
!
!-----------------------------------------------------------
!
! determine interval
!
      tn = t - real(int(t))
      if (tn < 0.) tn = tn + 1.
      tn = s(n)*tn + s(1)
      im1 = n
      if (tn < s(n)) im1 = intrvl_1(tn, s, n)
      i = im1 + 1
      if (i > n) i = 1
!
! denormalize tension factor
!
      sigmap = abs(sigma)*real(n)/s(n)
!
! set up and perform interpolation
!
      si = s(i)
      if (im1 == n) si = s(n) + s(1)
      del1 = tn - s(im1)
      del2 = si - tn
      dels = si - s(im1)
      sumx = (x(i)*del1 + x(im1)*del2)/dels
      sumy = (y(i)*del1 + y(im1)*del2)/dels
      sumxt = s(n)*(x(i) - x(im1))/dels
      sumyt = s(n)*(y(i) - y(im1))/dels
      if (sigmap == 0.) then
         dels6 = 6.*dels
         d = del1*del2/dels6
         c1 = -(del1 + dels)*d
         c2 = -(del2 + dels)*d
         dels6 = dels6/s(n)
         ct1 = (2.*del1*del1 - del2*(del1 + dels))/dels6
         ct2 = -(2.*del2*del2 - del1*(del2 + dels))/dels6
         dels = dels/(s(n)*s(n))
         ctt1 = del1/dels
         ctt2 = del2/dels
      else
         sigdel = sigmap*dels
         call snhcsh_1(ss, dummy, sigdel, -1)
         call snhcsh_1(s1, co1, sigmap*del1, 0)
         call snhcsh_1(s2, co2, sigmap*del2, 0)
         d = sigdel*sigmap*(1.+ss)
         c1 = del1*(s1 - ss)/d
         c2 = del2*(s2 - ss)/d
         ct1 = (co1 - ss)*s(n)/d
         ct2 = -(co2 - ss)*s(n)/d
         ctt1 = del1*(1.+s1)*s(n)*s(n)/(dels*(1.+ss))
         ctt2 = del2*(1.+s2)*s(n)*s(n)/(dels*(1.+ss))
      end if
      xs = sumx + c1*xp(i) + c2*xp(im1)
      ys = sumy + c1*yp(i) + c2*yp(im1)
      xst = sumxt + ct1*xp(i) + ct2*xp(im1)
      yst = sumyt + ct1*yp(i) + ct2*yp(im1)
      xstt = ctt1*xp(i) + ctt2*xp(im1)
      ystt = ctt1*yp(i) + ctt2*yp(im1)
      return
   end subroutine kurvpd
   subroutine surf1(m, n, x, y, z, iz, zx1, zxm, zy1, zyn, zxy11,&
   &zxym1, zxy1n, zxymn, islpsw, zp, temp, sigma, ierr)

!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer m, n, iz, islpsw, ierr
      real(kind=rprec) zxy11, zxym1, zxy1n, zxymn, sigma
      real(kind=rprec), dimension(m) :: x
      real(kind=rprec), dimension(n) :: y
      real(kind=rprec), dimension(iz, n) :: z
      real(kind=rprec), dimension(n) :: zx1, zxm
      real(kind=rprec), dimension(m) :: zy1, zyn
      real(kind=rprec), dimension(m, n, 3) :: zp
      real(kind=rprec), dimension(1) :: temp
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: mm1, mp1, nm1, np1, npm, i, npi, j, npmpj, jm1, jp1,&
      &jbak, jbakp1, im1, ip1, ibak, ibakp1, npibak
      real(kind=rprec) :: sigmay, dely1, dely2, c1, c2, c3,&
      &delyn, delynm, sigmax,&
      &delx1, delx2, zxy1ns, delxm, delxmm, zxymns, del1, deli, diag1&
      &, sdiag1, diagi, del2, diag2, sdiag2, diagin, t
!-----------------------------------------------
!
!
!                                 coded by alan kaylor cline
!                           from fitpack -- january 26, 1987
!                        a curve and surface fitting package
!                      a product of pleasant valley software
!                  8603 altus cove, austin, texas 78759, usa
!
! this subroutine determines the parameters necessary to
! compute an interpolatory surface passing through a rect-
! angular grid of functional values. the surface determined
! can be represented as the tensor product of splines under
! tension. the x- and y-partial derivatives around the
! boundary and the x-y-partial derivatives at the four
! corners may be specified or omitted. for actual mapping
! of points onto the surface it is necessary to call the
! function surf2.
!
! on input--
!
!   m is the number of grid lines in the x-direction, i. e.
!   lines parallel to the y-axis (m .ge. 2).
!
!   n is the number of grid lines in the y-direction, i. e.
!   lines parallel to the x-axis (n .ge. 2).
!
!   x is an array of the m x-coordinates of the grid lines
!   in the x-direction. these should be strictly increasing.
!
!   y is an array of the n y-coordinates of the grid lines
!   in the y-direction. these should be strictly increasing.
!
!   z is an array of the m * n functional values at the grid
!   points, i. e. z(i,j) contains the functional value at
!   (x(i),y(j)) for i = 1,...,m and j = 1,...,n.
!
!   iz is the row dimension of the matrix z used in the
!   calling program (iz .ge. m).
!
!   zx1 and zxm are arrays of the m x-partial derivatives
!   of the function along the x(1) and x(m) grid lines,
!   respectively. thus zx1(j) and zxm(j) contain the x-part-
!   ial derivatives at the points (x(1),y(j)) and
!   (x(m),y(j)), respectively, for j = 1,...,n. either of
!   these parameters will be ignored (and approximations
!   supplied internally) if islpsw so indicates.
!
!   zy1 and zyn are arrays of the n y-partial derivatives
!   of the function along the y(1) and y(n) grid lines,
!   respectively. thus zy1(i) and zyn(i) contain the y-part-
!   ial derivatives at the points (x(i),y(1)) and
!   (x(i),y(n)), respectively, for i = 1,...,m. either of
!   these parameters will be ignored (and estimations
!   supplied internally) if islpsw so indicates.
!
!   zxy11, zxym1, zxy1n, and zxymn are the x-y-partial
!   derivatives of the function at the four corners,
!   (x(1),y(1)), (x(m),y(1)), (x(1),y(n)), and (x(m),y(n)),
!   respectively. any of the parameters will be ignored (and
!   estimations supplied internally) if islpsw so indicates.
!
!   islpsw contains a switch indicating which boundary
!   derivative information is user-supplied and which
!   should be estimated by this subroutine. to determine
!   islpsw, let
!        i1 = 0 if zx1 is user-supplied (and = 1 otherwise),
!        i2 = 0 if zxm is user-supplied (and = 1 otherwise),
!        i3 = 0 if zy1 is user-supplied (and = 1 otherwise),
!        i4 = 0 if zyn is user-supplied (and = 1 otherwise),
!        i5 = 0 if zxy11 is user-supplied
!                                       (and = 1 otherwise),
!        i6 = 0 if zxym1 is user-supplied
!                                       (and = 1 otherwise),
!        i7 = 0 if zxy1n is user-supplied
!                                       (and = 1 otherwise),
!        i8 = 0 if zxymn is user-supplied
!                                       (and = 1 otherwise),
!   then islpsw = i1 + 2*i2 + 4*i3 + 8*i4 + 16*i5 + 32*i6
!                   + 64*i7 + 128*i8
!   thus islpsw = 0 indicates all derivative information is
!   user-supplied and islpsw = 255 indicates no derivative
!   information is user-supplied. any value between these
!   limits is valid.
!
!   zp is an array of at least 3*m*n locations.
!
!   temp is an array of at least n+n+m locations which is
!   used for scratch storage.
!
! and
!
!   sigma contains the tension factor. this value indicates
!   the curviness desired. if abs(sigma) is nearly zero
!   (e. g. .001) the resulting surface is approximately the
!   tensor product of cubic splines. if abs(sigma) is large
!   (e. g. 50.) the resulting surface is approximately
!   bi-linear. if sigma equals zero tensor products of
!   cubic splines result. a standard value for sigma is
!   approximately 1. in absolute value.
!
! on output--
!
!   zp contains the values of the xx-, yy-, and xxyy-partial
!   derivatives of the surface at the given nodes.
!
!   ierr contains an error flag,
!        = 0 for normal return,
!        = 1 if n is less than 2 or m is less than 2,
!        = 2 if the x-values or y-values are not strictly
!            increasing.
!
! and
!
!   m, n, x, y, z, iz, zx1, zxm, zy1, zyn, zxy11, zxym1,
!   zxy1n, zxymn, islpsw, and sigma are unaltered.
!
! this subroutine references package modules ceez, terms_1,
! and snhcsh_1.
!
!-----------------------------------------------------------
!
      mm1 = m - 1
      mp1 = m + 1
      nm1 = n - 1
      np1 = n + 1
      npm = n + m
      ierr = 0
      if (n > 1 .and. m > 1) then
         if (y(n) <= y(1)) go to 47
!
! denormalize tension factor in y-direction
!
         sigmay = abs(sigma)*real(n - 1)/(y(n) - y(1))
!
! obtain y-partial derivatives along y = y(1)
!
         if ((islpsw/8)*2 == islpsw/4) then
            zp(:m, 1, 1) = zy1(:m)
         else
            dely1 = y(2) - y(1)
            dely2 = dely1 + dely1
            if (n > 2) dely2 = y(3) - y(1)
            if (dely1 <= 0. .or. dely2 <= dely1) go to 47
            call ceez(dely1, dely2, sigmay, c1, c2, c3, n)
            zp(:m, 1, 1) = c1*z(:m, 1) + c2*z(:m, 2)
            if (n /= 2) then
               zp(:m, 1, 1) = zp(:m, 1, 1) + c3*z(:m, 3)
!
! obtain y-partial derivatives along y = y(n)
!
            end if
         end if
         if ((islpsw/16)*2 == islpsw/8) then
            do i = 1, m
               npi = n + i
               temp(npi) = zyn(i)
            end do
         else
            delyn = y(n) - y(nm1)
            delynm = delyn + delyn
            if (n > 2) delynm = y(n) - y(n - 2)
            if (delyn <= 0. .or. delynm <= delyn) go to 47
            call ceez((-delyn), (-delynm), sigmay, c1, c2, c3, n)
            do i = 1, m
               npi = n + i
               temp(npi) = c1*z(i, n) + c2*z(i, nm1)
            end do
            if (n /= 2) then
               do i = 1, m
                  npi = n + i
                  temp(npi) = temp(npi) + c3*z(i, n - 2)
               end do
            end if
         end if
         if (x(m) <= x(1)) go to 47
!
! denormalize tension factor in x-direction
!
         sigmax = abs(sigma)*real(m - 1)/(x(m) - x(1))
!
! obtain x-partial derivatives along x = x(1)
!
         if ((islpsw/2)*2 /= islpsw) go to 12
         zp(1, :n, 2) = zx1(:n)
         if ((islpsw/32)*2 == islpsw/16 .and. (islpsw/128)*2 == islpsw/64)&
         &go to 15
12       continue
         delx1 = x(2) - x(1)
         delx2 = delx1 + delx1
         if (m > 2) delx2 = x(3) - x(1)
         if (delx1 <= 0. .or. delx2 <= delx1) go to 47
         call ceez(delx1, delx2, sigmax, c1, c2, c3, m)
         if ((islpsw/2)*2 /= islpsw) then
            zp(1, :n, 2) = c1*z(1, :n) + c2*z(2, :n)
            if (m /= 2) then
               zp(1, :n, 2) = zp(1, :n, 2) + c3*z(3, :n)
!
! obtain x-y-partial derivative at (x(1),y(1))
!
            end if
         end if
15       continue
         if ((islpsw/32)*2 == islpsw/16) then
            zp(1, 1, 3) = zxy11
         else
            zp(1, 1, 3) = c1*zp(1, 1, 1) + c2*zp(2, 1, 1)
            if (m > 2) zp(1, 1, 3) = zp(1, 1, 3) + c3*zp(3, 1, 1)
!
! obtain x-y-partial derivative at (x(1),y(n))
!
         end if
         if ((islpsw/128)*2 == islpsw/64) then
            zxy1ns = zxy1n
         else
            zxy1ns = c1*temp(n + 1) + c2*temp(n + 2)
            if (m > 2) zxy1ns = zxy1ns + c3*temp(n + 3)
!
! obtain x-partial derivative along x = x(m)
!
         end if
         if ((islpsw/4)*2 /= islpsw/2) go to 21
         do j = 1, n
            npmpj = npm + j
            temp(npmpj) = zxm(j)
         end do
         if ((islpsw/64)*2 == islpsw/32 .and. (islpsw/256)*2 == islpsw/128)&
         &go to 24
21       continue
         delxm = x(m) - x(mm1)
         delxmm = delxm + delxm
         if (m > 2) delxmm = x(m) - x(m - 2)
         if (delxm <= 0. .or. delxmm <= delxm) go to 47
         call ceez((-delxm), (-delxmm), sigmax, c1, c2, c3, m)
         if ((islpsw/4)*2 /= islpsw/2) then
            do j = 1, n
               npmpj = npm + j
               temp(npmpj) = c1*z(m, j) + c2*z(mm1, j)
            end do
            if (m /= 2) then
               do j = 1, n
                  npmpj = npm + j
                  temp(npmpj) = temp(npmpj) + c3*z(m - 2, j)
               end do
!
! obtain x-y-partial derivative at (x(m),y(1))
!
            end if
         end if
24       continue
         if ((islpsw/64)*2 == islpsw/32) then
            zp(m, 1, 3) = zxym1
         else
            zp(m, 1, 3) = c1*zp(m, 1, 1) + c2*zp(mm1, 1, 1)
            if (m > 2) zp(m, 1, 3) = zp(m, 1, 3) + c3*zp(m - 2, 1, 1)
!
! obtain x-y-partial derivative at (x(m),y(n))
!
         end if
         if ((islpsw/256)*2 == islpsw/128) then
            zxymns = zxymn
         else
            zxymns = c1*temp(npm) + c2*temp(npm - 1)
            if (m > 2) zxymns = zxymns + c3*temp(npm - 2)
!
! set up right hand sides and tridiagonal system for y-grid
! perform forward elimination
!
         end if
         del1 = y(2) - y(1)
         if (del1 <= 0.) go to 47
         deli = 1./del1
         zp(:m, 2, 1) = deli*(z(:m, 2) - z(:m, 1))
         zp(1, 2, 3) = deli*(zp(1, 2, 2) - zp(1, 1, 2))
         zp(m, 2, 3) = deli*(temp(npm + 2) - temp(npm + 1))
         call terms_1(diag1, sdiag1, sigmay, del1)
         diagi = 1./diag1
         zp(:m, 1, 1) = diagi*(zp(:m, 2, 1) - zp(:m, 1, 1))
         zp(1, 1, 3) = diagi*(zp(1, 2, 3) - zp(1, 1, 3))
         zp(m, 1, 3) = diagi*(zp(m, 2, 3) - zp(m, 1, 3))
         temp(1) = diagi*sdiag1
         if (n /= 2) then
            do j = 2, nm1
               jm1 = j - 1
               jp1 = j + 1
               npmpj = npm + j
               del2 = y(jp1) - y(j)
               if (del2 <= 0.) go to 47
               deli = 1./del2
               zp(:m, jp1, 1) = deli*(z(:m, jp1) - z(:m, j))
               zp(1, jp1, 3) = deli*(zp(1, jp1, 2) - zp(1, j, 2))
               zp(m, jp1, 3) = deli*(temp(npmpj + 1) - temp(npmpj))
               call terms_1(diag2, sdiag2, sigmay, del2)
               diagin = 1./(diag1 + diag2 - sdiag1*temp(jm1))
               zp(:m, j, 1) = diagin*(zp(:m, jp1, 1) - zp(:m, j, 1) - sdiag1*zp(:m&
               &, jm1, 1))
               zp(1, j, 3) = diagin*(zp(1, jp1, 3) - zp(1, j, 3) - sdiag1*zp(1, jm1&
               &, 3))
               zp(m, j, 3) = diagin*(zp(m, jp1, 3) - zp(m, j, 3) - sdiag1*zp(m, jm1&
               &, 3))
               temp(j) = diagin*sdiag2
               diag1 = diag2
               sdiag1 = sdiag2
            end do
         end if
         diagin = 1./(diag1 - sdiag1*temp(nm1))
         do i = 1, m
            npi = n + i
            zp(i, n, 1) = diagin*(temp(npi) - zp(i, n, 1) - sdiag1*zp(i, nm1, 1))
         end do
         zp(1, n, 3) = diagin*(zxy1ns - zp(1, n, 3) - sdiag1*zp(1, nm1, 3))
         temp(n) = diagin*(zxymns - zp(m, n, 3) - sdiag1*zp(m, nm1, 3))
!
! perform back substitution
!
         do j = 2, n
            jbak = np1 - j
            jbakp1 = jbak + 1
            t = temp(jbak)
            zp(:m, jbak, 1) = zp(:m, jbak, 1) - t*zp(:m, jbakp1, 1)
            zp(1, jbak, 3) = zp(1, jbak, 3) - t*zp(1, jbakp1, 3)
            temp(jbak) = zp(m, jbak, 3) - t*temp(jbakp1)
         end do
!
! set up right hand sides and tridiagonal system for x-grid
! perform forward elimination
!
         del1 = x(2) - x(1)
         if (del1 <= 0.) go to 47
         deli = 1./del1
         zp(2, :n, 2) = deli*(z(2, :n) - z(1, :n))
         zp(2, :n, 3) = deli*(zp(2, :n, 1) - zp(1, :n, 1))
         call terms_1(diag1, sdiag1, sigmax, del1)
         diagi = 1./diag1
         zp(1, :n, 2) = diagi*(zp(2, :n, 2) - zp(1, :n, 2))
         zp(1, :n, 3) = diagi*(zp(2, :n, 3) - zp(1, :n, 3))
         temp(n + 1) = diagi*sdiag1
         if (m /= 2) then
            do i = 2, mm1
               im1 = i - 1
               ip1 = i + 1
               npi = n + i
               del2 = x(ip1) - x(i)
               if (del2 <= 0.) go to 47
               deli = 1./del2
               zp(ip1, :n, 2) = deli*(z(ip1, :n) - z(i, :n))
               zp(ip1, :n, 3) = deli*(zp(ip1, :n, 1) - zp(i, :n, 1))
               call terms_1(diag2, sdiag2, sigmax, del2)
               diagin = 1./(diag1 + diag2 - sdiag1*temp(npi - 1))
               zp(i, :n, 2) = diagin*(zp(ip1, :n, 2) - zp(i, :n, 2) - sdiag1*zp(&
               &im1, :n, 2))
               zp(i, :n, 3) = diagin*(zp(ip1, :n, 3) - zp(i, :n, 3) - sdiag1*zp(&
               &im1, :n, 3))
               temp(npi) = diagin*sdiag2
               diag1 = diag2
               sdiag1 = sdiag2
            end do
         end if
         diagin = 1./(diag1 - sdiag1*temp(npm - 1))
         do j = 1, n
            npmpj = npm + j
            zp(m, j, 2) = diagin*(temp(npmpj) - zp(m, j, 2) - sdiag1*zp(mm1, j, 2))
            zp(m, j, 3) = diagin*(temp(j) - zp(m, j, 3) - sdiag1*zp(mm1, j, 3))
         end do
!
! perform back substitution
!
         do i = 2, m
            ibak = mp1 - i
            ibakp1 = ibak + 1
            npibak = n + ibak
            t = temp(npibak)
            zp(ibak, :n, 2) = zp(ibak, :n, 2) - t*zp(ibakp1, :n, 2)
            zp(ibak, :n, 3) = zp(ibak, :n, 3) - t*zp(ibakp1, :n, 3)
         end do
         return
!
! too few points
!
      end if
      ierr = 1
      return
!
! points not strictly increasing
!
47    continue
      ierr = 2
      return
   end subroutine surf1
   function surf2(xx, yy, m, n, x, y, z, iz, zp, sigma)

!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer m, n, iz
      real(kind=rprec) xx, yy, sigma, surf2
      real(kind=rprec), dimension(m) :: x
      real(kind=rprec), dimension(n) :: y
      real(kind=rprec), dimension(iz, n) :: z
      real(kind=rprec), dimension(m, n, 3) :: zp
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: jm1, j, im1, i
      real(kind=rprec) :: del1, del2, dels, f1, f2, fp1, fp2,&
      &sigmap, sinhm1, sinhms&
      &, sinhm2, sigmax, sigmay, zim1, zi, zxxim1, zxxi, dummy
!-----------------------------------------------
!   E x t e r n a l   F u n c t i o n s
!-----------------------------------------------
!-----------------------------------------------
!
!
! inline one dimensional spline under tension interpolation
!
!
! denormalize tension factor in x and y direction
!
      sigmax = abs(sigma)*real(m - 1)/(x(m) - x(1))
      sigmay = abs(sigma)*real(n - 1)/(y(n) - y(1))
!
! determine y interval
!
      jm1 = intrvl_1(yy, y, n)
      j = jm1 + 1
!
! determine x interval
!
      im1 = intrvl_1(xx, x, m)
      i = im1 + 1
      del1 = yy - y(jm1)
      del2 = y(j) - yy
      dels = y(j) - y(jm1)
      if (sigmay == 0.) then
!
! perform four interpolations in y-direction
!
         zim1 = hermz(z(i - 1, j - 1), z(i - 1, j), zp(i - 1, j - 1, 1), zp(i - 1, j, 1))
         zi = hermz(z(i, j - 1), z(i, j), zp(i, j - 1, 1), zp(i, j, 1))
         zxxim1 = hermz(zp(i - 1, j - 1, 2), zp(i - 1, j, 2), zp(i - 1, j - 1, 3), zp(i - 1, j&
         &, 3))
         zxxi = hermz(zp(i, j - 1, 2), zp(i, j, 2), zp(i, j - 1, 3), zp(i, j, 3))
      else
         call snhcsh_1(sinhm1, dummy, sigmay*del1, -1)
         call snhcsh_1(sinhm2, dummy, sigmay*del2, -1)
         call snhcsh_1(sinhms, dummy, sigmay*dels, -1)
         zim1 = hermnz(z(i - 1, j - 1), z(i - 1, j), zp(i - 1, j - 1, 1), zp(i - 1, j, 1),&
         &sigmay)
         zi = hermnz(z(i, j - 1), z(i, j), zp(i, j - 1, 1), zp(i, j, 1), sigmay)
         zxxim1 = hermnz(zp(i - 1, j - 1, 2), zp(i - 1, j, 2), zp(i - 1, j - 1, 3), zp(i - 1,&
         &j, 3), sigmay)
         zxxi = hermnz(zp(i, j - 1, 2), zp(i, j, 2), zp(i, j - 1, 3), zp(i, j, 3),&
         &sigmay)
!
! perform final interpolation in x-direction
!
      end if
      del1 = xx - x(im1)
      del2 = x(i) - xx
      dels = x(i) - x(im1)
      if (sigmax == 0.) then
         surf2 = hermz(zim1, zi, zxxim1, zxxi)
         return
      end if
      call snhcsh_1(sinhm1, dummy, sigmax*del1, -1)
      call snhcsh_1(sinhm2, dummy, sigmax*del2, -1)
      call snhcsh_1(sinhms, dummy, sigmax*dels, -1)
      surf2 = hermnz(zim1, zi, zxxim1, zxxi, sigmax)
      return
   contains
      function hermz(f1, f2, fp1, fp2)
         real(kind=rprec) f1, hermz
         real(kind=rprec) f2
         real(kind=rprec) fp1
         real(kind=rprec) fp2
         hermz = (f2*del1 + f1*del2)/dels - del1*del2*(fp2*(del1 + dels) +&
         &fp1*(del2 + dels))/(6.*dels)
         return
      end function hermz
      function hermnz(f1, f2, fp1, fp2, sigmap)
         real(kind=rprec) f1, hermnz
         real(kind=rprec) f2
         real(kind=rprec) fp1
         real(kind=rprec) fp2
         real(kind=rprec) sigmap
         hermnz = (f2*del1 + f1*del2)/dels + (fp2*del1*(sinhm1 - sinhms) +&
         &fp1*del2*(sinhm2 - sinhms))/(sigmap*sigmap*dels*(1.+sinhms))
         return
      end function hermnz
   end function surf2
   subroutine ceez(del1, del2, sigma, c1, c2, c3, n)

!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer n
      real(kind=rprec) del1, del2, sigma, c1, c2, c3
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      real(kind=rprec)::del, dummy, coshm1, coshm2, delp, delm,&
      &sinhmp, sinhmm, denom
!-----------------------------------------------
!
!
!                                 coded by alan kaylor cline
!                           from fitpack -- january 26, 1987
!                        a curve and surface fitting package
!                      a product of pleasant valley software
!                  8603 altus cove, austin, texas 78759, usa
!
! this subroutine determines the coefficients c1, c2, and c3
! used to determine endpoint slopes. specifically, if
! function values y1, y2, and y3 are given at points x1, x2,
! and x3, respectively, the quantity c1*y1 + c2*y2 + c3*y3
! is the value of the derivative at x1 of a spline under
! tension (with tension factor sigma) passing through the
! three points and having third derivative equal to zero at
! x1. optionally, only two values, c1 and c2 are determined.
!
! on input--
!
!   del1 is x2-x1 (.gt. 0.).
!
!   del2 is x3-x1 (.gt. 0.). if n .eq. 2, this parameter is
!   ignored.
!
!   sigma is the tension factor.
!
! and
!
!   n is a switch indicating the number of coefficients to
!   be returned. if n .eq. 2 only two coefficients are
!   returned. otherwise all three are returned.
!
! on output--
!
!   c1, c2, and c3 contain the coefficients.
!
! none of the input parameters are altered.
!
! this subroutine references package module snhcsh_1.
!
!-----------------------------------------------------------
!
      if (n /= 2) then
         if (sigma == 0.) then
            del = del2 - del1
!
! tension .eq. 0.
!
            c1 = -(del1 + del2)/(del1*del2)
            c2 = del2/(del1*del)
            c3 = -del1/(del2*del)
            return
!
! tension .ne. 0.
!
         end if
         call snhcsh_1(dummy, coshm1, sigma*del1, 1)
         call snhcsh_1(dummy, coshm2, sigma*del2, 1)
         delp = sigma*(del2 + del1)/2.
         delm = sigma*(del2 - del1)/2.
         call snhcsh_1(sinhmp, dummy, delp, -1)
         call snhcsh_1(sinhmm, dummy, delm, -1)
         denom = coshm1*(del2 - del1) - 2.*del1*delp*delm*(1.+sinhmp)*&
         &(1.+sinhmm)
         c1 = 2.*delp*delm*(1.+sinhmp)*(1.+sinhmm)/denom
         c2 = -coshm2/denom
         c3 = coshm1/denom
         return
!
! two coefficients
!
      end if
      c1 = -1./del1
      c2 = -c1
      return
   end subroutine ceez
   subroutine curvpp(n, x, y, p, d, isw, s, eps, ys, ysp, sigma, td,&
   &tsd1, hd, hsd1, hsd2, rd, rsd1, rsd2, rnm1, rn, v, ierr)

!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer n, isw, ierr
      real(kind=rprec) p, s, eps, sigma
      real(kind=rprec), dimension(n) :: x, y, d, ys, ysp, td,&
      &tsd1, hd, hsd1, hsd2,&
      &rd, rsd1, rsd2, rnm1, rn, v
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: nm1, nm2, nm3, i, ip1, im1, ibak
      real(kind=rprec) :: q, sigmap, delxi1, delyi1, dim1, delxi,&
      &delyi, di, hsd11,&
      &sl, su, betapp, betap, alphap, sumd, sumy, disq, alpha, hsd1ip&
      &, beta, con, sum, hsd1p, hdim1, hdi, yspnm1, yspn, rsd1i,&
      &sumnm1, sum2, sumn, rsd2i, rnm1t, rnt, rdn, f, g, rnm1sm, rnsm&
      &, wim2, wim1, tui, wi, h, step
!-----------------------------------------------
!
!
!                                 coded by alan kaylor cline
!                           from fitpack -- january 26, 1987
!                        a curve and surface fitting package
!                      a product of pleasant valley software
!                  8603 altus cove, austin, texas 78759, usa
!
! this subroutine determines the parameters necessary to
! compute a periodic smoothing spline under tension. for a
! given increasing sequence of abscissae (x(i)), i = 1,...,n
! and associated ordinates (y(i)), i = 1,...,n, letting p be
! the period, x(n+1) = x(1)+p, and y(n+1) = y(1), the
! function determined minimizes the summation from i = 1 to
! n of the square of the second derivative of f plus sigma
! squared times the difference of the first derivative of f
! and (f(x(i+1))-f(x(i)))/(x(i+1)-x(i)) squared, over all
! functions f with period p and two continuous derivatives
! such that the summation of the square of
! (f(x(i))-y(i))/d(i) is less than or equal to a given
! constant s, where (d(i)), i = 1,...,n are a given set of
! observation weights. the function determined is a periodic
! spline under tension with third derivative discontinuities
! at (x(i)) i = 1,...,n (and all periodic translations of
! these values). for actual computation of points on the
! curve it is necessary to call the function curvp2.
!
! on input--
!
!   n is the number of values to be smoothed (n.ge.2).
!
!   x is an array of the n increasing abscissae of the
!   values to be smoothed.
!
!   y is an array of the n ordinates of the values to be
!   smoothed, (i. e. y(k) is the functional value
!   corresponding to x(k) ).
!
!   p is the period (p .gt. x(n)-x(1)).
!
!   d is a parameter containing the observation weights.
!   this may either be an array of length n or a scalar
!   (interpreted as a constant). the value of d
!   corresponding to the observation (x(k),y(k)) should
!   be an approximation to the standard deviation of error.
!
!   isw contains a switch indicating whether the parameter
!   d is to be considered a vector or a scalar,
!          = 0 if d is an array of length n,
!          = 1 if d is a scalar.
!
!   s contains the value controlling the smoothing. this
!   must be non-negative. for s equal to zero, the
!   subroutine does interpolation, larger values lead to
!   smoother funtions. if parameter d contains standard
!   deviation estimates, a reasonable value for s is
!   real(n).
!
!   eps contains a tolerance on the relative precision to
!   which s is to be interpreted. this must be greater than
!   or equal to zero and less than equal or equal to one. a
!   reasonable value for eps is sqrt(2./real(n)).
!
!   ys is an array of length at least n.
!
!   ysp is an array of length at least n.
!
!   sigma contains the tension factor. this value indicates
!   the degree to which the first derivative part of the
!   smoothing functional is emphasized. if sigma is nearly
!   zero (e. g. .001) the resulting curve is approximately a
!   cubic spline. if sigma is large (e. g. 50.) the
!   resulting curve is nearly a polygonal line. if sigma
!   equals zero a cubic spline results. a standard value for
!   sigma is approximately 1.
!
! and
!
!   td, tsd1, hd, hsd1, hsd2, rd, rsd1, rsd2, rnm1, rn, and
!   v are arrays of length at least n which are used for
!   scratch storage.
!
! on output--
!
!   ys contains the smoothed ordinate values.
!
!   ysp contains the values of the second derivative of the
!   smoothed curve at the given nodes.
!
!   ierr contains an error flag,
!        = 0 for normal return,
!        = 1 if n is less than 2,
!        = 2 if s is negative,
!        = 3 if eps is negative or greater than one,
!        = 4 if x-values are not strictly increasing,
!        = 5 if a d-value is non-positive,
!        = 6 if p is less than or equal to x(n)-x(1).
!
! and
!
!   n, x, y, d, isw, s, eps, and sigma are unaltered.
!
! this subroutine references package modules terms_1 and
! snhcsh_1.
!
!-----------------------------------------------------------
!
      if (n >= 2) then
         if (s < 0.) go to 26
         if (eps < 0. .or. eps > 1.) go to 27
         if (p <= x(n) - x(1)) go to 30
         ierr = 0
         q = 0.
         rsd1(1) = 0.
         rsd2(1) = 0.
         rsd2(2) = 0.
         rsd1(n - 1) = 0.
         rsd2(n - 1) = 0.
         rsd2(n) = 0.
!
! denormalize tension factor
!
         sigmap = abs(sigma)*real(n)/p
!
! form t matrix and second differences of y into ys
!
         nm1 = n - 1
         nm2 = n - 2
         nm3 = n - 3
         delxi1 = x(1) + p - x(n)
         delyi1 = (y(1) - y(n))/delxi1
         call terms_1(dim1, tsd1(1), sigmap, delxi1)
         hsd1(1) = 1./delxi1
         do i = 1, n
            ip1 = i + 1
            if (i == n) ip1 = 1
            delxi = x(ip1) - x(i)
            if (i == n) delxi = x(1) + p - x(n)
            if (delxi <= 0.) go to 28
            delyi = (y(ip1) - y(i))/delxi
            ys(i) = delyi - delyi1
            call terms_1(di, tsd1(ip1), sigmap, delxi)
            td(i) = di + dim1
            hd(i) = -(1./delxi + 1./delxi1)
            hsd1(ip1) = 1./delxi
            delxi1 = delxi
            delyi1 = delyi
            dim1 = di
         end do
         hsd11 = hsd1(1)
         if (n < 3) then
            tsd1(2) = tsd1(1) + tsd1(2)
            tsd1(1) = 0.
            hsd1(2) = hsd1(1) + hsd1(2)
            hsd1(1) = 0.
!
! calculate lower and upper tolerances
!
         end if
         sl = s*(1.-eps)
         su = s*(1.+eps)
         if (d(1) <= 0.) go to 29
         if (isw /= 1) then
!
! form h matrix - d array
!
            betapp = hsd1(n)*d(n)*d(n)
            betap = hsd1(1)*d(1)*d(1)
            alphap = hd(n)*d(n)*d(n)
            im1 = n
            sumd = 0.
            sumy = 0.
            do i = 1, n
               disq = d(i)*d(i)
               sumd = sumd + 1./disq
               sumy = sumy + y(i)/disq
               ip1 = i + 1
               if (i == n) ip1 = 1
               alpha = hd(i)*disq
               if (d(ip1) <= 0.) go to 29
               hsd1ip = hsd1(ip1)
               if (i == n) hsd1ip = hsd11
               beta = hsd1ip*d(ip1)*d(ip1)
               hd(i) = (hsd1(i)*d(im1))**2 + alpha*hd(i) + beta*hsd1ip
               hsd2(i) = hsd1(i)*betapp
               hsd1(i) = hsd1(i)*(alpha + alphap)
               im1 = i
               alphap = alpha
               betapp = betap
               betap = beta
            end do
            if (n == 3) hsd1(3) = hsd1(3) + hsd2(2)
!
! test for straight line fit
!
            con = sumy/sumd
            sum = 0.
            do i = 1, n
               sum = sum + ((y(i) - con)/d(i))**2
            end do
            if (sum <= su) go to 23
         else
!
! form h matrix - d constant
!
            sl = d(1)*d(1)*sl
            su = d(1)*d(1)*su
            hsd1p = hsd1(n)
            hdim1 = hd(n)
            sumy = 0.
            do i = 1, n
               sumy = sumy + y(i)
               hsd1ip = hsd11
               if (i < n) hsd1ip = hsd1(i + 1)
               hdi = hd(i)
               hd(i) = hsd1(i)*hsd1(i) + hdi*hdi + hsd1ip*hsd1ip
               hsd2(i) = hsd1(i)*hsd1p
               hsd1p = hsd1(i)
               hsd1(i) = hsd1p*(hdi + hdim1)
               hdim1 = hdi
            end do
            if (n == 3) hsd1(3) = hsd1(3) + hsd2(2)
!
! test for straight line fit
!
            con = sumy/real(n)
            sum = 0.
            do i = 1, n
               sum = sum + (y(i) - con)**2
            end do
            if (sum <= su) go to 23
!
! top of iteration
! cholesky factorization of q*t+h into r
!
!
! i = 1
!
         end if
8        continue
         rd(1) = 1./(q*td(1) + hd(1))
         rnm1(1) = hsd2(1)
         yspnm1 = ys(nm1)
         rn(1) = q*tsd1(1) + hsd1(1)
         yspn = ys(n)
         ysp(1) = ys(1)
         rsd1i = q*tsd1(2) + hsd1(2)
         rsd1(2) = rsd1i*rd(1)
         sumnm1 = 0.
         sum2 = 0.
         sumn = 0.
         if (n /= 3) then
            if (n == 2) go to 12
!
! i = 2
!
            rd(2) = 1./(q*td(2) + hd(2) - rsd1i*rsd1(2))
            rnm1(2) = -rnm1(1)*rsd1(2)
            rn(2) = hsd2(2) - rn(1)*rsd1(2)
            ysp(2) = ys(2) - rsd1(2)*ysp(1)
            if (n /= 4) then
               do i = 3, nm2
                  rsd2i = hsd2(i)
                  rsd1i = q*tsd1(i) + hsd1(i) - rsd2i*rsd1(i - 1)
                  rsd2(i) = rsd2i*rd(i - 2)
                  rsd1(i) = rsd1i*rd(i - 1)
                  rd(i) = 1./(q*td(i) + hd(i) - rsd1i*rsd1(i) - rsd2i*rsd2(i))
                  rnm1(i) = (-rnm1(i - 2)*rsd2(i)) - rnm1(i - 1)*rsd1(i)
                  rnm1t = rnm1(i - 2)*rd(i - 2)
                  sumnm1 = sumnm1 + rnm1t*rnm1(i - 2)
                  rnm1(i - 2) = rnm1t
                  sum2 = sum2 + rnm1t*rn(i - 2)
                  yspnm1 = yspnm1 - rnm1t*ysp(i - 2)
                  rn(i) = (-rn(i - 2)*rsd2(i)) - rn(i - 1)*rsd1(i)
                  rnt = rn(i - 2)*rd(i - 2)
                  sumn = sumn + rnt*rn(i - 2)
                  rn(i - 2) = rnt
                  yspn = yspn - rnt*ysp(i - 2)
                  ysp(i) = ys(i) - rsd1(i)*ysp(i - 1) - rsd2(i)*ysp(i - 2)
               end do
!
! i = n-3
!
            end if
            rnm1(nm3) = hsd2(nm1) + rnm1(nm3)
            rnm1(nm2) = rnm1(nm2) - hsd2(nm1)*rsd1(nm2)
            rnm1t = rnm1(nm3)*rd(nm3)
            sumnm1 = sumnm1 + rnm1t*rnm1(nm3)
            rnm1(nm3) = rnm1t
            sum2 = sum2 + rnm1t*rn(nm3)
            yspnm1 = yspnm1 - rnm1t*ysp(nm3)
            rnt = rn(nm3)*rd(nm3)
            sumn = sumn + rnt*rn(nm3)
            rn(nm3) = rnt
            yspn = yspn - rnt*ysp(nm3)
!
! i = n-2
!
         end if
         rnm1(nm2) = q*tsd1(nm1) + hsd1(nm1) + rnm1(nm2)
         rnm1t = rnm1(nm2)*rd(nm2)
         sumnm1 = sumnm1 + rnm1t*rnm1(nm2)
         rnm1(nm2) = rnm1t
         rn(nm2) = hsd2(n) + rn(nm2)
         sum2 = sum2 + rnm1t*rn(nm2)
         yspnm1 = yspnm1 - rnm1t*ysp(nm2)
         rnt = rn(nm2)*rd(nm2)
         sumn = sumn + rnt*rn(nm2)
         rn(nm2) = rnt
         yspn = yspn - rnt*ysp(nm2)
!
! i = n-1
!
12       continue
         rd(nm1) = 1./(q*td(nm1) + hd(nm1) - sumnm1)
         ysp(nm1) = yspnm1
         rn(nm1) = q*tsd1(n) + hsd1(n) - sum2
         rnt = rn(nm1)*rd(nm1)
         sumn = sumn + rnt*rn(nm1)
         rn(nm1) = rnt
         yspn = yspn - rnt*ysp(nm1)
!
! i = n
!
         rdn = q*td(n) + hd(n) - sumn
         rd(n) = 0.
         if (rdn > 0.) rd(n) = 1./rdn
         ysp(n) = yspn
!
! back solve of r(transpose)* r * ysp = ys
!
         ysp(n) = rd(n)*ysp(n)
         ysp(nm1) = rd(nm1)*ysp(nm1) - rn(nm1)*ysp(n)
         if (n /= 2) then
            yspn = ysp(n)
            yspnm1 = ysp(nm1)
            do ibak = 1, nm2
               i = nm1 - ibak
               ysp(i) = rd(i)*ysp(i) - rsd1(i + 1)*ysp(i + 1) - rsd2(i + 2)*&
               &ysp(i + 2) - rnm1(i)*yspnm1 - rn(i)*yspn
            end do
         end if
         sum = 0.
         delyi1 = (ysp(1) - ysp(n))/(x(1) + p - x(n))
         if (isw /= 1) then
!
! calculation of residual norm
!  - d array
!
            do i = 1, nm1
               delyi = (ysp(i + 1) - ysp(i))/(x(i + 1) - x(i))
               v(i) = (delyi - delyi1)*d(i)*d(i)
               sum = sum + v(i)*(delyi - delyi1)
               delyi1 = delyi
            end do
            delyi = (ysp(1) - ysp(n))/(x(1) + p - x(n))
            v(n) = (delyi - delyi1)*d(n)*d(n)
         else
!
! calculation of residual norm
!  - d constant
!
            do i = 1, nm1
               delyi = (ysp(i + 1) - ysp(i))/(x(i + 1) - x(i))
               v(i) = delyi - delyi1
               sum = sum + v(i)*(delyi - delyi1)
               delyi1 = delyi
            end do
            delyi = (ysp(1) - ysp(n))/(x(1) + p - x(n))
            v(n) = delyi - delyi1
         end if
         sum = sum + v(n)*(delyi - delyi1)
!
! test for convergence
!
         if (sum <= su .and. sum >= sl .and. q > 0.) go to 21
!
! calculation of newton correction
!
         f = 0.
         g = 0.
         rnm1sm = 0.
         rnsm = 0.
         im1 = n
         if (n /= 2) then
            wim2 = 0.
            wim1 = 0.
            do i = 1, nm2
               tui = tsd1(i)*ysp(im1) + td(i)*ysp(i) + tsd1(i + 1)*ysp(i + 1)
               wi = tui - rsd1(i)*wim1 - rsd2(i)*wim2
               rnm1sm = rnm1sm - rnm1(i)*wi
               rnsm = rnsm - rn(i)*wi
               f = f + tui*ysp(i)
               g = g + wi*wi*rd(i)
               im1 = i
               wim2 = wim1
               wim1 = wi
            end do
         end if
         tui = tsd1(nm1)*ysp(im1) + td(nm1)*ysp(nm1) + tsd1(n)*ysp(n)
         wi = tui + rnm1sm
         f = f + tui*ysp(nm1)
         g = g + wi*wi*rd(nm1)
         tui = tsd1(n)*ysp(nm1) + td(n)*ysp(n) + tsd1(1)*ysp(1)
         wi = tui + rnsm - rn(nm1)*wi
         f = f + tui*ysp(n)
         g = g + wi*wi*rd(n)
         h = f - q*g
         if (h <= 0. .and. q > 0.) go to 21
!
! update q - newton step
!
         step = (sum - sqrt(sum*sl))/h
         if (sl /= 0.) step = step*sqrt(sum/sl)
         q = q + step
         go to 8
!
! store smoothed y-values and second derivatives
!
21       continue
         ys = y - v
         ysp = q*ysp
         return
!
! store constant ys and zero ysp
!
23       continue
         ys = con
         ysp = 0.
         return
!
! n less than 2
!
      end if
      ierr = 1
      return
!
! s negative
!
26    continue
      ierr = 2
      return
!
! eps negative or greater than 1
!
27    continue
      ierr = 3
      return
!
! x-values not strictly increasing
!
28    continue
      ierr = 4
      return
!
! weight non-positive
!
29    continue
      ierr = 5
      return
!
! incorrect period
!
30    continue
      ierr = 6
      return
   end subroutine curvpp
   subroutine curvss(n, x, y, d, isw, s, eps, ys, ysp, sigma, td,&
   &tsd1, hd, hsd1, hsd2, rd, rsd1, rsd2, v, ierr)

!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer n, isw, ierr
      real(kind=rprec) s, eps, sigma
      real(kind=rprec), dimension(n) :: x, y, d, ys, ysp,&
      &td, tsd1, hd, hsd1, hsd2,&
      &rd, rsd1, rsd2, v
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: nm1, nm3, i, ibak
      real(kind=rprec) :: p, rdim1, yspim2, sigmap, delxi1,&
      &delyi1, dim1, delxi,&
      &delyi, di, sl, su, betapp, betap, alphap, alpha, beta, hsd1p,&
      &hdim1, hdi, rsd2i, rsd1i, sum, f, g, wim2, wim1, tui, wi, h,&
      &step
!-----------------------------------------------
!
!
!                                 coded by alan kaylor cline
!                           from fitpack -- january 26, 1987
!                        a curve and surface fitting package
!                      a product of pleasant valley software
!                  8603 altus cove, austin, texas 78759, usa
!
! this subroutine determines the parameters necessary to
! compute a smoothing spline under tension. for a given
! increasing sequence of abscissae (x(i)), i = 1,..., n and
! associated ordinates (y(i)), i = 1,..., n, the function
! determined minimizes the summation from i = 1 to n-1 of
! the square of the second derivative of f plus sigma
! squared times the difference of the first derivative of f
! and (f(x(i+1))-f(x(i)))/(x(i+1)-x(i)) squared, over all
! functions f with two continuous derivatives such that the
! summation of the square of (f(x(i))-y(i))/d(i) is less
! than or equal to a given constant s, where (d(i)), i = 1,
! ..., n are a given set of observation weights. the
! function determined is a spline under tension with third
! derivative discontinuities at (x(i)), i = 2,..., n-1. for
! actual computation of points on the curve it is necessary
! to call the function curv2.
!
! on input--
!
!   n is the number of values to be smoothed (n.ge.2).
!
!   x is an array of the n increasing abscissae of the
!   values to be smoothed.
!
!   y is an array of the n ordinates of the values to be
!   smoothed, (i. e. y(k) is the functional value
!   corresponding to x(k) ).
!
!   d is a parameter containing the observation weights.
!   this may either be an array of length n or a scalar
!   (interpreted as a constant). the value of d
!   corresponding to the observation (x(k),y(k)) should
!   be an approximation to the standard deviation of error.
!
!   isw contains a switch indicating whether the parameter
!   d is to be considered a vector or a scalar,
!          = 0 if d is an array of length n,
!          = 1 if d is a scalar.
!
!   s contains the value controlling the smoothing. this
!   must be non-negative. for s equal to zero, the
!   subroutine does interpolation, larger values lead to
!   smoother funtions. if parameter d contains standard
!   deviation estimates, a reasonable value for s is
!   real(n).
!
!   eps contains a tolerance on the relative precision to
!   which s is to be interpreted. this must be greater than
!   or equal to zero and less than equal or equal to one. a
!   reasonable value for eps is sqrt(2./real(n)).
!
!   ys is an array of length at least n.
!
!   ysp is an array of length at least n.
!
!   sigma contains the tension factor. this value indicates
!   the degree to which the first derivative part of the
!   smoothing functional is emphasized. if sigma is nearly
!   zero (e. g. .001) the resulting curve is approximately a
!   cubic spline. if sigma is large (e. g. 50.) the
!   resulting curve is nearly a polygonal line. if sigma
!   equals zero a cubic spline results. a standard value for
!   sigma is approximately 1.
!
! and
!
!   td, tsd1, hd, hsd1, hsd2, rd, rsd1, rsd2, and v are
!   arrays of length at least n which are used for scratch
!   storage.
!
! on output--
!
!   ys contains the smoothed ordinate values.
!
!   ysp contains the values of the second derivative of the
!   smoothed curve at the given nodes.
!
!   ierr contains an error flag,
!        = 0 for normal return,
!        = 1 if n is less than 2,
!        = 2 if s is negative,
!        = 3 if eps is negative or greater than one,
!        = 4 if x-values are not strictly increasing,
!        = 5 if a d-value is non-positive.
!
! and
!
!   n, x, y, d, isw, s, eps, and sigma are unaltered.
!
! this subroutine references package modules terms_1 and
! snhcsh_1.
!
!-----------------------------------------------------------
!
      if (n >= 2) then
         if (s < 0.) go to 17
         if (eps < 0. .or. eps > 1.) go to 18
         ierr = 0
         p = 0.
         v(1) = 0.
         v(n) = 0.
         ysp(1) = 0.
         ysp(n) = 0.
         if (n /= 2) then
            rsd1(1) = 0.
            rd(1) = 0.
            rsd2(n) = 0.
            rdim1 = 0.
            yspim2 = 0.
!
! denormalize tension factor
!
            sigmap = abs(sigma)*real(n - 1)/(x(n) - x(1))
!
! form t matrix and second differences of y into ys
!
            nm1 = n - 1
            nm3 = n - 3
            delxi1 = 1.
            delyi1 = 0.
            dim1 = 0.
            do i = 1, nm1
               delxi = x(i + 1) - x(i)
               if (delxi <= 0.) go to 19
               delyi = (y(i + 1) - y(i))/delxi
               ys(i) = delyi - delyi1
               call terms_1(di, tsd1(i + 1), sigmap, delxi)
               td(i) = di + dim1
               hd(i) = -(1./delxi + 1./delxi1)
               hsd1(i + 1) = 1./delxi
               delxi1 = delxi
               delyi1 = delyi
               dim1 = di
            end do
!
! calculate lower and upper tolerances
!
            sl = s*(1.-eps)
            su = s*(1.+eps)
            if (isw /= 1) then
!
! form h matrix - d array
!
               if (d(1) <= 0. .or. d(2) <= 0.) go to 20
               betapp = 0.
               betap = 0.
               alphap = 0.
               do i = 2, nm1
                  alpha = hd(i)*d(i)*d(i)
                  if (d(i + 1) <= 0.) go to 20
                  beta = hsd1(i + 1)*d(i + 1)*d(i + 1)
                  hd(i) = (hsd1(i)*d(i - 1))**2 + alpha*hd(i) + beta*hsd1(i + 1)
                  hsd2(i) = hsd1(i)*betapp
                  hsd1(i) = hsd1(i)*(alpha + alphap)
                  alphap = alpha
                  betapp = betap
                  betap = beta
               end do
            else
!
! form h matrix - d constant
!
               if (d(1) <= 0.) go to 20
               sl = d(1)*d(1)*sl
               su = d(1)*d(1)*su
               hsd1p = 0.
               hdim1 = 0.
               do i = 2, nm1
                  hdi = hd(i)
                  hd(i) = hsd1(i)*hsd1(i) + hdi*hdi + hsd1(i + 1)*hsd1(i + 1)
                  hsd2(i) = hsd1(i)*hsd1p
                  hsd1p = hsd1(i)
                  hsd1(i) = hsd1p*(hdi + hdim1)
                  hdim1 = hdi
               end do
!
! top of iteration
! cholesky factorization of p*t+h into r
!
            end if
5           continue
            do i = 2, nm1
               rsd2i = hsd2(i)
               rsd1i = p*tsd1(i) + hsd1(i) - rsd2i*rsd1(i - 1)
               rsd2(i) = rsd2i*rdim1
               rdim1 = rd(i - 1)
               rsd1(i) = rsd1i*rdim1
               rd(i) = 1./(p*td(i) + hd(i) - rsd1i*rsd1(i) - rsd2i*rsd2(i))
               ysp(i) = ys(i) - rsd1(i)*ysp(i - 1) - rsd2(i)*yspim2
               yspim2 = ysp(i - 1)
            end do
!
! back solve of r(transpose)* r * ysp = ys
!
            ysp(nm1) = rd(nm1)*ysp(nm1)
            if (n /= 3) then
               do ibak = 1, nm3
                  i = nm1 - ibak
                  ysp(i) = rd(i)*ysp(i) - rsd1(i + 1)*ysp(i + 1) - rsd2(i + 2)&
                  &*ysp(i + 2)
               end do
            end if
            sum = 0.
            delyi1 = 0.
            if (isw /= 1) then
!
! calculation of residual norm
!  - d array
!
               do i = 1, nm1
                  delyi = (ysp(i + 1) - ysp(i))/(x(i + 1) - x(i))
                  v(i) = (delyi - delyi1)*d(i)*d(i)
                  sum = sum + v(i)*(delyi - delyi1)
                  delyi1 = delyi
               end do
               v(n) = -delyi1*d(n)*d(n)
            else
!
! calculation of residual norm
!  - d constant
!
               do i = 1, nm1
                  delyi = (ysp(i + 1) - ysp(i))/(x(i + 1) - x(i))
                  v(i) = delyi - delyi1
                  sum = sum + v(i)*(delyi - delyi1)
                  delyi1 = delyi
               end do
               v(n) = -delyi1
            end if
            sum = sum - v(n)*delyi1
!
! test for convergence
!
            if (sum <= su) go to 14
!
! calculation of newton correction
!
            f = 0.
            g = 0.
            wim2 = 0.
            wim1 = 0.
            do i = 2, nm1
               tui = tsd1(i)*ysp(i - 1) + td(i)*ysp(i) + tsd1(i + 1)*ysp(i + 1)
               wi = tui - rsd1(i)*wim1 - rsd2(i)*wim2
               f = f + tui*ysp(i)
               g = g + wi*wi*rd(i)
               wim2 = wim1
               wim1 = wi
            end do
            h = f - p*g
            if (h <= 0.) go to 14
!
! update p - newton step
!
            step = (sum - sqrt(sum*sl))/h
            if (sl /= 0.) step = step*sqrt(sum/sl)
            p = p + step
            go to 5
!
! store smoothed y-values and second derivatives
!
         end if
14       continue
         ys = y - v
         ysp = p*ysp
         return
!
! n less than 2
!
      end if
      ierr = 1
      return
!
! s negative
!
17    continue
      ierr = 2
      return
!
! eps negative or greater than 1
!
18    continue
      ierr = 3
      return
!
! x-values not strictly increasing
!
19    continue
      ierr = 4
      return
!
! weight non-positive
!
20    continue
      ierr = 5
      return
   end subroutine curvss
   integer function intrvl_1(t, x, n)

!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer n
      real(kind=rprec) t
      real(kind=rprec), dimension(n) :: x
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: i, il, ih
      real(kind=rprec) :: tt

      save i
!-----------------------------------------------
!
!
!                                 coded by alan kaylor cline
!                           from fitpack -- january 26, 1987
!                        a curve and surface fitting package
!                      a product of pleasant valley software
!                  8603 altus cove, austin, texas 78759, usa
!
! this function determines the index of the interval
! (determined by a given increasing sequence) in which
! a given value lies.
!
! on input--
!
!   t is the given value.
!
!   x is a vector of strictly increasing values.
!
! and
!
!   n is the length of x (n .ge. 2).
!
! on output--
!
!   intrvl_1 returns an integer i such that
!
!          i =  1       if         e   t .lt. x(2)  ,
!          i =  n-1     if x(n-1) .le. t            ,
!          otherwise       x(i)  .le. t .le. x(i+1),
!
! none of the input parameters are altered.
!
!-----------------------------------------------------------
!
      data i/1/
!
      tt = t
!
! check for illegal i
!
      if (i >= n) i = n/2
!
! check old interval and extremes
!
      if (tt < x(i)) then
         if (tt <= x(2)) then
            i = 1
            intrvl_1 = 1
            return
         else
            il = 2
            ih = i
         end if
      else if (tt <= x(i + 1)) then
         intrvl_1 = i
         return
      else if (tt >= x(n - 1)) then
         i = n - 1
         intrvl_1 = n - 1
         return
      else
         il = i + 1
         ih = n - 1
      end if
!
! binary search loop
!
1     continue
      i = (il + ih)/2
      if (tt < x(i)) then
         ih = i
      else if (tt > x(i + 1)) then
         il = i + 1
      else
         intrvl_1 = i
         return
      end if
      go to 1
   end function intrvl_1
   integer function intrvp(t, x, n, p, tp)

!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer n
      real(kind=rprec) t, p, tp
      real(kind=rprec), dimension(n) :: x
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: i, nper, il, ih
      real(kind=rprec) :: tt

      save i
!-----------------------------------------------
!
!
!                                 coded by alan kaylor cline
!                           from fitpack -- january 26, 1987
!                        a curve and surface fitting package
!                      a product of pleasant valley software
!                  8603 altus cove, austin, texas 78759, usa
!
! this function determines the index of the interval
! (determined by a given increasing sequence) in which a
! given value lies, after translating the value to within
! the correct period.  it also returns this translated value.
!
! on input--
!
!   t is the given value.
!
!   x is a vector of strictly increasing values.
!
!   n is the length of x (n .ge. 2).
!
! and
!
!   p contains the period.
!
! on output--
!
!   tp contains a translated value of t (i. e. x(1) .le. tp,
!   tp .lt. x(1)+p, and tp = t + k*p for some integer k).
!
!   intrvl_1 returns an integer i such that
!
!          i = 1       if             tp .lt. x(2)  ,
!          i = n       if   x(n) .le. tp            ,
!          otherwise       x(i)  .le. tp .lt. x(i+1),
!
! none of the input parameters are altered.
!
!-----------------------------------------------------------
!
      data i/1/
!
      nper = (t - x(1))/p
      tp = t - real(nper)*p
      if (tp < x(1)) tp = tp + p
      tt = tp
!
! check for illegal i
!
      if (i >= n) i = n/2
!
! check old interval and extremes
!
      if (tt < x(i)) then
         if (tt <= x(2)) then
            i = 1
            intrvp = 1
            return
         else
            il = 2
            ih = i
         end if
      else if (tt <= x(i + 1)) then
         intrvp = i
         return
      else if (tt >= x(n)) then
         i = n
         intrvp = n
         return
      else
         il = i + 1
         ih = n
      end if
!
! binary search loop
!
1     continue
      i = (il + ih)/2
      if (tt < x(i)) then
         ih = i
      else if (tt > x(i + 1)) then
         il = i + 1
      else
         intrvp = i
         return
      end if
      go to 1
   end function intrvp
   subroutine snhcsh_1(sinhm, coshm, x, isw)

!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer isw
      real(kind=rprec) sinhm, coshm, x
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      real(kind=rprec) :: sp14, sp13, sp12, sp11, sq12, sq11,&
      &sq10, sp25, sp24, sp23&
      &, sp22, sp21, sq22, sq21, sq20, sp35, sp34, sp33, sp32, sp31,&
      &sq32, sq31, sq30, sp45, sp44, sp43, sp42, sp41, sq42, sq41,&
      &sq40, cp5, cp4, cp3, cp2, cp1, cq2, cq1, cq0, zp4, zp3, zp2,&
      &zp1, zq2, zq1, zq0, ax, xs, expx
!-----------------------------------------------
!
!
!                                 coded by alan kaylor cline
!                           from fitpack -- january 26, 1987
!                        a curve and surface fitting package
!                      a product of pleasant valley software
!                  8603 altus cove, austin, texas 78759, usa
!
! this subroutine returns approximations to
!       sinhm(x) = sinh(x)/x-1
!       coshm(x) = cosh(x)-1
! and
!       coshmm(x) = (cosh(x)-1-x*x/2)/(x*x)
! with relative error less than 4.0e-14.
!
! on input--
!
!   x contains the value of the independent variable.
!
!   isw indicates the function desired
!           = -1 if only sinhm is desired,
!           =  0 if both sinhm and coshm are desired,
!           =  1 if only coshm is desired,
!           =  2 if only coshmm is desired,
!           =  3 if both sinhm and coshmm are desired.
!
! on output--
!
!   sinhm contains the value of sinhm(x) if isw .le. 0 or
!   isw .eq. 3 (sinhm is unaltered if isw .eq.1 or isw .eq.
!   2).
!
!   coshm contains the value of coshm(x) if isw .eq. 0 or
!   isw .eq. 1 and contains the value of coshmm(x) if isw
!   .ge. 2 (coshm is unaltered if isw .eq. -1).
!
! and
!
!   x and isw are unaltered.
!
!-----------------------------------------------------------
!
      data sp14/.227581660976348E-7/
      data sp13/.612189863171694E-5/
      data sp12/.715314759211209E-3/
      data sp11/.398088289992973E-1/
      data sq12/.206382701413725E-3/
      data sq11/-.611470260009508E-1/
      data sq10/.599999999999986E+1/
      data sp25/.129094158037272E-9/
      data sp24/.473731823101666E-7/
      data sp23/.849213455598455E-5/
      data sp22/.833264803327242E-3/
      data sp21/.425024142813226E-1/
      data sq22/.106008515744821E-3/
      data sq21/-.449855169512505E-1/
      data sq20/.600000000268619E+1/
      data sp35/.155193945864942E-9/
      data sp34/.511529451668737E-7/
      data sp33/.884775635776784E-5/
      data sp32/.850447617691392E-3/
      data sp31/.428888148791777E-1/
      data sq32/.933128831061610E-4/
      data sq31/-.426677570538507E-1/
      data sq30/.600000145086489E+1/
      data sp45/.188070632058331E-9/
      data sp44/.545792817714192E-7/
      data sp43/.920119535795222E-5/
      data sp42/.866559391672985E-3/
      data sp41/.432535234960858E-1/
      data sq42/.824891748820670E-4/
      data sq41/-.404938841672262E-1/
      data sq40/.600005006283834E+1/
      data cp5/.552200614584744E-9/
      data cp4/.181666923620944E-6/
      data cp3/.270540125846525E-4/
      data cp2/.206270719503934E-2/
      data cp1/.744437205569040E-1/
      data cq2/.514609638642689E-4/
      data cq1/-.177792255528382E-1/
      data cq0/.200000000000000E+1/
      data zp4/.664418805876835E-8/
      data zp3/.218274535686385E-5/
      data zp2/.324851059327161E-3/
      data zp1/.244515150174258E-1/
      data zq2/.616165782306621E-3/
      data zq1/-.213163639579425E0/
      data zq0/.240000000000000E+2/
!
      ax = abs(x)
      if (isw < 0) then
!
! sinhm approximation
!
         if (ax <= 3.9) then
            xs = ax*ax
            if (ax <= 2.2) then
!
! sinhm approximation on (0.,2.2)
!
               sinhm = xs*((((sp14*xs + sp13)*xs + sp12)*xs + sp11)*xs&
               &+ 1.)/((sq12*xs + sq11)*xs + sq10)
               return
!
! sinhm approximation on (2.2,3.9)
!
            end if
            sinhm = xs*(((((sp25*xs + sp24)*xs + sp23)*xs + sp22)*xs +&
            &sp21)*xs + 1.)/((sq22*xs + sq21)*xs + sq20)
            return
         end if
         if (ax <= 5.1) then
!
! sinhm approximation on (3.9,5.1)
!
            xs = ax*ax
            sinhm = xs*(((((sp35*xs + sp34)*xs + sp33)*xs + sp32)*xs +&
            &sp31)*xs + 1.)/((sq32*xs + sq31)*xs + sq30)
            return
         end if
         if (ax <= 6.1) then
!
! sinhm approximation on (5.1,6.1)
!
            xs = ax*ax
            sinhm = xs*(((((sp45*xs + sp44)*xs + sp43)*xs + sp42)*xs +&
            &sp41)*xs + 1.)/((sq42*xs + sq41)*xs + sq40)
            return
!
! sinhm approximation above 6.1
!
         end if
         expx = exp(ax)
         sinhm = (expx - 1./expx)/(ax + ax) - 1.
         return
!
! coshm and (possibly) sinhm approximation
!
      end if
      if (isw < 2) then
         if (ax <= 2.2) then
            xs = ax*ax
            coshm = xs*(((((cp5*xs + cp4)*xs + cp3)*xs + cp2)*xs + cp1)*&
            &xs + 1.)/((cq2*xs + cq1)*xs + cq0)
            if (isw == 0) sinhm = xs*((((sp14*xs + sp13)*xs + sp12)*xs&
            &+ sp11)*xs + 1.)/((sq12*xs + sq11)*xs + sq10)
            return
         end if
         expx = exp(ax)
         coshm = (expx + 1./expx)/2.-1.
         if (isw == 0) sinhm = (expx - 1./expx)/(ax + ax) - 1.
         return
!
! coshmm and (possibly) sinhm approximation
!
      end if
      xs = ax*ax
      if (ax <= 2.2) then
         coshm = xs*((((zp4*xs + zp3)*xs + zp2)*xs + zp1)*xs + 1.)/((zq2&
         &*xs + zq1)*xs + zq0)
         if (isw == 3) sinhm = xs*((((sp14*xs + sp13)*xs + sp12)*xs +&
         &sp11)*xs + 1.)/((sq12*xs + sq11)*xs + sq10)
         return
      end if
      expx = exp(ax)
      coshm = ((expx + 1./expx - xs)/2.-1.)/xs
      if (isw == 3) sinhm = (expx - 1./expx)/(ax + ax) - 1.
      return
   end subroutine snhcsh_1
   subroutine terms_1(diag, sdiag, sigma, del)

!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      real(kind=rprec) diag, sdiag, sigma, del
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      real(kind=rprec) :: sigdel, sinhm, coshm, denom
!-----------------------------------------------
!
!
!                                 coded by alan kaylor cline
!                           from fitpack -- january 26, 1987
!                        a curve and surface fitting package
!                      a product of pleasant valley software
!                  8603 altus cove, austin, texas 78759, usa
!
! this subroutine computes the diagonal and superdiagonal
! terms of the tridiagonal linear system associated with
! spline under tension interpolation.
!
! on input--
!
!   sigma contains the tension factor.
!
! and
!
!   del contains the step size.
!
! on output--
!
!                sigma*del*cosh(sigma*del) - sinh(sigma*del)
!   diag = del*--------------------------------------------.
!                     (sigma*del)**2 * sinh(sigma*del)
!
!                   sinh(sigma*del) - sigma*del
!   sdiag = del*----------------------------------.
!                (sigma*del)**2 * sinh(sigma*del)
!
! and
!
!   sigma and del are unaltered.
!
! this subroutine references package module snhcsh_1.
!
!-----------------------------------------------------------
!
      if (sigma == 0.) then
         diag = del/3.
         sdiag = del/6.
         return
      end if
      sigdel = sigma*del
      call snhcsh_1(sinhm, coshm, sigdel, 0)
      denom = sigma*sigdel*(1.+sinhm)
      diag = (coshm - sinhm)/denom
      sdiag = sinhm/denom
      return
   end subroutine terms_1

end module fitpack