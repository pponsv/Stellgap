subroutine toFourier
   use fourier_lib
   use fitpack
   implicit none
!
!   Do Fourier transform integrations needed to convert data on a
!   theta, zeta grid [stored in array f(i=1,nznt)] to a set
!   of Fourier amplitudes [stored in array fnm(mn=1,mnmx)].
!   Typically, the number of grid points in each direction needs
!   to be > 3*number of modes used in each direction to avoid
!   aliasing errors(implies nznt > 9*mnmx).
!
!      do mn=1,mnmx     !loop over Fourier modes
!       fnm(mn) = 0.
!      do i=1,nznt      !loop over theta,zeta grid
!       if(sin_type .eq. 1 .and. cos_type .eq. 0) then
!         fnm(mn) = fnm(mn) + f(i)*sin_toF(i,mn)
!       else if(sin_type .eq. 0 .and. cos_type .eq. 1) then
!         fnm(mn) = fnm(mn) + f(i)*cos_toF(i,mn)
!       endif
!      end do
!      end do
!
   if (sin_type .eq. 1 .and. cos_type .eq. 0) then
      fnm = matmul(f, sin_toF)
   else if (sin_type .eq. 0 .and. cos_type .eq. 1) then
      fnm = matmul(f, cos_toF)
   end if
!
   return
end
!
subroutine old_toFourier
   use fourier_lib
   implicit none
   real(kind=rprec) dum, dnorm
!
!   Do Fourier transform integrations needed to convert data on a
!   theta, zeta grid [stored in array f(i=1,nznt)] to a set
!   of Fourier amplitudes [stored in array fnm(mn=1,mnmx)].
!   Typically, the number of grid points in each direction needs
!   to be > 3*number of modes used in each direction to avoid
!   aliasing errors(implies nznt > 9*mnmx).
!
   do mn = 1, mnmx     !loop over Fourier modes
      fnm(mn) = 0.
!       dnorm = 2.*real(nfp)/real(nznt)
      dnorm = 2./real(nznt)
      dum = abs(rn(mn)) + abs(rm(mn))
      if (nint(dum) .eq. 0) dnorm = .5*dnorm
      do i = 1, nznt      !loop over theta,zeta grid
!       arg = -rn(mn)*ztgrd(i) + rm(mn)*thtgrd(i)
         if (sin_type .eq. 1 .and. cos_type .eq. 0) then
            fnm(mn) = fnm(mn) + f(i)*sin_ar(i, mn)*dnorm
         else if (sin_type .eq. 0 .and. cos_type .eq. 1) then
            fnm(mn) = fnm(mn) + f(i)*cos_ar(i, mn)*dnorm
         end if
      end do
   end do
!
   return
end
!
subroutine toReal
!
!    Convert Fourier mode representation [stored in array anm(mn=1,mnmx)]
!    to values of function on a regularly spaced 2D grid
!    [stored in array f(i=1,nznt)].
!
   use fourier_lib
   implicit none
   do i = 1, nznt
      f(i) = 0.
      do mn = 1, mnmx
         if (sin_type .eq. 1 .and. cos_type .eq. 0) then
            f(i) = f(i) + anm(mn)*sin_ar(i, mn)
         else if (sin_type .eq. 0 .and. cos_type .eq. 1) then
            f(i) = f(i) + anm(mn)*cos_ar(i, mn)
         end if
      end do
   end do
end
!
!
!
subroutine dbydth
!
!    Take the theta derivative of the input Fourier amplitude array, fnm
!    and place the result in the output Fourier amplitude array, anm.
!    Changes to the sin/cos parity are reflected through the sin_type and
!    cos_type variables.
!
   use fourier_lib
   implicit none
   do i = 1, mnmx
      if (sin_type .eq. 1 .and. cos_type .eq. 0) then
         anm(i) = rm(i)*fnm(i)
      else if (sin_type .eq. 0 .and. cos_type .eq. 1) then
         anm(i) = -rm(i)*fnm(i)
      end if
   end do
   if (sin_type .eq. 1 .and. cos_type .eq. 0) then
      sin_type = 0; cos_type = 1
   else if (sin_type .eq. 0 .and. cos_type .eq. 1) then
      sin_type = 1; cos_type = 0
   end if
end
!
!
!
subroutine dbydzt
!
!    Take the zeta derivative of the input Fourier amplitude array, fnm
!    and place the result in the output Fourier amplitude array, anm.
!    Changes to the sin/cos parity are reflected through the sin_type and
!    cos_type variables.
!
   use fourier_lib
   implicit none
   do i = 1, mnmx
      if (sin_type .eq. 1 .and. cos_type .eq. 0) then
         anm(i) = -rn(i)*fnm(i)
      else if (sin_type .eq. 0 .and. cos_type .eq. 1) then
         anm(i) = rn(i)*fnm(i)
      end if
   end do
   if (sin_type .eq. 1 .and. cos_type .eq. 0) then
      sin_type = 0; cos_type = 1
   else if (sin_type .eq. 0 .and. cos_type .eq. 1) then
      sin_type = 1; cos_type = 0
   end if
end
