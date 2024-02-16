module fourier_lib_convolve

   use constants, only: r8
   ! use fourier_lib
   ! use fitpack
   ! use globals

   implicit none

contains

   function toFourier_new(fin, trigtype) result(fout)
      use globals, only: cos_toF, sin_toF

      real(r8), intent(in) :: fin(:,:)
      character*1, intent(in) :: trigtype
      real(r8), allocatable, dimension(:) :: fout, ftemp

      ftemp = reshape(transpose(fin), [size(fin)])
      select case (trigtype)
       case('c')
         fout = matmul(ftemp, cos_toF)
       case('s')
         fout = matmul(ftemp, sin_toF)
      end select

   end function toFourier_new


!    subroutine toFourier
! !
! !   Do Fourier transform integrations needed to convert data on a
! !   theta, zeta grid [stored in array f(i=1,nznt)] to a set
! !   of Fourier amplitudes [stored in array fnm(mn=1,mnmx)].
! !   Typically, the number of grid points in each direction needs
! !   to be > 3*number of modes used in each direction to avoid
! !   aliasing errors(implies nznt > 9*mnmx).
! !
! !      do mn=1,mnmx     !loop over Fourier modes
! !       fnm(mn) = 0.
! !      do i=1,nznt      !loop over theta,zeta grid
! !       if(sin_type .eq. 1 .and. cos_type .eq. 0) then
! !         fnm(mn) = fnm(mn) + f(i)*sin_toF(i,mn)
! !       else if(sin_type .eq. 0 .and. cos_type .eq. 1) then
! !         fnm(mn) = fnm(mn) + f(i)*cos_toF(i,mn)
! !       endif
! !      end do
! !      end do
! !
!       if (sin_type .eq. 1 .and. cos_type .eq. 0) then
!          fnm = matmul(f, sin_toF)
!       else if (sin_type .eq. 0 .and. cos_type .eq. 1) then
!          fnm = matmul(f, cos_toF)
!       end if
! !
!       return
!    end subroutine toFourier
!
!    subroutine old_toFourier

!       real(r8) :: dum, dnorm
! !
! !   Do Fourier transform integrations needed to convert data on a
! !   theta, zeta grid [stored in array f(i=1,nznt)] to a set
! !   of Fourier amplitudes [stored in array fnm(mn=1,mnmx)].
! !   Typically, the number of grid points in each direction needs
! !   to be > 3*number of modes used in each direction to avoid
! !   aliasing errors(implies nznt > 9*mnmx).
! !
!       do mn = 1, mnmx     !loop over Fourier modes
!          fnm(mn) = 0.
! !       dnorm = 2.*real(nfp)/real(nznt)
!          dnorm = 2./real(ith*izt)
!          dum = abs(rn(mn)) + abs(rm(mn))
!          if (nint(dum) .eq. 0) dnorm = .5*dnorm
!          do i = 1, ith*izt      !loop over theta,zeta grid
! !       arg = -rn(mn)*ztgrd(i) + rm(mn)*thtgrd(i)
!             if (sin_type .eq. 1 .and. cos_type .eq. 0) then
!                fnm(mn) = fnm(mn) + f(i)*sin_ar(i, mn)*dnorm
!             else if (sin_type .eq. 0 .and. cos_type .eq. 1) then
!                fnm(mn) = fnm(mn) + f(i)*cos_ar(i, mn)*dnorm
!             end if
!          end do
!       end do
! !
!       return
!    end subroutine old_toFourier
!
!    subroutine toReal
! !
! !    Convert Fourier mode representation [stored in array anm(mn=1,mnmx)]
! !    to values of function on a regularly spaced 2D grid
! !    [stored in array f(i=1,nznt)].
! !
!       integer :: i, mn

!       do i = 1, ith*izt
!          f(i) = 0.
!          do mn = 1, mnmx
!             if (sin_type .eq. 1 .and. cos_type .eq. 0) then
!                f(i) = f(i) + anm(mn)*sin_ar(i, mn)
!             else if (sin_type .eq. 0 .and. cos_type .eq. 1) then
!                f(i) = f(i) + anm(mn)*cos_ar(i, mn)
!             end if
!          end do
!       end do
!    end subroutine toReal
!
!
!

   subroutine dbydth(sin_type, cos_type)
      use globals, only: mnmx, m_fourier, anm, fnm
      integer, intent(inout) :: sin_type, cos_type
!
!    Take the theta derivative of the input Fourier amplitude array, fnm
!    and place the result in the output Fourier amplitude array, anm.
!    Changes to the sin/cos parity are reflected through the sin_type and
!    cos_type variables.
!
      integer :: i

      do i = 1, mnmx
         if (sin_type .eq. 1 .and. cos_type .eq. 0) then
            anm(i) = m_fourier(i)*fnm(i)
         else if (sin_type .eq. 0 .and. cos_type .eq. 1) then
            anm(i) = -m_fourier(i)*fnm(i)
         end if
      end do
      if (sin_type .eq. 1 .and. cos_type .eq. 0) then
         sin_type = 0; cos_type = 1
      else if (sin_type .eq. 0 .and. cos_type .eq. 1) then
         sin_type = 1; cos_type = 0
      end if
   end subroutine dbydth
!
!
!
   subroutine dbydzt(sin_type, cos_type)
      use globals, only: mnmx, n_fourier, anm, fnm
      integer, intent(inout) :: sin_type, cos_type
!
!    Take the zeta derivative of the input Fourier amplitude array, fnm
!    and place the result in the output Fourier amplitude array, anm.
!    Changes to the sin/cos parity are reflected through the sin_type and
!    cos_type variables.
!
      integer :: i
      do i = 1, mnmx
         if (sin_type .eq. 1 .and. cos_type .eq. 0) then
            anm(i) = -n_fourier(i)*fnm(i)
         else if (sin_type .eq. 0 .and. cos_type .eq. 1) then
            anm(i) = n_fourier(i)*fnm(i)
         end if
      end do
      if (sin_type .eq. 1 .and. cos_type .eq. 0) then
         sin_type = 0; cos_type = 1
      else if (sin_type .eq. 0 .and. cos_type .eq. 1) then
         sin_type = 1; cos_type = 0
      end if
   end subroutine dbydzt

end module fourier_lib_convolve