module eig_solver

   use globals, only: ipos_def_sym, mn_col, subset_eq, iopt, jobz
   use kind_spec, only: r8

   implicit none

   private

   public :: initialize_solver, eigenvalue_solver
   public :: amat, bmat, alfr, vr, alfi, omega, beta
   ! public :: dggev, dsygv, dsygvx

   external :: dggev, dsygv, dsygvx ! LAPACK subroutines

   integer :: lwork, m_red
   integer, allocatable, dimension(:)    :: ifail, iwork
   
   real(r8), allocatable, dimension(:,:) :: amat, bmat, vr, vl
   real(r8), allocatable, dimension(:)   :: alfr, alfi, work, omega, beta


contains

   subroutine initialize_solver
      integer :: istat

      lwork = 20 * mn_col

      allocate (amat(mn_col, mn_col), stat = istat)
      allocate (bmat(mn_col, mn_col), stat = istat)
      allocate (alfr(mn_col), stat = istat)
      allocate (alfi(mn_col), stat = istat)
      allocate (work(lwork), stat = istat)
      allocate (ifail(mn_col), stat = istat)
      allocate (iwork(lwork), stat = istat)
      allocate (omega(mn_col), stat = istat)
      allocate (vl(mn_col, mn_col), stat = istat)
      allocate (vr(mn_col, mn_col), stat = istat)
      allocate (beta(mn_col), stat = istat)
 

   end subroutine initialize_solver

   subroutine eigenvalue_solver
      ! real(r8) :: beta
      real(r8) :: egl, egu, abstol
      integer :: il, iu, info
      egl = 1.d-2
      egu = 0.6d0
      abstol = 1.d-8
      il = 0
      iu = 0

      if (ipos_def_sym .eqv. .false.) then
         call dggev('N', 'V', mn_col, amat, mn_col, bmat, mn_col, alfr, alfi, beta, vl, mn_col, vr, mn_col, work, lwork, info)
      else
         if (subset_eq .eqv. .false.) then
            call dsygv(iopt, jobz, 'L', mn_col, amat, mn_col, bmat, mn_col, omega, work, lwork, info)
         else
            call dsygvx(iopt, 'V', 'V', 'L', mn_col, amat, mn_col, bmat, mn_col, egl, egu, il, iu, abstol, m_red, omega, vr, &
               mn_col, work, lwork, iwork, ifail, info)
         end if
      end if

      if (info .ne. 0) write (*, '("info = ",i8)') info
   end subroutine eigenvalue_solver


end module eig_solver
