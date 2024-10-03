module eig_solver

   use constants, only: r8

   implicit none

   private

   public :: initialize_solver, eigenvalue_solver
   public :: amat, bmat, alfr, vr, alfi, omega, beta

   external :: dggev, dsygv, dsygvx ! LAPACK subroutines

   integer :: lwork, m_red

   integer, allocatable, dimension(:)    :: ifail, iwork
   real(r8), allocatable, dimension(:,:) :: amat, bmat, vr, vl
   real(r8), allocatable, dimension(:)   :: alfr, alfi, work, omega, beta


contains

   subroutine initialize_solver
      use globals, only: mn_col
      integer :: istat

      lwork = 100 * mn_col

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
      use globals, only: ipos_def_sym, mn_col, subset_eq, iopt, jobz
      real(r8) :: egl, egu, abstol
      integer :: il, iu, info
      egl = 1.d-2
      egu = 0.6d0
      abstol = 1.d-8
      il = 0
      iu = 0

      if (ipos_def_sym .eqv. .false.) then
         ! Solve the generalized eigenvalue problem
         ! A * x = lambda * B * x
         ! where A and B are symmetric matrices and B is positive definite
         ! x is the right eigenvector, lambda is the eigenvalue
         call dggev('N', & ! No left eigenvectors
            'V', &      ! Compute right eigenvectors
            mn_col, &   ! Number of rows in matrices
            amat, &     ! Matrix A
            mn_col, &   ! Leading dimension of A (number of rows)
            bmat, &     ! Matrix B
            mn_col, &   ! Leading dimension of B (number of rows)
            alfr, &     ! Real part of eigenvalues
            alfi, &     ! Imaginary part of eigenvalues
            beta, &     ! Beta - eigenvalues are (alfr + i*alfi) / beta
            vl, &       ! Left eigenvectors - not computed
            mn_col, &   ! Leading dimension of VL (number of rows)
            vr, &       ! Right eigenvectors - computed, stored in columns
            mn_col, &   ! Leading dimension of VR (number of rows)
            work, &     ! Workspace array of size lwork
            lwork, &    ! Length of work array
            info, &     ! Output: 0 if successful
            )
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
