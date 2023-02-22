module globals
!  Module for global variables

   use kind_spec

   implicit none

   !  PARAMETERS
   
   character(len=1), parameter :: jobz = 'V'
   logical, parameter :: subset_eq = .false.    !  Request a bracketed subset of eigenvalues. Not implemented?
   logical, parameter :: ipos_def_sym = .false. !  TODO???
   logical, parameter :: lrfp = .false.         ! If true: use reversed field pinch settings
   logical, parameter :: cyl = .false.          ! If true: cylindrical geometry???
   integer, parameter :: iopt = 1               !  Specifies problem to be solved by dggev, dsygv, dsygvx
   real(r8), parameter :: scale_khz = (1.d+3 * 2 * PI)**2


   !  Command line input
   integer :: irads, ir_fine_scl

   !  read_fourier_dat
   integer :: ith, izt, mpol, ntor, mnmx, ntors, nfp
   integer, allocatable, dimension(:) :: nw, mwl, mwu

   !  read_plasma_dat
   real(r8) :: aion, bion, cion, mass_ion, ion_density_0, nion(10), telec(10)
   integer :: ion_profile

   !  read_tae_data_boozer
   real(r8) :: bfavg

   !  Global variables

   integer :: mn_col
   integer, allocatable, dimension(:) :: im_col, in_col

   real(r8), allocatable, dimension(:) :: rho, rho_fine, rn, rm, rn_col, rm_col
   real(r8), allocatable, dimension(:) :: iotac, phipc, theta_tae, zeta_tae
   real(r8), allocatable, dimension(:) :: mu0_rho_ion, ion_density, iota_r, iota_r_inv
   
   !  Fourier variables
   real(r8), allocatable, dimension(:) :: fnm, f, anm
   real(r8), allocatable, dimension(:, :) :: cos_toF, sin_toF !cos_ar, sin_ar

   real(r8), allocatable, dimension(:,:,:) :: bfield, rjacob, gsssup
   real(r8), allocatable, dimension(:,:,:) :: bfield_lrg, rjacob_lrg, gsssup_lrg


end module globals
