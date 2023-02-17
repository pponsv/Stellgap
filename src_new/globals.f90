module globals
!  Module for global variables

   use kind_spec

   implicit none

   !    The subset_eq flag allows one to just request a bracketed subset of
   !     eigenvalues, rather than all eigenvalues. This option has not been
   !     developed beyond the call. Initial tests have not indicated it
   !     speeds things up.
   logical, parameter :: subset_eq = .false.

   logical, parameter :: ipos_def_sym = .false.
   logical :: lrfp ! If true: use reversed field pinch settings
   logical :: cyl  ! If true: cylindrical geometry???

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
   real(r8) :: scale_khz

   !  Global variables

   integer :: mn_col
   integer, allocatable, dimension(:) :: im_col, in_col

   real(r8), allocatable, dimension(:) :: rho, rho_fine, rn, rm, rn_col, rm_col
   real(r8), allocatable, dimension(:) :: iotac, phipc, theta_tae, zeta_tae
   real(r8), allocatable, dimension(:) :: mu0_rho_ion, ion_density, iota_r, iota_r_inv

   real(r8), allocatable, dimension(:,:,:) :: bfield, rjacob, gsssup
   real(r8), allocatable, dimension(:,:,:) :: bfield_lrg, rjacob_lrg, gsssup_lrg


end module globals
