module globals
!  Module for global variables

   use kind_spec

   implicit none

   logical :: lrfp ! If true: use reversed field pinch settings
   logical :: cyl  ! If true: cylindrical geometry???

   !  Command line input
   integer :: irads, ir_fine_scl

   !  read_fourier_dat
   integer :: ith, izt, mpol, ntor, nznt, mnmx, ntors, nfp, mode_family
   integer, allocatable, dimension(:) :: nw, mwl, mwu

   !  read_plasma_dat
   real(r8) :: aion, bion, cion, mass_ion, ion_density_0, nion(10), telec(10)
   integer :: ion_profile

   !  read_tae_data_boozer
   real(r8) :: bfavg
   real(r8) :: scale_khz

   real(r8), allocatable, dimension(:) :: rho, rho_fine
   real(r8), allocatable, dimension(:) :: iotac, phipc, theta_tae, zeta_tae
   real(r8), allocatable, dimension(:) :: mu0_rho_ion, ion_density, iota_r, iota_r_inv

   real(r8), allocatable, dimension(:,:,:) :: bfield, rjacob, gsssup
   real(r8), allocatable, dimension(:,:,:) :: bfield_lrg, rjacob_lrg, gsssup_lrg


end module globals
