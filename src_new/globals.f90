module globals
!  Module for global variables

   use kind_spec

   implicit none

   integer :: irads, ir_fine_scl, irad3

   logical :: lrfp ! If true: use reversed field pinch settings

   !  read_fourier_dat
   integer :: ith, izt, mpol, ntor, nznt, mnmx, ntors, nfp, mode_family
   integer, allocatable :: nw(:), mwl(:), mwu(:)

   !  read_plasma_dat
   real(r8) :: aion, bion, cion, mass_ion, ion_density_0, &
      nion(10), telec(10)
   integer :: ion_profile

   !  read_tae_data_boozer
   real(r8) :: bfavg
   real(r8), allocatable, dimension(:) :: iotac, phipc, theta_tae, zeta_tae
   real(r8), allocatable, dimension(:,:,:) :: bfield, rjacob, gsssup
   real(r8), allocatable :: rho(:)



end module globals
