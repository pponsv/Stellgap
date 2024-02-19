module globals

   use constants

   implicit none

   !  Parameters

   real(r8), parameter :: viz_flux = 0.5   !selects surface for AVS data - sync with metric_element_create.f
   logical, parameter :: viz = .true.
   logical, parameter :: test_jacob = .false.
   logical, parameter :: test_upr_lowr = .false.
   logical, parameter :: make_stellgap_data = .true.
   logical, parameter :: make_full_torus = .true.
   logical, parameter :: surf_compute = .true.

   !  Variables
   integer :: itheta = 80
   integer :: izeta = 80

   !  Booz_xform parameters
   integer :: nfp, nsd, mnboz
   real(r8) :: amin, r0, beta0

   real(r8), dimension(:), allocatable :: hiota, hpres, hphip, hjtor, hjpol
   real(r8), dimension(:), ALLOCATABLE :: jprl_coef0, jprl_coef1, jprl_coef2
   real(r8), dimension(:), allocatable :: xm, xn
   real(r8), dimension(:, :), allocatable :: rmncbh, zmnsbh, pmnsbh, bmncbh

contains

end module globals
