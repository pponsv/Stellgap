module read_input

   use constants
   use read_boozer_mod

   implicit none

contains

   subroutine read_args_and_boozer
      use globals, only : nfp, nsd, amin, r0, beta0, mnboz
      use globals, only : hiota, hpres, hjpol, hjtor, hphip
      use globals, only : xm, xn
      use globals, only : rmncbh, zmnsbh, pmnsbh, bmncbh
      integer :: ierr
      character(len=128) :: boozer_extension

      ! Check the number of command line arguments
      if (iargc() .ne. 1) then
         print *, ' MUST ENTER FILE SUFFIX ON COMMAND LINE'
         stop
      end if

      ! Read the file suffix from the command line
      call getarg(1, boozer_extension)

      ! Read the boozer file
      call read_boozer_file(boozer_extension, ierr)

      ! Check for errors
      if (ierr .ne. 0) then
         print *, ' ERROR READING BOOZER FILE'
         stop
      end if

      nfp = nfp_b ! Nunber of field periods
      nsd = ns_b ! Number of radial surfaces
      r0 = (rmax_b + rmin_b) / 2.0_r8 ! Mean radius?
      amin = r0 / aspect_b
      beta0 = betaxis_b
      mnboz = mnboz_b

      call allocate_nsd_mnboz
!
!...   IOTA, PRES, PHIP (= -PHIP_VMEC), JTOR (=I = -BUCO_VMEC) and
!         JPOL (=J = BVCO_VMEC) are ALL on HALF-MESH!!!

      hiota = iota_b
      hpres = mu_0 * pres_b  ! in VMEC versions > 6.00 pressure is given in pascals
      hphip = -phip_b        ! toroidal fluxes have REVERSED sign respect to VMEC!!
      hjpol = bvco_b
      hjtor = -buco_b        ! toroidal fluxes have REVERSED sign respect to VMEC!!

!       if (beta0 .le. ZERO) stop 'Beta0 <= 0'
!       b0 = sqrt((TWO/beta0)*(1.5_dp*hpres(2)-.5_dp*hpres(3)))

      xm = ixm_b
      xn = ixn_b

      ! RMN, ZMN, PMN and BMN are ALL on HALF-MESH
      bmncbh = bmnc_b
      rmncbh = rmnc_b
      zmnsbh = zmns_b
      pmnsbh = pmns_b

      call read_boozer_deallocate

   end subroutine read_args_and_boozer

   subroutine allocate_nsd_mnboz
      use globals, only : nsd, mnboz
      use globals, only : hiota, hpres, hjpol, hjtor, hphip, jprl_coef0, &
         jprl_coef1, jprl_coef2
      use globals, only: rmncbh, zmnsbh, pmnsbh, bmncbh
      use globals, only : xm, xn

      ! allocate (hiota(nsd))
      ! allocate (hpres(nsd))
      ! allocate (hjpol(nsd))
      ! allocate (hjtor(nsd))
      ! allocate (hphip(nsd))
      allocate (jprl_coef0(nsd))
      allocate (jprl_coef1(nsd))
      allocate (jprl_coef2(nsd))

      ! allocate (xm(mnboz))
      ! allocate (xn(mnboz))

      ! allocate (rmncbh(mnboz, nsd))
      ! allocate (zmnsbh(mnboz, nsd))
      ! allocate (pmnsbh(mnboz, nsd))
      ! allocate (bmncbh(mnboz, nsd))

   end subroutine allocate_nsd_mnboz

end module read_input
