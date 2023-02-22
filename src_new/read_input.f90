module input

   use kind_spec

   implicit none

   public

contains

   subroutine read_tae_data_boozer
      use globals, only: izt, ith, irads, lrfp, rho, rho_fine, iotac, theta_tae, zeta_tae,&
         rjacob, bfield, gsssup, bfavg, ir_fine_scl

      integer :: ir, i, j, nn(irads), istat
      real(r8) :: dum1, dum2, dm1, dm2, dm3, dm4, dm5
      real(r8) :: phipc(irads)

      character(len=*), parameter :: fmt_irad = '(1x,i3,4(2x,e15.7))'
      character(len=*), parameter :: fmt_ithz = '(1x,4(e24.12,2x),e24.12)'

      allocate (bfield(izt, ith, irads), stat = istat)
      allocate (gsssup(izt, ith, irads), stat = istat)
      allocate (rjacob(izt, ith, irads), stat = istat)
      allocate (theta_tae(ith), stat = istat)
      allocate (zeta_tae(izt), stat = istat)

      open (unit=20, file="tae_data_boozer", status="old")

      select case (lrfp)
       case (.true.)
         write (*, '("Using RFP settings")')
       case (.false.)
         write (*, '("Using tokamak/stellarator settings")')
      end select

      do ir = 1, irads
         read (20, fmt_irad) nn(ir), iotac(ir), phipc(ir), dum1, dum2
         do i = 1, izt
            do j = 1, ith
               read (20, fmt_ithz) theta_tae(j), zeta_tae(i), bfield(i, j, ir), gsssup(i, j, ir), dm1
               read (20, fmt_ithz) dm2, dm3, dm4, dm5, rjacob(i, j, ir)
            end do
         end do
      end do

      if (.not. lrfp) then
         do ir=1, irads
            rjacob(:, :, ir) = rjacob(:, :, ir)/phipc(ir)
         end do
      end if

      rho = real(nn, kind=r8) / real(nn(irads), kind=r8)
      rho_fine = rho(1) + (rho(irads) - rho(1)) * [(i, i=0,ir_fine_scl-1)] / (ir_fine_scl - 1)
      bfavg = sum(bfield) / real(irads*ith*izt, kind=r8)

   end subroutine read_tae_data_boozer


   subroutine read_args
      use globals, only : irads, ir_fine_scl

      character*10 :: arg1, arg2

      irads = 41
      ir_fine_scl = 128  !default values

      call getarg(1, arg1)
      call getarg(2, arg2)

      read (arg1, '(i3)') irads
      read (arg2, '(i4)') ir_fine_scl

      ! write (*, '("irads = ",i4," irads3 = ",i4," ir_fine_scl = ",i5)') irads, 3*irads, ir_fine_scl

   end subroutine read_args


   subroutine read_plasma_dat
      use globals, only: mass_ion, ion_density_0, ion_profile, nion, telec, aion, bion, cion

      real(r8) :: ion_to_proton_mass
      logical :: jdqz_data       ! Useless, for compatibility with AE3D
      character*4 :: egnout_form ! Useless, for compatibility with AE3D

      namelist /plasma_input/ ion_to_proton_mass, ion_density_0, ion_profile, jdqz_data, &
         egnout_form, nion, telec, aion, bion, cion

      !  Initialize to avoid surprises
      nion = 0.
      telec = 0.

      open (unit=4, file="plasma.dat", status="old")
      read (4, plasma_input)
      close(4)

      mass_ion = mass_proton * ion_to_proton_mass

   end subroutine read_plasma_dat


   subroutine read_fourier_dat
      use globals, only: nfp, ith, izt, ntors, nw, mwl, mwu, mpol, ntor, mnmx
      integer :: i, istat, mode_family

      open (unit=20, file="fourier.dat", status="old")
      read (20, *) nfp, ith, izt, mode_family
      read (20, *) ntors

      allocate (nw(ntors), stat=istat)
      allocate (mwl(ntors), stat=istat)
      allocate (mwu(ntors), stat=istat)
      
      do i = 1, ntors
         read (20, *) nw(i), mwl(i), mwu(i)
      end do
      close (unit=20)

      mpol = ith*2/5
      ntor = izt*2/5 ! Done to avoid aliasing issues

      ! nznt = ith*izt
      mnmx = (2*ntor + 1)*mpol - ntor

   end subroutine read_fourier_dat



end module input
