module input

   use kind_spec

   implicit none

   !  read_args
   character*10 :: arg1, arg2
   integer :: irads, ir_fine_scl, irad3

   !  read_plasma_dat
   real(r8) :: ion_to_proton_mass, ion_density_0, aion, bion, cion, &
      nion(10), telec(10)
   integer :: ion_profile

   !  read_fourier_dat
   integer :: ith, izt, mpol, ntor, nznt, mnmx, ntors, nfp, mode_family
   integer, allocatable :: nw(:), mwl(:), mwu(:)

   !  read_tae_data_boozer
   logical :: lrfp
   integer :: nn
   real(r8) :: bfavg
   real(r8), allocatable, dimension(:) :: iotac, phipc, iotac_inv, nsurf, b2, &
      theta_tae, zeta_tae
   real(r8), allocatable, dimension(:,:,:) :: bfield, rjacob, gsssup


contains

   subroutine read_tae_data_boozer
      integer :: ir, i, j
      real(r8) :: dum1, dum2
      real(r8) :: dm1, dm2, dm3, dm4, dm5

      open (unit=20, file="tae_data_boozer", status="old")

      if (lrfp) write (*, '("Using RFP settings")')
      if (.not. lrfp) write (*, '("Using tokamak/stellarator settings")')

      do ir = 1, irads
         read (20, '(1x,i3,4(2x,e15.7))') nn, iotac(ir), phipc(ir), dum1, dum2
         iotac_inv(ir) = 1.d0/iotac(ir)
         nsurf(ir) = dble(nn)
         b2(ir) = 0.d0
         do i = 1, izt
            do j = 1, ith
               read (20, '(1x,4(e24.12,2x),e24.12)') theta_tae(j), zeta_tae(i), &
                  bfield(i, j, ir), gsssup(i, j, ir), dm1
               b2(ir) = b2(ir) + bfield(i, j, ir)/(dble(izt*ith))
               read (20, '(1x,4(e24.12,2x),e24.12)') dm2,&
                  dm3, dm4, dm5, rjacob(i, j, ir)

               if (.not. lrfp) then
                  rjacob(i, j, ir) = rjacob(i, j, ir)/phipc(ir)
               end if
            end do
         end do
         !       write(*,*) ir, bfield(izt,ith,ir), rjacob(izt,ith,ir)
      end do

      bfavg = sum(b2)/(dble(irads))

   end subroutine read_tae_data_boozer


   subroutine read_args
      irads = 41; ir_fine_scl = 128  !default values
      call getarg(1, arg1)
      call getarg(2, arg2)

      read (arg1, '(i3)') irads
      read (arg2, '(i4)') ir_fine_scl
      irad3 = 3*irads

      write (*, '("irads = ",i4," irads3 = ",i4," ir_fine_scl = ",i5)') irads, irad3, ir_fine_scl

   end subroutine read_args


   subroutine read_plasma_dat
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

   end subroutine read_plasma_dat


   subroutine read_fourier_dat
      !  Dummy
      integer :: i, istat

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
      nznt = ith*izt
      mnmx = (2*ntor + 1)*mpol - ntor

   end subroutine read_fourier_dat



end module input
