!                              D   I   S   C   L   A   I   M   E   R
!
!       You are using a BETA version of the program stellgap.f, which is currently
!       under development by D. A. Spong of the Fusion Energy Division,
!       Oak Ridge National Laboratory.  Please report any problems or comments
!       to him.  As a BETA version, this program is subject to periodic change
!       and improvement without notice.
!
!
!  Note: this code uses normalized toroidal flux (variables: rho, r_pt) as
!        an independent variable. If a radius-like variable is desired
!        (i.e., sqrt(tor_flux) or radius based on sqrt(area/pi)), this must be
!        constructed in subsequent codes that use the AE3D eigenfunction data
!
!   Before running, check that the irads variable is set consistently with the
!   number of surfaces in the tae_boozer_data file. irads should be the number
!   of surfaces in the original VMEC run minus 1 or 2 to avoid edge problems.
!
!  Modified 2/23/2013 to properly treat RFPs where there is a reversal surface
!    (q -> 0, iota -> infinity). This is handled through options controlled by
!    the logical variable lrfp. lrfp is read from the first line of tae_data_boozer,
!    which is written by xmetric. For analysis of systems with lrfp = .true. both
!    updated versions of metric_element_create.f and stellagap.f are required.
!    If using an older tae_data_boozer file, a read error will result since the
!    first line will not have the lrfp logical. One can edit tae_data_boozer
!    and add "F" in the first line to resolve this. However, to treat cases with
!    lrfp = .true. it is necessary to have also used the upgraded xmetric. By
!    default lrfp is set to .false.
!
!  Modified 7/6/2012 making ir_fine_scl and irads to be command-line arguments
!   (if not given, default values of 128 and 41 are used)
!
!  Modified 12/2008 so that ion density profile, ion mass and several other
!   parameters are read in through the plasma.dat file rather than being
!   hardwired into the code. (see explanation below for variables contained
!   in plasma.dat file).
!
!  Modified on 9/18/2008 to offer both parallel or serial options, depending on
!   use of the precompiler.
!
!  This is the stellgap code for calculating Alfven contina in stellarators. The
!  associated publication that describes this calculation is:
!
!   D. A. Spong, R. Sanchez, A. Weller, "Shear Alfven Continua in Stellarators,"
!     Phys. of Plasmas, Vol. 10, pg. 3217, August, 2003.
!
!
!  This version of stellgap.f has been modified to be compatible with ae_solve.f.
!  The frequencies are output in kHz. The mode selection is done using the fourier.dat
!  file. The equilibrium data comes from the tae_data_boozer file. Following the
!  stellgap run, a post-processing code, post_process.f must be run. The gap plot
!  is then made using the results in the alfven_post file.
!  The central ion density, ion density profile, and ion mass are set from the
!  plasma.dat input file.
!  Two other variables that the user should set are irads and ir_fine_scl. These
!  are set througth paramter statements in the beginning of program tae_continua:
!      integer, parameter :: ir_fine_scl = 1000
!      integer, parameter :: irads = 39, irad3 = 3*irads
!  irads = # of flux surfaces in the initial equilibrium (wout file) - 2.
!  The outermost and innermost surfaces from the initial equilibrium
!  are left off to avoid inaccuracies that are sometimes present near
!  these regions.
!  ir_fine_scl = # of surfaces desired in the continuum output. Generally, to
!  resolve fine scale features in the Alfven continua, more surfaces are
!  desireable than in the original equilibrium data. The code performs
!  interpolations to allow this.
!
!
!  Input parameters:   (fourier.dat file).
!
!   ith, izt = number of theta and zeta grid points in the original mesh data
!     (in file tae_data_vmec) calculated by the xmetric code.
!    This file contains data
!    for iota, phipc, |B|, gssup (the psi-psi contra-variant metric element,
!    the Jacobian, and the contravariant poloidal and toroidal components of
!    the magnetic field)
!   nfp = number of field periods (for a tokamak run, set nfp = 1)
!   mpol, ntor = number of poloidal and toroidal modes used in the Fourier
!    expansion of the above gridded quantities suplied by the tae_data_vmec
!    data file.  To avoid anti-aliasing errors, mpol < 0.5*ith and
!    ntor < 0.5*izt.
!   mp_col = number of m,n mode pairs used in the representation of the
!    eigebfunction
!   nt_col = number of toroidal modes used in the representation of the
!    eigenfunction (should be an integral multiple of nfp)
!   mode_family = toroidal mode number about which the toroidal eigenfunction
!     spectrum is built. i.e., for m = 0,
!     n = mode_family, mode_family + nfp, mode_family + 2*nfp,
!         ..., mode_family + nt_col
!      while for m .ne. 0,
!      n = -nt_col + mode_family, -nt_col + nfp + mode_family, -nt_col
!            + 2*nfp + mode_family, ..., nt_col + mode_family
!
!   Input parameters:   (plasma.dat file).
!
!   ion_to_proton_mass = m_ion/m_proton
!   ion_density_0 = n_ion(0) = ion density (in m**-3) at magnetic axis
!   ion_profile = integer parameter that determines form of n_ion(r)/n_ion(0) fit
!      for ion_profile = 0 ion_density = [iota(rho)/iota(0)]**2
!      for ion_profile = 1 ion_density = polynomial fit
!                          = nion(1) + nion(2)*rho + nion(3)*(rho**2)
!                           + nion(4)*(rho**3) + ... + nion(9)*(rho**8) + nion(10)*(rho**9)
!               (here rho = normalized toroidal flux)
!      for ion_profile = 2 ion_density = ion_density = constant = 1.
!      for ion_profile = 3 ion_density = [1 - aion*(rho**bion)]**cion
!
!   nion(10) = array of polynomial fit coefficients for ion_profile = 1
!   aion, bion, cion = parameters for ion_profile = 3 fit
!   jdqz_data = logical variable that is used in ae3d, but not in stellgap.f
!                (included here only so that same plasma.dat file can be used)
!   egnout_form = "binr" or "asci" - like jdqz_data this variable is used in
!                  ae3d, but not stellgap.f
!
! compile: sh bld_stellgap
!
!  This code is can be run in serial or parallel mode, depending on whether
!   is it run through the precompiler with SERIAL or PARALLEL set.
!
!

!-----------------------------------------------------------------
program tae_continua
   use fourier_lib
   use Fourier_lib_convolve
   use kind_spec
   use postprocess
   use fitpack
   implicit none
#if defined (PARALLEL)
   include 'mpif.h'
#endif
!
!      This code is a test of a Fourier xform of a 2D function
!      (f) vs. theta, zeta.  f has cos/sin stellarator-like symmetry.
!      The Fourier xform is done by subroutine toFourier using
!      cos/sin expansions. (use with Fourier_lib.f)
!
!
   external dggev, dsygv, dsygvx
   integer :: ir_fine_scl, irads, irad3
   integer, parameter :: iopt = 1
   real(kind=rprec), parameter :: R0 = 1.0
   real(kind=rprec), parameter :: mass_proton = 1.67d-27
   real(kind=rprec) :: ion_density_0, mass_ion, ion_to_proton_mass,&
   &aion, bion, cion
!     Telec does nothing, added for compat. with sound version
   real(kind=rprec) :: nion(10), telec(10)
   integer :: ir, irr, il, iu, ion_profile
   character*4 :: egnout_form
   logical :: jdqz_data
   real(kind=rprec), dimension(:, :, :), allocatable :: bfield,&
   &gsssup, rjacob, bsupth, bsupzt, bfield_lrg,&
   &gsssup_lrg, bsupth_lrg, bsupzt_lrg, rjacob_lrg
   real(kind=rprec), dimension(:), allocatable :: theta_tae,&
   &zeta_tae, f1_nm, f3a_nm, f3b_nm, f3c_nm, bsupth_tmp, bsupzt_tmp
   real(kind=rprec), dimension(:, :), allocatable :: f1, f3a,&
   &f3b, f3c

   real(kind=rprec), dimension(:), allocatable :: bavg, mu0_rho_ion,&
   &ion_density, iota_r, iota_r_inv

   real(kind=rprec), allocatable :: beta(:), aux(:), eig_vect(:)
   real(kind=rprec), allocatable :: amat(:, :),&
   &bmat(:, :), z_r(:, :), vr(:, :), vl(:, :)
   real(kind=rprec), allocatable :: bgrad2(:, :, :)
   real(kind=rprec), allocatable :: omega(:), work(:), alfr(:),&
   &alfi(:)
   complex*16, allocatable :: alpha(:)

   real(kind=rprec), dimension(:), allocatable :: nsurf, iotac,&
   &phipc, rho, sp_fit, b2, iotac_inv

   real(kind=rprec) :: sp1, sp2, sigma_spl, r_pt,&
   &iotac_r, dum1, dum2, dum3, dum4

   real(kind=rprec), dimension(:), allocatable :: yp, temp,&
   &ypi, tempi

   real(kind=rprec) :: amat1, amat2, amat3, amat4, amat5
   real(kind=rprec) :: f1_coef, f2a_coef, f3a_coef,&
   &f3b_coef, f3c_coef
   real(kind=rprec) :: dm1, dm2, dm3, dm4, dm5, bfavg, va
   real(kind=rprec) :: f1max, eig_max, mu0, scale_khz,&
   &f3amax, f3bmax, f3cmax, ccci, scsi, test,&
   &egl, egu, abstol, f1_avg, f3_avg
   real(kind=rprec) :: ra, ra2, ra3, ra4, ra5, ra6
   integer :: ispl_opt, ierr_spl, num_eq, naux
   integer :: ier, npes, mype, numrads
   integer :: ni, nj, mi, mj, ieq, meq, neq, ii, jj
   integer :: j_max_index, ios
   integer :: m_emax, n_emax, isym_opt, in, im
   integer :: ldvl, ldvr, lwork, info, nn, m_red
   integer, allocatable :: ifail(:), iwork(:)
   logical :: cyl, lrfp
   character*20 outfile
   character*3 procnum
   character*1 jobz
   character*10 date, time, zone
   character*10 arg1, arg2
   integer values(8)
   namelist /plasma_input/ ion_to_proton_mass, ion_density_0,&
   &ion_profile, jdqz_data, egnout_form, nion, telec,&
   &aion, bion, cion
   cyl = .false.; lrfp = .false.
   call getarg(1, arg1)
   call getarg(2, arg2)
   irads = 41; ir_fine_scl = 128  !default values
   read (arg1, '(i3)') irads
   read (arg2, '(i4)') ir_fine_scl
   irad3 = 3*irads

   write (*, '("irads = ",i4," irads3 = ",i4," ir_fine_scl = ",i5)')&
   &irads, irad3, ir_fine_scl

   allocate (bavg(ir_fine_scl), mu0_rho_ion(ir_fine_scl),&
   &ion_density(ir_fine_scl), iota_r(ir_fine_scl),&
   &iota_r_inv(ir_fine_scl), stat=istat)
   allocate (nsurf(irads), iotac(irads), iotac_inv(irads),&
   &phipc(irads), rho(irads), sp_fit(irads),&
   &b2(irads), stat=istat)
   allocate (yp(irad3), temp(irad3), ypi(irad3),&
   &tempi(irad3), stat=istat)

   nion = 0.   ! Initialize to avoid random values if input too short
   open (unit=4, file="plasma.dat", status="old")
   read (4, plasma_input)
   print *, "NION", nion
   close (unit=4)
   mass_ion = mass_proton*ion_to_proton_mass
#if defined (PARALLEL)
! Get NPES and MYPE.  Requires initialization of MPI.
   call MPI_INIT(ier)
   if (ier .ne. 0) then
      write (*, '("mpi_init returned ier =")') ier
      stop
   end if
   call MPI_COMM_SIZE(MPI_COMM_WORLD, npes, ier)
   if (ier .ne. 0) then
      write (*, '("mpi_comm_size returned ier =")') ier
      stop
   end if
   call MPI_COMM_RANK(MPI_COMM_WORLD, mype, ier)
   if (ier .ne. 0) then
      write (*, '("mpi_comm_rank returned ier =")') ier
      stop
   end if
#endif
#if defined (SERIAL)
   mype = 0
   npes = 1
#endif
!
   jobz = 'V'
!
#if defined (PARALLEL)
   write (procnum, '(i3)') mype
   write (*, *) procnum
   outfile = "alfven_spec"//trim(adjustl(procnum))
#endif
#if defined (SERIAL)
   outfile = "alfven_spec"
#endif
   if (mype .eq. 0) then
      write (*, *) mype, npes
      write (*, *) procnum
      write (*, *) outfile
      call date_and_time(date, time, zone, values)
      write (*, *) time
   end if
   if (ipos_def_sym) isym_opt = 1
   if (.NOT. ipos_def_sym) isym_opt = 0
!
!    Generate Fourier arrays
!
   call readin
   call trig_array
   call convolution_array
   num_eq = mnmx
   naux = 10*mn_col
   mu0 = 2.d-7*twopi
   scale_khz = (1.d+3*twopi)**2
   allocate (aux(naux), stat=istat)
   allocate (alpha(mn_col), stat=istat)
   allocate (beta(mn_col), stat=istat)
   allocate (eig_vect(mn_col), stat=istat)
   allocate (z_r(mn_col, mn_col), stat=istat)
   allocate (omega(mn_col), stat=istat)
   allocate (amat(mn_col, mn_col), stat=istat)
   allocate (bmat(mn_col, mn_col), stat=istat)
   allocate (ifail(mn_col), stat=istat)
   ldvl = mn_col; ldvr = mn_col; lwork = 20*mn_col
   allocate (work(lwork), stat=istat)
   allocate (iwork(lwork), stat=istat)
   allocate (alfr(mn_col), stat=istat)
   allocate (alfi(mn_col), stat=istat)
   allocate (vl(mn_col, mn_col), stat=istat)
   allocate (vr(mn_col, mn_col), stat=istat)
   if (mype .eq. 0)&
   &write (*, *) izt, ith, irads, mnmx, nznt, mn_col
   allocate (bfield(izt, ith, irads), stat=istat)
   allocate (gsssup(izt, ith, irads), stat=istat)
   allocate (rjacob(izt, ith, irads), stat=istat)
   allocate (bsupth(izt, ith, irads), stat=istat)
   allocate (bsupzt(izt, ith, irads), stat=istat)
   allocate (bfield_lrg(izt, ith, ir_fine_scl), stat=istat)
   allocate (gsssup_lrg(izt, ith, ir_fine_scl), stat=istat)
   allocate (bsupth_lrg(izt, ith, ir_fine_scl), stat=istat)
   allocate (bsupzt_lrg(izt, ith, ir_fine_scl), stat=istat)
   allocate (rjacob_lrg(izt, ith, ir_fine_scl), stat=istat)
   allocate (theta_tae(ith), stat=istat)
   allocate (zeta_tae(izt), stat=istat)
   allocate (f1_nm(mnmx), stat=istat)
   allocate (f3a_nm(mnmx), stat=istat)
   allocate (f3b_nm(mnmx), stat=istat)
   allocate (f3c_nm(mnmx), stat=istat)
   allocate (bsupth_tmp(nznt), stat=istat)
   allocate (bsupzt_tmp(nznt), stat=istat)
   allocate (f1(izt, ith), stat=istat)
   allocate (f3a(izt, ith), stat=istat)
   allocate (f3b(izt, ith), stat=istat)
   allocate (f3c(izt, ith), stat=istat)
!
!   Open files for output
!
!      write(*,*) trim(adjustl(outfile))

   open (unit=21, file=trim(adjustl(outfile)), status="unknown")
   open (unit=8, file="coef_arrays", status="unknown")
   if (mype .eq. 0) then
      open (unit=7, file="data_post", status="unknown")
      open (unit=9, file="ion_profile", status="unknown")
   end if
!
!    Boozer coordinates input - new ae-mode-structure input
!
   open (unit=20, file="tae_data_boozer", status="old")
!      read(20,'(L)') lrfp    !uncomment and remove next line if lrfp added to tae_data_boozer
   lrfp = .false.
   if (lrfp) write (*, '("Using RFP settings")')
   if (.not. lrfp) write (*, '("Using tokamak/stellarator settings")')
   do ir = 1, irads
      read (20, '(1x,i3,4(2x,e15.7))') nn,&
      &iotac(ir), phipc(ir), dum1, dum2
      iotac_inv(ir) = 1.d0/iotac(ir)
      nsurf(ir) = dble(nn)
      b2(ir) = 0.d0
      do i = 1, izt
         do j = 1, ith
            read (20, '(1x,4(e24.12,2x),e24.12)') theta_tae(j),&
            &zeta_tae(i),&
            &bfield(i, j, ir), gsssup(i, j, ir),&
            &dm1
            b2(ir) = b2(ir) + bfield(i, j, ir)/(dble(izt*ith))
            read (20, '(1x,4(e24.12,2x),e24.12)') dm2,&
            &dm3, dm4, dm5,&
            &rjacob(i, j, ir)

            if (lrfp) then
               rjacob(i, j, ir) = rjacob(i, j, ir)
            else
               rjacob(i, j, ir) = rjacob(i, j, ir)/phipc(ir)
            end if
         end do
      end do
!       write(*,*) ir, bfield(izt,ith,ir), rjacob(izt,ith,ir)
   end do
   bfavg = sum(b2)/(dble(irads))
!
!    Record equilibrium and eigenfunction Fourier mode list
!
   if (mype .eq. 0) then
      open (unit=22, file="modes", status="unknown")
      write (22, '("Equilibrium modes:",/)')
      write (22, '("meq    neq")')
      do i = 1, mnmx
         meq = rm(i)
         neq = rn(i)
         write (22, '(i3,3x,i3)') meq, neq
      end do
      write (22, '(///,"Eigenvector modes:",/)')
      write (22, '("m    n")')
      do i = 1, mn_col
         write (22, '(i3,3x,i3)') im_col(i), in_col(i)
      end do
      close (unit=22)
   end if
!
!   Make spline fits and fill in fine_radius_scale arrays
!
   do ir = 1, irads
      rho(ir) = nsurf(ir)/nsurf(irads)
#if defined (SERIAL)
!         write(*,*) rho(ir)
#endif
   end do
   ispl_opt = 3; ierr_spl = 0; sigma_spl = 0.
   call curv1(irads, rho, iotac, sp1, sp2, ispl_opt, ypi,&
   &tempi, sigma_spl, ierr_spl)
   do irr = 1, ir_fine_scl
      r_pt = rho(1) + real(irr - 1)*(rho(irads) - rho(1))&
      &/real(ir_fine_scl - 1)
      ra = sqrt(r_pt); ra2 = ra**2; ra3 = ra**3; ra4 = ra**4
      ra5 = ra**5; ra6 = ra**6
      iota_r(irr) = curv2(r_pt, irads, rho, iotac, ypi, sigma_spl)

      if (ion_profile .eq. 0) then
         ion_density(irr) = (iota_r(irr)/iotac(1))**2   !profile that lines up gaps
      else if (ion_profile .eq. 1) then
         ion_density(irr) = nion(1) + r_pt*nion(2) + nion(3)*(r_pt**2)&
         &+ nion(4)*(r_pt**3) + nion(5)*(r_pt**4) + nion(6)*(r_pt**5)&
         &+ nion(7)*(r_pt**6) + nion(8)*(r_pt**7) + nion(9)*(r_pt**8)&
         &+ nion(10)*(r_pt**9)
      else if (ion_profile .eq. 2) then
         ion_density(irr) = 1.d0
      else if (ion_profile .eq. 3) then
         ion_density(irr) = (1.-aion*(r_pt**bion))**cion
      end if

      mu0_rho_ion(irr) = mu0*mass_ion*ion_density_0&
                        &*ion_density(irr)*scale_khz
      va = sqrt(bfavg**2/(mu0_rho_ion(irr)/scale_khz))
      if (mype .eq. 0)&
      &write (9, 67) r_pt, ion_density_0*ion_density(irr),&
      &iota_r(irr), va
   end do     !irr = 1,ir_fine_scl
!   Interpolate 1/iota for lrfp = true option
   ispl_opt = 3; ierr_spl = 0; sigma_spl = 0.
   call curv1(irads, rho, iotac_inv, sp1, sp2, ispl_opt, ypi,&
   &tempi, sigma_spl, ierr_spl)
   do irr = 1, ir_fine_scl
      r_pt = rho(1) + real(irr - 1)*(rho(irads) - rho(1))&
      &/real(ir_fine_scl - 1)
      iota_r_inv(irr) = curv2(r_pt, irads, rho, iotac_inv, ypi, sigma_spl)
   end do     !irr = 1,ir_fine_scl

   if (ierr_spl .ne. 0) write (*, '("spline error 1",i3)') ierr_spl
67 format(e12.5, 3(3x, e12.5))
   do i = 1, izt
      do j = 1, ith
!
!       Bfield array
!
         do ir = 1, irads
            sp_fit(ir) = bfield(i, j, ir)
         end do      !ir=1,irads
         ispl_opt = 3; ierr_spl = 0; sigma_spl = 0
         call curv1(irads, rho, sp_fit, sp1, sp2, ispl_opt, yp,&
         &temp, sigma_spl, ierr_spl)
         if (ierr_spl .ne. 0) write (*, '("spline error 2",i3)') ierr_spl
         do irr = 1, ir_fine_scl
            r_pt = rho(1) + real(irr - 1)*(rho(irads) - rho(1))&
            &/real(ir_fine_scl - 1)
            bfield_lrg(i, j, irr) = curv2(r_pt, irads, rho, sp_fit, yp, sigma_spl)
         end do     !irr = 1,ir_fine_scl
!
!
!       Jacobian array
!
         do ir = 1, irads
            sp_fit(ir) = rjacob(i, j, ir)
         end do      !ir=1,irads
         ispl_opt = 3; ierr_spl = 0; sigma_spl = 0
         call curv1(irads, rho, sp_fit, sp1, sp2, ispl_opt, yp,&
         &temp, sigma_spl, ierr_spl)
         if (ierr_spl .ne. 0) write (*, '("spline error 2",i3)') ierr_spl
         do irr = 1, ir_fine_scl
            r_pt = rho(1) + real(irr - 1)*(rho(irads) - rho(1))&
            &/real(ir_fine_scl - 1)
            rjacob_lrg(i, j, irr) = curv2(r_pt, irads, rho, sp_fit, yp, sigma_spl)
         end do     !irr = 1,ir_fine_scl
!
!       Gsssup array
!
         do ir = 1, irads
            sp_fit(ir) = gsssup(i, j, ir)
         end do
         ispl_opt = 3; ierr_spl = 0; sigma_spl = 0
         call curv1(irads, rho, sp_fit, sp1, sp2, ispl_opt, yp,&
         &temp, sigma_spl, ierr_spl)
         if (ierr_spl .ne. 0) write (*, '("spline error 3",i3)') ierr_spl
         do irr = 1, ir_fine_scl
            r_pt = rho(1) + real(irr - 1)*(rho(irads) - rho(1))&
            &/real(ir_fine_scl - 1)
            gsssup_lrg(i, j, irr) = curv2(r_pt, irads, rho, sp_fit, yp, sigma_spl)
         end do     !irr = 1,ir_fine_scl
      end do     !j=1,ith
   end do     !i=1,izt
   close (unit=9)
!
!
!     Fine_scale arrays and spline fits finished
!
   numrads = ir_fine_scl/npes
!       write(*,*) numrads
   do ir = (mype*numrads + 1), (mype*numrads + numrads)
!       write(*,*) ir
!      do ir=1,ir_fine_scl
      r_pt = rho(1) + real(ir - 1)*(rho(irads) - rho(1))&
      &/real(ir_fine_scl - 1)
      f1_avg = 0.; f3_avg = 0.
      do i = 1, izt
         do j = 1, ith
            f1(i, j) = gsssup_lrg(i, j, ir)*rjacob_lrg(i, j, ir)/&
            &(bfield_lrg(i, j, ir)**2)
            f1_avg = f1_avg + f1(i, j)/real(izt*ith)
            if (.not. lrfp) then
               f3c(i, j) = gsssup_lrg(i, j, ir)/&
               &(rjacob_lrg(i, j, ir)*(bfield_lrg(i, j, ir)**2))
               f3_avg = f3_avg + f3c(i, j)/real(izt*ith)
               f3b(i, j) = iota_r(ir)*f3c(i, j)
               f3a(i, j) = iota_r(ir)*f3b(i, j)
            else if (lrfp) then
               f3a(i, j) = gsssup_lrg(i, j, ir)/&
               &(rjacob_lrg(i, j, ir)*(bfield_lrg(i, j, ir)**2))
               f3b(i, j) = f3a(i, j)*iota_r_inv(ir)
               f3c(i, j) = f3b(i, j)*iota_r_inv(ir)
               f3_avg = f3_avg + f3a(i, j)/real(izt*ith)
            end if
         end do
      end do

      if (cyl) then
         if (.not. lrfp) then
            do i = 1, izt
               do j = 1, ith
                  f1(i, j) = f1_avg
!         f3a(i,j) = f3_avg*iota_r(ir)*iota_r(ir)
!         f3b(i,j) = f3_avg*iota_r(ir)
!         f3c(i,j) = f3_avg
               end do
            end do
         else if (lrfp) then
            do i = 1, izt
               do j = 1, ith
                  f1(i, j) = f1_avg
!         f3a(i,j) = f3_avg
!         f3b(i,j) = f3_avg*iota_r_inv(ir)
!         f3c(i,j) = f3_avg*iota_r_inv(ir)**iota_r_inv(ir)
               end do
            end do
         end if
      end if
!
!
!     Generate Fourier spectra of above gridded coefficients
!
      lg = 0
      do i = 1, izt
         do j = 1, ith
            lg = lg + 1
            f(lg) = f1(i, j)
         end do
      end do
      sin_type = 0; cos_type = 1
      call toFourier
      f1_nm(:) = fnm(:)
!
      lg = 0
      do i = 1, izt
         do j = 1, ith
            lg = lg + 1
            f(lg) = f3a(i, j)
         end do
      end do
      sin_type = 0; cos_type = 1
      call toFourier
      f3a_nm(:) = fnm(:)
!
      lg = 0
      do i = 1, izt
         do j = 1, ith
            lg = lg + 1
            f(lg) = f3b(i, j)
         end do
      end do
      sin_type = 0; cos_type = 1
      call toFourier
      f3b_nm(:) = fnm(:)
!
      lg = 0
      do i = 1, izt
         do j = 1, ith
            lg = lg + 1
            f(lg) = f3c(i, j)
         end do
      end do
      sin_type = 0; cos_type = 1
      call toFourier
      f3c_nm(:) = fnm(:)
!
!
!     Write out coefficient spectra if half way out in flux
!
      if (ir .eq. ir_fine_scl/2) then
         do mn = 1, mnmx
            write (8, '(f6.1,2x,f6.1,4(2x,e15.7))') rm(mn), rn(mn),&
            &f1_nm(mn), f3a_nm(mn),&
            &f3b_nm(mn), f3c_nm(mn)
         end do
      end if
!
!     Build A and B matrices
!
      do i = 1, mn_col
         do j = 1, mn_col
            bmat(i, j) = 0
            amat(i, j) = 0
            if (i .lt. j .and. ipos_def_sym) cycle     !i.e. symmetric storage mode: only need bottom half
            ni = in_col(i)
            nj = in_col(j)
            mi = im_col(i)
            mj = im_col(j)
            do ieq = 1, mnmx
               meq = rm(ieq)
               neq = rn(ieq)
               call ccc_convolve(ccci, mi, ni, mj, nj, meq, neq)
               call scs_convolve(scsi, mi, ni, mj, nj, meq, neq)
               bmat(i, j) = bmat(i, j) - ccci*f1_nm(ieq)*mu0_rho_ion(ir)
               amat(i, j) = amat(i, j)&
               &- (scsi*(f3a_nm(ieq)*rm_col(i)*rm_col(j)&
               &- f3b_nm(ieq)*rm_col(j)*rn_col(i)&
               &- f3b_nm(ieq)*rn_col(j)*rm_col(i)&
               &+ f3c_nm(ieq)*rn_col(j)*rn_col(i)))
            end do              !ieq = 1,mnmx
         end do                 !j=1,mn_col
      end do                  !i=1,mn_col
!
!     Call matrix eigenvalue solver
!

      egl = 1.d-2; egu = 0.6d0; abstol = 1.d-8
      il = 0; iu = 0

      if (.NOT. ipos_def_sym) then
         call dggev('N', 'V', mn_col, amat, mn_col, bmat, mn_col, alfr, alfi,&
         &beta, vl, ldvl, vr, ldvr, work, lwork, info)
      else if (ipos_def_sym .and. .not. subset_eq) then
         call dsygv(iopt, jobz, 'L', mn_col, amat, mn_col, bmat, mn_col,&
         &omega, work, lwork, info)
         if (info .ne. 0) write (*, '("info = ",i8)') info

      else if (ipos_def_sym .and. subset_eq) then
         call dsygvx(iopt, 'V', 'V', 'L', mn_col, amat, mn_col, bmat, mn_col,&
         &egl, egu, il, iu, abstol, m_red, omega, vr, mn_col,&
         &work, lwork, iwork, ifail, info)
         if (info .ne. 0) write (*, '("info = ",i8)') info

      end if
!
      do i = 1, mn_col
         if (iopt .eq. 1) then
            do j = 1, mn_col
               if (ipos_def_sym .and. jobz .eq. 'V' .and. .not. subset_eq)&
               &eig_vect(j) = abs(amat(j, i))
               if (.NOT. ipos_def_sym) eig_vect(j) = abs(vr(j, i))
!        if(ir .eq. 5 .and. i .eq. mn_col/2)
!     >   write(*,'(e15.7,2x,i4,2x,i4)') eig_vect(j),
!     >   im_col(j), in_col(j)
            end do

            eig_max = -1.d+30
            do j = 1, mn_col
               if (eig_vect(j) .gt. eig_max) then
                  eig_max = eig_vect(j)
                  j_max_index = j
               end if
            end do
            m_emax = im_col(j_max_index)
            n_emax = in_col(j_max_index)

!       if(ir .eq. 5 .and. i .eq. mn_col/2)
!     >   write(*,'(e15.7,2x,i4,2(2x,i4))') eig_max,m_emax,n_emax,
!     >    j_max_index

            if (ipos_def_sym) then
               write (21, '(2(e15.7,2x),i4,2x,i4)') r_pt, sqrt(abs(omega(i))),&
               &m_emax, n_emax
            else if (.NOT. ipos_def_sym) then
               write (21, '(4(e15.7,2x),i4,2x,i4)') r_pt, alfr(i),&
               &alfi(i), beta(i), m_emax, n_emax
            end if
         else if (iopt .eq. 0) then
            if (ipos_def_sym) then
               write (21, '(e15.7,2x,e15.7)') r_pt, sqrt(abs(omega(i)))
            else if (.NOT. ipos_def_sym) then
               write (21, '(e15.7,3(2x,e15.7))') r_pt, alfr(i),&
               &alfi(i), beta(i)
            end if
         end if
      end do          !do i=1,mn_col

   end do       !ir=1,ir_fine_scl
!
!
   close (unit=20)
   close (unit=21)
   close (unit=8)
   if (mype .eq. 0) write (*,&
   &'("modes = ",i5,2x,"no. of radial points = ",i5)')&
   &mn_col, ir_fine_scl
   if (mype .eq. 0) write (7, '(i5,3(2x,i5))') iopt, mn_col, ir_fine_scl,&
   &isym_opt
   if (mype .eq. 0) close (unit=7)
   call trg_deallocate
   deallocate (alpha, beta, aux)
   if (mype .eq. 0) then
      call date_and_time(date, time, zone, values)
      write (*, *) time
   end if
#if defined (PARALLEL)
   call MPI_FINALIZE(ier)
#endif
   call post_process
   stop
end program tae_continua
