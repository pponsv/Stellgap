!                              D   I   S   C   L   A   I   M   E   R
!
!       You are using a BETA version of the program metric_element_create.f, which is currently
!       under development by D. A. Spong of the Fusion Energy Division,
!       Oak Ridge National Laboratory.  Please report any problems or comments
!       to him.  As a BETA version, this program is subject to periodic change
!       and improvement without notice.
!
!   7/29/2010: A test was included to see if all the surfaces to the edge were included in
!              the xbooz_xform calculation. i.e., if(rmncbh(1,nsd) .eq. 0.)... This can also
!              be seen from running ncdump on the boozmn file and comparing ns_b to comput_surfs.
!              If the last surface is not included, then there can be problems with the outer
!              surface in boozmn since averages between surfaces are made. This particularly
!              shows up when xmetric is used for the outer surface metric calculation. In this
!              case, the outer surface area can be too low by a factor of ~4 due to the outer R
!              and z being down by 1/2 (from averaging over a 0 plus the correct value).
!
use read_boozer_mod, iota_bw => iota_b, bmnc_bw => bmnc_b, &
&rmnc_bw => rmnc_b, zmns_bw => zmns_b, gmn_bw => gmnc_b, &
&pmns_bw => pmns_b, pres_bw => pres_b, phip_bw => phip_b, &
&bvco_bw => bvco_b, buco_bw => buco_b, ixm_bw => ixm_b, &
&ixn_bw => ixn_b
use constants
!
!
implicit none

!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
real(r8), parameter :: p5 = 0.5_dp, one = 1.0_dp
!REAL(r8), PARAMETER :: twopi = 6.28318530717958623_dp
! REAL(r8), PARAMETER :: mu_0 = 2.0e-7_dp*twopi
real(r8), parameter :: zero = 0.0_dp, two = 2.0_dp
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
REAL(r8), DIMENSION(:), ALLOCATABLE :: hiota, hpres, &
&hphip, hjtor, hjpol
REAL(r8), DIMENSION(:), ALLOCATABLE :: iotaf, jpolf, &
&jtorf, phipf, presf
REAL(r8), DIMENSION(:), ALLOCATABLE :: iotapf, jpolpf, &
&jtorpf, phippf, prespf
REAL(r8), DIMENSION(:), ALLOCATABLE :: xmb, xnb, xm, xn
REAL(r8), DIMENSION(:), ALLOCATABLE :: radii, radii_flux
REAL(r8), DIMENSION(:, :), ALLOCATABLE :: rmncbf, zmnsbf, &
&pmnsbf, bmncbf
REAL(r8), DIMENSION(:, :), ALLOCATABLE :: rmncpbf, zmnspbf, &
&pmnspbf, bmncpbf
REAL(r8), DIMENSION(:, :), ALLOCATABLE :: rmncbh, zmnsbh, &
&pmnsbh, bmncbh
REAL(r8), DIMENSION(:, :), ALLOCATABLE :: bfield, &
&bfields, bfieldze, bfieldth, rjacob, rjacobze, rjacobth, &
&rjacobs, gsssup, gttsup, gzzsup, gstsup, &
&gszsup, gtzsup, brho
REAL(r8), DIMENSION(:), ALLOCATABLE :: jprl_coef0, &
&jprl_coef1, jprl_coef2
REAL(r8), DIMENSION(3, 3) :: mat_lowr
REAL(r8), DIMENSION(3, 3) :: mat_upr
REAL(r8), DIMENSION(3, 3) :: mat_diag
integer :: nsurf
real(r8) :: phipc, iotac, jtorc, jpolc, &
&presc, prespc, jtorpc, jpolpc, phippc, iotapc
real(r8) :: phis, phize, phith, rboo, rs, &
&rze, rth, zs, zze, zth, arg, ccosi, ssine, num1, &
&num2, ibf2, ibf3, rboo2, rjac2, rjac2i, gtssub, gstsub, &
&gzssub, gszsub, gtzsub, gztsub, gttsub, gsssub, &
&gzzsub, gaasup, gsasup, gtssup, gzssup, gztsup, &
&den1, cka, beta, t1, t2, t3, t4, t5, cks, &
&thetang, zetang, det, error, rjac2_avg, det_avg, &
&zboo, phiboo, xr, yr, zr
INTEGER :: nfp, nsd, mnboz, i, j, k, mn, kz, kt, ks
REAL(r8) :: r0, b0, amin, beta0, r0max, r0min, &
&aspect, betaxis, ohs2, mat_test_diag, mat_test_offdiag
character arg1*40, warg1*45, bozout*45
character*1 tb
integer nargs, numargs, itheta, izeta, nznt
integer iargc, unit_no, istat, ierr, ig, lf, is
real viz_flux   !for plotting interior flux surfaces in AVS
real surf_area_element, surf_area_total
logical lasym, viz, test_jacob, test_upr_lowr, &
&make_stellgap_data, make_full_torus, surf_compute
!-----------------------------------------------
tb = char(9); itheta = 80; izeta = 80; viz = .true.
viz_flux = 0.5   !selects surface for AVS data - sync with metric_element_create.f
test_jacob = .false.; test_upr_lowr = .false.
make_stellgap_data = .true.
make_full_torus = .true.
surf_compute = .true.
numargs = iargc()
call getarg(1, arg1)
if (numargs .ne. 1) then
   print *, ' MUST ENTER FILE SUFFIX ON COMMAND LINE'
   stop
end if
if (surf_compute) then
   open (unit = 47, file = "surf_area_elements", status = "unknown")
end if
if (make_stellgap_data) then
   open (unit = 20, file = "tae_data_boozer", status = "unknown")
end if
if (viz) then
   open (unit = 12, file = "cart_coords", status = "unknown")
   open (unit = 14, file = "metrics", status = "unknown")
end if
!       open(unit=15,file="ae_metric.dat",status="unknown")
open (unit = 15, file = "ae_metric.dat", &
&status = "unknown")
!
!
warg1 = arg1
!      write(*,*) warg1
!      call read_wout_file(warg1,ierr)
call read_boozer_file(warg1, ierr)
if (istat .ne. 0) stop 22

nfp = nfp_b
nsd = ns_b
aspect = aspect_b
r0max = rmax_b
r0min = rmin_b
betaxis = betaxis_b
mnboz = mnboz_b
!
if (surf_compute) then
   write (47, *) nfp, izeta, itheta
end if
!
!...   IOTA, PRES, PHIP (= -PHIP_VMEC), JTOR (=I = -BUCO_VMEC) and
!         JPOL (=J = BVCO_VMEC) are ALL on HALF-MESH!!!

allocate (hiota(nsd), hpres(nsd), hjpol(nsd), hjtor(nsd), &
&hphip(nsd), jprl_coef0(nsd), jprl_coef1(nsd), &
&jprl_coef2(nsd), stat = istat)
if (istat .ne. 0) stop 23
do k = 1, nsd
   hiota(k) = iota_bw(k)
   hpres(k) = mu_0 * pres_bw(k)  ! in VMEC versions > 6.00 pressure is given in pascals
   hphip(k) = -phip_bw(k)        ! toroidal fluxes have REVERSED sign respect to VMEC!!
   hjpol(k) = bvco_bw(k)
   hjtor(k) = -buco_bw(k)        ! toroidal fluxes have REVERSED sign respect to VMEC!!
   !       write(*,'(i4,5(2x,e12.5))') k,hiota(k),hpres(k),hphip(k),
   !     1     hjpol(k),hjtor(k)
end do
r0 = (r0max + r0min) / two
amin = r0 / aspect
beta0 = betaxis
!       if (beta0 .le. zero) stop 'Beta0 <= 0'
!       b0 = sqrt((two/beta0)*(1.5_dp*hpres(2)-.5_dp*hpres(3)))
write (*, '(///)')
write (*, *) nfp, nsd, aspect
write (*, '(/)')
write (*, *) r0max, r0min, betaxis, mnboz
write (*, '(/)')
write (*, *) r0, amin, beta0

allocate (xm(mnboz), xn(mnboz), stat = istat)
if (istat .ne. 0) stop 24
do mn = 1, mnboz
   xm(mn) = ixm_bw(mn)
   xn(mn) = ixn_bw(mn)
   !       write(*,*) mn,xm(mn),xn(mn)
end do
!...   RMN, ZMN, PMN and BMN are ALL on HALF-MESH

allocate (rmncbh(mnboz, nsd), zmnsbh(mnboz, nsd), &
&pmnsbh(mnboz, nsd), bmncbh(mnboz, nsd), stat = istat)
if (istat .ne. 0) stop 25

zmnsbh = zero; rmncbh = zero; pmnsbh = zero; bmncbh = zero

do k = 1, nsd
   do mn = 1, mnboz
      bmncbh(mn, k) = bmnc_bw(mn, k)
      rmncbh(mn, k) = rmnc_bw(mn, k)
      zmnsbh(mn, k) = zmns_bw(mn, k)
      pmnsbh(mn, k) = pmns_bw(mn, k)
   end do
end do
call read_boozer_deallocate
allocate (radii(nsd), radii_flux(nsd), stat = istat)
do i = 1, nsd
   radii(i) = sqrt(real(i - 1, r8) / (nsd - 1))     !  r/a= (radii_flux)**1/2
   radii_flux(i) = real(i - 1, r8) / (nsd - 1)    !  s = radii_flux
end do

if (rmncbh(1, nsd) .eq. 0.) then
   write (*, '("May need to add last surface to xbooz calculation")')
   !         stop 55
end if
!=================
!    MIGRATE TO FULL MESH
!=================

!...  Store surface quantites on RADIAL full mesh

allocate (iotaf(nsd), presf(nsd), jtorf(nsd), jpolf(nsd), &
&phipf(nsd), stat = k)
if (k .ne. 0) stop 'Allocation error'

iotaf = zero; presf = zero; jtorf = zero; jpolf = zero; phipf = zero

iotaf(2:nsd - 1) = p5 * (hiota(2:nsd - 1) + hiota(3:nsd))
presf(2:nsd - 1) = p5 * (hpres(2:nsd - 1) + hpres(3:nsd))
jtorf(2:nsd - 1) = p5 * (hjtor(2:nsd - 1) + hjtor(3:nsd))
jpolf(2:nsd - 1) = p5 * (hjpol(2:nsd - 1) + hjpol(3:nsd))
phipf(2:nsd - 1) = p5 * (hphip(2:nsd - 1) + hphip(3:nsd))

!...  Evaluate and store surface quantities derivatives on RADIAL full mesh

allocate (iotapf(nsd), jpolpf(nsd), jtorpf(nsd), &
&phippf(nsd), prespf(nsd), stat = k)
if (k .ne. 0) stop 'Allocation error in get_ballooning_grate'

iotapf = zero; prespf = zero; jtorpf = zero; jpolpf = zero; phippf = zero

!      ohs2 = 2.0_dp*dble(nsd-1)                       ! ds to differentiate on RADIAL half mesh
ohs2 = dble(nsd - 1)                                            ! 2 comes from interpolating to radial FULL mesh
iotapf(2:nsd - 1) = ohs2 * (hiota(3:nsd) - hiota(2:nsd - 1))
prespf(2:nsd - 1) = ohs2 * (hpres(3:nsd) - hpres(2:nsd - 1))
jtorpf(2:nsd - 1) = ohs2 * (hjtor(3:nsd) - hjtor(2:nsd - 1))
jpolpf(2:nsd - 1) = ohs2 * (hjpol(3:nsd) - hjpol(2:nsd - 1))
phippf(2:nsd - 1) = ohs2 * (hphip(3:nsd) - hphip(2:nsd - 1))

deallocate (hiota, hpres, hjtor, hjpol, hphip)

!...  store (m,n)-descriptors

allocate (xnb(mnboz), xmb(mnboz), stat = k)
if (k .ne. 0) stop 'Allocation error in get_ballooning_grate'
xnb = xn
xmb = xm

!...  EVALUATE AND STORE Boozer Fourier coefficients and
!           their derivatives on RADIAL full mesh

allocate (rmncbf(mnboz, nsd), zmnsbf(mnboz, nsd), pmnsbf(mnboz, nsd), &
&bmncbf(mnboz, nsd), rmncpbf(mnboz, nsd), &
&zmnspbf(mnboz, nsd), pmnspbf(mnboz, nsd), &
&bmncpbf(mnboz, nsd), stat = k)
if (k .ne. 0) stop 'Allocation error'

rmncbf = zero; zmnsbf = zero; pmnsbf = zero; bmncbf = zero
rmncpbf = zero; zmnspbf = zero; pmnspbf = zero; bmncpbf = zero
do mn = 1, mnboz
   do k = 1, nsd - 1

      !...   Boozer Fourier coefficients on RADIAL full mesh

      rmncbf(mn, k) = p5 * (rmncbh(mn, k + 1)&
      & + rmncbh(mn, k))
      zmnsbf(mn, k) = p5 * (zmnsbh(mn, k + 1)&
      & + zmnsbh(mn, k))
      pmnsbf(mn, k) = p5 * (pmnsbh(mn, k + 1)&
      & + pmnsbh(mn, k))
      bmncbf(mn, k) = p5 * (bmncbh(mn, k + 1)&
      & + bmncbh(mn, k))

      !...   Boozer Fourier coefficients radial derivatives on RADIAL full mesh

      rmncpbf(mn, k) = ohs2 * (rmncbh(mn, k + 1)&
      & - rmncbh(mn, k))
      zmnspbf(mn, k) = ohs2 * (zmnsbh(mn, k + 1)&
      & - zmnsbh(mn, k))
      pmnspbf(mn, k) = ohs2 * (pmnsbh(mn, k + 1)&
      & - pmnsbh(mn, k))
      bmncpbf(mn, k) = ohs2 * (bmncbh(mn, k + 1)&
      & - bmncbh(mn, k))

   end do
end do

deallocate (xn, xm, rmncbh, zmnsbh, pmnsbh, bmncbh, stat = k)
ig = itheta * izeta
allocate (bfield(ig, nsd), bfields(ig, nsd), bfieldze(ig, nsd), &
&bfieldth(ig, nsd), rjacob(ig, nsd), &
&rjacobze(ig, nsd), rjacobth(ig, nsd), &
&rjacobs(ig, nsd), gsssup(ig, nsd), gttsup(ig, nsd), gzzsup(ig, nsd), &
&gstsup(ig, nsd), gszsup(ig, nsd), &
&gtzsup(ig, nsd), brho(ig, nsd), stat = istat)

flux_surface : do ks = 2, nsd - 1

   !...   values of surface quantities at current surface

   phipc = phipf(ks)                            ! toroidal magnetic flux
   iotac = iotaf(ks)                            ! iota
   iotapc = iotapf(ks)                          ! radial iota derivative
   jtorc = jtorf(ks)                            ! torodial current density
   jpolc = jpolf(ks)                            ! poloidal current density
   presc = presf(ks)                            ! pressure
   prespc = prespf(ks)                          ! radial pressure gradient
   jtorpc = jtorpf(ks)                          ! radial derivative toroidal current density
   jpolpc = jpolpf(ks)                          ! radial derivative poloidal current density
   phippc = phippf(ks)                          ! radial derivative toroidal flux
   jprl_coef0(ks) = (jpolc * jtorpc - jpolpc * jtorc)&
   & / ((jpolc + jtorc * iotac) * phipc)
   jprl_coef1(ks) = jtorc / ((jpolc + jtorc * iotac) * phipc)
   jprl_coef2(ks) = -jpolc / ((jpolc + jtorc * iotac) * phipc)
   if (make_stellgap_data) then
      write (20, 21) ks, iotac, phipc, jtorc, jpolc
21    format(1x, i3, 4(2x, e15.7))
   end if
   !=============
   !   BEGIN FOURIER INVERSION
   !=============
   lf = 0
   error = 0.; rjac2_avg = 0.; det_avg = 0.
   surf_area_total = 0.
   toroidal : do kz = 1, izeta
      poloidal : do kt = 1, itheta
         lf = lf + 1
         thetang = twopi * dble(kt - 1) / dble(itheta - 1)
         zetang = twopi * dble(kz - 1) / (dble(nfp) * dble(izeta - 1))
         !...   initialize for Fourier inversion
         bfield(lf, ks) = zero
         bfieldze(lf, ks) = zero
         bfields(lf, ks) = zero
         bfieldth(lf, ks) = zero
         phis = zero
         phize = one              ! since Phi_cyl  =  Phi_VMEC  = Phi_BOOZ + Pmn
         phith = zero
         rboo = zero; zboo = zero; phiboo = zero
         rs = zero
         rze = zero
         rth = zero
         zs = zero
         zze = zero
         zth = zero

         fourier : do j = 1, mnboz             ! Fourier invert B, R, Z, Phi and derivatives

            arg = xmb(j) * thetang - zetang * xnb(j)
            ccosi = cos(arg)
            ssine = sin(arg)
            bfield(lf, ks) = bfield(lf, ks) + bmncbf(j, ks) * ccosi                !magnetic field magnitude
            bfields(lf, ks) = bfields(lf, ks) + bmncpbf(j, ks) * ccosi             ! ..... radial derivative
            bfieldze(lf, ks) = bfieldze(lf, ks) + xnb(j) * bmncbf(j, ks) * ssine       ! ..... zeta derivative
            bfieldth(lf, ks) = bfieldth(lf, ks) - xmb(j) * bmncbf(j, ks) * ssine       ! ..... theta derivative
            rboo = rboo + rmncbf(j, ks) * ccosi                    ! cylindrical R
            zboo = zboo + zmnsbf(j, ks) * ssine                    ! cylindrical Z
            phiboo = phiboo + pmnsbf(j, ks) * ssine                ! cylindrical Phi
            rth = rth - rmncbf(j, ks) * xmb(j) * ssine               ! ..... theta derivative
            rze = rze + rmncbf(j, ks) * xnb(j) * ssine               ! ..... zeta derivative
            rs = rs + rmncpbf(j, ks) * ccosi                       ! ..... radial derivative
            zth = zth + zmnsbf(j, ks) * xmb(j) * ccosi               ! cylindrical Z: theta derivative
            zze = zze - zmnsbf(j, ks) * xnb(j) * ccosi               ! ..... zeta derivative
            zs = zs + zmnspbf(j, ks) * ssine                       ! ..... radial derivative
            phith = phith + pmnsbf(j, ks) * xmb(j) * ccosi           ! cylindrical Phi: theta derivative
            phize = phize - pmnsbf(j, ks) * xnb(j) * ccosi           ! ..... zeta derivative
            phis = phis + pmnspbf(j, ks) * ssine                   ! ..... radial derivative

         end do fourier

         phiboo = phiboo + zetang
         !============
         !   END FOURIER INVERSION
         !============

         num1 = (iotac * jtorc - jpolc) * phipc
         num2 = (iotac * jtorc - jpolc) * phippc + (iotapc * jtorc + &
         &iotac * jtorpc - jpolpc) * phipc

         ibf2 = one / bfield(lf, ks)**2
         ibf3 = -2.0_dp / bfield(lf, ks)**3

         !...   Jacobian from cylindrical to Boozer  !note: this is sqrt(g), not 1/sqrt(g)

         rjacob(lf, ks) = num1 * ibf2                                 ! jacobian from cyl. to Boozer
         rjacobze(lf, ks) = num1 * bfieldze(lf, ks) * ibf3              ! ..... zeta derivative
         rjacobth(lf, ks) = num1 * bfieldth(lf, ks) * ibf3              ! ..... theta derivative
         rjacobs(lf, ks) = num1 * bfields(lf, ks) * ibf3 + num2 * ibf2      ! ..... radial derivative

         rboo2 = rboo**2
         rjac2 = rjacob(lf, ks)**2
         rjac2i = one / rjac2

         !...   Boozer lower metric elements

         gtssub = rth * rs + zth * zs + rboo2 * phith * phis
         gstsub = gtssub
         gzssub = rze * rs + zze * zs + rboo2 * phize * phis
         gszsub = gzssub
         gtzsub = rth * rze + zth * zze + rboo2 * phith * phize
         gztsub = gtzsub
         gttsub = (rth**2) + (zth**2) + rboo2 * (phith**2)
         gsssub = (rs**2) + (zs**2) + rboo2 * (phis**2)
         gzzsub = (rze**2) + (zze**2) + rboo2 * (phize**2)
         if (test_upr_lowr) then
            mat_lowr(1, 1) = gsssub; mat_lowr(1, 2) = gstsub; mat_lowr(1, 3) = gszsub
            mat_lowr(2, 1) = gtssub; mat_lowr(2, 2) = gttsub; mat_lowr(2, 3) = gtzsub
            mat_lowr(3, 1) = gzssub; mat_lowr(3, 2) = gztsub; mat_lowr(3, 3) = gzzsub
         end if
         if (ks .eq. nsd - 1 .and. surf_compute) then
            surf_area_element = sqrt(gttsub * gzzsub - gtzsub * gtzsub)
            surf_area_total = surf_area_total + surf_area_element&
            & * (twopi**2) / (dble(itheta - 1) * dble(izeta - 1))
            write (47, '(3e16.8)') zetang, thetang, surf_area_element
         end if

         !...   Consistency test:
         if (test_jacob) then
            det = gsssub * (gttsub * gzzsub - gztsub * gtzsub)&
            & + gstsub * (gtzsub * gzssub - gtssub * gzzsub)&
            & + gszsub * (gtssub * gztsub - gttsub * gzssub)
            det_avg = det_avg + det / dble(ig)
            rjac2_avg = rjac2_avg + rjac2 / dble(ig)
            error = error + (det - rjac2) / (rjac2 * dble(ig))
         end if

         !...   Boozer upper metric elements

         gsssup(lf, ks) = (gttsub * gzzsub - gtzsub * gztsub) * rjac2i
         gttsup(lf, ks) = (gsssub * gzzsub - gszsub * gzssub) * rjac2i
         gzzsup(lf, ks) = (gsssub * gttsub - gstsub * gtssub) * rjac2i
         gstsup(lf, ks) = (gztsub * gzssub - gstsub * gzzsub) * rjac2i
         gtssup = gstsup(lf, ks)
         gszsup(lf, ks) = (gtssub * gtzsub - gttsub * gzssub) * rjac2i
         gzssup = gszsup(lf, ks)
         gtzsup(lf, ks) = (gstsub * gzssub - gsssub * gztsub) * rjac2i
         gztsup = gtzsup(lf, ks)
         brho(lf, ks) = -(gstsup(lf, ks) * jtorc + gszsup(lf, ks) * jpolc)&
         & / gsssup(lf, ks)
         if (test_upr_lowr) then
            mat_upr(1, 1) = gsssup(lf, ks); mat_upr(1, 2) = gstsup(lf, ks)
            mat_upr(1, 3) = gszsup(lf, ks)
            mat_upr(2, 1) = gtssup; mat_upr(2, 2) = gttsup(lf, ks)
            mat_upr(2, 3) = gtzsup(lf, ks)
            mat_upr(3, 1) = gzssup; mat_upr(3, 2) = gztsup
            mat_upr(3, 3) = gzzsup(lf, ks)
            mat_diag = matmul(mat_lowr, mat_upr)
            mat_test_diag = (mat_diag(1, 1) + mat_diag(2, 2) + mat_diag(3, 3)) / 3.
            mat_test_offdiag = (mat_diag(1, 2) + mat_diag(1, 3) + mat_diag(2, 1)&
            & + mat_diag(2, 3) + mat_diag(3, 1) + mat_diag(3, 2)) / 6.
            if (abs(mat_test_diag - 1.) .gt. .04) then
               write (*, '(e15.7,2x,e15.7)') mat_test_diag, mat_test_offdiag
            end if
         end if

         !...   Upper metric elements in (s,alpha,phi)-coord.

         gaasup = gttsup(lf, ks) + (iotac**2) * gzzsup(lf, ks)&
         & - 2. * iotac * gtzsup(lf, ks)
         gsasup = gstsup(lf, ks) - iotac * gszsup(lf, ks)

         !...   Boozer geodesic curvature

         den1 = 2.0_dp * (iotac * jtorc - jpolc) * rjacob(lf, ks)
         cka = (jpolc * rjacobth(lf, ks) + jtorc * rjacobze(lf, ks)) / den1

         !...   Boozer normal curvature

         beta = (gstsup(lf, ks) * jtorc - gszsup(lf, ks) * jpolc) / gsssup(lf, ks)
         t1 = iotapc * jtorc + iotac * jtorpc - jpolpc + (phippc / phipc)&
         & * (iotac * jtorc - jpolc)
         t2 = 2.0_dp * prespc / phipc
         t3 = -(iotac * jtorc - jpolc) * rjacobs(lf, ks)
         t4 = -beta * iotac * rjacobth(lf, ks)
         t5 = -beta * rjacobze(lf, ks)
         cks = (rjacob(lf, ks) * t1 + rjac2 * t2 + t3 + t4 + t5) / den1
         !
         !      Section to write out 3D data for AVS:
         !
         if (ks .eq. (nsd - 1) .and. viz) then
            xr = rboo * cos(phiboo)
            yr = rboo * sin(phiboo)
            zr = zboo
            write (12, '(e15.7,2(2x,e15.7))') xr, yr, zr
            write (14, '(e15.7,7(2x,e15.7))') gsssup(lf, ks), gttsup(lf, ks), &
            &gzzsup(lf, ks), gstsup(lf, ks), gszsup(lf, ks), gtzsup(lf, ks), &
            &bfield(lf, ks)
         end if    !ks .eq. nsd-1
         !
         if (make_stellgap_data) then
            write (20, '(1x,4(e24.12,2x),e24.12)') thetang, &
            &zetang, bfield(lf, ks), gsssup(lf, ks), &
            &gtzsup(lf, ks)

            write (20, '(1x,4(e24.12,2x),e24.12)') gttsup(lf, ks), &
            &gzzsup(lf, ks), gstsup(lf, ks), gszsup(lf, ks), &
            &rjacob(lf, ks)
         end if
      end do poloidal
   end do toroidal
   if (ks .eq. 3) nznt = lf
   if (test_jacob) then
      write (*, '(i3,3(2x,e15.7))') ks, error, det_avg, rjac2_avg
   end if
end do flux_surface
!
!
!
if (surf_compute) then
   write (*, '("Total surface area = ",e15.7,"(m**2)")')&
   &surf_area_total
end if
write (15, *) (nsd - 2), izeta, itheta, nznt
do ks = 2, nsd - 1
   write (15, 19) iotaf(ks), iotapf(ks), jpolf(ks), jpolpf(ks), &
   &jtorf(ks), jtorpf(ks), phipf(ks), phippf(ks)
   !        write(*,*) presf(ks)
end do
do ks = 2, nsd - 1
   do lf = 1, nznt
      !
      write (15, 18) rjacob(lf, ks), bfield(lf, ks), &
      &gsssup(lf, ks), &
      &gttsup(lf, ks), gzzsup(lf, ks), &
      &gstsup(lf, ks), &
      &gszsup(lf, ks), gtzsup(lf, ks), &
      &bfields(lf, ks), &
      &bfieldth(lf, ks), bfieldze(lf, ks)
   end do
end do
!
!  Write out coefficients and B_rho that will be used subsequently
!  in AE3D to form J_prl/B:
!
do ks = 2, nsd - 1
   write (15, 48) jprl_coef0(ks), jprl_coef1(ks), &
   &jprl_coef2(ks), prespf(ks)
end do
!
!
do ks = 2, nsd - 1
   do lf = 1, nznt
      write (15, 49) brho(lf, ks)
   end do
end do
!
!
18 format(e15.7, 10(2x, e15.7))
19 format(e15.7, 7(2x, e15.7))
48 format(e15.7, 3(2x, e15.7))
49 format(e15.7)

if (make_full_torus) then
   itheta = 300; izeta = 300
   !       ks = nsd - 2
   ks = (nsd - 1) * viz_flux
   !       ks = 0.2*nsd
   open (unit = 17, file = "full_torus_coords", status = "unknown")
   do kz = 1, izeta
      do kt = 1, itheta
         thetang = twopi * dble(kt - 1) / dble(itheta - 1)
         !       zetang = twopi*dble(kz-1)/dble(izeta-1)
         zetang = 0.9 * twopi * dble(kz - 1) / dble(izeta - 1)
         rboo = 0.d0; zboo = 0.d0; phiboo = 0.d0
         do j = 1, mnboz             ! Fourier invert B, R, Z, Phi and derivatives
            arg = xmb(j) * thetang - zetang * xnb(j)
            ccosi = cos(arg)
            ssine = sin(arg)
            rboo = rboo + rmncbf(j, ks) * ccosi                    ! cylindrical R
            zboo = zboo + zmnsbf(j, ks) * ssine                    ! cylindrical Z
            phiboo = phiboo + pmnsbf(j, ks) * ssine                ! cylindrical Phi
         end do
         phiboo = phiboo + zetang
         xr = rboo * cos(phiboo)
         yr = rboo * sin(phiboo)
         zr = zboo
         write (17, '(e15.7,2(2x,e15.7))') xr, yr, zr
      end do
   end do

   open (unit = 18, file = "torus_slice_coords", status = "unknown")
   write (*, *) nsd - 2, itheta
   !       do is = 2,6
   !        ks = is*is
   !        if(is .eq. 6) ks = nsd-2
   !       do ks = 1,nsd-2
   !
   !       zetang = 0.d0
   !       zetang = twopi/50.
   zetang = 0.9 * twopi
   do ks = 2, int((nsd - 1) * viz_flux)
      !       do ks = 1,nsd-1
      do kt = 1, itheta
         thetang = twopi * dble(kt - 1) / dble(itheta - 1)
         rboo = 0.d0; zboo = 0.d0; phiboo = 0.d0
         do j = 1, mnboz             ! Fourier invert B, R, Z, Phi and derivatives
            arg = xmb(j) * thetang - zetang * xnb(j)
            ccosi = cos(arg)
            ssine = sin(arg)
            !         if(ks .ne. 1) then
            rboo = rboo + rmncbf(j, ks) * ccosi                    ! cylindrical R
            zboo = zboo + zmnsbf(j, ks) * ssine                    ! cylindrical Z
            phiboo = phiboo + pmnsbf(j, ks) * ssine                ! cylindrical Phi
            !         else if(ks .eq. 1) then
            !          rboo = rboo+rmncbf(j,ks+1)*ccosi                    ! cylindrical R
            !          zboo = 0.
            !          phiboo = 0.
            !         end if
         end do
         phiboo = phiboo + zetang
         xr = rboo * cos(phiboo)
         yr = rboo * sin(phiboo)
         zr = zboo
         write (18, '(e16.7,2(2x,e15.7))') xr, yr, zr
      end do
   end do
end if
stop
end
