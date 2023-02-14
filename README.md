# Stellgap

Stellgap calculates the shear Alfvén gap structure for 3D configurations (stellarators, RFPs, 3D tokamaks)

These codes are used to calculate shear Alfven continua for 3D configurations, both with and without sound wave coupling effects. The associated paper is D. A. Spong, R. Sanchez, A. Weller, "Shear Alfvén continua in stellarators," Phys. Plasmas 10 (2003) 3217–3224.
The workflow is as follows:

(a) prepare a VMEC equilibrium for the case of interest

(b) Run xbooz_xform (from Stellopt code suite) to convert from VMEC coordinates to Boozer coordinates. This run needs the in_booz.\* file, which tells it what range of m/n modes to use and the selected surfaces to use from the VMEC run. Typically, the first and last of the VMEC surfaces are removed, because they can sometimes have noisy data.

(c) Using the boozmn._ file produced in (b), run xmetric_ver_ to extract needed data and calculate the metric elements used in the continuum calculation. This produces a file called tae_data_boozer which is used as input for xstgap.

(d) The continuum calculation is done by xstgap*. There are two versions: one with sound wave couplings (xstgap_snd_ver*) and one without (xstgap). In addition to ae_data_boozer, these use the input files fourier.dat and plasma.dat. Fourier.dat specifies the set of fourier modes used to represent the continuum eigenmodes and plasma.dat contains information/profiles for the the plasma. The first line of fourier.dat gives the field periods and the surface grid parameters ith and izt used for calculating theta/zeta dependent coefficients of the continuum equation. It is important that these values for ith and izt are exactly the same as were used in the metric_element_create.f code.

(e) After running xstgap, the post-processing code (either post_process.f or post_process_snd.f) should be run. This produces a text file called alfven_post which contains columns of radial coordinate, frequency, dominant m and dominant n. Typically, the frequency vs. radius is plotted as a scatter plot. The code stelgp_to_silo.f is provided to convert data from alfven_post to a Visit silo file. This allows plotting the continua using the Visit software with the option to color code the points with either the dominant m or n.

Compilation scripts are provided for the different codes in the files whose name begins with “bld”. These will need to be edited, depending on what compiler is being used and where the needed libraries are located. Stellgap is constructed so that it can be compiled either as a serial version or a parallel version using precompilation flags. When running in parallel the number of surfaces requested will be divided by the number of processors and each group of surfaces allocated to a different processor. At the end of the run, all processors will write their results out to separate files, with names containing the processor number. These must be concatenated together in the post-processing step (e.g., as done in the post_process_snd.f code).

# Stellgap_ver5:

                            D   I   S   C   L   A   I   M   E   R

    You are using a BETA version of the program stellgap.f, which is currently
    under development by D. A. Spong of the Fusion Energy Division,
    Oak Ridge National Laboratory.  Please report any problems or comments
    to him.  As a BETA version, this program is subject to periodic change
    and improvement without notice.


    Note: this code uses normalized toroidal flux (variables: rho, r_pt) as
    an independent variable. If a radius-like variable is desired
    (i.e., sqrt(tor_flux) or radius based on sqrt(area/pi)), this must be
    constructed in subsequent codes that use the AE3D eigenfunction data

    Before running, check that the irads variable is set consistently with the
    number of surfaces in the tae_boozer_data file. irads should be the number
    of surfaces in the original VMEC run minus 1 or 2 to avoid edge problems.

    Modified 2/23/2013 to properly treat RFPs where there is a reversal surface
    (q -> 0, iota -> infinity). This is handled through options controlled by
    the logical variable lrfp. lrfp is read from the first line of tae_data_boozer,
    which is written by xmetric. For analysis of systems with lrfp = .true. both
    updated versions of metric_element_create.f and stellagap.f are required.
    If using an older tae_data_boozer file, a read error will result since the
    first line will not have the lrfp logical. One can edit tae_data_boozer
    and add "F" in the first line to resolve this. However, to treat cases with
    lrfp = .true. it is necessary to have also used the upgraded xmetric. By
    default lrfp is set to .false.

    Modified 7/6/2012 making ir_fine_scl and irads to be command-line arguments
    (if not given, default values of 128 and 41 are used)

    Modified 12/2008 so that ion density profile, ion mass and several other
    parameters are read in through the plasma.dat file rather than being
    hardwired into the code. (see explanation below for variables contained
    in plasma.dat file).

    Modified on 9/18/2008 to offer both parallel or serial options, depending on
    use of the precompiler.

    This is the stellgap code for calculating Alfven contina in stellarators. The
    associated publication that describes this calculation is:

    D. A. Spong, R. Sanchez, A. Weller, "Shear Alfven Continua in Stellarators,"
    Phys. of Plasmas, Vol. 10, pg. 3217, August, 2003.

    This version of stellgap.f has been modified to be compatible with ae_solve.f.
    The frequencies are output in kHz. The mode selection is done using the fourier.dat
    file. The equilibrium data comes from the tae_data_boozer file. Following the
    stellgap run, a post-processing code, post_process.f must be run. The gap plot
    is then made using the results in the alfven_post file.
    The central ion density, ion density profile, and ion mass are set from the
    plasma.dat input file.
    Two other variables that the user should set are irads and ir_fine_scl. These
    are set througth paramter statements in the beginning of program tae_continua:
    integer, parameter :: ir_fine_scl = 1000
    integer, parameter :: irads = 39, irad3 = 3\*irads
    irads = # of flux surfaces in the initial equilibrium (wout file) - 2.
    The outermost and innermost surfaces from the initial equilibrium
    are left off to avoid inaccuracies that are sometimes present near
    these regions.
    ir_fine_scl = # of surfaces desired in the continuum output. Generally, to
    resolve fine scale features in the Alfven continua, more surfaces are
    desireable than in the original equilibrium data. The code performs
    interpolations to allow this.

    Input parameters: (fourier.dat file).

    ith, izt = number of theta and zeta grid points in the original mesh data
    (in file tae_data_vmec) calculated by the xmetric code.
    This file contains data
    for iota, phipc, |B|, gssup (the psi-psi contra-variant metric element,
    the Jacobian, and the contravariant poloidal and toroidal components of
    the magnetic field)
    nfp = number of field periods (for a tokamak run, set nfp = 1)
    mpol, ntor = number of poloidal and toroidal modes used in the Fourier
    expansion of the above gridded quantities suplied by the tae_data_vmec
    data file. To avoid anti-aliasing errors, mpol < 0.5*ith and
    ntor < 0.5*izt.
    mp_col = number of m,n mode pairs used in the representation of the
    eigebfunction
    nt_col = number of toroidal modes used in the representation of the
    eigenfunction (should be an integral multiple of nfp)
    mode_family = toroidal mode number about which the toroidal eigenfunction
    spectrum is built. i.e., for m = 0,
    n = mode_family, mode_family + nfp, mode_family + 2*nfp,
    ..., mode_family + nt_col
    while for m .ne. 0,
    n = -nt_col + mode_family, -nt_col + nfp + mode_family, -nt_col + 2*nfp + mode_family, ..., nt_col + mode_family

    Input parameters: (plasma.dat file).

    ion_to_proton_mass = m_ion/m_proton
    ion_density_0 = n_ion(0) = ion density (in m**-3) at magnetic axis
    ion_profile = integer parameter that determines form of n_ion(r)/n_ion(0) fit
    for ion_profile = 0 ion_density = [iota(rho)/iota(0)]**2
    for ion_profile = 1 ion_density = polynomial fit
    = nion(1) + nion(2)_rho + nion(3)_(rho**2) + nion(4)\*(rho**3) + ... + nion(9)_(rho\*\*8) + nion(10)_(rho**9)
    (here rho = normalized toroidal flux)
    for ion_profile = 2 ion_density = ion_density = constant = 1.
    for ion_profile = 3 ion_density = [1 - aion\*(rho**bion)]\*\*cion

    nion(10) = array of polynomial fit coefficients for ion_profile = 1
    aion, bion, cion = parameters for ion_profile = 3 fit
    jdqz_data = logical variable that is used in ae3d, but not in stellgap.f
    (included here only so that same plasma.dat file can be used)
    egnout_form = "binr" or "asci" - like jdqz_data this variable is used in
    ae3d, but not stellgap.f

    compile: sh bld_stellgap

    This code is can be run in serial or parallel mode, depending on whether
    is it run through the precompiler with SERIAL or PARALLEL set.
