      module aorsa3din_mod

      implicit none

!     --------------------------------------------------------
!     Declarations and defaults for aorsa3d.in input variables
!     --------------------------------------------------------
      integer :: ndiste = 0               !-----ndist:  if (ndist .eq. 0) Maxwellian is used in sigmad_stix
                                          !-----        if (ndist .eq. 1) non-Maxwellian is used in sigmad_stix
      integer :: ndisti1 = 0
      integer :: ndisti2 = 0
      integer :: ndisti3 = 0

      integer :: nmodesx = 32             !-----number of modes used in the x direction
      integer :: nmodesy = 32             !-----number of modes used in the y direction
      integer :: nmodesphi = 8            !-----number of modes used in the phi direction
      integer :: nwdot = 0                !-----number of radial modes used in wdot and flow (fy) calculation

      integer :: lmax = 5                 !-----highest order Bessel function kept in plasma conductivity
      integer :: ibessel = 1              !-----flag determining whether or not to expand ion Bessel functions
c           if(ibessel.eq.1) Full Bessel functions and exponential are used
c           if(ibessel.eq.2) Exact 2nd order expansion from finite difference
c                            code is used
c           if(ibessel.eq.3) 2nd order Larmor radius expansion of Bessel
c                            functions is used with exponential = 1.0
      integer :: inu = 0                  !-----if(inu.eq.0) real collisions are left out
      integer :: iprint = 50              !-----output is printed every iprint grid points
      integer :: iexact = 1               !-----if(iexact.eq.1) full sixth order equation is solved
                                          !-----if(iexact.eq.0) approximate second order equation is solved
      integer :: iroot = 2                !-----decides which of the two fast wave roots to follow
      integer :: iequat = 1               !-----if(iequat.eq.1) complete equations are solved
                                          !-----if(iequat.eq.2) Fukayama's equations are solved
      integer :: iflag_gammab = 1         !-----flag controlling the sign of gammab (parallel gradient parameter)
c           if(iflag_gammab .eq. 0) parallel gradient of B is set to zero
c           if(iflag_gammab .eq. 1) parallel gradient of B is non-zero(default)
      integer :: nnodecx = 32             !-----number of radial mesh points used for wdot calculation
      integer :: nnodecy = 32             !-----number of vertical mesh points used for wdot calculation
      integer :: nnodecphi = 8

      real :: delta0 = 0.0000E+00         !-----numerical damping for Bernstein wave:  about 1.e-04 (dimensionless)
      real :: xwall = 0.0000E+00          !-----not used
      real :: xnwall = 0.0000E+00         !-----density of metal put on last mesh point
      real :: epszet = 1.0000E-07         !-----value of kx and ky at nx and ny = 0.
      real :: phistart = 0.0              !-----starting position in phi for VMEC magnetic field
      integer :: iqx = 4
c-----iqx:  if(iqx.eq.1) Vaclavik's  kinetic flux is used
c-----      if(iqx.eq.2) Romero's  kinetic flux is used
c-----      if(iqx.eq.3) Jaeger's  kinetic flux is used
c-----      if(iqx.eq.4) Batchelor's  kinetic flux is used (reduces to WKB)
      integer :: izfunc = 1
c-----izfunc: if(izfunc.eq.1) not used
c-----        if(izfunc.eq.2) not used
      integer :: iez = 0
c-----iez:  if(iez.eq.0) Ez is calculated from complete equation
c-----      if(iez.eq.1) Ez is set to zero

      real :: amu1 = 2.0000E+00           !-----ratio of majority ion to hydrogen ion mass
      real :: amu2 = 1.0000E+00           !-----ratio of minority ion to hydrogen ion mass
      real :: amu3 = 1.2000E+01           !-----ratio of third ion mass to hydrogen ion mass
      real :: z1 = 1.0000E+00             !-----ratio of majority ion charge to hydrogen ion charge
      real :: z2 = 1.0000E+00             !-----ratio of minority ion charge to hydrogen ion charge
      real :: z3 = 6.0000E+00             !-----ratio of third ion charge to hydrogen ion charge
      real :: eta = 4.5000E-01            !-----ratio of minority ion density to electron density
      real :: eta3 = 4.6600E-02           !-----ratio of third ion density to electron density
      real :: xn0 = 3.1100E+19            !-----electron density at x=0
      real :: xnlim = 0.0000E+00          !-----electron density in scrape-off region (x>aplasm)
      real :: te0 = 4.2900E+03            !-----central value of eletron temperature in eV
      real :: ti0 = 7.0700E+03            !-----central value of ion temperature in eV
      real :: telim = 0.0000E+00          !-----electron temperature in scrape-off region (x>aplasm)
      real :: tilim = 0.0000E+00          !-----ion temperature in scrape-off region (x>aplasm)
      integer :: iprofile = 1             !-----iprofile:  if (iprofile .eq. 1) generic profiles (Gaussian) (default)
                                          !-----           if (iprofile .eq. 2) generic profiles (parabolas)
                                          !-----           if (iprofile .eq. 3) fits of form (1 - rho**beta)**alpha
                                          !-----           if (iprofile .eq. 5) numerical profiles from namelist

      real :: alphan = 1.0
      real :: alphate = 1.0
      real :: alphati = 1.0
      real :: betan = 1.0
      real :: betate = 1.0
      real :: betati = 1.0

      real :: b0 = 2.0                    !-----value of magnetic field at x=0 in Tesla
      real :: q0 = 1.0                    !-----value of inverse rotational transform on axis
      real :: rt = 2.1                    !-----major radius of torus
      real :: ekappa = 1.5                !-----elongation
      real :: rwleft = .70                !-----major radius of the left conducting wall
      real :: rwright = 2.5               !-----major radius of the right conducting wall
      real :: awally = 7.0000E-01         !-----vertical location of the conducting wall
      real :: ymax = 0.0                  !-----radius in vertical (y) direction- in default it is set to awallx
      real :: aplasm = 7.0000E-01         !-----location of the plasma-scrape-off interface
      real :: alim = 100.0                !-----location of turning point in step function density function

      real :: grad = 0.0
c-----grad = 0.0 ignors gradients in Wdot (default)
c-----grad = 1.0 includes gradients in Wdot

      real :: flat = 0.0000E+00
c-----flat=0.0 gives parabolic profiles
c-----flat=1.0 gives flat profiles
      real :: b1rat = 7.0000E-01          !-----low field value in step function magnetic field (igeom=3)
      real :: b2rat = 1.3000E+00          !-----high field value in step function magnetic field (igeom=3)

      real :: xnurf = 3.2000E+07          !-----rf frequency in Hertz
      real :: prfin = 0.0                 !-----total applied RF power
      real :: curdnx = 0.0000E+00         !-----Amps/meter of toroidal length of antenna in the x direction
      real :: curdny = 1.0                !-----Amps/meter of toroidal length of antenna in the y direction
      real :: curdnz = 0.0000E+00         !-----Amps/meter of toroidal length of antenna in the z direction

      integer :: icurve = 0
      real :: rant = 1.6                  !-----major radius of antenna in meters
      real :: yant = .1                   !-----half height of antenna in meters
      real :: dpsiant0 = .02
      real :: theta_ant = 0.0             !-----poloidal angle at which the antenna sits (degrees)
c           if(theta_ant .eq. 0.) antenna is on low field side or right (default)
c           if(theta_ant .eq. 180.) antenna is on high field side or left
      real :: dthetant0 = 40.
      real :: dphiant0 = 5.0              !-----not used
      integer :: nstrap = 1
      integer :: nphiant = 24             !-----not used
      real :: strap_width = 0.1524
      real :: strap_separ = 0.4572
      real :: phase_diff = 180.0
      real :: amplt(2) = 1.0

      integer :: igeom = 2
c-----igeom: if(igeom.eq.1)Solovev flux surfaces
c-----       if(igeom.eq.2)Stellarator flux surfaces
c-----       if(igeom.eq.3)Stellarator (expanded) flux surfaces (right handed)
c-----       if(igeom.eq.4)Stellarator (expanded) flux surfaces (left  handed)
c-----       if(igeom.eq.5)Stellarator flux surfaces from VMEC

      integer :: nfp = 12
      integer :: mcap(100) = 0            !-----number of the field period calculated (0 .le. mcap .le. np-1)
      integer :: mcap_number = 1

      integer :: nstep  = 16              !-----determines steepness of step function magnetic field in option igeom=3
      integer :: nabs = 2                 !-----polynomial coefficient which determines slope of the absorber
      real :: xnuabs = 0.0000E+00         !-----magnitude of absorber used to stop wall reflections in benchmark case.
      real :: xbnch = 0.0000E+00          !-----abs(x)>xbnch is the artificial absorber region to stop wall reflections
                                          !     in benchmark case.
                                          !-----if (xbnch.eq.0.0) it is ignored.
      real :: xleft = -7.000E-01          !-----left boundary for energy integrals and outgoing energy flux
      real :: xright = 7.0000E-01         !-----right boundary for energy integrals and incoming energy flux
      integer :: iabsorb = 2
c-----if(iabsorb.eq.1)electron absorption from simple Landau formula
c-----if(iabsorb.eq.2)electron absorption from  .5 * real(J* dot E) i.e. ECH
      integer :: isigma = 1
c-----if(isigma.eq.0) cold plasma conductivity is used.
c-----if(isigma.eq.1) hot  plasma conductivity is used (default).
      integer :: nzfun = 1
c-----nzfun:  if(nzfun.eq.0) Simple Z function is used from ZFUN
c-----        if(nzfun.eq.1) Generalized Z function of Brambilla is used (default).
c-----        if(nzfun.eq.2) Z function of Smithe is used by do numerical integrals.
c-----        if(nzfun.eq.3) Z function table lookup of Smithe is used
      real :: qavg0 = 1.0                 !-----qavg0 is the rotational transform on axis

      real :: xnuomg = 0.0                !-----xnuomg is the collision rate used in hot and cold plasma dielectrics
      real :: xnuead = 0.0000E+00         !-----ad hoc collision frequency for electron in sec-1
      real :: xnu1ad = 0.0000E+00         !-----ad hoc collision frequency for majority ions in sec-1
      real :: xnu2ad = 0.0000E+00         !-----ad hoc collision frequency for minority ions in sec-1
      real :: xnu3ad = 0.0000E+00         !-----ad hoc collision frequency for 3rd ion species

      integer :: itemp = 0
c-----if(itemp.eq.0)use Gaussian temperature-finite at edge-use with nlim=0 at edge
c-----if(itemp.eq.1)use Fukuyama's profile for temperature with c=-2 and d=1 -use
c          with nlim=finite at edge
c-----if(itemp.eq.2)use Gaussian times 1-r**2/xant**2

      integer :: nfreqm = 1               !-----number of frequencies run
      real :: dfreq = 0.0000E+00          !-----frequency increment
      integer :: nkzm = 1                 !-----number of kz's run
      real :: dkz = 0.0000E+00            !-----kz increment
      integer :: idens = 0
c-----if(idens.eq.0)use parabolic density profile
c-----if(idens.eq.1)use D'Ippolito density profile
      real :: r0 = 2.3000E+00             !-----r0 is parameter r0 in D'Ippolito density profile
      real :: xnudip = 2.5000E+00         !-----xnudip is parameter nu in D'Ippolito density profile
      real :: adip = 0.0000E+00           !-----adip is parameter nu in D'Ippolito density profile
      real :: efold = 0.0000E+00          !-----efold is the number of e-foldings in density between
                                          !     plasma edge and awall = awallx
      real :: xdelta = 5.5000E-01         !-----center of Gaussian for numerical damping of IBW
      real :: wdelta = 0.0000E+00         !-----width of Gaussian for numerical damping of IBW
      real :: xdelt2 = -7.000E-02         !-----second center of Gaussian for numerical damping of IBW
      real :: wdelt2 = 0.0000E+00         !-----second width of Gaussian for numerical damping of IBW
      real :: zeffcd = 2.5000E+00         !-----Zeff for Ehst-Karney current drive calculation
      real :: rzoom1 = 0.0                !-----R on left side of zoomed plot
      real :: rzoom2 = 0.0                !-----R on right side of zoomed plot
      integer :: ibackground = 1          !-----ibackground controls color of plotting window and labels
c-----    For white background set ibackground = 0 (box is black)
c-----    For black background set ibackground = 1 (box is red)

      integer :: lhel = 2
      integer :: mhel = 12

      real :: acoil = 0.46

      integer :: kplot = 1
      integer :: idiag = 5
      integer :: jdiag = 5
      integer :: kdiag = 5
      integer :: nplot = 1
      integer :: mplot = 1
      integer :: nphiplot = 4

      integer :: iexpnd = 1

      real :: damping = 0.0               !  enhancement factor (default = 0.0) for the electron conductivity (sig3) 
      real :: xkperp_cutoff = 0.75        !  applied above the fractional value of xkperp (xkperp_cutoff) to 
                                          !  short out noise in E_parallel 

      real :: signbz = 1.0000E+00

      real :: psilim = 1.00
      real :: psiant = 0.95
      real :: psiplasm = 0.60

      integer :: nboundary = 1            !-----nboundary: if(nboundary .eq. 1)flux surface boundary (default)
                                          !-----           if(nboundary .eq. 0)square boundary

      integer :: upshift = 1              !-----upshift: if (upshift .ne.  0) upshift is turned on (default)
                                          !-----if (upshift .eq. -1) upshift is turned off for xkperp > xkperp_cutoff
                                          !-----if (upshift .eq.  0) upshift is turned off always

      integer :: nprow = 8
      integer :: npcol = 8

      integer :: nuper = 65
      integer :: nupar = 129
      integer :: nkperp = 201             !-----nkperp: number of kperp values used in Lee's interpolation version of the 
                                          !     non-Maxwellian sigma (default = 201: interpolates on 201 points)
                                          !     if (nkperp .eq. 0) there is no interpolation

      integer :: nzeta_wdot = 51          !-----nzeta_wdot:  if (nzeta_wdot .eq. 0) no wdot calculation
                                          !-----             if (nzeta_wdot .eq. 1) wdot is calculated without interpolation
                                          !-----             if (nzeta_wdot .ge. 2) wdot is calculated with interpolation
                                          !-----                 over nzeta_wdot grid points (default is 51)
      integer :: ftrap = 1                !-----ftrap = integer flag determining whether trapped particles effect current drive
                                          !        if(ftrap.eq.0) no trapped particles
                                          !        if(ftrap.ne.0) include trapped particles (default)

      integer :: i_write = 0              !-----i_write: if (i_write .eq. 0) 4-D ORBIT_RF file is NOT written (default)
                                          !-----         if (i_write .ne. 0) 4-D ORBIT_RF file IS written
      integer :: n_bin = 2


      namelist/aorsa3din/nmodesx, nmodesy, nmodesphi, nwdot, lmax,
     .    ibessel, mhel, lhel, nfp, mcap, mcap_number,
     .    acoil, kplot, nplot, mplot, nphiplot, idiag, jdiag, kdiag,
     .    dthetant0, dphiant0, dpsiant0,
     .    ti0, xnuead, xnu1ad, xnu2ad, rant, te0, yant,
     .    inu, iprint, iexact, delta0, xwall, xnwall,
     .    iroot, iequat, igeom, epszet, phistart,
     .    iqx, izfunc, iez, nprow, npcol,
     .    amu1, amu2, z1, z2, eta,
     .    b0, rt, awally, xnurf, aplasm, xnlim, signbz,
     .    xn0, flat, b1rat, b2rat, curdnx, curdny, curdnz,
     .    nstep, nabs, xnuabs, xbnch, xleft, xright,
     .    isigma, itemp, telim, tilim,
     .    nfreqm, dfreq, nkzm, dkz,
     .    idens, r0, xnudip, adip, efold,
     .    amu3, z3, eta3, xnu3ad,
     .    xdelta, wdelta, xdelt2, wdelt2, zeffcd,
     .    rzoom1, rzoom2, ibackground, iabsorb, q0, prfin,
     .    nzfun, alim, grad, qavg0, nnodecx, nnodecy, nnodecphi, ymax,
     .    alphan, alphate, alphati, betan, betate, betati,
     .    ekappa, rwleft, rwright,
     .    nphiant, iexpnd, icurve, xnuomg, psilim, psiant, psiplasm,
     .    nstrap, iflag_gammab, theta_ant, strap_width, strap_separ,
     .    phase_diff, amplt, damping, xkperp_cutoff, iprofile,
     .    nuper, nupar, nkperp, nzeta_wdot, i_write, n_bin, nboundary,
     .    upshift, ndiste, ndisti1, ndisti2, ndisti3, ftrap

      end module aorsa3din_mod
