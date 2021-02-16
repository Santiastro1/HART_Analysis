      CHARACTER*54      path1 ! you need to change the number of characters depending on the path length
      CHARACTER*63      path !  Path to analysis files
      character*37      gashead !
      character*98      filegas ! path1 + gashead + 7
      CHARACTER*45      HEADER
      CHARACTER*1 simclass,vopt

      PARAMETER (G= 6.672e-8) ! cgs units
      PARAMETER (pi = dacos(-1.d0))
      PARAMETER (deg = 2.d0*pi/360.d0)
      PARAMETER (densc=135.975) !Mo/Kpc (H=70Km/s/Mpc)
      PARAMETER (denscv=densc*93.0) 
      PARAMETER (densc200=densc*200.00)

      PARAMETER (Boxh  =20.0) ! box size Mpc/h, only for ART case

c If inputfiles.dat file is provided (Default), the next 2 parameters 
c are not needed

      PARAMETER (aini = 1.001) ! 1.000 initial scale factor for the analysis of ART/HART simulations
      PARAMETER (asteps = 0.200) ! time step for the analysis of ART/HART simulations


c      PARAMETER (numstep = 1)   ! number of snapshots to analyse
      
      PARAMETER (nbyteword=4)   ! defines length of direct-access record, 4 para HART, 1 para el resto

c      PARAMETER (NROW =1024, NGRID =128, ng=ngrid, NPAGE=NROW**2) !Usual in HART
      PARAMETER (NROW =2048, NGRID =128, ng=ngrid, NPAGE=NROW**2) !Usual in HART
   
      PARAMETER (NMAX=NGRID/2)
      PARAMETER (NRECL = NPAGE*6, NARR =MAX(2*NROW+1,NGRID+1))
      PARAMETER (NF67 = MAX(NROW,NGRID/2))
      PARAMETER (Nmaxpart = 5.2e7)
      parameter (nstarmax = 15.0e6) 
      PARAMETER (npartm = Nmaxpart)


      PARAMETER (path = 
     +'/Users/sroca/Postdoc/Simulations/MW_003/
     +RUN2.1_lowres/analysis/') ! 63
c     +'/Users/sroca/Postdoc/Simulations/MW_003/
c     +RUN2.1_highres/analysis/') ! 64
c     +'/Users/sroca/Postdoc/Simulations/VelaKlypin/
c     +astronomy.nmsu.edu/aklypin/HYDRO/Vela4_11/analysis/') ! 95
c     +'/Users/sroca/Postdoc/Simulations/MW_003/RUN2.1/analysis/') ! 56
      PARAMETER (path1 = 
     +'/Users/sroca/Postdoc/Simulations/MW_003/
     +RUN2.1_lowres/') ! 54
c     +'/Users/sroca/Postdoc/Simulations/MW_003/
c     +RUN2.1_highres/') ! 55
c     +'/Users/sroca/Postdoc/Simulations/VelaKlypin/
c     +astronomy.nmsu.edu/aklypin/HYDRO/Vela4_11/') ! 86
c     +'/Users/sroca/Postdoc/Simulations/MW_003/RUN2.2/') ! 47

      PARAMETER (steparm =  0.2 )     ! dR in kpc
      PARAMETER (rmax = 25.0) ! Last R to analyse in Kpc     
      PARAMETER (rmin = 0.1) ! First R to analyse, in some analysis is 0 by default in Kpc
      PARAMETER (nbins=int((rmax-rmin)/steparm)) ! Number of radial bins
      PARAMETER (narmbinsmax = 1000)
      PARAMETER (dZ_max  = 3.5)     ! max dZ for search (Kpc)
      PARAMETER (Rvir = 230.1) ! Rvir (kpc)

c  HART variables

      integer ncell0 , ng2 , narrg , nf67g ,  nctot
      integer nneigh , neighb , nchild , moct

      PARAMETER( gashead=
c     +'20MpcBox_HartGalMWgood_a') ! 24
c     +'20MpcBox_HartGalMW2012_a') !24
c     +'20MpcBox_HartGalMW003_RUN2_a') !28
c     +'20MpcBox_HartGalMW003_RUN2.1_a') !30
c     +'20MpcBox_HartGalMW003_RUN2.2_a') !30
c     +'20MpcBox_HartGalMW003_RUN2.3_a') !30
     +'20MpcBox_HartGalMW003_RUN2.1_lowres_a') !37
c     +'20MpcBox_HartGalMW2012_6_a') !26
c     +'Vela4_11.1_a') ! 12
      parameter(nchmax= 100 000 000)
      parameter(nbinmax= 200 )
      parameter ( ndim     = 3          ) ! # of spatial dimensions
      parameter ( mcell    = 100 000 000  ) ! # of refinement cells
      parameter ( nrowg    = 128       ) !# cells lev 0 gas
      parameter ( MinLevel = 0          ) ! minimum allowed level
      parameter ( MaxLevel = 11         ) ! maximum allowed level
      parameter ( nchem    = 2          ) ! # of chemical species
      parameter ( nvar     = 3          ) ! # of grav. variables(rho,phi1,phi2)
      parameter ( nhvar    = 8 + nchem  ) ! # of hydro variables
      parameter ( nspec    = 6          ) ! # of particle species 1, 2, 3, 6 (dm+1)
      parameter ( ncell0 = ng**3        ) ! # of zero level cells
      parameter ( nclmax = mcell       ) ! # max # of level cells  (>=ncell0)
      parameter ( ng2    = ng**2        ) ! # of cells in a grid layer
      parameter ( narrg   = ng + 1       ) ! FFT parameter
      parameter ( nf67g   = ng/2         ) ! FFT parameter
      double precision xng , yng
      parameter ( xng     = 1.d0 + 1.d0*ng) ! boundaries
      parameter ( yng     = 1.d0*ng - 1.0d-6    )
      parameter ( nctot  = mcell ) ! total number of cells
      parameter ( nneigh = 2*ndim       ) ! # of neighbors
      parameter ( neighb = 2*ndim       ) ! # of neighbors
      parameter ( nchild = 2**ndim      ) ! # of children
      parameter ( moct   = (mcell-ncell0)/nchild ) ! # of octs (old nky)

      parameter (optbigfile = 0 ) ! # If set to 1 only DM sp 1 will be readed

c Parameters when finding the center

      parameter (optcent=0 ) ! 0 if you look for the main resimulated galaxy,
c                            ! 1 if lookin for other system.
c     if optcent=1 the following initial values are requied
      parameter (xcentini=8.854) ! in Mpc 
      parameter (ycentini=25.635) ! in Mpc
      parameter (zcentini=0.8702) ! in Mpc
 
c Parameters when finding the galactic plane

      parameter (tolerance= 1d-6)
      parameter (maxiter= 100)


c Parameters for disk analysis as function of time

      PARAMETER (trmin = 5.0 ) ! Min rad for computing mean values of disk variables at t (in Kpc)
      PARAMETER (trmax = 6.0 ) ! Max rad for computing mean values of disk variables, at t (in Kpc)


c Parameters for bar analysis

      PARAMETER (angLimit = 7.) ! limit on deviation of bar position angle
      PARAMETER (bar_max = 8.0)      ! max bar radius in Kpc
      PARAMETER (bar_min = 0.5e-30)      ! min bar radius in Kpc

c  Select one or more stellar components

      PARAMETER (compo = 1) ! 0=gas; stars: 1=all, 2=disk, 3=bulge, 4=halo, 5=thin disk, 6=thick disk, 7=disk+bulge, 8=DM

c Necesary for vradrd computation (spherical or cylindrical shells)

      parameter (shell = 1) !1 for cyl. shells (disk), 2 for sph. (halo) (if shell = 2 select compo=8)

c Necesary for vcirc computation (binning)

      parameter( rinf= -1.0 )    ! galaxy halo
      parameter( width=20.0 )
      parameter( nbmax= 500 )

c Necesary for grid analysis computations

      PARAMETER( stepgrid=0.5 ) !grid smallest cell side length, in Kpc
      PARAMETER( ngridxy=150 ) !number of "half" grid elements in 1-D, x,y directions
      PARAMETER( ngridz=10 ) !number of "half" grid elements in 1-D, z direction
c                           !We use "half" grid elements because of symmetry between + and - axes              

      PARAMETER( steptheta=pi/60.) !Only used in the polar computation
      PARAMETER( nsup=1 ) !superposition factor. 2*nsup+1 consecutive cells will be superposed
c                         !except in the borders where the superpositon is only of nsup+1 cells


c Necesary for lv computation

c      PARAMETER( Vopt='V' ) !T if you want tangential velocity (without substractinb VC), V if you want galactocentric velocity




c
c     Needed variables but without changes to apply when changing the simulation
c
c

      DIMENSION RECDAT(NRECL),wspecies(10),lspecies(10)
      COMMON / CONTROL/ AEXPN,AEXP0,AMPLT,ASTEP,ISTEP,PARTW,
     +                  TINTG,EKIN,EKIN1,EKIN2,AU0,AEU0,
     +                  NROWC,NGRIDC,Nspecies,Nseed,Om0,Oml0,hubble,Wp5
     +                  ,Ocurv,extras(100)
      COMMON / ROW /    XPAR(NPAGE),YPAR(NPAGE),ZPAR(NPAGE),
     +                  VX(NPAGE),VY(NPAGE),VZ(NPAGE)
      COMMON / BWEIG/ iWeight(Nmaxpart),RWeight(Nmaxpart)
      COMMON /FOURAR/Zf(NARR),Yf(NARR)
      COMMON /F67COM/
     +                 IBC,      IP,       ISL,     L1,     N2,
     +                 N3,       N4,        N7,
     +                 SI(NF67),    INDEX(NF67) 

      EQUIVALENCE    (RECDAT(1),XPAR(1))
     +                               ,(wspecies(1),extras(1)),
     +                               (lspecies(1),extras(11)),
     +                               (Id,extras(90)),
     +                               (Norb,extras(91)),
     +                               (Rsnfw,extras(92)),
     +                               (diskmass,extras(93)),
     +                               (halomass,extras(94)),
     +                               (Rdisk,extras(95)),
     +                               (Cnfw,extras(96)),
     +                               (Nbulbo,extras(97)),
     +                               (Qt,extras(98)),
     +                               (Rtrunc,extras(99)),
     +                               (Caja,extras(100))



      parameter ( nil     = 0        )     ! integer zero

      parameter ( MaxL1     = MaxLevel + 1           )
      parameter ( isize     = 2**MaxL1               )
      parameter ( d_x       = 1.0 / 2**MaxL1         ) ! differs from Lesha!
      parameter ( nshift    = nchild - 1             )
      parameter ( nbshift   = nshift - ncell0        ) ! big shift; see Tools
      parameter ( mbshift   = -nbshift               )
      parameter ( ncell01   = -ncell0 - 1            )

c
c     ..........oct information..............
c
      common /TREE01/ iOctPs(ndim,0:moct)       ! coordinates of Oct centers
      common /TREE02/ iOctNb(nneigh,0:moct)     ! neighbouring cells
      common /TREE03/ iOctPr(0:moct)            ! parents/linked list index
      common /TREE04/ iOctLv(0:moct)            ! Level
      common /TREE05/ iOctCh(0:nctot)           ! children Octs

c
c    ........linked list of octs............
c
      common /LIST01/ iOctLL1(0:moct)              ! doubly linked list of octs
      common /LIST02/ iOctLL2(0:moct)
      common /LIST03/ iHOLL(MinLevel:MaxLevel+1)  ! linked list header
      common /LIST04/ iNOLL(MinLevel:MaxLevel+1)  ! # of ll entries at a Level

      real wmu
      parameter ( Y_p = 0.245 ) ! He mass fraction
      parameter ( wmu = 4.0 / (8.0 - 5.0 * Y_p) ) ! mol weight


      real var(nvar,nctot)
      real hvar(nhvar,nctot)

c     --------------------------------
c     Control parameters and variables
c     --------------------------------
c
      parameter ( nextra = 2 )  ! number of additional parameters
      real extra(nextra)
      character*256 lextra(nextra)

      real*8 r0, rho0, v0, t0, T_0, P0, S_0, aM0, E_0
      real*8 AL_SD, AL_Comp

c     ....... starformation parameters ........
c
      real*8 alpha_SF, C_SFR, eps_SF, dtmin_SF, dm_star_min,
     &       rho_SF, rho_SF_fact, T_SF, a_IMF, aM_stl, aM_stu, aM_SNII,
     &       aM_SNIa1,aM_SNIa2, t_SNIa, t_SNIai, C_SNIa, RIaf, ejM_SNIa,
     &       E_51, t_fb, C_fb, C_fbIa, fmass_met, c0_ml, T0_ml
c
c     indices of hvar for metal products of SN type II & Ia
c
      integer izII, izIa
      parameter ( izII = 9 , izIa = 10 )
c
c
c     ..................time...................

      real*8 t, dtime0, dtmin
      integer istep2
      real*8 tl(MinLevel:MaxLevel), tlold(MinLevel:MaxLevel)
      real*8 dtl(MinLevel:MaxLevel), dtlold(MinLevel:MaxLevel)
      real*8 aexp(MinLevel:MaxLevel), aexpold(MinLevel:MaxLevel)
      integer iTimeBin(MinLevel:MaxLevel)
      integer iSO(MinLevel:MaxLevel)  ! sweep order

c
c
c -----------------cooling parameters-------------------------
      integer nlt, nld, nlz, nrs
      real*8 tlmin, tlmax, dlt
      real*8 dlmin, dlmax, dld
      real*8 Zlmin, Zlmax, dlz
      real*8 rsmin, rsmax, drs
      real*8 dlti, dldi, dlzi, drsi
      common / CLCOOL01 / tlmin, tlmax, dlt, nlt
      common / CLCOOL02 / dlmin, dlmax, dld, nld
      common / CLCOOL03 / Zlmin, Zlmax, dlz, nlz
      common / CLCOOL04 / rsmin, rsmax, drs, nrs
      common / CLCOOL05 / dlti, dldi, dlzi, drsi

      real*8 smallrate
      parameter ( smallrate = 1.d-30 )
      parameter ( nltmax = 71 )
      parameter ( nldmax = 11 )
      parameter ( nlzmax = 9  )
      parameter ( nrsmax = 21 )
      real*8 coolcl(nltmax,nldmax,nlzmax,nrsmax)
      real*8 ccl_rs(nltmax,nldmax,nlzmax)
      real*8 f_ion(nltmax,nldmax,nlzmax,nrsmax)
      common / CLCOOL05 / coolcl
      common / CLCOOL06 / ccl_rs
      common / CLCOOL07 / f_ion

