! $Id: cppdefs.h 1628 2015-01-10 13:53:00Z marchesiello $
!
!======================================================================
! CROCO is a branch of ROMS developped at IRD and INRIA, in France
! The two other branches from UCLA (Shchepetkin et al) 
! and Rutgers University (Arango et al) are under MIT/X style license.
! CROCO specific routines (nesting) are under CeCILL-C license.
! 
! CROCO website : http://www.croco-ocean.org
!======================================================================
!
/*
   This is "cppdefs.h": MODEL CONFIGURATION FILE
   ==== == ============ ===== ============= ====
*/
#undef  BASIN           /* Basin Example */
#undef  CANYON_A        /* Canyon_A Example */
#undef  CANYON_B        /* Canyon_B Example */
#undef  EQUATOR         /* Equator Example  */
#undef  INNERSHELF      /* Inner Shelf Example */
#undef  RIVER           /* River run-off Example */
#undef  OVERFLOW        /* Graviational/Overflow Example */
#undef  SEAMOUNT        /* Seamount Example */
#undef  SHELFRONT       /* Shelf Front Example */
#undef  SOLITON         /* Equatorial Rossby Wave Example */
#undef  THACKER         /* Thacker wetting-drying Example */
#undef  UPWELLING       /* Upwelling Example */
#undef  VORTEX          /* Baroclinic Vortex Example */
#undef  INTERNAL        /* Internal Tide Example */
#undef  IGW             /* COMODO Internal Tide Example */
#undef  JET             /* Baroclinic Jet Example */
#undef  SHOREFACE       /* Shoreface Test Case on a Planar Beach */
#undef  RIP             /* Rip Current Test Case */
#undef  FLUME           /* Bar-generating Flume Example */
#undef  SWASH           /* Swash Test Case on a Planar Beach */
#undef  TANK            /* Tank Example */
#undef  ACOUSTIC        /* Acoustic wave test case */
#undef  GRAV_ADJ        /* Graviational Adjustment Example */
#undef  KH_INST         /* Kelvin-Helmholtz Instability Example */
#undef  S2DV            /* S2DV sections */ 
#undef  MILES            /* NBQ MILES Applications */ 
#undef REGIONAL        /* REGIONAL Applications */

! aponte (
#define  ITBALT             /* Baroclinic Jet Example, new setup */
! aponte )

! aponte jet (
#if defined ITBALT
/*
!   Baroclinic JET Example, new configuration: sept. 2013
!   Builds on the JET configuration and the roms_ird_fev2011 CANYON_A config.
!                       ========== === =======
*/
! to tag OSI parameter in output.mpi file
# define OSI 'OSI: '
! perturbation of the initial jet
# undef JET_PERTURB
! resolution (dx in km)
# define RESOLUTION 4
! stratification profiles
# define JET_CONFIG 1
! narrow channel for internal tide test case
# undef JET_NARROW
! type of wavemaker
# define WMAKER_TEST 5
! turn beta on off
# undef BETAON
! turn jet_decay on off
# undef JET_DECAY
# undef JET_EDDY
# undef ITIDEP
# define FSTURB
# ifdef FSTURB
#   define VMODES
# endif

# define XIOS
# undef FERMI_MEM

# undef  FLOATS
# undef  SIGMA_FLOATS
# undef DEPTH_FLOATS
# undef  FLOAT_DEBUG

! old but necessary options
# define  ANA_JET
# undef ITIDE         /* turns itide on/off*/
# undef FLAT_OCEAN    /* flattens the stratification on/off*/

/* reads initial profiles from a file, define ANA_INITIAL */
# define ANA_INITIAL
/* one tracer eos_linear */
# if (1)
#    undef ANA_DENSITY_FILE
#    undef NONLIN_EOS 
# else
/* two tracer non linear eos from file */
#    define ANA_DENSITY_FILE
#    define NONLIN_EOS 
#    define SALINITY
# endif

!
# undef  ETALON_CHECK
# define  MPI
# define PARALLEL_FILES

!                   General flags
# define SOLVE3D
# define UV_COR
# define UV_ADV

!                   Grid and periodicity, bdy conditions
# define EW_PERIODIC
# ifdef ANA_JET
#  define ANA_GRID
#  ifdef ANA_GRID
#    undef ANA_DEPTH_FILE
#  endif
# endif
# undef Y_STRETCH_GRID
# undef Y_STRETCH_RELAX

!                   Horizontal advection for tracer (default UP3)
# define  TS_HADV_UP5

!                   Surface and bottom fluxes
# define ANA_STFLUX
# define ANA_BTFLUX
# define ANA_SMFLUX
# ifdef SALINITY
#   define ANA_BSFLUX
#   define ANA_SSFLUX
# endif

!                   Climatology and nudging
# define CLIMATOLOGY
# ifdef CLIMATOLOGY
#  undef ZONAL_NUDGING
#  define M2CLIMATOLOGY
#  define M3CLIMATOLOGY
#  define TCLIMATOLOGY
#  define ZCLIMATOLOGY
#  ifdef ANA_JET
#   define ANA_M2CLIMA
#   define ANA_M3CLIMA
#   define ANA_TCLIMA
#   define ANA_SSH
#  endif
#  define M3NUDGING
#  define M2NUDGING
#  define TNUDGING
#  define ZNUDGING
#  undef ROBUST_DIAG  /* defined in JET */
#  define SPONGE
# endif
!                   Vertical mixing
# define ANA_VMIX
# define LMD_MIXING
# ifdef  LMD_MIXING
#  undef  ANA_VMIX
#  define ANA_SRFLUX
#  define   LMD_SKPP
#  undef    LMD_BKPP
#  define LMD_RIMIX
#  define LMD_CONVEC
#  undef  LMD_KPP /* this flag seems deprecated !? */
# endif
!
!                   Lateral mixing
# if (false)
#  define UV_VIS2
#  define TS_DIF2
#  undef UV_VIS4
#  undef TS_DIF4
# else
#  undef UV_VIS2
#  undef TS_DIF2
#  define UV_VIS4
#  define TS_DIF4
# endif
# define UV_MIX_S
# define TS_MIX_S
# undef SMAGORINSKY
# ifdef  UV_VIS4
#   undef  SMAGORINSKY
# endif


!
! Itide - Open boundary conditions
# ifdef ITIDE
#   define IWMAKER
#   define VMODES
# endif

# define NS_PERIODIC
# undef OBC_SOUTH
# undef OBC_NORTH
# define ANA_BRY
#define ANA_BRY_WKB


# define FRC_BRY
# ifdef FRC_BRY
#  define Z_FRC_BRY
#  define M2_FRC_BRY
#  define M3_FRC_BRY
#  define T_FRC_BRY
# endif

!                    Diagnostics
# define TIDAL_DIAGS
# undef DIAGNOSTICS_TS
# undef DIAGNOSTICS_UV
# undef DIAGS_UV_SPEC /* not implemented */
# undef DIAGS_TS_SPEC /* not implemented */


! aponte jet )



#elif defined REGIONAL
/*
!====================================================================
!               REGIONAL (realistic) Configurations
!==================================================================== 
!
!----------------------
! BASIC OPTIONS
!----------------------
!
*/
                      /* Configuration Name */
# define BENGUELA_LR
                      /* Parallelization */
# undef  OPENMP
# undef  MPI
                      /* I/O server */
# undef  XIOS
                      /* Non-hydrostatic option */
# undef  NBQ
                      /* Nesting */
# undef  AGRIF
# undef  AGRIF_2WAY
                      /* OA and OW Coupling via OASIS (MPI) */
# undef  OA_COUPLING
# undef  OW_COUPLING
                      /* Wave-current interactions */
# undef  MRL_WCI
                      /* Open Boundary Conditions */
# undef  TIDES
# define OBC_EAST
# define OBC_WEST
# define OBC_NORTH
# define OBC_SOUTH
                      /* Applications */
# undef  BIOLOGY
# undef  FLOATS
# undef  STATIONS
# undef  PASSIVE_TRACER
# undef  SEDIMENT
# undef  BBL
                      /* dedicated croco.log file */
# undef  LOGFILE

/*!
!-------------------------------------------------
! PRE-SELECTED OPTIONS
!
! ADVANCED OPTIONS ARE IN CPPDEFS_DEV.H
!-------------------------------------------------
*/
                      /* Parallelization */
# ifdef MPI
#  undef  PARALLEL_FILES
# endif
# undef  AUTOTILING
# undef  ETALON_CHECK
                      /* Non-hydrostatic options */
# ifdef NBQ
#  define W_HADV_TVD
#  define W_VADV_TVD
# endif
                      /* Grid configuration */
# define CURVGRID
# define SPHERICAL
# define MASKING
# undef  WET_DRY
# define NEW_S_COORD
                      /* Model dynamics */
# define SOLVE3D
# define UV_COR
# define UV_ADV
                      /* Equation of State */
# define SALINITY
# define NONLIN_EOS
                      /* Lateral Momentum Advection (default UP3) */
# define UV_HADV_UP3
# undef  UV_HADV_UP5
# undef  UV_HADV_WENO5
# undef  UV_HADV_TVD
                      /* Lateral Explicit Momentum Mixing */
# undef  UV_VIS2
# ifdef UV_VIS2
#  define UV_VIS_SMAGO
# endif
                      /* Vertical Momentum Advection  */
# define UV_VADV_SPLINES
# undef  UV_VADV_WENO5
# undef  UV_VADV_TVD
                      /* Lateral Tracer Advection (default UP3) */
# undef  TS_HADV_UP3
# define TS_HADV_RSUP3
# undef  TS_HADV_UP5
# undef  TS_HADV_WENO5
                      /* Lateral Explicit Tracer Mixing */
# undef  TS_DIF2
# undef  TS_DIF4
# undef  TS_MIX_S
                      /* Vertical Tracer Advection  */
# undef  TS_VADV_SPLINES
# define TS_VADV_AKIMA
# undef  TS_VADV_WENO5
                      /* Sponge layers for UV and TS */
# define SPONGE
                      /* Semi-implicit Vertical Tracer/Mom Advection */
# undef  VADV_ADAPT_IMP
                      /* Vertical Mixing */
# undef  BODYFORCE
# undef  BVF_MIXING
# define LMD_MIXING
# undef  GLS_MIXING
# undef  GLS_MIX2017  /* <--- Warning: under testing */

# ifdef LMD_MIXING
#  define LMD_SKPP
#  define LMD_BKPP
#  define LMD_RIMIX
#  define LMD_CONVEC
#  undef  LMD_DDMIX
#  define LMD_NONLOCAL
# endif
# ifdef GLS_MIXING
#  define GLS_KKL
#  undef  GLS_KOMEGA
#  undef  GLS_KEPSILON
#  undef  GLS_GEN
#  undef  KANTHA_CLAYSON
#  undef  CRAIG_BANNER
#  undef  CANUTO_A
#  undef  ZOS_HSIG
# endif
# ifdef GLS_MIX2017
#  undef  GLS_KOMEGA
#  define GLS_KEPSILON
#  undef  GLS_GEN
#  define CANUTO_A
#  undef  GibLau_78
#  undef  MelYam_82
#  undef  KanCla_94
#  undef  Luyten_96
#  undef  CANUTO_B 
#  undef  Cheng_02
# endif
                      /* Surface Forcing */
# undef BULK_FLUX
# ifdef BULK_FLUX
#  define BULK_FAIRALL
#  define BULK_LW
#  define BULK_EP
#  define BULK_SMFLUX
#  undef  SST_SKIN
#  undef  ANA_DIURNAL_SW
#  undef  ONLINE
#  undef  ERA_ECMWF
#  undef  RELATIVE_WIND
# else
#  define QCORRECTION
#  define SFLX_CORR
#  define ANA_DIURNAL_SW
# endif
                      /* Wave-current interactions */
# ifdef OW_COUPLING
#  define MRL_WCI
#  define BBL
# endif
# ifdef MRL_WCI
#  ifndef OW_COUPLING
#   define WAVE_OFFLINE
#   undef  WKB_WWAVE
#  endif
#  undef  WAVE_ROLLER
#  define WAVE_STREAMING
#  define WAVE_FRICTION
#  define WAVE_RAMP
#  ifdef WKB_WWAVE
#   undef  WKB_OBC_NORTH
#   undef  WKB_OBC_SOUTH
#   undef  WKB_OBC_WEST
#   undef  WKB_OBC_EAST
#  endif
# endif
                      /* Lateral Forcing */
# define CLIMATOLOGY
# ifdef CLIMATOLOGY
#  define ZCLIMATOLOGY
#  define M2CLIMATOLOGY
#  define M3CLIMATOLOGY
#  define TCLIMATOLOGY

#  define ZNUDGING
#  define M2NUDGING
#  define M3NUDGING
#  define TNUDGING
#  undef  ROBUST_DIAG
# endif

# undef  FRC_BRY
# ifdef FRC_BRY
#  define Z_FRC_BRY
#  define M2_FRC_BRY
#  define M3_FRC_BRY
#  define T_FRC_BRY
# endif
                      /* Bottom Forcing */
# define ANA_BSFLUX
# define ANA_BTFLUX
                      /* Point Sources - Rivers */
# undef PSOURCE
# undef PSOURCE_NCFILE
# ifdef PSOURCE_NCFILE                    
#  undef PSOURCE_NCFILE_TS
# endif
                      /* Open Boundary Conditions */
# ifdef TIDES
#  define SSH_TIDES
#  define UV_TIDES
#  define POT_TIDES
#  define TIDERAMP
# endif
# define OBC_M2CHARACT
# undef  OBC_M2ORLANSKI
# define OBC_M3ORLANSKI
# define OBC_TORLANSKI
# undef  OBC_M2SPECIFIED
# undef  OBC_M3SPECIFIED
# undef  OBC_TSPECIFIED
                      /* Input/Output & Diagnostics */
# define AVERAGES
# define AVERAGES_K
# undef  DIAGNOSTICS_TS
# undef  DIAGNOSTICS_UV
# ifdef DIAGNOSTICS_TS
#  undef  DIAGNOSTICS_TS_ADV
#  define DIAGNOSTICS_TS_MLD
# endif
/*
!           Applications:
!---------------------------------
! Biology, floats, Stations, 
! Passive tracer, Sediments, BBL
!---------------------------------
!
   Quasi-monotone lateral advection scheme (WENO5)
   for passive/biology/sediment tracers 
*/
# if defined PASSIVE_TRACER || defined BIOLOGY || defined SEDIMENT
#  define BIO_HADV_WENO5
# endif
                      /*   Choice of Biology models   */
# ifdef BIOLOGY
#  undef  PISCES
#  undef BIO_NChlPZD
#  undef  BIO_N2ChlPZD2
#  define BIO_BioEBUS
                      /*   Biology options    */
#  ifdef PISCES
#   undef  DIURNAL_INPUT_SRFLX
#   define key_pisces
#  endif
#  ifdef BIO_NChlPZD
#   define  OXYGEN
#  endif
#  ifdef BIO_BioEBUS
#   define NITROUS_OXIDE
#  endif
                      /*   Biology diagnostics    */
#  define DIAGNOSTICS_BIO
#  if defined DIAGNOSTICS_BIO && defined PISCES
#   define key_trc_diaadd
#   define key_trc_dia3d
#   define key_iomput
#  endif
# endif
                      /*   Lagrangian floats model    */
# ifdef FLOATS
#  undef  FLOATS_GLOBAL_ATTRIBUTES
#  undef  IBM
#  undef  RANDOM_WALK
#  ifdef RANDOM_WALK
#   define DIEL_MIGRATION
#   define RANDOM_VERTICAL
#   define RANDOM_HORIZONTAL
#  endif
# endif
                      /*   Stations recording    */
# ifdef STATIONS
#  define ALL_SIGMA
# endif
                      /*   Sediment dynamics model     */
# ifdef SEDIMENT
#  define ANA_SEDIMENT
#  undef  ANA_SPFLUX
#  undef  ANA_BPFLUX
# endif
/*
!
!==========================================================
!              IDEALIZED CONFIGURATIONS
!==========================================================
!
*/
#elif defined BASIN
/*
!                       Basin Example
!                       ===== =======
*/
# undef  OPENMP
# undef  MPI
# define UV_ADV
# define UV_COR
# define UV_VIS2
# define SOLVE3D
# define TS_DIF2
# define BODYFORCE
# define ANA_GRID
# define ANA_INITIAL
# define ANA_SMFLUX
# define ANA_STFLUX
# define ANA_BTFLUX
# define NO_FRCFILE

#elif defined CANYON_A
/*
!                       First Canyon Example
!                       ===== ====== =======
*/
# undef  OPENMP
# undef  MPI
# define UV_ADV
# define UV_COR
# define SOLVE3D
# define EW_PERIODIC
# define ANA_GRID
# define ANA_INITIAL
# define ANA_SMFLUX
# define ANA_STFLUX
# define ANA_BTFLUX
# define NO_FRCFILE

#elif defined CANYON_B
/*
!                       Second Canyon Example
!                       ====== ====== =======
*/
# undef  OPENMP
# undef  MPI
# define UV_ADV
# define UV_COR
# define SOLVE3D
# define EW_PERIODIC
# define ANA_GRID
# define ANA_INITIAL
# define ANA_SMFLUX
# define ANA_STFLUX
# define ANA_BTFLUX
# define ANA_VMIX
# define NO_FRCFILE

#elif defined EQUATOR
/*
!                       Equator Example
!                       ======= =======
! Boccaletti, G., R.C. Pacanowski, G.H. Philander and A.V. Fedorov, 2004,
! The Thermal Structure of the Upper Ocean, J.Phys.Oceanogr., 34, 888-902.
*/
# undef  OPENMP
# undef  MPI
# define UV_ADV
# define UV_COR
# define UV_VIS2
# define SOLVE3D
# define SALINITY
# define TS_DIF2
# define ANA_GRID
# define ANA_INITIAL
# define ANA_SMFLUX
# define ANA_STFLUX
# define ANA_SRFLUX
# define ANA_SSFLUX
# define ANA_BTFLUX
# define ANA_BSFLUX
# define QCORRECTION
# define ANA_SST
# define LMD_SKPP /* problem with MPI in Xi direction */
# define LMD_MIXING
# define LMD_RIMIX
# define LMD_CONVEC
# define NO_FRCFILE

#elif defined INNERSHELF
/*
!                       Inner Shelf Example
!                       ===== ===== =======
*/
# undef  OPENMP
# undef  MPI
# define INNERSHELF_EKMAN
# define INNERSHELF_APG
# define SOLVE3D
# define UV_COR
# define ANA_GRID
# define ANA_INITIAL
# define AVERAGES
# define ANA_SSFLUX
# define ANA_SRFLUX
# define ANA_STFLUX
# define ANA_BSFLUX
# define ANA_BTFLUX
# define ANA_SMFLUX
# define NS_PERIODIC
# define OBC_WEST
# define SPONGE
# ifndef INNERSHELF_EKMAN
#  define UV_ADV
#  define SALINITY
#  define NONLIN_EOS
#  define LMD_MIXING
#  ifdef LMD_MIXING
#   define LMD_SKPP
#   define LMD_BKPP
#   define LMD_RIMIX
#   define LMD_CONVEC
#  else
#   define GLS_MIXING
#   define GLS_KKL
#  endif
# endif
# define NO_FRCFILE

#elif defined INTERNAL
/*
!                       Internal Tide Example
!                       ======== ==== =======
!
! Di Lorenzo, E, W.R. Young and S.L. Smith, 2006, Numerical and anlytical estimates of M2
! tidal conversion at steep oceanic ridges, J. Phys. Oceanogr., 36, 1072-1084.  
*/
# undef  OPENMP
# undef  MPI
# define SOLVE3D
# define UV_COR
# define UV_ADV
# define BODYTIDE
# define ANA_GRID
# undef  INTERNALSHELF
# define ANA_INITIAL
# define ANA_BTFLUX
# define ANA_SMFLUX
# define ANA_SRFLUX
# define ANA_STFLUX
# define ANA_VMIX
# define EW_PERIODIC
# define NS_PERIODIC
# undef  OBC_EAST
# undef  OBC_WEST
# undef  SPONGE
# undef  ANA_SSH
# undef  ANA_M2CLIMA
# undef  ANA_M3CLIMA
# undef  ANA_TCLIMA
# undef  ZCLIMATOLOGY
# undef  M2CLIMATOLOGY
# undef  M3CLIMATOLOGY
# undef  TCLIMATOLOGY
# undef  M2NUDGING
# undef  M3NUDGING
# undef  TNUDGING
# define NO_FRCFILE

#elif defined IGW
/*
!                  COMODO Internal Tide Example
!                  ====== ======== ==== =======
!
! Pichon, A., 2007: Tests academiques de maree, Rapport interne n 21 du 19 octobre 2007, 
! Service Hydrographique et Oceanographique de la Marine. 
*/

# define EXPERIMENT3
# undef  OPENMP
# undef  MPI
# undef  NBQ
# define NEW_S_COORD
# define TIDES
# define TIDERAMP
# define SSH_TIDES
# define UV_TIDES
# define SOLVE3D 
# define UV_ADV
# define UV_COR
# define UV_VIS2
# undef  VADV_ADAPT_IMP
# define SPHERICAL
# define CURVGRID
# define ANA_INITIAL
# define ANA_VMIX
# define ANA_SMFLUX
# define ANA_STFLUX
# define ANA_SRFLUX
# define ANA_SSFLUX
# define ANA_BTFLUX
# define ANA_BSFLUX
# define NS_PERIODIC
# define OBC_EAST
# define OBC_WEST
# undef  SPONGE
# define ANA_SSH
# define ANA_M2CLIMA
# define ANA_M3CLIMA
# define ANA_TCLIMA
# define ZCLIMATOLOGY
# define M2CLIMATOLOGY
# define M3CLIMATOLOGY
# define TCLIMATOLOGY
# define M2NUDGING
# define M3NUDGING
# define TNUDGING
# undef  ONLINE_ANALYSIS

#elif defined RIVER
/*
!                       River run-off test problem
!                       ==========================
*/
# undef  OPENMP
# undef  MPI
# define SOLVE3D
# define UV_ADV
# define UV_COR
# define M2FILTER_FLAT
# define NONLIN_EOS
# define SALINITY
# define ANA_GRID
# define MASKING
# define ANA_INITIAL
# define ANA_SMFLUX
# define ANA_STFLUX
# define ANA_SSFLUX
# define ANA_SRFLUX
# define ANA_BTFLUX
# define ANA_BSFLUX
# define LMD_MIXING
# define LMD_SKPP
# define LMD_BKPP
# define LMD_RIMIX
# define LMD_CONVEC
# define PSOURCE
# define ANA_PSOURCE
# define NS_PERIODIC
# define FLOATS
# ifdef FLOATS
#   define RANDOM_WALK
#   ifdef RANDOM_WALK
#      define DIEL_MIGRATION
#      define RANDOM_VERTICAL
#      define RANDOM_HORIZONTAL
#   endif
# endif
# define NO_FRCFILE

#elif defined SEAMOUNT
/*
!                       Seamount Example
!                       ======== =======
*/
# undef  OPENMP
# undef  MPI
# define UV_ADV
# define UV_COR
# define SOLVE3D
# define SALINITY
# define NONLIN_EOS
# define ANA_GRID
# define ANA_INITIAL
# define ANA_SMFLUX
# define ANA_STFLUX
# define ANA_SSFLUX
# define ANA_SRFLUX
# define ANA_BTFLUX
# define ANA_BSFLUX
# define NO_FRCFILE

# elif defined SHELFRONT
/*
!                       Shelf Front Example
!                       ===== ===== =======
*/
# undef  OPENMP
# undef  MPI
# define UV_ADV
# define UV_COR
# define SOLVE3D
# define SALINITY
# define ANA_GRID
# define ANA_INITIAL
# define ANA_SMFLUX
# define ANA_STFLUX
# define ANA_SSFLUX
# define ANA_SRFLUX
# define ANA_BTFLUX
# define ANA_BSFLUX
# define EW_PERIODIC
# define NO_FRCFILE

#elif defined SOLITON
/*
!                       Equatorial Rossby Wave Example
!                       ========== ====== ==== =======
*/
# undef  OPENMP
# undef  MPI
# define UV_COR
# define UV_ADV
# define ANA_GRID
# define ANA_INITIAL
# define AVERAGES
# define EW_PERIODIC
# define ANA_SMFLUX
# define NO_FRCFILE

#elif defined THACKER
/*
!                       Thacker Example
!                       ======= =======
!
! Thacker, W., (1981), Some exact solutions to the nonlinear 
! shallow-water wave equations. 
! J. Fluid Mech., 107, 499â508.
*/
# undef  OPENMP
# undef  MPI
# define THACKER_2DV
# define SOLVE3D
# define UV_COR
# define UV_ADV
# undef  UV_VIS2
# define WET_DRY
# define NEW_S_COORD
# define ANA_GRID
# define ANA_INITIAL
# define ANA_BTFLUX
# define ANA_SMFLUX
# define ANA_SRFLUX
# define ANA_STFLUX
# define NO_FRCFILE

# elif defined OVERFLOW
/*
!                       Gravitational/Overflow Example
!                       ====================== =======
*/
# undef  OPENMP
# undef  MPI
# define UV_ADV
# define UV_COR
# define UV_VIS2
# define TS_DIF2
# define TS_MIX_GEO
# define SOLVE3D
# define ANA_GRID
# define ANA_INITIAL
# define ANA_SMFLUX
# define ANA_STFLUX
# define ANA_BTFLUX
# define NO_FRCFILE
/*
!                       Plume Example
!                       ===== =======
*/
#elif defined PLUME
# undef  OPENMP
# undef  MPI
# define NEW_S_COORD
# define ANA_GRID
# define ANA_INITIAL
# define ANA_SMFLUX
# define ANA_STFLUX
# define ANA_SSFLUX
# define ANA_SRFLUX
# define ANA_BSFLUX
# define ANA_BTFLUX
# define SOLVE3D
# define UV_COR
# define UV_ADV
# define SALINITY
# define NONLIN_EOS
# define SPLIT_EOS
# define TS_HADV_UP3
# define SPONGE
# undef  LMD_MIXING
# define GLS_MIX2017
# ifdef LMD_MIXING
#  define LMD_SKPP
#  define LMD_BKPP
#  define LMD_RIMIX
#  define LMD_CONVEC
#  undef  LMD_DDMIX
#  define LMD_NONLOCAL
# endif
# ifdef GLS_MIX2017
#  undef  GLS_KOMEGA
#  define GLS_KEPSILON
#  undef  GLS_GEN
#  define CANUTO_A
#  undef  GibLau_78
#  undef  MelYam_82
#  undef  KanCla_94
#  undef  Luyten_96
#  undef  CANUTO_B
#  undef  Cheng_02
# endif
# define NO_FRCFILE

#elif defined UPWELLING
/*
!                       Upwelling Example
!                       ========= =======
*/
# undef  OPENMP
# undef  MPI
# define SOLVE3D
# define UV_COR
# define UV_ADV
# define ANA_GRID
# define ANA_INITIAL
# define AVERAGES
# define SALINITY
# define NONLIN_EOS
# define ANA_SSFLUX
# define ANA_SRFLUX
# define ANA_STFLUX
# define ANA_BSFLUX
# define ANA_BTFLUX
# define ANA_SMFLUX
# define LMD_MIXING
# define LMD_SKPP
# define LMD_BKPP
# define LMD_RIMIX
# define LMD_CONVEC
# define EW_PERIODIC
# define NO_FRCFILE

#elif defined VORTEX
/*
!                       Baroclinic Vortex Example (TEST AGRIF)
!                       ========== ====== ======= ===== ======
*/
# undef  OPENMP
# undef  MPI
# define AGRIF
# define AGRIF_2WAY
# undef  NBQ
# define SOLVE3D
# define UV_COR
# define UV_ADV
# define ANA_STFLUX
# define ANA_SMFLUX
# define ANA_BSFLUX
# define ANA_BTFLUX
# define ANA_VMIX
# define OBC_EAST
# define OBC_WEST
# define OBC_NORTH
# define OBC_SOUTH
# define SPONGE
# define ZCLIMATOLOGY
# define M2CLIMATOLOGY
# define M3CLIMATOLOGY
# define TCLIMATOLOGY
# define ZNUDGING
# define M2NUDGING
# define M3NUDGING
# define TNUDGING
# define NO_FRCFILE

#elif defined JET
/*
!                       Baroclinic JET Example
!                       ========== === =======
*/
# define ANA_JET
# define MPI
# undef  NBQ
# undef  AGRIF
# undef  AGRIF_2WAY
# define SOLVE3D
# define UV_COR
# define UV_ADV
# define UV_VIS2
# ifdef ANA_JET
#  define ANA_GRID
#  define ANA_INITIAL
# endif
# define ANA_STFLUX
# define ANA_SMFLUX
# define ANA_BSFLUX
# define ANA_BTFLUX
# define ANA_VMIX
# define EW_PERIODIC
# define CLIMATOLOGY
# ifdef CLIMATOLOGY
#  define ZCLIMATOLOGY
#  define M2CLIMATOLOGY
#  define M3CLIMATOLOGY
#  define TCLIMATOLOGY
#  define ZNUDGING
#  define M2NUDGING
#  define M3NUDGING
#  define TNUDGING
#  define ROBUST_DIAG
#  define ZONAL_NUDGING
#  ifdef ANA_JET
#   define ANA_SSH
#   define ANA_M2CLIMA
#   define ANA_M3CLIMA
#   define ANA_TCLIMA
#  endif
# endif
# define LMD_MIXING 
# ifdef  LMD_MIXING
#  undef  ANA_VMIX
#  define ANA_SRFLUX
#  undef  LMD_KPP
#  define LMD_RIMIX
#  define LMD_CONVEC
# endif 
# define NO_FRCFILE

#elif defined SHOREFACE
/*
!                       PLANAR BEACH Example
!                       ====== ===== =======
!
!   Uchiyama, Y., McWilliams, J.C. and Shchepetkin, A.F. (2010): 
!      Wave-current interaction in an oceanic circulation model with a 
!      vortex force formalism: Application to the surf zone.
!      Ocean Modelling Vol. 34:1-2, pp.16-35.
*/
# undef  OPENMP
# undef  MPI
# define SOLVE3D
# define UV_ADV
# undef  MASKING
# define WET_DRY
# define NEW_S_COORD
# define ANA_GRID
# define ANA_INITIAL
# define ANA_SMFLUX
# define ANA_STFLUX
# define ANA_SSFLUX
# define ANA_SRFLUX
# define ANA_SST
# define ANA_BTFLUX
# define NS_PERIODIC
# define OBC_WEST
# define SPONGE
# define MRL_WCI
# ifdef MRL_WCI
#  undef  WAVE_OFFLINE
#  ifndef WAVE_OFFLINE
#   define WKB_WWAVE
#   define WKB_OBC_WEST
#   define WAVE_FRICTION
#   undef  WAVE_ROLLER
#   undef  MRL_CEW
#  endif
# endif
# define LMD_MIXING
# define LMD_SKPP
# define LMD_BKPP
# undef  BBL
# undef  SEDIMENT
# ifdef SEDIMENT
#  define TCLIMATOLOGY
#  define TNUDGING
#  define ANA_TCLIMA
# endif

#elif defined FLUME
/*
!                       FLUME Example
!                       ===== =======
!
!   Roelvink, J.A. and Stive, M.J.F., 1989: Bar-generating cross-shore 
!       flow mechanisms on a beach. Journal of Geophysical Research
*/
# undef  OPENMP
# undef  MPI
# define SOLVE3D
# define UV_ADV
# define NEW_S_COORD
# define ANA_GRID
# define ANA_INITIAL
# define ANA_SMFLUX
# define ANA_STFLUX
# define ANA_SSFLUX
# define ANA_SRFLUX
# define ANA_SST
# define ANA_BTFLUX
# define NS_PERIODIC
# define OBC_WEST
# define SPONGE
# define WET_DRY
# define MRL_WCI
# ifdef MRL_WCI
#  define WKB_WWAVE
#  define MRL_CEW
#  define WKB_OBC_WEST
#  define WAVE_ROLLER
#  define WAVE_FRICTION
#  define WAVE_BREAK_SWASH
#  undef  WAVE_STREAMING
#  undef  WAVE_RAMP
# endif
# define LMD_MIXING
# define LMD_SKPP
# define LMD_BKPP
# define LMD_VMIX_SWASH
# undef  GLS_MIX2017
# ifdef GLS_MIX2017
#  undef  GLS_KOMEGA
#  define GLS_KEPSILON
#  undef  GLS_GEN
#  define CANUTO_A
#  undef GibLau_78
#  undef MelYam_82
#  undef KanCla_94
#  undef Luyten_96
#  undef CANUTO_B 
#  undef Cheng_02
# endif

# define BBL
# define SEDIMENT
# ifdef SEDIMENT
#  define TCLIMATOLOGY
#  define TNUDGING
#  define ANA_TCLIMA
#  define MORPHODYN
# endif
# define NO_FRCFILE

#elif defined RIP
/*
!                       Rip Current Example
!                       === ======= =======
!
!   Weir, B., Uchiyama, Y.. (2010): 
!      A vortex force analysis of the interaction of rip 
!      currents and surface gravity wave
!      JGR Vol. 116
!
!  Default is idealized Duck Beach with 3D topography
!  RIP_TOPO_2D: Logshore uniform topography
!  BISCA: realistic case with Input files
!  GRANDPOPO: idealized Grand Popo Beach in Benin, 
!              longshore uniform
*/
# undef  BISCA
# undef  RIP_TOPO_2D
# undef  GRANDPOPO
# ifdef GRANDPOPO
#  define RIP_TOPO_2D
# endif
# undef ANA_TIDES
!
# undef  OPENMP
# undef  MPI
# define SOLVE3D
# define NEW_S_COORD
# define UV_ADV
# undef  VADV_ADAPT_IMP
# define UV_VIS2
# define UV_VIS_SMAGO
# define LMD_MIXING
# define LMD_SKPP
# define LMD_BKPP
# ifndef BISCA
#  define ANA_GRID
# endif
# define ANA_INITIAL
# define ANA_SMFLUX
# define ANA_STFLUX
# define ANA_SSFLUX
# define ANA_SRFLUX
# define ANA_SST
# define ANA_BTFLUX
# if !defined BISCA && !defined ANA_TIDES
#  define NS_PERIODIC
# else
#  define OBC_NORTH
#  define OBC_SOUTH
# endif
# define OBC_WEST
# define SPONGE
# ifdef ANA_TIDES
#  define ANA_SSH
#  define ANA_M2CLIMA
#  define ANA_M3CLIMA
#  define ZCLIMATOLOGY
#  define M2CLIMATOLOGY
#  define M3CLIMATOLOGY
#  define M2NUDGING
#  define M3NUDGING
# endif
# define WET_DRY
# define MRL_WCI
# ifdef MRL_WCI
#  define WKB_WWAVE
#  define WKB_OBC_WEST
#  define WAVE_ROLLER
#  define WAVE_FRICTION
#  define WAVE_STREAMING
#  define MRL_CEW
#  ifdef RIP_TOPO_2D
#   define WAVE_RAMP
#  endif
# endif
# ifdef BISCA
#  define BBL
# endif
# undef SEDIMENT
# ifdef SEDIMENT
#  define ANA_SEDIMENT
#  undef  ANA_SPFLUX
#  undef  ANA_BPFLUX
# endif
# undef  DIAGNOSTICS_UV

#elif defined SWASH
/*
!                       SWASH PLANAR BEACH Example
!                       ===== ====== ===== =======
!
*/
# undef  OPENMP
# undef  MPI
# define SOLVE3D
# undef  NBQ
# define UV_ADV
# define UV_VIS2
# define UV_VIS_SMAGO
# define MASKING
# define NEW_S_COORD
# define ANA_GRID
# define ANA_INITIAL
# define ANA_SMFLUX
# define ANA_STFLUX
# define ANA_SSFLUX
# define ANA_SRFLUX
# define ANA_SST
# define ANA_BTFLUX
# define OBC_WEST
# define SPONGE
# define FRC_BRY
# define ANA_BRY
# define Z_FRC_BRY
# define M2_FRC_BRY
# define M3_FRC_BRY
# define T_FRC_BRY
# define WET_DRY
# define NO_FRCFILE

#elif defined TANK
/*
!                       Tank Example
!                       ==== =======
!
! Chen, X.J., 2003. A fully hydrodynamic model for three-dimensional,
! free-surface flows.
! Int. J. Numer. Methods Fluids 42, 929â~@~S952.
*/
# undef  TANKY
# undef  MPI
# define NBQ
# ifdef NBQ
#  undef  NBQ_PRECISE
#  define NBQ_PERF
# endif
# define SOLVE3D
# undef  UV_ADV
# define NEW_S_COORD
# define ANA_GRID
# define ANA_INITIAL
# define ANA_BTFLUX
# define ANA_SMFLUX
# define ANA_SRFLUX
# define ANA_STFLUX
# define NO_FRCFILE

#elif defined ACOUSTIC 
/*
!                       ACOUSTIC WAVE TESTCASE 
!                       ======== ==== ========
*/
# undef  MPI
# define NBQ
# ifdef NBQ
#  undef  NBQ_PRECISE
#  define NBQ_PERF
# endif
# undef  UV_VIS2
# define SOLVE3D
# define NEW_S_COORD
# define ANA_GRID
# define ANA_INITIAL
# define ANA_SMFLUX
# define ANA_STFLUX
# define ANA_BTFLUX
# define NO_FRCFILE

#elif defined GRAV_ADJ
/*
!                       Gravitational Adjustment Example
!                       ============= ========== =======
!
!  Soliton case GRAV_ADJ_SOLITON (non-hydro test) is setup from:
!  Horn, D.A., J. Imberger, & G.N. Ivey, (2001). 
!  The degeneration of large-scale interfacial gravity waves in lakes. 
!  J. Fluid Mech., 434:181-207. 
!
*/
# undef  OPENMP
# undef  MPI
# undef  NBQ
# undef  XIOS 
# ifdef NBQ
#  define GRAV_ADJ_SOLITON
#  undef  NBQ_PRECISE
#  define NBQ_PERF
# endif
# define UV_VIS2
# define SOLVE3D
# define NEW_S_COORD
# define UV_ADV
# define TS_HADV_WENO5
# define TS_VADV_WENO5
# define ANA_GRID
# define ANA_INITIAL
# define ANA_SMFLUX
# define ANA_STFLUX
# define ANA_BTFLUX
# undef PASSIVE_TRACER
# define NO_FRCFILE

#elif defined KH_INST 
/*
!                       Kelvin-Helmholtz Instability Example
!                       ================ =========== =======
!
*/
# undef  KH_INST2D
# undef  KH_INST3D
# undef  KH_INSTY
# undef  MPI
# define NBQ
# ifdef NBQ
#  define NBQ_PERF
#  undef  NBQ_PRECISE
#  define NBQ_OBC
# endif
# undef  XIOS
# define SOLVE3D
# define NEW_S_COORD

# define UV_VIS2
# define UV_MIX_S
# define UV_VIS_SMAGO
# define SALINITY
# undef  PASSIVE_TRACER

# define UV_ADV
# define TS_HADV_WENO5
# define TS_VADV_WENO5
# define TS_DIF2
# define UV_VADV_C2
# define UV_HADV_C2
# undef  UV_VADV_TVD
# undef  UV_HADV_TVD
# undef  W_VADV_TVD
# undef  W_HADV_TVD
# define RESET_RHO0

# define ANA_GRID
# define ANA_INITIAL
# define ANA_SMFLUX
# define ANA_STFLUX
# undef  ANA_SRFLUX
# define ANA_BTFLUX
# define ANA_SSFLUX
# define ANA_BSFLUX
# ifndef KH_INSTY
#  define EW_PERIODIC
# else
#  define NS_PERIODIC
# endif

#elif defined S2DV 
/*
!                  2DV Sections 
!                  ============
*/
# undef  EXPERIMENT1
# define NBQ
# define MPI
# define XIOS
# define SPHERICAL
# define CURVGRID
# undef  MASKING
# define NEW_S_COORD
# undef  WET_DRY
# ifdef NBQ
#  undef  NBQ_PRECISE
#  define NBQ_PERF
# endif
# define SOLVE3D 
# define UV_ADV
# define UV_COR
# define SALINITY
# undef  VADV_ADAPT_IMP
# define SUPERBEE 
# define UV_HADV_TVD
# define UV_VADV_TVD
# undef  UV_VADV_WENO5
# undef  UV_HADV_WENO5
# define W_HADV_TVD
# define W_VADV_TVD
# define TS_VADV_WENO5
# define TS_HADV_WENO5
# define UV_VIS2
# define UV_VIS_SMAGO
# undef  TIDES
# ifdef TIDES
#  define TIDERAMP
#  define SSH_TIDES
#  define UV_TIDES 
# endif
# define OBC_EAST
# define OBC_WEST
# define NS_PERIODIC
# undef  ANA_INITIAL
# undef  ANA_VMIX
# define ANA_SMFLUX
# define ANA_STFLUX
# define ANA_SRFLUX 
# define ANA_SSFLUX
# define ANA_BTFLUX
# define ANA_BSFLUX
# undef  PASSIVE_TRACER

#elif defined MILES 
/*
!                NBQ MILES APPLICATIONS 
!                === ===== ============
!
*/
# undef  EXPERIMENT1
# undef  EXPERIMENT3
# define NBQ
# define MPI
# define XIOS
# define SOLVE3D 
# define UV_ADV
# define UV_COR
# ifdef NBQ
#  undef  NBQ_PRECISE
#  define NBQ_PERF
#  define NBQ_FREESLIP
# endif
# define SPHERICAL
# define CURVGRID
# define MASKING
# define NEW_S_COORD
# undef  WET_DRY
# undef  VADV_ADAPT_IMP
# undef  UV_VADV_C2
# undef  UV_HADV_C2
# define UV_VADV_TVD
# define UV_HADV_TVD
# undef  SUPERBEE
# define W_VADV_TVD
# define W_HADV_TVD
# define UV_VIS2
# undef  UV_VIS_SMAGO
# define SALINITY
# define TS_HADV_WENO5
# define TS_VADV_WENO5
# undef  PASSIVE_TRACER
# undef  TIDES
# ifdef TIDES
#  undef  TIDERAMP
#  define SSH_TIDES
#  define POT_TIDES
#  define UV_TIDES 
# endif
# define OBC_EAST
# define OBC_WEST
# undef  OBC_SOUTH
# undef  OBC_NORTH
# define SPONGE
# undef  LMD_MIXING
# undef  LMD_BKPP
# undef  GLS_MIX2017
# undef  ANA_INITIAL
# define ANA_SMFLUX
# define ANA_STFLUX
# define ANA_SRFLUX
# define ANA_SSFLUX
# define ANA_BTFLUX
# define ANA_BSFLUX


#endif /* END OF CONFIGURATION CHOICE */

#include "cppdefs_dev.h"
#include "set_global_definitions.h"

