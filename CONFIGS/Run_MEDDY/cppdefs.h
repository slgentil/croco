! $Id: cppdefs.h 1628 2015-01-10 13:53:00Z marchesiello $
!
!======================================================================
! ROMS_AGRIF is a branch of ROMS developped at IRD and INRIA, in France
! The two other branches from UCLA (Shchepetkin et al)
! and Rutgers University (Arango et al) are under MIT/X style license.
! ROMS_AGRIF specific routines (nesting) are under CeCILL-C license.
!
! ROMS_AGRIF website : http://www.romsagrif.org
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
#undef  GRAV_ADJ        /* Graviational Adjustment Example */
#undef  INNERSHELF      /* Inner Shelf Example */
#undef  RIVER           /* River run-off Example */
#undef  OVERFLOW        /* Graviational/Overflow Example */
#undef  SEAMOUNT        /* Seamount Example */
#undef  SHELFRONT       /* Shelf Front Example */
#undef  SOLITON         /* Equatorial Rossby Wave Example */
#undef  UPWELLING       /* Upwelling Example */
#undef  VORTEX          /* Baroclinic Vortex Example */
#undef  INTERNAL        /* Internal Tide Example */
#undef  IGW             /* COMODO Internal Tide Example */
#undef  JET             /* Baroclinic Jet Example */
#undef  RIP             /* Rip Current Test Case */
#undef  SHOREFACE       /* Shoreface Test Case on a Planar Beach */
#undef  SWASH           /* Swash Test Case on a Planar Beach */
#undef  THACKER         /* Thacker wetting-drying Example */
#undef  TANK            /* Tank Example */
! aponte (
#undef  JETN             /* Baroclinic Jet Example, new setup */
! aponte )
#define  GO_MEDDY        /* Meddy Example */
#undef REGIONAL        /* REGIONAL Applications */


/*
!
!==========================================================
!              IDEALIZED CONFIGURATIONS
!==========================================================
!
*/
/*
!
!==========================================================
!              JETN CONFIGURATION
!==========================================================
!
*/
! aponte jet (
#if defined JETN
/*
!   Baroclinic JET Example, new configuration: sept. 2013
!   Builds on the JET configuration and the roms_ird_fev2011 CANYON_A config.
!                       ========== === =======
*/
! perturbation of the initial jet
# define JET_PERTURB
! resolution (dx in km)
# define RESOLUTION 8
! stratification profiles
# define JET_CONFIG 1
! narrow channel for internal tide test case
# undef JET_NARROW
! type of wavemaker
# define WMAKER_TEST 5
! turn beta on off
# define BETAON
! turn jet_decay on off
# undef JET_DECAY

# define XIOS

# undef  FLOATS
# undef  FLOAT_DEBUG

! old but necessary options
# define  ANA_JET
# define ITIDE         /* turns itide on/off*/
# undef FLAT_OCEAN    /* flattens the stratification on/off*/

/* reads initial profiles from a file, define ANA_INITIAL */
# define ANA_INITIAL
/* one tracer eos_linear */
# if (1)
#    undef ANA_DENSITY_FILE
#    undef NONLIN_EOS
#    undef SALINITY
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
# undef TS_HADV_RSUP3
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
#  define ZONAL_NUDGING
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
# undef  SMAGORINSKY
# endif


!
! Itide - Open boundary conditions
# ifdef ITIDE
# define IWMAKER
# define VMODES
# endif

# define OBC_SOUTH
# define OBC_NORTH
# define ANA_BRY

! 0 = closed / 1 = specified / 2 = radiation
# define OBCS_BTCONFIG0
# define OBCS_BCCONFIG0
# define OBCN_BTCONFIG0
# define OBCN_BCCONFIG0

!!!!!!  BAROTROPIC obcs

!!! south boundary
!# define OBCS_BTCONFIG1
# ifdef OBCS_BTCONFIG0
! closed
#  undef OBC_SOUTH
#  define OBCS_Z_GRADIENT
#  define OBCS_M2N_SPECIFIED
#  define OBCS_M2T_GRADIENT
# elif defined OBCS_BTCONFIG1
! specified
#  define OBCS_M2N_SPECIFIED
#  define OBCS_M2T_SPECIFIED
#  define OBCS_Z_SPECIFIED
# elif defined OBCS_BTCONFIG2
! radiation (method of characteristics)
#  define OBCS_M2N_GRADIENT  /* shouldn't do much in fact */
#  define OBCS_M2T_GRADIENT
# endif

!!! north boundary
! closed:
!# define OBCN_BTCONFIG0
# ifdef OBCN_BTCONFIG0
! closed
#  undef OBC_NORTH
#  define OBCN_Z_GRADIENT
#  define OBCN_M2N_SPECIFIED
#  define OBCN_M2T_GRADIENT
# elif defined OBCN_BTCONFIG2
! radiation (method of characteristics)
#  define OBCN_M2N_GRADIENT  /* shouldn't do much in fact */
#  define OBCN_M2T_GRADIENT
# endif


!!!!!! BAROCLINIC obcs

!!! south boundary
!# define OBCS_BCCONFIG1
# ifdef OBCS_BCCONFIG0
! closed
#  undef OBC_SOUTH
#  define OBCS_M3N_SPECIFIED
#  define OBCS_M3T_GRADIENT
#  define OBCS_T_GRADIENT
# elif defined OBCS_BCCONFIG1
! specified
#  define OBCS_T_SPECIFIED
#  define OBCS_M3N_SPECIFIED
#  define OBCS_M3T_SPECIFIED
# elif defined OBCS_BCCONFIG2
! radiation (method of characteristics)
#  define OBCS_SHIFT
#  define OBCS_VMODES_CHAR
#  undef OBCS_VMODES_CHAR_BTCLOSED  /* gradient for bt mode */
#  define OBCS_M3N_GRADIENT  /* shouldn't do much in fact */
#  define OBCS_M3T_GRADIENT
# endif

!!! north boundary
!# define OBCN_BCCONFIG0
# ifdef OBCN_BCCONFIG0
! closed:
#  undef OBC_NORTH
#  define OBCN_M3N_SPECIFIED
#  define OBCN_M3T_GRADIENT
#  define OBCN_T_GRADIENT
# elif defined OBCN_BCCONFIG2
! radiation (method of characteristics)
#  define OBCN_SHIFT
#  define OBCN_VMODES_CHAR
#  undef OBCN_VMODES_CHAR_BTCLOSED  /* gradient for bt mode */
#  define OBCN_M3N_GRADIENT  /* shouldn't do much in fact */
#  define OBCN_M3T_GRADIENT
# endif

# define FRC_BRY
# ifdef FRC_BRY
#  define Z_FRC_BRY
#  define M2_FRC_BRY
#  define M3_FRC_BRY
#  define T_FRC_BRY
# endif

!                    Diagnostics
# define TIDAL_DIAGS
!# define XIOS
# undef DIAGNOSTICS_TS
# undef DIAGNOSTICS_UV
# undef DIAGS_UV_SPEC /* not implemented */
# undef DIAGS_TS_SPEC /* not implemented */


! aponte jet )

#elif defined GO_MEDDY
/*
!
!==========================================================
!              GO_MEDDY CONFIGURATION
!==========================================================
!
*/
# undef OPENMP
# define MPI
# define UV_ADV
# define UV_COR
# define SOLVE3D
# define ANA_GRID
# define ANA_INITIAL
# define SALINITY
# undef  NBQ

!         Vertical Mixing
# undef ANA_VMIX
# undef BVF_MIXING
# define LMD_MIXING
# ifdef  LMD_MIXING
#  undef  ANA_VMIX
#  undef  BVF_MIXING
#  undef  LMD_SKPP
#  undef  LMD_BKPP
#  define LMD_RIMIX
#  define LMD_CONVEC
#  undef  LMD_KPP /* this flag seems deprecated !? */
# endif


!	   MOMENTUM LATERAL advection-diffusion scheme (defaut UP3)

# undef  UV_HADV_C4        /* 4th-order centered lateral advection */
# undef  UV_HADV_C2        /* 2nd-order centered lateral advection */
# undef  UV_HADV_UP5	   /* 5th-order upstream lateral advection */
# undef  UV_HADV_C6	   /* 6th-order centered lateral advection */
# undef  UV_HADV_WENO5	   /* 5th-order WENOZ    lateral advection */

!	Tracer lateral advection scheme (default 3rd-order upstream UP3)

# define TS_HADV_RSUP3		/* Rotated-Split UP3   */
# undef TS_HADV_UP5		/* 5th-order upstream  */
# undef TS_HADV_WENO5		/* 5th-order WENO      */

!	Momentum vertical advection scheme (default 4th-order Splines)
# undef UV_VADV_COMPACT	/* 4th-order compact scheme (a.k.a. Parabolic splines) */

!	Tracer vertical advection scheme (default 4th-order  AKIMA)

# define  TS_VADV_SPLINES   	/* Splines vertical advection            */
# undef  TS_VADV_WENO5     	/* 5th-order WENOZ vertical advection    */
# undef TS_VADV_COMPACT    /* Splines vertical advection            */

# undef  VADV_ADAPT_IMP	/* Semi-implicit vertical advection TS and UV */

!	Horizontal momentum dissipation

# undef UV_VIS2       		/*    Laplacian Dissipation             */
# undef  UV_VIS4       		/*    Biharmonic Dissipation            */

!	Horizontal Tracer diffusion

# undef TS_DIF2       		/*    Laplacian Diffusion             */
# undef TS_DIF4       		/*    Biharmonic Diffusion            */

!       Resolution = 64,128,200,256,400,512 or 1024
# define RESOLUTION 64


		/* analytical forcing */
# define ANA_SMFLUX
# define ANA_STFLUX
# define ANA_SSFLUX
# define ANA_BTFLUX
# define ANA_BSFLUX

		/* Boundary Conditions */
# define PERIODIC
# ifdef PERIODIC
#   define EW_PERIODIC
#   define NS_PERIODIC
# endif

# define PARALLEL_FILES
# define XIOS

 		/* Diagnostics */ 
# define  DIAG_SPEC             /* save spectral diagnostics */
# undef  DIAG_SPEC_KT             /* save spectral diagnostics */
# if defined DIAG_SPEC || defined DIAG_SPEC_KT
#   define  DIAGNOSTICS_UV
#   define  DIAGNOSTICS_TS
# endif


#endif /* END OF CONFIGURATION CHOICE */

#include "cppdefs_dev.h"
#include "set_global_definitions.h"

