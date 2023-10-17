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
#define  SLOPE        /* Meddy Example */


/*
!
!==========================================================
!              IDEALIZED CONFIGURATIONS
!==========================================================
!
*/

#if defined SLOPE
/*
!
!==========================================================
!              SLOPE CONFIGURATION
!==========================================================
!
*/
! to tag OSI parameter in output.mpi file
# define OSI 'OSI: '

# define JET_SLOPE    /* BELL, JET_SLOPE */

# undef OPENMP
# define MPI
# define NEW_S_COORD
# define UV_ADV
# define UV_COR
# define SOLVE3D
# define ANA_GRID
# define ANA_INITIAL
# undef SALINITY
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
#  undef LMD_RIMIX
#  undef LMD_CONVEC
#  undef  LMD_KPP /* this flag seems deprecated !? */
# endif


!	   MOMENTUM LATERAL advection-diffusion scheme (defaut UP3)

# undef  UV_HADV_C4        /* 4th-order centered lateral advection */
# undef  UV_HADV_C2        /* 2nd-order centered lateral advection */
# define  UV_HADV_UP5	   /* 5th-order upstream lateral advection */
# undef  UV_HADV_C6	   /* 6th-order centered lateral advection */
# undef  UV_HADV_WENO5	   /* 5th-order WENOZ    lateral advection */

!	Tracer lateral advection scheme (default 3rd-order upstream UP3)

# undef TS_HADV_RSUP3		/* Rotated-Split UP3   */
# undef TS_HADV_RSUP5		/* Pseudo R-Split UP5 lateral advection */
# define TS_HADV_UP5		/* 5th-order upstream  */
# undef TS_HADV_C6
# undef TS_HADV_WENO5		/* 5th-order WENO      */

!	Momentum vertical advection scheme (default 4th-order Splines)

# undef UV_VADV_WENO5   /* 5th-order  WENOZ   vertical advection  */
# undef UV_VADV_COMPACT	/* 4th-order compact scheme (a.k.a. Parabolic splines) */
# undef UV_VADV_C2    /* Splines vertical advection            */

!	Tracer vertical advection scheme (default 4th-order  AKIMA)

# define  TS_VADV_SPLINES   	/* Splines vertical advection            */
# undef  TS_VADV_WENO5     	/* 5th-order WENOZ vertical advection    */
# undef TS_VADV_COMPACT    /* Splines vertical advection            */
# undef TS_VADV_C2    /* Splines vertical advection            */

# undef  VADV_ADAPT_IMP	/* Semi-implicit vertical advection TS and UV */

!	Horizontal momentum dissipation

# undef UV_VIS2       		/*    Laplacian Dissipation             */
# undef  UV_VIS4       		/*    Biharmonic Dissipation            */

!	Horizontal Tracer diffusion

# undef TS_DIF2       		/*    Laplacian Diffusion             */
# undef TS_DIF4       		/*    Biharmonic Diffusion            */



		/* analytical forcing */
# define ANA_SMFLUX
# define ANA_STFLUX
# define ANA_SSFLUX
# define ANA_BTFLUX
# define ANA_BSFLUX
# define NO_FRCFILE

		/* Boundary Conditions */

# ifdef JET_SLOPE
#   undef OBC_EAST
#   undef OBC_WEST
#   define EW_PERIODIC
# else  /* BELL */
#   undef OBC_NORTH
#   undef OBC_SOUTH
#   define NS_PERIODIC
# endif
/*
# if defined OBC_EAST || defined OBC_WEST
#   undef  OBC_M2SPECIFIED
#   undef  OBC_M2FLATHER
#   define OBC_M2CHARACT
#   undef  OBC_M2ORLANSKI
#   define OBC_M3ORLANSKI
#   define OBC_TORLANSKI
#   undef  OBC_M3SPECIFIED
#   undef  OBC_TSPECIFIED
# endif
*/

# undef SPONGE
# undef FRC_BRY
# undef ANA_BRY
# undef Z_FRC_BRY
# undef M2_FRC_BRY
# undef M3_FRC_BRY
# undef T_FRC_BRY

# undef PARALLEL_FILES
# define XIOS
# define USE_CALENDAR

 		/* Diagnostics */ 
# undef  DIAG_SPEC             /* save spectral diagnostics */
# undef  DIAG_SPEC_KT             /* save spectral diagnostics */
# undef DIAG_SPEC_ISORHO
# if defined DIAG_SPEC || defined DIAG_SPEC_KT
#   define  DIAGNOSTICS_UV
#   define  DIAGNOSTICS_TS
# endif


#endif /* END OF CONFIGURATION CHOICE */

#include "cppdefs_dev.h"
#include "set_global_definitions.h"

