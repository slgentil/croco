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
#undef  MILES           /* NBQ MILES Applications */ 
#undef  REGIONAL        /* REGIONAL Applications */
#define BASIN_EQ        /* Bassin rectangulaire à l'équateur */

#if defined BASIN_EQ



!                      Bassin rectangulaire à l''équateur  
!                      ==================================
!
!                       Parallelization
# undef OPENMP
# define MPI
# undef PARALLEL_FILES
# define XIOS

!                       Model dynamics
# define SOLVE3D
# define UV_ADV
# define UV_COR

!    MOMENTUM LATERAL advection-diffusion scheme (defaut UP3)

# undef  UV_HADV_C4        /* 4th-order centered lateral advection */
# undef  UV_HADV_C2        /* 2nd-order centered lateral advection */
# define  UV_HADV_UP5    /* 5th-order upstream lateral advection */
# undef  UV_HADV_C6    /* 6th-order centered lateral advection */
# undef  UV_HADV_WENO5     /* 5th-order WENOZ    lateral advection */

! Tracer lateral advection scheme (default 3rd-order upstream UP3)

# undef TS_HADV_RSUP3   /* Rotated-Split UP3   */
# define TS_HADV_UP5    /* 5th-order upstream  */
# undef TS_HADV_WENO5   /* 5th-order WENO      */

! Tracer vertical advection scheme (default 4th-order  AKIMA)

# define  TS_VADV_SPLINES     /* Splines vertical advection            */
# undef  TS_VADV_WENO5      /* 5th-order WENOZ vertical advection    */
# undef TS_VADV_COMPACT    /* Splines vertical advection            */

# undef  VADV_ADAPT_IMP /* Semi-implicit vertical advection TS and UV */

!                       Lateral Mixing
# undef UV_VIS2
# undef UV_VIS4
# undef TS_DIF2
# undef TS_DIF4
# undef SMAGORINSKY

!                       Grid configuration and Initialisation
# define ANA_GRID
# define DEGREE
# define ANA_INITIAL
# undef BV_VAR
!                       Forcing
# define ANA_SMFLUX
# define ANA_STFLUX
# define ANA_BTFLUX
# define NO_FRCFILE

!                       Equation of State
# undef SALINITY
# ifdef SALINITY
# define SALINITY_INI
# define ANA_SSFLUX
# define ANA_BSFLUX
# endif


#define VTFORC_INT
#ifdef  VTFORC_INT
! toutes les clefs de rappel et dissipation
#define SPONGE_LAYER
#define NORTH_SPONGE
#define WEST_SPONGE
#define TNU_MAT         /* allow tnu to be space dependant */
#define UVNU_MAT         /* allow uvnu to be space dependant */
#endif    /* VTFORC_INT */
!
/*# ifdef SPONGE_LAYER
# define ANA_STFLUX
# endif*/

!                       Vertical Mixing
# define  LMD_MIXING /* Activate KPP mixing */
# ifdef  LMD_MIXING
#  undef  ANA_VMIX
#  undef  BVF_MIXING
#  undef LMD_SKPP
#  undef LMD_BKPP
#  define LMD_RIMIX
#  define LMD_CONVEC
# endif


#endif /* END OF CONFIGURATION CHOICE */

#include "cppdefs_dev.h"
#include "set_global_definitions.h"

