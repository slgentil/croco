CROCO v1.0
-----------
Released date : 26 June 2018

Previous release : ROMS_AGRIF v3.1.1 (July 2014)

New in v1.0 :
=============
CROCO sources and CROCO_TOOLS (the follow-on of ROMS_TOOLS) are now distributed separately (for croco_tools releases, see associated tab at  https://www.croco-ocean.org/download/croco-project/ ).  

CROCO has now a new architecture tree. The OCEAN directory contains the sources. 

New in CROCO v1.0 (associated cppkey names and dedicated configuration in cppdefs.h are presented):

- Non-hydrostatic kernel: a (3D) non-Boussinesq “fast-mode” is available.  
#define NBQ 
Dedicated configurations : #define TANK, ACOUSTIC, GRAV_ADJ (NBQ), KH_INST, S2DV, MILES
croco.in: set CSOUND_NBQ (sound speed) to a minimum of 5 times the max external gravity wave phase speed (sqrt(gh)). Max is real sound speed of 1500 m/s. Lower sound speed allows larger time steps.
cppdefs.h: choice of the full scheme (NBQ_PRECISE) or for faster integration a simplified one (NBQ_PERF). In this case, the vertical grid and some non-dominant terms of the fast mode equations are only updated at the internal (slow) step. NBQ_PERF is the default option and is recommended for submesoscale applications (>10m resolution), while NBQ_PRECISE should be used for LES applications (e.g., KH_INST)
cppdefs.h: choice of advection scheme for vertical velocity w includes TVD and WENO (see below); default choice is W_HADV_C4 and W_VADV_SPLINES for horizontal and vertical advection respectively (see cppdefs_dev.h for default choices).
AGRIF nesting (1-WAY for now) is possible with NBQ and for NBQ integration in child grid only: #define NBQ_CHILD_ONLY

- Updated Ocean-Wave-Atmosphere coupling using the generic coupler OASIS3-MCT  
#define OA_COUPLING and/or #define OW_COUPLING
See specific documentation for coupling in the Documentation section of our website ; https://www.croco-ocean.org/documentation/  ; dedicated Coupling_tools are also available in croco_tools

- Updated Wave-Current interactions (McWilliams et al., JFM 2004)  
#define MRL_WCI 
Dedicated configuration : SHOREFACE

- Built-in wave propagation model (WKB) for nearshore applications 
#define WKB_WAVE

- Updates in sediments modules based on Warner et al. (2008): bedload transport, hydro-morphodynamics coupling and morphological acceleration factor 
#define SEDIMENT and associated cppkeys
Dedicated configuration :  FLUME, RIP

- Additional IO module (XIOS) providing more flexibility and better performances for HPC 
# define XIOS

- Choice of monotonic horizontal and vertical advection schemes for all variables (WENO5 &  TVD)
WENO5:  #define UV_HADV_ADV ;  UV_VADV_TVD ; W_HADV_TVD ;  W_VADV_TVD ; TS_VADV_WENO5 ; TS_HADV_WENO5
TVD:  # define UV_HADV_TVD ; UV_VADV_TVD ; W_HADV_TVD ; W_VADV_TVD

- Semi-implicit vertical advection for avoiding CFL limitation associated with “hot spots” (Shchepetkin, OM 2015) 
# define  VADV_ADAPT_IMP 
GLS turbulent closure sub-model (in addition to KPP)  
# define GLS_MIXING_2017 

- Wave‐induced (non breaking) vertical mixing in KPP modified according to Wang et al. (JGR 2010) based on Qiao et al (JPO 2004)  
# define WAVE_NONBRK_VMIX_QIAO in mrl_wci.F

- Basic 3D Smagorinsky model for LES-type applications   
# define UV_VIS_SMAGO_3D (with NBQ)

- Built-in diffusion in barotropic time stepping as an alternative to fast mode filtering (used in NBQ applications)  
# define M2FILTER_NONE

- Updated PISCES version (Aumont et Bopp, GBC 2006)  
# define PISCES 

- New vertical coordinate (not dependent on minimum depth), by default,  suited for decreasing pressure gradient errors in the thermocline above steep topography (Shchepetkin and McWilliams, 2009)
 # define NEW_S_COORD

- Dedicated log output file croco.log 
# define LOGFILE

ROMS_AGRIF  is not maintained anymore and we strongly encourage ROMS_AGRIF users to switch  to CROCO. 
