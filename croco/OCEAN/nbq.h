! $Id:$
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
#ifdef NBQ
!**********************************************************************
      logical WEST_INTER_NBQ         
      common /nbq_log1/ WEST_INTER_NBQ
      logical EAST_INTER_NBQ         
      common /nbq_log2/ EAST_INTER_NBQ                                         
      logical SOUTH_INTER_NBQ       
      common /nbq_log3/ SOUTH_INTER_NBQ
      logical NORTH_INTER_NBQ       
      common /nbq_log4/ NORTH_INTER_NBQ

!********************************************************************** 
      integer iteration_nbq_max  
      common /nbq_var1/ iteration_nbq_max             
      integer iteration_nbq      
      common /nbq_var2/ iteration_nbq      
      integer ifl_nbq     
      common /nbq_var3/ ifl_nbq     
      integer slip_nbq   
      common /nbq_var4/ slip_nbq    

!**********************************************************************
      real soundspeed_nbq(GLOBAL_2D_ARRAY) 
      common /nbq_param1/ soundspeed_nbq                              
      real soundspeed2_nbq(GLOBAL_2D_ARRAY) 
      common /nbq_param2/ soundspeed2_nbq

      double precision time_nbq               
      common /nbq_param3/ time_nbq
      double precision csvisc1_nbq                
      common /nbq_param4/ csvisc1_nbq 
      double precision csvisc2_nbq            
      common /nbq_param5/ csvisc2_nbq 
      double precision cw_int_nbq       
      common /nbq_param6/ cw_int_nbq
      double precision ifl_imp_nbq
      common /nbq_param7/ ifl_imp_nbq

!**********************************************************************
      integer ndtnbq
      common /time_nbq1/ ndtnbq
      real dtnbq
      common /time_nbq2/ dtnbq 
      real csound_nbq
      common /nbq_csound/ csound_nbq
      real visc2_nbq
      common /nbq_visc2/ visc2_nbq

      real dtgrid_nbq
      common /nbq_dtgrid/ dtgrid_nbq

!**********************************************************************
      real qdmu_nbq(GLOBAL_2D_ARRAY,N)
      common /nbq_qdmu_nbq/ qdmu_nbq
      real qdmv_nbq(GLOBAL_2D_ARRAY,N)
      common /nbq_qdmv_nbq/ qdmv_nbq
      real qdmw_nbq(GLOBAL_2D_ARRAY,0:N)
      common /nbq_qdmw_nbq/ qdmw_nbq

!**********************************************************************       
# ifdef MPI       
      real thetadiv_nbq(GLOBAL_2D_ARRAY,0:N)
# else
      real thetadiv_nbq(GLOBAL_2D_ARRAY_EXT_NBQ,0:N)
# endif      
      common /nbq_thetadiv_nbq/ thetadiv_nbq
      real thetadiv2_nbq(GLOBAL_2D_ARRAY,0:N)
      common /nbq_thetadiv2_nbq/ thetadiv2_nbq
      real thetadiv3_nbq(GLOBAL_2D_ARRAY,0:N)
      common /nbq_thetadiv3_nbq/ thetadiv3_nbq

!**********************************************************************
      real ru_int_nbq(GLOBAL_2D_ARRAY,N)
      common /nbq_ruint/ ru_int_nbq
      real rv_int_nbq(GLOBAL_2D_ARRAY,N)
      common /nbq_rvint/ rv_int_nbq
      real rw_int_nbq(GLOBAL_2D_ARRAY,0:N)
      common /nbq_rwint/ rw_int_nbq   

      real rho_nbq(GLOBAL_2D_ARRAY,N)
      common/nbq_rho_nbq/rho_nbq

      real ru_nbq(GLOBAL_2D_ARRAY,N)
      common /nbq_ru/ ru_nbq
      real rv_nbq(GLOBAL_2D_ARRAY,N)
      common /nbq_rv/ rv_nbq
      real rw_nbq(GLOBAL_2D_ARRAY,0:N)
      common /nbq_rw/ rw_nbq

      real ru_nbq_avg2(GLOBAL_2D_ARRAY,N)
      common /avg2_runbq/ ru_nbq_avg2
      real rv_nbq_avg2(GLOBAL_2D_ARRAY,N)
      common /avg2_rvnbq/ rv_nbq_avg2
      real rw_nbq_avg2(GLOBAL_2D_ARRAY,0:N)
      common /avg2_rwnbq/ rw_nbq_avg2

!**********************************************************************
      real DU_nbq(GLOBAL_2D_ARRAY)
      common /nbq_DU_nbq/ DU_nbq
      real DV_nbq(GLOBAL_2D_ARRAY)
      common /nbq_DV_nbq/ DV_nbq

# ifdef NBQ_MASS
      real rho_nbq_avg1(GLOBAL_2D_ARRAY,0:N)
      common /avg1_rhonbq/ rho_nbq_avg1
      real rhobar_nbq(GLOBAL_2D_ARRAY,4)
      common /nbq_rhobar/ rhobar_nbq
      real rhobar_nbq_avg1(GLOBAL_2D_ARRAY)
      common /nbq_rhobar_AVG1/ rhobar_nbq_avg1
#  ifndef NBQ_ZETAW
      real rhobar_nbq_int(GLOBAL_2D_ARRAY)
      common /nbq_rhobar_int/ rhobar_nbq_int
#  endif
# endif

# ifdef NBQ_ZETAW 
      real rho_bak(GLOBAL_2D_ARRAY,N)
      common/nbq_rho_bak/rho_bak
# endif

      real ru_int_nbq_2d (GLOBAL_2D_ARRAY)  
      common /nbq_ruint_2d/ru_int_nbq_2d
      real rv_int_nbq_2d (GLOBAL_2D_ARRAY)  
      common /nbq_rvint_2d/rv_int_nbq_2d

      real rubar_nbq(GLOBAL_2D_ARRAY)
      common /nbq_rubar/ rubar_nbq
      real rvbar_nbq(GLOBAL_2D_ARRAY)
      common /nbq_rvbar/ rvbar_nbq

!**********************************************************************
      real e(GLOBAL_2D_ARRAY)
      common /nbq_e/ e
      real eomn(GLOBAL_2D_ARRAY)
      common /nbq_eomn/ eomn
      real cosa(GLOBAL_2D_ARRAY)
      common /nbq_cosa/ cosa
      real sina(GLOBAL_2D_ARRAY)
      common /nbq_sina/ sina

!**********************************************************************
      real Hzw_half_nbq(GLOBAL_2D_ARRAY,0:N)
      common /grid_Hzw_half_nbq/ Hzw_half_nbq
# ifndef NBQ_ZETAW
      real zr_half_nbq(GLOBAL_2D_ARRAY,N)
      common /grid_zr_half_nbq/ zr_half_nbq
      real zw_half_nbq(GLOBAL_2D_ARRAY,0:N)
      common /grid_zw_half_nbq/ zw_half_nbq     
# endif    

# if defined NBQ_MASS && !defined NBQ_ZETAW 
       real Hzr_half_nbq(GLOBAL_2D_ARRAY,N)
       common /grid_Hzr_half_nbq/ Hzr_half_nbq
# endif

      real zr_nbq(GLOBAL_2D_ARRAY,N,4)
      common /nbq_zr/ zr_nbq
      real zw_nbq(GLOBAL_2D_ARRAY,0:N,4)
      common /nbq_zw/ zw_nbq

      real Hz_bak2(GLOBAL_2D_ARRAY,1:N)
      common /nbq_H_bak2/ Hz_bak2

# ifdef NBQ_HZCORRECT
       real Hz_correct(GLOBAL_2D_ARRAY,N)
       common /grid_Hz_correct/ Hz_correct
#  ifdef NBQ_HZCORR_DEBUG
      real  Hz_corr(GLOBAL_2D_ARRAY,N) 
      common/corr_Hz/Hz_corr
#  endif
# endif
  
!**********************************************************************
# ifdef NBQ_PERF
      real dthetadiv_nbqdz(GLOBAL_2D_ARRAY,0:N,2)
      common /nbq_nods3/ dthetadiv_nbqdz
      real dZdxq_w(GLOBAL_2D_ARRAY,0:N+1)
      common /nbq_nods5/ dZdxq_w
      real dZdyq_w(GLOBAL_2D_ARRAY,0:N+1)
      common /nbq_nods7/ dZdyq_w
# else
#  ifdef NBQ_ZETAW
      real dthetadiv_nbqdz(GLOBAL_2D_ARRAY)
      common /nbq_nods3/ dthetadiv_nbqdz
      real dZdxq_w(GLOBAL_2D_ARRAY,0:N+1)
      common /nbq_nods5/ dZdxq_w
      real dZdyq_w(GLOBAL_2D_ARRAY,0:N+1)
      common /nbq_nods7/ dZdyq_w
#  endif
# endif /* NBQ_PERF */

!**********************************************************************
# ifdef NBQ_ZETAW
      real wsurf_nbq(GLOBAL_2D_ARRAY,4)
      common /nbq_wsurf/wsurf_nbq
      real usurf_nbq(GLOBAL_2D_ARRAY)
      common /nbq_usurf/usurf_nbq
      real vsurf_nbq(GLOBAL_2D_ARRAY)
      common /nbq_vsurf/vsurf_nbq

      integer :: kbak2,kold2
      common /gridext1/kbak2,kold2

      real cff4
      common /com_cff4/ cff4
      real cff5
      common /com_cff5/ cff5
      real cff6
      common /com_cff6/ cff6
      real cff7
      common /com_cff7/ cff7
      real cff8
      common /com_cff8/ cff8
      real cff9
      common /com_cff9/ cff9
      real cff10
      common /com_cff10/ cff10

      integer IstrU2,JstrV2
      integer IstrR2,IendR2
      integer JstrR2,JendR2
# endif /* NBQ_ZETAW */

!**********************************************************************
# if defined OBC_NBQ && defined OBC_NBQORLANSKI
#  ifdef OBC_WEST
      real qdmu_nbq_west(GLOBAL_1D_ARRAYETA,N,2)
      common /bry_unbq_west/ qdmu_nbq_west
      real qdmv_nbq_west(GLOBAL_1D_ARRAYETA,N,2)
      common /bry_vnbq_west/ qdmv_nbq_west
      real qdmw_nbq_west(GLOBAL_1D_ARRAYETA,0:N,2)
      common /bry_wnbq_west/ qdmw_nbq_west
      real  rho_nbq_west(GLOBAL_1D_ARRAYETA,N,2)
      common /bry_rnbq_west/ rho_nbq_west
#  endif
#  ifdef OBC_EAST
      real qdmu_nbq_east(GLOBAL_1D_ARRAYETA,N,2)
      common /bry_unbq_east/ qdmu_nbq_east
      real qdmv_nbq_east(GLOBAL_1D_ARRAYETA,N,2)
      common /bry_vnbq_east/ qdmv_nbq_east
      real qdmw_nbq_east(GLOBAL_1D_ARRAYETA,0:N,2)
      common /bry_wnbq_east/ qdmw_nbq_east
      real  rho_nbq_east(GLOBAL_1D_ARRAYETA,N,2)
      common /bry_rnbq_east/ rho_nbq_east
#  endif
#  ifdef OBC_SOUTH
      real qdmu_nbq_south(GLOBAL_1D_ARRAYXI,N,2)
      common /bry_unbq_south/ qdmu_nbq_south
      real qdmv_nbq_south(GLOBAL_1D_ARRAYXI,N,2)
      common /bry_vnbq_south/ qdmv_nbq_south
      real qdmw_nbq_south(GLOBAL_1D_ARRAYXI,0:N,2)
      common /bry_wnbq_south/ qdmw_nbq_south
      real  rho_nbq_south(GLOBAL_1D_ARRAYXI,N,2)
      common /bry_rnbq_south/ rho_nbq_south
#  endif
#  ifdef OBC_NORTH
      real qdmu_nbq_north(GLOBAL_1D_ARRAYXI,N,2)
      common /bry_unbq_north/ qdmu_nbq_north
      real qdmv_nbq_north(GLOBAL_1D_ARRAYXI,N,2)
      common /bry_vnbq_north/ qdmv_nbq_north
      real qdmw_nbq_north(GLOBAL_1D_ARRAYXI,0:N,2)
      common /bry_wnbq_north/ qdmw_nbq_north
      real  rho_nbq_north(GLOBAL_1D_ARRAYXI,N,2)
      common /bry_rnbq_north/ rho_nbq_north
#  endif
# endif     

!**********************************************************************
# ifdef NBQ_NUDGING
      real NBQnudgcof(GLOBAL_2D_ARRAY)
      common /nbq_nudg/ NBQnudgcof
# endif

!**********************************************************************
# ifdef ACOUSTIC
      real  period_exp  
      common/ACOUS1/period_exp
      real  for_a_exp   
      common/ACOUS2/for_a_exp
      real  dg_exp     
      common/ACOUS3/dg_exp 
      real  hmax_exp    
      common/ACOUS4/hmax_exp
      real  amp_exp
      common/ACOUS4/amp_exp
# endif

#endif /* NBQ */

  
