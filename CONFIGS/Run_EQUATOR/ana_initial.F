! $Id: ana_initial.F 1620 2015-01-08 10:47:13Z marchesiello $
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
#include "cppdefs.h"

      subroutine ana_initial (tile)
      implicit none
#include "param.h"
       integer tile, trd
C$    integer omp_get_thread_num
#include "compute_tile_bounds.h"
      call ana_initial_tile (Istr,Iend,Jstr,Jend)
      return
      end
!
      subroutine ana_initial_tile (Istr,Iend,Jstr,Jend)
!
!----------------------------------------------------------------------
!  Set initial conditions for momentum and tracer variables using
!  analytical expressions.
!----------------------------------------------------------------------
!
#ifdef AGRIF
      use Agrif_Util
#endif
      implicit none
#include "param.h"
#include "grid.h"
#include "ocean2d.h"
#include "ocean3d.h"
#include "scalars.h"
#ifdef SEDIMENT
# include "sediment.h"
#endif
      integer Istr,Iend,Jstr,Jend, i,j,k, itrc
#if defined NBQ && ! defined MPI
      integer :: mynode
#endif
#ifdef SOLITON
      real fac, x,y, x0,y0, cff1,cff2,cff3
#endif
#ifdef KH_INST 
      real rho1_exp,rho2_exp,h1_exp, Lt
      real d(0:N)
# endif
#ifdef GRAV_ADJ_SOLITON
      real rho1_exp,rho2_exp,h1_exp,amp_exp, Lt
      real d(0:N)
#endif
#ifdef BASIN
      real cff1,cff2
#endif
#ifdef INTERNAL
      real cff
#endif
#if defined IGW || defined S2DV || defined MILES 
      real    cff,sig,rhom,depth
#endif 
#ifdef BIOLOGY
      real temp, SiO4
#endif
#ifdef THACKER
      real ETA0, omega, depth, Lt
#endif
#ifdef TANK
      real ETA0, depth, Lt,cff
#endif
#ifdef RIP
      integer iper
      real xper,yper,xx,yy,xs,db,eps
#endif
#ifdef KH_INST 
      real zu0,zr0,zint,depth,dz,gam,del,hint,du,drho
      real hu0,hup0,htemp0,hsalt0,eps3D
#endif
#ifdef GRAV_ADJ_SOLITON
      real zint,depth,dz,gam,del
#endif
!
#include "compute_auxiliary_bounds.h"
!
#if defined EW_PERIODIC && !defined MPI
# define IR_RANGE Istr,Iend
# define IU_RANGE Istr,Iend
#else
# define IR_RANGE IstrR,IendR
# define IU_RANGE Istr,IendR
#endif
#if defined NS_PERIODIC && !defined MPI
# define JR_RANGE Jstr,Jend
# define JV_RANGE Jstr,Jend
#else
# define JR_RANGE JstrR,JendR
# define JV_RANGE Jstr,JendR
#endif

#ifdef ANA_INITIAL
!
!----------------------------------------------------------------------
!  Initial conditions for free surface and 2D momentum components.
!----------------------------------------------------------------------
!
      if (nrrec.eq.0) then
!
# if defined BASIN_EQ

        do j=JR_RANGE
          do i=IR_RANGE
            zeta(i,j,1)=0.
            zeta(i,j,2)=0.
            ubar(i,j,1)=0.
            vbar(i,j,1)=0.
            ubar(i,j,2)=0.
            vbar(i,j,2)=0.
          enddo
        enddo
# else
        do j=JR_RANGE
          do i=IR_RANGE
            zeta(i,j,1)=0.
            zeta(i,j,2)=0.
            ubar(i,j,1)=0.
            vbar(i,j,1)=0.
            ubar(i,j,2)=0.
            vbar(i,j,2)=0.
          enddo
        enddo
# endif

      else  !  nrrec.ne.0.


      endif  !  nrrec.eq.0.

# ifdef SOLVE3D
!
!----------------------------------------------------------------------
!  Initial conditions for momentum components [m/s].
!----------------------------------------------------------------------
!
      if (nrrec.eq.0) then
!
#  if defined BASIN_EQ

        do k=1,N
          do j=JR_RANGE
            do i=IR_RANGE
              u(i,j,k,1)=0.
              v(i,j,k,1)=0.
              u(i,j,k,2)=0.
              v(i,j,k,2)=0.
            enddo
          enddo
        enddo
#  else
        do k=1,N
          do j=JR_RANGE
            do i=IR_RANGE
              u(i,j,k,1)=0.
              v(i,j,k,1)=0.
              u(i,j,k,2)=0.
              v(i,j,k,2)=0.
            enddo
          enddo
        enddo
#  endif

      else   ! nrrec.ne.0.
!
      endif   !nrrec.eq.0.
!
!----------------------------------------------------------------------
!  Set initial conditions for potential temperature [degC] and
!  salinity [PSU] and passive tracer
!----------------------------------------------------------------------
!
#  if defined BASIN_EQ

      call init_rhobar(Istr,Iend,Jstr,Jend)

#  endif

      if (nrrec.eq.0) then
!
#  if defined BASIN_EQ

c      print *,'apres init_rhobar  xnconst=',xnconst
      
      do k=1,N
        do j=JR_RANGE
          do i=IR_RANGE
            t(i,j,k,1,itemp)=tbar(i,j,k)
            t(i,j,k,2,itemp)=tbar(i,j,k)
          enddo
        enddo
      enddo

#  else
        do k=1,N
          do j=JR_RANGE
            do i=IR_RANGE
              t(i,j,k,1,itemp)=20.
              t(i,j,k,2,itemp)=20.
#   ifdef SALINITY
              t(i,j,k,1,isalt)=35.
              t(i,j,k,2,isalt)=35.
#   endif
            enddo
          enddo
        enddo
#  endif
!
      endif   !  nrrec.eq.0.
!
# endif /*SOLVE3D */
#endif /* ANA_INITIAL */

#if defined BASIN_EQ

#   ifdef SPONGE_LAYER
        call init_sponge(Istr,Iend,Jstr,Jend)
#   endif

#   ifdef VTFORC_INT
      if (iic.eq.ntstart) then
         call vtforc(Istr,Iend,Jstr,Jend)
      endif
#   endif  /* VTFORC_INT */ 

# endif  
!!
!--------------------------------------------------------------------
!  Exchange boundary information
!--------------------------------------------------------------------
!

#if defined EW_PERIODIC || defined NS_PERIODIC || defined  MPI
# ifdef ANA_INITIAL
      call exchange_r2d_tile (Istr,Iend,Jstr,Jend,
     &                        zeta(START_2D_ARRAY,1))
      call exchange_u2d_tile (Istr,Iend,Jstr,Jend,
     &                        ubar(START_2D_ARRAY,1))
      call exchange_v2d_tile (Istr,Iend,Jstr,Jend,
     &                        vbar(START_2D_ARRAY,1))
#  ifdef SOLVE3D
#   ifdef THREE_GHOST_POINTS_UV
      call exchange_u3d_3pts_tile (Istr,Iend,Jstr,Jend,
     &                             u(START_2D_ARRAY,1,1))
      call exchange_v3d_3pts_tile (Istr,Iend,Jstr,Jend,
     &                             v(START_2D_ARRAY,1,1))
#   else
      call exchange_u3d_tile (Istr,Iend,Jstr,Jend,
     &                        u(START_2D_ARRAY,1,1))
      call exchange_v3d_tile (Istr,Iend,Jstr,Jend,
     &                        v(START_2D_ARRAY,1,1))
#   endif
#  endif /* SOLVE3D */
# endif /* ANA_INITIAL */
# ifdef SOLVE3D
      do itrc=1,NT
        if (.not.got_tini(itrc)) then
#  ifdef THREE_GHOST_POINTS_TS
          call exchange_r3d_3pts_tile (Istr,Iend,Jstr,Jend,
     &                                 t(START_2D_ARRAY,1,1,itrc))
#  else
          call exchange_r3d_tile (Istr,Iend,Jstr,Jend,
     &                            t(START_2D_ARRAY,1,1,itrc))
#  endif
        endif
      enddo
# endif
#endif
#undef IR_RANGE
#undef IU_RANGE
#undef JR_RANGE
#undef JV_RANGE

      return
      end


