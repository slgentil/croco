! $Id: nbq_bry_store 1619 2015-01-07 13:53:03Z marchesiello $
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
!
#if defined OBC_NBQ && defined OBC_NBQORLANSKI
!
!======================================================================
!                   subroutine nbq_bry_store
!======================================================================
!
      subroutine nbq_bry_store_tile (Istr,Iend,Jstr,Jend) 
!
!----------------------------------------------------------------------
!  Store boundary values of nbq variables 
!  for use in radiative conditions
!----------------------------------------------------------------------
!
      implicit none          
      integer Istr,Iend,Jstr,Jend, i,j,k, itrc
# include "param.h"
# include "scalars.h"
# include "nbq.h"
!
# include "compute_auxiliary_bounds.h"
!
# ifdef OBC_WEST
      if (WESTERN_EDGE) then
        do k=1,N
          do j=Jstr-1,Jend+1
            qdmu_nbq_west(j,k,1)=qdmu_nbq(Istr  ,j,k)
            qdmu_nbq_west(j,k,2)=qdmu_nbq(Istr+1,j,k)
            qdmv_nbq_west(j,k,1)=qdmv_nbq(Istr-1,j,k)
            qdmv_nbq_west(j,k,2)=qdmv_nbq(Istr  ,j,k)
            qdmw_nbq_west(j,k,1)=qdmw_nbq(Istr-1,j,k)
            qdmw_nbq_west(j,k,2)=qdmw_nbq(Istr  ,j,k)
             rho_nbq_west(j,k,1)= rho_nbq(Istr-1,j,k)
             rho_nbq_west(j,k,2)= rho_nbq(Istr  ,j,k)
          enddo
        enddo
        do j=Jstr-1,Jend+1
          qdmw_nbq_west(j,0,1)=qdmw_nbq(Istr-1,j,0)
          qdmw_nbq_west(j,0,2)=qdmw_nbq(Istr  ,j,0)
        enddo
      endif
# endif /* OBC_WEST */
!
# ifdef OBC_EAST
      if (EASTERN_EDGE) then
        do k=1,N
          do j=Jstr-1,Jend+1
            qdmu_nbq_east(j,k,1)=qdmu_nbq(Iend+1,j,k)
            qdmu_nbq_east(j,k,2)=qdmu_nbq(Iend  ,j,k)
            qdmv_nbq_east(j,k,1)=qdmv_nbq(Iend+1,j,k)
            qdmv_nbq_east(j,k,2)=qdmv_nbq(Iend  ,j,k)
            qdmw_nbq_east(j,k,1)=qdmw_nbq(Iend+1,j,k)
            qdmw_nbq_east(j,k,2)=qdmw_nbq(Iend  ,j,k)
             rho_nbq_east(j,k,1)= rho_nbq(Iend+1,j,k)
             rho_nbq_east(j,k,2)= rho_nbq(Iend  ,j,k)
          enddo
        enddo
        do j=Jstr-1,Jend+1
          qdmw_nbq_east(j,0,1)=qdmw_nbq(Iend+1,j,0)
          qdmw_nbq_east(j,0,2)=qdmw_nbq(Iend  ,j,0)
        enddo
      endif
# endif /* OBC_EAST */
!
# ifdef OBC_SOUTH
      if (SOUTHERN_EDGE) then
        do k=1,N
          do i=Istr-1,Iend+1
            qdmu_nbq_south(i,k,1)=qdmu_nbq(i,Jstr-1,k)
            qdmu_nbq_south(i,k,2)=qdmu_nbq(i,Jstr  ,k)
            qdmv_nbq_south(i,k,1)=qdmv_nbq(i,Jstr  ,k)
            qdmv_nbq_south(i,k,2)=qdmv_nbq(i,Jstr+1,k)
            qdmw_nbq_south(i,k,1)=qdmw_nbq(i,Jstr-1,k)
            qdmw_nbq_south(i,k,2)=qdmw_nbq(i,Jstr  ,k)
             rho_nbq_south(i,k,1)= rho_nbq(i,Jstr-1,k)
             rho_nbq_south(i,k,2)= rho_nbq(i,Jstr  ,k)
          enddo
        enddo
        do i=Istr-1,Iend+1
          qdmw_nbq_south(i,0,1)=qdmw_nbq(i,Jstr-1,0)
          qdmw_nbq_south(i,0,2)=qdmw_nbq(i,Jstr  ,0)
        enddo
      endif
# endif /* OBC_SOUTH */
!
# ifdef OBC_NORTH
      if (NORTHERN_EDGE) then
        do k=1,N
          do i=Istr-1,Iend+1
            qdmu_nbq_north(i,k,1)=qdmu_nbq(i,Jend+1,k)
            qdmu_nbq_north(i,k,2)=qdmu_nbq(i,Jend  ,k)
            qdmv_nbq_north(i,k,1)=qdmv_nbq(i,Jend+1,k)
            qdmv_nbq_north(i,k,2)=qdmv_nbq(i,Jend  ,k)
            qdmw_nbq_north(i,k,1)=qdmw_nbq(i,Jend+1,k)
            qdmw_nbq_north(i,k,2)=qdmw_nbq(i,Jend  ,k)
             rho_nbq_north(i,k,1)= rho_nbq(i,Jend+1,k)
             rho_nbq_north(i,k,2)= rho_nbq(i,Jend  ,k)
          enddo
        enddo
        do i=Istr-1,Iend+1
          qdmw_nbq_north(i,0,1)=qdmw_nbq(i,Jend+1,0)
          qdmw_nbq_north(i,0,2)=qdmw_nbq(i,Jend  ,0)
        enddo
      endif
# endif /* OBC_NORTH */
!
      return
      end
#endif /* OBC_NBQ && OBC_NBQORLANSKI */

