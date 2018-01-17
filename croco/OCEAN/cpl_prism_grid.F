! $Id: cpl_prism_grid.F Joris Pianezze 9 aout 2016 $
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
      SUBROUTINE cpl_prism_grid(id_partition)
#include "cppdefs.h"
!     IMPLICIT NONE
#if defined OA_COUPLING || defined OW_COUPLING
#ifdef OA_MCT
      USE mod_prism             ! OASIS3 prism module
      INCLUDE 'mpif.h'
#endif

#include "param.h"
#include "grid.h"
#include "scalars.h"
#include "mpi_roms.h"
#include "set_global_definitions.h"
!
      integer, intent(in) :: id_partition
!
      integer error_flg
!
      real xbox, ybox
      real clonr(GLOBAL_2D_ARRAY, 4)
      real clonu(GLOBAL_2D_ARRAY, 4)
      real clonv(GLOBAL_2D_ARRAY, 4)
      real clatr(GLOBAL_2D_ARRAY, 4)
      real clatu(GLOBAL_2D_ARRAY, 4)
      real clatv(GLOBAL_2D_ARRAY, 4)
      real arear(GLOBAL_2D_ARRAY)
      real areau(GLOBAL_2D_ARRAY)
      real areav(GLOBAL_2D_ARRAY)
      integer rmask_inv(GLOBAL_2D_ARRAY)
      integer umask_inv(GLOBAL_2D_ARRAY)
      integer vmask_inv(GLOBAL_2D_ARRAY)
!
!
!     !-------------------------------------------------------------------
!     !             ***  ROUTINE cpl_prism_grid  ***
!     !
!     ! ** Purpose :   Define grids/areas/corners/masks information for ocean
!     !                exchange between ROMS and COUPLER. (OASIS3-MCT software)
!     !
!     ! ** Method  :   OASIS3 MPI communication
!     !--------------------------------------------------------------------

!------------------------------------------------------------------
!     Set the ROMS grid definition to PRISM system
!------------------------------------------------------------------
!
!    Need to use oasis3_mct_v3 !
!
      CALL prism_start_grids_writing(error_flg)

! Model grid longitudes and latitudes
!------------------------------------------------------------------
!     -- t-grid --
      CALL prism_write_grid('ocnt',LLm,MMm,
     &           lonr(1:Lm,1:Mm),latr(1:Lm,1:Mm),
     &           id_partition)
!     -- u-grid --
      CALL prism_write_grid('ocnu',LLm,MMm,
     &           lonu(1:Lm,1:Mm),latu(1:Lm,1:Mm),
     &           id_partition)
!     -- v-grid --
      CALL prism_write_grid('ocnv',LLm,MMm,
     &           lonv(1:Lm,1:Mm),latv(1:Lm,1:Mm),
     &           id_partition)

! Model grid cell corner longitudes and latitudes
!------------------------------------------------------------------
!     -- t-grid --
      xbox=lonr(1,2)-lonr(1,1)
      clonr(:,:,1)=lonr(:,:)+xbox/2.
      clonr(:,:,2)=lonr(:,:)-xbox/2.
      clonr(:,:,3)=lonr(:,:)-xbox/2.
      clonr(:,:,4)=lonr(:,:)+xbox/2.  
      ybox=latr(2,1)-latr(1,1)
      clatr(:,:,1)=latr(:,:)+ybox/2.
      clatr(:,:,2)=latr(:,:)+ybox/2.
      clatr(:,:,3)=latr(:,:)-ybox/2.
      clatr(:,:,4)=latr(:,:)-ybox/2.
      CALL prism_write_corner('ocnt', LLm, MMm, 4,
     &             clonr(1:Lm,1:Mm,:), clatr(1:Lm,1:Mm,:),
     &             id_partition)

!     -- u-grid --
      xbox=lonu(1,2)-lonu(1,1)
      clonu(:,:,1)=lonu(:,:)+xbox/2.
      clonu(:,:,2)=lonu(:,:)-xbox/2.
      clonu(:,:,3)=lonu(:,:)-xbox/2.
      clonu(:,:,4)=lonu(:,:)+xbox/2.  
      ybox=latu(2,1)-latu(1,1)
      clatu(:,:,1)=latu(:,:)+ybox/2.
      clatu(:,:,2)=latu(:,:)+ybox/2.
      clatu(:,:,3)=latu(:,:)-ybox/2.
      clatu(:,:,4)=latu(:,:)-ybox/2.
      CALL prism_write_corner('ocnu', LLm, MMm, 4,
     &             clonu(1:Lm,1:Mm,:), clatu(1:Lm,1:Mm,:),
     &             id_partition)

!     -- v-grid --
      xbox=lonv(1,2)-lonv(1,1)
      clonv(:,:,1)=lonv(:,:)+xbox/2.
      clonv(:,:,2)=lonv(:,:)-xbox/2.
      clonv(:,:,3)=lonv(:,:)-xbox/2.
      clonv(:,:,4)=lonv(:,:)+xbox/2.  
      ybox=latv(2,1)-latv(1,1)
      clatv(:,:,1)=latv(:,:)+ybox/2.
      clatv(:,:,2)=latv(:,:)+ybox/2.
      clatv(:,:,3)=latv(:,:)-ybox/2.
      clatv(:,:,4)=latv(:,:)-ybox/2.
      CALL prism_write_corner('ocnv', LLm, MMm, 4,
     &             clonv(1:Lm,1:Mm,:), clatv(1:Lm,1:Mm,:),
     &             id_partition)

! Model grid cell area
!------------------------------------------------------------------
!     -- t-grid --
      arear(:,:)=1./(pm(:,:)*pn(:,:))
      CALL prism_write_area('ocnt', LLm, MMm,
     &                         arear(1:Lm,1:Mm),
     &                         id_partition)

!     -- u-grid --
      areau(:,:)=1./(pm(:,:)*pn(:,:))
      CALL prism_write_area('ocnu', LLm, MMm,
     &                         areau(1:Lm,1:Mm),
     &                         id_partition)

!     -- v-grid --
      areav(:,:)=1./(pm(:,:)*pn(:,:))
      CALL prism_write_area('ocnv', LLm, MMm,
     &                         areav(1:Lm,1:Mm),
     &                         id_partition)

! Model grid mask : 0=sea / 1=land
!------------------------------------------------------------------
      rmask_inv=0
      umask_inv=0
      vmask_inv=0
      DO i=1, Lm
        DO j=1, Mm
          IF (rmask(i,j)==0) rmask_inv(i,j)=1
          IF (umask(i,j)==0) umask_inv(i,j)=1
          IF (vmask(i,j)==0) vmask_inv(i,j)=1
        ENDDO
      ENDDO

!     -- t-grid --
      CALL prism_write_mask('ocnt', LLm, MMm,
     &                     rmask_inv(1:Lm,1:Mm),
     &                     id_partition)

!     -- u-grid --
      CALL prism_write_mask('ocnu', LLm, MMm,
     &                     umask_inv(1:Lm,1:Mm),
     &                     id_partition)

!     -- v-grid --
      CALL prism_write_mask('ocnv', LLm, MMm,
     &                     vmask_inv(1:Lm,1:Mm),
     &                     id_partition)
      
! Terminate grid writing
!------------------------------------------------------------------
      CALL prism_terminate_grids_writing()

      CALL MPI_Barrier(ocean_grid_comm, error_flg)

!------------------------------------------------------------------
!     End of grid definition
!------------------------------------------------------------------

#else 
!---------------------------------------------------------------------------------------
!     Dummy subroutine                                                           NO Coupling
!---------------------------------------------------------------------------------------
#endif 

      END SUBROUTINE cpl_prism_grid
