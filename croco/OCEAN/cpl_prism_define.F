! $Id: cpl_prism_define.F 1274 2013-09-16 09:24:57Z gcambon $
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
      SUBROUTINE cpl_prism_define
#include "cppdefs.h"
#if defined OA_COUPLING || defined OW_COUPLING
      USE mod_prism             ! OASIS3 prism module

      IMPLICIT NONE

#include "param.h"
#include "grid.h"
#include "scalars.h"
#include "mpi_roms.h"
#include "ncscrum.h"
      integer i, j, ji, i0, j0, ja
      integer nproc, myproc, error_flg, ierr
      integer start_x, start_y
      real delta_lon, delta_lat
      integer ncid, varid, nf_fread
      INTEGER    :: level
      CHARACTER*1 clevel
      CHARACTER*5 catm
      character fname*64
#include "netcdf.inc"

!     !-------------------------------------------------------------------
!     !             ***  ROUTINE cpl_prism_define  ***
!     !
!     ! ** Purpose :   Define grid and field information for ocean
!     !                exchange between ROMS and COUPLER. (OASIS3-MCT software)
!     !
!     ! ** Method  :   OASIS3 MPI communication
!     !--------------------------------------------------------------------

!------------------------------------------------------------------
!     Set the ROMS grid definition to PRISM system
!------------------------------------------------------------------

      start_x =  iminmpi
      start_y =  jminmpi
      
!     Set the ROMS partition definition to PRISM system
!------------------------------------------------------------------

      oasis_ig_paral(1) = 2     ! Box partition
      oasis_ig_paral(2) = (LLm)*(start_y-1) + (start_x-1)
      oasis_ig_paral(3) = Lmmpi
      oasis_ig_paral(4) = Mmmpi
      oasis_ig_paral(5) = LLm


      CALL prism_def_partition_proto(oasis_part_id, 
     &     oasis_ig_paral, comp_ierror)
      IF ( comp_ierror /= PRISM_Ok ) THEN
        CALL prism_abort_proto ( comp_id, 'cpl_prism_grid',
     &       'Failure in prism_def_partition')
      END IF

!------------------------------------------------------------------
!     Set the ROMS grid definition to PRISM system
!------------------------------------------------------------------

      CALL cpl_prism_grid(oasis_part_id)

!------------------------------------------------------------------
!     Set the WRF exchange variables definition to PRISM system
!------------------------------------------------------------------

      oasis_var_nodims(1) = 2   ! Rank of the field array is 2
      oasis_var_nodims(2) = 1   ! Bundles always 1 for OASIS3

      oasis_var_type = PRISM_Real

      oasis_var_shape(1) = 1
      oasis_var_shape(2) = Lmmpi !Lm+2+padd_X + 1  !Lm
      oasis_var_shape(3) = 1
      oasis_var_shape(4) = Mmmpi !Mm+2+padd_E + 1  !Mm

! default definiton
      ssnd_clname(:) = 'not defined' 
      srcv_clname(:) = 'not defined' 
#ifdef AGRIF
      level=Agrif_Fixed()
#else
      level=0
#endif
      write(clevel,'(i1.1)') level

!------------------------------------------------------------------
!     Sent
!------------------------------------------------------------------
!     SST
      ssnd_clname(1) = 'SRMSSTV'//clevel
!     UOCE
      ssnd_clname(2) = 'SRMUOCE'//clevel
!     VOCE
      ssnd_clname(3) = 'SRMVOCE'//clevel
!------------------------------------------------------------------
!     Received
!------------------------------------------------------------------
!     Evaporation-Precipitation [kg/m^2/s]
      srcv_clname(1) = 'RRMEVPR'//clevel
!     net surface SW radiation [W/m2]
      srcv_clname(2) = 'RRMSRFL'//clevel
!     total non-solar heat flux [W/m2]
      srcv_clname(3) = 'RRMSTFL'//clevel
!     zonal wind tress [N/m2]
      srcv_clname(4) = 'RRMTAUX'//clevel
!     meridional wind tress [N/m2]
      srcv_clname(5) = 'RRMTAUY'//clevel
!     wind tress module [N/m2]
      srcv_clname(6) = 'RRMTAUM'//clevel
      
      DO ja = 0, nmaxatm

         IF ( ja == 0 ) THEN   ;   catm = ""
         ELSE                  ;   WRITE(catm,'(a,i1.1)') "_atm", ja
         ENDIF
!     
!     ... Announce sent variables.
      DO ji = 1, nmaxfld
          CALL prism_def_var_proto ( ssnd_nid(ja,ji),
     &       ssnd_clname(ji)//trim(catm),
     &       oasis_part_id, oasis_var_nodims,
     &       PRISM_Out, oasis_var_shape,
     &       oasis_var_type, comp_ierror)
          IF ( comp_ierror /= PRISM_Ok ) THEN
             CALL prism_abort_proto( ssnd_nid(ja,ji),'cpl_prism_define',
     &        'Failure in prism_def_var: '//ssnd_clname(ji)//trim(catm))
          ELSE
             IF ( ssnd_nid(ja,ji) /= -1 ) THEN
                MPI_master_only write(*,*) 
     & 'CPL-ROMS: sent ',ssnd_clname(ji)//trim(catm),ssnd_nid(ja,ji)
             ENDIF
          ENDIF
      END DO

!     ... Announce received variables.
      DO ji = 1, nmaxfld
         CALL prism_def_var_proto ( srcv_nid(ja,ji),
     &       srcv_clname(ji)//trim(catm),
     &       oasis_part_id, oasis_var_nodims, 
     &       PRISM_In, oasis_var_shape, 
     &       oasis_var_type, comp_ierror)
          IF ( comp_ierror /= PRISM_Ok ) THEN
             CALL prism_abort_proto( srcv_nid(ja,ji),'cpl_prism_define',
     &        'Failure in prism_def_var: '//srcv_clname(ji)//trim(catm))
          ELSE
             IF ( srcv_nid(ja,ji) /= -1 ) THEN
                MPI_master_only write(*,*) 
     & 'CPL-ROMS: received ',srcv_clname(ji)//trim(catm),srcv_nid(ja,ji)
             ENDIF
          ENDIF
      END DO

      END DO

!------------------------------------------------------------------
!     End of definition phase
!------------------------------------------------------------------
      CALL MPI_Barrier( MPI_COMM_WORLD, comp_ierror ) 
#ifdef AGRIF
      if (.not.Agrif_Root()) then
         if (Agrif_fixed() == Agrif_Nb_Fine_Grids()) then
#endif
            CALL prism_enddef_proto( comp_ierror )
            IF ( comp_ierror /= PRISM_Ok ) THEN
               CALL prism_abort_proto ( comp_id, 'cpl_prism_enddef', 
     &              'Failure in prism_enddef')
            ELSE
               MPI_master_only write(*,*) 
     &              'CPL-ROMS: ENDDEF OK'
            ENDIF
#ifdef AGRIF
         endif
      endif
#endif

      cplmsk(:,:,0) = 1.            ! default use, no coupling mask
      cplmsk(:,:,1:nmaxatm) = 0.    ! default definition when using coupling mask
      fname = 'coupling_masks'//clevel//'.nc'
# if defined MPI && defined PARALLEL_FILES
      ierr=0
      call insert_node (fname, len_trim(fname), mynode, NNODES, ierr)
# endif
      ierr=nf_open(trim(fname), nf_nowrite, ncid)
      if (ierr == NF_NOERR) then
          MPI_master_only write(*,*) 
     &    'CPL-ROMS: open coupling_masks'//clevel//'.nc ok'
          DO ja = 1, nmaxatm
             WRITE(catm,'(a,i1.1)') "_atm", ja
             ierr=nf_inq_varid (ncid, 'cplmsk'//catm, varid)
             if (ierr == NF_NOERR) then
      ierr=nf_fread (cplmsk(START_2D_ARRAY,ja), ncid, varid, 0, r2dvar)
                MPI_master_only write(*,*) 
     &               'CPL-ROMS: read cplmsk'//catm//' ok',
     &               minval(cplmsk(1:Lmmpi,1:Mmmpi,ja)),
     &               maxval(cplmsk(1:Lmmpi,1:Mmmpi,ja))
             endif
          enddo
          ierr=nf_close(ncid)
       else
          MPI_master_only write(*,*) 
     &         'CPL-ROMS: coupling_masks'//clevel//'.nc not found'
       endif


#else 
!---------------------------------------------------------------------------------------
!     Dummy subroutine                                                           NO Coupling
!---------------------------------------------------------------------------------------
#endif 

      END SUBROUTINE cpl_prism_define


