! $Id: cpl_prism_put.F 1274 2013-09-16 09:24:57Z gcambon $
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
      SUBROUTINE cpl_prism_put (sndtime)
!     !---------------------------------------------------------------------
!     !              ***  ROUTINE cpl_prism_put  ***
!     !
!     ! ** Purpose : - 
!     !           
!     !
!     !----------------------------------------------------------------------
#include "cppdefs.h"
#if defined OA_COUPLING || defined OW_COUPLING
      USE mod_prism             ! OASIS3 prism module
# ifdef AGRIF
      use Agrif_Util
# endif
      IMPLICIT NONE

# include "param.h"
# include "scalars.h"
# include "ocean3d.h"
# include "ocean2d.h"
# include "mpi_cpl.h"
      
      INTEGER    :: sndtime, i, level, ja
      DOUBLE PRECISION, DIMENSION(Lmmpi,Mmmpi)  :: snd_field 
      CHARACTER*1 clevel
      
# ifdef AGRIF
      level=Agrif_Fixed()
# else
      level=0
# endif
      write(clevel,'(i1.1)') level
!
!!$      if (level==0) then
!!$         MPI_master_only 
!!$     & write(*,*) 'CPL-CROCO : put for parent at time', sndtime
!!$      else		
!!$         MPI_master_only 
!!$     & write(*,*) 'CPL-CROCO : put for child:', level,
!!$     &  'at time', sndtime
!!$      endif
!      
      DO i = 1, nmaxfld
      
         if (ssnd_clname(i)=='SRMSSTV'//clevel) then
            snd_field = t(1:Lmmpi,1:Mmmpi,N,nnew,itemp) + 273.15    

      DO ja = 0, nmaxatm
         if ( ssnd_nid(ja,i) /= -1 ) then
            call prism_put_proto(ssnd_nid(ja,i), sndtime, 
     &        snd_field, comp_ierror)
            if ( (comp_ierror == OASIS_Sent   ) .or.
     &           (comp_ierror == OASIS_ToRest ) .or.
     &           (comp_ierror == OASIS_SentOut) .or.
     &           (comp_ierror == OASIS_ToRestOut) ) then
           MPI_master_only write(*,*) 'CPL-CROCO : ',ja,
     & ssnd_clname(i),' sent at time ',sndtime
            endif
         endif  
      enddo
         
         endif

         if (ssnd_clname(i)=='SRMSSHV'//clevel) then
            snd_field = zeta(1:Lmmpi,1:Mmmpi,nnew)
            	
      DO ja = 0, nmaxatm
         if ( ssnd_nid(ja,i) /= -1 ) then             
            call prism_put_proto(ssnd_nid(ja,i), sndtime,
     &        snd_field, comp_ierror)
                 if ( (comp_ierror == OASIS_Sent   ) .or.
     &           (comp_ierror == OASIS_ToRest ) .or.
     &           (comp_ierror == OASIS_SentOut) .or.
     &           (comp_ierror == OASIS_ToRestOut) ) then
           MPI_master_only write(*,*) 'CPL-CROCO : ',ja,
     & ssnd_clname(i),' sent at time ',sndtime
            endif
         endif
      enddo

         endif
      
         if (ssnd_clname(i)=='SRMUOCE'//clevel) then
            snd_field = 0.5 * (u(1:Lmmpi  ,1:Mmmpi,N,nnew)
     &                        +u(2:Lmmpi+1,1:Mmmpi,N,nnew))

      DO ja = 0, nmaxatm
         if ( ssnd_nid(ja,i) /= -1 ) then
            call prism_put_proto(ssnd_nid(ja,i), sndtime, 
     &        snd_field, comp_ierror)
            if ( (comp_ierror == OASIS_Sent   ) .or.
     &           (comp_ierror == OASIS_ToRest ) .or.
     &           (comp_ierror == OASIS_SentOut) .or.
     &           (comp_ierror == OASIS_ToRestOut) ) then
           MPI_master_only write(*,*) 'CPL-CROCO : ',ja,
     & ssnd_clname(i),' sent at time ',sndtime
            endif
         endif
      enddo

         endif
         
         if (ssnd_clname(i)=='SRMVOCE'//clevel) then
            snd_field = 0.5 * (v(1:Lmmpi,1:Mmmpi  ,N,nnew)
     &                        +v(1:Lmmpi,2:Mmmpi+1,N,nnew))

      DO ja = 0, nmaxatm
         if ( ssnd_nid(ja,i) /= -1 ) then
            call prism_put_proto(ssnd_nid(ja,i), sndtime, 
     &        snd_field, comp_ierror)
            if ( (comp_ierror == OASIS_Sent   ) .or.
     &           (comp_ierror == OASIS_ToRest ) .or.
     &           (comp_ierror == OASIS_SentOut) .or.
     &           (comp_ierror == OASIS_ToRestOut) ) then
           MPI_master_only write(*,*) 'CPL-CROCO : ',ja,
     & ssnd_clname(i),' sent at time ',sndtime
            endif
         endif
      enddo

         endif
      END DO
#else 
!---------------------------------------------------------------------------------------
!     Dummy subroutine                                                           NO Coupling
!---------------------------------------------------------------------------------------
#endif 
      END SUBROUTINE cpl_prism_put
      

      
