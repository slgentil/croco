! $Id: zoombc_3Dfast.F 1458 2014-02-03 15:01:25Z gcambon $
!
!======================================================================
! CROCO is a branch of ROMS developped at IRD and INRIA, in France
! The two other branches from UCLA (Shchepetkin et al) 
! and Rutgers University (Arango et al) are under MIT/X style license.
! CROCO specific routines (nesting) are under CeCILL-C license.
! 
! This routine belongs to the specific CROCO package.
! 
! CROCO website : http://www.croco-ocean.org
!======================================================================
!
!
#include "cppdefs.h"

#if defined AGRIF && defined NBQ
!
!====================================================================
!                   subroutine unbq_bc_interp_tile
!====================================================================
!
      subroutine unbq_bc_interp_tile(Istr,Iend,Jstr,Jend)
      use AGRIF_Util
      implicit none
# include "param.h"
# include "boundary.h"  
# include "climat.h"
# include "grid.h"
# include "scalars.h"
# include "zoom.h"
# include "ocean3d.h"
# include "nbq.h"

      integer Istr,Iend,Jstr,Jend, i,j,k
      real tinterp
      INTEGER :: nbstep3dparent
      External :: u3dinterp
      integer :: parentnbstep
      !
# include "compute_auxiliary_bounds.h"
!
# ifdef MPI
#  define LOCALLM Lmmpi
#  define LOCALMM Mmmpi
# else
#  define LOCALLM Lm
#  define LOCALMM Mm
# endif

      call u3dbc_interp_tile(Istr,Iend,Jstr,Jend)


      parentnbstep=Agrif_Parent(iif)
      if (U3DFastTimeindex .NE. parentnbstep) then ! <<< - 3Dfast time index ok?
        dUinterp = 0.
        tinterp=1.
# ifdef MASKING
        Agrif_UseSpecialValue = .true.
# endif
        Agrif_SpecialValue = 0.

        common_index = 1
!        Call Agrif_Bc_variable(qdmunbqid,             ! <<< - qdmunbq OK?
!     &    calledweight=tinterp,procname = u3dinterp)  ! <<< - implement 3dfast == 3d for future?

        Agrif_UseSpecialvalue=.false.
        U3DFastTimeindex = parentnbstep
        U3DFastTimeindex2 = agrif_nb_step()
      endif                     ! end if (U3DFastTimeindex .NE. parentnbstep)


      if (agrif_nb_step() .EQ. U3DFastTimeindex2) then

# ifdef AGRIF_OBC_SOUTH
! ---------------------------------
           if (SOUTHERN_EDGE) then
              do k=1,N      
                 do i=Istr,IendR
                    qdmu_nbq(i,Jstr-1,k)=u(i,Jstr-1,k,nstp)     !  *0.01
#  ifdef MASKING
     &                   *umask(i,Jstr-1)
#  endif
                    
                 enddo
              enddo
           endif         
# endif

# ifdef AGRIF_OBC_NORTH
! ---------------------------------
           if (NORTHERN_EDGE) then      
              do k=1,N      
                 do i=Istr,IendR
                    qdmu_nbq(i,Jend+1,k)=u(i,Jend+1,k,nstp)       !   *0.01 
#  ifdef MASKING
     &                   *umask(i,Jend+1)
#  endif
                 enddo
              enddo
           endif      
# endif
           

# ifdef AGRIF_OBC_WEST
! ---------------------------------
           if (WESTERN_EDGE) then      
              do k=1,N      
                 do j=JstrR,JendR
                    qdmu_nbq(Istr,j,k)=u(Istr,j,k,nstp)    !   *0.01 
#  ifdef MASKING
     &                   *umask(Istr,j)
#  endif
                 enddo
              enddo
           endif      
# endif
           
# ifdef AGRIF_OBC_EAST
! ---------------------------------
           if (EASTERN_EDGE) then      
              do k=1,N      
                 do j=JstrR,JendR
                    qdmu_nbq(Iend+1,j,k)=u(Iend+1,j,k,nstp)       !  *0.01 
#  ifdef MASKING
     &                   *umask(Iend+1,j)
#  endif
                 enddo
              enddo
           endif
# endif


      endif ! end if (agrif_nb_step() .EQ. U2DTimeindex2
      
!
! \TODO To be continued ???
! 3d fast interp?
!

      return
      end       
!
!====================================================================
!                   subroutine vnbq_bc_interp_tile
!====================================================================
!
      subroutine vnbq_bc_interp_tile(Istr,Iend,Jstr,Jend)
      use AGRIF_Util
!     implicit none
# include "param.h"
# include "grid.h"
# include "ocean2d.h"
# include "scalars.h"
# include "climat.h"
# include "boundary.h"
# include "zoom.h"
# include "coupling.h"
#  include "nbq.h"
      integer Istr,Iend,Jstr,Jend, i,j,k
      external :: v3dinterp
!$AGRIF_DO_NOT_TREAT
      integer :: indinterp,nbgrid
      common/interp2d/indinterp,nbgrid
!$AGRIF_END_DO_NOT_TREAT      
  
# ifdef MPI
       include 'mpif.h'
#  define LOCALLM Lmmpi
#  define LOCALMM Mmmpi
# else
#  define LOCALLM Lm
#  define LOCALMM Mm
# endif   
!          
# include "compute_auxiliary_bounds.h"
!

      parentnbstep=Agrif_Parent(iif)
      
      if (V2DTimeindex .NE. parentnbstep) then
        dVinterp = 0.
        tinterp=1.

# ifdef MASKING        
        Agrif_UseSpecialValue = .true.
# endif
        Agrif_SpecialValue = 0.
        
        nbgrid = Agrif_Fixed()
        common_index = 1
!        Call Agrif_Bc_variable(qdmvnbqid,
!     &   calledweight=tinterp,procname = v3dinterp)  ! <<< - implement 3dfast == 3d for future?

        Agrif_UseSpecialvalue=.false.
        V2DTimeindex = parentnbstep
        V2DTimeindex2 = agrif_nb_step()
      endif

      if (agrif_nb_step() .EQ. V2DTimeindex2) then


# ifdef AGRIF_OBC_SOUTH
! ---------------------------------
      if (SOUTHERN_EDGE) then
         do k=1,N      
            do i=IstrR,IendR
               qdmv_nbq(i,Jstr,k)=qdmv_nbq(i,Jstr,k)       *0.01 
#  ifdef MASKING
     &              *vmask(i,Jstr) 
#  endif
            enddo
         enddo
      endif
# endif

# ifdef AGRIF_OBC_NORTH
! ---------------------------------
      if (NORTHERN_EDGE) then      
         do k=1,N      
            do i=IstrR,IendR
               qdmv_nbq(i,Jend+1,k)=qdmv_nbq(i,Jend+1,k)       *0.01    
#  ifdef MASKING
     &              *vmask(i,Jend+1)
#  endif
            enddo
         enddo
      endif      
# endif

# ifdef AGRIF_OBC_WEST
! ---------------------------------
      if (WESTERN_EDGE) then      
         do k=1,N      
            do j=Jstr,JendR
               qdmv_nbq(Istr-1,j,k)=qdmv_nbq(Istr-1,j,k)       *0.01     
#  ifdef MASKING
     &              *vmask(Istr-1,j)
#  endif
            enddo
         enddo
      endif      
# endif

# ifdef AGRIF_OBC_EAST
! ---------------------------------
      if (EASTERN_EDGE) then      
         do k=1,N      
            do j=Jstr,JendR
               qdmv_nbq(Iend+1,j,k)=qdmv_nbq(Iend+1,j,k)      *0.01 
#  ifdef MASKING
     &           *vmask(Iend+1,j)
#  endif
         enddo
      enddo
      endif
# endif

      endif !         if (agrif_nb_step() .EQ. V2DTimeindex2) then

      
      return
      end       
!
!====================================================================
!                   subroutine wnbq_bc_interp_tile
!====================================================================
!
      subroutine wnbq_bc_interp_tile(Istr,Iend,Jstr,Jend)
      use AGRIF_Util
!     implicit none
# include "param.h"
# include "grid.h"
# include "ocean2d.h"
# include "scalars.h"
# include "climat.h"
# include "boundary.h"
# include "zoom.h"
# include "coupling.h"
#  include "nbq.h"
      integer Istr,Iend,Jstr,Jend,i,j,k
      integer :: parentnbstep
      INTEGER :: nbstep3dparent
      real tinterp
      external t3dinterp
!
#  include "compute_auxiliary_bounds.h"
!
#  ifdef MPI
#   define LOCALLM Lmmpi
#   define LOCALMM Mmmpi
#  else
#   define LOCALLM Lm
#   define LOCALMM Mm
#  endif   

      nbstep3dparent=Agrif_Parent(nbstep3d)
      parentnbstep=Agrif_Parent_Nb_Step()

      if (TTimeindex.NE.parentnbstep) then
        tinterp=1.
# ifdef MASKING        
        Agrif_UseSpecialvalue=.true.
# endif
        Agrif_Specialvalue=0.
        Call Agrif_Set_bc(tid,(/-1,0/),
     &     InterpolationShouldbemade=.TRUE.)
!        Call Agrif_Bc_variable(tid,calledweight=tinterp,
!     &       procname = t3dinterp)
        
        Agrif_UseSpecialValue=.false.
      endif


# ifdef AGRIF_OBC_SOUTH
! ------------------------
         if (SOUTHERN_EDGE) then
            do k=1,N      
               do i=IstrR,IendR
                  qdmw_nbq(i,Jstr-1,k)=qdmw_nbq(i,Jstr-1,k)  *0.01 
#  ifdef MASKING     
     &        *rmask(i,Jstr-1)
#  endif
          enddo
        enddo
      endif
# endif

# ifdef AGRIF_OBC_NORTH
! ------------------------
      if (NORTHERN_EDGE) then      
        do k=1,N      
          do i=IstrR,IendR
             qdmw_nbq(i,Jend+1,k)=qdmw_nbq(i,Jend+1,k)   *0.01     
#  ifdef MASKING     
     &          *rmask(i,Jend+1)
#  endif
          enddo
        enddo
      endif      
# endif

# ifdef AGRIF_OBC_WEST
! -----------------------
      if (WESTERN_EDGE) then      
        do k=1,N      
          do j=JstrR,JendR
            qdmw_nbq(Istr-1,j,k)=qdmw_nbq(Istr-1,j,k)        *0.01 
#  ifdef MASKING     
     &         *rmask(Istr-1,j)
#  endif
          enddo
        enddo
      endif      
# endif

# ifdef AGRIF_OBC_EAST
! ------------------------
      if (EASTERN_EDGE) then      
         do k=1,N      
            do j=JstrR,JendR
               qdmw_nbq(Iend+1,j,k)=qdmw_nbq(Iend+1,j,k)     *0.01     
#  ifdef MASKING
     &              *rmask(Iend+1,j)
#  endif
            enddo
         enddo
      endif
# endif
!
! \TODO To be continued
!
!

      return
      end     
!
!====================================================================
!                   subroutine rnbq_bc_interp_tile
!====================================================================
!
      subroutine rnbq_bc_interp_tile(Istr,Iend,Jstr,Jend)
!
! \TODO To be continued
!
!
      integer indx,itrc
      indx=indxR              ! \TODO check index
      itrc=1                  ! \TODO set index
      !call t3dbc_interp_tile(Istr,Iend,Jstr,Jend,indx,itrc)
      return
      end       
!
#else
      subroutine zoombc_3Dfast_empty()
      return
      end
#endif
