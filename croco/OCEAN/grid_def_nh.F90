#include "cppdefs.h"
#ifdef NBQ

      subroutine grid_def_nh

!******************************************************************************************
!                   Numbering of Mass Points
!
!******************************************************************************************

      use module_nh
      use module_nbq

      implicit none

# include "param_F90.h"
# include "scalars_F90.h"
# include "grid.h"

# include "def_bounds.h"

      integer :: i1,i2,j1,j2,i,j

!******************************************************************************************
! Inner MPI domain:
! cf: http://poc.obs-mip.fr/auclair/WOcean.fr/SNH/Restricted/NH-NBQ/Html_pages/Algebrique_Numerotation_Base.htm
!******************************************************************************************

!*******************************************************************
!     Definitions: *_INTER_NBQ logical variables
!*******************************************************************


# ifdef MPI 

       WEST_INTER_NBQ  = WEST_INTER
       EAST_INTER_NBQ  = EAST_INTER
       SOUTH_INTER_NBQ = SOUTH_INTER
       NORTH_INTER_NBQ = NORTH_INTER

#  ifdef EW_PERIODIC
       if (WESTERN_EDGE) then
!         WESTERN_EDGE   = .FALSE.
          WEST_INTER_NBQ = .TRUE.
       endif
       if (EASTERN_EDGE) then
!         EASTERN_EDGE   = .FALSE.
          EAST_INTER_NBQ = .TRUE.
       endif
#  endif

#  ifdef NS_PERIODIC
       if (SOUTHERN_EDGE) then
!         SOUTHERN_EDGE   = .FALSE.
          SOUTH_INTER_NBQ = .TRUE.
       endif
       if (NORTHERN_EDGE) then
!         NORTHERN_EDGE   = .FALSE.
          NORTH_INTER_NBQ = .TRUE.
       endif
#  endif

# endif 


!*******************************************************************
!     NH domain:  mass grid-points
!*******************************************************************

       istr_nh  = 1
       iend_nh  = LOCALLM

       jstr_nh  = 1
       jend_nh  = LOCALMM 

!......The following coef. are updated in nump_nh:
       istrq_nh = istr_nh
       jstrq_nh = jstr_nh
       iendq_nh = iend_nh 
       jendq_nh = jend_nh

!*******************************************************************
!     NH domain:  U-V velocity grid-points 
!         (coef. updated in numuvw_nh)
!*******************************************************************

       istru_nh = 2
       iendu_nh = LOCALLM 
       jstru_nh = jstr_nh
       jendu_nh = jend_nh
       istrv_nh = istr_nh
       iendv_nh = iend_nh 
       jstrv_nh = 2
       jendv_nh = LOCALMM 

# ifdef MPI 
      if (WEST_INTER_NBQ) then
       istru_nh = 1
      endif

      if (SOUTH_INTER_NBQ) then
       jstrv_nh = 1
      endif
# endif

# ifdef MASKING
!*******************************************************************
! Grid-mask for NBQ-use
!*******************************************************************

      rmask_nbq(Istr_nh:Iend_nh  ,Jstr_nh:Jend_nh)   = rmask(Istr_nh:Iend_nh  ,Jstr_nh:Jend_nh)
      umask_nbq(Istru_nh:Iendu_nh,Jstru_nh:Jendu_nh) = umask(Istru_nh:Iendu_nh,Jstru_nh:Jendu_nh)
      vmask_nbq(Istrv_nh:Iendv_nh,Jstrv_nh:Jendv_nh) = vmask(Istrv_nh:Iendv_nh,Jstrv_nh:Jendv_nh)

!  Exchange Grid
#  ifdef MPI
      call exchange_r2d_tile (Istr_nh,Iend_nh,Jstr_nh,Jend_nh,  &
                        rmask_nbq(START_2D_ARRAY))
      call exchange_u2d_tile (Istru_nh,Iendu_nh,Jstru_nh,Jendu_nh,  &
                        umask_nbq(START_2D_ARRAY))
      call exchange_v2d_tile (Istrv_nh,Iendv_nh,Jstrv_nh,Jendv_nh,  &
                        vmask_nbq(START_2D_ARRAY))
#  endif

      i1 = istr_nh
      i2 = iend_nh
      do j=jstr_nh-1,jend_nh+1
#  ifdef MPI
         if (WESTERN_EDGE) rmask_nbq(i1-1,j) = rmask_nbq(i1,j)
         if (EASTERN_EDGE) rmask_nbq(i2+1,j) = rmask_nbq(i2,j)
#  else
#   if defined OBC_WEST
            rmask_nbq(i1-1,j) = rmask_nbq(i1,j)
#   endif
#   if defined OBC_EAST
            rmask_nbq(i2+1,j) = rmask_nbq(i2,j)
#   endif
#  endif
      enddo

      j1 = jstr_nh
      j2 = jend_nh
      do i=istr_nh-1,iend_nh+1
#  ifdef MPI
        if (SOUTHERN_EDGE) rmask_nbq(i,j1-1) = rmask_nbq(i,j1)
        if (NORTHERN_EDGE) rmask_nbq(i,j2+1) = rmask_nbq(i,j2)
#  else
#   if defined OBC_SOUTH
           rmask_nbq(i,j1-1) = rmask_nbq(i,j1)
#   endif
#   if defined OBC_NORTH
           rmask_nbq(i,j2+1) = rmask_nbq(i,j2)
#   endif
#  endif
      enddo

      i1 = istru_nh
      i2 = iendu_nh
      do j=jstru_nh-1,jendu_nh+1
#  ifdef MPI
         if (WESTERN_EDGE) umask_nbq(i1-1,j) = umask_nbq(i1,j)
         if (EASTERN_EDGE) umask_nbq(i2+1,j) = umask_nbq(i2,j)
#  else
#   if defined OBC_WEST
            umask_nbq(i1-1,j) = umask_nbq(i1,j)
#   endif
#   if defined OBC_EAST
            umask_nbq(i2+1,j) = umask_nbq(i2,j)
#   endif
#  endif
      enddo

      j1 = jstru_nh
      j2 = jendu_nh
      do i=istru_nh-1,iendu_nh+1
#  ifdef MPI
        if (SOUTHERN_EDGE) umask_nbq(i,j1-1) = umask_nbq(i,j1)
        if (NORTHERN_EDGE) umask_nbq(i,j2+1) = umask_nbq(i,j2)  
#  else
#   if defined OBC_SOUTH
           umask_nbq(i,j1-1) = umask_nbq(i,j1)
#   endif
#   if defined OBC_NORTH
           umask_nbq(i,j2+1) = umask_nbq(i,j2)
#   endif
#  endif
      enddo

      i1 = istrv_nh
      i2 = iendv_nh
      do j=jstrv_nh-1,jendv_nh+1
#  ifdef MPI
         if (WESTERN_EDGE) vmask_nbq(i1-1,j) = vmask_nbq(i1,j)
         if (EASTERN_EDGE) vmask_nbq(i2+1,j) = vmask_nbq(i2,j)
#  else
#   if defined OBC_WEST
            vmask_nbq(i1-1,j) = vmask_nbq(i1,j)
#   endif
#   if defined OBC_EAST
            vmask_nbq(i2+1,j) = vmask_nbq(i2,j)
#   endif
#  endif
      enddo

      j1 = jstrv_nh
      j2 = jendv_nh
      do i=istrv_nh-1,iendv_nh+1
#  ifdef MPI
        if (SOUTHERN_EDGE) vmask_nbq(i,j1-1) = vmask_nbq(i,j1)
        if (NORTHERN_EDGE) vmask_nbq(i,j2+1) = vmask_nbq(i,j2)
#  else
#   if defined OBC_SOUTH
           vmask_nbq(i,j1-1) = vmask_nbq(i,j1)
#   endif
#   if defined OBC_NORTH
           vmask_nbq(i,j2+1) = vmask_nbq(i,j2)
#   endif
#  endif
      enddo
# endif

      return

      end subroutine grid_def_nh
 
#else
        subroutine grid_def_nh_empty
        return
        end 
#endif
