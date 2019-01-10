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

      subroutine windinit(Istr,Iend,Jstr,Jend)
!
!----------------------------------------------------------------------
!  Set wind profile.
!----------------------------------------------------------------------
!
      implicit none
#include "param.h"
#include "grid.h"
#include "scalars.h"

      integer Istr,Iend,Jstr,Jend, i,j,k
      real xmid, xwide, xtap, ymid, ywide, ytap, x1, x2, y1,y2
      real(kind=8) :: localsum  , globalsum, nr_maxofall

#include "compute_auxiliary_bounds.h"

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

! xmid : middle of the bump
! xwide : space between both centered tanh
! xtap : length of the slope 
! unit = degree

      xmid = xl/3.
      xwide = 20.
      xtap = 5.
      ymid = 0.
      ywide = 6.
      ytap = 5.

      x1 = xmid - 0.5*xwide
      x2 = xmid + 0.5*xwide
      y1 = ymid - 0.5*ywide
      y2 = ymid + 0.5*ywide

      do i=IR_RANGE
         prof_taux(i) = tanh((xr(i,1)-x1)/xtap) - tanh((xr(i,1)-x2)/xtap)
      enddo
      localsum = maxval(prof_taux)
      globalsum = nr_maxofall(localsum)
      prof_taux = prof_taux / globalsum
      do j=JR_RANGE
!        prof_tauy(j) = tanh((yr(1,j)-y1)/ytap) - tanh((yr(1,j)-y2)/ytap)
!        Blackman window
         prof_tauy(j) = 0.42 - 0.5*cos(2.*pi*(yr(1,j)+el/2.)/el) + 
     &                         0.08*cos(4.*pi*(yr(1,j)+el/2.)/el)
      enddo
      localsum = maxval(prof_tauy)
      globalsum = nr_maxofall(localsum)
      prof_tauy = prof_tauy / globalsum

      return
      end
