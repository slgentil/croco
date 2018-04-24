#include "cppdefs.h"

c *******************************************************************
c  Initialize a bunch of stuff.  This is called before writing out
c  the history header file or else after it is read in.  This allows
c  you to change the values of variable such as dt.
c *******************************************************************

      subroutine init_sponge(Istr,Iend,Jstr,Jend)

      implicit none
# include "param.h"
# include "grid.h"
# include "forces.h"
# include "mixing.h"
# include "scalars.h"
# include "ocean3d.h"

      integer Istr,Iend,Jstr,Jend, i, j, k, isp, itrc, ibnd

#ifdef SPONGE_LAYER
      real         yeq_pv_nord, ymd_pv_nord, ynr_pv_nord
      real         yeq_pv_sud, ymd_pv_sud, ynr_pv_sud, dzforc
      real         zforc_eq_nord, zforc_md_nord, zforc_nr_nord
      real         zforc_eq_sud, zforc_md_sud, zforc_nr_sud
      real         yspng, dyspng
      real         xspng, dxspng
      real         rdmp_bnd, rdmp
#endif  /* SPONGE_LAYER */

# include "compute_extended_bounds.h"

c  nudging time constant in units of sec**-1.

#if defined  SPONGE_LAYER
c ai rajoute la possibilite de mettre la couche de melange que sur une fraction
c septentrionale du domaine
      rdmp = 1./(10.*86400.)
      rdmp_bnd = 1./(10.*86400)
      yspng = 35.
      yspng = 12.      
      dyspng = 1.
      xspng = 48.
      dxspng = 1.
      dzforc = 0.0001
      

      if (iic.eq.ntstart) then
      
      MPI_master_only print *,' '
      MPI_master_only print *,'parameters for PV-FORCED LAYER'
      MPI_master_only print *,'rdmp = ',rdmp
      MPI_master_only print *,'rdmp-bnd = ',rdmp_bnd
      MPI_master_only print *,'yspng = ',yspng
      MPI_master_only print *,'dyspng = ',dyspng 
      MPI_master_only print *,'xspng = ',xspng
      MPI_master_only print *,'dxspng = ',dxspng
      MPI_master_only print *,'dzforc = ',dzforc
      MPI_master_only print *,' '
      
      endif


c rdmp_mxlyr for the northern region

      do j = JstrR,JendR
         if (     (yr(1,j).ge. yspng+dyspng) 
     &        .or.(yr(1,j).le.-yspng-dyspng) ) then
            rdmp_mxlyr(:,j,:) = rdmp_bnd 
         elseif (     (yr(1,j).le. yspng-dyspng) 
     &           .and.(yr(1,j).ge.-yspng+dyspng)) then
            rdmp_mxlyr(:,j,:) = 0.
         elseif (     (yr(1,j).gt. yspng-dyspng) 
     &           .and.(yr(1,j).lt. yspng+dyspng) ) then
            rdmp_mxlyr(:,j,:) = rdmp_bnd * 0.5 *
     &        (1.+sin(0.5*pi*(yr(1,j)-yspng)/dyspng))
         elseif (     (yr(1,j).gt.-yspng-dyspng) 
     &           .and.(yr(1,j).lt.-yspng+dyspng) ) then
            rdmp_mxlyr(:,j,:) = rdmp_bnd * 0.5 *
     &        (1.-sin(0.5*pi*(yr(1,j)+yspng)/dyspng))
         else
            print *,'Probleme rdmp_mxlyr'
            stop
         endif
      enddo
 
#endif /* SPONGE_LAYER */


#ifdef VTFORC_INT
      if (iic.eq.ntstart) then
         call vtforc(Istr,Iend,Jstr,Jend)
      endif
#endif  /* VTFORC_INT */   


      return
      end
