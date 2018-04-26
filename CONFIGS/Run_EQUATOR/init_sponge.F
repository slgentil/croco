#include "cppdefs.h"

c *******************************************************************
c  Initialize a bunch of stuff.  This is called before writing out
c  the history header file or else after it is read in.  This allows
c  you to change the values of variable such as dt.
c *******************************************************************

      subroutine init_sponge(Istr,Iend,Jstr,Jend)

#ifdef VTFORC_CORNER

      implicit none
# include "param.h"
# include "grid.h"
# include "forces.h"
# include "mixing.h"
# include "scalars.h"
# include "ocean3d.h"

      integer Istr,Iend,Jstr,Jend, i, j, k, isp, itrc, ibnd

      real         yspng, dyspng
      real         rdmp_bnd, rdmp

# include "compute_extended_bounds.h"

c  nudging time constant in units of sec**-1.

      rdmp = 1./(10.*86400.)
      rdmp_bnd = 1./(10.*86400)
      yspng = 12.      
      dyspng = 1.
      

      if (iic.eq.ntstart) then    
      MPI_master_only print *,' '
      MPI_master_only print *,'parameters for PV-FORCED LAYER'
      MPI_master_only print *,'rdmp = ',rdmp
      MPI_master_only print *,'rdmp-bnd =',1./(rdmp_bnd*86400),' days'
      MPI_master_only print *,'yspng = ',yspng
      MPI_master_only print *,'dyspng = ',dyspng 
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


      if (iic.eq.ntstart) then
         call vtforc(Istr,Iend,Jstr,Jend)
      endif

#endif  /* VTFORC_CORNER */   


      return
      end
