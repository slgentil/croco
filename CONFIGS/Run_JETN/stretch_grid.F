#include "cppdefs.h"

      subroutine stretch_grid(dymin,dy,nbpts_tot)

#ifdef Y_STRETCH_GRID
	
      implicit none

#include "param.h"
#include "scalars.h"
#include "grid.h"
#include "mpi_roms.h"

        integer :: nbpts_L1, nbpts_L2, nbpts_tot,j,ierr
        real :: dymin
        real :: a1,b1,a3,b3,dummy,mapping,coeff_stretch
        real, dimension(-1:nbpts_tot+2) :: dy_tot, yp_tot
        real, dimension(GLOBAL_1D_ARRAYETA) :: dy

	nbpts_L2 = jet_domain_width / dymin
	nbpts_L1 = (nbpts_tot - nbpts_L2 ) / 2       
        
       ! init stretched coeff for south and north grid
#if RESOLUTION == 8	
       coeff_stretch = 1.08009
#elif RESOLUTION == 4	
!      coeff_stretch = 1.09741    ! 576 points
       coeff_stretch = 1.03973    ! 640 points
!      coeff_stretch = 1.02172    ! 704 points
!      coeff_stretch = 1.00874    ! 832 points
#elif RESOLUTION == 2
       coeff_stretch = 1.04802
#elif RESOLUTION == 1
       coeff_stretch = 1.04982
#endif

       ! compute dy : dy(n+1) = coeff_stretch*dy(n)
       do j=-1,nbpts_tot+2			  ! middle regular grid
           dy_tot(j)=dymin
       enddo
       do j=nbpts_L1+1,-1,-1			  ! south stretched grid
           dy_tot(j)=dy_tot(j+1)*coeff_stretch
       enddo
       do j=nbpts_tot-nbpts_L1,nbpts_tot+2	  ! north stretched grid
           dy_tot(j)=dy_tot(j-1)*coeff_stretch
       enddo

       ! compute y 
       yp_tot(1)=0
       do j=2,nbpts_tot+2			  ! total domain
	      yp_tot(j) = yp_tot(j-1)+dy_tot(j-1)
       enddo
       yp_tot(0) = yp_tot(1) - dy_tot(0)
       yp_tot(-1) = yp_tot(0) - dy_tot(-1)
       do j=-1,Mm+2				  ! local subdomain
           yp(:,j) = yp_tot(j+jminmpi-1)
           dy(j) = dy_tot(j+jminmpi-1)
       enddo  
       el = yp_tot(nbpts_tot+1)      

       if (mynode==0) then
 	  print*,'dy_max (m) = ',dy_tot(1)
 	  print*,'Taille du domaine (m) = ',yp_tot(nbpts_tot)
       endif
       
#endif 
       return
       end

