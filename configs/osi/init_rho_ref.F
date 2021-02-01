#include "cppdefs.h"

      subroutine init_rho_ref()
      
      implicit none

# include "param.h"
# include "ocean3d.h"
# include "diagnostics.h"

      integer::k
      real(kind=8) :: localsum  , globalsum, nr_sumofall

! Initialisation d'un rho référence: rho moyenné par tranche horizontal    
  
      do k=1,N
         localsum = sum(rho(1:Lm,1:Mm,k))/Lm/Mm
         globalsum = nr_sumofall(localsum)
         rho_r(k) = globalsum/NNODES
      enddo          
      rho_r = -1. * rho_r
      call qsort(rho_r,N)
      rho_r =  -1. * rho_r 

      do k=2,N-1
         drhordz(k) = 0.5*(rho_r(k+1) - rho_r(k-1)) / Hz(1,1,k)
      enddo
      drhordz(N) = drhordz(N-1)
      drhordz(1) = drhordz(2)

      end subroutine init_rho_ref

