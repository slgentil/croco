!
! $Id: set_dissip.F, 2005/04/07 $
!
#include "cppdefs.h"

      subroutine set_dissip(Istr,Iend,Jstr,Jend)
!
! Computes 3D horizontal dissipation coefficients for u,v and t
!
      implicit none
      
#include "param.h"
#include "grid.h"
#include "ocean3d.h"
#include "mixing.h"
#include "scalars.h"

      integer Istr,Iend,Jstr,Jend, i,j
#ifdef SOLVE3D
     &       ,k,itrc
#endif
      real uvnu_spng,tnu_spng
      real yspng,dyspng,xspng_west,dxspng_west,dum
      
# ifdef MPI
      include 'mpif.h'
      integer ierr
# endif
!
# include "compute_extended_bounds.h"

#if defined UVNU_MAT || defined TNU_MAT 

      do k=1,N
        do j=JstrR,JendR
          do i=IstrR,IendR
            uvnu(i,j,k) = visc2
            enddo
         enddo
      enddo

      do itrc=1,NT
      do k=1,N
        do j=JstrR,JendR
          do i=IstrR,IendR
            tnu(i,j,k,itrc) =  tnu2(itrc)
            enddo
         enddo
      enddo
      enddo 
                                                                                                                                                
#if defined  NS_SPONGE
      uvnu_spng = 5000.
      tnu_spng = 5000.     
      yspng = 12.
      dyspng = 1.
 
      do k=1,N
        do j=JstrR,JendR
          do i=IstrR,IendR
           uvnu(i,j,k) = max( uvnu(i,j,k), uvnu_spng *
     &      ( 2. + tanh((yr(i,j)-yspng)/dyspng)
     &           - tanh((yr(i,j)+yspng)/dyspng) ) / 2. )
            enddo
         enddo
      enddo

      do itrc=1,NT
      do k=1,N
        do j=JstrR,JendR
          do i=IstrR,IendR
           tnu(i,j,k,itrc) = max( tnu(i,j,k,itrc), tnu_spng *
     &      ( 2. + tanh((yr(i,j)-yspng)/dyspng)
     &           - tanh((yr(i,j)+yspng)/dyspng) ) / 2. )
          enddo
        enddo
      enddo
      enddo

      if (iic.eq.ntstart) then
        MPI_master_only print*,'--------------------------------'
        MPI_master_only print*,' northern boundary sponge layer '
        MPI_master_only print*,'--------------------------------'
        MPI_master_only print *,'uvnu_spng = ',uvnu_spng
        MPI_master_only print *,'yspng = ',yspng
        MPI_master_only print *,'dyspng = ',dyspng
        MPI_master_only print *,'tnu_spng = ',tnu_spng
      endif

#endif   /* NS_SPONGE */

#ifdef  WEST_SPONGE
      uvnu_spng = 1000.
      tnu_spng = 1000.
      xspng_west = 2.5
      dxspng_west = 1.
        
      do k=1,N
        do j=JstrR,JendR
          do i=IstrR,IendR
           uvnu(i,j,k) = max( uvnu(i,j,k),  uvnu_spng *
     &      (1.-tanh((xr(i,j)-xspng_west)/dxspng_west)) / 2. )
            enddo
         enddo
      enddo

       do itrc=1,NT
       do k=1,N
        do j=JstrR,JendR
          do i=IstrR,IendR
           tnu(i,j,k,itrc) = max( tnu(i,j,k,itrc),  tnu_spng *
     &      (1.-tanh((xr(i,j)-xspng_west)/dxspng_west)) / 2. )
          enddo
        enddo
      enddo
      enddo

      if (iic.eq.ntstart) then
        MPI_master_only print*,'--------------------------------'  
        MPI_master_only print *,' Western Boundary sponge on '
        MPI_master_only print*,'--------------------------------' 
        MPI_master_only print *,'uvnu_spng=',uvnu_spng
        MPI_master_only print *,'xspng_west = ',xspng_west
        MPI_master_only print *,'dxspng_west = ',dxspng_west
        MPI_master_only print *,' tnu_spng=',tnu_spng
#ifdef  SALINITY
        MPI_master_only print *,'!!!  on salt along NS boundaries'
#endif
      endif

#endif   /* WEST_SPONGE */

        MPI_master_only print *,'max(uvnu)=',maxval(uvnu)
        dum=dt * (pi**2) * maxval(uvnu) * pm(1,1)**2
        MPI_master_only print *,'CFL diffusion horizontale  laplacien =',dum

#endif /* UVNU_MAT or TNU_MAT */

      return
      end

