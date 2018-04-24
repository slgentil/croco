#include "cppdefs.h"
 
      subroutine init_rhobar(istr,iend,jstr,jend)
!
! Set initial conditions for momentum and tracer variables using
!  analytical expressions.
!
      implicit none
      integer tile,istr,iend,jstr,jend
      integer i,j,k, itrc
      integer L,M
      real x,y
      
# include "param.h"
# include "forces.h"
# include "grid.h"
# include "ocean2d.h"
# include "ocean3d.h"
# include "coupling.h"
# include "scalars.h"
# include "work.h"
# include "private_scratch.h"

#if !defined MPI
      integer ii,jj
      parameter (ii=0, jj=0)
#endif 

!-- for parallelization (integration of zeta in eta direction)
#ifdef MPI
      include 'mpif.h'
#endif

c# ifdef BASIN
c      real zze, dze, drho, xnconst, xn0, xn1, ztherm, gorho
c# endif

# ifdef BV_VAR
      real xn2, dz1, z1, z2, drhotc, dz2, aa, bb
# endif
      

# include "compute_auxiliary_bounds.h"
!

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

! Initial conditions for free surface and 2D momentum components.
!======== ========== === ==== ======= === == ======== ===========

# if defined BASIN_EQ

!--------------------------------------------------------------------
!  Set initial conditions for potential temperature [degC] and
!  salinity [PSU].
!
# ifdef BV_VAR
      if (.false.) then
! Fred's formula
      zze = 150.
      dze = 50.
      drho = 2.3
      xn0 = 2.2e-3
      xn1 = 0.2e-3
      ztherm=1000.
      gorho = g/rho0
                                                                          
      do k=1,N
        do j=JR_RANGE
          do i=IR_RANGE
c            rhobar(i,j,k) = -(xnconst)**2. * (z_r(i,j,k)) / gorho
           rhobar(i,j,k)=drho * (1. - tanh((z_r(i,j,k)+zze)/dze)) / 2.
     &       - xn0**2 *(exp(z_r(i,j,k)/ztherm))   * rho0 * ztherm / g
     &       - xn1**2 * rho0 * z_r(i,j,k) / g
          enddo
        enddo
      enddo

      MPI_master_only print *,'=============================='
      MPI_master_only print *, '   '
      MPI_master_only print *, 'zze=',zze,'  dze=',dze,'  xn0=',xn0,
     &        ' ztherm=',ztherm, ' xn1=',xn1
      MPI_master_only print *, '   '
      MPI_master_only print *,'=================================='

      else
! Marc's formula
      print *,'inside init_rhobar Marc'
      gorho = g/rho0
      xn1=sqrt(10.)*1e-3                  
      xn2=sqrt(6e-7)
      dz1=300
      z1=1200
      z2=150
      drhotc=2.9
      dz2=70.      
      aa=(xn2**2-xn1**2)/2/gorho
      bb=(xn2**2+xn1**2)/2/gorho

      MPI_master_only print *,'=============================='
      MPI_master_only print *,'xn1,xn2,dz1,z1,z2,drhotc,dz2,aa,bb'
      MPI_master_only print *, xn1,xn2,dz1,z1,z2,drhotc,dz2,aa,bb       
      MPI_master_only print *,'=============================='

      do k=1,N
        do j=JR_RANGE
          do i=IR_RANGE
      	  rhobar(i,j,k) = aa*dz1*log(cosh(-(-z_r(i,j,k)-z1)/dz1))
     &                +bb*(-z_r(i,j,k)-z1)
     &                -(aa*dz1*log(cosh(-( -z1)/dz1))+bb*( -z1))
     &                + drhotc*(1+tanh((-z_r(i,j,k)-z2)/dz2))/2.
          enddo
        enddo
      enddo

      MPI_master_only print*,'=============================='
      MPI_master_only print *, '   '
      MPI_master_only print *, 
     &          'Variable stratification (Marc) mean BV=1.8e-3'
      MPI_master_only print *, '   '
      MPI_master_only print *,'=================================='
      endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# else
      xnconst = 2.e-3
      gorho = g/rho0
                                                                          
      do k=1,N
        do j=JR_RANGE
          do i=IR_RANGE
            rhobar(i,j,k) = -(xnconst)**2. * (z_r(i,j,k)) / gorho
            tbar(i,j,k) = - (rhobar(i,j,k)) / Tcoef + T0
           enddo
        enddo
      enddo

      MPI_master_only print *,'xnconst=',xnconst

# endif  /* BV_VAR */
# endif /* BASIN_EQ */


# undef IR_RANGE
# undef IU_RANGE
# undef JR_RANGE
# undef JV_RANGE
      return
      end
