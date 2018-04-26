#include "cppdefs.h"

      subroutine  vtforc(Istr,Iend,Jstr,Jend)
       
c  forcage temporellement variable  de v et rho aux frontieres Nord et Sud
c  attention le choix est que l'integrale en z de vtmp n'est pas nulle donc
c  il faut injecter et faire sortir la meme masse aux frontieres Nord et Sud

      implicit none
# include "param.h"
# include "grid.h"
# include "forces.h"
# include "mixing.h"
# include "scalars.h"
# include "ocean3d.h"

      integer Istr,Iend,Jstr,Jend, i, j, k

      real h0,zmbnd,njours,njours_p, dnjours, sumamp
      real xtmp0,vtmp0,vtmpmax
      
      real onedeg,onedeg0,beta0
      integer L0,M0

      integer ni

#ifdef VTFORC_CORNER
      
      h0=h(Istr,Jstr)
c      zmbnd = pi/h0			!mde 1
      zmbnd = 2*pi/h0		!mde 2
c      zmbnd = 4*pi/h0		!mde 4
      gorho = g/rho0
      
c      njours = 24.      
      njours = 30    
c      njours = 40.
c      njours = 57.  
c      njours = 74.      

      njours_p = 57.    ! superieure a njours
      dnjours= (njours_p -njours)/(nb_ni -1)

      xtmp0 = 1.         ! en degrÃ©
      vtmp0 = 0.2      ! en m/sec
c      vtmp0 = 0.1      ! en m/sec
c      vtmp0 = 0.05      ! en m/sec
        do ni=1,nb_ni
	   freqtmp(ni)=2.*pi/((   njours+(ni-1.)*dnjours)*86400.)
	   phitmp(ni)=2.*pi*ni/nb_ni
        enddo
        sumamp=0.
        do ni=1,nb_ni
	   sumamp= sumamp+ freqtmp(ni)**2.
        enddo
        do ni=1,nb_ni
           ampltmp(ni)=freqtmp(ni)/sqrt(sumamp)
        enddo
      if (iic.eq.ntstart) then
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       MPI_master_only print *, '   '
       MPI_master_only print*,'======================================='
       MPI_master_only print *,'FORCAGE v et rho aux bords Nord et Sud'
       MPI_master_only print *,' dnjours ', dnjours
       MPI_master_only print *,' periode ', 2.*pi/freqtmp(:)/86400.
       MPI_master_only print *,' periode njours min=',
     &   2.*pi/freqtmp(1)/86400.
       MPI_master_only print *,' periode njours max=',
     &  2.*pi/freqtmp(nb_ni)/86400.
       MPI_master_only print *,' nbre echantillonnage freq=',nb_ni
       MPI_master_only print *,' largeur courant xtmp0=',xtmp0      
       MPI_master_only print *,' amp vitesse meridienne vtmp0=',vtmp0
       MPI_master_only print*,'================================='
       MPI_master_only print *, '   '
       MPI_master_only print*,'============================='
       MPI_master_only print *, '!FORCAGE cos(mz) 2*pi/m =',2*pi/zmbnd
       MPI_master_only print*,'================================'
       MPI_master_only print *, '   '
      endif
     
      do k=1,N
        do j=Jstr-1,Jend+1
          do i=Istr-1,Iend+1
            vtmp(i,j,k)=xr(i,j)/xtmp0*exp(-0.5*(xr(i,j)/xtmp0)**2.)  
     &     *cos(zmbnd*z_r(i,j,k))
          enddo
        enddo
      enddo
      
!      vtmpmax=exp(-0.5)
      vtmpmax=  0.5971881681267380
      
#ifndef MPI      
      print *,' vtmpmax=',vtmpmax 
      L0=LLm+1
      M0=MMm+1
      print *,'L0 M0=',L0,M0
      print *,'xr(1,  et 2)=',xr(1,M0),xr(2,M0)
      print *,'xr(L  et Lm)=',xr(L0,M0),xr(L0-1,M0)
      print *,'yr(1,  et 2)=',yr(1,1),yr(1,2)
      print *,'yr(M  et Mm)=',yr(L0,M0),yr(L0,M0-1)

      print *,'xp(1,  et 2)=',xp(1,M0),xp(2,M0)
      print *,'xp(L  et Lm)=',xp(L0,M0),xp(L0-1,M0)
      print *,'yp(1,  et 2)=',yp(1,1),yp(1,2)
      print *,'yp(M  et Mm)=',yp(L0,M0),yp(L0,M0-1)
    
#endif /* MPI */
      
      vtmp(:,:,:)=vtmp0*vtmp(:,:,:)/vtmpmax
      onedeg0=1./180 *pi*6.37e6    ! valeur codée dans spem     
      beta0=2.2829e-11        ! valeur codée dans spem
      onedeg=pi*Eradius/180
      
      do k=1,N
         do j=Jstr-1,Jend+1
            do i=Istr-1,Iend+1
               ttmp(i,j,k) = - (     	! sert pour rappel bandes NordSud
     &      beta0*yr(i,j)*onedeg   		! changement signe vs Marc
     &     *1/gorho
     &     *vtmp0/vtmpmax
     &     *xtmp0*onedeg*exp(-0.5*(xr(i,j)/xtmp0)**2.)  
     &     *(-1.)*zmbnd*sin(zmbnd*z_r(i,j,k))  ) / Tcoef + T0
            enddo
         enddo
      enddo
      
#ifndef MPI
      print *,'onedeg= ', onedeg,'onedeg0= ', onedeg0
      print *,' beta0= ',beta0, ' beta=',beta
      print *,' f(1,M)=', f(1,M0)
      print *,' f(1,1)=', f(1,1)
      
      print *,' maxval(ttmp(1,:,:)= ', maxval(ttmp(1,:,:))

      vtmpmax=maxval(vtmp(:,2,:))

!      print *,'vtmp(:,35)=',vtmp(:,2,35)
!      print *,'vtmp(4,:)= ',vtmp(4,2,:)
#endif /* MPI */
           
     
#endif /* VTFORC_CORNER */
      return
      end
