#include "cppdefs.h"
!=====================================================================
       subroutine matrice(xlmb,w2)
!=====================================================================
#ifdef GO_MEDDY
    !
!!$c history
!!$c double-checked on 5/2/2010 that barotropic eigenvalue ~e-13 for 120 levels
!!$c for gaussian lense of initpsi.F (icase=0) arrays akm and akp have been tested
!!$c for the analytical value of the ratio of vortex stretching to relative vorticity 
!!$c wich should be given by 1/(2*Burger)
!!$c beware that xkrad1 = Lx/(2*pi*xlambda)
!!$c where xlambda = first radius of deformation = Nfreq*Lz/(pi*f0)
!!$c the 2*pi factor in xkrad1 originates from integer values of wavenumber used fo
!!$r 
!!$c derivation and integration in horizontal space
!!$c
!!$c thus apart from suppression of prints output no line has been changed in this routine
!!$c the burden of appropriate definition of xkrad1 is taken care of in initpsi

       implicit none  
# include "param.h"
# include "grid.h"
# include "ocean3d.h"
# include "scalars.h"
       include "mpif.h"

       integer,parameter :: lwork=3*N
       integer,parameter :: liwork = 2*N  
       integer i,i1,j,j1,k,mde,ilo,ihi,ier,zl       
       integer,dimension(N) :: ior 
       integer,dimension(liwork) :: iwork       
       real(8) :: abnrm, xmin, xlambda, xkrad1, w2
       real(8),dimension(0:N+1) :: cofh, cofn
       real(8),dimension(N,N) :: zr,cd_tot
       real(8),dimension(N) :: wr,wii,cr,gp,gi,scal,rconde,rcondv, xlmb
       real(8),dimension(lwork) :: work
       
       
! N levels or modes

       do i=0,N+1
         cofh(i)=1./float(N)
         cofn(i)=1.
       enddo

       if (mynode==0) 
     &    print *,'top and bottom value cofn(N), cofn(1) =',cofn(N), cofn(1)

!!$c gp is the relative gravity at each interface starting at ocean surface
!!$c cofh(i) relative depth of layer i
!!$c cofn(i) Brunt-Vaissala value between layer i and i+1 normalized by Nfreq

       do i=0,N-1
         gp(i+1) = 0.5*(cofh(i)+cofh(i+1)) * cofn(i)**2
       enddo
       gi(1)=0.
       do i=2,N
         gi(i)=1./gp(i)
       enddo
    
!!$c---------------------------------------------
!!$c    CONSTRUCTION DE LA MATRICE DU VORTEX STRETCHING [b]
!!$c---------------------------------------------
       do j=1,N
         do i=1,N
           cd_tot(i,j)=0.
         enddo
       enddo

       do i=1,N
         cd_tot(i,i)=-gi(i)/cofh(i)
       enddo
       do i=2,N
         cd_tot(i-1,i)=-cd_tot(i,i)
         cd_tot(i,i-1)=gi(i)/cofh(i-1)
         cd_tot(i-1,i-1)=cd_tot(i-1,i-1)-cd_tot(i,i-1)
       enddo


!!$c---------------------------------------------
!!$c    Appel subroutine eig de Lapack (monoproc) 
!!$c---------------------------------------------

       call dgeevx('N','N','V','N',N,cd_tot,N,wr,wii, 
     &    zl,1,zr,N,ilo,ihi,scal, 
     &    abnrm,rconde,rcondv,work,lwork,iwork,ier)
       if (ier.ne.0) then
         print *,'erreur dans dgeevx  ier=',ier
         stop
       endif
       do  i=1,N
         cr(i)=abs(wr(i))
       enddo

!!$c------------------------------------------------------
!!$c    CLASSEMENTS DE VALEURS PROPRES PAR ORDRE CROISSANT
!!$c------------------------------------------------------

       xmin=1.0e+10
       j1=0
       do 31 k=1,N
         do 33 i=1,N
           if (k.eq.1) go to 35
           do 34 j=1,j1
             if (i.eq.ior(j)) go to 33
34         continue
35         continue
           xmin=min(cr(i),xmin)
33       continue
         do 37 i=1,N
            i1=i
            if (cr(i).eq.xmin) go to 38
37       continue
38       continue
         j1=j1+1
         ior(j1)=i1
         wii(j1)=cr(i1)
         xmin=1.0e+10
31     continue

              if (mynode==0) then
                print *,'barotropic eigenvalue=', wii(1)/wii(2)
                print *,'first 4 eig prior to norm.',wii(1),wii(2),wii(3),wii(4)
              endif
              cr(1)=0.
              w2 = wii(2)
              do  j=2,N
                 cr(j)=wii(j)/w2
              enddo

!!$c---------------------------------
!!$c radius of deformation scaling
!!$c---------------------------------

!                xlambda=sqrt(Nfreq2)*(Hz(1,1,1)*(N-1))/(f(1,1)*sqrt(w2) )
                xlambda=sqrt(Nfreq2)*(Hz(1,1,N/2)*(N-1))/(f(1,1)*pi )
                xkrad1=((LLm-1)/pm(1,1))/(xlambda*2*pi)
                if (mynode==0) then
                   print *,'~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
                   print *,' resolution dx=',1./pm(1,1),' dz=',Hz(1,1,N/2)
                   print *,' Coriolis parameter ', f(1,1), 
     &                     ' Nfreq=', sqrt(Nfreq2)
                   print *,' radius of deformation (km)',  xlambda/1e3
                   print *,'~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
                   print  *,' xkrad1 ', xkrad1
                endif


!!$c---------------------------------
!!$c normalisation des valeurs propres pour le code spectral
!!$c---------------------------------
                do mde=1,N
                   xlmb(mde)=cr(mde)*xkrad1**2
                enddo

                if (mynode == 0)  print *,' xlmb=',minval(xlmb),maxval(xlmb)
#endif /*  GO_MEDDY */
       end subroutine matrice


!=====================================================================
       subroutine init_spectral(planf,planb,rvar,cvar)
!=====================================================================

# include "param.h"
       include "fftw3.f"      

       integer*8 :: planf, planb
       real(kind=8),dimension(1:LLm,1:MMm) :: rvar    
       complex(kind=8),dimension(LLm/2+1,1:MMm) :: cvar
    !     
    !     To transform a two-dimensional real array, out of place, you might use the following:
    !     
    !     double precision in
    !     dimension in(M,N)
    !     double complex out
    !     dimension out(M/2 + 1, N)
    !     integer*8 plan
    !     
    !     call dfftw_plan_dft_r2c_2d(plan,M,N,in,out,FFTW_ESTIMATE)
    !     call dfftw_execute_dft_r2c(plan, in, out)
    !     call dfftw_destroy_plan(plan)
       call dfftw_plan_dft_r2c_2d(planf,LLm,MMm, 
     &                            rvar,cvar,FFTW_ESTIMATE)
       call dfftw_plan_dft_c2r_2d(planb,LLm,MMm, 
     &                            cvar,rvar,FFTW_ESTIMATE)
       end subroutine init_spectral


!=====================================================================
       real(8) function rand_spec(xx)
!=====================================================================

       real(8) xx

!   definition of kinetic energy spectrum
!      rand_spec = 1./((1.+ (xx*10000.)**2)**(1.)) ! mid
       rand_spec = (1.+xx)**(-100.)
!   rand_spec = 1./((1.+ (xx*100.)**4)**(1.))   ! petite
!   rand_spec = 1./((1.+ (xx*1.)**(-9.))**(1.)) ! grand  
!   randspec = 1./((1.+ (xx/4.)**2)**(4))

       end function rand_spec

