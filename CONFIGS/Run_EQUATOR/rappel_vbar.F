

#if defined RAPPEL_VBAR ||  defined RAPPEL_VBARCLINE

      if ((iic.eq.ntstart).and.(FIRST_2D_STEP)) then
        rappel_rate = 1./(2* 86400)
        xrmid=0.
        rhalfw=4.5
        xrdecay=1.5
        rhalfn=4.5
        ypdecay=3.
#if !defined RAPPEL_VBARCLINE   
c cas forcage barotrope seul
	vbaramp = 0.20				!  amplitude in v (m/sec)
#else       /* RAPPEL_VBARCLINE */
	vbaramp = 0.2241588296622976		!!!!!!!!!!!! pas propre du tout!!!
#endif        /* RAPPEL_VBARCLINE */


        MPI_master_only print *, '   '
        MPI_master_only print*,'======================================'
        MPI_master_only print *, 'FORCAGE OUEST'
        MPI_master_only print *, 'vbaramp   ', vbaramp
        MPI_master_only print *,
     &                   'restoring rate in v only(days) RAPPEL_V',
     &                        1./(86400*rappel_rate)
        MPI_master_only print *,'localized plateau middle ',xrmid
        MPI_master_only print *,'half-width=',rhalfw,' decay',xrdecay
        MPI_master_only print*,'======================================'
        MPI_master_only print *, '   '

         mz = pi/h(0,0)      	                ! 1rst baroclinic mode (m^(-1))
	 onedeg = pi*Eradius/180                ! in meters
         lrossby=sqrt(xnconst/(mz*beta))/onedeg  ! in degrees of latitude
c         print *,' lrossby ', lrossby, 'beta ',beta
c         lx = 2.47                         ! width in degrees by convention; 73j
c         lx = 3.15                          ! 57j
         lx = 3.6			    ! 50j
c         lx = 6.                           ! 30j
         kx= -2*pi/(lx*onedeg )
         omegaY=-beta/kx
c         print *,'lx, kx, onedeg, beta, omegaY'
c         print *,lx, kx, onedeg, beta, omegaY

c 7jan08: ltourb=rayon du tourbillon de bord Ouest; ytourb=latitude du centre du tourbillon
        ltourb=3.
c        ytourb=2.
        ytourb=0.			! mieux vaut pas le decentrer
! dependance en y du champ vers lequel on rappelle
        do j=JstrV,Jend
          do i=Istr,Iend 
                yshape_tourb(i,j)=
     &               exp(-0.5*((yp(1,j)- ytourb)/ltourb)**2.)	
              enddo
           enddo
      vrap=0.
      raptab=0.
! champ vers lequel on rappelle	   
        do j=JstrV,Jend
          do i=Istr,Iend 
                vrap(i,j)= yshape_tourb(i,j)
     &              *cos(kx*xr(i,1)*onedeg)		
              enddo
           enddo
! coefficient de rappel
        do j=JstrV,Jend
          do i=Istr,Iend
            raptab(i,j)=rappel_rate*
     &              0.5*(tanh((xr(i,1)-xrmid+rhalfw)/xrdecay)
     &                   -tanh((xr(i,1)-xrmid-rhalfw)/xrdecay))  
     &             *0.5*(tanh((yp(1,j)-ytourb+rhalfn)/ypdecay)
     &                   -tanh((yp(1,j)-ytourb-rhalfn)/ypdecay))
          enddo
        enddo
c	 print *,'inside bar JstrV,Jend',JstrV,Jend
c	 print *,'           Istr,Iend',Istr,Iend

        MPI_master_only print *, '   '
        MPI_master_only print*,'==================================='
        MPI_master_only print *,' ytourb:',ytourb
        MPI_master_only print *,' ltourb:',ltourb
	 MPI_master_only print *,
     &   '  zonal wavelength (degrees) : ',(2*pi/kx)/onedeg
        MPI_master_only print *,'  period (days)',
     &           2.*pi/omegaY/86400
	 MPI_master_only print *,'group propagation speed(m/sec) ',
     &           beta/kx**2
        MPI_master_only print*,'==================================='
        MPI_master_only print *, '   '
      endif		!iic=nstart

#ifdef RAPPEL_VBAR
! rappel a tous les pas de temps barotropes
        do j=JstrV,Jend
          do i=Istr,Iend
            vbar(i,j,knew)=vbar(i,j,knew)
     &            +dtfast*(vrap(i,j)*sin(omegaY*time)*vbaramp
     &                            -vbar(i,j,knew))
     &           * raptab(i,j)   
          enddo
        enddo
# endif             /* defined RAPPEL_VBAR */
        
# endif             /* defined RAPPEL_VBAR || RAPPEL_VBARCLINE */
