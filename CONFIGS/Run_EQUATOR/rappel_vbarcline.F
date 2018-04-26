#ifdef RAPPEL_VBARCLINE
      if (iic.eq.ntstart) then
        zcent_tourb=-1250.
c        zcent_tourb=-850.
        zdepth_tourb=1000.

c        zcent_tourb=-1000.
c        zdepth_tourb=1500.


        amp_tourb=0.50
        
        do k=1,N
           zshape_tourb(k)= amp_tourb*
c     &         (((z_r(1,1,k)-zcent_tourb)/zdepth_tourb)**2 -1) *
     &          exp(-0.5*( (z_r(1,1,k)-zcent_tourb)/zdepth_tourb)**2 )
        enddo
        vbaramp=0.
        do k=1,N
        vbaramp= vbaramp +
     &            zshape_tourb(k)*Hz(1,1,k)		! prise en compte grille verticale variable
        enddo
        vbaramp=vbaramp/h(0,0)

        MPI_master_only print *, '   '
        MPI_master_only print*,'======================================'
        MPI_master_only print *, 'zcent_tourb ', zcent_tourb 
        MPI_master_only print *, 'zdepth_tourb ', zdepth_tourb
#if  !defined RAPPEL_VBAR
        MPI_master_only print *, 'no barotropic forcing:'
        MPI_master_only print *, ' vbaramp removed '
        zshape_tourb(:) = zshape_tourb(:) - vbaramp
#endif
        MPI_master_only print*,'======================================'
        MPI_master_only print *, '   '


#if  !defined BV_VAR
c projection sur cosinus
        amp_clin1=0
        amp_clin2=0
        h0=h(1,1)
        do k=1,N
        amp_clin1=amp_clin1+ zshape_tourb(k)*cos(pi*z_r(1,1,k)/h0)
        amp_clin2=amp_clin2+ zshape_tourb(k)*cos(2*pi*z_r(1,1,k)/h0)
        enddo
        amp_clin1=2*amp_clin1/N
        amp_clin2=2*amp_clin2/N
        MPI_master_only print *, 'amp_clin 1 & 2', amp_clin1, amp_clin2
# endif /* BV_VAR */
      endif             !iic=nstart

! rappel à tous les pas de temps baroclines
c      if (iic.gt.0) then
      if (.false.) then
      print *,'iic=',iic
      print *,'JV_RANGE', JV_RANGE
      print *,'IV_RANGE', IV_RANGE
      print *,'zshape_tourb', zshape_tourb
      endif
      do k=1,N
        do j=JV_RANGE
          do i=IV_RANGE
            v(i,j,k,nnew)=v(i,j,k,nnew)
     &            +dt*(vrap(i,j)*sin(omegaY*time)*zshape_tourb(k)   
     &                            -v(i,j,k,nnew))
     &           * raptab(i,j)

          enddo
        enddo
        enddo
c	 endif

        
# endif             /* defined RAPPEL_VBARCLINE */

