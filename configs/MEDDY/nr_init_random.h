#ifdef GO_MEDDY

!....................................................................
!============================= 
!  	random perturbations 
!=============================   

! initialisation 
      rhopr_rand = 0.  


      if (mynode.eq.0) write(*,*) '--------------------------------------------'
      if (mynode.eq.0) write(*,*) 'random seeds in vertical modes/Fourier space'
      if (mynode.eq.0) write(*,*) '--------------------------------------------'

      call init_spectral(planf,planb,rvar,cvar)    
      call matrice(xlmb,w2)

!     xlambda=sqrt(Nfreq2)*(CELL_DZ*(KM-KL))/(bgf_o(1,1)*sqrt(w2) )
!     xlambda=sqrt(Nfreq2)*(CELL_DZ*(KM-KL))/(bgf_o(1,1)*pi )
!     xkrad1=(CELL_DX*(LLm-0))/(xlambda*2*pi)

       do j=1,32	! (64+1)/2
	  do  i=1,32  	! (64+1)/2
                mk = float(i)/float(64)
                nk = float(j)/float(64)
		kk(i,j) = mk**2 + nk**2
		kk(i,(64+1)-j) = mk**2 + nk**2
		kk((64+1)-i,j) = mk**2 + nk**2
                kk((64+1)-i,(64+1)-j) = mk**2 + nk**2
          enddo
       enddo  

!      valmax = nr_maxofall(maxval(abs(kk)))
!      if (mynode.eq.0) print*,'XXXXXX',valmax,maxval(xlmb),1/w2*xkrad1**2,w2

       etot=0.0
       allocate(prd_spec(1:64,1:64,nmdef))
       
      do mde=1, nmdef
        ekin(mde)=0.0
        do j=1,64
          do  i=1,64
!            xk=sqrt(kk(i,j) + xlmb(mde))
             xk=sqrt(kk(i,j))
             aa= sqrt(rand_spec(xk))
             call random_number(xx)
             phase=2.*pi*xx
             prd_spec(i,j,mde)=cmplx( aa*cos(phase) , aa*sin(phase))
             alp = real( prd_spec(i,j,mde)*conjg(prd_spec(i,j,mde)) )
             ekin(mde)=ekin(mde) + kk(i,j)*alp
          enddo
        enddo
        etot = etot + ekin(mde) 		! kinetic energy only
      enddo


!...... renormalization
      energy   = energy / ((2.*pi)/((LLm+1)/pm(1,1)))**2	! energy in input file  in MKSA
      amp      = sqrt(energy/etot)
      prd_spec = amp * prd_spec
      

      allocate(prd(1:LLm,1:MMm,nmdef))
      prd=0.
      nj = mynode / NP_XI
      ni = mynode - nj*NP_XI   
      
!...... perturbation ds l'espace spectral aux dimensions de la simu
      do mde   = 1,nmdef
        cvar = 0.
        cvar(1:64/2+1,1:32)  = 
     &                    prd_spec(1:64/2+1,1:32,mde)  
        cvar(1:64/2+1,MMm-32+1:MMm)  = 
     &                    prd_spec(1:64/2+1,33:64,mde)   
        
!...... perturbation ds l'espace physique
        call dfftw_execute_dft_c2r(planb,cvar,rvar) 
        prd(:,:,mde) = rvar/(64**2)       
      enddo
      deallocate(prd_spec) 

         
!...... perturbation ds chaque proc
      do k = 1, N                        !    
        do j = 1, Mm
          do i = 1, Lm   
             iofset = i+ni*Lm
             jofset = j+nj*Mm
!            x = xp(i,j)
!            y = yp(i,j) 
!            d2 = ((x-x1)/(2.5*thick_y))**2 +((y-y1)/(2.5*thick_y))**2 
             rhopr_rand(i,j,k) =  prd(iofset,jofset,nmdef)      !*exp(-d2**2)             
          end do
        end do
      end do
      deallocate(prd)

      call exchange_r3d_tile(Istr,Iend,Jstr,Jend,
     &                          rhopr_rand(START_2D_ARRAY,1))


      if (mynode == 0) print *,'max rhopr_rand ', maxval(abs(rhopr_rand))
     
      if (random_pert_urot) then	! coding velocity streamfunction

      if (mynode == 0) then
        write(STDOUT,*) '--------------------------------------------------'
        write(STDOUT,*) 
     & 'random perturbations added to streamfunction field amp=',factor_amp_urot
        write(STDOUT,*) '--------------------------------------------------'
      endif


      endif !random_pert_urot

! end random perturbations
!..........................................................................

#endif /* GO_MEDDY */

