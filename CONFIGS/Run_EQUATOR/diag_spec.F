
#include "cppdefs.h"
  
      module diagspec
   
#ifdef DIAG_SPEC

      implicit none
      
      integer :: LLmspec,MMmspec
      integer :: Nax
      integer, parameter :: nspec=18  ! nb of cospectra
      integer, parameter :: nhat=31   ! nb of variables to Fourierize
      integer :: navgspec, nwrtspec, nrecspec
      
      integer :: nint, ncidspec, ktot
      integer*8,dimension(:,:),allocatable :: lookup_index
      integer*8 :: plan
      integer*8,dimension(:),allocatable ::reqid
      integer ,dimension(:), allocatable :: index_count
      
      real(kind=8),dimension(:,:),allocatable :: zglo
      real(kind=8),dimension(:,:),allocatable :: zpar
      real(kind=8),dimension(:,:,:),allocatable :: zexchang
      real(kind=8),dimension(:),allocatable::kx
      real(kind=8),dimension(:),allocatable :: ky
      real(kind=8),dimension(:,:),allocatable :: kk
      real(kind=8),dimension(:),allocatable :: rho_r,drhordz  
      real(kind=8),dimension(:,:),allocatable :: infft 
         
      character*24,dimension(:),allocatable:: name_cospec
      character *80 :: specname
      
      complex(kind=8),dimension(:,:,:),allocatable :: cospec
      complex(kind=8),dimension(:,:,:),allocatable :: hat
      complex(kind=8),dimension(:,:),allocatable :: outfft      


      contains


  !*************************************************************************
 
       subroutine init_diags_spec
      

! Allocation et initialisation des tableaux
! Calcul des indices pour cumuler sur des bandes spectrales (lookup_index)
! Initialisation des FFT
! Calcul d'un rho de référence (rho_r)

# include "param.h"
# include "scalars.h"
# include "ocean3d.h"
# include "grid.h"

      include 'mpif.h'
      include 'fftw3.f'

      !local
      real(kind=8)::kh,coef,ni,nj
      integer::i,j,k,l,kint,iproc,iret,i0,j0,mpierr
      integer,dimension(MPI_STATUS_SIZE):: istat
      real(kind=8) :: localsum  , globalsum, nr_sumofall 

      LLmspec=(LLm)/2+1
      MMmspec=MMm
      if ( mod(N,NNODES) == 0 ) then
         nint=N/NNODES
      else
         nint=N/NNODES + 1
      endif

! allocation dynamique des tableaux

      allocate(lookup_index(LLmspec,MMmspec))
      allocate(reqid(0:NNODES-1))
      allocate(zglo(LLm,MMm))
      allocate(zpar(Lm,Mm))
      allocate(zexchang(Lm,Mm,0:NNODES-1))
      allocate(kx(LLmspec))
      allocate(ky(MMmspec))
      allocate(kk(LLmspec,MMmspec))
      allocate(hat(LLmspec,MMmspec,nhat)) 
      allocate(rho_r(N))
      allocate(drhordz(N)) 
!     allocate(infft(LLm-1,MMm-1))
      allocate(infft(LLm,MMm))
      allocate(outfft(LLmspec,MMmspec))

! Initialisation des noms de variables pour le fichier netcdf

      allocate(name_cospec(nspec))       
      name_cospec(1)="VdVdt "
      name_cospec(2)="hor_adv"
      name_cospec(3)="ver_adv"
      name_cospec(4)="hor_advC"
      name_cospec(5)="cor"     
      name_cospec(6)="hor_pres"
      name_cospec(7)="ver_pres"
      name_cospec(8)="hor_vis"
      name_cospec(9)="ver_vis"
      name_cospec(10)="ke"
      name_cospec(11)="rhodrhodt"
      name_cospec(12)="rho_hadv"
      name_cospec(13)="rho_vadv"  
      name_cospec(14)="rho_hadvC"          
      name_cospec(15)="rho_hdiff"
      name_cospec(16)="rho_vdiff"       
      name_cospec(17)="rho_forc"
      name_cospec(18)="ape"          
            

!     Make wave number mask

      call get_kxky          
      
      do j = 1,MMmspec
         do i = 1,LLmspec
            lookup_index(i,j)=floor(kk(i,j))+1
         enddo
      enddo 
      Nax = maxval(lookup_index)  
      
      allocate(index_count(Nax))
      allocate(cospec(Nax,nint,nspec))
      cospec=0.  
      nrecspec=0

! output paramètres diags spectraux
 
      if(mynode.eq.0)then
         write(*,*)
         write(*,*)"=========== Initialisation diags spectraux ================"
         write(*,*)
         write(*,*)" Nombre de niveaux traités par un processeur,nint  = ",nint  
         write(*,*)" Nombre d ondes, Nax =",Nax  
         write(*,*)" Nombre de variables fourierisées, nhat  = ",nhat           
         write(*,*)" Nombre de co-spectres calculés, nspec = ",nspec
         write(*,*)" Fichier netcdf = diag_spec.nc"
         write(*,*)" Co-spectres calculés = "   
         do l=1,nspec
            write(*,*)"           ",trim(name_cospec(l))          
         enddo      
         write(*,*)
      endif
           
!
!To transform a two-dimensional real array, out of place, you might use the following:
!
!             double precision in
!             dimension in(M,N)
!             double complex out
!             dimension out(M/2 + 1, N)
!             integer*8 plan
!     
!             call dfftw_plan_dft_r2c_2d(plan,M,N,in,out,FFTW_ESTIMATE)
!             call dfftw_execute_dft_r2c(plan, in, out)
!             call dfftw_destroy_plan(plan)

!     call dfftw_plan_dft_r2c_2d(plan,LLm-1,MMm-1,infft,outfft,
!    &                           FFTW_ESTIMATE)   
      call dfftw_plan_dft_r2c_2d(plan,LLm,MMm,infft,outfft,
     &                           FFTW_ESTIMATE)       

      
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

      end subroutine init_diags_spec
      
!========================================

      subroutine get_kxky
      
# include "param.h"

      real(kind=8):: pi, coefx, coefy
      integer::i,j,k

      pi=acos(-1.)
      kx=0.
      ky=0.
      kk=0.
      
      kx = (/ (i, i= 0,LLmspec-1) /)               
      ky = (/ (j, j=1,MMmspec/2),(j, j=-MMmspec/2,0) /) 
      do j = 1,MMmspec
         do i = 1,LLmspec
            kk(i,j)=sqrt(kx(i)*kx(i) + ky(j)*ky(j))
         enddo
      enddo
      
      end subroutine get_kxky

!========================================

      subroutine comp_diags_spec
      
       
! Calcul des co-spectres -> cospec(kx,k,:)
! Variables spectrales:
!	1 = u(n+1/2)
!	2 = v(n+1/2)
!	3 = horizontal advection of u
!	4 = vertical advection of u
!	5 = horizontal advection of u C4 scheme
!	6 = horizontal advection of v
!	7 = vertical advection of v
!	8 = horizontal advection of v C4 scheme
!	9 = Coriolis part of u
!	10 = Coriolis part of v
!	11 = Pressure gradient of u
!	12 = Pressure gradient of v
!	13 = Pressure gradient of w
!	14 = Horizontal mixing of u
!	15 = vertical mixing of u
!	16 = horizontal mixing of v 
!	17 = vertical mixing of v
!	18 = dudt
!	19 = dvdt
!	20 = (u(n+1)+u(n))/2
!	21 = (v(n+1)+v(n))/2
!	22 = w
!	23 = rho(n+1/2)
!	24 = horizontal advection of rho
!	25 = vertical advection of rho
!	26 = horizontal advection of rho C4 or C6 scheme
!	27 = horizontal mixing
!	28 = vertical mixing
!	29 = forcing of rho
!	30 = drhodt
!	31 = (rho(n+1)+rho(n)/2
!
! co-spectres:
!	1 = ududt + vdvdt
!	2 = horizontal advection
!	3 = vertical advection
!	4 = horizontal advection C4 scheme
!	5 = Coriolis
!	6 = Horizontal pressure
!	7 = Vertical pressure
!	8 = horizontal viscosity
!	9 = vertical viscosity
!       10 = ke
!	11 = rhodrhodt
!	12 = horizontal advection of rho 	
!	13 = vertical advection of rho
!	14 = horizontal advection of rho C4 or C6 scheme	
!	15 = horizontal diffusion of rho
!	16 = vertical diffusion of rho
!	17 = Forcing of rho
!       18 = ape qg

# include "param.h"
# include "scalars.h"

	include 'mpif.h'

 
      !local
      integer :: ik,k,kint,kv


      ! échanges aux bords des sous-domaines
      call exchange_diags_spec

      do kint=1,nint  
      
         k=(mynode*nint)+kint
         
         do kv=1,nhat
            call comp_hat(kv,kint) 
         enddo

         if ((mynode*nint + kint) <= N) then 
                                 
! momentum equation, total kinetic energy, u*dudt  * v*dvdt  
            call comp_cospec(20,18,1,kint) 
            call comp_cospec(21,19,1,kint)                                
   
! momentum equation, horizontal advection : u*(u*gradu) + v*(v*gradv)
            call comp_cospec(1,3,2,kint) 
            call comp_cospec(2,6,2,kint)

! momentum equation, vertical advection : u*(u*gradu) + v*(v*gradv) 
            call comp_cospec(1,4,3,kint) 
            call comp_cospec(2,7,3,kint) 

! momentum equation, horizontal advection C4 scheme : u*(u*gradu) + v*(v*gradv)
            call comp_cospec(1,5,4,kint) 
            call comp_cospec(2,8,4,kint)  
            
! momentum equation, coriolis : u*(fv) - v*(fu)
            call comp_cospec(1,9,5,kint) 
            call comp_cospec(2,10,5,kint) 

! Momentum equation, gradient de pression horizontal : -u*dpdx - v*dpdy                        
            call comp_cospec(1,11,6,kint)  
            call comp_cospec(2,12,6,kint) 

! Momentum equation, gradient de pression vertical : -w*dpdz                  
            call comp_cospec(22,13,7,kint) 

! momentum equation, horizontal viscosity         
            call comp_cospec(1,14,8,kint)  
            call comp_cospec(2,16,8,kint) 
                          
! Momentum equation, vertical viscosity 
            call comp_cospec(1,15,9,kint) 
            call comp_cospec(2,17,9,kint)                        
              
! Momentum equation, ke, (u*u + v*v)
            call comp_cospec(1,1,10,kint) 
            call comp_cospec(2,2,10,kint) 

! density,-g/2/drhordz(k,nt)/rho0* rho*drhodt   
            call comp_cospec(31,30,11,kint)                            
                                 
! density horizontal advection, rho*(v*grad.rho)         
            call comp_cospec(23,24,12,kint)                       

! density vertical advection, rho*(v*grad.rho)          
            call comp_cospec(23,25,13,kint) 
            
! density horizontal advection C4 or C6 scheme, rho*(v*grad.rho)         
            call comp_cospec(23,26,14,kint)                                
                  
! density, horizontal diffusion, rho*(diffusion horizontale)
            call comp_cospec(23,27,15,kint) 
                
! density, vertical diffusion, rho*(diffusion horizontale)
            call comp_cospec(23,28,16,kint)                 

! density, forcing, rho*(forcing)
            call comp_cospec(23,29,17,kint)                   

! density, ape QG (rho*rho)
            call comp_cospec(23,23,18,kint) 

         endif

      enddo
      
      end subroutine  comp_diags_spec
     
!========================================

      subroutine comp_cospec(kv1,kv2,ksp,kint)

# include "param.h"
# include "grid.h"
     
      integer :: kv1,kv2,ksp,kint      
      


      ! compute the ksp-th cospectrum (azimuthal averaging)
      ! between hat(kv1) and hat(kv2) for vertical segment kint      

           

	include 'mpif.h'
        include 'fftw3.f'      
      
      !local
      integer ::i,j,ik
      complex(kind=8)::cov
      
      index_count=0
      do j = 1,MMmspec     
         do i = 1,LLmspec
            ik = lookup_index(i,j)
            index_count(ik) = index_count(ik)+1
            if (i==1) then
              cov=conjg(hat(i,j,kv1))*hat(i,j,kv2)
            else
              cov=conjg(hat(i,j,kv1))*hat(i,j,kv2) * 2. !  Symetrie espace spectral
            endif
            cospec(ik,kint,ksp) = cospec(ik,kint,ksp) + cov
         enddo
      enddo  
      
      end subroutine comp_cospec
      
!========================================

      subroutine comp_hat(kv,kint)
 
# include "param.h" 
# include "scalars.h"
     
      integer::kv,kint   
      real :: rho_mean 
        
      
      ! compute complex amplitude for variable kv and level kint

	include 'mpif.h'
        include 'fftw3.f'

      !local
      integer::iproc,i,j,k
      real :: mean

      k=(mynode*nint)+kint
      if (k <= N) then
         call launch_rcv
      endif
      do iproc=0,NNODES-1      
         ktot=(iproc*nint)+kint
         if (ktot <= N) then
            call compute_var(kv,ktot)
            call launch_send(iproc)
         endif
      enddo
      if (k <= N) then
         call launch_wait   
      endif       
      
             
      if (k <= N) then    
       
! fourierize
         infft = zglo  
         call dfftw_execute_dft_r2c(plan,infft,outfft)
         hat(:,:,kv)=outfft/(LLm)/(MMm)

      endif



      end subroutine comp_hat


!========================================
      subroutine compute_var(var,k)
     
      integer::var,k            
        
# include "param.h"
# include "scalars.h"
# include "ocean3d.h"
# include "diagnostics.h"

      include 'mpif.h'


      !local
      integer::i,j,tile
      real :: rhodiag, mean
      real, dimension(GLOBAL_2D_ARRAY,1:N) :: workr
         

      select case(var)

!---------------------------------------------------------------
! Momentum diagnostic variables.
!---------------------------------------------------------------

      case(1)
!     compute u at level k
         do j=1,Mm
            do i=1,Lm
               zpar(i,j)=MV(i,j,k,1)
            enddo
         enddo
         
      case(2)
!     compute v at level k
         do j=1,Mm
            do i=1,Lm
               zpar(i,j)=MV(i,j,k,2)
            enddo
         enddo
         
      case(3)
!     compute horizontal advection of u at level kl
         do j=1,Mm
            do i=1,Lm
               zpar(i,j)=MXadv(i,j,k,1)+MYadv(i,j,k,1)
            enddo
         enddo
         
      case(4)
!     compute vertical advection of u at level kl
         do j=1,Mm
            do i=1,Lm
               zpar(i,j)=MVadv(i,j,k,1)
            enddo
         enddo
 
      case(5)
!     compute horizontal advection of u C4 scheme at level kl
         do j=1,Mm
            do i=1,Lm
               zpar(i,j)=MXadvC(i,j,k,1)+MYadvC(i,j,k,1)
            enddo
         enddo 
                
      case(6)
!     compute horizontal advection of v at level kl
         do j=1,Mm
            do i=1,Lm
               zpar(i,j)=MXadv(i,j,k,2)+MYadv(i,j,k,2)
            enddo
         enddo   
      
      case(7)
!     compute vertical advection of v at level k
         do j=1,Mm
            do i=1,Lm
               zpar(i,j)=MVadv(i,j,k,2)
            enddo
         enddo
         
      case(8)
!     compute horizontal advection of v C4 scheme at level kl
         do j=1,Mm
            do i=1,Lm
               zpar(i,j)=MXadvC(i,j,k,2)+MYadvC(i,j,k,2)
            enddo
         enddo   
         
      case(9)
!     compute Coriolis part of u at level k
         do j=1,Mm
            do i=1,Lm
               zpar(i,j)=MCor(i,j,k,1)
            enddo
         enddo
         
      case(10)
!     compute coriolis part of v at level k
         do j=1,Mm
            do i=1,Lm
               zpar(i,j)=MCor(i,j,k,2)
            enddo
         enddo  
                       
      case(11)
!     compute Pressure gradient of u at level k
         do j=1,Mm
            do i=1,Lm
               zpar(i,j)=MPrsgrd(i,j,k,1)
            enddo
         enddo
         
         
      case(12)
!     compute Pressure gradient of v at level k
         do j=1,Mm
            do i=1,Lm
               zpar(i,j)=MPrsgrd(i,j,k,2)
            enddo
         enddo 
         
      case(13)        
!     compute Pressure gradient of w at level k
         do j=1,Mm
            do i=1,Lm
               zpar(i,j)= g / rho0 * rho(i,j,k)
            enddo
         enddo  
                    
      case(14)
!     compute  horizontal mixing of u at level k

         do j=1,Mm
            do i=1,Lm
               zpar(i,j)=MHmix(i,j,k,1)
            enddo
         enddo
         
      case(15)
!     compute  vertical mixing of u at level k

         do j=1,Mm
            do i=1,Lm
               zpar(i,j)=MVmix(i,j,k,1)
            enddo
         enddo
                 
      case(16)
!     compute  horizontal mixing of v at level k
         do j=1,Mm
            do i=1,Lm
               zpar(i,j)=MHmix(i,j,k,2)
            enddo
         enddo 
           
      case(17)
!     compute  vertical mixing of v at level k
         do j=1,Mm
            do i=1,Lm
               zpar(i,j)=MVmix(i,j,k,2)
            enddo
         enddo           
 
      case(18)
!     compute dudt at level k      
         do j=1,Mm
            do i=1,Lm
               zpar(i,j)=Mrate(i,j,k,1)
            enddo
         enddo             
                         
      case(19)
!     compute dvdt at level k
         do j=1,Mm
            do i=1,Lm
               zpar(i,j)=Mrate(i,j,k,2)
            enddo
         enddo

      case(20)
!     compute (u(n+1)+u(n))/2 at level k
         do j=1,Mm
            do i=1,Lm
               zpar(i,j)=MVmoy(i,j,k,1)
            enddo
         enddo   
         
      case(21)
!     compute (v(n+1)+v(n))/2 at level k
         do j=1,Mm
            do i=1,Lm
               zpar(i,j)=MVmoy(i,j,k,2)
            enddo
         enddo 
                         
      case(22)
!     compute w at level k
         do tile=0,NSUB_X*NSUB_E-1
            call Wvlcty (tile, workr)
         enddo
         do j=1,Mm
            do i=1,Lm
               zpar(i,j)=workr(i,j,k)
            enddo
         enddo
         
!---------------------------------------------------------------
! Tracer diagnostic variables.
!---------------------------------------------------------------
         
      case(23)         
!     compute rho at level k
         do j=1,Mm
            do i=1,Lm
               rhodiag= R0 - Tcoef*(TT(i,j,k,itemp)-T0)
#  ifdef SALINITY
     &                     + Scoef*(TT(i,j,k,isalt)-S0)
#  endif     
               zpar(i,j)=rhodiag - rho_r(k)  
            enddo
         enddo   
            
      case(24)
!     compute horizontal advection of rho at level k
         do j=1,Mm
            do i=1,Lm
               zpar(i,j)= 					  				
     &                 - (TXadv(i,j,k,itemp) +			 
     &                    TYadv(i,j,k,itemp) )*Tcoef		 
# ifdef SALINITY                      
     &                 + (TXadv(i,j,k,isalt) +			 
     &                    TYadv(i,j,k,isalt) )*Scoef 		 
# endif
            enddo
         enddo

      case(25)
!     compute vertical advection of rho at level k
         do j=1,Mm
            do i=1,Lm
               zpar(i,j)= 					  				
     &                 - TVadv(i,j,k,itemp)*Tcoef		 
# ifdef SALINITY                      
     &                 + TVadv(i,j,k,isalt)*Scoef 		 
# endif
            enddo
         enddo
 
      case(26)
!     compute horizontal advection of t & s at level k C4 or C6 scheme
         do j=1,Mm
            do i=1,Lm
               zpar(i,j)= 					  				
     &                 - (TXadvC(i,j,k,itemp) +			 
     &                    TYadvC(i,j,k,itemp) )*Tcoef		 
# ifdef SALINITY                      
     &                 + (TXadvC(i,j,k,isalt) +			 
     &                    TYadvC(i,j,k,isalt) )*Scoef 		 
# endif
            enddo
         enddo
                      
         
      case(27)
!     compute horizontal mixing of rho at level k
         do j=1,Mm
            do i=1,Lm
               zpar(i,j)= 					  				
     &                 - THmix(i,j,k,itemp)*Tcoef		 
# ifdef SALINITY                      
     &                 + THmix(i,j,k,isalt)*Scoef 		 
# endif
            enddo
         enddo
         
      case(28)
!     compute vertical mixing of rho at level k
         do j=1,Mm
            do i=1,Lm
               zpar(i,j)= 					  				
     &                 - TVmix(i,j,k,itemp)*Tcoef		 
# ifdef SALINITY                      
     &                 + TVmix(i,j,k,isalt)*Scoef 		  
# endif
            enddo
         enddo

      case(29)
!     compute forcing of rho at level k
         do j=1,Mm
            do i=1,Lm
               zpar(i,j)= 					  				
     &                 - TForc(i,j,k,itemp)*Tcoef		 
# ifdef SALINITY                      
     &                 + TForc(i,j,k,isalt)*Scoef 		 
# endif
            enddo
         enddo
   
      case(30)
!     compute drhodt at level k
         do j=1,Mm
            do i=1,Lm
               zpar(i,j)= 					  				
     &                 - Trate(i,j,k,itemp)*Tcoef		 
# ifdef SALINITY                      
     &                 + Trate(i,j,k,isalt)*Scoef 		 
# endif
            enddo
         enddo
 
      case(31)         
!     compute (rho(N)+rho(N+1))/2 at level k
         do j=1,Mm
            do i=1,Lm
               rhodiag = R0 - Tcoef*(TTmoy(i,j,k,itemp)-T0)
#  ifdef SALINITY
     &                     + Scoef*(TTmoy(i,j,k,isalt)-S0)
#  endif     
               zpar(i,j) = rhodiag - rho_r(k)  
            enddo
         enddo                                    
      end select
      
      end subroutine compute_var

!========================================

      subroutine launch_rcv
                 
# include "param.h"
# include "scalars.h"
# include "mpi_roms.h"

	include 'mpif.h'

	

      integer::iproc,i0,j0,npts,mpierr
      
!     launch all the nonblocking receive requests
      npts = Lm*Mm
      zexchang = 0.
      do iproc=0,NNODES-1
         if(mynode.ne.iproc)then
            call MPI_IRECV(zexchang(:,:,iproc),npts,  
     &            MPI_DOUBLE_PRECISION,iproc,iproc,MPI_COMM_WORLD,
     &            reqid(iproc),mpierr)
         endif
      enddo

      
      end subroutine launch_rcv

!========================================

      subroutine launch_wait
                 
# include "param.h"
# include "scalars.h"

	include 'mpif.h'


      integer::iproc,i0,j0,i,j,mpierr,i_offset,j_offset
      integer,dimension(MPI_STATUS_SIZE):: istat
      
!     launch all the wait
      do iproc=0,NNODES-1
         if(mynode.ne.iproc)then
            call MPI_WAIT(reqid(iproc),istat,mpierr)
         endif
      enddo
      do iproc=0,NNODES-1 
         if(mynode.ne.iproc)then         
           do j=1,Mm
             do i=1,Lm
               i_offset=mod(iproc,NP_XI)*Lm
               j_offset=(iproc/NP_XI)*Mm 
               i0 = i+i_offset
               j0 = j+j_offset
               zglo(i0,j0) = zexchang(i,j,iproc)    
             enddo
           enddo 
         endif           
      enddo
      
      end subroutine launch_wait

!========================================

      subroutine launch_send(iproc)

      integer::iproc
      
# include "param.h"
# include "scalars.h"
# include "mpi_roms.h"

	include 'mpif.h'


      !local
      integer::i,j,i0,j0,npts,mpierr,i_offset,j_offset
      
      if(mynode.eq.iproc)then
!copy zpar into the receive buffer
         do j=1,Mm
            do i=1,Lm
               i_offset=mod(iproc,NP_XI)*Lm
               j_offset=(iproc/NP_XI)*Mm 
               i0 = i + i_offset
               j0 = j + j_offset
               zglo(i0,j0)=zpar(i,j)
            enddo
         enddo            
      else
!send to iproc
         npts = Lm*Mm
         call MPI_SEND(zpar,npts,MPI_DOUBLE_PRECISION,iproc,mynode,  
     &                  MPI_COMM_WORLD,mpierr)
      endif


      end subroutine launch_send

!======================================== 

      subroutine wrt_diags_spec(iflag)
      
!     use netcdf      
      
      integer:: iflag           

# include "param.h"
# include "scalars.h"
# include "ocean3d.h"
# include "grid.h"
# include "netcdf.inc"
# include "mpi_roms.h"

      include 'mpif.h'


      !local
      integer::ierr
      integer::iproc,i,j,k,kint,l,Nin,kcount,mpierr
      integer,save:: nbtime,timeid
      integer,save::dummyid
      integer,dimension(40),save:: dimid,varid
      character*64 :: name
      integer,dimension(MPI_STATUS_SIZE):: istat

      select case(iflag)
      case(0) 		! initialization (open netcdf file)
      if (mynode.eq.0) then

         
! Create the netCDF file. The NF_NETCDF4 flag causes a
! HDF5/netCDF-4 file to be created. The comm and info parameters
! cause parallel I/O to be enabled.

         specname = 'diag_spec.nc' 
         call check(nf_create(trim(specname), NF_CLOBBER, ncidspec)  )   

         call check(nf_def_dim(ncidspec, 'kk',Nax,dimid(1)) )       
         call check(nf_def_dim(ncidspec, 'level',N,dimid(2)) )     
         nbtime = int((ntimes-ntstart)/nwrtspec)+1
         call check(nf_def_dim(ncidspec, 'time',NF_UNLIMITED,dimid(3)) )         
     
! put global attributes

         call check(nf_put_att_int(ncidspec,NF_GLOBAL,'nx',NF_INT,1,(LLm)))
         call check(nf_put_att_int(ncidspec,NF_GLOBAL,'ny',NF_INT,1,(MMm)))        
         call check(nf_put_att_double(ncidspec,NF_GLOBAL,'dx',
     &              NF_DOUBLE,1,1./pm(1,1)))  
         call check(nf_put_att_double(ncidspec,NF_GLOBAL,'dy',
     &              NF_DOUBLE,1,1./pn(1,1)))
         call check(nf_put_att_double(ncidspec,NF_GLOBAL,'dz',
     &              NF_DOUBLE,1,Hz(1,1,1)))

! define variables

         call check(nf_def_var(ncidspec,'time' ,
     &                    nf_double,1,(/dimid(3)/),timeid))
         do l=1,nspec
            name=name_cospec(l)            
            call check(nf_def_var(ncidspec,trim(name),NF_DOUBLE,3,  
     &                       (/dimid(1),dimid(2),dimid(3)/),varid(l)))                        
         enddo        
    
         call check(nf_enddef(ncidspec) )
          
      endif 
              
      case(1) ! write spec

! each processor write it's own cospec in the parallel netcdf file if needed
          
          nrecspec=nrecspec+1
          
          cospec(:,:,10) = 0.5 * cospec(:,:,10) 
          cospec(:,:,18) = 0.5 * cospec(:,:,18) 
                   
          do l=11,nspec             
             do k=1,nint
                ktot=mynode*nint+k               
                if (ktot <= N) then             
                     cospec(:,k,l) = - g/drhordz(ktot)/rho0*cospec(:,k,l)
                endif
             enddo
          enddo

          cospec = cospec / float(nwrtspec/navgspec)      

          if (mynode==0)  call check(nf_put_vara_double(ncidspec,timeid,  
     &                     (/nrecspec/),(/1/),time*sec2day) ) ! time in days      

          ! reads/writes to happen together on every processor.
          do l=1,nspec

             if (mynode==0) then
                   Nin=1
                   kcount=nint 

                   call check(nf_put_vara_double(ncidspec, varid(l),  
     &                  (/1,Nin,nrecspec/),  
     &                  (/Nax,kcount,1/),
     &                  real(cospec(1:Nax,1:kcount,l))) )  
             else
                    Nin=mynode*nint+1
                    if (Nin <= N) then
                       call MPI_SEND(cospec,Nax*nint*nspec,                
     &                       MPI_DOUBLE_COMPLEX,0,mynode,  
     &                       MPI_COMM_WORLD,mpierr)  
                    endif
             endif                     
               
             if (mynode==0) then
                  do iproc=1,NNODES-1
                     Nin=iproc*nint+1           
                     if (Nin <= N) then                    
                        call MPI_IRECV(cospec,  
     &                      Nax*nint*nspec,  
     &                      MPI_DOUBLE_COMPLEX,iproc,iproc,  
     &                      MPI_COMM_WORLD,reqid(iproc),mpierr)
                        call MPI_WAIT(reqid(iproc),istat,mpierr)
                        kcount = min(nint,N-Nin+1)  
 	                                             
                        call check(nf_put_vara_double(ncidspec, varid(l),  
     &                     (/1,Nin,nrecspec/),  
     &                     (/Nax,kcount,1/),
     &                     real(cospec(1:Nax,1:kcount,l)))  )     		   
                     endif
                  enddo
               endif
            enddo
                              
            if (nrecspec==nbtime) then
               ierr = nf_CLOSE(ncidspec)
            else
               ierr = nf_sync(ncidspec)
            endif
            cospec = 0. 
                 
      end select

      
      end subroutine wrt_diags_spec  
      
!======================================== 
       
      subroutine check(status)
      
!     use netcdf
# include "netcdf.inc"   
  
      integer, intent ( in) :: status
      integer :: err
    
      if(status /= nf_noerr) then 
        print *, trim(nf_strerror(status))
        call mpi_finalize(err)
        stop
      endif
      end subroutine check 
       
#endif  /* DIAG_SPEC */  
    
      end module diagspec
       
!======================================== 

      subroutine exchange_diags_spec()    
# include "param.h"
# include "scalars.h"
# include "ocean3d.h"
# include "diagnostics.h"

      include 'mpif.h'


#ifdef DIAG_SPEC
    
      !local
      integer::i,j,tile
      real, dimension(GLOBAL_2D_ARRAY,1:N) :: workr  
      
# include "compute_tile_bounds.h"       
    
      call exchange_u3d_tile (Istr,Iend,Jstr,Jend,
     &                        Mrate(START_2D_ARRAY,1,1))
      call exchange_u3d_tile (Istr,Iend,Jstr,Jend,
     &                        MXadv(START_2D_ARRAY,1,1)) 
      call exchange_u3d_tile (Istr,Iend,Jstr,Jend,
     &                        MYadv(START_2D_ARRAY,1,1))
      call exchange_u3d_tile (Istr,Iend,Jstr,Jend,
     &                        MXadvC(START_2D_ARRAY,1,1)) 
      call exchange_u3d_tile (Istr,Iend,Jstr,Jend,
     &                        MYadvC(START_2D_ARRAY,1,1))    
      call exchange_u3d_tile (Istr,Iend,Jstr,Jend,
     &                        MVadv(START_2D_ARRAY,1,1))
      call exchange_u3d_tile (Istr,Iend,Jstr,Jend,
     &                        MCor(START_2D_ARRAY,1,1))
      call exchange_u3d_tile (Istr,Iend,Jstr,Jend,
     &                        MPrsgrd(START_2D_ARRAY,1,1))   
      call exchange_u3d_tile (Istr,Iend,Jstr,Jend,
     &                        MHmix(START_2D_ARRAY,1,1))
      call exchange_u3d_tile (Istr,Iend,Jstr,Jend,
     &                        MVmix(START_2D_ARRAY,1,1))
      call exchange_u3d_tile (Istr,Iend,Jstr,Jend,
     &                        MV(START_2D_ARRAY,1,1))    
      call exchange_u3d_tile (Istr,Iend,Jstr,Jend,
     &                        MVmoy(START_2D_ARRAY,1,1))  
       
      call exchange_v3d_tile (Istr,Iend,Jstr,Jend,
     &                        Mrate(START_2D_ARRAY,1,2))
      call exchange_v3d_tile (Istr,Iend,Jstr,Jend,
     &                        MXadv(START_2D_ARRAY,1,2)) 
      call exchange_v3d_tile (Istr,Iend,Jstr,Jend,
     &                        MYadv(START_2D_ARRAY,1,2))
      call exchange_v3d_tile (Istr,Iend,Jstr,Jend,
     &                        MXadvC(START_2D_ARRAY,1,2)) 
      call exchange_v3d_tile (Istr,Iend,Jstr,Jend,
     &                        MYadvC(START_2D_ARRAY,1,2))    
      call exchange_v3d_tile (Istr,Iend,Jstr,Jend,
     &                        MVadv(START_2D_ARRAY,1,2))
      call exchange_v3d_tile (Istr,Iend,Jstr,Jend,
     &                        MCor(START_2D_ARRAY,1,2))
      call exchange_v3d_tile (Istr,Iend,Jstr,Jend,
     &                        MPrsgrd(START_2D_ARRAY,1,2))   
      call exchange_v3d_tile (Istr,Iend,Jstr,Jend,
     &                        MHmix(START_2D_ARRAY,1,2))
      call exchange_v3d_tile (Istr,Iend,Jstr,Jend,
     &                        MVmix(START_2D_ARRAY,1,2))
      call exchange_v3d_tile (Istr,Iend,Jstr,Jend,
     &                        MV(START_2D_ARRAY,1,2))    
      call exchange_v3d_tile (Istr,Iend,Jstr,Jend,
     &                        MVmoy(START_2D_ARRAY,1,2))  
       
      do itrc=1,NT
        call exchange_r3d_tile (Istr,Iend,Jstr,Jend,
     &                          TXadv(START_2D_ARRAY,1,itrc))
        call exchange_r3d_tile (Istr,Iend,Jstr,Jend,
     &                          TYadv(START_2D_ARRAY,1,itrc))
        call exchange_r3d_tile (Istr,Iend,Jstr,Jend,
     &                          TVadv(START_2D_ARRAY,1,itrc))
        call exchange_r3d_tile (Istr,Iend,Jstr,Jend,
     &                          TVmix(START_2D_ARRAY,1,itrc))
        call exchange_r3d_tile (Istr,Iend,Jstr,Jend,
     &                          TForc(START_2D_ARRAY,1,itrc))
        call exchange_r3d_tile (Istr,Iend,Jstr,Jend,
     &                          Trate(START_2D_ARRAY,1,itrc))
        call exchange_r3d_tile (Istr,Iend,Jstr,Jend,
     &                          THmix(START_2D_ARRAY,1,itrc))
        call exchange_r3d_tile (Istr,Iend,Jstr,Jend,
     &                          TT(START_2D_ARRAY,1,itrc))            
        call exchange_r3d_tile (Istr,Iend,Jstr,Jend,
     &                          TTmoy(START_2D_ARRAY,1,itrc))
      enddo   
       
#endif  /* DIAG_SPEC */  

      end subroutine exchange_diags_spec


