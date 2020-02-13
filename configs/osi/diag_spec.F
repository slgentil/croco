
#include "cppdefs.h"
  
      module diagspec
   
#ifdef DIAG_SPEC

      implicit none
      
      integer :: LLmglo,MMmglo,LLmspec,MMmspec
      integer :: ind_lon_min_glo, ind_lon_max_glo
      integer :: ind_lat_min_glo, ind_lat_max_glo
      integer :: ind_lon_min_loc, ind_lon_max_loc
      integer :: ind_lat_min_loc, ind_lat_max_loc
      integer :: in_proc_loc
      integer :: nbpts_spec_x, nbpts_spec_y
      integer :: Nax
      integer, parameter :: nspec=19  ! nb of cospectra
      integer, parameter :: nhat=32   ! nb of variables to Fourierize
      integer :: navgspec, nwrtspec, nrecspec
      integer :: coefewper, coefnsper
      
      integer :: nint, ncidspec, ktot
      integer*8,dimension(:,:),allocatable :: lookup_index
      integer*8 :: plan
      integer*8,dimension(:),allocatable ::reqid
      integer ,dimension(:), allocatable :: index_count
      
      real(kind=8) :: dnan, zero=0.
      real(kind=8) :: lon_min, lon_max, lat_min, lat_max
      real(kind=8), parameter :: land_factor=0.05

      real(kind=8),dimension(:,:),allocatable :: zglo, zsub
      real(kind=8),dimension(:,:),allocatable :: xr_glo, yr_glo
      real(kind=8),dimension(:,:),allocatable :: zpar
      real(kind=8),dimension(:,:,:),allocatable :: zexchang
      real(kind=8),dimension(:),allocatable::kx
      real(kind=8),dimension(:),allocatable :: ky
      real(kind=8),dimension(:,:),allocatable :: kk
      real(kind=8),dimension(:),allocatable :: rho_r, drhordz
      real(kind=8),dimension(:),allocatable :: rho_r_lin, drhordz_lin
      real(kind=8),dimension(:),allocatable :: drhordz_interp
      real(kind=8),dimension(:),allocatable :: drhordz_interp_lin
      real(kind=8),dimension(:,:),allocatable :: infft 
      real(kind=8),dimension(:),allocatable :: hann_ew, hann_ns

#ifdef MASKING
      real(kind=8),dimension(:,:),allocatable :: rmask_glo
#endif
        
      character*24,dimension(:),allocatable:: name_cospec
      character*80 :: specname
      
      complex(kind=8),dimension(:,:,:),allocatable :: cospec
      complex(kind=8),dimension(:,:,:),allocatable :: hat
      complex(kind=8),dimension(:,:),allocatable :: outfft      


      contains


  !*************************************************************************
 
       subroutine init_diags_spec
     
       use isorho 

! Allocation et initialisation des tableaux
! Calcul des indices pour cumuler sur des bandes spectrales (lookup_index)
! Initialisation des FFT
! Calcul d'un rho de référence (rho_r)

# include "param.h"
# include "scalars.h"
# include "ocean3d.h"
# include "grid.h"
# include "work.h"

      include 'mpif.h'
      include 'fftw3.f'

      !local
      real(kind=8):: kh,coef,ni,nj
      integer:: i,j,k,l,kint,iproc,iret,i0,j0,mpierr,i_offset,j_offset
      integer:: nbofnan
      integer,dimension(MPI_STATUS_SIZE):: istat
      real(kind=8) :: localsum, globalsum, nr_sumofall
      
      dnan=zero/zero


         
# ifdef DIAG_SPEC_ISORHO
! Find rho min and max over the whole domain and initialize the new rho axes for future interpolation
      call init_isorho 
# endif     

# ifdef DIAG_SPEC_ISORHO
      if ( mod(Ninterp,NNODES) == 0 ) then
         nint=Ninterp/NNODES
      else
         nint=Ninterp/NNODES + 1
      endif
# else     
      if ( mod(N,NNODES) == 0 ) then
         nint=N/NNODES
      else
         nint=N/NNODES + 1
      endif
# endif


! allocation dynamique des tableaux physiques

      allocate(reqid(0:NNODES-1))
      allocate(zglo(LLm,MMm))
      allocate(xr_glo(LLm,MMm))
      allocate(yr_glo(LLm,MMm))
      allocate(zpar(Lm,Mm))
      allocate(zexchang(Lm,Mm,0:NNODES-1))
#ifdef MASKING
      allocate(rmask_glo(LLm,MMm))
#endif

! exchange xr and yr to have a global coordinate for horizontal interpolation

      call launch_rcv
      do iproc=0,NNODES-1    
# ifdef SPHERICAL
         zpar=lonr(1:Lm,1:Mm)
#else  
         zpar=xr(1:Lm,1:Mm)
#endif
         call launch_send(iproc)
      enddo
      call launch_wait   
      xr_glo=zglo(1:LLm,1:MMm)
      
      call launch_rcv
      do iproc=0,NNODES-1   
# ifdef SPHERICAL
         zpar=latr(1:Lm,1:Mm)
#else      
         zpar=yr(1:Lm,1:Mm)
#endif
         call launch_send(iproc)
      enddo
      call launch_wait   
      yr_glo=zglo(1:LLm,1:MMm)

! exchange land mask to have a global domain for horizontal interpolation

# ifdef MASKING
      call launch_rcv
      do iproc=0,NNODES-1      
         zpar=rmask(1:Lm,1:Mm)
         call launch_send(iproc)
      enddo
      call launch_wait   
      rmask_glo=zglo(1:LLm,1:MMm)

# endif  

! find index, local and global,, corresponding to lon_min, lon_max, lat_min, lat_max
      ind_lon_min_glo = minloc(abs(xr_glo(:,1)-lon_min),DIM=1)
      ind_lon_max_glo = minloc(abs(xr_glo(:,1)-lon_max),DIM=1)
      ind_lat_min_glo = minloc(abs(yr_glo(1,:)-lat_min),DIM=1)
      ind_lat_max_glo = minloc(abs(yr_glo(1,:)-lat_max),DIM=1)
      i_offset=mod(mynode,NP_XI)*Lm
      j_offset=(mynode/NP_XI)*Mm 
      in_proc_loc=1
      if (ind_lon_min_glo>imaxmpi .or. ind_lon_max_glo< iminmpi) then
         in_proc_loc=0
      else
         ind_lon_min_loc = max(1,ind_lon_min_glo - iminmpi +1)
         ind_lon_max_loc = min(Lm,ind_lon_max_glo - iminmpi +1)
      endif
      if (ind_lat_min_glo>jmaxmpi .or. ind_lat_max_glo< jminmpi) then
         in_proc_loc=0
      else
         ind_lat_min_loc = max(1,ind_lat_min_glo - jminmpi +1)
         ind_lat_max_loc = min(Mm,ind_lat_max_glo - jminmpi +1)
      endif


! Nouvelle taille de domaine et fenêtre de Hanning si domaine si non périodique    
# ifdef EW_PERIODIC     
      LLmglo=LLm
      coefewper=1.
# else
      nbpts_spec_x = ind_lon_max_glo-ind_lon_min_glo+1
      LLmglo=2*nbpts_spec_x - 1
      coefewper=2.
      allocate(hann_ew(nbpts_spec_x))
      do i=1,nbpts_spec_x
         hann_ew(i) = 0.5 - 0.5*cos(2.*pi*float(i)/float(nbpts_spec_x)) 
      enddo
# endif
# ifdef NS_PERIODIC     
      MMmglo=MMm
      coefnsper=1.
# else
      nbpts_spec_y = ind_lat_max_glo-ind_lat_min_glo+1
      MMmglo=2*nbpts_spec_y - 1
      coefnsper=2.
      allocate(hann_ns(nbpts_spec_y))      
      do j=1,nbpts_spec_y
         hann_ns(j) = 0.5 - 0.5*cos(2.*pi*float(j)/float(nbpts_spec_y)) 
      enddo
# endif
      LLmspec=(LLmglo)/2+1
      MMmspec=MMmglo 

! allocation dynamique des tableaux spectraux

      allocate(zsub(LLmglo,MMmglo))
      allocate(lookup_index(LLmspec,MMmspec))
      allocate(kx(LLmspec))
      allocate(ky(MMmspec))
      allocate(kk(LLmspec,MMmspec))
      allocate(hat(LLmspec,MMmspec,nhat)) 
      allocate(rho_r(N))
      allocate(drhordz(N)) 
      allocate(drhordz_interp(Ninterp))
      allocate(rho_r_lin(N))
      allocate(drhordz_lin(N)) 
      allocate(drhordz_interp_lin(Ninterp))
      allocate(infft(LLmglo,MMmglo))
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
      name_cospec(19)="ape_lin"        


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
  
      call dfftw_plan_dft_r2c_2d(plan,LLmglo,MMmglo,infft,outfft,
     &                           FFTW_ESTIMATE)       

      
! Initialisation d'un rho référence (pour rho du code): rho moyenné par tranche horizontal

      do k=1,N


#ifdef MASKING  
         ! Somme de rho non masqués par niveau   
         work2d=0
         where (rmask /= 0.) work2d = rho(:,:,k)
         if (in_proc_loc==1) then
            localsum=sum(work2d(ind_lon_min_loc:ind_lon_max_loc,
     &                          ind_lat_min_loc:ind_lat_max_loc))
         else
            localsum=0
         endif
         globalsum = nr_sumofall(localsum) 
         ! Nb de points non masqués
         work2d=0
         where (rmask /= 0.) work2d = dnan
         if (in_proc_loc==1) then
            localsum=count(IsNan(work2d(ind_lon_min_loc:ind_lon_max_loc,
     &                               ind_lat_min_loc:ind_lat_max_loc)))
         else
            localsum=0
         endif
         nbofnan=nr_sumofall(localsum)          
         rho_r(k) = globalsum/nbofnan
 
#else     
         localsum = sum(rho(1:Lm,1:Mm,k))/Lm/Mm
         globalsum = nr_sumofall(localsum)
         rho_r(k) = globalsum/NNODES
#endif


      enddo   
      rho_r = -1. * rho_r
      call qsort(rho_r,N)
      rho_r =  -1. * rho_r 

      do k=2,N-1
         drhordz(k) = 0.5*(rho_r(k+1) - rho_r(k-1)) / Hz(1,1,k)
      enddo
      drhordz(N) = drhordz(N-1)
      drhordz(1) = drhordz(2)  

      ! interpolate drhordz on iso-density axis
            
# ifdef DIAG_SPEC_ISORHO
      call interp1d_isorho("rho", drhordz, drhordz_interp,1,1,N,Ninterp)
# endif

      
! Initialisation d'un rho référence (pour rho linéaire): rho moyenné par tranche horizontal

      do k=1,N


#ifdef MASKING  
         ! Somme de rho non masqués par niveau   
         work2d=0
         where (rmask /= 0.) 
     &         work2d = R0 - Tcoef*(t(:,:,k,nrhs,itemp)-T0)
#  ifdef SALINITY
     &                     + Scoef*(t(:,:,k,nrhs,isalt)-S0)
#  endif
         if (in_proc_loc==1) then
            localsum=sum(work2d(ind_lon_min_loc:ind_lon_max_loc,
     &                          ind_lat_min_loc:ind_lat_max_loc))
         else
            localsum=0
         endif
         globalsum = nr_sumofall(localsum) 
         ! Nb de points non masqués
         work2d=0
         where (rmask /= 0.) work2d = dnan
         if (in_proc_loc==1) then
            localsum=count(IsNan(work2d(ind_lon_min_loc:ind_lon_max_loc,
     &                               ind_lat_min_loc:ind_lat_max_loc)))
         else
            localsum=0
         endif
         nbofnan=nr_sumofall(localsum)  
         rho_r_lin(k) = globalsum/float(nbofnan)
 
#else  /* MASKING */   
         work2d = R0 - Tcoef*(t(:,:,k,nrhs,itemp)-T0)
#  ifdef SALINITY
     &               + Scoef*(t(:,:,k,nrhs,isalt)-S0)
#  endif
         localsum = sum(work2d)/Lm/Mm
         globalsum = nr_sumofall(localsum)
         rho_r_lin(k) = globalsum/NNODES
#endif /* MASKING */  


      enddo        
      rho_r_lin = -1. * rho_r_lin
      call qsort(rho_r_lin,N)
      rho_r_lin =  -1. * rho_r_lin 

      do k=2,N-1
         drhordz_lin(k) = 0.5*(rho_r_lin(k+1) - rho_r_lin(k-1)) / Hz(1,1,k)
      enddo
      drhordz_lin(N) = drhordz_lin(N-1)
      drhordz_lin(1) = drhordz_lin(2)  

      ! interpolate drhordz on iso-density axis
            
# ifdef DIAG_SPEC_ISORHO
      call interp1d_isorho("rho", drhordz_lin, drhordz_interp_lin,1,1,
     &                      N,Ninterp)
# endif

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
!  18 = ape qg

#ifdef  DIAG_SPEC_ISORHO
       use isorho
#endif

# include "param.h"
# include "scalars.h"
# include "ocean3d.h"
# include "grid.h"
# include "work.h"

      include 'mpif.h'

 
      !local
      integer :: ik,k,kint,kv
      real(kind=8) :: nr_maxofall,nr_minofall
      
      ! échanges aux bords des sous-domaines
      call exchange_diags_spec
      
      ! find maximum of rho surface to keep only complete isorho
      ! find index on rho graduation

#ifdef MASKING     
      work2d=0
      where (rmask /= 0.) work2d = rho(:,:,N)
      if (in_proc_loc==1) then

         rhomax_surf = maxval(work2d(
     &         ind_lon_min_loc:ind_lon_max_loc,
     &         ind_lat_min_loc:ind_lat_max_loc))
      else
         rhomax_surf = dnan
      endif
#else
      rhomax_surf = maxval(rho(1:Lm,1:Mm,N)) 
#endif
      rhomax_surf = nr_maxofall(rhomax_surf)

      call find_rho_index_down(rhomax_surf,rhomax_surf_index)

      do kint=1,nint  
      
         k=(mynode*nint)+kint
         
         do kv=1,nhat
            call comp_hat(kv,kint) 
         enddo
         
        if ( any(IsNaN(real(hat(:,:,1)))) .or.
     &        k > rhomax_surf_index    ) then
           cospec(:,kint,:) = dnan
           cycle
            
         else
                                 
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
            call comp_cospec(32,32,18,kint) 

! density, linear ape QG (rho_lin*rho_lin)
            call comp_cospec(23,23,19,kint) 

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

#ifdef  DIAG_SPEC_ISORHO
       use isorho
#endif     
 
# include "param.h" 
# include "scalars.h"
     
      integer::kv,kint   
      
      ! compute complex amplitude for variable kv and level kint

	include 'mpif.h'
        include 'fftw3.f'

      !local
      integer::iproc,i,j,k, nbofnan
      real :: mean, rho_mean
      
      k=(mynode*nint)+kint
      if (k <= rhomax_surf_index) then
         call launch_rcv
      endif
      do iproc=0,NNODES-1      
         ktot=(iproc*nint)+kint
         if (ktot <= rhomax_surf_index) then
            call compute_var(kv,ktot)
            call launch_send(iproc)
         endif
      enddo
      if (k <= rhomax_surf_index) then
         call launch_wait   
      endif 

!     array is Nan in place of the land mask
#ifdef MASKING
      where (rmask_glo == 0.) zglo(1:LLm,1:MMm) = dnan
#endif

! extract sub-domain for spectral diagnostics
      zsub(1:nbpts_spec_x,1:nbpts_spec_y) = 
     &             zglo(ind_lon_min_glo:ind_lon_max_glo,
     &                  ind_lat_min_glo:ind_lat_max_glo)

!     2d horizontal interpolation if pourcentage of land acceptable
      nbofnan = count(IsNan(zsub(1:nbpts_spec_x,1:nbpts_spec_y)))   
      if (nbofnan/=0 .and. nbofnan < land_factor*LLm*MMm) then
!        call fillnan(zglo(1:LLM,1:MMm),xr_glo(1:LLm,1),yr_glo(1,1:MMm),
!    &                 LLm,MMm)
         call fillnan(zsub(1:nbpts_spec_x,1:nbpts_spec_y),
     &                xr_glo(ind_lon_min_glo:ind_lon_max_glo,1),
     &                yr_glo(1,ind_lat_min_glo:ind_lat_max_glo),
     &                 nbpts_spec_x,nbpts_spec_y)

      endif
     
! doesn't treat the level with outcropping density interfaces
      if ((k <=rhomax_surf_index) .and.
     &    (.not. any(IsNaN(zsub(1:nbpts_spec_x,1:nbpts_spec_y)))) ) then
      
! double le domaine si nécessaire pour le rendre périodique  
! applique une fenêtre de hanning pour atténuer les bords
   
# ifndef EW_PERIODIC
         do j=1,MMmglo
            zsub(1:nbpts_spec_x,j) = zsub(1:nbpts_spec_x,j)*hann_ew(:)
         enddo
         do i=1,nbpts_spec_x-1
            zsub(nbpts_spec_x+i,:) = zsub(nbpts_spec_x-i,:)
         enddo
# endif
# ifndef NS_PERIODIC
         do i=1,LLmglo
            zsub(i,1:nbpts_spec_y) = zsub(i,1:nbpts_spec_y)*hann_ns(:)
         enddo
         do j=1,nbpts_spec_y-1
            zsub(:,nbpts_spec_y+j) = zsub(:,nbpts_spec_y-j)
         enddo
# endif 
       
! fourierize
         infft = zsub  
         call dfftw_execute_dft_r2c(plan,infft,outfft)
         hat(:,:,kv)=outfft/(LLmglo)/(MMmglo)
      
      else
         hat(:,:,kv)=dnan
      endif
  
      end subroutine comp_hat


!========================================
      subroutine compute_var(var,k)
      
      use isorho
     
      integer::var,k            
        
# include "param.h"
# include "scalars.h"
# include "ocean3d.h"
# include "diagnostics.h"
# include "grid.h"
 
      include 'mpif.h'


      !local
      integer::i,j,tile,kk
      real :: rhodiag, mean
      real, dimension(N) :: val_in
      real, dimension(GLOBAL_2D_ARRAY,1:N) :: workr
         

      select case(var)

!---------------------------------------------------------------
! Momentum diagnostic variables.
!---------------------------------------------------------------

      case(1)
!     compute u at level k        
         do j=1,Mm
            do i=1,Lm
#ifdef DIAG_SPEC_ISORHO  
               val_in = MV(i,j,:,1) 
               rho_axis(:) = 0.5*(rho(i-1,j,:)+rho(i,j,:)) 
               call interp0d(val_in,rho_axis,zpar(i,j),
     &              rho_axis_interp(k),method,N)
!              call  interp0d_isorho("u",val_in,zpar(i,j),i,j,k,N,Ninterp)                      
#else
               zpar(i,j)=MV(i,j,k,1)
#endif
            enddo
         enddo
         
      case(2)
!     compute v at level k
         do j=1,Mm
            do i=1,Lm
#ifdef DIAG_SPEC_ISORHO 
               val_in = MV(i,j,:,2)   
               rho_axis(:) = 0.5*(rho(i,j-1,:)+rho(i,j,:)) 
               call interp0d(val_in,rho_axis,zpar(i,j),
     &              rho_axis_interp(k),method,N)
!              call  interp0d_isorho("v",val_in,zpar(i,j),i,j,k,N,Ninterp)            
#else
               zpar(i,j)=MV(i,j,k,2)
#endif
            enddo
         enddo
         
      case(3)
!     compute horizontal advection of u at level kl
         do j=1,Mm
            do i=1,Lm
#ifdef DIAG_SPEC_ISORHO  
               val_in =  MXadv(i,j,:,1)+MYadv(i,j,:,1) 
               rho_axis(:) = 0.5*(rho(i-1,j,:)+rho(i,j,:))
               call interp0d(val_in,rho_axis,zpar(i,j),
     &              rho_axis_interp(k),method,N)
!              call  interp0d_isorho("u",val_in,zpar(i,j),i,j,k,N,Ninterp)             
#else
               zpar(i,j)=MXadv(i,j,k,1)+MYadv(i,j,k,1)
#endif
            enddo
         enddo
         
      case(4)
!     compute vertical advection of u at level kl
         do j=1,Mm
            do i=1,Lm
#ifdef DIAG_SPEC_ISORHO  
               val_in = MVadv(i,j,:,1)
               rho_axis(:) = 0.5*(rho(i-1,j,:)+rho(i,j,:)) 
               call interp0d(val_in,rho_axis,zpar(i,j),
     &              rho_axis_interp(k),method,N)  
!              call  interp0d_isorho("u",val_in,zpar(i,j),i,j,k,N,Ninterp)          
#else
               zpar(i,j)=MVadv(i,j,k,1)
#endif
            enddo
         enddo
 
      case(5)
!     compute horizontal advection of u C4 scheme at level kl
         do j=1,Mm
            do i=1,Lm
#ifdef DIAG_SPEC_ISORHO  
               val_in =  MXadvC(i,j,:,1)+MYadvC(i,j,:,1)
               rho_axis(:) = 0.5*(rho(i-1,j,:)+rho(i,j,:))
               call interp0d(val_in,rho_axis,zpar(i,j),
     &              rho_axis_interp(k),method,N)
!              call  interp0d_isorho("u",val_in,zpar(i,j),i,j,k,N,Ninterp)     
#else
               zpar(i,j)=MXadvC(i,j,k,1)+MYadvC(i,j,k,1)
#endif
            enddo
         enddo 
                
      case(6)
!     compute horizontal advection of v at level kl
         do j=1,Mm
            do i=1,Lm
#ifdef DIAG_SPEC_ISORHO 
               val_in =  MXadv(i,j,:,2)+MYadv(i,j,:,2)    
               rho_axis(:) = 0.5*(rho(i,j-1,:)+rho(i,j,:)) 
               call interp0d(val_in,rho_axis,zpar(i,j),
     &              rho_axis_interp(k),method,N) 
!              call  interp0d_isorho("v",val_in,zpar(i,j),i,j,k,N,Ninterp)            
#else
               zpar(i,j)=MXadv(i,j,k,2)+MYadv(i,j,k,2)
#endif
            enddo
         enddo   
      
      case(7)
!     compute vertical advection of v at level k
         do j=1,Mm
            do i=1,Lm
#ifdef DIAG_SPEC_ISORHO 
               val_in = MVadv(i,j,:,2)    
               rho_axis(:) = 0.5*(rho(i,j-1,:)+rho(i,j,:)) 
               call interp0d(val_in,rho_axis,zpar(i,j),
     &              rho_axis_interp(k),method,N)
!              call  interp0d_isorho("v",val_in,zpar(i,j),i,j,k,N,Ninterp)             
#else
               zpar(i,j)=MVadv(i,j,k,2)
#endif
            enddo
         enddo
         
      case(8)
!     compute horizontal advection of v C4 scheme at level kl
         do j=1,Mm
            do i=1,Lm
#ifdef DIAG_SPEC_ISORHO 
               val_in =  MXadvC(i,j,:,2)+MYadvC(i,j,:,2)     
               rho_axis(:) = 0.5*(rho(i,j-1,:)+rho(i,j,:)) 
               call interp0d(val_in,rho_axis,zpar(i,j),
     &              rho_axis_interp(k),method,N)
!              call  interp0d_isorho("v",val_in,zpar(i,j),i,j,k,N,Ninterp)           
#else
               zpar(i,j)=MXadvC(i,j,k,2)+MYadvC(i,j,k,2)
#endif
            enddo
         enddo   
         
      case(9)
!     compute Coriolis part of u at level k
         do j=1,Mm
            do i=1,Lm
#ifdef DIAG_SPEC_ISORHO  
               val_in = MCor(i,j,:,1) 
               rho_axis(:) = 0.5*(rho(i-1,j,:)+rho(i,j,:)) 
               call interp0d(val_in,rho_axis,zpar(i,j),
     &              rho_axis_interp(k),method,N) 
!              call  interp0d_isorho("u",val_in,zpar(i,j),i,j,k,N,Ninterp) 
#else
               zpar(i,j)=MCor(i,j,k,1)
#endif
            enddo
         enddo
         
      case(10)
!     compute coriolis part of v at level k
         do j=1,Mm
            do i=1,Lm
#ifdef DIAG_SPEC_ISORHO 
               val_in = MCor(i,j,:,2)    
               rho_axis(:) = 0.5*(rho(i,j-1,:)+rho(i,j,:)) 
               call interp0d(val_in,rho_axis,zpar(i,j),
     &              rho_axis_interp(k),method,N)    
!              call  interp0d_isorho("v",val_in,zpar(i,j),i,j,k,N,Ninterp)   
#else
               zpar(i,j)=MCor(i,j,k,2)
#endif
            enddo
         enddo  
                       
      case(11)
!     compute Pressure gradient of u at level k
         do j=1,Mm
            do i=1,Lm
#ifdef DIAG_SPEC_ISORHO  
               val_in =   MPrsgrd(i,j,:,1) 
               rho_axis(:) = 0.5*(rho(i-1,j,:)+rho(i,j,:))
               call interp0d(val_in,rho_axis,zpar(i,j),
     &              rho_axis_interp(k),method,N)
!              call  interp0d_isorho("u",val_in,zpar(i,j),i,j,k,N,Ninterp)
#else
               zpar(i,j)=MPrsgrd(i,j,k,1)
#endif
            enddo
         enddo
         
         
      case(12)
!     compute Pressure gradient of v at level k
         do j=1,Mm
            do i=1,Lm
#ifdef DIAG_SPEC_ISORHO 
               val_in = MPrsgrd(i,j,:,2)   
               rho_axis(:) = 0.5*(rho(i,j-1,:)+rho(i,j,:)) 
               call interp0d(val_in,rho_axis,zpar(i,j),
     &              rho_axis_interp(k),method,N) 
!              call  interp0d_isorho("v",val_in,zpar(i,j),i,j,k,N,Ninterp)
#else
               zpar(i,j)=MPrsgrd(i,j,k,2)
#endif
            enddo
         enddo 
         
      case(13)        
!     compute Pressure gradient of w at level k
         do j=1,Mm
            do i=1,Lm
#ifdef DIAG_SPEC_ISORHO 
               val_in = g / rho0 * rho(i,j,:)
               rho_axis(:) = rho(i,j,:)
               call interp0d(val_in,rho_axis,zpar(i,j),
     &              rho_axis_interp(k),method,N)
!              call  interp0d_isorho("rho",val_in,zpar(i,j),i,j,k,N,Ninterp)
#else
               zpar(i,j)= g / rho0 * rho(i,j,k)
#endif
            enddo
         enddo  
                    
      case(14)
!     compute  horizontal mixing of u at level k

         do j=1,Mm
            do i=1,Lm
#ifdef DIAG_SPEC_ISORHO  
               val_in = MHmix(i,j,:,1)
               rho_axis(:) = 0.5*(rho(i-1,j,:)+rho(i,j,:))
               call interp0d(val_in,rho_axis,zpar(i,j),
     &              rho_axis_interp(k),method,N)
!              call  interp0d_isorho("u",val_in,zpar(i,j),i,j,k,N,Ninterp)
#else
               zpar(i,j)=MHmix(i,j,k,1)
#endif
            enddo
         enddo
         
      case(15)
!     compute  vertical mixing of u at level k

         do j=1,Mm
            do i=1,Lm
#ifdef DIAG_SPEC_ISORHO  
               val_in = MVmix(i,j,:,1)
               rho_axis(:) = 0.5*(rho(i-1,j,:)+rho(i,j,:))
               call interp0d(val_in,rho_axis,zpar(i,j),
     &              rho_axis_interp(k),method,N)
!              call  interp0d_isorho("u",val_in,zpar(i,j),i,j,k,N,Ninterp)
#else
               zpar(i,j)=MVmix(i,j,k,1)
#endif
            enddo
         enddo
                 
      case(16)
!     compute  horizontal mixing of v at level k
         do j=1,Mm
            do i=1,Lm
#ifdef DIAG_SPEC_ISORHO 
               val_in =  MHmix(i,j,:,2)   
               rho_axis(:) = 0.5*(rho(i,j-1,:)+rho(i,j,:))
               call interp0d(val_in,rho_axis,zpar(i,j),
     &              rho_axis_interp(k),method,N)
!              call  interp0d_isorho("v",val_in,zpar(i,j),i,j,k,N,Ninterp)
#else
               zpar(i,j)=MHmix(i,j,k,2)
#endif
            enddo
         enddo 
           
      case(17)
!     compute  vertical mixing of v at level k
         do j=1,Mm
            do i=1,Lm
#ifdef DIAG_SPEC_ISORHO 
               val_in =  MVmix(i,j,:,2)   
               rho_axis(:) = 0.5*(rho(i,j-1,:)+rho(i,j,:))
               call interp0d(val_in,rho_axis,zpar(i,j),
     &              rho_axis_interp(k),method,N)
!              call  interp0d_isorho("v",val_in,zpar(i,j),i,j,k,N,Ninterp)
#else
               zpar(i,j)=MVmix(i,j,k,2)
#endif
            enddo
         enddo           
 
      case(18)
!     compute dudt at level k      
         do j=1,Mm
            do i=1,Lm
#ifdef DIAG_SPEC_ISORHO  
               val_in = Mrate(i,j,:,1)
               rho_axis(:) = 0.5*(rho(i-1,j,:)+rho(i,j,:))
               call interp0d(val_in,rho_axis,zpar(i,j),
     &              rho_axis_interp(k),method,N)
!              call  interp0d_isorho("u",val_in,zpar(i,j),i,j,k,N,Ninterp)
#else
               zpar(i,j)=Mrate(i,j,k,1)
#endif
            enddo
         enddo             
                         
      case(19)
!     compute dvdt at level k
         do j=1,Mm
            do i=1,Lm
#ifdef DIAG_SPEC_ISORHO 
               val_in = Mrate(i,j,:,2)   
               rho_axis(:) = 0.5*(rho(i,j-1,:)+rho(i,j,:))
               call interp0d(val_in,rho_axis,zpar(i,j),
     &              rho_axis_interp(k),method,N)
!              call  interp0d_isorho("v",val_in,zpar(i,j),i,j,k,N,Ninterp)
#else
               zpar(i,j)=Mrate(i,j,k,2)
#endif
            enddo
         enddo

      case(20)
!     compute (u(n+1)+u(n))/2 at level k
         do j=1,Mm
            do i=1,Lm
#ifdef DIAG_SPEC_ISORHO  
               val_in = MVmoy(i,j,:,1)
               rho_axis(:) = 0.5*(rho(i-1,j,:)+rho(i,j,:))
               call interp0d(val_in,rho_axis,zpar(i,j),
     &              rho_axis_interp(k),method,N)
!              call  interp0d_isorho("u",val_in,zpar(i,j),i,j,k,N,Ninterp)
#else
               zpar(i,j)=MVmoy(i,j,k,1)
#endif
            enddo
         enddo   
         
      case(21)
!     compute (v(n+1)+v(n))/2 at level k
         do j=1,Mm
            do i=1,Lm
#ifdef DIAG_SPEC_ISORHO 
               val_in = MVmoy(i,j,:,2)   
               rho_axis(:) = 0.5*(rho(i,j-1,:)+rho(i,j,:))
               call interp0d(val_in,rho_axis,zpar(i,j),
     &              rho_axis_interp(k),method,N)
!              call  interp0d_isorho("v",val_in,zpar(i,j),i,j,k,N,Ninterp)
#else
               zpar(i,j)=MVmoy(i,j,k,2)
#endif
            enddo
         enddo 
                         
      case(22)
!     compute w at level k
         do tile=0,NSUB_X*NSUB_E-1
            call Wvlcty (tile, workr)
         enddo
         do j=1,Mm
            do i=1,Lm
#ifdef DIAG_SPEC_ISORHO 
               val_in = workr(i,j,:)
               rho_axis(:) = rho(i,j,:)
               call interp0d(val_in,rho_axis,zpar(i,j),
     &              rho_axis_interp(k),method,N)
!              call  interp0d_isorho("rho",val_in,zpar(i,j),i,j,k,N,Ninterp)
#else
               zpar(i,j)=workr(i,j,k)
#endif
            enddo
         enddo
         
!---------------------------------------------------------------
! Tracer diagnostic variables.
!---------------------------------------------------------------
         
      case(23)         
!     compute rho at level k
         do j=1,Mm
            do i=1,Lm
#ifdef DIAG_SPEC_ISORHO 
               val_in = rho(i,j,:) - rho_r(:)   
               rho_axis(:) = rho(i,j,:)
               call interp0d(val_in,rho_axis,zpar(i,j),
     &              rho_axis_interp(k),method,N)
!              call  interp0d_isorho("rho",val_in,zpar(i,j),i,j,k,N,Ninterp)  
#else
               rhodiag= rho(i,j,k)    
               zpar(i,j)=rhodiag - rho_r(k)  
#endif
            enddo
         enddo
            
      case(24)
!     compute horizontal advection of rho at level k
         do j=1,Mm
            do i=1,Lm
#ifdef DIAG_SPEC_ISORHO 
               val_in =  					  				
     &                 - (TXadv(i,j,:,itemp) +			 
     &                    TYadv(i,j,:,itemp) )*Tcoef		 
# ifdef SALINITY                      
     &                 + (TXadv(i,j,:,isalt) +			 
     &                    TYadv(i,j,:,isalt) )*Scoef 		 
# endif   
               rho_axis(:) = rho(i,j,:)
               call interp0d(val_in,rho_axis,zpar(i,j),
     &              rho_axis_interp(k),method,N)
!              call  interp0d_isorho("rho",val_in,zpar(i,j),i,j,k,N,Ninterp) 
#else
               zpar(i,j)= 					  				
     &                 - (TXadv(i,j,k,itemp) +			 
     &                    TYadv(i,j,k,itemp) )*Tcoef		 
# ifdef SALINITY                      
     &                 + (TXadv(i,j,k,isalt) +			 
     &                    TYadv(i,j,k,isalt) )*Scoef 		 
# endif
#endif
            enddo
         enddo

      case(25)
!     compute vertical advection of rho at level k
         do j=1,Mm
            do i=1,Lm
#ifdef DIAG_SPEC_ISORHO 
               val_in =  					  				
     &                 - TVadv(i,j,:,itemp)*Tcoef		 
# ifdef SALINITY                      
     &                 + TVadv(i,j,:,isalt)*Scoef 		 
# endif	 
               rho_axis(:) = rho(i,j,:)
               call interp0d(val_in,rho_axis,zpar(i,j),
     &              rho_axis_interp(k),method,N)
!              call  interp0d_isorho("rho",val_in,zpar(i,j),i,j,k,N,Ninterp)
#else
               zpar(i,j)= 					  				
     &                 - TVadv(i,j,:,itemp)*Tcoef		 
# ifdef SALINITY                      
     &                 + TVadv(i,j,:,isalt)*Scoef 		 
# endif
#endif
            enddo
         enddo
 
      case(26)
!     compute horizontal advection of t & s at level k C4 or C6 scheme
         do j=1,Mm
            do i=1,Lm
#ifdef DIAG_SPEC_ISORHO 
               val_in =  					  				
     &                 - (TXadvC(i,j,:,itemp) +			 
     &                    TYadvC(i,j,:,itemp) )*Tcoef		 
# ifdef SALINITY                      
     &                 + (TXadvC(i,j,:,isalt) +			 
     &                    TYadvC(i,j,:,isalt) )*Scoef 		 
# endif
               rho_axis(:) = rho(i,j,:)
               call interp0d(val_in,rho_axis,zpar(i,j),
     &              rho_axis_interp(k),method,N)
!              call  interp0d_isorho("rho",val_in,zpar(i,j),i,j,k,N,Ninterp)
#else
               zpar(i,j)= 					  				
     &                 - (TXadvC(i,j,k,itemp) +			 
     &                    TYadvC(i,j,k,itemp) )*Tcoef		 
# ifdef SALINITY                      
     &                 + (TXadvC(i,j,k,isalt) +			 
     &                    TYadvC(i,j,k,isalt) )*Scoef 		 
# endif
#endif
            enddo
         enddo
                      
         
      case(27)
!     compute horizontal mixing of rho at level k
         do j=1,Mm
            do i=1,Lm
#ifdef DIAG_SPEC_ISORHO 
               val_in =  					  				
     &                 - THmix(i,j,:,itemp)*Tcoef		 
# ifdef SALINITY                      
     &                 + THmix(i,j,:,isalt)*Scoef 		 
# endif
               rho_axis(:) = rho(i,j,:)
               call interp0d(val_in,rho_axis,zpar(i,j),
     &              rho_axis_interp(k),method,N)
!              call  interp0d_isorho("rho",val_in,zpar(i,j),i,j,k,N,Ninterp)
#else
               zpar(i,j)= 					  				
     &                 - THmix(i,j,k,itemp)*Tcoef		 
# ifdef SALINITY                      
     &                 + THmix(i,j,k,isalt)*Scoef 		 
# endif
#endif
            enddo
         enddo
         
      case(28)
!     compute vertical mixing of rho at level k
         do j=1,Mm
            do i=1,Lm
#ifdef DIAG_SPEC_ISORHO 
               val_in =  					  				
     &                 - TVmix(i,j,:,itemp)*Tcoef		 
# ifdef SALINITY                      
     &                 + TVmix(i,j,:,isalt)*Scoef 		  
# endif
               rho_axis(:) = rho(i,j,:)
               call interp0d(val_in,rho_axis,zpar(i,j),
     &              rho_axis_interp(k),method,N)
!              call  interp0d_isorho("rho",val_in,zpar(i,j),i,j,k,N,Ninterp)
#else
               zpar(i,j)= 					  				
     &                 - TVmix(i,j,k,itemp)*Tcoef		 
# ifdef SALINITY                      
     &                 + TVmix(i,j,k,isalt)*Scoef 		  
# endif
#endif
            enddo
         enddo

      case(29)
!     compute forcing of rho at level k
         do j=1,Mm
            do i=1,Lm
#ifdef DIAG_SPEC_ISORHO 
               val_in =  					  				
     &                 - TForc(i,j,:,itemp)*Tcoef		 
# ifdef SALINITY                      
     &                 + TForc(i,j,:,isalt)*Scoef 		 
# endif
               rho_axis(:) = rho(i,j,:)
               call interp0d(val_in,rho_axis,zpar(i,j),
     &              rho_axis_interp(k),method,N)
!              call  interp0d_isorho("rho",val_in,zpar(i,j),i,j,k,N,Ninterp)
#else
               zpar(i,j)= 					  				
     &                 - TForc(i,j,k,itemp)*Tcoef		 
# ifdef SALINITY                      
     &                 + TForc(i,j,k,isalt)*Scoef 		 
# endif
#endif
            enddo
         enddo
   
      case(30)
!     compute drhodt at level k
         do j=1,Mm
            do i=1,Lm
#ifdef DIAG_SPEC_ISORHO 
               val_in =  					  				
     &                 - Trate(i,j,:,itemp)*Tcoef		 
# ifdef SALINITY                      
     &                 + Trate(i,j,:,isalt)*Scoef 		 
# endif
               rho_axis(:) = rho(i,j,:)
               call interp0d(val_in,rho_axis,zpar(i,j),
     &              rho_axis_interp(k),method,N)
!              call  interp0d_isorho("rho",val_in,zpar(i,j),i,j,k,N,Ninterp)
#else
               zpar(i,j)= 					  				
     &                 - Trate(i,j,k,itemp)*Tcoef		 
# ifdef SALINITY                      
     &                 + Trate(i,j,k,isalt)*Scoef 		 
# endif
#endif
            enddo
         enddo
 
      case(31)         
!     compute (rho(N)+rho(N+1))/2 at level k
         do j=1,Mm
            do i=1,Lm
#ifdef DIAG_SPEC_ISORHO 
               val_in = R0 - Tcoef*(TTmoy(i,j,:,itemp)-T0)
#  ifdef SALINITY
     &                     + Scoef*(TTmoy(i,j,:,isalt)-S0)
#  endif  
     &                     - rho_r(:)
               rho_axis(:) = rho(i,j,:)
               call interp0d(val_in,rho_axis,zpar(i,j),
     &              rho_axis_interp(k),method,N)
!              call  interp0d_isorho("rho",val_in,zpar(i,j),i,j,k,N,Ninterp)
#else
               rhodiag = R0 - Tcoef*(TTmoy(i,j,k,itemp)-T0)
#  ifdef SALINITY
     &                     + Scoef*(TTmoy(i,j,k,isalt)-S0)
#  endif      
               zpar(i,j) = rhodiag - rho_r(k)  
#endif
            enddo
         enddo  

      case(32)         
!     compute rho_lin at level k
         do j=1,Mm
            do i=1,Lm
#ifdef DIAG_SPEC_ISORHO 
               val_in = R0 - Tcoef*(TT(i,j,:,itemp)-T0)
#  ifdef SALINITY
     &                     + Scoef*(TT(i,j,:,isalt)-S0)
#  endif
     &                	   - rho_r_lin(:)   
               rho_axis(:) = rho(i,j,:)
               call interp0d(val_in,rho_axis,zpar(i,j),
     &              rho_axis_interp(k),method,N)
!              call  interp0d_isorho("rho",val_in,zpar(i,j),i,j,k,N,Ninterp)  
#else
               rhodiag= R0 - Tcoef*(TT(i,j,k,itemp)-T0)
#  ifdef SALINITY
     &                     + Scoef*(TT(i,j,k,isalt)-S0)
#  endif     
               zpar(i,j)=rhodiag - rho_r_lin(k)  
#endif
            enddo
         enddo                                 
      end select
      
      end subroutine compute_var

!========================================

      subroutine launch_rcv
                 
# include "param.h"
# include "scalars.h"
# include "mpi_cpl.h"

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
# include "mpi_cpl.h"

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
      use isorho    
      
      integer:: iflag           

# include "param.h"
# include "scalars.h"
# include "ocean3d.h"
# include "grid.h"
# include "netcdf.inc"
# include "mpi_cpl.h"

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
         call check(nf_def_dim(ncidspec, 'level',Ninterp,dimid(2)) )
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
          cospec(:,:,19) = 0.5 * cospec(:,:,19) 
                   
     
          do l=11,nspec-1            
             do k=1,nint
                ktot=mynode*nint+k               
                if (ktot <= Ninterp) then             
                     cospec(:,k,l) = - g/drhordz_interp(ktot)/rho0*
     &                               cospec(:,k,l)
                endif
             enddo
          enddo       
          do k=1,nint
             ktot=mynode*nint+k               
             if (ktot <= Ninterp) then             
                  cospec(:,k,19) = - g/drhordz_interp_lin(ktot)/rho0*
     &                            cospec(:,k,19)
             endif
          enddo

! cospec divided because of averaging, duplicate domain EW, duplicate domain NS          

          cospec = cospec / float(nwrtspec/navgspec) / coefewper / coefnsper     

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
                    if (Nin <= Ninterp) then                    
!                      call MPI_SEND(cospec,Nax*nint*nspec,
                       call MPI_SEND(cospec(:,:,l),Nax*nint,                
     &                       MPI_DOUBLE_COMPLEX,0,mynode,  
     &                       MPI_COMM_WORLD,mpierr)  
                    endif
             endif                     
               
             if (mynode==0) then
                  do iproc=1,NNODES-1
                     Nin=iproc*nint+1   
                     if (Nin <= Ninterp) then                    
!                       call MPI_IRECV(cospec,  
!    &                      Nax*nint*nspec, 
                        call MPI_IRECV(cospec(:,:,l),  
     &                      Nax*nint,  
     &                      MPI_DOUBLE_COMPLEX,iproc,iproc,  
     &                      MPI_COMM_WORLD,reqid(iproc),mpierr)
                        call MPI_WAIT(reqid(iproc),istat,mpierr)
                        kcount = min(nint,Ninterp-Nin+1)                         
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
       
       
!======================================== 

      subroutine exchange_diags_spec()    


# include "param.h"
# include "scalars.h"
# include "ocean3d.h"
# include "diagnostics.h"

      include 'mpif.h'


      !local
      integer::i,j,itrc,tile
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
       
      end subroutine exchange_diags_spec


#endif  /* DIAG_SPEC */  
    
      end module diagspec
