
#include "cppdefs.h"
  
      module diagspec_kt
   
#ifdef DIAG_SPEC_KT

      implicit none
      
      integer :: LLmspec,MMmspec
      integer, parameter :: nhat_kt=4   ! nb of variables to Fourierize
      integer, parameter :: Nmax_kt=10  ! Max number of level to save
      integer :: navgspec_kt, nwrtspec_kt, nrecspec_kt
      integer :: level_kt(Nmax_kt)
      
      integer :: nint, ncidspec_kt, ktot, N_kt
      integer*8 :: plan_kt
      integer*8,dimension(:),allocatable ::reqid
      
      real(kind=8),dimension(Nmax_kt) :: depth_kt
      real(kind=8),dimension(:,:),allocatable :: zglo
      real(kind=8),dimension(:,:),allocatable :: zpar
      real(kind=8),dimension(:,:,:),allocatable :: zexchang
      real(kind=8),dimension(:),allocatable::kx
      real(kind=8),dimension(:),allocatable :: ky
      real(kind=8),dimension(:,:),allocatable :: kk
      real(kind=8),dimension(:),allocatable :: rho_r_kt
      real(kind=8),dimension(:,:),allocatable :: infft 
         
      character*24,dimension(:),allocatable:: name_hat_kt
      character *80 :: specname_kt
      
      complex(kind=8),dimension(:,:,:,:),allocatable :: hat_kt
      complex(kind=8),dimension(:,:),allocatable :: outfft      


      contains


  !*************************************************************************
 
       subroutine init_diags_spec_kt
      

! Allocation et initialisation des tableaux
! Initialisation des FFT
! Calcul d'un rho de référence (rho_r_kt)

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
      if ( mod(N_kt,NNODES) == 0 ) then
         nint=N_kt/NNODES
      else
         nint=N_kt/NNODES + 1
      endif
!     Find indexes corresponding to nearest depths compute
      do k=1,N_kt
         level_kt(k) = minloc(abs(z_r(1,1,:) - depth_kt(k)),1)
      enddo


! allocation dynamique des tableaux

      allocate(reqid(0:NNODES-1))
      allocate(zglo(LLm,MMm))
      allocate(zpar(Lm,Mm))
      allocate(zexchang(Lm,Mm,0:NNODES-1))
      allocate(kx(LLmspec))
      allocate(ky(MMmspec))
      allocate(kk(LLmspec,MMmspec))
      allocate(hat_kt(LLmspec,MMmspec,N_kt,nhat_kt)) 
      allocate(rho_r_kt(N))
      allocate(infft(LLm,MMm))
      allocate(outfft(LLmspec,MMmspec)) 

! Initialisation des noms de variables pour le fichier netcdf      
      allocate(name_hat_kt(nhat_kt))       
      name_hat_kt(1)="u "
      name_hat_kt(2)="v"
      name_hat_kt(3)="w"
      name_hat_kt(4)="rho"  

!     Make wave number mask

      call get_kxky_kt         
      
      nrecspec_kt=0
      hat_kt= 0.

! output paramètres diags spectraux
 
      if(mynode.eq.0)then
         write(*,*)
         write(*,*)"====== Diags spectraux temporels ======"
         write(*,*)
         write(*,*)" Nombre de niveaux traités par processeur, nint  = ",nint  
         write(*,*)" Nombre de variables fourierisées, nhat_kt  = ",nhat_kt           
         write(*,*)" Fichier netcdf = diag_spec_kt.nc"
         write(*,*)" Variables fourierisées = "   
         do l=1,nhat_kt
            write(*,*)"           ",trim(name_hat_kt(l))          
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
 
      call dfftw_plan_dft_r2c_2d(plan_kt,LLm,MMm,infft,outfft,
     &                           FFTW_ESTIMATE)       

      
! Initialisation d'un rho référence: rho moyenné par tranche horizontal    
  
      do k=1,N
         localsum = sum(rho(1:Lm,1:Mm,k))/Lm/Mm
         globalsum = nr_sumofall(localsum)
         rho_r_kt(k) = globalsum/NNODES
      enddo          
      rho_r_kt = -1. * rho_r_kt
      call qsort(rho_r_kt,N)
      rho_r_kt =  -1. * rho_r_kt 


      end subroutine init_diags_spec_kt
      
!========================================

      subroutine get_kxky_kt
      
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
      
      end subroutine get_kxky_kt

!========================================

      subroutine comp_diags_spec_kt
      
       
! Calcul des variables spectrales:
!	1 = u(n+1/2)
!	2 = v(n+1/2)
!	3 = w
!	4 = rho(n+1/2)


# include "param.h"
# include "scalars.h"

	include 'mpif.h'

 
      !local
      integer :: ik,k,kint,kv


      ! échanges aux bords des sous-domaines
      call exchange_diags_spec_kt

      do kint=1,nint  
      
         do kv=1,nhat_kt
            call comp_hat_kt(kv,kint) 
         enddo

      enddo
      
      end subroutine  comp_diags_spec_kt
      
!========================================

      subroutine comp_hat_kt(kv,kint)
 
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
      if (k <= N_kt) then
         call launch_rcv_kt
      endif
      do iproc=0,NNODES-1      
         ktot=(iproc*nint)+kint
         if (ktot <= N_kt) then
            call compute_var_kt(kv,ktot)
            call launch_send_kt(iproc)
         endif
      enddo
      if (k <= N_kt) then
         call launch_wait_kt   
      endif       
      
             
      if (k <= N_kt) then    
       
! fourierize
         infft = zglo  
         call dfftw_execute_dft_r2c(plan_kt,infft,outfft)
         hat_kt(:,:,kint,kv)=hat_kt(:,:,kint,kv)+outfft/(LLm)/(MMm)

      endif



      end subroutine comp_hat_kt


!========================================
      subroutine compute_var_kt(var,ik)
     
      integer::var,ik            
        
# include "param.h"
# include "scalars.h"
# include "ocean3d.h"
# include "diagnostics.h"

      include 'mpif.h'


      !local
      integer::i,j,tile,k
      real :: rhodiag, mean
      real, dimension(GLOBAL_2D_ARRAY,1:N) :: workr
         

      k = level_kt(ik)

      select case(var)

!---------------------------------------------------------------
! Momentum diagnostic variables.
!---------------------------------------------------------------

      case(1)
!     compute (u(n+1)+u(n))/2 at level k
         do j=1,Mm
            do i=1,Lm
               zpar(i,j)=MVmoy(i,j,k,1)
            enddo
         enddo   
         
      case(2)
!     compute (v(n+1)+v(n))/2 at level k
         do j=1,Mm
            do i=1,Lm
               zpar(i,j)=MVmoy(i,j,k,2)
            enddo
         enddo 
                           
      case(3)
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
         
      case(4)         
         do j=1,Mm
            do i=1,Lm
               rhodiag = R0 - Tcoef*(TTmoy(i,j,k,itemp)-T0)
#  ifdef SALINITY
     &                     + Scoef*(TTmoy(i,j,k,isalt)-S0)
#  endif     
               zpar(i,j) = rhodiag - rho_r_kt(k)  
            enddo
         enddo 

      end select
      
      end subroutine compute_var_kt

!========================================

      subroutine launch_rcv_kt
                 
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

      
      end subroutine launch_rcv_kt

!========================================

      subroutine launch_wait_kt
                 
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
      
      end subroutine launch_wait_kt

!========================================

      subroutine launch_send_kt(iproc)

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


      end subroutine launch_send_kt

!======================================== 

      subroutine wrt_diags_spec_kt(iflag)
      
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

         specname_kt = 'diag_spec_kt.nc' 
         call check(nf_create(trim(specname_kt), NF_CLOBBER, ncidspec_kt)  )    

         call check(nf_def_dim(ncidspec_kt, 'kx',LLmspec,dimid(1)) )   
         call check(nf_def_dim(ncidspec_kt, 'ky',MMmspec,dimid(2)) )      
         call check(nf_def_dim(ncidspec_kt, 'level',N_kt,dimid(3)) )     
         nbtime = int((ntimes-ntstart)/nwrtspec_kt)+1
         call check(nf_def_dim(ncidspec_kt, 'time',NF_UNLIMITED,dimid(4)) )         
     
! put global attributes

         call check(nf_put_att_int(ncidspec_kt,NF_GLOBAL,'nx',NF_INT,1,(LLm)))
         call check(nf_put_att_int(ncidspec_kt,NF_GLOBAL,'ny',NF_INT,1,(MMm)))        
         call check(nf_put_att_double(ncidspec_kt,NF_GLOBAL,'dx',
     &              NF_DOUBLE,1,1./pm(1,1)))  
         call check(nf_put_att_double(ncidspec_kt,NF_GLOBAL,'dy',
     &              NF_DOUBLE,1,1./pn(1,1)))
         call check(nf_put_att_double(ncidspec_kt,NF_GLOBAL,'dz',
     &              NF_DOUBLE,1,Hz(1,1,1)))
         call check(nf_put_att_double(ncidspec_kt,NF_GLOBAL,'levels',
     &              NF_DOUBLE,N_kt,depth_kt(1:N_kt)))

! define variables

         call check(nf_def_var(ncidspec_kt,'time' ,
     &                    nf_double,1,(/dimid(4)/),timeid))
         do l=1,nhat_kt
            name=trim(name_hat_kt(l))//'_real'            
            call check(nf_def_var(ncidspec_kt,trim(name),NF_DOUBLE,4,  
     &            (/dimid(1),dimid(2),dimid(3),dimid(4)/),varid(l))) 
            name=trim(name_hat_kt(l))//'_imag'            
            call check(nf_def_var(ncidspec_kt,trim(name),NF_DOUBLE,4,  
     &            (/dimid(1),dimid(2),dimid(3),dimid(4)/),varid(l+nhat_kt)))                       
         enddo        
    
         call check(nf_enddef(ncidspec_kt) )
          
      endif 
              
      case(1) ! write spec

! each processor write it's own cospec in the parallel netcdf file if needed
          
          nrecspec_kt=nrecspec_kt+1

          hat_kt = hat_kt / float(nwrtspec_kt/navgspec_kt)      

          if (mynode==0)  call check(nf_put_vara_double(ncidspec_kt,timeid,  
     &                     (/nrecspec_kt/),(/1/),time*sec2day) ) ! time in days      

          ! reads/writes to happen together on every processor.
          do l=1,nhat_kt

             if (mynode==0) then
                   Nin=1
                   kcount=nint 

                   call check(nf_put_vara_double(ncidspec_kt, varid(l),  
     &                  (/1,1,Nin,nrecspec_kt/),  
     &                  (/LLmspec,MMmspec,kcount,1/),
     &                  real(hat_kt(:,:,1:kcount,l))) )  
                   call check(nf_put_vara_double(ncidspec_kt, varid(l+nhat_kt),  
     &                  (/1,1,Nin,nrecspec_kt/),  
     &                  (/LLmspec,MMmspec,kcount,1/),
     &                  imag(hat_kt(:,:,1:kcount,l))) )  
             else
                    Nin=mynode*nint+1
                    if (Nin <= N_kt) then
                       call MPI_SEND(hat_kt,LLmspec*MMmspec*nint*nhat_kt,                
     &                       MPI_DOUBLE_COMPLEX,0,mynode,  
     &                       MPI_COMM_WORLD,mpierr)  
                    endif
             endif                     
               
             if (mynode==0) then
                  do iproc=1,NNODES-1
                     Nin=iproc*nint+1           
                     if (Nin <= N_kt) then                    
                        call MPI_IRECV(hat_kt,  
     &                      LLmspec*MMmspec*nint*nhat_kt,  
     &                      MPI_DOUBLE_COMPLEX,iproc,iproc,  
     &                      MPI_COMM_WORLD,reqid(iproc),mpierr)
                        call MPI_WAIT(reqid(iproc),istat,mpierr)
                        kcount = min(nint,N_kt-Nin+1)  
 	                                             
                        call check(nf_put_vara_double(ncidspec_kt, 
     &                     varid(l),  
     &                     (/1,1,Nin,nrecspec_kt/),  
     &                     (/LLmspec,MMmspec,kcount,1/),
     &                     real(hat_kt(1:LLmspec,1:MMmspec,1:kcount,l)))  )  
                        call check(nf_put_vara_double(ncidspec_kt, 
     &                     varid(l+nhat_kt),  
     &                     (/1,1,Nin,nrecspec_kt/),  
     &                     (/LLmspec,MMmspec,kcount,1/),
     &                     imag(hat_kt(1:LLmspec,1:MMmspec,1:kcount,l)))  )     		   
                     endif
                  enddo
               endif
            enddo
                              
            if (nrecspec_kt==nbtime) then
               ierr = nf_CLOSE(ncidspec_kt)
            else
               ierr = nf_sync(ncidspec_kt)
            endif
            hat_kt = 0. 
                 
      end select

      
      end subroutine wrt_diags_spec_kt 
      
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

      subroutine exchange_diags_spec_kt()    
# include "param.h"
# include "scalars.h"
# include "ocean3d.h"
# include "diagnostics.h"

      include 'mpif.h'


      !local
      integer::i,j,tile
      real, dimension(GLOBAL_2D_ARRAY,1:N) :: workr  
      
# include "compute_tile_bounds.h"       

      call exchange_u3d_tile (Istr,Iend,Jstr,Jend,
     &                        MV(START_2D_ARRAY,1,1))    
      call exchange_v3d_tile (Istr,Iend,Jstr,Jend,
     &                        MV(START_2D_ARRAY,1,2))    
      do itrc=1,NT
        call exchange_r3d_tile (Istr,Iend,Jstr,Jend,
     &                          TT(START_2D_ARRAY,1,itrc))            
        call exchange_r3d_tile (Istr,Iend,Jstr,Jend,
     &                          TTmoy(START_2D_ARRAY,1,itrc))
      enddo   
       
      end subroutine exchange_diags_spec_kt



#endif  /* DIAG_SPEC_KT */  
    
      end module diagspec_kt