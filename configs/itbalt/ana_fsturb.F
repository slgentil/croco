! 
!  for each level k
!     Ffst(i,j,k,t) = FTfst * Wfst(i,j) * outfst(i,j)
!  
!  
!
!

#include "cppdefs.h"
       
      module fsturb

#ifdef FSTURB

      implicit none

! For the time stepping
!     indices of the old and new fsturb time step
      integer :: nfstold, nfstnew
!     time step of the computation of the forcing (second)
      real(kind=8), parameter :: dtfst = 10.*86400.
!     time after which forcing is killed
      real(kind=8), parameter :: fst_tdecay = 900.*86400.
!     time of computation corresponding to the nfstnew indice
      real(kind=8)            :: timefst

! Number of modes to keep, amplitude of each mode
      integer, parameter      :: Nmodefst=5
      real(kind=8), dimension(Nmodefst), parameter :: 
     &                               afst = (/0.1,0.1,0.1,0.1,0.1/)

! For the real form of the forcing
      real(kind=8), parameter :: xfstmid  = 512.e+3, xfstwid = 200.e+3
      real(kind=8), parameter :: yfstmid  = 1440.e+3, yfstwid = 500.e+3
!     amplitude of the forcing
      real(kind=8), parameter :: FTfst = 1.e3
!     gaussian form of the forcing
      real(kind=8),dimension(:,:), allocatable :: Wfst
!     total forcing
      real(kind=8),dimension(:,:,:,:), allocatable :: Ffst

! For the exchanges
      integer   :: nint, ktot
      integer*8,dimension(:), allocatable        ::reqid
      real(kind=8),dimension(:,:), allocatable   :: zpar
      real(kind=8),dimension(:,:,:), allocatable :: zexchang

! For the complex part of the forcing
      integer   :: Lspec, M2
!     Radius and width for the horizontal part of the spectral initialization      
      real(kind=8)            :: kFfst, dkFfst
!     wave length   
      real(kind=8),dimension(:), allocatable     ::kxfst
      real(kind=8),dimension(:), allocatable     :: kyfst
!     horizontal part of the spectral initialization   
      complex(kind=8),dimension(:,:), allocatable   :: psifst

! For the FFT
      integer*8 :: plan
      complex(kind=8),dimension(:,:), allocatable   :: infst 
      real(kind=8),dimension(:,:), allocatable   :: outfst

      contains

!----------------------------------------------------------------------

      subroutine get_kxky
      
# include "param.h"
# include "scalars.h"
      
      integer :: i,j
      real :: dx,dy

      dx = xl/(LLm)
      dy = el/(MMm)
      kxfst = (/ (i*dx, i= 0,Lspec-1) /)               
      kyfst = (/ (j*dy, j=1,M2/2),(j*dy, j=-M2/2,0) /)  
      
      end subroutine get_kxky

!----------------------------------------------------------------------

      subroutine init_fsturb

      use params

# include "param.h"
# include "scalars.h"
# include "grid.h"
        include 'fftw3.f'

!     include 'mpif.h'

!  This routine init the forcing for turbulence


      character(len=32) :: my_fmt
      integer i,j,k
      real(kind=8) :: kx, ky
!
      kFfst=xl/16.
!      dkFfst = xl/16.
      dkFfst = xl/32.

      if ( mod(N,NNODES) == 0 ) then
         nint=N/NNODES
      else
         nint=N/NNODES + 1
      endif

! allocation dynamique des tableaux

      Lspec = LLm/2+1
      M2=(MMm/2)*2+1
      allocate(Wfst(GLOBAL_2D_ARRAY))
      allocate(Ffst(GLOBAL_2D_ARRAY,N,2))
      allocate(reqid(0:NNODES-1))
      allocate(zpar(Lm,Mm))
      allocate(zexchang(Lm,Mm,0:NNODES-1))
      allocate(kxfst(Lspec))
      allocate(kyfst(M2))
      allocate(infst(Lspec,MMm))
      allocate(outfst(LLm,MMm))
      allocate(psifst(Lspec,MMm))

! Make wave number mask
      call get_kxky  
!
! Initialization FFT complex to real, out of place:
!
      call dfftw_plan_dft_c2r_2d(plan,LLm,MMm,infst,outfst,
     &                           FFTW_ESTIMATE)   


     ! Init the form of the forcing
     ! xfstmid,xfstwid,yfstmid,yfstwid initialized as parameters
      do j=1,Mm
        do i=1,Lm
          Wfst(i,j) = exp(-((yr(i,j)-yfstmid) / (0.25*yfstwid))**2)
!    &              * exp(-((xr(i,j)-xfstmid) / (0.25*xfstwid))**2)
        enddo
      enddo
      if (mynode==0) call var2cdf(real(Wfst),size(shape(real(Wfst))), 
     &                 shape(real(Wfst)),'Wfst.nc')
  
      ! compute psi
      do j=1,MMm
        ky = kyfst(j)
        do i=1,Lspec
            kx = kxfst(i)
            psifst(i,j) = cmplx(exp(-(kx**2 + ky**2 - kFfst**2)
     &                              /(dkFfst**2)),0.)
        enddo
      enddo
      if (mynode==0) call var2cdf(real(psifst),size(shape(real(psifst))), 
     &                 shape(real(psifst)),'psifst.nc')

      timefst = time
      nfstold = 1
      nfstnew = 1
      call comp_fsturb(time)
      nfstnew = 2
      call comp_fsturb(time+dtfst)

!  Output parameters

      call params_output('fst_dt',dtfst)
      call params_output('fst_tdecay',fst_tdecay)
      call params_output('fst_Nmode',Nmodefst)
      call params_output('fst_a',afst)
      call params_output('fst_xmid',xfstmid)
      call params_output('fst_xwid',xfstwid)
      call params_output('fst_ymid',yfstmid)
      call params_output('fst_ywid',yfstwid)
      call params_output('fst_FT',FTfst)
      call params_output('fst_k',kFfst)
      call params_output('fst_dk',dkFfst)

      end subroutine init_fsturb

!----------------------------------------------------------------------

      subroutine comp_fsturb(curtime)
      
       

# include "param.h"
# include "scalars.h"

      include 'mpif.h'

      !local
      integer :: ik,k,kint
      real(kind=8) :: curtime
      real(kind=8) :: decay

      if (mynode==0) print*,"compute turbulence forcing at time=",curtime
      do kint=1,nint  
        call spec_to_phys(kint,curtime)
      enddo

      if ( curtime <= fst_tdecay ) then
        decay = 1.
      else
        decay = 0.
      endif

      do k=1,N 
         Ffst(:,:,k,nfstnew) = Ffst(:,:,k,nfstnew) * FTfst * Wfst(:,:)
     &                          * decay
      enddo

      end subroutine comp_fsturb

!----------------------------------------------------------------------

      subroutine spec_to_phys(kint,curtime)
 
# include "param.h" 
# include "scalars.h"
     
      integer::kv,kint   
      real(kind=8) :: rho_mean, curtime
        
      
      ! compute complex amplitude for variable kv and level kint

      include 'mpif.h'
      include 'fftw3.f'

      !local
      integer::iproc,i,j,k
      real :: mean

      k=(mynode*nint)+kint
      if (k <= N) then
        call compute_spectral(k,curtime)
      endif
      call launch_rcv_fst(kint)
      if (k <= N) then
        do iproc=0,NNODES-1    
            call launch_send_fst(iproc,k)
        enddo
      endif
      call launch_wait_fst(kint)  


      end subroutine spec_to_phys

!----------------------------------------------------------------------

      subroutine launch_rcv_fst(kint)
                 
# include "param.h"
# include "scalars.h"
# include "mpi_cpl.h"

      include 'mpif.h'

    

      integer::iproc,i0,j0,npts,mpierr,k,kint
      
!     launch all the nonblocking receive requests
      npts = Lm*Mm
      zexchang = 0.
      do iproc=0,NNODES-1
         k=(iproc*nint)+kint
         if(mynode.ne.iproc .and. k <= N) then
            call MPI_IRECV(zexchang(:,:,iproc),npts,  
     &            MPI_DOUBLE_PRECISION,iproc,iproc,MPI_COMM_WORLD,
     &            reqid(iproc),mpierr)
         endif
      enddo

      
      end subroutine launch_rcv_fst

!----------------------------------------------------------------------

      subroutine launch_wait_fst(kint)
              
      integer :: kint

# include "param.h"
# include "scalars.h"

      include 'mpif.h'


      integer::iproc,i0,j0,i,j,k,mpierr,i_offset,j_offset
      integer,dimension(MPI_STATUS_SIZE):: istat
      
!     launch all the wait
      do iproc=0,NNODES-1
         k=(iproc*nint)+kint
         if (mynode.ne.iproc .and. k <= N) then
           call MPI_WAIT(reqid(iproc),istat,mpierr)
           do j=1,Mm
             do i=1,Lm
               Ffst(i,j,k,nfstnew) = zexchang(i,j,iproc)  
             enddo
           enddo 
         endif
      enddo
      
      end subroutine launch_wait_fst

!----------------------------------------------------------------------

      subroutine launch_send_fst(iproc,k)

      integer::iproc,k
      
# include "param.h"
# include "scalars.h"
# include "mpi_cpl.h"

      include 'mpif.h'


      !local
      integer::i,j,i0,j0,npts,mpierr,i_offset,j_offset
      
      if (mynode.eq.iproc) then
         ! copy directly in the array 
         do j=1,Mm
            do i=1,Lm
               i_offset=mod(iproc,NP_XI)*Lm
               j_offset=(iproc/NP_XI)*Mm 
               i0 = i + i_offset
               j0 = j + j_offset
               Ffst(i,j,k,nfstnew)=outfst(i0,j0)
            enddo
         enddo            
      else
         ! send to iproc
         npts = Lm*Mm
         do j=1,Mm
            do i=1,Lm
               i_offset=mod(iproc,NP_XI)*Lm
               j_offset=(iproc/NP_XI)*Mm 
               i0 = i + i_offset
               j0 = j + j_offset
               zpar(i,j)=outfst(i0,j0)
            enddo
         enddo            
         call MPI_SEND(zpar,npts,MPI_DOUBLE_PRECISION,iproc,mynode,  
     &                  MPI_COMM_WORLD,mpierr)
      endif
     

      end subroutine launch_send_fst

!----------------------------------------------------------------------

      subroutine compute_spectral(ktot,curtime)
!
# include "param.h"
# include "scalars.h"

      integer :: ktot
      real(kind=8) :: curtime

      ! local
      integer :: nbseed
      integer, allocatable :: seed(:)
      integer :: i, j, k
      real(kind=8),dimension(Lspec,MMm) :: rNoisefst, iNoisefst

      infst = cmplx(0.,0.)

!  Get the random noise Noisefst

      ! restart the state of the pseudorandom number generator
      call random_seed(size=nbseed)
      allocate(seed(nbseed))


      do k=1,Nmodefst

         ! compute N
         ! assign information in values to seed
         seed = int(curtime*k)
         ! seed the random number generator
         call random_seed(put=seed)
         call random_number(rNoisefst)
         seed = int(curtime*k+1)
         ! seed the random number generator
         call random_seed(put=seed)
         call random_number(iNoisefst)


         infst = infst + cmplx(rNoisefst,iNoisefst)* afst(k)*psifst
     &                   *phir(ktot,k)

      enddo

      if (ktot==N) then
        call var2cdf(real(infst),size(shape(real(infst))),
     &               shape(real(infst)),'infst.nc')
      endif

! fourierize
      call dfftw_execute(plan)
      outfst = outfst / LLm / MMm

      if (ktot==N) then
        call var2cdf(outfst,size(shape(outfst)),shape(outfst),
     &              'outfst.nc')
      endif

      end subroutine compute_spectral

#endif /* FSTURB */

      end module fsturb
