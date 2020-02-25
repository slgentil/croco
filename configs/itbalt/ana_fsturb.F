! 
!  for each level k
!     fst_F(i,j,k,t) = fst_FT * fst_W(i,j) * fst_out(i,j)
!  
! Doc about fftw:
!  http://www.fftw.org/fftw3_doc/Multi_002ddimensional-Transforms.html#Multi_002ddimensional-Transforms  
!
! ... further details about array sizes in spectral complex space ...
!

#include "cppdefs.h"
       
      module fsturb

#ifdef FSTURB

      implicit none

! For the time stepping
!     indices of the old and new fsturb time step
      integer :: nfstold, nfstnew
!     time step of the computation of the forcing (second)
      real(kind=8), parameter :: fst_dt = 10.*86400.
!     time after which forcing is killed
      real(kind=8), parameter :: fst_tdecay = 900.*86400.
!     time of computation corresponding to the nfstnew indice
      real(kind=8)            :: fst_time

! Number of modes to keep, amplitude of each mode
      integer, parameter      :: fst_Nmode=5
      real(kind=8), dimension(fst_Nmode), parameter :: 
     &                               fst_a = (/0.0,0.1,0.0,0.0,0.0/)

! For the real form of the forcing
      real(kind=8), parameter :: fst_xmid  = 512.e+3, fst_xwid = 200.e+3
      real(kind=8), parameter :: fst_ymid  = 1440.e+3, fst_ywid = 500.e+3
!     amplitude of the forcing
      real(kind=8), parameter :: fst_FT = 1.e-1
!     gaussian form of the forcing
      real(kind=8),dimension(:,:), allocatable :: fst_W
!     total forcing
      real(kind=8),dimension(:,:,:,:), allocatable :: fst_F

! For the exchanges
      integer   :: nint, ktot
      integer*8,dimension(:), allocatable        ::reqid
      real(kind=8),dimension(:,:), allocatable   :: zpar
      real(kind=8),dimension(:,:,:), allocatable :: zexchang

! For the complex part of the forcing
      integer   :: L2, M2
!     Radius and width for the horizontal part of the spectral initialization      
      real(kind=8)            :: fst_kF, fst_dkF
!     wave length   
      real(kind=8),dimension(:), allocatable     ::fst_kx
      real(kind=8),dimension(:), allocatable     :: fst_ky
!     horizontal part of the spectral initialization   
      complex(kind=8),dimension(:,:), allocatable   :: fst_psi

! For the FFT
      integer*8 :: plan
      complex(kind=8),dimension(:,:), allocatable   :: fst_in 
      real(kind=8),dimension(:,:), allocatable   :: fst_out

      contains

!----------------------------------------------------------------------

      subroutine get_kxky
      
# include "param.h"
# include "scalars.h"
      
      integer :: i,j
      real :: dkx,dky

      dkx = 1./xl
      dky = 1./el
      fst_kx = (/ (i*dkx, i=0,L2) /)
      fst_ky = (/ (j*dky, j=0,M2) , (j*dky, j=-M2+1,-1) /)
      
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
      fst_kF = 16./xl
      fst_dkF = fst_kF/2.
!      fst_dkF = xl/32.

      if ( mod(N,NNODES) == 0 ) then
         nint=N/NNODES
      else
         nint=N/NNODES + 1
      endif

! allocation dynamique des tableaux

      L2 = LLm/2
      M2 = MMm/2
      allocate(fst_W(GLOBAL_2D_ARRAY))
      allocate(fst_F(GLOBAL_2D_ARRAY,N,2))
      allocate(reqid(0:NNODES-1))
      allocate(zpar(Lm,Mm))
      allocate(zexchang(Lm,Mm,0:NNODES-1))
      allocate(fst_kx(L2+1))
      allocate(fst_ky(MMm))
      allocate(fst_in(L2+1,MMm))
      allocate(fst_out(LLm,MMm))
      allocate(fst_psi(L2+1,MMm))

! Make wave number mask
      call get_kxky  
!
! Initialization FFT complex to real, out of place:
!
      call dfftw_plan_dft_c2r_2d(plan,LLm,MMm,fst_in,fst_out,
     &                           FFTW_ESTIMATE)   


     ! Init the form of the forcing
     ! fst_xmid,fst_xwid,yfstmid,fst_ywid initialized as parameters
      do j=1,Mm
        do i=1,Lm
          fst_W(i,j) = exp(-((yr(i,j)-fst_ymid) / (0.25*fst_ywid))**2)
!    &              * exp(-((xr(i,j)-fst_xmid) / (0.25*fst_xwid))**2)
        enddo
      enddo
      if (mynode==0) call var2cdf(real(fst_W),size(shape(real(fst_W))), 
     &                 shape(real(fst_W)),'fst_W.nc')
  
      ! compute psi
      do j=1,MMm
        ky = fst_ky(j)
        do i=1,L2+1
            kx = fst_kx(i)
            fst_psi(i,j) = cmplx(exp(-(kx**2 + ky**2 - fst_kF**2)
     &                              /(fst_dkF**2)),0.)
        enddo
      enddo
      if (mynode==0) call var2cdf(real(fst_psi),size(shape(real(fst_psi))), 
     &                 shape(real(fst_psi)),'fst_psi.nc')

      fst_time = time
      nfstold = 1
      nfstnew = 1
      call comp_fsturb(time)
      nfstnew = 2
      call comp_fsturb(time+fst_dt)

!  Output parameters

      call params_output('fst_dt',fst_dt)
      call params_output('fst_tdecay',fst_tdecay)
      call params_output('fst_Nmode',fst_Nmode)
      call params_output('fst_a',fst_a)
      call params_output('fst_xmid',fst_xmid)
      call params_output('fst_xwid',fst_xwid)
      call params_output('fst_ymid',fst_ymid)
      call params_output('fst_ywid',fst_ywid)
      call params_output('fst_FT',fst_FT)
      call params_output('fst_k',fst_kF)
      call params_output('fst_dk',fst_dkF)

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
         fst_F(:,:,k,nfstnew) = fst_F(:,:,k,nfstnew) * fst_FT * fst_W(:,:)
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
               fst_F(i,j,k,nfstnew) = zexchang(i,j,iproc)  
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
               fst_F(i,j,k,nfstnew)=fst_out(i0,j0)
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
               zpar(i,j)=fst_out(i0,j0)
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
!      real(kind=8),dimension(Lspec,MMm) :: rNoisefst, iNoisefst
      real(kind=8),dimension(LLm,MMm) :: rNoisefst, iNoisefst

      fst_in = cmplx(0.,0.)

!  Get the random noise Noisefst

      ! restart the state of the pseudorandom number generator
      call random_seed(size=nbseed)
      allocate(seed(nbseed))


      do k=1,fst_Nmode

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


         fst_in = fst_in + cmplx(rNoisefst,iNoisefst)* fst_a(k)*fst_psi
     &                   *phir(ktot,k)

      enddo

      if (ktot==N) then
        call var2cdf(real(fst_in),size(shape(real(fst_in))),
     &               shape(real(fst_in)),'fst_in.nc')
      endif

! fourierize
      call dfftw_execute(plan)
      fst_out = fst_out / LLm / MMm

      if (ktot==N) then
        call var2cdf(fst_out,size(shape(fst_out)),shape(fst_out),
     &              'fst_out.nc')
      endif

      end subroutine compute_spectral

#endif /* FSTURB */

      end module fsturb
