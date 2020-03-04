! 
!  for each level k
!     fst_F_t(i,j,k,t) = fst_FT * fst_W(i,j) * fst_out(i,j)
!  
! Doc about fftw:
!  http://www.fftw.org/fftw3_doc/Multi_002ddimensional-Transforms.html#Multi_002ddimensional-Transforms  
!
! ... further details about array sizes in spectral complex space ...
!
! this lines activates copy of the present file to the backup folder: aponte

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
     &                               fst_a = (/1.0,0.0,0.0,0.0,0.0/)

! For the real form of the forcing
      real(kind=8), parameter :: fst_xmid  = 512.e+3, fst_xwid = 200.e+3
      real(kind=8), parameter :: fst_ymid  = 1440.e+3, fst_ywid = 500.e+3
!     amplitude of the forcing
!      real(kind=8), parameter :: fst_FT = 1.e-1
      real(kind=8) :: fst_FT
!     gaussian form of the forcing
      real(kind=8),dimension(:,:), allocatable :: fst_W
!     total forcing
      real(kind=8),dimension(:,:,:,:), allocatable :: fst_F_t
      real(kind=8),dimension(:,:,:), allocatable :: fst_F_zeta

! For the exchanges
      integer   :: nint
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
      real(kind=8) :: kx, ky, uF

      ! forcing wavenumbers
      fst_kF = 16./xl
      fst_dkF = fst_kF/4.
      ! forcing amplitude
      !uF = 0.1
      !fst_FT = uF*1.0313E-4 /fst_kF/g *rho0/abs(Tcoef) /fst_dt
      !  ~ 1e-4
      ! above is not pretty: too many hardcoded parameters:
      ! f0 = 1.0313E-4
      fst_FT = 0.1

      ! nint = number of vertical levels per procs
      if ( mod(N,NNODES) == 0 ) then
         nint=N/NNODES
      else
         nint=N/NNODES + 1
      endif

! allocation dynamique des tableaux

      L2 = LLm/2
      M2 = MMm/2
      allocate(fst_W(GLOBAL_2D_ARRAY))
      allocate(fst_F_t(GLOBAL_2D_ARRAY,N,2))
      allocate(fst_F_zeta(GLOBAL_2D_ARRAY,2))
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
            fst_psi(i,j) = cmplx(exp(-(sqrt(kx**2 + ky**2) - fst_kF)**2
     &                              /(fst_dkF**2)),0.)
        enddo
      enddo
      ! normalize fst_psi
      fst_psi = fst_psi/sqrt(sum(abs(fst_psi)**2)) /1.42
      ! 1.42 = std of a (-1,1) uniform random variable
      !
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

      if (mynode==0) print*,"compute turbulence forcing at time=", curtime
      
      ! loop around vertical levels for this proc, fill fst_F_t
      do kint=1,nint
        call spec_to_phys(kint, curtime)
      enddo
      ! for zeta, fill fst_F_zeta
      call spec_to_phys(0, curtime)

      if ( curtime <= fst_tdecay ) then
        decay = 1.
      else
        decay = 0.
      endif

      do k=1,N
         fst_F_t(:,:,k,nfstnew) = fst_F_t(:,:,k,nfstnew) * fst_FT 
     &                          * fst_W(:,:)
     &                          * decay
      enddo

      fst_F_zeta(:,:,nfstnew) = fst_F_zeta(:,:,nfstnew) * fst_FT 
     &                          * fst_W(:,:)
     &                          * decay

      end subroutine comp_fsturb

!----------------------------------------------------------------------

      subroutine spec_to_phys(kint, curtime)

# include "param.h" 
# include "scalars.h"

      integer :: kv,kint   
      !real(kind=8) :: rho_mean, curtime
      real(kind=8) :: curtime
              
      ! compute complex amplitude for variable kv and level kint

      include 'mpif.h'
      include 'fftw3.f'

      !local
      integer::iproc,i,j,k
      !real :: mean

      ! computes the true vertical level k corresponding to kint
      if ( kint > 0 ) then
        k = (mynode*nint)+kint
      else
        k = 0
      endif
      if (k <= N) then
        call compute_spectral(k, curtime)
      endif
      
      call launch_rcv_fst(kint)
      if (k <= N) then
        do iproc=0,NNODES-1
            call launch_send_fst(iproc, k)
        enddo
      endif
      call launch_wait_fst(kint)

      end subroutine spec_to_phys

!----------------------------------------------------------------------

      subroutine compute_spectral(k, curtime)
  !
  # include "param.h"
  # include "scalars.h"

      ! k
      integer :: k
      real(kind=8) :: curtime

      ! local
      integer :: nbseed
      integer, allocatable :: seed(:)
      integer :: i, j, m
      real(kind=8),dimension(L2+1,MMm) :: rNoise, iNoise
      real(kind=8),dimension(L2+1,MMm) :: Ftmp

      fst_in = cmplx(0.,0.)

      !  Get the random noise Noisefst

      ! restart the state of the pseudorandom number generator
      call random_seed(size=nbseed)
      allocate(seed(nbseed))

      do m=1,fst_Nmode

         ! compute N
         ! assign information in values to seed
         seed = int(curtime*m)
         ! seed the random number generator
         call random_seed(put=seed)
         call random_number(rNoise)
         seed = int(curtime*m+1)
         ! seed the random number generator
         call random_seed(put=seed)
         call random_number(iNoise)
         
  !         fst_in = cmplx(rNoisefst*2.-1.,iNoisefst*2.-1.)*fst_psi
         Ftmp = cmplx(rNoise*2.-1.,iNoise*2.-1.) 
     &                     *fst_a(k)*fst_psi
         ! reset 0 and Nyquist frequencies to real number
         Ftmp(1,1) = real(Ftmp(1,1))
         Ftmp(L2+1,MMm) = real(Ftmp(L2+1,MMm))

         if (k>0) then
           ! fst_F_t
           fst_in = fst_in + Ftmp*phir(k,m)
         elseif (mynode==0) then
           ! fst_F_zeta
           fst_in = fst_in - Ftmp*phiuv(N,m)/g/(-Tcoef)
         endif
      enddo

      if (k==N) then
        call var2cdf(real(fst_in),size(shape(real(fst_in))),
     &               shape(real(fst_in)),'fst_in.nc')
      endif

  ! inverse fourier transform
      call dfftw_execute(plan)

      if (k==N) then
        call var2cdf(fst_out,size(shape(fst_out)),shape(fst_out),
     &              'fst_out.nc')
      endif

      end subroutine compute_spectral

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
      if ( kint>0 ) then
        do iproc=0,NNODES-1
             k = (iproc*nint)+kint
             if( mynode.ne.iproc .and. k <= N) then
                call MPI_IRECV(zexchang(:,:,iproc),npts,  
     &            MPI_DOUBLE_PRECISION,iproc,iproc,MPI_COMM_WORLD,
     &            reqid(iproc),mpierr)
             endif
        enddo
      elseif ( mynode>0 ) then
        iproc = 0
        call MPI_IRECV(zexchang(:,:,iproc),npts,  
     &            MPI_DOUBLE_PRECISION,iproc,iproc,MPI_COMM_WORLD,
     &            reqid(iproc),mpierr)
      endif
      
      end subroutine launch_rcv_fst

!----------------------------------------------------------------------

      subroutine launch_send_fst(iproc,k)

      integer::iproc,k
      
# include "param.h"
# include "scalars.h"
# include "mpi_cpl.h"

      include 'mpif.h'

      !local
      integer::i,j,i0,j0,npts,mpierr,i_offset,j_offset
      
      npts = Lm*Mm
              
      if ( k>0 .and. mynode==iproc ) then
       ! copy directly in the array 
       ! fst_F_t
       do j=1,Mm
          do i=1,Lm
             i_offset=mod(iproc,NP_XI)*Lm
             j_offset=(iproc/NP_XI)*Mm 
             i0 = i + i_offset
             j0 = j + j_offset
             fst_F_t(i,j,k,nfstnew) = fst_out(i0,j0)
          enddo
       enddo            
      elseif ( k>0 .or. (k==0 .and. mynode==0) ) then
       ! send to iproc
       ! fst_F_t and fst_F_zeta
       do j=1,Mm
          do i=1,Lm
             i_offset=mod(iproc,NP_XI)*Lm
             j_offset=(iproc/NP_XI)*Mm 
             i0 = i + i_offset
             j0 = j + j_offset
             zpar(i,j) = fst_out(i0,j0)
          enddo
       enddo            
       call MPI_SEND(zpar,npts,MPI_DOUBLE_PRECISION,iproc,mynode,  
      &                  MPI_COMM_WORLD,mpierr)
      endif

      end subroutine launch_send_fst

!----------------------------------------------------------------------

      subroutine launch_wait_fst(kint)
              
      integer :: kint

# include "param.h"
# include "scalars.h"

      include 'mpif.h'

      integer::iproc,i0,j0,i,j,k,mpierr,i_offset,j_offset
      integer,dimension(MPI_STATUS_SIZE):: istat
      
!     launch all the wait
      if ( kint>0 ) then
        do iproc=0,NNODES-1
           k=(iproc*nint)+kint
           if (mynode.ne.iproc .and. k <= N) then
             call MPI_WAIT(reqid(iproc),istat,mpierr)
             do j=1,Mm
               do i=1,Lm
                 fst_F_t(i,j,k,nfstnew) = zexchang(i,j,iproc)  
               enddo
             enddo 
           endif
        enddo
      elseif ( mynode>0 ) then
        iproc = 0
        call MPI_WAIT(reqid(iproc),istat,mpierr)
        do j=1,Mm
          do i=1,Lm
            fst_F_zeta(i,j,nfstnew) = zexchang(i,j,iproc)  
          enddo
        enddo         
      endif
      
      end subroutine launch_wait_fst

#endif /* FSTURB */

      end module fsturb
