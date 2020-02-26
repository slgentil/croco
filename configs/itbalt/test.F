! Ce programme calcule la transformée de Fourier en deux dimensions d'un produit de cosinus. Il s'agit de la transformée réelle
!
      Program test
      implicit none

        include 'fftw3.f'

      real, parameter               :: xl=1024e+3, el=2880e+3  
      integer,parameter             :: L=256, M=720        
C       real, parameter               :: xl=1024e+3, el=xl  
C       integer,parameter             :: L=256, M=L        
      integer,parameter             :: Lspec=L/2+1, M2=(M/2)*2+1
      real(kind=8)                  :: fst_kF, fst_dkF, fst_ymid ,fst_ywid
      integer                       :: i, j
      integer*8                     :: plan

      real :: dx,dy
      real(kind=8),dimension(Lspec) :: fst_kx
      real(kind=8),dimension(M2) :: fst_ky
      real(kind=8)                  :: kx, ky
!
      integer :: nbseed
      integer, allocatable :: seed(:)
      real(kind=8), parameter :: fst_a=0.1, fst_FT = 1.
      real(kind=8),dimension(L,M) :: fst_W, yr
      real(kind=8),dimension(Lspec,M) :: fst_rNoise
      real(kind=8),dimension(Lspec,M) :: fst_iNoise
      complex(kind=8) :: cdata(Lspec,M)
      real(kind=8) :: rdata(L,M)
      real(kind=8) :: dummy

! Création des plans pour la transformée de Fourier
      call dfftw_plan_dft_c2r_2d(plan,L,M,cdata,rdata,
     &                           FFTW_ESTIMATE)       


      fst_ymid=real(M)/2.
      fst_ywid=real(M)/2.
      do i=1,L
        yr(i,:) = (/ (j, j=1,M) /)
      enddo
      do j=1,M
        do i=1,L
          fst_W(i,j) = exp(-((yr(i,j)-fst_ymid) / (0.25*fst_ywid))**2)
        enddo
      enddo
      call var2cdf(fst_W,size(shape(fst_W)),shape(fst_W),'fst_W.nc')

C       kFfst=L/4.
C       dkFfst = L/16.
C       kFfst=xl/16.
C       dkFfst = xl/16.
      fst_kF=xl/16.
      fst_dkF = xl/16.

!     Make wave number mask
      dx = 1./xl
      dy = 1./el
      fst_kx = (/ (i*dx, i= 0,Lspec-1) /)               
      fst_ky = (/ (j*dy, j=1,M2/2),(j*dy, j=-M2/2,0) /)   

      do j=1,M
        ky = fst_ky(j)
        do i=1,Lspec
            kx = fst_kx(i)
            cdata(i,j) = CMPLX(exp(-(kx**2 + ky**2 - fst_kF**2)
     &                             /(fst_dkF**2)),0.)
        enddo
      enddo

C       cdata = CMPLX(0.,0.)
C       cdata (Lspec/2,M/2)= CMPLX(1.,0.)
      call var2cdf(real(cdata),size(shape(real(cdata))),
     &             shape(real(cdata)),'fst_cdata.nc')

      ! restart the state of the pseudorandom number generator
      call random_seed(size=nbseed)
      allocate(seed(nbseed))
      ! assign information in values to seed
      seed = 10
      ! seed the random number generator
      call random_seed(put=seed)
      call random_number(fst_rNoise)
      ! assign information in values to seed
      seed = 11
      ! seed the random number generator
      call random_seed(put=seed)
      call random_number(fst_iNoise)
      call var2cdf(fst_rNoise,size(shape(fst_rNoise)),shape(fst_rNoise),
     &               'fst_Noise.nc')

      cdata = cdata * cmplx(fst_rNoise,fst_iNoise) * fst_a

C       call var2cdf(real(cdata),size(shape(real(cdata))),
C      &             shape(real(cdata)),'cdata.nc')
!
!  Appel des fonctions de transformation
!
      call dfftw_execute(plan)
      rdata = rdata/L/M

      call var2cdf(rdata,size(shape(rdata)),shape(rdata),'fst_rdata1.nc')

      rdata = rdata * fst_W * fst_FT

      call var2cdf(rdata,size(shape(rdata)),shape(rdata),'fst_rdata2.nc')
!
! Destruction des plans
!
      call dfftw_destroy_plan(plan)
      end
