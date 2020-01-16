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
      real(kind=8)                  :: kFfst, dkFfst, yfstmid ,yfstwid
      integer                       :: i, j
      integer*8                     :: plan

      real :: dx,dy
      real(kind=8),dimension(Lspec) :: kxfst
      real(kind=8),dimension(M2) :: kyfst
      real(kind=8)                  :: kx, ky
!
      integer :: nbseed
      integer, allocatable :: seed(:)
      real(kind=8), parameter :: afst=0.1, FTfst = 1.
      real(kind=8),dimension(L,M) :: Wfst, yr
      real(kind=8),dimension(Lspec,M) :: rNoisefst
      real(kind=8),dimension(Lspec,M) :: iNoisefst
      complex(kind=8) :: cdata(Lspec,M)
      real(kind=8) :: rdata(L,M)
      real(kind=8) :: dummy

! Création des plans pour la transformée de Fourier
      call dfftw_plan_dft_c2r_2d(plan,L,M,cdata,rdata,
     &                           FFTW_ESTIMATE)       


      yfstmid=real(M)/2.
      yfstwid=real(M)/2.
      do i=1,L
        yr(i,:) = (/ (j, j=1,M) /)
      enddo
      do j=1,M
        do i=1,L
          Wfst(i,j) = exp(-((yr(i,j)-yfstmid) / (0.25*yfstwid))**2)
        enddo
      enddo
      call var2cdf(Wfst,size(shape(Wfst)),shape(Wfst),'Wfst.nc')

C       kFfst=L/4.
C       dkFfst = L/16.
C       kFfst=xl/16.
C       dkFfst = xl/16.
      kFfst=xl/16.
      dkFfst = xl/16.

!     Make wave number mask
      dx = xl/(L)
      dy = el/(M)
      kxfst = (/ (i*dx, i= 0,Lspec-1) /)               
      kyfst = (/ (j*dy, j=1,M2/2),(j*dy, j=-M2/2,0) /)   

      do j=1,M
        ky = kyfst(j)
        do i=1,Lspec
            kx = kxfst(i)
            cdata(i,j) = CMPLX(exp(-(kx**2 + ky**2 - kFfst**2)
     &                             /(dkFfst**2)),0.)
        enddo
      enddo

C       cdata = CMPLX(0.,0.)
C       cdata (Lspec/2,M/2)= CMPLX(1.,0.)
      call var2cdf(real(cdata),size(shape(real(cdata))),
     &             shape(real(cdata)),'cdata.nc')

      ! restart the state of the pseudorandom number generator
      call random_seed(size=nbseed)
      allocate(seed(nbseed))
      ! assign information in values to seed
      seed = 10
      ! seed the random number generator
      call random_seed(put=seed)
      call random_number(rNoisefst)
      ! assign information in values to seed
      seed = 11
      ! seed the random number generator
      call random_seed(put=seed)
      call random_number(iNoisefst)
      call var2cdf(rNoisefst,size(shape(rNoisefst)),shape(rNoisefst),
     &               'Noisefst.nc')

      cdata = cdata * cmplx(rNoisefst,iNoisefst) * afst

C       call var2cdf(real(cdata),size(shape(real(cdata))),
C      &             shape(real(cdata)),'cdata.nc')
!
!  Appel des fonctions de transformation
!
      call dfftw_execute(plan)
      rdata = rdata/L/M

      call var2cdf(rdata,size(shape(rdata)),shape(rdata),'rdata1.nc')

      rdata = rdata * Wfst * FTfst

      call var2cdf(rdata,size(shape(rdata)),shape(rdata),'rdata2.nc')
!
! Destruction des plans
!
      call dfftw_destroy_plan(plan)
      end