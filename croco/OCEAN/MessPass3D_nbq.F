!
! $Id: MessPass3D.F 1143 2013-05-17 08:17:23Z serena $
!
#include "cppdefs.h"
#if defined MPI && defined SOLVE3D && defined NBQ
# ifndef MP_3PTS
      subroutine MessPass3D_tile_init_nbq (iv,Istr,Iend,Jstr,Jend,
     &	A, kst, nmax)
# else
      subroutine MessPass3D_3pts_tile_init_nbq (iv,Istr,Iend,Jstr,Jend,
     &	A,kst, nmax)
# endif
!
!======================================================================
!
! This subroutine is designed for ROMS-MPI code. It exchanges domain
! boundary information, including 2 (or 3) ghost-cells in each
! direction.
!
! Ping Wang 9/15/99.
! Patrick Marchesiello 2012: generic number of ghost-cells Npts
!
!======================================================================
!
!     implicit none
# include "param.h"
# include "scalars.h"
# include "mpi_cpl.h"
      include 'mpif.h'
!
! Nb of boundary points involved in communication
!
      integer Npts,ipts,jpts
# ifndef MP_3PTS
      parameter (Npts=2)
# else
      parameter (Npts=3)
# endif
      integer :: maxvar
      parameter(maxvar=10)
      integer kst, nmax
      real A(GLOBAL_2D_ARRAY,kst:nmax)
CSDISTRIBUTE_RESHAPE A(BLOCK_PATTERN) BLOCK_CLAUSE
      integer Istr,Iend,Jstr,Jend, i,j,k, isize,jsize,ksize,
     &        req(18), status(MPI_STATUS_SIZE,18), ierr, iv
      integer iter, mdii, mdjj
      integer sub_X,size_X, sub_E,size_E, size_Z

      
      integer,dimension(3) :: oldsize, newsize, starts
      integer :: WS_type, ES_type
      integer :: NS_type, SS_type
      integer :: SWS_type, NES_type
      integer :: SES_type, NWS_type
      integer :: WR_type, ER_type
      integer :: NR_type, SR_type
      integer :: SWR_type, NER_type
      integer :: SER_type, NWR_type
      integer,dimension(maxvar)  :: nbreq
      integer,dimension(16,maxvar) :: treq

      common /nbq_exchange/
     &	nbreq,		! 1->U
     &	treq

      
# if ! defined AGRIF

# else
      real, dimension(:), allocatable ::
     &     buf_snd4,  ibuf_sndN,  buf_snd2,
     &     buf_rev4,  ibuf_revN,  buf_rev2,
     &    jbuf_sndW,             jbuf_sndE,
     &    jbuf_revW,             jbuf_revE,
     &     buf_snd1,  ibuf_sndS,  buf_snd3,
     &     buf_rev1,  ibuf_revS,  buf_rev3
# endif

c**
c      common /buffers_3D/
c     &     buf_snd4,     ibuf_sndN,     buf_snd2,
c     &     buf_rev4,     ibuf_revN,     buf_rev2,
c
c     &    jbuf_sndW,                    jbuf_sndE,
c     &    jbuf_revW,                    jbuf_revE,
c
c     &     buf_snd1,     ibuf_sndS,     buf_snd3,
c     &     buf_rev1,     ibuf_revS,     buf_rev3
c**
!
# include "compute_message_bounds.h"

# ifdef AGRIF
      size_Z=Npts*Npts*(N+1)
      sub_X=(Lm+NSUB_X-1)/NSUB_X
      size_X=Npts*(N+1)*(4+sub_X)
      sub_E=(Mm+NSUB_E-1)/NSUB_E
      size_E=Npts*(N+1)*(4+sub_E)

      Allocate(buf_snd4(size_Z),  ibuf_sndN(size_X),  buf_snd2(size_Z),
     &         buf_rev4(size_Z),  ibuf_revN(size_X),  buf_rev2(size_Z),
     &        jbuf_sndW(size_E),                     jbuf_sndE(size_E),
     &        jbuf_revW(size_E),                     jbuf_revE(size_E),
     &         buf_snd1(size_Z),  ibuf_sndS(size_X),  buf_snd3(size_Z),
     &         buf_rev1(size_Z),  ibuf_revS(size_X),  buf_rev3(size_Z))
# endif
!
      ksize=Npts*Npts*nmax                    ! message sizes for
      isize=Npts*ishft*nmax              ! corner messages and sides
      jsize=Npts*jshft*nmax              ! in XI and ETA directions 
!
! Prepare to receive and send: sides....
!
                            !  Message passing split into two stages
                            !  in order to optimize Send-Recv pairing
                            !  in such a way that if one subdomain
!       do iter=0,1           !  sends message to, say, its WESTERN
!         mdii=mod(ii+iter,2) !  neighbor, that neighbor is preparing
!         mdjj=mod(jj+iter,2) !  to receive this message first (i.e.
                            !  message coming from its EASTERN side),
                            !  rather than send his WEST
                            !  bound message, similarly to the first
                            !  subdomain.
!         write(200+mynode,*) "A=",lbound(A)
!         write(200+mynode,*) "  ",ubound(A)
!         write(200+mynode,*) "shape(A)=",shape(A)

!              iv=1
	     ierr=0
	     nbreq(iv)=0
	     oldsize=shape(A)
!
          if (WEST_INTER) then
	     ! Define size of exchange subarray
	     newsize=(/ Npts, jmax-jmin+1, nmax /)
!	     ! Send bounds Carreful index begin in 0 even in fortran
	     starts = (/ 1, jmin, 1 /) - lbound(A)
	    ! Create and commit mpi type
	     nbreq(iv) = nbreq(iv)+1
	     call mpi_type_create_subarray(3,oldsize,newsize,
     &			starts,MPI_ORDER_FORTRAN, 
     & 			MPI_DOUBLE_PRECISION,
     &			WS_type,ierr)
	     call mpi_type_commit(WS_type,ierr)
	  ! Déclare persistant Comm
	     call mpi_ssend_init(A, 1, WS_type, p_W, 11, 
     &			MPI_COMM_WORLD,treq(nbreq(iv),iv), ierr)

	     newsize=(/ Npts, jmax-jmin+1, nmax /)
	     ! Recv
	     starts = (/ 1-Npts, jmin, 1 /) - lbound(A)
	     nbreq(iv) = nbreq(iv)+1
	     call mpi_type_create_subarray(3,oldsize,newsize,
     &			starts,MPI_ORDER_FORTRAN, 
     & 			MPI_DOUBLE_PRECISION,
     &			WR_type,ierr)
	     call mpi_type_commit(WR_type,ierr)
	     call mpi_recv_init(A, 1, WR_type, p_W, 12, 
     &			MPI_COMM_WORLD,treq(nbreq(iv),iv), ierr)

          endif
          
          if (EAST_INTER) then
	     ! Define size of exchange subarray
	     newsize=(/ Npts, jmax-jmin+1, nmax /)
	     ! Send bounds Carreful index begin in 0 even in fortran
	     starts = (/ Lmmpi-Npts+1, jmin, 1 /) - lbound(A)
	     nbreq(iv) = nbreq(iv)+1
	     call mpi_type_create_subarray(3,oldsize,newsize,
     &			starts,MPI_ORDER_FORTRAN, 
     & 			MPI_DOUBLE_PRECISION,
     &			ES_type,ierr)
	     call mpi_type_commit(ES_type,ierr)
	     call mpi_ssend_init(A, 1, ES_type, p_E, 12, 
     &			MPI_COMM_WORLD,treq(nbreq(iv),iv), ierr)

	     newsize=(/ Npts, jmax-jmin+1, nmax /)
! 	     newsize=(/ 1, jmax-jmin+1, nmax /)
	     ! Recv
	     starts = (/ Lmmpi+1, jmin, 1 /) - lbound(A)
	     nbreq(iv) = nbreq(iv)+1
	     call mpi_type_create_subarray(3,oldsize,newsize,
     &			starts,MPI_ORDER_FORTRAN, 
     & 			MPI_DOUBLE_PRECISION,
     &			ER_type,ierr)
	     call mpi_type_commit(ER_type,ierr)
	     call mpi_recv_init(A, 1, ER_type, p_E, 11, 
     &			MPI_COMM_WORLD,treq(nbreq(iv),iv), ierr)
 
          endif

          if (SOUTH_INTER) then
	    ! Define size of exchange subarray
	    newsize=(/ imax-imin+1, Npts, nmax /)
	    ! Send bounds Carreful index begin in 0 even in fortran
	    starts = (/ imin, 1, 1 /) - lbound(A)
	    ! Create and commit mpi type
	    nbreq(iv) = nbreq(iv)+1 ! increment comm number
	    call mpi_type_create_subarray(3,oldsize,newsize,
     &			starts,MPI_ORDER_FORTRAN, 
     & 			MPI_DOUBLE_PRECISION,
     &			SS_type,ierr)	  
	  call mpi_type_commit(SS_type,ierr)
	  ! Déclare persistant Comm
	  call mpi_ssend_init(A, 1, SS_type, p_S, 13, 
     &			MPI_COMM_WORLD,treq(nbreq(iv),iv), ierr)
	  !! Same couscous with Receive
	  ! Define size of exchange subarray
	  newsize=(/ imax-imin+1, Npts, nmax /)
	  ! Recv bounds Carreful index begin in 0 even in fortran
	  starts = (/ imin, 1-Npts, 1 /) - lbound(A)
	  ! Create and commit mpi type
	  nbreq(iv) = nbreq(iv)+1 ! increment comm number
	  call mpi_type_create_subarray(3,oldsize,newsize,
     &			starts,MPI_ORDER_FORTRAN, 
     & 			MPI_DOUBLE_PRECISION,
     &			SR_type,ierr)	  
	  call mpi_type_commit(SR_type,ierr)
	  ! Déclare persistant Comm
	  call mpi_recv_init(A, 1, SR_type, p_S, 14, 
     &			MPI_COMM_WORLD,treq(nbreq(iv),iv), ierr)
          endif
          
          if (NORTH_INTER) then
	  ! Define size of exchange subarray
	  newsize=(/ imax-imin+1, Npts, nmax /)
	  ! Send bounds Carreful index begin in 0 even in fortran
	  starts = (/ imin, Mmmpi-Npts+1, 1 /) - lbound(A)
	  ! Create and commit mpi type
	  nbreq(iv) = nbreq(iv)+1 ! increment comm number
	  call mpi_type_create_subarray(3,oldsize,newsize,
     &			starts,MPI_ORDER_FORTRAN, 
     & 			MPI_DOUBLE_PRECISION,
     &			NS_type,ierr)	  
	  call mpi_type_commit(NS_type,ierr)
	  ! Déclare persistant Comm
	  call mpi_ssend_init(A, 1, NS_type, p_N, 14, 
     &			MPI_COMM_WORLD,treq(nbreq(iv),iv), ierr)
	  !! Same couscous with Receive
	  ! Define size of exchange subarray
	  newsize=(/ imax-imin+1, Npts, nmax /)
	  ! Recv bounds Carreful index begin in 0 even in fortran
	  starts = (/ imin, Mmmpi+1, 1 /) - lbound(A)
	  ! Create and commit mpi type
	  nbreq(iv) = nbreq(iv)+1 ! increment comm number
	  call mpi_type_create_subarray(3,oldsize,newsize,
     &			starts,MPI_ORDER_FORTRAN, 
     & 			MPI_DOUBLE_PRECISION,
     &			NR_type,ierr)	  
	  call mpi_type_commit(NR_type,ierr)
	  ! Déclare persistant Comm
	  call mpi_recv_init(A, 1, NR_type, p_N, 13, 
     &			MPI_COMM_WORLD,treq(nbreq(iv),iv), ierr)
	  endif	 
	  
	  ! 
! ! ...corners:

          if (SOUTH_INTER .and. WEST_INTER) then
	     ! Define size of exchange subarray
	     newsize=(/ Npts, Npts, nmax /)
!	     ! Send bounds Carreful index begin in 0 even in fortran
	     starts = (/ 1, 1, 1 /) - lbound(A)
	    ! Create and commit mpi type
	     nbreq(iv) = nbreq(iv)+1
	     call mpi_type_create_subarray(3,oldsize,newsize,
     &			starts,MPI_ORDER_FORTRAN, 
     & 			MPI_DOUBLE_PRECISION,
     &			SWS_type,ierr)
	     call mpi_type_commit(SWS_type,ierr)
	     ! Déclare persistant Comm
	     call mpi_ssend_init(A, 1, SWS_type, p_SW, 16, 
     &			MPI_COMM_WORLD,treq(nbreq(iv),iv), ierr)
!	     ! Send bounds Carreful index begin in 0 even in fortran
	     starts = (/ 1-Npts,1-Npts,1 /) - lbound(A)
	    ! Create and commit mpi type
	     nbreq(iv) = nbreq(iv)+1
	     call mpi_type_create_subarray(3,oldsize,newsize,
     &			starts,MPI_ORDER_FORTRAN, 
     & 			MPI_DOUBLE_PRECISION,
     &			SWR_type,ierr)
	     call mpi_type_commit(SWR_type,ierr)
	     ! Déclare persistant Comm
	     call mpi_recv_init(A, 1, SWR_type, p_SW, 15, 
     &			MPI_COMM_WORLD,treq(nbreq(iv),iv), ierr)
          endif
          
          if (NORTH_INTER .and. EAST_INTER) then
	     ! Define size of exchange subarray
	     newsize=(/ Npts, Npts, nmax /)
!	     ! Send bounds Carreful index begin in 0 even in fortran
	     starts = (/ Lmmpi+1-Npts,Mmmpi+1-Npts,1 /) - lbound(A)
	    ! Create and commit mpi type
	     nbreq(iv) = nbreq(iv)+1
	     call mpi_type_create_subarray(3,oldsize,newsize,
     &			starts,MPI_ORDER_FORTRAN, 
     & 			MPI_DOUBLE_PRECISION,
     &			NES_type,ierr)
	     call mpi_type_commit(NES_type,ierr)
	     ! Déclare persistant Comm
	     call mpi_ssend_init(A, 1, NES_type, p_NE, 15, 
     &			MPI_COMM_WORLD,treq(nbreq(iv),iv), ierr)
     
!	     ! Send bounds Carreful index begin in 0 even in fortran
	     starts = (/ Lmmpi+1,Mmmpi+1,1 /) - lbound(A)
	    ! Create and commit mpi type
	     nbreq(iv) = nbreq(iv)+1
	     call mpi_type_create_subarray(3,oldsize,newsize,
     &			starts,MPI_ORDER_FORTRAN, 
     & 			MPI_DOUBLE_PRECISION,
     &			NER_type,ierr)
	     call mpi_type_commit(NER_type,ierr)
	     ! Déclare persistant Comm
	     call mpi_recv_init(A, 1, NER_type, p_NE, 16, 
     &			MPI_COMM_WORLD,treq(nbreq(iv),iv), ierr)
          endif

          if (SOUTH_INTER .and. EAST_INTER) then
	     ! Define size of exchange subarray
	     newsize=(/ Npts, Npts, nmax /)
!	     ! Send bounds Carreful index begin in 0 even in fortran
	     starts = (/ Lmmpi+1-Npts,1,1 /) - lbound(A)
	    ! Create and commit mpi type
	     nbreq(iv) = nbreq(iv)+1
	     call mpi_type_create_subarray(3,oldsize,newsize,
     &			starts,MPI_ORDER_FORTRAN, 
     & 			MPI_DOUBLE_PRECISION,
     &			SES_type,ierr)
	     call mpi_type_commit(SES_type,ierr)
	     ! Déclare persistant Comm
	     call mpi_ssend_init(A, 1, SES_type, p_SE, 18, 
     &			MPI_COMM_WORLD,treq(nbreq(iv),iv), ierr)
!	     ! Send bounds Carreful index begin in 0 even in fortran
	     starts = (/ Lmmpi+1,1-Npts,1 /) - lbound(A)
	    ! Create and commit mpi type
	     nbreq(iv) = nbreq(iv)+1
	     call mpi_type_create_subarray(3,oldsize,newsize,
     &			starts,MPI_ORDER_FORTRAN, 
     & 			MPI_DOUBLE_PRECISION,
     &			SER_type,ierr)
	     call mpi_type_commit(SER_type,ierr)
	     ! Déclare persistant Comm
	     call mpi_recv_init(A, 1, SER_type, p_SE, 17, 
     &			MPI_COMM_WORLD,treq(nbreq(iv),iv), ierr)
          endif
          
          if (NORTH_INTER .and. WEST_INTER) then
	     ! Define size of exchange subarray
	     newsize=(/ Npts, Npts, nmax /)
!	     ! Send bounds Carreful index begin in 0 even in fortran
	     starts = (/ 1,Mmmpi+1-Npts,1 /) - lbound(A)
	    ! Create and commit mpi type
	     nbreq(iv) = nbreq(iv)+1
	     call mpi_type_create_subarray(3,oldsize,newsize,
     &			starts,MPI_ORDER_FORTRAN, 
     & 			MPI_DOUBLE_PRECISION,
     &			NWS_type,ierr)
	     call mpi_type_commit(NWS_type,ierr)
	     ! Déclare persistant Comm
	     call mpi_ssend_init(A, 1, NWS_type, p_NW, 17, 
     &			MPI_COMM_WORLD,treq(nbreq(iv),iv), ierr)
!	     ! Send bounds Carreful index begin in 0 even in fortran
	     starts = (/ 1-Npts,Mmmpi+1,1 /) - lbound(A)
	    ! Create and commit mpi type
	     nbreq(iv) = nbreq(iv)+1
	     call mpi_type_create_subarray(3,oldsize,newsize,
     &			starts,MPI_ORDER_FORTRAN, 
     & 			MPI_DOUBLE_PRECISION,
     &			NWR_type,ierr)
	     call mpi_type_commit(NWR_type,ierr)
	     ! Déclare persistant Comm
	     call mpi_recv_init(A, 1, NWR_type, p_NW, 18, 
     &			MPI_COMM_WORLD,treq(nbreq(iv),iv), ierr)
          endif

!       print *,mynode, "(1) treq  =","nbreq(",iv,")=",
!      & nbreq(iv),treq(1:nbreq(iv),iv),ierr
!
!


# if defined AGRIF
        DeAllocate(buf_snd4,  ibuf_sndN,  buf_snd2,
     &     buf_rev4,  ibuf_revN,  buf_rev2,
     &    jbuf_sndW,                      jbuf_sndE,
     &    jbuf_revW,                      jbuf_revE,
     &     buf_snd1,  ibuf_sndS,  buf_snd3,
     &     buf_rev1,  ibuf_revS,  buf_rev3)
# endif

      return
      end

      
      
      !--------------------------------------------------
# ifndef MP_3PTS
      subroutine MessPass3D_tile_nbq (iv,Istr,Iend,Jstr,Jend, A, nmax)
# else
      subroutine MessPass3D_3pts_tile_nbq (iv,Istr,Iend,Jstr,Jend,A,nmax)
# endif
!
!======================================================================
!
! This subroutine is designed for ROMS-MPI code. It exchanges domain
! boundary information, including 2 (or 3) ghost-cells in each
! direction.
!
! Ping Wang 9/15/99.
! Patrick Marchesiello 2012: generic number of ghost-cells Npts
!
!======================================================================
!
!     implicit none
# include "param.h"
# include "scalars.h"
# include "mpi_cpl.h"
      include 'mpif.h'
!
! Nb of boundary points involved in communication
!
      integer Npts,ipts,jpts
# ifndef MP_3PTS
      parameter (Npts=2)
# else
      parameter (Npts=3)
# endif
      integer :: maxvar
      parameter(maxvar=10)
      integer nmax
      real A(GLOBAL_2D_ARRAY,nmax)
CSDISTRIBUTE_RESHAPE A(BLOCK_PATTERN) BLOCK_CLAUSE
      integer Istr,Iend,Jstr,Jend, i,j,k, isize,jsize,ksize,
     &        req(18), status(MPI_STATUS_SIZE,18), ierr, iv
      integer iter, mdii, mdjj
      integer sub_X,size_X, sub_E,size_E, size_Z

      
      integer,dimension(3) :: oldsize, newsize, starts
      integer :: WS_type, ES_type
      integer :: NS_type, SS_type
      integer :: SWS_type, NES_type
      integer :: SES_type, NWS_type
      integer :: WR_type, ER_type
      integer :: NR_type, SR_type
      integer :: SWR_type, NER_type
      integer :: SER_type, NWR_type
      integer,dimension(maxvar)  :: nbreq
      integer,dimension(16,maxvar) :: treq

      common /nbq_exchange/
     &	nbreq,		! 1->U
     &	treq

      
# if ! defined AGRIF

# else
      real, dimension(:), allocatable ::
     &     buf_snd4,  ibuf_sndN,  buf_snd2,
     &     buf_rev4,  ibuf_revN,  buf_rev2,
     &    jbuf_sndW,             jbuf_sndE,
     &    jbuf_revW,             jbuf_revE,
     &     buf_snd1,  ibuf_sndS,  buf_snd3,
     &     buf_rev1,  ibuf_revS,  buf_rev3
# endif

c**
c      common /buffers_3D/
c     &     buf_snd4,     ibuf_sndN,     buf_snd2,
c     &     buf_rev4,     ibuf_revN,     buf_rev2,
c
c     &    jbuf_sndW,                    jbuf_sndE,
c     &    jbuf_revW,                    jbuf_revE,
c
c     &     buf_snd1,     ibuf_sndS,     buf_snd3,
c     &     buf_rev1,     ibuf_revS,     buf_rev3
c**
!"
# include "compute_message_bounds.h"

# ifdef AGRIF
      size_Z=Npts*Npts*(N+1)
      sub_X=(Lm+NSUB_X-1)/NSUB_X
      size_X=Npts*(N+1)*(4+sub_X)
      sub_E=(Mm+NSUB_E-1)/NSUB_E
      size_E=Npts*(N+1)*(4+sub_E)

      Allocate(buf_snd4(size_Z),  ibuf_sndN(size_X),  buf_snd2(size_Z),
     &         buf_rev4(size_Z),  ibuf_revN(size_X),  buf_rev2(size_Z),
     &        jbuf_sndW(size_E),                     jbuf_sndE(size_E),
     &        jbuf_revW(size_E),                     jbuf_revE(size_E),
     &         buf_snd1(size_Z),  ibuf_sndS(size_X),  buf_snd3(size_Z),
     &         buf_rev1(size_Z),  ibuf_revS(size_X),  buf_rev3(size_Z))
# endif
!
      ksize=Npts*Npts*nmax                    ! message sizes for
      isize=Npts*ishft*nmax              ! corner messages and sides
      jsize=Npts*jshft*nmax              ! in XI and ETA directions 
      
            
      ierr = 0
!       write(200+mynode,*) "nbreq(",iv,")=",nbreq(iv)
!       write(200+mynode,*) treq(1:nbreq(iv),iv)
!       print *,mynode, "     treq =","nbreq(",iv,")=",nbreq(iv),
!      &	treq(1:nbreq(iv),iv),ierr
!       cff=1.0
!       print *,mynode,WEST_INTER,EAST_INTER,SOUTH_INTER,NORTH_INTER
!       print *,mynode," nbreq(",iv,")=",nbreq(iv),"---"

      ! Send-Recv
      call MPI_STARTALL(nbreq(iv),treq(1:nbreq(iv),iv), ierr)
!       print *,mynode, " call MPI_STARTAL OK",ierr
       ! Waitall
      call MPI_WAITALL(nbreq(iv), 
     &	treq(1:nbreq(iv),iv),
     &	status(:,1:nbreq(iv)), ierr)

     
!       print *,mynode," : ",status(:,1)
!       print *,mynode," : ",status(:,2)
!       print *,mynode," : ",status(:,3)
!       print *,mynode," : ",status(:,4)
!       print *,mynode,"nbreq(",iv,")=",nbreq(iv)
!    
      return
      end
      
# ifndef MP_3PTS
#  define MP_3PTS
#  include "MessPass3D_nbq.F"
#  undef MP_3PTS
# endif

#else
      subroutine MessPass3D_init_nbq_empty
      return
      end
#endif

