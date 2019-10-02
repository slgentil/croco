!
! $Id: MessPassFloats.F 1143 2013-05-17 08:17:23Z serena $
!
#include "cppdefs.h"

! aponte flt (
#if defined FLOATS && defined MPI

      subroutine exchange_floats ()
!
!======================================================================
!
! exchange float data: track, trackaux, trackinfo
! Could be made much faster if buffer where exchanged only if necessary
! also won't be performant with a larger number of floats
!
!======================================================================
!     implicit none
# include "param.h"
# include "scalars.h"
# include "mpi_roms.h"
# include "floats.h"
      include 'mpif.h'
!
! Nb of boundary points involved in communication
!
      integer Nptsmax,ipts,jpts
      integer string_len
      character*60  string_error
      parameter (Nptsmax=1000)
CSDISTRIBUTE_RESHAPE A(BLOCK_PATTERN) BLOCK_CLAUSE
!      integer Istr,Iend,Jstr,Jend, i,j, isize,jsize,ksize,
      integer iter,
     &        req(8), status(MPI_STATUS_SIZE,8), ierr, mdii,mdjj
      integer iflt, ifltmax, Npts, i0, i, j, k
      integer node_float
! buffers for North,South,East,West     
      real ibuf_sndN(Nptsmax), ibuf_revN(Nptsmax), 
     &     ibuf_sndS(Nptsmax), ibuf_revS(Nptsmax),
     &     ibuf_sndW(Nptsmax), ibuf_revW(Nptsmax),
     &     ibuf_sndE(Nptsmax), ibuf_revE(Nptsmax)
! buffers for corner SW, NE, SE, NW
      real ibuf_sndSW(Nptsmax), ibuf_sndNE(Nptsmax), 
     &     ibuf_revSW(Nptsmax), ibuf_revNE(Nptsmax),
     &     ibuf_sndSE(Nptsmax), ibuf_sndNW(Nptsmax), 
     &     ibuf_revSE(Nptsmax), ibuf_revNW(Nptsmax)
     
! init buffer arrays to 0
      ibuf_sndN(:)=0.0D0
      ibuf_revN(:)=0.0D0
      ibuf_sndS(:)=0.0D0
      ibuf_revS(:)=0.0D0
      ibuf_sndW(:)=0.0D0
      ibuf_revW(:)=0.0D0
      ibuf_sndE(:)=0.0D0
      ibuf_revE(:)=0.0D0  
      ibuf_sndSW(:)=0.0D0
      ibuf_revSW(:)=0.0D0  
      ibuf_sndNE(:)=0.0D0
      ibuf_revNE(:)=0.0D0
      ibuf_sndSE(:)=0.0D0
      ibuf_revSE(:)=0.0D0
      ibuf_sndNW(:)=0.0D0
      ibuf_revNW(:)=0.0D0

                            !  Message passing split into two stages
                            !  in order to optimize Send-Recv pairing
                            !  in such a way that if one subdomain
      do iter=0,1           !  sends message to, say, its WESTERN
        mdii=mod(ii+iter,2) !  neighbor, that neighbor is preparing
        mdjj=mod(jj+iter,2) !  to exchange_floatsreceive this message first (i.e.
                            !  message coming from its EASTERN side),
                            !  rather than send his WEST
                            !  bound message, similarly to the first
                            !  subdomain.

! Number of data transferred by floats (iflt+track+trackaux+trackinfo)      
        Npts = 1 + NFV*(NFT+1) + NDIAGS + 3
!
! Prepare to receive and send: sides....
        if (mdii.eq.0) then
          if (WEST_INTER) then
            ifltmax=0
            i0=1

            do iflt=1,nfloats   
              if ( trackinfo(1,iflt) .eq. 1 .and.
     &             trackinfo(2,iflt) .eq. mynode ) then 
                node_float= modulo(floor((track(ixgrd,nf,iflt)-0.5) 
     &                          /dble(Lm)),NP_XI)
                node_float= node_float
     &                 + floor((track(iygrd,nf,iflt)-0.5)
     &                       /dble(Mm ))*NP_XI
                if ( node_float .eq. p_W ) then
C                   write(*,'(A,i4,A,i4,A,i4)') 'exchange float number ',
C      &                 iflt,' from proc ',mynode,' to proc ',p_W            
                  if (i0+Npts .ge. Nptsmax) then
                   write(*,*)'exchange_floats:increase Nptsmax>',i0+Npts
                    call MPI_Finalize (ierr)
                    stop
                  endif
                  ! iflt
                  ibuf_sndW(i0+1) = dble(iflt)
                  i0=i0+1
                  ! track EW periodic
                  if (ii == 0 .and. track(ixgrd,nf,iflt) < iminmpi) then
                    do i=0,NFT 
                     track(ixgrd,i,iflt)=track(ixgrd,i,iflt)+dble(LLm-1)
                    enddo 
                  endif         
                  do j=1,NFV
                    do i=0,NFT
                      ibuf_sndW(i0+j+i*NFV) = track(j,i,iflt)
                    enddo
                  enddo
                  i0=i0+NFV*(NFT+1)
                  ! trackaux
                  do j=1,NDIAGS
                    ibuf_sndW(i0+j) = trackaux(j,iflt)
                  enddo
                  i0=i0+NDIAGS
                  ! update node number
                  trackinfo(2,iflt)=node_float
                  ! trackinfo
                  do j=1,3
                    ibuf_sndW(i0+j) = dble(trackinfo(j,iflt))
                  enddo
                  i0=i0+3
                  ! number of floats exchanged
                  ifltmax=ifltmax+1

                endif
              endif
            enddo                 
          
            ibuf_sndW(1)=dble(ifltmax)
            call MPI_Irecv (ibuf_revW, Nptsmax, MPI_DOUBLE_PRECISION,
     &                        p_W, 2, MPI_COMM_WORLD, req(1), ierr)  
            if (ierr .ne. 0) then
               call MPI_Error_string(ierr,string_error,string_len)
               print*,'---> MessPassFloats: 1.1 (recv):',string_error
               call MPI_Finalize (ierr)
               stop
            endif      
            call MPI_Send  (ibuf_sndW, Nptsmax, MPI_DOUBLE_PRECISION,
     &                        p_W, 1, MPI_COMM_WORLD,         ierr)  
            if (ierr .ne. 0) then
               call MPI_Error_string(ierr,string_error,string_len)
               print*,'---> MessPassFloats: 1.1 (send):',string_error
               call MPI_Finalize (ierr)
               stop
            endif       
          endif
        else
          if (EAST_INTER) then
            ifltmax=0
            i0=1
            do iflt=1,nfloats          
              if ( trackinfo(1,iflt) .eq. 1 .and.
     &             trackinfo(2,iflt) .eq. mynode ) then
                node_float= mod(int((track(ixgrd,nf,iflt)-0.5)
     &                          /dble(Lm)),NP_XI)
                node_float= node_float
     &                 + int((track(iygrd,nf,iflt)-0.5)
     &                       /dble(Mm))*NP_XI
                if ( node_float .eq. p_E ) then
C                   write(*,'(A,i4,A,i4,A,i4)') 'exchange float number ',
C      &                 iflt,' from proc ',mynode,' to proc ',p_E    
                  if (i0+Npts .ge. Nptsmax) then
                   write(*,*)'exchange_floats:increase Nptsmax>',i0+Npts
                    call MPI_Finalize (ierr)
                    stop
                  endif
                  ! iflt
                  ibuf_sndE(i0+1) = dble(iflt)
                  i0=i0+1
                  ! track EW periodic
                  if (ii == NP_XI-1 .and. track(ixgrd,nf,iflt)>imaxmpi) then 
                    do i=0,NFT 
                     track(ixgrd,i,iflt)=track(ixgrd,i,iflt)-dble(LLm-1)
                    enddo    
                  endif         
                  do j=1,NFV
                    do i=0,NFT
                      ibuf_sndE(i0+j+i*NFV) = track(j,i,iflt)
                    enddo
                  enddo
                  i0=i0+NFV*(NFT+1)
                  ! trackaux
                  do j=1,NDIAGS
                    ibuf_sndE(i0+j) = trackaux(j,iflt)
                  enddo
                  i0=i0+NDIAGS
                  ! update node number
                  trackinfo(2,iflt)=node_float
                  ! trackinfo
                  do j=1,3
                    ibuf_sndE(i0+j) = dble(trackinfo(j,iflt))
                  enddo
                  i0=i0+3
                  ! number of floats exchanged
                  ifltmax=ifltmax+1

                endif
              endif
            enddo
            ibuf_sndE(1)=dble(ifltmax)
            call MPI_Irecv (ibuf_revE, Nptsmax, MPI_DOUBLE_PRECISION,
     &                        p_E, 1, MPI_COMM_WORLD, req(2), ierr) 
            if (ierr .ne. 0) then
              call MPI_Error_string(ierr,string_error,string_len)
              print*,'---> MessPassFloats: 1.2 (send):',string_error
              call MPI_Finalize (ierr)
              stop
            endif             
            call MPI_Send  (ibuf_sndE, Nptsmax, MPI_DOUBLE_PRECISION,
     &                         p_E, 2, MPI_COMM_WORLD,         ierr)  
            if (ierr .ne. 0) then
              call MPI_Error_string(ierr,string_error,string_len)
              print*,'---> MessPassFloats: 1.2 (send):',string_error
              call MPI_Finalize (ierr)
              stop
            endif   
          endif
        endif

        if (mdjj.eq.0) then
          if (SOUTH_INTER) then
            ifltmax=0
            i0=1
            do iflt=1,nfloats
              if ( trackinfo(1,iflt) .eq. 1 .and.
     &             trackinfo(2,iflt) .eq. mynode ) then
                node_float= mod(int((track(ixgrd,nf,iflt)-0.5)
     &                          /dble(Lm)),NP_XI)
                node_float= node_float
     &                 + int((track(iygrd,nf,iflt)-0.5)
     &                       /dble(Mm))*NP_XI
                if ( node_float .eq. p_S ) then
C                   write(*,'(A,i4,A,i4,A,i4)') 'exchange float number ',
C      &                 iflt,' from proc ',mynode,' to proc ',p_S    
                  if (i0+Npts .ge. Nptsmax) then
                   write(*,*)'exchange_floats:increase Nptsmax>',i0+Npts
                    call MPI_Finalize (ierr)
                    stop
                  endif             
                  ! iflt
                  ibuf_sndS(i0+1) = dble(iflt)
                  i0=i0+1
                  ! track NS periodic
                  if (jj == 0 .and. track(iygrd,nf,iflt) < jminmpi) then
                    do i=0,NFT 
                     track(iygrd,i,iflt)=track(iygrd,i,iflt)+dble(MMm-1)
                    enddo 
                  endif         
            
                  do j=1,NFV
                    do i=0,NFT
                      ibuf_sndS(i0+j+i*NFV) = track(j,i,iflt)
                    enddo
                  enddo
                  i0=i0+NFV*(NFT+1)
                  ! trackaux
                  do j=1,NDIAGS
                    ibuf_sndS(i0+j) = trackaux(j,iflt)
                  enddo
                  i0=i0+NDIAGS
                  ! update node number
                  trackinfo(2,iflt)=node_float
                  ! trackinfo
                  do j=1,3
                    ibuf_sndS(i0+j) = dble(trackinfo(j,iflt))
                  enddo
                  i0=i0+3
                  ! number of floats exchanged
                  ifltmax=ifltmax+1

                endif
              endif
            enddo
            ibuf_sndS(1)=dble(ifltmax)
            call MPI_Irecv (ibuf_revS, Nptsmax, MPI_DOUBLE_PRECISION,
     &                        p_S, 4, MPI_COMM_WORLD, req(3), ierr) 
            if (ierr .ne. 0) then
              call MPI_Error_string(ierr,string_error,string_len)
              print*,'---> MessPassFloats: 1.3 (send):',string_error
              call MPI_Finalize (ierr)
              stop
            endif             
            call MPI_Send  (ibuf_sndS, Nptsmax, MPI_DOUBLE_PRECISION,
     &                        p_S, 3, MPI_COMM_WORLD,         ierr)
            if (ierr .ne. 0) then
              call MPI_Error_string(ierr,string_error,string_len)
              print*,'---> MessPassFloats: 1.3 (send):',string_error
              call MPI_Finalize (ierr)
              stop
            endif  
          endif
        else
          if (NORTH_INTER) then
            ifltmax=0
            i0=1
            do iflt=1,nfloats
              if ( trackinfo(1,iflt) .eq. 1 .and.
     &             trackinfo(2,iflt) .eq. mynode ) then
                node_float= mod(int((track(ixgrd,nf,iflt)-0.5)
     &                          /dble(Lm)),NP_XI)
                node_float= node_float
     &                 + int((track(iygrd,nf,iflt)-0.5)
     &                       /dble(Mm))*NP_XI
                if ( node_float .eq. p_N ) then
C                   write(*,'(A,i4,A,i4,A,i4)') 'exchange float number ',
C      &                 iflt,' from proc ',mynode,' to proc ',p_N    
                  if (i0+Npts .ge. Nptsmax) then
                   write(*,*)'exchange_floats:increase Nptsmax>',i0+Npts
                    call MPI_Finalize (ierr)
                    stop
                  endif             
                  ! iflt
                  ibuf_sndN(i0+1) = dble(iflt)
                  i0=i0+1
                  ! track NS periodic
                  if (jj == NP_ETA-1 .and. track(iygrd,nf,iflt)>jmaxmpi) then 
                    do i=0,NFT 
                     track(iygrd,i,iflt)=track(iygrd,i,iflt)-dble(MMm-1)
                    enddo    
                  endif         
                  do j=1,NFV
                    do i=0,NFT
                      ibuf_sndN(i0+j+i*NFV) = track(j,i,iflt)
                    enddo
                  enddo
                  i0=i0+NFV*(NFT+1)
                  ! trackaux
                  do j=1,NDIAGS
                    ibuf_sndN(i0+j) = trackaux(j,iflt)
                  enddo
                  i0=i0+NDIAGS
                  ! update node number
                  trackinfo(2,iflt)=node_float
                  ! trackinfo
                  do j=1,3
                    ibuf_sndN(i0+j) = dble(trackinfo(j,iflt))
                  enddo
                  i0=i0+3
                  ! number of floats exchanged
                  ifltmax=ifltmax+1

                endif
              endif
            enddo
            ibuf_sndN(1)=dble(ifltmax)
            call MPI_Irecv (ibuf_revN, Nptsmax, MPI_DOUBLE_PRECISION,
     &                         p_N, 3, MPI_COMM_WORLD, req(4), ierr) 
            if (ierr .ne. 0) then
              call MPI_Error_string(ierr,string_error,string_len)
              print*,'---> MessPassFloats: 1.4 (send):',string_error
              call MPI_Finalize (ierr)
              stop
            endif          
            call MPI_Send  (ibuf_sndN, Nptsmax, MPI_DOUBLE_PRECISION,
     &                         p_N, 4, MPI_COMM_WORLD,         ierr)
            if (ierr .ne. 0) then
              call MPI_Error_string(ierr,string_error,string_len)
              print*,'---> MessPassFloats: 1.4 (send):',string_error
              call MPI_Finalize (ierr)
              stop
            endif   
          endif
        endif

!
! ...corners:
!
        if (mdii.eq.0) then
          if (SOUTH_INTER .and. WEST_INTER) then

            ifltmax=0
            i0=1
            do iflt=1,nfloats
              if ( trackinfo(1,iflt) .eq. 1 .and.
     &             trackinfo(2,iflt) .eq. mynode ) then
                node_float= mod(int((track(ixgrd,nf,iflt)-0.5)
     &                          /dble(Lm)),NP_XI)
                node_float= node_float
     &                 + int((track(iygrd,nf,iflt)-0.5)
     &                       /dble(Mm))*NP_XI
                if ( node_float .eq. p_SW ) then
C                   write(*,'(A,i4,A,i4,A,i4)') 'exchange float number ',
C      &                 iflt,' from proc ',mynode,' to proc ',p_SW    
                  if (i0+Npts .ge. Nptsmax) then
                   write(*,*)'exchange_floats:increase Nptsmax>',i0+Npts
                     call MPI_Finalize (ierr)
                     stop
                  endif             
                  ! iflt
                  ibuf_sndSW(i0+1) = dble(iflt)
                  i0=i0+1

                  ! track EW periodic
                  if (ii == 0 .and. track(ixgrd,nf,iflt) < iminmpi) then
                    do i=0,NFT 
                     track(ixgrd,i,iflt)=track(ixgrd,i,iflt)+dble(LLm-1)
                    enddo 
                  endif  
                  ! track NS periodic
                  if (jj == 0 .and. track(iygrd,nf,iflt) < jminmpi) then
                    do i=0,NFT 
                     track(iygrd,i,iflt)=track(iygrd,i,iflt)+dble(MMm-1)
                    enddo 
                  endif          
   
                  do j=1,NFV
                     do i=0,NFT
                       ibuf_sndSW(i0+j+i*NFV) = track(j,i,iflt)
                     enddo
                  enddo
                  i0=i0+NFV*(NFT+1)
                  ! trackaux
                  do j=1,NDIAGS
                    ibuf_sndSW(i0+j) = trackaux(j,iflt)
                  enddo
                  i0=i0+NDIAGS
                  ! update node number
                  trackinfo(2,iflt)=node_float
                  ! trackinfo
                  do j=1,3
                    ibuf_sndSW(i0+j) = dble(trackinfo(j,iflt))
                  enddo
                  i0=i0+3
                  ! number of floats exchanged
                  ifltmax=ifltmax+1

                endif
              endif
            enddo
            ibuf_sndSW(1)=dble(ifltmax)

            call MPI_Irecv (ibuf_revSW, Nptsmax, MPI_DOUBLE_PRECISION,
     &                         p_SW, 6, MPI_COMM_WORLD, req(5), ierr) 
            if (ierr .ne. 0) then
              call MPI_Error_string(ierr,string_error,string_len)
              print*,'---> MessPassFloats: 1.5 (send):',string_error
              call MPI_Finalize (ierr)
              stop
            endif          
            call MPI_Send  (ibuf_sndSW, Nptsmax, MPI_DOUBLE_PRECISION,
     &                         p_SW, 5, MPI_COMM_WORLD,         ierr)
            if (ierr .ne. 0) then
              call MPI_Error_string(ierr,string_error,string_len)
              print*,'---> MessPassFloats: 1.5 (send):',string_error
              call MPI_Finalize (ierr)
              stop
            endif 
          endif
        else
          if (NORTH_INTER .and. EAST_INTER) then
            ifltmax=0
            i0=1
            do iflt=1,nfloats
              if ( trackinfo(1,iflt) .eq. 1 .and.
     &             trackinfo(2,iflt) .eq. mynode ) then
                node_float= mod(int((track(ixgrd,nf,iflt)-0.5)
     &                          /dble(Lm)),NP_XI)
                node_float= node_float
     &                 + int((track(iygrd,nf,iflt)-0.5)
     &                       /dble(Mm))*NP_XI
                if ( node_float .eq. p_NE ) then
C                   write(*,'(A,i4,A,i4,A,i4)') 'exchange float number ',
C      &                 iflt,' from proc ',mynode,' to proc ',p_NE    
                  if (i0+Npts .ge. Nptsmax) then
                   write(*,*)'exchange_floats:increase Nptsmax>',i0+Npts
                     call MPI_Finalize (ierr)
                     stop
                  endif             
                  ! iflt
                  ibuf_sndNE(i0+1) = dble(iflt)
                  i0=i0+1
                  ! track EW periodic
                  if (ii == NP_XI-1 .and. track(ixgrd,nf,iflt)>imaxmpi) then 
                    do i=0,NFT 
                     track(ixgrd,i,iflt)=track(ixgrd,i,iflt)-dble(LLm-1)
                    enddo    
                  endif         
                  ! track NS periodic
                  if (jj == NP_ETA-1 .and. track(iygrd,nf,iflt)>jmaxmpi) then 
                    do i=0,NFT 
                     track(iygrd,i,iflt)=track(iygrd,i,iflt)-dble(MMm-1)
                    enddo    
                  endif         
    
                  do j=1,NFV
                     do i=0,NFT
                       ibuf_sndNE(i0+j+i*NFV) = track(j,i,iflt)
                     enddo
                  enddo
                  i0=i0+NFV*(NFT+1)
                  ! trackaux
                  do j=1,NDIAGS
                    ibuf_sndNE(i0+j) = trackaux(j,iflt)
                  enddo
                  i0=i0+NDIAGS
                  ! update node number
                  trackinfo(2,iflt)=node_float
                  ! trackinfo
                  do j=1,3
                    ibuf_sndNE(i0+j) = dble(trackinfo(j,iflt))
                  enddo
                  i0=i0+3
                  ! number of floats exchanged
                  ifltmax=ifltmax+1

                endif
              endif
            enddo
            ibuf_sndNE(1)=dble(ifltmax)

            call MPI_Irecv (ibuf_revNE, Nptsmax, MPI_DOUBLE_PRECISION,
     &                         p_NE, 5, MPI_COMM_WORLD, req(6), ierr) 
            if (ierr .ne. 0) then
              call MPI_Error_string(ierr,string_error,string_len)
              print*,'---> MessPassFloats: 1.6 (send):',string_error
              call MPI_Finalize (ierr)
              stop
            endif          
            call MPI_Send  (ibuf_sndNE, Nptsmax, MPI_DOUBLE_PRECISION,
     &                         p_NE, 6, MPI_COMM_WORLD,         ierr)
            if (ierr .ne. 0) then
              call MPI_Error_string(ierr,string_error,string_len)
              print*,'---> MessPassFloats: 1.6 (send):',string_error
              call MPI_Finalize (ierr)
              stop
            endif  

          endif
        endif

        if (mdii.eq.1) then
          if (SOUTH_INTER .and. EAST_INTER) then

            ifltmax=0
            i0=1
            do iflt=1,nfloats
              if ( trackinfo(1,iflt) .eq. 1 .and.
     &             trackinfo(2,iflt) .eq. mynode ) then
                node_float= mod(int((track(ixgrd,nf,iflt)-0.5)
     &                          /dble(Lm)),NP_XI)
                node_float= node_float
     &                 + int((track(iygrd,nf,iflt)-0.5)
     &                       /dble(Mm))*NP_XI
                if ( node_float .eq. p_SE ) then
C                   write(*,'(A,i4,A,i4,A,i4)') 'exchange float number ',
C      &                 iflt,' from proc ',mynode,' to proc ',p_SE    
                  if (i0+Npts .ge. Nptsmax) then
                   write(*,*)'exchange_floats:increase Nptsmax>',i0+Npts
                     call MPI_Finalize (ierr)
                     stop
                  endif             
                  ! iflt
                  ibuf_sndSE(i0+1) = dble(iflt)
                  i0=i0+1
                  ! track NS periodic
                  if (jj == 0 .and. track(iygrd,nf,iflt) < jminmpi) then
                    do i=0,NFT 
                     track(iygrd,i,iflt)=track(iygrd,i,iflt)+dble(MMm-1)
                    enddo 
                  endif      
                  ! track EW periodic
                  if (ii == NP_XI-1 .and. track(ixgrd,nf,iflt)>imaxmpi) then 
                    do i=0,NFT 
                     track(ixgrd,i,iflt)=track(ixgrd,i,iflt)-dble(LLm-1)
                    enddo    
                  endif            
  
                  do j=1,NFV
                     do i=0,NFT
                       ibuf_sndSE(i0+j+i*NFV) = track(j,i,iflt)
                     enddo
                  enddo
                  i0=i0+NFV*(NFT+1)
                  ! trackaux
                  do j=1,NDIAGS
                    ibuf_sndSE(i0+j) = trackaux(j,iflt)
                  enddo
                  i0=i0+NDIAGS
                  ! update node number
                  trackinfo(2,iflt)=node_float
                  ! trackinfo
                  do j=1,3
                    ibuf_sndSE(i0+j) = dble(trackinfo(j,iflt))
                  enddo
                  i0=i0+3
                  ! number of floats exchanged
                  ifltmax=ifltmax+1

                endif
              endif
            enddo
            ibuf_sndSE(1)=dble(ifltmax)

            call MPI_Irecv (ibuf_revSE, Nptsmax, MPI_DOUBLE_PRECISION,
     &                         p_SE, 8, MPI_COMM_WORLD, req(7), ierr) 
            if (ierr .ne. 0) then
              call MPI_Error_string(ierr,string_error,string_len)
              print*,'---> MessPassFloats: 1.7 (send):',string_error
              call MPI_Finalize (ierr)
              stop
            endif          
            call MPI_Send  (ibuf_sndSE, Nptsmax, MPI_DOUBLE_PRECISION,
     &                         p_SE, 7, MPI_COMM_WORLD,         ierr)
            if (ierr .ne. 0) then
              call MPI_Error_string(ierr,string_error,string_len)
              print*,'---> MessPassFloats: 1.7 (send):',string_error
              call MPI_Finalize (ierr)
              stop
            endif 
          endif
        else
          if (NORTH_INTER .and. WEST_INTER) then

            ifltmax=0
            i0=1
            do iflt=1,nfloats
              if ( trackinfo(1,iflt) .eq. 1 .and.
     &             trackinfo(2,iflt) .eq. mynode ) then
                node_float= mod(int((track(ixgrd,nf,iflt)-0.5)
     &                          /dble(Lm)),NP_XI)
                node_float= node_float
     &                 + int((track(iygrd,nf,iflt)-0.5)
     &                       /dble(Mm))*NP_XI
                if ( node_float .eq. p_NW ) then
C                   write(*,'(A,i4,A,i4,A,i4)') 'exchange float number ',
C      &                 iflt,' from proc ',mynode,' to proc ',p_NW    
                  if (i0+Npts .ge. Nptsmax) then
                   write(*,*)'exchange_floats:increase Nptsmax>',i0+Npts
                     call MPI_Finalize (ierr)
                     stop
                  endif             
                  ! iflt
                  ibuf_sndNW(i0+1) = dble(iflt)
                  i0=i0+1
                  ! track NS periodic
                  if (jj == NP_ETA-1 .and. track(iygrd,nf,iflt)>jmaxmpi) then 
                    do i=0,NFT 
                     track(iygrd,i,iflt)=track(iygrd,i,iflt)-dble(MMm-1)
                    enddo    
                  endif         
                  ! track EW periodic
                  if (ii == 0 .and. track(ixgrd,nf,iflt) < iminmpi) then
                    do i=0,NFT 
                     track(ixgrd,i,iflt)=track(ixgrd,i,iflt)+dble(LLm-1)
                    enddo 
                  endif         
  
                  do j=1,NFV
                     do i=0,NFT
                       ibuf_sndNW(i0+j+i*NFV) = track(j,i,iflt)
                     enddo
                  enddo
                  i0=i0+NFV*(NFT+1)
                  ! trackaux
                  do j=1,NDIAGS
                    ibuf_sndNW(i0+j) = trackaux(j,iflt)
                  enddo
                  i0=i0+NDIAGS
                  ! update node number
                  trackinfo(2,iflt)=node_float
                  ! trackinfo
                  do j=1,3
                    ibuf_sndNW(i0+j) = dble(trackinfo(j,iflt))
                  enddo
                  i0=i0+3
                  ! number of floats exchanged
                  ifltmax=ifltmax+1

                endif
              endif
            enddo
            ibuf_sndNW(1)=dble(ifltmax)

            call MPI_Irecv (ibuf_revNW, Nptsmax, MPI_DOUBLE_PRECISION,
     &                         p_NW, 7, MPI_COMM_WORLD, req(8), ierr) 
            if (ierr .ne. 0) then
              call MPI_Error_string(ierr,string_error,string_len)
              print*,'---> MessPassFloats: 1.8 (send):',string_error
              call MPI_Finalize (ierr)
              stop
            endif          
            call MPI_Send  (ibuf_sndNW, Nptsmax, MPI_DOUBLE_PRECISION,
     &                         p_NW, 8, MPI_COMM_WORLD,         ierr)
            if (ierr .ne. 0) then
              call MPI_Error_string(ierr,string_error,string_len)
              print*,'---> MessPassFloats: 1.8 (send):',string_error
              call MPI_Finalize (ierr)
              stop
            endif 
          endif
        endif

      enddo   !<-- iter
!
! Wait for completion of receive and fill ghost points: sides...
!
      if (WEST_INTER) then
        call MPI_Wait (req(1),status(1,1),ierr)
        if (ierr .ne. 0) then
          call MPI_Error_string(ierr,string_error,string_len)
          print*,'---> MessPassFloats: 2.1:',string_error
          call MPI_Finalize (ierr)
          stop
        endif
        ifltmax=nint(ibuf_revW(1))        
        if (ifltmax .ge. 1) then
C           write(*,*) 'MessPassFloats: 2.1', mynode, p_W
          i0=1
          do k=1,ifltmax
            ! iflt
            iflt=nint(ibuf_revW(i0+1))
            i0=i0+1
            ! track
            do j=1,NFV
              do i=0,NFT
                track(j,i,iflt) = ibuf_revW(i0+j+i*NFV)
              enddo
            enddo
            i0=i0+NFV*(NFT+1)
            ! trackaux
            do j=1,NDIAGS
              trackaux(j,iflt) = ibuf_revW(i0+j)
            enddo
            i0=i0+NDIAGS
            ! trackinfo
            do j=1,3
            trackinfo(j,iflt) = nint(ibuf_revW(i0+j))
          enddo
          i0=i0+3
        enddo
       endif
      endif

      if (EAST_INTER) then
        call MPI_Wait (req(2),status(1,2),ierr)
        if (ierr .ne. 0) then
          call MPI_Error_string(ierr,string_error,string_len)
          print*,'---> MessPassFloats: 2.2:',string_error
          call MPI_Finalize (ierr)
          stop
        endif                
        ifltmax=nint(ibuf_revE(1))
        if (ifltmax .ge. 1) then
C           write(*,*) 'MessPassFloats: 2.2', mynode, p_E
          i0=1
          do k=1,ifltmax
            ! iflt
            iflt=nint(ibuf_revE(i0+1))
            i0=i0+1
            ! track
            do j=1,NFV
              do i=0,NFT
                track(j,i,iflt) = ibuf_revE(i0+j+i*NFV)
              enddo
            enddo
            i0=i0+NFV*(NFT+1)
            ! trackaux
            do j=1,NDIAGS
              trackaux(j,iflt) = ibuf_revE(i0+j)
            enddo
            i0=i0+NDIAGS
            ! trackinfo
            do j=1,3
              trackinfo(j,iflt) = nint(ibuf_revE(i0+j))
            enddo
          enddo
        endif
      endif

      if (SOUTH_INTER) then
        call MPI_Wait (req(3),status(1,3),ierr) 
        if (ierr .ne. 0) then
          call MPI_Error_string(ierr,string_error,string_len)
          print*,'---> MessPassFloats: 2.3:',string_error
          call MPI_Finalize (ierr)
          stop
        endif     
        ifltmax=nint(ibuf_revS(1))
        if (ifltmax .ge. 1) then
C           write(*,*) 'MessPassFloats: 2.3', mynode, p_S
          i0=1
          do k=1,ifltmax
            ! iflt
            iflt=nint(ibuf_revS(i0+1))
            i0=i0+1
            ! track
            do j=1,NFV
              do i=0,NFT
                track(j,i,iflt) = ibuf_revS(i0+j+i*NFV)
              enddo
            enddo
            i0=i0+NFV*(NFT+1)
            ! trackaux
            do j=1,NDIAGS
              trackaux(j,iflt) = ibuf_revS(i0+j)
            enddo
            i0=i0+NDIAGS
            ! trackinfo
            do j=1,3
              trackinfo(j,iflt) = nint(ibuf_revS(i0+j))
            enddo
            i0=i0+3
          enddo
        endif
      endif

      if (NORTH_INTER) then
        call MPI_Wait (req(4),status(1,4),ierr) 
        if (ierr .ne. 0) then
          call MPI_Error_string(ierr,string_error,string_len)
          print*,'---> MessPassFloats: 2.4:',string_error
          call MPI_Finalize (ierr)
          stop
        endif            
        ifltmax=nint(ibuf_revN(1))
        if (ifltmax .ge. 1) then
C           write(*,*) 'MessPassFloats: 2.4', mynode, p_N
          i0=1
          do k=1,ifltmax
            ! iflt
            iflt=nint(ibuf_revN(i0+1))
            i0=i0+1
            ! track
            do j=1,NFV
              do i=0,NFT
                track(j,i,iflt) = ibuf_revN(i0+j+i*NFV)
              enddo
            enddo
            i0=i0+NFV*(NFT+1)
            ! trackaux
            do j=1,NDIAGS
              trackaux(j,iflt) = ibuf_revN(i0+j)
            enddo
            i0=i0+NDIAGS
            ! trackinfo
            do j=1,3
              trackinfo(j,iflt) = nint(ibuf_revN(i0+j))
            enddo
            i0=i0+3
          enddo
        endif
      endif

!
! ...corners:
!
      if (SOUTH_INTER .and. WEST_INTER) then
        call MPI_Wait (req(5),status(1,5),ierr) 
        if (ierr .ne. 0) then
          call MPI_Error_string(ierr,string_error,string_len)
          print*,'---> MessPassFloats: 2.5:',string_error
          call MPI_Finalize (ierr)
          stop
        endif            
        ifltmax=nint(ibuf_revSW(1))
        if (ifltmax .ge. 1) then
C           write(*,*) 'MessPassFloats: 2.5', mynode, p_SW
          i0=1
          do k=1,ifltmax
            ! iflt
            iflt=nint(ibuf_revSW(i0+1))
            i0=i0+1
            ! track
            do j=1,NFV
              do i=0,NFT
                track(j,i,iflt) = ibuf_revSW(i0+j+i*NFV)
              enddo
            enddo
            i0=i0+NFV*(NFT+1)
            ! trackaux
            do j=1,NDIAGS
              trackaux(j,iflt) = ibuf_revSW(i0+j)
            enddo
            i0=i0+NDIAGS
            ! trackinfo
            do j=1,3
              trackinfo(j,iflt) = nint(ibuf_revSW(i0+j))
            enddo
            i0=i0+3
          enddo
        endif
      endif

      if (NORTH_INTER .and. EAST_INTER) then
        call MPI_Wait (req(6),status(1,6),ierr) 
        if (ierr .ne. 0) then
          call MPI_Error_string(ierr,string_error,string_len)
          print*,'---> MessPassFloats: 2.6:',string_error
          call MPI_Finalize (ierr)
          stop
        endif            
        ifltmax=nint(ibuf_revNE(1))
        if (ifltmax .ge. 1) then
C           write(*,*) 'MessPassFloats: 2.6', mynode, p_NE
          i0=1
          do k=1,ifltmax
            ! iflt
            iflt=nint(ibuf_revNE(i0+1))
            i0=i0+1
            ! track
            do j=1,NFV
              do i=0,NFT
                track(j,i,iflt) = ibuf_revNE(i0+j+i*NFV)
              enddo
            enddo
            i0=i0+NFV*(NFT+1)
            ! trackaux
            do j=1,NDIAGS
              trackaux(j,iflt) = ibuf_revNE(i0+j)
            enddo
            i0=i0+NDIAGS
            ! trackinfo
            do j=1,3
              trackinfo(j,iflt) = nint(ibuf_revNE(i0+j))
            enddo
            i0=i0+3
          enddo
        endif
      endif

      if (SOUTH_INTER .and. EAST_INTER) then
        call MPI_Wait (req(7),status(1,7),ierr) 
        if (ierr .ne. 0) then
          call MPI_Error_string(ierr,string_error,string_len)
          print*,'---> MessPassFloats: 2.7:',string_error
          call MPI_Finalize (ierr)
          stop
        endif            
        ifltmax=nint(ibuf_revSE(1))
        if (ifltmax .ge. 1) then
C           write(*,*) 'MessPassFloats: 2.7', mynode, p_SE
          i0=1
          do k=1,ifltmax
            ! iflt
            iflt=nint(ibuf_revSE(i0+1))
            i0=i0+1
            ! track
            do j=1,NFV
              do i=0,NFT
                track(j,i,iflt) = ibuf_revSE(i0+j+i*NFV)
              enddo
            enddo
            i0=i0+NFV*(NFT+1)
            ! trackaux
            do j=1,NDIAGS
              trackaux(j,iflt) = ibuf_revSE(i0+j)
            enddo
            i0=i0+NDIAGS
            ! trackinfo
            do j=1,3
              trackinfo(j,iflt) = nint(ibuf_revSE(i0+j))
            enddo
            i0=i0+3
          enddo
        endif
      endif

      if (NORTH_INTER .and. WEST_INTER) then
        call MPI_Wait (req(8),status(1,8),ierr) 
        if (ierr .ne. 0) then
          call MPI_Error_string(ierr,string_error,string_len)
          print*,'---> MessPassFloats: 2.8:',string_error
          call MPI_Finalize (ierr)
          stop
        endif            
        ifltmax=nint(ibuf_revNW(1))
        if (ifltmax .ge. 1) then
C           write(*,*) 'MessPassFloats: 2.8', mynode, p_NW
          i0=1
          do k=1,ifltmax
            ! iflt
            iflt=nint(ibuf_revNW(i0+1))
            i0=i0+1
            ! track
            do j=1,NFV
              do i=0,NFT
                track(j,i,iflt) = ibuf_revNW(i0+j+i*NFV)
              enddo
            enddo
            i0=i0+NFV*(NFT+1)
            ! trackaux
            do j=1,NDIAGS
              trackaux(j,iflt) = ibuf_revNW(i0+j)
            enddo
            i0=i0+NDIAGS
            ! trackinfo
            do j=1,3
              trackinfo(j,iflt) = nint(ibuf_revNW(i0+j))
            enddo
            i0=i0+3
          enddo
        endif
      endif

      call MPI_Barrier(MPI_COMM_WORLD,ierr)

#else
      subroutine exchange_floats_empty
#endif
! aponte flt )
      return
      end

