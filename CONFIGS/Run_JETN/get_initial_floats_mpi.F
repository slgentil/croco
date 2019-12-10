! $Id: get_initial_floats.F 1458 2014-02-03 15:01:25Z gcambon $
!
!======================================================================
! ROMS_AGRIF is a branch of ROMS developped at IRD and INRIA, in France
! The two other branches from UCLA (Shchepetkin et al) 
! and Rutgers University (Arango et al) are under MIT/X style license.
! ROMS_AGRIF specific routines (nesting) are under CeCILL-C license.
! 
! ROMS_AGRIF website : http://www.romsagrif.org
!======================================================================
!
#include "cppdefs.h"
! aponte flt (
#if defined FLOATS && defined MPI
                                         ! Read initial conditions for 
      subroutine get_initial_floats_mpi  ! the floats variables from 
                                         ! ascii initialization file.
C       implicit none
# include "param.h"
# include "scalars.h"
# include "floats.h"
# include "mpi_cpl.h"

      integer i, j, k, iflt, nfltmax, ierr
! ascii rst file
      character*20 filename
      character*20 rowfmt
      integer buff0(1:3)
      real buff1(1:NFV,0:NFT)
      real buff2(1:NDIAGS)
      real buff3(inode:izgrd)
      real(kind=8) :: nr_sumofall 





     ! ascii initialization file
      write (filename, "(A10,I0.3)") "float.rst.", mynode
      open(unit=12, file= trim(filename), action="read"
     &            , status="old")

      read(12,'(I4)') nfltmax

      !debug
      write(*,*) "Reads ", filename, nfltmax
 
      do i=1,nfltmax

       read(12,2)
     &         iflt
     &         , (buff0(j),j=1,3)
     &         , (buff1(j,0),j=1,NFV)
     &         , (buff1(j,1),j=1,NFV)
     &         , (buff1(j,2),j=1,NFV)
     &         , (buff1(j,3),j=1,NFV)
     &         , (buff2(j),j=1,NDIAGS)
     &         , (buff3(j),j=inode,izgrd)

       ! put data into float variables
       do j=1,3
        trackinfo(j,iflt)=buff0(j)
!        write(*,*) filename, 'trackinfo', iflt,  trackinfo(j,iflt)
       enddo
       do j=1,NFV
        do k=0,NFT
         track(j,k,iflt)=buff1(j,k)
!         write(*,*) filename, 'track', iflt, j, k, buff1(j,k)
        enddo
       enddo
       do j=1,NDIAGS
        trackaux(j,iflt)=buff2(j)
!        write(*,*) filename, 'trackaux', iflt, j, buff2(j)
       enddo
       do j=inode,izgrd
        Tinfo(j,iflt)=buff3(j)
!        write(*,*) filename, 'Tinfo', iflt, j, buff3(j)
       enddo


      enddo

      ! close output file
      close(12)
      
      nfloats = sumofall(real(nfltmax),MPI_COMM_WORLD)
      if ( mynode .eq. 0 )
     &    write(stdout,'(A44,I4)') 
     &      'ii Number of float trajectories to compute:',nfloats
   

      ! ascii output format
!   2  format(I4,(1+NFV*4+NDIAGS)(1x,e16.10))
! 1 + 6*4 + 10 + 6 = 41
   2  format(I4,3(1x,I4),41(1x,e16.10))


      ! Need to broadcast information in order to turn off floats
      ! that will be launched

#else 

      subroutine get_initial_floats_mpi_empty
#endif  /* FLOATS */
! aponte flt )
      return
      end



