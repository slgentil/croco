#include "cppdefs.h"

      subroutine stop_run()

#ifdef XIOS
          USE xios          ! XIOS module
#endif

      implicit none

#ifdef MPI
      include 'mpif.h'
      integer :: ierr
#endif
# include "mpi_cpl.h"


#ifdef MPI
      call MPI_Barrier(MPI_COMM_WORLD, ierr)  ! XIOS

# if defined OA_COUPLING || defined OW_COUPLING
      call prism_terminate_proto(ierr)       ! Finalize OASIS3
# elif defined XIOS
                             ! MPI finalize done by XIOS, XIOS+OASIS not coded
      call xios_finalize()   ! end mpp communications with xios
#  if !defined AGRIF      
      call MPI_Finalize (ierr)
#  endif 
# else
      call MPI_Finalize (ierr)
# endif     
#endif

      stop
      end

