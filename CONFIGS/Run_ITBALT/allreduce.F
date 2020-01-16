
       real function maxofall(maxperproc,communicator)

       implicit none
       include "mpif.h"
      
!      get the max of all proc

       integer :: code, communicator
       real :: maxperproc,result
      
       call mpi_allreduce(maxperproc,result,1,MPI_DOUBLE_PRECISION,
     &                   MPI_MAX,communicator,code)
       maxofall = result
      
       end function maxofall

!------------------------------------------------
  
       real function minofall(minperproc,communicator)

       implicit none
       include "mpif.h"
      
!      get the max of all proc

       integer :: code, communicator
       real :: minperproc,result
      
       call mpi_allreduce(minperproc,result,1,MPI_DOUBLE_PRECISION,
     &                   MPI_MIN,communicator,code)
       minofall = result
      
       end function minofall
  
!------------------------------------------------
    
       real function sumofall(sumperproc,communicator)

       implicit none
       include "mpif.h"
      
!      get the max of all proc

       integer :: code, communicator
       real :: sumperproc,result
      
       call mpi_allreduce(sumperproc,result,1,MPI_DOUBLE_PRECISION,
     &                   MPI_SUM,communicator,code)
       sumofall = result
      
       end function sumofall
