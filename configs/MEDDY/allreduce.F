#include "cppdefs.h"

       real(8) function nr_maxofall(maxperproc)

       implicit none
#include "param.h"
#include "mpi_roms.h"

       include "mpif.h"
  
! get the max of all proc

       integer :: code
       real(8) :: maxperproc,result
  
       call mpi_allreduce(maxperproc,result,1,MPI_DOUBLE_PRECISION,MPI_MAX,
     &                    MPI_COMM_WORLD,code)
       nr_maxofall = result
  
       end function nr_maxofall
       
!-------------------------------------------------------------------------------
       real(8) function nr_minofall(minperproc)

       implicit none
#include "mpi_roms.h"
       include "mpif.h"
  
! get the max of all proc

       integer :: code
       real(8) :: minperproc,result
  
       call mpi_allreduce(minperproc,result,1,MPI_DOUBLE_PRECISION,MPI_MIN,
     &                    MPI_COMM_WORLD,code)
       nr_minofall = result
  
       end function nr_minofall
  
!------------------------------------------------------------------------------- 
       real(8) function nr_sumofall(sumperproc)

       implicit none
#include "mpi_roms.h"      
       include "mpif.h"
  
!      get the max of all proc

       integer :: code
       real(8) :: sumperproc,result
  
       call mpi_allreduce(sumperproc,result,1,MPI_DOUBLE_PRECISION,MPI_SUM,
     &                    MPI_COMM_WORLD,code)
       nr_sumofall = result
  
       end function nr_sumofall
