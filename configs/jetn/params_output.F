#include "cppdefs.h"

      module params

          implicit none

          interface params_output
          module procedure params_output_scalar, params_output_array
          end interface params_output

          contains

          subroutine params_output_scalar(varname,vardata)

# include "param.h"
# include "scalars.h"

              character(len=*), intent(in) :: varname
              class(*), intent(in) :: vardata
              select type(vardata)
              type is (real)
                  MPI_master_only  write(*,'(A,A,A,2X,E17.10)') OSI, 
     &                             trim(varname), '=' , vardata 
              type is (integer)
                  MPI_master_only  write(*,'(A,A,A,2X,I2)') OSI, 
     &                             trim(varname), '='  , vardata
              end select
          end subroutine params_output_scalar

          subroutine params_output_array(varname,vardata)

# include "param.h"
# include "scalars.h"

              character(len=*), intent(in) :: varname
              class(*), dimension(:), intent(in) :: vardata
              integer :: i, length
              character(len=30) :: my_format


              length = size(vardata)
              select type(vardata)
              type is (real)
                  write(my_format, '(A, i0, A)') '(A,A,A,2X,', length,
     &                                           'E17.10,2X,A)'
                  MPI_master_only write(*,my_format) OSI,trim(varname),
     &                                           '=(' ,vardata,')'
              type is (integer)
                  write(my_format, '(A, i0, A)') '(A,A,A,2X,', length, 
     &                                            'I10,2X,A)'
                  MPI_master_only write(*,my_format) OSI,trim(varname),
     &                                           '=(' ,vardata,')'
              end select
          end subroutine params_output_array

      end module params
