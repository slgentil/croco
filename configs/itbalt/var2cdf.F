        subroutine var2cdf(array,nbdims,shape_array,filename)

        implicit none
        real(kind=8), intent(IN) :: array(*)
        character*(*),intent(IN) :: filename 
        integer, intent(IN) :: nbdims
        integer, intent(IN) :: shape_array(*)

        integer :: file_id, xdim_id, ydim_id
        integer :: array_id
        integer, dimension(:), allocatable :: arrdims
        !    character(len=*) :: argtitle = Flag_in

        integer :: i, j
        integer :: ierr
        character*4 :: string 

#include "netcdf.inc"

        allocate(arrdims(nbdims))

        ! create the file
        ierr = nf_create(filename, NF_CLOBBER, file_id)

        ! define the dimensions
        do i=1,nbdims
            write(string,'(A,I1)') 'dim',i
            ierr = nf_def_dim(file_id, trim(string), shape_array(i), 
     &                       xdim_id)
            arrdims(i) = xdim_id
        enddo

        ! definethe variable
        ierr = nf_def_var(file_id, 'array',  NF_DOUBLE, nbdims, arrdims, 
     &                     array_id)
        ! done defining
        ierr = nf_enddef(file_id)

        ! Write out the values
        ierr = nf_put_var(file_id, array_id, array)

        ! close; done
        ierr = nf_close(file_id)
        return
        end 