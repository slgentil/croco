! $Id: ANA_GET_DEPTH.F 2015-06-10 slgentil $
!
!
!aponte bathy(
#include "cppdefs.h"
                              
      subroutine ana_get_depth  ! Read depth from a NetCDF initialization file.
                               
      implicit none
# include "param.h"
# include "scalars.h"
# include "ncscrum.h"
# include "grid.h"

      integer ncid, indx, varid,  ierr, lstr, lenstr,
     &        nf_fread, checkdims
      character bathy_name*64
      
#include "netcdf.inc"
!
! Open bathymetry netCDF file for reading. Check that all
! spatial dimensions in that file are consistent with the model arrays,
! and read the depth data.
!
      bathy_name = 'bathy.nc'
      lstr=lenstr(bathy_name)
# if defined MPI && defined PARALLEL_FILES
      call insert_node (bathy_name, lstr, mynode, NNODES, ierr)
# endif
      ierr=nf_open(bathy_name(1:lstr), nf_nowrite, ncid)
      if (ierr .eq. nf_noerr) then
        ierr=checkdims (ncid, ininame, lstr, indx)
        if (ierr. ne. nf_noerr) then
          goto 99
        endif
      else
        write(stdout,'(/1x,2A/15x,3A)') 'ANA_GET_DEPTH ERROR: Cannot ',
     &               'open netCDF file', '''', bathy_name(1:lstr) ,'''.'
        goto 99                                           !--> ERROR
      endif
!
! Read depth variable:
! ---- -- -------- ----- ----------
!
      ierr=nf_inq_varid (ncid, 'depth', varid)
      if (ierr .eq. nf_noerr) then
        ierr=nf_fread (h(START_2D_ARRAY), ncid, varid, indx, r2dvar)
        if (ierr .ne. nf_noerr) then
          write(stdout,2) 'depth',bathy_name(1:lstr)
          goto 99                                         !--> ERROR
        endif
      else
        write(stdout,1) 'depth', bathy_name(1:lstr)
        goto 99                                           !--> ERROR
      endif

!
!  Close input NetCDF file.
!
      ierr=nf_close(ncid)

  1   format(/1x,'ANA_GET_DEPTH - unable to find variable:',    1x,A,
     &                            /15x,'in input NetCDF file:',1x,A/)
  2   format(/1x,'ANA_GET_DEPTH - error while reading variable:',1x, A,
     &    2x,'in input NetCDF file:',1x,A/)
  3   format(/1x,'ANA_GET_DEPTH - unable to find variable:',    1x,A,
     &                            /15x,'in input NetCDF file:',1x,A,
     &    1x,'-> analytical value'/)
     
 99   return
 
!aponte bathy)
      end
