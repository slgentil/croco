      ini_jet_name = "ini_jet_slope.nc"
      lstr=lenstr(ini_jet_name)
      ierr=nf_open(ini_jet_name(1:lstr), nf_nowrite, ncid)

      imin=iminmpi-3+3
      imax=imaxmpi+3+3
      jmin=jminmpi-3+3
      jmax=jmaxmpi+3+3
      icount=imax-imin+1
      jcount=jmax-jmin+1

!
! Topography. (h)
!
      ierr=nf_inq_varid (ncid, "h", varid)
      if (ierr .eq. nf_noerr) then
        ierr=nf_fread (h(START_2D_ARRAY), ncid, varid,
     &                                         0, r2dvar)
C         ierr = nf_get_vars_double(ncid,varid,(/imin,jmin/),
C      &        (/icount,jcount/), (/1,1/),h)
        if (ierr .ne. nf_noerr) then
          MPI_master_only write(stdout,4) "h "
     &                                ,ini_jet_name
          stop                                         !--> ERROR
        endif
      else
        MPI_master_only  write(stdout,3) "h"
     &                                ,ini_jet_name
        stop                                           !--> ERROR
      endif

  3   format(/1x,'ini_jet_slope - unable to find variable:',    1x,A,
     &                            /15x,'in input NetCDF file:',1x,A/)
  4   format(/1x,'ini_jet_slope - error while reading variable:',1x, A,
     &                            /15x,'in input NetCDF file:',1x,A/)
  
!
!  Close input NetCDF file.
!
       ierr=nf_close(ncid)

