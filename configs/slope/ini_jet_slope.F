C # ifdef THREE_GHOST_POINTS
C #   undef THREE_GHOST_POINTS
C #   define SLOPE_TGH
C # endif

# ifdef PARALLEL_FILES
#   undef PARALLEL_FILES
#   define  SLOPE_PF 
# endif

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
! Free-surface. (zeta)
!
      lvar=lenstr(vname(1,indxZ))
      ierr=nf_inq_varid (ncid, vname(1,indxZ)(1:lvar), varid)
      if (ierr .eq. nf_noerr) then
        ierr=nf_fread (zeta(START_2D_ARRAY,1), ncid, varid,
     &                                         indx, r2dvar)
C         ierr = nf_get_vars_double(ncid,varid,(/imin,jmin/),
C      &        (/icount,jcount/), (/1,1/),zeta(-2,-2,1))
        if (ierr .ne. nf_noerr) then
          MPI_master_only write(stdout,2) vname(1,indxZ)(1:lvar)
     &                                ,ininame(1:lstr)
          stop                                         !--> ERROR
        endif
      else
        MPI_master_only  write(stdout,1) vname(1,indxZ)(1:lvar)
     &                                ,ininame(1:lstr)
        stop                                           !--> ERROR
      endif
!
! 2D momentum component in the XI-direction. (ubar)
!
      lvar=lenstr(vname(1,indxUb))
      ierr=nf_inq_varid (ncid, vname(1,indxUb)(1:lvar), varid)
      if (ierr .eq. nf_noerr) then
        ierr=nf_fread (ubar(START_2D_ARRAY,1), ncid, varid,
     &                                         indx, u2dvar)
C         ierr = nf_get_vars_double(ncid,varid,(/imin,jmin/),
C      &        (/icount,jcount/), (/1,1/),ubar(-2,-2,1))
        if (ierr .ne. nf_noerr) then
          MPI_master_only write(stdout,2) vname(1,indxUb)(1:lvar)
     &                                , ininame(1:lstr)
          stop                                        !--> ERROR
        endif
      else
        MPI_master_only write(stdout,1) vname(1,indxUb)(1:lvar)
     &                                ,ininame(1:lstr)
        stop                                           !--> ERROR
      endif
!
! 2D momentum component in the ETA-direction (vbar)
!
      lvar=lenstr(vname(1,indxVb))
      ierr=nf_inq_varid (ncid, vname(1,indxVb)(1:lvar), varid)
      if (ierr .eq. nf_noerr) then
        ierr=nf_fread (vbar(START_2D_ARRAY,1), ncid, varid,
     &                                         indx, v2dvar)
C         ierr = nf_get_vars_double(ncid,varid,(/imin,jmin/),
C      &        (/icount,jcount/), (/1,1/),vbar(-2,-2,1))
        if (ierr .ne. nf_noerr) then
          MPI_master_only write(stdout,2) vname(1,indxVb)(1:lvar)
     &                                , ininame(1:lstr)
          stop                                          !--> ERROR
        endif
      else
        MPI_master_only write(stdout,1) vname(1,indxVb)(1:lvar)
     &                                ,ininame(1:lstr)
        stop                                           !--> ERROR
      endif
!
! 3D momentum component in the XI-direction. (u)
!
      lvar=lenstr(vname(1,indxU))
      ierr=nf_inq_varid (ncid, vname(1,indxU)(1:lvar), varid)
      if (ierr .eq. nf_noerr) then
        ierr=nf_fread (u(START_2D_ARRAY,1,1), ncid, varid,
     &                                        indx, u3dvar)
C         ierr = nf_get_vars_double(ncid,varid,(/imin,jmin,1/),
C      &        (/icount,jcount,N/), (/1,1,1/),u(-2,-2,1,1))
        if (ierr .ne. nf_noerr) then
          MPI_master_only write(stdout,2) vname(1,indxU)(1:lvar)
     &                                , ininame(1:lstr)
          stop                                         !--> ERROR
        endif
      else
        MPI_master_only write(stdout,1) vname(1,indxU)(1:lvar)
     &                                ,ininame(1:lstr)
        stop                                           !--> ERROR
      endif
!
! 3D momentum component in the ETA-direction. (v)
!
      lvar=lenstr(vname(1,indxV))
      ierr=nf_inq_varid (ncid, vname(1,indxV)(1:lvar), varid)
      if (ierr .eq. nf_noerr) then
        ierr=nf_fread (v(START_2D_ARRAY,1,1), ncid, varid,
     &                                        indx, v3dvar)
C         ierr = nf_get_vars_double(ncid,varid,(/imin,jmin,1/),
C      &        (/icount,jcount,N/), (/1,1,1/),v(-2,-2,1,1))
        if (ierr .ne. nf_noerr) then
          MPI_master_only write(stdout,2) vname(1,indxV)(1:lvar)
     &                                , ininame(1:lstr)
          stop                                        !--> ERROR
        endif
      else
        MPI_master_only write(stdout,1) vname(1,indxV)(1:lvar)
     &                                ,ininame(1:lstr)
        stop                                           !--> ERROR
      endif
!
! Tracer variables. (temp, salt)
!
      do itrc=1,NT
        lvar=lenstr(vname(1,indxT+itrc-1))
        ierr=nf_inq_varid (ncid, vname(1,indxT+itrc-1)(1:lvar),
     &                                                    varid)
        if (ierr .eq. nf_noerr) then
          ierr=nf_fread (t(START_2D_ARRAY,1,1,itrc), ncid,  varid,
     &                                               indx, r3dvar)
C         ierr = nf_get_vars_double(ncid,varid,(/imin,jmin,1/),
C      &        (/icount,jcount,N/), (/1,1,1/),t(-2,-2,1,1,itrc))

          if (ierr .ne. nf_noerr) then
        MPI_master_only write(stdout,2) vname(1,indxT+itrc-1)(1:lvar)
     &                                , ininame(1:lstr)
            stop                                       !--> ERROR
          endif
        else
        MPI_master_only write(stdout,1) vname(1,indxT+itrc-1)(1:lvar),
     &                                    ininame(1:lstr)
        endif
      enddo


  1   format(/1x,'ini_jet_slope - unable to find variable:',    1x,A,
     &                            /15x,'in input NetCDF file:',1x,A/)
  2   format(/1x,'ini_jet_slope - error while reading variable:',1x, A,
     &                            /15x,'in input NetCDF file:',1x,A/)


!
!  Close input NetCDF file.
!
      ierr=nf_close(ncid)

C # ifdef SLOPE_TGH
C #   define THREE_GHOST_POINTS
C # endif

# ifdef SLOPE_PF
#   define PARALLEL_FILES
# endif
