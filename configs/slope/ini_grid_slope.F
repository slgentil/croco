      ini_jet_name = "ini_jet_slope.nc"
      lstr=lenstr(ini_jet_name)
      ierr=nf_open(ini_jet_name(1:lstr), nf_nowrite, ncid)

      imin=iminmpi-3+3
      imax=imaxmpi+3+3
      jmin=jminmpi-3+3
      jmax=jmaxmpi+3+3
      icount=imax-imin+1
      jcount=jmax-jmin+1
C       print*,'mynode,imin,imax,jmin,jmax=',mynode,imin,imax,jmin,jmax

!
! Vertical S-coordinate surface stretching parameter: theta_s
!
      ierr=nf_inq_varid (ncid, "theta_s", varid)
      if (ierr .eq. nf_noerr) then
         ierr = nf_get_var1_double(ncid, varid, (/1/), theta_s)
         if (ierr .ne. nf_noerr) then
            MPI_master_only write(stdout,2) "theta_s ",ini_jet_name
          stop                                         !--> ERROR
        endif
      else
        MPI_master_only  write(stdout,1) "theta_s",ini_jet_name
        stop                                           !--> ERROR
      endif
      MPI_master_only  write(stdout,*) "theta_s done"
!
! Vertical S-coordinate bottom stretching parameter: theta_b
!
      ierr=nf_inq_varid (ncid, "theta_b", varid)
      if (ierr .eq. nf_noerr) then
         ierr = nf_get_var1_double(ncid, varid, (/1/), theta_b)
         if (ierr .ne. nf_noerr) then
            MPI_master_only write(stdout,2) "theta_b ",ini_jet_name
          stop                                         !--> ERROR
        endif
      else
        MPI_master_only  write(stdout,1) "theta_b",ini_jet_name
        stop                                           !--> ERROR
      endif
      MPI_master_only  write(stdout,*) "theta_b done"
!
! Vertical S-coordinate Hc parameter: hc
!
      ierr=nf_inq_varid (ncid, "hc", varid)
      if (ierr .eq. nf_noerr) then
         ierr = nf_get_var1_double(ncid, varid, (/1/), hc)
         if (ierr .ne. nf_noerr) then
            MPI_master_only write(stdout,2) "hc ",ini_jet_name
          stop                                         !--> ERROR
        endif
      else
        MPI_master_only  write(stdout,1) "hc",ini_jet_name
        stop                                           !--> ERROR
      endif
      MPI_master_only  write(stdout,*) "hc done"
!
! background density : rho0
!
      ierr=nf_inq_varid (ncid, "rho0", varid)
      if (ierr .eq. nf_noerr) then
         ierr = nf_get_var1_double(ncid, varid, (/1/), rho0)
         if (ierr .ne. nf_noerr) then
            MPI_master_only write(stdout,2) "rho0 ",ini_jet_name
          stop                                         !--> ERROR
        endif
      else
        MPI_master_only  write(stdout,1) "rho0",ini_jet_name
        stop                                           !--> ERROR
      endif
      MPI_master_only  write(stdout,*) "rho0 done"
!
! Size of the domain in the xi direction : xl
!
      ierr=nf_inq_varid (ncid, "xl", varid)
      if (ierr .eq. nf_noerr) then
         ierr = nf_get_var1_double(ncid, varid, (/1/), xl)
         if (ierr .ne. nf_noerr) then
            MPI_master_only write(stdout,2) "xl ",ini_jet_name
          stop                                         !--> ERROR
        endif
      else
        MPI_master_only  write(stdout,1) "xl",ini_jet_name
        stop                                           !--> ERROR
      endif
      MPI_master_only  write(stdout,*) "xl done"
!
! Size of the domain in the eta direction : el
!
      ierr=nf_inq_varid (ncid, "el", varid)
      if (ierr .eq. nf_noerr) then
         ierr = nf_get_var1_double(ncid, varid, (/1/), el)
         if (ierr .ne. nf_noerr) then
            MPI_master_only write(stdout,2) "el ",ini_jet_name
          stop                                         !--> ERROR
        endif
      else
        MPI_master_only  write(stdout,1) "el",ini_jet_name
        stop                                           !--> ERROR
      endif
      MPI_master_only  write(stdout,*) "el done"
!
! Size of a horizontal cell in the xi direction : dx
!
      ierr=nf_inq_varid (ncid, "dx", varid)
      if (ierr .eq. nf_noerr) then
         ierr = nf_get_var1_double(ncid, varid, (/1/), dx)
         if (ierr .ne. nf_noerr) then
            MPI_master_only write(stdout,2) "dx ",ini_jet_name
          stop                                         !--> ERROR
        endif
      else
        MPI_master_only  write(stdout,1) "dx",ini_jet_name
        stop                                           !--> ERROR
      endif
      MPI_master_only  write(stdout,*) "dx done"
!
! Size of a horizontal cell in the eta direction : dy
!
      ierr=nf_inq_varid (ncid, "dy", varid)
      if (ierr .eq. nf_noerr) then
         ierr = nf_get_var1_double(ncid, varid, (/1/), dy)
         if (ierr .ne. nf_noerr) then
            MPI_master_only write(stdout,2) "dy ",ini_jet_name
          stop                                         !--> ERROR
        endif
      else
        MPI_master_only  write(stdout,1) "dy",ini_jet_name
        stop                                           !--> ERROR
      endif
      MPI_master_only  write(stdout,*) "dy done"
!
! x at rho point. (xr)
!
      ierr=nf_inq_varid (ncid, "xr", varid)
      if (ierr .eq. nf_noerr) then
        ierr=nf_fread (xr(START_2D_ARRAY), ncid, varid,
     &                                         0, r2dvar)
C         ierr = nf_get_vars_double(ncid,varid,(/imin,jmin/),
C      &        (/icount,jcount/), (/1,1/),xr)
        if (ierr .ne. nf_noerr) then
          MPI_master_only write(stdout,2) "xr ",ini_jet_name
          stop                                         !--> ERROR
        endif
      else
        MPI_master_only  write(stdout,1) "xr",ini_jet_name
        stop                                           !--> ERROR
      endif
      MPI_master_only  write(stdout,*) "xr done"

!
! y at rho point. (yr)
!
      ierr=nf_inq_varid (ncid, "yr", varid)
      if (ierr .eq. nf_noerr) then
        ierr=nf_fread (yr(START_2D_ARRAY), ncid, varid,
     &                                         0, r2dvar)
C         ierr = nf_get_vars_double(ncid,varid,(/imin,jmin/),
C      &        (/icount,jcount/), (/1,1/),yr)
        if (ierr .ne. nf_noerr) then
          MPI_master_only write(stdout,2) "yr ",ini_jet_name
          stop                                         !--> ERROR
        endif
      else
        MPI_master_only  write(stdout,1) "yr",ini_jet_name
        stop                                           !--> ERROR
      endif
      MPI_master_only  write(stdout,*) "yr done"

!
! x at psi point. (xp)
!
      ierr=nf_inq_varid (ncid, "xp", varid)
      if (ierr .eq. nf_noerr) then
        ierr=nf_fread (xp(START_2D_ARRAY), ncid, varid,
     &                                         0, p2dvar)
C         ierr = nf_get_vars_double(ncid,varid,(/imin,jmin/),
C      &        (/icount,jcount/), (/1,1/),xp)
        if (ierr .ne. nf_noerr) then
          MPI_master_only write(stdout,2) "xp ",ini_jet_name
          stop                                         !--> ERROR
        endif
      else
        MPI_master_only  write(stdout,1) "xp",ini_jet_name
        stop                                           !--> ERROR
      endif
      MPI_master_only  write(stdout,*) "xp done"

!
! y at psi point. (yp)
!
      ierr=nf_inq_varid (ncid, "yp", varid)
      if (ierr .eq. nf_noerr) then
        ierr=nf_fread (yp(START_2D_ARRAY), ncid, varid,
     &                                         0, p2dvar)
C         ierr = nf_get_vars_double(ncid,varid,(/imin,jmin/),
C      &        (/icount,jcount/), (/1,1/),yp)
        if (ierr .ne. nf_noerr) then
          MPI_master_only write(stdout,2) "yp ",ini_jet_name
          stop                                         !--> ERROR
        endif
      else
        MPI_master_only  write(stdout,1) "yp",ini_jet_name
        stop                                           !--> ERROR
      endif
      MPI_master_only  write(stdout,*) "yp done"

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
          MPI_master_only write(stdout,2) "h "
     &                                ,ini_jet_name
          stop                                         !--> ERROR
        endif
      else
        MPI_master_only  write(stdout,1) "h"
     &                                ,ini_jet_name
        stop                                           !--> ERROR
      endif
      MPI_master_only  write(stdout,*) "h done"
!
!  Coriolis parameter (f)
!
      ierr=nf_inq_varid (ncid, "f", varid)
      if (ierr .eq. nf_noerr) then
        ierr=nf_fread (f(START_2D_ARRAY), ncid, varid,
     &                                         0, r2dvar)
C         ierr = nf_get_vars_double(ncid,varid,(/imin,jmin/),
C      &        (/icount,jcount/), (/1,1/),f)
        if (ierr .ne. nf_noerr) then
          MPI_master_only write(stdout,2) "f "
     &                                , ini_jet_name
          stop                                        !--> ERROR
        endif
      else
        MPI_master_only write(stdout,1) "f "
     &                                ,ini_jet_name
        stop                                           !--> ERROR
      endif
      MPI_master_only  write(stdout,*) "f done"

  1   format(/1x,'ini_jet_slope - unable to find variable:',    1x,A,
     &                            /15x,'in input NetCDF file:',1x,A/)
  2   format(/1x,'ini_jet_slope - error while reading variable:',1x, A,
     &                            /15x,'in input NetCDF file:',1x,A/)
  

      MPI_master_only write(stdout,*) "theta_s,theta_b,hc=",theta_s,theta_b,hc
      MPI_master_only write(stdout,*) "rh0 = ",rho0
      MPI_master_only write(stdout,*) "xl,el = ",xl,el
      MPI_master_only write(stdout,*) "dx,dy = ",dx,dy
      MPI_master_only write(stdout,*) "xr = ",xr(-2:1,-2)
      MPI_master_only write(stdout,*) "xp = ",xp(-2:1,-2)
      MPI_master_only write(stdout,*) "yr = ",yr(-2,-2:1)
      MPI_master_only write(stdout,*) "yp = ",yp(-2,-2:1)
      MPI_master_only write(stdout,*) "f = ",f(-2:1,-2)


!
!  Close input NetCDF file.
!
       ierr=nf_close(ncid)

