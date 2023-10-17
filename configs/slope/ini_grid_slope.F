      ini_jet_name = "ini_jet_slope.nc"
      lstr=lenstr(ini_jet_name)
      ierr=nf_open(ini_jet_name(1:lstr), nf_nowrite, ncid)

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
!
! Size of the domain in the xi direction : xl
!
      ierr=nf_inq_varid (ncid, "xl", varid)
      if (ierr .eq. nf_noerr) then
         ierr = nf_get_var1_double(ncid, varid, (/1/), Length_XI)
         if (ierr .ne. nf_noerr) then
            MPI_master_only write(stdout,2) "xl ",ini_jet_name
          stop                                         !--> ERROR
        endif
      else
        MPI_master_only  write(stdout,1) "xl",ini_jet_name
        stop                                           !--> ERROR
      endif
!
! Size of the domain in the eta direction : el
!
      ierr=nf_inq_varid (ncid, "el", varid)
      if (ierr .eq. nf_noerr) then
         ierr = nf_get_var1_double(ncid, varid, (/1/), Length_ETA)
         if (ierr .ne. nf_noerr) then
            MPI_master_only write(stdout,2) "el ",ini_jet_name
          stop                                         !--> ERROR
        endif
      else
        MPI_master_only  write(stdout,1) "el",ini_jet_name
        stop                                           !--> ERROR
      endif
!
!  Max topography depth
!
      ierr=nf_inq_varid (ncid, "depth", varid)
      if (ierr .eq. nf_noerr) then
         ierr = nf_get_var1_double(ncid, varid, (/1/), depth)
         if (ierr .ne. nf_noerr) then
            MPI_master_only write(stdout,2) "depth ",ini_jet_name
          stop                                         !--> ERROR
        endif
      else
        MPI_master_only  write(stdout,1) "depth",ini_jet_name
        stop                                           !--> ERROR
      endif
!
!
!  Coriolis parameter f0
!
      ierr=nf_inq_varid (ncid, "f0", varid)
      if (ierr .eq. nf_noerr) then
         ierr = nf_get_var1_double(ncid, varid, (/1/), f0)
         if (ierr .ne. nf_noerr) then
            MPI_master_only write(stdout,2) "f0 ",ini_jet_name
          stop                                         !--> ERROR
        endif
      else
        MPI_master_only  write(stdout,1) "f0",ini_jet_name
        stop                                           !--> ERROR
      endif
!
! Coriolis parameter beta
!
      ierr=nf_inq_varid (ncid, "beta", varid)
      if (ierr .eq. nf_noerr) then
         ierr = nf_get_var1_double(ncid, varid, (/1/), beta)
         if (ierr .ne. nf_noerr) then
            MPI_master_only write(stdout,2) "beta ",ini_jet_name
          stop                                         !--> ERROR
        endif
      else
        MPI_master_only  write(stdout,1) "beta",ini_jet_name
        stop                                           !--> ERROR
      endif

  1   format(/1x,'ini_jet_slope - unable to find variable:',    1x,A,
     &                            /15x,'in input NetCDF file:',1x,A/)
  2   format(/1x,'ini_jet_slope - error while reading variable:',1x, A,
     &                            /15x,'in input NetCDF file:',1x,A/)
  
!  Output 

      MPI_master_only  write(stdout,*) " "
      call params_output('theta_s',theta_s)
      call params_output('theta_b',theta_b)
      call params_output('hc',hc)
      call params_output('Length_XI',Length_XI)
      call params_output('Length_ETA',Length_ETA)
      call params_output('depth',depth)
      call params_output('rho0',rho0)
      call params_output('f0',f0)
      call params_output('beta',beta)
      MPI_master_only  write(stdout,*) " "


!
!  Close input NetCDF file.
!
       ierr=nf_close(ncid)

