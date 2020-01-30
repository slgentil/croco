! $Id: ANA_GET_DENSITY.F 2015-06-15 slgentil $
!
!
!aponte density(
#include "cppdefs.h"
                              
      subroutine ana_get_density (Istr,Iend,Jstr,Jend)  
      
! Read depth from a NetCDF initialization file.
                               
      implicit none
# include "param.h"
# include "grid.h"
# include "scalars.h"
# include "ncscrum.h"
# include "ocean2d.h"
# include "ocean3d.h"
# include "climat.h"
# include "netcdf.inc"
# include "private_scratch.h"

      integer Istr,Iend,Jstr,Jend, i,j,k, itrc      
      integer trd, omp_get_thread_num
      integer ncid, varid,  ierr, lstr, nf_fread
      character density_name*64
      
!
! Open densitymetry netCDF file for reading. Check that all
! spatial dimensions in that file are consistent with the model arrays,
! and read the depth data.
!

      ierr = 0
      density_name = "inifield.nc"   
      lstr=len_trim(density_name)
#if defined MPI && defined PARALLEL_FILES

      call insert_node (density_name, lstr, mynode, NNODES, ierr)

#endif
      ierr=nf_open(density_name(1:lstr), nf_nowrite, ncid) 
      if (ierr .ne. nf_noerr) then
        write(stdout,'(/1x,2A/15x,3A)') 'ANA_GET_DENSITY ERROR: Cannot ',
     &               'open netCDF file', '''', density_name(1:lstr)
        goto 99                                           !--> ERROR
      endif

#ifdef CLIMATOLOGY       

! Read climatology variables, ssh,tclm,sclm,uclm,vclm

! Free-surface.

         ierr=nf_inq_varid (ncid, 'ssh', varid)
         if (ierr .eq. nf_noerr) then
           ierr=nf_fread (ssh(START_2D_ARRAY), ncid, varid,
     &                                         0, r2dvar)
           if (ierr .ne. nf_noerr) then
             write(stdout,2) 'ssh',density_name(1:lstr)
             goto 99                                         !--> ERROR
           endif
         else
           write(stdout,1) 'ssh', density_name(1:lstr)
           goto 99                                           !--> ERROR
         endif    
                 
!     read non perturbated forcing fields tclm

         ierr=nf_inq_varid (ncid, 'tclm', varid)
         if (ierr .eq. nf_noerr) then
           ierr=nf_fread (tclm(START_2D_ARRAY,1,1), ncid, varid, 0, r3dvar)
           if (ierr .ne. nf_noerr) then
             write(stdout,2) 'tclm',density_name(1:lstr)
             goto 99                                         !--> ERROR
           endif
         else
           write(stdout,1) 'tclm', density_name(1:lstr)
           goto 99                                           !--> ERROR
         endif
         
#ifdef SALINITY
!	sclm         

         ierr=nf_inq_varid (ncid, 'sclm', varid)
         if (ierr .eq. nf_noerr) then
           ierr=nf_fread (tclm(START_2D_ARRAY,1,2), ncid, varid, 0, r3dvar)
           if (ierr .ne. nf_noerr) then
             write(stdout,2) 'sclm',density_name(1:lstr)
             goto 99                                         !--> ERROR
           endif
         else
           write(stdout,1) 'sclm', density_name(1:lstr)
           goto 99                       !--> ERROR
         endif
#endif  /* SALINITY */
         
! 3D momentum component in the XI-direction.

         ierr=nf_inq_varid (ncid, 'uclm', varid)
         if (ierr .eq. nf_noerr) then
            ierr=nf_fread (uclm(START_2D_ARRAY,1), ncid, varid,
     &                                        0, u3dvar)
            if (ierr .ne. nf_noerr) then
               write(stdout,2) 'uclm', density_name(1:lstr)
               goto 99                                         !--> ERROR
            endif
         else
            write(stdout,1) 'uclm', density_name(1:lstr)
            goto 99                                           !--> ERROR
         endif
    
! 3D momentum component in the ETA-direction.

         ierr=nf_inq_varid (ncid, 'vclm', varid)
         if (ierr .eq. nf_noerr) then
            ierr=nf_fread (vclm(START_2D_ARRAY,1), ncid, varid,
     &                                        0, v3dvar)
            if (ierr .ne. nf_noerr) then
               write(stdout,2) 'vclm', density_name(1:lstr)
               goto 99                                         !--> ERROR
            endif
         else
            write(stdout,1) 'vclm', density_name(1:lstr)
            goto 99                                           !--> ERROR
         endif

#endif /* CLIMATOLOGY */

! Read initial fields, zeta,t,s,u,v
         
! Free-surface.

         ierr=nf_inq_varid (ncid, 'zeta', varid)
         if (ierr .eq. nf_noerr) then
           ierr=nf_fread (zeta(START_2D_ARRAY,1), ncid, varid,
     &                                         0, r2dvar)
           if (ierr .ne. nf_noerr) then
             write(stdout,2) 'zeta',density_name(1:lstr)
             goto 99                                         !--> ERROR
           endif
         else
           write(stdout,1) 'zeta', density_name(1:lstr)
           goto 99                                           !--> ERROR
         endif   
         
!     read perturbated initial fields t and s  
! use the second time index in order no to interfer with ana_initial.F
                     
         ierr=nf_inq_varid (ncid, 't', varid)
         if (ierr .eq. nf_noerr) then
           ierr=nf_fread (t(START_2D_ARRAY,1,3,1), ncid, varid, 0, r3dvar)
           if (ierr .ne. nf_noerr) then
             write(stdout,2) 't',density_name(1:lstr)
             goto 99                                         !--> ERROR
           endif
         else
           write(stdout,1) 't', density_name(1:lstr)
           goto 99                                           !--> ERROR
         endif
         
#ifdef SALINITY

         ierr=nf_inq_varid (ncid, 's', varid)
         if (ierr .eq. nf_noerr) then
           ierr=nf_fread (t(START_2D_ARRAY,1,3,2), ncid, varid, 0, r3dvar)
           if (ierr .ne. nf_noerr) then
             write(stdout,2) 's',density_name(1:lstr)
             goto 99                                         !--> ERROR
           endif
         else
           write(stdout,1) 's', density_name(1:lstr)
           goto 99                                           !--> ERROR
         endif
#endif   /* SALINITY */ 
 
! 3D momentum component in the XI-direction.

         ierr=nf_inq_varid (ncid, 'u', varid)
         if (ierr .eq. nf_noerr) then
            ierr=nf_fread (u(START_2D_ARRAY,1,1), ncid, varid,
     &                                        0, u3dvar)
            if (ierr .ne. nf_noerr) then
               write(stdout,2) 'u', density_name(1:lstr)
               goto 99                                         !--> ERROR
            endif
         else
            write(stdout,1) 'u', density_name(1:lstr)
            goto 99                                           !--> ERROR
         endif
    
! 3D momentum component in the ETA-direction.

         ierr=nf_inq_varid (ncid, 'v', varid)
         if (ierr .eq. nf_noerr) then
            ierr=nf_fread (v(START_2D_ARRAY,1,1), ncid, varid,
     &                                        0, v3dvar)
            if (ierr .ne. nf_noerr) then
               write(stdout,2) 'v', density_name(1:lstr)
               goto 99                                         !--> ERROR
            endif
         else
            write(stdout,1) 'v', density_name(1:lstr)
            goto 99                                           !--> ERROR
         endif     

!
!  Close input NetCDF file.
!
      ierr=nf_close(ncid)

! --- Compute barotropic u velocity ---
!
      ubar(:,:,1)=0.
      do k=1,N
        do j=Jstr,Jend
          do i=Istr,Iend
            ubar(i,j,1)=ubar(i,j,1)+u(i,j,k,1)
     &                             *0.5*(Hz(i-1,j,k)+Hz(i,j,k))
          enddo
        enddo
      enddo
      do j=Jstr,Jend
        do i=Istr,Iend
          ubar(i,j,1) = ubar(i,j,1) /
     &                  ( 0.5*(zeta(i,j,1)+zeta(i-1,j,1))
     &                  + 0.5*(h(i,j)+h(i-1,j)) )
        enddo
      enddo
!
! --- Compute barotropic v velocity ---
!
      vbar(:,:,1)=0.
      do k=1,N
        do j=Jstr,Jend
          do i=Istr,Iend
            vbar(i,j,1)=vbar(i,j,1)+v(i,j,k,1)
     &                           *0.5*(Hz(i,j-1,k)+Hz(i,j,k))
          enddo
        enddo
      enddo
      do j=Jstr,Jend
        do i=Istr,Iend
          vbar(i,j,1) = vbar(i,j,1) / 
     &                ( 0.5*(zeta(i,j,1)+zeta(i,j-1,1))
     &                + 0.5*(h(i,j)+h(i,j-1)) )
        enddo
      enddo

#ifdef CLIMATOLOGY
!
! Fill climatological arrays. MPI exchanges are done in analytical
!
      ubclm(:,:)=0.
      do k=1,N
        do j=Jstr,Jend
          do i=Istr,Iend
            ubclm(i,j)=ubclm(i,j)+uclm(i,j,k)
     &                             *0.5*(Hz(i-1,j,k)+Hz(i,j,k))
          enddo
        enddo
      enddo
      do j=Jstr,Jend
        do i=Istr,Iend
          ubclm(i,j) = ubclm(i,j) /
     &                  ( 0.5*(zeta(i,j,1)+zeta(i-1,j,1))
     &                  + 0.5*(h(i,j)+h(i-1,j)) )
        enddo
      enddo
!
! --- Compute barotropic u velocity ---
!
      vbclm(:,:)=0.
      do k=1,N
        do j=Jstr,Jend
          do i=Istr,Iend
            vbclm(i,j)=vbclm(i,j)+vclm(i,j,k)
     &                           *0.5*(Hz(i,j-1,k)+Hz(i,j,k))
          enddo
        enddo
      enddo
      do j=Jstr,Jend
        do i=Istr,Iend
          vbclm(i,j) = vbclm(i,j) / 
     &                ( 0.5*(zeta(i,j,1)+zeta(i,j-1,1))
     &                + 0.5*(h(i,j)+h(i,j-1)) )
        enddo
      enddo
#endif  /* CLIMATOLOGY */

#if defined EW_PERIODIC || defined NS_PERIODIC || defined  MPI
      call exchange_u3d_tile (Istr,Iend,Jstr,Jend,
     &                        u(START_2D_ARRAY,1,1))
      call exchange_v3d_tile (Istr,Iend,Jstr,Jend,
     &                        v(START_2D_ARRAY,1,1))
      do itrc=1,NT
        call exchange_r3d_tile (Istr,Iend,Jstr,Jend,
     &                          t(START_2D_ARRAY,1,1,itrc))
      enddo
#endif

#if defined EW_PERIODIC || defined NS_PERIODIC || defined  MPI
      call exchange_r2d_tile (Istr,Iend,Jstr,Jend,
     &                        zeta(START_2D_ARRAY,1))
      call exchange_u2d_tile (Istr,Iend,Jstr,Jend,
     &                        ubar(START_2D_ARRAY,1))
      call exchange_v2d_tile (Istr,Iend,Jstr,Jend,
     &                        vbar(START_2D_ARRAY,1))
#endif

#ifdef CLIMATOLOGY

#if defined EW_PERIODIC || defined NS_PERIODIC || defined  MPI
      call exchange_u3d_tile (Istr,Iend,Jstr,Jend,
     &                        uclm(START_2D_ARRAY,1))
      call exchange_v3d_tile (Istr,Iend,Jstr,Jend,
     &                        vclm(START_2D_ARRAY,1))
      do itrc=1,NT
        call exchange_r3d_tile (Istr,Iend,Jstr,Jend,
     &                          tclm(START_2D_ARRAY,1,itrc))
      enddo
#endif

#if defined EW_PERIODIC || defined NS_PERIODIC || defined  MPI
      call exchange_r2d_tile (Istr,Iend,Jstr,Jend,
     &                        ssh(START_2D_ARRAY))
      call exchange_u2d_tile (Istr,Iend,Jstr,Jend,
     &                        ubclm(START_2D_ARRAY))
      call exchange_v2d_tile (Istr,Iend,Jstr,Jend,
     &                        vbclm(START_2D_ARRAY))
#endif

#endif  /* CLIMATOLOGY */

      trd=omp_get_thread_num()
      call rho_eos_tile (Istr,Iend,Jstr,Jend, A2d(1,1,trd), 
     &                                        A2d(1,2,trd))
      
  1   format(/1x,'ANA_GET_DENSITY - unable to find variable:',    1x,A,
     &                            /15x,'in input NetCDF file:',1x,A/)
  2   format(/1x,'ANA_GET_DENSITY - error while reading variable:',1x, A,
     &    2x,'in input NetCDF file:',1x,A/)
  3   format(/1x,'ANA_GET_DENSITY - unable to find variable:',    1x,A,
     &                            /15x,'in input NetCDF file:',1x,A,
     &    1x,'-> analytical value'/)
     
 99   return
 
!aponte density)
      end
