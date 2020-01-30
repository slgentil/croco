! $Id: get_initial.F 697 2011-04-11 12:35:17Z gcambon $
!
!======================================================================
! ROMS_AGRIF is a branch of ROMS developped at IRD and INRIA, in France
! The two other branches from UCLA (Shchepetkin et al) 
! and Rutgers University (Arango et al) are under MIT/X style license.
! ROMS_AGRIF specific routines (nesting) are under CeCILL-C license.
! 
! ROMS_AGRIF website : http://roms.mpl.ird.fr
!======================================================================
!
! aponte itide (
#include "cppdefs.h"
#if defined ITIDE

      subroutine get_vmodes  ! get the vertical structure of normal modes

      implicit none
# include "param.h"
# include "scalars.h"
# include "ncscrum.h"
# include "ocean2d.h"
# include "ocean3d.h"
      integer ncid, indx, varid,  ierr, lstr, lvar, latt, lenstr,
     &        ibuff(4),   nf_fread, checkdims
      character units*64
#include "netcdf.inc"
#include "buffer.h"
      integer type,   vert_type,  imin, imax,  start(3)
     &      , record, horiz_type, jmin, jmax,  count(3)
     &      , i,j,k, shift
      integer nmode



      if (may_day_flag.ne.0) return      !-->  EXIT
!      lstr=lenstr(vmodesname)
!      ierr=nf_open(ininame(1:lstr), nf_nowrite, ncid)
      ierr=nf_open('vmodes_in.nc', nf_nowrite, ncid)
      if (ierr .eq. nf_noerr) then
      else
        write(stdout,'(/1x,2A/15x,3A)') 'GET_VMODES ERROR: Cannot ',
     &               'open netCDF file', '''' ,'''.'
        goto 99                                           !--> ERROR
      endif

!
! Read rho modes
! ---- -- -------- ----- ----------
!
!      lvar=lenstr(vname(1,indxTime))
!      ierr=nf_inq_varid (ncid, vname(1,indxTime)(1:lvar), varid)
      ierr=nf_inq_varid (ncid, 'phi_r', varid)
      if (ierr .eq. nf_noerr) then

! $Id: nf_fread.F 697 2011-04-11 12:35:17Z gcambon $
! Read/Write a floating point array from/into an input/output
! NetCDF file.

      start(1)=1
      start(2)=1
      start(3)=1
      count(1)=N
      count(2)=5
      count(3)=1

!
! Read/Write array from the disk.
!===== ===== ===== ==== === =====
!
      ierr=nf_get_vara_FTYPE (ncid, varid, start, count, buff(1))
      nf_fread=ierr
      if (ierr .ne. nf_noerr) then
        write(*,'(/1x,2A,I5,1x,A,I4/)') 'NF_FREAD ERROR: ',
     &             'nf_get_vara netCDF error code =', ierr
     &              MYID
        return
      endif

      do nmode=1,count(2)
       do k=1,count(1)
!        do j=jmin,jmax
!          shift=1-imin+count(1)*(j-jmin+(k-1)*count(2))
!          do i=imin,imax
            phir(k,nmode)=buff(k+(nmode-1)*count(1))
!          enddo
!        enddo
       enddo
      enddo

! debug
!       do k=1,count(1)
!          write(*,*) phir(k,2)
!       enddo

      else
        write(stdout,1) 'phir', 'vmodes_in.nc'
        goto 99                                           !--> ERROR
      endif


!
! Read uv modes
! ---- -- -------- ----- ----------
!
!      lvar=lenstr(vname(1,indxTime))
!      ierr=nf_inq_varid (ncid, vname(1,indxTime)(1:lvar), varid)
      ierr=nf_inq_varid (ncid, 'phi_uv', varid)
      if (ierr .eq. nf_noerr) then

! $Id: nf_fread.F 697 2011-04-11 12:35:17Z gcambon $
! Read/Write a floating point array from/into an input/output
! NetCDF file.

!
! Read/Write array from the disk.
!===== ===== ===== ==== === =====
!
      ierr=nf_get_vara_FTYPE (ncid, varid, start, count, buff(1))
      nf_fread=ierr
      if (ierr .ne. nf_noerr) then
        write(*,'(/1x,2A,I5,1x,A,I4/)') 'NF_FREAD ERROR: ',
     &             'nf_get_vara netCDF error code =', ierr
     &              MYID
        return
      endif

      do nmode=1,count(2)
       do k=1,count(1)
          phiuv(k,nmode)=buff(k+(nmode-1)*count(1))
       enddo
      enddo



!
! Exchange periodic and computational margins (reader only).
!
#ifdef MPI
#define LOCALLM Lmmpi
#define LOCALMM Mmmpi
#else
#define LOCALLM Lm
#define LOCALMM Mm
#endif
#if defined EW_PERIODIC || defined NS_PERIODIC  || defined MPI
!      if (horiz_type.eq.0 .and. vert_type.eq.0) then
!        call exchange_r2d_tile (1,LOCALLM,1,LOCALMM, A)
!      elseif (horiz_type.eq.1 .and. vert_type.eq.0) then
!        call exchange_u2d_tile (1,LOCALLM,1,LOCALMM, A)
!      elseif (horiz_type.eq.2 .and. vert_type.eq.0) then
!        call exchange_v2d_tile (1,LOCALLM,1,LOCALMM, A)
!      elseif (horiz_type.eq.3 .and. vert_type.eq.0) then
!        call exchange_p2d_tile (1,LOCALLM,1,LOCALMM, A)
# ifdef SOLVE3D
!      elseif (horiz_type.eq.0 .and. vert_type.eq.1) then
!        call exchange_r3d_tile (1,LOCALLM,1,LOCALMM, A)
!      elseif (horiz_type.eq.1 .and. vert_type.eq.1) then
!        call exchange_u3d_tile (1,LOCALLM,1,LOCALMM, A)
!      elseif (horiz_type.eq.2 .and. vert_type.eq.1) then
!        call exchange_v3d_tile (1,LOCALLM,1,LOCALMM, A)
!      elseif (horiz_type.eq.3 .and. vert_type.eq.1) then
!        call exchange_p3d_tile (1,LOCALLM,1,LOCALMM, A)
!      elseif (horiz_type.eq.0 .and. vert_type.eq.2) then
!        call exchange_w3d_tile (1,LOCALLM,1,LOCALMM, A)
# endif
!      endif
#endif


      else
        write(stdout,1) 'phiuv', 'vmodes_in.nc'
        goto 99                                           !--> ERROR
      endif


!
!  Close input NetCDF file.
!
      ierr=nf_close(ncid)
  1   format(/1x,'GET_INITIAL - unable to find variable:',    1x,A,
     &                            /15x,'in input NetCDF file:',1x,A/)

      return
  99  may_day_flag=2
      return
      end

#else
      subroutine get_vmodes_empty
      return
      end
#endif /* ITIDE */


! aponte itide )
