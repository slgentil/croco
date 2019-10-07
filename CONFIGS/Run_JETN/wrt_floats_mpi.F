! $Id: wrt_floats.F 1143 2013-05-17 08:17:23Z serena $
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
#include "cppdefs.h"
#if defined FLOATS && defined MPI
                               ! Writes requested model
      subroutine wrt_floats_mpi    ! fields at requested levels 
                               ! into history netCDF file.
# ifdef AGRIF
      USE Agrif_Util 
# endif
      implicit none
# include "param.h"
# include "scalars.h"
# include "ncscrum.h"
# include "ncscrum_floats.h"
# include "grid.h"
# include "ocean2d.h"
# include "ocean3d.h"
# include "mixing.h"
# include "mpi_roms.h"
# include "floats.h"
# include "netcdf.inc"

      integer ierr,  record, lvar, lenstr, iflt, id
     &   , nfltrelmax, nfltunrelmax, indxrel(Mfloats)
     &   , indxunrel(Mfloats),Toutint(Mfloats)
     &   , start(2), count(2), ibuff(2), nf_fwrite
      logical flt_innode
      integer node_float
      logical newf, file_exists
      character*65 vinfo
      real Tout(Mfloats)
! ascii file
      character*20 filename
      character*20 rowfmt


! ----- --- -------- ----- ----------
!
! Save indices of released and non released floats 
!

! !!! Should only output for tiles that contain floats !!!

      if ( ANY( trackinfo(2,1:nfloats)==mynode ) ) then
      	where(trackinfo(2,1:nfloats)== mynode) 
     &        Tout(1:nfloats)=trackinfo(1,1:nfloats)
      	if (any(Tout(1:nfloats) == 1)) then

     ! ascii output file
      write (filename, "(A6,I0.4)") "float.", mynode
      
! If file does not exist, create it else append to it 
      inquire(file=filename, exist=file_exists)
      if ( file_exists) then 
! posision="append" is not 77 and may not work with some compilers
        open(unit=12, file= trim(filename), action="write", status="old",
     &              position="append")
      else
        open(unit=12, file= trim(filename), action="write")
      endif

      do iflt=1,nfloats

         ! check whether float is within mpi tile
         if ( trackinfo(1,iflt) .eq. 1 
     &      .and. trackinfo(2,iflt) .eq. mynode ) then

            write(12,2)
     &         iflt, tdays 
     &         , trackaux(ixgrd,iflt), trackaux(iygrd,iflt)
     &         , trackaux(izgrd,iflt), trackaux(ifdpt,iflt)
     &         , trackaux(iftem,iflt), trackaux(ifu,iflt) 
     &         , trackaux(ifv,iflt), trackaux(ifdudt,iflt) 
     &         , trackaux(ifdvdt,iflt), trackaux(ifpres,iflt) 

         endif

      enddo

      ! close output file
      close(12)

      ! ascii output format
   2  format(I4,11(1x,e16.8)) 

      	endif
      endif
      
      return
      end


#else
      subroutine wrt_floats_mpi_empty
      return
      end
#endif /* FLOATS */

