! $Id: init_floats.F 1458 2014-02-03 15:01:25Z gcambon $
!
!======================================================================
! ROMS_AGRIF is a branch of ROMS developped at IRD and INRIA, in France
! The two other branches from UCLA (Shchepetkin et al) 
! and Rutgers University (Arango et al) are under MIT/X style license.
! ROMS_AGRIF specific routines (nesting) are under CeCILL-C license.
! 
! ROMS_AGRIF website : http://www.romsagrif.org
!======================================================================
!
#include "cppdefs.h"

! aponte flt (
       subroutine init_floats_mpi
#if defined FLOATS && defined MPI
!
!================================================== John M. Klinck ===
!  Copyright (c) 2000 Rutgers/UCLA                                   !
!================================================ Hernan G. Arango ===
!                                                                    !
!  This routine reads in and process initial float locations from    !
!  input floats script.                                              !
!                                                                    !
!=====================================================================
!
      implicit none
# include "param.h"
# include "grid.h"
# include "ocean3d.h"
# include "ncscrum_floats.h"
# include "scalars.h"
# include "floats.h"
# include "init_floats.h"
!
      integer Ncount, i, icard, iunit, j, k, lstr, mc, nc
      integer index(Mfloats)
      integer lenstr
      real xfloat, yfloat, zfloat
      real Ip0(Mfloats), Jp0(Mfloats), lat(Mfloats), lon(Mfloats)
      character*35 frmt
      parameter (iunit=50)
! aponte flt (
#ifdef MPI
      integer node_float ! mpi node number
#endif

      real init
      parameter (init=0.)

! from init_arrays_floats.F

!      integer iflt, j, i
      integer iflt
      
!  Initialize floats trajectory array.
!
         do iflt=1,Mfloats
            fltgrd(iflt)=-1 ! -1 is the dead float value. 
            nrecvel(iflt)=init

            Tinfo(inode,iflt)=-1
            do j=igrd,izgrd  
               Tinfo(j,iflt)=init ! Float trajectory initial information.
            enddo
        
            do j=1,NDIAGS
               trackaux(j,iflt)=init ! Multivariate float data for writing 
            enddo

            do j=0,NFT  ! Number of float time levels (minus one). 
              do i=1,NFV  ! Number of float variables.
                track(i,j,iflt)=init
              enddo
            enddo
             
            trackinfo(1,iflt) = -1
            trackinfo(2,iflt) = -1
            trackinfo(3,iflt) = 0

         enddo

# ifdef AGRIF
         do j=1,maxgrids   ! consistent with initialization of nfp1
           floattindex(j)=0
         enddo
# endif

      if ( nrpfflt.gt.0 ) then
!        restart from float.rst.mynode       
         call get_initial_floats_mpi
      else

!        start from floats.in
!

! aponte flt )

!
!---------------------------------------------------------------------
!  Read in initial float locations.
!---------------------------------------------------------------------
!
         lstr=lenstr(fposnam)
         open(iunit,file=fposnam(1:lstr),form='formatted',
     &     status='old', err=195)
c     write(stdout,10)
c 10  format(/,' FLOATS processing parameters:',/)
!
!  Read input parameters according to their input card number.
!
         icard=0
         do while (icard.lt.99)
!
!  Read in floats identification title.
!
           if (icard.eq.1) then
             read(iunit,'(a)',err=70) Ftitle
             lstr=lenstr(Ftitle)
             if (mynode==0) write(stdout,20) Ftitle(1:lstr)
  20         format(1x,'(',a,')',/)
!
!  Read in initial floats location.
!
           elseif (icard.eq.2) then
             i=0
             nfloats=0
c         write(stdout,30)
c30       format(1x,'Fcoor',2x,'Ftype',2x,'Fcount',2x,'Ft0',
c    &        2x,'Fx0',2x,'Fy0',2x,'Fz0',2x,'Fgrd',2x,'Fdt',
c    &                           2x,'Fdx',2x,'Fdy',2x,'Fdz',/)
             do while (.true.)
               i=i+1
               read(iunit,*,err=40) Ft0(i), Fx0(i), Fy0(i), Fz0(i),
     &                    Fgrd(i),Fcoor(i), Ftype(i), Fcount(i),
     &                           Fdt(i), Fdx(i), Fdy(i), Fdz(i)
               nfloats=nfloats+Fcount(i)
c# if !defined SPHERICAL
c           if (Fcoor(i).eq.1) then
c             frmt='(i1,i2,i5,f10.4,2f8.0,f8.2,2x,i1,4f9.3)'
c           else
c            frmt='(i1,i2,i5,f10.4,3f8.2,2x,i1,4f9.3)'
c           endif
c# else 
c             frmt='(i1,i2,i5,f10.4,3f8.2,2x,i1,4f9.3)'
c# endif
c           write(stdout,frmt) Fcoor(i), Ftype(i), Fcount(i),
c    &                         Ft0(i), Fx0(i), Fy0(i), Fz0(i),
c    &                   Fgrd(i),Fdt(i), Fdx(i), Fdy(i), Fdz(i)
             enddo
  40         if (Ft0(i).ne.99.0) then
               write(stdout,50) icard, i, fposnam
  50           format(/,' INIT_FLOATS - error while reading input card: ',
     &             i2, ', floats location entry: ',i3,/,15x,
     &             'in input script: ',a)
               stop
             else
               Ncount=i-1
               i_floats=i-1
               goto 90
             endif
             if (Mfloats.lt.nfloats) then
               write(stdout,60) Mfloats, nfloats
  60           format(/,' INIT_FLOATS - too small dimension parameter,',
     &             ' Mfloats',2i6,/,15x,'change file  param.h  and',
     &             ' recompile.')
               stop
             endif
           endif
!
!  Read last input card ID.
!
           read(iunit,*,err=70) icard
         enddo
         goto 90
!
!  Error while reading input parameters.
!
  70     write(stdout,80) icard, fposnam
  80     format(/,' INIT_FLOATS - error while reading input card: ',
     &       i2,15x,'in input script: ',a)

  90     close(iunit)

         if ( mynode .eq. 0 ) then 
           write(stdout,100) nfloats
 100       format(/,2x,i6,4x,'nfloats',t26,
     &       'Number of float trajectories to compute.',/)
         endif

!
!---------------------------------------------------------------------
!  Set initial float location.
!---------------------------------------------------------------------
!
!  Set time of float release (seconds after model initialization) and
!  initial float horizontal positions (grid units).  Fill the initial
!  vertical level or depth position.
!
         mc=0
         nc=0
         do i=1,Ncount
           if (Fcount(i).eq.1) then
             nc=nc+1
             if (Ft0(i)<=0) then
                Tinfo(itstr,nc)=(tdays-Ft0(i))*day2sec ! <=0 relative time, >0 absolute time
             else
                Tinfo(itstr,nc)= Ft0(i) *day2sec
             endif
             Tinfo(izgrd,nc)=Fz0(i)
             Tinfo(igrd,i)=FLOAT(Fgrd(i)) ! each float belongs to this grid initially
             if (Fcoor(i).eq.0) then
                Tinfo(ixgrd,nc)=MIN(MAX(0.5,Fx0(i)),FLOAT(LLm)+0.5)
                track(ixgrd,nfp1,nc)=Tinfo(ixgrd,nc)    
                Tinfo(iygrd,nc)=MIN(MAX(0.5,Fy0(i)),FLOAT(MMm)+0.5)
                track(iygrd,nfp1,nc)=Tinfo(iygrd,nc)
             else
               mc=mc+1
               lon(mc)=Fx0(i)
               lat(mc)=Fy0(i)
               index(mc)=nc
             endif
           elseif (Fcount(i).gt.1) then   ! not verified
             do j=1,Fcount(i)
               nc=nc+1
               Tinfo(igrd,nc)=FLOAT(Fgrd(i)) ! each float belongs to this grid initially
! aponte flt (
               trackinfo(1,nc)=0
               trackinfo(3,nc)=Tinfo(igrd,nc)
! aponte flt )
               if (Fdt(i).gt.0.0) then
                 if (Ft0(i)<=0) then
                    Tinfo(itstr,nc)=(tdays-Ft0(i)+FLOAT(j-1)*Fdt(i))*day2sec ! <=0 relative time, >0 absolute time
                 else
                    Tinfo(itstr,nc)= (Ft0(i)+FLOAT(j-1)*Fdt(i)) *day2sec
                 endif
                 Tinfo(izgrd,nc)=Fz0(i)
                 if (Fcoor(i).eq.0) then
                   Tinfo(ixgrd,nc)=MIN(MAX(0.5,Fx0(i)),FLOAT(LLm)+0.5)
                   Tinfo(iygrd,nc)=MIN(MAX(0.5,Fy0(i)),FLOAT(MMm)+0.5)
                 else
                   mc=mc+1
                   lon(mc)=Fx0(i)
                   lat(mc)=Fy0(i)
                   index(mc)=nc
                 endif
               else
                 if (Ft0(i)<=0) then
                    Tinfo(itstr,nc)=(tdays-Ft0(i))*day2sec ! <=0 relative time, >0 absolute time
                 else
                    Tinfo(itstr,nc)= Ft0(i) *day2sec
                 endif
                 if (Fdz(i).eq.0.0) then
                   Tinfo(izgrd,nc)=Fz0(i)
                 else
                   if (Fz0(i).gt.0.0) then
                     zfloat=Fz0(i)+FLOAT(j-1)*Fdz(i)
                     Tinfo(izgrd,nc)=MIN(MAX(0.0,zfloat),FLOAT(N))
                   else
                     Tinfo(izgrd,nc)=Fz0(i)+FLOAT(j-1)*Fdz(i)
                   endif
                 endif
                 if (Fcoor(i).eq.0) then
                   xfloat=Fx0(i)+FLOAT(j-1)*Fdx(i)
                   Tinfo(ixgrd,nc)=MIN(MAX(0.5,xfloat),
     &                                FLOAT(LLm)+0.5)
                   yfloat=Fy0(i)+FLOAT(j-1)*Fdy(i)
                   Tinfo(iygrd,nc)=MIN(MAX(0.5,yfloat),
     &                                FLOAT(MMm)+0.5)
                 else
                   mc=mc+1
                   index(mc)=nc
                   lon(mc)=Fx0(i)+FLOAT(j-1)*Fdx(i)
                   lat(mc)=Fy0(i)+FLOAT(j-1)*Fdy(i)
                 endif
               endif
             enddo
           endif
         enddo

!
!  Set number of floats trajectories to compute.
!
         nfloats=nc
!
!  If applicable, convert floats initial (lon,lat) positions to grid
!  units.
!
! aponte flt this does not work at the moment
         if (mc.gt.0) then
# ifdef SPHERICAL
           call hindices (Ip0,Jp0,lon,lat,mc,lonr,latr)
# else
           call hindices (Ip0,Jp0,lon,lat,mc,xr,yr)
# endif
           do i=1,mc
             nc=index(i)
             Tinfo(ixgrd,nc)=MIN(MAX(0.5,Ip0(i)),FLOAT(LLm)+0.5)
             Tinfo(iygrd,nc)=MIN(MAX(0.5,Jp0(i)),FLOAT(MMm)+0.5)
           enddo
         endif
!
! determine node number

         do nc=1,Nfloats
           node_float= mod(int((Tinfo(ixgrd,nc)-0.5)/dble(Lm)),NP_XI)
           node_float= node_float 
     &                 + int((Tinfo(iygrd,nc)-0.5)/dble(Mm))*NP_XI
           Tinfo(inode,nc)=dble(node_float)
           trackinfo(2,nc)=node_float
         enddo
      
# ifdef SPHERICAL
         call interp_r2d_type_ini_mpi (lonr(START_2D_ARRAY), iflon)
         call interp_r2d_type_ini_mpi (latr(START_2D_ARRAY), iflat)
# else
         call interp_r2d_type_ini_mpi (  xr(START_2D_ARRAY), iflon)
         call interp_r2d_type_ini_mpi (  yr(START_2D_ARRAY), iflat)
# endif
# ifdef SOLVE3D
         call fill_ini ! fills in trackaux for ixgrd,iygrd,izgrd
                    ! and ifld (either izgrd or ifld is meaningful)
# endif

      endif /* nrfpflt */

      ! write initial position
      if (ldefflt) call wrt_floats_mpi
          
      return

 195  write(stdout,205) fposnam
 205  format(/,'float file ',A,/,' not found => float initialization', 
     &   ' using restart file')


C! Test de l'initialisation des flotteurs
C      do nc=1,nfloats
C         if (Tinfo(inode,nc)==mynode) then
C            write(stdout,53) mynode, nc, Tinfo(ixgrd,nc), Tinfo(iygrd,nc),
C     &      Tinfo(izgrd,nc), NINT(Tinfo(igrd,nc)), Tinfo(itstr,nc)*sec2day
C            format(/,'node ',i4,' float ',I4,5(1x,F6.1))
C         enddo
C      enddo
C! Fin test

#endif /* FLOATS && MPI */
! aponte flt )
      return
      end
