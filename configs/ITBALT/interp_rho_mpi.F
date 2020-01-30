! $Id: interp_rho.F 697 2011-04-11 12:35:17Z gcambon $
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
      subroutine interp_r3d_type_mpi (A, ifld, nfltmax, indx)
      implicit none
# include "param.h"
# include "grid.h"
# include "floats.h"
# include "scalars.h"
      integer ifld, nfltmax, indx(nfltmax)
      real A(GLOBAL_2D_ARRAY,N)
      integer id,iflt, i1,i2, j1,j2, k1,k2
      real p1,p2, q1,q2, r1,r2, cff1, cff2
      integer node_float

      integer node_offx
      integer node_offy
      integer i_offset, j_offset
      i_offset=mod(mynode,NP_XI)*Lm
      j_offset=(mynode/NP_XI)*Mm

      do id=1,nfltmax
        iflt=indx(id)

! check whether float is within mpi tile
       node_float= mod(int((track(ixgrd,nfp1,iflt)-0.5)/dble(Lm)),NP_XI)
       node_float= node_float
     &               + int((track(iygrd,nfp1,iflt)-0.5)/dble(Mm))*NP_XI

       if ( node_float .eq. mynode  ) then

        k1=int(track(izgrd,nfp1,iflt)+0.5)
        r2=track(izgrd,nfp1,iflt)+0.5 - float(k1)
        r1=1.-r2
        k1=max(k1,   1)
        k2=min(k1+1, N)

        i1=int(track(ixgrd,nfp1,iflt))
        i1=i1-i_offset
        i2=i1+1
!        p2=track(ixgrd,nfp1,iflt)-float(i1)
        p2=track(ixgrd,nfp1,iflt)-float(i1+i_offset)
        p1=1.-p2

        j1=int(track(iygrd,nfp1,iflt))
        j1=j1-j_offset
        j2=j1+1
!        q2=track(iygrd,nfp1,iflt)-float(j1)
        q2=track(iygrd,nfp1,iflt)-float(j1+j_offset)
        q1=1.0-q2

#ifdef MASKING
        cff1=p1*q1*rmask(i1,j1)*(r1*A(i1,j1,k1)+r2*A(i1,j1,k2))
     &      +p2*q1*rmask(i2,j1)*(r1*A(i2,j1,k1)+r2*A(i2,j1,k2))
     &      +p1*q2*rmask(i1,j2)*(r1*A(i1,j2,k1)+r2*A(i1,j2,k2))
     &      +p2*q2*rmask(i2,j2)*(r1*A(i2,j2,k1)+r2*A(i2,j2,k2))

        cff2=q1*(p1*rmask(i1,j1) + p2*rmask(i2,j1))
     &      +q2*(p1*rmask(i1,j2) + p2*rmask(i2,j2))

        if (cff2.gt.0.) then
          trackaux(ifld,iflt)=cff1/cff2
        else
          trackaux(ifld,iflt)=0.0
        endif
#else
        trackaux(ifld,iflt)=p1*q1*(r1*A(i1,j1,k1)+r2*A(i1,j1,k2))
     &                       +p2*q1*(r1*A(i2,j1,k1)+r2*A(i2,j1,k2))
     &                       +p1*q2*(r1*A(i1,j2,k1)+r2*A(i1,j2,k2))
     &                       +p2*q2*(r1*A(i2,j2,k1)+r2*A(i2,j2,k2))
#endif

       endif

      enddo
      return
      end


      subroutine interp_u3d_type_mpi (A, ifld, nfltmax, indx)
      implicit none
# include "param.h"
# include "grid.h"
# include "floats.h"
# include "scalars.h"
      integer ifld, nfltmax, indx(nfltmax)
      real A(GLOBAL_2D_ARRAY,N)
      integer id,iflt, i1,i2, j1,j2, k1,k2
      real p1,p2, q1,q2, r1,r2, cff1, cff2
      integer node_float

      integer node_offx
      integer node_offy
      integer i_offset, j_offset
      i_offset=mod(mynode,NP_XI)*Lm
      j_offset=(mynode/NP_XI)*Mm

      do id=1,nfltmax
        iflt=indx(id)

! check whether float is within mpi tile
       node_float= mod(int((track(ixgrd,nfp1,iflt)-0.5)/dble(Lm)),NP_XI)
       node_float= node_float
     &               + int((track(iygrd,nfp1,iflt)-0.5)/dble(Mm))*NP_XI

       if ( node_float .eq. mynode  ) then

        k1=int(track(izgrd,nfp1,iflt)+0.5)
        r2=track(izgrd,nfp1,iflt)+0.5 - float(k1)
        r1=1.-r2
        k1=max(k1,   1)
        k2=min(k1+1, N)

        i1=int(track(ixgrd,nfp1,iflt)-0.5)
        i1=i1-i_offset
        i2=i1+1
!        p2=track(ixgrd,nfp1,iflt)-float(i1)
        p2=track(ixgrd,nfp1,iflt)-float(i1+i_offset)
        p1=1.-p2

        j1=int(track(iygrd,nfp1,iflt))
        j1=j1-j_offset
        j2=j1+1
!        q2=track(iygrd,nfp1,iflt)-float(j1)
        q2=track(iygrd,nfp1,iflt)-float(j1+j_offset)
        q1=1.0-q2

#ifdef MASKING
        cff1=p1*q1*rmask(i1,j1)*(r1*A(i1,j1,k1)+r2*A(i1,j1,k2))
     &      +p2*q1*rmask(i2,j1)*(r1*A(i2,j1,k1)+r2*A(i2,j1,k2))
     &      +p1*q2*rmask(i1,j2)*(r1*A(i1,j2,k1)+r2*A(i1,j2,k2))
     &      +p2*q2*rmask(i2,j2)*(r1*A(i2,j2,k1)+r2*A(i2,j2,k2))

        cff2=q1*(p1*rmask(i1,j1) + p2*rmask(i2,j1))
     &      +q2*(p1*rmask(i1,j2) + p2*rmask(i2,j2))

        if (cff2.gt.0.) then
          trackaux(ifld,iflt)=cff1/cff2
        else
          trackaux(ifld,iflt)=0.0
        endif
#else
        trackaux(ifld,iflt)=p1*q1*(r1*A(i1,j1,k1)+r2*A(i1,j1,k2))
     &                       +p2*q1*(r1*A(i2,j1,k1)+r2*A(i2,j1,k2))
     &                       +p1*q2*(r1*A(i1,j2,k1)+r2*A(i1,j2,k2))
     &                       +p2*q2*(r1*A(i2,j2,k1)+r2*A(i2,j2,k2))
#endif

       endif

      enddo
      return
      end


      subroutine interp_v3d_type_mpi (A, ifld, nfltmax, indx)
      implicit none
# include "param.h"
# include "grid.h"
# include "floats.h"
# include "scalars.h"
      integer ifld, nfltmax, indx(nfltmax)
      real A(GLOBAL_2D_ARRAY,N)
      integer id,iflt, i1,i2, j1,j2, k1,k2
      real p1,p2, q1,q2, r1,r2, cff1, cff2
      integer node_float

      integer node_offx
      integer node_offy
      integer i_offset, j_offset
      i_offset=mod(mynode,NP_XI)*Lm
      j_offset=(mynode/NP_XI)*Mm

      do id=1,nfltmax
        iflt=indx(id)

! check whether float is within mpi tile
       node_float= mod(int((track(ixgrd,nfp1,iflt)-0.5)/dble(Lm)),NP_XI)
       node_float= node_float
     &               + int((track(iygrd,nfp1,iflt)-0.5)/dble(Mm))*NP_XI

       if ( node_float .eq. mynode  ) then

        k1=int(track(izgrd,nfp1,iflt)+0.5)
        r2=track(izgrd,nfp1,iflt)+0.5 - float(k1)
        r1=1.-r2
        k1=max(k1,   1)
        k2=min(k1+1, N)

        i1=int(track(ixgrd,nfp1,iflt))
        i1=i1-i_offset
        i2=i1+1
!        p2=track(ixgrd,nfp1,iflt)-float(i1)
        p2=track(ixgrd,nfp1,iflt)-float(i1+i_offset)
        p1=1.-p2

        j1=int(track(iygrd,nfp1,iflt)-0.5)
        j1=j1-j_offset
        j2=j1+1
!        q2=track(iygrd,nfp1,iflt)-float(j1)
        q2=track(iygrd,nfp1,iflt)-float(j1+j_offset)
        q1=1.0-q2

#ifdef MASKING
        cff1=p1*q1*rmask(i1,j1)*(r1*A(i1,j1,k1)+r2*A(i1,j1,k2))
     &      +p2*q1*rmask(i2,j1)*(r1*A(i2,j1,k1)+r2*A(i2,j1,k2))
     &      +p1*q2*rmask(i1,j2)*(r1*A(i1,j2,k1)+r2*A(i1,j2,k2))
     &      +p2*q2*rmask(i2,j2)*(r1*A(i2,j2,k1)+r2*A(i2,j2,k2))

        cff2=q1*(p1*rmask(i1,j1) + p2*rmask(i2,j1))
     &      +q2*(p1*rmask(i1,j2) + p2*rmask(i2,j2))

        if (cff2.gt.0.) then
          trackaux(ifld,iflt)=cff1/cff2
        else
          trackaux(ifld,iflt)=0.0
        endif
#else
        trackaux(ifld,iflt)=p1*q1*(r1*A(i1,j1,k1)+r2*A(i1,j1,k2))
     &                       +p2*q1*(r1*A(i2,j1,k1)+r2*A(i2,j1,k2))
     &                       +p1*q2*(r1*A(i1,j2,k1)+r2*A(i1,j2,k2))
     &                       +p2*q2*(r1*A(i2,j2,k1)+r2*A(i2,j2,k2))
#endif

       endif

      enddo
      return
      end


      subroutine interp_w3d_type_mpi (A, ifld, nfltmax,indx)
      implicit none
# include "param.h"
# include "floats.h"
# include "scalars.h"
      integer ifld, nfltmax, indx(nfltmax)
      real A(GLOBAL_2D_ARRAY,0:N)
      integer id,iflt, i1,i2, j1,j2, k1,k2
      real p1,p2, q1,q2, r1,r2
      integer node_float
      integer i_offset, j_offset
      i_offset=mod(mynode,NP_XI)*Lm
      j_offset=(mynode/NP_XI)*Mm

      do id=1,nfltmax
        iflt=indx(id)

! check whether float is within mpi tile
       node_float= mod(int((track(ixgrd,nfp1,iflt)-0.5)/dble(Lm)),NP_XI)
       node_float= node_float
     &               + int((track(iygrd,nfp1,iflt)-0.5)/dble(Mm))*NP_XI

       if ( node_float .eq. mynode  ) then

        k1=int(track(izgrd,nfp1,iflt))
        r2=track(izgrd,nfp1,iflt) - float(k1)
        r1=1.-r2
        k1=max(k1,   0)
        k2=min(k1+1, N)

        i1=int(track(ixgrd,nfp1,iflt))
        i1=i1-i_offset
        i2=i1+1
!        p2=track(ixgrd,nfp1,iflt)-float(i1)
        p2=track(ixgrd,nfp1,iflt)-float(i1+i_offset)
        p1=1.-p2

        j1=int(track(iygrd,nfp1,iflt))
        j1=j1-j_offset
        j2=j1+1
!        q2=track(iygrd,nfp1,iflt)-float(j1)
        q2=track(iygrd,nfp1,iflt)-float(j1+j_offset)
        q1=1.0-q2

        trackaux(ifld,iflt)=p1*q1*(r1*A(i1,j1,k1)+r2*A(i1,j1,k2))
     &                       +p2*q1*(r1*A(i2,j1,k1)+r2*A(i2,j1,k2))
     &                       +p1*q2*(r1*A(i1,j2,k1)+r2*A(i1,j2,k2))
     &                       +p2*q2*(r1*A(i2,j2,k1)+r2*A(i2,j2,k2))

       endif

      enddo
      return
      end


      subroutine interp_r2d_type_mpi (A, ifld, nfltmax, indx)
      implicit none
# include "param.h"
# include "floats.h"
# include "scalars.h"
      integer ifld, nfltmax, indx(nfltmax)
      real A(GLOBAL_2D_ARRAY)
      integer id,iflt, i1,i2, j1,j2
      real p1,p2, q1,q2
      integer node_float
      integer i_offset, j_offset
      i_offset=mod(mynode,NP_XI)*Lm
      j_offset=(mynode/NP_XI)*Mm

      do id=1,nfltmax
        iflt=indx(id)

! check whether float is within mpi tile
       node_float= mod(int((track(ixgrd,nfp1,iflt)-0.5)/dble(Lm)),NP_XI)
       node_float= node_float
     &               + int((track(iygrd,nfp1,iflt)-0.5)/dble(Mm))*NP_XI

       if ( node_float .eq. mynode  ) then

        i1=int(track(ixgrd,nfp1,iflt))
        i1=i1-i_offset
        i2=i1+1
!        p2=track(ixgrd,nfp1,iflt)-float(i1)
        p2=track(ixgrd,nfp1,iflt)-float(i1+i_offset)
        p1=1.-p2

        j1=int(track(iygrd,nfp1,iflt))
        j1=j1-j_offset
        j2=j1+1
!        q2=track(iygrd,nfp1,iflt)-float(j1)
        q2=track(iygrd,nfp1,iflt)-float(j1+j_offset)
        q1=1.0-q2

        trackaux(ifld,iflt)=q1*(p1*A(i1,j1) + p2*A(i2,j1))
     &                       +q2*(p1*A(i1,j2) + p2*A(i2,j2))

       endif

      enddo
      return
      end


      subroutine interp_r2d_type_ini_mpi (A, ifld)
      implicit none
# include "param.h"
# include "floats.h"
# include "scalars.h"
      integer ifld
      real A(GLOBAL_2D_ARRAY)
      integer iflt, i1,i2, j1,j2
      real p1,p2, q1,q2
      integer node_float
      integer i_offset, j_offset
      i_offset=mod(mynode,NP_XI)*Lm
      j_offset=(mynode/NP_XI)*Mm

      do iflt=1,nfloats       
!       node_float= mod(int((track(ixgrd,nfp1,iflt)-0.5)/dble(Lm)),NP_XI)
!       node_float= node_float
!    &               + int((track(iygrd,nfp1,iflt)-0.5)/dble(Mm))*NP_XI

!       if ( node_float .eq. mynode  ) then
        if (Tinfo(inode,iflt) == mynode) then
           i1=int(Tinfo(ixgrd,iflt))
           i1=i1-i_offset
           i2=i1+1
!          p2=Tinfo(ixgrd,iflt)-float(i1)
           p2=Tinfo(ixgrd,iflt)-float(i1+i_offset)
           p1=1.-p2

           j1=int(Tinfo(iygrd,iflt))
           j1=j1-j_offset
           j2=j1+1
!          q2=Tinfo(iygrd,iflt)-float(j1)
           q2=Tinfo(iygrd,iflt)-float(j1+j_offset)
           q1=1.0-q2

           trackaux(ifld,iflt)=q1*(p1*A(i1,j1) + p2*A(i2,j1))
     &                       +q2*(p1*A(i1,j2) + p2*A(i2,j2))
        endif
      enddo
      return
      end

#else
      subroutine interp_floats_mpi_empty
      end
#endif
