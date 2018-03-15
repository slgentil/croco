#ifndef MP_3PTS
      subroutine exchange_3d_tile_init_nbq (iv,Istr,Iend,Jstr,Jend, A)
#else
      subroutine exchange_3d_3pts_tile_init_nbq (iv,Istr,Iend,Jstr,Jend, A)
#endif
!
! Set periodic boundary conditions (if any) for a three-dimensional
! field A of RHO-, U-, V- or PSI-type. This file is designed to
! generate five different subroutines, by redefining (via CPP) the
! name of the subroutine exchange_2d_tile above and the starting
! indices ISTART = [Istr for U-,PSI-type; IstrR for V-,RHO-type]
! and JSTART = [Jstr for V-,PSI-type; JstrR for U-,RHO-type] below,
! as well as macro KSTART for the vertical RHO- and W-types. See
! also mounting file exchange.F
!
      implicit none
#include "param.h"
#include "scalars.h"
      integer Npts,ipts,jpts
# ifndef MP_3PTS
      parameter (Npts=2)
# else
      parameter (Npts=3)
# endif
      real A(GLOBAL_2D_ARRAY,KSTART:N)
      integer iv,Istr,Iend,Jstr,Jend, i,j,k
!
#include "compute_auxiliary_bounds.h"
!
              /* NS_PERIODIC */

#ifdef MPI
      k=N-KSTART+1
# ifndef MP_3PTS
      call MessPass3D_tile_init_nbq (iv,Istr,Iend,Jstr,Jend,  A,KSTART,N)
# else
      call MessPass3D_3pts_tile_init_nbq (iv,Istr,Iend,Jstr,Jend,  A,KSTART,N)
# endif
#endif
      return
      end

# ifndef MP_3PTS
#  define MP_3PTS
#  include "exchange_3d_tile_init_nbq.h"
#  undef MP_3PTS
# endif

