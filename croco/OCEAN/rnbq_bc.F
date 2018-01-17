! $Id: rnbq_bc.F 1458 2014-02-03 15:01:25Z gcambon $
!
!======================================================================
! CROCO is a branch of ROMS developped at IRD and INRIA, in France
! The two other branches from UCLA (Shchepetkin et al) 
! and Rutgers University (Arango et al) are under MIT/X style license.
! CROCO specific routines (nesting) are under CeCILL-C license.
! 
! CROCO website : http://www.croco-ocean.org
!======================================================================
!
#ifndef CHILD
!
# include "cppdefs.h"
# ifdef NBQ
      subroutine rnbq_bc_tile (Istr,Iend,Jstr,Jend,grad)
#  ifdef AGRIF      
      use AGRIF_Util
      integer Istr,Iend,Jstr,Jend
      real grad(PRIVATE_2D_SCRATCH_ARRAY)
      if (AGRIF_Root()) then
        call rnbq_bc_parent_tile (Istr,Iend,Jstr,Jend,grad)
      else
        call rnbq_bc_child_tile (Istr,Iend,Jstr,Jend,grad)
      endif
      return
      end
!
! PARENT
!
      subroutine rnbq_bc_parent_tile (Istr,Iend,Jstr,Jend,grad)
#  endif
!
! Set lateral boundary conditions for W-component momentum
! rhp_nbq_a(W_index) for the parent grid.
!
# endif /* NBQ */
#else
# ifdef NBQ
!
! CHILD
!
      subroutine rnbq_bc_child_tile (Istr,Iend,Jstr,Jend,grad)
!
! Set lateral boundary conditions for W-component momentum
! rhp_nbq_a(W_index) for the child grid.
!
# endif /* NBQ */
#endif /* CHILD */
#ifdef NBQ
!
! Common Code
!
# include "set_obc_definitions.h"
!
      USE module_nh
      USE module_nbq
!
      implicit none
# include "param.h"
# include "grid.h"
# include "ocean3d.h"
# include "climat.h"
# include "scalars.h"
# include "boundary.h"  
# ifdef MRL_WCI
#  include "forces.h"
# endif
      integer Istr,Iend,Jstr,Jend, i,j,k,
     &        l_b_j,l_b_jm1,l_bm1_j,l_bm1_jm1,l_bm2_j,
     &        l_b_i,l_b_im1,l_bm1_i,l_bm1_im1,l_bm2_i,
     &        l_c,l_c_im1,l_c_jm1
      real    grad(PRIVATE_2D_SCRATCH_ARRAY)
      real    cff,eps,
     &        cx,cy, dft,dfx,dfy, tau,tau_in,tau_out
      parameter (eps=1.E-20)
!
# include "compute_auxiliary_bounds.h"
!
! Interpolations of the parent values to get rnbqbry_east or rnbqclm
!
# ifdef CHILD
      call rnbq_bc_interp_tile(Istr,Iend,Jstr,Jend)
# endif
!
# if defined NBQ_FRC_BRY || defined NBQNUDGING
      tau_in=dt*tauT_in
      tau_out=dt*tauT_out
# endif
!
# define nnew 2
# define nstp 1
!
# ifndef EW_COM_PERIODIC
!
!====================================================================
!                            WESTERN BC
!====================================================================
      if (WESTERN_EDGE) then
#  if defined OBC_COM_WEST && defined OBC_COM_NBQORLANSKI
        do k=1,N                             ! Western edge radiation
          do j=Jstr,Jend+1                   ! ======= ==== =========

            l_b_j     = ijk2lq_nh(Istr-1,j  ,k)
            l_b_jm1   = ijk2lq_nh(Istr-1,j-1,k)
            l_bm1_j   = ijk2lq_nh(Istr  ,j  ,k)
            l_bm1_jm1 = ijk2lq_nh(Istr  ,j-1,k)

            grad(Istr-1,j)=(rhp_nbq_a(l_b_j  )
     &                     -rhp_nbq_a(l_b_jm1))
#   ifdef MASKING
     &                              *vmask(Istr-1,j)
#   endif
            grad(Istr  ,j)=(rhp_nbq_a(l_bm1_j  )
     &                     -rhp_nbq_a(l_bm1_jm1))
#   ifdef MASKING
     &                                  *vmask(Istr,j)
#   endif
          enddo
          do j=Jstr,Jend

            l_b_j     = ijk2lq_nh(Istr-1,j,k)
            l_bm1_j   = ijk2lq_nh(Istr  ,j,k)
            l_bm2_j   = ijk2lq_nh(Istr+1,j,k)

            dft=rhp_nbq_a(l_bm1_j)-rhp_nbq_a(l_bm1_j)
            dfx=rhp_nbq_a(l_bm1_j)-rhp_nbq_a(l_bm2_j)

            if (dfx*dft .lt. 0.) then
              dft=0.                       ! <-- cancel cx, if inflow
#   if defined NBQ_FRC_BRY || defined NBQNUDGING
              tau=tau_in
            else
              tau=tau_out
#   endif
            endif
 
            if (dft*(grad(Istr,j)+grad(Istr,j+1)) .gt. 0.) then
              dfy=grad(Istr,j)
            else
              dfy=grad(Istr,j+1)
            endif
 
#   ifdef OBC_COM_RAD_NORMAL
            dfy=0.
#   endif
            cff=max(dfx*dfx+dfy*dfy, eps)
            cx=dft*dfx
#   ifdef OBC_COM_RAD_NPO
            cy=0.
#   else
            cy=min(cff,max(dft*dfy,-cff))
#   endif
 
            rhp_nbq_a(l_b_j)=( cff*rhp_nbq_a(l_b_j  )
     &                              +cx*rhp_nbq_a(l_bm1_j)
     &                             -max(cy,0.)*grad(Istr-1,j  )
     &                             -min(cy,0.)*grad(Istr-1,j+1)
     &                                               )/(cff+cx)
#   if defined NBQ_FRC_BRY  || defined NBQNUDGING 
            rhp_nbq_a(l_b_j)=(1.-tau)*rhp_nbq_a(l_b_j)
#    ifdef NBQ_FRC_BRY
     &                                    +tau*rnbqbry_west(j,k)
#    else     
     &                                  +tau*rnbqclm(Istr-1,j,k)
#    endif
#   endif
#   ifdef MASKING
            rhp_nbq_a(l_b_j)=rhp_nbq_a(l_b_j)*rmask(Istr-1,j)
#   endif
          enddo
        enddo
!
#  elif defined OBC_COM_WEST && defined OBC_COM_NBQSPECIFIED
!                                            Western edge Specified BC
!                                            ======= ==== ========= ==
        do k=1,N
          do j=Jstr,Jend
            l_b_j     = ijk2lq_nh(Istr-1,j  ,k)
#   ifdef NBQ_FRC_BRY
            rhp_nbq_a(l_b_j)=rnbqbry_west(j,k)         ! specified
#   else
            rhp_nbq_a(l_b_j)=rnbqclm(Istr-1,j,k)
#   endif
#   ifdef MASKING
     &                               *rmask(Istr-1,j)
#   endif
          enddo
        enddo
!#  else
!                                            Western edge gradient BC
!                                            ======= ==== ======== ==
!        do k=1,N
!          do j=Jstr,Jend
!            l_b_j     = ijk2lq_nh(Istr-1,j  ,k)
!            l_bm1_j   = ijk2lq_nh(Istr  ,j  ,k)
!            rhp_nbq_a(l_b_j)=rhp_nbq_a(l_bm1_j) ! Gradient: default
!#   ifdef MASKING
!     &                                     *rmask(Istr,j)
!#   endif
!          enddo
!        enddo
# endif              /* OBC_COM_WEST */
      endif         !<-- WESTERN_EDGE
!
!====================================================================
!                            EASTERN BC
!====================================================================
      if (EASTERN_EDGE) then
#  if defined OBC_COM_EAST && defined OBC_COM_NBQORLANSKI
        do k=1,N                             ! Eastern edge radiation
          do j=Jstr,Jend+1                   ! ======= ==== =========

            l_b_j       = ijk2lq_nh(Iend+1,j  ,k)
            l_b_jm1     = ijk2lq_nh(Iend+1,j-1,k)
            l_bm1_j     = ijk2lq_nh(Iend  ,j  ,k)
            l_bm1_jm1   = ijk2lq_nh(Iend  ,j-1,k)

            grad(Iend  ,j)=(rhp_nbq_a(l_bm1_j  )
     &                     -rhp_nbq_a(l_bm1_jm1))
#   ifdef MASKING
     &                                *vmask(Iend  ,j)
#   endif
            grad(Iend+1,j)=(rhp_nbq_a(l_b_j  )
     &                     -rhp_nbq_a(l_b_jm1))
#   ifdef MASKING
     &                              *vmask(Iend+1,j)
#   endif
          enddo
          do j=Jstr,Jend

            l_b_j     = ijk2lq_nh(Iend+1,j,k)
            l_bm1_j   = ijk2lq_nh(Iend  ,j,k)
            l_bm2_j   = ijk2lq_nh(Iend-1,j,k)

            dft=rhp_nbq_a(l_bm1_j)-rhp_nbq_a(l_bm1_j)
            dfx=rhp_nbq_a(l_bm1_j)-rhp_nbq_a(l_bm2_j)
 
            if (dfx*dft .lt. 0.) then
              dft=0.                      ! <-- cancel cx, if inflow
#   if defined NBQ_FRC_BRY || defined NBQNUDGING
              tau=tau_in
            else
              tau=tau_out
#   endif
            endif
 
            if (dft*(grad(Iend,j)+grad(Iend,j+1)) .gt. 0.) then
              dfy=grad(Iend,j)
            else
              dfy=grad(Iend,j+1)
            endif
 
#   ifdef OBC_COM_RAD_NORMAL
            dfy=0.
#   endif
            cff=max(dfx*dfx+dfy*dfy, eps)
            cx=dft*dfx
#   ifdef OBC_COM_RAD_NPO
            cy=0.
#   else
            cy=min(cff,max(dft*dfy,-cff))
#   endif
 
            rhp_nbq_a(l_b_j)=( cff*rhp_nbq_a(l_b_j  )
     &                              +cx*rhp_nbq_a(l_bm1_j)
     &                             -max(cy,0.)*grad(Iend+1,j  )
     &                             -min(cy,0.)*grad(Iend+1,j+1)
     &                                               )/(cff+cx)
#   if defined NBQ_FRC_BRY  || defined NBQNUDGING 
            rhp_nbq_a(l_b_j)=(1.-tau)*rhp_nbq_a(l_b_j)
#    ifdef NBQ_FRC_BRY
     &                                    +tau*rnbqbry_east(j,k)
#    else     
     &                                  +tau*rnbqclm(Iend+1,j,k)
#    endif
#   endif
#   ifdef MASKING
            rhp_nbq_a(l_b_j)=rhp_nbq_a(l_b_j)*rmask(Iend+1,j)
#   endif
          enddo
        enddo
!
#  elif defined OBC_COM_EAST && defined OBC_COM_NBQSPECIFIED
!                                            Eastern edge Specified BC
!                                            ======= ==== ========= ==
        do k=1,N
          do j=Jstr,Jend
            l_b_j     = ijk2lq_nh(Iend+1,j,k)
#   ifdef NBQ_FRC_BRY
             rhp_nbq_a(l_b_j)=rnbqbry_east(j,k)       ! specified
#   else
             rhp_nbq_a(l_b_j)=rnbqclm(Iend+1,j,k)
#   endif
#   ifdef MASKING
     &                                *rmask(Iend+1,j)
#   endif
          enddo
        enddo
!#  else
!                                            Eastern edge gradient BC
!                                            ======= ==== ======== ==
!        do k=1,N
!          do j=Jstr,Jend
!            l_b_j     = ijk2lq_nh(Iend+1,j,k)
!            l_bm1_j   = ijk2lq_nh(Iend  ,j,k)
!            rhp_nbq_a(l_b_j)=rhp_nbq_a(l_bm1_j)  ! gradient (default)
!#   ifdef MASKING
!     &                                   *rmask(Iend+1,j)
!#   endif
!          enddo
!        enddo
#  endif            /* OBC_COM_EAST */
      endif         !<-- EASTERN_EDGE
# endif             /* !EW_COM_PERIODIC */
 
# ifndef NS_COM_PERIODIC
!
!====================================================================
!                            SOUTHERN BC
!====================================================================
      if (SOUTHERN_EDGE) then
#  if defined OBC_COM_SOUTH && defined OBC_COM_NBQORLANSKI
        do k=1,N                            ! Southern edge radiation
          do i=Istr,Iend+1                  ! ======== ==== =========

            l_b_i     = ijk2lq_nh(i  ,Jstr-1,k)
            l_b_im1   = ijk2lq_nh(i-1,Jstr-1,k)
            l_bm1_i   = ijk2lq_nh(i  ,Jstr  ,k)
            l_bm1_im1 = ijk2lq_nh(i-1,Jstr  ,k)

            grad(i,Jstr-1)=(rhp_nbq_a(l_b_i  )
     &                     -rhp_nbq_a(l_b_im1))
#   ifdef MASKING
     &                              *pmask(i  ,Jstr)
#   endif
            grad(i,Jstr  )=(rhp_nbq_a(l_bm1_i  )
     &                     -rhp_nbq_a(l_bm1_im1))
#   ifdef MASKING
     &                                *pmask(i,Jstr+1)
#   endif
          enddo
          do i=Istr,Iend

            l_b_i     = ijk2lq_nh(i  ,Jstr-1,k)
            l_bm1_i   = ijk2lq_nh(i  ,Jstr  ,k)
            l_bm2_i   = ijk2lq_nh(i  ,Jstr+1,k)

            dft=rhp_nbq_a(l_bm1_i)-rhp_nbq_a(l_bm1_i)
            dfx=rhp_nbq_a(l_bm1_i)-rhp_nbq_a(l_bm2_i)

            if (dfx*dft .lt. 0.) then
              dft=0.                      ! <-- cancel cx, if inflow
#   if defined NBQ_FRC_BRY || defined NBQNUDGING
              tau=tau_in
            else
              tau=tau_out
#   endif
            endif
 
            if (dft*(grad(i,Jstr)+grad(i+1,Jstr)) .gt. 0.) then
              dfy=grad(i,Jstr)
            else
              dfy=grad(i+1,Jstr)
            endif
 
#   ifdef OBC_COM_RAD_NORMAL
            dfy=0.
#   endif
            cff=max(dfx*dfx+dfy*dfy, eps)
            cx=dft*dfx
#   ifdef OBC_COM_RAD_NPO
            cy=0.
#   else
            cy=min(cff,max(dft*dfy,-cff))
#   endif
 
            rhp_nbq_a(l_b_i)=( cff*rhp_nbq_a(l_b_i  )
     &                              +cx*rhp_nbq_a(l_bm1_i)
     &                             -max(cy,0.)*grad(i  ,Jstr-1)
     &                             -min(cy,0.)*grad(i+1,Jstr-1)
     &                                               )/(cff+cx)
#   if defined NBQ_FRC_BRY  || defined NBQNUDGING 
            rhp_nbq_a(l_b_i)=(1.-tau)*rhp_nbq_a(l_b_i)
#    ifdef NBQ_FRC_BRY
     &                                   +tau*rnbqbry_south(i,k)
#    else     
     &                                  +tau*rnbqclm(i,Jstr-1,k)
#    endif
#   endif
#   ifdef MASKING
            rhp_nbq_a(l_b_i)=rhp_nbq_a(l_b_i)*rmask(i,Jstr-1)
#   endif
          enddo
        enddo
#  elif defined OBC_COM_SOUTH && defined OBC_COM_NBQSPECIFIED
!                                           Southern edge Specified BC
!                                           ======== ==== ========= ==
        do k=1,N
          do i=Istr,Iend
            l_b_i     = ijk2lq_nh(i  ,Jstr-1,k)
#   ifdef NBQ_FRC_BRY
            rhp_nbq_a(l_b_i)=rnbqbry_south(i,k)        ! specified
#   else
            rhp_nbq_a(l_b_i)=rnbqclm(i,Jstr-1,k)
#   endif
#   ifdef MASKING
     &                               *rmask(i,Jstr-1)
#   endif
          enddo
        enddo
!#  else
!                                           Southern edge gradient BC
!                                           ======== ==== ======== ==
!        do k=1,N
!          do i=Istr,Iend
!            l_b_i     = ijk2lq_nh(i  ,Jstr-1,k)
!            l_bm1_i   = ijk2lq_nh(i  ,Jstr  ,k)
!            rhp_nbq_a(l_b_i)=rhp_nbq_a(l_bm1_i)  ! gradient (default)
!#   ifdef MASKING
!     &                                   *rmask(i,Jstr-1)
!#   endif
!          enddo
!        enddo
#  endif              /* OBC_COM_SOUTH */
      endif         !<-- SOUTHERN_EDGE
!
!====================================================================
!                            NORTHERN BC
!====================================================================
      if (NORTHERN_EDGE) then
#  if defined OBC_COM_NORTH && defined OBC_COM_NBQORLANSKI
        do k=1,N                            ! Northern edge radiation
          do i=Istr,Iend+1                  ! ======== ==== =========

            l_b_i     = ijk2lq_nh(i  ,Jend+1,k)
            l_b_im1   = ijk2lq_nh(i-1,Jend+1,k)
            l_bm1_i   = ijk2lq_nh(i  ,Jend  ,k)
            l_bm1_im1 = ijk2lq_nh(i-1,Jend  ,k)

            grad(i,Jend  )=(rhp_nbq_a(l_bm1_j  )
     &                     -rhp_nbq_a(l_bm1_im1))
#   ifdef MASKING
     &                                  *umask(i,Jend)
#   endif
            grad(i,Jend+1)=(rhp_nbq_a(l_b_j  )
     &                     -rhp_nbq_a(l_b_im1))
#   ifdef MASKING
     &                              *umask(i,Jend+1)
#   endif
          enddo
          do i=Istr,Iend

            l_b_i     = ijk2lq_nh(i  ,Jend+1,k)
            l_bm1_i   = ijk2lq_nh(i  ,Jend  ,k)
            l_bm2_i   = ijk2lq_nh(i  ,Jend-1,k)

            dft=rhp_nbq_a(l_bm1_i)-rhp_nbq_a(l_bm1_i)
            dfx=rhp_nbq_a(l_bm1_i)-rhp_nbq_a(l_bm2_i)

            if (dfx*dft .lt. 0.) then
              dft=0.                       ! <-- cancel cx, if inflow
#   if defined NBQ_FRC_BRY || defined NBQNUDGING
              tau=tau_in
            else
              tau=tau_out
#   endif
            endif
 
            if (dft*(grad(i,Jend)+grad(i+1,Jend)) .gt. 0.) then
              dfy=grad(i,Jend)
            else
              dfy=grad(i+1,Jend)
            endif
 
#   ifdef OBC_COM_RAD_NORMAL
            dfy=0.
#   endif
            cff=max(dfx*dfx+dfy*dfy, eps)
            cx=dft*dfx
#   ifdef OBC_COM_RAD_NPO
            cy=0.
#   else
            cy=min(cff,max(dft*dfy,-cff))
#   endif
 
            rhp_nbq_a(l_b_i)=( cff*rhp_nbq_a(l_b_i)
     &                              +cx*rhp_nbq_a(l_bm1_i)
     &                             -max(cy,0.)*grad(i  ,Jend+1)
     &                             -min(cy,0.)*grad(i+1,Jend+1)
     &                                               )/(cff+cx)
#   if defined NBQ_FRC_BRY  || defined NBQNUDGING 
            rhp_nbq_a(l_b_i)=(1.-tau)*rhp_nbq_a(l_b_i)
#    ifdef NBQ_FRC_BRY
     &                                   +tau*rnbqbry_north(i,k)
#    else     
     &                                  +tau*rnbqclm(i,Jend+1,k)
#    endif
#   endif
#   ifdef MASKING
            rhp_nbq_a(l_b_i)=rhp_nbq_a(l_b_i)*rmask(i,Jend+1)
#   endif
          enddo
        enddo
!
#  elif defined OBC_COM_NORTH && defined OBC_COM_NBQSPECIFIED
!                                           Northern edge Specified BC
!                                           ======== ==== ========= ==
        do k=1,N
          do i=Istr,Iend
            l_b_i     = ijk2lq_nh(i  ,Jend+1,k)
#   ifdef NBQ_FRC_BRY
            rhp_nbq_a(l_b_i)=rnbqbry_north(i,k)      ! specified
#   else
            rhp_nbq_a(l_b_i)=rnbqclm(i,Jend+1,k)
#   endif
#   ifdef MASKING
     &                               *rmask(i,Jend+1)
#   endif
          enddo
        enddo
!#  else
!        do k=1,N
!          do i=Istr,Iend
!                                           Northern edge gradient BC
!                                           ======== ==== ======== ==
!            l_b_i     = ijk2lq_nh(i  ,Jend+1,k)
!            l_bm1_i   = ijk2lq_nh(i  ,Jend  ,k)
!            rhp_nbq_a(l_b_i)=rhp_nbq_a(l_bm1_i) ! gradient (default)
!#   ifdef MASKING
!     &                                   *rmask(i,Jend+1)
!#   endif
!          enddo
!        enddo
#  endif
      endif     !<--  NORTHERN_EDGE
# endif          /* !NS_COM_PERIODIC */
 
                           ! Corners between adjacent open boundaries
                           ! ======= ======= ======== ==== ==========
 
# if defined OBC_COM_SOUTH && defined OBC_COM_WEST
      if (WESTERN_EDGE .and. SOUTHERN_EDGE) then
        do k=1,N
          l_c     = ijk2lq_nh(Istr-1,Jstr-1,k)
          l_c_im1 = ijk2lq_nh(Istr  ,Jstr-1,k)
          l_c_jm1 = ijk2lq_nh(Istr-1,Jstr  ,k)
          rhp_nbq_a(l_c)=0.5*(rhp_nbq_a(l_c_im1)
     &                            +rhp_nbq_a(l_c_jm1))
#  ifdef MASKING
     &                                *rmask(Istr-1,Jstr-1)
#  endif
        enddo
      endif
# endif
# if defined OBC_COM_SOUTH && defined OBC_COM_EAST
      if (EASTERN_EDGE .and. SOUTHERN_EDGE) then
        do k=1,N
          l_c     = ijk2lq_nh(Iend+1,Jstr-1,k)
          l_c_im1 = ijk2lq_nh(Iend  ,Jstr-1,k)
          l_c_jm1 = ijk2lq_nh(Iend+1,Jstr  ,k)
          rhp_nbq_a(l_c)=0.5*(rhp_nbq_a(l_c_im1)
     &                            +rhp_nbq_a(l_c_jm1))
#  ifdef MASKING
     &                                *rmask(Iend+1,Jstr-1)
#  endif
        enddo
      endif
# endif
# if defined OBC_COM_NORTH && defined OBC_COM_WEST
      if (WESTERN_EDGE .and. NORTHERN_EDGE) then
        do k=1,N
          l_c     = ijk2lq_nh(Istr-1,Jend+1,k)
          l_c_im1 = ijk2lq_nh(Istr  ,Jend+1,k)
          l_c_jm1 = ijk2lq_nh(Istr-1,Jend  ,k)
          rhp_nbq_a(l_c)=0.5*(rhp_nbq_a(l_c_im1)
     &                            +rhp_nbq_a(l_c_jm1))
#  ifdef MASKING
     &                                *rmask(Istr-1,Jend+1)
#  endif
        enddo
      endif
# endif
# if defined OBC_COM_NORTH && defined OBC_COM_EAST
      if (EASTERN_EDGE .and. NORTHERN_EDGE) then
        do k=1,N
          l_c     = ijk2lq_nh(Iend+1,Jend+1,k)
          l_c_im1 = ijk2lq_nh(Iend  ,Jend+1,k)
          l_c_jm1 = ijk2lq_nh(Iend+1,Jend  ,k)
          rhp_nbq_a(l_c)=0.5*(rhp_nbq_a(l_c_im1)
     &                            +rhp_nbq_a(l_c_jm1))
#  ifdef MASKING
     &                                *rmask(Iend+1,Jend+1)
#  endif
        enddo
      endif
# endif

# undef nnew
# undef nstp

      return
      end
#else
# ifndef CHILD
      subroutine rnbq_bc_parent_empty
      end
# else
      subroutine rnbq_bc_child_empty
      end
# endif
#endif /* NBQ */
#ifndef CHILD
# define CHILD
# ifdef AGRIF
#  include "rnbq_bc.F"
# endif
# undef CHILD
#endif  /* !CHILD */
