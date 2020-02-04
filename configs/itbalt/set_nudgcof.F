! $Id: set_nudgcof.F 1458 2014-02-03 15:01:25Z gcambon $
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
#include "cppdefs.h"
#if defined ZNUDGING || defined M2NUDGING || defined M3NUDGING\
    || defined TNUDGING || defined NBQ_NUDGING || defined SPONGE

      subroutine set_nudgcof (tile)
      implicit none
      integer tile, trd, omp_get_thread_num
# include "param.h"
# include "private_scratch.h"
# include "compute_tile_bounds.h"
      trd=omp_get_thread_num()
      call set_nudgcof_tile (Istr,Iend,Jstr,Jend,A2d(1,1,trd))
      return
      end
!
      subroutine set_nudgcof_tile (Istr,Iend,Jstr,Jend, wrk)
!     implicit none
      use params
      integer ierr
# if defined MPI
      include 'mpif.h'
# endif
# include "param.h"
# include "grid.h"
# include "climat.h"
# include "mixing.h"
# include "scalars.h"
# include "ocean3d.h"
# include "mpi_cpl.h"
# ifdef NBQ_NUDGING
#  include "nbq.h"
# endif
      integer Istr,Iend,Jstr,Jend, i, j, k, isp, itrc, ibnd
      real    wrk(PRIVATE_2D_SCRATCH_ARRAY)
# if defined CANARY  || defined IBERIA
      real    lon0,lat0, rsponge, dx,dy,dr,cff,cff1,cff2
# endif

! aponte itide (
# if defined ITIDE
      real y
      real cff,cff1,cff2, cff3, cff3_itide, cff3_zano, cff4
      real invday
      parameter (invday=1.0/86400.0)
! complex sponge formulation
      integer N_spg
      parameter (N_spg=5)
      real y_spg(1:N_spg)
      real tau_spg(1:N_spg)
      integer i_spg
# endif
      real cff_kill_sponges
! aponte itide )
! aponte  (
#ifdef Y_STRETCH_RELAX
      real rl_start_south, rl_start_north
#endif
! aponte  )
!
# include "compute_extended_bounds.h"
!
!--------------------------------------------------------------------
!  Set up nudging towards climatology time-scale coeffcients [1/s]
!  and/or increase horizontal mixing in the sponge areas.
!--------------------------------------------------------------------
!
# if defined SPONGE_GRID || !defined SPONGE
#  ifdef IGW
      isp=155
#  else
      isp=10
#  endif
# else
#  ifdef MPI
      if (ii.eq.0.and.jj.eq.0.and.Istr.eq.1.and.Jstr.eq.1) then
        isp=int(x_sponge*pm(1,1))
      endif      
      call MPI_Bcast(isp, 1, MPI_INTEGER,
     &                    0, MPI_COMM_WORLD, ierr)      
#  else
      isp=int(x_sponge*pm(1,1)) ! number of points in layer     
#  endif
# endif /* SPONGE_GRID */

!
      do j=max(-1,JstrR-1),JendR
        do i=max(-1,IstrR-1),IendR
          ibnd=isp
# ifndef MPI
#  ifdef OBC_WEST
          ibnd=min(ibnd,i)
#  endif
#  ifdef OBC_EAST
          ibnd=min(ibnd,Lm+1-i)
#  endif
#  ifdef OBC_SOUTH
          ibnd=min(ibnd,j)
#  endif
#  ifdef OBC_NORTH
          ibnd=min(ibnd,Mm+1-j)
#  endif
# else
#  ifdef OBC_WEST
          ibnd=min(ibnd,i+iminmpi-1)          
#  endif
#  ifdef OBC_EAST
          ibnd=min(ibnd,LLm+1-(i+iminmpi-1))          
#  endif
#  ifdef OBC_SOUTH
          ibnd=min(ibnd,j+jminmpi-1)          
#  endif
#  ifdef OBC_NORTH
          ibnd=min(ibnd,MMm+1-(j+jminmpi-1))          
#  endif
# endif

# ifdef IGW
          if (ibnd .ge. isp) then
            wrk(i,j)=0.
          else
            wrk(i,j)=exp(-5.d-5*float(ibnd)*grdmax)       ! exp. profile
          endif
# else
          wrk(i,j)=.5*(cos(pi*float(ibnd)/float(isp))+1.) ! cosine profile
!         wrk(i,j)=float(isp-ibnd)/float(isp)             ! linear profile
# endif

! aponte jet (
! this piece of code should be commented and ROBUST_DIAG should be turned on 
# ifdef ITBALT
          wrk(i,j)=1.
# endif
! aponte jet )

# ifdef SPONGE_SED
          if (i.gt.0 .and. j.gt.0)
     &      cff_sponge(i,j)=float(ibnd)/float(isp)
# endif
        enddo
      enddo
!
!-------------------------------------------------------------------
! Compute nudging coefficients in nudging layers
!-------------------------------------------------------------------
!
      do j=JstrR,JendR
        do i=IstrR,IendR
# ifdef TNUDGING
          Tnudgcof(i,j,N,itemp)=tauT_out*wrk(i,j)
# endif
# ifdef ZNUDGING
          Znudgcof(i,j)=tauM_out*wrk(i,j)
# endif
# ifdef M2NUDGING
          M2nudgcof(i,j)=tauM_out*wrk(i,j)
# endif
# ifdef M3NUDGING
          M3nudgcof(i,j)=tauM_out*wrk(i,j)
# endif
# ifdef NBQ_NUDGING
          NBQnudgcof(i,j)=dtfast*tauM_out*wrk(i,j) ! nudg. to internal modes (dim. less) 
# endif
        enddo

! aponte itide (
# ifdef ITIDE
#  ifdef MPI
         y=j+jj*Mm
#  else
         y=j
#  endif
        nudgcof_itide(j) = 0.
        y=y/real(MMm)

!!!
#ifdef Y_STRETCH_RELAX

       ! compute relaxation for south and north stretched grid
        if ( rl_flag_south ) then
           rl_start_south = 0.5*el - rl_dist_south        ! start south from the middle
           cff = rl_L_south/el  			  ! nondim width of slope
           cff1 = invday/rl_tau_in_south  		  ! inverse relax time: inside value
           cff2 = invday/rl_tau_out_south  		  ! inverse relax time: outside value
           if ( yp(i,j)<=(rl_start_south-rl_L_south))then ! plateau
              nudgcof_itide(j) = cff2
           endif
           if ( yp(i,j)>(rl_start_south-rl_L_south) .and. 
     &          yp(i,j) <= rl_start_south ) then          ! slope  
              nudgcof_itide(j) = 
     &              -(cff2-cff1)/(cff) * (yp(i,j) - rl_start_south)/el + cff1         
           endif
        endif
 
        if ( rl_flag_north ) then
           rl_start_north = 0.5*el + rl_dist_north   	  ! start north from the middle     
           cff = rl_L_north/el  			  ! nondim width of slope
           cff1 = invday/rl_tau_in_north  		  ! inverse relax time: south
           cff2 = invday/rl_tau_out_north  		  ! inverse relax time: north 
           if ( yp(i,j)>= (rl_start_north) .and. 
     &          yp(i,j) < (rl_start_north+rl_L_north) )  then  !slope
              nudgcof_itide(j) =  
     &               (cff2-cff1)/(cff) * (yp(i,j) - rl_start_north)/el + cff1           
           endif
           if ( yp(i,j) >= (rl_start_north+rl_L_north)) then   ! plateau
              nudgcof_itide(j) = cff2 
           end if  

        endif

#else /* Y_STRETCH_RELAX */

! Flag that kills all sponges but not the zonal relaxation (0=no sponges)
        cff_kill_sponges = 0.0
        call params_output('flag_sponges', cff_kill_sponges)
!        if (mynode==0) write(*,*) 'cff_kill_sponges = ', cff_kill_sponges

        if ( rl_flag_south ) then

#   if 0
! original formulation with linear increase plus plateau
         cff = rl_L_south/el  ! nondim width
         cff1 = invday/rl_tau_in_south  ! inverse relax time: inside value
         cff2 = invday/rl_tau_out_south  ! inverse relax time: outside value

         if ( y .lt. cff ) then
          cff3_itide = cff1 * y/cff
     &         + 2.0*cff2 * (cff-y)/cff ! linear + plateau
!     &         + cff2 * (cff-y)/cff ! linear
          cff3_itide = max(cff3_itide,cff1)
          cff3_itide = min(cff3_itide,cff2)
!          nudgcof_itide(j) = nudgcof_itide(j) + cff3_itide
          nudgcof_itide(j) = cff3_itide
         endif
#   else
! damp zonal mean + anomalies 
! more complex formulation
        y_spg(1)=0.0D+3
        y_spg(2)=200.0D+3
        y_spg(3)=300.0D+3
        y_spg(4)=500.0D+3
        y_spg(5)=800.0D+3
        tau_spg(1)=0.5D+0
        tau_spg(2)=0.5D+0
        tau_spg(3)=10.0D+0
        tau_spg(4)=100.0D+0
        tau_spg(5)=100.0D+0
        do i_spg=1,N_spg
         y_spg(i_spg)=y_spg(i_spg)/el
         tau_spg(i_spg)=tau_spg(i_spg)/invday
        enddo

        cff3_itide=0.0
        do i_spg=1,N_spg-1
        if ( (y .ge. y_spg(i_spg)) .and. (y .le. y_spg(i_spg+1)) ) then
         cff = y_spg(i_spg+1)-y_spg(i_spg)  ! width
         cff1=1.0/tau_spg(i_spg)
         cff2=1.0/tau_spg(i_spg+1)
         cff3_itide=   cff1*(y_spg(i_spg+1)-y)/cff 
     &         + cff2*(y-y_spg(i_spg))/cff
        endif
        enddo
        nudgcof_itide(j) = cff3_itide
! kills efficiently this type of damping
        nudgcof_itide(j) = nudgcof_itide(j)*cff_kill_sponges


! damp zonal anomalies 
! more complex formulation
        y_spg(1)=0.0D+3
        y_spg(2)=100.0D+3
        y_spg(3)=300.0D+3
        y_spg(4)=600.0D+3
        y_spg(5)=800.0D+3
        tau_spg(1)=10.0D+0
        tau_spg(2)=10.0D+0
        tau_spg(3)=10.0D+0
        tau_spg(4)=10.0D+0
!       tau_spg(1)=1.0D+0
!       tau_spg(2)=1.0D+0
!       tau_spg(3)=1.0D+0
!       tau_spg(4)=1.0D+0
        tau_spg(5)=100.0D+0
        do i_spg=1,N_spg
         y_spg(i_spg)=y_spg(i_spg)/el
         tau_spg(i_spg)=tau_spg(i_spg)/invday
        enddo

        cff3_zano=0.0
        do i_spg=1,N_spg-1
        if ( (y .ge. y_spg(i_spg)) .and. (y .le. y_spg(i_spg+1)) ) then
         cff = y_spg(i_spg+1)-y_spg(i_spg)  ! width
         cff1=1.0/tau_spg(i_spg)
         cff2=1.0/tau_spg(i_spg+1)
         cff3_zano=   cff1*(y_spg(i_spg+1)-y)/cff
     &         + cff2*(y-y_spg(i_spg))/cff
        endif
        enddo
        nudgcof_zano_itide(j) = cff3_zano
! kills efficiently this type of damping
        nudgcof_zano_itide(j) = nudgcof_zano_itide(j)*cff_kill_sponges

#  endif

        endif

!!!
        if ( rl_flag_north ) then

#   if 1

         cff = rl_L_north/el  ! nondim width
         cff1 = invday/rl_tau_in_north  ! inverse relax time: south
         cff2 = invday/rl_tau_out_north  ! inverse relax time: north
!
         if ( y .gt. 1-cff ) then
          cff3_itide = cff1 * (1.0-y)/cff
     &         + 2.0*cff2 * (y-1.0+cff)/cff ! linear + plateau
!     &         + cff2 * (y-1.0+cff)/cff ! linear
          cff3_itide = max(cff3_itide,cff1)
          cff3_itide = min(cff3_itide,cff2)
!          nudgcof_itide(j) = nudgcof_itide(j) + cff3_itide
          nudgcof_itide(j) = cff3_itide
         endif
! kills efficiently this type of damping
        nudgcof_itide(j) = nudgcof_itide(j)*cff_kill_sponges 

#   else
! damp zonal mean + anomalies 
! more complex formulation
        y_spg(1)= el - 0.0D+3
        y_spg(2)= el - 200.0D+3
        y_spg(3)= el - 300.0D+3
        y_spg(4)= el - 500.0D+3
        y_spg(5)= el - 800.0D+3
        tau_spg(1)=0.5D+0
        tau_spg(2)=0.5D+0
        tau_spg(3)=10.0D+0
        tau_spg(4)=100.0D+0
        tau_spg(5)=100.0D+0
        do i_spg=1,N_spg
         y_spg(i_spg)=y_spg(i_spg)/el
         tau_spg(i_spg)=tau_spg(i_spg)/invday
        enddo

!       cff3=0.0
        do i_spg=1,N_spg-1
        if ( (y .le. y_spg(i_spg)) .and. (y .ge. y_spg(i_spg+1)) ) then
         cff = y_spg(i_spg+1)-y_spg(i_spg)  ! width
         cff1=1.0/tau_spg(i_spg)
         cff2=1.0/tau_spg(i_spg+1)
         cff3_itide=   cff1*(y_spg(i_spg+1)-y)/cff 
     &         + cff2*(y-y_spg(i_spg))/cff
        endif
        enddo
        nudgcof_itide(j) = cff3_itide
! kills efficiently this type of damping
        nudgcof_itide(j) = nudgcof_itide(j)*cff_kill_sponges


! damp zonal anomalies 
! more complex formulation
        y_spg(1)= el - 0.0D+3
        y_spg(2)= el - 100.0D+3
        y_spg(3)= el - 300.0D+3
        y_spg(4)= el - 600.0D+3
        y_spg(5)= el - 800.0D+3
        tau_spg(1)=10.0D+0
        tau_spg(2)=10.0D+0
        tau_spg(3)=10.0D+0
        tau_spg(4)=10.0D+0
!       tau_spg(1)=1.0D+0
!       tau_spg(2)=1.0D+0
!       tau_spg(3)=1.0D+0
!       tau_spg(4)=1.0D+0
        tau_spg(5)=100.0D+0
        do i_spg=1,N_spg
         y_spg(i_spg)=y_spg(i_spg)/el
         tau_spg(i_spg)=tau_spg(i_spg)/invday
        enddo

!       cff3=0.0
        do i_spg=1,N_spg-1
        if ( (y .le. y_spg(i_spg)) .and. (y .ge. y_spg(i_spg+1)) ) then
         cff = y_spg(i_spg+1)-y_spg(i_spg)  ! width
         cff1=1.0/tau_spg(i_spg)
         cff2=1.0/tau_spg(i_spg+1)
         cff3_zano=   cff1*(y_spg(i_spg+1)-y)/cff
     &         + cff2*(y-y_spg(i_spg))/cff
        endif
        enddo
        nudgcof_zano_itide(j) = cff3_zano
! kills efficiently this type of damping
        nudgcof_zano_itide(j) = nudgcof_zano_itide(j)*cff_kill_sponges

#  endif

        endif

#endif /* Y_STRETCH_RELAX */

!
!        y=y/real(MMm)
!!        cff = 0.3   ! nondim width
!        cff = 1. * 175000.0/el  ! nondim width
!        cff1 = 1.0/(86400.0*100.0)  ! inverse relax time: south
!        cff2 = 1.0/(86400.0*0.5)  ! inverse relax time: north 
!        nudgcof_itide(j) =
!     &          cff1 * (1.0-y)/cff
!     &        + cff2 * (y-1.0+cff)/cff
!        if ( nudgcof_itide(j) .lt. cff1 ) then
!         nudgcof_itide(j) = cff1
!        elseif ( nudgcof_itide(j) .gt. cff2 ) then
!         nudgcof_itide(j) = cff2
!        endif
!!        if ( (y .gt. (1.0-cff)) .and. (y .lt. 1.0) ) then
!!         nudgcof_itide(j) = 1.0 /(
!!     &          86400.0*100.0 * (1.0-y)/cff
!!     &        + 86400.0*1.0   * (y-1.0+cff)/cff
!!!     &         (1.0 + tanh((y-0.8)/0.2))*0.5
!!     &          )
!!        elseif (y .gt. 1.0) then
!!          nudgcof_itide(j) = 1.0 /( 86400.0*1.0 )
!!        else
!!          nudgcof_itide(j) = 1.0 /( 86400.0*100.0 )
!!        endif
!!        nudgcof_itide(j) = 0.0
!!        write(*,*) 'nudgcof_itide = ', nudgcof_itide(j)
# endif /* ITIDE */
! aponte itide )

      enddo
!
!  Apply nudging to other tracers and vertical levels
!
# ifdef TNUDGING
      do itrc=1,NT       ! includes BIOLOGY variables
        do k=1,N
          do j=JstrR,JendR
            do i=IstrR,IendR
              Tnudgcof(i,j,k,itrc)=Tnudgcof(i,j,N,itemp)
            enddo
          enddo
        enddo
      enddo
# endif
!
! Interior nudging
!
# ifdef ROBUST_DIAG
#  ifdef TNUDGING 
      do k=1,N 
        do j=JstrR,JendR
          do i=IstrR,IendR
            Tnudgcof(i,j,k,itemp)=tauT_out
#   ifdef SALINITY
            Tnudgcof(i,j,k,isalt)=tauT_out
#   endif
          enddo
        enddo
      enddo
#  endif /* TNUDGING */
#  ifdef JET
      do j=JstrR,JendR
        do i=IstrR,IendR
#   ifdef ZNUDGING
          Znudgcof(i,j)=tauM_out
#   endif
#   ifdef M2NUDGING
          M2nudgcof(i,j)=tauM_out
#   endif
#   ifdef M3NUDGING
          M3nudgcof(i,j)=tauM_out
#   endif
        enddo
      enddo
#   endif /* JET */
#  endif /* ROBUST_DIAG */

# ifdef SPONGE
!
!-------------------------------------------------------------------
!  Add Viscosity and Diffusivity in SPONGE layers
!-------------------------------------------------------------------
!
! aponte sponge (
#  if defined ITIDE
      do j=max(-1,JstrR-1),JendR
#   ifdef MPI
         y=j+jj*Mm
#   else
         y=j
#   endif
        cff4 = 0.
        y=y/real(MMm)
!
        if ( sp_flag_south ) then
         cff = sp_L_south/el  ! nondim width
         cff1 = invday/sp_tau_in_south  ! inverse relax time: inside value
         cff2 = invday/sp_tau_out_south  ! inverse relax time: outside value
!
         if ( y .lt. cff ) then
          cff3 = cff1 * y/cff
     &         + 2.0*cff2 * (cff-y)/cff ! linear + plateau
!     &         + cff2 * (cff-y)/cff ! linear
          cff3 = max(cff3,cff1)
          cff3 = min(cff3,cff2)
          cff4 = cff3
         endif
        endif
!
        if ( sp_flag_north ) then
         cff = sp_L_north/el  ! nondim width
         cff1 = invday/sp_tau_in_north  ! inverse relax time: south
         cff2 = invday/sp_tau_out_north  ! inverse relax time: north
!
         if ( y .gt. 1-cff ) then
          cff3 = cff1 * (1.0-y)/cff
     &         + 2.0*cff2 * (y-1.0+cff)/cff ! linear + plateau
!     &         + cff2 * (y-1.0+cff)/cff ! linear
          cff3 = max(cff3,cff1)
          cff3 = min(cff3,cff2)
          cff4 = cff3
         endif
        endif

        do i=max(-1,IstrR-1),IendR
! normalize by 1/k^2 to convert this relaxation time scale into
! a viscous variable
          wrk(i,j) = cff4/(2.0*3.1415/175000.0)**2
        enddo
      enddo
#ifdef ITBALT
! SPONGE_GRID is turned on automatically in set_global_definitions.h 
# undef SPONGE_GRID
      v_sponge=1.0
#endif

#  endif
! aponte sponge )
#  if defined UV_VIS2 || defined UV_VIS4 || defined SPONGE_VIS2
!
!  Add Viscosity at rho points 
!
      do j=JstrR,JendR
        do i=IstrR,IendR
#   ifdef SPONGE_GRID
          visc2_sponge_r(i,j)=(0.01 /(pm(i,j)*pn(i,j)*dt))*wrk(i,j)
#   else
          visc2_sponge_r(i,j)=v_sponge*wrk(i,j)
#   endif
#   ifdef UV_VIS4
          visc4_sponge_r(i,j)=visc2_sponge_r(i,j)/(pm(i,j)*pn(i,j))
#   endif
        enddo
      enddo
!
!  Interpolate Viscosity at psi points
!
      do j=Jstr,JendR
        do i=Istr,IendR
#   ifdef SPONGE_GRID
          visc2_sponge_p(i,j)=0.25*(0.01 /(pm(i,j)*pn(i,j)*dt))*
#   else
          visc2_sponge_p(i,j)=0.25*v_sponge*
#   endif
     &                              ( wrk(i,j  )+wrk(i-1,j  )
     &                               +wrk(i,j-1)+wrk(i-1,j-1) ) 
#   ifdef UV_VIS4
          visc4_sponge_p(i,j)=visc2_sponge_p(i,j)/(pm(i,j)*pn(i,j))
#   endif
        enddo
      enddo
#  endif /* UV_VIS2 || UV_VIS4 || SPONGE_VIS2 */

#  if defined TS_DIF2 || defined TS_DIF4 || defined SPONGE_DIF2
!
!  Add Diffusivity for all tracers
!
       do itrc=1,NT
        do j=JstrR,JendR
          do i=IstrR,IendR
#   ifdef SPONGE_GRID
            diff2_sponge(i,j)=(0.01 /(pm(i,j)*pn(i,j)*dt))*wrk(i,j)
#   else
            diff2_sponge(i,j)=v_sponge*wrk(i,j) 
#   endif
#   ifdef TS_DIF4
            diff4_sponge(i,j)=diff2_sponge(i,j)/(pm(i,j)*pn(i,j))
#   endif
        enddo
       enddo
      enddo
#  endif /* TS_DIF2 || TS_DIF4 || SPONGE_DIF2 */
# endif /* SPONGE */          
!
!-------------------------------------------------------------------
!  add configuration specific stuff
!-------------------------------------------------------------------
!
# if defined CANARY || defined IBERIA
      lon0=-8.2
      lat0=35.7
      rsponge=300.e3
      cff=1./(50.*86400.)
      cff1=Eradius*cos(lat0*deg2rad)*deg2rad
      cff2=Eradius*deg2rad
      do j=JstrR,JendR
        do i=IstrR,IendR
          dx=cff1*(lonr(i,j)-lon0)
          dy=cff2*(latr(i,j)-lat0)
          dr=sqrt(dx**2+dy**2)
          if (dr .lt. rsponge) then
            do k=1,N
              Tnudgcof(i,j,k,itemp)=.5*cff*(cos(pi*dr/rsponge)+1)
     &                    *(-atan((z_r(i,j,k)+750.)*2.e-2)/pi+.5)
              Tnudgcof(i,j,k,isalt)=Tnudgcof(i,j,k,itemp)
            enddo
          endif
        enddo
      enddo
# endif
#else
      subroutine set_nudgcof_empty
#endif /* TNUDGING || ZNUDGING || SPONGE */
      return
      end
