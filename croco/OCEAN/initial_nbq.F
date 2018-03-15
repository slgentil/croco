! $Id:$
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
#ifdef NBQ
!
      subroutine initial_nbq (tile, icall)
!
      implicit none
      integer tile, icall, trd
#include "param.h"
#include "private_scratch.h"
!$    integer omp_get_thread_num
#include "compute_tile_bounds.h"
      trd=0
!$    trd=omp_get_thread_num()

      call initial_nbq_tile (icall, 
     &                      Istr,Iend,Jstr,Jend
     &                     ,A3d(1,1,trd), A3d(1,2,trd)
     &                     ,A3d(1,3,trd), A3d(1,4,trd)
#ifdef NBQ_ZETAW
     &                     ,A3d(1,5,trd), A3d(1,6,trd)
#endif
#ifdef NONLIN_EOS
     &                     ,A2d(1,1,trd), A2d(1,2,trd)
#endif
     &                      )
  
      end subroutine initial_nbq
  
      subroutine initial_nbq_tile (icall, 
     &                            Istr,Iend,Jstr,Jend
     &                           ,Hzw_half_nbq_inv, Hzr_half_nbq_inv
     &                           ,Hzw_half_nbq_inv_u,Hzw_half_nbq_inv_v
#ifdef NBQ_ZETAW
     &                           ,Hzu_half_qdmu,Hzv_half_qdmv
#endif
#ifdef NONLIN_EOS
     &                           ,K_up, K_dw
#endif
     &                         )
!
!======================================================================
!
!                 NH and NBQ initialization
!
!======================================================================
!
      implicit none
      integer  Istr, Iend, Jstr, Jend
# ifdef MPI      
      include 'mpif.h'
# endif      
# include "param.h"
# include "scalars.h"
# include "private_scratch.h"
# include "nbq.h"
# include "work.h"
# include "grid.h"
# include "ocean2d.h"
# include "ocean3d.h"

# ifdef NBQ_ZETAW
       real Hzu_half_qdmu(PRIVATE_2D_SCRATCH_ARRAY,0:N)
       real Hzv_half_qdmv(PRIVATE_2D_SCRATCH_ARRAY,0:N)
# endif
       real Hzw_half_nbq_inv(PRIVATE_2D_SCRATCH_ARRAY,0:N)
       real Hzr_half_nbq_inv(PRIVATE_2D_SCRATCH_ARRAY,N)
       real Hzw_half_nbq_inv_u(PRIVATE_2D_SCRATCH_ARRAY,0:N)
       real Hzw_half_nbq_inv_v(PRIVATE_2D_SCRATCH_ARRAY,0:N)
# ifndef NBQ_ZETAW	   	   
       real work3d_nbq(PRIVATE_2D_SCRATCH_ARRAY,N,5)
# endif
# ifdef NONLIN_EOS      
      real K_up(PRIVATE_1D_SCRATCH_ARRAY,0:N)  ! work arrays for call
      real K_dw(PRIVATE_1D_SCRATCH_ARRAY,0:N)  ! to nonlinear EOS
# endif  

      integer i,j,k,ierr
      integer icall

# include "compute_extended_bounds.h"
!
# if defined EW_PERIODIC && !defined MPI
#  define IR_RANGE Istr,Iend
#  define IU_RANGE Istr,Iend
# else
#  define IR_RANGE IstrR,IendR
#  define IU_RANGE Istr,IendR
# endif
# if defined NS_PERIODIC && !defined MPI
#  define JR_RANGE Jstr,Jend
#  define JV_RANGE Jstr,Jend
# else
#  define JR_RANGE JstrR,JendR
#  define JV_RANGE Jstr,JendR
# endif

      if (icall.eq.1) then
!**********************************************************************
!
!                 Initializations (PART I)
!
!**********************************************************************

!----------------------------------------------------------------------
!........Semi-implicit scheme (0/1):
!----------------------------------------------------------------------
# ifdef NBQ_IMP
        ifl_imp_nbq = 1
        MPI_master_only write(6,*)
        MPI_master_only write(6,*) '--------------------------------'
        MPI_master_only write(6,*) ' NBQ: semi-implicit integration '
        MPI_master_only write(6,*) '--------------------------------'
        MPI_master_only write(6,*)
# else
        ifl_imp_nbq = 0
        MPI_master_only write(6,*)
        MPI_master_only write(6,*) '---------------------------'
        MPI_master_only write(6,*) ' NBQ: explicit integration '
        MPI_master_only write(6,*) '---------------------------'
        MPI_master_only write(6,*)
# endif
!
!----------------------------------------------------------------------
!  Initialize density perturbation and momentum arrays
!----------------------------------------------------------------------
!
        do k=1,N
          do j=JR_RANGE
            do i=IU_RANGE
              qdmu_nbq(i,j,k)=0.
            enddo
          enddo
        enddo
        do k=1,N
          do j=JV_RANGE
            do i=IR_RANGE
              qdmv_nbq(i,j,k)=0.
            enddo
          enddo
        enddo
        do k=1,N
          do j=JR_RANGE
            do i=IR_RANGE
              qdmw_nbq(i,j,k)=0.
              rho_nbq(i,j,k) =0.
            enddo
          enddo
        enddo
# ifdef NBQ_MASS
        do j=JR_RANGE
          do i=IR_RANGE
            rhobar_nbq(i,j,k) =1.
          enddo
        enddo
# endif
!
!----------------------------------------------------------------------
!  Initialize parameters: should be done in a NH-namelist
!
!   pseudoacoustic speed should be around 5 times 
!   external phase speed sqrt(g*h)
!----------------------------------------------------------------------
!
        ifl_nbq  = 1
        slip_nbq = 0
        iteration_nbq_max=ndtnbq
        soundspeed_nbq   =csound_nbq        ! pseudoacoustic speed
        soundspeed2_nbq  =csound_nbq**2
!
!----------------------------------------------------------------------
!... Initializes Grid-coef
!----------------------------------------------------------------------
!
        do k=1,N
          do j=JR_RANGE
            do i=IR_RANGE
              Hzw_half_nbq_inv(i,j,k)  =1.e-30
              Hzr_half_nbq_inv(i,j,k)  =1.e-30
              Hzw_half_nbq_inv_u(i,j,k)=1.e-30
              Hzw_half_nbq_inv_v(i,j,k)=1.e-30
            enddo
          enddo
        enddo
# ifdef NBQ_ZETAW
        do k=1,N
          do j=JR_RANGE
            do i=IU_RANGE
              Hzu_half_qdmu(i,j,k)=0.
            enddo
          enddo
        enddo
        do k=1,N
          do j=JV_RANGE
            do i=IR_RANGE
              Hzv_half_qdmv(i,j,k)=0.
            enddo
          enddo
        enddo
        call grid_nbq_tile(Istr,Iend,Jstr,Jend
     &                    ,Hzw_half_nbq_inv,Hzr_half_nbq_inv
     &                    ,Hzw_half_nbq_inv_u,Hzw_half_nbq_inv_v
     &                    ,Hzu_half_qdmu, Hzv_half_qdmv
     &                    )
# else
        call grid_nbq_tile(Istr,Iend,Jstr,Jend
     &                    ,Hzw_half_nbq_inv,Hzr_half_nbq_inv
     &                    ,Hzw_half_nbq_inv_u,Hzw_half_nbq_inv_v
     &                    ,work3d_nbq(PRIVATE_2D_SCRATCH_ARRAY,1,1)
     &                    ,work3d_nbq(PRIVATE_2D_SCRATCH_ARRAY,1,2)
     &                    )
# endif
!
!----------------------------------------------------------------------
!... Set second viscosity coefficient:
!----------------------------------------------------------------------
!
!       csvisc1_nbq  = dtnbq * soundspeed2_nbq + visc2_nbq
!       csvisc2_nbq  = dtnbq * soundspeed2_nbq / csvisc1_nbq 
!
      endif     ! icall == 1

      if (icall.eq.2) then

!**********************************************************************
!
!                 rhobar_nbq initializations (PART II)
!
!**********************************************************************
!.........Grid update:
# ifdef NBQ_PRECISE
        dtgrid_nbq = dtfast
# else
        dtgrid_nbq = dt
# endif

!.........EOS to compute rho (if not already done):
# ifdef NONLIN_EOS
        call rho_eos_tile(Istr,Iend,Jstr,Jend,K_up,K_dw)
# else
        call rho_eos_tile(Istr,Iend,Jstr,Jend)
# endif

# ifdef NBQ_MASS
!.........Initialize NBQ density field:
        do k=1,N
          do j=JR_RANGE
            do i=IR_RANGE
              rho_nbq_avg1(i,j,k)=(rho0+rho(i,j,k))/rho0
            enddo
          enddo
        enddo 
        do j=JR_RANGE
          do i=IR_RANGE
            work2d(i,j)      =0.
            rhobar_nbq(i,j,:)=0.
          enddo
        enddo
        do k=1,N
          do j=JR_RANGE
            do i=IR_RANGE
              work2d(i,j)      =work2d(i,j)+Hzr(i,j,k)
              rhobar_nbq(i,j,:)=rhobar_nbq(i,j,:)+
     &                          rho(i,j,k)*Hzr(i,j,k)/rho0
            enddo
          enddo
        enddo
 
!....... add rho0 for added precision:
          do j=JR_RANGE
            do i=IR_RANGE
            rhobar_nbq(i,j,:)   =rhobar_nbq(i,j,:)/work2d(i,j) + 1.
            rhobar_nbq_avg1(i,j)=rhobar_nbq(i,j,1) 
          enddo
        enddo
# else
!      rhobar_nbq     (:,:,:)=1.
!      rhobar_nbq_avg1(:,:  )=1.
# endif
        
!.......Remaining initializations:
        do k=1,N
          do j=JR_RANGE
            do i=IR_RANGE
              thetadiv_nbq(i,j,k)=0.
            enddo
          enddo
        enddo

# ifdef NBQ_MASS
#  if defined EW_PERIODIC || defined NS_PERIODIC || defined  MPI
        call exchange_r2d_tile (Istr,Iend,Jstr,Jend
     &                         ,rho_nbq_avg1(START_2D_ARRAY,1))
        call exchange_r2d_tile (Istr,Iend,Jstr,Jend
     &                         ,rhobar_nbq(START_2D_ARRAY,1))
        call exchange_r2d_tile (Istr,Iend,Jstr,Jend
     &                         ,rhobar_nbq(START_2D_ARRAY,2))
        call exchange_r2d_tile (Istr,Iend,Jstr,Jend
     &                         ,rhobar_nbq(START_2D_ARRAY,3))
        call exchange_r2d_tile (Istr,Iend,Jstr,Jend
     &                         ,rhobar_nbq(START_2D_ARRAY,4))
#  endif
# endif

      endif 

      if (icall == 3) then 
!**********************************************************************
!
!               NBQ initializations (PART III)
!
!**********************************************************************
        do k=1,N
          do j=JR_RANGE
            do i=IU_RANGE
              qdmu_nbq(i,j,k)=0.
            enddo
          enddo
        enddo
# ifdef NBQ_ZETAW
        do j=JR_RANGE
          do i=IU_RANGE
            usurf_nbq=0.
          enddo
        enddo
# endif
        do k=1,N
          do j=JV_RANGE
            do i=IR_RANGE
              qdmv_nbq(i,j,k)=0.
            enddo
          enddo
        enddo
# ifdef NBQ_ZETAW
        do j=JV_RANGE
          do i=IR_RANGE
            vsurf_nbq=0.
          enddo
        enddo
# endif
        do k=1,N
          do j=JR_RANGE
            do i=IR_RANGE
              qdmw_nbq(i,j,k)=0.
              rho_nbq(i,j,k) =0.
# ifndef NBQ_ZETAW
              zr_half_nbq(i,j,k)=0.
              zw_half_nbq(i,j,k)=0.
# endif
            enddo
          enddo
        enddo
# ifdef NBQ_ZETAW
        do j=JR_RANGE
          do i=IR_RANGE
            wsurf_nbq=0.
          enddo
        enddo
# endif

# ifdef KH_INST
        if (iic.le.1.and.iif.le.1) then
          do k=1,N
            do j=JR_RANGE
              do i=IU_RANGE
                qdmu_nbq(i,j,k)=0.5*u(i,j,k,nrhs)
# ifdef NBQ_MASS
     &                          *(Hz(i,j,k)+Hz(i-1,j,k))
     &                          *(1.+0.5*(rho(i,j,k)+rho(i-1,j,k))/rho0)
# else
     &                          *(Hz(i,j,k)+Hz(i-1,j,k))
# endif
              enddo
            enddo
          enddo

          do k=1,N-1
            do j=JR_RANGE
              do i=IR_RANGE
                qdmw_nbq(i,j,k)=0.5*wz(i,j,k,nrhs)
# ifdef NBQ_MASS
     &                          *(Hz(i,j,k)+Hz(i,j,k+1))
     &                          *(1.+0.5*(rho(i,j,k)+rho(i,j,k+1))/rho0)
# else
     &                          *(Hz(i,j,k)+Hz(i,j,k+1))
# endif
              enddo
            enddo
          enddo
          k=0 
          do j=JR_RANGE 
            do i=IR_RANGE
                qdmw_nbq(i,j,k)=0.5*wz(i,j,k,nrhs)
# ifdef NBQ_MASS
     &                          *(Hz(i,j,k)+Hz(i,j,k+1))
     &                          *(1.+0.5*(rho(i,j,k)+rho(i,j,k+1))/rho0)
# else
     &                          *(Hz(i,j,k)+Hz(i,j,k+1))
# endif            enddo
          enddo
          k=N 
          do j=JR_RANGE
            do i=IR_RANGE
                qdmw_nbq(i,j,k)=0.5*wz(i,j,k,nrhs)
# ifdef NBQ_MASS
     &                          *(Hz(i,j,k)+Hz(i,j,k+1))
     &                          *(1.+0.5*(rho(i,j,k)+rho(i,j,k+1))/rho0)
# else
     &                          *(Hz(i,j,k)+Hz(i,j,k+1))
# endif
            enddo
          enddo
        endif ! iic
# endif

      endif  ! icall
      
      return
      end subroutine initial_nbq_tile
#else
      subroutine initial_nbq_empty
      return
      end
#endif
