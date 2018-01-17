#include "cppdefs.h"
#ifdef NBQ

      subroutine density_nbq(icall)

!======================================================================
!
!                      Various Computations related to
!                            NBQ density
!
!> @note Main web documentation at:
! http://poc.obs-mip.fr/auclair/WOcean.fr/SNH/index_snh_home.htm
!
! DESCRIPTION: 
!
!> @brief Density anomaly driver 
!
!> @details 
! - icall=0 vertical averaging density to get the external 
!           mode density (rhpbar_nbq_t).
! - icall=1 internal mode density (rhpio_nbq_t).
! - icall=5 computes the NH pressure gradient for the internal mode.
! - icall=6 computes the divergence of the momentum (div_nbq_a) 
!      + used in the momentum equation (second viscosity term with 
!        gradient operator)
!      + used in the continuity equation
! - icall=7 time incrementation of NBQ mode density anomaly.
!
! REVISION HISTORY:
!
!> @authors
!> @date 2015 January
!> @todo
!
!======================================================================
!
!      use module_principal , only : kount0, iteration3d, rhp_t, rho, mask_t, &
!			    dz_t, iteration2d_max_now, hz_w, iteration2d,    &
!			    imax, jmax, kmax
!     use module_parallele                ! #MPI#
      use module_nh                       ! #NH#
      use module_nbq                      ! #NBQ#
      implicit none
      integer :: icall, i,j ,k, k1,indm_d

# include "param_F90.h"
# include "scalars_F90.h"
# include "work.h"
# include "grid.h"
# include "ocean2d.h"
# include "ocean3d.h"
# include "nbq.h"

      integer :: ncp
      real    :: dist_d

      double precision :: t1_d,t2_d

      if (icall.eq.20) then
#ifndef NBQ_VOL
!
!**********************************************************************
!  Transfer density field to i,j,k array 
!  and time filter, ready for external mode
!**********************************************************************

#ifdef old_version 
       
        rhobar_nbq(istrq_nh-1:iendq_nh+1,jstrq_nh-1:jendq_nh+1,knew)=0. 
        work2d    (istrq_nh-1:iendq_nh+1,jstrq_nh-1:jendq_nh+1)=0.

        do j=jstrq_nh-1,jendq_nh+1
        do k=1,N
        do i=istrq_nh-1,iendq_nh+1     
           l_nbq=ijk2lq_nh(i,j,k)
           rho_nbq_ext(i,j,k)  = 1.D0 + (rhp_nbq_a(l_nbq)+rho(i,j,k))/rho0     
           work2d(i,j)         = work2d(i,j)+Hzr_half_nbq(i,j,k)
           rhobar_nbq(i,j,knew)= 1./float(N)+ rhobar_nbq(i,j,knew) &
               +(rhp_nbq_a(l_nbq)+rho(i,j,k))/rho0*Hzr_half_nbq(i,j,k) 
         enddo
         enddo
         enddo
!
!.......Rho0 added subsequently for added precision 
!
         do j=jstrq_nh-1,jendq_nh+1
         do i=istrq_nh-1,iendq_nh+1
           rhobar_nbq(i,j,knew) = 1.D0 + rhobar_nbq(i,j,knew)/max(work2d(i,j),1.d-30)
         enddo
         enddo
#endif

        rhobar_nbq(istrq_nh-1:iendq_nh+1,jstrq_nh-1:jendq_nh+1,knew)=0.
!       work2d    (istrq_nh-1:iendq_nh+1,jstrq_nh-1:jendq_nh+1)=0.

         do l_nbq = 1 , neqcont_nh
           i     = l2iq_nh (l_nbq)
           j     = l2jq_nh (l_nbq)
           k     = l2kq_nh (l_nbq)                    
!          work2d(i,j)         = work2d(i,j)+Hzr_half_nbq(i,j,k)
           rhobar_nbq(i,j,knew)= rhobar_nbq(i,j,knew)                           &
                     +(rhp_nbq_a(l_nbq)+rho(i,j,k))/rho0 *Hzr_half_nbq(i,j,k)
           rho_nbq_ext(i,j,k)  = 1.+(rhp_nbq_a(l_nbq)+rho(i,j,k)  )/rho0     
         enddo
!
!.......Rho0 added subsequently for added precision 
!
!        rho_nbq_ext(:,:,:) = 1.D0 + rho_nbq_ext(:,:,:) 
         do j=jstrq_nh-1,jendq_nh+1
         do i=istrq_nh-1,iendq_nh+1
           rhobar_nbq(i,j,knew) = 1.D0 + rhobar_nbq(i,j,knew) & 
              / (zw_half_nbq(i,j,N)-zw_half_nbq(i,j,0))
         enddo
         enddo

#if defined EW_PERIODIC || defined NS_PERIODIC || defined  MPI
!        call exchange_r2d_tile (Istr_nh,Iend_nh,Jstr_nh,Jend_nh,rhobar_nbq(START_2D_ARRAY,knew))
!        call exchange_r3d_tile (Istr_nh,Iend_nh,Jstr_nh,Jend_nh &
!                               ,rho_nbq_ext(START_2D_ARRAY,1))
#endif

        

   !    rhobar_nbq=1.
   !    rho_nbq_ext=1.

#ifdef RVTK_DEBUG
      call check_tab3d(rho_nbq_ext(:,:,1:N),'rho_nbq_ext (density_nbq)','r')
      call check_tab2d(rhobar_nbq(:,:,knew),'rhobar_nbq (density_nbq)','r')
#endif    



#endif

      elseif (icall.eq.6) then
!
!**********************************************************************
!     Calcul de la divergence
!**********************************************************************
!
            call amux(                                                &
                  neqcont_nh                                          &
                 ,qdm_nbq_a(1)                              &
!                  ,div_nbq_a(1)                              &
                 ,div_nbq_a(1)                                        &
                 ,contv_nh (1)                                        &
                 ,contj_nh (1)                                        &
                 ,conti_nh (1)                                        &
                       )

      elseif (icall.eq.60) then
!
!**********************************************************************
!     Calcul de la divergence
!**********************************************************************
!
#ifndef NBQ_IMP
            call amux(                                                &
                  neqcont_nh                                          &
                 ,qdm_nbq_a(1:neqmom_nh(1)+neqmom_nh(2))             &
!                ,div_nbq_a(1,dnrhs_nbq)                              &
                 ,div_nbq_a(1)		                              &
                 ,contv_nh (1)                                        &
                 ,contj_nh (1)                                        &
                 ,conti_nh (1)                                        &
                       )         
#endif
 
            call amux(                                                &
                  neqcont_nh                                          &
                 ,qdm_nbq_a(1)     &
!                ,div_nbq_a(1,dnrhs_nbq)                              &
                 ,divz_nbq_a(1)		                              &
                 ,contzv_nh (1)                                        &
                 ,contzj_nh (1)                                        &
                 ,contzi_nh (1)                                        &
                       )

      div_nbq_a(1:)=div_nbq_a(1:)+divz_nbq_a(1:)

      elseif (icall.eq.61) then
!
!**********************************************************************
!     Calcul de la divergence
!**********************************************************************
!
            call amux(                                                &
                  neqcont_nh                                          &
                 ,qdm_nbq_a(1)                              &
!                  ,div_nbq_a(1)                             &
                 ,div_nbq_a(1)                                        &
                 ,contv_nh (1)                                        &
                 ,contj_nh (1)                                        &
                 ,conti_nh (1)                                        &
                       )

      elseif (icall.eq.62) then
!
!**********************************************************************
!     Calcul de la divergence
!**********************************************************************
!
            call amux(                                                &
                  neqcont_nh                                          &
                 ,qdm_nbq_a(1)                              &
!                  ,div_nbq_a(1)                              &
                 ,div_nbq_a(1)		                              &
                 ,contv_nh (1)                                        &
                 ,contj_nh (1)                                        &
                 ,conti_nh (1)                                        &
                       )

      elseif (icall.eq.7) then
        return
         
!
!*******************************************************************
!......Move forward: Masse
!*******************************************************************
!          
!   Remove Leap-Frog in step_NBQ Only FW-BW (FB)  scheme
!          ncp       = rnnew_nbq
!          rnnew_nbq = rnstp_nbq
!          rnstp_nbq = rnrhs_nbq
!          rnrhs_nbq = ncp
         ncp       = rnnew_nbq
	 rnnew_nbq = rnrhs_nbq 
	 rnrhs_nbq = ncp 

!          ncp       = dnnew_nbq
!          dnnew_nbq = dnstp_nbq
!          dnstp_nbq = dnrhs_nbq
!          dnrhs_nbq = ncp

# ifdef ACOUSTIC
      elseif (icall.eq.10) then
!
!*******************************************************************
!......Acoustic waves: Initialization
!*******************************************************************
          period_exp = 0.025/2.
          for_a_exp  = 2.5
          amp_exp = 1.e-3
          hmax_exp = 128.
          dg_exp = 2.

      elseif (icall.eq.11) then
!
!*******************************************************************
!......Acoustic waves: Forcing
!*******************************************************************
!
        time_nbq = time_nbq + 0.5*dtnbq

        do l_nbq = 1 , neqcont_nh 

          i=l2iq_nh(l_nbq)
          j=l2jq_nh(l_nbq)
          k=l2kq_nh(l_nbq)

          dist_d=sqrt((xr(i,j)-xl/2.)**2+(0.*(yr(i,j)-el/2.))**2       &
                              +(abs(z_r(i,j,k))-hmax_exp/2.)**2)
!         if (dist_d.le.for_a_exp(1)) then
!              div_nbq_a(l_nbq,dnrhs_nbq) = div_nbq_a(l_nbq,dnrhs_nbq)   &
             div_nbq_a(l_nbq) = div_nbq_a(l_nbq)   &
                        +amp_exp*sin(2*acos(-1.)*time_nbq/period_exp)  &
                                        *exp(-dist_d**2/for_a_exp**2)
!         endif
        enddo
!         write(6,*) 'ACOUSTIC',div_nbq_a(10,dnrhs_nbq)
# endif /* ACOUSTIC */

      endif  ! icall

        return
        end
#else
        subroutine density_nbq_empty
        return
        end
#endif
