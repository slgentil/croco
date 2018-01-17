#include "cppdefs.h"
#ifdef NBQ

      subroutine ru_nbq(icall)

!**********************************************************************
!
!                      Various Computations related to
!                            NBQ momentum
!
!> @note <a href="http://poc.obs-mip.fr/auclair/WOcean.fr/SNH/index_snh_home.htm"> Main web documentation </a>
!
! DESCRIPTION: 
!
!> @brief NBQ momentum routine.
!
! REVISION HISTORY:
!
!> @authors
!> @date 2015 January
!> @todo
!
!**********************************************************************

      use module_nh 
      use module_nbq
      implicit none

# include "param_F90.h"
# include "scalars_F90.h"
# include "grid.h"
# include "ocean3d.h"
# include "nbq.h"

      integer :: i,j,k
      integer :: icall,ncp
      real    :: cff


      if (icall.eq.1) then
!
!*******************************************************************
!  Fill external and internal forcing terms
!*******************************************************************
!
#if defined EW_PERIODIC || defined NS_PERIODIC || defined  MPI
      call exchange_u3d_tile (Istru_nh,Iendu_nh,Jstru_nh,Jendu_nh,  &
                                       ruint_nbq(START_2D_ARRAY,1))
      call exchange_u3d_tile (Istru_nh,Iendu_nh,Jstru_nh,Jendu_nh,  &
                        ruext_nbq(START_2D_ARRAY,1))
      call exchange_u3d_tile (Istrv_nh,Iendv_nh,Jstrv_nh,Jendv_nh,  &
                                       rvint_nbq(START_2D_ARRAY,1))
      call exchange_u3d_tile (Istrv_nh,Iendv_nh,Jstrv_nh,Jendv_nh,  &
                        rvext_nbq(START_2D_ARRAY,1))
      call exchange_u3d_tile (Istr_nh,Iend_nh,Jstr_nh,Jend_nh,  &
                                       rwint_nbq(START_2D_ARRAY,0))

#endif
!!       do l_nbq = nequ_nh(2)+1,nequ_nh(5)
         do l_nbq = 1,nequ_nh(7)
            i=l2imom_nh(l_nbq)
            j=l2jmom_nh(l_nbq)
            k=l2kmom_nh(l_nbq)
            dqdmdt_nbq_a(l_nbq)=rho0*(ruint_nbq(i,j,k)+ruext_nbq(i,j,k))
         enddo

!       do l_nbq = neqv_nh(2)+1,neqv_nh(5)  
        do l_nbq = nequ_nh(7)+1,neqv_nh(7)  
           i=l2imom_nh(l_nbq)
           j=l2jmom_nh(l_nbq)
           k=l2kmom_nh(l_nbq)
           dqdmdt_nbq_a(l_nbq)=rho0*(rvint_nbq(i,j,k)+rvext_nbq(i,j,k))
        enddo

        if (iif==1) then
!        do l_nbq = neqw_nh(2)+1,neqw_nh(5)
         do l_nbq = neqv_nh(7)+1,neqw_nh(7)
           i=l2imom_nh(l_nbq)
           j=l2jmom_nh(l_nbq)
           k=l2kmom_nh(l_nbq)
           dqdmdt_nbq_a(l_nbq)=rho0*rwint_nbq(i,j,k)
         enddo
        endif
       
      elseif (icall.eq.2) then
!
!*******************************************************************
!  Prepare feedback of NBQ rhs terms to external equations
!*******************************************************************
!
        cff=1/(rho0*real(ndtnbq))
!        
! X-direction:
!
        rubar_nbq(:,:)=0.
!       do l_nbq = nequ_nh(2)+1,nequ_nh(5)
!        do l_nbq = nequ_nh(1)+1,nequ_nh(6)
        do l_nbq = 1,nequ_nh(7)
           i=l2imom_nh(l_nbq)
           j=l2jmom_nh(l_nbq)
           k=l2kmom_nh(l_nbq)
           ru_nbq_ext(i,j,k)   = cff*rhssum_nbq_a(l_nbq)*on_u(i,j)*om_u(i,j)
           rhssum_nbq_a(l_nbq) = 0.
           rubar_nbq(i,j)      = rubar_nbq(i,j)+ru_nbq_ext(i,j,k)
        enddo
        
#if defined EW_PERIODIC || defined NS_PERIODIC || defined  MPI
      call exchange_u3d_tile (Istru_nh,Iendu_nh,Jstru_nh,Jendu_nh,  &
                                       ru_nbq_ext(START_2D_ARRAY,1))

      call exchange_u2d_tile (Istru_nh,Iendu_nh,Jstru_nh,Jendu_nh,  &
                        rubar_nbq(START_2D_ARRAY))
#endif
 !     rubar_nbq=0.
 !     ru_nbq_ext=0.

#ifdef RVTK_DEBUG
       call check_tab3d(ru_nbq_ext(:,:,1:N),'ru_nbq_ext (ru_nbq)','u')
       call check_tab2d(rubar_nbq(:,:),'rubar_nbq (ru_nbq)','u')
#endif    
!
! Y-direction:
!
        rvbar_nbq(:,:)=0.
!       do l_nbq = neqv_nh(2)+1,neqv_nh(5)  
!       do l_nbq = neqv_nh(1)+1,neqv_nh(6)  
        do l_nbq = nequ_nh(7)+1,neqv_nh(7)  
            i=l2imom_nh(l_nbq)
            j=l2jmom_nh(l_nbq)
            k=l2kmom_nh(l_nbq)
            rv_nbq_ext(i,j,k)   = cff*rhssum_nbq_a(l_nbq)*on_v(i,j)*om_v(i,j)
            rhssum_nbq_a(l_nbq) = 0.
            rvbar_nbq(i,j)      = rvbar_nbq(i,j)+rv_nbq_ext(i,j,k)
        enddo

#if defined EW_PERIODIC || defined NS_PERIODIC || defined  MPI
      call exchange_v3d_tile (Istrv_nh,Iendv_nh,Jstrv_nh,Jendv_nh,  &    
                                       rv_nbq_ext(START_2D_ARRAY,1)) 

      call exchange_v2d_tile (Istrv_nh,Iendv_nh,Jstrv_nh,Jendv_nh,  &
                        rvbar_nbq(START_2D_ARRAY))
#endif

#ifdef RVTK_DEBUG
       call check_tab3d(rv_nbq_ext(:,:,1:N),'rv_nbq_ext (ru_nbq)','v')
       call check_tab2d(rvbar_nbq(:,:),'rvbar_nbq (ru_nbq)','v')
#endif    

!    
! Z-direction:
!
!       do l_nbq = neqw_nh(2)+1,neqw_nh(5)
!       do l_nbq = neqw_nh(1)+1,neqw_nh(6)
        do l_nbq = neqv_nh(7)+1,neqw_nh(7)
            i = l2imom_nh (l_nbq)
            j = l2jmom_nh (l_nbq)
            k = l2kmom_nh (l_nbq)
            rw_nbq_ext(i,j,k)   = cff*rhssum_nbq_a(l_nbq)*on_r(i,j)*om_r(i,j)
            rhssum_nbq_a(l_nbq) = 0.
        enddo

#if defined EW_PERIODIC || defined NS_PERIODIC || defined  MPI
      call exchange_r3d_tile (Istr_nh,Iend_nh,Jstr_nh,Jend_nh,  &    
                                       rw_nbq_ext(START_2D_ARRAY,0)) 

#endif
#ifdef RVTK_DEBUG
       call check_tab3d(rw_nbq_ext(:,:,0:N),'rw_nbq_ext (ru_nbq)','v')
#endif    


      elseif (icall.eq.6) then
!
!*******************************************************************
!  Increment momentum :
!*******************************************************************
!
#ifndef NBQ_IMP
          div_nbq_a(1:neqcont_nh)=-visc2_nbq*div_nbq_a(1:neqcont_nh) /dtnbq     &
                                  +soundspeed2_nbq*rhp_nbq_a(1:neqcont_nh)
         call amux(                                                        &
                neqcorrt_nbq                                               &
               ,div_nbq_a(1:neqcont_nh)                                    & 
!              ,rhp_nbq_a(1:neqcont_nh)-rhp_bq_a(1:neqcont_nh)             & 
               ,rhs1_nbq (1)                                               &
               ,momv_nh (1)                                                &
               ,momj_nh  (1)                                               & 
               ,momi_nh  (1)                                               &
                 )
#else
          div_nbq_a(1:neqcont_nh)=-visc2_nbq*div_nbq_a(1:neqcont_nh)/dtnbq       &
                                  +soundspeed2_nbq*rhp_nbq_a(1:neqcont_nh)
         call amux(                                                        &
                neqmom_nh(1)+neqmom_nh(2)                                  &
               ,div_nbq_a(1:neqcont_nh)                                    & 
!              ,rhp_nbq_a(1:neqcont_nh)-rhp_bq_a(1:neqcont_nh)             & 
               ,rhs1_nbq (1)                                               &
               ,momv_nh (1)                                                &
               ,momj_nh  (1)                                               & 
               ,momi_nh  (1)                   &
                 )
#endif

      elseif (icall.eq.7) then
         return
       
!
!*******************************************************************
!  Move forward momentum
!*******************************************************************
!
!   Remove Leap-Frog in step_NBQ Only FW-BW (FB)  scheme
!           ncp       = vnnew_nbq
!           vnnew_nbq = vnstp_nbq
!           vnstp_nbq = vnrhs_nbq
!           vnrhs_nbq = ncp
	    ncp       = vnnew_nbq
	    vnnew_nbq = vnrhs_nbq
	    vnrhs_nbq = ncp
	    

!           elseif (icall.eq.8) then
! !
! !*******************************************************************
! !  Move forward momentum
! !*******************************************************************
! !
! 
!           ncp = vnrhs_nbq
!           vnstp_nbq = vnnew_nbq
!           vnrhs_nbq = vnstp_nbq
!           vnnew_nbq = ncp
          
       endif  ! icall

       return
       end
#else
      subroutine ru_nbq_empty 
      return
      end
#endif

