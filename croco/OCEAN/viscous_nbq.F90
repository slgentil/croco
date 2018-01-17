#include "cppdefs.h"
#ifdef NBQ

      subroutine viscous_nbq(icall)
!**********************************************************************
!
!                      Various Computations related to
!                              2nd viscosity
!
!> @note <a href="http://poc.obs-mip.fr/auclair/WOcean.fr/SNH/index_snh_home.htm"> Main web documentation </a>
!
! DESCRIPTION: 
!
!> @brief <a href="http://poc.obs-mip.fr/auclair/WOcean.fr/SNH/Restricted/NH-NBQ/Html_pages/Schemas_Filtrage.htm">
!! Filtering pseudo-acoustic waves</a>.
!! Second viscosity calculation. 
!> @details Divergence operator applied on the momentum, not on velocity.
!! The NH operators are reused.
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
# include "ocean2d.h"
# include "ocean3d.h"
# include "grid.h"
# include "nbq.h"

      integer :: i,j,k,icall,indm_v

      if (icall.eq.0) then
!*******************************************************************
!*******************************************************************
!       Initialisation de la seconde viscosite
!       calcul de la viscosite de reference...
!
!*******************************************************************
!*******************************************************************
       
        visc2_nbq_a(1:neqmom_nh(0)) = visc2_nbq
        csvisc1_nbq  = dtnbq * soundspeed2_nbq + visc2_nbq
        csvisc2_nbq  = dtnbq * soundspeed2_nbq / csvisc1_nbq 
       

      endif

      if (icall.eq.1) then
!*******************************************************************
!*******************************************************************
!               Seconde viscosite:
!
! Attention: telle que codee actuellement la seconde viscosite
!            est appliquee sur la qdm et non la vitesse comme cela
!            devrait etre le cas pour la version moleculaire !
!            (normalisation par rho dans main_nbq)
!           On "gagne" ainsi un produit MAT*VECT dans la boucle NBQ
!           On conserve toutefois la divergence calculee sur l'ancienne
!           grille... soit probleme potentiel!
!
!*******************************************************************
!*******************************************************************

         indm_v=momi_nh(neqcont_nh+1)+1
         rhsd2_nbq=0.
         call amux(                                                 &
              neqmom_nh(0)                                          &
!              ,div_nbq_a(1:neqcont_nh,dnstp_nbq)                     &    ! div decentree (0)
             ,div_nbq_a(1:neqcont_nh)                     &    ! div decentree (0)
             ,rhsd2_nbq(1:neqmom_nh(0))                             &
             ,momv_nh(1:indm_v)                                     &
             ,momj_nh(1:indm_v)                                     &
             ,momi_nh(1:neqcont_nh+1)                               &
                 )   
    
      endif

      if (icall.eq.2) then
!*******************************************************************
!*******************************************************************
!               Seconde viscosite:
!
! Attention: telle que codee actuellement la seconde viscosite
!            est appliquee sur la qdm et non la vitesse comme cela
!            devrait etre le cas pour la version moleculaire !
!            (normalisation par rho dans main_nbq)
!           On "gagne" ainsi un produit MAT*VECT dans la boucle NBQ
!           On conserve toutefois la divergence calculee sur l'ancienne
!           grille... soit probleme potentiel!
!
!*******************************************************************
!*******************************************************************

         indm_v=momi_nh(neqcont_nh+1)+1
!         rhsd2_nbq=0.D0

#ifndef NBQ_IMP
         call amux(                                                 &
              neqmom_nh(0)                                          &
!            ,div_nbq_a(1:neqcont_nh,dnrhs_nbq)                     &    ! div decentree (0)
             ,div_nbq_a(1:neqcont_nh)                     &    ! div decentree (0)
             ,rhsd2_nbq(1:neqmom_nh(0))                             &
             ,momv_nh(1:indm_v)                                     &
             ,momj_nh(1:indm_v)                                     &
             ,momi_nh(1:neqcont_nh+1)                               &
                 )  
#else
         call amux(                                                 &
              neqmom_nh(1)+neqmom_nh(2)                             &
!            ,div_nbq_a(1:neqcont_nh,dnrhs_nbq)                     &    ! div decentree (0)
             ,div_nbq_a(1:neqcont_nh)                     &    ! div decentree (0)
             ,rhsd2_nbq(1:neqmom_nh(0))                             &
             ,momv_nh(1:indm_v)                                     &
             ,momj_nh(1:indm_v)                                     &
             ,momi_nh(1:neqcont_nh+1)                               &
                 )  
#endif
      endif

      if (icall.eq.12) then
!*******************************************************************
!*******************************************************************
!               Seconde viscosite:
!
! Attention: telle que codee actuellement la seconde viscosite
!            est appliquee sur la qdm et non la vitesse comme cela
!            devrait etre le cas pour la version moleculaire !
!            (normalisation par rho dans main_nbq)
!           On "gagne" ainsi un produit MAT*VECT dans la boucle NBQ
!           On conserve toutefois la divergence calculee sur l'ancienne
!           grille... soit probleme potentiel!
!
!*******************************************************************
!*******************************************************************
         return
         indm_v=momi_nh(neqcont_nh+1)
         rhsd2_nbq=0.
         call amux(                                                 &
              neqmom_nh(0)                                          &
!              ,div_nbq_a(1:neqcont_nh,dnrhs_nbq)                     &    ! div decentree (0)
             ,div_nbq_a(1:neqcont_nh)                     &    ! div decentree (0)
             ,rhsd2_nbq(1:neqmom_nh(0))                             &
             ,momv_nh(1:indm_v)                                     &
             ,momj_nh(1:indm_v)                                     &
             ,momi_nh(1:neqcont_nh+1)                               &
                 )   
   !    rhsd2_nbq = 0.
    
      endif

      end subroutine viscous_nbq

#else
      subroutine viscous_nbq_empty
      end subroutine viscous_nbq_empty
#endif
