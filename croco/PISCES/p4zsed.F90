#include "cppdefs.h"

MODULE p4zsed
   !!======================================================================
   !!                         ***  MODULE p4sed  ***
   !! TOP :   PISCES Compute loss of organic matter in the sediments
   !!======================================================================
   !! History :   1.0  !  2004-03 (O. Aumont) Original code
   !!             2.0  !  2007-12  (C. Ethe, G. Madec)  F90
   !!----------------------------------------------------------------------
#if defined key_pisces
   !!----------------------------------------------------------------------
   !!   'key_pisces'                                       PISCES bio-model
   !!----------------------------------------------------------------------
   !!   p4z_sed        :  Compute loss of organic matter in the sediments
   !!   p4z_sbc        :  Read and interpolate time-varying nutrients fluxes
   !!   p4z_sed_init   :  Initialization of p4z_sed
   !!----------------------------------------------------------------------
   USE sms_pisces
   USE p4zint
   USE p4zopt
   USE p4zsink
   USE p4zrem
   USE p4zlim


   IMPLICIT NONE
   PRIVATE

   PUBLIC   p4z_sed   
   PUBLIC   p4z_sed_alloc   
   PUBLIC   p4z_sed_init   
   PUBLIC   p4z_sed_nam   

   !!* Substitution
#  include "ocean2pisces.h90"
#  include "top_substitute.h90"

   !! * Shared module variables
   LOGICAL, PUBLIC ::    &
     ln_dust     = .FALSE.      ,  &  !:
     ln_river    = .FALSE.      ,  &  !:
     ln_ndepo    = .FALSE.      ,  &  !:
     ln_sedinput = .FALSE.      ,  &  !:
     ln_sedmeta  = .FALSE.            !:

   REAL(wp), PUBLIC ::   &
     sedfeinput = 1.E-9   ,  &  !:
     dustsolub  = 0.014   ,  &  !:
     mfrac      = 0.035   ,  &  ! Fe mineral fraction of dust
     wdust      =  2.0    ,  &  ! Dust sinking speed 
     concfediaz = 1.E-10  ,  &  !: Diazotrophs half-saturation Cste for Iron
     diazolight = 50      ,  &  !: Diazotrophs sensitivity to light (W/m2)
     nitrfix    = 1.E-7         !:

   !! * Module variables
!   INTEGER ::                   &
!     ryyss,                     &  !: number of seconds per year
!     rmtss                         !: number of seconds per month

   REAL(wp), PUBLIC   ::  year2daydta                 
   INTEGER ::                   &
      numdust,                  &  !: logical unit for surface fluxes data
      nflx1 , nflx2,            &  !: first and second record used
      nflx11, nflx12      ! ???
   REAL(wp), DIMENSION(:,:,:), PUBLIC, ALLOCATABLE, SAVE ::    &  !:
     dustmo, no3depmo, nh4depmo, ironsed                                !: 2 consecutive set of dust fields 
   REAL(wp), DIMENSION(:,:,:), PUBLIC, ALLOCATABLE, SAVE ::    &  !:
     cmask
   REAL(wp), DIMENSION(:,:), ALLOCATABLE, SAVE   ::    &
     rivinp, cotdep, dust,  &
     no3dep, nh4dep, po4dep, sidep
   REAL(wp), DIMENSION(:,:,:), ALLOCATABLE, SAVE  ::   &
      nitrpot, irondep
   REAL(wp) :: sumdepsi, rivalkinput, rivpo4input, nitdepinput


CONTAINS

   SUBROUTINE p4z_sed(kt, jnt)
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE p4z_sed  ***
      !!
      !! ** Purpose :   Compute loss of organic matter in the sediments. This
      !!              is by no way a sediment model. The loss is simply 
      !!              computed to balance the inout from rivers and dust
      !!
      !! ** Method  : - ???
      !!---------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kt, jnt ! ocean time step
      INTEGER  ::   ji, jj, jk
!      INTEGER  ::   ikt
#if ! defined key_sed
      REAL(wp) ::   zsumsedsi, zsumsedpo4, zsumsedcal
      REAL(wp) ::   zrivsil, zrivalk, zrivno3
#endif
      REAL(wp) :: zconctmp, zws3, zws4, zwsc 
      REAL(wp) :: zlim, zconctmp2, zfact, zmsk
      REAL(wp) :: zrfact2
      REAL(wp) :: zo2, zno3, zflx, zpdenit, z1pdenit, zdenitt, zolimit
      REAL(wp) :: zsiloss, zcaloss, zfactcal, zdep, zwstpoc
      REAL(wp) :: ztrfer, ztrpo4, zwdust, zlight
      REAL(wp), DIMENSION(PRIV_2D_BIOARRAY) :: zdenit2d, zbureff
      REAL(wp), DIMENSION(PRIV_2D_BIOARRAY) :: zno3dep,znh4dep
      CHARACTER (len=25) :: charout
      !!---------------------------------------------------------------------

      IF( (jnt == 1) .and. ( ln_dust .OR. ln_ndepo ) )  CALL p4z_sbc( kt )

      DO jj = JRANGE
         DO ji = IRANGE
            zno3dep(ji,jj) = 0.
            znh4dep(ji,jj) = 0.
         END DO
      END DO
      ! Iron and Si deposition at the surface
      ! -------------------------------------
      IF( ln_dust ) THEN
         DO jj = JRANGE
            DO ji = IRANGE
               irondep(ji,jj,1) = ( dustsolub * mfrac * dust(ji,jj) &
                  &               / ( 55.85 * rmtss ) + 3.e-10 / ryyss )   &
                  &               * rfact2 / fse3t(ji,jj,KSURF)
                  !
               sidep  (ji,jj)   = 8.8 * 0.075  &
                  &              * mfrac * dust(ji,jj) / 28.1 / rmtss    &
                  &              * rfact2 / fse3t(ji,jj,KSURF) 
                  !
               po4dep (ji,jj)   = 0.1 * 0.021  &
                  &              * mfrac * dust(ji,jj) / 31. / rmtss    &
                  &              / po4r                                 &
                  &              * rfact2 / fse3t(ji,jj,KSURF) 
                  !
               tra(ji,jj,1,jpsil) = tra(ji,jj,1,jpsil) + sidep (ji,jj)
               tra(ji,jj,1,jppo4) = tra(ji,jj,1,jppo4) + po4dep(ji,jj)
               tra(ji,jj,1,jpfer) = tra(ji,jj,1,jpfer) + irondep(ji,jj,1)
            END DO
          END DO

          ! Iron solubilization of particles in the water column
          ! ----------------------------------------------------
          ! dust in kg/m2/s ---> 1/55.85 to put in mol/Fe
          !  wdust in m/j
          zwdust = 0.03 * rday / ( wdust * 55.85 ) / ( 270. * rday )
          DO jk = KRANGEL
             DO jj = JRANGE
                DO ji = IRANGE
                   irondep(ji,jj,jk) = dust(ji,jj) * mfrac * zwdust * rfact2 &
                   &                  * EXP( -fsdept(ji,jj,K) / 540. )
                   tra(ji,jj,jk,jpfer) = tra(ji,jj,jk,jpfer) + irondep(ji,jj,jk)
                END DO
             END DO
         END DO
      ENDIF

      ! Add the external input of nutrients, carbon and alkalinity
      ! ----------------------------------------------------------
      IF( ln_river ) THEN
         DO jj = JRANGE
            DO ji = IRANGE
               tra(ji,jj,1,jppo4) = tra(ji,jj,1,jppo4) + rivinp(ji,jj) * rfact2 
               tra(ji,jj,1,jpno3) = tra(ji,jj,1,jpno3) + rivinp(ji,jj) * rfact2
               tra(ji,jj,1,jpfer) = tra(ji,jj,1,jpfer) + rivinp(ji,jj) * 3.e-5 * rfact2
               tra(ji,jj,1,jpsil) = tra(ji,jj,1,jpsil) + cotdep(ji,jj)   * rfact2 / 6.
               tra(ji,jj,1,jpdic) = tra(ji,jj,1,jpdic) + rivinp(ji,jj) * 2.631 * rfact2
               tra(ji,jj,1,jptal) = tra(ji,jj,1,jptal) &
                  &                + ( cotdep(ji,jj) - rno3 * rivinp(ji,jj) )  * rfact2
            END DO
         END DO
      ENDIF

      IF( ln_ndepo ) THEN
         DO jj = JRANGE
            DO ji = IRANGE
               ! conversion from KgN/m2/month to molC/L/s
               zfact = rfact2 / rno3 / ( 14. * rmtss ) / fse3t(ji,jj,KSURF)
               zno3dep(ji,jj) =  zfact * no3dep(ji,jj) 
               znh4dep(ji,jj) =  zfact * nh4dep(ji,jj) 
               !
               tra(ji,jj,1,jpno3) = tra(ji,jj,1,jpno3) + zno3dep(ji,jj) 
               tra(ji,jj,1,jpnh4) = tra(ji,jj,1,jpnh4) + znh4dep(ji,jj) 
               tra(ji,jj,1,jptal) = tra(ji,jj,1,jptal) + rno3 * ( znh4dep(ji,jj) - zno3dep(ji,jj) )
            END DO
         END DO
      ENDIF

      ! Add the external input of iron which is 3D distributed
      ! (dust, river and sediment mobilization)
      ! ------------------------------------------------------

      IF( ln_sedinput ) THEN     
         DO jk = KRANGE
            DO jj = JRANGE
               DO ji = IRANGE
                  tra(ji,jj,jk,jpfer) = tra(ji,jj,jk,jpfer) + ironsed(ji,jj,jk) * rfact2
              END DO
            END DO
         END DO
      ENDIF


#if ! defined key_sed
      ! Loss of biogenic silicon, Caco3 organic carbon in the sediments. 
      ! First, the total loss is computed.
      ! The factor for calcite comes from the alkalinity effect
      ! -------------------------------------------------------------

      ! Initialisation of variables used to compute Sinking Speed
      zsumsedsi  = 0.e0
      zsumsedpo4 = 0.e0
      zsumsedcal = 0.e0
      DO jj = JRANGE
         DO ji = IRANGE
            zfact = e1t(ji,jj) * e2t(ji,jj) / rday * tmask_i(ji,jj)
            zws3  = zfact * wsbio3(ji,jj,ikt)
            zws4  = zfact * wsbio4(ji,jj,ikt)
            zwsc  = zfact * wscal (ji,jj,ikt)
# if defined key_kriest
            zsumsedsi  = zsumsedsi  + trn(ji,jj,KSED,jpdsi) * zwsc
            zsumsedpo4 = zsumsedpo4 + trn(ji,jj,KSED,jppoc) * zws3
# else
            zsumsedsi  = zsumsedsi  + trn(ji,jj,KSED,jpdsi) * zws4
            zsumsedpo4 = zsumsedpo4 + ( trn(ji,jj,KSED,jpgoc) * zws4   &
               &                      + trn(ji,jj,KSED,jppoc) * zws3 )
# endif
         END DO
      END DO
      IF( .NOT. ln_sedmeta ) THEN
         !
         DO jj = JRANGE
            DO ji = IRANGE
               zfact = e1t(ji,jj) * e2t(ji,jj) / rday * tmask_i(ji,jj)
               zwsc  = zfact * wscal(ji,jj,ikt)
               zsumsedcal = zsumsedcal + trn(ji,jj,KSED,jpcal) * zwsc * 2.e0
            END DO
         END DO
         !
      ELSE
         ! For calcite, burial efficiency is made a function of saturation
         DO jj = JRANGE
            DO ji = IRANGE
               zfact      = e1t(ji,jj) * e2t(ji,jj) / rday * tmask_i(ji,jj)
               zwsc       = zfact * wscal (ji,jj,ikt)
               zfactcal   = MIN( excess(ji,jj,ikt), 0.2 )
               zfactcal   = MIN( 1., 1.3 * ( 0.2 - zfactcal ) / ( 0.4 - zfactcal ) )
               zsumsedcal = zsumsedcal + trn(ji,jj,KSED,jpcal) * zwsc * 2.e0 * zfactcal
            END DO
         END DO
      ENDIF

      IF( lk_mpp ) THEN
         CALL mpp_sum( zsumsedsi  )   ! sums over the global domain
         CALL mpp_sum( zsumsedcal )   ! sums over the global domain
         CALL mpp_sum( zsumsedpo4 )   ! sums over the global domain
      ENDIF

#endif

      IF( .NOT. ln_sedmeta ) THEN
        ! Then this loss is scaled at each bottom grid cell for
        ! equilibrating the total budget of silica in the ocean.
        ! Thus, the amount of silica lost in the sediments equal
        ! the supply at the surface (dust+rivers)
        ! ------------------------------------------------------

        DO jj = JRANGE
           DO ji = IRANGE
              zconctmp = trn(ji,jj,KSED,jpdsi) * xstep / fse3t(ji,jj,KSED)   &
# if ! defined key_kriest
     &               * wscal (ji,jj,ikt)
# else
     &               * wsbio4(ji,jj,ikt)
# endif
              tra(ji,jj,ikt,jpdsi) = tra(ji,jj,ikt,jpdsi) - zconctmp

#if ! defined key_sed
              tra(ji,jj,ikt,jpsil) = tra(ji,jj,ikt,jpsil) + zconctmp   &
                &      * 0.98
!                &      * ( 1.- ( sumdepsi + rivalkinput / ryyss / 6. ) / zsumsedsi )
#endif
           END DO
        END DO

        DO jj = JRANGE
           DO ji = IRANGE
              zconctmp = trn(ji,jj,KSED,jpcal) * wscal(ji,jj,ikt) * xstep / fse3t(ji,jj,KSED)
              tra(ji,jj,ikt,jpcal) = tra(ji,jj,ikt,jpcal) - zconctmp

#if ! defined key_sed
              tra(ji,jj,ikt,jptal) = tra(ji,jj,ikt,jptal) + zconctmp   &
                 &   * 0.85 * 2.0
!                 &   * ( 1.- ( rivalkinput / ryyss ) / zsumsedcal ) * 2.e0
              tra(ji,jj,ikt,jpdic) = tra(ji,jj,ikt,jpdic) + zconctmp   &
                 &   * 0.85
!                 &   * ( 1.- ( rivalkinput / ryyss ) / zsumsedcal )
#endif
           END DO
        END DO

        DO jj = JRANGE
           DO ji = IRANGE
              zfact = xstep / fse3t(ji,jj,KSED)
              zws3 = wsbio3(ji,jj,ikt) * zfact
              zws4 = wsbio4(ji,jj,ikt) * zfact
# if ! defined key_kriest
              zconctmp  = trn(ji,jj,KSED,jpgoc)
              zconctmp2 = trn(ji,jj,KSED,jppoc)
              tra(ji,jj,ikt,jpgoc) = tra(ji,jj,ikt,jpgoc) - zconctmp  * zws4
              tra(ji,jj,ikt,jppoc) = tra(ji,jj,ikt,jppoc) - zconctmp2 * zws3
#if ! defined key_sed
              tra(ji,jj,ikt,jpdoc) = tra(ji,jj,ikt,jpdoc)    &
                &                  + ( zconctmp * zws4 + zconctmp2 * zws3 )  &
                &                    * 0.92
!              &                 * ( 1.- rivpo4input / (ryyss * zsumsedpo4 ) )
#endif
              tra(ji,jj,ikt,jpbfe) = tra(ji,jj,ikt,jpbfe) - trn(ji,jj,KSED,jpbfe) * zws4
              tra(ji,jj,ikt,jpsfe) = tra(ji,jj,ikt,jpsfe) - trn(ji,jj,KSED,jpsfe) * zws3

# else
              zconctmp  = trn(ji,jj,KSED,jpnum)
              zconctmp2 = trn(ji,jj,KSED,jppoc)
              tra(ji,jj,ikt,jpnum) = tra(ji,jj,ikt,jpnum) - zconctmp  * zws4
              tra(ji,jj,ikt,jppoc) = tra(ji,jj,ikt,jppoc) - zconctmp2 * zws3
#if ! defined key_sed
              tra(ji,jj,ikt,jpdoc) = tra(ji,jj,ikt,jpdoc) + zconctmp2 * zws3 &
                             &       * 0.92
!                             &      * ( 1.- rivpo4input / (ryyss * zsumsedpo4 ) )
#endif
              tra(ji,jj,ikt,jpsfe) = tra(ji,jj,ikt,jpsfe) - trn(ji,jj,KSED,jpsfe) * zws3

# endif
           END DO
        END DO

        ! Potential nitrogen fixation dependant on temperature and iron
        ! -------------------------------------------------------------

        DO jk = KRANGE
           DO jj = JRANGE
              DO ji = IRANGE
                 zlim = ( 1.- xnanono3(ji,jj,jk) - xnanonh4(ji,jj,jk) )
                 IF( zlim <= 0.2 )   zlim = 0.01
                 nitrpot(ji,jj,jk) = MAX( 0.e0, ( 0.6 * tgfunc(ji,jj,jk) - 2.15 ) / rday )   &
                 &                  * zlim * rfact2 * trn(ji,jj,K,jpfer)   &
                 &                  / ( conc3 + trn(ji,jj,K,jpfer) ) &
                 &                  * ( 1.- EXP( -etot(ji,jj,jk) / diazolight) )
              END DO
           END DO
        END DO


        DO jk = KRANGE
           DO jj = JRANGE
              DO ji = IRANGE
                 zfact = nitrpot(ji,jj,jk) * nitrfix
                 tra(ji,jj,jk,jpnh4) = tra(ji,jj,jk,jpnh4) + zfact
                 tra(ji,jj,jk,jpoxy) = tra(ji,jj,jk,jpoxy) + zfact   * o2nit
                 tra(ji,jj,jk,jppo4) = tra(ji,jj,jk,jppo4) + 30./ 46.* zfact
              END DO
           END DO
        END DO
       !
     ELSE

       ! Computation of the sediment denitrification proportion:
       ! The metamodel from midlleburg (2006) is being used
       ! Computation of the fraction of organic matter that is permanently buried from Dunne's model
       ! -------------------------------------------------------

       DO jj = JRANGE
          DO ji = IRANGE
# if defined key_kriest
             zflx =    trn(ji,jj,KSED,jppoc) * wsbio3(ji,jj,ikt)   * 1E3 * 1E6 / 1E4
# else
             zflx = (  trn(ji,jj,KSED,jpgoc) * wsbio4(ji,jj,ikt)   &
                &    + trn(ji,jj,KSED,jppoc) * wsbio3(ji,jj,ikt) )  * 1E3 * 1E6 / 1E4
#endif
             zflx  = LOG10( MAX( 1E-3, zflx ) )
             zo2   = LOG10( MAX( 10. , trn(ji,jj,KSED,jpoxy) * 1E6 ) )
             zno3  = LOG10( MAX( 1.  , trn(ji,jj,KSED,jpno3) * 1E6 * rno3 ) )
             zdep  = LOG10( fsdepw(ji,jj,ikt+1) )
             zdenit2d(ji,jj) = -2.2567 - 1.185 * zflx - 0.221 * zflx**2 &
                &              - 0.3995 * zno3 * zo2  + 1.25 * zno3    &
             &                 + 0.4721 * zo2 - 0.0996 * zdep + 0.4256 * zflx * zo2
             zdenit2d(ji,jj) = 10.0**( zdenit2d(ji,jj) )
             !
# if defined key_kriest
             zflx = (  trn(ji,jj,KSED,jppoc) * wsbio3(ji,jj,ikt) * 1E6
#else
             zflx = (  trn(ji,jj,KSED,jpgoc) * wsbio4(ji,jj,ikt)   &
               &     + trn(ji,jj,KSED,jppoc) * wsbio3(ji,jj,ikt) ) * 1E6
#endif
             zbureff(ji,jj) = 0.013 + 0.53 * zflx**2 / ( 7.0 + zflx )**2
          END DO
       END DO 


      ! This loss is scaled at each bottom grid cell for equilibrating the total budget of silica in the ocean.
      ! Thus, the amount of silica lost in the sediments equal the supply at the surface (dust+rivers)
      ! ------------------------------------------------------
#if ! defined key_sed
      zrivsil =  1. - ( sumdepsi + rivalkinput / ryyss ) &
          &     / ( zsumsedsi + rtrn )
#endif

      DO jj = JRANGE
         DO ji = IRANGE
            zdep = xstep / fse3t(ji,jj,KSED) 
            zws4 = wsbio4(ji,jj,ikt) * zdep
            zwsc = wscal (ji,jj,ikt) * zdep
# if defined key_kriest
            zsiloss = trn(ji,jj,KSED,jpdsi) * zws4
# else
            zsiloss = trn(ji,jj,KSED,jpdsi) * zwsc
# endif
            zcaloss = trn(ji,jj,KSED,jpcal) * zwsc
            !
            tra(ji,jj,ikt,jpdsi) = tra(ji,jj,ikt,jpdsi) - zsiloss
            tra(ji,jj,ikt,jpcal) = tra(ji,jj,ikt,jpcal) - zcaloss
#if ! defined key_sed
            tra(ji,jj,ikt,jpsil) = tra(ji,jj,ikt,jpsil) + zsiloss * zrivsil 
            zfactcal = MIN( excess(ji,jj,ikt), 0.2 )
            zfactcal = MIN( 1., 1.3 * ( 0.2 - zfactcal ) / ( 0.4 - zfactcal ) )
            zrivalk  =  1. - ( rivalkinput / ryyss ) * zfactcal / ( zsumsedcal + rtrn )
            tra(ji,jj,ikt,jptal) =  tra(ji,jj,ikt,jptal) + zcaloss * zrivalk * 2.0
            tra(ji,jj,ikt,jpdic) =  tra(ji,jj,ikt,jpdic) + zcaloss * zrivalk
#endif
         END DO
      END DO



      DO jj = JRANGE
         DO ji = IRANGE
            zdep = xstep / fse3t(ji,jj,KSED) 
            zws4 = wsbio4(ji,jj,ikt) * zdep
            zws3 = wsbio3(ji,jj,ikt) * zdep
            zrivno3 = 1. - zbureff(ji,jj)
# if ! defined key_kriest
            tra(ji,jj,ikt,jpgoc) = tra(ji,jj,ikt,jpgoc) - trn(ji,jj,KSED,jpgoc) * zws4 
            tra(ji,jj,ikt,jppoc) = tra(ji,jj,ikt,jppoc) - trn(ji,jj,KSED,jppoc) * zws3
            tra(ji,jj,ikt,jpbfe) = tra(ji,jj,ikt,jpbfe) - trn(ji,jj,KSED,jpbfe) * zws4
            tra(ji,jj,ikt,jpsfe) = tra(ji,jj,ikt,jpsfe) - trn(ji,jj,KSED,jpsfe) * zws3
            zwstpoc              = trn(ji,jj,KSED,jpgoc) * zws4 + trn(ji,jj,KSED,jppoc) * zws3
# else
            tra(ji,jj,ikt,jpnum) = tra(ji,jj,ikt,jpnum) - trn(ji,jj,KSED,jpnum) * zws4 
            tra(ji,jj,ikt,jppoc) = tra(ji,jj,ikt,jppoc) - trn(ji,jj,KSED,jppoc) * zws3
            tra(ji,jj,ikt,jpsfe) = tra(ji,jj,ikt,jpsfe) - trn(ji,jj,KSED,jpsfe) * zws3
            zwstpoc = trn(ji,jj,KSED,jppoc) * zws3 
# endif

#if ! defined key_sed
            ! The 0.5 factor in zpdenit and zdenitt is to avoid negative NO3 concentration after both denitrification
            ! in the sediments and just above the sediments. Not very clever, but simpliest option.
            zpdenit  = MIN( 0.5 * ( trn(ji,jj,KSED,jpno3) - rtrn ) / rdenit &
               &          , zdenit2d(ji,jj) * zwstpoc * zrivno3 )
            z1pdenit = zwstpoc * zrivno3 - zpdenit
            zolimit = MIN( ( trn(ji,jj,KSED,jpoxy) - rtrn ) / o2ut &
               &         , z1pdenit * ( 1.- nitrfac(ji,jj,ikt) ) )
            zdenitt = MIN(  0.5 * ( trn(ji,jj,KSED,jpno3) - rtrn ) / rdenit &
               &         , z1pdenit * nitrfac(ji,jj,ikt) )
            tra(ji,jj,ikt,jpdoc) = tra(ji,jj,ikt,jpdoc) + z1pdenit - zolimit - zdenitt
            tra(ji,jj,ikt,jppo4) = tra(ji,jj,ikt,jppo4) + zpdenit + zolimit + zdenitt
            tra(ji,jj,ikt,jpnh4) = tra(ji,jj,ikt,jpnh4) + zpdenit + zolimit + zdenitt
            tra(ji,jj,ikt,jpno3) = tra(ji,jj,ikt,jpno3) - rdenit * (zpdenit + zdenitt)
            tra(ji,jj,ikt,jpoxy) = tra(ji,jj,ikt,jpoxy) - zolimit * o2ut
            tra(ji,jj,ikt,jptal) = tra(ji,jj,ikt,jptal) + rno3 * (zolimit + (1.+rdenit) * (zpdenit + zdenitt) )
            tra(ji,jj,ikt,jpdic) = tra(ji,jj,ikt,jpdic) + zpdenit + zolimit + zdenitt
#endif
         END DO
      END DO


      ! Potential nitrogen fixation dependant on temperature and iron
      ! -------------------------------------------------------------

      DO jk = KRANGE
         DO jj = JRANGE
            DO ji = IRANGE
               zlim = ( 1.- xnanono3(ji,jj,jk) - xnanonh4(ji,jj,jk) )
               IF( zlim <= 0.2 )   zlim = 0.01
               zfact = zlim * rfact2
               ztrfer = trn(ji,jj,K,jpfer) / ( concfediaz + trn(ji,jj,K,jpfer) )
               ztrpo4 = trn(ji,jj,K,jppo4) / ( concnnh4   + trn(ji,jj,K,jppo4) ) 
               zlight =  ( 1.- EXP( -etot(ji,jj,jk) / diazolight ) ) 
               nitrpot(ji,jj,jk) =  MAX( 0.e0, ( 0.6 * tgfunc(ji,jj,jk) - 2.15 ) / rday )   &
                 &                *  zfact * MIN( ztrfer, ztrpo4 ) * zlight
            END DO
         END DO 
      END DO


      ! Nitrogen change due to nitrogen fixation
      ! ----------------------------------------

      DO jk = KRANGE
         DO jj = JRANGE
            DO ji = IRANGE
               zfact = nitrpot(ji,jj,jk) * nitrfix
               tra(ji,jj,jk,jpnh4) = tra(ji,jj,jk,jpnh4) +             zfact
               tra(ji,jj,jk,jptal) = tra(ji,jj,jk,jptal) + rno3      * zfact
               tra(ji,jj,jk,jpoxy) = tra(ji,jj,jk,jpoxy) + o2nit     * zfact 
               tra(ji,jj,jk,jppo4) = tra(ji,jj,jk,jppo4) + concdnh4 / ( concdnh4 + trn(ji,jj,K,jppo4) ) &
               &                     * 0.002 * trn(ji,jj,K,jpdoc) * xstep
            END DO
         END DO
      END DO


     ENDIF

#if defined key_trc_diaadd
        zrfact2 = 1.e+3 * rfact2r
        DO jj = JRANGE
           DO ji = IRANGE
              zmsk = zrfact2 * fse3t(ji,jj,KSURF) * tmask(ji,jj,KSURF)
              trc2d(ji,jj,jp_sildep)   = sidep(ji,jj)  * zmsk                ! iron deposition      
              trc2d(ji,jj,jp_po4dep)   = po4dep(ji,jj) * po4r * zmsk                ! iron deposition      
              trc2d(ji,jj,jp_no3dep )  = zno3dep(ji,jj) * rno3 * zmsk                ! iron deposition      
              trc2d(ji,jj,jp_nh4dep )  = znh4dep(ji,jj) * rno3 * zmsk                ! iron deposition      
              trc2d(ji,jj,jp_nfix   )  = nitrpot(ji,jj,1) * rno3 * zmsk * nitrfix       ! nitrogen fixation at surface
           END DO
        END DO

        DO jk = KRANGE
           DO jj = JRANGE
              DO ji = IRANGE
                 zmsk = zrfact2 * fse3t(ji,jj,K) * tmask(ji,jj,K)
                 trc3d(ji,jj,K,jp_nfixo2 )  = nitrpot(ji,jj,jk) * rno3 * zmsk * nitrfix * o2nit  ! O2 production by Nfix
                 trc3d(ji,jj,K,jp_irondep)  = irondep(ji,jj,jk) * zmsk                ! iron flux from dust
                 trc3d(ji,jj,K,jp_ironsed ) = ironsed(ji,jj,jk) * 1e+3 * tmask(ji,jj,K)  ! iron from  sediment
             END DO
           END DO
        ENDDO

# endif
      !
       IF(ln_ctl)   THEN  ! print mean trends (used for debugging)
         WRITE(charout, FMT="('sed ')")
         CALL prt_ctl_trc_info(charout)
         CALL prt_ctl_trc( charout, ltra='tra')
!         CALL prt_ctl_trc(tab4d=trn, mask=tmask, clinfo=ctrcnm)
       ENDIF

   END SUBROUTINE p4z_sed

   INTEGER FUNCTION p4z_sed_alloc()
      !!----------------------------------------------------------------------
      !!                     ***  ROUTINE p4z_sed_alloc  ***
      !!----------------------------------------------------------------------
      INTEGER ::   ierr(7)        ! Local variables

      ierr(:) = 0
      !
#ifdef NEMO
      ALLOCATE( dust(PRIV_2D_BIOARRAY), dustmo(PRIV_2D_BIOARRAY,2), STAT= ierr(1) )
#else
      ALLOCATE( dust(PRIV_2D_BIOARRAY), STAT= ierr(1) )
#endif
      ALLOCATE( irondep(PRIV_3D_BIOARRAY),po4dep(PRIV_2D_BIOARRAY), STAT= ierr(2) )
      ALLOCATE( sidep(PRIV_2D_BIOARRAY), STAT= ierr(3) )
      ALLOCATE( no3dep(PRIV_2D_BIOARRAY), nh4dep(PRIV_2D_BIOARRAY), STAT=ierr(4) )
      ALLOCATE( rivinp(PRIV_2D_BIOARRAY), cotdep(PRIV_2D_BIOARRAY), STAT=ierr(6) )
      ALLOCATE( nitrpot(PRIV_3D_BIOARRAY), ironsed(PRIV_3D_BIOARRAY), STAT=ierr(7) )
      !
      p4z_sed_alloc = MAXVAL( ierr )
      !
      IF( p4z_sed_alloc /= 0 )   CALL ctl_warn('p4z_sed_alloc: failed to allocate arrays')
      !
   END FUNCTION p4z_sed_alloc

# if defined NEMO

   SUBROUTINE p4z_sbc(kt)

      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE p4z_sbc  ***
      !!
      !! ** Purpose :   Read and interpolate the external sources of 
      !!                nutrients
      !!
      !! ** Method  :   Read the files and interpolate the appropriate variables
      !!
      !! ** input   :   external netcdf files
      !!
      !!----------------------------------------------------------------------
      !! * arguments
      INTEGER, INTENT( in  ) ::   kt   ! ocean time step
      !! * Local declarations
      INTEGER ::   &
         imois, imois2,       &  ! temporary integers
         i15  , iman             !    "          "
      REAL(wp) ::   &
         zxy                     !    "         "


      !!---------------------------------------------------------------------

      ! Initialization
      ! --------------

      i15 = nday / 16
      iman  = INT( raamo )
      imois = nmonth + i15 - 1
      IF( imois == 0 ) imois = iman
      imois2 = nmonth

      ! 1. first call kt=nit000
      ! -----------------------

      IF( kt == nit000 ) THEN
         ! initializations
         nflx1  = 0
         nflx11 = 0
         ! open the file
         IF(lwp) THEN
            WRITE(numout,*) ' '
            WRITE(numout,*) ' **** Routine p4z_sbc'
         ENDIF
         CALL iom_open ( 'dust.orca.nc', numdust )
      ENDIF


     ! Read monthly file
      ! ----------------

      IF( kt == nit000 .OR. imois /= nflx1 ) THEN

         ! Calendar computation

         ! nflx1 number of the first file record used in the simulation
         ! nflx2 number of the last  file record

         nflx1 = imois
         nflx2 = nflx1+1
         nflx1 = MOD( nflx1, iman )
         nflx2 = MOD( nflx2, iman )
         IF( nflx1 == 0 )   nflx1 = iman
         IF( nflx2 == 0 )   nflx2 = iman
         IF(lwp) WRITE(numout,*) 'first record file used nflx1 ',nflx1
         IF(lwp) WRITE(numout,*) 'last  record file used nflx2 ',nflx2

         ! Read monthly fluxes data

         ! humidity
         CALL iom_get ( numdust, jpdom_data, 'dust', dustmo(:,:,1), nflx1 )
         CALL iom_get ( numdust, jpdom_data, 'dust', dustmo(:,:,2), nflx2 )


      ENDIF

     ! 3. at every time step interpolation of fluxes
      ! ---------------------------------------------

      zxy = FLOAT( nday + 15 - 30 * i15 ) / 30
      dust(:,:) = ( (1.-zxy) * dustmo(:,:,1) + zxy * dustmo(:,:,2) )

      IF( kt == nitend ) CALL iom_close (numdust)

#else

   SUBROUTINE p4z_sbc(kt)

      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE p4z_sbc  ***
      !!
      !! ** Purpose :   Read and interpolate the external sources of 
      !!                nutrients
      !!
      !! ** Method  :   Read the files and interpolate the appropriate variables
      !!
      !! ** input   :   external netcdf files
      !!
      !!----------------------------------------------------------------------
      !! * arguments
      INTEGER, INTENT( in  ) ::   kt   ! ocean time step
      !
      INTEGER :: ji, jj, jk
      INTEGER, PARAMETER :: jpmois = 12
      INTEGER :: irec1, irec2, i15
      INTEGER :: nyear, nday, nmonth
      REAL    :: zpdtan, zpdtmo, zdemi, zt
      REAL    :: zxy, zjulian, zsec

      IF( kt == nit000 .AND. lwp ) THEN
        WRITE(numout,*) ' '
        WRITE(numout,*) ' Number of days per year in file year2daydta = ', year2daydta 
        WRITE(numout,*) ' '
      ENDIF

      zpdtan = ( year2daydta * day2sec ) / rdt
      zpdtmo = zpdtan / float( jpmois )
      zdemi  = zpdtmo / 2.
      zt     = ( float( kt ) + zdemi) / zpdtmo
      

      !  recherche de l'indice des enregistrements
      !  du modele dynamique encadrant le pas de temps kt.
      !  --------------------------------------------------
      irec1 = int( zt )
      irec2 = irec1 + 1
      irec1 = MOD( irec1, jpmois )
      IF ( irec1 == 0 ) irec1 = jpmois
      irec2 = MOD( irec2, jpmois )
      IF ( irec2 == 0 ) irec2 = jpmois

      zxy = zt - float(int ( zt ) )
      !
      IF( ln_dust ) THEN
         DO jj = JRANGE
            DO ji = IRANGE
               dust(ji,jj) = ( 1. - zxy ) * dustmo(ji,jj,irec1)  &
                  &               + zxy   * dustmo(ji,jj,irec2)
            END DO
          END DO
      ENDIF
      !
      IF( ln_ndepo ) THEN
         DO jj = JRANGE
            DO ji = IRANGE
               no3dep(ji,jj) = ( 1. - zxy ) * no3depmo(ji,jj,irec1)  &
                 &                  + zxy   * no3depmo(ji,jj,irec2)
               !
               nh4dep(ji,jj) = ( 1. - zxy ) * nh4depmo(ji,jj,irec1)  &
                 &                   + zxy  * nh4depmo(ji,jj,irec2)
            END DO
          END DO
      ENDIF
#endif

   END SUBROUTINE p4z_sbc

   SUBROUTINE p4z_sed_nam

      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE p4z_sed_init  ***
      !!
      !! ** Purpose :   Initialization of the external sources of nutrients
      !!
      !! ** Method  :   Read the files and compute the budget
      !!      called at the first timestep (nittrc000)
      !!
      !! ** input   :   external netcdf files
      !!
      !!----------------------------------------------------------------------
      NAMELIST/nampissed/ ln_dust, ln_river, ln_ndepo, ln_sedinput, &
         &                sedfeinput, dustsolub, wdust, mfrac, &
         &                nitrfix, diazolight, ln_sedmeta, concfediaz


      REWIND( numnatp )                     ! read numnatp
      READ  ( numnatp, nampissed )

      IF(lwp) THEN
         WRITE(numout,*) ' '
         WRITE(numout,*) ' Namelist : nampissed '
         WRITE(numout,*) ' ~~~~~~~~~~~~~~~~~ '
         WRITE(numout,*) '    Dust input from the atmosphere                         ln_dust  = ', ln_dust
         WRITE(numout,*) '    River input of nutrients                               ln_river    = ', ln_river
         WRITE(numout,*) '    Atmospheric deposition of N                            ln_ndepo    = ', ln_ndepo
         WRITE(numout,*) '    Fe input from sediments                                ln_sedinput = ', ln_sedinput
         WRITE(numout,*) '    Coastal release of Iron                                sedfeinput  = ', sedfeinput
         WRITE(numout,*) '    Solubility of the dust                                 dustsolub   = ', dustsolub
         WRITE(numout,*) '    Fe Mineral fraction of the dust                        mfrac       = ', mfrac
         WRITE(numout,*) '     Dust sinking speed                                    wdust       = ', wdust
         WRITE(numout,*) '    nitrogen fixation sensitivty to light                  diazolight  = ', diazolight
         WRITE(numout,*) '    nitrogen fixation rate                                 nitrfix     = ', nitrfix
         WRITE(numout,*) '    New sediment param : metamodel from midlleburg (y/n)   ln_sedmeta  = ', ln_sedmeta
         WRITE(numout,*) '    fe half-saturation cste for diazotrophs                concfediaz  = ', concfediaz
      ENDIF


   END SUBROUTINE p4z_sed_nam

   SUBROUTINE p4z_sed_init

      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE p4z_sed_init  ***
      !!
      !! ** Purpose :   Initialization of the external sources of nutrients
      !!
      !! ** Method  :   Read the files and compute the budget
      !!      called at the first timestep (nittrc000)
      !!
      !! ** input   :   external netcdf files
      !!
      !!----------------------------------------------------------------------

      INTEGER ::   ji, jj, jk, jm
      INTEGER , PARAMETER ::   jpmois = 12, jpan = 1
      INTEGER :: numriv, numbath, numdep


      REAL(wp) ::   zcoef
      REAL(wp) ::   expide, denitide,zmaskt
      REAL(wp) , DIMENSION(PRIV_2D_BIOARRAY)     ::   riverdoc, river, ndepo
      REAL(wp) , DIMENSION(PRIV_2D_BIOARRAY,jpmois)    ::   zdustmo

      DO jk = KRANGE
         DO jj = JRANGE
            DO ji = IRANGE
               irondep(ji,jj,jk) = 0.e0          
               ironsed(ji,jj,jk) = 0.e0          
               nitrpot(ji,jj,jk) = 0.e0          
            END DO
         END DO
     END DO
     DO jj = JRANGE
        DO ji = IRANGE
             sidep (ji,jj) = 0.e0          
             po4dep(ji,jj) = 0.e0          
             no3dep(ji,jj) = 0.e0
             nh4dep(ji,jj) = 0.e0
        END DO
     END DO

      ! Dust input from the atmosphere
      ! ------------------------------
      IF( ln_dust ) THEN 
!         IF(lwp) WRITE(numout,*) '    Initialize dust input from atmosphere '
!         IF(lwp) WRITE(numout,*) '    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ '
!         CALL iom_open ( 'dust.orca.nc', numdust )
!         DO jm = 1, jpmois
!            CALL iom_get( numdust, jpdom_data, 'dust', zdustmo(:,:,jm), jm )
!         END DO
!         CALL iom_close( numdust )
         !
         ! total atmospheric supply of Si
         ! ------------------------------
         sumdepsi = 0.e0
!         DO jm = 1, jpmois
!            DO jj = JRANGE
!               DO ji = IRANGE
!                  sumdepsi = sumdepsi + zdustmo(ji,jj,jm) / (12.*rmtss) * 8.8        &
!                     &     * 0.075/28.1 * e1t(ji,jj) * e2t(ji,jj) * tmask(ji,jj,KSURF) * tmask_i(ji,jj)
!              END DO
!            END DO
!         END DO
!        IF( lk_mpp )  CALL mpp_sum( sumdepsi )  ! sum over the global domain
      ELSE
         sumdepsi = 0.e0
      ENDIF

      ! Nutrient input from rivers
      ! --------------------------
      IF( ln_river ) THEN
!        IF(lwp) WRITE(numout,*) '    Initialize the nutrient input by rivers from river.orca.nc file'
!         IF(lwp) WRITE(numout,*) '    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
!         CALL iom_open ( 'river.orca.nc', numriv )
!         CALL iom_get  ( numriv, jpdom_data, 'riverdic', river   (:,:), jpan )
!         CALL iom_get  ( numriv, jpdom_data, 'riverdoc', riverdoc(:,:), jpan )
!         CALL iom_close( numriv )

        ! Number of seconds per year and per month
!      ryyss = nyear_len(1) * rday
!      rmtss = ryyss / raamo

         DO jj = JRANGE
            DO ji = IRANGE
               river   (ji,jj) = 0. 
               riverdoc(ji,jj) = 0.
            END DO
         END DO

         ! N/P and Si releases due to coastal rivers
         ! -----------------------------------------
         DO jj = JRANGE
            DO ji = IRANGE
               zcoef = ryyss * e1t(ji,jj) * e2t(ji,jj) * fse3t(ji,jj,KSURF) * tmask(ji,jj,KSURF) * tmask_i(ji,jj)
               cotdep(ji,jj) =  river(ji,jj)                  *1E9 / ( 12. * zcoef + rtrn )
               rivinp(ji,jj) = (river(ji,jj)+riverdoc(ji,jj)) *1E9 / ( 31.6* zcoef + rtrn )
            END DO
         END DO
         rivpo4input = 0.e0
         rivalkinput = 0.e0
!         DO jj = JRANGE
!            DO ji = IRANGE
!               zcoef = cvol(ji,jj,KSURF) * ryyss
!               rivpo4input = rivpo4input + rivinp(ji,jj) * zcoef
!               rivalkinput = rivalkinput + cotdep(ji,jj) * zcoef
!            END DO
!         END DO
!         IF( lk_mpp ) THEN
!            CALL mpp_sum( rivpo4input )  ! sum over the global domain
!            CALL mpp_sum( rivalkinput )  ! sum over the global domain
!         ENDIF
      ELSE
         rivpo4input = 0.e0
         rivalkinput = 0.e0
      ENDIF

      ! Nutrient input from dust
      ! ------------------------
      IF( ln_ndepo ) THEN
!        IF(lwp) WRITE(numout,*) '    Initialize the nutrient input by dust from ndeposition.orca.nc'
!        IF(lwp) WRITE(numout,*) '    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
!         CALL iom_open ( 'ndeposition.orca.nc', numdep )
!         CALL iom_get  ( numdep, jpdom_data, 'ndep', ndepo(:,:), jpan )
!         CALL iom_close( numdep )
!      
!         DO jj = JRANGE
!            DO ji = IRANGE
!               zcoef         = 14E6*ryyss*fse3t(ji,jj,KSURF) 
!               nitdep(ji,jj) = 7.6 * ndepo(ji,jj) / ( zcoef + rtrn )
!            END DO
!         END DO
         nitdepinput = 0.e0
!         DO jj = JRANGE
!            DO ji = IRANGE
!               zcoef = cvol(ji,jj,KSURF) * ryyss
!               nitdepinput = nitdepinput + nitdep(ji,jj) * zcoef
!            END DO
!         END DO
!         IF( lk_mpp ) CALL mpp_sum( nitdepinput )  ! sum over the global domain
      ELSE
         nitdepinput = 0.e0
      ENDIF

      ! Iron input from sediment
      IF( ln_sedinput ) THEN

         DO jj = JRANGE
            DO ji = IRANGE
               cmask(ji,jj,jpk) = 1
            ENDDO
         ENDDO
         DO jk = 2, N
            DO jj = JRANGE-1
               DO ji = IRANGE-1
                  IF( tmask(ji,jj,K) /= 0. ) THEN
                     zmaskt = tmask(ji+1,jj,K) * tmask(ji-1,jj  ,K) * tmask(ji,jj+1,K)    &
                        &                      * tmask(ji  ,jj-1,K) * tmask(ji,jj  ,K)
                     IF( zmaskt == 0. )   cmask(ji,jj,K ) = 0.1
                  ENDIF
               END DO
            END DO
         END DO
         DO jk = KRANGE
            DO jj = JRANGE
               DO ji = IRANGE
                  expide   = MIN( 8.,( fsdept(ji,jj,K) / 500. )**(-1.5) )
                  denitide = -0.9543 + 0.7662 * LOG( expide ) - 0.235 * LOG(expide )**2
                  cmask(ji,jj,jk) = cmask(ji,jj,jk) * MIN( 1., EXP( denitide ) / 0.5 )
               END DO
            END DO
         END DO
      
        ! Coastal supply of iron
        ! -------------------------
        DO jk = KRANGE
           DO jj = JRANGE
              DO ji = IRANGE
                 ironsed(ji,jj,jk) = sedfeinput * cmask(ji,jj,jk) / (fse3t(ji,jj,K) * rday )
              END DO
          END DO
         END DO
         !
      ENDIF

   END SUBROUTINE p4z_sed_init

#else
   !!======================================================================
   !!  Dummy module :                                   No PISCES bio-model
   !!======================================================================
CONTAINS
   SUBROUTINE p4z_sed                         ! Empty routine
   END SUBROUTINE p4z_sed
#endif 

   !!======================================================================
END MODULE  p4zsed
