#if defined RAPPEL_VBAR  || defined RAPPEL_VBARCLINE
!----- for short barotropic rossby wave nudging
      real            mz, kx, omegaY
      real            vbaramp, ltourb, ytourb, onedeg, cvit, lrossby
      real            vrap(GLOBAL_2D_ARRAY)
      real            raptab(GLOBAL_2D_ARRAY)
      real            yshape_tourb(GLOBAL_2D_ARRAY)
      real            lx
      real            rappel_rate,xrmid,rhalfw,xrdecay,rhalfn,ypdecay
      common /rappel_vbar/ vrap, raptab,omegaY, vbaramp

!#ifdef RAPPEL_VBARCLINE
      common /rappel_vbarcline/  zshape_tourb
      real           zshape_tourb(N), zdepth_tourb, zcent_tourb
      real           amp_tourb
      real           amp_clin1, amp_clin2, h0
!#endif        /* RAPPEL_VBARCLINE */

#endif        /* RAPPEL_VBAR || RAPPEL_VBARCLINE */

