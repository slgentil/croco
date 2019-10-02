!
!====================================================================
!           Wave maker for wave-resolving simulations: 
!              monochromatic, bichromatic or JONSMAP
!====================================================================
!
!--------------------------------------------------------------------
!  Configurations
!--------------------------------------------------------------------
!
!  default parameters
!
        wp=11.            ! period
        wa=0.4            ! amplitude
        wd=0.             ! incidence angle 
        wds=0.            ! directional spread
        gamma=3.3         ! JONSWAP peakedness parameter
!
!  Set configuration parameters
!
#   ifdef FLUME_WAVES
#    define WAVE_MAKER_SPECTRUM

#   elif defined RIP && !defined MRL_WCI
#    define WAVE_MAKER_JONSWAP
#    define WAVE_MAKER_OBLIQUE
#    define WAVE_MAKER_DSPREAD
        wp=11.            ! period
        wa=0.4            ! amplitude
        wd=-20.           ! incidence angle (deg)
        wds=15.           ! directional spread (deg)
                          !  -> crest length = wl/(2sin(wds*deg2rad))
#   elif defined SWASH
#    ifdef SWASH_GLOBEX_B2
#     define WAVE_MAKER_BICHROMATIC
        wf1=2*pi*0.42     ! GLOBEX B2
        wf2=2*pi*0.462
        wa1=0.09
        wa2=0.01
#    elif defined SWASH_GLOBEX_A3
#     define WAVE_MAKER_JONSWAP
        wp=2.25           ! period
        wa=0.0354         ! amplitude
        gamma=20.         ! JONSWAP peakedness parameter
#    else
        wp=2.
        wa=0.0442
#    endif
#   endif
        wf=2*pi/wp        ! frequency
!
!  Time & space origins
!
#   ifdef FLUME_WAVES
        x0=14.1
        time0=64.
#   else
        x0=0.
        y0=0.
        time0=0.
#   endif
!
!  Convert angles to rad
!
        wd=-wd*deg2rad    ! incidence angle 
        wds=wds*deg2rad   ! directional spread
!
!--------------------------------------------------------------------
!  Initialisation
!--------------------------------------------------------------------
!
#   ifdef FLUME_WAVES
!
!  Read file
!
        if (FIRST_TIME_STEP) then
!        open(117,file='datwaves_CORR1.txt',form='formatted',status='old')
          open(117,file='datwaves.txt',form='formatted',status='old')
          do k=1,320
            read(117,*) wa_bry(k), wf_bry(k),
     &                  wpha_bry(k), wk_bry(k)
            wa_bry(k)=wa_bry(k)*0.154/0.05
            !wpha_bry(k)=wpha_bry(k) + 1.5*pi
            khd=h(1,1)*wf_bry(k)**2/g   ! recompute wavenumber
            kh=sqrt( khd*khd+khd/(1.+khd*(K1+khd*(K2+khd*(K3+khd*(K4+
     &                                       khd*(K5+K6*khd)))))) )
            wk_bry(k)=kh/h(1,:1,$s/1)
          enddo
        endif
        ramp=tanh(dt/2.*float(iic-ntstart))

#   elif defined WAVE_MAKER_SPECTRUM
!
!  Build wave spectrum
!
        if (FIRST_TIME_STEP) then
          fmin=0.2*wf  ! frequency spread
          fmax=5.0*wf
          df=(fmax-fmin)/Nfrq
          cff2=0.
          do iw=1,Nfrq
            wf_bry(iw)=fmin+float(iw)*df
            khd=h(0,0)*wf_bry(iw)**2/g   ! wavenumber
            kh=sqrt( khd*khd+khd/(1.+khd*(K1+khd*(K2+khd*(K3+khd*(K4+
     &                                       khd*(K5+K6*khd)))))) )
            wk_bry(iw)=kh/h(0,0)
          enddo
#    ifdef WAVE_MAKER_JONSWAP
          do iw=1,Nfrq
            sigma=0.5*( 0.09*(1.+sign(1.,wf_bry(iw)-wf))+
     &                  0.07*(1.-sign(1.,wf_bry(iw)-wf)) )
            cff0=exp(-0.5*((wf_bry(iw)-wf)/(sigma*wf))**2)
            cff1=0.3119*(wf**4)/(wf_bry(iw)**5)
     &                 *exp(-1.25*(wf/wf_bry(iw))**4)*gamma**cff0
            cff2=16.*cff1*df  ! integral must be 1
            wa_bry(iw)=cff2
            cff3=cff3+cff2
          enddo
          do iw=1,Nfrq
            wa_bry(iw)=wa*sqrt(wa_bry(iw)/cff3) ! normalize
          enddo
#    elif defined WAVE_MAKER_GAUSSIAN
          cff2=0.
          do iw=1,Nfrq
            cff1=exp(-((wf_bry(iw)-wf)/0.1)**2)
            wa_bry(iw)=cff1
            cff2=cff2+cff1
          enddo
          do iw=1,Nfrq
            wa_bry(iw)=wa*sqrt(wa_bry(iw)/cff2) ! normalize
          enddo
#    endif
#    ifdef WAVE_MAKER_DSPREAD
          dmin=wd-25  ! directional spread
          dmax=wd+25
          dd=(dmax-dmin)/Ndir
          cff4=0.
          do jw=1,Ndir
            wd_bry(jw)=dmin+float(jw)*dd
            cff3=exp(-((wd_bry(jw)-wd)/max(wds,1.e-12))**2)
            wa_bry_d(jw)=cff3
            cff4=cff4+cff3
          enddo
          do jw=1,Ndir
            wa_bry_d(jw)=sqrt(wa_bry_d(jw)/cff4) ! normalize
          enddo
          call RANDOM_NUMBER(wpha_bry)  ! random phase
          do iw=1,Nfrq
            do jw=1,Ndir
              wpha_bry(iw,jw)=wpha_bry(iw,jw)*2.*pi
            enddo
          enddo
#    else
          call RANDOM_NUMBER(wpha_bry)  ! random phase
          do iw=1,Nfrq
            wpha_bry(iw)=wpha_bry(iw)*2.*pi
          enddo
        endif ! FIRST_TIME_STEP
        ramp=tanh(dt/wp*float(iic-ntstart))
#    endif /* WAVE_MAKER_DSPREAD */

#   elif defined WAVE_MAKER_BICHROMATIC
!
!  Bichromatic waves
!
        h0=h(0,1)
        ramp=tanh(dt/2.*float(iic-ntstart))
        wa1=ramp*wa1
        wa2=ramp*wa2
        khd=h0*wf1**2/g   ! compute wavenumber
        wk1=sqrt( khd*khd+khd/(1.+khd*(K1+khd*(K2+khd*(K3+khd*(K4+
     &                                    khd*(K5+K6*khd)))))) )/h0
        khd=h0*wf2**2/g
        wk2=sqrt( khd*khd+khd/(1.+khd*(K1+khd*(K2+khd*(K3+khd*(K4+
     &                                    khd*(K5+K6*khd)))))) )/h0
#   else
!
!  Monochromatic waves
!
        ramp=tanh(dt/wp*float(iic-ntstart))
        wa=wa*ramp
        wf=2*pi/wp
        khd=h0*wf**2/g   ! compute wavenumber
        wk=sqrt( khd*khd+khd/(1.+khd*(K1+khd*(K2+khd*(K3+khd*(K4+
     &                                   khd*(K5+K6*khd)))))) )/h0
#   endif /* FLUME ... */
!
!--------------------------------------------------------------------
!  Sea level zetabry
!--------------------------------------------------------------------
!
#   ifdef Z_FRC_BRY
        do j=JstrR,JendR
          h0=h(0,j)
#    if defined WAVE_MAKER_SPECTRUM
          zetabry_west(j)=0.
          do iw=1,Nfrq   ! frequency spread
#     if WAVE_MAKER_DSPREAD
            do jw=1,Ndir ! directional spread
              theta=(xr(0,j)-x0)*wk_bry(iw)*cos(wd_bry(jw))
     &             +(yr(0,j)-y0)*wk_bry(iw)*sin(wd_bry(jw))
     &             -(time-time0)*wf_bry(iw) 
     &                        -wpha_bry(iw,jw)
              zetabry_west(j)=zetabry_west(j) + 
     &                        ramp*wa_bry(iw)*wa_bry_d(jw)
     &                        *cos(theta)
            enddo
#     else
            theta=(xr(0,j)-x0)*wk_bry(iw)*cos(wd)
     &           +(yr(0,j)-y0)*wk_bry(iw)*sin(wd)
     &           -(time-time0)*wf_bry(iw) 
     &                      -wpha_bry(iw)
            zetabry_west(j)=zetabry_west(j) + 
     &                      ramp*wa_bry(iw)*cos(theta)
#     endif
          enddo

#    elif defined WAVE_MAKER_BICHROMATIC
           cff1=tanh(wk1*h0)
           cff2=tanh(wk2*h0)
           cff3=wa1*wa1*wk1*(3.-cff1**2)/(4.*cff1**3)
           cff4=wa2*wa2*wk2*(3.-cff2**2)/(4.*cff2**3)
           zetabry_west(j)=wa1*cos(wf1*time)   ! GLOBEX B2
     &                    +wa2*cos(wf2*time)
     &                 +cff3*cos(2*wf1*time)
     &                 +cff4*cos(2*wf2*time)
#    else
          time0=0.
          khd=h0*wf**2/g   ! compute wavenumber
          wk=sqrt( khd*khd+khd/(1.+khd*(K1+khd*(K2+khd*(K3+khd*(K4+
     &                                     khd*(K5+K6*khd)))))) )/h0
          sigma=tanh(wk*h0)
          theta=(xr(0,j)-x0)*wk*cos(wd)*cos(wds)
     &         +(yr(0,j)-y0)*wk*sin(wd)*cos(wds)
     &         -(time-time0)*wf
          cff_spread=cos(-(xr(0,j)-x0)*sin(wd)*sin(wds)*wk
     &                   +(yr(0,j)-y0)*cos(wd)*sin(wds)*wk)
          zetabry_west(j)=( wa*cos(theta)
#     ifdef STOKES_WAVES
     &            +wk*wa*wa*(3.-sigma**2)/(4.*sigma**3)*cos(2.*theta)
#     endif
     &                        )*cff_spread
#    endif /* FLUME ... */
        enddo  ! j loop
#   endif /* Z_FRC_BRY */
!
!--------------------------------------------------------------------
!  XI velocity components ubry and ubarbry
!--------------------------------------------------------------------
!
#   ifdef M3_FRC_BRY
        do j=JstrR,JendR
          h0=0.5*(h(0,j)+h(1,j))
          !Du=0.5*(z_w(0,j,N)+z_w(1,j,N)-z_w(0,j,0)-z_w(1,j,0))
          Du=h0
#    if defined FLUME_WAVES || \
     (defined WAVE_MAKER_SPECTRUM && !defined WAVE_MAKER_OBLIQUE \
                                  && !defined WAVE_MAKER_DSPREAD)
          do k=1,N
            ubry_west(j,k)=0.
            Zu=Du+0.5*(z_r(0,j,k)+z_r(1,j,k))
            do iw=1,Nfrq
              ubry_west(j,k)=ubry_west(j,k)+
     &                               ramp*wa_bry(iw)*wf_bry(iw)*
     &                                          cosh(wk_bry(iw)*Zu)/
     &                                          sinh(wk_bry(iw)*Du)
     &               *cos((0.5*(xr(0,j)+xr(1,j))-x0)*wk_bry(iw)
     &                                 -(time-time0)*wf_bry(iw)
     &                                            -wpha_bry(iw))
            enddo
          enddo

#    elif defined WAVE_MAKER_SPECTRUM
          khd=h0*wf**2/g   ! compute mean wavenumber
          wk=sqrt( khd*khd+khd/(1.+khd*(K1+khd*(K2+khd*(K3+khd*(K4+
     &                                     khd*(K5+K6*khd)))))) )/h0
          do k=1,N
            Zu=Du+0.5*(z_r(0,j,k)+z_r(1,j,k))
            ubry_west(j,k)=zetabry_west(j)*wf*cos(wd)
     &                                    *cosh(wk*Zu)
     &                                    /sinh(wk*Du)
          enddo

#    elif defined WAVE_MAKER_BICHROMATIC
          cff1=wa1*cos(wf1*time)*wf1/sinh(wk1*Du)
          cff2=wa2*cos(wf2*time)*wf2/sinh(wk2*Du)
          cff3=0.75*wa1*wa1*wk1*wf1/(sinh(wk1*Du))**4
     &                             *cos(2*wf1*time)
          cff4=0.75*wa2*wa2*wk2*wf2/(sinh(wk2*Du))**4
     &                             *cos(2*wf2*time)
          do k=1,N
            Zu=Du+0.5*(z_r(0,j,k)+z_r(1,j,k))
            ubry_west(j,k)=cff1*cosh(wk1*Zu)
     &                    +cff2*cosh(wk2*Zu)
     &                    +cff3*cosh(2*wk1*Zu)
     &                    +cff4*cosh(2*wk2*Zu)
          enddo
          cff1=0.5*g*wa1*wa1*wk1/(wf1*Du)
          cff2=0.5*g*wa2*wa2*wk2/(wf2*Du)
          !cff1=0.5*wa1*wa1*wk1*wf1/(sinh(wk1*Du))**2
          !cff2=0.5*wa2*wa2*wk2*wf2/(sinh(wk2*Du))**2
          do k=1,N
            ubry_west(j,k)=ubry_west(j,k)   ! compensation flow
     &                       -cff1 !*cosh(2*wk1*Zu)
     &                       -cff2 !*cosh(2*wk2*Zu)
          enddo
#    else
          xu=0.5*(xr(0,j)+xr(1,j))-x0
          yu=0.5*(yr(0,j)+yr(1,j))-y0
          theta=xu*cos(wd)*cos(wds)*wk
     &         +yu*sin(wd)*cos(wds)*wk
     &         -(time-time0)*wf
          cff_spread=cos(-xu*sin(wd)*sin(wds)*wk
     &                   +yu*cos(wd)*sin(wds)*wk)
          cff=wa*wf*cos(wd)/sinh(wk*Du)*cos(theta)*cff_spread
#     ifdef STOKES_WAVES
          cff2=cff*cos(2*theta)/cos(theta)*0.75*wk*wa/(sinh(wk*Du))**3
#     endif
          do k=1,N
            Zu=Du+0.5*(z_r(0,j,k)+z_r(1,j,k))
            ubry_west(j,k)=cff*cosh(wk*Zu)
#     ifdef STOKES_WAVES
     &                    +cff2*cosh(2*wk*Zu)
#     endif
          enddo
          cff1=0.5*g*wa*wa*wk/(wf*Du)
          do k=1,N
            ubry_west(j,k)=ubry_west(j,k)   ! compensation flow
     &                    -cff1
          enddo
#    endif /* FLUME_WAVES ... */
        enddo  ! j loop
#   endif /* M3_FRC_BRY */

#   ifdef M2_FRC_BRY
        do j=JstrR,JendR
          cff4=0.
          cff5=0.
          do k=1,N
            cff4=cff4+ubry_west(j,k)*(Hz(0,j,k)+Hz(1,j,k))
            cff5=cff5+(Hz(0,j,k)+Hz(1,j,k))
          enddo
          ubarbry_west(j)=cff4/cff5
        enddo
#   endif /* M2_FRC_BRY */
!
!--------------------------------------------------------------------
!  ETA velocity components vbry and vbarbry
!--------------------------------------------------------------------
!
#   ifdef M3_FRC_BRY
        do j=JstrV,JendR
          h0=0.5*(h(0,j)+h(0,j-1))
          !Dv=0.5*(z_w(0,j,N)+z_w(0,j-1,N)-z_w(0,j,0)-z_w(0,j-1,0))
          Dv=h0
#    if defined WAVE_MAKER_SPECTRUM && WAVE_MAKER_OBLIQUE
          khd=h0*wf**2/g   ! compute mean wavenumber
          wk=sqrt( khd*khd+khd/(1.+khd*(K1+khd*(K2+khd*(K3+khd*(K4+
     &                                     khd*(K5+K6*khd)))))) )/h0
          do k=1,N
            Zv=Dv+0.5*(z_r(0,j,k)+z_r(0,j-1,k))
            vbry_west(j,k)=zetabry_west(j)*wf*sin(wd)
     &                                    *cosh(wk*Zv)
     &                                    /sinh(wk*Dv)
          enddo

#    elif defined WAVE_MAKER_BICHROMATIC && WAVE_MAKER_OBLIQUE
          do k=1,N
              vbry_west(j,k)=0.
          enddo

#    elif defined WAVE_MAKER_OBLIQUE
          xv=0.5*(xr(0,j)+xr(0,j-1))-x0
          yv=0.5*(yr(0,j)+yr(0,j-1))-y0
          theta=xv*cos(wd)*cos(wds)*wk
     &         +yv*sin(wd)*cos(wds)*wk 
     &         -(time-time0)*wf
          cff_spread=cos(-xv*sin(wd)*sin(wds)*wk
     &                   +yv*cos(wd)*sin(wds)*wk)
          cff=wa*wf*sin(wd)/sinh(wk*Dv)*cos(theta)*cff_spread
#      ifdef STOKES_WAVES
          cff2=cff*cos(2*theta)/cos(theta)*0.75*wk*wa/(sinh(wk*Dv))**3
#      endif
          do k=1,N
            Zv=Dv+0.5*(z_r(0,j,k)+z_r(0,j-1,k))
            vbry_west(j,k)=cff*cosh(wk*Zv)
#      ifdef STOKES_WAVES
     &                    +cff2*cosh(2*wk*Zv)
#      endif
          enddo
#    else
          do k=1,N
              vbry_west(j,k)=0.
          enddo
#    endif /* WAVE_MAKER_SPECTRUM */
        enddo  ! j loop
#   endif /* M3_FRC_BRY */

#   ifdef M2_FRC_BRY
#    ifdef WAVE_MAKER_OBLIQUE
        do j=JstrV,JendR
          cff4=0.
          cff5=0.
          do k=1,N
            cff4=cff4+vbry_west(j,k)*(Hz(0,j,k)+Hz(0,j-1,k))
            cff5=cff5+(Hz(0,j,k)+Hz(0,j-1,k))
          enddo
          vbarbry_west(j)=cff4/cff5
        enddo
#    else
        do j=JstrV,JendR
          vbarbry_west(j)=0.
        enddo
#    endif
#   endif /* M2_FRC_BRY */
!
!--------------------------------------------------------------------
!  Z velocity component wbry
!--------------------------------------------------------------------
!
#   ifdef W_FRC_BRY
        do j=JstrR,JendR
          !Dr=z_w(0,j,N)-z_w(0,j,0)
          Dr=h(0,j)
#    ifdef WAVE_MAKER_SPECTRUM 
          do k=1,N
            wbry_west(j,k)=0.
            Zr=Dr+z_w(0,j,k)
            do iw=1,Nfrq
              wbry_west(j,k)=wbry_west(j,k)+
     &                               ramp*wa_bry(iw)*wf_bry(iw)*
     &                                          sinh(wk_bry(iw)*Zr)/
     &                                          sinh(wk_bry(iw)*Dr)
     &                             *sin((xr(0,j)-x0)*wk_bry(iw)
     &                                 -(time-time0)*wf_bry(iw)
     &                                            -wpha_bry(iw))
            enddo
          enddo

#    elif defined WAVE_MAKER_BICHROMATIC
          cff1=wa1*cos(wf1*time)*wf1/sinh(wk1*Dr)
          cff2=wa2*cos(wf2*time)*wf2/sinh(wk2*Dr)
          cff3=0.75*wa1*wa1*wk1*wf1/(sinh(wk1*Dr))**4
     &                               *cos(2*wf1*time)
          cff4=0.75*wa1*wa1*wk1*wf1/(sinh(wk2*Dr))**4
     &                               *cos(2*wf2*time)
          do k=0,N   
            Zr=Dr+z_w(0,j,k)     
            wbry_west(j,k)=cff1*sinh(wk1*Zr)
     &                    +cff2*sinh(wk1*Zr)
     &                  +cff3*sinh(2*wk1*Zr)
     &                  +cff4*sinh(2*wk1*Zr)
            wbry_west(j,k)=0.
          enddo
          do k=1,N-1
            wnbqbry_west(j,k)=0.5*wbry_west(j,k)*(Hz(0,j,k)+Hz(0,j,k+1))
          enddo
          wnbqbry_west(j,0)=0.5*wbry_west(j,0)*Hz(0,j,1)
          wnbqbry_west(j,N)=0.5*wbry_west(j,N)*Hz(0,j,N)

#    else
          theta=(xr(0,j)-x0)*cos(wd)*cos(wds)*wk
     &         +(yr(0,j)-y0)*sin(wd)*cos(wds)*wk - time*wf
          cff=wa*wf/sinh(wk*Dr)*cos(theta)
          cff2=0.75*cff*cos(2*theta)/cos(theta)*wk*wa/(sinh(wk*Dr))**3
          do k=1,N
            Zr=Dr+z_w(0,j,k)
            wbry_west(j,k)=cff*sinh(wk*Zr)
     &                 +cff2*sinh(2*wk*Zr)
          enddo
#    endif /* WAVE_MAKER_SPECTRUM ... */
        enddo  ! j loop
#   endif /* W_FRC_BRY */
!
!--------------------------------------------------------------------
!  TRACERS tbry
!--------------------------------------------------------------------
!
#   ifdef T_FRC_BRY
        if (FIRST_TIME_STEP) then
          do k=1,N
            do j=JstrR,JendR
              do itrc=1,NT
                tbry_west(j,k,itrc)=t(0,j,k,1,itrc)
              enddo
            enddo
          enddo
        endif
#   endif

