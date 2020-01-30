! $Id: ana_initial.F 1432 2014-01-20 15:24:14Z rblod $
!
!======================================================================
! ROMS_AGRIF is a branch of ROMS developped at IRD and INRIA, in France
! The two other branches from UCLA (Shchepetkin et al) 
! and Rutgers University (Arango et al) are under MIT/X style license.
! ROMS_AGRIF specific routines (nesting) are under CeCILL-C license.
! 
! ROMS_AGRIF website : http://roms.mpl.ird.fr
!======================================================================
!
#include "cppdefs.h"
#ifdef ANA_INITIAL
      subroutine ana_initial (tile)
      implicit none
# include "param.h"
       integer tile, trd
C$    integer omp_get_thread_num
# include "compute_tile_bounds.h"
      call ana_initial_tile (Istr,Iend,Jstr,Jend)
      return
      end
!*********************************************************************
      subroutine ana_initial_tile (Istr,Iend,Jstr,Jend)
!*********************************************************************
!
!--------------------------------------------------------------------
!  Set initial conditions for momentum and tracer variables using
!  analytical expressions.
!--------------------------------------------------------------------
!
# ifdef AGRIF
      use Agrif_Util
# endif
      use params
      implicit none
# include "param.h"
# include "grid.h"
# include "ocean2d.h"
# include "ocean3d.h"
# include "scalars.h"
      integer Istr,Iend,Jstr,Jend, i,j,k, itrc
# ifdef SOLITON
      real fac, x,y, x0,y0, cff1,cff2,cff3
# endif
# ifdef BASIN
      real cff1,cff2
# endif
# ifdef INTERNAL
      real cff
# endif
# ifdef BIOLOGY
      real temp, SiO4
# endif
# ifdef THACKER
      real ETA0, omega, depth, Lt
# endif
! aponte itide (
# ifdef ITIDE
      real lf, lk, ll, lsc, Ap_c, Ap_s
      integer nmode
!      integer nmode_ptb
!      real Acos_ptb, Asin_ptb
# endif
! aponte itide )

# include "compute_auxiliary_bounds.h"
!
# define IR_RANGE IstrR,IendR
# define IU_RANGE  Istr,IendR
# define JR_RANGE JstrR,JendR
# define JV_RANGE  Jstr,JendR
!
!  Initial conditions for free surface and 2D momentum components.
!------------------------------------------------------------------
!
# if defined BASIN      || defined CANYON_A   || defined CANYON_B  \
   || defined EQUATOR   || defined GRAV_ADJ   || defined OVERFLOW  \
   || defined SEAMOUNT  || defined SHELFRONT  || defined UPWELLING \
   || defined RIVER     || defined INTERNAL   || defined REGIONAL  \
   || defined JET       || defined INNERSHELF || defined SHOREFACE \
   || defined JETN || defined JETR || defined BTCHANNEL \
   || defined RIP
      do j=JR_RANGE
        do i=IR_RANGE
          zeta(i,j,1)=0.
          zeta(i,j,2)=0.
          ubar(i,j,1)=0.
          vbar(i,j,1)=0.
          ubar(i,j,2)=0.
          vbar(i,j,2)=0.
        enddo
      enddo
# elif defined SOLITON
#  ifdef AGRIF
      if (Agrif_Root()) then
#  endif
        x0=2.*xl/3.
        y0=el/2.
#  ifdef AGRIF
      else
        x0=2.*Agrif_Parent_xl/3.
        y0=Agrif_Parent_el/2.
      endif
#  endif
      cff1=0.395
      cff2=0.771*(cff1*cff1)
      do j=JR_RANGE
        do i=IR_RANGE
          x=xr(i,j)-x0
          y=yr(i,j)-y0
          cff3=exp(-cff1*x)
          fac=cff2*(2.*cff3/(1.+cff3*cff3))**2 
          zeta(i,j,1)=0.25*fac*(6.*y*y+3.)*exp(-0.5*y*y)
        enddo
      enddo
      do j=JR_RANGE
        do i=IU_RANGE
          x=0.5*(xr(i-1,j)+xr(i,j))-x0
          y=0.5*(yr(i-1,j)+yr(i,j))-y0
          cff3=exp(-cff1*x)
          fac=cff2 * (2.*cff3/(1.+cff3*cff3))**2
          ubar(i,j,1)=0.25*fac*(6.*y*y-9.)*exp(-0.5*y*y)
        enddo
      enddo
      do j=JV_RANGE
        do i=IR_RANGE
          x=0.5*(xr(i,j-1)+xr(i,j))-x0
          y=0.5*(yr(i,j-1)+yr(i,j))-y0
          cff3=exp(-cff1*x)
          fac=cff2 * (2.*cff3/(1.+cff3*cff3))**2
          vbar(i,j,1)=2.*fac*y*(-2.*cff1*tanh(cff1*x))
     &                                  *exp(-0.5*y*y)
        enddo
      enddo
# elif defined THACKER
      ETA0=0.1
      depth=10.
      Lt=80.e3
      omega=0.5*f(1,1) + sqrt( 0.25*f(1,1)**2 + 2*g*depth/Lt**2 )
      do j=JR_RANGE
        do i=IR_RANGE
          zeta(i,j,1)=2*ETA0*depth*(xr(i,j)/Lt - .5*ETA0/Lt)
          zeta(i,j,2)=zeta(i,j,1)
          ubar(i,j,1)=0.
          vbar(i,j,1)=-ETA0*omega*Lt
          ubar(i,j,2)=0
          vbar(i,j,2)=vbar(i,j,1)
        enddo
      enddo
# else
      ERROR: ANA_INITIAL: no initials for ubar,vbar,zeta
# endif

# ifdef SOLVE3D
!
!  Initial conditions for momentum components [m/s].
!--------------------------------------------------------------------
!
#  if defined BASIN     || defined CANYON_A   || defined CANYON_B  \
   || defined EQUATOR   || defined GRAV_ADJ   || defined OVERFLOW  \
   || defined SEAMOUNT  || defined SHELFRONT  || defined UPWELLING \
   || defined RIVER     || defined INTERNAL   || defined REGIONAL  \
   || defined JET       || defined INNERSHELF || defined SHOREFACE \
   || defined JETN || defined JETR || defined BTCHANNEL \
   || defined RIP
      do k=1,N
        do j=JR_RANGE
          do i=IR_RANGE
            u(i,j,k,1)=0.
            u(i,j,k,2)=0.
            v(i,j,k,1)=0.
            v(i,j,k,2)=0.
          enddo
        enddo
      enddo
#  elif defined THACKER
      do k=1,N
        do j=JR_RANGE
          do i=IR_RANGE
            u(i,j,k,1)=0.
            v(i,j,k,1)=-ETA0*omega*Lt
            u(i,j,k,2)=0.
            v(i,j,k,2)=v(i,j,k,1)
          enddo
        enddo
      enddo
#  else
        ERROR: ANA_INITIAL: no initial conditions for u,v.
#  endif
!
!  Initial conditions for tracer type variables.
!--------------------------------------------------------------------
!  Set initial conditions for potential temperature [degC] and
!  salinity [PSU].
!
#  ifdef BASIN
      cff1=(44.690/39.382)**2
      cff2=cff1*(rho0*800./g)*(5.0e-5/((42.689/44.690)**2))
      do k=1,N
        do j=JR_RANGE
          do i=IR_RANGE
            t(i,j,k,1,itemp)=cff2*exp(z_r(i,j,k)/800.)
     &                 *(0.6-0.4*tanh(z_r(i,j,k)/800.))
            t(i,j,k,2,itemp)=t(i,j,k,1,itemp)
          enddo
        enddo
      enddo
#  elif defined CANYON_A
      do k=1,N
        do j=JR_RANGE
          do i=IR_RANGE
            t(i,j,k,1,itemp)=T0
            t(i,j,k,2,itemp)=t(i,j,k,1,itemp)
          enddo
        enddo
      enddo
#  elif defined CANYON_B
      do k=1,N
        do j=JR_RANGE
          do i=IR_RANGE
            t(i,j,k,1,itemp)=3.488*exp(z_r(i,j,k)/800.)
     &        *(1.-0.666666666666*tanh(z_r(i,j,k)/800.))
            t(i,j,k,2,itemp)=t(i,j,k,1,itemp)
          enddo
        enddo
      enddo
#  elif defined SHOREFACE || defined RIP
      do k=1,N
        do j=jstrR,jendR
          do i=istrR,iendR
            t(i,j,k,1,itemp)=18.0
            t(i,j,k,2,itemp)=t(i,j,k,1,itemp)
#   ifdef SALINITY
            t(i,j,k,1,isalt)=32.0
            t(i,j,k,2,isalt)=t(i,j,k,1,isalt)
#   endif
          enddo
        enddo
      enddo
#   ifdef ANA_BRY
#    if defined OBC_WEST && (defined M3_FRC_BRY || defined T_FRC_BRY)
      if (WESTERN_EDGE) then
        do k=1,N
          do j=jstrR,jendR
#     ifdef M3_FRC_BRY
            u_west(j,k)=u(istr  ,j,k,1)
            v_west(j,k)=v(istr-1,j,k,1)
#     endif
#     ifdef T_FRC_BRY
            t_west(j,k,itemp)=t(istr-1,j,k,1,itemp)
#      ifdef SALINITY
            t_west(j,k,isalt)=t(istr-1,j,k,1,isalt)
#      endif
#     endif
          enddo
        enddo
      endif
#    endif
#   endif   /* ANA_BRY */
#  elif defined THACKER
      do k=1,N
        do j=JR_RANGE
          do i=IR_RANGE
            t(i,j,k,1,itemp)=20.
!            t(i,j,k,1,itemp)=20.-(zeta(i,j,1)-z_r(i,j,k))
            t(i,j,k,2,itemp)=t(i,j,k,1,itemp)
          enddo
        enddo
      enddo

#  elif defined EQUATOR
      do k=1,N
        do j=JstrR,JendR
          do i=IstrR,IendR
            t(i,j,k,1,itemp)=10.0
            t(i,j,k,2,itemp)=t(i,j,k,1,itemp)
#    ifdef SALINITY
            t(i,j,k,1,isalt)=35.0
            t(i,j,k,2,isalt)=t(i,j,k,1,isalt)
#    endif
          enddo
        enddo
      enddo
#  elif defined GRAV_ADJ
      do k=1,N
        do j=JstrR,JendR
          do i=IstrR,IendR
#   ifdef MPI
           if (i+iminmpi-1 .LE. LLm/2) then
#   else
           if (i .LE. LLm/2) then
#   endif            
            t(i,j,k,1,itemp)=5/Tcoef
            t(i,j,k,2,itemp)=5/Tcoef
            else
            t(i,j,k,1,itemp)=10/Tcoef
            t(i,j,k,2,itemp)=10/Tcoef           
            endif
          enddo
        enddo
      enddo
#  elif defined INNERSHELF
      do k=1,N
        do j=JR_RANGE
          do i=IR_RANGE
#   ifdef INNERSHELF_EKMAN
            t(i,j,k,1,itemp)=15.
#   else
#    define Z0 (-80.)
#    define THKNSS 50.
#    define Z1 0.
#    define STRAT 1000.
            t(i,j,k,1,itemp)=14.+5.*THKNSS*log(
     &                         cosh((z_w(i,j,k  )-Z0)/THKNSS)
     &                        /cosh((z_w(i,j,k-1)-Z0)/THKNSS)
     &                            )/(z_w(i,j,k)-z_w(i,j,k-1))
     &               +((z_w(i,j,k)+z_w(i,j,k-1))/2.-Z1)/STRAT
#    undef Z0
#    undef THKNSS
#    undef Z1
#    undef STRAT
#   endif
            t(i,j,k,2,itemp)=t(i,j,k,1,itemp)
#   ifdef SALINITY
            t(i,j,k,1,isalt)=35.
            t(i,j,k,2,isalt)=t(i,j,k,1,isalt)
#   endif
          enddo
        enddo
      enddo
#  elif defined  INTERNAL
#   define BVFREQ 2.e-3
      cff=BVFREQ*BVFREQ/g
      do k=1,N
        do j=JR_RANGE
          do i=IR_RANGE
            t(i,j,k,1,itemp)=rho0*(1-cff*z_r(i,j,k))
            t(i,j,k,2,itemp)=t(i,j,k,1,itemp)
          enddo
        enddo
      enddo
#  elif defined OVERFLOW
      do k=1,N
        do j=JR_RANGE
          do i=IR_RANGE

            t(i,j,k,1,itemp)=T0*(0.5-0.5*tanh( yr(i,j)/1000.-25.))

            t(i,j,k,2,itemp)=t(i,j,k,1,itemp)
          enddo
        enddo
      enddo
#  elif defined  RIVER
      do k=1,N
        do j=JR_RANGE
          do i=IR_RANGE
            t(i,j,k,1,itemp)=10.  !4.+10.*exp(z_r(i,j,k)/50.)
            t(i,j,k,1,isalt)=36.
            t(i,j,k,2,itemp)=t(i,j,k,1,itemp)
            t(i,j,k,2,isalt)=t(i,j,k,1,isalt)
          enddo
        enddo
      enddo
#  elif defined SEAMOUNT
      do k=1,N
        do j=JR_RANGE
          do i=IR_RANGE
            t(i,j,k,1,itemp)=5.*15.*exp(z_r(i,j,k)/100.)
            t(i,j,k,2,itemp)=t(i,j,k,1,itemp)
#   ifdef SALINITY
            t(i,j,k,1,isalt)=35.
            t(i,j,k,2,isalt)=t(i,j,k,1,isalt)
#   endif
          enddo
        enddo
      enddo
#  elif defined SHELFRONT
      do k=1,N
        do j=JR_RANGE
          do i=IR_RANGE
            t(i,j,k,1,itemp)=T0+2.5*tanh((yr(i,j)-50000.0)/20000.0)
            t(i,j,k,2,itemp)=t(i,j,k,1,itemp)
            t(i,j,k,1,isalt)=S0  +  tanh((yr(i,j)-50000.0)/20000.0)
            t(i,j,k,2,isalt)=t(i,j,k,1,isalt)
          enddo
        enddo
      enddo
#  elif  defined SED_TEST2 
      do k=1,N
        do j=JstrR,JendR
          do i=IstrR,IendR
            t(i,j,k,1,itemp)=20.0
            t(i,j,k,1,isalt)=0.0
            t(i,j,k,2,itemp)=t(i,j,k,1,itemp)
            t(i,j,k,2,isalt)=t(i,j,k,1,isalt)
          enddo
        enddo
      enddo
#  elif defined UPWELLING
      do k=1,N
        do j=JR_RANGE
          do i=IR_RANGE
#   define Z0 (-35.)
#   define THKNSS 6.5
#   define Z1 (-75.)
#   define STRAT 150.
           t(i,j,k,1,itemp)=14.+4.*THKNSS*log(
     &                          cosh((z_w(i,j,k )-Z0)/THKNSS)
     &                         /cosh((z_w(i,j,k-1)-Z0)/THKNSS)
     &                             )/(z_w(i,j,k)-z_w(i,j,k-1))
     &               +((z_w(i,j,k)+z_w(i,j,k-1))/2.-Z1)/STRAT
#   undef Z0
#   undef THKNSS
#   undef Z1
#   undef STRAT
            t(i,j,k,2,itemp)=t(i,j,k,1,itemp)
#   ifdef SALINITY
            t(i,j,k,1,isalt)=35.
            t(i,j,k,2,isalt)=t(i,j,k,1,isalt)
#   endif
          enddo
        enddo
      enddo
#  elif defined REGIONAL
      do k=1,N
        do j=JR_RANGE
          do i=IR_RANGE
#   define Z0 (-80.)
#   define THKNSS 50.
#   define Z1 0.
#   define STRAT 1000.
           t(i,j,k,1,itemp)=14.+5.*THKNSS*log(
     &                            cosh((z_w(i,j,k )-Z0)/THKNSS)
     &                           /cosh((z_w(i,j,k-1)-Z0)/THKNSS)
     &                               )/(z_w(i,j,k)-z_w(i,j,k-1))
     &               +((z_w(i,j,k)+z_w(i,j,k-1))/2.-Z1)/STRAT
#   undef Z0
#   undef THKNSS
#   undef Z1
#   undef STRAT
            t(i,j,k,2,itemp)=t(i,j,k,1,itemp)
#   ifdef SALINITY
            t(i,j,k,1,isalt)=35.
            t(i,j,k,2,isalt)=t(i,j,k,1,isalt)
#   endif
#   ifdef BIOLOGY
            temp=t(i,j,k,1,itemp)
            if (temp.lt.8.) then
               SiO4=30.
            elseif (temp.ge.8. .and. temp.le.11.) then
               SiO4=30.-((temp-8.)*(20./3.))
            elseif (temp.gt.11. .and. temp.le.13.) then
               SiO4=10.-((tdy=Length_ETA/float(MMm) emp-11.)*(8./2.))
            elseif (temp.gt.13. .and. temp.le.16.) then
               SiO4=2.-((temp-13.)*(2./3.))
            elseif (temp.gt.16.) then
              SiO4=0.
            endif
            t(i,j,k,1,iNO3_)=1.67+0.5873*SiO4+0.0144*SiO4**2
     &                               +0.0003099*SiO4**3
            t(i,j,k,2,iNO3_)=t(i,j,k,1,iNO3_)
#  ifdef PISCES
            t(i,j,k,1,iDIC_)=2150.
            t(i,j,k,2,iDIC_)=t(i,j,k,1,iDIC_)
            t(i,j,k,1,iTAL_)=2350.
            t(i,j,k,2,iTAL_)=t(i,j,k,1,iTAL_)
            t(i,j,k,1,iOXY_)=200.
            t(i,j,k,2,iOXY_)=t(i,j,k,1,iOXY_)
            t(i,j,k,1,iCAL_)=0.01
            t(i,j,k,2,iCAL_)=t(i,j,k,1,iCAL_)
            t(i,j,k,1,iPO4_)=t(i,j,k,1,iNO3_)/16.
            t(i,j,k,2,iPO4_)=t(i,j,k,1,iPO4_)
            t(i,j,k,1,iPOC_)=0.01
            t(i,j,k,2,iPOC_)=t(i,j,k,1,iPOC_)
            t(i,j,k,1,iSIL_)=91.51
            t(i,j,k,2,iSIL_)=t(i,j,k,1,iSIL_)
            t(i,j,k,1,iPHY_)=0.01
            t(i,j,k,2,iPHY_)=t(i,j,k,1,iPHY_)
            t(i,j,k,1,iZOO_)=0.01
            t(i,j,k,2,iZOO_)=t(i,j,k,1,iZOO_)
            t(i,j,k,1,iDOC_)=5.
            t(i,j,k,2,iDOC_)=t(i,j,k,1,iDOC_)
            t(i,j,k,1,iDIA_)=0.01
            t(i,j,k,2,iDIA_)=t(i,j,k,1,iDIA_)
            t(i,j,k,1,iMES_)=0.01
            t(i,j,k,2,iMES_)=t(i,j,k,1,iMES_)
            t(i,j,k,1,iBSI_)=1.5e-3
            t(i,j,k,2,iBSI_)=t(i,j,k,1, zeta(i,j,1)iBSI_)
            t(i,j,k,1,iFER_)=6.e-4
            t(i,j,k,2,iFER_)=t(i,j,k,1,iFER_)
            t(i,j,k,1,iBFE_)=1.E-2*5E-6
            t(i,j,k,2,iBFE_)=t(i,j,k,1,iBFE_)
            t(i,j,k,1,iGOC_)=0.01
            t(i,j,k,2,iGOC_)=t(i,j,k,1,iGOC_)
            t(i,j,k,1,iSFE_)=0.01*5.E-6
            t(i,j,k,2,iSFE_)=t(i,j,k,1,iSFE_)
            t(i,j,k,1,iDFE_)=0.01*5.E-6
            t(i,j,k,2,iDFE_)=t(i,j,k,1,iDFE_)
            t(i,j,k,1,iDSI_)=0.01*0.15
            t(i,j,k,2,iDSI_)=t(i,j,k,1,iDSI_)
            t(i,j,k,1,iNFE_)=0.01*5.E-6
            t(i,j,k,2,iNFE_)=t(i,j,k,1,iNFE_)
            t(i,j,k,1,iNCH_)=0.01*12./55.
            t(i,j,k,2,iNCH_)=t(i,j,k,1,iNCH_)
            t(i,j,k,1,iDCH_)=0.01*12./55.
            t(i,j,k,2,iDCH_)=t(i,j,k,1,iDCH_)
            t(i,j,k,1,iNH4_)=0.01
            t(i,j,k,2,iNH4_)=t(i,j,k,1,iNH4_)
#  elif defined BIO_NChlPZD
            t(i,j,k,1,iChla)=0.08
            t(i,j,k,2,iChla)=t(i,j,k,1,iChla)
            t(i,j,k,1,iPhy1)=0.1
            t(i,j,k,2,iPhy1)=t(i,j,k,1,iPhy1)
            t(i,j,k,1,iZoo1)=0.06
            t(i,j,k,2,iZoo1)=t(i,j,k,1,iZoo1)
            t(i,j,k,1,iDet1)=0.02
            t(i,j,k,2,iDet1)=t(i,j,k,1,iDet1)
#     ifdef OXYGEN
            t(i,j,k,1,iO2)=250.0
            t(i,j,k,2,iO2)=t(i,j,k,1,iO2)
#     endif
#    elif defined BIO_N2ChlPZD2
            t(i,j,k,1,iNH4_)=0.1
            t(i,j,k,2,iNH4_)=t(i,j,k,1,iNH4_)
            t(i,j,k,1,iPhy1)=0.06
            t(i,j,k,2,iPhy1)=t(i,j,k,1,iPhy1)
            t(i,j,k,1,iZoo1)=0.04
            t(i,j,k,2,iZoo1)=t(i,j,k,1,iZoo1)
            t(i,j,k,1,iDet1)=0.02
            t(i,j,k,2,iDet1)=t(i,j,k,1,iDet1)
            t(i,j,k,1,iDet2)=0.02
            t(i,j,k,2,iDet2)=t(i,j,k,1,iDet2)
#    elif defined BIO_BioEBUS
            t(i,j,k,1,iNO2_)=0.05*exp(z_r(i,j,k)/100.)
            t(i,j,k,2,iNO2_)=t(i,j,k,1,iNO2_)
            t(i,j,k,1,iNH4_)=0.1*exp(z_r(i,j,k)/100.)
            t(i,j,k,2,iNH4_)=t(i,j,k,1,iNH4_)
            t(i,j,k,1,iPhy1)=0.04*exp(z_r(i,j,k)/50.)
            t(i,j,k,2,iPhy1)=t(i,j,k,1,iPhy1)
            t(i,j,k,1,iPhy2)=0.06*exp(z_r(i,j,k)/50.)
            t(i,j,k,2,iPhy2)=t(i,j,k,1,iPhy2)
            t(i,j,k,1,iZoo1)=0.04*exp(z_r(i,j,k)/100.)
            t(i,j,k,2,iZoo1)=t(i,j,k,1,iZoo1)
            t(i,j,k,1,iZoo2)=0.04*exp(z_r(i,j,k)/100.)
            t(i,j,k,2,iZoo2)=t(i,j,k,1,iZoo2)
            t(i,j,k,1,iDet1)=0.02
            t(i,j,k,2,iDet1)=t(i,j,k,1,iDet1)
            t(i,j,k,1,iDet2)=0.02
            t(i,j,k,2,iDet2)=t(i,j,k,1,iDet2)
            t(i,j,k,1,iDON) =0.5*exp(z_r(i,j,k)/100.)    
            t(i,j,k,2,iDON) =t(i,j,k,1,iDON)
            t(i,j,k,2,iO2) = 250.*exp(z_r(i,j,k)/100.)
            t(i,j,k,2,iO2) =t(i,j,k,1,iO2)
#     ifdef NITROUS_OXIDE
            t(i,j,k,1,iN2O)=-(0.008*exp(z_r(i,j,k)/100.)-0.01)
            if (t(i,j,k,1,iN2O).lt.0.) then
              t(i,j,k,1,iN2O)=0.
            endif
            t(i,j,k,2,iN2O)=t(i,j,k,1,iN2O)  
#     endif
#    endif
#   endif /* BIOLOGY */
          enddo
        enddo
      enddo

#  elif defined JET
!
! Compute analytical t and geostrophically adjusted zeta,ub,vb,u,v fields
! with no zonal perturbation first
!  --> store in climatological arrays
!
      call ana_jet_tile (Istr,Iend,Jstr,Jend,0) 
!
! Compute same fields with zonal perturbation 
! --> store in initial arrays
!
      call ana_jet_tile (Istr,Iend,Jstr,Jend,1)



! aponte jet (
#  elif defined JETN

! same but with different parameters

! aponte density (
# ifdef ANA_DENSITY_FILE

      call ana_get_density (Istr,Iend,Jstr,Jend)
      call ana_jetn_tile (Istr,Iend,Jstr,Jend,0)
      call ana_jetn_tile (Istr,Iend,Jstr,Jend,1)

# else   /* ANA_DENSITY_FILE  */

      call ana_jetn_tile (Istr,Iend,Jstr,Jend,0)
      call vmodes(Istr,Iend,Jstr,Jend)
      call ana_jetn_tile (Istr,Iend,Jstr,Jend,1)

      
# endif /* ANA_DENSITY_FILE  */
! aponte density )
! aponte itide (
! aponte vmodes (
! load the structure of vertical modes
!       call get_vmodes
! compute vertical modes
!      call vmodes(Istr,Iend,Jstr,Jend)
! aponte vmodes )
! aponte itide

! aponte jet )
#    if defined ITIDE
! tidal forcing at the southern boundary
! will eventually go into the input file
!       omega_itide = 1.402D-4
!       nmode_itide = 2
!      start_time_itide = 0.0D0


! new formulation
!       lf = f(Istr+1,Jstr+1)
!       ll = sqrt(omega_itide**2-lf**2)/c_vmodes(nmode_itide)
       Arc_itide = 0.0D+0
       Ars_itide = 0.0D+0
       Auc_itide = 0.0D+0
       Aus_itide = 0.0D+0
! scaling follows reponse to a dirac
! should use speeds that include rotation
! should produce a 1 cm tidal sea level amplitude for mode 1
       Avc_itide = Az_itide/dely_itide
     &       *c_vmodes(nmode_itide)**2/c_vmodes(2)
     &       *2.27D+0

# ifdef ITIDEP
! alternative: mulitple amplitudes can be generated
!!
!! truth :
!!
! barotropic mode:
       nmode = 1
       Atide_cos_mode(nmode) = 0.0
       Atide_sin_mode(nmode) = 0.0
       dAtide_cos_mode(nmode) = 0.0
       dAtide_sin_mode(nmode) = 0.0
! first baroclinic mode
       nmode = 2
       Atide_cos_mode(nmode) = Az_itide/dely_itide
     &       *c_vmodes(nmode)**2/c_vmodes(2)
     &       *2.27D+0
       Atide_sin_mode(nmode) = 0.0
       dAtide_cos_mode(nmode) = 0.0
       dAtide_sin_mode(nmode) = 0.0
!       MPI_master_only write( stdout,'(3x,A,I2,A,1PE15.8)') 
!     &    'Truth Atide_cos_mode(',nmode,') = ',
!     &    Atide_cos_mode(nmode)
!       MPI_master_only write( stdout,'(3x,A,I2,A,1PE15.8)')
!     &    'Truth Atide_sin_mode(',nmode,') = ',
!     &    Atide_sin_mode(nmode)
! higher modes
       do nmode=3,10
        Atide_cos_mode(nmode) = 0.0 
        Atide_sin_mode(nmode) = 0.0
        dAtide_cos_mode(nmode) = 0.0
        dAtide_sin_mode(nmode) = 0.0
!        MPI_master_only write( stdout,'(3x,A,I2,A,1PE15.8)')
!     &    'Truth Atide_cos_mode(',nmode,') = ',
!     &    Atide_cos_mode(nmode)
!        MPI_master_only write( stdout,'(3x,A,I2,A,1PE15.8)')
!     &    'Truth Atide_sin_mode(',nmode,') = ',
!     &    Atide_sin_mode(nmode)
       enddo
       call params_output('itideptb_truth_Acos_mode', Atide_cos_mode)
       call params_output('itideptb_truth_Asin_mode', Atide_sin_mode)
!!
!! first guess
!!
       if ( flag_ptb .ge. 1 ) then
! first guess is 0 for now
        do nmode=1,10
         Atide_cos_mode(nmode) = 0.0 
         Atide_sin_mode(nmode) = 0.0 
!         MPI_master_only write( stdout,'(3x,A,I2,A,1PE15.8)')
!     &    'First guess Atide_cos_mode(',nmode,') = ',
!     &    Atide_cos_mode(nmode)
!         MPI_master_only write( stdout,'(3x,A,I2,A,1PE15.8)')
!     &    'First guess Atide_sin_mode(',nmode,') = ',
!     &    Atide_sin_mode(nmode)
        enddo
        call params_output('itideptb_fguess_Acos_mode', Atide_cos_mode)
        call params_output('itideptb_fguess_Asin_mode', Atide_sin_mode)        
       endif
!!
!! perturbation on first guess
!! 
       if ( flag_ptb .ge. 2 ) then
!       nmode_ptb=4
!        MPI_master_only write(*,*) 'nmode_ptb = ', nmode_ptb
!        MPI_master_only write( stdout,'(3x,A,I2,A,1PE15.8)')
!     &    'Pertubation dcos(',nmode_ptb,') = ', dcos_ptb
!         MPI_master_only write( stdout,'(3x,A,I2,A,1PE15.8)')
!     &    'Perturbation dsin_ptb(',nmode_ptb,') = ', dsin_ptb
!        Atide_cos_mode(nmode_ptb) = Atide_cos_mode(nmode_ptb) + dcos_ptb
!        Atide_sin_mode(nmode_ptb) = Atide_sin_mode(nmode_ptb) + dsin_ptb
!        MPI_master_only write( stdout,'(3x,A,I2,A,1PE15.8)')
!     &    'Atide_cos_mode(',nmode_ptb,') = ', Atide_cos_mode(nmode_ptb)
!         MPI_master_only write( stdout,'(3x,A,I2,A,1PE15.8)')
!     &    'Atide_sin_mode(',nmode_ptb,') = ', Atide_sin_mode(nmode_ptb)
        call params_output('itideptb_nmode', nmode_ptb)
        call params_output('itideptb_dcos', dcos_ptb)
        call params_output('itideptb_dsin', dsin_ptb)
        call params_output('itideptb_full_Acos_mode', Atide_cos_mode)
        call params_output('itideptb_full_Asin_mode', Atide_sin_mode)        
       endif
! ....
# endif /* ITIDEP */

       Avs_itide = 0.0D+0
       Azc_itide = 0.0D+0
       Azs_itide = 0.0D+0
!       Arc_itide = 6.5372914D0
!       Ars_itide = 0.0D0
!       Auc_itide = 0.0D0
!       Aus_itide = -2.3844872D-02
!       Avc_itide = -3.1044374D-02
!       Avs_itide = 0.0D0
!       Azc_itide = 1.6314420D-02
!       Azs_itide = 0.0D0

! input is sea level displacement
! we divide it by nmode_itide^2 in order to maintain
! the Froude number approximately constant
!       Azc_itide = 1.0D-2
!       Azc_itide = 0.0D0
!       Azs_itide = 0.0D0
!       Azs_itide = Az_itide
!       Azs_itide = Az_itide/max(1.0,real(nmode_itide)-1.0)**2
!Azeta = Ap * phi(end,im) /grd.g/grd.rho0;
!       Ap_c = g*Azc_itide/phiuv(N,nmode_itide)
!       Ap_s = g*Azs_itide/phiuv(N,nmode_itide)
!Ap = -grd.g*Arho;
!       Arc_itide = -Ap_c*rho0/g
!       Ars_itide = -Ap_s*rho0/g
!% polarization relations
!Au = (omega*k-sqrt(-1)*f*l)/(omega^2-f^2) *Ap/grd.rho0;
!Av = (omega*l+sqrt(-1)*f*k)/(omega^2-f^2) *Ap/grd.rho0;
!
! from sage worksheet: https://cloud.sagemath.com/projects
! ucos = (omega*k*pcos-f*l*psin)/(omega^2-f^2)/grd.rho0;
! usin = (f*l*pcos+omega*k*psin)/(omega^2-f^2)/grd.rho0;
! vcos = (omega*l*pcos+f*k*psin)/(omega^2-f^2)/grd.rho0;
! vsin = (-f*k*pcos+omega*l*psin)/(omega^2-f^2)/grd.rho0;
!
!       lf = f(Istr+1,Jstr+1)
!       lk = 0.
!       ll = sqrt(omega_itide**2-lf**2)/c_vmodes(nmode_itide)
!       lsc = 1./(omega_itide**2-lf**2)
!
!       Auc_itide = (omega_itide*lk*Ap_c-lf*ll*Ap_s)*lsc
!       Aus_itide = (omega_itide*lk*Ap_s+lf*ll*Ap_c)*lsc
!       Avc_itide = (omega_itide*ll*Ap_c+lf*lk*Ap_s)*lsc
!       Avs_itide = (omega_itide*ll*Ap_s-lf*lk*Ap_c)*lsc
!
!   eta0 = 1.0D-2
!       Awc_itide = 0.0D0
!       Aws_itide = 0.0D0
!       Aws_itide = 1.0D-2*omega_itide*g/phiuv(N,nmode_itide) *1.0D0 
!       Aws_itide = Az_itide*omega_itide*g/phiuv(N,nmode_itide) 
!       MPI_master_only write(*,*) 'Awc_itide = ', Awc_itide
!       MPI_master_only write(*,*) 'Aws_itide = ', Aws_itide
!
!% real( A*exp(k*x+l*y-omega*t) ) = 
!%   Ar*cos()-Ai*sin()

       if ( .not. flag_itide ) then
        Arc_itide = 0.0D0
        Ars_itide = 0.0D0
        Auc_itide = 0.0D0
        Aus_itide = 0.0D0
        Avc_itide = 0.0D0
        Avs_itide = 0.0D0
        Azc_itide = 0.0D0
        Azs_itide = 0.0D0
        Awc_itide = 0.0D0
        Aws_itide = 0.0D0
       endif
       
       call params_output('itide_Ar_cos', Arc_itide)
       call params_output('itide_Ar_sin', Ars_itide)       
       call params_output('itide_Au_cos', Auc_itide)
       call params_output('itide_Au_sin', Aus_itide)       
       call params_output('itide_Av_cos', Avc_itide)
       call params_output('itide_Av_sin', Avs_itide)       
       call params_output('itide_Az_cos', Azc_itide)
       call params_output('itide_Az_sin', Azs_itide)       
       call params_output('itide_Aw_cos', Awc_itide)
       call params_output('itide_Aw_sin', Aws_itide)

#    endif
! aponte jetn )


! aponte jetr (
#  elif defined JETR

! same but with different parameters
      call ana_jetr_tile (Istr,Iend,Jstr,Jend,0)
      call ana_jetr_tile (Istr,Iend,Jstr,Jend,1)

! aponte jetr )

! aponte btc (
#   elif defined BTCHANNEL

#    if defined ITIDE

       start_time_itide = 0.D0
       Arc_itide = 0.D0
       Ars_itide = 0.D0
       Auc_itide = 0.D0
       Aus_itide = 0.D0
       Avc_itide = 0.D0
       Avs_itide = 0.D0
       Azc_itide = 1.0D-02
       Azs_itide = 0.0D0

       if ( .not. flag_itide ) then
        Auc_itide = 0.0D0
        Aus_itide = 0.0D0
        Avc_itide = 0.0D0
        Avs_itide = 0.0D0
        Azc_itide = 0.0D0
        Azs_itide = 0.0D0
       endif

#     endif

      do k=1,N
        do j=JR_RANGE
          do i=IR_RANGE
            t(i,j,k,1,itemp)=T0
            t(i,j,k,2,itemp)=T0
!            t(i,j,k,1,isalt)=???
!            t(i,j,k,2,isalt)=???
          enddo
        enddo
      enddo
! aponte btc )

#  else
      do k=1,N
        do j=JR_RANGE
          do i=IR_RANGE
            t(i,j,k,1,itemp)=???
            t(i,j,k,2,itemp)=???
            t(i,j,k,1,isalt)=???
            t(i,j,k,2,isalt)=???
          enddo
        enddo
      enddo
#  endif

#  if defined EW_PERIODIC || defined NS_PERIODIC || defined  MPI
      call exchange_u3d_tile (Istr,Iend,Jstr,Jend,
     &                        u(START_2D_ARRAY,1,1))
      call exchange_v3d_tile (Istr,Iend,Jstr,Jend,
     &                        v(START_2D_ARRAY,1,1))
      do itrc=1,NT
        call exchange_r3d_tile (Istr,Iend,Jstr,Jend,
     &                          t(START_2D_ARRAY,1,1,itrc))
      enddo
#  endif

# endif /* SOLVE3D */

# if defined EW_PERIODIC || defined NS_PERIODIC || defined  MPI
      call exchange_r2d_tile (Istr,Iend,Jstr,Jend,
     &                        zeta(START_2D_ARRAY,1))
      call exchange_u2d_tile (Istr,Iend,Jstr,Jend,
     &                        ubar(START_2D_ARRAY,1))
      call exchange_v2d_tile (Istr,Iend,Jstr,Jend,
     &                        vbar(START_2D_ARRAY,1))
# endif

# undef IR_RANGE
# undef IU_RANGE
# undef JR_RANGE
# undef JV_RANGE
      return
      end

# ifdef JET
!*********************************************************************
      subroutine ana_jet_tile (Istr,Iend,Jstr,Jend, jet_perturb)
!*********************************************************************
!
!--------------------------------------------------------------------
!  Set initial conditions for momentum and tracer variables using
!  analytical expressions for the JET test case
!--------------------------------------------------------------------
!
#  ifdef AGRIF
      use Agrif_Util
#  endif
      implicit none
#  include "param.h"
#  include "grid.h"
#  include "ocean2d.h"
#  include "ocean3d.h"
#  include "scalars.h"
#  include "mpi_cpl.h"
#  include "private_scratch.h"
#  ifdef CLIMATOLOGY
#   include "climat.h"
#  endif
      integer Istr,Iend,Jstr,Jend, i,j,k, itrc, 
     &        jet_perturb, tndx
      real cff,x,y,Fyz,
     &     rhomaxs,bgdrhos,zs1,dzs,drhos,drhosfs,
     &     rhomaxn,bgdrhon,zn1,dzn,drhon,drhosfn,
     &     drhosf,z0,z0p,dzs_a,dzn_a,zk,
     &     Ljet,cff_perturb,ierr,tanh
      real*QUAD mssh,area,cff1
! aponte quad (
      real mssh_dbl,area_dbl
! aponte quad )
      real rhoprof(0:N,2)
      real ru(PRIVATE_2D_SCRATCH_ARRAY,N),   
     &     rv(PRIVATE_2D_SCRATCH_ARRAY,N),   
     &      P(PRIVATE_2D_SCRATCH_ARRAY,N), 
     &     FC(PRIVATE_2D_SCRATCH_ARRAY),
     &    dZx(PRIVATE_2D_SCRATCH_ARRAY),
     &     rx(PRIVATE_2D_SCRATCH_ARRAY),
     &    dRx(PRIVATE_2D_SCRATCH_ARRAY)
#  ifdef MPI
      real*QUAD allmssh(1,NNODES),allarea(1,NNODES)
! aponte quad (
      real allmssh_dbl(1,NNODES),allarea_dbl(1,NNODES)
! aponte quad )
      include 'mpif.h'
#  endif
      real asymth,zz,zasym,dzasym
      asymth(zz,zasym,dzasym) = 
     &         (zz-zasym)
     &         *sqrt(1+0.5*((zz-zasym)+abs(zz-zasym))**2/dzasym**2)
     &         +zasym
!
#  include "compute_auxiliary_bounds.h"
!
! ---------------------------------------------------------------------
! BAROCLINIC JET TEST CASE
! ---------------------------------------------------------------------
!
! --- Jet and density profile parameters ---
!
      Ljet        = 1600.e3;  ! Jet width in km
!
! southern profile
      rhomaxs=27.75     ! max density
      bgdrhos=9.8e-6    ! background stratification
      zs1=-1000
      dzs=700
      drhos=1.4         ! gives strength of Phillips mode (0 to 1.4)
! northern profile 
      rhomaxn=27.7573   ! max density
      bgdrhon=9.8e-6    ! background stratification
      zn1=-400 
      dzn=300
! surface Charney mode 
      drhosfs=1.5       ! northern anomaly  (0 for Phillips mode only)
      drhosfn=0.0       ! southern anomaly  (0 for Phillips mode only)
      z0=-300.          ! vertical penetration (-600 for strong Charney mode)
      drhosf=0.00       ! weakens Charney mode (canceled for 1.05)
      z0p=-110.         ! only used for drhosf non zero
!
! zonal perturbation
!
      if (jet_perturb .eq. 0) then
        cff_perturb = 0.0    ! zonal perturbation coefficient
        tndx=2
      else
        cff_perturb = 0.02
        tndx=1
        zeta(:,:,tndx)=0.
        call set_depth_tile (Istr,Iend,Jstr,Jend)
      endif
!
! --- Build northern and southern density profiles ---
!
! First compute background density profiles associated with
! a gentle stratification that does not depend 
! on jet side. It is there to ensure static stability. 
      do k=1,N
        zk = z_r(istr,jstrR,k)
        rhoprof(k,1)=rhomaxs-bgdrhos*(zk+h(istr,jstrR))
        rhoprof(k,2)=rhomaxn-bgdrhon*(zk+h(istr,jstrR))
      enddo
!
! Second, get main north/south contrast with a distorded 
! tanh shape of variable amplitude.Distorsion of the tanh 
! is done with the asym function that increases 
! stratification in the upper ocean
      dzs_a=1.3*dzs
      do k=1,N
        zk=asymth(z_r(istr,jstrR,k),zs1,dzs_a)
        rhoprof(k,1)=rhoprof(k,1)-drhos*(0.5+0.5*tanh((zk-zs1)/dzs))
      enddo
      dzn_a=1.3*dzn
      drhon=-(rhoprof(N,1)-rhoprof(N,2))/(0.5+0.5*tanh((zk-zn1)/dzn))
      do k=1,N
        zk=asymth(z_r(istr,jstrR,k),zn1,dzn_a)
        rhoprof(k,2)=rhoprof(k,2)-drhon*(0.5+0.5*tanh((zk-zn1)/dzn))
      enddo
      do k=1,N
        zk = z_r(istr,jstrR,k)
        rhoprof(k,1)=rhoprof(k,1)
     &             -drhosf*(exp((zk-z0p)/abs(z0p)))/(exp(1.))
        rhoprof(k,1)=rhoprof(k,1)
     &             -drhosfs*0.5*(1+tanh((zk-z0)/abs(z0)))/tanh(1.)
        rhoprof(k,2)=rhoprof(k,2)
     &             -drhosf*(exp((zk-z0p)/abs(z0p)))/(exp(1.))
        rhoprof(k,2)=rhoprof(k,2)
     &             -drhosfn*0.5*(1+tanh((zk-z0)/abs(z0)))/tanh(1.)
      enddo
!
! --- Fit southern and northern profiles ---
!
      do k=1,N
        do j=JstrR,JendR
          do i=IstrR,IendR
#  ifdef MPI
            y=j-0.5+jj*Mm
            x=i-0.5+ii*Lm
#  else
            y=j-0.5
            x=i-0.5
#  endif
            y=y/real(MMm)-0.5
            x=x/real(LLm)
!
! Add vertical attenuation and XI perturbation
!
            y = y + cff_perturb*exp(z_r(i,j,k)/1000.)
     &                          *exp(-(x-0.5)**2/0.05)
     &                          *( 0.5*sin(2.*pi*x)
     &                            +0.5*sin(6.*pi*x) )
!
! Set Jet structure and place its center at y=pi/2
!
            y = y*pi*el/Ljet + pi/2 

            if (y.lt.0) then
              Fyz=1.0
            elseif (y.gt.pi) then
              Fyz=0.0
            else
              Fyz=1.0-(y-sin(y)*cos(y))/pi
            endif
            rho(i,j,k)=Fyz*rhoprof(k,1)+(1.-Fyz)*rhoprof(k,2)
          enddo
        enddo
      enddo

      do k=1,N
        do j=JstrR,JendR
          do i=IstrR,IendR
            t(i,j,k,tndx,1)=rho(i,j,k)
          enddo
        enddo
      enddo

#  if defined EW_PERIODIC || defined NS_PERIODIC || defined  MPI
      call exchange_r3d_tile (istr,iend,jstr,jend,
     &                        t(START_2D_ARRAY,1,tndx,1))
#  endif
      rho(:,:,:)=t(:,:,:,tndx,1)
!
! --- Compute zeta ---
!
! Get pressure using pressure gradient routine
!
      call prsgrd_tile (Istr,Iend,Jstr,Jend, ru,rv, P,
     &                  A2d(1,3,0),A2d(1,4,0), FC,dZx,rx,dRx)
!
! Compute zeta so as to cancel bottom pressure perturbation
!
      do j=Jstr,Jend
        do i=Istr,Iend
          zeta(i,j,tndx)=-P(i,j,1)/g
        enddo
      enddo
!
! Remove mean sea level (mssh)  !!! MPI DEBUG problems here !!!
!
      mssh=QuadZero
      area=QuadZero
      do j=Jstr,Jend
        do i=Istr,Iend
           cff1=1./(pm(i,j)*pn(i,j))
           mssh=mssh+cff1*zeta(i,j,tndx)
           area=area+cff1
        enddo
      enddo
      mssh=mssh/area
! aponte quad (
      mssh_dbl=dble(mssh)
      area_dbl=dble(area)
! aponte quad )
#  ifdef MPI
! xios (
!      call MPI_ALLGATHER(mssh,1,MPI_DOUBLE_PRECISION,
!     &                allmssh,1,MPI_DOUBLE_PRECISION,
!     &                          MPI_COMM_WORLD,ierr)
!      call MPI_ALLGATHER(area,1,MPI_DOUBLE_PRECISION,
!     &                allarea,1,MPI_DOUBLE_PRECISION,
!     &                          MPI_COMM_WORLD,ierr)
! aponte quad (
!      call MPI_ALLGATHER(mssh,1,MPI_DOUBLE_PRECISION,
!     &                allmssh,1,MPI_DOUBLE_PRECISION,
!     &                          ocean_grid_comm,ierr)
!      call MPI_ALLGATHER(area,1,MPI_DOUBLE_PRECISION,
!     &                allarea,1,MPI_DOUBLE_PRECISION,
!     &                          ocean_grid_comm,ierr)
      call MPI_ALLGATHER(mssh_dbl,1,MPI_DOUBLE_PRECISION,
     &                allmssh_dbl,1,MPI_DOUBLE_PRECISION,
     &                          ocean_grid_comm,ierr)
      call MPI_ALLGATHER(area_dbl,1,MPI_DOUBLE_PRECISION,
     &                allarea_dbl,1,MPI_DOUBLE_PRECISION,
     &                          ocean_grid_comm,ierr)
! aponte quad )
! xios )
! aponte quad (
!      mssh=QuadZero
!      area=QuadZero
!      do i=1,NNODES
!        mssh=mssh+allmssh(1,i)*allarea(1,i)
!        area=area+allarea(1,i)
!      enddo
!      mssh=mssh/area
      mssh_dbl=0.0D+0
      area_dbl=0.0D+0
      do i=1,NNODES
        mssh_dbl=mssh_dbl+allmssh_dbl(1,i)*allarea_dbl(1,i)
        area_dbl=area_dbl+allarea_dbl(1,i)
      enddo
      mssh_dbl=mssh_dbl/area_dbl
! aponte quad )
#  endif
      do j=Jstr,Jend
        do i=Istr,Iend
! aponte quad (
!          zeta(i,j,tndx)=zeta(i,j,tndx)-mssh
          zeta(i,j,tndx)=zeta(i,j,tndx)-mssh_dbl
! aponte quad )
        enddo
      enddo
!
! Apply zeta boundary conditions
!
      if (SOUTHERN_EDGE) then
        do i=Istr,Iend                     ! Southern edge gradient BC
          zeta(i,Jstr-1,tndx)=zeta(i,Jstr,tndx)
          enddo
      endif
      if (NORTHERN_EDGE) then
        do i=istr,iend                     ! Northern edge gradient BC
          zeta(i,Jend+1,tndx)=zeta(i,Jend,tndx)
        enddo
      endif
#  if defined EW_PERIODIC || defined NS_PERIODIC || defined  MPI
      call exchange_r2d_tile (istr,iend,jstr,jend,
     &                        zeta(START_2D_ARRAY,tndx))
#  endif
!
! --- Get gressure gradient (ru,rv) using updated zeta field ---
!
      call set_depth_tile (Istr,Iend,Jstr,Jend)
      call prsgrd_tile (Istr,Iend,Jstr,Jend, ru,rv, P,
     &                  A2d(1,3,0),A2d(1,4,0), FC,dZx,rx,dRx)

#  define WORK rho1
! Bug : rv déclaré PRIVATE_2D_SCRATCH_ARRAY plus petit que GLOBAL_2D_ARRAY
!      WORK(:,:,:) = rv(GLOBAL_2D_ARRAY,:)
       WORK(PRIVATE_2D_SCRATCH_ARRAY,:) = rv(PRIVATE_2D_SCRATCH_ARRAY,:)   
#  if defined EW_PERIODIC || defined NS_PERIODIC || defined  MPI
      call exchange_r3d_tile (Istr,Iend,Jstr,Jend,
     &                        WORK(START_2D_ARRAY,1))
#  endif
      if (SOUTHERN_EDGE) then
        WORK(:,0,:)=WORK(:,1,:)
      endif
      if (NORTHERN_EDGE) then
        WORK(:,Mm+1,:)=WORK(:,Mm,:)
      endif
!     rv(GLOBAL_2D_ARRAY,:) = WORK(:,:,:)
      rv(PRIVATE_2D_SCRATCH_ARRAY,:) = WORK(PRIVATE_2D_SCRATCH_ARRAY,:)            
!
! Bug : ru déclaré PRIVATE_2D_SCRATCH_ARRAY plus petit que GLOBAL_2D_ARRAY
!     WORK(:,:,:) = ru(GLOBAL_2D_ARRAY,:)
      WORK(PRIVATE_2D_SCRATCH_ARRAY,:) = ru(PRIVATE_2D_SCRATCH_ARRAY,:)      
#  if defined EW_PERIODIC || defined NS_PERIODIC || defined  MPI
      call exchange_r3d_tile (istr,iend,jstr,jend,
     &                        WORK(START_2D_ARRAY,1))
#  endif
      if (SOUTHERN_EDGE) then
        WORK(:,0,:)=0.
      endif
      if (NORTHERN_EDGE) then
        WORK(:,Mm+1,:)=0.
      endif
!     ru(GLOBAL_2D_ARRAY,:) = WORK(:,:,:)
      ru(PRIVATE_2D_SCRATCH_ARRAY,:) = WORK(PRIVATE_2D_SCRATCH_ARRAY,:)            
#  undef WORK
!
! --- Compute geostrophic u,v velocities ---
!
      do k=1,N
        do j=Jstr,Jend
          do i=Istr,Iend
            cff=0.5*(pm(i,j)+pm(i,j-1))
     &             *(pn(i,j)+pn(i,j-1))
     &            /(Hz(i,j,k)+Hz(i,j-1,k))/f(i,j)
            u(i,j,k,tndx) = cff*rv(i,j,k)
          enddo
        enddo
      enddo
      do k=1,N
        do j=Jstr,Jend
          do i=Istr,Iend
            cff=0.5*(pm(i,j)+pm(i-1,j))
     &                *(pn(i,j)+pn(i-1,j))
     &                /(Hz(i,j,k)+Hz(i-1,j,k))/f(i,j)
            v(i,j,k,tndx) =-cff*ru(i,j,k)
          enddo
        enddo
      enddo
!
! --- Compute barotropic u velocity ---
!
      ubar(:,:,tndx)=0.
      do k=1,N
        do j=Jstr,Jend
          do i=Istr,Iend
            ubar(i,j,tndx)=ubar(i,j,tndx)+u(i,j,k,tndx)
     &                             *0.5*(Hz(i-1,j,k)+Hz(i,j,k))
          enddo
        enddo
      enddo
      do j=Jstr,Jend
        do i=Istr,Iend
          ubar(i,j,tndx) = ubar(i,j,tndx) /
     &                  ( 0.5*(zeta(i,j,tndx)+zeta(i-1,j,tndx))
     &                  + 0.5*(h(i,j)+h(i-1,j)) )
        enddo
      enddo
!
! --- Compute barotropic u velocity ---
!
      vbar(:,:,tndx)=0.
      do k=1,N
        do j=Jstr,Jend
          do i=Istr,Iend
            vbar(i,j,tndx)=vbar(i,j,tndx)+v(i,j,k,tndx)
     &                           *0.5*(Hz(i,j-1,k)+Hz(i,j,k))
          enddo
        enddo
      enddo
      do j=Jstr,Jend
        do i=Istr,Iend
          vbar(i,j,tndx) = vbar(i,j,tndx) / 
     &                ( 0.5*(zeta(i,j,tndx)+zeta(i,j-1,tndx))
     &                + 0.5*(h(i,j)+h(i,j-1)) )
        enddo
      enddo

#  ifdef CLIMATOLOGY
      if (jet_perturb .eq. 0) then
!
! Fill climatological arrays. MPI exchanges are done in analytical
!
        do k=1,N
          do j=Jstr,Jend
            do i=Istr,Iend
              uclm(i,j,k) = u(i,j,k,tndx)
              vclm(i,j,k) = v(i,j,k,tndx)
              tclm(i,j,k,1) = t(i,j,k,tndx,1)
            enddo
          enddo
        enddo
        do j=Jstr,Jend
          do i=Istr,Iend
            ssh(i,j) = zeta(i,j,tndx)
          enddo
        enddo
        do j=Jstr,Jend
          do i=Istr,Iend
            ubclm(i,j) = ubar(i,j,tndx)
          enddo
        enddo
        do j=Jstr,Jend
          do i=Istr,Iend
            vbclm(i,j) = vbar(i,j,tndx)
          enddo
        enddo
      endif
#  endif

      return
      end
# endif /* JET */








! aponte jet (
# ifdef JETN
!*********************************************************************
      subroutine ana_jetn_tile (Istr,Iend,Jstr,Jend, jet_perturb)
!*********************************************************************
!
!--------------------------------------------------------------------
!  Set initial conditions for momentum and tracer variables using
!  analytical expressions for the JET test case
!--------------------------------------------------------------------
!
#  ifdef AGRIF
      use Agrif_Util
#  endif
      implicit none
#  include "param.h"
#  include "grid.h"
#  include "ocean2d.h"
#  include "ocean3d.h"
#  include "scalars.h"
#  include "mpi_cpl.h"
#  include "private_scratch.h"
#  ifdef CLIMATOLOGY
#   include "climat.h"
#  endif
      integer Istr,Iend,Jstr,Jend, i,j,k, itrc, 
     &        jet_perturb, tndx
      real cff,x,y,Fyz,
     &     rhomaxs,bgdrhos,zs1,dzs,drhos,drhosfs,
     &     rhomaxn,bgdrhon,zn1,dzn,drhon,drhosfn,
     &     drhosf,z0,z0p,dzs_a,dzn_a,zk,
     &     Ljet,cff_perturb,ierr,tanh
! aponte (
      integer trd, omp_get_thread_num
      real alpha1
      real alpha2
      real drhosfs1, drhosfs2, drhosfn1, drhosfn2
! aponte )
      real*QUAD mssh,area,cff1
! aponte quad (
      real mssh_dbl,area_dbl
! aponte quad )
      real rhoprof(0:N,2)
      real ru(PRIVATE_2D_SCRATCH_ARRAY,N),   
     &     rv(PRIVATE_2D_SCRATCH_ARRAY,N),   
     &      P(PRIVATE_2D_SCRATCH_ARRAY,N), 
     &     FC(PRIVATE_2D_SCRATCH_ARRAY),
     &    dZx(PRIVATE_2D_SCRATCH_ARRAY),
     &     rx(PRIVATE_2D_SCRATCH_ARRAY),
     &    dRx(PRIVATE_2D_SCRATCH_ARRAY)
#  ifdef MPI
      real*QUAD allmssh(1,NNODES),allarea(1,NNODES)
! aponte quad (
      real allmssh_dbl(1,NNODES),allarea_dbl(1,NNODES)
! aponte quad )
! aponte JET_EDDY (
# ifdef JET_EDDY
      real nbeddy, keddy, xx

# endif
! aponte JET_EDDY )
      include 'mpif.h'
#  endif
      real asymth,zz,zasym,dzasym
! aponte (
      asymth(zz,zasym,dzasym) =
     &         (zz-zasym)
     &         *(1+0.5*((zz-zasym)+abs(zz-zasym))/dzasym)
     &         +zasym
!      asymth(zz,zasym,dzasym) = 
!     &         (zz-zasym)
!     &         *sqrt(1+0.5*((zz-zasym)+abs(zz-zasym))**2/dzasym**2)
!     &         +zasym
! aponte )
!
#  include "compute_auxiliary_bounds.h"

# ifndef ANA_DENSITY_FILE

!
! ---------------------------------------------------------------------
! BAROCLINIC JET TEST CASE
! ---------------------------------------------------------------------
!
! --- Jet and density profile parameters ---
!
! aponte jet (
!      Ljet        = 1600.e3  ! Jet width in km
      Ljet        = jet_ywidth  ! Jet width in km
! aponte jet )
!
! southern profile
      rhomaxs=27.75     ! max density
      bgdrhos=9.8e-6    ! background stratification
      zs1=-1000
      dzs=700
!      drhos=1.4         ! gives strength of Phillips mode (0 to 1.4)
! northern profile 
      rhomaxn=27.7573   ! max density
      bgdrhon=9.8e-6    ! background stratification
      zn1=-400 
      dzn=300
! surface Charney mode 
!      drhosfs=1.5       ! northern anomaly  (0 for Phillips mode only)
!      drhosfn=0.0       ! southern anomaly  (0 for Phillips mode only)
!      z0=-300.          ! vertical penetration (-600 for strong Charney mode)
!      drhosf=0.00       ! weakens Charney mode (canceled for 1.05)
!      z0p=-110.         ! only used for drhosf non zero
! aponte jet (
!      z0p=z0 ! for consistency issues with roms_ird_fev2011
#  if JET_CONFIG == 1
      drhos=1.4         ! gives strength of Phillips mode (0 to 1.4)
      drhosfs=1.5       ! northern anomaly  (0 for Phillips mode only)
      drhosfn=0.0       ! southern anomaly  (0 for Phillips mode only)
      z0=-300.0          ! vertical penetration (-600 for strong Charney mode)
      z0p=z0 ! for consistency issues with roms_ird_fev2011
      alpha1=0.0075   !   weak charney mode => 1 otherwise 0 (or 0.0075 to avoid inversion)
      alpha2=1.0   ! strong charney mode => 1 otherwise 0
!
      drhosfs1=drhosfs*alpha1
      drhosfs2=drhosfs*alpha2
      drhosfn1=drhosfn*alpha1
      drhosfn2=drhosfn*alpha2
#  elif JET_CONFIG == 2
      drhos=1.4         ! gives strength of Phillips mode (0 to 1.4)
      drhosfs=1.5       ! northern anomaly  (0 for Phillips mode only)
      drhosfn=0.0       ! southern anomaly  (0 for Phillips mode only)
      z0=-300.0          ! vertical penetration (-600 for strong Charney mode)
      z0p=z0 ! for consistency issues with roms_ird_fev2011
      alpha1=0.0075   !   weak charney mode => 1 otherwise 0 (or 0.0075 to avoid inversion)
      alpha2=0.0   ! strong charney mode => 1 otherwise 0
!
      drhosfs1=drhosfs*alpha1
      drhosfs2=drhosfs*alpha2
      drhosfn1=drhosfn*alpha1
      drhosfn2=drhosfn*alpha2
#  elif JET_CONFIG == 3
      ! setup 1: Phillips > Charney
      drhos=1.4         ! gives strength of Phillips mode (0 to 1.4)
      drhosfs1=0.0
      drhosfn1=0.0
      drhosfs2=1.5
      drhosfn2=0.0
      z0=-300.0          ! vertical penetration (-600 for strong Charney mode)
      z0p=-110.0
#  elif JET_CONFIG == 4
      ! setup 2: Phillips
      drhos=1.4         ! gives strength of Phillips mode (0 to 1.4)
      drhosfs1=0.0
      drhosfn1=0.0
      drhosfs2=0.0
      drhosfn2=0.0
      z0=-300.0          ! vertical penetration (-600 for strong Charney mode)
      z0p=-110.0
#  elif JET_CONFIG == 5
      ! setup 4: Phillips + surface anomaly
      drhos=1.4         ! gives strength of Phillips mode (0 to 1.4)
      drhosfs1=0.71
      drhosfn1=0.71
      drhosfs2=1.5
      drhosfn2=0.0
      z0=-300.0          ! vertical penetration (-600 for strong Charney mode)
      z0p=-110.0
#  elif JET_CONFIG == 6
      ! setup 5
      drhos=0.0         ! gives strength of Phillips mode (0 to 1.4)
      drhosfs1=0.0
      drhosfn1=0.0
      drhosfs2=2.5
      drhosfn2=1.0
      z0=-600.0          ! vertical penetration (-600 for strong Charney mode)
      z0p=-110.0
#  endif
! A parameter that weakens the jet by averaging north and south profiles with
! a weight in order to obtain the final profile in the north
!      jet_weight = 0.5
! aponte jet )
!
! zonal perturbation
!
      if (jet_perturb .eq. 0) then
        cff_perturb = 0.0    ! zonal perturbation coefficient
        tndx=2
      else
        cff_perturb = 0.02
        tndx=1
        zeta(:,:,tndx)=0.
        call set_depth_tile (Istr,Iend,Jstr,Jend)
      endif
      
!     do j=JstrR,JendR
!       do i=IstrR,IendR      
!
! --- Build northern and southern density profiles ---
!
! First compute background density profiles associated with
! a gentle stratification that does not depend 
! on jet side. It is there to ensure static stability. 
      do k=1,N
        zk = z_r(istr,jstrR,k)
        rhoprof(k,1)=rhomaxs-bgdrhos*(zk+h(istr,jstrR))
        rhoprof(k,2)=rhomaxn-bgdrhon*(zk+h(istr,jstrR))
      enddo
!
! Second, get main north/south contrast with a distorded 
! tanh shape of variable amplitude.Distorsion of the tanh 
! is done with the asym function that increases 
! stratification in the upper ocean
      dzs_a=1.3*dzs
      do k=1,N
        zk=asymth(z_r(istr,jstrR,k),zs1,dzs_a)
        rhoprof(k,1)=rhoprof(k,1)-drhos*(0.5+0.5*tanh((zk-zs1)/dzs))
      enddo
      dzn_a=1.3*dzn
      drhon=-(rhoprof(N,1)-rhoprof(N,2))/(0.5+0.5*tanh((zk-zn1)/dzn))
      do k=1,N
        zk=asymth(z_r(istr,jstrR,k),zn1,dzn_a)
        rhoprof(k,2)=rhoprof(k,2)-drhon*(0.5+0.5*tanh((zk-zn1)/dzn))
      enddo
      do k=1,N
        zk = z_r(istr,jstrR,k)
        rhoprof(k,1)=rhoprof(k,1)
! aponte jet (
     &             -drhosfs1*(exp((zk-z0p)/abs(z0p)))/(exp(1.))
!     &             -drhosfs*alpha1*(exp((zk-z0p)/abs(z0p)))/(exp(1.))
!     &             -drhosf*(exp((zk-z0p)/abs(z0p)))/(exp(1.))
! aponte jet )
        rhoprof(k,1)=rhoprof(k,1)
! aponte jet (
     &             -drhosfs2*0.5*(1+tanh((zk-z0)/abs(z0)))/tanh(1.)
!     &             -drhosfs*alpha2*0.5*(1+tanh((zk-z0)/abs(z0)))/tanh(1.)
!     &             -drhosfs*0.5*(1+tanh((zk-z0)/abs(z0)))/tanh(1.)
! aponte jet )
        rhoprof(k,2)=rhoprof(k,2)
! aponte jet (
     &             -drhosfn1*(exp((zk-z0p)/abs(z0p)))/(exp(1.))
!     &             -drhosfn*alpha1*(exp((zk-z0p)/abs(z0p)))/(exp(1.))
!     &             -drhosf*(exp((zk-z0p)/abs(z0p)))/(exp(1.))
! aponte jet )
        rhoprof(k,2)=rhoprof(k,2)
! aponte jet (
     &             -drhosfn2*0.5*(1+tanh((zk-z0)/abs(z0)))/tanh(1.)
!     &             -drhosfn*alpha2*0.5*(1+tanh((zk-z0)/abs(z0)))/tanh(1.)
!     &             -drhosfn*0.5*(1+tanh((zk-z0)/abs(z0)))/tanh(1.)
! aponte jet )
! aponte jet (
! averages the profiles in order to have smaller difference
        rhoprof(k,2) = jet_weight * rhoprof(k,2) 
     &                 + (1.0 - jet_weight) * rhoprof(k,1)
! aponte jet )
      enddo
!
! --- Fit southern and northern profiles ---
!
      do k=1,N
        do j=JstrR,JendR
          do i=IstrR,IendR
#  ifdef MPI
            y=j-0.5+jj*Mm
            x=i-0.5+ii*Lm
#  else
            y=j-0.5
            x=i-0.5
#  endif
            y=y/real(MMm)-0.5
            x=x/real(LLm)

# ifdef JET_EDDY
!
! Set Eddy structure
!
             nbeddy = 10
             keddy = nbeddy * pi / xl
             x = x*xl
             y = y*el
             
!             if (nrrec==0) then
             if ( jet_perturb ) then
                call random_number(xx)
                !xx=0.
                Fyz = cos(keddy*x + xx/10.)*
     &                cos(keddy*y + xx/10.)
             else
                Fyz = cos(keddy*x)*
     &                cos(keddy*y)            
             endif
         
             Fyz = 1. - Fyz*(1+tanh((y+Ljet/2.)/50000.))*
     &                 (1-tanh((y-Ljet/2.)/50000.))/4.
#   ifdef FLAT_OCEAN
             Fyz=1.0
#   endif
! aponte jet )
#else
! aponte jet (
! I keep this part of the code for now even though it is different from fev2011 code
#   ifdef JET_PERTURB
!
! Add vertical attenuation and XI perturbation
!
        if ( flag_jet_perturb ) then
            y = y + cff_perturb*exp(z_r(i,j,k)/1000.)
     &                          *exp(-(x-0.5)**2/0.05)
     &                          *( 0.5*sin(2.*pi*x)
     &                            +0.5*sin(6.*pi*x) )
        endif
#   endif
! aponte jet )
!
! Set Jet structure and place its center at y=pi/2
!
             y = y*pi*el/Ljet + pi/2 

             if(y.lt.0)then
               Fyz=0.5
             elseif(y.gt.pi)then
               Fyz=-0.5
             else
               Fyz=0.5-(y-sin(y)*cos(y))/pi
             endif
! aponte jet (
             Fyz=Fyz+0.5
             if ( .not. flag_jet ) then
              Fyz=1.0
             endif
#   ifdef FLAT_OCEAN
             Fyz=1.0
#   endif
! aponte jet )

# endif /* JET_EDDY */
             rho(i,j,k)=Fyz*rhoprof(k,1)+(1.-Fyz)*rhoprof(k,2)
           enddo
         enddo
       enddo

       do k=1,N
         do j=JstrR,JendR
           do i=IstrR,IendR
             t(i,j,k,tndx,1)=rho(i,j,k)
           enddo
         enddo
       enddo

#  if defined EW_PERIODIC || defined NS_PERIODIC || defined  MPI
      call exchange_r3d_tile (istr,iend,jstr,jend,
     &                        t(START_2D_ARRAY,1,tndx,1))
#  endif
      rho(:,:,:)=t(:,:,:,tndx,1)

#else

      if (jet_perturb .eq. 0) then     
! climatology
       t(:,:,:,nrhs,itemp) = tclm(:,:,:,itemp)
#ifdef SALINITY
       t(:,:,:,nrhs,isalt) = tclm(:,:,:,isalt)
#endif
       tndx=2
      else
! initial conditions (perturbed clim) 
       t(:,:,:,nrhs,itemp) = t(:,:,:,3,itemp)
#ifdef SALINITY
       t(:,:,:,nrhs,isalt) = t(:,:,:,3,isalt)
#endif
! init depth profiles
       tndx=1
       zeta(:,:,tndx)=0.
       call set_depth_tile (Istr,Iend,Jstr,Jend)
      endif
! tndx=1,2 depending on clim or ini
      t(:,:,:,tndx,itemp) = t(:,:,:,nrhs,itemp)
#ifdef SALINITY
      t(:,:,:,tndx,isalt) = t(:,:,:,nrhs,isalt)
#endif

      trd=omp_get_thread_num()
      call rho_eos_tile (Istr,Iend,Jstr,Jend, A2d(1,1,trd),
     &                                        A2d(1,2,trd))

#endif

! aponte density )    

!     if (jet_perturb .eq. 0) t(:,:,:,nrhs,1) = t(:,:,:,tndx,1)
!     call rho_eos_tile (Istr,Iend,Jstr,Jend)
!
! --- Compute zeta ---
!
! Get pressure using pressure gradient routine
!
      call prsgrd_tile (Istr,Iend,Jstr,Jend, ru,rv, P,
     &                  A2d(1,3,0),A2d(1,4,0), FC,dZx,rx,dRx)
!
! Compute zeta so as to cancel bottom pressure perturbation
!
      do j=Jstr,Jend
        do i=Istr,Iend
          zeta(i,j,tndx)=-P(i,j,1)/g
!         zeta(i,j,tndx)= 0.         
        enddo
      enddo
!
! Remove mean sea level (mssh)  !!! MPI DEBUG problems here !!!
!
      mssh=QuadZero
      area=QuadZero
      do j=Jstr,Jend
        do i=Istr,Iend
           cff1=1./(pm(i,j)*pn(i,j))
           mssh=mssh+cff1*zeta(i,j,tndx)
           area=area+cff1
        enddo
      enddo
      mssh=mssh/area
! aponte quad (
      mssh_dbl=dble(mssh)
      area_dbl=dble(area)
! aponte quad )
#  ifdef MPI
! xios (
!      call MPI_ALLGATHER(mssh,1,MPI_DOUBLE_PRECISION,
!     &                allmssh,1,MPI_DOUBLE_PRECISION,
!     &                          MPI_COMM_WORLD,ierr)
!      call MPI_ALLGATHER(area,1,MPI_DOUBLE_PRECISION,
!     &                allarea,1,MPI_DOUBLE_PRECISION,
!     &                          MPI_COMM_WORLD,ierr)
! aponte quad (
!      call MPI_ALLGATHER(mssh,1,MPI_DOUBLE_PRECISION,
!     &                allmssh,1,MPI_DOUBLE_PRECISION,
!     &                          ocean_grid_comm,ierr)
!      call MPI_ALLGATHER(area,1,MPI_DOUBLE_PRECISION,
!     &                allarea,1,MPI_DOUBLE_PRECISION,
!     &                          ocean_grid_comm,ierr)
      call MPI_ALLGATHER(mssh_dbl,1,MPI_DOUBLE_PRECISION,
     &                allmssh_dbl,1,MPI_DOUBLE_PRECISION,
     &                          ocean_grid_comm,ierr)
      call MPI_ALLGATHER(area_dbl,1,MPI_DOUBLE_PRECISION,
     &                allarea_dbl,1,MPI_DOUBLE_PRECISION,
     &                          ocean_grid_comm,ierr)
! aponte quad )
! xios )
! aponte quad (
!      mssh=QuadZero
!      area=QuadZero
!      do i=1,NNODES
!        mssh=mssh+allmssh(1,i)*allarea(1,i)
!        area=area+allarea(1,i)
!      enddo
!      mssh=mssh/area
      mssh_dbl=0.0D+0
      area_dbl=0.0D+0
      do i=1,NNODES
        mssh_dbl=mssh_dbl+allmssh_dbl(1,i)*allarea_dbl(1,i)
        area_dbl=area_dbl+allarea_dbl(1,i)
      enddo
      mssh_dbl=mssh_dbl/area_dbl
! aponte quad )
#  endif
      do j=Jstr,Jend
        do i=Istr,Iend
! aponte quad (
!          zeta(i,j,tndx)=zeta(i,j,tndx)-mssh
          zeta(i,j,tndx)=zeta(i,j,tndx)-mssh_dbl
! aponte quad )
        enddo
      enddo
!
! Apply zeta boundary conditions
!
      if (SOUTHERN_EDGE) then
        do i=Istr,Iend                     ! Southern edge gradient BC
          zeta(i,Jstr-1,tndx)=zeta(i,Jstr,tndx)
          enddo
      endif
      if (NORTHERN_EDGE) then
        do i=istr,iend                     ! Northern edge gradient BC
          zeta(i,Jend+1,tndx)=zeta(i,Jend,tndx)
        enddo
      endif
#  if defined EW_PERIODIC || defined NS_PERIODIC || defined  MPI
      call exchange_r2d_tile (istr,iend,jstr,jend,
     &                        zeta(START_2D_ARRAY,tndx))
#  endif
!
! --- Get gressure gradient (ru,rv) using updated zeta field ---
!
      call set_depth_tile (Istr,Iend,Jstr,Jend)
      call prsgrd_tile (Istr,Iend,Jstr,Jend, ru,rv, P,
     &                  A2d(1,3,0),A2d(1,4,0), FC,dZx,rx,dRx)

#  define ANA_INI_GEO_APONTE

! new implementation, should decrease amplitude of transients

#  if defined ANA_INI_GEO_APONTE
! exchange pressure gradients
#  define LWORK rho1
!
      LWORK(PRIVATE_2D_SCRATCH_ARRAY,:) = ru(:,:,:)
      call exchange_u3d_tile (Istr,Iend,Jstr,Jend,
     &                                 LWORK(START_2D_ARRAY,1))
      ru(:,:,:) = LWORK(PRIVATE_2D_SCRATCH_ARRAY,:)
!
      LWORK(PRIVATE_2D_SCRATCH_ARRAY,:) = rv(:,:,:)
      call exchange_v3d_tile (Istr,Iend,Jstr,Jend,
     &                                 rho1(START_2D_ARRAY,1))
      rv(:,:,:) = LWORK(PRIVATE_2D_SCRATCH_ARRAY,:)
!
#  undef LWORK

!
! --- Compute geostrophic u,v velocities ---
!
      do k=1,N
        do j=Jstr-1,Jend+1
          do i=Istr-1,Iend+1
            cff=0.5*(pm(i,j)+pm(i-1,j))
     &                *(pn(i,j)+pn(i-1,j))
     &                /(Hz(i,j,k)+Hz(i-1,j,k))/f(i,j)
            ru(i,j,k) =cff*ru(i,j,k)
!
            cff=0.5*(pm(i,j)+pm(i,j-1))
     &             *(pn(i,j)+pn(i,j-1))
     &            /(Hz(i,j,k)+Hz(i,j-1,k))/f(i,j)
            rv(i,j,k) = cff*rv(i,j,k)
          enddo
        enddo
      enddo
      do k=1,N
        do j=Jstr,Jend
          do i=Istr,Iend
            u(i,j,k,tndx) = ( rv(i,j,k)+rv(i,j+1,k)
     &                       +rv(i-1,j,k)+rv(i-1,j+1,k) ) * 0.25
            v(i,j,k,tndx) =-( ru(i,j,k)+ru(i+1,j,k)
     &                       +ru(i,j-1,k)+ru(i+1,j-1,k) ) * 0.25
          enddo
        enddo
      enddo

#  else

! old implementation

#  define WORK rho1
! Bug : rv déclaré PRIVATE_2D_SCRATCH_ARRAY plus petit que GLOBAL_2D_ARRAY
!      WORK(:,:,:) = rv(GLOBAL_2D_ARRAY,:)
       WORK(PRIVATE_2D_SCRATCH_ARRAY,:) = rv(PRIVATE_2D_SCRATCH_ARRAY,:)
#  if defined EW_PERIODIC || defined NS_PERIODIC || defined  MPI
      call exchange_r3d_tile (Istr,Iend,Jstr,Jend,
     &                        WORK(START_2D_ARRAY,1))
#  endif
      if (SOUTHERN_EDGE) then
        WORK(:,0,:)=WORK(:,1,:)
      endif
      if (NORTHERN_EDGE) then
        WORK(:,Mm+1,:)=WORK(:,Mm,:)
      endif
!     rv(GLOBAL_2D_ARRAY,:) = WORK(:,:,:)
      rv(PRIVATE_2D_SCRATCH_ARRAY,:) = WORK(PRIVATE_2D_SCRATCH_ARRAY,:)
!
! Bug : ru déclaré PRIVATE_2D_SCRATCH_ARRAY plus petit que GLOBAL_2D_ARRAY
!     WORK(:,:,:) = ru(GLOBAL_2D_ARRAY,:)
      WORK(PRIVATE_2D_SCRATCH_ARRAY,:) = ru(PRIVATE_2D_SCRATCH_ARRAY,:)
#  if defined EW_PERIODIC || defined NS_PERIODIC || defined  MPI
      call exchange_r3d_tile (istr,iend,jstr,jend,
     &                        WORK(START_2D_ARRAY,1))
#  endif
      if (SOUTHERN_EDGE) then
        WORK(:,0,:)=0.
      endif
      if (NORTHERN_EDGE) then
        WORK(:,Mm+1,:)=0.
      endif
!     ru(GLOBAL_2D_ARRAY,:) = WORK(:,:,:)
      ru(PRIVATE_2D_SCRATCH_ARRAY,:) = WORK(PRIVATE_2D_SCRATCH_ARRAY,:)
#  undef WORK
!
! --- Compute geostrophic u,v velocities ---
!
      do k=1,N
        do j=Jstr,Jend
          do i=Istr,Iend
            cff=0.5*(pm(i,j)+pm(i,j-1))
     &             *(pn(i,j)+pn(i,j-1))
     &            /(Hz(i,j,k)+Hz(i,j-1,k))/f(i,j)
            u(i,j,k,tndx) = cff*rv(i,j,k)
          enddo
        enddo
      enddo
      do k=1,N
        do j=Jstr,Jend
          do i=Istr,Iend
            cff=0.5*(pm(i,j)+pm(i-1,j))
     &                *(pn(i,j)+pn(i-1,j))
     &                /(Hz(i,j,k)+Hz(i-1,j,k))/f(i,j)
            v(i,j,k,tndx) =-cff*ru(i,j,k)
          enddo
        enddo
      enddo

#  endif
!
! --- Compute barotropic u velocity ---
!
      ubar(:,:,tndx)=0.
      do k=1,N
        do j=Jstr,Jend
          do i=Istr,Iend
            ubar(i,j,tndx)=ubar(i,j,tndx)+u(i,j,k,tndx)
     &                             *0.5*(Hz(i-1,j,k)+Hz(i,j,k))
          enddo
        enddo
      enddo
      do j=Jstr,Jend
        do i=Istr,Iend
          ubar(i,j,tndx) = ubar(i,j,tndx) /
     &                  ( 0.5*(zeta(i,j,tndx)+zeta(i-1,j,tndx))
     &                  + 0.5*(h(i,j)+h(i-1,j)) )
        enddo
      enddo
!
! --- Compute barotropic u velocity ---
!
      vbar(:,:,tndx)=0.
      do k=1,N
        do j=Jstr,Jend
          do i=Istr,Iend
            vbar(i,j,tndx)=vbar(i,j,tndx)+v(i,j,k,tndx)
     &                           *0.5*(Hz(i,j-1,k)+Hz(i,j,k))
          enddo
        enddo
      enddo
      do j=Jstr,Jend
        do i=Istr,Iend
          vbar(i,j,tndx) = vbar(i,j,tndx) / 
     &                ( 0.5*(zeta(i,j,tndx)+zeta(i,j-1,tndx))
     &                + 0.5*(h(i,j)+h(i,j-1)) )
        enddo
      enddo

#  ifdef CLIMATOLOGY
      if (jet_perturb .eq. 0) then
!
! Fill climatological arrays. MPI exchanges are done in analytical
!
! aponte jet (  there are issues with the indices here
!        do k=1,N
!          do j=Jstr,Jend
!            do i=Istr,Iend
!              uclm(i,j,k) = u(i,j,k,tndx)
!              vclm(i,j,k) = v(i,j,k,tndx)
!              tclm(i,j,k,1) = t(i,j,k,tndx,1)
!            enddo
!          enddo
!        enddo
        do k=1,N
          do j=JstrR,JendR
            do i=IstrR,IendR
              uclm(i,j,k) = u(i,j,k,tndx)
              vclm(i,j,k) = v(i,j,k,tndx)
              tclm(i,j,k,itemp) = t(i,j,k,tndx,itemp)
#ifdef SALINITY
              tclm(i,j,k,isalt) = t(i,j,k,tndx,isalt)
#endif
            enddo
          enddo
        enddo
!        do j=Jstr,Jend
!          do i=Istr,Iend
!            ssh(i,j) = zeta(i,j,tndx)
!          enddo
!        enddo
!        do j=Jstr,Jend
!          do i=Istr,Iend
!            ubclm(i,j) = ubar(i,j,tndx)
!          enddo
!        enddo
!        do j=Jstr,Jend
!          do i=Istr,Iend
!            vbclm(i,j) = vbar(i,j,tndx)
!          enddo
!        enddo
        do j=JstrR,JendR
          do i=IstrR,IendR
            ssh(i,j) = zeta(i,j,tndx)
            ubclm(i,j) = ubar(i,j,tndx)
            vbclm(i,j) = vbar(i,j,tndx)
          enddo
        enddo
      endif
! aponte )
#  endif  /* CLIMATOLOGY */



      return
      end
# endif /* JETN */
! aponte jet )





! aponte jetr (
# ifdef JETR
!*********************************************************************
      subroutine ana_jetr_tile (Istr,Iend,Jstr,Jend, jet_perturb)
!*********************************************************************
!
!--------------------------------------------------------------------
!  Set initial conditions for momentum and tracer variables using
!  analytical expressions for the JET test case
!--------------------------------------------------------------------
!
#  ifdef AGRIF
      use Agrif_Util
#  endif
      implicit none
#  include "param.h"
#  include "grid.h"
#  include "ocean2d.h"
#  include "ocean3d.h"
#  include "scalars.h"
#  include "mpi_cpl.h"
#  include "private_scratch.h"
#  ifdef CLIMATOLOGY
#   include "climat.h"
#  endif
      integer Istr,Iend,Jstr,Jend, i,j,k, itrc, 
     &        jet_perturb, tndx
      real cff,x,y,Fyz,
     &     zk,
     &     rhomax,bgdrho,
     &     Ljet,cff_perturb,ierr,tanh
      real*QUAD mssh,area,cff1
! aponte quad (
      real mssh_dbl,area_dbl
! aponte quad )
      real rhoprof(0:N,2)
      real ru(PRIVATE_2D_SCRATCH_ARRAY,N),   
     &     rv(PRIVATE_2D_SCRATCH_ARRAY,N),   
     &      P(PRIVATE_2D_SCRATCH_ARRAY,N), 
     &     FC(PRIVATE_2D_SCRATCH_ARRAY),
     &    dZx(PRIVATE_2D_SCRATCH_ARRAY),
     &     rx(PRIVATE_2D_SCRATCH_ARRAY),
     &    dRx(PRIVATE_2D_SCRATCH_ARRAY)
#  ifdef MPI
      real*QUAD allmssh(1,NNODES),allarea(1,NNODES)
! aponte quad (
      real allmssh_dbl(1,NNODES),allarea_dbl(1,NNODES)
! aponte quad )
      include 'mpif.h'
#  endif
! aponte )
!
#  include "compute_auxiliary_bounds.h"
!
! ---------------------------------------------------------------------
! BAROCLINIC JET TEST CASE
! ---------------------------------------------------------------------
!
! --- Jet and density profile parameters ---
!
! Background stratification
      bgdrho=1.0e-4    ! background stratification
      rhomax=27.75     ! max density
!
! depth, width and strength of the jet
      Ljet = width_jet  ! Jet width in km
!      H_jet
!      drho_jet

!
! zonal perturbation
!
      if (jet_perturb .eq. 0) then
        cff_perturb = 0.0    ! zonal perturbation coefficient
        tndx=2
      else
!        cff_perturb = 0.02
        cff_perturb = offy_jet/el
        tndx=1
        zeta(:,:,tndx)=0.
        call set_depth_tile (Istr,Iend,Jstr,Jend)
      endif
!
! --- Build northern and southern density profiles ---
!
! Linear stratification
      do k=1,N
        zk = z_r(istr,jstrR,k)
        rhoprof(k,1)=rhomax-bgdrho*(zk+h(istr,jstrR))
        rhoprof(k,2)=rhomax-bgdrho*(zk+h(istr,jstrR))
      enddo
! Add to the north the perturbation that will produce the jet
      do k=1,N
        zk = z_r(istr,jstrR,k)
        rhoprof(k,2)=rhoprof(k,2) - 
     &          drho_jet * exp(zk/H_jet)
      enddo

!
! --- Fit southern and northern profiles ---
!
      do k=1,N
        do j=JstrR,JendR
          do i=IstrR,IendR
#  ifdef MPI
            y=j-0.5+jj*Mm
            x=i-0.5+ii*Lm
#  else
            y=j-0.5
            x=i-0.5
#  endif
            y=y/real(MMm)-0.5
            x=x/real(LLm)

! aponte jet (
! I keep this part of the code for now even though it is different from fev2011 code
#   ifdef JET_PERTURB
!
! Add vertical attenuation and XI perturbation
!
        if ( flag_jet_perturb ) then
            y = y + cff_perturb
     &            *exp(-(x-0.5)**2/0.05)
!     &            *exp(z_r(i,j,k)/100.)
!     &                          *( 0.5*sin(2.*pi*x)
!     &                            +0.5*sin(6.*pi*x) )
        endif
#   endif
! aponte jet )
!
! Set Jet structure and place its center at y=pi/2
!
             y = y*pi*el/Ljet + pi/2 

             if(y.lt.0)then
               Fyz=0.5
             elseif(y.gt.pi)then
               Fyz=-0.5
             else
               Fyz=0.5-(y-sin(y)*cos(y))/pi
             endif
! aponte jet (
             Fyz=Fyz+0.5
             if ( .not. flag_jet ) then
              Fyz=1.0
             endif
#   ifdef FLAT_OCEAN
             Fyz=1.0
#   endif
! aponte jet )
             rho(i,j,k)=Fyz*rhoprof(k,1)+(1.-Fyz)*rhoprof(k,2)
           enddo
         enddo
       enddo

       do k=1,N
         do j=JstrR,JendR
           do i=IstrR,IendR
             t(i,j,k,tndx,1)=rho(i,j,k)
           enddo
         enddo
       enddo

#  if defined EW_PERIODIC || defined NS_PERIODIC || defined  MPI
      call exchange_r3d_tile (istr,iend,jstr,jend,
     &                        t(START_2D_ARRAY,1,tndx,1))
#  endif
      rho(:,:,:)=t(:,:,:,tndx,1)
!
! --- Compute zeta ---
!
! Get pressure using pressure gradient routine
!
      call prsgrd_tile (Istr,Iend,Jstr,Jend, ru,rv, P,
     &                  A2d(1,3,0),A2d(1,4,0), FC,dZx,rx,dRx)
!
! Compute zeta so as to cancel bottom pressure perturbation
!
      do j=Jstr,Jend
        do i=Istr,Iend
          zeta(i,j,tndx)=-P(i,j,1)/g
        enddo
      enddo
!
! Remove mean sea level (mssh)  !!! MPI DEBUG problems here !!!
!
      mssh=QuadZero
      area=QuadZero
      do j=Jstr,Jend
        do i=Istr,Iend
           cff1=1./(pm(i,j)*pn(i,j))
           mssh=mssh+cff1*zeta(i,j,tndx)
           area=area+cff1
        enddo
      enddo
      mssh=mssh/area
! aponte quad (
      mssh_dbl=dble(mssh)
      area_dbl=dble(area)
! aponte quad )
#  ifdef MPI
! xios (
!      call MPI_ALLGATHER(mssh,1,MPI_DOUBLE_PRECISION,
!     &                allmssh,1,MPI_DOUBLE_PRECISION,
!     &                          MPI_COMM_WORLD,ierr)
!      call MPI_ALLGATHER(area,1,MPI_DOUBLE_PRECISION,
!     &                allarea,1,MPI_DOUBLE_PRECISION,
!     &                          MPI_COMM_WORLD,ierr)
! aponte quad (
!      call MPI_ALLGATHER(mssh,1,MPI_DOUBLE_PRECISION,
!     &                allmssh,1,MPI_DOUBLE_PRECISION,
!     &                          ocean_grid_comm,ierr)
!      call MPI_ALLGATHER(area,1,MPI_DOUBLE_PRECISION,
!     &                allarea,1,MPI_DOUBLE_PRECISION,
!     &                          ocean_grid_comm,ierr)
      call MPI_ALLGATHER(mssh_dbl,1,MPI_DOUBLE_PRECISION,
     &                allmssh_dbl,1,MPI_DOUBLE_PRECISION,
     &                          ocean_grid_comm,ierr)
      call MPI_ALLGATHER(area_dbl,1,MPI_DOUBLE_PRECISION,
     &                allarea_dbl,1,MPI_DOUBLE_PRECISION,
     &                          ocean_grid_comm,ierr)
! aponte quad )
! xios )
! aponte quad (
!      mssh=QuadZero
!      area=QuadZero
!      do i=1,NNODES
!        mssh=mssh+allmssh(1,i)*allarea(1,i)
!        area=area+allarea(1,i)
!      enddo
!      mssh=mssh/area
      mssh_dbl=0.0D+0
      area_dbl=0.0D+0
      do i=1,NNODES
        mssh_dbl=mssh_dbl+allmssh_dbl(1,i)*allarea_dbl(1,i)
        area_dbl=area_dbl+allarea_dbl(1,i)
      enddo
      mssh_dbl=mssh_dbl/area_dbl
! aponte quad )
#  endif
      do j=Jstr,Jend
        do i=Istr,Iend
! aponte quad (
!          zeta(i,j,tndx)=zeta(i,j,tndx)-mssh
          zeta(i,j,tndx)=zeta(i,j,tndx)-mssh_dbl
! aponte quad )
        enddo
      enddo
!
! Apply zeta boundary conditions
!
      if (SOUTHERN_EDGE) then
        do i=Istr,Iend                     ! Southern edge gradient BC
          zeta(i,Jstr-1,tndx)=zeta(i,Jstr,tndx)
          enddo
      endif
      if (NORTHERN_EDGE) then
        do i=istr,iend                     ! Northern edge gradient BC
          zeta(i,Jend+1,tndx)=zeta(i,Jend,tndx)
        enddo
      endif
#  if defined EW_PERIODIC || defined NS_PERIODIC || defined  MPI
      call exchange_r2d_tile (istr,iend,jstr,jend,
     &                        zeta(START_2D_ARRAY,tndx))
#  endif
!
! --- Get gressure gradient (ru,rv) using updated zeta field ---
!
      call set_depth_tile (Istr,Iend,Jstr,Jend)
      call prsgrd_tile (Istr,Iend,Jstr,Jend, ru,rv, P,
     &                  A2d(1,3,0),A2d(1,4,0), FC,dZx,rx,dRx)

#  define ANA_INI_GEO_APONTE

! new implementation, should decrease amplitude of transients

#  if defined ANA_INI_GEO_APONTE
! exchange pressure gradients
#  define LWORK rho1
!
      LWORK(PRIVATE_2D_SCRATCH_ARRAY,:) = ru(:,:,:)
      call exchange_u3d_tile (Istr,Iend,Jstr,Jend,
     &                                 LWORK(START_2D_ARRAY,1))
      ru(:,:,:) = LWORK(PRIVATE_2D_SCRATCH_ARRAY,:)
!
      LWORK(PRIVATE_2D_SCRATCH_ARRAY,:) = rv(:,:,:)
      call exchange_v3d_tile (Istr,Iend,Jstr,Jend,
     &                                 rho1(START_2D_ARRAY,1))
      rv(:,:,:) = LWORK(PRIVATE_2D_SCRATCH_ARRAY,:)
!
#  undef LWORK

!
! --- Compute geostrophic u,v velocities ---
!
      do k=1,N
        do j=Jstr-1,Jend+1
          do i=Istr-1,Iend+1
            cff=0.5*(pm(i,j)+pm(i-1,j))
     &                *(pn(i,j)+pn(i-1,j))
     &                /(Hz(i,j,k)+Hz(i-1,j,k))/f(i,j)
            ru(i,j,k) =cff*ru(i,j,k)
!
            cff=0.5*(pm(i,j)+pm(i,j-1))
     &             *(pn(i,j)+pn(i,j-1))
     &            /(Hz(i,j,k)+Hz(i,j-1,k))/f(i,j)
            rv(i,j,k) = cff*rv(i,j,k)
          enddo
        enddo
      enddo
      do k=1,N
        do j=Jstr,Jend
          do i=Istr,Iend
            u(i,j,k,tndx) = ( rv(i,j,k)+rv(i,j+1,k)
     &                       +rv(i-1,j,k)+rv(i-1,j+1,k) ) * 0.25
            v(i,j,k,tndx) =-( ru(i,j,k)+ru(i+1,j,k)
     &                       +ru(i,j-1,k)+ru(i+1,j-1,k) ) * 0.25
          enddo
        enddo
      enddo

#  else

! old implementation

#  define WORK rho1
! Bug : rv déclaré PRIVATE_2D_SCRATCH_ARRAY plus petit que GLOBAL_2D_ARRAY
!      WORK(:,:,:) = rv(GLOBAL_2D_ARRAY,:)
       WORK(PRIVATE_2D_SCRATCH_ARRAY,:) = rv(PRIVATE_2D_SCRATCH_ARRAY,:)         
#  if defined EW_PERIODIC || defined NS_PERIODIC || defined  MPI
      call exchange_r3d_tile (Istr,Iend,Jstr,Jend,
     &                        WORK(START_2D_ARRAY,1))
#  endif
      if (SOUTHERN_EDGE) then
        WORK(:,0,:)=WORK(:,1,:)
      endif
      if (NORTHERN_EDGE) then
        WORK(:,Mm+1,:)=WORK(:,Mm,:)
      endif
!     rv(GLOBAL_2D_ARRAY,:) = WORK(:,:,:)
      rv(PRIVATE_2D_SCRATCH_ARRAY,:) = WORK(PRIVATE_2D_SCRATCH_ARRAY,:)            

!
! Bug : ru déclaré PRIVATE_2D_SCRATCH_ARRAY plus petit que GLOBAL_2D_ARRAY
!     WORK(:,:,:) = ru(GLOBAL_2D_ARRAY,:)
      WORK(PRIVATE_2D_SCRATCH_ARRAY,:) = ru(PRIVATE_2D_SCRATCH_ARRAY,:)   
#  if defined EW_PERIODIC || defined NS_PERIODIC || defined  MPI
      call exchange_r3d_tile (istr,iend,jstr,jend,
     &                        WORK(START_2D_ARRAY,1))
#  endif
      if (SOUTHERN_EDGE) then
        WORK(:,0,:)=0.
      endif
      if (NORTHERN_EDGE) then
        WORK(:,Mm+1,:)=0.
      endif
!     ru(GLOBAL_2D_ARRAY,:) = WORK(:,:,:)
      ru(PRIVATE_2D_SCRATCH_ARRAY,:) = WORK(PRIVATE_2D_SCRATCH_ARRAY,:)          
#  undef WORK
!
! --- Compute geostrophic u,v velocities ---
!
      do k=1,N
        do j=Jstr,Jend
          do i=Istr,Iend
            cff=0.5*(pm(i,j)+pm(i,j-1))
     &             *(pn(i,j)+pn(i,j-1))
     &            /(Hz(i,j,k)+Hz(i,j-1,k))/f(i,j)
            u(i,j,k,tndx) = cff*rv(i,j,k)
          enddo
        enddo
      enddo
      do k=1,N
        do j=Jstr,Jend
          do i=Istr,Iend
            cff=0.5*(pm(i,j)+pm(i-1,j))
     &                *(pn(i,j)+pn(i-1,j))
     &                /(Hz(i,j,k)+Hz(i-1,j,k))/f(i,j)
            v(i,j,k,tndx) =-cff*ru(i,j,k)
          enddo
        enddo
      enddo

#  endif
!
! --- Compute barotropic u velocity ---
!
      ubar(:,:,tndx)=0.
      do k=1,N
        do j=Jstr,Jend
          do i=Istr,Iend
            ubar(i,j,tndx)=ubar(i,j,tndx)+u(i,j,k,tndx)
     &                             *0.5*(Hz(i-1,j,k)+Hz(i,j,k))
          enddo
        enddo
      enddo
      do j=Jstr,Jend
        do i=Istr,Iend
          ubar(i,j,tndx) = ubar(i,j,tndx) /
     &                  ( 0.5*(zeta(i,j,tndx)+zeta(i-1,j,tndx))
     &                  + 0.5*(h(i,j)+h(i-1,j)) )
        enddo
      enddo
!
! --- Compute barotropic u velocity ---
!
      vbar(:,:,tndx)=0.
      do k=1,N
        do j=Jstr,Jend
          do i=Istr,Iend
            vbar(i,j,tndx)=vbar(i,j,tndx)+v(i,j,k,tndx)
     &                           *0.5*(Hz(i,j-1,k)+Hz(i,j,k))
          enddo
        enddo
      enddo
      do j=Jstr,Jend
        do i=Istr,Iend
          vbar(i,j,tndx) = vbar(i,j,tndx) / 
     &                ( 0.5*(zeta(i,j,tndx)+zeta(i,j-1,tndx))
     &                + 0.5*(h(i,j)+h(i,j-1)) )
        enddo
      enddo

#  ifdef CLIMATOLOGY
      if (jet_perturb .eq. 0) then
!
! Fill climatological arrays. MPI exchanges are done in analytical
!
! aponte jet (  there are issues with the indices here
!        do k=1,N
!          do j=Jstr,Jend
!            do i=Istr,Iend
!              uclm(i,j,k) = u(i,j,k,tndx)
!              vclm(i,j,k) = v(i,j,k,tndx)
!              tclm(i,j,k,1) = t(i,j,k,tndx,1)
!            enddo
!          enddo
!        enddo
        do k=1,N
          do j=JstrR,JendR
            do i=IstrR,IendR
              uclm(i,j,k) = u(i,j,k,tndx)
              vclm(i,j,k) = v(i,j,k,tndx)
              tclm(i,j,k,1) = t(i,j,k,tndx,1)
            enddo
          enddo
        enddo
!        do j=Jstr,Jend
!          do i=Istr,Iend
!            ssh(i,j) = zeta(i,j,tndx)
!          enddo
!        enddo
!        do j=Jstr,Jend
!          do i=Istr,Iend
!            ubclm(i,j) = ubar(i,j,tndx)
!          enddo
!        enddo
!        do j=Jstr,Jend
!          do i=Istr,Iend
!            vbclm(i,j) = vbar(i,j,tndx)
!          enddo
!        enddo
        do j=JstrR,JendR
          do i=IstrR,IendR
            ssh(i,j) = zeta(i,j,tndx)
            ubclm(i,j) = ubar(i,j,tndx)
            vbclm(i,j) = vbar(i,j,tndx)
          enddo
        enddo
      endif
! aponte )
#  endif



      return
      end
# endif /* JETR */
! aponte jetr )





#else
      subroutine ana_initial_empty
      return
      end
#endif /* ANA_INITIAL */
