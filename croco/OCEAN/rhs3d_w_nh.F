#include "cppdefs.h"

#define W_HADV
#define W_VADV

#ifndef NBQ_FREESLIP
# define NBQ_WNOSLIP
#else
# undef  NBQ_WNOSLIP
#endif

#if defined SOLVE3D && defined NBQ

      SUBROUTINE rhs3d_w_nh(tile)

!======================================================================
!                       ***  Subroutine RHS3D  ***
! NBQ mode :  compute right-hand-side for the vertical velocity wz 
!======================================================================
! History :   2016-11  (F. Lemarié)  Original code
!----------------------------------------------------------------------

      implicit none
      integer tile, trd, omp_get_thread_num
# include "param.h"
# include "private_scratch.h"
# include "compute_tile_bounds.h"
      trd=omp_get_thread_num()
      CALL rhs3d_w_tile (Istr,Iend,Jstr,Jend,  A3d(1,3,trd),    !<-- A3d(1,1,trd), A3d(1,2,trd) contain ru,rv 
     &             A2d(1,1,trd), A2d(1,2,trd), A2d(1,3,trd),    !<-- A3d(1,3,trd) contains rw 
     &             A2d(1,4,trd), A2d(1,5,trd), A2d(1,6,trd), 
     &             A2d(1,7,trd), A3d(1,4,trd), A3d(1,5,trd)     !<-- A3d(1,3,trd) and A3d(1,4,trd) should 
     &                                                     )    !<-- be available at this point ...
      return
      END
      
!======================================================================
      
      SUBROUTINE rhs3d_w_tile (istr, iend,   jstr,   jend, rw ,
     &                           CF,   FC,  We_r ,    lap, Hzw, 
     &                          WFx,  WFe, HUon_w, HVom_w       )      
      IMPLICIT NONE
      !!
      INTEGER  :: Istr,Iend,Jstr,Jend
      INTEGER  :: i,j,k,kp
      INTEGER  :: imin,imax,jmin,jmax
      real    rr,rrm,rrp,limiter2,Cr,cdt,cu
      integer kp2,kp1,km1
      !!
# include "param.h"
      REAL     :: rw     (PRIVATE_2D_SCRATCH_ARRAY,0:N  )
      REAL     :: CF     (PRIVATE_1D_SCRATCH_ARRAY,0:N  )
      REAL     :: FC     (PRIVATE_1D_SCRATCH_ARRAY,0:N+1)  
      REAL     :: We_r   (PRIVATE_1D_SCRATCH_ARRAY,0:N+1)        
      REAL     :: lap    (PRIVATE_2D_SCRATCH_ARRAY      )
      REAL     :: Hzw    (PRIVATE_1D_SCRATCH_ARRAY,0:N  )      
      REAL     :: WFx    (PRIVATE_2D_SCRATCH_ARRAY      ) 
      REAL     :: WFe    (PRIVATE_2D_SCRATCH_ARRAY      )
      REAL     :: HUon_w (PRIVATE_2D_SCRATCH_ARRAY,0:N  )
      REAL     :: HVom_w (PRIVATE_2D_SCRATCH_ARRAY,0:N  )
      !!    
      REAL     :: gamma, epsil, cff, cff1, cff2, cff3, Omeg_r
      REAL     :: wz_bot,wz_sfc
      PARAMETER (gamma=0.25  )
      PARAMETER (epsil=1.E-16)
      !!
# include "grid.h"
# include "ocean3d.h"
# include "coupling.h"
# include "forces.h"
# include "scalars.h"  
# include "nbq.h"
# include "compute_auxiliary_bounds.h"
      !!
# ifdef EW_PERIODIC 
#   define IU_RANGE Istr,Iend
#   define IV_RANGE Istr,Iend
# else
#   define IU_RANGE Istr,IendR
#   define IV_RANGE IstrR,IendR
# endif

# ifdef NS_PERIODIC 
#   define JU_RANGE Jstr,Jend
#   define JV_RANGE Jstr,Jend
# else
#   define JU_RANGE JstrR,JendR
#   define JV_RANGE Jstr,JendR
# endif

      DO j = jstr,jend
        DO i = istr,iend+1
          HUon_w(i,j,0)=0.5*HUon(i,j,1)                
          HUon_w(i,j,N)=0.5*HUon(i,j,N) 
        ENDDO
      ENDDO
      !
      DO j = jstr,jend+1
        DO i = istr,iend
          HVom_w(i,j,0)=0.5*HVom(i,j,1)             
          HVom_w(i,j,N)=0.5*HVom(i,j,N) 
        ENDDO
      ENDDO 
      !
      DO k = 1,N-1
        DO j = jstr,jend
          DO i = istr,iend+1
            HUon_w(i,j,k)=0.5*( HUon(i,j,k)+HUon(i,j,k+1) )
          ENDDO
        ENDDO
         
        DO j = jstr,jend+1
          DO i = istr,iend
            HVom_w(i,j,k)=0.5*( HVom(i,j,k)+HVom(i,j,k+1) )
          ENDDO
        ENDDO         
      ENDDO

# if defined UV_COR
!
!-----------------------------------------------------------------------
!  Add in Coriolis  FWcor = +e (U cos(angle) - V sin(angle) ) 
!  with e = 2 Omega cos Phi
!-----------------------------------------------------------------------
!
      DO k = 0, N
        DO j = jstr,jend   
          DO i = istr,iend        
            rw(i,j,k) = 0.5*eomn(i,j)*(
     &          cosa(i,j)*pn(i,j)*(HUon_w(i,j,k)+HUon_w(i+1,j,k)) 
     &        - sina(i,j)*pm(i,j)*(HVom_w(i,j,k)+HVom_w(i,j+1,k))
     &                                 )             
          ENDDO
        ENDDO
      ENDDO      
# else
      DO k = 0, N
         DO j = jstr,jend   
            DO i = istr,iend        
              rw(i,j,k) = 0.            
            ENDDO
         ENDDO
      ENDDO 
# endif

# if defined UV_ADV && defined W_HADV
!
!-----------------------------------------------------------------------      
!  Horizontal advection
!----------------------------------------------------------------------- 
!
      if(nrhs.eq.3) then
        cdt = dt                     !<-- Corrector
      elseif(FIRST_TIME_STEP) then
        cdt = 0.5*dt
      else
        cdt = (1.-2*gamma)*dt        !<-- Predictor
      endif

#  ifndef EW_PERIODIC
      IF (WESTERN_EDGE) THEN       
         imin=istr                  
      ELSE                        
         imin=istr-1             
      ENDIF                      
      IF (EASTERN_EDGE) THEN      
         imax=iend                 
      ELSE                  
         imax=iend+1           
      ENDIF
#  else
      imin=istr-1
      imax=iend+1
#  endif

#  ifndef NS_PERIODIC
      IF (SOUTHERN_EDGE) THEN
         jmin=jstr
      ELSE
         jmin=jstr-1
      ENDIF
      IF (NORTHERN_EDGE) then
         jmax=jend
      ELSE
         jmax=jend+1
      ENDIF
#  else
      jmin=jstr-1
      jmax=jend+1
#  endif
 
      DO k=0,N             

#  ifdef W_HADV_C2
!
! === C2 horizontal advection scheme ===
!
        DO j=Jstr,Jend
          DO i=Istr,Iend+1
            WFx(i,j)=0.5*( wz(i,j,k,nrhs)+wz(i-1,j,k,nrhs)
     &                                     )*Huon_w(i,j,k)
          ENDDO    
        ENDDO
        DO j=jstr,jend+1
          DO i=istr,iend
            WFe(i,j)=0.5*( Wz(i,j,k,nrhs)+Wz(i,j-1,k,nrhs)
     &                                     )*Hvom_w(i,j,k)
          ENDDO
        ENDDO

#  elif defined W_HADV_TVD
!
! === TVD horizontal advection scheme ===
!
        DO j=Jstr,Jend
          DO i=Istr,Iend+1
            WFx(i,j)=0.5*( wz(i,j,k,nrhs)+wz(i-1,j,k,nrhs)
     &                                     )*Huon_w(i,j,k)
          ENDDO    
        ENDDO

        DO j=Jstr,Jend
          DO i=Istr,Iend+1
            if (k==0) then
              Hzw (i,k) = 0.25 * (HZR(i-1,j,k+1) + HZR(i,j,k+1))
            elseif (k==N) then
              Hzw (i,k) = 0.25 * (HZR(i-1,j,k) + HZR(i,j,k))
            else
              Hzw (i,k) = 0.25 * ( HZR(i  ,j,k) + HZR(i  ,j,k+1)
     &                           + HZR(i-1,j,k) + HZR(i-1,j,k+1) )
            endif
            cff = Huon_w(i,j,k)
            cu = pn_u(i,j)*cdt*cff
     &           /Hzw(i,k)*pm(i,j)

            rrp= (wz(i+1,j,k,nrhs) -wz(i,j,k,nrhs))
#   ifdef MASKING
     &                                *umask(i+1,j)
#   endif
            rr = (wz(i,j,k,nrhs) -wz(i-1 ,j,k,nrhs))
#   ifdef MASKING
     &                                *umask(i,j)
#   endif
            rrm= (wz(i-1,j,k,nrhs) -wz(i-2 ,j,k,nrhs))
#   ifdef MASKING
     &                                *umask(i-1,j)
#   endif

            cff1=(cff*(wz(i-1,j,k,nrhs) + wz(i,j,k,nrhs) )
     &                                   -ABS(cff)*rr)*0.5

            Cr=limiter2(cu,WFx(i,j),cff1,rrm,rr,rrp)
            WFx(i,j) = (1-cr)*cff1 + Cr* WFx(i,j)
          ENDDO
        ENDDO
!---
        DO j=jstr,jend+1
          DO i=istr,iend
            WFe(i,j)=0.5*( Wz(i,j,k,nrhs)+Wz(i,j-1,k,nrhs)
     &                                     )*Hvom_w(i,j,k)
          ENDDO
        ENDDO 

        DO j=JstrV-1,Jend+1
          DO i=Istr,Iend
            if (k==0) then
              Hzw(i,k) = 0.25 * (HZR(i,j,k+1) + HZR(i,j-1,k+1))
            elseif (k==N) then
              Hzw(i,k) = 0.25 * (HZR(i,j,k  ) + HZR(i,j-1,k  ))
            else
              Hzw(i,k) = 0.25 * (HZR(i,j  ,k) + HZR(i,j  ,k+1)
     &                         + HZR(i,j-1,k) + HZR(i,j-1,k+1) )
            endif
            cff = Hvom_w(i,j,k)
            cu = pm_v(i,j)*cdt*cff
     &           /Hzw(i,k)*pn(i,j)

            rrp= (wz(i,j+1,k,nrhs) -wz(i,j,k,nrhs))
#   ifdef MASKING
     &                                *vmask(i,j+1)
#   endif
            rr = (wz(i,j,k,nrhs) -wz(i,j-1,k,nrhs))
#   ifdef MASKING
     &                                *vmask(i,j)
#   endif
            rrm= (wz(i,j-1,k,nrhs) -wz(i,j-2,k,nrhs))
#   ifdef MASKING
     &                                *vmask(i,j-1)
#   endif
            cff1=(cff*(wz(i,j,k,nrhs) + wz(i,j-1,k,nrhs) )
     &                                  -ABS(cff)*rr )*0.5
            Cr=limiter2(cu,WFe(i,j),cff1,rrm,rr,rrp)
            WFe(i,j) = (1-Cr)*cff1  + Cr* WFe(i,j)
          ENDDO
        ENDDO

#  else /* W_HADV_C4 */
!
! === TVD horizontal advection scheme ===
!
        DO j=jstr,jend
          DO i=imin,imax+1
            WFx(i,j) = ( wz(i,j,k,nrhs) - wz(i-1,j,k,nrhs) )
#   ifdef MASKING
     &                                           *umask(i,j)
#   endif
          ENDDO
        ENDDO  

#   ifndef EW_PERIODIC
        IF (WESTERN_EDGE) THEN
          DO j=jstr,jend
            WFx(istr-1,j) = WFx(istr  ,j)
          ENDDO
        ENDIF
        IF (EASTERN_EDGE) then
          DO j=jstr,jend
            WFx(iend+2,j) = WFx(iend+1,j)
          ENDDO
        ENDIF
#   endif

        DO j=Jstr,Jend
          DO i=Istr-1,Iend+1          
            lap(i,j) = WFx(i+1,j) - WFx(i,j)
          ENDDO
        ENDDO  

        DO j=Jstr,Jend
          DO i=Istr,Iend+1
            IF ( Huon_w(i,j,k).gt. 0. ) THEN
              cff = lap( i-1, j )
            ELSE
              cff = lap( i  , j )
            ENDIF            
            WFx(i,j)=0.5*( wz(i,j,k,nrhs)+wz(i-1,j,k,nrhs)
     &                           -0.25*cff )*Huon_w(i,j,k)
          ENDDO    
        ENDDO
         
        DO j=jmin,jmax+1
          DO i=istr,iend
            WFe(i,j) = ( Wz(i,j,k,nrhs)-Wz(i,j-1,k,nrhs) ) 
#   ifdef MASKING
     &                                         *vmask(i,j)
#   endif
          ENDDO 
        ENDDO         
#   ifndef NS_PERIODIC
        IF (SOUTHERN_EDGE) then
          DO i=istr,iend
            WFe(i,jstr-1) = WFe(i,jstr)
          ENDDO
        ENDIF
        IF (NORTHERN_EDGE) then
          DO i=istr,iend
            WFe(i,jend+2) = WFe(i,jend+1)
          ENDDO
        ENDIF
#   endif         
 
        DO j=jstr-1,jend+1
          DO i=istr,iend
            lap(i,j) = WFe(i,j+1) - WFe(i,j  )
          ENDDO
        ENDDO  

        DO j=jstr,jend+1
          DO i=istr,iend
            IF ( HVom_w(i,j,k).gt. 0. ) THEN
               cff = lap( i , j-1 )
            ELSE
               cff = lap( i , j   )
            ENDIF 
            WFe(i,j)=0.5*( Wz(i,j,k,nrhs)+Wz(i,j-1,k,nrhs)
     &                           -0.25*cff )*Hvom_w(i,j,k)
          ENDDO
        ENDDO 
#  endif /* W_HADV_C2 */
!
!  Finalize horizontal advection: compute flux divergence
!
        DO j=Jstr,Jend
          DO i=Istr,Iend
            rw(i,j,k)= rw(i,j,k) - WFx(i+1,j) + WFx(i,j)
     &                           - WFe(i,j+1) + WFe(i,j)    
          ENDDO
        ENDDO

      ENDDO          !<-- outer loop k

# endif /* W_HADV */


# if defined UV_ADV && defined W_VADV
!
!-----------------------------------------------------------------------      
!  Vertical advection
!----------------------------------------------------------------------- 
!
#  define udH_x WFx
#  define vdH_y WFe  
#  define DUon(i,j) Huon_w(i,j,0)  
#  define DVom(i,j) Hvom_w(i,j,0)
#  define udZ_x(i,j) Huon_w(i,j,1)  
#  define vdZ_y(i,j) Hvom_w(i,j,1)
#  define FC1 lap
!
!----------------------------------------------
!  First, set bottom boundary condition for wz:
!  free-slip or no-slip
!----------------------------------------------
!
      IF(nrhs.eq.3) then    !<-- corrector
        DO j=Jstr,Jend 
          DO i=Istr,Iend+1 
#  ifndef NBQ_WNOSLIP
            udH_x(i,j) = - u(i,j,1,nrhs) * pm_u(i,j) 
     &                                   * ( z_w(i  ,j,0)-z_w(i-1,j,0) )
#  else
            udH_x(i,j) = 0. 
#  endif
#  ifdef NBQ_MASS
            DUon(i,j) = 2.*DU_avg2(i,j) 
     &             / ( rhobar_nbq_avg1(i,j) + rhobar_nbq_avg1(i-1,j) ) 
#  else
            DUon(i,j) = DU_avg2(i,j) 
#  endif
            udZ_x(i,j) = u(i,j,N,nrhs) * pm_u(i,j) 
     &                                 * ( z_w(i  ,j,N)-z_w(i-1,j,N) )
          ENDDO
        ENDDO

        DO j=Jstr,Jend+1 
          DO i=Istr,Iend 
#  ifndef NBQ_WNOSLIP
            vdH_y(i,j) = - v(i,j,1,nrhs) * pn_v(i,j) 
     &                                   * ( z_w(i,j,0)-z_w(i,j-1,0) )
#  else
            vdH_y(i,j) = 0. 
#  endif
#  ifdef NBQ_MASS
            DVom(i,j)  = 2.*DV_avg2(i,j)
     &           / ( rhobar_nbq_avg1(i,j) + rhobar_nbq_avg1(i,j-1) )
#  else
            DVom(i,j)  = DV_avg2(i,j)
#  endif
            vdZ_y(i,j) = v(i,j,N,nrhs) * pn_v(i,j) 
     &                                 * ( z_w(i,j,N)-z_w(i,j-1,N) )
          ENDDO
        ENDDO
      ELSE                  !<-- predictor
        DO j=Jstr,Jend 
          DO i=Istr,Iend+1 
#  ifndef NBQ_WNOSLIP
            udH_x(i,j) = - u(i,j,1,nstp) * pm_u(i,j) 
     &                                   * ( z_w(i  ,j,0)-z_w(i-1,j,0) )
#  else
            udH_x(i,j) = 0.
#  endif
#  ifdef NBQ_MASS
            DUon(i,j) = 2.*DU_avg1(i,j,nstp) 
     &             / ( rhobar_nbq_avg1(i,j) + rhobar_nbq_avg1(i-1,j) )  
#  else
            DUon(i,j) = DU_avg1(i,j,nstp) 
#  endif
            udZ_x(i,j) = u(i,j,N,nstp) * pm_u(i,j) 
     &                                 * ( z_w(i  ,j,N)-z_w(i-1,j,N) )
          ENDDO
        ENDDO

        DO j=Jstr,Jend+1 
          DO i=Istr,Iend 
#  ifndef NBQ_WNOSLIP
            vdH_y(i,j) = - v(i,j,1,nstp) * pn_v(i,j) 
     &                                   * ( z_w(i,j,0)-z_w(i,j-1,0) )
#  else
            vdH_y(i,j) = 0.
#  endif
#  ifdef NBQ_MASS
            DVom(i,j) = 2.*DV_avg1(i,j,nstp)
     &           / ( rhobar_nbq_avg1(i,j) + rhobar_nbq_avg1(i,j-1) )
#  else
            DVom(i,j) = DV_avg1(i,j,nstp)
#  endif
            vdZ_y(i,j) = v(i,j,N,nstp) * pn_v(i,j) 
     &                                 * ( z_w(i,j,N)-z_w(i,j-1,N) )
          ENDDO
        ENDDO
      ENDIF    

#  ifdef NBQ_FREESLIP
      DO j=Jstr,Jend
        DO i=Istr,Iend
          wz(i,j,0,nnew) = - 0.5d0 * ( udH_x(i,j) + udH_x(i+1,j)
     &                             +   vdH_y(i,j) + vdH_y(i,j+1) )
     &                             *0.5*Hz(i,j,1)
        ENDDO
      ENDDO
#  else
      DO j=Jstr,Jend
        DO i=Istr,Iend
          wz(i,j,0,nnew) = 0.
        ENDDO
      ENDDO
#  endif
!
!--------------------------------------------
!  Second, compute vertical advection
!--------------------------------------------
!
      DO j=Jstr,Jend 
 
        DO i=Istr,Iend          
          Hzw (i,0  ) = 0.5 * HZR(i,j,1) 
          Hzw (i,N  ) = 0.5 * HZR(i,j,N)  
          We_r(i,0  ) =        We(i,j,0)
          We_r(i,N  ) = 0.5 * (We(i,j,N)+We(i,j,N-1))             
          We_r(i,N+1) =        We(i,j,N)        
        ENDDO
        DO k=1,N-1
          DO i=Istr,Iend
            Hzw (i,k) = 0.5d0 * ( HZR(i,j,k) + HZR(i,j,k+1) )
            We_r(i,k) = 0.5d0 * ( We (i,j,k) + We (i,j,k-1) ) 
          ENDDO        
        ENDDO

#  ifdef W_VADV_SPLINES          
!
! === SPLINES vertical advection scheme ===
!
        DO i=Istr,Iend      
          wz_bot  = - 0.5d0 * ( udH_x(i,j) + udH_x(i+1,j)
     &                      +   vdH_y(i,j) + vdH_y(i,j+1) )  
#   ifdef MOVING_BATHY
          wz_bot  =  wz_bot - dh(i,j)     
#   endif  
          cff     =   0.5 / ( Hzw (i,1) + Hzw (i,0) )  !<--- 1 / b(1)
          CF(i,1) =   cff * Hzw(i,0)                 !<--- q(1) = c(1) / b(1)
          FC(i,1) =   cff * (  3.*Hzw(i,0)*wz(i,j,1,nrhs)     
     &                       + 3.*Hzw(i,1)*wz(i,j,0,nrhs)
     &                       -    Hzw(i,1)*wz_bot      !<--- f(1) / b(1) 
     &                      )     
        ENDDO          

        DO k=2,N-1
          DO i=Istr,Iend          
            cff     = 1./( 2.*Hzw(i,k-1)+Hzw(i,k)*(2.-CF(i,k-1)) )   !<--- p = 1 / ( b(k)+a(k)*q(k-1) )
            CF(i,k) = cff* Hzw(i,k-1)                                !<--- c(k) * p
            FC(i,k) = cff*( 3.*( Hzw(i,k-1)*wz(i,j,k  ,nrhs)         !<--- f(k)=( f(k)-a(k)*f(k-1) )*p
     &                         + Hzw(i,k  )*wz(i,j,k-1,nrhs) )
     &                         - Hzw(i,k  )*FC(i,k-1)  )
          ENDDO
        ENDDO                   
          
        DO i=Istr,Iend           
          wz_sfc  = pn(i,j)*pm(i,j)*( DUon(i,j)-DUon(i+1,j  )
     &                            +   DVom(i,j)-DVom(i  ,j+1) )    !
     &            + 0.5d0 * ( udZ_x(i,j) + udZ_x(i+1,j)
     &                    +   vdZ_y(i,j) + vdZ_y(i,j+1) )

          cff     = 1./( 2.*Hzw(i,N-1)+Hzw(i,N)*(2.-CF(i,N-1)) )   !<--- p = 1 / ( b(k)+a(k)*q(k-1) )
          FC(i,N) = cff*( 3.*( Hzw(i,N-1)*wz(i,j,N  ,nrhs)         !<--- f(k)=( f(k)-a(k)*f(k-1) )*p
     &                       + Hzw(i,N  )*wz(i,j,N-1,nrhs) )
     &                       - Hzw(i,N-1)*wz_sfc
     &                       - Hzw(i,N  )*FC(i,N-1)  )    
        ENDDO          

        DO k=N-1,1,-1
          DO i=Istr,Iend 
            FC(i,k  )=FC(i,k)-CF(i,k)*FC(i,k+1)
            FC(i,k+1)=FC(i,k+1)*We_r(i,k+1)  
          ENDDO
        ENDDO

        DO i=Istr,Iend 
          FC(i,1  )=FC(i,1  )*We_r(i,1)
          FC(i,0  )=0. 
          FC(i,N+1)=0.
        ENDDO

#  else
!
! === C2 vertical advection scheme ===
!
        DO k=2,N-1
          DO i=Istr,Iend
            FC(i,k) = 0.5*(wz(i,j,k,nrhs) + wz(i,j,k-1,nrhs))*We_r(i,k)
#   ifdef W_VADV_TVD
            FC1(i,k)=(We_r(i,k) * (wz(i,j,k,nrhs) + wz(i,j,k-1,nrhs) )
     &           -abs(We_r(i,k))* (wz(i,j,k,nrhs) - wz(i,j,k-1,nrhs) ) 
     &               ) *0.5
#   endif
          ENDDO
        ENDDO

        DO i=Istr,Iend        
!         wz_bot  = - 0.5d0 * ( udH_x(i,j) + udH_x(i+1,j)
!     &                     +   vdH_y(i,j) + vdH_y(i,j+1) )  
!#   ifdef MOVING_BATHY
!         wz_bot  =  wz_bot - dh(i,j)     
!#   endif 
        ! FC(i,1) = 0.5*( wz(i,j,1,nrhs) + wz_bot  )*We_r(i,1)
          FC(i,1) = 0.5*( wz(i,j,1,nrhs) + wz(i,j,0,nrhs)  )*We_r(i,1)
#   ifdef W_VADV_TVD
          FC1(i,1)=(We_r(i,1) * (wz(i,j,1,nrhs) + wz(i,j,0,nrhs) )
     &         -abs(We_r(i,1))* (wz(i,j,1,nrhs) - wz(i,j,0,nrhs) ) )*0.5
#   endif
!         wz_sfc  = pn(i,j)*pm(i,j)*( DUon(i,j)-DUon(i+1,j  )
!     &                           +   DVom(i,j)-DVom(i  ,j+1) )    !
!     &           + 0.5d0 * ( udZ_x(i,j) + udZ_x(i+1,j)
!     &                   +   vdZ_y(i,j) + vdZ_y(i,j+1) )
!        FC(i,N) = 0.5*( wz_sfc + wz(i,j,N-1,nrhs))*We_r(i,N)
          FC(i,N) = 0.5*( wz(i,j,N,nrhs) + wz(i,j,N-1,nrhs))*We_r(i,N)
#   ifdef W_VADV_TVD
          FC1(i,N)=(We_r(i,N) * (wz(i,j,N,nrhs) + wz(i,j,N-1,nrhs) )
     &         -abs(We_r(i,N))* (wz(i,j,N,nrhs) - wz(i,j,N-1,nrhs) ) )*0.5
#   endif
        ENDDO          

        DO i=Istr,Iend
          FC(i,0  )=0.
          FC(i,N+1)=0.
        ENDDO          

#   ifdef W_VADV_TVD
        DO k=1,N
          kp2=min(N+1,k+2)
          kp1=min(N+1,k+1)
          km1=max(1,k-1)
          DO i=IstrU,Iend
            cff= We_r(i,k)
            cu = cff*pm(i,j)*pn(i,j)
     &            *cdt/Hzw(i,k)
            rrp= (wz(i,j,min(kp2,N),nrhs) -wz(i,j,min(kp1,N),nrhs))
            rr = (wz(i,j,min(kp1,N),nrhs) -wz(i ,j,k,nrhs))
            rrm= (wz(i,j,k  ,nrhs) -wz(i ,j,km1,nrhs))
            Cr = limiter2(cu,FC(i,k),cff1,rrm,rr,rrp)
            FC(i,k) = (1-Cr)*FC1(i,k) + Cr*FC(i,k)
          ENDDO
        ENDDO
#   endif

#  endif /* W_VADV_SPLINES */

#  undef udH_x
#  undef vdH_y   
#  undef DUon 
#  undef DVom 
#  undef udZ_x  
#  undef vdZ_y
#  undef FC1
!
!-----------------------------------------------------
! Finalize vertical advection: compite flux divergence
!-----------------------------------------------------
! 
        DO k=1,N
          DO i=Istr,Iend
            rw(i,j,k) = rw(i,j,k)  - FC(i,k+1) + FC(i,k)
#  ifdef MASKING              
            rw(i,j,k) = rw(i,j,k) * rmask(i,j)
#  endif             
          ENDDO
        ENDDO

!       DO i=Istr,Iend
!         rw(i,j,0) = rw(i,j,0)  + FC(i,0)
!#  ifdef MASKING             
!         rw(i,j,0) = rw(i,j,0) * rmask(i,j)
!#  endif             
!       ENDDO

      ENDDO  !<-- j loop
                     
# endif /* UV_ADV && W_VADV */

      return
      end
!
!=======================================================================
!
# if defined W_VADV_TVD || defined W_VADV_TVD

      function limiter2(nu,hf,lf,rrm,rr,rrp)

      implicit none
      real    :: cff,limiter2,Rj,nu,hf,lf,rrm,rr,rrp
!
!  Upwind        Limiter(Cr)=0.
!  Lax-Wendroff  Limiter(Cr)=1.
!  Min-Mod       Limiter(Cr)=max(0.,min(1.,Cr))
!  Suberbee      Limiter(Cr)=max(0.,max(min(1.,2*Cr*(1-nu)/nu),min(2.,Cr)))
!
      if (abs(rr).gt.1.e-20) THEN
        if (nu .gt. 1.e-20) THEN
          Rj=rrm/rr
        else
          Rj=rrp/rr
        endif
      else
        if (nu.gt.1.e-20) then
          Rj=Rrm*1.E20
        else
          Rj=Rrp*1.E20
        endif
      endif
       
!     limiter=max(0., min(1.,(hf-lf)/rr,(1-abs(nu))/abs(nu)*Rj) )
      cff=(1-abs(nu))/abs(nu)
!  Add stability condition/tvd Rj => Rj*cff

#  ifdef SUPERBEE
       limiter2=max(0.,max(min(1.,2.*Rj),
     &                     min(2.,Rj)))
#  elif defined MINMOD
       limiter2=max(0.D0,min(1.,Rj))
#  else
       limiter2=(Rj+abs(Rj))/(1+abs(Rj))
#  endif       
      end function

# endif /* W_VADV_TVD || W_VADV_TVD */
!
!=======================================================================
!
#else
      subroutine rhs3d_w_empty
      end
#endif /* SOLVE3D && NBQ */       
