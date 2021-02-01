! $Id: ana_initial.F 1620 2015-01-08 10:47:13Z marchesiello $
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

      subroutine ana_initial (tile)

      implicit none
#include "param.h"
       integer tile, trd
C$    integer omp_get_thread_num
#include "compute_tile_bounds.h"
#ifdef GO_MEDDY
      call ana_initial_tile (Istr,Iend,Jstr,Jend,tile)
#else     
      call ana_initial_tile (Istr,Iend,Jstr,Jend) 
#endif          
      return
      end
!
#ifdef GO_MEDDY
      subroutine ana_initial_tile (Istr,Iend,Jstr,Jend,tile)
#else
      subroutine ana_initial_tile (Istr,Iend,Jstr,Jend)
#endif
!
!----------------------------------------------------------------------
!  Set initial conditions for momentum and tracer variables using
!  analytical expressions.
!----------------------------------------------------------------------
!
      use params
#ifdef AGRIF
      use Agrif_Util
#endif
      implicit none
#include "param.h"
#include "grid.h"
#include "ocean2d.h"
#include "ocean3d.h"
#include "scalars.h"
#ifdef SEDIMENT
# include "sediment.h"
#endif      
#ifdef GO_MEDDY
      include "mpif.h"
#endif
      integer Istr,Iend,Jstr,Jend, i,j,k, itrc
      
#if defined NBQ && ! defined MPI
      integer :: mynode
#endif


#ifdef GO_MEDDY
      
      integer :: icase, tile , ierr

      real(8) :: valmax, valmin, nr_maxofall,nr_minofall, rand_spec
      real(8) :: thick_y, thick_z
      real(8) :: Burger, Rossby, stretch
      real(8) ::  f0    
      real(8) :: grd_z,buoyancy,pres_dif
      real(8) :: pres_up(GLOBAL_2D_ARRAY)
      real(8) :: cnt
      real(8),dimension(GLOBAL_2D_ARRAY,N) :: ph  
      real(8) :: x,y,z
      real(8) ::  x1,y1,z1
      
! tableau pour rhopr et gradient wind initialisation
      real(8),dimension(GLOBAL_2D_ARRAY,0:N) :: denom
      real(8) :: z2, d2, gss
      logical :: linear_balance
      real(8) :: u_geo(GLOBAL_2D_ARRAY,N)
      real(8) :: v_geo(GLOBAL_2D_ARRAY,N)
  
! Nvar
      integer ::  power
      real(8) :: dhh,shift,deltaN,zth1,zth2,delta1,delta2,wid1,wid2
      real(8) :: Fp, dNfreq_bk, Bu0
      real(8) :: zz,Nfreq2_bk_m1,Nfreq2_bk_p1,dNfreq2_bk, rhobk(N), 
     &           rhoin(N)
      real(8),dimension(1:N) :: rho_bk_sav
      real(8),dimension(GLOBAL_2D_ARRAY,N) :: gaussian

! random_seed
      logical   :: random_init, random_pert_urot
      integer   :: mde,nmdef,ni,nj,iofset,jofset
      integer*8 :: planf, planb
      real(8)   :: facsolpert, factor_amp_urot, w2, xlambda, xkrad1, mk, nk
      
      real(kind=8),dimension(1:64,1:64) :: kk
      real(8),dimension(N) :: xlmb , ekin
      real(kind=8)          :: e,etot,energy,xk,aa,xx,phase,alp,amp
      real(kind=8),dimension(1:LLm,1:MMm) :: rvar    
      complex, allocatable  :: prd_spec(:,:,:)
      complex, allocatable  :: prd(:,:,:)
      complex(kind=8),dimension(LLm/2+1,1:MMm) :: cvar
      
#endif
!
#include "compute_auxiliary_bounds.h"
!
#if defined EW_PERIODIC && !defined MPI
# define IR_RANGE Istr,Iend
# define IU_RANGE Istr,Iend
#else
# define IR_RANGE IstrR,IendR
# define IU_RANGE Istr,IendR
#endif
#if defined NS_PERIODIC && !defined MPI
# define JR_RANGE Jstr,Jend
# define JV_RANGE Jstr,Jend
#else
# define JR_RANGE JstrR,JendR
# define JV_RANGE Jstr,JendR
#endif

#ifdef ANA_INITIAL
!
!----------------------------------------------------------------------
!  Initial conditions for free surface and 2D momentum components.
!----------------------------------------------------------------------
!
      if (nrrec.eq.0) then
!
        do j=JR_RANGE
          do i=IR_RANGE
            zeta(i,j,1)=0.
            zeta(i,j,2)=0.
            ubar(i,j,1)=0.
            ubar(i,j,2)=0.
            vbar(i,j,1)=0.
            vbar(i,j,2)=0.
          enddo
        enddo


# ifdef NBQ
        do j=JR_RANGE
          do i=IR_RANGE
            zeta(i,j,1)=h(i,j)+zeta(i,j,1)
            zeta(i,j,2)=zeta(i,j,1)
          enddo
        enddo
# endif

      endif  !  nrrec.eq.0.


!
!  Initial conditions for tracer type variables.
!--------------------------------------------------------------------
!  Set initial conditions for potential temperature [degC] and
!  salinity [PSU].
!
# ifdef GO_MEDDY

      if (nrrec.eq.0) then
!
!====================================
!  input parameters of various icase
!==================================== 
#include "nr_init_icases.h"
!

      if (mynode == 0) then
         write(STDOUT,*) '--------------------------------------------------'
         write(STDOUT,*) ' initialisation vortex lens '
         write(STDOUT,*) '--------------------------------------------------'
         write(STDOUT,'(A,F7.1)')' largest time step (s) horizontal',
     &                  (1./pm(1,1))/(abs(2*z1*Nfreq/pi)+0.5)/4
         write(STDOUT,'(A,F7.1)')' largest time step (s) vertical', 
     &                  Hz(1,1,1)/(0.1)/4
         write(STDOUT,*)' max coriolis, min', f(1,1),  f(1,Mm)
         write(STDOUT,*) '--------------------------------------------------'
      endif
      call params_output('icase', icase)
      call params_output('Nfreq2',Nfreq2)
      call params_output('thick_y',thick_y)
      call params_output('thick_z',thick_z)
      call params_output('Burger',Burger)
      call params_output('Rossby',Rossby)
      call params_output('stretch',stretch)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% debut du background en z

! background rho(z)
!------------------
      if (icase == 50) then           ! Go meddy stratification   
        open(unit=12, file='exp.txt', form='formatted')
        read(12,*) rhobk
        read(12,*) rhoin
        close(12)
      else
        do k = 1, N
          do j = JR_RANGE
            do i = IR_RANGE
               z = z_r(i,j,k)
               zz = (z-zmid)/dhh
               rhobk(k) = (- Nfreq2*( (1-shift)**2*dhh*zz  
     &          + 2*(1-shift)*deltaN*dhh*log(cosh(zz))  
     &          + dhh*deltaN**2*(zz-tanh(zz)) )  
     &          - delta1*atan(exp(  wid1*(z+zth1)  ))  
     &          - delta2*atan(exp(  wid2*(z+zth2)  ))   )
            enddo
          enddo
        enddo
      endif

! background t,s
!-----------------

      do k = 1, N
        do j = JR_RANGE
          do i = IR_RANGE
            t(i,j,k,1,itemp) = T0 - FRACTS/Tcoef*rho0/g*rhobk(k)
#ifdef SALINITY
            t(i,j,k,1,isalt) = S0 + (1-FRACTS)/Scoef*rho0/g*rhobk(k)
#endif
          end do
        end do
      end do

! background rho 
!---------------

      rho = 0.0
      call rho_eos(tile)

      rho_bk_sav = rho(1,1,:)   !!! invariant sur l'horizontal

      do k = 2, N-1                            
         Nfreq2_bk(k) = -g/rho0*   
     &                         0.5*(rho(1,1,k+1)-rho(1,1,k-1))/Hz(1,1,k)
      end do

      Nfreq2_bk(N) = Nfreq2_bk(N-1)
      Nfreq2_bk(1) = Nfreq2_bk(2)

           
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% fin du background en z

!!! rhopr = vortex density profile/rho0*g
!----------------------------------------
           if (icase.eq.2) then
              if (linear_balance) then   
                if (mynode == 0) 
     &             write(STDOUT,*) 'initialisation linear thermal wind balance'
              else
                if (mynode == 0) 
     &             write(STDOUT,*) 'initialisation gradient wind balance'
             endif    ! linear balance
           endif   ! icase = 2
 
      do k = 1, N                           
        do j = JR_RANGE
          do i = IR_RANGE
            x = xr(i,j)
            y = yr(i,j)
            z = z_r(i,j,k)

            d2 = ((x-x1)/thick_y)**2 +((y-y1)/thick_y)**2 

            if (icase.gt.100) then    
              z2 = ((z-z1)/thick_z)**2 
              distance2    = d2 + Nfreq2_bk(k)/Nfreq2*z2
              denom(i,j,k) = 1 + distance2 + distance2**2 
              dNfreq_bk    = 0.5*(sqrt(Nfreq2_bk(k+1))-sqrt(Nfreq2_bk(k-1)))
     &                       /Hz(i,j,k)
              Fp           = -Rossby/4*(1+2*distance2)/denom(i,j,k)**(power+1)
              rhopr(i,j,k) =  2.*(z-z1) *(dNfreq_bk*(z-z1)+sqrt(Nfreq2_bk(k)))
     &                        /Nfreq * Fp*(1+2*Fp)*Nfreq2/(Burger)

           elseif (icase.eq.1) then
              z2 = ((z-z1)/thick_z)**2 
              distance2    = d2 + z2
              denom(i,j,k) = 1 + distance2 + distance2**2 
              rhopr(i,j,k) = -Rossby/2.*Nfreq2/Burger   
     &                       *(z-z1)*(1+2*distance2)/denom(i,j,k)**(power+1)
           elseif (icase.eq.2) then
              z2 = ((z-z1)/thick_z)**2 
              distance2    = d2 + z2
              if (linear_balance) then   
              rhopr(i,j,k) = -Rossby/2.*Nfreq2/Burger   
     &                     *(z-z1)*exp(-distance2)
              else
              rhopr(i,j,k) = -Rossby/2.*Nfreq2/Burger   
     &                     *(z-z1)*exp(-distance2)
     &                     *(1+Rossby/2.*exp(-distance2))
              endif    ! linear balance
           elseif (icase.eq.50) then        ! Go meddy data
              rhopr(i,j,k) = rhoin(k)*exp(-d2)
           endif
         end do
       end do
      end do

      valmax = nr_maxofall(maxval(abs(rhopr)))
      if (mynode == 0) print *,'-------------------------------------------'
      if (mynode == 0) print *,'max rhopr =', valmax
      if (mynode == 0) print *,'-------------------------------------------'


! vortex density profile : t,s,rho
!---------------------------------

      t(:,:,:,1,itemp) = T0 - FRACTS/Tcoef*rho0/g*rhopr
#ifdef SALINITY
      t(:,:,:,1,isalt) = S0 + (1-FRACTS)/Scoef*rho0/g*rhopr
#endif

      call rho_eos(tile)


! stratification totale: N2_bk+N2_pr
!-----------------------------------
      do k = 2, N-1                            
        do j = JR_RANGE
          do i = IR_RANGE
            Nfreq2_tot(i,j,k) = Nfreq2_bk(k)-g/rho0*   
     &                           0.5d0*(rho(i,j,k+1)-rho(i,j,k-1))/Hz(i,j,k)
          end do
        end do
      end do
      Nfreq2_tot(:,:,N) = Nfreq2_tot(:,:,N-1)
      Nfreq2_tot(:,:,1) = Nfreq2_tot(:,:,2)

      rho = 0.0
      rho(:,:,1) = rho_bk_sav(1)

! Calcul rho 3d
!--------------
      do k = 2, N                            
        do j = JR_RANGE
          do i = IR_RANGE
            rho(i,j,k) = rho(i,j,k-1) - rho0/g*0.5*  
     &                        (Nfreq2_tot(i,j,k)+Nfreq2_tot(i,j,k-1))*Hz(i,j,k)
          end do
        end do
      end do
!      rho(:,:,1) = rho(:,:,2)

! inversion de T et S
!---------------------
      t(:,:,:,1,itemp) = T0-FRACTS/Tcoef*(rho-rho0)
#ifdef SALINITY
      t(:,:,:,1,isalt) = S0+(1-FRACTS)/Scoef*(rho-rho0)
#endif


!!! + la partie compensée !!!
!----------------------------
      do k = 1, N                            
        do j = JR_RANGE
          do i = IR_RANGE
            x = xr(i,j)
            y = yr(i,j)
            z = z_r(i,j,k)
            d2 = ((x-x1)/thick_y)**2 +((y-y1)/thick_y)**2 
            if (icase.gt.100) then
              gaussian(i,j,k)=  (-3.5d-3/4.*log(-z/1.5+10.)+4.d-3) +  
     &           exp(-d2)*((-tanh((z+670)/200)+tanh((z+1520)/100))*1.d-3)
            elseif (icase.eq.50) then
              gaussian(i,j,k)=  (-3.0e-3/4.*log(-z/1.5+10.)+4.d-3) +  
     &           exp(-d2)*((-tanh((z+700)/200)+tanh((z+1520)/100))*0.8e-3)
            else
              gaussian(i,j,k)=  (0.8d-6*z-4d-4) + 
     &           exp(-d2)*(exp(-((z-z1)/450)**4)*1.5d-3)
            endif
          end do
        end do
      end do

      t(:,:,:,1,itemp) = t(:,:,:,1,itemp)  + (FRACTS/Tcoef*rho0)*gaussian
#ifdef SALINITY
      t(:,:,:,1,isalt) = t(:,:,:,1,isalt) + (FRACTS/Scoef*rho0)*gaussian
#endif

      valmax = nr_maxofall(maxval(abs(t(:,:,:,1,itemp))))
      if (mynode == 0)  print *,'max t =', valmax
      valmax = nr_maxofall(maxval(abs(t(:,:,:,1,isalt))))
      if (mynode == 0)  print *,'max s =', valmax

      valmax = nr_maxofall(maxval(rho))
      if (mynode==0) print *,'max potend_inittracer', valmax
      valmin = nr_minofall(minval(rho))
      if (mynode==0) print *,'min potend_inittracer', valmin


!....................................................................
!============================= 
!  	random perturbations 
!=============================   

        if (random_init) then  
#include "nr_init_random.h"
        endif


        if (nmdef.gt.1) then
! following option of balanced density perturbations IS included in random_pert_urot
          if (mynode == 0) then
            write(STDOUT,*) '--------------------------------------------------'
            write(STDOUT,*)
     &     ' random balanced perturbations added to density field '
            write(STDOUT,*) '--------------------------------------------------'
          endif
      
! adding to temperature field

          do k = 1, N                 !loop excluding surface/bottom boundaries
            do j = JR_RANGE
              do i = IR_RANGE
                f0 = (f(i,j)+f(i,j+1))*0.5d0
                t(i,j,k,1,itemp)= t(i,j,k,1,itemp)* facsolpert 
     &                   -0.5*25.0/4.0*f0*(rho0/g)* 
!    &                    -0.5*25.0/4.0*rho0*Rossby*Nfreq2/g* 
     &                  (rhopr_rand(i,j,k+1) -rhopr_rand(i,j,k-1))/(Hz(1,1,1)*2)
              end do
            end do
          end do
          if (mynode == 0)  
     &      print *,'random_pert_urot max t ', maxval(abs(t(:,:,:,1,itemp)))

! adding to salinity field
          do k = 1, N                  !loop excluding surface/bottom boundaries
            do j = JR_RANGE
              do i = IR_RANGE
                f0 = (f(i,j)+f(i,j+1))*0.5d0
                t(i,j,k,1,isalt)= t (i,j,k,1,isalt)* facsolpert 
     &                   +(1-0.5)*5.0/4.0*f0*(rho0/g)* 
!    &                    +(1-0.5)*5.0/4.0*rho0*Rossby*Nfreq2/g* 
     &                 (rhopr_rand(i,j,k+1) -rhopr_rand(i,j,k-1))/(Hz(i,j,k)*2)
              end do
            end do
          end do
          if (mynode == 0)  
     &      print *,'random_pert_urot max s ',maxval(abs(t(:,:,:,1,isalt)))
        endif ! if nmdef>1
         

!     call mpi_barrier(MPI_COMM_WORLD,ierr)  


#if defined EW_PERIODIC || defined NS_PERIODIC || defined  MPI
# ifdef SOLVE3D
      do itrc=1,NT
        if (.not.got_tini(itrc)) then

#  ifdef THREE_GHOST_POINTS_TS
          call exchange_r3d_3pts_tile (Istr,Iend,Jstr,Jend,
     &                            t(START_2D_ARRAY,1,1,itrc))
#  else
          call exchange_r3d_tile (Istr,Iend,Jstr,Jend,
     &                            t(START_2D_ARRAY,1,1,itrc))
#  endif

        endif
      enddo

# endif
#endif

!
!  Initial conditions for momentum components [m/s].
!--------------------------------------------------------------------
!



      ph=0.
      pres_up = 0.0
      u_geo = 0.
      v_geo = 0.

      if (linear_balance) then   
      !!!! geostrophic component
      !-------------------------

	      call rho_eos(tile)

	      do k = 1, N
		do j = JstrR-1,JendR+1
		  do i = IstrR-1 ,IendR+1
		    buoyancy = (rho(i,j,k))*g*(-1.0d0)/rho0
		    grd_z = buoyancy
		    pres_dif = Hz(i,j,k)*(grd_z)       ! valid under uniform grid           
		    ph(i,j,k) = (pres_up(i,j) + 0.5d0*pres_dif)
		    pres_up(i,j) = pres_up(i,j) + pres_dif
		  enddo
		enddo
	      enddo

       else

	      do k = 1, N
		do j = JstrR-1,JendR+1
		  do i = IstrR-1 ,IendR+1
                      x = xr(i,j)
                      y = yr(i,j)
                      d2 = ((x-x1)/thick_y)**2 +((y-y1)/thick_y)**2
                      z = z_r(i,j,k)
                      z2 = ((z-z1)/thick_z)**2
                      distance2    = d2 + z2
                      f0 = f(i,j)
                      ph(i,j,k) = - f0*thick_y**2*Rossby/4.*exp(-distance2)
		  end do
		end do
	      end do

       endif   ! linear_balance

	      do k = 1, N
		do j = 1, Mm
		  do i = 1, Lm
                    if (linear_balance) then   
		        f0 = (f(i,j)+f(i+1,j))*0.5d0
                    else
                        f0 = 1.
                    endif
		    cnt = 4.
		    u_geo(i,j,k) = ((ph(i-1,j+1,k)-ph(i-1,j,k)) + 
     &       	           (ph(i-1,j,k)-ph(i-1,j-1,k)) +   
     &         	         (ph(i,j+1,k)-ph(i,j,k)) +       
     &         	         (ph(i,j,k)-ph(i,j-1,k))           
     &         	         )*(-1.0d0)*pn(i,j)/(f0*cnt)
		  end do
		end do
	      end do

	      do k = 1, N
		do j = 1, Mm
		  do i = 1, Lm
                    if (linear_balance) then   
		        f0 = (f(i,j)+f(i,j+1))*0.5d0
                    else
                        f0 = 1.
                    endif
		    cnt = 4.
		    v_geo(i,j,k) = ((ph(i+1,j-1,k)-ph(i,j-1,k)) + 
     &  	                (ph(i,j-1,k)-ph(i-1,j-1,k)) +    
     &  	                (ph(i+1,j,k)-ph(i,j,k)) +       
     &  	                (ph(i,j,k)-ph(i-1,j,k))            
     & 		                 )*pm(i,j)/(f0*cnt)
		  end do
		end do
	      end do

	      u(:,:,:,1) = u_geo
	      v(:,:,:,1) = v_geo


      valmax = nr_maxofall(maxval(abs(u(:,:,:,1))))
      if (mynode == 0) print *,'max u =', valmax
      valmax = nr_maxofall(maxval(abs(v(:,:,:,1))))
      if (mynode == 0) print *,'max v =', valmax
!.....................................................................

!....................................................................
!============================= 
!  	random perturbations 
!=============================  
        if (random_init) then  
!#include "nr_init_random.h"

!...........
      do k = 1, N
        do j = 1, Mm
          do i = 1, Lm
            cnt = 4.
            u(i,j,k,1) = u(i,j,k,1) * facsolpert + 
     &            ((rhopr_rand(i-1,j+1,k)-rhopr_rand(i-1,j,k)) + 
     &             (rhopr_rand(i-1,j,k)-rhopr_rand(i-1,j-1,k)) +   
     &             (rhopr_rand(i,j+1,k)-rhopr_rand(i,j,k)) +       
     &             (rhopr_rand(i,j,k)-rhopr_rand(i,j-1,k))           
     &             )*(-1.0d0)*pn(i,j)/(cnt)

          end do
        end do
      end do
      if (mynode == 0)  
     &   print *,'random_pert_urot max sol_nu ', maxval(abs(u(:,:,:,1)))

      do k = 1, N
        do j = 1, Mm
          do i = 1, Lm
            cnt = 4.
            v(i,j,k,1) = v(i,j,k,1) * facsolpert + 
     &            ((rhopr_rand(i+1,j-1,k)-rhopr_rand(i,j-1,k)) + 
     &             (rhopr_rand(i,j-1,k)-rhopr_rand(i-1,j-1,k)) +    
     &             (rhopr_rand(i+1,j,k)-rhopr_rand(i,j,k)) +       
     &             (rhopr_rand(i,j,k)-rhopr_rand(i-1,j,k))            
     &            )*pm(i,j)/(cnt)
          end do
        end do
      end do
      if (mynode == 0)  
     &   print *,'random_pert_urot max sol_nv ', maxval(abs(v(:,:,:,1)))
     
        endif  ! init_random

!  Initial conditions for 2D momentum components.
!  barotropic (depth-averaged) velocity

      do k = 1, N
        do j = 1, Mm
          do i = 1, Lm
            ubar(i,j,1) = ubar(i,j,1) + u(i,j,k,1)*Hz(i,j,k)
          enddo
        enddo
      enddo
      do j = 1, Mm
        do i = 1, Lm
          ubar(i,j,1) = ubar(i,j,1)*hinv(i,j)
        enddo
      enddo       

      do k = 1, N
        do j = 1, Mm
          do i = 1, Lm
            vbar(i,j,1) = vbar(i,j,1) + v(i,j,k,1)*Hz(i,j,k)
          enddo
        enddo
      enddo
      do j = 1, Mm
        do i = 1, Lm
          vbar(i,j,1) = vbar(i,j,1)*hinv(i,j)
        enddo
      enddo
      
      endif  ! nrrec     
            
# else
      ERROR: ANA_INITIAL: no initials for ubar and vbar
# endif /* GO_MEDDY */

!
!--------------------------------------------------------------------
!  Exchange boundary information
!--------------------------------------------------------------------
!
#if defined EW_PERIODIC || defined NS_PERIODIC || defined  MPI
      call exchange_r2d_tile (Istr,Iend,Jstr,Jend,
     &                        zeta(START_2D_ARRAY,1))
      call exchange_u2d_tile (Istr,Iend,Jstr,Jend,
     &                        ubar(START_2D_ARRAY,1))
      call exchange_v2d_tile (Istr,Iend,Jstr,Jend,
     &                        vbar(START_2D_ARRAY,1))
#  ifdef SOLVE3D
#   ifdef THREE_GHOST_POINTS_UV
      call exchange_u3d_3pts_tile (Istr,Iend,Jstr,Jend,
     &                        u(START_2D_ARRAY,1,1))
      call exchange_v3d_3pts_tile (Istr,Iend,Jstr,Jend,
     &                        v(START_2D_ARRAY,1,1))
#   else
      call exchange_u3d_tile (Istr,Iend,Jstr,Jend,
     &                        u(START_2D_ARRAY,1,1))
      call exchange_v3d_tile (Istr,Iend,Jstr,Jend,
     &                        v(START_2D_ARRAY,1,1))
#   endif
#  endif /* SOLVE3D */

      ubar(:,:,2) = ubar(:,:,1)
      vbar(:,:,2) = vbar(:,:,1)
      u(:,:,:,2) = u(:,:,:,1)
      v(:,:,:,2) = v(:,:,:,1)
      t(:,:,:,2,:)=t(:,:,:,1,:)
#endif

#endif /* ANA_INITIAL */

#undef IR_RANGE
#undef IU_RANGE
#undef JR_RANGE
#undef JV_RANGE

      return
      end



