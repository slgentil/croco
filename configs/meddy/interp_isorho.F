
#include "cppdefs.h"
  
      module isorho

      implicit none

      integer :: Ninterp, method
      integer :: rhomax_surf_index
      real, parameter :: Ninterp_fact = 2.
      
      real(kind=8) :: rhomax_surf, rhomin_bottom
      real(kind=8),dimension(:),allocatable :: rho_axis, rho_axis_interp
      real(kind=8),dimension(:,:,:),allocatable :: zpar_interp
      
      parameter(method=2)
      
      contains

!*************************************************************************
      subroutine init_isorho
      
# include "param.h"
# include "scalars.h"
# include "ocean3d.h"

      integer::i,j,k
      real(kind=8) :: rhomin, rhomax, nr_maxofall, nr_minofall, drho  
      
      rhomin = minval(rho(1:Lm,1:Mm,:))
      rhomin = nr_minofall(rhomin)
      rhomax = maxval(rho(1:Lm,1:Mm,:))
      rhomax = nr_maxofall(rhomax)
      Ninterp = floor(N*Ninterp_fact)
      allocate(rho_axis(N))
      allocate(rho_axis_interp(Ninterp))
      allocate(zpar_interp(GLOBAL_2D_ARRAY,Ninterp))
      drho = (rhomax - rhomin)/(Ninterp-1)
      do k=1,Ninterp
         rho_axis_interp(k) = rhomax - (k-1)*drho
      enddo

      end subroutine init_isorho
      
      
!*************************************************************************
     
      subroutine interp0d_isorho(str,val_in, val_out, i, j, k, 
     &            Nin , Nout)

! Interpolation d'un point sur la verticale en iso-densité

# include "param.h"
# include "scalars.h"
# include "ocean3d.h"
# include "grid.h"
# include "work.h"

      !external
      character(len=*), intent(in) :: str
      real(kind=8), dimension(Nin), intent(in) :: val_in
      real(kind=8), intent(out) :: val_out
      integer, intent(in) ::i,j,k, Nin, Nout

      select case(trim(str))

      case("u")
              rho_axis(:) = 0.5*(rho(i-1,j,:)+rho(i,j,:))
         
      case("v")
              rho_axis(:) = 0.5*(rho(i,j-1,:)+rho(i,j,:))
         
      case("rho")
              rho_axis(:) = rho(i,j,:)
         
      end select   
      call interp0d(val_in(:),rho_axis,val_out,rho_axis_interp(k),
     &              method,Nin )
     
      end subroutine interp0d_isorho        
!*********************************************************************************************

      subroutine interp0d(vari, yi, varo, yo, method, nyi)
      
    ! Interpolation along the axis (y)
    !
    ! - vari: input vector variable
    ! - varo: output scalar variable
    ! - yi: input y vector axis
    ! - yo: output y scalar axis
    ! - method: 0 = nearest, 1 = linear, 2 = cubic, 3 = hermit
    !

      implicit none

    ! Extrernal
      integer, intent(in) :: method, nyi
      real(kind=8), dimension(nyi), intent(in) :: vari, yi
      real(kind=8), intent(in) :: yo
      real(kind=8), intent(out) :: varo

    ! Internal
      integer :: iyi, iyo
      real(kind=8) :: dy0,dy1,yi0,yi1,mu
      real(kind=8) :: tension,bias,a0,a1,a2,a3 ! Hermit
      real(kind=8),dimension(nyi) :: zyi
      real(kind=8) :: vc0,vc1
      real(kind=8) :: mv, zero
      logical :: bmask(nyi)

    ! Treat only monotonically decreasing grids
      zyi = yi
      if (yi(nyi)>yi(1)) zyi = -yi

    ! Initialisation
      zero = 0.
      varo = zero / zero
      bias = 0.
      tension = 0.
      bmask = abs(vari-mv)<=abs(epsilon(0d0)*1.1*mv)

    ! Loop on input grid
      iyo = 1
      do iyi = 1, nyi-1

        yi0 = zyi(iyi)
        yi1 = zyi(iyi+1)

        if (yi1>yo) cycle
        if (yi0<yo) exit

        ! Distances to neighbours
        dy0 = yi0-yo
        dy1 = yo-yi1

        ! Interpolation
        if (dy0==0.) then
           varo = vari(iyi)
        else if (dy1==0.) then
           varo = vari(iyi+1)

        elseif (method==0) then

                ! Nearest neighbour
                if (dy0 < dy1) then
                    varo = vari(iyi)
                else
                    varo = vari(iyi+1)
                endif

        elseif (method==1)then

                ! Linear
                varo = (vari(iyi)*dy1 + vari(iyi+1)*dy0) / (dy0+dy1)

        else

                ! Cubic and Hermit
                !
                if (iyi==1)then ! y0
                    vc0 = 2*vari(iyi)-vari(iyi+1)
                else
                    vc0 = vari(iyi-1)
                endif
                if (iyi==nyi-1)then ! y3
                    vc1 = 2*vari(iyi+1)-vari(iyi)
                else
                    vc1 = vari(iyi+2)
                endif
                mu = dy0/(dy0+dy1)

                if (method==2)then

                    ! Cubic
                    !   mu2 = mu*mu;
                    !   a0 = y3 - y2 - y0 + y1;
                    !   a1 = y0 - y1 - a0;
                    !   a2 = y2 - y0;
                    !   a3 = y1;
                    !   return (a0*mu*mu2+a1*mu2+a2*mu+a3);

                    varo = vc1 - vari(iyi+1) - vc0 + vari(iyi) !a0
                    varo = mu**3*varo + mu**2*(vc0-vari(iyi)-varo) ! a0*mu^3 + a1*mu
                    varo = varo + mu*(vari(iyi+1)-vc0)
                    varo = varo + vari(iyi)

                else

                    ! Hermit
                    !   mu2 = mu * mu;
                    !   mu3 = mu2 * mu;
                    !   a0 =  2*mu3 - 3*mu2 + 1;
                    !   a1 =    mu3 - 2*mu2 + mu;
                    !   a2 =    mu3 -   mu2;
                    !   a3 = -2*mu3 + 3*mu2;
                    !   m0  = (y1-y0)*(1+bias)*(1-tension)/2;
                    !   m0 += (y2-y1)*(1-bias)*(1-tension)/2;
                    !   m1  = (y2-y1)*(1+bias)*(1-tension)/2;
                    !   m1 += (y3-y2)*(1-bias)*(1-tension)/2;
                    !   return(a0*y1 + a1*m0 + a2*m1 + a3*y2);
                    a0 = 2*mu**3 - 3*mu**2 + 1
                    a1 =    mu**3 - 2*mu**2 + mu
                    a2 =    mu**3 -   mu**2
                    a3 = -2*mu**3 + 3*mu**2
                    varo = a0*vari(iyi)
                    varo = varo + a1*(  
     &                   (vari(iyi)-vc0)        *(1+bias)*(1-tension)/2 +  
     &                   (vari(iyi+1)-vari(iyi))*(1-bias)*(1-tension)/2)
                    varo = varo + a2*(
     &                   (vari(iyi+1)-vari(iyi))*(1+bias)*(1-tension)/2 + 
     &                   (vc1-vari(iyi+1))      *(1-bias)*(1-tension)/2)
                    varo = varo + a3*vari(iyi+1)
                endif

                ! Mask
                varo = merge(mv, varo,  
     &                 any(bmask( max(iyi-1,1):min(iyi+2,nyi))) )
    
        endif
        iyo = iyo + 1

      end do
    
      end subroutine interp0d  
         
!*************************************************************************
     
      subroutine interp1d_isorho(str, val_in, val_out, i, j, Nin , Nout)

! Interpolation d'un vecteur sur la verticale en iso-densité

# include "param.h"
# include "scalars.h"
# include "ocean3d.h"
# include "grid.h"
# include "work.h"

      !external
      character(len=*), intent(in) :: str
      real(kind=8), dimension(Nin), intent(in) :: val_in
      real(kind=8), dimension(Nout), intent(out) :: val_out
      integer, intent(in) :: i,j, Nin, Nout
      
     ! for interpolation
      integer :: ier, nderiv
      integer, parameter :: nknots=4
      integer, parameter :: nwrk=nknots*(nknots+3)
      real(kind=8) :: ssmth , wts
      real(kind=8), dimension(nwrk) :: winterp
      
      ! local
      real(kind=8) :: zero,dnan
      
      zero=0. ; dnan=zero/zero
!     val_out = dnan
      
      select case(trim(str))

!---------------------------------------------------------------
! Momentum diagnostic variables.
!---------------------------------------------------------------

      case("u")
              rho_axis(:) = 0.5*(rho(i-1,j,:)+rho(i,j,:)) 
         
      case("v")
              rho_axis(:) = 0.5*(rho(i,j-1,:)+rho(i,j,:)) 
         
      case("rho")
              rho_axis(:) = rho(i,j,:)
         
      end select 

!     ssmth = 0.
!     nderiv = 0
!     wts = -1.
!     call csa1xd (Nin,rho_axis,val_in,wts,nknots,ssmth, 
!    &             nderiv,Nout,rho_axis_interp,val_out,
!    &                    nwrk,winterp,ier)     
      call interp1d(val_in(:),rho_axis,val_out(:),
     &                      rho_axis_interp,dnan,1,1,Nin,Nout,2)   

      end subroutine interp1d_isorho  
            
!*************************************************************************
     
      subroutine interp3d_isorho(str,val_in, val_out)

! Interpolation d'une matrice(i,j,k) sur la verticale en iso-densité

# include "param.h"
# include "scalars.h"
# include "ocean3d.h"
# include "grid.h"
# include "work.h"

      include 'mpif.h'

      ! external
      character(len=*), intent(in) :: str
      real(kind=8), dimension(GLOBAL_2D_ARRAY,N), 
     &                               intent(in) :: val_in
      real(kind=8), dimension(GLOBAL_2D_ARRAY,Ninterp), 
     &                              intent(out) :: val_out
      ! local
      integer::i,j,k
      real(kind=8) :: zero,dnan

      
     ! for interpolation
      integer :: ier, nderiv
      integer, parameter :: nknots=4
      integer, parameter :: nwrk=nknots*(nknots+3)
      real(kind=8) :: ssmth , wts
      real(kind=8), dimension(nwrk) :: winterp
      
      zero=0. ; dnan=zero/zero
!     val_out = dnan
      do j=1,Mm
         do i=1,Lm        
            select case(trim(str))
            case("u")
                rho_axis(:) = 0.5*(rho(i-1,j,:)+rho(i,j,:)) 
         
            case("v")
               rho_axis(:) = 0.5*(rho(i,j-1,:)+rho(i,j,:)) 
         
            case("rho")
               rho_axis(:) = rho(i,j,:)
         
            end select 
            
!           ssmth = 0.
!           nderiv = 0
!           wts = -1.
!           call csa1xd (N,rho_axis,val_in(i,j,:),wts,nknots,ssmth, 
!    &                    nderiv,Ninterp,rho_axis_interp,val_out(i,j,:),
!    &                    nwrk,winterp,ier) 
            call interp1d(val_in(i,j,:),rho_axis,val_out(i,j,:),
     &                    rho_axis_interp,dnan,method,1,N,Ninterp,0)
         enddo
      enddo
      
       
      end subroutine interp3d_isorho   
      
!*************************************************************************
      
      subroutine find_rho_index_down(value,index_down)
      
# include "param.h"


      integer::k, index_down
      real(kind=8), intent(in) ::  value 
      
! Loop on input grid

      index_down=1
      do k = 2, Ninterp
         if (value>rho_axis_interp(k)) exit
         index_down = k
      enddo   
      end subroutine find_rho_index_down

    
!*************************************************************************
      
      subroutine find_rho_index_up(value,index_up)
      
# include "param.h"


      integer::k, index_up
      real(kind=8), intent(in) ::  value 
      
! Loop on input grid

      index_up = Ninterp
      do k = Ninterp-1, 1, -1
         if (value<rho_axis_interp(k)) exit
         index_up = k         
      enddo   
      end subroutine find_rho_index_up     
      
      end module
