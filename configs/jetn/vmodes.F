#include "cppdefs.h"
!
! aponte jet (
#if defined VMODES
! aponte vmodes (
! start with a wrapper 
      subroutine vmodes(Istr,Iend,Jstr,Jend)

# include "set_obc_definitions.h"

      implicit none
# include "param.h"
# include "scalars.h"
# include "ocean3d.h"
!# include "ncscrum.h"
# include "private_scratch.h"

#ifdef MPI
      include 'mpif.h'
# include "mpi_cpl.h"
!      integer size, step, status(MPI_STATUS_SIZE), ierr
      integer status(MPI_STATUS_SIZE), ierr
      integer rank_root
      real buff(6*(N+1)*nvmodes)
#endif

! local variable
      integer Istr, Iend, Jstr, Jend
      real lNf(0:N)
      real lzc(1:N)
      real lzf(0:N)
!      real cvmodes(nvmodes)
      integer i,j,k,im
      character*20 filename
      character*20 rowfmt
! for pressure
      real ru(PRIVATE_2D_SCRATCH_ARRAY,N),
     &     rv(PRIVATE_2D_SCRATCH_ARRAY,N),
     &      P(PRIVATE_2D_SCRATCH_ARRAY,N),
     &     FC(PRIVATE_2D_SCRATCH_ARRAY),
     &    dZx(PRIVATE_2D_SCRATCH_ARRAY),
     &     rx(PRIVATE_2D_SCRATCH_ARRAY),
     &    dRx(PRIVATE_2D_SCRATCH_ARRAY)

#include "compute_auxiliary_bounds.h"


!      if (SOUTHERN_EDGE .and. EASTERN_EDGE) then
      if (SOUTHERN_EDGE) then
       i=Istr+1
       j=Jstr+1
      elseif (NORTHERN_EDGE) then
       i=Istr+1
       j=Jend
      else
       i=Istr+Lm/2
       j=Jstr+Mm/2
      endif
! store indices
      i_vmodes = i
      j_vmodes = j

      
      do k=1,N-1
! N is computed at k w-grid points 

        lNf(k) = sqrt( -g*(rho(i,j,k+1)-rho(i,j,k))/rho0
     &                  /(z_r(i,j,k+1)-z_r(i,j,k)) )
        zr_vmodes(k) = z_r(i,j,k)
        zw_vmodes(k) = z_w(i,j,k)
      enddo
      lNf(N)=lNf(N-1)
      lNf(0)=lNf(1)
      zw_vmodes(0) = z_w(i,j,0)
      zw_vmodes(N) = z_w(i,j,N)
      zr_vmodes(N) = z_r(i,j,N)
      call vmodes_core(lNf,zr_vmodes,zw_vmodes,c_vmodes,phiuv,phir)

! store reference stratification
      N2_vmodes(0) = lNf(0)**2
      do k=1,N
       T0_vmodes(k) = rho(i,j,k)
       N2_vmodes(k) = lNf(k)**2
      enddo

! computes scaling values for rho modes (for projections)
      do im=1,nvmodes
       scr_vmodes(im) = 0.0
       do k=1,N
        scr_vmodes(im) = scr_vmodes(im)
     &             + phir(k,im)**2 
     &               *0.5/(N2_vmodes(k)+N2_vmodes(k-1))
     &               *(zw_vmodes(k)-zw_vmodes(k-1))
       enddo
       scr_vmodes(im) = (zw_vmodes(N)-zw_vmodes(0))/scr_vmodes(im)
      enddo
       
! computes the vertical velocity modal structure
! as well as scaling values (for projections)
      do im=1,nvmodes
       scw_vmodes(im) = 0.0
!
       k=0
       phiw(k,im) = phir(k,im) /N2_vmodes(k)
       scw_vmodes(im) = phiw(k,im)**2 * N2_vmodes(k) 
     &              * (zr_vmodes(k+1)-zw_vmodes(k))
!
       do k=1,N-1
        phiw(k,im) = phir(k,im) /N2_vmodes(k)
        scw_vmodes(im) = scw_vmodes(im)
     &             + phiw(k,im)**2 * N2_vmodes(k)
     &               * (zr_vmodes(k+1)-zr_vmodes(k))
       enddo
!
       k=N
       phiw(k,im) = phir(k,im) /N2_vmodes(k)
       scw_vmodes(im) = scw_vmodes(im) 
     &            + phiw(k,im)**2 * N2_vmodes(k)
     &              * (zw_vmodes(k)-zr_vmodes(k))

       scw_vmodes(im) = (zw_vmodes(N)-zw_vmodes(0))/scw_vmodes(im)
      enddo

! get and store pressure at reference point
      call prsgrd_tile (Istr,Iend,Jstr,Jend, ru,rv, P,
     &                  A2d(1,3,0),A2d(1,4,0), FC,dZx,rx,dRx)
      do k=1,N
       P0_vmodes(k) = P(i_vmodes,j_vmodes,k)
      enddo

!      if (EASTERN_EDGE .and. NORTHERN_EDGE) then
!       write(*,*) 'c_vmodes = ', c_vmodes
!       write(*,*) 'phiuv ', phiuv(:,2)
!       write(*,*) 'phir ', phir(:,2)
!       write(*,*) 'phiw ', phiw(:,2)
!      endif

      if (mynode==0) then
        call var2cdf(zr_vmodes,size(shape(zr_vmodes)), 
     &                 shape(zr_vmodes),'vmodes_zr.nc')
        call var2cdf(zw_vmodes,size(shape(zw_vmodes)), 
     &                 shape(zw_vmodes),'vmodes_zw.nc')
        call var2cdf(c_vmodes,size(shape(c_vmodes)), 
     &                 shape(c_vmodes),'vmodes_c.nc')
        call var2cdf(phiuv,size(shape(phiuv)), 
     &                 shape(phiuv),'vmodes_phiuv.nc')
        call var2cdf(phir,size(shape(phir)), 
     &                 shape(phir),'vmodes_phir.nc')
        call var2cdf(phiw,size(shape(phiw)), 
     &                 shape(phiw),'vmodes_phiw.nc')
        call var2cdf(N2_vmodes,size(shape(N2_vmodes)), 
     &                 shape(N2_vmodes),'vmodes_N2.nc')
        call var2cdf(T0_vmodes,size(shape(T0_vmodes)), 
     &                 shape(T0_vmodes),'vmodes_T0.nc')
       endif
       
! write c and modes to ascii file
       if (.false.) then
          write (filename, "(A7,I0.3)") "vmodes.", mynode
          open(unit=12, file= trim(filename), action="write") 

          write(12,'(1X,e16.8)',advance='no') 0.0
          write(rowfmt,'(A,I4,A)') '(', N,'(1X,e16.8))'
          write(12,rowfmt,advance='no') (T0_vmodes(k), k=1,N)
          write(12,rowfmt,advance='no') (zr_vmodes(k), k=1,N)
          write(rowfmt,'(A,I4,A)') '(', N+1,'(1X,e16.8))'
          write(12,rowfmt) (zw_vmodes(k), k=0,N)

          do im=1,nvmodes
           write(12,'(1X,e16.8)',advance='no') c_vmodes(im)
           write(rowfmt,'(A,I4,A)') '(', N,'(1X,e16.8))'
           write(12,rowfmt,advance='no') (phiuv(k,im), k=1,N)
           write(12,rowfmt,advance='no') (phir(k,im), k=1,N)
           write(rowfmt,'(A,I4,A)') '(', N+1,'(1X,e16.8))'
           write(12,rowfmt) (phiw(k,im), k=0,N)
          enddo

          close(12)
       endif

#ifdef MPI
! The profile of one proc is used over the entire domain
      rank_root=0
!      if ( mynode .eq. rank_root ) then
!        buff( ...
!      endif
      call MPI_Bcast(c_vmodes, nvmodes, MPI_DOUBLE_PRECISION,
     &                       rank_root, ocean_grid_comm, ierr)
      call MPI_Bcast(phiuv, N*nvmodes, MPI_DOUBLE_PRECISION,
     &                       rank_root, ocean_grid_comm, ierr)
      call MPI_Bcast(phir, N*nvmodes, MPI_DOUBLE_PRECISION,
     &                       rank_root, ocean_grid_comm, ierr)
      call MPI_Bcast(phiw, (N+1)*nvmodes, MPI_DOUBLE_PRECISION,
     &                       rank_root, ocean_grid_comm, ierr)
      call MPI_Bcast(zr_vmodes, N, MPI_DOUBLE_PRECISION,
     &                       rank_root, ocean_grid_comm, ierr)
      call MPI_Bcast(zw_vmodes, N+1, MPI_DOUBLE_PRECISION,
     &                       rank_root, ocean_grid_comm, ierr)
      call MPI_Bcast(scr_vmodes, nvmodes, MPI_DOUBLE_PRECISION,
     &                       rank_root, ocean_grid_comm, ierr)
      call MPI_Bcast(scw_vmodes, nvmodes, MPI_DOUBLE_PRECISION,
     &                       rank_root, ocean_grid_comm, ierr)
      call MPI_Bcast(T0_vmodes, N, MPI_DOUBLE_PRECISION,
     &                       rank_root, ocean_grid_comm, ierr)
      call MPI_Bcast(N2_vmodes, N+1, MPI_DOUBLE_PRECISION,
     &                       rank_root, ocean_grid_comm, ierr)
!      call MPI_Bcast(P0_vmodes, N, MPI_DOUBLE_PRECISION,
!     &                       rank_root, ocean_grid_comm, ierr)
#endif     

!              call MPI_Send (buff,  14, MPI_DOUBLE_PRECISION,
!     &             mynode-step, 17, ocean_grid_comm,      ierr)
!            elseif (mynode .lt. size-step) then
!              call MPI_Recv (buff,  14, MPI_DOUBLE_PRECISION,
!     &             mynode+step, 17, ocean_grid_comm, status, ierr)


      MPI_master_only write(*,*) ' End of vmodes.F '


      return
      end


!=====================================================================
!      program vmodes
      subroutine vmodes_core(Nf,zc,zf,cr,phi,dphidz)
!=====================================================================
!
! assumes N is at cell edges (w levels)
!
! normalization is such that int(phi^2 dz) /h = 1
!
! Gill p160
! Kloosterziel 1995
! (phat of gill = phi of kloorsterziel)
!
! phi is the vertical structure function of the (u,v,p) and satisfies:
! d(dphi/dz*1/N^2)/dz + 1/c^2 * phi = 0
! g dphi/dz + N^2 phi = 0 at z = 0
! dphi/dz = 0 at z = -h
! c is the eigenvalue to be determined
!
! (u,v,p) = ( uh , vh , ph ) *phi
! (rho,w) = ( rhoh, wh/N(z)^2 ) * dphidz
!
! time dependent functions satisfy:
! duh/dt - f vh = -dph/dx /rho0
! dvh/dt + f uh = -dph/dy /rho0
! drhoh/dt = rho0/g * wh
! -wh/c^2 + duh/dx + dvh/dy = 0
!
! along with : ph = - g rhoh

!!$ C.Menesguen & S.Le Gentil (october 2013)
!!$ calculates eigen values and eigen modes for the current stratification
!!$ Nfreq2_bk     : profil de brunt-vaisala
!!$ nmde  : number of modes to calculate
!!$ vm_speed    : celerity en m/s
!!$ vm_baseu,vm_basev    : normalized eigen vectors
!!$

      implicit none  
!      include "mpif.h"

!!!! User input: you need to provide the number of vertical points and bounds
!!!! also contains the number of vertical modes nvmodes
#include "param.h"
#define nz N
! center cell first and last indices
#define ib 1
#define it N
! zf(i+ioff) is just below zc(i)
!#define ioff 0
#define ioff 0-1
! cell edges first and last indices
!#define ibf 0
!#define itf N
!!!!

! test, otherwise need to be in param.h
!      integer nz, ib, it, nvmodes
!      parameter (nz=100, ib=1, it=100, nvmodes=5)

! input
      real Nf(ib+ioff:it+1+ioff)
      real zc(ib:it) 
      real zf(ib+ioff:it) 

! choice of upper boundary condition
! bc_style = 0   is   w=0
! bc_style = 1   is   w=deta/dt
      integer bc_style
      parameter (bc_style=1)

! output
      real cr(nvmodes)
      real phi(ib:it,nvmodes)
      real dphidz(ib+ioff:it+1+ioff,nvmodes)

! local variables
      real lphi(ib:it,ib:it)
      real dzc(ib:it-1)
      real dzf(ib:it)
      real b(ib:it,ib:it)
      real g
      parameter (g=9.81)
      integer i,j,k,im
      integer imap(nvmodes)
      real sc(nvmodes)
! lapack stuff
      integer ilo,ihi,ierr,zl
      real wr(ib:it)
      real wi(ib:it)
      real scal(ib:it)
      real rconde(ib:it)
      real rcondv(ib:it)
      integer lwork, liwork
      parameter (lwork=3*nz, liwork=2*nz)
      integer iwork(liwork)
      real work(lwork)
      real abnrm, normv

! test
!      Nf(:)=1.D-3
!      do i=ib,it
!       Nf(i)=1.D-3
!       zf(i) = -4000.+40.*(i-1)
!       zc(i) = zf(i) + 20.
!      enddo
!      zf(it+1)=0.


! compute vertical grid spacing
      do i=ib,it-1
        dzc(i)=zc(i+1)-zc(i)
        dzf(i)=zf(i+1+ioff)-zf(i+ioff)
      enddo
      dzf(it)=zf(it+1+ioff)-zf(it+ioff)


!!! build b operator
      b(:,:)=0.
      do i=ib+1,it-1
       b(i,i+1) = -1./dzc(i)/dzf(i)/Nf(i+ioff+1)**2
       b(i,i) =  1./dzc(i)/dzf(i)/Nf(i+ioff+1)**2
     &          +1./dzc(i-1)/dzf(i)/Nf(i+ioff)**2
       b(i,i-1) = -1./dzc(i-1)/dzf(i)/Nf(i+ioff)**2
      enddo

! bottom boundary condition: dphi/dz=0 (w=0)
      i=ib
      b(i,i+1) = -1./dzc(i)/dzf(i)/Nf(i+1+ioff)**2
      b(i,i) = 1./dzc(i)/dzf(i)/Nf(i+1+ioff)**2
! top boundary condition
      i=it
      b(i,i) = 1/g/dzf(i)
     &        +1./dzc(i-1)/dzf(i)/Nf(i+ioff)**2
      b(i,i-1) = -1./dzc(i-1)/dzf(i)/Nf(i+ioff)**2

      if ( bc_style .eq. 0) then
       b(i,i) = 1./dzc(i-1)/dzf(i)/Nf(i+ioff)**2
      endif


!!! solve the boundary
!      call sgeevx('N','N','V','N',nz,b,nz,wr,wi,
      call dgeevx('N','N','V','N',nz,b,nz,wr,wi,
     &   zl,1,lphi,nz,ilo,ihi,scal,
     &   abnrm,rconde,rcondv,work,lwork,iwork,ierr)

      if (ierr.ne.0) then
       print *,'erreur dans dgeevx  ierr=',ierr
       stop
      endif

      do  i=ib,it
        wr(i)=1./sqrt(abs(wr(i)))
      enddo

!!! sort in descending eigenvalues
      im=1
      imap(im)=ib
      do i=ib,it
       if ( wr(imap(im)) .lt. wr(i) ) then 
        imap(im)=i
       endif
      enddo
      do im=2,nvmodes
        imap(im)=ib
        do i=ib,it
         if ( wr(imap(im)).lt.wr(i) .and. wr(i).lt.wr(imap(im-1)) ) then
          imap(im)=i
         endif
        enddo
      enddo

! store modes eigenvalues
      do im=1,nvmodes
       cr(im)=wr(imap(im))
      enddo

!      write(*,*) 'cvmode = ', cr

! store and scale vertical structures
! scaling is such that int_-h^0 phi^2 dz = h
      sc(:) = 0.
      do im=1,nvmodes
       do k=ib,it
        phi(k,im) = lphi(k,imap(im))
        sc(im) = sc(im) + phi(k,im)**2 *dzf(k)
       enddo

       sc(im) = sqrt(sc(im)/(zf(it+1+ioff)-zf(ib+ioff)))
       if ( phi(it,im) .lt. 0. )  sc(im)=-sc(im)

       do k=ib,it
        phi(k,im) = phi(k,im) / sc(im)
       enddo
      enddo

! compute the vertical derivative of the mode
      do im=1,nvmodes
       do k=ib+1,it
        dphidz(k+ioff,im) = (phi(k,im)-phi(k-1,im))/dzc(k-1)
       enddo
! top boundary
       k=it+1
       dphidz(k+ioff,im) = -Nf(k-1+ioff)**2/g*phi(k-1,im)
       if ( bc_style.eq.0 ) then
        dphidz(k+ioff,im) = 0.
       endif
! bottom boundary
       k=ib
       dphidz(k+ioff,im) = 0.
      enddo

!      write(*,*) 'phi(2)', phi(:,2)
!      write(*,*) 'dphidz(2)', dphidz(:,2)

      end subroutine vmodes_core
! aponte vmodes )
# endif /* VMODES */
! aponte jet )
