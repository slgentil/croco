      subroutine fillnan(datain,xr,yr,nx,ny)
      
      implicit none
      
!  Do a horizontal 2D cubic spline interpolation to fill Nan value of
!  input array
 
!  external 
      integer, intent(in) :: nx, ny
      real(kind=8), dimension(nx,ny), intent(inout) :: datain 	! sparse array with Nan
      real(kind=8), dimension(nx), intent(in) :: xr		! x coordinate
      real(kind=8), dimension(ny), intent(in) :: yr		! y coordinate
      
! Local      
!
      integer :: i,j,indx,nbofnan,ier
      integer, parameter :: ndim=2 	!  Dimension of the problem.
      integer :: ndata			!  The number of non-Nan points.
      integer, parameter :: n1=4, n2=4	!  Number of knots for splines
      integer, parameter :: ncf=n1*n2, nwrk=ncf*(ncf+3)
      integer, dimension(ndim) :: knots, nderiv      

      real(kind=8), dimension(:,:), allocatable :: xydata	! coordinates of defined values
      real(kind=8), dimension(:), allocatable :: vdata		! defined data
      real(kind=8), dimension(nwrk) :: work

      
!
!  Define error file, Fortran unit number, workstation type,
!  and workstation ID.
!
!     parameter (ierrf=6, lunit=2, iwtype=1, iwkid=1)

!
!  Find number of nan in input data 
!  Allocate and fill array with defined data only
!
      nbofnan = count(IsNan(datain(1:nx,1:ny)))
      
      ndata = (nx*ny) - nbofnan
      allocate(xydata(ndim,ndata))
      allocate(vdata(ndata))
      indx=0
      do j=1,ny
         do i=1,nx  
            if (.not. IsNaN(datain(i,j))) then
               indx = indx+1
               vdata(indx) = datain(i,j)
               xydata(1,indx) = xr(i)
               xydata(2,indx) = yr(j)
            endif
         enddo
      enddo

 
!  Specify the numbers of knots in each coordinate direction.
 
      knots(1) = n1
      knots(2) = n2
 
!  Interpolate sparse array
 
      call csa2d(ndata,xydata,vdata,knots,nx,ny,xr,yr,datain,nwrk,work,ier)

      end subroutine fillnan

