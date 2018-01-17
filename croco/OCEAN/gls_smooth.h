! $Id: gls_smooth.h 1160 2013-06-11 09:48:32Z gcambon $
!
!======================================================================
! CROCO is a branch of ROMS developped at IRD and INRIA, in France
! The two other branches from UCLA (Shchepetkin et al)
! and Rutgers University (Arango et al) are under MIT/X style license.
! CROCO specific routines (nesting) are under CeCILL-C license.
!
! CROCO website : http://www.croco-ocean.org/
!======================================================================
!
! PART OF KPP2005 (Shchepetkin et al. 2005)
!
#   ifndef EW_PERIODIC
      if (WESTERN_EDGE) then
        do j=J_EXT_RANGE
          hwrk(Istr-1,j,k,nnew,ig)=hwrk(Istr,j,k,nnew,ig)
        enddo
      endif
      if (EASTERN_EDGE) then
        do j=J_EXT_RANGE
          hwrk(Iend+1,j,k,nnew,ig)=hwrk(Iend,j,k,nnew,ig)
        enddo
      endif
#   endif
#   ifndef NS_PERIODIC
      if (SOUTHERN_EDGE) then
        do i=I_EXT_RANGE
          hwrk(i,Jstr-1,k,nnew,ig)=hwrk(i,Jstr,k,nnew,ig)
        enddo
      endif
      if (NORTHERN_EDGE) then
        do i=I_EXT_RANGE
          hwrk(i,Jend+1,k,nnew,ig)=hwrk(i,Jend,k,nnew,ig)
        enddo
      endif
#    ifndef EW_PERIODIC
      if (WESTERN_EDGE.and.SOUTHERN_EDGE) then
        hwrk(Istr-1,Jstr-1,k,nnew,ig)=hwrk(Istr,Jstr,k,nnew,ig)
      endif
      if (WESTERN_EDGE.and.NORTHERN_EDGE) then
        hwrk(Istr-1,Jend+1,k,nnew,ig)=hwrk(Istr,Jend,k,nnew,ig)
      endif
      if (EASTERN_EDGE.and.SOUTHERN_EDGE) then
        hwrk(Iend+1,Jstr-1,k,nnew,ig)=hwrk(Iend,Jstr,k,nnew,ig)
      endif
      if (EASTERN_EDGE.and.NORTHERN_EDGE) then
        hwrk(Iend+1,Jend+1,k,nnew,ig)=hwrk(Iend,Jend,k,nnew,ig)
      endif
#    endif
#   endif

      do j=Jstr,Jend+1
        do i=Istr,Iend+1
          wrk(i,j)=0.25*( hwrk(i  ,j  ,k,nnew,ig)  
     &                  + hwrk(i-1,j  ,k,nnew,ig)
     &                  + hwrk(i  ,j-1,k,nnew,ig)
     &                  + hwrk(i-1,j-1,k,nnew,ig) )
#   ifdef MASKING
     &                   *pmask2(i,j)
#   endif
        enddo
      enddo
      do j=Jstr,Jend
        do i=Istr,Iend
#   ifdef MASKING
          cff=0.25*(pmask2(i,j)   +pmask2(i+1,j)
     &             +pmask2(i,j+1) +pmask2(i+1,j+1))
#   else
          cff=1.
#   endif
          hwrk(i,j,k,nnew,ig)=(1.-cff)*hwrk(i,j,k,nnew,ig)+
     &              0.25*(wrk(i,j)  +wrk(i+1,j)
     &                   +wrk(i,j+1)+wrk(i+1,j+1))
#   ifdef MASKING
          hwrk(i,j,k,nnew,ig)=hwrk(i,j,k,nnew,ig)*rmask(i,j)
#   endif
        enddo
      enddo
