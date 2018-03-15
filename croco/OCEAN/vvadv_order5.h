!----------------------------------------------------------
! Compute vertical advective fluxes
! using 5th-order WENO scheme
!----------------------------------------------------------
!
          do k=3,N-3
            do i=Istr,Iend
              cff=0.5*(We(i  ,j,k)+We(i,j-1,k))
              FC(i,k)=cff*
     &              FLUX5(
     &             v(i,j,k-2,nrhs), v(i,j,k-1,nrhs),
     &             v(i,j,k  ,nrhs), v(i,j,k+1,nrhs),
     &             v(i,j,k+2,nrhs), v(i,j,k+3,nrhs), cff)
            enddo
          enddo
            do i=Istr,Iend
              cff=0.5*(We(i  ,j,2)+We(i,j-1,2))
              FC(i,2)=cff*
     &               FLUX3(
     &           v(i,j,1,nrhs), v(i,j,2,nrhs),
     &           v(i,j,3,nrhs), v(i,j,4,nrhs), cff)
            cff=0.5*(We(i  ,j,N-2)+We(i,j-1,N-2))
            FC(i,N-2)=cff*
     &               FLUX3(
     &           v(i,j,N-3,nrhs), v(i,j,N-2,nrhs),
     &           v(i,j,N-1,nrhs), v(i,j,N  ,nrhs),cff)

            FC(i,1  )=0.25*(We(i,j,  1)+We(i,j-1,  1))*( v(i,j,1 ,nrhs)
     &                               +  v(i,j,2  ,nrhs))
            FC(i,N-1)=0.25*(We(i,j,N-1)+We(i,j-1,N-1))*( v(i,j,N-1,nrhs)
     &                               +  v(i,j,N,  nrhs))

#  ifdef MOVING_BATHY
            FC(i,0)=0.5*v(i,j,1,nrhs)*
     &            (We(i,j  ,0)+
     &             We(i,j-1,0))
#  else
            FC(i,0)=0.
#  endif
            FC(i,N )=0.
          enddo
