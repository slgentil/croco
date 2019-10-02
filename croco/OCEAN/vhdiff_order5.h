!
!===============================================================
!
! Compute diffusive part of UP5
!
!===============================================================
!
#  ifdef NS_PERIODIC
          jmin=0
          jmax=LOCALMM+1
#  else
#   ifdef MPI
          if (SOUTH_INTER) then
            jmin=0
          else
            jmin=3
          endif
          if (NORTH_INTER) then
            jmax=Mmmpi+1
          else
            jmax=Mmmpi-1
          endif
#   else
          jmin=3
          jmax=Mm-1
#   endif
#  endif
#  ifdef EW_PERIODIC
          imin=1
          imax=LOCALLM+1
#  else
#   ifdef MPI
          if (WEST_INTER) then
            imin=1
          else
            imin=3
          endif
          if (EAST_INTER) then
            imax=Lmmpi+1
          else
            imax=Lmmpi-1
          endif
#   else
          imin=3
          imax=Lm-1
#   endif
#  endif
!
!----------------------------------------------------------------------
!  j loop: VFe4
!----------------------------------------------------------------------
!
          DO j = JstrV-1,Jend  !j_loop_y_flux_5
                                                  !
            IF ( j.ge.jmin .and. j.le.jmax ) THEN ! use full stencil
                                                  !
              DO i = Istr,Iend
                vel = 0.5*(Hvom(i,j,k)+Hvom(i,j+1,k))
                flx5 = vel*flux6(
     &             v(i,j-2,k,nrhs), v(i,j-1,k,nrhs), 
     &             v(i,j  ,k,nrhs), v(i,j+1,k,nrhs),
     &             v(i,j+2,k,nrhs), v(i,j+3,k,nrhs), vel )
#  ifdef MASKING 
                flx3 = vel*flux4(
     &             v(i,j-1,k,nrhs), v(i,j  ,k,nrhs),
     &             v(i,j+1,k,nrhs), v(i,j+2,k,nrhs), vel ) 
                flx2 = vel*flux2(
     &             v(i,j  ,k,nrhs), v(i,j+1,k,nrhs), vel, cdif)
#   ifdef UP5_MASKING
                mask0=vmask(i,j)*vmask(i,j+1)
                mask2=vmask(i,j-1)*mask0*vmask(i,j+2)
                IF (vel.gt.0) THEN
                  mask1=vmask(i,j-1)*mask0
                  mask3=vmask(i,j-2)*mask2          
                ELSE
                  mask1=vmask(i,j+2)*mask0
                  mask3=vmask(i,j+3)*mask2
                ENDIF
                VFe4(i,j)=mask3*flx5+(1-mask3)*mask1*flx3+
     &                              (1-mask3)*(1-mask1)*mask0*flx2
#   else
                mask1=vmask(i,j-1)*vmask(i,j+2)
                mask2=vmask(i,j-2)*vmask(i,j+3)
                mask0=mask1*mask2
                VFe4(i,j)=mask0*flx5+(1-mask0)*mask1*flx3+
     &                              (1-mask0)*(1-mask1)*flx2
#   endif /* UP5_MASKING */
#  else
                VFe4(i,j)=flx5
#  endif /* MASKING */
              ENDDO

            ELSE IF ( j.eq.jmin-2 ) THEN   ! 2nd order flux next to south
                                           ! boundary
              DO i = Istr,Iend
                vel = 0.5*(Hvom(i,j,k)+Hvom(i,j+1,k))
                VFe4(i,j) = vel*flux2(
     &             v(i,j,k,nrhs), v(i,j+1,k,nrhs), vel, cdif)
              ENDDO
                                                             !
            ELSE IF ( j.eq.jmin-1 .and. jmax.ge.jmin ) THEN  ! 3rd of 4th order flux 2 in
                                                             ! from south boundary
              DO i = Istr,Iend
                vel = 0.5*(Hvom(i,j,k)+Hvom(i,j+1,k))
                flx3 = vel*flux4(
     &             v(i,j-1,k,nrhs), v(i,j  ,k,nrhs),
     &             v(i,j+1,k,nrhs), v(i,j+2,k,nrhs), vel )
#  ifdef MASKING
                flx2 = vel*flux2(
     &             v(i,j  ,k,nrhs), v(i,j+1,k,nrhs), vel, cdif)
                mask1=vmask(i,j-1)*vmask(i,j+2)
                VFe4(i,j)=mask1*flx3+(1-mask1)*flx2
#  else
                VFe4(i,j)=flx3
#  endif
              ENDDO
                                          !
            ELSE IF ( j.eq.jmax+2 ) THEN  ! 2nd order flux next to north
                                          ! boundary
              DO i = Istr,Iend
                vel = 0.5*(Hvom(i,j,k)+Hvom(i,j+1,k))
                VFe4(i,j) = vel*flux2(
     &             v(i,j,k,nrhs), v(i,j+1,k,nrhs), vel, cdif)
              ENDDO
                                          !
            ELSE IF ( j.eq.jmax+1 ) THEN  ! 3rd or 4th order flux 2 in from
                                          ! north boundary
              DO i = Istr,Iend
                vel = 0.5*(Hvom(i,j,k)+Hvom(i,j+1,k))
                flx3 = vel*flux4(
     &             v(i,j-1,k,nrhs), v(i,j  ,k,nrhs),
     &             v(i,j+1,k,nrhs), v(i,j+2,k,nrhs), vel )
#  ifdef MASKING
                flx2 = vel*flux2(
     &             v(i,j  ,k,nrhs), v(i,j+1,k,nrhs), vel, cdif)
                mask1=vmask(i,j-1)*vmask(i,j+2)
                VFe4(i,j)=mask1*flx3+(1-mask1)*flx2
#  else
                VFe4(i,j)=flx3
#  endif
              ENDDO
            ENDIF
          ENDDO ! j_loop_y_flux_5
!
!----------------------------------------------------------------------
!  i loop: VFx4
!----------------------------------------------------------------------
!
          DO i = Istr,Iend+1  !i_loop_x_flux_5
                                                  !
            IF ( i.ge.imin .and. i.le.imax ) THEN ! use full stencil
                                                  !
              DO j = JstrV,Jend
                vel = 0.5*(Huon(i,j-1,k)+Huon(i,j,k))
                flx5 = vel*flux6(
     &             v(i-3,j,k,nrhs), v(i-2,j,k,nrhs),
     &             v(i-1,j,k,nrhs), v(i  ,j,k,nrhs),
     &             v(i+1,j,k,nrhs), v(i+2,j,k,nrhs), vel )
#  ifdef MASKING
                flx3 = vel*flux4(
     &             v(i-2,j,k,nrhs), v(i-1,j,k,nrhs),
     &             v(i  ,j,k,nrhs), v(i+1,j,k,nrhs), vel )
                flx2 = vel*flux2(
     &             v(i-1,j,k,nrhs), v(i  ,j,k,nrhs), vel, cdif)
#   ifdef UP5_MASKING
                mask0=vmask(i-1,j)*vmask(i,j)
                mask2=vmask(i-2,j)*mask0*vmask(i+1,j)
                IF (vel.gt.0) THEN
                  mask1=vmask(i-2,j)*mask0
                  mask3=vmask(i-3,j)*mask2          
                ELSE
                  mask1=vmask(i+1,j)*mask0
                  mask3=vmask(i+2,j)*mask2
                ENDIF
                VFx4(i,j)=mask3*flx5+(1-mask3)*mask1*flx3+
     &                              (1-mask3)*(1-mask1)*mask0*flx2
#   else
                mask1=vmask(i-2,j)*vmask(i+1,j)
                mask2=vmask(i-3,j)*vmask(i+2,j)
                mask0=mask1*mask2
                VFx4(i,j)=mask0*flx5+(1-mask0)*mask1*flx3+
     &                              (1-mask0)*(1-mask1)*flx2
#   endif /* UP5_MASKING */
#  else
                VFx4(i,j)=flx5
#  endif /* MASKING */
              ENDDO
                                           !
            ELSE IF ( i.eq.imin-2 ) THEN   ! 2nd order flux next to south
                                           ! boundary
              DO j = JstrV,Jend
                vel = 0.5*(Huon(i,j-1,k)+Huon(i,j,k))
                VFx4(i,j) = vel*flux2(
     &             v(i-1,j,k,nrhs), v(i,j,k,nrhs), vel, cdif)
              ENDDO
                                                             !
            ELSE IF ( i.eq.imin-1 .and. imax.ge.imin ) THEN  ! 3rd of 4th order flux 2 in
                                                             ! from south boundary
              DO j = JstrV,Jend
                vel = 0.5*(Huon(i,j-1,k)+ Huon(i,j,k))
                flx3 = vel*flux4(
     &             v(i-2,j,k,nrhs), v(i-1,j,k,nrhs),
     &             v(i  ,j,k,nrhs), v(i+1,j,k,nrhs), vel )
#  ifdef MASKING
                flx2 = vel*flux2(
     &             v(i-1,j,k,nrhs), v(i  ,j,k,nrhs), vel, cdif)
                mask1=vmask(i-2,j)*vmask(i+1,j)
                VFx4(i,j)=mask1*flx3+(1-mask1)*flx2
#  else
                VFx4(i,j)=flx3
#  endif
              ENDDO
                                          !
            ELSE IF ( i.eq.imax+2 ) THEN  ! 2nd order flux next to north
                                          ! boundary
              DO j = JstrV,Jend
                vel = 0.5*(Huon(i,j-1,k)+ Huon(i,j,k))
                VFx4(i,j) = vel*flux2(
     &             v(i-1,j,k,nrhs), v(i,j,k,nrhs), vel, cdif)
              ENDDO
                                          !
            ELSE IF ( i.eq.imax+1 ) THEN  ! 3rd or 4th order flux 2 in from
                                          ! north boundary
              DO j = JstrV,Jend
                vel = 0.5*(Huon(i,j-1,k)+ Huon(i,j,k))
                flx3 = vel*flux4(
     &             v(i-2,j,k,nrhs), v(i-1,j,k,nrhs),
     &             v(i  ,j,k,nrhs), v(i+1,j,k,nrhs),  vel )
#  ifdef MASKING
                flx2 = vel*flux2(
     &             v(i-1,j,k,nrhs), v(i,j,k,nrhs), vel, cdif)
                mask1=vmask(i-2,j)*vmask(i+1,j)
                VFx4(i,j)=mask1*flx3+(1-mask1)*flx2
#  else
                VFx4(i,j)=flx3
#  endif
              ENDDO
            ENDIF
          ENDDO ! i_loop_x_flux_5

