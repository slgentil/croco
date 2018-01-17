#include "cppdefs.h"
#if defined NBQ && defined NBQ_IMPIJK
!
!======================================================================
!                      NBQ-Mode for NH-modeling
!                            Main Routine
!======================================================================
!
!> @note Main web documentation at:
! http://poc.obs-mip.fr/auclair/WOcean.fr/SNH/index_snh_home.htm
!
! DESCRIPTION: 
!
!> @brief SNBQ driver : Non-hydrostatic algorithm with the 
!                       Non-boussinesq solver.
!
!> @details Loops on the NBQ time step. See SNBQ web pages :
! REVISION HISTORY:
!
!> @authors
!> @date 2015 October
!> @todo
!
!======================================================================
!
      subroutine implicit_ijk_nbq ( icall )


!      use module_principal  
!      use module_parallele
       use module_nh 
       use module_nbq
!      use module_exp 
# ifdef TRACETXT
      use module_tracetxt_out
# endif
      implicit none

# include "param_F90.h"
# include "scalars_F90.h"
# include "ocean3d.h"
# include "grid.h"
# include "nbq.h"

# ifdef MPI
  include "mpif.h"
# endif

      integer :: i,j,k
      integer :: icall,flag_i,ierr_i,ifl_ascii,nd_i
      integer :: rank,couleur,clef

!.....Flag pour les sorties ASCII liees aux calculs implicites:
      ifl_ascii = 0

      if (icall.eq.0) then
!*******************************************************************************
! Initialisations
!
!*******************************************************************************

          theta_imp_nbq = 1. 

!.........Pre-process AMUB ==> nnzimp_nbq

          call amubdg(                                                &
                neqmom_nh(0)                                          &
               ,neqcont_nh                                            &
               ,neqmom_nh(0)                                          &
               ,momj_nh                                               &
               ,momi_nh                                               &
               ,contj_nh                                              &
               ,conti_nh                                              &
               ,iwork1_nbq                                            &
               ,nnzimp_nbq                                            &
               ,iwork2_nbq       )

      endif ! icall == 0

      if (icall.eq.1) then
!*******************************************************************************
!   Construction du systeme: I-alpha*G*D
!     + amub normal
!     + calcul structure au premier pas de temps 
!     + passage en mom * cont
!     + factorisation ou pre-conditionnement du solver
!*******************************************************************************

          flag_i = 1

!.........Process AMUB  ==> "imp" CRC-sparse matrix: alpha*G*D

          call amub3(                                                 &
                neqmom_nh(0)                                          &
               ,neqcont_nh                                            &
               ,flag_i                                                &
               ,dtnbq * ( soundspeed_nbq**2 * dtnbq + visc2_nbq )     &
                      * theta_imp_nbq                                 &
               ,momvg_nh(1:)                                          &  ! GGGG
               ,momj_nh                                               &
               ,momi_nh                                               &
               ,contv_nh(1:)                                          &
               ,contj_nh                                              &
               ,conti_nh                                              &
               ,impv_nbq(1:)                                          &
               ,impj_nbq                                              &
               ,impi_nbq                                              &
               ,nnzimp_nbq                                            &
               ,iwork2_nbq                                            &
               ,ierr_i      )

!........Computes "implicit" matrix: I + alpha*G*D

          call aplsca1(                                               &
                neqmom_nh(0)                                          &
               ,impv_nbq(1:)                                          &
               ,impj_nbq                                              &
               ,impi_nbq                                              &
               ,iwork1_nbq )

!         nnzimp_nbq = impi_nbq(neqmom_nh(0)+1)-1

      endif  ! icall == 1

      if (icall.eq.2) then
!*******************************************************************************
!   RHS 
!*******************************************************************************
       
!.......RHS = MOMxRHS(cont) 
        call amux(                                                  &
               neqmom_nh(0)                                         &
              ,rhp_nbq_a (1:neqcont_nh,rnrhs_nbq)                   &
              ,rhsimp_nbq(1)                                        &
              ,momvg_nh  (1)                                        &
              ,momj_nh   (1)                                        &
              ,momi_nh   (1)       )  

!.....Sauvegarde du second membre &
!     Ajouts des deux termes de l'equation QDM  du RHS: MOM * RHS(cont) + RHS(mom)
!     Equation de continuite:

        if (theta_imp_nbq.ne.1) then
          rhsimp_nbq    (1:neqmom_nh(0)) =                                         &
           qdm_nbq_a   (1:neqmom_nh(0),vnrhs_nbq)                                  &   
         + dtnbq*(  theta_imp_nbq*rhsimp_nbq  (1:neqmom_nh(0))*soundspeed_nbq**2   &
                   +(1.-theta_imp_nbq)*( soundspeed_nbq**2   * rhs1_nbq(l_nbq)     &
                                        - visc2_nbq_a(l_nbq) * rhsd2_nbq(l_nbq) )  &
                   + dqdmdt_nbq_a(1:neqmom_nh(0)) )
        else                        
          rhsimp_nbq    (1:neqmom_nh(0)) =                                         &
           qdm_nbq_a   (1:neqmom_nh(0),vnrhs_nbq)                                  &   
         + dtnbq*(  rhsimp_nbq  (1:neqmom_nh(0))*soundspeed_nbq**2                 &
                   + dqdmdt_nbq_a(1:neqmom_nh(0)) )
        endif
     
!*******************************************************************************
!   ASCII OUTPUT 
!*******************************************************************************
!         open(unit=10,file="matprod.dat")
!         open(unit=11,file="rhs.dat")
!         do l_nh=1,neqmom_nh(0)
!         write(11,*) rhsimp_nbq(l_nh)
!         do l1_nh=impi_nbq(l_nh),impi_nbq(l_nh+1)-1
!            write(10,*) l_nh,impj_nbq(l1_nh),impv_nbq(l1_nh)
!         enddo
!         enddo
!         close(10)
!         stop 'coucou'

!*******************************************************************************
!*******************************************************************************
!*******************************************************************************
! Resolution du systeme implicite
!*******************************************************************************
!*******************************************************************************
!*******************************************************************************


! XXX : resolution avec solver super-mega efficace!!!

#ifdef MPI
! !.....Initialization: MPI
!       CALL MPI_INIT(ierr_i)
	call MPI_COMM_RANK( MPI_COMM_WORLD,rank,ierr_i)

!.... Define a communicator for the package.
!     MUMPS in sequential need a communcitor with one process
      if ( (iic ==1) .and. (iif ==1) .and. iteration_nbq==1 ) then
	couleur = rank
	clef = -1 
	call MPI_COMM_SPLIT( MPI_COMM_WORLD,couleur,clef,mumps_comm,ierr_i)
	mumps_par%COMM = mumps_comm
      endif
#endif

      if (iif==1 .and. iic==1 .and. iteration_nbq==1 ) then

!     Initialize an instance of the package
!     for L U factorization (sym = 0, with working host)

      mumps_par%JOB = -1
      mumps_par%SYM = 0
      mumps_par%PAR = 1

      CALL DMUMPS(mumps_par)

      IF (mumps_par%INFOG(1).LT.0) THEN
       WRITE(6,'(A,A,I6,A,I9)') " ERROR RETURN: ",             &
                 "  mumps_par%INFOG(1)= ", mumps_par%INFOG(1), &
                 "  mumps_par%INFOG(2)= ", mumps_par%INFOG(2)
       stop 'err1'
       GOTO 500
      END IF
!      endif
!  Define problem on the host (processor 0)
! 	print *,iif,iic

	mumps_par%N  = neqmom_nh(0)
        mumps_par%NZ = impi_nbq(neqmom_nh(0)+1)-1 
        mumps_par%ICNTL(2) = 0 
        mumps_par%ICNTL(3) = 0 

        ALLOCATE( mumps_par%IRN ( mumps_par%NZ ) )
        ALLOCATE( mumps_par%JCN ( mumps_par%NZ ) )
        ALLOCATE( mumps_par%A   ( mumps_par%NZ ) )
        ALLOCATE( mumps_par%RHS ( mumps_par%N  ) )

        do l_nh  = 1 , neqmom_nh(0)
          do l1_nh = impi_nbq(l_nh),impi_nbq(l_nh+1)-1 
            mumps_par%IRN(l1_nh) = l_nh
            mumps_par%JCN(l1_nh) = impj_nbq(l1_nh) 
          enddo
        enddo  
	  ! Facto Symbolique
	mumps_par%JOB = 1
 	CALL DMUMPS(mumps_par)
      endif !if ( (iic == 1) .and (iif ==1))

!.....Matrix (impv_nbq) : external-mode time-step
      if ( iif==1 .and. iteration_nbq==1 ) then
        do l_nh  = 1 , neqmom_nh(0)
           do l1_nh = impi_nbq(l_nh),impi_nbq(l_nh+1)-1 
              mumps_par%A  (l1_nh) = impv_nbq(l1_nh) 
           enddo
        enddo  
     endif

      if (iif==1 .and. iteration_nbq==1 ) then
!     if (iif==1 .and. iic==1 .and. iteration_nbq==1 ) then
        ! Facto numerique
 	mumps_par%JOB = 2
 	CALL DMUMPS(mumps_par)

     endif

!....RHS (NBQ time-step)
     do l_nh  = 1 , neqmom_nh(0)
           mumps_par%RHS(l_nh) = rhsimp_nbq(l_nh)
      enddo  

!  Call package for solution
      mumps_par%JOB = 3
      CALL DMUMPS(mumps_par)

      IF (mumps_par%INFOG(1).LT.0) THEN
       WRITE(6,'(A,A,I6,A,I9)') " ERROR RETURN: ",             &
                 "  mumps_par%INFOG(1)= ", mumps_par%INFOG(1), &
                 "  mumps_par%INFOG(2)= ", mumps_par%INFOG(2)
       stop 'err2'
       GOTO 500
      END IF

!*******************************************************************************
!*******************************************************************************
!*******************************************************************************
!  Solver: end
!*******************************************************************************
!*******************************************************************************
!*******************************************************************************
!.....Recopie de la solution implicite pour w (update) et calcul rhs (sum)
      
        qdm_nbq_a(1:neqmom_nh(0),vnnew_nbq) = mumps_par%RHS(1:neqmom_nh(0))
        rhssum_nbq_a(1:neqmom_nh(0)) = rhssum_nbq_a(1:neqmom_nh(0))              &
             + ( qdm_nbq_a (1:neqmom_nh(0),vnnew_nbq)                            &
                -qdm_nbq_a (1:neqmom_nh(0),vnrhs_nbq)  ) / dtnbq                 &
             - dqdmdt_nbq_a(1:neqmom_nh(0)) 
 
!  Deallocate user data
!        IF ( mumps_par%MNID .eq. 0 )THEN
!          DEALLOCATE( mumps_par%IRN )
!          DEALLOCATE( mumps_par%JCN )
!          DEALLOCATE( mumps_par%A   )
!          DEALLOCATE( mumps_par%RHS )
!        END IF
! 
! 
! !     Destroy the instance (deallocate internal data structures)
!        mumps_par%JOB = -2
!        CALL DMUMPS(mumps_par)
!        IF (mumps_par%INFOG(1).LT.0) THEN
!         WRITE(6,'(A,A,I6,A,I9)') " ERROR RETURN: ",                &
!                   "  mumps_par%INFOG(1)= ", mumps_par%INFOG(1),    &
!                   "  mumps_par%INFOG(2)= ", mumps_par%INFOG(2)
!         stop 'err3'
!         GOTO 500
!        END IF

      endif
      return

!500  CALL MPI_FINALIZE(IERR)
 500  continue 
      stop 'passage 500'
      end subroutine implicit_ijk_nbq

#else
      subroutine implicit_ijk_nbq_empty
      end subroutine implicit_ijk_nbq_empty
#endif
