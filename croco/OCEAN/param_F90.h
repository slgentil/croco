! $Id: param.h 1619 2015-01-07 13:53:03Z marchesiello $
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
!
# define F90CODE
# include "param.h"

# ifdef BIOLOGY
#  ifdef PISCES
#   ifdef DIAGNOSTICS_BIO 
#    ifdef key_trc_diaadd
     integer  Nhi,Nco3,Naksp,Netot,Nprorca   &
     &          , Nprorcad,Npronew,Npronewd  &
     &          , Nprobsi,Nprofed,Nprofen    &
     &          , Ngrapoc,Ngrapoc2           &
     &          , Nmico2,Nmeso2              &
     &          , Nnitrifo2,Nfixo2,Nremino2  &
     &          , Npronewo2,Nprorego2        &
     &          , Nfld,Nflu16,Nkgco2,Natcco2,Nsinking &
     &          , Nsinkfer,Nsinksil,Nironsed        &
     &          , Nsinkcal,Nheup,Nnitrpot           &
     &          , Nirondep,Nsildep,Npo4dep          &
     &          , Nno3dep,Nnh4dep                   &
#    endif
     &          , NumFluxTerms,NumVSinkTerms,NumGasExcTerms
#   endif

#   ifdef key_trc_diaadd
       parameter (Nhi       = 1, &
     &            Nco3      = 2,&
     &            Naksp     = 3,&
     &            Netot     = 4,&
     &            Nprorca   = 5,&
     &            Nprorcad  = 6,&
     &            Npronew   = 7,&
     &            Npronewd  = 8,&
     &            Nprobsi   = 9,&
     &            Nprofen   = 10,&
     &            Nprofed   = 11,&
     &            Npronewo2 = 12,&
     &            Nprorego2 = 13,&
     &            Ngrapoc   = 14,&
     &            Ngrapoc2  = 15,&
     &            Nmico2    = 16,&
     &            Nmeso2    = 17,&
     &            Nnitrifo2 = 18,&
     &            Nremino2  = 19,&
     &            Nfixo2    = 20,&
     &            Nirondep  = 21,&
     &            Nironsed  = 22,&
     &            NumFluxTerms = Nironsed)

       parameter (Nfld      = 1,&
     &            Nflu16    = 2,&
     &            Nkgco2    = 3,&
     &            Natcco2   = 4,&
     &            Nsinking  = 5,&
     &            Nsinkfer  = 6,&
     &            Nsinksil  = 7,&
     &            Nsinkcal  = 8,&
     &            Nheup     = 9,&
     &            Nsildep   = 10,&
     &            Npo4dep   = 11,&
     &            Nno3dep   = 12,&
     &            Nnh4dep   = 13,&
     &            Nnitrpot  = 14,&
     &            NumGasExcTerms = 0,&
     &            NumVSinkTerms = Nnitrpot)
#   else
     integer NumFluxTerms,NumVSinkTerms,NumGasExcTerms
     parameter (NumFluxTerms = 0)
     parameter (NumGasExcTerms = 0, NumVSinkTerms = 0)
#   endif
#endif
#endif

# undef  F90CODE
 



