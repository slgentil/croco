! $Id: init_arrays.F 1554 2014-06-18 16:40:09Z ericm $
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

      subroutine init_arrays (tile)
      implicit none
#include "param.h"
      integer tile, trd
C$    integer omp_get_thread_num
#include "compute_tile_bounds.h"
      call init_arrays_tile (Istr,Iend,Jstr,Jend)
      return
      end

      subroutine init_arrays_tile (Istr,Iend,Jstr,Jend)
!
! Initialize (first touch) globally accessable arrays. Most of them
! are assigned to zeros, vertical mixing coefficients are assinged
! to their background values. These will remain unchenged if no
! vertical mixing scheme is applied. Because of the "first touch"
! default data distribution policy, this operation actually performs
! distribution of the shared arrays accross the cluster, unless
! another distribution policy is specified to override the default.
!
!     implicit none
#include "param.h"
#include "ncscrum.h"
#include "grid.h"
#include "ocean2d.h"
#include "ocean3d.h"
#include "coupling.h"
#include "averages.h"
#include "forces.h"
#include "mixing.h"
#include "climat.h"
#include "scalars.h"
#include "tides.h"
#include "sediment.h"
#include "bbl.h"
#include "diagnostics.h"
#include "boundary.h"
#include "sources.h"
#include "nbq.h"

#ifdef DIAGNOSTICS_VRT
# include "diags_vrt.h"
#endif

#ifdef DIAGNOSTICS_EK
# include "diags_ek.h"
#endif

#ifdef DIAGNOSTICS_PV
# include "diags_pv.h"
#endif

#ifdef DIAGNOSTICS_EDDY
# include "diags_eddy.h"
#endif

#if defined OUTPUTS_SURFACE && defined AVERAGES
# include "surf.h"
#endif
      integer Istr,Iend,Jstr,Jend, i,j
#ifdef SOLVE3D
     &       ,k,itrc
# ifdef DIAGNOSTICS_BIO
     &       ,iflux
# endif
#endif
#ifdef SEDIMENT
     &       ,ilay
#endif
#if defined SSH_TIDES || defined UV_TIDES
     &       ,itide
#endif
      real init    !!!!  0xFFFA5A5A ==> NaN
      parameter (init=0.)
#ifdef MPI
      include 'mpif.h'
      integer ierr
#endif
!
#include "compute_extended_bounds.h"
!
!  Initialize 2-D primitive variables.
!
      do j=JstrR,JendR
        do i=IstrR,IendR
          zeta(i,j,1)=0. ! init
          zeta(i,j,2)=init
          zeta(i,j,3)=init
          zeta(i,j,4)=init
!
          ubar(i,j,1)=init
          ubar(i,j,2)=init
          ubar(i,j,3)=init
          ubar(i,j,4)=init
!
          vbar(i,j,1)=init
          vbar(i,j,2)=init
          vbar(i,j,3)=init
          vbar(i,j,4)=init
#ifdef AVERAGES
          zeta_avg(i,j)=init 
          ubar_avg(i,j)=init
          vbar_avg(i,j)=init
#endif
#ifdef MRL_WCI
          sup(i,j)=init 
          ust2d(i,j)=init
          vst2d(i,j)=init
# ifdef AVERAGES
          sup_avg(i,j)=init 
          ust2d_avg(i,j)=init
          vst2d_avg(i,j)=init
# endif
#endif
#ifdef WET_DRY
          rmask_wet(i,j)=1
          umask_wet(i,j)=1
          vmask_wet(i,j)=1
          pmask_wet(i,j)=1
#endif
#ifdef SOLVE3D
          rufrc(i,j)=init
          rufrc(i,j)=init
!
          Zt_avg1(i,j)=0. !init
          DU_avg1(i,j,1)=0. !init
          DV_avg1(i,j,1)=0. !init
          DU_avg1(i,j,2)=0. !init
          DV_avg1(i,j,2)=0. !init
          DU_avg2(i,j)=0. !init
          DV_avg2(i,j)=0. !init
#endif

#if defined OUTPUTS_SURFACE && defined AVERAGES
          surft_avg(i,j)=0. ! init
          surfs_avg(i,j)=0. ! init
          surfz_avg(i,j)=0. ! init
          surfu_avg(i,j)=0. ! init
          surfv_avg(i,j)=0. ! init
#endif
        enddo
      enddo
           
#ifdef SOLVE3D
!
!  Initialize 3-D primitive variables.
!
! aponte jet (
# if defined JETN
! aponte itide (
#  if defined ITIDE
! aponte vmodes (
      do k=1,nvmodes
        do i=IstrR,IendR
          v_vmodes(i,k)=init
        enddo
      enddo
! aponte vmodes )
! aponte relax (
      nudgcof_itide = init
      nudgcof_zano_itide = init
! aponte relax )
      Atide_cos_mode = 0
      Atide_sin_mode = 0
! aponte itide )
#  endif
# endif /* JETN */
! aponte jet )

      do k=1,N
        do j=JstrR,JendR
          do i=IstrR,IendR
            u(i,j,k,1)=init
            u(i,j,k,2)=init
!
            v(i,j,k,1)=init
            v(i,j,k,2)=init
!
            rho(i,j,k) =init
# ifdef AVERAGES
            rho_avg(i,j,k)=init
            u_avg(i,j,k)=init
            v_avg(i,j,k)=init
            w_avg(i,j,k)=init
# endif
# ifdef MRL_WCI
            ust(i,j,k)=init
            vst(i,j,k)=init
            wst(i,j,k)=init
#  ifdef AVERAGES
            ust_avg(i,j,k)=init
            vst_avg(i,j,k)=init
            wst_avg(i,j,k)=init
#  endif
# endif
!
# if defined DIAGNOSTICS_UV || defined DIAGNOSTICS_EK_FULL
            MVadv(i,j,k,1)=init
            MVadv(i,j,k,2)=init
            MCor(i,j,k,1)=init
            MCor(i,j,k,2)=init
            MPrsgrd(i,j,k,1)=init
            MPrsgrd(i,j,k,2)=init
# if defined DIAGNOSTICS_BARO
            MBaro(i,j,k,1)=init
            MBaro(i,j,k,2)=init
# endif
# endif
# if defined DIAGNOSTICS_UV || defined DIAGNOSTICS_EK_FULL || defined DIAGNOSTICS_PV
            MXadv(i,j,k,1)=init ! indice 1 for u_momentum
            MXadv(i,j,k,2)=init ! indice 2 for v_momentum
            MYadv(i,j,k,1)=init ! ...
            MYadv(i,j,k,2)=init
            MHmix(i,j,k,1,1)=init
            MHmix(i,j,k,2,1)=init
            MHmix(i,j,k,1,2)=init
            MHmix(i,j,k,2,2)=init
            MVmix(i,j,k,1)=init
            MVmix(i,j,k,2)=init
            MHdiff(i,j,k,1)=init
            MHdiff(i,j,k,2)=init
# endif

# if defined DIAGNOSTICS_UV || defined DIAGNOSTICS_PV
            MVmix2(i,j,k,1)=init
            MVmix2(i,j,k,2)=init
            Mrate(i,j,k,1)=init
            Mrate(i,j,k,2)=init
# endif


# ifdef DIAGNOSTICS_UV
#  ifdef MRL_WCI
            Mvf(i,j,k,1)=init
            Mvf(i,j,k,2)=init
            Mbrk(i,j,k,1)=init
            Mbrk(i,j,k,2)=init
            MStCo(i,j,k,1)=init
            MStCo(i,j,k,2)=init
            MVvf(i,j,k,1)=init
            MVvf(i,j,k,2)=init
            MPrscrt(i,j,k,1)=init
            MPrscrt(i,j,k,2)=init
            Msbk(i,j,k,1)=init
            Msbk(i,j,k,2)=init
            Mbwf(i,j,k,1)=init
            Mbwf(i,j,k,2)=init
            Mfrc(i,j,k,1)=init
            Mfrc(i,j,k,2)=init
#  endif             
#  ifdef AVERAGES
            MXadv_avg(i,j,k,1)=init
            MXadv_avg(i,j,k,2)=init
            MYadv_avg(i,j,k,1)=init
            MYadv_avg(i,j,k,2)=init
            MVadv_avg(i,j,k,1)=init
            MVadv_avg(i,j,k,2)=init
            MCor_avg(i,j,k,1)=init
            MCor_avg(i,j,k,2)=init
            MPrsgrd_avg(i,j,k,1)=init
            MPrsgrd_avg(i,j,k,2)=init
# if defined DIAGNOSTICS_BARO
            MBaro_avg(i,j,k,1)=init
            MBaro_avg(i,j,k,2)=init
#  endif
            MHmix_avg(i,j,k,1)=init
            MHmix_avg(i,j,k,2)=init
            MVmix_avg(i,j,k,1)=init
            MVmix_avg(i,j,k,2)=init
            MVmix2_avg(i,j,k,1)=init
            MVmix2_avg(i,j,k,2)=init
            MHdiff_avg(i,j,k,1)=init
            MHdiff_avg(i,j,k,2)=init
            Mrate_avg(i,j,k,1)=init
            Mrate_avg(i,j,k,2)=init
#   ifdef MRL_WCI
            Mvf_avg(i,j,k,1)=init
            Mvf_avg(i,j,k,2)=init
            Mbrk_avg(i,j,k,1)=init
            Mbrk_avg(i,j,k,2)=init
            MStCo_avg(i,j,k,1)=init
            MStCo_avg(i,j,k,2)=init
            MVvf_avg(i,j,k,1)=init
            MVvf_avg(i,j,k,2)=init
            MPrscrt_avg(i,j,k,1)=init
            MPrscrt_avg(i,j,k,2)=init
            Msbk_avg(i,j,k,1)=init
            Msbk_avg(i,j,k,2)=init
            Mbwf_avg(i,j,k,1)=init
            Mbwf_avg(i,j,k,2)=init
            Mfrc_avg(i,j,k,1)=init
            Mfrc_avg(i,j,k,2)=init
#   endif             
#  endif /* AVERAGES */
# endif /* DIAGNOSTICS_UV */
          enddo
        enddo
      enddo
      do k=0,N
        do j=JstrR,JendR
          do i=IstrR,IendR
             We(i,j,k)=init
# ifdef VADV_ADAPT_IMP
             Wi(i,j,k)=init
# endif
# ifdef AVERAGES
            omega_avg(i,j,k)=init
# if defined ANA_VMIX || defined BVF_MIXING \
  || defined LMD_MIXING || defined LMD_SKPP || defined LMD_BKPP \
  || defined GLS_MIX2017 || defined GLS_MIXING
            bvf_avg(i,j,k)=init
# endif
# endif
          enddo
        enddo
      enddo
      do itrc=1,NT
        got_tini(itrc)=.false.
        do k=1,N
          do j=JstrR,JendR
            do i=IstrR,IendR
              t(i,j,k,1,itrc)=init
              t(i,j,k,2,itrc)=init
# ifdef AVERAGES
              t_avg(i,j,k,itrc)=init
# endif
            enddo
          enddo
        enddo
      enddo

# if defined DIAGNOSTICS_TS || defined DIAGNOSTICS_PV
      do itrc=1,NT
        do k=1,N
          do j=JstrR,JendR
            do i=IstrR,IendR
              THmix(i,j,k,itrc)=init
              TVmix(i,j,k,itrc)=init
              TForc(i,j,k,itrc)=init
              Trate(i,j,k,itrc)=init
            enddo
          enddo
        enddo
      enddo
# endif /*DIAGNOSTICS_TS || DIAGNOSTICS_PV*/

# ifdef DIAGNOSTICS_TS
      do itrc=1,NT
        do k=1,N
          do j=JstrR,JendR
            do i=IstrR,IendR
              TXadv(i,j,k,itrc)=init
              TYadv(i,j,k,itrc)=init
              TVadv(i,j,k,itrc)=init
              THmix(i,j,k,itrc)=init
              TVmix(i,j,k,itrc)=init
              TForc(i,j,k,itrc)=init
              Trate(i,j,k,itrc)=init
#  ifdef AVERAGES
              TXadv_avg(i,j,k,itrc)=init
              TYadv_avg(i,j,k,itrc)=init
              TVadv_avg(i,j,k,itrc)=init
              THmix_avg(i,j,k,itrc)=init
              TVmix_avg(i,j,k,itrc)=init
              TForc_avg(i,j,k,itrc)=init
              Trate_avg(i,j,k,itrc)=init
#  endif /* AVERAGES*/
            enddo
          enddo
        enddo
      enddo
#  ifdef DIAGNOSTICS_TS_MLD
      do itrc=1,NT
        do j=JstrR,JendR
          do i=IstrR,IendR
            TXadv_mld(i,j,itrc)=init
            TYadv_mld(i,j,itrc)=init
            TVadv_mld(i,j,itrc)=init
            THmix_mld(i,j,itrc)=init
            TVmix_mld(i,j,itrc)=init
            TForc_mld(i,j,itrc)=init
            Trate_mld(i,j,itrc)=init
            Tentr_mld(i,j,itrc)=init
#  ifdef AVERAGES
            TXadv_mld_avg(i,j,itrc)=init
            TYadv_mld_avg(i,j,itrc)=init
            TVadv_mld_avg(i,j,itrc)=init
            THmix_mld_avg(i,j,itrc)=init
            TVmix_mld_avg(i,j,itrc)=init
            TForc_mld_avg(i,j,itrc)=init
            Trate_mld_avg(i,j,itrc)=init
            Tentr_mld_avg(i,j,itrc)=init
#  endif /* AVERAGES*/


# ifdef DIAGNOSTICS_EDDY
#  ifdef AVERAGES
            eddyuu_avg(i,j,k)=init
            eddyvv_avg(i,j,k)=init
            eddyuv_avg(i,j,k)=init
            eddyub_avg(i,j,k)=init
            eddyvb_avg(i,j,k)=init
            eddywb_avg(i,j,k)=init
            eddyuw_avg(i,j,k)=init
            eddyvw_avg(i,j,k)=init
#  endif /* AVERAGES */
# endif /* DIAGNOSTICS_EDDY */


          enddo
        enddo
      enddo
#  endif
# endif /*DIAGNOSTICS_TS*/



# ifdef DIAGNOSTICS_VRT
        do j=jstrR,jendR             ! 3-D primitive
          do i=istrR,iendR           ! variables.

            vrtXadv(i,j)=init 
            vrtYadv(i,j)=init
            vrtHdiff(i,j)=init
            vrtCor(i,j)=init
            vrtPrsgrd(i,j)=init
            vrtHmix(i,j)=init
            vrtVmix(i,j)=init
            vrtrate(i,j)=init
            vrtVmix2(i,j)=init
            vrtWind(i,j)=init
            vrtDrag(i,j)=init
# if defined DIAGNOSTICS_BARO
            vrtBaro(i,j)=init
# endif
#  ifdef AVERAGES
            vrtXadv_avg(i,j)=init
            vrtYadv_avg(i,j)=init
            vrtHdiff_avg(i,j)=init
            vrtCor_avg(i,j)=init
            vrtPrsgrd_avg(i,j)=init
            vrtHmix_avg(i,j)=init
            vrtVmix_avg(i,j)=init
            vrtrate_avg(i,j)=init
            vrtVmix2_avg(i,j)=init
            vrtWind_avg(i,j)=init
            vrtDrag_avg(i,j)=init
# if defined DIAGNOSTICS_BARO
            vrtBaro_avg(i,j)=init
# endif
#  endif /* AVERAGES */
          enddo
        enddo
# endif /* DIAGNOSTICS_VRT */


# ifdef DIAGNOSTICS_EK
        do j=jstrR,jendR             ! 3-D primitive
          do i=istrR,iendR           ! variables.

            ekHadv(i,j)=init 
            ekHdiff(i,j)=init
            ekVadv(i,j)=init            
            ekCor(i,j)=init
            ekPrsgrd(i,j)=init
            ekHmix(i,j)=init
            ekVmix(i,j)=init
            ekrate(i,j)=init
            ekvol(i,j)=init
            ekVmix2(i,j)=init
            ekWind(i,j)=init
            ekDrag(i,j)=init
# if defined DIAGNOSTICS_BARO
            ekBaro(i,j)=init
# endif
#  ifdef AVERAGES
            ekHadv_avg(i,j)=init
            ekHdiff_avg(i,j)=init
            ekVadv_avg(i,j)=init           
            ekCor_avg(i,j)=init
            ekPrsgrd_avg(i,j)=init
            ekHmix_avg(i,j)=init
            ekVmix_avg(i,j)=init
            ekrate_avg(i,j)=init
            ekvol_avg(i,j)=init
            ekVmix2_avg(i,j)=init
            ekWind_avg(i,j)=init
            ekDrag_avg(i,j)=init
# if defined DIAGNOSTICS_BARO
            ekBaro_avg(i,j)=init
# endif
#  endif /* AVERAGES */

# ifdef DIAGNOSTICS_EK_MLD
            ekHadv_mld(i,j)=init 
            ekHdiff_mld(i,j)=init
            ekVadv_mld(i,j)=init            
            ekCor_mld(i,j)=init
            ekPrsgrd_mld(i,j)=init
            ekHmix_mld(i,j)=init
            ekVmix_mld(i,j)=init
            ekrate_mld(i,j)=init
            ekvol_mld(i,j)=init
            ekVmix2_mld(i,j)=init
# if defined DIAGNOSTICS_BARO
            ekBaro_mld(i,j)=init
# endif
#  ifdef AVERAGES
            ekHadv_mld_avg(i,j)=init
            ekHdiff_mld_avg(i,j)=init
            ekVadv_mld_avg(i,j)=init           
            ekCor_mld_avg(i,j)=init
            ekPrsgrd_mld_avg(i,j)=init
            ekHmix_mld_avg(i,j)=init
            ekVmix_mld_avg(i,j)=init
            ekrate_mld_avg(i,j)=init
            ekvol_mld_avg(i,j)=init
            ekVmix2_mld_avg(i,j)=init
# if defined DIAGNOSTICS_BARO
            ekBaro_mld_avg(i,j)=init
# endif
#  endif /* AVERAGES */
# endif /* DIAGNOSTICS_EK */

          enddo
        enddo
# endif /* DIAGNOSTICS_EK */

# ifdef DIAGNOSTICS_PV
# if defined DIAGNOSTICS_PV_FULL
      do k=1,N
        do j=jstrR,jendR             ! 3-D primitive
          do i=istrR,iendR           ! variables.
            pv(i,j,k)=init
            pvd(i,j,k)=init
#  ifdef AVERAGES
            pv_avg(i,j,k)=init
            pvd_avg(i,j,k)=init
#  endif /* AVERAGES */
          enddo
        enddo
      enddo
# endif
      do itrc=1,NT
        do k=1,N
          do j=JstrR,JendR
            do i=IstrR,IendR
              Trhs(i,j,k,itrc)=init
#  ifdef AVERAGES
              Trhs_avg(i,j,k,itrc)=init
#  endif /* AVERAGES*/
            enddo
          enddo
        enddo
      enddo

      do itrc=1,2
        do k=1,N
          do j=JstrR,JendR
            do i=IstrR,IendR
              Mrhs(i,j,k,itrc)=init
#  ifdef AVERAGES
              Mrhs_avg(i,j,k,itrc)=init
#  endif /* AVERAGES*/
            enddo
          enddo
        enddo
      enddo

# endif /* DIAGNOSTICS_PV */

# ifdef DIAGNOSTICS_BIO
      do k=1,N
        do j=JstrR,JendR
          do i=IstrR,IendR
#  ifdef PISCES
#   ifdef key_trc_diaadd
            do iflux=1,NumFluxTerms
              bioFlux(i,j,k,iflux) = init
            enddo
            do iflux=1,NumVSinkTerms
              bioVSink(i,j,iflux)= init
            enddo
#   endif
#  else
            do iflux=1,NumFluxTerms
              bioFlux(i,j,k,iflux) = init
            enddo
            do iflux=1,NumVSinkTerms
              bioVSink(i,j,k,iflux)= init
            enddo
#  endif
          enddo
        enddo
      enddo
#  if (defined BIO_NChlPZD && defined OXYGEN) || defined BIO_BioEBUS 
        do j=JstrR,JendR
          do i=IstrR,IendR
            do iflux=1,NumGasExcTerms
              GasExcFlux(i,j,iflux)= init
            enddo        
          enddo
        enddo
#  endif
# endif /*DIAGNOSTICS_BIO*/

# ifdef BIO_BioEBUS
          do j=JstrR,JendR
            do i=IstrR,IendR
               wind10(i,j) =init
#   ifdef AVERAGES
               wind10_avg(i,j) =init
#   endif     	    
            enddo
          enddo

        do k=1,N  
          do j=JstrR,JendR
            do i=IstrR,IendR
               AOU(i,j,k) =init
#   ifdef AVERAGES
               AOU_avg(i,j,k) =init
#   endif 
            enddo
          enddo
        enddo
# endif 

# ifdef SEDIMENT
!
!  Initialize sediment module arrays.
!
      do j=JstrR,JendR
        do i=IstrR,IendR
           Ssize(i,j)=init
           Sdens(i,j)=init
           taucb(i,j)=init
           w_set(i,j)=init
           Hripple(i,j)=init
           Lripple(i,j)=init
          do ilay=1,NLAY
            bed_thick(i,j,ilay)=init
            bed_poros(i,j,ilay)=init
            bed_age(i,j,ilay)=init
            do itrc=1,NST
              bed_frac(i,j,ilay,itrc)=init
              bed_mass(i,j,ilay,1,itrc)=init
              bed_mass(i,j,ilay,2,itrc)=init
            enddo
          enddo
        enddo
      enddo
# endif /* SEDIMENT */
# ifdef BBL
      do j=JstrR,JendR
        do i=IstrR,IendR
          Abed(i,j)=init
          Hripple(i,j)=init
          Lripple(i,j)=init
          Zbnot(i,j)=init
          Zbapp(i,j)=init
          Bustrc(i,j)=init
          Bvstrc(i,j)=init
          Bustrw(i,j)=init
          Bvstrw(i,j)=init
          Bustrcwmax(i,j)=init
          Bvstrcwmax(i,j)=init
#ifdef MOVING_BATHY 
          dh(i,j)=init
#endif
        enddo
      enddo
# endif /* BBL */
# ifdef MRL_WCI
      do k=0,N
        do j=JstrR,JendR
          do i=IstrR,IendR
            Akb(i,j,k)=init 
            Akw(i,j,k)=init 
            E_pre(i,j,k)=init
# ifdef AVERAGES
            Akb_avg(i,j,k)=init 
            Akw_avg(i,j,k)=init 
# endif
          enddo
        enddo
      enddo
      do k=1,N
        do j=JstrR,JendR
          do i=IstrR,IendR
            kvf(i,j,k)=init
# ifdef AVERAGES
            kvf_avg(i,j,k)=init
# endif
          enddo
        enddo
      enddo
      do j=JstrR,JendR
        do i=IstrR,IendR
          calP(i,j)=init 
          Kapsrf(i,j)=init 
# ifdef AVERAGES
          calP_avg(i,j)=init 
          Kaps_avg(i,j)=init 
# endif
          enddo
        enddo
# endif
#endif /* SOLVE3D */
!
!  Initialize forcing arrays (see "forces.h").
!
      do j=JstrR,JendR
        do i=IstrR,IendR
          sustr(i,j)=init
          svstr(i,j)=init
# ifdef OA_COUPLING
          smstr(i,j)=init
# endif
#ifndef ANA_SMFLUX
          sustrg(i,j,1)=init
          svstrg(i,j,1)=init
          sustrg(i,j,2)=init
          svstrg(i,j,2)=init
#endif
          bustr(i,j)=init
          bvstr(i,j)=init
#ifdef AVERAGES
          bostr_avg(i,j)=init
          wstr_avg(i,j)=init
          sustr_avg(i,j)=init
          svstr_avg(i,j)=init
#endif
#ifndef ANA_BMFLUX
          bustrg(i,j,1)=init
          bvstrg(i,j,1)=init
          bustrg(i,j,2)=init
          bvstrg(i,j,2)=init
#endif
        enddo
      enddo
#ifdef SOLVE3D
      do itrc=1,NT
        do j=JstrR,JendR
          do i=IstrR,IendR
            stflx(i,j,itrc)=init
# ifdef AVERAGES
            stflx_avg(i,j,itrc)=init
# endif
# if !defined ANA_STFLUX || !defined ANA_SSFLUX
            stflxg(i,j,1,itrc)=init
            stflxg(i,j,2,itrc)=init
# endif
            btflx(i,j,itrc)=init
# ifndef ANA_BTFLUX
            btflxg(i,j,1,itrc)=init
            btflxg(i,j,2,itrc)=init
# endif
          enddo
        enddo
      enddo
      do j=JstrR,JendR
        do i=IstrR,IendR
# ifdef QCORRECTION
          dqdt(i,j)=init
          sst (i,j)=init
#  ifndef ANA_SST
          dqdtg(i,j,1)=init
          sstg (i,j,1)=init
          dqdtg(i,j,2)=init
          sstg (i,j,2)=init
#  endif
#  if defined SALINITY && defined SFLX_CORR 
          sss (i,j)=init
#    ifndef ANA_SSS
          sssg (i,j,1)=init
          sssg (i,j,2)=init
#    endif
#  endif
# endif
!
# ifdef BULK_FLUX
          tair (i,j)=init
          rhum (i,j)=init
          prate (i,j)=init
          radlw (i,j)=init
          radsw (i,j)=init
          wspd(i,j)=init
#  ifdef BULK_SM_UPDATE
          uwnd(i,j)=init
          vwnd(i,j)=init
#  endif
# ifdef DIURNAL_INPUT_SRFLX
          radswbio (i,j)=init
# endif

          tairg(i,j,1)=init
          rhumg (i,j,1)=init
          prateg (i,j,1)=init
          radlwg (i,j,1)=init
          radswg (i,j,1)=init
          tairg(i,j,2)=init
          rhumg (i,j,2)=init
          prateg (i,j,2)=init
          radlwg (i,j,2)=init
          radswg (i,j,2)=init
          wspdg(i,j,1)=init
          wspdg(i,j,2)=init
#  ifdef BULK_SM_UPDATE
          uwndg(i,j,1)=init
          vwndg(i,j,1)=init
          uwndg(i,j,2)=init
          vwndg(i,j,2)=init
#  endif
# ifdef DIURNAL_INPUT_SRFLX
          radswbiog (i,j,1)=init
          radswbiog (i,j,2)=init
# endif
          shflx_rsw(i,j)=init
          shflx_rlw(i,j)=init
          shflx_lat(i,j)=init
          shflx_sen(i,j)=init
#  ifdef AVERAGES
          shflx_rsw_avg(i,j)=init
          shflx_rlw_avg(i,j)=init
          shflx_lat_avg(i,j)=init
          shflx_sen_avg(i,j)=init
#  endif
# endif /* BULK_FLUX */
# ifdef SST_SKIN
          sst_skin=init
#  ifdef AVERAGES
          sst_skin_avg=init
#  endif
# endif

          srflx(i,j)=init
# ifdef AVERAGES
          srflx_avg(i,j)=init
# endif
# ifndef ANA_SRFLUX
          srflxg(i,j,1)=init
          srflxg(i,j,2)=init
# endif
# if (defined BBL || defined MRL_WCI) && defined WAVE_OFFLINE
          wwag(i,j,1)=init
          wwdg(i,j,1)=init
          wwpg(i,j,1)=init
          wwag(i,j,2)=init
          wwdg(i,j,2)=init
          wwpg(i,j,2)=init
#  ifdef WAVE_OFFLINE_BREAKING
          wweb(i,j,1)=init
          wweb(i,j,2)=init
#  endif
# endif
        enddo
      enddo
# if defined PSOURCE_NCFILE
          do i=1,Msrc
           qbar(i)=init
           qbarg(i,1)=init
           qbarg(i,2)=init
          enddo
# endif
#endif /* SOLVE3D */
!
! Initialize climatology arrays (see "climat.h").
!
#if defined ZNUDGING && defined ZCLIMATOLOGY
      do j=JstrR,JendR
        do i=IstrR,IendR
          Znudgcof(i,j)=init
        enddo
      enddo
#endif
#if defined ZCLIMATOLOGY || (defined AGRIF && !defined Z_FRC_BRY)
      do j=JstrR,JendR
        do i=IstrR,IendR
          ssh(i,j)=init
        enddo
      enddo
#endif
#if defined ZCLIMATOLOGY && !defined ANA_SSH
      do j=JstrR,JendR
        do i=IstrR,IendR
          sshg(i,j,1)=init
          sshg(i,j,2)=init
        enddo
      enddo
#endif
#ifdef SOLVE3D
# if defined TNUDGING && defined TCLIMATOLOGY
      do itrc=1,NT
        do k=1,N
          do j=JstrR,JendR
            do i=IstrR,IendR
              Tnudgcof(i,j,k,itrc)=init
            enddo
          enddo
        enddo
      enddo
# endif
# if defined TCLIMATOLOGY || (defined AGRIF && !defined T_FRC_BRY)
      do itrc=1,NT
        do k=1,N
          do j=JstrR,JendR
            do i=IstrR,IendR
              tclm(i,j,k,itrc)=init
            enddo
          enddo
        enddo
      enddo
# endif
# if defined TCLIMATOLOGY && !defined ANA_TCLIMA
      do itrc=1,NT
        do k=1,N
          do j=JstrR,JendR
            do i=IstrR,IendR
              tclima(i,j,k,1,itrc)=init
              tclima(i,j,k,2,itrc)=init
            enddo
          enddo
        enddo
      enddo
# endif
#endif /* SOLVE3D */

#if defined M2CLIMATOLOGY || (defined AGRIF && !defined M2_FRC_BRY)
      do j=JstrR,JendR
        do i=IstrR,IendR
          ubclm(i,j)=init
          vbclm(i,j)=init
        enddo
      enddo
#endif
#if defined M2CLIMATOLOGY && !defined ANA_M2CLIMA
      do j=JstrR,JendR
        do i=IstrR,IendR
          ubclima(i,j,1)=init
          ubclima(i,j,2)=init
          vbclima(i,j,1)=init
          vbclima(i,j,2)=init
        enddo
      enddo
#endif
#if defined SOLVE3D
# if defined M3CLIMATOLOGY && !defined ANA_M3CLIMA
      do k=1,N
        do j=JstrR,JendR
          do i=IstrR,IendR
            uclima(i,j,k,1)=init
            uclima(i,j,k,2)=init
            vclima(i,j,k,1)=init
            vclima(i,j,k,2)=init
          enddo
        enddo
      enddo
# endif
# if defined M3CLIMATOLOGY || \
      (defined AGRIF && !defined M3_FRC_BRY)
      do k=1,N
        do j=JstrR,JendR
          do i=IstrR,IendR
            uclm(i,j,k)=init
            vclm(i,j,k)=init
          enddo
        enddo
      enddo
# endif
#endif
#if defined T_FRC_BRY || defined M2_FRC_BRY || \
      defined M3_FRC_BRY || defined Z_FRC_BRY
!
! Initialize boundary arrays (see "boundary.h").
!
# ifdef OBC_WEST
      if (WESTERN_EDGE) then
#  ifdef Z_FRC_BRY
        do j=JstrR,JendR
          zetabry_west(j)=init
        enddo
#  endif
#  ifdef M2_FRC_BRY
        do j=JstrR,JendR
          ubarbry_west(j)=init
          vbarbry_west(j)=init
        enddo
#  endif
#  if defined SOLVE3D && (defined M3_FRC_BRY || defined T_FRC_BRY)
        do k=1,N
          do j=JstrR,JendR
#   ifdef M3_FRC_BRY
            ubry_west(j,k)=init
            vbry_west(j,k)=init
#   endif
#   ifdef T_FRC_BRY
            do itrc=1,NT
              tbry_west(j,k,itrc)=init
            enddo
#   endif
          enddo
        enddo
#  endif    /* SOLVE3D && (M3_FRC_BRY || T_FRC_BRY)*/
      endif
# endif /* OBC_WEST */
!
# ifdef OBC_EAST
      if (EASTERN_EDGE) then
#  ifdef Z_FRC_BRY
        do j=JstrR,JendR
          zetabry_east(j)=init
        enddo
#  endif
#  ifdef M2_FRC_BRY
        do j=JstrR,JendR
          ubarbry_east(j)=init
          vbarbry_east(j)=init
        enddo
#  endif
#  if defined SOLVE3D && (defined M3_FRC_BRY || defined T_FRC_BRY)
        do k=1,N
          do j=JstrR,JendR
#   ifdef M3_FRC_BRY
            ubry_east(j,k)=init
            vbry_east(j,k)=init
#   endif
#   ifdef T_FRC_BRY
            do itrc=1,NT
              tbry_east(j,k,itrc)=init
            enddo
#   endif
          enddo
        enddo
#  endif    /* SOLVE3D && (M3_FRC_BRY || T_FRC_BRY)*/
      endif
# endif /* OBC_EAST */
!
# ifdef OBC_SOUTH
      if (SOUTHERN_EDGE) then
#  ifdef Z_FRC_BRY
        do i=IstrR,IendR
          zetabry_south(i)=init
        enddo
#  endif
#  ifdef M2_FRC_BRY
        do i=IstrR,IendR
          ubarbry_south(i)=init
          vbarbry_south(i)=init
        enddo
#  endif
#  if defined SOLVE3D && (defined M3_FRC_BRY || defined T_FRC_BRY)
        do k=1,N
          do i=IstrR,IendR
#   ifdef M3_FRC_BRY
            ubry_south(i,k)=init
            vbry_south(i,k)=init
#   endif
#   ifdef T_FRC_BRY
            do itrc=1,NT
              tbry_south(i,k,itrc)=init
            enddo
#   endif
          enddo
        enddo
#  endif    /* SOLVE3D && (M3_FRC_BRY || T_FRC_BRY)*/
      endif
# endif /* OBC_SOUTH */
!
# ifdef OBC_NORTH
      if (NORTHERN_EDGE) then
#  ifdef Z_FRC_BRY
        do i=IstrR,IendR
          zetabry_north(i)=init
        enddo
#  endif
#  ifdef M2_FRC_BRY
        do i=IstrR,IendR
          ubarbry_north(i)=init
          vbarbry_north(i)=init
        enddo
#  endif
#  if defined SOLVE3D && (defined M3_FRC_BRY || defined T_FRC_BRY)
        do k=1,N
          do i=IstrR,IendR
#   ifdef M3_FRC_BRY
            ubry_north(i,k)=init
            vbry_north(i,k)=init
#   endif
#   ifdef T_FRC_BRY
            do itrc=1,NT
              tbry_north(i,k,itrc)=init
            enddo
#   endif
          enddo
        enddo
#  endif    /* SOLVE3D && (M3_FRC_BRY || T_FRC_BRY)*/
      endif
# endif /* OBC_NORTH */
#endif
#if defined SSH_TIDES || defined UV_TIDES
!
! Initialize array for tidal elevation forcing.
!
      do itide=1,Ntides
        Tperiod(itide)=init
        do j=JstrR,JendR
          do i=IstrR,IendR
# ifdef SSH_TIDES
            SSH_Tamp(i,j,itide)=init
            SSH_Tphase(i,j,itide)=init
# endif
# ifdef UV_TIDES
            UV_Tangle(i,j,itide)=init
            UV_Tphase(i,j,itide)=init
            UV_Tmajor(i,j,itide)=init
            UV_Tminor(i,j,itide)=init
# endif
          enddo
        enddo
      enddo
#endif
!
! Set variable horizontal viscosities and tracer diffusion
! coefficients (see "mixing.h") to their background values.
!
#if defined UV_VIS2 || defined SPONGE_VIS2
        do j=JstrR,JendR
          do i=IstrR,IendR
            visc2_r(i,j)=visc2
            visc2_p(i,j)=visc2
            visc2_sponge_r(i,j)=init
            visc2_sponge_p(i,j)=init
          enddo
        enddo
#endif
#ifdef UV_VIS4
        do j=JstrR,JendR
          do i=IstrR,IendR
            visc4_r(i,j)=visc4
            visc4_p(i,j)=visc4
            visc4_sponge_r(i,j)=init
            visc4_sponge_p(i,j)=init
          enddo
        enddo
#endif
#ifdef VIS_COEF_3D
        do k=1,N
          do j=JstrR,JendR
            do i=IstrR,IendR
              visc3d_r(i,j,k)=init
              visc3d_p(i,j,k)=init
            enddo
          enddo
        enddo
#endif
#if defined TS_DIF2 || defined SPONGE_DIF2
        do itrc=1,NT
          do j=JstrR,JendR
            do i=IstrR,IendR
              diff2(i,j,itrc)=tnu2(itrc)
            enddo
          enddo
        enddo
          do j=JstrR,JendR
            do i=IstrR,IendR
              diff2_sponge(i,j)=init
            enddo
          enddo
#endif
#ifdef TS_DIF4
        do itrc=1,NT
          do j=JstrR,JendR
            do i=IstrR,IendR
              diff4(i,j,itrc)=tnu4(itrc)
            enddo
          enddo
        enddo
          do j=JstrR,JendR
            do i=IstrR,IendR
              diff4_sponge(i,j)=init
            enddo
          enddo
#endif
#ifdef DIF_COEF_3D
        do k=1,N
          do j=JstrR,JendR
            do i=IstrR,IendR
              diff3d_u(i,j,k)=init
              diff3d_v(i,j,k)=init
#   ifdef TS_DIF_SMAGO
              diff3d_r(i,j,k)=init
#   endif
            enddo
          enddo
        enddo
#endif
#ifdef TS_MIX_ISO
        do k=1,N-1
          do j=JstrR,JendR
            do i=IstrR,IendR
              dRdx(i,j,k)=0.
              dRde(i,j,k)=0.
            enddo
          enddo
        enddo
#endif
#ifdef SOLVE3D
!
!  Initialize vertical mixing coefficients (see "mixing.h") to their
!  background values. If no vertical closure scheme is selected, the
!  vertical mixing coefficients are those specified by the background
!  values.
!
      do k=0,N
        do j=JstrR,JendR
          do i=IstrR,IendR
# if !defined LMD_MIXING && !defined BVF_MIXING && !defined GLS_MIX2017 && !defined GLS_MIXING
            Akv(i,j,k)=Akv_bak
# else
            Akv(i,j,k)=init
# endif
# if defined BVF_MIXING \
      || defined LMD_MIXING || defined LMD_SKPP || defined LMD_BKPP \
      || defined GLS_MIX2017 || defined GLS_MIXING
            bvf(i,j,k)=init
# endif
          enddo
        enddo
        do j=JstrR,JendR
          do i=IstrR,IendR
# if !defined LMD_MIXING && !defined BVF_MIXING && !defined GLS_MIX2017 && !defined GLS_MIXING
              Akt(i,j,k,itemp)=Akt_bak(itemp)
#  ifdef SALINITY
              Akt(i,j,k,isalt)=Akt_bak(isalt)
#  endif
# else
              Akt(i,j,k,itemp)=init
#  ifdef SALINITY
              Akt(i,j,k,isalt)=init
#  endif
# endif
# ifdef GLS_MIXING
              Akk(i,j,k)=Akk_bak
              Akp(i,j,k)=Akp_bak
# endif
          enddo
        enddo
      enddo
!
# if (defined LMD_SKPP || defined LMD_BKPP) && defined LMD_NONLOCAL
      do k=1,N
        do j=JstrR,JendR
          do i=IstrR,IendR
            ghats(i,j,k)=init
          enddo
        enddo
      enddo
# endif
# ifdef LMD_SKPP
!
!  Initialize depth of surface boundary layer.
!
      do j=JstrR,JendR
        do i=IstrR,IendR
#ifdef LMD_SKPP2005
          hbls(i,j,1)=init          
          hbls(i,j,2)=init 
#else        
          hbl(i,j)=init
#endif         
          kbl(i,j)=init
#  ifdef AVERAGES
          hbl_avg(i,j)=init
#  endif
#  if defined DIAGNOSTICS_TS && defined DIAGNOSTICS_TS_MLD
          kbl_nstp(i,j)=init
#  endif
        enddo
      enddo
# endif /* LMD_SKPP */
# ifdef LMD_BKPP
!
!  Initialize depth of bottom boundary layer.
!
      do j=JstrR,JendR
        do i=IstrR,IendR
          hbbl(i,j)=init
#ifdef AVERAGES
          hbbl_avg(i,j)=init
#endif
        enddo
      enddo
# endif /* LMD_BKPP */
# ifdef LMD_SKPP
#  ifdef BIOLOGY
!
!  Initialize depth of the euphotic layer.
!
      do j=JstrR,JendR
        do i=IstrR,IendR
          hel(i,j)=init
#  if (defined BIO_NChlPZD && defined OXYGEN)
          u10(i,j)=init
          Kv_O2(i,j)=init
          O2satu(i,j)=init
#  endif
# ifdef AVERAGES
          hel_avg(i,j)=init
#  if (defined BIO_NChlPZD && defined OXYGEN)
          u10_avg(i,j)=init
          Kv_O2_avg(i,j)=init
          O2satu_avg(i,j)=init
#  endif
# endif
        enddo
      enddo
#  endif /* BIOLOGY */
# endif /* LMD_SKPP */

# ifdef GLS_MIXING
      do k=0,N
        do j=JstrR,JendR
          do i=IstrR,IendR
            tke(i,j,k,1)=gls_Kmin
            tke(i,j,k,2)=gls_Kmin
            gls(i,j,k,1)=gls_Pmin
            gls(i,j,k,2)=gls_Pmin
            Lscale(i,j,k)=init
          enddo
        enddo
      enddo
# endif

# ifdef NBQ
      do j=JstrR,JendR
        do i=IstrR,IendR
#  ifdef NBQ_MASS
          rhobar_nbq(i,j,:)   =1.
          rhobar_nbq_avg1(i,j)=1.
#  endif
          rubar_nbq(i,j)    =init
          rvbar_nbq(i,j)    =init
          rufrc_bak(i,j,1:2)=init
          rvfrc_bak(i,j,1:2)=init
          
        enddo
      enddo
      do k=1,N
        do j=JstrR,JendR
          do i=IstrR,IendR
# ifdef NBQ_MASS
            rho_nbq(i,j,k)     =1.
            rho_nbq_avg1(i,j,k)=1.
!            rho_nbq_avg2(i,j,k)=1.
# endif
            ru_nbq(i,j,k)     =init
            rv_nbq(i,j,k)     =init
#if !defined M2FILTER_NONE
            ru_nbq_avg1(i,j,k)=init
            rv_nbq_avg1(i,j,k)=init
#endif
            ru_nbq_avg2(i,j,k)=init
            rv_nbq_avg2(i,j,k)=init
            ru_int_nbq(i,j,k) =init
            rv_int_nbq(i,j,k) =init
          enddo
        enddo
      enddo
      do k=0,N
        do j=JstrR,JendR
          do i=IstrR,IendR
            wz(i,j,k,1)  =init
            wz(i,j,k,2)  =init
            rw_nbq(i,j,k)=init
#if !defined M2FILTER_NONE
            rw_nbq_avg1(i,j,k)=init
#endif
            rw_nbq_avg2(i,j,k)=init
            rw_int_nbq(i,j,k) =init
          enddo
        enddo
      enddo
# endif

#endif /* SOLVE3D */

      return
      end
