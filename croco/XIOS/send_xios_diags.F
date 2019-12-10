!     $Id:$
!     
!======================================================================
!     ROMS_AGRIF is a branch of ROMS developped at IRD and INRIA, in France
!     The two other branches from UCLA (Shchepetkin et al) 
!     and Rutgers University (Arango et al) are under MIT/X style license.
!     ROMS_AGRIF specific routines (nesting) are under CeCILL-C license.
!     
!     ROMS_AGRIF website : http://roms.mpl.ird.fr
!======================================================================
!     
#include "cppdefs.h"

#ifdef XIOS

      MODULE xios_module
      USE xios
      IMPLICIT NONE

      PUBLIC set_grid
      PUBLIC iom_set_domain_attr
      PUBLIC iom_set_axis_attr
      PUBLIC iom_set_field_attr
      PUBLIC iom_set_file_attr
      PUBLIC iom_set_file_global_attr_char
      PUBLIC iom_set_file_global_attr_float
      PUBLIC iom_get_file_attr
      PUBLIC set_scalar
      




      CONTAINS
!#include "cppdefs.h"


      SUBROUTINE set_grid( cdgrd, plon, plat)
# include "param.h"
# include "ncscrum.h"
# include "scalars.h"
      !!----------------------------------------------------------------------
      !!                     ***  ROUTINE set_grid  ***
      !!
      !! ** Purpose :   set domain definition to XIOS library
      !!
      !! ** Method  :   
      !!
      !! ** History :   S. Theetten , adaptation for XIOS-2.5 (May, 2019)
      !!                A. Ponte    , adaptation for CROCO
      !!                S. Masson   , original version (NEMO)  
      !!
      !!----------------------------------------------------------------------
      CHARACTER(LEN=*)            , INTENT(in) ::   cdgrd
      REAL plon(GLOBAL_2D_ARRAY)
      REAL plat(GLOBAL_2D_ARRAY)
      !
      INTEGER vtype, horiz_type, vert_type
      INTEGER imin, imax, jmin, jmax, ni, nj
      INTEGER nx_glo, ny_glo 
      INTEGER istart, jstart
      INTEGER data_ni,data_nj 
      INTEGER :: ji
      !!----------------------------------------------------------------------
      !
      ! See ncscrum.h for r2dvar ...
      SELECT CASE (TRIM(cdgrd))
      CASE ("rho")
         vtype=r2dvar
         nx_glo=LLm+2
         ny_glo=MMm+2
      CASE ("b")
         vtype=r2dvar
         nx_glo=LLm+2
         ny_glo=MMm+2
      CASE ("u")
         vtype=u2dvar
         nx_glo=LLm+1
         ny_glo=MMm+2
      CASE ("v")
         vtype=v2dvar
         nx_glo=LLm+2
         ny_glo=MMm+1
      CASE ("w")
         vtype=r2dvar
         nx_glo=LLm+2
         ny_glo=MMm+2
      CASE ("psi")
         vtype=p2dvar
         nx_glo=LLm+1
         ny_glo=MMm+1
      END SELECT

      ! From nf_fread.F
      vert_type=vtype/4         ! vertical and horizontal
      horiz_type=vtype-4*vert_type ! grid types, then calculate
      jmin=horiz_type/2         ! starting indices indices
      imin=horiz_type-2*jmin    ! in horizontal directions.
      !     
      istart=1
      jstart=1

# ifdef MPI
      IF (ii.gt.0) then
         istart=1-imin+iminmpi
         imin=1
      ENDIF
      IF (ii.eq.NP_XI-1) then
         imax=Lmmpi+1
      ELSE
         imax=Lmmpi
      ENDIF
      IF (jj.gt.0) then
         jstart=1-jmin+jminmpi
         jmin=1
      ENDIF
      IF (jj.eq.NP_ETA-1) then
         jmax=Mmmpi+1
      ELSE
         jmax=Mmmpi
      ENDIF
# else
      imin=1
      imax=Lm+1
      jmin=1
      jmax=Mm+1
# endif

      ! there is probably already a variable for ni and nj 
      ni=imax-imin+1
      nj=jmax-jmin+1

      ! from set_global_definitions.h
# ifdef THREE_GHOST_POINTS
#  ifdef MPI
      data_ni = Lm+3+padd_X+2+1
      data_nj = Mm+3+padd_E+2+1
#  endif
# else
#  ifdef MPI
      data_ni = Lm+2+padd_X+1+1
      data_nj = Mm+2+padd_E+1+1
#  endif
# endif

      ! notes iminmpi, jminmpi are computed in MPI_Setup.F
# if ! defined XIOS2
      CALL iom_set_domain_attr(TRIM(cdgrd),ni_glo=nx_glo, nj_glo=ny_glo,
     &     ibegin=istart, ni=ni,jbegin=jstart,nj=nj)
# else
      CALL iom_set_domain_attr(TRIM(cdgrd),ni_glo=nx_glo, nj_glo=ny_glo,
     &     ibegin=istart-1, ni=ni,jbegin=jstart-1,nj=nj)
# endif
      ! to do: below should have be rewritten to match all configs (periodic, number of ghost points, ...)
      ! see set_global_definitions.h
      CALL iom_set_domain_attr(TRIM(cdgrd),data_dim=2,
     &     data_ibegin=START_1D_ARRAYXI-imin,
     &     data_ni=data_ni,
     &     data_jbegin=START_1D_ARRAYETA-jmin, 
     &     data_nj=data_nj)

      CALL iom_set_domain_attr(TRIM(cdgrd),
     &     lonvalue=RESHAPE(plon(imin:imax,jmin:jmax), (/ni*nj/)),
     &     latvalue=RESHAPE(plat(imin:imax,jmin:jmax), (/ni*nj/)))
      !      
      END SUBROUTINE set_grid

#ifdef MASKING
      SUBROUTINE set_mask(cdgrd)
# include "param.h"
# include "ncscrum.h"
# include "scalars.h"
# include "grid.h"
      !!----------------------------------------------------------------------
      !!                     ***  ROUTINE set_mask  ***
      !!
      !! ** Purpose :   set mask in domain definition to XIOS library
      !!
      !! ** Method  :   
      !!
      !! ** History :   S. Petton   , adapt from set_grid
      !!                S. Theetten , adaptation for XIOS-2.5 (May, 2019)
      !!                A. Ponte    , adaptation for CROCO
      !!                S. Masson   , original version (NEMO)  
      !!
      !!----------------------------------------------------------------------
      CHARACTER(LEN=*)            , INTENT(in) ::   cdgrd
      LOGICAL zmask(GLOBAL_2D_ARRAY)
      !
      INTEGER vtype, horiz_type, vert_type
      INTEGER imin, imax, jmin, jmax, ni, nj
      !!----------------------------------------------------------------------
      !
      zmask(:,:) = .false.
      SELECT CASE (TRIM(cdgrd))
      CASE ("rho")
         vtype=r2dvar
         where(rmask(:,:) > 0) zmask(:,:) = .true.
      CASE ("b")
         vtype=r2dvar
         where(rmask(:,:) > 0) zmask(:,:) = .true.
      CASE ("u")
         vtype=u2dvar
         where(umask(:,:) > 0) zmask(:,:) = .true.
      CASE ("v")
         vtype=v2dvar
         where(vmask(:,:) > 0) zmask(:,:) = .true.
      CASE ("w")
         vtype=r2dvar
         where(rmask(:,:) > 0) zmask(:,:) = .true.
      CASE ("psi")
         vtype=p2dvar
         where(pmask(:,:) > 0) zmask(:,:) = .true.
      END SELECT

      ! From nf_fread.F
      vert_type=vtype/4            ! vertical and horizontal
      horiz_type=vtype-4*vert_type ! grid types, then calculate
      jmin=horiz_type/2            ! starting indices indices
      imin=horiz_type-2*jmin       ! in horizontal directions.

# ifdef MPI
      IF (ii.gt.0) then
         imin=1
      ENDIF
      IF (ii.eq.NP_XI-1) then
         imax=Lmmpi+1
      ELSE
         imax=Lmmpi
      ENDIF
      IF (jj.gt.0) then
         jmin=1
      ENDIF
      IF (jj.eq.NP_ETA-1) then
         jmax=Mmmpi+1
      ELSE
         jmax=Mmmpi
      ENDIF
# else
      imin=1
      imax=Lm+1
      jmin=1
      jmax=Mm+1
# endif

      ! there is probably already a variable for ni and nj 
      ni=imax-imin+1
      nj=jmax-jmin+1

      CALL iom_set_domain_attr(TRIM(cdgrd),
     &     mask=RESHAPE(zmask(imin:imax,jmin:jmax), (/ni*nj/)))
      
      END SUBROUTINE set_mask
#endif

      SUBROUTINE iom_set_domain_attr( cdid, ni_glo, nj_glo, ibegin, jbegin,
     & ni, nj, zoom_ibegin, zoom_jbegin, zoom_ni, zoom_nj,
     & data_dim, data_ibegin, data_ni, data_jbegin, data_nj, 
     & lonvalue, latvalue, mask,nvertex, bounds_lon, bounds_lat, area )
      !!----------------------------------------------------------------------
      !!                     ***  ROUTINE iom_set_domain_attr  ***
      !!
      !! ** Purpose :   set domain definition to XIOS library
      !!
      !! ** Method  :   
      !!
      !! ** History :   S. Theetten , adaptation for XIOS-2.5 (May, 2019)
      !!                R. Benshila , adaptation for CROCO
      !!                S. Masson   , original version (NEMO)  
      !!
      !!----------------------------------------------------------------------
      CHARACTER(LEN=*) , INTENT(in) :: cdid
      INTEGER , OPTIONAL, INTENT(in) :: ni_glo, nj_glo, ibegin, jbegin 
      INTEGER , OPTIONAL, INTENT(in) :: ni, nj
      INTEGER , OPTIONAL, INTENT(in) :: data_dim, data_ibegin, data_ni
      INTEGER , OPTIONAL, INTENT(in) :: data_jbegin, data_nj
      INTEGER , OPTIONAL, INTENT(in) :: zoom_ibegin, zoom_jbegin
      INTEGER , OPTIONAL, INTENT(in) :: zoom_ni, zoom_nj, nvertex
      REAL, DIMENSION(:)   , OPTIONAL, INTENT(in) ::   lonvalue, latvalue
      REAL, DIMENSION(:,:) , OPTIONAL, INTENT(in) ::   bounds_lon, 
     & bounds_lat, area
# if ! defined XIOS2
      LOGICAL,  DIMENSION(:,:) , OPTIONAL, INTENT(in) ::   mask
# else
      LOGICAL,  DIMENSION(:) , OPTIONAL, INTENT(in) ::   mask
# endif
      !!----------------------------------------------------------------------
      !
# if ! defined XIOS2
      IF ( xios_is_valid_domain     (cdid) ) THEN
         CALL xios_set_domain_attr     ( cdid, ni_glo=ni_glo, nj_glo=nj_glo, 
     &    ibegin=ibegin, jbegin=jbegin, ni=ni, nj=nj,  
     &    data_dim=data_dim, data_ibegin=data_ibegin, data_ni=data_ni, 
     &    data_jbegin=data_jbegin, data_nj=data_nj ,  
     &    zoom_ibegin=zoom_ibegin, zoom_jbegin=zoom_jbegin, 
     &    zoom_ni=zoom_ni, zoom_nj=zoom_nj,                       
     &    lonvalue=lonvalue, latvalue=latvalue, mask=mask,
     &    nvertex=nvertex, bounds_lon=bounds_lon,                 
     &    bounds_lat=bounds_lat)!RB, area=area )
      ENDIF
      IF ( xios_is_valid_domaingroup(cdid) ) THEN
         CALL xios_set_domaingroup_attr( cdid, ni_glo=ni_glo, nj_glo=nj_glo,
     &    ibegin=ibegin, jbegin=jbegin, ni=ni, nj=nj,
     &    data_dim=data_dim, data_ibegin=data_ibegin, data_ni=data_ni, 
     &    data_jbegin=data_jbegin, data_nj=data_nj ,   
     &    zoom_ibegin=zoom_ibegin, zoom_jbegin=zoom_jbegin, 
     &    zoom_ni=zoom_ni, zoom_nj=zoom_nj,                       
     &    lonvalue=lonvalue, latvalue=latvalue, mask=mask, 
     &    nvertex=nvertex, bounds_lon=bounds_lon,                  
     &    bounds_lat=bounds_lat)!RB, area=area )
      ENDIF

# else

      IF ( xios_is_valid_domain(cdid) ) THEN
         CALL xios_set_domain_attr( cdid, ni_glo=ni_glo, nj_glo=nj_glo, 
     &    ibegin=ibegin, jbegin=jbegin, ni=ni, nj=nj,
     &    data_dim=data_dim, data_ibegin=data_ibegin, data_ni=data_ni,
     &    data_jbegin=data_jbegin, data_nj=data_nj,   
     &    lonvalue_1D=lonvalue, latvalue_1D=latvalue, mask_1D=mask,
     &    nvertex=nvertex, bounds_lon_1D=bounds_lon,       
     &    bounds_lat_1D=bounds_lat, area=area, type='curvilinear')
      ENDIF
      !
      IF ( xios_is_valid_domaingroup(cdid) ) THEN
         CALL xios_set_domaingroup_attr( cdid, ni_glo=ni_glo, nj_glo=nj_glo, 
     &    ibegin=ibegin, jbegin=jbegin, ni=ni, nj=nj,
     &    data_dim=data_dim, data_ibegin=data_ibegin, data_ni=data_ni, 
     &    data_jbegin=data_jbegin, data_nj=data_nj ,   
     &    lonvalue_1D=lonvalue, latvalue_1D=latvalue, mask_1D=mask,
     & nvertex=nvertex, bounds_lon_1D=bounds_lon,                  
     &    bounds_lat_1D=bounds_lat, area=area, type='curvilinear')
      ENDIF
# endif
      !
      CALL xios_solve_inheritance()
      ! 
      END SUBROUTINE iom_set_domain_attr

# if defined XIOS2
      SUBROUTINE iom_set_zoom_domain_attr( cdid, ibegin, jbegin, ni, nj)
      CHARACTER(LEN=*)          , INTENT(in) ::   cdid
      INTEGER         , OPTIONAL, INTENT(in) ::   ibegin, jbegin, ni, nj

      IF ( xios_is_valid_zoom_domain     (cdid) ) THEN
         CALL xios_set_zoom_domain_attr     ( cdid, ibegin=ibegin-1, 
     &        jbegin=jbegin-1, ni=ni, nj=nj)
      ENDIF
      END SUBROUTINE iom_set_zoom_domain_attr
# endif

      
      SUBROUTINE iom_set_axis_attr( cdid, paxis,bounds )
      !!----------------------------------------------------------------------
      !!                     ***  ROUTINE iom_set_axis_attr  ***
      !!
      !! ** Purpose :   set axis attribute for XIOS
      !!
      !! ** Method  :   
      !!
      !! ** History :   S. Theetten , adaptation for XIOS-2.5 (May, 2019)
      !!                A. Ponte    , adaptation for CROCO
      !!                S. Masson   , original version (NEMO)
      !!
      !!----------------------------------------------------------------------
      CHARACTER(LEN=*)      , INTENT(in) ::   cdid
      REAL, DIMENSION(:)  , OPTIONAL, INTENT(in) ::   paxis
      REAL, DIMENSION(:,:), OPTIONAL, INTENT(in) ::   bounds
      !!----------------------------------------------------------------------
      !
      IF ( PRESENT(paxis) ) THEN
# if ! defined XIOS2
      IF ( xios_is_valid_axis (cdid) ) then
       CALL xios_set_axis_attr(cdid, size=size(paxis),value=paxis )
      ENDIF
      !
      IF ( xios_is_valid_axisgroup(cdid) ) then
       CALL xios_set_axisgroup_attr(cdid, size=size(paxis),value=paxis )
      ENDIF
# else
      IF ( xios_is_valid_axis (cdid) ) then
       CALL xios_set_axis_attr(cdid,n_glo=SIZE(paxis), value=paxis )
      ENDIF
      !
      IF ( xios_is_valid_axisgroup(cdid) ) then
       CALL xios_set_axisgroup_attr(cdid, n_glo=SIZE(paxis),value=paxis)
      ENDIF
# endif
      !
      ENDIF
      IF ( xios_is_valid_axis     (cdid) )  then 
         CALL xios_set_axis_attr     ( cdid, bounds=bounds )
      ENDIF
      IF ( xios_is_valid_axisgroup(cdid) )  then
         CALL xios_set_axisgroup_attr( cdid, bounds=bounds )
      ENDIF
      
      CALL xios_solve_inheritance()
      !
      END SUBROUTINE iom_set_axis_attr
 
      SUBROUTINE iom_set_field_attr( cdid, freq_op, freq_offset )
      !!----------------------------------------------------------------------
      !!                     ***  ROUTINE iom_set_field_attr  ***
      !!
      !! ** Purpose :   set output field attribute for XIOS
      !!
      !! ** Method  :   
      !!
      !! ** History :   S. Theetten , adaptation for XIOS-2.5 (May, 2019)
      !!                A. Ponte    , adaptation for CROCO 
      !!                S. Masson   , original version (NEMO)  
      !!
      !!----------------------------------------------------------------------
      CHARACTER(LEN=*)          , INTENT(in) ::   cdid
      TYPE(xios_file) :: file_hdl
      TYPE(xios_variable) :: var_hdl
      LOGICAL ok
# if ! defined XIOS2
      CHARACTER(LEN=*),OPTIONAL , INTENT(in) ::   freq_op
      CHARACTER(LEN=*),OPTIONAL , INTENT(in) ::   freq_offset
# else
      TYPE(xios_duration),OPTIONAL , INTENT(in) ::   freq_op
      TYPE(xios_duration),OPTIONAL , INTENT(in) ::   freq_offset
# endif      
      !!----------------------------------------------------------------------
      !
      IF ( xios_is_valid_field     (cdid) ) 
     &   CALL xios_set_field_attr   
     &  ( cdid, freq_op=freq_op, freq_offset=freq_offset )
      !
      IF ( xios_is_valid_fieldgroup(cdid) )   
     & CALL xios_set_fieldgroup_attr( 
     &     cdid, freq_op=freq_op, freq_offset=freq_offset )
      !
      CALL xios_solve_inheritance()
      !
      END SUBROUTINE iom_set_field_attr

      SUBROUTINE iom_set_file_global_attr_char(cdid, type, name, value)
# include "param.h"
# include "scalars.h"
      !!----------------------------------------------------------------------
      !!                     ***  ROUTINE iom_set_file_global_attr_char  ***
      !!
      !! ** Purpose :   set file global attribute for XIOS
      !!
      !! ** Method  :   
      !!
      !! ** History :   S. Theetten : add global attributes , May, 2019

      !!
      !!----------------------------------------------------------------------
      CHARACTER(LEN=*)          , INTENT(in) ::   cdid
      CHARACTER(LEN=*)          , INTENT(in) ::   type
      CHARACTER(LEN=*),OPTIONAL , INTENT(in) ::   name, value
      LOGICAL :: ok
      TYPE(xios_file) :: file_hdl
      TYPE(xios_variable) :: var_hdl
      !!----------------------------------------------------------------------

      CALL xios_solve_inheritance()

      IF ( xios_is_valid_file     (cdid) )   THEN
          CALL xios_get_handle(cdid,file_hdl)

          CALL xios_add_child(file_hdl, var_hdl, trim(name))
          CALL xios_set_variable_attr(trim(name), name=name, type=type)
          ok = xios_setvar(trim(name),value)
      ENDIF

      CALL xios_solve_inheritance()
      !
      END SUBROUTINE iom_set_file_global_attr_char


      SUBROUTINE iom_set_file_global_attr_float(cdid, type, name, value)
# include "param.h"
# include "scalars.h"
      !!----------------------------------------------------------------------
      !!                     ***  ROUTINE iom_set_file_global_attr_char  ***
      !!
      !! ** Purpose :   set file global attribute for XIOS
      !!
      !! ** Method  :   
      !!
      !! ** History :   S. Theetten : add global attributes , May, 2019

      !!
      !!----------------------------------------------------------------------
      CHARACTER(LEN=*)          , INTENT(in) ::   cdid
      CHARACTER(LEN=*)          , INTENT(in) ::   type
      CHARACTER(LEN=*),OPTIONAL , INTENT(in) ::   name
	  REAL                                   ::   value
      LOGICAL :: ok
      TYPE(xios_file) :: file_hdl
      TYPE(xios_variable) :: var_hdl
      !!----------------------------------------------------------------------

      CALL xios_solve_inheritance()

      IF ( xios_is_valid_file     (cdid) )   THEN
          CALL xios_get_handle(cdid,file_hdl)

          CALL xios_add_child(file_hdl, var_hdl, trim(name))
          CALL xios_set_variable_attr(trim(name), name=name, type=type)
          ok = xios_setvar(trim(name),value)
      ENDIF

      CALL xios_solve_inheritance()
      !
      END SUBROUTINE iom_set_file_global_attr_float
	  
      SUBROUTINE iom_set_file_global_attr_int(cdid, type, name, value)
# include "param.h"
# include "scalars.h"
      !!----------------------------------------------------------------------
      !!                     ***  ROUTINE iom_set_file_global_attr_int  ***
      !!
      !! ** Purpose :   set file global attribute for XIOS
      !!
      !! ** Method  :   
      !!
      !! ** History :   S. Theetten : add global attributes , May, 2019

      !!
      !!----------------------------------------------------------------------
      CHARACTER(LEN=*)          , INTENT(in) ::   cdid
      CHARACTER(LEN=*)          , INTENT(in) ::   type
      CHARACTER(LEN=*),OPTIONAL , INTENT(in) ::   name
      INTEGER                                ::   value
      LOGICAL :: ok
      TYPE(xios_file) :: file_hdl
      TYPE(xios_variable) :: var_hdl
      !!----------------------------------------------------------------------

      CALL xios_solve_inheritance()

      IF ( xios_is_valid_file     (cdid) )   THEN
          CALL xios_get_handle(cdid,file_hdl)

          CALL xios_add_child(file_hdl, var_hdl, trim(name))
          CALL xios_set_variable_attr(trim(name), name=name, type=type)
          ok = xios_setvar(trim(name),value)
      ENDIF

      CALL xios_solve_inheritance()
      !
      END SUBROUTINE iom_set_file_global_attr_int
	  
	  
      SUBROUTINE iom_set_file_attr( cdid, name, name_suffix )
# include "param.h"
# include "scalars.h"
      !!----------------------------------------------------------------------
      !!                     ***  ROUTINE iom_set_file_attr  ***
      !!
      !! ** Purpose :   set file attribute for XIOS
      !!
      !! ** Method  :   
      !!
      !! ** History :   A. Ponte   , adaptation for CROCO 
      !!                S. Masson  , original version (NEMO)  
      !!
      !!----------------------------------------------------------------------
      CHARACTER(LEN=*)          , INTENT(in) ::   cdid
      CHARACTER(LEN=*),OPTIONAL , INTENT(in) ::   name, name_suffix
      !!----------------------------------------------------------------------
      IF ( xios_is_valid_file     (cdid) )   CALL xios_set_file_attr  
     &    ( cdid, name=name, name_suffix=name_suffix )
      !
      IF ( xios_is_valid_filegroup(cdid) )   CALL xios_set_filegroup_attr
     & ( cdid, name=name, name_suffix=name_suffix )
      !
      CALL xios_solve_inheritance()
      !
      END SUBROUTINE iom_set_file_attr

      
      SUBROUTINE iom_get_file_attr( cdid, name, name_suffix, output_freq )
# include "param.h"
# include "scalars.h"
# include "ncscrum.h"
      !!----------------------------------------------------------------------
      !!                     ***  ROUTINE iom_get_file_attr  ***
      !!
      !! ** Purpose :   get attribute from xml file
      !!
      !! ** Method  :   
      !!
      !! ** History :   R. Benshila, adaptation for CROCO 
      !!                S. Masson  , original version (NEMO)  
      !!
      !!----------------------------------------------------------------------
      CHARACTER(LEN=*)          , INTENT(in ) :: 
     &   cdid
      CHARACTER(LEN=*),OPTIONAL , INTENT(out) :: 
     &   name, name_suffix
# if ! defined XIOS2
      CHARACTER(LEN=*),OPTIONAL , INTENT(out) ::    output_freq
# else
      TYPE(xios_duration)   ,OPTIONAL , INTENT(out) :: output_freq
# endif 
      LOGICAL                                 :: 
     &   llexist1,llexist2,llexist3
      !---------------------------------------------------------------------
      !
      IF( PRESENT( name        ) )   name = ''          ! default values
      IF( PRESENT( name_suffix ) )   name_suffix = ''
# if ! defined XIOS2
      IF( PRESENT( output_freq ) )   output_freq = ''
# else
      IF( PRESENT( output_freq ) )   output_freq = xios_duration(0,0,0,0,0,0)
# endif
       IF ( xios_is_valid_file     (cdid) ) THEN
         CALL xios_solve_inheritance()
         CALL xios_is_defined_file_attr   
     &      ( cdid, name = llexist1, name_suffix = llexist2,
     &        output_freq = llexist3)
         IF(llexist1)   CALL xios_get_file_attr   
     &      ( cdid, name = name )
         IF(llexist2)   CALL xios_get_file_attr   
     &      ( cdid, name_suffix = name_suffix )
         IF(llexist3)   CALL xios_get_file_attr   
     &      ( cdid, output_freq = output_freq )
      ENDIF

      IF ( xios_is_valid_filegroup(cdid) ) THEN
         CALL xios_solve_inheritance()
         CALL xios_is_defined_filegroup_attr( cdid, name = llexist1, 
     &    name_suffix = llexist2, output_freq = llexist3)
         IF(llexist1)   CALL xios_get_filegroup_attr( cdid, 
     &    name = name )
         IF(llexist2)   CALL xios_get_filegroup_attr( cdid, 
     &    name_suffix = name_suffix )
         IF(llexist3)   CALL xios_get_filegroup_attr( cdid, 
     &    output_freq = output_freq )
      ENDIF
      !
      END SUBROUTINE iom_get_file_attr

      
      SUBROUTINE set_scalar
# include "param.h"
# include "scalars.h"
      !!----------------------------------------------------------------------
      !!                     ***  ROUTINE set_scalar  ***
      !!
      !! ** Purpose :   fake axis for scalar output
      !!
      !! ** Method  :   
      !! ** History :   S. Theetten , adaptation for XIOS-2.5 (May, 2019)
      !!                R. Benshila , adaptation for CROCO 
      !!                S. Masson   , original version (NEMO)  
      !!
      !!----------------------------------------------------------------------
      REAL, DIMENSION(1) ::   zz = 1.
      !!----------------------------------------------------------------------
      ! 
# if ! defined XIOS2
      CALL iom_set_domain_attr('scalarpoint', ni_glo=NNODES, nj_glo=1, 
     &  ibegin=mynode+1, jbegin=1, ni=1, nj=1)
!      CALL iom_set_domain_attr('scalarpoint', data_dim=2,data_ibegin = 1, 
!     &                         data_ni = 1, data_jbegin = 1, data_nj = 1)
# else
      CALL iom_set_domain_attr('scalarpoint', ni_glo=NNODES, nj_glo=1,
     &  ibegin=mynode, jbegin=0, ni=1, nj=1)
!      CALL iom_set_domain_attr('scalarpoint', data_dim=2,data_ibegin = 0,
!     &                         data_ni = 1, data_jbegin = 0, data_nj = 1)
# endif
      CALL iom_set_domain_attr('scalarpoint', data_dim=2,data_ibegin = 1, 
     &                         data_ni = 1, data_jbegin = 1, data_nj = 1)
      zz=REAL(mynode+1)
      CALL iom_set_domain_attr('scalarpoint', lonvalue=zz, 
     &  latvalue=zz)
      !
      END SUBROUTINE set_scalar

      SUBROUTINE set_xml_global_att(cdid)
# include "param.h"
# include "scalars.h"
# include "ncscrum.h"
# include "strings.h"
      !!----------------------------------------------------------------------
      !!                     ***  ROUTINE set_xml_global_att  ***
      !!
      !! ** Purpose :   automatic definitions of global attributs
      !!
      !! ** Method  :
      !!
      !! ** History :   S. Theetten, (May, 2019)

      !!
      !!----------------------------------------------------------------------
      CHARACTER(LEN=*)          , INTENT(in) ::   cdid
      CHARACTER(LEN=40)  ::   char_value

      char_value=trim(title)
      CALL iom_set_file_global_attr_char   (cdid
     &     ,type='string',name='title', value=char_value)
#ifdef START_DATE
C start_date
#endif
C       char_value=trim(date_str)
C        CALL iom_set_file_global_attr   (cdid
C       &          ,type='string',name='date', value=char_value)
      char_value=trim(rstname)
      CALL iom_set_file_global_attr_char   (cdid
     &     ,type='string',name='rst_file', value=char_value)
#ifndef ANA_GRID
      char_value=trim(grdname)
      CALL iom_set_file_global_attr_char   (cdid
     &     ,type='string',name='grd_file', value=char_value)
#endif
#ifndef ANA_INITIAL
      char_value=trim(ininame)
      CALL iom_set_file_global_attr_char   (cdid
     &     ,type='string',name='ini_file', value=char_value)
#endif
#if !defined ANA_SMFLUX || !defined ANA_STFLUX \
      || !defined ANA_BTFLUX \
      || (defined BBL && !defined ANA_BSEDIM)  \
      || (defined BBL && !defined ANA_WWAVE)   \
      || (defined SALINITY && !defined ANA_SSFLUX)  \
      || ((defined LMD_SKPP || defined LMD_BKPP)  && !defined ANA_SRFLUX)
      char_value=trim(frcname)
      CALL iom_set_file_global_attr_char   (cdid
     &     ,type='string',name='frc_file', value=char_value)
#endif
#ifdef PSOURCE_NCFILE
C qbar_file
#endif
#ifdef ASSIMILATION
C ass_file
C apar_file
#endif
#ifdef SOLVE3D
# ifdef  NEW_S_COORD
      CALL iom_set_file_global_attr_char   (cdid
     &     ,type='string',name='VertCoordType', value='NEW')
# endif
# ifdef  LMD_SKPP2005
      CALL iom_set_file_global_attr_char   (cdid
     &     ,type='string',name='skpp', value='2005')
# endif
      CALL iom_set_file_global_attr_float   (cdid
     &     ,type='float',name='theta_s', value=theta_s)
      CALL iom_set_file_global_attr_char   (cdid
     &     ,type='string',name='theta_s_expl'
     &     ,value='S-coordinate surface control parameter')
      CALL iom_set_file_global_attr_float   (cdid
     &     ,type='float',name='theta_b', value=theta_b)
      CALL iom_set_file_global_attr_char   (cdid
     &     ,type='string',name='theta_b_expl'
     &     , value='S-coordinate bottom control parameter')
      CALL iom_set_file_global_attr_float   (cdid
     &     ,type='float',name='Tcline', value=Tcline)
      CALL iom_set_file_global_attr_char   (cdid
     &     ,type='string',name='Tcline_expl'
     &     , value='S-coordinate surface/bottom layer width')
      CALL iom_set_file_global_attr_char   (cdid
     &     ,type='string',name='Tcline_units', value='meter')
      CALL iom_set_file_global_attr_float   (cdid
     &      ,type='float',name='hc', value=hc)
      CALL iom_set_file_global_attr_char   (cdid
     &     ,type='string',name='hc_expl'
     &     , value='S-coordinate parameter, critical depth')
      CALL iom_set_file_global_attr_char   (cdid
     &     ,type='string',name='hc_units', value='meter')
      CALL iom_set_file_global_attr_char   (cdid
     &     ,type='string',name='sc_w'
     &     , value='have a look at variable sc_w in this file')
      CALL iom_set_file_global_attr_char   (cdid
     &     ,type='string',name='Cs_w'
     &     , value='have a look at variable Cs_w in this file')
      CALL iom_set_file_global_attr_char   (cdid
     &     ,type='string',name='sc_r'
     &     , value='have a look at variable sc_r in this file')
      CALL iom_set_file_global_attr_char   (cdid
     &     ,type='string',name='Cs_r'
     &     , value='have a look at variable Cs_r in this file')
#endif
      CALL iom_set_file_global_attr_int   (cdid
     &     ,type='int',name='ntimes', value=ntimes)
      CALL iom_set_file_global_attr_int   (cdid
     &     ,type='int',name='ndtfast', value=ndtfast)
      CALL iom_set_file_global_attr_float   (cdid
     &     ,type='float',name='dt', value=dt)
      CALL iom_set_file_global_attr_float   (cdid
     &     ,type='float',name='dtfast', value=dtfast)
      CALL iom_set_file_global_attr_int   (cdid
     &      ,type='int',name='nwrt', value=nwrt)
#ifdef UV_VIS2
#endif
#ifdef UV_VIS4
#endif
#ifdef SOLVE3D
# ifdef TS_DIF2
# endif
# ifdef TS_DIF4
C       CALL iom_set_file_global_attr_float   (cdid
C      &     ,type='float',name='tnu4', value=tnu4)
      CALL iom_set_file_global_attr_char   (cdid
     &     ,type='string',name='tnu4_expl'
     &     , value='biharmonic mixing coefficient for tracers')
      CALL iom_set_file_global_attr_char   (cdid
     &     ,type='string',name='units', value='meter4 second-1')
# endif
# if !defined LMD_MIXING && !defined BVF_MIXING
!
! Background vertical viscosity and tracer mixing coefficients.
!
# endif
#endif
!
! Bottom drag coefficients.
!
C      if (maxval(Zob).ne.0.) then
C 		  elseif (rdrg2.gt.0.) then
C 			  elseif (rdrg.ne.0) then
 
      CALL iom_set_file_global_attr_float   (cdid
     &     ,type='float',name='rdrg', value=rdrg)
      CALL iom_set_file_global_attr_char   (cdid
     &     ,type='string',name='rdrg_expl'
     &     , value='linear drag coefficient')
      CALL iom_set_file_global_attr_char   (cdid
     &     ,type='string',name='rdrg_units'
     &     , value='meter second-1')
C 				  endif
#ifdef SOLVE3D
!
! Equation of State parameters.
!
      CALL iom_set_file_global_attr_float   (cdid
     &     ,type='float',name='rho0', value=rho0)
      CALL iom_set_file_global_attr_char   (cdid
     &     ,type='string',name='rho0_expl'
     &     ,value='Mean density used in Boussinesq approximation')
      CALL iom_set_file_global_attr_char   (cdid
     &     ,type='string',name='rho0_units'
     &     ,value='kilogram meter-3')
# ifndef NONLIN_EOS
# endif

! Various parameters.
!
# ifdef BODYFORCE
# endif
#endif /* SOLVE3D */
!
! Slipperiness parameters.
!
      CALL iom_set_file_global_attr_float   (cdid
     &     ,type='float',name='gamma2', value=gamma2)
      CALL iom_set_file_global_attr_char   (cdid
     &     ,type='string',name='gamma2_expl', value=
     &          'Slipperiness parameter')
# ifdef SPONGE
!
! Sponge parameters
!
      CALL iom_set_file_global_attr_float   (cdid
     &      ,type='float',name='x_sponge', value=x_sponge)
      CALL iom_set_file_global_attr_float   (cdid
     &     ,type='float',name='v_sponge', value=v_sponge)
      CALL iom_set_file_global_attr_char   (cdid,
     &     type='string',name='sponge_expl',value=
     &     'Sponge parameters : extent (m) & viscosity (m2.s-1)')
# endif
# ifdef SEDIMENT
!
! Sediment parameters
!
# endif
!
! List of Source Codes and Activated CPP-switches
!
      CALL iom_set_file_global_attr_char   (cdid
     &     ,type='string',name='SRCS', value=srcs)
      CALL iom_set_file_global_attr_char   (cdid
     &     ,type='string',name='CPP-options', value=Coptions)

      END SUBROUTINE set_xml_global_att

      SUBROUTINE set_xmlatt
# include "param.h"
# include "scalars.h"

      !!----------------------------------------------------------------------
      !!                     ***  ROUTINE set_xmlatt  ***
      !!
      !! ** Purpose :   automatic definitions of some of the xml attributs...
      !!
      !! ** Method  :   
      !!
      !! ** History :   S. Theetten, adaptation for XIOS-2.5 (May, 2019)
      !!                R. Benshila, adaptation for CROCO 
      !!                S. Masson  , original version (NEMO)  
      !!
      !!----------------------------------------------------------------------
      CHARACTER(len=1),DIMENSION( 3) ::   clgrd                    ! suffix name
      CHARACTER(len=256)             ::   clsuff                   ! suffix name
      CHARACTER(len=1)               ::   cl1                      ! 1 character
      CHARACTER(len=2)               ::   cl2                      ! 2 characters
      CHARACTER(len=3)               ::   cl3                      ! 3 characters
      CHARACTER(LEN=20)              ::   char_value
      INTEGER                        ::   ji, jg                   ! loop counters
      INTEGER                        ::   ix, iy                   ! i-,j- index
      REAL       ,DIMENSION(11) ::   zlontao                  ! longitudes of tao    moorings
      REAL      ,DIMENSION( 7) ::   zlattao                  ! latitudes  of tao    moorings
      REAL        ,DIMENSION( 4) ::   zlonrama                 ! longitudes of rama   moorings
      REAL       ,DIMENSION(11) ::   zlatrama                 ! latitudes  of rama   moorings
      REAL        ,DIMENSION( 3) ::   zlonpira                 ! longitudes of pirata moorings
      REAL     ,DIMENSION( 9) ::   zlatpira                 ! latitudes  of pirata moorings
!!!!!!!!!!!!!!!!#if  defined XIOS2
!!!!!!!!!!!!!!!!      TYPE(xios_duration)            ::   f_op, f_of
!!!!!!!!!!!!!!!!#endif
      !!----------------------------------------------------------------------
      !
      ! frequency of the call of iom_put (attribut: freq_op)
      WRITE(cl1,'(i1)') 1 
!!!!!!!!!!!!!!!!!!      CALL iom_set_field_attr('field_definition', 
!!!!!!!!!!!!!!!!!!     & freq_op = cl1//'ts', freq_offset='0ts')

      ! output file names (attribut: name)
      DO ji = 1, 9
         WRITE(cl1,'(i1)') ji
         CALL iom_update_file_name('file'//cl1)
      END DO
      DO ji = 1, 99
         WRITE(cl2,'(i2.2)') ji
         CALL iom_update_file_name('file'//cl2)
      END DO
      DO ji = 1, 999
         WRITE(cl3,'(i3.3)') ji
         CALL iom_update_file_name('file'//cl3)
      END DO

      ! Zooms...
      clgrd = (/ 'T', 'U', 'W' /)
      DO jg = 1, SIZE(clgrd)     ! grid type
         cl1 = clgrd(jg)
         ! Equatorial section (attributs: jbegin, ni, name_suffix)
         CALL dom_ngb( 0., 0., ix, iy, cl1 )
         CALL iom_set_domain_attr ('Eq'//cl1, 
     &    zoom_jbegin=iy, zoom_ni=LLm0)
         CALL iom_get_file_attr   ('Eq'//cl1, 
     &    name_suffix = clsuff             )
         CALL iom_set_file_attr   ('Eq'//cl1, 
     &    name_suffix = TRIM(clsuff)//'_Eq')
         CALL iom_update_file_name('Eq'//cl1)
      END DO
      ! TAO moorings (attributs: ibegin, jbegin, name_suffix)
      zlontao = (/ 137.0, 147.0, 156.0, 165.0, -180.0, -170.0,
     &  -155.0, -140.0, -125.0, -110.0, -95.0 /)
      zlattao = (/  -8.0,  -5.0,  -2.0,   0.0,    2.0,    5.0,
     &     8.0 /)
      CALL set_mooring( zlontao, zlattao )
      ! RAMA moorings (attributs: ibegin, jbegin, name_suffix)
      zlonrama = (/  55.0,  67.0, 80.5, 90.0 /)
      zlatrama = (/ -16.0, -12.0, -8.0, -4.0, -1.5, 0.0, 1.5,
     &  4.0, 8.0, 12.0, 15.0 /)
      CALL set_mooring( zlonrama, zlatrama )
      ! PIRATA moorings (attributs: ibegin, jbegin, name_suffix)
      zlonpira = (/ -38.0, -23.0, -10.0 /)
      zlatpira = (/ -19.0, -14.0,  -8.0, 0.0, 4.0, 8.0, 12.0, 15.0, 20.0 /)
      CALL set_mooring( zlonpira, zlatpira )




      !
      END SUBROUTINE set_xmlatt


      SUBROUTINE set_mooring( plon, plat)
      !!----------------------------------------------------------------------
      !!                     ***  ROUTINE set_mooring  ***
      !!
      !! ** Purpose :   automatic definitions of moorings xml attributs...
      !!
      !! ** Method  :   provide a list of coordinates
      !!
      !! ** History :   R. Benshila, adaptation for CROCO 
      !!                S. Masson  , original version (NEMO)  
      !!
      !!----------------------------------------------------------------------
      REAL, DIMENSION(:), INTENT(in) ::  plon, plat           ! longitudes/latitudes oft the mooring
      !
      CHARACTER(len=1),DIMENSION(1) ::   clgrd = (/ 'T' /)        ! suffix name
      CHARACTER(len=256)            ::   clname                   ! file name
      CHARACTER(len=256)            ::   clsuff                   ! suffix name
      CHARACTER(len=1)              ::   cl1                      ! 1 character
      CHARACTER(len=6)              ::   clon,clat                ! name of longitude, latitude
      INTEGER                       ::   ji, jj, jg               ! loop counters
      INTEGER                       ::   ix, iy                   ! i-,j- index
      REAL                          ::   zlon, zlat
      !!----------------------------------------------------------------------
      !
      DO jg = 1, SIZE(clgrd)
         cl1 = clgrd(jg)
         DO ji = 1, SIZE(plon)
            DO jj = 1, SIZE(plat)
               zlon = plon(ji)
               zlat = plat(jj)
               ! modifications for RAMA moorings
               IF( zlon ==  67. .AND. zlat ==  15. )   zlon =  65.
               IF( zlon ==  90. .AND. zlat <=  -4. )   zlon =  95.
               IF( zlon ==  95. .AND. zlat ==  -4. )   zlat =  -5.
               ! modifications for PIRATA moorings
               IF( zlon == -38. .AND. zlat == -19. )   zlon = -34.
               IF( zlon == -38. .AND. zlat == -14. )   zlon = -32.
               IF( zlon == -38. .AND. zlat ==  -8. )   zlon = -30.
               IF( zlon == -38. .AND. zlat ==   0. )   zlon = -35.
               IF( zlon == -23. .AND. zlat ==  20. )   zlat =  21.
               IF( zlon == -10. .AND. zlat == -14. )   zlat = -10.
               IF( zlon == -10. .AND. zlat ==  -8. )   zlat =  -6.
               IF( zlon == -10. .AND. zlat ==   4. ) THEN  
                   zlon = 0.   ;   zlat = 0.   
               ENDIF
               CALL dom_ngb( zlon, zlat, ix, iy, cl1 )
               IF( zlon >= 0. ) THEN  
                  IF( zlon == REAL(NINT(zlon)) ) THEN   
                  WRITE(clon, '(i3,  a)') NINT( zlon), 'e'
                  ELSE                                      
                  WRITE(clon, '(f5.1,a)')       zlon , 'e'
                  ENDIF
               ELSE             
                  IF( zlon == REAL(NINT(zlon)) ) THEN  
                  WRITE(clon, '(i3,  a)') NINT(-zlon), 'w'
                  ELSE                                     
                  WRITE(clon, '(f5.1,a)')      -zlon , 'w'
                  ENDIF
               ENDIF
               IF( zlat >= 0. ) THEN  
                  IF( zlat == REAL(NINT(zlat)) ) THEN  
                  WRITE(clat, '(i2,  a)') NINT( zlat), 'n'
                  ELSE                                      
                  WRITE(clat, '(f4.1,a)')       zlat , 'n'
                  ENDIF
               ELSE             
                  IF( zlat == REAL(NINT(zlat)) ) THEN   
                  WRITE(clat, '(i2,  a)') NINT(-zlat), 's'
                  ELSE                                      
                  WRITE(clat, '(f4.1,a)')      -zlat , 's'
                  ENDIF
               ENDIF
               clname = TRIM(ADJUSTL(clat))//TRIM(ADJUSTL(clon))
               CALL iom_set_domain_attr 
     &           (TRIM(clname)//cl1, zoom_ibegin= ix, zoom_jbegin= iy)
               CALL iom_get_file_attr   
     &           (TRIM(clname)//cl1, name_suffix = clsuff)
               CALL iom_set_file_attr   
     &           (TRIM(clname)//cl1, name_suffix = TRIM(clsuff)
     &            //'_'//TRIM(clname))
               CALL iom_update_file_name
     &           (TRIM(clname)//cl1)
            END DO
         END DO
      END DO
      ! 
      END SUBROUTINE set_mooring


      SUBROUTINE dom_ngb( plon, plat, kii, kjj, cdgrid )
# include "param.h"
# include "grid.h"
# ifdef MPI
#      include "mpi_cpl.h"
       include 'mpif.h'
# endif
      !!----------------------------------------------------------------------
      !!                    ***  ROUTINE dom_ngb  ***
      !!
      !! ** Purpose :   find the closest grid point from a given lon/lat position
      !!
      !! ** Method  :   look for minimum distance in cylindrical projection 
      !!                -> not good if located at too high latitude...
      !!
      !!** History  : R. Benshila, adaptation for CROCO 
      !!              S. Masson  , original version (NEMO)  
      !!
      !----------------------------------------------------------------------
      REAL            , INTENT(in   ) ::   plon, plat   ! longitude,latitude of the point
      INTEGER         , INTENT(  out) ::   kii, kjj     ! i-,j-index of the closes grid point
      CHARACTER(len=1), INTENT(in   ) ::   cdgrid       ! grid name 'T', 'U', 'V', 'W'
      !
      INTEGER , DIMENSION(2) ::   iloc
      REAL                   ::   zlon, zmini
      REAL,DIMENSION(GLOBAL_2D_ARRAY) ::  zglam, zgphi, 
     &                                    zmask, zdist
      REAL, DIMENSION(2,1) ::   zain,zaout
      INTEGER :: ierror
      !!--------------------------------------------------------------------
      !
      zmask(:,:) = 0.
# ifdef SPHERICAL
      SELECT CASE( cdgrid )
      CASE( 'U' ) 
      zglam(:,:) = lonu(:,:) 
      zgphi(:,:) = latu(:,:)
#  ifdef MASKING
      zmask(:,:) = umask(:,:)
#  endif
      CASE( 'V' ) 
      zglam(:,:) = lonv(:,:) 
      zgphi(:,:) = latv(:,:) 
#  ifdef MASKING
      zmask(:,:) = vmask(:,:)
#  endif
      CASE DEFAULT 
      zglam(:,:) = lonr(:,:) 
      zgphi(:,:) = latr(:,:) ; 
#  ifdef MASKING
      zmask(:,:) = rmask(:,:)
#  endif
      END SELECT
# else
      SELECT CASE( cdgrid )
      CASE( 'U' ) 
      zglam(:,:) = yr(:,:) 
      zgphi(:,:) = xp(:,:)
#  ifdef MASKING
      zmask(:,:) = umask(:,:)
#  endif
      CASE( 'V' ) 
      zglam(:,:) = yp(:,:) 
      zgphi(:,:) = xr(:,:) 
# ifdef MASKING
      zmask(:,:) = vmask(:,:)
# endif
      CASE DEFAULT 
      zglam(:,:) = yr(:,:) 
      zgphi(:,:) = xr(:,:) ; 
#  ifdef MASKING
      zmask(:,:) = rmask(:,:)
#  endif
      END SELECT
# endif

# ifdef SPHERICAL
      zlon       = MOD( plon       + 720., 360. )                                     ! plon between    0 and 360
      zglam(:,:) = MOD( zglam(:,:) + 720., 360. )                                     ! glam between    0 and 360
      IF( zlon > 270. )   zlon = zlon - 360.                                          ! zlon between  -90 and 270
      IF( zlon <  90. )   WHERE( zglam(:,:) > 180. ) zglam(:,:) = zglam(:,:) - 360.   ! glam between -180 and 180
# endif
      zglam(:,:) = zglam(:,:) - zlon
      zgphi(:,:) = zgphi(:,:) - plat
      zdist(:,:) = zglam(:,:) * zglam(:,:) + zgphi(:,:) * zgphi(:,:)
      
# ifdef MPI
      zmini = MINVAL( zdist(:,:) , mask= zmask == 1.e0 )
      iloc = MINLOC( zdist(:,:) , mask= zmask == 1.e0 )
      !
      kii = iloc(1) + iminmpi - 1
      kjj = iloc(2) + jminmpi - 1
      !
      zain(1,:)=zmini
      zain(2,:)=kii+10000.*kjj
      !
      CALL MPI_ALLREDUCE( zain,zaout, 1, MPI_2DOUBLE_PRECISION,
     &                    MPI_MINLOC,MPI_COMM_WORLD,ierror)
      !
      zmini = zaout(1,1)
      kii = INT(zaout(2,1)/10000.)
      kjj = INT(zaout(2,1) - 10000.*kjj )
# else
         iloc(:) = MINLOC( zdist(:,:), mask = zmask(:,:) == 1.e0 )
         kii = iloc(1) - 1
         kjj = iloc(2) - 1
# endif
      !
      END SUBROUTINE dom_ngb


      SUBROUTINE iom_update_file_name( cdid )
# include "param.h"
# include "scalars.h"
# include "ncscrum.h"
# include "strings.h"
      !!----------------------------------------------------------------------
      !!                    ***  ROUTINE iom_update_file_name  ***
      !!
      !! ** Purpose :   update output file name from model information
      !!
      !! ** Method  :   request xml attribute and modify it
      !!
      !! ** History :   S. Theetten, adaptation for XIOS-2.5 (May, 2019)
      !!                R. Benshila, adaptation for CROCO 
      !!                S. Masson  , original version (NEMO)  
      !!
      !!---------------------------------------------------------------------
      CHARACTER(LEN=*)          , INTENT(in) ::   cdid
      !
      CHARACTER(LEN=256) ::   clname
      CHARACTER(LEN=20)  ::   clfreq
      CHARACTER(LEN=20)  ::   cldate
      CHARACTER(LEN=20)  ::   type, name, value, nom, valeur
      CHARACTER(LEN=20)  ::   char_value
      INTEGER            ::   idx
      INTEGER            ::   jn
      INTEGER            ::   itrlen
      INTEGER            ::   lvar,lenstr
      INTEGER            ::   iyear, imonth, iday, isec
      REAL               ::   zsec, jday0,rvalue
      LOGICAL            ::   llexist
# if  defined XIOS2
      TYPE(xios_duration)   ::   output_freq
# endif     
      !!----------------------------------------------------------------------
      !
# ifdef START_DATE
      ! start_date with format dd-mm-yyyy
       READ(start_date(1:2),fmt='(i2)') iday
       READ(start_date(4:5),fmt='(i2)') imonth
       READ(start_date(7:10),fmt='(i4)') iyear
       CALL ymds2ju( iyear, imonth,iday, 0., jday0 )  
# elif defined USE_CALENDAR
       CALL ymds2ju(1900, 1, 1, 0., jday0 )  
# else
       CALL ymds2ju( 1, 1, 1, 0., jday0 )  
# endif
      DO jn = 1,2

# if ! defined XIOS2
         IF( jn == 1 )   CALL iom_get_file_attr( cdid, name 
     &           = clname, output_freq = clfreq )
# else
         IF( jn == 1 )   CALL iom_get_file_attr( cdid, name 
     &           = clname, output_freq = output_freq )
# endif
         IF( jn == 2 )   CALL iom_get_file_attr( cdid, name_suffix = clname )

         IF ( TRIM(clname) /= '' ) THEN 

            idx = INDEX(clname,'@expname@') + INDEX(clname,'@EXPNAME@')
            DO WHILE ( idx /= 0 ) 
               clname = clname(1:idx-1)//TRIM(replace_blank(TRIM(title),"_"))
     &          //clname(idx+9:LEN_TRIM(clname))
               idx = INDEX(clname,'@expname@') + INDEX(clname,'@EXPNAME@')
            END DO

# if ! defined XIOS2
            idx = INDEX(clname,'@freq@') + INDEX(clname,'@FREQ@')
            DO WHILE ( idx /= 0 ) 
               IF ( TRIM(clfreq) /= '' ) THEN
                  itrlen = LEN_TRIM(clfreq)
                  IF ( clfreq(itrlen-1:itrlen) == 'mo' ) THEN
                  clfreq = clfreq(1:itrlen-1)
                  ENDIF
                  clname = clname(1:idx-1)//TRIM(clfreq)
     &                    //clname(idx+6:LEN_TRIM(clname))
               ENDIF
               idx = INDEX(clname,'@freq@') + INDEX(clname,'@FREQ@')
            END DO
# else
            idx = INDEX(clname,'@freq@') + INDEX(clname,'@FREQ@')
            DO WHILE ( idx /= 0 )
              IF ( output_freq%timestep /= 0) THEN
                  WRITE(clfreq,'(I18,A2)')INT(output_freq%timestep),'ts'
                  itrlen = LEN_TRIM(ADJUSTL(clfreq))
              ELSE IF ( output_freq%hour /= 0 ) THEN
                  WRITE(clfreq,'(I19,A1)')INT(output_freq%hour),'h'
                  itrlen = LEN_TRIM(ADJUSTL(clfreq))
              ELSE IF ( output_freq%day /= 0 ) THEN
                  WRITE(clfreq,'(I19,A1)')INT(output_freq%day),'d'
                  itrlen = LEN_TRIM(ADJUSTL(clfreq))
              ELSE IF ( output_freq%month /= 0 ) THEN
                  WRITE(clfreq,'(I19,A1)')INT(output_freq%month),'m'
                  itrlen = LEN_TRIM(ADJUSTL(clfreq))
              ELSE IF ( output_freq%year /= 0 ) THEN
                  WRITE(clfreq,'(I19,A1)')INT(output_freq%year),'y'
                  itrlen = LEN_TRIM(ADJUSTL(clfreq))
              ELSE
                  if(mynode .eq. 0 ) print*,'error in the name of file id '
     &           //TRIM(cdid),
     & ' attribute output_freq is undefined -> cannot replace @freq@ in '
     & //TRIM(clname) 
				  STOP
!                  CALL ctl_stop('error in the name of file id '
!     &			  //TRIM(cdid),   
!     & ' attribute output_freq is undefined -> cannot replace @freq@ in '
!     & //TRIM(clname) )
              ENDIF
              clname = clname(1:idx-1)//TRIM(ADJUSTL(clfreq))//
     &        clname(idx+6:LEN_TRIM(clname))
              idx = INDEX(clname,'@freq@') + INDEX(clname,'@FREQ@')
            END DO
# endif
            idx = INDEX(clname,'@startdate@') + INDEX(clname,'@STARTDATE@')
            DO WHILE ( idx /= 0 ) 
               cldate = iom_sdate( tdays+jday0  )
               clname = clname(1:idx-1)//TRIM(cldate)//
     &          clname(idx+11:LEN_TRIM(clname))
               idx = INDEX(clname,'@startdate@') + INDEX(clname,'@STARTDATE@')
            END DO

            idx = INDEX(clname,'@startdatefull@') 
     &       + INDEX(clname,'@STARTDATEFULL@')
            DO WHILE ( idx /= 0 ) 
               cldate = iom_sdate(tdays+jday0 , ldfull = .TRUE. )
               clname = clname(1:idx-1)//TRIM(cldate)//
     &          clname(idx+15:LEN_TRIM(clname))
               idx = INDEX(clname,'@startdatefull@')
     &           + INDEX(clname,'@STARTDATEFULL@')
            END DO

            idx = INDEX(clname,'@enddate@') + INDEX(clname,'@ENDDATE@')
            DO WHILE ( idx /= 0 ) 
               cldate = iom_sdate( tdays+jday0 + dt
     &          / day2sec * REAL( ntimes - ntstart ), ld24 = .TRUE. )
               clname = clname(1:idx-1)//TRIM(cldate)
     &          //clname(idx+9:LEN_TRIM(clname))
               idx = INDEX(clname,'@enddate@') + INDEX(clname,'@ENDDATE@')
            END DO

            idx = INDEX(clname,'@enddatefull@')
     &        + INDEX(clname,'@ENDDATEFULL@')
            DO WHILE ( idx /= 0 ) 
               cldate = iom_sdate( tdays+jday0 + dt  
     &          /day2sec * REAL( ntimes - ntstart ), ld24 = .TRUE., 
     &           ldfull = .TRUE. )
               clname = clname(1:idx-1)//TRIM(cldate)
     &          //clname(idx+13:LEN_TRIM(clname))
               idx = INDEX(clname,'@enddatefull@')
     &           + INDEX(clname,'@ENDDATEFULL@')
            END DO
            IF( jn == 1 )   CALL iom_set_file_attr( cdid, name
     &               = clname )
            IF( jn == 2 )   CALL iom_set_file_attr( cdid, name_suffix
     &               = clname )
	 
	        CALL set_xml_global_att(cdid)

         ENDIF

      END DO
      !
      END SUBROUTINE iom_update_file_name


      FUNCTION iom_sdate( pjday, ld24, ldfull )
      !!----------------------------------------------------------------------
      !!                     ***  ROUTINE iom_sdate  ***
      !!
      !! ** Purpose :   send back the date corresponding to the given julian day
      !!
      !! ** Method  :   
      !!
      !! ** History :   R. Benshila, adaptation for CROCO 
      !!                S. Masson  , original version (NEMO)  
      !!
      !!---------------------------------------------------------------------
      REAL    , INTENT(in   )           ::   pjday         ! julian day
      LOGICAL , INTENT(in   ), OPTIONAL ::   ld24          ! true to force 24:00 instead of 00:00
      LOGICAL , INTENT(in   ), OPTIONAL ::   ldfull        ! true to get the compleate date: yyyymmdd_hh:mm:ss
      !
      CHARACTER(LEN=20) ::   iom_sdate
      CHARACTER(LEN=50) ::   clfmt                         !  format used to write the date
      INTEGER           ::   iyear, imonth, iday, ihour, iminute, isec
      REAL              ::   zsec
      LOGICAL           ::   ll24, llfull
      !!---------------------------------------------------------------------
      !
      IF( PRESENT(ld24) ) THEN   ;   ll24 = ld24
      ELSE                       ;   ll24 = .FALSE.
      ENDIF

      IF( PRESENT(ldfull) ) THEN   ;   llfull = ldfull
      ELSE                         ;   llfull = .FALSE.
      ENDIF

      CALL ju2ymds( pjday, iyear, imonth, iday, zsec )
      isec = NINT(zsec)

      IF( ll24 .AND. isec == 0 ) THEN   ! 00:00 of the next day -> move to 24:00 of the current day
         CALL ju2ymds( pjday - 1., iyear, imonth, iday, zsec )
         isec = 86400
      ENDIF

      IF( iyear < 10000 ) THEN  
      clfmt = "i4.4,2i2.2"                ! format used to write the date
      ELSE                   
      WRITE(clfmt, "('i',i1,',2i2.2')") INT(LOG10(REAL(iyear))) + 1
      ENDIF

      IF( llfull ) THEN
         clfmt = TRIM(clfmt)//",'_',3i2.2"
         ihour   = isec / 3600
         isec    = MOD(isec, 3600)
         iminute = isec / 60
         isec    = MOD(isec, 60)
         WRITE(iom_sdate, '('//TRIM(clfmt)//')') 
     &    iyear, imonth, iday, ihour, iminute, isec    ! date of the end of run
      ELSE
         WRITE(iom_sdate, '('//TRIM(clfmt)//')')
     &     iyear, imonth, iday                          ! date of the end of run
      ENDIF

      END FUNCTION iom_sdate

 
      SUBROUTINE ymds2ju (year,month,day,sec,julian)
      !!----------------------------------------------------------------------
      !!                     ***  ROUTINE ymds2ju  ***
      !!
      !! ** Purpose :   send back the julian day corresponding to the given date
      !!
      !! ** Method  :   
      !!
      !! ** History :   R. Benshila, adaptation for CROCO 
      !!                IPSL       , original version  
      !!
      !!---------------------------------------------------------------------
      IMPLICIT NONE
      !
      INTEGER,INTENT(IN) :: year,month,day
      REAL,INTENT(IN)    :: sec
      REAL,INTENT(OUT) :: julian
      !
      INTEGER :: julian_day
      REAL    :: julian_sec
      REAL,PARAMETER :: one_day = 86400.0    
      !---------------------------------------------------------------------
      !
      CALL ymds2ju_internal (year,month,day,sec,julian_day,julian_sec)
      julian = julian_day+julian_sec/one_day
      !
      END SUBROUTINE ymds2ju

      SUBROUTINE ymds2ju_internal (year,month,day,sec,julian_day,julian_sec)
      !---------------------------------------------------------------------
      !- Converts year, month, day and seconds into a julian day
      !-
      !- In 1968 in a letter to the editor of Communications of the ACM
      !- (CACM, volume 11, number 10, October 1968, p.657) Henry F. Fliegel
      !- and Thomas C. Van Flandern presented such an algorithm.
      !-
      !- See also : http://www.magnet.ch/serendipity/hermetic/cal_stud/jdn.htm
      !-
      !- In the case of the Gregorian calendar we have chosen to use
      !- the Lilian day numbers. This is the day counter which starts
      !- on the 15th October 1582.
      !- This is the day at which Pope Gregory XIII introduced the
      !- Gregorian calendar.
      !- Compared to the true Julian calendar, which starts some
      !- 7980 years ago, the Lilian days are smaler and are dealt with
      !- easily on 32 bit machines. With the true Julian days you can only
      !- the fraction of the day in the real part to a precision of
      !- a 1/4 of a day with 32 bits.
      !---------------------------------------------------------------------
      IMPLICIT NONE
      !
      INTEGER,INTENT(IN) :: year,month,day
      REAL,INTENT(IN)    :: sec
      INTEGER,INTENT(OUT) :: julian_day
      REAL,INTENT(OUT)    :: julian_sec
      !
      REAL,PARAMETER :: one_day = 86400.0    
      REAL,PARAMETER :: one_year = 365.2425
      INTEGER :: mon_len(12)=(/31,28,31,30,3
     &     1,30,31,31,30,31,30,31/)
      INTEGER :: jd,m,y,d,ml
      !---------------------------------------------------------------------
      !
      m = month
      y = year
      d = day
      !
      ! We deduce the calendar from the length of the year as it
      ! is faster than an INDEX on the calendar variable.
      !
      IF ( (one_year > 365.0).AND.(one_year < 366.0) ) THEN
      !-- "Gregorian"
       jd = (1461*(y+4800+INT((m-14)/12)))/4 
     &      +(367*(m-2-12*(INT((m-14)/12))))/12 
     &      -(3*((y+4900+INT((m-14)/12))/100))/4 
     &      +d-32075
       jd = jd-2299160
      ELSE IF (    (ABS(one_year-365.0) <= EPSILON(one_year))  
     &  .OR.(ABS(one_year-366.0) <= EPSILON(one_year)) ) THEN
      !-- "No leap" or "All leap"
      ml = SUM(mon_len(1:m-1))
      jd = y*NINT(one_year)+ml+(d-1)
      ELSE
      !-- Calendar with regular month
       ml = NINT(one_year/12.)
       jd = y*NINT(one_year)+(m-1)*ml+(d-1)
      ENDIF
      !
      julian_day = jd
      julian_sec = sec
      !
      END SUBROUTINE ymds2ju_internal

     
      SUBROUTINE ju2ymds (julian,year,month,day,sec)
      !!----------------------------------------------------------------------
      !!                     ***  ROUTINE ymds2ju  ***
      !!
      !! ** Purpose :   send back the date corresponding to the given julian day
      !!
      !! ** Method  :   
      !!
      !! ** History :   R. Benshila, adaptation for CROCO 
      !!                IPSL       , original version  
      !!
      !!---------------------------------------------------------------------
	  IMPLICIT NONE
	  !
	  REAL,INTENT(IN) :: julian
      INTEGER,INTENT(OUT) :: year,month,day
	  REAL,INTENT(OUT)    :: sec
	  !
	  INTEGER :: julian_day
	  REAL    :: julian_sec
      REAL,PARAMETER :: one_day = 86400.0	
      !---------------------------------------------------------------------
	  ! 
	  julian_day = INT(julian)
	  julian_sec = (julian-julian_day)*one_day
      !
      CALL ju2ymds_internal(julian_day,julian_sec,year,month,day,sec)
      !
	  END SUBROUTINE ju2ymds


      SUBROUTINE ju2ymds_internal (julian_day,julian_sec,year,
     &  month,day,sec)
	  !---------------------------------------------------------------------
	  !- This subroutine computes from the julian day the year,
	  !- month, day and seconds
      !-
	  !- In 1968 in a letter to the editor of Communications of the ACM
	  !- (CACM, volume 11, number 10, October 1968, p.657) Henry F. Fliegel
	  !- and Thomas C. Van Flandern presented such an algorithm.
	  !-
	  !- See also : http://www.magnet.ch/serendipity/hermetic/cal_stud/jdn.htm
	  !-
	  !- In the case of the Gregorian calendar we have chosen to use
	  !- the Lilian day numbers. This is the day counter which starts
      !- on the 15th October 1582. This is the day at which Pope
	  !- Gregory XIII introduced the Gregorian calendar.
	  !- Compared to the true Julian calendar, which starts some 7980
	  !- years ago, the Lilian days are smaler and are dealt with easily
      !- on 32 bit machines. With the true Julian days you can only the
	  !- fraction of the day in the real part to a precision of a 1/4 of
	  !- a day with 32 bits.
	  !---------------------------------------------------------------------
	  IMPLICIT NONE
      !
	  INTEGER,INTENT(IN) :: julian_day
	  REAL,INTENT(IN)    :: julian_sec
	  INTEGER,INTENT(OUT) :: year,month,day
	  REAL,INTENT(OUT)    :: sec
	  !
	  INTEGER :: l,n,i,jd,j,d,m,y,ml
	  INTEGER :: add_day
	  REAL :: eps_day
      REAL,PARAMETER :: one_day = 86400.0	
      REAL,PARAMETER :: one_year = 365.2425
      INTEGER :: mon_len(12)=(/31,28,31,30,3
     &     1,30,31,31,30,31,30,31/)
      !---------------------------------------------------------------------
	  !	 
	  eps_day = SPACING(one_day)
	  !
	  jd = julian_day
	  sec = julian_sec
	  IF (sec > (one_day-eps_day)) THEN
	    add_day = INT(sec/one_day)
	    sec = sec-add_day*one_day
	    jd = jd+add_day
	  ENDIF
	  IF (sec < -eps_day) THEN
		  sec = sec+one_day
		  jd = jd-1
	  ENDIF
	  !
      IF ( (one_year > 365.0).AND.(one_year < 366.0) ) THEN
	      !-- Gregorian
	     jd = jd+2299160
	     !
	     l = jd+68569
         n = (4*l)/146097
         l = l-(146097*n+3)/4
	     i = (4000*(l+1))/1461001
	     l = l-(1461*i)/4+31
	     j = (80*l)/2447
	     d = l-(2447*j)/80
	     l = j/11
		 m = j+2-(12*l)
	     y = 100*(n-49)+i+l
	  ELSE IF (    (ABS(one_year-365.0) <= EPSILON(one_year)) 
     &  .OR.(ABS(one_year-366.0) <= EPSILON(one_year)) ) THEN
	     !-- No leap or All leap
	     y = jd/NINT(one_year)
		 l = jd-y*NINT(one_year)
         m = 1
	     ml = 0
	     DO WHILE (ml+mon_len(m) <= l)
		    ml = ml+mon_len(m)
		    m = m+1
	     ENDDO
		 d = l-ml+1
      ELSE
	     !-- others
		 ml = NINT(one_year/12.)
		 y = jd/NINT(one_year)
		 l = jd-y*NINT(one_year)
		 m = (l/ml)+1
		 d = l-(m-1)*ml+1
      ENDIF
	  !
	  day = d
      month = m
      year = y
      !
	  END SUBROUTINE ju2ymds_internal

   
      FUNCTION replace_blank(string,rep) RESULT(outs)
      !!----------------------------------------------------------------------
      !!                     ***  FUNCTION replace_blank  ***
      !!
      !! ** Purpose :   replace blank strin by a given pattern
      !!
      !! ** Method  :   
      !!
      !! ** History :   R. Benshila
      !!
      !!---------------------------------------------------------------------
      !
      CHARACTER(len=*) :: string, rep
      CHARACTER(len=100) :: outs
      !!---------------------------------------------------------------------
      INTEGER :: stringLen 
      INTEGER :: last, actual
      !!---------------------------------------------------------------------
      ! 
      outs=string
      stringLen = LEN(string)
      last = 1
      actual = 1

      DO WHILE (actual < stringLen)
        IF (string(last:last) == ' ') THEN
            actual = actual + 1
            outs(last:last) = '_'
            last=actual
        ELSE
            last = last + 1
            IF (actual < last) actual = last
        ENDIF
      END DO
      !
      END FUNCTION replace_blank

		    
      END MODULE xios_module

#include "cppdefs.h"
       SUBROUTINE init_xios(tile)

      IMPLICIT NONE
# include "param.h"
      INTEGER tile
# include "compute_tile_bounds.h"
      
      CALL init_xios_tile (Istr,Iend,Jstr,Jend)
      !
      END SUBROUTINE init_xios

      SUBROUTINE init_xios_tile(Istr,Iend,Jstr,Jend)
      !!----------------------------------------------------------------------
      !!                     ***  ROUTINE iom_sdate  ***
      !!
      !! ** Purpose :   definition of grid for output
      !!
      !! ** Method  :   
      !!
      !! ** History :   S. Theetten, adaptation for XIOS-2.5 (May, 2019)
      !!                R. Benshila, improvement for CROCO 
      !!                A. Ponte   , adaption for ROMS  
      !!                S. Masson  , original version (NEMO)  
      !!
      !!---------------------------------------------------------------------
      !
      USE xios     ! XIOS library
      USE xios_module
      !    
      IMPLICIT NONE
      !
# include "param.h"
# include "grid.h"
# include "ocean3d.h"
# include "scalars.h"
# include "ncscrum.h"
# include "mpi_cpl.h"
      !
      INTEGER Istr,Iend,Jstr,Jend
      !
      INTEGER tile,ji
      INTEGER ilocal_comm
# if !  defined XIOS2
      TYPE(xios_time) :: dtime  = xios_time(0, 0, 0, 0, 0, 0)
# else
      TYPE(xios_duration) :: dtime    = xios_duration(0,
     &  0, 0, 0, 0, 0)
!      TYPE(xios_date)     :: start_date
# endif     
      TYPE(xios_context) :: ctx_hdl
      TYPE(xios_field) :: field_hdl
      TYPE(xios_fieldgroup) :: fieldgroup_hdl
      TYPE(xios_file) :: file_hdl
      TYPE(xios_variable) :: var_hdl
      LOGICAL ok
      !
      CHARACTER(len=19) :: cldate
      CHARACTER(len=10) :: clname
      CHARACTER(len=10) :: cdname
      INTEGER iday,imonth,iyear,ihour,iminute,isec
#ifdef USE_CALENDAR
      CHARACTER (len=19)    :: tool_sectodat
#endif

      !!---------------------------------------------------------------------
      !  
# include "compute_auxiliary_bounds.h"
      !
# ifdef EW_PERIODIC
#  define IU_RANGE Istr,Iend
#  define IV_RANGE Istr,Iend
# else
#  define IU_RANGE Istr,IendR
#  define IV_RANGE IstrR,IendR
# endif

# ifdef NS_PERIODIC
#  define JU_RANGE Jstr,Jend
#  define JV_RANGE Jstr,Jend
# else
#  define JU_RANGE JstrR,JendR
#  define JV_RANGE Jstr,JendR
# endif
      !
      cdname='roms'
      clname = cdname
# ifdef AGRIF
      IF( TRIM(Agrif_CFixed()) /= '0' )
     & clname =TRIM(cdname)//"."//TRIM(Agrif_CFixed())
# endif
# if defined MPI
      CALL xios_context_initialize(TRIM(clname), MPI_COMM_WORLD)
# else
      CALL xios_context_initialize(TRIM(clname), 0)
# endif
      CALL iom_swap( cdname )

      ! calendar parameters
# if ! defined XIOS2	  
      CALL xios_set_context_attr(TRIM(clname), calendar_type= 
     &      "Gregorian")
# endif

# ifdef START_DATE
      ! start_date with format dd-mm-yyyy
      READ(start_date(1:2),fmt='(i2)') iday
      READ(start_date(4:5),fmt='(i2)') imonth
      READ(start_date(7:10),fmt='(i4)') iyear
      ihour=0
      iminute=0
      isec=0
# elif defined USE_CALENDAR
      cldate = tool_sectodat(time)
      CALL tool_decompdate(cldate,iday,imonth,iyear,ihour,iminute,isec)
# else
      iday=1
      imonth=1       
      iyear=1
      ihour=0
      iminute=0
      isec=0
# endif
# if ! defined XIOS2
      WRITE(cldate,"(i4.4,'-',i2.2,'-',i2.2,' 00:00:00')")
     & iyear,imonth,iday
      CALL xios_set_context_attr(TRIM(clname), 
     &        start_date=cldate )
# else
      CALL xios_define_calendar( TYPE = "Gregorian", 
     & time_origin = xios_date(1979,01,01,00,00,00),
     & start_date = xios_date(iyear,imonth,iday,ihour,iminute,isec))
#endif
      ! horizontal grid definition
      CALL set_scalar

# ifdef SPHERICAL
      CALL set_grid( "rho", lonr, latr )
      CALL set_grid( "b", lonr, latr )
      CALL set_grid( "u", lonu, latu )
      CALL set_grid( "v", lonv, latv )
      CALL set_grid( "w", lonr, latr )
      CALL set_grid( "psi", lonu, latv )
# else
      CALL set_grid( "rho", xr, yr )
      CALL set_grid( "b", xr, yr )
      CALL set_grid( "u", xp, yr )
      CALL set_grid( "v", xr, yp )
      CALL set_grid( "w", xr, yr )
      CALL set_grid( "psi", xp, yp )
# endif

# ifdef MASKING
      CALL set_mask("rho")
      CALL set_mask("b")
      CALL set_mask("u")
      CALL set_mask("v")
      CALL set_mask("w")
      CALL set_mask("psi")
# endif
 
      ! vertical grid definition
      CALL iom_set_axis_attr( "s_rho", sc_r )
      CALL iom_set_axis_attr( "s_w", sc_w )
# ifdef SEDIMENT
      CALL iom_set_axis_attr( "s_b",(/ (REAL(ji), ji=1,NLAY) /))
# else
      ! fake axis to avoid chaning xml file 
      CALL iom_set_axis_attr( "s_b",(/ (REAL(ji), ji=1,1) /))
# endif
      
      ! time axis (for ROMS tools)

      ! automatic definitions of some of the xml attributs
      CALL set_xmlatt

      ! end file definition

      ! Accademic test with dt> 1second are non comptabile with
      ! XIOS, set xios time step unit to 1 model time_step
      if (dt < 1) then
        dtime%second =  1
      else
        dtime%second = dt
      endif
      CALL xios_set_timestep(dtime)

      ! end definition phase
      CALL xios_close_context_definition()
      !
      CALL xios_update_calendar(0)

      END SUBROUTINE init_xios_tile
      
      
      SUBROUTINE iom_swap( cdname )
      !!---------------------------------------------------------------------
      !!                   ***  SUBROUTINE  iom_swap  ***
      !!
      !! ** Purpose :  swap context between different agrif grid for xmlio_server
      !!---------------------------------------------------------------------
      USE xios
      !
      CHARACTER(len=*), INTENT(in) :: cdname
      TYPE(xios_context) :: roms_hdl
# ifdef AGRIF
      IF( TRIM(Agrif_CFixed()) == '0' ) THEN
        CALL xios_get_handle(TRIM(cdname),roms_hdl)
      ELSE
        CALL xios_get_handle(TRIM(cdname)//"."//TRIM(Agrif_CFixed()),roms_hdl)
      ENDIF
# else
      CALL xios_get_handle(TRIM(cdname),roms_hdl)

# endif
      !
      CALL xios_set_current_context(roms_hdl)
      !
      END SUBROUTINE iom_swap


      SUBROUTINE iom_context_finalize( cdname )
      !!---------------------------------------------------------------------
      !!                   ***  SUBROUTINE  iom_context_finalize  ***
      !!
      !! ** Purpose :  finalization of output definition
      !!---------------------------------------------------------------------
      USE xios
      !
      CHARACTER(LEN=*), INTENT(in) :: cdname
      TYPE(xios_context) :: roms_hdl
      !
      IF( xios_is_valid_context(cdname) ) THEN
         CALL iom_swap( cdname )   ! swap to cdname context
         CALL xios_context_finalize() ! finalize the context
         IF( cdname /= "roms" ) CALL iom_swap( "roms" )   ! return back to roms context
      ENDIF
      !
      END SUBROUTINE iom_context_finalize


      subroutine send_xios_diags (tile)

      implicit none

     
# include "param.h"
# ifdef SOLVE3D
#  include "work.h"
#  include "ncscrum.h"
# endif
      integer tile
# include "compute_tile_bounds.h"

      call send_xios_diags_tile (tile,Istr,Iend,Jstr,Jend)

      return
      end


      subroutine send_xios_diags_tile (tile,Istr,Iend,Jstr,Jend)

!
!================================================== John M. Klinck ===
!  Copyright (c) 2000 Rutgers/UCLA                                   !
!================================================ Hernan G. Arango ===
!                                                                    !
!  This routine sends every diagnostics fields to Xios server         
!
!  ** History :   S. Theetten , adaptation for XIOS-2.5 (May, 2019)
!                                                                    !
!=====================================================================
!

          USE xios     ! XIOS io
          USE xios_module
#if defined BIOLOGY && defined PISCES
          USE trcwri_pisces
#endif
          IMPLICIT NONE

# include "param.h"
# include "ncscrum.h"
# include "grid.h"
# include "ocean2d.h"
# include "ocean3d.h"
# include "forces.h"
# include "mixing.h"
# include "coupling.h"
# ifdef SEDIMENT
#  include "sediment.h"
# endif
# ifdef BBL
#  include "bbl.h"
# endif
# ifdef WKB_WWAVE
#  include "wkb_wwave.h"
# endif
# ifdef DIAGNOSTICS_UV
#  include "diagnostics.h"
# endif
# if defined DIAGNOSTICS_VRT
#  include "diags_vrt.h"
# endif
# ifdef DIAGNOSTICS_EK
#  include "diags_ek.h"
# endif
# ifdef DIAGNOSTICS_PV
#  include "diags_pv.h"
# endif
# ifdef MASKING
#  define SWITCH *
# else
#  define SWITCH !
# endif
# ifdef NBQ
# include "nbq.h"
# endif
# include "scalars.h"
# include "work.h"
# include "private_scratch.h"
!
#ifdef DIAGNOSTICS_EDDY
      real workru(GLOBAL_2D_ARRAY,1:N)
      common /work3d_r/ workru
      real workrv(GLOBAL_2D_ARRAY,1:N)
      common /work3d_r/ workrv
      real GRho
#endif
      integer tile
      integer Istr,Iend,Jstr,Jend
      integer i,j,k
      integer ilc, itrc
      character*20 nametrc
      real eps
      parameter (eps=1.D-20)
      real stf_cff
      real zz(GLOBAL_2D_ARRAY,0:N)
      parameter(stf_cff=86400/0.01)
!
# include "compute_auxiliary_bounds.h"
!
#  define I_EXT_RANGE Istr-1,Iend+1
#  define J_EXT_RANGE Jstr-1,Jend+1

#  ifdef EW_PERIODIC 
#   define IU_RANGE Istr,Iend
#   define IV_RANGE Istr,Iend
#  else
#   define IU_RANGE Istr,IendR
#   define IV_RANGE IstrR,IendR
#  endif

#  ifdef NS_PERIODIC 
#   define JU_RANGE Jstr,Jend
#   define JV_RANGE Jstr,Jend
#  else
#   define JU_RANGE JstrR,JendR
#   define JV_RANGE Jstr,JendR
#  endif


!       ilc=iic   ! number of time step since start
       ilc=1+iic-ntstart   ! number of time step since restart
       ! To be conform with clasic (no xios) croco file first 
       ! time step must be writen
       if (ilc == 0) ilc=1
       CALL xios_update_calendar(ilc)

!       call xios_send_field("scrum_time",time)
       call xios_send_field("time",time)

! static attributes
# ifdef SOLVE3D
!
! Variable dimension on vertical : s_rho and s_w
!
! done automatically
!
! S-coordinate independent variables "sc_w", "sc_r" and stretching
! curves "Cs_w", "Cs_r" at W- and RHO-points.
!
!
! Ugly fix to write scalars until XIOS 2.0
!
# if ! defined XIOS2
      do k=1,N
      zz(:,:,k)=Cs_r(k)
      end do
      call xios_send_field("Cs_r",zz(:,:,1:N) )
      do k=0,N
      zz(:,:,k)=Cs_w(k)
      end do
      call xios_send_field("Cs_w",zz)
      do k=1,N
      zz(:,:,k)=sc_r(k)
      end do
      call xios_send_field("sc_r",zz(:,:,1:N) )
      do k=0,N
      zz(:,:,k)=sc_w(k)
      end do
      call xios_send_field("sc_w",zz)
      zz(:,:,1)=hc
      call xios_send_field("hc",zz(:,:,1))
      zz(:,:,1)=theta_s
      call xios_send_field("theta_s",zz(:,:,1))
      zz(:,:,1)=theta_b
      call xios_send_field("theta_b",zz(:,:,1))
      zz(:,:,1)=Tcline
      call xios_send_field("Tcline",zz(:,:,1))
      zz(:,:,1)=REAL(Vtransform)
      call xios_send_field("Vtransform",zz(:,:,1))
# else
! Beautiful  way to write scalars with XIOS-2.x
      call xios_send_field("Cs_r",Cs_r(1:N))
      call xios_send_field("Cs_w",Cs_w(0:N))
      call xios_send_field("sc_r",sc_r(1:N))
      call xios_send_field("sc_w",sc_w(0:N))
      call xios_send_field("hc",hc)
      call xios_send_field("theta_s",theta_s)
      call xios_send_field("theta_b",theta_b)
      call xios_send_field("Tcline",Tcline)
      call xios_send_field("Vtransform",REAL(Vtransform))
# endif

      call xios_send_field("levels_rho",z_r(:,:,:))
      call xios_send_field("levels_w"  ,z_w(:,:,:))
      do k=1,N
        do j=JstrR,JendR
          do i=Istr,IendR
            workr(i,j,k)=0.5*(z_r(i,j,k)+z_r(i-1,j,k))
          enddo
        enddo
      enddo
      call xios_send_field("levels_u"  ,workr(:,:,:))
      do k=1,N
       do j=JV_RANGE
          do i=IV_RANGE
            workr(i,j,k)=0.5*(z_r(i,j,k)+z_r(i,j-1,k))
          enddo
       enddo
      enddo
      call xios_send_field("levels_v"  ,workr(:,:,:))
      
# endif
!
! Bathymetry.
!
      call xios_send_field("h",h)
!
! Coriolis parameter.
!
      call xios_send_field("f",f)
!
! Curvilinear transformation metrics.
!
      call xios_send_field("pm",pm)
      call xios_send_field("pn",pn)
!
! Longitude/latitude or cartezian coordinates of RHO-points
!
# ifdef SPHERICAL
      call xios_send_field("lon_rho",lonr)
      call xios_send_field("lat_rho",latr)
      call xios_send_field("lon_u",lonu)
      call xios_send_field("lat_u",latu)
      call xios_send_field("lon_v",lonv)
      call xios_send_field("lat_v",latv)
# else  
      call xios_send_field("x_rho",xr)
      call xios_send_field("y_rho",yr)
# endif
# ifdef CURVGRID
      call xios_send_field("angle",angler)
# endif
# ifdef MASKING
      call xios_send_field("mask_rho",rmask)
# endif

!
! Critical Depth for Drying cells
!
# ifdef WET_DRY
      call xios_send_field("Dcrit",Dcrit)
# endif

      work2d=zeta(:,:,fast_indx_out)
      call xios_send_field("zeta",work2d)
      call xios_send_field("ubar",ubar(:,:,fast_indx_out))
      call xios_send_field("vbar",vbar(:,:,fast_indx_out))

# ifdef MOVING_BATHY
      call xios_send_field("hmorph",h)
# endif
# ifdef WET_DRY
      call xios_send_field("rmask_wet",rmask_wet)
      call xios_send_field("umask_wet",umask_wet)
      call xios_send_field("vmask_wet",vmask_wet)
#endif

      do j=0,Mm
         do i=0,Lm
            work2d(i,j)=0.5*sqrt((bustr(i,j)+bustr(i+1,j))**2
     &                           +(bvstr(i,j)+bvstr(i,j+1))**2)
     &                                                    *rho0
         enddo
      enddo
      call xios_send_field("bustr",bustr*rho0)
      call xios_send_field("bvstr",bvstr*rho0)
      call xios_send_field("bostr",work2d)

      do j=1,Mm
         do i=1,Lm
            work2d2(i,j)=0.5*sqrt((sustr(i,j)+sustr(i+1,j))**2
     &                           +(svstr(i,j)+svstr(i,j+1))**2)
     &                                                    *rho0
         enddo
      enddo
      call xios_send_field("wstr",work2d)

      work2d=sustr*rho0
      call xios_send_field("sustr",work2d)
      work2d=svstr*rho0 
      call xios_send_field("svstr",work2d)
#ifdef BULK_SM_UPDATE
      call xios_send_field("uwnd",uwnd)
      call xios_send_field("vwnd",vwnd)
#endif

# ifdef WKB_WWAVE
      call xios_send_field("hrm",hrm(:,:,wstp))
      call xios_send_field("frq",frq(:,:,wstp))
      call xios_send_field("wac",wac(:,:,wstp))
      call xios_send_field("wkx",wkx(:,:,wstp))
      call xios_send_field("wke",wke(:,:,wstp))
      call xios_send_field("epb",wdsp)
      call xios_send_field("epd",wdrg)
#  ifdef WAVE_ROLLER
      call xios_send_field("war",war(:,:,wstp))
      call xios_send_field("epr",rdsp)
#  endif
# endif

# ifdef MRL_WCI
      call xios_send_field("sup",sup)
      call xios_send_field("ust2d",ust2d)
      call xios_send_field("vst2d",vst2d)
#  ifdef SOLVE3D
      call xios_send_field("ust",ust)
      call xios_send_field("vst",vst)
      call xios_send_field("wst",wst)
      call xios_send_field("Akb",Akb)
      call xios_send_field("Akw",Akw)
      call xios_send_field("kvf",kvf)
      call xios_send_field("calP",calP)
      call xios_send_field("Kapsrf",Kapsrf)
#  endif
# endif

! ------ 3D variables
! nstp vs nnew : present choice follows set_avg.F
#  ifdef SOLVE3D
      call xios_send_field("u",u(:,:,:,nstp))
      call xios_send_field("v",v(:,:,:,nstp)) 
      call xios_send_field("u_surf",u(:,:,N,nstp)) 
      call xios_send_field("v_surf",v(:,:,N,nstp)) 
# ifdef NBQ
      call xios_send_field("temp",t(:,:,:,nstp,itemp)) 
#   else
      call xios_send_field("temp",t(:,:,:,nstp,itemp)) 
      call xios_send_field("temp_surf",t(:,:,N,nstp,itemp)) 
# endif
# ifdef SALINITY
      call xios_send_field("salt",t(:,:,:,nstp,isalt))
      call xios_send_field("salt_surf",t(:,:,N,nstp,isalt))
# endif /* SALINITY */
! WARNING here passive tracers are missing
! we could easily add them in using all the vname stuff
! but I found it quiet ugly

      call xios_send_field("rho",rho)

# if defined ANA_VMIX || defined BVF_MIXING \
  || defined LMD_MIXING || defined LMD_SKPP || defined LMD_BKPP \
  || defined GLS_MIXING || defined UV_VIS_SMAGO_3D
      call xios_send_field("bvf",bvf)
# endif

      do k=0,N
         do j=0,Mm+1
            do i=0,Lm+1
               work(i,j,k)= ( We(i,j,k)
#  ifdef VADV_ADAPT_IMP
     &                    +   Wi(i,j,k)
#  endif                  
     &                           ) *pm(i,j)*pn(i,j) 
#  ifdef NBQ_MASS
     &                         /rho_nbq_avg1(i,j,k)
#  endif
            enddo
         enddo
      enddo
      call xios_send_field("omega",work)

#  ifdef NBQ
      call xios_send_field("w_nbq",wz(:,:,:,nstp))
      call xios_send_field("rho_nbq",rho_nbq)
#  endif

      call Wvlcty (tile, workr) 
      call xios_send_field("w",workr)

#  ifdef DIAGNOSTICS_EDDY

         GRho=-1*g/rho0

          do k=1,N
            do j=JstrR,JendR
              do i=IstrR,IendR

                  GRho=-1*g/rho0

                  workru = 0.5*(u(i+1,j,k,nstp)+u(i,j,k,nstp))
                  workrv = 0.5*(v(i,j+1,k,nstp)+v(i,j,k,nstp))

                enddo
              enddo
           enddo

          call xios_send_field("uu", workru**2)
          call xios_send_field("uv", workru*workrv)
          call xios_send_field("vv", workrv**2)
          call xios_send_field("ub", workru*GRho*rho)
          call xios_send_field("vb", workrv*GRho*rho)
          call xios_send_field("wb", workr *GRho*rho)
          call xios_send_field("uT", workru*t(:,:,:,nstp,itemp))
          call xios_send_field("vT", workrv*t(:,:,:,nstp,itemp))
          call xios_send_field("wT", workr *t(:,:,:,nstp,itemp))
          call xios_send_field("uS", workru*t(:,:,:,nstp,isalt))
          call xios_send_field("vS", workrv*t(:,:,:,nstp,isalt))
          call xios_send_field("wS", workr *t(:,:,:,nstp,isalt))

#  endif


        
#  ifdef VIS_COEF_3D
      call xios_send_field("visc3d",visc3d_r)
#  endif

#  ifdef DIF_COEF_3D
      do k=1,N
         do j=1,Mm
            do i=1,Lm
               workr(i,j,k)=
#   ifdef TS_DIF2
     &              diff2(i,j,itemp)
#    ifdef TS_DIF_SMAGO
     &              +diff3d_r(i,j,k)
#    endif
#   elif defined TS_DIF4
     &              diff4(i,j,itemp)
#    ifdef TS_DIF_SMAGO
     &              +diff3d_r(i,j,k)*om_r(i,j)*on_r(i,j)
#    endif
     &              +0.25*(diff3d_u(i,j,k)+diff3d_u(i+1,j,k)
     &              +diff3d_v(i,j,k)+diff3d_v(i,j+1,k))
#   endif
            enddo
         enddo
      enddo
      call xios_send_field("diff3d",workr)
#  endif

      call xios_send_field("AKv",Akv)
      work=Akt(:,:,:,itemp)
      call xios_send_field("AKt",work)
#  ifdef SALINITY      
      work=Akt(:,:,:,isalt)
      call xios_send_field("AKs",work)
#  endif
#  ifdef GLS_MIXING
!      call xios_send_field("AKk",Akk)
!      call xios_send_field("AKp",Akp)
#  endif
#  ifdef LMD_SKPP
#   ifdef LMD_SKPP2005
!      work2d=hbls(:,:,nstp)
#   else
!      work2d=hbl
#   endif
!      call xios_send_field("hbl",work2d)
#  elif defined GLS_MIXING
!      work2d=hbl
!      call xios_send_field("hbl",work2d)
#  endif
#  ifdef LMD_BKPP
      call xios_send_field("hbbl",work2d)
#  endif
#  ifdef GLS_MIXING
!      call xios_send_field("tke",tke(:,:,:,nstp))
!      call xios_send_field("gls",gls(:,:,:,nstp))
      call xios_send_field("Lscale",Lscale)
#  endif

      work2d=stflx(:,:,itemp)*rho0*Cp
     &       SWITCH rmask
      call xios_send_field("shflx",work2d)
#  ifdef SALINITY    
         do j=0,Mm+1
            do i=0,Lm+1
               work2d(i,j)=stf_cff*stflx(i,j,isalt)/
     &              ( max(eps,t(i,j,N,nstp,isalt))) 
     &              SWITCH rmask(i,j)
            enddo
        enddo
      call xios_send_field("swflx",work2d)
#  endif

#  ifdef BULK_FLUX
      work2d=shflx_rsw*rho0*Cp
      call xios_send_field("radsw",work2d)
#  else 
      work2d=srflx*rho0*Cp
      call xios_send_field("swrad",work2d)
#  endif
      
#  ifdef BULK_FLUX
      work2d=shflx_rlw*rho0*Cp
      call xios_send_field("shflx_rlw",work2d)

      work2d=shflx_lat*rho0*Cp
      call xios_send_field("shflx_lat",work2d)

      work2d=shflx_sen*rho0*Cp
      call xios_send_field("shflx_sen",work2d)
#  endif
#  ifdef SST_SKIN
      call xios_send_field("sst_skin",sst_skin)
#  endif
     

#  ifdef DIAGNOSTICS_UV

      itrc = 1
          call xios_send_field("u_rate", Mrate(:,:,:,itrc))
          call xios_send_field("u_adv", MXadv(:,:,:,itrc)
     &           + MYadv(:,:,:,itrc) + MVadv(:,:,:,itrc))
          call xios_send_field("u_Cor",MCor(:,:,:,itrc))
          call xios_send_field("u_Prsgrd",MPrsgrd(:,:,:,itrc))
          call xios_send_field("u_Hmix",MHmix(:,:,:,itrc,nstp))
          call xios_send_field("u_Hdiff",MHdiff(:,:,:,itrc))
          call xios_send_field("u_Vmix",MVmix(:,:,:,itrc))
          call xios_send_field("u_Vmix2",MVmix2(:,:,:,itrc))
#  ifdef DIAGNOSTICS_BARO
          call xios_send_field("u_Baro",MBaro(:,:,:,itrc))
#  endif
#  ifdef M3FAST
          call xios_send_field("u_fast",Mfast(:,:,:,itrc))
#  endif
      itrc = 2
          call xios_send_field("v_rate", Mrate(:,:,:,itrc))
          call xios_send_field("v_adv", MXadv(:,:,:,itrc)
     &         + MYadv(:,:,:,itrc) + MVadv(:,:,:,itrc))
          call xios_send_field("v_Cor",MCor(:,:,:,itrc))
          call xios_send_field("v_Prsgrd",MPrsgrd(:,:,:,itrc))
          call xios_send_field("v_Hmix",MHmix(:,:,:,itrc,nstp))
          call xios_send_field("v_Hdiff",MHdiff(:,:,:,itrc))
          call xios_send_field("v_Vmix",MVmix(:,:,:,itrc))
          call xios_send_field("v_Vmix2",MVmix2(:,:,:,itrc))
#  ifdef DIAGNOSTICS_BARO
          call xios_send_field("v_Baro",MBaro(:,:,:,itrc))
#  endif
#  ifdef M3FAST
          call xios_send_field("v_fast",Mfast(:,:,:,itrc))
#  endif
#  endif


#  ifdef DIAGNOSTICS_VRT
          call xios_send_field("vrtrate", vrtrate(:,:))
          call xios_send_field("vrtadv", vrtXadv(:,:) + vrtYadv(:,:))
          call xios_send_field("vrtCor",vrtCor(:,:))
          call xios_send_field("vrtPrsgrd",vrtPrsgrd(:,:))
          call xios_send_field("vrtHmix",vrtHmix(:,:))
          call xios_send_field("vrtHdiff",vrtHdiff(:,:))
!          call xios_send_field("vrtVmix",vrtVmix(:,:))
          call xios_send_field("vrtVmix2",vrtVmix2(:,:))
          call xios_send_field("vrtBaro",vrtBaro(:,:))
          call xios_send_field("vrtDrag",vrtDrag(:,:))
          call xios_send_field("vrtWind",vrtWind(:,:))
#  ifdef M3FAST
          call xios_send_field("vrtfast",vrtfast(:,:))
#  endif
#  endif


#  ifdef DIAGNOSTICS_EK

          call xios_send_field("ekrate", ekrate(:,:))
          call xios_send_field("ekadv", ekHadv(:,:) + ekVadv(:,:))
          call xios_send_field("ekCor",ekCor(:,:))
          call xios_send_field("ekPrsgrd",ekPrsgrd(:,:))
          call xios_send_field("ekHmix",ekHmix(:,:))
          call xios_send_field("ekHdiff",ekHdiff(:,:))
          call xios_send_field("ekVmix",ekVmix(:,:))
          call xios_send_field("ekVmix2",ekVmix2(:,:))
          call xios_send_field("ekBaro",ekBaro(:,:))
          call xios_send_field("ekvol",ekvol(:,:))
          call xios_send_field("ekDrag",ekDrag(:,:))
          call xios_send_field("ekWind",ekWind(:,:))
#  ifdef M3FAST
          call xios_send_field("ekfast",ekfast(:,:))
#  endif
#  endif






#  if defined BIOLOGY
#   if defined PISCES
      CALL trc_wri_pisces
#   else
      call xios_send_field("hel",hel)     
#   if (defined BIO_NChlPZD ||  defined BIO_N2ChlPZD2)
      call xios_send_field("theta",theta)
#    ifdef OXYGEN
      call xios_send_field("U10",u10)
      call xios_send_field("KvO2",Kv_O2)
      call xios_send_field("O2sat",O2satu)
#    endif /* OXYGEN */
#   elif defined BIO_BioEBUS 
      call xios_send_field("AOU",AOU)
      call xios_send_field("wind10",wind10)
#   endif
#   endif
#  endif /* BIOLOGY */


#  ifdef SEDIMENT
      do i=1,NGRAV
         write(nametrc,'(A,I1)') 'gravel_',i
         call xios_send_field(TRIM(nametrc),
     &    t(:,:,:,nstp,1+ntrc_salt+ntrc_pas+ntrc_bio+i))
      enddo      
      do i=1,NSAND
         write(nametrc,'(A,I1)') 'sand_',i
         call xios_send_field(TRIM(nametrc),
     &    t(:,:,:,nstp,1+ntrc_salt+ntrc_pas+ntrc_bio+ngrav+i)) 
      enddo      
      do i=1,NMUD
         write(nametrc,'(A,I1)') 'mud_',i
         call xios_send_field(TRIM(nametrc),
     &    t(:,:,:,nstp,1+ntrc_salt+ntrc_pas+ntrc_bio+ngrav+nsand+i)) 
      enddo      

      call xios_send_field("bed_thick",bed_thick) 
      call xios_send_field("bed_poros",bed_poros) 

      do i=1,NGRAV
         write(nametrc,'(A,I1)') 'bed_frac_grav_',i
         call xios_send_field(TRIM(nametrc),
     &    bed_frac(:,:,:,i))
         write(nametrc,'(A,I1)') 'bed_mass_grav_',i
         call xios_send_field(TRIM(nametrc),
     &    bed_mass(:,:,:,nstp,i))
      enddo      
      do i=1,NSAND
         write(nametrc,'(A,I1)') 'bed_frac_sand_',i
         call xios_send_field(TRIM(nametrc),
     &    bed_frac(:,:,:,NGRAV+i))
         write(nametrc,'(A,I1)') 'bed_mass_sand_',i
         call xios_send_field(TRIM(nametrc),
     &    bed_mass(:,:,:,nstp,NGRAV+i))
      enddo      
      do i=1,NMUD
         write(nametrc,'(A,I1)') 'bed_frac_mud_',i
         call xios_send_field(TRIM(nametrc),
     &    bed_frac(:,:,:,NGRAV+NSAND+i))
         write(nametrc,'(A,I1)') 'bed_mass_mud_',i
         call xios_send_field(TRIM(nametrc),
     &    bed_mass(:,:,:,nstp,NGRAV+NSAND+i))
      enddo      

#   ifdef SUSPLOAD
      do i=1,NGRAV
         write(nametrc,'(A,I1)') 'dflux_grav_',i
         call xios_send_field(TRIM(nametrc),
     &    settling_flux(:,:,i))
         write(nametrc,'(A,I1)') 'eflux_grav_',i
         call xios_send_field(TRIM(nametrc),
     &    ero_flux(:,:,i))
      enddo      
      do i=1,NSAND
         write(nametrc,'(A,I1)') 'dflux_sand_',i
         call xios_send_field(TRIM(nametrc),
     &    settling_flux(:,:,NGRAV+i))
         write(nametrc,'(A,I1)') 'eflux_sand_',i
         call xios_send_field(TRIM(nametrc),
     &    ero_flux(:,:,NGRAV+i))
      enddo      
      do i=1,NMUD
         write(nametrc,'(A,I1)') 'dflux_mud_',i
         call xios_send_field(TRIM(nametrc),
     &    settling_flux(:,:,NGRAV+NSAND+i))
         write(nametrc,'(A,I1)') 'edflux_mud_',i
         call xios_send_field(TRIM(nametrc),
     &    ero_flux(:,:,NGRAV+NSAND+i))
      enddo      
#   endif

#   ifdef BEDLOAD
      do i=1,NGRAV
         write(nametrc,'(A,I1)') 'bdlu_gravel_',i
         call xios_send_field(TRIM(nametrc),
     &    bedldu(:,:,i))
         write(nametrc,'(A,I1)') 'bdlv_gravel_',i
         call xios_send_field(TRIM(nametrc),
     &    bedldv(:,:,i))
      enddo      
      do i=1,NSAND
         write(nametrc,'(A,I1)') 'bdlu_sand_',i
         call xios_send_field(TRIM(nametrc),
     &    bedldu(:,:,NGRAV+i))
         write(nametrc,'(A,I1)') 'bdlv_sand_',i
         call xios_send_field(TRIM(nametrc),
     &    bedldv(:,:,NGRAV+i))
      enddo      
      do i=1,NMUD
         write(nametrc,'(A,I1)') 'bdlu_mud_',i
         call xios_send_field(TRIM(nametrc),
     &    bedldu(:,:,NGRAV+NSAND+i))
         write(nametrc,'(A,I1)') 'bdlv_mud_',i
         call xios_send_field(TRIM(nametrc),
     &    bedldv(:,:,NGRAV+NSAND+i))
      enddo      
#   endif
#  endif /* SEDIMENT */
     
#  ifdef BBL
         call xios_send_field("Abed",Abed) 
         call xios_send_field("Hripple",Hripple)
         call xios_send_field("Lripple",Lripple) 
         call xios_send_field("Zbnot",Zbnot)
         call xios_send_field("Zbapp",Zbapp)

         work2d=sqrt(bustrw**2
     &                             +bvstrw**2)
     &                                   *rho0
         call xios_send_field("bostrw",work2d)
         call xios_send_field("bustrc",bustrc*rho0)
         call xios_send_field("bvstrc",bvstrc*rho0)
         call xios_send_field("bustrw",bustrw*rho0)
         call xios_send_field("bvstrw",bvstrw*rho0)
         call xios_send_field("bustrcwmax",bustrcwmax*rho0)
         call xios_send_field("bvstrcwmax",bvstrcwmax*rho0)
!     
#  endif /* BBL */
# endif /* SOLVE3D */
  

      return
      end

#else   /* XIOS */
 
      subroutine send_xios_empty
      return
      end

#endif /* XIOS */
