!=====================================================================
!
!               S p e c f e m 3 D  V e r s i o n  2 . 1
!               ---------------------------------------
!
!          Main authors: Dimitri Komatitsch and Jeroen Tromp
!    Princeton University, USA and CNRS / INRIA / University of Pau
! (c) Princeton University / California Institute of Technology and CNRS / INRIA / University of Pau
!                             July 2012
!
! This program is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 2 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License along
! with this program; if not, write to the Free Software Foundation, Inc.,
! 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
!
!=====================================================================

!--------------------------------------------------------------------------------------------------
!
! Interface to USGS Bay Area 3-D seismic velocity model.
! http://earthquake.usgs.gov/data/3dgeologic/
!
! Brad Aagaard, USGS
!
!--------------------------------------------------------------------------------------------------


  module usgscencalvm_model

    ! pointers (must match sizeof(void*) in C)
    integer(kind=8) :: query
    integer(kind=8) :: errHandler

    ! Hardwired parameters
    character(len=64), parameter :: filenameLog = "cencalvm.log"
    character(len=64), parameter :: filenameDB = "USGSBayAreaVM-08.3.0.etree"
    character(len=64), parameter :: filenameDBExt = "USGSBayAreaVMExt-08.3.0.etree"
    
    integer, parameter :: queryType = 2 ! waveres
    double precision, parameter :: queryRes = 1.0
    double precision, parameter :: minVs = 500.0
    integer, parameter :: cacheSize = 512
!---
!
! ADD YOUR MODEL HERE
!
!---

  ! only here to illustrate an example
  !  type model_external_variables
  !    sequence
  !    double precision dvs(0:dummy_size)
  !  end type model_external_variables
  !  type (model_external_variables) MEXT_V

  end module usgscencalvm_model

!
!-------------------------------------------------------------------------------------------------
!

  subroutine model_usgscencalvm_open(myrank)

    ! setup model

    use usgscencalvm_model
    use generate_databases_par, only: TOMOGRAPHY_PATH

    implicit none

    include "constants.h"

    ! args
    integer :: myrank

    ! local variables
    integer :: ok
    character(len=256) :: filenameLog_full
    character(len=256) :: filenameDB_full
    character(len=256) :: filenamdDBExt_full

    ! Create query
    call cencalvm_createquery_f(query)
    if(query == 0) stop "USGS CenCalVM module: Could not create query."

    ! Get handle to error handler
    call cencalvm_errorhandler_f(query, errHandler)
    if(errHandler == 0) stop "USGS CenCalVM module: Could not get error handler."

    ! Set log file in error handler
    filenameLog_full = TOMOGRAPHY_PATH(1:len_trim(TOMOGRAPHY_PATH))//trim(filenameLog)
    call cencalvm_error_logfilename_f(errHandler, filenameLog_full, ok)
    call model_usgscencalvm_checkerr(ok)

    ! Set database filename
    filenameDB_full = TOMOGRAPHY_PATH(1:len_trim(TOMOGRAPHY_PATH))//trim(filenameDB)
    call cencalvm_filename_f(query, filenameDB_full, ok)
    call model_usgscencalvm_checkerr(ok)

    ! Set extended database filename
    !filenameDBExt_full = TOMOGRAPHY_PATH(1:len_trim(TOMOGRAPHY_PATH))//trim(filenameDBExt)
    !call cencalvm_filename_f(query, filenameDBExt_full, ok)
    !call model_usgscencalvm_checkerr(ok)

    call cencalvm_cachesize_f(query, cacheSize, ok)
    call model_usgscencalvm_checkerr(ok)

    call cencalvm_open_f(query, ok)
    call model_usgscencalvm_checkerr(ok)

    ! Set query type
    call cencalvm_querytype_f(query, queryType, ok)
    call model_usgscencalvm_checkerr(ok)

    ! Set query resolution
    call cencalvm_queryres_f(query, queryRes, ok)
    call model_usgscencalvm_checkerr(ok)

  end subroutine model_usgscencalvm_open

!
!-------------------------------------------------------------------------------------------------
!

  subroutine model_usgscencalvm_close(myrank)

    ! routine to close model

    use usgscencalvm_model
    
    implicit none
    
    include "constants.h"

    ! args
    integer :: myrank

    ! local variables
    integer :: ok

    call cencalvm_close_f(query, ok)
    call model_usgscencalvm_checkerr(ok)

    call cencalvm_destroyquery_f(query, ok)
    call model_usgscencalvm_checkerr(ok)

  end subroutine model_usgscencalvm_close

!
!-------------------------------------------------------------------------------------------------
!

  subroutine model_usgscencalvm_checkerr(ok)

    use usgscencalvm_model

    ! args
    integer :: ok

    ! local variables
    character(len=256) errorMsg

    call cencalvm_error_message_f(errHandler, errorMsg, ok)
    write(6,*) "USGS CenCalVM module: ERROR - "//errorMsg

    ! if query generated an error, then bail out, otherwise reset status
    if(ok == 2) stop "USGS CenCalVM module: FATAL ERROR"
    call cencalvm_error_resetstatus_f(errHandler, ok)

  end subroutine model_usgscencalvm_checkerr
  
!
!-------------------------------------------------------------------------------------------------
!


  subroutine model_usgscencalvm_values(xmesh,ymesh,zmesh,rho,vp,vs,qkappa_atten,qmu_atten,iflag_aniso,idomain_id )

    ! given a GLL point, returns super-imposed velocity model values

    use generate_databases_par,only: &
         nspec => NSPEC_AB,ibool,UTM_PROJECTION_ZONE,SUPPRESS_UTM_PROJECTION

    use create_regions_mesh_ext_par

    use usgscencalvm_model

    implicit none

    ! GLL point
    double precision, intent(in) :: xmesh,ymesh,zmesh

    ! density, Vp and Vs
    real(kind=CUSTOM_REAL) :: vp,vs,rho

    ! attenuation flag
    real(kind=CUSTOM_REAL) :: qkappa_atten,qmu_atten

    ! anisotropy flag
    integer :: iflag_aniso

    ! acoustic/elastic/.. domain flag ( 1 = acoustic / 2 = elastic / ... )
    integer :: idomain_id

    ! local variables
    double precision :: lon
    double precision :: lat
    double precision :: elev

    integer :: ok

    integer, parameter :: numVals = 9
    double precision :: vals(9)

    call utm_geo(lon,lat,xmesh,ymesh,UTM_PROJECTION_ZONE,IUTM2LONGLAT,SUPPRESS_UTM_PROJECTION)

    call cencalvm_query_f(query,vals,numVals,lon,lat,elev,ok)
    call model_usgscencalvm_checkerr(ok)

    vp = vals(1)
    vs = vals(2)
    rho = vals(3)

    ! attenuation: Qs
    qmu_atten = vals(4)

    ! attenuation Qp
    qkappa_atten = 9999.

    ! no anisotropy
    iflag_aniso = 0

    ! elastic material
    idomain_id = IDOMAIN_ELASTIC

  end subroutine model_usgscencalvm_values
