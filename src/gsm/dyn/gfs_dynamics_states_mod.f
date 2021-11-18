!#include "../../../ESMFVersionDefine.h"

       MODULE gfs_dynamics_states_mod

! 
!  June 2005 	 weiyu yang  initial code.
!  february 2007 hann-ming henry juang  
!		for gfs dynamics and gaussian grid DATA.
!  March 2009           Weiyu Yang, modified for the ensemble NEMS run.
!  2009/10/05           Sarah Lu, grid_gr unfolded to 3D
!  Jun Wang    2009-11-09 set difital filter variables to export state
!  Sarah Lu    2010-02-09 set cpi, ri attributes to export state
!  Sarah Lu    2010-02-11 add set_dynexp_attribute subprogram for ri/cpi
!  Henry Juang 2011-01-03 change routine name back to gfs_dynamics_ and for NDSL
!  Weiyu Yang  2011-02    Updated to use both the ESMF 4.0.0rp2 library,
!                         ESMF 5 series library and the the
!                         ESMF 3.1.0rp2 library.
!  Weiyu Yang  2011-09    Updated to use the ESMF 5.2.0r library,
!  Philip Pegion 2014-05  Added Stochastic parameterization capability
!  S   Moorthi 2015-11    Update for digital filter with semi-Lagrangian
!                         gridded tracers
!  S   Moorthi 2015-12    Removed digital filtering from p, pd and dpdt
!                         and added rqtk when semilag=.true.
!
!!USEs:
!
      USE ESMF

! the derived TYPE of the internal state.
!----------------------------------------
      USE gfs_dynamics_internal_state_mod
      USE gfs_dynamics_namelist_mod         ! Dynamics configuration

! routines which can be used to add a fortran array to 
! an esmf state and to get a fortran array from an esmf state
!-----------------------------------------------------
      USE gfs_dynamics_grid_create_mod
      USE gfs_dynamics_add_get_state_mod
      USE gfs_dynamics_err_msg_mod
      USE gfs_dyn_tracer_config, ONLY : gfs_dyn_tracer
      use layout_grid_tracers ,  only : rgt_h,xhalo,yhalo
      use gfs_dyn_resol_def,     only : adiabatic
      use namelist_dynamics_def, only : semilag, gg_tracers
      use gfs_dyn_layout1,       only : ipt_lats_node_a,me,lats_node_a

      IMPLICIT none

      TYPE(GFS_Dyn_State_Namelist) :: cf

! !REVISION HISTORY:
!  da Silva  2009-01-22  First version
!  Sarah Lu  2009-01-26  Revised to include all 3d atmos fields
!  Sarah Lu  2009-08-06  Add gfs_dynamics_import2internal_mgrid and
!                        gfs_dynamics_internal2export_mgrid
!  Sarah Lu  2009-10-12  Port to the latest trunk
!  Sarah Lu  2009-10-17  Tracer bundle added; (shum, oz, cld) removed
!  S Moorthi 2015-11-10  Adding semilag+gg_tracer option for gridded tracers
!
      CONTAINS

      SUBROUTINE gfs_dynamics_import2internal(imp_gfs_dyn, int_state, rc, exp_gfs_dyn)

! this subroutine can be used to update the initial condition 
! fields of the internal state from the esmf inport state.
!------------------------------------------------------------

! every possible import state has its own turn-on/turn-off switch flag
! which can be used to fit different interface requirement from different
! outside grid component systems.
!------------------------------------------------------------------------

      TYPE(esmf_state)                                          :: imp_gfs_dyn  
      TYPE(gfs_dynamics_internal_state), POINTER, INTENT(inout) :: int_state 
      INTEGER,          OPTIONAL,                 INTENT(out)   :: rc     
      TYPE(ESMF_State), OPTIONAL,                 INTENT(in)    :: exp_gfs_dyn

      TYPE(ESMF_Field)                  :: Field
      TYPE(ESMF_FieldBundle)            :: Bundle

      INTEGER                           :: rc1, rcfinal, ylan, lat, lons_lat,    &
                                           i, lan, k, n, kk, kstr, kend, mstr, mend

      REAL, DIMENSION(:, :),    POINTER :: FArr2D
      REAL, DIMENSION(:, :, :), POINTER :: FArr3D

! initialize the error signal variables.
!---------------------------------------
      rc1     = esmf_success
      rcfinal = esmf_success

      CALL esmf_logwrite("update internal state with esmf import state", &
                          ESMF_LOGMSG_INFO, rc = rc1)

      cf = int_state%esmf_sta_list


! check the internal state, if it is ready to run in the start_step
! then we don't need to do import to internal state
!-----------------------------------------------------------------------
      IF( int_state%start_step ) THEN
        PRINT *,' It is starting, so no need for import_state2internal '
        RETURN
!     ELSE
!        PRINT *,' do import state to internal state '
      END IF

! get the surface orography array from the esmf import state.
!------------------------------------------------------------
      IF(cf%z_import == 1) THEN
        kstr = int_state%g_gz
        IF(ASSOCIATED(FArr2D)) NULLIFY(FArr2D)
        IF(.NOT. PRESENT(exp_gfs_dyn)) THEN
          CALL getf90arrayfromstate(imp_gfs_dyn, 'hs',     FArr2D, 0, rc = rc1)
        ELSE
          CALL getf90arrayfromstate(exp_gfs_dyn, 'hs_dfi', FArr2D, 0, rc = rc1)
        END IF

        CALL gfs_dynamics_err_msg(rc1,"gete esmf state - hs_im",rcfinal)

        int_state%grid_gr(:,:,kstr) = FArr2D
      END IF

! get the surface pressure array from the esmf import state.
! for the detailed comments for every computational steps
! please refer to the surface orography array code.
!-----------------------------------------------------------
      IF(cf%ps_import == 1) THEN
        IF(ASSOCIATED(FArr2D)) NULLIFY(FArr2D)
        IF(.NOT. PRESENT(exp_gfs_dyn)) THEN
          IF(int_state%ENS .AND. int_state%Cpl_flag) THEN
            mstr = int_state%g_q
            CALL GetF90ArrayFromState(imp_gfs_dyn, 'pps', FArr2D, 0, rc = rc1)
            int_state%grid_gr(:,:,mstr) = FArr2D
            mstr = int_state%g_qm
            CALL GetF90ArrayFromState(imp_gfs_dyn, 'psm', FArr2D, 0, rc = rc1)
            int_state%grid_gr(:,:,mstr) = FArr2D
          ELSE
            kstr = int_state%g_zqp
            CALL getf90arrayfromstate(imp_gfs_dyn, 'ps', FArr2D, 0, rc = rc1)
            int_state%grid_gr(:,:,kstr) = FArr2D
!           enddo
          END IF
        ELSE
          kstr = int_state%g_zqp
          CALL getf90arrayfromstate(exp_gfs_dyn, 'ps_dfi', FArr2D, 0, rc = rc1)
          int_state%grid_gr(:,:,kstr) = FArr2D
        END IF

        CALL gfs_dynamics_err_msg(rc1,"gete esmf state - ps_im",rcfinal)
      END IF

!     if (gg_tracers) then
!     IF(cf%rqtk_import == 1) THEN

      if (semilag .and. gg_tracers .and. .not. adiabatic) then
        IF(ASSOCIATED(FArr2D)) NULLIFY(FArr2D)
!       IF(.NOT. PRESENT(exp_gfs_dyn)) THEN
!         IF(int_state%ENS .AND. int_state%Cpl_flag) THEN
!           mstr = int_state%g_q
!           CALL GetF90ArrayFromState(imp_gfs_dyn, 'pps', FArr2D, 0, rc = rc1)
!           int_state%grid_gr(:,:,mstr) = FArr2D
!           mstr = int_state%g_qm
!           CALL GetF90ArrayFromState(imp_gfs_dyn, 'psm', FArr2D, 0, rc = rc1)
!           int_state%grid_gr(:,:,mstr) = FArr2D
!         ELSE
          if(.not. int_state%adiabatic) then
            kstr = int_state%g_rqtk
            CALL getf90arrayfromstate(imp_gfs_dyn, 'rqtk', FArr2D, 0, rc = rc1)
            int_state%grid_gr(:,:,kstr) = FArr2D
          endif
!         END IF
!       ELSE
!         kstr = int_state%g_rqtk
!         CALL getf90arrayfromstate(exp_gfs_dyn, 'rqtk_dfi', FArr2D, 0, rc = rc1)
!         int_state%grid_gr(:,:,kstr) = FArr2D
!       END IF

        CALL gfs_dynamics_err_msg(rc1,"gete esmf state - rqtk_im",rcfinal)
      endif

! get the temperature array from the esmf import state.
!------------------------------------------------------
      IF(cf%temp_import == 1) THEN
        IF(ASSOCIATED(FArr3D)) NULLIFY(FArr3D)
        IF(.NOT. PRESENT(exp_gfs_dyn)) THEN
          IF(int_state%ENS .AND. int_state%Cpl_flag) THEN
            mstr = int_state%g_tt
            mend = mstr + int_state%levs - 1
            CALL GetF90ArrayFromState(imp_gfs_dyn, 'tt', FArr3D, 0, rc = rc1)
            int_state%grid_gr(:,:,mstr:mend) = FArr3D
            mstr = int_state%g_ttm
            mend = mstr + int_state%levs - 1
            CALL GetF90ArrayFromState(imp_gfs_dyn, 'tm', FArr3D, 0, rc = rc1)
            int_state%grid_gr(:,:,mstr:mend) = FArr3D
          ELSE
            kstr = int_state%g_ttp
            kend = kstr + int_state%levs - 1
            CALL getf90arrayfromstate(imp_gfs_dyn, 't', FArr3D, 0, rc = rc1)
            int_state%grid_gr(:,:,kstr:kend) = FArr3D
          END IF
        ELSE
          kstr = int_state%g_ttp
          kend = kstr + int_state%levs - 1
          CALL getf90arrayfromstate(exp_gfs_dyn, 't_dfi', FArr3D, 0, rc = rc1)
          int_state%grid_gr(:,:,kstr:kend) = FArr3D
        END IF

        CALL gfs_dynamics_err_msg(rc1,"gete esmf state - t_im",rcfinal)
      END IF

! get the zonal-wind array from the esmf import state.
! for detailed line by line comments please refer to 
! the temperature code.
!-----------------------------------------------------
      IF(cf%u_import == 1) THEN
        IF(ASSOCIATED(FArr3D)) NULLIFY(FArr3D)
        IF(.NOT. PRESENT(exp_gfs_dyn)) THEN
          IF(int_state%ENS .AND. int_state%Cpl_flag) THEN
            mstr = int_state%g_uu
            mend = mstr + int_state%levs - 1
            CALL GetF90ArrayFromState(imp_gfs_dyn, 'uu', FArr3D, 0, rc = rc1)
            int_state%grid_gr(:,:,mstr:mend) = FArr3D
            mstr = int_state%g_uum
            mend = mstr + int_state%levs - 1
            CALL GetF90ArrayFromState(imp_gfs_dyn, 'um', FArr3D, 0, rc = rc1)
            int_state%grid_gr(:,:,mstr:mend) = FArr3D
          ELSE
            CALL getf90arrayfromstate(imp_gfs_dyn, 'u', FArr3D, 0, rc = rc1)
            kstr = int_state%g_uup
            kend = kstr + int_state%levs - 1
            int_state%grid_gr(:,:,kstr:kend) = FArr3D
          END IF
        ELSE
          CALL getf90arrayfromstate(exp_gfs_dyn, 'u_dfi', FArr3D, 0, rc = rc1)
          kstr = int_state%g_uup
          kend = kstr + int_state%levs - 1
          int_state%grid_gr(:,:,kstr:kend) = FArr3D
        END IF

        CALL gfs_dynamics_err_msg(rc1,"gete esmf state - u_im",rcfinal)
      END IF

! get the meridian-wind array from the esmf import state.
!-----------------------------------------------------
      IF(cf%v_import == 1) THEN
        IF(ASSOCIATED(FArr3D)) NULLIFY(FArr3D)
        IF(.NOT. PRESENT(exp_gfs_dyn)) THEN
          IF(int_state%ENS .AND. int_state%Cpl_flag) THEN
            mstr = int_state%g_vv
            mend = mstr + int_state%levs - 1
            CALL GetF90ArrayFromState(imp_gfs_dyn, 'vv', FArr3D, 0, rc = rc1)
            int_state%grid_gr(:,:,mstr:mend) = FArr3D
            mstr = int_state%g_vvm
            mend = mstr + int_state%levs - 1
            CALL GetF90ArrayFromState(imp_gfs_dyn, 'vm', FArr3D, 0, rc = rc1)
            int_state%grid_gr(:,:,mstr:mend) = FArr3D
          ELSE
            CALL getf90arrayfromstate(imp_gfs_dyn, 'v', FArr3D, 0, rc = rc1)
            kstr = int_state%g_vvp
            kend = kstr + int_state%levs - 1
            int_state%grid_gr(:,:,kstr:kend) = FArr3D
          END IF
        ELSE
          CALL getf90arrayfromstate(exp_gfs_dyn, 'v_dfi', FArr3D, 0, rc = rc1)
          kstr = int_state%g_vvp
          kend = kstr + int_state%levs - 1
          int_state%grid_gr(:,:,kstr:kend) = FArr3D
        END IF

        CALL gfs_dynamics_err_msg(rc1,"gete esmf state - v_im",rcfinal)
      END IF

! get the tracer array from the esmf import state.
!-------------------------------------------------
      IF(cf%tracer_import == 1) THEN
        IF(.NOT. PRESENT(exp_gfs_dyn)) THEN
          IF(int_state%ENS .AND. int_state%Cpl_flag) THEN
            CALL ESMF_StateGet(imp_gfs_dyn, 'tracers_ENS', Bundle, rc = rc1 )
            CALL gfs_dynamics_err_msg(rc1, 'retrieve Ebundle from state', rcfinal)
            DO k = 1, int_state%ntrac
              mstr = int_state%g_rq + (k - 1) * int_state%levs
              mend = mstr + int_state%levs - 1
              IF(ASSOCIATED(FArr3D)) NULLIFY(FArr3D)
              CALL ESMF_FieldBundleGet(Bundle,                               &
!                                      trim(int_state%gfs_dyn_tracer%vname(k, 2)), &
                                       trim(gfs_dyn_tracer%vname(k, 2)),     &
                                       field = Field,                        &
                                       rc = rc1)
              CALL gfs_dynamics_err_msg(rc1, 'retrieve Efield from bundle', rcfinal)

              CALL ESMF_FieldGet(Field, farrayPtr = FArr3D, localDE = 0, rc = rc1)

              CALL gfs_dynamics_err_msg(rc1, 'retrieve Farray from field', rcfinal)
              int_state%grid_gr(:, :, mstr:mend) = FArr3D
            END DO
            DO k = 1, int_state%ntrac
              mstr = int_state%g_rm + (k - 1) * int_state%levs
              mend = mstr + int_state%levs - 1
              NULLIFY(FArr3D)
              CALL ESMF_FieldBundleGet(Bundle,                               &
!                                      trim(int_state%gfs_dyn_tracer%vname(k, 3)), &
                                       trim(gfs_dyn_tracer%vname(k, 3)),     &
                                       field = Field,                        &
                                       rc = rc1)
              CALL gfs_dynamics_err_msg(rc1, 'retrieve Efield from bundle', rcfinal)

              CALL ESMF_FieldGet(Field, farrayPtr = FArr3D, localDE = 0, rc = rc1)

              CALL gfs_dynamics_err_msg(rc1, 'retrieve Farray from field', rcfinal)
              int_state%grid_gr(:, :, mstr:mend) = FArr3D
            END DO
          ELSE
            CALL ESMF_StateGet(imp_gfs_dyn, 'tracers', Bundle, rc = rc1 )
            CALL gfs_dynamics_err_msg(rc1, 'retrieve Ebundle from state', rcfinal)
            DO k = 1, int_state%ntrac
              kstr = int_state%g_rqp + (k - 1) * int_state%levs
              kend = kstr + int_state%levs - 1
              IF(ASSOCIATED(FArr3D)) NULLIFY(FArr3D)
              CALL ESMF_FieldBundleGet(Bundle,                               &
!                                      trim(int_state%gfs_dyn_tracer%vname(k, 1)), &
                                       trim(gfs_dyn_tracer%vname(k, 1)),     &
                                       field = Field,                        &
                                       rc = rc1)
              CALL gfs_dynamics_err_msg(rc1, 'retrieve Efield from bundle', rcfinal)
                  
              CALL ESMF_FieldGet(Field, farrayPtr = FArr3D, localDE = 0, rc = rc1)

              CALL gfs_dynamics_err_msg(rc1, 'retrieve Farray from field', rcfinal)
              int_state%grid_gr(:, :, kstr:kend) = FArr3D
            END DO
          END IF
        ELSE
!         if (semilag .and. gg_tracers .and. .not. adiabatic) then
!           DO n = 1, int_state%ntrac
!             CALL getf90arrayfromstate(exp_gfs_dyn,                  &
!             TRIM(gfs_dyn_tracer%vname(n, 1))//'_dfi', FArr3D, 0, rc=rc1)
!             do k=1,int_state%levs
!               kk = n*int_state%levs - k + 1
!               kk = int_state%levs - k + 1
!!$omp parallel do private(i,lan,ylan,lat,lons_lat)
!               do lan=1,int_state%lats_node_a
!                 ylan = yhalo + lan
!                 ylan = int_state%lats_node_a + 1 - lan + yhalo
!                 lat  = int_state%global_lats_a(ipt_lats_node_a+int_state%lats_node_a-lan)
!                 lons_lat = int_state%lonsperlat(lat)
!                 do i=1,lons_lat
!                   rgt_h(xhalo+i,k,ylan,n) = FArr3D(i,lan,kk)
!                 enddo
!               enddo
!             enddo
!             CALL gfs_dynamics_err_msg(rc1, 'retrieve Farray from field', rcfinal)
!           END DO
!         else
            DO k = 1, int_state%ntrac
              CALL getf90arrayfromstate(exp_gfs_dyn,                  &
!             TRIM(int_state%gfs_dyn_tracer%vname(k, 1))//'_dfi',     &
              TRIM(gfs_dyn_tracer%vname(k, 1))//'_dfi', FArr3D, 0, rc = rc1)
              kstr = int_state%g_rqp + (k - 1) * int_state%levs
              kend = kstr + int_state%levs - 1
              int_state%grid_gr(:, :, kstr:kend) = FArr3D
              CALL gfs_dynamics_err_msg(rc1, 'retrieve Farray from field', rcfinal)

            END DO
!         endif
        END IF
      END IF

! get the pressure array from the esmf import state.
!-------------------------------------------------------------
      IF(cf%p_import == 1) THEN
        IF(.NOT. PRESENT(exp_gfs_dyn)) THEN
          CALL getf90arrayfromstate(imp_gfs_dyn, 'p',     FArr3D, 0, rc = rc1)
        ELSE
!       ELSEif (.not. semilag) then
          CALL getf90arrayfromstate(exp_gfs_dyn, 'p_dfi', FArr3D, 0, rc = rc1)
        END IF
        kstr = int_state%g_p
        kend = kstr + int_state%levs - 1
        int_state%grid_gr(:, :, kstr:kend) = FArr3D

        CALL gfs_dynamics_err_msg(rc1,"gete esmf state - p_im",rcfinal)
      END IF

! get the pressure layer depth (dp) array from the esmf import state.
!-------------------------------------------------------------
      IF(cf%dp_import == 1) THEN
        IF(.NOT. PRESENT(exp_gfs_dyn)) THEN
          CALL getf90arrayfromstate(imp_gfs_dyn, 'dp',     FArr3D, 0, rc = rc1)
        ELSE
!       ELSEif (.not. semilag) then
          CALL getf90arrayfromstate(exp_gfs_dyn, 'dp_dfi', FArr3D, 0, rc = rc1)
        END IF
        kstr = int_state%g_dpp
        kend = kstr + int_state%levs - 1
        int_state%grid_gr(:, :, kstr:kend) = FArr3D
        CALL gfs_dynamics_err_msg(rc1,"gete esmf state - dp_im",rcfinal)
      END IF

!
! get the omega (dpdt) array from the esmf import state.
!-------------------------------------------------------------

      IF(cf%dpdt_import == 1) THEN
        IF(.NOT. PRESENT(exp_gfs_dyn)) THEN
          CALL getf90arrayfromstate(imp_gfs_dyn, 'dpdt',     FArr3D, 0, rc = rc1)
        ELSE
!       ELSEif (.not. semilag) then
          CALL getf90arrayfromstate(exp_gfs_dyn, 'dpdt_dfi', FArr3D, 0, rc = rc1)
        END IF
        kstr = int_state%g_dpdt
        kend = kstr + int_state%levs - 1
        int_state%grid_gr(:, :, kstr:kend) = FArr3D

        CALL gfs_dynamics_err_msg(rc1,"gete esmf state - dpdt_im",rcfinal)
      END IF

!
!
! PRINT out the final error signal message and put it to rc.
!-----------------------------------------------------------
      IF(PRESENT(rc)) CALL gfs_dynamics_err_msg_final(rcfinal,          &
                                     "gfs_dynamics_import2internal",rc)

      END SUBROUTINE gfs_dynamics_import2internal


! =========================================================================


      SUBROUTINE gfs_dynamics_internal2export(int_state, exp_gfs_dyn, rc)

! 
! this subroutine will be changed that all export
! esmf states will get data directly from the gfs internal structure data arrays
! 

!
!!uses:
!
!
! !input/output variables and parameters:
!----------------------------------------
      TYPE(esmf_state),                           INTENT(inout) :: exp_gfs_dyn
      TYPE(gfs_dynamics_internal_state), POINTER, INTENT(inout) :: int_state
      INTEGER, OPTIONAL,                          INTENT(out)   :: rc

! local array size parameter of the esmf export state arrays.
!------------------------------------------------------------
      INTEGER               :: rc1     ! error signal work variable.
      INTEGER               :: rcfinal ! the final error signal variable.
      INTEGER               :: k, kstr, kend, mstr, mend

      TYPE(ESMF_Field)                              :: Field
      TYPE(ESMF_FieldBundle), SAVE                  :: Bundle

      REAL(KIND = kind_evod), DIMENSION(:, :),    POINTER ::    &
                   hold_z,    hold_ps,     hold_pps,  hold_psm, &
                   hold_pps6, hold_psm6
!                  hold_pps6, hold_psm6,   hold_rqtk

      REAL(KIND = kind_evod), DIMENSION(:, :, :), POINTER ::    &
                   hold_temp,   hold_u,    hold_v,              &
                   hold_p,      hold_dp,   hold_dpdt,           &
                   hold_ttemp,  hold_uu,   hold_vv,             &
                   hold_tempm,  hold_um,   hold_vm,             &
                   hold_ttemp6, hold_uu6,  hold_vv6,            &
                   hold_tempm6, hold_um6,  hold_vm6,            &
                   hold_shum_wts, hold_sppt_wts, hold_skebu_wts,&
                   hold_skebv_wts,hold_vcu_wts,hold_vcv_wts

      REAL, DIMENSION(:, :, :), POINTER :: FArr3D

      LOGICAL, SAVE :: first
      DATA first/.true./
      SAVE hold_z,    hold_ps,     hold_temp, hold_u, hold_v, &
           hold_p,    hold_dp,     hold_dpdt,                 &
           hold_psm,  hold_tempm,  hold_um,   hold_vm,        &
           hold_pps,  hold_ttemp,  hold_uu,   hold_vv,        &
           hold_psm6, hold_tempm6, hold_um6,  hold_vm6,       &
           hold_pps6, hold_ttemp6, hold_uu6,  hold_vv6,       &
           hold_shum_wts, hold_sppt_wts, hold_skebu_wts,      &
           hold_skebv_wts,hold_vcu_wts,hold_vcv_wts

! initialize the error signal variables.
!---------------------------------------
      rc1     = esmf_success
      rcfinal = esmf_success

      IF(.NOT. first) RETURN

      cf = int_state%esmf_sta_list

! PRINT out the information.
!-----------------------------------------------

      CALL esmf_logwrite("begining to put the esmf export state.", &
                          ESMF_LOGMSG_INFO, rc = rc1)

! orography field. gaussian grid
!----------------------------------------------------------------
      IF(cf%z_export == 1) THEN
        kstr   = int_state%g_gz
        hold_z => int_state%grid_gr(:, :, kstr)
        CALL addf90arraytostate(exp_gfs_dyn, mgrid, 'hs', hold_z, rc = rc1)
        CALL gfs_dynamics_err_msg(rc1, "done z_export.", rcfinal)
        IF(int_state%ndfi > 0 .AND. cf%z_import == 1) THEN
          hold_z => int_state%grid_gr_dfi%hs(:, :, 1)
          CALL addf90arraytostate(exp_gfs_dyn, mgrid, 'hs_dfi', hold_z, rc = rc1)
          CALL gfs_dynamics_err_msg(rc1,"add to esmf export state -z-dfi",rcfinal)
        END IF
      END IF

! surface pressure
!----------------------------------------------------------------

      IF(cf%ps_export == 1) THEN
        IF(int_state%ENS) THEN
          mstr      =  int_state%g_q
          hold_pps  => int_state%grid_gr (:,  :, mstr)
          hold_pps6 => int_state%grid_gr6(:, :, mstr - 1)
          mstr      =  int_state%g_qm
          hold_psm  => int_state%grid_gr (:,  :, mstr)
          hold_psm6 => int_state%grid_gr6(:, :, mstr - 1)
          CALL addf90arraytostate(exp_gfs_dyn, mgrid, 'pps',  hold_pps,  rc = rc1)
          CALL addf90arraytostate(exp_gfs_dyn, mgrid, 'psm',  hold_psm,  rc = rc1)
          CALL addf90arraytostate(exp_gfs_dyn, mgrid, 'pps6', hold_pps6, rc = rc1)
          CALL addf90arraytostate(exp_gfs_dyn, mgrid, 'psm6', hold_psm6, rc = rc1)
        END IF
        kstr    = int_state%g_zqp
        hold_ps => int_state%grid_gr(:, :, kstr)
        CALL addf90arraytostate(exp_gfs_dyn, mgrid, 'ps', hold_ps, rc = rc1)
        CALL gfs_dynamics_err_msg(rc1, "put to esmf state - ps_ex", rcfinal)
        IF(int_state%ndfi > 0 .AND. cf%ps_import == 1) THEN
          hold_ps => int_state%grid_gr_dfi%ps(:, :, 1)
          CALL addf90arraytostate(exp_gfs_dyn, mgrid, 'ps_dfi', hold_ps, rc = rc1)
          CALL gfs_dynamics_err_msg(rc1,"add to esmf export state -ps-dfi",rcfinal)
        END IF
      END IF
!     if (gg_tracers) then
!       kstr    = int_state%g_rqtk
!       hold_rqtk => int_state%grid_gr(:, :, kstr)
!       CALL addf90arraytostate(exp_gfs_dyn, mgrid, 'rqtk',hold_rqtk, rc = rc1)
!       CALL gfs_dynamics_err_msg(rc1, "put to esmf state - rqtk_ex", rcfinal)
!       IF(int_state%ndfi > 0 ) THEN
!         hold_rqtk => int_state%grid_gr_dfi%rqtk(:, :)
!         CALL addf90arraytostate(exp_gfs_dyn, mgrid, 'rqtk_dfi', hold_rqtk, rc = rc1)
!         CALL gfs_dynamics_err_msg(rc1,"add to esmf export state -rqtk-dfi",rcfinal)
!       endif
!     endif

! add the temperature fileds into the esmf export state.
!-------------------------------------------------------

      IF(cf%temp_export == 1) THEN
        IF(int_state%ENS) THEN
          mstr        = int_state%g_tt
          mend        = mstr + int_state%levs - 1
          hold_ttemp  => int_state%grid_gr (:, :, mstr     : mend)
          hold_ttemp6 => int_state%grid_gr6(:, :, mstr - 1 : mend - 1)
          mstr        = int_state%g_ttm
          mend        = mstr + int_state%levs - 1
          hold_tempm  => int_state%grid_gr (:, :, mstr     : mend)
          hold_tempm6 => int_state%grid_gr6(:, :, mstr - 1 : mend - 1)

          CALL addf90arraytostate(exp_gfs_dyn, mgrid, 'tt',  hold_ttemp,  rc = rc1)
          CALL addf90arraytostate(exp_gfs_dyn, mgrid, 'tm',  hold_tempm,  rc = rc1)
          CALL addf90arraytostate(exp_gfs_dyn, mgrid, 'tt6', hold_ttemp6, rc = rc1)
          CALL addf90arraytostate(exp_gfs_dyn, mgrid, 'tm6', hold_tempm6, rc = rc1)
        END IF
        kstr = int_state%g_ttp
        kend = kstr + int_state%levs - 1
        hold_temp => int_state%grid_gr(:, :, kstr : kend)
        CALL addf90arraytostate(exp_gfs_dyn, mgrid, 't', hold_temp, rc = rc1)
        CALL gfs_dynamics_err_msg(rc1, "put to esmf state - t_ex", rcfinal)
        IF(int_state%ndfi > 0 .AND. cf%temp_import == 1) THEN
          hold_temp => int_state%grid_gr_dfi%t(:, :, :)
          CALL addf90arraytostate(exp_gfs_dyn, mgrid, 't_dfi', hold_temp, rc = rc1)
          CALL gfs_dynamics_err_msg(rc1,"add to esmf export state -t-dfi",rcfinal)
        END IF
      END IF

! to add the u field into the esmf export state.  for the detailed
! description comments please refer to the temperature field.
!--------------------------------------------------------------------------

      IF(cf%u_export == 1) THEN
        IF(int_state%ENS) THEN
          mstr        = int_state%g_uu
          mend        = mstr + int_state%levs - 1
          hold_uu     => int_state%grid_gr (:, :, mstr     : mend)
          hold_uu6    => int_state%grid_gr6(:, :, mstr - 1 : mend - 1)
          mstr        = int_state%g_uum
          mend        = mstr + int_state%levs - 1
          hold_um     => int_state%grid_gr (:, :, mstr     : mend)
          hold_um6    => int_state%grid_gr6(:, :, mstr - 1 : mend - 1)

          CALL addf90arraytostate(exp_gfs_dyn, mgrid, 'uu',  hold_uu,  rc = rc1)
          CALL addf90arraytostate(exp_gfs_dyn, mgrid, 'um',  hold_um,  rc = rc1)
          CALL addf90arraytostate(exp_gfs_dyn, mgrid, 'uu6', hold_uu6, rc = rc1)
          CALL addf90arraytostate(exp_gfs_dyn, mgrid, 'um6', hold_um6, rc = rc1)
        END IF
        kstr   = int_state%g_uup
        kend   = kstr + int_state%levs - 1
        hold_u => int_state%grid_gr(:, :, kstr : kend)
        CALL addf90arraytostate(exp_gfs_dyn, mgrid, 'u', hold_u, rc = rc1)
        CALL gfs_dynamics_err_msg(rc1, "put to esmf state - u_ex", rcfinal)
        IF(int_state%ndfi > 0 .AND. cf%u_import == 1) THEN
          hold_u => int_state%grid_gr_dfi%u(:, :, :)
          CALL addf90arraytostate(exp_gfs_dyn, mgrid, 'u_dfi', hold_u, rc = rc1)
          CALL gfs_dynamics_err_msg(rc1,"add to esmf export state -u-dfi",rcfinal)
        END IF
      END IF

! add the v field into the esmf export state.
!----------------------------------------------------

      IF(cf%v_export == 1) THEN
        IF(int_state%ENS) THEN
          mstr        = int_state%g_vv
          mend        = mstr + int_state%levs - 1
          hold_vv     => int_state%grid_gr (:, :, mstr     : mend)
          hold_vv6    => int_state%grid_gr6(:, :, mstr - 1 : mend - 1)
          mstr        = int_state%g_vvm
          mend        = mstr + int_state%levs - 1
          hold_vm     => int_state%grid_gr (:, :, mstr     : mend)
          hold_vm6    => int_state%grid_gr6(:, :, mstr - 1 : mend - 1)

          CALL addf90arraytostate(exp_gfs_dyn, mgrid, 'vv',  hold_vv,  rc = rc1)
          CALL addf90arraytostate(exp_gfs_dyn, mgrid, 'vm',  hold_vm,  rc = rc1)
          CALL addf90arraytostate(exp_gfs_dyn, mgrid, 'vv6', hold_vv6, rc = rc1)
          CALL addf90arraytostate(exp_gfs_dyn, mgrid, 'vm6', hold_vm6, rc = rc1)
        END IF
        kstr   = int_state%g_vvp
        kend   = kstr + int_state%levs - 1
        hold_v => int_state%grid_gr(:, :, kstr : kend)
        CALL addf90arraytostate(exp_gfs_dyn, mgrid, 'v', hold_v, rc = rc1)
        CALL gfs_dynamics_err_msg(rc1, "put to esmf state - v_ex", rcfinal)
        IF(int_state%ndfi > 0 .AND. cf%v_import == 1) THEN
          hold_v => int_state%grid_gr_dfi%v(:, :, :)
          CALL addf90arraytostate(exp_gfs_dyn, mgrid, 'v_dfi', hold_v, rc = rc1)
          CALL gfs_dynamics_err_msg(rc1,"add to esmf export state -v_dfi",rcfinal)
        END IF
      END IF

! add the tracer fields into the esmf export state.
!---------------------------------------------------
      IF(cf%tracer_export == 1) THEN
        Bundle = ESMF_FieldBundleCreate(name = 'tracers', rc = rc1)
        CALL gfs_dynamics_err_msg(rc1, "create empty fieldbundle", rcfinal)

        IF(int_state%ENS) THEN
          DO k = 1, int_state%ntrac
            mstr = int_state%g_rq + (k - 1) * int_state%levs
            mend = mstr + int_state%levs - 1
            NULLIFY(FArr3D)
            FArr3D => int_state%grid_gr(:, :, mstr : mend)
            Field  = ESMF_FieldCreate(mgrid, FArr3D,                    &
                                      name=trim(gfs_dyn_tracer%vname(k,2)), rc=rc1)
!           name   = trim(int_state%gfs_dyn_tracer%vname(k, 2)), rc = rc1)
            CALL ESMF_FieldBundleAdd(Bundle, (/Field/), rc = rc1)
            CALL gfs_dynamics_err_msg(rc1, "add field to bundle, vname=2", rcfinal)

            NULLIFY(FArr3D)
            FArr3D => int_state%grid_gr6(:, :, mstr - 1 : mend - 1)
            Field  = ESMF_FieldCreate(mgrid, FArr3D,                    &
                                      name=trim(gfs_dyn_tracer%vname(k,4)), rc=rc1)
!           name   = trim(int_state%gfs_dyn_tracer%vname(k, 4)), rc = rc1)
            CALL ESMF_FieldBundleAdd(Bundle, (/Field/), rc = rc1)
            CALL gfs_dynamics_err_msg(rc1, "add field to bundle, vname=4", rcfinal)
          END DO
          DO k = 1, int_state%ntrac
            mstr = int_state%g_rm + (k - 1) * int_state%levs
            mend = mstr + int_state%levs - 1
            NULLIFY(FArr3D)
            FArr3D => int_state%grid_gr(:, :, mstr : mend)
            Field  = ESMF_FieldCreate(mgrid, FArr3D,                    &
                                      name=trim(gfs_dyn_tracer%vname(k,3)), rc=rc1)
!           name = trim(int_state%gfs_dyn_tracer%vname(k, 3)), rc = rc1)
            CALL ESMF_FieldBundleAdd(Bundle, (/Field/), rc = rc1)
            CALL gfs_dynamics_err_msg(rc1, "add field to bundle, vname=3", rcfinal)

            NULLIFY(FArr3D)
            FArr3D => int_state%grid_gr6(:, :, mstr - 1 : mend - 1)
            Field  = ESMF_FieldCreate(mgrid, FArr3D,                    &
                                      name=trim(gfs_dyn_tracer%vname(k,5)), rc=rc1)
!           name   = trim(int_state%gfs_dyn_tracer%vname(k, 5)), rc = rc1)
            CALL ESMF_FieldBundleAdd(Bundle, (/Field/), rc = rc1)
            CALL gfs_dynamics_err_msg(rc1, "add field to bundle, vname=5", rcfinal)
          END DO
        END IF

        DO k = 1, int_state%ntrac
          kstr = int_state%g_rqp + (k - 1) * int_state%levs
          kend = kstr + int_state%levs - 1
          NULLIFY(FArr3D)
          FArr3D => int_state%grid_gr(:, :, kstr : kend)
          Field  = ESMF_FieldCreate(mgrid, FArr3D,                    &
                                    name=trim(gfs_dyn_tracer%vname(k,1)), rc=rc1)
!         name = trim(int_state%gfs_dyn_tracer%vname(k, 1)), rc = rc1)
          CALL ESMF_FieldBundleAdd(Bundle, (/Field/), rc = rc1)
          CALL gfs_dynamics_err_msg(rc1, "add field to bundle, vname=1", rcfinal)
        END DO

!  --- set attributes for ri and cpi to be exported to physics
        call set_dynexp_attribute
!  ---  inputs:  (in scope variables)
!  ---  outputs: (in scope variables)

        CALL ESMF_StateAddReplace(exp_gfs_dyn, (/Bundle/), rc = rc1)
        CALL gfs_dynamics_err_msg(rc1, "add to esmf state - tracer", rcfinal)

        IF(int_state%ndfi > 0 .AND. cf%tracer_import == 1) THEN
          DO k = 1, int_state%ntrac
            kstr = (k - 1) * int_state%levs + 1
            kend = kstr + int_state%levs - 1
            NULLIFY(FArr3D)
            FArr3D => int_state%grid_gr_dfi%tracer(:, :, kstr : kend)
            Field  = ESMF_FieldCreate(mgrid, FArr3D,                                &
                                      name=trim(gfs_dyn_tracer%vname(k,1))//'_dfi', &
                                      rc = rc1)
!           name   = trim(int_state%gfs_dyn_tracer%vname(k,1))//'_dfi', rc = rc1)
            CALL ESMF_StateAddReplace(exp_gfs_dyn, (/Field/), rc = rc1)
            CALL gfs_dynamics_err_msg(rc1, "add field to state, vname=1", rcfinal)
          END DO
        END IF
      END IF
! add layer pressure (pp) fileds into the esmf export state.
!-------------------------------------------------------

      IF(cf%p_export == 1) THEN
          kstr = int_state%g_p
          kend = kstr + int_state%levs - 1
          hold_p => int_state%grid_gr(:, :, kstr : kend)
          CALL addf90arraytostate(exp_gfs_dyn, mgrid, 'p',       & 
                                  hold_p, rc = rc1)
          CALL gfs_dynamics_err_msg(rc1, "put to esmf state - p_ex",rcfinal)
!         IF(int_state%ndfi > 0 .AND. cf%p_import == 1 .and. .not.  semilag) THEN
          IF(int_state%ndfi > 0 .AND. cf%p_import == 1) THEN
              hold_p => int_state%grid_gr_dfi%p(:, :, :)
              CALL addf90arraytostate(exp_gfs_dyn, mgrid, 'p_dfi',             & 
                                      hold_p, rc = rc1)
              CALL gfs_dynamics_err_msg(rc1,"add to esmf export state -p-dfi",rcfinal)
          END IF
      END IF

! add pressure depth (dp) fileds into the esmf export state.
!-------------------------------------------------------

      IF(cf%dp_export == 1) THEN
          kstr = int_state%g_dpp
          kend = kstr + int_state%levs - 1
          hold_dp => int_state%grid_gr(:, :, kstr : kend)
          CALL addf90arraytostate(exp_gfs_dyn, mgrid, 'dp',          & 
                                  hold_dp, rc = rc1)
          CALL gfs_dynamics_err_msg(rc1, "put to esmf state - dp_ex",rcfinal)
!         IF(int_state%ndfi > 0 .AND. cf%dp_import == 1 .and. .not.  semilag) THEN
          IF(int_state%ndfi > 0 .AND. cf%dp_import == 1) THEN
              hold_dp => int_state%grid_gr_dfi%dp(:, :, :)
              CALL addf90arraytostate(exp_gfs_dyn, mgrid, 'dp_dfi',             & 
                                      hold_dp, rc = rc1)
              CALL gfs_dynamics_err_msg(rc1,"add to esmf export state -dp-dfi",rcfinal)
          END IF
      END IF


! add omega (dpdt) fileds into the esmf export state.
!-------------------------------------------------------

      IF(cf%dpdt_export == 1) THEN
          kstr = int_state%g_dpdt
          kend = kstr + int_state%levs - 1
          hold_dpdt => int_state%grid_gr(:, :, kstr : kend)
          CALL addf90arraytostate(exp_gfs_dyn, mgrid, 'dpdt',        & 
                                  hold_dpdt, rc = rc1)
          CALL gfs_dynamics_err_msg(rc1, "put to esmf state - dpdt_ex",rcfinal)
!         IF(int_state%ndfi > 0 .AND. cf%dpdt_import == 1 .and. .not.  semilag) THEN
          IF(int_state%ndfi > 0 .AND. cf%dpdt_import == 1) THEN
              hold_dpdt => int_state%grid_gr_dfi%dpdt(:, :, :)
              CALL addf90arraytostate(exp_gfs_dyn, mgrid, 'dpdt_dfi',             & 
                                      hold_dpdt, rc = rc1)
              CALL gfs_dynamics_err_msg(rc1,"add to esmf export state -dpdt-dfi",rcfinal)
          END IF
      END IF
!
! add stochastic humidity perturbations (shum_wts) fileds into the esmf export state.
!-------------------------------------------------------
      IF(cf%shum_wts_export == 1) THEN
          hold_shum_wts => int_state%shum_wts(:, :, : )
          CALL addf90arraytostate(exp_gfs_dyn, mgrid, 'shum_wts',    & 
                                  hold_shum_wts, rc = rc1)
          CALL gfs_dynamics_err_msg(rc1,                             &
                  "put to esmf state - shum_wts_ex",rcfinal)
      END IF
!
!
! add SPPT perturbations (sppt_wts) fileds into the esmf export state.
!-------------------------------------------------------

      IF(cf%sppt_wts_export == 1) THEN
          hold_sppt_wts => int_state%sppt_wts(:, :, : )
          CALL addf90arraytostate(exp_gfs_dyn, mgrid, 'sppt_wts',    & 
                                  hold_sppt_wts, rc = rc1)
          CALL gfs_dynamics_err_msg(rc1,                             &
                  "put to esmf state - sppt_wts_ex",rcfinal)
      END IF
!
!
! add SKEB perturbations (skebu_wts and skebv_wts) fileds into the esmf export state.
!-------------------------------------------------------

      IF(cf%skeb_wts_export == 1) THEN
          hold_skebu_wts => int_state%skebu_wts(:, :, : )
          hold_skebv_wts => int_state%skebv_wts(:, :, : )
          CALL addf90arraytostate(exp_gfs_dyn, mgrid, 'skebu_wts',    & 
                                  hold_skebu_wts, rc = rc1)
          CALL gfs_dynamics_err_msg(rc1,                             &
                  "put to esmf state - skebu_wts_ex",rcfinal)
          CALL addf90arraytostate(exp_gfs_dyn, mgrid, 'skebv_wts',    & 
                                  hold_skebv_wts, rc = rc1)
          CALL gfs_dynamics_err_msg(rc1,                             &
                  "put to esmf state - skebv_wts_ex",rcfinal)
      END IF
!
! add SKEB perturbations (vcu_wts and vcv_wts) fileds into the esmf export state.
!-------------------------------------------------------

      IF(cf%vc_wts_export == 1) THEN
          hold_vcu_wts => int_state%vcu_wts(:, :, : )
          hold_vcv_wts => int_state%vcv_wts(:, :, : )
          CALL addf90arraytostate(exp_gfs_dyn, mgrid, 'vcu_wts',    & 
                                  hold_vcu_wts, rc = rc1)
          CALL gfs_dynamics_err_msg(rc1,                             &
                  "put to esmf state - vcu_wts_ex",rcfinal)
          CALL addf90arraytostate(exp_gfs_dyn, mgrid, 'vcv_wts',    & 
                                  hold_vcv_wts, rc = rc1)
          CALL gfs_dynamics_err_msg(rc1,                             &
                  "put to esmf state - vcv_wts_ex",rcfinal)
      END IF
!
!! print out the final error signal information and put it to the rc.
!!-------------------------------------------------------------------
      IF(PRESENT(rc)) CALL gfs_dynamics_err_msg_final(rcfinal,       &
                                   "gfs_dynamics_internal2export", rc)

      first = .false.

      contains
! =================


!-----------------------------
      subroutine set_dynexp_attribute
!.............................
!  ---  inputs:  (in scope variables)
!  ---  outputs: (in scope variables)

       TYPE(ESMF_Info) :: info

! ---  Set attributes for ri and cpi

       CALL ESMF_InfoGetFromHost(Bundle, info, rc=RC)
       call gfs_dynamics_err_msg(rc,"retrieve info from dyn_exp",rcfinal)

       CALL ESMF_InfoSet(info                                   &  !<-- Info handle from dyn export state tracer bundle
               ,key   ='cpi_dryair'                             &  !<-- Name of the attribute to insert
               ,value = gfs_dyn_tracer%cpi(0)                   &  !<-- Value of the attribute
               ,rc    = RC)
       call gfs_dynamics_err_msg(rc,"insert cpi(0) attribute to dyn_exp",rcfinal)

       CALL ESMF_InfoSet(info                                   &  !<-- Info handle from dyn export state tracer bundle
               ,key   ='ri_dryair'                              &  !<-- Name of the attribute to insert
               ,value = gfs_dyn_tracer%ri(0)                    &  !<-- Value of the attribute
               ,rc    = RC)
       call gfs_dynamics_err_msg(rc,"insert ri(0) attribute to dyn_exp",rcfinal)

       CALL ESMF_InfoSet(info                                   &  !<-- Info handle from dyn export state tracer bundle
               ,key    ='cpi'                                   &  !<-- Name of the attribute array
               ,values = gfs_dyn_tracer%cpi(1:int_state%ntrac)  &  !<-- The array being inserted
               ,rc     = RC)
       call gfs_dynamics_err_msg(rc,"insert cpi(:) attribute to dyn_exp",rcfinal)

       CALL ESMF_InfoSet(info                                   &  !<-- Info handle from dyn export state tracer bundle
               ,key    ='ri'                                    &  !<-- Name of the attribute array
               ,values = gfs_dyn_tracer%ri(1:int_state%ntrac)   &!<-- The array being inserted
               ,rc     = RC)
       call gfs_dynamics_err_msg(rc,"insert ri(:) attribute to dyn_exp",rcfinal)

      RETURN
      end subroutine set_dynexp_attribute

      END SUBROUTINE gfs_dynamics_internal2export


! ==========================================================================

      END MODULE gfs_dynamics_states_mod
