#include "../../../ESMFVersionDefine.h"

!
      module atmos_phy_chem_cpl_comp_mod

!-----------------------------------------------------------------------
!
!** This module holds the phy-to-chem coupler's register and run routines
!** -> setservices (registers init and run) is called by GOCART_SETUP
!** -> init associates tracer bundle in phy export state with iAERO bundle 
!**    in chem import state;
!** -> run transfers/converts data from phy export state to chem import state
!
!! Code Revision:
!! 11Nov 2009     Sarah Lu,   First Crack
!! 18Nov 2009     Sarah Lu,   Revise coupler run to do data copy
!! 29Dec 2009     Sarah Lu,   Comments added for clarification
!! 01Feb 2010     Sarah Lu,   Extend to include all 2d/3d fields needed by
!!                            GOCART
!! 07Feb 2010     Sarah Lu,   Add getrh subroutine
!! 11Feb 2010     Sarah Lu,   Add get_attribute subroutine to retrieve tracer
!!                            specification 
!! 12Feb 2010     Sarah Lu,   Include chemical tracers
!! 06Mar 2010     Sarah Lu,   Flip vertical profile index from bottom-up to
!!                            top-down; add Init routine
!! 12Mar 2010     Sarah Lu,   Clean up phy2chem run
!! 13Mar 2010     Sarah Lu,   Use m_chars: set statename to lower case
!! 16Mar 2010     Sarah Lu,   Add FillDefault_ and patch_; Make dimension and
!!                            tracer specification public (to chem2phy coupler)
!! 23Mar 2010     Sarah Lu,   Call rtc to track run-time-clock; debug print optional
!! 09Apr 2010     Sarah Lu,   Clean up the code
!! 09May 2010     Sarah Lu,   Revise species name for SU, OC, and BC
!! 13May 2010     Sarah Lu,   Change GetPointer_3D_ from private to public
!! 09Jun 2010     Sarah Lu,   Remove ref to m_chars/lowercase; add g_fixchar
!! 10Jun 2010     Sarah Lu,   Remove aerosol tarcer pointer for iAERO
!! 30Jun 2010     Sarah Lu,   Revise g_fixchar to allow longer char string
!! 04Aug 2010     Sarah Lu,   run_DU[SU,SS,OC,BC] are determined from chemReg
!! 09Sep 2010     Sarah Lu,   wet1 is exported from phys, no longer calculated
!!                            in the phy-to-chem coupler; correct how cn_prcp
!!                            and ncn_prcp are computed
!! 10Oct 2010     Sarah Lu,   pass g2d_fld%met from chem_imp to phys_exp
!! 15Oct 2010     Sarah Lu,   pass fscav from chem_imp to phys_exp
!! 16Oct 2010     Sarah Lu,   change g2d_fld%met from instant to accumulated 
!! 08Nov 2010     Sarah Lu,   set zle floor values to hs; aer_diag fields
!!                            are modified
!! 14Nov 2010     Sarah Lu,   pass deltim from phy_exp to chem_imp in init
!! 29Dec 2010     Sarah Lu,   Fields not used by GOCART are removed from diag
!!   Feb 2011     Weiyu Yang, Updated to use both the ESMF 4.0.0rp2 library,
!!                            ESMF 5 library and the the ESMF 3.1.0rp2 library.
!! 12May 2011     Weiyu Yang, Modified for using the ESMF 5.2.0r_beta_snapshot_07.
!! 27Nov 2011     Sarah Lu,   Modified to pass kdt from phys_exp to chem_imp;
!!                            specify i,j dimension in mapping aerosol array
!! 06Feb 2012     Weiyu Yang, Modified for using the ESMF 5.2.0r library.
!! 30Sep 2014     Sarah Lu,   Remove the code passing fscav from chem_imp to phys_exp
!! 02Oct 2014     Sarah Lu,   Add additional 2D and 3D fields as GOCART import state
!!                            has changed
!! 11Nov 2014     S Moorthi   Some cleanup and optimization and bug fix
!!                            threading is yet to be done
!! 02Dec 2014     Jun Wang,   remove O3 and update to esmf6.3
!! 17Dec 2014     Jun Wang,   fix oout of boundary problem when doing aerosol flipflop 
!! 26Dec 2014     Jun Wang,   set addtional fields for GOCART import
!!                            state in NGACv.2
!------------------------------------------------------------------------------

      USE ESMF

      USE MODULE_ERR_MSG, ONLY: ERR_MSG, MESSAGE_CHECK
      use MODULE_gfs_machine,  ONLY: kind_phys
      use MODULE_gfs_physcons, ONLY: con_rd,  con_fvirt, con_g, &
                                     con_eps, con_epsm1, con_rerth, &
                                     con_pi
      use MODULE_gfs_tropp,    ONLY: tpause
      use MODULE_gfs_funcphys

      USE Chem_RegistryMod

!-----------------------------------------------------------------------
!
      implicit none
      SAVE
!
!-----------------------------------------------------------------------
!
!     Gaussian grid 
      integer(ESMF_KIND_I4),allocatable,public :: lonsperlar_r(:)
      real(ESMF_KIND_R8),allocatable,public    :: slat_r(:), dlat_r(:)
      integer, public               :: lonr, lats_node_r, lats_node_r_max
      integer, public               :: im, jm, km

!     Tracer specification
      integer, public               :: ntrac
      logical, public               :: run_DU, run_SU, run_SS, run_OC, run_BC
      character(10), allocatable    :: spec(:)

! --- public interface
      public::  SetServices, GetPointer_tracer_, CkPointer_, GetPointer_3D_, &
                GetPointer_diag_

      private
      TYPE(Chem_Registry)          :: chemReg      !<-- The GOCART Chem_Registry
      logical, parameter           :: lckprnt = .false.

      contains
!
!-----------------------------------------------------------------------
!#######################################################################
!-----------------------------------------------------------------------
!
      subroutine setservices(GC, RC_REG)
!
!-----------------------------------------------------------------------
!!
!! This routine register the coupler component's init and run routines    
!!
!! Code Revision:
!! 11Nov 2009     Sarah Lu, First Crack
!! 06Mar 2010     Sarah Lu, Register init routine
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!***  argument variables
!-----------------------------------------------------------------------
!
      type(ESMF_cplcomp)               :: gc         ! coupler component
!
      integer,intent(out) :: rc_reg                  ! return code for register
!
!-----------------------------------------------------------------------
!***  local variables
!-----------------------------------------------------------------------
!
      integer :: rc=ESMF_success                     ! the error signal variable

!-----------------------------------------------------------------------
!***  register the coupler component's init routine
!-----------------------------------------------------------------------
!
      MESSAGE_CHECK="Set Entry Point for phy2chem coupler init"

      call ESMF_CplCompSetEntryPoint(GC                        & !<-- The gridded component
                                    ,ESMF_METHOD_INITIALIZE    & !<-- Predefined subroutine type
                                    ,INIT                      & !<-- User's subroutineName
                                    ,rc=RC)

      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_REG)

!-----------------------------------------------------------------------
!***  register the coupler component's run routine
!-----------------------------------------------------------------------
!
      MESSAGE_CHECK="Set Entry Point for phy2chem coupler run"

      call ESMF_CplCompSetEntryPoint(GC                        & !<-- The gridded component
                                    ,ESMF_METHOD_RUN           & !<-- Predefined subroutine type
                                    ,RUN                       & !<-- User's subroutineName
                                    ,rc=RC)

      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_REG)

!-----------------------------------------------------------------------
!***  Check the final error signal variable 
!-----------------------------------------------------------------------
!
      IF(RC_REG==ESMF_SUCCESS)THEN
        WRITE(0,*)'PHY2CHEM CPL SET_SERVICES SUCCEEDED'
      ELSE
        WRITE(0,*)'PHY2CHEM CPL SET_SERVICES FAILED RC_REG=',RC_REG
      ENDIF

!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
      end subroutine setservices

!-----------------------------------------------------------------------
!#######################################################################
!-----------------------------------------------------------------------
!
      subroutine init(GC, PHY_EXP_STATE, CHEM_IMP_STATE, CLOCK, RC_CPL)
!
!-----------------------------------------------------------------------
!!
!! This routine associates tracer arrays in chem import state (iAERO bundle)
!! to phy export state (tracers bundle)
!!
!! Code Revision:
!! 06Mar 2010     Sarah Lu, First Crack
!! 04Aug 2010     Sarah Lu, Determine run_DU[SU,SS,OC,BC] from ChemRegistry
!! 14Nov 2010     Sarah Lu, Pass deltim from phy_exp to chem_imp
!-----------------------------------------------------------------------

      implicit none
!
!-----------------------------------------------------------------------
!***  argument variables
!-----------------------------------------------------------------------

      type(ESMF_cplcomp)               :: GC
      type(ESMF_state)                 :: PHY_EXP_STATE, CHEM_IMP_STATE
      type(ESMF_clock)                 :: CLOCK
!
      integer,           intent(out)   :: RC_CPL
!
!-----------------------------------------------------------------------
!***  local variables
!-----------------------------------------------------------------------

      integer                       :: rc=ESMF_success  ! the error signal variable
      integer                       :: IERR, N, L
      type(ESMF_FieldBundle)        :: Bundle, iBundle
      type(ESMF_Field)              :: Field
      character(len=ESMF_MAXSTR)    :: vname
      real                          :: deltim
!
!-----------------------------------------------------------------------
!***  Read Chem_Registry to retrive tracer name
!-----------------------------------------------------------------------

      print *, 'PHY2CHEM_INIT: get ChemReg '
      chemReg = Chem_RegistryCreate ( IERR )             !<-- read Chem_Registry

!-----------------------------------------------------------------------
!***  Pass tracer bundle from physics export to chemistry import
!-----------------------------------------------------------------------

      MESSAGE_CHECK="PHY2CHEM_INIT: get tracer bundle from phy exp"
      call ESMF_StateGet(PHY_EXP_STATE, 'tracers', iBundle, RC=RC)
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_CPL)

      MESSAGE_CHECK="PHY2CHEM_INIT: get iAERO bundle from chem imp"
      call ESMF_StateGet(CHEM_IMP_STATE, 'iAERO', Bundle, RC=RC)
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_CPL)

      allocate ( spec(chemReg%n_GOCART) )
      do L = 1, chemReg%n_GOCART

         N = chemReg%i_GOCART + L - 1

         vname = chemReg%vname(N)
         spec(L) = vname
         MESSAGE_CHECK="PHY2CHEM_INIT: get field from tracers: "//vname
         call ESMF_FieldBundleGet(iBundle, vname, FIELD=Field, rc = RC )

         CALL ERR_MSG(RC,MESSAGE_CHECK,RC_CPL)

         MESSAGE_CHECK="PHY2CHEM_INIT: add field to iAero: "//vname
         call ESMF_FieldBundleAdd(Bundle, (/Field/), rc    =RC)
         CALL ERR_MSG(RC,MESSAGE_CHECK,RC_CPL)

      end do
!
!-----------------------------------------------------------------------
!***  Determine what module(s) to be included in the simulation
!-----------------------------------------------------------------------

      run_DU = chemReg%doing_DU
      run_SU = chemReg%doing_SU
      run_OC = chemReg%doing_OC
      run_BC = chemReg%doing_BC
      run_SS = chemReg%doing_SS

!
!-----------------------------------------------------------------------
!***  Destroy Chem_Registry 
!-----------------------------------------------------------------------

      call Chem_RegistryDestroy ( chemReg, IERR )


!-----------------------------------------------------------------------
!***  Pass time-step from phy_exp to chem_imp
!-----------------------------------------------------------------------

      MESSAGE_CHECK="PHY2CHEM_INIT: get deltim from phy_exp"
      CALL ESMF_AttributeGet(PHY_EXP_STATE, name = 'deltim',  &
                             value = deltim , rc=RC)
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_CPL)

      MESSAGE_CHECK="PHY2CHEM_INIT: add deltim to chem_imp"
      CALL ESMF_AttributeSet(CHEM_IMP_STATE, name = 'deltim',  &
                             value = deltim , rc=RC)
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_CPL)

!-----------------------------------------------------------------------
!***  Check the final error signal variable
!-----------------------------------------------------------------------
!
      IF(RC_CPL==ESMF_SUCCESS)THEN
        WRITE(0,*)'PHY2CHEM CPL INIT SUCCEEDED'
      ELSE
        WRITE(0,*)'PHY2CHEM CPL INIT FAILED RC_CPL=',RC_CPL
      ENDIF
!
      end subroutine init


!-----------------------------------------------------------------------
!#######################################################################
!-----------------------------------------------------------------------
!
      subroutine run(GC, PHY_EXP_STATE, CHEM_IMP_STATE, CLOCK, RC_CPL)
!
!-----------------------------------------------------------------------
!!
!! This routine transfer/convert data from phy_exp state to chem_imp state
!!
!! Code Revision:
!! 18Nov 2009     Sarah Lu, First Crack
!! 29Dec 2009     Sarah Lu, Comments added for clarification
!! 01Feb 2010     Sarah Lu, Transfer/convert fields from phy export state 
!!                          to GOCART import state
!! 11Feb 2010     Sarah Lu, Add get_attribute
!! 06Mar 2010     Sarah Lu, Flip vertical profile index from bottom-up to
!!                          top-down
!! 12Mar 2010     Sarah Lu, Code clean up
!! 16Mar 2010     Sarah Lu, Add FillDefault_ and patch_; correct how vertical
!!                          index is flipped for 3D arrays
!! 23Mar 2010     Sarah Lu, Track run-time-clock; debug print optional
!! 09Api 2010     Sarah Lu, Remove some rtc tracking print
!! 25Aug 2012     Sarah Lu, Remove item_name
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!***  argument variables
!-----------------------------------------------------------------------

      type(ESMF_cplcomp)               :: GC
      type(ESMF_state)                 :: PHY_EXP_STATE    ! coupler import state
      type(ESMF_state)                 :: CHEM_IMP_STATE   ! coupler export state
      type(ESMF_clock)                 :: CLOCK
!
      integer,           intent(out)   :: RC_CPL

!-----------------------------------------------------------------------
!***  local variables
!-----------------------------------------------------------------------
!
      real, save                       :: deltim
!
      integer                          :: rc=ESMF_success  ! the error signal variable
      integer                          :: item_count_phys, item_count_chem
      logical, save                    :: first =  .true.
      real                             :: dx, tem2
      real(ESMF_KIND_R8), pointer      :: Array(:,:,:)
      type(ESMF_Field)                 :: Field
!
      integer                          :: kdt

!
      integer, parameter               :: nfld_2d  = 18          
!jw   integer, parameter               :: nfld_3d  = 10         
      integer, parameter               :: nfld_3d  = 9         
      integer, parameter               :: nfld_trc = 5           
      character(8)                     :: vname_2d (nfld_2d)    
      character(8)                     :: vname_3d (nfld_3d)   
      character(8)                     :: vname_trc(nfld_trc) 

      real(kind=8)                     :: rtc, t1, t1_i, t1_o, t2, t2_i, t2_o

! Fortran array for phy export state
      real(ESMF_KIND_R8), pointer, dimension(:,:) ::                     &
               p_slmsk, p_hpbl,  p_wet1,  p_stype,  p_vtype, p_vfrac,    &
               p_rain,  p_rainc, p_dtsfci,p_tsea,   p_stc1,  p_u10m,     &
               p_v10m,  p_ustar, p_zorl,  p_hs,     p_ps,    p_fice

      real (ESMF_KIND_R8), pointer, dimension(:,:,:) ::                  &
               p_t, p_u, p_v, p_p, p_dp, p_fcld, p_dqdt,                 &
               p_cnv_mfc, p_cnv_mfd, p_cnv_qc

      real (ESMF_KIND_R8), pointer, dimension(:,:,:) ::                  &
               p_spfh,                                                   & ! met tracers
!jw               p_spfh, p_o3mr,                                           & ! met tracers
               p_du001, p_du002, p_du003, p_du004, p_du005,              & ! DU
               p_ss001, p_ss002, p_ss003, p_ss004, p_ss005,              & ! SS
               p_msa,   p_so4,   p_so2,   p_dms,                         & ! SU
               p_ocphobic, p_ocphilic, p_bcphobic, p_bcphilic              ! OC/BC

! Fortran array for chem import state
      real(ESMF_KIND_R8), pointer, dimension(:,:) ::                     &
               c_lwi,   c_zpbl, c_frlake,  c_fraci,  c_wet1, c_lai,      &
               c_grn,   c_cn_prcp, c_ncn_prcp, c_sh, c_ta,   c_tsoil1,   & 
               c_u10m,  c_v10m,  c_ustar,  c_z0h,  c_tropp,  c_ps

      real (ESMF_KIND_R8), pointer, dimension(:,:,:) :: c_ple, c_zle,    &
               c_airdens, c_t, c_u, c_v, c_fcld, c_dqdt

      real (ESMF_KIND_R8), pointer, dimension(:,:,:) ::  c_o3, c_rh2       ! met tracers

! Additional Fortran arrays for chem import stat (for GOCART v2)
      real(ESMF_KIND_R8), pointer, dimension(:,:) ::          &
               c_area, c_cldtt, c_frocean, c_swndsrf,         &
               c_ts, c_u10n, c_v10n

      real (ESMF_KIND_R8), pointer, dimension(:,:,:) ::       &
!jw             c_ch4, c_h2o2, c_no3, c_oh,                   &
               c_cnv_mfc, c_cnv_mfd, c_cnv_qc, c_dqrl 
!
!---  Add the following for 2d aerosol diag fields ---
!  Fortran data pointer for phy export state
      real(ESMF_KIND_R8), pointer, dimension(:,:) ::  p_diag

!  Fortran data pointer for chem import state
      real(ESMF_KIND_R8), pointer, dimension(:,:) ::  c_diag

      TYPE(ESMF_FieldBundle)            :: Bundle, iBundle
      character*10                      :: tag, vname, BundleName, FieldName
      integer                           :: kcount, lons_lat, i, j, k, kk, n
      integer, save                     :: nfld_met
      character*10, dimension(30)       :: name_lst
      character*10,  allocatable, save  :: name_met(:)
!
! local variables for conversion
      real       :: ptp, utp, vtp, ttp, htp, shrtp, tv1, dz, tem
      real (ESMF_KIND_R8), allocatable, dimension(:), save ::           &
                              prsln, sh, rh, shs, rho, pi, h
!
      real(kind=kind_phys), parameter :: rovg = con_rd / con_g, qmin = 1.0e-10
      real(ESMF_KIND_R8),   parameter :: f_one = 1.0, f_zero = 0.0
!
      data vname_2d /'TROPP', 'LWI', 'ZPBL', 'FRLAKE',     &       
                     'FRACI', 'WET1', 'LAI', 'GRN', 'TA',  &     
                     'CN_PRCP', 'NCN_PRCP', 'PS', 'SH',    &    
                     'TSOIL1', 'U10M',  'V10M','USTAR','Z0H'/ 

      data vname_3d /'PLE', 'ZLE' , 'AIRDENS', 'FCLD', 'DQDT', &   
                     'T', 'U', 'V', 'RH2'  /               
!jw                     'T', 'U', 'V', 'O3', 'RH2'  /               

      data vname_trc/'du001', 'du002', 'du003', 'du004', 'du005'/            

!---------------------------------------------
!* Determine dimension/tracer and allocate local arrays
!---------------------------------------------
!
      IF ( FIRST ) THEN
!
!  --- Retrieve attributes (lat/lon and tracer specification) 
        MESSAGE_CHECK = 'PHY2CHEM_RUN: Retrive phy_exp attributes'
        call get_attribute(PHY_EXP_STATE, RC)
        CALL ERR_MSG(RC, MESSAGE_CHECK, RC_CPL)
!
        MESSAGE_CHECK = 'PHY2CHEM_RUN: Retrive field t from phy export'
        call ESMF_StateGet(PHY_EXP_STATE , 't', Field, rc=RC)
        CALL ERR_MSG(RC, MESSAGE_CHECK, RC_CPL)
!
        nullify(Array)
        MESSAGE_CHECK = 'PHY2CHEM_RUN: Get Fortran data pointer from t'

        CALL ESMF_FieldGet(field=Field, localDe=0, farrayPtr=Array, rc = RC)

        CALL ERR_MSG(RC, MESSAGE_CHECK, RC_CPL)

!       Determine dimension in x-, y-, and z-direction
        im = size(Array, dim=1)
        jm = size(Array, dim=2)
        km = size(Array, dim=3)

!  --- Retrieve deltim
        MESSAGE_CHECK="PHY2CHEM_RUN: get deltim from phy_exp"
        CALL ESMF_AttributeGet(PHY_EXP_STATE, name = 'deltim',  &
                               value = deltim , rc=RC)
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_CPL)


!  ---  allocate arrays at mid layers
        allocate (                               &
                     sh   (km),                  &
                     rh   (km),                  &
                     rho  (km),                  &
                     shs  (km)                   &
                    ) 
!  ---  allocate arrays at interfaces
        allocate (                               &
                     prsln(0:km),                &
                     pi   (0:km),                &
                     h    (0:km)                 &
                    )
!
!  ---  Compute all physics function tables
!
        call gfuncphys      
!
!  ---  Retrive g2d_fld%met from the bundle attribute
!
        FieldName = 'met_nfld'
        MESSAGE_CHECK="PHY2CHEM_RUN: get attribute from phy_exp"
        CALL ESMF_AttributeGet(PHY_EXP_STATE, name = FieldName,  &
                               value = kcount , rc=RC)
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_CPL)

        if ( kcount > 0 ) then

          MESSAGE_CHECK="PHY2CHEM_RUN: get bundle from phy_exp"
          BundleName='dgmet'
          CALL ESMF_StateGet(PHY_EXP_STATE, BundleName, &
                             Bundle, rc=RC)
          CALL ERR_MSG(RC,MESSAGE_CHECK,RC_CPL)

          do i = 1, kcount
            write(tag, '(i2.2)') i
            FieldName = trim(BundleName)//'_'//trim(tag)
            MESSAGE_CHECK="PHY2CHEM_RUN: get attribute from bundle"
            CALL ESMF_AttributeGet(Bundle, name=FieldName,   &
                                   value=vname, rc=RC)
            CALL ERR_MSG(RC,MESSAGE_CHECK,RC_CPL)
            name_lst(i) = trim(vname)
          enddo

          nfld_met = kcount
          allocate(name_met(kcount) )
          name_met(1:kcount) = name_lst(1:kcount)

        endif

!
! ---   Debug print (optional)
        if ( lckprnt ) then
          print *, 'PHY2CHEM_RUN: im, jm, km=', im, jm, km
          print *, 'PHY2CHEM_RUN: ntrac    =', ntrac
     	  print *, 'PHY2CHEM_RUN: lonr            =', lonr
  	  print *, 'PHY2CHEM_RUN: lats_node_r     =', lats_node_r
  	  print *, 'PHY2CHEM_RUN: lats_node_r_max =', lats_node_r_max
	  print *, 'PHY2CHEM_RUN: lonsperlar_r =', lonsperlar_r(:)
	print *, 'PHY2CHEM_RUN: slat_r       =', slat_r(:)
	print *, 'PHY2CHEM_RUN: dlat_r       =', dlat_r(:)
        endif

      ENDIF
!
!---------------------------------------------
!* Get Fortran array from phy export state
!---------------------------------------------

!  --- transfer kdt
      MESSAGE_CHECK="PHY2CHEM_RUN: Get kdt from phy export state"
      CALL ESMF_AttributeGet(PHY_EXP_STATE, name = 'kdt',  &
                             value = kdt , rc=RC)
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_CPL)
      MESSAGE_CHECK="PHY2CHEM_RUN: Add kdt to chem import state"
      CALL ESMF_AttributeSet(CHEM_IMP_STATE, name = 'kdt',  &
                             value = kdt , rc=RC)
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_CPL)


      t1_i=rtc()
      MESSAGE_CHECK="PHY2CHEM_RUN: Get ItemCount from phy export state"
      call ESMF_StateGet(PHY_EXP_STATE                    &
                        ,itemcount = item_count_phys      &
                        ,rc   =rc)
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_CPL)

      IF ( RC == ESMF_SUCCESS ) THEN

        if (item_count_phys == 0 ) then

          print *, 'PHY2CHEM_RUN: Empty phy export state; Fortran array not filled'

        else

          call GetPointer_(PHY_EXP_STATE, 'slmsk', p_slmsk, rc)
          call GetPointer_(PHY_EXP_STATE, 'hpbl',  p_hpbl , rc)
          call GetPointer_(PHY_EXP_STATE, 'wet1',  p_wet1 , rc)
          call GetPointer_(PHY_EXP_STATE, 'stype', p_stype, rc)
          call GetPointer_(PHY_EXP_STATE, 'vtype', p_vtype, rc)
          call GetPointer_(PHY_EXP_STATE, 'vfrac', p_vfrac, rc)
          call GetPointer_(PHY_EXP_STATE, 'rain',  p_rain , rc)
          call GetPointer_(PHY_EXP_STATE, 'rainc', p_rainc, rc)
          call GetPointer_(PHY_EXP_STATE, 'dtsfci',p_dtsfci,rc)
          call GetPointer_(PHY_EXP_STATE, 'tsea',  p_tsea , rc)
          call GetPointer_(PHY_EXP_STATE, 'stc1',  p_stc1 , rc)
          call GetPointer_(PHY_EXP_STATE, 'u10m',  p_u10m , rc)
          call GetPointer_(PHY_EXP_STATE, 'v10m',  p_v10m , rc)
          call GetPointer_(PHY_EXP_STATE, 'ustar', p_ustar, rc)
          call GetPointer_(PHY_EXP_STATE, 'zorl',  p_zorl , rc)
          call GetPointer_(PHY_EXP_STATE, 'hs'  ,  p_hs   , rc)
          call GetPointer_(PHY_EXP_STATE, 'ps'  ,  p_ps   , rc)
!ngac v2
          call GetPointer_(PHY_EXP_STATE, 'fice' ,  p_fice, rc)

! for GFS, vertical index is from surface to toa
          call GetPointer_3D_(PHY_EXP_STATE, 't' ,  p_t   , rc)
          call GetPointer_3D_(PHY_EXP_STATE, 'u' ,  p_u   , rc)
          call GetPointer_3D_(PHY_EXP_STATE, 'v' ,  p_v   , rc)
          call GetPointer_3D_(PHY_EXP_STATE, 'p' ,  p_p   , rc)
          call GetPointer_3D_(PHY_EXP_STATE, 'dp',  p_dp  , rc)
          call GetPointer_3D_(PHY_EXP_STATE, 'fcld',p_fcld , rc)
          call GetPointer_3D_(PHY_EXP_STATE, 'dqdt',p_dqdt , rc)
          call GetPointer_3D_(PHY_EXP_STATE, 'cnv_mfc',  p_cnv_mfc, rc)
          call GetPointer_3D_(PHY_EXP_STATE, 'cnv_mfd',  p_cnv_mfd, rc)
          call GetPointer_3D_(PHY_EXP_STATE, 'cnv_qc' ,  p_cnv_qc, rc)

! get met + chem tracers
          call GetPointer_tracer_(PHY_EXP_STATE,'spfh', p_spfh, rc)
!jw       call GetPointer_tracer_(PHY_EXP_STATE,'o3mr', p_o3mr, rc)

          if ( run_DU ) then
            call GetPointer_tracer_(PHY_EXP_STATE,'du001', p_du001, rc)
            call GetPointer_tracer_(PHY_EXP_STATE,'du002', p_du002, rc)
            call GetPointer_tracer_(PHY_EXP_STATE,'du003', p_du003, rc)
            call GetPointer_tracer_(PHY_EXP_STATE,'du004', p_du004, rc)
            call GetPointer_tracer_(PHY_EXP_STATE,'du005', p_du005, rc)
          endif 

          if ( run_SS ) then
            call GetPointer_tracer_(PHY_EXP_STATE,'ss001', p_ss001, rc)
            call GetPointer_tracer_(PHY_EXP_STATE,'ss002', p_ss002, rc)
            call GetPointer_tracer_(PHY_EXP_STATE,'ss003', p_ss003, rc)
            call GetPointer_tracer_(PHY_EXP_STATE,'ss004', p_ss004, rc)
            call GetPointer_tracer_(PHY_EXP_STATE,'ss005', p_ss005, rc)
          endif

          if ( run_SU ) then
            call GetPointer_tracer_(PHY_EXP_STATE,'MSA', p_msa, rc)
            call GetPointer_tracer_(PHY_EXP_STATE,'SO4', p_so4, rc)
            call GetPointer_tracer_(PHY_EXP_STATE,'SO2', p_so2, rc)
            call GetPointer_tracer_(PHY_EXP_STATE,'DMS', p_dms, rc)
          endif

          if ( run_OC ) then
            call GetPointer_tracer_(PHY_EXP_STATE,'OCphobic', p_ocphobic, rc)
            call GetPointer_tracer_(PHY_EXP_STATE,'OCphilic', p_ocphilic, rc)
          endif

           if ( run_BC ) then
            call GetPointer_tracer_(PHY_EXP_STATE,'BCphobic', p_bcphobic, rc)
            call GetPointer_tracer_(PHY_EXP_STATE,'BCphilic', p_bcphilic, rc)
          endif

        endif

      ELSE

         print *, 'PHY2CHEM_RUN: phy export state ItemCount failed,', RC

      ENDIF

!---------------------------------------------
!* Get Fortran array from gocart import state
!---------------------------------------------

      MESSAGE_CHECK="PHY2CHEM_RUN: Get ItemCount from chem import state"
      call ESMF_StateGet(CHEM_IMP_STATE                   &
                        ,itemcount = item_count_chem      &
                        ,rc= RC)
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_CPL)

      IF ( RC == ESMF_SUCCESS ) THEN

        if (item_count_chem == 0 ) then

          print *, 'PHY2CHEM_RUN: Empty chem import state; Fortran array not filled'

        else

          call GetPointer_(CHEM_IMP_STATE, 'LWI'   ,c_lwi  , rc)
          call GetPointer_(CHEM_IMP_STATE, 'ZPBL'  ,c_zpbl , rc)
          call GetPointer_(CHEM_IMP_STATE, 'FRLAKE',c_frlake,rc)
          call GetPointer_(CHEM_IMP_STATE, 'FRACI' ,c_fraci, rc)
          call GetPointer_(CHEM_IMP_STATE, 'WET1'  ,c_wet1 , rc)
          call GetPointer_(CHEM_IMP_STATE, 'LAI'   ,c_lai  , rc)
          call GetPointer_(CHEM_IMP_STATE, 'GRN'   ,c_grn  , rc)
          call GetPointer_(CHEM_IMP_STATE, 'CN_PRCP',c_cn_prcp,rc)
          call GetPointer_(CHEM_IMP_STATE, 'NCN_PRCP',c_ncn_prcp,rc)
          call GetPointer_(CHEM_IMP_STATE, 'SH'    ,c_sh   , rc)
          call GetPointer_(CHEM_IMP_STATE, 'TA'    ,c_ta   , rc)
          call GetPointer_(CHEM_IMP_STATE, 'TSOIL1',c_tsoil1,rc)
          call GetPointer_(CHEM_IMP_STATE, 'U10M'  ,c_u10m , rc)
          call GetPointer_(CHEM_IMP_STATE, 'V10M'  ,c_v10m , rc)
          call GetPointer_(CHEM_IMP_STATE, 'USTAR' ,c_ustar, rc)
          call GetPointer_(CHEM_IMP_STATE, 'Z0H'   ,c_z0h  , rc)
          call GetPointer_(CHEM_IMP_STATE, 'TROPP' ,c_tropp, rc)
          call GetPointer_(CHEM_IMP_STATE, 'PS'    ,c_ps   , rc)

! for GOCART v2
          call GetPointer_(CHEM_IMP_STATE, 'AREA'  ,c_area , rc)
          call GetPointer_(CHEM_IMP_STATE, 'CLDTT' ,c_cldtt, rc)
          call GetPointer_(CHEM_IMP_STATE, 'FROCEAN',c_frocean, rc)
          call GetPointer_(CHEM_IMP_STATE, 'SWNDSRF',c_swndsrf, rc)
          call GetPointer_(CHEM_IMP_STATE, 'TS'     ,c_ts, rc)
          call GetPointer_(CHEM_IMP_STATE, 'U10N'   ,c_u10n, rc)
          call GetPointer_(CHEM_IMP_STATE, 'V10N'   ,c_v10n, rc)

! for GOCART, vertical index is fom toa to surface
          call GetPointer_3D_(CHEM_IMP_STATE,'PLE', c_ple , rc)
          call GetPointer_3D_(CHEM_IMP_STATE,'ZLE', c_zle , rc)
          call GetPointer_3D_(CHEM_IMP_STATE,'AIRDENS',  c_airdens, rc)
          call GetPointer_3D_(CHEM_IMP_STATE,'T'  , c_t   , rc)
          call GetPointer_3D_(CHEM_IMP_STATE,'U'  , c_u   , rc)
          call GetPointer_3D_(CHEM_IMP_STATE,'V'  , c_v   , rc)
          call GetPointer_3D_(CHEM_IMP_STATE,'FCLD', c_fcld , rc)
          call GetPointer_3D_(CHEM_IMP_STATE,'DQDT', c_dqdt , rc)

! for GOCART v2
!jw       call GetPointer_3D_(CHEM_IMP_STATE,'CH4' , c_ch4 , rc)
!jw       call GetPointer_3D_(CHEM_IMP_STATE,'H2O2', c_h2o2 , rc)
!jw       call GetPointer_3D_(CHEM_IMP_STATE,'NO3' , c_no3 , rc)
!jw       call GetPointer_3D_(CHEM_IMP_STATE,'OH'  , c_oh , rc)
          call GetPointer_3D_(CHEM_IMP_STATE,'CNV_MFC', c_cnv_mfc , rc)
          call GetPointer_3D_(CHEM_IMP_STATE,'CNV_MFD', c_cnv_mfd , rc)
          call GetPointer_3D_(CHEM_IMP_STATE,'CNV_QC' , c_cnv_qc , rc)
          call GetPointer_3D_(CHEM_IMP_STATE,'DQRL'   , c_dqrl , rc)

! get met tracers
!jw       call GetPointer_3D_(CHEM_IMP_STATE,'O3' , c_o3  , rc)
          call GetPointer_3D_(CHEM_IMP_STATE,'RH2' , c_rh2 , rc)

! chem tracers already pointed to phy export state (tracer bundle)

        endif

      ELSE

         print *, 'PHY2CHEM_RUN: chem import state ItemCount failed,', RC

      ENDIF

      t1_o=rtc()
      t1 = t1_o - t1_i
!
!---------------------------------------------
!* Do the actual coupling 
!---------------------------------------------
! 
      t2_i = rtc()
      IF ( RC_CPL == ESMF_SUCCESS ) THEN

        if (item_count_chem==0 .or. item_count_phys==0) then

          print *, 'PHY2CHEM_RUN: Empty state; Coupling phy_exp with chem_imp skipped'

        else


! --- fill in default values
          call FillDefault_

!
! --- Aerosol tracer fields: flip from bottom-up to top-down
!
          if ( run_DU ) then
            do k=1,km/2
              kk = km + 1 - k
              do j=1,jm
                lons_lat = lonsperlar_r(j)
                do i=1,lons_lat
                  tem = p_du001(i,j,kk)
                  p_du001(i,j,kk) = p_du001(i,j,k)
                  p_du001(i,j,k)  = tem

                  tem = p_du002(i,j,kk)
                  p_du002(i,j,kk) = p_du002(i,j,k)
                  p_du002(i,j,k)  = tem

                  tem = p_du003(i,j,kk)
                  p_du003(i,j,kk) = p_du003(i,j,k)
                  p_du003(i,j,k)  = tem

                  tem = p_du004(i,j,kk)
                  p_du004(i,j,kk) = p_du004(i,j,k)
                  p_du004(i,j,k)  = tem

                  tem = p_du005(i,j,kk)
                  p_du005(i,j,kk) = p_du005(i,j,k)
                  p_du005(i,j,k)  = tem
                enddo
              enddo
            enddo
          endif

          if ( run_SS ) then
            do k=1,km/2
              kk = km + 1 - k
              do j=1,jm
                lons_lat = lonsperlar_r(j)
                do i=1,lons_lat
                  tem = p_ss001(i,j,kk)
                  p_ss001(i,j,kk) = p_ss001(i,j,k)
                  p_ss001(i,j,k)  = tem

                  tem = p_ss002(i,j,kk)
                  p_ss002(i,j,kk) = p_ss002(i,j,k)
                  p_ss002(i,j,k)  = tem

                  tem = p_ss003(i,j,kk)
                  p_ss003(i,j,kk) = p_ss003(i,j,k)
                  p_ss003(i,j,k)  = tem

                  tem = p_ss004(i,j,kk)
                  p_ss004(i,j,kk) = p_ss004(i,j,k)
                  p_ss004(i,j,k)  = tem

                  tem = p_ss005(i,j,kk)
                  p_ss005(i,j,kk) = p_ss005(i,j,k)
                  p_ss005(i,j,k)  = tem
                enddo
              enddo
            enddo
          endif

          if ( run_OC ) then
            do k=1,km/2
              kk = km + 1 - k
              do j=1,jm
                lons_lat = lonsperlar_r(j)
                do i=1,lons_lat
                  tem = p_ocphobic(i,j,kk)
                  p_ocphobic(i,j,kk) = p_ocphobic(i,j,k)
                  p_ocphobic(i,j,k)  = tem


                  tem = p_ocphilic(i,j,kk)
                  p_ocphilic(i,j,kk) = p_ocphilic(i,j,k)
                  p_ocphilic(i,j,k)  = tem
                enddo
              enddo
            enddo
          endif

          if ( run_BC ) then
            do k=1,km/2
              kk = km + 1 - k
              do j=1,jm
                lons_lat = lonsperlar_r(j)
                do i=1,lons_lat
                  tem = p_bcphobic(i,j,kk)
                  p_bcphobic(i,j,kk) = p_bcphobic(i,j,k)
                  p_bcphobic(i,j,k)  = tem

                  tem = p_bcphilic(i,j,kk)
                  p_bcphilic(i,j,kk) = p_bcphilic(i,j,k)
                  p_bcphilic(i,j,k)  = tem
                enddo
              enddo
            enddo
          endif

          if ( run_SU ) then
            do k=1,km/2
              kk = km + 1 - k
              do j=1,jm
                lons_lat = lonsperlar_r(j)
                do i=1,lons_lat
                  tem = p_msa(i,j,kk)
                  p_msa(i,j,kk) = p_msa(i,j,k)
                  p_msa(i,j,k)  = tem

                  tem = p_so2(i,j,kk)
                  p_so2(i,j,kk) = p_so2(i,j,k)
                  p_so2(i,j,k)  = tem

                  tem = p_so4(i,j,kk)
                  p_so4(i,j,kk) = p_so4(i,j,k)
                  p_so4(i,j,k)  = tem

                  tem = p_dms(i,j,kk)
                  p_dms(i,j,kk) = p_dms(i,j,k)
                  p_dms(i,j,k)  = tem
                enddo
              enddo
            enddo
          endif

! --- 2D array: data copy
          do j=1,jm
            lons_lat = lonsperlar_r(j)
            do i=1,lons_lat
              c_frlake(i,j)  = 0.                            ! fraction_of_lake (1)
              c_lai(i,j)     = 3.                            ! leaf_area_index (1)
!
              c_zpbl(i,j)    = p_hpbl(i,j)                   ! boundary layer height (m)
              c_grn(i,j)     = p_vfrac(i,j)                  !  greeness_fraction (1)
              c_sh(i,j)      = p_dtsfci(i,j)                 ! sensible heat flux (W/m^2)
              c_ta(i,j)      = p_tsea(i,j)                   ! surface air Temperature (K)
              c_tsoil1(i,j)  = p_stc1(i,j)                   ! soil temperatures layer_1 (k)
              c_u10m(i,j)    = p_u10m(i,j)                   ! 10-meter eastward_wind (m s-1)
              c_v10m(i,j)    = p_v10m(i,j)                   ! 10-meter northward_wind (m s-1)
              c_ustar(i,j)   = p_ustar(i,j)                  ! surface velocity scale (m s-1)
              c_lwi(i,j)     = p_slmsk(i,j)                  !  land-ocean-ice mask  (1)
              c_ps(i,j)      = p_ps(i,j)                     ! surface pressure (Pa)
              c_wet1(i,j)    = p_wet1(i,j)                   ! soil wetness (1)
!
! --- 2D array: data copy with unit conversion
              c_cn_prcp(i,j) = 1.E3*p_rainc(i,j) /deltim     ! surface conv. rain flux (kg/m^2/s)
              c_z0h(i,j)     = p_zorl(i,j) / 1.E2            ! surface roughness (m)
!
! --- 2D array for ngacv2
              c_fraci(i,j)   = p_slmsk(i,j)*p_fice(i,j)      ! ice_covered_fraction_of_tile (1)
              c_frocean(i,j) = p_slmsk(i,j)                  ! ocean_covered_fraction_of_tile (1)
              c_u10n(i,j)    = 0.                            ! equivalent neutral u-wind at 10m
              c_v10n(i,j)    = 0.                            ! equivalent neutral v-wind at 10m

            enddo
          enddo
!
! --- get 2D array: area
          do j = 1, jm
            lons_lat = lonsperlar_r(j)
            if( lons_lat > 0 ) then 
              dx = con_rerth*(con_pi+con_pi)*sqrt(1.-slat_r(j)*slat_r(j))/lons_lat
              tem2 = con_rerth*dlat_r(j)
              do i = 1, lons_lat
                c_area(i,j) = dx*tem2
              enddo 
              print *,'in atm,phys-chem,j=',j,'dx=',dx,'c_area=',c_area(1,j),'glat=',dlat_r(j),slat_r(j)
            endif
          enddo 

! --- 3D array: filp vertical index from bottom-up to top-down
          do k=1,km
            kk = km + 1 - 1
            do j=1,jm
              lons_lat = lonsperlar_r(j)
              do i=1,lons_lat
                c_t (i,j,k)      = p_t (i,j,kk)       ! air temp at mid-layer (K)
                c_u (i,j,k)      = p_u (i,j,kk)       ! zonal wind at mid-layer (m/s)
                c_v (i,j,k)      = p_v (i,j,kk)       ! meridian wind at mid-layer (m/s)
!jw             c_o3(i,j,k)      = p_o3mr(i,j,kk)     ! ozone mixing ratio at mid-layer (kg/kg)
                c_fcld(i,j,k)    = p_fcld(i,j,kk)     ! cloud cover  (1)
                c_dqdt(i,j,k)    = p_dqdt(i,j,kk)     ! total moisture tendency (kg/kg/s)
                c_dqrl(i,j,k)    = 0.                 ! large scale rainwater source
                c_cnv_mfc(i,j,k) = p_cnv_mfc(i,j,kk)  ! cumulative mass flux (kg/m2/s)
                c_cnv_mfd(i,j,k) = p_cnv_mfd(i,j,kk)  ! detrained mass flux (kg/m2/s)
                c_cnv_qc(i,j,k)  = p_cnv_qc(i,j,kk)   ! total cpnvective condensat (kg/kg)

              enddo
            enddo
          enddo
!
! --- Compute ple, zle, airdens, rh2, tropp, wet1
          do j = 1, jm
            lons_lat = lonsperlar_r(j)
            do i = 1, im
!
              c_ncn_prcp(i,j) = 1.E3*(p_rain(i,j)-p_rainc(i,j))/deltim ! Non-conv. precip rate (kg/m^2/s)

!         local array
              sh(1:km) = max( p_spfh(i,j,1:km), qmin )

!         compute air pressure at interface (bottom-up)
              pi(0) = p_ps(i,j)
              do k=1,km-1
                pi(k) = pi(k-1) - p_dp(i,j,k)
              enddo
              pi(km) = f_zero

!         compute prsln (for thickness computation)
              do k = 0, km-1                    ! from SFC to TOA
                prsln(k) = log(0.01*pi(k))      ! convert from Pa to mb
              enddo
              prsln(km) = log(0.01*p_p(i,j,km)) ! convert from Pa to mb

!         compute rho (layer air density in kg/m^3), h (interface height in m)
              h(0) =  max ( p_hs(i,j), 0.0 )
!             h(0) =  f_zero
              do k = 1, km                                       ! from SFC to TOA
                tv1 = p_t(i,j,k) * (f_one + con_fvirt * sh(k))   ! virtual temp (k)
                rho(k) = p_p(i,j,k) /(con_rd * tv1)              ! air density (kg/m3)
                dz = rovg * (prsln(k-1)-prsln(k)) * tv1          ! thickness (m)
                if ( k == km ) dz = dz + dz
                h(k) = h(k-1) + dz                               ! geopotential height at interface (m)
              enddo

!         compute rh2 (relative humidity (in percent)
              call getrh(km,p_p(i,j,:),sh,p_t(i,j,:),shs,rh)

!         compute tropp (ptp: tropopause pressure in Pa)
              call tpause(km,p_p(i,j,:),p_u(i,j,:),p_v(i,j,:),p_t(i,j,:),h, &
                          ptp,utp,vtp,ttp,htp,shrtp)

!         pass local array to chem import
              c_tropp(i,j)    = ptp
              do k=0,km
                c_ple(i,j,k) = pi(km-k)                 ! air pressure at interface (Pa)
                c_zle(i,j,k) = h (km-k)                 ! geopotential height at interface (m)
              enddo
              do k=1,km
                c_airdens(i,j,k) = rho(km-k+1)          ! air density at mid-layer (kg/m3)
                c_rh2(i,j,k)     = rh(km-k+1) * 0.01    ! relative humidity in precent
              enddo

            enddo
          enddo

        endif

      ELSE

        print *, 'PHY2CHEM_RUN: Coupling chem_imp with phy_exp skipped'

      ENDIF
      t2_o=rtc()
      t2 = t2_o - t2_i
!
! --- now let's take care 2d aer_diag fields 
!
      BundleName='dgmet'
      do i = 1, nfld_met

        vname = name_met(i)
        nullify(p_diag)
        MESSAGE_CHECK = "Phys2Chem CPL_RUN: Get Farray from Phy_Exp-"//vname
        call GetPointer_diag_(PHY_EXP_STATE, BundleName, vname, p_diag, rc)
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_CPL)

        select case ( vname )
          case ('xU10M')
            p_diag = p_diag + deltim * c_u10m
          case ('xV10M')
            p_diag = p_diag + deltim * c_v10m
          case ('xUUSTAR')
            p_diag = p_diag + deltim * c_ustar
          case ('xZ0H')
            p_diag = p_diag + deltim * c_z0h
          case ('xLWI')
            p_diag = p_diag + deltim * c_lwi
          case ('xZPBL')
            p_diag = p_diag + deltim * c_zpbl
          case ('xWET1')
            p_diag = p_diag + deltim * c_wet1
          case ('xSH')
            p_diag = p_diag + deltim * c_sh
          case ('xCNPRCP')
            p_diag = p_diag + deltim * c_cn_prcp
          case ('xNCNPRCP')
            p_diag = p_diag + deltim * c_ncn_prcp

! change the following from accumulation to instant values;
! dqdt is scaled by 1e6

          case ('xZLE01')
            p_diag = c_zle(:,:,1)
          case ('xAIRDEN01')
            p_diag = c_airdens(:,:,1)
          case ('xT01')
            p_diag = c_t(:,:,1)
          case ('xU01')
            p_diag = c_u(:,:,1)
          case ('xV01')
            p_diag = c_v(:,:,1)
          case ('xFCLD01')
            p_diag = c_fcld(:,:,1)
          case ('xDQDT01')
            p_diag = c_dqdt(:,:,1) * 1.e6

          case ('xZLE64')
            p_diag = c_zle(:,:,64)
          case ('xAIRDEN64')
            p_diag = c_airdens(:,:,64)
          case ('xT64')
            p_diag = c_t(:,:,64)
          case ('xU64')
            p_diag = c_u(:,:,64)
          case ('xV64')
            p_diag = c_v(:,:,64)
          case ('xFCLD64')
            p_diag = c_fcld(:,:,64)
          case ('xDQDT64')
            p_diag = c_dqdt(:,:,64) * 1.e6

        end select   

      enddo    ! kcount-loop

!
!
!---------------------------------------------
!** Patch the array with valid values for chemistry import
!---------------------------------------------
        call patch_
!  ---  inputs:  (in scope variables)
!  ---  outputs: (in scope variables)

!! debug print
      if ( lckprnt ) then
       print *, 'RTC for Get Array', t1, t1_i, t1_o
       print *, 'RTC for Computations:', t2, t2_i, t2_o

       print *, 'PHY2CHEM_RUN: check chem_imp before exit PHY2CHEM_RUN'

       DO I = 1, nfld_2d
        call CkPointer_ (CHEM_IMP_STATE, vname_2d(I) , '2D', rc)
       ENDDO

       DO I = 1, nfld_3d
        call CkPointer_ (CHEM_IMP_STATE, vname_3d(I) , '3D',  rc)
       ENDDO

       DO I = 1, nfld_trc
        call CkPointer_ (CHEM_IMP_STATE, vname_trc(I) , 'TR', rc)
       ENDDO
      endif
!
!     Reset first flag
      FIRST = .False.
!
!-----------------------------------------------------------------------
!***  Check the final error signal variable 
!-----------------------------------------------------------------------
!
      IF(RC_CPL==ESMF_SUCCESS)THEN
        WRITE(0,*)'PHY2CHEM CPL RUN SUCCEEDED'
      ELSE
        WRITE(0,*)'PHY2CHEM CPL RUN FAILED RC_CPL=',RC_CPL
      ENDIF

!! 
      contains 

! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      subroutine FillDefault_

!  ---  inputs:  (in scope variables)
!  ---  outputs: (in scope variables)

        do j = 1, jm
          do i = 1, im

            c_PLE(i,j,:) = (/ 1, 2, 3, 4, 6, 8, 11, 15, 21, 27, 36, 47, 61, 79, 101, 130,       &
                           165, 208, 262, 327, 407, 504, 621, 761, 929, 1127, 1364, 1645,  &
                           1979, 2373, 2836, 3381, 4017, 4764, 5638, 6660, 7851, 9236,     &
                           10866, 12783, 15039, 17693, 20792, 24398, 28606, 33388, 37003,  &
                           40612, 44214, 47816, 51405, 54997, 58584, 62170, 65769, 68147,  &
                           70540, 72931, 75313, 77711, 79623, 81046, 82485, 83906, 85344,  &
                           86765, 88201, 89636, 91071, 92516, 93921, 95376 /)       

            c_ZLE(i,j,:) = (/ 78676, 74222, 71032, 68578, 66390, 64345, 62371, 60419, 58455, &
                            56469, 54463, 52449, 50446, 48476, 46563, 44718, 42946, 41256, &
                            39651, 38123, 36656, 35234, 33847, 32499, 31199, 29940, 28704, &
                            27494, 26310, 25151, 24017, 22905, 21815, 20745, 19691, 18656, &
                            17629, 16609, 15589, 14559, 13514, 12470, 11475, 10487, 9469, &
                            8438, 7731, 7076, 6463, 5889, 5348, 4838, 4355, 3898, 3464, &
                            3187, 2918, 2656, 2403, 2155, 1963, 1821, 1682, 1546, 1412, &
                            280, 1149, 1022, 896, 773, 654, 535, 417 /)

            c_AIRDENS(i,j,:) = (/ 2.27987766266e-05, 4.03523445129e-05, 6.19888305664e-05, 8.63075256348e-05, &
                                0.000117659568787, 0.000159025192261, 0.000209808349609, 0.000270366668701, &
                                0.000345230102539, 0.000439167022705, 0.00055980682373, 0.000717163085938, &
                                0.000923156738281, 0.00120162963867, 0.00156402587891, 0.00202178955078, &
                                0.00262451171875, 0.00339889526367, 0.00437164306641, 0.00555419921875, &
                                0.00694274902344, 0.00857543945312, 0.0105895996094, 0.0131225585938, &
                                0.0160827636719, 0.0195617675781, 0.0237731933594, 0.0287780761719, &
                                0.0347290039062, 0.0416870117188, 0.0499267578125, 0.0596313476562, &
                                0.0711669921875, 0.084716796875, 0.100830078125, 0.11865234375, 0.138671875, &
                                0.1630859375, 0.190185546875, 0.22021484375, 0.25927734375, 0.318359375, &
                                0.3720703125, 0.42138671875, 0.47265625, 0.521484375, 0.5615234375, &
                                0.6005859375, 0.638671875, 0.677734375, 0.71875, 0.759765625, 0.8017578125, &
                                0.8447265625, 0.8798828125, 0.90625, 0.9326171875, 0.958984375, 0.986328125, &
                                1.013671875, 1.03515625, 1.052734375, 1.072265625, 1.08984375, 1.10546875, &
                                1.123046875, 1.140625, 1.162109375, 1.1953125, 1.21875, 1.234375, 1.25 /)

            c_FCLD(i,j,:) = 1e-2 * &
                          (/ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
                             0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
                             0, 0, 0, 0, 0, 0, 16, 21, 26, 28, 0, 0, 0, 0, 0, 0, 0, &
                             0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0 /)

            c_DQDT(i,j,:) = 1e-12 * &
                          (/ 9, 11, -3, -3, -2, -18, -10, 2, 0, -3, -6, -5, -3, -1, 1, &
                             1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, &
                             1, 1, 0, 0, -5, -22, -33, 95, 474, 348, 177, 3377, 11045, &
                             11788, -5267, -7756, -17491, -19790, -10884, -6082, 8120, 4381, &
                             -10346, 8033, 69151, 77650, 61351, 46508, 33936, 23022, 15658, &
                             11598, 6469, 4861, -846, -7974, -30500, -20663, -14930 /)

            c_T(i,j,:) = (/ 219, 221, 223, 228, 230, 230, 232, 238, 245, 253, 259, 263, &
                          264, 262, 258, 253, 247, 239, 233, 229, 227, 227, 226, 223, &
                          222, 221, 220, 219, 218, 217, 216, 215, 214, 213, 212, 212, &
                          214, 214, 216, 219, 219, 210, 210, 218, 227, 234, 240, 245, &
                          250, 254, 257, 260, 262, 263, 265, 266, 267, 268, 269, 270, &
                          270, 270, 270, 270, 271, 271, 271, 270, 267, 265, 266, 266 /)


            c_U(i,j,:) = (/ -18, -13, 0, 10, 26, 36, 39, 40, 38, 37, 36, 35, 32, 28,    &
                           23, 16, 6, -2, -9, -13, -15, -16, -14, -14, -12, -12, -11, &
                          -10, -9, -5, -3, -2, 0, 1, 3, 5, 9, 13, 17, 22, 24, 26,     &
                           25, 26, 26, 22, 19, 17, 14, 12, 12, 11, 11, 11, 11, 10, 9, &
                            8, 6, 4, 3, 2, 1, 0, -1, -2, -3, -4, -5, -6, -6, -6 /)

            c_V(i,j,:) = (/ 20, 13, 9, 4, -1, -9, -20, -24, -25, -27, -28, -28, -26,    &
                         -25, -27, -28, -28, -28, -27, -27, -25, -23, -19, -15, -11,  &
                         -10, -9, -8, -7, -7, -8, -9, -10, -12, -14, -15, -16, -18,   &
                         -21, -22, -22, -25, -29, -25, -23, -23, -22, -20, -17, -13,  &
                         -9, -6, -4, -4, -4, -3, -2, -1, 0, 0, 0, 1, 1, 1, 2, 2,      &
                          3, 3, 3, 4, 4, 3 /)

!jw            c_O3(i,j,:) = 1.E-9 * & 
!jw                        (/ 16182, 9700, 7294, 5781, 4164, 3017, 2440, 2287, 2324, 2514,  &
!jw                           2838, 3304, 4030, 4924, 5915, 7033, 8434, 9894, 11101, 11414, &
!jw                           10475, 9745, 10058, 9119, 8538, 9238, 9164, 10028, 10132, 10237, &
!jw                           9447, 7972, 7174, 5222, 4008, 3296, 2231, 1320, 768, 628, 685, &
!jw                           676, 202, 122, 96, 88, 86, 83, 83, 84, 84, 83, 82, 81, 79, &
!jw                           79, 77, 76, 77, 80, 84, 87, 89, 90, 89, 88, 83, 76, 69, 65, &
!jw                           64, 64 /)
 
            c_RH2(i,j,:) = 1e-6 * &
                         (/ 1, 2, 2, 2, 3, 4, 4, 3, 4, 4, 4, 4, 4, 4, 4, 6, 18, 51,          &
                            129, 267, 394, 502, 682, 1135, 1603, 2076, 2820, 3792, 5120,     &
                            6806, 8912, 11597, 15397, 20386, 28168, 29755, 28748, 33875,     &
                            34058, 28657, 43458, 401856, 947266, 932618, 902344, 657227,     &
                            371583, 203370, 235108, 317872, 413086, 511719, 691407, 686524,  &
                            601563, 456055, 475098, 626954, 590821, 483399, 380860, 297852,  &
                            230958, 183594, 144288, 111084, 96558, 136963, 369629, 770508,   &
                            793946, 799805 /)
!           2D
!           --
            c_TROPP(i,j)    = 20363.5
            c_LWI(i,j)      = 1.
            c_ZPBL(i,j)     = 59.
            c_FRLAKE(i,j)   = 0.
            c_FRACI(i,j)    = 0.
            c_WET1(i,j)     = 0.0
            c_LAI(i,j)      = 0.280273
            c_GRN(i,j)      = 0.5
            c_CN_PRCP(i,j)  = 0.0
            c_NCN_PRCP(i,j) = 3.18323e-10
            c_PS(i,j)       = 96825.3
            c_SH(i,j)       = -28.548
            c_TSOIL1(i,j)   = 260.014
            c_U10M(i,j)     = -3.5
            c_V10M(i,j)     = 2.8
            c_USTAR(i,j)    = 0.29
            c_Z0H(i,j)      = 0.02005
            c_TA(i,j)       = 270.014

        end do
      end do

      end subroutine FillDefault_

      subroutine patch_

!  ---  inputs:  (in scope variables)
!  ---  outputs: (in scope variables)


      do j = 1, lats_node_r_max
        lons_lat = lonsperlar_r(j)
        if (lons_lat  < lonr ) then
          do i=lons_lat+1,lonr           ! fill in 2D array except for frlake, fraci, lai
            c_lwi  (i,j)    = c_lwi   (1,1)
            c_zpbl (i,j)    = c_zpbl  (1,1)
            c_wet1 (i,j)    = c_wet1  (1,1)
            c_grn  (i,j)    = c_grn   (1,1)
            c_cn_prcp(i,j)  = c_cn_prcp(1,1)
            c_ncn_prcp(i,j) = c_ncn_prcp(1,1)
            c_sh   (i,j)    = c_sh    (1,1)
            c_ta   (i,j)    = c_ta    (1,1)
            c_tsoil1(i,j)   = c_tsoil1(1,1)
            c_u10m (i,j)    = c_u10m  (1,1)
            c_v10m (i,j)    = c_v10m  (1,1)
            c_ustar(i,j)    = c_ustar (1,1)
            c_z0h  (i,j)    = c_z0h   (1,1)
            c_tropp(i,j)    = c_tropp (1,1)
            c_ps   (i,j)    = c_ps    (1,1)
!
            c_frocean   (i,j)    = c_frocean    (1,1)
            c_fraci     (i,j)    = c_fraci      (1,1)
            c_area      (i,j)    = c_area       (1,1)
            c_u10n      (i,j)    = 0. 
            c_v10n      (i,j)    = 0. 
          enddo
        endif
      enddo

      do k = 0, km                       ! fill in 3D array at interface
        do j = 1, lats_node_r_max
          lons_lat = lonsperlar_r(j)
          if (lons_lat  < lonr ) then
            do i=lons_lat+1,lonr
              c_ple(i,j,k) = c_ple(1,1,k)
              c_zle(i,j,k) = c_zle(1,1,k)
             enddo
          endif
        enddo
      enddo
      do k = 1, km                       ! fill in 3D array at mid-layer
        do j = 1, lats_node_r_max
          lons_lat = lonsperlar_r(j)
          if (lons_lat  < lonr ) then
            do i=lons_lat+1,lonr
              c_airdens(i,j,k) = c_airdens(1,1,k)
              c_t      (i,j,k) = c_t      (1,1,k)
              c_u      (i,j,k) = c_u      (1,1,k)
              c_v      (i,j,k) = c_v      (1,1,k)
              c_fcld   (i,j,k) = c_fcld   (1,1,k)
              c_dqdt   (i,j,k) = c_dqdt   (1,1,k)
              c_cnv_mfc(i,j,k) = c_cnv_mfc(1,1,k)
              c_cnv_mfd(i,j,k) = c_cnv_mfd(1,1,k)
              c_cnv_qc (i,j,k) = c_cnv_qc (1,1,k)
              c_dqrl   (i,j,k) = 0.              

!jw              c_o3     (i,j,k)  = c_o3    (1,1,k)
              c_rh2    (i,j,k)  = c_rh2   (1,1,k)
            enddo
          endif
        enddo
      enddo

      RETURN

      end subroutine patch_

      END subroutine run

!=========================
      subroutine get_attribute(STATE, RC_CPL)

!  ---  input
        type(ESMF_State)                     :: STATE

!  ---  output
        integer, intent(out)                 :: RC_CPL

!  ---  locals:
        type(ESMF_FieldBundle)   :: Bundle
        integer                  :: STATUS, RC
        character(esmf_maxstr)   :: statename

! ---   Retrieve statename
        MESSAGE_CHECK='Retrive state name'
        call ESMF_StateGet(state=STATE, name=statename, rc=RC )
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_CPL)

!  ---  Retrieve lat/lon info
        MESSAGE_CHECK = 'Extract lonr attribute from '//trim(statename)
        CALL ESMF_AttributeGet(state=STATE, name='lonr'  &
                              ,value=lonr, rc=RC)
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_CPL)

        MESSAGE_CHECK = 'Extract lats_node_r attribute from '//trim(statename)
        CALL ESMF_AttributeGet(state=STATE, name='lats_node_r'  &  
                             ,value=lats_node_r, rc=RC)
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_CPL)

        MESSAGE_CHECK = 'Extract lats_node_r_max attribute from '//trim(statename)
        CALL ESMF_AttributeGet(state=STATE, name='lats_node_r_max'&
                             ,value=lats_node_r_max, rc=RC)
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_CPL)

        if ( .not. allocated (lonsperlar_r)) then
          allocate ( lonsperlar_r(lats_node_r_max))
        endif

        MESSAGE_CHECK = 'Extract lonsperlar_r attribute from '//trim(statename)

        CALL ESMF_AttributeGet(state  =STATE               &  !<-- Name of the state
                           ,name      ='lonsperlar_r'      &  !<-- Name of the attribute to retrieve
                           ,itemCount = lats_node_r_max    &  !<-- Number of values in the attribute
                           ,valueList =lonsperlar_r        &  !<-- Value of the attribute
                           ,rc        =RC)

        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_CPL)

        if ( .not. allocated (slat_r)) then
          allocate ( slat_r(lats_node_r_max))
          allocate ( dlat_r(lats_node_r_max))
        endif

        MESSAGE_CHECK = 'Extract slat_r attribute from '//trim(statename)

        CALL ESMF_AttributeGet(state  =STATE               &  !<-- Name of the state
                           ,name      ='slat_r'            &  !<-- Name of the attribute to retrieve
                           ,itemCount = lats_node_r_max    &  !<-- Number of values in the attribute
                           ,valueList =slat_r              &  !<-- Value of the attribute
                           ,rc        =RC)
!        print *,'in phy_chem,get slat_r,slat_r=',slat_r

        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_CPL)

        MESSAGE_CHECK = 'Extract dlat_r attribute from '//trim(statename)

        CALL ESMF_AttributeGet(state  =STATE               &  !<-- Name Cof the state
                           ,name      ='dlat_r'            &  !<-- Name of the attribute to retrieve
                           ,itemCount = lats_node_r_max    &  !<-- Number of values in the attribute
                           ,valueList =dlat_r              &  !<-- Value of the attribute
                           ,rc        =RC)
!       print *,'in phy_chem,get dlat_r,dlat_r=',dlat_r

        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_CPL)


! ---   Retrieve tracer specification
        MESSAGE_CHECK = 'Extract tracer bundle from '//trim(statename)
        call ESMF_StateGet(state=STATE, ItemName='tracers',    &
                           fieldbundle=Bundle, rc = rc )
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_CPL)

        MESSAGE_CHECK = 'Extract ntrac attribute from '//trim(statename)
        CALL ESMF_AttributeGet(Bundle, name  ='ntrac'  & !<-- Name of the attribute to retrieve
                              ,value = ntrac, rc = RC)   !<-- Value of the attribute
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_CPL)

       RETURN
      end subroutine get_attribute

!!! ---------------- ! ------------------ ! ---------------- !----------------!

      subroutine CkPointer_(STATE,NAME,TAG,rc)

! --- input/output arguments
        type(ESMF_State)                :: STATE
        character(len=*)                :: NAME
        character(len=*)                :: TAG
        integer, intent (OUT)           :: rc

! --- locals
        integer                         :: rc1, ii, jj, kk
        type(ESMF_Field)                :: Field
        type(ESMF_FieldBundle)          :: Bundle
        character(esmf_maxstr)          :: StateName, BundleName, LName
        real(ESMF_KIND_R8), pointer     :: Array2D(:,:), Array3D(:,:,:)
!
        MESSAGE_CHECK='Retrive state name'
        call ESMF_StateGet(state=State, name=statename, rc=rc1)
        CALL ERR_MSG(rc1, MESSAGE_CHECK, rc)

!       select case ( lowercase(statename))
        call g_fixchar(StateName, LName, 2)
        select case ( Lname  )
        case ( 'physics export') 
           BundleName='tracers'
        case ( 'chemistry import') 
           BundleName='iAERO'
        case ( 'chemistry export') 
           BundleName='AERO'
        end select
!
        if ( TAG == '2D' ) then
          MESSAGE_CHECK = 'Extract '//NAME//' from '//statename
          call ESMF_StateGet(state=STATE, itemName=NAME, field=Field, rc=rc1)
          CALL ERR_MSG(rc1, MESSAGE_CHECK, rc)

          nullify(Array2D)
          MESSAGE_CHECK = 'Get 2d Fortran data pointer from '//NAME

          CALL ESMF_FieldGet(field=Field, localDe=0, farrayPtr=Array2D, rc=rc1)

          CALL ERR_MSG(rc1, MESSAGE_CHECK, rc)
!
          ii = size(Array2D, dim=1)
          jj = size(Array2D, dim=2)
!
!         print *, trim(statename), '/', NAME,':', Array2D(1,1),Array2D(ii,jj), &
!                           minval(Array2D),maxval(Array2D)

        ELSEIF (TAG == '3D' ) then
          MESSAGE_CHECK = 'Extract '//NAME//' from '//statename
          call ESMF_StateGet(state=STATE, itemName=NAME, field=Field, rc=rc1)
          CALL ERR_MSG(rc1, MESSAGE_CHECK, rc)

          nullify(Array3D)
          MESSAGE_CHECK = 'Get 3d Fortran data pointer from '//NAME

          CALL ESMF_FieldGet(field=Field, localDe=0, farrayPtr=Array3D, rc=rc1)

          CALL ERR_MSG(rc1, MESSAGE_CHECK, rc)
!
          ii = size(Array3D, dim=1)
          jj = size(Array3D, dim=2)
          kk = size(Array3D, dim=3)
!
!         print *, trim(statename), '/', NAME,':', Array3D(1,1,1),Array3D(ii,jj,kk), &
!                           minval(Array3D),maxval(Array3D)

        ELSEIF (TAG == 'TR' ) then
          MESSAGE_CHECK = 'Extract '//trim(BundleName)//' from '//statename
          call ESMF_StateGet(state=State, ItemName=BundleName, fieldbundle=Bundle, rc=rc1)
          CALL ERR_MSG(rc1, MESSAGE_CHECK, rc)
!
          MESSAGE_CHECK = 'Extract '//NAME//' from '//trim(BundleName)
          CALL ESMF_FieldBundleGet(Bundle, name, field=Field, rc=rc1)
          CALL ERR_MSG(rc1, MESSAGE_CHECK, rc)

!
!
          nullify(Array3D)
          MESSAGE_CHECK = 'Get 3d Fortran data pointer from '//NAME

          CALL ESMF_FieldGet(field=Field, localDe=0, farrayPtr=Array3D, rc=rc1)

          CALL ERR_MSG(rc1, MESSAGE_CHECK, rc)
!
          ii = size(Array3D, dim=1)
          jj = size(Array3D, dim=2)
          kk = size(Array3D, dim=3)
!
!         print*, trim(statename), '/', NAME,':', Array3D(1,1,1),Array3D(ii,jj,kk), &
!                           minval(Array3D),maxval(Array3D)

        ELSE

          print *, 'PHY2CHEM_RUN: Unknown data type, abort now'
          call abort

        ENDIF


        return
!
      end subroutine CkPointer_

!!! ---------------- ! ------------------ ! ---------------- !----------------!

        subroutine GetPointer_ (STATE, NAME, Array, RC)

! --- input/output arguments
        type(ESMF_State)                :: State
        character(len=*), intent(IN)    :: NAME
        real(ESMF_KIND_R8), pointer, intent(OUT)  :: Array(:,:)
        integer, intent (OUT)           :: rc

! --- locals
        type(ESMF_Field)                :: Field
        integer                         :: rc1
        character(esmf_maxstr)          :: statename
!
!===>  ...  begin here
!
        MESSAGE_CHECK = 'PHY2CHEM_RUN: Retrive statename'
        call ESMF_StateGet(state=State, name=statename, rc=rc1)
        CALL ERR_MSG(rc1, MESSAGE_CHECK, rc)

        MESSAGE_CHECK = 'PHY2CHEM_RUN: Extract '//NAME//' from '//statename
        call ESMF_StateGet(state=STATE, itemName=NAME, field=Field, rc=rc1)
        CALL ERR_MSG(rc1, MESSAGE_CHECK, rc)
!
        nullify(Array)
        MESSAGE_CHECK = 'PHY2CHEM_RUN: Get Fortran data pointer from '//NAME

        CALL ESMF_FieldGet(field=Field, localDe=0, farrayPtr=Array, rc=rc1)

        CALL ERR_MSG(rc1, MESSAGE_CHECK, rc)
!
!       check i- and j-dimension
        if ( lonr .ne. size(Array, dim=1) )  print *, 'ERROR !',   &
             'Invalid lonr:',  lonr, size(Array,dim=1)
        if (lats_node_r_max .ne. size(Array, dim=2)) print *, 'ERROR !', &
             'Invalid lats_node_r_max', lats_node_r_max, size(Array,dim=2)

       end subroutine GetPointer_

!!! ---------------- ! ------------------ ! ---------------- !----------------!

        subroutine GetPointer_3D_ (STATE, NAME, Array, RC)

! --- input/output arguments
        type(ESMF_State)                :: State
        character(len=*), intent(IN)    :: NAME
        real(ESMF_KIND_R8), pointer, intent(OUT)  :: Array(:,:,:)
        integer, intent (OUT)           :: rc

! --- locals
        type(ESMF_Field)                :: Field
        integer                         :: rc1
        character(esmf_maxstr)          :: statename
!
!===>  ...  begin here
!
        MESSAGE_CHECK = 'PHY2CHEM_RUN: Retrive statename'
        call ESMF_StateGet(state=State, name=statename, rc=rc1 )
        CALL ERR_MSG(rc1, MESSAGE_CHECK, rc)

        MESSAGE_CHECK = 'PHY2CHEM_RUN: Extract '//NAME//' from '//statename
        call ESMF_StateGet(state=STATE, itemName=NAME, field=Field, rc=rc1)
        CALL ERR_MSG(rc1, MESSAGE_CHECK, rc)
!
        nullify(Array)
        MESSAGE_CHECK = 'PHY2CHEM_RUN: Get Fortran data pointer from '//NAME

        CALL ESMF_FieldGet(field=Field, localDe=0, farrayPtr=Array, rc=rc1)

        CALL ERR_MSG(rc1, MESSAGE_CHECK, rc)
!
!       check i- and j-dimension
        if ( lonr .ne. size(Array, dim=1) )  print *, 'ERROR !',   &
             'Invalid lonr:',  lonr, size(Array,dim=1)
        if (lats_node_r_max .ne. size(Array, dim=2)) print *, 'ERROR !', &
             'Invalid lats_node_r_max', lats_node_r_max, size(Array,dim=2)

       end subroutine GetPointer_3D_

!!! ---------------- ! ------------------ ! ---------------- !----------------!

        subroutine GetPointer_tracer_ (STATE, NAME, Array, RC)

! --- input/output arguments
        type(ESMF_State)                :: State
        character(len=*), intent(IN)    :: NAME
        real(ESMF_KIND_R8), pointer, intent(OUT) :: Array(:,:,:)
        integer, intent (OUT)           :: rc

! --- locals
        type(ESMF_Field)                :: Field
        type(ESMF_FieldBundle)          :: Bundle
        integer                         :: rc1
        character(esmf_maxstr)          :: statename, BundleName, Lname
!
        character*8                     :: FldName(8)
        integer                         :: nameCount, i
!===>  ...  begin here

        MESSAGE_CHECK = 'PHY2CHEM_RUN: Retrive statename'
        call ESMF_StateGet(state=State, name=statename, rc=rc1 )
        CALL ERR_MSG(rc1, MESSAGE_CHECK, rc)

!       select case ( lowercase(statename))
        call g_fixchar(StateName, LName, 2)
        select case ( Lname  )
        case ( 'physics export') 
           BundleName='tracers'
        case ( 'chemistry import') 
           BundleName='iAERO'
        case ( 'chemistry export') 
           BundleName='AERO'
        end select

        MESSAGE_CHECK = 'PHY2CHEM_RUN: Extract '//BundleName//' from '//statename
        call ESMF_StateGet(state=State, ItemName=BundleName, fieldbundle=Bundle, rc=rc1)
        CALL ERR_MSG(rc1, MESSAGE_CHECK, rc)

        MESSAGE_CHECK = 'PHY2CHEM_RUN: Extract '//NAME//' from '//BundleName
!jw        CALL ESMF_FieldBundleGet(Bundle, FIELDNAME=name, field=Field, rc=rc1)
        CALL ESMF_FieldBundleGet(Bundle, name, field=Field, rc=rc1)
        CALL ERR_MSG(rc1, MESSAGE_CHECK, rc)

        nullify(Array)
        MESSAGE_CHECK = 'PHY2CHEM_RUN: Get Fortran data pointer from '//NAME

        CALL ESMF_FieldGet(field=Field, localDe=0, farrayPtr=Array, rc=rc1)

        CALL ERR_MSG(rc1, MESSAGE_CHECK, rc)
!
!       check i- and j-dimension
        if ( lonr .ne. size(Array, dim=1) )  print *, 'ERROR !',   &
             'Invalid lonr:',  lonr, size(Array,dim=1)
        if (lats_node_r_max .ne. size(Array, dim=2)) print *, 'ERROR !', &
             'Invalid lats_node_r_max', lats_node_r_max, size(Array,dim=2)

       end subroutine GetPointer_tracer_

!!! ---------------- ! ------------------ ! ---------------- !----------------!

        subroutine GetPointer_diag_ (STATE, BUNDLENAME, NAME, ARRAY, RC)

! --- input/output arguments
        type(ESMF_State)                :: STATE
        character(len=*), intent(IN)    :: BUNDLENAME
        character(len=*), intent(IN)    :: NAME
        real(ESMF_KIND_R8), pointer, intent(OUT) :: ARRAY(:,:)
        integer, intent (OUT)           :: RC

! --- locals
        type(ESMF_Field)                :: Field
        type(ESMF_FieldBundle)          :: Bundle
        integer                         :: rc1
!
!===>  ...  begin here

        IF ( BundleName(1:4) == 'xxxx') then
         MESSAGE_CHECK = 'GetPointer_diag: Extract Field '//NAME
         CALL ESMF_StateGet(state=State, ItemName=NAME, field=Field, rc=rc1)
         CALL ERR_MSG(rc1, MESSAGE_CHECK, rc)
        ELSE
         MESSAGE_CHECK = 'GetPointer_diag: Extract Bundle '//BundleName
         call ESMF_StateGet(state=State, ItemName=BundleName,  &
                           fieldbundle=Bundle, rc=rc1)
         CALL ERR_MSG(rc1, MESSAGE_CHECK, rc)

         MESSAGE_CHECK = 'GetPointer_diag:: Extract Field '//NAME
         CALL ESMF_FieldBundleGet(Bundle, name, field=Field, rc=rc1)
         CALL ERR_MSG(rc1, MESSAGE_CHECK, rc)
        ENDIF

        nullify(Array)
        MESSAGE_CHECK = 'GetPointer_diag:: Get Fortran data pointer from '//NAME

        CALL ESMF_FieldGet(field=Field, localDe=0, farrayPtr=Array, rc=rc1)

        CALL ERR_MSG(rc1, MESSAGE_CHECK, rc)
!
       end subroutine GetPointer_diag_

!! adopt getrh routine from /nwprod/sorc/global_postgp.fd/postgp.f

      subroutine getrh(km,p,sh,t,shs,rh)
!
! Subprogram: getrh      Compute saturation humidity and relative humidity
!   Prgmmr: Iredell      Org: np23        Date: 1999-10-18
!
! Abstract: This subprogram computes the saturation specific humidity and the
!           relative humidity.  The relative humidity is constrained to be
!           between 0 and 100.
!
! Program history log:
!   1999-10-18  Mark Iredell
!
! Usage:  call getrh(km,p,sh,t,shs,rh)
!   Input argument list:
!     km       integer number of levels
!     p        real (km) pressure (Pa)
!     sh       real (km) specific humidity (kg/kg)
!     t        real (km) temperature (K)
!   Output argument list:
!     shs      real (km) saturation specific humidity (kg/kg)
!     rh       real (km) relative humidity (percent)
!
! Modules used:
!   funcphys       Physical functions
!
! Files included:
!   physcons.h     Physical constants
!
! Subprograms called:
!   fpvs           compute saturation vapor pressure
!
! Attributes:
!   Language: Fortran 90
!
!$$$
       implicit none
       integer,intent(in):: km
       real,intent(in):: p(km),sh(km),t(km)
       real,intent(out):: shs(km),rh(km)
       real(krealfp) pr,tr,es
       integer k
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       do k=1,km
         pr     = p(k)
         tr     = t(k)
         es     = fpvs(tr)
         es     = min(es,pr)
         shs(k) = con_eps*es/(pr+con_epsm1*es)
         rh(k)  = 1.e2*min(max(sh(k)/shs(k),0.),1.)
       enddo
       end subroutine
!
      subroutine g_fixchar(name_in, name_out, option)
      implicit none

      character*(*), intent(in)   ::  name_in
      character*(*), intent(out)  ::  name_out
      integer, intent(in)         ::  option

      character*30                :: temp
      integer                     :: i, ic

      name_out= '          '
      temp = trim(adjustl(name_in))
      do i = 1, len_trim(temp)
        ic = IACHAR(temp(i:i))
        if(option == 1 ) then             !<--- convert to upper case
          if(ic .ge. 97 .and. ic .le. 122) then
            name_out(i:i) = CHAR( IC-32 )
          else
            name_out(i:i) = temp(i:i)
          endif
        endif
        if(option == 2 ) then             !<--- convert to lower case
          if(ic .ge. 65 .and. ic .le. 90) then
            name_out(i:i) = CHAR( IC+32 )
          else
            name_out(i:i) = temp(i:i)
          endif
        endif

      enddo
      name_out = trim(name_out)
      return

      end subroutine g_fixchar
!
!
      END module atmos_phy_chem_cpl_comp_mod



