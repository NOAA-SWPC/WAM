#include "../../../ESMFVersionDefine.h"

!
      module atmos_chem_phy_cpl_comp_mod

!-----------------------------------------------------------------------
!
!** This module holds the chem-to-phys coupler's register and run routines
!** setservices (only registers run step) is called by GOCART_SETUP
!** run transfer/convert data from chem export state to phy export state
!
!! Code Revision:
!! 24Feb 2010     Sarah Lu,   First Crack
!! 12Mar 2010     Sarah Lu,   Use routines from phy2chem coupler
!! 16Mar 2010     Sarah Lu,   Dimension and tracer specification are passed 
!!                            in from phy2chem coupler module; flip the 
!!                            vertical index from top-down to bottom-up
!! 09May 2010     Sarah Lu,   Revise species name for SU, OC, and BC
!! 13May 2010     Sarah Lu,   Gaseous species (DMS, SO2, MSA) are extracted
!!                            from GOCART export state (not from AERO bundle)
!! 10Jun 2010     Sarah Lu,   Gaseous species are taken from AERO bundle (as
!!                            GOCART grid component is revised)
!! 06Aug 2010     Sarah Lu,   Modify phy2chem run routine to pass g2d_fld
!!                            from chem_exp to phys_exp 
!! 10Aug 2010     Sarah Lu,   Modify chem2phy run routine to accumulate g2d_fld
!! 10Oct 2010     Sarah Lu,   Move GetPointer_diag_ to phy2chem coupler
!! 2011-05-11     Weiyu Yang, Modified for using the ESMF 5.2.0r_beta_snapshot_07.
!! 2011-11-27     Sarah Lu,   Specify i, j dimension in mapping aerosol arrays
!! 2012-02-06     Weiyu Yang, Modified for using the ESMF 5.2.0r library.
!! 2014-11-03     S Moorthi   Some cleanup and optimization
!! 17Dec 2014     Jun Wang,   retrieve aerosol pointer once and do aerosol 
!!                            flipflop on the pointer itself
!! 18Dec 2014     Jun Wang,   reset diagnostic fields
!------------------------------------------------------------------------------

      USE ESMF

      USE MODULE_ERR_MSG, ONLY: ERR_MSG, MESSAGE_CHECK
      use MODULE_gfs_machine,  ONLY: kind_phys
!
      use atmos_phy_chem_cpl_comp_mod, only:                  &
                          GetPointer_tracer_, CkPointer_,     &
                          GetPointer_diag_,                   &
                          lonr, lats_node_r, lats_node_r_max, &
                          lonsperlar_r, im, jm, km, ntrac,    &
                          run_DU, run_SU, run_SS, run_OC, run_BC

!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
      public :: setservices

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
!! This routine register the coupler component's run routine        
!!
!! Code Revision:
!! 24Feb 2009     Sarah Lu, First Crack
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
!***  register the coupler component's run routine
!-----------------------------------------------------------------------
!
      MESSAGE_CHECK="Set Entry Point for chem2phy coupler run"

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
        WRITE(0,*)'CHEM2PHY CPL SET_SERVICES SUCCEEDED'
      ELSE
        WRITE(0,*)'CHEM2PHY CPL SET_SERVICES FAILED RC_REG=',RC_REG
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
      subroutine run(GC, CHEM_EXP_STATE, PHY_EXP_STATE, CLOCK, RC_CPL)
!
!-----------------------------------------------------------------------
!!
!! This routine transfer/convert data from chem_exp state to phy_exp state
!!
!! Code Revision:
!! 24Feb 2010     Sarah Lu, First Crack
!! 24Aug 2012     Sarah Lu, Remove item_name
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!***  argument variables
!-----------------------------------------------------------------------

      type(ESMF_cplcomp)               :: GC
      type(ESMF_state)                 :: CHEM_EXP_STATE   ! coupler import state
      type(ESMF_state)                 :: PHY_EXP_STATE    ! coupler export state
      type(ESMF_clock)                 :: CLOCK
!
      integer,           intent(out)   :: RC_CPL

!-----------------------------------------------------------------------
!***  local variables
!-----------------------------------------------------------------------
!
      integer                          :: rc=ESMF_success  ! the error signal variable
      integer                          :: item_count_phys, item_count_chem
      logical, save                    :: first =  .true.

! Fortran array for phy export state
      real (ESMF_KIND_R8), pointer, dimension(:,:,:) ::                  &
               p_du001, p_du002, p_du003, p_du004, p_du005,              & ! DU
               p_ss001, p_ss002, p_ss003, p_ss004, p_ss005,              & ! SS
               p_msa,   p_so4,   p_so2,   p_dms,                         & ! SU
               p_ocphobic, p_ocphilic, p_bcphobic, p_bcphilic              ! OC/BC
           
! Fortran array for chem export state
      real (ESMF_KIND_R8), pointer, dimension(:,:,:) ::                  &
               c_du001, c_du002, c_du003, c_du004, c_du005,              & ! DU
               c_ss001, c_ss002, c_ss003, c_ss004, c_ss005,              & ! SS
               c_msa,   c_so4,   c_so2,   c_dms,                         & ! SU
               c_ocphobic, c_ocphilic, c_bcphobic, c_bcphilic              ! OC/BC
!

!---  Add the following for 2d aerosol diag fields ---

!  Fortran data pointer for phy export state
      real(ESMF_KIND_R8), pointer, dimension(:,:) ::  p_diag

!  Fortran data pointer for chem export state
      real(ESMF_KIND_R8), pointer, dimension(:,:) ::  c_diag

      logical                 :: get_attribute
      integer                 :: i, j, k, kk, kcount, lons_lat
      character*10            :: aerosol_list(5), aerosol, tag, vname
      character*10            :: BundleName, FieldName
      TYPE(ESMF_FieldBundle)  :: Bundle
      character*10, dimension(30)  :: name_lst
      real                       tem
      real, save              :: deltim
      integer,save            :: nfld_du, nfld_ss, nfld_su, nfld_oc, nfld_bc
      character*10,  allocatable, save  :: name_du(:), name_ss(:), &
                                           name_oc(:), name_bc(:), name_su(:)

      data aerosol_list / 'du', 'su', 'ss', 'oc', 'bc'/


!---------------------------------------------
!* Determine dimension and allocate local array 
!---------------------------------------------
!
      IF ( first ) THEN
!
! dimension/tracer setting is passed on from phy2chem coupler
!
        print *, 'CKS=>CHEM2PHY_RUN: im, jm, km=', im, jm, km
        print *, 'CKS=>CHEM2PHY_RUN: ntrac =', ntrac
        print *, 'CKS=>CHEM2PHY_RUN: doing_DU =', run_DU
        print *, 'CKS=>CHEM2PHY_RUN: doing_SU =', run_SU
        print *, 'CKS=>CHEM2PHY_RUN: doing_SS =', run_SS
        print *, 'CKS=>CHEM2PHY_RUN: doing_OC =', run_OC
        print *, 'CKS=>CHEM2PHY_RUN: doing_BC =', run_BC
        print *, 'CKS=>CHEM2PHY_RUN: lonr =', lonr
        print *, 'CKS=>CHEM2PHY_RUN: lats_node_r =', lats_node_r
        print *, 'CKS=>CHEM2PHY_RUN: lats_node_r_max =', lats_node_r_max
        print *, 'CKS=>CHEM2PHY_RUN: lonsperlar_r =', lonsperlar_r(:)

!
! determine 2d aer_diag fields from the bundle attribute
!
        name_lst(:) = 'xxxxx'
        kcount      = 0

        MESSAGE_CHECK="CHEM2PHY_RUN: get deltim from phy_exp"
        CALL ESMF_AttributeGet(PHY_EXP_STATE, name = 'deltim',  &
                               value = deltim , rc=RC)
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_CPL)

        lab_setup: DO k = 1, 5
          aerosol = aerosol_list(k)
  
          get_attribute = .False. 
          if (aerosol=='du' .and. run_DU ) get_attribute = .True.
          if (aerosol=='su' .and. run_SU ) get_attribute = .True.
          if (aerosol=='ss' .and. run_SS ) get_attribute = .True.
          if (aerosol=='oc' .and. run_OC ) get_attribute = .True.
          if (aerosol=='bc' .and. run_BC ) get_attribute = .True.

          lab_get_attribute: IF ( get_attribute ) then
          FieldName = trim(aerosol)//'_nfld'
          MESSAGE_CHECK="CHEM2PHY_RUN: get attribute from phy_exp"
          CALL ESMF_AttributeGet(PHY_EXP_STATE, name = FieldName,  &  
                                 value = kcount , rc=RC)
          CALL ERR_MSG(RC,MESSAGE_CHECK,RC_CPL)

          IF ( kcount > 0 ) then

            MESSAGE_CHECK="CHEM2PHY_RUN: get bundle from phy_exp"
            BundleName='dg'//trim(aerosol)
            CALL ESMF_StateGet(PHY_EXP_STATE, BundleName, &
                               Bundle, rc=RC)
            CALL ERR_MSG(RC,MESSAGE_CHECK,RC_CPL)

            do i = 1, kcount
              write(tag, '(i2.2)') i
              FieldName = trim(BundleName)//'_'//trim(tag)
              MESSAGE_CHECK="CHEM2PHY_RUN: get attribute from bundle"
              CALL ESMF_AttributeGet(Bundle, name=FieldName,   &
                                     value=vname, rc=RC)
              CALL ERR_MSG(RC,MESSAGE_CHECK,RC_CPL)
              name_lst(i) = trim(vname)
            enddo
          ENDIF

          select case ( aerosol )
          case ( 'du' )
            nfld_du = kcount
            allocate(name_du(kcount) )
            name_du(1:kcount) = name_lst(1:kcount)
          case ( 'ss' )
            nfld_ss = kcount
            allocate(name_ss(kcount) )
            name_ss(1:kcount) = name_lst(1:kcount)
          case ( 'su' )
            nfld_su = kcount
            allocate(name_su(kcount) )
            name_su(1:kcount) = name_lst(1:kcount)
          case ( 'oc' )
            nfld_oc = kcount
            allocate(name_oc(kcount) )
            name_oc(1:kcount) = name_lst(1:kcount)
          case ( 'bc' )
            nfld_bc = kcount
            allocate(name_bc(kcount) )
            name_bc(1:kcount) = name_lst(1:kcount)
          end select
          ENDIF  lab_get_attribute

        ENDDO  lab_setup

!       reset first flag
        first = .false.

      ENDIF

!---------------------------------------------
!* Get Fortran array from phy export state
!---------------------------------------------

      MESSAGE_CHECK="CHEM2PHY_RUN: Get ItemCount from phy export state"
      call ESMF_StateGet(PHY_EXP_STATE                    &
                        ,itemcount = item_count_phys      &
                        ,rc   =rc)
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_CPL)

      IF ( RC == ESMF_SUCCESS ) THEN

      if (item_count_phys == 0 ) then

       print *, 'CHEM2PHY_RUN: Empty phy export state; Fortran array not filled'

      else

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

       print *, 'CHEM2PHY_RUN: phy export state ItemCount failed,', RC

      ENDIF

!---------------------------------------------
!* Get Fortran array from gocart export state
!---------------------------------------------

      MESSAGE_CHECK="CHEM2PHY_RUN: Get ItemCount from chem export state"
      call ESMF_StateGet(CHEM_EXP_STATE                   &
                        ,itemcount = item_count_chem      &
                        ,rc   =rc)
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_CPL)

      IF ( RC == ESMF_SUCCESS ) THEN

      if (item_count_chem == 0 ) then

       print *, 'CHEM2PHY_RUN: Empty chem export state; Fortran array not filled'

      endif

      ELSE

       print *, 'CHEM2PHY_RUN: chem export state ItemCount failed,', RC

      ENDIF


!---------------------------------------------
!* Do the actual coupling 
!---------------------------------------------
! 
      IF ( RC_CPL == ESMF_SUCCESS ) THEN

       if (item_count_chem==0 .or. item_count_phys==0) then

         print *, 'CHEM2PHY_RUN: Empty state; Coupling phy_exp with chem_exp skipped'

       else

!
! ---  data copy between phy export state and chem export state

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

        if ( run_BC ) then
          do k=1,km/2
            kk = km + 1 - k
            do j=1,jm
              lons_lat = lonsperlar_r(j)
              do i=1,lons_lat
                tem = p_bcphilic(i,j,kk)
                p_bcphilic(i,j,kk) = p_bcphilic(i,j,k)
                p_bcphilic(i,j,k)  = tem

                tem = p_bcphobic(i,j,kk)
                p_bcphobic(i,j,kk) = p_bcphobic(i,j,k)
                p_bcphobic(i,j,k)  = tem
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
                tem = p_ocphilic(i,j,kk)
                p_ocphilic(i,j,kk) = p_ocphilic(i,j,k)
                p_ocphilic(i,j,k)  = tem

                tem = p_ocphobic(i,j,kk)
                p_ocphobic(i,j,kk) = p_ocphobic(i,j,k)
                p_ocphobic(i,j,k) = tem
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

       endif

      ELSE

       print *, 'CHEM2PHY_RUN: Coupling phy_exp with chem_exp skipped'

      ENDIF

!
! --- now let's take care 2d aer_diag fields
!
      lab_copy: DO k = 1, 5
        aerosol = aerosol_list(k)
  
        get_attribute = .False. 
        if (aerosol=='du' .and. run_DU ) get_attribute = .True.
        if (aerosol=='su' .and. run_SU ) get_attribute = .True.
        if (aerosol=='ss' .and. run_SS ) get_attribute = .True.
        if (aerosol=='oc' .and. run_OC ) get_attribute = .True.
        if (aerosol=='bc' .and. run_BC ) get_attribute = .True.

        lab_get_attribute2: IF ( get_attribute ) then

        select case ( aerosol )
          case ( 'du' )
            kcount = nfld_du 
            name_lst(1:kcount) = name_du(1:kcount) 
          case ( 'ss' )
            kcount = nfld_ss
            name_lst(1:kcount) = name_ss(1:kcount) 
          case ( 'su' )
            kcount = nfld_su
            name_lst(1:kcount) = name_su(1:kcount) 
          case ( 'oc' )
            kcount = nfld_oc
            name_lst(1:kcount) = name_oc(1:kcount) 
          case ( 'bc' )
            kcount = nfld_bc
            name_lst(1:kcount) = name_bc(1:kcount) 
        end select

        BundleName='dg'//trim(aerosol)
        do i = 1, kcount

          vname = name_lst(i) 
          nullify(p_diag)
          MESSAGE_CHECK = "Chem2Phys CPL_RUN: Get Farray from Phy_Exp-"//vname
          call GetPointer_diag_(PHY_EXP_STATE, BundleName, vname, p_diag, rc)
          CALL ERR_MSG(RC,MESSAGE_CHECK,RC_CPL)

          nullify(c_diag)
          MESSAGE_CHECK = "Chem2Phys CPL_RUN: Get Farray from Chem_Exp-"//vname
          call GetPointer_diag_(CHEM_EXP_STATE, 'xxxx', vname, c_diag, rc)
!         print *,' inchem_phys, run, be get chem_exp,i=',i,'vname=',trim(vname),'rc=',rc
          CALL ERR_MSG(RC,MESSAGE_CHECK,RC_CPL)

!*        p_diag(:,:) = c_diag(:,:)
          p_diag(:,:) = p_diag(:,:) + c_diag(:,:)*deltim
        enddo    ! kcount-loop

        ENDIF lab_get_attribute2
      ENDDO lab_copy

!
!-----------------------------------------------------------------------
!***  Check the final error signal variable 
!-----------------------------------------------------------------------
!
      IF(RC_CPL==ESMF_SUCCESS)THEN
        WRITE(0,*)'CHEM2PHY CPL RUN SUCCEEDED'
      ELSE
        WRITE(0,*)'CHEM2PHY CPL RUN FAILED RC_CPL=',RC_CPL
      ENDIF
!! 
      END subroutine run

      END module atmos_chem_phy_cpl_comp_mod
