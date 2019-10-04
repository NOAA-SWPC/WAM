!#include "../../../ESMFVersionDefine.h"

!-----------------------------------------------------------------------
!
      module gfs_dyn_phy_cpl_comp_mod
!
!-----------------------------------------------------------------------
!
!***  this module holds the coupler's register, init, run, and finalize 
!***  routines.  they are called from the main gridded component
!***  in module_main_grid_comp.f.
!
!***  the coupler provides 2-way coupling between the dynamics and
!***  physics gridded components by transfering their export and
!***  import states between the two.
!
!! Code Revision:
!! Oct 12 2009        Sarah Lu, atm_cpl_run modified to move Fields and
!!                    FieldBundle between import and export states
!! Oct 16 2009        Sarah Lu, move tracer bundle between states
!! Mar 05 2010        Sarah Lu, modify init routine (associate export state
!!                    to import state) 
!! May 31 2010        Sarah Lu, remove ref to NMM_B and StatePrint call
!! Aug 17 2010        Sarah Lu, debug print only for master PE and the first
!!                    run step
!! Mar 30 2011        Weiyu Yang, modified the code to avoid the ESMF 
!!                    log error.
!! May 12 2011        Weiyu Yang, modified for using the ESMF 5.2.0r_beta_snapshot_07.
!! May 27 2011        Weiyu Yang, modified for using the ESMF 5.2.0r_ library.
!! Jan 13 2016        Jun Wang, fix of the mis-matching of the number of fields for adiabatic
!-----------------------------------------------------------------------
!
      USE ESMF
      use module_export_import_data
      use module_err_msg
!
!
!
!-----------------------------------------------------------------------
!
      implicit none
!     integer 	lm
!
!-----------------------------------------------------------------------
!
      private

      INTEGER, save         :: MYPE             !<-- Each MPI task ID
      logical, save         :: adiabatic
      integer, save         :: n2dflddiff,n3dflddiff
!
      public :: gfs_cpl_setservices
!
!-----------------------------------------------------------------------
!
      contains
!
!-----------------------------------------------------------------------
!#######################################################################
!-----------------------------------------------------------------------
!
      subroutine gfs_cpl_setservices(gc_gfs_cpl,rc_reg)
!
!-----------------------------------------------------------------------
!***  register the coupler component's initialize, run, and finalize
!***  routines.
!-----------------------------------------------------------------------
!
      type(esmf_cplcomp)               :: gc_gfs_cpl 	! coupler component
!
      integer,intent(out) :: rc_reg               	! return code for register
!
!-----------------------------------------------------------------------
!***  local variables
!-----------------------------------------------------------------------
!
      integer :: rc = esmf_success
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
      rc_reg = esmf_success  ! the error signal variable

!-----------------------------------------------------------------------
!***  register the coupler initialize subroutine.  since it is just one
!***  subroutine, no need to specify phase.  the second argument is
!***  a pre-defined subroutine type, such as esmf_setinit, esmf_setrun,
!***  or esmf_setfinal.
!-----------------------------------------------------------------------
!
      call esmf_logwrite("set entry point for coupler initialize"       &
                        ,ESMF_LOGMSG_INFO,rc=rc)
!
      call esmf_cplcompsetentrypoint(gc_gfs_cpl                         &  
                                    ,ESMF_METHOD_INITIALIZE             &  
                                    ,gfs_cpl_initialize                 &  
                                    ,rc=rc)
!
      call err_msg(rc,'set entry point for coupler initialize',rc_reg)
!
!-----------------------------------------------------------------------
!***  register the coupler run subroutine.
!-----------------------------------------------------------------------
!
      call esmf_logwrite("set entry point for coupler run"              &
                        ,ESMF_LOGMSG_INFO,rc=rc)
!
      call esmf_cplcompsetentrypoint(gc_gfs_cpl                         &  
                                    ,ESMF_METHOD_RUN                    &  
                                    ,gfs_cpl_run                        &  
                                    ,rc=rc)
!
      call err_msg(rc,'set entry point for coupler run',rc_reg)
!
!-----------------------------------------------------------------------
!***  register the coupler finalize subroutine.
!-----------------------------------------------------------------------
!
      call esmf_logwrite("set entry point for coupler finalize"         &
                        ,ESMF_LOGMSG_INFO,rc=rc)
!
      call esmf_cplcompsetentrypoint(gc_gfs_cpl                         &  
                                    ,ESMF_METHOD_FINALIZE               &  
                                    ,gfs_cpl_finalize                   &  
                                    ,rc=rc)
!
      call err_msg(rc,'set entry point for coupler finalize',rc_reg)
!
!-----------------------------------------------------------------------
!***  check the error signal variable.
!-----------------------------------------------------------------------
!
      call err_msg(rc,'set coupler services step',rc_reg)
!
!-----------------------------------------------------------------------
!
      end subroutine gfs_cpl_setservices
!
!
!-----------------------------------------------------------------------
!#######################################################################
!-----------------------------------------------------------------------
      subroutine gfs_cpl_initialize(gc_gfs_cpl,imp_state,exp_state      &
                                   ,clock,rc_cpl)
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!***
!***  The init routine was an no-opt place holder before the GOCART plug-in
!***  To enable GOCART initialization, the init routine now associates
!***  fields between import state and export state  (Mar 2010, Sarah Lu)
!***  
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!***  argument variables
!-----------------------------------------------------------------------
!
      type(esmf_cplcomp)               :: gc_gfs_cpl
      type(esmf_state)                 :: imp_state, exp_state
      type(esmf_clock)                 :: clock
!
      integer,           intent(out)   :: rc_cpl
!
!-----------------------------------------------------------------------
!***  local variables
!-----------------------------------------------------------------------
!
      integer                        :: rc,n
      integer                        :: ndata2i,ndata3i,ndata2o,ndata3o
!
      character(esmf_maxstr)         :: import_statename
      character(esmf_maxstr)         :: export_statename
!
      character(20)                  :: array_name
      character(50)                  :: msg
!
      TYPE(ESMF_Field)               :: Field      
      TYPE(ESMF_FieldBundle)         :: Bundle      
!
      TYPE(ESMF_VM)                  :: VM

!-----------------------------------------------------------------------
!***  initialize the error signal variables.
!-----------------------------------------------------------------------
!
      rc = esmf_success
!
!-----------------------------------------------------------------------
!***  determine MPI task ID
!-----------------------------------------------------------------------
      CALL ESMF_CplCompGet(CplComp=gc_gfs_cpl, vm=VM, rc=RC)

      CALL ESMF_VMGet     (vm       = VM               &  !<-- The virtual machine
                          ,localpet = MYPE             &  !<-- Each MPI task ID
                          ,rc = RC)
!
!-----------------------------------------------------------------------
!***  determine the direction of the transfer by extracting
!***  the statename from the import state.
!-----------------------------------------------------------------------
!
      call esmf_logwrite("retrieve state name in coupler"               &
                        ,ESMF_LOGMSG_INFO,rc=rc)
!
      msg = 'retrieve state name from cpl import state'
      call esmf_stateget(imp_state, name=import_statename, rc=rc)
      call err_msg(rc,msg,rc_cpl)

      msg = 'retrieve adiabatic from cpl import state'
      call ESMF_AttributeGet(imp_state,"adiabatic",adiabatic,rc=rc)
      call err_msg(rc,msg,rc_cpl)
      n2dflddiff = max(ndata_2d_phy_exp-ndata_2d_phy_imp,0)
      n3dflddiff = max(ndata_3d_phy_exp-ndata_3d_phy_imp,0)

      msg = 'retrieve state name from cpl export state'
      call esmf_stateget(exp_state, name=export_statename, rc=rc)
      call err_msg(rc,msg,rc_cpl)
!
      call err_msg(rc,msg,rc_cpl)
      if ( MYPE == 0) then
         print *,'GFS_CPL INIT is to associate data from '              &
             ,' (',trim(import_statename),') with '                     &
             ,' (',trim(export_statename),') ,adiabatic=',adiabatic
      endif
!
!-----------------------------------------------------------------------
!***  the number of fields transferred from the dynamics to
!***  the physics may not equal the number of fields transferred
!***  in the other direction.  these values are specified in
!***  module_export_import_data.
!-----------------------------------------------------------------------
!
! get from dynamics import state
      if(trim(import_statename) == "GFS dynamics import") then
        ndata3i = ndata_3d_dyn_imp
        ndata2i = ndata_2d_dyn_imp
        if(adiabatic) then
          ndata2i = ndata_2d_dyn_imp - n2dflddiff
          ndata3i = ndata_3d_dyn_imp - n3dflddiff
        endif
!
! get from dynamics export state
      elseif(trim(import_statename) == "GFS dynamics export") then
        ndata3i = ndata_3d_dyn_exp
        ndata2i = ndata_2d_dyn_exp
!
! get from physics import state
      elseif(trim(import_statename) == "physics import") then
        ndata3i = ndata_3d_phy_imp
        ndata2i = ndata_2d_phy_imp
!
! get from physics export state
      elseif(trim(import_statename) == "physics export") then
        ndata3i = ndata_3d_phy_exp
        ndata2i = ndata_2d_phy_exp
        if(adiabatic) then
          ndata2i = ndata_2d_phy_exp - n2dflddiff
          ndata3i = ndata_3d_phy_exp - n3dflddiff
        endif
!
      else
        print *,' Error: no state name match, state_name='         &
               , trim(import_statename)
      endif
      if ( MYPE == 0) then
        print *,'GFS_CPL INIT import state is ',trim(import_statename)
        print *,'GFS_CPL INIT ndata2i ndata3i are ',ndata2i,ndata3i
      endif
!
! ---------------------------------------------------------------------
! put to dynamics import state
      if(trim(export_statename) == "GFS dynamics import") then
        ndata3o = ndata_3d_dyn_imp
        ndata2o = ndata_2d_dyn_imp
        if(adiabatic) then
          ndata2o = ndata_2d_dyn_imp - n2dflddiff
          ndata3o = ndata_3d_dyn_imp - n3dflddiff
        endif
!
! put to dynamics export state
      elseif(trim(export_statename) == "GFS dynamics export") then
        ndata3o = ndata_3d_dyn_exp
        ndata2o = ndata_2d_dyn_exp
!
! put to physics import state
      elseif(trim(export_statename) == "physics import") then
        ndata3o = ndata_3d_phy_imp
        ndata2o = ndata_2d_phy_imp
!
! put to physics export state
      elseif(trim(export_statename) == "physics export") then
        ndata3o = ndata_3d_phy_exp
        ndata2o = ndata_2d_phy_exp
        if(adiabatic) then
          ndata2o = ndata_2d_phy_exp - n2dflddiff
          ndata3o = ndata_3d_phy_exp - n3dflddiff
        endif
!
      else
        print *,' Error: no state name match, state_name='         &
               , trim(export_statename)
      endif
      if ( MYPE == 0) then
        print *,'GFS_CPL INIT export state is ',trim(export_statename)
        print *,'GFS_CPL INIT ndata2o ndata3o are ',ndata2o,ndata3o,'adiabatic=',adiabatic
      endif
!
! --  check item member
      if ( ndata2o > ndata2i .or. ndata3o > ndata3i ) then
!       print *,' ERROR: import data is too few for export data '
        print *,' ndata2o=',ndata2o,' ndata2i=',ndata2i,           &
                ' ndata3o=',ndata3o,' ndata3i=',ndata3i
!       call abort
      endif
!
!-----------------------------------------------------------------------
!***  loop through the field data names, extract those fields from the
!***  import state, and add them to the export state.
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!***  do 2-d arrays.
!-----------------------------------------------------------------------
!
      data_2d: do n=1,ndata2o

        array_name=trim(datanames_2d(n))
        if(MYPE==0)print *, 'GFS_CPL INIT transfer ', array_name
!
        msg = 'retrieve field '//array_name//' from cpl import'
        call ESMF_StateGet(imp_state, array_name, Field, rc=rc)
        call err_msg(rc,msg,rc_cpl)    

        msg = '1, add field to cpl export'
        call ESMF_StateAddReplace(exp_state, (/Field/), rc=rc)
        call err_msg(rc,msg,rc_cpl)    
!
      enddo data_2d
!
!-----------------------------------------------------------------------
!***  do 3-d arrays.
!-----------------------------------------------------------------------
!
      data_3d: do n=1,ndata3o

        array_name=trim(datanames_3d(n))
        if(MYPE==0)print *, 'GFS_CPL INIT transfer ', array_name
!
        msg = 'retrieve field '//array_name//' from cpl import'
        call ESMF_StateGet(imp_state, array_name, Field, rc=rc)
        call err_msg(rc,msg,rc_cpl)    

        msg = '2, add field to cpl export'
        call ESMF_StateAddReplace(exp_state, (/Field/), rc=rc)
        call err_msg(rc,msg,rc_cpl)    
!
      enddo data_3d

!-----------------------------------------------------------------------
!***  do tracer arrays.
!-----------------------------------------------------------------------
!
      if(MYPE==0)print *, 'GFS_CPL INIT transfer tracer bundle'
      msg = 'retrieve tracer bundle from cpl import'
      call ESMF_StateGet(imp_state, 'tracers', Bundle, rc=rc)
      call err_msg(rc,msg,rc_cpl)    

      msg = 'add bundle to cpl export'
      call ESMF_StateAddReplace(exp_state, (/Bundle/), rc=rc)
      call err_msg(rc,msg,rc_cpl)    

      end subroutine gfs_cpl_initialize
!
!-----------------------------------------------------------------------
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!-----------------------------------------------------------------------
!
      subroutine gfs_cpl_run(gc_gfs_cpl, imp_state, exp_state, clock, rc_cpl)
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!***  run the coupler to transfer data between the gridded components.
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!***  argument variables
!-----------------------------------------------------------------------
!
      type(esmf_cplcomp)               :: gc_gfs_cpl
      type(esmf_state)                 :: imp_state, exp_state
      type(esmf_clock)                 :: clock
!
      integer,           intent(out)   :: rc_cpl
!
!-----------------------------------------------------------------------
!***  local variables
!-----------------------------------------------------------------------
!
      integer :: rc, rcfinal, l, n, imp_item, exp_item
      integer :: ndata1i,ndata2i,ndata3i,ndata1o,ndata2o,ndata3o
!
      character(esmf_maxstr) :: import_statename, export_statename
!
      type(esmf_array) :: hold_array
!     type(esmf_routehandle) :: routehandle
!
      character(3)  :: model_level
      character(6)  :: format
      character(20) :: array_name
      character(20) :: imp_item_name(100), exp_item_name(100)
!
! Additional variables are added, allowing the couple to handle array,
! field, and/or field bundle  (Sarah Lu)
!
      TYPE(ESMF_Field)               :: ESMFField        
      TYPE(ESMF_FieldBundle)         :: ESMFBundle      
      integer                        :: rc2           
      logical, save :: cpl_array  = .false., cpl_field  = .false.    &
                     , cpl_bundle = .false., first      = .true.     &
                     , from_exp_dyn_to_imp_phy = .false.            &
                     , from_exp_phy_to_imp_dyn = .false.            &
                     , from_imp_dyn_to_exp_dyn = .false.            &
                     , from_imp_phy_to_exp_phy = .false.

!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!***  initialize the error signal variables.
!-----------------------------------------------------------------------
!
      rc      = esmf_success
      rcfinal = esmf_success
!
!-----------------------------------------------------------------------
!***  determine the direction of the transfer by extracting
!***  the statename from the import state.
!-----------------------------------------------------------------------
!
      call esmf_logwrite("retrieve state name in coupler"               &
                        ,ESMF_LOGMSG_INFO,rc=rc)
!
      call esmf_stateget(imp_state, name=import_statename, rc=rc)
      call esmf_stateget(exp_state, name=export_statename, rc=rc)
!
      IF ( first .and. (MYPE == 0)  ) then
           print *,' coupler is to move data from '                     &
             ,' (',trim(import_statename),') to '                       &
             ,' (',trim(export_statename),') '
      ENDIF
!
! if the pointer has already set, return
!
      if( trim(import_statename) == 'GFS dynamics export' .and.         &
          trim(export_statename) == 'physics import'  .and.             &
          from_exp_dyn_to_imp_phy ) return

      if( trim(import_statename) == 'physics export'  .and.             &
          trim(export_statename) == 'GFS dynamics import' .and.         &
          from_exp_phy_to_imp_dyn ) return

      if( trim(import_statename) == 'GFS dynamics import' .and.         &
          trim(export_statename) == 'GFS dynamics export' .and.         &
          from_imp_dyn_to_exp_dyn ) return

      if( trim(import_statename) == 'physics import' .and.              &
          trim(export_statename) == 'physics export' .and.              &
          from_imp_phy_to_exp_phy ) return

      call esmf_stateget(imp_state                                      &
                        ,itemcount = imp_item                           &
                        ,itemnamelist = imp_item_name                   &
                        ,rc   =rc)
      call esmf_stateget(exp_state                                      &
                        ,itemcount = exp_item                           &
                        ,itemnamelist = exp_item_name                   &
                        ,rc   =rc)
!
!      print *,' import item count is ',imp_item
!      print *,' export item count is ',exp_item
!
!      print *,' import item name ',(imp_item_name(n),n=1,imp_item)
!      print *,' export item name ',(exp_item_name(n),n=1,exp_item)
!
      call err_msg(rc,'retrieve state name in coupler',rcfinal)
!
      rcfinal = rc

!
!-----------------------------------------------------------------------
! check the data class in the imp/exp state (Sarah Lu)

      if ( first ) then                                               
      array_name = trim(datanames_2d(1))                               
!---> query hs as an ESMF Field
      CALL ESMF_StateGet(state     = imp_state                      & 
                        ,itemName  = array_name                     &
                        ,field     = ESMFField                      & 
                        ,rc        = rc2)                            
      if ( rc2 == esmf_success ) then                               
         if(MYPE == 0) print *,' LU_CPL: ESMFField found in imp state'
         cpl_field = .true.                                     
      else                                                     
         if(MYPE == 0) print *,' LU_CPL: ESMFField not found in imp state'
         cpl_field = .false.                                 
      endif                                           
!---> query tracer bundle
      CALL ESMF_StateGet(state       = imp_state                    & 
                        ,itemName    = 'tracers'                    & 
                        ,fieldbundle = ESMFBundle                   & 
                        ,rc          = rc2)                          
      if ( rc2 == esmf_success ) then                          
         if(MYPE == 0) print *,' LU_CPL: ESMFBundle found in imp state'
         cpl_bundle = .true.                                       
      else                                                        
         if(MYPE == 0) print *,' LU_CPL: ESMFBundle not found in imp state'
         cpl_bundle = .false.                                  
      endif                                              

!     first = .false.              
      endif    
!
!-----------------------------------------------------------------------
!***  the number of fields transferred from the dynamics to
!***  the physics may not equal the number of fields transferred
!***  in the other direction.  these values are specified in
!***  module_export_import_data.
!-----------------------------------------------------------------------
!
! get from dynamics import state
      if(trim(import_statename) == "GFS dynamics import") then
        ndata3i = ndata_3d_dyn_imp
        ndata2i = ndata_2d_dyn_imp
        ndata1i = ndata_1d_dyn_imp
        if(adiabatic) then
          ndata2i = ndata_2d_dyn_imp - n2dflddiff
          ndata3i = ndata_3d_dyn_imp - n3dflddiff
        endif
!       write(0,*)' import state is from dynamics_import_state '
!       write(0,*)' ndata1i ndata2i ndata3i are ',ndata1i,ndata2i,ndata3i
!
! get from dynamics export state
      elseif(trim(import_statename) == "GFS dynamics export") then
        ndata3i = ndata_3d_dyn_exp
        ndata2i = ndata_2d_dyn_exp
        ndata1i = ndata_1d_dyn_exp
!       write(0,*)' import state is from dynamics_export_state '
!       write(0,*)' ndata1i ndata2i ndata3i are ',ndata1i,ndata2i,ndata3i
!
! get from physics import state
      elseif(trim(import_statename) == "physics import") then
        ndata3i = ndata_3d_phy_imp
        ndata2i = ndata_2d_phy_imp
        ndata1i = ndata_1d_phy_imp
!       write(0,*)' import state is from physics_import_state '
!       write(0,*)' ndata1i ndata2i ndata3i are ',ndata1i,ndata2i,ndata3i
!
! get from physics export state
      elseif(trim(import_statename) == "physics export") then
        ndata3i = ndata_3d_phy_exp
        ndata2i = ndata_2d_phy_exp
        ndata1i = ndata_1d_phy_exp

!       write(0,*)' import state is from physics_export_state '
!       write(0,*)' ndata1i ndata2i ndata3i are ',ndata1i,ndata2i,ndata3i
        if(adiabatic) then
          ndata2i = ndata_2d_phy_exp - n2dflddiff
          ndata3i = ndata_3d_phy_exp - n3dflddiff
        endif!
      else
        print *,' Error: no state name match, state_name='         &
               , trim(import_statename)
      endif
!
! ---------------------------------------------------------------------
! put to dynamics import state
      if(trim(export_statename) == "GFS dynamics import") then
        ndata3o = ndata_3d_dyn_imp
        ndata2o = ndata_2d_dyn_imp
        ndata1o = ndata_1d_dyn_imp
        if(adiabatic) then
          ndata2o = ndata_2d_dyn_imp - n2dflddiff
          ndata3o = ndata_3d_dyn_imp - n3dflddiff
        endif
!       write(0,*)' export state is for dynamics_import_state '
!       write(0,*)' ndata1o ndata2o ndata3o are ',ndata1o,ndata2o,ndata3o
!
! put to dynamics export state
      elseif(trim(export_statename) == "GFS dynamics export") then
        ndata3o = ndata_3d_dyn_exp
        ndata2o = ndata_2d_dyn_exp
        ndata1o = ndata_1d_dyn_exp
!       write(0,*)' export state is for dynamics_export_state '
!       write(0,*)' ndata1o ndata2o ndata3o are ',ndata1o,ndata2o,ndata3o
!
! put to physics import state
      elseif(trim(export_statename) == "physics import") then
        ndata3o = ndata_3d_phy_imp
        ndata2o = ndata_2d_phy_imp
        ndata1o = ndata_1d_phy_imp
!       write(0,*)' export state is for physics_import_state '
!       write(0,*)' ndata1o ndata2o ndata3o are ',ndata1o,ndata2o,ndata3o
!
! put to physics export state
      elseif(trim(export_statename) == "physics export") then
        ndata3o = ndata_3d_phy_exp
        ndata2o = ndata_2d_phy_exp
        ndata1o = ndata_1d_phy_exp
        if(adiabatic) then
          ndata2o = ndata_2d_phy_exp - n2dflddiff
          ndata3o = ndata_3d_phy_exp - n3dflddiff
        endif
!       write(0,*)' export state is for physics_export_state '
!       write(0,*)' ndata1o ndata2o ndata3o are ',ndata1o,ndata2o,ndata3o
!
      else
        print *,' Error: no state name match, state_name='         &
               , trim(export_statename)
      endif
!
!-----------------------------------------------------------------------
!***  loop through the field data names, extract those fields from the
!***  import state, and add them to the export state.
!-----------------------------------------------------------------------
!***  loop through the field data names, extract those fields from the
!***  import state, and add them to the export state.
!-----------------------------------------------------------------------
!
      if( ndata1o > ndata1i .or.                                       &
          ndata2o > ndata2i .or.                                       &
          ndata3o > ndata3i ) then
!       print *,' ERROR: import data is too few for export data '
        write(0,*)' ndata1o=',ndata1o,' ndata1i=',ndata1i,           &
                  ' ndata2o=',ndata2o,' ndata2i=',ndata2i,           &
                  ' ndata3o=',ndata3o,' ndata3i=',ndata3i
!       call abort
      endif
!
!-----------------------------------------------------------------------
!***  do the 1-d arrays.
!-----------------------------------------------------------------------

!    
  lab_if_1d : if ( cpl_array ) then
      data_1d: do n=1,ndata1o
!
!-----------------------------------------------------------------------
!
        array_name = trim(datanames_1d(n))
!
!       print *,' get ',array_name
!
        CALL ESMF_StateGet(imp_state, array_name, hold_array, rc = rc)
!
!       call err_msg(rc,'retrieve array from cpl import',rcfinal)
!
!-----------------------------------------------------------------------
!
!       print *,' add ',array_name
!
        CALL ESMF_StateAddReplace(exp_state, (/hold_array/), rc = rc)
!
!       call err_msg(rc,'add array to cpl export',rcfinal)
!
!-----------------------------------------------------------------------
!
      enddo data_1d
  endif  lab_if_1d
!
!-----------------------------------------------------------------------
!***  now do the 2-d arrays.
!-----------------------------------------------------------------------
!
      data_2d: do n=1,ndata2o
!

!-----------------------------------------------------------------------
!
        array_name = trim(datanames_2d(n))
!
!       print *,' get ',array_name
!
        lab_if_2darray : if ( cpl_array ) then     ! loop through 2d array

        CALL ESMF_StateGet(imp_state, array_name, hold_array, rc = rc)

        CALL ESMF_StateAddReplace(exp_state, (/hold_array/), rc = rc)
!
!-----------------------------------------------------------------------
!
        endif lab_if_2darray

        lab_if_2dfield : if ( cpl_field ) then     ! loop through field

          call ESMF_StateGet(state     = imp_state                   & 
                           ,itemName   = array_name                  & 
                           ,field      = ESMFField                   & 
                           ,rc   =rc)                                 
          call err_msg(rc,'retrieve field from cpl import',rcfinal)    

          call ESMF_StateAddReplace(exp_state, (/ESMFField/),rc   =rc)
          call err_msg(rc,'6, add field to cpl export',rcfinal)         

        endif lab_if_2dfield                                       

      enddo data_2d
!
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!***  do 3-d arrays.
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!
        data_3d: do n=1,ndata3o
!
!-----------------------------------------------------------------------
!
          array_name = trim(datanames_3d(n))
!
!         print *,' get ',array_name
!

         lab_if_3darray : if ( cpl_array ) then   ! loop through 3d array

          CALL ESMF_StateGet(imp_state, array_name, hold_array, rc = rc)
!
          CALL ESMF_StateAddReplace(exp_state, (/hold_array/), rc = rc)
!
         endif lab_if_3darray   

         lab_if_3dfield : if ( cpl_field ) then    ! loop through 3d fields
          call ESMF_StateGet(state      = imp_state                    &
                             ,itemName  = array_name                   &
                             ,field     = ESMFField                    &
                             ,rc   =rc)                                 
          call err_msg(rc,'retrieve field from cpl import',rcfinal)   

          call ESMF_StateAddReplace(exp_state, (/ESMFField/),rc = rc)
          call err_msg(rc,'3, add field to cpl export',rcfinal)        
!
         endif lab_if_3dfield                               
!-----------------------------------------------------------------------
!
        enddo data_3d

!
!-----------------------------------------------------------------------
!
!!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!! Extract tracer field bundle from the import state and add them
!! to the export state  (Sarah Lu)
!! 
         lab_if_bundle : if ( cpl_bundle ) then

          call ESMF_StateGet(state      = imp_state                    &
                             ,itemName  = 'tracers'                    &
                             ,fieldbundle = ESMFBundle                 &
                             ,rc   =rc)                                
          call err_msg(rc,'retrieve bundle from cpl import',rcfinal)

          call ESMF_StateAddReplace(exp_state, (/ESMFBundle/),rc = rc)
          call err_msg(rc,'add bundle to cpl export',rcfinal)         

         endif lab_if_bundle                                   
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
!
! check the export state  (optional)
!

      IF ( first ) THEN
        first = .false.                         !<-- reset first flag

        if ( MYPE==0 ) THEN
        call esmf_stateget(exp_state                                    &
                        ,name =export_statename                         &
                        ,itemcount = exp_item                           &
                        ,itemnamelist = exp_item_name                   &
                        ,rc   =rc)

        print *,' coupler is done for expor state '                     &
             ,' (',trim(export_statename),') with item =',exp_item      &
             ,' and item name is ',(exp_item_name(n),n=1,exp_item)
!
        rcfinal=rc
        ENDIF

      ENDIF
!
! make sure to run once
!
!      if( trim(import_statename).eq.'GFS dynamics export' .and.          &
!          trim(export_statename).eq.'physics import'  )                  &
!          from_exp_dyn_to_imp_phy = .true.
!
!      if( trim(import_statename).eq.'physics export'  .and.             &
!          trim(export_statename).eq.'GFS dynamics import' )             &
!          from_exp_phy_to_imp_dyn = .true.
!
!      if( trim(import_statename).eq.'GFS dynamics import' .and.         &
!          trim(export_statename).eq.'GFS dynamics export' )             &
!          from_imp_dyn_to_exp_dyn = .true.
!
!      if( trim(import_statename).eq.'physics import' .and.              &
!          trim(export_statename).eq.'physics export' )                  &
!          from_imp_phy_to_exp_phy = .true.

!-----------------------------------------------------------------------
!
      call err_msg(rc,'main cpl run step',rcfinal)
!
      rc_cpl=rcfinal
!
!-----------------------------------------------------------------------
!
      end subroutine gfs_cpl_run
!
!-----------------------------------------------------------------------
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!-----------------------------------------------------------------------
!
      subroutine gfs_cpl_finalize(gc_gfs_cpl,imp_state,exp_state, clock,rc_cpl)
!
!-----------------------------------------------------------------------
!***  finalize the coupler.
!-----------------------------------------------------------------------
!
!
      type(esmf_cplcomp)               :: gc_gfs_cpl
      type(esmf_state)                 :: imp_state, exp_state
      type(esmf_clock)                 :: clock
!
      integer,           intent(out)   :: rc_cpl
!      
!-----------------------------------------------------------------------
!***  local variables
!-----------------------------------------------------------------------
!
      integer :: rc,rcfinal
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
      rc=esmf_success
      rcfinal=esmf_success
!
!-----------------------------------------------------------------------
!
      call err_msg(rc,'main cpl finalize step',rcfinal)
!
      rc_cpl = rcfinal
!
!-----------------------------------------------------------------------
!
      end subroutine gfs_cpl_finalize
!
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!
      end module gfs_dyn_phy_cpl_comp_mod
!
!-----------------------------------------------------------------------
