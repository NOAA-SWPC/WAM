!#include "../../../ESMFVersionDefine.h"

!-----------------------------------------------------------------------
!
      MODULE gfs_dynamics_output
!
!-----------------------------------------------------------------------
!***  LIST THE QUANTITIES FROM THE DYNAMICS INTERNAL STATE THAT
!***  CAN BE SELECTED FOR HISTORY OUTPUT AND ASSOCIATE THEM WITH
!***  UNIQUE INTEGERS.
!***  THE USER WILL PLACE AN 'H' IN THE 3rd ELEMENT OF THE FOLLOWING
!***  LIST OF QUANTITIES THAT ARE IN THE DYNAMICS INTERNAL STATE
!***  IF THAT QUANTITY IS TO BE WRITTEN TO THE HISTORY FILES.
!-----------------------------------------------------------------------
!***  WHEN NEW QUANTITIES ARE ADDED TO THE DYNAMICS INTERNAL STATE
!***  THAT ARE CANDIDATES FOR HISTORY OUTPUT THEN THEY NEED TO BE
!***  ADDED IN TWO PLACES BELOW: 
!*    (1) THE APPROPRIATE DATA LIST PRECEDING THE 'CONTAINS' STATEMENT
!*    (2) 'THE DYNAMICS INTERNAL STATE POINTER BLOCK'
!*        IN SUBROUTINE POINT_DYNAMICS_OUPUT
!
!***  REVISION LOG:
!***  Nov 23 2009, Sarah Lu   Add gocart species to DYN_INT_STATE_3D_R_DIAB
!***  Feb 20 2011, Henry Juang  use gfs_dyn_write_state for more related name
!***  February 2011 Weiyu Yang, Updated to use both the ESMF 4.0.0rp2 library,
!***                            ESMF 5 library and the the ESMF 3.1.0rp2 library.
!***  Aug 2011      Jun Wang   remove hardcode of adiabetic tracer when setting
!***                           diabetic tracer names
!***  Jan 2016      S Moorthi  Added option not to output pres and dpres for semilag
!***                           and some cleanup
!***
!-----------------------------------------------------------------------
!
      USE ESMF
      USE gfs_dyn_machine
      USE gfs_dynamics_internal_state_mod,ONLY: gfs_dynamics_internal_state 
      use namelist_dynamics_def,          only: semilag
      use gfs_dynamics_err_msg_mod
      use gfs_dyn_tracer_config,          ONLY: gfs_dyn_tracer      ! generalized tracer
!
      IMPLICIT NONE
!
      PRIVATE
!
      PUBLIC :: POINT_DYNAMICS_OUTPUT_GFS,                              &
                DYN_INT_STATE_ISCALAR,    DYN_INT_STATE_RSCALAR,        &
                DYN_INT_STATE_LSCALAR,    DYN_INT_STATE_1D_I,           &
                DYN_INT_STATE_2D_I,       DYN_INT_STATE_1D_R,           &
                DYN_INT_STATE_2D_R,       DYN_INT_STATE_3D_R_DIAB,      &
                DYN_INT_STATE_3D_R_ADIAB, DYN_INT_STATE_4D_R
!
!-----------------------------------------------------------------------
!
      INTEGER,PARAMETER :: MAX_KOUNT=100
!
!-----------------------------------------------------------------------
!***  LIST THE VARIOUS QUANTITIES IN THE DYNAMICS INTERNAL STATE 
!***  THAT ARE CANDIDATES FOR HISTORY OUTPUT.
!***  GROUP THEM BY TYPE (VS_OUTPUTARIOUS SCALARS, VARIOUS ARRAYS) SINCE
!***  WE WILL USE A WORKING POINTER FOR EACH TYPE.
!-----------------------------------------------------------------------
!
!-------------------------
!***  INTEGER SCALARS  ***
!-------------------------
!
      CHARACTER(13),DIMENSION(3,MAX_KOUNT) :: DYN_INT_STATE_ISCALAR   &
!
       =RESHAPE((/ 'latg         ', 'OGFS_SIG     ', 'R            '  &
                  ,'lonf         ', 'OGFS_SIG     ', 'R            '  &
                  ,'levs         ', 'OGFS_SIG     ', 'R            '  &
                  ,'jcap         ', 'OGFS_SIG     ', 'R            '  &
                  ,'ntoz         ', 'OGFS_SIG     ', 'R            '  &
                  ,'ntcw         ', 'OGFS_SIG     ', 'R            '  &
                  ,'ncld         ', 'OGFS_SIG     ', 'R            '  &
                  ,'ntke         ', 'OGFS_SIG     ', 'R            '  &
                  ,'ntrac        ', 'OGFS_SIG     ', 'R            '  &
                  ,'vertcoord_id ', 'OGFS_SIG     ', 'R            '  &
                  ,'thermodyn_id ', 'OGFS_SIG     ', 'R            '  &
                  ,'sfcpress_id  ', 'OGFS_SIG     ', 'R            '  &
                  ,'ienst        ', 'OGFS_SIG     ', 'R            '  &
                  ,'iensi        ', 'OGFS_SIG     ', 'R            '  &
                  ,'itrun        ', 'OGFS_SIG     ', 'R            '  &
                  ,'icen2        ', 'OGFS_SIG     ', 'R            '  &
                 /)                                                   &
               ,(/3,MAX_KOUNT/)                                       &
               ,(/'*            ', '*            ', '*            '/))
!
!----------------------
!***  REAL SCALARS  ***
!----------------------
!
      CHARACTER(15),DIMENSION(3,MAX_KOUNT) :: DYN_INT_STATE_RSCALAR       &
!
       =RESHAPE((/ 'pdryini        ', 'OGFS_SIG       ', 'R              '&
                 /)                                                       &
               ,(/3,MAX_KOUNT/)                                           &
               ,(/'*              ', '*              ', '*              '/))
!
!----------------------
!***  LOGICAL SCALARS  ***
!----------------------
!
      CHARACTER(16),DIMENSION(3,MAX_KOUNT) :: DYN_INT_STATE_LSCALAR           &
!
       =RESHAPE((/ 'HYBRID          ', 'OGFS_SIG        ', 'R               ' &
                  ,'GEN_COORD_HYBRID', 'OGFS_SIG        ', 'R               ' &
                  ,'ADIABATIC       ', 'OGFS_SIG        ', 'R               ' &
                 /)                                                           &
               ,(/3,MAX_KOUNT/)                                               &
               ,(/'*               ', '*               ', '*               '/))
!
!----------------------------
!***  INTEGER 1-D ARRAYS  ***
!----------------------------
!
      CHARACTER(10),DIMENSION(3,MAX_KOUNT) :: DYN_INT_STATE_1D_I        &
!
       =RESHAPE((/ 'IDAT      ', 'OGFS_SIG  ', 'R         '             &
                 /)                                                     &
               ,(/3,MAX_KOUNT/)                                         &
               ,(/'**********', '**********', '**********'/))
!
!----------------------------
!***  INTEGER 2-D ARRAYS  ***
!----------------------------
!
      CHARACTER(10),DIMENSION(3,MAX_KOUNT) :: DYN_INT_STATE_2D_I        &
!
       =RESHAPE((/ '-         ', '-         ', '-         '             &
                 /)                                                     &
               ,(/3,MAX_KOUNT/)                                         &
               ,(/'**********', '**********', '**********'/))
!
!-------------------------
!***  REAL 1-D ARRAYS  *** 
!-------------------------
!
      CHARACTER(13),DIMENSION(3,MAX_KOUNT) :: DYN_INT_STATE_1D_R        &
!
       =RESHAPE((/ 'AK5          ', 'OGFS_SIG     ', 'R            '    &
                  ,'BK5          ', 'OGFS_SIG     ', 'R            '    &
                  ,'CK5          ', 'OGFS_SIG     ', 'R            '    &
                  ,'SI           ', 'OGFS_SIG     ', 'R            '    &
                  ,'CPI          ', 'OGFS_SIG     ', 'R            '    &
                  ,'RI           ', 'OGFS_SIG     ', 'R            '    &
                 /)                                                     &
               ,(/3,MAX_KOUNT/)                                         &
               ,(/'*            ', '*            ', '*            '/))
!
!-------------------------
!***  REAL 2-D ARRAYS  ***
!
      CHARACTER(10),DIMENSION(3,MAX_KOUNT) :: DYN_INT_STATE_2D_R        &
!
       =RESHAPE((/ 'hgt       ', 'OGFS_SIG  ', 'R         '             &
                  ,'pres      ', 'OGFS_SIG  ', 'R         '             &
                 /)                                                     &
               ,(/3,MAX_KOUNT/)                                         &
               ,(/'*         ', '*         ', '*         '/))
!
!-------------------------
!***  REAL 3-D ARRAYS  ***
!-------------------------
!
!------------------------------
!***  Diabatic forecast output
!------------------------------
!
!  generalized tracer: 
!  DYN_INT_STATE_3D_R_DIAB is determined in run-time
!
      CHARACTER(10),DIMENSION(3,MAX_KOUNT),TARGET ::                    &
!
                                             DYN_INT_STATE_3D_R_DIAB
!                                            DYN_INT_STATE_3D_R_DIAB     &  
!!
!       =RESHAPE((/ 'dpres     ', 'OGFS_SIG  ', 'levs      ' &  !<-- The physics counterparts of these variables
!                  ,'pres      ', 'OGFS_SIG  ', 'levs      ' &  !<-- The physics counterparts of these variables
!                  ,'ugrd      ', 'OGFS_SIG  ', 'levs      ' &  !    are being designated for history output.
!                  ,'vgrd      ', 'OGFS_SIG  ', 'levs      ' &  !    This assumes that history output always
!                  ,'tmp       ', 'OGFS_SIG  ', 'levs      ' &  !    immediately follows a call to the Physics.
!                  ,'spfh      ', 'OGFS_SIG  ', 'levs      ' &  !
!                  ,'o3mr      ', 'OGFS_SIG  ', 'levs      ' &  !<--
!                  ,'clwmr     ', 'OGFS_SIG  ', 'levs      ' &
!                  ,'du001     ', 'OGFS_SIG  ', 'levs      ' &  !<-- gocart species
!                  ,'du002     ', 'OGFS_SIG  ', 'levs      ' &
!                  ,'du003     ', 'OGFS_SIG  ', 'levs      ' &
!                  ,'du004     ', 'OGFS_SIG  ', 'levs      ' &
!                  ,'du005     ', 'OGFS_SIG  ', 'levs      ' &
!
!                 /)                                         &
!               ,(/3,MAX_KOUNT/)                             &
!               ,(/'**********', '**********', '**********'/))

!
!-------------------------------
!***  Adiabatic forecast output
!-------------------------------
!
      CHARACTER(10),DIMENSION(3,MAX_KOUNT),TARGET ::        &
                                DYN_INT_STATE_3D_R_ADIAB    &  
!
       =RESHAPE((/ 'dpres     ', 'OGFS_SIG  ', 'levs      ' &  !<-- The physics counterparts of these variables
                  ,'pres      ', 'OGFS_SIG  ', 'levs      ' &  !<-- The physics counterparts of these variables
                  ,'ugrd      ', 'OGFS_SIG  ', 'levs      ' &  !    are being designated for history output.
                  ,'vgrd      ', 'OGFS_SIG  ', 'levs      ' &  !    This assumes that history output always
                  ,'tmp       ', 'OGFS_SIG  ', 'levs      ' &  !    immediately follows a call to the Physics.
                  ,'spfh      ', 'OGFS_SIG  ', 'levs      ' &  !
                  ,'o3mr      ', 'OGFS_SIG  ', 'levs      ' &  !<--
                  ,'clwmr     ', 'OGFS_SIG  ', 'levs      ' &
                  ,'tke       ', 'OGFS_SIG  ', 'levs      ' &
                 /)                                         &
               ,(/3,MAX_KOUNT/)                             &
               ,(/'**********', '**********', '**********'/))
!
!-------------------------
!***  REAL 4-D ARRAYS  ***
!-------------------------
!
      CHARACTER(12),DIMENSION(3,MAX_KOUNT) :: DYN_INT_STATE_4D_R        &
!
       =RESHAPE((/ 'TRACERS     ', 'OGFS_SIG    ', 'levs        '       &
                 /)                                                     &
               ,(/3,MAX_KOUNT/)                                         &
               ,(/'*           ', '*           ', '*           '/))
!
!-----------------------------------------------------------------------
!***  FIRST WE MUST PROVIDE POINTERS INTO THE DYNAMICS INTERNAL STATE.
!-----------------------------------------------------------------------
      TYPE DYN_ISC
        INTEGER(KIND=kind_io4),               POINTER :: NAME  !<-- Pointer for integer scalars
      END TYPE DYN_ISC
!-----------------------------------------------------------------------
      TYPE DYN_RSC
        REAL(KIND=kind_grid),                 POINTER :: NAME  !<-- Pointer for real scalars
      END TYPE DYN_RSC
!-----------------------------------------------------------------------
      TYPE DYN_LSC
        LOGICAL,                              POINTER :: NAME  !<-- Pointer for real scalars
      END TYPE DYN_LSC
!-----------------------------------------------------------------------
      TYPE DYN_I1D
        INTEGER(KIND=kind_io4),DIMENSION(:),  POINTER :: NAME  !<-- Pointer for 1D integer arrays
      END TYPE DYN_I1D
!-----------------------------------------------------------------------
      TYPE DYN_I2D
        INTEGER(KIND=kind_io4),DIMENSION(:,:),POINTER :: NAME   !<-- Pointer for 2D integer arrays
      END TYPE DYN_I2D
!-----------------------------------------------------------------------
      TYPE DYN_R1D
        REAL(KIND=kind_evod),DIMENSION(:),    POINTER :: NAME   !<-- Pointer for 1D real arrays
      END TYPE DYN_R1D
!-----------------------------------------------------------------------
      TYPE DYN_R2D
        REAL(KIND=kind_io4),DIMENSION(:,:),   POINTER :: NAME   !<-- Pointer for 2D real arrays
      END TYPE DYN_R2D
!-----------------------------------------------------------------------
      TYPE DYN_R3D
        REAL(KIND=kind_io4),DIMENSION(:,:,:), POINTER :: NAME   !<-- Pointer for 3D real arrays
      END TYPE DYN_R3D
!-----------------------------------------------------------------------
!
      CONTAINS
!-----------------------------------------------------------------------
!
      SUBROUTINE POINT_DYNAMICS_OUTPUT_GFS(INT_STATE,IMP_STATE_WRITE)
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
       use gfs_dyn_mpi_def
       use gfs_dyn_resol_def,      only : thermodyn_id,sfcpress_id,ngrids_gg
       use gfs_dyn_coordinate_def, only : vertcoord_id,AK5,BK5,CK5
       use gfs_dyn_vert_def,       only : si
       use gfs_dyn_date_def,       only : fhour,idate
       use gfs_dyn_write_state,    only : buff_mult_pieceg
       use gfs_dyn_io_header,      only : icen2,ienst,iensi,itrun
       use gfs_dyn_tracer_const,   only : cpi,ri
       use namelist_dynamics_def,  only : hybrid,gen_coord_hybrid
!
!-----------------------------------------------------------------------
!***  THIS ROUTINE TAKES THE USER'S SELECTIONS FOR OUTPUT QUANTITIES,
!***  POINTS AT THEM, AND INSERTS THOSE POINTERS INTO THE IMPORT STATE
!***  OF THE WRITE COMPONENTS.
!-----------------------------------------------------------------------
!
      TYPE(ESMF_State)    ,INTENT(INOUT)        :: IMP_STATE_WRITE      !<-- Import state for the Write components
      TYPE(gfs_dynamics_internal_state),POINTER :: INT_STATE            !<-- The dynamics internal state
!
!-----------------------------------------------------------------------
!
      TYPE(ESMF_Grid)             :: GRID           !<-- Create The ESMF Grid for write bundle
      INTEGER                     :: MYPE,RC,RC_DYN_OUT
      INTEGER                     :: minidx(2),maxidx(2)
      TYPE(ESMF_DistGrid)         :: DISTGRID
      TYPE(ESMF_FieldBundle),SAVE :: GFS_DYN_BUNDLE
      character(esmf_maxstr)      :: MESSAGE_CHECK

      INTEGER , DIMENSION(:, :), POINTER :: i2
!
!-----------------------------------------------------------------------
!***  ESMF VERSIONS OF THE LOGICALS IN THE DYNAMICS INTERNAL STATE
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!***  ARRAYS OF POINTERS OF THE ABOVE TYPES
!-----------------------------------------------------------------------
!
      TYPE(DYN_ISC),DIMENSION(MAX_KOUNT) :: I_SC
      TYPE(DYN_RSC),DIMENSION(MAX_KOUNT) :: R_SC
      TYPE(DYN_LSC),DIMENSION(MAX_KOUNT) :: L_SC
      TYPE(DYN_I1D),DIMENSION(MAX_KOUNT) :: I_1D
      TYPE(DYN_R1D),DIMENSION(MAX_KOUNT) :: R_1D
      TYPE(DYN_R2D),DIMENSION(MAX_KOUNT) :: R_2D
      TYPE(DYN_R3D),DIMENSION(MAX_KOUNT) :: R_3D
!
!-----------------------------------------------------------------------
!
      MYPE=int_state%ME
!
!-----------------------------------------------------------------------
!*** create grid for write bundle on belted domain
!-----------------------------------------------------------------------
!
!*** fortran array buff_mult_pieceg is on dimension 
!*** (int_state%lonf,int_state%lats_node_a_max)
     minidx(1) = 1
     minidx(2) = 1
     maxidx(1) = size(buff_mult_pieceg,1)*int_state%nodes
     maxidx(2) = size(buff_mult_pieceg,2)

!    if (mype == 0) write(0,*)'lonf=',int_state%lonf,'lats_node_a_max=',      &
!                   int_state%lats_node_a_max,'minidx=',minidx,'maxidx=',maxidx
!
!create dist  grid
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Create DISGRID for write Bundle"
      CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      DistGrid = ESMF_DistGridCreate(minidx, maxidx, rc=rc)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL gfs_dynamics_err_msg(RC,MESSAGE_CHECK,RC_DYN_OUT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!jw double check the grid index with buff_mult_pieceg

      allocate(i2(2,int_state%nodes))
      i2 = 0
      CALL ESMF_DistGridGet(DistGrid, indexCountPDe=i2, rc=rc)
!     print *,'dist grid dimension,i2=',i2(1,:), '2d i2=',i2(2,:)
!
!create grid
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Create GRID from DistGrid for write Bundle"
      CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      grid = ESMF_GridCreate(name="gridwrt", distgrid=DistGrid, rc=rc)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL gfs_dynamics_err_msg(RC,MESSAGE_CHECK,RC_DYN_OUT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  THE DYNAMICS INTERNAL STATE POINTER BLOCK
!-----------------------------------------------------------------------
!***  POINT AT ALL VARIABLES IN THE DYNAMICS INTERNAL STATE 
!***  THAT COULD BE WRITTEN TO HISTORY OUTPUT, I.E., THOSE
!***  LISTED AT THE TOP OF THIS MODULE.
!-----------------------------------------------------------------------
!
!***  INTEGER SCALARS
!
      I_SC(1)%NAME  => int_state%latg
      I_SC(2)%NAME  => int_state%lonf
      I_SC(3)%NAME  => int_state%levs
      I_SC(4)%NAME  => int_state%jcap
      I_SC(5)%NAME  => int_state%ntoz
      I_SC(6)%NAME  => int_state%ntcw
      I_SC(7)%NAME  => int_state%ncld
      I_SC(8)%NAME  => int_state%ntke
      I_SC(9)%NAME  => int_state%ntrac
      I_SC(10)%NAME => vertcoord_id
      I_SC(11)%NAME => thermodyn_id
      I_SC(12)%NAME => sfcpress_id
      I_SC(13)%NAME => ienst
      I_SC(14)%NAME => iensi
      I_SC(15)%NAME => itrun
      I_SC(16)%NAME => icen2
!        
!***  REAL SCALARS
!     ------------
!
      R_SC(1)%NAME => int_state%pdryini
      R_SC(2)%NAME => fhour
!        
!***  LOGICAL SCALARS
!     ---------------
!
      L_SC(1)%NAME => hybrid
      L_SC(2)%NAME => gen_coord_hybrid
      L_SC(3)%NAME => int_state%adiabatic
!        
!***  1D INTEGER ARRAYS
!     -----------------
!
      I_1D(1)%NAME => IDATE
!        
!***  1D REAL ARRAYS
!       ------------
!
      R_1D(1)%NAME => ak5
      R_1D(2)%NAME => bk5
      R_1D(3)%NAME => ck5
      R_1D(4)%NAME => si
      R_1D(5)%NAME => cpi
      R_1D(6)%NAME => ri
!        
!***  2D REAL ARRAYS
!
!***  3D REAL ARRAYS
!     --------------
!
      R_3D(1)%NAME => buff_mult_pieceg
!
!-----------------------------------------------------------------------
!##jw
!***  use addf90arraytostate to add output fields into import_write_state
!***  since the fields may be in different grid, then add 
!***  import_write_state into export_dyn_state
!##jw
!***  CREATE AN ESMF Bundle THAT WILL HOLD HISTORY OUTPUT DATA
!***  AND NOTHING ELSE.  THIS WILL SERVE TO ISOLATE THE OUTPUT
!***  DATA FROM EVERYTHING ELSE INSIDE THE WRITE COMPONENT'S IMPORT STATE.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK = "Create History Data Bundle"
      CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      GFS_DYN_BUNDLE = ESMF_FieldBundleCreate(name=int_state%filename_base(1) &  !<-- The Bundle's name
                                           ,rc  =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL gfs_dynamics_err_msg(RC,MESSAGE_CHECK,RC_DYN_OUT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  FIRST ADD THE LOCAL SUBDOMAIN LIMITS TO THE WRITE COMPONENT'S
!***  IMPORT STATE AS ATTRIBUTES ALONG WITH THE GLOBAL/REGIONAL MODE.
!***  THIS INFORMATION IS NEEDED FOR QUILTING THE LOCAL DOMAIN DATA
!***  INTO FULL DOMAIN FIELDS.
!***  THE LOCAL DOMAIN LIMITS GO DIRECTLY INTO THE WRITE COMPONENT'S
!***  IMPORT STATE TO KEEP THEM SEPARATE FROM THE HISTORY DATA THAT 
!***  WILL BE INSERTED INTO A Bundle.
!
!***  DO THE SAME WITH THE NUMBER OF FCST TASKS (INPESxJNPES) SINCE
!***  THE WRITE COMPONENT ALSO NEEDS THAT INFORMATION AS WELL AS
!***  THE HALO DEPTHS.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK = "Add Local Subdomain Limits to the Write Import State"
      CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_AttributeSet(state    =IMP_STATE_WRITE                  &  !<-- The Write component import state
                            ,name     ='im'                             &  !<-- Name of the integer array
                            ,value    =int_state%lonf                   &  !<-- The array being inserted into the import state
                            ,rc       =RC)

      CALL ESMF_AttributeSet(state    =IMP_STATE_WRITE                  &  !<-- The Write component import state
                            ,name     ='jm'                             &!<-- Name of the integer array
                            ,value    =int_state%latg                   &  !<-- The array being inserted into the import state
                            ,rc       =RC)
!
      CALL ESMF_AttributeSet(state    =IMP_STATE_WRITE                  &  !<-- The Write component import state
                            ,name     ='lats_node_a'                    &  !<-- Name of the integer array
                            ,value    =int_state%lats_node_a            &  !<-- The array being inserted into the import state
                            ,rc       =RC)
!
      CALL ESMF_AttributeSet(state    =IMP_STATE_WRITE                  &  !<-- The Write component import state
                            ,name     ='ipt_lats_node_a'                &  !<-- Name of the integer array
                            ,value    =int_state%ipt_lats_node_a        &  !<-- The array being inserted into the import state
                            ,rc       =RC)

      CALL ESMF_AttributeSet(state    =IMP_STATE_WRITE                  &  !<-- The Write component import state
                            ,name     ='global_lats_a'                  &  !<-- Name of the integer array
                            ,itemCount=int_state%latg                   &  !<-- Length of array being inserted into the import state
                            ,valueList=int_state%global_lats_a          &  !<-- The array being inserted into the import state
                            ,rc       =RC)

      CALL ESMF_AttributeSet(state    =IMP_STATE_WRITE                  &  !<-- The Write component import state
                            ,name     ='zhour'                          &  !<-- Name of the integer array
                            ,value    =int_state%zhour                  &  !<-- The array being inserted into the import state
                            ,rc       =RC)
!
      CALL ESMF_AttributeSet(state    =IMP_STATE_WRITE                  &  !<-- The Write component import state
                            ,name     ='pdryini'                        &  !<-- Name of the integer array
                            ,value    =int_state%pdryini                &  !<-- The array being inserted into the import state
                            ,rc       =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL gfs_dynamics_err_msg(RC,MESSAGE_CHECK,RC_DYN_OUT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  NOW INSERT INTO THE WRITE COMPONENTS' IMPORT STATE THE POINTERS
!***  OF ONLY THOSE QUANTITIES THAT ARE SPECIFIED BY THE USER FOR
!***  HISTORY OUTPUT.  THE DATA IS PLACED INTO AN ESMF Bundle WHICH
!***  ITSELF WILL BE PLACED INTO THE IMPORT STATE AT THE END OF THE ROUTINE.
!-----------------------------------------------------------------------
!
      CALL ADD_BUNDLE_TO_WRITE(int_state, imp_state_write, semilag,     &
                               grid, I_SC, R_SC, L_SC, I_1D,            &
                               R_1D, R_2D, R_3D,            &
                               2,'OGFS_DYN','OGFS_SIG',GFS_DYN_BUNDLE)
!
!-----------------------------------------------------------------------
!
      END SUBROUTINE POINT_DYNAMICS_OUTPUT_GFS
!
!-----------------------------------------------------------------------
!
      SUBROUTINE ADD_BUNDLE_TO_WRITE(int_state, imp_state_write, semilag,   &
                                     GRID, I_SC, R_SC, L_SC, I_1D,          &
                                     R_1D, R_2D, R_3D,                      &
                                     file_index,file_gen,file_id,file_bundle)
!
!-----------------------------------------------------------------------
!
      type(gfs_dynamics_internal_state),intent(in) :: int_state
      type(ESMF_STATE),intent(inout)       :: imp_state_write
      type(ESMF_GRID),intent(in)           :: GRID
      character(*),intent(in)              :: file_id,file_gen
      integer,intent(in)                   :: file_index
      type(ESMF_FieldBUNDLE),intent(inout) :: file_bundle
      TYPE(DYN_ISC),intent(in) :: I_SC(:)
      TYPE(DYN_RSC),intent(in) :: R_SC(:)
      TYPE(DYN_LSC),intent(in) :: L_SC(:)
      TYPE(DYN_I1D),intent(in) :: I_1D(:)
      TYPE(DYN_R1D),intent(in) :: R_1D(:)
      TYPE(DYN_R2D),intent(in) :: R_2D(:)
      TYPE(DYN_R3D),intent(in) :: R_3D(:)
      logical,      intent(in) :: semilag
!
!--- local variables
!
      INTEGER                     :: K,LENGTH,i,ii,ij,ik,mype,nskp      &
                                    ,N,NDIM3,NFIND,NUM_2D_FIELDS        &
                                    ,LONF,LEV,NTRACS,RC,RC_DYN_OUT      &
                                    ,LDIM1, LDIM2, UDIM1, UDIM2

      REAL(KIND=kind_io4),DIMENSION(:,:),POINTER            :: TEMP_R2D
      REAL(KIND=kind_io4),DIMENSION(:,:),allocatable,target :: mytmp
!
      CHARACTER(3)                :: MODEL_LEVEL
      CHARACTER(2)                :: TRACERS_KIND
      CHARACTER(6)                :: FMT='(I3.3)'
      CHARACTER(ESMF_MAXSTR)      :: VBL_NAME
!
      CHARACTER(10),DIMENSION(:,:),POINTER :: DYN_INT_STATE_3D_R
!
      TYPE(ESMF_Field)         :: FIELD
      TYPE(ESMF_DataCopy_Flag) :: COPYFLAG=ESMF_DATACOPY_REFERENCE

      LOGICAL                     :: LOG_ESMF

      character(esmf_maxstr)      :: MESSAGE_CHECK
!
!-----------------------------------------------------------------------
!***  BEGIN WITH THE INTEGER SCALARS.
!-----------------------------------------------------------------------
!
      MYPE = int_state%ME
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK = "Insert Dynamics Integer Scalars into Hist/Res Bundles"
      CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      DO NFIND=1,MAX_KOUNT
        IF(trim(DYN_INT_STATE_ISCALAR(2,NFIND)) == 'OGFS_SIG') THEN        !<-- Take integer scalar data specified for history output
          VBL_NAME = TRIM(DYN_INT_STATE_ISCALAR(1,NFIND))

          CALL ESMF_AttributeSet(FIELDBUNDLE=file_bundle                &  !<-- The Write component output history Bundle
                                ,name       =VBL_NAME                   &  !<-- Name of the integer scalar
                                ,value      =I_SC(NFIND)%NAME           &  !<-- The scalar being inserted into the import state
                                ,rc         =RC)

!          if(MYPE==0) write(0,*)'dyn_output,after VBL_NAME=',VBL_NAME,'NFIND=',NFIND,  &
!           'value=',I_SC(NFIND)%NAME
!
        ENDIF
!
        IF(DYN_INT_STATE_ISCALAR(2,NFIND) == '*' .AND.                  &
           DYN_INT_STATE_ISCALAR(3,NFIND) == '*'      )THEN                !<-- End of the integer scalar list
          EXIT
        ENDIF
!
      ENDDO
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL gfs_dynamics_err_msg(RC,MESSAGE_CHECK,RC_DYN_OUT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  THE REAL SCALARS.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK = "Insert Dynamics Real Scalars into History/Res Bundles"
      CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      DO NFIND=1,MAX_KOUNT
        IF(trim(DYN_INT_STATE_RSCALAR(2,NFIND)) == 'OGFS_SIG') THEN        !<-- Take real scalar data specified for history output
          VBL_NAME = TRIM(DYN_INT_STATE_RSCALAR(1,NFIND))
!
!          write(0,*)'dyn_output,VBL_NAME=',VBL_NAME,'NFIND=',NFIND,     &
!           'value=',R_SC(NFIND)%NAME

          CALL ESMF_AttributeSet(FIELDBUNDLE=file_bundle                &  !<-- The Write component output history Bundle
                                ,name       =VBL_NAME                   &  !<-- Name of the integer scalar
                                ,value      =R_SC(NFIND)%NAME           &  !<-- The scalar being inserted into the import state
                                ,rc         =RC)
        ENDIF
!
        IF(DYN_INT_STATE_RSCALAR(2,NFIND) == '*' .AND.                  &
           DYN_INT_STATE_RSCALAR(3,NFIND) == '*'      )THEN                !<-- End of the real scalar list
          EXIT
        ENDIF
!
      ENDDO
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL gfs_dynamics_err_msg(RC,MESSAGE_CHECK,RC_DYN_OUT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  THE LOGICAL SCALARS.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK = "Insert Dynamics Logical Scalars into History/Res Bundles"
      CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      DO NFIND=1,MAX_KOUNT
!
        IF(trim(DYN_INT_STATE_LSCALAR(2,NFIND)) == 'OGFS_SIG') THEN        !<-- Take real scalar data specified for history output
          VBL_NAME = TRIM(DYN_INT_STATE_LSCALAR(1,NFIND))
!
!***  ESMF VERSION OF LOGICALS NEEDED FOR THEIR INSERTION
!***  INTO THE HISTORY OUTPUT Bundle OF THE WRITE COMPONENT'S
!***  IMPORT STATE.
!
!          write(0,*)'dyn_output,log VBL_NAME=',VBL_NAME,'NFIND=',NFIND

          LOG_ESMF = .false.
          if(L_SC(NFIND)%NAME) LOG_ESMF = .true.
          CALL ESMF_AttributeSet(FIELDBUNDLE=file_bundle                &  !<-- The Write component output history Bundle
                                ,name       =VBL_NAME                   &  !<-- Name of the integer scalar
                                ,value      =LOG_ESMF                   &  !<-- The scalar being inserted into the import state
                                ,rc         =RC)
        ENDIF
!
        IF(DYN_INT_STATE_LSCALAR(2,NFIND) == '*' .AND.                  &
           DYN_INT_STATE_LSCALAR(3,NFIND) == '*'      ) THEN               !<-- End of the real scalar list
          EXIT
        ENDIF
!
      ENDDO
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL gfs_dynamics_err_msg(RC,MESSAGE_CHECK,RC_DYN_OUT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  THE 1D INTEGER ARRAYS
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK = "Insert Dynamics 1-D Integer Arrays into Hist/Res Bundles"
      CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      DO NFIND=1,MAX_KOUNT
        IF(trim(DYN_INT_STATE_1D_I(2,NFIND)) == 'OGFS_SIG')THEN            !<-- Take 1D integer array data specified for history output
          VBL_NAME = TRIM(DYN_INT_STATE_1D_I(1,NFIND))
          LENGTH   = SIZE(I_1D(NFIND)%NAME)
!
!          write(0,*)'dyn_output,1d int VBL_NAME=',VBL_NAME,'NFIND=',NFIND,     &
!           'LENGTH=',LENGTH,'value=',I_1D(NFIND)%NAME(1:LENGTH)

          CALL ESMF_AttributeSet(FIELDBUNDLE =file_bundle               &  !<-- The Write component output history Bundle
                                ,name        =VBL_NAME                  &  !<-- Name of the integer scalar
                                ,itemCount   =LENGTH                    &  !<-- # of elements in this attribute
                                ,valueList   =I_1D(NFIND)%NAME          &  !<-- The 1D integer being inserted into the import state
                                ,rc          =RC)
        ENDIF
!
        IF(DYN_INT_STATE_1D_I(2,NFIND) == '*' .AND.                     &
           DYN_INT_STATE_1D_I(3,NFIND) == '*'      )THEN                   !<-- End of the 1D integer array list
          EXIT
        ENDIF
!
      ENDDO
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL gfs_dynamics_err_msg(RC,MESSAGE_CHECK,RC_DYN_OUT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  THE 1D REAL ARRAYS.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK = "Insert Dynamics 1-D Real Arrays into Hist/Res Bundles"
      CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      DO NFIND=1,MAX_KOUNT
        IF(trim(DYN_INT_STATE_1D_R(2,NFIND)) == 'OGFS_SIG') THEN           !<-- Take 1D real array data specified for history output
          VBL_NAME = TRIM(DYN_INT_STATE_1D_R(1,NFIND))
          LENGTH   = SIZE(R_1D(NFIND)%NAME)
          if(VBL_NAME == 'CPI' .or. VBL_NAME == 'RI') LENGTH = int_state%ntrac+1
!
!          write(0,*)'dyn_output,1d real VBL_NAME=',VBL_NAME,'NFIND=',NFIND,     &
!           'LENGTH=',LENGTH,'value=',R_1D(NFIND)%NAME(1:LENGTH)

          CALL ESMF_AttributeSet(FIELDBUNDLE =file_bundle               &  !<-- The Write component output history Bundle
                                ,name        =VBL_NAME                  &  !<-- Name of the integer scalar
                                ,itemCount   =LENGTH                    &  !<-- # of elements in this attribute
                                ,valueList   =R_1D(NFIND)%NAME          &  !<-- The 1D real being inserted into the import state
                                ,rc          =RC)
        ENDIF
!
        IF(DYN_INT_STATE_1D_R(2,NFIND) == '*' .AND.                     &
           DYN_INT_STATE_1D_R(3,NFIND) == '*'      ) THEN                  !<-- End of the 1D real array list
          EXIT
        ENDIF
!
      ENDDO
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL gfs_dynamics_err_msg(RC,MESSAGE_CHECK,RC_DYN_OUT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  DEREFERENCE THE LOCAL MEMORY LIMTS.  THESE DETERMINE THE TOTAL
!***  NUMBER OF WORDS LOADED PER SUBDOMAIN.
!-----------------------------------------------------------------------
!
      lonf = int_state%lonf
      lev  = 1
!
!-----------------------------------------------------------------------
!  THE 2D REAL ARRAYS on GRID3.
!-----------------------------------------------------------------------
!
      NUM_2D_FIELDS = 0
!
      DO NFIND=1,MAX_KOUNT
        IF(trim(DYN_INT_STATE_2D_R(2,NFIND)) == 'OGFS_SIG') THEN           !<-- Take 2D real array data specified for history output
          VBL_NAME = TRIM(DYN_INT_STATE_2D_R(1,NFIND))
          TEMP_R2D => R_3D(1)%NAME(1:int_state%lonf,1:int_state%lats_node_a_max,NFIND)
!
!         write(0,*)'dyn_output,2d real VBL_NAME=',VBL_NAME,'NFIND=',NFIND, &
!            maxval(temp_r2d),maxval(temp_r2d),' me=',mype
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          MESSAGE_CHECK = "Insert Dynamics 2-D Real Data into Field"
          CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          FIELD = ESMF_FieldCreate(            GRID                     &  !<-- The ESMF grid
                                ,              TEMP_R2D                 &  !<-- The 2D real array being inserted into the import state
                                ,datacopyflag =COPYFLAG                 &
                                ,name         =VBL_NAME                 &  !<-- Name of the 2D real array
                                ,indexFlag=ESMF_INDEX_DELOCAL           &
                                ,rc           =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          CALL gfs_dynamics_err_msg(RC,MESSAGE_CHECK,RC_DYN_OUT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          MESSAGE_CHECK = "Insert Dynamics 2-D Real Field into Hist/Res Bundles"
          CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          CALL ESMF_FieldBundleAdd(file_bundle &  !<-- The Write component output history Bundle
                                  ,(/FIELD/)   &  !<-- ESMF Field holding the 2D real array
                                  ,rc    =RC)
          CALL gfs_dynamics_err_msg(RC,MESSAGE_CHECK,RC_DYN_OUT)
!
          NUM_2D_FIELDS = NUM_2D_FIELDS + 1
!
        ENDIF
!
        IF(DYN_INT_STATE_2D_R(2,NFIND) == '*' .AND.                     &
           DYN_INT_STATE_2D_R(3,NFIND) == '*'      ) THEN                  !<-- End of the 2D real array list
          EXIT
        ENDIF
!
      ENDDO
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL gfs_dynamics_err_msg(RC,MESSAGE_CHECK,RC_DYN_OUT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  THE 3D REAL ARRAYS  on GRID3 too
!***  WE ARE WORKING WITH 3D ARRAYS BUT THEY ARE LOADED LAYER BY LAYER
!***  INTO 2D Fields.
!-----------------------------------------------------------------------
!

!! generalized tracer:
!! set up DYN_INT_STATE_3D_R_DIAB(1:3,1:MAX_KOUNT)

      DO II =1,MAX_KOUNT
         DYN_INT_STATE_3D_R_DIAB(1:3,II) = '**********'
      ENDDO
!
!! find first tracer, !!! first tracer will always be "spfh'
      DO II =1,MAX_KOUNT
         if(trim(DYN_INT_STATE_3D_R_ADIAB(1,II)) == 'spfh') then
           NTRACS = II
           exit
         endif
      ENDDO
!  fill in met fields
      DO II=1,NTRACS-1
         DYN_INT_STATE_3D_R_DIAB(1,II) = DYN_INT_STATE_3D_R_ADIAB(1,II)
      ENDDO

!  fill in tracer fields
      II = NTRACS-1
      DO IJ=1, gfs_dyn_tracer%ntrac
        DYN_INT_STATE_3D_R_DIAB(1,II+IJ) = gfs_dyn_tracer%vname(IJ, 1)
      ENDDO                                        
      IK = II  + gfs_dyn_tracer%ntrac     

!  ntke
      if( int_state%ntke > 0 ) then
        DYN_INT_STATE_3D_R_DIAB(1,IK+1) = 'tke       '
        IK = IK + 1
      endif
!  fill in file type and level type
      DO II = 1, IK
        DYN_INT_STATE_3D_R_DIAB(2,II) = 'OGFS_SIG  '
        DYN_INT_STATE_3D_R_DIAB(3,II) = 'levs      '
      ENDDO
!    

!      IF(.NOT. int_state%ADIABATIC)THEN
      DYN_INT_STATE_3D_R => DYN_INT_STATE_3D_R_DIAB(1:3,1:MAX_KOUNT)      !<-- Select diabatic output
!     ELSE
!       DYN_INT_STATE_3D_R => DYN_INT_STATE_3D_R_ADIAB(1:3,1:MAX_KOUNT)     !<-- Select adiabatic output
!     ENDIF
!
      NULLIFY(TEMP_R2D)
      LDIM1 = LBOUND(R_3D(1)%NAME,1)
      UDIM1 = UBOUND(R_3D(1)%NAME,1)
      LDIM2 = LBOUND(R_3D(1)%NAME,2)
      UDIM2 = UBOUND(R_3D(1)%NAME,2)

!     write(0,*)'LDIM1=',LDIM1,'UDIM1=',UDIM1,'LDIM2=',LDIM2,'UDIM2=',UDIM2
!
      nskp = 0
      DO NFIND=1,MAX_KOUNT
!
        IF(trim(DYN_INT_STATE_3D_R(2,NFIND)) == 'OGFS_SIG') THEN            !<-- Take 3D real array data specified for history output

          if(trim(DYN_INT_STATE_3D_R(3,NFIND)) == 'levs')then
            NdiM3 = int_state%levs
          elseif(trim(DYN_INT_STATE_3D_R(3,NFIND)) == 'levsp1')then
            NDIM3 = int_state%levs+1
          endif
          if (semilag .and. (trim(DYN_INT_STATE_3D_R(1,NFIND)) == 'dpres' &
                        .or. trim(DYN_INT_STATE_3D_R(1,NFIND)) == 'pres')) then
            nskp = nskp + ndim3
            cycle
          endif
!
          DO K=1,NDIM3
            WRITE(MODEL_LEVEL,FMT)K
            VBL_NAME = TRIM(DYN_INT_STATE_3D_R(1,NFIND))//'_'//MODEL_LEVEL//'_2D'
            TEMP_R2D => R_3D(1)%NAME(LDIM1:UDIM1,LDIM2:UDIM2,NUM_2D_FIELDS+K+nskp)
!
!           write(0,*)'dyn_output,3d real VBL_NAME=',VBL_NAME,'NFIND=',NFIND,     &
!           'value=',maxval(R_3D(1)%NAME(LDIM1:UDIM1,LDIM2:UDIM2,NUM_2D_FIELDS+K+nskp))&
!           ,minval(R_3D(1)%NAME(LDIM1:UDIM1,LDIM2:UDIM2,NUM_2D_FIELDS+K+nskp)),' me=',mype
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
            MESSAGE_CHECK = "Fill 2-D Fields with Each Level of Dynamics 3-D Data"
            CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
            FIELD = ESMF_FieldCreate(              GRID                 &  !<-- The ESMF grid
                                   ,              TEMP_R2D              &  !<-- Level K of 3D real array being inserted into the import state
                                   ,datacopyflag =COPYFLAG              &
                                   ,name         =VBL_NAME              &  !<-- Name of this level of the 3D real array
                                   ,rc           =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
            CALL gfs_dynamics_err_msg(RC,MESSAGE_CHECK,RC_DYN_OUT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
            MESSAGE_CHECK = "Insert Dynamics 3-D Data into History Bundle"
            CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
            CALL ESMF_FieldBundleAdd(file_bundle &  !<-- The Write component output history Bundle
                                    ,(/FIELD/)   &  !<-- ESMF Field holding the 1D real array
                                    ,rc    =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
            CALL gfs_dynamics_err_msg(RC,MESSAGE_CHECK,RC_DYN_OUT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          ENDDO
!
          NUM_2D_FIELDS = NUM_2D_FIELDS + NDIM3
!
        ENDIF
!
        IF(DYN_INT_STATE_3D_R(2,NFIND)=='*' .AND.                       &
           DYN_INT_STATE_3D_R(3,NFIND)=='*'      )THEN                     !<-- End of the 3D real array list
          EXIT
        ENDIF
!
      ENDDO
!
!-----------------------------------------------------------------------
!***  THE 4D REAL ARRAYS.
!***  WE ARE WORKING WITH 4D ARRAYS BUT THEY ARE LOADED LAYER BY LAYER
!***  INTO 2D Fields.
!-----------------------------------------------------------------------
!
!
!jw      DO NFIND=1,MAX_KOUNT
!
!jw        IF(trim(DYN_INT_STATE_4D_R(2,NFIND))=='OGFS_SIG')THEN           !<-- Take 4D real array data specified for history output
!
!jw         if(trim(DYN_INT_STATE_4D_R(3,NFIND))=='levs')then
!jw           NDIM3=levs
!jw         elseif(trim(DYN_INT_STATE_4D_R(3,NFIND))=='levsp1')then
!jw           NDIM3=levs+1
!jw         endif
!jw         LDIM1=LBOUND(R_4D(NFIND)%NAME,1)
!jw         UDIM1=UBOUND(R_4D(NFIND)%NAME,1)

!jw         DO N=int_state%ntcw+1,ntrac                                    !<-- Loop through the tracers kind
!jw          DO K=1,NDIM3
!
!jw            WRITE(TRACERS_KIND,FMT)N
!jw            WRITE(MODEL_LEVEL,FMT)K
!
!jw            VBL_NAME=TRIM(DYN_INT_STATE_4D_R(1,NFIND))//'_'//TRACERS_KIND//'_'//MODEL_LEVEL//'_2D'
!jw            NULLIFY(TEMP_R2D)
!jw            TEMP_R2D=>R_2D(NFIND)%NAME(LDIM1:UDIM1,NUM_2D_FIELDS)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!jw            MESSAGE_CHECK="Fill 2-D Fields with Each Level of Physics 4-D Data"
!      CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!jw            FIELD=ESMF_FieldCreate(grid         =GRID3               &  !<-- The ESMF grid
!jw                                  ,farray       =TEMP_R2D            &  !<-- Level K of 4D real array being inserted into the data Bundle
!jw                                  ,copyflag     =COPYFLAG            &
!jw                                  ,name         =VBL_NAME            &  !<-- Name of this level of
!jwthe 4D real array
!jw                                  ,rc           =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!jw            CALL gfs_dynamics_err_msg(RC,MESSAGE_CHECK,RC_DYN_OUT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!jw            MESSAGE_CHECK="Insert Physics 4-D Data into History Bundle"
!      CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!jw            CALL ESMF_FieldBundleAdd(bundle=GFS_DYN_BUNDLE           &  !<-- The write component's output data Bundle
!jw                                    ,field =FIELD                    &  !<-- ESMF Field holding the 1D real array
!jw                                    ,rc    =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!jw            CALL gfs_dynamics_err_msg(RC,MESSAGE_CHECK,RC_DYN_OUT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!jw            NUM_2D_FIELDS=NUM_2D_FIELDS+1                              !<-- Continue adding up all levels of 4D Fields
!
!jw          ENDDO
!jw         ENDDO
!
!jw        ENDIF
!
!jw        IF(DYN_INT_STATE_4D_R(2,NFIND)=='*' .AND.                    &
!jw           DYN_INT_STATE_4D_R(3,NFIND)=='*'      )THEN                  !<-- End of the 4D real array list
!jw          EXIT
!jw        ENDIF
!
!jw      ENDDO
!-----------------------------------------------------------------------
!
      IF(MYPE == 0) WRITE(0,*)' Exit DYNAMICS_OUTPUT num_2d_fields=',NUM_2D_FIELDS
!
!-----------------------------------------------------------------------
!***  INSERT THE HISTORY DATA Bundle INTO THE WRITE COMPONENT'S
!***  IMPORT STATE.
!***  SINCE DYNAMICS IS CALLED BEFORE PHYSICS, WE WILL INSERT THE
!***  BUNDLE NOW AND SIMPLY USE IT IN POINT_PHYSICS_OUTPUT.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Dynamics: Insert History Bundle into the Write Import State"
      CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_StateAddReplace(IMP_STATE_WRITE &  !<-- The Write component's import state
                               ,(/file_bundle/) &  !<-- The ESMF Bundle holding all Dynamics history data
                               ,rc         =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL gfs_dynamics_err_msg(RC,MESSAGE_CHECK,RC_DYN_OUT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!
      END SUBROUTINE ADD_BUNDLE_TO_WRITE
!
!-----------------------------------------------------------------------
!
      END MODULE gfs_dynamics_output
