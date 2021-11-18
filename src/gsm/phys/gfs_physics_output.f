!#include "../../../ESMFVersionDefine.h"

!-----------------------------------------------------------------------
!
      MODULE gfs_physics_output
!
!-----------------------------------------------------------------------
!***  LIST THE QUANTITIES FROM THE PHYSICS INTERNAL STATE THAT
!***  CAN BE SELECTED FOR HISTORY OUTPUT AND ASSOCIATE THEM WITH
!***  UNIQUE INTEGERS.
!***  THE USER WILL PLACE AN 'H' IN THE 3rd ELEMENT OF THE FOLLOWING
!***  LIST OF QUANTITIES THAT ARE IN THE PHYSICS INTERNAL STATE
!***  IF THAT QUANTITY IS TO BE WRITTEN TO THE HISTORY FILES.
!
!***  THE ORDER OF NAMES IN THE CHARACTER SELECTION LIST MUST BE
!***  THE SAME AS THE ORDER OF THE ACTUAL DATA BEING TARGETED
!***  BY POINTERS.
!-----------------------------------------------------------------------
!***  WHEN NEW QUANTITIES ARE ADDED TO THE PHYSICS INTERNAL STATE
!***  THAT ARE CANDIDATES FOR HISTORY OUTPUT THEN THEY NEED TO BE
!***  ADDED IN TWO PLACES BELOW: 
!*    (1) THE APPROPRIATE DATA LIST PRECEDING THE 'CONTAINS' STATEMENT
!*    (2) 'THE PHYSICS INTERNAL STATE POINTER BLOCK'
!*        IN SUBROUTINE POINT_PHYSICS_OUPUT
!
!***  REVISION LOG:
!***  Oct    2008  J. WANG    create gfs physics output
!***  Jul 23 2010, Sarah Lu   write out 2d aer_diag (g2d_fld)
!***  Jul    2010  S. Moorthi Updated for new physics and extra output
!***  Dec 23 2010, Sarah Lu   add g2d_fld%met to 2d aer_diag
!***  Feb    2011  Weiyu Yang, Updated to use both the ESMF 4.0.0rp2 library,
!***                         ESMF 5 library and the the ESMF 3.1.0rp2 library.
!***  Mar    2011  J. WANG    output lpl,zsoil for sfc file
!***  MAR    2011  Weiyu Yang, Use the fieldbundle list to avoid the ESMF
!***                         log error. Fix bug.
!***  Jun    2011  J. Wang   change field names in flx file to standard grib
!***                         names
!***  Sep    2011  Weiyu Yang, Updated to use both the ESMF 5.2.0r library.
!***  Apr    2013  J. Wang   add ivegsrc and cu_physics to gfs flx file
!***  Nov 22 2013  Sarah Lu  Updated PHY_INT_STATE_2D_R_FLX to match re-organized
!***                         flx files with 5 new fields
!***  Aug 21 2015 Xu Li      change nst_fcst to be nstf_name
!***  Jan    2016  S Moorthi Add tke related changes and cleanup
!
!-----------------------------------------------------------------------
!
      USE ESMF
      USE module_gfs_machine
      USE gfs_physics_internal_state_mod,ONLY: gfs_physics_internal_state
      use gfs_physics_err_msg_mod,       ONLY: gfs_physics_err_msg
!
      IMPLICIT NONE
      include 'mpif.h'
!
!-----------------------------------------------------------------------
!
      PRIVATE
!
      PUBLIC :: POINT_PHYSICS_OUTPUT_GFS,                               &
     &          PHY_INT_STATE_ISCALAR,  PHY_INT_STATE_RSCALAR,          &
     &          PHY_INT_STATE_1D_I,     PHY_INT_STATE_1D_R,             &
     &          PHY_INT_STATE_2D_R_SFC, PHY_INT_STATE_3D_R,             &
     &          PHY_INT_STATE_2D_R_FLX, PHY_INT_STATE_2D_R_AER,         &
     &          PHY_INT_STATE_2D_R_NST, PHY_INT_STATE_2D_R_FLX_CLM
!
!-----------------------------------------------------------------------
!
!     INTEGER,PARAMETER      :: MAX_KOUNT=100
      INTEGER,PARAMETER      :: MAX_KOUNT=150      
      character(esmf_maxstr) :: MESSAGE_CHECK
!
!-----------------------------------------------------------------------
!***  LIST THE VARIOUS QUANTITIES IN THE PHYSICS INTERNAL STATE 
!***  THAT ARE CANDIDATES FOR HISTORY OUTPUT.
!***  GROUP THEM BY TYPE (VARIOUS SCALARS, VARIOUS ARRAYS) SINCE
!***  WE WILL USE A WORKING POINTER FOR EACH TYPE.
!-----------------------------------------------------------------------
!
!-------------------------
!***  INTEGER SCALARS  ***
!-------------------------
!
      CHARACTER(15),DIMENSION(3,MAX_KOUNT) :: PHY_INT_STATE_ISCALAR        &
!
       =RESHAPE((/ 'latr           ', 'OGFS_PHY       ', 'R              ' &
                  ,'lonr           ', 'OGFS_PHY       ', 'R              ' &
                  ,'levs           ', 'OGFS_PHY       ', 'R              ' &
                  ,'ntoz           ', 'OGFS_PHY       ', 'R              ' &
                  ,'ntcw           ', 'OGFS_PHY       ', 'R              ' &
                  ,'ncld           ', 'OGFS_PHY       ', 'R              ' &
                  ,'ntke           ', 'OGFS_PHY       ', 'R              ' &
                  ,'ntrac          ', 'OGFS_PHY       ', 'R              ' &
                  ,'thermodyn_id   ', 'OGFS_PHY       ', 'R              ' &
                  ,'sfcpress_id    ', 'OGFS_PHY       ', 'R              ' &
                  ,'lsoil          ', 'OGFS_PHY       ', 'R              ' &
                  ,'idrt           ', 'OGFS_PHY       ', 'R              ' &
                  ,'ivssfc         ', 'OGFS_SFC       ', 'R              ' &
                  ,'ivegsrc        ', 'OGFS_FLX       ', 'R              ' &
                  ,'cu_physics     ', 'OGFS_FLX       ', 'R              ' &
                  ,'ivsnst         ', 'OGFS_NST       ', 'R              ' &
                 /)                                                        &
                ,(/3,MAX_KOUNT/)                                           &
                ,(/'*              ', '*              ', '*              '/))
!----------------------
!***  REAL SCALARS  ***
!----------------------
!
      CHARACTER(10),DIMENSION(3,MAX_KOUNT) :: PHY_INT_STATE_RSCALAR     &
!
       =RESHAPE((/ 'zhour     ', 'OGFS_FLX  ', 'R         '             &
!jw                           '-         ', '-         ', '-         '  &
                 /)                                                     &
               ,(/3,MAX_KOUNT/)                                         &
               ,(/'**********', '**********', '**********'/))
!
!----------------------------
!***  INTEGER 1-D ARRAYS  ***
!----------------------------
!
      CHARACTER(10),DIMENSION(3,MAX_KOUNT) :: PHY_INT_STATE_1D_I        &
!
       =RESHAPE((/ 'IDAT      ', 'OGFS_PHY  ', '-         '             &
                   ,'LPL       ', 'OGFS_SFC  ', '-         '            &
                 /)                                                     &
               ,(/3,MAX_KOUNT/)                                         &
               ,(/'**********', '**********', '**********'/))
!
!----------------------------
!***  INTEGER 2-D ARRAYS  ***
!----------------------------
!
      CHARACTER(10),DIMENSION(3,MAX_KOUNT) :: PHY_INT_STATE_2D_I        &
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
      CHARACTER(13),DIMENSION(3,MAX_KOUNT) :: PHY_INT_STATE_1D_R        &
!
       =RESHAPE((/ 'zsoil        ', 'OGFS_SFC     ', '             '    &
!jw               ,'BK5          ', 'OGFS_PHY     ', 'R            '    &
!jw               ,'CK5          ', 'OGFS_PHY     ', 'R            '    &
!
                 /)                                                     &
               ,(/3,MAX_KOUNT/)                                         &
               ,(/'*            ', '*            ', '*            '/))
!
!-------------------------
!***  REAL 2-D ARRAYS  sfc***
!-------------------------
!
      CHARACTER(16),DIMENSION(3,MAX_KOUNT),TARGET :: PHY_INT_STATE_2D_R_SFC    & 
!
       =RESHAPE((/ 'tmp             ', 'OGFS_SFC        ', 'sfc             ' &
                  ,'weasd           ', 'OGFS_SFC        ', 'sfc             ' &
!                 ,'tcdc            ', 'OGFS_SFC        ', 'cnvt cld layer  ' &
!                 ,'pres            ', 'OGFS_SFC        ', 'cnvt cld bot    ' &
!                 ,'pres            ', 'OGFS_SFC        ', 'cnvt cld top    ' &
                  ,'tg3             ', 'OGFS_SFC        ', 'sfc             ' &
                  ,'sfcr            ', 'OGFS_SFC        ', 'sfc             ' &
                  ,'alvsf           ', 'OGFS_SFC        ', 'sfc             ' &
                  ,'alvwf           ', 'OGFS_SFC        ', 'sfc             ' &
                  ,'alnsf           ', 'OGFS_SFC        ', 'sfc             ' &
                  ,'alnwf           ', 'OGFS_SFC        ', 'sfc             ' &
                  ,'land            ', 'OGFS_SFC        ', 'sfc             ' &
                  ,'veg             ', 'OGFS_SFC        ', 'sfc             ' &
                  ,'cnwat           ', 'OGFS_SFC        ', 'sfc             ' &
                  ,'f10m            ', 'OGFS_SFC        ', '10 m above gnd  ' &
                  ,'tmp             ', 'OGFS_SFC        ', '2 m above gnd   ' &
                  ,'spfh            ', 'OGFS_SFC        ', '2 m above gnd   ' &
                  ,'vtype           ', 'OGFS_SFC        ', 'sfc             ' &
                  ,'sotyp           ', 'OGFS_SFC        ', 'sfc             ' &
                  ,'facsf           ', 'OGFS_SFC        ', 'sfc             ' &
                  ,'facwf           ', 'OGFS_SFC        ', 'sfc             ' &
                  ,'fricv           ', 'OGFS_SFC        ', 'sfc             ' &
                  ,'ffmm            ', 'OGFS_SFC        ', 'sfc             ' &
                  ,'ffhh            ', 'OGFS_SFC        ', 'sfc             ' &
                  ,'icetk           ', 'OGFS_SFC        ', 'sfc             ' &
                  ,'icec            ', 'OGFS_SFC        ', 'sfc             ' &
                  ,'tisfc           ', 'OGFS_SFC        ', 'sfc             ' &
                  ,'tprcp           ', 'OGFS_SFC        ', 'sfc             ' &
                  ,'crain           ', 'OGFS_SFC        ', 'sfc             ' &
                  ,'snod            ', 'OGFS_SFC        ', 'sfc             ' &
                  ,'shdmin          ', 'OGFS_SFC        ', 'sfc             ' &
                  ,'shdmax          ', 'OGFS_SFC        ', 'sfc             ' &
                  ,'sltyp           ', 'OGFS_SFC        ', 'sfc             ' &
                  ,'salbd           ', 'OGFS_SFC        ', 'sfc             ' &
                  ,'orog            ', 'OGFS_SFC        ', 'sfc             ' &
                 /)                                                            &
               ,(/3,MAX_KOUNT/)                                               &
               ,(/'*               ', '*               ', '*               '/))
!
!-------------------------
!***  REAL 3-D ARRAYS  ***
!-------------------------
!
!<-- Append "_ikj" to any 3D variables with IKJ storage order

      CHARACTER(16),DIMENSION(4,MAX_KOUNT) :: PHY_INT_STATE_3D_R                                & 
       =RESHAPE((/ 'smc             ', 'OGFS_SFC        ','soil layer      ', 'lsoil           ' &
                  ,'stc             ', 'OGFS_SFC        ','soil layer      ', 'lsoil           ' &
                  ,'slc             ', 'OGFS_SFC        ','soil layer      ', 'lsoil           ' &
                 /)                                                                              &
               ,(/4,MAX_KOUNT/)                                                                  &
               ,(/'*               ', '*               ', '*               ','*               '/))
!
!-------------------------
!***  REAL 4-D ARRAYS  ***
!-------------------------
!
      CHARACTER(12),DIMENSION(4,MAX_KOUNT) :: PHY_INT_STATE_4D_R                 &
!
       =RESHAPE((/ '-           ', '-           ', '-           ','-           ' &
                 /)                                                              &
               ,(/4,MAX_KOUNT/)                                                  &
               ,(/'*           ', '*           ', '*           ','*           '/))
!
!-------------------------
!***  REAL 2-D ARRAYS  flx***
!

      CHARACTER(16),DIMENSION(3,MAX_KOUNT),TARGET :: PHY_INT_STATE_2D_R_FLX        &
!
       =RESHAPE((/ 'uflx_ave        ', 'OGFS_FLX        ', 'sfc             ' &   !1
                  ,'vflx_ave        ', 'OGFS_FLX        ', 'sfc             ' &   !2
                  ,'shtfl_ave       ', 'OGFS_FLX        ', 'sfc             ' &   !3
                  ,'lhtfl_ave       ', 'OGFS_FLX        ', 'sfc             ' &   !4
                  ,'tmp             ', 'OGFS_FLX        ', 'sfc             ' &   !5
                  ,'soilw           ', 'OGFS_FLX        ', '0-10 cm down    ' &   !6
                  ,'soilw           ', 'OGFS_FLX        ', '10-40 cm down   ' &   !7
                  ,'tmp             ', 'OGFS_FLX        ', '0-10 cm down    ' &   !8
                  ,'tmp             ', 'OGFS_FLX        ', '10-40 cm down   ' &   !9
                  ,'weasd           ', 'OGFS_FLX        ', 'sfc             ' &   !10
                  ,'ulwrf_ave       ', 'OGFS_FLX        ', 'nom. top        ' &   !11
                  ,'uswrf_ave       ', 'OGFS_FLX        ', 'nom. top        ' &   !12
                  ,'uswrf_ave       ', 'OGFS_FLX        ', 'sfc             ' &   !13
                  ,'dswrf_ave       ', 'OGFS_FLX        ', 'sfc             ' &   !14
                  ,'tcdc_ave        ', 'OGFS_FLX        ', 'high cld lay    ' &   !15
                  ,'pres_ave        ', 'OGFS_FLX        ', 'high cld top    ' &   !16
                  ,'pres_ave        ', 'OGFS_FLX        ', 'high cld bot    ' &   !17
                  ,'tmp_ave         ', 'OGFS_FLX        ', 'high cld top    ' &   !18
                  ,'tcdc_ave        ', 'OGFS_FLX        ', 'mid cld lay     ' &   !19
                  ,'pres_ave        ', 'OGFS_FLX        ', 'mid cld top     ' &   !20
                  ,'pres_ave        ', 'OGFS_FLX        ', 'mid cld bot     ' &   !21
                  ,'tmp_ave         ', 'OGFS_FLX        ', 'mid cld top     ' &   !22
                  ,'tcdc_ave        ', 'OGFS_FLX        ', 'low cld lay     ' &   !23
                  ,'pres_ave        ', 'OGFS_FLX        ', 'low cld top     ' &   !24
                  ,'pres_ave        ', 'OGFS_FLX        ', 'low cld bot     ' &   !25
                  ,'tmp_ave         ', 'OGFS_FLX        ', 'low cld top     ' &   !26
                  ,'tcdc_ave        ', 'OGFS_FLX        ', 'atmos col       ' &   !27
                  ,'tcdc_ave        ', 'OGFS_FLX        ', 'bndary-layer cld' &   !28
                  ,'dlwrf_ave       ', 'OGFS_FLX        ', 'sfc             ' &   !29
                  ,'ulwrf_ave       ', 'OGFS_FLX        ', 'sfc             ' &   !30
                  ,'duvb_ave        ', 'OGFS_FLX        ', 'sfc             ' &   !31
                  ,'cduvb_ave       ', 'OGFS_FLX        ', 'sfc             ' &   !32
                  ,'dswrf           ', 'OGFS_FLX        ', 'nom. top        ' &   !33
                  ,'vbdsf_ave       ', 'OGFS_FLX        ', 'sfc             ' &   !34
                  ,'vddsf_ave       ', 'OGFS_FLX        ', 'sfc             ' &   !35
                  ,'nbdsf_ave       ', 'OGFS_FLX        ', 'sfc             ' &   !36
                  ,'nddsf_ave       ', 'OGFS_FLX        ', 'sfc             ' &   !37
                  ,'csulf           ', 'OGFS_FLX        ', 'nom. top        ' &   !38
                  ,'csusf           ', 'OGFS_FLX        ', 'nom. top        ' &   !39
                  ,'csdlf           ', 'OGFS_FLX        ', 'sfc             ' &   !40
                  ,'csusf           ', 'OGFS_FLX        ', 'sfc             ' &   !41
                  ,'csdsf           ', 'OGFS_FLX        ', 'sfc             ' &   !42
                  ,'csulf           ', 'OGFS_FLX        ', 'sfc             ' &   !43
                  ,'albdo_ave       ', 'OGFS_FLX        ', 'sfc             ' &   !44
                  ,'prate_ave       ', 'OGFS_FLX        ', 'sfc             ' &   !45
                  ,'cprat_ave       ', 'OGFS_FLX        ', 'sfc             ' &   !46
                  ,'gflux_ave       ', 'OGFS_FLX        ', 'sfc             ' &   !47
                  ,'land            ', 'OGFS_FLX        ', 'sfc             ' &   !48
                  ,'icec            ', 'OGFS_FLX        ', 'sfc             ' &   !49
                  ,'ugrd            ', 'OGFS_FLX        ', '10 m above gnd  ' &   !50
                  ,'vgrd            ', 'OGFS_FLX        ', '10 m above gnd  ' &   !51
                  ,'tmp             ', 'OGFS_FLX        ', '2 m above gnd   ' &   !52
                  ,'spfh            ', 'OGFS_FLX        ', '2 m above gnd   ' &   !53
                  ,'pres            ', 'OGFS_FLX        ', 'sfc             ' &   !54
                  ,'tmax            ', 'OGFS_FLX        ', '2 m above gnd   ' &   !55
                  ,'tmin            ', 'OGFS_FLX        ', '2 m above gnd   ' &   !56
                  ,'spfhmax         ', 'OGFS_FLX        ', '2 m above gnd   ' &   !57
                  ,'spfhmin         ', 'OGFS_FLX        ', '2 m above gnd   ' &   !58
                  ,'watr_acc        ', 'OGFS_FLX        ', 'sfc             ' &   !59
                  ,'pevpr_ave       ', 'OGFS_FLX        ', 'sfc             ' &   !60
                  ,'cwork_ave       ', 'OGFS_FLX        ', 'atmos col       ' &   !61
                  ,'u-gwd_ave       ', 'OGFS_FLX        ', 'sfc             ' &   !62
                  ,'v-gwd_ave       ', 'OGFS_FLX        ', 'sfc             ' &   !63
                  ,'hpbl            ', 'OGFS_FLX        ', 'sfc             ' &   !64
                  ,'pwat            ', 'OGFS_FLX        ', 'atmos col       ' &   !65
                  ,'tcdc            ', 'OGFS_FLX        ', 'convect-cld laye' &   !66
                  ,'pres            ', 'OGFS_FLX        ', 'convect-cld top ' &   !67
                  ,'pres            ', 'OGFS_FLX        ', 'convect-cld bot ' &   !68
                  ,'icetk           ', 'OGFS_FLX        ', 'sfc             ' &   !69
                  ,'soilw           ', 'OGFS_FLX        ', '40-100 cm down  ' &   !70
                  ,'soilw           ', 'OGFS_FLX        ', '100-200 cm down ' &   !71
                  ,'tmp             ', 'OGFS_FLX        ', '40-100 cm down  ' &   !72
                  ,'tmp             ', 'OGFS_FLX        ', '100-200 cm down ' &   !73
                  ,'soill           ', 'OGFS_FLX        ', '0-10 cm down    ' &   !74
                  ,'soill           ', 'OGFS_FLX        ', '10-40 cm down   ' &   !75
                  ,'soill           ', 'OGFS_FLX        ', '40-100 cm down  ' &   !76
                  ,'soill           ', 'OGFS_FLX        ', '100-200 cm down ' &   !77
                  ,'snod            ', 'OGFS_FLX        ', 'sfc             ' &   !78
                  ,'cnwat           ', 'OGFS_FLX        ', 'sfc             ' &   !79
                  ,'sfcr            ', 'OGFS_FLX        ', 'sfc             ' &   !80
                  ,'veg             ', 'OGFS_FLX        ', 'sfc             ' &   !81
                  ,'vgtyp           ', 'OGFS_FLX        ', 'sfc             ' &   !82
                  ,'sotyp           ', 'OGFS_FLX        ', 'sfc             ' &   !83
                  ,'sltyp           ', 'OGFS_FLX        ', 'sfc             ' &   !84
                  ,'fricv           ', 'OGFS_FLX        ', 'sfc             ' &   !85
                  ,'hgt             ', 'OGFS_FLX        ', 'sfc             ' &   !86
                  ,'crain           ', 'OGFS_FLX        ', 'sfc             ' &   !87
                  ,'sfexc           ', 'OGFS_FLX        ', 'sfc             ' &   !88
                  ,'acond           ', 'OGFS_FLX        ', 'sfc             ' &   !89
                  ,'pevpr           ', 'OGFS_FLX        ', 'sfc             ' &   !90
                  ,'dlwrf           ', 'OGFS_FLX        ', 'sfc             ' &   !91
                  ,'ulwrf           ', 'OGFS_FLX        ', 'sfc             ' &   !92
                  ,'uswrf           ', 'OGFS_FLX        ', 'sfc             ' &   !93
                  ,'dswrf           ', 'OGFS_FLX        ', 'sfc             ' &   !94
                  ,'shtfl           ', 'OGFS_FLX        ', 'sfc             ' &   !95
                  ,'lhtfl           ', 'OGFS_FLX        ', 'sfc             ' &   !96
                  ,'gflux           ', 'OGFS_FLX        ', 'sfc             ' &   !97
                  ,'ssrun_acc       ', 'OGFS_FLX        ', 'sfc             ' &   !98
                  ,'tmp             ', 'OGFS_FLX        ', 'hybrid lev 1    ' &   !99
                  ,'spfh            ', 'OGFS_FLX        ', 'hybrid lev 1    ' &   !100
                  ,'ugrd            ', 'OGFS_FLX        ', 'hybrid lev 1    ' &   !101
                  ,'vgrd            ', 'OGFS_FLX        ', 'hybrid lev 1    ' &   !102
                  ,'hgt             ', 'OGFS_FLX        ', 'hybrid lev 1    ' &   !103
                  ,'evbs_ave        ', 'OGFS_FLX        ', 'sfc             ' &   !104
                  ,'evcw_ave        ', 'OGFS_FLX        ', 'sfc             ' &   !105
                  ,'trans_ave       ', 'OGFS_FLX        ', 'sfc             ' &   !106
                  ,'sbsno_ave       ', 'OGFS_FLX        ', 'sfc             ' &   !107
                  ,'snowc_ave       ', 'OGFS_FLX        ', 'sfc             ' &   !108
                  ,'soilm           ', 'OGFS_FLX        ', '0-200 cm down   ' &   !109
                  ,'snohf           ', 'OGFS_FLX        ', 'sfc             ' &   !110
                  ,'wilt            ', 'OGFS_FLX        ', 'sfc             ' &   !111
                  ,'fldcp           ', 'OGFS_FLX        ', 'sfc             ' &   !112
                  ,'sunsd_acc       ', 'OGFS_FLX        ', 'sfc             ' &   !113
                  ,'cpofp           ', 'OGFS_FLX        ', 'sfc             ' &   !114

                 /)                                                           &
               ,(/3,MAX_KOUNT/)                                               &
               ,(/'*               ', '*               ','*               '/))


      CHARACTER(16),DIMENSION(3,MAX_KOUNT),TARGET :: PHY_INT_STATE_2D_R_FLX_CLM    &
!
       =RESHAPE((/ 'uflx_ave        ', 'OGFS_FLX        ', 'sfc             ' &   !1
                  ,'vflx_ave        ', 'OGFS_FLX        ', 'sfc             ' &   !2
                  ,'shtfl_ave       ', 'OGFS_FLX        ', 'sfc             ' &   !3
                  ,'lhtfl_ave       ', 'OGFS_FLX        ', 'sfc             ' &   !4
                  ,'tmp             ', 'OGFS_FLX        ', 'sfc             ' &   !5
                  ,'soilw           ', 'OGFS_FLX        ', '0-10 cm down    ' &   !6
                  ,'soilw           ', 'OGFS_FLX        ', '10-40 cm down   ' &   !7
                  ,'tmp             ', 'OGFS_FLX        ', '0-10 cm down    ' &   !8
                  ,'tmp             ', 'OGFS_FLX        ', '10-40 cm down   ' &   !9
                  ,'weasd           ', 'OGFS_FLX        ', 'sfc             ' &   !10
                  ,'ulwrf_ave       ', 'OGFS_FLX        ', 'nom. top        ' &   !11
                  ,'uswrf_ave       ', 'OGFS_FLX        ', 'nom. top        ' &   !12
                  ,'uswrf_ave       ', 'OGFS_FLX        ', 'sfc             ' &   !13
                  ,'dswrf_ave       ', 'OGFS_FLX        ', 'sfc             ' &   !14
                  ,'tcdc_ave        ', 'OGFS_FLX        ', 'high cld lay    ' &   !15
                  ,'pres_ave        ', 'OGFS_FLX        ', 'high cld top    ' &   !16
                  ,'pres_ave        ', 'OGFS_FLX        ', 'high cld bot    ' &   !17
                  ,'tmp_ave         ', 'OGFS_FLX        ', 'high cld top    ' &   !18
                  ,'tcdc_ave        ', 'OGFS_FLX        ', 'mid cld lay     ' &   !19
                  ,'pres_ave        ', 'OGFS_FLX        ', 'mid cld top     ' &   !20
                  ,'pres_ave        ', 'OGFS_FLX        ', 'mid cld bot     ' &   !21
                  ,'tmp_ave         ', 'OGFS_FLX        ', 'mid cld top     ' &   !22
                  ,'tcdc_ave        ', 'OGFS_FLX        ', 'low cld lay     ' &   !23
                  ,'pres_ave        ', 'OGFS_FLX        ', 'low cld top     ' &   !24
                  ,'pres_ave        ', 'OGFS_FLX        ', 'low cld bot     ' &   !25
                  ,'tmp_ave         ', 'OGFS_FLX        ', 'low cld top     ' &   !26
                  ,'tcdc_ave        ', 'OGFS_FLX        ', 'atmos col       ' &   !27
                  ,'tcdc_ave        ', 'OGFS_FLX        ', 'bndary-layer cld' &   !28
                  ,'dlwrf_ave       ', 'OGFS_FLX        ', 'sfc             ' &   !29
                  ,'ulwrf_ave       ', 'OGFS_FLX        ', 'sfc             ' &   !30
                  ,'duvb_ave        ', 'OGFS_FLX        ', 'sfc             ' &   !31
                  ,'cduvb_ave       ', 'OGFS_FLX        ', 'sfc             ' &   !32
                  ,'dswrf           ', 'OGFS_FLX        ', 'nom. top        ' &   !33
                  ,'vbdsf_ave       ', 'OGFS_FLX        ', 'sfc             ' &   !34
                  ,'vddsf_ave       ', 'OGFS_FLX        ', 'sfc             ' &   !35
                  ,'nbdsf_ave       ', 'OGFS_FLX        ', 'sfc             ' &   !36
                  ,'nddsf_ave       ', 'OGFS_FLX        ', 'sfc             ' &   !37
                  ,'csulf           ', 'OGFS_FLX        ', 'nom. top        ' &   !38
                  ,'csusf           ', 'OGFS_FLX        ', 'nom. top        ' &   !39
                  ,'csdlf           ', 'OGFS_FLX        ', 'sfc             ' &   !40
                  ,'csusf           ', 'OGFS_FLX        ', 'sfc             ' &   !41
                  ,'csdsf           ', 'OGFS_FLX        ', 'sfc             ' &   !42
                  ,'csulf           ', 'OGFS_FLX        ', 'sfc             ' &   !43
                  ,'albdo_ave       ', 'OGFS_FLX        ', 'sfc             ' &   !44
                  ,'prate_ave       ', 'OGFS_FLX        ', 'sfc             ' &   !45
                  ,'cprat_ave       ', 'OGFS_FLX        ', 'sfc             ' &   !46
                  ,'gflux_ave       ', 'OGFS_FLX        ', 'sfc             ' &   !47
                  ,'land            ', 'OGFS_FLX        ', 'sfc             ' &   !48
                  ,'icec            ', 'OGFS_FLX        ', 'sfc             ' &   !49
                  ,'ugrd            ', 'OGFS_FLX        ', '10 m above gnd  ' &   !50
                  ,'vgrd            ', 'OGFS_FLX        ', '10 m above gnd  ' &   !51
                  ,'tmp             ', 'OGFS_FLX        ', '2 m above gnd   ' &   !52
                  ,'spfh            ', 'OGFS_FLX        ', '2 m above gnd   ' &   !53
                  ,'pres            ', 'OGFS_FLX        ', 'sfc             ' &   !54
                  ,'tmax            ', 'OGFS_FLX        ', '2 m above gnd   ' &   !55
                  ,'tmin            ', 'OGFS_FLX        ', '2 m above gnd   ' &   !56
                  ,'spfhmax         ', 'OGFS_FLX        ', '2 m above gnd   ' &   !57
                  ,'spfhmin         ', 'OGFS_FLX        ', '2 m above gnd   ' &   !58
                  ,'watr_acc        ', 'OGFS_FLX        ', 'sfc             ' &   !59
                  ,'pevpr_ave       ', 'OGFS_FLX        ', 'sfc             ' &   !60
                  ,'cwork_ave       ', 'OGFS_FLX        ', 'atmos col       ' &   !61
                  ,'u-gwd_ave       ', 'OGFS_FLX        ', 'sfc             ' &   !62
                  ,'v-gwd_ave       ', 'OGFS_FLX        ', 'sfc             ' &   !63
                  ,'hpbl            ', 'OGFS_FLX        ', 'sfc             ' &   !64
                  ,'pwat            ', 'OGFS_FLX        ', 'atmos col       ' &   !65
                  ,'tcdc            ', 'OGFS_FLX        ', 'convect-cld laye' &   !66
                  ,'pres            ', 'OGFS_FLX        ', 'convect-cld top ' &   !67
                  ,'pres            ', 'OGFS_FLX        ', 'convect-cld bot ' &   !68
                  ,'icetk           ', 'OGFS_FLX        ', 'sfc             ' &   !69
                  ,'soilw           ', 'OGFS_FLX        ', '40-100 cm down  ' &   !70
                  ,'soilw           ', 'OGFS_FLX        ', '100-200 cm down ' &   !71
                  ,'tmp             ', 'OGFS_FLX        ', '40-100 cm down  ' &   !72
                  ,'tmp             ', 'OGFS_FLX        ', '100-200 cm down ' &   !73
                  ,'soill           ', 'OGFS_FLX        ', '0-10 cm down    ' &   !74
                  ,'soill           ', 'OGFS_FLX        ', '10-40 cm down   ' &   !75
                  ,'soill           ', 'OGFS_FLX        ', '40-100 cm down  ' &   !76
                  ,'soill           ', 'OGFS_FLX        ', '100-200 cm down ' &   !77
                  ,'snod            ', 'OGFS_FLX        ', 'sfc             ' &   !78
                  ,'cnwat           ', 'OGFS_FLX        ', 'sfc             ' &   !79
                  ,'sfcr            ', 'OGFS_FLX        ', 'sfc             ' &   !80
                  ,'veg             ', 'OGFS_FLX        ', 'sfc             ' &   !81
                  ,'vgtyp           ', 'OGFS_FLX        ', 'sfc             ' &   !82
                  ,'sotyp           ', 'OGFS_FLX        ', 'sfc             ' &   !83
                  ,'sltyp           ', 'OGFS_FLX        ', 'sfc             ' &   !84
                  ,'fricv           ', 'OGFS_FLX        ', 'sfc             ' &   !85
                  ,'hgt             ', 'OGFS_FLX        ', 'sfc             ' &   !86
                  ,'crain           ', 'OGFS_FLX        ', 'sfc             ' &   !87
                  ,'sfexc           ', 'OGFS_FLX        ', 'sfc             ' &   !88
                  ,'acond           ', 'OGFS_FLX        ', 'sfc             ' &   !89
                  ,'pevpr           ', 'OGFS_FLX        ', 'sfc             ' &   !90
! -skip 91-97 for climate mode-
                  ,'ssrun_acc       ', 'OGFS_FLX        ', 'sfc             ' &   !98
                  ,'tmp             ', 'OGFS_FLX        ', 'hybrid lev 1    ' &   !99
                  ,'spfh            ', 'OGFS_FLX        ', 'hybrid lev 1    ' &   !100
                  ,'ugrd            ', 'OGFS_FLX        ', 'hybrid lev 1    ' &   !101
                  ,'vgrd            ', 'OGFS_FLX        ', 'hybrid lev 1    ' &   !102
                  ,'hgt             ', 'OGFS_FLX        ', 'hybrid lev 1    ' &   !103
                  ,'evbs_ave        ', 'OGFS_FLX        ', 'sfc             ' &   !104
                  ,'evcw_ave        ', 'OGFS_FLX        ', 'sfc             ' &   !105
                  ,'trans_ave       ', 'OGFS_FLX        ', 'sfc             ' &   !106
                  ,'sbsno_ave       ', 'OGFS_FLX        ', 'sfc             ' &   !107
                  ,'snowc_ave       ', 'OGFS_FLX        ', 'sfc             ' &   !108
                  ,'soilm           ', 'OGFS_FLX        ', '0-200 cm down   ' &   !109
                  ,'snohf           ', 'OGFS_FLX        ', 'sfc             ' &   !110
                  ,'wilt            ', 'OGFS_FLX        ', 'sfc             ' &   !111
                  ,'fldcp           ', 'OGFS_FLX        ', 'sfc             ' &   !112
                  ,'sunsd_acc       ', 'OGFS_FLX        ', 'sfc             ' &   !113
                  ,'cpofp           ', 'OGFS_FLX        ', 'sfc             ' &   !114

                 /)                                                           &
               ,(/3,MAX_KOUNT/)                                               &
               ,(/'*               ', '*               ','*               '/))
!
!-------------------------
!***  REAL 2-D ARRAYS aer***                                     (for g2d_fld)
!-------------------------
!
      CHARACTER(16),DIMENSION(3,MAX_KOUNT),TARGET :: PHY_INT_STATE_2D_R_AER
!
!-------------------------
!***  REAL 2-D ARRAYS  nst**
!-------------------------
!
      CHARACTER(16),DIMENSION(3,MAX_KOUNT),TARGET :: PHY_INT_STATE_2D_R_NST        &
!
       =RESHAPE((/ 'slmsk           ', 'OGFS_NST        ', 'sfc             ' &   !1
                  ,'xt              ', 'OGFS_NST        ', 'sfc             ' &   !2
                  ,'xs              ', 'OGFS_NST        ', 'sfc             ' &   !3
                  ,'xu              ', 'OGFS_NST        ', 'sfc             ' &   !4
                  ,'xv              ', 'OGFS_NST        ', 'sfc             ' &   !5
                  ,'xz              ', 'OGFS_NST        ', 'sfc             ' &   !6
                  ,'zm              ', 'OGFS_NST        ', 'sfc             ' &   !7
                  ,'xtts            ', 'OGFS_NST        ', 'sfc             ' &   !8
                  ,'xzts            ', 'OGFS_NST        ', 'sfc             ' &   !9
                  ,'dtCool          ', 'OGFS_NST        ', 'sfc             ' &   !10
                  ,'zc              ', 'OGFS_NST        ', 'sfc             ' &   !11
                  ,'c0              ', 'OGFS_NST        ', 'sfc             ' &   !12
                  ,'cd              ', 'OGFS_NST        ', 'sfc             ' &   !13
                  ,'w0              ', 'OGFS_NST        ', 'sfc             ' &   !14
                  ,'wd              ', 'OGFS_NST        ', 'sfc             ' &   !15
                  ,'dconv           ', 'OGFS_NST        ', 'sfc             ' &   !16
                  ,'ifd             ', 'OGFS_NST        ', 'sfc             ' &   !17
                  ,'tref            ', 'OGFS_NST        ', 'sfc             ' &   !18
                  ,'qrain           ', 'OGFS_NST        ', 'sfc             ' &   !19
                 /)                                                           &
               ,(/3,MAX_KOUNT/)                                               &
               ,(/'*               ', '*               ','*               '/))
!
!-----------------------------------------------------------------------
!***  FIRST WE MUST PROVIDE POINTERS INTO THE PHYSICS INTERNAL STATE.
!-----------------------------------------------------------------------

      TYPE PHY_ISC
        INTEGER(KIND=kind_io4),                POINTER :: NAME    !<-- Pointer for integer scalars
      END TYPE PHY_ISC
!-----------------------------------------------------------------------
      TYPE PHY_RSC
        REAL(KIND=kind_grid),                  POINTER :: NAME    !<-- Pointer for real scalars
      END TYPE PHY_RSC
!-----------------------------------------------------------------------
      TYPE PHY_I1D
        INTEGER(KIND=kind_io4),DIMENSION(:),   POINTER :: NAME    !<-- Pointer for 1D integer arrays
      END TYPE PHY_I1D
!-----------------------------------------------------------------------
      TYPE PHY_I2D
        INTEGER(KIND=kind_io4),DIMENSION(:,:), POINTER :: NAME    !<-- Pointer for 2D integer arrays
      END TYPE PHY_I2D
!-----------------------------------------------------------------------
      TYPE PHY_R1D
        REAL(KIND=kind_io4),DIMENSION(:),      POINTER :: NAME    !<-- Pointer for 1D real arrays
      END TYPE PHY_R1D
!-----------------------------------------------------------------------
      TYPE PHY_R2D
        REAL(KIND=kind_io4),DIMENSION(:,:),    POINTER :: NAME    !<-- Pointer for 2D real arrays
      END TYPE PHY_R2D
!-----------------------------------------------------------------------
      TYPE PHY_R3D
        REAL(KIND=kind_io4),DIMENSION(:,:,:),  POINTER :: NAME    !<-- Pointer for 3D real arrays
      END TYPE PHY_R3D
!-----------------------------------------------------------------------
      TYPE PHY_R4D
        REAL(KIND=kind_io4),DIMENSION(:,:,:,:),POINTER :: NAME    !<-- Pointer for 4D real arrays
      END TYPE PHY_R4D
!-----------------------------------------------------------------------
!
      CONTAINS
!
!-----------------------------------------------------------------------
!
      SUBROUTINE POINT_PHYSICS_OUTPUT_GFS(INT_STATE,IMP_WRITE_STATE)
!
!-----------------------------------------------------------------------
!
       use machine,              only: kind_io4
       use resol_def,            only: thermodyn_id,sfcpress_id,ntrac,ntke, &
                                       ngrids_gg,ivssfc,ivsnst
       use coordinate_def,       only: vertcoord_id,AK5,BK5,CK5
       use date_def,             only: fhour,idate,zhour
       use namelist_physics_def, only: nstf_name,ivegsrc,cu_physics
!
      use mod_state, only:  buff_mult_piecea2d                & !for sfc 2d file
                           ,buff_mult_piecea3d                & !for sfc 3d file
                           ,buff_mult_piecef                  & !for flx file
                           ,buff_mult_pieceg                  & !for aer file
                           ,buff_mult_piecenst                  !for nst file

!-----------------------------------------------------------------------
!***  THIS ROUTINE TAKES THE USER'S SELECTIONS FOR OUTPUT QUANTITIES,
!***  POINTS AT THEM, AND INSERTS THOSE POINTERS INTO THE IMPORT STATE
!***  OF THE WRITE COMPONENTS.
!-----------------------------------------------------------------------
!
      TYPE(ESMF_State), INTENT(INOUT)          :: IMP_WRITE_STATE  !<-- Import state for the write gridded components
      TYPE(gfs_physics_internal_state),POINTER :: INT_STATE        !<-- The physics internal state
!
!-----------------------------------------------------------------------
!***  LOCAL VARIABLES
!-----------------------------------------------------------------------
!
      INTEGER                     :: mype,RC,RC_PHY_OUT
      INTEGER                     :: MINIDX(2),MAXIDX(2)
      INTEGER                     :: II, IJ, K,n,nst, ndg
      INTEGER                     :: GFS_BUNDLE_LIST_size, GFS_BUNDLE_LIST_aer, GFS_BUNDLE_LIST_nst
      INTEGER,allocatable         :: I2(:,:)
!     real(kind=kind_io4),target  :: fhourtmp
      real(kind=kind_io4),dimension(:),allocatable,target  :: r4tmp
!
!-----------------------------------------------------------------------
!***  ARRAYS OF POINTERS OF THE ABOVE TYPES
!-----------------------------------------------------------------------
!
! Add the ESMF field bundle list array to 
! avoid the ESMF log error.  WY.
! GFS_BUNDLE_LIST(1) -------- GFS_SFC_BUNDLE.
! GFS_BUNDLE_LIST(2) -------- GFS_FLX_BUNDLE.
! GFS_BUNDLE_LIST(GFS_BUNDLE_LIST_aer) -------- GFS_AER_BUNDLE.
! GFS_BUNDLE_LIST(GFS_BUNDLE_LIST_nst) -------- GFS_NST_BUNDLE.
!--------------------------------------------
      TYPE(ESMF_FieldBundle),ALLOCATABLE, SAVE :: GFS_BUNDLE_LIST(:)

      TYPE(ESMF_DISTGRID) :: DistGrid
      TYPE(ESMF_GRID)     :: Grid
      TYPE(ESMF_Info)     :: info
!
      TYPE(PHY_ISC),DIMENSION(MAX_KOUNT) :: I_SC
      TYPE(PHY_RSC),DIMENSION(MAX_KOUNT) :: R_SC
      TYPE(PHY_I1D),DIMENSION(MAX_KOUNT) :: I_1D
      TYPE(PHY_I2D),DIMENSION(MAX_KOUNT) :: I_2D
      TYPE(PHY_R1D),DIMENSION(MAX_KOUNT) :: R_1D
      TYPE(PHY_R2D),DIMENSION(MAX_KOUNT) :: R_2D
      TYPE(PHY_R3D),DIMENSION(MAX_KOUNT) :: R_3D
      TYPE(PHY_R4D),DIMENSION(MAX_KOUNT) :: R_4D
!
!-----------------------------------------------------------------------
!
      MYPE = int_state%ME
!
!-----------------------------------------------------------------------
!*** create dist  grid
!-----------------------------------------------------------------------
!*** fortran array buff_mult_piecea2d is on dimension (buff_mult_piece3d
!*** has same dimension as buff_mult_piece3d )
!*** (int_state%lonf,int_state%lats_node_a_max)
!
      minidx(1) = 1
      minidx(2) = 1
      maxidx(1) = size(buff_mult_piecea2d,1)*int_state%nodes
      maxidx(2) = size(buff_mult_piecea2d,2)

!     write(0,*)'in phys_output,lonr=',int_state%lonr,'lats_node_a_max=',int_state%lats_node_r_max, &
!       'minidx=',minidx,'maxidx=',maxidx

! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK = "Create DISGRID for write Bundle"
      CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      DistGrid = ESMF_DistGridCreate(minidx, maxidx, rc=rc)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL gfs_physics_err_msg(RC,MESSAGE_CHECK,RC_PHY_OUT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!jw double check the grid index with buff_mult_pieceg
      allocate(i2(2,int_state%nodes))
      i2 = 0
      CALL ESMF_DistGridGet(DistGrid, indexCountPDe=i2, rc=rc)
!
!-----------------------------------------------------------------------
!*** create grid
!-----------------------------------------------------------------------
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK = "Create GRID from DistGrid for write Bundle"
      CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      grid = ESMF_GridCreate(name="gridwrt", distgrid=DistGrid, rc=rc)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL gfs_physics_err_msg(RC,MESSAGE_CHECK,RC_PHY_OUT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!Set up GFS_BUNDLE_LIST_size, allocate GFS_BUNDLE_LIST.
!------------------------------------------------------
      GFS_BUNDLE_LIST_size = 2
      IF(int_state%lgocart) THEN
          GFS_BUNDLE_LIST_size = GFS_BUNDLE_LIST_size + 1
          GFS_BUNDLE_LIST_aer  = 3
          IF(nstf_name(1) > 0) THEN
              GFS_BUNDLE_LIST_nst = 4
              GFS_BUNDLE_LIST_size = GFS_BUNDLE_LIST_size + 1
          END IF
       ELSE
           IF(nstf_name(1) > 0) THEN
               GFS_BUNDLE_LIST_size = GFS_BUNDLE_LIST_size + 1
               GFS_BUNDLE_LIST_nst  = 3
           END IF
       END IF
       ALLOCATE(GFS_BUNDLE_LIST(GFS_BUNDLE_LIST_size))
!
!-----------------------------------------------------------------------
!***  THE PHYSICS INTERNAL STATE POINTER BLOCK
!-----------------------------------------------------------------------
!***  POINT AT ALL VARIABLES IN THE PHYSICS INTERNAL STATE 
!***  THAT COULD BE WRITTEN TO HISTORY OUTPUT, I.E., THOSE
!***  LISTED AT THE TOP OF THIS MODULE.
!-----------------------------------------------------------------------
!
!---------------------
!***  INTEGER SCALARS
!---------------------
!
      I_SC( 1)%NAME => int_state%latr
      I_SC( 2)%NAME => int_state%lonr
      I_SC( 3)%NAME => int_state%levs
      I_SC( 4)%NAME => int_state%ntoz
      I_SC( 5)%NAME => int_state%ntcw
      I_SC( 6)%NAME => int_state%ncld
      I_SC( 7)%NAME => int_state%ntke
      I_SC( 8)%NAME => int_state%ntrac
      I_SC( 9)%NAME => thermodyn_id
      I_SC(10)%NAME => sfcpress_id
      I_SC(11)%NAME => int_state%lsoil
      I_SC(12)%NAME => int_state%idrt
      I_SC(13)%NAME => ivssfc
      I_SC(14)%NAME => ivegsrc
      I_SC(15)%NAME => cu_physics
      I_SC(16)%NAME => ivsnst
!        
!------------------
!***  REAL SCALARS
!------------------
!
      R_SC(1)%NAME => fhour
      zhour = int_state%zhour
      R_SC(2)%NAME => zhour
!        
!-----------------------
!***  1D INTEGER ARRAYS
!-----------------------
!
      I_1D(1)%NAME => idate
      I_1D(2)%NAME => int_state%lonsperlar
!        
!--------------------
!***  1D REAL ARRAYS
!--------------------
!
      R_1D(1)%NAME => int_state%zsoil
!     R_1D(2)%NAME => bk5
!     R_1D(3)%NAME => ck5
!        
!-----------------------
!***  2D INTEGER ARRAYS
!-----------------------
!
!--------------------
!***  2D REAL ARRAYS, save in 3d array buff_mult_piecea2D
!--------------------
!
!jw   R_2D( 1)%NAME => buff_mult_piecea2D
!        
!--------------------
!***  3D REAL ARRAYS
!--------------------
!
!      write(0,*)'buff_mult_piecea2D=',maxval(buff_mult_piecea2D), &
!       minval(buff_mult_piecea2D),'buff_mult_piecea3D=',   &
!       maxval(buff_mult_piecea3d),minval(buff_mult_piecea3D)

      R_3D(1)%NAME => buff_mult_piecea2D
      R_3D(2)%NAME => buff_mult_piecea3D
!flx
      R_3D(3)%NAME => buff_mult_piecef
!aer
      R_3D(4)%NAME => buff_mult_pieceg                 ! for g2d_fld
!nst
      R_3D(5)%NAME => buff_mult_piecenst
!
!
!-----------------------------------------------------------------------
!##jw
!***  CREATE AN ESMF Bundle THAT WILL HOLD HISTORY OUTPUT DATA
!***  AND NOTHING ELSE.  THIS WILL SERVE TO ISOLATE THE OUTPUT
!***  DATA FROM EVERYTHING ELSE INSIDE THE WRITE COMPONENT'S
!***  IMPORT STATE.
!-----------------------------------------------------------------------
!***  SINCE PHY output separate file, do not need to:
!***  EXTRACT THE HISTORY OUTPUT Bundle FROM THE WRITE COMPONENT'S
!***  IMPORT STATE.  IT ALREADY CONTAINS OUTPUT VARIABLES FROM
!***  THE DYNAMICS.  WE ARE PREPARING TO ADD HISTORY VARIABLES
!***  FROM THE PHYSICS.
!-----------------------------------------------------------------------
!! for sfc file
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK = "Extract History Data Bundle from the Write Import State"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      GFS_BUNDLE_LIST(1) = ESMF_FieldBundleCreate(name=int_state%filename_base(2) &  !<-- The Bundle's name
                                           ,rc  =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL gfs_physics_err_msg(RC,MESSAGE_CHECK,RC_PHY_OUT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!*** add im, jm into import state
!
      CALL ESMF_InfoGetFromHost(IMP_WRITE_STATE, info, rc=RC)

      CALL ESMF_InfoSet(info                                            &  !<-- The Write component import state's info handle
                        ,key      ='im'                                 &  !<-- Name of the integer attribute
                        ,value    =int_state%lonr                       &  !<-- The value being inserted into the import state
                        ,rc       =RC)

      CALL ESMF_InfoSet(info                                            &  !<-- The Write component import state's info handle
                        ,key      ='jm'                                 &  !<-- Name of the integer attribute
                        ,value    =int_state%latr                       &  !<-- The value being inserted into the import state
                        ,rc       =RC)
!
      CALL ESMF_InfoSet(info                                            &  !<-- The Write component import state's info handle
                        ,key      ='zhour'                              &  !<-- Name of the integer attribute
                        ,value    =int_state%zhour                      &  !<-- The value being inserted into the import state
                        ,rc       =RC)
!
!-----------------------------------------------------------------------
!***  INSERT INTO THE WRITE COMPONENTS' IMPORT STATE THE POINTERS
!***  OF ONLY THOSE QUANTITIES THAT ARE SPECIFIED BY THE USER
!***  FOR HISTORY OUTPUT.
!***  PUT THE HISTORY DATA INTO THE ESMF Bundle RESIDING IN THE
!***  WRITE COMPONENT'S IMPORT STATE.
!-----------------------------------------------------------------------
!
      CALL ADD_BUNDLE_TO_WRITE(int_state,  imp_write_state,   &
                               grid, I_SC, R_SC, I_1D, I_2D,  &
                               R_1D, R_2D, R_3D, R_4D,        &
                               2,'OGFS_PHY','OGFS_SFC',GFS_BUNDLE_LIST(1))
!
!-----------------------------------------------------------------------
!! for flx file
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK = "Extract History Data Bundle from the Write Import State"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      GFS_BUNDLE_LIST(2) = ESMF_FieldBundleCreate(name=int_state%filename_base(3) &  !<-- The Bundle's name
                                           ,rc  =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL gfs_physics_err_msg(RC,MESSAGE_CHECK,RC_PHY_OUT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
! 
!-----------------------------------------------------------------------
!***  INSERT INTO THE WRITE COMPONENTS' IMPORT STATE THE POINTERS
!***  OF ONLY THOSE QUANTITIES THAT ARE SPECIFIED BY THE USER
!***  FOR HISTORY OUTPUT.
!***  PUT THE HISTORY DATA INTO THE ESMF Bundle RESIDING IN THE
!***  WRITE COMPONENT'S IMPORT STATE.
!-----------------------------------------------------------------------
!
      CALL ADD_BUNDLE_TO_WRITE(int_state,  imp_write_state,   &
                               grid, I_SC, R_SC, I_1D, I_2D,  &
                               R_1D, R_2D, R_3D, R_4D,        &
                               2,'OGFS_PHY','OGFS_FLX',GFS_BUNDLE_LIST(2))
!-----------------------------------------------------------------------
!! for aer file (optional)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

! SARAH: this is the logic from the trunk version (num_file is hardwired to 4)
!      lab_if_gocart : IF ( int_state%lgocart .and.                     &
!                          (int_state%num_file .eq. 4) ) THEN
! SARAH: now I change it to be consisten with Jun's nst revision

       lab_if_gocart : IF ( int_state%lgocart ) THEN
        ndg = 0
        do n=1,int_state%num_file
         if (trim(int_state%filename_base(n))=='AER.F' ) ndg = n
        enddo
        if (ndg == 0) then
         write(0,*)"ERROR!!! gocart forecast without aer_diag output! "
         stop 555
        endif


! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
! set up PHY_INT_STATE_2D_R_AER
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!  fill in default values
      DO II =1,MAX_KOUNT
        PHY_INT_STATE_2D_R_AER(1:3,II) = '**********'
      ENDDO

      II = 0
!  fill in variable name
      if ( int_state%g2d_fld%du%nfld > 0 ) then   ! fill in du diag
        do  k = 1, int_state%g2d_fld%du%nfld
         II = II + 1
         PHY_INT_STATE_2D_R_AER(1,II) = int_state%g2d_fld%du%diag(k)%name
        enddo
      endif

      if ( int_state%g2d_fld%su%nfld > 0 ) then   ! fill in su diag
        do  k = 1, int_state%g2d_fld%su%nfld
         II = II + 1
         PHY_INT_STATE_2D_R_AER(1,II) = int_state%g2d_fld%su%diag(k)%name
        enddo
      endif

      if ( int_state%g2d_fld%ss%nfld > 0 ) then   ! fill in ss diag
        do  k = 1, int_state%g2d_fld%ss%nfld
         II = II + 1
         PHY_INT_STATE_2D_R_AER(1,II) = int_state%g2d_fld%ss%diag(k)%name
        enddo
      endif

      if ( int_state%g2d_fld%oc%nfld > 0 ) then   ! fill in oc diag
        do  k = 1, int_state%g2d_fld%oc%nfld
         II = II + 1
         PHY_INT_STATE_2D_R_AER(1,II) = int_state%g2d_fld%oc%diag(k)%name
        enddo
      endif

      if ( int_state%g2d_fld%bc%nfld > 0 ) then   ! fill in bc diag
        do  k = 1, int_state%g2d_fld%bc%nfld
         II = II + 1
         PHY_INT_STATE_2D_R_AER(1,II) = int_state%g2d_fld%bc%diag(k)%name
        enddo
      endif

      if ( int_state%g2d_fld%met%nfld > 0 ) then   ! fill in met diag
        do  k = 1, int_state%g2d_fld%met%nfld
         II = II + 1
         PHY_INT_STATE_2D_R_AER(1,II) = int_state%g2d_fld%met%diag(k)%name
        enddo
      endif

!  fill in file type and level type
      DO K = 1, II
        PHY_INT_STATE_2D_R_AER(2,K) = 'OGFS_AER  '
        PHY_INT_STATE_2D_R_AER(3,K) = 'atmos col '
      ENDDO
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK = "Extract History Data Bundle from the Write Import State"
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      GFS_BUNDLE_LIST(GFS_BUNDLE_LIST_aer) = ESMF_FieldBundleCreate(          &
!                                         name=int_state%filename_base(4),    &  !<-- The Bundle's name
                                          name=int_state%filename_base(ndg),  &  !<-- The Bundle's name
                                          rc  =RC)
!
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL gfs_physics_err_msg(RC,MESSAGE_CHECK,RC_PHY_OUT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  INSERT INTO THE WRITE COMPONENTS' IMPORT STATE THE POINTERS
!***  OF ONLY THOSE QUANTITIES THAT ARE SPECIFIED BY THE USER
!***  FOR HISTORY OUTPUT.
!***  PUT THE HISTORY DATA INTO THE ESMF Bundle RESIDING IN THE
!***  WRITE COMPONENT'S IMPORT STATE.
!-----------------------------------------------------------------------
!
      CALL ADD_BUNDLE_TO_WRITE(int_state,  imp_write_state,   &
                               grid, I_SC, R_SC, I_1D, I_2D,  &
                               R_1D, R_2D, R_3D, R_4D,        &
                               2,'OGFS_PHY','OGFS_AER',GFS_BUNDLE_LIST(GFS_BUNDLE_LIST_aer))

      ENDIF lab_if_gocart
!
!
!-----------------------------------------------------------------------
!! for nst file (optional)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      nstf_name_if:    IF ( nstf_name(1) > 0)  THEN
!
       nst = 0
       do n=1,int_state%num_file
        if (trim(int_state%filename_base(n)) == 'NST.F' ) nst = n
       enddo
       if(nst == 0) then
         write(0,*)"WRONG!!! nst forecast without nst output! "
       endif
!       print *,'in phys_output,nst=',nst,'filename_base=',trim(int_state%filename_base(nst))
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK = "Extract History Data Bundle from the Write Import State"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      GFS_BUNDLE_LIST(GFS_BUNDLE_LIST_nst) = ESMF_FieldBundleCreate(         &
                                           name=int_state%filename_base(nst) &  !<-- The Bundle's name
                                           ,rc  =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL gfs_physics_err_msg(RC,MESSAGE_CHECK,RC_PHY_OUT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ADD_BUNDLE_TO_WRITE(int_state,  imp_write_state,     &
                               grid, I_SC, R_SC, I_1D, I_2D,    &
                               R_1D, R_2D, R_3D, R_4D,          &
                               2,'OGFS_PHY','OGFS_NST',GFS_BUNDLE_LIST(GFS_BUNDLE_LIST_nst))

!
      endif nstf_name_if

!-----------------------------------------------------------------------
!***  INSERT THE OUTPUT DATA Bundle INTO THE WRITE COMPONENT'S IMPORT STATE.
!-----------------------------------------------------------------------
!
      MESSAGE_CHECK = "Phsyics: Insert History Bundle into the Write Import State"
      CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_StateAddReplace(state    = IMP_WRITE_STATE              &
                        ,fieldbundleList = GFS_BUNDLE_LIST              &
                        ,rc              = RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL gfs_physics_err_msg(RC,MESSAGE_CHECK,RC_PHY_OUT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!
      END SUBROUTINE POINT_PHYSICS_OUTPUT_GFS
!
!-----------------------------------------------------------------------
!
      SUBROUTINE ADD_BUNDLE_TO_WRITE(int_state, imp_write_state,        &
                                     GRID, I_SC, R_SC, I_1D, I_2D,      &
                                     R_1D, R_2D, R_3D, R_4D,            &
                                     file_index,file_gen,file_id,file_bundle)
!
!-----------------------------------------------------------------------
!
      type(gfs_physics_internal_state),intent(in) :: int_state
      type(ESMF_STATE),intent(inout)              :: imp_write_state
      type(ESMF_GRID),intent(in)                  :: GRID
      character(*),intent(in)                     :: file_id,file_gen
      integer,intent(in)                          :: file_index
      type(ESMF_FieldBUNDLE),intent(inout)        :: file_bundle
      TYPE(PHY_ISC),intent(in)                    :: I_SC(:)
      TYPE(PHY_RSC),intent(in)                    :: R_SC(:)
      TYPE(PHY_I1D),intent(in)                    :: I_1D(:)
      TYPE(PHY_I2D),intent(in)                    :: I_2D(:)
      TYPE(PHY_R1D),intent(in)                    :: R_1D(:)
      TYPE(PHY_R2D),intent(in)                    :: R_2D(:)
      TYPE(PHY_R3D),intent(in)                    :: R_3D(:)
      TYPE(PHY_R4D),intent(in)                    :: R_4D(:)
!
!--- local variables
!
      INTEGER                      :: INDEX_IKJ,K,LENGTH                &
                                     ,MP_PHYSICS,MYPE,N,NDIM3,NFIND     &
                                     ,NUM_2D_FIELDS_I,NUM_2D_FIELDS_R   &
                                     ,RC,RC_PHY_OUT,NLEVS               &
                                     ,LDIM1,LDIM2,LDIM3,LDIM4           &
                                     ,UDIM1,UDIM2,UDIM3,UDIM4,NDIM4
!
      INTEGER(ESMF_KIND_I4),DIMENSION(:,:),POINTER :: TEMP_I2D
      REAL(ESMF_KIND_R4)   ,DIMENSION(:,:),POINTER :: TEMP_R2D
      CHARACTER(16),DIMENSION(:,:),POINTER :: PHY_INT_STATE_2D_R
!
      CHARACTER(2)             :: MODEL_LEVEL,TRACERS_KIND
      CHARACTER(6)             :: FMT='(I2.2)'
      CHARACTER(ESMF_MAXSTR)   :: VBL_NAME,VBL_NAME_X
!
      TYPE(ESMF_Field)         :: FIELD
      TYPE(ESMF_DataCopy_Flag) :: COPYFLAG=ESMF_DATACOPY_REFERENCE
      TYPE(ESMF_Info)          :: info

      MYPE = int_state%ME
!
!-----------------------------------------------------------------------
!***  THE MICROPHYSICS SCHEME SPECIFICATION IS NEEDED IN THE OUTPUT

!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!***  INSERT INTO THE WRITE COMPONENTS' IMPORT STATE THE POINTERS
!***  OF ONLY THOSE QUANTITIES THAT ARE SPECIFIED BY THE USER
!***  FOR HISTORY OUTPUT.
!***  PUT THE HISTORY DATA INTO THE ESMF Bundle RESIDING IN THE
!***  WRITE COMPONENT'S IMPORT STATE.
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!***  BEGIN WITH THE INTEGER SCALARS.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK = "Insert integer scalars into file bundle"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_InfoGetFromHost(file_bundle, info, rc=RC)

      DO NFIND=1,MAX_KOUNT
        IF(trim(PHY_INT_STATE_ISCALAR(file_index,NFIND)) == trim(file_id).or. &
           trim(PHY_INT_STATE_ISCALAR(file_index,NFIND)) == trim(file_gen) ) THEN !<-- Take integer scalar data specified for history output
          VBL_NAME = TRIM(PHY_INT_STATE_ISCALAR(1,NFIND))

!     write(0,*)' VBL_NAME =',VBL_NAME,' nfind=',nfind,' vali=',I_SC(NFIND)%NAME,' me=',mype

          CALL ESMF_InfoSet(info                                        &  !<-- The info handle for the Write component output history Bundle
                            ,key        =VBL_NAME                       &  !<-- Name of the integer scalar
                            ,value      =I_SC(NFIND)%NAME               &  !<-- The scalar being inserted into the output data Bundle
                            ,rc         =RC)
        ENDIF
!
        IF(PHY_INT_STATE_ISCALAR(file_index,NFIND)=='**********' ) THEN           !<-- End of the integer scalar list
          EXIT
        ENDIF
!
      ENDDO
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL gfs_physics_err_msg(RC,MESSAGE_CHECK,RC_PHY_OUT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  THE REAL SCALARS.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK = "Insert Physics Real Scalars into History Bundle"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      DO NFIND=1,MAX_KOUNT
        IF(trim(PHY_INT_STATE_RSCALAR(file_index,NFIND)) == trim(file_id).or. &
           trim(PHY_INT_STATE_RSCALAR(file_index,NFIND)) == trim(file_gen) ) THEN    !<-- Take real scalar data specified for history output
          VBL_NAME = TRIM(PHY_INT_STATE_RSCALAR(1,NFIND))
!     write(0,*)' VBL_NAME =',VBL_NAME,' nfind=',nfind,' val=',R_SC(NFIND)%NAME,' me=',mype

          CALL ESMF_InfoSet(info                                        &  !<-- The info handle for the Write component output history Bundle
                            ,key        =VBL_NAME                       &  !<-- Name of the integer scalar
                            ,value      =R_SC(NFIND)%NAME               &  !<-- The scalar being inserted into the output data Bundle
                            ,rc         =RC)

        ENDIF
!
        IF(PHY_INT_STATE_RSCALAR(file_index,NFIND)=='**********' ) THEN    !<-- End of the real scalar list
          EXIT
        ENDIF
!
      ENDDO
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL gfs_physics_err_msg(RC,MESSAGE_CHECK,RC_PHY_OUT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  THE 1D INTEGER ARRAYS
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK = "Insert Physics 1-D Integer Arrays into History Bundle"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      DO NFIND=1,MAX_KOUNT
        IF(trim(PHY_INT_STATE_1D_I(file_index,NFIND)) == trim(file_id).or.  &
           trim(PHY_INT_STATE_1D_I(file_index,NFIND)) == trim(file_gen) ) THEN       !<-- Take 1D integer array data specified for history output
          VBL_NAME = TRIM(PHY_INT_STATE_1D_I(1,NFIND))
          LENGTH   = SIZE(I_1D(NFIND)%NAME)
          if(VBL_NAME == "LPL") LENGTH = SIZE(I_1D(NFIND)%NAME)/2
!
!         write(0,*)'phy_out,I1D NFIND=',NFIND,'VBL_NAME=',trim(VBL_NAME),  &
!            'length=',length,'value=',I_1D(NFIND)%NAME

          CALL ESMF_InfoSet(info                                        &  !<-- The info handle for the Write component output history Bundle
                            ,key        =VBL_NAME                       &  !<-- Name of the integer array
                            ,values     =I_1D(NFIND)%NAME               &  !<-- The 1D integer being inserted into the output data Bundle
                            ,rc         =RC)

        ENDIF
!

        IF(PHY_INT_STATE_1D_I(file_index,NFIND)=='**********' ) THEN   !<-- End of the 1D integer array list
          EXIT
        ENDIF
!
      ENDDO
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL gfs_physics_err_msg(RC,MESSAGE_CHECK,RC_PHY_OUT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  THE 1D REAL ARRAYS.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK = "Insert Physics 1-D Real Arrays into History Bundle"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      DO NFIND=1,MAX_KOUNT
        IF(trim(PHY_INT_STATE_1D_R(file_index,NFIND)) == trim(file_id).or. &
           trim(PHY_INT_STATE_1D_R(file_index,NFIND)) == trim(file_gen) ) THEN  !<-- Take 1D real array data specified for history output
          VBL_NAME = TRIM(PHY_INT_STATE_1D_R(1,NFIND))
          LENGTH   = SIZE(R_1D(NFIND)%NAME)
!
!          write(0,*)'phy_out,R1D NFIND=',NFIND,'VBL_NAME=',trim(VBL_NAME),  &
!             'length=',length,'value=',R_1D(NFIND)%NAME

          CALL ESMF_InfoSet(info                                        &  !<-- The info handle for the Write component output history Bundle
                            ,key        =VBL_NAME                       &  !<-- Name of the real array
                            ,values     =R_1D(NFIND)%NAME               &  !<-- The 1D real being inserted into the output data Bundle
                            ,rc         =RC)

        ENDIF
!
        IF(PHY_INT_STATE_1D_R(file_index,NFIND)=='**********' ) THEN                !<-- End of the 1D real array list
          EXIT
        ENDIF
!
      ENDDO
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL gfs_physics_err_msg(RC,MESSAGE_CHECK,RC_PHY_OUT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  THE 2D INTEGER ARRAYS.
!-----------------------------------------------------------------------
!
      NUM_2D_FIELDS_I = 0
      NULLIFY(TEMP_I2D)
!
      DO NFIND=1,MAX_KOUNT
        IF(trim(PHY_INT_STATE_2D_I(file_index,NFIND)) == trim(file_id)) THEN    !<-- Take 2D integer array data specified for history output
          LDIM1 = LBOUND(I_2D(NFIND)%NAME,1)
          UDIM1 = UBOUND(I_2D(NFIND)%NAME,1)
          LDIM2 = LBOUND(I_2D(NFIND)%NAME,2)
          UDIM2 = UBOUND(I_2D(NFIND)%NAME,2)
!
          VBL_NAME = TRIM(PHY_INT_STATE_2D_I(1,NFIND))
          TEMP_I2D => I_2D(NFIND)%NAME(LDIM1:UDIM1,LDIM2:UDIM2)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!          write(0,*)'phy_out,I2D NFIND=',NFIND,'VBL_NAME=',trim(VBL_NAME),  &
!             'value=',maxval(I_2D(NFIND)%NAME(LDIM1:UDIM1,LDIM2:UDIM2)),    &
!              minval(I_2D(NFIND)%NAME(LDIM1:UDIM1,LDIM2:UDIM2))
!
          MESSAGE_CHECK = "Insert Physics 2-D Integer Data into Field"
!         CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          FIELD=ESMF_FieldCreate(grid         =GRID                     &  !<-- The ESMF grid
                                ,farray       =TEMP_I2D                 &  !<-- The 2D integer array being inserted into output data Bundle
                                ,datacopyflag =COPYFLAG                 &
                                ,name         =VBL_NAME                 &  !<-- Name of the 2D integer array
                                ,indexFlag=ESMF_INDEX_DELOCAL           &
                                ,rc           =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          CALL gfs_physics_err_msg(RC,MESSAGE_CHECK,RC_PHY_OUT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          MESSAGE_CHECK = "Insert Physics 2-D Integer Field into History Bundle"
!         CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          CALL ESMF_FieldBundleAdd(FILE_BUNDLE  &  !<-- The write component's output data Bundle
                                  ,(/FIELD/)    &  !<-- ESMF Field holding the 2D integer array
                                  ,rc    =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          CALL gfs_physics_err_msg(RC,MESSAGE_CHECK,RC_PHY_OUT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          NUM_2D_FIELDS_I = NUM_2D_FIELDS_I + 1                            !<-- Add upt the number of 2D integer Fields
!
        ENDIF
!
        IF(PHY_INT_STATE_2D_I(file_index,NFIND)=='**********' ) THEN       !<-- End of the 2D integer array list
          EXIT
        ENDIF
!
      ENDDO
!      IF(MYPE==0)WRITE(0,*)' PHYSICS_OUTPUT: Number of 2-D Integer Fields=',NUM_2D_FIELDS_I
!
!-----------------------------------------------------------------------
!***  THE 2D REAL ARRAYS.
!-----------------------------------------------------------------------
!
      NUM_2D_FIELDS_R = 0
      LDIM1 = LBOUND(R_3D(1)%NAME,1)
      UDIM1 = UBOUND(R_3D(1)%NAME,1)
      LDIM2 = LBOUND(R_3D(1)%NAME,2)
      UDIM2 = UBOUND(R_3D(1)%NAME,2)
      NULLIFY(TEMP_R2D)
!
      if(file_id == 'OGFS_SFC') then
        PHY_INT_STATE_2D_R => PHY_INT_STATE_2D_R_SFC
      elseif(file_id == 'OGFS_FLX') then
        if ( int_state%climate) then
         PHY_INT_STATE_2D_R => PHY_INT_STATE_2D_R_FLX_CLM
        else
         PHY_INT_STATE_2D_R => PHY_INT_STATE_2D_R_FLX
        endif
      elseif(file_id == 'OGFS_AER') then                  
        PHY_INT_STATE_2D_R => PHY_INT_STATE_2D_R_AER
      elseif(file_id == 'OGFS_NST') then                  
        PHY_INT_STATE_2D_R => PHY_INT_STATE_2D_R_NST
      endif 
!
      DO NFIND=1,MAX_KOUNT
        IF(trim(PHY_INT_STATE_2D_R(file_index,NFIND)) == trim(file_id)) THEN    !<-- Take 2D real array data specified for history output
          VBL_NAME = TRIM(PHY_INT_STATE_2D_R(1,NFIND))//'_'//trim(PHY_INT_STATE_2D_R(3,NFIND))
          if(file_id == 'OGFS_SFC') then
            TEMP_R2D => R_3D(1)%NAME(LDIM1:UDIM1,LDIM2:UDIM2,NUM_2D_FIELDS_R+1)
          elseif(file_id == 'OGFS_FLX') then
            TEMP_R2D => R_3D(3)%NAME(LDIM1:UDIM1,LDIM2:UDIM2,NUM_2D_FIELDS_R+1)
          elseif(file_id == 'OGFS_AER') then
            TEMP_R2D => R_3D(4)%NAME(LDIM1:UDIM1,LDIM2:UDIM2,NUM_2D_FIELDS_R+1)
          elseif(file_id == 'OGFS_NST') then
            TEMP_R2D => R_3D(5)%NAME(LDIM1:UDIM1,LDIM2:UDIM2,NUM_2D_FIELDS_R+1)
          endif
!
!          write(0,*)'phy_out,R2D NFIND=',NFIND,'VBL_NAME=',trim(VBL_NAME),  &
!             'value=',maxval(R_3D(1)%NAME(LDIM1:UDIM1,LDIM2:UDIM2,NUM_2D_FIELDS_R+1)),    &
!              minval(R_3D(1)%NAME(LDIM1:UDIM1,LDIM2:UDIM2,NUM_2D_FIELDS_R+1))
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          MESSAGE_CHECK = "Insert Physics 2-D Real Data into Field"
!         CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          FIELD = ESMF_FieldCreate(grid         =GRID                   &  !<-- The ESMF grid
                                  ,farray       =TEMP_R2D               &  !<-- The 2D real array being inserted into the output data Bundle
                                  ,datacopyflag =COPYFLAG               &
                                  ,name         =VBL_NAME               &  !<-- Name of the 2D real array
                                  ,indexFlag=ESMF_INDEX_DELOCAL         &
                                  ,rc           =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          CALL gfs_physics_err_msg(RC,MESSAGE_CHECK,RC_PHY_OUT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          MESSAGE_CHECK = "Insert Physics 2-D Real Field into History Bundle"
!         CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          CALL ESMF_FieldBundleAdd(FILE_BUNDLE  &  !<-- The write component's output data Bundle
                                  ,(/FIELD/)    &  !<-- ESMF Field holding the 2D real array
                                  ,rc    =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          CALL gfs_physics_err_msg(RC,MESSAGE_CHECK,RC_PHY_OUT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          NUM_2D_FIELDS_R = NUM_2D_FIELDS_R + 1                             !<-- Add up the number of 2D real Fields
!
        ENDIF
!
        IF(PHY_INT_STATE_2D_R(file_index,NFIND) == '**********' ) THEN            !<-- End of the 2D real array list
          EXIT
        ENDIF
!
      ENDDO
!
!-----------------------------------------------------------------------
!***  THE 3D REAL ARRAYS.
!***  WE ARE WORKING WITH 3D ARRAYS BUT THEY ARE LOADED LAYER BY LAYER
!***  INTO 2D Fields.
!-----------------------------------------------------------------------
!
!
!################****!
!### all the 3D sfc variables are saved in 2D array. buff_bult_pieces
!###
      LDIM1 = LBOUND(R_3D(2)%NAME,1)
      UDIM1 = UBOUND(R_3D(2)%NAME,1)
      LDIM2 = LBOUND(R_3D(2)%NAME,2)
      UDIM2 = UBOUND(R_3D(2)%NAME,2)
      NLEVS = 0
      DO NFIND=1,MAX_KOUNT
!
        IF(trim(PHY_INT_STATE_3D_R(file_index,NFIND)) == trim(file_id)) THEN      !<-- Take 3D real array data specified for history output
!
         if(trim(PHY_INT_STATE_3D_R(4,NFIND)) == 'levs') then
           NDIM3 = int_state%levs
         elseif(trim(PHY_INT_STATE_3D_R(4,NFIND)) == 'levsp1') then
           NDIM3 = int_state%levs + 1
         elseif(trim(PHY_INT_STATE_3D_R(4,NFIND)) == 'lsoil') then
           NDIM3 = int_state%lsoil
         else
           write(0,*)'WRONG: The vertical level for 3D variable ',  &
             trim(PHY_INT_STATE_3D_R(file_index,NFIND)),'has not been specified!'
         endif
!
          DO K=1,NDIM3                                                     !<-- Loop through the levels of the array
            WRITE(MODEL_LEVEL,FMT)K
!
            VBL_NAME = TRIM(PHY_INT_STATE_3D_R(1,NFIND))//'_'//TRIM(PHY_INT_STATE_3D_R(3,NFIND)) &
     &         //'_'//MODEL_LEVEL //'_2D'
            NULLIFY(TEMP_R2D)
            TEMP_R2D => R_3D(2)%NAME(LDIM1:UDIM1,LDIM2:UDIM2,NLEVS+K)      !<-- Point at appropriate section of this IJK level 
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
            MESSAGE_CHECK = "Fill 2-D Fields with Each Level of Physics 3-D Data"
!           CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
            FIELD = ESMF_FieldCreate(grid         =GRID                 &  !<-- The ESMF grid
                                    ,farray       =TEMP_R2D             &  !<-- Level K of 3D real array being inserted into the data Bundle
                                    ,datacopyflag =COPYFLAG             &
                                    ,name         =VBL_NAME             &  !<-- Name of this level of the 3D real array
                                    ,indexFlag=ESMF_INDEX_DELOCAL         &
                                    ,rc           =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
            CALL gfs_physics_err_msg(RC,MESSAGE_CHECK,RC_PHY_OUT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
            MESSAGE_CHECK = "Insert Physics 3-D Data into History Bundle"
!           CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
             CALL ESMF_FieldBundleAdd(FILE_BUNDLE &  !<-- The write component's output data Bundle
                                    ,(/FIELD/)    &  !<-- ESMF Field holding the 1D real array
                                    ,rc    =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
            CALL gfs_physics_err_msg(RC,MESSAGE_CHECK,RC_PHY_OUT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
            NUM_2D_FIELDS_R=NUM_2D_FIELDS_R+1                              !<-- Continue adding up all levels of 3D Fields
!
          ENDDO
          NLEVS = NLEVS + NDIM3
!
        ENDIF
!
        IF(PHY_INT_STATE_3D_R(2,NFIND) == '**********' .AND.            &
           PHY_INT_STATE_3D_R(3,NFIND) == '**********'      ) THEN     !<-- End of the 3D real array list
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
      DO NFIND=1,MAX_KOUNT
!
        IF(trim(PHY_INT_STATE_4D_R(file_index,NFIND)) == 'file_id') THEN   !<-- Take 4D real array data specified for history output
!
         LDIM1 = LBOUND(R_4D(NFIND)%NAME,1)
         UDIM1 = UBOUND(R_4D(NFIND)%NAME,1)
         LDIM2 = LBOUND(R_4D(NFIND)%NAME,2)
         UDIM2 = UBOUND(R_4D(NFIND)%NAME,2)
         NDIM3 = 1
         NDIM4 = 1
!
         if(trim(PHY_INT_STATE_3D_R(4,NFIND)) == 'levs') then
           NDIM3 = int_state%levs
         elseif(trim(PHY_INT_STATE_3D_R(4,NFIND)) == 'levsp1') then
           NDIM3 = int_state%levs + 1
         endif
!
         if(trim(PHY_INT_STATE_3D_R(4,NFIND)) == 'ntrac') then
           NDIM4 = int_state%ntrac
         endif
!
         DO N=1,NDIM4                                                  !<-- Loop through the tracers kind
          DO K=1,NDIM3                                                 !<-- Loop through the levels of the array
            WRITE(TRACERS_KIND,FMT) N
            WRITE(MODEL_LEVEL, FMT) K
!
              VBL_NAME = TRIM(PHY_INT_STATE_4D_R(1,NFIND))//'_'//TRACERS_KIND//'_'//MODEL_LEVEL//'_2D'
              NULLIFY(TEMP_R2D)
              TEMP_R2D => R_4D(NFIND)%NAME(LDIM1:UDIM1,LDIM2:UDIM2,K,N)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
            MESSAGE_CHECK = "Fill 2-D Fields with Each Level of Physics 4-D Data"
!           CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
            FIELD=ESMF_FieldCreate(grid         =GRID                   &  !<-- The ESMF grid
                                  ,farray       =TEMP_R2D               &  !<-- Level K of 4D real array being inserted into the data Bundle
                                  ,datacopyflag =COPYFLAG               &
                                  ,name         =VBL_NAME               &  !<-- Name of this level of the 4D real array
                                  ,indexFlag=ESMF_INDEX_DELOCAL           &
                                  ,rc           =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
            CALL gfs_physics_err_msg(RC,MESSAGE_CHECK,RC_PHY_OUT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
            MESSAGE_CHECK = "Insert Physics 4-D Data into History Bundle"
!           CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
            CALL ESMF_FieldBundleAdd(FILE_BUNDLE  &  !<-- The write component's output data Bundle
                                    ,(/FIELD/)    &  !<-- ESMF Field holding the 1D real array
                                    ,rc    =RC)

!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
            CALL gfs_physics_err_msg(RC,MESSAGE_CHECK,RC_PHY_OUT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
            NUM_2D_FIELDS_R = NUM_2D_FIELDS_R + 1                          !<-- Continue adding up all levels of 4D Fields
!
          ENDDO
         ENDDO
!
        ENDIF
!
        IF(PHY_INT_STATE_4D_R(2,NFIND) == '**********' .AND.            &
           PHY_INT_STATE_4D_R(3,NFIND) == '**********'      ) THEN         !<-- End of the 4D real array list
          EXIT
        ENDIF
!
      ENDDO
!-----------------------------------------------------------------------
!
!WY bug fixed.
      MYPE = int_state%ME

      IF(MYPE == 0) WRITE(0,*)' PHYSICS_OUTPUT: Number of 2-D Real Fields=',NUM_2D_FIELDS_R
!
      END SUBROUTINE ADD_BUNDLE_TO_WRITE
!
!-----------------------------------------------------------------------
      END MODULE gfs_physics_output
!-----------------------------------------------------------------------
