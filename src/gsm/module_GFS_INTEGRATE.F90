!#include "../../ESMFVersionDefine.h"

!-----------------------------------------------------------------------
!
      MODULE MODULE_GFS_INTEGRATE
!
!-----------------------------------------------------------------------
!
!***  THIS MODULE HOLDS THE PRIMARY INTEGRATION RUNSTREAM OF THE GFS
!***  WITHIN SUBROUTINE GFS_INTEGRATE.
!
!-----------------------------------------------------------------------
!
! PROGRAM HISTORY LOG:
!   2009-12-23  Lu        - GFS_INTEGRATE modified to loop thru dyn, phy, &
!                           chem gridded component
!   2010-02-04  Lu        - GOCART_INTEGRATE added
!   2010-02-05  WANG      - change alarm set up for restart option of GFS
!   2010-03-09  Lu        - Add CHEM2PHY CPL
!   2010-08-18  WANG      - output filtered fields at HALFDFITIME
!   2010-11-10  WANG      - reset integration loop for dfi
!   2011-02     W Yang    - Updated to use both the ESMF 4.0.0rp2 library,
!                           ESMF 5 series library and the the
!                           ESMF 3.1.0rp2 library.
!   2011-03     W Yang    - Modified the digiter filter code for turning off
!                           the digiter filter case.
!   2011-10-01  Wang/Lu   - MYPE added to GOCART_INTEGRATE argument
!   2011-10     W Yang    - Modified for using the ESMF 5.2.0r library.
!   2015-08-17  S Moorthi - turn off history writing when restarting from
!                           history file
!   2015-10-22  S Moorthi - Modify digital filter to work for restart step
!                           e.g. second segment of forecast
!   2016-01-12  J Wang    - restore dyn-phys coupler for adiabatic option
!   2016-01-13  P Tripp   - NUOPC/GSM merge - changed NUOPC_ClockPrintCurrTime
!                           to ESMF_ClockPrint
!   2016-01-30  J Wang    - add high frequency output
!   2016-02-26  S Moorthi - fix for grid-point digital filter for semilag option
!   2016-03-05  S Moorthi - revise grid-point digital filter logic flow
!-----------------------------------------------------------------------

      USE ESMF
      use NUOPC
      USE MODULE_ERR_MSG
      USE MODULE_INCLUDE

      USE MODULE_DIGITAL_FILTER_GFS
      USE MODULE_GFS_WRITE,        ONLY: WRITE_ASYNC_GFS
      USE MODULE_GOCART_ROUTINES,  ONLY: GOCART_INTEGRATE
!
!-----------------------------------------------------------------------
!
      IMPLICIT NONE
!
!-----------------------------------------------------------------------
!
      PRIVATE
!
!
      PUBLIC :: GFS_INTEGRATE
!
!
      character(len=160) :: nuopcMsg

      CONTAINS
!
!-----------------------------------------------------------------------

      SUBROUTINE GFS_INTEGRATE(gc_gfs_dyn                            &
                              ,gc_gfs_phy                            &
                              ,GC_GFS_CHEM                           &
                              ,gc_gfs_cpl                            &
                              ,GC_PHY2CHEM_CPL                       &
                              ,GC_CHEM2PHY_CPL                       &
                              ,wrt_comps                             &
                              ,imp_gfs_dyn                           &
                              ,exp_gfs_dyn                           &
                              ,imp_gfs_phy                           &
                              ,exp_gfs_phy                           &
                              ,IMP_GFS_CHEM                          &
                              ,EXP_GFS_CHEM                          &
                              ,imp_gfs_wrt                           &
                              ,exp_gfs_wrt                           &
                              ,CLOCK_GFS                             &
                              ,OUTPUT_INTERVAL                       &
                              ,OUTPUT_INTERVAL_HF                    &
                              ,quilting                              &
                              ,WRITE_GROUP_READY_TO_GO               &
                              ,CURRTIME                              &
                              ,STARTTIME                             &
                              ,RESTRTTIME                            &
                              ,NTIMESTEP                             &
                              ,TIMESTEP                              &
                              ,NFHMAX_HF                             &
                              ,DFIHR                                 &
                              ,DFILEVS                               &
                              ,LDFI_GRD                              &
                              ,LDFIFLTO                              &
                              ,LWRTGRDCMP                            &
                              ,MYPE                                  &
                              ,PHYSICS_ON                            &
                              ,CHEMISTRY_ON, VM_local)
!
!-----------------------------------------------------------------------
!
      TYPE(ESMF_GridComp),INTENT(INOUT)      :: gc_gfs_dyn
      TYPE(ESMF_GridComp),INTENT(INOUT)      :: gc_gfs_phy &
                                               ,GC_GFS_CHEM        !<-- The Chemistry component
      TYPE(ESMF_CplComp), INTENT(INOUT)      :: gc_gfs_cpl
      TYPE(ESMF_CplComp), INTENT(INOUT)      :: GC_PHY2CHEM_CPL    !<-- The Phy-to-Chem coupler component
      TYPE(ESMF_CplComp), INTENT(INOUT)      :: GC_CHEM2PHY_CPL    !<-- The Chem-to-Phy coupler component

!jw
      TYPE(ESMF_GridComp),INTENT(INOUT)      :: wrt_comps(:)
      TYPE(ESMF_State),   INTENT(INOUT)      :: imp_gfs_dyn,exp_gfs_dyn
      TYPE(ESMF_State),   INTENT(INOUT)      :: imp_gfs_phy,exp_gfs_phy  &
                                               ,IMP_GFS_CHEM,EXP_GFS_CHEM  !<-- import/export states for Chemistry component

!jw
      TYPE(ESMF_State),   INTENT(INOUT)      :: imp_gfs_wrt,exp_gfs_wrt
      TYPE(ESMF_Clock),   INTENT(INOUT)      :: CLOCK_GFS                  !<-- GFS Component's ESMF Clock
      TYPE(ESMF_Time),    INTENT(INOUT)      :: CURRTIME               &   !<-- current forecast time
                                              , STARTTIME              &
                                              , RESTRTTIME
      INTEGER(KIND=KINT), INTENT(INOUT)      :: DFIHR, NTIMESTEP
      INTEGER(KIND=KINT), INTENT(IN)         :: DFILEVS                  !<-- level above which no DFI(WAM)
      INTEGER(KIND=KINT), INTENT(IN)         :: MYPE
      TYPE(ESMF_TimeInterval),INTENT(INout)  :: TIMESTEP                 !<-- ESMF timestep (s)
      TYPE(ESMF_TimeInterval)                :: donetime 
      TYPE(ESMF_Logical),INTENT(IN)          :: PHYSICS_ON               !<-- physics on (true) or off (false) switch
      TYPE(ESMF_Logical),INTENT(IN)          :: CHEMISTRY_ON             !<-- chemistry on (true) or off (false) switch
!jw
      TYPE(ESMF_TimeInterval),INTENT(INOUT)  :: OUTPUT_INTERVAL, OUTPUT_INTERVAL_HF
      LOGICAL,INTENT(IN)                     :: QUILTING, LDFIFLTO, LWRTGRDCMP, LDFI_GRD
      INTEGER(KIND=KINT),INTENT(IN)          :: NFHMAX_HF
      INTEGER(KIND=KINT),INTENT(INOUT)       :: WRITE_GROUP_READY_TO_GO
      TYPE(ESMF_VM), INTENT(IN)              :: VM_local
!
!---------------------
!***  Local variables
!---------------------
!
      INTEGER(KIND=KINT)                     :: RC,RC_LOOP,I
!     INTEGER(kind=ESMF_KIND_I8)             :: NTIMESTEP_ESMF         & !<-- current forecast timestep (ESMF_INT)
      INTEGER(kind=ESMF_KIND_I8)             :: NTIMESTEP_ESMF                   &
                                               ,NTIMESTEPH               !<-- timestep at fdfi

      INTEGER(KIND=KINT), save               :: NTIMESTEPB
      INTEGER(KIND=KINT)                     :: NDFISTEP

      TYPE(ESMF_Time)                        :: HALFDFITIME, DFITIME
      TYPE(ESMF_TimeInterval)                :: HALFDFIINTVAL
!jw
      TYPE(ESMF_Time)                        :: ALARM_OUTPUT_RING
      TYPE(ESMF_Time)                        :: ALARM_OUTPUT_HF_RING
      TYPE(ESMF_Time),  save                 :: ALARM_OUTPUT_HF_STOP
      TYPE(ESMF_Alarm), SAVE                 :: ALARM_OUTPUT
      TYPE(ESMF_Alarm), SAVE                 :: ALARM_OUTPUT_HF
      TYPE(ESMF_TimeInterval), save          :: OUTPUT_HFMAX
      LOGICAL                                :: Cpl_flag
      TYPE(ESMF_LOGICAL)                     :: Cpl_flag_ESMF
      LOGICAL, SAVE                          :: write_flag = .true.
      LOGICAL, SAVE                          :: first      = .true.
      LOGICAL, SAVE                          :: first_dfi  = .true.
      LOGICAL, SAVE                          :: end_dfi    = .false.
      LOGICAL, SAVE                          :: LDFI       = .false.
      LOGICAL                                :: LSPC, LALARM
      LOGICAL                                :: lalm1,lalm2,lskip
      INTEGER                                :: YY, MM, DD, H, M, S
!
      TYPE(ESMF_Field)                       :: FIELD
      TYPE(ESMF_info)                        :: info
!     real(ESMF_KIND_R8), dimension(:,:), pointer :: tmp_ptr2d
!-----------------------------------------------------------------------
!***  Set up alarm for output,alarm starts from current time
!-----------------------------------------------------------------------
!
      IF (first) THEN
        IF (NFHMAX_HF > 0) then
          CALL ESMF_TimeIntervalSet(OUTPUT_HFMAX                           &  !<-- ESMF time interval between GFS history output
                                   ,h           =NFHMAX_HF                 &  !<-- Hours between GFS history output
                                   ,m           =0                         &  !<-- Minutes between GFS history output
                                   ,s           =0                         &  !<-- Seconds between GFS history output
                                   ,rc          =RC)


          ALARM_OUTPUT_HF_STOP = starttime + OUTPUT_HFMAX + OUTPUT_INTERVAL_HF
          if (currtime <= starttime+output_hfmax) then
            ALARM_OUTPUT_HF_RING = CURRTIME  + OUTPUT_INTERVAL_HF
            ALARM_OUTPUT_RING = STARTTIME + OUTPUT_HFMAX + OUTPUT_INTERVAL

!
            ALARM_OUTPUT_HF = ESMF_AlarmCreate(name             ='ALARM_OUTPUT_HF'    &
                                              ,clock            =CLOCK_GFS            &  !<-- GFS Clock
                                              ,ringTime         =ALARM_OUTPUT_HF_RING &  !<-- Forecast/Restart start time (ESMF)
                                              ,ringInterval     =OUTPUT_INTERVAL_HF   &  !<-- Time interval between
                                              ,stoptime         =ALARM_OUTPUT_HF_STOP &  !<-- Time interval between
                                              ,ringTimeStepCount=1                    &  !<-- The Alarm rings for this many timesteps
                                              ,sticky           =.false.              &  !<-- Alarm does not ring until turned off
                                              ,rc               =RC)
          else
            ALARM_OUTPUT_RING = CURRTIME + OUTPUT_INTERVAL
          endif

        ELSE
          ALARM_OUTPUT_RING = CURRTIME + OUTPUT_INTERVAL
        ENDIF

!        CALL ESMF_TimeGet(time = ALARM_OUTPUT_RING                         &
!                         ,yy   = YY                                        &
!                         ,mm   = MM                                        &
!                         ,dd   = DD                                        &
!                         ,h    = H                                         &
!                         ,m    = M                                         &
!                         ,s    = S                                         &
!                         ,rc   = RC)
!
!         if (mype == 0) write(0,*)'set up alarm_output_ring,H=',H,'m=',m,'s=',s, &
!           'YYMMDD=',YY,MM,DD,'NFHMAX_HF=',NFHMAX_HF
!
         ALARM_OUTPUT = ESMF_AlarmCreate(name             ='ALARM_OUTPUT'     &
                                        ,clock            =CLOCK_GFS          &  !<-- GFS Clock
                                        ,ringTime         =ALARM_OUTPUT_RING  &  !<-- Forecast/Restart start time (ESMF)
                                        ,ringInterval     =OUTPUT_INTERVAL    &  !<-- Time interval between
                                        ,ringTimeStepCount=1                  &  !<-- The Alarm rings for this many timesteps
                                        ,sticky           =.false.            &  !<-- Alarm does not ring until turned off
                                        ,rc               =RC)

        call esmf_clockget(clock_gfs, timestep=timestep, starttime=starttime, &
                           currtime=currtime, rc=rc)
        donetime   = currtime - starttime
        NTIMESTEPB = nint(donetime/timeStep)

!       if (mype == 0) print *,' NTIMESTEPB =', NTIMESTEPB
        CALL ESMF_TimeGet(time = CURRTIME                         &
                         ,yy   = YY                               &
                         ,mm   = MM                               &
                         ,dd   = DD                               &
                         ,h    = H                                &
                         ,m    = M                                &
                         ,s    = S                                &
                         ,rc   = RC)
        if (mype == 0) write(0,*)' curr time=',yy,mm,dd,h,m,s,' first=',first
        first = .false.
      END IF
!
      IF (DFIHR > 0) THEN
        CALL ESMF_TimeIntervalSet(timeinterval=HALFDFIINTVAL, h=DFIHR, rc=RC)
        NDFISTEP    = HALFDFIINTVAL / TIMESTEP
        HALFDFITIME = RESTRTTIME    + HALFDFIINTVAL
        DFITIME     = HALFDFITIME   + HALFDFIINTVAL
        LDFI        = .true.
        NTIMESTEPH  = NDFISTEP
      ELSE
        NTIMESTEPH  = -1
        HALFDFITIME = RESTRTTIME
        DFITIME     = RESTRTTIME - timestep
      END IF
      if ((currtime-restrttime)/timestep > 2*ntimesteph) then
        ldfi = .false.
      endif
!
!-----------------------------------------------------------------------
!***  Execute the Run step of the Dynamics component
!-----------------------------------------------------------------------
!
      call ESMF_ClockPrint(CLOCK_GFS, options="currTime", &
      preString="Before entering GFS_Integrate time-loop - CLOCK current:                                   ", &
      unit=nuopcMsg)
      call ESMF_LogWrite(nuopcMsg, ESMF_LOGMSG_INFO)

      integrate: DO WHILE(.NOT.ESMF_ClockIsStopTime(CLOCK_GFS, rc = RC))
!                         .OR. (LDFI .and. DFIHR > 0) )

        CALL ESMF_LogWrite("Execute GFS Dynamics",ESMF_LOGMSG_INFO,rc=RC)
!
        call ESMF_ClockPrint(CLOCK_GFS, options="currTime", &
        preString="Right inside GFS_Integrate time-loop - CLOCK current:                                      ", &
        unit=nuopcMsg)
        call ESMF_LogWrite(nuopcMsg, ESMF_LOGMSG_INFO)

        CALL ESMF_GridCompRun(gridcomp    = GC_GFS_DYN                  &
                             ,importstate = IMP_GFS_DYN                 &
                             ,exportstate = EXP_GFS_DYN                 &
                             ,clock       = CLOCK_GFS                   &
                             ,rc          = RC)
!
        CALL ERR_MSG(RC,'execute dynamics',RC_LOOP)
!
        CALL ESMF_LogWrite("after dyn run, couple dyn_exp-to-phy_imp"   &
                          ,ESMF_LOGMSG_INFO,rc=rc)
        call ESMF_ClockPrint(CLOCK_GFS, options="currTime", &
        preString="Right after dyn.Run() - CLOCK current:                                                     ", &
          unit=nuopcMsg)
        call ESMF_LogWrite(nuopcMsg, ESMF_LOGMSG_INFO)
!
        lskip = .false.

        call esmf_clockget(clock_gfs, timestep=timestep, starttime=starttime, &
                           currtime=currtime, rc=rc)
        donetime  = currtime - starttime
        NTIMESTEP = nint(donetime/timeStep)


!       if (mype == 0) print *,' NTIMESTEP=',NTIMESTEP,' before disable alarm'  &
!                             ,' NTIMESTEPB=',NTIMESTEPB,'DFIHR=',DFIHR
!
!*** decide when to disable alarm

        IF(DFIHR > 0) THEN
!         LALM1 = .not.LDFIFLTO .and. CURRTIME > HALFDFITIME     &
!                 .and. CURRTIME <= DFITIME .and. first_dfi

          LALM2 = LDFIFLTO .and. CURRTIME >= HALFDFITIME .and. first_dfi

!         if (mype == 0) print *,'LALM2=',LALM2,'LDFIFLTO=',LDFIFLTO,  &
!           CURRTIME >= HALFDFITIME,'first_dfi=',first_dfi
!         IF(LALM1 .or. LALM2) THEN
          IF(LALM2) THEN
            if (nfhmax_hf > 0) then
              if (CURRTIME <= ALARM_OUTPUT_HF_STOP ) then
                CALL ESMF_AlarmDisable(alarm=ALARM_OUTPUT_HF, rc=RC)
              else
                call ESMF_AlarmDisable(ALARM_OUTPUT,rc=rc)
              endif
            else
              call ESMF_AlarmDisable(ALARM_OUTPUT,rc=rc)
            endif
            first_dfi = .false.
          ENDIF

!-----------------------------------------------------------------------
!***  Digital filter
!-----------------------------------------------------------------------
!
          filter_block: IF (LDFI) THEN
!
            IF(CURRTIME == RESTRTTIME) THEN

!***                                          Filter's first stage
!                                             --------------------
              if(LDFI_GRD) then
                CALL DIGITAL_FILTER_DYN_INIT_GFS(EXP_GFS_DYN,NDFISTEP,DFILEVS)
              endif

!***                                           The initial summation
!                                              ---------------------
              IF(PHYSICS_ON == ESMF_True) THEN
                CALL DIGITAL_FILTER_PHY_INIT_GFS(imp_gfs_phy)
              ENDIF

            else
!***                                            The summation stage
!                                               -------------------
              if(LDFI_GRD) then
!               if (mype == 0) write(0,*)' calling DIGITAL_FILTER_DYN_SUM_GFS', &
!                 ' at ntimestep=',ntimestep
                CALL DIGITAL_FILTER_DYN_SUM_GFS(EXP_GFS_DYN)
              endif
!
              IF(PHYSICS_ON == ESMF_True) THEN
                IF(CURRTIME == HALFDFITIME) THEN
                  CALL DIGITAL_FILTER_PHY_SAVE_GFS(IMP_GFS_PHY)
                ENDIF
              ENDIF
            ENDIF
!
!***                                            The final stage
!                                               ---------------
            CALL ESMF_TimeGet(time = CURRTIME                         &
                             ,yy   = YY                               &
                             ,mm   = MM                               &
                             ,dd   = DD                               &
                             ,h    = H                                &
                             ,m    = M                                &
                             ,s    = S                                &
                             ,rc   = RC)
            if (mype == 0) write(0,*)' curr time=',yy,mm,dd,h,m,s
!

            IF(CURRTIME == DFITIME) THEN
              if(LDFI_GRD) then
!               if (mype == 0) write(0,*)' calling DIGITAL_FILTER_DYN_AVERAGE_GFS',&
!                ' at ntimestep=',ntimestep
                CALL DIGITAL_FILTER_DYN_AVERAGE_GFS(exp_gfs_dyn)
              endif
!
              IF(PHYSICS_ON == ESMF_True) THEN
                CALL DIGITAL_FILTER_PHY_RESTORE_GFS(imp_gfs_phy)
              ENDIF
!
!             if (mype == 0) write(0,*)' setting clock back after digital filter'

              CALL ESMF_ClockSet(clock        = CLOCK_GFS             &
                                ,currtime     = HALFDFITIME           &
                                ,rc           = RC)

              call esmf_clockget(clock_gfs, timestep=timestep, starttime=starttime,&
                                 currtime=currtime, rc=rc)

              donetime = currtime - starttime
              NTIMESTEP = nint(donetime/timeStep)

              if (mype == 0) print *,' aft setting clock back NTIMESTEP=',NTIMESTEP

!
              call ESMF_ClockPrint(CLOCK_GFS, options="currTime", &
              preString="Right after ESMF_ClockSet() inside (CURRTIME == DFITIME) - CLOCK current:                  ", &
                      unit=nuopcMsg)
              call ESMF_LogWrite(nuopcMsg, ESMF_LOGMSG_INFO)

              DFITIME = RESTRTTIME
              DFIHR   = 0
              end_dfi = .true.
              lskip   = .true.
              ldfi    = .false.

            END IF
          ENDIF  filter_block
        END IF

!-----------------------------------------------------------------------
!*** enable alarm
!
        LALM1 = end_dfi .and. LDFIFLTO .and. currTime >= HALFDFITIME

!       if (mype == 0) print *,' LALM1=',LALM1,'ntimestep=',ntimestep, & 
!        'end_dfi=',end_dfi,' LDFIFLTO=',LDFIFLTO,CURRTIME <= ALARM_OUTPUT_HF_STOP
!
        IF (lalm1) then
          if (nfhmax_hf > 0) then
            if(CURRTIME <= ALARM_OUTPUT_HF_STOP ) then
              CALL ESMF_AlarmEnable(alarm = ALARM_OUTPUT_HF, rc=RC)
            else
              CALL ESMF_AlarmEnable(alarm = ALARM_OUTPUT, rc=RC)
            endif
          else 
            CALL ESMF_AlarmEnable(alarm = ALARM_OUTPUT, rc=RC)
          endif
          end_dfi = .false.
        endif
!
!-----------------------------------------------------------------------
!***  Call the Write component if it is time.
!-----------------------------------------------------------------------
!
        LSPC = (NTIMESTEP == 1 .and. currtime-timestep == starttime) .OR. &
               (LDFIFLTO .and. CurrTime == halfdfitime)
        write_flag = .false.
        LALARM = .false.
        if (nfhmax_hf > 0) then
          if(currtime <= starttime+output_hfmax) then
            if(ESMF_AlarmIsEnabled(alarm = ALARM_OUTPUT_HF, rc = RC)) then
              if( ESMF_AlarmIsRinging(alarm = ALARM_OUTPUT_HF,rc = Rc)) then
                LALARM = .true.
              endif
            endif
          else
            if(ESMF_AlarmIsEnabled(alarm = ALARM_OUTPUT, rc = RC)) then
              if(ESMF_AlarmIsRinging(alarm = ALARM_OUTPUT,rc = Rc)) then
                LALARM = .true.
              endif
            endif
          endif
        endif
        if(ESMF_AlarmIsEnabled(alarm = ALARM_OUTPUT, rc = RC)) then
          if(ESMF_AlarmIsRinging(alarm = ALARM_OUTPUT,rc = Rc)) then
            LALARM = .true.
          endif
        endif
!        if (mype == 0) print *,' LALARM=',LALARM,' LSPC=',LSPC,' LWRTGRDCMP=',LWRTGRDCMP, &
!           'lskip=',Lskip, 'halfdfitime?=',CurrTime == halfdfitime, &
!           'NTIMESTEP=',NTIMESTEP,currtime-timestep == starttime,'LDFIFLTO=',LDFIFLTO, &
!           'regular alarm=',ESMF_AlarmIsEnabled(alarm = ALARM_OUTPUT, rc = RC), &
!           ESMF_AlarmIsRinging(alarm = ALARM_OUTPUT,rc = Rc)

        outputdyn: IF((LALARM .or. LSPC) .AND. LWRTGRDCMP) THEN
          CALL WRITE_ASYNC_GFS(WRT_COMPs,exp_gfs_dyn               &
                              ,imp_gfs_wrt,exp_gfs_wrt             &
                              ,CLOCK_GFS                           &
                              ,MYPE                                &
                              ,WRITE_GROUP_READY_TO_GO)
        END IF outputdyn
!
!
        lskipif: if(.not. LSKIP) then                    ! if true, skip physics
!
!-----------------------------------------------------------------------
!***  import dynamics export data into physics-dynamics coupler and export to Physics
!-----------------------------------------------------------------------
!
          call esmf_cplcomprun(cplcomp     = gc_gfs_cpl          &
                              ,importstate = exp_gfs_dyn         &
                              ,exportstate = imp_gfs_phy         &
                              ,clock       = CLOCK_GFS           &
                              ,rc          = RC)
          call err_msg(RC,'couple dyn-to-phy',RC_LOOP)
!
          call ESMF_ClockPrint(CLOCK_GFS, options="currTime", &
               preString="Right after gc_gfs_cpl.Run() - CLOCK current:                                              ", &
               unit=nuopcMsg)
          call ESMF_LogWrite(nuopcMsg, ESMF_LOGMSG_INFO)
          IF (PHYSICS_ON == ESMF_True) THEN
!
!-----------------------------------------------------------------------
!***              Execute the Run step of the Physics Component
!-----------------------------------------------------------------------
!
            call esmf_logwrite("execute physics",ESMF_LOGMSG_INFO,rc=rc)
            call esmf_gridcomprun(gridcomp    = gc_gfs_phy            &
                                 ,importstate = imp_gfs_phy           &
                                 ,exportstate = exp_gfs_phy           &
                                 ,clock       = CLOCK_GFS             &
                                 ,rc          = RC)
            call err_msg(RC,'execute physics',RC_LOOP)
!check time step
            CALL ESMF_ClockGet(clock        = CLOCK_GFS               &
                              ,advanceCount = NTIMESTEP_ESMF          &  !<-- # of times the clock has advanced
                              ,rc           = RC)
!
            call ESMF_ClockPrint(CLOCK_GFS, options="currTime", &
                 preString="Right after gc_gfs_phy.Run() - CLOCK current:                                              ", &
                unit=nuopcMsg)
            call ESMF_LogWrite(nuopcMsg, ESMF_LOGMSG_INFO)
            NTIMESTEP = NTIMESTEP_ESMF
!
!-----------------------------------------------------------------------
!***              Invoke GOCART
!-----------------------------------------------------------------------
            IF (CHEMISTRY_ON == ESMF_True) THEN

              MESSAGE_CHECK = "Execute GOCART module"

              CALL GOCART_INTEGRATE(                                    &
                                   GC_GFS_CHEM,                         &
                                   GC_PHY2CHEM_CPL,                     &
                                   GC_CHEM2PHY_CPL,                     &
                                   EXP_GFS_PHY,                         &
                                   IMP_GFS_CHEM, EXP_GFS_CHEM,          &
                                   CLOCK_GFS, MYPE, RC )

              CALL ERR_MSG(RC,MESSAGE_CHECK,RC_LOOP)

              call ESMF_ClockPrint(CLOCK_GFS, options="currTime", &
                   preString="Right after GOCART_INTEGRATE() - CLOCK current:                                            ", &
                   unit=nuopcMsg)
              call ESMF_LogWrite(nuopcMsg, ESMF_LOGMSG_INFO)

            ENDIF
!
          ELSE       !***              Skip Physics if not turned on
!                                      -----------------------------

            call esmf_logwrite("pass phy_imp to phy_exp ",       &
                                ESMF_LOGMSG_INFO,rc=rc)
!
            call esmf_cplcomprun(              gc_gfs_cpl        &
                                ,importstate = imp_gfs_phy       &
                                ,exportstate = exp_gfs_phy       &
                                ,clock       = CLOCK_GFS         &
                                ,rc          = RC)
!
            call err_msg(RC,'pass phy_imp-to-phy_exp',RC_LOOP)

            call ESMF_ClockPrint(CLOCK_GFS, options="currTime",  &
                 preString="Right after gc_gfs_cpl.Run() WITHOUT PHYSICS - CLOCK current:                              ", &
                 unit=nuopcMsg)
            call ESMF_LogWrite(nuopcMsg, ESMF_LOGMSG_INFO)
          ENDIF

!
!-----------------------------------------------------------------------
!***  import physics export data into physics-dynamics coupler and export to Dynamics
!-----------------------------------------------------------------------
!
          call esmf_logwrite("couple phy_exp-to-dyn_imp",           &
                              ESMF_LOGMSG_INFO,rc=RC)
!
          call esmf_cplcomprun(cplcomp     = gc_gfs_cpl             &
                              ,importstate = exp_gfs_phy            &
                              ,exportstate = imp_gfs_dyn            &
                              ,clock       = CLOCK_GFS              &
                              ,rc          = RC)
!
          call err_msg(RC,'couple phy_exp-to-dyn_imp',RC_LOOP)
!
          call ESMF_ClockPrint(CLOCK_GFS, options="currTime",       &
               preString="Right after gc_gfs_cpl.Run() 2nd time - CLOCK current:                                     ", &
               unit=nuopcMsg)
          call ESMF_LogWrite(nuopcMsg, ESMF_LOGMSG_INFO)

!         Not sure if this cpl_flag stuff works - Moorthi
!         -----------------------------------------------
          CALL ESMF_InfoGetFromHost(imp_gfs_dyn, info, rc = rc)
          CALL ESMF_InfoGet(info, 'Cpl_flag', Cpl_flag, rc = rc)

          IF(.NOT. Cpl_flag) THEN

            CALL ESMF_ClockAdvance(clock = CLOCK_GFS, rc = RC)   ! advance ESMF clock
                                                                 ! ------------------

            call ESMF_ClockPrint(CLOCK_GFS, options="currTime", &
            preString="Right after ESMF_ClockAdvance() inside (.NOT. Cpl_flag) - CLOCK current: ", &
                    unit=nuopcMsg)
            call ESMF_LogWrite(nuopcMsg, ESMF_LOGMSG_INFO)
          END IF

          CALL ESMF_ClockGet(clock=CLOCK_GFS, currTime=CURRTIME, rc=RC)

        endif lskipif


      ENDDO integrate

      call ESMF_ClockPrint(CLOCK_GFS, options="currTime", &
      preString="Right after exiting GFS_Integrate time-loop - CLOCK current:                               ", &
              unit=nuopcMsg)
      call ESMF_LogWrite(nuopcMsg, ESMF_LOGMSG_INFO)


      call ESMF_ClockPrint(CLOCK_GFS, options="currTime", &
        preString="Right before returning from GFS_Integrate() - CLOCK current:                               ", &
        unit=nuopcMsg)
      call ESMF_LogWrite(nuopcMsg, ESMF_LOGMSG_INFO)

      END SUBROUTINE GFS_INTEGRATE
!
      END MODULE MODULE_GFS_INTEGRATE
