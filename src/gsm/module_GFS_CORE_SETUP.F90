!#include "../../ESMFVersionDefine.h"

!-----------------------------------------------------------------------
!
      MODULE MODULE_GFS_CORE_SETUP
!
!-----------------------------------------------------------------------
!
!***  THIS MODULE HOLDS THE DYNAMICS REGISTER, INIT, RUN, AND FINALIZE
!***  ROUTINES.  THEY ARE CALLED FROM THE ATM GRIDDED COMPONENT
!***  (ATM INITIALIZE CALLS DYNAMICS INITIALIZE, ETC.)
!***  IN MODULE_MAIN_GRID_COMP.F.
!
!
!-----------------------------------------------------------------------
!
      USE ESMF
      USE module_DM_PARALLEL_GFS  ,only: SETUP_SERVERS_GFS
      USE module_INCLUDE
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
      PUBLIC :: GFS_SETUP       !<-- An NMM-specific routine to set up parallelism and ESMF Grid
!
!
      CONTAINS

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
      SUBROUTINE GFS_SETUP(gc_atm,grid_atmos)
!
!-----------------------------------------------------------------------
!***  THIS ROUTINE CONTAINS NMM-SPECIFIC CODE FOR THE ATM COMPONENT:
!***    (1) SETTING UP DISTRIBUTED MEMORY PARALLELISM IN THE NMM;
!***    (2) CREATING THE ESMF Grid FOR THE ATM COMPONENT;
!***    (3) SHARING LOCAL SUBDOMAIN INDEX LIMITS AMONG TASKS.
!-----------------------------------------------------------------------
!
!
      USE MODULE_INCLUDE
      USE MODULE_GFS_MPI_DEF,ONLY : num_pes_fcst,last_fcst_pe          &
                                   ,first_fcst_pe                      &
                                   ,petlist_fcst                       &
                                   ,mpi_comm_inter                     &
                                   ,mc_comp,mpi_comm_comp              &
                                   ,quilting
      USE MODULE_IO_MPI_DEF,ONLY  : wrt_num_pes_fcst=>num_pes_fcst     &
                                   ,wrt_last_fcst_pe=>last_fcst_pe     &
                                   ,wrt_quilting=>quilting             &
                                   ,wrt_mpi_comm_comp=>mpi_comm_comp   &
                                   ,petlist_write,write_tasks_per_group&
                                   ,write_groups,max_inter_groups      &
                                   ,mpi_comm_inter_array 
!
      type(ESMF_gridcomp),intent(inout) :: gc_atm
      type(ESMF_grid),    intent(out)   :: grid_atmos  ! the ESMF grid for the integration attached to

!
!-----------------------------------------------------------------------
!***  SET UP THE ESMF GRID FOR THE NMM AND ESTABLISH THE
!***  DISTRIBUTED MEMORY PARALLELISM.
!-----------------------------------------------------------------------
!
      type(ESMF_config)            :: cf           ! the config object
!
      type(ESMF_DistGrid)          :: DistGrid_atmos
      type(ESMF_VM)                :: vm
!
      integer, dimension(2)        :: ncounts     ! parameter array to set up the
                                                  ! size of the 2-d ESMF grid.
      integer, dimension(2)        :: min,max     ! parameter arrays to set up the
                                                  ! start number and the end number of
                                                  ! the ESMF grid in each dimension.
      real(ESMF_kind_r8),dimension(ESMF_maxgriddim) :: mincoords,maxcoords
      integer,dimension(ESMF_maxgriddim) :: counts
      INTEGER , DIMENSION(2)             :: i1
      INTEGER , DIMENSION(:, :), POINTER :: i2
      INTEGER , DIMENSION(1)             :: mypelocal
      INTEGER , allocatable              :: petlistvm(:)
      integer                      :: num_pes_tot,num_pes,im,jm,lm
      integer                      :: mpi_intra,mpi_intra_b     ! the mpi intra-communicator
      integer                      :: rc,irtn,mype
      integer                      :: RC_RUN
      logical                      :: global
      character(50)                :: mode
      rc      = ESMF_success
      RC_RUN  = ESMF_success


!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!***  RETRIEVE THE VM (VIRTUAL MACHINE) OF THE ATM GRIDDED COMPONENT.
!***  CALL ESMF_GridCompGet TO RETRIEVE THE VM ANYWHERE YOU NEED IT.
!***  WE NEED VM NOW TO SET UP THE DE LAYOUT.
!-----------------------------------------------------------------------
!
!     CALL ESMF_logwrite("retrieve the config object and vm ",          &
!                        ESMF_logmsg_info,rc=RC)

!
      CALL ESMF_gridcompget(         gc_atm                             &
                           ,config  =cf                                 &
                           ,vm      =vm                                 &
                           ,rc      =RC)
!
!
!-----------------------------------------------------------------------
!***  set up parameters for mpi communications.
!***  use ESMF utility to get pe identification and total number of pes
!***  (referred to here as nodes.)
!-----------------------------------------------------------------------
!
!     CALL ESMF_logwrite("get mype and nodes from vm",ESMF_logmsg_info,rc=RC)
!
      CALL ESMF_vmget(vm                        &  !<-- the virtual machine
                     ,localpet=mype             &  !<-- local pe rank
                     ,petcount=num_pes          &  !<-- total # of tasks
                     ,rc      =RC)
!
      num_pes_tot  = num_pes
!
      mypelocal(1) = mype
      allocate(petlistvm(num_pes))
      CALL ESMF_VMAllGather(vm,                    &
                            sendData=mypelocal,    &
                            recvData=petlistvm,    &
                            count=1,               &
                            rc=RC)
!
!***  note: at this point, num_pes is the total number of mpi tasks,
!***        i.e., forecast tasks + quilt tasks.
!
!-----------------------------------------------------------------------
!***  establish the task layout including the quilt servers
!***  here in the main gridded component.  get the global
!***  mpi communicator and give it to setup_servers who will
!***  split it between forecast and quilt tasks.
!-----------------------------------------------------------------------
!
      CALL ESMF_vmget(vm                                  &
                     ,mpicommunicator=mpi_intra           &  !<-- the global communicator
                     ,rc             =RC)
!
      CALL mpi_comm_dup(mpi_intra,mpi_intra_b,rc)
!
      CALL ESMF_configgetattribute(cf                      &
                                  ,value =quilting         &  !<-- # of fcst tasks in j direction
                                  ,label ='quilting:'      &
                                  ,rc    =RC)
      wrt_quilting = quilting

!
!-----------------------------------------------------------------------
!***  set up quilt/write task specifications
!***  FIRST RETRIEVE THE TASK AND GROUP COUNTS FROM THE CONFIG FILE.
!-----------------------------------------------------------------------
      CALL ESMF_ConfigGetAttribute(CF                                   &  !<-- The configure file
                                  ,WRITE_GROUPS                         &  !<-- Number of write groups from config file
                                  ,label ='write_groups:'               &
                                  ,rc    =RC)
!
      CALL ESMF_ConfigGetAttribute(CF                                   &  !<-- The configure file
                                  ,WRITE_TASKS_PER_GROUP                &  !<-- Number of write tasks per group from config file
                                  ,label ='write_tasks_per_group:'      &
                                  ,rc    =RC)
!
      IF(quilting) THEN
          num_pes_fcst = num_pes_tot - WRITE_GROUPS * WRITE_TASKS_PER_GROUP
      ELSE
          num_pes_fcst = num_pes_tot
      END IF
      wrt_num_pes_fcst = num_pes_fcst
      
      allocate(petlist_fcst(num_pes_fcst))
      petlist_fcst(1:num_pes_fcst) = petlistvm(1:num_pes_fcst)
      last_fcst_pe     = maxval(petlist_fcst(1:num_pes_fcst) )
      wrt_last_fcst_pe = last_fcst_pe
      first_fcst_pe    = minval(petlist_fcst(1:num_pes_fcst) )

!      write(0,*)'gfs_setup,first_fcst_pe=',first_fcst_pe,'last_fcst_pe=', &
!        last_fcst_pe

      if(quilting) then

!-----------------------------------------------------------------------
!***  SEGREGATE THE FORECAST TASKS FROM THE QUILT/WRITE TASKS.
!-----------------------------------------------------------------------
!
        CALL SETUP_SERVERS_GFS(MYPE,NUM_PES,last_fcst_pe                  &
                              ,WRITE_GROUPS,WRITE_TASKS_PER_GROUP         &
                              ,mpi_intra_b,max_inter_groups               &
                              ,mpi_comm_inter_array)
        if (mype == 0)                                                    &
          write(0,*)'after setup_servers_gfs, write_groups=',write_groups,&
                    'WRITE_TASKS_PER_GROUP=', WRITE_TASKS_PER_GROUP,      &
                    'last_fcst_pe=', last_fcst_pe

      else            !if not quilt, for 1pe
        NUM_PES       = num_pes_fcst
        mc_comp       = mpi_intra_b
        mpi_comm_comp = mc_comp
      endif
      wrt_mpi_comm_comp = mpi_comm_comp
!***
!***  NOTE: At this point, NUM_PES is the number of Forecast tasks only.
!***
!
!-----------------------------------------------------------------------------
! Allocate the local index array i2 to store the local size information of the
! ditributed grid.  Information is based per dimension and per De.
!-----------------------------------------------------------------------------
      ALLOCATE(i2(2, num_pes_fcst))
!-----------------------------------------------------------------------
!***  create the ESMF grid.
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!***  the first dimension of ncounts is the i dimension for parallelization.
!***  the second dimension of ncounts is the j dimension.
!-----------------------------------------------------------------------
!
      CALL ESMF_configgetattribute(       cf                            &
                                  ,value =im                            &
                                  ,label ='im:'                         &
                                  ,rc    =RC)
!
      CALL ESMF_configgetattribute(       cf                            &
                                  ,value =jm                            &
                                  ,label ='jm:'                         &
                                  ,rc    =RC)
!
!---------------------------------------------------------
!***  if this is a global mode forecast, extend im and jm.
!***  retrieve the mode from the config file.
!---------------------------------------------------------
!
      CALL ESMF_configgetattribute(       cf                            &
                                  ,value =mode                          &
                                  ,label ='global:'                     &
                                  ,rc    =RC)
!
      if(trim(mode)=='.true.')then
        GLOBAL = .true.
        if (mype == 0) print *,'Run as a global model.'
        ncounts(1) = im + 2 !<-- global mode horizontal dimensions.
        ncounts(2) = jm + 2
      else
        GLOBAL = .false.
        if (mype == 0) print *,'Run as a regional model.'
        ncounts(1) = im     !<-- regional mode horizontal dimensions.
        ncounts(2) = jm
      endif
!
      max(1) = ncounts(1)
      max(2) = ncounts(2)
!
      min(1) = 1
      min(2) = 1
!
!-----------------------------------------------------------------------
!***  now create the main gridded component's ESMF grid.
!-----------------------------------------------------------------------
!
! Create the ESMF DistGrid_atmos.
!--------------------------------
!     CALL ESMF_LogWrite("Create DistGrid_atmos", ESMF_LOGMSG_INFO, rc = rc)

      DistGrid_atmos = ESMF_DistGridCreate(minIndex  = min,                 &
                                           maxIndex  = max,                 &
                                           regDecomp = (/1, num_pes_fcst/), &
                                           rc        = rc)


! Create the ESMF grid_atmos based on the created ESMF DistGrid_atmos information.
!---------------------------------------------------------------------------------
!     CALL ESMF_logwrite("create grid_atmos",ESMF_logmsg_info,rc=RC)
!
      grid_atmos = ESMF_GridCreate(name     = "grid_atmos",   &
                                   distgrid = DistGrid_atmos, &
                                   rc       = rc)
!
!-----------------------------------------------------------------------
!***  get the local array sizes for the main grid.
!***  again, only forecast tasks are relevant here.
!-----------------------------------------------------------------------
!
      if(mype<num_pes_fcst)then
          i2 = 0
          CALL ESMF_DistGridGet(DistGrid_atmos, indexCountPDe = i2, rc = rc)
      endif
!
!-----------------------------------------------------------------------
!***  USING 'computationalCount' FROM ARRAY I1 OBTAINED IN THE
!***  PREVIOUS CALL, GENERATE ALL OF THE LOCAL TASK INDEX LIMITS
!***  FOR ALL FORECAST TASKS.
!***  THE USER, NOT ESMF, DOES THIS WORK.
!------------------------------------------------- ----------------------
!
      END SUBROUTINE GFS_SETUP
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
      END MODULE MODULE_GFS_CORE_SETUP
!
!-----------------------------------------------------------------------

