!-----------------------------------------------------------------------
                        module module_dm_parallel
!-----------------------------------------------------------------------
!
!***  This module contains all codes directly related to distributed
!***  memory issues except for halo exchange although note that the
!***  halo widths must be set here.
!
!-----------------------------------------------------------------------
!
      use mpi
!
      use module_kinds
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!---domain decomposition info-------------------------------------------
!-----------------------------------------------------------------------
!
integer(kind=kint),parameter :: &
!
 ihalo=3 &                   ! halo width in I direction
,jhalo=3 &                   ! halo width in J direction
!
!***  Hardwire to 100 the maximum number of server groups allowed.
!***  This is far greater than should ever be needed.
!
,max_groups=100 &            ! max number of quilt server groups
!
!***  For now, set the number of threads here to 1.
!***  Clearly, this must be reconciled before actually using threading.
!
,num_tiles=1                 ! number of threads
!
!-----------------------------------------------------------------------
integer(kind=kint) :: &
 ide &                       ! ending data index, x direction
,ids &                       ! starting data index, x direction
,ims &                       ! the starting memory I for each task
,ime &                       ! the ending memory I for each task
,its &                       ! the starting integration I for each task
,ite &                       ! the ending integration I for each task
,jde &                       ! ending data index, y direction
,jds &                       ! starting data index, y direction
,jms &                       ! the starting memory J for each task
,jme &                       ! the ending memory J for each task
,jts &                       ! the starting integration J for each task
,jte &                       ! the ending integration J for each task
,lm  &                       ! the number of atmospheric model layers
,mpi_comm_comp &             ! local mpi communicator
,mpi_comm_inter &            ! intercommunicator for the quilt/write tasks
,mype_share &                ! my task ID to be seen by any USEs
,npes &                      ! total number of forecast tasks after SETUP_SERVERS 
!
,num_pts_max                 ! max points in any task's subdomain
!
integer(kind=kint) :: &
 its_b1 &                    ! its AND 1 point from global boundary
,its_b2 &                    ! its AND 2 points from global boundary
,its_h1 &                    ! its AND 1 point into halo
,its_h2 &                    ! its AND 2 points into halo
,its_b1_h1 &                 ! its AND _b1 AND _h1
,its_b1_h2 &                 ! its AND _b1 AND _h2
,ite_b1 &                    ! ite AND 1 point from global boundary
,ite_b2 &                    ! ite AND 2 points from global boundary
,ite_h1 &                    ! ite AND 1 point into halo
,ite_h2 &                    ! ite AND 2 points into halo
,ite_b1_h1 &                 ! ite AND _b1 AND _h1
,ite_b1_h2 &                 ! ite AND _b1 AND _h2
,jts_b1 &                    ! jts AND 1 point from global boundary
,jts_b2 &                    ! jts AND 2 points from global boundary
,jts_h1 &                    ! jts AND 1 point into halo
,jts_h2 &                    ! jts AND 2 points into halo
,jts_b1_h1 &                 ! jts AND _b1 AND _h1
,jts_b1_h2 &                 ! jts AND _b1 AND _h2
,jte_b1 &                    ! jte AND 1 point from global boundary
,jte_b2 &                    ! jte AND 2 points from global boundary
,jte_h1 &                    ! jte AND 1 point into halo
,jte_h2 &                    ! jte AND 2 points into halo
,jte_b1_h1 &                 ! jte AND _b1 AND _h1
,jte_b1_h2                   ! jte AND _b1 AND _h2
!
integer(kind=kint) :: &
 ide_m1 & 
,ide_m2 & 
,ids_p1 & 
,jde_m1 & 
,jde_m2 & 
,jds_p1 
integer(kind=kint),dimension(8) :: &
 my_neb                      ! my task's eight neighbors 
!
integer(kind=kint),dimension(max_groups) :: &
 mpi_intercomm_array &       ! intercommunicators for the integration tasks
,num_serv_per_grp            ! number of tasks in each group
!
integer(kind=kint),allocatable,dimension(:) :: &
 local_iend &
,local_istart &
,local_jend &
,local_jstart
!
integer :: mpi_intra
!-----------------------------------------------------------------------
!
       contains
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
      subroutine setup_servers(mype,inpes,jnpes,npes                    &
                              ,ngroups_write,write_tasks_per_group      &
                              ,mpi_intra,quilting)
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
!***  SETUP_SERVERS splits the communicator between integration 
!***  and output tasks.
!
!-----------------------------------------------------------------------
!
!   input argument list:
!    mype          - The fsct/quilt task rank in the domain intracommunicator.
!    inpes         - Number of mpi tasks in the X direction
!    jnpes         - Number of mpi tasks in the Y direction
!    npes          - Total number of mpi tasks provided to the job.  As input 
!                    to SETUP_SERVERS it includes the forecast tasks plus all
!                    write tasks in all groups of write tasks.  npes must at least
!                    equal the product of inpes*jnpes otherwise the integration
!                    cannot proceed. The difference between the product npes_fcst
!                    and npes is the number of mpi tasks that are available
!                    for i/o serving. This can be zero, in which case output will
!                    write a direct access file that can be separately quilted. 
!                    In order to skip the separate quilting step, make sure that
!                    the number of mpi tasks that the code is initiated with is at
!                    least one greater than npes_fcst.
!                    Later in the routine, npes is reset to the number of fcst tasks.
!    ngroups_write - Number of groups of write tasks.
!    mpi_intra     - The domain intracommunicator.
!    write_tasks_per_group - # of write tasks per write group
!
!-----------------------------------------------------------------------
!***  Argument variables.
!-----------------------------------------------------------------------
!
      integer(kind=kint),intent(in) :: &
 mype &                     ! each task's ID
,inpes &                    ! number of compute tasks in X direction
,jnpes &                    ! number of compute tasks in Y direction
,ngroups_write &            ! number of groups of write tasks
,write_tasks_per_group &    ! number of groups of write tasks per group
,mpi_intra                  ! global communicator
!
      logical(kind=klog),intent(in) :: &
quilting                    ! has output via quilting been specified?
!
      integer(kind=kint),intent(inout) :: &
 npes                       ! total number of tasks provided
                            ! then converted to the number of fcst tasks
!-----------------------------------------------------------------------
!***  Local variables.
!-----------------------------------------------------------------------
!
      integer(kind=kint) :: comdup,i,icc,icolor,iendq,iendxx,ierr &
                           ,igroup,igroup_x,iqserver,iquilt_group &
                           ,iss,issl &
                           ,istaq,istaxx,iworld,iworld_minus,ixx,jj,kk & 
                           ,lead_remote,npes_fcst,one
!
      integer(kind=kint),allocatable,dimension(:) :: irank
!
      logical :: include_mype
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!
!***  Let npes_fcst be the product of inpes and jnpes (namelist variables).
!***  This is the number of mpi tasks the executable has been built for.
!***  npes, returned from mpi_comm_size, must be at least this size
!***  otherwise the integration cannot proceed. The difference between
!***  npes_fcst and npes is the number of mpi tasks that are available
!***  for quilt/write serving. This can be zero, in which case output will
!***  write a direct access file that can be separately quilted.
!***  In order to skip the separate quilting step, make sure that
!***  the number of mpi tasks that the code is initiated with is at
!***  least one greater than npes_fcst.
!
!-----------------------------------------------------------------------
!
!!!   mype=mype_share
      npes_fcst=inpes*jnpes
      mpi_comm_comp=mpi_intra                                            !<-- Set mpi_comm_comp to the global communicator
!
      if(.not.quilting)then
        return                                                           !<-- If no output then nothing else to do.
      endif
!
!-----------------------------------------------------------------------
!
!***  At this point npes is the total number of mpi tasks. We will
!***  reset this at the end of the subroutine to the number of mpi
!***  tasks that are working on the model integration.
!
!***  First, however, we need to make sure that a sufficient number
!***  of mpi tasks have been initiated. If not, we will stop.
!
!-----------------------------------------------------------------------
!***  Compare the total number of MPI tasks (npes) to the number 
!***  specified for the forecast integration (npes_fcst=inpes*jnpes).
!***  Obviously the total number cannot be less than the number
!***  used for the forecast.
!-----------------------------------------------------------------------
!
      if(npes<npes_fcst)then
         write(0,*)' ***********************************************'
         write(0,*)' ***********************************************'
         write(0,*)' *************MAJOR PROBLEM*********************'
         write(0,*)' *************MAJOR PROBLEM*********************'
         write(0,*)' *************MAJOR PROBLEM*********************'
         write(0,*)' *************MAJOR PROBLEM*********************'
         write(0,*)' '
         write(0,*)' THERE ARE INSUFFICIENT MPI TASKS TO CONTINUE'
         write(0,*)' YOU MUST SPECIFY AT LEAST ',npes_fcst,' TASKS'
         write(0,*)' STOPPING NOW'
         write(0,*)' '
         write(0,*)' *************MAJOR PROBLEM*********************'
         write(0,*)' *************MAJOR PROBLEM*********************'
         write(0,*)' *************MAJOR PROBLEM*********************'
         write(0,*)' *************MAJOR PROBLEM*********************'
         write(0,*)' ***********************************************'
         write(0,*)' ***********************************************'
!!!      call mpi_abort(mpi_comm_world,1,ierr)
         call mpi_abort(mpi_intra,1,ierr)
      endif
!-----------------------------------------------------------------------
!
!***  OK, we have a sufficient number of mpi tasks to continue.
!
!***  How many groups of write tasks?  The default is 1 group.
!
!-----------------------------------------------------------------------
      one=1
!     iquilt_group=ngroups_read+ngroups_write
      iquilt_group=ngroups_write
      iquilt_group=max(iquilt_group,one)
!-----------------------------------------------------------------------
!
!***  ERROR CHECK FOR NUMBER OF GROUPS - MAXIMUM IS 100 - THAT IS A LOT!
!
!-----------------------------------------------------------------------
      if(iquilt_group>100)then
        write(0,*)' ***** IQUILT_GROUP IS GREATER THAN 100'
        write(0,*)' ***** DO YOU REALLY WANT THIS ?'
        write(0,*)' ***** IF SO THEN INCREASE SIZE IN mpp.h'
        write(0,*)' ***** ALSO, CHANGE IF CHECK IN SETUP_SERVERS'
        write(0,*)' ***** RESETTING THE NUMBER OF SERVER GROUPS TO 100'
        write(0,*)' ***** WE ARE CONTINUING ....   '
        iquilt_group=max_groups
      endif
!
      if(mype==0)then
        write(0,*)' Number of Server Groups: ',iquilt_group
      endif
!-----------------------------------------------------------------------
!
!***  Compute the number of quilt tasks per group.
!***  All mpi tasks beyond npes_fcst will be quilt tasks.
!***  If the number of tasks is not equally divisible by
!***  the number of groups of tasks then some groups may have
!***  more tasks then others.  This is fine.
!***  Note that we require at least one task per group.
!***  We may need to reduce the number of groups if
!***  it exceeds the number of tasks.
!
!-----------------------------------------------------------------------
!!!!  iqserver=npes-npes_fcst                                              !<-- Total # of quilt tasks in all groups
!     iqserver=ngroups_read*read_tasks_per_group                        &  !<-- Total # of quilt tasks in all groups
!             +ngroups_write*write_tasks_per_group
      iqserver=ngroups_write*write_tasks_per_group                         !<-- Total # of quilt tasks in all groups
!
      if(iqserver==0)then
        if(mype==0)then
          write(0,*)' *** You specified 0 Write tasks '
          write(0,*)' *** Output will write a direct access file'
        endif
        iquilt_group=0
      endif
!
      if(iquilt_group>iqserver)then
        iquilt_group=iqserver
        write(0,*)' ***** NOT ENOUGH WRITE/QUILT TASKS'
        write(0,*)' ***** WE NEED TO REDUCE THE NUMBER OF WRITE GROUPS'
        write(0,*)' ***** NUMBER OF GROUPS IS ',iquilt_group
      endif
!
      do i=0,iquilt_group-1
!!!!    call para_range(one,iqserver,iquilt_group,i,istaq,iendq)
!!!!    num_serv_per_grp(i+1)=iendq-istaq+1                              !<-- Store the # of tasks per group
        num_serv_per_grp(i+1)=write_tasks_per_group                      !<-- Store the # of tasks per group
        if(mype==0)then
          write(0,*)' Number of tasks for Group ',i+1,' is ',num_serv_per_grp(i+1)
        endif
      enddo
!
!-----------------------------------------------------------------------
!***  If there are more tasks executing this job that there are 
!***  forecast tasks plus all quilt tasks then warn the user!!!!
!-----------------------------------------------------------------------
!
      if(npes>npes_fcst+iqserver)then
!       write(0,*)' ABORTING IN SETUP_SERVERS'
        write(0,*)' MORE TASKS ARE EXECUTING THIS JOB THAN'             &
                 ,' THERE ARE FORECAST TASKS PLUS QUILT TASKS'
        write(0,*)' npes=',npes,' npes_fcst=',npes_fcst,' iqserver=',iqserver
!       call mpi_abort(mpi_intra,1,ierr)
      endif
!
!-----------------------------------------------------------------------
!***  Set up the "color" for mpi_comm_split.
!***  Those tasks which will do model integration will be color 0.
!***  The quilt tasks will have the color of the group number to
!***  which they will belong.
!-----------------------------------------------------------------------
!
      if(mype<npes_fcst)then
        icolor=0
      elseif(mype<npes)then
        istaxx=npes_fcst
        do i=1,iquilt_group
          iendxx=istaxx+num_serv_per_grp(i)-1
          if(mype>=istaxx.and.mype<=iendxx)then
            icolor=i
          endif
          istaxx=iendxx+1
        enddo
      endif
!
!-----------------------------------------------------------------------
!***  Split the communicator - The new intracommunicator for all tasks
!***  is mpi_comm_comp. mpi_comm_world (mpi_intra) is still available but it 
!***  refers to all the mpi tasks (model integration AND quilt tasks).
!-----------------------------------------------------------------------
!        
      call mpi_comm_dup(mpi_intra,comdup,ierr)
      call mpi_comm_split(comdup,icolor,mype,mpi_comm_comp,ierr)
!
!-----------------------------------------------------------------------
!***  At this point we have a new communicator, mpi_comm_comp,
!***  that can be used by the forecast tasks and the quilt tasks
!***  for their internal communications. On to the intercommunicators ...
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!***  Now we must create the intercommunicators for use between the mpi
!***  tasks doing the model integration and the mpi tasks for each 
!***  quilt group.  The first step is to exclude the tasks that do not
!***  belong.  We will do this for each quilt group by excluding the
!***  tasks from all of the other quilt groups.
!-----------------------------------------------------------------------
!
      allocate(irank(iqserver))                                          !<-- Dimension irank to the total # of quilt tasks
!
      if(iqserver>0)then
        do i=1,iqserver
          irank(i)=-1                                                    !<-- Initialize irank to meaningless values
        enddo
      endif
!
      ixx=npes_fcst                                                      !<-- Let ixx be the # of forecast tasks
!
!-----------------------------------------------------------------------
!
      inter_comm : do i=1,iquilt_group                                   !<-- Create intercommunicators between the set of fcst tasks
                                                                         !    and each quilt group individually.
!-----------------------------------------------------------------------
!
        include_mype=.true.
!
        if(mype<npes_fcst)then
          lead_remote=ixx                                                !<-- All fcst tasks set lead_remote to total # of fcst tasks.
                                                                         !    This is the rank of the lead quilt task in each
                                                                         !    quilt group as seen by the fcst tasks.  In other
                                                                         !    words it is the rank of the final fcst task plus 1.
        else
          lead_remote=0                                                  !<-- All quilt tasks set lead_remote to zero.
                                                                         !<-- This is the rank of the lead fcst task as seen by
                                                                         !<-- all quilt tasks.
        endif
!
        icc=0
        iss=npes_fcst
!
!***  This is the first possible task id that could be excluded.
!
        do jj=1,iquilt_group
          if(jj/=i)then
            issl=iss
            do kk=1,num_serv_per_grp(jj)
              icc=icc+1                               !<-- Add up number of excluded quilt tasks (those outside of current quilt group)
              irank(icc)=issl                         !<-- Save global IDs of the excluded quilt tasks
              if(mype==issl)include_mype=.false.      !<-- If MYPE is an excluded task, it sets its flag to .FALSE.
              issl=issl+1
            enddo
          endif
          iss=iss+num_serv_per_grp(jj)
        enddo
!
!-----------------------------------------------------------------------
!***  At this point we have an array, irank, with quilt task ids
!***  to exclude.  There are icc excluded tasks.
!***  Create a new group with all tasks in the other quilt groups
!***  excluded and then create a new intracommunicator (iworld_minus)
!***  that contains only the mpi tasks doing the model integration and
!***  the tasks that belong to the quilt group i we are considering.
!-----------------------------------------------------------------------
!
!!!     iworld=mpi_comm_world
        iworld=mpi_intra                                                 !<-- The global communicator
        call mpi_comm_group(iworld,igroup,ierr)                          !<-- igroup contains ALL tasks
        call mpi_group_excl(igroup,icc,irank,igroup_x,ierr)              !<-- igroup_x contains fcst tasks plus quilt tasks in group i
        call mpi_comm_create(iworld,igroup_x,iworld_minus,ierr)          !<-- iworld_minus is new intracommunicator for tasks in igroup_x
        call mpi_group_free(igroup,ierr)                                 !<-- Clear (set to null) igroup
        call mpi_group_free(igroup_x,ierr)                               !<-- Clear (set to null) igroup_x
!
!-----------------------------------------------------------------------
!***  At this point we have an intracommunicator that excludes the tasks we don't want.
!***  Create an intercommunicator for use between the mpi tasks doing the model
!***  integration and quilt group i that we are considering. This process is
!***  a collective routine that can only be done by the tasks that have not
!***  been excluded. Save this new intercommunicator in mpi_comm_inter for use
!***  by the tasks that belong to the quilt group that we are considering. The
!***  tasks that are performing the model integration will reference
!***  mpi_intercomm_array() since they will need to select which quilt
!***  group they wish to communicate with.
!-----------------------------------------------------------------------
!
        if(include_mype)then                                              !<-- Select all fcst tasks plus the quilt tasks in group i
          call mpi_intercomm_create(mpi_comm_comp                       & !<-- Each task's local intracommunicator
                                   ,0                                   & !<-- Rank of lead task in each local intracommunicator
                                   ,iworld_minus                        & !<-- Intracommunicator between fcst tasks and
                                                                          !    quilt tasks in group i
                                   ,lead_remote                         & !<-- Rank of leader in the remote group in iworld_minus
                                   ,0                                   & !<-- A tag
                                   ,mpi_intercomm_array(i)              & !<-- The new intercommunicator between the fcst tasks and
                                                                          !    the tasks in quilt group i
                                   ,ierr)
           mpi_comm_inter=mpi_intercomm_array(i)
        endif
!
        call mpi_barrier(mpi_intra,ierr)
!
!-----------------------------------------------------------------------
      enddo  inter_comm
!-----------------------------------------------------------------------
!***
!***  Set npes to the number of tasks working on the model integration.
!***
      npes=npes-iqserver
!
      if(mype==0)then
        write(0,*)' Number of Integration Tasks: ',npes
        write(0,*)' Total Number of Quilt Servers: ',iqserver
        write(0,*)' exit SETUP_SERVERS npes=',npes
      endif
!***
      deallocate (irank)
!-----------------------------------------------------------------------
!
      end subroutine setup_servers
!
!-----------------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!-----------------------------------------------------------------------
!
      subroutine para_range &
!
!!!   (n1,n2,nprocs,irank,ista,iend)
      (jnpes,ntasks_per_group,n_write_task,jrow_first,jrow_last)
!
!-----------------------------------------------------------------------
!***  Each write task will receive history data from a subset of the
!***  forecast tasks.  Entire rows (not partial rows) of forecast
!***  tasks will send to a write task.  Determine the range of rows
!***  of forecast tasks that will be sending data to each write
!***  task.  Row 1 is the row of forecast tasks along the domain's
!***  southern boundary.
!-----------------------------------------------------------------------
!
!   input argument list:
!     jnpes            - # of forecast tasks in the J direction
!     ntasks_per_group - # of write tasks per write group
!     n_write_task     - the "index" of the write task being
!                        considered given that every write group
!                        contains tasks 1-->ntasks_per_group
!
!   output argument list:
!     jrow_first - the first row of forecast tasks to send history
!                  data to this write task
!     jrow_last  - the last row of forecast tasks to send history
!                  data to this write task
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
      integer(kind=kint),intent(in) ::                                  &
        jnpes                                                           &
       ,n_write_task                                                    &
       ,ntasks_per_group
!
      integer(kind=kint),intent(out) ::                                 &
        jrow_first                                                      &
       ,jrow_last
!
!-----------------------------------------------------------------------
!***  Local variables
!-----------------------------------------------------------------------
!
      integer(kind=kint) ::                                             &
        iadd                                                            &
       ,n_remain                                                        &
       ,ntask                                                           &
       ,num_base   
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
      num_base=jnpes/ntasks_per_group
      n_remain=jnpes-ntasks_per_group*num_base
      jrow_last=0
!
      iadd=1
      do ntask=1,n_write_task
        if(ntask>n_remain)iadd=0
        jrow_first=jrow_last+1
        jrow_last=jrow_first+num_base+iadd-1
      enddo
!
!-----------------------------------------------------------------------
!
      end subroutine para_range
!
!-----------------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!-----------------------------------------------------------------------
!
      subroutine decomp &
!
      (mype,inpes,jnpes,npes_fcst,im,jm,lmx,global,ijcount)
!
!-----------------------------------------------------------------------
!
!***  DECOMP specifies the domain decomposition.
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!***  Argument variables.
!-----------------------------------------------------------------------
!
      integer(kind=kint),intent(in) :: im &
                                      ,jm &
                                      ,lmx &
                                      ,inpes &
                                      ,jnpes &
                                      ,mype &
                                      ,npes_fcst
!
      logical,intent(in) :: global
!
      integer(kind=kint),dimension(2),intent(in) :: ijcount
!
!-----------------------------------------------------------------------
!***  Local variables.
!-----------------------------------------------------------------------
!
      integer(kind=kint) :: i &
                           ,i_add &
                           ,icol &
                           ,iend &
                           ,ierr &
                           ,iguess &
                           ,ipe &
                           ,irecv &
                           ,iremain &
                           ,irtn &
                           ,isend &
                           ,istart &
                           ,istat &
                           ,j &
                           ,j_add  &
                           ,jend &
                           ,jguess &
                           ,jremain &
                           ,jrow &
                           ,jstart &
                           ,k2 &
                           ,l_remain &
                           ,lyr_frac &
                           ,my_e &
                           ,my_n &
                           ,my_ne &
                           ,my_nw &
                           ,my_s &
                           ,my_se &
                           ,my_sw &
                           ,my_w &
                           ,myi &
                           ,myj &
                           ,n &
                           ,npe &
                           ,num_pts 
!
      integer(kind=kint),dimension(4) :: limits
!
      integer(kind=kint),dimension(mpi_status_size) :: jstat
!
      integer(kind=kint),allocatable,dimension(:,:) :: ijcount_all      &
                                                      ,itemp
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!***  LMX is the number of model layers sent to this routine from the
!***  MAIN_GRID_COMP module where it was retrieved from the config file.
!***  The variable LM is available to all other modules from this one.
!***  Give LM the value of LMX so that any other module may use it.
!-----------------------------------------------------------------------
!
      lm=lmx
!
!-----------------------------------------------------------------------
!
      mype_share=mype  ! The gather/scatter routines and other parallel
                       ! routines will get mype from here and thus
                       ! remain ESMF-neutral themselves
!
!-----------------------------------------------------------------------
!***
!***  Compute the index limits within each MPI task.
!***  We will divide the number of points in the I and J directions
!***  evenly and then give as many of the initial tasks in each
!***  direction single additional points until any remainders are
!***  used up.
!***  The task IDs will start with 0 in the lower left corner and
!***  increase in the I direction and then in the J direction.
!***
!-----------------------------------------------------------------------
!***  The full dimensions of the integration domain.
!-----------------------------------------------------------------------
!
      if(global)then
        ids=1
        ide=im+2
        jds=1
        jde=jm+2
      else
        ids=1
        ide=im
        jds=1
        jde=jm
      endif
!
!-----------------------------------------------------------------------
!***  Find the remainders of points in each direction that will be
!***  incrementally added to each of the final tasks in each direction.
!-----------------------------------------------------------------------
!
      iguess=(ide-ids+1)/inpes
      iremain=(ide-ids+1)-iguess*inpes
      jguess=(jde-jds+1)/jnpes
      jremain=(jde-jds+1)-jguess*jnpes
!
!-----------------------------------------------------------------------
!***  Let every task know where all other tasks start and end
!***  on the full grid.
!***  Each task will save its own start/end values.
!-----------------------------------------------------------------------
!
      if(allocated(local_istart))then
        deallocate(local_istart)
        deallocate(local_iend)
        deallocate(local_jstart)
        deallocate(local_jend)
      endif
!
      npes=npes_fcst
      allocate(local_istart(0:npes-1),stat=istat)
      allocate(local_iend(0:npes-1),stat=istat)
      allocate(local_jstart(0:npes-1),stat=istat)
      allocate(local_jend(0:npes-1),stat=istat)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      if(mype==0)then
        if(allocated(ijcount_all))then
          deallocate(ijcount_all)
        endif
!
        allocate(ijcount_all(2,0:npes-1),stat=istat)
        ijcount_all(1,0)=ijcount(1)
        ijcount_all(2,0)=ijcount(2)
        local_istart(0)=ids
        local_iend(0)=ids+ijcount(1)-1
        local_jstart(0)=jds
        local_jend(0)=jds+ijcount(2)-1
!
        do npe=1,npes-1
          call mpi_recv(ijcount_all(1,npe),2,mpi_integer,npe,npe        &
                       ,mpi_comm_comp,jstat,irecv)
!
!-----------------------------------------------------------------------
!***  Find the start and end I values on this task's 
!***  primary integration region (inside the haloes).
!-----------------------------------------------------------------------
!
          icol=mod(npe,inpes)+1
          if(icol==1)then
            local_istart(npe)=ids
          else
            local_istart(npe)=local_istart(npe-1)+ijcount_all(1,npe-1)
          endif
          local_iend(npe)=local_istart(npe)+ijcount_all(1,npe)-1
!
!-----------------------------------------------------------------------
!***  Find the start and end J values on this task's 
!***  primary integration region (inside the haloes).
!-----------------------------------------------------------------------
!
          jrow=npe/inpes+1
          if(jrow==1)then
            local_jstart(npe)=jds
          else
            local_jstart(npe)=local_jstart(npe-inpes)+ijcount_all(2,npe-inpes)
          endif
          local_jend(npe)=local_jstart(npe)+ijcount_all(2,npe)-1
        enddo
!
      else
        call mpi_send(ijcount,2,mpi_integer,0,mype,mpi_comm_comp,isend)
      endif
!
      call mpi_bcast(local_istart,npes,mpi_integer,0,mpi_comm_comp,ierr)
      call mpi_bcast(local_iend,npes,mpi_integer,0,mpi_comm_comp,ierr)
      call mpi_bcast(local_jstart,npes,mpi_integer,0,mpi_comm_comp,ierr)
      call mpi_bcast(local_jend,npes,mpi_integer,0,mpi_comm_comp,ierr)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!-----------------------------------------------------------------------
!
      local_ij: do npe=0,npes-1
!
        if(mype==npe)then
          its=local_istart(npe)
          ite=local_iend(npe)
          jts=local_jstart(npe)
          jte=local_jend(npe)
!     write(0,*)' after bcast its=',its,' ite=',ite,' jts=',jts,' jte=',jte
        endif
!
      enddo local_ij
!
!-----------------------------------------------------------------------
!***  The memory or storage dimensions include the halo.
!-----------------------------------------------------------------------
!
!     ims=max(its-ihalo,ids)
!     ime=min(ite+ihalo,ide)
!     jms=max(jts-jhalo,jds)
!     jme=min(jte+jhalo,jde)
      ims=its-ihalo
      ime=ite+ihalo
      jms=jts-jhalo
      jme=jte+jhalo
!
!-----------------------------------------------------------------------
!***  Additional loop limits regarding global boundary and haloes.
!***  If "_bN" is appended to a start/end limit then that means
!***  that the MPI tasks along the global boundary must stay N
!***  points away from that boundary.
!***  If "_hN" is appended to a start/end limit then that means
!***  that each task will compute N points into the halo unless
!***  being blocked by the global boundary.
!-----------------------------------------------------------------------
!
      ids_p1=max(its,ids+1)
      ide_m1=min(ite,ide-1)
      ide_m2=min(ite,ide-2)
      jds_p1=max(jts,jds+1)
      jde_m1=min(jte,jde-1)
      jde_m2=min(jte,jde-2)
!
      its_b1=max(its,ids+1)
      ite_b1=min(ite,ide-1)
      its_b2=max(its,ids+2)
      ite_b2=min(ite,ide-2)
      jts_b1=max(jts,jds+1)
      jte_b1=min(jte,jde-1)
      jts_b2=max(jts,jds+2)
      jte_b2=min(jte,jde-2)
!
      its_h1=max(its-1,ids)
      ite_h1=min(ite+1,ide)
      its_h2=max(its-2,ids)
      ite_h2=min(ite+2,ide)
      jts_h1=max(jts-1,jds)
      jte_h1=min(jte+1,jde)
      jts_h2=max(jts-2,jds)
      jte_h2=min(jte+2,jde)
!
      its_b1_h1=max(its-1,ids+1)
      ite_b1_h1=min(ite+1,ide-1)
      ite_b1_h2=min(ite+2,ide-1)
      jts_b1_h1=max(jts-1,jds+1)
      jte_b1_h1=min(jte+1,jde-1)
      jte_b1_h2=min(jte+2,jde-1)
!
      if(mype==0)then
        write(0,*)' ids=',ids,' ide=',ide,' jds=',jds,' jde=',jde
      endif
!
      do npe=0,npes-1
        if(mype==npe)then
!!!       write(0,*)' PE=',mype
          write(0,*)' inpes=',inpes,' jnpes=',jnpes
          write(0,*)' its=',its,' ite=',ite,' jts=',jts,' jte=',jte
          write(0,*)' ims=',ims,' ime=',ime,' jms=',jms,' jme=',jme
          write(0,*)' ids=',ids,' ide=',ide,' jds=',jds,' jde=',jde
        endif
        call mpi_barrier(mpi_comm_comp,irtn)
      enddo
!-----------------------------------------------------------------------
!***  Find the maximum horizontal size of each task's subdomain
!***  since task 0 will need that in subroutine DSTRB.
!-----------------------------------------------------------------------
!
      num_pts_max=0
!
      if(mype==0)then
        do ipe=1,npes-1
          call mpi_recv(limits,4,mpi_integer,ipe,ipe,mpi_comm_comp      &
     &,                 jstat,irecv)
!
          istart=limits(1)
          iend=limits(2)
          jstart=limits(3)
          jend=limits(4)
!
          num_pts=(iend-istart+1)*(jend-jstart+1)
          if(num_pts>num_pts_max)then
            num_pts_max=num_pts
          endif
        enddo
!
      else
!
        limits(1)=its
        limits(2)=ite
        limits(3)=jts
        limits(4)=jte
!
        call mpi_send(limits,4,mpi_integer,0,mype,mpi_comm_comp,isend)
      endif
!
!-----------------------------------------------------------------------
!***  Let each task determine who its eight neighbors are because we
!***  will need to know that for the halo exchanges.  The direction
!***  to each neighbor will be designated by the following integers:
!     
!***      north: 1
!***       east: 2
!***      south: 3
!***       west: 4
!***  northeast: 5
!***  southeast: 6
!***  southwest: 7
!***  northwest: 8
!
!***  If a task has no neighbor in a particular direction because of
!***  the presence of the global domain boundary then that element
!***  of my_neb is set to -1.
!-----------------------------------------------------------------------
!
      if(allocated(itemp))then
        deallocate(itemp)
      endif
!
      allocate(itemp(inpes,jnpes),stat=istat)
      ipe=0
!
      do j=1,jnpes
      do i=1,inpes
        itemp(i,j)=ipe
        if(ipe==mype)then
          myi=i
          myj=j
        endif
        ipe=ipe+1
      enddo
      enddo
!
      my_n=-1
      if(myj+1<=jnpes)my_n=itemp(myi,myj+1)
!
      my_e=-1
      if(myi+1<=inpes)my_e=itemp(myi+1,myj)
!
      my_s=-1
      if(myj-1>=1)my_s=itemp(myi,myj-1)
!
      my_w=-1
      if(myi-1>=1)my_w=itemp(myi-1,myj)
!
      my_ne=-1
      if((myi+1<=inpes).and.(myj+1<=jnpes)) &
         my_ne=itemp(myi+1,myj+1)
!
      my_se=-1
      if((myi+1<=inpes).and.(myj-1>=1)) &
         my_se=itemp(myi+1,myj-1)
!
      my_sw=-1
      if((myi-1>=1).and.(myj-1>=1)) &
         my_sw=itemp(myi-1,myj-1)
!
      my_nw=-1
      if((myi-1>=1).and.(myj+1<=jnpes)) &
         my_nw=itemp(myi-1,myj+1)
!
      my_neb(1)=my_n
      my_neb(2)=my_e
      my_neb(3)=my_s
      my_neb(4)=my_w
      my_neb(5)=my_ne
      my_neb(6)=my_se
      my_neb(7)=my_sw
      my_neb(8)=my_nw
!
      deallocate(itemp)
!
!!!   write(0,*)' Exiting DECOMP'
!!!   write(0,*)' its=',its,' ite=',ite,' jts=',jts,' jte=',jte
!!!   write(0,*)' ims=',ims,' ime=',ime,' jms=',jms,' jme=',jme
!-----------------------------------------------------------------------
!
      end subroutine decomp
!
!-----------------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!-----------------------------------------------------------------------
                        subroutine dstrb &
      (arrayg,arrayl,lgs,lge,lls,lle,l1,mype,mpi_comm_comp)
!-----------------------------------------------------------------------
!     DSTRB distributes the elements of real global array arrayg to the
!     real local array arrayl. 
!-----------------------------------------------------------------------
!     input argument list:
!       arrayg - global array
!       lgs    - starting vertical index of global array
!       lge    - ending vertical index of global array
!       lls    - starting vertical index of local array
!       lle    - ending vertical index of local array
!       l1     - vertical level of arrayl being filled in this call
!                (used only when lge=1 and lle>1, i.e. when the global
!                 array is actually just one level of a multi_level
!                 array)
!       mype   - task rank
!       mpi_comm_comp - the local intracommunicator
!
!     output argument list:
!       arrayl - local array
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!***
!***  argument variables
!***
      integer(kind=kint),intent(in) :: l1,lge,lgs,lle,lls
!
      integer(kind=kint),intent(in) :: mype,mpi_comm_comp
!
      real(kind=kfpt),dimension(ids:ide,jds:jde,lgs:lge),intent(in) :: &
                                                                  arrayg
      real(kind=kfpt),dimension(ims:ime,jms:jme,lls:lle),intent(out) :: &
                                                                  arrayl
!
!-----------------------------------------------------------------------
!***
!***  local variables
!***
!
      integer(kind=kint) :: i,iend,ipe,irecv,irtn,isend,istart,j,jend &
                           ,jstart,knt,l,numvals
      integer,dimension(4) :: limits
      integer,dimension(mpi_status_size) :: jstat
!
      real(kind=kfpt),allocatable,dimension(:) :: arrayx
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!***  Task 0 fills its own local domain then parcels out all the other 
!***  pieces to the other tasks.
!-----------------------------------------------------------------------
!
      tasks : if(mype==0)then
!
        if(lge==lgs)then
          do j=jts,jte
          do i=its,ite
            arrayl(i,j,l1)=arrayg(i,j,lgs)
          enddo
          enddo
!
        else
!
          do l=lgs,lge
            do j=jts,jte
            do i=its,ite
              arrayl(i,j,l)=arrayg(i,j,l)
            enddo
            enddo
          enddo
        endif
!
!***  Task 0 creates an array to hold all the values from the other
!***  tasks' pieces of the global array and then sends those pieces 
!***  out to the appropriate tasks.
!
        numvals=num_pts_max*(lge-lgs+1)
        allocate(arrayx(numvals),stat=i)
!
        do ipe=1,npes-1
!
          call mpi_recv(limits,4,mpi_integer,ipe,ipe,mpi_comm_comp &
                       ,jstat,irecv)
!
          istart=limits(1)
          iend=limits(2)
          jstart=limits(3)
          jend=limits(4)
          knt=0
!
          do l=lgs,lge
            do j=jstart,jend
            do i=istart,iend
              knt=knt+1
              arrayx(knt)=arrayg(i,j,l)
            enddo
            enddo
          enddo
!
          call mpi_send(arrayx,knt,mpi_real,ipe,ipe,mpi_comm_comp,isend)
!
        enddo
!
        deallocate(arrayx)
!
!-----------------------------------------------------------------------
!***  All other tasks tell task 0 what their horizontal limits are and
!***  receive their piece of the global array from task 0.
!-----------------------------------------------------------------------
!
      else
!
        limits(1)=its
        limits(2)=ite
        limits(3)=jts
        limits(4)=jte
!
        call mpi_send(limits,4,mpi_integer,0,mype,mpi_comm_comp,isend)
!
        knt=(ite-its+1)*(jte-jts+1)*(lge-lgs+1)
        allocate(arrayx(knt),stat=i)
!
        call mpi_recv(arrayx,knt,mpi_real,0,mype,mpi_comm_comp &
                     ,jstat,irecv)
!
        knt=0
        if(lge==lgs)then
          do j=jts,jte
          do i=its,ite
            knt=knt+1
            arrayl(i,j,l1)=arrayx(knt)
          enddo
          enddo
        else
          do l=lgs,lge
            do j=jts,jte
            do i=its,ite
              knt=knt+1
              arrayl(i,j,l)=arrayx(knt)
            enddo
            enddo
          enddo
        endif
!
        deallocate(arrayx)
!
!-----------------------------------------------------------------------
!
      endif tasks
!
!-----------------------------------------------------------------------
      call mpi_barrier(mpi_comm_comp,irtn)
!-----------------------------------------------------------------------
!
      end subroutine dstrb
!
!-----------------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!-----------------------------------------------------------------------
                        subroutine idstrb &
      (iarrayg,iarrayl,mype,mpi_comm_comp)
!-----------------------------------------------------------------------
!     IDSTRB distributes the elements of integer global array iarrayg
!     to the integer local array iarrayl. 
!-----------------------------------------------------------------------
!     input argument list:
!       iarrayg - global array
!
!     output argument list:
!       iarrayl - local array
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!***
!***  argument variables
!***
      integer(kind=kint),intent(in) :: mype,mpi_comm_comp
!
      integer(kind=kint),dimension(ids:ide,jds:jde),intent(in) :: &
                                                                 iarrayg
      integer(kind=kint),dimension(ims:ime,jms:jme),intent(out) :: &
                                                                 iarrayl
!-----------------------------------------------------------------------
!***
!***  local variables
!***
!
      integer(kind=kint) :: i,iend,ipe,irecv,irtn,isend,istart,j,jend &
!xxx                       ,jstart,knt,l,mype,numvals
                           ,jstart,knt,l,numvals
      integer,dimension(4) :: limits
      integer,dimension(mpi_status_size) :: jstat
!
      integer(kind=kint),allocatable,dimension(:) :: iarrayx
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
!***  Initialize the output array.
!
      do j=jms,jme
      do i=ims,ime
        iarrayl(i,j)=0.
      enddo
      enddo
!
!-----------------------------------------------------------------------
!***  Task 0 fills its own local domain then parcels out all the other 
!***  pieces to the other tasks.
!-----------------------------------------------------------------------
!
      tasks : if(mype==0)then
!
        do j=jts,jte
        do i=its,ite
          iarrayl(i,j)=iarrayg(i,j)
        enddo
        enddo
!
!***  Task 0 creates an array to hold all the values from the other
!***  tasks' pieces of the global array and then sends those pieces 
!***  out to the appropriate tasks.
!
        numvals=num_pts_max
        allocate(iarrayx(numvals),stat=i)
!
        do ipe=1,npes-1
!
          call mpi_recv(limits,4,mpi_integer,ipe,ipe,mpi_comm_comp &
                       ,jstat,irecv)
!
          istart=limits(1)
          iend=limits(2)
          jstart=limits(3)
          jend=limits(4)
          knt=0
!
          do j=jstart,jend
          do i=istart,iend
            knt=knt+1
            iarrayx(knt)=iarrayg(i,j)
          enddo
          enddo
!
          call mpi_send(iarrayx,knt,mpi_integer,ipe,ipe,mpi_comm_comp &
                       ,isend)
!
        enddo
!
        deallocate(iarrayx)
!
!-----------------------------------------------------------------------
!***  All other tasks tell task 0 what their horizontal limits are and
!***  receive their piece of the global array from task 0.
!-----------------------------------------------------------------------
!
      else
!
        limits(1)=its
        limits(2)=ite
        limits(3)=jts
        limits(4)=jte
!
        call mpi_send(limits,4,mpi_integer,0,mype,mpi_comm_comp,isend)
!
        knt=(ite-its+1)*(jte-jts+1)
        allocate(iarrayx(knt),stat=i)
!
        call mpi_recv(iarrayx,knt,mpi_integer,0,mype,mpi_comm_comp &
                     ,jstat,irecv)
!
        knt=0
        do j=jts,jte
        do i=its,ite
          knt=knt+1
          iarrayl(i,j)=iarrayx(knt)
        enddo
        enddo
!
        deallocate(iarrayx)
!
!-----------------------------------------------------------------------
!
      endif tasks
!
!-----------------------------------------------------------------------
      call mpi_barrier(mpi_comm_comp,irtn)
!-----------------------------------------------------------------------
!
      end subroutine idstrb
!
!-----------------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!-----------------------------------------------------------------------
                        subroutine dstrb_soil &
      (arrayg,arrayl,lgs,lge,lls,lle)
!-----------------------------------------------------------------------
!     DSTRB distributes the elements of real global array arrayg to the
!     real local array arrayl. 
!-----------------------------------------------------------------------
!     input argument list:
!       arrayg - global soil array
!       lgs    - starting vertical index of global array
!       lge    - ending vertical index of global array
!       lls    - starting vertical index of local array
!       lle    - ending vertical index of local array
!
!     output argument list:
!       arrayl - local soil array
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!***
!***  argument variables
!***
      integer(kind=kint),intent(in) :: lge,lgs,lle,lls
!
      real(kind=kfpt),dimension(lgs:lge,ids:ide,jds:jde),intent(in) :: &
                                                                  arrayg
      real(kind=kfpt),dimension(lls:lle,ims:ime,jms:jme),intent(out) :: &
                                                                  arrayl
!-----------------------------------------------------------------------
!***
!***  local variables
!***
!
      integer(kind=kint) :: i,iend,ipe,irecv,irtn,isend,istart,j,jend &
                           ,jstart,knt,l,mype,numvals
      integer,dimension(4) :: limits
      integer,dimension(mpi_status_size) :: jstat
!
      real(kind=kfpt),allocatable,dimension(:) :: arrayx
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
      mype=mype_share
!
!***  Initialize the output array.
!
      do j=jms,jme
      do i=ims,ime
      do l=lls,lle
        arrayl(l,i,j)=0.
      enddo
      enddo
      enddo
!
!-----------------------------------------------------------------------
!***  Task 0 fills its own local domain then parcels out all the other 
!***  pieces to the other tasks.
!-----------------------------------------------------------------------
!
      tasks : if(mype==0)then
!
        do j=jts,jte
        do i=its,ite
          do l=lgs,lge
            arrayl(l,i,j)=arrayg(l,i,j)
          enddo
        enddo
        enddo
!
!***  Task 0 creates an array to hold all the values from the other
!***  tasks' pieces of the global array and then sends those pieces 
!***  out to the appropriate tasks.
!
        numvals=num_pts_max*(lge-lgs+1)
        allocate(arrayx(numvals),stat=i)
!
        do ipe=1,npes-1
!
          call mpi_recv(limits,4,mpi_integer,ipe,ipe,mpi_comm_comp &
                       ,jstat,irecv)
!
          istart=limits(1)
          iend=limits(2)
          jstart=limits(3)
          jend=limits(4)
          knt=0
!
          do j=jstart,jend
          do i=istart,iend
            do l=lgs,lge
              knt=knt+1
              arrayx(knt)=arrayg(l,i,j)
            enddo
          enddo
          enddo
!
          call mpi_send(arrayx,knt,mpi_real,ipe,ipe,mpi_comm_comp,isend)
!
        enddo
!
        deallocate(arrayx)
!
!-----------------------------------------------------------------------
!***  All other tasks tell task 0 what their horizontal limits are and
!***  receive their piece of the global array from task 0.
!-----------------------------------------------------------------------
!
      else
!
        limits(1)=its
        limits(2)=ite
        limits(3)=jts
        limits(4)=jte
!
        call mpi_send(limits,4,mpi_integer,0,mype,mpi_comm_comp,isend)
!
        knt=(ite-its+1)*(jte-jts+1)*(lge-lgs+1)
        allocate(arrayx(knt),stat=i)
!
        call mpi_recv(arrayx,knt,mpi_real,0,mype,mpi_comm_comp &
                     ,jstat,irecv)
!
        knt=0
        do j=jts,jte
        do i=its,ite
          do l=lgs,lge
            knt=knt+1
            arrayl(l,i,j)=arrayx(knt)
          enddo
        enddo
        enddo
!
        deallocate(arrayx)
!
!-----------------------------------------------------------------------
!
      endif tasks
!
!-----------------------------------------------------------------------
      call mpi_barrier(mpi_comm_comp,irtn)
!-----------------------------------------------------------------------
!
      end subroutine dstrb_soil
!
!-----------------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!-----------------------------------------------------------------------
                        subroutine gather_layers &
      (field,lm,npes,msize_dummy_fft &
      ,lm_fft,k1_fft,k2_fft &
      ,local_istart,local_iend &
      ,local_jstart,local_jend &
      ,jstart_fft,jend_fft &
      ,my_jrow_start,my_jrow_end &
      ,ipe_start,ipe_end &
      ,my_domain_has_fft_lats &
      ,mype,mpi_comm_comp &
      ,its,ite,ims,ime,ids,ide,jts,jte,jms,jme &
      ,array_lyrs)
!-----------------------------------------------------------------------
!***  GATHER_LAYERS distributes all the elements of field between layers
!***  k1 and k2 inclusive to the appropriate task for subsequent 
!***  application of FFTs.
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
integer(kind=kint),intent(in) :: &
 its &
,ite &
,ims &
,ime &
,ids &
,ide &
,jts &
,jte &
,jms &
,jme
!
integer(kind=kint),intent(in) :: &
 ipe_end &
,ipe_start &
,jstart_fft &
,jend_fft &
,lm &
,lm_fft &
,mpi_comm_comp &
,msize_dummy_fft &
,mype &
,npes 
!
integer(kind=kint),dimension(0:npes-1),intent(in) :: &
 k1_fft &
,k2_fft &
,local_iend &
,local_istart &
,local_jend &
,local_jstart &
,my_jrow_start &
,my_jrow_end 
!
real(kind=kfpt),dimension(ims:ime,jms:jme,lm),intent(in) :: &
 field
!
real(kind=kfpt),dimension(ids:ide,jstart_fft:jend_fft,1:lm_fft),intent(out) :: &
 array_lyrs
!
logical(kind=klog),dimension(0:npes-1),intent(in) :: &
 my_domain_has_fft_lats
!
!-----------------------------------------------------------------------
!***  Local Variables
!-----------------------------------------------------------------------
!
integer(kind=kint) :: &
 i &
,i_extent &
,iend &
,ierr &
,istart &
,j &
,j_extent &
,jend_fft_local &
,jstart_fft_local &
,k &
,k1 &
,k2 &
,k_extent &
!!!,mype &
,n &
,nbuf &
,npe 
!
integer(kind=kint),dimension(mpi_status_size) :: jstat
integer(kind=kint),dimension(:),allocatable,target,save :: handle
!
real(kind=kfpt),dimension(:),allocatable,save :: &
 dummy_recv 
real(kind=kfpt),dimension(:,:),allocatable,save :: &
 dummy_send
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!***  Handles are needed for the mpi_waits associated with mpi_isend.
!***  Allocate the handle array and initialize it.
!-----------------------------------------------------------------------
!
      if(.not.allocated(handle))then
        allocate(handle(0:npes-1))
        do n=0,npes-1
          handle(n)=mpi_request_null
        enddo
        allocate(dummy_recv(msize_dummy_fft))
        if(my_domain_has_fft_lats(mype))then
          allocate(dummy_send(msize_dummy_fft,ipe_start:ipe_end))
        endif
      endif
!
!-----------------------------------------------------------------------
!***  Each task is responsible for multiple full layers (the entire
!***  horizontal expanse of points in a given model layer in a given
!***  hemisphere that are on latitude circles where FFTs are applied),
!***  a single layer, or even a partial layer if there are more MPI
!***  compute tasks than there are layers (by the definition above,
!***  there are 2*LM layers).
!***  Each task must gather full latitude circles for its designated
!***  layers and rows and send to the other tasks its pieces of the
!***  latitude circles they need for the model layers they will be
!***  handling.
!
!***  k1_fft and k2_fft provide the vertical limits for each task's
!***  group of model layers.
!
!***  We gather into array array_lyrs which will hold full latitude
!***  circles of data in each task's own set of assigned model layers.
!
!-----------------------------------------------------------------------
!***  Compute the number of words sent by mype to everyone else.
!***  A task will send only if it has FFT latitude circles
!***  within a receiving task's designated region for computing FFTs.
!***  The sender uses the receiver's K extent, its own I extent,
!***  and the J extent for which its subdomain has FFT latitudes.
!-----------------------------------------------------------------------
!
      if(my_domain_has_fft_lats(mype))then
!
        send_to_npe: do npe=ipe_start,ipe_end     !<--- Send subsets of my FFT points to all
                                                  !     other tasks in this hemisphere.
!
          nbuf=0                                  !<--- This counts the words we are inserting
                                                  !     into the send buffer.
          jstart_fft_local=max(jts,my_jrow_start(npe))
          jend_fft_local=min(jte,my_jrow_end(npe))
!
          if(jstart_fft_local<=jend_fft_local)then
!
            k1=k1_fft(npe)
            k2=k2_fft(npe)
            do k=k1,k2
              do j=jstart_fft_local,jend_fft_local
              do i=its,ite
                nbuf=nbuf+1                       !<-- This counts words going to all tasks.
                dummy_send(nbuf,npe)=field(i,j,k)
              enddo
              enddo
            enddo
!
!           call mpi_issend(dummy_send(1,npe),nbuf,mpi_real,     &
            call mpi_isend(dummy_send(1,npe),nbuf,mpi_real,      &
                            npe,mype,mpi_comm_comp,              &
                            handle(npe),ierr)
!
          endif
!
        enddo send_to_npe
!
      endif

!-----------------------------------------------------------------------
!***  Compute the number of words received by mype from everyone.
!***  Only those tasks containing FFT latitude circles will be
!***  sending more than zero words.
!***  The receiver uses its K extent, the sender's full I extent,
!***  and the sender's J extent over which its subdomain has
!***  FFT latitude circles.
!***  Everyone receives since FFT work is shared by all.
!
!------------------------------------------------------------------------
!***  Loop through the tasks in each hemisphere from which mype receives.
!------------------------------------------------------------------------
!
!***  We need to specify which points each task will receive from tasks with FFT latitudes.
!
      k1=k1_fft(mype)
      k2=k2_fft(mype)
      k_extent=k2-k1+1
!
      recv_from_npe: do npe=ipe_start,ipe_end    !<--- Tasks with FFT lats in this task's
                                                 !     hemisphere will send
!
        from_senders: if(my_domain_has_fft_lats(npe))then
          jstart_fft_local=max(local_jstart(npe),my_jrow_start(mype))
          jend_fft_local=min(local_jend(npe),my_jrow_end(mype))
          j_extent=jend_fft_local-jstart_fft_local+1
          if(j_extent<=0)cycle
          istart=local_istart(npe)
          iend=local_iend(npe)
          i_extent=iend-istart+1
          n=j_extent*k_extent*i_extent      !<--- Total # of 3-D FFT points coming from remote
                                            !     task
!
!-----------------------------------------------------------------------
!***  Receive data only from tasks with FFT latitudes
!-----------------------------------------------------------------------
!
          call mpi_recv(dummy_recv,n,mpi_real &
                       ,npe,npe &
                       ,mpi_comm_comp &
                       ,jstat &
                       ,ierr)
!
!-----------------------------------------------------------------------
!***  Fill in the working array.
!***  Use the horizontal domain information from each source task that
!***  contained the original pieces of FFT circles since what we are
!***  doing here is combining all those pieces for particular layers
!***  into full circles.
!-----------------------------------------------------------------------
!
!
          nbuf=0      !<---  This counts the words we are receiving from the recv buffer.
!
          n=0
          do k=k1,k2
             n=n+1
             do j=jstart_fft_local,jend_fft_local
             do i=istart,iend
                nbuf=nbuf+1
                array_lyrs(i,j,n)=dummy_recv(nbuf)
             enddo
             enddo
          enddo
!
        endif from_senders
!
      enddo recv_from_npe
!
!-----------------------------------------------------------------------
!***  Clear the ISend request handles.
!-----------------------------------------------------------------------
!
      if(my_domain_has_fft_lats(mype))then
!
        do npe=ipe_start,ipe_end  
          call mpi_wait(handle(npe),jstat,ierr)
        enddo
!
      endif
!
!-----------------------------------------------------------------------
!
      end subroutine gather_layers
!
!-----------------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!-----------------------------------------------------------------------
                        subroutine scatter_layers &
      (array_lyrs,lm,npes,msize_dummy_fft &
      ,lm_fft,k1_fft,k2_fft &
      ,local_istart,local_iend &
      ,local_jstart,local_jend &
      ,jstart_fft,jend_fft &
      ,my_jrow_start,my_jrow_end &
      ,ipe_start,ipe_end &
      ,my_domain_has_fft_lats &
      ,mype,mpi_comm_comp &
      ,its,ite,ims,ime,ids,ide,jts,jte,jms,jme &
      ,field)
!-----------------------------------------------------------------------
!***  SCATTER_LAYERS distributes the elements of array_lyrs between 
!***  layers k1 and k2 inclusive to the appropriate tasks that actually
!***  own the FFT latitude rows.
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
integer(kind=kint),intent(in) :: &
 its &
,ite &
,ims &
,ime &
,ids &
,ide &
,jts &
,jte &
,jms &
,jme
!
integer(kind=kint),intent(in) :: &
 ipe_end &
,ipe_start &
,jstart_fft &
,jend_fft &
,lm &
,lm_fft &
,mpi_comm_comp &
,msize_dummy_fft &
,mype &
,npes 
!
integer(kind=kint),dimension(0:npes-1),intent(in) :: &
 k1_fft &
,k2_fft &
,local_iend &
,local_istart &
,local_jend &
,local_jstart &
,my_jrow_start &
,my_jrow_end
!
real(kind=kfpt),dimension(ids:ide,jstart_fft:jend_fft,1:lm_fft),intent(in) :: &
 array_lyrs
!
real(kind=kfpt),dimension(ims:ime,jms:jme,lm),intent(out) :: &
 field
!
logical(kind=klog),dimension(0:npes-1),intent(in) :: &
 my_domain_has_fft_lats
!
!-----------------------------------------------------------------------
!***  Local Variables
!-----------------------------------------------------------------------
!
integer(kind=kint) :: &
 i &
,i_extent &
,iend &
,ierr &
,iprod &
,istart &
,j &
,j_extent &
,jend_fft_local &
,jstart_fft_local &
,k &
,k_extent &
,k1 &
,k2 &
!!!,mype &
,n &
,ncount_recv &
,nbuf &
,nn &
,npe
!
integer(kind=kint),dimension(mpi_status_size) :: jstat
integer(kind=kint),dimension(:),allocatable,target,save :: handle
!
real(kind=kfpt),dimension(:),allocatable,save :: &
 dummy_recv 
real(kind=kfpt),dimension(:,:),allocatable,save :: &
 dummy_send
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
!xxx  mype=mype_share
!
!-----------------------------------------------------------------------
!***  Handles are needed for the mpi_waits associated with mpi_isend.
!***  Allocate the handle array and initialize it.
!-----------------------------------------------------------------------
!
      if(.not.allocated(handle))then
        allocate(handle(0:npes-1))
        do n=0,npes-1
          handle(n)=mpi_request_null
        enddo
!
        allocate(dummy_recv(msize_dummy_fft))
        allocate(dummy_send(msize_dummy_fft,ipe_start:ipe_end))
      endif
!
!-----------------------------------------------------------------------
!***  Each task holds full latitude circles of data within its own
!***  subset of model layers that it was assigned.  We want to
!***  redistribute the data back to individual tasks whose domains
!***  actually contain the particular points.
!
!***  k1_fft and k2_fft provide the vertical limits for each task's
!***  group of model layers that it will use for FFT computation.
!***  In the general case, these layers and the points in them will
!***  come from other tasks.
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!***  Compute the number of words sent by mype to everyone.  A task
!***  will only send to tasks that have FFT latitude circles that
!***  were specified to be handled by the sending task.
!***  The sender uses its own K extent and the receiver's I extent
!***  and the receiver's J extent that spans the appropriate
!***  FFT latitude rows computed by the sender.
!-----------------------------------------------------------------------
!
      k1=k1_fft(mype)
      k2=k2_fft(mype)
      k_extent=k2-k1+1
!
      send_to_npe: do npe=ipe_start,ipe_end       !<-- Tasks with FFT lats in this task's
                                                  !    hemisphere will receive
!
        if(my_domain_has_fft_lats(npe))then       !<--- Send back to only those tasks with FFT
                                                  !     lats in this hemisphere
          istart=local_istart(npe)
          iend=local_iend(npe)
          i_extent=iend-istart+1
          jstart_fft_local=max(local_jstart(npe),my_jrow_start(mype))
          jend_fft_local=min(local_jend(npe),my_jrow_end(mype))
          j_extent=jend_fft_local-jstart_fft_local+1
!
          nn=0                                    !<--- Counter for number of words sent to
                                                  !     each task
          if(j_extent>0)then
!
            n=0
            nbuf=0                                !<--- Counter for number of words inserted 
                                                  !     into send buffer
!
            do k=k1,k2
              n=n+1
              do j=jstart_fft_local,jend_fft_local
              do i=istart,iend
                nn=nn+1
                nbuf=nbuf+1
                dummy_send(nbuf,npe)=array_lyrs(i,j,n)
              enddo
              enddo
            enddo
!
!           call mpi_issend(dummy_send(1,npe),nbuf,mpi_real             &
            call mpi_isend(dummy_send(1,npe),nbuf,mpi_real              &
                           ,npe,npe,mpi_comm_comp                       &
                           ,handle(npe),ierr)
          endif
!
        endif
!
      enddo send_to_npe
!
!-----------------------------------------------------------------------
!***  Each task with FFT latitudes now fills in its output array.
!***  The data received from remote tasks is sized to the remote
!***  tasks' FFT computation layers and the local tasks' I extent
!***  and J extent that lie within FFT latitudes.
!-----------------------------------------------------------------------
!
      if(my_domain_has_fft_lats(mype))then
!
!-----------------------------------------------------------------------
!
        i_extent=ite-its+1
!
        recv_from_npe: do npe=ipe_start,ipe_end       !<---  Tasks with FFT lats recv from 
                                                      ! evryone in the hemisphere
                                                      ! who has computed any FFTs within 
                                                      ! this task's FFT lats.
          jstart_fft_local=max(jts,my_jrow_start(npe))
          jend_fft_local=min(jte,my_jrow_end(npe))
          if(jstart_fft_local>jend_fft_local)cycle
          j_extent=jend_fft_local-jstart_fft_local+1
!
          iprod=max(i_extent*j_extent,0)
!
          k1=k1_fft(npe)
          k2=k2_fft(npe)
          k_extent=k2-k1+1
!
          ncount_recv=iprod*k_extent               !<--  # of words to receive from task npe
!
          nbuf=0                                   !<---  Counter for words in recv buffer
!
          call mpi_recv(dummy_recv,ncount_recv,mpi_real               &
                       ,npe,mype                                      &
                       ,mpi_comm_comp                                 &
                       ,jstat                                         &
                       ,ierr)
!
          do k=k1,k2
            do j=jstart_fft_local,jend_fft_local
            do i=its,ite
              nbuf=nbuf+1
              field(i,j,k)=dummy_recv(nbuf)
            enddo
            enddo
          enddo
!
        enddo recv_from_npe
!
      endif
!
!-----------------------------------------------------------------------
!***  Clear the ISend request handles.
!-----------------------------------------------------------------------
!
      do npe=ipe_start,ipe_end
!
        if(my_domain_has_fft_lats(npe))then       !<--- Send back to only those tasks with FFT
          call mpi_wait(handle(npe),jstat,ierr)
        endif
!
      enddo
!
!-----------------------------------------------------------------------
!
      end subroutine scatter_layers
!
!-----------------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!-----------------------------------------------------------------------
!
!===============================================================
!  From Jim Abeles.
!  Block partition the loop bounds (lb...ub) -> (i1...i2).
!  The number of tasks is ntasks;  taskid = 0, 1, ..., ntasks-1.
!  The first nt1 tasks get a chunk one bigger than the rest.
!  The counts and displacements arrays range from 1 to ntasks.
!===============================================================
!
      subroutine looplimits(taskid, ntasks, lb, ub, i1, i2)
      implicit none
      integer taskid, ntasks, lb, ub, i1, i2
      integer chunk, nwork, nt1, nt2
      integer itask, netdisp
      integer counts(ntasks), displacements(ntasks)

      nwork = ub - lb + 1
      chunk = nwork/ntasks
      nt1 = nwork - ntasks*chunk
      nt2 = ntasks - nt1

      netdisp = lb
      do itask = 1, nt1
         counts(itask) = chunk + 1
         displacements(itask) = netdisp  
         netdisp = min(ub,netdisp+chunk+1)
      end do
      do itask = nt1 + 1 , ntasks
         counts(itask) = chunk
         displacements(itask) = netdisp  
         netdisp = min(ub,netdisp+chunk)
      end do

      i1 = displacements(taskid+1)
      i2 = min(ub,i1+counts(taskid+1)-1)

      return
      end subroutine looplimits
!
!-----------------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!-----------------------------------------------------------------------
!
                       end module module_dm_parallel
!
!-----------------------------------------------------------------------
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!-----------------------------------------------------------------------
