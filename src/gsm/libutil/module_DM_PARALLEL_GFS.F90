!-----------------------------------------------------------------------
                        module module_dm_parallel_gfs
!-----------------------------------------------------------------------
!
! 09/2009  Weiyu Yang   -  Ensemble GEFS.
!----------------------------------------
!
!***  This module contains all codes directly related to distributed
!***  memory issues except for halo exchange although note that the
!***  halo widths must be set here.
!
!-----------------------------------------------------------------------
!
      use module_gfs_machine
      use module_gfs_mpi_def,only : mc_comp                      &
                                   ,mpi_comm_comp                &
                                   ,mpi_comm_inter               &
                                   ,num_pes_fcst
!
!-----------------------------------------------------------------------
!
      implicit none
      public setup_servers_gfs
!
!
!-----------------------------------------------------------------------
!---domain decomposition info-------------------------------------------
!-----------------------------------------------------------------------
!
       contains
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
      subroutine setup_servers_gfs(mype,npes,last_fcst_pe               &
                              ,ngroups_write,write_tasks_per_group      &
                              ,mpi_intra,max_inter_groups               &
                              ,mpi_comm_inter_array)
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
!    mype          - My task ID
!    npes          - Total number of mpi tasks provided to the job.  As input 
!                    to SETUP_SERVERS it includes the forecast tasks plus all
!                    write tasks in all groups of write tasks.  npes must at least
!                    equal the forecast pes otherwise the integration
!                    cannot proceed. The difference between the product npes_fcst
!                    and npes is the number of mpi tasks that are available
!                    for i/o serving. This can be zero, in which case output will
!                    write a direct access file that can be separately quilted. 
!                    In order to skip the separate quilting step, make sure that
!                    the number of mpi tasks that the code is initiated with is at
!                    least one greater than npes_fcst.
!                    Later in the routine, npes is reset to the number of fcst tasks.
!    ngroups_write - Number of groups of write tasks.
!    mpi_intra     - The global communicator.
!    write_tasks_per_group - # of write tasks per write group
!
!-----------------------------------------------------------------------
!***  Argument variables.
!-----------------------------------------------------------------------
!
      integer(kind=kind_io4),intent(in) :: &
 mype &                     ! each task's ID
,last_fcst_pe  &            ! last fcst pe ID
,ngroups_write &            ! number of groups of write tasks
,write_tasks_per_group &    ! number of groups of write tasks per group
,mpi_intra     &            ! global communicator
,max_inter_groups           ! max number of quilt server groups
!
      INTEGER,dimension(max_inter_groups),intent(inout)  ::            &
 mpi_comm_inter_array       ! quilt server group mpi communicators
      integer(kind=kind_io4),intent(inout) :: &
 npes                       ! total number of tasks provided
                            ! then converted to the number of fcst tasks
!-----------------------------------------------------------------------
!***  Local variables.
!-----------------------------------------------------------------------
!
      integer(kind=kind_io4) :: comdup,i,icc,icolor,iendq,iendxx,ierr &
                           ,igroup,igroup_x,iqserver,iquilt_group &
                           ,iss,issl &
                           ,istaq,istaxx,iworld,iworld_minus,ixx,jj,kk & 
                           ,lead_remote,npes_fcst,one,last_qserver
!
      integer(kind=kind_io4),allocatable,dimension(:) :: irank         &
                           ,num_serv_per_grp
!
      logical :: include_mype
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!
!***  This is the number of mpi tasks the executable has been built for.
!***  npes, returned from mpi_comm_size, must be at least the npes_fcst size
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
      npes_fcst=num_pes_fcst
      mc_comp=mpi_intra                                            !<-- Set mpi_comm_comp to the global communicator
      mpi_comm_comp=mpi_intra 
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
!***  specified for the forecast integration (npes_fcst).
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
      allocate(num_serv_per_grp(iquilt_group))
!-----------------------------------------------------------------------
!
!***  ERROR CHECK FOR NUMBER OF GROUPS - MAXIMUM IS 100 - THAT IS ALOT!
!
!-----------------------------------------------------------------------
      if(iquilt_group>100)then
        write(0,*)' ***** IQUILT_GROUP IS GREATER THAN 100'
        write(0,*)' ***** DO YOU REALLY WANT THIS ?'
        write(0,*)' ***** IF SO THEN INCREASE SIZE IN mpp.h'
        write(0,*)' ***** ALSO, CHANGE IF CHECK IN SETUP_SERVERS'
        write(0,*)' ***** RESETTING THE NUMBER OF SERVER GROUPS TO 100'
        write(0,*)' ***** WE ARE CONTINUING ....   '
        iquilt_group=max_inter_groups
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
!junwang
       write(0,*)'in setup_gfs,last_fcst_pe=',last_fcst_pe,'ngroup_write=', &
         ngroups_write,'iqserver=',iqserver
      last_qserver=last_fcst_pe+iqserver
!
      if(iqserver==0)then
        if(mype==0)then
          write(0,*)' *** You specified 0 Write tasks '
          write(0,*)' *** Output will write a direct access file'
        endif
        iquilt_group=0
      endif
!
      if(iqserver/=0 .and. iquilt_group>iqserver)then
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
      if(mype<=last_fcst_pe)then
        icolor=0
      elseif(mype<=last_qserver)then
        istaxx=last_fcst_pe+1
        do i=1,iquilt_group
          iendxx=istaxx+num_serv_per_grp(i)-1
          if(mype>=istaxx.and.mype<=iendxx)then
            icolor=i
          endif
          istaxx=iendxx+1
        enddo
      endif
      write(0,*)'my color=',icolor
!
!-----------------------------------------------------------------------
!***  Split the communicator - The new intracommunicator for all tasks
!***  is mpi_comm_comp. mpi_comm_world (mpi_intra) is still available but it 
!***  refers to all the mpi tasks (model integration AND quilt tasks).
!-----------------------------------------------------------------------
!        
!!!   call mpi_comm_dup(mpi_comm_world,comdup,ierr)
      call mpi_comm_dup(mpi_intra,comdup,ierr)
!
!-----------------------------------------------------------------------
!### jw
!***  if ntasks >0 then mmunicator - The new intracommunicator for all tasks
!###
!-----------------------------------------------------------------------
!
      call mpi_comm_split(comdup,icolor,mype,mc_comp,ierr)
      mpi_comm_comp=mc_comp
      write(0,*)'after mpi_comm_split'
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
                                                                         !<-- Dimension irank to the total # of quilt tasks
      if(iqserver>0)then
        do i=1,iqserver
          irank(i)=-1                                                    !<-- Initialize irank to meaningless values
        enddo
      endif
!
!junwang      ixx=npes_fcst                                                      !<-- Let ixx be the # of forecast tasks
      ixx=last_fcst_pe+1                                                    !<-- Let ixx be the # of forecast tasks
!
!-----------------------------------------------------------------------
!
      inter_comm : do i=1,iquilt_group                                   !<-- Create intercommunicators between the set of fcst tasks
                                                                         !    and each quilt group individually.
!-----------------------------------------------------------------------
!
        include_mype=.true.
!
!junwang        if(mype<npes_fcst)then
        if(mype<=last_fcst_pe)then
          lead_remote=ixx                                                !<-- All fcst tasks set lead_remote to total # of fcst tasks.
!!!       lead_remote=ixx+(i-1)*write_tasks_per_group                    !<-- All fcst tasks set lead_remote to total # of fcst tasks.
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
!junwang        iss=npes_fcst
        iss=last_fcst_pe+1
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
        write(0,*)'icc=',icc,'irank=',irank,'lead_remote=',lead_remote,  &
         'include_mype=',include_mype,'npes_fcst=',npes_fcst
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
        write(0,*)'after mpi_comm_create'
!
!-----------------------------------------------------------------------
!***  At this point we have an intracommunicator that excludes the tasks we don't want.
!***  Create an intercommunicator for use between the mpi tasks doing the model
!***  integration and quilt group i that we are considering. This process is
!***  a collective routine that can only be done by the tasks that have not
!***  been excluded. Save this new intercommunicator in mpi_comm_inter for use
!***  by the tasks that belong to the quilt group that we are considering. The
!***  tasks that are performing the model integration will reference
!***  mpi_comm_inter_array() since they will need to select which quilt
!***  group they wish to communicate with.
!-----------------------------------------------------------------------
!
        if(include_mype.and.npes_fcst>0)then                                              !<-- Select all fcst tasks plus the quilt tasks in group i
          write(0,*)'before mpi_intercomm_create'
          call mpi_intercomm_create(mc_comp                             & !<-- Each task's local intracommunicator
                                   ,0                                   & !<-- Rank of lead task in each local intracommunicator
                                   ,iworld_minus                        & !<-- Intracommunicator between fcst tasks and
                                                                          !    quilt tasks in group i
                                   ,lead_remote                         & !<-- Rank of leader in the remote group in iworld_minus
                                   ,0                                   & !<-- A tag
                                   ,mpi_comm_inter_array(i)             & !<-- The new intercommunicator between the fcst tasks and
                                                                          !    the tasks in quilt group i
                                   ,ierr)
          write(0,*)'after mpi_intercomm_create,mpi_comm(',i,')=',   &
            mpi_comm_inter_array(i)
           mpi_comm_inter=mpi_comm_inter_array(i)
        endif
!
        call mpi_barrier(mpi_intra,ierr)
!
!-----------------------------------------------------------------------
      enddo  inter_comm
!------------------------------------------------------------------------
!
      deallocate (irank)
!-----------------------------------------------------------------------
!### jw
!***  if ntasks >0 then mmunicator - The new intracommunicator for all tasks
!###
!-----------------------------------------------------------------------
!
!***
!***  Set npes to the number of tasks working on the model integration.
!***
      npes=npes-iqserver
      write(0,*)'npes=',npes,'last_fcst_pe=',last_fcst_pe
!
!!!   if(mype==0)then
        write(0,*)' Number of Integration Tasks: ',npes
        write(0,*)' Total Number of Quilt Servers: ',iqserver
        write(0,*)' exit setup_servers npes=',npes
!!!   endif
!***
      deallocate (num_serv_per_grp)
!-----------------------------------------------------------------------
!
      end subroutine setup_servers_gfs
!
!-----------------------------------------------------------------------
!
      end module module_dm_parallel_gfs
