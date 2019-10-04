module module_core_setup

!JR Left undesirable name "module_core_setup" as is because svn external
!JR icosio/ files "use" this module, and changing the name would require
!JR separate commits to FIM, icosio, and NIM repositories.
!JR Prefer "use mpi" to "include 'mpif.h'" for argument checking, but that 
!JR implementation needs to wait till SMS can provide a declaration that says 
!JR it doesn't need to look at the module file.
!JR Moved mpi use/include to top of module instead of inside included routines
!JR so can easily swap between "use" and "include" as desired. 

!  use mpi
  use module_wtinfo,only:wtinfo

  implicit none
  save
  private

  include 'mpif.h'

  integer, parameter :: badret = 1 ! Return code for MPI_Abort to pass to environment

! Public stuff

  public :: core_setup_fim, iam_compute_root, iam_write_root
  integer, public :: my_comm=MPI_COMM_NULL          ! MPI intra-communicator for FIM compute or write tasks.
  integer, public :: output_intercomm=MPI_COMM_NULL ! MPI inter-communicator between FIM compute and write tasks.
  integer, public :: total_comm=MPI_COMM_NULL       ! MPI intra-communicator for all available MPI tasks.
  integer, public :: world_rank                     ! rank of this task in mpi_comm_in or mpi_comm_world, public for debugging only
  integer, public :: nct=-1                         ! number of compute tasks
  integer, public :: nwt=-1                         ! number of write tasks

! TODO Ensure use_write_tasks, iam_fim_task, iam_write task dont get used before
! core_setup_fim is called
  logical, public :: use_write_tasks = .false. ! .true. iff write tasks are enabled (init to false)
  logical, public :: iam_fim_task = .true.     ! .true. iff I am a FIM compute task (init to true)
  logical, public :: iam_write_task = .false.  ! .true. iff i am a FIM write task (init to false)
  logical, public :: core_setup_done = .false. ! .true. iff core_setup_fim has finished without error

CONTAINS


! public interfaces


!*********************************************************************
  subroutine core_setup_fim (mpi_comm_in,tasklist_compute,tasklist_write)
!
! SUMMARY:  
!       When MPI is used, set up communicators for FIM compute tasks and 
!       optional write tasks.  
!       When MPI is not used, return immediately.  
! NOTE:  Executable SMS directives must not be placed before 
! NOTE:  this routine is called!  
! NOTE:  This includes writes or prints without !SMS$ignore because they 
! NOTE:  cause SMS to generate code.
!
! ARGUMENTS:  
!       Optional argument mpi_comm_in is an MPI communicator that may be 
!       passed in to restrict FIM to a subset of MPI_COMM_WORLD. 
!
!       If present, optional argument tasklist_compute must be an 
!       unassociated pointer to a rank-1 integer array.  It will be 
!       allocated and filled with ranks of compute tasks in MPI_COMM_WORLD 
!       or, if present, mpi_comm_in.  It is the responsibility of the 
!       caller to deallocate tasklist_compute.  
!
!       If present, optional argument tasklist_write must be an 
!       unassociated pointer to a rank-1 integer array.  It will be 
!       allocated and filled with ranks of write tasks in MPI_COMM_WORLD 
!       or, if present, mpi_comm_in.  If the number of write tasks requested 
!       is 0, tasklist_write will be allocated with size 0.  It is the 
!       responsibility of the caller to deallocate tasklist_write.   
!
! DETAILS:  
!       Split the MPI communicator and create intercommunicators.  
!     MPI tasks may be split into two or three groups depending upon 
!     settings in FIMnamelist.  The first group will contain the 
!     "compute" tasks responsible for model computations.  The 
!     optional second group will contain "write" tasks responsible 
!     for writing output to disk allowing overlap of computation with 
!     output.  If needed, a third group will contain "do-nothing" tasks 
!     that are neither compute tasks nor write tasks.  The use of 
!     "do-nothing" tasks allows the compute root to optionally live 
!     on its own node (if more memory is needed).  The "do-nothing" 
!     tasks also allow write tasks to be mapped to nodes optionally 
!     leaving some cores unused (if more memory is needed or to allow 
!     concurrent writes of output files from multiple nodes which is 
!     very fast on some machines).  Use of "do-nothing" tasks avoids 
!     dependence upon site-specific features of batch queuing systems 
!     to map MPI tasks onto nodes (e.g. mpich $MACHINE_FILE).  
!       Some run time options and their associated FIMnamelist settings 
!     are listed below.  When specified in FIMnamelist, "cpn" is the 
!     number of cores per node.  
!       -------------------------------------------------------------
!       Case 1:  ComputeTasks = N
!                root_own_node = .false.
!                num_write_tasks = 0
!         0                               compute "root"
!         1..N-1                          compute tasks
!       -------------------------------------------------------------
!       Case 2:  ComputeTasks = N
!                root_own_node = .true.
!                cpn = cn
!                num_write_tasks = 0
!         0                               compute "root"
!         1..cn-1                         do-nothing tasks
!         cn..cn-1+N-1                    compute tasks
!       -------------------------------------------------------------
!       Case 3:  ComputeTasks = N
!                root_own_node = .false.
!                cpn = cn
!                num_write_tasks = 1
!         0                               compute "root"
!         1..cn-1                         compute tasks
!         cn                              write task
!         cn+1..(2*cn)-1                  do-nothing tasks
!         (2*cn)..cn+N-1                  compute tasks
!       -------------------------------------------------------------
!       Case 4:  ComputeTasks = N
!                root_own_node = .true.
!                cpn = cn
!                num_write_tasks = 1
!         0                               compute "root"
!         1..cn-1                         do-nothing tasks
!         cn                              write task
!         cn+1..(2*cn)-1                  do-nothing tasks
!         (2*cn)..(2*cn)-1+N-1            compute tasks
!       -------------------------------------------------------------
!       Case 5:  ComputeTasks = N
!                root_own_node = .true.
!                cpn = cn
!                num_write_tasks = 21
!                max_write_tasks_per_node = 3
!         0                               compute "root"
!         1..cn-1                         do-nothing tasks
!         cn                              write "root"
!         cn+1..cn+(3-1)                  write tasks
!         cn+3..(2*cn)-1                  do-nothing tasks
!         2*cn..(2*cn)+(3-1)              write tasks
!         (2*cn)+3..(3*cn)-1              do-nothing tasks
!         3*cn..(3*cn)+(3-1)              write tasks
!         (3*cn)+3..(4*cn)-1              do-nothing tasks
!         4*cn..(4*cn)+(3-1)              write tasks
!         (4*cn)+3..(5*cn)-1              do-nothing tasks
!         5*cn..(5*cn)+(3-1)              write tasks
!         (5*cn)+3..(6*cn)-1              do-nothing tasks
!         6*cn..(6*cn)+(3-1)              write tasks
!         (6*cn)+3..(7*cn)-1              do-nothing tasks
!         7*cn..(7*cn)+(3-1)              write tasks
!         (7*cn)+3..(8*cn)-1              do-nothing tasks
!         (8*cn)..(8*cn)-1+N-1            compute tasks
!       -------------------------------------------------------------
!
!TODO:  We put compute root on the first node to avoid problems with file 
!TODO:  system sync on some machines.  For example, if compute root is not 
!TODO:  on node zero, then file system changes made by prep, which does run 
!TODO:  on node zero, may not be up to date by the time the compute root runs 
!TODO:  on another node.  This approach may need fine-tuning as it may cause 
!TODO:  the "node affinity" features on other machines to go haywire.  
!
!*********************************************************************

!TODO:  At present, this routine ignores serial vs. parallel build and 
!TODO:  always calls MPI routines.  In practice this is not a problem 
!TODO:  because MPI is always linked in even for serial builds due to use of 
!TODO:  the MPI timers.  If this becomes a problem, wrap this code in cpp 
!TODO:  directives.  

!TODO:  MPI communicators created here are never freed.  Fix.  

    IMPLICIT NONE

! Arguments
    INTEGER, OPTIONAL, INTENT(IN) :: mpi_comm_in
    INTEGER, OPTIONAL, POINTER :: tasklist_compute(:)
    INTEGER, OPTIONAL, POINTER :: tasklist_write(:)

! Local declarations
    logical, parameter :: debugon=.false.    ! set to .true. to turn on debug prints
    integer, parameter :: tag = 998          ! tag for MPI_Recv
    integer, parameter :: intercomm_tag = 15 ! tag for intercommnicator creation

    REAL*8  :: t0,t1=0.0d0     ! local timers
    LOGICAL :: initialized     ! whether mpi_init has been called
    INTEGER :: my_rank         ! rank of this task in total_comm
    INTEGER :: num_tasks       ! total number of MPI tasks available
    INTEGER :: num_tasks_world ! total number of MPI tasks in mpi_comm_in or mpi_comm_world
    INTEGER :: mwtpn           ! max write tasks per node
    INTEGER :: color           ! input arg to mpi_comm_split()--must be non-negative

    integer :: cpn             ! cores per node
    integer :: status(mpi_status_size) ! returned from mpi_recv
    integer :: ignore          ! return code from mpi routines (ignored)
    integer :: ierr            ! return code from mpi routines (not ignored)
    integer :: coresleft       ! number of cores to be allocated after root, write_tasks
    integer :: i               ! loop index over ranks in input communicator
    integer :: ii              ! loop index over ranks in fimcomm_tweaked
    integer :: i2              ! another loop index
    integer :: inc             ! increment
    integer :: looplim         ! ending loop index
    integer :: begwriterank = 999999  ! rank of first write task in fimcomm_tweaked
    integer :: endwriterank = -999999 ! rank of last write task in fimcomm_tweaked
    integer :: n               ! index over write tasks
    integer :: iamnew          ! Rank in new communicator

    logical :: root_own_node   ! whether root has a node to himself
! abort_on_bad_task_distrib controls whether to abort the FIM run if something fishy is 
! encountered w.r.t. node names associated with MPI tasks
    logical :: abort_on_bad_task_distrib ! whether to die if cpn is wrong

    character(len=mpi_max_processor_name) :: mynode ! node name
    character(len=mpi_max_error_string)   :: string ! error string
    integer :: resultlen        ! length return from MPI routines
    integer :: tmp_comm_in      ! input communicator (MPI_COMM_WORLD or derivative)
    integer :: tmp_comm         ! dup of tmp_comm_in
    integer :: fimcomm_tweaked  ! tmp_comm modified to contain only MPI ranks which map to 
! compute tasks or write tasks, i.e. what "tweak_hostfile" 
! used to do--leave out cores which will not participate

! Dynamic arrays

! proclist will contain the list of MPI ranks from the input communicator which will 
! define a new communicator (fimcomm_tweaked). This new communicator eliminates "holes"
! to allow for such things as the root MPI task having a node to himself, and write tasks
! starting on a node boundary.
    integer, allocatable :: proclist(:)
    integer :: lastproclist       ! last valid index of proclist
    integer :: numcomputetasks    ! number of compute tasks
    logical :: debugmsg_on        ! write-task debug message control
    character(len=mpi_max_processor_name), allocatable :: nodes(:) ! node names

! Return immediately iff this is a serial run (i.e. MPI_INIT has not been called)

    call mpi_initialized (initialized, ierr)
    if (.not.initialized) then
! this is a serial run
      nct = 1
! TODO:  remove duplication
      iam_fim_task = .true.
      iam_write_task = .false.
      core_setup_done = .true.
      return
    end if

    CALL StartTimer(t0)

!TODO:  For NEMS, read from config file instead of namelist.  Keep details 
!TODO:  inside get_write_task_info().  

    call wtinfo(cpn,nwt,mwtpn,root_own_node,abort_on_bad_task_distrib,debugmsg_on)
    if (nwt > 0 .and. mwtpn > cpn) then
!SMS$IGNORE BEGIN
      print *,'ERROR in core_setup_fim: wtinfo returned max write tasks per node=', &
        mwtpn, ' > cores per node=', cpn
      call mpi_abort (mpi_comm_world, badret, ignore)
      stop
!SMS$IGNORE END
    end if

    if (cpn < 1) then
!SMS$IGNORE BEGIN
      print *,'ERROR in core_setup_fim: wtinfo returned non-positive cpn=', cpn
      print *,'This MUST be set in /writetasknamelist/'
      call mpi_abort (mpi_comm_world, badret, ignore)
      stop
!SMS$IGNORE END
    end if

    IF (nwt < 0) THEN
!SMS$IGNORE BEGIN
      PRINT *,'ERROR in core_setup_fim: wtinfo returned negative num_write_tasks'
      CALL MPI_ABORT (MPI_COMM_WORLD, badret, ignore)
      STOP
!SMS$IGNORE END
    ENDIF

! use passed-in communicator if present, otherwise use MPI_COMM_WORLD
    tmp_comm_in = MPI_COMM_WORLD
    IF ( PRESENT( mpi_comm_in ) ) THEN
      tmp_comm_in = mpi_comm_in
    ENDIF
! dup is needed to make sends and receives below safe
    CALL MPI_COMM_DUP (tmp_comm_in, tmp_comm, ierr)
    IF (ierr /= 0) THEN
!SMS$IGNORE BEGIN
      PRINT *,'ERROR in core_setup_fim:  MPI_COMM_DUP(tmp_comm_in) returned ',ierr
      CALL MPI_ABORT (MPI_COMM_WORLD, badret, ignore)
      STOP
!SMS$IGNORE END
    ENDIF

!JR Define a new communicator which places root and write tasks appropriately
!JR i.e. maybe root on a node by himself, and likewise write tasks on their
!JR own node(s). 
!JR Check that cpn is consistent with the return from MPI_GET_PROCESSOR_NAME
!JR The checking code ASSUMES that tweak_hostfile or equivalent is NOT in effect

    mynode = ' '
    call mpi_get_processor_name (mynode, resultlen, ierr)
    if (ierr /= 0) then
!SMS$IGNORE BEGIN
      call mpi_error_string (ierr, string, resultlen, ignore)
      write(6,*)'Error from mpi_get_processor_name:', string(1:resultlen)
      call mpi_abort (mpi_comm_world, badret, ignore)
      stop
!SMS$IGNORE END
    end if

    call mpi_comm_size (tmp_comm, num_tasks_world, ignore)
    call mpi_comm_rank (tmp_comm, world_rank, ignore)

    allocate (nodes(0:2*cpn-1))                  ! Check 1st 2 nodes
    allocate (proclist(0:num_tasks_world-1))     ! Ranks start at 0
    proclist(:)=-1
    nodes(:) = ' '

! Gather node name info for 1st 2 nodes to check conformance with 
! expected node allocation

    looplim = min (num_tasks_world-1, 2*cpn-1)
    if (world_rank == 0) then             ! get node name info from other tasks
      nodes(0) = mynode
      do i=1,looplim
        call mpi_recv (nodes(i), mpi_max_processor_name, mpi_character, &
          i, tag, tmp_comm, status, ignore)
      end do
    else if (world_rank <= looplim) then  ! send my info to root
      call mpi_send (mynode, mpi_max_processor_name, mpi_character, &
        0, tag, tmp_comm, ignore)
    end if

! Root task ensures that node name changes for second node

    if (world_rank == 0 .and. num_tasks_world > cpn) then
      if (nodes(0) == nodes(cpn)) then
!SMS$IGNORE BEGIN
        write(6,*)'core_setup_fim: MPI task 0 name=', trim (nodes(0)), &
          ' matches task ',cpn, ' name=', trim (nodes(cpn)),' cpn=', cpn
        write(6,*)'Perhaps namelist value for cpn=', cpn, ' is incorrect?'
        if (abort_on_bad_task_distrib) then
          write(6,*)'abort_on_bad_task_distrib is true so model is aborting...'
          call mpi_abort (mpi_comm_world, badret, ignore)
          stop
        end if
!SMS$IGNORE END
      end if
    end if

! When task distribution leaves "holes" (for example when 
! root_own_task is set on systems with more than one CPU per 
! node) some MPI tasks will not be used.  These relaxed and 
! unmotivated tasks are known as "do-nothing" tasks.  

! First: Set proclist for first node. Contents depend on whether root 
! has a node to himself

! I am a "do-nothing" task unless modified below
    iam_fim_task=.false.
    iam_write_task=.false.
    proclist(0) = 0  ! Root task always same in both communictors
    if (world_rank==0) iam_fim_task = .true.  ! compute task
    if (root_own_node) then
      i = cpn         ! Move pointer to 1st core of next node
      ii = 1          ! Set proclist index for next entry
    else              ! Fill remainder of node 0 with compute tasks
! TODO: Fix this for pathological case where root_own_node=.false. and nct < cpn
! TODO: Currently utils/get_num_cores guards against this case.
! TODO: Could be done here by calling GetNprocs, but that routine is INSANELY expensive because
! TODO: it opens 2 files and reads namelists from them in order to get the value.
      looplim = min (num_tasks_world-1, cpn-1)
      ii = 1          ! Set proclist index for next entry
      do i=1,looplim
        proclist(ii) = i
        if (world_rank==i) iam_fim_task = .true.  ! compute task
        if (world_rank == 0 .and. nodes(i-1) /= nodes(i)) then
!SMS$IGNORE BEGIN
          write(6,*)'core_setup_fim: MPI task ',i-1,' name=',trim(nodes(i-1)), &
            ' does not match task ',i,' name=',trim(nodes(i)),' cpn=',cpn
          if (abort_on_bad_task_distrib) then
            write(6,*)'abort_on_bad_task_distrib is true so model is aborting...'
            call mpi_abort (mpi_comm_world, badret, ignore)
            stop
          end if
!SMS$IGNORE END
        end if
        ii = ii + 1
      end do
      i = cpn
    end if

! Add write tasks to proclist. "i" currently points to 1st core of 2nd node

    if (nwt > 0) then
      begwriterank = ii
    end if

    do n=1,nwt
      if (i >= num_tasks_world) then
!SMS$IGNORE BEGIN
        if (world_rank == 0) then
          write(6,*)'Ran out of MPI tasks assigning write task=',n,'. i=',i
        end if
        call mpi_abort (mpi_comm_world, badret, ignore)
        stop
!SMS$IGNORE END
      end if
      proclist(ii) = i
      if (world_rank==i) iam_write_task = .true.  ! write task
      endwriterank = ii
      if (mod (n, mwtpn) == 0) then  ! Move pointer to start of next node
        inc = cpn - mwtpn + 1
      else                           ! Haven't filled up the node with write tasks yet
        inc = 1
      end if
      i = i + inc
      ii = ii + 1
    end do

! Compute tasks after write tasks start at next empty node

    if (mod (i, cpn) /= 0) then
      coresleft = cpn - mod (i, cpn)
      i = i + coresleft
    end if

! Remainder of proclist is compute tasks 

    do while (i < num_tasks_world)
      proclist(ii) = i
      if (world_rank==i) iam_fim_task = .true.  ! compute task
      i = i + 1
      ii = ii + 1
    end do
    lastproclist = ii - 1

    if (nwt == 0) then
      use_write_tasks = .false.
    else
      use_write_tasks = .true.
    end if

! Split "do-nothing" tasks from other tasks.  
! Yes, two splits are needed so intercommunicators can be correctly 
! constructed.  
! Define the new communicator fimcomm_tweaked.  
! tasks with color==1 will be in the "do-nothing" group
    color = 1
    if (iam_fim_task.or.iam_write_task) color = 0
    CALL MPI_COMM_SPLIT (tmp_comm, color, world_rank, fimcomm_tweaked, ierr)
    if (ierr /= 0) then
!SMS$IGNORE BEGIN
      call mpi_error_string (ierr, string, resultlen, ignore)
      write(6,*)'Error in mpi_comm_split:', string(1:resultlen)
      call flush (6)
      call mpi_abort (mpi_comm_world, badret, ignore)
!SMS$IGNORE END
    endif

    if (iam_fim_task.or.iam_write_task) then

      call mpi_comm_rank (fimcomm_tweaked, iamnew, ignore)
!SMS$IGNORE BEGIN
      if (iamnew == 0) then
        write(6,'(a,a)') 'Compute task: WORLD rank 0 = fimcomm_tweaked rank 0 is running on ', &
          trim(mynode)
      else if (iamnew < begwriterank .or. iamnew > endwriterank) then
        write(6,'(a,i0,a,i0,a,a)') 'Compute task: WORLD rank ',world_rank, &
          ' = fimcomm_tweaked rank ', iamnew, ' is running on ', trim(mynode)
      else 
        write(6,'(a,i0,a,i0,a,a)') 'Write task: WORLD rank ',world_rank, &
          ' = fimcomm_tweaked rank ', iamnew, ' is running on ', trim(mynode)
      end if
!SMS$IGNORE END

    else   ! we are a "do-nothing" task

      my_comm = fimcomm_tweaked   ! used later in !SMS$SET_COMMUNICATOR
!SMS$IGNORE BEGIN
      write(6,'(a,i0,a,a)') 'Do-nothing task: WORLD rank ',world_rank, &
        ' is running on ', trim(mynode)
!SMS$IGNORE END
! NOTE:  The do-nothing tasks must return from this routine for NEMS runs 
! NOTE:  due to the ESMF requirement that all tasks in a VM participate 
! NOTE:  in component creation, even tasks that are not used.  

    end if

! TODO:  This is the opposite of DRY!  
! TODO:  The ordering of compute and write tasks is implicit in THREE 
! TODO:  PLACES.  DRY this out!!!  

! build task list(s) using world_rank ranks if argument(s) are present
! NOTE:  do-nothing tasks actually do do this, but calling them 
! NOTE:  "do-almost-nothing tasks" is not really all that illuminating
    if (present(tasklist_compute)) then
      numcomputetasks = lastproclist + 1 - nwt
      allocate(tasklist_compute(numcomputetasks))
      if (nwt > 0) then
        i2 = 1
        do i = 0,begwriterank-1
          tasklist_compute(i2) = proclist(i)
          i2 = i2 + 1
        enddo
        do i = endwriterank+1,lastproclist
          tasklist_compute(i2) = proclist(i)
          i2 = i2 + 1
        enddo
      else
        do i = 0,numcomputetasks-1
          tasklist_compute(i+1) = proclist(i)
        enddo
      endif
!SMS$IGNORE BEGIN
      if (debugon) then
        do i=0,num_tasks_world-1
          if (i == world_rank) then
            write (6,'(a,i0,a,i0)') 'DEBUG0c: ',world_rank,' numcomputetasks = ',numcomputetasks
            do n=1,numcomputetasks
              write (6,'(a,i0,a,i0,a,i0)') 'DEBUG0c: ',world_rank,' tasklist_compute(',n,') = ',tasklist_compute(n)
            enddo
          endif
          call flush(6)
          CALL MPI_BARRIER (tmp_comm, ierr)
          IF (ierr /= 0) THEN
            PRINT *,'ERROR in DEBUG0c:  MPI_BARRIER returned ',ierr
            call flush(6)
            CALL MPI_ABORT (MPI_COMM_WORLD, badret, ignore)
            STOP
          ENDIF
        enddo
      endif
!SMS$IGNORE END
    endif
    if (present(tasklist_write)) then
! allocate size=0 if nwt==0
      allocate(tasklist_write(nwt))
      if (nwt > 1) then
        i2 = 1
        do i = begwriterank,endwriterank
          tasklist_write(i2) = proclist(i)
          i2 = i2 + 1
        enddo
      endif
!SMS$IGNORE BEGIN
      if (debugon) then
        do i=0,num_tasks_world-1
          if (i == world_rank) then
            write (6,'(a,i0,a,i0)') 'DEBUG0w: ',world_rank,' nwt = ',nwt
            do n=1,nwt
              write (6,'(a,i0,a,i0,a,i0)') 'DEBUG0w: ',world_rank,' tasklist_write(',n,') = ',tasklist_write(n)
            enddo
          endif
          call flush(6)
          CALL MPI_BARRIER (tmp_comm, ierr)
          IF (ierr /= 0) THEN
            PRINT *,'ERROR in DEBUG0w:  MPI_BARRIER returned ',ierr
            call flush(6)
            CALL MPI_ABORT (MPI_COMM_WORLD, badret, ignore)
            STOP
          ENDIF
        enddo
      endif
!SMS$IGNORE END
    endif

    deallocate(proclist)
    deallocate(nodes)

! Split the remaining MPI communicator into "compute" and "write" 
! sections and create intercommunicators.  
! Not surprisingly, do-nothing tasks skip this.  
    finish_core_setup: if (iam_fim_task.or.iam_write_task) then

      CALL MPI_COMM_DUP (fimcomm_tweaked, total_comm, ierr)
      IF (ierr /= 0) THEN
!SMS$IGNORE BEGIN
        PRINT *,'ERROR in core_setup_fim:  MPI_COMM_DUP returned ',ierr
        CALL MPI_ABORT (MPI_COMM_WORLD, badret, ignore)
        STOP
!SMS$IGNORE END
      ENDIF

      CALL MPI_COMM_SIZE (total_comm, num_tasks, ierr)
      IF (ierr /= 0) THEN
!SMS$IGNORE BEGIN
        PRINT *,'ERROR in core_setup_fim:  MPI_COMM_SIZE returned ',ierr
        CALL MPI_ABORT (MPI_COMM_WORLD, badret, ignore)
        STOP
!SMS$IGNORE END
      ENDIF

      CALL MPI_COMM_RANK (total_comm, my_rank, ierr)
      IF (ierr /= 0) THEN
!SMS$IGNORE BEGIN
        PRINT *,'ERROR in core_setup_fim:  MPI_COMM_RANK returned ',ierr
        CALL MPI_ABORT (MPI_COMM_WORLD, badret, ignore)
        STOP
!SMS$IGNORE END
      ENDIF

! set up MPI communicators for FIM compute and optional write tasks
      IF (nwt == 0) THEN
        nct = num_tasks
! my_comm = total_comm
        CALL MPI_COMM_DUP (total_comm, my_comm, ierr)
        IF (ierr /= 0) THEN
!SMS$IGNORE BEGIN
          PRINT *,'ERROR in core_setup_fim:  MPI_COMM_DUP returned ',ierr
          CALL MPI_ABORT (MPI_COMM_WORLD, badret, ignore)
          STOP
!SMS$IGNORE END
        ENDIF
      ELSE   ! nwt > 0: nwt < 0 already checked above
        nct = num_tasks - nwt
        IF (nct <= 0) THEN
!SMS$IGNORE BEGIN
          PRINT *,'ERROR in core_setup_fim:  number of MPI tasks (',num_tasks, &
            ') <= num_write_tasks (',nwt,')'
          CALL MPI_ABORT (MPI_COMM_WORLD, badret, ignore)
          STOP
!SMS$IGNORE END
        ENDIF
! Split the MPI communicator and create intercommunicators for 
! FIM compute and write tasks.  
! See comments above for details.  
        IF (iam_fim_task) THEN
          color = 0
        ELSE
          color = 1
        ENDIF

        CALL MPI_COMM_SPLIT (total_comm, color, my_rank, my_comm, ierr)
        IF (ierr /= 0) THEN
!SMS$IGNORE BEGIN
          PRINT *,'ERROR in core_setup_fim:  MPI_COMM_SPLIT returned ',ierr
          CALL MPI_ABORT (MPI_COMM_WORLD, badret, ignore)
          STOP
!SMS$IGNORE END
        ENDIF
        IF (iam_fim_task) THEN
! call from FIM compute tasks: "remote leader" is begwriterank
          CALL MPI_INTERCOMM_CREATE (my_comm, 0, total_comm, begwriterank, intercomm_tag, &
            output_intercomm, ierr)
        ELSE IF (iam_write_task) THEN
! call from write tasks: "remote leader" is root (0)
!TODO:  extend to multiple groups of write tasks
          CALL MPI_INTERCOMM_CREATE (my_comm, 0, total_comm, 0,            intercomm_tag, &
            output_intercomm, ierr)
        ELSE
!SMS$IGNORE BEGIN
          PRINT *,'ERROR in core_setup_fim: entered no-mans land where I am neither ', &
            'fim task nor write task nor do-nothing task'
          CALL MPI_ABORT (MPI_COMM_WORLD, badret, ignore)
          STOP
!SMS$IGNORE END
        ENDIF

        IF (ierr /= 0) THEN
!SMS$IGNORE BEGIN
          PRINT *,'ERROR in core_setup_fim:  MPI_INTERCOMM_CREATE returned ',ierr
          CALL MPI_ABORT (MPI_COMM_WORLD, badret, ignore)
          STOP
!SMS$IGNORE END
        ENDIF
      ENDIF

    endif finish_core_setup

! Execution of directive on all tasks ensures that SMS_INIT is called 
! properly on the write tasks and on the do-nothing tasks.  This allows 
! these tasks to call NNT_STOP->NNT_EXIT without assertion errors.  As 
! a side-effect, duplicate messages are printed by the roots of the 
! write and do-nothing groups during NNT_EXIT().  
! Finally, tell SMS what communicator to use.  
! All tasks do this passing in their respective intra-communicators.  
! NOTE:  Executable SMS directives must not be placed before 
! NOTE:  this directive!  
! NOTE:  This includes writes or prints without !SMS$ignore because they 
! NOTE:  cause SMS to generate code.
!SMS$SET_COMMUNICATOR( my_comm )

!JR Print rank and node name. Useful for debugging assignment of tasks
!JR SMS$IGNORE keeps SMS from wrapping iam_root around the prints.

!SMS$IGNORE BEGIN
    if (iam_fim_task) then
      write(6,'(a,i0,a,a)')'FIM Compute task rank ', my_rank, ' is running on node ', trim(mynode)
    else if (iam_write_task) then
      write(6,'(a,i0,a,a)')'FIM Write task rank ', my_rank, ' is running on node ', trim(mynode)
    end if
!SMS$IGNORE END

    CALL IncrementTimer(t0,t1)
! print time only from FIM compute root
    IF (iam_fim_task) THEN
      PRINT *,'core_setup_fim time =',t1
      PRINT"(' Number of Write Tasks:'               ,I24,' processors')",nwt
    ENDIF

    core_setup_done = .true.

  end subroutine core_setup_fim

  logical function iam_compute_root()
    logical, save :: first_call = .true.
    logical, save :: save_result
    logical::initialized
    integer       :: my_rank, ierr, ignore
    if (first_call) then
      save_result = .false.
      call mpi_initialized(initialized,ierr)
      if (.not.initialized) then ! serial run
        save_result=.true.
      else
        if (iam_fim_task) then
! figure out who is the "root" of the compute tasks
          CALL MPI_COMM_RANK (my_comm, my_rank, ierr)
          IF (ierr /= 0) THEN
!SMS$IGNORE BEGIN
            PRINT *,'ERROR in iam_compute_root:  MPI_COMM_RANK returned ',ierr
            CALL MPI_ABORT (MPI_COMM_WORLD, badret, ignore)
            STOP
!SMS$IGNORE END
          ENDIF
          save_result = (my_rank == 0)
        endif
      endif
      first_call = .false.
    endif
    iam_compute_root = save_result
  end function iam_compute_root

!TODO:  remove duplication with iam_compute_root
  logical function iam_write_root()
    logical, save :: first_call = .true.
    logical, save :: save_result
    integer       :: my_rank, ierr, ignore
    if (first_call) then
      save_result = .false.
      if (.NOT.iam_fim_task) then
! figure out who is the "root" of the write tasks
        CALL MPI_COMM_RANK (my_comm, my_rank, ierr)
        IF (ierr /= 0) THEN
!SMS$IGNORE BEGIN
          PRINT *,'ERROR in iam_write_root:  MPI_COMM_RANK returned ',ierr
          CALL MPI_ABORT (MPI_COMM_WORLD, badret, ignore)
          STOP
!SMS$IGNORE END
        ENDIF
        save_result = (my_rank == 0)
      endif
      first_call = .false.
    endif
    iam_write_root = save_result
  end function iam_write_root

end module module_core_setup
