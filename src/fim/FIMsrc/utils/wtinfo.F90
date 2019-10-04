module module_wtinfo

  implicit none

contains

  subroutine wtinfo(cpn_out,num_write_tasks_out,mwtpn_out,root_own_node_out,&
    abort_on_bad_task_distrib_out,debugmsg_on_out,comm_in)

! Read namelist and return cores per node (cpn), number of write tasks
! (num_write_tasks), max write tasks per node (mwtpn), and whether root has node
! to himself (root_own_node).
!
! If this module is built in an SMS-parallelized context, insert and compile code
! to read the namelist only on the root task, package the write-task info in a
! user-defined type, and broadcast these values to the other tasks via an MPI
! derived type.

    implicit none

!sms$insert include 'mpif.h'

    type wtconf
      sequence
      integer::cpn
      integer::mwtpn
      integer::nwt
      logical::aobtd
      logical::do
      logical::ron
    end type wtconf

    integer, intent(out) :: cpn_out             ! number of cores per node
    integer, intent(out) :: mwtpn_out           ! max write tasks per node
    integer, intent(out) :: num_write_tasks_out ! number of write tasks
    integer, parameter   :: lun = 11            ! unit number
    logical, intent(out) :: abort_on_bad_task_distrib_out
    logical, intent(out) :: debugmsg_on_out     ! write-task debug msg control
    logical, intent(out) :: root_own_node_out   ! rank 0 has node to himself?
    integer,intent(in),optional :: comm_in      ! an MPI intracommunicator

    integer       :: comm,istatus,me=0,comm_wtconf,size_i,size_l
    integer       :: count(2),offset(2),type(2)
    type(wtconf)  :: conf

! Namelist variables

! If something fishy is encountered w.r.t. node names associated with MPI tasks,
! default is to abort the model

    integer, save :: cpn = 0 ! init to bad value (user MUST specify in namelist)
    integer, save :: max_write_tasks_per_node = 7
    integer, save :: num_write_tasks = 0    ! default is no write tasks
    logical, save :: abort_on_bad_task_distrib = .true.
    logical, save :: debugmsg_on=.false.    ! write-task debug message control
    logical, save :: root_own_node = .true. ! Put rank 0 on node by himself
    logical, save :: alreadyReadWriteTaskInfo = .false.

    namelist /WRITETASKnamelist/ abort_on_bad_task_distrib,cpn,debugmsg_on,&
      max_write_tasks_per_node,num_write_tasks,root_own_node

! Prefer the passed-in MPI communicator, if provided.

!sms$insert if (present(comm_in)) then
!sms$insert comm=comm_in
!sms$insert else
!sms$insert comm=mpi_comm_world
!sms$insert endif

    if (.not.alreadyReadWriteTaskInfo) then

!sms$insert call mpi_comm_rank(comm,me,istatus)

!sms$insert if (me.eq.0) then

!sms$ignore begin
      open (lun,file='FIMnamelist',status='old',action='read',iostat=istatus)

      if (istatus.ne.0) then
        write (*,'(a,i0,a,i0)') 'wtinfo: task ',me,&
          ' failed to open namelist file on unit ',lun
        call flush(6)
        stop
      endif

      read (lun,WRITETASKnamelist,iostat=istatus)
      if (istatus.ne.0) then
        write (*,'(a,i0,a,i0)') 'wtinfo: task ',me,' failed to read WRITETASKnamelist on unit ',lun
        call flush(6)
        stop
      endif

      close(lun)
!sms$ignore end

! Populate conf with individual values.

!sms$insert conf%aobtd = abort_on_bad_task_distrib
!sms$insert conf%cpn   = cpn
!sms$insert conf%do    = debugmsg_on
!sms$insert conf%mwtpn = max_write_tasks_per_node
!sms$insert conf%nwt   = num_write_tasks
!sms$insert conf%ron   = root_own_node

!sms$insert endif ! me.eq.0

! Determine the size in bytes of the MPI integer type

!sms$insert call mpi_type_extent(mpi_integer,size_i,istatus)
!sms$insert if (istatus.ne.0) then
!sms$ignore begin
!sms$insert write (*,'(a,i0)') 'wtinfo: MPI_Type_extent returned ',istatus
!sms$insert call mpi_abort(comm,istatus)
!sms$insert stop
!sms$ignore end
!sms$insert endif

! Determine the size in bytes of the MPI logical type

!sms$insert call mpi_type_extent(mpi_logical,size_l,istatus)
!sms$insert if (istatus.ne.0) then
!sms$ignore begin
!sms$insert write (*,'(a,i0)') 'wtinfo: MPI_Type_extent returned ',istatus
!sms$insert call mpi_abort(comm,istatus)
!sms$insert stop
!sms$ignore end
!sms$insert endif

! Populate the arrays defining the MPI derived type

!sms$insert count(1)=3
!sms$insert offset(1)=0
!sms$insert type(1)=mpi_integer
!sms$insert count(2)=3
!sms$insert offset(2)=offset(1)+(count(1)*size_i)
!sms$insert type(2)=mpi_logical

! Create the MPI derived-type struct

!sms$insert call mpi_type_struct(2,count,offset,type,comm_wtconf,istatus)
!sms$insert if (istatus.ne.0) then
!sms$ignore begin
!sms$insert write (*,'(a,i0)') 'wtinfo: MPI_Type_struct returned ',istatus
!sms$insert call mpi_abort(comm,istatus)
!sms$insert stop
!sms$ignore end
!sms$insert endif

! Commit the MPI derived type

!sms$insert call mpi_type_commit(comm_wtconf,istatus)
!sms$insert if (istatus.ne.0) then
!sms$ignore begin
!sms$insert write (*,'(a,i0)') 'wtinfo: MPI_Type_commit returned ',istatus
!sms$insert call mpi_abort(comm,istatus)
!sms$insert stop
!sms$ignore end
!sms$insert endif

! Broadcast the writetask configuration info

!sms$insert call mpi_bcast(conf,1,comm_wtconf,0,comm,istatus)
!sms$insert if (istatus.ne.0) then
!sms$ignore begin
!sms$insert write (*,'(a,i0)') 'wtinfo: MPI_Bcast returned ',istatus
!sms$insert call mpi_abort(comm,istatus)
!sms$insert stop
!sms$ignore end
!sms$insert endif

! Non-root tasks unpack conf to individual values.

!sms$insert if (me.ne.0) then
!sms$insert abort_on_bad_task_distrib = conf%aobtd
!sms$insert cpn                       = conf%cpn
!sms$insert debugmsg_on               = conf%do
!sms$insert max_write_tasks_per_node  = conf%mwtpn
!sms$insert num_write_tasks           = conf%nwt
!sms$insert root_own_node             = conf%ron
!sms$insert endif ! me.ne.0

      alreadyReadWriteTaskInfo=.true.

    endif ! .not.alreadyReadWriteTaskInfo

    num_write_tasks_out           = num_write_tasks
    mwtpn_out                     = max_write_tasks_per_node
    cpn_out                       = cpn
    root_own_node_out             = root_own_node
    abort_on_bad_task_distrib_out = abort_on_bad_task_distrib
    debugmsg_on_out               = debugmsg_on

  end subroutine wtinfo

end module module_wtinfo
