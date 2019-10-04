module icosio

! NOTE Two cpp tokens, SERIAL and NOGRIB, guard portions of this code. When
! NOTE SERIAL is set, code sections that introduce a dependence on MPI are
! NOTE removed by cpp. Therefore, if SERIAL is *not* set, icosio must be linked
! NOTE with MPI. Likewise, when NOGRIB is set, code sections relating to grib
! NOTE output file production are removed by cpp. For SERIAL, the following
! NOTE convention has been followed: If a routine is unneeded in serial mode but
! NOTE introduces no MPI dependence, it is left alone; if it is unneeded but
! NOTE does introduce an MPI dependence, it is completely removed; and if it
! NOTE is needed in serial mode, only the lines introducing MPI dependencies
! NOTE are removed. If NOGRIB is not set, the three routines post_init_file,
! NOTE post_write_field, and post_finalize_file must be available at link time.

! TODO Ideally, the necessary grib-production code would be moved into and built
! TODO as part of icosio.

  implicit none
  save
  private

#ifndef SERIAL

  include 'mpif.h'

! Public subroutines defining the icosio API: for parallel builds only

  public icosio_prep,icosio_run

#endif /* SERIAL */

! Public subroutines defining the icosio API: for parallel and serial builds

  public icosio_end_frame,icosio_out,icosio_set_inv_perm,icosio_setup

! Module parameters (for serial and parallel builds)
  
  integer,parameter::max_filename_len=80 ! max file name length
  integer,parameter::max_header_size=1024 ! max header length
  integer,parameter::max_output_files=100
  integer,parameter::max_varname_len=64 ! max variable name length

! Module variables (for serial and parallel builds)

  character(len=max_filename_len),allocatable::filename_list(:)
  character*1000::msg ! buffer for messages
  character*2::tasktype='__' ! set to 'ct' or 'wt' in write_share_init()
  integer,allocatable::interior_sizes(:)
  integer,pointer::inv_perm_global(:)=>null()
  integer,pointer::inv_perm_local(:)=>null()
  integer,private::size_c,size_i,size_l,size_r ! intrinsic sizes per MPI
  integer::comm_framecmd ! mpi derived type for frame commands
  integer::comm_varmeta ! mpi derived type for variable metadata
  integer::intercomm ! intercommunicator between compute and write tasks
  integer::interior_size=0 ! number of interior points on this compute task
  integer::intracomm ! communicator for my group
  integer::istatus
  integer::me ! 0..(n-1) MPI rank
  integer::nct=-1 ! number of compute tasks
  integer::nwt=-1 ! number of write tasks
  integer::outfiles=0 ! number of assigned write tasks
  integer::wtindex=0 ! tracks which write task is next to receive data
  logical::i_am_compute_root ! am I serial or root of a compute intracommunicator?
  logical::i_am_compute_task ! am I a compute task?
  logical::i_am_write_root ! am I root of a write intracommunicator?
  logical::icosio_setup_called=.false. ! has icosio_setup() been called?
  logical::serial ! am I running serial?
  logical::single ! am I serial or the only task in my intracommunicator?

! MPI comm tags

  integer,parameter::tag_cmd=100
  integer,parameter::tag_collect_inv_perm_control=101
  integer,parameter::tag_collect_inv_perm_segment=102
  integer,parameter::tag_collect_var_bounds=103
  integer,parameter::tag_collect_var_segment=104
  integer,parameter::tag_data=105
  integer,parameter::tag_interior_sizes=106
  integer,parameter::tag_inv_perm_control=107
  integer,parameter::tag_inv_perm_data=108
  integer,parameter::tag_metadata=109

! These shared variables, set in icosio_setup(), are set directly by the
! calling model code. See icosio_setup() for more details.

  character(len=12)::yyyymmddhhmm='____________' ! Date
  integer::comm                                  ! An MPI communicator
  integer::filename_len=-1                       ! Length of filenames
  integer::glvl=-1                               ! Grid level
  integer::header_size                           ! Product of header cols x rows
  integer::ime=-1                                ! Bounds: upper outer
  integer::ims=-1                                ! Bounds: lower outer
  integer::ipe=-1                                ! Bounds: upper inner
  integer::ips=-1                                ! Bounds: lower inner
  integer::lunout=-1                             ! Use this lun for disk writes
  integer::nip=-1                                ! Number of icosahedral points
  integer::nts=-1                                ! Number of time steps
  integer::nvl=-1                                ! Number of vertical levels
  integer::varname_len=-1                        ! Length of variable names
  logical::binout=.true.                         ! Write non-grib history files?
  logical::client_server_io                      ! Use client/server io?
  logical::debugmsg_on=.false.                   ! Print verbose status messages?
  logical::gribout=.false.                       ! Write grib files?
  logical::i_am_write_task                       ! Am I a write task?
  logical::print_diags                           ! Only print diagnostics?
  logical::using_write_tasks                     ! Are we using write tasks?
  real::dt=-1.0                                  ! Time step length

! Packet types for compute-task/write-task communications

  type framecmd
    sequence
    integer::its=-1      ! Time step
    integer::segments=-1 ! # of segments to expect from each compute task
  end type framecmd

  type varmeta
    sequence
    real::scalefactor=1.     ! Scale factor for grib output
    integer::accum_start=-1  ! Time accumulation factor for grib output
    integer::levels=-1       ! Number of vertical levels for this variable
    integer::segment_size=-1 ! Bytes in this variable segment
    integer::time=-1         ! Time this variable represents
    integer::filename_len=-1 ! Length of the supplied filename
    character(len=max_varname_len)::varname=''   ! Name of the variable
    character(len=max_filename_len)::filename='' ! Name of the output file
    character::header(max_header_size) ! Header for native binary output
! Padding may be required here if data isn't aligned to 8-byte boundary
  end type varmeta

! Types for buffering data & metadata

  type buffer
    integer::segments=0
    type(buffer_node),pointer::current=>null()
    type(buffer_node),pointer::head=>null()
  end type buffer

  type buffer_node
    real,pointer::segment(:)=>null()
    type(buffer_node),pointer::next=>null()
    type(varmeta)::vm
  end type buffer_node

  type(buffer),pointer::buffers(:) ! read/write buffers

contains

!--------------------------------------------------------------------------------
  subroutine append_to_list(writetask,varname,levels,segment_size,segment,&
    filename,header,time,scalefactor,accum_start)
!--------------------------------------------------------------------------------

! Attach a new node to the linked list of buffer nodes for the specified write
! task. If the optional variable name is not specified an empty node is appended;
! otherwise, it is assumed that the rest of the optional arguments are also
! present, and the variable metadata and field data are filled in with the
! supplied values. This is a private routine -- not part of the icosio API. It is
! called by both compute and write tasks.

    character(len=*),intent(in),optional::filename,varname
    character,intent(in),optional::header(header_size)
    integer,intent(in),optional::levels,segment_size,time,accum_start
    integer,intent(in)::writetask
    real,intent(in),optional::scalefactor

    character(len=14)::this='append_to_list'
    real,pointer,optional::segment(:)
    type(buffer),pointer::bp
    type(buffer_node),pointer::newnode

! Allocate a new node, set its variable metadata and point it at the passed-in
! variable data segment.

    allocate(newnode,stat=istatus)
    write(msg,'(a,a,a,a,a,i0,a)') 'ERROR: ',this,' (',tasktype,' ',me,&
      '): Failed to allocate new buffer_node.'
    call die_if(istatus.ne.0,istatus,msg)

! "segment" may be large, so use a pointer instead of copying. "segment" is
! deallocated in clear_list().

    if (present(varname)) then
      newnode%vm%levels=levels
      newnode%vm%segment_size=segment_size
      newnode%vm%varname(1:varname_len)=varname(1:varname_len)
      newnode%vm%filename(1:filename_len)=filename(1:filename_len)
      newnode%vm%header(1:header_size)=header(1:header_size)
      newnode%vm%time=time
      newnode%vm%scalefactor=scalefactor
      newnode%vm%accum_start=accum_start
      newnode%segment=>segment
      write(msg,'(a,a,a,a,i0,a,a,a)') this,' (',tasktype,' ',me,&
        '): Set varmeta for ',varname,'.'
      call debugmsg
    endif

! Append the new node to the appropriate linked list.

    bp=>buffers(writetask)
    if (associated(bp%head)) then
      bp%current%next=>newnode
    else
      bp%head=>newnode
    endif
    bp%current=>newnode
    bp%segments=bp%segments+1

    if (present(varname)) then
      write(msg,'(a,a,a,a,i0,a,a,a)') this,' (',tasktype,' ',me,&
        '): Appended ',varname,' to list.'
      call debugmsg
    else
      write(msg,'(a,a,a,a,i0,a)') this,' (',tasktype,' ',me,&
        '): Appended blank node to list.'
      call debugmsg
    endif

  end subroutine append_to_list

!--------------------------------------------------------------------------------
  subroutine buffer_var(its,varname,var,levels,filename,header,time,scalefactor,&
    accum_start)
!--------------------------------------------------------------------------------

! Buffer data for eventual transmission to write task(s), which is triggered by a
! icosio_flush() call. This is a private routine -- not part of the icosio API.
! It is called only by compute tasks.

    character(len=*),intent(in)::filename,varname
    character,intent(in)::header(header_size)
    integer,intent(in)::accum_start,its,levels,time
    real,intent(in)::scalefactor
    real,intent(in)::var(levels,ims:ime)

    character(len=10)::this='buffer_var'
    integer::i,ipn,ivl,offset,segment_size,writetask
    real,pointer::segment(:)
    type(buffer_node),pointer::node

    write(msg,'(a,a,a,a,i0,a,i0,a,a,a)') this,' (',tasktype,' ',me,'): its=',&
      its,' varname=',varname(1:varname_len),' entry'
    call debugmsg

! Lay out variable data in a 1D buffer.

    segment_size=levels*interior_size
    allocate(segment(segment_size),stat=istatus)
    write(msg,'(a,a,a,a,a,i0,a)') 'ERROR: ',this,' (',tasktype,' ',me,&
      '): Failed to allocate segment.'
    call die_if(istatus.ne.0,istatus,msg)
    offset=0
    do ipn=ips,ipe
      do ivl=1,levels
        offset=offset+1
        segment(offset)=var(ivl,ipn)
      enddo
    enddo

    write(msg,'(a,a,a,a,i0,a,a,a)') this,' (',tasktype,' ',me,&
      '): Laid out ',varname(1:varname_len),' in 1D buffer.'
    call debugmsg

! Determine which write task will receive this variable. mod() ensures
! round-robin assignment of responsibility for output files to write tasks.
! When there are fewer write tasks than output files, multiple variables will be
! assigned to some write task(s). When there are more write tasks than output
! files, write tasks not used in one output frame may be used in the next. In a
! future enhancement, it may be possible to split a single variable across two or
! more write tasks for output, for scalability.

! Search the buffers for a matching filename. If we find one, the same
! write task should receive this variable segment as well.

    writetask=0
    do i=1,nwt
      node=>buffers(i)%head
      do while (associated(node).and.writetask.eq.0)
        if (node%vm%filename(1:filename_len).eq.filename(1:filename_len))&
          writetask=i
        node=>node%next
      enddo
      if (writetask.ne.0) exit
    enddo

! If no matching filename was found, do round-robin assignment.

    if (writetask.eq.0) then
      wtindex=mod(wtindex,nwt)
      outfiles=outfiles+1
      wtindex=wtindex+1
      writetask=wtindex
      write(msg,'(a,a,a,a,i0,a,a,a,i0,a)') this,' (',tasktype,' ',me,'): ',&
        varname(1:varname_len),' assigned to wt ',writetask,' round-robin.'
      call debugmsg
    else
      write(msg,'(a,a,a,a,i0,a,a,a,i0,a)') this,' (',tasktype,' ',me,'): ',&
        varname(1:varname_len),' assigned to wt ',writetask,'.'
      call debugmsg
    endif

! Buffer data to send later.

    call append_to_list(writetask,varname,levels,segment_size,segment,filename,&
      header,time,scalefactor,accum_start)

    write(msg,'(a,a,a,a,i0,a,i0,a,a,a)') this,' (',tasktype,' ',me,'): its=',&
      its,' varname=',varname(1:varname_len),' exit'
    call debugmsg

  end subroutine buffer_var

!--------------------------------------------------------------------------------
  subroutine check_setup_called(caller)
!--------------------------------------------------------------------------------

! Return with error if icosio_setup() has not been called.

    character(len=*),intent(in)::caller

    write(msg,'(a,a,a,a,a,i0,a)') 'ERROR: ',caller,' (',tasktype,' ',me,&
      '): icosio_setup has not been called.'
    call die_if(.not.icosio_setup_called,istatus,msg)

  end subroutine check_setup_called

!--------------------------------------------------------------------------------
  subroutine clear_list
!--------------------------------------------------------------------------------

! Walk the list of data buffers for the specified write task and deallocate
! dynamic buffers. Reset the buffer-head pointers & counter. This is a private
! routine -- not part of the icosio API. It is called by both compute and write
! tasks.

    character(len=10)::this='clear_list'
    integer::i
    type(buffer_node),pointer::node,next

    write(msg,'(a,a,a,a,i0,a)') this,' (',tasktype,' ',me,'): entry'
    call debugmsg

    if (associated(buffers)) then
      do i=1,size(buffers)
        if (associated(buffers(i)%head)) then
          node=>buffers(i)%head
          do while (associated(node))
            next=>node%next
            deallocate(node%segment)
            node%segment=>null()
            deallocate(node)
            node=>next
          enddo
        endif
        buffers(i)%segments=0
        buffers(i)%current=>null()
        buffers(i)%head=>null()
      enddo
    else
      write(msg,'(a,a,a,a,a,i0,a)') 'ERROR: ',this,' (',tasktype,' ',me,&
        '): Called but buffers not allocated!'
      call die(istatus,msg)
    endif

    write(msg,'(a,a,a,a,i0,a)') this,' (',tasktype,' ',me,'): exit'
    call debugmsg

  end subroutine clear_list

!--------------------------------------------------------------------------------
  subroutine collect_inv_perm
!--------------------------------------------------------------------------------

! Collect the various segments of the inverse grid permutation array on the
! compute root. This is a private routine -- not parts of the icosio API. It is
! called only by compute tasks, and only when write tasks are not in use.

! TODO Consider using MPI 2 intercomm collective operations (probably MPI_Gather
! TODO and MPI_Gatherv) to get interior sizes and then global inv_perm.

    character(len=16)::this='collect_inv_perm'
    integer::bound(2),ct,inv_perm_global_size,segment_size

    write(msg,'(a,a,a,a,i0,a)') this,' (',tasktype,' ',me,'): entry'
    call debugmsg

! Do some initial error checking.

    write(msg,'(a,a,a,a,a,i0,a,a,a)') 'ERROR: ',this,' (',tasktype,' ',me,&
      '): Calling ',this,' twice is a bug.'
    call die_if(associated(inv_perm_global),istatus,msg)

    write(msg,'(a,a,a,a,a,i0,a,a,a)') 'ERROR: ',this,' (',tasktype,' ',me,&
      '): A write task calling ',this,' is a bug.'
    call die_if(associated(inv_perm_global),istatus,msg)

! The association status of inv_perm_global is used to signal the the compute
! root has already collected the necessary global-size inv_perm, so both root
! and non-root compute tasks need to allocate it. So, allocate the compute root's
! inv_perm_global as global size, and allocate the other compute tasks' as a
! single integer.

    if (i_am_compute_root) then
      inv_perm_global_size=nip
    else
      inv_perm_global_size=1
    endif

    allocate(inv_perm_global(inv_perm_global_size),stat=istatus)
    write(msg,'(a,a,a,a,a,i0,a)') 'ERROR: ',this,' (',tasktype,' ',me,&
      '): Failed to allocate inv_perm_global.'
    call die_if(istatus.ne.0,istatus,msg)

    if (i_am_compute_root) then

! The compute root already has its own local segment: Copy it into the correct
! location in inv_perm_global.

      inv_perm_global(ips:ipe)=inv_perm_local(ips:ipe)

! If running in single mode (serial or one-task MPI), inv_perm_global=inv_perm
! on the one and only task, and the previous assignment statement is sufficient,
! so simply return here.

      if (single) then
        write(msg,'(a,a,a,a,i0,a)') this,' (',tasktype,' ',me,&
          '): Single task returning...'
        call debugmsg
        return
      endif

#ifndef SERIAL
! Otherwise, the compute root collects the local inv_perm segments from each
! compute task to assemble a global-size inv_perm.

      do ct=1,nct-1 ! MPI tasks start from 0

        write(msg,'(a,a,a,a,i0,a,i0,a)') this,' (',tasktype,' ',me,&
          '): Receiving inv_perm bounds from ',ct,'.'
        call debugmsg

! Receive the start and end indicies of inv_perm_global into which to copy the
! segment about to be received.

        call mpi_recv(bound,2,mpi_integer,ct,tag_collect_inv_perm_control,&
          intracomm,mpi_status_ignore,istatus)
        write(msg,'(a,a,a,a,a,i0,a,i0)') 'ERROR: ',this,' (',tasktype,' ',me,&
          '): MPI_Recv returned ',istatus
        call die_if(istatus.ne.mpi_success,istatus,msg)

        segment_size=bound(2)-bound(1)+1

! Receive a segment into the appropriate location in the global array.

        write(msg,'(a,a,a,a,i0,a,i0,a)') this,' (',tasktype,' ',me,&
          '): Receiving inv_perm segment from ',ct,'.'
        call debugmsg

        call mpi_recv(inv_perm_global(bound(1):bound(2)),segment_size,&
          mpi_integer,ct,tag_collect_inv_perm_segment,intracomm,&
          mpi_status_ignore,istatus)
        write(msg,'(a,a,a,a,a,i0,a,i0)') 'ERROR: ',this,' (',tasktype,' ',me,&
          '): MPI_Recv returned ',istatus
        call die_if(istatus.ne.mpi_success,istatus,msg)

      enddo

    else ! I am a non-root compute task...

      inv_perm_global(1)=-1  ! A default value.

      segment_size=ipe-ips+1

      bound(1)=ips
      bound(2)=ipe

      write(msg,'(a,a,a,a,i0,a)') this,' (',tasktype,' ',me,&
        '): Sending inv_perm bounds to root.'
      call debugmsg

! Inform the compute root of the inv_perm_global bounds into which this task's
! segment should be written.

      call mpi_send(bound,2,mpi_integer,0,tag_collect_inv_perm_control,intracomm,&
        istatus)
      write(msg,'(a,a,a,a,a,i0,a,i0)') 'ERROR: ',this,' (',tasktype,' ',me,&
        '): MPI_Send returned ',istatus
      call die_if(istatus.ne.mpi_success,istatus,msg)

! Send this task's segment to the compute root.

      write(msg,'(a,a,a,a,i0,a)') this,' (',tasktype,' ',me,&
        '): Sending inv_perm segment to root.'
      call debugmsg

      call mpi_send(inv_perm_local(ips:ipe),segment_size,mpi_integer,0,&
        tag_collect_inv_perm_segment,intracomm,istatus)
      write(msg,'(a,a,a,a,a,i0,a,i0)') 'ERROR: ',this,' (',tasktype,' ',me,&
        '): MPI_Send returned ',istatus
      call die_if(istatus.ne.mpi_success,istatus,msg)
#endif /* SERIAL */

    endif ! (i_am_compute_root)

    write(msg,'(a,a,a,a,i0,a)') this,' (',tasktype,' ',me,'): exit'
    call debugmsg

  end subroutine collect_inv_perm

!--------------------------------------------------------------------------------
  subroutine collect_var(var,levels,glbvar)
!--------------------------------------------------------------------------------

! Collect the various segments of a global model field array on the compute root.
! This is a private routine -- not parts of the icosio API. It is called only by
! compute tasks, and only when write tasks are not in use.

! TODO Consider using MPI 2 intercomm collective operations (probably MPI_Gather
! TODO and MPI_Gatherv) to get interior sizes and then global fields.

    integer,intent(in)::levels
    real,intent(in)::var(levels,ims:ime)
    real,pointer::glbvar(:,:)

    character(len=11)::this='collect_var'
    integer::bound(2),ct,segment_size

    if (i_am_compute_root) then

! Allocate space for a global-size variable.

      allocate(glbvar(levels,nip),stat=istatus)
      write(msg,'(a,a,a,a,a,i0,a)') 'ERROR: ',this,' (',tasktype,' ',me,&
        '): Failed to allocate glbvar.'
      call die_if(istatus.ne.0,istatus,msg)

! Copy my local segment into the global array.

      glbvar(:,ips:ipe)=var(:,ips:ipe)

#ifndef SERIAL
! If nct=1 (i.e. we are 'single') this loop will not execute.

      do ct=1,nct-1 ! MPI tasks start from 0

! Receive the start and end indicies of inv_perm_global into which to copy the
! segment about to be received.

        write(msg,'(a,a,a,a,i0,a,i0,a)') this,' (',tasktype,' ',me,&
          '): Receiving var bounds from ',ct,'.'
        call debugmsg

        call mpi_recv(bound,2,mpi_integer,ct,tag_collect_var_bounds,intracomm,&
          mpi_status_ignore,istatus)
        write(msg,'(a,a,a,a,a,i0,a,i0)') 'ERROR: ',this,' (',tasktype,' ',me,&
          '): MPI_Recv returned ',istatus
        call die_if(istatus.ne.mpi_success,istatus,msg)

        segment_size=levels*(bound(2)-bound(1)+1)

! Receive a segment into the appropriate location in the global array.

        write(msg,'(a,a,a,a,i0,a,i0,a)') this,' (',tasktype,' ',me,&
          '): Receiving var segment from ',ct,'.'
        call debugmsg

        call mpi_recv(glbvar(:,bound(1):bound(2)),segment_size,mpi_real,ct,&
          tag_collect_var_segment,intracomm,mpi_status_ignore,istatus)
        write(msg,'(a,a,a,a,a,i0,a,i0)') 'ERROR: ',this,' (',tasktype,' ',me,&
          '): MPI_Recv returned ',istatus
        call die_if(istatus.ne.mpi_success,istatus,msg)

      enddo

    else ! I am a non-root compute task...

      segment_size=levels*(ipe-ips+1)

      bound(1)=ips
      bound(2)=ipe

! Inform the compute root of the inv_perm_global bounds into which this task's
! segment should be written.

      write(msg,'(a,a,a,a,i0,a)') this,' (',tasktype,' ',me,&
        '): Sending var bounds to root.'
      call debugmsg

      call mpi_send(bound,2,mpi_integer,0,tag_collect_var_bounds,intracomm,&
        istatus)
      write(msg,'(a,a,a,a,a,i0,a,i0)') 'ERROR: ',this,' (',tasktype,' ',me,&
        '): MPI_Send returned ',istatus
      call die_if(istatus.ne.mpi_success,istatus,msg)

! Send this task's segment to the compute root.

      write(msg,'(a,a,a,a,i0,a)') this,' (',tasktype,' ',me,&
        '): Sending var segment to root.'
      call debugmsg

      call mpi_send(var(:,ips:ipe),segment_size,mpi_real,0,&
        tag_collect_var_segment,intracomm,istatus)
      write(msg,'(a,a,a,a,a,i0,a,i0)') 'ERROR: ',this,' (',tasktype,' ',me,&
        '): MPI_Send returned ',istatus
      call die_if(istatus.ne.mpi_success,istatus,msg)
#endif /* SERIAL */

    endif ! (i_am_compute_root)

  end subroutine collect_var

#ifndef SERIAL
!--------------------------------------------------------------------------------
  subroutine comm_type_common(n,counts,offsets,types,comm_type)
!--------------------------------------------------------------------------------

! Common steps for creating MPI derived types. This is a private routine -- not
! part of the icosio API. It is called by both compute and write tasks.

    integer,intent(in)::n,counts(:),offsets(:),types(:)
    integer,intent(out)::comm_type

    call mpi_type_struct(n,counts,offsets,types,comm_type,istatus)
    call die_if(istatus.ne.mpi_success,istatus,'ERROR calling MPI_Type_struct.')
    call mpi_type_commit(comm_type,istatus)
    call die_if(istatus.ne.mpi_success,istatus,'ERROR calling MPI_Type_commit.')

  end subroutine comm_type_common
#endif /* SERIAL */

#ifndef SERIAL
!--------------------------------------------------------------------------------
  subroutine comm_type_framecmd(comm_framecmd)
!--------------------------------------------------------------------------------

! MPI derived type to communicate type framecmd (frame command). This is a
! private routine -- not part of the icosio API. It is called by both compute and
! write tasks.

    integer,intent(out)::comm_framecmd

    integer,parameter::n=1
    integer::counts(n),offsets(n),types(n)

    call comm_type_sizes

! For: its, segments
    counts(1)=2
    offsets(1)=0
    types(1)=mpi_integer

    call comm_type_common(n,counts,offsets,types,comm_framecmd)

  end subroutine comm_type_framecmd
#endif /* SERIAL*/

#ifndef SERIAL
!--------------------------------------------------------------------------------
  subroutine comm_type_setup
!--------------------------------------------------------------------------------

! Create the necessary MPI derived types. This is a private routine -- not part
! of the icosio API. It is called by both compute and write tasks.

    call comm_type_framecmd(comm_framecmd)
    call comm_type_varmeta(comm_varmeta)

  end subroutine comm_type_setup
#endif /* SERIAL */

#ifndef SERIAL
!--------------------------------------------------------------------------------
  subroutine comm_type_sizes
!--------------------------------------------------------------------------------

! Query MPI for the sizes of intrinsic types. This is a private routine -- not
! part of the icosio API. It is called by both compute and write tasks.

    logical::set=.false.

    if (.not.set) then
      call mpi_type_extent(mpi_character,size_c,istatus)
      call die_if(istatus.ne.mpi_success,istatus,&
        'ERROR calling MPI_Type_extent for MPI_CHARACTER.')
      call mpi_type_extent(mpi_integer,size_i,istatus)
      call die_if(istatus.ne.mpi_success,istatus,&
        'ERROR calling MPI_Type_extent for MPI_INTEGER.')
      call mpi_type_extent(mpi_logical,size_l,istatus)
      call die_if(istatus.ne.mpi_success,istatus,&
        'ERROR calling MPI_Type_extent for MPI_LOGICAL.')
      call mpi_type_extent(mpi_real,size_r,istatus)
      call die_if(istatus.ne.mpi_success,istatus,&
        'ERROR calling MPI_Type_extent for MPI_REAL.')
      set=.true.
    endif

  end subroutine comm_type_sizes
#endif /* SERIAL */

#ifndef SERIAL
!--------------------------------------------------------------------------------
  subroutine comm_type_varmeta(comm_varmeta)
!--------------------------------------------------------------------------------

! MPI derived type to communicate type varmeta (variable metadata). This is a
! private routine -- not part of the icosio API. It is called by both compute and
! write tasks.

    integer,intent(out)::comm_varmeta

    integer,parameter::n=3  ! number of MPI kinds
    integer::counts(n),offsets(n),types(n)

    call comm_type_sizes

! For: scalefactor
    counts(1)=1
    offsets(1)=0
    types(1)=mpi_real

! For: accum_start, levels, segment_size, time, filename_len
    counts(2)=5
    offsets(2)=offsets(1) + counts(1)*size_r
    types(2)=mpi_integer

! For: varname, filename, header
    counts(3)=max_varname_len+max_filename_len+max_header_size
    offsets(3)=offsets(2) + counts(2)*size_i
    types(3)=mpi_character

    call comm_type_common(n,counts,offsets,types,comm_varmeta)

  end subroutine comm_type_varmeta
#endif /* SERIAL */

!--------------------------------------------------------------------------------
  subroutine debugmsg
!--------------------------------------------------------------------------------

! Print a message to stdout if verbose debugging is enabled. This is a private
! routine -- not part of the icosio API. It is called by both compute and write
! tasks.

    if (debugmsg_on) then
      write (6,'(a,a)') 'DEBUGMSG: ',trim(msg)
      call flush(6)
    endif

  end subroutine debugmsg

!--------------------------------------------------------------------------------
  subroutine die(i,message,oldstatus)
!--------------------------------------------------------------------------------

! End program after printing the specified message to stdout. This is a private
! routine -- not part of the icosio API. It is called by both compute and write
! tasks.

    integer,intent(inout)::i
    character(len=*),intent(in)::message
    integer,intent(in),optional::oldstatus

    integer::ignore

    if (present(oldstatus)) then
      write (*,"(a,i0)") trim(message),oldstatus
    else
      write (*,*) trim(message)
    endif
    call flush(6)
#ifndef SERIAL
    call mpi_abort(mpi_comm_world,i,ignore)
#endif /* SERIAL */
    stop

  end subroutine die

!--------------------------------------------------------------------------------
  subroutine die_if(condition,i,message,oldstatus)
!--------------------------------------------------------------------------------

! End program and print message if condition is true. This is a private
! routine -- not part of the icosio API. It is called by both compute and write
! tasks.

    logical,intent(in)::condition
    integer,intent(inout)::i
    character(len=*),intent(in)::message
    integer,intent(in),optional::oldstatus

    if (condition) then
      if (present(oldstatus)) then
        call die(i,message,oldstatus)
      else
        call die(i,message)
      endif
    endif

  end subroutine die_if

!--------------------------------------------------------------------------------
  subroutine icosio_end_frame(its)
!--------------------------------------------------------------------------------

! Reset the list of filenames written to during this output interval. If write
! tasks are not in use, simply return; otherwise, call icosio_flush(). This is a
! public routine -- part of the icosio API. It is called by both compute and
! write tasks.

    integer,intent(in)::its

    character(len=16)::this='icosio_end_frame'

    call check_setup_called(this)
    if (allocated(filename_list)) deallocate(filename_list)
    if (nwt.eq.0) return
#ifndef SERIAL
    call icosio_flush(its)
#endif /* SERIAL */

  end subroutine icosio_end_frame

#ifndef SERIAL
!--------------------------------------------------------------------------------
  subroutine flush_one(its,writetask)
!--------------------------------------------------------------------------------

! Send the metadata and field-data segments accumulated in one write-task linked-
! list buffer to the appropriate write task. This is a private routine -- not
! part of the icosio API. It is called only by compute tasks.

    integer,intent(in)::its,writetask

    character(len=9)::this='flush_one'
    integer::wt
    type(buffer_node),pointer::head,node
    type(framecmd)::cmd

    write(msg,'(a,a,a,a,i0,a,i0,a,i0,a)') this,' (',tasktype,' ',me,&
      '): its=',its,' writetask=',writetask,' entry'
    call debugmsg

    head=>buffers(writetask)%head

    if (associated(head)) then

      wt=writetask-1 ! MPI tasks start from 0

      if ((writetask.lt.1).or.(writetask.gt.nwt)) then
        write(msg,'(a,a,a,a,a,i0,a,i0,a)') 'ERROR: ',this,' (',tasktype,' ',me,&
          '): (writetask=',writetask,') Out of range.'
        call die(istatus,msg)
      endif

! Send framecmd. A frame command with value > 0 informs the write task that a
! frame of output metadata and field data follows, and its value indicates how
! many different variables' segments will be sent.

      if (i_am_compute_root) then
        cmd%its=its
        cmd%segments=buffers(writetask)%segments
        call mpi_send(cmd,1,comm_framecmd,wt,tag_cmd,intercomm,istatus)
        write(msg,'(a,a,a,a,a,i0,a,i0,a,i0,a)') 'ERROR: ',this,' (',tasktype,&
          ' ',me,'): MPI_Send of go cmd to wt ',writetask,' returned ',istatus,&
          '.'
        call die_if(istatus.ne.mpi_success,istatus,msg)
      endif

! Send metadata: Walk the linked list and send the metadata packets.

      node=>head
      do while (associated(node))
        write(msg,'(a,a,a,a,i0,a,i0,a,a,a)') this,' (',tasktype,' ',me,&
          '): sending wt ',writetask,' ',trim(node%vm%varname),' metadata...'
        call debugmsg
        call mpi_send(node%vm,1,comm_varmeta,wt,tag_metadata,intercomm,istatus)
        write(msg,'(a,a,a,a,a,i0,a,i0,a,i0,a)') 'ERROR: ',this,' (',tasktype,&
          ' ',me,'): MPI_Send var meta to ',writetask,' returned ',istatus,'.'
        call die_if(istatus.ne.mpi_success,istatus,msg)
        node=>node%next
      enddo

! Send variable data: Walk the linked list and send the field data packets.

      node=>head
      do while (associated(node))
        call mpi_send(node%segment,node%vm%segment_size,mpi_real,wt,tag_data,&
          intercomm,istatus)
        write(msg,'(a,a,a,a,a,i0,a,i0,a,i0,a)') 'ERROR: ',this,' (',tasktype,&
          ' ',me,'): MPI_Send var data to ',writetask,' returned ',istatus,'.'
        call die_if(istatus.ne.mpi_success,istatus,msg)
        node=>node%next
      enddo

    endif ! associated(head)

    write(msg,'(a,a,a,a,i0,a,i0,a,i0,a)') this,' (',tasktype,' ',me,'): its=',&
      its,' writetask=',writetask,' exit'
    call debugmsg

  end subroutine flush_one
#endif /* SERIAL */

#ifndef SERIAL
!--------------------------------------------------------------------------------
  subroutine icosio_flush(its)
!--------------------------------------------------------------------------------

! Send accumulated metadata and field-data segments to all write tasks involed in
! the output of data for this frame. This is a private routine -- not part of
! the icosio API. If is called only by compute tasks.

    integer,intent(in)::its

    character(len=12)::this='icosio_flush'
    integer::i,wt,j,k
    logical::warning_given=.false.
    type(framecmd)::cmd

! Return if buffers have not accumulated any output data

    if (buffers(1)%segments.eq.0) return

! Print one-per-file write task layout instructions.

    if (i_am_compute_root.and..not.warning_given) then
      if (nwt.eq.1) then
        write(*,'(a,i0,a)') 'NOTE: Using one write task for ',outfiles,&
          ' output files.'
      else if (outfiles.ne.nwt) then
        write(*,'(a,i0,a)') 'NOTE: Set num_write_tasks=',outfiles,&
          ' in namelist file for a one-per-file write-task layout.'
      endif
      warning_given=.true.
    endif

! Flush the buffers of all write tasks holding data from the current frame.

    i=mod(outfiles-wtindex,nwt)
    if (outfiles.lt.nwt) then
      j=outfiles-1
    else
      j=nwt-1
    endif
    do k=i,i+j
      call flush_one(its,mod(k,nwt)+1)
    enddo

    call clear_list

! If this was the terminal timestep, tell the write tasks we're done by sending a
! framecmd packet with a zero value...

    if (its.eq.nts) then
      if (i_am_compute_root) then
        cmd%segments=0.
        do i=1,nwt
          wt=i-1
          call mpi_send(cmd,1,comm_framecmd,wt,tag_cmd,intercomm,istatus)
          write(msg,'(a,a,a,a,a,i0,a,i0,a,i0,a)') 'ERROR: ',this,' (',tasktype,&
            ' ',me,'): MPI_Send terminal cmd to wt ',i,' returned ',istatus,'.'
          call die_if(istatus.ne.mpi_success,istatus,msg)
        enddo
      endif
    else ! ...otherwise reset for next frame.
      outfiles=0
    endif

  end subroutine icosio_flush
#endif /* SERIAL */

#ifndef SERIAL
!--------------------------------------------------------------------------------
  subroutine icosio_prep
!--------------------------------------------------------------------------------

! Send some one-time initialization data to write tasks. Allocates some arrays
! needed for communication with write tasks. This is a public routine -- part of
! the icosio API. It is called only by compute tasks.

! TODO Consider using MPI 2 intercomm collective operations here to broadcast
! TODO interior size and inv_perm segments.

    character(len=11)::this='isocio_prep'
    integer::wt
    logical::send_inv_perm

    write(msg,'(a,a,a,a,i0,a,i0,a)') this,' (',tasktype,' ',me,'): entry'
    call debugmsg

    call check_setup_called(this)

! Set up MPI derived types.

    call comm_type_setup

! Allocate dynamic buffers.

    allocate(buffers(nwt),stat=istatus)
    write(msg,'(a,a,a,a,a,i0,a)') 'ERROR: ',this,' (',tasktype,' ',me,&
      '): Failed to allocate write buffers.'
    call die_if(istatus.ne.0,istatus,msg)

! Calculate interior points owned by this compute task.

    interior_size=ipe-ips+1

! Send interior size to write task(s).

    write(msg,'(a,a,a,a,i0,a)') this,' (',tasktype,' ',me,&
      '): Sending interior_size to write tasks.'
    call debugmsg
    do wt=0,nwt-1 ! MPI tasks start from 0
      call mpi_send(interior_size,1,mpi_integer,wt,tag_interior_sizes,&
        intercomm,istatus)
      write(msg,'(a,a,a,a,a,i0,a,i0,a)') 'ERROR: ',this,' (',tasktype,' ',me,&
        '): MPI_Send interior_size returned ',istatus,'.'
      call die_if(istatus.ne.mpi_success,istatus,msg)
      write(msg,'(a,a,a,a,i0,a,i0,a,i0,a)') this,' (',tasktype,' ',me,&
        '): Sent interior_size ',interior_size,' to wt ',wt+1,'.'
      call debugmsg
    enddo

! If inv_perm_local is associated -- by a call to icosio_set_inv_perm() -- then
! it must be sent to the write task(s). Decide, then let the write task(s) know
! whether or not inv_perm segments will be sent.

    send_inv_perm=.false.
    if (associated(inv_perm_local)) then
      send_inv_perm=.true.
    endif

    write(msg,'(a,a,a,a,i0,a)') this,' (',tasktype,' ',me,&
      '): Sending inv_perm control bit to write tasks.'
    call debugmsg

    do wt=0,nwt-1 ! MPI tasks start from 0
      call mpi_send(send_inv_perm,1,mpi_logical,wt,tag_inv_perm_control,&
        intercomm,istatus)
      write(msg,'(a,a,a,a,a,i0,a,i0,a)') 'ERROR: ',this,' (',tasktype,' ',me,&
        '): MPI_Send send_inv_perm returned ',istatus,'.'
      call die_if(istatus.ne.mpi_success,istatus,msg)
      write(msg,'(a,a,a,a,i0,a,i0,a)') this,' (',tasktype,' ',me,&
        '): Sent inv_perm control bit to wt ',wt+1,'.'
      call debugmsg
    enddo

! Send local inv_perm segment to write task(s), if needed.

    if (send_inv_perm) then

      write(msg,'(a,a,a,a,i0,a)') this,' (',tasktype,' ',me,&
        '): Sending inv_perm segment to write tasks.'
      call debugmsg

      do wt=0,nwt-1 ! MPI tasks start from 0
        call mpi_send(inv_perm_local(ips:ipe),interior_size,mpi_integer,wt,&
          tag_inv_perm_data,intercomm,istatus)
        write(msg,'(a,a,a,a,a,i0,a,i0,a)') 'ERROR: ',this,' (',tasktype,' ',me,&
          '): MPI_Send inv_perm returned ',istatus,'.'
        call die_if(istatus.ne.mpi_success,istatus,msg)
        write(msg,'(a,a,a,a,i0,a,i0,a)') this,' (',tasktype,' ',me,&
          '): Sent inv_perm segment to wt ',wt+1,'.'
        call debugmsg
      enddo

    endif

    write(msg,'(a,a,a,a,i0,a)') this,' (',tasktype,' ',me,'): exit'
    call debugmsg

  end subroutine icosio_prep
#endif /* SERIAL */

!--------------------------------------------------------------------------------
  subroutine icosio_out(its,time,varname,var,levels,filename,header,scalefactor,&
    accum_start)
!--------------------------------------------------------------------------------

! Receive a field-data array and either write it to disk directly, or buffer it
! for a later transmisstion to write task(s). This is a public routine -- part of
! the icosio API. It is called only by compute tasks.

    character(len=*),intent(in)::filename,varname
    character,intent(in)::header(header_size)
    integer,intent(in),optional::accum_start
    integer,intent(in)::levels,its,time
    real,intent(in),optional::scalefactor
    real,intent(in)::var(levels,ims:ime)

    character(len=10)::this='icosio_out'
    integer::accum_start_local
    real,pointer::glbvar(:,:)
    real::scalefactor_local

    call check_setup_called(this)

! Set default values for scalefactor (factor to scale input by for grib output)
! and accum_start (accumulation start) if not provided by caller.

    scalefactor_local=1.
    if (present(scalefactor)) scalefactor_local=scalefactor

    accum_start_local=-1
    if (present(accum_start)) accum_start_local=accum_start

! If write tasks are in use, buffer data for eventual overlapped write;
! otherwise, write data directly.

    if (using_write_tasks) then

      call buffer_var(its,varname,var,levels,filename,header,time,&
        scalefactor_local,accum_start_local)

    else

! Collect inv_perm on the compute root if fixed grid-order reordering is needed.

      if (i_am_compute_task) then
        if (associated(inv_perm_local)) then
          if (.not.associated(inv_perm_global)) then
            call collect_inv_perm
          endif
        endif
      endif

! If running in single mode (serial or one-task MPI), write the array immediately
! to disk. Otherwise, collect the distributed array onto the compute root and
! then write to disk.

      if (single) then
        call var_to_disk(accum_start_local,ips,ims,ipe,ime,filename,header,its,&
          levels,scalefactor_local,time,var,varname)
      else
        call collect_var(var,levels,glbvar)
        if (i_am_compute_root) then
          call var_to_disk(accum_start_local,1,1,nip,nip,filename,header,its,&
            levels,scalefactor_local,time,glbvar,varname)
          deallocate(glbvar)
        endif
      endif

    endif

  end subroutine icosio_out

#ifndef SERIAL
!--------------------------------------------------------------------------------
  subroutine icosio_run
!--------------------------------------------------------------------------------

! Main routine for write tasks. Allocate dynamic buffers, receive some one-time
! data, then enter a loop waiting for and responding to frame commands, which
! trigger either the receipt and output of a frame of output data, or stop. This
! is a public routine -- part of the icosio API. It is called only by write
! tasks.

! TODO Consider using MPI 2 intercomm collective operations here to broadcast
! TODO (receive, in ths case) interior size and inv_perm segments.

    character(len=10)::this='icosio_run'
    integer::ct,offset
    logical::recv_inv_perm,running=.true.
    type(framecmd)::cmd

    call check_setup_called(this)

! If the model is running in an "only print diagnostics" mode, either stop or
! return, depending on the client/server io settting.

    if (print_diags) then
      if (client_server_io) then
        call mpi_finalize(istatus)
        stop
      else
        return
      endif
    endif

! Set up MPI derived types.

    call comm_type_setup

! Allocate dynamic buffers.

    allocate(buffers(nct),stat=istatus)
    write(msg,'(a,a,a,a,a,i0,a)') 'ERROR: ',this,' (',tasktype,' ',me,&
      '): Failed to allocate buffers.'
    call die_if(istatus.ne.0,istatus,msg)

! Allocate and receive interior sizes from compute task(s).

    write(msg,'(a,a,a,a,i0,a)') this,' (',tasktype,' ',me,&
      '): Receiving interior sizes from compute tasks.'
    call debugmsg

    allocate(interior_sizes(0:nct-1),stat=istatus)
    write(msg,'(a,a,a,a,a,i0,a)') 'ERROR: ',this,' (',tasktype,' ',me,&
      '): Failed to allocate interior_sizes.'
    call die_if(istatus.ne.0,istatus,msg)

    do ct=0,nct-1
      call mpi_recv(interior_sizes(ct),1,mpi_integer,ct,tag_interior_sizes,&
        intercomm,mpi_status_ignore,istatus)
      write(msg,'(a,a,a,a,a,i0,a,i0,a)') 'ERROR: ',this,' (',tasktype,' ',me,&
        '): MPI_Recv interior_sizes returned ',istatus,'.'
      call die_if(istatus.ne.mpi_success,istatus,msg)
      write(msg,'(a,a,a,a,i0,a,i0,a,i0,a)') this,' (',tasktype,' ',me,&
        '): Received interior size ',interior_sizes(ct),' from ct ',ct,'.'
      call debugmsg
    enddo

! Allocate and receive from compute task(s) inv_perm, if needed.

    call mpi_recv(recv_inv_perm,1,mpi_logical,0,tag_inv_perm_control,&
      intercomm,mpi_status_ignore,istatus)
    write(msg,'(a,a,a,a,a,i0,a,i0,a)') 'ERROR: ',this,' (',tasktype,' ',me,&
      '): MPI_Recv recv_inv_perm returned ',istatus,'.'
    call die_if(istatus.ne.mpi_success,istatus,msg)
    write(msg,'(a,a,a,a,i0,a,i0,a)') this,' (',tasktype,' ',me,&
      '): Received inv_perm control bit from ct ',ct,'.'
    call debugmsg

    if (recv_inv_perm) then

      allocate(inv_perm_global(nip),stat=istatus)
      write(msg,'(a,a,a,a,a,i0,a)') 'ERROR: ',this,' (',tasktype,' ',me,&
        '): Failed to allocate inv_perm_global.'
      call die_if(istatus.ne.0,istatus,msg)

      offset=1
      do ct=0,nct-1
        call mpi_recv(inv_perm_global(offset:offset-1+interior_sizes(ct)),&
          interior_sizes(ct),mpi_integer,ct,tag_inv_perm_data,intercomm,&
          mpi_status_ignore,istatus)
        write(msg,'(a,a,a,a,a,i0,a,i0,a)') 'ERROR: ',this,' (',tasktype,' ',me,&
          '): MPI_Recv inv_perm_global returned ',istatus,'.'
        call die_if(istatus.ne.mpi_success,istatus,msg)
        write(msg,'(a,a,a,a,i0,a,i0,a)') this,' (',tasktype,' ',me,&
          '): Received inv_perm segment from ct ',ct,'.'
        call debugmsg
        offset=offset+interior_sizes(ct)
      enddo

    endif ! (recv_inv_perm)

! Main loop

    do while (running)

! Receve a frame command from the compute root.

      call mpi_recv(cmd,1,comm_framecmd,0,tag_cmd,intercomm,&
        mpi_status_ignore,istatus)
      write(msg,'(a,a,a,a,a,i0,a,i0,a)') 'ERROR: ',this,' (',tasktype,' ',me,&
        '): MPI_Recv cmd returned ',istatus,'.'
      call die_if(istatus.ne.mpi_success,istatus,msg)
      call die_if(cmd%segments.lt.0,istatus,'framecmd segment count < 0.')

! If a zero frame command value was sent, stop in the appriate manner...

      if (cmd%segments.eq.0) then
        write(msg,'(a,a,a,a,i0,a)') this,' (',tasktype,' ',me,&
          '): Stop command received.'
        call debugmsg
        running=.false.
        if (client_server_io) then
          write(msg,'(a,a,a,a,i0,a)') this,' (',tasktype,' ',me,&
            '): Client-server mode enabled, stopping...'
          call debugmsg
          call mpi_finalize(istatus)
          stop
        else
          return
        endif
      else ! ...otherwise, process the incoming frame of output data.
        write(msg,'(a,a,a,a,i0,a,i0,a)') this,' (',tasktype,' ',me,&
          '): its=',cmd%its,' Go command received.'
        call debugmsg
        call process_frame(cmd%its,cmd%segments)
      endif
    enddo

  end subroutine icosio_run
#endif /* SERIAL */

!--------------------------------------------------------------------------------
  subroutine icosio_set_inv_perm(inv_perm_in)
!--------------------------------------------------------------------------------

! Point inv_perm_local to the supplied inv_perm segment. If this routine is
! called, icosio assumes that grid reordering (via a reassembled global-size
! inverse grid permutation array) should be applied to variables before they are
! written to disk. Otherwise, arrays will be written as-is. This is a public
! routine -- part of the icosio inteface. It is called only by compute tasks.

    character(len=19)::this='icosio_set_inv_perm'
    integer,allocatable,target,intent(in)::inv_perm_in(:)

    call check_setup_called(this)

    inv_perm_local=>inv_perm_in

  end subroutine icosio_set_inv_perm

!--------------------------------------------------------------------------------
  subroutine icosio_setup(binout_in,client_server_io_in,comm_in,debugmsg_on_in,&
    dt_in,filename_len_in,glvl_in,gribout_in,header_size_in,i_am_write_task_in,&
    ips_in,ims_in,ipe_in,ime_in,lunout_in,nip_in,nts_in,nvl_in,print_diags_in,&
    using_write_tasks_in,varname_len_in,yyyymmddhhmm_in)
!--------------------------------------------------------------------------------

! Sets values known by the calling model and required by icosio. This is a public
! routine -- part of the icosio API. It is called by both compute and write
! tasks.

    character(len=12),intent(in)::yyyymmddhhmm_in
    integer,intent(in)::comm_in,filename_len_in,glvl_in,header_size_in,ims_in,&
      ips_in,ipe_in,ime_in,lunout_in,nip_in,nts_in,nvl_in,varname_len_in
    logical,intent(in)::binout_in,client_server_io_in,debugmsg_on_in,gribout_in,&
      i_am_write_task_in,print_diags_in,using_write_tasks_in
    real,intent(in)::dt_in

    character(len=12)::this='icosio_setup'
    integer::group
    logical::have_intercomm

! Set module variables.
!
! A number of necessary parameters are deduced from three settings: comm,
! i_am_write_task, and using_write_tasks. If write tasks are in use, comm must
! be an intercommunicator between the compute and write tasks; otherwise, it must
! be an intracommunicator for the compute tasks. The two logical values have the
! obvious meanings.

    binout=binout_in
    client_server_io=client_server_io_in
    comm=comm_in
    debugmsg_on=debugmsg_on_in
    dt=dt_in
    filename_len=filename_len_in
    glvl=glvl_in
    gribout=gribout_in
    header_size=header_size_in
    i_am_write_task=i_am_write_task_in
    ime=ime_in
    ims=ims_in
    ipe=ipe_in
    ips=ips_in
    lunout=lunout_in
    nip=nip_in
    nts=nts_in
    nvl=nvl_in
    print_diags=print_diags_in
    using_write_tasks=using_write_tasks_in
    varname_len=varname_len_in
    yyyymmddhhmm=yyyymmddhhmm_in

! Set serial non-MPI defaults.

    i_am_compute_root=.true.
    i_am_compute_task=.true.
    i_am_write_root=.false.
    me=0
    nct=1
    nwt=0
    serial=.true.
    single=.true.
    tasktype='ct'

#ifndef SERIAL
    intercomm=mpi_comm_null ! default value
    intracomm=mpi_comm_null ! default value

! Override defaults.

    if (using_write_tasks) then
      call mpi_comm_test_inter(comm,have_intercomm,istatus)
      write(msg,'(a,a,a,a,a,i0,a,i0)') 'ERROR: ',this,' (',tasktype,' ',me,&
        '): MPI_Comm_test_inter returned ',istatus
      call die_if(istatus.ne.mpi_success,istatus,msg)
      write(msg,'(a,a,a,a,a,i0,a)') 'ERROR: ',this,' (',tasktype,' ',me,&
        '): Pass icosio an intercommunicator when write tasks are enabled.'
      call die_if(.not.have_intercomm,istatus,msg)
      call mpi_comm_group(comm,group,istatus)
      write(msg,'(a,a,a,a,a,i0,a,i0)') 'ERROR: ',this,' (',tasktype,' ',me,&
        '): MPI_Comm_group returned ',istatus
      call die_if(istatus.ne.mpi_success,istatus,msg)
      call mpi_comm_create(comm,group,intracomm,istatus)
      write(msg,'(a,a,a,a,a,i0,a,i0)') 'ERROR: ',this,' (',tasktype,' ',me,&
        '): MPI_Comm_create returned ',istatus
      call die_if(istatus.ne.mpi_success,istatus,msg)
      intercomm=comm
    else if (comm.ne.mpi_comm_null) then
      call mpi_comm_test_inter(comm,have_intercomm,istatus)
      write(msg,'(a,a,a,a,a,i0,a,i0)') 'ERROR: ',this,' (',tasktype,' ',me,&
        '): MPI_Comm_test_inter returned ',istatus
      call die_if(istatus.ne.mpi_success,istatus,msg)
      write(msg,'(a,a,a,a,a,i0,a)') 'ERROR: ',this,' (',tasktype,' ',me,&
        '): Pass icosio an intracommunicator when write tasks are disabled.'
      call die_if(have_intercomm,istatus,msg)
      intracomm=comm
    endif

    if (i_am_write_task) then
      call mpi_comm_rank(intracomm,me,istatus)
      write(msg,'(a,a,a,a,a,i0,a,i0)') 'ERROR: ',this,' (',tasktype,' ',me,&
        '): MPI_Comm_rank returned ',istatus
      call die_if(istatus.ne.mpi_success,istatus,msg)
      call mpi_comm_remote_size(intercomm,nct,istatus)
      write(msg,'(a,a,a,a,a,i0,a,i0)') 'ERROR: ',this,' (',tasktype,' ',me,&
        '): MPI_Comm_remote_size returned ',istatus
      call die_if(istatus.ne.mpi_success,istatus,msg)
      call mpi_comm_size(intracomm,nwt,istatus)
      write(msg,'(a,a,a,a,a,i0,a,i0)') 'ERROR: ',this,' (',tasktype,' ',me,&
        '): MPI_Comm_size returned ',istatus
      call die_if(istatus.ne.mpi_success,istatus,msg)
      i_am_compute_root=.false.
      i_am_compute_task=.false.
      serial=.false.
      if (nwt.gt.1) single=.false.
      tasktype='wt'
      if (me.eq.0) i_am_write_root=.true.
    else
      if (intracomm.ne.mpi_comm_null) then
        call mpi_comm_rank(intracomm,me,istatus)
        write(msg,'(a,a,a,a,a,i0,a,i0)') 'ERROR: ',this,' (',tasktype,' ',me,&
          '): MPI_Comm_rank returned ',istatus
        call die_if(istatus.ne.mpi_success,istatus,msg)
        call mpi_comm_size(intracomm,nct,istatus)
        write(msg,'(a,a,a,a,a,i0,a,i0)') 'ERROR: ',this,' (',tasktype,' ',me,&
          '): MPI_Comm_size returned ',istatus
        call die_if(istatus.ne.mpi_success,istatus,msg)
        if (me.ne.mpi_success) i_am_compute_root=.false.
        serial=.false.
        if (nct.ne.1) single=.false.
      endif
      if (intercomm.ne.mpi_comm_null) then
        call mpi_comm_remote_size(intercomm,nwt,istatus)
        write(msg,'(a,a,a,a,a,i0,a,i0)') 'ERROR: ',this,' (',tasktype,' ',me,&
          '): MPI_Comm_remote_size returned ',istatus
        call die_if(istatus.ne.mpi_success,istatus,msg)
      endif
    endif
#endif /* SERIAL */

! Check for problems.

    write(msg,'(a,a,a,a,a,i0,a,i0,a,i0,a)') 'ERROR: ',this,' (',tasktype,' ',me,&
      '): model filename length (',filename_len,') exceeds icosio max (',&
      max_filename_len,'). Adjust "max_filename_len" parameter in icosio.'
    call die_if(filename_len.gt.max_filename_len,istatus,msg)

    write(msg,'(a,a,a,a,a,i0,a,i0,a,i0,a)') 'ERROR: ',this,' (',tasktype,' ',me,&
      '): model varname length (',varname_len,') exceeds icosio max (',&
      max_varname_len,'). Adjust "max_varname_len" parameter in icosio.'
    call die_if(varname_len.gt.max_varname_len,istatus,msg)

    write(msg,'(a,a,a,a,a,i0,a,i0,a,i0,a)') 'ERROR: ',this,' (',tasktype,' ',me,&
      '): model header size (',header_size,') exceeds icosio max (',&
      max_header_size,'). Adjust "max_header_size" parameter in icosio.'
    call die_if(varname_len.gt.max_varname_len,istatus,msg)

! Record that icosio_setup has run.

    icosio_setup_called=.true.

#ifndef SERIAL
! If I am a write task, enter client-server mode if selected.

    if (i_am_write_task.and.client_server_io) call icosio_run
#endif

  end subroutine icosio_setup

#ifndef SERIAL
!--------------------------------------------------------------------------------
  subroutine process_frame(its,segments)
!--------------------------------------------------------------------------------

! Receive metadata and field data from compute tasks, then write to disk. This is
! a public routine -- part of the icosio API. It is called only by write tasks.

    integer,intent(in)::its,segments

    character(len=13)::this='process_frame'
    integer::ct,ctbuffer,isegment
    type(buffer_node),pointer::node

    write(msg,'(a,a,a,a,i0,a,i0,a,i0,a)') this,' (',tasktype,' ',me,'): its=',&
      its,' segments=',segments,' entry'
    call debugmsg

! Receive metadata.

    write(msg,'(a,a,a,a,i0,a,i0,a)') this,' (',tasktype,' ',me,&
      '): receiving variable metadata for its=',its,'.'
    call debugmsg
    do ct=0,nct-1
      ctbuffer=ct+1
      do isegment=1,segments
        call append_to_list(ctbuffer) ! append an empty buffer node
        node=>buffers(ctbuffer)%current
        call mpi_recv(node%vm,1,comm_varmeta,ct,tag_metadata,intercomm,&
          mpi_status_ignore,istatus)
        write(msg,'(a,a,a,a,a,i0,a,i0,a)') 'ERROR: ',this,' (',tasktype,' ',me,&
          '): MPI_Recv metadata returned ',istatus,'.'
        call die_if(istatus.ne.mpi_success,istatus,msg)
        write(msg,'(a,a,a,a,i0,a,i0,a,a,a,i0,a)') this,' (',tasktype,' ',me,&
          '): its=',its,' received ',node%vm%varname(1:varname_len),&
          ' metadata: allocating ',node%vm%segment_size,&
          ' bytes for variable data.'
        call debugmsg
        allocate(node%segment(node%vm%segment_size),stat=istatus)
        write(msg,'(a,a,a,a,a,i0,a,a,a)') 'ERROR: ',this,' (',tasktype,' ',me,&
          '): Allocate failed for ',node%vm%varname(1:varname_len),' segment.'
        call die_if(istatus.ne.0,istatus,msg)
      enddo
    enddo
    write(msg,'(a,a,a,a,i0,a,i0,a)') this,' (',tasktype,' ',me,&
      '): received variable metadata for its=',its,'.'
    call debugmsg

! Receive variable data.

    write(msg,'(a,a,a,a,i0,a,i0,a)') this,' (',tasktype,' ',me,&
      '): receiving variable data for its=',its,'.'
    call debugmsg
    do ct=0,nct-1
      ctbuffer=ct+1
      node=>buffers(ctbuffer)%head
      do while (associated(node))
        call mpi_recv(node%segment,node%vm%segment_size,mpi_real,ct,tag_data,&
          intercomm,mpi_status_ignore,istatus)
        write (msg,'(a,a,a,a,a,i0,a,i0,a)') 'ERROR: ',this,' (',tasktype,' ',me,&
          '): MPI_Recv variable data returned ',istatus,'.'
        call die_if(istatus.ne.mpi_success,istatus,msg)
        write(msg,'(a,a,a,a,i0,a,i0,a,a,a)') this,' (',tasktype,' ',me,&
          '): its=',its,' received ',node%vm%varname(1:varname_len),&
          ' variable data.'
        call debugmsg
        node=>node%next
      enddo
    enddo
    write(msg,'(a,a,a,a,i0,a,i0,a)') this,' (',tasktype,' ',me,&
      '): received variable data for its=',its,'.'
    call debugmsg
    call write_vars_to_disk(its)
    call clear_list
    write(msg,'(a,a,a,a,i0,a,i0,a,i0,a)') this,' (',tasktype,' ',me,'): its=',&
      its,' segments=',segments,' exit'
    call debugmsg

  end subroutine process_frame
#endif /* SERIAL */

!--------------------------------------------------------------------------------
  subroutine var_to_disk(accum_start,ips,ims,ipe,ime,filename,header,its,levels,&
    scalefactor,time,var,varname)
!--------------------------------------------------------------------------------

! Write a global-size model field array to disk. This is a private routine -- not
! part of the icosio API. It is called by both compute and write tasks.

    character(len=*),intent(in)::filename,varname
    character,intent(in)::header(header_size)
    integer,intent(in)::accum_start,ips,ims,ipe,ime,its,levels,time
    real,intent(in)::scalefactor,var(levels,ims:ime)

    character(len=11)::this='var_to_disk'
    integer::index,ipn,ivl
    integer,save::count
    logical::append
    real,allocatable::reordered_var(:,:)

! If the list of filenames that have already been written to during this output
! frame has not already been allocated, do that now and zero its associated
! counter.
    
    if (.not.allocated(filename_list)) then
      allocate(filename_list(max_output_files),stat=istatus)
      write(msg,'(a,a,a,a,a,i0,a)') 'ERROR: ',this,' (',tasktype,' ',me,&
        '): failed to allocate filename_list.'
      call die_if(istatus.ne.0,istatus,msg)
      count=0
    endif

! Determine whether or not we need to append to the output file: If the current
! filename appears in the list of files we've already written to this frame,
! we will append; otherwise we will create a new file, potentially overwriting
! any existing file (as may be the case for restart runs).

    append=.false.
    do index=1,count
      if (trim(filename).eq.trim(filename_list(index))) then
        append=.true.
        exit
      endif
    enddo

! If we've decided not to append, check that we're not exceeding the size of
! our list of already-written-to filenames. If that's ok, increment the counter
! and record the current filename.

    if (.not.append) then
      count=count+1
      write(msg,'(a,a,a,a,a,i0,a)') 'ERROR: ',this,' (',tasktype,' ',me,&
        '): max_output_files limit exceeded.'
      call die_if(count.gt.max_output_files,istatus,msg)
      filename_list(count)=filename
    endif

! If non-grib binary history output is enabled:

    if (binout) then
      if (append) then
        open (lunout,file=filename(1:filename_len),form="unformatted",&
          iostat=istatus,position="append")
      else
        open (lunout,file=filename(1:filename_len),form="unformatted",&
          iostat=istatus)
      endif
      write(msg,'(a,a,a,a,a,i0,a,a,a)') 'ERROR: ',this,' (',tasktype,' ',me,&
        '): failed to open file ',filename(1:filename_len),'.'
      call die_if(istatus.ne.0,istatus,msg)
      write(msg,'(a,a,a,a,a)') 'ERROR: ',this,': Could not open file ',&
        filename(1:filename_len),'.'
      write (lunout) header
      if (associated(inv_perm_global)) then
        allocate(reordered_var(levels,ipe-ips+1))
        do ipn=ips,ipe
          do ivl=1,levels
            reordered_var(ivl,ipn)=var(ivl,inv_perm_global(ipn))
          enddo
        enddo
        write (lunout) reordered_var
        deallocate(reordered_var)
      else
        write (lunout) var(:,ips:ipe)
      endif
      write(msg,'(a,a,a,a,i0,a,a,a)') this,' (',tasktype,' ',me,&
        '): Wrote binary FIM field ',varname(1:varname_len),'.'
      call debugmsg
      close (lunout)
    endif

#ifndef NOGRIB
! If grib output is enabled:

    if (gribout) then
      call post_write_field(var,varname(1:varname_len),scalefactor,accum_start,&
        istatus)
      write(msg,'(a,a,a,a,a,i0,a)') 'ERROR: ',this,' (',tasktype,' ',me,&
        '): Bad return from post_write_field.'
      call die_if(istatus.ne.0,istatus,msg)
    endif
#endif /* NOGRIB */

  end subroutine var_to_disk

!--------------------------------------------------------------------------------
  subroutine write_vars_to_disk(its)
!--------------------------------------------------------------------------------

! Write to disk one output frame of all the variables for which this write task
! is responsible. If grib output is enabled, one one write task is allowed: The
! single write task open the grib file, writes all necessary data to it, then
! closes the file. This is a private routine -- not part of the icosio API. It
! is called only by write tasks.

! Array "glbvar" is allocated and deallocated inside of a while loop over
! variables, so its size neveer exceeds that of a single 3D global variable.

    integer,intent(in)::its

    type buffer_node_ptr
      type(buffer_node),pointer::ptr=>null()
    end type buffer_node_ptr

    character(len=18)::this='write_vars_to_disk'
    character(len=max_filename_len)::filename
    character(len=max_varname_len)::varname
    integer::accum_start,inode,ipn,ipns,ipnstart,istatus,levels,offset,ret,time
    real,pointer::glbvar(:,:),segment(:)
    real::scalefactor
    type(buffer_node),pointer::master
    type(buffer_node_ptr),pointer::nodes(:)

    write(msg,'(a,a,a,a,i0,a,i0,a)') this,' (',tasktype,' ',me,'): its=',its,&
      ' entry'
    call debugmsg

! Allocate buffer node pointers, one per compute task.

    allocate(nodes(nct),stat=istatus)
    write(msg,'(a,a,a,a,a,i0,a)') 'ERROR: ',this,' (',tasktype,' ',me,&
      '): Failed to allocate nodes.'
    call die_if(istatus.ne.0,istatus,msg)

! Point the buffer node pointers to the list heads.

    do inode=1,size(nodes)
      nodes(inode)%ptr=>buffers(inode)%head
    enddo

    master=>nodes(1)%ptr ! Set the master pointer

#ifndef NOGRIB
! If gribout enabled, open the grib file.

    if (gribout) then
      time=master%vm%time
      call post_init_file(time,istatus)
      write (msg,'(a,a,a,a,a,i0,a)') 'ERROR: ',this,' (',tasktype,' ',me,&
        '): post_init_file failed.'
      call die_if (istatus.ne.0,istatus,msg)

      write (msg,'(a,a,a,a,i0,a,i0,a)') this,' (',tasktype,' ',me,&
        '): post_init_file opened grib file time=',time,'.'
      call debugmsg
    endif
#endif /* NOGRIB */

! Allocate & populate the global variable buffer and have it written to disk.

    do while (associated(master)) ! One iteration per variable

      accum_start=master%vm%accum_start
      filename(1:filename_len)=master%vm%filename(1:filename_len)
      levels=master%vm%levels
      scalefactor=master%vm%scalefactor
      time=master%vm%time
      varname(1:varname_len)=master%vm%varname(1:varname_len)

      allocate(glbvar(levels,nip),stat=istatus)
      write(msg,'(a,a,a,a,a,i0,a)') 'ERROR: ',this,' (',tasktype,' ',me,&
        '): Failed to allocate glbvar.'
      call die_if(istatus.ne.0,istatus,msg)

! Copy data from buffer nodes into global buffer.

      ipnstart=1
      do inode=1,size(nodes)
        segment=>nodes(inode)%ptr%segment
        ipns=nodes(inode)%ptr%vm%segment_size/levels
        offset=1
        do ipn=1,ipns
          glbvar(1:levels,ipnstart-1+ipn)=segment(offset:offset-1+levels)
          offset=offset+levels
        enddo
        ipnstart=ipnstart+ipns
      enddo

      call var_to_disk(accum_start,1,1,nip,nip,filename,&
        master%vm%header(1:header_size),its,levels,scalefactor,time,glbvar,&
        varname)

! Step forward in buffer node list to handle another variable.

      do inode=1,size(nodes)
        nodes(inode)%ptr=>nodes(inode)%ptr%next
      enddo
      master=>nodes(1)%ptr

      deallocate(glbvar)

    enddo

! Reset the list of filenames written to during this output interval.

    if (allocated(filename_list)) deallocate(filename_list)

    write(msg,'(a,a,a,a,i0,a,i0,a)') this,' (',tasktype,' ',me,'): its=',its,&
      ' exit'
    call debugmsg

#ifndef NOGRIB
! Close the output grib file.

    if (gribout) then
      call post_finalize_file(ret)
      write (msg,'(a,a,a,a,i0,a,i0,a)') this,' (',tasktype,' ',me,&
        '): post_finalize_file closed grib file time=',time,'.'
      call debugmsg
    endif
#endif

  end subroutine write_vars_to_disk

end module icosio
