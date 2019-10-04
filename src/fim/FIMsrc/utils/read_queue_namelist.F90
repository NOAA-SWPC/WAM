module read_queue_namelist

  use module_initial_chem_namelists
  implicit none
  save

  character(8)   :: ComputeTasks = '10'            ! Number of compute tasks for FIM; 'S' means Serial
  character(8)   :: MaxQueueTime = '00:05:00'      ! Run time for the complete job (HH:MM:SS)
  character(160) :: SRCDIR  = '../FIMsrc'          ! Location of the FIM source
  character(160) :: PREPDIR = 'nodir'              ! If exists, use for prep otherwise calculate prep
  character(160) :: FIMDIR  = 'nodir'              ! If exists, use for FIM otherwise calculate FIM
!JR Changed paths to something non-existent, so namelist value is always used
  character(160) :: DATADIR = '/no_such_path'      ! Location of gfsltln and global_mtnvar files
  character(160) :: DATADR2 = '/no_such_path'      ! Location of the sanl file and the sfcanl file
  character(160) :: chem_datadir = '/no_such_path' ! Location of the chemistry data files
  integer        :: glvl                           ! The grid level defined in the Makefile
  integer        :: SubdivNum(20)                  ! Subdivision numbers for each recursive refinement(2: bisection, 3: trisection, etc.)
  integer        :: nvl                            ! Number of vertical native levels

contains

  subroutine ReadQUEUEnamelist
    integer::ierr
    logical,save::alreadyReadNamelists = .false.
    namelist /QUEUEnamelist/ ComputeTasks,MaxQueueTime,SRCDIR,PREPDIR,FIMDIR,DATADIR,DATADR2,chem_datadir
    namelist /CNTLnamelist/ glvl, SubdivNum, nvl
! TODO:  Remove duplicate namelist declarations via use of include 
! TODO:  or use association.  Or, better yet, just call FIM routines...  
! Note:  REWIND required by IBM!  
! TODO:  Using open-read-close in place of REWIND until SMS is updated
!JR add status and action to ensure 0-byte file doesn't get created by accident!
    if (.not.alreadyReadNamelists) then
      OPEN (10,file="FIMnamelist",status='old',action='read',iostat=ierr)
      if (ierr.ne.0) then
        write (0,'(a)') 'ERROR: Could not open FIMnamelist for read.'
        stop
      endif
      READ (10,NML=QUEUEnamelist,iostat=ierr)
      if (ierr.ne.0) then
        write (0,'(a)') 'ERROR: Could not read QUEUEnamelist namelist.'
        stop
      endif
      CLOSE(10)
      OPEN  (10,file="FIMnamelist",status='old',action='read',iostat=ierr)
      if (ierr.ne.0) then
        write (0,'(a)') 'ERROR: Could not open FIMnamelist for read.'
        stop
      endif
      READ (10,NML=CNTLnamelist,iostat=ierr)
      if (ierr.ne.0) then
        write (0,'(a)') 'ERROR: Could not read CNTLnamelist namelist.'
        stop
      endif
      CLOSE(10)
      OPEN  (10,file="FIMnamelist",status='old',action='read',iostat=ierr)
      if (ierr.ne.0) then
        write (0,'(a)') 'ERROR: Could not open FIMnamelist for read.'
        stop
      endif
      READ (10,NML=chemwrf,iostat=ierr)
      if (ierr.ne.0) then
        write (0,'(a)') 'NOTE: chemwrf namelist not read, continuing with chemistry disabled...'
      endif
      CLOSE(10)
      OPEN  (10,file="FIMnamelist",status='old',action='read',iostat=ierr)
      if (ierr.ne.0) then
        write (0,'(a)') 'ERROR: Could not open FIMnamelist for read.'
        stop
      endif
      READ (10,NML=wrfphysics,iostat=ierr)
      if (ierr.ne.0) then
        write (0,'(a)') 'NOTE: wrfphysics namelist not read, continuing...'
      endif
      close(10)
      alreadyReadNamelists = .true.
    endif
  end subroutine ReadQUEUEnamelist

  subroutine GetNprocs  (nprocs)
    integer,intent(OUT) :: nprocs
    call ReadQUEUEnamelist
    if(ComputeTasks=='S'.or.ComputeTasks=='s') then
      nprocs=1
    else
      read(ComputeTasks,*) nprocs
    endif
  end subroutine GetNprocs

  subroutine GetMaxQueueTime (QueueTime)
    character(8),intent(OUT) :: QueueTime
    call ReadQUEUEnamelist
    QueueTime = MaxQueueTime
  end subroutine GetMaxQueueTime

  subroutine ReturnGLVL (glvlout)
    integer,intent(OUT) :: glvlout
    call ReadQUEUEnamelist
    glvlout = glvl
  end subroutine ReturnGLVL

  subroutine ReturnNIP (nipout)
    integer,intent(OUT) :: nipout
    integer i
    call ReadQUEUEnamelist
    nipout = 1
    do i = 1, glvl
      nipout = nipout * SubdivNum(i)
    enddo
    nipout = 10 * nipout * nipout + 2
  end subroutine ReturnNIP

  subroutine ReturnSubdivNum(SubdivNumout)
    integer,intent(OUT) :: SubdivNumout(20)
    integer i
    call ReadQUEUEnamelist
    SubdivNumout = 0
    do i = 1, glvl
      SubdivNumout(i) = SubdivNum(i)
    enddo
  end subroutine ReturnSubdivNum

  subroutine ReturnNVL  (nvlout)
    integer,intent(OUT) :: nvlout
    call ReadQUEUEnamelist
    nvlout = nvl
  end subroutine ReturnNVL

  subroutine ReturnDT(dtout)
    REAL,intent(OUT) :: dtout
    integer i, rsl
    call ReadQUEUEnamelist
    rsl = 1
    do i = 1, glvl
      rsl = rsl * SubdivNum(i)
    enddo
    rsl =  INT(REAL(rsl) / 2.0)
    dtout = 1.2*4800./REAL(rsl)
  end subroutine ReturnDT

  subroutine GetChemOn (chem_on)
    logical,intent(OUT) :: chem_on
    call ReadQUEUEnamelist
!TODO:  need a better way to maintain this...  
    chem_on = (chem_opt /= 0)
  end subroutine GetChemOn

  subroutine GetWRFcuOn (wrfcu_on)
    logical,intent(OUT) :: wrfcu_on
    call ReadQUEUEnamelist
!TODO:  need a better way to maintain this...  
    wrfcu_on = (cu_physics /= 0)
  end subroutine GetWRFcuOn

  subroutine GetWRFmpOn (wrfmp_on)
    logical,intent(OUT) :: wrfmp_on
    call ReadQUEUEnamelist
!TODO:  need a better way to maintain this...  
    wrfmp_on = (mp_physics /= 0)
  end subroutine GetWRFmpOn

! returns .true. iff any WRF physics or chemistry is turned on
  subroutine GetWRFOn (wrf_on)
    logical,intent(OUT) :: wrf_on
    logical :: chem_on, wrfcu_on, wrfmp_on
    call GetChemOn(chem_on)
    call GetWRFcuOn(wrfcu_on)
    call GetWRFmpOn(wrfmp_on)
!TODO:  need a better way to maintain this...  
    wrf_on = (chem_on .OR. wrfcu_on .OR. wrfmp_on)
  end subroutine GetWRFOn

end module read_queue_namelist
