module gfs_dyn_iau_module

! module for iau data, IO and time interpolation.
 use ESMF
 use gfs_dyn_layout1, only: me, len_trie_ls, len_trio_ls, ls_dim, nodes, &
                    lats_node_a, me_l_0,lats_node_a_max
 use gfs_dyn_machine, only: kind_evod, kind_grid,kind_io4,kind_io8,kind_phys
 use namelist_dynamics_def, only: iaufiles_fg,iaufiles_anl,iaufhrs,iau,iau_delthrs,ens_nam,nemsio_in,semilag
 use gfs_dyn_resol_def
 use gfs_dyn_mpi_def, only: mc_comp,mpi_sum,mpi_real4,mpi_complex,mpi_r_io_r,mpi_real8
 use gfs_dynamics_states_mod, only: gfs_dynamics_internal_state
 USE gfs_dynamics_namelist_mod, only: GFS_Dyn_State_Namelist
 USE gfs_dynamics_grid_create_mod
 use gfs_dyn_tracer_config, only: gfs_dyn_tracer
 USE gfs_dynamics_add_get_state_mod
 USE gfs_dynamics_err_msg_mod
 USE gfs_dyn_coordinate_def, ONLY: ak5, bk5


! iaufiles_fg: filenames for first-guess fields.
! iaufiles_anl: filenames for analysis fields.
! iaufhrs: forecast hours for input files.
! iau_delthrs: length of IAU window (in hours).

 implicit none
 private

 public :: init_iau, destroy_iau, getiauforcing,applyiauforcing

 real(kind_evod), dimension(:,:,:,:),allocatable ::  grid_gr_iauall
 integer, public :: nfiles,gq_iau,gu_iau,gv_iau,gt_iau,grq_iau
 logical, public :: iau_initialized = .false.,iauinc=.false.



 contains

 subroutine init_iau(int_state)

! read in first-guess and analysis files, compute and store increments in
! spectral space.
 
   type(gfs_dynamics_internal_state),intent(in)  :: int_state   
   integer, allocatable, dimension(:) :: idt
   integer n,nfilesall,idate(4),ierr

   integer lan,k,lat,i,lon_dim,lons_lat

   gu_iau=1            ! u wind 
   gv_iau=gu_iau+levs  ! v wind
   gt_iau=gv_iau+levs  ! virtual temperature
   grq_iau=gt_iau+levs ! tracers
   gq_iau=grq_iau+levh    ! surface pressure
   if (me.EQ.me_l_0) print*,'gu_iau,gv_iau,gt_iau=',gu_iau,gv_iau,gt_iau
   if (me.EQ.me_l_0) print*,'grq_iau,gq_iau=',grq_iau,gq_iau

   iau_initialized = .true.
   iauinc=.false.
   if (trim(iaufiles_fg(1)) .eq. '' .and. trim(iaufiles_anl(1)) .ne. '') then
      iauinc=.true.
   endif

   nfilesall = size(iaufiles_anl)
   nfiles = 0
   do n=1,nfilesall
      if (trim(iaufiles_anl(n)) .eq. '' .or. iaufhrs(n) .lt. 0) exit
      if (me .eq. me_l_0) then
         print *,n,trim(adjustl(iaufiles_anl(n)))
          if (.not. iauinc) print *,n,trim(adjustl(iaufiles_fg(n)))
      endif
      nfiles = nfiles + 1
   enddo
   if (me .eq. me_l_0) print *,'nfiles = ',nfiles
   call mpi_barrier(mc_comp,ierr)
   if (nfiles < 1) then
      print *,'must be at least one file in iaufiles_fg and iaufiles_anal'
      call mpi_quit(9999)
   endif
   allocate(idt(nfiles-1))
   idt = iaufhrs(2:nfiles)-iaufhrs(1:nfiles-1)
   do n=1,nfiles-1
      if (idt(n) .ne. iaufhrs(2)-iaufhrs(1)) then
        print *,'forecast intervals in iaufhrs must be constant'
        call mpi_quit(9999)
      endif
   enddo
   if (me .eq. me_l_0) print *,'iau interval = ',iau_delthrs,' hours'
   deallocate(idt)
   allocate(grid_gr_iauall(lonf,lats_node_a_max,gq_iau,nfiles))
   if (nemsio_in) then
       call read_iau_nemsio(int_state)
   else
       call read_iau_sigio(int_state)
   endif
 end subroutine init_iau
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine read_iau_nemsio(int_state)
! read in first-guess and analysis files, compute and store increments in
! grid space.
   implicit none 
   type(gfs_dynamics_internal_state),intent(in)  :: int_state   
   real(kind_grid), dimension(:,:,:),allocatable ::  grid_gr_fg
   character(len=120) filename
   integer n,iprint,idate(4),ierr

   REAL(KIND=KIND_GRID) ptrc   (lonf,lats_node_a,ntrac)          !glbsum

   iprint = 1
   if (me.eq.me_l_0) print*,'reading iau_nemsio',iauinc,me
   if (.not. iauinc) allocate(grid_gr_fg(lonf,lats_node_a_max,gq_iau))
   if (me.eq.me_l_0) print*,'allocated',lonf,lats_node_a_max,gq_iau
   do n=1,nfiles
      if (.not. iauinc) then ! skip reading first guess and only read in increment file
         filename = iaufiles_fg(n)
         if (int_state%ens) then
           filename = trim(filename) // ens_nam
         endif
         if (me .eq. me_l_0) print *,'reading fg ',trim(filename)
         CALL TREADEO_nemsio_iau(filename,IDATE,                 &
              grid_gr_fg(:,1:lats_node_a,gq_iau),                  &
              grid_gr_fg(:,1:lats_node_a,gu_iau:gu_iau+levs-1),    &
              grid_gr_fg(:,1:lats_node_a,gv_iau:gv_iau+levs-1),    &
              grid_gr_fg(:,1:lats_node_a,gt_iau:gt_iau+levs-1),    &
              grid_gr_fg(:,1:lats_node_a,grq_iau:grq_iau+levh-1),      &
              int_state%LS_NODE,int_state%LS_NODES,int_state%MAX_LS_NODES,IPRINT,               &
              int_state%global_lats_a,int_state%lats_nodes_a,int_state%lonsperlat)
         if (me .eq. me_l_0)print*,'fg grid_gr ps',me, grid_gr_fg(1,1,gq_iau),grid_gr_fg(2,lats_node_a,gq_iau)
         if(iaulnp.EQ.1) grid_gr_fg(:,:,gq_iau)=alog(grid_gr_fg(:,:,gq_iau))
      endif ! not iauinc
      filename = iaufiles_anl(n)   ! analysis or increment file
      if (int_state%ens) then
         filename = trim(filename) // ens_nam
      endif
      if (me .eq. me_l_0) print *,'reading anl ',trim(filename)
      CALL TREADEO_nemsio_iau(filename,IDATE,                    &
           grid_gr_iauall(:,1:lats_node_a,gq_iau,n),               &
           grid_gr_iauall(:,1:lats_node_a,gu_iau:gu_iau+levs-1,n), &
           grid_gr_iauall(:,1:lats_node_a,gv_iau:gv_iau+levs-1,n), &
           grid_gr_iauall(:,1:lats_node_a,gt_iau:gt_iau+levs-1,n), &
           grid_gr_iauall(:,1:lats_node_a,grq_iau:grq_iau+levh-1,n),   &
           int_state%LS_NODE,int_state%LS_NODES,int_state%MAX_LS_NODES,IPRINT,                 &
           int_state%global_lats_a,int_state%lats_nodes_a,int_state%lonsperlat)
      if(iaulnp.EQ.1) grid_gr_iauall(:,:,gq_iau,n)=alog(grid_gr_iauall(:,:,gq_iau,n))
      if (me.EQ.me_l_0) print*,'anl grid_gr ps',me, grid_gr_iauall(1,1,gq_iau,n),grid_gr_iauall(2,lats_node_a,gq_iau,n)
!!$omp workshare     
      if (.not. iauinc)grid_gr_iauall(:,:,:,n)=grid_gr_iauall(:,:,:,n)-grid_gr_fg(:,:,:)
!!$omp end workshare      
  enddo
  if (me.EQ.me_l_0) print*,'iau grid_gr ps',me, grid_gr_iauall(1,1,gq_iau,1),grid_gr_iauall(2,lats_node_a,gq_iau,1)
  call mpi_barrier(mc_comp,ierr)
  if (.not.iauinc) deallocate(grid_gr_fg)
 end subroutine read_iau_nemsio
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 subroutine read_iau_sigio(int_state)

! read in first-guess and analysis files, compute and store increments in
! spectral space.
 
   implicit none 
   type(gfs_dynamics_internal_state),intent(in)  :: int_state   
   real(kind_grid), dimension(:,:,:),allocatable ::  grid_gr_fg
   real(kind_grid), dimension(:,:,:),allocatable ::  grid_gr_anl
   real(kind_grid), dimension(:,:)  ,allocatable ::  grid_gr_ps
   real(kind_grid), dimension(:,:,:),allocatable ::  grid_gr_tr
   real(kind_grid), dimension(:,:,:),allocatable ::  grid_gr_dummy
   character(len=120) filename
   integer n,iprint,idate(4),ierr

   integer N1,lan,k,lat,i,lon_dim,lons_lat,j

   iprint=1
   if (.not. iauinc) allocate(grid_gr_fg(lonf,lats_node_a_max,gq_iau))
   allocate(grid_gr_anl(lonf,lats_node_a_max,gq_iau))
   allocate(grid_gr_tr(lonf,lats_node_a_max,levh))
   allocate(grid_gr_ps(lonf,lats_node_a_max))
   allocate(grid_gr_dummy(lonf,lats_node_a_max,levs))
   do n=1,nfiles
      if (.not. iauinc) then ! skip reading first guess and only read in increment file
         filename = iaufiles_fg(n)
         if (int_state%ens) then
           filename = trim(filename) // ens_nam
         endif
         if (me .eq. me_l_0) print *,'reading ',trim(filename)
         CALL TREADEO_IAU(IDATE,grid_gr_fg,              &
                      int_state%LS_NODE,int_state%LS_NODES,int_state%MAX_LS_NODES,IPRINT,    &
                      int_state%global_lats_a,int_state%lats_nodes_a,int_state%lonsperlat, filename,   &
                      int_state%epse, int_state%epso, int_state%epsedn,  &
                      int_state%epsodn, int_state%plnew_a, int_state%plnow_a,                &
                      int_state%plnev_a, int_state%plnod_a,int_state%snnp1ev, int_state%snnp1od)
         if (me .eq. me_l_0) print *,'done reading ',trim(filename)
         if (me .eq. me_l_0) print *,'calling model_to_common_vars'
         call model_to_common_vars (grid_gr_fg(1,1,gq_iau),  &
                                    grid_gr_fg(1,1,gt_iau),  &
                                    grid_gr_fg(1,1,grq_iau), &
                                    grid_gr_fg(1,1,gu_iau),  &
                                    grid_gr_fg(1,1,gv_iau),  &
                                    grid_gr_dummy(1,1,1),    &
                                    grid_gr_dummy(1,1,1),    &
                                    grid_gr_dummy(1,1,1),    &
                                    int_state%global_lats_a,          &
                                    int_state%lonsperlat,1)
         if(iaulnp.EQ.1) grid_gr_fg(:,:,gq_iau)=alog(grid_gr_fg(:,:,gq_iau))
         if (me .eq. me_l_0) print *,'back from model_to_common_vars'
         call mpi_barrier(mc_comp,i)
      endif
      filename = iaufiles_anl(n)   ! analysis or increment file
      if (int_state%ens) then
         filename = trim(filename) // ens_nam
      endif
      if (me.eq.me_l_0) print*,'calling treadeo',trim(filename)
      CALL TREADEO_IAU(IDATE,grid_gr_anl,                                               &
                   int_state%LS_NODE,int_state%LS_NODES,int_state%MAX_LS_NODES,IPRINT,  &
                   int_state%global_lats_a,int_state%lats_nodes_a,int_state%lonsperlat, filename,   &
                   int_state%epse, int_state%epso, int_state%epsedn,                                &
                   int_state%epsodn, int_state%plnew_a, int_state%plnow_a,                          &
                   int_state%plnev_a, int_state%plnod_a,int_state%snnp1ev, int_state%snnp1od)
      if (me.eq.me_l_0) print*,'back from treadeo',trim(filename)
      if (iauinc) then
        grid_gr_tr=int_state%grid_gr(:,:,g_rq:g_rq+  levh-1)
      else
        grid_gr_tr=grid_gr_anl(:,:,grq_iau:grq_iau+  levh-1)
      endif
      grid_gr_ps=grid_gr_anl(:,:,gq_iau)
      call model_to_common_vars (grid_gr_ps(1,1),          &
                                 grid_gr_anl(1,1,gt_iau),  &
                                 grid_gr_tr(1,1,1),        &
                                 grid_gr_anl(1,1,gu_iau),  &
                                 grid_gr_anl(1,1,gv_iau),  &
                                 grid_gr_dummy(1,1,1),     &
                                 grid_gr_dummy(1,1,1),     &
                                 grid_gr_dummy(1,1,1),     &
                                 int_state%global_lats_a,  &
                                 int_state%lonsperlat,1)
     ! copy back converted pressure if not increments
     if (.not. iauinc) grid_gr_anl(:,:,gq_iau)=grid_gr_ps
     if(iaulnp.EQ.1) grid_gr_anl(:,:,gq_iau)=alog(grid_gr_anl(:,:,gq_iau))
     print*,'b4 grid_gr ps',me, grid_gr_anl(1,1,gq_iau),grid_gr_fg(1,1,gq_iau)
!!$omp workshare     
     if (iauinc) then
         grid_gr_iauall(:,:,gq_iau,n)=grid_gr_anl(:,:,gq_iau)
         grid_gr_iauall(:,:,gu_iau:gu_iau+levs-1,n)=grid_gr_anl(:,:,gu_iau:gu_iau+levs-1)
         grid_gr_iauall(:,:,gv_iau:gv_iau+levs-1,n)=grid_gr_anl(:,:,gv_iau:gv_iau+levs-1)
         grid_gr_iauall(:,:,gt_iau:gt_iau+levs-1,n)=grid_gr_anl(:,:,gt_iau:gt_iau+levs-1)
         grid_gr_iauall(:,:,grq_iau:grq_iau+levh-1,n)=grid_gr_anl(:,:,grq_iau:grq_iau+levh-1)
     else
         grid_gr_iauall(:,:,gq_iau,n)=grid_gr_anl(:,:,gq_iau)-grid_gr_fg(:,:,gq_iau)
         grid_gr_iauall(:,:,gu_iau:gu_iau+levs-1,n)=grid_gr_anl(:,:,gu_iau:gu_iau+levs-1)-grid_gr_fg(:,:,gu_iau:gu_iau+levs-1)
         grid_gr_iauall(:,:,gv_iau:gv_iau+levs-1,n)=grid_gr_anl(:,:,gv_iau:gv_iau+levs-1)-grid_gr_fg(:,:,gv_iau:gv_iau+levs-1)
         grid_gr_iauall(:,:,gt_iau:gt_iau+levs-1,n)=grid_gr_anl(:,:,gt_iau:gt_iau+levs-1)-grid_gr_fg(:,:,gt_iau:gt_iau+levs-1)
         grid_gr_iauall(:,:,grq_iau:grq_iau+levh-1,n)=grid_gr_anl(:,:,grq_iau:grq_iau+levh-1)-grid_gr_fg(:,:,grq_iau:grq_iau+levh-1)
     endif
!!$omp end workshare      
  enddo
  if (me.EQ.me_l_0) print*,'iau grid_gr ps',me, grid_gr_iauall(1,1,gq_iau,1),grid_gr_iauall(2,lats_node_a,gq_iau,1)
  if (me.EQ.me_l_0) print*,'iau grid_gr t',me, grid_gr_iauall(1,1,gt_iau+30,1),grid_gr_iauall(2,lats_node_a,gt_iau+30,1)
  if (me.EQ.me_l_0) print*,'iau grid_gr u',me, grid_gr_iauall(1,1,gu_iau+30,1),grid_gr_iauall(2,lats_node_a,gu_iau+30,1)
  if (me.EQ.me_l_0) print*,'iau grid_gr q',me, grid_gr_iauall(1,1,grq_iau+30,1),grid_gr_iauall(2,lats_node_a,grq_iau+30,1)
  if (.not.iauinc) deallocate(grid_gr_fg)
  deallocate(grid_gr_anl)
  deallocate(grid_gr_ps,grid_gr_tr)
 end subroutine read_iau_sigio
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 subroutine getiauforcing(grid_gr_iau,t,int_state,rc)
        
! compute an IAU forcing by interpolating increments to model time set
! and dividing by length of IAU window (in seconds).
      
   implicit none 
   REAL(KIND=KIND_GRID),intent(out) :: grid_gr_iau(lonf,lats_node_a_max,gq_iau)
   type(gfs_dynamics_internal_state),intent(in)  :: int_state   
   real(kind_evod), intent(in) :: t
   INTEGER,intent(out) :: rc
   real(kind_evod) delt, dt
   integer n,t1,t2
   integer i,j,k,lat
   !real(kind_io4), allocatable, dimension(:,:) :: workg,workg_out
   !real(kind_phys), allocatable, dimension(:,:) :: grid_gr_tmp
   !real (kind=kind_io8)   glolal(lonf,lats_node_a)
   !integer kmsk0(lonf,lats_node_a)
   !kmsk0 = 0
  
   grid_gr_iau = 0.
   dt = iau_delthrs*3600.
   if (me .eq. me_l_0) print *,'in getiauforcing1',nfiles,iaufhrs(1:nfiles)
   if (me .eq. me_l_0) print *,'in getiauforcing2',t, iaufhrs(1)*3600,iaufhrs(nfiles)*3600.
   ! set forcing to zero and return if outside iau window.
   if ( nfiles > 1) then  ! IAU forcing files bookend interval
      if (t <= iaufhrs(1)*3600. .or. t > iaufhrs(nfiles)*3600.) then
         if (me .eq. me_l_0) print *,'no iau forcing'
         rc=1
         return
      endif
   else  ! single file at middle of window
      t1=iaufhrs(1)*3600 - dt*0.5
      t2=iaufhrs(1)*3600 + dt*0.5
      if ( t <= t1 .or. t > t2 ) then
         if (me .eq. me_l_0) print *,'no iau forcing'
         rc=1
         return
      endif
   endif
   if (nfiles > 1) then
      if (t .eq. 3600.*iaufhrs(nfiles)) then
!!$omp workshare
         grid_gr_iau = grid_gr_iauall(:,:,:,nfiles)/dt
!!$omp end workshare
         return
      else if (t .eq. 3600.*iaufhrs(1)) then
!!$omp workshare
         grid_gr_iau = grid_gr_iauall(:,:,:,1)/dt
!!$omp end workshare
         return
      endif
      do n=1,nfiles
         if (iaufhrs(n)*3600. > t) exit
      enddo
      if (me .eq. me_l_0) print *,'n,t,to',n,t/3600.,iaufhrs(n)
      delt = (iaufhrs(n)-(t/3600.))/(iaufhrs(n)-iaufhrs(n-1))
!!$omp workshare
      grid_gr_iau = ((1.-delt)*grid_gr_iauall(:,:,:,n) + delt*grid_gr_iauall(:,:,:,n-1))/dt
!!$omp end workshare
      if (me .eq. me_l_0) print *,'getiauforcing:',t/3600.,1.-delt,n,iaufhrs(n),delt,n-1,iaufhrs(n-1)
   else
      grid_gr_iau = grid_gr_iauall(:,:,:,1)/dt
   endif
   if (me .eq. me_l_0) print *,'have iau forcing',grid_gr_iau(1,1,gq_iau)
   rc=0
   !allocate(workg(lonf,latg))
   !allocate(workg_out(lonf,latg))
   !allocate(grid_gr_tmp(lonf,lats_node_a_max))
! write out data
   !if (me.eq.me_l_0) open(77,file='../iau_forcing.dat',form='unformatted')
   !do k=1,gq_iau
   !  workg= 0.
   !  grid_gr_tmp(:,:)=grid_gr_iauall(:,:,k,1)
   !  CALL uninterpred(2,kmsk0,glolal,grid_gr_tmp,&
   !                   int_state%global_lats_a,int_state%lonsperlat)
   !  do j=1,lats_node_a
   !     lat=int_state%global_lats_a(ipt_lats_node_a-1+j)
   !     do i=1,lonf
   !        workg(i,lat) = glolal(i,j)
   !     enddo
   !  enddo
   !  call mpi_reduce(workg,workg_out,lonf*latg,&
   !                mpi_real4,mpi_sum,me_l_0,mc_comp,i)
   !  if (me .eq. me_l_0) then
   !     write(77) workg_out
   !     print*,'array a k',grid_gr_iauall(1,1,k,1),glolal(1,1)
   !     print*,'writing k',workg(1,1),workg_out(1,1)
   !  endif
   !enddo
   !close(77)
   !deallocate(workg,workg_out)
   !call mpi_barrier(mc_comp,i)
   !call mpi_quit(9999)
! end write
 end subroutine getiauforcing

SUBROUTINE applyiauforcing(grid_gr_iau,imp_gfs_dyn,int_state,deltim,rc)

! this subroutine added the iau forcing term to the gfs import state
!-----------------------------------------------------------

      implicit none 
      type(esmf_state),    intent(inout) :: imp_gfs_dyn 
      type(gfs_dynamics_internal_state),intent(in)  :: int_state   
      REAL(KIND=KIND_GRID),intent(in) :: grid_gr_iau(lonf,lats_node_a_max,gq_iau)
      REAL(kind_evod),intent(in)     :: deltim
      INTEGER,          OPTIONAL,                 INTENT(out)   :: rc     

      TYPE(ESMF_Field)         :: Field
      TYPE(ESMF_FieldBundle)   :: Bundle
      TYPE(GFS_Dyn_State_Namelist) :: cf

      INTEGER                  :: rc1, rcfinal,irec
      INTEGER                  :: k, kstr, kend,i,j
      INTEGER                  :: mstr, mend,lat,ierr

      REAL, DIMENSION(:, :),    POINTER :: FArr2D
      REAL, DIMENSION(:, :, :), POINTER :: FArr3D
  ! real(kind_io4), allocatable, dimension(:,:) :: workg,workg_out
  ! real(kind_phys), allocatable, dimension(:,:) :: grid_gr_tmp
  ! real (kind=kind_io8)   glolal(lonf,lats_node_a)
  ! integer kmsk0(lonf,lats_node_a)
  ! character*3,fstr
  ! kmsk0 = 0

! initialize the error signal variables.
!---------------------------------------
      rc1     = esmf_success
      rcfinal = esmf_success

      CALL esmf_logwrite(						&
           " update import state with IAU forcing", 	&
            ESMF_LOGMSG_INFO, rc = rc1)

      cf = int_state%esmf_sta_list

! get the surface pressure array from the esmf import state.
! for the detailed comments for every computational steps
! please refer to the surface orography array code.
!-----------------------------------------------------------
      !write(fstr,FMT='(I3.3)') int_state%kdt
      !if (me.eq.me_l_0) open(78,file='../psiau_'//fstr//'.dat',form='unformatted')
      !allocate(workg(lonf,latg))
      !allocate(workg_out(lonf,latg))
      !allocate(grid_gr_tmp(lonf,lats_node_a_max))
! write out data
     !workg= 0.
     !grid_gr_tmp(:,:)=grid_gr_iau(:,:,gq_iau)
     !CALL uninterpred(2,kmsk0,glolal,grid_gr_tmp,&
     !                 int_state%global_lats_a,int_state%lonsperlat)
     !  do j=1,lats_node_a
     !     lat=int_state%global_lats_a(ipt_lats_node_a-1+j)
     !     do i=1,lonf
     !        workg(i,lat) = glolal(i,j)
     !     enddo
     !  enddo
     !  call mpi_reduce(workg,workg_out,lonf*latg,&
     !              mpi_real4,mpi_sum,me_l_0,mc_comp,i)
     !  if (me .eq. me_l_0) then
     !     write(78) workg_out
     ! endif
         IF(cf%ps_import == 1) THEN
            CALL getf90arrayfromstate(imp_gfs_dyn, 'ps', FArr2D, 0, rc = rc1)
            CALL gfs_dynamics_err_msg(rc1,"get esmf state - ps",rcfinal)
! write out data
     !    workg= 0.
     !    grid_gr_tmp(:,:)=FArr2D
     !    CALL uninterpred(2,kmsk0,glolal,grid_gr_tmp,&
     !                     int_state%global_lats_a,int_state%lonsperlat)
     !    do j=1,lats_node_a
     !       lat=int_state%global_lats_a(ipt_lats_node_a-1+j)
     !       do i=1,lonf
     !          workg(i,lat) = glolal(i,j)
     !       enddo
     !    enddo
     !call mpi_reduce(workg,workg_out,lonf*latg,&
     !              mpi_real4,mpi_sum,me_l_0,mc_comp,i)
     !if (me .eq. me_l_0) then
     !   write(78) workg_out
     !endif
!! add iau forcing
             if (me .eq. me_l_0) print *,'before forcing ps:',farr2d(4,4)
             If (iaulnp.GT.0) THEN
             If (iaulnp.EQ.1) THEN
                 FArr2D = exp(alog(FArr2D) + grid_gr_iau(:,:,gq_iau) * deltim)
             ELSE
                IF (iauinc .and. .not.semilag) THEN
                   FArr2D = FArr2D * (1.0 + grid_gr_iau(:,:,gq_iau) * deltim) 
!             
                ELSE
                    FArr2D = FArr2D + grid_gr_iau(:,:,gq_iau) * deltim 
                ENDIF
             ENDIF
             ENDIF
             if (me .eq. me_l_0) print *,'applyiauforcing ps:',grid_gr_iau(4,4,gq_iau),farr2d(4,4)
          END IF
! write out data
     !workg= 0.
     !grid_gr_tmp(:,:)=FArr2D
     !CALL uninterpred(2,kmsk0,glolal,grid_gr_tmp,&
     !                 int_state%global_lats_a,int_state%lonsperlat)
     !do j=1,lats_node_a
     !   lat=int_state%global_lats_a(ipt_lats_node_a-1+j)
     !   do i=1,lonf
     !      workg(i,lat) = glolal(i,j)
     !   enddo
     !enddo
     !call mpi_reduce(workg,workg_out,lonf*latg,&
     !              mpi_real4,mpi_sum,me_l_0,mc_comp,i)
     !if (me .eq. me_l_0) then
     !   write(78) workg_out
     !   close(78)
     !endif
   !deallocate(workg,workg_out)

! get the temperature array from the esmf import state.
!------------------------------------------------------
          IF(cf%temp_import == 1) THEN
             CALL getf90arrayfromstate(imp_gfs_dyn, 't', FArr3D, 0, rc = rc1)
             if (me .eq. me_l_0) print *,'before forcing t:',farr3d(4,4,4)
             CALL gfs_dynamics_err_msg(rc1,"get esmf state - t",rcfinal)
             FArr3D = FArr3D + grid_gr_iau(:,:,gt_iau:gt_iau+levs-1) * deltim
             if (me .eq. me_l_0) print *,'applyiauforcing t:',grid_gr_iau(4,4,gt_iau+3),farr3d(4,4,4)
          END IF

! get the zonal-wind array from the esmf import state.
! for detailed line by line comments please refer to 
! the temperature code.
!-----------------------------------------------------
          IF(cf%u_import == 1) THEN
             CALL getf90arrayfromstate(imp_gfs_dyn, 'u', FArr3D, 0, rc = rc1)
             if (me .eq. me_l_0) print *,'before forcing u:',farr3d(4,4,4)
             CALL gfs_dynamics_err_msg(rc1,"get esmf state - u",rcfinal)
             FArr3D = FArr3D + grid_gr_iau(:,:,gu_iau:gu_iau+levs-1) * deltim
             if (me .eq. me_l_0) print *,'applyiauforcing u:',grid_gr_iau(4,4,gu_iau+3),farr3d(4,4,4)
          END IF

! get the meridian-wind array from the esmf import state.
!-----------------------------------------------------
          IF(cf%v_import == 1) THEN
             CALL getf90arrayfromstate(imp_gfs_dyn, 'v', FArr3D, 0, rc = rc1)
             if (me .eq. me_l_0) print *,'before forcing v:',farr3d(4,4,4)
             CALL gfs_dynamics_err_msg(rc1,"get esmf state - v",rcfinal)
             FArr3D = FArr3D + grid_gr_iau(:,:,gv_iau:gv_iau+levs-1) * deltim
             if (me .eq. me_l_0) print *,'applyiauforcing v:',grid_gr_iau(4,4,gv_iau+3),farr3d(4,4,4)
          END IF

! get the tracer array from the esmf import state.
!-------------------------------------------------
          IF(cf%tracer_import == 1) THEN
             CALL ESMF_StateGet(imp_gfs_dyn, 'tracers', Bundle, rc = rc1 )
             CALL gfs_dynamics_err_msg(rc1, 'retrieve Ebundle from state', rcfinal)
             DO k = 1, int_state%ntrac
                 IF(ASSOCIATED(FArr3D)) NULLIFY(FArr3D)
                 CALL ESMF_FieldBundleGet(Bundle,                               &
                                          trim(gfs_dyn_tracer%vname(k, 1)), &
                                          field = Field,                        &
                                          rc = rc1)
                 CALL gfs_dynamics_err_msg(rc1, 'retrieve Efield from bundle', rcfinal)
                         
                 CALL ESMF_FieldGet(Field, farrayPtr = FArr3D, localDE = 0, rc = rc1)
    
                 CALL gfs_dynamics_err_msg(rc1, 'retrieve Farray from field', rcfinal)
                 kstr = grq_iau + (k-1)*levs
                 kend = kstr + levs - 1
                 if (me .eq. me_l_0) print *,'before forcing tr:',k,farr3d(4,4,4)
                 FArr3D = FArr3D + grid_gr_iau(:,:,kstr:kend) * deltim
                 if (me .eq. me_l_0) print *,'applyiauforcing tr:',grid_gr_iau(4,4,kstr+3),farr3d(4,4,4)
             END DO
         END IF

!
!
! PRINT out the final error signal message and put it to rc.
!-----------------------------------------------------------
      IF(PRESENT(rc)) CALL gfs_dynamics_err_msg_final(rcfinal, &
          "applyiauforcing",rc)

      END SUBROUTINE applyiauforcing


 subroutine destroy_iau()

   implicit none 
! deallocate array.
   deallocate(grid_gr_iauall)

 end subroutine destroy_iau

end module gfs_dyn_iau_module
