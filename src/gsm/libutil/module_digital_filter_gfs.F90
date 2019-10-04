!#include "../../../ESMFVersionDefine.h"

      module module_digital_filter_gfs
!
! a generic digital filter for any model under ESMF 
!
! March 2007	Hann-Ming Henry Juang
! February 2008 Weiyu Yang, updated to use the ESMF 3.1.0 library.
! November 2009 Jun Wang, digital filter is done on n+1 time step.
! February 2011 Weiyu Yang, Updated to use both the ESMF 4.0.0rp2 library,
!                           ESMF 5 library and the the ESMF 3.1.0rp2 library.
! May      2011 Weiyu Yang, Modified for using the ESMF 5.2.0r_beta_snapshot_07.
! Sep      2011 Weiyu Yang, Modified for using the ESMF 5.2.0r library.
! Oct 22   2015 S Moorthi   cleanup and optimize (vectorization)
! Nov ??   2015 S Moorthi   update for digital filter for semi-Lagrangian
!                           with gridded tracers
! Dec 15   2015 S Moorthi   Many changes to save 2d and 3d fields separately
!                           (thus reducing memory and fix a small memory leak)
!----------------------------------------------------------------------------
!
      USE ESMF

      implicit none

! ---------
! dynamics
! ---------
      real(8), allocatable,       save :: dyn_save_3d(:,:,:,:)
      real(8), allocatable,       save :: dyn_save_2d(:,:,:)
      real(8)             ,       save :: totalsum
      character(20), allocatable, save :: dyn_name(:)
      integer,       allocatable, save :: dyn_dim(:,:)
      integer,       allocatable, save :: dyn_items_ord2d(:),dyn_items_ord3d(:)
      integer,                    save :: kstep, nstep,           &
                                          lvl_change, items_2d, items_3d
      character(20)                       state_name
! ---------
! physics
! ---------
      type(esmf_state) ,          save :: phy_state_save
      character(20), allocatable, save :: phy_name(:)
      integer                             phy_items

      contains

! ---------------------------------------------------------------
! subroutine for dynamics
! ---------------------------------------------------------------
      subroutine digital_filter_dyn_init_gfs(dyn_state,ndfistep,dfilevs)
!
      implicit none
      type(esmf_state), intent(in) :: dyn_state
      integer,          intent(in) :: ndfistep, dfilevs
!
!     TYPE(ESMF_VM)                :: vm
      TYPE(ESMF_Field)             :: FIELD
      type(esmf_Grid)              :: GRID
      type(esmf_DistGrid)          :: distGRID
      type(ESMF_DELayout)          :: LAYOUT
      integer                      :: gridRank, dim_max(3)
      integer                      :: m,n,rc,ierr,dyn_items_tmp,  &
                                      me, nodes,nDes,deId,        &
                                      deList(1),i,j,k,nn

      INTEGER, DIMENSION(:, :), POINTER :: AL, AU
!
!      CALL ESMF_VMGetCurrent(vm, rc = rc)
!      CALL ESMF_VMGet(vm, localpet = me, petcount = nodes, rc = rc)

      nstep      = ndfistep 
      kstep      = - nstep -1
      lvl_change = dfilevs
      items_2d   = 0
      items_3d   = 0

      call esmf_stateget(state=dyn_state, name=state_name,       &
                         itemcount=dyn_items_tmp, rc=rc)

      allocate(dyn_name(dyn_items_tmp))
      allocate(dyn_items_ord2d(dyn_items_tmp), dyn_items_ord3d(dyn_items_tmp))
      allocate(dyn_dim(3,dyn_items_tmp))

      call esmf_stateget(state=dyn_state, itemnamelist=dyn_name, rc=rc)

      dim_max   = 0
      do n=1,dyn_items_tmp
        if(index(trim(dyn_name(n)),"_dfi") > 0) then

          CALL ESMF_StateGet(dyn_state, ItemName=trim(dyn_name(n)), field=Field, rc=rc)
          call esmf_Fieldget(FIELD, grid=grid, rc=rc)
!
          call ESMF_GridGet    (GRID, dimCount=gridRank, distGrid=distGrid, rc=rc)
          call ESMF_DistGridGet(distGRID, delayout=layout, rc=rc)
!         call ESMF_DELayoutGet(layout, deCount=nDEs, localDeList=deList, rc=rc )
          call ESMF_DELayoutGet(layout, deCount=nDEs, localDeToDeMap=deList, rc=rc )

          deId = deList(1)
          if(gridRank /= 3) then
            write(0,*)'WARNING: the rank of digital filter variables'  &
                   ,' is not 3! gridRank= ',gridRank
            stop
          endif
!
          allocate (AL(gridRank,0:nDEs-1), stat = ierr )
          allocate (AU(gridRank,0:nDEs-1), stat = ierr )
!
          call ESMF_DistGridGet(distgrid, minIndexPDe=AL,           &
                                maxIndexPDe=AU, rc=rc )

          do m=1,gridRank
            dyn_dim(m,n) = AU(m, deId) - AL(m, deId) + 1
            dim_max(m)   = max(dyn_dim(m,n),dim_max(m))
          enddo
          deallocate(AU, AL)

          if(trim(dyn_name(n)) == 'hs_dfi' .or.              &
             trim(dyn_name(n)) == 'ps_dfi') dyn_dim(3,n) = 1

          if (dyn_dim(3,n) > 1) then
            items_3d = items_3d + 1
            dyn_items_ord3d(items_3d) = n
          else
            items_2d = items_2d + 1
            dyn_items_ord2d(items_2d) = n
          endif
        endif
      enddo

      totalsum = 0.0

      if(minval(dim_max) > 1) then
        if (items_2d > 0) then
          allocate(dyn_save_2d(dim_max(1),dim_max(2),items_2d))
          do n=1,items_2d
!$omp parallel do private(i,j)
            do j=1,dim_max(2)
              do i=1,dim_max(1)
                dyn_save_2d(i,j,n) = 0.0
              enddo
            enddo
          enddo
        endif
        if (items_3d > 0) then
          allocate(dyn_save_3d(dim_max(1),dim_max(2),dim_max(3),items_3d))
          do n=1,items_3d
!$omp parallel do private(i,j,k)
            do k=1,dim_max(3)
              do j=1,dim_max(2)
                do i=1,dim_max(1)
                  dyn_save_3d(i,j,k,n) = 0.0
                enddo
              enddo
            enddo
          enddo
        endif
      endif
!      
      end subroutine digital_filter_dyn_init_gfs

! ---------------------------------------------------------------
      subroutine digital_filter_dyn_sum_gfs(dyn_state)
!
      implicit none
      type(esmf_state), intent(in)  :: dyn_state
!
!     TYPE(ESMF_VM)                                 :: vm
      TYPE(ESMF_Field)                              :: Field
      real(ESMF_KIND_R8), dimension(:,:),   pointer :: tmp_ptr2d
      real(ESMF_KIND_R8), dimension(:,:,:), pointer :: tmp_ptr
      real(8), parameter                            :: pi=acos(-1.0)
      real(8)                                       :: sx, wx, digfil
      integer                                       :: n, i, j, k, rc,item, &
                                                       dim1,dim2,dim3,nn    &
                                                      ,me, nodes
!     CALL ESMF_VMGetCurrent(vm, rc = rc)
!     CALL ESMF_VMGet(vm, localpet = me, petcount = nodes, rc = rc)

      kstep  = kstep + 1
      sx     = (pi*kstep) / nstep
      wx     = (pi*kstep) / (nstep+1)

      if( kstep /= 0 ) then
          digfil = (sin(wx)*sin(sx)) / (wx*sx)
      else
          digfil = 1
      endif 

      totalsum = totalsum + digfil

      if (items_2d > 0) then
        do n=1,items_2d
          item = dyn_items_ord2d(n)

          CALL ESMF_StateGet(dyn_state, ItemName=dyn_name(item), Field=FIELD, rc=rc)
          dim1 = dyn_dim(1,item)
          dim2 = dyn_dim(2,item)
          if (associated(tmp_ptr2D)) nullify(tmp_ptr2D)

          CALL ESMF_FieldGet(FIELD, localDe=0, farrayPtr=tmp_ptr2D, rc=rc)
!$omp parallel do private(i,j)
          do j=1,dim2
            do i=1,dim1
              dyn_save_2d(i,j,n) = dyn_save_2d(i,j,n) + digfil*tmp_ptr2D(i,j)
            enddo
          enddo

        enddo
      endif

      if (items_3d > 0) then
       do n=1,items_3d
          item = dyn_items_ord3d(n)

          CALL ESMF_StateGet(dyn_state, ItemName=dyn_name(item), Field=FIELD, rc=rc)
          dim1 = dyn_dim(1,item)
          dim2 = dyn_dim(2,item)
          dim3 = dyn_dim(3,item)

          if (associated(tmp_ptr)) nullify(tmp_ptr)

          CALL ESMF_FieldGet(FIELD, localDe=0, farrayPtr=tmp_ptr, rc=rc) 

          if(lvl_change /= dim3) then
            if(kstep == 0) then
              do k=lvl_change+1, dim3
!$omp parallel do private(i,j)
                do j=1,dim2
                  do i=1,dim1
                    dyn_save_3d(i,j,k,n) = tmp_ptr(i,j,k)
                  enddo
                enddo
              enddo
            endif
            do k=1,lvl_change
!$omp parallel do private(i,j)
              do j=1,dim2
                do i=1,dim1
                  dyn_save_3d(i,j,k,n) = dyn_save_3d(i,j,k,n)         &
                                       + digfil * tmp_ptr(i,j,k)
                enddo
              enddo
            enddo
          else
              
!$omp parallel do private(i,j,k)
            do k=1,dim3
              do j=1,dim2
                do i=1,dim1
                  dyn_save_3d(i,j,k,n) = dyn_save_3d(i,j,k,n)          &
                                       + digfil * tmp_ptr(i,j,k)
                enddo
              enddo
            enddo

          endif
        enddo
      endif
                                                    !     For tracers
      end subroutine digital_filter_dyn_sum_gfs

! ---------------------------------------------------------------
      subroutine digital_filter_dyn_average_gfs(dyn_state)
!
      implicit none
      type(esmf_state), intent(inout)               :: dyn_state
!
!     TYPE(ESMF_VM)                                 :: vm
      TYPE(ESMF_Field)                              :: FIELD
      real(ESMF_KIND_R8), dimension(:,:,:), pointer :: tmp_ptr
      real(ESMF_KIND_R8), dimension(:,:),   pointer :: tmp_ptr2D
!     real(ESMF_KIND_R8), dimension(:,:,:), pointer :: tmp_ptr1
      real(8)                                       :: totalsumi
      integer                                       :: n, i, j, k, rc, item &
                                                     , dim1,dim2,dim3,nn    &
                                                     , me, nodes
!
!     CALL ESMF_VMGetCurrent(vm, rc = rc)
!     CALL ESMF_VMGet(vm, localpet = me, petcount = nodes, rc = rc)

      totalsumi = 1.0 / totalsum

      if (items_2d > 0) then
        do n=1,items_2d
!
          item = dyn_items_ord2d(n)
          CALL ESMF_StateGet(dyn_state, ItemName=dyn_name(item), Field=FIELD, rc = rc)

          dim1 = dyn_dim(1,item)
          dim2 = dyn_dim(2,item)

          if (associated(tmp_ptr2D)) nullify(tmp_ptr2D)

          CALL ESMF_FieldGet(FIELD, localDe=0, farrayPtr=tmp_ptr2D, rc=rc)

!$omp parallel do private(i,j)
          do j=1,dim2
            do i=1,dim1
              tmp_ptr2D(i,j) = dyn_save_2d(i,j,n)*totalsumi
            enddo
          enddo
        enddo
      endif

      if (items_3d > 0) then
        do n=1,items_3d
!
          item = dyn_items_ord3d(n)
          CALL ESMF_StateGet(dyn_state, ItemName=dyn_name(item), Field=FIELD, rc = rc)

          dim1 = dyn_dim(1,item)
          dim2 = dyn_dim(2,item)
          dim3 = dyn_dim(3,item)

          if (associated(tmp_ptr)) nullify(tmp_ptr)

          CALL ESMF_FieldGet(FIELD, localDe=0, farrayPtr=tmp_ptr, rc=rc)

!       print *,'dfi dyn save field ',trim(dyn_name(item)),'=',  &
!       maxval(tmp_ptr),minval(tmp_ptr)

          if (lvl_change == dim3) then
!$omp parallel do private(i,j,k)
            do k=1,dim3
              do j=1,dim2
                do i=1,dim1
                  tmp_ptr(i,j,k) = dyn_save_3d(i,j,k,n) * totalsumi
                enddo
              enddo
            enddo
          else
            do k=1,lvl_change
!$omp parallel do private(i,j)
              do j=1,dim2
                do i=1,dim1
                  tmp_ptr(i,j,k) = dyn_save_3d(i,j,k,n) * totalsumi
                enddo
              enddo
            enddo
            do k=lvl_change+1,dim3
!$omp parallel do private(i,j)
              do j=1,dim2
                do i=1,dim1
                  tmp_ptr(i,j,k) = dyn_save_3d(i,j,k,n)
                enddo
              enddo
            enddo
          endif
!
!jun testing
!       CALL ESMF_FieldGet(FIELD, localDe=0, farrayPtr=tmp_ptr1, rc=rc) 
!       do i=1,dim3
!         write(2000+me,*)'field is updated lvl=',i,' data=',maxval(tmp_ptr1(1:dim1,1:dim2,i)),minval(tmp_ptr1(1:dim1,1:dim2,i))
!       enddo
!
        enddo
      endif

      if(allocated(dyn_name))        deallocate(dyn_name)
      if(allocated(dyn_dim))         deallocate(dyn_dim)
      if(allocated(dyn_save_2d))     deallocate(dyn_save_2d)
      if(allocated(dyn_save_3d))     deallocate(dyn_save_3d)
      if(allocated(dyn_items_ord2d)) deallocate(dyn_items_ord2d)
      if(allocated(dyn_items_ord3d)) deallocate(dyn_items_ord3d)

      end subroutine digital_filter_dyn_average_gfs

! ---------------------------------------------------------------
! subroutines for physics
! ---------------------------------------------------------------
      subroutine digital_filter_phy_init_gfs(phy_state)
!
      implicit none
      type(esmf_state), intent(in) :: phy_state
!
      integer rc

      call esmf_stateget(state=phy_state, itemcount=phy_items, rc=rc)
      allocate(phy_name(phy_items))
      call esmf_stateget(state=phy_state, itemnamelist=phy_name, rc=rc)
!       print *,'dfi phy init, phy_items=',phy_items,'names=',phy_name(1:phy_items)

      phy_state_save = esmf_statecreate(name  ="digital filter phy"         &
                                       ,stateintent=ESMF_STATEINTENT_UNSPECIFIED &
                                       ,rc         =rc)
!
      end subroutine digital_filter_phy_init_gfs

! ---------------------------------------------------------------
      subroutine digital_filter_phy_save_gfs(phy_state)
!
      implicit none
      type(esmf_state), intent(in) :: phy_state
!
!jw   TYPE(ESMF_Array)             :: tmp_array
      TYPE(ESMF_Field)             :: tmp_field
      TYPE(ESMF_FieldBUNDLE)       :: tmp_bundle
      TYPE(ESMF_STATE)             :: tmp_state
      type(ESMF_StateItem_Flag)    :: itemtype

      integer                      :: n, rc
!
      do n=1,phy_items

        CALL ESMF_StateGet(phy_state, phy_name(n),itemtype,rc=rc)
!        print *,'restor data,phy_name=',phy_name(n),'itemtype=',itemtype,'rc=',rc

        if(itemtype == ESMF_STATEITEM_FIELDBUNDLE) then
          CALL ESMF_StateGet(phy_state, phy_name(n), tmp_bundle, rc=rc)
          CALL ESMF_StateAddReplace(phy_state_save, (/tmp_bundle/), rc=rc)
!        print *,'save bundle data,phy_name=',phy_name(n),'rc=',rc
        else if(itemtype == ESMF_STATEITEM_FIELD) then
          CALL ESMF_StateGet(phy_state, phy_name(n), tmp_field, rc=rc)
          CALL ESMF_StateAddReplace(phy_state_save, (/tmp_field/), rc=rc)
!        print *,'save field data,phy_name=',phy_name(n),'rc=',rc
        else if(itemtype == ESMF_STATEITEM_STATE) then
          CALL ESMF_StateGet(phy_state, phy_name(n), tmp_state, rc=rc)
          CALL ESMF_StateAddReplace(phy_state_save, (/tmp_state/), rc=rc)
!        print *,'save state data,phy_name=',phy_name(n),'rc=',rc
        endif

      enddo
      end subroutine digital_filter_phy_save_gfs

! ---------------------------------------------------------------
      subroutine digital_filter_phy_restore_gfs(phy_state)
!
      implicit none
      type(esmf_state), intent(inout) :: phy_state
!
      TYPE(ESMF_field)                :: tmp_field
      TYPE(ESMF_FieldBundle)          :: tmp_bundle
      TYPE(ESMF_STATE)                :: tmp_state
      type(ESMF_StateItem_Flag)       :: itemtype
      integer                         :: n, rc
!
      do n=1,phy_items

        CALL ESMF_StateGet(phy_state_save, phy_name(n),itemtype,rc=rc)
!        print *,'restor data,phy_name=',phy_name(n),'itemtype=',itemtype,'rc=',rc

        if(itemtype == ESMF_STATEITEM_FIELDBUNDLE) then
          CALL ESMF_StateGet(phy_state_save, phy_name(n), tmp_bundle, rc=rc)
          CALL ESMF_StateAddReplace(phy_state, (/tmp_bundle/), rc=rc)
!          print *,'restor bundle, ',trim(phy_name(n)),'rc=',rc
        else if(itemtype == ESMF_STATEITEM_FIELD) then
          CALL ESMF_StateGet(phy_state_save, phy_name(n), tmp_field, rc = rc)
          CALL ESMF_StateAddReplace(phy_state, (/tmp_field/), rc=rc)
!          print *,'restor field, ',trim(phy_name(n)),'rc=',rc
        else if(itemtype == ESMF_STATEITEM_STATE) then
          CALL ESMF_StateGet(phy_state_save, phy_name(n), tmp_state, rc=rc)
          CALL ESMF_StateAddReplace(phy_state, (/tmp_state/), rc=rc)
!          print *,'restor state, ',trim(phy_name(n)),'rc=',rc
        endif
      enddo

      call esmf_statedestroy(phy_state_save,rc=rc)
      if (allocated(phy_name)) deallocate(phy_name)

      end subroutine digital_filter_phy_restore_gfs


      end module module_digital_filter_gfs
