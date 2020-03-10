#include "wam_defs.h"

module module_CPLFIELDS

  !-----------------------------------------------------------------------------
  ! ATM Coupling Fields: export and import
  !
  !-----------------------------------------------------------------------------

  use ESMF
  use NUOPC
  use netcdf

  implicit none
  logical, public :: wam_ipe_cpl_rst_input, wam_ipe_cpl_rst_output 
  
  private
  
  integer, parameter :: MAXNAMELEN = 128

  ! private internal state to keep instance data
  integer, dimension(:,:), allocatable :: localNodeToIndexMap
  integer :: localnodes

  real(kind=ESMF_KIND_R8),parameter :: Rearth=6376000.  ! copied from atmos/share/module_CONSTANTS.F90

  ! Regular (non-reduced) Gaussian Grid ---------------
  public            :: gauss2d, wam2dmesh, wamlevels
  type(ESMF_Grid)   :: gauss2d
  type(ESMF_Mesh)   :: wam2dmesh
  integer           :: wamlevels
  
  ! Export Fields ----------------------------------------
  integer, public, parameter :: NexportFields = 56
  type(ESMF_Field), public   :: exportFields(NexportFields)
  character(len=40), public, parameter :: exportFieldsList(NexportFields) = (/ &
      "mean_zonal_moment_flx                  ", &
      "mean_merid_moment_flx                  ", &
      "mean_sensi_heat_flx                    ", &
      "mean_laten_heat_flx                    ", &
      "mean_down_lw_flx                       ", &
      "mean_down_sw_flx                       ", &
      "mean_prec_rate                         ", &
      "inst_zonal_moment_flx                  ", &
      "inst_merid_moment_flx                  ", &
      "inst_sensi_heat_flx                    ", &
      "inst_laten_heat_flx                    ", &
      "inst_down_lw_flx                       ", &
      "inst_down_sw_flx                       ", &
      "inst_temp_height2m                     ", &
      "inst_spec_humid_height2m               ", &
      "inst_zonal_wind_height10m              ", &
      "inst_merid_wind_height10m              ", &
      "inst_temp_height_surface               ", &
      "inst_pres_height_surface               ", &
      "inst_surface_height                    ", &
      "mean_net_lw_flx                        ", &
      "mean_net_sw_flx                        ", &
      "inst_net_lw_flx                        ", &
      "inst_net_sw_flx                        ", &
      "mean_down_sw_ir_dir_flx                ", &
      "mean_down_sw_ir_dif_flx                ", &
      "mean_down_sw_vis_dir_flx               ", &
      "mean_down_sw_vis_dif_flx               ", &
      "inst_down_sw_ir_dir_flx                ", &
      "inst_down_sw_ir_dif_flx                ", &
      "inst_down_sw_vis_dir_flx               ", &
      "inst_down_sw_vis_dif_flx               ", &
      "mean_net_sw_ir_dir_flx                 ", &
      "mean_net_sw_ir_dif_flx                 ", &
      "mean_net_sw_vis_dir_flx                ", &
      "mean_net_sw_vis_dif_flx                ", &
      "inst_net_sw_ir_dir_flx                 ", &
      "inst_net_sw_ir_dif_flx                 ", &
      "inst_net_sw_vis_dir_flx                ", &
      "inst_net_sw_vis_dif_flx                ", &
!     "inst_ir_dir_albedo                     ", &
!     "inst_ir_dif_albedo                     ", &
!     "inst_vis_dir_albedo                    ", &
!     "inst_vis_dif_albedo                    ", &
      "inst_land_sea_mask                     ", &
      "inst_temp_height_lowest                ", &
      "inst_spec_humid_height_lowest          ", &
      "inst_zonal_wind_height_lowest          ", &
      "inst_merid_wind_height_lowest          ", &
      "inst_pres_height_lowest                ", &
      "inst_height_lowest                     ", &
      "mean_fprec_rate                        ", &
      "northward_wind_neutral                 ", &
      "eastward_wind_neutral                  ", &
      "upward_wind_neutral                    ", &
      "temp_neutral                           ", &
      "O_Density                              ", &
      "O2_Density                             ", &
      "N2_Density                             ", &
      "height                                 "  &
  /)

  ! Import Fields ----------------------------------------
  integer, public, parameter :: NimportFields = 16
  type(ESMF_Field), public   :: importFields(NimportFields)
  logical, public            :: importFieldsValid(NimportFields)
  character(len=40), public, parameter :: importFieldsList(NimportFields) = (/ &
      "land_mask                              ", &
      "surface_temperature                    ", &
      "sea_surface_temperature                ", &
      "ice_fraction                           ", &
      "inst_ice_ir_dif_albedo                 ", &
      "inst_ice_ir_dir_albedo                 ", &
      "inst_ice_vis_dif_albedo                ", &
      "inst_ice_vis_dir_albedo                ", &
      "mean_up_lw_flx                         ", &
      "mean_laten_heat_flx                    ", &
      "mean_sensi_heat_flx                    ", &
      "mean_evap_rate                         ", &
      "mean_zonal_moment_flx                  ", &
      "mean_merid_moment_flx                  ", &
      "mean_ice_volume                        ", &
      "mean_snow_volume                       "  /)
  
  ! Utility GSM members ----------------------------------
  public            :: global_lats_ptr
  integer, pointer  :: global_lats_ptr(:)
  public            :: lonsperlat_ptr
  integer, pointer  :: lonsperlat_ptr(:)

  ! Methods
  public fillExportFields
  public queryFieldList
  public setupGauss2d
  public fillWAMFields
  public MeshCreateReducedGaussian
  
  !-----------------------------------------------------------------------------
  contains
  !-----------------------------------------------------------------------------
  
  subroutine fillExportFields(data_a2oi, lonr, latr, rootPet, rc)
    real(kind=ESMF_KIND_R8), target, intent(in) :: data_a2oi(:,:,:)
    integer, intent(in)                         :: lonr, latr, rootPet
    integer, intent(out), optional              :: rc
    
    integer           :: n
    !-----
    ! Fill updated data into the export Fields.
    !-----
    
    if (present(rc)) rc=ESMF_SUCCESS
    
    do n=1, size(exportFields)
      if (ESMF_FieldIsCreated(exportFields(n))) then
        call ESMF_FieldScatter(exportFields(n), data_a2oi(:,:,n), &
          rootPet=rootPet, rc=rc)
        ESMF_ERR_RETURN(rc,rc)
      endif
    enddo

    ESMF_ERR_RETURN(rc,rc)

  end subroutine
  
  !-----------------------------------------------------------------------------

  function MeshCreateReducedGaussian(ipt_lats_node_a, lats_node_a, &
    lonsperlat, global_lats_a, colrad_a, vm, rc) result (mesh)

    integer,                       intent(in)  :: ipt_lats_node_a
    integer,                       intent(in)  :: lats_node_a
    integer,                       intent(in)  :: lonsperlat(:)
    integer,                       intent(in)  :: global_lats_a(:)
    real(ESMF_KIND_R8),            intent(in)  :: colrad_a(:)
    type(ESMF_VM),       optional, intent(in)  :: vm
    integer,             optional, intent(out) :: rc

    type(ESMF_Mesh) :: mesh

    ! -- local variables
    integer :: localrc, stat
    integer :: localPet, petCount
    integer :: latg, latg2, long
    integer :: i, iHemi, ip1, j, jb, jm1, js, k, kk, kp1, l, l1, m, n
    integer :: id, id1, id2, id3, id4
    integer :: lats_nodes, lats_start
    integer :: lnumNodes, numNodes, numElems, numQuads, numTris
    logical :: isVMCreated
    logical, dimension(:),   allocatable :: localNodes
    integer, dimension(:),   allocatable :: nodeIds, nodeOwners
    integer, dimension(:),   allocatable :: globalToLocalIdMap
    integer, dimension(:),   allocatable :: lnodeIds, lnodeOwners
    integer, dimension(:),   allocatable :: elemIds, elemType, elemConn
    integer, dimension(:),   allocatable :: nodeToPetMap
    integer, dimension(:,:), allocatable :: indexToIdMap
    integer, dimension(:,:), allocatable :: indexToLocalIdMap
    real(ESMF_KIND_R8) :: dx, x1, x2, x3, x4
    real(ESMF_KIND_R8), dimension(:),   allocatable :: nodeCoords
    real(ESMF_KIND_R8), dimension(:),   allocatable :: lnodeCoords
    real(ESMF_KIND_R8), dimension(:),   allocatable :: y
    real(ESMF_KIND_R8), dimension(:,:), allocatable :: x
    type(ESMF_VM) :: localVM

    ! -- local parameters
    real(ESMF_KIND_R8), parameter :: rad2deg = &
      57.29577951308232087721_ESMF_KIND_R8

    ! begin
    if (present(rc)) rc = ESMF_SUCCESS

    ! perform sanity check on input array sizes
    latg = size(lonsperlat)
    latg2 = latg / 2

    if (size(global_lats_a) /= latg) then
      call ESMF_LogSetError(ESMF_RC_ARG_SIZE, &
        msg="sizes of global lats and lonsperlats arrays must be the same", &
        line=__LINE__,  &
        file=__FILE__,  &
        rcToReturn=rc)
      return
    end if

    if (size(colrad_a) /= latg2) then
      call ESMF_LogSetError(ESMF_RC_ARG_SIZE, &
        msg="size of colatitude array is inconsistent", &
        line=__LINE__,  &
        file=__FILE__,  &
        rcToReturn=rc)
      return
    end if

    ! get information on parallel environment
    if (present(vm)) then
      isVMCreated = ESMF_VMIsCreated(vm, rc=localrc)
      if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__,  &
        file=__FILE__,  &
        rcToReturn=rc)) &
        return
      if (isVMCreated) then
        ! - use provided VM
        localVM = vm
      else
        call ESMF_LogSetError(ESMF_RC_ARG_BAD, &
          msg="provided VM was not created", &
          line=__LINE__,  &
          file=__FILE__,  &
          rcToReturn=rc)
        return
      end if
    else
      ! - retrieve VM from context
      call ESMF_VMGetCurrent(localVM, rc=localrc)
      if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__,  &
        file=__FILE__,  &
        rcToReturn=rc)) &
        return
    end if

    ! - retrieve grid decomposition across PETs
    call ESMF_VMGet(localVM, localPet=localPet, petCount=petCount, rc=localrc)
    if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__,  &
      file=__FILE__,  &
      rcToReturn=rc)) &
      return

    ! - build node-to-PET map on all PETs
    allocate(nodeToPetMap(2*petCount), stat=stat)
    if (ESMF_LogFoundAllocError(statusToCheck=stat, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__,  &
      file=__FILE__,  &
      rcToReturn=rc)) &

    nodeToPetMap = 0

    call ESMF_VMAllGather(localVM, (/ ipt_lats_node_a, lats_node_a /), &
      nodeToPetMap, 2, rc=localrc)
    if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__,  &
      file=__FILE__,  &
      rcToReturn=rc)) &
      return

    ! pre-compute global mesh coordinates, including 2-point halo region
    long = maxval(lonsperlat)

    n = long + 2
    allocate(x(n,latg), y(latg), stat=stat)
    if (ESMF_LogFoundAllocError(statusToCheck=stat, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__,  &
      file=__FILE__,  &
      rcToReturn=rc)) &
      return

    x = 0._ESMF_KIND_R8
    y = 0._ESMF_KIND_R8
    do j = 1, latg2
      y(j) =  90._ESMF_KIND_R8 - rad2deg * colrad_a(j)
      dx    = 360._ESMF_KIND_R8 / lonsperlat(j)
      do i = 1, lonsperlat(j) + 2
        x(i,j) = (i-1)*dx
      end do
    end do
    x(:,latg:latg2+1:-1) =  x(:,1:latg2)
    y(  latg:latg2+1:-1) = -y(  1:latg2)

    ! compute global nodes
    numNodes = sum(lonsperlat)

    ! - allocate workspace
    allocate(indexToIdMap(long,latg), nodeIds(numNodes), &
      nodeCoords(2*numNodes), nodeOwners(numNodes), stat=stat)
    if (ESMF_LogFoundAllocError(statusToCheck=stat, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__,  &
      file=__FILE__,  &
      rcToReturn=rc)) &
      return

    ! - define global node ids and coordinates
    k = 0
    l = 0
    do j = 1, latg
      do i = 1, lonsperlat(j)
        k = k + 1
        indexToIdMap(i,j) = k
        nodeIds(k) = k
        nodeCoords(l+1) = x(i,j)
        nodeCoords(l+2) = y(j)
        l = l + 2
      end do
    end do

    ! - assign node ownership
    k = 0
    do n = 0, petCount-1
      lats_start = nodeToPetMap(k+1)
      lats_nodes = nodeToPetMap(k+2)
      k = k + 2
      do j = 1, lats_nodes
        l = global_lats_a(lats_start+j-1)
        do i = 1, lonsperlat(l)
          nodeOwners(indexToIdMap(i,l)) = n
        end do
      end do
    end do

    deallocate(nodeToPetMap, stat=stat)
    if (ESMF_LogFoundDeallocError(statusToCheck=stat, &
      msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__,  &
      file=__FILE__,  &
      rcToReturn=rc)) &
      return

    allocate(localNodes(numNodes), &
      indexToLocalIdMap(long,latg), stat=stat)
    if (ESMF_LogFoundAllocError(statusToCheck=stat, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__,  &
      file=__FILE__,  &
      rcToReturn=rc)) &
      return

    ! - assign local nodes based on node ownership first
    do i = 1, numNodes
      id = nodeIds(i)
      localNodes(id) = (nodeOwners(id) == localPet)
    end do

    ! - map longitude and latitude indices to local node id
    ! - this is used to build a local node id to local (i,j)
    !   map to correctly populate fields built on this mesh
    k = 0
    do j = 1, latg
      do i = 1, lonsperlat(j)
        id = indexToIdMap(i,j)
        if (localNodes(id)) then
          k = k + 1
          indexToLocalIdMap(i,j) = k
        end if
      end do
    end do

    ! compute global number and type of elements
    numElems = 0
    numQuads = 0
    numTris = 0
    do iHemi = 0, 1
      jb = (latg-1)*iHemi + 2*(1-iHemi)
      js = 1 - 2*iHemi
      do j = jb, latg2, js
        k = 1
        jm1 = j-js
        l  = lonsperlat(j)
        l1 = lonsperlat(jm1)
        do i = 1, l
          x1 = x(i  ,j)
          x2 = x(i+1,j)
          x3 = x(k+1,jm1)
          x4 = x(k  ,jm1)
          ip1 = mod(i  ,l)+1
          kk  = mod(k-1,l1)+1
          kp1 = mod(k  ,l1)+1
          id1 = indexToIdMap(i  ,j)
          id2 = indexToIdMap(ip1,j)
          id3 = indexToIdMap(kp1,jm1)
          id4 = indexToIdMap(kk ,jm1)
          if (x3 > x2) then
            ! - elements are triangles
            if ((nodeOwners(id1) == localPet) .or. &
                (nodeOwners(id2) == localPet) .or. &
                (nodeOwners(id4) == localPet)) then
              ! - assign all element's nodes to the local PET
              numElems = numElems + 1
              numTris = numTris + 1
              localNodes(id1) = .true.
              localNodes(id2) = .true.
              localNodes(id4) = .true.
            end if
            if (x4 < x2) then
              k = k + 1
              if ((nodeOwners(id2) == localPet) .or. &
                  (nodeOwners(id3) == localPet) .or. &
                  (nodeOwners(id4) == localPet)) then
                numElems = numElems + 1
                numTris = numTris + 1
                localNodes(id2) = .true.
                localNodes(id3) = .true.
                localNodes(id4) = .true.
              end if
            end if
          else
            ! - elements are quadrilaters
            k = k + 1
            if ((nodeOwners(id1) == localPet) .or. &
                (nodeOwners(id2) == localPet) .or. &
                (nodeOwners(id3) == localPet) .or. &
                (nodeOwners(id4) == localPet)) then
              numElems = numElems + 1
              numQuads = numQuads + 1
              ! - assign all element's nodes to the local PET
              localNodes(id1) = .true.
              localNodes(id2) = .true.
              localNodes(id3) = .true.
              localNodes(id4) = .true.
            end if
          end if
        end do
      end do
    end do

    ! now compute local nodes
    lnumNodes = count(localNodes)

    ! - allocate local node arrays and global to local node id map
    allocate(lnodeIds(lnumNodes), lnodeCoords(2*lnumNodes), &
      lnodeOwners(lnumNodes), globalToLocalIdMap(numNodes), stat=stat)
    if (ESMF_LogFoundAllocError(statusToCheck=stat, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__,  &
      file=__FILE__,  &
      rcToReturn=rc)) &
      return

    globalToLocalIdMap = 0

    k = 0
    l = 0
    n = 0
    do i = 1, numNodes
      id = nodeIds(i)
      n = 2*(id-1)
      if (localNodes(id)) then
        k = k + 1
        globalToLocalIdMap(id) = k
        lnodeIds(k) = id
        lnodeOwners(k) = nodeOwners(id)
        lnodeCoords(l+1) = nodeCoords(n+1)
        lnodeCoords(l+2) = nodeCoords(n+2)
        l = l + 2
      end if
    end do

    ! free up memory used by nodes
    deallocate(nodeIds, nodeCoords, stat=stat)
    if (ESMF_LogFoundDeallocError(statusToCheck=stat, &
      msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__,  &
      file=__FILE__,  &
      rcToReturn=rc)) &
      return

    ! create Mesh object
    mesh = ESMF_MeshCreate(parametricDim=2, spatialDim=2, &
      coordSys=ESMF_COORDSYS_SPH_DEG, rc=localrc)
    if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__,  &
      file=__FILE__,  &
      rcToReturn=rc)) &
      return

    ! add local nodes
    call ESMF_MeshAddNodes(mesh, nodeIds=lnodeIds, &
      nodeCoords=lnodeCoords, nodeOwners=lnodeOwners, rc=localrc)
    if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__,  &
      file=__FILE__,  &
      rcToReturn=rc)) &
      return

    ! free up memory used by nodes
    deallocate(lnodeIds, lnodeCoords, lnodeOwners, stat=stat)
    if (ESMF_LogFoundDeallocError(statusToCheck=stat, &
      msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__,  &
      file=__FILE__,  &
      rcToReturn=rc)) &
      return

    ! define local elements
    ! - allocate work arrays for local elements
    allocate(elemIds(numElems), elemType(numElems), &
      elemConn(3*numTris+4*numQuads), stat=stat)
    if (ESMF_LogFoundAllocError(statusToCheck=stat, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__,  &
      file=__FILE__,  &
      rcToReturn=rc)) &
      return

   ! -- define Mesh elementa and connectivity
    m = 0
    n = 0
    id = 0
    do iHemi = 0, 1
      jb = (latg-1)*iHemi + 2*(1-iHemi)
      js = 1 - 2*iHemi
      do j = jb, latg2, js
        k = 1
        jm1 = j-js
        l  = lonsperlat(j)
        l1 = lonsperlat(jm1)
        do i = 1, l
          x1 = x(i  ,j)
          x2 = x(i+1,j)
          x3 = x(k+1,jm1)
          x4 = x(k  ,jm1)
          ip1 = mod(i  ,l)+1
          kk  = mod(k-1,l1)+1
          kp1 = mod(k  ,l1)+1
          id1 = indexToIdMap(i  ,j)
          id2 = indexToIdMap(ip1,j)
          id3 = indexToIdMap(kp1,jm1)
          id4 = indexToIdMap(kk ,jm1)
          if (x3 > x2) then
            ! - create triangles
            id = id + 1
            if ((nodeOwners(id1) == localPet) .or. &
                (nodeOwners(id2) == localPet) .or. &
                (nodeOwners(id4) == localPet)) then
              n = n + 1
              elemIds(n) = id
              elemType(n) = ESMF_MESHELEMTYPE_TRI
              elemConn(m + 1) = id1
              elemConn(m + 2) = id2
              elemConn(m + 3) = id4
              m = m + 3
            end if
            if (x4 < x2) then
              k = k + 1
              id = id + 1
              if ((nodeOwners(id2) == localPet) .or. &
                  (nodeOwners(id3) == localPet) .or. &
                  (nodeOwners(id4) == localPet)) then
                n = n + 1
                elemIds(n) = id
                elemType(n) = ESMF_MESHELEMTYPE_TRI
                elemConn(m + 1) = id4
                elemConn(m + 2) = id2
                elemConn(m + 3) = id3
                m = m + 3
              end if
            end if
          else
            ! - create quadrilaters
            k = k + 1
            id = id + 1
            if ((nodeOwners(id1) == localPet) .or. &
                (nodeOwners(id2) == localPet) .or. &
                (nodeOwners(id3) == localPet) .or. &
                (nodeOwners(id4) == localPet)) then
              n = n + 1
              elemIds(n) = id
              elemType(n) = ESMF_MESHELEMTYPE_QUAD
              elemConn(m + 1) = id1
              elemConn(m + 2) = id2
              elemConn(m + 3) = id3
              elemConn(m + 4) = id4
              m = m + 4
            end if
          end if
        end do
      end do
    end do

    ! - convert local element id to local element id
    do i = 1, size(elemConn)
      elemConn(i) = globalToLocalIdMap(elemConn(i))
    end do

    ! - free up memory used for elements definition
    deallocate(globalToLocalIdMap, indexToIdMap, localNodes, &
      nodeOwners, x, y, stat=stat)
    if (ESMF_LogFoundDeallocError(statusToCheck=stat, &
      msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__,  &
      file=__FILE__,  &
      rcToReturn=rc)) &
      return

    ! -- add elements to Mesh object
    call ESMF_MeshAddElements(mesh, elementIds=elemIds, &
      elementTypes=elemType, elementConn=elemConn, rc=localrc)
    if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__,  &
      file=__FILE__,  &
      rcToReturn=rc)) &
      return

    ! - free up memory used for local element arrays
    deallocate(elemIds, elemType, elemConn, stat=stat)
    if (ESMF_LogFoundDeallocError(statusToCheck=stat, &
      msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__,  &
      file=__FILE__,  &
      rcToReturn=rc)) &
      return

    ! create id map required to correctly populate fields on this mesh
    ! map local node ids to local longitude and latitude indices
    allocate(localNodeToIndexMap(lnumNodes,2), stat=stat)
    if (ESMF_LogFoundAllocError(statusToCheck=stat, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__,  &
      file=__FILE__,  &
      rcToReturn=rc)) &
      return

    localNodeToIndexMap = 0

    do j = 1, lats_node_a
      l = global_lats_a(ipt_lats_node_a+j-1)
      do i = 1, lonsperlat(l)
        id = indexToLocalIdMap(i,l)
        localNodeToIndexMap(id,1) = i
        localNodeToIndexMap(id,2) = j
      end do
    end do

    deallocate(indexToLocalIdMap, stat=stat)
    if (ESMF_LogFoundDeallocError(statusToCheck=stat, &
      msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__,  &
      file=__FILE__,  &
      rcToReturn=rc)) &
      return

  end function MeshCreateReducedGaussian

  !-----------------------------------------------------------------------------

  subroutine setupGauss2d(lonr, latr, pi, colrad_a, lats_node_a, &
    global_lats_a, lonsperlat, rc)
    integer, intent(in)                         :: lonr, latr 
    real(kind=ESMF_KIND_R8), intent(in)         :: pi, colrad_a(:)
    integer, intent(in)                         :: lats_node_a
    integer, intent(in), target                 :: global_lats_a(:)
    integer, intent(in), target                 :: lonsperlat(:)
    integer, intent(out), optional              :: rc
    
    !-----
    ! Create a regular (non-reduced) Gaussian Grid according to NEMS parameters.
    !-----

    integer                                     :: i, j
    real(kind=ESMF_KIND_R8), pointer            :: lonPtr(:,:), latPtr(:,:)
    real(kind=ESMF_KIND_R8), pointer            :: lonCorPtr(:,:), latCorPtr(:,:)
    real(kind=ESMF_KIND_R8), pointer            :: areaPtr(:,:)
    integer(kind=ESMF_KIND_I4), pointer         :: maskPtr(:,:)
    real(kind=ESMF_KIND_R8)                     :: latCorjp1
    character(len=256)                          :: tmpstr
    type(ESMF_VM)                               :: vm
    integer                                     :: petCount
    integer, allocatable                        :: latCounts(:)

    if (present(rc)) rc=ESMF_SUCCESS
    
    call ESMF_VMGetCurrent(vm, rc=rc)
    ESMF_ERR_RETURN(rc,rc)
    
    call ESMF_VMGet(vm, petCount=petCount, rc=rc)
    ESMF_ERR_RETURN(rc,rc)
    
    allocate(latCounts(petCount))

    ! gather the latitude counts on all PETs as an array
    call ESMF_VMAllGather(vm, (/lats_node_a/), latCounts, count=1, rc=rc)
    ESMF_ERR_RETURN(rc,rc)

    ! Create a global spherical grid that is decomposed along latitude dim
    ! the same way that GSM decomposes the Grid.
    gauss2d = ESMF_GridCreate1PeriDim(minIndex=(/1,1/), &
      countsPerDEDim1=(/lonr/),&! 1 DE along "i", i.e. longitude, w/ all longit.
      countsPerDEDim2=latCounts,&! petCount DEs along "j", i.e. latitude w/ cnts
      indexflag=ESMF_INDEX_GLOBAL, coordSys=ESMF_COORDSYS_SPH_DEG, rc=rc)
    ESMF_ERR_RETURN(rc,rc)
    
    ! add coordinates    
    call ESMF_GridAddCoord(gauss2d, staggerLoc=ESMF_STAGGERLOC_CENTER, rc=rc)
    ESMF_ERR_RETURN(rc,rc)

    call ESMF_GridAddCoord(gauss2d, staggerLoc=ESMF_STAGGERLOC_CORNER, rc=rc)
    ESMF_ERR_RETURN(rc,rc)
    
    call ESMF_GridAddItem(gauss2d, itemFlag=ESMF_GRIDITEM_MASK, itemTypeKind=ESMF_TYPEKIND_I4, staggerLoc=ESMF_STAGGERLOC_CENTER, rc=rc)
    ESMF_ERR_RETURN(rc,rc)

    call ESMF_GridAddItem(gauss2d, itemFlag=ESMF_GRIDITEM_AREA, itemTypeKind=ESMF_TYPEKIND_R8, staggerLoc=ESMF_STAGGERLOC_CENTER, rc=rc)
    ESMF_ERR_RETURN(rc,rc)

    !--- CORNERS ---

    call ESMF_GridGetCoord(gauss2d, coordDim=1, staggerloc=ESMF_STAGGERLOC_CORNER, farrayPtr=lonCorPtr, rc=rc)
    ESMF_ERR_RETURN(rc,rc)

    call ESMF_GridGetCoord(gauss2d, coordDim=2, staggerloc=ESMF_STAGGERLOC_CORNER, farrayPtr=latCorPtr, rc=rc)
    ESMF_ERR_RETURN(rc,rc)
    
    write(tmpstr,'(a,4i8)') 'gsm gauss2d corner ',lbound(lonCorPtr,1),ubound(lonCorPtr,1),lbound(lonCorPtr,2),ubound(lonCorPtr,2)
    call ESMF_LogWrite(trim(tmpstr), ESMF_LOGMSG_INFO, line=__LINE__, file=__FILE__, rc=rc)

    ! fill coordinate arrays the same way GSM sets up a non-reduced Gaussian
    do j=lbound(lonCorPtr,2),ubound(lonCorPtr,2)
    do i=lbound(lonCorPtr,1),ubound(lonCorPtr,1)
      lonCorPtr(i,j) = 360./real(lonr) * (real(i)-1.5)
      if (j == 1) then
        latCorPtr(i,j) = 90.
      elseif (j == latr+1) then
        latCorPtr(i,j) = -90.
      elseif (j == latr/2+1) then
        latCorPtr(i,j) = 0.
      elseif (j < latr/2+1) then
        latCorPtr(i,j) = 90. - 180./pi * 0.5*(colrad_a(j)+colrad_a(j-1))
      else
        latCorPtr(i,j) = 180./pi * 0.5*(colrad_a(latr+1-j)+colrad_a(latr+1-j+1)) - 90.
      endif
    enddo
    enddo

    !--- CENTERS ---

    call ESMF_GridGetCoord(gauss2d, coordDim=1, farrayPtr=lonPtr, rc=rc)
    ESMF_ERR_RETURN(rc,rc)

    call ESMF_GridGetCoord(gauss2d, coordDim=2, farrayPtr=latPtr, rc=rc)
    ESMF_ERR_RETURN(rc,rc)
    
    call ESMF_GridGetItem(gauss2d, itemflag=ESMF_GRIDITEM_MASK, farrayPtr=maskPtr, rc=rc)
    ESMF_ERR_RETURN(rc,rc)

    call ESMF_GridGetItem(gauss2d, itemflag=ESMF_GRIDITEM_AREA, farrayPtr=areaPtr, rc=rc)
    ESMF_ERR_RETURN(rc,rc)

    write(tmpstr,'(a,4i8)') 'gsm gauss2d center ',lbound(lonPtr,1),ubound(lonPtr,1),lbound(lonPtr,2),ubound(lonPtr,2)
    call ESMF_LogWrite(trim(tmpstr), ESMF_LOGMSG_INFO, line=__LINE__, file=__FILE__, rc=rc)

    ! fill coordinate arrays the same way GSM sets up a non-reduced Gaussian
    ! tcraig, this is not correct, but is a starting point.
    do j=lbound(lonPtr,2),ubound(lonPtr,2)
      if (j+1 == 1) then
        latCorjp1 = 90.
      elseif (j+1 == latr+1) then
        latCorjp1 = -90.
      elseif (j+1 == latr/2+1) then
        latCorjp1 = 0.
      elseif (j+1 < latr/2+1) then
        latCorjp1 = 90. - 180./pi * 0.5*(colrad_a(j)+colrad_a(j+1))
      else
        latCorjp1 = 180./pi * 0.5*(colrad_a(latr+1-j)+colrad_a(latr+1-j-1)) - 90.
      endif
    do i=lbound(lonPtr,1),ubound(lonPtr,1)
      lonPtr(i,j) = 360./real(lonr) * (i-1)
      if (j <= latr/2) then
        latPtr(i,j) = 90. - 180./pi * colrad_a(j)
      else
        latPtr(i,j) = 180./pi * colrad_a(latr+1-j) - 90.
      endif
      maskPtr(i,j) = 1
!      areaPtr(i,j) = abs(2.*pi/real(lonr) * cos(latPtr(i,j)*pi/180.) * pi/real(latr) * Rearth * Rearth)
      areaPtr(i,j) = abs(2.*pi/real(lonr) * cos(latPtr(i,j)*pi/180.) * pi/180.*(latCorjp1-latCorPtr(i,j)) * Rearth * Rearth)
    enddo
    enddo
    
    ! store GSM members for easier access
    global_lats_ptr => global_lats_a
    lonsperlat_ptr => lonsperlat
    
    deallocate(latCounts)

  end subroutine

  ! Create analytical fields for the 2D WAM built on a ESMF_Mesh
  subroutine fillWAMFields(uug, vvg, wwg, ttg, zzg, n2g, rqg, ipt_lats_node_a, global_lats_a, rc)
    
    real(ESMF_KIND_R8), target :: uug(:,:,:)
    real(ESMF_KIND_R8), target :: vvg(:,:,:)
    real(ESMF_KIND_R8), target :: wwg(:,:,:)
    real(ESMF_KIND_R8), target :: ttg(:,:,:)
    real(ESMF_KIND_R8), target :: zzg(:,:,:)
    real(ESMF_KIND_R8), target :: n2g(:,:,:)
    real(ESMF_KIND_R8), intent(in) :: rqg(:,:,:)
    integer, optional :: ipt_lats_node_a
    integer, optional, intent(in) :: global_lats_a(:)
    integer, optional :: rc

    real(ESMF_KIND_R8), pointer   :: fptr(:,:), infptr(:,:,:)
    integer :: i, j, n, k, kk, levels
    integer :: PetNo, PetCnt
    type(ESMF_VM) :: vm
    character(len=128):: fieldName
    character(len=128):: fileName
    real(ESMF_KIND_R8), pointer   :: varbuf(:,:,:)
    integer                         :: nc, varid, status
    integer                         :: start3(3), count3(3)  

    if (present(rc)) rc=ESMF_SUCCESS

    !------------------------------------------------------------------------
    ! get global vm information
    !
    call ESMF_VMGetCurrent(vm, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    ! set up local pet info
    call ESMF_VMGet(vm, localPet=PetNo, petCount=PetCnt, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    do n=1, NexportFields
       if (ESMF_FieldIsCreated(exportFields(n))) then
         call ESMF_FieldGet(exportFields(n), name=fieldName, rc=rc)
         if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, &
            file=__FILE__)) &
            return  ! bail out
         call ESMF_FieldGet(exportFields(n), farrayPtr=fptr, rc=rc)
         if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, &
            file=__FILE__)) &
            return  ! bail out
         print *, 'bound wwg', ubound(wwg), lbound(wwg)
         print *, trim(fieldName), ' wwg size and field size', size(wwg,1),size(wwg,2),&
	       size(wwg,3), size(fptr,1),size(fptr,2)
        
         levels = size(fptr,2)
         if (trim(fieldName) == "O_Density" .or. trim(fieldName)=="O2_Density") then
            if (trim(fieldName) == "O_Density") then
               kk = 3
            else
               kk = 4
            endif
            allocate(varbuf(size(rqg,1),size(rqg,2),levels))
            do i=1,size(rqg,1)
              do j=1,size(rqg,2)
                 do k=1,levels
                   varbuf(i,j,k)=rqg(i,j,levels*kk+k)
                 enddo
              enddo
            enddo
            infptr=>varbuf
         elseif (trim(fieldName) == "northward_wind_neutral") then
            infptr=>vvg
         elseif (trim(fieldName) == "eastward_wind_neutral") then
            infptr=>uug
         elseif (trim(fieldName) == "upward_wind_neutral") then
            infptr=>wwg
         elseif (trim(fieldName) == "temp_neutral") then
            infptr=>ttg
         elseif (trim(fieldName) == "N2_Density") then
            infptr=>n2g
         elseif (trim(fieldName) == "height") then
            infptr=>zzg
         endif
         do i=1,size(fptr,1)
           do j=1,levels
             fptr(i,j)=infptr(localNodeToIndexMap(i,1),&
                              localNodeToIndexMap(i,2),j)
           enddo
         enddo
         if (trim(fieldName) == "O_Density" .or. trim(fieldName)=="O2_Density") then
	    deallocate(varbuf)
         endif
	 print *, trim(fieldName), ' min/max ', minval(fptr), maxval(fptr)
      endif
   enddo

  end subroutine

  integer function queryFieldList(fieldlist, fieldname, abortflag, rc)
    ! returns integer index of first found fieldname in fieldlist
    ! by default, will abort if field not found, set abortflag to false 
    !   to turn off the abort.
    ! return value of < 1 means the field was not found
    character(len=*),intent(in) :: fieldlist(:)
    character(len=*),intent(in) :: fieldname
    logical, optional :: abortflag
    integer, optional :: rc

    integer :: n
    logical :: labort

    labort = .true.
    if (present(abortflag)) then
      labort = abortflag
    endif

    queryFieldList = 0
    n = 1
    do while (queryFieldList < 1 .and. n <= size(fieldlist))  
      if (trim(fieldlist(n)) == trim(fieldname)) then
        queryFieldList = n
      else
        n = n + 1
      endif
    enddo

    if (labort .and. queryFieldList < 1) then
      call ESMF_LogWrite('queryFieldList ABORT on fieldname '//trim(fieldname), ESMF_LOGMSG_INFO, line=__LINE__, file=__FILE__, rc=rc)
      CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    endif
  end function queryFieldList

  !-----------------------------------------------------------------------------

end module
