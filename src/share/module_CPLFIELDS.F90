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
  integer, pointer :: lonind(:),latind(:)
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
  public createWAMGrid
  public fillWAMFields
  
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

  ! Create 2D WAM as a ESMF_Mesh with only distgrid (no coordinates)
  subroutine createWAMGrid(long, latg, levs, ipt_lats_node_a, lats_node_a, global_lats_a, &
  	     	           lonsperlat, xlon, xlat, rc)

     integer                        :: long, latg, levs ! grid dimension (192x94x150)
     integer                        :: ipt_lats_node_a  ! starting lat index for the local processor
     integer                        :: lats_node_a      ! number of latitues in the local processor
     integer(ESMF_KIND_I4), target  :: global_lats_a(:) ! array holds the random shuffle order of latitude index
     integer(ESMF_KIND_I4), target  :: lonsperlat(:)    ! number of longitude points per lat 
     real(ESMF_KIND_R8), target     :: xlon(:,:)        ! local longitude array
     real(ESMF_KIND_R8), target     :: xlat(:,:)        ! local latitude array
     integer, optional              :: rc

     integer             :: i, j, k, ind1
     type(ESMF_DistGrid) :: distgrid
     integer(ESMF_KIND_I4), pointer :: indList(:)
     integer             :: ind
     character(len=128):: fileName
     integer :: PetNo, PetCnt
     type(ESMF_VM) :: vm
     integer           :: nc, varid, status
     integer           :: start2(2), count2(2)  
     real(ESMF_KIND_R8) :: rad2deg
     real(ESMF_KIND_R8), allocatable :: lonbuf(:,:), latbuf(:,:)

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

      ! find the total number of nodes in each processor and create local index table
      localnodes=0
      ind1 = 0
      do i=1,lats_node_a
         ind=global_lats_a(ipt_lats_node_a+i-1)
         localnodes=localnodes + lonsperlat(ind)
      enddo

      write(0,*) 'start lat, count and totalnodes: ', ipt_lats_node_a, lats_node_a, localnodes

      ! Create a distgrid using a collapsed 1D index array based on the local row index
      ! WAM's 2D grid has the latitude from North to South.  Need to create an index that 
      ! use the south to north order.
      allocate(indList(localnodes),lonind(localnodes),latind(localnodes))
      k=1
      do i=1,lats_node_a
        ind=global_lats_a(ipt_lats_node_a+i-1)
        do j=1,lonsperlat(ind)
           ind1 = latg-ind+1
           indList(k)=long*(ind1-1)+j
           lonind(k)=j
           latind(k)=i
           k=k+1
        enddo
      enddo
      distgrid = ESMF_DistGridCreate(indList, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, &
            file=__FILE__)) &
            return  ! bail out
      ! Create mesh using the distgrid as the nodaldistgrid,  no elemdistgrid available
      ! just use nodeldistgrid for both
      wam2dmesh = ESMF_MeshCreate(distgrid,distgrid,rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, &
            file=__FILE__)) &
            return  ! bail out
      wamlevels = levs
      deallocate(indList)

      ! Write out xlon and xlat from every PET to wam3dgridnew2.nc
      filename = 'wam3dgridnew2.nc'
      rad2deg = 180.0/3.1415926
      allocate(lonbuf(size(xlon,1),size(xlon,2)), latbuf(size(xlon,1),size(xlon,2)))
      lonbuf = -999.0
      latbuf = -999.0
      do i=1,lats_node_a
        ind=global_lats_a(ipt_lats_node_a+i-1)
        do j=1,lonsperlat(ind)
           lonbuf(j,i)=xlon(j,i)*rad2deg
           latbuf(j,i)=xlat(j,i)*rad2deg
        enddo
      enddo
      do i=0, PetCnt-1
        if (PetNo == i) then
                 status = nf90_open(filename, NF90_WRITE, nc)
                 call CheckNCError(status, filename)
                 status = nf90_inq_varId(nc, 'lons', varid)
                 call CheckNCError(status, 'lons')
                 start2(2)=ipt_lats_node_a
                 start2(1)=1
                 count2(2)=size(xlon,2)
                 count2(1)=size(xlon,1)
	         status = nf90_put_var(nc, varid, lonbuf, & 
		    start2, count2)
                 call CheckNCError(status, 'lons')
                 status = nf90_inq_varId(nc, 'lats', varid)
                 call CheckNCError(status, 'lats')
	         status = nf90_put_var(nc, varid, latbuf, & 
		    start2, count2)
                 call CheckNCError(status, 'lats')
                 status = nf90_close(nc)
                 call CheckNCError(status, filename)
         endif
         call ESMF_VMBarrier(vm)
      enddo  !i=0,PetCnt-1
      deallocate(lonbuf, latbuf)
      ESMF_ERR_RETURN(rc,rc)

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
    integer, save                   :: slice=1
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
	     fptr(i,j)=infptr(lonind(i),latind(i),j)
           enddo
         enddo
	 ! For validation purpose, write the field out
         if (slice==1 .and. present(global_lats_a)) then
            filename = 'wam3dgridnew2.nc'
            ! write shuffle order global_lats_a(:) first
            if (PetNo == 0) then
                 status = nf90_open(filename, NF90_WRITE, nc)
                 call CheckNCError(status, filename)
                 status = nf90_inq_varId(nc, "ShuffleOrder", varid)
                 call CheckNCError(status, "ShuffleOrder")
	         status = nf90_put_var(nc, varid, global_lats_a) 
                 call CheckNCError(status, 'ShuffleOrder')
                 status = nf90_close(nc)
                 call CheckNCError(status, filename)
	    endif    
            ! The variable in the output file is named "heights" because there is 
            ! another variable named "height" already
            if (trim(fieldName) == "height") fieldName = "heights"
            if( present( ipt_lats_node_a ) )then

              do i=0, PetCnt-1
                if (PetNo == i) then
                   status = nf90_open(filename, NF90_WRITE, nc)
                   call CheckNCError(status, filename)
                   status = nf90_inq_varId(nc, trim(fieldName), varid)
                   call CheckNCError(status, trim(fieldName))
                   start3(3)=1
                   start3(2)=ipt_lats_node_a
                   start3(1)=1
                   count3(3)=size(infptr,3)
                   count3(2)=size(infptr,2)
                   count3(1)=size(infptr,1)
  	         status = nf90_put_var(nc, varid, infptr, & 
  		    start3, count3)
                   call CheckNCError(status, trim(fieldName))
                   status = nf90_close(nc)
                   call CheckNCError(status, filename)
                endif
                call ESMF_VMBarrier(vm)
              enddo  !i=0,PetCnt-1
 
            endif
            if (PetNo == 0) then
       	       print *, 'Write out ', trim(fieldName)
            endif
         endif !slice==2
         if (trim(fieldName) == "O_Density" .or. trim(fieldName)=="O2_Density") then
	    deallocate(varbuf)
         endif
	 print *, trim(fieldName), ' min/max ', minval(fptr), maxval(fptr)
      endif
   enddo
   slice = slice+1

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

!------------------------------------------------------------------------------
!
!  check CDF file error code
!
#undef  ESMF_METHOD
#define ESMF_METHOD "CheckNCError"
subroutine CheckNCError (ncStatus, errmsg)

    integer,          intent(in)  :: ncStatus
    character(len=*), intent(in)  :: errmsg

    character(len=256) :: msg
    integer, parameter :: nf90_noerror = 0
    if ( ncStatus .ne. nf90_noerror) then
        write(msg, '("NetCDF Error: ", A, " : ", A)') &
    		trim(errmsg),trim(nf90_strerror(ncStatus))
        call ESMF_LogSetError(ESMF_FAILURE, &
	      msg=msg, &
	      line=__LINE__, &
              file=__FILE__)
        !bail out 
        call ESMF_Finalize(endflag=ESMF_END_ABORT)
    end if
    return
end subroutine CheckNCError

  !-----------------------------------------------------------------------------

end module
