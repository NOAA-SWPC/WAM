!-----------------------------------------------------------------------
!
      MODULE MODULE_DOMAIN_NUOPC_SET
!
!-----------------------------------------------------------------------
!
!***  Create/generate/set various NUOPC-related items for an NMM domain.
!
!-----------------------------------------------------------------------
!
      USE ESMF
      USE NUOPC
!
      USE MODULE_KINDS 
!
      USE MODULE_CONSTANTS,ONLY : A
!
      USE module_NMM_INTERNAL_STATE,ONLY: NMM_INTERNAL_STATE            &
                                         ,WRAP_NMM_INTERNAL_STATE
!
      USE module_DOMAIN_INTERNAL_STATE,ONLY: DOMAIN_INTERNAL_STATE      &
                                            ,WRAP_DOMAIN_INTERNAL_STATE
!
      USE module_SOLVER_INTERNAL_STATE,ONLY: SOLVER_INTERNAL_STATE      &
                                            ,WRAP_SOLVER_INT_STATE
!
      USE module_DOMAIN_TASK_SPECS,ONLY: DOMAIN_TASK_SPECS
!
      USE module_CPLFIELDS,ONLY: exportFieldsList,importFieldsList      &
                                ,nExportFields,nImportFields            &
                                ,queryFieldList
!
!-----------------------------------------------------------------------
!
      IMPLICIT NONE
!
!-----------------------------------------------------------------------
!
      PRIVATE
!
      PUBLIC :: CONNECT_EXPORT_FIELDS                                   &
               ,CONNECT_IMPORT_FIELDS                                   &
               ,DOMAIN_DESCRIPTORS                                      &
               ,I_AM_ROOT                                               &
               ,I_AM_PET                                                &
               ,DUMP_DOMAIN_DESCRIPTOR                                  &
               ,DUMP_DOMAIN_DESCRIPTORS                                 &
               ,NMMB_CreateDomainFields                                 &
               ,NMMB_CreateRouteHandle                                  &
               ,NMMB_GridCreate                                         &
               ,NMMB_GridUpdate                                         &
               ,NMMB_RegridExport                                       &
               ,NMMB_RegridImport                                       &
               ,ROTANGLE_CELLAREA_SEAMASK
!
!-----------------------------------------------------------------------
!
      TYPE(DOMAIN_TASK_SPECS),DIMENSION(:),POINTER,SAVE ::               &
                                                   DOMAIN_DESCRIPTORS      !<-- Object holding basic task info for coupling
!
!-----------------------------------------------------------------------
!
      CONTAINS
!
!-----------------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!-----------------------------------------------------------------------
!
      FUNCTION NMMB_GridCreate(nestingDomainGridcompIndex               &
                              ,DOMAIN_GRID_COMP                         &
                              ,DOMAIN_DESCRIPTORS,TPH0D,TLM0D           &
                              ,RC )
!
!-----------------------------------------------------------------------
! The NMMB component always runs on the super set of processors of all
! the domains including the parent. ESMF supports a way to create a grid
! on the NMMB component that spans the entire processor list but only has
! data allocation on the given domain processor list, more specifically
! only on the computational PETs of the given domain's processor list.
! In ESMF terms, the NMMB grid for a given domain only has domain elements
! (DEs) on the given domain's processor list while the number of DE is 0
! on all the other processors.
!-----------------------------------------------------------------------
!
    ! arguments
    INTEGER, INTENT(IN)                :: nestingDomainGridcompIndex
    TYPE(ESMF_GridComp),INTENT(IN)     :: DOMAIN_GRID_COMP
    REAL(kind=KFPT),INTENT(IN)         :: TLM0D,TPH0D
    TYPE(DOMAIN_TASK_SPECS),POINTER    :: DOMAIN_DESCRIPTORS(:)
    INTEGER, INTENT(OUT)               :: rc

    ! return
    type(ESMF_Grid)                                  :: NMMB_GridCreate

    ! local variables
    !TYPE(NMM_INTERNAL_STATE),POINTER,SAVE :: NMM_INT_STATE               !<-- The NMM component internal state pointer
    TYPE(WRAP_NMM_INTERNAL_STATE) :: WRAP_NMM                            !<-- The F90 wrap of the NMM internal state
    type(ESMF_GRID)    :: PGRID
    type(ESMF_TypeKind_Flag) :: ctk
    integer :: tcount, ldecount, dimcount, slcount, I, J, lpet, plpet
    integer :: iend, istart, jend, jstart
    integer :: elb(2), eub(2), clb(2), cub(2)
    integer(ESMF_KIND_I4),dimension(:,:), pointer :: fptr_mask
    real(kind=kdbl) :: deg2rad,dlm,dph,j_center,lam0,phi0,pi,rad2deg
    real(ESMF_KIND_R8),dimension(:,:), pointer :: fptr_lat,fptr_lon
    real(ESMF_KIND_R8),dimension(:,:), pointer :: fptr_area
    TYPE(DOMAIN_INTERNAL_STATE),POINTER :: DOMAIN_INT_STATE
    TYPE(WRAP_DOMAIN_INTERNAL_STATE) :: WRAP_DOMAIN
    TYPE(SOLVER_INTERNAL_STATE),POINTER :: SOLVER_INT_STATE
    TYPE(WRAP_SOLVER_INT_STATE) :: WRAP_SOLVER
    type(ESMF_Field)            :: sstField, sstField_nmmb
    type(ESMF_State)            :: pgridState

    integer                     :: ID_X, nElem, pGridIndexCount, deCount, DE, offset
    integer                     :: ITS,ITE,JTS,JTE
    integer                     :: INPES, JNPES
    integer                     :: minIndex(2), maxIndex(2)
    integer                     :: lDetoDeMap(0:0)        ! 1LocalDE/DE/PET
    integer, allocatable        :: pGridIndexCountPDE(:,:)
    integer, allocatable        :: IndicesD1(:), IndicesD2(:)
    integer, allocatable        :: deBlockList(:,:,:)
    integer, allocatable        :: pGridArbIndexList(:,:)
    type(ESMF_DistGrid)         :: pdistGrid, pDistGrid_nmmb
    integer, allocatable        :: petMap(:)
    type(ESMF_DELayout)         :: delayout
    type(ESMF_DistGrid)         :: distgrid_nmmb
    type(ESMF_Grid)             :: grid_nmmb
    character(4096)             :: tmpstr

!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
    ! Implementation
    RC = ESMF_SUCCESS

    INPES = DOMAIN_DESCRIPTORS(nestingDomainGridcompIndex)%INPES
    JNPES = DOMAIN_DESCRIPTORS(nestingDomainGridcompIndex)%JNPES

    minIndex = DOMAIN_DESCRIPTORS(nestingDomainGridcompIndex)%INDX_MIN
    maxIndex = DOMAIN_DESCRIPTORS(nestingDomainGridcompIndex)%INDX_MAX

    delayout = ESMF_DElayOutCreate(DOMAIN_DESCRIPTORS(nestingDomainGridcompIndex)%PET_MAP, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    distgrid_nmmb = ESMF_DistGridCreate(minIndex, maxIndex, &
      regDecomp=(/INPES, JNPES/), delayout=delayout, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    !NMMB_GridCreate=ESMF_GRIDCREATE(distgrid_nmmb, coordSys=ESMF_COORDSYS_SPH_RAD, &
    NMMB_GridCreate=ESMF_GRIDCREATE(distgrid_nmmb, &
      indexflag=ESMF_INDEX_GLOBAL, rc=RC)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    ! Attach coordinates to the NMMB_GridCreate
    ! Attach center coordinate to parent Grid
    call ESMF_GridAddCoord(NMMB_GridCreate, staggerLoc=ESMF_STAGGERLOC_CENTER, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    ! Attach corner coordinate to parent Grid
    call ESMF_GridAddCoord(NMMB_GridCreate, staggerLoc=ESMF_STAGGERLOC_CORNER, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    ! Attach cell area to NMMB Grid
    call ESMF_GridAddItem(NMMB_GridCreate, itemflag=ESMF_GRIDITEM_AREA, &
      staggerLoc=ESMF_STAGGERLOC_CENTER, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return ! bail out
 
    ! Attach cell area to NMMB Grid
    call ESMF_GridAddItem(NMMB_GridCreate, itemflag=ESMF_GRIDITEM_MASK, &
      staggerLoc=ESMF_STAGGERLOC_CENTER, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return ! bail out
 
    active: if(DOMAIN_DESCRIPTORS(nestingDomainGridcompIndex)%TASK_ACTIVE) then

      CALL ESMF_GridCompGetInternalState(DOMAIN_GRID_COMP &
                                        ,WRAP_DOMAIN      &
                                        ,RC)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out

      DOMAIN_INT_STATE=>wrap_domain%DOMAIN_INT_STATE

      CALL ESMF_GridCompGetInternalState(domain_int_state%SOLVER_GRID_COMP, &
        WRAP_SOLVER, RC)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out

      SOLVER_INT_STATE=>wrap_solver%INT_STATE

!-----------------------------------------------------------------------
!***  Compute/save the wind rotation angles (grid to earth), grid cell
!***  areas, and sea mask for the given domain.
!-----------------------------------------------------------------------
!
      ITS=solver_int_state%ITS
      ITE=solver_int_state%ITE
      JTS=solver_int_state%JTS
      JTE=solver_int_state%JTE
!
      ID_X=nestingDomainGridcompIndex
!
      CALL ROTANGLE_CELLAREA_SEAMASK(SOLVER_INT_STATE                   &
                                    ,ITS,ITE,JTS,JTE                    &
                                    ,DOMAIN_DESCRIPTORS(ID_X)%ROT_ANGLE &
                                    ,DOMAIN_DESCRIPTORS(ID_X)%CELL_AREA &
                                    ,DOMAIN_DESCRIPTORS(ID_X)%SEA_MASK )

      ! Define cell area
      call ESMF_GridGetItem(NMMB_GridCreate, itemflag=ESMF_GRIDITEM_AREA, &
        staggerLoc=ESMF_STAGGERLOC_CENTER, &
        localDe=0, farrayPtr=fptr_area, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out

      DOMAIN_DESCRIPTORS(nestingDomainGridcompIndex)%NEST_CELL_AREA => fptr_area
      fptr_area = DOMAIN_DESCRIPTORS(nestingDomainGridcompIndex)%CELL_AREA

      ! Define cell mask
      call ESMF_GridGetItem(NMMB_GridCreate, itemflag=ESMF_GRIDITEM_MASK, &
        staggerLoc=ESMF_STAGGERLOC_CENTER, &
        localDe=0, farrayPtr=fptr_mask, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out

      DOMAIN_DESCRIPTORS(nestingDomainGridcompIndex)%NEST_SEA_MASK => fptr_mask
      fptr_mask = DOMAIN_DESCRIPTORS(nestingDomainGridcompIndex)%SEA_MASK

      offset = 0

      ! Define Center longitude
      call ESMF_GridGetCoord(NMMB_GridCreate, coordDim=1, &
        staggerLoc=ESMF_STAGGERLOC_CENTER, &
        localDe=0, farrayPtr=fptr_lon, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out

      ! Save a reference so moving nest can update the coordinate values
      DOMAIN_DESCRIPTORS(nestingDomainGridcompIndex)%NEST_GLON => fptr_lon

      write(tmpstr, *) 'minIndex= ', minIndex
      call ESMF_LogWrite(trim(tmpstr), ESMF_LOGMSG_INFO, rc=rc)
      write(tmpstr, *) 'maxIndex= ', maxIndex
      call ESMF_LogWrite(trim(tmpstr), ESMF_LOGMSG_INFO, rc=rc)
      write(tmpstr, *) 'ITS= ', solver_int_state%ITS
      call ESMF_LogWrite(trim(tmpstr), ESMF_LOGMSG_INFO, rc=rc)
      write(tmpstr, *) 'ITE= ', solver_int_state%ITE
      call ESMF_LogWrite(trim(tmpstr), ESMF_LOGMSG_INFO, rc=rc)
      write(tmpstr, *) 'JTS= ', solver_int_state%JTS
      call ESMF_LogWrite(trim(tmpstr), ESMF_LOGMSG_INFO, rc=rc)
      write(tmpstr, *) 'JTE= ', solver_int_state%JTE
      call ESMF_LogWrite(trim(tmpstr), ESMF_LOGMSG_INFO, rc=rc)
      write(tmpstr, *) 'shape fptr_lon: ', lbound(fptr_lon), ubound(fptr_lon)
      !write(tmpstr, *) 'value fptr_lon: ', fptr_lon
      call ESMF_LogWrite(trim(tmpstr), ESMF_LOGMSG_INFO, rc=rc)
      write(tmpstr, *) 'shape glon: ', shape(solver_int_state%glon)
      !write(tmpstr, *) 'value glon: ', solver_int_state%glon
      call ESMF_LogWrite(trim(tmpstr), ESMF_LOGMSG_INFO, rc=rc)
      write(tmpstr, *) 'shape glat: ', shape(solver_int_state%glat)
      call ESMF_LogWrite(trim(tmpstr), ESMF_LOGMSG_INFO, rc=rc)
      write(tmpstr, *) 'shape vlon: ', shape(solver_int_state%vlon)
      call ESMF_LogWrite(trim(tmpstr), ESMF_LOGMSG_INFO, rc=rc)
      write(tmpstr, *) 'shape vlat: ', shape(solver_int_state%vlat)
      call ESMF_LogWrite(trim(tmpstr), ESMF_LOGMSG_INFO, rc=rc)
      call ESMF_LogFlush()

      pi=acos(-1._kdbl)
      rad2deg=180._kdbl/pi

      do J = solver_int_state%JTS, solver_int_state%JTE
        do I = solver_int_state%ITS, solver_int_state%ITE
          fptr_lon(I,J) = solver_int_state%GLON(I+offset,J+offset)*rad2deg
!     if((i>=1.and.i<=5.and.j>=1.and.j<=2).or.(i>=1.and.i<=2.and.j>=1.and.j<=5))then
!       write(0,24240)i,j,fptr_lon(i,j)
24240   format(' center i=',i3,' j=',i3,' lon=',e12.5)
!     endif
!     if((i==1.or.i==solver_int_state%IDE).and.(j==1.or.j==solver_int_state%JDE))then
!       write(0,44421)i,j,fptr_lon(i,j)
44421   format(' i=',i3,' j=',i3,' center lon=',e13.6)
!     endif
        enddo
      enddo

      ! Define Center latitude
      call ESMF_GridGetCoord(NMMB_GridCreate, coordDim=2, &
        staggerLoc=ESMF_STAGGERLOC_CENTER, &
        localDe=0, farrayPtr=fptr_lat, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out

      DOMAIN_DESCRIPTORS(nestingDomainGridcompIndex)%NEST_GLAT => fptr_lat

      do J = solver_int_state%JTS, solver_int_state%JTE
        do I = solver_int_state%ITS, solver_int_state%ITE
          fptr_lat(I,J) = solver_int_state%GLAT(I+offset,J+offset)*rad2deg
!     if((i>=1.and.i<=5.and.j>=1.and.j<=2).or.(i>=1.and.i<=2.and.j>=1.and.j<=5))then
!       write(0,24241)i,j,fptr_lat(i,j)
24241   format(' center i=',i3,' j=',i3,' lat=',e12.5)
!     endif
!     if((i==1.or.i==solver_int_state%IDE).and.(j==1.or.j==solver_int_state%JDE))then
!       write(0,44422)i,j,fptr_lat(i,j)
44422   format(' i=',i3,' j=',i3,' center lat=',e13.6)
!     endif
        enddo
      enddo

      ! Define Corner longitude
      call ESMF_GridGetCoord(NMMB_GridCreate, coordDim=1, &
        staggerLoc=ESMF_STAGGERLOC_CORNER, &
        localDe=0, farrayPtr=fptr_lon, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out

      DOMAIN_DESCRIPTORS(nestingDomainGridcompIndex)%NEST_VLON => fptr_lon

!-----------------------------------------------------------------------
!***  On the B grid V points lie to the northeast of H points with the
!***  same indices.  So first fill all the lower left grid cell corners
!***  with V-point coordinates then fill in those corner values along
!***  the south and east domain edges by computing the correct
!***  coordinates.
!-----------------------------------------------------------------------

      offset=-1

      istart=max(2,solver_int_state%ITS)
      iend=solver_int_state%ITE
      if(solver_int_state%ITE==solver_int_state%IDE)THEN
        iend=solver_int_state%IDE+1
      endif

      jstart=max(2,solver_int_state%JTS)
      jend=solver_int_state%JTE
      if(solver_int_state%JTE==solver_int_state%JDE)THEN
        jend=solver_int_state%JDE+1
      endif

      do J = jstart, jend
        do I = istart, iend
          fptr_lon(I,J) = solver_int_state%VLON(I+offset,J+offset)
!     if(i>=solver_int_state%ide-1.or.j>=solver_int_state%jde-1)then
!       write(0,45451)i,j,i+offset,j+offset,solver_int_state%VLAT(i+offset,J+offset)*rad2deg &
!                                          ,solver_int_state%VLON(I+offset,J+offset)*rad2deg
45451   format(' NMMB_GridCreate i=',i3,' j=',i3,' i+offset=',i3,' j+offset=',i3,' vlat=',e12.5,' vlon(offset)=',e12.5)
!     endif
        enddo
      enddo

      ! Define Corner latitude
      call ESMF_GridGetCoord(NMMB_GridCreate, coordDim=2, &
        staggerLoc=ESMF_STAGGERLOC_CORNER, &
        localDe=0, farrayPtr=fptr_lat, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out

      DOMAIN_DESCRIPTORS(nestingDomainGridcompIndex)%NEST_VLAT => fptr_lat

      do J = jstart, jend
        do I = istart, iend
          fptr_lat(I,J) = solver_int_state%VLAT(I+offset,J+offset)
        enddo
      enddo

!***  The dimensions of the arrays holding the cell corners depends on
!***  whether the cell lies on the domain boundary.  For cells on the
!***  north boundary the J dimension is JTE+1 while for cells on the
!***  east boundary the I dimension is ITE+1.

      call cell_corners_on_bndry(fptr_lat,fptr_lon,iend,jend                 &
                                ,solver_int_state%GLAT,solver_int_state%GLON &
                                ,solver_int_state%ITS,solver_int_state%ITE   &
                                ,solver_int_state%JTS,solver_int_state%JTE   &
                                ,solver_int_state%IMS,solver_int_state%IME   &
                                ,solver_int_state%JMS,solver_int_state%JME   &
                                ,solver_int_state%IDS,solver_int_state%IDE   &
                                ,solver_int_state%JDS,solver_int_state%JDE )

!***  Convert corner coordinates from radians to degrees for ESMF regrid.

      do J = solver_int_state%JTS, jend
      do I = solver_int_state%ITS, iend
        fptr_lat(I,J) = fptr_lat(I,J) * rad2deg
        fptr_lon(I,J) = fptr_lon(I,J) * rad2deg
!     if((i>=1.and.i<=5.and.j>=1.and.j<=2).or.(i>=1.and.i<=2.and.j>=1.and.j<=5))then
!       write(0,24242)i,j,fptr_lat(i,j),fptr_lon(i,j)
24242   format(' corner i=',i3,' j=',i3,' lat=',e12.5,' lon=',e12.5)
!     endif
!     if((i<=2.or.i>=solver_int_state%IDE).and.(j<=2.or.j>=solver_int_state%JDE))then
!       write(0,44401)i,j,fptr_lat(I,J),fptr_lon(I,J)
44401   format(' i=',i3,' j=',i3,' corner glat=',e13.6,' corner glon=',e13.6)
!     endif
      enddo
      enddo

!-----------------------------------------------------------------------

    endif active

!
!-----------------------------------------------------------------------
!
    contains
!
!-----------------------------------------------------------------------
!
      subroutine cell_corners_on_bndry(vlat,vlon,ilim,jlim              &
                                      ,glat,glon                        &
                                      ,its,ite,jts,jte                  &
                                      ,ims,ime,jms,jme                  &
                                      ,ids,ide,jds,jde)

!-----------------------------------------------------------------------
!***  Compute the geographic lat/lon of the lower left corners of
!***  NMMB grid cells that lie on the south and west boundary of
!***  the domain.
!-----------------------------------------------------------------------

      implicit none

      integer(kind=KINT),intent(in) :: ilim,jlim                        &  !<-- Upper limits of array of corners on this task
                                      ,ids,ide,jds,jde                  &  !<-- Domain limits
                                      ,ims,ime,jms,jme                  &  !<-- Memory limits of MPI task subdomain
                                      ,its,ite,jts,jte                     !<-- Integration limits of task subdomain

      real(ESMF_KIND_R8),dimension(its:ilim,jts:jlim),intent(inout) ::  &
                                                               vlat,vlon   !<-- Geographic lat/lon (rad) of V points on NMMB grid
      real(kind=KFPT),dimension(ims:ime,jms:jme),intent(in) :: glat,glon   !<-- Geographic lat/lon (rad) of H points on NMMB grid
      integer(kind=KINT) :: i,j
      real(kind=KDBL) :: arg1,arg2,deg2rad,dlm,dph,fctr                 &
                        ,pi,tlat,tlm0,tlon,tph0,x,y,z
      real(kind=KDBL),dimension(its:its+1,jts:jts+1) :: hlat,hlon

      pi=acos(-1._kdbl)
      deg2rad=pi/180._kdbl
      tph0=tph0d*deg2rad
      tlm0=tlm0d*deg2rad
!     write(0,40001)tph0d,tlm0d
40001 format(' corners tph0d=',e12.5,' tlm0d=',e12.5)

!***  We need to know the grid increments in terms of rotated lat/lon.

      do j=jts,jts+1
      do i=its,its+1
        x=cos(tph0)*cos(glat(i,j))*cos(glon(i,j)-tlm0)+sin(tph0)*sin(glat(i,j))
        y= cos(glat(i,j))*sin(glon(i,j)-tlm0)
        z=cos(tph0)*sin(glat(i,j))-sin(tph0)*cos(glat(i,j))*cos(glon(i,j)-tlm0)
        hlat(i,j)=atan(z/sqrt(x*x+y*y))                                    !<-- Rotated lat (radians) of SW corner H points
        hlon(i,j)=atan(y/x)                                                !<-- Rotated lon (radians) of SW corner H points
!     if(its==1.and.jts==1)then
!       write(0,32321)i,j,hlat(i,j)/deg2rad,hlon(i,j)/deg2rad
32321   format(' rot h lat i=',i3,' j=',i3,' lat=',e12.5,' lon=',e12.5)
!     endif
      enddo
      enddo

      dph=hlat(its,jts+1)-hlat(its,jts)                                    !<-- The grid increment of rotated latitude (radians)
      dlm=hlon(its+1,jts)-hlon(its,jts)                                    !<-- The grid increment of rotated longitude (radians)
!     if(its==1.and.jts==1)then
!       write(0,32322)dph/deg2rad,dlm/deg2rad
32322   format(' dphd=',e17.10,' dlmd=',e17.10)
!     endif

!***  Using the rotated lat/lon of phantom V points just south and
!***  west of the domain boundary find their geographic lat/lon.
!***  First the south side.

      tlat=hlat(its,jts)-dph*0.5_kdbl                                      !<-- Rotated lat (deg) of phantom V points south of boundary
      tlon=hlon(its,jts)-dlm*1.5_kdbl

      if(jts==1)then                                                       !<-- Select the tasks on the south boundary.
        j=1
        do i=its,ilim
          tlon=tlon+dlm
!     if((i<=2.or.i>=ide))then
!       write(0,55581)i,j,tlat/deg2rad,tlon/deg2rad
55581   format(' corners south i=',i3,' j=',i3,' tlat=',e12.5,' tlon=',e12.5)
!     endif
          arg1=sin(tlat)*cos(tph0)+cos(tlat)*sin(tph0)*cos(tlon)
          vlat(i,j)=asin(arg1)
          arg2=cos(tlat)*cos(tlon)/(cos(vlat(i,j))*cos(tph0))-          &
               tan(vlat(i,j))*tan(tph0)
          if(abs(arg2)>1.)arg2=abs(arg2)/arg2
          fctr=-1._kdbl
          if(tlon>0.)fctr=1._kdbl
          if(tlon>pi)fctr=-1._kdbl
!
          vlon(i,j)=tlm0+fctr*acos(arg2)
!     if((i<=2.or.i>=ide))then
!       write(0,55582)i,j,vlat(i,j)/deg2rad,vlon(i,j)/deg2rad
55582   format(' corners south i=',i3,' j=',i3,' glat=',e12.5,' glon=',e12.5)
!     endif
          if(vlon(i,j)<-pi)vlon(i,j)=vlon(i,j)+pi*2._kdbl
          if(vlon(i,j)> pi)vlon(i,j)=vlon(i,j)-pi*2._kdbl
        enddo
      endif

!***  Then the west side.

      tlon=hlon(its,jts)-dlm*0.5_kdbl                                      !<-- Rotated lon (deg) of phantom V points west of boundary
      tlat=hlat(its,jts)-dph*1.5_kdbl

      if(its==1)then                                                       !<-- Select the tasks on the west boundary.
        i=1
        do j=jts,jlim
          tlat=tlat+dph
!     if((j<=2.or.j>=jde))then
!       write(0,55583)i,j,tlat/deg2rad,tlon/deg2rad
55583   format(' corners west i=',i3,' j=',i3,' tlat=',e12.5,' tlon=',e12.5)
!     endif
          arg1=sin(tlat)*cos(tph0)+cos(tlat)*sin(tph0)*cos(tlon)
          vlat(i,j)=asin(arg1)
          arg2=cos(tlat)*cos(tlon)/(cos(vlat(i,j))*cos(tph0))-          &
               tan(vlat(i,j))*tan(tph0)
          if(abs(arg2)>1.)arg2=abs(arg2)/arg2
          fctr=-1._kdbl
          if(tlon>0.)fctr=1._kdbl
          if(tlon>pi)fctr=-1._kdbl
!
          vlon(i,j)=tlm0+fctr*acos(arg2)
          if(vlon(i,j)<-pi)vlon(i,j)=vlon(i,j)+pi*2._kdbl
          if(vlon(i,j)> pi)vlon(i,j)=vlon(i,j)-pi*2._kdbl
!     if((j<=2.or.j>=jde))then
!       write(0,55584)i,j,vlat(i,j)/deg2rad,vlon(i,j)/deg2rad
55584   format(' corners west i=',i3,' j=',i3,' glat=',e12.5,' glon=',e12.5)
!     endif
        enddo
      endif

!-----------------------------------------------------------------------

      end subroutine cell_corners_on_bndry

!-----------------------------------------------------------------------
!
      END FUNCTION NMMB_GridCreate
!
!-----------------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!-----------------------------------------------------------------------
!
      SUBROUTINE ROTANGLE_CELLAREA_SEAMASK(SOLVER_INT_STATE             &
                                          ,ITS,ITE,JTS,JTE              &
                                          ,ROT_ANGLE                    &
                                          ,CELL_AREA                    &
                                          ,SEA_MASK )
!
!-----------------------------------------------------------------------
!***  Compute/save the wind rotation angles (grid to earth), grid cell
!***  areas, and sea mask for the given domain.
!-----------------------------------------------------------------------
!
!------------------------
!***  Argument variables
!------------------------
!
      INTEGER(kind=KINT) :: ITS,ITE,JTS,JTE                                !<-- Subdomain integration limits for the given MPI task.
!
      TYPE(SOLVER_INTERNAL_STATE),POINTER,INTENT(IN) :: SOLVER_INT_STATE   !<-- The Solver component's internal state.
!
      REAL(kind=KDBL),DIMENSION(ITS:ITE,JTS:JTE),INTENT(OUT) :: ROT_ANGLE &  !<-- Wind rotations angles (radians)
                                                               ,CELL_AREA &  !<-- Grid cell areas (m**2)
                                                               ,SEA_MASK     !<-- Sea mask (1.0->water; 0.0->land)
!
!---------------------
!***  Local variables
!---------------------
!
      INTEGER(kind=KINT) :: I,J
!
      REAL(kind=KDBL) :: ARG,COS_LAM,COS_PHI,DEG2RAD,DLM,DPH,J_CENTER   &
                        ,LAM0,LAMBDA,PHI0,PI,ROT_LAT,SIN_LAM,SIN_PHI,X,Y
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
      PI=ACOS(-1._kdbl)
      DEG2RAD=PI/180._kdbl
      PHI0=solver_int_state%TPH0D*DEG2RAD
      LAM0=solver_int_state%TLM0D*DEG2RAD
!
      DPH=2._kdbl*solver_int_state%SBD*DEG2RAD/                         &  !<-- Increment of cell in rotated latitude (rad)
           (solver_int_state%JDE-solver_int_state%JDS)
      DLM=2._kdbl*solver_int_state%WBD*DEG2RAD/                         &  !<-- Increment of cell in rotated longitude (rad)
           (solver_int_state%IDE-solver_int_state%IDS)
!
      J_CENTER=0.5_kdbl*(solver_int_state%JDS+solver_int_state%JDE)        !<-- J index of this domain's central location
!
!
!-----------------------------------------------------------------------
!
      DO J=JTS,JTE
      DO I=ITS,ITE
!
!-----------------------------------------------------------------------
!***  Compute the rotation angle for winds (grid to earth).
!-----------------------------------------------------------------------
!
        COS_PHI=COS(solver_int_state%GLAT(I,J))
        SIN_PHI=SIN(solver_int_state%GLAT(I,J))
        COS_LAM=COS(solver_int_state%GLON(I,J)-LAM0)
        SIN_LAM=SIN(solver_int_state%GLON(I,J)-LAM0)
!
        X=COS(PHI0)*COS_PHI*COS_LAM+SIN(PHI0)*SIN_PHI
        Y=COS_PHI*SIN_LAM
        LAMBDA=ATAN(Y/X)
!
        ARG=SIN(PHI0)*SIN(LAMBDA)/COS_PHI
        IF(ABS(ARG)>1.)THEN
          ARG=SIGN(1.,ARG)
        ENDIF
!
        ROT_ANGLE(I,J)=ASIN(ARG)
!
!-----------------------------------------------------------------------
!***  Compute the grid cell areas.
!-----------------------------------------------------------------------
!
        ROT_LAT=(J-J_CENTER)*DPH                                           !<-- Rotated latitude (rad) of center of grid cell
!
        CELL_AREA(I,J)=2._kdbl*A**2*DLM*COS(ROT_LAT)*SIN(0.5_kdbl*DPH)     !<-- Grid cell area (m**2)
!
!-----------------------------------------------------------------------
!***  Store the sea masks for each domain (1.0->water; 0.0->land).
!-----------------------------------------------------------------------
!
        SEA_MASK(I,J)=solver_int_state%SM(I,J)
!
!-----------------------------------------------------------------------
!
      ENDDO
      ENDDO
!
!-----------------------------------------------------------------------
!
      END SUBROUTINE ROTANGLE_CELLAREA_SEAMASK
!
!-----------------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!-----------------------------------------------------------------------
!
      SUBROUTINE NMMB_CreateRouteHandle(N, DOMAIN_DESCRIPTORS, RC)
!
      INTEGER, INTENT(IN)                       :: N
      TYPE(DOMAIN_TASK_SPECS), POINTER          :: DOMAIN_DESCRIPTORS(:)
      INTEGER, INTENT(OUT), OPTIONAL            :: RC

      TYPE(ESMF_Field)                          :: parentField, nestField
      INTEGER                                   :: parent_domain_id
      INTEGER                                   :: srcTermProcessing_Value = 0
      CHARACTER(len=2)                          :: msg
      INTEGER                                   :: nthgrid = 1

!-----------------------------------------------------------------------

      if(PRESENT(RC)) RC = ESMF_SUCCESS
      if(N == DOMAIN_DESCRIPTORS(N)%PARENT_DOMAIN_ID) return !<-- Upper parent will not regrid to itself.

      parent_domain_id = DOMAIN_DESCRIPTORS(N)%PARENT_DOMAIN_ID
      parentField = ESMF_FieldCreate(DOMAIN_DESCRIPTORS(parent_domain_id)%GRID, typekind=ESMF_TYPEKIND_R8, RC=RC)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out
      nestField = ESMF_FieldCreate(DOMAIN_DESCRIPTORS(N)%GRID, typekind=ESMF_TYPEKIND_R8, RC=RC)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out

      nthgrid = nthgrid + 1
      if(nthgrid == 10) then
        write(msg, '(I2.2)') N
        call Grid_Write(DOMAIN_DESCRIPTORS(N)%GRID, "NMMB_Nest_"//trim(msg), rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
      endif

      call ESMF_FieldRegridStore(parentField, nestField, &
        routehandle=DOMAIN_DESCRIPTORS(N)%PARENT2SELF, &
        regridmethod=ESMF_REGRIDMETHOD_BILINEAR, &
        polemethod=ESMF_POLEMETHOD_ALLAVG, &
        srcTermProcessing=srcTermProcessing_Value, &      !<-- no partial sum on src side
        ignoreDegenerate=.true., &
        unmappedaction=ESMF_UNMAPPEDACTION_IGNORE, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out
      call ESMF_FieldRegridStore(nestField, parentField, &
        routehandle=DOMAIN_DESCRIPTORS(N)%SELF2PARENT, &
        regridmethod=ESMF_REGRIDMETHOD_BILINEAR, &
        polemethod=ESMF_POLEMETHOD_ALLAVG, &
        srcTermProcessing=srcTermProcessing_Value, &      !<-- no partial sum on src side
        ignoreDegenerate=.true., &
        unmappedaction=ESMF_UNMAPPEDACTION_IGNORE, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out

      END SUBROUTINE NMMB_CreateRouteHandle

!-----------------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!-----------------------------------------------------------------------
!
     SUBROUTINE NMMB_CreateDomainFields(N, DOMAIN_DESCRIPTORS, RC)
!
!-----------------------------------------------------------------------

      INTEGER, INTENT(IN)                       :: N
      TYPE(DOMAIN_TASK_SPECS), POINTER          :: DOMAIN_DESCRIPTORS(:)
      INTEGER, INTENT(OUT), OPTIONAL            :: RC

      TYPE(ESMF_Field)                          :: parentField, nestField
      INTEGER                                   :: fieldIdx

!-----------------------------------------------------------------------

      if(PRESENT(RC)) RC = ESMF_SUCCESS
      ! Allocate for all domains
      if(.not.allocated(DOMAIN_DESCRIPTORS(N)%ImportFieldsList))then
        allocate(DOMAIN_DESCRIPTORS(N)%ImportFieldsList(nImportFields))
        allocate(DOMAIN_DESCRIPTORS(N)%ExportFieldsList(nExportFields))
      endif

      !Print *, 'End of NMMB_CreateDomainFields: ', 'nImportFields = ', nImportFields, &
      !' nExportFields_NMMB = ', nExportFields, ' N = ', N, &
      !' sizeIm = ', size(DOMAIN_DESCRIPTORS(N)%ImportFieldsList), &
      !' sizeEx = ', size(DOMAIN_DESCRIPTORS(N)%ExportFieldsList)

      if(N == DOMAIN_DESCRIPTORS(N)%PARENT_DOMAIN_ID) return !<-- CAP will create the Fields for parent domain

      do fieldIdx = 1,nImportFields
        DOMAIN_DESCRIPTORS(N)%ImportFieldsList(fieldIdx) = ESMF_FieldCreate( &
          DOMAIN_DESCRIPTORS(N)%GRID, typekind=ESMF_TYPEKIND_R8, &
          name=trim(importFieldsList(fieldIdx)), RC=RC)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
      enddo

      do fieldIdx = 1,nExportFields
        call ESMF_LogWrite('NMMB_CreateDomainFields: '//trim(exportFieldsList(fieldIdx)), ESMF_LOGMSG_INFO, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out

        DOMAIN_DESCRIPTORS(N)%ExportFieldsList(fieldIdx) = ESMF_FieldCreate( &
          DOMAIN_DESCRIPTORS(N)%GRID, typekind=ESMF_TYPEKIND_R8, &
          name=trim(exportFieldsList(fieldIdx)), RC=RC)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
      enddo

!-----------------------------------------------------------------------

      END SUBROUTINE NMMB_CreateDomainFields

!-----------------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!-----------------------------------------------------------------------
!
      SUBROUTINE CONNECT_IMPORT_FIELDS(DOMAIN_ID,DOMAIN_DESCRIPTORS     &
                                      ,NMM_GRID_COMP,RC)
!
!-----------------------------------------------------------------------
!***  The ocean's SST has been interpolated onto the NMM-B's grid.
!***  Point at that SST field with a domain-dependent pointer.  The
!***  values will be copied to the Solver's internal state in
!***  DOMAIN_RUN immediately before the call to SOLVER_RUN.
!-----------------------------------------------------------------------
!
      IMPLICIT NONE
!
!-----------------------------------------------------------------------

      INTEGER,INTENT(IN)                    :: DOMAIN_ID
      TYPE(DOMAIN_TASK_SPECS),POINTER       :: DOMAIN_DESCRIPTORS(:)
      TYPE(ESMF_GridComp),INTENT(IN)        :: NMM_GRID_COMP
!
      INTEGER, INTENT(OUT), OPTIONAL        :: RC
!
!---------------------
!***  Local variables
!---------------------
!
      TYPE(NMM_INTERNAL_STATE),POINTER,SAVE :: NMM_INT_STATE               !<-- The NMM component internal state pointer
      TYPE(WRAP_NMM_INTERNAL_STATE) :: WRAP_NMM                            !<-- The F90 wrap of the NMM internal state
      TYPE(DOMAIN_INTERNAL_STATE),POINTER :: DOMAIN_INT_STATE              !<-- The Domain component internal state pointer
      TYPE(WRAP_DOMAIN_INTERNAL_STATE) :: WRAP_DOMAIN                      !<-- The F90 wrap of the Domain internal state
      REAL(ESMF_KIND_R8),POINTER :: SST_FPTR(:,:)
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!     write(0,36611)
36611 format(' enter CONNECT')
!
      IF(PRESENT(RC)) RC= ESMF_SUCCESS
!
!-----------------------------------------------------------------------
!***  Extract the NMM component's internal state.
!-----------------------------------------------------------------------
!
      CALL ESMF_GridCompGetInternalState(NMM_GRID_COMP                  &
                                        ,WRAP_NMM                       &
                                        ,RC )
!
      IF(RC/=0)THEN
        WRITE(0,*)' CONNECT_IMPORT_FIELDS failed to get NMM int state' 
        WRITE(0,*)' RC=',RC
        CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
      ENDIF
!
      NMM_INT_STATE=>wrap_nmm%NMM_INT_STATE
!
!-----------------------------------------------------------------------
!
!      DO DOMAIN_ID=1,NUM_DOMAINS_TOTAL

!     write(0,40021)domain_id,domain_descriptors(domain_id)%task_active
40021 format(' CONNECT domain_id=',i2,' task_active=',l1)
        IF(DOMAIN_DESCRIPTORS(DOMAIN_ID)%TASK_ACTIVE) THEN
          CALL ESMF_GridCompGetInternalState(&
                 nmm_int_state%DOMAIN_GRID_COMP(DOMAIN_ID) &
                ,WRAP_DOMAIN                          &
                ,RC)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, &
            file=__FILE__)) &
            return  ! bail out
!
          DOMAIN_INT_STATE=>wrap_domain%DOMAIN_INT_STATE
!
!         fieldIdx = queryFieldList(importFieldsList, 'sea_surface_temperature', rc=rc)
          CALL ESMF_FieldGet(DOMAIN_DESCRIPTORS(DOMAIN_ID)%importFieldsList(3), farrayPtr=sst_fptr, rc=rc)
!
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, &
            file=__FILE__)) &
            return  ! bail out
!
          domain_int_state%SST_COUPLED => sst_fptr
!     write(0,36612)
36612 format(' CONNECT pointed at sst_fptr')
!
        ENDIF
!
!     write(0,36613)
36613 format(' exit CONNECT')
!-----------------------------------------------------------------------
!
      END SUBROUTINE CONNECT_IMPORT_FIELDS
!
!-----------------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!-----------------------------------------------------------------------
!
      SUBROUTINE CONNECT_EXPORT_FIELDS(DOMAIN_ID, DOMAIN_DESCRIPTORS    &
                                      ,NMM_GRID_COMP, RC)
!
!-----------------------------------------------------------------------
!***  Connect fluxes generated inside the NMM with pointers in the
!***  cap in order to tranfer the fields to the ocean model.
!-----------------------------------------------------------------------
!
      IMPLICIT NONE
!
!-----------------------------------------------------------------------
!
!-----------------------
!*** Argument variables
!-----------------------
!
      INTEGER, INTENT(IN)                   :: DOMAIN_ID
      TYPE(DOMAIN_TASK_SPECS),POINTER       :: DOMAIN_DESCRIPTORS(:)
      TYPE(ESMF_GridComp),INTENT(IN)        :: NMM_GRID_COMP
      INTEGER, INTENT(OUT), OPTIONAL        :: RC
!
!---------------------
!***  Local variables
!---------------------
!
      TYPE(NMM_INTERNAL_STATE),POINTER,SAVE :: NMM_INT_STATE               !<-- The NMM component internal state pointer
      TYPE(WRAP_NMM_INTERNAL_STATE) :: WRAP_NMM                            !<-- The F90 wrap of the NMM internal state
      TYPE(DOMAIN_INTERNAL_STATE),POINTER :: DOMAIN_INT_STATE              !<-- The Domain component internal state pointer
      TYPE(WRAP_DOMAIN_INTERNAL_STATE) :: WRAP_DOMAIN                      !<-- The F90 wrap of the Domain internal state
      REAL(ESMF_KIND_R8),POINTER :: inst_sfc_pressure_fptr(:,:)         &  !<-- Instantaneous surface pressure pointer in cap
                                   ,inst_latent_htflx_fptr(:,:)         &  !<-- Instantaneous latent heat flux pointer in cap
                                   ,inst_sensible_htflx_fptr(:,:)       &  !<-- Instantaneous sensible heat flux pointer in cap
                                   ,inst_net_lwflx_fptr(:,:)            &  !<-- Instantaneous net lw flux pointer in cap
                                   ,inst_net_swflx_fptr(:,:)            &  !<-- Instantaneous net sw flux pointer in cap
                                   ,mean_zonal_momflx_fptr(:,:)         &  !<-- Mean zonal momentum flux pointer in cap
                                   ,mean_merid_momflx_fptr(:,:)         &  !<-- Mean meridional momentum flux pointer in cap
                                   ,mean_prec_rate_fptr(:,:)               !<-- Mean precipitation rate pointer in cap
!
      INTEGER :: fieldIdx                                                  !<-- Field ID from the Field List
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!     write(0,36611)
36611 format(' enter CONNECT export')
!
      IF(PRESENT(RC)) RC= ESMF_SUCCESS
!
!-----------------------------------------------------------------------
!***  Extract the NMM component's internal state.
!-----------------------------------------------------------------------
!
!
      CALL ESMF_GridCompGetInternalState(NMM_GRID_COMP                  &
                                        ,WRAP_NMM                       &
                                        ,RC )
!
      IF(RC/=0)THEN
        WRITE(0,*)' CONNECT_IMPORT_FIELDS failed to get NMM int state'
        WRITE(0,*)' RC=',RC
        CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
      ENDIF
!
      NMM_INT_STATE=>wrap_nmm%NMM_INT_STATE
!
!-----------------------------------------------------------------------
!

!     write(0,40021)domain_id,domain_descriptors(domain_id)%task_active
40021 format(' CONNECT domain_id=',i2,' task_active=',l1)
      IF(DOMAIN_DESCRIPTORS(DOMAIN_ID)%TASK_ACTIVE) THEN
!
        CALL ESMF_GridCompGetInternalState(&
               nmm_int_state%DOMAIN_GRID_COMP(DOMAIN_ID) &
              ,WRAP_DOMAIN                          &
              ,RC)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
!
        DOMAIN_INT_STATE=>wrap_domain%DOMAIN_INT_STATE
!
!------------------------------------------------------------
!***  Connect instantaneous sensible heat flux with the cap.
!------------------------------------------------------------
!
! Fei for hycom compatibility
!       fieldIdx = queryFieldList(exportFieldsList, 'inst_sensi_heat_flx', rc=rc)
        fieldIdx = queryFieldList(exportFieldsList, 'mean_sensi_heat_flx', rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
        CALL ESMF_FieldGet(DOMAIN_DESCRIPTORS(DOMAIN_ID)%exportFieldsList(fieldIdx) &
                          ,farrayPtr=inst_sensible_htflx_fptr                       &
                          ,rc       =rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
!
        domain_int_state%INST_SENS_HT_FLX_COUPLED => inst_sensible_htflx_fptr   !<-- Connect inst sensible heat flx between atm and cap
!
!----------------------------------------------------------
!***  Connect instantaneous latent heat flux with the cap.
!----------------------------------------------------------
!
! Fei for hycom compatibility
!        fieldIdx = queryFieldList(exportFieldsList, 'inst_laten_heat_flx', rc=rc)
        fieldIdx = queryFieldList(exportFieldsList, 'mean_laten_heat_flx', rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
        CALL ESMF_FieldGet(DOMAIN_DESCRIPTORS(DOMAIN_ID)%exportFieldsList(fieldIdx) &
                          ,farrayPtr=inst_latent_htflx_fptr                         &
                          ,rc       =rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
!
        domain_int_state%INST_LAT_HT_FLX_COUPLED => inst_latent_htflx_fptr  !<-- Connect mean latent heat flx between atm and cap
!
!-----------------------------------------------------
!***  Connect instantaneous net LW flux with the cap.
!-----------------------------------------------------
!
!keep        fieldIdx = queryFieldList(exportFieldsList, 'inst_net_lw_flx', rc=rc)
        fieldIdx = queryFieldList(exportFieldsList, 'mean_net_lw_flx', rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
        CALL ESMF_FieldGet(DOMAIN_DESCRIPTORS(DOMAIN_ID)%exportFieldsList(fieldIdx) &
                          ,farrayPtr=inst_net_lwflx_fptr                            &
                          ,rc       =rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
!
        domain_int_state%INST_NET_LW_FLX_COUPLED => inst_net_lwflx_fptr    !<-- Connect instantaneous net lw flx between atm and cap
!
!-----------------------------------------------------
!***  Connect instantaneous net SW flux with the cap.
!-----------------------------------------------------
!
!keep   fieldIdx = queryFieldList(exportFieldsList, 'inst_net_sw_flx', rc=rc)
        fieldIdx = queryFieldList(exportFieldsList, 'mean_net_sw_flx', rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
        CALL ESMF_FieldGet(DOMAIN_DESCRIPTORS(DOMAIN_ID)%exportFieldsList(fieldIdx) &
                          ,farrayPtr=inst_net_swflx_fptr                            &
                          ,rc       =rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
        domain_int_state%INST_NET_SW_FLX_COUPLED => inst_net_swflx_fptr   !<-- Connect instantaneous net sw flx between atm and cap
!
!----------------------------------------------------
!***  Connect mean zonal momentum flux with the cap.
!----------------------------------------------------
!
        fieldIdx = queryFieldList(exportFieldsList, 'mean_zonal_moment_flx', rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
        CALL ESMF_FieldGet(DOMAIN_DESCRIPTORS(DOMAIN_ID)%exportFieldsList(fieldIdx) &
                          ,farrayPtr=mean_zonal_momflx_fptr                            &
                          ,rc       =rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
!
        domain_int_state%MEAN_ZONAL_MOM_FLX_COUPLED => mean_zonal_momflx_fptr   !<-- Connect mean zonal mom flux between atm and cap
!
!---------------------------------------------------------
!***  Connect mean meridional momentum flux with the cap.
!---------------------------------------------------------
!
        fieldIdx = queryFieldList(exportFieldsList, 'mean_merid_moment_flx', rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
        CALL ESMF_FieldGet(DOMAIN_DESCRIPTORS(DOMAIN_ID)%exportFieldsList(fieldIdx) &
                          ,farrayPtr=mean_merid_momflx_fptr                         &
                          ,rc       =rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
!
        domain_int_state%MEAN_MERID_MOM_FLX_COUPLED => mean_merid_momflx_fptr   !<-- Connect mean meridionalal mom flux between atm and cap
!
!---------------------------------------------------
!***  Connect mean precipitation rate with the cap.
!---------------------------------------------------
!
        fieldIdx = queryFieldList(exportFieldsList, 'mean_prec_rate', rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
        CALL ESMF_FieldGet(DOMAIN_DESCRIPTORS(DOMAIN_ID)%exportFieldsList(fieldIdx) &
                          ,farrayPtr=mean_prec_rate_fptr                            &
                          ,rc       =rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
!
        domain_int_state%MEAN_PREC_RATE_COUPLED => mean_prec_rate_fptr     !<-- Connect mean precipitation rate between atm and cap
!
!----------------------------------------------------------
!***  Connect instantaneous surface pressure with the cap.
!----------------------------------------------------------
!
        fieldIdx = queryFieldList(exportFieldsList, 'inst_pres_height_surface', rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
        CALL ESMF_FieldGet(DOMAIN_DESCRIPTORS(DOMAIN_ID)%exportFieldsList(fieldIdx) &
                          ,farrayPtr=inst_sfc_pressure_fptr                         &
                          ,rc       =rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
!
        domain_int_state%INST_SFC_PRESSURE_COUPLED => inst_sfc_pressure_fptr   !<-- Connect instantaneous sfc pressure between atm and cap
!
!-----------------------------------------------------------------------
!
      ENDIF
!
!     write(0,36613)
36613 format(' exit CONNECT export')
!-----------------------------------------------------------------------
!
      END SUBROUTINE CONNECT_EXPORT_FIELDS
!
!-----------------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!-----------------------------------------------------------------------
!
      SUBROUTINE NMMB_RegridImport(slice, DOMAIN_DESCRIPTORS, NUM_DOMAINS_TOTAL, RC)
!
!-----------------------------------------------------------------------
!***  If there are nests interpolate the import field(s) from the
!***  upper parent's grid to those of the nests since for now only
!***  the parent is directly coupled to the external model(s).
!-----------------------------------------------------------------------

      INTEGER, INTENT(IN)                       :: slice                &
                                                  ,NUM_DOMAINS_TOTAL
      TYPE(DOMAIN_TASK_SPECS), POINTER          :: DOMAIN_DESCRIPTORS(:)
      INTEGER, INTENT(OUT), OPTIONAL            :: RC

      TYPE(ESMF_Field)                          :: parentField, nestField
      INTEGER                                   :: parent_domain_id, N, fieldIdx
      character(len=2)                          :: msg
      character(len=32)                          :: pname, nname

!-----------------------------------------------------------------------

      if(PRESENT(RC)) RC = ESMF_SUCCESS

      DO N = 2, NUM_DOMAINS_TOTAL
        parent_domain_id = DOMAIN_DESCRIPTORS(N)%PARENT_DOMAIN_ID
        !DO fieldIdx = 1, nImportFields
        DO fieldIdx = 3,3
          parentField=DOMAIN_DESCRIPTORS(parent_domain_id)%ImportFieldsList(fieldIdx)
          nestField  =DOMAIN_DESCRIPTORS(N)               %ImportFieldsList(fieldIdx)
          !print *, 'sizeIm = ', size(DOMAIN_DESCRIPTORS(1)%ImportFieldsList), &
          !         'sizeEx = ', size(DOMAIN_DESCRIPTORS(1)%ExportFieldsList), &
          !         fieldIdx, parent_domain_id, N
          call ESMF_FieldGet(parentField, name=pname, rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, &
            file=__FILE__)) &
            return  ! bail out
          call ESMF_FieldGet(nestField, name=nname, rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, &
            file=__FILE__)) &
            return  ! bail out
          !print *, 'pname = ', pname, ' nname = ', nname

          call ESMF_FieldRegrid(parentField, nestField, &
            routehandle=DOMAIN_DESCRIPTORS(N)%PARENT2SELF, &
            termorderflag=ESMF_TERMORDER_SRCSEQ, rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, &
            file=__FILE__)) &
            return  ! bail out
          write(msg, '(I2.2)') N
!          call ESMF_FieldWrite(nestField,'field_atm_import_sst_'//trim(msg)//'.nc',overwrite=.true.,&
!            timeslice=slice, rc=rc)
!          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
!            line=__LINE__, &
!            file=__FILE__)) &
!            return  ! bail out
        ENDDO
      ENDDO

!-----------------------------------------------------------------------

      END SUBROUTINE NMMB_RegridImport
!
!-----------------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!-----------------------------------------------------------------------
!
      SUBROUTINE NMMB_RegridExport(slice,DOMAIN_DESCRIPTORS             &
                                  ,NUM_DOMAINS_TOTAL,nExportFields_NMMB &
                                  ,EXPORT_FIELDS_INDX, RC)
!
!-----------------------------------------------------------------------
!***  If there are nests then blend the export fields from
!***  their grids onto that of the upper parent since
!***  for now only the parent is directly coupled to the
!***  external model(s).  The field 'mean_prec_rate' is
!***  only from the upper parent and contains no information
!***  from the nests.
!-----------------------------------------------------------------------

      INTEGER, INTENT(IN)                       :: slice                &
                                                  ,nExportFields_NMMB   &
                                                  ,NUM_DOMAINS_TOTAL
      INTEGER,DIMENSION(:),INTENT(IN)           :: EXPORT_FIELDS_INDX
      TYPE(DOMAIN_TASK_SPECS), POINTER          :: DOMAIN_DESCRIPTORS(:)
      INTEGER, INTENT(OUT), OPTIONAL            :: RC

      TYPE(ESMF_Field)                          :: nestField, parentField
      INTEGER                                   :: parent_domain_id, N, NN, fieldIdx
      character(len=2)                          :: msg

!-----------------------------------------------------------------------

      if(PRESENT(RC)) RC = ESMF_SUCCESS

      DO N = 2, NUM_DOMAINS_TOTAL
        parent_domain_id = DOMAIN_DESCRIPTORS(N)%PARENT_DOMAIN_ID
!       DO fieldIdx = 1, nExportFields
        DO NN = 1, nExportFields_NMMB
          fieldIdx = EXPORT_FIELDS_INDX(NN)
          nestField  =DOMAIN_DESCRIPTORS(N)               %ExportFieldsList(fieldIdx)
          parentField=DOMAIN_DESCRIPTORS(parent_domain_id)%ExportFieldsList(fieldIdx)
          call ESMF_FieldRegrid(nestField, parentField, &
            routehandle=DOMAIN_DESCRIPTORS(N)%SELF2PARENT, &
            zeroregion=ESMF_REGION_SELECT, &
            termorderflag=ESMF_TERMORDER_SRCSEQ, rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, &
            file=__FILE__)) &
            return  ! bail out
          write(msg, '(I2.2)') N
!         call ESMF_FieldWrite(nestField,'field_atm_import_sst_'//trim(msg)//'.nc',overwrite=.true.,&
!           timeslice=slice, rc=rc)
!         if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
!           line=__LINE__, &
!           file=__FILE__)) &
!           return  ! bail out
        ENDDO
      ENDDO

!-----------------------------------------------------------------------

      END SUBROUTINE NMMB_RegridExport
!
!-----------------------------------------------------------------------
!
      ! Dump the contents of one element of DOMAIN_DESCRIPTORS
      SUBROUTINE DUMP_DOMAIN_DESCRIPTOR(DOMAIN_DESCRIPTORS, ID_X)
        TYPE(DOMAIN_TASK_SPECS), POINTER             :: DOMAIN_DESCRIPTORS(:)
        INTEGER, INTENT(IN)                          :: ID_X

        INTEGER                                      :: N, RC
        CHARACTER(4096)                              :: TMPSTR
        CHARACTER(64)                                :: fname

        write(tmpstr, *) "DOMAIN ID_X = ", ID_X
        call ESMF_LogWrite(trim(tmpstr), ESMF_LOGMSG_INFO, rc=rc)
        write(tmpstr, *) "INPES = ", DOMAIN_DESCRIPTORS(ID_X)%INPES
        call ESMF_LogWrite(trim(tmpstr), ESMF_LOGMSG_INFO, rc=rc)
        write(tmpstr, *) "JNPES = ", DOMAIN_DESCRIPTORS(ID_X)%JNPES
        call ESMF_LogWrite(trim(tmpstr), ESMF_LOGMSG_INFO, rc=rc)
        write(tmpstr, *) "MIN_INDEX = ", DOMAIN_DESCRIPTORS(ID_X)%INDX_MIN
        call ESMF_LogWrite(trim(tmpstr), ESMF_LOGMSG_INFO, rc=rc)
        write(tmpstr, *) "MAX_INDEX = ", DOMAIN_DESCRIPTORS(ID_X)%INDX_MAX
        call ESMF_LogWrite(trim(tmpstr), ESMF_LOGMSG_INFO, rc=rc)
        write(tmpstr, *) "TASK_ACTIVE = ", DOMAIN_DESCRIPTORS(ID_X)%TASK_ACTIVE
        call ESMF_LogWrite(trim(tmpstr), ESMF_LOGMSG_INFO, rc=rc)
        write(tmpstr, *) "NUM_PETS = ", DOMAIN_DESCRIPTORS(ID_X)%NUM_PETS
        call ESMF_LogWrite(trim(tmpstr), ESMF_LOGMSG_INFO, rc=rc)
        write(tmpstr, *) "PET_MAP = ", size(DOMAIN_DESCRIPTORS(ID_X)%PET_MAP), DOMAIN_DESCRIPTORS(ID_X)%PET_MAP
        call ESMF_LogWrite(trim(tmpstr), ESMF_LOGMSG_INFO, rc=rc)

        write(tmpstr, *) "ImportFieldsList Size = ", SIZE(DOMAIN_DESCRIPTORS(ID_X)%ImportFieldsList)
        call ESMF_LogWrite(trim(tmpstr), ESMF_LOGMSG_INFO, rc=rc)
        DO N = 1, SIZE(DOMAIN_DESCRIPTORS(ID_X)%ImportFieldsList)
          Call ESMF_FieldGet(DOMAIN_DESCRIPTORS(ID_X)%ImportFieldsList(N), name=fname, RC=RC)
          write(tmpstr, *) "Import Field #", N, " NAME = ", TRIM(fname)
          call ESMF_LogWrite(trim(tmpstr), ESMF_LOGMSG_INFO, rc=rc)
        ENDDO

        write(tmpstr, *) "ExportFieldsList Size = ", SIZE(DOMAIN_DESCRIPTORS(ID_X)%ExportFieldsList)
        call ESMF_LogWrite(trim(tmpstr), ESMF_LOGMSG_INFO, rc=rc)
        DO N = 1, SIZE(DOMAIN_DESCRIPTORS(ID_X)%ExportFieldsList)
          Call ESMF_FieldGet(DOMAIN_DESCRIPTORS(ID_X)%ExportFieldsList(N), name=fname, RC=RC)
          write(tmpstr, *) "Export Field #", N, " NAME = ", TRIM(fname)
          call ESMF_LogWrite(trim(tmpstr), ESMF_LOGMSG_INFO, rc=rc)
        ENDDO
      END SUBROUTINE
!
!-----------------------------------------------------------------------
!
      ! Dump the contents of DOMAIN_DESCRIPTORS on a PET(optional)
      SUBROUTINE DUMP_DOMAIN_DESCRIPTORS(DOMAIN_DESCRIPTORS, PETNO)

        TYPE(DOMAIN_TASK_SPECS), POINTER             :: DOMAIN_DESCRIPTORS(:)
        INTEGER, INTENT(IN), OPTIONAL                :: PETNO

        INTEGER                                      :: ID_X

        IF(present(PETNO)) THEN
          IF(I_AM_PET(PETNO)) THEN
            DO ID_X = 1, SIZE(DOMAIN_DESCRIPTORS)
              CALL DUMP_DOMAIN_DESCRIPTOR(DOMAIN_DESCRIPTORS, ID_X)
            ENDDO
          ENDIF
        ELSE
          DO ID_X = 1, SIZE(DOMAIN_DESCRIPTORS)
            CALL DUMP_DOMAIN_DESCRIPTOR(DOMAIN_DESCRIPTORS, ID_X)
          ENDDO
        ENDIF
        call ESMF_LogFlush()

      END SUBROUTINE DUMP_DOMAIN_DESCRIPTORS
!
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
      FUNCTION I_AM_ROOT(RC)
!-----------------------------------------------------------------------

        integer, intent(out) :: RC
        logical :: I_AM_ROOT

        type(ESMF_VM) :: VM
        integer       :: lpet

        RC = ESMF_SUCCESS
        I_AM_ROOT = .FALSE.
        call ESMF_VMGetCurrent(VM, RC=RC)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
        call ESMF_VMGet(VM, localPet=lpet, RC=RC)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
        if(lpet == 0) I_AM_ROOT = .TRUE.

!-----------------------------------------------------------------------
      END FUNCTION I_AM_ROOT
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
      FUNCTION I_AM_PET(PETNO)
!-----------------------------------------------------------------------

        integer, intent(in),optional  :: PETNO
        logical :: I_AM_PET

        type(ESMF_VM) :: VM
        integer       :: lpet, PETNO_loc, RC

        RC = ESMF_SUCCESS
        I_AM_PET = .FALSE.
        call ESMF_VMGetCurrent(VM, RC=RC)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
        call ESMF_VMGet(VM, localPet=lpet, RC=RC)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
        if(present(PETNO)) THEN
          if(lpet == PETNO) I_AM_PET = .TRUE.
        else
          I_AM_PET = .TRUE.
        endif

!-----------------------------------------------------------------------
      END FUNCTION I_AM_PET
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
      SUBROUTINE NMMB_GridUpdate(nestingDomainGridcompIndex             &
                              ,DOMAIN_DESCRIPTORS                       &
                              ,AREA, MASK, GLON, GLAT, VLON, VLAT       &
                              ,RC )
!-----------------------------------------------------------------------
!     Update the pointer values saved in DOMAIN_DESCRIPTORS when
!     a nest moves. The pointer values are set up to point to memory
!     allocated in nest's Grid.
!-----------------------------------------------------------------------

      INTEGER, INTENT(IN)                :: nestingDomainGridcompIndex
      TYPE(DOMAIN_TASK_SPECS),POINTER    :: DOMAIN_DESCRIPTORS(:)
      REAL(kind=KDBL),DIMENSION(:,:),INTENT(IN) :: AREA, MASK
      REAL(kind=KDBL),DIMENSION(:,:),POINTER :: GLON, GLAT              &
                                               ,VLON, VLAT
      INTEGER, INTENT(OUT)               :: RC

      RC = ESMF_SUCCESS

      DOMAIN_DESCRIPTORS(nestingDomainGridcompIndex)%NEST_CELL_AREA = AREA
      DOMAIN_DESCRIPTORS(nestingDomainGridcompIndex)%NEST_SEA_MASK  = NINT(MASK)
      DOMAIN_DESCRIPTORS(nestingDomainGridcompIndex)%NEST_GLON = GLON
      DOMAIN_DESCRIPTORS(nestingDomainGridcompIndex)%NEST_GLAT = GLAT
      DOMAIN_DESCRIPTORS(nestingDomainGridcompIndex)%NEST_VLON = VLON
      DOMAIN_DESCRIPTORS(nestingDomainGridcompIndex)%NEST_VLAT = VLAT

      END SUBROUTINE NMMB_GridUpdate

!-----------------------------------------------------------------------
  subroutine Grid_Write(grid, string, rc)
    type(ESMF_Grid) ,intent(in)  :: grid
    character(len=*),intent(in)  :: string
    integer         ,intent(out) :: rc
  
    ! local 
    type(ESMF_Array)            :: array
    character(len=*),parameter  :: subname='(module_MEDIATOR:Grid_Write)'
    logical                     :: isPresent

    ! -- centers --

    rc = ESMF_SUCCESS

    ! -- centers --

    call ESMF_GridGetCoord(grid, staggerLoc=ESMF_STAGGERLOC_CENTER, isPresent=isPresent, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
    if (isPresent) then
      call ESMF_GridGetCoord(grid, coordDim=1, staggerLoc=ESMF_STAGGERLOC_CENTER, array=array, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
      call ESMF_ArraySet(array, name="lon_center", rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
      call ESMF_ArrayWrite(array, trim(string)//"_grid_coord1.nc", rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
      call ESMF_GridGetCoord(grid, coordDim=2, staggerLoc=ESMF_STAGGERLOC_CENTER, array=array, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
      call ESMF_ArraySet(array, name="lat_center", rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
      call ESMF_ArrayWrite(array, trim(string)//"_grid_coord2.nc", rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
    endif

    ! -- corners --

    call ESMF_GridGetCoord(grid, staggerLoc=ESMF_STAGGERLOC_CORNER, isPresent=isPresent, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
    if (isPresent) then
      call ESMF_GridGetCoord(grid, coordDim=1, staggerLoc=ESMF_STAGGERLOC_CORNER, array=array, rc=rc)
      if (.not. ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) then
        call ESMF_ArraySet(array, name="lon_corner", rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
        call ESMF_ArrayWrite(array, trim(string)//"_grid_corner1.nc", rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
      endif
      call ESMF_GridGetCoord(grid, coordDim=2, staggerLoc=ESMF_STAGGERLOC_CORNER, array=array, rc=rc)
      if (.not. ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) then
        call ESMF_ArraySet(array, name="lat_corner", rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
        call ESMF_ArrayWrite(array, trim(string)//"_grid_corner2.nc", rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
      endif
    endif


    ! -- mask --

    call ESMF_GridGetItem(grid, itemflag=ESMF_GRIDITEM_MASK, staggerLoc=ESMF_STAGGERLOC_CENTER, isPresent=isPresent, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
    if (isPresent) then
      call ESMF_GridGetItem(grid, staggerLoc=ESMF_STAGGERLOC_CENTER, itemflag=ESMF_GRIDITEM_MASK, array=array, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
      call ESMF_ArraySet(array, name="mask", rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
      call ESMF_ArrayWrite(array, trim(string)//"_grid_mask.nc", rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
    endif

    ! -- area --

    call ESMF_GridGetItem(grid, itemflag=ESMF_GRIDITEM_AREA, staggerLoc=ESMF_STAGGERLOC_CENTER, isPresent=isPresent, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
    if (isPresent) then
      call ESMF_GridGetItem(grid, staggerLoc=ESMF_STAGGERLOC_CENTER, itemflag=ESMF_GRIDITEM_AREA, array=array, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
      call ESMF_ArraySet(array, name="area", rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
      call ESMF_ArrayWrite(array, trim(string)//"_grid_area.nc", rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
    endif

  end subroutine Grid_Write

      END MODULE MODULE_DOMAIN_NUOPC_SET
!
!-----------------------------------------------------------------------
