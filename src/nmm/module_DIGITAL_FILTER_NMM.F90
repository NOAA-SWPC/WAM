#define ESMF_ERR_ABORT(rc) if (ESMF_LogFoundError(rc, msg="Aborting NMMB", line=__LINE__, file=__FILE__)) call ESMF_Finalize(endflag=ESMF_END_ABORT)
!-----------------------------------------------------------------------
!
      MODULE module_DIGITAL_FILTER_NMM
!
!-----------------------------------------------------------------------
!
! a generic digital filter for any model under ESMF
!
!-----------------------------------------------------------------------
! March    2007	Hann-Ming Henry Juang
! February 2008 Weiyu Yang, updated to use the ESMF 3.1.0 library.
! February 2011 Weiyu Yang, Updated to use both the ESMF 4.0.0rp2 library,
!                           ESMF 5 library and the the ESMF 3.1.0rp2 library.
! May      2011 Weiyu Yang, Modified for using the ESMF 5.2.0r_beta_snapshot_07.
! September2011 Weiyu Yang, Modified for using the ESMF 5.2.0r library.
! July     2012    T Black, Modified for generational task usage.
!----------------------------------------------------------------------------
!
      USE ESMF
      use module_kinds
      use module_exchange,only: halo_exch

      implicit none

      private

      public :: digital_filter_dyn_init_nmm,                            &
                digital_filter_dyn_sum_nmm,                             &
                digital_filter_dyn_average_nmm,                         &
                digital_filter_phy_init_nmm,                            &
                digital_filter_phy_save_nmm,                            &
                digital_filter_phy_restore_nmm

! ---------
! dynamics
! ---------
      character(20), allocatable, save :: name_save_2d(:)
      character(20), allocatable, save :: name_save_3d(:)

! ---------
! physics
! ---------
      character(20), allocatable, save :: name_save_2d_phys(:)
      character(20), allocatable, save :: name_save_3d_phys(:)

      contains

!-----------------------------------------------------------------------
!#######################################################################
!-----------------------------------------------------------------------
      subroutine digital_filter_dyn_init_nmm(filt_bundle                &
                                            ,ndfistep                   &
                                            ,dt_int,dt_num,dt_den       &
                                            ,its,ite,jts,jte,lm         &
                                            ,tot_rank_2d                &
                                            ,tot_rank_3d                &
                                            ,kstep,nstep                &
                                            ,totalsum                   &
                                            ,dolph_wgts                 &
                                            ,array_save_2d              &
                                            ,array_save_3d )
!-----------------------------------------------------------------------
!
      type(ESMF_FieldBundle),intent(inout) :: filt_bundle
      integer(kind=kint), intent(in)       :: ndfistep
      integer(kind=kint), intent(in)       :: dt_int,dt_num,dt_den
      integer(kind=kint), intent(in)       :: its,ite,jts,jte,lm
!
      integer(kind=kint), intent(out) :: kstep,nstep
      integer(kind=kint), intent(out) :: tot_rank_2d                    &
                                        ,tot_rank_3d
!
      real(kind=kfpt), dimension(:), pointer, intent(inout) :: dolph_wgts
      real(kind=kfpt), dimension(:,:,:), pointer, intent(inout) :: array_save_2d
      real(kind=kfpt), dimension(:,:,:,:), pointer, intent(inout) :: array_save_3d
!
      real(kind=kfpt), intent(out) :: totalsum
!
      integer(kind=kint)           :: tmp_rank
      integer(kind=kint)           :: istat,rc,n,NUM_FIELDS
      character(20)                :: field_name
      real(kind=kfpt)              :: taus, dt
      type(ESMF_Field)             :: tmpfield
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
      nstep = ndfistep
      kstep = - nstep -1

      if (associated(dolph_wgts)) deallocate(dolph_wgts)
      allocate(dolph_wgts(-nstep:nstep),stat=istat)
      if(istat/=0)then
        write(0,*)' DIGITAL_FILTER_DYN_INIT_NMM failed to allocate dolph_wgts stat=',istat
        write(0,*)' Aborting!!'
        CALL ESMF_Finalize(endflag=ESMF_END_ABORT,rc=RC); ESMF_ERR_ABORT(RC)
      endif

      dt=float(dt_int)+float(dt_num)/float(dt_den)

!     hardwiring cutoff frequency based on length of filtering window

      taus=float(2*ndfistep)*dt
!
      CALL dolph(dt,taus,nstep,dolph_wgts)
!

!     Retrieve the dynamical fields to be filtered from the bundle

      CALL ESMF_FieldBundleGet(FIELDBUNDLE      = FILT_BUNDLE           &  !<-- The ESMF Bundle of arrays to be filtered
                              ,fieldCount       = NUM_FIELDS            &  !<-- # of Fields in the Bundle
                              ,rc               = RC); ESMF_ERR_ABORT(RC)

      tot_rank_2d=0
      tot_rank_3d=0

      if (.not. allocated(name_save_2d))                                &
      allocate(name_save_2d(NUM_FIELDS))

      if (.not. allocated(name_save_3d))                                &
      allocate(name_save_3d(NUM_FIELDS))

      DO N=1,NUM_FIELDS

        CALL ESMF_FieldBundleGet(FIELDBUNDLE    = FILT_BUNDLE           &  !<-- The ESMF Bundle of arrays to be filtered
                                ,fieldIndex     = N                     &
                                ,field          = tmpfield              &
                                ,rc             = RC); ESMF_ERR_ABORT(RC)

        CALL ESMF_FieldGet(field=tmpfield, name=field_name, dimCount=tmp_rank, rc=rc); ESMF_ERR_ABORT(RC)
!
        IF (tmp_rank == 2) THEN
          tot_rank_2d=tot_rank_2d+1
          name_save_2d(tot_rank_2d)=field_name
        ENDIF
!
        IF (tmp_rank == 3) THEN
          tot_rank_3d=tot_rank_3d+1
          name_save_3d(tot_rank_3d)=field_name
        ENDIF
!
      ENDDO
!
      IF (tot_rank_2d > 0 .and. .not. associated(array_save_2d)) THEN
        allocate(array_save_2d(ITS:ITE,JTS:JTE,tot_rank_2d),stat=istat)
        if(istat/=0)then
          write(0,*)' DIGITAL_FILTER_DYN_INIT_NMM failed to allocate array_save_2d stat=',istat
          write(0,*)' Aborting!!'
          CALL ESMF_Finalize(endflag=ESMF_END_ABORT,rc=RC); ESMF_ERR_ABORT(RC)
        endif
        array_save_2d=0.
      ENDIF
!
      IF (tot_rank_3d > 0 .and. .not. associated(array_save_3d)) THEN
        allocate(array_save_3d(ITS:ITE,JTS:JTE,LM+1,tot_rank_3d),stat=istat)
        if(istat/=0)then
          write(0,*)' DIGITAL_FILTER_DYN_INIT_NMM failed to allocate array_save_3d stat=',istat
          write(0,*)' Aborting!!'
          CALL ESMF_Finalize(endflag=ESMF_END_ABORT,rc=RC); ESMF_ERR_ABORT(RC)
        endif
        array_save_3d=0.
      ENDIF
!
      totalsum=0.
!
!-----------------------------------------------------------------------
!
      end subroutine digital_filter_dyn_init_nmm

!-----------------------------------------------------------------------
!#######################################################################
!-----------------------------------------------------------------------
      subroutine digital_filter_dyn_sum_nmm(filt_bundle                 &
                                           ,mean_on                     &
                                           ,its,ite,jts,jte,lm          &
                                           ,tot_rank_2d                 &
                                           ,tot_rank_3d                 &
                                           ,kstep,nstep                 &
                                           ,totalsum                    &
                                           ,dolph_wgts                  &
                                           ,array_save_2d               &
                                           ,array_save_3d )
!-----------------------------------------------------------------------

      type(ESMF_FieldBundle),intent(inout) :: filt_bundle
      integer(kind=kint),intent(in) :: mean_on,nstep
      integer(kind=kint),intent(in) :: its,ite,jts,jte,lm
      integer(kind=kint), intent(in) :: tot_rank_2d,tot_rank_3d
!
      integer(kind=kint),intent(inout) :: kstep
!
      real(kind=kfpt), dimension(:), pointer, intent(in) :: dolph_wgts
!
      real(kind=kfpt), intent(inout) :: totalsum
!
      real(kind=kfpt), dimension(:,:,:), pointer, intent(inout) :: array_save_2d
      real(kind=kfpt), dimension(:,:,:,:), pointer, intent(inout) :: array_save_3d
!
      integer(kind=kint) :: i,j,l,n,rc,LDIM
      real(kind=kfpt) :: digfil,prod,sx,wx
      real(kind=kfpt),dimension(:,:)    ,pointer :: hold_2d
      real(kind=kfpt),dimension(:,:,:)  ,pointer :: hold_3d
      logical :: dolph
!
      type(ESMF_Field) :: hold_field
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
      rc    =esmf_success
!
      kstep = kstep + 1

!-----------------------------------------------------------------------
!     Future task: make this dolph switch logical a configure file item
!-----------------------------------------------------------------------

      dolph=.true.

      IF (dolph) THEN
        digfil=dolph_wgts(kstep)

      ELSE
        sx = acos(-1.)*kstep/nstep
        wx = acos(-1.)*kstep/(nstep+1)
        if( kstep/=0)then
          digfil = sin(wx)/wx*sin(sx)/sx
        else
          digfil=1.
        endif

        if(mean_on>0)then
          digfil=1.
        endif

      ENDIF
!
      totalsum = totalsum + digfil
!
      if(tot_rank_2d>0) then
        do n=1,tot_rank_2d
          NULLIFY(HOLD_2D)

          CALL ESMF_FieldBundleGet(FIELDBUNDLE = FILT_BUNDLE            &  !<-- The ESMF Bundle of arrays to be filtered
                                  ,fieldName   = name_save_2d(N)        &
                                  ,field       = HOLD_FIELD             &
                                  ,rc          = RC); ESMF_ERR_ABORT(RC)

          CALL ESMF_FieldGet(field     = HOLD_FIELD                     &  !<-- Field that holds the data pointer
                            ,localDe   = 0                              &
                            ,farrayPtr = HOLD_2D                        &  !<-- Put the pointer here
                            ,rc        = RC); ESMF_ERR_ABORT(RC)

          do j=jts,jte
          do i=its,ite
            array_save_2d(i,j,n)=array_save_2d(i,j,n)+digfil*hold_2d(i,j)
          enddo
          enddo
        enddo
      endif
!
      if(tot_rank_3d>0)then
        do n=1,tot_rank_3d
          NULLIFY(HOLD_3D)

          CALL ESMF_FieldBundleGet(FIELDBUNDLE = FILT_BUNDLE            &  !<-- The ESMF Bundle of arrays to be filtered
                                  ,fieldName   = name_save_3d(N)        &
                                  ,field       = HOLD_FIELD             &
                                  ,rc          = RC); ESMF_ERR_ABORT(RC)

          CALL ESMF_FieldGet(field     = HOLD_FIELD                     &  !<-- Field that holds the data pointer
                            ,localDe   = 0                              &
                            ,farrayPtr = HOLD_3D                        &  !<-- Put the pointer here
                            ,rc        = RC); ESMF_ERR_ABORT(RC)

          LDIM=size(HOLD_3D,dim=3)

          do l=1,LDIM
            do j=jts,jte
            do i=its,ite
              array_save_3d(i,j,l,n)=array_save_3d(i,j,l,n)+digfil*hold_3d(i,j,l)
            enddo
            enddo
          enddo
        enddo
      endif
!
!-----------------------------------------------------------------------
!
      end subroutine digital_filter_dyn_sum_nmm

!-----------------------------------------------------------------------
!#######################################################################
!-----------------------------------------------------------------------
      subroutine digital_filter_dyn_average_nmm(filt_bundle             &
                                               ,its,ite,jts,jte,lm      &
                                               ,tot_rank_2d             &
                                               ,tot_rank_3d             &
                                               ,kstep,nstep             &
                                               ,totalsum                &
                                               ,array_save_2d           &
                                               ,array_save_3d )
!
!-----------------------------------------------------------------------
!
      TYPE(ESMF_FieldBundle),INTENT(INOUT) :: FILT_BUNDLE
      INTEGER(kind=KINT),INTENT(IN) :: ITS,ITE,JTS,JTE,LM
!
      INTEGER(kind=KINT),INTENT(INOUT) :: TOT_RANK_2D                   &
                                         ,TOT_RANK_3D
!
      INTEGER(kind=KINT),INTENT(OUT) :: KSTEP,NSTEP
!
      REAL(kind=KFPT),INTENT(IN) :: TOTALSUM
!
      REAL(kind=KFPT),DIMENSION(:,:,:),POINTER,INTENT(INOUT) :: ARRAY_SAVE_2D
      REAL(kind=KFPT),DIMENSION(:,:,:,:),POINTER,INTENT(INOUT) :: ARRAY_SAVE_3D
!
      INTEGER(KIND=KINT) :: I,J,L,N,RC,LDIM
      REAL(KIND=KFPT),DIMENSION(:,:)    ,POINTER :: HOLD_2D
      REAL(KIND=KFPT),DIMENSION(:,:,:)  ,POINTER :: HOLD_3D
!
      CHARACTER(20)    :: FIELD_NAME
!
      TYPE(ESMF_Field) :: HOLD_FIELD
!
      REAL(KIND=KFPT)  :: totalsumi
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
      RC    =ESMF_SUCCESS
!
      totalsumi = 1.0 / totalsum

      IF (tot_rank_2d > 0) THEN
        DO N=1,tot_rank_2d
          DO J=JTS,JTE
          DO I=ITS,ITE
            array_save_2d(I,J,N)=totalsumi*array_save_2d(I,J,N)
          ENDDO
          ENDDO
        ENDDO
      ENDIF
!
      IF (tot_rank_3d > 0) THEN
        DO N=1,tot_rank_3d
          FIELD_NAME=name_save_3d(N)

          IF (FIELD_NAME == 'PINT') THEN
            LDIM=LM+1
          ELSE
            LDIM=LM
          ENDIF

          DO L=1,LDIM
            DO J=JTS,JTE
            DO I=ITS,ITE
              array_save_3d(I,J,L,N)=totalsumi*array_save_3d(I,J,L,N)
            ENDDO
            ENDDO
          ENDDO
        ENDDO
      ENDIF
!
      IF (tot_rank_2d > 0) THEN
        DO N=1,tot_rank_2d
          NULLIFY(HOLD_2D)

          CALL ESMF_FieldBundleGet(FIELDBUNDLE = FILT_BUNDLE            &  !<-- The ESMF Bundle of arrays to be filtered
                                  ,fieldName   = name_save_2d(N)        &
                                  ,field       = HOLD_FIELD             &
                                  ,rc          = RC); ESMF_ERR_ABORT(RC)

          CALL ESMF_FieldGet(field     = HOLD_FIELD                     &  !<-- Field that holds the data pointer
                            ,localDe   = 0                              &
                            ,farrayPtr = HOLD_2D                        &  !<-- Put the pointer here
                            ,rc        = RC); ESMF_ERR_ABORT(RC)

          DO J=JTS,JTE
          DO I=ITS,ITE
            HOLD_2D(I,J)=array_save_2d(I,J,N)
          ENDDO
          ENDDO

          CALL HALO_EXCH(hold_2d,1,2,2)
        ENDDO
      ENDIF
!
      IF (tot_rank_3d > 0) THEN
        DO N=1,tot_rank_3d
          NULLIFY(HOLD_3D)

          CALL ESMF_FieldBundleGet(FIELDBUNDLE = FILT_BUNDLE            &  !<-- The ESMF Bundle of arrays to be filtered
                                  ,fieldName   = name_save_3d(N)        &
                                  ,field       = HOLD_FIELD             &
                                  ,rc          = RC); ESMF_ERR_ABORT(RC)

          CALL ESMF_FieldGet(field     = HOLD_FIELD                     &  !<-- Field that holds the data pointer
                            ,localDe   = 0                              &
                            ,farrayPtr = HOLD_3D                        &  !<-- Put the pointer here
                            ,rc        = RC); ESMF_ERR_ABORT(RC)

          LDIM=size(HOLD_3D,dim=3)

          DO L=1,LDIM
            DO J=JTS,JTE
            DO I=ITS,ITE
              HOLD_3D(I,J,L)=array_save_3d(I,J,L,N)
            ENDDO
            ENDDO
          ENDDO

          CALL HALO_EXCH(hold_3d,LDIM,2,2)
        ENDDO
      ENDIF
!
      IF (tot_rank_2d > 0) THEN
        deallocate(array_save_2d)
      ENDIF

      IF (tot_rank_3d > 0) THEN
        deallocate(array_save_3d)
      ENDIF

      tot_rank_2d=0
      tot_rank_3d=0
      kstep=0
      nstep=0

!-----------------------------------------------------------------------
      end subroutine digital_filter_dyn_average_nmm
!-----------------------------------------------------------------------
!#######################################################################
!-----------------------------------------------------------------------
      subroutine digital_filter_phy_init_nmm(filt_bundle                &
                                            ,its,ite,jts,jte,lm         &
                                            ,tot_rank_2d_phys           &
                                            ,tot_rank_3d_phys           &
                                            ,array_save_2d_phys         &
                                            ,array_save_3d_phys )
!-----------------------------------------------------------------------
      TYPE(ESMF_FieldBundle),INTENT(INOUT) :: FILT_BUNDLE
      INTEGER(kind=KINT),INTENT(IN) :: ITS,ITE,JTS,JTE,LM
!
      REAL(kind=KFPT),DIMENSION(:,:,:),POINTER,INTENT(INOUT) :: ARRAY_SAVE_2D_PHYS
      REAL(kind=KFPT),DIMENSION(:,:,:,:),POINTER,INTENT(INOUT) :: ARRAY_SAVE_3D_PHYS
!
      INTEGER(kind=KINT),INTENT(OUT) :: TOT_RANK_2D_PHYS                &
                                       ,TOT_RANK_3D_PHYS
!
      TYPE(ESMF_Field) :: tmpfield
      integer :: istat, rc, NUM_FIELDS, n, tmp_rank
      character(len=20) :: field_name
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------

!     Retrieve the dynamical fields to be filtered from the bundle

      CALL ESMF_FieldBundleGet(FIELDBUNDLE    = FILT_BUNDLE             &  !<-- The ESMF Bundle of arrays to be filtered
                              ,fieldCount     = NUM_FIELDS              &  !<-- # of Fields in the Bundle
                              ,rc             = RC); ESMF_ERR_ABORT(RC)
!
      tot_rank_2d_phys=0
      tot_rank_3d_phys=0

      if (.not. allocated(name_save_2d_phys))                           &
      allocate(name_save_2d_phys(NUM_FIELDS))

      if (.not. allocated(name_save_3d_phys))                           &
      allocate(name_save_3d_phys(NUM_FIELDS))

      DO N=1,NUM_FIELDS

        CALL ESMF_FieldBundleGet(FIELDBUNDLE = FILT_BUNDLE              &  !<-- The ESMF Bundle of arrays to be filtered
                                ,fieldindex  = N                        &
                                ,field       = tmpfield                 &
                                ,rc          = RC); ESMF_ERR_ABORT(RC)

        CALL ESMF_FieldGet(field=tmpfield, name=field_name, dimCount=tmp_rank, rc=rc); ESMF_ERR_ABORT(RC)
!
        IF (tmp_rank == 2) THEN
          tot_rank_2d_phys=tot_rank_2d_phys+1
          name_save_2d_phys(tot_rank_2d_phys)=field_name
        ENDIF
!
        IF (tmp_rank == 3) THEN
          tot_rank_3d_phys=tot_rank_3d_phys+1
          name_save_3d_phys(tot_rank_3d_phys)=field_name
        ENDIF
!
      ENDDO
!
      IF (tot_rank_2d_phys > 0 .and. .not. associated(array_save_2d_phys)) THEN
        allocate(array_save_2d_phys(ITS:ITE,JTS:JTE,tot_rank_2d_phys),stat=istat)
        if(istat/=0)then
          write(0,*)' DIGITAL_FILTER_PHY_INIT_NMM failed to allocate array_save_2d_phys stat=',istat
          write(0,*)' Aborting!!'
          CALL ESMF_Finalize(endflag=ESMF_END_ABORT,rc=RC); ESMF_ERR_ABORT(RC)
        endif
        array_save_2d_phys=0.
      ENDIF
!
      IF (tot_rank_3d_phys > 0 .and. .not. associated(array_save_3d_phys)) THEN
        allocate(array_save_3d_phys(ITS:ITE,JTS:JTE,LM,tot_rank_3d_phys),stat=istat)
        if(istat/=0)then
          write(0,*)' DIGITAL_FILTER_PHY_INIT_NMM failed to allocate array_save_3d_phys stat=',istat
          write(0,*)' Aborting!!'
          CALL ESMF_Finalize(endflag=ESMF_END_ABORT,rc=RC); ESMF_ERR_ABORT(RC)
        endif
        array_save_3d_phys=0.
      ENDIF
!
!-----------------------------------------------------------------------
!
      end subroutine digital_filter_phy_init_nmm

!-----------------------------------------------------------------------
!#######################################################################
!-----------------------------------------------------------------------
      subroutine digital_filter_phy_save_nmm(filt_bundle                &
                                            ,its,ite,jts,jte            &
                                            ,tot_rank_2d_phys           &
                                            ,tot_rank_3d_phys           &
                                            ,array_save_2d_phys         &
                                            ,array_save_3d_phys )
!
!-----------------------------------------------------------------------
!
      TYPE(ESMF_FieldBundle),INTENT(INOUT) :: FILT_BUNDLE
      INTEGER(kind=KINT),INTENT(IN) :: ITS,ITE,JTS,JTE
      INTEGER(kind=KINT),INTENT(IN) :: TOT_RANK_2D_PHYS                 &
                                      ,TOT_RANK_3D_PHYS
!
      REAL(kind=KFPT),DIMENSION(:,:,:),POINTER,INTENT(INOUT) :: ARRAY_SAVE_2D_PHYS
      REAL(kind=KFPT),DIMENSION(:,:,:,:),POINTER,INTENT(INOUT) :: ARRAY_SAVE_3D_PHYS
!
      TYPE(ESMF_Field)             :: HOLD_FIELD
      integer                      :: n, rc, i,j,l, LDIM
      REAL(KIND=KFPT),DIMENSION(:,:)    ,POINTER :: HOLD_2D
      REAL(KIND=KFPT),DIMENSION(:,:,:)  ,POINTER :: HOLD_3D
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
      if(tot_rank_2d_phys>0) then
        do n=1,tot_rank_2d_phys
          nullify(hold_2d)

          CALL ESMF_FieldBundleGet(FIELDBUNDLE = FILT_BUNDLE            &  !<-- The ESMF Bundle of arrays to be filtered
                                  ,fieldName   = name_save_2d_phys(N)   &
                                  ,field       = HOLD_FIELD             &
                                  ,rc          = RC); ESMF_ERR_ABORT(RC)

          CALL ESMF_FieldGet(field     = HOLD_FIELD                     &  !<-- Field that holds the data pointer
                            ,localDe   = 0                              &
                            ,farrayPtr = HOLD_2D                        &  !<-- Put the pointer here
                            ,rc        = RC); ESMF_ERR_ABORT(RC)

          do j=jts,jte
          do i=its,ite
            array_save_2d_phys(i,j,n)=hold_2d(i,j)
          enddo
          enddo

        enddo
      endif

      if(tot_rank_3d_phys>0)then
        do n=1,tot_rank_3d_phys
          nullify(hold_3d)

          CALL ESMF_FieldBundleGet(FIELDBUNDLE = FILT_BUNDLE            &  !<-- The ESMF Bundle of arrays to be filtered
                                  ,fieldName   = name_save_3d_phys(N)   &
                                  ,field       = HOLD_FIELD             &
                                  ,rc          = RC); ESMF_ERR_ABORT(RC)

          CALL ESMF_FieldGet(field     = HOLD_FIELD                     &  !<-- Field that holds the data pointer
                            ,localDe   = 0                              &
                            ,farrayPtr = HOLD_3D                        &  !<-- Put the pointer here
                            ,rc        = RC); ESMF_ERR_ABORT(RC)

          LDIM=size(HOLD_3D,dim=3)

          do l=1, LDIM
            do j=jts,jte
            do i=its,ite
              array_save_3d_phys(i,j,l,n)=hold_3d(i,j,l)
            enddo
            enddo
          enddo
        enddo
!
      endif
!
!-----------------------------------------------------------------------
!
      end subroutine digital_filter_phy_save_nmm

!-----------------------------------------------------------------------
!#######################################################################
!-----------------------------------------------------------------------
      subroutine digital_filter_phy_restore_nmm( FILT_BUNDLE            &
                                                ,its,ite,jts,jte        &
                                                ,tot_rank_2d_phys       &
                                                ,tot_rank_3d_phys       &
                                                ,array_save_2d_phys     &
                                                ,array_save_3d_phys )
!
!-----------------------------------------------------------------------
!
      TYPE(ESMF_FieldBundle),INTENT(INOUT) :: FILT_BUNDLE
      INTEGER(kind=KINT),INTENT(IN) :: ITS,ITE,JTS,JTE
      INTEGER(kind=KINT),INTENT(IN) :: TOT_RANK_2D_PHYS                 &
                                      ,TOT_RANK_3D_PHYS
!
      REAL(kind=KFPT),DIMENSION(:,:,:),POINTER,INTENT(INOUT) :: ARRAY_SAVE_2D_PHYS
      REAL(kind=KFPT),DIMENSION(:,:,:,:),POINTER,INTENT(INOUT) :: ARRAY_SAVE_3D_PHYS
!
      TYPE(ESMF_Field) :: HOLD_FIELD
      INTEGER(KIND=KINT) :: I,J,L,N,RC,LDIM
      REAL(KIND=KFPT),DIMENSION(:,:)    ,POINTER :: HOLD_2D
      REAL(KIND=KFPT),DIMENSION(:,:,:)  ,POINTER :: HOLD_3D
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
      IF (tot_rank_2d_phys > 0) THEN
        DO N=1,tot_rank_2d_phys
          NULLIFY(HOLD_2D)

          CALL ESMF_FieldBundleGet(FIELDBUNDLE = FILT_BUNDLE            &  !<-- The ESMF Bundle of arrays to be filtered
                                  ,fieldName   = name_save_2d_phys(N)   &
                                  ,field       = HOLD_FIELD             &
                                  ,rc          = rc); ESMF_ERR_ABORT(RC)

          CALL ESMF_FieldGet(field     = HOLD_FIELD                     &  !<-- Field that holds the data pointer
                            ,localDe   = 0                              &
                            ,farrayPtr = HOLD_2D                        &  !<-- Put the pointer here
                            ,rc        = RC); ESMF_ERR_ABORT(RC)

          DO J=JTS,JTE
          DO I=ITS,ITE
            HOLD_2D(I,J)=array_save_2d_phys(I,J,N)
          ENDDO
          ENDDO

          CALL HALO_EXCH(HOLD_2D,1,2,2)

        ENDDO
      ENDIF

      IF (tot_rank_3d_phys > 0) THEN
        DO N=1,tot_rank_3d_phys
          NULLIFY(HOLD_3D)

          CALL ESMF_FieldBundleGet(FIELDBUNDLE = FILT_BUNDLE            &  !<-- The ESMF Bundle of arrays to be filtered
                                  ,fieldName   = name_save_3d_phys(N)   &
                                  ,field       = HOLD_FIELD             &
                                  ,rc          = rc); ESMF_ERR_ABORT(RC)

          CALL ESMF_FieldGet(field     = HOLD_FIELD                     &  !<-- Field that holds the data pointer
                            ,localDe   = 0                              &
                            ,farrayPtr = HOLD_3D                        &  !<-- Put the pointer here
                            ,rc        = RC); ESMF_ERR_ABORT(RC)
!
          LDIM=size(HOLD_3D,dim=3)
!
          DO L=1,LDIM
            DO J=JTS,JTE
            DO I=ITS,ITE
              HOLD_3D(I,J,L)=array_save_3d_phys(I,J,L,N)
            ENDDO
            ENDDO
          ENDDO

          CALL HALO_EXCH(HOLD_3D,LDIM,2,2)

        ENDDO
      ENDIF
!
!-----------------------------------------------------------------------
!
      end subroutine digital_filter_phy_restore_nmm

!-----------------------------------------------------------------------

   SUBROUTINE dolph(deltat, taus, m, window)

!     calculation of dolph-chebyshev window or, for short,
!     dolph window, using the expression in the reference:
!
!     antoniou, andreas, 1993: digital filters: analysis,
!     design and applications. mcgraw-hill, inc., 689pp.
!
!     the dolph window is optimal in the following sense:
!     for a given main-lobe width, the stop-band attenuation
!     is minimal; for a given stop-band level, the main-lobe
!     width is minimal.

      IMPLICIT NONE

      ! Arguments
      INTEGER, INTENT(IN)                  ::  m
      REAL, DIMENSION(0:2*M), INTENT(OUT)  ::  window
      REAL, INTENT(IN)                     :: deltat, taus

      ! local data
      integer, PARAMETER        :: NMAX = 5000
      REAL, dimension(0:NMAX)   :: t, w, time
      real, dimension(0:2*nmax) :: w2
      INTEGER                   :: NPRPE=0        ! no of pe

      real    :: pi, thetas, x0, term1, term2, rr, r,db, sum, arg, sumw
      integer :: n, nm1, i, nt

      PI = 4*ATAN(1.D0)

!      print *, 'in dfcoef, deltat = ', deltat, 'taus=',taus

      N = 2*M+1
      NM1 = N-1

      THETAS = 2*PI*ABS(DELTAT/TAUS)
      X0 = 1/COS(THETAS/2)
      TERM1 = (X0 + SQRT(X0**2-1))**(FLOAT(N-1))
      TERM2 = (X0 - SQRT(X0**2-1))**(FLOAT(N-1))
      RR = 0.5*(TERM1+TERM2)
      R = 1/RR
      DB = 20*LOG10(R)


!      WRITE(0,'(1X,''DOLPH: M,N='',2I8)')M,N
!      WRITE(0,'(1X,''DOLPH: THETAS (STOP-BAND EDGE)='',F10.3)')THETAS
!      WRITE(0,'(1X,''DOLPH: R,DB='',2F10.3)')R, DB

      DO NT=0,M
         SUM = 1
         DO I=1,M
            ARG = X0*COS(I*PI/N)
            CALL CHEBY(T,NM1,ARG)
            TERM1 = T(NM1)
            TERM2 = COS(2*NT*PI*I/N)
            SUM = SUM + R*2*TERM1*TERM2
         ENDDO
         W(NT) = SUM/N
         TIME(NT) = NT
      ENDDO
!     fill in the negative-time values by symmetry.
      DO NT=0,M
         W2(M+NT) = W(NT)
         W2(M-NT) = W(NT)
      ENDDO

!     fill up the array for return
      SUMW = 0.
      DO NT=0,2*M
         SUMW = SUMW + W2(NT)
      ENDDO
!      WRITE(0,'(1X,''DOLPH: SUM OF WEIGHTS W2='',F10.4)')SUMW

      DO NT=0,2*M
         WINDOW(NT) = W2(NT)
      ENDDO

      RETURN

   END SUBROUTINE dolph


   SUBROUTINE cheby(t, n, x)

!     calculate all chebyshev polynomials up to order n
!     for the argument value x.

!     reference: numerical recipes, page 184, recurrence
!         t_n(x) = 2xt_{n-1}(x) - t_{n-2}(x) ,  n>=2.

      IMPLICIT NONE

      ! Arguments
      INTEGER, INTENT(IN)  :: n
      REAL, INTENT(IN)     :: x
      REAL, DIMENSION(0:N) :: t

      integer  :: nn

      T(0) = 1
      T(1) = X
      IF(N.LT.2) RETURN
      DO NN=2,N
         T(NN) = 2*X*T(NN-1) - T(NN-2)
      ENDDO

      RETURN

   END SUBROUTINE cheby

      end module module_digital_filter_nmm
!
!-----------------------------------------------------------------------

