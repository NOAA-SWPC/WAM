!
! !MODULE: gfs_physics_g2d_mod  ---      Definition of the 2d aerosol
!                                        diag fields in the ESMF internal state.
!
! !DESCRIPTION: gfs_physics_g2d_mod ---    Define the 2d aerosol diag fields
!                                             in the ESMF internal state.
!---------------------------------------------------------------------------
! !REVISION HISTORY:
!
!  2010/07/14     Sarah Lu,  Initial code.
!  2010/08/10     Sarah Lu,  Add g2d_zerout routine
!  2010/09/11     Sarah Lu,  Modify g2d_zerout 
!  2010/10/10     Sarah Lu,  Add GFS forcing fields to g2d_fld
!  2010/10/15     Sarah Lu,  Revise g2d_zerout
!  2010/12/23     Sarah Lu,  Fields not used by GOCART are removed from g2d_fld;
!                            set doing_MET to F
!  2014/11/04     Jun Wang,  Fields duexttauv,not used by GOCART are removed from g2d_fld;
!
! !INTERFACE:
!
 MODULE gfs_physics_g2d_mod
!
 use machine, only: kind_phys
 use gfs_phy_tracer_config,      ONLY: gfs_phy_tracer_type
!
 IMPLICIT none

 INTEGER, PARAMETER        :: MAX_AER_DIAG=150

 TYPE AER_R2D
  real(kind=kind_phys), dimension(:,:), pointer :: flds =>null()
  character*10             :: name
 END TYPE AER_R2D

 TYPE AER_Diag_Data
  TYPE (AER_R2D), DIMENSION(MAX_AER_DIAG)       :: diag 
  integer                                       :: nfld 
 END TYPE AER_Diag_Data

 TYPE G2D_Var_Data
   TYPE (AER_Diag_Data)         :: DU
   TYPE (AER_Diag_Data)         :: SS
   TYPE (AER_Diag_Data)         :: SU
   TYPE (AER_Diag_Data)         :: OC 
   TYPE (AER_Diag_Data)         :: BC
   TYPE (AER_Diag_Data)         :: MET
 END TYPE G2D_Var_Data
!
 logical, public      :: doing_MET = .false.
!
 contains

! !IROUTINE: g2d_aldata ---

!---------------------------------------------------------------------------
    subroutine g2d_aldata (dim1, dim2, gfs_phy_tracer,  &
                           g2d_fld, iret)

    implicit none
! 
    integer, intent(in)                   :: dim1, dim2
    type(gfs_phy_tracer_type), intent(in) :: gfs_phy_tracer
    TYPE(G2D_Var_Data), INTENT(out)       :: g2d_fld
    integer, intent(out)                  :: iret
!
    integer n

!! .............................................................
!!
!! For GEOS-5, MAPL will allocate these diag fields on-the-fly
!! Since NCEP does not port MAPL, the allocation will be done
!! manually in phys grid component where these fields are outputted
!!
 character*10      ::   name_du(25), name_su(30)
 character*10      ::   name_oc(15), name_bc(14)
 character*10      ::   name_ss(29)
 character*10      ::   name_met(24)
!
 data name_du(1:25) /                                    &
  'DUEM001', 'DUEM002', 'DUEM003', 'DUEM004', 'DUEM005', &
  'DUSD001', 'DUSD002', 'DUSD003', 'DUSD004', 'DUSD005', &
  'DUDP001', 'DUDP002', 'DUDP003', 'DUDP004', 'DUDP005', &
  'DUWT001', 'DUWT002', 'DUWT003', 'DUWT004', 'DUWT005', &
  'DUSMASS', 'DUCMASS', 'DUSMASS25','DUCMASS25',         &
  'DUAERIDX'/

 data name_su(1:30) /                                    &
  'SUEM001', 'SUEM002', 'SUEM003', 'SUEM004',            &
  'SUDP001', 'SUDP002', 'SUDP003', 'SUDP004',            &
  'SUWT001', 'SUWT002', 'SUWT003', 'SUWT004',            &
  'SO2SMASS','SO2CMASS','SO4SMASS','SO4CMASS',           &
  'DMSSMASS','DMSCMASS','SUPSO2',                        &
  'SUPSO4g', 'SUPSO4aq','SUPSO4wt',                      &
  'SO4EMAN', 'SO2EMAN', 'SO2EMBB', 'SO2EMVN','SO2EMVE',  &
  'SUPMSA',  'SUEXTTAU','SUSCATAU' /

 data name_oc(1:15) /                                    &
  'OCEM001', 'OCEM002', 'OCDP001', 'OCDP002', 'OCWT001', &
  'OCWT002', 'OCHYPHIL','OCEMAN',  'OCEMBB',  'OCEMBF',  &
  'OCEMBG',  'OCSMASS', 'OCCMASS', 'OCEXTTAU','OCSCATAU'/

 data name_bc(1:14) /                                    &
  'BCEM001', 'BCEM002', 'BCDP001', 'BCDP002', 'BCWT001', &
  'BCWT002', 'BCHYPHIL','BCEMAN',  'BCEMBB', 'BCEMBF',   &
  'BCSMASS', 'BCCMASS', 'BCEXTTAU','BCSCATAU'/

 data name_ss(1:29) /                                    &
   'SSEM001', 'SSEM002', 'SSEM003', 'SSEM004', 'SSEM005',&
   'SSSD001', 'SSSD002', 'SSSD003', 'SSSD004', 'SSSD005',&
   'SSDP001', 'SSDP002', 'SSDP003', 'SSDP004', 'SSDP005',&
   'SSWT001', 'SSWT002', 'SSWT003', 'SSWT004', 'SSWT005',&
   'SSSMASS', 'SSCMASS', 'SSEXTTAU', 'SSSCATAU',         &
   'SSSMASS25','SSCMASS25','SSEXTT25','SSSCAT25',        &
   'SSAERIDX'/

 data name_met(1:24) /                                   &
   'U10M', 'V10M', 'UUSTAR', 'Z0H', 'LWI', 'ZPBL',       &
   'WET1', 'SH'  , 'CNPRCP', 'NCNPRCP' ,                 &
   'ZLE01' , 'AIRDEN01', 'T01'    , 'U01',               &
   'V01'   , 'FCLD01'   , 'DQDT01' ,                     &
   'ZLE64' , 'AIRDEN64', 'T64'    , 'U64',               &
   'V64'   , 'FCLD64'   , 'DQDT64'/
!
!! .............................................................
!
!   allocate g2d_fld for GFS met fields
    g2d_fld%met%nfld = 0
    if ( doing_MET ) then
      call aldata_ (dim1, dim2, name_met, g2d_fld%met, iret)
    endif

!   allocate 2d aer diag fields for DU module
    g2d_fld%du%nfld = 0
    if ( gfs_phy_tracer%doing_DU ) then
      call aldata_ (dim1, dim2, name_du, g2d_fld%du, iret)
    endif

!   allocate 2d aer diag fields for SU module
    g2d_fld%su%nfld = 0
    if ( gfs_phy_tracer%doing_SU ) then
      call aldata_ (dim1, dim2, name_su, g2d_fld%su, iret)
    endif
!
!   allocate 2d aer diag fields for SS module
    g2d_fld%ss%nfld = 0
    if ( gfs_phy_tracer%doing_SS ) then
      call aldata_ (dim1, dim2, name_ss, g2d_fld%ss, iret)
    endif

!   allocate 2d aer diag fields for OC module
    g2d_fld%oc%nfld = 0
    if ( gfs_phy_tracer%doing_OC ) then
      call aldata_ (dim1, dim2, name_oc, g2d_fld%oc, iret)
    endif

!   allocate 2d aer diag fields for BC module
    g2d_fld%bc%nfld = 0
    if ( gfs_phy_tracer%doing_BC ) then
      call aldata_ (dim1, dim2, name_bc, g2d_fld%bc, iret)
    endif

    return
    end subroutine

!!!
    subroutine aldata_ (dim1, dim2, vname, g2d, iret)

    implicit none
! 
    integer, intent(in)                :: dim1, dim2
    character*10, intent(in)           :: vname(:)
    TYPE (AER_Diag_Data),intent (out)  :: g2d
    integer, intent(out)               :: iret
!
! local
    integer                            :: n, num_vname
!
    num_vname = size (vname, dim=1)
!
!   allocate 2d aerosol diag fields: diag(:)%flds
!   fill in diag field name:  diag(:)%name
!   fill in diag field count: nfld

    do n = 1, num_vname
      allocate(g2d%diag(n)%flds(dim1,dim2), stat=iret)
      if(iret.ne.0) iret=-3
      g2d%diag(n)%name = vname(n)
    enddo
!
    g2d%nfld = num_vname
!
    return
    end subroutine


!---------------------------------------------------------------------------
! !IROUTINE: g2d_zerout ---

    subroutine g2d_zerout (g2d_fld, iret)

    implicit none
!
    TYPE(G2D_Var_Data), INTENT(inout)     :: g2d_fld
    integer, intent(out)                  :: iret
!
    integer k
!
    if (  g2d_fld%met%nfld .gt. 0 ) then
      do k =1, g2d_fld%met%nfld
        g2d_fld%met%diag(k)%flds(:,:) = 0.
      enddo
    endif
!
    if (  g2d_fld%du%nfld .gt. 0 ) then
      do k =1, g2d_fld%du%nfld
        g2d_fld%du%diag(k)%flds(:,:) = 0.
      enddo
    endif

    if (  g2d_fld%su%nfld .gt. 0 ) then
      do k =1, g2d_fld%su%nfld
        g2d_fld%su%diag(k)%flds(:,:) = 0.
      enddo
    endif

    if (  g2d_fld%ss%nfld .gt. 0 ) then
      do k =1, g2d_fld%ss%nfld
       g2d_fld%ss%diag(k)%flds(:,:) = 0.
      enddo
    endif

    if (  g2d_fld%oc%nfld .gt. 0 ) then
      do k =1, g2d_fld%oc%nfld
       g2d_fld%oc%diag(k)%flds(:,:) = 0.
      enddo
    endif

    if (  g2d_fld%bc%nfld .gt. 0 ) then
      do k =1, g2d_fld%bc%nfld
       g2d_fld%bc%diag(k)%flds(:,:) = 0.
      enddo
    endif

    iret = 0

    return
    end subroutine
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 END MODULE gfs_physics_g2d_mod
