!
!! ! Subroutine : gocart_tracer_config
!
! ! Description: initialize gocart chem tracers
!
! ! Revision history:
!   Aug 09 2011   Jun Wang, initial code from tracer_config_init
!   Sep 16 2011   Sarah Lu, pass Chem_Registry info to gfs_phy_tracer
! -------------------------------------------------------------------------
!
      subroutine gocart_tracer_config (me)
!      subroutine gocart_tracer_config (gfs_phy_tracer,me)
!
      use  gfs_phy_tracer_config
      use Chem_RegistryMod
!
      implicit none
!
! input
      integer, intent(in)    ::  me
! output
!      type (gfs_phy_tracer_type), intent(out)    ::  gfs_phy_tracer
!
! local
      integer                 :: i, status, ierr
      type(Chem_Registry)     :: reg

      if(me==0)print *,'LU_TRC: tracer_config is running'

! Read Chem_Registry
      reg = Chem_RegistryCreate ( ierr )

!!    if ( me == 0) call Chem_RegistryPrint (reg)

! ntrac_chem = number of chem tracers
      gfs_phy_tracer%ntrac_chem = reg%nq
      gfs_phy_tracer%doing_OC = reg%doing_OC
      gfs_phy_tracer%doing_BC = reg%doing_BC
      gfs_phy_tracer%doing_DU = reg%doing_DU
      gfs_phy_tracer%doing_SS = reg%doing_SS
      gfs_phy_tracer%doing_SU = reg%doing_SU
      gfs_phy_tracer%doing_GOCART = reg%doing_GOCART

!--- fill in chem tracers
      allocate(gfs_phy_tracer%chem_name(reg%nq),stat=status)
      if( status .ne. 0 ) go to 999

      do i = 1, reg%nq
       gfs_phy_tracer%chem_name(i)=reg%vname(i)
      enddo

! Destroy Chem_Registry
      call Chem_RegistryDestroy ( reg, ierr )

      return

999   print *,'TRAC_CONFIG: error in allocate gfs_phy_tracer :',
     &    status,me

      end subroutine gocart_tracer_config
