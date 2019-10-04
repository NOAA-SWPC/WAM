!
!! ! Subroutine : dyn_gocart_tracer_config
!
! ! Description: initialize GOCART chem tracers
!
! ! Revision history:
!   Aug 09 2011   Jun Wang, initial code from tracer_config_init
!   Sep 16 2011   Sarah Lu, pass Chem_Registry info to gfs_dyn_tracer
! -------------------------------------------------------------------------
!
      subroutine dyn_gocart_tracer_config (me)
!      subroutine dyn_gocart_tracer_config (gfs_dyn_tracer,me)
!
      use gfs_dyn_tracer_config
      use Chem_RegistryMod
!
      implicit none
! input
      integer, intent(in)    ::  me
! output
!      type (gfs_dyn_tracer_type), intent(out)    ::  gfs_dyn_tracer
!
! local
      integer                 :: i, status, ierr
      type(Chem_Registry)     :: reg

      if(me==0)print *,'LU_TRC: tracer_config is running'

! Read Chem_Registry
      reg = Chem_RegistryCreate ( ierr )

! ntrac_chem = number of chem tracers
      gfs_dyn_tracer%ntrac_chem = reg%nq
      gfs_dyn_tracer%doing_OC = reg%doing_OC
      gfs_dyn_tracer%doing_BC = reg%doing_BC
      gfs_dyn_tracer%doing_DU = reg%doing_DU
      gfs_dyn_tracer%doing_SS = reg%doing_SS
      gfs_dyn_tracer%doing_SU = reg%doing_SU
      gfs_dyn_tracer%doing_GOCART = reg%doing_GOCART

!--- fill in chem tracers
      allocate(gfs_dyn_tracer%chem_name(reg%nq),stat=status)
      if( status .ne. 0 ) go to 999

      do i = 1, reg%nq
       gfs_dyn_tracer%chem_name(i)=reg%vname(i)
      enddo

! Destroy Chem_Registry
      call Chem_RegistryDestroy ( reg, ierr )

      return

999   print *,'LU_TRC: error in allocate gfs_dyn_tracer :',status,me

      end subroutine dyn_gocart_tracer_config

