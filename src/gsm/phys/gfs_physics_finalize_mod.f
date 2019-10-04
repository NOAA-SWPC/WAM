!
! !module: gfs_physics_finalize_mod 
!          --- finalize module of the grided component of the gfs physics.
!
! !description: gfs finalize module.
!
! !revision history:
!
!  november 2004      weiyu yang initial code.
!  february 2006      shrinivas moorthi - removed some comments
!  january  2007      hann-ming henry juang -- modify to be dynamics only
!  july     2007      shrinivas moorthi -- modified for physics component
!  november 2007      hann-ming henry juang -- remove some calls
!
!
! !interface:
!
      module gfs_physics_finalize_mod
!
!!uses:
!
      use gfs_physics_internal_state_mod

      implicit none

      contains

      subroutine gfs_physics_finalize(gis_phy, rc)

      type(gfs_physics_internal_state)                :: gis_phy
      integer, optional,                intent(out)   :: rc

!
!***********************************************************************
!
      if(present(rc)) then
        rc = 0
      end if

      end subroutine gfs_physics_finalize

      end module gfs_physics_finalize_mod
