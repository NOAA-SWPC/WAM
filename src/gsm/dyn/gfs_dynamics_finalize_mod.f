!
! !module: gfs_dynamics_finalize_mod 
!          --- finalize module of the grided
!              component of the gfs dynamics system.
!
! !description: gfs finalize module.
!
! !revision history:
!
!  november 2004      weiyu yang initial code.
!  february 2006      shrinivas moorthi - removed some comments
!  january  2007      hann-ming henry juang -- modify to be dynamics only
!  Sept     2010      Jun Wang   change gfsio to nemsio
!  Dec      2010      Jun Wang   change to nemsio library
!
!
! !interface:
!
      module gfs_dynamics_finalize_mod
!
!!uses:
!
      use gfs_dynamics_internal_state_mod
      use nemsio_module , only : nemsio_finalize

      implicit none

      contains

      subroutine gfs_dynamics_finalize(gis_dyn, rc)

      type(gfs_dynamics_internal_state)                :: gis_dyn
      integer, optional,                 intent(out)   :: rc

!
!***********************************************************************
!
!     if(me.eq.0) then
!       call w3tage('gsm     ')
!     endif
!!
      if (nemsio_out .or. nemsio_in) then
        call nemsio_finalize()
      endif

      if(present(rc)) then
        rc = 0
      end if

      end subroutine gfs_dynamics_finalize

      end module gfs_dynamics_finalize_mod
