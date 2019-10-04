      module gfs_dyn_tracer_const
      use gfs_dyn_machine , only : kind_grid
      implicit none

! !revision history:
!
!  09Feb2010   Sarah Lu, ri/cpi dimension increased
!  17Aug2010   Sarah Lu, print tracer_const only for master PE
!  04Jan2012   Henry Juang, remove num_tracer

      integer, parameter :: max_num_tracer=50
!jw      real(kind=kind_grid) ri(0:20),cpi(0:20)
!     real(kind=kind_grid),target :: ri(0:20),cpi(0:20)
      real(kind=kind_grid),target ::  ri(0:max_num_tracer)
      real(kind=kind_grid),target :: cpi(0:max_num_tracer)
!hmhj integer, parameter :: num_tracer=3

      contains
! -------------------------------------------------------------------   
      subroutine get_tracer_const (ntrac,me,nlunit)
      use gfs_dyn_machine , only : kind_grid
      use gfs_dyn_physcons , only : rd => con_rd , cpd => con_cp
      implicit none
      integer ntrac,me,nlunit
      namelist /tracer_constant/ ri,cpi

!     print *,' enter get_tracer_const',ntrac,num_tracer
c
! Remark (09Feb2010) by Sarah Lu
! This routine reads ri/cpi for meteorological tracers;
! Number of met tracers (num_tracer) is hardwired to 3 

!hmhj if( ntrac.ne.num_tracer ) then
!hmhj   if( me.eq.0 ) then
!hmhj     write(0,*) ' Error ; inconsistent number of tracer '
!hmhj     write(0,*) ' ntrac=',ntrac,' num_tracer=',num_tracer
!hmhj   endif
!hmhj   call abort
!hmhj endif

      ri=0.0
      cpi=0.0
      ri(0)=rd
      cpi(0)=cpd

      rewind(nlunit)
      read(nlunit, tracer_constant)
      if(me==0) write(*, tracer_constant)

!     print *,' done get_tracer_const'

      return
      end subroutine get_tracer_const

      end module gfs_dyn_tracer_const
