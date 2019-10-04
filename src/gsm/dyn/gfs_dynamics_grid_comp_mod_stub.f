!#include "../../../ESMFVersionDefine.h"

! !module: gfs_dynamics_grid_comp_mod --- 
!                       esmf gridded component of gfs dynamics
!
! !description: gfs dynamics gridded component main module.
!
! !revision history:
!
!  january 2007     hann-ming henry juang
!                           
!
! !interface:
!
      module gfs_dynamics_grid_comp_mod
 
!!uses:
!------
      USE ESMF

      implicit none

      private   ! by default, data is private to this module

      public gfs_dyn_setservices	! only set service is public

!eop
!-------------------------------------------------------------------------


      contains


!----------------------------------------------------------------------
!bop
!
! !routine: gfs_dyn_setservices --- 
!           set services for gfs dynamics gridded component.
! 
! !interface:
!
      subroutine gfs_dyn_setservices (gc_gfs_dyn, rc)
 
! !arguments:
!------------

      type(esmf_gridcomp), intent(in)  :: gc_gfs_dyn 	! gridded component
      integer,             intent(out) :: rc    	! return code
     

      end subroutine gfs_dyn_setservices

!----------------------------------------------------------------------
!bop
! !routine:  gfs_dyn_initialize --- initialize routine to initialize 
!                                   and set up the gfs running job.
!
! !description: this subroutine initializes the gfs running before
!               the main running loop.
!
!
! !revision history:
!
!  november 2004     weiyu yang initial code.
!  may      2005     weiyu yang for the updated gfs version.
!  february 2006     moorthi
!  february 2007     h.-m. h. juang
!
! !interface:
!

! this argument list is a standard list for all the initialize,
! the run and finalize routines for an esmf system.
!--------------------------------------------------------------
      subroutine gfs_dyn_initialize(gc_gfs_dyn, 			&
                                   imp_gfs_dyn, exp_gfs_dyn, clock, rc)

! user code, for computations related to the esmf interface states.
!------------------------------------------------------------------
!
! !input/output variables and parameters:
!----------------------------------------

      type(esmf_gridcomp), intent(inout) :: gc_gfs_dyn 
      type(esmf_state),    intent(inout) :: imp_gfs_dyn
      type(esmf_state),    intent(inout) :: exp_gfs_dyn
      type(esmf_clock),    intent(inout) :: clock

!
! !output variables and parameters:
!----------------------------------

      integer, intent(out) :: rc  


      end subroutine gfs_dyn_initialize


!----------------------------------------------------------------------
!bop
!
! !routine: gfs_dyn_run --- 
!           main grid component routine to run the gfs dynamics.
!
! !description: this subroutine will run the most part computations 
!               of the gfs dynamics.
!
! !revision history:
!
!  november 2004     weiyu yang initial code.
!  may      2005     weiyu yang for the updated gfs version.
!  february 2006     moorthi
!  july     2007     hann-ming henry juang
!
! !interface:
!

      subroutine gfs_dyn_run(gc_gfs_dyn, 				&
                            imp_gfs_dyn, exp_gfs_dyn, clock, rc)

!
! !input variables and parameters:
!---------------------------------
      type(esmf_gridcomp), intent(inout) :: gc_gfs_dyn   
      type(esmf_state),    intent(in)    :: imp_gfs_dyn 
 
! !output variables and parameters:
!----------------------------------
      type(esmf_clock),    intent(inout) :: clock
      type(esmf_state),    intent(inout) :: exp_gfs_dyn
      integer,             intent(out)   :: rc   


      end subroutine gfs_dyn_run


!----------------------------------------------------------------------
!bop
!
! !routine: finalize --- finalizing routine to finish the 
!                        gfs running job.
!
! !description: this subroutine will finish the gfs computations,
! !             and will release the memory space.
!
! !revision history:
!
!  november 2004     weiyu yang initial code.
!  may      2005     weiyu yang for the updated gfs version.
!  february 2006     moorthi
!  february 2007     juang for dynamics only
!
! !interface:

      subroutine gfs_dyn_finalize(gc_gfs_dyn, 				&
                                 imp_gfs_dyn, exp_gfs_dyn, clock, rc)

!
! !input variables and parameters:
!---------------------------------
      type(esmf_gridcomp), intent(inout)  :: gc_gfs_dyn
      type(esmf_state),    intent(inout)  :: imp_gfs_dyn
      type(esmf_state),    intent(inout)  :: exp_gfs_dyn
      type(esmf_clock),    intent(inout)  :: clock

! !output variables and parameters:
!----------------------------------
      integer,             intent(out)    :: rc

      end subroutine gfs_dyn_finalize
      end module gfs_dynamics_grid_comp_mod
