!
! !module: gfs_dynamics_namelist_mod  ---       definition of the name list
!                                               in the esmf internal state.
!
! !description:   define the name list variables
!                 in the esmf internal state.
!---------------------------------------------------------------------------
! !revision history:
!
!  november 2004      weiyu yang initial code.
!  february 2006      took out model namelists
!  january  2007      hann-ming henry juang for gfs dynamics only
!  oct 2009           sarah lu, tracer added, (q, oz, cld) removed
!
! !interface:
!
      module gfs_dynamics_namelist_mod


      implicit none

      type nam_gfs_dyn_namelist
           integer                :: nlunit, total_member, member_id
           real                   :: deltim
           character(80)          :: gfs_dyn_namelist
           character(20)          :: grid_ini,grid_ini2,sig_ini, sig_ini2
      end type nam_gfs_dyn_namelist
!
      type gfs_dyn_state_namelist
!          LOGICAL                :: redgg_a
!
! for the sigma file.
!--------------------
           integer                :: idate1_import
           integer                :: z_import
           integer                :: ps_import
           integer                :: u_import
           integer                :: v_import
           integer                :: temp_import
           integer                :: tracer_import
           integer                :: p_import
           integer                :: dp_import
           integer                :: dpdt_import
     
           integer                :: idate1_export
           integer                :: z_export
           integer                :: ps_export
           integer                :: u_export
           integer                :: v_export
           integer                :: temp_export
           integer                :: tracer_export
           integer                :: p_export
           integer                :: dp_export
           integer                :: dpdt_export
           integer                :: shum_wts_export
           integer                :: sppt_wts_export
           integer                :: skeb_wts_export
           integer                :: vc_wts_export

      end type gfs_dyn_state_namelist

      end module gfs_dynamics_namelist_mod
