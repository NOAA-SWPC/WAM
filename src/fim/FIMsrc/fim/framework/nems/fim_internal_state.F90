      MODULE module_FIM_INTERNAL_STATE

      USE ESMF_MOD

      IMPLICIT NONE

      PRIVATE
      PUBLIC :: FIM_INTERNAL_STATE,WRAP_FIM_INTERNAL_STATE

      TYPE FIM_INTERNAL_STATE

        type(esmf_gridcomp) :: gc_fim_dyn
        type(esmf_gridcomp) :: gc_fim_phy
        type(esmf_state   ) :: imp_fim_dyn
        type(esmf_state   ) :: exp_fim_dyn      !<-- import/export states for fim dynamics
        type(esmf_state   ) :: imp_fim_phy
        type(esmf_state   ) :: exp_fim_phy      !<-- import/export states for fim physics
!       type(esmf_state   ) :: imp_fim_wrt
!       type(esmf_state   ) :: exp_fim_wrt      !<-- import/export states for fim write
        type(esmf_clock   ) :: clock_fim
        type(esmf_logical)  :: cpl_flag
!       type(esmf_logical)  :: chemistry_on                !<-- is chemistry active?
        type(esmf_cplcomp)  :: gc_fim_cpl
        integer             :: mype                        !<-- each mpi task id
!       integer             :: write_group_ready_to_go     !<-- the write group to use
        logical             :: quilting                    !<-- is asynchronous quilting specified?
        type(esmf_logical)  :: physics_on                  !<-- is physics active?
!       type(esmf_gridcomp), pointer :: wrt_comps(:)
        type(esmf_timeinterval) :: timeinterval_fim_output !<-- time interval between fim history output

        ! Task ID list of fcst tasks (for Dyn, Phy, and Cpl components)
        ! Task IDs are based on FIM component local VM
        integer, pointer :: petlist_fcst(:)
        ! Task ID list of all write tasks
        ! Task IDs are based on FIM component local VM
        integer, pointer :: petlist_write(:)

      END TYPE FIM_INTERNAL_STATE

      TYPE WRAP_FIM_INTERNAL_STATE
        TYPE(FIM_INTERNAL_STATE),POINTER :: FIM_INT_STATE
      END TYPE WRAP_FIM_INTERNAL_STATE

      END MODULE module_FIM_INTERNAL_STATE

