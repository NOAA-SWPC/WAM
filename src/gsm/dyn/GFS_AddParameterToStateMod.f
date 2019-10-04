!#include "../../../ESMFVersionDefine.h"

!----------------------------------------------------------------------
! !MODULE: GFS_AddParameterToStateMod
!        --- Add required parameters to the GFS ESMF export state
!            for the ensemble coupler to do the spectral transform
!            for the stochastic perturbation scheme, the second step.
!
! !DESCRIPTION: Add all required parameters to the GFS ESMF export state.
!
! !REVISION HISTORY:
!
!  May      2007     Weiyu Yang Initial code.
!  March    2009     Weiyu Yang Modified for the NEMS model.
!  February 2011     Weiyu Yang Updated to use both the ESMF 4.0.0rp2 library,
!                               ESMF 5 library and the the ESMF 3.1.0rp2 library.
!  May      2011     Weiyu yang, modified for using the ESMF 5.2.0r_beta_snapshot_07.
!
!
! !INTERFACE:
!

 MODULE GFS_AddParameterToStateMod

      USE ESMF

 USE gfs_dyn_resol_def
 USE gfs_dyn_layout1
 USE gfs_dyn_mpi_def
 USE gfs_dynamics_internal_state_mod

 IMPLICIT none

 CONTAINS

 SUBROUTINE AddParameterToState(State, Int_State, rc)

 TYPE(ESMF_State),                           INTENT(inout) :: State
 TYPE(gfs_dynamics_internal_state), POINTER, INTENT(in)    :: Int_State
 INTEGER, OPTIONAL,                          INTENT(out)   :: rc

 INTEGER                                         :: rc1, rcfinal

 rc1     = ESMF_SUCCESS
 rcfinal = ESMF_SUCCESS

! One by one add the parameters to the GFS ESMF export state.
!------------------------------------------------------------
 CALL ESMF_AttributeSet(State, 'NTRAC', ntrac, rc = rc1)
     IF(ESMF_LogFoundError(rc1, msg="Add ntrac to the GFS export state.")) THEN
         rcfinal = ESMF_FAILURE
         PRINT*, 'Error Happened When Adding ntrac to the GFS export state, rc = ', rc1
         rc1 = ESMF_SUCCESS
     END IF

 CALL ESMF_AttributeSet(State, 'MPI_R_MPI_R', MPI_R_MPI_R, rc = rc1)

     IF(ESMF_LogFoundError(rc1, msg="Add MPI_R_MPI_R to the GFS export state.")) THEN
         rcfinal = ESMF_FAILURE
         PRINT*, 'Error Happened When Adding MPI_R_MPI_R to the GFS export state, rc = ', rc1
         rc1 = ESMF_SUCCESS
     END IF

 CALL ESMF_AttributeSet(State, 'JCAP', jcap, rc = rc1)

     IF(ESMF_LogFoundError(rc1, msg="Add JCAP to the GFS export state.")) THEN
         rcfinal = ESMF_FAILURE
         PRINT*, 'Error Happened When Adding JCAP to the GFS export state, rc = ', rc1
         rc1 = ESMF_SUCCESS
     END IF

 CALL ESMF_AttributeSet(State, 'NODES_COMP', nodes, rc = rc1)

     IF(ESMF_LogFoundError(rc1, msg="Add NODES_COMP to the GFS export state.")) THEN
         rcfinal = ESMF_FAILURE
         PRINT*, 'Error Happened When Adding NODES_COMP to the GFS export state, rc = ', rc1
         rc1 = ESMF_SUCCESS
     END IF

 CALL ESMF_AttributeSet(State, 'ME_COMP', me, rc = rc1)

     IF(ESMF_LogFoundError(rc1, msg="Add ME_COMP to the GFS export state.")) THEN
         rcfinal = ESMF_FAILURE
         PRINT*, 'Error Happened When Adding ME_COMP to the GFS export state, rc = ', rc1
         rc1 = ESMF_SUCCESS
     END IF

 CALL ESMF_AttributeSet(State, 'MC_COMP', MC_COMP, rc = rc1)

     IF(ESMF_LogFoundError(rc1, msg="Add MC_COMP to the GFS export state.")) THEN
         rcfinal = ESMF_FAILURE
         PRINT*, 'Error Happened When Adding MC_COMP to the GFS export state, rc = ', rc1
         rc1 = ESMF_SUCCESS
     END IF

 CALL ESMF_AttributeSet(State, 'LATS_NODE_A', lats_node_a, rc = rc1)

     IF(ESMF_LogFoundError(rc1, msg="Add LATS_NODE_A to the GFS export state.")) THEN
         rcfinal = ESMF_FAILURE
         PRINT*, 'Error Happened When Adding LATS_NODE_A to the GFS export state, rc = ', rc1
         rc1 = ESMF_SUCCESS
     END IF

 CALL ESMF_AttributeSet(State, 'IPT_LATS_NODE_A', ipt_lats_node_a, rc = rc1)

     IF(ESMF_LogFoundError(rc1, msg="Add IPT_LATS_NODE_A to the GFS export state.")) THEN
         rcfinal = ESMF_FAILURE
         PRINT*, 'Error Happened When Adding IPT_LATS_NODE_A to the GFS export state, rc = ', rc1
         rc1 = ESMF_SUCCESS
     END IF

 CALL ESMF_AttributeSet(State, 'LONF', lonf, rc = rc1)

     IF(ESMF_LogFoundError(rc1, msg="Add LONF to the GFS export state.")) THEN
         rcfinal = ESMF_FAILURE
         PRINT*, 'Error Happened When Adding LONF to the GFS export state, rc = ', rc1
         rc1 = ESMF_SUCCESS
     END IF

 CALL ESMF_AttributeSet(State, 'LATG', latg, rc = rc1)

     IF(ESMF_LogFoundError(rc1, msg="Add LATG to the GFS export state.")) THEN
         rcfinal = ESMF_FAILURE
         PRINT*, 'Error Happened When Adding LATG to the GFS export state, rc = ', rc1
         rc1 = ESMF_SUCCESS
     END IF

 CALL ESMF_AttributeSet(State, 'GLOBAL_LATS_A', Int_State%global_lats_a, &
     itemCount = latg, rc = rc1)

     IF(ESMF_LogFoundError(rc1, msg="Add GLOBAL_LATS_A to the GFS export state.")) THEN
         rcfinal = ESMF_FAILURE
         PRINT*, 'Error Happened When Adding GLOBAL_LATS_A to the GFS export state, rc = ', rc1
         rc1 = ESMF_SUCCESS
     END IF

 CALL ESMF_AttributeSet(State, 'LONSPERLAT', Int_State%lonsperlat, &
     itemCount = latg, rc = rc1)

     IF(ESMF_LogFoundError(rc1, msg="Add LONSPERLAT to the GFS export state.")) THEN
         rcfinal = ESMF_FAILURE
         PRINT*, 'Error Happened When Adding LONSPERLAT to the GFS export state, rc = ', rc1
         rc1 = ESMF_SUCCESS
     END IF

 IF(rcfinal == ESMF_SUCCESS) THEN
     PRINT*, "PASS: GFS_AddParameterToStateMod.f"
 ELSE
     PRINT*, "FAIL: GFS_AddParameterToStateMod.f"
 END IF

 IF(PRESENT(rc)) THEN
     rc = rcfinal
 END IF

 END SUBROUTINE AddParameterToState

 END MODULE GFS_AddParameterToStateMod
