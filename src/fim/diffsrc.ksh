#!/bin/ksh

#diffcmd=xxdiff
diffcmd=diff

#fimsrcnems=/export/Tasks/FIMtoNEMSrepo/FIM_r1710_WORK/FIMsrc/fim/framework/nems
#fimsrcnems=FIMsrc/fim/framework/nems
fimsrcnems=/export/Tasks/FIMtoNEMSrepo/FIM_r1911/FIMsrc/fim/framework/nems

echo "Compare $fimsrcnems with ."

$diffcmd ${fimsrcnems}/fim_grid_comp.F90 fim_grid_comp.F90
$diffcmd ${fimsrcnems}/fim_internal_state.F90 fim_internal_state.F90
$diffcmd ${fimsrcnems}/module_DYNAMICS_GRID_COMP.F90 module_DYNAMICS_GRID_COMP.F90
$diffcmd ${fimsrcnems}/module_PHYSICS_GRID_COMP.F90 module_PHYSICS_GRID_COMP.F90
$diffcmd ${fimsrcnems}/module_FIM_INTEGRATE.F90 module_FIM_INTEGRATE.F90
$diffcmd ${fimsrcnems}/module_DYN_PHY_CPL_COMP.F90 module_DYN_PHY_CPL_COMP.F90

