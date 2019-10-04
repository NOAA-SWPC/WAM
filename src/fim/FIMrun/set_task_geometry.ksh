#!/bin/ksh
#
# Set "LSB_PJL_TASK_GEOMETRY" for NCAR LSF installation.  
# Set LoadLeveler "#@task_geometry" for standard IBM installations.  
#


# turn this on for debugging
#export VERBOSE="true"
export VERBOSE="false"


function fail
{
  print "$@"
  exit 1
}

# Distribute MPI tasks across nodes as evenly as possible using 
# LSF "LSB_PJL_TASK_GEOMETRY" and LoadLeveler "task_geometry" sub-syntax:  
#  (taskid,taskid,...)(taskid,taskid,...)...
# In this syntax each "()" encloses a node and each "taskid" is an MPI task 
# ID in MPI_COMM_WORLD.  
#
# Arguments are:  
#   num_tasks           Number of MPI tasks.  
#   max_tasks_per_node  Maximum number of tasks per node.  
#
# Global (instance) variables:  
#   $dist_tasks         Return string is stored therein.  
#   $current_task       Next unused MPI task ID is incremented.  
#
function DistributeTasks
{
  # local variables
  typeset num_tasks=$1; shift
  typeset max_tasks_per_node=$1; shift

  if [[ $VERBOSE == "true" ]] ; then
    print "Begin DistributeTasks:  current_task = $current_task"
    print "Begin DistributeTasks:  num_tasks = $num_tasks"
    print "Begin DistributeTasks:  max_tasks_per_node = $max_tasks_per_node"
  fi

  if (( max_tasks_per_node <= 0 )) ; then
    fail "DistributeTasks ERROR:  max_tasks_per_node must be positive"
  fi

  if (( num_tasks <= 0 )) ; then
    dist_tasks=""
  else
    typeset num_nodes=0
    typeset extra_tasks=0
    typeset base_tasks_per_node=0
    (( num_nodes = num_tasks / max_tasks_per_node ))
    (( extra_tasks = num_tasks % max_tasks_per_node ))
    if (( extra_tasks > 0 )) ; then
      (( num_nodes += 1 ))
    fi
    (( base_tasks_per_node = num_tasks / num_nodes ))
    (( extra_tasks = num_tasks % num_nodes ))
    if [[ $VERBOSE == "true" ]] ; then
      print "  DistributeTasks:  num_nodes = $num_nodes"
      print "  DistributeTasks:  extra_tasks = $extra_tasks"
      print "  DistributeTasks:  base_tasks_per_node = $base_tasks_per_node"
    fi
    dist_tasks="("
    typeset task_count=0
    typeset node_task_count=0
    typeset this_tasks_per_node=0
    if (( base_tasks_per_node > 0 )) ; then
# refactor to remove duplication
      if (( extra_tasks > 0 )) ; then
        (( this_tasks_per_node = base_tasks_per_node + 1 ))
        (( extra_tasks -= 1 ))
      else
        (( this_tasks_per_node = base_tasks_per_node ))
      fi
    else
      (( this_tasks_per_node = extra_tasks ))
    fi
    while (( task_count < num_tasks )) ; do
      if (( node_task_count == this_tasks_per_node )) ; then
        # start a new node
        dist_tasks="$dist_tasks)($current_task"
        node_task_count=0
        if (( extra_tasks > 0 )) ; then
          (( this_tasks_per_node = base_tasks_per_node + 1 ))
          (( extra_tasks -= 1 ))
        else
          (( this_tasks_per_node = base_tasks_per_node ))
        fi
      else
        if (( node_task_count > 0 )) ; then
          dist_tasks="$dist_tasks,"
        fi
        dist_tasks="$dist_tasks$current_task"
      fi
      (( node_task_count += 1 ))
      (( current_task += 1 ))
      (( task_count += 1 ))
    done
    dist_tasks="$dist_tasks)"
  fi
}


# Map MPI tasks to LSF IBM nodes using "LSB_PJL_TASK_GEOMETRY" syntax:  
#  export LSB_PJL_TASK_GEOMETRY="{(taskid,taskid,...)(taskid,taskid,...)... }"
# In this syntax each "()" encloses a node and each "taskid" is an MPI task ID 
# in MPI_COMM_WORLD.  
# This syntax also matches LoadLeveler's #@task_geometry keyword.  
#
# Nodes are filled as evenly as possible.  
#
# Arguments are:  
#   compute_tasks     Number of compute tasks.  
#   nwt               Number of write tasks.  
#   mwtpn             Maximum number of write tasks per node.  
#   mctpn             Maximum number of compute tasks per node.  
#
# Global (instance) variables used:  
#   $dist_tasks         Return value from function DistributeTasks. 
#   $current_task       Next unused MPI task ID, updated by function 
#                       DistributeTasks. 
#
# Return value is printed.  
#
function SetTaskGeometry
{
  # local variables
  typeset compute_tasks="$1"; shift
  typeset nwt=$1; shift
  typeset mwtpn=$1; shift
  typeset mctpn=$1; shift

  if [[ $VERBOSE == "true" ]] ; then
    print "Begin SetTaskGeometry:  compute_tasks = $compute_tasks"
    print "Begin SetTaskGeometry:  nwt = $nwt"
    print "Begin SetTaskGeometry:  mwtpn = $mwtpn"
    print "Begin SetTaskGeometry:  mctpn = $mctpn"
  fi

  # First MPI task is compute root which lives on its own separate node
  current_task=0
  DistributeTasks 1 $mctpn
  typeset task_geometry="$dist_tasks"
  # Next set of MPI tasks are write tasks:  
  DistributeTasks $nwt $mwtpn
  task_geometry=$task_geometry$dist_tasks
  # Remaining tasks are compute tasks, distribute among nodes as evenly as 
  # possible.  
  (( remaining_ct = compute_tasks - 1 ))
  DistributeTasks $remaining_ct $mctpn
  task_geometry=$task_geometry$dist_tasks
  print "{$task_geometry}"
}


#set -x

usagestr="Usage:  ${0} compute_tasks [nwt [mwtpn [mctpn]]]"

nwt=0
mctpn=32
(( mwtpn = mctpn ))
case $# in
  1) compute_tasks="$1";;
  2) compute_tasks="$1"; nwt=$2;;
  3) compute_tasks="$1"; nwt=$2; mwtpn=$3;;
  4) compute_tasks="$1"; nwt=$2; mwtpn=$3; mctpn=$4;;
  *) fail "${usagestr}";;
esac 

if (( compute_tasks <= 0 )) ; then
  fail "${usagestr}, compute_tasks must be a positive integer"
fi

  if [[ $VERBOSE == "true" ]] ; then
    print "  compute_tasks = $compute_tasks"
    print "  nwt = $nwt"
    print "  mwtpn = $mwtpn"
    print "  mctpn = $mctpn"
  fi

SetTaskGeometry $compute_tasks $nwt $mwtpn $mctpn

