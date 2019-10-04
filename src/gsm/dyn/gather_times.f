      subroutine gather_times_a(lats_nodes_a,global_lats_a,
     .                        global_times_a,global_time_a)

      use gfs_dyn_resol_def
      use gfs_dyn_layout1
      use gfs_dyn_mpi_def
      implicit none
      
      integer lan,ierr,ilat,lat,node,nsend
      integer              global_lats_a(latg)
      integer              lats_nodes_a(nodes)
      integer icount
 
      real(kind=kind_evod) global_times_a(latg,nodes)
      real(kind=kind_evod) global_time_a(latg)
      real(kind=kind_io8) tmps(lats_node_a_max,nodes)
      real(kind=kind_io8) tmpr(lats_node_a_max,nodes)
 
      tmps = -55.55
      if (nodes.ne.1) then
        icount = 0
        do node=1,nodes
          do lan=1,lats_node_a
             icount = icount + 1
           lat=global_lats_a(ipt_lats_node_a-1+lan)
           tmps(lan,node)=global_times_a(lat,me+1)
c$$$           tmps(lan,node) = node * 1000 + lan
          enddo
        enddo
!sela   write(600+me,*) ' icount = ',icount
c
!sela   do node=1,nodes
!sela        write(600+me,*) ' node = ',node
!sela        write(600+me,620) (tmps(lan,node),lan=1,lats_node_a_max)
 620          format(5(e13.5,3x))
!sela   enddo
!!
        nsend=lats_node_a_max
c$$$        print*,' nsend = ',nsend
 
        call mpi_alltoall(tmps,nsend,MPI_A_DEF,
     x                     tmpr,nsend,MPI_A_DEF,
     x                     MC_COMP,ierr)
!!
 
!sela   do node=1,nodes
!sela     write(700+me,*) ' node = ',node
!sela     write(700+me,620) (tmpr(lan,node),lan=1,lats_node_a_max)
!sela   enddo
 
        ilat=1
        do node=1,nodes
          do lan=1,lats_nodes_a(node)
             lat=global_lats_a(ilat)
             global_time_a(lat)=tmpr(lan,node)
             ilat=ilat+1
          enddo
        enddo
      ELSE
        do lan=1,latg
           global_time_a(lan) = global_times_a(lan,1)
        enddo
      ENDIF
      return
      end
c
