c***********************************************************************
c
      SUBROUTINE EXCHA(lats_nodes_a,global_lats_a,X1,X2,Y1,Y2)
c
c***********************************************************************
c
      use gfs_dyn_resol_def
      use gfs_dyn_layout1
      use gfs_dyn_mpi_def
      implicit none
                                                                                
      integer n,i,j,ierr,ilat,lat,node,nsend
      integer              global_lats_a(latg)
      integer              lats_nodes_a(nodes)
      real(kind=kind_io8) X1(lats_node_a),X2(lats_node_a)
      real(kind=kind_io8) Y1(latg),Y2(latg)
      real(kind=kind_io8) tmps(2,lats_node_a_max,nodes)
      real(kind=kind_io8) tmpr(2,lats_node_a_max,nodes)
c
c@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
c
      if (nodes.ne.1) then
        do node=1,nodes
          do i=1,lats_node_a
           lat=global_lats_a(ipt_lats_node_a-1+i)
           tmps(1,i,node)=X1(I)
           tmps(2,i,node)=X2(I)
          enddo
        enddo
!!
        nsend=2*lats_node_a_max
        call mpi_alltoall(tmps,nsend,mpi_a_def,
     x                     tmpr,nsend,mpi_a_def,
     x                     MC_COMP,ierr)
!!
        ilat=1
        do node=1,nodes
          do i=1,lats_nodes_a(node)
             lat=global_lats_a(ilat)
             Y1(lat)=tmpr(1,i,node)
             Y2(lat)=tmpr(2,i,node)
             ilat=ilat+1
          enddo
        enddo
!!
      else
        y1=x1
        y2=x2
      endif
!!
      return
      end

