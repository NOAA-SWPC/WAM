       subroutine if_shuff(global_lats_a_old,lats_nodes_a_old,
     .         global_lats_a,lats_nodes_a,timesum_a,kdt,ifshuff,
     .         shuffle_overhead)

c
      use gfs_dyn_resol_def
      use gfs_dyn_layout1
      use gfs_dyn_gg_def
      use gfs_dyn_vert_def
      use gfs_dyn_date_def
      use gfs_dyn_mpi_def
      implicit none
c
      integer               lats_nodes_a(nodes)
      integer              global_lats_a(latg)
c
      integer               lats_nodes_a_old(nodes)
      integer              global_lats_a_old(latg)
      integer kdt
      logical ifshuff
      real(kind=kind_evod) timesum_a(latg)
      real timesumold,timesumnew,timeold,timenew
      real shuffle_overhead
      integer node,ilat
      integer lat_old,lat_new
c
      timesumold = 0.0
      timeold = 0.0
      timesumnew = 0.0
      timenew = 0.0
      lat_old = 0
      lat_new = 0
      do node=1,nodes
cmy old
      timesumold = 0.0
      timesumnew = 0.0
         do ilat=1,lats_nodes_a_old(node)
            lat_old = lat_old + 1
       timesumold=timesumold + timesum_a(global_lats_a_old(lat_old))
c$$$         if (node .eq. nodes) print*,' old node, lat time = ',
c$$$     .     node,timesum_a(global_lats_a_old(lat_old)),
c$$$     .     global_lats_a_old(lat_old),timesumold
         enddo
           print*,' for node  timesumold = ',node,timesumold
 
          timeold = max(timeold,timesumold)
cmy new
         do ilat=1,lats_nodes_a(node)
            lat_new = lat_new + 1
       timesumnew=timesumnew + timesum_a(global_lats_a(lat_new))
c$$$         if (node .eq. nodes) print*,' new node, lat time = ',
c$$$     .     node,timesum_a(global_lats_a(lat_new)),
c$$$     .     global_lats_a(lat_new),timesumnew
         enddo
         print*,' for node timesumnew = ',node,timesumnew
          timenew = max(timenew,timesumnew)
      enddo
c
      if (timenew+shuffle_overhead .le. timeold) ifshuff = .true.
 
      print*,' from if_shuff kdt,new,old,overhead = ',
     .  kdt,timenew,timeold,shuffle_overhead
 
      return
      end
