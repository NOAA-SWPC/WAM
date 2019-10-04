      subroutine getcon_lag(lats_nodes_a,global_lats_a,
     &                      lats_nodes_h,global_lats_h_sn,
     &                      lonsperlat,xhalo,yhalo)
      use gfs_dyn_resol_def,     only : jcap,latg,latg2,lonf
      use gfs_dyn_layout1,       only : me,nodes

      use gfs_dyn_gg_def,        only : colrad_a,sinlat_a
      use namelist_dynamics_def, only : shuff_lats_a
      use layout_lag   , only : ipt_lats_node_h,lat1s_h,lats_dim_h,
     &                          lats_node_h,lats_node_h_max,lon_dim_h
      implicit none
!
      integer     yhalo,xhalo
!
      integer, dimension(nodes) :: lats_nodes_a, lats_nodes_h
      integer, dimension(latg)  :: lonsperlat,   global_lats_a

      integer, dimension(latg+2*yhalo*nodes) :: global_lats_h_sn
!
      integer  i,j,l,n,lat,i1,i2,node,nodesio
      integer, dimension(latg+2*yhalo*nodes) :: global_lats_h_ns
!     logical  shuffled
!
      if (me == 0) print 100, jcap, me
100   format ('getcon_h jcap= ',i4,2x,'me=',i3)

      do lat = 1, latg2
         lonsperlat(latg+1-lat) = lonsperlat(lat)
      end do
      nodesio=nodes

!     print*,'con_h me,nodes,nodesio = ',me,nodes,nodesio

      shuff_lats_a = .false.

      if (me == 0) print*,' shuff_lats_a = ',shuff_lats_a

!     shuffled = shuff_lats_a
!     print*,'  else in if shuff_lats_a=',shuff_lats_a

      call setlats_lag(lats_nodes_a,global_lats_a,
     &               lats_nodes_h,global_lats_h_ns,yhalo)

!  reverse order for use in set_halos

      i1 = 1
      i2 = 0
      do n=1,nodes
         j  = 0
         i2 = i2 + lats_nodes_h(n)
         do i=i1,i2
            j = j + 1
            global_lats_h_sn(i) = global_lats_h_ns(i2+1-j)
         enddo
         i1 = i2 + 1
      enddo

!sela  write(6,*) ' getcon after setlats_h global_lats_a = '
!sela  write(6,830)     global_lats_a

!sela  write(6,*) ' getcon after setlats_h global_lats_h = '
!sela  write(6,830)     global_lats_h

 830   format(10(i4,1x))
      lats_dim_h = 0
      do node=1,nodes
         lats_dim_h = max(lats_dim_h, lats_nodes_h(node))
      enddo
      lats_node_h     = lats_nodes_h(me+1)
      lats_node_h_max = 0
      do i=1,nodes
        lats_node_h_max  = max(lats_node_h_max, lats_nodes_h(i))
      enddo
      ipt_lats_node_h = 1
      if ( me > 0 ) then
         do node=1,me
            ipt_lats_node_h = ipt_lats_node_h + lats_nodes_h(node)
         enddo
      endif
      do j=1,latg2
        sinlat_a(j) = cos(colrad_a(j))
      enddo
      do l=0,jcap
         do lat = 1, latg2
            if ( l <= min(jcap,lonsperlat(lat)/2) ) then
               lat1s_h(l) = lat
               go to 200
            endif
         end do
  200    continue
      end do
ccmr
!mjr  allocate ( lon_dims_h_coef(lats_node_h) )
!mjr  allocate ( lon_dims_h_grid(lats_node_h) )
ccmr
!mjr  do j=1,lats_node_h
!mjr     lat = global_lats_h(ipt_lats_node_h-1+j)
!mjr     if ( lonsperlat(lat) .eq. lonf ) then
!mjr        lon_dims_h_coef(j) = lonf +2 + 2*xhalo !even
!mjr        lon_dims_h_grid(j) = lonf +1 + 2*xhalo !even/odd
!mjr     else
ccmr        lon_dims_h_coef(j) = lonsperlat(lat) +2 + 2*xhalo !even
ccmr        lon_dims_h_grid(j) = lonsperlat(lat) +1 + 2*xhalo !even/odd
!mjr        lon_dims_h_coef(j) = lonf +2 + 2*xhalo !even
!mjr        lon_dims_h_grid(j) = lonf +1 + 2*xhalo !even/odd
!mjr     endif
!mjr  enddo
!mjr
      lon_dim_h = lonf + 1 + xhalo + xhalo !even/odd
!mjr
      return
      end
