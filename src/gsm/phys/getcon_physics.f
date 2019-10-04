      subroutine getcon_physics(n3,n4,
     &                          lats_nodes_r,global_lats_r,
     &                          lonsperlar,
     &                          lats_nodes_ext,global_lats_ext,
     &                          colat1,idrt)
!
!---------------------------------------------------------------------
! revision
! Dec 2014  Jun Wang  add slat_r,dlat_r

      use resol_def,            ONLY: latr, jintmx, nypt, lonrx, lonr,
     &                                latr2                             
      use layout1,              ONLY: me, nodes, lon_dims_r, 
     &                                lon_dims_ext, ipt_lats_node_r, 
     &                                ipt_lats_node_ext,
     &                                lats_node_r_max, lats_node_ext,
     &                                lats_node_r, lats_dim_r, 
     &                                lats_dim_ext
      use gg_def,               ONLY: colrad_r, wgt_r, wgtcs_r, rcs2_r, 
     &                                sinlat_r, coslat_r,slat_r,dlat_r
      use namelist_physics_def, ONLY: shuff_lats_r,semilag
!jw   use mpi_def,              ONLY: icolor, liope
      use machine,              ONLY: kind_dbl_prec, kind_evod
      implicit none
!!
      integer              i,j,k,l,lat,lev
      integer              n,n3,n4,idrt
!
      integer               lats_nodes_r(nodes)
      integer              global_lats_r(latr)
      integer                 lonsperlar(latr)
!
      integer                lats_nodes_ext(nodes)
      integer        global_lats_ext(latr+2*jintmx+2*nypt*(nodes-1))
!
      real(kind=kind_dbl_prec) ,allocatable:: colrad_dp(:), wgt_dp(:),
     &                                        wgtcs_dp(:),  rcs2_dp(:)
!
      integer              iprint,locl,node,nodesio
      integer              len_trie_ls_nod, len_trio_ls_nod
!
      integer              indlsev,jbasev,indlsod,jbasod
!
      integer gl_lats_index
      integer global_time_sort_index_r(latr)
      integer nodes_tmp
!
      include 'function_indlsev'
      include 'function_indlsod'
!
      real(kind=kind_evod) global_time_r(latr)
!
!     logical shuffled
!
      real(kind=kind_evod) colat1
!
      real(kind=kind_evod), parameter :: cons0    = 0.d0
!    &                                  ,cons0p5  = 0.5d0,
!    &                                   cons0p92 = 0.92d0,
!    &                                   cons1    = 1.d0
!
      iprint = 0
!
!     print 100, jcap, levs
100   format (1h0,'getcon physics ',i3,i3,' created january 2008')
!
      do lat = 1, latr2
         lonsperlar(latr+1-lat) = lonsperlar(lat)
      end do

!     write(0,*)' in getcon physics: lonsperlar ',lonsperlar
!
!jw
      idrt=4                                 !INTEGER DATA REPRESENTATION TYPE :4 Gaussian ,0:LATLON
!jw      if (liope) then
!jw         if (icolor.eq.2) then
!jw           nodesio=1
!jw         else
!jw           nodesio=nodes
!jw         endif
!jw      else

         nodesio = nodes

!jw      endif
!
!jw      if (nodesio .eq. 1 .and. nodes .eq. 1
!jw     .    .or. (nodes .eq. 2 .and. nodesio .eq. 1) ) then

      if (nodesio == 1 .and. nodes == 1) then
         shuff_lats_r = .false.
      endif
!!
!     print *,' getcon shuff_lats_r ',shuff_lats_r

      if (shuff_lats_r) then
 
        gl_lats_index = 0
        global_lats_r = -1

        do lat = 1,latr
         global_time_r(lat) = lonsperlar(lat)
        enddo
!
!my sort the lat times in descending order
!
        call sortrx(latr,-global_time_r,global_time_sort_index_r)

        if (iprint == 1)
     &   print *,' getcon_physics after sortrx for r index = ',
     &    global_time_sort_index_r
 
!my input lat time index in descending order
!my output global_lats_r and lats_nodes_r (gl_lats_index temp)
!
        gl_lats_index = 0
!
        nodes_tmp = nodes

!       if (liope .and. icolor .eq. 2) then
!         nodes_tmp           = nodes - 1
!         lats_nodes_r(nodes) = 0  ! gwvbugfix initialize lats_nodes_r(iotask)
!                                  ! which is not set in the loop below
!       endif

        do node=1,nodes_tmp
!          print *,' node gl_lats_index ',gl_lats_index
           call get_lats_node_r(node-1, global_lats_r,
     &                          lats_nodes_r(node),
     &                          gl_lats_index,global_time_sort_index_r,
     &                          iprint)
            if (me+1 == node .and. iprint == 1)
     &     print *,' node lats_nodes_r(node) ',lats_nodes_r(node)
        enddo
        call setlats_r_ext_shuff(lats_nodes_r,lats_nodes_ext,
     &           global_lats_r, global_lats_ext,iprint,lonsperlar)
      else

         if (semilag) then
           call setlats_r_slg(lats_nodes_r,global_lats_r,iprint,
     &                         lonsperlar)
         else
           call setlats_r(lats_nodes_r,lats_nodes_ext,global_lats_r,
     &                    global_lats_ext,iprint,lonsperlar)
         endif

      endif ! shuff_lats_r

!     write(0,*)' in getcon physics: lonsperlar ',lonsperlar
!    &,' me=',me,' global_lats_r=',global_lats_r
!
!      print *,' getcon physics: lats_nodes_r',lats_nodes_r
      iprint = 0
!
      lats_dim_r = 0
      do node=1,nodes
         lats_dim_r = max(lats_dim_r,lats_nodes_r(node))
      enddo
!      print *,' getcon physics: lats_dim_r',lats_dim_r
!
      lats_dim_ext = 0
      do node=1,nodes
             lats_dim_ext =
     &   max(lats_dim_ext, lats_nodes_ext(node), lats_nodes_r(node))
      enddo
!      print *,' getcon physics: lats_dim_ext',lats_dim_ext
!
      lats_node_r = lats_nodes_r(me+1)
!      print *,' getcon physics: lats_node_r',lats_node_r
!
      lats_node_ext = lats_nodes_ext(me+1)
!      print *,' getcon physics: lats_node_ext',lats_node_ext
!
      lats_node_r_max = 0
      do i=1,nodes
        lats_node_r_max = max(lats_node_r_max,lats_nodes_r(i))
      enddo
!      print *,' getcon physics: lats_node_r_max',lats_node_r_max
!
!
      ipt_lats_node_r = 1
      ipt_lats_node_ext = 1
  
!     if ( .not. shuffled .and. me .gt. 0 ) then
      if ( .not. shuff_lats_r .and. me .gt. 0 ) then
         do node=1,me
          ipt_lats_node_ext = ipt_lats_node_ext + lats_nodes_ext(node)
         enddo
      endif
!
      if ( me .gt. 0 ) then
         do node=1,me
            ipt_lats_node_r = ipt_lats_node_r + lats_nodes_r(node)
         enddo
      endif
!     print *,' getcon physics: ipt_lats_node_ext',ipt_lats_node_ext
!     print *,' getcon physics: ipt_lats_node_r',ipt_lats_node_r
!
!
      n3    = 51
      n4    = 52
!
!
      iprint = 0
!     if ( me == 0 ) iprint = 1
!
      if ( kind_evod == 8 ) then !------------------------------------

           call glats_physics(latr2,colrad_r(1:latr2),wgt_r,wgtcs_r,
     &                        rcs2_r,iprint)
!!
           colat1 = colrad_r(1)
!!
           do i=latr2+1,latr
              colrad_r(i) = colrad_r(latr+1-i)
           enddo
!
      else !------------------------------------------------------------
           allocate  ( colrad_dp(latr2) )
           allocate  (    wgt_dp(latr2) )
           allocate  (  wgtcs_dp(latr2) )
           allocate  (   rcs2_dp(latr2) )
!
           call glats_physics(latr2,colrad_dp,wgt_dp,wgtcs_dp,rcs2_dp,
     &                        iprint)
!!
           colat1=colrad_dp(1)
!!
           do i=1,latr2
              colrad_r(i) = colrad_dp(i)
                 wgt_r(i) =    wgt_dp(i)
               wgtcs_r(i) =  wgtcs_dp(i)
                rcs2_r(i) =   rcs2_dp(i)
           enddo
!
           do i=latr2+1,latr
              colrad_r(i) = colrad_dp(latr+1-i)
           enddo
!
           deallocate  ( colrad_dp )
           deallocate  (    wgt_dp )
           deallocate  (  wgtcs_dp )
           deallocate  (   rcs2_dp )
!
      endif !-----------------------------------------------------------
!
!     print *,' getcon physics: colrad_r',colrad_r
!
!
      do j=1,latr
        if (j <= latr2) then
          sinlat_r(j) = cos(colrad_r(j))
        else
          sinlat_r(j) = -cos(colrad_r(j))
        endif
        coslat_r(j) = sqrt(1. E 0 -sinlat_r(j)*sinlat_r(j))
      enddo
      if(.not.allocated(slat_r)) allocate(slat_r(lats_node_r_max))
      if(.not.allocated(dlat_r)) allocate(dlat_r(lats_node_r_max))
      do j=1,lats_node_r
         lat = global_lats_r(ipt_lats_node_r-1+j)
         slat_r(j) = sinlat_r(lat)
         if(lat < latr) then
           dlat_r(j) = abs(asin(sinlat_r(lat))-asin(sinlat_r(lat+1)))
         else
           dlat_r(j) = abs(-90.-asin(sinlat_r(lat)))
         endif
      enddo
!
      if(iprint == 1 ) then
        print *,' getcon physics: sinlat_r',sinlat_r
        print *,' getcon physics: slat_r',slat_r(1:lats_node_r)
        print *,' getcon physics: dlat_r',dlat_r(1:lats_node_r),
     &  'dim=',lats_node_r_max,lats_node_r
        print *,' getcon physics: coslat_r',coslat_r
      endif
!
!
      do j=1,lats_node_r
         lat = global_lats_r(ipt_lats_node_r-1+j)
         if ( lonsperlar(lat) .eq. lonr ) then
            lon_dims_r(j) = lonrx
         else
            lon_dims_r(j) = lonsperlar(lat) + 2
         endif
      enddo
!
!     print *,' getcon physics: lon_dims_r',lon_dims_r
!
!     if (.not. shuffled) then
!     do j=1,lats_node_ext
!        lat = global_lats_ext(ipt_lats_node_ext-1+j)
!        if ( lonsperlar(lat) .eq. lonr ) then
!           lon_dims_ext(j) = lonrx
!        else
!           lon_dims_ext(j) = lonsperlar(lat) + 1+2*nxpt+1
!        endif
!     enddo
!     endif
!
!     print *,' end of getcon physics '
      return
      end
