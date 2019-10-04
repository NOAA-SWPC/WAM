        
        subroutine GFS_simple_scatter ( global, local )

        use resol_def,   ONLY: lonr, latr
     &,                        scatter_lats, scatter_lons
        use layout1,     ONLY: me, lats_node_r, 
     &                         ipt_lats_node_r, nodes
c    intro of implicit none > >
	implicit none
	integer :: i
	integer :: j
	integer :: kmsk
	integer :: lat
c	^^^^^^^^^^^^
        real  :: global (lonr, latr)
        real  :: dist   (lonr, lats_node_r)
        real  :: local  (lonr, lats_node_r)

        real*4, allocatable :: global_r4(:,:)

!        integer ::  kmsk(lonr,latr)

!	integer ::  i, j, k, lat

!=== Begin here

        kmsk = 0

! r8 to r4
        allocate(global_r4(ubound(global,1),ubound(global,2)))
        global_r4(:,:) = global(:,:)

! fixed emission strength for debuging testing
!        do i =1, lonr
!        do j =1, latr
!         if ( global_r4(i,j)>0.)  global_r4(i,j) = 0.5
!        enddo
!        enddo

! scatter 
        do j=1,lats_node_r
           lat=scatter_lats(ipt_lats_node_r-1+j)
           do i=1,lonr
             dist(i,j)= global_r4(i,lat)
           enddo
        enddo


! full vs reduced grids
        local(:,:) = dist(:,:)
        CALL interpred_phys(1,kmsk,dist,local,scatter_lats,scatter_lons)

! deallocate
        deallocate(global_r4)

        end subroutine GFS_simple_scatter


