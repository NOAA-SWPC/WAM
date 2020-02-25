      SUBROUTINE IDEA_IAU_SPLIT3d
     & (Uanl,nz, me, Ua, global_lats_r, lonsperlar)
   
!
! here  NZ =/=levs
!
      USE machine,   ONLY: kind_io4, kind_io8
      use layout1,   ONLY: lats_node_r, ipt_lats_node_r 
      use resol_def, ONLY: latr, lonr, levs
!
      use idea_iau_gmao, only : nxa, nya, nza
      implicit none
      integer, intent(in) ::  nz,  me
!                                                   
      integer, dimension(latr)       :: global_lats_r, lonsperlar
!
      real(kind=kind_io8), dimension(nxa,nya,nz) :: Uanl
      real(kind=kind_io8), dimension(lonr,lats_node_r, nz)   :: Ua
      integer                                  :: iop 
      real(kind=kind_io8), allocatable, dimension(:,:) :: buff2, buffo 
!
      integer kmsk0(lonr,lats_node_r)
      integer         :: i,j,k, lat
!     print *, ' VAY-me-split', me, ipt_lats_node_r, lats_node_r

      allocate ( buff2(lonr,lats_node_r))
      allocate ( buffo(lonr,lats_node_r))
      DO k=1, nz
!$omp parallel do private(i,j)
      do j=1,lats_node_r
          lat = global_lats_r(ipt_lats_node_r-1+j) 
       do i=1,lonr
         buffo(i,j)  = Uanl(i,lat,k)
        enddo     
       enddo
!
         CALL wam_spread(buffo, buff2, global_lats_r,lonsperlar)     
!
        do j=1,lats_node_r
          do i=1,lonr
            Ua(i,j,k) = buff2(i,j)
          enddo
        enddo
      ENDDO
!
! On each processor it takes slice all "lons" for selected lats
!
      deallocate(buff2, buffo)
      RETURN

      END SUBROUTINE IDEA_IAU_SPLIT3D
!
      SUBROUTINE IDEA_IAU_SPLIT2d(global_lats_r, lonsperlar, 
     & me, Uanl, Ua)
      use layout1,   ONLY: lats_node_r, ipt_lats_node_r     
      use resol_def, ONLY: latr, lonr
      use idea_iau_gmao, only : nxa, nya
      implicit none
      integer, intent(in) :: me
!                                            
      integer, dimension(latr)       :: global_lats_r, lonsperlar
      real, dimension(nxa,  nya)     :: Uanl
      real, dimension(lonr, lats_node_r)       :: Ua
      integer                             :: iop =0
      real,   allocatable, dimension(:,:) :: buff2, buffo 
      integer   :: j, i, lat
!
       allocate ( buffo(lonr,lats_node_r))
       allocate ( buff2(lonr,lats_node_r))
!      call split2d_phys_iau(Uanl, buffo, global_lats_r, iop)

      print *, ' VAY-me-split2', me, ipt_lats_node_r, lats_node_r
!   
      do j=1,lats_node_r
        lat = global_lats_r(ipt_lats_node_r-1+j) 
       do i=1,lonr
         buffo(i,j) =  Uanl(i,lat)
        enddo
      enddo
!
         CALL wam_spread(buffo, buff2, global_lats_r,lonsperlar)     
!
        do j=1,lats_node_r
          do i=1,lonr
            Ua(i,j) = buff2(i,j)
          enddo
        enddo

      deallocate(buff2, buffo)  ! buffo, buff2
     
      RETURN
      end SUBROUTINE IDEA_IAU_SPLIT2D
!
!***********************************************************************
