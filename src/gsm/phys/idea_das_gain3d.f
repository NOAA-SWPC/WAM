!
      MODULE IDEA_DAS_PE
!
!       use resol_def, only : lonr, latr, levs !lonf, latg
       use module_gfs_machine, only : kind_grid
       real(kind=kind_grid),
     &  allocatable, dimension(:, :, :) :: wT, wO3, wOP, wP3d
       real(kind=kind_grid),allocatable :: T_das(:,:,:), O3_das(:,:,:)
       real(kind=kind_grid),allocatable :: P_das(:,:,:), OP_das(:,:,:)
       real(kind=kind_grid),allocatable :: Tmp3d(:,:,:)
!
      END MODULE IDEA_DAS_PE
!
      SUBROUTINE idea_init_das
      use mpi_def,     ONLY: liope, info, mpi_comm_all, 
     &                 mc_comp, mpi_comm_null
      use IDEA_DAS_PE           
      use resol_def, only : lonr, latr, levs
      use layout1,only : n1 => lats_node_r, me
      IMPLICIT NONE

      integer :: ierr
!
        allocate(wp3d(lonr, n1, levs),  wT(lonr,n1,levs))
        allocate(wOP(lonr,  n1, levs), wO3(lonr,n1,levs))

           allocate(T_das(lonr, latr, levs))
           allocate(P_das(lonr, latr, levs))
           allocate(O3_das(lonr, latr, levs))
           allocate(OP_das(lonr, latr, levs)) 
           allocate(Tmp3d(lonr, latr, levs))       
!  
       call mpi_barrier(mpi_comm_all,ierr) 
!
      END SUBROUTINE idea_init_das
!
!
      SUBROUTINE idea_cp4_omf(grid_fld, ind_lats, ind_lons)
!
       use layout1, only : lats_node_r, lats_node_r_max,  me
       use resol_def, only : lonr, latr, levs             !lonf, latg
       use gfs_physics_gridgr_mod, only : Grid_Var_Data
      use mpi_def,     ONLY: liope, info, mpi_comm_all, 
     &                 mc_comp, mpi_comm_null
       use IDEA_DAS_PE, only : wt, wop, wo3, wp3d 
!  
       IMPLICIT NONE
          TYPE(Grid_Var_Data)     :: grid_fld
          integer, dimension(lats_node_r) :: ind_lats, ind_lons
          integer :: n1, n2, lons_lat 
          integer ::  i,j,k, ierr       
  
          n1 = lonr
          n2 = lats_node_r

      do k=1, levs
!
!$omp parallel do private(i,j)
      do j=1,n2
         lons_lat=ind_lons(j)
        do i=1,lons_lat
        wt(i,j,k)=grid_fld%T(i, j,k)
        wO3(i,j,k)=grid_fld%tracers(2)%flds(i,j,k)
        wOP(i,j,k)=grid_fld%tracers(4)%flds(i,j,k)
        wp3d(i,j,k)=grid_fld%p(i, j,k)

!        wu(i,j,k)=grid_fld%U(i, j,k)
!        wv(i,j,k)=grid_fld%V(i, j,k)

!        wQ(i,j,k)=grid_fld%tracers(1)%flds(i,j,k)

!        wO2(i,j,k)=grid_fld%tracers(5)%flds(i,j,k)
!        wdp(i,j,k)=grid_fld%dp(i, j,k)

        enddo
      enddo

       call mpi_barrier(mpi_comm_all,ierr) 

      ENDDO ! levels
      END SUBROUTINE idea_cp4_omf
!
      SUBROUTINE idea_cp4_anal(grid_fld,ind_lats, ind_lons)
       use layout1, only : lats_node_r, lats_node_r_max,  me
       use resol_def, only : lonr, latr, levs             !lonf, latg
       use gfs_physics_gridgr_mod, only : Grid_Var_Data
       use IDEA_DAS_PE, only : wt, wo3, wop
       use mpi_def,     ONLY: liope, info, mpi_comm_all, 
     &                 mc_comp, mpi_comm_null      
!  
       IMPLICIT NONE
!
          TYPE(Grid_Var_Data)     :: grid_fld
          integer, dimension(lats_node_r) :: ind_lats, ind_lons
          integer :: n1, n2, lons_lat
          integer ::  i,j,k, ierr   
          
  
          n1 = lonr
          n2 = lats_node_r

!$omp parallel do private(i,j,k)
      do k=1, levs
!
      do j=1,n2
          lons_lat=ind_lons(j)
        do i=1,lons_lat
        grid_fld%T(i, j,k) = wt(i,j,k)
        grid_fld%tracers(2)%flds(i,j,k)=wO3(i,j,k)
        grid_fld%tracers(4)%flds(i,j,k)=wOP(i,j,k)
!        wp3d(i,j,k)=grid_fld%p(i, j,k)
!        wu(i,j,k)=grid_fld%U(i, j,k)
!        wv(i,j,k)=grid_fld%V(i, j,k)

        enddo
      enddo

      call mpi_barrier(mpi_comm_all,ierr) 

      ENDDO ! levels
      END SUBROUTINE idea_cp4_anal
!
       SUBROUTINE idea_das_window(global_lats_r,  lonsperlar, 
     &              fhour, Curr_NC_WAMDAY, tw1, tw2)
!
! Prototype from wam_nc_output16.f
!
        use resol_def,   ONLY: latr, levs, levp1, lonr
        use layout1,     ONLY: me, nodes, ipt_lats_node_r,
     &  lats_node_r
        use mpi_def,     ONLY: liope, info, mpi_comm_all, 
     &    mpi_r_io_r, mc_comp, mpi_comm_null

        USE machine,   ONLY: kind_io4, kind_io8
!
        USE IDEA_DAS_PE, only :  T_das, wT
        USE IDEA_DAS_PE, only : OP_das, wOP
        USE IDEA_DAS_PE, only : O3_das, wO3
        USE IDEA_DAS_PE, only :  P_das, wP3d
        USE IDEA_DAS_PE, only : Tmp3d
!        use  idea_ncout_phys
!        use gg_def,         only : colrad_r      ! latr
        use coordinate_def, only : ak5,bk5       ! levs+1
        use netcdf
!=======================================================
          implicit none
!
          real, intent(in)        ::  tw1, tw2
          real, intent(in)        :: fhour
          integer, intent(in)     :: Curr_NC_WAMDAY  
!        
!
       integer, intent(in), dimension(latr)::global_lats_r, lonsperlar
!
!locals
!
          integer      :: IOPROC, node, lat, ierr

          real         :: Tmp_t(lonr,  latr) 
          real         :: Tmp_op(lonr, latr)
          real         :: Tmp_o3(lonr, latr)
          real         :: Buffo(lonr, lats_node_r)

          integer      :: nelements 
          integer      :: i,j, k
!=======================================================
! put names in module
!
!=======================================================

           IOPROC  = nodes-1   !0

!       call mpi_barrier(mpi_comm_all,ierr)
!
! Global 3D-fields
!
!
       if (me == 0) then
        print *, nodes, ioproc, ' VAY-das-ioproc'
        print *, tw1, tw2,      ' VAY-das-ioproc'
        print *, Curr_NC_WAMDAY, ' VAY-das-WAMDAY'
!        print *, 'VAY-b-SGAIN-Wt ', maxval(Wt), minval(Wt)    
       endif

!  
!        call IDEA_DAS_GAIN3D
!     & (ioproc, Tmp3d, levs, wT, GLOBAL_LATS_R,LONSPERLAR)      

       call IDEA_DAS_GAIN3D
     & (ioproc, T_das, levs, wT, GLOBAL_LATS_R,LONSPERLAR)

      call mpi_barrier(mpi_comm_all,ierr)

!
       call IDEA_DAS_GAIN3D
     & (ioproc, P_das, levs, wP3d, GLOBAL_LATS_R,LONSPERLAR)
      call mpi_barrier(mpi_comm_all,ierr)

       call IDEA_DAS_GAIN3D
     & (ioproc, O3_das, levs, wO3, GLOBAL_LATS_R,LONSPERLAR)
      call mpi_barrier(mpi_comm_all,ierr)

       call IDEA_DAS_GAIN3D
     & (ioproc, OP_das, levs, wOP, GLOBAL_LATS_R,LONSPERLAR)
      call mpi_barrier(mpi_comm_all,ierr)
!
      if (me == ioproc) then
!
! T_DAS collected on ioproc
!
       print *, 'VAY-a-SG_TW1=', tw1,  ' tw2 =',tw2
       print *, 'VAY-a-SG_TD ',     maxval(T_das), minval(T_das), me  
       print *, 'VAY-a-SGAIN_OP ',  maxval(OP_das), minval(OP_das)
       print *, 'VAY-a-SGAIN_O3 ',  maxval(O3_das), minval(O3_das)
       print *, 'VAY-a-SGAIN_P3 ',  maxval(P_das), minval(P_das)   


  
      call make_das(T_das,P_das,O3_das,OP_das, lonr,latr,levs,ioproc,
     & tw1, tw2, fhour,  Curr_NC_WAMDAY)  
!    

      endif  

       nelements = size(T_das)
 
          
!
! We are sending 2D-slices of 3D-array at each level to all PEs from IOPROC
!
!
       do k=1, levs


          if(me == ioproc) then
            Tmp_t(:,:) =   T_das(:,:,K)
            Tmp_o3(:,:) =  O3_das(:,:,K)
            Tmp_op(:,:) =  OP_das(:,:,K)
          endif  

       call mpi_barrier(mpi_comm_all,ierr)   

       call mpi_bcast(Tmp_t, lonr*latr, MPI_R_IO_R, ioproc,MPI_COMM_ALL,info) 
       call mpi_barrier(mpi_comm_all,ierr)

       call mpi_bcast(Tmp_o3, lonr*latr, MPI_R_IO_R, ioproc,MPI_COMM_ALL,info) 
       call mpi_barrier(mpi_comm_all,ierr)

       call mpi_bcast(Tmp_op, lonr*latr, MPI_R_IO_R, ioproc,MPI_COMM_ALL,info) 
       call mpi_barrier(mpi_comm_all,ierr)

!       Tmp_o3 = reshape(O3_das, (/nelements/) )
!       call mpi_bcast(O3_das, nelements, MPI_R_IO_R, ioproc,MPI_COMM_ALL,info) 
!       call mpi_barrier(mpi_comm_all,ierr)

!       Tmp_op = reshape(OP_das, (/nelements/) )

!       call mpi_bcast(Op_das, nelements, MPI_R_IO_R, ioproc,MPI_COMM_ALL,info) 
!       call mpi_barrier(mpi_comm_all,ierr)

        call IDEA_DAS_SPLIT2d(Tmp_t,   wt(1,1,k), global_lats_r, lonsperlar)
        call IDEA_DAS_SPLIT2d(Tmp_o3, wo3(1,1,k), global_lats_r, lonsperlar) 
        call IDEA_DAS_SPLIT2d(Tmp_op, wop(1,1,k), global_lats_r, lonsperlar)  
   
       enddo


      
!       call mpi_bcast(Tmp3d, nelements, MPI_R_IO_R, ioproc,MPI_COMM_ALL,info) 
!       call mpi_barrier(mpi_comm_all,ierr)
!       endif


!      if(me .ne. ioproc) then 
!       T_das =  reshape(Tmp_t, (/lonr, latr, levs/))
!       O3_das = reshape(Tmp_o3, (/lonr, latr, levs/))
!       Op_das = reshape(Tmp_op, (/lonr, latr, levs/))
!       P_das = Tmp3d-T_das
!      print *, 'VAY-mpi-TOMF ', info, me, minval(P_das), maxval(P_das)
!      endif
! 
!      call mpi_bcast(O3_das,lonr*latr*levs,MPI_R_IO_R,ioproc,MPI_COMM_ALL,info) 
!      call mpi_bcast(OP_das,lonr*latr*levs,MPI_R_IO_R,ioproc,MPI_COMM_ALL,info) 
!      print *, 'VAY-mpi-T-das ', me, maxval(T_das), minval(T_das)
!      print *, 'VAY-mpi_bcast O3-das ', me, maxval(O3_das)
!      print *, 'VAY-mpi_bcast OP-das ', me, maxval(OP_das)
!      CALL mpi_send(tmp(1,1),illen*lonr,MPI_R_IO,ioproc,msgtag,MPI_COMM_ALL,info)

      RETURN
!
! split results back to wT, wO3, wOP
!
!$omp parallel do private(k, i,j,lat)
      DO k=1, levs
      do j=1,lats_node_r
          lat = global_lats_r(ipt_lats_node_r-1+j) 
       do i=1,lonr
         wt(i,j,k)   = T_das(i,lat,k)
         wo3(i,j,k)  = O3_das(i,lat,k)
         woP(i,j,k)  = OP_das(i,lat,k)
        enddo
      enddo
      ENDDO
!
       call mpi_barrier(mpi_comm_all,ierr)


       RETURN
   
!
!
      return

      END SUBROUTINE idea_das_window

!
      SUBROUTINE IDEA_DAS_GAIN3D
     & (ioproc, t3d, levs, t3dpe, GLOBAL_LATS_R,LONSPERLAR)
!--------------------------------------------------------
!
      use resol_def,   ONLY: latr, lonr
      use layout1,     ONLY: me, lats_node_r, lats_node_r_max,
     &                       ipt_lats_node_r, nodes
      use mpi_def,     ONLY: info, mpi_comm_all, liope, mpi_r_io,
     &                       stat,mpi_r_io_r
      USE machine,     ONLY: kind_io4, kind_io8
     
!--------------------------------------------------------
      IMPLICIT NONE
!     
!
      integer :: ioproc
      integer :: levs
      integer :: global_lats_r(latr), lonsperlar(latr)
!!
      real (kind=kind_io8) t3d(lonr, latr,         levs)     !out
      real (kind=kind_io8) t3dpe(lonr,lats_node_r, levs)

      real (kind=kind_io8) xl(lonr,lats_node_r)  ! for in-copy
      real (kind=kind_io8) xf(lonr,lats_node_r)  ! for in-copy

      real (kind=kind_io8)  x(lonr,latr)         !out

!    
   
      integer :: K, i,j
!
! t3dpe => xl => xf => X => t3d
!  
      DO k=1, levs
!
!$omp parallel do private(i,j)
      do j=1,lats_node_r
      do i=1,lonr
        xl(i, j) = t3dpe(i,j,K)
      enddo
      enddo
!
!subroutine unsplit2d_phys(ioproc,x,xl,global_lats_r)

       CALL wam_fill(xl,xf,global_lats_r,lonsperlar)
       CALL wam_unsplit2d(ioproc,x,xf,global_lats_r)
!
        do j=1,latr
           do i=1,lonr
              t3d(i, j,k) = X(i, j)     !fill-in is needed only for IOPROC
          enddo
        enddo

       enddo

       RETURN
       END  SUBROUTINE IDEA_DAS_GAIN3D
!
      SUBROUTINE IDEA_DAS_SPLIT3d
     & (Uanl, nz, me,  Ua, global_lats_r, lonsperlar)  
!
! here  NZ =levs
!    T_das, levs, me, Wt, GLOBAL_LATS_R,LONSPERLAR)
!
      USE machine,   ONLY: kind_io4, kind_io8
      use layout1,   ONLY: lats_node_r, ipt_lats_node_r 
      use resol_def, ONLY: latr, lonr, levs
!
      implicit none
      integer, intent(in) ::  nz, me
!                                                   
      integer, dimension(latr)       :: global_lats_r, lonsperlar
!
      real(kind=kind_io8), dimension(lonr,       latr, nz)   :: Uanl
      real(kind=kind_io8), dimension(lonr,lats_node_r, nz)   :: Ua

      real(kind=kind_io8)     ::   tmp(lonr,lats_node_r)
      real(kind=kind_io8)     ::  ftmp(lonr,lats_node_r)
      integer                                                :: iop 
      integer         :: i,j,k, lat

      DO k=1, nz
!
!-- get subdomain of data from 2D-global Uanl
!$omp parallel do private(i,j)
      do j=1,lats_node_r
          lat = global_lats_r(ipt_lats_node_r-1+j) 
       do i=1,lonr
         Tmp(i,j)  = Uanl(i,lat,k)
        enddo
      enddo

      CALL wam_spread(Tmp ,fTmp, global_lats_r,lonsperlar)

      do j=1,lats_node_r
       do i=1,lonr
        Ua(i,j,k) =fTmp(i,j)
        enddo
      enddo

      ENDDO
!
! On each processor it takes slice all "lons" for selected lats
!
      END SUBROUTINE IDEA_DAS_SPLIT3d


      SUBROUTINE IDEA_DAS_SPLIT2d(X, buffo, global_lats_r, LONSPERLAR)
!==================================================================== 
! 
!        X (full) => buffo(domain)
!
!       from "lons-PE => "lonr-PE"
!call split2d_rst(buff1,sfc_fld%slmsk,fieldsize,global_lats_r,LONSPERLAR)
!
!
!====================================================================
      USE machine,   ONLY: kind_io8
      use layout1,   ONLY: lats_node_r, ipt_lats_node_r 
      use resol_def, ONLY: latr, lonr, levs
!
      implicit none
!                                                   
      integer, dimension(latr)  :: global_lats_r, LONSPERLAR
!
      real(kind=kind_io8), dimension(lonr,       latr)   ::   X
      real(kind=kind_io8), dimension(lonr,lats_node_r)   :: buffo

      real(kind=kind_io8)     ::   tmp(lonr,lats_node_r)
!     real(kind=kind_io8)     ::  ftmp(lonr,lats_node_r)
      integer                                                :: iop 
      integer         :: i,j,k, lat
!
!-- get subdomain of data from 2D-global Uanl

      do j=1,lats_node_r
          lat = global_lats_r(ipt_lats_node_r-1+j) 
       do i=1,lonr
         Tmp(i,j)  =X(i,lat)
        enddo
      enddo
      CALL wam_spread(Tmp, buffo, global_lats_r,lonsperlar)   
    

! On each processor it takes slice all "lons" for selected lats
!
      END SUBROUTINE IDEA_DAS_SPLIT2d
!
!
      subroutine idea_split2d_1D(x,xl,fieldsize,global_lats_r,lonsperlar)
!
      use resol_def,      ONLY: latr, lonr
      use layout1,        ONLY: me, lats_node_r, ipt_lats_node_r, nodes
      use mpi_def,        ONLY: liope, mpi_comm_all, info,mpi_r_io_r
      USE machine,        ONLY: kind_ior, kind_io8
      implicit none
!!
!!
      integer,intent(in) :: fieldsize,global_lats_r(latr),
     &                      lonsperlar(latr)
      real(kind =kind_io8), intent(in)    :: x(fieldsize)
      real (kind=kind_io8),intent(inout)  :: xl(lonr,lats_node_r)
      integer j,lat,i,lon
!
!--- get subdomain of data
!
       do j=1,lats_node_r
         lat = global_lats_r(ipt_lats_node_r-1+j)
         if(lat /= 1) then
           lon = sum(lonsperlar(1:lat-1))
         else
           lon = 0
         endif
!
         do i=1,lonsperlar(lat)    !  i=1,lons
           xl(i,j) = X(lon+i)      !  storage by from i+ [lon = sum(lonsperlar(1:lat-1))]
                                   ! like grd%fld(i, j, klev) => X( lats_node_r*lonsperlar(lat) )
         enddo
       enddo
!
!
      return
      end subroutine idea_split2d_1D
