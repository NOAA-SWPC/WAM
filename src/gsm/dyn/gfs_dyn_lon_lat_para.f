      SUBROUTINE GFS_DYN_LONLAT_PARA(global_lats_a,XLON,XLAT,lonsperlat)
!
! 2009/10   sarah lu, compute lat/lon coordinate
!

c***********************************************************************
!
      USE GFS_DYN_MACHINE , ONLY : kind_grid

      use gfs_dyn_resol_def
      use gfs_dyn_layout1
      use gfs_dyn_gg_def
      use gfs_dyn_physcons, pi => con_pi

      use module_CPLFIELDS, ONLY: setupGauss2d,
     &                            wam2dmesh, MeshCreateReducedGaussian,
     &                            wamlevels

      implicit none

      integer i,j,lat
      integer                 lonsperlat(latg)
      real (kind=kind_grid) tpi,hpi,bphi
      real (kind=kind_grid) a
      PARAMETER (TPI=2.E0*PI,HPI=0.5E0*PI)
      integer              global_lats_a(latg)
      REAL(KIND=kind_dbl_prec) colrad_a_tmp(latg)
      real (kind=kind_grid) XLON(lonf,lats_node_a)
      real (kind=kind_grid) XLAT(lonf,lats_node_a)
!
      colrad_a_tmp(1:latg2)=colrad_a(1:latg2)
      do i=latg2+1,latg
         colrad_a_tmp(i)=colrad_a(latg+1-i)
      enddo
!
      xlon=0.
      xlat=0.
 
      DO j=1,lats_node_a
        lat = global_lats_a(ipt_lats_node_a-1+j)
        BPHI = TPI/lonsperlat(lat)
        if (lat.le.latg2) then
          DO i=1,lonsperlat(lat)
            XLON(I,J) = (i-1) * BPHI
            XLAT(I,J) = HPI - colrad_a_tmp(lat)
          ENDDO
        else
          DO i=1,lonsperlat(lat)
            XLON(I,J) =  (i-1) * BPHI
            XLAT(I,J) = colrad_a_tmp(lat)-HPI
          ENDDO
        endif
      ENDDO
      
      call setupGauss2d(lonsperlat(latg2), 2*latg2, pi, colrad_a, 
     & lats_node_a, global_lats_a, lonsperlat)
 
      wam2dmesh = MeshCreateReducedGaussian(ipt_lats_node_a,
     &  lats_node_a, lonsperlat, global_lats_a, colrad_a)
      wamlevels = levs

      RETURN
      END
 
