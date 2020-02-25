       subroutine idea_writediag_wamday( global_lats_r,lonsperlar,
     &        fhour,idate,Curr_NC_WAMDAY)
!
! Prototype from wam_ncdiag_output17.f => GW-diagnostics
!
!

        use resol_def,   ONLY: latr, levs, levp1, lonr
        use layout1,     ONLY: me, nodes, lats_node_r
        use mpi_def,     ONLY: liope, info, mpi_comm_all, 
     &                                   mc_comp, mpi_comm_null
        USE machine,   ONLY: kind_io4, kind_io8
!
        use  idea_ncout_diag
        use gg_def,         only : colrad_r      ! latr
        use coordinate_def, only : ak5,bk5       ! levs+1
        use netcdf
        use idea_lat_gaus, only : lat_gaus_t62
!
	use idea_gdas_calendar, only : hist_diag, Dirnc_case



          implicit none

          real        :: fhour
          integer     :: idate(4)
          integer     :: Curr_NC_WAMDAY            !, Wam_daysteps

!
          integer, dimension(latr)   :: global_lats_r, lonsperlar
!
          integer         :: IOPROC, node, lat, ierr
!
      real, allocatable   :: tmp3d(:,:,:), tmp2d(:,:)
!=======================================================================
! need to be passed: hyam, hybm, hyai, hybi, lons_deg, lats_deg, levi
!=======================================================================
      character(len=8)    :: S8_DATE
      character(len=8)    :: file_wam='wamdiag_'
      character(len=21)   :: file_hist
      character(len=11)   :: File_day='20160125' 
      character(len=3)    :: Fend ='.nc'
      character(len=121)  :: file_hist_das
      character(len=100)  :: dir_case='/scratch3/NCEPDEV/swpc/scrub/Valery.Yudin/NC_DAS/'
!
!
! Coordinate arrays
!
!
      real, dimension(levs)   :: pmbm, hyam, hybm 
      real, dimension(levs+1) :: pmbi, hyai, hybi 
      real                    :: lonwam(lonr)
      real                    :: latwam(latr)
      real                    :: Pref
!
      integer    :: status, ncid, iernc    
      real       ::  time  
  
      integer :: NxDimID, NyDimID, NzDimID, NzIDimID, NtDimID
      integer :: scalDimID, Vid_sp0
      integer :: VidLat, VidLon, Vidlev, VidIlev
      integer :: VidHyam, VidHybm, VidHyai, VidHybi
      integer :: VidTime, VidDate, VidDsec
!
      integer :: VidHs, VidPs
      integer :: VidT, VidQ, VidV, VidU, VidZgkm
      integer :: Vidqo2, Vidqo3, Vidqo  

      integer :: Vid

      integer :: Datesec, Ymd

         integer :: start1(1), count1(1)
         integer :: start(4), count(4)
         integer :: nz, ny, nx
!
         integer :: start3(4), count3(4)
         integer :: ihr
!
         integer :: i, j,k
!
             ihr =1                        ! daily-record
          YMD = Curr_NC_WAMDAY
!
! define IOPROC
!
          IOPROC =nodes-1                  !like in SFC-RSTSR  wrtout_physics.f  
!
        call mpi_barrier(mpi_comm_all,ierr)
!      t3  = rtc()
!      call mpi_barrier(mpi_comm_all,ierr)
!      t4  = rtc()
!      tba = t4-t3


          nz =levs
          ny = latr
          nx = lonr

          start=(/1,1,1, ihr/)
          count=(/nz, ny, nx,1/)

          start1(1) =1
          count1(1) =1


      allocate ( tmp2d (lonr,latr), tmp3d (lonr,latr, levs) )

      if(me == ioproc) then   
!
! Create name of the HIST-file
!
    
      write(S8_DATE, fmt='(I8.8)') Curr_NC_WAMDAY

      print *, S8_DATE, ' NC_STR-date '
      File_hist=trim(file_wam)//trim(S8_DATE)//trim(Fend)

!      File_hist_das=trim(dir_case)//trim(file_wam)//trim(S8_DATE)//trim(Fend)
      File_hist_das=trim(dirnc_case)//trim(hist_diag)//trim(S8_DATE)//trim(Fend)      
      print *, File_hist_das, ' VAY DIAG-File_hist_das NC-file '
!
! open-File_hist
!23456
      ierNC=NF90_OPEN(trim(File_hist_das), nf90_write, ncid)
      status = nf90_create(trim(File_hist_das), nf90_clobber, ncid)
!
! Dimensions 
!
      status = nf90_def_dim(ncid, 'lon',  lonr,          NxDimID)
      status = nf90_def_dim(ncid, 'lat',  latr,          NyDimID)
      status = nf90_def_dim(ncid, 'time', NF90_UNLIMITED,NtDimID)
      status = nf90_def_dim(ncid, 'scalar', 1,         scalDimID)

      status = nf90_def_dim(ncid, 'lev',  levs,          NzDimID)
      status = nf90_def_dim(ncid, 'ilev', levs+1,        NzIDimID)
!
! Create vars and attributes .........
! Coordinates
!
      status = nf90_def_var(ncid, 'time',nf90_float, (/ NtDimID /), VidTime)
          status = nf90_put_att(ncid,VidTime, 'long_name', 'data-time')
          status = nf90_put_att(ncid,VidTime, 'units', 'days since 2012-01-01 00:00:00')
!
          status = nf90_def_var(ncid, 'date',     nf90_INT, (/ NtDimID /), VidDate)
          status = nf90_put_att(ncid,VidDate, 'long_name', 'current date (YYYYMMDD)')

          status = nf90_def_var(ncid, 'datesec',     nf90_INT, (/ NtDimID /), VidDsec)
          status = nf90_put_att(ncid,VidDsec, 'long_name', 'current seconds of current date')
          status = nf90_put_att(ncid,VidDsec, 'units', 'seconds')


!
          status = nf90_def_var(ncid, 'lat',     nf90_float, (/ NyDimID /),  VidLat)
          status = nf90_put_att(ncid,Vidlat, 'long_name', 'data-latitude')
          status = nf90_put_att(ncid,Vidlat, 'units', 'degrees_north')
!
          status = nf90_def_var(ncid, 'lon',     nf90_float, (/ NxDimID /), Vidlon)
          status = nf90_put_att(ncid,Vidlon, 'long_name', 'data-longitude')
          status = nf90_put_att(ncid,Vidlon, 'units', 'degrees_east')
!
          status = nf90_def_var(ncid, 'lev',     nf90_float, (/ NzDimID /), VidLev)
          status = nf90_put_att(ncid,Vidlev, 'units', 'level')
          status = nf90_put_att(ncid,Vidlev, 'standard_name', 'atmosphere_hybrid_sigma_pressure_coordinate')
          status = nf90_put_att(ncid,Vidlev, 'long_name', 'hybrid level at midpoint (1000*(A+B))')
!          status = nf90_put_att(ncid,Vidlev, 'units', 'mb')
          status = nf90_put_att(ncid,Vidlev, 'positive', 'down') !   ncdf_attput, id, ilev_id, 'positive', 'down'
          status = nf90_def_var(ncid, 'hyam',     nf90_float, (/ NzDimID /), VidHyam)
          status = nf90_put_att(ncid,VidHyam, 'long_name', 'hybrid level at midpoint A')
          status = nf90_put_att(ncid,VidHyam, 'units', ' dimensionless ')

          status = nf90_def_var(ncid, 'hybm',     nf90_float, (/ NzDimID /), VidHybm)
          status = nf90_put_att(ncid,VidHybm, 'long_name', 'hybrid level at midpoint B')
          status = nf90_put_att(ncid,VidHybm, 'units', ' dimensionless ')
!
          status = nf90_def_var(ncid, 'ilev',     nf90_float, (/ NzIDimID /), VidiLev)
          status = nf90_put_att(ncid,VidIlev, 'standard_name', 'atmosphere_hybrid_sigma_pressure_coordinate')
          status = nf90_put_att(ncid,VidIlev, 'long_name', 'hybrid level at interface (1000*(A+B))')
          status = nf90_put_att(ncid,VidIlev, 'units', 'lev')
          status = nf90_put_att(ncid,VidIlev, 'positive', 'down')
          status = nf90_def_var(ncid, 'hyai',     nf90_float, (/ NzIDimID /), VidHyaI)
          status = nf90_put_att(ncid,VidHyam, 'long_name', 'hybrid level at interface Ai')
          status = nf90_put_att(ncid,VidHyam, 'units', ' dimensionless ')

          status = nf90_def_var(ncid, 'hybi',     nf90_float, (/ NzIDimID /), VidHybI)
          status = nf90_put_att(ncid,VidHybI, 'long_name', 'hybrid level at interface Bi')
          status = nf90_put_att(ncid,VidHybI, 'units', ' dimensionless ')
!
          status = nf90_def_var(ncid, 'P0',     nf90_float, (/ scalDimID /), Vid_sp0)
          status = nf90_put_att(ncid,Vid_sp0, 'long_name', ' Reference surface pressure')
          status = nf90_put_att(ncid,Vid_sp0, 'units', ' Pa ')

!          status = nf90_def_var(ncid, 'day_obsed',   nf90_double, (/ NtDimID /),  VIDdate)
!============================================
!=========================================== 
!23456
      status = nf90_def_var(ncid, 'HS', nf90_float, (/ NxDimId, NyDimId, NtDimID /), VidHs)
      status = nf90_put_att(ncid,VidT, 'long_name', ' Surface elevation ')
      status = nf90_put_att(ncid,VidT, 'units', 'm')
!
      status = nf90_def_var(ncid, 'PS', nf90_float, (/ NxDimId, NyDimId, NtDimID /), VidHs)
      status = nf90_put_att(ncid,VidT, 'long_name', ' Surface pressure ')
      status = nf90_put_att(ncid,VidT, 'units', ' Pa ')
!
      status = nf90_def_var(ncid, 'U', nf90_float, (/ NxDimId, NyDimId, NzDimID,NtDimID /), VidU)
      status = nf90_put_att(ncid,VidU, 'long_name', ' Zonal wind ')
      status = nf90_put_att(ncid,VidU, 'units', 'm/s')
!
      status = nf90_def_var(ncid, 'V', nf90_float, (/ NxDimId, NyDimId, NzDimID,NtDimID /), VidV)
      status = nf90_put_att(ncid,VidV, 'long_name', ' Meridional wind ')
      status = nf90_put_att(ncid,VidV, 'units', 'm/s' )
!
      status = nf90_def_var(ncid, 'T', nf90_float, (/ NxDimId, NyDimId, NzDimID,NtDimID /), VidT)
      status = nf90_put_att(ncid,VidT, 'long_name', ' Temperature')
      status = nf90_put_att(ncid,VidT, 'units', 'K') 
!
      status = nf90_def_var(ncid, 'P3D', nf90_float, (/ NxDimId, NyDimId, NzDimID,NtDimID /), VidT)
      status = nf90_put_att(ncid,VidT, 'long_name', ' Pressure-3D')
      status = nf90_put_att(ncid,VidT, 'units', 'Pa')
!
      status = nf90_def_var(ncid, 'DP', nf90_float, (/ NxDimId, NyDimId, NzDimID,NtDimID /), VidT)
      status = nf90_put_att(ncid,VidT, 'long_name', ' P-thickness')
      status = nf90_put_att(ncid,VidT, 'units', 'Pa')
!
!      status = nf90_def_var(ncid, 'AXGW', nf90_float, (/ NxDimId, NyDimId, NzDimID,NtDimID /), VidQ)
!      status = nf90_put_att(ncid,VidQ, 'long_name', ' AX-GW, W-E GW Acc-n')
!      status = nf90_put_att(ncid,VidQ, 'units', ' m/s2 ') 

!      status = nf90_def_var(ncid, 'AYGW', nf90_float, (/ NxDimId, NyDimId, NzDimID,NtDimID /), VidQ)
!      status = nf90_put_att(ncid,VidQ, 'long_name', ' AX-GW, S-N GW Acc-n')
!      status = nf90_put_att(ncid,VidQ, 'units', ' m/s2 ') 

!      status = nf90_def_var(ncid, 'EPSGW', nf90_float, (/ NxDimId, NyDimId, NzDimID,NtDimID /), VidQ)
!      status = nf90_put_att(ncid,VidQ, 'long_name', ' Eps-GW Heat')
!      status = nf90_put_att(ncid,VidQ, 'units', ' K/day ') 

!      status = nf90_def_var(ncid, 'KEDGW', nf90_float, (/ NxDimId, NyDimId, NzDimID,NtDimID /), VidQ)
!      status = nf90_put_att(ncid,VidQ, 'long_name', ' GW Eddy Diffusion, m2/s')
!      status = nf90_put_att(ncid,VidQ, 'units', ' m2/s ') 
!
      status = nf90_def_var(ncid, 'ZGKM', nf90_float, (/ NxDimId, NyDimId, NzDimID,NtDimID /), VidZgkm)
      status = nf90_put_att(ncid,VidZgkm, 'long_name', ' Geo-Height')
      status = nf90_put_att(ncid,VidZgkm, 'units', ' km ') 
!
      status = nf90_enddef(ncid)

!========================================
! Put Ymd and Model Coordiantes
!========================================

      time = 24.0           ! daily mean
      Datesec = ihr * 3600     
!
! Grids:     see glats_physics.f
!
      do i=1, lonr
         lonwam(i) = 0.+ (i-1)*360./float(lonr)
      enddo

      do i=1, latr
!         latwam(i) = 90.- (i-1)*180./float(latr-1)
         latwam(i)= lat_gaus_t62(i)     !colrad_r(i)*45./atan(1.)      ! Gaussian Lats
      enddo
!
! need to define hyai/hybi on edges  with zeroes at the top
! VG-read/write in the dynamical core
!
      do k=1, levs+1
        hyai(k)=ak5(k)*1.e-2         !*1000/(Ps = 1.e5) = 1.e-2
        hybi(k)=bk5(k)
        pmbi(k) =1.e5*(hyai(k)+hybi(k))
      enddo

      do k=1, levs
        pmbm(k)=0.5*(pmbi(k)+pmbi(k+1))
        hyam(k)=0.5*(hyai(k)+hyai(k+1))
        hybm(k)=0.5*(hybi(k)+hybi(k+1))
      enddo
!
! Grids completed
!
      status = nf90_put_var(ncid,  VidDate,   YMD)      !, start = start1,  count = count1)
      status = nf90_put_var(ncid,  VidDsec,   Datesec) !
      status = nf90_put_var(ncid,  VidTime,   Time)
 
      status = nf90_put_var(ncid, VidLon, LonWam)
      status = nf90_put_var(ncid, VidLat, LatWam)

      status = nf90_put_var(ncid, VidTime, Time)
      status = nf90_put_var(ncid, VidLev, Pmbm)
      status = nf90_put_var(ncid, VidHyam, hyam)
      status = nf90_put_var(ncid, VidHybm, hybm)

      status = nf90_put_var(ncid, VidiLev, Pmbi)
      status = nf90_put_var(ncid, VidHyai, hyai)
      status = nf90_put_var(ncid, VidHybi, hybi)
      Pref = 1.e5
      status = nf90_put_var(ncid, Vid_sp0, Pref)


      ENDIF     ! (IOPROC)
!=================================================
! Now Gain 2D/3D and write on ME=IOPROC
!=================================================

        call IDEA_NCOUT_GAIN3D
     & (ioproc, tmp2d, 1, d_wPs, GLOBAL_LATS_R,LONSPERLAR)
!
       call mpi_barrier(mpi_comm_all,ierr)
      if(me == ioproc) then    
        iernc=nf90_inq_varid( ncid, 'PS', vid )
        iernc= nf90_put_var( ncid, vid, tmp2d)
      endif
        call IDEA_NCOUT_GAIN3D
     & (ioproc, tmp2d, 1, d_wPhis, GLOBAL_LATS_R,LONSPERLAR)
      if(me == ioproc) then    
        iernc=nf90_inq_varid( ncid, 'HS', vid )
        iernc= nf90_put_var( ncid, vid, tmp2d)
      endif
!U
      call IDEA_NCOUT_GAIN3D
     & (ioproc, tmp3d, levs, d_wU, GLOBAL_LATS_R,LONSPERLAR)

       call mpi_barrier(mpi_comm_all,ierr)

      if(me == ioproc) then 
        iernc=nf90_inq_varid( ncid, 'U', vid )
        iernc= nf90_put_var( ncid, vid, tmp3d)
       endif
!V
      call IDEA_NCOUT_GAIN3D
     & (ioproc, tmp3d, levs, d_wV, GLOBAL_LATS_R,LONSPERLAR)

       call mpi_barrier(mpi_comm_all,ierr)
      if(me == ioproc) then 
        iernc=nf90_inq_varid( ncid, 'V', vid )
        iernc= nf90_put_var( ncid, vid, tmp3d)
      endif
!T
      call IDEA_NCOUT_GAIN3D
     & (ioproc, tmp3d, levs, d_wT, GLOBAL_LATS_R,LONSPERLAR)
      call mpi_barrier(mpi_comm_all,ierr)
      if(me == ioproc) then 
        iernc=nf90_inq_varid( ncid, 'T', vid )
        iernc= nf90_put_var( ncid, vid, tmp3d)
      endif
!       print *, 'VAY-NCHIST-T ', maxval(tmp3d), minval(tmp3d)
!
! ADD P3D, DP, ZPHIL + TIDES
!
!P3d
      call IDEA_NCOUT_GAIN3D
     & (ioproc, tmp3d, levs, d_P3d, GLOBAL_LATS_R,LONSPERLAR)
      call mpi_barrier(mpi_comm_all,ierr)
      if(me == ioproc) then 
        iernc=nf90_inq_varid( ncid, 'P3D', vid )
        iernc= nf90_put_var( ncid, vid, tmp3d)
      endif
!DP
      call IDEA_NCOUT_GAIN3D
     & (ioproc, tmp3d, levs, d_Dp, GLOBAL_LATS_R,LONSPERLAR)
      call mpi_barrier(mpi_comm_all,ierr)

      if(me == ioproc) then 
        iernc=nf90_inq_varid( ncid, 'DP', vid )
        iernc= nf90_put_var( ncid, vid, tmp3d)
!
! Need Compute 'ZGKM' geopotential
!
        
        iernc=nf90_inq_varid( ncid, 'ZGKM', vid )
        iernc= nf90_put_var( ncid, vid, tmp3d )
      endif
!=================================================================
!GW-diagnostics:   wAxz, wAyz, wEps, wKed 
      goto 555

      call IDEA_NCOUT_GAIN3D
     & (ioproc, tmp3d, levs, wAxz, GLOBAL_LATS_R,LONSPERLAR)
      call mpi_barrier(mpi_comm_all,ierr)
      if(me == ioproc) then 
        iernc=nf90_inq_varid( ncid, 'AXGW', vid )
        iernc= nf90_put_var( ncid, vid, tmp3d)
      endif
!AY
       call IDEA_NCOUT_GAIN3D
     & (ioproc, tmp3d, levs, wAyz, GLOBAL_LATS_R,LONSPERLAR)
      call mpi_barrier(mpi_comm_all,ierr)
      if(me == ioproc) then 
        iernc=nf90_inq_varid( ncid, 'AYGW', vid )
        iernc= nf90_put_var( ncid, vid, tmp3d)
      endif     
!Eps
       call IDEA_NCOUT_GAIN3D
     & (ioproc, tmp3d, levs, wEps,  GLOBAL_LATS_R,LONSPERLAR)
      call mpi_barrier(mpi_comm_all,ierr)
      if(me == ioproc) then 
        iernc=nf90_inq_varid( ncid, 'EPSGW', vid )
        iernc= nf90_put_var( ncid, vid, tmp3d )
      endif
!Ked
       call IDEA_NCOUT_GAIN3D
     & (ioproc, tmp3d, levs, wKed, GLOBAL_LATS_R,LONSPERLAR)
      call mpi_barrier(mpi_comm_all,ierr)
      if(me == ioproc) then 
        iernc=nf90_inq_varid( ncid, 'KEDGW', vid )
        iernc= nf90_put_var( ncid, vid, tmp3d )
      ENDIF
!
! close file
!
 555    continue 
        if(me == ioproc) then
           status = nf90_close(ncid)

!
! data TYPES: NF90_BYTE, NF90_CHAR, NF90_SHORT, NF90_INT, NF90_FLOAT, and NF90_DOUBLE
!
           print *, ' NO-collects for NC-file '
        endif    ! IOPROC

        deallocate ( tmp2d, tmp3d )


         RETURN

 
!==============================================================================
          end subroutine idea_writediag_wamday
!         
