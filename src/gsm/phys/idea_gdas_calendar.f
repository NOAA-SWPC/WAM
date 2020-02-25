!23456
      module idea_gdas_calendar
!      
! contain appropriate names for hist_gdas hist_diag hist_wtide
!      
      character (len=100) Dirnc_case
!     $ ='/scratch3/NCEPDEV/swpc/scrub/Valery.Yudin/NC_DAS48/'
      character (len=20)  hist_gdas ! ='2016010100_inhr_'
      character (len=20)  hist_diag !='hist2016010100_diag_'
      character (len=20)  hist_wtide! ='hist2016010100_davr_'
      character (len=10)  str_idaygdas 
      character (len=9), parameter :: nml_gdas = 'gdasnc_in'   ! should be inside $RUNDIR
  
      integer           :: idaygdas
!      
      integer  :: irec_gdas
      integer  :: ihr
      integer  :: imm
      integer  :: idd     
      integer  :: iy4
      contains
!      
      subroutine init_idea_gdas_calendar(idate)
      
      implicit none
      
      integer idate(4)
      integer ierr
      
      namelist /gdasnc_nl/ Dirnc_case
      
      open(unit=155, file=trim(nml_gdas), status='old' )
      read(155, gdasnc_nl, iostat=ierr)  
      if (ierr .ne. 0) then
      print *, 'VAY-GDAS:-nml_gdas', trim(nml_gdas)
      print *, 'VAY-GDAS:-nml_gdas error ', ierr      
      write(6, gdasnc_nl)
      Dirnc_case='/scratch1/NCEPDEV/swpc/Svetlana.Karol/scrub/NC_DAS48/'
      endif 
      close(155)   

      iy4=idate(4)
      imm=idate(2)
      idd=idate(3)
      ihr=idate(1)  
      irec_gdas =1 
      idaygdas = iy4*1000000+imm*10000+idd*100+ihr
      write(STR_IDAYGDAS, fmt='(I10.10)') idaygdas
      hist_gdas ='hist'//STR_IDAYGDAS//'_inhr_'
      hist_diag ='hist'//STR_IDAYGDAS//'_diag_'
      hist_wtide='hist'//STR_IDAYGDAS//'_davr_'
!      if (ihr == 0 ) irec_gdas=1
      if (ihr == 6 )  irec_gdas=2
      if (ihr == 12 ) irec_gdas=3
      if (ihr == 18 ) irec_gdas=4                 
!       
! update     hist_gdas   hist_diag  hist_wtide  
! get/READ  Dirnc_case from Namelist ....
!
       print *, 'VAY-GDAS:', hist_gdas
       print *, 'VAY-GDAS:', hist_diag
       print *, 'VAY-GDAS:', hist_wtide
       print *, 'VAY-GDAS:', Dirnc_case    
!                    
      end subroutine init_idea_gdas_calendar  
      end module idea_gdas_calendar
