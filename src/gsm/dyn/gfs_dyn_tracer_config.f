!
!! ! Module: gfs_dyn_tracer_config
!
! ! Description: gfs dynamics tracer configuration module
!
! ! Revision history:
!   Aug 08 2009   Sarah Lu, initial code
!   Aug 09 2009   Sarah Lu, add ntrac_chem, ntrac_met
!   Aug 10 2009   Sarah Lu, gfs_dyn_tracer is determined from ChemRegistry
!   Oct 16 2009   Sarah Lu, remove ChemRegistry; hardwire tracer specification 
!                           for testing; port to the latest trunk
!   Nov 13, 2009  Weiyu Yang, modified for the ensemble GEFS code.
!   Nov 19 2009   Sarah Lu, chem tracer specified from ChemRegistry
!   Feb 09 2009   Sarah Lu, ri/cpi added to gfs_dyn_tracer_type
!   Aug 17 2010   Sarah Lu, remove debug print
!   Aug 30 2010   Sarah Lu, set glbsum default as F
!   Aug 08 2011   Jun Wang, compile gocart only when running GOCART
!   Sep 16 2011   Sarah Lu, revise chem tracer initialization 
!   Jan 06 2012   Henry Juang, revise meteorological tracers
!          2015   S Moorthi    Adding tke
! -------------------------------------------------------------------------
!
      module gfs_dyn_tracer_config
      use gfs_dyn_machine ,     only : kind_grid
      use gfs_dyn_tracer_const, only : cpi,ri

      implicit none
      SAVE
!
! tracer specification
!
      type    gfs_dyn_tracer_type
        character*20,         pointer      :: chem_name(:)   ! chem_tracer name
        character*20,         pointer      :: vname(:, :)    ! variable name
        real(kind=kind_grid), pointer      :: ri(:)
        real(kind=kind_grid), pointer      :: cpi(:)
        integer                  :: ntrac, ntrac_met, ntrac_chem
        logical                  :: doing_DU, doing_SU, doing_SS
     &,                             doing_OC, doing_BC, doing_GOCART
      endtype gfs_dyn_tracer_type

      type (gfs_dyn_tracer_type), save     ::  gfs_dyn_tracer
!
! misc tracer options
!
      logical, save                        :: glbsum  = .false.
!

! --- public interface
      public     tracer_config_init

      contains

! -------------------------------------------------------------------   
      subroutine tracer_config_init (ntrac,
     &                               ntoz,ntcw,ncld,ntke,me)

!  
!  This subprogram sets up gfs_dyn_tracer
! 
      implicit none
! input
      integer, intent(in)       ::  me, ntoz,ntcw,ncld,ntke
! output
!      type (gfs_dyn_tracer_type), intent(out)    ::  gfs_dyn_tracer
! input/output
      integer, intent(inout)    :: ntrac
! local
      integer                   :: i, j, status, ierr
      character*20              :: rgname

! initialize ntrac_chem (the default is no chemistry)
      gfs_dyn_tracer%ntrac_chem = 0
      gfs_dyn_tracer%doing_GOCART = .false.

! initialize chem tracers
      call dyn_gocart_tracer_config(me)
!     call dyn_gocart_tracer_config(gfs_dyn_tracer,me)

! ntrac_met = number of met tracers
!hmhj if ( ntoz < ntcw ) then                       
!hmhj   gfs_dyn_tracer%ntrac_met = ntcw + ncld - 1   
!hmhj else                                                           
!hmhj   gfs_dyn_tracer%ntrac_met = ntoz                              
!hmhj endif                                          
!hmhj if ( gfs_dyn_tracer%ntrac_met /= ntrac ) then
!hmhj   print *,'LU_TRC: ERROR ! inconsistency in ntrac:',
!hmhj&           ntrac, gfs_dyn_tracer%ntrac_met
!hmhj   stop  222
!hmhj endif
! input ntrac is meteorological tracers
      gfs_dyn_tracer%ntrac_met = ntrac

! update ntrac = total number of tracers
      gfs_dyn_tracer%ntrac = gfs_dyn_tracer%ntrac_met +     
     &                       gfs_dyn_tracer%ntrac_chem
      ntrac = gfs_dyn_tracer%ntrac

      if(me==0) then
       print *, 'LU_TRC: ntrac_met =',gfs_dyn_tracer%ntrac_met
       print *, 'LU_TRC: ntrac_chem=',gfs_dyn_tracer%ntrac_chem
       print *, 'LU_TRC: ntrac     =',gfs_dyn_tracer%ntrac
      endif

! Set up tracer name, cpi, and ri
      if ( gfs_dyn_tracer%ntrac > 0 ) then      
       allocate(gfs_dyn_tracer%vname(ntrac, 5), stat=status)
           if( status .ne. 0 ) go to 999         
       allocate(gfs_dyn_tracer%ri(0:ntrac),  stat=status)
           if( status .ne. 0 ) go to 999
       allocate(gfs_dyn_tracer%cpi(0:ntrac), stat=status)
           if( status /= 0 ) go to 999

!--- fill in met tracers
      gfs_dyn_tracer%vname(1,    1) = 'spfh'   
      gfs_dyn_tracer%vname(1,    2) = 'spfh_q'   
      gfs_dyn_tracer%vname(1,    3) = 'spfh_m'   
      gfs_dyn_tracer%vname(1,    4) = 'spfh_q6'   
      gfs_dyn_tracer%vname(1,    5) = 'spfh_m6'   
      if(ntoz > 1) then
        gfs_dyn_tracer%vname(ntoz, 1) = 'o3mr'  
        gfs_dyn_tracer%vname(ntoz, 2) = 'o3mr_q'  
        gfs_dyn_tracer%vname(ntoz, 3) = 'o3mr_m'  
        gfs_dyn_tracer%vname(ntoz, 4) = 'o3mr_q6'  
        gfs_dyn_tracer%vname(ntoz, 5) = 'o3mr_m6'  
      endif
      if(ntcw > 1) then
        gfs_dyn_tracer%vname(ntcw, 1) = 'clwmr' 
        gfs_dyn_tracer%vname(ntcw, 2) = 'clwmr_q' 
        gfs_dyn_tracer%vname(ntcw, 3) = 'clwmr_m' 
        gfs_dyn_tracer%vname(ntcw, 4) = 'clwmr_q6' 
        gfs_dyn_tracer%vname(ntcw, 5) = 'clwmr_m6' 
      endif
      if(ntke > 1) then
        gfs_dyn_tracer%vname(ntke, 1) = 'tke'
        gfs_dyn_tracer%vname(ntke, 2) = 'tke_q'
        gfs_dyn_tracer%vname(ntke, 3) = 'tke_m'
        gfs_dyn_tracer%vname(ntke, 4) = 'tke_q6'
        gfs_dyn_tracer%vname(ntke, 5) = 'tke_m6'
      endif
      if(gfs_dyn_tracer%ntrac_met == 5) then
        gfs_dyn_tracer%vname(4, 1) = 'o' 
        gfs_dyn_tracer%vname(4, 2) = 'o_q' 
        gfs_dyn_tracer%vname(4, 3) = 'o_m' 
        gfs_dyn_tracer%vname(4, 4) = 'o_q6' 
        gfs_dyn_tracer%vname(4, 5) = 'o_m6' 
        gfs_dyn_tracer%vname(5, 1) = 'o2' 
        gfs_dyn_tracer%vname(5, 2) = 'o2_q' 
        gfs_dyn_tracer%vname(5, 3) = 'o2_m' 
        gfs_dyn_tracer%vname(5, 4) = 'o2_q6' 
        gfs_dyn_tracer%vname(5, 5) = 'o2_m6' 
      endif

      gfs_dyn_tracer%cpi(0:gfs_dyn_tracer%ntrac_met) =
     &               cpi(0:gfs_dyn_tracer%ntrac_met)
      gfs_dyn_tracer%ri(0:gfs_dyn_tracer%ntrac_met) =
     &               ri(0:gfs_dyn_tracer%ntrac_met)

! 
!--- fill in chem tracers
      if ( gfs_dyn_tracer%ntrac_chem > 0 ) then      
        do i = 1,gfs_dyn_tracer%ntrac_chem
          j = i + gfs_dyn_tracer%ntrac_met
          rgname = trim(gfs_dyn_tracer%chem_name(i))
          if(me==0)print *, 'LU_TRC_dyn: vname=',j,rgname
          gfs_dyn_tracer%vname(j, 1) = rgname
          gfs_dyn_tracer%vname(j, 2) = rgname//'_q'
          gfs_dyn_tracer%vname(j, 3) = rgname//'_m'
          gfs_dyn_tracer%vname(j, 4) = rgname//'_q6'
          gfs_dyn_tracer%vname(j, 5) = rgname//'_m6'
          gfs_dyn_tracer%cpi(j) = 0.
          gfs_dyn_tracer%ri (j) = 0.
        enddo
      endif

      endif     !!

      return

999   print *,'LU_TRC: error in allocate gfs_dyn_tracer :',status,me

      end subroutine tracer_config_init

! ========================================================================= 

      end module gfs_dyn_tracer_config
