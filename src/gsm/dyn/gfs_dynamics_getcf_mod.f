!
! !module: gfs_dynamics_getcf_mod --- configure data module of the 
!                          gridded component of the gfs dynamics system.
!
! !description: this program uses the esmf configure class software
!               to get all parameters of the original namelists of
!               gfs dynamics and all trun-on/turn-off switch flags 
!               import state and the export state.
!
! !revision history:
!
!  november 2004 weiyu yang    initial code for esmf gfs model.
!  may      2005 weiyu yang    for the updated gfs version.
!  march    2006 s. moorthi    modified for the new gfs.
!  january  2007 h. juang      modified for the gfs dynamics.
!  May      2009 j. wang       modified for the gfs wrt grid comp
!  Oct 2009 Sarah Lu           tracer added; (q, oz, cld) removed
!  November 2009 weiyu yang    modified for the ensemble NEMS run.
!  Sep      2012 j. wang       add nemsio_in to specify input files
!
! !interface:
!
      module gfs_dynamics_getcf_mod

!
!!uses:
!

! the esmf internal state contents.
!----------------------------------
      use gfs_dynamics_internal_state_mod
      use gfs_dynamics_err_msg_mod

      implicit none

      contains

!------------------------------------------------------------------------
!bop
!
! !iroutine:
!
! !interface:
!
     subroutine gfs_dynamics_getcf(gc_gfs_dyn, int_state, rc)

      type(esmf_gridcomp),                intent(inout) :: gc_gfs_dyn 
      type(gfs_dynamics_internal_state),  intent(inout) :: int_state 
      type(esmf_vm)                                     :: vm           
      integer,                            intent(out)   :: rc  
!
! !description: load run parameters from the resource file into internal
!               state.

!
! !revision history:
!
!  november 3, 2004 weiyu yang, for ncep gfs model.
!  may         2005 weiyu yang  for the updated gfs version.
!  january     2007 hann-ming henry juang  for the  gfs dynamics version.
!
!eop
!------------------------------------------------------------------------

      type(esmf_config) :: cf
      integer           :: ndat, i, rc1, rcfinal
      character*20      :: clab, pelab
      integer           :: pe_member,j, i1, me, tasks, totpe

! initialize the error signal variables.
!---------------------------------------
      rc1     = esmf_success
      rcfinal = esmf_success

! get the esmf config from the gfs grid component.
!-------------------------------------------------
      call esmf_gridcompget(gc_gfs_dyn, config = cf, rc = rc1)

      call gfs_dynamics_err_msg(rc1,				&
               'grid component get configure',rcfinal)

! start to get the configure data.
!---------------------------------

! the original gfs dynamics has namelist file, "nam_gfs_dyn".
! besides to get the information from namelist file, the getcf
! process will get the "esmf_state_namelist" which contains all turn-on/turn-off
! switch flags for all possible fields in the esmf import and export states.
! the single esmf config file contains all information of the namelist files.
!-------------------------------------------------------------------------------

! for "nam_gfs_dyn".
!---------------

! "esmf_configgetattribute" is an esmf config utility routine which is used
! to get information from the config file.  
! the first argument is the gfs esmf config.
! the second argument is the destination variable 
! to which the value of the namelist 
! variable will be sent. the third one is the label 
! in the config file which is used
! to identify the required namelist variable.  
! the last one is the error signal variable.
! "esmf_configgetattribute" is a generic interface name 
! which can be used to get the 
! the different data types data from the config file, 
! such as real(4), real(8), 
! integer(4), integer(8), character(...), etc..
!------------------------------------------------------------------------------
      call esmf_vmgetglobal(vm, rc = rc1)
      call esmf_vmget(vm, localpet = me, pecount = tasks, rc = rc1)

      call esmf_configgetattribute(cf, 					&
                             int_state%nam_gfs_dyn%total_member,   	&
                             label = 'total_member:', rc    = rc1)

      IF(int_state%nam_gfs_dyn%total_member == 1) THEN
          int_state%ENS = .FALSE.
      ELSE
          int_state%ENS = .TRUE.
      END IF

      call gfs_dynamics_err_msg_var(rc1,                                &
                 'gfs dynamics getcf','total_member',rcfinal)

!      print *,' total_member=',int_state%nam_gfs_dyn%total_member
!      print *,' total tasks=',tasks

      call esmf_configgetattribute(cf, 					&
                             int_state%grib_inp,                	&
                             label = 'grib_input:', rc    = rc1)
      call gfs_dynamics_err_msg_var(rc1,                                &
                 'gfs dynamics getcf','grib_input',rcfinal)

      totpe = 0
      i1 = 0

!     print *,' before j loop up to ',1,int_state%nam_gfs_dyn%total_member

      do j=1,int_state%nam_gfs_dyn%total_member
        write(pelab,'("PE_MEMBER",i2.2,":")') j
        call esmf_configgetattribute(cf, 				&
                             pe_member, label = pelab, rc = rc1)
        call gfs_dynamics_err_msg_var(rc1,                              &
                 'gfs dynamics getcf',pelab,rcfinal)
        if (pe_member == 0) 						&
            pe_member = tasks / int_state%nam_gfs_dyn%total_member

!       print *,' pe_member=',pe_member
!       print *,' tasks=',tasks

        totpe = totpe + pe_member
        do i = 1, pe_member
          if (me == i1) then
             int_state%nam_gfs_dyn%member_id = j
          end if
          i1 = i1+1
        end do
      end do

      if (totpe /= tasks) then
       print *,' totpe=',totpe,' and tasks=',tasks, ' do not match'
       stop 9999
      endif

      call esmf_configgetattribute(cf, 					&
                              int_state%nam_gfs_dyn%nlunit,         	&
                              label = 'nlunit:',  rc = rc1)
      call gfs_dynamics_err_msg_var(rc1,				&
               'gfs dynamics getcf','nlunit',rcfinal)

!jws
      call esmf_configgetattribute(cf,                                  &
                              int_state%adiabatic,                      &
                              label = 'adiabatic:',  rc = rc1)
      call gfs_dynamics_err_msg_var(rc1,                                &
               'gfs dynamics getcf','adiabatic',rcfinal)
      call esmf_configgetattribute(cf,                                  &
                              int_state%restart_run,                    &
                              label = 'restart:',  rc = rc1)
      call esmf_configgetattribute(cf,                                  &
                              int_state%nemsio_in,                      &
                              label = 'nemsio_in:',  rc = rc1)
!jwe

      call esmf_configgetattribute(cf, 					&
                              int_state%nam_gfs_dyn%gfs_dyn_namelist, 	&
                              label = 'namelist:',  rc = rc1)
      call gfs_dynamics_err_msg_var(rc1,				&
                 'gfs dynamics getcf','namelist',rcfinal)

      call esmf_configgetattribute(cf, 					&
                              int_state%nam_gfs_dyn%deltim,         	&
                              label = 'deltim:', rc = rc1)
      int_state%deltim = int_state%nam_gfs_dyn%deltim

      call gfs_dynamics_err_msg_var(rc1,				&
               'gfs dynamics getcf','deltim',rcfinal)

      if (int_state%nam_gfs_dyn%total_member <= 1) then
        int_state%nam_gfs_dyn%grid_ini  = 'grid_ini'
        int_state%nam_gfs_dyn%grid_ini2 = 'grid_ini2'
        int_state%nam_gfs_dyn%sig_ini  = 'sig_ini'
        int_state%nam_gfs_dyn%sig_ini2 = 'sig_ini2'
      else
        write(int_state%nam_gfs_dyn%grid_ini,                          &
              '("grid_ini_",i2.2)') int_state%nam_gfs_dyn%member_id
        write(int_state%nam_gfs_dyn%grid_ini2,                         &
              '("grid_ini2_",i2.2)') int_state%nam_gfs_dyn%member_id
        write(int_state%nam_gfs_dyn%sig_ini, 				&
              '("sig_ini_",i2.2)') int_state%nam_gfs_dyn%member_id
        write(int_state%nam_gfs_dyn%sig_ini2,				&
              '("sig_ini2_",i2.2)') int_state%nam_gfs_dyn%member_id
      endif

!
!jws get output file name
!--------------------------
!
      call esmf_configgetattribute(cf,                                  &
                              int_state%num_file,                       &
                              label = 'num_file:',    rc = rc1)

      call gfs_dynamics_err_msg_var(rc1,                                &
               'gfs dynamics getcf','num_file',rcfinal)
!
      allocate(int_state%filename_base(int_state%num_file))
      call esmf_configFindLabel(CF,'filename_base:',rc=RC)
      Do I=1,int_state%num_file
        call esmf_configgetattribute(cf,                                &
                              int_state%filename_base(i),               &
                              rc = rc1)
      enddo
      call gfs_dynamics_err_msg_var(rc1,                                &
               'gfs dynamics getcf','filename_base',rcfinal)
!jwe

      call esmf_configgetattribute(cf, 					&
                              int_state%slg_flag, label = 'SLG_FLAG:',  &
                              rc = rc1)
      call gfs_dynamics_err_msg_var(rc1,				&
               'gfs dynamics getcf', 'slg_flag', rcfinal)
!
!
! for "esmf_state_namelist"
!--------------------------
!
      call esmf_configgetattribute(cf, 					&
                              int_state%esmf_sta_list%idate1_import,    &
                              label = 'idate1_import:',    rc = rc1)
      call gfs_dynamics_err_msg_var(rc1,				&
               'gfs dynamics getcf','idate1_import',rcfinal)
!
      call esmf_configgetattribute(cf, 					&
                              int_state%esmf_sta_list%z_import,         &
                              label = 'z_import:',    rc = rc1)
      call gfs_dynamics_err_msg_var(rc1,				&
               'gfs dynamics getcf','z_import',rcfinal)
!
      call esmf_configgetattribute(cf, 					&
                              int_state%esmf_sta_list%ps_import,	&
                              label = 'ps_import:',    rc = rc1)
      call gfs_dynamics_err_msg_var(rc1,				&
               'gfs dynamics getcf','ps_import',rcfinal)
!
      call esmf_configgetattribute(cf, 					&
                              int_state%esmf_sta_list%u_import,		&
                              label = 'u_import:',    rc = rc1)
      call gfs_dynamics_err_msg_var(rc1,				&
               'gfs dynamics getcf','u_import',rcfinal)
!
      call esmf_configgetattribute(cf, 					&
                              int_state%esmf_sta_list%v_import,		&
                              label = 'v_import:',    rc = rc1)
      call gfs_dynamics_err_msg_var(rc1,				&
               'gfs dynamics getcf','v_import',rcfinal)
!
      call esmf_configgetattribute(cf, 					&
                              int_state%esmf_sta_list%temp_import,	&
                              label = 'temp_import:',    rc = rc1)
      call gfs_dynamics_err_msg_var(rc1,				&
               'gfs dynamics getcf','temp_import',rcfinal)
!
      call esmf_configgetattribute(cf,                                  &
                              int_state%esmf_sta_list%tracer_import,    &
                              label = 'tracer_import:',    rc = rc1)
      call gfs_dynamics_err_msg_var(rc1,                                &
               'gfs dynamics getcf','tracer_import',rcfinal)
!
      call esmf_configgetattribute(cf, 					&
                              int_state%esmf_sta_list%p_import,	        &
                              label = 'p_import:',    rc = rc1)
      call gfs_dynamics_err_msg_var(rc1,				&
               'gfs dynamics getcf','p_import',rcfinal)
!
      call esmf_configgetattribute(cf, 					&
                              int_state%esmf_sta_list%dp_import,	&
                              label = 'dp_import:',    rc = rc1)
      call gfs_dynamics_err_msg_var(rc1,				&
               'gfs dynamics getcf','dp_import',rcfinal)
!
      call esmf_configgetattribute(cf, 					&
                              int_state%esmf_sta_list%dpdt_import,	&
                              label = 'dpdt_import:',    rc = rc1)
      call gfs_dynamics_err_msg_var(rc1,				&
               'gfs dynamics getcf','dpdt_import',rcfinal)

!----------------------------------
!
      call esmf_configgetattribute(cf, 					&
                              int_state%esmf_sta_list%idate1_export,	&
                              label = 'idate1_export:',    rc = rc1)
      call gfs_dynamics_err_msg_var(rc1,				&
               'gfs dynamics getcf','idate1_export',rcfinal)
!
      call esmf_configgetattribute(cf, 					&
                              int_state%esmf_sta_list%z_export, 	&
                              label = 'z_export:',    rc = rc1)
      call gfs_dynamics_err_msg_var(rc1,				&
               'gfs dynamics getcf','z_export',rcfinal)
!
      call esmf_configgetattribute(cf, 					&
                              int_state%esmf_sta_list%ps_export,	&
                              label = 'ps_export:',    rc = rc1)
      call gfs_dynamics_err_msg_var(rc1,				&
               'gfs dynamics getcf','ps_export',rcfinal)
!
      call esmf_configgetattribute(cf, 					&
                              int_state%esmf_sta_list%u_export,		&
                              label = 'u_export:',    rc = rc1)
      call gfs_dynamics_err_msg_var(rc1,				&
               'gfs dynamics getcf','u_export',rcfinal)
!
      call esmf_configgetattribute(cf, 					&
                              int_state%esmf_sta_list%v_export,		&
                              label = 'v_export:',    rc = rc1)
      call gfs_dynamics_err_msg_var(rc1,				&
               'gfs dynamics getcf','v_export',rcfinal)
!
      call esmf_configgetattribute(cf, 					&
                              int_state%esmf_sta_list%temp_export,	&
                              label = 'temp_export:',    rc = rc1)
      call gfs_dynamics_err_msg_var(rc1,				&
               'gfs dynamics getcf','temp_export',rcfinal)

      call esmf_configgetattribute(cf,                                  &
                              int_state%esmf_sta_list%tracer_export,    &
                              label = 'tracer_export:',    rc = rc1)
      call gfs_dynamics_err_msg_var(rc1,                                &
               'gfs dynamics getcf','tracer_export',rcfinal)

      call esmf_configgetattribute(cf, 					&
                              int_state%esmf_sta_list%p_export,	        &
                              label = 'p_export:',    rc = rc1)
      call gfs_dynamics_err_msg_var(rc1,				&
               'gfs dynamics getcf','p_export',rcfinal)

      call esmf_configgetattribute(cf, 					&
                              int_state%esmf_sta_list%dp_export,	&
                              label = 'dp_export:',    rc = rc1)
      call gfs_dynamics_err_msg_var(rc1,				&
               'gfs dynamics getcf','dp_export',rcfinal)

      call esmf_configgetattribute(cf, 					&
                              int_state%esmf_sta_list%dpdt_export,	&
                              label = 'dpdt_export:',    rc = rc1)
      call gfs_dynamics_err_msg_var(rc1,				&
               'gfs dynamics getcf','dpdt_export',rcfinal)

      call esmf_configgetattribute(cf, 					&
                              int_state%esmf_sta_list%shum_wts_export,	&
                              label = 'shum_wts_export:',    rc = rc1)
      call gfs_dynamics_err_msg_var(rc1,				&
               'gfs dynamics getcf','shum_wts_export',rcfinal)

      call esmf_configgetattribute(cf, 					&
                              int_state%esmf_sta_list%sppt_wts_export,	&
                              label = 'sppt_wts_export:',    rc = rc1)
      call gfs_dynamics_err_msg_var(rc1,				&
               'gfs dynamics getcf','sppt_wts_export',rcfinal)

      call esmf_configgetattribute(cf, 					&
                              int_state%esmf_sta_list%skeb_wts_export,	&
                              label = 'skeb_wts_export:',    rc = rc1)
      call gfs_dynamics_err_msg_var(rc1,				&
               'gfs dynamics getcf','skeb_wts_export',rcfinal)

      call esmf_configgetattribute(cf, 					&
                              int_state%esmf_sta_list%vc_wts_export,	&
                              label = 'vc_wts_export:',    rc = rc1)
      call gfs_dynamics_err_msg_var(rc1,				&
               'gfs dynamics getcf','vc_wts_export',rcfinal)

!----------------------------------

      call gfs_dynamics_err_msg_final(rcfinal,'gfs_dynamics_getcf',rc)

      end subroutine gfs_dynamics_getcf

      end module gfs_dynamics_getcf_mod
