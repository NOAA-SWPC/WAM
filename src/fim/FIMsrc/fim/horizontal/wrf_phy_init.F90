module module_wrf_phy_init

contains

  subroutine wrf_phy_init
!*********************************************************************
!       Loads the initial variables and constants for the WRF physics 
!       component.  
!*********************************************************************

    use module_control    ,only: nvl,nip
    use module_wrf_control,only: ids,ide,jds,jde,kds,kde, &
                                 ims,ime,jms,jme,kms,kme, &
                                 its,ite,jts,jte,kts,kte
    use module_wrfphysvars
    use module_wrfphys_alloc,only: wrfphys_alloc
    use module_wrf_share    ,only: wrf_set_array_bounds
    use module_wrf_variables,only: phys3dwrf,phys2dwrf,exch,pb2d

!TODO:  clean up duplication in these modules and add "only"
! for chemistry and WRF physics namelists:
    use module_chem_namelist_defaults
! contains config_flags
    use module_initial_chem_namelists  !, only: cu_physics, mp_physics
! TBH:  Ignore these so PPP doesn't have to translate them
!SMS$ignore begin
    use module_cu_gd        ,only: gdinit
    use module_species_decs
    use module_set_wrfphys
    use units, only: getunit, returnunit
!SMS$ignore end

    implicit none

! Local variables
    character(64) :: filename 
    real*8        :: t0,t1=0.0d0
    integer :: unitno             ! Unit number for I/O

!SMS$insert integer :: mype

!SMS$insert call nnt_me(mype)

    call StartTimer(t0)

    call wrf_set_array_bounds(nvl,nip,                 &
                              ids,ide,jds,jde,kds,kde, &
                              ims,ime,jms,jme,kms,kme, &
                              its,ite,jts,jte,kts,kte)
    print *,'start reading WRF physics namelists'
    filename='./FIMnamelist'
    print *,'set wrfphys namelists'
    call set_wrfphys_namelist_defaults
    unitno = getunit ()
    if (unitno < 0) then
      print*,'wrf_phy_init: getunit failed. Stopping'
      stop
    end if

    open (unitno, file=filename, form='formatted', action='read', err=70)
    read (unitno, wrfphysics, err=90)
    print *,'read wrfphys, cu_phys, mp_physics = ',cu_physics,mp_physics
    close(unitno)

    open (unitno, file=filename, form='formatted', action='read', err=70)
    read (unitno, chemwrf, err=90)
    print *,'read chem, chem_opt = ',chem_opt
    close(unitno)
    call returnunit (unitno)

    if (chem_opt == 0 .and. mp_physics == 0 .and. cu_physics == 0) return
    config_flags%mp_physics = mp_physics
    config_flags%cu_physics = cu_physics
    call set_wrfphys (mp_physics)
!
! we need the wrfphys variables for wrfchem
!
    if ((.not.mp_physics==0).or.(.not.cu_physics==0).or.(.not.chem_opt==0)) then
      write(0,*)'allocatewrfphys variables'
      call wrfphys_alloc
    end if   ! if ((.not.mp_physics==0).or.(.not.cu_physics==0)) then
!TODO:  move to wrfphys_alloc() ??  
!TODO:  avoid allocation when these variables are not used
    allocate( exch(nvl,nip)) ! 
    exch = 0.
    allocate( pb2d(nip)) ! 
    pb2d = 0.
    allocate( phys3dwrf(nvl,nip,11)) ! WRF Physics diagnostic variable to store tendencies for microphys and cu
                                     ! 1 = rqvcu
                                     ! 2 = rqvbl
                                     ! 3 = rqvf
                                     ! 4 = rthcu
                                     ! 5 = rthbl
                                     ! 6 = rthra
                                     ! 7 = rthf
                                     ! 8 = rqccu
                                     ! 9 = rqrcu
                                     ! 10 = rqscu
                                     ! 11 = rqicu
    allocate( phys2dwrf(nip,8))      ! Physics diagnostic variable

    phys3dwrf = 0.
    phys2dwrf = 0.

    if (.not.cu_physics==0) then
!TODO:  this really only applies to GD cumulus scheme -- improve "if" 
!TOD:   logic and generalize call
      CALL gdinit(RTHCUTEN,RQVCUTEN,RQCCUTEN,RQICUTEN,        &
                  MASS_FLUX,1004.6855,.false.,                &
                  0,0,0,                                      &
                  RTHFTEN, RQVFTEN,                           &
                  APR_GR,APR_W,APR_MC,APR_ST,APR_AS,          &
                  APR_CAPMA,APR_CAPME,APR_CAPMI,              &
                  .false.,                                    &
                  ids, ide, jds, jde, kds, kde,               &
                  ims, ime, jms, jme, kms, kme,               &
                  its, ite, jts, jte, kts, kte               )
    endif

    call IncrementTimer(t0,t1)
    print"(' WRF PHYSICS INIT time:',F10.0)",t1

    return

70  write(6,*)'wrf_phy_init: error opening unit=', unitno, '. Stopping'
    call flush(6)
    stop

90  write(6,*)'wrf_phy_init: error reading from unit=', unitno, '. Stopping'
    call flush(6)
    stop

  end subroutine wrf_phy_init
end module module_wrf_phy_init
