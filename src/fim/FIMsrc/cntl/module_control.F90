module module_control
!********************************************************************
!	This module specifies control variables for the global fim 
!  	A. E. MacDonald		October 11, 2004
!  	J. LEE         		September,  2005
!********************************************************************

implicit none
save

!  MODEL GRID LEVEL: glvl

!  Grid levels for the globe are:

!	glvl	Number of grid points	Linear Scale (km)
!	0		12                7,071
!	1		42                3,779 
!	2		162               1,174
!	3		642               891
!	4		2562              446
!	5		10,242            223
!	6		40,962            111
!	7		163,842           56
!	8		655,362           28
!	9		2,621,442         14
!	10		10,485,762        7
!	11		41,943,042        3.5
!	12		167,772,162       1.75

integer	, parameter :: nr=10            ! number of rhombi
integer	, parameter :: npp=6            ! number of proximity points(max)
integer	, parameter :: nd=2             ! number of directions (x,y)
integer , parameter :: nabl=3           ! number adams bashforth levels
integer , parameter :: VarNameLen=4     ! length of variable names for output routines
integer , parameter :: filename_len=80  ! max length of output filenames
!integer             :: ntr=4            ! # of tracers, 1=theta, 2=qv, 3=qw, 4=O3
integer             :: ntra=4           ! # of tracers advected on small dt, 1=theta, 2=qv, 3=qw, 4=O3
integer             :: ntrb=0           ! # of tracers advected on large dt (will include chemistry)
!integer , parameter :: nvarp=5          ! # of isobaric variables, 1=height, 2=temp, 3=RH,4=U,5=V
integer             :: nvarp=5          ! # of isobaric variables, 1=height, 2=temp, 3=RH,4=U,5=V
integer , parameter :: kbl=2            ! number of thin layers in surface boundary lyr
! integer , parameter :: nvlp = 40        ! number of isobaric levels - 1000-25 hPa
integer , parameter :: nvar2d = 6       ! number of extra 2d diagnostic variables for output
integer             :: stencl_frst=0    ! lowest layer for stencl diagnostics
integer             :: stencl_last=0    ! uppermost layer for stencl diagnostics
integer             :: stencl_step=0    ! layer increment in stencl diagnostics
integer             :: nvlp = 0         ! number of isobaric vertical levels - ex. 1000-25 hPa
integer, allocatable :: pres_hpa(:)     ! Holds isobaric vertical levels (hPa)
integer, allocatable :: pres_pa(:)      ! Holds isobaric vertical levels (Pa)

! namelist variables
!-------------------
integer       :: glvl                   ! the grid level defined in the Makefile
integer       :: nvl                    ! number of vertical native layers
integer       :: curve=1                ! 0: ij order; 1: Hilbert curve order; 2: ij block order, 3: Square Layout
integer       :: NumCacheBlocksPerPE=1  ! Number of cache blocks per processor. Used only for ij block order
logical       :: alt_topo=.false.       ! if true: use non-GFS surface height
character(2)  :: ArchvTimeUnit='hr'     ! ts:timestep; mi:minute; hr:hour; dy:day; mo:month
integer       :: TotalTime=5            ! total integration time in ArchvTimeUnit
integer       :: ArchvIntvl=6           ! archive interval in ArchvTimeUnit
integer       :: ArchvStep=1            ! archive interval in time steps
real          :: hrs_in_month=730.      ! length of month in hrs (=24*365/12)
integer       :: nts = -999             ! number of time steps to run: init to bad value

! OUTPUTnamelist items related to restart
integer       :: restart_freq = 999999  ! Interval (units=ArchvTimeUnit) to write restart file
integer       :: itsStart = 1           ! index of starting time step (is and always has been 1)
logical       :: readrestart = .false.  ! True means start by reading a restart file (default false)
! End of OUTPUTnamelist items related to restart

logical       :: UpdateSST      = .false. ! True means update SST field with interpolated monthly SST and Sea Ice Fraction (redef in FIMnamelist)
logical       :: EnKFAnl  = .false.     ! True means start by reading an EnKF analysis file (redef in FIMnamelist)
logical       :: EnKFIO   = .false.     ! True means special IO for EnKF analysis.
integer       :: PrintIpnDiag     = -1  ! Global ipn value at which to print diagnostics (-1 means no print)
integer       :: PrintDiagProgVars=  6  ! Hourly increment to print diagnostic prognosis variables (-1=>none, 0=>every step)
integer       :: PrintDiagNoise   =  1  ! Hourly increment to print diagnostic gravity wave noise  (-1=>none, 0=>every step)
logical       :: PrintDiags   = .false. ! True means print diagnostic messages (redef in FIMnamelist)
integer       :: PhysicsInterval  =  0  ! Interval in seconds to call physics, 0 => every time step
integer       :: RadiationInterval=3600 ! Interval in seconds to call radiation, 0 => every time step
integer       :: SSTInterval=-999       ! Interval in seconds to call sst, 0 => every time step (redef in FIMnamelist)
character(12) :: yyyymmddhhmm           ! Forecast initial time
logical       :: GravityWaveDrag=.true. ! True means calculate gravity wave drag
logical       :: digifilt=.false.       ! True means use digitial filter

logical       :: ras = .false.          ! false means call SAS
integer       :: num_p3d          =  4  ! 4 means call Zhao/Carr/Sundqvist Microphysics
character(80) :: EnKFFileName='enkf.anal' ! Name of the EnKF analysis file
character(80) :: FGFileName  ='fg.hybrid' ! Name of the EnKF analysis file
character(80) :: FGFileNameSig='fg.sigma' ! Name of the EnKF analysis file
character(80) :: sst_dat='HADISST_MONTHLY.1991-2009'  ! Name of SST file
character(80) :: ocean_bcs_ltln='ocean_bcs_ltln.360x180.dat'
logical       :: PrintMAXMINtimes=.true.! True means print MAX MIN routine times, false means print for each PE
logical       :: TimingBarriers=.false. ! True means insert timed barriers to measure task skew, will slow down model
logical       :: FixedGridOrder=.true.  ! True => always output in the same order (IJ), False => order determined by curve
!Control variables calculated in init.F90 from namelist variables
integer       :: CallPhysics            ! Timestep interval to call physics
integer       :: CallRadiation          ! Timestep interval to call radiation
integer       :: CallSST                ! Timestep interval to call radiation
integer       :: ipnDiagLocal     =  0  ! Local ipn value at which to print diagnostics
integer       :: ipnDiagPE        = -1  ! Processor on which ipnDiag resides.
integer       :: TestDiagProgVars       ! Calculated test value for printing diagnostic prognosis variables
integer       :: TestDiagNoise          ! Calculated test value for printing diagnostic gravity wave noise

integer       :: nip                    ! # of icosahedral points
integer       :: nvlp1                  ! # of vertical levels (= layers+1)
real          :: dt                     ! model time step (seconds)
integer       :: dtratio = 6            ! tracer B / tracer A timestep ratio
real          :: tfiltwin               ! length of digital filter window (secs)
integer       :: wts_type=3             ! type of digital filter window (1=Lanczos,2=Hamming,3=Dolph)
integer       :: numphr                 ! # of time steps/hr
integer	      :: nx	                ! rhombus x dimension
integer	      :: ny                     ! rhombus y dimension
real          :: rleigh_light = 0.      ! rayleigh damping time scale (days^-1) if top-layer wind < 100 m/s
real          :: rleigh_heavy = 0.2     ! rayleigh damping time scale (days^-1) if top-layer wind > 100 m/s
real          :: ptop= 10.		! pres at top of model domain (Pa)
real          :: thktop = 0.            ! min.thickness (Pa) of uppermost layer
integer       :: prev_date(4),next_date(4) ! time information for SST and ICE FRACTION slices held in memory
logical       :: have_next_sst=.false.  ! true when day = mid month, and sst's have been updated to next month


!TOPOnamelist
integer       :: toponpass=0         ! number of passes of shuman smoother-desmoother of input 5' wrf topo grid
real          :: toposmoothfact=1.25 ! radius of influence factor in grid lengths
character(120):: topodatfile='/no_such_file' ! path to topo dat file (set in FIMnamelist)
character(120):: topoglvldir='./'    ! dir of glvl.dat
character(len=80):: gfsltln_file ='NO_SUCH_FILE'
character(len=80):: mtnvar_file  ='NO_SUCH_FILE'
character(len=80):: aerosol_file ='NO_SUCH_FILE'
character(len=80):: co2_2008_file='NO_SUCH_FILE'
character(len=80):: co2_glb_file ='NO_SUCH_FILE'
character(len=80) :: isobaric_levels_file='NO_SUCH_FILE'

logical       :: pure_sig = .false.     ! if true, use pure sigma coord.
real          :: intfc_smooth = 50.     ! diffusivity (m/s) for intfc smoothing
real          :: slak = 0.5             ! intfc movement retardation coeff.
real          :: veldff_bkgnd = 0.      ! diffusion velocity (=diffusion/mesh size)
real          :: veldff_boost = 0.      ! veldff at model top (linear ramp-up over several layers)
real          :: dt_reducer_numerator   = 8. ! dt = dt * dt_reducer_numerator/dt_reducer_denominator
real          :: dt_reducer_denominator = 9.

integer       :: i = 0

contains

subroutine control(quiet_arg)
!SMS$ignore begin
use read_queue_namelist,only: ReturnGLVL,ReturnNVL,ReturnNIP,ReturnDT,GetWRFOn
use module_wrf_control, only: wrf_control   ! to re-compute ntra if needed
use units, only: getunit, returnunit
!SMS$ignore end
! arguments
logical, optional, intent(in) :: quiet_arg

! local variables
logical :: quiet
logical :: wrf_flag
integer :: unitno
!integer :: ioerr
! Define and read in the namelists
NAMELIST /PREPnamelist/  curve,NumCacheBLocksPerPE,alt_topo,gfsltln_file,mtnvar_file &
                        ,aerosol_file,co2_2008_file,co2_glb_file
NAMELIST /DIAGnamelist/  PrintIpnDiag,PrintDiagProgVars,PrintDiagNoise,PrintDiags
NAMELIST /MODELnamelist/ nts,	&
                         digifilt,wts_type,tfiltwin,   &
                         rleigh_light,rleigh_heavy,ptop,thktop,		&
                         pure_sig,intfc_smooth,slak,veldff_bkgnd,       &
                         veldff_boost,UpdateSST,dt_reducer_numerator,   &
                         dt_reducer_denominator,                        &
                         EnKFAnl,EnKFIO,sst_dat,ocean_bcs_ltln,dtratio
NAMELIST /TIMEnamelist/  yyyymmddhhmm
NAMELIST /OUTPUTnamelist/TotalTime,ArchvTimeUnit,ArchvIntvl,PrintMAXMINtimes, &
                         FixedGridOrder,TimingBarriers, &
                         restart_freq, readrestart, itsStart
                         
NAMELIST /ISOBARICnamelist/isobaric_levels_file
! namelist for calculating topo from wrf 5' data
namelist /TOPOnamelist/ toponpass,toposmoothfact,topodatfile,topoglvldir

quiet = .false.
if (present(quiet_arg)) then
  quiet = quiet_arg
end if

unitno = getunit ()
if (unitno < 0) then
  print*,'control: getunit failed for namelist files. Stopping'
  stop
end if

! Note:  REWIND required by IBM!  
! TODO:  Using open-read-close in place of REWIND until SMS is updated
!JR The following commented "open" call DOES NOT WORK! SMS puts a if(iam_root) around the 
!JR iostat=ioerr, meaning the test on slave nodes uses an uninitialized value! Instead use 
!JR antiquated "err=" feature.
!OPEN  (10, file="./FIMnamelist", status='old', action='read', iostat=ioerr)
!TODO: Fix the requirement to use antiquated f77 features

open  (unitno, file="./FIMnamelist", status='old', action='read', err=70)
write(6,*) 'control: successfully opened FIMnamelist'
read  (unitno, NML=PREPnamelist, err=90)
close (unitno)

!TODO: If we're reading a restart file, skip reading the other namelists and read from disk
open  (unitno, file="./FIMnamelist", status='old', action='read', err=70)
read  (unitno, NML=DIAGnamelist, err=90)
close (unitno)

open  (unitno, file="./FIMnamelist", status='old', action='read', err=70)
read  (unitno, NML=MODELnamelist, err=90)
close (unitno)

open  (unitno, file="./FIMnamelist", status='old', action='read', err=70)
read  (unitno, NML=TIMEnamelist, err=90)
close (unitno)

open  (unitno, file="./FIMnamelist", status='old', action='read', err=70)
read  (unitno, NML=OUTPUTnamelist, err=90)
close (unitno)

!JR Isnt it OK to enforce that toponamelist has to be there, but it can be empty?
open  (unitno, file="./FIMnamelist", status='old', action='read')
read  (unitno, NML=TOPOnamelist)  ! OK to use defaults if /TOPOnamelist/ not present
close (unitno)

open  (unitno, file="./FIMnamelist", status='old', action='read', err=70)
read  (unitno, NML=ISOBARICnamelist, err=75)
close (unitno)
! need to read file and do initialization here because this is needed by both fim and pop
open (unitno, file="./"//isobaric_levels_file, form="formatted", status="old", err=80)
print *, 'module_control: opened: ', isobaric_levels_file
read (unitno, *, err=81) nvlp

if (readrestart) then
  if (updatesst) then
    write(6,*)'control: restart does not yet work with updatesst=.true.'
    call flush(6)
    stop
  end if

  if (enkfanl) then
    write(6,*)'control: restart does not yet work with enkfanl=.true.'
    call flush(6)
    stop
  end if

  if (digifilt) then
    write(6,*)'control: restart does not yet work with digifilt=.true.'
    call flush(6)
    stop
  end if
end if

allocate (pres_hpa(nvlp))
allocate (pres_pa(nvlp))
read (unitno, *, err=82) pres_hpa
close (unitno)
call returnunit (unitno)

do i=1,nvlp
  pres_pa(i) = pres_hpa(i) * 100
end do

print *, '*** nvlp: ', nvlp
print '(a/(10i4))','pres (hPa):',pres_hpa
print '(a/(10i7))','pres (pa):',pres_pa

if (.not.quiet) then
  write(*, NML=PREPnamelist)
  write(*, NML=DIAGnamelist)
  write(*, NML=MODELnamelist)
  write(*, NML=TIMEnamelist)
  write(*, NML=OUTPUTnamelist)
  write(*, NML=TOPOnamelist)
  write(*, NML=ISOBARICnamelist)
end if

!SMS$serial begin
call ReturnGLVL(glvl)
call ReturnNVL ( nvl)
call ReturnNIP (nip)
call ReturnDT(dt)
call GetWRFOn(wrf_flag)
!SMS$serial end

if (glvl>=7) dt =  dt * dt_reducer_numerator / dt_reducer_denominator
nvlp1             = nvl+1                 ! # of vertical levels (= layers+1)
!dt               = 5760./2**(glvl-1)     ! model time step (seconds)
numphr            = nint(3600./dt)        ! # of time steps/hr
!dt                = 3600./float(numphr)
nx                = 2**glvl	          ! rhombus x dimension
ny                = 2**glvl	          ! rhombus y dimension

if (wrf_flag) then
  ! add WRF variables to ntr-dimensioned arrays, if needed
  call wrf_control(nvarp,ntrb)
end if
!ntra = ntr

if (curve == 0) then !The grid order is already IJ for curve=0.
  FixedGridOrder=.false.
end if

if (ArchvTimeUnit == 'ts') then
  nts = TotalTime
  ArchvStep = ArchvIntvl
else if (ArchvTimeUnit == 'mi') then
  nts = TotalTime*numphr/60
  ArchvStep = nint(60.*ArchvIntvl/dt)
  restart_freq = restart_freq*numphr/60
else if (ArchvTimeUnit == 'hr') then
  nts = TotalTime*numphr
  ArchvStep = nint(3600.*ArchvIntvl/dt)
  restart_freq = restart_freq*numphr
else if (ArchvTimeUnit == 'dy') then
  nts = TotalTime*24*numphr
  ArchvStep = nint(86400.*ArchvIntvl/dt)
  restart_freq = restart_freq*numphr*24
else if (ArchvTimeUnit == 'mo') then
  nts = TotalTime*hrs_in_month*numphr
  ArchvStep = nint(hrs_in_month*3600.*ArchvIntvl/dt)
  restart_freq = restart_freq*numphr*30*24
else
  write (*,'(a,a)') 'ERROR in module_control unrecognized output time unit: ',ArchvTimeUnit
  stop
end if

if      (PrintDiagProgVars == 0) then
  TestDiagProgVars = 1
else if (PrintDiagProgVars <  0) then
  TestDiagProgVars = 2*nts
else
  TestDiagProgVars = PrintDiagProgVars*numphr
end if

if      (PrintDiagNoise == 0) then
  TestDiagNoise = 1
else if (PrintDiagNoise <  0) then
  TestDiagNoise = 2*nts
else
  TestDiagNoise = PrintDiagNoise*numphr
end if

return

70 write(6,*)'control: error opening ./FIMnamelist'
call flush(6)
stop

75 write(6,*)'control: error reading ./ISOBARICnamelist'
call flush(6)
stop

80 write(6,*)'module_control: error reading isobaric_levels_file: COULD NOT OPEN: ', isobaric_levels_file, ' program aborted'
call flush(6)
stop

81 write(6,*)'module_control: error reading isobaric_levels_file - nvlp -  program aborted'
call flush(6)
stop

82 write(6,*)'module_control: error reading isobaric_levels_file - pres_ha - program aborted'
call flush(6)
stop

90 write(6,*)'control: error reading one of the namelists in ./FIMnamelist'
call flush(6)
stop

end subroutine control
end module module_control
