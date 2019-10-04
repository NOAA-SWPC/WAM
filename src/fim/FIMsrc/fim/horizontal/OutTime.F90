
!*********************************************************************
        module module_outtime_main
!       This module stores and prints elapsed wall-clock times for 
!       various parts of FIM.  
!       J. Middlecoff           April, 2008
!*********************************************************************

implicit none

real*8, save :: MainLoopTime=0.0d0

contains

subroutine OutTime(maxmin_times)

! iff .TRUE., print only max and min times
! otherwise, print times for many tasks
logical, optional, intent(in) :: maxmin_times

! local variables
logical :: maxmin_times_lcl
real :: MainLoopTimeMIN
real :: MainLoopTimeMAX

integer :: mype=0
!SMS$insert call nnt_me(mype)

maxmin_times_lcl = .false.
if (present(maxmin_times)) maxmin_times_lcl = maxmin_times

MainLoopTimeMIN = MainLoopTime
!SMS$reduce(MainLoopTimeMIN,MIN)
MainLoopTimeMAX = MainLoopTime
!SMS$reduce(MainLoopTimeMAX,MAX)

if (maxmin_times_lcl) then
  print"(' Main loop     ',2f15.3)",MainLoopTimeMIN,MainLoopTimeMAX
else
!SMS$ignore begin
  print"(' Main loop     ',f15.3,i10)",MainLoopTime,mype
!SMS$ignore end
endif

return
end subroutine OutTime 

end module module_outtime_main

!*********************************************************************
        module module_outtime_dyn
!       This module stores and prints elapsed wall-clock times for 
!       various parts of FIM dynamics.  
!       J. Middlecoff           April, 2008
!*********************************************************************

implicit none

real*8, save :: tdyn=0.0d0
real*8, save :: tout=0.0d0
real*8, save :: toutputBa=0.0d0
real*8, save :: thystat=0.0d0
real*8, save :: tedgvar=0.0d0
real*8, save :: tmomtum=0.0d0
real*8, save :: thybgen=0.0d0
real*8, save :: tprofout=0.0d0
real*8, save :: tabstart=0.0d0
real*8, save :: tcnuity=0.0d0
real*8, save :: ttrcadv=0.0d0
real*8, save :: ttransp=0.0d0
real*8, save :: tedgvarEx=0.0d0
real*8, save :: tedgvarBa=0.0d0
real*8, save :: tcnuityEx=0.0d0
real*8, save :: tcnuityBa=0.0d0
real*8, save :: ttrcadvEx=0.0d0
real*8, save :: ttrcadvBa=0.0d0
real*8, save :: ttranspEx=0.0d0
real*8, save :: ttranspBa=0.0d0
real*8, save :: tread_restart=0.0d0
real*8, save :: twrite_restart=0.0d0

contains

subroutine OutTime(maxmin_times,TimingBarriers,print_header)

! iff .TRUE., print only max and min times
! otherwise, print times for many tasks
logical, optional, intent(in) :: maxmin_times

! iff .TRUE., print barrier times
logical, optional, intent(in) :: TimingBarriers

! iff .TRUE., print header
logical, optional, intent(in) :: print_header

! local variables
logical :: maxmin_times_lcl
logical :: TimingBarriers_lcl
logical :: print_header_lcl

real :: tedgvarE,tcnuityE,ttrcadvE,ExchangeTime
real :: tedgvarB,tcnuityB,ttrcadvB,ExchBarrierTime
real :: tedgvarBMIN,tcnuityBMIN,ttrcadvBMIN,ExchBarrierTimeMIN
real :: tedgvarBMAX,tcnuityBMAX,ttrcadvBMAX,ExchBarrierTimeMAX

real :: tdynMIN,toutMIN,thybgenMIN,thystatMIN,tprofoutMIN
real :: tabstartMIN,tedgvarMIN,tcnuityMIN,ttrcadvMIN,tmomtumMIN
real :: tedgvarEMIN,ttrcadvEMIN,ExchangeTimeMIN
real :: tcnuityEMIN,toutputBMIN,tread_restartMIN,twrite_restartMIN

real :: tdynMAX,toutMAX,thybgenMAX,thystatMAX,tprofoutMAX
real :: tabstartMAX,tedgvarMAX,tcnuityMAX,ttrcadvMAX,tmomtumMAX
real :: tedgvarEMAX,ttrcadvEMAX,ExchangeTimeMAX
real :: tcnuityEMAX,toutputBMAX,tread_restartMAX,twrite_restartMAX

integer :: mype=0
!SMS$insert call nnt_me(mype)

maxmin_times_lcl = .false.
if (present(maxmin_times)) maxmin_times_lcl = maxmin_times
TimingBarriers_lcl = .false.
if (present(TimingBarriers)) TimingBarriers_lcl = TimingBarriers
print_header_lcl = .false.
if (present(print_header)) print_header_lcl = print_header

tedgvarE        = tedgvarEx
tcnuityE        = tcnuityEx
ttrcadvE        = ttrcadvEx
ExchangeTime    = tedgvarE+tcnuityE+ttrcadvE

tedgvarB        = tedgvarBa
tcnuityB        = tcnuityBa
ttrcadvB        = ttrcadvBa
ExchBarrierTime = tedgvarB+tcnuityB+ttrcadvB

tdynMIN         = tdyn
toutMIN         = tout
toutputBMIN     = toutputBa
thybgenMIN      = thybgen
thystatMIN      = thystat
tprofoutMIN     = tprofout
tabstartMIN     = tabstart
tedgvarMIN      = tedgvar
tcnuityMIN      = tcnuity
ttrcadvMIN      = ttrcadv
tmomtumMIN      = tmomtum
tedgvarEMIN     = tedgvarE
tcnuityEMIN     = tcnuityE
ttrcadvEMIN     = ttrcadvE
ExchangeTimeMIN = ExchangeTime
tedgvarBMIN     = tedgvarB
tcnuityBMIN     = tcnuityB
ttrcadvBMIN     = ttrcadvB
ExchBarrierTimeMIN = ExchBarrierTime
tread_restartMIN = tread_restart
twrite_restartMIN = twrite_restart

!SMS$reduce(tdynMIN,toutMIN,thybgenMIN,thystatMIN,toutputBMIN,
!SMS$>      tprofoutMIN,tabstartMIN,tedgvarMIN,tcnuityMIN,ttrcadvMIN,
!SMS$>      tmomtumMIN,tedgvarEMIN,tcnuityEMIN,
!SMS$>      ttrcadvEMIN,ExchangeTimeMIN,
!SMS$>      tedgvarBMIN,tcnuityBMIN,ttrcadvBMIN,ExchBarrierTimeMIN,
!SMS$>      tread_restartMIN,twrite_restartMIN,
!SMS$>      MIN)

tdynMAX         = tdyn
toutMAX         = tout
toutputBMAX     = toutputBa
thybgenMAX      = thybgen
thystatMAX      = thystat
tprofoutMAX     = tprofout
tabstartMAX     = tabstart
tedgvarMAX      = tedgvar
tcnuityMAX      = tcnuity
ttrcadvMAX      = ttrcadv
tmomtumMAX      = tmomtum
tedgvarEMAX     = tedgvarE
tcnuityEMAX     = tcnuityE
ttrcadvEMAX     = ttrcadvE
ExchangeTimeMAX = ExchangeTime
tedgvarBMAX     = tedgvarB
tcnuityBMAX     = tcnuityB
ttrcadvBMAX     = ttrcadvB
ExchBarrierTimeMAX = ExchBarrierTime
tread_restartMAX = tread_restart
twrite_restartMAX = twrite_restart

!SMS$reduce(tdynMAX,toutMAX,thybgenMAX,thystatMAX,toutputBMAX,
!SMS$>      tprofoutMAX,tabstartMAX,tedgvarMAX,tcnuityMAX,ttrcadvMAX,
!SMS$>      tmomtumMAX,tedgvarEMAX,tcnuityEMAX,
!SMS$>      ttrcadvEMAX,ExchangeTimeMAX,
!SMS$>      tedgvarBMAX,tcnuityBMAX,ttrcadvBMAX,ExchBarrierTimeMAX,
!SMS$>      tread_restartMAX,twrite_restartMAX,
!SMS$>      MAX)

if (maxmin_times_lcl) then

if (print_header_lcl) then
  if (TimingBarriers_lcl) then
    print"('                        MODULE TIME (sec)              EXCHANGE TIME (sec)            EXCH BARRIER TIME (sec)')"
    print"(' Module                  MIN            MAX             MIN            MAX             MIN            MAX')"
  else
    print"('                        MODULE TIME (sec)              EXCHANGE TIME (sec)')"
    print"(' Module                  MIN            MAX             MIN            MAX')"
  endif
endif
  print"(' Dynamics      ',2f15.3)",tdynMIN,tdynMAX
  print"(' Output        ',2f15.3)",toutMIN,toutMAX
  if (TimingBarriers_lcl) then
    print"(' Output Wait   ',2f15.3)",toutputBMIN,toutputBMAX
  endif
  print"(' hybgen        ',2f15.3)",thybgenMIN,thybgenMAX
  print"(' hystat        ',2f15.3)",thystatMIN,thystatMAX
  print"(' profout       ',2f15.3)",tprofoutMIN,tprofoutMAX
  print"(' abstart       ',2f15.3)",tabstartMIN,tabstartMAX
  if (TimingBarriers_lcl) then
    print"(' edgvar        ',6f15.3)",tedgvarMIN,tedgvarMAX,tedgvarEMIN,tedgvarEMAX,tedgvarBMIN,tedgvarBMAX
    print"(' cnuity        ',6f15.3)",tcnuityMIN,tcnuityMAX,tcnuityEMIN,tcnuityEMAX,tcnuityBMIN,tcnuityBMAX
    print"(' trcadv        ',6f15.3)",ttrcadvMIN,ttrcadvMAX,ttrcadvEMIN,ttrcadvEMAX,ttrcadvBMIN,ttrcadvBMAX
  else
    print"(' edgvar        ',4f15.3)",tedgvarMIN,tedgvarMAX,tedgvarEMIN,tedgvarEMAX
    print"(' cnuity        ',4f15.3)",tcnuityMIN,tcnuityMAX,tcnuityEMIN,tcnuityEMAX
    print"(' trcadv        ',4f15.3)",ttrcadvMIN,ttrcadvMAX,ttrcadvEMIN,ttrcadvEMAX
  endif
  print"(' momtum        ',2f15.3)",tmomtumMIN,tmomtumMAX
  print"(' Exchange      ',2f15.3)",ExchangeTimeMIN,ExchangeTimeMAX
  if (TimingBarriers_lcl) then
    print"(' ExchBarrier   ',2f15.3)",ExchBarrierTimeMIN,ExchBarrierTimeMAX
  endif
  print"(' read_restart  ',2f15.3)",tread_restartMIN,tread_restartMAX
  print"(' write_restart ',2f15.3)",twrite_restartMIN,twrite_restartMAX

else

!SMS$ignore begin
  print"(' Dynamics      ',f15.3,i10)",tdyn,mype
  print"(' Output        ',f15.3,i10)",tout,mype
  if (TimingBarriers_lcl) then
    print"(' Output_Wait   ',f15.3,i10)",toutputBa,mype
  endif
  print"(' hybgen        ',f15.3,i10)",thybgen,mype
  print"(' hystat        ',f15.3,i10)",thystat,mype
  print"(' profout       ',f15.3,i10)",tprofout,mype
  print"(' abstart       ',f15.3,i10)",tabstart,mype
  print"(' edgvar        ',f15.3,i10)",tedgvar,mype
  print"(' edgvarE       ',f15.3,i10)",tedgvarE,mype
  if (TimingBarriers_lcl) then
    print"(' edgvarB       ',f15.3,i10)",tedgvarB,mype
  endif
  print"(' cnuity        ',f15.3,i10)",tcnuity,mype
  print"(' cnuityE       ',f15.3,i10)",tcnuityE,mype
  if (TimingBarriers_lcl) then
    print"(' cnuityB       ',f15.3,i10)",tcnuityB,mype
  endif
  print"(' trcadv        ',f15.3,i10)",ttrcadv,mype
  print"(' trcadvE       ',f15.3,i10)",ttrcadvE,mype
  if (TimingBarriers_lcl) then
    print"(' trcadvB       ',f15.3,i10)",ttrcadvB,mype
  endif
  print"(' momtum        ',f15.3,i10)",tmomtum,mype
  print"(' Exchange      ',f15.3,i10)",ExchangeTime,mype
  if (TimingBarriers_lcl) then
    print"(' ExchBarrier   ',f15.3,i10)",ExchBarrierTime,mype
  endif
  print"(' read_restart  ',2f15.3)",tread_restart,mype
  print"(' write_restart ',2f15.3)",twrite_restart,mype
!SMS$ignore end

endif

return
end subroutine OutTime 

end module module_outtime_dyn



!*********************************************************************
        module module_outtime_phy
!       This module stores and prints elapsed wall-clock times for 
!       various parts of FIM physics.  
!       J. Middlecoff           April, 2008
!*********************************************************************

implicit none

real*8, save :: tphy=0.0d0

contains

subroutine OutTime(maxmin_times)

! iff .TRUE., print only max and min times
! otherwise, print times for many tasks
logical, optional, intent(in) :: maxmin_times

! local variables
logical :: maxmin_times_lcl
real :: tphyMIN
real :: tphyMAX

integer :: mype=0
!SMS$insert call nnt_me(mype)

maxmin_times_lcl = .false.
if (present(maxmin_times)) maxmin_times_lcl = maxmin_times

tphyMIN = tphy
!SMS$reduce(tphyMIN,MIN)
tphyMAX = tphy
!SMS$reduce(tphyMAX,MAX)

if (maxmin_times_lcl) then
  print"(' Physics       ',2f15.3)",tphyMIN,tphyMAX
else
!SMS$ignore begin
  print"(' Physics       ',f15.3,i10)",tphy,mype
!SMS$ignore end
endif

return
end subroutine OutTime 

end module module_outtime_phy




!*********************************************************************
        module module_outtime_chem
!       This module stores and prints elapsed wall-clock times for 
!       various parts of FIM chemistry.  
!       J. Middlecoff           April, 2008
!*********************************************************************

implicit none

real*8, save :: tchem=0.0d0

contains

subroutine OutTime(maxmin_times)

! iff .TRUE., print only max and min times
! otherwise, print times for many tasks
logical, optional, intent(in) :: maxmin_times

! local variables
logical :: maxmin_times_lcl
real :: tchemMIN
real :: tchemMAX

!SMS$insert integer :: mype
!SMS$insert call nnt_me(mype)

maxmin_times_lcl = .false.
if (present(maxmin_times)) maxmin_times_lcl = maxmin_times

tchemMIN = tchem
!SMS$reduce(tchemMIN,MIN)
tchemMAX = tchem
!SMS$reduce(tchemMAX,MAX)

if (maxmin_times_lcl) then
!SMS$insert print"(' Chemistry     ',2f15.3)",tchemMIN,tchemMAX
else
!SMS$ignore begin
!SMS$insert print"(' Chemistry     ',f15.3,i10)",tchem,mype
!SMS$ignore end
endif

return
end subroutine OutTime 

end module module_outtime_chem


!TODO:  DRY this out!  *Way* too much overlap with module_outtime_phy 
!TODO:  and others!  
!*********************************************************************
        module module_outtime_wrf_phy
!       This module stores and prints elapsed wall-clock times for 
!       various parts of WRF physics.  
!*********************************************************************

implicit none

real*8, save :: tphy=0.0d0

contains

subroutine OutTime(maxmin_times)

! iff .TRUE., print only max and min times
! otherwise, print times for many tasks
logical, optional, intent(in) :: maxmin_times

! local variables
logical :: maxmin_times_lcl
real :: tphyMIN
real :: tphyMAX

integer :: mype=0
!SMS$insert call nnt_me(mype)

maxmin_times_lcl = .false.
if (present(maxmin_times)) maxmin_times_lcl = maxmin_times

tphyMIN = tphy
!SMS$reduce(tphyMIN,MIN)
tphyMAX = tphy
!SMS$reduce(tphyMAX,MAX)

if (maxmin_times_lcl) then
  print"(' WRF Physics ',2f15.3)",tphyMIN,tphyMAX
else
!SMS$ignore begin
  print"(' WRF Physics ',f15.3,i10)",tphy,mype
!SMS$ignore end
endif

return
end subroutine OutTime 

end module module_outtime_wrf_phy




!*********************************************************************
        module module_outtime_cpl
!       This module stores and prints elapsed wall-clock times for 
!       various parts of FIM coupler.  
!       T. Henderson            March, 2009
!*********************************************************************

implicit none

real*8, save :: tcpl=0.0d0

contains

subroutine OutTime(maxmin_times,print_header)

! iff .TRUE., print only max and min times
! otherwise, print times for many tasks
logical, optional, intent(in) :: maxmin_times

! iff .TRUE., print header
logical, optional, intent(in) :: print_header

! local variables
logical :: maxmin_times_lcl
logical :: print_header_lcl
real :: tcplMIN
real :: tcplMAX

integer :: mype=0
!SMS$insert call nnt_me(mype)

maxmin_times_lcl = .false.
if (present(maxmin_times)) maxmin_times_lcl = maxmin_times

print_header_lcl = .false.
if (present(print_header)) print_header_lcl = print_header

tcplMIN = tcpl
!SMS$reduce(tcplMIN,MIN)
tcplMAX = tcpl
!SMS$reduce(tcplMAX,MAX)

if (maxmin_times_lcl) then
  if (print_header_lcl) then
    print"('                        MODULE TIME (sec)')"
    print"(' Module                  MIN            MAX')"
  endif
  print"(' Coupler       ',2f15.3)",tcplMIN,tcplMAX
else
!SMS$ignore begin
  print"(' Coupler       ',f15.3,i10)",tcpl,mype
!SMS$ignore end
endif

return

end subroutine OutTime

end module module_outtime_cpl
