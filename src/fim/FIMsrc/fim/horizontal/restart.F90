!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! module restart contains these public entities:
!   write_restart: writes a restart file
!   read_restart: reads a restart file
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module restart
  use module_control,     only: itsstart
  use units,              only: getunit, returnunit
  use module_outtime_dyn, only: twrite_restart, tread_restart
  
  implicit none
  
  private
  public :: write_restart, read_restart
  
contains

! write_restart: open the output restart file, then call write_restart_dyn and write_restart_phy
  subroutine write_restart (itsm1)
    integer, intent(in) :: itsm1     ! current time step to write to restart file

    integer           :: ivl         ! Index for vertical levels
    character(len=16) :: rfilename   ! Restart file name
    integer           :: unitno      ! Fortran unit number to write to
    integer           :: ret         ! function return
    real*8            :: t0          ! for timing
    character(len=8)  :: funcval     ! for Lahey

    character(len=8), external :: its2string ! Converts an integer to a string

#ifdef MANUALGPTL
#include <gptl.inc>
#endif

    write (6,*) 'Entered write_restart itsm1=', itsm1
!JR The following causes Lahey to fail! Looks like a compiler bug, so code around it:
!    write (rfilename, "('restart_',a8)") its2string (itsm1)
    funcval = its2string (itsm1)
    write (rfilename, "(a8,a8)") 'restart_', funcval

    unitno = getunit ()
    if (unitno < 0) then
      write(6,*)'write_restart: Bad return from getunit'
      stop
    end if

    call starttimer(t0)
!SMS$SERIAL BEGIN    
    open (unitno, file=rfilename, form='unformatted', action='write', err=70)
    write (unitno, err=90) itsm1
!SMS$SERIAL END

#ifdef MANUALGPTL
    ret = gptlprint_memusage ('before write_restart')
#endif
    call write_restart_dyn (unitno)
    call write_restart_phy (unitno)
#ifdef MANUALGPTL
    ret = gptlprint_memusage ('after write_restart')
#endif
    call incrementtimer (t0, twrite_restart)

    close (unitno)
    call returnunit (unitno)

#ifdef MANUALGPTL
    ret = gptlprint_memusage ('before system()')
#endif

!SMS$SERIAL (DEFAULT=IGNORE) BEGIN
    call system ('/bin/rm rpointer; ln -s ' // trim(rfilename) // ' rpointer')
    write(6,*)'restart: created link file rpointer -> ', trim(rfilename)
    call system ('/bin/ls -l rpointer ' // trim(rfilename))
!SMS$SERIAL END

#ifdef MANUALGPTL
    ret = gptlprint_memusage ('after system()')
#endif
    return

70  write(6,*)'write_restart: Error opening unit ', unitno, '. Stopping'
    stop

90  write(6,*)'write_restart: Error writing restart file. Stopping'
    stop
  end subroutine write_restart

! read_restart: open the input restart file, then call read_restart_dyn and read_restart_phy
  subroutine read_restart ()
    integer :: itsm1       ! time step read from restart file
    integer :: unitno      ! unit number for restart file
    integer :: ret         ! function return
    real*8  :: t0          ! for timing

#ifdef MANUALGPTL
#include <gptl.inc>
#endif

    unitno = getunit ()
    if (unitno < 0) then
      write(6,*)'read_restart: Bad return from getunit'
      stop
    end if
    
    call starttimer (t0)
!SMS$SERIAL BEGIN
    open (unitno, file='rpointer', form='unformatted', action='read', err=70)
    write(6,*)'read_restart: successfully opened restart file rpointer on unit ', unitno
    read (unitno, err=90) itsm1
!SMS$SERIAL END

    write (6,*) 'read_restart: starting up from itsm1=', itsm1

#ifdef MANUALGPTL
    ret = gptlprint_memusage ('before read_restart')
#endif
    call read_restart_dyn (unitno)
    call read_restart_phy (unitno)
#ifdef MANUALGPTL
    ret = gptlprint_memusage ('after read_restart')
#endif
    call incrementtimer (t0, tread_restart)

    close (unitno)
    call returnunit (unitno)
    write(6,*)'read_restart: returned unit ', unitno, ' to the list of available units'

    itsstart = itsm1 + 1     ! Set starting iteration counter for the restart run

    return

70  write(6,*)'read_restart: Error opening rpointer. Stopping'
    stop

90  write(6,*)'read_restart: Error reading restart file. Stopping'
    stop
  end subroutine read_restart
end module restart
