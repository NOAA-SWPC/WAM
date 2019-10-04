!#include "../../../ESMFVersionDefine.h"

      subroutine gfs_dynamics_start_time_get(				&
                 yy, mm, dd, hh, mns, sec, kfhour, n1, n2,    		&
                 grib_inp, cfile, cfile2, nemsio_in,me,restart_run,rc)

! this subroutine gets and calculates the start time from reading the
! sigma file information.

! !revision history:
!
!  march 2005      weiyu yang initial code.
!  Feb   2010      jun wang   read data from nemsio file
!  Sep   2010      jun wang   remove gfsio option
!  Dec   2010      jun wang   change to nemsio library
!  Feb   2011      sarah lu   change to read nfhour
!  Sep   2012      Jun Wang   add sigio option
!  Aug   2013      Moorthi/Lu changed to read nemsio files (zeus porting)
!  Dec   2013      Jun Wang   for restart, read nemsio files
!
!uses:
!
      USE ESMF, ONLY: esmf_success

      use gfs_dyn_machine,      only : kind_io4, kind_evod
      use gfs_dyn_date_def,     only : idate,idate7
      use nemsio_module
      use sigio_module
      use sigio_r_module

      implicit none

!
! arguments:
!-----------

      integer,                intent(in)  :: grib_inp, me
      character (len=*),intent(in)        :: cfile, cfile2
      logical,intent(in)                  :: nemsio_in, restart_run
      integer,                intent(out) :: yy, mm, dd, hh, mns, sec
      integer,                intent(out) :: n1, n2
      integer,                intent(out) :: kfhour
      integer,                intent(out) :: rc     ! return code

      integer                             :: rc1 = esmf_success
      integer                             :: nfhour 
      type(nemsio_gfile) :: nfile
      type(sigio_head)   :: head
      integer iret, khour

      n1    = 11
      n2    = 12
 
      if (me == 0) print *,' dyn_start_time,grib_inp=',grib_inp,' cfile=',cfile,  &
    &          'nemsio_in=',nemsio_in
!        write(0,*)'in dyn_start_time'

      if(nemsio_in .or. restart_run) then
        call nemsio_init()
        call nemsio_open(nfile,trim(cfile),'read',iret=iret)
!        print *,'nemsio open instart iret=',iret
        call nemsio_getfilehead(nfile,idate=idate7,nfhour=kfhour,   & 
     &     iret=iret)                                            
        call nemsio_close(nfile)
        call nemsio_finalize()
!
        idate(1)=idate7(4)
        idate(2:3)=idate7(2:3)
        idate(4)=idate7(1)
        yy     = idate7(1)
        mm     = idate7(2)
        dd     = idate7(3)
        hh     = idate7(4)
        mns    = idate7(5)
        if(idate7(7)/=0) then
          sec    = idate7(6)*1./idate7(7)
        else
          sec    = 0
        endif
       else
         if (me == 0) print *,'in start_time,sigio input option,n1=',n1
         call sigio_rropen(n1,cfile,iret)
         if (me == 0) print *,'open cfile=',trim(cfile),' iret=',iret
         call sigio_rrhead(n1,head,iret)
         if (me == 0) print *,'rhead cfile=',trim(cfile),' iret=',iret
!
        idate = head%idate
        yy    = idate(4)
        mm    = idate(2)
        dd    = idate(3)
        hh    = idate(1)
        mns = 0.
        sec = 0.
        kfhour  = head%fhour
        idate7  = 0
        idate7(1:4) = idate(1:4)
        idate7(1)   = yy
        idate7(4)   = hh
        idate7(7)   = 100
!
        call sigio_rclose(n1,iret)

       endif

!      print *,' fhour=',fhour,' idate=',idate7,' iret=',iret

      if (iret /= 0) call mpi_quit(5555)

      if (me == 0) print *,' idate=',idate,' kfhour=',kfhour

      rc = rc1

      end subroutine gfs_dynamics_start_time_get
