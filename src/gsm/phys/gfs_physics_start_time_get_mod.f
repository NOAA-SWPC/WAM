!#include "../../../ESMFVersionDefine.h"

      subroutine gfs_physics_start_time_get(yy, mm, dd, hh, mns, sec,  &
                                            kfhour, fhini,n1,cfile, rc)

! this subroutine gets and calculates the start time from reading the
! surface file information.

! !revision history:
!
!  march 2005      weiyu yang initial code.
!        2006      modified version 
!  november 2007   henry juang
!  Sep      2010   Jun Wang  change to nemsio file
!  Dec      2010   Jun Wang  change to nemsio library
!  Feb      2011   Sarah Lu  change to read nfhour from nemsio file header
!  sep      2011   Weiyu yang modified for using the ESMF 5.2.0r library.
!  Jan      2016   S Moorthi - ddd fhour in use date_def
!
!uses:
!
      use ESMF,         only: esmf_success
      use machine,      only: kind_io4, kind_evod
      use date_def,     only: idate,idate7,fhour
      use sfcio_module
      use nemsio_module

      implicit none

!
! arguments:
!-----------

      integer,                intent(out) :: yy, mm, dd, hh, mns, sec
      integer,                intent(out) :: n1, kfhour
      real(kind = kind_evod), intent(out) :: fhini
      integer,                intent(out) :: rc     ! return code
      character (len=*),intent(in)        :: cfile

      integer                        :: rc1 = esmf_success
!     real(kind = kind_evod)         :: fhour
      type(sfcio_head) head
      type(nemsio_gfile) nfile
      integer iret, khour

      n1 = 13
 
      call sfcio_sropen(n1,cfile,iret)
      call sfcio_srhead(n1,head,iret)
!     print *,'sfcio_srhead, iret=',iret
      if(iret == 0) then
        call sfcio_sclose(n1,iret)

        fhour  = head%fhour
        idate  = head%idate
        fhini  = fhour

        yy     = idate(4)
        mm     = idate(2)
        dd     = idate(3)
        hh     = idate(1)
        mns    = 0
        sec    = 0

        kfhour = nint(fhour) 

        rc = rc1
      else
        call sfcio_sclose(n1,iret)
!
        call nemsio_init()
        call nemsio_open(nfile,trim(cfile),'read',iret=iret)
!       print *,'in start time,after nemsio_open,iret=',iret

        call nemsio_getfilehead(nfile,idate=idate7,nfhour=kfhour,iret=iret)
        fhini = real(kfhour) 

!       write(0,*)'after nemsio,idate=',idate7,'nfhour=',kfhour,iret
        call nemsio_close(nfile)
        call nemsio_finalize()
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

      endif
!
!     write(0,*)' idate=',idate,' kfhour=',kfhour
      if (iret /= 0) call mpi_quit(5555)

      end subroutine gfs_physics_start_time_get
