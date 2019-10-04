   subroutine checkiauforcing(t,INIAUINTERVAL)
        
! check to see if model is in IAU window
   use namelist_physics_def,only : iaufhrs,iau_delthrs,iaufiles_anl
   use machine, only:kind_evod
   implicit none
   real(kind_evod), intent(in) :: t
   LOGICAL,intent(out) :: INIAUINTERVAL
   real(kind_evod) dt
   integer n,t1,t2,nfiles,nfilesall
  
   dt = iau_delthrs*3600.
   INIAUINTERVAL=.false.

   nfilesall = size(iaufiles_anl)
   nfiles = 0
   do n=1,nfilesall
      if (trim(iaufiles_anl(n)) .eq. '' .or. iaufhrs(n) .lt. 0) exit
      nfiles = nfiles + 1
   enddo
   ! set forcing to zero and return if outside iau window.
   if ( nfiles > 1) then  ! IAU forcing files bookend interval
      if (t <= iaufhrs(1)*3600. .or. t > iaufhrs(nfiles)*3600.) then
         return
      endif
   else  ! single file at middle of window
      t1=iaufhrs(1)*3600 - dt*0.5
      t2=iaufhrs(1)*3600 + dt*0.5
      if ( t <= t1 .or. t > t2 ) then
         return
      endif
   endif
   if (nfiles > 1) then
      if (t .eq. 3600.*iaufhrs(nfiles)) then
         return
      else if (t .eq. 3600.*iaufhrs(1)) then
         return
      endif
      do n=1,nfiles
         if (iaufhrs(n)*3600. > t) exit
      enddo
   endif
   INIAUINTERVAL=.true.
 end subroutine checkiauforcing
