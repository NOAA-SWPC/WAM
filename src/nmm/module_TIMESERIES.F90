!-----------------------------------------------------------------------
      module module_timeseries
!-----------------------------------------------------------------------

      use module_kinds
      use module_constants
      use module_my_domain_specs
      use module_vars

!-----------------------------------------------------------------------

      implicit none

!-----------------------------------------------------------------------

      private

      public :: timeseries_initialize, timeseries_run, timeseries_finalize

      integer, parameter :: max_point=20
      logical, save :: initialized=.false.
      integer, save :: npoints
      real, dimension(max_point), save :: points_lon, points_lat
      integer, dimension(max_point), save :: ipnt, jpnt
      real, dimension(max_point), save :: pnt_lon, pnt_lat

      real, parameter :: rtd=180.0/pi

      real, allocatable, dimension(:,:) :: zint   ! height of full levels
      real, allocatable, dimension(:,:) :: zmid   ! height of half levels

      integer, dimension(max_point), save :: tsunit

      integer :: var2d_number
      character(len=8), dimension(1000) :: var2d_name
      character(len=24), dimension(1000) :: var2d_units
      character(len=64), dimension(1000) :: var2d_description
      integer :: var3d_number
      character(len=8), dimension(1000) :: var3d_name
      character(len=24), dimension(1000) :: var3d_units
      character(len=64), dimension(1000) :: var3d_description
      character(len=4), dimension(1000) :: var3d_lvlind

      logical :: nml_exist

      integer, parameter :: max_fulllevel_vars = 10
      character(len=8), dimension(max_fulllevel_vars) :: fulllevel_vars &
           = reshape ( (/                                               &
                              'PINT    '                                &
                       /)                                               &
                      ,(/max_fulllevel_vars/)                           &
                      ,(/'********'/)                                   &
                     )
!-----------------------------------------------------------------------

      contains

!-----------------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!-----------------------------------------------------------------------

      subroutine timeseries_initialize(solver_state,ntsd,ierr)

      use module_solver_internal_state, only : solver_internal_state

      implicit none

      type(solver_internal_state),intent(in) :: solver_state
      integer, intent(in) :: ntsd
      integer, intent(out) :: ierr

      integer :: ios

      integer :: im, jm, lm
      real :: tlm0d,tph0d,dlmd,dphd,wbd,sbd
      real :: tlmd,tphd
      integer :: i,k,j, np, var, inrs,jnrs, nvar, n, l
      logical :: inside
      integer, dimension(8) :: modelstarttime
      character(len=3) :: trnum

      integer year, month, day, hour, minute, second, ten_thousandth
      integer :: ihr
      integer :: nf_hours,nf_minutes
      real    :: nf_seconds
      real    :: secfcst
      integer :: idatv(3),ihrv,idaywk

      character*128  :: filename

      namelist /ts_locations/ npoints, points_lon, points_lat


      ierr = 0

      if (initialized) return

      inquire ( file='ts_locations.nml', exist=nml_exist)

      if (.not.nml_exist) then
        if (mype == 0) then
          write(0,*) ' ts_locations.nml does not exist. will skip timeseries output'
        end if
        return
      end if

      open ( unit=80, file='ts_locations.nml', status='old', iostat = ios)
      if (ios /= 0) then
         write(0,*) 'error missing ts_locations.nml'
         ierr = -1
         return
      end if
      read  ( unit=80, nml=ts_locations,  iostat = ios)
      if (ios /= 0) then
         write(0,*) 'error reading ts_locations.nml'
         ierr = -1
         return
      end if
      close(unit=80,iostat=ios)
      if (ios /= 0) then
         write(0,*) 'error closing ts_locations.nml'
         ierr = -1
         return
      end if

      im = solver_state%im
      jm = solver_state%jm
      lm = solver_state%lm
      tlm0d = solver_state%tlm0d
      tph0d = solver_state%tph0d
      dlmd = solver_state%dlmd
      dphd = solver_state%dphd
      wbd = solver_state%wbd
      sbd = solver_state%sbd

!! calculate forecast time for this time step

      secfcst = solver_state%dt * ntsd
      nf_hours = int(secfcst/3600)
      ihr = nf_hours
      nf_minutes = int( mod(secfcst,3600.0)/60.0 )
      nf_seconds = (secfcst - nf_hours*3600.0) - nf_minutes*60.0
      if (nf_seconds< 0.000001) nf_seconds=0.0
      ten_thousandth = nint((secfcst-int(secfcst))*10000)

      call valid(solver_state%idat,solver_state%ihrst,nf_hours,idatv,ihrv,idaywk)

      year = idatv(3)
      month = idatv(2)
      day = idatv(1)
      hour = ihrv
      minute = nf_minutes
      second = nf_seconds

      modelstarttime(1) = year
      modelstarttime(2) = month
      modelstarttime(3) = day
      modelstarttime(4) = hour
      modelstarttime(5) = minute
      modelstarttime(6) = second
      modelstarttime(7) = 0
      modelstarttime(8) = ntsd

      allocate(zint(npoints,lm+1))
      allocate(zmid(npoints,lm))

      var2d_number = 0
      var3d_number = 0


!! loop over all solver state variables and find which ones are selected for timeseries output.
!! currently only 2d and 3d real variables are supported. this is due to limitations of used
!! timeseries binary file format

      do n=1,solver_state%num_vars
        if (solver_state%vars(n)%tseries) then
        select case(solver_state%vars(n)%tkr)
          case(tkr_r2d)
            var2d_number = var2d_number + 1
            var2d_name(var2d_number) = trim(solver_state%vars(n)%vbl_name)
            var2d_units(var2d_number) = ""
            var2d_description(var2d_number) = trim(solver_state%vars(n)%description)
          case(tkr_r3d)
            var3d_number = var3d_number + 1
            var3d_name(var3d_number) = trim(solver_state%vars(n)%vbl_name)
            var3d_units(var3d_number) = ""
            var3d_description(var3d_number) = trim(solver_state%vars(n)%description)
            var3d_lvlind(var3d_number) = 'H   '

            !! for 3d variables, level indicator is set to 'H' by default which means 'half' level
            !! variable, or variable located at the middle of the layer in vertical. if the variable name
            !! is included in the list of "full" level (or interface) variables then
            !! level indicator is reset to 'F'

            do l=1,max_fulllevel_vars
               if (trim(fulllevel_vars(l))==trim(var3d_name(var3d_number))) then
                  var3d_lvlind(var3d_number) = 'F   '
                  exit
               end if
            end do
          case default
            write(0,*)' unknown tkr = ', solver_state%vars(n)%tkr, trim(solver_state%vars(n)%vbl_name)
            ierr = -1
            return
        end select
        end if
      end do

!! loop over number of points specified in ts_locations namelist and calculate i,j
!! indexes of the nearest H point

      np_loop: do np=1,npoints

         call tll(points_lon(np),points_lat(np),tlmd,tphd,tph0d,tlm0d)
         call ijnrs(tlmd,tphd,dlmd,dphd,wbd,sbd,im,jm,inrs,jnrs,inside)
         if (inside) then
            ipnt(np)=inrs
            jpnt(np)=jnrs
         end if

!! if this point, with i,j indexes ipnt(np),jpnt(np) is inside the tile (its:ite,jts:jte)
!! located on this PE then proceed and create output file and write out file header

         inside_if: if (its <= ipnt(np) .and. ipnt(np) <= ite .and. jts <= jpnt(np) .and. jpnt(np) <= jte ) then

            i = ipnt(np)
            j = jpnt(np)
            zint(np,lm+1)=solver_state%fis(i,j)/g

            do k=lm,1,-1
                  zint(np,k)=zint(np,k+1)+solver_state%t(i,j,k)*(0.608*max(solver_state%q(i,j,k),epsq)+1.)*r_d        &
                            *(solver_state%pint(i,j,k+1)-solver_state%pint(i,j,k))                                     &
                            /((solver_state%sgml2(k)*solver_state%pd(i,j)+solver_state%psgml1(k))*g)
                  zmid(np,k)=0.5*(zint(np,k+1)+zint(np,k))
            end do

!! open the timeseries output file
            tsunit(np) = 90 + np
            write(filename,fmt='(a,i2.2,a,i2.2,a)') "ts_p",np,"_d",my_domain_id,".bin"
            write(0,*) ' open file:', np,tsunit(np),filename
            open(unit=tsunit(np), file=filename, form='unformatted', iostat=ios)
            if (ios /= 0) then
               write(0,*) 'error opening file '//trim(filename)//' in timeseries_initialize'
               ierr = -1
               return
            end if

!! header
            write(tsunit(np)) solver_state%dt
            write(tsunit(np)) modelstarttime
            write(tsunit(np)) 0  ! avgyn
            write(tsunit(np)) 0  ! avgfrq
            write(tsunit(np)) 0  ! avglen
            write(tsunit(np)) 0  ! avgfirst
            write(tsunit(np)) var2d_number
            do nvar=1,var2d_number
               write(tsunit(np)) nvar,var2d_name(nvar),var2d_units(nvar),var2d_description(nvar)
               write(0,*)        nvar,var2d_name(nvar),var2d_units(nvar),var2d_description(nvar)
            end do
            write(tsunit(np)) var3d_number
            do nvar=1,var3d_number
               write(tsunit(np)) nvar,var3d_name(nvar),var3d_lvlind(nvar),var3d_units(nvar),var3d_description(nvar)
               write(0,*)        nvar,var3d_name(nvar),var3d_lvlind(nvar),var3d_units(nvar),var3d_description(nvar)
            end do
            write(tsunit(np)) 1
            write(tsunit(np)) ipnt(np),jpnt(np), points_lat(np),points_lon(np)
            write(tsunit(np)) lm+1            ! number of full levels
            do k=1,lm+1
               write(tsunit(np)) zint(np,k)
            end do
            write(tsunit(np)) lm              ! number of half levels
            do k=1,lm
             write(tsunit(np)) zmid(np,k)
            end do
!! end of header

!! close the file
            close(unit=tsunit(np), iostat=ios)
            if (ios /= 0) then
               write(0,*) 'error closing '//filename
            end if

         end if inside_if

      end do np_loop

      initialized=.true.

      end subroutine timeseries_initialize

!-----------------------------------------------------------------------
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!-----------------------------------------------------------------------

      subroutine timeseries_run(solver_state,ntsd,ierr)

      use module_solver_internal_state, only : solver_internal_state

      implicit none

      type(solver_internal_state),intent(in) :: solver_state
      integer,intent(in) :: ntsd
      integer, intent(out) :: ierr

      integer i,j,k, lm, lmax, l
      integer n, np, var
      integer year, month, day, hour, minute, second, ten_thousandth
      real :: tlm0d,tph0d
      real :: pus,pvs

      integer :: ihr
      integer :: nf_hours,nf_minutes
      real    :: nf_seconds
      real    :: secfcst

      integer :: idatv(3),ihrv,idaywk

      integer :: ios
      character*128  :: filename

!!
      if (.not.nml_exist) then
        ierr = 0
        return
      end if

      tlm0d = solver_state%tlm0d
      tph0d = solver_state%tph0d
      lm = solver_state%lm

!! calculate forecast time for this time step

      secfcst = solver_state%dt * ntsd
      nf_hours = int(secfcst/3600)
      ihr = nf_hours
      nf_minutes = int( mod(secfcst,3600.0)/60.0 )
      nf_seconds = (secfcst - nf_hours*3600.0) - nf_minutes*60.0
      if (nf_seconds< 0.000001) nf_seconds=0.0
      ten_thousandth = nint((secfcst-int(secfcst))*10000)

      call valid(solver_state%idat,solver_state%ihrst,nf_hours,idatv,ihrv,idaywk)

      year = idatv(3)
      month = idatv(2)
      day = idatv(1)
      hour = ihrv
      minute = nf_minutes
      second = nf_seconds


      j_loop: do j = jts, jte
      i_loop: do i = its, ite

        np_loop: do np = 1, npoints

          if (i.eq.ipnt(np) .and. j.eq.jpnt(np)) then

            tsunit(np) = 90 + np
            write(filename,fmt='(a,i2.2,a,i2.2,a)') "ts_p",np,"_d",my_domain_id,".bin"
            write(0,*) ' open file:', np,tsunit(np),filename
            open(unit=tsunit(np), file=filename, form='unformatted', position='append', iostat=ios)
            if (ios /= 0) then
               write(0,*) 'error opening file '//trim(filename)//' in timeseries_run'
               ierr = -1
               return
            end if

            write(tsunit(np)) year,month,day,hour,minute,second,ten_thousandth,ntsd

! 2d
            do n=1,solver_state%num_vars
              if (solver_state%vars(n)%tseries .and. solver_state%vars(n)%tkr==tkr_r2d) then
                write(tsunit(np)) solver_state%vars(n)%r2d(i,j)
              end if
            end do

!3d
            do n=1,solver_state%num_vars
              if (solver_state%vars(n)%tseries .and. solver_state%vars(n)%tkr==tkr_r3d) then
                lmax = lm
                do l=1,max_fulllevel_vars
                  if (trim(fulllevel_vars(l))==trim(solver_state%vars(n)%vbl_name)) then
                    lmax = lm+1
                    exit
                  end if
                end do
                do k=1,lmax
                  write(tsunit(np)) solver_state%vars(n)%r3d(i,j,k)
                end do
                if(lmax==lm) write(tsunit(np)) 0.0
              end if
            end do

            close(unit=tsunit(np), iostat=ios)
            if (ios /= 0) then
               write(0,*) 'error closing '//filename
            end if

          end if

        end do np_loop

      end do i_loop
      end do j_loop

      end subroutine timeseries_run

!-----------------------------------------------------------------------
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!-----------------------------------------------------------------------

      subroutine timeseries_finalize
      end subroutine timeseries_finalize

!-----------------------------------------------------------------------
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!-----------------------------------------------------------------------

      subroutine tll(almd,aphd,tlmd,tphd,tph0d,tlm0d)
!-----------------------------------------------------------------------
      real, intent(in) :: almd, aphd
      real, intent(out) :: tlmd, tphd
      real, intent(in) :: tph0d, tlm0d
!-----------------------------------------------------------------------
      real, parameter :: pi=3.141592654
      real, parameter :: dtr=pi/180.0
!
      real :: tph0, ctph0, stph0, relm, srlm, crlm
      real :: aph, sph, cph, cc, anum, denom
!-----------------------------------------------------------------------

      tph0=tph0d*dtr
      ctph0=cos(tph0)
      stph0=sin(tph0)

      relm=(almd-tlm0d)*dtr
      srlm=sin(relm)
      crlm=cos(relm)
      aph=aphd*dtr
      sph=sin(aph)
      cph=cos(aph)
      cc=cph*crlm
      anum=cph*srlm
      denom=ctph0*cc+stph0*sph

      tlmd=atan2(anum,denom)/dtr
      tphd=asin(ctph0*sph-stph0*cc)/dtr

      return

      end subroutine tll

!-----------------------------------------------------------------------
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!-----------------------------------------------------------------------

      subroutine ijnrs(tlmpt,tphpt,dlmd,dphd,wbd,sbd,im,jm,inrs,jnrs,inside)

!     ******************************************************************
!     *                                                                *
!     *  code to:                                                      *
!     *  - find the i and j indices of nearest h point of the b        *
!     *    grid box containing point tlmpt,tphpt                       *
!     *                                                                *
!     *  ictp version - d.jovic                                        *
!     ******************************************************************

      real, intent(in) :: tlmpt,tphpt,dlmd,dphd,wbd,sbd
      integer, intent(in) :: im,jm
      integer, intent(out) :: inrs,jnrs
      logical, intent(out) :: inside

!-----------------------------------------------------------------------

      inrs=nint((tlmpt-wbd)/dlmd)+1
      jnrs=nint((tphpt-sbd)/dphd)+1

      if (jnrs >= 1 .and. jnrs <= jm .and. inrs >= 1 .and. inrs <= im ) then
         inside=.true.
      else
         inside=.false.
      endif

      return

      end subroutine ijnrs

!-----------------------------------------------------------------------
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!-----------------------------------------------------------------------

      subroutine tllwin(almd,aphd,tph0d,tlm0d,pus,pvs)

!     ******************************************************************
!     *                                                                *
!     *  ll to tll transformation of velocity                          *
!     *                                                                *
!     *  programer: z.janjic, yugoslav fed. hydromet. inst.,           *
!     *                 beograd, 1982                                  *
!     *                                                                *
!     ******************************************************************

      real, intent (in) :: almd, aphd
      real, intent (in) :: tph0d, tlm0d
      real, intent (inout) :: pus,pvs

      real, parameter :: pi=3.141592654
      real, parameter :: dtr=pi/180.0

      real :: relm,ctph0,stph0,tph0,srlm,crlm,ph,sph,cph,cc,tph
      real :: rctph,cray,dray
      real :: tpus,tpvs
      real :: rdenom

!-----------------------------------------------------------------------

      tph0=tph0d*dtr

      ctph0=cos(tph0)
      stph0=sin(tph0)

         relm=(almd-tlm0d)*dtr
         srlm=sin(relm)
         crlm=cos(relm)

         ph=aphd*dtr
         sph=sin(ph)
         cph=cos(ph)

         cc=cph*crlm
         tph=asin(ctph0*sph-stph0*cc)
         rctph=1.0/cos(tph)

         cray=stph0*srlm*rctph
         dray=(ctph0*cph+stph0*sph*crlm)*rctph

         tpus=pus
         tpvs=pvs
         rdenom=1.0/(dray*dray + cray*cray)
         pus=(dray*tpus+cray*tpvs)*rdenom
         pvs=(-cray*tpus+dray*tpvs)*rdenom

      return

      end subroutine tllwin

!-----------------------------------------------------------------------
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!-----------------------------------------------------------------------

      subroutine valid(idat,ihrst,ihr,idatv,ihrv,idaywk)
      integer idat(3),ihrst,ihr,idatv(3),ihrv,idaywk

      integer :: ijulian,iadd,iadday,idayyr

      ijulian=iw3jdn(idat(3),idat(2),idat(1))
      iadd=ihrst+ihr
      iadday=int((ihrst+ihr)/24)
      ihrv=iadd-24*iadday
      ijulian=ijulian+iadday
      call w3fs26(ijulian,idatv(3),idatv(2),idatv(1),idaywk,idayyr)
      return
      end subroutine valid

!-----------------------------------------------------------------------
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!-----------------------------------------------------------------------

      integer  function iw3jdn(iyear,month,iday)

      integer, intent(in) :: iyear,month,iday

      iw3jdn  = iday - 32075                                           &
              + 1461 * (iyear + 4800 + (month - 14) / 12) / 4          &
              + 367 * (month - 2 - (month -14) / 12 * 12) / 12         &
              - 3 * ((iyear + 4900 + (month - 14) / 12) / 100) / 4
      return
      end function iw3jdn

!-----------------------------------------------------------------------
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!-----------------------------------------------------------------------

      subroutine w3fs26(jldayn,iyear,month,iday,idaywk,idayyr)

      integer, intent(in) :: jldayn
      integer, intent(out) :: iyear, month, iday,idaywk,idayyr

      integer :: l,n,i,j

      l      = jldayn + 68569
      n      = 4 * l / 146097
      l      = l - (146097 * n + 3) / 4
      i      = 4000 * (l + 1) / 1461001
      l      = l - 1461 * i / 4 + 31
      j      = 80 * l / 2447
      iday   = l - 2447 * j / 80
      l      = j / 11
      month  = j + 2 - 12 * l
      iyear  = 100 * (n - 49) + i + l
      idaywk = mod((jldayn + 1),7) + 1
      idayyr = jldayn - (-31739 +1461 * (iyear+4799) / 4 - 3 * ((iyear+4899)/100)/4)
      return
      end subroutine w3fs26

!-----------------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!-----------------------------------------------------------------------

      end module module_timeseries

!-----------------------------------------------------------------------
