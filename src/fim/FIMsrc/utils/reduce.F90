program reduce
!*********************************************************************
!This program reads in one or more FIM pressure-level output files and 
! creates new files with one selected pressure level per file.
!The pressure levels and file names are set in the file REDUCEinput.
!There are only 5 pressure level variables - see PLvars below.
!Before output the files are changed to fixed ij grid order using inv_perm.
!Jacques Middlecoff December 2009.
!*********************************************************************
implicit none
integer,parameter   :: NPLvars = 5 !Number of pressure level variables
character(200)      :: FullName    !PathName+FileName
character(120)      :: PathName    !Path name where the pressure level files reside
character(80)       :: FileName    !Name of the pressure level file to be read
character(80)       :: header(10)  !FIM output file header
character(12)       :: yyyymmddhhmm! Forecast initial time
character(6)        :: FileTime    !FIM variable name read from the FIM header
character(4)        :: VarName     !FIM variable name read from the FIM header
character(4)        :: PLvars(NPLvars)=(/'up3P','vp3P','tmpP','rp3P','hgtP'/)
character(2)        :: Var(NPLvars)!Name of the FIM out variables to process 
integer             :: Nvars       !Number of FIM out variables to process (max 5)
integer             :: Ntimes      !Number of FIM output times to process
integer             :: NumLevels   !The number of pressure levels to output
integer             :: Dim1        !FIM vertical dimension read from the header
integer             :: nip         !FIM # of icosahedral points read from the header
integer             :: unit=50     !I/O unit
integer             :: n,level,ivl !Indexes
integer             :: ipn,t,v     !Indexes
integer,allocatable :: Time    (:) !The pressure levels (mb) for output
integer,allocatable :: Pressure(:) !The pressure levels (mb) for output
integer,allocatable :: inv_perm(:) !Permutation array for fixed grid order
real   ,allocatable :: invar (:,:) !Location for FIM variable to be read in
real   ,allocatable :: outvar(:,:) !Output variable before inv_perm
real   ,allocatable :: fixvar  (:) !Output variable after inv_perm

OPEN (10,file="REDUCEinput")
read (10,*) PathName
read (10,*) Nvars
read (10,*) (Var(n),n=1,Nvars)
read (10,*) Ntimes
allocate(   Time(Ntimes))
read (10,*) Time
read (10,*) NumLevels
allocate(   Pressure(NumLevels))
read (10,*) Pressure
close(10)
write(FileTime,"(i6.6)") Time(1)
FileName =  'fim_out_' // trim(Var(1)) // FileTime
FullName = trim(PathName) // trim(FileName)
open(unit,file=FullName,form="unformatted")
read(unit) header
read(header(3),"(5x,i2,6x,i10)") Dim1,nip
close(unit)
allocate(inv_perm(nip))
call GetInvPerm(nip,inv_perm)
do v=1,Nvars
  do t=1,Ntimes
    write(FileTime,"(i6.6)") Time(t)
    FileName =  'fim_out_' // trim(Var(v)) // FileTime
    FullName = trim(PathName) // trim(FileName)
    open(unit,file=FullName,form="unformatted")
    read(unit) header
    read(header,"(4x,A,37x,A)") VarName,yyyymmddhhmm
    read(header(3),"(5x,i2,6x,i10)") Dim1,nip
    if(NumLevels > Dim1) then
      print"('NumLevels=',I0,' > the number of levels Dim1=',I0)",NumLevels,Dim1
      print*,'In the file ',trim(Filename)
      print*,'Fatal error'
      stop
    endif
    do n=1,NPLvars
      if(VarName == PLvars(n)) then
        exit
      elseif(n==NPLvars) then
        print*,'Error: ',VarName,' is not a pressure level variable'
        stop
      endif
    enddo  
    allocate(invar(Dim1,nip))
    read (unit) invar
    close(unit)
    allocate(outvar(nip,NumLevels),fixvar(nip))
    level = 1
    do ivl=1,Dim1
      if( 1000-(ivl-1)*25 == Pressure(level) ) then
        outvar(:,level) = invar(ivl,:)
        level=level+1
        if(level > NumLevels) exit
      endif
    enddo
    if(level /= NumLevels+1) then
      print*,'Error: did not find correct pressure levels'
      stop
    endif
    do level = 1,NumLevels
      do ipn=1,nip
        fixvar(ipn) = outvar(inv_perm(ipn),level)
      enddo
      write(FullName,"(a10,'_',a6,'_',i4.4,'_',a8,'_',a2)") FileName(1:10),FileName(11:16),Pressure(level), &
           yyyymmddhhmm(1:8),yyyymmddhhmm(9:10)
      write(header(10),"('Pressure level = ',i4.4,' mb')") Pressure(level)
      open (unit,file=FullName,form="unformatted")
      write(unit) header
      write(unit) fixvar
    enddo
    deallocate(invar,outvar,fixvar)
  enddo !Ntimes
enddo !Nvars

stop
end program reduce

subroutine GetInvPerm(nip,inv_perm)
implicit none
integer,parameter   :: npp=6	        ! number of proximity points(max)
integer,intent(IN)  :: nip
integer,intent(OUT) :: inv_perm(nip)
real                ::   work2d(nip)
integer             ::  iwork2d(nip)
integer             :: isn,i2,i3
character(16)       :: header

OPEN(10, file='icos_grid_info_level.dat',form='unformatted')
READ(10) header
READ(10) header
READ(10) work2d,work2d
do isn = 1,npp
  READ(10) iwork2d
enddo
READ(10) iwork2d
do i3 = 1,2
  do i2 = 1,2
    do isn = 1,npp
      READ(10) work2d
    enddo
  enddo
enddo
READ(10) inv_perm
CLOSE(10)

return
end subroutine GetInvPerm
