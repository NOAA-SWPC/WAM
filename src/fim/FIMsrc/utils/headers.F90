!These are the four routines that write, and verify the headers in the *.dat file.
!The I/O in all four routines must match.

subroutine WriteGlvlHeader(unit,glvl)
integer     ,intent(IN ) :: unit
integer     ,intent(IN ) :: glvl
character(16)            :: header
write(header,"('glvl =',I2)") glvl
write(unit) header
return
end subroutine WriteGlvlHeader

subroutine WriteCurveHeader(unit,curve)
integer     ,intent(IN ) :: unit
integer     ,intent(IN ) :: curve
character(16)            :: header
write(header,"('curve =',I2)") curve
write(unit) header
return
end subroutine WriteCurveHeader

subroutine TestGlvlHeader(unit,FileName,RoutineName,glvl)
integer     ,intent(IN) :: unit
character(*),intent(IN) :: FileName
character(*),intent(IN) :: RoutineName
integer     ,intent(IN) :: glvl
integer                 :: glvlHeader
character(16)           :: header
character(80)           :: FMT
integer                 :: ioerr
FMT ="('Error in ',a,' unit ',i0,' glvl=',i0,' does not match header glvl=',i0)"
read(unit, iostat=ioerr) header
if (ioerr /= 0) then
  write(6,*) 'testglvlheader: bad attempt to read header info file=', trim(filename)
  stop
end if
read(header,"(6x,I2)") glvlHeader
if(glvl /= glvlHeader) then
  write(6,FMT) RoutineName,unit,glvl,glvlHeader
endif
end subroutine TestGlvlHeader

subroutine TestCurveHeader(unit,FileName,RoutineName,curve)
integer     ,intent(IN) :: unit
character(*),intent(IN) :: FileName
character(*),intent(IN) :: RoutineName
integer     ,intent(IN) :: curve
integer                 :: curveHeader
character(16)           :: header
character(80)           :: FMT
integer                 :: ioerr
FMT ="('Error in ',a,' unit ',i0,' curve=',i0,' does not match header curve=',i0)"
read(unit, iostat=ioerr) header
if (ioerr /= 0) then
  write(6,*) 'testcurveheader: bad attempt to read header info file=', trim(filename)
  stop
end if
read(header,"(7x,I2)") curveHeader
if(curve /= curveHeader) then
  write(6,FMT) RoutineName,unit,curve,curveHeader
endif
end subroutine TestCurveHeader
