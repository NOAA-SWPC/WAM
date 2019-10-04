program jdate
implicit none
CHARACTER(len=12) :: yyyymmddhhmm
CHARACTER(len=9) :: JulianDate 

call getarg(1,yyyymmddhhmm)
call GetJdate(yyyymmddhhmm,JulianDate)
print'(A9,$)',JulianDate
stop
end program jdate
