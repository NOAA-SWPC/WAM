subroutine datetime
character (24) :: DateT
character ( 8) :: date
character (10) :: time
character ( 5) :: zone
character ( 3) :: month(12)
character (80) :: FMT
integer                    :: values(8)
data month /'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'/

call date_and_time(date,time,zone,values)
FMT = "(a3,i3,',',i5,i3,':',i2.2,':',i2.2)"
write(DateT,FMT) month(values(2)),values(3),values(1),values(5),values(6),values(7)
print "(' DATE-TIME:  ',A)",DateT

end subroutine datetime

