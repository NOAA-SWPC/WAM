      subroutine GetJdate(yyyymmddhhmm,jdate)
      implicit none
      CHARACTER(len=12) :: yyyymmddhhmm
      CHARACTER(len=9) :: jdate 
      INTEGER year, month, day, hour, minute, jday, IW3JDN

      ! get date info from the date string
      READ(UNIT=yyyymmddhhmm(1:4), FMT='(I4)') year
      READ(UNIT=yyyymmddhhmm(5:6), FMT='(I2)') month
      READ(UNIT=yyyymmddhhmm(7:8), FMT='(I2)') day
      READ(UNIT=yyyymmddhhmm(9:10), FMT='(I2)') hour
      READ(UNIT=yyyymmddhhmm(11:12), FMT='(I2)') minute
     
      ! create the jdate string
      jday = IW3JDN(year,month,day) - IW3JDN(year,1, 1) + 1
      WRITE(UNIT=jdate(1:2), FMT='(I2.2)') MOD (year, 100) 
      WRITE(UNIT=jdate(3:5), FMT='(I3.3)') jday 
      WRITE(UNIT=jdate(6:7), FMT='(I2.2)') hour 
      WRITE(UNIT=jdate(8:9), FMT='(I2.2)') minute

      return
      end subroutine GetJdate
