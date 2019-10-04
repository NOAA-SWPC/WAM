program readcase
implicit none
integer , parameter :: glvl=5 ! the grid level
integer , parameter :: nip=10*(2**glvl)**2+2 ! # of icosahedral points
integer , parameter :: nvl=25	         ! number of vertical levels
integer , parameter :: nvlp1=nvl+1     ! number of vertical levels
real    :: us1(nvl  ,nip),vs1 (nvl,nip),dp1(nvl,nip)
real    :: pr1(nvlp1,nip),th1 (nvl,nip)
real    :: mp1(nvl  ,nip),tk1 (nvl,nip)
real    :: ph1(nvlp1,nip),vor1(nvl,nip),div1(nvl,nip)
real    :: qv1(nvl  ,nip),rh1 (nvl,nip),rn2d1(nip),pw2d1(nip)

real    :: us2(nvl  ,nip),vs2 (nvl,nip),dp2(nvl,nip)
real    :: pr2(nvlp1,nip),th2 (nvl,nip)
real    :: mp2(nvl  ,nip),tk2 (nvl,nip)
real    :: ph2(nvlp1,nip),vor2(nvl,nip),div2(nvl,nip)
real    :: qv2(nvl  ,nip),rh2 (nvl,nip),rn2d2(nip),pw2d2(nip)

integer :: lunin1=13
integer :: lunin2=14

read(lunin1)its1,us1,vs1,dp1,pr1,mp1,th1,vor1,ph1,qv1,rh1,rn2d1,pw2d1
read(lunin2)its2,us2,vs2,dp2,pr2,mp2,th2,vor2,ph2,qv2,rh2,rn2d2,pw2d2

print*,its1,its2
do ipn=1,nip
  do ivl=1,nvl
    if(us1(ivl,ipn) /= us2(ivl,ipn)) then
      print"(i8,2i5,1p2e20.7)",its1,ivl,ipn,us1(ivl,ipn),us2(ivl,ipn)
    endif
  enddo
enddo

