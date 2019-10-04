program readINI
  implicit none
  integer,parameter :: npp=6
  integer,parameter :: nd=2	! number of directions (x,y)
  integer,parameter :: glvl=5
  integer,parameter :: nip=10*(2**glvl)**2+2
  integer , parameter :: nvl=25	        ! number of vertical levels
  integer , parameter :: nvlp1=nvl+1      ! number of vertical levels plus one
  real us3d(nvl,nip)	! zonal wind (m/s), layer
  real vs3d(nvl,nip)	! meridional wind (m/s), layer
  real dp3d(nvl,nip)	! del p between coord levels (pascals)
  real mp3d(nvl,nip)	! Montgomery Potential (m**2/s**2)
  real th3d(nvl,nip)	! theta (k), layer
  real pr3d(nvlp1,nip)	! pressure (pascals), level
  real ex3d(nvlp1,nip)	! exner funciton, level
  real ph3d(nvlp1,nip)	! phi (=gz), m**2/s**2
  real phse(npp,nip)	! phi bottom interpolated to the edges
  real pref(nvlp1,nip)
  integer isn,ipn,ivl

  open(unit=28,file="fim_ini.dat",form="unformatted")
  read (28) us3d,vs3d,dp3d,pr3d,ex3d,mp3d,th3d,ph3d,phse,pref
  close(28)
  write(100,100) ((ivl,ipn,us3d(ivl,ipn),ivl=1,nvl),ipn=1,nip)
  write(101,100) ((ivl,ipn,vs3d(ivl,ipn),ivl=1,nvl),ipn=1,nip)
  write(102,100) ((ivl,ipn,dp3d(ivl,ipn),ivl=1,nvl),ipn=1,nip)
  write(103,100) ((ivl,ipn,pr3d(ivl,ipn),ivl=1,nvlp1),ipn=1,nip)
  write(104,100) ((ivl,ipn,ex3d(ivl,ipn),ivl=1,nvlp1),ipn=1,nip)
  write(105,100) ((ivl,ipn,mp3d(ivl,ipn),ivl=1,nvl),ipn=1,nip)
  write(106,100) ((ivl,ipn,th3d(ivl,ipn),ivl=1,nvl),ipn=1,nip)
  write(107,100) ((ivl,ipn,ph3d(ivl,ipn),ivl=1,nvlp1),ipn=1,nip)
  write(108,100) ((isn,ipn,phse(ivl,ipn),isn=1,npp),ipn=1,nip)
  write(109,100) ((ivl,ipn,pref(ivl,ipn),ivl=1,nvlp1),ipn=1,nip)
100 format(2i10,1pe15.7)

end program readINI
