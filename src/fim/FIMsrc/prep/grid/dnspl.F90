!=============================================================
!
! Ning Wang, March 2007 
!
!=============================================================
   program dnspl
   IMPLICIT NONE
 
! command line argument count
   INTEGER :: iargc
   CHARACTER(len=128) :: org_file
   CHARACTER(len=1) :: glvls
   INTEGER :: nip, glvl, i, seq
   REAL*8 lat, lon, d2r

   REAL*8 ll(2), min_dist
   INTEGER nn(3)

   REAL*8, ALLOCATABLE :: icos_grid_h(:,:) ! grid location (ll)
   REAL*8, ALLOCATABLE :: icos_grid(:,:) ! grid location (ll)
   INTEGER, ALLOCATABLE ::seqs(:) ! seq numbers for the grid

   IF (iargc() .NE. 2 ) THEN
     WRITE(0,*) 'Usage: dnspl [original grid file] [glvl:4-8]'
     STOP
   END IF
   
   CALL getarg(1,org_file)
   CALL getarg(2,glvls)

   d2r = 4.0 * atan(1.0) / 180.0

   glvl = 7
   nip = 10 * (2 ** (2 * glvl)) + 2 

! init kd tree for the original high resolution grid.
   CALL init_kd_tree2(org_file, nip, 3)

! allocate the data structure for original high resolution grid.
   ALLOCATE(icos_grid_h(nip, 2))	

! open and read in the original high resolution grid file.
   OPEN(10,file=org_file)
   READ(10,*) nip
   DO i = 1, nip 
     READ(10,*) seq, lat, lon
     icos_grid_h(i,1) = lat 
     icos_grid_h(i,2) = lon
   END DO
   CLOSE(10)

! allocate the data structure for low resolution grid.
   READ(glvls, *) glvl
   nip = 10 * (2 ** (2 * glvl)) + 2 
   ALLOCATE(icos_grid(nip, 2))	
   ALLOCATE(seqs(nip))	

! open and read in the low resolution grid file.
   OPEN               (10,"icos_grid_level.dat")
   call TestGlvlHeader(10,'icos_grid_level.dat','dnspl',glvl)
   READ(10,*) nip
   DO i = 1, nip 
     READ(10,*) seqs(i), lat, lon
     icos_grid(i,1) = lat 
     icos_grid(i,2) = lon
   END DO
   CLOSE(10)

   DO i = 1, nip
     ll(1) = icos_grid(i,1) * d2r 
     ll(2) = icos_grid(i,2) * d2r
     CALL knn_search(ll, nn, min_dist)
     icos_grid(i,1) = icos_grid_h(nn(1),1)
     icos_grid(i,2) = icos_grid_h(nn(1),2)
   END DO

! open and write in the down sampled low resolution grid file.
   OPEN(10,file="icos_grid_level.dat")
   call WriteGlvlHeader(10,glvls)
   WRITE(10,*) nip
   DO i = 1, nip 
     lat = icos_grid(i,1)  
     lon = icos_grid(i,2) 
     WRITE(10,*) seqs(i), lat, lon
   END DO
   CLOSE(10)

   END program dnspl
       

