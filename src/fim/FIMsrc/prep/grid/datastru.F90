MODULE DataStru

     INTEGER :: ml, sl, offset
     INTEGER :: init_a_sn, init_b_sn, init_c_sn
     INTEGER :: gc, xdim, ydim, orient

     TYPE GridPoint
       REAL*8 :: latlon(2)
       INTEGER :: seqnum
     END TYPE GridPoint

     TYPE GridPointWnb
       REAL*8 :: latlon(2)
       REAL*8 :: neighbor(2,6)
       INTEGER :: seqnum
     END TYPE GridPointWnb
  
     TYPE GridPointXYZ
       REAL*8 :: x
       REAL*8 :: y
       REAL*8 :: z
       INTEGER :: seqnum
     END TYPE GridPointXYZ
  
     TYPE Triangle
       TYPE(GridPoint) :: vertex(3)
       INTEGER start_idx, end_idx
     END TYPE Triangle

     TYPE(GridPoint), TARGET, ALLOCATABLE :: gp(:)
     INTEGER, ALLOCATABLE :: perm(:)
     INTEGER, ALLOCATABLE :: inv_perm(:)

END MODULE DataStru


