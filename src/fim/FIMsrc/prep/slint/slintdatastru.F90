MODULE SlintDataStru

     TYPE GRID
       INTEGER :: type
       INTEGER :: ngp, mx, my
       REAL*8, ALLOCATABLE :: latlon(:,:)
       REAL*8, ALLOCATABLE :: coeffs(:,:)
       REAL*4, ALLOCATABLE :: data(:)
       INTEGER, ALLOCATABLE :: nn(:,:)
       
     END TYPE GRID

     TYPE(GRID) src_grid, tgt_grid

END MODULE SlintDataStru


