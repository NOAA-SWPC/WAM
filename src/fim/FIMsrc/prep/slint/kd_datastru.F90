MODULE kd_datastru

     TYPE Node
       LOGICAL :: bucket   ! true if it is a bucket node
       INTEGER :: npoints  ! number of points of the node 
       INTEGER :: discrim  ! discrimator  
       REAL :: cutval      ! the cut value for the space
       INTEGER :: lopt     ! lower point  
       INTEGER :: hipt     ! high point  
       TYPE(Node) , POINTER :: loson 
       TYPE(Node) , POINTER :: hison 
     END TYPE Node

     INTEGER n, d, gc, bottom_level, dir, num_k

     TYPE(Node), TARGET, ALLOCATABLE :: nodes(:)
     INTEGER, ALLOCATABLE :: perm(:)
     INTEGER, TARGET, ALLOCATABLE :: nni(:)
     REAL, ALLOCATABLE :: nnd(:)
     REAL, ALLOCATABLE :: os(:)
     REAL, POINTER:: ppoints(:,:)
     TYPE(Node), POINTER :: root

     REAL :: curDistSq, cur_o_s 


END MODULE kd_datastru
