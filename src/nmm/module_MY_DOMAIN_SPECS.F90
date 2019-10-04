!-----------------------------------------------------------------------
!
      MODULE MODULE_MY_DOMAIN_SPECS
!
!-----------------------------------------------------------------------
!
!***  Set and hold key domain/subdomain dimensions, the forecast task
!***  intracommunicator, and task rank within that communicator for
!***  the currently active domain.  
!
!-----------------------------------------------------------------------
!
      USE MODULE_KINDS
!
!-----------------------------------------------------------------------
!
      IMPLICIT NONE
!
!-----------------------------------------------------------------------
!
      PUBLIC
!
      INTEGER(kind=KINT),SAVE :: IDS,IDE,JDS,JDE                        &
                                ,IMS,IME,JMS,JME                        &
                                ,ITS,ITE,JTS,JTE                        &
                                ,ITS_B1,ITE_B1,ITS_B2,ITE_B2            &
                                ,ITS_B1_H1,ITE_B1_H1,ITE_B1_H2          &
                                ,ITS_B1_H2                              &
                                ,ITS_H1,ITE_H1,ITS_H2,ITE_H2            &
                                ,JTS_B1,JTE_B1,JTS_B2,JTE_B2            &
                                ,JTS_B1_H1,JTE_B1_H1,JTE_B1_H2          &
                                ,JTS_B1_H2                              &
                                ,JTS_H1,JTE_H1,JTS_H2,JTE_H2
!
      INTEGER(kind=KINT),SAVE :: IHALO,JHALO                            &
                                ,MPI_COMM_COMP                          &
                                ,MY_DOMAIN_ID                           &
                                ,MYPE
!
      INTEGER(kind=KINT),DIMENSION(1:8),SAVE :: MY_NEB
!
      INTEGER(kind=KINT),DIMENSION(:),ALLOCATABLE,SAVE :: LOCAL_ISTART  &
                                                         ,LOCAL_IEND    &
                                                         ,LOCAL_JSTART  &
                                                         ,LOCAL_JEND
!
      LOGICAL(kind=KLOG),SAVE :: ADV_STANDARD                           &
                                ,ADV_UPSTREAM                           &
                                ,E_BDY                                  &
                                ,N_BDY                                  &
                                ,S_BDY                                  &
                                ,W_BDY
!
!-----------------------------------------------------------------------
!
      CONTAINS
!
!-----------------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!-----------------------------------------------------------------------
!
      SUBROUTINE SET_DOMAIN_SPECS(ITS_IN,ITE_IN,JTS_IN,JTE_IN           &
                                 ,IMS_IN,IME_IN,JMS_IN,JME_IN           &
                                 ,IDS_IN,IDE_IN,JDS_IN,JDE_IN           &
                                 ,IHALO_IN,JHALO_IN                     &
                                 ,MY_DOMAIN_ID_IN                       &
                                 ,MYPE_IN                               &
                                 ,MY_NEB_IN                             &
                                 ,MPI_COMM_COMP_IN                      &
                                 ,NUM_PES                               &
!
                                 ,LOCAL_ISTART_IN,LOCAL_IEND_IN         &  !    ^
                                 ,LOCAL_JSTART_IN,LOCAL_JEND_IN         &  !    |
                                 ,ADV_STANDARD_IN,ADV_UPSTREAM_IN       &  ! Optional arguments
                                 ,S_BDY_IN,N_BDY_IN                     &  !    |
                                 ,W_BDY_IN,E_BDY_IN                     &  !    v
                                   )
!
!-----------------------------------------------------------------------
!
!***  Set these key domain-related variables for the active domain.
!
!-----------------------------------------------------------------------
!
!------------------------
!***  Argument variables
!------------------------
!
      INTEGER(kind=KINT),INTENT(IN) :: ITS_IN,ITE_IN,JTS_IN,JTE_IN     &
                                      ,IMS_IN,IME_IN,JMS_IN,JME_IN     &
                                      ,IDS_IN,IDE_IN,JDS_IN,JDE_IN     &
                                      ,IHALO_IN,JHALO_IN               &
                                      ,MPI_COMM_COMP_IN                &
                                      ,MY_DOMAIN_ID_IN                 &
                                      ,MYPE_IN                         &
                                      ,NUM_PES
!
      INTEGER(KIND=KINT),DIMENSION(1:8),INTENT(IN) :: MY_NEB_IN
!
!-----------------------
!*** Optional arguments
!-----------------------
!
      INTEGER(kind=KINT),DIMENSION(0:NUM_PES-1),INTENT(IN),OPTIONAL :: &
                                                      LOCAL_ISTART_IN  &
                                                     ,LOCAL_IEND_IN    &
                                                     ,LOCAL_JSTART_IN  &
                                                     ,LOCAL_JEND_IN
!
      LOGICAL(kind=KLOG),INTENT(IN),OPTIONAL :: ADV_STANDARD_IN        &
                                               ,ADV_UPSTREAM_IN        &
                                               ,E_BDY_IN,N_BDY_IN      &
                                               ,S_BDY_IN,W_BDY_IN
!
!---------------------
!***  Local variables
!---------------------
!
      INTEGER(kind=KINT) :: N
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
      ITS=ITS_IN
      ITE=ITE_IN
      IMS=IMS_IN
      IME=IME_IN
      IDS=IDS_IN
      IDE=IDE_IN
!
      JTS=JTS_IN
      JTE=JTE_IN
      JMS=JMS_IN
      JME=JME_IN
      JDS=JDS_IN
      JDE=JDE_IN
!
      IHALO=IHALO_IN
      JHALO=JHALO_IN
!
      MPI_COMM_COMP=MPI_COMM_COMP_IN
      MY_DOMAIN_ID=MY_DOMAIN_ID_IN
      MYPE=MYPE_IN
!
      ITS_B1=MAX(ITS,IDS+1)
      ITE_B1=MIN(ITE,IDE-1)
      ITS_B2=MAX(ITS,IDS+2)
      ITE_B2=MIN(ITE,IDE-2)
      ITS_B1_H1=MAX(ITS-1,IDS+1)
      ITE_B1_H1=MIN(ITE+1,IDE-1)
      ITE_B1_H2=MIN(ITE+2,IDE-1)
      ITS_H1=MAX(ITS-1,IDS)
      ITE_H1=MIN(ITE+1,IDE)
      ITS_H2=MAX(ITS-2,IDS)
      ITE_H2=MIN(ITE+2,IDE)
      JTS_B1=MAX(JTS,JDS+1)
      JTE_B1=MIN(JTE,JDE-1)
      JTS_B2=MAX(JTS,JDS+2)
      JTE_B2=MIN(JTE,JDE-2)
      JTS_B1_H1=MAX(JTS-1,JDS+1)
      JTE_B1_H1=MIN(JTE+1,JDE-1)
      JTE_B1_H2=MIN(JTE+2,JDE-1)
      JTS_H1=MAX(JTS-1,JDS)
      JTE_H1=MIN(JTE+1,JDE)
      JTS_H2=MAX(JTS-2,JDS)
      JTE_H2=MIN(JTE+2,JDE)
!
      DO N=1,8
        MY_NEB(N)=MY_NEB_IN(N)
      ENDDO
!
      IF(PRESENT(ADV_STANDARD_IN))THEN
        ADV_STANDARD=ADV_STANDARD_IN
        ADV_UPSTREAM=ADV_UPSTREAM_IN
      ENDIF
!
      IF(PRESENT(S_BDY_IN))THEN
        S_BDY=S_BDY_IN
        N_BDY=N_BDY_IN
        W_BDY=W_BDY_IN
        E_BDY=E_BDY_IN
      ENDIF
!
      IF(PRESENT(LOCAL_ISTART_IN))THEN
        IF(ALLOCATED(LOCAL_ISTART))THEN
          DEALLOCATE(LOCAL_ISTART)
          DEALLOCATE(LOCAL_IEND)
          DEALLOCATE(LOCAL_JSTART)
          DEALLOCATE(LOCAL_JEND)
        ENDIF
!
        ALLOCATE(LOCAL_ISTART(0:NUM_PES-1))
        ALLOCATE(LOCAL_IEND(0:NUM_PES-1))
        ALLOCATE(LOCAL_JSTART(0:NUM_PES-1))
        ALLOCATE(LOCAL_JEND(0:NUM_PES-1))
!
        DO N=0,NUM_PES-1
          LOCAL_ISTART(N)=LOCAL_ISTART_IN(N)
          LOCAL_IEND(N)  =LOCAL_IEND_IN(N)
          LOCAL_JSTART(N)=LOCAL_JSTART_IN(N)
          LOCAL_JEND(N)  =LOCAL_JEND_IN(N)
        ENDDO
      ENDIF
!
!-----------------------------------------------------------------------
!
      END SUBROUTINE SET_DOMAIN_SPECS
!
!-----------------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!-----------------------------------------------------------------------
!
      END MODULE MODULE_MY_DOMAIN_SPECS
!
!-----------------------------------------------------------------------
