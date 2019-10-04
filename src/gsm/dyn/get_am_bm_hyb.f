      SUBROUTINE AM_BM_hyb

! program log:
! 20110220:  Henry Juang updated code for NDSL advection
!
      USE gfs_dyn_MACHINE , ONLY : kind_grid
      use gfs_dyn_resol_def
      use gfs_dyn_coordinate_def		! hmhj
      use namelist_dynamics_def , only : ref_temp,
     &                                   semi_implicit_temp_profile
      use gfs_dyn_layout1, only : me
      use gfs_dyn_physcons, rd => con_rd
     &                    , cp => con_cp, rearth => con_rerth
      IMPLICIT NONE 

      REAL(KIND=KIND_EVOD)
     &     PK5REF(LEVP1),BETA,DPKREF(LEVS),
     &     TREF(levs),PSREF,KAPPA,FACTOR,
     &     ALFAREF(LEVS),
     &     VECM(LEVS),   YECM(LEVS,LEVS),TECM(LEVS,LEVS)
      INTEGER K,KK,J,IROW,ICOL,ICOLBEG,ICOLEND
 
!     print *,' enter get_am_bm_hyb_fd ' 		! hmhj

      if( semi_implicit_temp_profile ) then
        if (me == 0)
     &    print *,' use layer mean temperature for semi-implicit '
      else
        if (me == 0)
     &  print *,' use constant temperature for semi-implicit '
        do k=1,levs
          TREF(k) = ref_temp
        enddo
      endif
      PSREF = 80.
      BETA  = 1.
      KAPPA = RD/CP
 
!print PRINT *,' BEGIN  AM_BM '
 
      DO K=1,LEVP1
         PK5REF(K) = AK5(K) + BK5(K)*PSREF
!print PRINT 100,K,AK5(K),BK5(K),PK5REF(K)
      ENDDO
 
100   FORMAT('K=',I2,2X,'AK5=',E10.3,2X,'BK5=',E10.3,2X,'PK5REF=',E10.3)
!printPRINT*,'-------------------------------------------------------'
 
      DO K=1,LEVS
        DPKREF(K)  = PK5REF(K+1) - PK5REF(K)
        tor_hyb(K) = BETA*RD*tref(k)/(rearth*rearth)
!print PRINT 110,K,DPKREF(K),tor_hyb(K)
      ENDDO
 
110   FORMAT('K=',I2,2X,' IN AM_BM DPKREF=',E11.4,2x,'tor_hyb=',E11.4)
!printPRINT*,'-------------------------------------------------------'
!printPRINT*,'-  CALCULATE ALFAREF  WATCH ALFAREF(1)     '
 
      ALFAREF(1) = LOG(2.) ! COULD ALSO BE=1.  but watch for layer values
 
!printPRINT*,'-------------------------------------------------------'
      DO K=2,LEVS
       ALFAREF(K) = 1.
     &            - (PK5REF(K)/DPKREF(K))*LOG(PK5REF(K+1)/PK5REF(K))

!printPRINT 210,K,K,K,K,K
210   FORMAT('ALFA(',I2,')=1.-(PK5(',I2,')/DPK(',I2,'))*LOG(PK5(',I2,
     & '+1)/PK5(',I2,'))')
      ENDDO
 
!sela PRINT 125,ALFAREF(1)
125   FORMAT('WORRY --- ALFAREF(1)=',E10.3)
      DO K=1,LEVS
!print PRINT 130,K,ALFAREF(K)
      ENDDO
130   FORMAT('K=',I2,2X,'ALFAREF',E16.8)
!printPRINT*,'---- BEGIN MATRICES COMPUTATION -----------'
 
!     PRINT*,'333333333333333333333333333333333333333333'
!printPRINT 144
144   FORMAT('BEGIN YECM COMPUTATION')

       YECM = 0.
       DO IROW=1,LEVS
         YECM(IROW,IROW)=ALFAREF(IROW)*RD
         ICOLBEG=IROW+1
         IF(ICOLBEG.LE.LEVS)THEN
           DO ICOL=ICOLBEG,LEVS
             YECM(IROW,ICOL) = RD*LOG( PK5REF(ICOL+1)/PK5REF(ICOL) )
           ENDDO
         ENDIF
       ENDDO

150    FORMAT('YECM(',I2,',',I2,')=RD*LOG( PK5REF(',I2,
     &        '+1)/PK5REF(',I2,'))')
c     PRINT*,'-----------------1234567------------------'
160    FORMAT('YECM=',4(1X,E10.3))
 
      TECM = 0.
 
      DO IROW=1,LEVS
!       PRINT*,' DOING ROW ...............................',IROW

         TECM(IROW,IROW) = KAPPA*TREF(irow)*ALFAREF(IROW)
         ICOLEND = IROW - 1
 
        DO ICOL=1,ICOLEND
          FACTOR = (KAPPA*TREF(irow)/
     &                   DPKREF(IROW))*LOG(PK5REF(IROW+1)/PK5REF(IROW))
          TECM(IROW,ICOL) = FACTOR*DPKREF(ICOL)
        ENDDO
      ENDDO

165    FORMAT('IROW=',I2,2X,'FACTOR=',E16.8,2X,'ICOLEND=',I2)
166    FORMAT('FACTOR=(KAPPA*TREF/DPKREF(',I2,'))*LOG(PK5REF(',I2,
     & '+1)/PK5REF(',I2,'))')
167    FORMAT('INNERLUP IROW=',I2,2X,'ICOL=',I2,2X,'TECM(IR,IC)=',E12.4)
!     PRINT*,'4444444  PRINT YECM      44444444444444444'
 
!      DO IROW=1,LEVS
!       PRINT*,'YECM ROW=',IROW,'LEVS=',LEVS
!       PRINT 1700,(YECM(IROW,J),J=1,LEVS/2)
!       PRINT 1701,(YECM(IROW,J),J=LEVS/2+1,LEVS)
!      ENDDO
1700   FORMAT('  A  ',10(1X,E10.3))
1701   FORMAT('  B  ',10(1X,E10.3))
 
!     PRINT*,'5555555  PRINT TECM      55555555555555555'
 
!      DO IROW=1,LEVS
!       PRINT*,'TECM ROW=',IROW,'LEVS=',LEVS
!       PRINT 1700,(TECM(IROW,J),J=1,LEVS/2)
!       PRINT 1701,(TECM(IROW,J),J=LEVS/2+1,LEVS)
!      ENDDO
 
!     PRINT*,'666666666666666666666666666666666666666666'
!printPRINT 171
171   FORMAT('BEGIN VVEC DCOMPUTATION')

       DO ICOL=1,LEVS
         VECM(ICOL) = DPKREF(ICOL)/PSREF
       ENDDO

!      DO ICOL=1,LEVS
!print  PRINT 175,ICOL,VECM(ICOL)
!      ENDDO
175    FORMAT('ICOL=',I2,2X,'VECM=',E16.8)
 
 
 
      DO J=1,LEVS
       SVHYB(J) = VECM(LEVS+1-J)
        DO K=1,LEVS
          AMHYB(K,J) = YECM(LEVS+1-K,LEVS+1-J)
          BMHYB(K,J) = TECM(LEVS+1-K,LEVS+1-J)
        ENDDO
      ENDDO
 
      DO J=1,LEVS
        DO K=1,LEVS
         AMHYB(K,J) = AMHYB(K,J)*BETA/(rearth*rearth)
        ENDDO
      ENDDO
180    FORMAT('AMHYB=',4(1X,E10.3))
!     PRINT*,'777777777777777777777777777777777777777777'
 
185    FORMAT('BMHYB=',4(1X,E10.3))
 
!     PRINT*,'888888888888888888888888888888888888888888'
!      DO k=1,LEVS
!print  PRINT 186,K,SVHYB(k)
!      ENDDO
186    FORMAT('K=',I2,2X,' IN AM_BM SV= SVHYB=',F7.4)
!print PRINT*,'FIN SHALOM AM_BM)'

!     print *,' leave get_am_bm_hyb_fd ' 		! hmhj
 
      RETURN
      END
