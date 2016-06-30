SUBROUTINE RESGKL_MAIN(MODE,ACCURACY,M,N,L,K,MATDESCRA,INDXA,PNTRBA,PNTREA,A)
  !$ use omp_lib
  IMPLICIT NONE

  INTEGER M, N, L, K, INDXA(*), PNTRBA(*), PNTREA(*)
  DOUBLE PRECISION A(*),te1,te2,NRM,DNRM2,STARTT,TMPT,MINT
  CHARACTER*6 MATDESCRA
  CHARACTER MODE

  INTEGER INFO,COU,SELEK,ACCURACY
  INTEGER J,I,W
  INTEGER ISEED(4),ITR
  DOUBLE PRECISION MIN_ERR,TMP_ERR,ERR,BK(K,K),VK(N,K),UK(M,K),VPLUS(N),VPLUS_ORIG(N)
  EXTERNAL ERR,DNRM2
  ISEED( 1 ) = 12
  ISEED( 2 ) = 14
  ISEED( 3 ) = 13
  ISEED( 4 ) = 733
  print *, "MODE = ",MODE
  WRITE (*,*)"M N K L"
  WRITE (*,*)M,N,K,L
  DO W = 1,1
     CALL DLARNV(1,ISEED,N,VPLUS_ORIG)
     NRM = DNRM2(N,VPLUS_ORIG,1)
     VPLUS_ORIG = VPLUS_ORIG / NRM
     DO SELEK = 1, 5 , 2
        !$ STARTT = omp_get_wtime()
        !CALL SYSTEM_CLOCK(STARTT)
        WRITE(*,*) 
        BK=0
        UK=0
        VK=0
        MIN_ERR = -1.0
        COU = 0
        ITR = 1
        J = 0
        VPLUS=VPLUS_ORIG
        DO WHILE (COU < 100)
           CALL RESGKL(J,mode,MATDESCRA,INDXA,PNTRBA,PNTREA,A,M,N,K,L,BK,VK,UK,VPLUS,INFO,SELEK)
   
           TMP_ERR=0.0D+0
           DO I = 1,L
              TMP_ERR=max(TMP_ERR,abs(bk(I,L+1)))
           END DO
           ITR=ITR+1
           !CALL SYSTEM_CLOCK(TMPT)
           !$ tmpt = omp_get_wtime()
           IF ( (TMP_ERR < MIN_ERR) .OR. (MIN_ERR .EQ. -1.0) ) THEN
              MIN_ERR = TMP_ERR
              MINT = TMPT
              COU = -1
           END IF
           !WRITE(*,*) SELEK,ITR,((TMPT-STARTT)+GKL_TIME)/DBLE(RATE),TMP_ERR
           WRITE(*,*) SELEK,ITR,(TMPT-STARTT),TMP_ERR
           IF(10**(0.0 - accuracy) > TMP_ERR) exit
           COU = COU + 1
        END DO
        te2 = 0.0D+0
        DO I = 1,L
           TE1 = ERR(mode,MATDESCRA,INDXA,PNTRBA,PNTREA,A,M,N,K,I,BK,VK,UK)
           te2 = te2 + te1
           write(*,*) TE1
        END DO
        write(*,*)
        !write(*,*) W,"SUM_OF_ERR_1-L" , SELEK, te2
        WRITE(*,*) W,"TIME", SELEK,ITR,(TMPT-STARTT),TMP_ERR,"AVE",te2/DBLE(L)
        
        !WRITE(*,*) SELEK,ITR-7,((MINT-STARTT)+GKL_TIME)/DBLE(RATE),MIN_ERR
        DO I = 1,3
           WRITE (*,*) "SingularValue(1-3)",I, BK(I,I)
        END DO
     END DO
  END DO

END SUBROUTINE RESGKL_MAIN
