SUBROUTINE RESGKL_MAIN(MODE,ACCURACY,M,N,L,K,MATDESCRA,INDXA,PNTRBA,PNTREA,A)
  !$ use omp_lib
  IMPLICIT NONE

  INTEGER M, N, L, K, INDXA(*), PNTRBA(*), PNTREA(*)
  DOUBLE PRECISION A(*),cpbeta,te1,te2,te3
  CHARACTER*6 MATDESCRA
  CHARACTER MODE

  INTEGER INFO,COU,SELEK,ACCURACY
  INTEGER J,I,STARTT,TMPT,MINT,RATE
  INTEGER ISEED(4),ITR
  DOUBLE PRECISION MIN_ERR,TMP_ERR,ERR,CPBK(K,K),CPVK(N,K),CPUK(M,K),CPVPLUS(N),GKL_TIME,BK(K,K),VK(N,K),UK(M,K),VPLUS(N),beta
  EXTERNAL GKL,ERR
  ISEED( 1 ) = 12
  ISEED( 2 ) = 14
  ISEED( 3 ) = 13
  ISEED( 4 ) = 733
  print *, ",,,,,MODE = ",MODE
  WRITE (*,*)",,,,,M,N,K,L"
  WRITE (*,*)",,,,,",M,",",N,",",K,",",L
  !CALL SYSTEM_CLOCK(STARTT,RATE,TMPT)
  !$ STARTT = omp_get_wtime()
  !  CALL GKL(A,M,N,K,BK,VK,UK,VPLUS)
  CALL GKL(mode,MATDESCRA,INDXA,PNTRBA,PNTREA,A,M,N,K,BK,VK,UK,VPLUS,beta)
  !CALL SYSTEM_CLOCK(TMPT)
  !$ TMPT = omp_get_wtime()
  GKL_TIME = TMPT - STARTT
  DO SELEK = 1, 6!, 2
     WRITE(*,*) 
     CPBK=BK
     CPUK=UK
     CPVK=VK
     CPVPLUS=VPLUS
     cpbeta = beta
     MIN_ERR = -1.0
     COU = 0
     ITR = 1
     !CALL SYSTEM_CLOCK(STARTT)
     !$ STARTT = omp_get_wtime()
     DO WHILE (COU < 10)
        CALL RESGKL(mode,MATDESCRA,INDXA,PNTRBA,PNTREA,A,M,N,K,L,CPBK,CPVK,CPUK,CPVPLUS,cpbeta,INFO,SELEK)

        TMP_ERR=0.0D+0
        DO I = 1,L
           TMP_ERR=max(TMP_ERR,abs(cpbk(I,L+1)))
        END DO
        TE1 = ERR(mode,MATDESCRA,INDXA,PNTRBA,PNTREA,A,M,N,K,1,CPBK,CPVK,CPUK)
        TE2 = ERR(mode,MATDESCRA,INDXA,PNTRBA,PNTREA,A,M,N,K,2,CPBK,CPVK,CPUK)
        TE3 = ERR(mode,MATDESCRA,INDXA,PNTRBA,PNTREA,A,M,N,K,3,CPBK,CPVK,CPUK)

        ITR=ITR+1
        !CALL SYSTEM_CLOCK(TMPT)
        !$ tmpt = omp_get_wtime()
        IF ( (TMP_ERR < MIN_ERR) .OR. (MIN_ERR .EQ. -1.0) ) THEN
           ! ボトムを更新出来ないかぎり、COUはリセットされず増える
           MIN_ERR = TMP_ERR
           MINT = TMPT
           COU = -1
        END IF
        !WRITE(*,*) SELEK,",",ITR,",",((TMPT-STARTT)+GKL_TIME)/DBLE(RATE),",",TMP_ERR
        WRITE(*,*) SELEK,",",ITR,",",(TMPT-STARTT)+GKL_TIME,",",TMP_ERR,",",te1,",",te2,",",te3
        WRITE(*,*) SELEK,",",ITR,",",(TMPT-STARTT)+GKL_TIME,",",TMP_ERR,",",cpbk(1,1),",",cpbk(2,2),",",cpbk(3,3)
        IF(10**(0.0 - accuracy) > TMP_ERR) exit
        COU = COU + 1
     END DO
     !WRITE(*,*) SELEK,",",ITR-7,", ",((MINT-STARTT)+GKL_TIME)/DBLE(RATE),",",MIN_ERR
     DO I = 1,3
        WRITE (*,*) I, CPBK(I,I)
     END DO
  END DO

END SUBROUTINE RESGKL_MAIN
