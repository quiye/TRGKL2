SUBROUTINE RESGKL_MAIN(MODE,ACCURACY,M,N,L,K,IAP,JA,A,WORK,LWORK)
  use omp_lib
  IMPLICIT NONE

  INTEGER M, N, L, K, IAP(*), JA(*), LWORK
  DOUBLE PRECISION A(*),TE0,STARTT,TMPT,MINT,WORK(*)
  CHARACTER MODE
  INTEGER INFO,COU,SELEK,ACCURACY
  INTEGER J,I,W
  INTEGER ISEED(4),ITR
  DOUBLE PRECISION MIN_ERR,TMP_ERR,BK(K,K),VK(N,K),UK(M,K),VPLUS(N),VPLUS_ORIG(N),CONST
  DOUBLE PRECISION ZERO,MINUSONE,TEN
  PARAMETER (ZERO = 0.0D+0,MINUSONE = -1.0D+0,TEN = 10.0D+0)
  DOUBLE PRECISION ERR,DNRM2
  EXTERNAL ERR,DNRM2

  ISEED( 1 ) = 12
  ISEED( 2 ) = 14
  ISEED( 3 ) = 13
  ISEED( 4 ) = 733
  print *, "MODE = ",MODE
  WRITE (*,*) "M N K L"
  WRITE (*,*) M,N,K,L

  CONST = TEN**(ZERO - accuracy)

  DO W = 1,1
     CALL DLARNV(1,ISEED,N,VPLUS_ORIG)
     VPLUS_ORIG = VPLUS_ORIG / DNRM2(N,VPLUS_ORIG,1)

     DO SELEK = 1, 7 , 1
        STARTT = omp_get_wtime()
        WRITE(*,*) 
        BK = ZERO
        UK = ZERO
        VK = ZERO
        MIN_ERR = MINUSONE
        COU = 0
        ITR = 1
        J = 0
        VPLUS = VPLUS_ORIG
        DO WHILE (COU < 100)
           CALL RESGKL(J,mode,IAP,JA,A,M,N,K,L,BK,VK,UK,VPLUS,INFO,SELEK,WORK,LWORK)
           TMP_ERR = MAXVAL(ABS(BK(1:L,L+1)))

           ITR = ITR+1
           TMPT = omp_get_wtime()
           IF ( (TMP_ERR < MIN_ERR) .OR. (MIN_ERR .EQ. MINUSONE) ) THEN
              MIN_ERR = TMP_ERR
              MINT = TMPT
              COU = -1
           END IF

           WRITE(*,*) SELEK,ITR,(TMPT-STARTT),TMP_ERR
           IF(TMP_ERR .LE. CONST) exit
           COU = COU + 1
        END DO

        TE0 = ERR(mode,IAP,JA,A,M,N,K,L,BK,VK,UK)
        write(*,*)
        WRITE(*,*) W,"TIME", SELEK,ITR,(TMPT-STARTT),TMP_ERR,"AVE",TE0
        
        DO I = 1,3
           WRITE (*,*) "SingularValue(1-3)",I, BK(I,I)
        END DO
     END DO
  END DO

END SUBROUTINE RESGKL_MAIN
