SUBROUTINE RESGKL_MAIN(start_row,ini,mtd,MODE,LS,ACCURACY,M,N,L,K,IAP,JA,A,WORK,LWORK)
  use omp_lib
  IMPLICIT NONE

  CHARACTER MODE,LS
  DOUBLE PRECISION A(*),WORK(*),BK(K,K),VK(N,K),UK(M,K),VPLUS(N),VPLUS_ORIG(N)
  DOUBLE PRECISION MIN_ERR,TMP_ERR,TAVE,STARTT,TMPT,MINT,CONST,ZERO,MINUSONE,TEN,ERRMAX,ERR,DNRM2,TMAX
  PARAMETER (ZERO = 0.0D+0,MINUSONE = -1.0D+0,TEN = 10.0D+0)
  INTEGER ISEED(4), IAP(*), JA(*),ini,mtd,start_row(*)
  INTEGER M,N,L,K,INFO,COU,SELEK,ACCURACY,I,W,ITR,LWORK

  ISEED( 1 ) = ini
  ISEED( 2 ) = 402+ini*13
  ISEED( 3 ) = 2+ini*3
  ISEED( 4 ) = 111+ini*7
  print *, "MODE = ",MODE
  WRITE (*,*) "M N K L"
  WRITE (*,*) M,N,K,L

  CONST = TEN**(ZERO - accuracy)
  if(min(M,N)<K) then
     print *,"! ERR L must be smaller than",min(M,N)/2
  end if
  DO W = 1,1
     CALL DLARNV(1,ISEED,N,VPLUS_ORIG)
     VPLUS_ORIG = VPLUS_ORIG / DNRM2(N,VPLUS_ORIG,1)

     !DO SELEK = 1, 1!, 2
        SELEK = mtd
        STARTT = omp_get_wtime()
        WRITE(*,*) 
        BK = ZERO
        UK = ZERO
        VK = ZERO
        MIN_ERR = MINUSONE
        COU = 0
        ITR = 0
        I = 0
        VPLUS = VPLUS_ORIG
        DO WHILE (COU < 5000)
           CALL RESGKL(start_row,I,mode,LS,IAP,JA,A,M,N,K,L,BK,VK,UK,VPLUS,INFO,SELEK,WORK,LWORK)
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

        TAVE = ERR(start_row,mode,IAP,JA,A,M,N,K,L,BK,VK,UK,WORK)
        TMAX = ERRMAX(start_row,mode,IAP,JA,A,M,N,K,L,BK,VK,UK,WORK)
        WRITE(*,*) W,"TIME", SELEK,ITR,(TMPT-STARTT),TMP_ERR,"AVE",TAVE,"MAX",TMAX
        
        DO I = 1,3
           WRITE (*,*) "SingularValue(1-3)",I, BK(I,I)
        END DO
     !END DO
  END DO

END SUBROUTINE RESGKL_MAIN
