      SUBROUTINE DOQDS1(UPLO,N0,A,B,SU,LDSU,SVT,LDSVT,WORK,WORK2,INFO)
      
      IMPLICIT NONE
      
      CHARACTER UPLO
      INTEGER N0, INFO, LDSU, LDSVT, ISUB, MAXITER, SIT
      DOUBLE PRECISION A(*), B(*), WORK(N0,*)
      DOUBLE PRECISION SU(LDSU, *), SVT(LDSVT, *), WORK2(*)
      INTEGER N, M, I, J, K, OLDM, OLDN, IDIR, M0
      INTEGER INDRV1, INDRV2, INDRV3, INDRV4, INDRV5, INDRV6
      DOUBLE PRECISION TMP1, TMP2, TMP3, TMP4
      DOUBLE PRECISION TAU, TAU2
      DOUBLE PRECISION SIGMA, SIGMA2, DESIG, T, DESIG0
      DOUBLE PRECISION C1, S1, C2, S2, SMIN, EPS, TOL
      
      DOUBLE PRECISION ONE, ZERO, HALF, TWO
      INTEGER MAXREC
      PARAMETER (ONE = 1.0D0, ZERO = 0.0D0, HALF = 0.5D0, TWO = 2.0D0)
      DOUBLE PRECISION HUNDRD 
      PARAMETER (HUNDRD = 100.0D0)
      PARAMETER (MAXREC = 2000000000)
*     
      EXTERNAL DLAMCH
      DOUBLE PRECISION DLAMCH
      EXTERNAL DFMA0
      DOUBLE PRECISION DFMA0
      LOGICAL LSAME
      EXTERNAL LSAME
*     
      INDRV1 = 0
      INDRV2 = INDRV1+2*N0
      INDRV3 = INDRV2+2*N0
      INDRV4 = INDRV3+N0
      INDRV5 = INDRV4+N0
      INDRV6 = INDRV5+N0
*     
      EPS = DLAMCH( 'Precision' )
      TOL = HUNDRD*EPS
*     
      IF ( LSAME( UPLO, 'U' ) ) THEN
         
         DO I = 1, N0-1
            IF (A(I) .LT. ZERO) THEN
               A(I) = -A(I)
               B(I) = -B(I)
               CALL DSCAL(N0,-ONE,SU(1,I),1)
            ENDIF
            IF (B(I) .LT. ZERO) THEN
               B(I) = -B(I)
               A(I+1) = -A(I+1)
               CALL DSCAL(N0,-ONE,SVT(I+1,1),LDSVT)
            ENDIF
         ENDDO
         IF (A(N0) .LT. ZERO) THEN
            A(N0) = -A(N0)
            CALL DSCAL(N0,-ONE,SU(1,N0),1)
         ENDIF
         
         OLDM = 1
         OLDN = N0

         IF ( A(1) .GE. A(N0) ) THEN

            IDIR = 1

            TMP1 = A(1)
            DO J = 1, N0-1
               CALL DLARTG2(TMP1,B(J),C1,S1,A(J))
               WORK2(INDRV3+J) = C1
               WORK2(INDRV4+J) = S1
               B(J) = S1*A(J+1)
               TMP1 = C1*A(J+1)
            ENDDO
            A(N0) = TMP1
*     
            DO J = 1, N0-1
               IF( ( WORK2(INDRV3+J).NE.ONE ) .OR.
     $              ( WORK2(INDRV4+J).NE.ZERO ) ) THEN
                  CALL DROT(N0,SVT(J,1),LDSVT,SVT(J+1,1),LDSVT,
     $                 WORK2(INDRV3+J),WORK2(INDRV4+J))
               ENDIF
            ENDDO

         ELSE

            IDIR = 2

            TMP1 = A(N0)
            DO J = N0-1,1,-1
               CALL DLARTG2(TMP1,B(J),C1,S1,A(J+1))
               WORK2(INDRV1+J+1) = C1
               WORK2(INDRV2+J+1) = S1
               B(J) = S1*A(J)
               TMP1 = C1*A(J)
            ENDDO
            A(1) = TMP1

            DO J = N0-1, 1, -1
               IF( ( WORK2(INDRV1+J+1).NE.ONE ) .OR.
     $              ( WORK2(INDRV2+J+1).NE.ZERO ) ) THEN
                  CALL DROT(N0,SU(1,J+1),1,SU(1,J),1,
     $                 WORK2(INDRV1+J+1),WORK2(INDRV2+J+1))
               ENDIF
            ENDDO
            
         ENDIF
*     
      ELSE
*     
         DO I = 1, N0-1
            IF (A(I) .LT. ZERO) THEN
               A(I) = -A(I)
               B(I) = -B(I)
               CALL DSCAL(N0,-ONE,SVT(I,1),LDSVT)
            ENDIF
            IF (B(I) .LT. ZERO) THEN
               B(I) = -B(I)
               A(I+1) = -A(I+1)
               CALL DSCAL(N0,-ONE,SU(1,I+1),1)
            ENDIF
         ENDDO
         IF (A(N0) .LT. ZERO) THEN
            A(N0) = -A(N0)
            CALL DSCAL(N0,-ONE,SVT(N0,1),LDSVT)
         ENDIF
*     
         OLDM = -1
         OLDN = -1
*     
      ENDIF
*     
      DO I = 1, N0
         DO J = 1, N0
            WORK(J,I)=ZERO
         ENDDO
      ENDDO
*     
      M = 1
      N = N0
*     
      DO J = M, N-3
         IF (B(J) .LE. ZERO) THEN
            B(J) = -ZERO
            WORK2(INDRV6+J) = -ZERO
            M = J+1
         ENDIF
      ENDDO
*     
      B(N) = ZERO
      WORK2(INDRV6+N) = ZERO
*     
 3000 SIGMA = -B(N)
      SIGMA2 = TOL*SIGMA
      DESIG = -WORK2(INDRV6+N)

      MAXITER = 30*(N-M+1)
      DO I = 1, MAXITER
*     
 15      IF (N-M+1 .EQ. 1) THEN
            
            CALL DLARTG2(A(N),SIGMA,C1,S1,TMP1)
            
            IF( ( C1.NE.ONE ) .OR. ( S1.NE.ZERO ) ) THEN
               CALL DROT(N0,SU(1,N),1,WORK(1,N),1,C1,S1)
            ENDIF
            
            CALL DLARTG7(SIGMA,DESIG,A(N),A(N),DESIG0)
            
            GO TO 700
            
         ENDIF
         
         IF (N-M+1 .EQ. 2) THEN
            
            CALL DLASV2(A(N-1), B(N-1), A(N), A(N), A(N-1), S1, C1, 
     $           S2, C2)
            
            IF( ( C2.NE.ONE ) .OR. ( S2.NE.ZERO ) ) THEN
               CALL DROT(N0,SVT(N-1,1),LDSVT,SVT(N,1),LDSVT,C2,S2)
            ENDIF
            
            IF( ( C1.NE.ONE ) .OR. ( S1.NE.ZERO ) ) THEN
               CALL DROT(N0,SU(1,N-1),1,SU(1,N),1,C1,S1)
            ENDIF
            
            IF( ( C2.NE.ONE ) .OR. ( S2.NE.ZERO ) ) THEN
               CALL DROT(N0,WORK(1,N-1),1,WORK(1,N),1,C2,S2)
            ENDIF
            
            CALL DLARTG2(A(N-1),SIGMA,C1,S1,TMP1)
            
            IF( ( C1.NE.ONE ) .OR. ( S1.NE.ZERO ) ) THEN
               CALL DROT(N0,SU(1,N-1),1,WORK(1,N-1),1,C1,S1)
            ENDIF
            
            CALL DLARTG7(SIGMA,DESIG,A(N-1),A(N-1),DESIG0)
            
            CALL DLARTG2(A(N),SIGMA,C1,S1,TMP1)
            
            IF( ( C1.NE.ONE ) .OR. ( S1.NE.ZERO ) ) THEN
               CALL DROT(N0,SU(1,N),1,WORK(1,N),1,C1,S1)
            ENDIF
            
            CALL DLARTG7(SIGMA,DESIG,A(N),A(N),DESIG0)
            
            GO TO 700
         ENDIF
         
         IF (B(N-1) .LE. SIGMA2) THEN
            
            CALL DLARTG2(A(N),SIGMA,C1,S1,TMP1)
            
            IF( ( C1.NE.ONE ) .OR. ( S1.NE.ZERO ) ) THEN
               CALL DROT(N0,SU(1,N),1,WORK(1,N),1,C1,S1)
            ENDIF
            
            CALL DLARTG7(SIGMA,DESIG,A(N),A(N),DESIG0)
            
            N = N - 1
            GO TO 15
            
         ENDIF
*     
         IF (B(N-2) .LE. SIGMA2) THEN
            
            CALL DLASV2(A(N-1), B(N-1), A(N), A(N), A(N-1), S1, 
     $           C1, S2, C2)
            
            IF( ( C2.NE.ONE ) .OR. ( S2.NE.ZERO ) ) THEN
               CALL DROT(N0,SVT(N-1,1),LDSVT,SVT(N,1),LDSVT,C2,S2)
            ENDIF
            
            IF( ( C1.NE.ONE ) .OR. ( S1.NE.ZERO ) ) THEN
               CALL DROT(N0,SU(1,N-1),1,SU(1,N),1,C1,S1)
            ENDIF
            
            IF( ( C2.NE.ONE ) .OR. ( S2.NE.ZERO ) ) THEN
               CALL DROT(N0,WORK(1,N-1),1,WORK(1,N),1,C2,S2)
            ENDIF
            
            CALL DLARTG2(A(N-1),SIGMA,C1,S1,TMP1)
            
            IF( ( C1.NE.ONE ) .OR. ( S1.NE.ZERO ) ) THEN
               CALL DROT(N0,SU(1,N-1),1,WORK(1,N-1),1,C1,S1)
            ENDIF
            
            CALL DLARTG7(SIGMA,DESIG,A(N-1),A(N-1),DESIG0)
            
            CALL DLARTG2(A(N),SIGMA,C1,S1,TMP1)
            
            IF( ( C1.NE.ONE ) .OR. ( S1.NE.ZERO ) ) THEN
               CALL DROT(N0,SU(1,N),1,WORK(1,N),1,C1,S1)
            ENDIF
            
            CALL DLARTG7(SIGMA,DESIG,A(N),A(N),DESIG0)
            
            N = N - 2
            GO TO 15
            
         ENDIF
*     
         IF (M .GT. OLDN .OR. N .LT. OLDM) THEN
            IF ( A(M) .GE. A(N) ) THEN
               IDIR = 1
            ELSE
               IDIR = 2
            ENDIF
         ENDIF
*     
         IF ( IDIR .EQ. 1) THEN
*     
            OLDM = M
            OLDN = N

            CALL DLAS2(A(N-1), B(N-1), A(N), TAU, TMP3)
            TAU = MIN(TAU,A(N))
            IF (TAU .EQ. ZERO) GO TO 350

            TMP2 = MIN(A(N),A(N-1))
            TMP3 = MAX(A(N),A(N-1))
            IF (TMP3 .GE. TWO*TMP2) THEN
               SIT = 1
            ELSE
               SIT = 0
            ENDIF

            CALL DLARTG7(SIGMA,DESIG,TAU,T,DESIG0)
            IF (T .LE. SIGMA .AND. SIT .EQ. 1) GO TO 350

            TAU2 = MINVAL(A(M:N-1))
            IF (TAU .GE. TAU2) THEN
               IF (TAU2 .EQ. ZERO) GO TO 350
               CALL DLARTG7(SIGMA,DESIG,TAU2,T,DESIG0)
               IF (T .LE. SIGMA .AND. SIT .EQ. 1) GO TO 350
               GO TO 160
            ENDIF
*     
            TMP4 = A(M)-TAU
            IF (TMP4 .LE. ZERO) THEN
               GO TO 160
            ELSE
               TMP3 = SQRT(TMP4)*SQRT(A(M)+TAU)
            ENDIF
            DO J = M, N-2
               CALL DLARTG4(TMP3,B(J),C1)
               TMP4 = DFMA0(C1,A(J+1),-TAU)
               IF (TMP4 .LE. ZERO) THEN
                  GO TO 160
               ELSE
                  TMP3 = SQRT(TMP4)*SQRT(DFMA0(C1,A(J+1),TAU))
               ENDIF
            ENDDO
            CALL DLARTG4(TMP3,B(N-1),C1)
            TMP4 = DFMA0(C1,A(N),-TAU)
            IF (TMP4 .LT. ZERO) THEN
               DO K=0,MAXREC
                  TAU2 = MAX(ONE-DBLE(K)*EPS,ZERO)*(C1*A(N))
                  IF (TAU2 .LT. TAU) EXIT
               ENDDO
               TAU = TAU2
            ENDIF
            IF (TAU .GT. ZERO) GO TO 125
*     
 160        WORK2(INDRV5+N) = ONE/A(N)
            TMP3 = WORK2(INDRV5+N)
            DO J = N-1,M,-1
               WORK2(INDRV1+J) = B(J)/A(J)
               WORK2(INDRV5+J) = ONE/A(J)+
     $              WORK2(INDRV1+J)*WORK2(INDRV5+J+1)
               TMP3 = MAX(TMP3,WORK2(INDRV5+J))
            ENDDO
            WORK2(INDRV5+M) = (WORK2(INDRV5+M)/TMP3)/A(M)
            TMP1 = WORK2(INDRV5+M)
            DO J = M+1,N
               WORK2(INDRV2+J) = B(J-1)/A(J)
               WORK2(INDRV5+J) = (WORK2(INDRV5+J)/TMP3)/A(J)+
     $              WORK2(INDRV2+J)*WORK2(INDRV5+J-1)
               TMP1 = MAX(TMP1,WORK2(INDRV5+J))
            ENDDO
            TAU = ZERO
            IF (TMP3 .GT. ZERO .AND. TMP1 .GT. ZERO) THEN
               TAU = MAX(TAU,ONE/SQRT(TMP1)/SQRT(TMP3))
            ENDIF
            
            WORK2(INDRV5+N) = WORK2(INDRV5+N)/TMP1
            WORK2(INDRV6+N) = WORK2(INDRV5+N)/A(N)
            TMP3 = WORK2(INDRV6+N)
            DO J = N-1,M,-1
               WORK2(INDRV5+J) = WORK2(INDRV5+J)/TMP1
               WORK2(INDRV6+J) = WORK2(INDRV5+J)/A(J)+
     $              WORK2(INDRV1+J)*WORK2(INDRV6+J+1)
               TMP3 = MAX(TMP3,WORK2(INDRV6+J))
            ENDDO
            TMP2 = (WORK2(INDRV6+M)/TMP3)/A(M)
            TMP1 = TMP2/WORK2(INDRV5+M)
            DO J = M+1,N
               TMP2 = (WORK2(INDRV6+J)/TMP3)/A(J)+WORK2(INDRV2+J)*TMP2
               TMP1 = MAX(TMP1,TMP2/WORK2(INDRV5+J))
            ENDDO
            IF (TMP3 .GT. ZERO .AND. TMP1 .GT. ZERO) THEN
               TAU = MAX(TAU,ONE/SQRT(TMP1)/SQRT(TMP3))
            ENDIF

            IF (TAU .GT. ZERO) GO TO 125
       
            TAU=A(N)-HALF*B(N-1)
            IF (TAU .LE. ZERO) GO TO 350
            DO J=N-1,M+1,-1
               TMP2=A(J)-HALF*(B(J)+B(J-1))
               IF (TMP2 .LE. ZERO) GO TO 350
               TAU=MIN(TAU,TMP2)
            ENDDO
            TMP2=A(M)-HALF*B(M)
            IF (TMP2 .LE. ZERO) GO TO 350
            TAU=MIN(TAU,TMP2)
*     
 125        IF (TAU .EQ. ZERO) GO TO 350
            CALL DLARTG7(SIGMA,DESIG,TAU,T,DESIG0)
            IF (T .LE. SIGMA .AND. SIT .EQ. 1) GO TO 350
*     
            IF (SIGMA .EQ. ZERO) THEN
*     
               TMP4 = A(M)-TAU
               IF (TMP4 .LE. ZERO) THEN
                  DO K=0,MAXREC
                     TAU2 = MAX(ONE-DBLE(K)*EPS,ZERO)*A(M)
                     IF (TAU2 .LT. TAU) EXIT
                  ENDDO
                  TAU = TAU2
                  GO TO 125
               ELSE
                  TMP3 = SQRT(TMP4)*SQRT(A(M)+TAU)
               ENDIF
               
               CALL DLARTG3(TMP3,-T,C1,S1)
               WORK2(INDRV1+2*M-1) = C1
               WORK2(INDRV2+2*M-1) = S1
               DO J = M, N-2
                  CALL DLARTG2(TMP3,B(J),C1,S1,WORK2(INDRV5+J))
                  WORK2(INDRV1+2*J) = C1
                  WORK2(INDRV2+2*J) = S1
                  WORK2(INDRV6+J) = S1*A(J+1)
                  
                  TMP4 = DFMA0(C1,A(J+1),-TAU)
                  IF (TMP4 .LE. ZERO) THEN
                     DO K=0,MAXREC
                        TAU2 = MAX(ONE-DBLE(K)*EPS,ZERO)*(C1*A(J+1))
                        IF (TAU2 .LT. TAU) EXIT
                     ENDDO
                     TAU = MAX(HALF*TAU,TAU2)
                     GO TO 125
                  ELSE
                     TMP3 = SQRT(TMP4)*SQRT(DFMA0(C1,A(J+1),TAU))
                  ENDIF
                  
                  CALL DLARTG3(TMP3,-T,C1,S1)
                  WORK2(INDRV1+2*J+1) = C1
                  WORK2(INDRV2+2*J+1) = S1
               ENDDO
               CALL DLARTG2(TMP3,B(N-1),C1,S1,WORK2(INDRV5+N-1))
               WORK2(INDRV1+2*N-2) = C1
               WORK2(INDRV2+2*N-2) = S1
               WORK2(INDRV6+N-1) = S1*A(N)
               
               TMP4 = DFMA0(C1,A(N),-TAU)
               IF(TMP4 .LT. ZERO .AND. TAU+TMP4 .EQ. TAU) TMP4=ZERO
               IF(TMP4 .LT. ZERO) THEN
                  DO K=0,MAXREC
                     TAU2 = MAX(ONE-DBLE(K)*EPS,ZERO)*(C1*A(N))
                     IF (TAU2 .LT. TAU) EXIT
                  ENDDO
                  TAU = MAX(HALF*TAU,TAU2)
                  GO TO 125
               ELSE
                  TMP3 = SQRT(TMP4)*SQRT(DFMA0(C1,A(N),TAU))
               ENDIF
               
               CALL DLARTG3(TMP3,-T,C1,S1)
               WORK2(INDRV1+2*N-1) = C1
               WORK2(INDRV2+2*N-1) = S1
               WORK2(INDRV5+N) = TMP3
*     
            ELSE
*     
               TMP4 = A(M)-TAU
               IF (TMP4 .LE. ZERO) THEN
                  DO K=0,MAXREC
                     TAU2 = MAX(ONE-DBLE(K)*EPS,ZERO)*A(M)
                     IF (TAU2 .LT. TAU) EXIT
                  ENDDO
                  TAU = TAU2
                  GO TO 125
               ELSE
                  TMP3 = SQRT(TMP4)*SQRT(A(M)+TAU)
               ENDIF
               
               CALL DLARTG6(A(M),TMP3,SIGMA,T,C1,S1)
               WORK2(INDRV1+2*M-1) = C1
               WORK2(INDRV2+2*M-1) = S1
               DO J = M, N-2
                  CALL DLARTG2(TMP3,B(J),C1,S1,WORK2(INDRV5+J))
                  WORK2(INDRV1+2*J) = C1
                  WORK2(INDRV2+2*J) = S1
                  WORK2(INDRV6+J) = S1*A(J+1)
                  
                  TMP4 = DFMA0(C1,A(J+1),-TAU)
                  IF (TMP4 .LE. ZERO) THEN
                     DO K=0,MAXREC
                        TAU2 = MAX(ONE-DBLE(K)*EPS,ZERO)*(C1*A(J+1))
                        IF (TAU2 .LT. TAU) EXIT
                     ENDDO
                     TAU = MAX(HALF*TAU,TAU2)
                     GO TO 125
                  ELSE
                     TMP3 = SQRT(TMP4)*SQRT(DFMA0(C1,A(J+1),TAU))
                  ENDIF
                  
                  CALL DLARTG6(C1*A(J+1),TMP3,SIGMA,T,C1,S1)
                  WORK2(INDRV1+2*J+1) = C1
                  WORK2(INDRV2+2*J+1) = S1
               ENDDO
               CALL DLARTG2(TMP3,B(N-1),C1,S1,WORK2(INDRV5+N-1))
               WORK2(INDRV1+2*N-2) = C1
               WORK2(INDRV2+2*N-2) = S1
               WORK2(INDRV6+N-1) = S1*A(N)
               
               TMP4 = DFMA0(C1,A(N),-TAU)
               IF(TMP4 .LT. ZERO .AND. TAU+TMP4 .EQ. TAU) TMP4=ZERO
               IF(TMP4 .LT. ZERO) THEN
                  DO K=0,MAXREC
                     TAU2 = MAX(ONE-DBLE(K)*EPS,ZERO)*(C1*A(N))
                     IF (TAU2 .LT. TAU) EXIT
                  ENDDO
                  TAU = MAX(HALF*TAU,TAU2)
                  GO TO 125
               ELSE
                  TMP3 = SQRT(TMP4)*SQRT(DFMA0(C1,A(N),TAU))
               ENDIF
               
               IF (TMP3 .EQ. ZERO) THEN
                  CALL DLARTG3(SIGMA,-C1*A(N),C1,S1)
               ELSE
                  CALL DLARTG6(C1*A(N),TMP3,SIGMA,T,C1,S1)
               ENDIF
               WORK2(INDRV1+2*N-1) = C1
               WORK2(INDRV2+2*N-1) = S1
               WORK2(INDRV5+N) = TMP3
*     
            ENDIF
*     
            SIGMA = T
            SIGMA2 = TOL*SIGMA
            DESIG = DESIG0
*     
            TMP1 = WORK2(INDRV5+M)
            DO J = M, N-1
               CALL DLARTG2(TMP1,WORK2(INDRV6+J),C1,S1,A(J))
               WORK2(INDRV3+J) = C1
               WORK2(INDRV4+J) = S1
               B(J) = S1*WORK2(INDRV5+J+1)
               TMP1 = C1*WORK2(INDRV5+J+1)
            ENDDO
            A(N) = TMP1
*     
            IF( ( WORK2(INDRV1+2*M-1).NE.ONE ) .OR. 
     $           ( WORK2(INDRV2+2*M-1).NE.ZERO ) ) THEN
               CALL DROT(N0,SU(1,M),1,WORK(1,M),1,
     $              WORK2(INDRV1+2*M-1),WORK2(INDRV2+2*M-1))
            ENDIF
            DO J = M, N-1
               IF( ( WORK2(INDRV1+2*J).NE.ONE ) .OR. 
     $              ( WORK2(INDRV2+2*J).NE.ZERO ) ) THEN
                  CALL DROT(N0,SU(1,J),1,SU(1,J+1),1,WORK2(INDRV1+2*J),
     $                 WORK2(INDRV2+2*J))
               ENDIF
               IF( ( WORK2(INDRV1+2*J+1).NE.ONE ) .OR. 
     $              ( WORK2(INDRV2+2*J+1).NE.ZERO ) ) 
     $              THEN
                  CALL DROT(N0,SU(1,J+1),1,WORK(1,J+1),1,
     $                 WORK2(INDRV1+2*J+1),WORK2(INDRV2+2*J+1))
               ENDIF
               IF( ( WORK2(INDRV3+J).NE.ONE ) .OR. 
     $              ( WORK2(INDRV4+J).NE.ZERO ) ) THEN
                  CALL DROT(N0,WORK(1,J),1,WORK(1,J+1),1,
     $                 WORK2(INDRV3+J),WORK2(INDRV4+J))
               ENDIF
            ENDDO
*     
            DO J = M, N-1
               IF( ( WORK2(INDRV3+J).NE.ONE ) .OR. 
     $              ( WORK2(INDRV4+J).NE.ZERO ) ) THEN
                  CALL DROT(N0,SVT(J,1),LDSVT,SVT(J+1,1),LDSVT,
     $                 WORK2(INDRV3+J),WORK2(INDRV4+J))
               ENDIF
            ENDDO
*     
            GO TO 400
*     
 350        TMP1 = A(M)
            DO J = M, N-1
               CALL DLARTG2(TMP1,B(J),C1,S1,WORK2(INDRV5+J))
               WORK2(INDRV1+J) = C1
               WORK2(INDRV2+J) = S1
               WORK2(INDRV6+J) = S1*A(J+1)
               TMP1 = C1*A(J+1)
            ENDDO
            WORK2(INDRV5+N) = TMP1
*     
            DO J = M, N-1
               IF( ( WORK2(INDRV1+J).NE.ONE ) .OR. 
     $              ( WORK2(INDRV2+J).NE.ZERO ) ) THEN
                  CALL DROT(N0,SU(1,J),1,SU(1,J+1),1,WORK2(INDRV1+J),
     $                 WORK2(INDRV2+J))
               ENDIF
            ENDDO
*
            GO TO 395
*     
 390        IF (N-M+1 .EQ. 1) THEN
            
               CALL DLARTG2(WORK2(INDRV5+N),SIGMA,C1,S1,TMP1)
               
               IF( ( C1.NE.ONE ) .OR. ( S1.NE.ZERO ) ) THEN
                  CALL DROT(N0,SU(1,N),1,WORK(1,N),1,C1,S1)
               ENDIF
               
               CALL DLARTG7(SIGMA,DESIG,WORK2(INDRV5+N),A(N),DESIG0)
               
               GO TO 700
               
            ENDIF
            
            IF (N-M+1 .EQ. 2) THEN
               
               CALL DLASV2(WORK2(INDRV5+N-1), WORK2(INDRV6+N-1),
     $              WORK2(INDRV5+N), WORK2(INDRV5+N), WORK2(INDRV5+N-1),
     $              S1, C1, S2, C2)
               
               IF( ( C1.NE.ONE ) .OR. ( S1.NE.ZERO ) ) THEN
                  CALL DROT(N0,SVT(N-1,1),LDSVT,SVT(N,1),LDSVT,C1,S1)
               ENDIF
               
               IF( ( C2.NE.ONE ) .OR. ( S2.NE.ZERO ) ) THEN
                  CALL DROT(N0,SU(1,N-1),1,SU(1,N),1,C2,S2)
               ENDIF
               
               IF( ( C1.NE.ONE ) .OR. ( S1.NE.ZERO ) ) THEN
                  CALL DROT(N0,WORK(1,N-1),1,WORK(1,N),1,C1,S1)
               ENDIF
               
               CALL DLARTG2(WORK2(INDRV5+N-1),SIGMA,C1,S1,TMP1)
               
               IF( ( C1.NE.ONE ) .OR. ( S1.NE.ZERO ) ) THEN
                  CALL DROT(N0,SU(1,N-1),1,WORK(1,N-1),1,C1,S1)
               ENDIF
               
               CALL DLARTG7(SIGMA,DESIG,WORK2(INDRV5+N-1),A(N-1),DESIG0)
               
               CALL DLARTG2(WORK2(INDRV5+N),SIGMA,C1,S1,TMP1)
               
               IF( ( C1.NE.ONE ) .OR. ( S1.NE.ZERO ) ) THEN
                  CALL DROT(N0,SU(1,N),1,WORK(1,N),1,C1,S1)
               ENDIF
               
               CALL DLARTG7(SIGMA,DESIG,WORK2(INDRV5+N),A(N),DESIG0)
               
               GO TO 700
            ENDIF
            
 395        IF (WORK2(INDRV6+N-1) .LE. SIGMA2) THEN
               
               CALL DLARTG2(WORK2(INDRV5+N),SIGMA,C1,S1,TMP1)
               
               IF( ( C1.NE.ONE ) .OR. ( S1.NE.ZERO ) ) THEN
                  CALL DROT(N0,SU(1,N),1,WORK(1,N),1,C1,S1)
               ENDIF
               
               CALL DLARTG7(SIGMA,DESIG,WORK2(INDRV5+N),A(N),DESIG0)
               
               N = N - 1
               GO TO 390
               
            ENDIF
*     
            IF (WORK2(INDRV6+N-2) .LE. SIGMA2) THEN
               
               CALL DLASV2(WORK2(INDRV5+N-1), WORK2(INDRV6+N-1),
     $              WORK2(INDRV5+N), WORK2(INDRV5+N), WORK2(INDRV5+N-1),
     $              S1, C1, S2, C2)
               
               IF( ( C1.NE.ONE ) .OR. ( S1.NE.ZERO ) ) THEN
                  CALL DROT(N0,SVT(N-1,1),LDSVT,SVT(N,1),LDSVT,C1,S1)
               ENDIF
               
               IF( ( C2.NE.ONE ) .OR. ( S2.NE.ZERO ) ) THEN
                  CALL DROT(N0,SU(1,N-1),1,SU(1,N),1,C2,S2)
               ENDIF
               
               IF( ( C1.NE.ONE ) .OR. ( S1.NE.ZERO ) ) THEN
                  CALL DROT(N0,WORK(1,N-1),1,WORK(1,N),1,C1,S1)
               ENDIF
               
               CALL DLARTG2(WORK2(INDRV5+N-1),SIGMA,C1,S1,TMP1)
               
               IF( ( C1.NE.ONE ) .OR. ( S1.NE.ZERO ) ) THEN
                  CALL DROT(N0,SU(1,N-1),1,WORK(1,N-1),1,C1,S1)
               ENDIF
               
               CALL DLARTG7(SIGMA,DESIG,WORK2(INDRV5+N-1),A(N-1),DESIG0)
               
               CALL DLARTG2(WORK2(INDRV5+N),SIGMA,C1,S1,TMP1)
               
               IF( ( C1.NE.ONE ) .OR. ( S1.NE.ZERO ) ) THEN
                  CALL DROT(N0,SU(1,N),1,WORK(1,N),1,C1,S1)
               ENDIF
               
               CALL DLARTG7(SIGMA,DESIG,WORK2(INDRV5+N),A(N),DESIG0)
               
               N = N - 2
               GO TO 390
               
            ENDIF
*
            TMP1 = WORK2(INDRV5+M)
            DO J = M, N-1
               CALL DLARTG2(TMP1,WORK2(INDRV6+J),C1,S1,A(J))
               WORK2(INDRV3+J) = C1
               WORK2(INDRV4+J) = S1
               B(J) = S1*WORK2(INDRV5+J+1)
               TMP1 = C1*WORK2(INDRV5+J+1)
            ENDDO
            A(N) = TMP1
*     
            DO J = M, N-1
               IF( ( WORK2(INDRV3+J).NE.ONE ) .OR. 
     $              ( WORK2(INDRV4+J).NE.ZERO ) ) THEN
                  CALL DROT(N0,WORK(1,J),1,WORK(1,J+1),1,
     $                 WORK2(INDRV3+J),WORK2(INDRV4+J))
               ENDIF
            ENDDO
*     
            DO J = M, N-1
               IF( ( WORK2(INDRV3+J).NE.ONE ) .OR. 
     $              ( WORK2(INDRV4+J).NE.ZERO ) ) THEN
                  CALL DROT(N0,SVT(J,1),LDSVT,SVT(J+1,1),LDSVT,
     $                 WORK2(INDRV3+J),WORK2(INDRV4+J))
               ENDIF
            ENDDO
*     
 400        TMP1 = A(M)
            DO J = M, N-3
               IF (B(J) .LE. SIGMA2 .OR. TMP1 .EQ. B(J)+TMP1) THEN
                  B(J) = -SIGMA
                  WORK2(INDRV6+J) = -DESIG
                  M = J+1
                  TMP1 = A(J+1)
               ELSE
                  TMP1 = A(J+1)*(TMP1/(TMP1+B(J)))
               ENDIF
            ENDDO

            IF (B(N-2) .LE. SIGMA2 .OR. TMP1 .EQ. B(N-2)+TMP1) THEN
               B(N-2) = ZERO
               TMP1 = A(N-1)
            ELSE
               TMP1 = A(N-1)*(TMP1/(TMP1+B(N-2)))
            ENDIF
 
            IF (B(N-1) .LE. SIGMA2 .OR. TMP1 .EQ. B(N-1)+TMP1) THEN
               B(N-1) = ZERO
            ENDIF
*     
         ELSE
*     
            OLDM = M
            OLDN = N

            TMP1 = A(N)
            DO J = N-1,M,-1
               CALL DLARTG2(TMP1,B(J),C1,S1,WORK2(INDRV5+J+1))
               WORK2(INDRV3+J+1) = C1
               WORK2(INDRV4+J+1) = S1
               WORK2(INDRV6+J) = S1*A(J)
               TMP1 = C1*A(J)
            ENDDO
            WORK2(INDRV5+M) = TMP1
*     
            CALL DLAS2(WORK2(INDRV5+M), WORK2(INDRV6+M), 
     $           WORK2(INDRV5+M+1), TAU, TMP3)
            TAU = MIN(TAU,WORK2(INDRV5+M))
            IF (TAU .EQ. ZERO) GO TO 1350

            TMP2 = MIN(WORK2(INDRV5+M),WORK2(INDRV5+M+1))
            TMP3 = MAX(WORK2(INDRV5+M),WORK2(INDRV5+M+1))
            IF (TMP3 .GE. TWO*TMP2) THEN
               SIT = 1
            ELSE
               SIT = 0
            ENDIF

            CALL DLARTG7(SIGMA,DESIG,TAU,T,DESIG0)
            IF (T .LE. SIGMA .AND. SIT .EQ. 1) GO TO 1350

            TAU2 = MINVAL(WORK2(INDRV5+M+1:INDRV5+N))
            IF (TAU .GE. TAU2) THEN
               IF (TAU2 .EQ. ZERO) GO TO 1350
               CALL DLARTG7(SIGMA,DESIG,TAU2,T,DESIG0)
               IF (T .LE. SIGMA .AND. SIT .EQ. 1) GO TO 1350
               GO TO 1160
            ENDIF
*     
            TMP4 = WORK2(INDRV5+N)-TAU
            IF (TMP4 .LE. ZERO) THEN
               GO TO 1160
            ELSE
               TMP3 = SQRT(TMP4)*SQRT(WORK2(INDRV5+N)+TAU)
            ENDIF
            DO J = N-1,M+1,-1
               CALL DLARTG4(TMP3,WORK2(INDRV6+J),C1)
               TMP4 = DFMA0(C1,WORK2(INDRV5+J),-TAU)
               IF (TMP4 .LE. ZERO) THEN
                  GO TO 1160
               ELSE
                  TMP3 = SQRT(TMP4)*SQRT(DFMA0(C1,WORK2(INDRV5+J),TAU))
               ENDIF
            ENDDO
            CALL DLARTG4(TMP3,WORK2(INDRV6+M),C1)
            TMP4 = DFMA0(C1,WORK2(INDRV5+M),-TAU)
            IF (TMP4 .LT. ZERO) THEN
               DO K=0,MAXREC
                  TAU2 = MAX(ONE-DBLE(K)*EPS,ZERO)*(C1*WORK2(INDRV5+M))
                  IF (TAU2 .LT. TAU) EXIT
               ENDDO
               TAU = TAU2
            ENDIF
            IF (TAU .GT. ZERO) GO TO 1125
*     
 1160       A(M) = ONE/WORK2(INDRV5+M)
            TMP3 = A(M)
            DO J = M+1,N
               WORK2(INDRV1+J) = WORK2(INDRV6+J-1)/WORK2(INDRV5+J)
               A(J) = ONE/WORK2(INDRV5+J)+WORK2(INDRV1+J)*A(J-1)
               TMP3 = MAX(TMP3,A(J))
            ENDDO
            A(N) = (A(N)/TMP3)/WORK2(INDRV5+N)
            TMP1 = A(N)
            DO J = N-1,M,-1
               WORK2(INDRV2+J) = WORK2(INDRV6+J)/WORK2(INDRV5+J)
               A(J) = (A(J)/TMP3)/WORK2(INDRV5+J)+WORK2(INDRV2+J)*A(J+1)
               TMP1 = MAX(TMP1,A(J))
            ENDDO
            TAU = ZERO
            IF (TMP3 .GT. ZERO .AND. TMP1 .GT. ZERO) THEN
               TAU = MAX(TAU,ONE/SQRT(TMP1)/SQRT(TMP3))
            ENDIF
            
            A(M) = A(M)/TMP1
            B(M) = A(M)/WORK2(INDRV5+M)
            TMP3 = B(M)
            DO J = M+1, N
               A(J) = A(J)/TMP1
               B(J) = A(J)/WORK2(INDRV5+J)+WORK2(INDRV1+J)*B(J-1)
               TMP3 = MAX(TMP3,B(J))
            ENDDO
            TMP2 = (B(N)/TMP3)/WORK2(INDRV5+N)
            TMP1 = TMP2/A(N)
            DO J = N-1,M,-1
               TMP2 = (B(J)/TMP3)/WORK2(INDRV5+J)+WORK2(INDRV2+J)*TMP2
               TMP1 = MAX(TMP1,TMP2/A(J))
            ENDDO
            IF (TMP3 .GT. ZERO .AND. TMP1 .GT. ZERO) THEN
               TAU = MAX(TAU,ONE/SQRT(TMP1)/SQRT(TMP3))
            ENDIF

            IF (TAU .GT. ZERO) GO TO 1125
            
            TAU=WORK2(INDRV5+M)-HALF*WORK2(INDRV6+M)
            IF (TAU .LE. ZERO) GO TO 1350
            DO J=M+1,N-1
               TMP2=WORK2(INDRV5+J)
     $              -HALF*(WORK2(INDRV6+J)+WORK2(INDRV6+J-1))
               IF (TMP2 .LE. ZERO) GO TO 1350
               TAU=MIN(TAU,TMP2)
            ENDDO
            TMP2=WORK2(INDRV5+N)-HALF*WORK2(INDRV6+N-1)
            IF (TMP2 .LE. ZERO) GO TO 1350
            TAU=MIN(TAU,TMP2)
*     
 1125       IF (TAU .EQ. ZERO) GO TO 1350
            CALL DLARTG7(SIGMA,DESIG,TAU,T,DESIG0)
            IF (T .LE. SIGMA .AND. SIT .EQ. 1) GO TO 1350
*     
            IF (SIGMA .EQ. ZERO) THEN
*     
               TMP4 = WORK2(INDRV5+N)-TAU
               IF (TMP4 .LE. ZERO) THEN
                  DO K=0,MAXREC
                     TAU2 = MAX(ONE-DBLE(K)*EPS,ZERO)*WORK2(INDRV5+N)
                     IF (TAU2 .LT. TAU) EXIT
                  ENDDO
                  TAU = TAU2
                  GO TO 1125
               ELSE
                  TMP3 = SQRT(TMP4)*SQRT(WORK2(INDRV5+N)+TAU)
               ENDIF
               
               CALL DLARTG3(TMP3,-T,C1,S1)
               WORK2(INDRV1+2*N) = C1
               WORK2(INDRV2+2*N) = S1
               DO J = N-1, M+1, -1
                  CALL DLARTG2(TMP3,WORK2(INDRV6+J),C1,S1,A(J+1))
                  WORK2(INDRV1+2*J+1) = C1
                  WORK2(INDRV2+2*J+1) = S1
                  B(J) = S1*WORK2(INDRV5+J)
                  
                  TMP4 = DFMA0(C1,WORK2(INDRV5+J),-TAU)
                  IF (TMP4 .LE. ZERO) THEN
                     DO K=0,MAXREC
                        TAU2 = MAX(ONE-DBLE(K)*EPS,ZERO)*
     $                       (C1*WORK2(INDRV5+J))
                        IF (TAU2 .LT. TAU) EXIT
                     ENDDO
                     TAU = MAX(HALF*TAU,TAU2)
                     GO TO 1125
                  ELSE
                     TMP3 = SQRT(TMP4)*
     $                    SQRT(DFMA0(C1,WORK2(INDRV5+J),TAU))
                  ENDIF
                  
                  CALL DLARTG3(TMP3,-T,C1,S1)
                  WORK2(INDRV1+2*J) = C1
                  WORK2(INDRV2+2*J) = S1
               ENDDO
               CALL DLARTG2(TMP3,WORK2(INDRV6+M),C1,S1,A(M+1))
               WORK2(INDRV1+2*M+1) = C1
               WORK2(INDRV2+2*M+1) = S1
               B(M) = S1*WORK2(INDRV5+M)
               
               TMP4 = DFMA0(C1,WORK2(INDRV5+M),-TAU)
               IF(TMP4 .LT. ZERO .AND. TAU+TMP4 .EQ. TAU) TMP4=ZERO
               IF(TMP4 .LT. ZERO) THEN
                  DO K=0,MAXREC
                     TAU2 = MAX(ONE-DBLE(K)*EPS,ZERO)*
     $                    (C1*WORK2(INDRV5+M))
                     IF (TAU2 .LT. TAU) EXIT
                  ENDDO
                  TAU = MAX(HALF*TAU,TAU2)
                  GO TO 1125
               ELSE
                  TMP3 = SQRT(TMP4)*SQRT(DFMA0(C1,WORK2(INDRV5+M),TAU))
               ENDIF
               
               CALL DLARTG3(TMP3,-T,C1,S1)
               WORK2(INDRV1+2*M) = C1
               WORK2(INDRV2+2*M) = S1
               A(M) = TMP3
*     
            ELSE
*     
               TMP4 = WORK2(INDRV5+N)-TAU
               IF (TMP4 .LE. ZERO) THEN
                  DO K=0,MAXREC
                     TAU2 = MAX(ONE-DBLE(K)*EPS,ZERO)*WORK2(INDRV5+N)
                     IF (TAU2 .LT. TAU) EXIT
                  ENDDO
                  TAU = TAU2
                  GO TO 1125
               ELSE
                  TMP3 = SQRT(TMP4)*SQRT(WORK2(INDRV5+N)+TAU)
               ENDIF
               
               CALL DLARTG6(WORK2(INDRV5+N),TMP3,SIGMA,T,C1,S1)
               WORK2(INDRV1+2*N) = C1
               WORK2(INDRV2+2*N) = S1
               DO J = N-1, M+1, -1
                  CALL DLARTG2(TMP3,WORK2(INDRV6+J),C1,S1,A(J+1))
                  WORK2(INDRV1+2*J+1) = C1
                  WORK2(INDRV2+2*J+1) = S1
                  B(J) = S1*WORK2(INDRV5+J)
                  
                  TMP4 = DFMA0(C1,WORK2(INDRV5+J),-TAU)
                  IF (TMP4 .LE. ZERO) THEN
                     DO K=0,MAXREC
                        TAU2 = MAX(ONE-DBLE(K)*EPS,ZERO)*
     $                       (C1*WORK2(INDRV5+J))
                        IF (TAU2 .LT. TAU) EXIT
                     ENDDO
                     TAU = MAX(HALF*TAU,TAU2)
                     GO TO 1125
                  ELSE
                     TMP3 = SQRT(TMP4)*
     $                    SQRT(DFMA0(C1,WORK2(INDRV5+J),TAU))
                  ENDIF
                  
                  CALL DLARTG6(C1*WORK2(INDRV5+J),TMP3,SIGMA,T,C1,S1)
                  WORK2(INDRV1+2*J) = C1
                  WORK2(INDRV2+2*J) = S1
               ENDDO
               CALL DLARTG2(TMP3,WORK2(INDRV6+M),C1,S1,A(M+1))
               WORK2(INDRV1+2*M+1) = C1
               WORK2(INDRV2+2*M+1) = S1
               B(M) = S1*WORK2(INDRV5+M)
               
               TMP4 = DFMA0(C1,WORK2(INDRV5+M),-TAU)
               IF(TMP4 .LT. ZERO .AND. TAU+TMP4 .EQ. TAU) TMP4=ZERO
               IF(TMP4 .LT. ZERO) THEN
                  DO K=0,MAXREC
                     TAU2 = MAX(ONE-DBLE(K)*EPS,ZERO)*
     $                    (C1*WORK2(INDRV5+M))
                     IF (TAU2 .LT. TAU) EXIT
                  ENDDO
                  TAU = MAX(HALF*TAU,TAU2)
                  GO TO 1125
               ELSE
                  TMP3 = SQRT(TMP4)*SQRT(DFMA0(C1,WORK2(INDRV5+M),TAU))
               ENDIF
               
               IF (TMP3 .EQ. ZERO) THEN
                  CALL DLARTG3(SIGMA,-C1*WORK2(INDRV5+M),C1,S1)
               ELSE
                  CALL DLARTG6(C1*WORK2(INDRV5+M),TMP3,SIGMA,T,C1,S1)
               ENDIF
               WORK2(INDRV1+2*M) = C1
               WORK2(INDRV2+2*M) = S1
               A(M) = TMP3
*     
            ENDIF
*     
            SIGMA = T
            SIGMA2 = TOL*SIGMA
            DESIG = DESIG0
*     
            DO J = N-1, M, -1
               IF( ( WORK2(INDRV3+J+1).NE.ONE ) .OR. 
     $              ( WORK2(INDRV4+J+1).NE.ZERO ) ) THEN
                  CALL DROT(N0,SVT(J+1,1),LDSVT,SVT(J,1),LDSVT,
     $                 WORK2(INDRV3+J+1),WORK2(INDRV4+J+1))
               ENDIF
            ENDDO
*     
            DO J = N-1, M, -1
               IF( ( WORK2(INDRV3+J+1).NE.ONE ) .OR. 
     $              ( WORK2(INDRV4+J+1).NE.ZERO ) ) THEN
                  CALL DROT(N0,WORK(1,J+1),1,WORK(1,J),1,
     $                 WORK2(INDRV3+J+1),WORK2(INDRV4+J+1))
               ENDIF
               IF( ( WORK2(INDRV1+2*J+2).NE.ONE ) .OR. 
     $              ( WORK2(INDRV2+2*J+2).NE.ZERO ) ) 
     $              THEN
                  CALL DROT(N0,SU(1,J+1),1,WORK(1,J+1),1,
     $                 WORK2(INDRV1+2*J+2),WORK2(INDRV2+2*J+2))
               ENDIF
               IF( ( WORK2(INDRV1+2*J+1).NE.ONE ) .OR. 
     $              ( WORK2(INDRV2+2*J+1).NE.ZERO ) ) THEN
                  CALL DROT(N0,SU(1,J+1),1,SU(1,J),1,
     $                 WORK2(INDRV1+2*J+1),WORK2(INDRV2+2*J+1))
               ENDIF
            ENDDO
            IF( ( WORK2(INDRV1+2*M).NE.ONE ) .OR. 
     $           ( WORK2(INDRV2+2*M).NE.ZERO ) ) THEN
               CALL DROT(N0,SU(1,M),1,WORK(1,M),1,
     $              WORK2(INDRV1+2*M),WORK2(INDRV2+2*M))
            ENDIF
*     
            GO TO 410
*     
 1350       DO J = N-1, M, -1
               IF( ( WORK2(INDRV3+J+1).NE.ONE ) .OR. 
     $              ( WORK2(INDRV4+J+1).NE.ZERO ) ) THEN
                  CALL DROT(N0,SVT(J+1,1),LDSVT,SVT(J,1),LDSVT,
     $                 WORK2(INDRV3+J+1),WORK2(INDRV4+J+1))
               ENDIF
            ENDDO
*     
            DO J = N-1, M, -1
               IF( ( WORK2(INDRV3+J+1).NE.ONE ) .OR. 
     $              ( WORK2(INDRV4+J+1).NE.ZERO ) ) THEN
                  CALL DROT(N0,WORK(1,J+1),1,WORK(1,J),1,
     $                 WORK2(INDRV3+J+1),WORK2(INDRV4+J+1))
               ENDIF
            ENDDO

            GO TO 1395
            
 1390       IF (N-M+1 .EQ. 1) THEN
            
               CALL DLARTG2(WORK2(INDRV5+N),SIGMA,C1,S1,TMP1)
               
               IF( ( C1.NE.ONE ) .OR. ( S1.NE.ZERO ) ) THEN
                  CALL DROT(N0,SU(1,N),1,WORK(1,N),1,C1,S1)
               ENDIF
               
               CALL DLARTG7(SIGMA,DESIG,WORK2(INDRV5+N),A(N),DESIG0)
               
               GO TO 700
               
            ENDIF
            
            IF (N-M+1 .EQ. 2) THEN
               
               CALL DLASV2(WORK2(INDRV5+N-1), WORK2(INDRV6+N-1),
     $              WORK2(INDRV5+N), WORK2(INDRV5+N), WORK2(INDRV5+N-1),
     $              S1, C1, S2, C2)
               
               IF( ( C1.NE.ONE ) .OR. ( S1.NE.ZERO ) ) THEN
                  CALL DROT(N0,SVT(N-1,1),LDSVT,SVT(N,1),LDSVT,C1,S1)
               ENDIF
               
               IF( ( C2.NE.ONE ) .OR. ( S2.NE.ZERO ) ) THEN
                  CALL DROT(N0,SU(1,N-1),1,SU(1,N),1,C2,S2)
               ENDIF
               
               IF( ( C1.NE.ONE ) .OR. ( S1.NE.ZERO ) ) THEN
                  CALL DROT(N0,WORK(1,N-1),1,WORK(1,N),1,C1,S1)
               ENDIF
               
               CALL DLARTG2(WORK2(INDRV5+N-1),SIGMA,C1,S1,TMP1)
               
               IF( ( C1.NE.ONE ) .OR. ( S1.NE.ZERO ) ) THEN
                  CALL DROT(N0,SU(1,N-1),1,WORK(1,N-1),1,C1,S1)
               ENDIF
               
               CALL DLARTG7(SIGMA,DESIG,WORK2(INDRV5+N-1),A(N-1),DESIG0)
               
               CALL DLARTG2(WORK2(INDRV5+N),SIGMA,C1,S1,TMP1)
               
               IF( ( C1.NE.ONE ) .OR. ( S1.NE.ZERO ) ) THEN
                  CALL DROT(N0,SU(1,N),1,WORK(1,N),1,C1,S1)
               ENDIF
               
               CALL DLARTG7(SIGMA,DESIG,WORK2(INDRV5+N),A(N),DESIG0)
               
               GO TO 700
            ENDIF
            
 1395       IF (WORK2(INDRV6+N-1) .LE. SIGMA2) THEN
               
               CALL DLARTG2(WORK2(INDRV5+N),SIGMA,C1,S1,TMP1)
               
               IF( ( C1.NE.ONE ) .OR. ( S1.NE.ZERO ) ) THEN
                  CALL DROT(N0,SU(1,N),1,WORK(1,N),1,C1,S1)
               ENDIF
               
               CALL DLARTG7(SIGMA,DESIG,WORK2(INDRV5+N),A(N),DESIG0)
               
               N = N - 1
               GO TO 1390
               
            ENDIF
*     
            IF (WORK2(INDRV6+N-2) .LE. SIGMA2) THEN
               
               CALL DLASV2(WORK2(INDRV5+N-1), WORK2(INDRV6+N-1),
     $              WORK2(INDRV5+N), WORK2(INDRV5+N), WORK2(INDRV5+N-1),
     $              S1, C1, S2, C2)
               
               IF( ( C1.NE.ONE ) .OR. ( S1.NE.ZERO ) ) THEN
                  CALL DROT(N0,SVT(N-1,1),LDSVT,SVT(N,1),LDSVT,C1,S1)
               ENDIF
               
               IF( ( C2.NE.ONE ) .OR. ( S2.NE.ZERO ) ) THEN
                  CALL DROT(N0,SU(1,N-1),1,SU(1,N),1,C2,S2)
               ENDIF
               
               IF( ( C1.NE.ONE ) .OR. ( S1.NE.ZERO ) ) THEN
                  CALL DROT(N0,WORK(1,N-1),1,WORK(1,N),1,C1,S1)
               ENDIF
               
               CALL DLARTG2(WORK2(INDRV5+N-1),SIGMA,C1,S1,TMP1)
               
               IF( ( C1.NE.ONE ) .OR. ( S1.NE.ZERO ) ) THEN
                  CALL DROT(N0,SU(1,N-1),1,WORK(1,N-1),1,C1,S1)
               ENDIF
               
               CALL DLARTG7(SIGMA,DESIG,WORK2(INDRV5+N-1),A(N-1),DESIG0)
               
               CALL DLARTG2(WORK2(INDRV5+N),SIGMA,C1,S1,TMP1)
               
               IF( ( C1.NE.ONE ) .OR. ( S1.NE.ZERO ) ) THEN
                  CALL DROT(N0,SU(1,N),1,WORK(1,N),1,C1,S1)
               ENDIF
               
               CALL DLARTG7(SIGMA,DESIG,WORK2(INDRV5+N),A(N),DESIG0)
               
               N = N - 2
               GO TO 1390
               
            ENDIF

            TMP1 = WORK2(INDRV5+N)
            DO J = N-1, M, -1
               CALL DLARTG2(TMP1,WORK2(INDRV6+J),C1,S1,A(J+1))
               WORK2(INDRV1+J+1) = C1
               WORK2(INDRV2+J+1) = S1
               B(J) = S1*WORK2(INDRV5+J)
               TMP1 = C1*WORK2(INDRV5+J)
            ENDDO
            A(M) = TMP1
*     
            DO J = N-1, M, -1
               IF( ( WORK2(INDRV1+J+1).NE.ONE ) .OR. 
     $              ( WORK2(INDRV2+J+1).NE.ZERO ) ) THEN
                  CALL DROT(N0,SU(1,J+1),1,SU(1,J),1,
     $                 WORK2(INDRV1+J+1),WORK2(INDRV2+J+1))
               ENDIF
            ENDDO
*     
 410        TMP1 = A(N)
            IF (B(N-1) .LE. SIGMA2 .OR. TMP1 .EQ. B(N-1)+TMP1) THEN
               B(N-1) = ZERO
               TMP1 = A(N-1)
            ELSE
               TMP1 = A(N-1)*(TMP1/(TMP1+B(N-1)))
            ENDIF

            IF (B(N-2) .LE. SIGMA2 .OR. TMP1 .EQ. B(N-2)+TMP1) THEN
               B(N-2) = ZERO
               TMP1 = A(N-2)
            ELSE
               TMP1 = A(N-2)*(TMP1/(TMP1+B(N-2)))
            ENDIF

            M0 = M
            DO J = N-3, M, -1
               IF (B(J) .LE. SIGMA2 .OR. TMP1 .EQ. B(J)+TMP1) THEN
                  B(J) = -SIGMA
                  WORK2(INDRV6+J) = -DESIG
                  IF (M0 .EQ. M) M0 = J+1
                  TMP1 = A(J)
               ELSE
                  TMP1 = A(J)*(TMP1/(TMP1+B(J)))
               ENDIF
            ENDDO
            M = M0
*     
         ENDIF
*     
      ENDDO
      INFO = 2
      RETURN
      
 700  N = M-1
      IF (N .LE. 0) THEN
         GO TO 4000
      ENDIF
*     
      DO J = M-2,1,-1
         IF (B(J) .LE. ZERO) THEN
            M = J+1
            GO TO 3000
         ENDIF
      ENDDO
      M = 1
      GO TO 3000
*     
 4000 DO I = 1, N0 - 1
*     
*     Scan for smallest D(I)
*     
         ISUB = 1
         SMIN = A( 1 )
         DO J = 2, N0 + 1 - I
            IF( A( J ).LE.SMIN ) THEN
               ISUB = J
               SMIN = A( J )
            END IF
         ENDDO
         IF( ISUB.NE.N0+1-I ) THEN
*     
*     Swap singular values and vectors
*     
            A( ISUB ) = A( N0+1-I )
            A( N0+1-I ) = SMIN
            CALL DSWAP( N0, SVT( ISUB, 1 ), LDSVT, SVT( N0+1-I, 1 ), 
     $           LDSVT )
            CALL DSWAP( N0, SU( 1, ISUB ), 1, SU( 1, N0+1-I ), 1 )
         END IF
      ENDDO
*     
      INFO = 0
      RETURN
*     
      END
