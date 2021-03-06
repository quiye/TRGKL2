      SUBROUTINE DOQDS1(UPLO,N0,A,B,SU,LDSU,SVT,LDSVT,WORK,WORK2,
     $     INFO)
      
      IMPLICIT NONE
      
      CHARACTER UPLO
      INTEGER N0, INFO, LDSU, LDSVT, ISUB, MAXITER, SIT, M0, FLAG
      DOUBLE PRECISION A(*), B(*), WORK(N0,*)
      DOUBLE PRECISION SU(LDSU, *), SVT(LDSVT, *), WORK2(*)
      INTEGER N, M, I, J, K, OLDM, OLDN
      INTEGER INDRV1, INDRV2, INDRV3, INDRV4, INDRV5, INDRV6, INDRV7,
     $     INDRV8
      DOUBLE PRECISION TMP1, TMP2, TMP3, TMP4, TMP5
      DOUBLE PRECISION TAU, TAU1, TAU2, TAU3
      DOUBLE PRECISION SIGMA, SIGMA2, DESIG, T, DESIG0, S
      DOUBLE PRECISION C1, S1, C2, S2, SMIN, EPS, TOL
      
      DOUBLE PRECISION ONE, ZERO, HALF, TWO, CONST
      INTEGER MAXREC
      PARAMETER (ONE = 1.0D0, ZERO = 0.0D0, HALF = 0.5D0, TWO = 2.0D0)
      PARAMETER (CONST = 0.75D0)
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
      INDRV2 = INDRV1+N0
      INDRV3 = INDRV2+N0
      INDRV4 = INDRV3+N0
      INDRV5 = INDRV4+N0
      INDRV6 = INDRV5+N0
      INDRV7 = INDRV6+N0
      INDRV8 = INDRV7+N0
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
         
         K=N0/2
         DO J=1,K
            TMP1=A(J)
            A(J)=A(N0+1-J)
            A(N0+1-J)=TMP1
            TMP1=B(J)
            B(J)=B(N0-J)
            B(N0-J)=TMP1
            CALL DSWAP(N0,SU(1,J),1,SU(1,N0+1-J),1)
            CALL DSWAP(N0,SVT(J,1),LDSVT,SVT(N0+1-J,1),LDSVT)
         ENDDO
            
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
      ENDIF
*
      OLDM = -1
      OLDN = -1
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
      IF (N-M+1 .GE. 3) THEN

         IF (A(M) .GE. A(N)) THEN
         
            TMP1 = A(M)
            DO J = M, N-3
               IF (B(J) .LE. TOL*TMP1) THEN
                  B(J) = -ZERO
                  WORK2(INDRV6+J) = -ZERO
                  M = J+1
                  TMP1 = A(J+1)
               ELSE
                  TMP1 = A(J+1)*(TMP1/(TMP1+B(J)))
               ENDIF
            ENDDO
            
            IF (B(N-2) .LE. TOL*TMP1) THEN
               B(N-2) = ZERO
               TMP1 = A(N-1)
            ELSE
               TMP1 = A(N-1)*(TMP1/(TMP1+B(N-2)))
            ENDIF
            
            IF (B(N-1) .LE. TOL*TMP1) THEN
               B(N-1) = ZERO
            ENDIF
            
         ELSE
            
            TMP1 = A(N)
            IF (B(N-1) .LE. TOL*TMP1) THEN
               B(N-1) = ZERO
               TMP1 = A(N-1)
            ELSE
               TMP1 = A(N-1)*(TMP1/(TMP1+B(N-1)))
            ENDIF
            
            IF (B(N-2) .LE. TOL*TMP1) THEN
               B(N-2) = ZERO
               TMP1 = A(N-2)
            ELSE
               TMP1 = A(N-2)*(TMP1/(TMP1+B(N-2)))
            ENDIF
            
            M0 = M
            DO J = N-3, M, -1
               IF (B(J) .LE. TOL*TMP1) THEN
                  B(J) = -ZERO
                  WORK2(INDRV6+J) = -ZERO
                  IF (M0 .EQ. M) M0 = J+1
                  TMP1 = A(J)
               ELSE
                  TMP1 = A(J)*(TMP1/(TMP1+B(J)))
               ENDIF
            ENDDO
            M = M0

         ENDIF
            
      ENDIF
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
            
            CALL DLARTG(A(N),SIGMA,C1,S1,TMP1)
            
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
            
            CALL DLARTG(A(N-1),SIGMA,C1,S1,TMP1)
            
            IF( ( C1.NE.ONE ) .OR. ( S1.NE.ZERO ) ) THEN
               CALL DROT(N0,SU(1,N-1),1,WORK(1,N-1),1,C1,S1)
            ENDIF
            
            CALL DLARTG7(SIGMA,DESIG,A(N-1),A(N-1),DESIG0)
            
            CALL DLARTG(A(N),SIGMA,C1,S1,TMP1)
            
            IF( ( C1.NE.ONE ) .OR. ( S1.NE.ZERO ) ) THEN
               CALL DROT(N0,SU(1,N),1,WORK(1,N),1,C1,S1)
            ENDIF
            
            CALL DLARTG7(SIGMA,DESIG,A(N),A(N),DESIG0)
            
            GO TO 700
         ENDIF

         CALL DLARTG(A(N),SIGMA,C1,S1,T)

         IF (B(N-1) .LE. TOL*T) THEN
         
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
            
            CALL DLARTG(A(N-1),SIGMA,C1,S1,TMP1)
            
            IF( ( C1.NE.ONE ) .OR. ( S1.NE.ZERO ) ) THEN
               CALL DROT(N0,SU(1,N-1),1,WORK(1,N-1),1,C1,S1)
            ENDIF
            
            CALL DLARTG7(SIGMA,DESIG,A(N-1),A(N-1),DESIG0)
            
            CALL DLARTG(A(N),SIGMA,C1,S1,TMP1)
            
            IF( ( C1.NE.ONE ) .OR. ( S1.NE.ZERO ) ) THEN
               CALL DROT(N0,SU(1,N),1,WORK(1,N),1,C1,S1)
            ENDIF
            
            CALL DLARTG7(SIGMA,DESIG,A(N),A(N),DESIG0)
            
            N = N - 2
            GO TO 15
            
         ENDIF
*
         IF (M .GT. OLDN .OR. N .LT. OLDM) THEN
            IF ( A(M) .LT. A(N) ) THEN

               IF (SIGMA .GT. ZERO) THEN
                  
                  TMP1 = A(N)
                  CALL DLARTG(SIGMA,TMP1,C1,S1,T)
                  IF (B(N-1)+T .EQ. T) THEN
                     B(N-1) = ZERO
                     GO TO 15
                  ELSE
                     TMP1 = A(N-1)*(TMP1/(TMP1+B(N-1)))
                  ENDIF
                  
                  CALL DLARTG(SIGMA,TMP1,C1,S1,T)
                  IF (B(N-2)+T .EQ. T) THEN
                     B(N-2) = ZERO
                     GO TO 15
                  ELSE
                     TMP1 = A(N-2)*(TMP1/(TMP1+B(N-2)))
                  ENDIF
                  
                  M0 = M
                  DO J = N-3, M, -1
                     CALL DLARTG(SIGMA,TMP1,C1,S1,T)
                     IF (B(J)+T .EQ. T) THEN
                        B(J) = -SIGMA
                        WORK2(INDRV6+J) = -DESIG
                        IF (M0 .EQ. M) M0 = J+1
                        TMP1 = A(J)
                     ELSE
                        TMP1 = A(J)*(TMP1/(TMP1+B(J)))
                     ENDIF
                  ENDDO

               ELSE
                  
                  TMP1 = A(N)
                  IF (B(N-1) .LE. TOL*TMP1) THEN
                     B(N-1) = ZERO
                     GO TO 15
                  ELSE
                     TMP1 = A(N-1)*(TMP1/(TMP1+B(N-1)))
                  ENDIF
                  
                  IF (B(N-2) .LE. TOL*TMP1) THEN
                     B(N-2) = ZERO
                     GO TO 15
                  ELSE
                     TMP1 = A(N-2)*(TMP1/(TMP1+B(N-2)))
                  ENDIF
                  
                  M0 = M
                  DO J = N-3, M, -1
                     IF (B(J) .LE. TOL*TMP1) THEN
                        B(J) = -SIGMA
                        WORK2(INDRV6+J) = -DESIG
                        IF (M0 .EQ. M) M0 = J+1
                        TMP1 = A(J)
                     ELSE
                        TMP1 = A(J)*(TMP1/(TMP1+B(J)))
                     ENDIF
                  ENDDO

               ENDIF
               
               IF (M0 .NE. M) THEN
                  M = M0
                  GO TO 15
               ENDIF

               K=(N-M+1)/2+M-1
               DO J=M,K
                  TMP1=A(J)
                  A(J)=A(N+M-J)
                  A(N+M-J)=TMP1
                  TMP1=B(J)
                  B(J)=B(N-1+M-J)
                  B(N-1+M-J)=TMP1
                  CALL DSWAP(N0,SU(1,J),1,SU(1,N+M-J),1)
                  CALL DSWAP(N0,WORK(1,J),1,WORK(1,N+M-J),1)
                  CALL DSWAP(N0,SVT(J,1),LDSVT,SVT(N+M-J,1),LDSVT)
               ENDDO

               TMP1 = A(M)
               DO J = M, N-1
                  CALL DLARTG(TMP1,B(J),C1,S1,A(J))
                  WORK2(INDRV3+J) = C1
                  WORK2(INDRV4+J) = S1
                  B(J) = S1*A(J+1)
                  TMP1 = C1*A(J+1)
               ENDDO
               A(N) = TMP1
*     
               CALL DLASR( 'R', 'V', 'F', N0, N-M+1,
     $              WORK2(INDRV3+M), WORK2(INDRV4+M),
     $              WORK(1,M), N0 )
*     
               CALL DLASR( 'L', 'V', 'F', N-M+1, N0,
     $              WORK2( INDRV3+M ),
     $              WORK2( INDRV4+M ), SVT( M, 1 ), N0 )
               
               GO TO 400
            ENDIF
         ENDIF
*     
         CALL DLARTG(SIGMA,A(N),C1,S1,T)
         IF (T .LE. SIGMA) GO TO 350

         TAU2 = MINVAL(A(M:N-1))
         CALL DLARTG(SIGMA,TAU2,C1,S1,T)
         IF (T .LE. SIGMA) GO TO 350
*
         CALL DLAS2(A(N-1), B(N-1), A(N), TAU, TMP3)
         TAU = MIN(TAU,A(N))
         IF (TAU .EQ. ZERO) GO TO 350
         
         IF (TMP3 .GE. TWO*TAU) THEN
            SIT = 1
         ELSE
            SIT = 0
         ENDIF

         CALL DLARTG(SIGMA,TAU,C1,S1,T)
         IF (T .LE. SIGMA .AND. SIT .EQ. 1) GO TO 350
         CALL DLARTG7(SIGMA,DESIG,TAU,T,DESIG0)
         IF (T .LE. SIGMA .AND. SIT .EQ. 1) GO TO 350
         
         IF (TAU2 .LE. TAU) THEN
            CALL DLARTG7(SIGMA,DESIG,TAU2,T,DESIG0)
            IF (T .LE. SIGMA .AND. SIT .EQ. 1) GO TO 350
            TAU1 = TAU2
            GO TO 160
         ELSE
            TAU1 = TAU
         ENDIF
*     
         TMP4 = A(M)-TAU
         TMP5 = A(M)+TAU
         IF (TMP4 .LE. ZERO) THEN
            GO TO 160
         ELSE
            TMP3 = SQRT(TMP4)*SQRT(TMP5)
         ENDIF
         DO J = M, N-2
            CALL DLARTG(TMP3,B(J),C1,S1,T)
            TMP4 = DFMA0(C1,A(J+1),-TAU)
            TMP5 = DFMA0(C1,A(J+1),TAU)
            IF (TMP4 .LE. ZERO) THEN
               GO TO 160
            ELSE
               TMP3 = SQRT(TMP4)*SQRT(TMP5)
            ENDIF
         ENDDO
         CALL DLARTG(TMP3,B(N-1),C1,S1,T)
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
 160     WORK2(INDRV5+N) = ONE/A(N)
         TMP3 = WORK2(INDRV5+N)
         DO J = N-1,M,-1
            WORK2(INDRV3+J) = B(J)/A(J)
            WORK2(INDRV5+J) = ONE/A(J)+
     $           WORK2(INDRV3+J)*WORK2(INDRV5+J+1)
            TMP3 = MAX(TMP3,WORK2(INDRV5+J))
         ENDDO
         WORK2(INDRV5+M) = (WORK2(INDRV5+M)/TMP3)/A(M)
         TMP1 = WORK2(INDRV5+M)
         DO J = M+1,N
            WORK2(INDRV4+J) = B(J-1)/A(J)
            WORK2(INDRV5+J) = (WORK2(INDRV5+J)/TMP3)/A(J)+
     $           WORK2(INDRV4+J)*WORK2(INDRV5+J-1)
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
     $           WORK2(INDRV3+J)*WORK2(INDRV6+J+1)
            TMP3 = MAX(TMP3,WORK2(INDRV6+J))
         ENDDO
         TMP2 = (WORK2(INDRV6+M)/TMP3)/A(M)
         TMP1 = TMP2/WORK2(INDRV5+M)
         DO J = M+1,N
            TMP2 = (WORK2(INDRV6+J)/TMP3)/A(J)+WORK2(INDRV4+J)*TMP2
            TMP1 = MAX(TMP1,TMP2/WORK2(INDRV5+J))
         ENDDO
         IF (TMP3 .GT. ZERO .AND. TMP1 .GT. ZERO) THEN
            TAU = MAX(TAU,ONE/SQRT(TMP1)/SQRT(TMP3))
         ENDIF
         
         TAU = MIN(TAU,TAU1)
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
 125     IF (TAU .EQ. ZERO) GO TO 350
         CALL DLARTG(SIGMA,TAU,C1,S1,T)
         IF (T .LE. SIGMA .AND. SIT .EQ. 1) GO TO 350
         CALL DLARTG7(SIGMA,DESIG,TAU,T,DESIG0)
         IF (T .LE. SIGMA .AND. SIT .EQ. 1) GO TO 350
*     
         IF (SIGMA .EQ. ZERO) THEN
*     
            TMP4 = A(M)-TAU
            TMP5 = A(M)+TAU
            IF (TMP4 .LE. ZERO) THEN
               DO K=0,MAXREC
                  TAU2 = MAX(ONE-DBLE(K)*EPS,ZERO)*A(M)
                  IF (TAU2 .LT. TAU) EXIT
               ENDDO
               TAU = TAU2
               GO TO 125
            ELSE
               TMP3 = SQRT(TMP4)*SQRT(TMP5)
            ENDIF
            
            CALL DLARTG(TMP3,T,C1,S1,S)
            WORK2(INDRV1+M) = C1
            WORK2(INDRV2+M) = -S1
            DO J = M, N-2
               CALL DLARTG(TMP3,B(J),C1,S1,WORK2(INDRV5+J))
               WORK2(INDRV7+J) = C1
               WORK2(INDRV8+J) = S1
               WORK2(INDRV6+J) = S1*A(J+1)
               
               TMP4 = DFMA0(C1,A(J+1),-TAU)
               TMP5 = DFMA0(C1,A(J+1),TAU)
               IF (TMP4 .LE. ZERO) THEN
                  DO K=0,MAXREC
                     TAU2 = MAX(ONE-DBLE(K)*EPS,ZERO)*(C1*A(J+1))
                     IF (TAU2 .LT. TAU) EXIT
                  ENDDO
                  TAU3 = CONST*TAU
                  IF (TAU3 .GE. TAU) TAU3 = HALF*TAU
                  TAU = MAX(TAU3,TAU2)
                  GO TO 125
               ELSE
                  TMP3 = SQRT(TMP4)*SQRT(TMP5)
               ENDIF
               
               CALL DLARTG(TMP3,T,C1,S1,S)
               WORK2(INDRV1+J+1) = C1
               WORK2(INDRV2+J+1) = -S1
            ENDDO
            CALL DLARTG(TMP3,B(N-1),C1,S1,WORK2(INDRV5+N-1))
            WORK2(INDRV7+N-1) = C1
            WORK2(INDRV8+N-1) = S1
            WORK2(INDRV6+N-1) = S1*A(N)
            
            TMP4 = DFMA0(C1,A(N),-TAU)
            TMP5 = DFMA0(C1,A(N),TAU)
            IF (TMP4 .LT. ZERO .AND. TAU+TMP4 .EQ. TAU) TMP4 = ZERO
            IF (TMP4 .LT. ZERO) THEN
               DO K=0,MAXREC
                  TAU2 = MAX(ONE-DBLE(K)*EPS,ZERO)*(C1*A(N))
                  IF (TAU2 .LT. TAU) EXIT
               ENDDO
               TAU3 = CONST*TAU
               IF (TAU3 .GE. TAU) TAU3 = HALF*TAU
               TAU = MAX(TAU3,TAU2)
               GO TO 125
            ELSE
               TMP3 = SQRT(TMP4)*SQRT(TMP5)
            ENDIF
            
            CALL DLARTG(TMP3,T,C1,S1,S)
            WORK2(INDRV1+N) = C1
            WORK2(INDRV2+N) = -S1
            WORK2(INDRV5+N) = TMP3
*     
         ELSE
*     
            TMP4 = A(M)-TAU
            TMP5 = A(M)+TAU
            IF (TMP4 .LE. ZERO) THEN
               DO K=0,MAXREC
                  TAU2 = MAX(ONE-DBLE(K)*EPS,ZERO)*A(M)
                  IF (TAU2 .LT. TAU) EXIT
               ENDDO
               TAU = TAU2
               GO TO 125
            ELSE
               TMP3 = SQRT(TMP4)*SQRT(TMP5)
            ENDIF
            
            CALL DLARTG6(A(M),TMP3,SIGMA,T,C1,S1)
            WORK2(INDRV1+M) = C1
            WORK2(INDRV2+M) = S1
            DO J = M, N-2
               CALL DLARTG(TMP3,B(J),C1,S1,WORK2(INDRV5+J))
               WORK2(INDRV7+J) = C1
               WORK2(INDRV8+J) = S1
               WORK2(INDRV6+J) = S1*A(J+1)
               
               TMP4 = DFMA0(C1,A(J+1),-TAU)
               TMP5 = DFMA0(C1,A(J+1),TAU)
               IF (TMP4 .LE. ZERO) THEN
                  DO K=0,MAXREC
                     TAU2 = MAX(ONE-DBLE(K)*EPS,ZERO)*(C1*A(J+1))
                     IF (TAU2 .LT. TAU) EXIT
                  ENDDO
                  TAU3 = CONST*TAU
                  IF (TAU3 .GE. TAU) TAU3 = HALF*TAU
                  TAU = MAX(TAU3,TAU2)
                  GO TO 125
               ELSE
                  TMP3 = SQRT(TMP4)*SQRT(TMP5)
               ENDIF
               
               CALL DLARTG6(C1*A(J+1),TMP3,SIGMA,T,C1,S1)
               WORK2(INDRV1+J+1) = C1
               WORK2(INDRV2+J+1) = S1
            ENDDO
            CALL DLARTG(TMP3,B(N-1),C1,S1,WORK2(INDRV5+N-1))
            WORK2(INDRV7+N-1) = C1
            WORK2(INDRV8+N-1) = S1
            WORK2(INDRV6+N-1) = S1*A(N)
            
            TMP4 = DFMA0(C1,A(N),-TAU)
            TMP5 = DFMA0(C1,A(N),TAU)
            IF (TMP4 .LT. ZERO .AND. TAU+TMP4 .EQ. TAU) TMP4 = ZERO
            IF (TMP4 .LT. ZERO) THEN
               DO K=0,MAXREC
                  TAU2 = MAX(ONE-DBLE(K)*EPS,ZERO)*(C1*A(N))
                  IF (TAU2 .LT. TAU) EXIT
               ENDDO
               TAU3 = CONST*TAU
               IF (TAU3 .GE. TAU) TAU3 = HALF*TAU
               TAU = MAX(TAU3,TAU2)
               GO TO 125
            ELSE
               TMP3 = SQRT(TMP4)*SQRT(TMP5)
            ENDIF
            
            IF (TMP3 .EQ. ZERO) THEN
               CALL DLARTG(SIGMA,C1*A(N),C1,S1,S)
               S1 = -S1
            ELSE
               CALL DLARTG6(C1*A(N),TMP3,SIGMA,T,C1,S1)
            ENDIF
            WORK2(INDRV1+N) = C1
            WORK2(INDRV2+N) = S1
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
            CALL DLARTG(TMP1,WORK2(INDRV6+J),C1,S1,A(J))
            WORK2(INDRV3+J) = C1
            WORK2(INDRV4+J) = S1
            B(J) = S1*WORK2(INDRV5+J+1)
            TMP1 = C1*WORK2(INDRV5+J+1)
         ENDDO
         A(N) = TMP1
*     
         IF( ( WORK2(INDRV1+M).NE.ONE ) .OR. 
     $        ( WORK2(INDRV2+M).NE.ZERO ) ) THEN
            CALL DROT(N0,SU(1,M),1,WORK(1,M),1,
     $           WORK2(INDRV1+M),WORK2(INDRV2+M))
         ENDIF
         DO J = M, N-1
            IF( ( WORK2(INDRV7+J).NE.ONE ) .OR. 
     $           ( WORK2(INDRV8+J).NE.ZERO ) ) THEN
               CALL DROT(N0,SU(1,J),1,SU(1,J+1),1,WORK2(INDRV7+J),
     $              WORK2(INDRV8+J))
            ENDIF
            IF( ( WORK2(INDRV1+J+1).NE.ONE ) .OR. 
     $           ( WORK2(INDRV2+J+1).NE.ZERO ) ) 
     $           THEN
               CALL DROT(N0,SU(1,J+1),1,WORK(1,J+1),1,
     $              WORK2(INDRV1+J+1),WORK2(INDRV2+J+1))
            ENDIF
            IF( ( WORK2(INDRV3+J).NE.ONE ) .OR. 
     $              ( WORK2(INDRV4+J).NE.ZERO ) ) THEN
               CALL DROT(N0,WORK(1,J),1,WORK(1,J+1),1,
     $              WORK2(INDRV3+J),WORK2(INDRV4+J))
            ENDIF
         ENDDO
*     
         CALL DLASR( 'L', 'V', 'F', N-M+1, N0,
     $        WORK2( INDRV3+M ),
     $        WORK2( INDRV4+M ), SVT( M, 1 ), N0 )
*
         GO TO 400
*     
 350     FLAG = -1
         TMP1 = A(M)
         CALL DLARTG(SIGMA,TMP1,C1,S1,T)
         IF (T .LE. SIGMA) THEN
            FLAG = M-1
            WORK2(INDRV1+M) = C1
            WORK2(INDRV2+M) = -S1
            TMP1 = ZERO
         ELSE
            WORK2(INDRV1+M) = ONE
            WORK2(INDRV2+M) = ZERO
         ENDIF
         DO J = M, N-1
            CALL DLARTG(TMP1,B(J),C1,S1,A(J))
            WORK2(INDRV7+J) = C1
            WORK2(INDRV8+J) = S1
            B(J) = S1*A(J+1)
            TMP1 = C1*A(J+1)
            CALL DLARTG(SIGMA,TMP1,C1,S1,T)
            IF (T .LE. SIGMA) THEN
               IF (FLAG .EQ. -1) FLAG = J
               WORK2(INDRV1+J+1) = C1
               WORK2(INDRV2+J+1) = -S1
               TMP1 = ZERO
            ELSE
               WORK2(INDRV1+J+1) = ONE
               WORK2(INDRV2+J+1) = ZERO
            ENDIF
         ENDDO
         A(N) = TMP1
*
         TMP1 = A(M)
         DO J = M, N-1
            CALL DLARTG(TMP1,B(J),C1,S1,A(J))
            WORK2(INDRV3+J) = C1
            WORK2(INDRV4+J) = S1
            B(J) = S1*A(J+1)
            TMP1 = C1*A(J+1)
         ENDDO
         A(N) = TMP1
*
         IF (FLAG .EQ. -1) THEN
            
            CALL DLASR( 'R', 'V', 'F', N0, N-M+1,
     $           WORK2(INDRV7+M), WORK2(INDRV8+M),
     $           SU(1,M), N0 )
            
         ELSE

            CALL DLASR( 'R', 'V', 'F', N0, (FLAG+1)-M+1,
     $           WORK2(INDRV7+M), WORK2(INDRV8+M),
     $           SU(1,M), N0 )
            
            IF( ( WORK2(INDRV1+FLAG+1).NE.ONE ) .OR.
     $           ( WORK2(INDRV2+FLAG+1).NE.ZERO ) )
     $           THEN
               CALL DROT(N0,SU(1,FLAG+1),1,WORK(1,FLAG+1),1,
     $              WORK2(INDRV1+FLAG+1),WORK2(INDRV2+FLAG+1))
            ENDIF

            CALL DLASR( 'R', 'V', 'F', N0, N-(FLAG+1)+1,
     $           WORK2(INDRV7+FLAG+1), WORK2(INDRV8+FLAG+1),
     $           SU(1,FLAG+1), N0 )

         ENDIF
*     
         CALL DLASR( 'R', 'V', 'F', N0, N-M+1,
     $        WORK2(INDRV3+M), WORK2(INDRV4+M),
     $        WORK(1,M), N0 )

*
         CALL DLASR( 'L', 'V', 'F', N-M+1, N0,
     $        WORK2( INDRV3+M ),
     $        WORK2( INDRV4+M ), SVT( M, 1 ), N0 )
*     
 400     OLDM = M
         OLDN = N
*
         IF (SIGMA .GT. ZERO) THEN

            TMP1 = A(M)
            DO J = M, N-3
               CALL DLARTG(SIGMA,TMP1,C1,S1,T)
               IF (B(J) .LE. SIGMA2 .OR. B(J)+T .EQ. T) THEN
                  B(J) = -SIGMA
                  WORK2(INDRV6+J) = -DESIG
                  M = J+1
                  TMP1 = A(J+1)
               ELSE
                  TMP1 = A(J+1)*(TMP1/(TMP1+B(J)))
               ENDIF
            ENDDO
            
            CALL DLARTG(SIGMA,TMP1,C1,S1,T)
            IF (B(N-2) .LE. SIGMA2 .OR. B(N-2)+T .EQ. T) THEN
               B(N-2) = ZERO
               TMP1 = A(N-1)
            ELSE
               TMP1 = A(N-1)*(TMP1/(TMP1+B(N-2)))
            ENDIF
            
            CALL DLARTG(SIGMA,TMP1,C1,S1,T)
            IF (B(N-1) .LE. SIGMA2 .OR. B(N-1)+T .EQ. T) THEN
               B(N-1) = ZERO
            ENDIF

         ELSE
            
            TMP1 = A(M)
            DO J = M, N-3
               IF (B(J) .LE. TOL*TMP1) THEN
                  B(J) = -SIGMA
                  WORK2(INDRV6+J) = -DESIG
                  M = J+1
                  TMP1 = A(J+1)
               ELSE
                  TMP1 = A(J+1)*(TMP1/(TMP1+B(J)))
               ENDIF
            ENDDO
            
            IF (B(N-2) .LE. TOL*TMP1) THEN
               B(N-2) = ZERO
               TMP1 = A(N-1)
            ELSE
               TMP1 = A(N-1)*(TMP1/(TMP1+B(N-2)))
            ENDIF
            
            IF (B(N-1) .LE. TOL*TMP1) THEN
               B(N-1) = ZERO
            ENDIF

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
