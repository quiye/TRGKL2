SUBROUTINE RESGKL(J,MODE,MATDESCRA,INDXA,PNTRBA,PNTREA,A,M,N,K,L,BK,VK,UK,VPLUS,INFO,SELEK)
  IMPLICIT NONE

  INTEGER INDXA(*), PNTRBA(*), PNTREA(*)
  DOUBLE PRECISION A(*)
  CHARACTER*6 MATDESCRA
  CHARACTER MODE

  INTEGER M,N,K,L,I,J,LWORK,IINFO,INFO,SELEK
  DOUBLE PRECISION ONE,ZERO,MINUSONE
  PARAMETER(ONE=1.0D+0,ZERO=0.0D+0,MINUSONE=-1.0D+0)
  DOUBLE PRECISION NRM,CDUMMY(1,1),beta
  DOUBLE PRECISION BK(K,K),CPBK(K,K),VK(N,K),UK(M,K),VPLUS(N),VMTEMP(M),BD(K),BE(K)
  DOUBLE PRECISION VNTEMP(N),VL(K,K),TMPKK(K,K),KK(K,K)
  DOUBLE PRECISION DNRM2,VM(N,K),UM(M,K),HIGE(L),DLAMCH,WORK2(K*8)
  DOUBLE PRECISION Q(K,K),P(K,K),TAU(L),Q_arr(K)
  DOUBLE PRECISION,ALLOCATABLE :: WORK(:)
  EXTERNAL CGS2,DLAMCH,DOQDS2,DGEBRD_2

  ! LWORK = 5*K
  LWORK =MAX(5*K,MAX(M,N),K*K)
  ALLOCATE (WORK(LWORK))
  WORK = ONE
  WORK(1) = ONE+2.0D0*DLAMCH('E')
  
  VK(: ,J+1) = VPLUS
  NRM = 0.0D+0
  DO WHILE(J < K)
     IF(MODE=='s') THEN
        CALL MKL_DCSRMV( 'N', M, N, ONE, MATDESCRA, A, INDXA, &
             PNTRBA,PNTREA, VK(:,J+1), ZERO, VMTEMP(1))
     ELSE IF(MODE=='d') THEN
        CALL DGEMV('N',M,N,ONE,A,M,VK(:,J+1),1,ZERO,VMTEMP,1)
     END IF
     CALL CGS2(VMTEMP,UK,M,K,J)
     NRM = DNRM2(M,VMTEMP,1)
     VMTEMP = VMTEMP / NRM
     BK(J+1,J+1) = NRM
     UK(:,J+1) = VMTEMP
     IF(MODE=='s') THEN
        CALL MKL_DCSRMV( 'T', M, N, ONE, MATDESCRA, A, INDXA, &
             PNTRBA,PNTREA, UK(:,J+1), ZERO, VNTEMP(1))
     ELSE IF(MODE=='d') THEN
        CALL DGEMV('T',M,N,ONE,A,M,UK(:,J+1),1,ZERO,VNTEMP,1)
     END IF
     CALL CGS2(VNTEMP,VK,N,K,J+1)
     NRM = DNRM2(N,VNTEMP,1)
     VNTEMP = VNTEMP / NRM
     IF(J < K - 1) THEN
        BK(J+1,J+2) = NRM
        VK(: ,J+2) = VNTEMP
     END IF
     J = J+1
  END DO
  VPLUS = VNTEMP
  beta = NRM
  J = L

  cpbk=bk
  Q=0
  DO I = 1,K
     Q(I,I)=1
  END DO
  P=Q
  Q_arr=0
  Q_arr(k)=1.0D+0
   IF (SELEK == 1) THEN
     ! with lapack 1.0 QR
     !  WRITE(*,*) "GIVENS回転(DGEBRDG_LP1)+QR法(DBDSQRU)+両側(DGEMM)"
     CALL DGEBRDG_4_BISIDE(L+1,BK,K,Q,P)
     !上下左右逆転実装スべし
     DO I = 1,K
        BD(I)=BK(I,I)
     END DO
     DO I = 1,K-1
        BE(I)=BK(I,I+1)
     END DO
     BK=0
     CALL DBDSQRU( 'U',K,K,K,0,BD,BE,P,K,Q,K,CDUMMY,K,WORK,IINFO )
     DO I =1 ,L
        BK(I,I) = BD(k+1-i)
     END DO
     DO I = 1,L
        HIGE(I) = beta * Q(k,k+1-I)
     END DO
     DO I = 1, L
        Q_arr = Q(1:K,K+1-I)
        Q(1:K,K+1-I) = Q(1:K,I)
        Q(1:K,I) = Q_arr
        Q_arr = P(K+1-I,1:K)
        P(K+1-I,1:K)=P(I,1:K)
        P(I,1:K) = Q_arr
     END DO
     CALL DGEMM('N','N',M,L,K,ONE,UK,M,Q,K,ZERO,UM,M)
     CALL DGEMM('N','T',N,L,K,ONE,VK,N,P,K,ZERO,VM,N)
  ELSE IF ( SELEK==2 ) THEN
     ! with lapack 1.0 QR
     !  WRITE(*,*) "GIVENS回転(DGEBRDG_LP1)+QR法(DBDSQRU)+片側(DORMQR)"

     CALL DGEBRDG_LP1(L+1,BK,K,Q,P)

     DO I = 1,K
        BD(I)=BK(I,I)
     END DO

     DO I = 1,K-1
        BE(I)=BK(I,I+1)
     END DO

     CALL DBDSQRU('U',K,K,0,0,BD,BE,P,K,Q,K,CDUMMY,K,WORK,IINFO )
     VL=TRANSPOSE(P)
     CALL DGEQRF(K,L,VL,K,TAU,WORK,LWORK,IINFO)
     CALL DORMQR('R','N',K,K,L,VL,K,TAU,CPBK,K,WORK,LWORK,IINFO )
     CALL DORMQR('R','N',N,K,L,VL,K,TAU,VK,N,WORK,LWORK,IINFO )
     VM(:,1:L)=VK(:,1:L)
     CALL DGEQRF(K,L,CPBK,K,TAU,WORK,LWORK,IINFO)
     CALL DORMQR('R','N',M,K,L,CPBK,K,TAU,UK,M,WORK,LWORK,IINFO )
     UM(:,1:L)=UK(:,1:L)
     BK=0
     DO I = 1,L
        BK(I,I:L)=CPBK(I,I:L)
     END DO
     CALL DORMQR('R','N',1,K,L,CPBK,K,TAU,Q_arr,1,WORK,LWORK,IINFO )
     DO I = 1,L
        HIGE(I) = beta * Q_arr(i)
     END DO
  ELSE IF ( SELEK==6) THEN
     ! そもそも、初回はDGEBRDする必要がない
     ! 下のELSEから始まるものと比較して、ロバスト性が高くなってる。

     !  WRITE(*,*) "ヤコビ法(DGESVJ)+片側(DORMQR)"
     I=M
     CALL DGESVJ('U','C','V',K,K,BK,K,VNTEMP,I,VL,K,WORK,LWORK,IINFO)
     CALL DGEQRF(K,L,VL,K,TAU,WORK,LWORK,IINFO)
     CALL DORMQR('R','N',K,K,L,VL,K,TAU,CPBK,K,WORK,LWORK,IINFO )
     CALL DORMQR('R','N',N,K,L,VL,K,TAU,VK,N,WORK,LWORK,IINFO )
     VM=VK
     CALL DGEQRF(K,L,CPBK,K,TAU,WORK,LWORK,IINFO)
     CALL DORMQR('R','N',M,K,L,CPBK,K,TAU,UK,M,WORK,LWORK,IINFO )
     UM=UK
     BK=0
     DO I = 1,L
        BK(I,I:L)=CPBK(I,I:L)
     END DO
     CALL DORMQR('R','N',1,K,L,CPBK,K,TAU,Q_arr,1,WORK,LWORK,IINFO )
     DO I = 1,L
        HIGE(I) = beta * Q_arr(i)
     END DO

  ELSE IF (SELEK==5) THEN
     !  WRITE(*,*) "ヤコビ法(DGESVJ)+両側(DGEMM)"
     I=M
     CALL DGESVJ('U','C','V',K,K,BK,K,VNTEMP,I,VL,K,WORK,LWORK,IINFO)
     IF ( NINT(WORK(2)) .LT. L ) THEN
        INFO = 2
        RETURN
     END IF
     DO I = 1,L
        HIGE(I) = beta * BK(k,k+1-I)
     END DO

     DO I = 1,L
         Q_arr = VL(1:K,K+1-I)
         VL(1:K,K+1-I) = VL(1:K,I)
         VL(1:K,I) = Q_arr
         Q_arr = BK(1:K,K+1-I)
         BK(1:K,K+1-I) = BK(1:K,I)
         BK(1:K,I) = Q_arr
     END DO

     CALL DGEMM('N','N',N,L,K,ONE,VK,N,VL,K,ZERO,VM,N)
     CALL DGEMM('N','N',M,L,K,ONE,UK,M,BK,K,ZERO,UM,M)
     BK=0
     DO I =1 ,L
        BK(I,I) = VNTEMP(k+1-I)
     END DO
 
  ELSE IF(SELEK==4) THEN
     !  WRITE(*,*) "GIVENS回転(DGEBRDG_LP1)+OQDS2法(DOQDS2)+片側(DORMQR)"

     CALL DGEBRDG_LP1(L+1,BK,K,Q,P) ! OUT PUT P IS TRANSPOSED

     DO I = 1,K
        BD(I)=BK(I,I)
     END DO

     DO I = 1,K-1
        BE(I)=BK(I,I+1)
     END DO

     P=TRANSPOSE(P)
     CALL DOQDS2('L',K,BD,BE,P,K,WORK,WORK2,INFO)
     CALL DGEQRF(K,L,P,K,TAU,WORK,LWORK,IINFO)
     CALL DORMQR('R','N',K,K,L,P,K,TAU,CPBK,K,WORK,LWORK,IINFO )
     CALL DORMQR('R','N',N,K,L,P,K,TAU,VK,N,WORK,LWORK,IINFO )
     VM=VK
     CALL DGEQRF(K,L,CPBK,K,TAU,WORK,LWORK,IINFO)
     CALL DORMQR('R','N',M,K,L,CPBK,K,TAU,UK,M,WORK,LWORK,IINFO )
     UM=UK
     BK=0
     DO I = 1,L
        BK(I,I:L)=CPBK(I,I:L)
     END DO
     CALL DORMQR('R','N',1,K,L,CPBK,K,TAU,Q_arr,1,WORK,LWORK,IINFO )
     DO I = 1,L
        HIGE(I) = beta * Q_arr(i)
     END DO
  ELSE IF(SELEK==3) THEN
     !WRITE(*,*) "GIVENS回転(DGEBRDG_LP1)+OQDS1法(DOQDS1)+両側(DGEMM)"

     CALL DGEBRDG_4_BISIDE(L+1,BK,K,Q,P)
     P=TRANSPOSE(P)
     Q=TRANSPOSE(Q)
     DO I = 1,K
        BD(I)=BK(I,I)
     END DO
     DO I = 1,K-1
        BE(I)=BK(I,I+1)
     END DO
     CALL DOQDS1('L',K,BD,BE,P,K,Q,K,WORK,WORK2,INFO)
     BK=0
     DO I =1 ,L
        BK(I,I) = BD(k-I+1)
     END DO
     DO I = 1,L
        HIGE(I) = beta * Q(K-I+1,K)
     END DO
     DO I = 1,L
        Q_arr = Q(K+1-I,1:K)
        Q(K+1-I,1:K) = Q(I,1:K)
        Q(I,1:K) = Q_arr
        Q_arr = P(1:K,K+1-I)
        P(1:K,K+1-I) = P(1:K,I)
        P(1:K,I) = Q_arr
     END DO
     CALL DGEMM('N','T',M,L,K,ONE,UK,M,Q,K,ZERO,UM,M)
     CALL DGEMM('N','N',N,L,K,ONE,VK,N,P,K,ZERO,VM,N)
  ELSE IF(SELEK==7) THEN
     !WRITE(*,*) "GIVENS回転(DGEBRDG_LP1)+OQDS1法(DOQDS1)+両側(DGEMM)"

     CALL DGEBRDG_4_BISIDE(L+1,BK,K,Q,P)
     TMPKK=0
     DO I = 1, K
        TMPKK(I,K-I+1) = 1.0D+0
     END DO
     !CALL DGEMM('N','N',K,K,K,ONE,BK,K,TMPKK,K,ZERO,KK,K)
     !CALL DGEMM('N','N',K,K,K,ONE,TMPKK,K,KK,K,ZERO,BK,K)
     CALL DGEMM('N','N',K,K,K,ONE,Q,K,TMPKK,K,ZERO,KK,K)
     Q=KK
     CALL DGEMM('N','N',K,K,K,ONE,TMPKK,K,P,K,ZERO,KK,K)
     P=KK
     !P=TRANSPOSE(P)
     !Q=TRANSPOSE(Q)
     DO I = 1,K
        BD(I)=BK(K+1-I,K+1-I)
     END DO
     DO I = 1,K-1
        BE(I)=BK(K-I,K+1-I)
     END DO
     CALL DOQDS1('L',K,BD,BE,Q,K,P,K,WORK,WORK2,INFO)
     P=TRANSPOSE(P)
     Q=TRANSPOSE(Q)
     BK=0
     DO I =1 ,L
        BK(I,I) = BD(k-I+1)
     END DO
     DO I = 1,L
        HIGE(I) = beta * Q(K-I+1,K)
     END DO
     DO I = 1,L
        Q_arr = Q(K+1-I,1:K)
        Q(K+1-I,1:K) = Q(I,1:K)
        Q(I,1:K) = Q_arr
        Q_arr = P(1:K,K+1-I)
        P(1:K,K+1-I) = P(1:K,I)
        P(1:K,I) = Q_arr
     END DO
     CALL DGEMM('N','T',M,L,K,ONE,UK,M,Q,K,ZERO,UM,M)
     CALL DGEMM('N','N',N,L,K,ONE,VK,N,P,K,ZERO,VM,N)
   END IF
  !IF(MODE=='s') THEN
     !CALL MKL_DCSRMV( 'N', M, N, ONE, MATDESCRA, A, INDXA, &
          !PNTRBA,PNTREA, VPLUS(1), ZERO, VMTEMP(1))
  !ELSE IF(MODE=='d') THEN
     !CALL DGEMV('N',M,N,ONE,A,M,VPLUS,1,ZERO,VMTEMP,1)
  !END IF
  !CALL DGEMV('T',M,L,ONE,UM,M,VMTEMP,1,ZERO,HIGE,1)
  BK(1: L,L+1) = HIGE
  VK = VM
  UK = UM
  RETURN
END SUBROUTINE RESGKL
