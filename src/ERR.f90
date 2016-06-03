DOUBLE PRECISION FUNCTION ERR(mode,matdescra,indxA,pntrbA,pntreA,A,M,N,K,L,BK,VK,UK)
  IMPLICIT NONE

  INTEGER indxA(*), pntrbA(*), pntreA(*)
  DOUBLE PRECISION A(*)
  CHARACTER*6 matdescra
  CHARACTER MODE

  INTEGER M,N,K,L,DM,DMM,S
  DOUBLE PRECISION ONE,ZERO,SUM
  PARAMETER(ONE=1.0D+0,ZERO=0.0D+0)
  DOUBLE PRECISION DL
  DOUBLE PRECISION BK(K,K),VK(N,K),UK(M,K),VMTEMP(M)
  DOUBLE PRECISION VNTEMP(N)
  DOUBLE PRECISION DNRM2
  EXTERNAL DNRM2
  DM = 10 * K
  DMM = K
  DL = L

  SUM=0
  DO S=0,L-1
     IF(MODE=='s') THEN
        CALL MKL_DCSRMV( 'N', M, N, ONE, matdescra, A, indxA, &
             pntrbA,pntreA, VK(:,S+1), ZERO, VMTEMP(1))
     ELSE IF(MODE=='d') THEN
        CALL DGEMV('N',M,N,ONE,A,M,VK(:,S+1),1,ZERO,VMTEMP,1)
     END IF
     VMTEMP = VMTEMP - BK(S+1,S+1)*UK(:,S+1)
     IF(MODE=='s') THEN
        CALL MKL_DCSRMV( 'T', M, N, ONE, matdescra, A, indxA, &
             pntrbA,pntreA, UK(:,S+1), ZERO, VNTEMP(1))
     ELSE IF(MODE=='d') THEN
        CALL DGEMV('T',M,N,ONE,A,M,UK(:,S+1),1,ZERO,VNTEMP,1)
     END IF
     VNTEMP = VNTEMP - BK(S+1,S+1)*VK(:,S+1)
     SUM = SUM + DNRM2(M,VMTEMP,1) !収束が早い
     SUM = SUM + DNRM2(N,VNTEMP,1) !収束が遅い
  END DO
  ERR = SUM / DL
  RETURN
END FUNCTION ERR
