DOUBLE PRECISION FUNCTION ERR(mode,IAP,JA,A,M,N,K,L,BK,VK,UK,VTEMP)
  IMPLICIT NONE
  CHARACTER MODE
  INTEGER M,N,K,L,I
  INTEGER IAP(*),JA(*)
  DOUBLE PRECISION ONE,ZERO,DNRM2
  PARAMETER(ONE = 1.0D+0,ZERO = 0.0D+0)
  DOUBLE PRECISION A(*),BK(K,K),VK(N,K),UK(M,K),VTEMP(*)
  ERR = ZERO
  DO I = 1, L
     IF(MODE=='s') THEN
        CALL AV(M,IAP,JA,A,VK(:,I), VTEMP)
     ELSE IF(MODE=='d') THEN
        CALL DGEMV('N',M,N,ONE,A,M,VK(:,I),1,ZERO,VTEMP,1)
     END IF
     CALL DAXPY(M, -BK(I,I), UK(1,I), 1, VTEMP,1)
     ERR = ERR + DNRM2(M,VTEMP,1) / SQRT(DBLE(M))

     IF(MODE=='s') THEN
        CALL ATV(M,N,IAP,JA,A,UK(:,I), VTEMP)
     ELSE IF(MODE=='d') THEN
        CALL DGEMV('T',M,N,ONE,A,M,UK(:,I),1,ZERO,VTEMP,1)
     END IF
     CALL DAXPY(N, -BK(I,I), VK(1,I), 1, VTEMP,1)
     ERR = ERR + DNRM2(N,VTEMP,1) / SQRT(DBLE(N))
  ENDDO
  RETURN
END FUNCTION ERR
