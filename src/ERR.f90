DOUBLE PRECISION FUNCTION ERR(start_row,mode,IAP,JA,A,M,N,K,L,BK,VK,UK,VTEMP,errmax)
  IMPLICIT NONE
  CHARACTER MODE
  INTEGER M,N,K,L,I
  INTEGER IAP(*),JA(*)
  INTEGER start_row(*)
  DOUBLE PRECISION ONE,ZERO,DNRM2,TMP,DDOT
  PARAMETER(ONE = 1.0D+0,ZERO = 0.0D+0)
  DOUBLE PRECISION A(*),BK(K,K),VK(N,K),UK(M,K),VTEMP(*),hoge,errmax
  errmax = ZERO
  ERR = ZERO
  hoge = ZERO
  DO I = 1, L
     TMP = ZERO

     IF(MODE=='s') THEN
        CALL AV(hoge,start_row,IAP,JA,A,VK(:,I), VTEMP)
     ELSE IF(MODE=='d') THEN
        CALL DGEMV('N',M,N,ONE,A,M,VK(:,I),1,ZERO,VTEMP,1)
     END IF
     CALL DAXPY(M, -BK(I,I), UK(1,I), 1, VTEMP,1)
     TMP = TMP + DDOT(M,VTEMP,1,VTEMP,1)

     IF(MODE=='s') THEN
        CALL ATV(hoge,start_row,N,IAP,JA,A,UK(:,I), VTEMP)
     ELSE IF(MODE=='d') THEN
        CALL DGEMV('T',M,N,ONE,A,M,UK(:,I),1,ZERO,VTEMP,1)
     END IF
     CALL DAXPY(N, -BK(I,I), VK(1,I), 1, VTEMP,1)
     TMP = TMP + DDOT(N,VTEMP,1,VTEMP,1)

     TMP=SQRT(TMP/2.0D+0)
     errmax = max(errmax,tmp)
     ERR = ERR + TMP / DBLE(L)
     print *,TMP
  ENDDO
  RETURN
END FUNCTION ERR
