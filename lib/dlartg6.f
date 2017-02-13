      SUBROUTINE DLARTG6( K, L, M, N, CS, SN )
      IMPLICIT NONE
*
*     .. Scalar Arguments ..
      DOUBLE PRECISION   K, L, M, N
      DOUBLE PRECISION   CS, R, SN
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO
      PARAMETER          ( ZERO = 0.0D0 )
*     ..
*     .. Local Scalars ..
      DOUBLE PRECISION   F, G
*     ..
*     .. External Functions ..
      DOUBLE PRECISION   DFMA0
      EXTERNAL           DFMA0
*     ..
*     .. Executable Statements ..
*
      IF (N .GE. K) THEN

         F = DFMA0(K/N,L,M)
         G = DFMA0(M/N,L,-K)

      ELSE

         F = DFMA0(M,N/K,L)
         G = DFMA0(M,L/K,-N)

      ENDIF
      
      IF( G .GE. ZERO ) THEN
         CALL DLARTG(F,G,CS,SN,R)
      ELSE
         CALL DLARTG(F,-G,CS,SN,R)
         SN = -SN
      END IF
      RETURN
*     
*     End of DLARTG6
*     
      END
