      SUBROUTINE DLARTG4( F, G, CS )
      IMPLICIT NONE
*
*     .. Scalar Arguments ..
      DOUBLE PRECISION   CS, F, G, R
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO
      PARAMETER          ( ZERO = 0.0D0 )
      DOUBLE PRECISION   ONE
      PARAMETER          ( ONE = 1.0D0 )
      DOUBLE PRECISION   TWO
      PARAMETER          ( TWO = 2.0D0 )
*     ..
*     .. Local Scalars ..
*     LOGICAL            FIRST
c      INTEGER            COUNT
c      INTEGER            I
      DOUBLE PRECISION   EPS, F1, G1, SAFMIN, SAFMN2, SAFMX2, SCALE
*     ..
*     .. External Functions ..
      DOUBLE PRECISION   DLAMCH
      EXTERNAL           DLAMCH
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, INT, LOG, MAX, SQRT
*     ..
*     .. Save statement ..
*     SAVE               FIRST, SAFMX2, SAFMIN, SAFMN2
*     ..
*     .. Data statements ..
*     DATA               FIRST / .TRUE. /
*     ..
*     .. Executable Statements ..
*
*     IF( FIRST ) THEN
         SAFMIN = DLAMCH( 'S' )
         EPS = DLAMCH( 'E' )
         SAFMN2 = DLAMCH( 'B' )**INT( LOG( SAFMIN / EPS ) /
     $            LOG( DLAMCH( 'B' ) ) / TWO )
         SAFMX2 = ONE / SAFMN2
*        FIRST = .FALSE.
*     END IF
      IF( G.EQ.ZERO ) THEN
         CS = ONE
c         SN = ZERO
c         R = F
      ELSE IF( F.EQ.ZERO ) THEN
         CS = ZERO
c         SN = ONE
c         R = G
      ELSE
         F1 = F
         G1 = G
         SCALE = MAX( F1, G1 )
         IF( SCALE.GE.SAFMX2 ) THEN
c            COUNT = 0
   10       CONTINUE
c            COUNT = COUNT + 1
            F1 = F1*SAFMN2
            G1 = G1*SAFMN2
            SCALE = MAX( F1, G1 )
            IF( SCALE.GE.SAFMX2 )
     $         GO TO 10
            R = SQRT( F1**2+G1**2 )
            CS = F1 / R
c            SN = G1 / R
c            DO 20 I = 1, COUNT
c               R = R*SAFMX2
c   20       CONTINUE
         ELSE IF( SCALE.LE.SAFMN2 ) THEN
c            COUNT = 0
   30       CONTINUE
c            COUNT = COUNT + 1
            F1 = F1*SAFMX2
            G1 = G1*SAFMX2
            SCALE = MAX( F1, G1 )
            IF( SCALE.LE.SAFMN2 )
     $         GO TO 30
            R = SQRT( F1**2+G1**2 )
            CS = F1 / R
c            SN = G1 / R
c            DO 40 I = 1, COUNT
c               R = R*SAFMN2
c   40       CONTINUE
         ELSE
            R = SQRT( F1**2+G1**2 )
            CS = F1 / R
c            SN = G1 / R
         END IF
c         IF( ABS( F ).GT.ABS( G ) .AND. CS.LT.ZERO ) THEN
c            CS = -CS
c            SN = -SN
c            R = -R
c         END IF
      END IF
      RETURN
*
*     End of DLARTG4
*
      END
