SUBROUTINE DGEBRDG_LP1(N,A,LDA,Q,P)
  IMPLICIT NONE
  INTEGER N,LDA,i,j
  DOUBLE PRECISION  ONE, ZERO,cs,sn
  PARAMETER ( ONE = 1.0D+0, ZERO=0.0D+0 )
  DOUBLE PRECISION A(LDA,*),Q(LDA,*),P(LDA,*)
  EXTERNAL dlartg,drot

  IF ( n .LE. 2 ) THEN
     RETURN
  END IF
  DO i = N, 3, -1
     DO j = 1, i -2
        ! WRITE(*,*) "i=",i,"j=",j,"step1"
        ! DO k = 1, n
        !    WRITE (*,*) A(k,1:n)
        ! END DO
        CALL DLARTG( a(j+1,i), a(j,i), CS, SN,a(j+1,i))
        a(j,i)=0
        IF (( cs .NE. 1) .OR. (sn .NE. 0)  ) THEN
           CALL drot (i-j+1-1,A(j+1,j),LDA,A(j,j),LDA,CS,SN)
           CALL drot (LDA,Q(1,j+1),1,Q(1,j),1,CS,SN)
        END IF
        ! WRITE(*,*) "step2"
        ! DO k = 1, n
        !    WRITE (*,*) A(k,1:n)
        ! END DO
        CALL DLARTG( a(j+1,j+1),a(j+1,j),CS, SN,a(j+1,j+1))
        a(j+1,j)=0
        IF (( cs .NE. 1) .OR. (sn .NE. 0)  ) THEN
           CALL drot (j+1-1,A(1,j+1),1,A(1,j),1,CS,SN)
           CALL drot (LDA,P(j+1,1),LDA,P(j,1),LDA,CS,SN)
        END IF

        ! WRITE(*,*) "step3"
        ! DO k = 1, n
        !    WRITE (*,*) A(k,1:n)
        ! END DO
        ! WRITE(*,*)
     END DO
  END DO
  RETURN
END SUBROUTINE DGEBRDG_LP1
