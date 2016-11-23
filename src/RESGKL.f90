SUBROUTINE RESGKL(J,MODE,LS,IAP,JA,A,M,N,K,L,BK,VK,UK,VPLUS,INFO,SELEK,WORK,LWORK)
  IMPLICIT NONE

  INTEGER IAP(*), JA(*), LWORK
  DOUBLE PRECISION A(*), WORK(*)
  CHARACTER MODE,LS

  INTEGER M,N,K,L,I,J,IINFO,INFO,SELEK
  DOUBLE PRECISION ONE,ZERO,TWO,MINUSONE
  PARAMETER(ONE=1.0D+0,ZERO=0.0D+0,TWO=2.0D+0,MINUSONE=-1.0D+0)
  DOUBLE PRECISION CDUMMY(1,1),BETA
  DOUBLE PRECISION BK(K,K),CPBK(K,K),VK(N,K),UK(M,K),VPLUS(N),BD(K),BE(K)
  DOUBLE PRECISION VM(N,L),UM(M,L),WORK2(K*8)
  DOUBLE PRECISION Q(K,K),P(K,K)
  DOUBLE PRECISION DNRM2,DLAMCH

  VK(1:N,J+1) = VPLUS(1:N)
  DO WHILE(J < K)
     IF(MODE=='s') THEN
        CALL av(M,IAP,JA,A, VK(1:N,J+1), UM(1:M,1))
     ELSE IF(MODE=='d') THEN
        CALL DGEMV('N',M,N,ONE,A,M,VK(1:N,J+1),1,ZERO,UM(1:M,1),1)
     END IF
     CALL CGS2(UM(1:M,1),UK,M,J,WORK)
     BK(J+1,J+1) = DNRM2(M,UM(1:M,1),1)
     UK(1:M,J+1) = UM(1:M,1) / BK(J+1,J+1)
     IF(MODE=='s') THEN
        CALL atv(M,N,IAP,JA,A, UK(1:M,J+1), VPLUS)
     ELSE IF(MODE=='d') THEN
        CALL DGEMV('T',M,N,ONE,A,M,UK(1:M,J+1),1,ZERO,VPLUS,1)
     END IF
     CALL CGS2(VPLUS,VK,N,J+1,WORK)
     IF(J < K - 1) THEN
        BK(J+1,J+2) = DNRM2(N,VPLUS,1)
        VK(1:N,J+2) = VPLUS(1:N) / BK(J+1,J+2)
     ELSE
        BETA = DNRM2(N,VPLUS,1)
        VPLUS(1:N) = VPLUS(1:N) / BETA
        !call dscal(N,1/BETA,VPLUS,1)
     END IF
     J = J+1
  END DO
  J = L

  Q = ZERO
  DO I = 1,K
     Q(I,I) = ONE
  END DO
  P = ZERO
  DO I = 1,K
     P(I,I) = ONE
  END DO

  IF (SELEK == 1) THEN
     ! with lapack 1.0 QR
     !  WRITE(*,*) "GIVENS回転(DGEBRDG_LP1)+QR法(DBDSQRU)+両側(DGEMM)"
     CALL DGEBRDG_4_BISIDE(L+1,BK,K,Q,P)

     DO I = 1,K
        BD(I) = BK(I,I)
     END DO
     DO I = 1,K-1
        BE(I) = BK(I,I+1)
     END DO
     CALL DBDSQRU( 'U',K,K,K,0,BD,BE,P,K,Q,K,CDUMMY,K,WORK,IINFO )
     IF(LS == "s") then
     DO I = 1, K/2
        CALL DSWAP(K,Q(1,K+1-I),1,Q(1,I),1) !Q(1:K,I) = Q(1:K,K+1-I)
     END DO
     DO I = 1, K/2
        CALL DSWAP(K,P(K+1-I,1),K,P(I,1),K) !P(I,1:K) = P(K+1-I,1:K)
     END DO
     end if
     CALL DGEMM('N','N',M,L,K,ONE,UK,M,Q,K,ZERO,UM,M)
     CALL DCOPY(M*L,UM,1,UK,1)
     CALL DGEMM('N','T',N,L,K,ONE,VK,N,P,K,ZERO,VM,N)
     CALL DCOPY(N*L,VM,1,VK,1)

     BK = ZERO
     IF(LS == "s") then
     DO I =1 ,L
        BK(I,I) = BD(K-I+1)
     END DO
     else
     DO I =1 ,L
        BK(I,I) = BD(I)
     END DO
     end if
     CALL DAXPY(L,BETA,Q(K,1),K,BK(1,L+1),1)

  ELSE IF ( SELEK==2 ) THEN
     ! with lapack 1.0 QR
     !  WRITE(*,*) "GIVENS回転(DGEBRDG_LP1)+QR法(DBDSQRU)+片側(DORMQR)"
     CALL DCOPY(K*K,BK,1,CPBK,1)
     CALL DGEBRDG_LP1(L+1,BK,K,Q,P)

     DO I = 1,K
        BD(I) = BK(I,I)
     END DO

     DO I = 1,K-1
        BE(I) = BK(I,I+1)
     END DO

     CALL DBDSQRU('U',K,K,0,0,BD,BE,P,K,Q,K,CDUMMY,K,WORK,IINFO )
     IF(LS == "s") then
     DO I = 1, K/2
        CALL DSWAP(K,P(K+1-I,1),K,P(I,1),K) !P(I,1:K) = P(K+1-I,1:K)
     END DO
     end if
     P = TRANSPOSE(P)
     CALL DGEQRF(K,L,P,K,BD,WORK,LWORK,IINFO)
     CALL DORMQR('R','N',K,K,L,P,K,BD,CPBK,K,WORK,LWORK,IINFO )
     CALL DORMQR('R','N',N,K,L,P,K,BD,VK,N,WORK,LWORK,IINFO )
     CALL DGEQRF(K,L,CPBK,K,BD,WORK,LWORK,IINFO)
     CALL DORMQR('R','N',M,K,L,CPBK,K,BD,UK,M,WORK,LWORK,IINFO )

     BK = ZERO
     DO I = 1,L
        CALL DCOPY(I,CPBK(1,I),1,BK(1,I),1)
     END DO

     BE(1:K) = ZERO
     BE(K) = ONE
     CALL DORMQR('R','N',1,K,L,CPBK,K,BD,BE,1,WORK,LWORK,IINFO ) 
     CALL DAXPY(L,BETA,BE,1,BK(1,L+1),1)

  ELSE IF(SELEK==3) THEN
     !WRITE(*,*) "GIVENS回転(DGEBRDG_LP1)+OQDS1法(DOQDS1)+両側(DGEMM)"
     CALL DGEBRDG_4_BISIDE(L+1,BK,K,Q,P)
     P = TRANSPOSE(P)
     Q = TRANSPOSE(Q)
     DO I = 1,K
        BD(I) = BK(I,I)
     END DO
     DO I = 1,K-1
        BE(I) = BK(I,I+1)
     END DO
     CALL DOQDS1('L',K,BD,BE,P,K,Q,K,WORK,WORK2,INFO)
     IF(LS == "s") then
     DO I = 1, K/2
        CALL DSWAP(K,P(1,K+1-I),1,P(1,I),1) !P(1:K,I) = P(1:K,K+1-I)
     END DO
     DO I = 1, K/2
        CALL DSWAP(K,Q(K+1-I,1),K,Q(I,1),K) !Q(I,1:K) = Q(K+1-I,1:K)
     END DO
     end if
     CALL DGEMM('N','T',M,L,K,ONE,UK,M,Q,K,ZERO,UM,M)
     CALL DCOPY(M*L,UM,1,UK,1)
     CALL DGEMM('N','N',N,L,K,ONE,VK,N,P,K,ZERO,VM,N)
     CALL DCOPY(N*L,VM,1,VK,1)

     BK = ZERO
     IF(LS == "s") then
     DO I =1 ,L
        BK(I,I) = BD(K-I+1)
     END DO
     else
     DO I =1 ,L
        BK(I,I) = BD(I)
     END DO
     end if
     CALL DAXPY(L,BETA,Q(1,K),1,BK(1,L+1),1)
  
  ELSE IF(SELEK==4) THEN
     !  WRITE(*,*) "GIVENS回転(DGEBRDG_LP1)+OQDS2法(DOQDS2)+片側(DORMQR)"

     CALL DCOPY(K*K,BK,1,CPBK,1)
     CALL DGEBRDG_LP1(L+1,BK,K,Q,P)

     DO I = 1,K
        BD(I) = BK(I,I)
     END DO

     DO I = 1,K-1
        BE(I) = BK(I,I+1)
     END DO

     CALL DOQDS3('U',K,BD,BE,P,K,WORK,IINFO )
     IF(LS == "s") then
     DO I = 1, K/2
        CALL DSWAP(K,P(K+1-I,1),K,P(I,1),K) !P(I,1:K) = P(K+1-I,1:K)
     END DO
     end if
     P = TRANSPOSE(P)
     CALL DGEQRF(K,L,P,K,BD,WORK,LWORK,IINFO)
     CALL DORMQR('R','N',K,K,L,P,K,BD,CPBK,K,WORK,LWORK,IINFO )
     CALL DORMQR('R','N',N,K,L,P,K,BD,VK,N,WORK,LWORK,IINFO )
     CALL DGEQRF(K,L,CPBK,K,BD,WORK,LWORK,IINFO)
     CALL DORMQR('R','N',M,K,L,CPBK,K,BD,UK,M,WORK,LWORK,IINFO )

     BK = ZERO
     DO I = 1,L
        CALL DCOPY(I,CPBK(1,I),1,BK(1,I),1)
     END DO

     BE(1:K) = ZERO
     BE(K) = ONE
     CALL DORMQR('R','N',1,K,L,CPBK,K,BD,BE,1,WORK,LWORK,IINFO ) 
     CALL DAXPY(L,BETA,BE,1,BK(1,L+1),1)

  ELSE IF (SELEK==5) THEN
     !  WRITE(*,*) "ヤコビ法(DGESVJ)+両側(DGEMM)"

     WORK(1) = ONE+TWO*DLAMCH('E')
     CALL DGESVJ('U','C','V',K,K,BK,K,BD,M,P,K,WORK,LWORK,IINFO)
     IF ( NINT(WORK(2)) .LT. L ) THEN
        INFO = 2
        RETURN
     END IF

     CALL DGEMM('N','N',N,L,K,ONE,VK,N,P,K,ZERO,VM,N)
     CALL DCOPY(N*L,VM,1,VK,1)
     CALL DGEMM('N','N',M,L,K,ONE,UK,M,BK,K,ZERO,UM,M)
     CALL DCOPY(M*L,UM,1,UK,1)

     CALL DCOPY(L,BK(K,1),K,BE,1)

     BK = ZERO
     DO I =1 ,L
        BK(I,I) = BD(I)
     END DO
     CALL DAXPY(L,BETA,BE,1,BK(1,L+1),1)

  ELSE IF ( SELEK==6) THEN
     !  WRITE(*,*) "ヤコビ法(DGESVJ)+片側(DORMQR)"
     CALL DCOPY(K*K,BK,1,CPBK,1)

     WORK(1) = ONE+TWO*DLAMCH('E')
     CALL DGESVJ('U','C','V',K,K,BK,K,BD,M,P,K,WORK,LWORK,IINFO)
     CALL DGEQRF(K,L,P,K,BD,WORK,LWORK,IINFO)
     CALL DORMQR('R','N',K,K,L,P,K,BD,CPBK,K,WORK,LWORK,IINFO )
     CALL DORMQR('R','N',N,K,L,P,K,BD,VK,N,WORK,LWORK,IINFO )
     CALL DGEQRF(K,L,CPBK,K,BD,WORK,LWORK,IINFO)
     CALL DORMQR('R','N',M,K,L,CPBK,K,BD,UK,M,WORK,LWORK,IINFO )

     BK = ZERO
     DO I = 1,L
        CALL DCOPY(I,CPBK(1,I),1,BK(1,I),1)
     END DO

     BE(1:K) = ZERO
     BE(K) = ONE
     CALL DORMQR('R','N',1,K,L,CPBK,K,BD,BE,1,WORK,LWORK,IINFO )
     CALL DAXPY(L,BETA,BE,1,BK(1,L+1),1)

  ELSE IF (SELEK==7) THEN
     !  WRITE(*,*) "GIVENSrotation+ヤコビ法(DGESVJ)+両側(DGEMM)"
     ! 作成中
     CALL DGEBRDG_4_BISIDE(L+1,BK,K,Q,P)
     P = TRANSPOSE(P)
     WORK(1) = ONE+TWO*DLAMCH('E')
     CALL DGESVJ('U','C','A',K,K,BK,K,BD,K,P,K,WORK,LWORK,IINFO)

     CALL DGEMM('N','N',N,L,K,ONE,VK,N,P,K,ZERO,VM,N)
     CALL DCOPY(N*L,VM,1,VK,1)

     CALL DGEMM('N','N',K,L,K,ONE,Q,K,BK,K,ZERO,P,K)
     CALL DGEMM('N','N',M,L,K,ONE,UK,M,P,K,ZERO,UM,M)
     CALL DCOPY(M*L,UM,1,UK,1)

     BK = ZERO
     DO I =1 ,L
        BK(I,I) = BD(I)
     END DO
     CALL DAXPY(L,BETA,P(K,1),K,BK(1,L+1),1)

  END IF
  RETURN
END SUBROUTINE RESGKL

subroutine av (M, IAP, JA ,A, P, AP)
      
  IMPLICIT NONE
  include 'omp_lib.h'
  integer M
  integer IAP(*),JA(*)
  double precision A(*),P(*),AP(*),zero
  parameter (zero=0.0d0)
  integer i,j
  
  !$OMP PARALLEL DO PRIVATE(j)
  do i=1,M
     AP(i)=zero
     do j=IAP(i),IAP(i+1)-1
        AP(i)=AP(i)+(A(j))*P(JA(j))
     enddo
  enddo
  !$OMP END PARALLEL DO
  
  return
end subroutine av

subroutine atv (M, N, IAP, JA, A, Q, AQ)
  
  IMPLICIT NONE
  include 'omp_lib.h'
  integer M,N
  integer IAP(*),JA(*)
  double precision A(*),Q(*),AQ(N),zero
  parameter (zero=0.0d0)
  
  integer i,j
  
  do i=1,N
     AQ(i)=zero
  enddo
  
  !$OMP PARALLEL DO PRIVATE(j) REDUCTION(+:AQ)
  do i=1,M
     do j=IAP(i),IAP(i+1)-1
        AQ(JA(j))=AQ(JA(j))+(A(j))*Q(i)
     enddo
  enddo
  !$OMP END PARALLEL DO
  
  return
end subroutine atv

!SUBROUTINE CGS3(vntemp,V,n,j,WORK)
!  IMPLICIT NONE
!  include 'omp_lib.h'
!  INTEGER n,j,i
!  DOUBLE PRECISION  ONE, ZERO, MINUSONE
!  PARAMETER ( ONE = 1.0D+0, ZERO=0.0D+0,MINUSONE = -1.0D+0 )
!  DOUBLE PRECISION V(n,j),vntemp(n),tmp,WORK(n)
!  do i=1,N
!     WORK(i)=zero
!  enddo
!  !$OMP PARALLEL DO PRIVATE(tmp) REDUCTION(+:WORK)
!  do i = 1,j
!     CALL DDOT(n,V(1,i),1,vntemp,1,tmp)
!     WORK =WORK+ tmp * V(1:n,i)
!  end do
!  !$OMP END PARALLEL DO
!  vntemp = vntemp - WORK
!  do i=1,N
!     WORK(i)=zero
!  enddo
!  !$OMP PARALLEL DO PRIVATE(tmp) REDUCTION(+:WORK)
!  do i = 1,j
!     CALL DDOT(n,V(1,i),1,vntemp,1,tmp)
!     WORK =WORK+ tmp * V(1:n,i)
!  end do
!  !$OMP END PARALLEL DO
!  vntemp = vntemp - WORK
!
!END SUBROUTINE CGS3
