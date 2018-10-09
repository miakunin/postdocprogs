MODULE NUMERICS

use LAKE_DATATYPES, only : ireals, iintegers
use NUMERIC_PARAMS, only : vector_length

contains
SUBROUTINE MATRIXSUM(a,b,c,k)
 implicit none
 
!MATRIXES: C=A+B
 real(kind=ireals), dimension (vector_length,2,2)::a,b,c
 integer(kind=iintegers) k,j,i

 do j=1,2
  do i=1,2
   c(k,i,j)=a(k,i,j)+b(k,i,j)
  enddo
 enddo
 
 END SUBROUTINE MATRIXSUM


 SUBROUTINE MATRIXMULT(a,b,c,k)
 implicit none

!MATRIXES: C=A*B      

 real(kind=ireals), dimension (vector_length,2,2)::a,b,c
 integer(kind=iintegers) k

 c(k,1,1)=a(k,1,1)*b(k,1,1)+a(k,1,2)*b(k,2,1)
 c(k,1,2)=a(k,1,1)*b(k,1,2)+a(k,1,2)*b(k,2,2)
 c(k,2,1)=a(k,2,1)*b(k,1,1)+a(k,2,2)*b(k,2,1)
 c(k,2,2)=a(k,2,1)*b(k,1,2)+a(k,2,2)*b(k,2,2)

 END SUBROUTINE MATRIXMULT


 SUBROUTINE MATRIXMULTVECTOR(a,g,f,k)
 implicit none

!MATRIX A, VECTORS g, f: Ag=f

 real(kind=ireals) a(vector_length,2,2),f(vector_length,2), &
 & g(vector_length,2)
 integer(kind=iintegers) k

 f(k,1)=a(k,1,1)*g(k,1)+a(k,1,2)*g(k,2)
 f(k,2)=a(k,2,1)*g(k,1)+a(k,2,2)*g(k,2)

 return
 END SUBROUTINE MATRIXMULTVECTOR


 SUBROUTINE VECTORSUM(a,b,c,k)
 implicit none

!VECTORS: C=A+B

 real(kind=ireals), dimension(vector_length,2)::a,b,c
 integer(kind=iintegers) k
 c(k,1)=a(k,1)+b(k,1)
 c(k,2)=a(k,2)+b(k,2)
 END SUBROUTINE VECTORSUM


 SUBROUTINE INVERSMATRIX(a,a1,k)
 implicit none 
 
!MATRIXES: A1=A*(-1)

 real(kind=ireals), dimension(vector_length,2,2)::a,a1
 integer(kind=iintegers) k
 
 a1(k,1,1)=a(k,2,2)/(a(k,1,1)*a(k,2,2)-a(k,1,2)*a(k,2,1)) 
 a1(k,1,2)=-a(k,1,2)/(a(k,1,1)*a(k,2,2)-a(k,1,2)*a(k,2,1))
 a1(k,2,1)=-a(k,2,1)/(a(k,1,1)*a(k,2,2)-a(k,1,2)*a(k,2,1))
 a1(k,2,2)=a(k,1,1)/(a(k,1,1)*a(k,2,2)-a(k,1,2)*a(k,2,1))  

 END SUBROUTINE INVERSMATRIX

 
 SUBROUTINE MATRIXPROGONKA(a,b,c,d,y,N)
 
!MATRIXPROGONKA solves the set of MATRIX three-point diference equations 
 
 implicit none

 real(kind=ireals), dimension(vector_length,2,2):: a,b,c,x3,x32,x31,x4,alpha
 real(kind=ireals), dimension(vector_length,2):: y,d,x2,beta,x21
 integer(kind=iintegers) N,i,j,k

 call INVERSMATRIX(c,x32,1)
 call MATRIXMULT(x32,b,x3,1)
 do j = 1, 2
   do i = 1, 2
     alpha(2,i,j) = x3(1,i,j)
   enddo
 enddo
 call MATRIXMULTVECTOR(x32,d,x2,1)
 do i = 1, 2
   beta(2,i) = x2(1,i)
 enddo

 do k = 3, N
  CALL MATRIXMULT(A,ALPHA,X3,k-1)
  X4(k-1,1:2,1:2) = - X3(k-1,1:2,1:2)
  CALL MATRIXSUM(C,X4,X31,k-1)
  CALL INVERSMATRIX(X31,X32,k-1)
  CALL MATRIXMULT(X32,B,X3,k-1)
  do j = 1, 2
    do i = 1, 2
      alpha(k,i,j) = X3(k-1,i,j)
    enddo
  enddo
  !call matrixmult(x3,x31,x33,k-1)
  !call matrixsum(-b,x33,x3,k-1)
  CALL MATRIXMULTVECTOR(A,BETA,X2,K-1)
  CALL VECTORSUM(D,X2,X21,K-1)
  CALL MATRIXMULTVECTOR(X32,X21,X2,K-1)
  do i = 1, 2
    beta(k,i) = X2(k-1,i)
  enddo
 enddo

 CALL MATRIXMULT(A,ALPHA,X3,N)
 X4(N,1:2,1:2) = - X3(N,1:2,1:2)
 CALL MATRIXSUM(C,X4,X31,N)
 CALL INVERSMATRIX(X31,X32,N)
 CALL MATRIXMULTVECTOR(A,BETA,X2,N)
 CALL VECTORSUM(D,X2,X21,N)
 CALL MATRIXMULTVECTOR(X32,X21,X2,N)
 do i = 1, 2
   Y(N,i) = X2(N,i)
 enddo

 do k = N-1, 1, -1
  CALL MATRIXMULTVECTOR(ALPHA,Y,X2,K+1)
  CALL VECTORSUM(X2,BETA,X21,K+1)
  Y(K,1) = X21(K+1,1)
  Y(K,2) = X21(K+1,2)
 enddo

 return

 END SUBROUTINE MATRIXPROGONKA


 real(kind=ireals) FUNCTION KRON(i,j)
 implicit none
 integer(kind=iintegers) i,j
 kron=0.
 if(i==j) kron=1.
 END FUNCTION 


 SUBROUTINE IND_STAB_FACT_DB (a,b,c,N,M,ind_stab,ind_bound)

 implicit none

 real(kind=ireals), dimension(1:vector_length):: a,b,c
 integer(kind=iintegers) M,i,N
 logical ind_stab, ind_bound 

 SAVE

 ind_stab=.true.
 if (ind_bound .eqv. .true.) then 
  if (abs(b(N))>=abs(c(N)).or.abs(a(M))>=abs(c(M))) then
   ind_stab=.false.
   RETURN
  endif
 endif
 do i=N+1,M-1
  if (abs(a(i))+abs(b(i))>=abs(c(i))) then
   ind_stab=.false.
   RETURN
  endif
 enddo
 
 END SUBROUTINE IND_STAB_FACT_DB


 LOGICAL FUNCTION CHECK_PROGONKA(N,a,b,c,d,y)

!Function CHECK_PROGONKA checks the accuracy
!of tridiagonal matrix system solution

 implicit none

 integer(kind=iintegers), intent(in):: N
 real(kind=ireals), intent(in):: a(1:N)
 real(kind=ireals), intent(in):: b(1:N)
 real(kind=ireals), intent(in):: c(1:N)
 real(kind=ireals), intent(in):: d(1:N)
 real(kind=ireals), intent(in):: y(1:N)

 real(kind=ireals), parameter:: del0 = 1.0d-13
 real(kind=ireals) del

 integer(kind=iintegers) i

 del = 0.d0
 del = max(c(1)*y(1)-b(1)*y(2)-d(1),del)
 del = max(c(N)*y(N)-a(N)*y(N-1)-d(N),del)
 do i = 2, N-1
   del = max(-a(i)*y(i-1)+c(i)*y(i)-b(i)*y(i+1)-d(i),del)
 enddo

 CHECK_PROGONKA = del < del0

 END FUNCTION CHECK_PROGONKA



 SUBROUTINE PROGONKA(a, b, c, f, y, K, N)
 implicit none
!FACTORIZATION METHOD FOR THE FOLLOWING SYSTEM OF LINEAR EQUATIONS:
!-a(i)*y(i-1)+c(i)*y(i)-b(i)*y(i+1)=f(i) i=K+1,N-1
! c(K)*y(K)-b(K)*y(K+1)=f(K)
!-a(N)*y(N-1)+c(N)*y(N)=f(N)
!

 integer(kind=iintegers), intent(in) :: K, N
 real(kind=ireals), intent(in) :: a(vector_length), b(vector_length), &
 & c(vector_length), f(vector_length) 
 real(kind=ireals), intent(out) :: y(vector_length)

 real(kind=ireals) :: alpha(vector_length+2), beta(vector_length+2) 
 integer(kind=iintegers) :: i

 SAVE
             
 alpha(K+1) = b(K)/c(K)
 beta(K+1) = f(K)/c(K)
 do i = K+2, N+1 
   alpha(i) = b(i-1)/(c(i-1)-a(i-1)*alpha(i-1))
   beta(i) = (f(i-1)+a(i-1)*beta(i-1))/ &
   & (c(i-1)-a(i-1)*alpha(i-1))
 end do
 y(N) = beta(N+1)
 do i = N-1, K, -1
   y(i) = alpha(i+1)*y(i+1)+beta(i+1)
 end do
  
 END SUBROUTINE PROGONKA


FUNCTION STEP(x)
! Heavyside (step) function of x
implicit none

real(kind=ireals) :: STEP

real(kind=ireals), intent(in) :: x

STEP = 0.5*(sign(1._ireals,x) + 1.)

END FUNCTION STEP


!>Function ACCUMSUM updates an accumulated mean
FUNCTION ACCUMM(n,sum_nm1,xn)

implicit none

real(kind=ireals), intent(in) :: sum_nm1 !> mean over n-1 values
real(kind=ireals), intent(in) :: xn !> n-th value of time series
integer(kind=iintegers), intent(in) :: n 

real(kind=ireals) :: ACCUMM

ACCUMM = ((n-1)*sum_nm1 + xn)/real(n)

END FUNCTION ACCUMM

REAL FUNCTION FindDet(matrix, n)
    IMPLICIT NONE
    REAL, DIMENSION(n,n) :: matrix
    INTEGER, INTENT(IN) :: n
    REAL :: m, temp
    INTEGER :: i, j, k, l
    LOGICAL :: DetExists = .TRUE.
    l = 1
    !Convert to upper triangular form
    DO k = 1, n-1
        IF (matrix(k,k) == 0) THEN
            DetExists = .FALSE.
            DO i = k+1, n
                IF (matrix(i,k) /= 0) THEN
                    DO j = 1, n
                        temp = matrix(i,j)
                        matrix(i,j)= matrix(k,j)
                        matrix(k,j) = temp
                    END DO
                    DetExists = .TRUE.
                    l=-l
                    EXIT
                ENDIF
            END DO
            IF (DetExists .EQV. .FALSE.) THEN
                FindDet = 0
                return
            END IF
        ENDIF
        DO j = k+1, n
            m = matrix(j,k)/matrix(k,k)
            DO i = k+1, n
                matrix(j,i) = matrix(j,i) - m*matrix(k,i)
            END DO
        END DO
    END DO
   
    !Calculate determinant by finding product of diagonal elements
    FindDet = l
    DO i = 1, n
        FindDet = FindDet * matrix(i,i)
    END DO
   
END FUNCTION FindDet

END MODULE NUMERICS
