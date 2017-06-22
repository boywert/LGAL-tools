      PROGRAM   MAIN
      IMPLICIT NONE
      include 'mkl_blas.fi'
      include 'mkl_vml.f90'
      DOUBLE PRECISION ALPHA, BETA, PI
      INTEGER          M, K, N, I, J
      PARAMETER        (M=3, K=3, N=100000)
      DOUBLE PRECISION A(M,K), B(K,N), C(M,N), D(M,N)
      DOUBLE PRECISION theta 
      PI = 4*ATAN(1.0)
      theta = PI/2. 
      PRINT *, "This example computes real matrix C=alpha*A*B+beta*C"
      PRINT *, "using Intel(R) MKL function dgemm, where A, B, and C"
      PRINT *, "are matrices and alpha and beta are double precision "
      PRINT *, "scalars"
      PRINT *, ""

      PRINT *, "Initializing data for matrix multiplication C=A*B for "
      PRINT 10, " matrix A(",M," x",K, ") and matrix B(", K," x", N, ")"
 10   FORMAT(a,I5,a,I5,a,I5,a,I5,a)
      PRINT *, ""
      ALPHA = 1.0 
      BETA = 0.0

      PRINT *, "Intializing matrix data"
      PRINT *, ""
      
      PRINT *, "Rotation Matrix"
      A(:,:) = 0.0
      A(1,1) = 1.0
      A(2,2) = COS(theta)
      A(3,3) = COS(theta)
      A(2,3) = -1.*SIN(theta)
      A(3,2) = SIN(theta)

      B(1,:) = 1.0
      B(2,:) = 1.0
      B(3,:) = 1.0
      
      C(:,:) = 0.0

      PRINT *, "Computing matrix product using Intel(R) MKL DGEMM "
      PRINT *, "subroutine"
      CALL DGEMM('N','N',M,N,K,ALPHA,A,M,B,K,BETA,C,M)
      PRINT *, "Computations completed."
      PRINT *, ""
      D(1,:) = dsqrt(ddot(3, C(1:3,:),1, C(1:3,:), 1))
      call vdatan( N, C(3,:)/D(1,:), D(2,:) )
      call vdatan( N, C(2,:)/C(1,:), D(3,:) )
      D(2,:) = D(2,:)/PI*180
      D(3,:) = D(3,:)/PI*180
      PRINT *, "Top left corner of matrix A:"
      PRINT 20, ((A(I,J), J = 1,MIN(K,3)), I = 1,MIN(M,3))
      PRINT *, ""

      PRINT *, "Top left corner of matrix B:"
      PRINT 20, ((B(I,J),J = 1,MIN(N,3)), I = 1,MIN(K,3))
      PRINT *, ""

 20   FORMAT(3(F12.0,1x))

      PRINT *, "Top left corner of matrix C:"
      PRINT 30, ((C(I,J), J = 1,MIN(N,3)), I = 1,MIN(M,3))
      PRINT *, ""


      PRINT *, "Top left corner of matrix D:"
      PRINT 30, ((D(I,J), J = 1,MIN(N,3)), I = 1,MIN(M,3))
      PRINT *, ""

 30   FORMAT(3(ES12.4,1x))
      

      PRINT *, "Example completed."
      STOP 

      END
