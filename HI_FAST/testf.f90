subroutine blas_3dvsdot(N,A,B) bind (c,name='blas_3dvsdot')
  use iso_c_binding
  include 'mkl_blas.fi'
  include 'mkl_vml.f90'
  integer (c_int), intent(in), value :: N
  real (c_double), intent(IN):: A(3,N)
  real (c_double) :: AA(3,N)
  real (c_double), intent(OUT):: B(3,N)
  call vdabs( 3, A(1:3,1), B(1) )
  call vdabs( 3, A(1:3,2), B(2) ) 

  B(1,:) = vdsqrt(A(1,1:N)*A(1,1:N)+A(2,1:N)*A(2,1:N)+A(3,1:N)*A(3,1:N))

end subroutine blas_3dvsdot
