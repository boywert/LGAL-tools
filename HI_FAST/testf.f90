subroutine blas_3dvsdot(N,A,B) bind (c,name='blas_3dvsdot')
  use iso_c_binding
  include 'mkl_blas.fi'
  include 'mkl_vml.f90'
  integer (c_int), intent(in), value :: N
  real (c_double), intent(IN):: A(3,N)
  real (c_double) :: AA(3,N)
  real (c_double), intent(OUT):: B(N)
  call vdabs( 1, A(1:3,1), B(1) )
  call vdabs( 1, A(1:3,2), B(2) ) 
  print *, A(1:3,1:2)
  print *, ""
  print *, B(1:2)
end subroutine blas_3dvsdot
