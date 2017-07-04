subroutine blas_3dvsdot(N,A,B) bind (c,name='blas_3dvsdot')
  use iso_c_binding
  include 'mkl_blas.fi'
  include 'mkl_vml.f90'
  integer (c_int), intent(in), value :: N
  real (c_double), intent(IN):: A(N,3)
  real (c_double) :: AA(N,3)
  real (c_double), intent(OUT):: B(N)
  call vdabs( N, A(:,1:3), B(:) ) 
  print *, A(1:10,:)
  print *, B(1:10,:)
end subroutine blas_3dvsdot
