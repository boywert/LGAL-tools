subroutine blas_3dvsdot(N,A,B) bind (c,name='blas_3dvsdot')
  use iso_c_binding
  include 'mkl_blas.fi'
  include 'mkl_vml.f90'
  integer (c_int), intent(in), value :: N
  real (c_float), intent(IN):: A(3,N)
  real (c_float), intent(OUT):: B(N)
  B(:) =  sqrt(sdot(3, A(1:3,:),1, A(1:3,:), 1))
  print *, A
end subroutine blas_3dvsdot
