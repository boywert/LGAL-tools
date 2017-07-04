subroutine blas_3dvsdot(N,A,B) bind (c,name='blas_3dvsdot')
  use iso_c_binding
  include 'mkl_blas.fi'
  include 'mkl_vml.f90'
  integer (c_int), intent(in), value :: N
  real (c_double), intent(IN):: A(N,3)
  real (c_double) :: AA(3,N)
  real (c_double), intent(OUT):: B(N)
  AA = transpose(A)
  B(:) =  ddot(3, AA(1:3,:),1, AA(1:3,:), 1)
  print *, B(1:10)
end subroutine blas_3dvsdot
