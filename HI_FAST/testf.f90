subroutine blas_3dvsdot(N,A,B)
  implicit none
  include 'mkl_blas.fi'
  include 'mkl_vml.f90'
  integer :: N
  real (kind=4), intent(IN):: A(3,N),B(N)
  real (kind=4), intent(OUT):: B(N)
  B(:) =  sqrt(sdot(3, A(1:3,:),1, A(1:3,:), 1))
end subroutine blas_3dvsdot
