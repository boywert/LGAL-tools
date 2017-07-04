subroutine blas_3dvsdot(N,A,B) bind (c,name='blas_3dvsdot')
  use iso_c_binding
  include 'mkl_blas.fi'
  include 'mkl_vml.f90'
  integer :: i 
  integer (c_int), intent(in), value :: N
  real (c_double), intent(IN):: A(3,N)
  real (c_double) :: AA(3,N)
  real (c_double), intent(OUT):: B(3,N)
  do i = 1,100
     B(1,:) = dsqrt(A(1,1:N)*A(1,1:N)+A(2,1:N)*A(2,1:N)+A(3,1:N)*A(3,1:N))
     B(2,:) = acos(A(3,:)/B(1,:))
     B(3,:) = atan(A(2,:)/A(1,:))
  end do
end subroutine blas_3dvsdot
