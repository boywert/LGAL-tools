subroutine cart2sphere(N,A,B) bind (c,name='cart2sphere')
  use iso_c_binding
end subroutine cart2sphere

subroutine cart2sphere1(N,A,B) bind (c,name='cart2sphere1')
  use iso_c_binding
  include 'mkl_blas.fi'
  include 'mkl_vml.f90'
  integer (c_int), intent(in), value :: N
  real (c_float), intent(IN):: A(3,N)
  real (c_float), intent(OUT):: B(3,N)
  B(1,:) = sqrt(A(1,1:N)*A(1,1:N)+A(2,1:N)*A(2,1:N)+A(3,1:N)*A(3,1:N))
  B(2,:) = acos(A(3,:)/B(1,:))
  B(3,:) = atan(A(2,:)/A(1,:))
end subroutine cart2sphere1

subroutine cart2sphere2(N,A,B) bind (c,name='cart2sphere2')
  use iso_c_binding
  include 'mkl_blas.fi'
  include 'mkl_vml.f90'
  integer (c_int), intent(in), value :: N
  real (c_float), intent(IN):: A(3,N)
  real (c_float), intent(OUT):: B(3,N)
  call vssqrt(N,A(1,1:N)*A(1,1:N)+A(2,1:N)*A(2,1:N)+A(3,1:N)*A(3,1:N),B(1,:))
  call vsacos(N,A(3,:)/B(1,:),B(2,:))
  call vsatan(N,A(2,:)/A(1,:),B(3,:))
end subroutine cart2sphere2

