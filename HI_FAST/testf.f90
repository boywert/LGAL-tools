module test
contains

! subroutine make_sphere(N,boxsize,A,B) bind (c,name='make_sphere')
!   use iso_c_binding
!   !implicit none
!   integer (c_int), intent(in), value :: N
!   real (c_float), intent(IN) :: boxsize
!   real (c_float), intent(IN):: A(3,N)
!   !real (c_float), allocatable :: AC(:,:)
!   real (c_float), intent(OUT):: B(3,N)
!   !integer  :: i,j,k,l,index

!   !allocate(AC(3,N))
!   ! do all 8 quadrants
!   ! do i=1,2
!   !    do j=1,2
!   !       do k=1,2
!   !          index = (i-1)*2*2 + (j-1)*2 + k - 1
!   !          do l=1,N
!   !             AC(1:3,index*N+l) = A(1:3,l) - (/ (i-1), (j-1), (k-1) /)*boxsize
!   !          end do
!   !       end do
!   !    end do
!   ! end do

!   !AC = A
!   B(1,:) = sqrt(A(1,:)*A(1,:)+A(2,:)*A(2,:)+A(3,:)*A(3,:))
!   B(2,:) = acos(A(3,:)/B(1,:))
!   B(3,:) = atan(A(2,:)/A(1,:))
!   !deallocate(AC)
!   print *, "ready"
! end subroutine make_sphere

subroutine make_sphere(N,boxsize,A1,A2,B1,B2) bind (c,name='make_sphere')
  use iso_c_binding
  integer (c_int), intent(in), value :: N
  real (c_float), intent(IN), value :: boxsize
  real (c_float), intent(IN):: A1(3,N),A2(3,N)
  real (c_float), intent(OUT):: B1(3,8*N)
  real (c_float), intent(OUT):: B2(3,8*N)
  real (c_float), allocatable :: AC(:,:)
  real (c_float) :: PI 
  integer :: i,j,k,l
  PI = 4.0*atan(1.0)
  print *, N,boxsize

  allocate(AC(3,N))
  do i=1,2
     do j=1,2
        do k=1,2
           index = (i-1)*2*2 + (j-1)*2 + k - 1
           AC(1,:) = A1(1,:) - (i-1)*boxsize
           AC(2,:) = A1(2,:) - (j-1)*boxsize
           AC(3,:) = A1(3,:) - (k-1)*boxsize
           B1(1,index*N+1:index*N+N) = sqrt(AC(1,1:N)*AC(1,1:N) + AC(2,1:N)*AC(2,1:N) + AC(3,1:N)*AC(3,1:N))
           B1(2,index*N+1:index*N+N) = acos(max(-1.0,min(1.0,AC(3,1:N)/B1(1,index*N+1:index*N+N))))
           B1(3,index*N+1:index*N+N) = mod((atan2(AC(2,1:N),AC(1,1:N)) + 2*PI),2*PI)
           B2(1,index*N+1:index*N+N) = (AC(1,1:N)*A2(1,1:N) + AC(2,1:N)*A2(2,1:N) + AC(3,1:N)*A2(3,1:N))/B1(1,index*N+1:index*N+N)
           B2(2,index*N+1:index*N+N) = A2(1,1:N)*cos(B1(3,index*N+1:index*N+N))*cos(B1(2,index*N+1:index*N+N))*-1.0 + A2(2,1:N)*sin(B1(3,index*N+1:index*N+N))*cos(B1(2,index*N+1:index*N+N)) - A2(3,1:N)*sin(B1(2,index*N+1:index*N+N))
           
           !B2(2,index*N+1:index*N+N) = (A2(1,1:N)*cos(B1(2,index*N+1:index*N+N))*cos(B1(3,index*N+1:index*N+N)) + AC(1:N))*sin(B1(2,index*N+1:index*N+N))*cos(B1(3,index*N+1:index*N+N)) - AC(3,1:N)*sin(B1(3,index*N+1:index*N+N)))/B1(1,index*N+1:index*N+N)
           !B2(3,index*N+1:index*N+N) = (-1*sin(B1(3,index*N+1:index*N+N))*A2(1,:) + cos(B1(3,index*N+1:index*N+N))*A2(2,:))/sqrt(A2(1,:)*A2(1,:) + A2(2,:)*A2(2,:) + A2(3,:)*A2(3,:))
        end do
     end do
  end do
  deallocate(AC)
end subroutine make_sphere

subroutine cart2sphere1(N,A,B) bind (c,name='cart2sphere1')
  use iso_c_binding
  integer (c_int), intent(in), value :: N
  real (c_float), intent(IN):: A(3,N)
  real (c_float), intent(OUT):: B(3,N)
  B(1,:) = sqrt(A(1,1:N)*A(1,1:N)+A(2,1:N)*A(2,1:N)+A(3,1:N)*A(3,1:N))
  B(2,:) = acos(A(3,:)/B(1,:))
  B(3,:) = atan(A(2,:)/A(1,:))
end subroutine cart2sphere1

  ! subroutine cart2sphere2(N,A,B) bind (c,name='cart2sphere2')
  !   use iso_c_binding
  !   include 'mkl_blas.fi'
  !   include 'mkl_vml.f90'
  !   integer (c_int), intent(in), value :: N
  !   real (c_float), intent(IN):: A(3,N)
  !   real (c_float), intent(OUT):: B(3,N)
  !   call vssqrt(N,A(1,1:N)*A(1,1:N)+A(2,1:N)*A(2,1:N)+A(3,1:N)*A(3,1:N),B(1,:))
  !   call vsacos(N,A(3,:)/B(1,:),B(2,:))
  !   call vsatan(N,A(2,:)/A(1,:),B(3,:))
  ! end subroutine cart2sphere2

end module test
