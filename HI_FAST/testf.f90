! module test
! contains

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

subroutine make_sphere(N,boxsize,A,B) bind (c,name='make_sphere')
  use iso_c_binding
  integer (c_int), intent(in), value :: N
  real (c_float), intent(IN), value :: boxsize
  real (c_float), intent(IN):: A(3,N)
  real (c_float), intent(OUT):: B(N,3)
  real (c_float), allocatable :: AC(:,:)
  integer :: i,j,k,l
  print *, N,boxsize

  allocate(AC(N,3))
  do i=1,2
     do j=1,2
        do k=1,2
           index = (i-1)*2*2 + (j-1)*2 + k - 1
           if ((i .eq. 1) .and. (j .eq. 1) .and. (k .eq. 1)) then
              call cart2sphere1(N,AC,B(1:3,index*N+1:index*N+N))
           else
              AC(1,:) = A(1,:) - (i-1)*boxsize
              AC(2,:) = A(2,:) - (j-1)*boxsize
              AC(3,:) = A(3,:) - (k-1)*boxsize
              call cart2sphere1(N,AC,B(1:3,index*N+1:index*N+N))
           end if
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

! end module test                 
