INCLUDE 'test_mod.f03'
      program test
      USE test_mod
      implicit none
      real(kind=real64),dimension(:,:),allocatable::A,B,C
!
      allocate(A(3,2),B(3,6),C(2,6))
!
      A(1,1) = float(1)
      A(2,1) = float(4)
      A(3,1) = float(2)
      A(1,2) = float(5)
      A(2,2) = float(1)
      A(3,2) = float(0)
      A(1,3) = float(3)
      A(2,3) = float(3)
      A(3,3) = float(6)
!
      B(1,1) = float(2)
      B(2,1) = float(1)
      B(3,1) = float(1)
      B(1,2) = float(3)
      B(2,2) = float(1)
      B(3,2) = float(4)
      B(1,3) = float(4)
      B(2,3) = float(2)
      B(3,3) = float(2)
      B(1,4) = float(1)
      B(2,4) = float(2)
      B(3,4) = float(3)
      B(1,5) = float(3)
      B(2,5) = float(5)
      B(3,5) = float(5)
      B(1,6) = float(1)
      B(2,6) = float(4)
      B(3,6) = float(3)
!
      C = MatMul(Transpose(A),B)
      call mqc_print(iOut,C,header='matrix c from MatMul')
!
      deallocate(C)
      allocate(C(2,6))
      call dgemm('t','n',2,6,  &
        3,float(1),A,3,B,3,float(0),  &
        C,2)
      call mqc_print(iOut,C,header='matrix c from DGEMM')
!
      end program test
