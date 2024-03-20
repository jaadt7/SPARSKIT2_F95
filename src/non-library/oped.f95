!*==oped.f90 processed by SPAG 8.04RA 12:07 23 Feb 2024
!!SPAG Open source Personal, Educational or Academic User Clemson University  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
!-------
SUBROUTINE oped(N,X,Y,Diag,Ioff,Ndiag)
!======================================================
! this kernel performs a matrix by vector multiplication
! for a diagonally structured  matrix stored in diagonal
! format
!======================================================
   USE ISO_FORTRAN_ENV                 
   IMPLICIT NONE
!
! COMMON variable declarations rewritten by SPAG
!
   INTEGER :: Nmvec , Nope
   COMMON Nope , Nmvec
!
! Dummy argument declarations rewritten by SPAG
!
   INTEGER , INTENT(IN) :: N
   INTEGER , INTENT(IN) :: Ndiag
   REAL(REAL64) , INTENT(IN) , DIMENSION(N) :: X
   REAL(REAL64) , INTENT(INOUT) , DIMENSION(N) :: Y
   REAL(REAL64) , INTENT(IN) , DIMENSION(N,Ndiag) :: Diag
   INTEGER , INTENT(IN) , DIMENSION(Ndiag) :: Ioff
!
! Local variable declarations rewritten by SPAG
!
   INTEGER :: i1 , i2 , io , j , k
!
! End of declarations rewritten by SPAG
!
!DIR$ IVDEP
   DO j = 1 , N
      Y(j) = 0.00
   ENDDO
!
   DO j = 1 , Ndiag
      io = Ioff(j)
      i1 = max0(1,1-io)
      i2 = min0(N,N-io)
!DIR$ IVDEP
      DO k = i1 , i2
         Y(k) = Y(k) + Diag(k,j)*X(k+io)
      ENDDO
   ENDDO
   Nmvec = Nmvec + 1
   Nope = Nope + 2*Ndiag*N
END SUBROUTINE oped
