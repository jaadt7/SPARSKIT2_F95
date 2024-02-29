!*==phitest.f90 processed by SPAG 8.04RA 12:07 23 Feb 2024
!!SPAG Open source Personal, Educational or Academic User Clemson University  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
!
PROGRAM phitest
!-------------------------------------------------------------------
! Test program for exponential propagator using Arnoldi approach
! This main program is a very simple test using diagonal matrices
! (Krylov subspace methods are blind to the structure of the matrix
! except for symmetry). This provides a good way of testing the
! accuracy of the method as well as the error estimates. This test
! program tests the phi-function variant instead of the exponential.
!
! The subroutine phipro which is tested here computes an approximation
! to the vector
!
!   	w(tn) = w(t0) + tn *  phi( - A * tn ) * (r - A w(t0))
!
! where phi(z) = (1-exp(z)) / z
!
! i.e. it solves dw/dt = -A w + r in [t0,t0+ tn] (returns only w(t0+tn))
!
! In this program t0=0, tn is input by the user and A is a simple
! diagonal matrix.
!
!-------------------------------------------------------------------
   USE ISO_FORTRAN_ENV                 
   IMPLICIT NONE
!
! PARAMETER definitions rewritten by SPAG
!
   INTEGER , PARAMETER :: NMAX = 400 , IH0 = 60 , NZMAX = 7*NMAX
!
! Local variable declarations rewritten by SPAG
!
   REAL(REAL64) , DIMENSION(NZMAX) :: a
   REAL(REAL64) , SAVE :: a0 , b0 , eps , epsmac
   REAL(REAL64) , EXTERNAL :: ddot
   REAL(REAL64) :: h , t , tn
   INTEGER , DIMENSION(10) :: ioff
   INTEGER , SAVE :: iout
   INTEGER :: j , k , m , n , ndiag
   REAL(REAL64) , DIMENSION(NMAX) :: r , w , w1 , x , y
   REAL(REAL64) , DIMENSION(IH0*NMAX) :: u
   LOGICAL verboz
   EXTERNAL phiprod
!
! End of declarations rewritten by SPAG
!
   DATA iout/6/ , a0/0.0/ , b0/1.0/ , epsmac/1.D-10/ , eps/1.D-10/
!
!     dimendion of matrix
!
   n = 100
!---------------------------define matrix -----------------------------
!     A is a single diagonal matrix (ndiag = 1 and ioff(1) = 0 )
!-----------------------------------------------------------------------
   ndiag = 1
   ioff(1) = 0
!
!     entriesin the diagonal are uniformly distributed.
!
   h = 1.0D0/real(n+1)
   DO j = 1 , n
      a(j) = real(j+1)*h
   ENDDO
!--------
   WRITE (6,'(10hEnter tn: ,$)')
   READ (5,*) tn
!
   WRITE (6,'(36hEpsilon (desired relative accuracy): ,$)')
   READ (5,*) eps
!-------
   WRITE (6,'(36h m (= dimension of Krylov subspace): ,$)')
   READ (5,*) m
!
!     define initial conditions to be (1,1,1,1..1)^T
!
   DO j = 1 , n
      w(j) = 1.0 - 1.0/real(j+10)
      r(j) = 1.0
      w1(j) = w(j)
   ENDDO
!
!
!     HERE IT IS DONE IN 10 SUBSTEPS
!     nsteps = 10
!     tnh = tn/real(nsteps)
!     MARCHING LOOP
!     do jj =1, nsteps
!     call phiprod (n, m, eps, tnh, u, w, r, x, y, a, ioff, ndiag)
   verboz = .true.
   CALL phiprod(n,m,eps,tn,u,w,r,x,y,a,ioff,ndiag,verboz)
!     enddo
   PRINT * , ' final answer '
   PRINT * , (w(k),k=1,20)
!
!     exact solution
!
   DO k = 1 , n
      w1(k) = w1(k) + ((1.0-dexp(-a(k)*tn))/a(k))*(r(k)-a(k)*w1(k))
   ENDDO
!
!     exact solution
!
   PRINT * , ' exact answer '
   PRINT * , (w1(k),k=1,20)
!
!     computing actual 2-norm of error
!
   t = 0.0D0
   DO k = 1 , n
      t = t + (w1(k)-w(k))**2
   ENDDO
   t = dsqrt(t/ddot(n,w,1,w,1))
!
   WRITE (6,*) ' final error' , t
!-----------------------------------------------------------------------
END PROGRAM phitest
!*==oped.f90 processed by SPAG 8.04RA 12:07 23 Feb 2024
!!SPAG Open source Personal, Educational or Academic User Clemson University  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
!-----------------------------------------------------------------------
SUBROUTINE oped(N,X,Y,Diag,Ioff,Ndiag)
!======================================================
! this kernel performs a matrix by vector multiplication
! for a diagonally structured  matrix stored in diagonal format
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
!
