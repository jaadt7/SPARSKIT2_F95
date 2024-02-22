!*==amux.f90 processed by SPAG 8.04RA 15:47 21 Feb 2024
!!SPAG Open source Personal, Educational or Academic User Clemson University  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
!----------------------------------------------------------------------c
!                          S P A R S K I T                             c
!----------------------------------------------------------------------c
!          BASIC MATRIX-VECTOR OPERATIONS - MATVEC MODULE              c
!         Matrix-vector Mulitiplications and Triang. Solves            c
!----------------------------------------------------------------------c
! contents: (as of Nov 18, 1991)                                       c
!----------                                                            c
! 1) Matrix-vector products:                                           c
!---------------------------                                           c
! amux  : A times a vector. Compressed Sparse Row (CSR) format.        c
! amuxms: A times a vector. Modified Compress Sparse Row format.       c
! atmux : Transp(A) times a vector. CSR format.                        c
! atmuxr: Transp(A) times a vector. CSR format. A rectangular.         c
! amuxe : A times a vector. Ellpack/Itpack (ELL) format.               c
! amuxd : A times a vector. Diagonal (DIA) format.                     c
! amuxj : A times a vector. Jagged Diagonal (JAD) format.              c
! vbrmv : Sparse matrix-full vector product, in VBR format             c
!                                                                      c
! 2) Triangular system solutions:                                      c
!-------------------------------                                       c
! lsol  : Unit Lower Triang. solve. Compressed Sparse Row (CSR) format.c
! ldsol : Lower Triang. solve.  Modified Sparse Row (MSR) format.      c
! lsolc : Unit Lower Triang. solve. Comp. Sparse Column (CSC) format.  c
! ldsolc: Lower Triang. solve. Modified Sparse Column (MSC) format.    c
! ldsoll: Lower Triang. solve with level scheduling. MSR format.       c
! usol  : Unit Upper Triang. solve. Compressed Sparse Row (CSR) format.c
! udsol : Upper Triang. solve.  Modified Sparse Row (MSR) format.      c
! usolc : Unit Upper Triang. solve. Comp. Sparse Column (CSC) format.  c
! udsolc: Upper Triang. solve.  Modified Sparse Column (MSC) format.   c
!----------------------------------------------------------------------c
! 1)     M A T R I X    B Y    V E C T O R     P R O D U C T S         c
!----------------------------------------------------------------------c
SUBROUTINE amux(N,X,Y,A,Ja,Ia)
   USE ISO_FORTRAN_ENV                 
   IMPLICIT NONE
!
! Dummy argument declarations rewritten by SPAG
!
   INTEGER , INTENT(IN) :: N
   REAL(REAL64) , INTENT(IN) , DIMENSION(*) :: X
   REAL(REAL64) , INTENT(OUT) , DIMENSION(*) :: Y
   REAL(REAL64) , INTENT(IN) , DIMENSION(*) :: A
   INTEGER , INTENT(IN) , DIMENSION(*) :: Ja
   INTEGER , INTENT(IN) , DIMENSION(*) :: Ia
!
! Local variable declarations rewritten by SPAG
!
   INTEGER :: i , k
   REAL(REAL64) :: t
!
! End of declarations rewritten by SPAG
!
!-----------------------------------------------------------------------
!         A times a vector
!-----------------------------------------------------------------------
! multiplies a matrix by a vector using the dot product form
! Matrix A is stored in compressed sparse row storage.
!
! on entry:
!----------
! n     = row dimension of A
! x     = real array of length equal to the column dimension of
!         the A matrix.
! a, ja,
!    ia = input matrix in compressed sparse row format.
!
! on return:
!-----------
! y     = real array of length n, containing the product y=Ax
!
!-----------------------------------------------------------------------
! local variables
!
!-----------------------------------------------------------------------
   DO i = 1 , N
!
!     compute the inner product of row i with vector x
!
      t = 0.0D0
      DO k = Ia(i) , Ia(i+1) - 1
         t = t + A(k)*X(Ja(k))
      ENDDO
!
!     store result in y(i)
!
      Y(i) = t
   ENDDO
!
!---------end-of-amux---------------------------------------------------
!-----------------------------------------------------------------------
END SUBROUTINE amux
!*==amuxms.f90 processed by SPAG 8.04RA 15:47 21 Feb 2024
!!SPAG Open source Personal, Educational or Academic User Clemson University  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
!-----------------------------------------------------------------------
SUBROUTINE amuxms(N,X,Y,A,Ja)
   USE ISO_FORTRAN_ENV                 
   IMPLICIT NONE
!
! Dummy argument declarations rewritten by SPAG
!
   INTEGER , INTENT(IN) :: N
   REAL(REAL64) , INTENT(IN) , DIMENSION(*) :: X
   REAL(REAL64) , INTENT(INOUT) , DIMENSION(*) :: Y
   REAL(REAL64) , INTENT(IN) , DIMENSION(*) :: A
   INTEGER , INTENT(IN) , DIMENSION(*) :: Ja
!
! Local variable declarations rewritten by SPAG
!
   INTEGER :: i , k
!
! End of declarations rewritten by SPAG
!
!-----------------------------------------------------------------------
!         A times a vector in MSR format
!-----------------------------------------------------------------------
! multiplies a matrix by a vector using the dot product form
! Matrix A is stored in Modified Sparse Row storage.
!
! on entry:
!----------
! n     = row dimension of A
! x     = real array of length equal to the column dimension of
!         the A matrix.
! a, ja,= input matrix in modified compressed sparse row format.
!
! on return:
!-----------
! y     = real array of length n, containing the product y=Ax
!
!-----------------------------------------------------------------------
! local variables
!
!-----------------------------------------------------------------------
   DO i = 1 , N
      Y(i) = A(i)*X(i)
   ENDDO
   DO i = 1 , N
!
!     compute the inner product of row i with vector x
!
      DO k = Ja(i) , Ja(i+1) - 1
         Y(i) = Y(i) + A(k)*X(Ja(k))
      ENDDO
   ENDDO
!
!---------end-of-amuxm--------------------------------------------------
!-----------------------------------------------------------------------
END SUBROUTINE amuxms
!*==atmux.f90 processed by SPAG 8.04RA 15:47 21 Feb 2024
!!SPAG Open source Personal, Educational or Academic User Clemson University  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
!-----------------------------------------------------------------------
SUBROUTINE atmux(N,X,Y,A,Ja,Ia)
   USE ISO_FORTRAN_ENV                 
   IMPLICIT NONE
!
! Dummy argument declarations rewritten by SPAG
!
   INTEGER , INTENT(IN) :: N
   REAL(REAL64) , INTENT(IN) , DIMENSION(*) :: X
   REAL(REAL64) , INTENT(INOUT) , DIMENSION(*) :: Y
   REAL(REAL64) , INTENT(IN) , DIMENSION(*) :: A
   INTEGER , INTENT(IN) , DIMENSION(*) :: Ja
   INTEGER , INTENT(IN) , DIMENSION(*) :: Ia
!
! Local variable declarations rewritten by SPAG
!
   INTEGER :: i , k
!
! End of declarations rewritten by SPAG
!
!-----------------------------------------------------------------------
!         transp( A ) times a vector
!-----------------------------------------------------------------------
! multiplies the transpose of a matrix by a vector when the original
! matrix is stored in compressed sparse row storage. Can also be
! viewed as the product of a matrix by a vector when the original
! matrix is stored in the compressed sparse column format.
!-----------------------------------------------------------------------
!
! on entry:
!----------
! n     = row dimension of A
! x     = real array of length equal to the column dimension of
!         the A matrix.
! a, ja,
!    ia = input matrix in compressed sparse row format.
!
! on return:
!-----------
! y     = real array of length n, containing the product y=transp(A)*x
!
!-----------------------------------------------------------------------
!     local variables
!
!-----------------------------------------------------------------------
!
!     zero out output vector
!
   DO i = 1 , N
      Y(i) = 0.0
   ENDDO
!
! loop over the rows
!
   DO i = 1 , N
      DO k = Ia(i) , Ia(i+1) - 1
         Y(Ja(k)) = Y(Ja(k)) + X(i)*A(k)
      ENDDO
   ENDDO
!
!-------------end-of-atmux----------------------------------------------
!-----------------------------------------------------------------------
END SUBROUTINE atmux
!*==atmuxr.f90 processed by SPAG 8.04RA 15:47 21 Feb 2024
!!SPAG Open source Personal, Educational or Academic User Clemson University  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
!-----------------------------------------------------------------------
SUBROUTINE atmuxr(M,N,X,Y,A,Ja,Ia)
   USE ISO_FORTRAN_ENV                 
   IMPLICIT NONE
!
! Dummy argument declarations rewritten by SPAG
!
   INTEGER , INTENT(IN) :: M
   INTEGER , INTENT(IN) :: N
   REAL(REAL64) , INTENT(IN) , DIMENSION(*) :: X
   REAL(REAL64) , INTENT(INOUT) , DIMENSION(*) :: Y
   REAL(REAL64) , INTENT(IN) , DIMENSION(*) :: A
   INTEGER , INTENT(IN) , DIMENSION(*) :: Ja
   INTEGER , INTENT(IN) , DIMENSION(*) :: Ia
!
! Local variable declarations rewritten by SPAG
!
   INTEGER :: i , k
!
! End of declarations rewritten by SPAG
!
!-----------------------------------------------------------------------
!         transp( A ) times a vector, A can be rectangular
!-----------------------------------------------------------------------
! See also atmux.  The essential difference is how the solution vector
! is initially zeroed.  If using this to multiply rectangular CSC
! matrices by a vector, m number of rows, n is number of columns.
!-----------------------------------------------------------------------
!
! on entry:
!----------
! m     = column dimension of A
! n     = row dimension of A
! x     = real array of length equal to the column dimension of
!         the A matrix.
! a, ja,
!    ia = input matrix in compressed sparse row format.
!
! on return:
!-----------
! y     = real array of length n, containing the product y=transp(A)*x
!
!-----------------------------------------------------------------------
!     local variables
!
!-----------------------------------------------------------------------
!
!     zero out output vector
!
   DO i = 1 , M
      Y(i) = 0.0
   ENDDO
!
! loop over the rows
!
   DO i = 1 , N
      DO k = Ia(i) , Ia(i+1) - 1
         Y(Ja(k)) = Y(Ja(k)) + X(i)*A(k)
      ENDDO
   ENDDO
!
!-------------end-of-atmuxr---------------------------------------------
!-----------------------------------------------------------------------
END SUBROUTINE atmuxr
!*==amuxe.f90 processed by SPAG 8.04RA 15:47 21 Feb 2024
!!SPAG Open source Personal, Educational or Academic User Clemson University  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
!-----------------------------------------------------------------------
SUBROUTINE amuxe(N,X,Y,Na,Ncol,A,Ja)
   USE ISO_FORTRAN_ENV                 
   IMPLICIT NONE
!
! Dummy argument declarations rewritten by SPAG
!
   INTEGER , INTENT(IN) :: N
   INTEGER , INTENT(IN) :: Na
   REAL(REAL64) , INTENT(IN) , DIMENSION(N) :: X
   REAL(REAL64) , INTENT(INOUT) , DIMENSION(N) :: Y
   INTEGER , INTENT(IN) :: Ncol
   REAL(REAL64) , INTENT(IN) , DIMENSION(Na,*) :: A
   INTEGER , INTENT(IN) , DIMENSION(Na,*) :: Ja
!
! Local variable declarations rewritten by SPAG
!
   INTEGER :: i , j
!
! End of declarations rewritten by SPAG
!
!-----------------------------------------------------------------------
!        A times a vector in Ellpack Itpack format (ELL)
!-----------------------------------------------------------------------
! multiplies a matrix by a vector when the original matrix is stored
! in the ellpack-itpack sparse format.
!-----------------------------------------------------------------------
!
! on entry:
!----------
! n     = row dimension of A
! x     = real array of length equal to the column dimension of
!         the A matrix.
! na    = integer. The first dimension of arrays a and ja
!         as declared by the calling program.
! ncol  = integer. The number of active columns in array a.
!         (i.e., the number of generalized diagonals in matrix.)
! a, ja = the real and integer arrays of the itpack format
!         (a(i,k),k=1,ncol contains the elements of row i in matrix
!          ja(i,k),k=1,ncol contains their column numbers)
!
! on return:
!-----------
! y     = real array of length n, containing the product y=y=A*x
!
!-----------------------------------------------------------------------
! local variables
!
!-----------------------------------------------------------------------
   DO i = 1 , N
      Y(i) = 0.0
   ENDDO
   DO j = 1 , Ncol
      DO i = 1 , N
         Y(i) = Y(i) + A(i,j)*X(Ja(i,j))
      ENDDO
   ENDDO
!
!--------end-of-amuxe---------------------------------------------------
!-----------------------------------------------------------------------
END SUBROUTINE amuxe
!*==amuxd.f90 processed by SPAG 8.04RA 15:47 21 Feb 2024
!!SPAG Open source Personal, Educational or Academic User Clemson University  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
!-----------------------------------------------------------------------
SUBROUTINE amuxd(N,X,Y,Diag,Ndiag,Idiag,Ioff)
   USE ISO_FORTRAN_ENV                 
   IMPLICIT NONE
!
! Dummy argument declarations rewritten by SPAG
!
   INTEGER , INTENT(IN) :: N
   INTEGER , INTENT(IN) :: Ndiag
   INTEGER , INTENT(IN) :: Idiag
   REAL(REAL64) , INTENT(IN) , DIMENSION(N) :: X
   REAL(REAL64) , INTENT(INOUT) , DIMENSION(N) :: Y
   REAL(REAL64) , INTENT(IN) , DIMENSION(Ndiag,Idiag) :: Diag
   INTEGER , INTENT(IN) , DIMENSION(Idiag) :: Ioff
!
! Local variable declarations rewritten by SPAG
!
   INTEGER :: i1 , i2 , io , j , k
!
! End of declarations rewritten by SPAG
!
!-----------------------------------------------------------------------
!        A times a vector in Diagonal storage format (DIA)
!-----------------------------------------------------------------------
! multiplies a matrix by a vector when the original matrix is stored
! in the diagonal storage format.
!-----------------------------------------------------------------------
!
! on entry:
!----------
! n     = row dimension of A
! x     = real array of length equal to the column dimension of
!         the A matrix.
! ndiag  = integer. The first dimension of array adiag as declared in
!         the calling program.
! idiag  = integer. The number of diagonals in the matrix.
! diag   = real array containing the diagonals stored of A.
! idiag  = number of diagonals in matrix.
! diag   = real array of size (ndiag x idiag) containing the diagonals
!
! ioff   = integer array of length idiag, containing the offsets of the
!   	   diagonals of the matrix:
!          diag(i,k) contains the element a(i,i+ioff(k)) of the matrix.
!
! on return:
!-----------
! y     = real array of length n, containing the product y=A*x
!
!-----------------------------------------------------------------------
! local variables
!
!-----------------------------------------------------------------------
   DO j = 1 , N
      Y(j) = 0.0D0
   ENDDO
   DO j = 1 , Idiag
      io = Ioff(j)
      i1 = max0(1,1-io)
      i2 = min0(N,N-io)
      DO k = i1 , i2
         Y(k) = Y(k) + Diag(k,j)*X(k+io)
      ENDDO
   ENDDO
!
!----------end-of-amuxd-------------------------------------------------
!-----------------------------------------------------------------------
END SUBROUTINE amuxd
!*==amuxj.f90 processed by SPAG 8.04RA 15:47 21 Feb 2024
!!SPAG Open source Personal, Educational or Academic User Clemson University  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
!-----------------------------------------------------------------------
SUBROUTINE amuxj(N,X,Y,Jdiag,A,Ja,Ia)
   USE ISO_FORTRAN_ENV                 
   IMPLICIT NONE
!
! Dummy argument declarations rewritten by SPAG
!
   INTEGER , INTENT(IN) :: N
   REAL(REAL64) , INTENT(IN) , DIMENSION(N) :: X
   REAL(REAL64) , INTENT(INOUT) , DIMENSION(N) :: Y
   INTEGER , INTENT(IN) :: Jdiag
   REAL(REAL64) , INTENT(IN) , DIMENSION(*) :: A
   INTEGER , INTENT(IN) , DIMENSION(*) :: Ja
   INTEGER , INTENT(IN) , DIMENSION(*) :: Ia
!
! Local variable declarations rewritten by SPAG
!
   INTEGER :: i , ii , j , k1 , len
!
! End of declarations rewritten by SPAG
!
!-----------------------------------------------------------------------
!        A times a vector in Jagged-Diagonal storage format (JAD)
!-----------------------------------------------------------------------
! multiplies a matrix by a vector when the original matrix is stored
! in the jagged diagonal storage format.
!-----------------------------------------------------------------------
!
! on entry:
!----------
! n      = row dimension of A
! x      = real array of length equal to the column dimension of
!         the A matrix.
! jdiag  = integer. The number of jadded-diagonals in the data-structure.
! a      = real array containing the jadded diagonals of A stored
!          in succession (in decreasing lengths)
! j      = integer array containing the colum indices of the
!          corresponding elements in a.
! ia     = integer array containing the lengths of the  jagged diagonals
!
! on return:
!-----------
! y      = real array of length n, containing the product y=A*x
!
! Note:
!-------
! Permutation related to the JAD format is not performed.
! this can be done by:
!     call permvec (n,y,y,iperm)
! after the call to amuxj, where iperm is the permutation produced
! by csrjad.
!-----------------------------------------------------------------------
! local variables
!
!-----------------------------------------------------------------------
   DO i = 1 , N
      Y(i) = 0.0D0
   ENDDO
   DO ii = 1 , Jdiag
      k1 = Ia(ii) - 1
      len = Ia(ii+1) - k1 - 1
      DO j = 1 , len
         Y(j) = Y(j) + A(k1+j)*X(Ja(k1+j))
      ENDDO
   ENDDO
!
!----------end-of-amuxj-------------------------------------------------
!-----------------------------------------------------------------------
END SUBROUTINE amuxj
!*==vbrmv.f90 processed by SPAG 8.04RA 15:47 21 Feb 2024
!!SPAG Open source Personal, Educational or Academic User Clemson University  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
!-----------------------------------------------------------------------
SUBROUTINE vbrmv(Nr,Nc,Ia,Ja,A,Kvstr,Kvstc,X,B)
   USE ISO_FORTRAN_ENV                 
   IMPLICIT NONE
!
! Dummy argument declarations rewritten by SPAG
!
   INTEGER , INTENT(IN) :: Nr
   INTEGER , INTENT(IN) :: Nc
   INTEGER , INTENT(IN) , DIMENSION(Nr+1) :: Ia
   INTEGER , INTENT(IN) , DIMENSION(*) :: Ja
   REAL(REAL64) , INTENT(IN) , DIMENSION(*) :: A
   INTEGER , INTENT(IN) , DIMENSION(Nr+1) :: Kvstr
   INTEGER , INTENT(IN) , DIMENSION(*) :: Kvstc
   REAL(REAL64) , INTENT(IN) , DIMENSION(*) :: X
   REAL(REAL64) , INTENT(INOUT) , DIMENSION(*) :: B
!
! Local variable declarations rewritten by SPAG
!
   INTEGER :: i , ii , istart , istop , j , jj , k , n
   REAL(REAL64) :: xjj
!
! End of declarations rewritten by SPAG
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!     Sparse matrix-full vector product, in VBR format.
!-----------------------------------------------------------------------
!     On entry:
!--------------
!     nr, nc  = number of block rows and columns in matrix A
!     ia,ja,(),a,kvstr,kvstc = matrix A in variable block row format
!     x       = multiplier vector in full format
!
!     On return:
!---------------
!     b = product of matrix A times vector x in full format
!
!     Algorithm:
!---------------
!     Perform multiplication by traversing a in order.
!
!-----------------------------------------------------------------------
!-----local variables
!---------------------------------
   n = Kvstc(Nc+1) - 1
   DO i = 1 , n
      B(i) = 0.D0
   ENDDO
!---------------------------------
   k = 1
   DO i = 1 , Nr
      istart = Kvstr(i)
      istop = Kvstr(i+1) - 1
      DO j = Ia(i) , Ia(i+1) - 1
         DO jj = Kvstc(Ja(j)) , Kvstc(Ja(j)+1) - 1
            xjj = X(jj)
            DO ii = istart , istop
               B(ii) = B(ii) + xjj*A(k)
               k = k + 1
            ENDDO
         ENDDO
      ENDDO
   ENDDO
!---------------------------------
END SUBROUTINE vbrmv
!*==lsol.f90 processed by SPAG 8.04RA 15:47 21 Feb 2024
!!SPAG Open source Personal, Educational or Academic User Clemson University  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
!-----------------------------------------------------------------------
!----------------------end-of-vbrmv-------------------------------------
!-----------------------------------------------------------------------
!----------------------------------------------------------------------c
! 2)     T R I A N G U L A R    S Y S T E M    S O L U T I O N S       c
!----------------------------------------------------------------------c
SUBROUTINE lsol(N,X,Y,Al,Jal,Ial)
   USE ISO_FORTRAN_ENV                 
   IMPLICIT NONE
!
! Dummy argument declarations rewritten by SPAG
!
   INTEGER , INTENT(IN) :: N
   REAL(REAL64) , INTENT(INOUT) , DIMENSION(N) :: X
   REAL(REAL64) , INTENT(IN) , DIMENSION(N) :: Y
   REAL(REAL64) , INTENT(IN) , DIMENSION(*) :: Al
   INTEGER , INTENT(IN) , DIMENSION(*) :: Jal
   INTEGER , INTENT(IN) , DIMENSION(N+1) :: Ial
!
! Local variable declarations rewritten by SPAG
!
   INTEGER :: j , k
   REAL(REAL64) :: t
!
! End of declarations rewritten by SPAG
!
!-----------------------------------------------------------------------
!   solves    L x = y ; L = lower unit triang. /  CSR format
!-----------------------------------------------------------------------
! solves a unit lower triangular system by standard (sequential )
! forward elimination - matrix stored in CSR format.
!-----------------------------------------------------------------------
!
! On entry:
!----------
! n      = integer. dimension of problem.
! y      = real array containg the right side.
!
! al,
! jal,
! ial,    = Lower triangular matrix stored in compressed sparse row
!          format.
!
! On return:
!-----------
!	x  = The solution of  L x  = y.
!--------------------------------------------------------------------
! local variables
!
!-----------------------------------------------------------------------
   X(1) = Y(1)
   DO k = 2 , N
      t = Y(k)
      DO j = Ial(k) , Ial(k+1) - 1
         t = t - Al(j)*X(Jal(j))
      ENDDO
      X(k) = t
   ENDDO
!
!----------end-of-lsol--------------------------------------------------
!-----------------------------------------------------------------------
END SUBROUTINE lsol
!*==ldsol.f90 processed by SPAG 8.04RA 15:47 21 Feb 2024
!!SPAG Open source Personal, Educational or Academic User Clemson University  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
!-----------------------------------------------------------------------
SUBROUTINE ldsol(N,X,Y,Al,Jal)
   USE ISO_FORTRAN_ENV                 
   IMPLICIT NONE
!
! Dummy argument declarations rewritten by SPAG
!
   INTEGER , INTENT(IN) :: N
   REAL(REAL64) , INTENT(INOUT) , DIMENSION(N) :: X
   REAL(REAL64) , INTENT(IN) , DIMENSION(N) :: Y
   REAL(REAL64) , INTENT(IN) , DIMENSION(*) :: Al
   INTEGER , INTENT(IN) , DIMENSION(*) :: Jal
!
! Local variable declarations rewritten by SPAG
!
   INTEGER :: j , k
   REAL(REAL64) :: t
!
! End of declarations rewritten by SPAG
!
!-----------------------------------------------------------------------
!     Solves L x = y    L = triangular. MSR format
!-----------------------------------------------------------------------
! solves a (non-unit) lower triangular system by standard (sequential)
! forward elimination - matrix stored in MSR format
! with diagonal elements already inverted (otherwise do inversion,
! al(1:n) = 1.0/al(1:n),  before calling ldsol).
!-----------------------------------------------------------------------
!
! On entry:
!----------
! n      = integer. dimension of problem.
! y      = real array containg the right hand side.
!
! al,
! jal,   = Lower triangular matrix stored in Modified Sparse Row
!          format.
!
! On return:
!-----------
!	x = The solution of  L x = y .
!--------------------------------------------------------------------
! local variables
!
!-----------------------------------------------------------------------
   X(1) = Y(1)*Al(1)
   DO k = 2 , N
      t = Y(k)
      DO j = Jal(k) , Jal(k+1) - 1
         t = t - Al(j)*X(Jal(j))
      ENDDO
      X(k) = Al(k)*t
   ENDDO
!----------end-of-ldsol-------------------------------------------------
!-----------------------------------------------------------------------
END SUBROUTINE ldsol
!*==lsolc.f90 processed by SPAG 8.04RA 15:47 21 Feb 2024
!!SPAG Open source Personal, Educational or Academic User Clemson University  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
!-----------------------------------------------------------------------
SUBROUTINE lsolc(N,X,Y,Al,Jal,Ial)
   USE ISO_FORTRAN_ENV                 
   IMPLICIT NONE
!
! Dummy argument declarations rewritten by SPAG
!
   INTEGER , INTENT(IN) :: N
   REAL(REAL64) , INTENT(INOUT) , DIMENSION(N) :: X
   REAL(REAL64) , INTENT(IN) , DIMENSION(N) :: Y
   REAL(REAL64) , INTENT(IN) , DIMENSION(*) :: Al
   INTEGER , INTENT(IN) , DIMENSION(*) :: Jal
   INTEGER , INTENT(IN) , DIMENSION(*) :: Ial
!
! Local variable declarations rewritten by SPAG
!
   INTEGER :: j , k
   REAL(REAL64) :: t
!
! End of declarations rewritten by SPAG
!
!-----------------------------------------------------------------------
!       SOLVES     L x = y ;    where L = unit lower trang. CSC format
!-----------------------------------------------------------------------
! solves a unit lower triangular system by standard (sequential )
! forward elimination - matrix stored in CSC format.
!-----------------------------------------------------------------------
!
! On entry:
!----------
! n      = integer. dimension of problem.
! y      = real*8 array containg the right side.
!
! al,
! jal,
! ial,    = Lower triangular matrix stored in compressed sparse column
!          format.
!
! On return:
!-----------
!	x  = The solution of  L x  = y.
!-----------------------------------------------------------------------
! local variables
!
!-----------------------------------------------------------------------
   DO k = 1 , N
      X(k) = Y(k)
   ENDDO
   DO k = 1 , N - 1
      t = X(k)
      DO j = Ial(k) , Ial(k+1) - 1
         X(Jal(j)) = X(Jal(j)) - t*Al(j)
      ENDDO
   ENDDO
!
!----------end-of-lsolc-------------------------------------------------
!-----------------------------------------------------------------------
END SUBROUTINE lsolc
!*==ldsolc.f90 processed by SPAG 8.04RA 15:47 21 Feb 2024
!!SPAG Open source Personal, Educational or Academic User Clemson University  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
!-----------------------------------------------------------------------
SUBROUTINE ldsolc(N,X,Y,Al,Jal)
   USE ISO_FORTRAN_ENV                 
   IMPLICIT NONE
!
! Dummy argument declarations rewritten by SPAG
!
   INTEGER , INTENT(IN) :: N
   REAL(REAL64) , INTENT(INOUT) , DIMENSION(N) :: X
   REAL(REAL64) , INTENT(IN) , DIMENSION(N) :: Y
   REAL(REAL64) , INTENT(IN) , DIMENSION(*) :: Al
   INTEGER , INTENT(IN) , DIMENSION(*) :: Jal
!
! Local variable declarations rewritten by SPAG
!
   INTEGER :: j , k
   REAL(REAL64) :: t
!
! End of declarations rewritten by SPAG
!
!-----------------------------------------------------------------------
!    Solves     L x = y ;    L = nonunit Low. Triang. MSC format
!-----------------------------------------------------------------------
! solves a (non-unit) lower triangular system by standard (sequential)
! forward elimination - matrix stored in Modified Sparse Column format
! with diagonal elements already inverted (otherwise do inversion,
! al(1:n) = 1.0/al(1:n),  before calling ldsol).
!-----------------------------------------------------------------------
!
! On entry:
!----------
! n      = integer. dimension of problem.
! y      = real array containg the right hand side.
!
! al,
! jal,
! ial,    = Lower triangular matrix stored in Modified Sparse Column
!           format.
!
! On return:
!-----------
!	x = The solution of  L x = y .
!--------------------------------------------------------------------
! local variables
!
!-----------------------------------------------------------------------
   DO k = 1 , N
      X(k) = Y(k)
   ENDDO
   DO k = 1 , N
      X(k) = X(k)*Al(k)
      t = X(k)
      DO j = Jal(k) , Jal(k+1) - 1
         X(Jal(j)) = X(Jal(j)) - t*Al(j)
      ENDDO
   ENDDO
!
!----------end-of-lsolc------------------------------------------------
!-----------------------------------------------------------------------
END SUBROUTINE ldsolc
!*==ldsoll.f90 processed by SPAG 8.04RA 15:47 21 Feb 2024
!!SPAG Open source Personal, Educational or Academic User Clemson University  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
!-----------------------------------------------------------------------
SUBROUTINE ldsoll(N,X,Y,Al,Jal,Nlev,Lev,Ilev)
   USE ISO_FORTRAN_ENV                 
   IMPLICIT NONE
!
! Dummy argument declarations rewritten by SPAG
!
   INTEGER , INTENT(IN) :: N
   INTEGER , INTENT(IN) :: Nlev
   REAL(REAL64) , INTENT(INOUT) , DIMENSION(N) :: X
   REAL(REAL64) , INTENT(IN) , DIMENSION(N) :: Y
   REAL(REAL64) , INTENT(IN) , DIMENSION(*) :: Al
   INTEGER , INTENT(IN) , DIMENSION(*) :: Jal
   INTEGER , INTENT(IN) , DIMENSION(N) :: Lev
   INTEGER , INTENT(IN) , DIMENSION(Nlev+1) :: Ilev
!
! Local variable declarations rewritten by SPAG
!
   INTEGER :: i , ii , jrow , k
   REAL(REAL64) :: t
!
! End of declarations rewritten by SPAG
!
!-----------------------------------------------------------------------
!    Solves L x = y    L = triangular. Uses LEVEL SCHEDULING/MSR format
!-----------------------------------------------------------------------
!
! On entry:
!----------
! n      = integer. dimension of problem.
! y      = real array containg the right hand side.
!
! al,
! jal,   = Lower triangular matrix stored in Modified Sparse Row
!          format.
! nlev   = number of levels in matrix
! lev    = integer array of length n, containing the permutation
!          that defines the levels in the level scheduling ordering.
! ilev   = pointer to beginning of levels in lev.
!          the numbers lev(i) to lev(i+1)-1 contain the row numbers
!          that belong to level number i, in the level shcheduling
!          ordering.
!
! On return:
!-----------
!	x = The solution of  L x = y .
!--------------------------------------------------------------------
!
!     outer loop goes through the levels. (SEQUENTIAL loop)
!
   DO ii = 1 , Nlev
!
!     next loop executes within the same level. PARALLEL loop
!
      DO i = Ilev(ii) , Ilev(ii+1) - 1
         jrow = Lev(i)
!
! compute inner product of row jrow with x
!
         t = Y(jrow)
         DO k = Jal(jrow) , Jal(jrow+1) - 1
            t = t - Al(k)*X(Jal(k))
         ENDDO
         X(jrow) = t*Al(jrow)
      ENDDO
   ENDDO
!-----------------------------------------------------------------------
END SUBROUTINE ldsoll
!*==usol.f90 processed by SPAG 8.04RA 15:47 21 Feb 2024
!!SPAG Open source Personal, Educational or Academic User Clemson University  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
!-----------------------------------------------------------------------
SUBROUTINE usol(N,X,Y,Au,Jau,Iau)
   USE ISO_FORTRAN_ENV                 
   IMPLICIT NONE
!
! Dummy argument declarations rewritten by SPAG
!
   INTEGER , INTENT(IN) :: N
   REAL(REAL64) , INTENT(INOUT) , DIMENSION(N) :: X
   REAL(REAL64) , INTENT(IN) , DIMENSION(N) :: Y
   REAL(REAL64) , INTENT(IN) , DIMENSION(*) :: Au
   INTEGER , INTENT(IN) , DIMENSION(*) :: Jau
   INTEGER , INTENT(IN) , DIMENSION(N+1) :: Iau
!
! Local variable declarations rewritten by SPAG
!
   INTEGER :: j , k
   REAL(REAL64) :: t
!
! End of declarations rewritten by SPAG
!
!-----------------------------------------------------------------------
!             Solves   U x = y    U = unit upper triangular.
!-----------------------------------------------------------------------
! solves a unit upper triangular system by standard (sequential )
! backward elimination - matrix stored in CSR format.
!-----------------------------------------------------------------------
!
! On entry:
!----------
! n      = integer. dimension of problem.
! y      = real array containg the right side.
!
! au,
! jau,
! iau,    = Lower triangular matrix stored in compressed sparse row
!          format.
!
! On return:
!-----------
!	x = The solution of  U x = y .
!--------------------------------------------------------------------
! local variables
!
!-----------------------------------------------------------------------
   X(N) = Y(N)
   DO k = N - 1 , 1 , -1
      t = Y(k)
      DO j = Iau(k) , Iau(k+1) - 1
         t = t - Au(j)*X(Jau(j))
      ENDDO
      X(k) = t
   ENDDO
!
!----------end-of-usol--------------------------------------------------
!-----------------------------------------------------------------------
END SUBROUTINE usol
!*==udsol.f90 processed by SPAG 8.04RA 15:47 21 Feb 2024
!!SPAG Open source Personal, Educational or Academic User Clemson University  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
!-----------------------------------------------------------------------
SUBROUTINE udsol(N,X,Y,Au,Jau)
   USE ISO_FORTRAN_ENV                 
   IMPLICIT NONE
!
! Dummy argument declarations rewritten by SPAG
!
   INTEGER , INTENT(IN) :: N
   REAL(REAL64) , INTENT(INOUT) , DIMENSION(N) :: X
   REAL(REAL64) , INTENT(IN) , DIMENSION(N) :: Y
   REAL(REAL64) , INTENT(IN) , DIMENSION(*) :: Au
   INTEGER , INTENT(IN) , DIMENSION(*) :: Jau
!
! Local variable declarations rewritten by SPAG
!
   INTEGER :: j , k
   REAL(REAL64) :: t
!
! End of declarations rewritten by SPAG
!
!-----------------------------------------------------------------------
!             Solves   U x = y  ;   U = upper triangular in MSR format
!-----------------------------------------------------------------------
! solves a non-unit upper triangular matrix by standard (sequential )
! backward elimination - matrix stored in MSR format.
! with diagonal elements already inverted (otherwise do inversion,
! au(1:n) = 1.0/au(1:n),  before calling).
!-----------------------------------------------------------------------
!
! On entry:
!----------
! n      = integer. dimension of problem.
! y      = real array containg the right side.
!
! au,
! jau,    = Lower triangular matrix stored in modified sparse row
!          format.
!
! On return:
!-----------
!	x = The solution of  U x = y .
!--------------------------------------------------------------------
! local variables
!
!-----------------------------------------------------------------------
   X(N) = Y(N)*Au(N)
   DO k = N - 1 , 1 , -1
      t = Y(k)
      DO j = Jau(k) , Jau(k+1) - 1
         t = t - Au(j)*X(Jau(j))
      ENDDO
      X(k) = Au(k)*t
   ENDDO
!
!----------end-of-udsol-------------------------------------------------
!-----------------------------------------------------------------------
END SUBROUTINE udsol
!*==usolc.f90 processed by SPAG 8.04RA 15:47 21 Feb 2024
!!SPAG Open source Personal, Educational or Academic User Clemson University  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
!-----------------------------------------------------------------------
SUBROUTINE usolc(N,X,Y,Au,Jau,Iau)
   USE ISO_FORTRAN_ENV                 
   IMPLICIT NONE
!
! Dummy argument declarations rewritten by SPAG
!
   INTEGER , INTENT(IN) :: N
   REAL(REAL64) , INTENT(INOUT) , DIMENSION(*) :: X
   REAL(REAL64) , INTENT(IN) , DIMENSION(*) :: Y
   REAL(REAL64) , INTENT(IN) , DIMENSION(*) :: Au
   INTEGER , INTENT(IN) , DIMENSION(*) :: Jau
   INTEGER , INTENT(IN) , DIMENSION(*) :: Iau
!
! Local variable declarations rewritten by SPAG
!
   INTEGER :: j , k
   REAL(REAL64) :: t
!
! End of declarations rewritten by SPAG
!
!-----------------------------------------------------------------------
!       SOUVES     U x = y ;    where U = unit upper trang. CSC format
!-----------------------------------------------------------------------
! solves a unit upper triangular system by standard (sequential )
! forward elimination - matrix stored in CSC format.
!-----------------------------------------------------------------------
!
! On entry:
!----------
! n      = integer. dimension of problem.
! y      = real*8 array containg the right side.
!
! au,
! jau,
! iau,    = Uower triangular matrix stored in compressed sparse column
!          format.
!
! On return:
!-----------
!	x  = The solution of  U x  = y.
!-----------------------------------------------------------------------
! local variables
!
!-----------------------------------------------------------------------
   DO k = 1 , N
      X(k) = Y(k)
   ENDDO
   DO k = N , 1 , -1
      t = X(k)
      DO j = Iau(k) , Iau(k+1) - 1
         X(Jau(j)) = X(Jau(j)) - t*Au(j)
      ENDDO
   ENDDO
!
!----------end-of-usolc-------------------------------------------------
!-----------------------------------------------------------------------
END SUBROUTINE usolc
!*==udsolc.f90 processed by SPAG 8.04RA 15:47 21 Feb 2024
!!SPAG Open source Personal, Educational or Academic User Clemson University  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
!-----------------------------------------------------------------------
SUBROUTINE udsolc(N,X,Y,Au,Jau)
   USE ISO_FORTRAN_ENV                 
   IMPLICIT NONE
!
! Dummy argument declarations rewritten by SPAG
!
   INTEGER , INTENT(IN) :: N
   REAL(REAL64) , INTENT(INOUT) , DIMENSION(N) :: X
   REAL(REAL64) , INTENT(IN) , DIMENSION(N) :: Y
   REAL(REAL64) , INTENT(IN) , DIMENSION(*) :: Au
   INTEGER , INTENT(IN) , DIMENSION(*) :: Jau
!
! Local variable declarations rewritten by SPAG
!
   INTEGER :: j , k
   REAL(REAL64) :: t
!
! End of declarations rewritten by SPAG
!
!-----------------------------------------------------------------------
!    Solves     U x = y ;    U = nonunit Up. Triang. MSC format
!-----------------------------------------------------------------------
! solves a (non-unit) upper triangular system by standard (sequential)
! forward elimination - matrix stored in Modified Sparse Column format
! with diagonal elements already inverted (otherwise do inversion,
! auuuul(1:n) = 1.0/au(1:n),  before calling ldsol).
!-----------------------------------------------------------------------
!
! On entry:
!----------
! n      = integer. dimension of problem.
! y      = real*8 array containg the right hand side.
!
! au,
! jau,   = Upper triangular matrix stored in Modified Sparse Column
!          format.
!
! On return:
!-----------
!	x = The solution of  U x = y .
!--------------------------------------------------------------------
! local variables
!
!-----------------------------------------------------------------------
   DO k = 1 , N
      X(k) = Y(k)
   ENDDO
   DO k = N , 1 , -1
      X(k) = X(k)*Au(k)
      t = X(k)
      DO j = Jau(k) , Jau(k+1) - 1
         X(Jau(j)) = X(Jau(j)) - t*Au(j)
      ENDDO
   ENDDO
!
!----------end-of-udsolc------------------------------------------------
!-----------------------------------------------------------------------
END SUBROUTINE udsolc
