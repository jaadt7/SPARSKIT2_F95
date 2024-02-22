!*==sobel.f90 processed by SPAG 8.04RA 15:47 21 Feb 2024
!!SPAG Open source Personal, Educational or Academic User Clemson University  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
SUBROUTINE sobel(N,Nrowc,Ncolc,C,Jc,Ic,A,Ja,Ia,B,Jb,Ib,Nzmax,Ierr)
   USE ISO_FORTRAN_ENV                 
   IMPLICIT NONE
!
! Dummy argument declarations rewritten by SPAG
!
   INTEGER , INTENT(IN) :: N
   INTEGER :: Nrowc
   INTEGER :: Ncolc
   REAL(REAL64) , DIMENSION(*) :: C
   INTEGER , DIMENSION(*) :: Jc
   INTEGER , DIMENSION(*) :: Ic
   REAL(REAL64) , DIMENSION(*) :: A
   INTEGER , DIMENSION(*) :: Ja
   INTEGER , DIMENSION(*) :: Ia
   REAL(REAL64) , DIMENSION(*) :: B
   INTEGER , DIMENSION(*) :: Jb
   INTEGER , DIMENSION(*) :: Ib
   INTEGER :: Nzmax
   INTEGER :: Ierr
!
! Local variable declarations rewritten by SPAG
!
   INTEGER :: i , ipos , ncola , ncolb , nrowa , nrowb , offset
   EXTERNAL addblk , copmat , diagblk , leftblk , rightblk
!
! End of declarations rewritten by SPAG
!
!-----------------------------------------------------------------------
!     This subroutine generates a matrix used in the statistical problem
!     presented by Prof. Sobel. The matrix is formed by a series of
!     small submatrix on or adjancent to the diagonal. The submatrix on
!     the diagonal is square and the size goes like 1, 1, 2, 2, 3, 3,...
!     Each of the diagonal block is a triadiagonal matrix, each of the
!     off-diagonal block is a bidiagonal block. The values of elements
!     in the off-diagonal block are all -1. So are the values of the
!     elements on the sub- and super-diagonal of the blocks on the
!     diagonal. The first element(1,1) of the diagonal block is alternating
!     between 3 and 5, the rest of the diagonal elements (of the block
!     on the diagonal) are 6.
!-----------------------------------------------------------------------
!     This subroutine calls following subroutines to generate the three
!     thypes of submatrices:
!     diagblk -- generates diagonal block.
!     leftblk -- generates the block left of the diagonal one.
!     rightblk-- generates the block right of the diagonal one.
!-----------------------------------------------------------------------
   IF ( N<2 ) RETURN
 
   ipos = 1
   offset = 1
   CALL diagblk(1,Nrowc,Ncolc,C,Jc,Ic)
   DO i = 2 , N - 2
      nrowa = Nrowc
      ncola = Ncolc
      CALL copmat(Nrowc,C,Jc,Ic,A,Ja,Ia,1,1)
      CALL rightblk(i-1,nrowb,ncolb,B,Jb,Ib)
      CALL addblk(nrowa,ncola,A,Ja,Ia,ipos,ipos+offset,1,nrowb,ncolb,B,Jb,Ib,Nrowc,Ncolc,C,Jc,Ic,Nzmax,Ierr)
      CALL leftblk(i,nrowb,ncolb,B,Jb,Ib)
      CALL addblk(Nrowc,Ncolc,C,Jc,Ic,ipos+offset,ipos,1,nrowb,ncolb,B,Jb,Ib,nrowa,ncola,A,Ja,Ia,Nzmax,Ierr)
      ipos = ipos + offset
      CALL diagblk(i,nrowb,ncolb,B,Jb,Ib)
      CALL addblk(nrowa,ncola,A,Ja,Ia,ipos,ipos,1,nrowb,ncolb,B,Jb,Ib,Nrowc,Ncolc,C,Jc,Ic,Nzmax,Ierr)
      offset = 1 + (i-1)/2
   ENDDO
END SUBROUTINE sobel
!*==diagblk.f90 processed by SPAG 8.04RA 15:47 21 Feb 2024
!!SPAG Open source Personal, Educational or Academic User Clemson University  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
!-----------------------------------------------------------------------
SUBROUTINE diagblk(N,Nrow,Ncol,A,Ja,Ia)
   USE ISO_FORTRAN_ENV                 
   IMPLICIT NONE
!
! Dummy argument declarations rewritten by SPAG
!
   INTEGER , INTENT(IN) :: N
   INTEGER , INTENT(INOUT) :: Nrow
   INTEGER , INTENT(INOUT) :: Ncol
   REAL(REAL64) , INTENT(OUT) , DIMENSION(1:*) :: A
   INTEGER , INTENT(OUT) , DIMENSION(1:*) :: Ja
   INTEGER , INTENT(OUT) , DIMENSION(1:*) :: Ia
!
! Local variable declarations rewritten by SPAG
!
   INTEGER :: i , k
!
! End of declarations rewritten by SPAG
!
!-----------------------------------------------------------------------
!     generates the diagonal block for the given problem.
!-----------------------------------------------------------------------
   Nrow = 1 + (N-1)/2
   Ncol = Nrow
   k = 1
   Ia(1) = 1
   Ja(1) = 1
   IF ( mod(N,2)==1 ) THEN
      A(1) = 3
   ELSE
      A(1) = 5
   ENDIF
   k = k + 1
   IF ( Ncol>1 ) THEN
      Ja(2) = 2
      A(2) = -1.0
      k = k + 1
   ENDIF
 
   DO i = 2 , Nrow
      Ia(i) = k
      Ja(k) = i - 1
      A(k) = -1.0
      k = k + 1
      Ja(k) = i
      A(k) = 6.0
      k = k + 1
      IF ( i<Nrow ) THEN
         Ja(k) = i + 1
         A(k) = -1.0
         k = k + 1
      ENDIF
   ENDDO
   Ia(Nrow+1) = k
!---------end-of-diagblk------------------------------------------------
END SUBROUTINE diagblk
!*==leftblk.f90 processed by SPAG 8.04RA 15:47 21 Feb 2024
!!SPAG Open source Personal, Educational or Academic User Clemson University  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
!-----------------------------------------------------------------------
SUBROUTINE leftblk(N,Nrow,Ncol,A,Ja,Ia)
   USE ISO_FORTRAN_ENV                 
   IMPLICIT NONE
!
! Dummy argument declarations rewritten by SPAG
!
   INTEGER , INTENT(IN) :: N
   INTEGER , INTENT(INOUT) :: Nrow
   INTEGER , INTENT(INOUT) :: Ncol
   REAL(REAL64) , INTENT(OUT) , DIMENSION(1:*) :: A
   INTEGER , INTENT(OUT) , DIMENSION(1:*) :: Ja
   INTEGER , INTENT(OUT) , DIMENSION(1:*) :: Ia
!
! Local variable declarations rewritten by SPAG
!
   INTEGER :: i , k
!
! End of declarations rewritten by SPAG
!
!-----------------------------------------------------------------------
!     Generate the subdiagonal block for the problem given.
!-----------------------------------------------------------------------
   Nrow = 1 + (N-1)/2
   Ncol = N/2
   k = 1
   DO i = 1 , Nrow
      Ia(i) = k
      IF ( Nrow/=Ncol ) THEN
         IF ( i>1 ) THEN
            Ja(k) = i - 1
            A(k) = -1.0
            k = k + 1
         ENDIF
      ENDIF
      IF ( i<=Ncol ) THEN
         Ja(k) = i
         A(k) = -1.0
         k = k + 1
      ENDIF
      IF ( Nrow==Ncol ) THEN
         IF ( i<Ncol ) THEN
            Ja(k) = i + 1
            A(k) = -1.0
            k = k + 1
         ENDIF
      ENDIF
   ENDDO
   Ia(Nrow+1) = k
!---------end-of-leftblk------------------------------------------------
END SUBROUTINE leftblk
!*==rightblk.f90 processed by SPAG 8.04RA 15:47 21 Feb 2024
!!SPAG Open source Personal, Educational or Academic User Clemson University  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
!-----------------------------------------------------------------------
SUBROUTINE rightblk(N,Nrow,Ncol,A,Ja,Ia)
   USE ISO_FORTRAN_ENV                 
   IMPLICIT NONE
!
! Dummy argument declarations rewritten by SPAG
!
   INTEGER , INTENT(IN) :: N
   INTEGER , INTENT(INOUT) :: Nrow
   INTEGER , INTENT(INOUT) :: Ncol
   REAL(REAL64) , INTENT(OUT) , DIMENSION(1:*) :: A
   INTEGER , INTENT(OUT) , DIMENSION(1:*) :: Ja
   INTEGER , INTENT(OUT) , DIMENSION(1:*) :: Ia
!
! Local variable declarations rewritten by SPAG
!
   INTEGER :: i , k
!
! End of declarations rewritten by SPAG
!
   Nrow = 1 + (N-1)/2
   Ncol = 1 + N/2
   k = 1
   DO i = 1 , Nrow
      Ia(i) = k
      IF ( Nrow==Ncol ) THEN
         IF ( i>1 ) THEN
            Ja(k) = i - 1
            A(k) = -1.0
            k = k + 1
         ENDIF
      ENDIF
      Ja(k) = i
      A(k) = -1.0
      k = k + 1
      IF ( Nrow/=Ncol ) THEN
         Ja(k) = i + 1
         A(k) = -1.0
         k = k + 1
      ENDIF
   ENDDO
   Ia(Nrow+1) = k
!---------end-of-rightblk-----------------------------------------------
END SUBROUTINE rightblk
