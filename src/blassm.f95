!*==amub.f90 processed by SPAG 8.04RA 15:47 21 Feb 2024
!!SPAG Open source Personal, Educational or Academic User Clemson University  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
!----------------------------------------------------------------------c
!                          S P A R S K I T                             c
!----------------------------------------------------------------------c
!        BASIC LINEAR ALGEBRA FOR SPARSE MATRICES. BLASSM MODULE       c
!----------------------------------------------------------------------c
! amub   :   computes     C = A*B                                      c
! aplb   :   computes     C = A+B                                      c
! aplb1  :   computes     C = A+B  [Sorted version: A, B, C sorted]    c
! aplsb  :   computes     C = A + s B                                  c
! aplsb1 :   computes     C = A+sB  [Sorted version: A, B, C sorted]   c
! apmbt  :   Computes     C = A +/- transp(B)                          c
! aplsbt :   Computes     C = A + s * transp(B)                        c
! diamua :   Computes     C = Diag * A                                 c
! amudia :   Computes     C = A* Diag                                  c
! aplsca :   Computes     A:= A + s I    (s = scalar)                  c
! apldia :   Computes     C = A + Diag.                                c
!----------------------------------------------------------------------c
! Note: this module still incomplete.                                  c
!----------------------------------------------------------------------c
SUBROUTINE amub(Nrow,Ncol,Job,A,Ja,Ia,B,Jb,Ib,C,Jc,Ic,Nzmax,Iw,Ierr)
   USE ISO_FORTRAN_ENV                 
   IMPLICIT NONE
!
! Dummy argument declarations rewritten by SPAG
!
   INTEGER , INTENT(IN) :: Nrow
   INTEGER , INTENT(IN) :: Ncol
   INTEGER , INTENT(IN) :: Job
   REAL(REAL64) , INTENT(IN) , DIMENSION(*) :: A
   INTEGER , INTENT(IN) , DIMENSION(*) :: Ja
   INTEGER , INTENT(IN) , DIMENSION(Nrow+1) :: Ia
   REAL(REAL64) , INTENT(IN) , DIMENSION(*) :: B
   INTEGER , INTENT(IN) , DIMENSION(*) :: Jb
   INTEGER , INTENT(IN) , DIMENSION(*) :: Ib
   REAL(REAL64) , INTENT(INOUT) , DIMENSION(*) :: C
   INTEGER , INTENT(INOUT) , DIMENSION(*) :: Jc
   INTEGER , INTENT(INOUT) , DIMENSION(*) :: Ic
   INTEGER , INTENT(IN) :: Nzmax
   INTEGER , INTENT(INOUT) , DIMENSION(Ncol) :: Iw
   INTEGER , INTENT(OUT) :: Ierr
!
! Local variable declarations rewritten by SPAG
!
   INTEGER :: ii , j , jcol , jj , jpos , k , ka , kb , len
   REAL(REAL64) :: scal
   LOGICAL :: values
!
! End of declarations rewritten by SPAG
!
!-----------------------------------------------------------------------
! performs the matrix by matrix product C = A B
!-----------------------------------------------------------------------
! on entry:
! ---------
! nrow  = integer. The row dimension of A = row dimension of C
! ncol  = integer. The column dimension of B = column dimension of C
! job   = integer. Job indicator. When job = 0, only the structure
!                  (i.e. the arrays jc, ic) is computed and the
!                  real values are ignored
!
! a,
! ja,
! ia   = Matrix A in compressed sparse row format.
!
! b,
! jb,
! ib    =  Matrix B in compressed sparse row format.
!
! nzmax = integer. The  length of the arrays c and jc.
!         amub will stop if the result matrix C  has a number
!         of elements that exceeds exceeds nzmax. See ierr.
!
! on return:
!----------
! c,
! jc,
! ic    = resulting matrix C in compressed sparse row sparse format.
!
! ierr  = integer. serving as error message.
!         ierr = 0 means normal return,
!         ierr .gt. 0 means that amub stopped while computing the
!         i-th row  of C with i=ierr, because the number
!         of elements in C exceeds nzmax.
!
! work arrays:
!------------
! iw    = integer work array of length equal to the number of
!         columns in A.
! Note:
!-------
!   The row dimension of B is not needed. However there is no checking
!   on the condition that ncol(A) = nrow(B).
!
!-----------------------------------------------------------------------
   values = (Job/=0)
   len = 0
   Ic(1) = 1
   Ierr = 0
!     initialize array iw.
   DO j = 1 , Ncol
      Iw(j) = 0
   ENDDO
!
   DO ii = 1 , Nrow
!     row i
      DO ka = Ia(ii) , Ia(ii+1) - 1
         IF ( values ) scal = A(ka)
         jj = Ja(ka)
         DO kb = Ib(jj) , Ib(jj+1) - 1
            jcol = Jb(kb)
            jpos = Iw(jcol)
            IF ( jpos==0 ) THEN
               len = len + 1
               IF ( len>Nzmax ) THEN
                  Ierr = ii
                  RETURN
               ENDIF
               Jc(len) = jcol
               Iw(jcol) = len
               IF ( values ) C(len) = scal*B(kb)
            ELSE
               IF ( values ) C(jpos) = C(jpos) + scal*B(kb)
            ENDIF
         ENDDO
      ENDDO
      DO k = Ic(ii) , len
         Iw(Jc(k)) = 0
      ENDDO
      Ic(ii+1) = len + 1
   ENDDO
!-------------end-of-amub-----------------------------------------------
!-----------------------------------------------------------------------
END SUBROUTINE amub
!*==aplb.f90 processed by SPAG 8.04RA 15:47 21 Feb 2024
!!SPAG Open source Personal, Educational or Academic User Clemson University  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
!-----------------------------------------------------------------------
SUBROUTINE aplb(Nrow,Ncol,Job,A,Ja,Ia,B,Jb,Ib,C,Jc,Ic,Nzmax,Iw,Ierr)
   USE ISO_FORTRAN_ENV                 
   IMPLICIT NONE
!
! Dummy argument declarations rewritten by SPAG
!
   INTEGER , INTENT(IN) :: Nrow
   INTEGER , INTENT(IN) :: Ncol
   INTEGER , INTENT(IN) :: Job
   REAL(REAL64) , INTENT(IN) , DIMENSION(*) :: A
   INTEGER , INTENT(IN) , DIMENSION(*) :: Ja
   INTEGER , INTENT(IN) , DIMENSION(Nrow+1) :: Ia
   REAL(REAL64) , INTENT(IN) , DIMENSION(*) :: B
   INTEGER , INTENT(IN) , DIMENSION(*) :: Jb
   INTEGER , INTENT(IN) , DIMENSION(Nrow+1) :: Ib
   REAL(REAL64) , INTENT(INOUT) , DIMENSION(*) :: C
   INTEGER , INTENT(INOUT) , DIMENSION(*) :: Jc
   INTEGER , INTENT(INOUT) , DIMENSION(Nrow+1) :: Ic
   INTEGER , INTENT(IN) :: Nzmax
   INTEGER , INTENT(INOUT) , DIMENSION(Ncol) :: Iw
   INTEGER , INTENT(OUT) :: Ierr
!
! Local variable declarations rewritten by SPAG
!
   INTEGER :: ii , j , jcol , jpos , k , ka , kb , len
   LOGICAL :: values
!
! End of declarations rewritten by SPAG
!
!-----------------------------------------------------------------------
! performs the matrix sum  C = A+B.
!-----------------------------------------------------------------------
! on entry:
! ---------
! nrow     = integer. The row dimension of A and B
! ncol  = integer. The column dimension of A and B.
! job   = integer. Job indicator. When job = 0, only the structure
!                  (i.e. the arrays jc, ic) is computed and the
!                  real values are ignored.
!
! a,
! ja,
! ia   = Matrix A in compressed sparse row format.
!
! b,
! jb,
! ib     =  Matrix B in compressed sparse row format.
!
! nzmax	= integer. The  length of the arrays c and jc.
!         amub will stop if the result matrix C  has a number
!         of elements that exceeds exceeds nzmax. See ierr.
!
! on return:
!----------
! c,
! jc,
! ic	= resulting matrix C in compressed sparse row sparse format.
!
! ierr	= integer. serving as error message.
!         ierr = 0 means normal return,
!         ierr .gt. 0 means that amub stopped while computing the
!         i-th row  of C with i=ierr, because the number
!         of elements in C exceeds nzmax.
!
! work arrays:
!------------
! iw	= integer work array of length equal to the number of
!         columns in A.
!
!-----------------------------------------------------------------------
   values = (Job/=0)
   Ierr = 0
   len = 0
   Ic(1) = 1
   DO j = 1 , Ncol
      Iw(j) = 0
   ENDDO
!
   DO ii = 1 , Nrow
!     row i
      DO ka = Ia(ii) , Ia(ii+1) - 1
         len = len + 1
         jcol = Ja(ka)
         IF ( len>Nzmax ) THEN
            CALL spag_block_1
            RETURN
         ENDIF
         Jc(len) = jcol
         IF ( values ) C(len) = A(ka)
         Iw(jcol) = len
      ENDDO
!
      DO kb = Ib(ii) , Ib(ii+1) - 1
         jcol = Jb(kb)
         jpos = Iw(jcol)
         IF ( jpos==0 ) THEN
            len = len + 1
            IF ( len>Nzmax ) THEN
               CALL spag_block_1
               RETURN
            ENDIF
            Jc(len) = jcol
            IF ( values ) C(len) = B(kb)
            Iw(jcol) = len
         ELSE
            IF ( values ) C(jpos) = C(jpos) + B(kb)
         ENDIF
      ENDDO
      DO k = Ic(ii) , len
         Iw(Jc(k)) = 0
      ENDDO
      Ic(ii+1) = len + 1
   ENDDO
   RETURN
CONTAINS
   SUBROUTINE spag_block_1
      Ierr = ii
   END SUBROUTINE spag_block_1
!------------end of aplb -----------------------------------------------
!-----------------------------------------------------------------------
END SUBROUTINE aplb
!*==aplb1.f90 processed by SPAG 8.04RA 15:47 21 Feb 2024
!!SPAG Open source Personal, Educational or Academic User Clemson University  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
!-----------------------------------------------------------------------
SUBROUTINE aplb1(Nrow,Ncol,Job,A,Ja,Ia,B,Jb,Ib,C,Jc,Ic,Nzmax,Ierr)
   USE ISO_FORTRAN_ENV                 
   IMPLICIT NONE
!
! Dummy argument declarations rewritten by SPAG
!
   INTEGER , INTENT(IN) :: Nrow
   INTEGER , INTENT(IN) :: Ncol
   INTEGER , INTENT(IN) :: Job
   REAL(REAL64) , INTENT(IN) , DIMENSION(*) :: A
   INTEGER , INTENT(IN) , DIMENSION(*) :: Ja
   INTEGER , INTENT(IN) , DIMENSION(Nrow+1) :: Ia
   REAL(REAL64) , INTENT(IN) , DIMENSION(*) :: B
   INTEGER , INTENT(IN) , DIMENSION(*) :: Jb
   INTEGER , INTENT(IN) , DIMENSION(Nrow+1) :: Ib
   REAL(REAL64) , INTENT(OUT) , DIMENSION(*) :: C
   INTEGER , INTENT(OUT) , DIMENSION(*) :: Jc
   INTEGER , INTENT(OUT) , DIMENSION(Nrow+1) :: Ic
   INTEGER , INTENT(IN) :: Nzmax
   INTEGER , INTENT(OUT) :: Ierr
!
! Local variable declarations rewritten by SPAG
!
   INTEGER :: i , j1 , j2 , ka , kamax , kb , kbmax , kc
   LOGICAL :: values
!
! End of declarations rewritten by SPAG
!
!-----------------------------------------------------------------------
! performs the matrix sum  C = A+B for matrices in sorted CSR format.
! the difference with aplb  is that the resulting matrix is such that
! the elements of each row are sorted with increasing column indices in
! each row, provided the original matrices are sorted in the same way.
!-----------------------------------------------------------------------
! on entry:
! ---------
! nrow	= integer. The row dimension of A and B
! ncol  = integer. The column dimension of A and B.
! job   = integer. Job indicator. When job = 0, only the structure
!                  (i.e. the arrays jc, ic) is computed and the
!                  real values are ignored.
!
! a,
! ja,
! ia   = Matrix A in compressed sparse row format with entries sorted
!
! b,
! jb,
! ib	=  Matrix B in compressed sparse row format with entries sorted
!        ascendly in each row
!
! nzmax	= integer. The  length of the arrays c and jc.
!         amub will stop if the result matrix C  has a number
!         of elements that exceeds exceeds nzmax. See ierr.
!
! on return:
!----------
! c,
! jc,
! ic	= resulting matrix C in compressed sparse row sparse format
!         with entries sorted ascendly in each row.
!
! ierr	= integer. serving as error message.
!         ierr = 0 means normal return,
!         ierr .gt. 0 means that amub stopped while computing the
!         i-th row  of C with i=ierr, because the number
!         of elements in C exceeds nzmax.
!
! Notes:
!-------
!     this will not work if any of the two input matrices is not sorted
!-----------------------------------------------------------------------
   INTEGER :: spag_nextblock_1
   spag_nextblock_1 = 1
   SPAG_DispatchLoop_1: DO
      SELECT CASE (spag_nextblock_1)
      CASE (1)
         values = (Job/=0)
         Ierr = 0
         kc = 1
         Ic(1) = kc
!
         DO i = 1 , Nrow
            ka = Ia(i)
            kb = Ib(i)
            kamax = Ia(i+1) - 1
            kbmax = Ib(i+1) - 1
            SPAG_Loop_3_1: DO
               IF ( ka<=kamax ) THEN
                  j1 = Ja(ka)
               ELSE
                  j1 = Ncol + 1
               ENDIF
               IF ( kb<=kbmax ) THEN
                  j2 = Jb(kb)
               ELSE
                  j2 = Ncol + 1
               ENDIF
!
!     three cases
!
               IF ( kc>Nzmax ) THEN
                  spag_nextblock_1 = 2
                  CYCLE SPAG_DispatchLoop_1
               ENDIF
               IF ( j1==j2 ) THEN
                  IF ( values ) C(kc) = A(ka) + B(kb)
                  Jc(kc) = j1
                  ka = ka + 1
                  kb = kb + 1
                  kc = kc + 1
               ELSEIF ( j1<j2 ) THEN
                  Jc(kc) = j1
                  IF ( values ) C(kc) = A(ka)
                  ka = ka + 1
                  kc = kc + 1
               ELSEIF ( j1>j2 ) THEN
                  Jc(kc) = j2
                  IF ( values ) C(kc) = B(kb)
                  kb = kb + 1
                  kc = kc + 1
               ENDIF
               IF ( ka>kamax .AND. kb>kbmax ) THEN
                  Ic(i+1) = kc
                  EXIT SPAG_Loop_3_1
               ENDIF
            ENDDO SPAG_Loop_3_1
         ENDDO
         RETURN
      CASE (2)
         Ierr = i
         EXIT SPAG_DispatchLoop_1
      END SELECT
   ENDDO SPAG_DispatchLoop_1
!------------end-of-aplb1-----------------------------------------------
!-----------------------------------------------------------------------
END SUBROUTINE aplb1
!*==aplsb.f90 processed by SPAG 8.04RA 15:47 21 Feb 2024
!!SPAG Open source Personal, Educational or Academic User Clemson University  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
!-----------------------------------------------------------------------
SUBROUTINE aplsb(Nrow,Ncol,A,Ja,Ia,S,B,Jb,Ib,C,Jc,Ic,Nzmax,Iw,Ierr)
   USE ISO_FORTRAN_ENV                 
   IMPLICIT NONE
!
! Dummy argument declarations rewritten by SPAG
!
   INTEGER , INTENT(IN) :: Nrow
   INTEGER , INTENT(IN) :: Ncol
   REAL(REAL64) , INTENT(IN) , DIMENSION(*) :: A
   INTEGER , INTENT(IN) , DIMENSION(*) :: Ja
   INTEGER , INTENT(IN) , DIMENSION(Nrow+1) :: Ia
   REAL(REAL64) , INTENT(IN) :: S
   REAL(REAL64) , INTENT(IN) , DIMENSION(*) :: B
   INTEGER , INTENT(IN) , DIMENSION(*) :: Jb
   INTEGER , INTENT(IN) , DIMENSION(Nrow+1) :: Ib
   REAL(REAL64) , INTENT(INOUT) , DIMENSION(*) :: C
   INTEGER , INTENT(INOUT) , DIMENSION(*) :: Jc
   INTEGER , INTENT(INOUT) , DIMENSION(Nrow+1) :: Ic
   INTEGER , INTENT(IN) :: Nzmax
   INTEGER , INTENT(INOUT) , DIMENSION(Ncol) :: Iw
   INTEGER , INTENT(OUT) :: Ierr
!
! Local variable declarations rewritten by SPAG
!
   INTEGER :: ii , j , jcol , jpos , k , ka , kb , len
!
! End of declarations rewritten by SPAG
!
!-----------------------------------------------------------------------
! performs the matrix sum  C = A+s*B.
!-----------------------------------------------------------------------
! on entry:
! ---------
! nrow     = integer. The row dimension of A and B
! ncol  = integer. The column dimension of A and B.
! job   = integer. Job indicator. When job = 0, only the structure
!                  (i.e. the arrays jc, ic) is computed and the
!                  real values are ignored.
!
! a,
! ja,
! ia   = Matrix A in compressed sparse row format.
!
! s    = real*8 - coefficient that multiplies B.
! b,
! jb,
! ib     =  Matrix B in compressed sparse row format.
!
! nzmax	= integer. The  length of the arrays c and jc.
!         amub will stop if the result matrix C  has a number
!         of elements that exceeds exceeds nzmax. See ierr.
!
! on return:
!----------
! c,
! jc,
! ic	= resulting matrix C in compressed sparse row sparse format.
!
! ierr	= integer. serving as error message.
!         ierr = 0 means normal return,
!         ierr .gt. 0 means that aplsb1 stopped while computing the
!         i-th row  of C with i=ierr, because the number
!         of elements in C exceeds nzmax.
!
! work arrays:
!------------
! iw	= integer work array of length equal to the number of
!         columns in A.
! note: expanded  row implementation. Does not require column indices to
!       be sorted.
!-----------------------------------------------------------------------
   Ierr = 0
   len = 0
   Ic(1) = 1
   DO j = 1 , Ncol
      Iw(j) = 0
   ENDDO
!
   DO ii = 1 , Nrow
!     copy row ii  to C
      DO ka = Ia(ii) , Ia(ii+1) - 1
         len = len + 1
         jcol = Ja(ka)
         IF ( len>Nzmax ) THEN
            CALL spag_block_1
            RETURN
         ENDIF
         Jc(len) = jcol
         C(len) = A(ka)
         Iw(jcol) = len
      ENDDO
!
      DO kb = Ib(ii) , Ib(ii+1) - 1
         jcol = Jb(kb)
         jpos = Iw(jcol)
         IF ( jpos==0 ) THEN
            len = len + 1
            IF ( len>Nzmax ) THEN
               CALL spag_block_1
               RETURN
            ENDIF
            Jc(len) = jcol
            C(len) = S*B(kb)
            Iw(jcol) = len
         ELSE
            C(jpos) = C(jpos) + S*B(kb)
         ENDIF
      ENDDO
      DO k = Ic(ii) , len
         Iw(Jc(k)) = 0
      ENDDO
      Ic(ii+1) = len + 1
   ENDDO
   RETURN
CONTAINS
   SUBROUTINE spag_block_1
      Ierr = ii
   END SUBROUTINE spag_block_1
!------------end of aplsb ----------------------------------------------
!-----------------------------------------------------------------------
END SUBROUTINE aplsb
!*==aplsb1.f90 processed by SPAG 8.04RA 15:47 21 Feb 2024
!!SPAG Open source Personal, Educational or Academic User Clemson University  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
!-----------------------------------------------------------------------
SUBROUTINE aplsb1(Nrow,Ncol,A,Ja,Ia,S,B,Jb,Ib,C,Jc,Ic,Nzmax,Ierr)
   USE ISO_FORTRAN_ENV                 
   IMPLICIT NONE
!
! Dummy argument declarations rewritten by SPAG
!
   INTEGER , INTENT(IN) :: Nrow
   INTEGER , INTENT(IN) :: Ncol
   REAL(REAL64) , INTENT(IN) , DIMENSION(*) :: A
   INTEGER , INTENT(IN) , DIMENSION(*) :: Ja
   INTEGER , INTENT(IN) , DIMENSION(Nrow+1) :: Ia
   REAL(REAL64) , INTENT(IN) :: S
   REAL(REAL64) , INTENT(IN) , DIMENSION(*) :: B
   INTEGER , INTENT(IN) , DIMENSION(*) :: Jb
   INTEGER , INTENT(IN) , DIMENSION(Nrow+1) :: Ib
   REAL(REAL64) , INTENT(OUT) , DIMENSION(*) :: C
   INTEGER , INTENT(OUT) , DIMENSION(*) :: Jc
   INTEGER , INTENT(OUT) , DIMENSION(Nrow+1) :: Ic
   INTEGER , INTENT(IN) :: Nzmax
   INTEGER , INTENT(OUT) :: Ierr
!
! Local variable declarations rewritten by SPAG
!
   INTEGER :: i , j1 , j2 , ka , kamax , kb , kbmax , kc
!
! End of declarations rewritten by SPAG
!
   INTEGER :: spag_nextblock_1
   spag_nextblock_1 = 1
   SPAG_DispatchLoop_1: DO
      SELECT CASE (spag_nextblock_1)
      CASE (1)
!-----------------------------------------------------------------------
! performs the operation C = A+s B for matrices in sorted CSR format.
! the difference with aplsb is that the resulting matrix is such that
! the elements of each row are sorted with increasing column indices in
! each row, provided the original matrices are sorted in the same way.
!-----------------------------------------------------------------------
! on entry:
! ---------
! nrow	= integer. The row dimension of A and B
! ncol  = integer. The column dimension of A and B.
!
! a,
! ja,
! ia   = Matrix A in compressed sparse row format with entries sorted
!
! s	= real. scalar factor for B.
!
! b,
! jb,
! ib	=  Matrix B in compressed sparse row format with entries sorted
!        ascendly in each row
!
! nzmax	= integer. The  length of the arrays c and jc.
!         amub will stop if the result matrix C  has a number
!         of elements that exceeds exceeds nzmax. See ierr.
!
! on return:
!----------
! c,
! jc,
! ic	= resulting matrix C in compressed sparse row sparse format
!         with entries sorted ascendly in each row.
!
! ierr	= integer. serving as error message.
!         ierr = 0 means normal return,
!         ierr .gt. 0 means that amub stopped while computing the
!         i-th row  of C with i=ierr, because the number
!         of elements in C exceeds nzmax.
!
! Notes:
!-------
!     this will not work if any of the two input matrices is not sorted
!-----------------------------------------------------------------------
         Ierr = 0
         kc = 1
         Ic(1) = kc
!
!     the following loop does a merge of two sparse rows + adds  them.
!
         DO i = 1 , Nrow
            ka = Ia(i)
            kb = Ib(i)
            kamax = Ia(i+1) - 1
            kbmax = Ib(i+1) - 1
            SPAG_Loop_3_1: DO
!
!     this is a while  -- do loop --
!
               IF ( ka<=kamax .OR. kb<=kbmax ) THEN
!
                  IF ( ka<=kamax ) THEN
                     j1 = Ja(ka)
                  ELSE
!     take j1 large enough  that always j2 .lt. j1
                     j1 = Ncol + 1
                  ENDIF
                  IF ( kb<=kbmax ) THEN
                     j2 = Jb(kb)
                  ELSE
!     similarly take j2 large enough  that always j1 .lt. j2
                     j2 = Ncol + 1
                  ENDIF
!
!     three cases
!
                  IF ( j1==j2 ) THEN
                     C(kc) = A(ka) + S*B(kb)
                     Jc(kc) = j1
                     ka = ka + 1
                     kb = kb + 1
                     kc = kc + 1
                  ELSEIF ( j1<j2 ) THEN
                     Jc(kc) = j1
                     C(kc) = A(ka)
                     ka = ka + 1
                     kc = kc + 1
                  ELSEIF ( j1>j2 ) THEN
                     Jc(kc) = j2
                     C(kc) = S*B(kb)
                     kb = kb + 1
                     kc = kc + 1
                  ENDIF
                  IF ( kc<=Nzmax ) CYCLE
                  spag_nextblock_1 = 2
                  CYCLE SPAG_DispatchLoop_1
!
!     end while loop
!
               ENDIF
               Ic(i+1) = kc
               EXIT SPAG_Loop_3_1
            ENDDO SPAG_Loop_3_1
         ENDDO
         RETURN
      CASE (2)
         Ierr = i
         EXIT SPAG_DispatchLoop_1
      END SELECT
   ENDDO SPAG_DispatchLoop_1
!------------end-of-aplsb1 ---------------------------------------------
!-----------------------------------------------------------------------
END SUBROUTINE aplsb1
!*==apmbt.f90 processed by SPAG 8.04RA 15:47 21 Feb 2024
!!SPAG Open source Personal, Educational or Academic User Clemson University  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
!-----------------------------------------------------------------------
SUBROUTINE apmbt(Nrow,Ncol,Job,A,Ja,Ia,B,Jb,Ib,C,Jc,Ic,Nzmax,Iw,Ierr)
   USE ISO_FORTRAN_ENV                 
   IMPLICIT NONE
!
! Dummy argument declarations rewritten by SPAG
!
   INTEGER :: Nrow
   INTEGER :: Ncol
   INTEGER , INTENT(IN) :: Job
   REAL(REAL64) , INTENT(IN) , DIMENSION(*) :: A
   INTEGER , INTENT(IN) , DIMENSION(*) :: Ja
   INTEGER , INTENT(IN) , DIMENSION(Nrow+1) :: Ia
   REAL(REAL64) , DIMENSION(*) :: B
   INTEGER , DIMENSION(*) :: Jb
   INTEGER , DIMENSION(Ncol+1) :: Ib
   REAL(REAL64) , INTENT(INOUT) , DIMENSION(*) :: C
   INTEGER , INTENT(INOUT) , DIMENSION(*) :: Jc
   INTEGER , INTENT(INOUT) , DIMENSION(*) :: Ic
   INTEGER , INTENT(IN) :: Nzmax
   INTEGER , INTENT(INOUT) , DIMENSION(*) :: Iw
   INTEGER :: Ierr
!
! Local variable declarations rewritten by SPAG
!
   INTEGER :: i , ii , ipos , j , jcol , jpos , k , ka , len , ljob , nnza , nnzb
   LOGICAL :: values
   EXTERNAL coicsr , csrcoo , csrcsc
!
! End of declarations rewritten by SPAG
!
!-----------------------------------------------------------------------
! performs the matrix sum  C = A + transp(B) or C = A - transp(B)
!-----------------------------------------------------------------------
! on entry:
! ---------
! nrow	= integer. The row dimension of A and transp(B)
! ncol  = integer. The column dimension of A. Also the row
!                  dimension of B.
!
! job	= integer. if job = -1, apmbt will compute C= A - transp(B)
!         (structure + values)
!         if (job .eq. 1)  it will compute C=A+transp(A)
!         (structure+ values)
!         if (job .eq. 0) it will compute the structure of
!         C= A+/-transp(B) only (ignoring all real values).
!         any other value of job will be treated as  job=1
! a,
! ja,
! ia    = Matrix A in compressed sparse row format.
!
! b,
! jb,
! ib	=  Matrix B in compressed sparse row format.
!
! nzmax	= integer. The  length of the arrays c, jc, and ic.
!         amub will stop if the result matrix C  has a number
!         of elements that exceeds exceeds nzmax. See ierr.
!
! on return:
!----------
! c,
! jc,
! ic	= resulting matrix C in compressed sparse row format.
!
! ierr	= integer. serving as error message.
!         ierr = 0 means normal return.
!         ierr = -1 means that nzmax was .lt. either the number of
!         nonzero elements of A or the number of nonzero elements in B.
!         ierr .gt. 0 means that amub stopped while computing the
!         i-th row  of C with i=ierr, because the number
!         of elements in C exceeds nzmax.
!
! work arrays:
!------------
! iw	= integer work array of length at least max(ncol,nrow)
!
! Notes:
!------- It is important to note that here all of three arrays c, ic,
!        and jc are assumed to be of length nnz(c). This is because
!        the matrix is internally converted in coordinate format.
!
!-----------------------------------------------------------------------
   values = (Job/=0)
!
   Ierr = 0
   DO j = 1 , Ncol
      Iw(j) = 0
   ENDDO
!
   nnza = Ia(Nrow+1) - 1
   nnzb = Ib(Ncol+1) - 1
   len = nnzb
   IF ( Nzmax<nnzb .OR. Nzmax<nnza ) THEN
      Ierr = -1
      RETURN
   ENDIF
!
! trasnpose matrix b into c
!
   ljob = 0
   IF ( values ) ljob = 1
   ipos = 1
   CALL csrcsc(Ncol,ljob,ipos,B,Jb,Ib,C,Jc,Ic)
!-----------------------------------------------------------------------
   IF ( Job==-1 ) THEN
      DO k = 1 , len
         C(k) = -C(k)
      ENDDO
   ENDIF
!
!--------------- main loop --------------------------------------------
!
   DO ii = 1 , Nrow
      DO k = Ic(ii) , Ic(ii+1) - 1
         Iw(Jc(k)) = k
      ENDDO
!-----------------------------------------------------------------------
      DO ka = Ia(ii) , Ia(ii+1) - 1
         jcol = Ja(ka)
         jpos = Iw(jcol)
         IF ( jpos==0 ) THEN
!
!     if fill-in append in coordinate format to matrix.
!
            len = len + 1
            IF ( len>Nzmax ) THEN
               CALL spag_block_1
               RETURN
            ENDIF
            Jc(len) = jcol
 
            Ic(len) = ii
            IF ( values ) C(len) = A(ka)
         ELSE
!     else do addition.
            IF ( values ) C(jpos) = C(jpos) + A(ka)
         ENDIF
      ENDDO
      DO k = Ic(ii) , Ic(ii+1) - 1
         Iw(Jc(k)) = 0
      ENDDO
   ENDDO
!
!     convert first part of matrix (without fill-ins) into coo format
!
   ljob = 2
   IF ( values ) ljob = 3
   DO i = 1 , Nrow + 1
      Iw(i) = Ic(i)
   ENDDO
   CALL csrcoo(Nrow,ljob,nnzb,C,Jc,Iw,nnzb,C,Ic,Jc,Ierr)
!
!     convert the whole thing back to csr format.
!
   ljob = 0
   IF ( values ) ljob = 1
   CALL coicsr(Nrow,len,ljob,C,Jc,Ic,Iw)
   RETURN
CONTAINS
   SUBROUTINE spag_block_1
      Ierr = ii
   END SUBROUTINE spag_block_1
!--------end-of-apmbt---------------------------------------------------
!-----------------------------------------------------------------------
END SUBROUTINE apmbt
!*==aplsbt.f90 processed by SPAG 8.04RA 15:47 21 Feb 2024
!!SPAG Open source Personal, Educational or Academic User Clemson University  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
!-----------------------------------------------------------------------
SUBROUTINE aplsbt(Nrow,Ncol,A,Ja,Ia,S,B,Jb,Ib,C,Jc,Ic,Nzmax,Iw,Ierr)
   USE ISO_FORTRAN_ENV                 
   IMPLICIT NONE
!
! Dummy argument declarations rewritten by SPAG
!
   INTEGER :: Nrow
   INTEGER :: Ncol
   REAL(REAL64) , INTENT(IN) , DIMENSION(*) :: A
   INTEGER , INTENT(IN) , DIMENSION(*) :: Ja
   INTEGER , INTENT(IN) , DIMENSION(Nrow+1) :: Ia
   REAL(REAL64) , INTENT(IN) :: S
   REAL(REAL64) , DIMENSION(*) :: B
   INTEGER , DIMENSION(*) :: Jb
   INTEGER , DIMENSION(Ncol+1) :: Ib
   REAL(REAL64) , INTENT(INOUT) , DIMENSION(*) :: C
   INTEGER , INTENT(INOUT) , DIMENSION(*) :: Jc
   INTEGER , INTENT(INOUT) , DIMENSION(*) :: Ic
   INTEGER , INTENT(IN) :: Nzmax
   INTEGER , INTENT(INOUT) , DIMENSION(*) :: Iw
   INTEGER :: Ierr
!
! Local variable declarations rewritten by SPAG
!
   INTEGER :: i , ii , ipos , j , jcol , jpos , k , ka , len , ljob , nnza , nnzb
   EXTERNAL coicsr , csrcoo , csrcsc
!
! End of declarations rewritten by SPAG
!
!-----------------------------------------------------------------------
! performs the matrix sum  C = A + transp(B).
!-----------------------------------------------------------------------
! on entry:
! ---------
! nrow	= integer. The row dimension of A and transp(B)
! ncol  = integer. The column dimension of A. Also the row
!                  dimension of B.
!
! a,
! ja,
! ia    = Matrix A in compressed sparse row format.
!
! s	= real. scalar factor for B.
!
!
! b,
! jb,
! ib	=  Matrix B in compressed sparse row format.
!
! nzmax	= integer. The  length of the arrays c, jc, and ic.
!         amub will stop if the result matrix C  has a number
!         of elements that exceeds exceeds nzmax. See ierr.
!
! on return:
!----------
! c,
! jc,
! ic	= resulting matrix C in compressed sparse row format.
!
! ierr	= integer. serving as error message.
!         ierr = 0 means normal return.
!         ierr = -1 means that nzmax was .lt. either the number of
!         nonzero elements of A or the number of nonzero elements in B.
!         ierr .gt. 0 means that amub stopped while computing the
!         i-th row  of C with i=ierr, because the number
!         of elements in C exceeds nzmax.
!
! work arrays:
!------------
! iw	= integer work array of length at least max(nrow,ncol)
!
! Notes:
!------- It is important to note that here all of three arrays c, ic,
!        and jc are assumed to be of length nnz(c). This is because
!        the matrix is internally converted in coordinate format.
!
!-----------------------------------------------------------------------
   Ierr = 0
   DO j = 1 , Ncol
      Iw(j) = 0
   ENDDO
!
   nnza = Ia(Nrow+1) - 1
   nnzb = Ib(Ncol+1) - 1
   len = nnzb
   IF ( Nzmax<nnzb .OR. Nzmax<nnza ) THEN
      Ierr = -1
      RETURN
   ENDIF
!
!     transpose matrix b into c
!
   ljob = 1
   ipos = 1
   CALL csrcsc(Ncol,ljob,ipos,B,Jb,Ib,C,Jc,Ic)
   DO k = 1 , len
      C(k) = C(k)*S
   ENDDO
!
!     main loop. add rows from ii = 1 to nrow.
!
   DO ii = 1 , Nrow
!     iw is used as a system to recognize whether there
!     was a nonzero element in c.
      DO k = Ic(ii) , Ic(ii+1) - 1
         Iw(Jc(k)) = k
      ENDDO
!
      DO ka = Ia(ii) , Ia(ii+1) - 1
         jcol = Ja(ka)
         jpos = Iw(jcol)
         IF ( jpos==0 ) THEN
!
!     if fill-in append in coordinate format to matrix.
!
            len = len + 1
            IF ( len>Nzmax ) THEN
               CALL spag_block_1
               RETURN
            ENDIF
            Jc(len) = jcol
            Ic(len) = ii
            C(len) = A(ka)
         ELSE
!     else do addition.
            C(jpos) = C(jpos) + A(ka)
         ENDIF
      ENDDO
      DO k = Ic(ii) , Ic(ii+1) - 1
         Iw(Jc(k)) = 0
      ENDDO
   ENDDO
!
!     convert first part of matrix (without fill-ins) into coo format
!
   ljob = 3
   DO i = 1 , Nrow + 1
      Iw(i) = Ic(i)
   ENDDO
   CALL csrcoo(Nrow,ljob,nnzb,C,Jc,Iw,nnzb,C,Ic,Jc,Ierr)
!
!     convert the whole thing back to csr format.
!
   ljob = 1
   CALL coicsr(Nrow,len,ljob,C,Jc,Ic,Iw)
   RETURN
CONTAINS
   SUBROUTINE spag_block_1
      Ierr = ii
   END SUBROUTINE spag_block_1
!--------end-of-aplsbt--------------------------------------------------
!-----------------------------------------------------------------------
END SUBROUTINE aplsbt
!*==diamua.f90 processed by SPAG 8.04RA 15:47 21 Feb 2024
!!SPAG Open source Personal, Educational or Academic User Clemson University  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
!-----------------------------------------------------------------------
SUBROUTINE diamua(Nrow,Job,A,Ja,Ia,Diag,B,Jb,Ib)
   USE ISO_FORTRAN_ENV                 
   IMPLICIT NONE
!
! Dummy argument declarations rewritten by SPAG
!
   INTEGER , INTENT(IN) :: Nrow
   INTEGER , INTENT(IN) :: Job
   REAL(REAL64) , INTENT(IN) , DIMENSION(*) :: A
   INTEGER , INTENT(IN) , DIMENSION(*) :: Ja
   INTEGER , INTENT(IN) , DIMENSION(Nrow+1) :: Ia
   REAL(REAL64) , INTENT(IN) , DIMENSION(Nrow) :: Diag
   REAL(REAL64) , INTENT(OUT) , DIMENSION(*) :: B
   INTEGER , INTENT(OUT) , DIMENSION(*) :: Jb
   INTEGER , INTENT(OUT) , DIMENSION(Nrow+1) :: Ib
!
! Local variable declarations rewritten by SPAG
!
   INTEGER :: ii , k , k1 , k2
   REAL(REAL64) :: scal
!
! End of declarations rewritten by SPAG
!
!-----------------------------------------------------------------------
! performs the matrix by matrix product B = Diag * A  (in place)
!-----------------------------------------------------------------------
! on entry:
! ---------
! nrow	= integer. The row dimension of A
!
! job   = integer. job indicator. Job=0 means get array b only
!         job = 1 means get b, and the integer arrays ib, jb.
!
! a,
! ja,
! ia   = Matrix A in compressed sparse row format.
!
! diag = diagonal matrix stored as a vector dig(1:n)
!
! on return:
!----------
!
! b,
! jb,
! ib	= resulting matrix B in compressed sparse row sparse format.
!
! Notes:
!-------
! 1)        The column dimension of A is not needed.
! 2)        algorithm in place (B can take the place of A).
!           in this case use job=0.
!-----------------------------------------------------------------
   DO ii = 1 , Nrow
!
!     normalize each row
!
      k1 = Ia(ii)
      k2 = Ia(ii+1) - 1
      scal = Diag(ii)
      DO k = k1 , k2
         B(k) = A(k)*scal
      ENDDO
   ENDDO
!
   IF ( Job==0 ) RETURN
!
   DO ii = 1 , Nrow + 1
      Ib(ii) = Ia(ii)
   ENDDO
   DO k = Ia(1) , Ia(Nrow+1) - 1
      Jb(k) = Ja(k)
   ENDDO
!----------end-of-diamua------------------------------------------------
!-----------------------------------------------------------------------
END SUBROUTINE diamua
!*==amudia.f90 processed by SPAG 8.04RA 15:47 21 Feb 2024
!!SPAG Open source Personal, Educational or Academic User Clemson University  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
!-----------------------------------------------------------------------
SUBROUTINE amudia(Nrow,Job,A,Ja,Ia,Diag,B,Jb,Ib)
   USE ISO_FORTRAN_ENV                 
   IMPLICIT NONE
!
! Dummy argument declarations rewritten by SPAG
!
   INTEGER , INTENT(IN) :: Nrow
   INTEGER , INTENT(IN) :: Job
   REAL(REAL64) , INTENT(IN) , DIMENSION(*) :: A
   INTEGER , INTENT(IN) , DIMENSION(*) :: Ja
   INTEGER , INTENT(IN) , DIMENSION(Nrow+1) :: Ia
   REAL(REAL64) , INTENT(IN) , DIMENSION(Nrow) :: Diag
   REAL(REAL64) , INTENT(OUT) , DIMENSION(*) :: B
   INTEGER , INTENT(OUT) , DIMENSION(*) :: Jb
   INTEGER , INTENT(OUT) , DIMENSION(Nrow+1) :: Ib
!
! Local variable declarations rewritten by SPAG
!
   INTEGER :: ii , k , k1 , k2
!
! End of declarations rewritten by SPAG
!
!-----------------------------------------------------------------------
! performs the matrix by matrix product B = A * Diag  (in place)
!-----------------------------------------------------------------------
! on entry:
! ---------
! nrow	= integer. The row dimension of A
!
! job   = integer. job indicator. Job=0 means get array b only
!         job = 1 means get b, and the integer arrays ib, jb.
!
! a,
! ja,
! ia   = Matrix A in compressed sparse row format.
!
! diag = diagonal matrix stored as a vector dig(1:n)
!
! on return:
!----------
!
! b,
! jb,
! ib	= resulting matrix B in compressed sparse row sparse format.
!
! Notes:
!-------
! 1)        The column dimension of A is not needed.
! 2)        algorithm in place (B can take the place of A).
!-----------------------------------------------------------------
   DO ii = 1 , Nrow
!
!     scale each element
!
      k1 = Ia(ii)
      k2 = Ia(ii+1) - 1
      DO k = k1 , k2
         B(k) = A(k)*Diag(Ja(k))
      ENDDO
   ENDDO
!
   IF ( Job==0 ) RETURN
!
   DO ii = 1 , Nrow + 1
      Ib(ii) = Ia(ii)
   ENDDO
   DO k = Ia(1) , Ia(Nrow+1) - 1
      Jb(k) = Ja(k)
   ENDDO
!-----------------------------------------------------------------------
!-----------end-of-amudiag----------------------------------------------
END SUBROUTINE amudia
!*==aplsca.f90 processed by SPAG 8.04RA 15:47 21 Feb 2024
!!SPAG Open source Personal, Educational or Academic User Clemson University  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
!-----------------------------------------------------------------------
SUBROUTINE aplsca(Nrow,A,Ja,Ia,Scal,Iw)
   USE ISO_FORTRAN_ENV                 
   IMPLICIT NONE
!
! Dummy argument declarations rewritten by SPAG
!
   INTEGER :: Nrow
   REAL(REAL64) , INTENT(INOUT) , DIMENSION(*) :: A
   INTEGER , INTENT(INOUT) , DIMENSION(*) :: Ja
   INTEGER , INTENT(INOUT) , DIMENSION(Nrow+1) :: Ia
   REAL(REAL64) , INTENT(IN) :: Scal
   INTEGER , INTENT(INOUT) , DIMENSION(*) :: Iw
!
! Local variable declarations rewritten by SPAG
!
   INTEGER :: icount , ii , j , k , k1 , k2 , ko
   LOGICAL :: test
   EXTERNAL diapos
!
! End of declarations rewritten by SPAG
!
!-----------------------------------------------------------------------
! Adds a scalar to the diagonal entries of a sparse matrix A :=A + s I
!-----------------------------------------------------------------------
! on entry:
! ---------
! nrow	= integer. The row dimension of A
!
! a,
! ja,
! ia    = Matrix A in compressed sparse row format.
!
! scal  = real. scalar to add to the diagonal entries.
!
! on return:
!----------
!
! a,
! ja,
! ia	= matrix A with diagonal elements shifted (or created).
!
! iw    = integer work array of length n. On return iw will
!         contain  the positions of the diagonal entries in the
!         output matrix. (i.e., a(iw(k)), ja(iw(k)), k=1,...n,
!         are the values/column indices of the diagonal elements
!         of the output matrix. ).
!
! Notes:
!-------
!     The column dimension of A is not needed.
!     important: the matrix a may be expanded slightly to allow for
!     additions of nonzero elements to previously nonexisting diagonals.
!     The is no checking as to whether there is enough space appended
!     to the arrays a and ja. if not sure allow for n additional
!     elemnts.
!     coded by Y. Saad. Latest version July, 19, 1990
!-----------------------------------------------------------------------
!
   CALL diapos(Nrow,Ja,Ia,Iw)
   icount = 0
   DO j = 1 , Nrow
      IF ( Iw(j)==0 ) THEN
         icount = icount + 1
      ELSE
         A(Iw(j)) = A(Iw(j)) + Scal
      ENDIF
   ENDDO
!
!     if no diagonal elements to insert in data structure return.
!
   IF ( icount==0 ) RETURN
!
! shift the nonzero elements if needed, to allow for created
! diagonal elements.
!
   ko = Ia(Nrow+1) + icount
!
!     copy rows backward
!
   DO ii = Nrow , 1 , -1
!
!     go through  row ii
!
      k1 = Ia(ii)
      k2 = Ia(ii+1) - 1
      Ia(ii+1) = ko
      test = (Iw(ii)==0)
      DO k = k2 , k1 , -1
         j = Ja(k)
         IF ( test .AND. (j<ii) ) THEN
            test = .FALSE.
            ko = ko - 1
            A(ko) = Scal
            Ja(ko) = ii
            Iw(ii) = ko
         ENDIF
         ko = ko - 1
         A(ko) = A(k)
         Ja(ko) = j
      ENDDO
!     diagonal element has not been added yet.
      IF ( test ) THEN
         ko = ko - 1
         A(ko) = Scal
         Ja(ko) = ii
         Iw(ii) = ko
      ENDIF
   ENDDO
   Ia(1) = ko
!-----------------------------------------------------------------------
!----------end-of-aplsca------------------------------------------------
END SUBROUTINE aplsca
!*==apldia.f90 processed by SPAG 8.04RA 15:47 21 Feb 2024
!!SPAG Open source Personal, Educational or Academic User Clemson University  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
!-----------------------------------------------------------------------
SUBROUTINE apldia(Nrow,Job,A,Ja,Ia,Diag,B,Jb,Ib,Iw)
   USE ISO_FORTRAN_ENV                 
   IMPLICIT NONE
!
! Dummy argument declarations rewritten by SPAG
!
   INTEGER :: Nrow
   INTEGER , INTENT(IN) :: Job
   REAL(REAL64) , INTENT(IN) , DIMENSION(*) :: A
   INTEGER , DIMENSION(*) :: Ja
   INTEGER , DIMENSION(Nrow+1) :: Ia
   REAL(REAL64) , INTENT(IN) , DIMENSION(Nrow) :: Diag
   REAL(REAL64) , INTENT(OUT) , DIMENSION(*) :: B
   INTEGER , INTENT(INOUT) , DIMENSION(*) :: Jb
   INTEGER , INTENT(INOUT) , DIMENSION(Nrow+1) :: Ib
   INTEGER , INTENT(INOUT) , DIMENSION(*) :: Iw
!
! Local variable declarations rewritten by SPAG
!
   INTEGER :: icount , ii , j , k , k1 , k2 , ko , nnz
   LOGICAL :: test
   EXTERNAL diapos
!
! End of declarations rewritten by SPAG
!
!-----------------------------------------------------------------------
! Adds a diagonal matrix to a general sparse matrix:  B = A + Diag
!-----------------------------------------------------------------------
! on entry:
! ---------
! nrow	= integer. The row dimension of A
!
! job   = integer. job indicator. Job=0 means get array b only
!         (i.e. assume that a has already been copied into array b,
!         or that algorithm is used in place. ) For all practical
!         purposes enter job=0 for an in-place call and job=1 otherwise
!
!         Note: in case there are missing diagonal elements in A,
!         then the option job =0 will be ignored, since the algorithm
!         must modify the data structure (i.e. jb, ib) in this
!         situation.
!
! a,
! ja,
! ia   = Matrix A in compressed sparse row format.
!
! diag = diagonal matrix stored as a vector dig(1:n)
!
! on return:
!----------
!
! b,
! jb,
! ib	= resulting matrix B in compressed sparse row sparse format.
!
!
! iw    = integer work array of length n. On return iw will
!         contain  the positions of the diagonal entries in the
!         output matrix. (i.e., a(iw(k)), ja(iw(k)), k=1,...n,
!         are the values/column indices of the diagonal elements
!         of the output matrix. ).
!
! Notes:
!-------
! 1)        The column dimension of A is not needed.
! 2)        algorithm in place (b, jb, ib, can be the same as
!           a, ja, ia, on entry). See comments for parameter job.
!
! coded by Y. Saad. Latest version July, 19, 1990
!-----------------------------------------------------------------
!
!     copy integer arrays into b's data structure if required
!
   IF ( Job/=0 ) THEN
      nnz = Ia(Nrow+1) - 1
      DO k = 1 , nnz
         Jb(k) = Ja(k)
         B(k) = A(k)
      ENDDO
      DO k = 1 , Nrow + 1
         Ib(k) = Ia(k)
      ENDDO
   ENDIF
!
!     get positions of diagonal elements in data structure.
!
   CALL diapos(Nrow,Ja,Ia,Iw)
!
!     count number of holes in diagonal and add diag(*) elements to
!     valid diagonal entries.
!
   icount = 0
   DO j = 1 , Nrow
      IF ( Iw(j)==0 ) THEN
         icount = icount + 1
      ELSE
         B(Iw(j)) = A(Iw(j)) + Diag(j)
      ENDIF
   ENDDO
!
!     if no diagonal elements to insert return
!
   IF ( icount==0 ) RETURN
!
!     shift the nonzero elements if needed, to allow for created
!     diagonal elements.
!
   ko = Ib(Nrow+1) + icount
!
!     copy rows backward
!
   DO ii = Nrow , 1 , -1
!
!     go through  row ii
!
      k1 = Ib(ii)
      k2 = Ib(ii+1) - 1
      Ib(ii+1) = ko
      test = (Iw(ii)==0)
      DO k = k2 , k1 , -1
         j = Jb(k)
         IF ( test .AND. (j<ii) ) THEN
            test = .FALSE.
            ko = ko - 1
            B(ko) = Diag(ii)
            Jb(ko) = ii
            Iw(ii) = ko
         ENDIF
         ko = ko - 1
         B(ko) = A(k)
         Jb(ko) = j
      ENDDO
!     diagonal element has not been added yet.
      IF ( test ) THEN
         ko = ko - 1
         B(ko) = Diag(ii)
         Jb(ko) = ii
         Iw(ii) = ko
      ENDIF
   ENDDO
   Ib(1) = ko
!-----------------------------------------------------------------------
!------------end-of-apldiag---------------------------------------------
END SUBROUTINE apldia
