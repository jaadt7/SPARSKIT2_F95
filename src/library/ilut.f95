!*==ilut.f90 processed by SPAG 8.04RA 15:47 21 Feb 2024
!!SPAG Open source Personal, Educational or Academic User Clemson University  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
!-----end of lctcsr-----------------------------------------------------
 
 
!----------------------------------------------------------------------c
!                          S P A R S K I T                             c
!----------------------------------------------------------------------c
!                   ITERATIVE SOLVERS MODULE                           c
!----------------------------------------------------------------------c
! This Version Dated: August 13, 1996. Warning: meaning of some        c
! ============ arguments have changed w.r.t. earlier versions. Some    c
!              Calling sequences may also have changed                 c
!----------------------------------------------------------------------c
! Contents:                                                            c
!-------------------------preconditioners------------------------------c
!                                                                      c
! ILUT    : Incomplete LU factorization with dual truncation strategy  c
! ILUTP   : ILUT with column  pivoting                                 c
! ILUD    : ILU with single dropping + diagonal compensation (~MILUT)  c
! ILUDP   : ILUD with column pivoting                                  c
! ILUK    : level-k ILU                                                c
! ILU0    : simple ILU(0) preconditioning                              c
! MILU0   : MILU(0) preconditioning                                    c
!                                                                      c
!----------sample-accelerator-and-LU-solvers---------------------------c
!                                                                      c
! PGMRES  : preconditioned GMRES solver                                c
! LUSOL   : forward followed by backward triangular solve (Precond.)   c
! LUTSOL  : solving v = (LU)^{-T} u (used for preconditioning)         c
!                                                                      c
!-------------------------utility-routine------------------------------c
!                                                                      c
! QSPLIT  : quick split routine used by ilut to sort out the k largest c
!           elements in absolute value                                 c
!                                                                      c
!----------------------------------------------------------------------c
!                                                                      c
! Note: all preconditioners are preprocessors to pgmres.               c
! usage: call preconditioner then call pgmres                          c
!                                                                      c
!----------------------------------------------------------------------c
SUBROUTINE ilut(N,A,Ja,Ia,Lfil,Droptol,Alu,Jlu,Ju,Iwk,W,Jw,Ierr)
!-----------------------------------------------------------------------
   USE ISO_FORTRAN_ENV                 
   IMPLICIT NONE
!
! Dummy argument declarations rewritten by SPAG
!
   INTEGER , INTENT(IN) :: N
   REAL(REAL64) , INTENT(IN) , DIMENSION(*) :: A
   INTEGER , INTENT(IN) , DIMENSION(*) :: Ja
   INTEGER , INTENT(IN) , DIMENSION(N+1) :: Ia
   INTEGER , INTENT(IN) :: Lfil
   REAL(REAL64) , INTENT(IN) :: Droptol
   REAL(REAL64) , INTENT(INOUT) , DIMENSION(*) :: Alu
   INTEGER , INTENT(INOUT) , DIMENSION(*) :: Jlu
   INTEGER , INTENT(INOUT) , DIMENSION(N) :: Ju
   INTEGER , INTENT(IN) :: Iwk
   REAL(REAL64) , INTENT(INOUT) , DIMENSION(N+1) :: W
   INTEGER , INTENT(INOUT) , DIMENSION(2*N) :: Jw
   INTEGER , INTENT(OUT) :: Ierr
!
! Local variable declarations rewritten by SPAG
!
   REAL(REAL64) :: fact , s , t , tnorm
   INTEGER :: i , ii , j , j1 , j2 , jj , jpos , jrow , ju0 , k , len , lenl , lenu
   EXTERNAL qsplit
!
! End of declarations rewritten by SPAG
!
!----------------------------------------------------------------------*
!                      *** ILUT preconditioner ***                     *
!      incomplete LU factorization with dual truncation mechanism      *
!----------------------------------------------------------------------*
!     Author: Yousef Saad *May, 5, 1990, Latest revision, August 1996  *
!----------------------------------------------------------------------*
! PARAMETERS
!-----------
!
! on entry:
!==========
! n       = integer. The row dimension of the matrix A. The matrix
!
! a,ja,ia = matrix stored in Compressed Sparse Row format.
!
! lfil    = integer. The fill-in parameter. Each row of L and each row
!           of U will have a maximum of lfil elements (excluding the
!           diagonal element). lfil must be .ge. 0.
!           ** WARNING: THE MEANING OF LFIL HAS CHANGED WITH RESPECT TO
!           EARLIER VERSIONS.
!
! droptol = real*8. Sets the threshold for dropping small terms in the
!           factorization. See below for details on dropping strategy.
!
!
! iwk     = integer. The lengths of arrays alu and jlu. If the arrays
!           are not big enough to store the ILU factorizations, ilut
!           will stop with an error message.
!
! On return:
!===========
!
! alu,jlu = matrix stored in Modified Sparse Row (MSR) format containing
!           the L and U factors together. The diagonal (stored in
!           alu(1:n) ) is inverted. Each i-th row of the alu,jlu matrix
!           contains the i-th row of L (excluding the diagonal entry=1)
!           followed by the i-th row of U.
!
! ju      = integer array of length n containing the pointers to
!           the beginning of each row of U in the matrix alu,jlu.
!
! ierr    = integer. Error message with the following meaning.
!           ierr  = 0    --> successful return.
!           ierr .gt. 0  --> zero pivot encountered at step number ierr.
!           ierr  = -1   --> Error. input matrix may be wrong.
!                            (The elimination process has generated a
!                            row in L or U whose length is .gt.  n.)
!           ierr  = -2   --> The matrix L overflows the array al.
!           ierr  = -3   --> The matrix U overflows the array alu.
!           ierr  = -4   --> Illegal value for lfil.
!           ierr  = -5   --> zero row encountered.
!
! work arrays:
!=============
! jw      = integer work array of length 2*n.
! w       = real work array of length n+1.
!
!----------------------------------------------------------------------
! w, ju (1:n) store the working array [1:ii-1 = L-part, ii:n = u]
! jw(n+1:2n)  stores nonzero indicators
!
! Notes:
! ------
! The diagonal elements of the input matrix must be  nonzero (at least
! 'structurally').
!
!----------------------------------------------------------------------*
!---- Dual drop strategy works as follows.                             *
!                                                                      *
!     1) Theresholding in L and U as set by droptol. Any element whose *
!        magnitude is less than some tolerance (relative to the abs    *
!        value of diagonal element in u) is dropped.                   *
!                                                                      *
!     2) Keeping only the largest lfil elements in the i-th row of L   *
!        and the largest lfil elements in the i-th row of U (excluding *
!        diagonal elements).                                           *
!                                                                      *
! Flexibility: one  can use  droptol=0  to get  a strategy  based on   *
! keeping  the largest  elements in  each row  of L  and U.   Taking   *
! droptol .ne.  0 but lfil=n will give  the usual threshold strategy   *
! (however, fill-in is then mpredictible).                             *
!----------------------------------------------------------------------*
!     locals
   INTEGER :: spag_nextblock_1
   spag_nextblock_1 = 1
   SPAG_DispatchLoop_1: DO
      SELECT CASE (spag_nextblock_1)
      CASE (1)
         IF ( Lfil<0 ) THEN
!
!     illegal lfil entered.
!
            Ierr = -4
            RETURN
         ELSE
!-----------------------------------------------------------------------
!     initialize ju0 (points to next element to be added to alu,jlu)
!     and pointer array.
!-----------------------------------------------------------------------
            ju0 = N + 2
            Jlu(1) = ju0
!
!     initialize nonzero indicator array.
!
            DO j = 1 , N
               Jw(N+j) = 0
            ENDDO
!-----------------------------------------------------------------------
!     beginning of main loop.
!-----------------------------------------------------------------------
            DO ii = 1 , N
               j1 = Ia(ii)
               j2 = Ia(ii+1) - 1
               tnorm = 0.0D0
               DO k = j1 , j2
                  tnorm = tnorm + abs(A(k))
               ENDDO
               IF ( tnorm==0.0 ) THEN
                  spag_nextblock_1 = 5
                  CYCLE SPAG_DispatchLoop_1
               ENDIF
               tnorm = tnorm/real(j2-j1+1)
!
!     unpack L-part and U-part of row of A in arrays w
!
               lenu = 1
               lenl = 0
               Jw(ii) = ii
               W(ii) = 0.0
               Jw(N+ii) = ii
!
               DO j = j1 , j2
                  k = Ja(j)
                  t = A(j)
                  IF ( k<ii ) THEN
                     lenl = lenl + 1
                     Jw(lenl) = k
                     W(lenl) = t
                     Jw(N+k) = lenl
                  ELSEIF ( k==ii ) THEN
                     W(ii) = t
                  ELSE
                     lenu = lenu + 1
                     jpos = ii + lenu - 1
                     Jw(jpos) = k
                     W(jpos) = t
                     Jw(N+k) = jpos
                  ENDIF
               ENDDO
               jj = 0
               len = 0
               SPAG_Loop_3_1: DO
!
!     eliminate previous rows
!
                  jj = jj + 1
                  IF ( jj>lenl ) THEN
!
!     reset double-pointer to zero (U-part)
!
                     DO k = 1 , lenu
                        Jw(N+Jw(ii+k-1)) = 0
                     ENDDO
!
!     update L-matrix
!
                     lenl = len
                     len = min0(lenl,Lfil)
!
!     sort by quick-split
!
                     CALL qsplit(W,Jw,lenl,len)
!
!     store L-part
!
                     DO k = 1 , len
                        IF ( ju0>Iwk ) THEN
                           spag_nextblock_1 = 3
                           CYCLE SPAG_DispatchLoop_1
                        ENDIF
                        Alu(ju0) = W(k)
                        Jlu(ju0) = Jw(k)
                        ju0 = ju0 + 1
                     ENDDO
!
!     save pointer to beginning of row ii of U
!
                     Ju(ii) = ju0
!
!     update U-matrix -- first apply dropping strategy
!
                     len = 0
                     DO k = 1 , lenu - 1
                        IF ( abs(W(ii+k))>Droptol*tnorm ) THEN
                           len = len + 1
                           W(ii+len) = W(ii+k)
                           Jw(ii+len) = Jw(ii+k)
                        ENDIF
                     ENDDO
                     lenu = len + 1
                     len = min0(lenu,Lfil)
!
                     CALL qsplit(W(ii+1),Jw(ii+1),lenu-1,len)
!
!     copy
!
                     t = abs(W(ii))
                     IF ( len+ju0>Iwk ) THEN
                        spag_nextblock_1 = 4
                        CYCLE SPAG_DispatchLoop_1
                     ENDIF
                     DO k = ii + 1 , ii + len - 1
                        Jlu(ju0) = Jw(k)
                        Alu(ju0) = W(k)
                        t = t + abs(W(k))
                        ju0 = ju0 + 1
                     ENDDO
!
!     store inverse of diagonal element of u
!
                     IF ( W(ii)==0.0 ) W(ii) = (0.0001+Droptol)*tnorm
!
                     Alu(ii) = 1.0D0/W(ii)
!
!     update pointer to beginning of next row of U.
!
                     Jlu(ii+1) = ju0
                     EXIT SPAG_Loop_3_1
                  ELSE
!-----------------------------------------------------------------------
!     in order to do the elimination in the correct order we must select
!     the smallest column index among jw(k), k=jj+1, ..., lenl.
!-----------------------------------------------------------------------
                     jrow = Jw(jj)
                     k = jj
!
!     determine smallest column index
!
                     DO j = jj + 1 , lenl
                        IF ( Jw(j)<jrow ) THEN
                           jrow = Jw(j)
                           k = j
                        ENDIF
                     ENDDO
!
                     IF ( k/=jj ) THEN
!     exchange in jw
                        j = Jw(jj)
                        Jw(jj) = Jw(k)
                        Jw(k) = j
!     exchange in jr
                        Jw(N+jrow) = jj
                        Jw(N+j) = k
!     exchange in w
                        s = W(jj)
                        W(jj) = W(k)
                        W(k) = s
                     ENDIF
!
!     zero out element in row by setting jw(n+jrow) to zero.
!
                     Jw(N+jrow) = 0
!
!     get the multiplier for row to be eliminated (jrow).
!
                     fact = W(jj)*Alu(jrow)
                     IF ( abs(fact)>Droptol ) THEN
!
!     combine current row and row jrow
!
                        DO k = Ju(jrow) , Jlu(jrow+1) - 1
                           s = fact*Alu(k)
                           j = Jlu(k)
                           jpos = Jw(N+j)
                           IF ( j>=ii ) THEN
!
!     dealing with upper part.
!
                              IF ( jpos==0 ) THEN
!
!     this is a fill-in element
!
                                 lenu = lenu + 1
                                 IF ( lenu>N ) THEN
                                    spag_nextblock_1 = 2
                                    CYCLE SPAG_DispatchLoop_1
                                 ENDIF
                                 i = ii + lenu - 1
                                 Jw(i) = j
                                 Jw(N+j) = i
                                 W(i) = -s
                              ELSE
!
!     this is not a fill-in element
!
                                 W(jpos) = W(jpos) - s
 
                              ENDIF
!
!     dealing  with lower part.
!
                           ELSEIF ( jpos==0 ) THEN
!
!     this is a fill-in element
!
                              lenl = lenl + 1
                              IF ( lenl>N ) THEN
                                 spag_nextblock_1 = 2
                                 CYCLE SPAG_DispatchLoop_1
                              ENDIF
                              Jw(lenl) = j
                              Jw(N+j) = lenl
                              W(lenl) = -s
                           ELSE
!
!     this is not a fill-in element
!
                              W(jpos) = W(jpos) - s
                           ENDIF
                        ENDDO
!
!     store this pivot element -- (from left to right -- no danger of
!     overlap with the working elements in L (pivots).
!
                        len = len + 1
                        W(len) = fact
                        Jw(len) = jrow
                     ENDIF
                  ENDIF
               ENDDO SPAG_Loop_3_1
!-----------------------------------------------------------------------
!     end main loop
!-----------------------------------------------------------------------
            ENDDO
            Ierr = 0
            RETURN
         ENDIF
      CASE (2)
!
!     incomprehensible error. Matrix must be wrong.
!
         Ierr = -1
         RETURN
      CASE (3)
!
!     insufficient storage in L.
!
         Ierr = -2
         RETURN
      CASE (4)
!
!     insufficient storage in U.
!
         Ierr = -3
         RETURN
      CASE (5)
!
!     zero row encountered
!
         Ierr = -5
         EXIT SPAG_DispatchLoop_1
      END SELECT
   ENDDO SPAG_DispatchLoop_1
!----------------end-of-ilut--------------------------------------------
!-----------------------------------------------------------------------
END SUBROUTINE ilut
!*==ilutp.f90 processed by SPAG 8.04RA 15:47 21 Feb 2024
!!SPAG Open source Personal, Educational or Academic User Clemson University  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
!----------------------------------------------------------------------
SUBROUTINE ilutp(N,A,Ja,Ia,Lfil,Droptol,Permtol,Mbloc,Alu,Jlu,Ju,Iwk,W,Jw,Iperm,Ierr)
   USE ISO_FORTRAN_ENV                 
   IMPLICIT NONE
!
! Dummy argument declarations rewritten by SPAG
!
   INTEGER , INTENT(IN) :: N
   REAL(REAL64) , INTENT(IN) , DIMENSION(*) :: A
   INTEGER , INTENT(INOUT) , DIMENSION(*) :: Ja
   INTEGER , INTENT(IN) , DIMENSION(N+1) :: Ia
   INTEGER , INTENT(IN) :: Lfil
   REAL(REAL64) , INTENT(IN) :: Droptol
   REAL(REAL64) , INTENT(IN) :: Permtol
   INTEGER , INTENT(IN) :: Mbloc
   REAL(REAL64) , INTENT(INOUT) , DIMENSION(*) :: Alu
   INTEGER , INTENT(INOUT) , DIMENSION(*) :: Jlu
   INTEGER , INTENT(INOUT) , DIMENSION(N) :: Ju
   INTEGER , INTENT(IN) :: Iwk
   REAL(REAL64) , INTENT(INOUT) , DIMENSION(N+1) :: W
   INTEGER , INTENT(INOUT) , DIMENSION(2*N) :: Jw
   INTEGER , INTENT(INOUT) , DIMENSION(2*N) :: Iperm
   INTEGER , INTENT(OUT) :: Ierr
!
! Local variable declarations rewritten by SPAG
!
   REAL(REAL64) :: fact , s , t , tmp , tnorm , xmax , xmax0
   INTEGER :: i , icut , ii , imax , j , j1 , j2 , jj , jpos , jrow , ju0 , k , len , lenl , lenu
   EXTERNAL qsplit
!
! End of declarations rewritten by SPAG
!
!-----------------------------------------------------------------------
!     implicit none
!----------------------------------------------------------------------*
!       *** ILUTP preconditioner -- ILUT with pivoting  ***            *
!      incomplete LU factorization with dual truncation mechanism      *
!----------------------------------------------------------------------*
! author Yousef Saad *Sep 8, 1993 -- Latest revision, August 1996.     *
!----------------------------------------------------------------------*
! on entry:
!==========
! n       = integer. The dimension of the matrix A.
!
! a,ja,ia = matrix stored in Compressed Sparse Row format.
!           ON RETURN THE COLUMNS OF A ARE PERMUTED. SEE BELOW FOR
!           DETAILS.
!
! lfil    = integer. The fill-in parameter. Each row of L and each row
!           of U will have a maximum of lfil elements (excluding the
!           diagonal element). lfil must be .ge. 0.
!           ** WARNING: THE MEANING OF LFIL HAS CHANGED WITH RESPECT TO
!           EARLIER VERSIONS.
!
! droptol = real*8. Sets the threshold for dropping small terms in the
!           factorization. See below for details on dropping strategy.
!
! lfil    = integer. The fill-in parameter. Each row of L and
!           each row of U will have a maximum of lfil elements.
!           WARNING: THE MEANING OF LFIL HAS CHANGED WITH RESPECT TO
!           EARLIER VERSIONS.
!           lfil must be .ge. 0.
!
! permtol = tolerance ratio used to  determne whether or not to permute
!           two columns.  At step i columns i and j are permuted when
!
!                     abs(a(i,j))*permtol .gt. abs(a(i,i))
!
!           [0 --> never permute; good values 0.1 to 0.01]
!
! mbloc   = if desired, permuting can be done only within the diagonal
!           blocks of size mbloc. Useful for PDE problems with several
!           degrees of freedom.. If feature not wanted take mbloc=n.
!
!
! iwk     = integer. The lengths of arrays alu and jlu. If the arrays
!           are not big enough to store the ILU factorizations, ilut
!           will stop with an error message.
!
! On return:
!===========
!
! alu,jlu = matrix stored in Modified Sparse Row (MSR) format containing
!           the L and U factors together. The diagonal (stored in
!           alu(1:n) ) is inverted. Each i-th row of the alu,jlu matrix
!           contains the i-th row of L (excluding the diagonal entry=1)
!           followed by the i-th row of U.
!
! ju      = integer array of length n containing the pointers to
!           the beginning of each row of U in the matrix alu,jlu.
!
! iperm   = contains the permutation arrays.
!           iperm(1:n) = old numbers of unknowns
!           iperm(n+1:2*n) = reverse permutation = new unknowns.
!
! ierr    = integer. Error message with the following meaning.
!           ierr  = 0    --> successful return.
!           ierr .gt. 0  --> zero pivot encountered at step number ierr.
!           ierr  = -1   --> Error. input matrix may be wrong.
!                            (The elimination process has generated a
!                            row in L or U whose length is .gt.  n.)
!           ierr  = -2   --> The matrix L overflows the array al.
!           ierr  = -3   --> The matrix U overflows the array alu.
!           ierr  = -4   --> Illegal value for lfil.
!           ierr  = -5   --> zero row encountered.
!
! work arrays:
!=============
! jw      = integer work array of length 2*n.
! w       = real work array of length n
!
! IMPORTANR NOTE:
! --------------
! TO AVOID PERMUTING THE SOLUTION VECTORS ARRAYS FOR EACH LU-SOLVE,
! THE MATRIX A IS PERMUTED ON RETURN. [all column indices are
! changed]. SIMILARLY FOR THE U MATRIX.
! To permute the matrix back to its original state use the loop:
!
!      do k=ia(1), ia(n+1)-1
!         ja(k) = iperm(ja(k))
!      enddo
!
!-----------------------------------------------------------------------
!     local variables
!
   INTEGER :: spag_nextblock_1
   spag_nextblock_1 = 1
   SPAG_DispatchLoop_1: DO
      SELECT CASE (spag_nextblock_1)
      CASE (1)
!
         IF ( Lfil<0 ) THEN
!
!     illegal lfil entered.
!
            Ierr = -4
            RETURN
         ELSE
!-----------------------------------------------------------------------
!     initialize ju0 (points to next element to be added to alu,jlu)
!     and pointer array.
!-----------------------------------------------------------------------
            ju0 = N + 2
            Jlu(1) = ju0
!
!  integer double pointer array.
!
            DO j = 1 , N
               Jw(N+j) = 0
               Iperm(j) = j
               Iperm(N+j) = j
            ENDDO
!-----------------------------------------------------------------------
!     beginning of main loop.
!-----------------------------------------------------------------------
            DO ii = 1 , N
               j1 = Ia(ii)
               j2 = Ia(ii+1) - 1
               tnorm = 0.0D0
               DO k = j1 , j2
                  tnorm = tnorm + abs(A(k))
               ENDDO
               IF ( tnorm==0.0 ) THEN
                  spag_nextblock_1 = 5
                  CYCLE SPAG_DispatchLoop_1
               ENDIF
               tnorm = tnorm/(j2-j1+1)
!
!     unpack L-part and U-part of row of A in arrays  w  --
!
               lenu = 1
               lenl = 0
               Jw(ii) = ii
               W(ii) = 0.0
               Jw(N+ii) = ii
!
               DO j = j1 , j2
                  k = Iperm(N+Ja(j))
                  t = A(j)
                  IF ( k<ii ) THEN
                     lenl = lenl + 1
                     Jw(lenl) = k
                     W(lenl) = t
                     Jw(N+k) = lenl
                  ELSEIF ( k==ii ) THEN
                     W(ii) = t
                  ELSE
                     lenu = lenu + 1
                     jpos = ii + lenu - 1
                     Jw(jpos) = k
                     W(jpos) = t
                     Jw(N+k) = jpos
                  ENDIF
               ENDDO
               jj = 0
               len = 0
               SPAG_Loop_3_1: DO
!
!     eliminate previous rows
!
                  jj = jj + 1
                  IF ( jj>lenl ) THEN
!
!     reset double-pointer to zero (U-part)
!
                     DO k = 1 , lenu
                        Jw(N+Jw(ii+k-1)) = 0
                     ENDDO
!
!     update L-matrix
!
                     lenl = len
                     len = min0(lenl,Lfil)
!
!     sort by quick-split
!
                     CALL qsplit(W,Jw,lenl,len)
!
!     store L-part -- in original coordinates ..
!
                     DO k = 1 , len
                        IF ( ju0>Iwk ) THEN
                           spag_nextblock_1 = 3
                           CYCLE SPAG_DispatchLoop_1
                        ENDIF
                        Alu(ju0) = W(k)
                        Jlu(ju0) = Iperm(Jw(k))
                        ju0 = ju0 + 1
                     ENDDO
!
!     save pointer to beginning of row ii of U
!
                     Ju(ii) = ju0
!
!     update U-matrix -- first apply dropping strategy
!
                     len = 0
                     DO k = 1 , lenu - 1
                        IF ( abs(W(ii+k))>Droptol*tnorm ) THEN
                           len = len + 1
                           W(ii+len) = W(ii+k)
                           Jw(ii+len) = Jw(ii+k)
                        ENDIF
                     ENDDO
                     lenu = len + 1
                     len = min0(lenu,Lfil)
                     CALL qsplit(W(ii+1),Jw(ii+1),lenu-1,len)
!
!     determine next pivot --
!
                     imax = ii
                     xmax = abs(W(imax))
                     xmax0 = xmax
                     icut = ii - 1 + Mbloc - mod(ii-1,Mbloc)
                     DO k = ii + 1 , ii + len - 1
                        t = abs(W(k))
                        IF ( t>xmax .AND. t*Permtol>xmax0 .AND. Jw(k)<=icut ) THEN
                           imax = k
                           xmax = t
                        ENDIF
                     ENDDO
!
!     exchange w's
!
                     tmp = W(ii)
                     W(ii) = W(imax)
                     W(imax) = tmp
!
!     update iperm and reverse iperm
!
                     j = Jw(imax)
                     i = Iperm(ii)
                     Iperm(ii) = Iperm(j)
                     Iperm(j) = i
!
!     reverse iperm
!
                     Iperm(N+Iperm(ii)) = ii
                     Iperm(N+Iperm(j)) = j
!-----------------------------------------------------------------------
!
                     IF ( len+ju0>Iwk ) THEN
                        spag_nextblock_1 = 4
                        CYCLE SPAG_DispatchLoop_1
                     ENDIF
!
!     copy U-part in original coordinates
!
                     DO k = ii + 1 , ii + len - 1
                        Jlu(ju0) = Iperm(Jw(k))
                        Alu(ju0) = W(k)
                        ju0 = ju0 + 1
                     ENDDO
!
!     store inverse of diagonal element of u
!
                     IF ( W(ii)==0.0 ) W(ii) = (1.0D-4+Droptol)*tnorm
                     Alu(ii) = 1.0D0/W(ii)
!
!     update pointer to beginning of next row of U.
!
                     Jlu(ii+1) = ju0
                     EXIT SPAG_Loop_3_1
                  ELSE
!-----------------------------------------------------------------------
!     in order to do the elimination in the correct order we must select
!     the smallest column index among jw(k), k=jj+1, ..., lenl.
!-----------------------------------------------------------------------
                     jrow = Jw(jj)
                     k = jj
!
!     determine smallest column index
!
                     DO j = jj + 1 , lenl
                        IF ( Jw(j)<jrow ) THEN
                           jrow = Jw(j)
                           k = j
                        ENDIF
                     ENDDO
!
                     IF ( k/=jj ) THEN
!     exchange in jw
                        j = Jw(jj)
                        Jw(jj) = Jw(k)
                        Jw(k) = j
!     exchange in jr
                        Jw(N+jrow) = jj
                        Jw(N+j) = k
!     exchange in w
                        s = W(jj)
                        W(jj) = W(k)
                        W(k) = s
                     ENDIF
!
!     zero out element in row by resetting jw(n+jrow) to zero.
!
                     Jw(N+jrow) = 0
!
!     get the multiplier for row to be eliminated: jrow
!
                     fact = W(jj)*Alu(jrow)
!
!     drop term if small
!
                     IF ( abs(fact)>Droptol ) THEN
!
!     combine current row and row jrow
!
                        DO k = Ju(jrow) , Jlu(jrow+1) - 1
                           s = fact*Alu(k)
!     new column number
                           j = Iperm(N+Jlu(k))
                           jpos = Jw(N+j)
                           IF ( j>=ii ) THEN
!
!     dealing with upper part.
!
                              IF ( jpos==0 ) THEN
!
!     this is a fill-in element
!
                                 lenu = lenu + 1
                                 i = ii + lenu - 1
                                 IF ( lenu>N ) THEN
                                    spag_nextblock_1 = 2
                                    CYCLE SPAG_DispatchLoop_1
                                 ENDIF
                                 Jw(i) = j
                                 Jw(N+j) = i
                                 W(i) = -s
                              ELSE
!     no fill-in element --
                                 W(jpos) = W(jpos) - s
                              ENDIF
!
!     dealing with lower part.
!
                           ELSEIF ( jpos==0 ) THEN
!
!     this is a fill-in element
!
                              lenl = lenl + 1
                              IF ( lenl>N ) THEN
                                 spag_nextblock_1 = 2
                                 CYCLE SPAG_DispatchLoop_1
                              ENDIF
                              Jw(lenl) = j
                              Jw(N+j) = lenl
                              W(lenl) = -s
                           ELSE
!
!     this is not a fill-in element
!
                              W(jpos) = W(jpos) - s
                           ENDIF
                        ENDDO
!
!     store this pivot element -- (from left to right -- no danger of
!     overlap with the working elements in L (pivots).
!
                        len = len + 1
                        W(len) = fact
                        Jw(len) = jrow
                     ENDIF
                  ENDIF
               ENDDO SPAG_Loop_3_1
!-----------------------------------------------------------------------
!     end main loop
!-----------------------------------------------------------------------
            ENDDO
!
!     permute all column indices of LU ...
!
            DO k = Jlu(1) , Jlu(N+1) - 1
               Jlu(k) = Iperm(N+Jlu(k))
            ENDDO
!
!     ...and of A
!
            DO k = Ia(1) , Ia(N+1) - 1
               Ja(k) = Iperm(N+Ja(k))
            ENDDO
!
            Ierr = 0
            RETURN
         ENDIF
      CASE (2)
!
!     incomprehensible error. Matrix must be wrong.
!
         Ierr = -1
         RETURN
      CASE (3)
!
!     insufficient storage in L.
!
         Ierr = -2
         RETURN
      CASE (4)
!
!     insufficient storage in U.
!
         Ierr = -3
         RETURN
      CASE (5)
!
!     zero row encountered
!
         Ierr = -5
         EXIT SPAG_DispatchLoop_1
      END SELECT
   ENDDO SPAG_DispatchLoop_1
!----------------end-of-ilutp-------------------------------------------
!-----------------------------------------------------------------------
END SUBROUTINE ilutp
!*==ilud.f90 processed by SPAG 8.04RA 15:47 21 Feb 2024
!!SPAG Open source Personal, Educational or Academic User Clemson University  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
!-----------------------------------------------------------------------
SUBROUTINE ilud(N,A,Ja,Ia,Alph,Tol,Alu,Jlu,Ju,Iwk,W,Jw,Ierr)
!-----------------------------------------------------------------------
   USE ISO_FORTRAN_ENV                 
   IMPLICIT NONE
!
! Dummy argument declarations rewritten by SPAG
!
   INTEGER , INTENT(IN) :: N
   REAL(REAL64) , INTENT(IN) , DIMENSION(*) :: A
   INTEGER , INTENT(IN) , DIMENSION(*) :: Ja
   INTEGER , INTENT(IN) , DIMENSION(N+1) :: Ia
   REAL(REAL64) , INTENT(IN) :: Alph
   REAL(REAL64) , INTENT(IN) :: Tol
   REAL(REAL64) , INTENT(INOUT) , DIMENSION(*) :: Alu
   INTEGER , INTENT(INOUT) , DIMENSION(*) :: Jlu
   INTEGER , INTENT(INOUT) , DIMENSION(N) :: Ju
   INTEGER , INTENT(IN) :: Iwk
   REAL(REAL64) , INTENT(INOUT) , DIMENSION(2*N) :: W
   INTEGER , INTENT(INOUT) , DIMENSION(2*N) :: Jw
   INTEGER , INTENT(OUT) :: Ierr
!
! Local variable declarations rewritten by SPAG
!
   REAL(REAL64) :: dropsum , fact , s , t , tnorm
   INTEGER :: i , ii , j , j1 , j2 , jj , jpos , jrow , ju0 , k , len , lenl , lenu
!
! End of declarations rewritten by SPAG
!
!----------------------------------------------------------------------*
!                     *** ILUD preconditioner ***                      *
!    incomplete LU factorization with standard droppoing strategy      *
!----------------------------------------------------------------------*
! Author: Yousef Saad * Aug. 1995 --                                   *
!----------------------------------------------------------------------*
! This routine computes the ILU factorization with standard threshold  *
! dropping: at i-th step of elimination, an element a(i,j) in row i is *
! dropped  if it satisfies the criterion:                              *
!                                                                      *
!  abs(a(i,j)) < tol * [average magnitude of elements in row i of A]   *
!                                                                      *
! There is no control on memory size required for the factors as is    *
! done in ILUT. This routines computes also various diagonal compensa- *
! tion ILU's such MILU. These are defined through the parameter alph   *
!----------------------------------------------------------------------*
! on entry:
!==========
! n       = integer. The row dimension of the matrix A. The matrix
!
! a,ja,ia = matrix stored in Compressed Sparse Row format
!
! alph    = diagonal compensation parameter -- the term:
!
!           alph*(sum of all dropped out elements in a given row)
!
!           is added to the diagonal element of U of the factorization
!           Thus: alph = 0 ---> ~ ILU with threshold,
!                 alph = 1 ---> ~ MILU with threshold.
!
! tol     = Threshold parameter for dropping small terms in the
!           factorization. During the elimination, a term a(i,j) is
!           dropped whenever abs(a(i,j)) .lt. tol * [weighted norm of
!           row i]. Here weighted norm = 1-norm / number of nnz
!           elements in the row.
!
! iwk     = The length of arrays alu and jlu -- this routine will stop
!           if storage for the factors L and U is not sufficient
!
! On return:
!===========
!
! alu,jlu = matrix stored in Modified Sparse Row (MSR) format containing
!           the L and U factors together. The diagonal (stored in
!           alu(1:n) ) is inverted. Each i-th row of the alu,jlu matrix
!           contains the i-th row of L (excluding the diagonal entry=1)
!           followed by the i-th row of U.
!
! ju      = integer array of length n containing the pointers to
!           the beginning of each row of U in the matrix alu,jlu.
!
! ierr    = integer. Error message with the following meaning.
!           ierr  = 0    --> successful return.
!           ierr .gt. 0  --> zero pivot encountered at step number ierr.
!           ierr  = -1   --> Error. input matrix may be wrong.
!                            (The elimination process has generated a
!                            row in L or U whose length is .gt.  n.)
!           ierr  = -2   --> Insufficient storage for the LU factors --
!                            arrays alu/ jalu are  overflowed.
!           ierr  = -3   --> Zero row encountered.
!
! Work Arrays:
!=============
! jw      = integer work array of length 2*n.
! w       = real work array of length n
!
!----------------------------------------------------------------------
!
! w, ju (1:n) store the working array [1:ii-1 = L-part, ii:n = u]
! jw(n+1:2n)  stores the nonzero indicator.
!
! Notes:
! ------
! All diagonal elements of the input matrix must be  nonzero.
!
!-----------------------------------------------------------------------
!     locals
   INTEGER :: spag_nextblock_1
   spag_nextblock_1 = 1
   SPAG_DispatchLoop_1: DO
      SELECT CASE (spag_nextblock_1)
      CASE (1)
!-----------------------------------------------------------------------
!     initialize ju0 (points to next element to be added to alu,jlu)
!     and pointer array.
!-----------------------------------------------------------------------
         ju0 = N + 2
         Jlu(1) = ju0
!
!     initialize nonzero indicator array.
!
         DO j = 1 , N
            Jw(N+j) = 0
         ENDDO
!-----------------------------------------------------------------------
!     beginning of main loop.
!-----------------------------------------------------------------------
         DO ii = 1 , N
            j1 = Ia(ii)
            j2 = Ia(ii+1) - 1
            dropsum = 0.0D0
            tnorm = 0.0D0
            DO k = j1 , j2
               tnorm = tnorm + abs(A(k))
            ENDDO
            IF ( tnorm==0.0 ) THEN
               spag_nextblock_1 = 4
               CYCLE SPAG_DispatchLoop_1
            ENDIF
            tnorm = tnorm/real(j2-j1+1)
!
!     unpack L-part and U-part of row of A in arrays w
!
            lenu = 1
            lenl = 0
            Jw(ii) = ii
            W(ii) = 0.0
            Jw(N+ii) = ii
!
            DO j = j1 , j2
               k = Ja(j)
               t = A(j)
               IF ( k<ii ) THEN
                  lenl = lenl + 1
                  Jw(lenl) = k
                  W(lenl) = t
                  Jw(N+k) = lenl
               ELSEIF ( k==ii ) THEN
                  W(ii) = t
               ELSE
                  lenu = lenu + 1
                  jpos = ii + lenu - 1
                  Jw(jpos) = k
                  W(jpos) = t
                  Jw(N+k) = jpos
               ENDIF
            ENDDO
            jj = 0
            len = 0
            SPAG_Loop_3_1: DO
!
!     eliminate previous rows
!
               jj = jj + 1
               IF ( jj>lenl ) THEN
!
!     reset double-pointer to zero (For U-part only)
!
                  DO k = 1 , lenu
                     Jw(N+Jw(ii+k-1)) = 0
                  ENDDO
!
!     update l-matrix
!
                  DO k = 1 , len
                     IF ( ju0>Iwk ) THEN
                        spag_nextblock_1 = 3
                        CYCLE SPAG_DispatchLoop_1
                     ENDIF
                     Alu(ju0) = W(k)
                     Jlu(ju0) = Jw(k)
                     ju0 = ju0 + 1
                  ENDDO
!
!     save pointer to beginning of row ii of U
!
                  Ju(ii) = ju0
!
!     go through elements in U-part of w to determine elements to keep
!
                  len = 0
                  DO k = 1 , lenu - 1
!            if (abs(w(ii+k)) .gt. tnorm*tol) then
                     IF ( abs(W(ii+k))>abs(W(ii))*Tol ) THEN
                        len = len + 1
                        W(ii+len) = W(ii+k)
                        Jw(ii+len) = Jw(ii+k)
                     ELSE
                        dropsum = dropsum + W(ii+k)
                     ENDIF
                  ENDDO
!
!     now update u-matrix
!
                  IF ( ju0+len-1>Iwk ) THEN
                     spag_nextblock_1 = 3
                     CYCLE SPAG_DispatchLoop_1
                  ENDIF
                  DO k = ii + 1 , ii + len
                     Jlu(ju0) = Jw(k)
                     Alu(ju0) = W(k)
                     ju0 = ju0 + 1
                  ENDDO
!
!     define diagonal element
!
                  W(ii) = W(ii) + Alph*dropsum
!
!     store inverse of diagonal element of u
!
                  IF ( W(ii)==0.0 ) W(ii) = (0.0001+Tol)*tnorm
!
                  Alu(ii) = 1.0D0/W(ii)
!
!     update pointer to beginning of next row of U.
!
                  Jlu(ii+1) = ju0
                  EXIT SPAG_Loop_3_1
               ELSE
!-----------------------------------------------------------------------
!     in order to do the elimination in the correct order we must select
!     the smallest column index among jw(k), k=jj+1, ..., lenl.
!-----------------------------------------------------------------------
                  jrow = Jw(jj)
                  k = jj
!
!     determine smallest column index
!
                  DO j = jj + 1 , lenl
                     IF ( Jw(j)<jrow ) THEN
                        jrow = Jw(j)
                        k = j
                     ENDIF
                  ENDDO
!
                  IF ( k/=jj ) THEN
!     exchange in jw
                     j = Jw(jj)
                     Jw(jj) = Jw(k)
                     Jw(k) = j
!     exchange in jr
                     Jw(N+jrow) = jj
                     Jw(N+j) = k
!     exchange in w
                     s = W(jj)
                     W(jj) = W(k)
                     W(k) = s
                  ENDIF
!
!     zero out element in row by setting resetting jw(n+jrow) to zero.
!
                  Jw(N+jrow) = 0
!
!     drop term if small
!
!         if (abs(w(jj)) .le. tol*tnorm) then
!            dropsum = dropsum + w(jj)
!            goto 150
!         endif
!
!     get the multiplier for row to be eliminated (jrow).
!
                  fact = W(jj)*Alu(jrow)
!
!     drop term if small
!
                  IF ( abs(fact)<=Tol ) THEN
                     dropsum = dropsum + W(jj)
                     CYCLE
                  ENDIF
!
!     combine current row and row jrow
!
                  DO k = Ju(jrow) , Jlu(jrow+1) - 1
                     s = fact*Alu(k)
                     j = Jlu(k)
                     jpos = Jw(N+j)
                     IF ( j>=ii ) THEN
!
!     dealing with upper part.
!
                        IF ( jpos==0 ) THEN
!
!     this is a fill-in element
!
                           lenu = lenu + 1
                           IF ( lenu>N ) THEN
                              spag_nextblock_1 = 2
                              CYCLE SPAG_DispatchLoop_1
                           ENDIF
                           i = ii + lenu - 1
                           Jw(i) = j
                           Jw(N+j) = i
                           W(i) = -s
                        ELSE
!
!     this is not a fill-in element
!
                           W(jpos) = W(jpos) - s
                        ENDIF
!
!     dealing with lower part.
!
                     ELSEIF ( jpos==0 ) THEN
!
!     this is a fill-in element
!
                        lenl = lenl + 1
                        IF ( lenl>N ) THEN
                           spag_nextblock_1 = 2
                           CYCLE SPAG_DispatchLoop_1
                        ENDIF
                        Jw(lenl) = j
                        Jw(N+j) = lenl
                        W(lenl) = -s
                     ELSE
!
!     this is not a fill-in element
!
                        W(jpos) = W(jpos) - s
                     ENDIF
                  ENDDO
                  len = len + 1
                  W(len) = fact
                  Jw(len) = jrow
               ENDIF
            ENDDO SPAG_Loop_3_1
!-----------------------------------------------------------------------
!     end main loop
!-----------------------------------------------------------------------
         ENDDO
         Ierr = 0
         RETURN
      CASE (2)
!
!     incomprehensible error. Matrix must be wrong.
!
         Ierr = -1
         RETURN
      CASE (3)
!
!     insufficient storage in alu/ jlu arrays for  L / U factors
!
         Ierr = -2
         RETURN
      CASE (4)
!
!     zero row encountered
!
         Ierr = -3
         EXIT SPAG_DispatchLoop_1
      END SELECT
   ENDDO SPAG_DispatchLoop_1
!----------------end-of-ilud  ------------------------------------------
!-----------------------------------------------------------------------
END SUBROUTINE ilud
!*==iludp.f90 processed by SPAG 8.04RA 15:47 21 Feb 2024
!!SPAG Open source Personal, Educational or Academic User Clemson University  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
!----------------------------------------------------------------------
SUBROUTINE iludp(N,A,Ja,Ia,Alph,Droptol,Permtol,Mbloc,Alu,Jlu,Ju,Iwk,W,Jw,Iperm,Ierr)
!-----------------------------------------------------------------------
   USE ISO_FORTRAN_ENV                 
   IMPLICIT NONE
!
! Dummy argument declarations rewritten by SPAG
!
   INTEGER , INTENT(IN) :: N
   REAL(REAL64) , INTENT(IN) , DIMENSION(*) :: A
   INTEGER , INTENT(INOUT) , DIMENSION(*) :: Ja
   INTEGER , INTENT(IN) , DIMENSION(N+1) :: Ia
   REAL(REAL64) , INTENT(IN) :: Alph
   REAL(REAL64) , INTENT(IN) :: Droptol
   REAL(REAL64) , INTENT(IN) :: Permtol
   INTEGER , INTENT(IN) :: Mbloc
   REAL(REAL64) , INTENT(INOUT) , DIMENSION(*) :: Alu
   INTEGER , INTENT(INOUT) , DIMENSION(*) :: Jlu
   INTEGER , INTENT(INOUT) , DIMENSION(N) :: Ju
   INTEGER , INTENT(IN) :: Iwk
   REAL(REAL64) , INTENT(INOUT) , DIMENSION(2*N) :: W
   INTEGER , INTENT(INOUT) , DIMENSION(2*N) :: Jw
   INTEGER , INTENT(INOUT) , DIMENSION(2*N) :: Iperm
   INTEGER , INTENT(OUT) :: Ierr
!
! Local variable declarations rewritten by SPAG
!
   REAL(REAL64) :: dropsum , fact , s , t , tmp , tnorm , xmax , xmax0
   INTEGER :: i , icut , ii , imax , j , j1 , j2 , jj , jpos , jrow , ju0 , k , len , lenl , lenu
!
! End of declarations rewritten by SPAG
!
!----------------------------------------------------------------------*
!                     *** ILUDP preconditioner ***                     *
!    incomplete LU factorization with standard droppoing strategy      *
!    and column pivoting                                               *
!----------------------------------------------------------------------*
! author Yousef Saad -- Aug 1995.                                      *
!----------------------------------------------------------------------*
! on entry:
!==========
! n       = integer. The dimension of the matrix A.
!
! a,ja,ia = matrix stored in Compressed Sparse Row format.
!           ON RETURN THE COLUMNS OF A ARE PERMUTED.
!
! alph    = diagonal compensation parameter -- the term:
!
!           alph*(sum of all dropped out elements in a given row)
!
!           is added to the diagonal element of U of the factorization
!           Thus: alph = 0 ---> ~ ILU with threshold,
!                 alph = 1 ---> ~ MILU with threshold.
!
! droptol = tolerance used for dropping elements in L and U.
!           elements are dropped if they are .lt. norm(row) x droptol
!           row = row being eliminated
!
! permtol = tolerance ratio used for determning whether to permute
!           two columns.  Two columns are permuted only when
!           abs(a(i,j))*permtol .gt. abs(a(i,i))
!           [0 --> never permute; good values 0.1 to 0.01]
!
! mbloc   = if desired, permuting can be done only within the diagonal
!           blocks of size mbloc. Useful for PDE problems with several
!           degrees of freedom.. If feature not wanted take mbloc=n.
!
! iwk     = integer. The declared lengths of arrays alu and jlu
!           if iwk is not large enough the code will stop prematurely
!           with ierr = -2 or ierr = -3 (see below).
!
! On return:
!===========
!
! alu,jlu = matrix stored in Modified Sparse Row (MSR) format containing
!           the L and U factors together. The diagonal (stored in
!           alu(1:n) ) is inverted. Each i-th row of the alu,jlu matrix
!           contains the i-th row of L (excluding the diagonal entry=1)
!           followed by the i-th row of U.
!
! ju      = integer array of length n containing the pointers to
!           the beginning of each row of U in the matrix alu,jlu.
! iperm   = contains the permutation arrays ..
!           iperm(1:n) = old numbers of unknowns
!           iperm(n+1:2*n) = reverse permutation = new unknowns.
!
! ierr    = integer. Error message with the following meaning.
!           ierr  = 0    --> successful return.
!           ierr .gt. 0  --> zero pivot encountered at step number ierr.
!           ierr  = -1   --> Error. input matrix may be wrong.
!                            (The elimination process has generated a
!                            row in L or U whose length is .gt.  n.)
!           ierr  = -2   --> The L/U matrix overflows the arrays alu,jlu
!           ierr  = -3   --> zero row encountered.
!
! work arrays:
!=============
! jw      = integer work array of length 2*n.
! w       = real work array of length 2*n
!
! Notes:
! ------
! IMPORTANT: TO AVOID PERMUTING THE SOLUTION VECTORS ARRAYS FOR EACH
! LU-SOLVE, THE MATRIX A IS PERMUTED ON RETURN. [all column indices are
! changed]. SIMILARLY FOR THE U MATRIX.
! To permute the matrix back to its original state use the loop:
!
!      do k=ia(1), ia(n+1)-1
!         ja(k) = perm(ja(k))
!      enddo
!
!-----------------------------------------------------------------------
!     local variables
!
   INTEGER :: spag_nextblock_1
   spag_nextblock_1 = 1
   SPAG_DispatchLoop_1: DO
      SELECT CASE (spag_nextblock_1)
      CASE (1)
!-----------------------------------------------------------------------
!     initialize ju0 (points to next element to be added to alu,jlu)
!     and pointer array.
!-----------------------------------------------------------------------
         ju0 = N + 2
         Jlu(1) = ju0
!
!  integer double pointer array.
!
         DO j = 1 , N
            Jw(N+j) = 0
            Iperm(j) = j
            Iperm(N+j) = j
         ENDDO
!-----------------------------------------------------------------------
!     beginning of main loop.
!-----------------------------------------------------------------------
         DO ii = 1 , N
            j1 = Ia(ii)
            j2 = Ia(ii+1) - 1
            dropsum = 0.0D0
            tnorm = 0.0D0
            DO k = j1 , j2
               tnorm = tnorm + abs(A(k))
            ENDDO
            IF ( tnorm==0.0 ) THEN
               spag_nextblock_1 = 4
               CYCLE SPAG_DispatchLoop_1
            ENDIF
            tnorm = tnorm/(j2-j1+1)
!
!     unpack L-part and U-part of row of A in arrays  w  --
!
            lenu = 1
            lenl = 0
            Jw(ii) = ii
            W(ii) = 0.0
            Jw(N+ii) = ii
!
            DO j = j1 , j2
               k = Iperm(N+Ja(j))
               t = A(j)
               IF ( k<ii ) THEN
                  lenl = lenl + 1
                  Jw(lenl) = k
                  W(lenl) = t
                  Jw(N+k) = lenl
               ELSEIF ( k==ii ) THEN
                  W(ii) = t
               ELSE
                  lenu = lenu + 1
                  jpos = ii + lenu - 1
                  Jw(jpos) = k
                  W(jpos) = t
                  Jw(N+k) = jpos
               ENDIF
            ENDDO
            jj = 0
            len = 0
            SPAG_Loop_3_1: DO
!
!     eliminate previous rows
!
               jj = jj + 1
               IF ( jj>lenl ) THEN
!
!     reset double-pointer to zero (U-part)
!
                  DO k = 1 , lenu
                     Jw(N+Jw(ii+k-1)) = 0
                  ENDDO
!
!     update L-matrix
!
                  DO k = 1 , len
                     IF ( ju0>Iwk ) THEN
                        spag_nextblock_1 = 3
                        CYCLE SPAG_DispatchLoop_1
                     ENDIF
                     Alu(ju0) = W(k)
                     Jlu(ju0) = Iperm(Jw(k))
                     ju0 = ju0 + 1
                  ENDDO
!
!     save pointer to beginning of row ii of U
!
                  Ju(ii) = ju0
!
!     update u-matrix -- first apply dropping strategy
!
                  len = 0
                  DO k = 1 , lenu - 1
                     IF ( abs(W(ii+k))>tnorm*Droptol ) THEN
                        len = len + 1
                        W(ii+len) = W(ii+k)
                        Jw(ii+len) = Jw(ii+k)
                     ELSE
                        dropsum = dropsum + W(ii+k)
                     ENDIF
                  ENDDO
!
                  imax = ii
                  xmax = abs(W(imax))
                  xmax0 = xmax
                  icut = ii - 1 + Mbloc - mod(ii-1,Mbloc)
!
!     determine next pivot --
!
                  DO k = ii + 1 , ii + len
                     t = abs(W(k))
                     IF ( t>xmax .AND. t*Permtol>xmax0 .AND. Jw(k)<=icut ) THEN
                        imax = k
                        xmax = t
                     ENDIF
                  ENDDO
!
!     exchange w's
!
                  tmp = W(ii)
                  W(ii) = W(imax)
                  W(imax) = tmp
!
!     update iperm and reverse iperm
!
                  j = Jw(imax)
                  i = Iperm(ii)
                  Iperm(ii) = Iperm(j)
                  Iperm(j) = i
!     reverse iperm
                  Iperm(N+Iperm(ii)) = ii
                  Iperm(N+Iperm(j)) = j
!-----------------------------------------------------------------------
                  IF ( len+ju0-1>Iwk ) THEN
                     spag_nextblock_1 = 3
                     CYCLE SPAG_DispatchLoop_1
                  ENDIF
!
!     copy U-part in original coordinates
!
                  DO k = ii + 1 , ii + len
                     Jlu(ju0) = Iperm(Jw(k))
                     Alu(ju0) = W(k)
                     ju0 = ju0 + 1
                  ENDDO
!
!     define diagonal element
!
                  W(ii) = W(ii) + Alph*dropsum
!
!     store inverse of diagonal element of u
!
                  IF ( W(ii)==0.0 ) W(ii) = (1.0D-4+Droptol)*tnorm
!
                  Alu(ii) = 1.0D0/W(ii)
!
!     update pointer to beginning of next row of U.
!
                  Jlu(ii+1) = ju0
                  EXIT SPAG_Loop_3_1
               ELSE
!-----------------------------------------------------------------------
!     in order to do the elimination in the correct order we must select
!     the smallest column index among jw(k), k=jj+1, ..., lenl.
!-----------------------------------------------------------------------
                  jrow = Jw(jj)
                  k = jj
!
!     determine smallest column index
!
                  DO j = jj + 1 , lenl
                     IF ( Jw(j)<jrow ) THEN
                        jrow = Jw(j)
                        k = j
                     ENDIF
                  ENDDO
!
                  IF ( k/=jj ) THEN
!     exchange in jw
                     j = Jw(jj)
                     Jw(jj) = Jw(k)
                     Jw(k) = j
!     exchange in jr
                     Jw(N+jrow) = jj
                     Jw(N+j) = k
!     exchange in w
                     s = W(jj)
                     W(jj) = W(k)
                     W(k) = s
                  ENDIF
!
!     zero out element in row by resetting jw(n+jrow) to zero.
!
                  Jw(N+jrow) = 0
!
!     drop term if small
!
                  IF ( abs(W(jj))<=Droptol*tnorm ) THEN
                     dropsum = dropsum + W(jj)
                     CYCLE
                  ENDIF
!
!     get the multiplier for row to be eliminated: jrow
!
                  fact = W(jj)*Alu(jrow)
!
!     combine current row and row jrow
!
                  DO k = Ju(jrow) , Jlu(jrow+1) - 1
                     s = fact*Alu(k)
!     new column number
                     j = Iperm(N+Jlu(k))
                     jpos = Jw(N+j)
!
!     if fill-in element is small then disregard:
!
                     IF ( j>=ii ) THEN
!
!     dealing with upper part.
!
                        IF ( jpos==0 ) THEN
!     this is a fill-in element
                           lenu = lenu + 1
                           i = ii + lenu - 1
                           IF ( lenu>N ) THEN
                              spag_nextblock_1 = 2
                              CYCLE SPAG_DispatchLoop_1
                           ENDIF
                           Jw(i) = j
                           Jw(N+j) = i
                           W(i) = -s
                        ELSE
!     no fill-in element --
                           W(jpos) = W(jpos) - s
                        ENDIF
!
!     dealing with lower part.
!
                     ELSEIF ( jpos==0 ) THEN
!     this is a fill-in element
                        lenl = lenl + 1
                        IF ( lenl>N ) THEN
                           spag_nextblock_1 = 2
                           CYCLE SPAG_DispatchLoop_1
                        ENDIF
                        Jw(lenl) = j
                        Jw(N+j) = lenl
                        W(lenl) = -s
                     ELSE
!     no fill-in element --
                        W(jpos) = W(jpos) - s
                     ENDIF
                  ENDDO
                  len = len + 1
                  W(len) = fact
                  Jw(len) = jrow
               ENDIF
            ENDDO SPAG_Loop_3_1
!-----------------------------------------------------------------------
!     end main loop
!-----------------------------------------------------------------------
         ENDDO
!
!     permute all column indices of LU ...
!
         DO k = Jlu(1) , Jlu(N+1) - 1
            Jlu(k) = Iperm(N+Jlu(k))
         ENDDO
!
!     ...and of A
!
         DO k = Ia(1) , Ia(N+1) - 1
            Ja(k) = Iperm(N+Ja(k))
         ENDDO
!
         Ierr = 0
         RETURN
      CASE (2)
!
!     incomprehensible error. Matrix must be wrong.
!
         Ierr = -1
         RETURN
      CASE (3)
!
!     insufficient storage in arrays alu, jlu to store factors
!
         Ierr = -2
         RETURN
      CASE (4)
!
!     zero row encountered
!
         Ierr = -3
         EXIT SPAG_DispatchLoop_1
      END SELECT
   ENDDO SPAG_DispatchLoop_1
!----------------end-of-iludp---------------------------!----------------
!-----------------------------------------------------------------------
END SUBROUTINE iludp
!*==iluk.f90 processed by SPAG 8.04RA 15:47 21 Feb 2024
!!SPAG Open source Personal, Educational or Academic User Clemson University  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
!-----------------------------------------------------------------------
SUBROUTINE iluk(N,A,Ja,Ia,Lfil,Alu,Jlu,Ju,Levs,Iwk,W,Jw,Ierr)
   USE ISO_FORTRAN_ENV                 
   IMPLICIT NONE
!
! Dummy argument declarations rewritten by SPAG
!
   INTEGER , INTENT(IN) :: N
   REAL(REAL64) , INTENT(IN) , DIMENSION(*) :: A
   INTEGER , INTENT(IN) , DIMENSION(*) :: Ja
   INTEGER , INTENT(IN) , DIMENSION(N+1) :: Ia
   INTEGER , INTENT(IN) :: Lfil
   REAL(REAL64) , INTENT(INOUT) , DIMENSION(*) :: Alu
   INTEGER , INTENT(INOUT) , DIMENSION(*) :: Jlu
   INTEGER , INTENT(INOUT) , DIMENSION(N) :: Ju
   INTEGER , INTENT(INOUT) , DIMENSION(*) :: Levs
   INTEGER , INTENT(IN) :: Iwk
   REAL(REAL64) , INTENT(INOUT) , DIMENSION(N) :: W
   INTEGER , INTENT(INOUT) , DIMENSION(3*N) :: Jw
   INTEGER , INTENT(OUT) :: Ierr
!
! Local variable declarations rewritten by SPAG
!
   REAL(REAL64) :: fact , s , t
   INTEGER :: i , ii , j , j1 , j2 , jj , jlev , jpos , jrow , ju0 , k , lenl , lenu , n2
!
! End of declarations rewritten by SPAG
!
!----------------------------------------------------------------------*
!     SPARSKIT ROUTINE ILUK -- ILU WITH LEVEL OF FILL-IN OF K (ILU(k)) *
!----------------------------------------------------------------------*
!
! on entry:
!==========
! n       = integer. The row dimension of the matrix A. The matrix
!
! a,ja,ia = matrix stored in Compressed Sparse Row format.
!
! lfil    = integer. The fill-in parameter. Each element whose
!           leve-of-fill exceeds lfil during the ILU process is dropped.
!           lfil must be .ge. 0
!
! tol     = real*8. Sets the threshold for dropping small terms in the
!           factorization. See below for details on dropping strategy.
!
! iwk     = integer. The minimum length of arrays alu, jlu, and levs.
!
! On return:
!===========
!
! alu,jlu = matrix stored in Modified Sparse Row (MSR) format containing
!           the L and U factors together. The diagonal (stored in
!           alu(1:n) ) is inverted. Each i-th row of the alu,jlu matrix
!           contains the i-th row of L (excluding the diagonal entry=1)
!           followed by the i-th row of U.
!
! ju      = integer array of length n containing the pointers to
!           the beginning of each row of U in the matrix alu,jlu.
!
! levs    = integer (work) array of size iwk -- which contains the
!           levels of each element in alu, jlu.
!
! ierr    = integer. Error message with the following meaning.
!           ierr  = 0    --> successful return.
!           ierr .gt. 0  --> zero pivot encountered at step number ierr.
!           ierr  = -1   --> Error. input matrix may be wrong.
!                            (The elimination process has generated a
!                            row in L or U whose length is .gt.  n.)
!           ierr  = -2   --> The matrix L overflows the array al.
!           ierr  = -3   --> The matrix U overflows the array alu.
!           ierr  = -4   --> Illegal value for lfil.
!           ierr  = -5   --> zero row encountered in A or U.
!
! work arrays:
!=============
! jw      = integer work array of length 3*n.
! w       = real work array of length n
!
! Notes/known bugs: This is not implemented efficiently storage-wise.
!       For example: Only the part of the array levs(*) associated with
!       the U-matrix is needed in the routine.. So some storage can
!       be saved if needed. The levels of fills in the LU matrix are
!       output for information only -- they are not needed by LU-solve.
!
!----------------------------------------------------------------------
! w, ju (1:n) store the working array [1:ii-1 = L-part, ii:n = u]
! jw(n+1:2n)  stores the nonzero indicator.
!
! Notes:
! ------
! All the diagonal elements of the input matrix must be  nonzero.
!
!----------------------------------------------------------------------*
!     locals
   INTEGER :: spag_nextblock_1
   spag_nextblock_1 = 1
   SPAG_DispatchLoop_1: DO
      SELECT CASE (spag_nextblock_1)
      CASE (1)
         IF ( Lfil<0 ) THEN
!
!     illegal lfil entered.
!
            Ierr = -4
            RETURN
         ELSE
!-----------------------------------------------------------------------
!     initialize ju0 (points to next element to be added to alu,jlu)
!     and pointer array.
!-----------------------------------------------------------------------
            n2 = N + N
            ju0 = N + 2
            Jlu(1) = ju0
!
!     initialize nonzero indicator array + levs array --
!
            DO j = 1 , 2*N
               Jw(j) = 0
            ENDDO
!-----------------------------------------------------------------------
!     beginning of main loop.
!-----------------------------------------------------------------------
            DO ii = 1 , N
               j1 = Ia(ii)
               j2 = Ia(ii+1) - 1
!
!     unpack L-part and U-part of row of A in arrays w
!
               lenu = 1
               lenl = 0
               Jw(ii) = ii
               W(ii) = 0.0
               Jw(N+ii) = ii
!
               DO j = j1 , j2
                  k = Ja(j)
                  t = A(j)
                  IF ( t/=0.0 ) THEN
                     IF ( k<ii ) THEN
                        lenl = lenl + 1
                        Jw(lenl) = k
                        W(lenl) = t
                        Jw(n2+lenl) = 0
                        Jw(N+k) = lenl
                     ELSEIF ( k==ii ) THEN
                        W(ii) = t
                        Jw(n2+ii) = 0
                     ELSE
                        lenu = lenu + 1
                        jpos = ii + lenu - 1
                        Jw(jpos) = k
                        W(jpos) = t
                        Jw(n2+jpos) = 0
                        Jw(N+k) = jpos
                     ENDIF
                  ENDIF
               ENDDO
!
               jj = 0
               SPAG_Loop_3_1: DO
!
!     eliminate previous rows
!
                  jj = jj + 1
                  IF ( jj>lenl ) THEN
!
!     reset double-pointer to zero (U-part)
!
                     DO k = 1 , lenu
                        Jw(N+Jw(ii+k-1)) = 0
                     ENDDO
!
!     update l-matrix
!
                     DO k = 1 , lenl
                        IF ( ju0>Iwk ) THEN
                           spag_nextblock_1 = 3
                           CYCLE SPAG_DispatchLoop_1
                        ENDIF
                        IF ( Jw(n2+k)<=Lfil ) THEN
                           Alu(ju0) = W(k)
                           Jlu(ju0) = Jw(k)
                           ju0 = ju0 + 1
                        ENDIF
                     ENDDO
!
!     save pointer to beginning of row ii of U
!
                     Ju(ii) = ju0
!
!     update u-matrix
!
                     DO k = ii + 1 , ii + lenu - 1
                        IF ( ju0>Iwk ) THEN
                           spag_nextblock_1 = 4
                           CYCLE SPAG_DispatchLoop_1
                        ENDIF
                        IF ( Jw(n2+k)<=Lfil ) THEN
                           Jlu(ju0) = Jw(k)
                           Alu(ju0) = W(k)
                           Levs(ju0) = Jw(n2+k)
                           ju0 = ju0 + 1
                        ENDIF
                     ENDDO
 
                     IF ( W(ii)==0.0 ) THEN
                        spag_nextblock_1 = 5
                        CYCLE SPAG_DispatchLoop_1
                     ENDIF
!
                     Alu(ii) = 1.0D0/W(ii)
!
!     update pointer to beginning of next row of U.
!
                     Jlu(ii+1) = ju0
                     EXIT SPAG_Loop_3_1
                  ELSE
!-----------------------------------------------------------------------
!     in order to do the elimination in the correct order we must select
!     the smallest column index among jw(k), k=jj+1, ..., lenl.
!-----------------------------------------------------------------------
                     jrow = Jw(jj)
                     k = jj
!
!     determine smallest column index
!
                     DO j = jj + 1 , lenl
                        IF ( Jw(j)<jrow ) THEN
                           jrow = Jw(j)
                           k = j
                        ENDIF
                     ENDDO
!
                     IF ( k/=jj ) THEN
!     exchange in jw
                        j = Jw(jj)
                        Jw(jj) = Jw(k)
                        Jw(k) = j
!     exchange in jw(n+  (pointers/ nonzero indicator).
                        Jw(N+jrow) = jj
                        Jw(N+j) = k
!     exchange in jw(n2+  (levels)
                        j = Jw(n2+jj)
                        Jw(n2+jj) = Jw(n2+k)
                        Jw(n2+k) = j
!     exchange in w
                        s = W(jj)
                        W(jj) = W(k)
                        W(k) = s
                     ENDIF
!
!     zero out element in row by resetting jw(n+jrow) to zero.
!
                     Jw(N+jrow) = 0
!
!     get the multiplier for row to be eliminated (jrow) + its level
!
                     fact = W(jj)*Alu(jrow)
                     jlev = Jw(n2+jj)
                     IF ( jlev<=Lfil ) THEN
!
!     combine current row and row jrow
!
                        DO k = Ju(jrow) , Jlu(jrow+1) - 1
                           s = fact*Alu(k)
                           j = Jlu(k)
                           jpos = Jw(N+j)
                           IF ( j>=ii ) THEN
!
!     dealing with upper part.
!
                              IF ( jpos==0 ) THEN
!
!     this is a fill-in element
!
                                 lenu = lenu + 1
                                 IF ( lenu>N ) THEN
                                    spag_nextblock_1 = 2
                                    CYCLE SPAG_DispatchLoop_1
                                 ENDIF
                                 i = ii + lenu - 1
                                 Jw(i) = j
                                 Jw(N+j) = i
                                 W(i) = -s
                                 Jw(n2+i) = jlev + Levs(k) + 1
                              ELSE
!
!     this is not a fill-in element
!
                                 W(jpos) = W(jpos) - s
                                 Jw(n2+jpos) = min(Jw(n2+jpos),jlev+Levs(k)+1)
                              ENDIF
!
!     dealing with lower part.
!
                           ELSEIF ( jpos==0 ) THEN
!
!     this is a fill-in element
!
                              lenl = lenl + 1
                              IF ( lenl>N ) THEN
                                 spag_nextblock_1 = 2
                                 CYCLE SPAG_DispatchLoop_1
                              ENDIF
                              Jw(lenl) = j
                              Jw(N+j) = lenl
                              W(lenl) = -s
                              Jw(n2+lenl) = jlev + Levs(k) + 1
                           ELSE
!
!     this is not a fill-in element
!
                              W(jpos) = W(jpos) - s
                              Jw(n2+jpos) = min(Jw(n2+jpos),jlev+Levs(k)+1)
                           ENDIF
                        ENDDO
                        W(jj) = fact
                        Jw(jj) = jrow
                     ENDIF
                  ENDIF
               ENDDO SPAG_Loop_3_1
!-----------------------------------------------------------------------
!     end main loop
!-----------------------------------------------------------------------
            ENDDO
            Ierr = 0
            RETURN
         ENDIF
      CASE (2)
!
!     incomprehensible error. Matrix must be wrong.
!
         Ierr = -1
         RETURN
      CASE (3)
!
!     insufficient storage in L.
!
         Ierr = -2
         RETURN
      CASE (4)
!
!     insufficient storage in U.
!
         Ierr = -3
         RETURN
      CASE (5)
!
!     zero row encountered in A or U.
!
         Ierr = -5
         EXIT SPAG_DispatchLoop_1
      END SELECT
   ENDDO SPAG_DispatchLoop_1
!----------------end-of-iluk--------------------------------------------
!-----------------------------------------------------------------------
END SUBROUTINE iluk
!*==ilu0.f90 processed by SPAG 8.04RA 15:47 21 Feb 2024
!!SPAG Open source Personal, Educational or Academic User Clemson University  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
!----------------------------------------------------------------------
SUBROUTINE ilu0(N,A,Ja,Ia,Alu,Jlu,Ju,Iw,Ierr)
   USE ISO_FORTRAN_ENV                 
   IMPLICIT NONE
!
! Dummy argument declarations rewritten by SPAG
!
   INTEGER , INTENT(IN) :: N
   REAL(REAL64) , INTENT(IN) , DIMENSION(*) :: A
   INTEGER , INTENT(IN) , DIMENSION(*) :: Ja
   INTEGER , INTENT(IN) , DIMENSION(*) :: Ia
   REAL(REAL64) , INTENT(INOUT) , DIMENSION(*) :: Alu
   INTEGER , INTENT(INOUT) , DIMENSION(*) :: Jlu
   INTEGER , INTENT(INOUT) , DIMENSION(*) :: Ju
   INTEGER , INTENT(INOUT) , DIMENSION(*) :: Iw
   INTEGER , INTENT(OUT) :: Ierr
!
! Local variable declarations rewritten by SPAG
!
   INTEGER :: i , ii , j , jcol , jf , jj , jm , jrow , js , ju0 , jw
   REAL(REAL64) :: tl
!
! End of declarations rewritten by SPAG
!
!------------------ right preconditioner ------------------------------*
!                    ***   ilu(0) preconditioner.   ***                *
!----------------------------------------------------------------------*
! Note that this has been coded in such a way that it can be used
! with pgmres. Normally, since the data structure of the L+U matrix is
! the same as that the A matrix, savings can be made. In fact with
! some definitions (not correct for general sparse matrices) all we
! need in addition to a, ja, ia is an additional diagonal.
! ILU0 is not recommended for serious problems. It is only provided
! here for comparison purposes.
!-----------------------------------------------------------------------
!
! on entry:
!---------
! n       = dimension of matrix
! a, ja,
! ia      = original matrix in compressed sparse row storage.
!
! on return:
!-----------
! alu,jlu = matrix stored in Modified Sparse Row (MSR) format containing
!           the L and U factors together. The diagonal (stored in
!           alu(1:n) ) is inverted. Each i-th row of the alu,jlu matrix
!           contains the i-th row of L (excluding the diagonal entry=1)
!           followed by the i-th row of U.
!
! ju         = pointer to the diagonal elements in alu, jlu.
!
! ierr         = integer indicating error code on return
!            ierr = 0 --> normal return
!            ierr = k --> code encountered a zero pivot at step k.
! work arrays:
!-------------
! iw           = integer work array of length n.
!------------
! IMPORTANT
!-----------
! it is assumed that the the elements in the input matrix are stored
!    in such a way that in each row the lower part comes first and
!    then the upper part. To get the correct ILU factorization, it is
!    also necessary to have the elements of L sorted by increasing
!    column number. It may therefore be necessary to sort the
!    elements of a, ja, ia prior to calling ilu0. This can be
!    achieved by transposing the matrix twice using csrcsc.
!
!-----------------------------------------------------------------------
   ju0 = N + 2
   Jlu(1) = ju0
!
! initialize work vector to zero's
!
   DO i = 1 , N
      Iw(i) = 0
   ENDDO
!
! main loop
!
   DO ii = 1 , N
      js = ju0
!
! generating row number ii of L and U.
!
      DO j = Ia(ii) , Ia(ii+1) - 1
!
!     copy row ii of a, ja, ia into row ii of alu, jlu (L/U) matrix.
!
         jcol = Ja(j)
         IF ( jcol==ii ) THEN
            Alu(ii) = A(j)
            Iw(jcol) = ii
            Ju(ii) = ju0
         ELSE
            Alu(ju0) = A(j)
            Jlu(ju0) = Ja(j)
            Iw(jcol) = ju0
            ju0 = ju0 + 1
         ENDIF
      ENDDO
      Jlu(ii+1) = ju0
      jf = ju0 - 1
      jm = Ju(ii) - 1
!
!     exit if diagonal element is reached.
!
      DO j = js , jm
         jrow = Jlu(j)
         tl = Alu(j)*Alu(jrow)
         Alu(j) = tl
!
!     perform  linear combination
!
         DO jj = Ju(jrow) , Jlu(jrow+1) - 1
            jw = Iw(Jlu(jj))
            IF ( jw/=0 ) Alu(jw) = Alu(jw) - tl*Alu(jj)
         ENDDO
      ENDDO
!
!     invert  and store diagonal element.
!
      IF ( Alu(ii)==0.0D0 ) THEN
         CALL spag_block_1
         RETURN
      ENDIF
      Alu(ii) = 1.0D0/Alu(ii)
!
!     reset pointer iw to zero
!
      Iw(ii) = 0
      DO i = js , jf
         Iw(Jlu(i)) = 0
      ENDDO
   ENDDO
   Ierr = 0
   RETURN
CONTAINS
   SUBROUTINE spag_block_1
!
!     zero pivot :
!
      Ierr = ii
   END SUBROUTINE spag_block_1
!
!------- end-of-ilu0 ---------------------------------------------------
!-----------------------------------------------------------------------
END SUBROUTINE ilu0
!*==milu0.f90 processed by SPAG 8.04RA 15:47 21 Feb 2024
!!SPAG Open source Personal, Educational or Academic User Clemson University  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
!----------------------------------------------------------------------
SUBROUTINE milu0(N,A,Ja,Ia,Alu,Jlu,Ju,Iw,Ierr)
   USE ISO_FORTRAN_ENV                 
   IMPLICIT NONE
!
! Dummy argument declarations rewritten by SPAG
!
   INTEGER , INTENT(IN) :: N
   REAL(REAL64) , INTENT(IN) , DIMENSION(*) :: A
   INTEGER , INTENT(IN) , DIMENSION(*) :: Ja
   INTEGER , INTENT(IN) , DIMENSION(*) :: Ia
   REAL(REAL64) , INTENT(INOUT) , DIMENSION(*) :: Alu
   INTEGER , INTENT(INOUT) , DIMENSION(*) :: Jlu
   INTEGER , INTENT(INOUT) , DIMENSION(*) :: Ju
   INTEGER , INTENT(INOUT) , DIMENSION(*) :: Iw
   INTEGER , INTENT(OUT) :: Ierr
!
! Local variable declarations rewritten by SPAG
!
   INTEGER :: i , ii , j , jcol , jf , jj , jm , jrow , js , ju0 , jw
   REAL(REAL64) :: s , tl
!
! End of declarations rewritten by SPAG
!
!----------------------------------------------------------------------*
!                *** simple milu(0) preconditioner. ***                *
!----------------------------------------------------------------------*
! Note that this has been coded in such a way that it can be used
! with pgmres. Normally, since the data structure of a, ja, ia is
! the same as that of a, ja, ia, savings can be made. In fact with
! some definitions (not correct for general sparse matrices) all we
! need in addition to a, ja, ia is an additional diagonal.
! Ilu0 is not recommended for serious problems. It is only provided
! here for comparison purposes.
!-----------------------------------------------------------------------
!
! on entry:
!----------
! n       = dimension of matrix
! a, ja,
! ia      = original matrix in compressed sparse row storage.
!
! on return:
!----------
! alu,jlu = matrix stored in Modified Sparse Row (MSR) format containing
!           the L and U factors together. The diagonal (stored in
!           alu(1:n) ) is inverted. Each i-th row of the alu,jlu matrix
!           contains the i-th row of L (excluding the diagonal entry=1)
!           followed by the i-th row of U.
!
! ju         = pointer to the diagonal elements in alu, jlu.
!
! ierr         = integer indicating error code on return
!            ierr = 0 --> normal return
!            ierr = k --> code encountered a zero pivot at step k.
! work arrays:
!-------------
! iw           = integer work array of length n.
!------------
! Note (IMPORTANT):
!-----------
! it is assumed that the the elements in the input matrix are ordered
!    in such a way that in each row the lower part comes first and
!    then the upper part. To get the correct ILU factorization, it is
!    also necessary to have the elements of L ordered by increasing
!    column number. It may therefore be necessary to sort the
!    elements of a, ja, ia prior to calling milu0. This can be
!    achieved by transposing the matrix twice using csrcsc.
!-----------------------------------------------------------
   ju0 = N + 2
   Jlu(1) = ju0
! initialize work vector to zero's
   DO i = 1 , N
      Iw(i) = 0
   ENDDO
!
!-------------- MAIN LOOP ----------------------------------
!
   DO ii = 1 , N
      js = ju0
!
! generating row number ii or L and U.
!
      DO j = Ia(ii) , Ia(ii+1) - 1
!
!     copy row ii of a, ja, ia into row ii of alu, jlu (L/U) matrix.
!
         jcol = Ja(j)
         IF ( jcol==ii ) THEN
            Alu(ii) = A(j)
            Iw(jcol) = ii
            Ju(ii) = ju0
         ELSE
            Alu(ju0) = A(j)
            Jlu(ju0) = Ja(j)
            Iw(jcol) = ju0
            ju0 = ju0 + 1
         ENDIF
      ENDDO
      Jlu(ii+1) = ju0
      jf = ju0 - 1
      jm = Ju(ii) - 1
!     s accumulates fill-in values
      s = 0.0D0
      DO j = js , jm
         jrow = Jlu(j)
         tl = Alu(j)*Alu(jrow)
         Alu(j) = tl
!-----------------------perform linear combination --------
         DO jj = Ju(jrow) , Jlu(jrow+1) - 1
            jw = Iw(Jlu(jj))
            IF ( jw/=0 ) THEN
               Alu(jw) = Alu(jw) - tl*Alu(jj)
            ELSE
               s = s + tl*Alu(jj)
            ENDIF
         ENDDO
      ENDDO
!----------------------- invert and store diagonal element.
      Alu(ii) = Alu(ii) - s
      IF ( Alu(ii)==0.0D0 ) THEN
         CALL spag_block_1
         RETURN
      ENDIF
      Alu(ii) = 1.0D0/Alu(ii)
!----------------------- reset pointer iw to zero
      Iw(ii) = 0
      DO i = js , jf
         Iw(Jlu(i)) = 0
      ENDDO
   ENDDO
   Ierr = 0
   RETURN
CONTAINS
   SUBROUTINE spag_block_1
!     zero pivot :
      Ierr = ii
   END SUBROUTINE spag_block_1
!------- end-of-milu0 --------------------------------------------------
!-----------------------------------------------------------------------
END SUBROUTINE milu0
!*==pgmres.f90 processed by SPAG 8.04RA 15:47 21 Feb 2024
!!SPAG Open source Personal, Educational or Academic User Clemson University  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
!-----------------------------------------------------------------------
SUBROUTINE pgmres(N,Im,Rhs,Sol,Vv,Eps,Maxits,Iout,Aa,Ja,Ia,Alu,Jlu,Ju,Ierr)
!-----------------------------------------------------------------------
   USE ISO_FORTRAN_ENV                 
   IMPLICIT NONE
!
! PARAMETER definitions rewritten by SPAG
!
   INTEGER , PARAMETER :: KMAX = 50
!
! Dummy argument declarations rewritten by SPAG
!
   INTEGER :: N
   INTEGER , INTENT(IN) :: Im
   REAL(REAL64) , INTENT(INOUT) , DIMENSION(N) :: Rhs
   REAL(REAL64) , INTENT(INOUT) , DIMENSION(N) :: Sol
   REAL(REAL64) , INTENT(INOUT) , DIMENSION(N,*) :: Vv
   REAL(REAL64) , INTENT(IN) :: Eps
   INTEGER , INTENT(IN) :: Maxits
   INTEGER , INTENT(IN) :: Iout
   REAL(REAL64) , DIMENSION(*) :: Aa
   INTEGER , DIMENSION(*) :: Ja
   INTEGER , DIMENSION(N+1) :: Ia
   REAL(REAL64) , DIMENSION(*) :: Alu
   INTEGER , DIMENSION(*) :: Jlu
   INTEGER , DIMENSION(N) :: Ju
   INTEGER , INTENT(OUT) :: Ierr
!
! Local variable declarations rewritten by SPAG
!
   REAL(REAL64) , DIMENSION(KMAX) :: c , s
   REAL(REAL64) , EXTERNAL :: ddot , dnrm2
   REAL(REAL64) :: eps1 , gam , ro
   REAL(REAL64) , SAVE :: epsmac
   REAL(REAL64) , DIMENSION(KMAX+1,KMAX) :: hh
   INTEGER :: i , i1 , ii , its , j , jj , k , k1
   REAL(REAL64) , DIMENSION(KMAX+1) :: rs
   REAL(REAL64) :: t
   EXTERNAL amux , daxpy , lusol
!
! End of declarations rewritten by SPAG
!
!----------------------------------------------------------------------*
!                                                                      *
!                 *** ILUT - Preconditioned GMRES ***                  *
!                                                                      *
!----------------------------------------------------------------------*
! This is a simple version of the ILUT preconditioned GMRES algorithm. *
! The ILUT preconditioner uses a dual strategy for dropping elements   *
! instead  of the usual level of-fill-in approach. See details in ILUT *
! subroutine documentation. PGMRES uses the L and U matrices generated *
! from the subroutine ILUT to precondition the GMRES algorithm.        *
! The preconditioning is applied to the right. The stopping criterion  *
! utilized is based simply on reducing the residual norm by epsilon.   *
! This preconditioning is more reliable than ilu0 but requires more    *
! storage. It seems to be much less prone to difficulties related to   *
! strong nonsymmetries in the matrix. We recommend using a nonzero tol *
! (tol=.005 or .001 usually give good results) in ILUT. Use a large    *
! lfil whenever possible (e.g. lfil = 5 to 10). The higher lfil the    *
! more reliable the code is. Efficiency may also be much improved.     *
! Note that lfil=n and tol=0.0 in ILUT  will yield the same factors as *
! Gaussian elimination without pivoting.                               *
!                                                                      *
! ILU(0) and MILU(0) are also provided for comparison purposes         *
! USAGE: first call ILUT or ILU0 or MILU0 to set up preconditioner and *
! then call pgmres.                                                    *
!----------------------------------------------------------------------*
! Coded by Y. Saad - This version dated May, 7, 1990.                  *
!----------------------------------------------------------------------*
! parameters                                                           *
!-----------                                                           *
! on entry:                                                            *
!==========                                                            *
!                                                                      *
! n     == integer. The dimension of the matrix.                       *
! im    == size of krylov subspace:  should not exceed 50 in this      *
!          version (can be reset by changing parameter command for     *
!          kmax below)                                                 *
! rhs   == real vector of length n containing the right hand side.     *
!          Destroyed on return.                                        *
! sol   == real vector of length n containing an initial guess to the  *
!          solution on input. approximate solution on output           *
! eps   == tolerance for stopping criterion. process is stopped        *
!          as soon as ( ||.|| is the euclidean norm):                  *
!          || current residual||/||initial residual|| <= eps           *
! maxits== maximum number of iterations allowed                        *
! iout  == output unit number number for printing intermediate results *
!          if (iout .le. 0) nothing is printed out.                    *
!                                                                      *
! aa, ja,                                                              *
! ia    == the input matrix in compressed sparse row format:           *
!          aa(1:nnz)  = nonzero elements of A stored row-wise in order *
!          ja(1:nnz) = corresponding column indices.                   *
!          ia(1:n+1) = pointer to beginning of each row in aa and ja.  *
!          here nnz = number of nonzero elements in A = ia(n+1)-ia(1)  *
!                                                                      *
! alu,jlu== A matrix stored in Modified Sparse Row format containing   *
!           the L and U factors, as computed by subroutine ilut.       *
!                                                                      *
! ju     == integer array of length n containing the pointers to       *
!           the beginning of each row of U in alu, jlu as computed     *
!           by subroutine ILUT.                                        *
!                                                                      *
! on return:                                                           *
!==========                                                            *
! sol   == contains an approximate solution (upon successful return).  *
! ierr  == integer. Error message with the following meaning.          *
!          ierr = 0 --> successful return.                             *
!          ierr = 1 --> convergence not achieved in itmax iterations.  *
!          ierr =-1 --> the initial guess seems to be the exact        *
!                       solution (initial residual computed was zero)  *
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
! work arrays:                                                         *
!=============                                                         *
! vv    == work array of length  n x (im+1) (used to store the Arnoli  *
!          basis)                                                      *
!----------------------------------------------------------------------*
! subroutines called :                                                 *
! amux   : SPARSKIT routine to do the matrix by vector multiplication  *
!          delivers y=Ax, given x  -- see SPARSKIT/BLASSM/amux         *
! lusol : combined forward and backward solves (Preconditioning ope.) *
! BLAS1  routines.                                                     *
!----------------------------------------------------------------------*
   INTEGER :: spag_nextblock_1
!-------------------------------------------------------------
! arnoldi size should not exceed kmax=50 in this version..
! to reset modify paramter kmax accordingly.
!-------------------------------------------------------------
   DATA epsmac/1.D-16/
   spag_nextblock_1 = 1
   SPAG_DispatchLoop_1: DO
      SELECT CASE (spag_nextblock_1)
      CASE (1)
         eps1 = 0.0
         its = 0
!-------------------------------------------------------------
! outer loop starts here..
!-------------- compute initial residual vector --------------
         CALL amux(N,Sol,Vv,Aa,Ja,Ia)
         DO j = 1 , N
            Vv(j,1) = Rhs(j) - Vv(j,1)
         ENDDO
         DO
!-------------------------------------------------------------
            ro = dnrm2(N,Vv,1)
            IF ( Iout>0 .AND. its==0 ) WRITE (Iout,99001) its , ro
            IF ( ro==0.0D0 ) THEN
               Ierr = -1
               EXIT SPAG_DispatchLoop_1
            ELSE
               t = 1.0D0/ro
               DO j = 1 , N
                  Vv(j,1) = Vv(j,1)*t
               ENDDO
               IF ( its==0 ) eps1 = Eps*ro
!     ** initialize 1-st term  of rhs of hessenberg system..
               rs(1) = ro
               i = 0
               SPAG_Loop_3_1: DO
                  i = i + 1
                  its = its + 1
                  i1 = i + 1
                  CALL lusol(N,Vv(1,i),Rhs,Alu,Jlu,Ju)
                  CALL amux(N,Rhs,Vv(1,i1),Aa,Ja,Ia)
!-----------------------------------------
!     modified gram - schmidt...
!-----------------------------------------
                  DO j = 1 , i
                     t = ddot(N,Vv(1,j),1,Vv(1,i1),1)
                     hh(j,i) = t
                     CALL daxpy(N,-t,Vv(1,j),1,Vv(1,i1),1)
                  ENDDO
                  t = dnrm2(N,Vv(1,i1),1)
                  hh(i1,i) = t
                  IF ( t/=0.0D0 ) THEN
                     t = 1.0D0/t
                     DO k = 1 , N
                        Vv(k,i1) = Vv(k,i1)*t
                     ENDDO
                  ENDIF
!
!     done with modified gram schimd and arnoldi step..
!     now  update factorization of hh
!
                  IF ( i/=1 ) THEN
!--------perfrom previous transformations  on i-th column of h
                     DO k = 2 , i
                        k1 = k - 1
                        t = hh(k1,i)
                        hh(k1,i) = c(k1)*t + s(k1)*hh(k,i)
                        hh(k,i) = -s(k1)*t + c(k1)*hh(k,i)
                     ENDDO
                  ENDIF
                  gam = sqrt(hh(i,i)**2+hh(i1,i)**2)
!
!     if gamma is zero then any small value will do...
!     will affect only residual estimate
!
                  IF ( gam==0.0D0 ) gam = epsmac
!
!     get  next plane rotation
!
                  c(i) = hh(i,i)/gam
                  s(i) = hh(i1,i)/gam
                  rs(i1) = -s(i)*rs(i)
                  rs(i) = c(i)*rs(i)
!
!     detrermine residual norm and test for convergence-
!
                  hh(i,i) = c(i)*hh(i,i) + s(i)*hh(i1,i)
                  ro = abs(rs(i1))
! 131   format(1h ,2e14.4)
                  IF ( Iout>0 ) WRITE (Iout,99001) its , ro
                  IF ( i>=Im .OR. (ro<=eps1) ) THEN
!
!     now compute solution. first solve upper triangular system.
!
                     rs(i) = rs(i)/hh(i,i)
                     DO ii = 2 , i
                        k = i - ii + 1
                        k1 = k + 1
                        t = rs(k)
                        DO j = k1 , i
                           t = t - hh(k,j)*rs(j)
                        ENDDO
                        rs(k) = t/hh(k,k)
                     ENDDO
!
!     form linear combination of v(*,i)'s to get solution
!
                     t = rs(1)
                     DO k = 1 , N
                        Rhs(k) = Vv(k,1)*t
                     ENDDO
                     DO j = 2 , i
                        t = rs(j)
                        DO k = 1 , N
                           Rhs(k) = Rhs(k) + t*Vv(k,j)
                        ENDDO
                     ENDDO
!
!     call preconditioner.
!
                     CALL lusol(N,Rhs,Rhs,Alu,Jlu,Ju)
                     DO k = 1 , N
                        Sol(k) = Sol(k) + Rhs(k)
                     ENDDO
!
!     restart outer loop  when necessary
!
                     IF ( ro<=eps1 ) THEN
                        Ierr = 0
                        RETURN
                     ELSE
                        IF ( its>=Maxits ) THEN
                           spag_nextblock_1 = 2
                           CYCLE SPAG_DispatchLoop_1
                        ENDIF
!
!     else compute residual vector and continue..
!
                        DO j = 1 , i
                           jj = i1 - j + 1
                           rs(jj-1) = -s(jj-1)*rs(jj)
                           rs(jj) = c(jj-1)*rs(jj)
                        ENDDO
                        DO j = 1 , i1
                           t = rs(j)
                           IF ( j==1 ) t = t - 1.0D0
                           CALL daxpy(N,t,Vv(1,j),1,Vv,1)
                        ENDDO
!     restart outer loop.
                        EXIT SPAG_Loop_3_1
                     ENDIF
                  ENDIF
               ENDDO SPAG_Loop_3_1
            ENDIF
         ENDDO
         spag_nextblock_1 = 2
      CASE (2)
         Ierr = 1
         RETURN
      END SELECT
   ENDDO SPAG_DispatchLoop_1
99001 FORMAT ('   its =',i4,' res. norm =',d20.6)
!-----------------end of pgmres ---------------------------------------
!-----------------------------------------------------------------------
END SUBROUTINE pgmres
!*==lusol.f90 processed by SPAG 8.04RA 15:47 21 Feb 2024
!!SPAG Open source Personal, Educational or Academic User Clemson University  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
!-----------------------------------------------------------------------
SUBROUTINE lusol(N,Y,X,Alu,Jlu,Ju)
   USE ISO_FORTRAN_ENV                 
   IMPLICIT NONE
!
! Dummy argument declarations rewritten by SPAG
!
   INTEGER , INTENT(IN) :: N
   REAL(REAL64) , INTENT(IN) , DIMENSION(N) :: Y
   REAL(REAL64) , INTENT(INOUT) , DIMENSION(N) :: X
   REAL(REAL64) , INTENT(IN) , DIMENSION(*) :: Alu
   INTEGER , INTENT(IN) , DIMENSION(*) :: Jlu
   INTEGER , INTENT(IN) , DIMENSION(*) :: Ju
!
! Local variable declarations rewritten by SPAG
!
   INTEGER :: i , k
!
! End of declarations rewritten by SPAG
!
!-----------------------------------------------------------------------
!
! This routine solves the system (LU) x = y,
! given an LU decomposition of a matrix stored in (alu, jlu, ju)
! modified sparse row format
!
!-----------------------------------------------------------------------
! on entry:
! n   = dimension of system
! y   = the right-hand-side vector
! alu, jlu, ju
!     = the LU matrix as provided from the ILU routines.
!
! on return
! x   = solution of LU x = y.
!-----------------------------------------------------------------------
!
! Note: routine is in place: call lusol (n, x, x, alu, jlu, ju)
!       will solve the system with rhs x and overwrite the result on x .
!
!-----------------------------------------------------------------------
! local variables
!
!
! forward solve
!
   DO i = 1 , N
      X(i) = Y(i)
      DO k = Jlu(i) , Ju(i) - 1
         X(i) = X(i) - Alu(k)*X(Jlu(k))
      ENDDO
   ENDDO
!
!     backward solve.
!
   DO i = N , 1 , -1
      DO k = Ju(i) , Jlu(i+1) - 1
         X(i) = X(i) - Alu(k)*X(Jlu(k))
      ENDDO
      X(i) = Alu(i)*X(i)
   ENDDO
!
!----------------end of lusol ------------------------------------------
!-----------------------------------------------------------------------
END SUBROUTINE lusol
!*==lutsol.f90 processed by SPAG 8.04RA 15:47 21 Feb 2024
!!SPAG Open source Personal, Educational or Academic User Clemson University  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
!-----------------------------------------------------------------------
SUBROUTINE lutsol(N,Y,X,Alu,Jlu,Ju)
   USE ISO_FORTRAN_ENV                 
   IMPLICIT NONE
!
! Dummy argument declarations rewritten by SPAG
!
   INTEGER , INTENT(IN) :: N
   REAL(REAL64) , INTENT(IN) , DIMENSION(N) :: Y
   REAL(REAL64) , INTENT(INOUT) , DIMENSION(N) :: X
   REAL(REAL64) , INTENT(IN) , DIMENSION(*) :: Alu
   INTEGER , INTENT(IN) , DIMENSION(*) :: Jlu
   INTEGER , INTENT(IN) , DIMENSION(*) :: Ju
!
! Local variable declarations rewritten by SPAG
!
   INTEGER :: i , k
!
! End of declarations rewritten by SPAG
!
!-----------------------------------------------------------------------
!
! This routine solves the system  Transp(LU) x = y,
! given an LU decomposition of a matrix stored in (alu, jlu, ju)
! modified sparse row format. Transp(M) is the transpose of M.
!-----------------------------------------------------------------------
! on entry:
! n   = dimension of system
! y   = the right-hand-side vector
! alu, jlu, ju
!     = the LU matrix as provided from the ILU routines.
!
! on return
! x   = solution of transp(LU) x = y.
!-----------------------------------------------------------------------
!
! Note: routine is in place: call lutsol (n, x, x, alu, jlu, ju)
!       will solve the system with rhs x and overwrite the result on x .
!
!-----------------------------------------------------------------------
! local variables
!
!
   DO i = 1 , N
      X(i) = Y(i)
   ENDDO
!
! forward solve (with U^T)
!
   DO i = 1 , N
      X(i) = X(i)*Alu(i)
      DO k = Ju(i) , Jlu(i+1) - 1
         X(Jlu(k)) = X(Jlu(k)) - Alu(k)*X(i)
      ENDDO
   ENDDO
!
!     backward solve (with L^T)
!
   DO i = N , 1 , -1
      DO k = Jlu(i) , Ju(i) - 1
         X(Jlu(k)) = X(Jlu(k)) - Alu(k)*X(i)
      ENDDO
   ENDDO
!
!----------------end of lutsol -----------------------------------------
!-----------------------------------------------------------------------
END SUBROUTINE lutsol
!*==qsplit.f90 processed by SPAG 8.04RA 15:47 21 Feb 2024
!!SPAG Open source Personal, Educational or Academic User Clemson University  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
!-----------------------------------------------------------------------
SUBROUTINE qsplit(A,Ind,N,Ncut)
   USE ISO_FORTRAN_ENV                 
   IMPLICIT NONE
!
! Dummy argument declarations rewritten by SPAG
!
   INTEGER , INTENT(IN) :: N
   REAL(REAL64) , INTENT(INOUT) , DIMENSION(N) :: A
   INTEGER , INTENT(INOUT) , DIMENSION(N) :: Ind
   INTEGER , INTENT(IN) :: Ncut
!
! Local variable declarations rewritten by SPAG
!
   REAL(REAL64) :: abskey , tmp
   INTEGER :: first , itmp , j , last , mid
!
! End of declarations rewritten by SPAG
!
!-----------------------------------------------------------------------
!     does a quick-sort split of a real array.
!     on input a(1:n). is a real array
!     on output a(1:n) is permuted such that its elements satisfy:
!
!     abs(a(i)) .ge. abs(a(ncut)) for i .lt. ncut and
!     abs(a(i)) .le. abs(a(ncut)) for i .gt. ncut
!
!     ind(1:n) is an integer array which permuted in the same way as a(*).
!-----------------------------------------------------------------------
!-----
   first = 1
   last = N
   IF ( Ncut<first .OR. Ncut>last ) RETURN
   DO
!
!     outer loop -- while mid .ne. ncut do
!
      mid = first
      abskey = abs(A(mid))
      DO j = first + 1 , last
         IF ( abs(A(j))>abskey ) THEN
            mid = mid + 1
!     interchange
            tmp = A(mid)
            itmp = Ind(mid)
            A(mid) = A(j)
            Ind(mid) = Ind(j)
            A(j) = tmp
            Ind(j) = itmp
         ENDIF
      ENDDO
!
!     interchange
!
      tmp = A(mid)
      A(mid) = A(first)
      A(first) = tmp
!
      itmp = Ind(mid)
      Ind(mid) = Ind(first)
      Ind(first) = itmp
!
!     test for while loop
!
      IF ( mid==Ncut ) RETURN
      IF ( mid>Ncut ) THEN
         last = mid - 1
      ELSE
         first = mid + 1
      ENDIF
   ENDDO
!----------------end-of-qsplit------------------------------------------
!-----------------------------------------------------------------------
END SUBROUTINE qsplit
