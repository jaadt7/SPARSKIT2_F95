!*==submat.f90 processed by SPAG 8.04RA 15:47 21 Feb 2024
!!SPAG Open source Personal, Educational or Academic User Clemson University  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
!----------------------------------------------------------------------c
!                          S P A R S K I T                             c
!----------------------------------------------------------------------c
!                     UNARY SUBROUTINES MODULE                         c
!----------------------------------------------------------------------c
! contents:                                                            c
!----------                                                            c
! submat : extracts a submatrix from a sparse matrix.                  c
! filter : filters elements from a matrix according to their magnitude.c
! filterm: same as above, but for the MSR format                       c
! csort  : sorts the elements in increasing order of columns           c
! clncsr : clean up the CSR format matrix, remove duplicate entry, etc c
! transp : in-place transposition routine (see also csrcsc in formats) c
! copmat : copy of a matrix into another matrix (both stored csr)      c
! msrcop : copies a matrix in MSR format into a matrix in MSR format   c
! getelm : returns a(i,j) for any (i,j) from a CSR-stored matrix.      c
! getdia : extracts a specified diagonal from a matrix.                c
! getl   : extracts lower triangular part                              c
! getu   : extracts upper triangular part                              c
! levels : gets the level scheduling structure for lower triangular    c
!          matrices.                                                   c
! amask  : extracts     C = A mask M                                   c
! rperm  : permutes the rows of a matrix (B = P A)                     c
! cperm  : permutes the columns of a matrix (B = A Q)                  c
! dperm  : permutes both the rows and columns of a matrix (B = P A Q ) c
! dperm1 : general extractiob routine (extracts arbitrary rows)        c
! dperm2 : general submatrix permutation/extraction routine            c
! dmperm : symmetric permutation of row and column (B=PAP') in MSR fmt c
! dvperm : permutes a real vector (in-place)                           c
! ivperm : permutes an integer vector (in-place)                       c
! retmx  : returns the max absolute value in each row of the matrix    c
! diapos : returns the positions of the diagonal elements in A.        c
! extbdg : extracts the main diagonal blocks of a matrix.              c
! getbwd : returns the bandwidth information on a matrix.              c
! blkfnd : finds the block-size of a matrix.                           c
! blkchk : checks whether a given integer is the block size of A.      c
! infdia : obtains information on the diagonals of A.                  c
! amubdg : gets number of nonzeros in each row of A*B (as well as NNZ) c
! aplbdg : gets number of nonzeros in each row of A+B (as well as NNZ) c
! rnrms  : computes the norms of the rows of A                         c
! cnrms  : computes the norms of the columns of A                      c
! roscal : scales the rows of a matrix by their norms.                 c
! coscal : scales the columns of a matrix by their norms.              c
! addblk : Adds a matrix B into a block of A.                          c
! get1up : Collects the first elements of each row of the upper        c
!          triangular portion of the matrix.                           c
! xtrows : extracts given rows from a matrix in CSR format.            c
! csrkvstr:  Finds block row partitioning of matrix in CSR format      c
! csrkvstc:  Finds block column partitioning of matrix in CSR format   c
! kvstmerge: Merges block partitionings, for conformal row/col pattern c
!----------------------------------------------------------------------c
SUBROUTINE submat(Job,I1,I2,J1,J2,A,Ja,Ia,Nr,Nc,Ao,Jao,Iao)
   USE ISO_FORTRAN_ENV                 
   IMPLICIT NONE
!
! Dummy argument declarations rewritten by SPAG
!
   INTEGER , INTENT(IN) :: Job
   INTEGER , INTENT(IN) :: I1
   INTEGER , INTENT(IN) :: I2
   INTEGER , INTENT(IN) :: J1
   INTEGER , INTENT(IN) :: J2
   REAL(REAL64) , INTENT(IN) , DIMENSION(*) :: A
   INTEGER , INTENT(IN) , DIMENSION(*) :: Ja
   INTEGER , INTENT(IN) , DIMENSION(*) :: Ia
   INTEGER , INTENT(INOUT) :: Nr
   INTEGER , INTENT(INOUT) :: Nc
   REAL(REAL64) , INTENT(OUT) , DIMENSION(*) :: Ao
   INTEGER , INTENT(OUT) , DIMENSION(*) :: Jao
   INTEGER , INTENT(OUT) , DIMENSION(*) :: Iao
!
! Local variable declarations rewritten by SPAG
!
   INTEGER :: i , ii , j , k , k1 , k2 , klen
!
! End of declarations rewritten by SPAG
!
!-----------------------------------------------------------------------
! extracts the submatrix A(i1:i2,j1:j2) and puts the result in
! matrix ao,iao,jao
!---- In place: ao,jao,iao may be the same as a,ja,ia.
!--------------
! on input
!---------
! n	= row dimension of the matrix
! i1,i2 = two integers with i2 .ge. i1 indicating the range of rows to be
!          extracted.
! j1,j2 = two integers with j2 .ge. j1 indicating the range of columns
!         to be extracted.
!         * There is no checking whether the input values for i1, i2, j1,
!           j2 are between 1 and n.
! a,
! ja,
! ia    = matrix in compressed sparse row format.
!
! job	= job indicator: if job .ne. 1 then the real values in a are NOT
!         extracted, only the column indices (i.e. data structure) are.
!         otherwise values as well as column indices are extracted...
!
! on output
!--------------
! nr	= number of rows of submatrix
! nc	= number of columns of submatrix
!	  * if either of nr or nc is nonpositive the code will quit.
!
! ao,
! jao,iao = extracted matrix in general sparse format with jao containing
!	the column indices,and iao being the pointer to the beginning
!	of the row,in arrays a,ja.
!----------------------------------------------------------------------c
!           Y. Saad, Sep. 21 1989                                      c
!----------------------------------------------------------------------c
   Nr = I2 - I1 + 1
   Nc = J2 - J1 + 1
!
   IF ( Nr<=0 .OR. Nc<=0 ) RETURN
!
   klen = 0
!
!     simple procedure. proceeds row-wise...
!
   DO i = 1 , Nr
      ii = I1 + i - 1
      k1 = Ia(ii)
      k2 = Ia(ii+1) - 1
      Iao(i) = klen + 1
!-----------------------------------------------------------------------
      DO k = k1 , k2
         j = Ja(k)
         IF ( j>=J1 .AND. j<=J2 ) THEN
            klen = klen + 1
            IF ( Job==1 ) Ao(klen) = A(k)
            Jao(klen) = j - J1 + 1
         ENDIF
      ENDDO
   ENDDO
   Iao(Nr+1) = klen + 1
!------------end-of submat----------------------------------------------
!-----------------------------------------------------------------------
END SUBROUTINE submat
!*==filter.f90 processed by SPAG 8.04RA 15:47 21 Feb 2024
!!SPAG Open source Personal, Educational or Academic User Clemson University  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
!-----------------------------------------------------------------------
SUBROUTINE filter(N,Job,Drptol,A,Ja,Ia,B,Jb,Ib,Len,Ierr)
   USE ISO_FORTRAN_ENV                 
   IMPLICIT NONE
!
! Dummy argument declarations rewritten by SPAG
!
   INTEGER , INTENT(IN) :: N
   INTEGER , INTENT(IN) :: Job
   REAL(REAL64) , INTENT(IN) :: Drptol
   REAL(REAL64) , INTENT(IN) , DIMENSION(*) :: A
   INTEGER , INTENT(IN) , DIMENSION(*) :: Ja
   INTEGER , INTENT(IN) , DIMENSION(*) :: Ia
   REAL(REAL64) , INTENT(OUT) , DIMENSION(*) :: B
   INTEGER , INTENT(OUT) , DIMENSION(*) :: Jb
   INTEGER , INTENT(OUT) , DIMENSION(*) :: Ib
   INTEGER , INTENT(IN) :: Len
   INTEGER , INTENT(OUT) :: Ierr
!
! Local variable declarations rewritten by SPAG
!
   INTEGER :: index , k , k1 , k2 , row
   REAL(REAL64) :: loctol , norm
!
! End of declarations rewritten by SPAG
!
!-----------------------------------------------------------------------
!     This module removes any elements whose absolute value
!     is small from an input matrix A and puts the resulting
!     matrix in B.  The input parameter job selects a definition
!     of small.
!-----------------------------------------------------------------------
! on entry:
!---------
!  n	 = integer. row dimension of matrix
!  job   = integer. used to determine strategy chosen by caller to
!         drop elements from matrix A.
!          job = 1
!              Elements whose absolute value is less than the
!              drop tolerance are removed.
!          job = 2
!              Elements whose absolute value is less than the
!              product of the drop tolerance and the Euclidean
!              norm of the row are removed.
!          job = 3
!              Elements whose absolute value is less that the
!              product of the drop tolerance and the largest
!              element in the row are removed.
!
! drptol = real. drop tolerance used for dropping strategy.
! a
! ja
! ia     = input matrix in compressed sparse format
! len	 = integer. the amount of space available in arrays b and jb.
!
! on return:
!----------
! b
! jb
! ib    = resulting matrix in compressed sparse format.
!
! ierr	= integer. containing error message.
!         ierr .eq. 0 indicates normal return
!         ierr .gt. 0 indicates that there is'nt enough
!         space is a and ja to store the resulting matrix.
!         ierr then contains the row number where filter stopped.
! note:
!------ This module is in place. (b,jb,ib can ne the same as
!       a, ja, ia in which case the result will be overwritten).
!----------------------------------------------------------------------c
!           contributed by David Day,  Sep 19, 1989.                   c
!----------------------------------------------------------------------c
! local variables
!
   index = 1
   DO row = 1 , N
      k1 = Ia(row)
      k2 = Ia(row+1) - 1
      Ib(row) = index
      IF ( Job==2 ) THEN
         norm = 0.0D0
         DO k = k1 , k2
            norm = norm + A(k)*A(k)
         ENDDO
         norm = sqrt(norm)
      ELSEIF ( Job==3 ) THEN
         norm = 0.0D0
         DO k = k1 , k2
            IF ( abs(A(k))>norm ) norm = abs(A(k))
         ENDDO
      ELSE
         norm = 1.0D0
      ENDIF
      loctol = Drptol*norm
      DO k = k1 , k2
         IF ( abs(A(k))>loctol ) THEN
            IF ( index>Len ) THEN
               Ierr = row
               RETURN
            ENDIF
            B(index) = A(k)
            Jb(index) = Ja(k)
            index = index + 1
         ENDIF
      ENDDO
   ENDDO
   Ib(N+1) = index
!--------------------end-of-filter -------------------------------------
!-----------------------------------------------------------------------
END SUBROUTINE filter
!*==filterm.f90 processed by SPAG 8.04RA 15:47 21 Feb 2024
!!SPAG Open source Personal, Educational or Academic User Clemson University  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
!-----------------------------------------------------------------------
SUBROUTINE filterm(N,Job,Drop,A,Ja,B,Jb,Len,Ierr)
   USE ISO_FORTRAN_ENV                 
   IMPLICIT NONE
!
! Dummy argument declarations rewritten by SPAG
!
   INTEGER , INTENT(IN) :: N
   INTEGER , INTENT(IN) :: Job
   REAL(REAL64) , INTENT(IN) :: Drop
   REAL(REAL64) , INTENT(IN) , DIMENSION(*) :: A
   INTEGER , INTENT(IN) , DIMENSION(*) :: Ja
   REAL(REAL64) , INTENT(OUT) , DIMENSION(*) :: B
   INTEGER , INTENT(OUT) , DIMENSION(*) :: Jb
   INTEGER , INTENT(IN) :: Len
   INTEGER , INTENT(OUT) :: Ierr
!
! Local variable declarations rewritten by SPAG
!
   INTEGER :: index , k , k1 , k2 , row
   REAL(REAL64) :: loctol , norm
!
! End of declarations rewritten by SPAG
!
!-----------------------------------------------------------------------
!     This subroutine removes any elements whose absolute value
!     is small from an input matrix A. Same as filter but
!     uses the MSR format.
!-----------------------------------------------------------------------
! on entry:
!---------
!  n	 = integer. row dimension of matrix
!  job   = integer. used to determine strategy chosen by caller to
!         drop elements from matrix A.
!          job = 1
!              Elements whose absolute value is less than the
!              drop tolerance are removed.
!          job = 2
!              Elements whose absolute value is less than the
!              product of the drop tolerance and the Euclidean
!              norm of the row are removed.
!          job = 3
!              Elements whose absolute value is less that the
!              product of the drop tolerance and the largest
!              element in the row are removed.
!
! drop = real. drop tolerance used for dropping strategy.
! a
! ja     = input matrix in Modifief Sparse Row format
! len	 = integer. the amount of space in arrays b and jb.
!
! on return:
!----------
!
! b, jb = resulting matrix in Modifief Sparse Row format
!
! ierr	= integer. containing error message.
!         ierr .eq. 0 indicates normal return
!         ierr .gt. 0 indicates that there is'nt enough
!         space is a and ja to store the resulting matrix.
!         ierr then contains the row number where filter stopped.
! note:
!------ This module is in place. (b,jb can ne the same as
!       a, ja in which case the result will be overwritten).
!----------------------------------------------------------------------c
!           contributed by David Day,  Sep 19, 1989.                   c
!----------------------------------------------------------------------c
! local variables
!
!
   index = N + 2
   DO row = 1 , N
      k1 = Ja(row)
      k2 = Ja(row+1) - 1
      Jb(row) = index
      IF ( Job==2 ) THEN
         norm = A(row)**2
         DO k = k1 , k2
            norm = norm + A(k)*A(k)
         ENDDO
         norm = sqrt(norm)
      ELSEIF ( Job==3 ) THEN
         norm = abs(A(row))
         DO k = k1 , k2
            norm = max(abs(A(k)),norm)
         ENDDO
      ELSE
         norm = 1.0D0
      ENDIF
      loctol = Drop*norm
      DO k = k1 , k2
         IF ( abs(A(k))>loctol ) THEN
            IF ( index>Len ) THEN
               Ierr = row
               RETURN
            ENDIF
            B(index) = A(k)
            Jb(index) = Ja(k)
            index = index + 1
         ENDIF
      ENDDO
   ENDDO
   Jb(N+1) = index
!--------------------end-of-filterm-------------------------------------
!-----------------------------------------------------------------------
END SUBROUTINE filterm
!*==csort.f90 processed by SPAG 8.04RA 15:47 21 Feb 2024
!!SPAG Open source Personal, Educational or Academic User Clemson University  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
SUBROUTINE csort(N,A,Ja,Ia,Values)
   USE ISO_FORTRAN_ENV                 
   IMPLICIT NONE
!
! Dummy argument declarations rewritten by SPAG
!
   INTEGER , INTENT(IN) :: N
   REAL(REAL64) , INTENT(INOUT) , DIMENSION(*) :: A
   INTEGER , INTENT(INOUT) , DIMENSION(*) :: Ja
   INTEGER , INTENT(IN) , DIMENSION(N+1) :: Ia
   LOGICAL , INTENT(IN) :: Values
!
! Local variable declarations rewritten by SPAG
!
   INTEGER :: i , j , k , row
   REAL(REAL64) :: rj
!
! End of declarations rewritten by SPAG
!
!-----------------------------------------------------------------------
! This routine sorts the elements of  a matrix (stored in Compressed
! Sparse Row Format) in increasing order of their column indices within
! each row. It uses insertion sort
!-----------------------------------------------------------------------
! on entry:
!---------
! n     = the row dimension of the matrix
! a     = the matrix A in compressed sparse row format.
! ja    = the array of column indices of the elements in array a.
! ia    = the array of pointers to the rows.
! values= logical indicating whether or not the real values a(*) must
!         also be permuted. if (.not. values) then the array a is not
!         touched by csort and can be a dummy array.
!
! on return:
!----------
! the matrix stored in the structure a, ja, ia is permuted in such a
! way that the column indices are in increasing order within each row.
! iwork(1:nnz) contains the permutation used  to rearrange the elements.
!-----------------------------------------------------------------------
! Y. Saad - recoded Dec. 20th -  2017 -
!-----------------------------------------------------------------------
! local variables
   INTEGER :: spag_nextblock_1
   rj = 0.0
!
! for each row
!
   DO row = 1 , N
!
! scan row and do an insertion sort
!
      DO k = Ia(row) + 1 , Ia(row+1) - 1
         spag_nextblock_1 = 1
         SPAG_DispatchLoop_1: DO
            SELECT CASE (spag_nextblock_1)
            CASE (1)
               j = Ja(k)
               IF ( Values ) rj = A(k)
               i = k - 1
               DO WHILE ( (i>=Ia(row)) .AND. (j<Ja(i)) )
                  Ja(i+1) = Ja(i)
                  IF ( Values ) A(i+1) = A(i)
                  i = i - 1
                  IF ( i==0 ) THEN
                     spag_nextblock_1 = 2
                     CYCLE SPAG_DispatchLoop_1
                  ENDIF
               ENDDO
               spag_nextblock_1 = 2
            CASE (2)
               Ja(i+1) = j
               IF ( Values ) A(i+1) = rj
               EXIT SPAG_DispatchLoop_1
            END SELECT
         ENDDO SPAG_DispatchLoop_1
      ENDDO
   ENDDO
!---------------end-of-csort--------------------------------------------
!-----------------------------------------------------------------------
END SUBROUTINE csort
!*==clncsr.f90 processed by SPAG 8.04RA 15:47 21 Feb 2024
!!SPAG Open source Personal, Educational or Academic User Clemson University  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
!-----------------------------------------------------------------------
SUBROUTINE clncsr(Job,Value2,Nrow,A,Ja,Ia,Indu,Iwk)
   USE ISO_FORTRAN_ENV                 
   IMPLICIT NONE
!
! Dummy argument declarations rewritten by SPAG
!
   INTEGER , INTENT(IN) :: Nrow
   INTEGER , INTENT(IN) :: Job
   INTEGER , INTENT(IN) :: Value2
   REAL(REAL64) , INTENT(INOUT) , DIMENSION(*) :: A
   INTEGER , INTENT(INOUT) , DIMENSION(*) :: Ja
   INTEGER , INTENT(INOUT) , DIMENSION(Nrow+1) :: Ia
   INTEGER , INTENT(INOUT) , DIMENSION(Nrow) :: Indu
   INTEGER , INTENT(INOUT) , DIMENSION(Nrow+1) :: Iwk
!
! Local variable declarations rewritten by SPAG
!
   INTEGER :: i , ipos , j , k , kfirst , klast , ko
   REAL(REAL64) :: tmp
!
! End of declarations rewritten by SPAG
!
!     .. Scalar Arguments ..
!     ..
!     .. Array Arguments ..
!     ..
!
!     This routine performs two tasks to clean up a CSR matrix
!     -- remove duplicate/zero entries,
!     -- perform a partial ordering, new order lower triangular part,
!        main diagonal, upper triangular part.
!
!     On entry:
!
!     job   = options
!         0 -- nothing is done
!         1 -- eliminate duplicate entries, zero entries.
!         2 -- eliminate duplicate entries and perform partial ordering.
!         3 -- eliminate duplicate entries, sort the entries in the
!              increasing order of clumn indices.
!
!     value2  -- 0 the matrix is pattern only (a is not touched)
!                1 matrix has values too.
!     nrow    -- row dimension of the matrix
!     a,ja,ia -- input matrix in CSR format
!
!     On return:
!     a,ja,ia -- cleaned matrix
!     indu    -- pointers to the beginning of the upper triangular
!                portion if job > 1
!
!     Work space:
!     iwk     -- integer work space of size nrow+1
!
!     .. Local Scalars ..
!     ..
!
   IF ( Job<=0 ) RETURN
!
!     .. eliminate duplicate entries --
!     array INDU is used as marker for existing indices, it is also the
!     location of the entry.
!     IWK is used to stored the old IA array.
!     matrix is copied to squeeze out the space taken by the duplicated
!     entries.
!
   DO i = 1 , Nrow
      Indu(i) = 0
      Iwk(i) = Ia(i)
   ENDDO
   Iwk(Nrow+1) = Ia(Nrow+1)
   k = 1
   DO i = 1 , Nrow
      Ia(i) = k
      ipos = Iwk(i)
      klast = Iwk(i+1)
      SPAG_Loop_2_1: DO
         IF ( ipos<klast ) THEN
            j = Ja(ipos)
            IF ( Indu(j)==0 ) THEN
!     .. new entry ..
               IF ( Value2==0 ) THEN
                  Indu(j) = k
                  Ja(k) = Ja(ipos)
                  k = k + 1
               ELSEIF ( A(ipos)/=0.0D0 ) THEN
                  Indu(j) = k
                  Ja(k) = Ja(ipos)
                  A(k) = A(ipos)
                  k = k + 1
               ENDIF
            ELSEIF ( Value2/=0 ) THEN
!     .. duplicate entry ..
               A(Indu(j)) = A(Indu(j)) + A(ipos)
            ENDIF
            ipos = ipos + 1
            CYCLE
         ENDIF
!     .. remove marks before working on the next row ..
         DO ipos = Ia(i) , k - 1
            Indu(Ja(ipos)) = 0
         ENDDO
         EXIT SPAG_Loop_2_1
      ENDDO SPAG_Loop_2_1
   ENDDO
   Ia(Nrow+1) = k
   IF ( Job<=1 ) RETURN
!
!     .. partial ordering ..
!     split the matrix into strict upper/lower triangular
!     parts, INDU points to the the beginning of the upper part.
!
   DO i = 1 , Nrow
      klast = Ia(i+1) - 1
      kfirst = Ia(i)
      SPAG_Loop_2_2: DO
         IF ( klast>kfirst ) THEN
            IF ( Ja(klast)<i .AND. Ja(kfirst)>=i ) THEN
!     .. swap klast with kfirst ..
               j = Ja(klast)
               Ja(klast) = Ja(kfirst)
               Ja(kfirst) = j
               IF ( Value2/=0 ) THEN
                  tmp = A(klast)
                  A(klast) = A(kfirst)
                  A(kfirst) = tmp
               ENDIF
            ENDIF
            IF ( Ja(klast)>=i ) klast = klast - 1
            IF ( Ja(kfirst)<i ) kfirst = kfirst + 1
            CYCLE
         ENDIF
!
         IF ( Ja(klast)<i ) THEN
            Indu(i) = klast + 1
         ELSE
            Indu(i) = klast
         ENDIF
         EXIT SPAG_Loop_2_2
      ENDDO SPAG_Loop_2_2
   ENDDO
   IF ( Job<=2 ) RETURN
!
!     .. order the entries according to column indices
!     burble-sort is used
!
   DO i = 1 , Nrow
      DO ipos = Ia(i) , Indu(i) - 1
         DO j = Indu(i) - 1 , ipos + 1 , -1
            k = j - 1
            IF ( Ja(k)>Ja(j) ) THEN
               ko = Ja(k)
               Ja(k) = Ja(j)
               Ja(j) = ko
               IF ( Value2/=0 ) THEN
                  tmp = A(k)
                  A(k) = A(j)
                  A(j) = tmp
               ENDIF
            ENDIF
         ENDDO
      ENDDO
      DO ipos = Indu(i) , Ia(i+1) - 1
         DO j = Ia(i+1) - 1 , ipos + 1 , -1
            k = j - 1
            IF ( Ja(k)>Ja(j) ) THEN
               ko = Ja(k)
               Ja(k) = Ja(j)
               Ja(j) = ko
               IF ( Value2/=0 ) THEN
                  tmp = A(k)
                  A(k) = A(j)
                  A(j) = tmp
               ENDIF
            ENDIF
         ENDDO
      ENDDO
   ENDDO
!---- end of clncsr ----------------------------------------------------
!-----------------------------------------------------------------------
END SUBROUTINE clncsr
!*==copmat.f90 processed by SPAG 8.04RA 15:47 21 Feb 2024
!!SPAG Open source Personal, Educational or Academic User Clemson University  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
!-----------------------------------------------------------------------
SUBROUTINE copmat(Nrow,A,Ja,Ia,Ao,Jao,Iao,Ipos,Job)
   USE ISO_FORTRAN_ENV                 
   IMPLICIT NONE
!
! Dummy argument declarations rewritten by SPAG
!
   INTEGER , INTENT(IN) :: Nrow
   REAL(REAL64) , INTENT(IN) , DIMENSION(*) :: A
   INTEGER , INTENT(IN) , DIMENSION(*) :: Ja
   INTEGER , INTENT(IN) , DIMENSION(*) :: Ia
   REAL(REAL64) , INTENT(OUT) , DIMENSION(*) :: Ao
   INTEGER , INTENT(OUT) , DIMENSION(*) :: Jao
   INTEGER , INTENT(OUT) , DIMENSION(*) :: Iao
   INTEGER , INTENT(IN) :: Ipos
   INTEGER , INTENT(IN) :: Job
!
! Local variable declarations rewritten by SPAG
!
   INTEGER :: i , k , kst
!
! End of declarations rewritten by SPAG
!
!----------------------------------------------------------------------
! copies the matrix a, ja, ia, into the matrix ao, jao, iao.
!----------------------------------------------------------------------
! on entry:
!---------
! nrow	= row dimension of the matrix
! a,
! ja,
! ia    = input matrix in compressed sparse row format.
! ipos  = integer. indicates the position in the array ao, jao
!         where the first element should be copied. Thus
!         iao(1) = ipos on return.
! job   = job indicator. if (job .ne. 1) the values are not copies
!         (i.e., pattern only is copied in the form of arrays ja, ia).
!
! on return:
!----------
! ao,
! jao,
! iao   = output matrix containing the same data as a, ja, ia.
!-----------------------------------------------------------------------
!           Y. Saad, March 1990.
!-----------------------------------------------------------------------
! local variables
!
   kst = Ipos - Ia(1)
   DO i = 1 , Nrow + 1
      Iao(i) = Ia(i) + kst
   ENDDO
!
   DO k = Ia(1) , Ia(Nrow+1) - 1
      Jao(kst+k) = Ja(k)
   ENDDO
!
   IF ( Job/=1 ) RETURN
   DO k = Ia(1) , Ia(Nrow+1) - 1
      Ao(kst+k) = A(k)
   ENDDO
!
!--------end-of-copmat -------------------------------------------------
!-----------------------------------------------------------------------
END SUBROUTINE copmat
!*==msrcop.f90 processed by SPAG 8.04RA 15:47 21 Feb 2024
!!SPAG Open source Personal, Educational or Academic User Clemson University  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
!-----------------------------------------------------------------------
SUBROUTINE msrcop(Nrow,A,Ja,Ao,Jao,Job)
   USE ISO_FORTRAN_ENV                 
   IMPLICIT NONE
!
! Dummy argument declarations rewritten by SPAG
!
   INTEGER , INTENT(IN) :: Nrow
   REAL(REAL64) , INTENT(IN) , DIMENSION(*) :: A
   INTEGER , INTENT(IN) , DIMENSION(*) :: Ja
   REAL(REAL64) , INTENT(OUT) , DIMENSION(*) :: Ao
   INTEGER , INTENT(OUT) , DIMENSION(*) :: Jao
   INTEGER , INTENT(IN) :: Job
!
! Local variable declarations rewritten by SPAG
!
   INTEGER :: i , k
!
! End of declarations rewritten by SPAG
!
!----------------------------------------------------------------------
! copies the MSR matrix a, ja, into the MSR matrix ao, jao
!----------------------------------------------------------------------
! on entry:
!---------
! nrow	= row dimension of the matrix
! a,ja  = input matrix in Modified compressed sparse row format.
! job   = job indicator. Values are not copied if job .ne. 1
!
! on return:
!----------
! ao, jao   = output matrix containing the same data as a, ja.
!-----------------------------------------------------------------------
!           Y. Saad,
!-----------------------------------------------------------------------
! local variables
!
   DO i = 1 , Nrow + 1
      Jao(i) = Ja(i)
   ENDDO
!
   DO k = Ja(1) , Ja(Nrow+1) - 1
      Jao(k) = Ja(k)
   ENDDO
!
   IF ( Job/=1 ) RETURN
   DO k = Ja(1) , Ja(Nrow+1) - 1
      Ao(k) = A(k)
   ENDDO
   DO k = 1 , Nrow
      Ao(k) = A(k)
   ENDDO
!
!--------end-of-msrcop -------------------------------------------------
!-----------------------------------------------------------------------
END SUBROUTINE msrcop
!*==getelm.f90 processed by SPAG 8.04RA 15:47 21 Feb 2024
!!SPAG Open source Personal, Educational or Academic User Clemson University  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
!-----------------------------------------------------------------------
FUNCTION getelm(I,J,A,Ja,Ia,Iadd,Sorted)
   USE ISO_FORTRAN_ENV                 
   IMPLICIT NONE
!
! Function and Dummy argument declarations rewritten by SPAG
!
   REAL(REAL64) :: getelm
   INTEGER , INTENT(IN) :: I
   INTEGER , INTENT(IN) :: J
   REAL(REAL64) , INTENT(IN) , DIMENSION(*) :: A
   INTEGER , INTENT(IN) , DIMENSION(*) :: Ja
   INTEGER , INTENT(IN) , DIMENSION(*) :: Ia
   INTEGER , INTENT(INOUT) :: Iadd
   LOGICAL , INTENT(IN) :: Sorted
!
! Local variable declarations rewritten by SPAG
!
   INTEGER :: ibeg , iend , imid , k
!
! End of declarations rewritten by SPAG
!
!-----------------------------------------------------------------------
!     purpose:
!     --------
!     this function returns the element a(i,j) of a matrix a,
!     for any pair (i,j).  the matrix is assumed to be stored
!     in compressed sparse row (csr) format. getelm performs a
!     binary search in the case where it is known that the elements
!     are sorted so that the column indices are in increasing order.
!     also returns (in iadd) the address of the element a(i,j) in
!     arrays a and ja when the search is successsful (zero if not).
!-----
!     first contributed by noel nachtigal (mit).
!     recoded jan. 20, 1991, by y. saad [in particular
!     added handling of the non-sorted case + the iadd output]
!-----------------------------------------------------------------------
!     parameters:
!     -----------
! on entry:
!----------
!     i      = the row index of the element sought (input).
!     j      = the column index of the element sought (input).
!     a      = the matrix a in compressed sparse row format (input).
!     ja     = the array of column indices (input).
!     ia     = the array of pointers to the rows' data (input).
!     sorted = logical indicating whether the matrix is knonw to
!              have its column indices sorted in increasing order
!              (sorted=.true.) or not (sorted=.false.).
!              (input).
! on return:
!-----------
!     getelm = value of a(i,j).
!     iadd   = address of element a(i,j) in arrays a, ja if found,
!              zero if not found. (output)
!
!     note: the inputs i and j are not checked for validity.
!-----------------------------------------------------------------------
!     noel m. nachtigal october 28, 1990 -- youcef saad jan 20, 1991.
!-----------------------------------------------------------------------
!
!     local variables.
!
!
!     initialization
!
   Iadd = 0
   getelm = 0.0
   ibeg = Ia(I)
   iend = Ia(I+1) - 1
!
!     case where matrix is not necessarily sorted
!
   IF ( .NOT.Sorted ) THEN
!
! scan the row - exit as soon as a(i,j) is found
!
      DO k = ibeg , iend
         IF ( Ja(k)==J ) THEN
            Iadd = k
            CALL spag_block_1
            RETURN
         ENDIF
      ENDDO
!
!     end unsorted case. begin sorted case
!
   ELSE
      DO
!
!     begin binary search.   compute the middle index.
!
         imid = (ibeg+iend)/2
!
!     test if  found
!
         IF ( Ja(imid)==J ) THEN
            Iadd = imid
            CALL spag_block_1
            RETURN
         ENDIF
         IF ( ibeg>=iend ) THEN
            CALL spag_block_1
            RETURN
         ENDIF
!
!     else     update the interval bounds.
!
         IF ( Ja(imid)>J ) THEN
            iend = imid - 1
         ELSE
            ibeg = imid + 1
         ENDIF
      ENDDO
!
!     end both cases
!
   ENDIF
   CALL spag_block_1
CONTAINS
   SUBROUTINE spag_block_1
!
      IF ( Iadd/=0 ) getelm = A(Iadd)
   END SUBROUTINE spag_block_1
!
!--------end-of-getelm--------------------------------------------------
!-----------------------------------------------------------------------
END FUNCTION getelm
!*==getdia.f90 processed by SPAG 8.04RA 15:47 21 Feb 2024
!!SPAG Open source Personal, Educational or Academic User Clemson University  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
!-----------------------------------------------------------------------
SUBROUTINE getdia(Nrow,Ncol,Job,A,Ja,Ia,Len,Diag,Idiag,Ioff)
   USE ISO_FORTRAN_ENV                 
   IMPLICIT NONE
!
! Dummy argument declarations rewritten by SPAG
!
   INTEGER , INTENT(IN) :: Nrow
   INTEGER , INTENT(IN) :: Ncol
   INTEGER , INTENT(IN) :: Job
   REAL(REAL64) , INTENT(INOUT) , DIMENSION(*) :: A
   INTEGER , INTENT(INOUT) , DIMENSION(*) :: Ja
   INTEGER , INTENT(INOUT) , DIMENSION(*) :: Ia
   INTEGER , INTENT(INOUT) :: Len
   REAL(REAL64) , INTENT(OUT) , DIMENSION(*) :: Diag
   INTEGER , INTENT(INOUT) , DIMENSION(*) :: Idiag
   INTEGER , INTENT(IN) :: Ioff
!
! Local variable declarations rewritten by SPAG
!
   INTEGER :: i , iend , istart , k , kdiag , ko , kold
!
! End of declarations rewritten by SPAG
!
!-----------------------------------------------------------------------
! this subroutine extracts a given diagonal from a matrix stored in csr
! format. the output matrix may be transformed with the diagonal removed
! from it if desired (as indicated by job.)
!-----------------------------------------------------------------------
! our definition of a diagonal of matrix is a vector of length nrow
! (always) which contains the elements in rows 1 to nrow of
! the matrix that are contained in the diagonal offset by ioff
! with respect to the main diagonal. if the diagonal element
! falls outside the matrix then it is defined as a zero entry.
! thus the proper definition of diag(*) with offset ioff is
!
!     diag(i) = a(i,ioff+i) i=1,2,...,nrow
!     with elements falling outside the matrix being defined as zero.
!
!-----------------------------------------------------------------------
!
! on entry:
!----------
!
! nrow	= integer. the row dimension of the matrix a.
! ncol	= integer. the column dimension of the matrix a.
! job   = integer. job indicator.  if job = 0 then
!         the matrix a, ja, ia, is not altered on return.
!         if job.ne.0  then getdia will remove the entries
!         collected in diag from the original matrix.
!         this is done in place.
!
! a,ja,
!    ia = matrix stored in compressed sparse row a,ja,ia,format
! ioff  = integer,containing the offset of the wanted diagonal
!	  the diagonal extracted is the one corresponding to the
!	  entries a(i,j) with j-i = ioff.
!	  thus ioff = 0 means the main diagonal
!
! on return:
!-----------
! len   = number of nonzero elements found in diag.
!         (len .le. min(nrow,ncol-ioff)-max(1,1-ioff) + 1 )
!
! diag  = real*8 array of length nrow containing the wanted diagonal.
!	  diag contains the diagonal (a(i,j),j-i = ioff ) as defined
!         above.
!
! idiag = integer array of  length len, containing the poisitions
!         in the original arrays a and ja of the diagonal elements
!         collected in diag. a zero entry in idiag(i) means that
!         there was no entry found in row i belonging to the diagonal.
!
! a, ja,
!    ia = if job .ne. 0 the matrix is unchanged. otherwise the nonzero
!         diagonal entries collected in diag are removed from the
!         matrix and therefore the arrays a, ja, ia will change.
!	  (the matrix a, ja, ia will contain len fewer elements)
!
!----------------------------------------------------------------------c
!     Y. Saad, sep. 21 1989 - modified and retested Feb 17, 1996.      c
!----------------------------------------------------------------------c
!     local variables
!
   istart = max(0,-Ioff)
   iend = min(Nrow,Ncol-Ioff)
   Len = 0
   DO i = 1 , Nrow
      Idiag(i) = 0
      Diag(i) = 0.0D0
   ENDDO
!
!     extract  diagonal elements
!
   DO i = istart + 1 , iend
      SPAG_Loop_2_1: DO k = Ia(i) , Ia(i+1) - 1
         IF ( Ja(k)-i==Ioff ) THEN
            Diag(i) = A(k)
            Idiag(i) = k
            Len = Len + 1
            EXIT SPAG_Loop_2_1
         ENDIF
      ENDDO SPAG_Loop_2_1
   ENDDO
   IF ( Job==0 .OR. Len==0 ) RETURN
!
!     remove diagonal elements and rewind structure
!
   ko = 0
   DO i = 1 , Nrow
      kold = ko
      kdiag = Idiag(i)
      DO k = Ia(i) , Ia(i+1) - 1
         IF ( k/=kdiag ) THEN
            ko = ko + 1
            A(ko) = A(k)
            Ja(ko) = Ja(k)
         ENDIF
      ENDDO
      Ia(i) = kold + 1
   ENDDO
!
!     redefine ia(nrow+1)
!
   Ia(Nrow+1) = ko + 1
!------------end-of-getdia----------------------------------------------
!-----------------------------------------------------------------------
END SUBROUTINE getdia
!*==transp.f90 processed by SPAG 8.04RA 15:47 21 Feb 2024
!!SPAG Open source Personal, Educational or Academic User Clemson University  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
!-----------------------------------------------------------------------
SUBROUTINE transp(Nrow,Ncol,A,Ja,Ia,Iwk,Ierr)
   USE ISO_FORTRAN_ENV                 
   IMPLICIT NONE
!
! Dummy argument declarations rewritten by SPAG
!
   INTEGER , INTENT(IN) :: Nrow
   INTEGER , INTENT(INOUT) :: Ncol
   REAL(REAL64) , INTENT(INOUT) , DIMENSION(*) :: A
   INTEGER , INTENT(INOUT) , DIMENSION(*) :: Ja
   INTEGER , INTENT(INOUT) , DIMENSION(*) :: Ia
   INTEGER , INTENT(INOUT) , DIMENSION(*) :: Iwk
   INTEGER , INTENT(OUT) :: Ierr
!
! Local variable declarations rewritten by SPAG
!
   INTEGER :: i , inext , init , j , jcol , k , l , nnz
   REAL(REAL64) :: t , t1
!
! End of declarations rewritten by SPAG
!
!------------------------------------------------------------------------
! In-place transposition routine.
!------------------------------------------------------------------------
! this subroutine transposes a matrix stored in compressed sparse row
! format. the transposition is done in place in that the arrays a,ja,ia
! of the transpose are overwritten onto the original arrays.
!------------------------------------------------------------------------
! on entry:
!---------
! nrow	= integer. The row dimension of A.
! ncol	= integer. The column dimension of A.
! a	= real array of size nnz (number of nonzero elements in A).
!         containing the nonzero elements
! ja	= integer array of length nnz containing the column positions
! 	  of the corresponding elements in a.
! ia	= integer of size n+1, where n = max(nrow,ncol). On entry
!         ia(k) contains the position in a,ja of  the beginning of
!         the k-th row.
!
! iwk	= integer work array of same length as ja.
!
! on return:
!----------
!
! ncol	= actual row dimension of the transpose of the input matrix.
!         Note that this may be .le. the input value for ncol, in
!         case some of the last columns of the input matrix are zero
!         columns. In the case where the actual number of rows found
!         in transp(A) exceeds the input value of ncol, transp will
!         return without completing the transposition. see ierr.
! a,
! ja,
! ia	= contains the transposed matrix in compressed sparse
!         row format. The row dimension of a, ja, ia is now ncol.
!
! ierr	= integer. error message. If the number of rows for the
!         transposed matrix exceeds the input value of ncol,
!         then ierr is  set to that number and transp quits.
!         Otherwise ierr is set to 0 (normal return).
!
! Note:
!----- 1) If you do not need the transposition to be done in place
!         it is preferrable to use the conversion routine csrcsc
!         (see conversion routines in formats).
!      2) the entries of the output matrix are not sorted (the column
!         indices in each are not in increasing order) use csrcsc
!         if you want them sorted.
!----------------------------------------------------------------------c
!           Y. Saad, Sep. 21 1989                                      c
!  modified Oct. 11, 1989.                                             c
!----------------------------------------------------------------------c
! local variables
   INTEGER :: spag_nextblock_1
   spag_nextblock_1 = 1
   SPAG_DispatchLoop_1: DO
      SELECT CASE (spag_nextblock_1)
      CASE (1)
         Ierr = 0
         nnz = Ia(Nrow+1) - 1
!
!     determine column dimension
!
         jcol = 0
         DO k = 1 , nnz
            jcol = max(jcol,Ja(k))
         ENDDO
         IF ( jcol>Ncol ) THEN
            Ierr = jcol
            RETURN
         ENDIF
!
!     convert to coordinate format. use iwk for row indices.
!
         Ncol = jcol
!
         DO i = 1 , Nrow
            DO k = Ia(i) , Ia(i+1) - 1
               Iwk(k) = i
            ENDDO
         ENDDO
!     find pointer array for transpose.
         DO i = 1 , Ncol + 1
            Ia(i) = 0
         ENDDO
         DO k = 1 , nnz
            i = Ja(k)
            Ia(i+1) = Ia(i+1) + 1
         ENDDO
         Ia(1) = 1
!------------------------------------------------------------------------
         DO i = 1 , Ncol
            Ia(i+1) = Ia(i) + Ia(i+1)
         ENDDO
!
!     loop for a cycle in chasing process.
!
         init = 1
         k = 0
         spag_nextblock_1 = 2
      CASE (2)
         t = A(init)
         i = Ja(init)
         j = Iwk(init)
         Iwk(init) = -1
         DO
!------------------------------------------------------------------------
            k = k + 1
!     current row number is i.  determine  where to go.
            l = Ia(i)
!     save the chased element.
            t1 = A(l)
            inext = Ja(l)
!     then occupy its location.
            A(l) = t
            Ja(l) = j
!     update pointer information for next element to be put in row i.
            Ia(i) = l + 1
!     determine  next element to be chased
            IF ( Iwk(l)<0 ) THEN
               DO
                  init = init + 1
                  IF ( init>nnz ) THEN
                     spag_nextblock_1 = 3
                     CYCLE SPAG_DispatchLoop_1
                  ENDIF
                  IF ( Iwk(init)>=0 ) THEN
                     spag_nextblock_1 = 2
                     CYCLE SPAG_DispatchLoop_1
!     restart chasing --
                  ENDIF
               ENDDO
            ELSE
               t = t1
               i = inext
               j = Iwk(l)
               Iwk(l) = -1
               IF ( k>=nnz ) THEN
                  spag_nextblock_1 = 3
                  CYCLE SPAG_DispatchLoop_1
               ENDIF
            ENDIF
         ENDDO
         spag_nextblock_1 = 3
      CASE (3)
         DO i = Ncol , 1 , -1
            Ia(i+1) = Ia(i)
         ENDDO
         Ia(1) = 1
         EXIT SPAG_DispatchLoop_1
      END SELECT
   ENDDO SPAG_DispatchLoop_1
!
!------------------end-of-transp ----------------------------------------
!------------------------------------------------------------------------
END SUBROUTINE transp
!*==getl.f90 processed by SPAG 8.04RA 15:47 21 Feb 2024
!!SPAG Open source Personal, Educational or Academic User Clemson University  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
!------------------------------------------------------------------------
SUBROUTINE getl(N,A,Ja,Ia,Ao,Jao,Iao)
   USE ISO_FORTRAN_ENV                 
   IMPLICIT NONE
!
! Dummy argument declarations rewritten by SPAG
!
   INTEGER , INTENT(IN) :: N
   REAL(REAL64) , INTENT(IN) , DIMENSION(*) :: A
   INTEGER , INTENT(IN) , DIMENSION(*) :: Ja
   INTEGER , INTENT(IN) , DIMENSION(*) :: Ia
   REAL(REAL64) , INTENT(INOUT) , DIMENSION(*) :: Ao
   INTEGER , INTENT(INOUT) , DIMENSION(*) :: Jao
   INTEGER , INTENT(OUT) , DIMENSION(*) :: Iao
!
! Local variable declarations rewritten by SPAG
!
   INTEGER :: i , k , kdiag , ko , kold
   REAL(REAL64) :: t
!
! End of declarations rewritten by SPAG
!
!------------------------------------------------------------------------
! this subroutine extracts the lower triangular part of a matrix
! and writes the result ao, jao, iao. The routine is in place in
! that ao, jao, iao can be the same as a, ja, ia if desired.
!-----------
! on input:
!
! n     = dimension of the matrix a.
! a, ja,
!    ia = matrix stored in compressed sparse row format.
! On return:
! ao, jao,
!    iao = lower triangular matrix (lower part of a)
!	stored in a, ja, ia, format
! note: the diagonal element is the last element in each row.
! i.e. in  a(ia(i+1)-1 )
! ao, jao, iao may be the same as a, ja, ia on entry -- in which case
! getl will overwrite the result on a, ja, ia.
!
!------------------------------------------------------------------------
! local variables
!
! inititialize ko (pointer for output matrix)
!
   ko = 0
   DO i = 1 , N
      kold = ko
      kdiag = 0
      DO k = Ia(i) , Ia(i+1) - 1
         IF ( Ja(k)<=i ) THEN
            ko = ko + 1
            Ao(ko) = A(k)
            Jao(ko) = Ja(k)
            IF ( Ja(k)==i ) kdiag = ko
         ENDIF
      ENDDO
      IF ( kdiag/=0 .AND. kdiag/=ko ) THEN
!
!     exchange
!
         t = Ao(kdiag)
         Ao(kdiag) = Ao(ko)
         Ao(ko) = t
!
         k = Jao(kdiag)
         Jao(kdiag) = Jao(ko)
         Jao(ko) = k
      ENDIF
      Iao(i) = kold + 1
   ENDDO
!     redefine iao(n+1)
   Iao(N+1) = ko + 1
!----------end-of-getl -------------------------------------------------
!-----------------------------------------------------------------------
END SUBROUTINE getl
!*==getu.f90 processed by SPAG 8.04RA 15:47 21 Feb 2024
!!SPAG Open source Personal, Educational or Academic User Clemson University  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
!-----------------------------------------------------------------------
SUBROUTINE getu(N,A,Ja,Ia,Ao,Jao,Iao)
   USE ISO_FORTRAN_ENV                 
   IMPLICIT NONE
!
! Dummy argument declarations rewritten by SPAG
!
   INTEGER , INTENT(IN) :: N
   REAL(REAL64) , INTENT(IN) , DIMENSION(*) :: A
   INTEGER , INTENT(IN) , DIMENSION(*) :: Ja
   INTEGER , INTENT(IN) , DIMENSION(*) :: Ia
   REAL(REAL64) , INTENT(INOUT) , DIMENSION(*) :: Ao
   INTEGER , INTENT(INOUT) , DIMENSION(*) :: Jao
   INTEGER , INTENT(OUT) , DIMENSION(*) :: Iao
!
! Local variable declarations rewritten by SPAG
!
   INTEGER :: i , k , kdiag , kfirst , ko
   REAL(REAL64) :: t
!
! End of declarations rewritten by SPAG
!
!------------------------------------------------------------------------
! this subroutine extracts the upper triangular part of a matrix
! and writes the result ao, jao, iao. The routine is in place in
! that ao, jao, iao can be the same as a, ja, ia if desired.
!-----------
! on input:
!
! n     = dimension of the matrix a.
! a, ja,
!    ia = matrix stored in a, ja, ia, format
! On return:
! ao, jao,
!    iao = upper triangular matrix (upper part of a)
!	stored in compressed sparse row format
! note: the diagonal element is the last element in each row.
! i.e. in  a(ia(i+1)-1 )
! ao, jao, iao may be the same as a, ja, ia on entry -- in which case
! getu will overwrite the result on a, ja, ia.
!
!------------------------------------------------------------------------
! local variables
   ko = 0
   DO i = 1 , N
      kfirst = ko + 1
      kdiag = 0
      DO k = Ia(i) , Ia(i+1) - 1
         IF ( Ja(k)>=i ) THEN
            ko = ko + 1
            Ao(ko) = A(k)
            Jao(ko) = Ja(k)
            IF ( Ja(k)==i ) kdiag = ko
         ENDIF
      ENDDO
      IF ( kdiag/=0 .AND. kdiag/=kfirst ) THEN
!     exchange
         t = Ao(kdiag)
         Ao(kdiag) = Ao(kfirst)
         Ao(kfirst) = t
!
         k = Jao(kdiag)
         Jao(kdiag) = Jao(kfirst)
         Jao(kfirst) = k
      ENDIF
      Iao(i) = kfirst
   ENDDO
!     redefine iao(n+1)
   Iao(N+1) = ko + 1
!----------end-of-getu -------------------------------------------------
!-----------------------------------------------------------------------
END SUBROUTINE getu
!*==levels.f90 processed by SPAG 8.04RA 15:47 21 Feb 2024
!!SPAG Open source Personal, Educational or Academic User Clemson University  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
!-----------------------------------------------------------------------
SUBROUTINE levels(N,Jal,Ial,Nlev,Lev,Ilev,Levnum)
   IMPLICIT NONE
!
! Dummy argument declarations rewritten by SPAG
!
   INTEGER , INTENT(IN) :: N
   INTEGER , INTENT(IN) , DIMENSION(*) :: Jal
   INTEGER , INTENT(IN) , DIMENSION(*) :: Ial
   INTEGER , INTENT(INOUT) :: Nlev
   INTEGER , INTENT(OUT) , DIMENSION(*) :: Lev
   INTEGER , INTENT(INOUT) , DIMENSION(*) :: Ilev
   INTEGER , INTENT(INOUT) , DIMENSION(*) :: Levnum
!
! Local variable declarations rewritten by SPAG
!
   INTEGER :: i , j , levi
!
! End of declarations rewritten by SPAG
!
!-----------------------------------------------------------------------
! levels gets the level structure of a lower triangular matrix
! for level scheduling in the parallel solution of triangular systems
! strict lower matrices (e.g. unit) as well matrices with their main
! diagonal are accepted.
!-----------------------------------------------------------------------
! on entry:
!----------
! n        = integer. The row dimension of the matrix
! jal, ial =
!
! on return:
!-----------
! nlev     = integer. number of levels found
! lev      = integer array of length n containing the level
!            scheduling permutation.
! ilev     = integer array. pointer to beginning of levels in lev.
!            the numbers lev(i) to lev(i+1)-1 contain the row numbers
!            that belong to level number i, in the level scheduling
!            ordering. The equations of the same level can be solved
!            in parallel, once those of all the previous levels have
!            been solved.
! work arrays:
!-------------
! levnum   = integer array of length n (containing the level numbers
!            of each unknown on return)
!-----------------------------------------------------------------------
   DO i = 1 , N
      Levnum(i) = 0
   ENDDO
!
!     compute level of each node --
!
   Nlev = 0
   DO i = 1 , N
      levi = 0
      DO j = Ial(i) , Ial(i+1) - 1
         levi = max(levi,Levnum(Jal(j)))
      ENDDO
      levi = levi + 1
      Levnum(i) = levi
      Nlev = max(Nlev,levi)
   ENDDO
!-------------set data structure  --------------------------------------
   DO j = 1 , Nlev + 1
      Ilev(j) = 0
   ENDDO
!------count  number   of elements in each level -----------------------
   DO j = 1 , N
      i = Levnum(j) + 1
      Ilev(i) = Ilev(i) + 1
   ENDDO
!---- set up pointer for  each  level ----------------------------------
   Ilev(1) = 1
   DO j = 1 , Nlev
      Ilev(j+1) = Ilev(j) + Ilev(j+1)
   ENDDO
!-----determine elements of each level --------------------------------
   DO j = 1 , N
      i = Levnum(j)
      Lev(Ilev(i)) = j
      Ilev(i) = Ilev(i) + 1
   ENDDO
!     reset pointers backwards
   DO j = Nlev , 1 , -1
      Ilev(j+1) = Ilev(j)
   ENDDO
   Ilev(1) = 1
!----------end-of-levels------------------------------------------------
!-----------------------------------------------------------------------
END SUBROUTINE levels
!*==amask.f90 processed by SPAG 8.04RA 15:47 21 Feb 2024
!!SPAG Open source Personal, Educational or Academic User Clemson University  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
!-----------------------------------------------------------------------
SUBROUTINE amask(Nrow,Ncol,A,Ja,Ia,Jmask,Imask,C,Jc,Ic,Iw,Nzmax,Ierr)
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
   INTEGER , INTENT(IN) , DIMENSION(*) :: Jmask
   INTEGER , INTENT(IN) , DIMENSION(Nrow+1) :: Imask
   REAL(REAL64) , INTENT(OUT) , DIMENSION(*) :: C
   INTEGER , INTENT(OUT) , DIMENSION(*) :: Jc
   INTEGER , INTENT(OUT) , DIMENSION(Nrow+1) :: Ic
   LOGICAL , INTENT(INOUT) , DIMENSION(Ncol) :: Iw
   INTEGER , INTENT(IN) :: Nzmax
   INTEGER , INTENT(OUT) :: Ierr
!
! Local variable declarations rewritten by SPAG
!
   INTEGER :: ii , j , k , k1 , k2 , len
!
! End of declarations rewritten by SPAG
!
!---------------------------------------------------------------------
!-----------------------------------------------------------------------
! This subroutine builds a sparse matrix from an input matrix by
! extracting only elements in positions defined by the mask jmask, imask
!-----------------------------------------------------------------------
! On entry:
!---------
! nrow  = integer. row dimension of input matrix
! ncol	= integer. Column dimension of input matrix.
!
! a,
! ja,
! ia	= matrix in Compressed Sparse Row format
!
! jmask,
! imask = matrix defining mask (pattern only) stored in compressed
!         sparse row format.
!
! nzmax = length of arrays c and jc. see ierr.
!
! On return:
!-----------
!
! a, ja, ia and jmask, imask are unchanged.
!
! c
! jc,
! ic	= the output matrix in Compressed Sparse Row format.
!
! ierr  = integer. serving as error message.c
!         ierr = 1  means normal return
!         ierr .gt. 1 means that amask stopped when processing
!         row number ierr, because there was not enough space in
!         c, jc according to the value of nzmax.
!
! work arrays:
!-------------
! iw	= logical work array of length ncol.
!
! note:
!------ the  algorithm is in place: c, jc, ic can be the same as
! a, ja, ia in which cas the code will overwrite the matrix c
! on a, ja, ia
!
!-----------------------------------------------------------------------
   Ierr = 0
   len = 0
   DO j = 1 , Ncol
      Iw(j) = .FALSE.
   ENDDO
!     unpack the mask for row ii in iw
   DO ii = 1 , Nrow
!     save pointer in order to be able to do things in place
      DO k = Imask(ii) , Imask(ii+1) - 1
         Iw(Jmask(k)) = .TRUE.
      ENDDO
!     add umasked elemnts of row ii
      k1 = Ia(ii)
      k2 = Ia(ii+1) - 1
      Ic(ii) = len + 1
      DO k = k1 , k2
         j = Ja(k)
         IF ( Iw(j) ) THEN
            len = len + 1
            IF ( len>Nzmax ) THEN
               Ierr = ii
               RETURN
            ENDIF
            Jc(len) = j
            C(len) = A(k)
         ENDIF
      ENDDO
!
      DO k = Imask(ii) , Imask(ii+1) - 1
         Iw(Jmask(k)) = .FALSE.
      ENDDO
   ENDDO
   Ic(Nrow+1) = len + 1
!
!-----end-of-amask -----------------------------------------------------
!-----------------------------------------------------------------------
END SUBROUTINE amask
!*==rperm.f90 processed by SPAG 8.04RA 15:47 21 Feb 2024
!!SPAG Open source Personal, Educational or Academic User Clemson University  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
!-----------------------------------------------------------------------
SUBROUTINE rperm(Nrow,A,Ja,Ia,Ao,Jao,Iao,Perm,Job)
   USE ISO_FORTRAN_ENV                 
   IMPLICIT NONE
!
! Dummy argument declarations rewritten by SPAG
!
   INTEGER , INTENT(IN) :: Nrow
   REAL(REAL64) , INTENT(IN) , DIMENSION(*) :: A
   INTEGER , INTENT(IN) , DIMENSION(*) :: Ja
   INTEGER , INTENT(IN) , DIMENSION(Nrow+1) :: Ia
   REAL(REAL64) , INTENT(OUT) , DIMENSION(*) :: Ao
   INTEGER , INTENT(OUT) , DIMENSION(*) :: Jao
   INTEGER , INTENT(INOUT) , DIMENSION(Nrow+1) :: Iao
   INTEGER , INTENT(IN) , DIMENSION(Nrow) :: Perm
   INTEGER , INTENT(IN) :: Job
!
! Local variable declarations rewritten by SPAG
!
   INTEGER :: i , ii , j , k , ko
   LOGICAL :: values
!
! End of declarations rewritten by SPAG
!
!-----------------------------------------------------------------------
! this subroutine permutes the rows of a matrix in CSR format.
! rperm  computes B = P A  where P is a permutation matrix.
! the permutation P is defined through the array perm: for each j,
! perm(j) represents the destination row number of row number j.
! Youcef Saad -- recoded Jan 28, 1991.
!-----------------------------------------------------------------------
! on entry:
!----------
! n 	= dimension of the matrix
! a, ja, ia = input matrix in csr format
! perm 	= integer array of length nrow containing the permutation arrays
!	  for the rows: perm(i) is the destination of row i in the
!         permuted matrix.
!         ---> a(i,j) in the original matrix becomes a(perm(i),j)
!         in the output  matrix.
!
! job	= integer indicating the work to be done:
! 		job = 1	permute a, ja, ia into ao, jao, iao
!                       (including the copying of real values ao and
!                       the array iao).
! 		job .ne. 1 :  ignore real values.
!                     (in which case arrays a and ao are not needed nor
!                      used).
!
!------------
! on return:
!------------
! ao, jao, iao = input matrix in a, ja, ia format
! note :
!        if (job.ne.1)  then the arrays a and ao are not used.
!----------------------------------------------------------------------c
!           Y. Saad, May  2, 1990                                      c
!----------------------------------------------------------------------c
   values = (Job==1)
!
!     determine pointers for output matix.
!
   DO j = 1 , Nrow
      i = Perm(j)
      Iao(i+1) = Ia(j+1) - Ia(j)
   ENDDO
!
! get pointers from lengths
!
   Iao(1) = 1
   DO j = 1 , Nrow
      Iao(j+1) = Iao(j+1) + Iao(j)
   ENDDO
!
! copying
!
   DO ii = 1 , Nrow
!
! old row = ii  -- new row = iperm(ii) -- ko = new pointer
!
      ko = Iao(Perm(ii))
      DO k = Ia(ii) , Ia(ii+1) - 1
         Jao(ko) = Ja(k)
         IF ( values ) Ao(ko) = A(k)
         ko = ko + 1
      ENDDO
   ENDDO
!
!---------end-of-rperm -------------------------------------------------
!-----------------------------------------------------------------------
END SUBROUTINE rperm
!*==cperm.f90 processed by SPAG 8.04RA 15:47 21 Feb 2024
!!SPAG Open source Personal, Educational or Academic User Clemson University  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
!-----------------------------------------------------------------------
SUBROUTINE cperm(Nrow,A,Ja,Ia,Ao,Jao,Iao,Perm,Job)
   USE ISO_FORTRAN_ENV                 
   IMPLICIT NONE
!
! Dummy argument declarations rewritten by SPAG
!
   INTEGER , INTENT(IN) :: Nrow
   REAL(REAL64) , INTENT(IN) , DIMENSION(*) :: A
   INTEGER , INTENT(IN) , DIMENSION(*) :: Ja
   INTEGER , INTENT(IN) , DIMENSION(Nrow+1) :: Ia
   REAL(REAL64) , INTENT(OUT) , DIMENSION(*) :: Ao
   INTEGER , INTENT(OUT) , DIMENSION(*) :: Jao
   INTEGER , INTENT(OUT) , DIMENSION(Nrow+1) :: Iao
   INTEGER , INTENT(IN) , DIMENSION(*) :: Perm
   INTEGER , INTENT(IN) :: Job
!
! Local variable declarations rewritten by SPAG
!
   INTEGER :: i , k , nnz
!
! End of declarations rewritten by SPAG
!
!-----------------------------------------------------------------------
! this subroutine permutes the columns of a matrix a, ja, ia.
! the result is written in the output matrix  ao, jao, iao.
! cperm computes B = A P, where  P is a permutation matrix
! that maps column j into column perm(j), i.e., on return
!      a(i,j) becomes a(i,perm(j)) in new matrix
! Y. Saad, May 2, 1990 / modified Jan. 28, 1991.
!-----------------------------------------------------------------------
! on entry:
!----------
! nrow 	= row dimension of the matrix
!
! a, ja, ia = input matrix in csr format.
!
! perm	= integer array of length ncol (number of columns of A
!         containing the permutation array  the columns:
!         a(i,j) in the original matrix becomes a(i,perm(j))
!         in the output matrix.
!
! job	= integer indicating the work to be done:
! 		job = 1	permute a, ja, ia into ao, jao, iao
!                       (including the copying of real values ao and
!                       the array iao).
! 		job .ne. 1 :  ignore real values ao and ignore iao.
!
!------------
! on return:
!------------
! ao, jao, iao = input matrix in a, ja, ia format (array ao not needed)
!
! Notes:
!-------
! 1. if job=1 then ao, iao are not used.
! 2. This routine is in place: ja, jao can be the same.
! 3. If the matrix is initially sorted (by increasing column number)
!    then ao,jao,iao  may not be on return.
!
!----------------------------------------------------------------------c
! local parameters:
!
   nnz = Ia(Nrow+1) - 1
   DO k = 1 , nnz
      Jao(k) = Perm(Ja(k))
   ENDDO
!
!     done with ja array. return if no need to touch values.
!
   IF ( Job/=1 ) RETURN
!
! else get new pointers -- and copy values too.
!
   DO i = 1 , Nrow + 1
      Iao(i) = Ia(i)
   ENDDO
!
   DO k = 1 , nnz
      Ao(k) = A(k)
   ENDDO
!
!---------end-of-cperm--------------------------------------------------
!-----------------------------------------------------------------------
END SUBROUTINE cperm
!*==dperm.f90 processed by SPAG 8.04RA 15:47 21 Feb 2024
!!SPAG Open source Personal, Educational or Academic User Clemson University  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
!-----------------------------------------------------------------------
SUBROUTINE dperm(Nrow,A,Ja,Ia,Ao,Jao,Iao,Perm,Qperm,Job)
   USE ISO_FORTRAN_ENV                 
   IMPLICIT NONE
!
! Dummy argument declarations rewritten by SPAG
!
   INTEGER :: Nrow
   REAL(REAL64) , DIMENSION(*) :: A
   INTEGER , DIMENSION(*) :: Ja
   INTEGER , DIMENSION(Nrow+1) :: Ia
   REAL(REAL64) , DIMENSION(*) :: Ao
   INTEGER , DIMENSION(*) :: Jao
   INTEGER , DIMENSION(Nrow+1) :: Iao
   INTEGER , DIMENSION(Nrow) :: Perm
   INTEGER , DIMENSION(*) :: Qperm
   INTEGER , INTENT(IN) :: Job
!
! Local variable declarations rewritten by SPAG
!
   INTEGER :: locjob
   EXTERNAL cperm , rperm
!
! End of declarations rewritten by SPAG
!
!-----------------------------------------------------------------------
! This routine permutes the rows and columns of a matrix stored in CSR
! format. i.e., it computes P A Q, where P, Q are permutation matrices.
! P maps row i into row perm(i) and Q maps column j into column qperm(j):
!      a(i,j)    becomes   a(perm(i),qperm(j)) in new matrix
! In the particular case where Q is the transpose of P (symmetric
! permutation of A) then qperm is not needed.
! note that qperm should be of length ncol (number of columns) but this
! is not checked.
!-----------------------------------------------------------------------
! Y. Saad, Sep. 21 1989 / recoded Jan. 28 1991.
!-----------------------------------------------------------------------
! on entry:
!----------
! n 	= dimension of the matrix
! a, ja,
!    ia = input matrix in a, ja, ia format
! perm 	= integer array of length n containing the permutation arrays
!	  for the rows: perm(i) is the destination of row i in the
!         permuted matrix -- also the destination of column i in case
!         permutation is symmetric (job .le. 2)
!
! qperm	= same thing for the columns. This should be provided only
!         if job=3 or job=4, i.e., only in the case of a nonsymmetric
!	  permutation of rows and columns. Otherwise qperm is a dummy
!
! job	= integer indicating the work to be done:
! * job = 1,2 permutation is symmetric  Ao :== P * A * transp(P)
! 		job = 1	permute a, ja, ia into ao, jao, iao
! 		job = 2 permute matrix ignoring real values.
! * job = 3,4 permutation is non-symmetric  Ao :== P * A * Q
! 		job = 3	permute a, ja, ia into ao, jao, iao
! 		job = 4 permute matrix ignoring real values.
!
! on return:
!-----------
! ao, jao, iao = input matrix in a, ja, ia format
!
! in case job .eq. 2 or job .eq. 4, a and ao are never referred to
! and can be dummy arguments.
! Notes:
!-------
!  1) algorithm is in place
!  2) column indices may not be sorted on return even  though they may be
!     on entry.
!----------------------------------------------------------------------c
! local variables
!
!     locjob indicates whether or not real values must be copied.
!
   locjob = mod(Job,2)
!
! permute rows first
!
   CALL rperm(Nrow,A,Ja,Ia,Ao,Jao,Iao,Perm,locjob)
!
! then permute columns
!
   locjob = 0
!
   IF ( Job<=2 ) THEN
      CALL cperm(Nrow,Ao,Jao,Iao,Ao,Jao,Iao,Perm,locjob)
   ELSE
      CALL cperm(Nrow,Ao,Jao,Iao,Ao,Jao,Iao,Qperm,locjob)
   ENDIF
!
!-------end-of-dperm----------------------------------------------------
!-----------------------------------------------------------------------
END SUBROUTINE dperm
!*==dperm1.f90 processed by SPAG 8.04RA 15:47 21 Feb 2024
!!SPAG Open source Personal, Educational or Academic User Clemson University  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
!-----------------------------------------------------------------------
SUBROUTINE dperm1(I1,I2,A,Ja,Ia,B,Jb,Ib,Perm,Ipos,Job)
   USE ISO_FORTRAN_ENV                 
   IMPLICIT NONE
!
! Dummy argument declarations rewritten by SPAG
!
   INTEGER , INTENT(IN) :: I1
   INTEGER , INTENT(IN) :: I2
   REAL(REAL64) , INTENT(IN) , DIMENSION(*) :: A
   INTEGER , INTENT(IN) , DIMENSION(*) :: Ja
   INTEGER , INTENT(IN) , DIMENSION(*) :: Ia
   REAL(REAL64) , INTENT(OUT) , DIMENSION(*) :: B
   INTEGER , INTENT(OUT) , DIMENSION(*) :: Jb
   INTEGER , INTENT(OUT) , DIMENSION(*) :: Ib
   INTEGER , INTENT(IN) , DIMENSION(*) :: Perm
   INTEGER , INTENT(IN) :: Ipos
   INTEGER , INTENT(IN) :: Job
!
! Local variable declarations rewritten by SPAG
!
   INTEGER :: i , irow , k , ko
   LOGICAL :: values
!
! End of declarations rewritten by SPAG
!
!-----------------------------------------------------------------------
!     general submatrix extraction routine.
!-----------------------------------------------------------------------
!     extracts rows perm(i1), perm(i1+1), ..., perm(i2) (in this order)
!     from a matrix (doing nothing in the column indices.) The resulting
!     submatrix is constructed in b, jb, ib. A pointer ipos to the
!     beginning of arrays b,jb,is also allowed (i.e., nonzero elements
!     are accumulated starting in position ipos of b, jb).
!-----------------------------------------------------------------------
! Y. Saad,Sep. 21 1989 / recoded Jan. 28 1991 / modified for PSPARSLIB
! Sept. 1997..
!-----------------------------------------------------------------------
! on entry:
!----------
! n 	= dimension of the matrix
! a,ja,
!   ia  = input matrix in CSR format
! perm 	= integer array of length n containing the indices of the rows
!         to be extracted.
!
! job   = job indicator. if (job .ne.1) values are not copied (i.e.,
!         only pattern is copied).
!
! on return:
!-----------
! b,ja,
! ib   = matrix in csr format. b(ipos:ipos+nnz-1),jb(ipos:ipos+nnz-1)
!     contain the value and column indices respectively of the nnz
!     nonzero elements of the permuted matrix. thus ib(1)=ipos.
!
! Notes:
!-------
!  algorithm is NOT in place
!-----------------------------------------------------------------------
! local variables
!
!-----------------------------------------------------------------------
   values = (Job==1)
   ko = Ipos
   Ib(1) = ko
   DO i = I1 , I2
      irow = Perm(i)
      DO k = Ia(irow) , Ia(irow+1) - 1
         IF ( values ) B(ko) = A(k)
         Jb(ko) = Ja(k)
         ko = ko + 1
      ENDDO
      Ib(i-I1+2) = ko
   ENDDO
!--------end-of-dperm1--------------------------------------------------
!-----------------------------------------------------------------------
END SUBROUTINE dperm1
!*==dperm2.f90 processed by SPAG 8.04RA 15:47 21 Feb 2024
!!SPAG Open source Personal, Educational or Academic User Clemson University  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
!-----------------------------------------------------------------------
SUBROUTINE dperm2(I1,I2,A,Ja,Ia,B,Jb,Ib,Cperm,Rperm,Istart,Ipos,Job)
   USE ISO_FORTRAN_ENV                 
   IMPLICIT NONE
!
! Dummy argument declarations rewritten by SPAG
!
   INTEGER , INTENT(IN) :: I1
   INTEGER , INTENT(IN) :: I2
   REAL(REAL64) , INTENT(IN) , DIMENSION(*) :: A
   INTEGER , INTENT(IN) , DIMENSION(*) :: Ja
   INTEGER , INTENT(IN) , DIMENSION(*) :: Ia
   REAL(REAL64) , INTENT(OUT) , DIMENSION(*) :: B
   INTEGER , INTENT(OUT) , DIMENSION(*) :: Jb
   INTEGER , INTENT(OUT) , DIMENSION(*) :: Ib
   INTEGER , INTENT(IN) , DIMENSION(*) :: Cperm
   INTEGER , INTENT(IN) , DIMENSION(*) :: Rperm
   INTEGER , INTENT(IN) :: Istart
   INTEGER , INTENT(IN) :: Ipos
   INTEGER , INTENT(IN) :: Job
!
! Local variable declarations rewritten by SPAG
!
   INTEGER :: i , irow , k , ko
   LOGICAL :: values
!
! End of declarations rewritten by SPAG
!
!-----------------------------------------------------------------------
!     general submatrix permutation/ extraction routine.
!-----------------------------------------------------------------------
!     extracts rows rperm(i1), rperm(i1+1), ..., rperm(i2) and does an
!     associated column permutation (using array cperm). The resulting
!     submatrix is constructed in b, jb, ib. For added flexibility, the
!     extracted elements are put in sequence starting from row 'istart'
!     of B. In addition a pointer ipos to the beginning of arrays b,jb,
!     is also allowed (i.e., nonzero elements are accumulated starting in
!     position ipos of b, jb). In most applications istart and ipos are
!     equal to one. However, the generality adds substantial flexiblity.
!     EXPLE: (1) to permute msr to msr (excluding diagonals)
!     call dperm2 (1,n,a,ja,ja,b,jb,jb,rperm,rperm,1,n+2)
!            (2) To extract rows 1 to 10: define rperm and cperm to be
!     identity permutations (rperm(i)=i, i=1,n) and then
!            call dperm2 (1,10,a,ja,ia,b,jb,ib,rperm,rperm,1,1)
!            (3) to achieve a symmetric permutation as defined by perm:
!            call dperm2 (1,10,a,ja,ia,b,jb,ib,perm,perm,1,1)
!            (4) to get a symmetric permutation of A and append the
!            resulting data structure to A's data structure (useful!)
!            call dperm2 (1,10,a,ja,ia,a,ja,ia(n+1),perm,perm,1,ia(n+1))
!-----------------------------------------------------------------------
! Y. Saad,Sep. 21 1989 / recoded Jan. 28 1991.
!-----------------------------------------------------------------------
! on entry:
!----------
! n 	= dimension of the matrix
! i1,i2 = extract rows rperm(i1) to rperm(i2) of A, with i1<i2.
!
! a,ja,
!   ia  = input matrix in CSR format
! cperm = integer array of length n containing the permutation arrays
!	  for the columns: cperm(i) is the destination of column j,
!         i.e., any column index ja(k) is transformed into cperm(ja(k))
!
! rperm	=  permutation array for the rows. rperm(i) = origin (in A) of
!          row i in B. This is the reverse permutation relative to the
!          ones used in routines cperm, dperm,....
!          rows rperm(i1), rperm(i1)+1, ... rperm(i2) are
!          extracted from A and stacked into B, starting in row istart
!          of B.
! istart= starting row for B where extracted matrix is to be added.
!         this is also only a pointer of the be beginning address for
!         ib , on return.
! ipos  = beginning position in arrays b and jb where to start copying
!         elements. Thus, ib(istart) = ipos.
!
! job   = job indicator. if (job .ne.1) values are not copied (i.e.,
!         only pattern is copied).
!
! on return:
!-----------
! b,ja,
! ib   = matrix in csr format. positions 1,2,...,istart-1 of ib
!     are not touched. b(ipos:ipos+nnz-1),jb(ipos:ipos+nnz-1)
!     contain the value and column indices respectively of the nnz
!     nonzero elements of the permuted matrix. thus ib(istart)=ipos.
!
! Notes:
!-------
!  1) algorithm is NOT in place
!  2) column indices may not be sorted on return even  though they
!     may be on entry.
!-----------------------------------------------------------------------
! local variables
!
!-----------------------------------------------------------------------
   values = (Job==1)
   ko = Ipos
   Ib(Istart) = ko
   DO i = I1 , I2
      irow = Rperm(i)
      DO k = Ia(irow) , Ia(irow+1) - 1
         IF ( values ) B(ko) = A(k)
         Jb(ko) = Cperm(Ja(k))
         ko = ko + 1
      ENDDO
      Ib(Istart+i-I1+1) = ko
   ENDDO
!--------end-of-dperm2--------------------------------------------------
!-----------------------------------------------------------------------
END SUBROUTINE dperm2
!*==dmperm.f90 processed by SPAG 8.04RA 15:47 21 Feb 2024
!!SPAG Open source Personal, Educational or Academic User Clemson University  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
!-----------------------------------------------------------------------
SUBROUTINE dmperm(Nrow,A,Ja,Ao,Jao,Perm,Job)
   USE ISO_FORTRAN_ENV                 
   IMPLICIT NONE
!
! Dummy argument declarations rewritten by SPAG
!
   INTEGER :: Nrow
   REAL(REAL64) , DIMENSION(*) :: A
   INTEGER , DIMENSION(*) :: Ja
   REAL(REAL64) , DIMENSION(*) :: Ao
   INTEGER , INTENT(INOUT) , DIMENSION(*) :: Jao
   INTEGER , DIMENSION(Nrow) :: Perm
   INTEGER :: Job
!
! Local variable declarations rewritten by SPAG
!
   INTEGER :: j , n1 , n2
   EXTERNAL dperm
!
! End of declarations rewritten by SPAG
!
!-----------------------------------------------------------------------
! This routine performs a symmetric permutation of the rows and
! columns of a matrix stored in MSR format. i.e., it computes
! B = P A transp(P), where P, is  a permutation matrix.
! P maps row i into row perm(i) and column j into column perm(j):
!      a(i,j)    becomes   a(perm(i),perm(j)) in new matrix
! (i.e.  ao(perm(i),perm(j)) = a(i,j) )
! calls dperm.
!-----------------------------------------------------------------------
! Y. Saad, Nov 15, 1991.
!-----------------------------------------------------------------------
! on entry:
!----------
! n 	= dimension of the matrix
! a, ja = input matrix in MSR format.
! perm 	= integer array of length n containing the permutation arrays
!	  for the rows: perm(i) is the destination of row i in the
!         permuted matrix -- also the destination of column i in case
!         permutation is symmetric (job .le. 2)
!
! job	= integer indicating the work to be done:
! 		job = 1	permute a, ja, ia into ao, jao, iao
! 		job = 2 permute matrix ignoring real values.
!
! on return:
!-----------
! ao, jao = output matrix in MSR.
!
! in case job .eq. 2 a and ao are never referred to and can be dummy
! arguments.
!
! Notes:
!-------
!  1) algorithm is NOT in place
!  2) column indices may not be sorted on return even  though they may be
!     on entry.
!----------------------------------------------------------------------c
!     local variables
!
   n1 = Nrow + 1
   n2 = n1 + 1
!
   CALL dperm(Nrow,A,Ja,Ja,Ao(n2),Jao(n2),Jao,Perm,Perm,Job)
!
   Jao(1) = n2
   DO j = 1 , Nrow
      Ao(Perm(j)) = A(j)
      Jao(j+1) = Jao(j+1) + n1
   ENDDO
!
! done
!
!-----------------------------------------------------------------------
!--------end-of-dmperm--------------------------------------------------
END SUBROUTINE dmperm
!*==dvperm.f90 processed by SPAG 8.04RA 15:47 21 Feb 2024
!!SPAG Open source Personal, Educational or Academic User Clemson University  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
!-----------------------------------------------------------------------
SUBROUTINE dvperm(N,X,Perm)
   USE ISO_FORTRAN_ENV                 
   IMPLICIT NONE
!
! Dummy argument declarations rewritten by SPAG
!
   INTEGER , INTENT(IN) :: N
   REAL(REAL64) , INTENT(INOUT) , DIMENSION(N) :: X
   INTEGER , INTENT(INOUT) , DIMENSION(N) :: Perm
!
! Local variable declarations rewritten by SPAG
!
   INTEGER :: ii , init , j , k , next
   REAL(REAL64) :: tmp , tmp1
!
! End of declarations rewritten by SPAG
!
!-----------------------------------------------------------------------
! this subroutine performs an in-place permutation of a real vector x
! according to the permutation array perm(*), i.e., on return,
! the vector x satisfies,
!
!	x(perm(j)) :== x(j), j=1,2,.., n
!
!-----------------------------------------------------------------------
! on entry:
!---------
! n 	= length of vector x.
! perm 	= integer array of length n containing the permutation  array.
! x	= input vector
!
! on return:
!----------
! x	= vector x permuted according to x(perm(*)) :=  x(*)
!
!----------------------------------------------------------------------c
!           Y. Saad, Sep. 21 1989                                      c
!----------------------------------------------------------------------c
! local variables
   INTEGER :: spag_nextblock_1
   spag_nextblock_1 = 1
   SPAG_DispatchLoop_1: DO
      SELECT CASE (spag_nextblock_1)
      CASE (1)
!
         init = 1
         tmp = X(init)
         ii = Perm(init)
         Perm(init) = -Perm(init)
         k = 0
         DO
!
! loop
!
            k = k + 1
!
! save the chased element --
!
            tmp1 = X(ii)
            X(ii) = tmp
            next = Perm(ii)
            IF ( next<0 ) THEN
               SPAG_Loop_3_1: DO
!
! reinitilaize cycle --
!
                  init = init + 1
                  IF ( init>N ) THEN
                     spag_nextblock_1 = 2
                     CYCLE SPAG_DispatchLoop_1
                  ENDIF
                  IF ( Perm(init)>=0 ) THEN
                     tmp = X(init)
                     ii = Perm(init)
                     Perm(init) = -Perm(init)
                     EXIT SPAG_Loop_3_1
                  ENDIF
               ENDDO SPAG_Loop_3_1
            ELSE
!
! test for end
!
               IF ( k>N ) THEN
                  spag_nextblock_1 = 2
                  CYCLE SPAG_DispatchLoop_1
               ENDIF
               tmp = tmp1
               Perm(ii) = -Perm(ii)
!
! end loop
!
               ii = next
            ENDIF
         ENDDO
         spag_nextblock_1 = 2
      CASE (2)
!
         DO j = 1 , N
            Perm(j) = -Perm(j)
         ENDDO
         EXIT SPAG_DispatchLoop_1
      END SELECT
   ENDDO SPAG_DispatchLoop_1
!
!-------------------end-of-dvperm---------------------------------------
!-----------------------------------------------------------------------
END SUBROUTINE dvperm
!*==ivperm.f90 processed by SPAG 8.04RA 15:47 21 Feb 2024
!!SPAG Open source Personal, Educational or Academic User Clemson University  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
!-----------------------------------------------------------------------
SUBROUTINE ivperm(N,Ix,Perm)
   IMPLICIT NONE
!
! Dummy argument declarations rewritten by SPAG
!
   INTEGER , INTENT(IN) :: N
   INTEGER , INTENT(INOUT) , DIMENSION(N) :: Ix
   INTEGER , INTENT(INOUT) , DIMENSION(N) :: Perm
!
! Local variable declarations rewritten by SPAG
!
   INTEGER :: ii , init , j , k , next , tmp , tmp1
!
! End of declarations rewritten by SPAG
!
!-----------------------------------------------------------------------
! this subroutine performs an in-place permutation of an integer vector
! ix according to the permutation array perm(*), i.e., on return,
! the vector x satisfies,
!
!	ix(perm(j)) :== ix(j), j=1,2,.., n
!
!-----------------------------------------------------------------------
! on entry:
!---------
! n 	= length of vector x.
! perm 	= integer array of length n containing the permutation  array.
! ix	= input vector
!
! on return:
!----------
! ix	= vector x permuted according to ix(perm(*)) :=  ix(*)
!
!----------------------------------------------------------------------c
!           Y. Saad, Sep. 21 1989                                      c
!----------------------------------------------------------------------c
! local variables
   INTEGER :: spag_nextblock_1
   spag_nextblock_1 = 1
   SPAG_DispatchLoop_1: DO
      SELECT CASE (spag_nextblock_1)
      CASE (1)
!
         init = 1
         tmp = Ix(init)
         ii = Perm(init)
         Perm(init) = -Perm(init)
         k = 0
         DO
!
! loop
!
            k = k + 1
!
! save the chased element --
!
            tmp1 = Ix(ii)
            Ix(ii) = tmp
            next = Perm(ii)
            IF ( next<0 ) THEN
               SPAG_Loop_3_1: DO
!
! reinitilaize cycle --
!
                  init = init + 1
                  IF ( init>N ) THEN
                     spag_nextblock_1 = 2
                     CYCLE SPAG_DispatchLoop_1
                  ENDIF
                  IF ( Perm(init)>=0 ) THEN
                     tmp = Ix(init)
                     ii = Perm(init)
                     Perm(init) = -Perm(init)
                     EXIT SPAG_Loop_3_1
                  ENDIF
               ENDDO SPAG_Loop_3_1
            ELSE
!
! test for end
!
               IF ( k>N ) THEN
                  spag_nextblock_1 = 2
                  CYCLE SPAG_DispatchLoop_1
               ENDIF
               tmp = tmp1
               Perm(ii) = -Perm(ii)
!
! end loop
!
               ii = next
            ENDIF
         ENDDO
         spag_nextblock_1 = 2
      CASE (2)
!
         DO j = 1 , N
            Perm(j) = -Perm(j)
         ENDDO
         EXIT SPAG_DispatchLoop_1
      END SELECT
   ENDDO SPAG_DispatchLoop_1
!
!-------------------end-of-ivperm---------------------------------------
!-----------------------------------------------------------------------
END SUBROUTINE ivperm
!*==retmx.f90 processed by SPAG 8.04RA 15:47 21 Feb 2024
!!SPAG Open source Personal, Educational or Academic User Clemson University  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
!-----------------------------------------------------------------------
SUBROUTINE retmx(N,A,Ja,Ia,Dd)
   USE ISO_FORTRAN_ENV                 
   IMPLICIT NONE
!
! Dummy argument declarations rewritten by SPAG
!
   INTEGER , INTENT(IN) :: N
   REAL(REAL64) , INTENT(IN) , DIMENSION(*) :: A
   INTEGER , INTENT(IN) , DIMENSION(*) :: Ja
   INTEGER , INTENT(IN) , DIMENSION(*) :: Ia
   REAL(REAL64) , INTENT(OUT) , DIMENSION(*) :: Dd
!
! Local variable declarations rewritten by SPAG
!
   INTEGER :: i , k , k1 , k2
   REAL(REAL64) :: t , t1 , t2
!
! End of declarations rewritten by SPAG
!
!-----------------------------------------------------------------------
! returns in dd(*) the max absolute value of elements in row *.
! used for scaling purposes. superseded by rnrms  .
!
! on entry:
! n	= dimension of A
! a,ja,ia
!	= matrix stored in compressed sparse row format
! dd	= real*8 array of length n. On output,entry dd(i) contains
!	  the element of row i that has the largest absolute value.
!	  Moreover the sign of dd is modified such that it is the
!	  same as that of the diagonal element in row i.
!----------------------------------------------------------------------c
!           Y. Saad, Sep. 21 1989                                      c
!----------------------------------------------------------------------c
! local variables
   t = 0.0
   t2 = 0.0
   t1 = 0.0
!
! initialize
!
   k2 = 1
   DO i = 1 , N
      k1 = k2
      k2 = Ia(i+1) - 1
      t = 0.0D0
      DO k = k1 , k2
         t1 = abs(A(k))
         IF ( t1>t ) t = t1
         IF ( Ja(k)==i ) THEN
            IF ( A(k)>=0.0 ) THEN
               t2 = A(k)
            ELSE
               t2 = -A(k)
            ENDIF
         ENDIF
      ENDDO
      Dd(i) = t2*t
!     we do not invert diag
   ENDDO
!---------end of retmx -------------------------------------------------
!-----------------------------------------------------------------------
END SUBROUTINE retmx
!*==diapos.f90 processed by SPAG 8.04RA 15:47 21 Feb 2024
!!SPAG Open source Personal, Educational or Academic User Clemson University  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
!-----------------------------------------------------------------------
SUBROUTINE diapos(N,Ja,Ia,Idiag)
   IMPLICIT NONE
!
! Dummy argument declarations rewritten by SPAG
!
   INTEGER , INTENT(IN) :: N
   INTEGER , INTENT(IN) , DIMENSION(*) :: Ja
   INTEGER , INTENT(IN) , DIMENSION(N+1) :: Ia
   INTEGER , INTENT(OUT) , DIMENSION(N) :: Idiag
!
! Local variable declarations rewritten by SPAG
!
   INTEGER :: i , k
!
! End of declarations rewritten by SPAG
!
!-----------------------------------------------------------------------
! this subroutine returns the positions of the diagonal elements of a
! sparse matrix a, ja, ia, in the array idiag.
!-----------------------------------------------------------------------
! on entry:
!----------
!
! n	= integer. row dimension of the matrix a.
! a,ja,
!    ia = matrix stored compressed sparse row format. a array skipped.
!
! on return:
!-----------
! idiag  = integer array of length n. The i-th entry of idiag
!          points to the diagonal element a(i,i) in the arrays
!          a, ja. (i.e., a(idiag(i)) = element A(i,i) of matrix A)
!          if no diagonal element is found the entry is set to 0.
!----------------------------------------------------------------------c
!           Y. Saad, March, 1990
!----------------------------------------------------------------------c
   DO i = 1 , N
      Idiag(i) = 0
   ENDDO
!
!     sweep through data structure.
!
   DO i = 1 , N
      DO k = Ia(i) , Ia(i+1) - 1
         IF ( Ja(k)==i ) Idiag(i) = k
      ENDDO
   ENDDO
!----------- -end-of-diapos---------------------------------------------
!-----------------------------------------------------------------------
END SUBROUTINE diapos
!*==dscaldg.f90 processed by SPAG 8.04RA 15:47 21 Feb 2024
!!SPAG Open source Personal, Educational or Academic User Clemson University  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
!-----------------------------------------------------------------------
SUBROUTINE dscaldg(N,A,Ja,Ia,Diag,Job)
   USE ISO_FORTRAN_ENV                 
   IMPLICIT NONE
!
! Dummy argument declarations rewritten by SPAG
!
   INTEGER :: N
   REAL(REAL64) , INTENT(INOUT) , DIMENSION(*) :: A
   INTEGER , DIMENSION(*) :: Ja
   INTEGER , DIMENSION(*) :: Ia
   REAL(REAL64) , INTENT(INOUT) , DIMENSION(*) :: Diag
   INTEGER , INTENT(IN) :: Job
!
! Local variable declarations rewritten by SPAG
!
   INTEGER :: i , j , k , k1 , k2
   REAL(REAL64) :: t
   EXTERNAL retmx
!
! End of declarations rewritten by SPAG
!
!-----------------------------------------------------------------------
! scales rows by diag where diag is either given (job=0)
! or to be computed:
!  job = 1 ,scale row i by  by  +/- max |a(i,j) | and put inverse of
!       scaling factor in diag(i),where +/- is the sign of a(i,i).
!  job = 2 scale by 2-norm of each row..
! if diag(i) = 0,then diag(i) is replaced by one
! (no scaling)..
!----------------------------------------------------------------------c
!           Y. Saad, Sep. 21 1989                                      c
!----------------------------------------------------------------------c
   IF ( Job+1==1 ) THEN
   ELSEIF ( Job+1==2 ) THEN
      CALL retmx(N,A,Ja,Ia,Diag)
   ELSE
      DO j = 1 , N
         k1 = Ia(j)
         k2 = Ia(j+1) - 1
         t = 0.0D0
         DO k = k1 , k2
            t = t + A(k)*A(k)
         ENDDO
         Diag(j) = sqrt(t)
      ENDDO
   ENDIF
!------
   DO j = 1 , N
      IF ( Diag(j)/=0.0D0 ) THEN
         Diag(j) = 1.0D0/Diag(j)
      ELSE
         Diag(j) = 1.0D0
      ENDIF
   ENDDO
   DO i = 1 , N
      t = Diag(i)
      DO k = Ia(i) , Ia(i+1) - 1
         A(k) = A(k)*t
      ENDDO
   ENDDO
!--------end of dscaldg -----------------------------------------------
!-----------------------------------------------------------------------
END SUBROUTINE dscaldg
!*==extbdg.f90 processed by SPAG 8.04RA 15:47 21 Feb 2024
!!SPAG Open source Personal, Educational or Academic User Clemson University  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
!-----------------------------------------------------------------------
SUBROUTINE extbdg(N,A,Ja,Ia,Bdiag,Nblk,Ao,Jao,Iao)
   USE ISO_FORTRAN_ENV                 
   IMPLICIT NONE
!
! Dummy argument declarations rewritten by SPAG
!
   INTEGER , INTENT(IN) :: N
   REAL(REAL64) , INTENT(IN) , DIMENSION(*) :: A
   INTEGER , INTENT(IN) , DIMENSION(*) :: Ja
   INTEGER , INTENT(IN) , DIMENSION(*) :: Ia
   REAL(REAL64) , INTENT(OUT) , DIMENSION(*) :: Bdiag
   INTEGER , INTENT(IN) :: Nblk
   REAL(REAL64) , INTENT(OUT) , DIMENSION(*) :: Ao
   INTEGER , INTENT(OUT) , DIMENSION(*) :: Jao
   INTEGER , INTENT(OUT) , DIMENSION(*) :: Iao
!
! Local variable declarations rewritten by SPAG
!
   INTEGER :: i , j , j1 , j2 , jj , k , kb , ko , l , ltr , m
!
! End of declarations rewritten by SPAG
!
!-----------------------------------------------------------------------
! this subroutine extracts the main diagonal blocks of a
! matrix stored in compressed sparse row format and puts the result
! into the array bdiag and the remainder in ao,jao,iao.
!-----------------------------------------------------------------------
! on entry:
!----------
! n	= integer. The row dimension of the matrix a.
! a,
! ja,
! ia    = matrix stored in csr format
! nblk  = dimension of each diagonal block. The diagonal blocks are
!         stored in compressed format rowwise,i.e.,we store in
!	  succession the i nonzeros of the i-th row after those of
!	  row number i-1..
!
! on return:
!----------
! bdiag = real*8 array of size (n x nblk) containing the diagonal
!	  blocks of A on return
! ao,
! jao,
! iao   = remainder of the matrix stored in csr format.
!----------------------------------------------------------------------c
!           Y. Saad, Sep. 21 1989                                      c
!----------------------------------------------------------------------c
   m = 1 + (N-1)/Nblk
! this version is sequential -- there is a more parallel version
! that goes through the structure twice ....
   ltr = ((Nblk-1)*Nblk)/2
   l = m*ltr
   DO i = 1 , l
      Bdiag(i) = 0.0D0
   ENDDO
   ko = 0
   kb = 1
   Iao(1) = 1
!-------------------------
   DO jj = 1 , m
      j1 = (jj-1)*Nblk + 1
      j2 = min0(N,j1+Nblk-1)
      DO j = j1 , j2
         DO i = Ia(j) , Ia(j+1) - 1
            k = Ja(i)
            IF ( k<j1 ) THEN
               ko = ko + 1
               Ao(ko) = A(i)
               Jao(ko) = k
            ELSEIF ( k<j ) THEN
!     kb = (jj-1)*ltr+((j-j1)*(j-j1-1))/2+k-j1+1
!     bdiag(kb) = a(i)
               Bdiag(kb+k-j1) = A(i)
            ENDIF
         ENDDO
         kb = kb + j - j1
         Iao(j+1) = ko + 1
      ENDDO
   ENDDO
!---------end-of-extbdg-------------------------------------------------
!-----------------------------------------------------------------------
END SUBROUTINE extbdg
!*==getbwd.f90 processed by SPAG 8.04RA 15:47 21 Feb 2024
!!SPAG Open source Personal, Educational or Academic User Clemson University  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
!-----------------------------------------------------------------------
SUBROUTINE getbwd(N,Ja,Ia,Ml,Mu)
   IMPLICIT NONE
!
! Dummy argument declarations rewritten by SPAG
!
   INTEGER , INTENT(IN) :: N
   INTEGER , INTENT(IN) , DIMENSION(*) :: Ja
   INTEGER , INTENT(IN) , DIMENSION(N+1) :: Ia
   INTEGER , INTENT(INOUT) :: Ml
   INTEGER , INTENT(INOUT) :: Mu
!
! Local variable declarations rewritten by SPAG
!
   INTEGER :: i , k , ldist
!
! End of declarations rewritten by SPAG
!
!-----------------------------------------------------------------------
! gets the bandwidth of lower part and upper part of A.
! does not assume that A is sorted.
!-----------------------------------------------------------------------
! on entry:
!----------
! n	= integer = the row dimension of the matrix
! a, ja,
!    ia = matrix in compressed sparse row format.
!        Only arrays ja, ia are used and passed
! on return:
!-----------
! ml	= integer. The bandwidth of the strict lower part of A
! mu	= integer. The bandwidth of the strict upper part of A
!
! Notes:
! ===== ml and mu are allowed to be negative or return. This may be
!       useful since it will tell us whether a band is confined
!       in the strict  upper/lower triangular part.
!       indeed the definitions of ml and mu are
!
!       ml = max ( (i-j)  s.t. a(i,j) .ne. 0  )
!       mu = max ( (j-i)  s.t. a(i,j) .ne. 0  )
!----------------------------------------------------------------------c
! Y. Saad, Sep. 21 1989                                                c
!----------------------------------------------------------------------c
   Ml = -N
   Mu = -N
   DO i = 1 , N
      DO k = Ia(i) , Ia(i+1) - 1
         ldist = i - Ja(k)
         Ml = max(Ml,ldist)
         Mu = max(Mu,-ldist)
      ENDDO
   ENDDO
!---------------end-of-getbwd ------------------------------------------
!-----------------------------------------------------------------------
END SUBROUTINE getbwd
!*==blkfnd.f90 processed by SPAG 8.04RA 15:47 21 Feb 2024
!!SPAG Open source Personal, Educational or Academic User Clemson University  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
!-----------------------------------------------------------------------
SUBROUTINE blkfnd(Nrow,Ja,Ia,Nblk)
   IMPLICIT NONE
!
! Dummy argument declarations rewritten by SPAG
!
   INTEGER :: Nrow
   INTEGER , DIMENSION(*) :: Ja
   INTEGER , DIMENSION(Nrow+1) :: Ia
   INTEGER , INTENT(INOUT) :: Nblk
!
! Local variable declarations rewritten by SPAG
!
   INTEGER :: i , i1 , i2 , iblk , imsg , irow , jf , jfirst , jl , jlast , jrow , len , len0 , minlen
   EXTERNAL blkchk
!
! End of declarations rewritten by SPAG
!
!-----------------------------------------------------------------------
! This routine attemptps to determine whether or not  the input
! matrix has a block structure and finds the blocks size
! if it does. A block matrix is one which is
! comprised of small square dense blocks. If there are zero
! elements within the square blocks and the original data structure
! takes these zeros into account then blkchk may fail to find the
! correct block size.
!-----------------------------------------------------------------------
! on entry
!---------
! nrow	= integer equal to the row dimension of the matrix.
! ja    = integer array containing the column indices of the entries
!         nonzero entries of the matrix stored by row.
! ia    = integer array of length nrow + 1 containing the pointers
!         beginning of each row in array ja.
!
! nblk  = integer containing the assumed value of nblk if job = 0
!
! on return
!----------
! nblk  = integer containing the value found for nblk when job = 1.
!         if imsg .ne. 0 this value is meaningless however.
!
!----------------------------------------------------------------------c
!           Y. Saad, Sep. 21 1989                                      c
!----------------------------------------------------------------------c
!-----------------------------------------------------------------------
! first part of code will find candidate block sizes.
! criterion used here is a simple one: scan rows and  determine groups
! of rows that have the same length and such that the first column
! number and the last column number are identical.
!-----------------------------------------------------------------------
   minlen = Ia(2) - Ia(1)
   irow = 1
   DO i = 2 , Nrow
      len = Ia(i+1) - Ia(i)
      IF ( len<minlen ) THEN
         minlen = len
         irow = i
      ENDIF
   ENDDO
!
!     ---- candidates are all dividers of minlen
!
   Nblk = 1
   IF ( minlen<=1 ) RETURN
!
   SPAG_Loop_1_1: DO iblk = minlen , 1 , -1
      IF ( mod(minlen,iblk)==0 ) THEN
         len = Ia(2) - Ia(1)
         len0 = len
         jfirst = Ja(1)
         jlast = Ja(Ia(2)-1)
         DO jrow = irow + 1 , irow + Nblk - 1
            i1 = Ia(jrow)
            i2 = Ia(jrow+1) - 1
            len = i2 + 1 - i1
            jf = Ja(i1)
            jl = Ja(i2)
            IF ( len/=len0 .OR. jf/=jfirst .OR. jl/=jlast ) CYCLE SPAG_Loop_1_1
         ENDDO
!
!     check for this candidate ----
!
         CALL blkchk(Nrow,Ja,Ia,iblk,imsg)
         IF ( imsg==0 ) THEN
!
!     block size found
!
            Nblk = iblk
            RETURN
         ENDIF
      ENDIF
   ENDDO SPAG_Loop_1_1
!--------end-of-blkfnd -------------------------------------------------
!-----------------------------------------------------------------------
END SUBROUTINE blkfnd
!*==blkchk.f90 processed by SPAG 8.04RA 15:47 21 Feb 2024
!!SPAG Open source Personal, Educational or Academic User Clemson University  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
!-----------------------------------------------------------------------
SUBROUTINE blkchk(Nrow,Ja,Ia,Nblk,Imsg)
   IMPLICIT NONE
!
! Dummy argument declarations rewritten by SPAG
!
   INTEGER , INTENT(IN) :: Nrow
   INTEGER , INTENT(IN) , DIMENSION(*) :: Ja
   INTEGER , INTENT(IN) , DIMENSION(Nrow+1) :: Ia
   INTEGER , INTENT(IN) :: Nblk
   INTEGER , INTENT(OUT) :: Imsg
!
! Local variable declarations rewritten by SPAG
!
   INTEGER :: i , i1 , ii , irow , j , j2 , jstart , k , len , lena , nr
!
! End of declarations rewritten by SPAG
!
!-----------------------------------------------------------------------
! This routine checks whether the input matrix is a block
! matrix with block size of nblk. A block matrix is one which is
! comprised of small square dense blocks. If there are zero
! elements within the square blocks and the data structure
! takes them into account then blkchk may fail to find the
! correct block size.
!-----------------------------------------------------------------------
! on entry
!---------
! nrow	= integer equal to the row dimension of the matrix.
! ja    = integer array containing the column indices of the entries
!         nonzero entries of the matrix stored by row.
! ia    = integer array of length nrow + 1 containing the pointers
!         beginning of each row in array ja.
!
! nblk  = integer containing the value of nblk to be checked.
!
! on return
!----------
!
! imsg  = integer containing a message  with the following meaning.
!          imsg = 0 means that the output value of nblk is a correct
!                   block size. nblk .lt. 0 means nblk not correct
!                   block size.
!          imsg = -1 : nblk does not divide nrow
!          imsg = -2 : a starting element in a row is at wrong position
!             (j .ne. mult*nblk +1 )
!          imsg = -3 : nblk does divide a row length -
!          imsg = -4 : an element is isolated outside a block or
!             two rows in same group have different lengths
!----------------------------------------------------------------------c
!           Y. Saad, Sep. 21 1989                                      c
!----------------------------------------------------------------------c
!----------------------------------------------------------------------
! first part of code will find candidate block sizes.
! this is not guaranteed to work . so a check is done at the end
! the criterion used here is a simple one:
! scan rows and determine groups of rows that have the same length
! and such that the first column number and the last column number
! are identical.
!----------------------------------------------------------------------
   Imsg = 0
   IF ( Nblk<=1 ) RETURN
   nr = Nrow/Nblk
   IF ( nr*Nblk/=Nrow ) THEN
      Imsg = -1
      RETURN
   ELSE
!--   main loop ---------------------------------------------------------
      irow = 1
      DO ii = 1 , nr
!     i1= starting position for group of nblk rows in original matrix
         i1 = Ia(irow)
         j2 = i1
!     lena = length of each row in that group  in the original matrix
         lena = Ia(irow+1) - i1
!     len = length of each block-row in that group in the output matrix
         len = lena/Nblk
         IF ( len*Nblk/=lena ) THEN
            CALL spag_block_2
            RETURN
         ENDIF
!
!     for each row
!
         DO i = 1 , Nblk
            irow = irow + 1
            IF ( Ia(irow)-Ia(irow-1)/=lena ) THEN
               CALL spag_block_3
               RETURN
            ENDIF
!
!     for each block
!
            DO k = 0 , len - 1
               jstart = Ja(i1+Nblk*k) - 1
               IF ( (jstart/Nblk)*Nblk/=jstart ) THEN
                  CALL spag_block_1
                  RETURN
               ENDIF
!
!     for each column
!
               DO j = 1 , Nblk
                  IF ( jstart+j/=Ja(j2) ) THEN
                     CALL spag_block_3
                     RETURN
                  ENDIF
                  j2 = j2 + 1
               ENDDO
            ENDDO
         ENDDO
      ENDDO
!     went through all loops successfully:
      RETURN
   ENDIF
CONTAINS
   SUBROUTINE spag_block_1
      Imsg = -2
      RETURN
   END SUBROUTINE spag_block_1
   SUBROUTINE spag_block_2
      Imsg = -3
      RETURN
   END SUBROUTINE spag_block_2
   SUBROUTINE spag_block_3
      Imsg = -4
   END SUBROUTINE spag_block_3
!----------------end of chkblk -----------------------------------------
!-----------------------------------------------------------------------
END SUBROUTINE blkchk
!*==infdia.f90 processed by SPAG 8.04RA 15:47 21 Feb 2024
!!SPAG Open source Personal, Educational or Academic User Clemson University  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
!-----------------------------------------------------------------------
SUBROUTINE infdia(N,Ja,Ia,Ind,Idiag)
   IMPLICIT NONE
!
! Dummy argument declarations rewritten by SPAG
!
   INTEGER , INTENT(IN) :: N
   INTEGER , INTENT(IN) , DIMENSION(*) :: Ja
   INTEGER , INTENT(IN) , DIMENSION(*) :: Ia
   INTEGER , INTENT(INOUT) , DIMENSION(*) :: Ind
   INTEGER , INTENT(INOUT) :: Idiag
!
! Local variable declarations rewritten by SPAG
!
   INTEGER :: i , j , k , n2
!
! End of declarations rewritten by SPAG
!
!-----------------------------------------------------------------------
!     obtains information on the diagonals of A.
!-----------------------------------------------------------------------
! this subroutine finds the lengths of each of the 2*n-1 diagonals of A
! it also outputs the number of nonzero diagonals found.
!-----------------------------------------------------------------------
! on entry:
!----------
! n	= dimension of the matrix a.
!
! a,    ..... not needed here.
! ja,
! ia    = matrix stored in csr format
!
! on return:
!-----------
!
! idiag = integer. number of nonzero diagonals found.
!
! ind   = integer array of length at least 2*n-1. The k-th entry in
!         ind contains the number of nonzero elements in the diagonal
!         number k, the numbering beeing from the lowermost diagonal
!         (bottom-left). In other words ind(k) = length of diagonal
!         whose offset wrt the main diagonal is = - n + k.
!----------------------------------------------------------------------c
!           Y. Saad, Sep. 21 1989                                      c
!----------------------------------------------------------------------c
   n2 = N + N - 1
   DO i = 1 , n2
      Ind(i) = 0
   ENDDO
   DO i = 1 , N
      DO k = Ia(i) , Ia(i+1) - 1
         j = Ja(k)
         Ind(N+j-i) = Ind(N+j-i) + 1
      ENDDO
   ENDDO
!     count the nonzero ones.
   Idiag = 0
   DO k = 1 , n2
      IF ( Ind(k)/=0 ) Idiag = Idiag + 1
   ENDDO
! done
!------end-of-infdia ---------------------------------------------------
!-----------------------------------------------------------------------
END SUBROUTINE infdia
!*==amubdg.f90 processed by SPAG 8.04RA 15:47 21 Feb 2024
!!SPAG Open source Personal, Educational or Academic User Clemson University  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
!-----------------------------------------------------------------------
SUBROUTINE amubdg(Nrow,Ncol,Ncolb,Ja,Ia,Jb,Ib,Ndegr,Nnz,Iw)
   IMPLICIT NONE
!
! Dummy argument declarations rewritten by SPAG
!
   INTEGER , INTENT(IN) :: Nrow
   INTEGER , INTENT(IN) :: Ncol
   INTEGER , INTENT(IN) :: Ncolb
   INTEGER , INTENT(IN) , DIMENSION(*) :: Ja
   INTEGER , INTENT(IN) , DIMENSION(Nrow+1) :: Ia
   INTEGER , INTENT(IN) , DIMENSION(*) :: Jb
   INTEGER , INTENT(IN) , DIMENSION(Ncol+1) :: Ib
   INTEGER , INTENT(INOUT) , DIMENSION(Nrow) :: Ndegr
   INTEGER , INTENT(INOUT) :: Nnz
   INTEGER , INTENT(INOUT) , DIMENSION(Ncolb) :: Iw
!
! Local variable declarations rewritten by SPAG
!
   INTEGER :: ii , j , jc , jr , k , last , ldg
!
! End of declarations rewritten by SPAG
!
!-----------------------------------------------------------------------
! gets the number of nonzero elements in each row of A*B and the total
! number of nonzero elements in A*B.
!-----------------------------------------------------------------------
! on entry:
! --------
!
! nrow  = integer.  row dimension of matrix A
! ncol  = integer.  column dimension of matrix A = row dimension of
!                   matrix B.
! ncolb = integer. the colum dimension of the matrix B.
!
! ja, ia= row structure of input matrix A: ja = column indices of
!         the nonzero elements of A stored by rows.
!         ia = pointer to beginning of each row  in ja.
!
! jb, ib= row structure of input matrix B: jb = column indices of
!         the nonzero elements of A stored by rows.
!         ib = pointer to beginning of each row  in jb.
!
! on return:
! ---------
! ndegr	= integer array of length nrow containing the degrees (i.e.,
!         the number of nonzeros in  each row of the matrix A * B
!
! nnz   = total number of nonzero elements found in A * B
!
! work arrays:
!-------------
! iw	= integer work array of length ncolb.
!-----------------------------------------------------------------------
   DO k = 1 , Ncolb
      Iw(k) = 0
   ENDDO
 
   DO k = 1 , Nrow
      Ndegr(k) = 0
   ENDDO
!
!     method used: Transp(A) * A = sum [over i=1, nrow]  a(i)^T a(i)
!     where a(i) = i-th row of  A. We must be careful not to add  the
!     elements already accounted for.
!
!
   DO ii = 1 , Nrow
!
!     for each row of A
!
      ldg = 0
!
!    end-of-linked list
!
      last = -1
      DO j = Ia(ii) , Ia(ii+1) - 1
!
!     row number to be added:
!
         jr = Ja(j)
         DO k = Ib(jr) , Ib(jr+1) - 1
            jc = Jb(k)
            IF ( Iw(jc)==0 ) THEN
!
!     add one element to the linked list
!
               ldg = ldg + 1
               Iw(jc) = last
               last = jc
            ENDIF
         ENDDO
      ENDDO
      Ndegr(ii) = ldg
!
!     reset iw to zero
!
      DO k = 1 , ldg
         j = Iw(last)
         Iw(last) = 0
         last = j
      ENDDO
!-----------------------------------------------------------------------
   ENDDO
!
   Nnz = 0
   DO ii = 1 , Nrow
      Nnz = Nnz + Ndegr(ii)
   ENDDO
!
!---------------end-of-amubdg ------------------------------------------
!-----------------------------------------------------------------------
END SUBROUTINE amubdg
!*==aplbdg.f90 processed by SPAG 8.04RA 15:47 21 Feb 2024
!!SPAG Open source Personal, Educational or Academic User Clemson University  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
!-----------------------------------------------------------------------
SUBROUTINE aplbdg(Nrow,Ncol,Ja,Ia,Jb,Ib,Ndegr,Nnz,Iw)
   IMPLICIT NONE
!
! Dummy argument declarations rewritten by SPAG
!
   INTEGER , INTENT(IN) :: Nrow
   INTEGER , INTENT(IN) :: Ncol
   INTEGER , INTENT(IN) , DIMENSION(*) :: Ja
   INTEGER , INTENT(IN) , DIMENSION(Nrow+1) :: Ia
   INTEGER , INTENT(IN) , DIMENSION(*) :: Jb
   INTEGER , INTENT(IN) , DIMENSION(Nrow+1) :: Ib
   INTEGER , INTENT(INOUT) , DIMENSION(Nrow) :: Ndegr
   INTEGER , INTENT(INOUT) :: Nnz
   INTEGER , INTENT(INOUT) , DIMENSION(Ncol) :: Iw
!
! Local variable declarations rewritten by SPAG
!
   INTEGER :: ii , j , jc , jr , k , last , ldg
!
! End of declarations rewritten by SPAG
!
!-----------------------------------------------------------------------
! gets the number of nonzero elements in each row of A+B and the total
! number of nonzero elements in A+B.
!-----------------------------------------------------------------------
! on entry:
! ---------
! nrow	= integer. The row dimension of A and B
! ncol  = integer. The column dimension of A and B.
!
! a,
! ja,
! ia   = Matrix A in compressed sparse row format.
!
! b,
! jb,
! ib	=  Matrix B in compressed sparse row format.
!
! on return:
!----------
! ndegr	= integer array of length nrow containing the degrees (i.e.,
!         the number of nonzeros in  each row of the matrix A + B.
!
! nnz   = total number of nonzero elements found in A * B
!
! work arrays:
!------------
! iw	= integer work array of length equal to ncol.
!
!-----------------------------------------------------------------------
   DO k = 1 , Ncol
      Iw(k) = 0
   ENDDO
!
   DO k = 1 , Nrow
      Ndegr(k) = 0
   ENDDO
!
   DO ii = 1 , Nrow
      ldg = 0
!
!    end-of-linked list
!
      last = -1
!
!     row of A
!
      DO j = Ia(ii) , Ia(ii+1) - 1
         jr = Ja(j)
!
!     add element to the linked list
!
         ldg = ldg + 1
         Iw(jr) = last
         last = jr
      ENDDO
!
!     row of B
!
      DO j = Ib(ii) , Ib(ii+1) - 1
         jc = Jb(j)
         IF ( Iw(jc)==0 ) THEN
!
!     add one element to the linked list
!
            ldg = ldg + 1
            Iw(jc) = last
            last = jc
         ENDIF
      ENDDO
!     done with row ii.
      Ndegr(ii) = ldg
!
!     reset iw to zero
!
      DO k = 1 , ldg
         j = Iw(last)
         Iw(last) = 0
         last = j
      ENDDO
!-----------------------------------------------------------------------
   ENDDO
!
   Nnz = 0
   DO ii = 1 , Nrow
      Nnz = Nnz + Ndegr(ii)
   ENDDO
!----------------end-of-aplbdg -----------------------------------------
!-----------------------------------------------------------------------
END SUBROUTINE aplbdg
!*==rnrms.f90 processed by SPAG 8.04RA 15:47 21 Feb 2024
!!SPAG Open source Personal, Educational or Academic User Clemson University  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
!-----------------------------------------------------------------------
SUBROUTINE rnrms(Nrow,Nrm,A,Ia,Diag)
   USE ISO_FORTRAN_ENV                 
   IMPLICIT NONE
!
! Dummy argument declarations rewritten by SPAG
!
   INTEGER , INTENT(IN) :: Nrow
   INTEGER , INTENT(IN) :: Nrm
   REAL(REAL64) , INTENT(IN) , DIMENSION(*) :: A
   INTEGER , INTENT(IN) , DIMENSION(Nrow+1) :: Ia
   REAL(REAL64) , INTENT(OUT) , DIMENSION(Nrow) :: Diag
!
! Local variable declarations rewritten by SPAG
!
   INTEGER :: ii , k , k1 , k2
   REAL(REAL64) :: scal
!
! End of declarations rewritten by SPAG
!
!-----------------------------------------------------------------------
! gets the norms of each row of A. (choice of three norms)
!-----------------------------------------------------------------------
! on entry:
! ---------
! nrow	= integer. The row dimension of A
!
! nrm   = integer. norm indicator. nrm = 1, means 1-norm, nrm =2
!                  means the 2-nrm, nrm = 0 means max norm
!
! a,
! ja,
! ia   = Matrix A in compressed sparse row format.
!        ja is not used - not passed.
! on return:
!----------
!
! diag = real vector of length nrow containing the norms
!
!-----------------------------------------------------------------
   DO ii = 1 , Nrow
!
!     compute the norm if each element.
!
      scal = 0.0D0
      k1 = Ia(ii)
      k2 = Ia(ii+1) - 1
      IF ( Nrm==0 ) THEN
         DO k = k1 , k2
            scal = max(scal,abs(A(k)))
         ENDDO
      ELSEIF ( Nrm==1 ) THEN
         DO k = k1 , k2
            scal = scal + abs(A(k))
         ENDDO
      ELSE
         DO k = k1 , k2
            scal = scal + A(k)**2
         ENDDO
      ENDIF
      IF ( Nrm==2 ) scal = sqrt(scal)
      Diag(ii) = scal
   ENDDO
!-----------------------------------------------------------------------
!-------------end-of-rnrms----------------------------------------------
END SUBROUTINE rnrms
!*==cnrms.f90 processed by SPAG 8.04RA 15:47 21 Feb 2024
!!SPAG Open source Personal, Educational or Academic User Clemson University  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
!-----------------------------------------------------------------------
SUBROUTINE cnrms(Nrow,Nrm,A,Ja,Ia,Diag)
   USE ISO_FORTRAN_ENV                 
   IMPLICIT NONE
!
! Dummy argument declarations rewritten by SPAG
!
   INTEGER , INTENT(IN) :: Nrow
   INTEGER , INTENT(IN) :: Nrm
   REAL(REAL64) , INTENT(IN) , DIMENSION(*) :: A
   INTEGER , INTENT(IN) , DIMENSION(*) :: Ja
   INTEGER , INTENT(IN) , DIMENSION(Nrow+1) :: Ia
   REAL(REAL64) , INTENT(INOUT) , DIMENSION(Nrow) :: Diag
!
! Local variable declarations rewritten by SPAG
!
   INTEGER :: ii , j , k , k1 , k2
!
! End of declarations rewritten by SPAG
!
!-----------------------------------------------------------------------
! gets the norms of each column of A. (choice of three norms)
!-----------------------------------------------------------------------
! on entry:
! ---------
! nrow	= integer. The row dimension of A
!
! nrm   = integer. norm indicator. nrm = 1, means 1-norm, nrm =2
!                  means the 2-nrm, nrm = 0 means max norm
!
! a,
! ja,
! ia   = Matrix A in compressed sparse row format.
!
! on return:
!----------
!
! diag = real vector of length nrow containing the norms
! NOTE: this is designed for square matrices. For the case
! nrow .ne. ncol -- diag must be of size max(nrow,ncol) even
! though only the ncol first entries will be filled..
! [report E. Canot 10/20/05 ]
!-----------------------------------------------------------------
   DO k = 1 , Nrow
      Diag(k) = 0.0D0
   ENDDO
   DO ii = 1 , Nrow
      k1 = Ia(ii)
      k2 = Ia(ii+1) - 1
      DO k = k1 , k2
         j = Ja(k)
!     update the norm of each column
         IF ( Nrm==0 ) THEN
            Diag(j) = max(Diag(j),abs(A(k)))
         ELSEIF ( Nrm==1 ) THEN
            Diag(j) = Diag(j) + abs(A(k))
         ELSE
            Diag(j) = Diag(j) + A(k)**2
         ENDIF
      ENDDO
   ENDDO
   IF ( Nrm/=2 ) RETURN
   DO k = 1 , Nrow
      Diag(k) = sqrt(Diag(k))
   ENDDO
!-----------------------------------------------------------------------
!------------end-of-cnrms-----------------------------------------------
END SUBROUTINE cnrms
!*==roscal.f90 processed by SPAG 8.04RA 15:47 21 Feb 2024
!!SPAG Open source Personal, Educational or Academic User Clemson University  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
!-----------------------------------------------------------------------
SUBROUTINE roscal(Nrow,Job,Nrm,A,Ja,Ia,Diag,B,Jb,Ib,Ierr)
   USE ISO_FORTRAN_ENV                 
   IMPLICIT NONE
!
! Dummy argument declarations rewritten by SPAG
!
   INTEGER :: Nrow
   INTEGER :: Job
   INTEGER :: Nrm
   REAL(REAL64) , DIMENSION(*) :: A
   INTEGER , DIMENSION(*) :: Ja
   INTEGER , DIMENSION(Nrow+1) :: Ia
   REAL(REAL64) , INTENT(INOUT) , DIMENSION(Nrow) :: Diag
   REAL(REAL64) , DIMENSION(*) :: B
   INTEGER , DIMENSION(*) :: Jb
   INTEGER , DIMENSION(Nrow+1) :: Ib
   INTEGER , INTENT(OUT) :: Ierr
!
! Local variable declarations rewritten by SPAG
!
   INTEGER :: j
   EXTERNAL diamua , rnrms
!
! End of declarations rewritten by SPAG
!
!-----------------------------------------------------------------------
! scales the rows of A such that their norms are one on return
! 3 choices of norms: 1-norm, 2-norm, max-norm.
!-----------------------------------------------------------------------
! on entry:
! ---------
! nrow	= integer. The row dimension of A
!
! job   = integer. job indicator. Job=0 means get array b only
!         job = 1 means get b, and the integer arrays ib, jb.
!
! nrm   = integer. norm indicator. nrm = 1, means 1-norm, nrm =2
!                  means the 2-nrm, nrm = 0 means max norm
!
! a,
! ja,
! ia   = Matrix A in compressed sparse row format.
!
! on return:
!----------
!
! diag = diagonal matrix stored as a vector containing the matrix
!        by which the rows have been scaled, i.e., on return
!        we have B = Diag*A.
!
! b,
! jb,
! ib	= resulting matrix B in compressed sparse row sparse format.
!
! ierr  = error message. ierr=0     : Normal return
!                        ierr=i > 0 : Row number i is a zero row.
! Notes:
!-------
! 1)        The column dimension of A is not needed.
! 2)        algorithm in place (B can take the place of A).
!-----------------------------------------------------------------
   CALL rnrms(Nrow,Nrm,A,Ia,Diag)
   Ierr = 0
   DO j = 1 , Nrow
      IF ( Diag(j)==0.0D0 ) THEN
         Ierr = j
         RETURN
      ELSE
         Diag(j) = 1.0D0/Diag(j)
      ENDIF
   ENDDO
   CALL diamua(Nrow,Job,A,Ja,Ia,Diag,B,Jb,Ib)
!-------end-of-roscal---------------------------------------------------
!-----------------------------------------------------------------------
END SUBROUTINE roscal
!*==coscal.f90 processed by SPAG 8.04RA 15:47 21 Feb 2024
!!SPAG Open source Personal, Educational or Academic User Clemson University  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
!-----------------------------------------------------------------------
SUBROUTINE coscal(Nrow,Job,Nrm,A,Ja,Ia,Diag,B,Jb,Ib,Ierr)
   USE ISO_FORTRAN_ENV                 
   IMPLICIT NONE
!
! Dummy argument declarations rewritten by SPAG
!
   INTEGER :: Nrow
   INTEGER :: Job
   INTEGER :: Nrm
   REAL(REAL64) , DIMENSION(*) :: A
   INTEGER , DIMENSION(*) :: Ja
   INTEGER , DIMENSION(Nrow+1) :: Ia
   REAL(REAL64) , INTENT(INOUT) , DIMENSION(Nrow) :: Diag
   REAL(REAL64) , DIMENSION(*) :: B
   INTEGER , DIMENSION(*) :: Jb
   INTEGER , DIMENSION(Nrow+1) :: Ib
   INTEGER , INTENT(OUT) :: Ierr
!
! Local variable declarations rewritten by SPAG
!
   INTEGER :: j
   EXTERNAL amudia , cnrms
!
! End of declarations rewritten by SPAG
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
! scales the columns of A such that their norms are one on return
! result matrix written on b, or overwritten on A.
! 3 choices of norms: 1-norm, 2-norm, max-norm. in place.
!-----------------------------------------------------------------------
! on entry:
! ---------
! nrow	= integer. The row dimension of A
!
! job   = integer. job indicator. Job=0 means get array b only
!         job = 1 means get b, and the integer arrays ib, jb.
!
! nrm   = integer. norm indicator. nrm = 1, means 1-norm, nrm =2
!                  means the 2-nrm, nrm = 0 means max norm
!
! a,
! ja,
! ia   = Matrix A in compressed sparse row format.
!
! on return:
!----------
!
! diag = diagonal matrix stored as a vector containing the matrix
!        by which the columns have been scaled, i.e., on return
!        we have B = A * Diag
!
! b,
! jb,
! ib	= resulting matrix B in compressed sparse row sparse format.
!
! ierr  = error message. ierr=0     : Normal return
!                        ierr=i > 0 : Column number i is a zero row.
! Notes:
!-------
! 1)     The column dimension of A is not needed.
! 2)     algorithm in place (B can take the place of A).
!-----------------------------------------------------------------
   CALL cnrms(Nrow,Nrm,A,Ja,Ia,Diag)
   Ierr = 0
   DO j = 1 , Nrow
      IF ( Diag(j)==0.0 ) THEN
         Ierr = j
         RETURN
      ELSE
         Diag(j) = 1.0D0/Diag(j)
      ENDIF
   ENDDO
   CALL amudia(Nrow,Job,A,Ja,Ia,Diag,B,Jb,Ib)
!--------end-of-coscal--------------------------------------------------
!-----------------------------------------------------------------------
END SUBROUTINE coscal
!*==addblk.f90 processed by SPAG 8.04RA 15:47 21 Feb 2024
!!SPAG Open source Personal, Educational or Academic User Clemson University  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
!-----------------------------------------------------------------------
SUBROUTINE addblk(Nrowa,Ncola,A,Ja,Ia,Ipos,Jpos,Job,Nrowb,Ncolb,B,Jb,Ib,Nrowc,Ncolc,C,Jc,Ic,Nzmx,Ierr)
   USE ISO_FORTRAN_ENV                 
   IMPLICIT NONE
!
! Dummy argument declarations rewritten by SPAG
!
   INTEGER , INTENT(IN) :: Nrowa
   INTEGER , INTENT(IN) :: Ncola
   REAL(REAL64) , INTENT(IN) , DIMENSION(1:*) :: A
   INTEGER , INTENT(IN) , DIMENSION(1:*) :: Ja
   INTEGER , INTENT(IN) , DIMENSION(1:*) :: Ia
   INTEGER , INTENT(IN) :: Ipos
   INTEGER , INTENT(IN) :: Jpos
   INTEGER , INTENT(IN) :: Job
   INTEGER , INTENT(IN) :: Nrowb
   INTEGER , INTENT(IN) :: Ncolb
   REAL(REAL64) , INTENT(IN) , DIMENSION(1:*) :: B
   INTEGER , INTENT(IN) , DIMENSION(1:*) :: Jb
   INTEGER , INTENT(IN) , DIMENSION(1:*) :: Ib
   INTEGER , INTENT(INOUT) :: Nrowc
   INTEGER , INTENT(INOUT) :: Ncolc
   REAL(REAL64) , INTENT(OUT) , DIMENSION(1:*) :: C
   INTEGER , INTENT(OUT) , DIMENSION(1:*) :: Jc
   INTEGER , INTENT(OUT) , DIMENSION(1:*) :: Ic
   INTEGER , INTENT(IN) :: Nzmx
   INTEGER , INTENT(OUT) :: Ierr
!
! Local variable declarations rewritten by SPAG
!
   INTEGER :: i , j1 , j2 , ka , kamax , kb , kbmax , kc
   LOGICAL :: values
!
! End of declarations rewritten by SPAG
!
!      implicit none
!-----------------------------------------------------------------------
!     This subroutine adds a matrix B into a submatrix of A whose
!     (1,1) element is located in the starting position (ipos, jpos).
!     The resulting matrix is allowed to be larger than A (and B),
!     and the resulting dimensions nrowc, ncolc will be redefined
!     accordingly upon return.
!     The input matrices are assumed to be sorted, i.e. in each row
!     the column indices appear in ascending order in the CSR format.
!-----------------------------------------------------------------------
! on entry:
! ---------
! nrowa    = number of rows in A.
! bcola    = number of columns in A.
! a,ja,ia  = Matrix A in compressed sparse row format with entries sorted
! nrowb    = number of rows in B.
! ncolb    = number of columns in B.
! b,jb,ib  = Matrix B in compressed sparse row format with entries sorted
!
! nzmax	   = integer. The  length of the arrays c and jc. addblk will
!            stop if the number of nonzero elements in the matrix C
!            exceeds nzmax. See ierr.
!
! on return:
!----------
! nrowc    = number of rows in C.
! ncolc    = number of columns in C.
! c,jc,ic  = resulting matrix C in compressed sparse row sparse format
!            with entries sorted ascendly in each row.
!
! ierr	   = integer. serving as error message.
!         ierr = 0 means normal return,
!         ierr .gt. 0 means that addblk stopped while computing the
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
         kamax = 0
 
         values = (Job/=0)
         Ierr = 0
         Nrowc = max(Nrowa,Nrowb+Ipos-1)
         Ncolc = max(Ncola,Ncolb+Jpos-1)
         kc = 1
         kbmax = 0
         Ic(1) = kc
!
         DO i = 1 , Nrowc
            IF ( i<=Nrowa ) THEN
               ka = Ia(i)
               kamax = Ia(i+1) - 1
            ELSE
               ka = Ia(Nrowa+1)
            ENDIF
            IF ( (i>=Ipos) .AND. ((i-Ipos)<=Nrowb) ) THEN
               kb = Ib(i-Ipos+1)
               kbmax = Ib(i-Ipos+2) - 1
            ELSE
               kb = Ib(Nrowb+1)
            ENDIF
            SPAG_Loop_3_1: DO
!
!     a do-while type loop -- goes through all the elements in a row.
!
               IF ( ka<=kamax ) THEN
                  j1 = Ja(ka)
               ELSE
                  j1 = Ncolc + 1
               ENDIF
               IF ( kb<=kbmax ) THEN
                  j2 = Jb(kb) + Jpos - 1
               ELSE
                  j2 = Ncolc + 1
               ENDIF
!
!     if there are more elements to be added.
!
               IF ( (ka<=kamax .OR. kb<=kbmax) .AND. (j1<=Ncolc .OR. j2<=Ncolc) ) THEN
!
!     three cases
!
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
                  IF ( kc<=Nzmx ) CYCLE
                  spag_nextblock_1 = 2
                  CYCLE SPAG_DispatchLoop_1
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
!---------end-of-addblk-------------------------------------------------
END SUBROUTINE addblk
!*==get1up.f90 processed by SPAG 8.04RA 15:47 21 Feb 2024
!!SPAG Open source Personal, Educational or Academic User Clemson University  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
!-----------------------------------------------------------------------
SUBROUTINE get1up(N,Ja,Ia,Ju)
   IMPLICIT NONE
!
! Dummy argument declarations rewritten by SPAG
!
   INTEGER , INTENT(IN) :: N
   INTEGER , INTENT(IN) , DIMENSION(*) :: Ja
   INTEGER , INTENT(IN) , DIMENSION(*) :: Ia
   INTEGER , INTENT(OUT) , DIMENSION(*) :: Ju
!
! Local variable declarations rewritten by SPAG
!
   INTEGER :: i , k
!
! End of declarations rewritten by SPAG
!
!----------------------------------------------------------------------
! obtains the first element of each row of the upper triangular part
! of a matrix. Assumes that the matrix is already sorted.
!-----------------------------------------------------------------------
! parameters
! input
! -----
! ja      = integer array containing the column indices of aij
! ia      = pointer array. ia(j) contains the position of the
!           beginning of row j in ja
!
! output
! ------
! ju      = integer array of length n. ju(i) is the address in ja
!           of the first element of the uper triangular part of
!           of A (including rthe diagonal. Thus if row i does have
!           a nonzero diagonal element then ju(i) will point to it.
!           This is a more general version of diapos.
!-----------------------------------------------------------------------
! local vAriables
!
   DO i = 1 , N
      Ju(i) = 0
      k = Ia(i)
      SPAG_Loop_2_1: DO
!
         IF ( Ja(k)>=i ) THEN
            Ju(i) = k
         ELSEIF ( k<Ia(i+1)-1 ) THEN
            k = k + 1
!
! go try next element in row
!
            CYCLE
         ENDIF
         EXIT SPAG_Loop_2_1
      ENDDO SPAG_Loop_2_1
   ENDDO
!-----end-of-get1up-----------------------------------------------------
END SUBROUTINE get1up
!*==xtrows.f90 processed by SPAG 8.04RA 15:47 21 Feb 2024
!!SPAG Open source Personal, Educational or Academic User Clemson University  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
!----------------------------------------------------------------------
SUBROUTINE xtrows(I1,I2,A,Ja,Ia,Ao,Jao,Iao,Iperm,Job)
   USE ISO_FORTRAN_ENV                 
   IMPLICIT NONE
!
! Dummy argument declarations rewritten by SPAG
!
   INTEGER , INTENT(IN) :: I1
   INTEGER , INTENT(IN) :: I2
   REAL(REAL64) , INTENT(IN) , DIMENSION(*) :: A
   INTEGER , INTENT(IN) , DIMENSION(*) :: Ja
   INTEGER , INTENT(IN) , DIMENSION(*) :: Ia
   REAL(REAL64) , INTENT(OUT) , DIMENSION(*) :: Ao
   INTEGER , INTENT(OUT) , DIMENSION(*) :: Jao
   INTEGER , INTENT(OUT) , DIMENSION(*) :: Iao
   INTEGER , INTENT(IN) , DIMENSION(*) :: Iperm
   INTEGER , INTENT(IN) :: Job
!
! Local variable declarations rewritten by SPAG
!
   INTEGER :: ii , j , k , ko
   LOGICAL :: values
!
! End of declarations rewritten by SPAG
!
!-----------------------------------------------------------------------
! this subroutine extracts given rows from a matrix in CSR format.
! Specifically, rows number iperm(i1), iperm(i1+1), ...., iperm(i2)
! are extracted and put in the output matrix ao, jao, iao, in CSR
! format.  NOT in place.
! Youcef Saad -- coded Feb 15, 1992.
!-----------------------------------------------------------------------
! on entry:
!----------
! i1,i2   = two integers indicating the rows to be extracted.
!           xtrows will extract rows iperm(i1), iperm(i1+1),..,iperm(i2),
!           from original matrix and stack them in output matrix
!           ao, jao, iao in csr format
!
! a, ja, ia = input matrix in csr format
!
! iperm	= integer array of length nrow containing the reverse permutation
!         array for the rows. row number iperm(j) in permuted matrix PA
!         used to be row number j in unpermuted matrix.
!         ---> a(i,j) in the permuted matrix was a(iperm(i),j)
!         in the inout matrix.
!
! job	= integer indicating the work to be done:
! 		job .ne. 1 : get structure only of output matrix,,
!               i.e., ignore real values. (in which case arrays a
!               and ao are not used nor accessed).
! 		job = 1	get complete data structure of output matrix.
!               (i.e., including arrays ao and iao).
!------------
! on return:
!------------
! ao, jao, iao = input matrix in a, ja, ia format
! note :
!        if (job.ne.1)  then the arrays a and ao are not used.
!----------------------------------------------------------------------c
!           Y. Saad, revised May  2, 1990                              c
!----------------------------------------------------------------------c
   values = (Job==1)
!
! copying
!
   ko = 1
   Iao(1) = ko
   DO j = I1 , I2
!
! ii=iperm(j) is the index of old row to be copied.
!
      ii = Iperm(j)
      DO k = Ia(ii) , Ia(ii+1) - 1
         Jao(ko) = Ja(k)
         IF ( values ) Ao(ko) = A(k)
         ko = ko + 1
      ENDDO
      Iao(j-I1+2) = ko
   ENDDO
!
!---------end-of-xtrows-------------------------------------------------
!-----------------------------------------------------------------------
END SUBROUTINE xtrows
!*==csrkvstr.f90 processed by SPAG 8.04RA 15:47 21 Feb 2024
!!SPAG Open source Personal, Educational or Academic User Clemson University  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
!-----------------------------------------------------------------------
SUBROUTINE csrkvstr(N,Ia,Ja,Nr,Kvstr)
   IMPLICIT NONE
!
! Dummy argument declarations rewritten by SPAG
!
   INTEGER , INTENT(IN) :: N
   INTEGER , INTENT(IN) , DIMENSION(N+1) :: Ia
   INTEGER , INTENT(IN) , DIMENSION(*) :: Ja
   INTEGER , INTENT(INOUT) :: Nr
   INTEGER , INTENT(OUT) , DIMENSION(*) :: Kvstr
!
! Local variable declarations rewritten by SPAG
!
   INTEGER :: i , j , jdiff
!
! End of declarations rewritten by SPAG
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!     Finds block row partitioning of matrix in CSR format.
!-----------------------------------------------------------------------
!     On entry:
!--------------
!     n       = number of matrix scalar rows
!     ia,ja   = input matrix sparsity structure in CSR format
!
!     On return:
!---------------
!     nr      = number of block rows
!     kvstr   = first row number for each block row
!
!     Notes:
!-----------
!     Assumes that the matrix is sorted by columns.
!     This routine does not need any workspace.
!
!-----------------------------------------------------------------------
!     local variables
!-----------------------------------------------------------------------
   Nr = 1
   Kvstr(1) = 1
!---------------------------------
   DO i = 2 , N
      jdiff = Ia(i+1) - Ia(i)
      IF ( jdiff==Ia(i)-Ia(i-1) ) THEN
         SPAG_Loop_2_1: DO j = Ia(i) , Ia(i+1) - 1
            IF ( Ja(j)/=Ja(j-jdiff) ) THEN
               Nr = Nr + 1
               Kvstr(Nr) = i
               EXIT SPAG_Loop_2_1
            ENDIF
         ENDDO SPAG_Loop_2_1
      ELSE
         Nr = Nr + 1
         Kvstr(Nr) = i
      ENDIF
   ENDDO
   Kvstr(Nr+1) = N + 1
!---------------------------------
END SUBROUTINE csrkvstr
!*==csrkvstc.f90 processed by SPAG 8.04RA 15:47 21 Feb 2024
!!SPAG Open source Personal, Educational or Academic User Clemson University  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
!-----------------------------------------------------------------------
!------------------------end-of-csrkvstr--------------------------------
SUBROUTINE csrkvstc(N,Ia,Ja,Nc,Kvstc,Iwk)
   IMPLICIT NONE
!
! Dummy argument declarations rewritten by SPAG
!
   INTEGER , INTENT(IN) :: N
   INTEGER , INTENT(IN) , DIMENSION(N+1) :: Ia
   INTEGER , INTENT(IN) , DIMENSION(*) :: Ja
   INTEGER , INTENT(INOUT) :: Nc
   INTEGER , INTENT(OUT) , DIMENSION(*) :: Kvstc
   INTEGER , INTENT(INOUT) , DIMENSION(*) :: Iwk
!
! Local variable declarations rewritten by SPAG
!
   INTEGER :: i , j , k , ncol
!
! End of declarations rewritten by SPAG
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!     Finds block column partitioning of matrix in CSR format.
!-----------------------------------------------------------------------
!     On entry:
!--------------
!     n       = number of matrix scalar rows
!     ia,ja   = input matrix sparsity structure in CSR format
!
!     On return:
!---------------
!     nc      = number of block columns
!     kvstc   = first column number for each block column
!
!     Work space:
!----------------
!     iwk(*) of size equal to the number of scalar columns plus one.
!        Assumed initialized to 0, and left initialized on return.
!
!     Notes:
!-----------
!     Assumes that the matrix is sorted by columns.
!
!-----------------------------------------------------------------------
!     local variables
!
!-----------------------------------------------------------------------
!-----use ncol to find maximum scalar column number
   ncol = 0
!-----mark the beginning position of the blocks in iwk
   DO i = 1 , N
      IF ( Ia(i)<Ia(i+1) ) THEN
         j = Ja(Ia(i))
         Iwk(j) = 1
         DO k = Ia(i) + 1 , Ia(i+1) - 1
            j = Ja(k)
            IF ( Ja(k-1)/=j-1 ) THEN
               Iwk(j) = 1
               Iwk(Ja(k-1)+1) = 1
            ENDIF
         ENDDO
         Iwk(j+1) = 1
         ncol = max0(ncol,j)
      ENDIF
   ENDDO
!---------------------------------
   Nc = 1
   Kvstc(1) = 1
   DO i = 2 , ncol + 1
      IF ( Iwk(i)/=0 ) THEN
         Nc = Nc + 1
         Kvstc(Nc) = i
         Iwk(i) = 0
      ENDIF
   ENDDO
   Nc = Nc - 1
!---------------------------------
END SUBROUTINE csrkvstc
!*==kvstmerge.f90 processed by SPAG 8.04RA 15:47 21 Feb 2024
!!SPAG Open source Personal, Educational or Academic User Clemson University  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
!-----------------------------------------------------------------------
!------------------------end-of-csrkvstc--------------------------------
!-----------------------------------------------------------------------
SUBROUTINE kvstmerge(Nr,Kvstr,Nc,Kvstc,N,Kvst)
   IMPLICIT NONE
!
! Dummy argument declarations rewritten by SPAG
!
   INTEGER , INTENT(IN) :: Nr
   INTEGER , INTENT(IN) :: Nc
   INTEGER , INTENT(IN) , DIMENSION(Nr+1) :: Kvstr
   INTEGER , INTENT(IN) , DIMENSION(Nc+1) :: Kvstc
   INTEGER , INTENT(INOUT) :: N
   INTEGER , INTENT(OUT) , DIMENSION(*) :: Kvst
!
! Local variable declarations rewritten by SPAG
!
   INTEGER :: i , j
!
! End of declarations rewritten by SPAG
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!     Merges block partitionings, for conformal row/col pattern.
!-----------------------------------------------------------------------
!     On entry:
!--------------
!     nr,nc   = matrix block row and block column dimension
!     kvstr   = first row number for each block row
!     kvstc   = first column number for each block column
!
!     On return:
!---------------
!     n       = conformal row/col matrix block dimension
!     kvst    = conformal row/col block partitioning
!
!     Notes:
!-----------
!     If matrix is not square, this routine returns without warning.
!
!-----------------------------------------------------------------------
!-----local variables
!---------------------------------
   IF ( Kvstr(Nr+1)/=Kvstc(Nc+1) ) RETURN
   i = 1
   j = 1
   N = 1
   SPAG_Loop_1_1: DO
      IF ( i>Nr+1 ) THEN
         Kvst(N) = Kvstc(j)
         j = j + 1
      ELSEIF ( j>Nc+1 ) THEN
         Kvst(N) = Kvstr(i)
         i = i + 1
      ELSEIF ( Kvstc(j)==Kvstr(i) ) THEN
         Kvst(N) = Kvstc(j)
         j = j + 1
         i = i + 1
      ELSEIF ( Kvstc(j)<Kvstr(i) ) THEN
         Kvst(N) = Kvstc(j)
         j = j + 1
      ELSE
         Kvst(N) = Kvstr(i)
         i = i + 1
      ENDIF
      N = N + 1
      IF ( i>Nr+1 .AND. j>Nc+1 ) THEN
         N = N - 2
         EXIT SPAG_Loop_1_1
      ENDIF
   ENDDO SPAG_Loop_1_1
!---------------------------------
!------------------------end-of-kvstmerge-------------------------------
END SUBROUTINE kvstmerge
