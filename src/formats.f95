!*==test.f90 processed by SPAG 8.04RA 15:47 21 Feb 2024
!!SPAG Open source Personal, Educational or Academic User Clemson University  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
!*==csrdns.f90 processed by SPAG 8.04RA 15:47 21 Feb 2024
!!SPAG Open source Personal, Educational or Academic User Clemson University  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
!----------------------------------------------------------------------c
!                          S P A R S K I T                             c
!----------------------------------------------------------------------c
!                    FORMAT CONVERSION MODULE                          c
!----------------------------------------------------------------------c
! contents:                                                            c
!----------                                                            c
! csrdns  : converts a row-stored sparse matrix into the dense format. c
! dnscsr  : converts a dense matrix to a sparse storage format.        c
! coocsr  : converts coordinate to  to csr format                      c
! coicsr  : in-place conversion of coordinate to csr format            c
! csrcoo  : converts compressed sparse row to coordinate.              c
! csrssr  : converts compressed sparse row to symmetric sparse row     c
! ssrcsr  : converts symmetric sparse row to compressed sparse row     c
! csrell  : converts compressed sparse row to ellpack format           c
! ellcsr  : converts ellpack format to compressed sparse row format    c
! csrmsr  : converts compressed sparse row format to modified sparse   c
!           row format                                                 c
! msrcsr  : converts modified sparse row format to compressed sparse   c
!           row format.                                                c
! csrcsc  : converts compressed sparse row format to compressed sparse c
!           column format (transposition)                              c
! csrcsc2 : rectangular version of csrcsc                              c
! csrlnk  : converts compressed sparse row to linked list format       c
! lnkcsr  : converts linked list format to compressed sparse row fmt   c
! csrdia  : converts a compressed sparse row format into a diagonal    c
!           format.                                                    c
! diacsr  : converts a diagonal format into a compressed sparse row    c
!           format.                                                    c
! bsrcsr  : converts a block-row sparse format into a compressed       c
!           sparse row format.                                         c
! csrbsr  : converts a compressed sparse row format into a block-row   c
!           sparse format.                                             c
! csrbnd  : converts a compressed sparse row format into a banded      c
!           format (linpack style).                                    c
! bndcsr  : converts a banded format (linpack style) into a compressed c
!           sparse row storage.                                        c
! csrssk  : converts the compressed sparse row format to the symmetric c
!           skyline format                                             c
! sskssr  : converts symmetric skyline format to symmetric  sparse row c
!           format.                                                    c
! csrjad  : converts the csr format into the jagged diagonal format    c
! jadcsr  : converts the jagged-diagonal format into the csr format    c
! csruss  : Compressed Sparse Row to Unsymmetric Sparse Skyline        c
!           format                                                     c
! usscsr  : Unsymmetric Sparse Skyline format to Compressed Sparse Row c
! csrsss  : Compressed Sparse Row to Symmetric Sparse Skyline format   c
! ssscsr  : Symmetric Sparse Skyline format to Compressed Sparse Row   c
! csrvbr  : Converts compressed sparse row to var block row format     c
! vbrcsr  : Converts var block row to compressed sparse row format     c
! csorted : Checks if matrix in CSR format is sorted by columns        c
!--------- miscalleneous additions not involving the csr format--------c
! cooell  : converts coordinate to Ellpack/Itpack format               c
! dcsort  : sorting routine used by crsjad                             c
!----------------------------------------------------------------------c
SUBROUTINE csrdns(Nrow,Ncol,A,Ja,Ia,Dns,Ndns,Ierr)
   USE ISO_FORTRAN_ENV                 
   IMPLICIT NONE
!
! Dummy argument declarations rewritten by SPAG
!
   INTEGER , INTENT(IN) :: Ndns
   INTEGER , INTENT(IN) :: Nrow
   INTEGER , INTENT(IN) :: Ncol
   REAL(REAL64) , INTENT(IN) , DIMENSION(*) :: A
   INTEGER , INTENT(IN) , DIMENSION(*) :: Ja
   INTEGER , INTENT(IN) , DIMENSION(*) :: Ia
   REAL(REAL64) , INTENT(OUT) , DIMENSION(Ndns,*) :: Dns
   INTEGER , INTENT(OUT) :: Ierr
!
! Local variable declarations rewritten by SPAG
!
   INTEGER :: i , j , k
!
! End of declarations rewritten by SPAG
!
!-----------------------------------------------------------------------
! Compressed Sparse Row    to    Dense
!-----------------------------------------------------------------------
!
! converts a row-stored sparse matrix into a densely stored one
!
! On entry:
!----------
!
! nrow	= row-dimension of a
! ncol	= column dimension of a
! a,
! ja,
! ia    = input matrix in compressed sparse row format.
!         (a=value array, ja=column array, ia=pointer array)
! dns   = array where to store dense matrix
! ndns	= first dimension of array dns
!
! on return:
!-----------
! dns   = the sparse matrix a, ja, ia has been stored in dns(ndns,*)
!
! ierr  = integer error indicator.
!         ierr .eq. 0  means normal return
!         ierr .eq. i  means that the code has stopped when processing
!         row number i, because it found a column number .gt. ncol.
!
!-----------------------------------------------------------------------
   Ierr = 0
   DO i = 1 , Nrow
      DO j = 1 , Ncol
         Dns(i,j) = 0.0D0
      ENDDO
   ENDDO
!
   DO i = 1 , Nrow
      DO k = Ia(i) , Ia(i+1) - 1
         j = Ja(k)
         IF ( j>Ncol ) THEN
            Ierr = i
            RETURN
         ENDIF
         Dns(i,j) = A(k)
      ENDDO
   ENDDO
!---- end of csrdns ----------------------------------------------------
!-----------------------------------------------------------------------
END SUBROUTINE csrdns
!*==dnscsr.f90 processed by SPAG 8.04RA 15:47 21 Feb 2024
!!SPAG Open source Personal, Educational or Academic User Clemson University  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
!-----------------------------------------------------------------------
SUBROUTINE dnscsr(Nrow,Ncol,Nzmax,Dns,Ndns,A,Ja,Ia,Ierr)
   USE ISO_FORTRAN_ENV                 
   IMPLICIT NONE
!
! Dummy argument declarations rewritten by SPAG
!
   INTEGER , INTENT(IN) :: Ndns
   INTEGER , INTENT(IN) :: Nrow
   INTEGER , INTENT(IN) :: Ncol
   INTEGER , INTENT(IN) :: Nzmax
   REAL(REAL64) , INTENT(IN) , DIMENSION(Ndns,*) :: Dns
   REAL(REAL64) , INTENT(OUT) , DIMENSION(*) :: A
   INTEGER , INTENT(OUT) , DIMENSION(*) :: Ja
   INTEGER , INTENT(OUT) , DIMENSION(*) :: Ia
   INTEGER , INTENT(OUT) :: Ierr
!
! Local variable declarations rewritten by SPAG
!
   INTEGER :: i , j , next
!
! End of declarations rewritten by SPAG
!
!-----------------------------------------------------------------------
! Dense		to    Compressed Row Sparse
!-----------------------------------------------------------------------
!
! converts a densely stored matrix into a row orientied
! compactly sparse matrix. ( reverse of csrdns )
! Note: this routine does not check whether an element
! is small. It considers that a(i,j) is zero if it is exactly
! equal to zero: see test below.
!-----------------------------------------------------------------------
! on entry:
!---------
!
! nrow	= row-dimension of a
! ncol	= column dimension of a
! nzmax = maximum number of nonzero elements allowed. This
!         should be set to be the lengths of the arrays a and ja.
! dns   = input nrow x ncol (dense) matrix.
! ndns	= first dimension of dns.
!
! on return:
!----------
!
! a, ja, ia = value, column, pointer  arrays for output matrix
!
! ierr	= integer error indicator:
!         ierr .eq. 0 means normal retur
!         ierr .eq. i means that the the code stopped while
!         processing row number i, because there was no space left in
!         a, and ja (as defined by parameter nzmax).
!-----------------------------------------------------------------------
   Ierr = 0
   next = 1
   Ia(1) = 1
   DO i = 1 , Nrow
      DO j = 1 , Ncol
         IF ( Dns(i,j)/=0.0D0 ) THEN
            IF ( next>Nzmax ) THEN
               Ierr = i
               RETURN
            ENDIF
            Ja(next) = j
            A(next) = Dns(i,j)
            next = next + 1
         ENDIF
      ENDDO
      Ia(i+1) = next
   ENDDO
!---- end of dnscsr ----------------------------------------------------
!-----------------------------------------------------------------------
END SUBROUTINE dnscsr
!*==coocsr.f90 processed by SPAG 8.04RA 15:47 21 Feb 2024
!!SPAG Open source Personal, Educational or Academic User Clemson University  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
!-----------------------------------------------------------------------
SUBROUTINE coocsr(Nrow,Nnz,A,Ir,Jc,Ao,Jao,Iao)
   USE ISO_FORTRAN_ENV                 
   IMPLICIT NONE
!
! Dummy argument declarations rewritten by SPAG
!
   INTEGER , INTENT(IN) :: Nrow
   INTEGER , INTENT(IN) :: Nnz
   REAL(REAL64) , INTENT(IN) , DIMENSION(*) :: A
   INTEGER , INTENT(IN) , DIMENSION(*) :: Ir
   INTEGER , INTENT(IN) , DIMENSION(*) :: Jc
   REAL(REAL64) , INTENT(OUT) , DIMENSION(*) :: Ao
   INTEGER , INTENT(OUT) , DIMENSION(*) :: Jao
   INTEGER , INTENT(INOUT) , DIMENSION(*) :: Iao
!
! Local variable declarations rewritten by SPAG
!
   INTEGER :: i , iad , j , k , k0
   REAL(REAL64) :: x
!
! End of declarations rewritten by SPAG
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!  Coordinate     to   Compressed Sparse Row
!-----------------------------------------------------------------------
! converts a matrix that is stored in coordinate format
!  a, ir, jc into a row general sparse ao, jao, iao format.
!
! on entry:
!---------
! nrow	= dimension of the matrix
! nnz	= number of nonzero elements in matrix
! a,
! ir,
! jc    = matrix in coordinate format. a(k), ir(k), jc(k) store the nnz
!         nonzero elements of the matrix with a(k) = actual real value of
! 	  the elements, ir(k) = its row number and jc(k) = its column
!	  number. The order of the elements is arbitrary.
!
! on return:
!-----------
! ir 	is destroyed
!
! ao, jao, iao = matrix in general sparse matrix format with ao
! 	continung the real values, jao containing the column indices,
!	and iao being the pointer to the beginning of the row,
!	in arrays ao, jao.
!
! Notes:
!------ This routine is NOT in place.  See coicsr
!       On return the entries  of each row are NOT sorted by increasing
!       column number
!
!------------------------------------------------------------------------
   DO k = 1 , Nrow + 1
      Iao(k) = 0
   ENDDO
! determine row-lengths.
   DO k = 1 , Nnz
      Iao(Ir(k)) = Iao(Ir(k)) + 1
   ENDDO
! starting position of each row..
   k = 1
   DO j = 1 , Nrow + 1
      k0 = Iao(j)
      Iao(j) = k
      k = k + k0
   ENDDO
! go through the structure  once more. Fill in output matrix.
   DO k = 1 , Nnz
      i = Ir(k)
      j = Jc(k)
      x = A(k)
      iad = Iao(i)
      Ao(iad) = x
      Jao(iad) = j
      Iao(i) = iad + 1
   ENDDO
! shift back iao
   DO j = Nrow , 1 , -1
      Iao(j+1) = Iao(j)
   ENDDO
   Iao(1) = 1
!------------- end of coocsr -------------------------------------------
!-----------------------------------------------------------------------
END SUBROUTINE coocsr
!*==coicsr.f90 processed by SPAG 8.04RA 15:47 21 Feb 2024
!!SPAG Open source Personal, Educational or Academic User Clemson University  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
!-----------------------------------------------------------------------
SUBROUTINE coicsr(N,Nnz,Job,A,Ja,Ia,Iwk)
   USE ISO_FORTRAN_ENV                 
   IMPLICIT NONE
!
! Dummy argument declarations rewritten by SPAG
!
   INTEGER , INTENT(IN) :: N
   INTEGER , INTENT(IN) :: Nnz
   INTEGER , INTENT(IN) :: Job
   REAL(REAL64) , INTENT(INOUT) , DIMENSION(*) :: A
   INTEGER , INTENT(INOUT) , DIMENSION(Nnz) :: Ja
   INTEGER , INTENT(INOUT) , DIMENSION(*) :: Ia
   INTEGER , INTENT(INOUT) , DIMENSION(N+1) :: Iwk
!
! Local variable declarations rewritten by SPAG
!
   INTEGER :: i , inext , init , ipos , j , jnext , k
   REAL(REAL64) :: t , tnext
   LOGICAL :: values
!
! End of declarations rewritten by SPAG
!
!------------------------------------------------------------------------
! IN-PLACE coo-csr conversion routine.
!------------------------------------------------------------------------
! this subroutine converts a matrix stored in coordinate format into
! the csr format. The conversion is done in place in that the arrays
! a,ja,ia of the result are overwritten onto the original arrays.
!------------------------------------------------------------------------
! on entry:
!---------
! n	= integer. row dimension of A.
! nnz	= integer. number of nonzero elements in A.
! job   = integer. Job indicator. when job=1, the real values in a are
!         filled. Otherwise a is not touched and the structure of the
!         array only (i.e. ja, ia)  is obtained.
! a	= real array of size nnz (number of nonzero elements in A)
!         containing the nonzero elements
! ja	= integer array of length nnz containing the column positions
! 	  of the corresponding elements in a.
! ia	= integer array of length max(nnz,n+1) containing the row
! 	  positions of the corresponding elements in a.
! iwk	= integer work array of length n+1
! on return:
!----------
! a
! ja
! ia	= contains the compressed sparse row data structure for the
!         resulting matrix.
! Note:
!-------
!         the entries of the output matrix are not sorted (the column
!         indices in each are not in increasing order) use coocsr
!         if you want them sorted. Note also that ia has to have at
!         least n+1 locations
!----------------------------------------------------------------------c
!  Coded by Y. Saad, Sep. 26 1989                                      c
!----------------------------------------------------------------------c
   INTEGER :: spag_nextblock_1
   spag_nextblock_1 = 1
   SPAG_DispatchLoop_1: DO
      SELECT CASE (spag_nextblock_1)
      CASE (1)
!-----------------------------------------------------------------------
         t = 0.0
         tnext = 0.0
         values = (Job==1)
! find pointer array for resulting matrix.
         DO i = 1 , N + 1
            Iwk(i) = 0
         ENDDO
         DO k = 1 , Nnz
            i = Ia(k)
            Iwk(i+1) = Iwk(i+1) + 1
         ENDDO
!------------------------------------------------------------------------
         Iwk(1) = 1
         DO i = 2 , N
            Iwk(i) = Iwk(i-1) + Iwk(i)
         ENDDO
!
!     loop for a cycle in chasing process.
!
         init = 1
         k = 0
         spag_nextblock_1 = 2
      CASE (2)
         IF ( values ) t = A(init)
         i = Ia(init)
         j = Ja(init)
         Ia(init) = -1
         DO
!------------------------------------------------------------------------
            k = k + 1
!     current row number is i.  determine  where to go.
            ipos = Iwk(i)
!     save the chased element.
            IF ( values ) tnext = A(ipos)
            inext = Ia(ipos)
            jnext = Ja(ipos)
!     then occupy its location.
            IF ( values ) A(ipos) = t
            Ja(ipos) = j
!     update pointer information for next element to come in row i.
            Iwk(i) = ipos + 1
!     determine  next element to be chased,
            IF ( Ia(ipos)<0 ) THEN
               DO
                  init = init + 1
                  IF ( init>Nnz ) THEN
                     spag_nextblock_1 = 3
                     CYCLE SPAG_DispatchLoop_1
                  ENDIF
                  IF ( Ia(init)>=0 ) THEN
                     spag_nextblock_1 = 2
                     CYCLE SPAG_DispatchLoop_1
!     restart chasing --
                  ENDIF
               ENDDO
            ELSE
               t = tnext
               i = inext
               j = jnext
               Ia(ipos) = -1
               IF ( k>=Nnz ) THEN
                  spag_nextblock_1 = 3
                  CYCLE SPAG_DispatchLoop_1
               ENDIF
            ENDIF
         ENDDO
         spag_nextblock_1 = 3
      CASE (3)
         DO i = 1 , N
            Ia(i+1) = Iwk(i)
         ENDDO
         Ia(1) = 1
         EXIT SPAG_DispatchLoop_1
      END SELECT
   ENDDO SPAG_DispatchLoop_1
!----------------- end of coicsr ----------------------------------------
!------------------------------------------------------------------------
END SUBROUTINE coicsr
!*==csrcoo.f90 processed by SPAG 8.04RA 15:47 21 Feb 2024
!!SPAG Open source Personal, Educational or Academic User Clemson University  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
!-----------------------------------------------------------------------
SUBROUTINE csrcoo(Nrow,Job,Nzmax,A,Ja,Ia,Nnz,Ao,Ir,Jc,Ierr)
   USE ISO_FORTRAN_ENV                 
   IMPLICIT NONE
!
! Dummy argument declarations rewritten by SPAG
!
   INTEGER , INTENT(IN) :: Nrow
   INTEGER , INTENT(IN) :: Job
   INTEGER , INTENT(IN) :: Nzmax
   REAL(REAL64) , INTENT(IN) , DIMENSION(*) :: A
   INTEGER , INTENT(IN) , DIMENSION(*) :: Ja
   INTEGER , INTENT(IN) , DIMENSION(Nrow+1) :: Ia
   INTEGER , INTENT(INOUT) :: Nnz
   REAL(REAL64) , INTENT(OUT) , DIMENSION(*) :: Ao
   INTEGER , INTENT(OUT) , DIMENSION(*) :: Ir
   INTEGER , INTENT(OUT) , DIMENSION(*) :: Jc
   INTEGER , INTENT(OUT) :: Ierr
!
! Local variable declarations rewritten by SPAG
!
   INTEGER :: i , k , k1 , k2
!
! End of declarations rewritten by SPAG
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!  Compressed Sparse Row      to      Coordinate
!-----------------------------------------------------------------------
! converts a matrix that is stored in coordinate format
!  a, ir, jc into a row general sparse ao, jao, iao format.
!
! on entry:
!---------
! nrow	= dimension of the matrix.
! job   = integer serving as a job indicator.
!         if job = 1 fill in only the array ir, ignore jc, and ao.
!         if job = 2 fill in ir, and jc but not ao
!         if job = 3 fill in everything.
!         The reason why these options are provided is that on return
!         ao and jc are the same as a, ja. So when job = 3, a and ja are
!         simply copied into ao, jc.  When job=2, only jc and ir are
!         returned. With job=1 only the array ir is returned. Moreover,
!         the algorithm is in place:
!	     call csrcoo (nrow,1,nzmax,a,ja,ia,nnz,a,ia,ja,ierr)
!         will write the output matrix in coordinate format on a, ja,ia.
!
! a,
! ja,
! ia    = matrix in compressed sparse row format.
! nzmax = length of space available in ao, ir, jc.
!         the code will stop immediatly if the number of
!         nonzero elements found in input matrix exceeds nzmax.
!
! on return:
!-----------
! ao, ir, jc = matrix in coordinate format.
!
! nnz        = number of nonzero elements in matrix.
! ierr       = integer error indicator.
!         ierr .eq. 0 means normal retur
!         ierr .eq. 1 means that the the code stopped
!         because there was no space in ao, ir, jc
!         (according to the value of  nzmax).
!
! NOTES: 1)This routine is PARTIALLY in place: csrcoo can be called with
!         ao being the same array as as a, and jc the same array as ja.
!         but ir CANNOT be the same as ia.
!         2) note the order in the output arrays,
!------------------------------------------------------------------------
   Ierr = 0
   Nnz = Ia(Nrow+1) - 1
   IF ( Nnz>Nzmax ) THEN
      Ierr = 1
      RETURN
   ENDIF
!------------------------------------------------------------------------
   IF ( Job==1 ) THEN
      CALL spag_block_1
      RETURN
   ENDIF
   IF ( Job/=2 ) THEN
      DO k = 1 , Nnz
         Ao(k) = A(k)
      ENDDO
   ENDIF
   DO k = 1 , Nnz
      Jc(k) = Ja(k)
   ENDDO
   CALL spag_block_1
CONTAINS
   SUBROUTINE spag_block_1
!
!     copy backward to allow for in-place processing.
!
      DO i = Nrow , 1 , -1
         k1 = Ia(i+1) - 1
         k2 = Ia(i)
         DO k = k1 , k2 , -1
            Ir(k) = i
         ENDDO
      ENDDO
   END SUBROUTINE spag_block_1
!------------- end-of-csrcoo -------------------------------------------
!-----------------------------------------------------------------------
END SUBROUTINE csrcoo
!*==csrssr.f90 processed by SPAG 8.04RA 15:47 21 Feb 2024
!!SPAG Open source Personal, Educational or Academic User Clemson University  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
!-----------------------------------------------------------------------
SUBROUTINE csrssr(Nrow,A,Ja,Ia,Nzmax,Ao,Jao,Iao,Ierr)
   USE ISO_FORTRAN_ENV                 
   IMPLICIT NONE
!
! Dummy argument declarations rewritten by SPAG
!
   INTEGER , INTENT(IN) :: Nrow
   REAL(REAL64) , INTENT(IN) , DIMENSION(*) :: A
   INTEGER , INTENT(IN) , DIMENSION(*) :: Ja
   INTEGER , INTENT(IN) , DIMENSION(*) :: Ia
   INTEGER , INTENT(IN) :: Nzmax
   REAL(REAL64) , INTENT(INOUT) , DIMENSION(*) :: Ao
   INTEGER , INTENT(INOUT) , DIMENSION(*) :: Jao
   INTEGER , INTENT(OUT) , DIMENSION(*) :: Iao
   INTEGER , INTENT(OUT) :: Ierr
!
! Local variable declarations rewritten by SPAG
!
   INTEGER :: i , k , kdiag , ko , kold
   REAL(REAL64) :: t
!
! End of declarations rewritten by SPAG
!
!-----------------------------------------------------------------------
! Compressed Sparse Row     to     Symmetric Sparse Row
!-----------------------------------------------------------------------
! this subroutine extracts the lower triangular part of a matrix.
! It can used as a means for converting a symmetric matrix for
! which all the entries are stored in sparse format into one
! in which only the lower part is stored. The routine is in place in
! that the output matrix ao, jao, iao can be overwritten on
! the  input matrix  a, ja, ia if desired. Csrssr has been coded to
! put the diagonal elements of the matrix in the last position in
! each row (i.e. in position  ao(ia(i+1)-1   of ao and jao)
!-----------------------------------------------------------------------
! On entry
!-----------
! nrow  = dimension of the matrix a.
! a, ja,
!    ia = matrix stored in compressed row sparse format
!
! nzmax = length of arrays ao,  and jao.
!
! On return:
!-----------
! ao, jao,
!     iao = lower part of input matrix (a,ja,ia) stored in compressed sparse
!          row format format.
!
! ierr   = integer error indicator.
!          ierr .eq. 0  means normal return
!          ierr .eq. i  means that the code has stopped when processing
!          row number i, because there is not enough space in ao, jao
!          (according to the value of nzmax)
!
!-----------------------------------------------------------------------
   Ierr = 0
   ko = 0
!-----------------------------------------------------------------------
   DO i = 1 , Nrow
      kold = ko
      kdiag = 0
      DO k = Ia(i) , Ia(i+1) - 1
         IF ( Ja(k)<=i ) THEN
            ko = ko + 1
            IF ( ko>Nzmax ) THEN
               Ierr = i
               RETURN
            ENDIF
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
   Iao(Nrow+1) = ko + 1
!--------- end of csrssr -----------------------------------------------
!-----------------------------------------------------------------------
END SUBROUTINE csrssr
!*==ssrcsr.f90 processed by SPAG 8.04RA 15:47 21 Feb 2024
!!SPAG Open source Personal, Educational or Academic User Clemson University  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
!-----------------------------------------------------------------------
SUBROUTINE ssrcsr(Job,Value2,Nrow,A,Ja,Ia,Nzmax,Ao,Jao,Iao,Indu,Iwk,Ierr)
   USE ISO_FORTRAN_ENV                 
   IMPLICIT NONE
!
! Dummy argument declarations rewritten by SPAG
!
   INTEGER , INTENT(IN) :: Nrow
   INTEGER , INTENT(IN) :: Nzmax
   INTEGER , INTENT(IN) :: Job
   INTEGER , INTENT(IN) :: Value2
   REAL(REAL64) , INTENT(IN) , DIMENSION(*) :: A
   INTEGER , INTENT(IN) , DIMENSION(*) :: Ja
   INTEGER , INTENT(IN) , DIMENSION(Nrow+1) :: Ia
   REAL(REAL64) , INTENT(INOUT) , DIMENSION(Nzmax) :: Ao
   INTEGER , INTENT(INOUT) , DIMENSION(Nzmax) :: Jao
   INTEGER , INTENT(INOUT) , DIMENSION(Nrow+1) :: Iao
   INTEGER , INTENT(INOUT) , DIMENSION(Nrow) :: Indu
   INTEGER , INTENT(INOUT) , DIMENSION(Nrow+1) :: Iwk
   INTEGER , INTENT(OUT) :: Ierr
!
! Local variable declarations rewritten by SPAG
!
   INTEGER :: i , ipos , j , k , kfirst , klast , ko , kosav , nnz
   REAL(REAL64) :: tmp
!
! End of declarations rewritten by SPAG
!
!     .. Scalar Arguments ..
!     ..
!     .. Array Arguments ..
!     ..
!-----------------------------------------------------------------------
!     Symmetric Sparse Row to Compressed Sparse Row format
!-----------------------------------------------------------------------
!     This subroutine converts a given matrix in SSR format to regular
!     CSR format by computing Ao = A + A' - diag(A), where A' is A
!     transpose.
!
!     Typically this routine is used to expand the SSR matrix of
!     Harwell Boeing matrices, or to obtain a symmetrized graph of
!     unsymmetric matrices.
!
!     This routine is inplace, i.e., (Ao,jao,iao) may be same as
!     (a,ja,ia).
!
!     It is possible to input an arbitrary CSR matrix to this routine,
!     since there is no syntactical difference between CSR and SSR
!     format. It also removes duplicate entries and perform a partial
!     ordering. The output matrix has an order of lower half, main
!     diagonal and upper half after the partial ordering.
!-----------------------------------------------------------------------
! on entry:
!---------
!
! job   = options
!         0 -- duplicate entries are not removed. If the input matrix is
!              SSR (not an arbitary CSR) matrix, no duplicate entry should
!              arise from this routine.
!         1 -- eliminate duplicate entries, zero entries.
!         2 -- eliminate duplicate entries and perform partial ordering.
!         3 -- eliminate duplicate entries, sort the entries in the
!              increasing order of clumn indices.
!
! value2= will the values of A be copied?
!         0 -- only expand the graph (a, ao are not touched)
!         1 -- expand the matrix with the values.
!
! nrow  = column dimension of inout matrix
! a,
! ia,
! ja    = matrix in compressed sparse row format.
!
! nzmax = size of arrays ao and jao. SSRCSR will abort if the storage
!          provided in ao, jao is not sufficient to store A. See ierr.
!
! on return:
!----------
! ao, jao, iao
!       = output matrix in compressed sparse row format. The resulting
!         matrix is symmetric and is equal to A+A'-D. ao, jao, iao,
!         can be the same as a, ja, ia in the calling sequence.
!
! indu  = integer array of length nrow. INDU will contain pointers
!         to the beginning of upper traigular part if job > 1.
!         Otherwise it is also used as a work array (size nrow).
!
! iwk   = integer work space (size nrow+1).
!
! ierr  = integer. Serving as error message. If the length of the arrays
!         ao, jao exceeds nzmax, ierr returns the minimum value
!         needed for nzmax. otherwise ierr=0 (normal return).
!
!-----------------------------------------------------------------------
!     .. Local Scalars ..
!     ..
!     .. Executable Statements ..
   Ierr = 0
   DO i = 1 , Nrow
      Indu(i) = 0
      Iwk(i) = 0
   ENDDO
   Iwk(Nrow+1) = 0
!
!     .. compute number of elements in each row of (A'-D)
!     put result in iwk(i+1)  for row i.
!
   DO i = 1 , Nrow
      DO k = Ia(i) , Ia(i+1) - 1
         j = Ja(k)
         IF ( j/=i ) Iwk(j+1) = Iwk(j+1) + 1
      ENDDO
   ENDDO
!
!     .. find addresses of first elements of ouput matrix. result in iwk
!
   Iwk(1) = 1
   DO i = 1 , Nrow
      Indu(i) = Iwk(i) + Ia(i+1) - Ia(i)
      Iwk(i+1) = Iwk(i+1) + Indu(i)
      Indu(i) = Indu(i) - 1
   ENDDO
!.....Have we been given enough storage in ao, jao ?
   nnz = Iwk(Nrow+1) - 1
   IF ( nnz>Nzmax ) THEN
      Ierr = nnz
      RETURN
   ENDIF
!
!     .. copy the existing matrix (backwards).
!
   kosav = Iwk(Nrow+1)
   DO i = Nrow , 1 , -1
      klast = Ia(i+1) - 1
      kfirst = Ia(i)
      Iao(i+1) = kosav
      kosav = Iwk(i)
      ko = Iwk(i) - kfirst
      Iwk(i) = ko + klast + 1
      DO k = klast , kfirst , -1
         IF ( Value2/=0 ) Ao(k+ko) = A(k)
         Jao(k+ko) = Ja(k)
      ENDDO
   ENDDO
   Iao(1) = 1
!
!     now copy (A'-D). Go through the structure of ao, jao, iao
!     that has already been copied. iwk(i) is the address
!     of the next free location in row i for ao, jao.
!
   DO i = 1 , Nrow
      DO k = Iao(i) , Indu(i)
         j = Jao(k)
         IF ( j/=i ) THEN
            ipos = Iwk(j)
            IF ( Value2/=0 ) Ao(ipos) = Ao(k)
            Jao(ipos) = i
            Iwk(j) = ipos + 1
         ENDIF
      ENDDO
   ENDDO
   IF ( Job<=0 ) RETURN
!
!     .. eliminate duplicate entries --
!     array INDU is used as marker for existing indices, it is also the
!     location of the entry.
!     IWK is used to stored the old IAO array.
!     matrix is copied to squeeze out the space taken by the duplicated
!     entries.
!
   DO i = 1 , Nrow
      Indu(i) = 0
      Iwk(i) = Iao(i)
   ENDDO
   Iwk(Nrow+1) = Iao(Nrow+1)
   k = 1
   DO i = 1 , Nrow
      Iao(i) = k
      ipos = Iwk(i)
      klast = Iwk(i+1)
      SPAG_Loop_2_1: DO
         IF ( ipos<klast ) THEN
            j = Jao(ipos)
            IF ( Indu(j)==0 ) THEN
!     .. new entry ..
               IF ( Value2==0 ) THEN
                  Indu(j) = k
                  Jao(k) = Jao(ipos)
                  k = k + 1
               ELSEIF ( Ao(ipos)/=0.0D0 ) THEN
                  Indu(j) = k
                  Jao(k) = Jao(ipos)
                  Ao(k) = Ao(ipos)
                  k = k + 1
               ENDIF
            ELSEIF ( Value2/=0 ) THEN
!     .. duplicate entry ..
               Ao(Indu(j)) = Ao(Indu(j)) + Ao(ipos)
            ENDIF
            ipos = ipos + 1
            CYCLE
         ENDIF
!     .. remove marks before working on the next row ..
         DO ipos = Iao(i) , k - 1
            Indu(Jao(ipos)) = 0
         ENDDO
         EXIT SPAG_Loop_2_1
      ENDDO SPAG_Loop_2_1
   ENDDO
   Iao(Nrow+1) = k
   IF ( Job<=1 ) RETURN
!
!     .. partial ordering ..
!     split the matrix into strict upper/lower triangular
!     parts, INDU points to the the beginning of the strict upper part.
!
   DO i = 1 , Nrow
      klast = Iao(i+1) - 1
      kfirst = Iao(i)
      SPAG_Loop_2_2: DO
         IF ( klast>kfirst ) THEN
            IF ( Jao(klast)<i .AND. Jao(kfirst)>=i ) THEN
!     .. swap klast with kfirst ..
               j = Jao(klast)
               Jao(klast) = Jao(kfirst)
               Jao(kfirst) = j
               IF ( Value2/=0 ) THEN
                  tmp = Ao(klast)
                  Ao(klast) = Ao(kfirst)
                  Ao(kfirst) = tmp
               ENDIF
            ENDIF
            IF ( Jao(klast)>=i ) klast = klast - 1
            IF ( Jao(kfirst)<i ) kfirst = kfirst + 1
            CYCLE
         ENDIF
!
         IF ( Jao(klast)<i ) THEN
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
!     bubble-sort is used
!
   DO i = 1 , Nrow
      DO ipos = Iao(i) , Indu(i) - 1
         DO j = Indu(i) - 1 , ipos + 1 , -1
            k = j - 1
            IF ( Jao(k)>Jao(j) ) THEN
               ko = Jao(k)
               Jao(k) = Jao(j)
               Jao(j) = ko
               IF ( Value2/=0 ) THEN
                  tmp = Ao(k)
                  Ao(k) = Ao(j)
                  Ao(j) = tmp
               ENDIF
            ENDIF
         ENDDO
      ENDDO
      DO ipos = Indu(i) , Iao(i+1) - 1
         DO j = Iao(i+1) - 1 , ipos + 1 , -1
            k = j - 1
            IF ( Jao(k)>Jao(j) ) THEN
               ko = Jao(k)
               Jao(k) = Jao(j)
               Jao(j) = ko
               IF ( Value2/=0 ) THEN
                  tmp = Ao(k)
                  Ao(k) = Ao(j)
                  Ao(j) = tmp
               ENDIF
            ENDIF
         ENDDO
      ENDDO
   ENDDO
!
!---- end of ssrcsr ----------------------------------------------------
!-----------------------------------------------------------------------
END SUBROUTINE ssrcsr
!*==xssrcsr.f90 processed by SPAG 8.04RA 15:47 21 Feb 2024
!!SPAG Open source Personal, Educational or Academic User Clemson University  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
!-----------------------------------------------------------------------
SUBROUTINE xssrcsr(Nrow,A,Ja,Ia,Nzmax,Ao,Jao,Iao,Indu,Ierr)
   USE ISO_FORTRAN_ENV                 
   IMPLICIT NONE
!
! Dummy argument declarations rewritten by SPAG
!
   INTEGER , INTENT(IN) :: Nrow
   INTEGER , INTENT(IN) :: Nzmax
   REAL(REAL64) , INTENT(IN) , DIMENSION(*) :: A
   INTEGER , INTENT(IN) , DIMENSION(*) :: Ja
   INTEGER , INTENT(IN) , DIMENSION(Nrow+1) :: Ia
   REAL(REAL64) , INTENT(INOUT) , DIMENSION(Nzmax) :: Ao
   INTEGER , INTENT(INOUT) , DIMENSION(Nzmax) :: Jao
   INTEGER , INTENT(INOUT) , DIMENSION(Nrow+1) :: Iao
   INTEGER , INTENT(INOUT) , DIMENSION(Nrow+1) :: Indu
   INTEGER , INTENT(OUT) :: Ierr
!
! Local variable declarations rewritten by SPAG
!
   INTEGER :: i , ipos , j , k , kfirst , klast , ko , kosav , lenrow , nnz
!
! End of declarations rewritten by SPAG
!
!-----------------------------------------------------------------------
! Symmetric Sparse Row   to    (regular) Compressed Sparse Row
!-----------------------------------------------------------------------
! this subroutine converts  a symmetric  matrix in which only the lower
! part is  stored in compressed sparse row format, i.e.,
! a matrix stored in symmetric sparse format, into a fully stored matrix
! i.e., a matrix where both the lower and upper parts are stored in
! compressed sparse row format. the algorithm is in place (i.e. result
! may be overwritten onto the input matrix a, ja, ia ----- ).
! the output matrix delivered by ssrcsr is such that each row starts with
! the elements of the lower part followed by those of the upper part.
!-----------------------------------------------------------------------
! on entry:
!---------
!
! nrow  = row dimension of inout matrix
! a,
! ia,
! ja    = matrix in compressed sparse row format. This is assumed to be
!         a lower triangular matrix.
!
! nzmax	= size of arrays ao and jao. ssrcsr will abort if the storage
!	   provided in a, ja is not sufficient to store A. See ierr.
!
! on return:
!----------
! ao, iao,
!   jao = output matrix in compressed sparse row format. The resulting
!         matrix is symmetric and is equal to A+A**T - D, if
!         A is the original lower triangular matrix. ao, jao, iao,
!         can be the same as a, ja, ia in the calling sequence.
!
! indu  = integer array of length nrow+1. If the input matrix is such
!         that the last element in each row is its diagonal element then
!         on return, indu will contain the pointers to the diagonal
!         element in each row of the output matrix. Otherwise used as
!         work array.
! ierr  = integer. Serving as error message. If the length of the arrays
!         ao, jao exceeds nzmax, ierr returns the minimum value
!         needed for nzmax. otherwise ierr=0 (normal return).
!
!-----------------------------------------------------------------------
   Ierr = 0
   DO i = 1 , Nrow + 1
      Indu(i) = 0
   ENDDO
!
!     compute  number of elements in each row of strict upper part.
!     put result in indu(i+1)  for row i.
!
   DO i = 1 , Nrow
      DO k = Ia(i) , Ia(i+1) - 1
         j = Ja(k)
         IF ( j<i ) Indu(j+1) = Indu(j+1) + 1
      ENDDO
   ENDDO
!-----------
!     find addresses of first elements of ouput matrix. result in indu
!-----------
   Indu(1) = 1
   DO i = 1 , Nrow
      lenrow = Ia(i+1) - Ia(i)
      Indu(i+1) = Indu(i) + Indu(i+1) + lenrow
   ENDDO
!--------------------- enough storage in a, ja ? --------
   nnz = Indu(Nrow+1) - 1
   IF ( nnz>Nzmax ) THEN
      Ierr = nnz
      RETURN
   ENDIF
!
!     now copy lower part (backwards).
!
   kosav = Indu(Nrow+1)
   DO i = Nrow , 1 , -1
      klast = Ia(i+1) - 1
      kfirst = Ia(i)
      Iao(i+1) = kosav
      ko = Indu(i)
      kosav = ko
      DO k = kfirst , klast
         Ao(ko) = A(k)
         Jao(ko) = Ja(k)
         ko = ko + 1
      ENDDO
      Indu(i) = ko
   ENDDO
   Iao(1) = 1
!
!     now copy upper part. Go through the structure of ao, jao, iao
!     that has already been copied (lower part). indu(i) is the address
!     of the next free location in row i for ao, jao.
!
   DO i = 1 , Nrow
!     i-th row is now in ao, jao, iao structure -- lower half part
      SPAG_Loop_2_1: DO k = Iao(i) , Iao(i+1) - 1
         j = Jao(k)
         IF ( j>=i ) EXIT SPAG_Loop_2_1
         ipos = Indu(j)
         Ao(ipos) = Ao(k)
         Jao(ipos) = i
         Indu(j) = Indu(j) + 1
      ENDDO SPAG_Loop_2_1
   ENDDO
!----- end of xssrcsr --------------------------------------------------
!-----------------------------------------------------------------------
END SUBROUTINE xssrcsr
!*==csrell.f90 processed by SPAG 8.04RA 15:47 21 Feb 2024
!!SPAG Open source Personal, Educational or Academic User Clemson University  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
!-----------------------------------------------------------------------
SUBROUTINE csrell(Nrow,A,Ja,Ia,Maxcol,Coef,Jcoef,Ncoef,Ndiag,Ierr)
   USE ISO_FORTRAN_ENV                 
   IMPLICIT NONE
!
! Dummy argument declarations rewritten by SPAG
!
   INTEGER , INTENT(IN) :: Nrow
   INTEGER , INTENT(IN) :: Ncoef
   REAL(REAL64) , INTENT(IN) , DIMENSION(*) :: A
   INTEGER , INTENT(IN) , DIMENSION(*) :: Ja
   INTEGER , INTENT(IN) , DIMENSION(Nrow+1) :: Ia
   INTEGER , INTENT(IN) :: Maxcol
   REAL(REAL64) , INTENT(OUT) , DIMENSION(Ncoef,1) :: Coef
   INTEGER , INTENT(OUT) , DIMENSION(Ncoef,1) :: Jcoef
   INTEGER , INTENT(INOUT) :: Ndiag
   INTEGER , INTENT(OUT) :: Ierr
!
! Local variable declarations rewritten by SPAG
!
   INTEGER :: i , j , k , k1 , k2
!
! End of declarations rewritten by SPAG
!
!-----------------------------------------------------------------------
! Compressed Sparse Row	    to    Ellpack - Itpack format
!-----------------------------------------------------------------------
! this subroutine converts  matrix stored in the general a, ja, ia
! format into the coef, jcoef itpack format.
!
!-----------------------------------------------------------------------
! on entry:
!----------
! nrow 	  = row dimension of the matrix A.
!
! a,
! ia,
! ja      = input matrix in compressed sparse row format.
!
! ncoef  = first dimension of arrays coef, and jcoef.
!
! maxcol = integer equal to the number of columns available in coef.
!
! on return:
!----------
! coef	= real array containing the values of the matrix A in
!         itpack-ellpack format.
! jcoef = integer array containing the column indices of coef(i,j)
!         in A.
! ndiag = number of active 'diagonals' found.
!
! ierr 	= error message. 0 = correct return. If ierr .ne. 0 on
!	  return this means that the number of diagonals found
!         (ndiag) exceeds maxcol.
!
!-----------------------------------------------------------------------
! first determine the length of each row of lower-part-of(A)
   Ierr = 0
   Ndiag = 0
   DO i = 1 , Nrow
      k = Ia(i+1) - Ia(i)
      Ndiag = max0(Ndiag,k)
   ENDDO
!----- check whether sufficient columns are available. -----------------
   IF ( Ndiag>Maxcol ) THEN
      Ierr = 1
      RETURN
   ENDIF
!
! fill coef with zero elements and jcoef with row numbers.------------
!
   DO j = 1 , Ndiag
      DO i = 1 , Nrow
         Coef(i,j) = 0.0D0
         Jcoef(i,j) = i
      ENDDO
   ENDDO
!
!------- copy elements row by row.--------------------------------------
!
   DO i = 1 , Nrow
      k1 = Ia(i)
      k2 = Ia(i+1) - 1
      DO k = k1 , k2
         Coef(i,k-k1+1) = A(k)
         Jcoef(i,k-k1+1) = Ja(k)
      ENDDO
   ENDDO
!--- end of csrell------------------------------------------------------
!-----------------------------------------------------------------------
END SUBROUTINE csrell
!*==ellcsr.f90 processed by SPAG 8.04RA 15:47 21 Feb 2024
!!SPAG Open source Personal, Educational or Academic User Clemson University  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
!-----------------------------------------------------------------------
SUBROUTINE ellcsr(Nrow,Coef,Jcoef,Ncoef,Ndiag,A,Ja,Ia,Nzmax,Ierr)
   USE ISO_FORTRAN_ENV                 
   IMPLICIT NONE
!
! Dummy argument declarations rewritten by SPAG
!
   INTEGER , INTENT(IN) :: Nrow
   INTEGER , INTENT(IN) :: Ncoef
   REAL(REAL64) , INTENT(IN) , DIMENSION(Ncoef,1) :: Coef
   INTEGER , INTENT(IN) , DIMENSION(Ncoef,1) :: Jcoef
   INTEGER , INTENT(IN) :: Ndiag
   REAL(REAL64) , INTENT(OUT) , DIMENSION(*) :: A
   INTEGER , INTENT(OUT) , DIMENSION(*) :: Ja
   INTEGER , INTENT(OUT) , DIMENSION(Nrow+1) :: Ia
   INTEGER , INTENT(IN) :: Nzmax
   INTEGER , INTENT(OUT) :: Ierr
!
! Local variable declarations rewritten by SPAG
!
   INTEGER :: i , k , kpos
!
! End of declarations rewritten by SPAG
!
!-----------------------------------------------------------------------
!  Ellpack - Itpack format  to  Compressed Sparse Row
!-----------------------------------------------------------------------
! this subroutine converts a matrix stored in ellpack-itpack format
! coef-jcoef into the compressed sparse row format. It actually checks
! whether an entry in the input matrix is a nonzero element before
! putting it in the output matrix. The test does not account for small
! values but only for exact zeros.
!-----------------------------------------------------------------------
! on entry:
!----------
!
! nrow 	= row dimension of the matrix A.
! coef	= array containing the values of the matrix A in ellpack format.
! jcoef = integer arraycontains the column indices of coef(i,j) in A.
! ncoef = first dimension of arrays coef, and jcoef.
! ndiag = number of active columns in coef, jcoef.
!
! ndiag = on entry the number of columns made available in coef.
!
! on return:
!----------
! a, ia,
!    ja = matrix in a, ia, ja format where.
!
! nzmax	= size of arrays a and ja. ellcsr will abort if the storage
!	   provided in a, ja is not sufficient to store A. See ierr.
!
! ierr 	= integer. serves are output error message.
!         ierr = 0 means normal return.
!         ierr = 1 means that there is not enough space in
!         a and ja to store output matrix.
!-----------------------------------------------------------------------
! first determine the length of each row of lower-part-of(A)
   Ierr = 0
!-----check whether sufficient columns are available. -----------------
!
!------- copy elements row by row.--------------------------------------
   kpos = 1
   Ia(1) = kpos
   DO i = 1 , Nrow
      DO k = 1 , Ndiag
         IF ( Coef(i,k)/=0.0D0 ) THEN
            IF ( kpos>Nzmax ) THEN
               Ierr = kpos
               RETURN
            ENDIF
            A(kpos) = Coef(i,k)
            Ja(kpos) = Jcoef(i,k)
            kpos = kpos + 1
         ENDIF
      ENDDO
      Ia(i+1) = kpos
   ENDDO
!--- end of ellcsr -----------------------------------------------------
!-----------------------------------------------------------------------
END SUBROUTINE ellcsr
!*==csrmsr.f90 processed by SPAG 8.04RA 15:47 21 Feb 2024
!!SPAG Open source Personal, Educational or Academic User Clemson University  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
!-----------------------------------------------------------------------
SUBROUTINE csrmsr(N,A,Ja,Ia,Ao,Jao,Wk,Iwk)
   USE ISO_FORTRAN_ENV                 
   IMPLICIT NONE
!
! Dummy argument declarations rewritten by SPAG
!
   INTEGER , INTENT(IN) :: N
   REAL(REAL64) , INTENT(IN) , DIMENSION(*) :: A
   INTEGER , INTENT(IN) , DIMENSION(*) :: Ja
   INTEGER , INTENT(IN) , DIMENSION(N+1) :: Ia
   REAL(REAL64) , INTENT(OUT) , DIMENSION(*) :: Ao
   INTEGER , INTENT(INOUT) , DIMENSION(*) :: Jao
   REAL(REAL64) , INTENT(INOUT) , DIMENSION(N) :: Wk
   INTEGER , INTENT(INOUT) , DIMENSION(N+1) :: Iwk
!
! Local variable declarations rewritten by SPAG
!
   INTEGER :: i , icount , ii , iptr , j , k
!
! End of declarations rewritten by SPAG
!
!-----------------------------------------------------------------------
! Compressed Sparse Row   to      Modified - Sparse Row
!                                 Sparse row with separate main diagonal
!-----------------------------------------------------------------------
! converts a general sparse matrix a, ja, ia into
! a compressed matrix using a separated diagonal (referred to as
! the bell-labs format as it is used by bell labs semi conductor
! group. We refer to it here as the modified sparse row format.
! Note: this has been coded in such a way that one can overwrite
! the output matrix onto the input matrix if desired by a call of
! the form
!
!     call csrmsr (n, a, ja, ia, a, ja, wk,iwk)
!
! In case ao, jao, are different from a, ja, then one can
! use ao, jao as the work arrays in the calling sequence:
!
!     call csrmsr (n, a, ja, ia, ao, jao, ao,jao)
!
!-----------------------------------------------------------------------
!
! on entry :
!---------
! a, ja, ia = matrix in csr format. note that the
!	     algorithm is in place: ao, jao can be the same
!            as a, ja, in which case it will be overwritten on it
!            upon return.
!
! on return :
!-----------
!
! ao, jao  = sparse matrix in modified sparse row storage format:
!	   +  ao(1:n) contains the diagonal of the matrix.
!	   +  ao(n+2:nnz) contains the nondiagonal elements of the
!             matrix, stored rowwise.
!	   +  jao(n+2:nnz) : their column indices
!	   +  jao(1:n+1) contains the pointer array for the nondiagonal
!             elements in ao(n+1:nnz) and jao(n+2:nnz).
!             i.e., for i .le. n+1 jao(i) points to beginning of row i
!	      in arrays ao, jao.
!	       here nnz = number of nonzero elements+1
! work arrays:
!------------
! wk	= real work array of length n
! iwk   = integer work array of length n+1
!
! notes:
!-------
!        Algorithm is in place.  i.e. both:
!
!          call csrmsr (n, a, ja, ia, ao, jao, ao,jao)
!          (in which  ao, jao, are different from a, ja)
!           and
!          call csrmsr (n, a, ja, ia, a, ja, wk,iwk)
!          (in which  wk, jwk, are different from a, ja)
!        are OK.
!--------
! coded by Y. Saad Sep. 1989. Rechecked Feb 27, 1990.
!-----------------------------------------------------------------------
   icount = 0
!
! store away diagonal elements and count nonzero diagonal elements.
!
   DO i = 1 , N
      Wk(i) = 0.0D0
      Iwk(i+1) = Ia(i+1) - Ia(i)
      DO k = Ia(i) , Ia(i+1) - 1
         IF ( Ja(k)==i ) THEN
            Wk(i) = A(k)
            icount = icount + 1
            Iwk(i+1) = Iwk(i+1) - 1
         ENDIF
      ENDDO
   ENDDO
!
! compute total length
!
   iptr = N + Ia(N+1) - icount
!
!     copy backwards (to avoid collisions)
!
   DO ii = N , 1 , -1
      DO k = Ia(ii+1) - 1 , Ia(ii) , -1
         j = Ja(k)
         IF ( j/=ii ) THEN
            Ao(iptr) = A(k)
            Jao(iptr) = j
            iptr = iptr - 1
         ENDIF
      ENDDO
   ENDDO
!
! compute pointer values and copy wk(*)
!
   Jao(1) = N + 2
   DO i = 1 , N
      Ao(i) = Wk(i)
      Jao(i+1) = Jao(i) + Iwk(i+1)
   ENDDO
!------------ end of subroutine csrmsr ---------------------------------
!-----------------------------------------------------------------------
END SUBROUTINE csrmsr
!*==msrcsr.f90 processed by SPAG 8.04RA 15:47 21 Feb 2024
!!SPAG Open source Personal, Educational or Academic User Clemson University  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
!-----------------------------------------------------------------------
SUBROUTINE msrcsr(N,A,Ja,Ao,Jao,Iao,Wk,Iwk)
   USE ISO_FORTRAN_ENV                 
   IMPLICIT NONE
!
! Dummy argument declarations rewritten by SPAG
!
   INTEGER , INTENT(IN) :: N
   REAL(REAL64) , INTENT(IN) , DIMENSION(*) :: A
   INTEGER , INTENT(IN) , DIMENSION(*) :: Ja
   REAL(REAL64) , INTENT(OUT) , DIMENSION(*) :: Ao
   INTEGER , INTENT(OUT) , DIMENSION(*) :: Jao
   INTEGER , INTENT(OUT) , DIMENSION(N+1) :: Iao
   REAL(REAL64) , INTENT(INOUT) , DIMENSION(N) :: Wk
   INTEGER , INTENT(INOUT) , DIMENSION(N+1) :: Iwk
!
! Local variable declarations rewritten by SPAG
!
   LOGICAL :: added
   INTEGER :: i , idiag , ii , iptr , j , k
!
! End of declarations rewritten by SPAG
!
!-----------------------------------------------------------------------
!       Modified - Sparse Row  to   Compressed Sparse Row
!
!-----------------------------------------------------------------------
! converts a compressed matrix using a separated diagonal
! (modified sparse row format) in the Compressed Sparse Row
! format.
! does not check for zero elements in the diagonal.
!
!
! on entry :
!---------
! n          = row dimension of matrix
! a, ja      = sparse matrix in msr sparse storage format
!              see routine csrmsr for details on data structure
!
! on return :
!-----------
!
! ao,jao,iao = output matrix in csr format.
!
! work arrays:
!------------
! wk       = real work array of length n
! iwk      = integer work array of length n+1
!
! notes:
!   The original version of this was NOT in place, but has
!   been modified by adding the vector iwk to be in place.
!   The original version had ja instead of iwk everywhere in
!   loop 500.  Modified  Sun 29 May 1994 by R. Bramley (Indiana).
!
!-----------------------------------------------------------------------
   DO i = 1 , N
      Wk(i) = A(i)
      Iwk(i) = Ja(i)
   ENDDO
   Iwk(N+1) = Ja(N+1)
   Iao(1) = 1
   iptr = 1
!---------
   DO ii = 1 , N
      added = .FALSE.
      idiag = iptr + (Iwk(ii+1)-Iwk(ii))
      DO k = Iwk(ii) , Iwk(ii+1) - 1
         j = Ja(k)
         IF ( j<ii ) THEN
            Ao(iptr) = A(k)
            Jao(iptr) = j
            iptr = iptr + 1
         ELSEIF ( added ) THEN
            Ao(iptr) = A(k)
            Jao(iptr) = j
            iptr = iptr + 1
         ELSE
! add diag element - only reserve a position for it.
            idiag = iptr
            iptr = iptr + 1
            added = .TRUE.
!     then other element
            Ao(iptr) = A(k)
            Jao(iptr) = j
            iptr = iptr + 1
         ENDIF
      ENDDO
      Ao(idiag) = Wk(ii)
      Jao(idiag) = ii
      IF ( .NOT.added ) iptr = iptr + 1
      Iao(ii+1) = iptr
   ENDDO
!------------ end of subroutine msrcsr ---------------------------------
!-----------------------------------------------------------------------
END SUBROUTINE msrcsr
!*==csrcsc.f90 processed by SPAG 8.04RA 15:47 21 Feb 2024
!!SPAG Open source Personal, Educational or Academic User Clemson University  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
!-----------------------------------------------------------------------
SUBROUTINE csrcsc(N,Job,Ipos,A,Ja,Ia,Ao,Jao,Iao)
   USE ISO_FORTRAN_ENV                 
   IMPLICIT NONE
!
! Dummy argument declarations rewritten by SPAG
!
   INTEGER :: N
   INTEGER :: Job
   INTEGER :: Ipos
   REAL(REAL64) , DIMENSION(*) :: A
   INTEGER , DIMENSION(*) :: Ja
   INTEGER , DIMENSION(N+1) :: Ia
   REAL(REAL64) , DIMENSION(*) :: Ao
   INTEGER , DIMENSION(*) :: Jao
   INTEGER , DIMENSION(N+1) :: Iao
   EXTERNAL csrcsc2
!
! End of declarations rewritten by SPAG
!
!-----------------------------------------------------------------------
! Compressed Sparse Row     to      Compressed Sparse Column
!
! (transposition operation)   Not in place.
!-----------------------------------------------------------------------
! -- not in place --
! this subroutine transposes a matrix stored in a, ja, ia format.
! ---------------
! on entry:
!----------
! n	= dimension of A.
! job	= integer to indicate whether to fill the values (job.eq.1) of the
!         matrix ao or only the pattern., i.e.,ia, and ja (job .ne.1)
!
! ipos  = starting position in ao, jao of the transposed matrix.
!         the iao array takes this into account (thus iao(1) is set to ipos.)
!         Note: this may be useful if one needs to append the data structure
!         of the transpose to that of A. In this case use for example
!                call csrcsc (n,1,ia(n+1),a,ja,ia,a,ja,ia(n+2))
!	  for any other normal usage, enter ipos=1.
! a	= real array of length nnz (nnz=number of nonzero elements in input
!         matrix) containing the nonzero elements.
! ja	= integer array of length nnz containing the column positions
! 	  of the corresponding elements in a.
! ia	= integer of size n+1. ia(k) contains the position in a, ja of
!	  the beginning of the k-th row.
!
! on return:
! ----------
! output arguments:
! ao	= real array of size nzz containing the "a" part of the transpose
! jao	= integer array of size nnz containing the column indices.
! iao	= integer array of size n+1 containing the "ia" index array of
!	  the transpose.
!
!-----------------------------------------------------------------------
   CALL csrcsc2(N,N,Job,Ipos,A,Ja,Ia,Ao,Jao,Iao)
END SUBROUTINE csrcsc
!*==csrcsc2.f90 processed by SPAG 8.04RA 15:47 21 Feb 2024
!!SPAG Open source Personal, Educational or Academic User Clemson University  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
!-----------------------------------------------------------------------
SUBROUTINE csrcsc2(N,N2,Job,Ipos,A,Ja,Ia,Ao,Jao,Iao)
   USE ISO_FORTRAN_ENV                 
   IMPLICIT NONE
!
! Dummy argument declarations rewritten by SPAG
!
   INTEGER , INTENT(IN) :: N
   INTEGER , INTENT(IN) :: N2
   INTEGER , INTENT(IN) :: Job
   INTEGER , INTENT(IN) :: Ipos
   REAL(REAL64) , INTENT(IN) , DIMENSION(*) :: A
   INTEGER , INTENT(IN) , DIMENSION(*) :: Ja
   INTEGER , INTENT(IN) , DIMENSION(N+1) :: Ia
   REAL(REAL64) , INTENT(OUT) , DIMENSION(*) :: Ao
   INTEGER , INTENT(OUT) , DIMENSION(*) :: Jao
   INTEGER , INTENT(INOUT) , DIMENSION(N2+1) :: Iao
!
! Local variable declarations rewritten by SPAG
!
   INTEGER :: i , j , k , next
!
! End of declarations rewritten by SPAG
!
!-----------------------------------------------------------------------
! Compressed Sparse Row     to      Compressed Sparse Column
!
! (transposition operation)   Not in place.
!-----------------------------------------------------------------------
! Rectangular version.  n is number of rows of CSR matrix,
!                       n2 (input) is number of columns of CSC matrix.
!-----------------------------------------------------------------------
! -- not in place --
! this subroutine transposes a matrix stored in a, ja, ia format.
! ---------------
! on entry:
!----------
! n	= number of rows of CSR matrix.
! n2    = number of columns of CSC matrix.
! job	= integer to indicate whether to fill the values (job.eq.1) of the
!         matrix ao or only the pattern., i.e.,ia, and ja (job .ne.1)
!
! ipos  = starting position in ao, jao of the transposed matrix.
!         the iao array takes this into account (thus iao(1) is set to ipos.)
!         Note: this may be useful if one needs to append the data structure
!         of the transpose to that of A. In this case use for example
!                call csrcsc2 (n,n,1,ia(n+1),a,ja,ia,a,ja,ia(n+2))
!	  for any other normal usage, enter ipos=1.
! a	= real array of length nnz (nnz=number of nonzero elements in input
!         matrix) containing the nonzero elements.
! ja	= integer array of length nnz containing the column positions
! 	  of the corresponding elements in a.
! ia	= integer of size n+1. ia(k) contains the position in a, ja of
!	  the beginning of the k-th row.
!
! on return:
! ----------
! output arguments:
! ao	= real array of size nzz containing the "a" part of the transpose
! jao	= integer array of size nnz containing the column indices.
! iao	= integer array of size n+1 containing the "ia" index array of
!	  the transpose.
!
!-----------------------------------------------------------------------
!----------------- compute lengths of rows of transp(A) ----------------
   DO i = 1 , N2 + 1
      Iao(i) = 0
   ENDDO
   DO i = 1 , N
      DO k = Ia(i) , Ia(i+1) - 1
         j = Ja(k) + 1
         Iao(j) = Iao(j) + 1
      ENDDO
   ENDDO
!---------- compute pointers from lengths ------------------------------
   Iao(1) = Ipos
   DO i = 1 , N2
      Iao(i+1) = Iao(i) + Iao(i+1)
   ENDDO
!--------------- now do the actual copying -----------------------------
   DO i = 1 , N
      DO k = Ia(i) , Ia(i+1) - 1
         j = Ja(k)
         next = Iao(j)
         IF ( Job==1 ) Ao(next) = A(k)
         Jao(next) = i
         Iao(j) = next + 1
      ENDDO
   ENDDO
!-------------------------- reshift iao and leave ----------------------
   DO i = N2 , 1 , -1
      Iao(i+1) = Iao(i)
   ENDDO
   Iao(1) = Ipos
!--------------- end of csrcsc2 ----------------------------------------
!-----------------------------------------------------------------------
END SUBROUTINE csrcsc2
!*==csrlnk.f90 processed by SPAG 8.04RA 15:47 21 Feb 2024
!!SPAG Open source Personal, Educational or Academic User Clemson University  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
!-----------------------------------------------------------------------
SUBROUTINE csrlnk(N,Ia,Link)
   IMPLICIT NONE
!
! Dummy argument declarations rewritten by SPAG
!
   INTEGER , INTENT(IN) :: N
   INTEGER , INTENT(INOUT) , DIMENSION(N+1) :: Ia
   INTEGER , INTENT(OUT) , DIMENSION(*) :: Link
!
! Local variable declarations rewritten by SPAG
!
   INTEGER :: i , iend , istart , k
!
! End of declarations rewritten by SPAG
!
!-----------------------------------------------------------------------
!      Compressed Sparse Row         to    Linked storage format.
!-----------------------------------------------------------------------
! this subroutine translates a matrix stored in compressed sparse
! row into one with a linked list storage format. Only the link
! array needs to be obtained since the arrays a, ja, and ia may
! be unchanged and  carry the same meaning for the output matrix.
! in  other words a, ja, ia, link   is the output linked list data
! structure with a, ja, unchanged from input, and ia possibly
! altered (in case therea re null rows in matrix). Details on
! the output array link are given below.
!-----------------------------------------------------------------------
! Coded by Y. Saad, Feb 21, 1991.
!-----------------------------------------------------------------------
!
! on entry:
!----------
! n	= integer equal to the dimension of A.
! a, ja, ia == matrix in CSR format  -- only ia is used.
! a	= real array of size nna containing the nonzero elements
! ja	= integer array of size	nnz containing the column positions
! 	  of the corresponding elements in a.
! ia	= integer of size n+1 containing the pointers to the beginning
!         of each row. ia(k) contains the position in a, ja of the
!         beginning of the k-th row.
!
! on return:
!----------
! a, ja, are not changed.
! ia    may be changed if there are null rows.
!
! a     = nonzero elements.
! ja    = column positions.
! ia    = ia(i) points to the first element of row i in linked structure.
! link	= integer array of size containing the linked list information.
!         link(k) points to the next element of the row after element
!         a(k), ja(k). if link(k) = 0, then there is no next element,
!         i.e., a(k), jcol(k) is the last element of the current row.
!
!  Thus row number i can be accessed as follows:
!     next = ia(i)
!     while(next .ne. 0) do
!          value = a(next)      ! value a(i,j)
!          jcol  = ja(next)     ! column index j
!          next  = link(next)   ! address of next element in row
!     endwhile
! notes:
! ------ ia may be altered on return.
!-----------------------------------------------------------------------
! local variables
!
! loop through all rows
!
   DO i = 1 , N
      istart = Ia(i)
      iend = Ia(i+1) - 1
      IF ( iend>istart ) THEN
         DO k = istart , iend - 1
            Link(k) = k + 1
         ENDDO
         Link(iend) = 0
      ELSE
         Ia(i) = 0
      ENDIF
   ENDDO
!
!-------------end-of-csrlnk --------------------------------------------
!-----------------------------------------------------------------------
END SUBROUTINE csrlnk
!*==lnkcsr.f90 processed by SPAG 8.04RA 15:47 21 Feb 2024
!!SPAG Open source Personal, Educational or Academic User Clemson University  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
!-----------------------------------------------------------------------
SUBROUTINE lnkcsr(N,A,Jcol,Istart,Link,Ao,Jao,Iao)
   USE ISO_FORTRAN_ENV                 
   IMPLICIT NONE
!
! Dummy argument declarations rewritten by SPAG
!
   INTEGER , INTENT(IN) :: N
   REAL(REAL64) , INTENT(IN) , DIMENSION(*) :: A
   INTEGER , INTENT(IN) , DIMENSION(*) :: Jcol
   INTEGER , INTENT(IN) , DIMENSION(N) :: Istart
   INTEGER , INTENT(IN) , DIMENSION(*) :: Link
   REAL(REAL64) , INTENT(OUT) , DIMENSION(*) :: Ao
   INTEGER , INTENT(OUT) , DIMENSION(*) :: Jao
   INTEGER , INTENT(OUT) , DIMENSION(*) :: Iao
!
! Local variable declarations rewritten by SPAG
!
   INTEGER :: ipos , irow , next
!
! End of declarations rewritten by SPAG
!
!-----------------------------------------------------------------------
!     Linked list storage format   to      Compressed Sparse Row  format
!-----------------------------------------------------------------------
! this subroutine translates a matrix stored in linked list storage
! format into the compressed sparse row format.
!-----------------------------------------------------------------------
! Coded by Y. Saad, Feb 21, 1991.
!-----------------------------------------------------------------------
!
! on entry:
!----------
! n	= integer equal to the dimension of A.
!
! a	= real array of size nna containing the nonzero elements
! jcol	= integer array of size	nnz containing the column positions
! 	  of the corresponding elements in a.
! istart= integer array of size n poiting to the beginning of the rows.
!         istart(i) contains the position of the first element of
!         row i in data structure. (a, jcol, link).
!         if a row is empty istart(i) must be zero.
! link	= integer array of size nnz containing the links in the linked
!         list data structure. link(k) points to the next element
!         of the row after element ao(k), jcol(k). if link(k) = 0,
!         then there is no next element, i.e., ao(k), jcol(k) is
!         the last element of the current row.
!
! on return:
!-----------
! ao, jao, iao = matrix stored in csr format:
!
! ao    = real array containing the values of the nonzero elements of
!         the matrix stored row-wise.
! jao	= integer array of size nnz containing the column indices.
! iao	= integer array of size n+1 containing the pointers array to the
!         beginning of each row. iao(i) is the address in ao,jao of
!         first element of row i.
!
!-----------------------------------------------------------------------
! first determine individial bandwidths and pointers.
!-----------------------------------------------------------------------
! local variables
!-----------------------------------------------------------------------
   ipos = 1
   Iao(1) = ipos
!
!     loop through all rows
!
   DO irow = 1 , N
!
!     unroll i-th row.
!
      next = Istart(irow)
      DO WHILE ( next/=0 )
         Jao(ipos) = Jcol(next)
         Ao(ipos) = A(next)
         ipos = ipos + 1
         next = Link(next)
      ENDDO
      Iao(irow+1) = ipos
   ENDDO
!
!-------------end-of-lnkcsr -------------------------------------------
!-----------------------------------------------------------------------
END SUBROUTINE lnkcsr
!*==csrdia.f90 processed by SPAG 8.04RA 15:47 21 Feb 2024
!!SPAG Open source Personal, Educational or Academic User Clemson University  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
!-----------------------------------------------------------------------
SUBROUTINE csrdia(N,Idiag,Job,A,Ja,Ia,Ndiag,Diag,Ioff,Ao,Jao,Iao,Ind)
   USE ISO_FORTRAN_ENV                 
   IMPLICIT NONE
!
! Dummy argument declarations rewritten by SPAG
!
   INTEGER , INTENT(INOUT) :: Idiag
   INTEGER , INTENT(IN) :: Ndiag
   INTEGER :: N
   INTEGER , INTENT(IN) :: Job
   REAL(REAL64) , INTENT(IN) , DIMENSION(*) :: A
   INTEGER , DIMENSION(*) :: Ja
   INTEGER , DIMENSION(*) :: Ia
   REAL(REAL64) , INTENT(OUT) , DIMENSION(Ndiag,Idiag) :: Diag
   INTEGER , INTENT(INOUT) , DIMENSION(*) :: Ioff
   REAL(REAL64) , INTENT(OUT) , DIMENSION(*) :: Ao
   INTEGER , INTENT(OUT) , DIMENSION(*) :: Jao
   INTEGER , INTENT(OUT) , DIMENSION(*) :: Iao
   INTEGER , INTENT(INOUT) , DIMENSION(*) :: Ind
!
! Local variable declarations rewritten by SPAG
!
   INTEGER :: i , idum , ii , j , jmax , job1 , job2 , k , ko , l , n2
   EXTERNAL infdia
!
! End of declarations rewritten by SPAG
!
!-----------------------------------------------------------------------
! Compressed sparse row     to    diagonal format
!-----------------------------------------------------------------------
! this subroutine extracts  idiag diagonals  from the  input matrix a,
! a, ia, and puts the rest of  the matrix  in the  output matrix ao,
! jao, iao.  The diagonals to be extracted depend  on the  value of job
! (see below for details.)  In  the first  case, the  diagonals to be
! extracted are simply identified by  their offsets  provided in ioff
! by the caller.  In the second case, the  code internally determines
! the idiag most significant diagonals, i.e., those  diagonals of the
! matrix which  have  the  largest  number  of  nonzero elements, and
! extracts them.
!-----------------------------------------------------------------------
! on entry:
!----------
! n	= dimension of the matrix a.
! idiag = integer equal to the number of diagonals to be extracted.
!         Note: on return idiag may be modified.
! a, ja,
!    ia = matrix stored in a, ja, ia, format
! job	= integer. serves as a job indicator.  Job is better thought
!         of as a two-digit number job=xy. If the first (x) digit
!         is one on entry then the diagonals to be extracted are
!         internally determined. In this case csrdia exctracts the
!         idiag most important diagonals, i.e. those having the largest
!         number on nonzero elements. If the first digit is zero
!         then csrdia assumes that ioff(*) contains the offsets
!         of the diagonals to be extracted. there is no verification
!         that ioff(*) contains valid entries.
!         The second (y) digit of job determines whether or not
!         the remainder of the matrix is to be written on ao,jao,iao.
!         If it is zero  then ao, jao, iao is not filled, i.e.,
!         the diagonals are found  and put in array diag and the rest is
!         is discarded. if it is one, ao, jao, iao contains matrix
!         of the remaining elements.
!         Thus:
!         job= 0 means do not select diagonals internally (pick those
!                defined by ioff) and do not fill ao,jao,iao
!         job= 1 means do not select diagonals internally
!                      and fill ao,jao,iao
!         job=10 means  select diagonals internally
!                      and do not fill ao,jao,iao
!         job=11 means select diagonals internally
!                      and fill ao,jao,iao
!
! ndiag = integer equal to the first dimension of array diag.
!
! on return:
!-----------
!
! idiag = number of diagonals found. This may be smaller than its value
!         on entry.
! diag  = real array of size (ndiag x idiag) containing the diagonals
!         of A on return
!
! ioff  = integer array of length idiag, containing the offsets of the
!   	  diagonals to be extracted.
! ao, jao
!  iao  = remainder of the matrix in a, ja, ia format.
! work arrays:
!------------
! ind   = integer array of length 2*n-1 used as integer work space.
!         needed only when job.ge.10 i.e., in case the diagonals are to
!         be selected internally.
!
! Notes:
!-------
!    1) The algorithm is in place: ao, jao, iao can be overwritten on
!       a, ja, ia if desired
!    2) When the code is required to select the diagonals (job .ge. 10)
!       the selection of the diagonals is done from left to right
!       as a result if several diagonals have the same weight (number
!       of nonzero elemnts) the leftmost one is selected first.
!-----------------------------------------------------------------------
   job1 = Job/10
   job2 = Job - job1*10
   IF ( job1/=0 ) THEN
      n2 = N + N - 1
      CALL infdia(N,Ja,Ia,Ind,idum)
!----------- determine diagonals to  accept.----------------------------
!-----------------------------------------------------------------------
      ii = 0
      SPAG_Loop_1_1: DO
         ii = ii + 1
         jmax = 0
         DO k = 1 , n2
            j = Ind(k)
            IF ( j>jmax ) THEN
               i = k
               jmax = j
            ENDIF
         ENDDO
         IF ( jmax<=0 ) THEN
            ii = ii - 1
            EXIT SPAG_Loop_1_1
         ENDIF
         Ioff(ii) = i - N
         Ind(i) = -jmax
         IF ( ii>=Idiag ) EXIT SPAG_Loop_1_1
      ENDDO SPAG_Loop_1_1
      Idiag = ii
   ENDIF
!---------------- initialize diago to zero -----------------------------
   DO j = 1 , Idiag
      DO i = 1 , N
         Diag(i,j) = 0.0D0
      ENDDO
   ENDDO
!-----------------------------------------------------------------------
   ko = 1
!-----------------------------------------------------------------------
! extract diagonals and accumulate remaining matrix.
!-----------------------------------------------------------------------
   DO i = 1 , N
      SPAG_Loop_2_2: DO k = Ia(i) , Ia(i+1) - 1
         j = Ja(k)
         DO l = 1 , Idiag
            IF ( j-i==Ioff(l) ) THEN
               Diag(i,l) = A(k)
               CYCLE SPAG_Loop_2_2
            ENDIF
         ENDDO
!--------------- append element not in any diagonal to ao,jao,iao -----
         IF ( job2/=0 ) THEN
            Ao(ko) = A(k)
            Jao(ko) = j
            ko = ko + 1
         ENDIF
      ENDDO SPAG_Loop_2_2
      IF ( job2/=0 ) Ind(i+1) = ko
   ENDDO
   IF ( job2==0 ) RETURN
!     finish with iao
   Iao(1) = 1
   DO i = 2 , N + 1
      Iao(i) = Ind(i)
   ENDDO
!----------- end of csrdia ---------------------------------------------
!-----------------------------------------------------------------------
END SUBROUTINE csrdia
!*==diacsr.f90 processed by SPAG 8.04RA 15:47 21 Feb 2024
!!SPAG Open source Personal, Educational or Academic User Clemson University  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
!-----------------------------------------------------------------------
SUBROUTINE diacsr(N,Job,Idiag,Diag,Ndiag,Ioff,A,Ja,Ia)
   USE ISO_FORTRAN_ENV                 
   IMPLICIT NONE
!
! Dummy argument declarations rewritten by SPAG
!
   INTEGER , INTENT(IN) :: Idiag
   INTEGER , INTENT(IN) :: Ndiag
   INTEGER , INTENT(IN) :: N
   INTEGER , INTENT(IN) :: Job
   REAL(REAL64) , INTENT(IN) , DIMENSION(Ndiag,Idiag) :: Diag
   INTEGER , INTENT(IN) , DIMENSION(*) :: Ioff
   REAL(REAL64) , INTENT(OUT) , DIMENSION(*) :: A
   INTEGER , INTENT(OUT) , DIMENSION(*) :: Ja
   INTEGER , INTENT(OUT) , DIMENSION(*) :: Ia
!
! Local variable declarations rewritten by SPAG
!
   INTEGER :: i , j , jj , ko
   REAL(REAL64) :: t
!
! End of declarations rewritten by SPAG
!
!-----------------------------------------------------------------------
!    diagonal format     to     compressed sparse row
!-----------------------------------------------------------------------
! this subroutine extract the idiag most important diagonals from the
! input matrix a, ja, ia, i.e, those diagonals of the matrix which have
! the largest number of nonzero elements. If requested (see job),
! the rest of the matrix is put in a the output matrix ao, jao, iao
!-----------------------------------------------------------------------
! on entry:
!----------
! n	= integer. dimension of the matrix a.
! job	= integer. job indicator with the following meaning.
!         if (job .eq. 0) then check for each entry in diag
!         whether this entry is zero. If it is then do not include
!         in the output matrix. Note that the test is a test for
!         an exact arithmetic zero. Be sure that the zeros are
!         actual zeros in double precision otherwise this would not
!         work.
!
! idiag = integer equal to the number of diagonals to be extracted.
!         Note: on return idiag may be modified.
!
! diag  = real array of size (ndiag x idiag) containing the diagonals
!         of A on return.
!
! ndiag = integer equal to the first dimension of array diag.
!
! ioff  = integer array of length idiag, containing the offsets of the
!   	  diagonals to be extracted.
!
! on return:
!-----------
! a,
! ja,
! ia    = matrix stored in a, ja, ia, format
!
! Note:
! ----- the arrays a and ja should be of length n*idiag.
!
!-----------------------------------------------------------------------
   Ia(1) = 1
   ko = 1
   DO i = 1 , N
      DO jj = 1 , Idiag
         j = i + Ioff(jj)
         IF ( j>=1 .AND. j<=N ) THEN
            t = Diag(i,jj)
            IF ( Job/=0 .OR. t/=0.0D0 ) THEN
               A(ko) = t
               Ja(ko) = j
               ko = ko + 1
            ENDIF
         ENDIF
      ENDDO
      Ia(i+1) = ko
   ENDDO
!----------- end of diacsr ---------------------------------------------
!-----------------------------------------------------------------------
END SUBROUTINE diacsr
!*==bsrcsr.f90 processed by SPAG 8.04RA 15:47 21 Feb 2024
!!SPAG Open source Personal, Educational or Academic User Clemson University  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
!-----------------------------------------------------------------------
SUBROUTINE bsrcsr(Job,N,M,Na,A,Ja,Ia,Ao,Jao,Iao)
   USE ISO_FORTRAN_ENV                 
   IMPLICIT NONE
!
! Dummy argument declarations rewritten by SPAG
!
   INTEGER , INTENT(IN) :: N
   INTEGER , INTENT(IN) :: Na
   INTEGER , INTENT(IN) :: Job
   INTEGER , INTENT(IN) :: M
   REAL(REAL64) , INTENT(IN) , DIMENSION(Na,*) :: A
   INTEGER , INTENT(IN) , DIMENSION(*) :: Ja
   INTEGER , INTENT(IN) , DIMENSION(*) :: Ia
   REAL(REAL64) , INTENT(OUT) , DIMENSION(*) :: Ao
   INTEGER , INTENT(OUT) , DIMENSION(*) :: Jao
   INTEGER , INTENT(OUT) , DIMENSION(N+1) :: Iao
!
! Local variable declarations rewritten by SPAG
!
   INTEGER :: i , i1 , i2 , ii , ij , irow , j , jstart , k , krow
   LOGICAL :: val
!
! End of declarations rewritten by SPAG
!
!-----------------------------------------------------------------------
!             Block Sparse Row  to Compressed Sparse Row.
!-----------------------------------------------------------------------
! NOTE: ** meanings of parameters may have changed wrt earlier versions
! FORMAT DEFINITION HAS CHANGED WRT TO EARLIER VERSIONS...
!-----------------------------------------------------------------------
!
! converts a  matrix stored in block-reduced   a, ja, ia  format to the
! general  sparse row a,  ja, ia  format.  A matrix   that has  a block
! structure is a matrix whose entries are blocks  of the same size m
! (e.g.  3 x 3).   Then it is often preferred  to work with the reduced
! graph of the matrix. Instead of storing one element at a time one can
! store a whole block at a time.  In this storage scheme  an entry is a
! square array holding the m**2 elements of a block.
!
!-----------------------------------------------------------------------
! on entry:
!----------
! job   = if job.eq.0 on entry, values are not copied (pattern only)
!
! n	= the block row dimension of the matrix.
!
! m     = the dimension of each block. Thus, the actual row dimension
!         of A is n x m.
!
! na	= first dimension of array a as declared in calling program.
!         This should be .ge. m**2.
!
! a	= real array containing the real entries of the matrix. Recall
!         that each entry is in fact an m x m block. These entries
!         are stored column-wise in locations a(1:m*m,k) for each k-th
!         entry. See details below.
!
! ja	= integer array of length n. ja(k) contains the column index
!         of the leading element, i.e., the element (1,1) of the block
!         that is held in the column a(*,k) of the value array.
!
! ia    = integer array of length n+1. ia(i) points to the beginning
!         of block row number i in the arrays a and ja.
!
! on return:
!-----------
! ao, jao,
! iao   = matrix stored in compressed sparse row format. The number of
!         rows in the new matrix is n x m.
!
! Notes: THIS CODE IS NOT IN PLACE.
!
!-----------------------------------------------------------------------
! BSR FORMAT.
!----------
! Each row of A contains the m x m block matrix unpacked column-
! wise (this allows the user to declare the array a as a(m,m,*) on entry
! if desired). The block rows are stored in sequence just as for the
! compressed sparse row format.
!
!-----------------------------------------------------------------------
!     example  with m = 2:
!                                                       1  2 3
!    +-------|--------|--------+                       +-------+
!    | 1   2 |  0   0 |  3   4 |     Block             | x 0 x | 1
!    | 5   6 |  0   0 |  7   8 |     Representation:   | 0 x x | 2
!    +-------+--------+--------+                       | x 0 0 | 3
!    | 0   0 |  9  10 | 11  12 |                       +-------+
!    | 0   0 | 13  14 | 15  16 |
!    +-------+--------+--------+
!    | 17 18 |  0   0 |  0   0 |
!    | 22 23 |  0   0 |  0   0 |
!    +-------+--------+--------+
!
!    For this matrix:     n    = 3
!                         m    = 2
!                         nnz  = 5
!-----------------------------------------------------------------------
! Data structure in Block Sparse Row format:
!-------------------------------------------
! Array A:
!-------------------------
!     1   3   9   11   17   <<--each m x m block is stored column-wise
!     5   7   13  15   22       in a  column of the array A.
!     2   4   10  12   18
!     6   8   14  16   23
!-------------------------
! JA  1   3   2    3    1   <<-- column indices for each block. Note that
!-------------------------       these indices are wrt block matrix.
! IA  1   3   5    6        <<-- pointers to beginning of each block row
!-------------------------       in arrays A and JA.
!-----------------------------------------------------------------------
! locals
!
!
   val = (Job/=0)
   irow = 1
   krow = 1
   Iao(irow) = 1
!-----------------------------------------------------------------------
   DO ii = 1 , N
!
!     recall: n is the block-row dimension
!
      i1 = Ia(ii)
      i2 = Ia(ii+1) - 1
!
!     create m rows for each block row -- i.e., each k.
!
      DO i = 1 , M
         DO k = i1 , i2
            jstart = M*(Ja(k)-1)
            DO j = 1 , M
               ij = (j-1)*M + i
               IF ( val ) Ao(krow) = A(ij,k)
               Jao(krow) = jstart + j
               krow = krow + 1
            ENDDO
         ENDDO
         irow = irow + 1
         Iao(irow) = krow
      ENDDO
   ENDDO
!-------------end-of-bsrcsr --------------------------------------------
!-----------------------------------------------------------------------
END SUBROUTINE bsrcsr
!*==csrbsr.f90 processed by SPAG 8.04RA 15:47 21 Feb 2024
!!SPAG Open source Personal, Educational or Academic User Clemson University  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
!-----------------------------------------------------------------------
SUBROUTINE csrbsr(Job,Nrow,M,Na,A,Ja,Ia,Ao,Jao,Iao,Iw,Ierr)
   USE ISO_FORTRAN_ENV                 
   IMPLICIT NONE
!
! Dummy argument declarations rewritten by SPAG
!
   INTEGER , INTENT(IN) :: Nrow
   INTEGER , INTENT(IN) :: Na
   INTEGER , INTENT(IN) :: Job
   INTEGER , INTENT(IN) :: M
   REAL(REAL64) , INTENT(IN) , DIMENSION(*) :: A
   INTEGER , INTENT(IN) , DIMENSION(*) :: Ja
   INTEGER , INTENT(IN) , DIMENSION(Nrow+1) :: Ia
   REAL(REAL64) , INTENT(OUT) , DIMENSION(Na,*) :: Ao
   INTEGER , INTENT(INOUT) , DIMENSION(Na) :: Jao
   INTEGER , INTENT(INOUT) , DIMENSION(*) :: Iao
   INTEGER , INTENT(INOUT) , DIMENSION(*) :: Iw
   INTEGER , INTENT(INOUT) :: Ierr
!
! Local variable declarations rewritten by SPAG
!
   INTEGER :: i , ii , ij , io , irow , j , jpos , jr , k , ko , len , m2 , nr
   LOGICAL :: vals
!
! End of declarations rewritten by SPAG
!
!-----------------------------------------------------------------------
!     Compressed Sparse Row  to    Block Sparse Row
!-----------------------------------------------------------------------
!
! This  subroutine converts a matrix stored  in a general compressed a,
! ja, ia format into a a block  sparse row format a(m,m,*),ja(*),ia(*).
! See routine  bsrcsr  for more  details on  data   structure for block
! matrices.
!
! NOTES: 1) the initial matrix does not have to have a block structure.
! zero padding is done for general sparse matrices.
!        2) For most practical purposes, na should be the same as m*m.
!
!-----------------------------------------------------------------------
!
! In what follows nr=1+(nrow-1)/m = block-row dimension of output matrix
!
! on entry:
!----------
!
! job   =  job indicator.
!          job =  0 -> only the pattern of output matrix is generated
!          job >  0 -> both pattern and values are generated.
!          job = -1 -> iao(1) will return the number of nonzero blocks,
!            in the output matrix. In this case jao(1:nr) is used as
!            workspace, ao is untouched, iao is untouched except iao(1)
!
! nrow	= integer, the actual row dimension of the matrix.
!
! m     = integer equal to the dimension of each block. m should be > 0.
!
! na	= first dimension of array ao as declared in calling program.
!         na should be .ge. m*m.
!
! a, ja,
!    ia = input matrix stored in compressed sparse row format.
!
! on return:
!-----------
!
! ao    = real  array containing the  values of the matrix. For details
!         on the format  see below. Each  row of  a contains the  m x m
!         block matrix  unpacked column-wise (this  allows the  user to
!         declare the  array a as ao(m,m,*) on  entry if desired).  The
!         block rows are stored in sequence  just as for the compressed
!         sparse row format. The block  dimension  of the output matrix
!         is  nr = 1 + (nrow-1) / m.
!
! jao   = integer array. containing the block-column indices of the
!         block-matrix. Each jao(k) is an integer between 1 and nr
!         containing the block column index of the block ao(*,k).
!
! iao   = integer array of length nr+1. iao(i) points to the beginning
!         of block row number i in the arrays ao and jao. When job=-1
!         iao(1) contains the number of nonzero blocks of the output
!         matrix and the rest of iao is unused. This is useful for
!         determining the lengths of ao and jao.
!
! ierr  = integer, error code.
!              0 -- normal termination
!              1 -- m is equal to zero
!              2 -- NA too small to hold the blocks (should be .ge. m**2)
!
! Work arrays:
!-------------
! iw    = integer work array of dimension  nr = 1 + (nrow-1) / m
!
! NOTES:
!-------
!     1) this code is not in place.
!     2) see routine bsrcsr for details on data sctructure for block
!        sparse row format.
!
!-----------------------------------------------------------------------
!     nr is the block-dimension of the output matrix.
!
!-----
   Ierr = 0
   IF ( M*M>Na ) Ierr = 2
   IF ( M==0 ) Ierr = 1
   IF ( Ierr/=0 ) RETURN
!-----------------------------------------------------------------------
   vals = (Job>0)
   nr = 1 + (Nrow-1)/M
   m2 = M*M
   ko = 1
   io = 1
   Iao(io) = 1
   len = 0
!
!     iw determines structure of block-row (nonzero indicator)
!
   DO j = 1 , nr
      Iw(j) = 0
   ENDDO
!
!     big loop -- leap by m rows each time.
!
   DO ii = 1 , Nrow , M
      irow = 0
!
!     go through next m rows -- make sure not to go beyond nrow.
!
      DO WHILE ( ii+irow<=Nrow .AND. irow<=M-1 )
         DO k = Ia(ii+irow) , Ia(ii+irow+1) - 1
!
!     block column index = (scalar column index -1) / m + 1
!
            j = Ja(k) - 1
            jr = j/M + 1
            j = j - (jr-1)*M
            jpos = Iw(jr)
            IF ( jpos==0 ) THEN
!
!     create a new block
!
               Iw(jr) = ko
               Jao(ko) = jr
               IF ( vals ) THEN
!
!     initialize new block to zero -- then copy nonzero element
!
                  DO i = 1 , m2
                     Ao(i,ko) = 0.0D0
                  ENDDO
                  ij = j*M + irow + 1
                  Ao(ij,ko) = A(k)
               ENDIF
               ko = ko + 1
            ELSE
!
!     copy column index and nonzero element
!
               Jao(jpos) = jr
               ij = j*M + irow + 1
               IF ( vals ) Ao(ij,jpos) = A(k)
            ENDIF
         ENDDO
         irow = irow + 1
      ENDDO
!
!     refresh iw
!
      DO j = Iao(io) , ko - 1
         Iw(Jao(j)) = 0
      ENDDO
      IF ( Job==-1 ) THEN
         len = len + ko - 1
         ko = 1
      ELSE
         io = io + 1
         Iao(io) = ko
      ENDIF
   ENDDO
   IF ( Job==-1 ) Iao(1) = len
!
!--------------end-of-csrbsr--------------------------------------------
!-----------------------------------------------------------------------
END SUBROUTINE csrbsr
!*==csrbnd.f90 processed by SPAG 8.04RA 15:47 21 Feb 2024
!!SPAG Open source Personal, Educational or Academic User Clemson University  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
!-----------------------------------------------------------------------
SUBROUTINE csrbnd(N,A,Ja,Ia,Job,Abd,Nabd,Lowd,Ml,Mu,Ierr)
   USE ISO_FORTRAN_ENV                 
   IMPLICIT NONE
!
! Dummy argument declarations rewritten by SPAG
!
   INTEGER :: N
   INTEGER , INTENT(IN) :: Nabd
   REAL(REAL64) , INTENT(IN) , DIMENSION(*) :: A
   INTEGER , DIMENSION(*) :: Ja
   INTEGER , DIMENSION(N+1) :: Ia
   INTEGER , INTENT(IN) :: Job
   REAL(REAL64) , INTENT(OUT) , DIMENSION(Nabd,N) :: Abd
   INTEGER , INTENT(INOUT) :: Lowd
   INTEGER :: Ml
   INTEGER :: Mu
   INTEGER , INTENT(INOUT) :: Ierr
!
! Local variable declarations rewritten by SPAG
!
   INTEGER :: i , ii , j , k , m , mdiag
   EXTERNAL getbwd
!
! End of declarations rewritten by SPAG
!
!-----------------------------------------------------------------------
!   Compressed Sparse Row  to  Banded (Linpack ) format.
!-----------------------------------------------------------------------
! this subroutine converts a general sparse matrix stored in
! compressed sparse row format into the banded format. for the
! banded format,the Linpack conventions are assumed (see below).
!-----------------------------------------------------------------------
! on entry:
!----------
! n	= integer,the actual row dimension of the matrix.
!
! a,
! ja,
! ia    = input matrix stored in compressed sparse row format.
!
! job	= integer. if job=1 then the values of the lower bandwith ml
!         and the upper bandwidth mu are determined internally.
!         otherwise it is assumed that the values of ml and mu
!         are the correct bandwidths on input. See ml and mu below.
!
! nabd  = integer. first dimension of array abd.
!
! lowd  = integer. this should be set to the row number in abd where
!         the lowest diagonal (leftmost) of A is located.
!         lowd should be  ( 1  .le.  lowd  .le. nabd).
!         if it is not known in advance what lowd should be
!         enter lowd = 0 and the default value lowd = ml+mu+1
!         will be chosen. Alternative: call routine getbwd from unary
!         first to detrermione ml and mu then define lowd accordingly.
!         (Note: the banded solvers in linpack use lowd=2*ml+mu+1. )
!
! ml	= integer. equal to the bandwidth of the strict lower part of A
! mu	= integer. equal to the bandwidth of the strict upper part of A
!         thus the total bandwidth of A is ml+mu+1.
!         if ml+mu+1 is found to be larger than lowd then an error
!         flag is raised (unless lowd = 0). see ierr.
!
! note:   ml and mu are assumed to have	 the correct bandwidth values
!         as defined above if job is set to zero on entry.
!
! on return:
!-----------
!
! abd   = real array of dimension abd(nabd,n).
!         on return contains the values of the matrix stored in
!         banded form. The j-th column of abd contains the elements
!         of the j-th column of  the original matrix comprised in the
!         band ( i in (j-ml,j+mu) ) with the lowest diagonal at
!         the bottom row (row lowd). See details below for this format.
!
! ml	= integer. equal to the bandwidth of the strict lower part of A
! mu	= integer. equal to the bandwidth of the strict upper part of A
!         if job=1 on entry then these two values are internally computed.
!
! lowd  = integer. row number in abd where the lowest diagonal
!         (leftmost) of A is located on return. In case lowd = 0
!         on return, then it is defined to ml+mu+1 on return and the
!         lowd will contain this value on return. `
!
! ierr  = integer. used for error messages. On return:
!         ierr .eq. 0  :means normal return
!         ierr .eq. -1 : means invalid value for lowd. (either .lt. 0
!         or larger than nabd).
!         ierr .eq. -2 : means that lowd is not large enough and as
!         result the matrix cannot be stored in array abd.
!         lowd should be at least ml+mu+1, where ml and mu are as
!         provided on output.
!
!----------------------------------------------------------------------*
! Additional details on banded format.  (this closely follows the      *
! format used in linpack. may be useful for converting a matrix into   *
! this storage format in order to use the linpack  banded solvers).    *
!----------------------------------------------------------------------*
!             ---  band storage format  for matrix abd ---             *
! uses ml+mu+1 rows of abd(nabd,*) to store the diagonals of           *
! a in rows of abd starting from the lowest (sub)-diagonal  which  is  *
! stored in row number lowd of abd. the minimum number of rows needed  *
! in abd is ml+mu+1, i.e., the minimum value for lowd is ml+mu+1. the  *
! j-th  column  of  abd contains the elements of the j-th column of a, *
! from bottom to top: the element a(j+ml,j) is stored in  position     *
! abd(lowd,j), then a(j+ml-1,j) in position abd(lowd-1,j) and so on.   *
! Generally, the element a(j+k,j) of original matrix a is stored in    *
! position abd(lowd+k-ml,j), for k=ml,ml-1,..,0,-1, -mu.               *
! The first dimension nabd of abd must be .ge. lowd                    *
!                                                                      *
!     example [from linpack ]:   if the original matrix is             *
!                                                                      *
!              11 12 13  0  0  0                                       *
!              21 22 23 24  0  0                                       *
!               0 32 33 34 35  0     original banded matrix            *
!               0  0 43 44 45 46                                       *
!               0  0  0 54 55 56                                       *
!               0  0  0  0 65 66                                       *
!                                                                      *
! then  n = 6, ml = 1, mu = 2. lowd should be .ge. 4 (=ml+mu+1)  and   *
! if lowd = 5 for example, abd  should be:                             *
!                                                                      *
! untouched --> x  x  x  x  x  x                                       *
!               *  * 13 24 35 46                                       *
!               * 12 23 34 45 56    resulting abd matrix in banded     *
!              11 22 33 44 55 66    format                             *
!  row lowd--> 21 32 43 54 65  *                                       *
!                                                                      *
! * = not used                                                         *
!
!
!----------------------------------------------------------------------*
! first determine ml and mu.
!-----------------------------------------------------------------------
   Ierr = 0
!-----------
   IF ( Job==1 ) CALL getbwd(N,Ja,Ia,Ml,Mu)
   m = Ml + Mu + 1
   IF ( Lowd==0 ) Lowd = m
   IF ( m>Lowd ) Ierr = -2
   IF ( Lowd>Nabd .OR. Lowd<0 ) Ierr = -1
   IF ( Ierr<0 ) RETURN
!------------
   DO i = 1 , m
      ii = Lowd - i + 1
      DO j = 1 , N
         Abd(ii,j) = 0.0D0
      ENDDO
   ENDDO
!---------------------------------------------------------------------
   mdiag = Lowd - Ml
   DO i = 1 , N
      DO k = Ia(i) , Ia(i+1) - 1
         j = Ja(k)
         Abd(i-j+mdiag,j) = A(k)
      ENDDO
   ENDDO
!------------- end of csrbnd -------------------------------------------
!-----------------------------------------------------------------------
END SUBROUTINE csrbnd
!*==bndcsr.f90 processed by SPAG 8.04RA 15:47 21 Feb 2024
!!SPAG Open source Personal, Educational or Academic User Clemson University  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
!-----------------------------------------------------------------------
SUBROUTINE bndcsr(N,Abd,Nabd,Lowd,Ml,Mu,A,Ja,Ia,Len,Ierr)
   USE ISO_FORTRAN_ENV                 
   IMPLICIT NONE
!
! Dummy argument declarations rewritten by SPAG
!
   INTEGER , INTENT(IN) :: N
   INTEGER , INTENT(IN) :: Nabd
   REAL(REAL64) , INTENT(IN) , DIMENSION(Nabd,*) :: Abd
   INTEGER , INTENT(IN) :: Lowd
   INTEGER , INTENT(IN) :: Ml
   INTEGER , INTENT(IN) :: Mu
   REAL(REAL64) , INTENT(OUT) , DIMENSION(*) :: A
   INTEGER , INTENT(OUT) , DIMENSION(*) :: Ja
   INTEGER , INTENT(OUT) , DIMENSION(N+1) :: Ia
   INTEGER , INTENT(IN) :: Len
   INTEGER , INTENT(OUT) :: Ierr
!
! Local variable declarations rewritten by SPAG
!
   INTEGER :: i , irow , j , ko
   REAL(REAL64) :: t
!
! End of declarations rewritten by SPAG
!
   INTEGER :: spag_nextblock_1
!-----------------------------------------------------------------------
! Banded (Linpack ) format   to    Compressed Sparse Row  format.
!-----------------------------------------------------------------------
! on entry:
!----------
! n	= integer,the actual row dimension of the matrix.
!
! nabd  = first dimension of array abd.
!
! abd   = real array containing the values of the matrix stored in
!         banded form. The j-th column of abd contains the elements
!         of the j-th column of  the original matrix,comprised in the
!         band ( i in (j-ml,j+mu) ) with the lowest diagonal located
!         in row lowd (see below).
!
! lowd  = integer. this should be set to the row number in abd where
!         the lowest diagonal (leftmost) of A is located.
!         lowd should be s.t.  ( 1  .le.  lowd  .le. nabd).
!         The subroutines dgbco, ... of linpack use lowd=2*ml+mu+1.
!
! ml	= integer. equal to the bandwidth of the strict lower part of A
! mu	= integer. equal to the bandwidth of the strict upper part of A
!         thus the total bandwidth of A is ml+mu+1.
!         if ml+mu+1 is found to be larger than nabd then an error
!         message is set. see ierr.
!
! len   = integer. length of arrays a and ja. bndcsr will stop if the
!         length of the arrays a and ja is insufficient to store the
!         matrix. see ierr.
!
! on return:
!-----------
! a,
! ja,
! ia    = input matrix stored in compressed sparse row format.
!
! lowd  = if on entry lowd was zero then lowd is reset to the default
!         value ml+mu+l.
!
! ierr  = integer. used for error message output.
!         ierr .eq. 0 :means normal return
!         ierr .eq. -1 : means invalid value for lowd.
!	  ierr .gt. 0 : means that there was not enough storage in a and ja
!         for storing the ourput matrix. The process ran out of space
!         (as indicated by len) while trying to fill row number ierr.
!         This should give an idea of much more storage might be required.
!         Moreover, the first irow-1 rows are correctly filled.
!
! notes:  the values in abd found to be equal to zero
! -----   (actual test: if (abd(...) .eq. 0.0d0) are removed.
!         The resulting may not be identical to a csr matrix
!         originally transformed to a bnd format.
!
!-----------------------------------------------------------------------
   Ierr = 0
!-----------
   IF ( Lowd>Nabd .OR. Lowd<=0 ) THEN
      Ierr = -1
      RETURN
   ENDIF
!-----------
   ko = 1
   Ia(1) = 1
   DO irow = 1 , N
      spag_nextblock_1 = 1
      SPAG_DispatchLoop_1: DO
         SELECT CASE (spag_nextblock_1)
         CASE (1)
!-----------------------------------------------------------------------
            i = Lowd
            DO j = irow - Ml , irow + Mu
               IF ( j>0 ) THEN
                  IF ( j>N ) THEN
                     spag_nextblock_1 = 2
                     CYCLE SPAG_DispatchLoop_1
                  ENDIF
                  t = Abd(i,j)
                  IF ( t/=0.0D0 ) THEN
                     IF ( ko>Len ) THEN
                        Ierr = irow
                        RETURN
                     ENDIF
                     A(ko) = t
                     Ja(ko) = j
                     ko = ko + 1
                  ENDIF
               ENDIF
               i = i - 1
            ENDDO
            spag_nextblock_1 = 2
         CASE (2)
!     end for row irow
            Ia(irow+1) = ko
            EXIT SPAG_DispatchLoop_1
         END SELECT
      ENDDO SPAG_DispatchLoop_1
   ENDDO
!------------- end of bndcsr -------------------------------------------
!-----------------------------------------------------------------------
END SUBROUTINE bndcsr
!*==csrssk.f90 processed by SPAG 8.04RA 15:47 21 Feb 2024
!!SPAG Open source Personal, Educational or Academic User Clemson University  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
!-----------------------------------------------------------------------
SUBROUTINE csrssk(N,Imod,A,Ja,Ia,Asky,Isky,Nzmax,Ierr)
   USE ISO_FORTRAN_ENV                 
   IMPLICIT NONE
!
! Dummy argument declarations rewritten by SPAG
!
   INTEGER , INTENT(IN) :: N
   INTEGER , INTENT(IN) :: Nzmax
   INTEGER , INTENT(IN) :: Imod
   REAL(REAL64) , INTENT(IN) , DIMENSION(*) :: A
   INTEGER , INTENT(IN) , DIMENSION(*) :: Ja
   INTEGER , INTENT(IN) , DIMENSION(N+1) :: Ia
   REAL(REAL64) , INTENT(OUT) , DIMENSION(Nzmax) :: Asky
   INTEGER , INTENT(INOUT) , DIMENSION(N+1) :: Isky
   INTEGER , INTENT(OUT) :: Ierr
!
! Local variable declarations rewritten by SPAG
!
   INTEGER :: i , j , k , kend , ml , nnz
!
! End of declarations rewritten by SPAG
!
!-----------------------------------------------------------------------
!      Compressed Sparse Row         to     Symmetric Skyline Format
!  or  Symmetric Sparse Row
!-----------------------------------------------------------------------
! this subroutine translates a compressed sparse row or a symmetric
! sparse row format into a symmetric skyline format.
! the input matrix can be in either compressed sparse row or the
! symmetric sparse row format. The output matrix is in a symmetric
! skyline format: a real array containing the (active portions) of the
! rows in  sequence and a pointer to the beginning of each row.
!
! This module is NOT  in place.
!-----------------------------------------------------------------------
! Coded by Y. Saad, Oct 5, 1989. Revised Feb. 18, 1991.
!-----------------------------------------------------------------------
!
! on entry:
!----------
! n	= integer equal to the dimension of A.
! imod  = integer indicating the variant of skyline format wanted:
!         imod = 0 means the pointer isky points to the `zeroth'
!         element of the row, i.e., to the position of the diagonal
!         element of previous row (for i=1, isky(1)= 0)
!         imod = 1 means that itpr points to the beginning of the row.
!         imod = 2 means that isky points to the end of the row (diagonal
!                  element)
!
! a	= real array of size nna containing the nonzero elements
! ja	= integer array of size	nnz containing the column positions
! 	  of the corresponding elements in a.
! ia	= integer of size n+1. ia(k) contains the position in a, ja of
!	  the beginning of the k-th row.
! nzmax = integer. must be set to the number of available locations
!         in the output array asky.
!
! on return:
!----------
!
! asky    = real array containing the values of the matrix stored in skyline
!         format. asky contains the sequence of active rows from
!         i=1, to n, an active row being the row of elemnts of
!         the matrix contained between the leftmost nonzero element
!         and the diagonal element.
! isky	= integer array of size n+1 containing the pointer array to
!         each row. The meaning of isky depends on the input value of
!         imod (see above).
! ierr  =  integer.  Error message. If the length of the
!         output array asky exceeds nzmax. ierr returns the minimum value
!         needed for nzmax. otherwise ierr=0 (normal return).
!
! Notes:
!         1) This module is NOT  in place.
!         2) even when imod = 2, length of  isky is  n+1, not n.
!
!-----------------------------------------------------------------------
! first determine individial bandwidths and pointers.
!-----------------------------------------------------------------------
   Ierr = 0
   Isky(1) = 0
   DO i = 1 , N
      ml = 0
      DO k = Ia(i) , Ia(i+1) - 1
         ml = max(ml,i-Ja(k)+1)
      ENDDO
      Isky(i+1) = Isky(i) + ml
   ENDDO
!
!     test if there is enough space  asky to do the copying.
!
   nnz = Isky(N+1)
   IF ( nnz>Nzmax ) THEN
      Ierr = nnz
      RETURN
   ENDIF
!
!   fill asky with zeros.
!
   DO k = 1 , nnz
      Asky(k) = 0.0D0
   ENDDO
!
!     copy nonzero elements.
!
   DO i = 1 , N
      kend = Isky(i+1)
      DO k = Ia(i) , Ia(i+1) - 1
         j = Ja(k)
         IF ( j<=i ) Asky(kend+j-i) = A(k)
      ENDDO
   ENDDO
!
! modify pointer according to imod if necessary.
!
   IF ( Imod==0 ) RETURN
   IF ( Imod==1 ) THEN
      DO k = 1 , N + 1
         Isky(k) = Isky(k) + 1
      ENDDO
   ENDIF
   IF ( Imod==2 ) THEN
      DO k = 1 , N
         Isky(k) = Isky(k+1)
      ENDDO
   ENDIF
!
!------------- end of csrssk -------------------------------------------
!-----------------------------------------------------------------------
END SUBROUTINE csrssk
!*==sskssr.f90 processed by SPAG 8.04RA 15:47 21 Feb 2024
!!SPAG Open source Personal, Educational or Academic User Clemson University  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
!-----------------------------------------------------------------------
SUBROUTINE sskssr(N,Imod,Asky,Isky,Ao,Jao,Iao,Nzmax,Ierr)
   USE ISO_FORTRAN_ENV                 
   IMPLICIT NONE
!
! Dummy argument declarations rewritten by SPAG
!
   INTEGER , INTENT(IN) :: N
   INTEGER , INTENT(IN) :: Nzmax
   INTEGER , INTENT(IN) :: Imod
   REAL(REAL64) , INTENT(IN) , DIMENSION(*) :: Asky
   INTEGER , INTENT(IN) , DIMENSION(N+1) :: Isky
   REAL(REAL64) , INTENT(OUT) , DIMENSION(Nzmax) :: Ao
   INTEGER , INTENT(OUT) , DIMENSION(Nzmax) :: Jao
   INTEGER , INTENT(OUT) , DIMENSION(N+1) :: Iao
   INTEGER , INTENT(OUT) :: Ierr
!
! Local variable declarations rewritten by SPAG
!
   INTEGER :: i , j , k , kend , kstart , next
!
! End of declarations rewritten by SPAG
!
!-----------------------------------------------------------------------
!     Symmetric Skyline Format  to  Symmetric Sparse Row format.
!-----------------------------------------------------------------------
!  tests for exact zeros in skyline matrix (and ignores them in
!  output matrix).  In place routine (a, isky :: ao, iao)
!-----------------------------------------------------------------------
! this subroutine translates a  symmetric skyline format into a
! symmetric sparse row format. Each element is tested to see if it is
! a zero element. Only the actual nonzero elements are retained. Note
! that the test used is simple and does take into account the smallness
! of a value. the subroutine filter (see unary module) can be used
! for this purpose.
!-----------------------------------------------------------------------
! Coded by Y. Saad, Oct 5, 1989. Revised Feb 18, 1991./
!-----------------------------------------------------------------------
!
! on entry:
!----------
! n	= integer equal to the dimension of A.
! imod  = integer indicating the variant of skyline format used:
!         imod = 0 means the pointer iao points to the `zeroth'
!         element of the row, i.e., to the position of the diagonal
!         element of previous row (for i=1, iao(1)= 0)
!         imod = 1 means that itpr points to the beginning of the row.
!         imod = 2 means that iao points to the end of the row
!                  (diagonal element)
! asky  = real array containing the values of the matrix. asky contains
!         the sequence of active rows from i=1, to n, an active row
!         being the row of elemnts of the matrix contained between the
!         leftmost nonzero element and the diagonal element.
! isky 	= integer array of size n+1 containing the pointer array to
!         each row. isky (k) contains the address of the beginning of the
!         k-th active row in the array asky.
! nzmax = integer. equal to the number of available locations in the
!         output array ao.
!
! on return:
! ----------
! ao	= real array of size nna containing the nonzero elements
! jao	= integer array of size	nnz containing the column positions
! 	  of the corresponding elements in a.
! iao	= integer of size n+1. iao(k) contains the position in a, ja of
!	  the beginning of the k-th row.
! ierr  = integer. Serving as error message. If the length of the
!         output arrays ao, jao exceeds nzmax then ierr returns
!         the row number where the algorithm stopped: rows
!         i, to ierr-1 have been processed succesfully.
!         ierr = 0 means normal return.
!         ierr = -1  : illegal value for imod
! Notes:
!-------
! This module is in place: ao and iao can be the same as asky, and isky.
!-----------------------------------------------------------------------
! local variables
   Ierr = 0
!
! check for validity of imod
!
   IF ( Imod/=0 .AND. Imod/=1 .AND. Imod/=2 ) THEN
      Ierr = -1
      RETURN
   ENDIF
!
! next  = pointer to next available position in output matrix
! kend  = pointer to end of current row in skyline matrix.
!
   next = 1
!
! set kend = start position -1 in  skyline matrix.
!
   kend = 0
   IF ( Imod==1 ) kend = Isky(1) - 1
   IF ( Imod==0 ) kend = Isky(1)
!
! loop through all rows
!
   DO i = 1 , N
!
! save value of pointer to ith row in output matrix
!
      Iao(i) = next
!
! get beginnning and end of skyline  row
!
      kstart = kend + 1
      IF ( Imod==0 ) kend = Isky(i+1)
      IF ( Imod==1 ) kend = Isky(i+1) - 1
      IF ( Imod==2 ) kend = Isky(i)
!
! copy element into output matrix unless it is a zero element.
!
      DO k = kstart , kend
         IF ( Asky(k)/=0.0D0 ) THEN
            j = i - (kend-k)
            Jao(next) = j
            Ao(next) = Asky(k)
            next = next + 1
            IF ( next>Nzmax+1 ) THEN
               Ierr = i
               RETURN
            ENDIF
         ENDIF
      ENDDO
   ENDDO
   Iao(N+1) = next
!-------------end-of-sskssr --------------------------------------------
!-----------------------------------------------------------------------
END SUBROUTINE sskssr
!*==csrjad.f90 processed by SPAG 8.04RA 15:47 21 Feb 2024
!!SPAG Open source Personal, Educational or Academic User Clemson University  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
!-----------------------------------------------------------------------
SUBROUTINE csrjad(Nrow,A,Ja,Ia,Idiag,Iperm,Ao,Jao,Iao)
   USE ISO_FORTRAN_ENV                 
   IMPLICIT NONE
!
! Dummy argument declarations rewritten by SPAG
!
   INTEGER :: Nrow
   REAL(REAL64) , INTENT(IN) , DIMENSION(*) :: A
   INTEGER , INTENT(IN) , DIMENSION(*) :: Ja
   INTEGER , INTENT(IN) , DIMENSION(Nrow+1) :: Ia
   INTEGER , INTENT(INOUT) :: Idiag
   INTEGER , INTENT(INOUT) , DIMENSION(Nrow) :: Iperm
   REAL(REAL64) , INTENT(OUT) , DIMENSION(*) :: Ao
   INTEGER , INTENT(INOUT) , DIMENSION(*) :: Jao
   INTEGER , INTENT(INOUT) , DIMENSION(Nrow) :: Iao
!
! Local variable declarations rewritten by SPAG
!
   INTEGER :: i , ilo , j , jj , k , k0 , k1 , len
   EXTERNAL dcsort
!
! End of declarations rewritten by SPAG
!
!-----------------------------------------------------------------------
!    Compressed Sparse Row  to   JAgged Diagonal storage.
!-----------------------------------------------------------------------
! this subroutine converts  matrix stored in the compressed sparse
! row format to the jagged diagonal format. The data structure
! for the JAD (Jagged Diagonal storage) is as follows. The rows of
! the matrix are (implicitly) permuted so that their lengths are in
! decreasing order. The real entries ao(*) and their column indices
! jao(*) are stored in succession. The number of such diagonals is idiag.
! the lengths of each of these diagonals is stored in iao(*).
! For more details see [E. Anderson and Y. Saad,
! ``Solving sparse triangular systems on parallel computers'' in
! Inter. J. of High Speed Computing, Vol 1, pp. 73-96 (1989).]
! or  [Y. Saad, ``Krylov Subspace Methods on Supercomputers''
! SIAM J. on  Stat. Scient. Comput., volume 10, pp. 1200-1232 (1989).]
!-----------------------------------------------------------------------
! on entry:
!----------
! nrow 	  = row dimension of the matrix A.
!
! a,
! ia,
! ja      = input matrix in compressed sparse row format.
!
! on return:
!----------
!
! idiag = integer. The number of jagged diagonals in the matrix.
!
! iperm = integer array of length nrow containing the permutation
!         of the rows that leads to a decreasing order of the
!         number of nonzero elements.
!
! ao    = real array containing the values of the matrix A in
!         jagged diagonal storage. The j-diagonals are stored
!         in ao in sequence.
!
! jao   = integer array containing the column indices of the
!         entries in ao.
!
! iao   = integer array containing pointers to the beginning
!         of each j-diagonal in ao, jao. iao is also used as
!         a work array and it should be of length n at least.
!
!-----------------------------------------------------------------------
!     ---- define initial iperm and get lengths of each row
!     ---- jao is used a work vector to store tehse lengths
!
   Idiag = 0
   ilo = Nrow
   DO j = 1 , Nrow
      Iperm(j) = j
      len = Ia(j+1) - Ia(j)
      ilo = min(ilo,len)
      Idiag = max(Idiag,len)
      Jao(j) = len
   ENDDO
!
!     call sorter to get permutation. use iao as work array.
!
   CALL dcsort(Jao,Nrow,Iao,Iperm,ilo,Idiag)
!
!     define output data structure. first lengths of j-diagonals
!
   DO j = 1 , Nrow
      Iao(j) = 0
   ENDDO
   DO k = 1 , Nrow
      len = Jao(Iperm(k))
      DO i = 1 , len
         Iao(i) = Iao(i) + 1
      ENDDO
   ENDDO
!
!     get the output matrix itself
!
   k1 = 1
   k0 = k1
   DO jj = 1 , Idiag
      len = Iao(jj)
      DO k = 1 , len
         i = Ia(Iperm(k)) + jj - 1
         Ao(k1) = A(i)
         Jao(k1) = Ja(i)
         k1 = k1 + 1
      ENDDO
      Iao(jj) = k0
      k0 = k1
   ENDDO
   Iao(Idiag+1) = k1
!----------end-of-csrjad------------------------------------------------
!-----------------------------------------------------------------------
END SUBROUTINE csrjad
!*==jadcsr.f90 processed by SPAG 8.04RA 15:47 21 Feb 2024
!!SPAG Open source Personal, Educational or Academic User Clemson University  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
!-----------------------------------------------------------------------
SUBROUTINE jadcsr(Nrow,Idiag,A,Ja,Ia,Iperm,Ao,Jao,Iao)
   USE ISO_FORTRAN_ENV                 
   IMPLICIT NONE
!
! Dummy argument declarations rewritten by SPAG
!
   INTEGER , INTENT(IN) :: Nrow
   INTEGER , INTENT(IN) :: Idiag
   REAL(REAL64) , INTENT(IN) , DIMENSION(*) :: A
   INTEGER , INTENT(IN) , DIMENSION(*) :: Ja
   INTEGER , INTENT(IN) , DIMENSION(Idiag+1) :: Ia
   INTEGER , INTENT(IN) , DIMENSION(Nrow) :: Iperm
   REAL(REAL64) , INTENT(OUT) , DIMENSION(*) :: Ao
   INTEGER , INTENT(INOUT) , DIMENSION(*) :: Jao
   INTEGER , INTENT(INOUT) , DIMENSION(Nrow+1) :: Iao
!
! Local variable declarations rewritten by SPAG
!
   INTEGER :: i , j , jj , k , k1 , kpos , len
!
! End of declarations rewritten by SPAG
!
!-----------------------------------------------------------------------
!     Jagged Diagonal Storage   to     Compressed Sparse Row
!-----------------------------------------------------------------------
! this subroutine converts a matrix stored in the jagged diagonal format
! to the compressed sparse row format.
!-----------------------------------------------------------------------
! on entry:
!----------
! nrow 	  = integer. the row dimension of the matrix A.
!
! idiag   = integer. The  number of jagged diagonals in the data
!           structure a, ja, ia.
!
! a,
! ja,
! ia      = input matrix in jagged diagonal format.
!
! iperm   = permutation of the rows used to obtain the JAD ordering.
!
! on return:
!----------
!
! ao, jao,
! iao     = matrix in CSR format.
!-----------------------------------------------------------------------
! determine first the pointers for output matrix. Go through the
! structure once:
!
   DO j = 1 , Nrow
      Jao(j) = 0
   ENDDO
!
!     compute the lengths of each row of output matrix -
!
   DO i = 1 , Idiag
      len = Ia(i+1) - Ia(i)
      DO k = 1 , len
         Jao(Iperm(k)) = Jao(Iperm(k)) + 1
      ENDDO
   ENDDO
!
!     remember to permute
!
   kpos = 1
   Iao(1) = 1
   DO i = 1 , Nrow
      kpos = kpos + Jao(i)
      Iao(i+1) = kpos
   ENDDO
!
!     copy elemnts one at a time.
!
   DO jj = 1 , Idiag
      k1 = Ia(jj) - 1
      len = Ia(jj+1) - k1 - 1
      DO k = 1 , len
         kpos = Iao(Iperm(k))
         Ao(kpos) = A(k1+k)
         Jao(kpos) = Ja(k1+k)
         Iao(Iperm(k)) = kpos + 1
      ENDDO
   ENDDO
!
!     rewind pointers
!
   DO j = Nrow , 1 , -1
      Iao(j+1) = Iao(j)
   ENDDO
   Iao(1) = 1
!----------end-of-jadcsr------------------------------------------------
!-----------------------------------------------------------------------
END SUBROUTINE jadcsr
!*==dcsort.f90 processed by SPAG 8.04RA 15:47 21 Feb 2024
!!SPAG Open source Personal, Educational or Academic User Clemson University  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
SUBROUTINE dcsort(Ival,N,Icnt,Index,Ilo,Ihi)
   IMPLICIT NONE
!
! Dummy argument declarations rewritten by SPAG
!
   INTEGER , INTENT(IN) :: N
   INTEGER , INTENT(IN) :: Ilo
   INTEGER , INTENT(IN) :: Ihi
   INTEGER , INTENT(IN) , DIMENSION(N) :: Ival
   INTEGER , INTENT(INOUT) , DIMENSION(Ilo:Ihi) :: Icnt
   INTEGER , INTENT(OUT) , DIMENSION(N) :: Index
!
! Local variable declarations rewritten by SPAG
!
   INTEGER :: i , ivalj , j
!
! End of declarations rewritten by SPAG
!
!-----------------------------------------------------------------------
!     Specifications for arguments:
!     ----------------------------
!-----------------------------------------------------------------------
!    This routine computes a permutation which, when applied to the
!    input vector ival, sorts the integers in ival in descending
!    order.  The permutation is represented by the vector index.  The
!    permuted ival can be interpreted as follows:
!      ival(index(i-1)) .ge. ival(index(i)) .ge. ival(index(i+1))
!
!    A specialized sort, the distribution counting sort, is used
!    which takes advantage of the knowledge that
!        1)  The values are in the (small) range [ ilo, ihi ]
!        2)  Values are likely to be repeated often
!
!    contributed to SPARSKIT by Mike Heroux. (Cray Research)
!    ---------------------------------------
!-----------------------------------------------------------------------
! Usage:
!------
!     call dcsort( ival, n, icnt, index, ilo, ihi )
!
! Arguments:
!-----------
!    ival  integer array (input)
!          On entry, ia is an n dimensional array that contains
!          the values to be sorted.  ival is unchanged on exit.
!
!    n     integer (input)
!          On entry, n is the number of elements in ival and index.
!
!    icnt  integer (work)
!          On entry, is an integer work vector of length
!          (ihi - ilo + 1).
!
!    index integer array (output)
!          On exit, index is an n-length integer vector containing
!          the permutation which sorts the vector ival.
!
!    ilo   integer (input)
!          On entry, ilo is .le. to the minimum value in ival.
!
!    ihi   integer (input)
!          On entry, ihi is .ge. to the maximum value in ival.
!
! Remarks:
!---------
!         The permutation is NOT applied to the vector ival.
!
!----------------------------------------------------------------
!
!    Other integer values are temporary indices.
!
! Author:
!--------
!    Michael Heroux
!    Sandra Carney
!       Mathematical Software Research Group
!       Cray Research, Inc.
!
! References:
!    Knuth, Donald E., "The Art of Computer Programming, Volume 3:
!    Sorting and Searching," Addison-Wesley, Reading, Massachusetts,
!    1973, pp. 78-79.
!
! Revision history:
!    05/09/90: Original implementation.  A variation of the
!              Distribution Counting Sort recommended by
!              Sandra Carney. (Mike Heroux)
!
!-----------------------------------------------------------------
!     ----------------------------------
!     Specifications for local variables
!     ----------------------------------
!
!     --------------------------
!     First executable statement
!     --------------------------
   DO i = Ilo , Ihi
      Icnt(i) = 0
   ENDDO
!
   DO i = 1 , N
      Icnt(Ival(i)) = Icnt(Ival(i)) + 1
   ENDDO
!
   DO i = Ihi - 1 , Ilo , -1
      Icnt(i) = Icnt(i) + Icnt(i+1)
   ENDDO
!
   DO j = N , 1 , -1
      ivalj = Ival(j)
      Index(Icnt(ivalj)) = j
      Icnt(ivalj) = Icnt(ivalj) - 1
   ENDDO
END SUBROUTINE dcsort
!*==cooell.f90 processed by SPAG 8.04RA 15:47 21 Feb 2024
!!SPAG Open source Personal, Educational or Academic User Clemson University  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
!-------end-of-dcsort---------------------------------------------------
!-----------------------------------------------------------------------
SUBROUTINE cooell(Job,N,Nnz,A,Ja,Ia,Ao,Jao,Lda,Ncmax,Nc,Ierr)
   USE ISO_FORTRAN_ENV                 
   IMPLICIT NONE
!
! PARAMETER definitions rewritten by SPAG
!
   REAL(REAL64) , PARAMETER :: ZERO = 0.0D0
!
! Dummy argument declarations rewritten by SPAG
!
   INTEGER , INTENT(IN) :: Nnz
   INTEGER , INTENT(IN) :: Lda
   INTEGER , INTENT(IN) :: Ncmax
   INTEGER , INTENT(IN) :: Job
   INTEGER , INTENT(IN) :: N
   REAL(REAL64) , INTENT(IN) , DIMENSION(Nnz) :: A
   INTEGER , INTENT(IN) , DIMENSION(Nnz) :: Ja
   INTEGER , INTENT(IN) , DIMENSION(Nnz) :: Ia
   REAL(REAL64) , INTENT(OUT) , DIMENSION(Lda,Ncmax) :: Ao
   INTEGER , INTENT(INOUT) , DIMENSION(Lda,Ncmax) :: Jao
   INTEGER , INTENT(INOUT) :: Nc
   INTEGER , INTENT(OUT) :: Ierr
!
! Local variable declarations rewritten by SPAG
!
   LOGICAL :: copyval
   INTEGER :: i , ip , j , k
!
! End of declarations rewritten by SPAG
!
!-----------------------------------------------------------------------
!     COOrdinate format to ELLpack format
!-----------------------------------------------------------------------
!     On entry:
!     job     -- 0 if only pattern is to be processed(AO is not touched)
!     n       -- number of rows in the matrix
!     a,ja,ia -- input matix in COO format
!     lda     -- leading dimension of array AO and JAO
!     ncmax   -- size of the second dimension of array AO and JAO
!
!     On exit:
!     ao,jao  -- the matrix in ELL format
!     nc      -- maximum number of nonzeros per row
!     ierr    -- 0 if convertion succeeded
!                -1 if LDA < N
!                nc if NC > ncmax
!
!     NOTE: the last column of JAO is used as work space!!
!-----------------------------------------------------------------------
!     .. first executable statement ..
   copyval = (Job/=0)
   IF ( Lda<N ) THEN
      Ierr = -1
      RETURN
   ENDIF
!     .. use the last column of JAO as workspace
!     .. initialize the work space
   DO i = 1 , N
      Jao(i,Ncmax) = 0
   ENDDO
   Nc = 0
!     .. go through ia and ja to find out number nonzero per row
   DO k = 1 , Nnz
      i = Ia(k)
      Jao(i,Ncmax) = Jao(i,Ncmax) + 1
   ENDDO
!     .. maximum number of nonzero per row
   Nc = 0
   DO i = 1 , N
      IF ( Nc<Jao(i,Ncmax) ) Nc = Jao(i,Ncmax)
      Jao(i,Ncmax) = 0
   ENDDO
!     .. if nc > ncmax retrun now
   IF ( Nc>Ncmax ) THEN
      Ierr = Nc
      RETURN
   ENDIF
!     .. go through ia and ja to copy the matrix to AO and JAO
   DO k = 1 , Nnz
      i = Ia(k)
      j = Ja(k)
      Jao(i,Ncmax) = Jao(i,Ncmax) + 1
      ip = Jao(i,Ncmax)
      IF ( ip>Nc ) Nc = ip
      IF ( copyval ) Ao(i,ip) = A(k)
      Jao(i,ip) = j
   ENDDO
!     .. fill the unspecified elements of AO and JAO with zero diagonals
   DO i = 1 , N
      DO j = Ia(i+1) - Ia(i) + 1 , Nc
         Jao(i,j) = i
         IF ( copyval ) Ao(i,j) = ZERO
      ENDDO
   ENDDO
   Ierr = 0
!
END SUBROUTINE cooell
!*==xcooell.f90 processed by SPAG 8.04RA 15:47 21 Feb 2024
!!SPAG Open source Personal, Educational or Academic User Clemson University  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
!-----end-of-cooell-----------------------------------------------------
!-----------------------------------------------------------------------
SUBROUTINE xcooell(N,Nnz,A,Ja,Ia,Ac,Jac,Nac,Ner,Ncmax,Ierr)
   USE ISO_FORTRAN_ENV                 
   IMPLICIT NONE
!
! Dummy argument declarations rewritten by SPAG
!
   INTEGER , INTENT(IN) :: Nnz
   INTEGER , INTENT(IN) :: Nac
   INTEGER , INTENT(IN) :: Ner
   INTEGER , INTENT(IN) :: N
   REAL(REAL64) , INTENT(IN) , DIMENSION(Nnz) :: A
   INTEGER , INTENT(IN) , DIMENSION(Nnz) :: Ja
   INTEGER , INTENT(IN) , DIMENSION(Nnz) :: Ia
   REAL(REAL64) , INTENT(INOUT) , DIMENSION(Nac,Ner) :: Ac
   INTEGER , INTENT(OUT) , DIMENSION(Nac,Ner) :: Jac
   INTEGER , INTENT(INOUT) :: Ncmax
   INTEGER , INTENT(OUT) :: Ierr
!
! Local variable declarations rewritten by SPAG
!
   INTEGER :: icount , ii , in , inn , innz , is , k
!
! End of declarations rewritten by SPAG
!
!-----------------------------------------------------------------------
!   coordinate format to ellpack format.
!-----------------------------------------------------------------------
!
!   DATE WRITTEN: June 4, 1989.
!
!   PURPOSE
!   -------
!  This subroutine takes a sparse matrix in coordinate format and
!  converts it into the Ellpack-Itpack storage.
!
!  Example:
!  -------
!       (   11   0   13    0     0     0  )
!       |   21  22    0   24     0     0  |
!       |    0  32   33    0    35     0  |
!   A = |    0   0   43   44     0    46  |
!       |   51   0    0   54    55     0  |
!       (   61  62    0    0    65    66  )
!
!   Coordinate storage scheme:
!
!    A  = (11,22,33,44,55,66,13,21,24,32,35,43,46,51,54,61,62,65)
!    IA = (1, 2, 3, 4, 5, 6, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6, 6 )
!    JA = ( 1, 2, 3, 4, 5, 6, 3, 1, 4, 2, 5, 3, 6, 1, 4, 1, 2, 5)
!
!   Ellpack-Itpack storage scheme:
!
!       (   11  13    0    0   )          (   1   3   *    *  )
!       |   22  21   24    0   |          |   2   1   4    *  |
!  AC = |   33  32   35    0   |    JAC = |   3   2   5    *  |
!       |   44  43   46    0   |          |   4   3   6    *  |
!       |   55  51   54    0   |          |   5   1   4    *  |
!       (   66  61   62   65   )          (   6   1   2    5  )
!
!   Note: * means that you can store values from 1 to 6 (1 to n, where
!         n is the order of the matrix) in that position in the array.
!
!   Contributed by:
!   ---------------
!   Ernest E. Rothman
!   Cornell Thoery Center/Cornell National Supercomputer Facility
!   e-mail address: BITNET:   EER@CORNELLF.BITNET
!                   INTERNET: eer@cornellf.tn.cornell.edu
!
!   checked and modified  04/13/90 Y.Saad.
!
!   REFERENCES
!   ----------
!   Kincaid, D. R.; Oppe, T. C.; Respess, J. R.; Young, D. M. 1984.
!   ITPACKV 2C User's Guide, CNA-191. Center for Numerical Analysis,
!   University of Texas at Austin.
!
!   "Engineering and Scientific Subroutine Library; Guide and
!   Reference; Release 3 (SC23-0184-3). Pp. 79-86.
!
!-----------------------------------------------------------------------
!
!   INPUT PARAMETERS
!   ----------------
!  N       - Integer. The size of the square matrix.
!
!  NNZ     - Integer. Must be greater than or equal to the number of
!            nonzero elements in the sparse matrix. Dimension of A, IA
!            and JA.
!
!  NCA     - Integer. First dimension of output arrays ca and jac.
!
!  A(NNZ)  - Real array. (Double precision)
!            Stored entries of the sparse matrix A.
!            NNZ is the number of nonzeros.
!
!  IA(NNZ) - Integer array.
!            Pointers to specify rows for the stored nonzero entries
!            in A.
!
!  JA(NNZ) - Integer array.
!            Pointers to specify columns for the stored nonzero
!            entries in A.
!
!  NER     - Integer. Must be set greater than or equal to the maximum
!            number of nonzeros in any row of the sparse matrix.
!
!  OUTPUT PARAMETERS
!  -----------------
!  AC(NAC,*)  - Real array. (Double precision)
!               Stored entries of the sparse matrix A in compressed
!               storage mode.
!
!  JAC(NAC,*) - Integer array.
!               Contains the column numbers of the sparse matrix
!               elements stored in the corresponding positions in
!               array AC.
!
!  NCMAX   -  Integer. Equals the maximum number of nonzeros in any
!             row of the sparse matrix.
!
!  IERR    - Error parameter is returned as zero on successful
!             execution of the subroutin<e.
!             Error diagnostics are given by means of positive values
!             of this parameter as follows:
!
!             IERR = -1   -  NER is too small and should be set equal
!                            to NCMAX. The array AC may not be large
!                            enough to accomodate all the non-zeros of
!                            of the sparse matrix.
!             IERR =  1   -  The array AC has a zero column. (Warning)
!             IERR =  2   -  The array AC has a zero row.    (Warning)
!
!---------------------------------------------------------------------
!
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!   Initial error parameter to zero:
!
   Ierr = 0
!
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!   Initial output arrays to zero:
!
   DO in = 1 , Ner
      DO innz = 1 , N
         Jac(innz,in) = N
         Ac(innz,in) = 0.0D0
      ENDDO
   ENDDO
!
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
!   Assign nonzero elements of the sparse matrix (stored in the one
!   dimensional array A to the two dimensional array AC.
!   Also, assign the correct values with information about their
!   column indices to the two dimensional array KA. And at the same
!   time count the number of nonzeros in each row so that the
!   parameter NCMAX equals the maximum number of nonzeros in any row
!   of the sparse matrix.
!
   Ncmax = 1
   DO is = 1 , N
      k = 0
      DO ii = 1 , Nnz
         IF ( Ia(ii)==is ) THEN
            k = k + 1
            IF ( k<=Ner ) THEN
               Ac(is,k) = A(ii)
               Jac(is,k) = Ja(ii)
            ENDIF
         ENDIF
      ENDDO
      IF ( k>=Ncmax ) Ncmax = k
   ENDDO
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
!     Perform some simple error checks:
!
!heck maximum number of nonzeros in each row:
   IF ( Ncmax==Ner ) Ierr = 0
   IF ( Ncmax>Ner ) THEN
      Ierr = -1
      RETURN
   ENDIF
!
!heck if there are any zero columns in AC:
!
   DO in = 1 , Ncmax
      icount = 0
      DO inn = 1 , N
         IF ( Ac(inn,in)/=0.0D0 ) icount = 1
      ENDDO
      IF ( icount==0 ) THEN
         Ierr = 1
         RETURN
      ENDIF
   ENDDO
!
!heck if there are any zero rows in AC:
!
   DO inn = 1 , N
      icount = 0
      DO in = 1 , Ncmax
         IF ( Ac(inn,in)/=0.0D0 ) icount = 1
      ENDDO
      IF ( icount==0 ) THEN
         Ierr = 2
         RETURN
      ENDIF
   ENDDO
!------------- end of xcooell -------------------------------------------
END SUBROUTINE xcooell
!*==csruss.f90 processed by SPAG 8.04RA 15:47 21 Feb 2024
!!SPAG Open source Personal, Educational or Academic User Clemson University  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
!-----------------------------------------------------------------------
SUBROUTINE csruss(Nrow,A,Ja,Ia,Diag,Al,Jal,Ial,Au,Jau,Iau)
   USE ISO_FORTRAN_ENV                 
   IMPLICIT NONE
!
! Dummy argument declarations rewritten by SPAG
!
   INTEGER , INTENT(IN) :: Nrow
   REAL(REAL64) , INTENT(IN) , DIMENSION(*) :: A
   INTEGER , INTENT(IN) , DIMENSION(*) :: Ja
   INTEGER , INTENT(IN) , DIMENSION(Nrow+1) :: Ia
   REAL(REAL64) , INTENT(OUT) , DIMENSION(*) :: Diag
   REAL(REAL64) , INTENT(OUT) , DIMENSION(*) :: Al
   INTEGER , INTENT(OUT) , DIMENSION(*) :: Jal
   INTEGER , INTENT(INOUT) , DIMENSION(Nrow+1) :: Ial
   REAL(REAL64) , INTENT(OUT) , DIMENSION(*) :: Au
   INTEGER , INTENT(OUT) , DIMENSION(*) :: Jau
   INTEGER , INTENT(INOUT) , DIMENSION(Nrow+1) :: Iau
!
! Local variable declarations rewritten by SPAG
!
   INTEGER :: i , j , k , kl , ku
!
! End of declarations rewritten by SPAG
!
!-----------------------------------------------------------------------
! Compressed Sparse Row     to     Unsymmetric Sparse Skyline format
!-----------------------------------------------------------------------
! this subroutine converts a matrix stored in csr format into a nonsym.
! sparse skyline format. This latter format does not assume
! that the matrix has a symmetric pattern and consists of the following
! * the diagonal of A stored separately in diag(*);
! * The strict lower part of A is stored  in CSR format in al,jal,ial
! * The strict upper part is stored in CSC format in au,jau,iau.
!-----------------------------------------------------------------------
! On entry
!---------
! nrow  = dimension of the matrix a.
! a     = real array containing the nonzero values of the matrix
!         stored rowwise.
! ja    = column indices of the values in array a
! ia    = integer array of length n+1 containing the pointers to
!         beginning of each row in arrays a, ja.
!
! On return
!----------
! diag  = array containing the diagonal entries of A
! al,jal,ial = matrix in CSR format storing the strict lower
!              trangular part of A.
! au,jau,iau = matrix in CSC format storing the strict upper
!              triangular part of A.
!-----------------------------------------------------------------------
!
! determine U's data structure first
!
   DO i = 1 , Nrow + 1
      Iau(i) = 0
   ENDDO
   DO i = 1 , Nrow
      DO k = Ia(i) , Ia(i+1) - 1
         j = Ja(k)
         IF ( j>i ) Iau(j+1) = Iau(j+1) + 1
      ENDDO
   ENDDO
!
!     compute pointers from lengths
!
   Iau(1) = 1
   DO i = 1 , Nrow
      Iau(i+1) = Iau(i) + Iau(i+1)
      Ial(i+1) = Ial(i) + Ial(i+1)
   ENDDO
!
!     now do the extractions. scan all rows.
!
   kl = 1
   Ial(1) = kl
   DO i = 1 , Nrow
!
!     scan all elements in a row
!
      DO k = Ia(i) , Ia(i+1) - 1
         j = Ja(k)
!
!     if in upper part, store in row j (of transp(U) )
!
         IF ( j>i ) THEN
            ku = Iau(j)
            Au(ku) = A(k)
            Jau(ku) = i
            Iau(j) = ku + 1
         ELSEIF ( j==i ) THEN
            Diag(i) = A(k)
         ELSEIF ( j<i ) THEN
            Al(kl) = A(k)
            Jal(kl) = j
            kl = kl + 1
         ENDIF
      ENDDO
      Ial(i+1) = kl
   ENDDO
!
! readjust iau
!
   DO i = Nrow , 1 , -1
      Iau(i+1) = Iau(i)
   ENDDO
   Iau(1) = 1
!--------------- end-of-csruss -----------------------------------------
!-----------------------------------------------------------------------
END SUBROUTINE csruss
!*==usscsr.f90 processed by SPAG 8.04RA 15:47 21 Feb 2024
!!SPAG Open source Personal, Educational or Academic User Clemson University  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
!-----------------------------------------------------------------------
SUBROUTINE usscsr(Nrow,A,Ja,Ia,Diag,Al,Jal,Ial,Au,Jau,Iau)
   USE ISO_FORTRAN_ENV                 
   IMPLICIT NONE
!
! Dummy argument declarations rewritten by SPAG
!
   INTEGER , INTENT(IN) :: Nrow
   REAL(REAL64) , INTENT(OUT) , DIMENSION(*) :: A
   INTEGER , INTENT(OUT) , DIMENSION(*) :: Ja
   INTEGER , INTENT(INOUT) , DIMENSION(Nrow+1) :: Ia
   REAL(REAL64) , INTENT(IN) , DIMENSION(*) :: Diag
   REAL(REAL64) , INTENT(IN) , DIMENSION(*) :: Al
   INTEGER , INTENT(IN) , DIMENSION(*) :: Jal
   INTEGER , INTENT(IN) , DIMENSION(Nrow+1) :: Ial
   REAL(REAL64) , INTENT(IN) , DIMENSION(*) :: Au
   INTEGER , INTENT(IN) , DIMENSION(*) :: Jau
   INTEGER , INTENT(IN) , DIMENSION(Nrow+1) :: Iau
!
! Local variable declarations rewritten by SPAG
!
   INTEGER :: i , j , jak , k , ka
!
! End of declarations rewritten by SPAG
!
!-----------------------------------------------------------------------
! Unsymmetric Sparse Skyline   format   to Compressed Sparse Row
!-----------------------------------------------------------------------
! this subroutine converts a matrix stored in nonsymmetric sparse
! skyline format into csr format. The sparse skyline format is
! described in routine csruss.
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
! On entry
!-----------------------------------------------------------------------
! nrow  = dimension of the matrix a.
! diag  = array containing the diagonal entries of A
! al,jal,ial = matrix in CSR format storing the strict lower
!              trangular part of A.
! au,jau,iau = matrix in CSC format storing the strict upper
!              trangular part of A.
! On return
! ---------
! a     = real array containing the nonzero values of the matrix
!         stored rowwise.
! ja    = column indices of the values in array a
! ia    = integer array of length n+1 containing the pointers to
!         beginning of each row in arrays a, ja.
!
!-----------------------------------------------------------------------
!
! count elements in lower part + diagonal
!
   DO i = 1 , Nrow
      Ia(i+1) = Ial(i+1) - Ial(i) + 1
   ENDDO
!
! count elements in upper part
!
   DO i = 1 , Nrow
      DO k = Iau(i) , Iau(i+1) - 1
         j = Jau(k)
         Ia(j+1) = Ia(j+1) + 1
      ENDDO
   ENDDO
!---------- compute pointers from lengths ------------------------------
   Ia(1) = 1
   DO i = 1 , Nrow
      Ia(i+1) = Ia(i) + Ia(i+1)
   ENDDO
!
! copy lower part + diagonal
!
   DO i = 1 , Nrow
      ka = Ia(i)
      DO k = Ial(i) , Ial(i+1) - 1
         A(ka) = Al(k)
         Ja(ka) = Jal(k)
         ka = ka + 1
      ENDDO
      A(ka) = Diag(i)
      Ja(ka) = i
      Ia(i) = ka + 1
   ENDDO
!
!     copy upper part
!
   DO i = 1 , Nrow
      DO k = Iau(i) , Iau(i+1) - 1
!
! row number
!
         jak = Jau(k)
!
! where element goes
!
         ka = Ia(jak)
         A(ka) = Au(k)
         Ja(ka) = i
         Ia(jak) = ka + 1
      ENDDO
   ENDDO
!
! readjust ia
!
   DO i = Nrow , 1 , -1
      Ia(i+1) = Ia(i)
   ENDDO
   Ia(1) = 1
!----------end-of-usscsr------------------------------------------------
END SUBROUTINE usscsr
!*==csrsss.f90 processed by SPAG 8.04RA 15:47 21 Feb 2024
!!SPAG Open source Personal, Educational or Academic User Clemson University  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
!-----------------------------------------------------------------------
SUBROUTINE csrsss(Nrow,A,Ja,Ia,Sorted,Diag,Al,Jal,Ial,Au)
   USE ISO_FORTRAN_ENV                 
   IMPLICIT NONE
!
! Dummy argument declarations rewritten by SPAG
!
   INTEGER :: Nrow
   REAL(REAL64) , INTENT(IN) , DIMENSION(*) :: A
   INTEGER , INTENT(IN) , DIMENSION(*) :: Ja
   INTEGER , INTENT(IN) , DIMENSION(Nrow+1) :: Ia
   LOGICAL , INTENT(IN) :: Sorted
   REAL(REAL64) , INTENT(OUT) , DIMENSION(*) :: Diag
   REAL(REAL64) , DIMENSION(*) :: Al
   INTEGER , DIMENSION(*) :: Jal
   INTEGER , INTENT(INOUT) , DIMENSION(Nrow+1) :: Ial
   REAL(REAL64) , INTENT(OUT) , DIMENSION(*) :: Au
!
! Local variable declarations rewritten by SPAG
!
   INTEGER :: i , jak , k , kl , ku
   EXTERNAL csort
!
! End of declarations rewritten by SPAG
!
!-----------------------------------------------------------------------
! Compressed Sparse Row     to     Symmetric Sparse Skyline   format
!-----------------------------------------------------------------------
! this subroutine converts a matrix stored in csr format into the
! Symmetric sparse skyline   format. This latter format assumes that
! that the matrix has a symmetric pattern. It consists of the following
! * the diagonal of A stored separately in diag(*);
! * The strict lower part of A is stored  in csr format in al,jal,ial
! * The values only of strict upper part as stored in csc format in au.
!-----------------------------------------------------------------------
! On entry
!-----------
! nrow  = dimension of the matrix a.
! a     = real array containing the nonzero values of the matrix
!         stored rowwise.
! ja    = column indices of the values in array a
! ia    = integer array of length n+1 containing the pointers to
!         beginning of each row in arrays a, ja.
! sorted= a logical indicating whether or not the elements in a,ja,ia
!         are sorted.
!
! On return
! ---------
! diag  = array containing the diagonal entries of A
! al,jal,ial = matrix in csr format storing the strict lower
!              trangular part of A.
! au    = values of the strict upper trangular part of A, column wise.
!-----------------------------------------------------------------------
!
!     extract lower part and diagonal.
!
   kl = 1
   Ial(1) = kl
   DO i = 1 , Nrow
!
! scan all elements in a row
!
      DO k = Ia(i) , Ia(i+1) - 1
         jak = Ja(k)
         IF ( jak==i ) THEN
            Diag(i) = A(k)
         ELSEIF ( jak<i ) THEN
            Al(kl) = A(k)
            Jal(kl) = jak
            kl = kl + 1
         ENDIF
      ENDDO
      Ial(i+1) = kl
   ENDDO
!
! sort if not sorted
!
!%%%%%---- incompatible arg list!
   IF ( .NOT.Sorted ) CALL csort(Nrow,Al,Jal,Ial,.TRUE.)
!
! copy u
!
   DO i = 1 , Nrow
!
! scan all elements in a row
!
      DO k = Ia(i) , Ia(i+1) - 1
         jak = Ja(k)
         IF ( jak>i ) THEN
            ku = Ial(jak)
            Au(ku) = A(k)
            Ial(jak) = ku + 1
         ENDIF
      ENDDO
   ENDDO
!
! readjust ial
!
   DO i = Nrow , 1 , -1
      Ial(i+1) = Ial(i)
   ENDDO
   Ial(1) = 1
!--------------- end-of-csrsss -----------------------------------------
!-----------------------------------------------------------------------
END SUBROUTINE csrsss
!*==ssscsr.f90 processed by SPAG 8.04RA 15:47 21 Feb 2024
!!SPAG Open source Personal, Educational or Academic User Clemson University  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
!
SUBROUTINE ssscsr(Nrow,A,Ja,Ia,Diag,Al,Jal,Ial,Au)
   USE ISO_FORTRAN_ENV                 
   IMPLICIT NONE
!
! Dummy argument declarations rewritten by SPAG
!
   INTEGER , INTENT(IN) :: Nrow
   REAL(REAL64) , INTENT(OUT) , DIMENSION(*) :: A
   INTEGER , INTENT(OUT) , DIMENSION(*) :: Ja
   INTEGER , INTENT(INOUT) , DIMENSION(Nrow+1) :: Ia
   REAL(REAL64) , INTENT(IN) , DIMENSION(*) :: Diag
   REAL(REAL64) , INTENT(IN) , DIMENSION(*) :: Al
   INTEGER , INTENT(IN) , DIMENSION(*) :: Jal
   INTEGER , INTENT(IN) , DIMENSION(Nrow+1) :: Ial
   REAL(REAL64) , INTENT(IN) , DIMENSION(*) :: Au
!
! Local variable declarations rewritten by SPAG
!
   INTEGER :: i , j , jak , k , ka
!
! End of declarations rewritten by SPAG
!
!-----------------------------------------------------------------------
! Unsymmetric Sparse Skyline   format   to Compressed Sparse Row
!-----------------------------------------------------------------------
! this subroutine converts a matrix stored in nonsymmetric sparse
! skyline format into csr format. The sparse skyline format is
! described in routine csruss.
!-----------------------------------------------------------------------
! On entry
!---------
! diag  = array containing the diagonal entries of A
! al,jal,ial = matrix in csr format storing the strict lower
!              trangular part of A.
! au    = values of strict upper part.
!
! On return
! ---------
! nrow  = dimension of the matrix a.
! a     = real array containing the nonzero values of the matrix
!         stored rowwise.
! ja    = column indices of the values in array a
! ia    = integer array of length n+1 containing the pointers to
!         beginning of each row in arrays a, ja.
!
!-----------------------------------------------------------------------
!
! count elements in lower part + diagonal
!
   DO i = 1 , Nrow
      Ia(i+1) = Ial(i+1) - Ial(i) + 1
   ENDDO
!
! count elements in upper part
!
   DO i = 1 , Nrow
      DO k = Ial(i) , Ial(i+1) - 1
         j = Jal(k)
         Ia(j+1) = Ia(j+1) + 1
      ENDDO
   ENDDO
!---------- compute pointers from lengths ------------------------------
   Ia(1) = 1
   DO i = 1 , Nrow
      Ia(i+1) = Ia(i) + Ia(i+1)
   ENDDO
!
! copy lower part + diagonal
!
   DO i = 1 , Nrow
      ka = Ia(i)
      DO k = Ial(i) , Ial(i+1) - 1
         A(ka) = Al(k)
         Ja(ka) = Jal(k)
         ka = ka + 1
      ENDDO
      A(ka) = Diag(i)
      Ia(i) = ka + 1
   ENDDO
!
!     copy upper part
!
   DO i = 1 , Nrow
      DO k = Ial(i) , Ial(i+1) - 1
!
! row number
!
         jak = Jal(k)
!
! where element goes
!
         ka = Ia(jak)
         A(ka) = Au(k)
         Ja(ka) = i
         Ia(jak) = ka + 1
      ENDDO
   ENDDO
!
! readjust ia
!
   DO i = Nrow , 1 , -1
      Ia(i+1) = Ia(i)
   ENDDO
   Ia(1) = 1
!----------end-of-ssscsr------------------------------------------------
END SUBROUTINE ssscsr
!*==csrvbr.f90 processed by SPAG 8.04RA 15:47 21 Feb 2024
!!SPAG Open source Personal, Educational or Academic User Clemson University  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
!-----------------------------------------------------------------------
SUBROUTINE csrvbr(N,Ia,Ja,A,Nr,Nc,Kvstr,Kvstc,Ib,Jb,Kb,B,Job,Iwk,Nkmax,Nzmax,Ierr)
   USE ISO_FORTRAN_ENV                 
   IMPLICIT NONE
!
! Dummy argument declarations rewritten by SPAG
!
   INTEGER :: N
   INTEGER , INTENT(IN) :: Nkmax
   INTEGER , INTENT(IN) :: Nzmax
   INTEGER , DIMENSION(N+1) :: Ia
   INTEGER , DIMENSION(*) :: Ja
   REAL(REAL64) , DIMENSION(*) :: A
   INTEGER , INTENT(INOUT) :: Nr
   INTEGER , INTENT(INOUT) :: Nc
   INTEGER , INTENT(INOUT) , DIMENSION(*) :: Kvstr
   INTEGER , INTENT(INOUT) , DIMENSION(*) :: Kvstc
   INTEGER , INTENT(INOUT) , DIMENSION(*) :: Ib
   INTEGER , INTENT(INOUT) , DIMENSION(Nkmax-1) :: Jb
   INTEGER , INTENT(INOUT) , DIMENSION(Nkmax) :: Kb
   REAL(REAL64) , INTENT(OUT) , DIMENSION(Nzmax) :: B
   INTEGER , INTENT(IN) :: Job
   INTEGER , INTENT(INOUT) , DIMENSION(*) :: Iwk
   INTEGER , INTENT(OUT) :: Ierr
!
! Local variable declarations rewritten by SPAG
!
   INTEGER :: a0 , b0 , b1 , i , ii , j , jj , jnew , k0 , nb , ncol , neqr , numc
   LOGICAL :: sorted
   EXTERNAL csort , csorted , csrkvstc , csrkvstr , kvstmerge
!
! End of declarations rewritten by SPAG
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!     Converts compressed sparse row to variable block row format.
!-----------------------------------------------------------------------
!     On entry:
!--------------
!     n       = number of matrix rows
!     ia,ja,a = input matrix in CSR format
!
!     job     = job indicator.
!               If job=0, kvstr and kvstc are used as supplied.
!               If job=1, kvstr and kvstc are determined by the code.
!               If job=2, a conformal row/col partitioning is found and
!               returned in both kvstr and kvstc.  In the latter two cases,
!               an optimized algorithm can be used to perform the
!               conversion because all blocks are full.
!
!     nkmax   = size of supplied jb and kb arrays
!     nzmax   = size of supplied b array
!
!     If job=0 then the following are input:
!     nr,nc   = matrix block row and block column dimension
!     kvstr   = first row number for each block row
!     kvstc   = first column number for each block column.
!               (kvstr and kvstc may be the same array)
!
!     On return:
!---------------
!
!     ib,jb,kb,b = output matrix in VBR format
!
!     ierr    = error message
!               ierr = 0 means normal return
!               ierr = 1 out of space in jb and/or kb arrays
!               ierr = 2 out of space in b array
!               ierr = 3 nonsquare matrix used with job=2
!
!     If job=1,2 then the following are output:
!     nr,nc   = matrix block row and block column dimension
!     kvstr   = first row number for each block row
!     kvstc   = first column number for each block column
!               If job=2, then kvstr and kvstc contain the same info.
!
!     Work space:
!----------------
!     iwk(1:ncol) = inverse kvstc array.  If job=1,2 then we also need:
!     iwk(ncol+1:ncol+nr) = used to help determine sparsity of each block row.
!     The workspace is not assumed to be initialized to zero, nor is it
!     left that way.
!
!     Algorithms:
!----------------
!     There are two conversion codes in this routine.  The first assumes
!     that all blocks are full (there is a nonzero in the CSR data
!     structure for each entry in the block), and is used if the routine
!     determines the block partitioning itself.  The second code makes
!     no assumptions about the block partitioning, and is used if the
!     caller provides the partitioning.  The second code is much less
!     efficient than the first code.
!
!     In the first code, the CSR data structure is traversed sequentially
!     and entries are placed into the VBR data structure with stride
!     equal to the row dimension of the block row.  The columns of the
!     CSR data structure are sorted first if necessary.
!
!     In the second code, the block sparsity pattern is first determined.
!     This is done by traversing the CSR data structure and using an
!     implied linked list to determine which blocks are nonzero.  Then
!     the VBR data structure is filled by mapping each individual entry
!     in the CSR data structure into the VBR data structure.  The columns
!     of the CSR data structure are sorted first if necessary.
!
!-----------------------------------------------------------------------
!     Local variables:
!---------------------
!
!     ncol = number of scalar columns in matrix
!     nb = number of blocks in conformal row/col partitioning
!     neqr = number of rows in block row
!     numc = number of nonzero columns in row
!     a0 = index for entries in CSR a array
!     b0 = index for entries in VBR b array
!     b1 = temp
!     k0 = index for entries in VBR kb array
!     i  = loop index for block rows
!     ii = loop index for scalar rows in block row
!     j  = loop index for block columns
!     jj = loop index for scalar columns in block column
!     jnew = block column number
!     sorted = used to indicate if matrix already sorted by columns
!
!-----------------------------------------------------------------------
   Ierr = 0
!-----sort matrix by column indices
   CALL csorted(N,Ia,Ja,sorted)
   IF ( .NOT.sorted ) CALL csort(N,A,Ja,Ia,.TRUE.)
   IF ( Job==1 .OR. Job==2 ) THEN
!--------need to zero workspace; first find ncol
      ncol = 0
      DO i = 2 , N
         ncol = max0(ncol,Ja(Ia(i)-1))
      ENDDO
      DO i = 1 , ncol
         Iwk(i) = 0
      ENDDO
      CALL csrkvstr(N,Ia,Ja,Nr,Kvstr)
      CALL csrkvstc(N,Ia,Ja,Nc,Kvstc,Iwk)
   ENDIF
!-----check if want conformal partitioning
   IF ( Job==2 ) THEN
      IF ( Kvstr(Nr+1)/=Kvstc(Nc+1) ) THEN
         Ierr = 3
         RETURN
      ENDIF
!        use iwk temporarily
      CALL kvstmerge(Nr,Kvstr,Nc,Kvstc,nb,Iwk)
      Nr = nb
      Nc = nb
      DO i = 1 , nb + 1
         Kvstr(i) = Iwk(i)
         Kvstc(i) = Iwk(i)
      ENDDO
   ENDIF
!-----------------------------------------------------------------------
!     inverse kvst (scalar col number) = block col number
!     stored in iwk(1:n)
!-----------------------------------------------------------------------
   DO i = 1 , Nc
      DO j = Kvstc(i) , Kvstc(i+1) - 1
         Iwk(j) = i
      ENDDO
   ENDDO
   ncol = Kvstc(Nc+1) - 1
!-----jump to conversion routine
   IF ( Job==0 ) THEN
!-----------------------------------------------------------------------
!     Conversion for user supplied block partitioning
!-----------------------------------------------------------------------
!-----initialize workspace for sparsity indicator
      DO i = ncol + 1 , ncol + Nc
         Iwk(i) = 0
      ENDDO
      k0 = 1
      Kb(1) = 1
!-----find sparsity of block rows
      DO i = 1 , Nr
         neqr = Kvstr(i+1) - Kvstr(i)
         numc = Ia(Kvstr(i)+1) - Ia(Kvstr(i))
         Ib(i) = k0
!--------loop on all the elements in the block row to determine block sparsity
         DO jj = Ia(Kvstr(i)) , Ia(Kvstr(i+1)) - 1
            Iwk(Iwk(Ja(jj))+ncol) = 1
         ENDDO
!--------use sparsity to set jb and kb arrays
         DO j = 1 , Nc
            IF ( Iwk(j+ncol)/=0 ) THEN
!--------------check there is enough space in kb and jb arrays
               IF ( k0+1>Nkmax ) THEN
                  Ierr = 1
                  WRITE (*,*) 'csrvbr: no space in kb for block row ' , i
                  RETURN
               ENDIF
               Kb(k0+1) = Kb(k0) + neqr*(Kvstc(j+1)-Kvstc(j))
               Jb(k0) = j
               k0 = k0 + 1
               Iwk(j+ncol) = 0
            ENDIF
         ENDDO
      ENDDO
      Ib(Nr+1) = k0
!-----Fill b with entries from a by traversing VBR data structure.
      a0 = 1
!-----loop on block rows
      DO i = 1 , Nr
         neqr = Kvstr(i+1) - Kvstr(i)
!--------loop on scalar rows in block row
         DO ii = 0 , neqr - 1
            b0 = Kb(Ib(i)) + ii
!-----------loop on block columns
            DO j = Ib(i) , Ib(i+1) - 1
!--------------loop on scalar columns within block column
               DO jj = Kvstc(Jb(j)) , Kvstc(Jb(j)+1) - 1
!-----------------check there is enough space in b array
                  IF ( b0>Nzmax ) THEN
                     Ierr = 2
                     WRITE (*,*) 'csrvbr: no space in b for blk row' , i
                     RETURN
                  ENDIF
                  IF ( a0>=Ia(Kvstr(i)+ii+1) ) THEN
                     B(b0) = 0.D0
                  ELSEIF ( jj==Ja(a0) ) THEN
                     B(b0) = A(a0)
                     a0 = a0 + 1
                  ELSE
                     B(b0) = 0.D0
                  ENDIF
                  b0 = b0 + neqr
!--------------endloop on scalar columns
               ENDDO
!-----------endloop on block columns
            ENDDO
         ENDDO
      ENDDO
      RETURN
   ENDIF
!-----------------------------------------------------------------------
!     Fast conversion for computed block partitioning
!-----------------------------------------------------------------------
   a0 = 1
   b0 = 1
   k0 = 1
   Kb(1) = 1
!-----loop on block rows
   DO i = 1 , Nr
      neqr = Kvstr(i+1) - Kvstr(i)
      numc = Ia(Kvstr(i)+1) - Ia(Kvstr(i))
      Ib(i) = k0
!--------loop on first row in block row to determine block sparsity
      j = 0
      DO jj = Ia(Kvstr(i)) , Ia(Kvstr(i)+1) - 1
         jnew = Iwk(Ja(jj))
         IF ( jnew/=j ) THEN
!--------------check there is enough space in kb and jb arrays
            IF ( k0+1>Nkmax ) THEN
               Ierr = 1
               WRITE (*,*) 'csrvbr: no space in kb for block row ' , i
               RETURN
            ENDIF
!--------------set entries for this block
            j = jnew
            b0 = b0 + neqr*(Kvstc(j+1)-Kvstc(j))
            Kb(k0+1) = b0
            Jb(k0) = j
            k0 = k0 + 1
         ENDIF
      ENDDO
!--------loop on scalar rows in block row
      DO ii = 0 , neqr - 1
         b1 = Kb(Ib(i)) + ii
!-----------loop on elements in a scalar row
         DO jj = 1 , numc
!--------------check there is enough space in b array
            IF ( b1>Nzmax ) THEN
               Ierr = 2
               WRITE (*,*) 'csrvbr: no space in b for block row ' , i
               RETURN
            ENDIF
            B(b1) = A(a0)
            b1 = b1 + neqr
            a0 = a0 + 1
         ENDDO
      ENDDO
   ENDDO
   Ib(Nr+1) = k0
   RETURN
END SUBROUTINE csrvbr
!*==vbrcsr.f90 processed by SPAG 8.04RA 15:47 21 Feb 2024
!!SPAG Open source Personal, Educational or Academic User Clemson University  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
!-----------------------------------------------------------------------
!----------------------------end-of-csrvbr------------------------------
!----------------------------------------------------------------------c
SUBROUTINE vbrcsr(Ia,Ja,A,Nr,Kvstr,Kvstc,Ib,Jb,Kb,B,Nzmax,Ierr)
   USE ISO_FORTRAN_ENV                 
   IMPLICIT NONE
!
! Dummy argument declarations rewritten by SPAG
!
   INTEGER , INTENT(IN) :: Nr
   INTEGER , INTENT(IN) :: Nzmax
   INTEGER , INTENT(OUT) , DIMENSION(*) :: Ia
   INTEGER , INTENT(INOUT) , DIMENSION(*) :: Ja
   REAL(REAL64) , INTENT(OUT) , DIMENSION(*) :: A
   INTEGER , INTENT(IN) , DIMENSION(Nr+1) :: Kvstr
   INTEGER , INTENT(IN) , DIMENSION(*) :: Kvstc
   INTEGER , INTENT(IN) , DIMENSION(Nr+1) :: Ib
   INTEGER , INTENT(IN) , DIMENSION(*) :: Jb
   INTEGER , INTENT(IN) , DIMENSION(*) :: Kb
   REAL(REAL64) , INTENT(IN) , DIMENSION(Nzmax) :: B
   INTEGER , INTENT(INOUT) :: Ierr
!
! Local variable declarations rewritten by SPAG
!
   INTEGER :: a0 , b0 , i , ii , j , jj , neqr , numc
!
! End of declarations rewritten by SPAG
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!     Converts variable block row to compressed sparse row format.
!-----------------------------------------------------------------------
!     On entry:
!--------------
!     nr      = number of block rows
!     kvstr   = first row number for each block row
!     kvstc   = first column number for each block column
!     ib,jb,kb,b = input matrix in VBR format
!     nzmax   = size of supplied ja and a arrays
!
!     On return:
!---------------
!     ia,ja,a = output matrix in CSR format
!
!     ierr    = error message
!               ierr = 0 means normal return
!               ierr = negative row number when out of space in
!                      ja and a arrays
!
!     Work space:
!----------------
!     None
!
!     Algorithm:
!---------------
!     The VBR data structure is traversed in the order that is required
!     to fill the CSR data structure.  In a given block row, consecutive
!     entries in the CSR data structure are entries in the VBR data
!     structure with stride equal to the row dimension of the block.
!     The VBR data structure is assumed to be sorted by block columns.
!
!-----------------------------------------------------------------------
!     Local variables:
!---------------------
!
!     neqr = number of rows in block row
!     numc = number of nonzero columns in row
!     a0 = index for entries in CSR a array
!     b0 = index for entries in VBR b array
!     i  = loop index for block rows
!     ii = loop index for scalar rows in block row
!     j  = loop index for block columns
!     jj = loop index for scalar columns in block column
!
!-----------------------------------------------------------------------
   Ierr = 0
   a0 = 1
   b0 = 1
!-----loop on block rows
   DO i = 1 , Nr
!--------set num of rows in block row, and num of nonzero cols in row
      neqr = Kvstr(i+1) - Kvstr(i)
      numc = (Kb(Ib(i+1))-Kb(Ib(i)))/neqr
!--------construct ja for a scalar row
      DO j = Ib(i) , Ib(i+1) - 1
         DO jj = Kvstc(Jb(j)) , Kvstc(Jb(j)+1) - 1
            Ja(a0) = jj
            a0 = a0 + 1
         ENDDO
      ENDDO
!--------construct neqr-1 additional copies of ja for the block row
      DO ii = 1 , neqr - 1
         DO j = 1 , numc
            Ja(a0) = Ja(a0-numc)
            a0 = a0 + 1
         ENDDO
      ENDDO
!--------reset a0 back to beginning of block row
      a0 = Kb(Ib(i))
!--------loop on scalar rows in block row
      DO ii = 0 , neqr - 1
         Ia(Kvstr(i)+ii) = a0
         b0 = Kb(Ib(i)) + ii
!-----------loop on elements in a scalar row
         DO jj = 1 , numc
!--------------check there is enough space in a array
            IF ( a0>Nzmax ) THEN
               Ierr = -(Kvstr(i)+ii)
               WRITE (*,*) 'vbrcsr: no space for row ' , -Ierr
               RETURN
            ENDIF
            A(a0) = B(b0)
            a0 = a0 + 1
            b0 = b0 + neqr
         ENDDO
      ENDDO
!-----endloop on block rows
   ENDDO
   Ia(Kvstr(Nr+1)) = a0
END SUBROUTINE vbrcsr
!*==csorted.f90 processed by SPAG 8.04RA 15:47 21 Feb 2024
!!SPAG Open source Personal, Educational or Academic User Clemson University  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
!-----------------------------------------------------------------------
!---------------------------end-of-vbrcsr-------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
SUBROUTINE csorted(N,Ia,Ja,Sorted)
   IMPLICIT NONE
!
! Dummy argument declarations rewritten by SPAG
!
   INTEGER , INTENT(IN) :: N
   INTEGER , INTENT(IN) , DIMENSION(N+1) :: Ia
   INTEGER , INTENT(IN) , DIMENSION(*) :: Ja
   LOGICAL , INTENT(OUT) :: Sorted
!
! Local variable declarations rewritten by SPAG
!
   INTEGER :: i , j
!
! End of declarations rewritten by SPAG
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!     Checks if matrix in CSR format is sorted by columns.
!-----------------------------------------------------------------------
!     On entry:
!--------------
!     n       = number of rows in matrix
!     ia, ja  = sparsity structure of matrix in CSR format
!
!     On return:
!---------------
!     sorted  = indicates if matrix is sorted by columns
!
!-----------------------------------------------------------------------
!-----local variables
!---------------------------------
   DO i = 1 , N
      DO j = Ia(i) + 1 , Ia(i+1) - 1
         IF ( Ja(j-1)>=Ja(j) ) THEN
            Sorted = .FALSE.
            RETURN
         ENDIF
      ENDDO
   ENDDO
   Sorted = .TRUE.
END SUBROUTINE csorted
