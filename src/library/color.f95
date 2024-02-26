!*==multic.f90 processed by SPAG 8.04RA 15:47 21 Feb 2024
!!SPAG Open source Personal, Educational or Academic User Clemson University  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
!----------------------------------------------------------------------c
!                          S P A R S K I T                             c
!----------------------------------------------------------------------c
!          REORDERING ROUTINES -- COLORING BASED ROUTINES              c
!----------------------------------------------------------------------c
! contents:                                                            c
!----------                                                            c
! multic  : greedy algorithm for multicoloring                         c
! indset0 : greedy algorithm for independent set ordering              c
! indset1 : independent set ordering using minimal degree traversal    c
! indset2 : independent set ordering with local minimization           c
! indset3 : independent set ordering by vertex cover algorithm         c
! HeapSort, FixHeap, HeapInsert, interchange, MoveBack, FiHeapM,       c
!           FixHeapM, HeapInsertM,indsetr,rndperm, are utility         c
!           routines for sorting, generating random permutations, etc. c
!----------------------------------------------------------------------c
SUBROUTINE multic(N,Ja,Ia,Ncol,Kolrs,Il,Iord,Maxcol,Ierr)
   IMPLICIT NONE
!
! Dummy argument declarations rewritten by SPAG
!
   INTEGER , INTENT(IN) :: N
   INTEGER , INTENT(IN) :: Maxcol
   INTEGER , INTENT(IN) , DIMENSION(*) :: Ja
   INTEGER , INTENT(IN) , DIMENSION(N+1) :: Ia
   INTEGER , INTENT(INOUT) :: Ncol
   INTEGER , INTENT(INOUT) , DIMENSION(N) :: Kolrs
   INTEGER , INTENT(INOUT) , DIMENSION(Maxcol+1) :: Il
   INTEGER , INTENT(INOUT) , DIMENSION(N) :: Iord
   INTEGER , INTENT(OUT) :: Ierr
!
! Local variable declarations rewritten by SPAG
!
   INTEGER :: i , icol , ii , j , k , kol , mcol , mycol
!
! End of declarations rewritten by SPAG
!
!-----------------------------------------------------------------------
!     multicoloring ordering -- greedy algorithm --
!     determines the coloring permutation and sets up
!     corresponding data structures for it.
!-----------------------------------------------------------------------
! on entry
! --------
! n     = row and column dimention of matrix
! ja    = column indices of nonzero elements of matrix, stored rowwise.
! ia    = pointer to beginning of each row in ja.
! maxcol= maximum number of colors allowed -- the size of il is
!         maxcol+1 at least. Note: the number of colors does not
!         exceed the maximum degree of each node +1.
! iord  = en entry iord gives the order of traversal of the nodes
!         in the multicoloring algorithm. If there is no preference
!         then set iord(j)=j for j=1,...,n
!
! on return
! ---------
! ncol  = number of colours found
! kolrs = integer array containing the color number assigned to each node
! il    = integer array containing the pointers to the
!         beginning of each color set. In the permuted matrix
!         the rows /columns il(kol) to il(kol+1)-1 have the same color.
! iord  = permutation array corresponding to the multicolor ordering.
!         row number i will become row nbumber iord(i) in permuted
!         matrix. (iord = destination permutation array).
! ierr  = integer. Error message. normal return ierr = 0. If ierr .eq.1
!         then the array il was overfilled.
!
!-----------------------------------------------------------------------
!
   INTEGER :: spag_nextblock_1
   INTEGER :: spag_nextblock_2
   spag_nextblock_1 = 1
   SPAG_DispatchLoop_1: DO
      SELECT CASE (spag_nextblock_1)
      CASE (1)
!
         Ierr = 0
         DO j = 1 , N
            Kolrs(j) = 0
         ENDDO
         DO j = 1 , Maxcol
            Il(j) = 0
         ENDDO
!
         Ncol = 0
!
!     scan all nodes
!
         DO ii = 1 , N
            spag_nextblock_2 = 1
            SPAG_DispatchLoop_2: DO
               SELECT CASE (spag_nextblock_2)
               CASE (1)
                  i = Iord(ii)
!
!     look at adjacent nodes to determine colors already assigned
!
                  mcol = 0
                  DO k = Ia(i) , Ia(i+1) - 1
                     j = Ja(k)
                     icol = Kolrs(j)
                     IF ( icol/=0 ) THEN
                        mcol = max(mcol,icol)
!
!     il used as temporary to record already assigned colors.
!
                        Il(icol) = 1
                     ENDIF
                  ENDDO
!
!     taken colors determined. scan il until a slot opens up.
!
                  mycol = 1
                  DO WHILE ( Il(mycol)==1 )
                     mycol = mycol + 1
                     IF ( mycol>Maxcol ) THEN
                        spag_nextblock_1 = 2
                        CYCLE SPAG_DispatchLoop_1
                     ENDIF
                     IF ( mycol>mcol ) THEN
                        spag_nextblock_2 = 2
                        CYCLE SPAG_DispatchLoop_2
                     ENDIF
                  ENDDO
                  spag_nextblock_2 = 2
               CASE (2)
!
!     reset il to zero for next nodes
!
                  DO j = 1 , mcol
                     Il(j) = 0
                  ENDDO
!
!     assign color and update number of colors so far
!
                  Kolrs(i) = mycol
                  Ncol = max(Ncol,mycol)
                  EXIT SPAG_DispatchLoop_2
               END SELECT
            ENDDO SPAG_DispatchLoop_2
         ENDDO
!
!     every node has now been colored. Count nodes of each color
!
         DO j = 1 , N
            kol = Kolrs(j) + 1
            Il(kol) = Il(kol) + 1
         ENDDO
!
!     set pointers il
!
         Il(1) = 1
         DO j = 1 , Ncol
            Il(j+1) = Il(j) + Il(j+1)
         ENDDO
!
!     set iord
!
         DO j = 1 , N
            kol = Kolrs(j)
            Iord(j) = Il(kol)
            Il(kol) = Il(kol) + 1
         ENDDO
!
!     shift il back
!
         DO j = Ncol , 1 , -1
            Il(j+1) = Il(j)
         ENDDO
         Il(1) = 1
!
         RETURN
      CASE (2)
         Ierr = 1
         EXIT SPAG_DispatchLoop_1
      END SELECT
   ENDDO SPAG_DispatchLoop_1
!----end-of-multic------------------------------------------------------
!-----------------------------------------------------------------------
END SUBROUTINE multic
!*==indset0.f90 processed by SPAG 8.04RA 15:47 21 Feb 2024
!!SPAG Open source Personal, Educational or Academic User Clemson University  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
!-----------------------------------------------------------------------
 
SUBROUTINE indset0(N,Ja,Ia,Nset,Iord,Riord,Sym,Iptr)
   IMPLICIT NONE
!
! Dummy argument declarations rewritten by SPAG
!
   INTEGER , INTENT(IN) :: N
   INTEGER , INTENT(IN) , DIMENSION(*) :: Ja
   INTEGER , INTENT(IN) , DIMENSION(*) :: Ia
   INTEGER , INTENT(INOUT) :: Nset
   INTEGER , INTENT(INOUT) , DIMENSION(*) :: Iord
   INTEGER , INTENT(OUT) , DIMENSION(*) :: Riord
   LOGICAL , INTENT(IN) :: Sym
   INTEGER , INTENT(IN) :: Iptr
!
! Local variable declarations rewritten by SPAG
!
   INTEGER :: ipos , j , k , k1 , k2 , mat , nod , nummat
!
! End of declarations rewritten by SPAG
!
!----------------------------------------------------------------------
! greedy algorithm for independent set ordering
!----------------------------------------------------------------------
! parameters:
! ----------
! n      = row dimension of matrix
! ja, ia = matrix pattern in CRS format
! nset   = (output) number of elements in the independent set
! iord   = permutation array corresponding to the independent set
!          ordering. Row number i will become row number iord(i) in
!          permuted matrix.
! riord  = reverse permutation array. Row number i in the permutated
!          matrix is row number riord(i) in original matrix.
!----------------------------------------------------------------------
! notes: works for CSR, MSR, and CSC formats but assumes that the
! matrix has a symmetric structure.
!----------------------------------------------------------------------
! local variables
!
   DO j = 1 , N
      Iord(j) = 0
   ENDDO
   nummat = 1
   IF ( .NOT.Sym ) nummat = 2
!
!     iord used as a marker
!
   Nset = 0
   DO nod = 1 , N
      IF ( Iord(nod)==0 ) THEN
         Nset = Nset + 1
         Iord(nod) = 1
!
!     visit all neighbors of current nod
!
         ipos = 0
         DO mat = 1 , nummat
            DO k = Ia(ipos+nod) , Ia(ipos+nod+1) - 1
               j = Ja(k)
               IF ( j/=nod ) Iord(j) = 2
            ENDDO
            ipos = Iptr - 1
         ENDDO
      ENDIF
   ENDDO
!
!     get permutation
!
   k1 = 0
   k2 = Nset
   DO j = 1 , N
      IF ( Iord(j)==1 ) THEN
         k1 = k1 + 1
         k = k1
      ELSE
         k2 = k2 + 1
         k = k2
      ENDIF
      Riord(k) = j
      Iord(j) = k
   ENDDO
!----------------------------------------------------------------------
END SUBROUTINE indset0
!*==indset1.f90 processed by SPAG 8.04RA 15:47 21 Feb 2024
!!SPAG Open source Personal, Educational or Academic User Clemson University  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
!----------------------------------------------------------------------
SUBROUTINE indset1(N,Ja,Ia,Nset,Iord,Riord,Iw,Sym,Iptr)
   IMPLICIT NONE
!
! Dummy argument declarations rewritten by SPAG
!
   INTEGER :: N
   INTEGER , INTENT(IN) , DIMENSION(*) :: Ja
   INTEGER , INTENT(IN) , DIMENSION(*) :: Ia
   INTEGER , INTENT(INOUT) :: Nset
   INTEGER , INTENT(INOUT) , DIMENSION(*) :: Iord
   INTEGER , DIMENSION(*) :: Riord
   INTEGER , INTENT(INOUT) , DIMENSION(*) :: Iw
   LOGICAL , INTENT(IN) :: Sym
   INTEGER , INTENT(IN) :: Iptr
!
! Local variable declarations rewritten by SPAG
!
   INTEGER :: ii , imat , ipos , iptrm1 , j , k , k1 , k2 , mat , nod , nummat
   EXTERNAL heapsort
!
! End of declarations rewritten by SPAG
!
!----------------------------------------------------------------------
! greedy algorithm for independent set ordering -- with intial
! order of traversal given by that of min degree.
!----------------------------------------------------------------------
! parameters:
! ----------
! n      = row dimension of matrix
! ja, ia = matrix pattern in CRS format
! nset   = (output) number of elements in the independent set
! iord   = permutation array corresponding to the independent set
!          ordering. Row number i will become row number iord(i) in
!          permuted matrix.
! riord  = reverse permutation array. Row number i in the permutated
!          matrix is row number riord(i) in original matrix.
!----------------------------------------------------------------------
! notes: works for CSR, MSR, and CSC formats but assumes that the
! matrix has a symmetric structure.
!----------------------------------------------------------------------
! local variables
!
!     nummat is the number of matrices to loop through (A in symmetric
!     pattern case (nummat=1) or A,and transp(A) otherwise (mummat=2)
!
   IF ( Sym ) THEN
      nummat = 1
   ELSE
      nummat = 2
   ENDIF
   iptrm1 = Iptr - 1
!
!     initialize arrays
!
   DO j = 1 , N
      Iord(j) = j
      Riord(j) = j
      Iw(j) = 0
   ENDDO
!
!     initialize degrees of all nodes
!
   ipos = 0
   DO imat = 1 , nummat
      DO j = 1 , N
         Iw(j) = Iw(j) + Ia(ipos+j+1) - Ia(ipos+j)
      ENDDO
      ipos = iptrm1
   ENDDO
!
!     call heapsort -- sorts nodes in increasing degree.
!
   CALL heapsort(Iw,Iord,Riord,N,N)
!
!     weights no longer needed -- use iw to store order of traversal.
!
   DO j = 1 , N
      Iw(N-j+1) = Iord(j)
      Iord(j) = 0
   ENDDO
!
!     iord used as a marker
!
   Nset = 0
   DO ii = 1 , N
      nod = Iw(ii)
      IF ( Iord(nod)==0 ) THEN
         Nset = Nset + 1
         Iord(nod) = 1
!
!     visit all neighbors of current nod
!
         ipos = 0
         DO mat = 1 , nummat
            DO k = Ia(ipos+nod) , Ia(ipos+nod+1) - 1
               j = Ja(k)
               IF ( j/=nod ) Iord(j) = 2
            ENDDO
            ipos = iptrm1
         ENDDO
      ENDIF
   ENDDO
!
!     get permutation
!
   k1 = 0
   k2 = Nset
   DO j = 1 , N
      IF ( Iord(j)==1 ) THEN
         k1 = k1 + 1
         k = k1
      ELSE
         k2 = k2 + 1
         k = k2
      ENDIF
      Riord(k) = j
      Iord(j) = k
   ENDDO
!----------------------------------------------------------------------
END SUBROUTINE indset1
!*==indset2.f90 processed by SPAG 8.04RA 15:47 21 Feb 2024
!!SPAG Open source Personal, Educational or Academic User Clemson University  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
!----------------------------------------------------------------------
SUBROUTINE indset2(N,Ja,Ia,Nset,Iord,Riord,Iw,Sym,Iptr)
   IMPLICIT NONE
!
! Dummy argument declarations rewritten by SPAG
!
   INTEGER :: N
   INTEGER , INTENT(IN) , DIMENSION(*) :: Ja
   INTEGER , INTENT(IN) , DIMENSION(*) :: Ia
   INTEGER , INTENT(INOUT) :: Nset
   INTEGER , INTENT(INOUT) , DIMENSION(N) :: Iord
   INTEGER , INTENT(INOUT) , DIMENSION(N) :: Riord
   INTEGER , INTENT(INOUT) , DIMENSION(N) :: Iw
   LOGICAL , INTENT(IN) :: Sym
   INTEGER , INTENT(IN) :: Iptr
!
! Local variable declarations rewritten by SPAG
!
   INTEGER :: i , ichild , imat , ipos , iptrm1 , j , jn , jnew , jo , jold , k , k1 , k2 , last , lastlast , nod , nummat
   EXTERNAL fixheap , heapinsert , moveback
!
! End of declarations rewritten by SPAG
!
!----------------------------------------------------------------------
! greedy algorithm for independent set ordering -- local minimization
! using heap strategy --
!----------------------------------------------------------------------
! This version for BOTH unsymmetric and symmetric patterns
!----------------------------------------------------------------------
! on entry
! --------
! n     = row and column dimension of matrix
! ja    = column indices of nonzero elements of matrix,stored rowwise.
! ia    = pointer to beginning of each row in ja.
! sym   = logical indicating whether the matrix has a symmetric pattern.
!         If not the transpose must also be provided -- appended to the
!         ja, ia structure -- see description of iptr next.
! iptr  = in case the matrix has an unsymmetric pattern,the transpose
!         is assumed to be stored in the same arrays ia,ja. iptr is the
!         location in ia of the pointer to the first row of transp(A).
!         more generally, ia(iptr),...,ia(iptr+n) are the pointers to
!         the beginnings of rows 1, 2, ...., n+1 (row n+1 is fictitious)
!         of the transpose of A in the array ja. For example,when using
!         the msr format,one can write:
!          iptr = ja(n+1)
!          ipos = iptr+n+2                ! get the transpose of A:
!          call csrcsc (n,0,ipos,a,ja,ja,a,ja,ja(iptr))    ! and then:
!          call indset(n,ja,ja,nset,iord,riord,iwk,.false.,iptr)
!
! iw    = work space of length n.
!
! on return:
!----------
! nset  = integer. The number of unknowns in the independent set.
! iord  = permutation array corresponding to the new ordering. The
!         first nset unknowns correspond to the independent set.
! riord = reverse permutation array.
!----------------------------------------------------------------------
! local variables --
!
!
!     nummat is the number of matrices to loop through (A in symmetric
!     pattern case (nummat=1) or A,and transp(A) otherwise (mummat=2)
!
   IF ( Sym ) THEN
      nummat = 1
   ELSE
      nummat = 2
   ENDIF
   iptrm1 = Iptr - 1
!
!     initialize arrays
!
   DO j = 1 , N
      Iord(j) = j
      Riord(j) = j
      Iw(j) = 0
   ENDDO
!
!     initialize degrees of all nodes
!
   ipos = 0
   DO imat = 1 , nummat
      DO j = 1 , N
         Iw(j) = Iw(j) + Ia(ipos+j+1) - Ia(ipos+j)
      ENDDO
      ipos = iptrm1
   ENDDO
!
! start by constructing a heap
!
   DO i = N/2 , 1 , -1
      j = i
      CALL fixheap(Iw,Iord,Riord,j,j,N)
   ENDDO
!
! main loop -- remove nodes one by one.
!
   last = N
   Nset = 0
   SPAG_Loop_1_1: DO
      lastlast = last
      nod = Iord(1)
!
!     move first element to end
!
      CALL moveback(Iw,Iord,Riord,last)
      last = last - 1
      Nset = Nset + 1
!
!     scan all neighbors of accepted node -- move them to back --
!
      ipos = 0
      DO imat = 1 , nummat
         DO k = Ia(ipos+nod) , Ia(ipos+nod+1) - 1
            jold = Ja(k)
            jnew = Riord(jold)
            IF ( jold/=nod .AND. jnew<=last ) THEN
               Iw(jnew) = -1
               CALL heapinsert(Iw,Iord,Riord,jnew,ichild,jnew)
               CALL moveback(Iw,Iord,Riord,last)
               last = last - 1
            ENDIF
         ENDDO
         ipos = iptrm1
      ENDDO
!
! update the degree of each edge
!
      DO k = last + 1 , lastlast - 1
         jold = Iord(k)
!
!     scan the neighbors of current node
!
         ipos = 0
         DO imat = 1 , nummat
            DO i = Ia(ipos+jold) , Ia(ipos+jold+1) - 1
               jo = Ja(i)
               jn = Riord(jo)
!
!     consider this node only if it has not been moved
!
               IF ( jn<=last ) THEN
!     update degree of this neighbor
                  Iw(jn) = Iw(jn) - 1
!     and fix the heap accordingly
                  CALL heapinsert(Iw,Iord,Riord,jn,ichild,jn)
               ENDIF
            ENDDO
            ipos = iptrm1
         ENDDO
      ENDDO
!
!     stopping test -- end main "while"loop
!
      IF ( last<=1 ) THEN
         Nset = Nset + last
!
!     rescan all nodes one more time to determine the permutations
!
         k1 = 0
         k2 = Nset
         DO j = N , 1 , -1
            IF ( Iw(j)>=0 ) THEN
               k1 = k1 + 1
               k = k1
            ELSE
               k2 = k2 + 1
               k = k2
            ENDIF
            Riord(k) = Iord(j)
         ENDDO
         DO j = 1 , N
            Iord(Riord(j)) = j
         ENDDO
         EXIT SPAG_Loop_1_1
      ENDIF
   ENDDO SPAG_Loop_1_1
!----------------------------------------------------------------------
END SUBROUTINE indset2
!*==indset3.f90 processed by SPAG 8.04RA 15:47 21 Feb 2024
!!SPAG Open source Personal, Educational or Academic User Clemson University  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
!----------------------------------------------------------------------
SUBROUTINE indset3(N,Ja,Ia,Nset,Iord,Riord,Iw,Sym,Iptr)
   IMPLICIT NONE
!
! Dummy argument declarations rewritten by SPAG
!
   INTEGER :: N
   INTEGER , INTENT(IN) , DIMENSION(*) :: Ja
   INTEGER , INTENT(IN) , DIMENSION(*) :: Ia
   INTEGER , INTENT(INOUT) :: Nset
   INTEGER , INTENT(INOUT) , DIMENSION(N) :: Iord
   INTEGER , INTENT(INOUT) , DIMENSION(N) :: Riord
   INTEGER , INTENT(INOUT) , DIMENSION(N) :: Iw
   LOGICAL , INTENT(IN) :: Sym
   INTEGER , INTENT(IN) :: Iptr
!
! Local variable declarations rewritten by SPAG
!
   INTEGER :: i , ideg , imat , ipos , iptrm1 , j , jnew , jold , k , nnz , nod , nummat
   EXTERNAL fixheapm , movebackm
!
! End of declarations rewritten by SPAG
!
!----------------------------------------------------------------------
! greedy algorithm for independent set ordering -- local minimization
! using heap strategy -- VERTEX COVER ALGORITHM --
! ASSUMES MSR FORMAT (no diagonal element) -- ADD A SWITCH FOR CSR --
!----------------------------------------------------------------------
! This version for BOTH unsymmetric and symmetric patterns
!----------------------------------------------------------------------
! on entry
! --------
! n     = row and column dimension of matrix
! ja    = column indices of nonzero elements of matrix,stored rowwise.
! ia    = pointer to beginning of each row in ja.
! sym   = logical indicating whether the matrix has a symmetric pattern.
!         If not the transpose must also be provided -- appended to the
!         ja, ia structure -- see description of iptr next.
! iptr  = in case the matrix has an unsymmetric pattern,the transpose
!         is assumed to be stored in the same arrays ia,ja. iptr is the
!         location in ia of the pointer to the first row of transp(A).
!         more generally, ia(iptr),...,ia(iptr+n) are the pointers to
!         the beginnings of rows 1, 2, ...., n+1 (row n+1 is fictitious)
!         of the transpose of A in the array ja. For example,when using
!         the msr format,one can write:
!          iptr = ja(n+1)
!          ipos = iptr+n+2                ! get the transpose of A:
!          call csrcsc (n,0,ipos,a,ja,ja,a,ja,ja(iptr))    ! and then:
!          call indset(n,ja,ja,nset,iord,riord,iwk,.false.,iptr)
!
! iw    = work space of length n.
!
! on return:
!----------
! nset  = integer. The number of unknowns in the independent set.
! iord  = permutation array corresponding to the new ordering. The
!         first nset unknowns correspond to the independent set.
! riord = reverse permutation array.
!----------------------------------------------------------------------
! local variables --
!
!
!     nummat is the number of matrices to loop through (A in symmetric
!     pattern case (nummat=1) or A,and transp(A) otherwise (mummat=2)
!
   IF ( Sym ) THEN
      nummat = 1
   ELSE
      nummat = 2
   ENDIF
   iptrm1 = Iptr - 1
!
!     initialize arrays
!
   DO j = 1 , N
      Riord(j) = j
      Iord(j) = j
      Iw(j) = 0
   ENDDO
!
!     initialize degrees of all nodes
!
   nnz = 0
   ipos = 0
   DO imat = 1 , nummat
      DO j = 1 , N
         ideg = Ia(ipos+j+1) - Ia(ipos+j)
         Iw(j) = Iw(j) + ideg
         nnz = nnz + ideg
      ENDDO
      ipos = iptrm1
   ENDDO
!
!     number of edges
!
   IF ( Sym ) nnz = 2*nnz
!
! start by constructing a Max heap
!
   DO i = N/2 , 1 , -1
      j = i
      CALL fixheapm(Iw,Riord,Iord,j,j,N)
   ENDDO
   Nset = N
   SPAG_Loop_1_1: DO
!----------------------------------------------------------------------
! main loop -- remove nodes one by one.
!----------------------------------------------------------------------
!      lastnset = nset
      nod = Riord(1)
!
!     move first element to end
!
      CALL movebackm(Iw,Riord,Iord,Nset)
      nnz = nnz - Iw(Nset)
      Nset = Nset - 1
!
!     scan all neighbors of accepted node --
!
      ipos = 0
      DO imat = 1 , nummat
         DO k = Ia(ipos+nod) , Ia(ipos+nod+1) - 1
            jold = Ja(k)
            jnew = Iord(jold)
            IF ( jold/=nod .AND. jnew<=Nset ) THEN
               Iw(jnew) = Iw(jnew) - 1
               nnz = nnz - 1
               CALL fixheapm(Iw,Riord,Iord,jnew,jnew,Nset)
            ENDIF
         ENDDO
         ipos = iptrm1
      ENDDO
!
      IF ( nnz<=0 ) EXIT SPAG_Loop_1_1
   ENDDO SPAG_Loop_1_1
!-----------------------------------------------------------------------
END SUBROUTINE indset3
!*==heapsort.f90 processed by SPAG 8.04RA 15:47 21 Feb 2024
!!SPAG Open source Personal, Educational or Academic User Clemson University  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
!-----------------------------------------------------------------------
SUBROUTINE heapsort(A,Ind,Rind,N,Ncut)
   IMPLICIT NONE
!
! Dummy argument declarations rewritten by SPAG
!
   INTEGER :: N
   INTEGER , DIMENSION(*) :: A
   INTEGER , DIMENSION(N) :: Ind
   INTEGER , DIMENSION(N) :: Rind
   INTEGER , INTENT(IN) :: Ncut
!
! Local variable declarations rewritten by SPAG
!
   INTEGER :: i , j , jlast , last
   EXTERNAL fixheap , moveback
!
! End of declarations rewritten by SPAG
!
!----------------------------------------------------------------------
! integer version -- min heap sorts decreasinly.
!----------------------------------------------------------------------
! sorts inger keys in array a increasingly and permutes the companion
! array ind rind accrodingly.
! n    = size of array
! ncut = integer indicating when to cut the process.the process is
!        stopped after ncut outer steps of the heap-sort algorithm.
!        The first ncut values are sorted and they are the smallest
!        ncut values of the array.
!----------------------------------------------------------------------
! local variables
!
!
!    Heap sort algorithm ---
!
!    build heap
   DO i = N/2 , 1 , -1
      j = i
      CALL fixheap(A,Ind,Rind,j,j,N)
   ENDDO
!
!   done -- now remove keys one by one
!
   jlast = max(2,N-Ncut+1)
   DO last = N , jlast , -1
      CALL moveback(A,Ind,Rind,last)
   ENDDO
END SUBROUTINE heapsort
!*==fixheap.f90 processed by SPAG 8.04RA 15:47 21 Feb 2024
!!SPAG Open source Personal, Educational or Academic User Clemson University  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
!----------------------------------------------------------------------
SUBROUTINE fixheap(A,Ind,Rind,Jkey,Vacant,Last)
   IMPLICIT NONE
!
! Dummy argument declarations rewritten by SPAG
!
   INTEGER , INTENT(INOUT) , DIMENSION(*) :: A
   INTEGER , INTENT(INOUT) , DIMENSION(*) :: Ind
   INTEGER , INTENT(OUT) , DIMENSION(*) :: Rind
   INTEGER , INTENT(IN) :: Jkey
   INTEGER , INTENT(INOUT) :: Vacant
   INTEGER , INTENT(IN) :: Last
!
! Local variable declarations rewritten by SPAG
!
   INTEGER :: child , ikey , lchild , rchild , xkey
!
! End of declarations rewritten by SPAG
!
!----------------------------------------------------------------------
!     inserts a key (key and companion index) at the vacant position
!     in a (min) heap -
! arguments
!     a(1:last)    = real array
!     ind(1:last)  = integer array -- permutation of initial data
!     rind(1:last) = integer array -- reverse permutation
!     jkey         = position of key to be inserted. a(jkey)
!                    will be inserted into the heap
!     vacant       = vacant where a key is to be inserted
!     last         = number of elements in the heap.
!----------------------------------------------------------------------
! local variables
!
   xkey = A(Jkey)
   ikey = Ind(Jkey)
   lchild = 2*Vacant
   DO
      rchild = lchild + 1
      child = lchild
      IF ( rchild<=Last .AND. A(rchild)<A(child) ) child = rchild
      IF ( xkey<=A(child) .OR. child>Last ) THEN
         CALL spag_block_1
         RETURN
      ENDIF
      A(Vacant) = A(child)
      Ind(Vacant) = Ind(child)
      Rind(Ind(Vacant)) = Vacant
      Vacant = child
      lchild = 2*Vacant
      IF ( lchild>Last ) THEN
         CALL spag_block_1
         RETURN
      ENDIF
   ENDDO
   CALL spag_block_1
CONTAINS
   SUBROUTINE spag_block_1
      A(Vacant) = xkey
      Ind(Vacant) = ikey
      Rind(ikey) = Vacant
   END SUBROUTINE spag_block_1
!----------------------------------------------------------------------
END SUBROUTINE fixheap
!*==heapinsert.f90 processed by SPAG 8.04RA 15:47 21 Feb 2024
!!SPAG Open source Personal, Educational or Academic User Clemson University  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
!----------------------------------------------------------------------
SUBROUTINE heapinsert(A,Ind,Rind,Jkey,Child,Node)
   IMPLICIT NONE
!
! Dummy argument declarations rewritten by SPAG
!
   INTEGER , INTENT(INOUT) , DIMENSION(*) :: A
   INTEGER , INTENT(INOUT) , DIMENSION(*) :: Ind
   INTEGER , DIMENSION(*) :: Rind
   INTEGER , INTENT(IN) :: Jkey
   INTEGER , INTENT(INOUT) :: Child
   INTEGER , INTENT(IN) :: Node
!
! Local variable declarations rewritten by SPAG
!
   INTEGER :: ikey , parent , xkey
   EXTERNAL interchange
!
! End of declarations rewritten by SPAG
!
!----------------------------------------------------------------------
! inserts a key to a heap from `node'. Checks values up
! only -- i.e.,assumes that the subtree (if any) whose root
! is node is such that the keys are all inferior to those
! to ge inserted.
!
! child is where the key ended up.
!----------------------------------------------------------------------
!---- local variables
   xkey = A(Jkey)
   ikey = Ind(Jkey)
!      node = node + 1
   A(Node) = xkey
   Ind(Node) = ikey
   Rind(ikey) = Node
   IF ( Node<=1 ) RETURN
   Child = Node
   SPAG_Loop_1_1: DO
      parent = Child/2
      IF ( A(parent)<=A(Child) ) EXIT SPAG_Loop_1_1
      CALL interchange(A,Ind,Rind,Child,parent)
      Child = parent
      IF ( Child<=1 ) EXIT SPAG_Loop_1_1
   ENDDO SPAG_Loop_1_1
END SUBROUTINE heapinsert
!*==interchange.f90 processed by SPAG 8.04RA 15:47 21 Feb 2024
!!SPAG Open source Personal, Educational or Academic User Clemson University  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
!-----------------------------------------------------------------------
SUBROUTINE interchange(A,Ind,Rind,I,J)
   IMPLICIT NONE
!
! Dummy argument declarations rewritten by SPAG
!
   INTEGER , INTENT(INOUT) , DIMENSION(*) :: A
   INTEGER , INTENT(INOUT) , DIMENSION(*) :: Ind
   INTEGER , INTENT(OUT) , DIMENSION(*) :: Rind
   INTEGER , INTENT(IN) :: I
   INTEGER , INTENT(IN) :: J
!
! Local variable declarations rewritten by SPAG
!
   INTEGER :: itmp , tmp
!
! End of declarations rewritten by SPAG
!
   tmp = A(I)
   itmp = Ind(I)
!
   A(I) = A(J)
   Ind(I) = Ind(J)
!
   A(J) = tmp
   Ind(J) = itmp
   Rind(Ind(J)) = J
   Rind(Ind(I)) = I
!
END SUBROUTINE interchange
!*==moveback.f90 processed by SPAG 8.04RA 15:47 21 Feb 2024
!!SPAG Open source Personal, Educational or Academic User Clemson University  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
!----------------------------------------------------------------------
SUBROUTINE moveback(A,Ind,Rind,Last)
   IMPLICIT NONE
!
! Dummy argument declarations rewritten by SPAG
!
   INTEGER , INTENT(INOUT) , DIMENSION(*) :: A
   INTEGER , INTENT(INOUT) , DIMENSION(*) :: Ind
   INTEGER , DIMENSION(*) :: Rind
   INTEGER :: Last
!
! Local variable declarations rewritten by SPAG
!
   INTEGER :: imin , vacant , xmin
   EXTERNAL fixheap
!
! End of declarations rewritten by SPAG
!
! moves the front key to the back and inserts the last
! one back in from the top --
!
! local variables
!
!
   vacant = 1
   xmin = A(vacant)
   imin = Ind(vacant)
   CALL fixheap(A,Ind,Rind,Last,vacant,Last-1)
   A(Last) = xmin
   Ind(Last) = imin
   Rind(Ind(Last)) = Last
!
END SUBROUTINE moveback
!*==fixheapm.f90 processed by SPAG 8.04RA 15:47 21 Feb 2024
!!SPAG Open source Personal, Educational or Academic User Clemson University  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
!----------------------------------------------------------------------
SUBROUTINE fixheapm(A,Ind,Rind,Jkey,Vacant,Last)
   IMPLICIT NONE
!
! Dummy argument declarations rewritten by SPAG
!
   INTEGER , INTENT(INOUT) , DIMENSION(*) :: A
   INTEGER , INTENT(INOUT) , DIMENSION(*) :: Ind
   INTEGER , INTENT(OUT) , DIMENSION(*) :: Rind
   INTEGER , INTENT(IN) :: Jkey
   INTEGER , INTENT(INOUT) :: Vacant
   INTEGER , INTENT(IN) :: Last
!
! Local variable declarations rewritten by SPAG
!
   INTEGER :: child , ikey , lchild , rchild , xkey
!
! End of declarations rewritten by SPAG
!
!----
!     inserts a key (key and companion index) at the vacant position
!     in a heap -  THIS IS A MAX HEAP VERSION
! arguments
!     a(1:last)    = real array
!     ind(1:last)  = integer array -- permutation of initial data
!     rind(1:last) = integer array -- reverse permutation
!     jkey         = position of key to be inserted. a(jkey)
!                    will be inserted into the heap
!     vacant       = vacant where a key is to be inserted
!     last         = number of elements in the heap.
!----
! local variables
!
   xkey = A(Jkey)
   ikey = Ind(Jkey)
   lchild = 2*Vacant
   DO
      rchild = lchild + 1
      child = lchild
      IF ( rchild<=Last .AND. A(rchild)>A(child) ) child = rchild
      IF ( xkey>=A(child) .OR. child>Last ) THEN
         CALL spag_block_1
         RETURN
      ENDIF
      A(Vacant) = A(child)
      Ind(Vacant) = Ind(child)
      Rind(Ind(Vacant)) = Vacant
      Vacant = child
      lchild = 2*Vacant
      IF ( lchild>Last ) THEN
         CALL spag_block_1
         RETURN
      ENDIF
   ENDDO
   CALL spag_block_1
CONTAINS
   SUBROUTINE spag_block_1
      A(Vacant) = xkey
      Ind(Vacant) = ikey
      Rind(ikey) = Vacant
   END SUBROUTINE spag_block_1
END SUBROUTINE fixheapm
!*==heapinsertm.f90 processed by SPAG 8.04RA 15:47 21 Feb 2024
!!SPAG Open source Personal, Educational or Academic User Clemson University  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
!
SUBROUTINE heapinsertm(A,Ind,Rind,Jkey,Child,Node)
   IMPLICIT NONE
!
! Dummy argument declarations rewritten by SPAG
!
   INTEGER , INTENT(INOUT) , DIMENSION(*) :: A
   INTEGER , INTENT(INOUT) , DIMENSION(*) :: Ind
   INTEGER , DIMENSION(*) :: Rind
   INTEGER , INTENT(IN) :: Jkey
   INTEGER , INTENT(INOUT) :: Child
   INTEGER , INTENT(IN) :: Node
!
! Local variable declarations rewritten by SPAG
!
   INTEGER :: ikey , parent , xkey
   EXTERNAL interchange
!
! End of declarations rewritten by SPAG
!
!----------------------------------------------------------------------
! inserts a key to a heap from `node'. Checks values up
! only -- i.e.,assumes that the subtree (if any) whose root
! is node is such that the keys are all inferior to those
! to ge inserted.
!
! child is where the key ended up.
!----------------------------------------------------------------------
!---- local variables
   xkey = A(Jkey)
   ikey = Ind(Jkey)
!      node = node + 1
   A(Node) = xkey
   Ind(Node) = ikey
   Rind(ikey) = Node
   IF ( Node<=1 ) RETURN
   Child = Node
   SPAG_Loop_1_1: DO
      parent = Child/2
      IF ( A(parent)>=A(Child) ) EXIT SPAG_Loop_1_1
      CALL interchange(A,Ind,Rind,Child,parent)
      Child = parent
      IF ( Child<=1 ) EXIT SPAG_Loop_1_1
   ENDDO SPAG_Loop_1_1
END SUBROUTINE heapinsertm
!*==movebackm.f90 processed by SPAG 8.04RA 15:47 21 Feb 2024
!!SPAG Open source Personal, Educational or Academic User Clemson University  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
!----------------------------------------------------------------------
SUBROUTINE movebackm(A,Ind,Rind,Last)
   IMPLICIT NONE
!
! Dummy argument declarations rewritten by SPAG
!
   INTEGER , INTENT(INOUT) , DIMENSION(*) :: A
   INTEGER , INTENT(INOUT) , DIMENSION(*) :: Ind
   INTEGER , DIMENSION(*) :: Rind
   INTEGER :: Last
!
! Local variable declarations rewritten by SPAG
!
   INTEGER :: imin , vacant , xmin
   EXTERNAL fixheapm
!
! End of declarations rewritten by SPAG
!
!----------------------------------------------------------------------
! moves the front key to the back and inserts the last
! one back in from the top --  MAX HEAP VERSION
!----------------------------------------------------------------------
!
! local variables
!
!
   vacant = 1
   xmin = A(vacant)
   imin = Ind(vacant)
   CALL fixheapm(A,Ind,Rind,Last,vacant,Last-1)
   A(Last) = xmin
   Ind(Last) = imin
   Rind(Ind(Last)) = Last
!----------------------------------------------------------------------
END SUBROUTINE movebackm
!*==indsetr.f90 processed by SPAG 8.04RA 15:47 21 Feb 2024
!!SPAG Open source Personal, Educational or Academic User Clemson University  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
!----------------------------------------------------------------------
SUBROUTINE indsetr(N,Ja,Ia,Nset,Iord,Riord,Sym,Iptr)
   IMPLICIT NONE
!
! Dummy argument declarations rewritten by SPAG
!
   INTEGER :: N
   INTEGER , INTENT(IN) , DIMENSION(*) :: Ja
   INTEGER , INTENT(IN) , DIMENSION(*) :: Ia
   INTEGER , INTENT(INOUT) :: Nset
   INTEGER , INTENT(INOUT) , DIMENSION(*) :: Iord
   INTEGER , INTENT(INOUT) , DIMENSION(*) :: Riord
   LOGICAL , INTENT(IN) :: Sym
   INTEGER , INTENT(IN) :: Iptr
!
! Local variable declarations rewritten by SPAG
!
   INTEGER :: ii , ipos , iseed , j , k , k1 , k2 , mat , nod , nummat
   EXTERNAL rndperm
!
! End of declarations rewritten by SPAG
!
!----------------------------------------------------------------------
! greedy algorithm for independent set ordering -- RANDOM TRAVERSAL --
!----------------------------------------------------------------------
! parameters:
! ----------
! n      = row dimension of matrix
! ja, ia = matrix pattern in CRS format
! nset   = (output) number of elements in the independent set
! iord   = permutation array corresponding to the independent set
!          ordering. Row number i will become row number iord(i) in
!          permuted matrix.
! riord  = reverse permutation array. Row number i in the permutated
!          matrix is row number riord(i) in original matrix.
!----------------------------------------------------------------------
! notes: works for CSR, MSR, and CSC formats but assumes that the
! matrix has a symmetric structure.
!----------------------------------------------------------------------
! local variables
!
   DO j = 1 , N
      Iord(j) = 0
   ENDDO
!
! generate random permutation
!
   iseed = 0
   CALL rndperm(N,Riord,iseed)
   WRITE (8,'(10i6)') (Riord(j),j=1,N)
!
   nummat = 1
   IF ( .NOT.Sym ) nummat = 2
!
! iord used as a marker
!
   Nset = 0
   DO ii = 1 , N
      nod = Riord(ii)
      IF ( Iord(nod)==0 ) THEN
         Nset = Nset + 1
         Iord(nod) = 1
!
! visit all neighbors of current nod
!
         ipos = 0
         DO mat = 1 , nummat
            DO k = Ia(ipos+nod) , Ia(ipos+nod+1) - 1
               j = Ja(k)
               IF ( j/=nod ) Iord(j) = 2
            ENDDO
            ipos = Iptr - 1
         ENDDO
      ENDIF
   ENDDO
!
! get permutation
!
   k1 = 0
   k2 = Nset
   DO j = 1 , N
      IF ( Iord(j)==1 ) THEN
         k1 = k1 + 1
         k = k1
      ELSE
         k2 = k2 + 1
         k = k2
      ENDIF
      Riord(k) = j
      Iord(j) = k
   ENDDO
!----------------------------------------------------------------------
END SUBROUTINE indsetr
!*==rndperm.f90 processed by SPAG 8.04RA 15:47 21 Feb 2024
!!SPAG Open source Personal, Educational or Academic User Clemson University  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
!----------------------------------------------------------------------
SUBROUTINE rndperm(N,Iord,Iseed)
   IMPLICIT NONE
!
! Dummy argument declarations rewritten by SPAG
!
   INTEGER , INTENT(IN) :: N
   INTEGER , INTENT(INOUT) , DIMENSION(N) :: Iord
   INTEGER , INTENT(IN) :: Iseed
!
! Local variable declarations rewritten by SPAG
!
   INTEGER :: i , itmp , j
   INTEGER , EXTERNAL :: irand
!
! End of declarations rewritten by SPAG
!
!----------------------------------------------------------------------
! this subroutine will generate a pseudo random permutation of the
! n integers 1,2, ...,n.
! iseed is the initial seed. any integer.
!----------------------------------------------------------------------
! local
!
!----------------------------------------------------------------------
   DO j = 1 , N
      Iord(j) = j
   ENDDO
!
   DO i = 1 , N
      j = mod(Iseed+irand(0),N) + 1
      itmp = Iord(i)
      Iord(i) = Iord(j)
      Iord(j) = itmp
   ENDDO
!----------------------------------------------------------------------
!----------------------------------------------------------------------
END SUBROUTINE rndperm
