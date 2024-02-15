subroutine ilut(n,a,ja,ia,lfil,droptol,alu,jlu,ju,iwk,w,jw,ierr)
! -----------------------------------------------------------------------
     integer, intent(In) :: n, lfil, iwk
     real(kind=8), dimension(*), intent(In) :: a
     real(kind=8), dimension(*), intent(InOut) :: alu
     real(kind=8), dimension(n + 1), intent(Out) :: w
     real(kind=8), intent(In) :: droptol
     integer, dimension(*), intent(In) :: ja
     integer, dimension(n + 1), intent(In) :: ia
     integer, dimension(*), intent(Out) :: jlu
     integer, dimension(n), intent(InOut) :: ju
     integer, dimension(2 * n), intent(Out) :: jw
     integer, intent(Out) :: ierr
     !      locals
     integer :: ju0, k, j1, j2, j, ii, i, lenl, lenu, jj, jrow, jpos, len
     real(kind=8) :: tnorm, abs, s, fact

! ----------------------------------------------------------------------*
!                       *** ILUT preconditioner ***                     *
!       incomplete LU factorization with dual truncation mechanism      *
! ----------------------------------------------------------------------*
!      Author: Yousef Saad *May, 5, 1990, Latest revision, August 1996  *
! ----------------------------------------------------------------------*
!  PARAMETERS                                                           
! -----------                                                           
! 
!  on entry:
! ========== 
!  n       = integer. The row dimension of the matrix A. The matrix 
! 
!  a,ja,ia = matrix stored in Compressed Sparse Row format.              
! 
!  lfil    = integer. The fill-in parameter. Each row of L and each row
!            of U will have a maximum of lfil elements (excluding the 
!            diagonal element). lfil must be .ge. 0.
!            ** WARNING: THE MEANING OF LFIL HAS CHANGED WITH RESPECT TO
!            EARLIER VERSIONS. 
! 
!  droptol = real*8. Sets the threshold for dropping small terms in the
!            factorization. See below for details on dropping strategy.
! 
!   
!  iwk     = integer. The lengths of arrays alu and jlu. If the arrays
!            are not big enough to store the ILU factorizations, ilut
!            will stop with an error message. 
! 
!  On return:
! ===========
! 
!  alu,jlu = matrix stored in Modified Sparse Row (MSR) format containing
!            the L and U factors together. The diagonal (stored in
!            alu(1:n) ) is inverted. Each i-th row of the alu,jlu matrix
!            contains the i-th row of L (excluding the diagonal entry=1)
!            followed by the i-th row of U.
! 
!  ju      = integer array of length n containing the pointers to
!            the beginning of each row of U in the matrix alu,jlu.
! 
!  ierr    = integer. Error message with the following meaning.
!            ierr  = 0    --> successful return.
!            ierr .gt. 0  --> zero pivot encountered at step number ierr.
!            ierr  = -1   --> Error. input matrix may be wrong.
!                             (The elimination process has generated a
!                             row in L or U whose length is .gt.  n.)
!            ierr  = -2   --> The matrix L overflows the array al.
!            ierr  = -3   --> The matrix U overflows the array alu.
!            ierr  = -4   --> Illegal value for lfil.
!            ierr  = -5   --> zero row encountered.
! 
!  work arrays:
! =============
!  jw      = integer work array of length 2*n.
!  w       = real work array of length n+1.
!   
! ----------------------------------------------------------------------
!  w, ju (1:n) store the working array [1:ii-1 = L-part, ii:n = u] 
!  jw(n+1:2n)  stores nonzero indicators
!  
!  Notes:
!  ------
!  The diagonal elements of the input matrix must be  nonzero (at least
!  'structurally'). 
! 
! ----------------------------------------------------------------------* 
! ---- Dual drop strategy works as follows.                             *
!                                                                       *
!      1) Theresholding in L and U as set by droptol. Any element whose *
!         magnitude is less than some tolerance (relative to the abs    *
!         value of diagonal element in u) is dropped.                   *
!                                                                       *
!      2) Keeping only the largest lfil elements in the i-th row of L   * 
!         and the largest lfil elements in the i-th row of U (excluding *
!         diagonal elements).                                           *
!                                                                       *
!  Flexibility: one  can use  droptol=0  to get  a strategy  based on   *
!  keeping  the largest  elements in  each row  of L  and U.   Taking   *
!  droptol .ne.  0 but lfil=n will give  the usual threshold strategy   *
!  (however, fill-in is then mpredictible).                             *
! ----------------------------------------------------------------------*

     if (lfil < 0) then
          ierr = -4
          return
     end if
! -----------------------------------------------------------------------
!      initialize ju0 (points to next element to be added to alu,jlu)
!      and pointer array.
! -----------------------------------------------------------------------
     ju0 = n+2
     jlu(1) = ju0
! 
!      initialize nonzero indicator array. 
! 
     do j=1,n
          jw(n+j) = 0
     end do
! -----------------------------------------------------------------------
!      beginning of main loop.
! -----------------------------------------------------------------------
     do ii = 1, n
          j1 = ia(ii)
          j2 = ia(ii+1) - 1
          tnorm = 0.0d0
          do k=j1,j2
               tnorm = tnorm+abs(a(k))
          end do
          if (tnorm == 0.0) then
               ierr = -5
               return
          end if
          tnorm = tnorm/real(j2-j1+1)
!      
!      unpack L-part and U-part of row of A in arrays w 
!      
          lenu = 1
          lenl = 0
          jw(ii) = ii
          w(ii) = 0.0
          jw(n+ii) = ii
! 
          do j = j1, j2
               k = ja(j)
               t = a(j)
               if (k < ii) then
                    lenl = lenl+1
                    jw(lenl) = k
                    w(lenl) = t
                    jw(n+k) = lenl
               else if (k == ii) then
                    w(ii) = t
               else
                    lenu = lenu+1
                    jpos = ii+lenu-1 
                    jw(jpos) = k
                    w(jpos) = t
                    jw(n+k) = jpos
               endif
          end do
          jj = 0
          len = 0 
!      
!      eliminate previous rows
!      
 150 jj = jj+1
 if (jj > lenl) goto 160
! -----------------------------------------------------------------------
!      in order to do the elimination in the correct order we must select
!      the smallest column index among jw(k), k=jj+1, ..., lenl.
! -----------------------------------------------------------------------
 jrow = jw(jj)
 k = jj
!      
!      determine smallest column index
!      
 do j=jj+1,lenl
 if (jw(j) < jrow) then
 jrow = jw(j)
 k = j
 endif
 end do
! 
 if (k /= jj) then
!      exchange in jw
 j = jw(jj)
 jw(jj) = jw(k)
 jw(k) = j
!      exchange in jr
 jw(n+jrow) = jj
 jw(n+j) = k
!      exchange in w
 s = w(jj)
 w(jj) = w(k)
 w(k) = s
 endif
! 
!      zero out element in row by setting jw(n+jrow) to zero.
!      
 jw(n+jrow) = 0
! 
!      get the multiplier for row to be eliminated (jrow).
!      
 fact = w(jj)*alu(jrow)
 if (abs(fact) <= droptol) goto 150
!      
!      combine current row and row jrow
! 
 do k = ju(jrow), jlu(jrow+1)-1
 s = fact*alu(k)
 j = jlu(k)
 jpos = jw(n+j)
 if (j >= ii) then
!      
!      dealing with upper part.
!      
 if (jpos == 0) then
! 
!      this is a fill-in element
!      
 lenu = lenu+1
 if (lenu > n) goto 995
 i = ii+lenu-1
 jw(i) = j
 jw(n+j) = i
 w(i) = - s
 else
! 
!      this is not a fill-in element 
! 
 w(jpos) = w(jpos) - s
 endif
 else
!      
!      dealing  with lower part.
!      
 if (jpos == 0) then
! 
!      this is a fill-in element
!      
 lenl = lenl+1
 if (lenl > n) goto 995
 jw(lenl) = j
 jw(n+j) = lenl
 w(lenl) = - s
 else
!      
!      this is not a fill-in element 
!      
 w(jpos) = w(jpos) - s
 endif
 endif
 end do
!      
!      store this pivot element -- (from left to right -- no danger of
!      overlap with the working elements in L (pivots). 
!      
 len = len+1 
 w(len) = fact
 jw(len) = jrow
 goto 150
 160 continue
!      
!      reset double-pointer to zero (U-part)
!      
 do k=1, lenu
 jw(n+jw(ii+k-1)) = 0
 end do
!      
!      update L-matrix
!      
 lenl = len 
 len = min0(lenl,lfil)
!      
!      sort by quick-split
! 
 w_qsplit_cr = reshape(w,shape(w_qsplit_cr))
 call qsplit(w_qsplit_cr,jw,lenl,len) 
 w = reshape(w_qsplit_cr, shape(w))
! 
!      store L-part
!  
 do k=1, len 
 if (ju0 > iwk) goto 996
 alu(ju0) = w(k)
 jlu(ju0) = jw(k)
 ju0 = ju0+1
 end do
!      
!      save pointer to beginning of row ii of U
!      
 ju(ii) = ju0
! 
!      update U-matrix -- first apply dropping strategy 
! 
 len = 0
 do k=1, lenu-1
 if (abs(w(ii+k)) > droptol*tnorm) then 
 len = len+1
 w(ii+len) = w(ii+k) 
 jw(ii+len) = jw(ii+k) 
 endif
 enddo
 lenu = len+1
 len = min0(lenu,lfil)
! 
 call qsplit(w(ii + 1),jw(ii + 1),lenu - 1,len) 
! 
!      copy
!  
 t = abs(w(ii))
 if (len + ju0 > iwk) goto 997
 do k=ii+1,ii+len-1 
 jlu(ju0) = jw(k)
 alu(ju0) = w(k)
 t = t + abs(w(k) )
 ju0 = ju0+1
 end do
!      
!      store inverse of diagonal element of u
!      
 if (w(ii) == 0.0) w(ii) = (0.0001 + droptol)*tnorm
!      
 alu(ii) = 1.0d0/ w(ii) 
!      
!      update pointer to beginning of next row of U.
!      
 jlu(ii+1) = ju0
! -----------------------------------------------------------------------
!      end main loop
! -----------------------------------------------------------------------
 end do
 ierr = 0
 return
! 
!      incomprehensible error. Matrix must be wrong.
!      
 995 ierr = -1
 return
!      
!      insufficient storage in L.
!      
 996 ierr = -2
 return
!      
!      insufficient storage in U.
!      
 997 ierr = -3
 return
!      
!      illegal lfil entered.
!      
 998 ierr = -4
 return
!      
!      zero row encountered
!      
 999 ierr = -5
 return
! ----------------end-of-ilut--------------------------------------------
! -----------------------------------------------------------------------
 end subroutine ilut
 subroutine ilutp(n,a,ja,ia,lfil,droptol,permtol,mbloc,alu,jlu,ju,iwk,w,jw,iperm,ierr)
! BEGIN new declarations
! END new declarations
 implicit none
 integer, intent(In) :: n
 integer, dimension(1:*), intent(Out) :: ja
 integer, dimension(1:n + 1), intent(In) :: ia
 integer, intent(In) :: lfil
 integer, dimension(1:*), intent(Out) :: jlu
 integer, dimension(1:n), intent(InOut) :: ju
 integer, dimension(1:2 * n), intent(Out) :: jw
 integer, intent(In) :: iwk
 integer, dimension(1:2 * n), intent(Out) :: iperm
 integer, intent(Out) :: ierr
 real(kind=8), dimension(1:*), intent(In) :: a
 real(kind=8), dimension(1:*), intent(InOut) :: alu
 real(kind=8), dimension(1:n + 1), intent(Out) :: w
 real(kind=8), intent(In) :: droptol
! -----------------------------------------------------------------------
!      implicit none
! ----------------------------------------------------------------------*
!        *** ILUTP preconditioner -- ILUT with pivoting  ***            *
!       incomplete LU factorization with dual truncation mechanism      *
! ----------------------------------------------------------------------*
!  author Yousef Saad *Sep 8, 1993 -- Latest revision, August 1996.     *
! ----------------------------------------------------------------------*
!  on entry:
! ==========
!  n       = integer. The dimension of the matrix A.
! 
!  a,ja,ia = matrix stored in Compressed Sparse Row format.
!            ON RETURN THE COLUMNS OF A ARE PERMUTED. SEE BELOW FOR 
!            DETAILS. 
! 
!  lfil    = integer. The fill-in parameter. Each row of L and each row
!            of U will have a maximum of lfil elements (excluding the 
!            diagonal element). lfil must be .ge. 0.
!            ** WARNING: THE MEANING OF LFIL HAS CHANGED WITH RESPECT TO
!            EARLIER VERSIONS. 
! 
!  droptol = real*8. Sets the threshold for dropping small terms in the
!            factorization. See below for details on dropping strategy.
! 
!  lfil    = integer. The fill-in parameter. Each row of L and
!            each row of U will have a maximum of lfil elements.
!            WARNING: THE MEANING OF LFIL HAS CHANGED WITH RESPECT TO
!            EARLIER VERSIONS. 
!            lfil must be .ge. 0.
! 
!  permtol = tolerance ratio used to  determne whether or not to permute
!            two columns.  At step i columns i and j are permuted when 
! 
!                      abs(a(i,j))*permtol .gt. abs(a(i,i))
! 
!            [0 --> never permute; good values 0.1 to 0.01]
! 
!  mbloc   = if desired, permuting can be done only within the diagonal
!            blocks of size mbloc. Useful for PDE problems with several
!            degrees of freedom.. If feature not wanted take mbloc=n.
! 
!   
!  iwk     = integer. The lengths of arrays alu and jlu. If the arrays
!            are not big enough to store the ILU factorizations, ilut
!            will stop with an error message. 
! 
!  On return:
! ===========
! 
!  alu,jlu = matrix stored in Modified Sparse Row (MSR) format containing
!            the L and U factors together. The diagonal (stored in
!            alu(1:n) ) is inverted. Each i-th row of the alu,jlu matrix
!            contains the i-th row of L (excluding the diagonal entry=1)
!            followed by the i-th row of U.
! 
!  ju      = integer array of length n containing the pointers to
!            the beginning of each row of U in the matrix alu,jlu.
! 
!  iperm   = contains the permutation arrays. 
!            iperm(1:n) = old numbers of unknowns
!            iperm(n+1:2*n) = reverse permutation = new unknowns.
! 
!  ierr    = integer. Error message with the following meaning.
!            ierr  = 0    --> successful return.
!            ierr .gt. 0  --> zero pivot encountered at step number ierr.
!            ierr  = -1   --> Error. input matrix may be wrong.
!                             (The elimination process has generated a
!                             row in L or U whose length is .gt.  n.)
!            ierr  = -2   --> The matrix L overflows the array al.
!            ierr  = -3   --> The matrix U overflows the array alu.
!            ierr  = -4   --> Illegal value for lfil.
!            ierr  = -5   --> zero row encountered.
! 
!  work arrays:
! =============
!  jw      = integer work array of length 2*n.
!  w       = real work array of length n 
! 
!  IMPORTANR NOTE:
!  --------------
!  TO AVOID PERMUTING THE SOLUTION VECTORS ARRAYS FOR EACH LU-SOLVE, 
!  THE MATRIX A IS PERMUTED ON RETURN. [all column indices are
!  changed]. SIMILARLY FOR THE U MATRIX. 
!  To permute the matrix back to its original state use the loop:
! 
!       do k=ia(1), ia(n+1)-1
!          ja(k) = iperm(ja(k)) 
!       enddo
!  
! -----------------------------------------------------------------------
!      local variables
! 
 integer :: k
 integer :: i
 integer :: j
 integer :: jrow
 integer :: ju0
 integer :: ii
 integer :: j1
 integer :: j2
 integer :: jpos
 integer :: len
 integer :: imax
 integer :: lenu
 integer :: lenl
 integer :: jj
 integer, intent(In) :: mbloc
 integer :: icut
 real(kind=8) :: s
 real(kind=8) :: tmp
 real(kind=8) :: tnorm
 real(kind=8) :: xmax
 real(kind=8) :: xmax0
 real(kind=8) :: fact
 real(kind=8) :: abs
 real(kind=8) :: t
 real(kind=8), intent(In) :: permtol
!      
 real(kind=8), dimension(1:n + 1) :: w_qsplit_cr
 if (lfil < 0) goto 998
! ----------------------------------------------------------------------- 
!      initialize ju0 (points to next element to be added to alu,jlu)
!      and pointer array.
! -----------------------------------------------------------------------
 ju0 = n+2
 jlu(1) = ju0
! 
!   integer double pointer array.
! 
 do j=1, n
 jw(n+j) = 0
 iperm(j) = j
 iperm(n+j) = j
 end do
! -----------------------------------------------------------------------
!      beginning of main loop.
! -----------------------------------------------------------------------
 do ii = 1, n
 j1 = ia(ii)
 j2 = ia(ii+1) - 1
 tnorm = 0.0d0
 do k=j1,j2
 tnorm = tnorm+abs(a(k))
 end do
 if (tnorm == 0.0) goto 999
 tnorm = tnorm/(j2-j1+1)
! 
!      unpack L-part and U-part of row of A in arrays  w  --
! 
 lenu = 1
 lenl = 0
 jw(ii) = ii
 w(ii) = 0.0
 jw(n+ii) = ii
! 
 do j = j1, j2
 k = iperm(n+ja(j))
 t = a(j)
 if (k < ii) then
 lenl = lenl+1
 jw(lenl) = k
 w(lenl) = t
 jw(n+k) = lenl
 else if (k == ii) then
 w(ii) = t
 else
 lenu = lenu+1
 jpos = ii+lenu-1 
 jw(jpos) = k
 w(jpos) = t
 jw(n+k) = jpos
 endif
 end do
 jj = 0
 len = 0 
! 
!      eliminate previous rows
! 
 150 jj = jj+1
 if (jj > lenl) goto 160
! -----------------------------------------------------------------------
!      in order to do the elimination in the correct order we must select
!      the smallest column index among jw(k), k=jj+1, ..., lenl.
! -----------------------------------------------------------------------
 jrow = jw(jj)
 k = jj
! 
!      determine smallest column index
! 
 do j=jj+1,lenl
 if (jw(j) < jrow) then
 jrow = jw(j)
 k = j
 endif
 end do
! 
 if (k /= jj) then
!      exchange in jw
 j = jw(jj)
 jw(jj) = jw(k)
 jw(k) = j
!      exchange in jr
 jw(n+jrow) = jj
 jw(n+j) = k
!      exchange in w
 s = w(jj)
 w(jj) = w(k)
 w(k) = s
 endif
! 
!      zero out element in row by resetting jw(n+jrow) to zero.
!      
 jw(n+jrow) = 0
! 
!      get the multiplier for row to be eliminated: jrow
! 
 fact = w(jj)*alu(jrow)
! 
!      drop term if small
!      
 if (abs(fact) <= droptol) goto 150
! 
!      combine current row and row jrow
! 
 do k = ju(jrow), jlu(jrow+1)-1
 s = fact*alu(k)
!      new column number
 j = iperm(n+jlu(k))
 jpos = jw(n+j)
 if (j >= ii) then
! 
!      dealing with upper part.
! 
 if (jpos == 0) then
! 
!      this is a fill-in element
! 
 lenu = lenu+1
 i = ii+lenu-1 
 if (lenu > n) goto 995
 jw(i) = j
 jw(n+j) = i 
 w(i) = - s
 else
!      no fill-in element --
 w(jpos) = w(jpos) - s
 endif
 else
! 
!      dealing with lower part.
! 
 if (jpos == 0) then
! 
!      this is a fill-in element
! 
 lenl = lenl+1
 if (lenl > n) goto 995
 jw(lenl) = j
 jw(n+j) = lenl
 w(lenl) = - s
 else
! 
!      this is not a fill-in element
! 
 w(jpos) = w(jpos) - s
 endif
 endif
 end do
!      
!      store this pivot element -- (from left to right -- no danger of
!      overlap with the working elements in L (pivots). 
!      
 len = len+1 
 w(len) = fact
 jw(len) = jrow
 goto 150
 160 continue
! 
!      reset double-pointer to zero (U-part)
!      
 do k=1, lenu
 jw(n+jw(ii+k-1)) = 0
 end do
! 
!      update L-matrix
! 
 lenl = len 
 len = min0(lenl,lfil)
!      
!      sort by quick-split
! 
 w_qsplit_cr = reshape(w,shape(w_qsplit_cr))
 call qsplit(w_qsplit_cr,jw,lenl,len) 
 w = reshape(w_qsplit_cr, shape(w))
! 
!      store L-part -- in original coordinates ..
! 
 do k=1, len
 if (ju0 > iwk) goto 996
 alu(ju0) = w(k) 
 jlu(ju0) = iperm(jw(k))
 ju0 = ju0+1
 end do
! 
!      save pointer to beginning of row ii of U
! 
 ju(ii) = ju0
! 
!      update U-matrix -- first apply dropping strategy 
! 
 len = 0
 do k=1, lenu-1
 if (abs(w(ii+k)) > droptol*tnorm) then 
 len = len+1
 w(ii+len) = w(ii+k) 
 jw(ii+len) = jw(ii+k) 
 endif
 enddo
 lenu = len+1
 len = min0(lenu,lfil)
 call qsplit(w(ii + 1),jw(ii + 1),lenu - 1,len) 
! 
!      determine next pivot -- 
! 
 imax = ii
 xmax = abs(w(imax))
 xmax0 = xmax
 icut = ii - 1 + mbloc - mod(ii-1,mbloc)
 do k=ii+1,ii+len-1
 t = abs(w(k))
 if (t > xmax .and. t*permtol > xmax0 .and. jw(k) <= icut) then
 imax = k
 xmax = t
 endif
 enddo
! 
!      exchange w's
! 
 tmp = w(ii)
 w(ii) = w(imax)
 w(imax) = tmp
! 
!      update iperm and reverse iperm
! 
 j = jw(imax)
 i = iperm(ii)
 iperm(ii) = iperm(j)
 iperm(j) = i
! 
!      reverse iperm
! 
 iperm(n+iperm(ii)) = ii
 iperm(n+iperm(j)) = j
! ----------------------------------------------------------------------- 
! 
 if (len + ju0 > iwk) goto 997
! 
!      copy U-part in original coordinates
!      
 do k=ii+1,ii+len-1 
 jlu(ju0) = iperm(jw(k))
 alu(ju0) = w(k)
 ju0 = ju0+1
 end do
! 
!      store inverse of diagonal element of u
! 
 if (w(ii) == 0.0) w(ii) = (1.0d-4 + droptol)*tnorm
 alu(ii) = 1.0d0/ w(ii) 
! 
!      update pointer to beginning of next row of U.
! 
 jlu(ii+1) = ju0
! -----------------------------------------------------------------------
!      end main loop
! -----------------------------------------------------------------------
 end do
! 
!      permute all column indices of LU ...
! 
 do k = jlu(1),jlu(n+1)-1
 jlu(k) = iperm(n+jlu(k))
 enddo
! 
!      ...and of A
! 
 do k=ia(1), ia(n+1)-1
 ja(k) = iperm(n+ja(k))
 enddo
! 
 ierr = 0
 return
! 
!      incomprehensible error. Matrix must be wrong.
! 
 995 ierr = -1
 return
! 
!      insufficient storage in L.
! 
 996 ierr = -2
 return
! 
!      insufficient storage in U.
! 
 997 ierr = -3
 return
! 
!      illegal lfil entered.
! 
 998 ierr = -4
 return
! 
!      zero row encountered
! 
 999 ierr = -5
 return
! ----------------end-of-ilutp-------------------------------------------
! -----------------------------------------------------------------------
 end subroutine ilutp
 subroutine ilud(n,a,ja,ia,alph,tol,alu,jlu,ju,iwk,w,jw,ierr)
! -----------------------------------------------------------------------
 implicit none 
! BEGIN new declarations
! END new declarations
 integer, intent(In) :: n
 real(kind=8), dimension(1:*), intent(In) :: a
 real(kind=8), dimension(1:*), intent(InOut) :: alu
 real(kind=8), dimension(1:2 * n), intent(Out) :: w
 real(kind=8), intent(In) :: tol
 real(kind=8), intent(In) :: alph
 integer, dimension(1:*), intent(In) :: ja
 integer, dimension(1:n + 1), intent(In) :: ia
 integer, dimension(1:*), intent(Out) :: jlu
 integer, dimension(1:n), intent(InOut) :: ju
 integer, dimension(1:2 * n), intent(Out) :: jw
 integer, intent(In) :: iwk
 integer, intent(Out) :: ierr
! ----------------------------------------------------------------------*
!                      *** ILUD preconditioner ***                      *
!     incomplete LU factorization with standard droppoing strategy      *
! ----------------------------------------------------------------------*
!  Author: Yousef Saad * Aug. 1995 --                                   * 
! ----------------------------------------------------------------------*
!  This routine computes the ILU factorization with standard threshold  *
!  dropping: at i-th step of elimination, an element a(i,j) in row i is *
!  dropped  if it satisfies the criterion:                              *
!                                                                       *
!   abs(a(i,j)) < tol * [average magnitude of elements in row i of A]   *
!                                                                       *
!  There is no control on memory size required for the factors as is    *
!  done in ILUT. This routines computes also various diagonal compensa- * 
!  tion ILU's such MILU. These are defined through the parameter alph   *
! ----------------------------------------------------------------------* 
!  on entry:
! ========== 
!  n       = integer. The row dimension of the matrix A. The matrix 
! 
!  a,ja,ia = matrix stored in Compressed Sparse Row format              
! 
!  alph    = diagonal compensation parameter -- the term: 
! 
!            alph*(sum of all dropped out elements in a given row) 
! 
!            is added to the diagonal element of U of the factorization 
!            Thus: alph = 0 ---> ~ ILU with threshold,
!                  alph = 1 ---> ~ MILU with threshold. 
!  
!  tol     = Threshold parameter for dropping small terms in the
!            factorization. During the elimination, a term a(i,j) is 
!            dropped whenever abs(a(i,j)) .lt. tol * [weighted norm of
!            row i]. Here weighted norm = 1-norm / number of nnz 
!            elements in the row. 
!   
!  iwk     = The length of arrays alu and jlu -- this routine will stop
!            if storage for the factors L and U is not sufficient 
! 
!  On return:
! =========== 
! 
!  alu,jlu = matrix stored in Modified Sparse Row (MSR) format containing
!            the L and U factors together. The diagonal (stored in
!            alu(1:n) ) is inverted. Each i-th row of the alu,jlu matrix
!            contains the i-th row of L (excluding the diagonal entry=1)
!            followed by the i-th row of U.
! 
!  ju      = integer array of length n containing the pointers to
!            the beginning of each row of U in the matrix alu,jlu.
! 
!  ierr    = integer. Error message with the following meaning.
!            ierr  = 0    --> successful return.
!            ierr .gt. 0  --> zero pivot encountered at step number ierr.
!            ierr  = -1   --> Error. input matrix may be wrong.
!                             (The elimination process has generated a
!                             row in L or U whose length is .gt.  n.)
!            ierr  = -2   --> Insufficient storage for the LU factors --
!                             arrays alu/ jalu are  overflowed. 
!            ierr  = -3   --> Zero row encountered.
! 
!  Work Arrays:
! =============
!  jw      = integer work array of length 2*n.
!  w       = real work array of length n 
!   
! ----------------------------------------------------------------------
! 
!  w, ju (1:n) store the working array [1:ii-1 = L-part, ii:n = u] 
!  jw(n+1:2n)  stores the nonzero indicator. 
!  
!  Notes:
!  ------
!  All diagonal elements of the input matrix must be  nonzero.
! 
! ----------------------------------------------------------------------- 
!      locals
 integer :: ju0
 integer :: k
 integer :: j1
 integer :: j2
 integer :: j
 integer :: ii
 integer :: i
 integer :: lenl
 integer :: lenu
 integer :: jj
 integer :: jrow
 integer :: jpos
 integer :: len
 real(kind=8) :: tnorm
 real(kind=8) :: t
 real(kind=8) :: abs
 real(kind=8) :: s
 real(kind=8) :: fact
 real(kind=8) :: dropsum
! -----------------------------------------------------------------------
!      initialize ju0 (points to next element to be added to alu,jlu)
!      and pointer array.
! -----------------------------------------------------------------------
 ju0 = n+2
 jlu(1) = ju0
! 
!      initialize nonzero indicator array. 
! 
 do j=1,n
 jw(n+j) = 0
 end do
! -----------------------------------------------------------------------
!      beginning of main loop.
! -----------------------------------------------------------------------
 do ii = 1, n
 j1 = ia(ii)
 j2 = ia(ii+1) - 1
 dropsum = 0.0d0 
 tnorm = 0.0d0
 do k=j1,j2
 tnorm = tnorm + abs(a(k)) 
 end do
 if (tnorm == 0.0) goto 997
 tnorm = tnorm / real(j2-j1+1) 
!      
!      unpack L-part and U-part of row of A in arrays w 
!      
 lenu = 1
 lenl = 0
 jw(ii) = ii
 w(ii) = 0.0
 jw(n+ii) = ii
! 
 do j = j1, j2
 k = ja(j)
 t = a(j)
 if (k < ii) then
 lenl = lenl+1
 jw(lenl) = k
 w(lenl) = t
 jw(n+k) = lenl
 else if (k == ii) then
 w(ii) = t
 else
 lenu = lenu+1
 jpos = ii+lenu-1 
 jw(jpos) = k
 w(jpos) = t
 jw(n+k) = jpos
 endif
 end do
 jj = 0
 len = 0 
!      
!      eliminate previous rows
!      
 150 jj = jj+1
 if (jj > lenl) goto 160
! -----------------------------------------------------------------------
!      in order to do the elimination in the correct order we must select
!      the smallest column index among jw(k), k=jj+1, ..., lenl.
! -----------------------------------------------------------------------
 jrow = jw(jj)
 k = jj
!      
!      determine smallest column index
!      
 do j=jj+1,lenl
 if (jw(j) < jrow) then
 jrow = jw(j)
 k = j
 endif
 end do
! 
 if (k /= jj) then
!      exchange in jw
 j = jw(jj)
 jw(jj) = jw(k)
 jw(k) = j
!      exchange in jr
 jw(n+jrow) = jj
 jw(n+j) = k
!      exchange in w
 s = w(jj)
 w(jj) = w(k)
 w(k) = s
 endif
! 
!      zero out element in row by setting resetting jw(n+jrow) to zero.
!      
 jw(n+jrow) = 0
! 
!      drop term if small
!      
!          if (abs(w(jj)) .le. tol*tnorm) then
!             dropsum = dropsum + w(jj) 
!             goto 150
!          endif
!      
!      get the multiplier for row to be eliminated (jrow).
!      
 fact = w(jj)*alu(jrow)
! 
!      drop term if small
!      
 if (abs(fact) <= tol) then
 dropsum = dropsum + w(jj) 
 goto 150
 endif
!      
!      combine current row and row jrow
! 
 do k = ju(jrow), jlu(jrow+1)-1
 s = fact*alu(k)
 j = jlu(k)
 jpos = jw(n+j)
 if (j >= ii) then
!      
!      dealing with upper part.
!      
 if (jpos == 0) then
! 
!      this is a fill-in element
!      
 lenu = lenu+1
 if (lenu > n) goto 995
 i = ii+lenu-1
 jw(i) = j
 jw(n+j) = i
 w(i) = - s
 else
! 
!      this is not a fill-in element 
! 
 w(jpos) = w(jpos) - s
 endif
 else
!      
!      dealing with lower part.
!      
 if (jpos == 0) then
! 
!      this is a fill-in element
! 
 lenl = lenl+1
 if (lenl > n) goto 995
 jw(lenl) = j
 jw(n+j) = lenl
 w(lenl) = - s
 else
! 
!      this is not a fill-in element 
! 
 w(jpos) = w(jpos) - s
 endif
 endif
 end do
 len = len+1 
 w(len) = fact
 jw(len) = jrow
 goto 150
 160 continue
!      
!      reset double-pointer to zero (For U-part only)
!      
 do k=1, lenu
 jw(n+jw(ii+k-1)) = 0
 end do
! 
!      update l-matrix
! 
 do k=1, len
 if (ju0 > iwk) goto 996
 alu(ju0) = w(k) 
 jlu(ju0) = jw(k)
 ju0 = ju0+1
 end do
!      
!      save pointer to beginning of row ii of U
!      
 ju(ii) = ju0
! 
!      go through elements in U-part of w to determine elements to keep
! 
 len = 0
 do k=1, lenu-1
!             if (abs(w(ii+k)) .gt. tnorm*tol) then 
 if (abs(w(ii+k)) > abs(w(ii))*tol) then 
 len = len+1
 w(ii+len) = w(ii+k) 
 jw(ii+len) = jw(ii+k)
 else
 dropsum = dropsum + w(ii+k) 
 endif
 enddo
! 
!      now update u-matrix
! 
 if (ju0 + len-1 > iwk) goto 996
 do k=ii+1,ii+len
 jlu(ju0) = jw(k)
 alu(ju0) = w(k)
 ju0 = ju0+1
 end do
! 
!      define diagonal element 
!  
 w(ii) = w(ii) + alph*dropsum 
! 
!      store inverse of diagonal element of u
!               
 if (w(ii) == 0.0) w(ii) = (0.0001 + tol)*tnorm
!      
 alu(ii) = 1.0d0/ w(ii) 
!      
!      update pointer to beginning of next row of U.
!      
 jlu(ii+1) = ju0
! -----------------------------------------------------------------------
!      end main loop
! -----------------------------------------------------------------------
 end do
 ierr = 0
 return
! 
!      incomprehensible error. Matrix must be wrong.
!      
 995 ierr = -1
 return
!      
!      insufficient storage in alu/ jlu arrays for  L / U factors 
!      
 996 ierr = -2
 return
!      
!      zero row encountered
!      
 997 ierr = -3 
 return
! ----------------end-of-ilud  ------------------------------------------
! -----------------------------------------------------------------------
 end subroutine ilud
 subroutine iludp(n,a,ja,ia,alph,droptol,permtol,mbloc,alu,jlu,ju,iwk,w,jw,iperm,ierr)
! -----------------------------------------------------------------------
 implicit none
! BEGIN new declarations
! END new declarations
 integer, intent(In) :: n
 integer, dimension(1:*), intent(Out) :: ja
 integer, dimension(1:n + 1), intent(In) :: ia
 integer, intent(In) :: mbloc
 integer, dimension(1:*), intent(Out) :: jlu
 integer, dimension(1:n), intent(InOut) :: ju
 integer, dimension(1:2 * n), intent(Out) :: jw
 integer, intent(In) :: iwk
 integer, dimension(1:2 * n), intent(Out) :: iperm
 integer, intent(Out) :: ierr
 real(kind=8), dimension(1:*), intent(In) :: a
 real(kind=8), dimension(1:*), intent(InOut) :: alu
 real(kind=8), dimension(1:2 * n), intent(Out) :: w
 real(kind=8), intent(In) :: alph
 real(kind=8), intent(In) :: droptol
 real(kind=8), intent(In) :: permtol
! ----------------------------------------------------------------------*
!                      *** ILUDP preconditioner ***                     *
!     incomplete LU factorization with standard droppoing strategy      *
!     and column pivoting                                               * 
! ----------------------------------------------------------------------*
!  author Yousef Saad -- Aug 1995.                                      *
! ----------------------------------------------------------------------*
!  on entry:
! ==========
!  n       = integer. The dimension of the matrix A.
! 
!  a,ja,ia = matrix stored in Compressed Sparse Row format.
!            ON RETURN THE COLUMNS OF A ARE PERMUTED.
! 
!  alph    = diagonal compensation parameter -- the term: 
! 
!            alph*(sum of all dropped out elements in a given row) 
! 
!            is added to the diagonal element of U of the factorization 
!            Thus: alph = 0 ---> ~ ILU with threshold,
!                  alph = 1 ---> ~ MILU with threshold. 
!  
!  droptol = tolerance used for dropping elements in L and U.
!            elements are dropped if they are .lt. norm(row) x droptol
!            row = row being eliminated
! 
!  permtol = tolerance ratio used for determning whether to permute
!            two columns.  Two columns are permuted only when 
!            abs(a(i,j))*permtol .gt. abs(a(i,i))
!            [0 --> never permute; good values 0.1 to 0.01]
! 
!  mbloc   = if desired, permuting can be done only within the diagonal
!            blocks of size mbloc. Useful for PDE problems with several
!            degrees of freedom.. If feature not wanted take mbloc=n.
! 
!  iwk     = integer. The declared lengths of arrays alu and jlu
!            if iwk is not large enough the code will stop prematurely
!            with ierr = -2 or ierr = -3 (see below).
! 
!  On return:
! ===========
! 
!  alu,jlu = matrix stored in Modified Sparse Row (MSR) format containing
!            the L and U factors together. The diagonal (stored in
!            alu(1:n) ) is inverted. Each i-th row of the alu,jlu matrix
!            contains the i-th row of L (excluding the diagonal entry=1)
!            followed by the i-th row of U.
! 
!  ju      = integer array of length n containing the pointers to
!            the beginning of each row of U in the matrix alu,jlu.
!  iperm   = contains the permutation arrays ..
!            iperm(1:n) = old numbers of unknowns
!            iperm(n+1:2*n) = reverse permutation = new unknowns.
! 
!  ierr    = integer. Error message with the following meaning.
!            ierr  = 0    --> successful return.
!            ierr .gt. 0  --> zero pivot encountered at step number ierr.
!            ierr  = -1   --> Error. input matrix may be wrong.
!                             (The elimination process has generated a
!                             row in L or U whose length is .gt.  n.)
!            ierr  = -2   --> The L/U matrix overflows the arrays alu,jlu
!            ierr  = -3   --> zero row encountered.
! 
!  work arrays:
! =============
!  jw      = integer work array of length 2*n.
!  w       = real work array of length 2*n 
! 
!  Notes:
!  ------
!  IMPORTANT: TO AVOID PERMUTING THE SOLUTION VECTORS ARRAYS FOR EACH 
!  LU-SOLVE, THE MATRIX A IS PERMUTED ON RETURN. [all column indices are
!  changed]. SIMILARLY FOR THE U MATRIX. 
!  To permute the matrix back to its original state use the loop:
! 
!       do k=ia(1), ia(n+1)-1
!          ja(k) = perm(ja(k)) 
!       enddo
!  
! -----------------------------------------------------------------------
!      local variables
! 
 integer :: k
 integer :: i
 integer :: j
 integer :: jrow
 integer :: ju0
 integer :: ii
 integer :: j1
 integer :: j2
 integer :: jpos
 integer :: len
 integer :: imax
 integer :: lenu
 integer :: lenl
 integer :: jj
 integer :: icut
 real(kind=8) :: s
 real(kind=8) :: tmp
 real(kind=8) :: tnorm
 real(kind=8) :: xmax
 real(kind=8) :: xmax0
 real(kind=8) :: fact
 real(kind=8) :: abs
 real(kind=8) :: t
 real(kind=8) :: dropsum
! ----------------------------------------------------------------------- 
!      initialize ju0 (points to next element to be added to alu,jlu)
!      and pointer array.
! -----------------------------------------------------------------------
 ju0 = n+2
 jlu(1) = ju0
! 
!   integer double pointer array.
! 
 do j=1,n
 jw(n+j) = 0
 iperm(j) = j
 iperm(n+j) = j
 end do
! -----------------------------------------------------------------------
!      beginning of main loop.
! -----------------------------------------------------------------------
 do ii = 1, n
 j1 = ia(ii)
 j2 = ia(ii+1) - 1
 dropsum = 0.0d0 
 tnorm = 0.0d0
 do k=j1,j2
 tnorm = tnorm+abs(a(k))
 end do
 if (tnorm == 0.0) goto 997
 tnorm = tnorm/(j2-j1+1)
! 
!      unpack L-part and U-part of row of A in arrays  w  --
! 
 lenu = 1
 lenl = 0
 jw(ii) = ii
 w(ii) = 0.0
 jw(n+ii) = ii
! 
 do j = j1, j2
 k = iperm(n+ja(j))
 t = a(j)
 if (k < ii) then
 lenl = lenl+1
 jw(lenl) = k
 w(lenl) = t
 jw(n+k) = lenl
 else if (k == ii) then
 w(ii) = t
 else
 lenu = lenu+1
 jpos = ii+lenu-1 
 jw(jpos) = k
 w(jpos) = t
 jw(n+k) = jpos
 endif
 end do
 jj = 0
 len = 0 
! 
!      eliminate previous rows
! 
 150 jj = jj+1
 if (jj > lenl) goto 160
! -----------------------------------------------------------------------
!      in order to do the elimination in the correct order we must select
!      the smallest column index among jw(k), k=jj+1, ..., lenl.
! -----------------------------------------------------------------------
 jrow = jw(jj)
 k = jj
! 
!      determine smallest column index
! 
 do j=jj+1,lenl
 if (jw(j) < jrow) then
 jrow = jw(j)
 k = j
 endif
 end do
! 
 if (k /= jj) then
!      exchange in jw
 j = jw(jj)
 jw(jj) = jw(k)
 jw(k) = j
!      exchange in jr
 jw(n+jrow) = jj
 jw(n+j) = k
!      exchange in w
 s = w(jj)
 w(jj) = w(k)
 w(k) = s
 endif
! 
!      zero out element in row by resetting jw(n+jrow) to zero.
!      
 jw(n+jrow) = 0
! 
!      drop term if small
!      
 if (abs(w(jj)) <= droptol*tnorm) then
 dropsum = dropsum + w(jj) 
 goto 150
 endif 
! 
!      get the multiplier for row to be eliminated: jrow
! 
 fact = w(jj)*alu(jrow)
! 
!      combine current row and row jrow
! 
 do k = ju(jrow), jlu(jrow+1)-1
 s = fact*alu(k)
!      new column number
 j = iperm(n+jlu(k))
 jpos = jw(n+j)
! 
!      if fill-in element is small then disregard:
!      
 if (j >= ii) then
! 
!      dealing with upper part.
! 
 if (jpos == 0) then
!      this is a fill-in element
 lenu = lenu+1
 i = ii+lenu-1 
 if (lenu > n) goto 995
 jw(i) = j
 jw(n+j) = i 
 w(i) = - s
 else
!      no fill-in element --
 w(jpos) = w(jpos) - s
 endif
 else
! 
!      dealing with lower part.
! 
 if (jpos == 0) then
!      this is a fill-in element
 lenl = lenl+1
 if (lenl > n) goto 995
 jw(lenl) = j
 jw(n+j) = lenl
 w(lenl) = - s
 else
!      no fill-in element --
 w(jpos) = w(jpos) - s
 endif
 endif
 end do
 len = len+1 
 w(len) = fact
 jw(len) = jrow
 goto 150
 160 continue
! 
!      reset double-pointer to zero (U-part)
!      
 do k=1, lenu
 jw(n+jw(ii+k-1)) = 0
 end do
! 
!      update L-matrix
! 
 do k=1, len
 if (ju0 > iwk) goto 996
 alu(ju0) = w(k)
 jlu(ju0) = iperm(jw(k))
 ju0 = ju0+1
 end do
! 
!      save pointer to beginning of row ii of U
! 
 ju(ii) = ju0
! 
!      update u-matrix -- first apply dropping strategy 
! 
 len = 0
 do k=1, lenu-1
 if (abs(w(ii+k)) > tnorm*droptol) then 
 len = len+1
 w(ii+len) = w(ii+k) 
 jw(ii+len) = jw(ii+k) 
 else
 dropsum = dropsum + w(ii+k) 
 endif
 enddo
! 
 imax = ii
 xmax = abs(w(imax))
 xmax0 = xmax
 icut = ii - 1 + mbloc - mod(ii-1,mbloc)
! 
!      determine next pivot -- 
!  
 do k=ii+1,ii+len 
 t = abs(w(k))
 if (t > xmax .and. t*permtol > xmax0 .and. jw(k) <= icut) then
 imax = k
 xmax = t
 endif
 enddo
! 
!      exchange w's
! 
 tmp = w(ii)
 w(ii) = w(imax)
 w(imax) = tmp
! 
!      update iperm and reverse iperm
! 
 j = jw(imax)
 i = iperm(ii)
 iperm(ii) = iperm(j)
 iperm(j) = i
!      reverse iperm
 iperm(n+iperm(ii)) = ii
 iperm(n+iperm(j)) = j
! ----------------------------------------------------------------------- 
 if (len + ju0-1 > iwk) goto 996
! 
!      copy U-part in original coordinates
!      
 do k=ii+1,ii+len
 jlu(ju0) = iperm(jw(k))
 alu(ju0) = w(k)
 ju0 = ju0+1
 end do
! 
!      define diagonal element 
!  
 w(ii) = w(ii) + alph*dropsum 
! 
!      store inverse of diagonal element of u
! 
 if (w(ii) == 0.0) w(ii) = (1.0d-4 + droptol)*tnorm
! 
 alu(ii) = 1.0d0/ w(ii) 
! 
!      update pointer to beginning of next row of U.
! 
 jlu(ii+1) = ju0
! -----------------------------------------------------------------------
!      end main loop
! -----------------------------------------------------------------------
 end do
! 
!      permute all column indices of LU ...
! 
 do k = jlu(1),jlu(n+1)-1
 jlu(k) = iperm(n+jlu(k))
 enddo
! 
!      ...and of A
! 
 do k=ia(1), ia(n+1)-1
 ja(k) = iperm(n+ja(k))
 enddo
! 
 ierr = 0
 return
! 
!      incomprehensible error. Matrix must be wrong.
! 
 995 ierr = -1
 return
! 
!      insufficient storage in arrays alu, jlu to store factors
! 
 996 ierr = -2
 return
! 
!      zero row encountered
! 
 997 ierr = -3 
 return
! ----------------end-of-iludp---------------------------!----------------
! -----------------------------------------------------------------------
 end subroutine iludp
 subroutine iluk(n,a,ja,ia,lfil,alu,jlu,ju,levs,iwk,w,jw,ierr)
 implicit none 
! BEGIN new declarations
! END new declarations
 integer, intent(In) :: n
 real(kind=8), dimension(1:*), intent(In) :: a
 real(kind=8), dimension(1:*), intent(InOut) :: alu
 real(kind=8), dimension(1:n), intent(Out) :: w
 integer, dimension(1:*), intent(In) :: ja
 integer, dimension(1:n + 1), intent(In) :: ia
 integer, dimension(1:*), intent(Out) :: jlu
 integer, dimension(1:n), intent(InOut) :: ju
 integer, dimension(1:*), intent(InOut) :: levs
 integer, dimension(1:3 * n), intent(Out) :: jw
 integer, intent(In) :: lfil
 integer, intent(In) :: iwk
 integer, intent(Out) :: ierr
! ----------------------------------------------------------------------* 
!      SPARSKIT ROUTINE ILUK -- ILU WITH LEVEL OF FILL-IN OF K (ILU(k)) *
! ----------------------------------------------------------------------*
! 
!  on entry:
! ========== 
!  n       = integer. The row dimension of the matrix A. The matrix 
! 
!  a,ja,ia = matrix stored in Compressed Sparse Row format.              
! 
!  lfil    = integer. The fill-in parameter. Each element whose
!            leve-of-fill exceeds lfil during the ILU process is dropped.
!            lfil must be .ge. 0 
! 
!  tol     = real*8. Sets the threshold for dropping small terms in the
!            factorization. See below for details on dropping strategy.
!   
!  iwk     = integer. The minimum length of arrays alu, jlu, and levs.
! 
!  On return:
! ===========
! 
!  alu,jlu = matrix stored in Modified Sparse Row (MSR) format containing
!            the L and U factors together. The diagonal (stored in
!            alu(1:n) ) is inverted. Each i-th row of the alu,jlu matrix
!            contains the i-th row of L (excluding the diagonal entry=1)
!            followed by the i-th row of U.
! 
!  ju      = integer array of length n containing the pointers to
!            the beginning of each row of U in the matrix alu,jlu.
! 
!  levs    = integer (work) array of size iwk -- which contains the 
!            levels of each element in alu, jlu.
! 
!  ierr    = integer. Error message with the following meaning.
!            ierr  = 0    --> successful return.
!            ierr .gt. 0  --> zero pivot encountered at step number ierr.
!            ierr  = -1   --> Error. input matrix may be wrong.
!                             (The elimination process has generated a
!                             row in L or U whose length is .gt.  n.)
!            ierr  = -2   --> The matrix L overflows the array al.
!            ierr  = -3   --> The matrix U overflows the array alu.
!            ierr  = -4   --> Illegal value for lfil.
!            ierr  = -5   --> zero row encountered in A or U.
! 
!  work arrays:
! =============
!  jw      = integer work array of length 3*n.
!  w       = real work array of length n 
! 
!  Notes/known bugs: This is not implemented efficiently storage-wise.
!        For example: Only the part of the array levs(*) associated with
!        the U-matrix is needed in the routine.. So some storage can 
!        be saved if needed. The levels of fills in the LU matrix are
!        output for information only -- they are not needed by LU-solve. 
!         
! ----------------------------------------------------------------------
!  w, ju (1:n) store the working array [1:ii-1 = L-part, ii:n = u] 
!  jw(n+1:2n)  stores the nonzero indicator. 
!  
!  Notes:
!  ------
!  All the diagonal elements of the input matrix must be  nonzero.
! 
! ----------------------------------------------------------------------* 
!      locals
 integer :: ju0
 integer :: k
 integer :: j1
 integer :: j2
 integer :: j
 integer :: ii
 integer :: i
 integer :: lenl
 integer :: lenu
 integer :: jj
 integer :: jrow
 integer :: jpos
 integer :: n2
 integer :: jlev
 integer :: min
 real(kind=8) :: t
 real(kind=8) :: s
 real(kind=8) :: fact
 if (lfil < 0) goto 998
! -----------------------------------------------------------------------
!      initialize ju0 (points to next element to be added to alu,jlu)
!      and pointer array.
! -----------------------------------------------------------------------
 n2 = n+n 
 ju0 = n+2
 jlu(1) = ju0
! 
!      initialize nonzero indicator array + levs array -- 
! 
 do j=1,2*n 
 jw(j) = 0
 end do
! -----------------------------------------------------------------------
!      beginning of main loop.
! -----------------------------------------------------------------------
 do ii = 1, n
 j1 = ia(ii)
 j2 = ia(ii+1) - 1
!      
!      unpack L-part and U-part of row of A in arrays w 
!      
 lenu = 1
 lenl = 0
 jw(ii) = ii
 w(ii) = 0.0
 jw(n+ii) = ii
! 
 do j = j1, j2
 k = ja(j)
 t = a(j)
 if (t == 0.0) goto 170 
 if (k < ii) then
 lenl = lenl+1
 jw(lenl) = k
 w(lenl) = t
 jw(n2+lenl) = 0 
 jw(n+k) = lenl
 else if (k == ii) then
 w(ii) = t
 jw(n2+ii) = 0 
 else
 lenu = lenu+1
 jpos = ii+lenu-1 
 jw(jpos) = k
 w(jpos) = t
 jw(n2+jpos) = 0 
 jw(n+k) = jpos
 endif
 170 continue
 end do
! 
 jj = 0
! 
!      eliminate previous rows
!      
 150 jj = jj+1
 if (jj > lenl) goto 160
! -----------------------------------------------------------------------
!      in order to do the elimination in the correct order we must select
!      the smallest column index among jw(k), k=jj+1, ..., lenl.
! -----------------------------------------------------------------------
 jrow = jw(jj)
 k = jj
!      
!      determine smallest column index
!      
 do j=jj+1,lenl
 if (jw(j) < jrow) then
 jrow = jw(j)
 k = j
 endif
 end do
! 
 if (k /= jj) then
!      exchange in jw
 j = jw(jj)
 jw(jj) = jw(k)
 jw(k) = j
!      exchange in jw(n+  (pointers/ nonzero indicator).
 jw(n+jrow) = jj
 jw(n+j) = k
!      exchange in jw(n2+  (levels) 
 j = jw(n2+jj) 
 jw(n2+jj) = jw(n2+k) 
 jw(n2+k) = j
!      exchange in w
 s = w(jj)
 w(jj) = w(k)
 w(k) = s
 endif
! 
!      zero out element in row by resetting jw(n+jrow) to zero.
!      
 jw(n+jrow) = 0
!      
!      get the multiplier for row to be eliminated (jrow) + its level
!      
 fact = w(jj)*alu(jrow)
 jlev = jw(n2+jj) 
 if (jlev > lfil) goto 150
! 
!      combine current row and row jrow
! 
 do k = ju(jrow), jlu(jrow+1)-1
 s = fact*alu(k)
 j = jlu(k)
 jpos = jw(n+j)
 if (j >= ii) then
!      
!      dealing with upper part.
!      
 if (jpos == 0) then
! 
!      this is a fill-in element
!      
 lenu = lenu+1
 if (lenu > n) goto 995
 i = ii+lenu-1
 jw(i) = j
 jw(n+j) = i
 w(i) = - s
 jw(n2+i) = jlev+levs(k)+1 
 else
! 
!      this is not a fill-in element 
! 
 w(jpos) = w(jpos) - s
 jw(n2+jpos) = min(jw(n2+jpos),jlev+levs(k)+1)
 endif
 else
!      
!      dealing with lower part.
!      
 if (jpos == 0) then
! 
!      this is a fill-in element
! 
 lenl = lenl+1
 if (lenl > n) goto 995
 jw(lenl) = j
 jw(n+j) = lenl
 w(lenl) = - s
 jw(n2+lenl) = jlev+levs(k)+1 
 else
! 
!      this is not a fill-in element 
! 
 w(jpos) = w(jpos) - s
 jw(n2+jpos) = min(jw(n2+jpos),jlev+levs(k)+1)
 endif
 endif
 end do
 w(jj) = fact
 jw(jj) = jrow
 goto 150 
 160 continue 
!      
!      reset double-pointer to zero (U-part) 
!      
 do k=1, lenu
 jw(n+jw(ii+k-1)) = 0
 end do
! 
!      update l-matrix
!          
 do k=1, lenl 
 if (ju0 > iwk) goto 996
 if (jw(n2+k) <= lfil) then
 alu(ju0) = w(k)
 jlu(ju0) = jw(k)
 ju0 = ju0+1
 endif
 end do
!      
!      save pointer to beginning of row ii of U
!      
 ju(ii) = ju0
! 
!      update u-matrix
! 
 do k=ii+1,ii+lenu-1 
 if (ju0 > iwk) goto 997
 if (jw(n2+k) <= lfil) then
 jlu(ju0) = jw(k)
 alu(ju0) = w(k)
 levs(ju0) = jw(n2+k) 
 ju0 = ju0+1
 endif
 end do
 if (w(ii) == 0.0) goto 999 
!      
 alu(ii) = 1.0d0/ w(ii) 
!      
!      update pointer to beginning of next row of U.
!      
 jlu(ii+1) = ju0
! -----------------------------------------------------------------------
!      end main loop
! -----------------------------------------------------------------------
 end do
 ierr = 0
 return
! 
!      incomprehensible error. Matrix must be wrong.
!      
 995 ierr = -1
 return
!      
!      insufficient storage in L.
!      
 996 ierr = -2
 return
!      
!      insufficient storage in U.
!      
 997 ierr = -3
 return
!      
!      illegal lfil entered.
!      
 998 ierr = -4
 return
!      
!      zero row encountered in A or U. 
!      
 999 ierr = -5
 return
! ----------------end-of-iluk--------------------------------------------
! -----------------------------------------------------------------------
 end subroutine iluk
 subroutine ilu0(n,a,ja,ia,alu,jlu,ju,iw,ierr)
 implicit none
! BEGIN new declarations
 integer, intent(In) :: n
 integer, intent(Out) :: ierr
 integer :: ju0
 integer :: i
 integer :: ii
 integer :: js
 integer :: j
 integer :: jcol
 integer :: jf
 integer :: jm
 integer :: jrow
 real :: tl
 integer :: jj
 integer :: jw
! END new declarations
 real(kind=8), dimension(1:*), intent(In) :: a
 real(kind=8), dimension(1:*), intent(Out) :: alu
 integer, dimension(1:*), intent(In) :: ja
 integer, dimension(1:*), intent(In) :: ia
 integer, dimension(1:*), intent(Out) :: ju
 integer, dimension(1:*), intent(Out) :: jlu
 integer, dimension(1:*), intent(Out) :: iw
! ------------------ right preconditioner ------------------------------*
!                     ***   ilu(0) preconditioner.   ***                *
! ----------------------------------------------------------------------*
!  Note that this has been coded in such a way that it can be used
!  with pgmres. Normally, since the data structure of the L+U matrix is
!  the same as that the A matrix, savings can be made. In fact with
!  some definitions (not correct for general sparse matrices) all we
!  need in addition to a, ja, ia is an additional diagonal.
!  ILU0 is not recommended for serious problems. It is only provided
!  here for comparison purposes.
! -----------------------------------------------------------------------
! 
!  on entry:
! ---------
!  n       = dimension of matrix
!  a, ja,
!  ia      = original matrix in compressed sparse row storage.
! 
!  on return:
! -----------
!  alu,jlu = matrix stored in Modified Sparse Row (MSR) format containing
!            the L and U factors together. The diagonal (stored in
!            alu(1:n) ) is inverted. Each i-th row of the alu,jlu matrix
!            contains the i-th row of L (excluding the diagonal entry=1)
!            followed by the i-th row of U.
! 
!  ju         = pointer to the diagonal elements in alu, jlu.
! 
!  ierr         = integer indicating error code on return
!             ierr = 0 --> normal return
!             ierr = k --> code encountered a zero pivot at step k.
!  work arrays:
! -------------
!  iw           = integer work array of length n.
! ------------
!  IMPORTANT
! -----------
!  it is assumed that the the elements in the input matrix are stored
!     in such a way that in each row the lower part comes first and
!     then the upper part. To get the correct ILU factorization, it is
!     also necessary to have the elements of L sorted by increasing
!     column number. It may therefore be necessary to sort the
!     elements of a, ja, ia prior to calling ilu0. This can be
!     achieved by transposing the matrix twice using csrcsc.
! 
! -----------------------------------------------------------------------
 ju0 = n+2
 jlu(1) = ju0
! 
!  initialize work vector to zero's
! 
 do i=1, n
 iw(i) = 0
 end do
! 
!  main loop
! 
 do ii = 1, n
 js = ju0
! 
!  generating row number ii of L and U.
! 
 do j=ia(ii),ia(ii+1)-1
! 
!      copy row ii of a, ja, ia into row ii of alu, jlu (L/U) matrix.
! 
 jcol = ja(j)
 if (jcol == ii) then
 alu(ii) = a(j)
 iw(jcol) = ii
 ju(ii) = ju0
 else
 alu(ju0) = a(j)
 jlu(ju0) = ja(j)
 iw(jcol) = ju0
 ju0 = ju0+1
 endif
 end do
 jlu(ii+1) = ju0
 jf = ju0-1
 jm = ju(ii)-1
! 
!      exit if diagonal element is reached.
! 
 do j=js, jm
 jrow = jlu(j)
 tl = alu(j)*alu(jrow)
 alu(j) = tl
! 
!      perform  linear combination
! 
 do jj = ju(jrow), jlu(jrow+1)-1
 jw = iw(jlu(jj))
 if (jw /= 0) alu(jw) = alu(jw) - tl*alu(jj)
 end do
 end do
! 
!      invert  and store diagonal element.
! 
 if (alu(ii) == 0.0d0) goto 600
 alu(ii) = 1.0d0/alu(ii)
! 
!      reset pointer iw to zero
! 
 iw(ii) = 0
 do i = js, jf
 201 iw(jlu(i)) = 0
 end do
 end do
 ierr = 0
 return
! 
!      zero pivot :
! 
 600 ierr = ii
! 
 return
! ------- end-of-ilu0 ---------------------------------------------------
! -----------------------------------------------------------------------
 end subroutine ilu0
 subroutine milu0(n,a,ja,ia,alu,jlu,ju,iw,ierr)
 implicit none
! BEGIN new declarations
 integer, intent(In) :: n
 integer, intent(Out) :: ierr
 integer :: ju0
 integer :: i
 integer :: ii
 integer :: js
 integer :: j
 integer :: jcol
 integer :: jf
 integer :: jm
 real :: s
 integer :: jrow
 real :: tl
 integer :: jj
 integer :: jw
! END new declarations
 real(kind=8), dimension(1:*), intent(In) :: a
 real(kind=8), dimension(1:*), intent(Out) :: alu
 integer, dimension(1:*), intent(In) :: ja
 integer, dimension(1:*), intent(In) :: ia
 integer, dimension(1:*), intent(Out) :: ju
 integer, dimension(1:*), intent(Out) :: jlu
 integer, dimension(1:*), intent(Out) :: iw
! ----------------------------------------------------------------------*
!                 *** simple milu(0) preconditioner. ***                *
! ----------------------------------------------------------------------*
!  Note that this has been coded in such a way that it can be used
!  with pgmres. Normally, since the data structure of a, ja, ia is
!  the same as that of a, ja, ia, savings can be made. In fact with
!  some definitions (not correct for general sparse matrices) all we
!  need in addition to a, ja, ia is an additional diagonal.
!  Ilu0 is not recommended for serious problems. It is only provided
!  here for comparison purposes.
! -----------------------------------------------------------------------
! 
!  on entry:
! ----------
!  n       = dimension of matrix
!  a, ja,
!  ia      = original matrix in compressed sparse row storage.
! 
!  on return:
! ----------
!  alu,jlu = matrix stored in Modified Sparse Row (MSR) format containing
!            the L and U factors together. The diagonal (stored in
!            alu(1:n) ) is inverted. Each i-th row of the alu,jlu matrix
!            contains the i-th row of L (excluding the diagonal entry=1)
!            followed by the i-th row of U.
! 
!  ju         = pointer to the diagonal elements in alu, jlu.
! 
!  ierr         = integer indicating error code on return
!             ierr = 0 --> normal return
!             ierr = k --> code encountered a zero pivot at step k.
!  work arrays:
! -------------
!  iw           = integer work array of length n.
! ------------
!  Note (IMPORTANT):
! -----------
!  it is assumed that the the elements in the input matrix are ordered
!     in such a way that in each row the lower part comes first and
!     then the upper part. To get the correct ILU factorization, it is
!     also necessary to have the elements of L ordered by increasing
!     column number. It may therefore be necessary to sort the
!     elements of a, ja, ia prior to calling milu0. This can be
!     achieved by transposing the matrix twice using csrcsc.
! -----------------------------------------------------------
 ju0 = n+2
 jlu(1) = ju0
!  initialize work vector to zero's
 do i=1, n
 31 iw(i) = 0
 end do
! 
! -------------- MAIN LOOP ----------------------------------
! 
 do ii = 1, n
 js = ju0
! 
!  generating row number ii or L and U.
! 
 do j=ia(ii),ia(ii+1)-1
! 
!      copy row ii of a, ja, ia into row ii of alu, jlu (L/U) matrix.
!      
 jcol = ja(j)
 if (jcol == ii) then
 alu(ii) = a(j)
 iw(jcol) = ii
 ju(ii) = ju0
 else
 alu(ju0) = a(j)
 jlu(ju0) = ja(j)
 iw(jcol) = ju0
 ju0 = ju0+1
 endif
 end do
 jlu(ii+1) = ju0
 jf = ju0-1
 jm = ju(ii)-1
!      s accumulates fill-in values
 s = 0.0d0
 do j=js, jm
 jrow = jlu(j)
 tl = alu(j)*alu(jrow)
 alu(j) = tl
! -----------------------perform linear combination --------
 do jj = ju(jrow), jlu(jrow+1)-1
 jw = iw(jlu(jj))
 if (jw /= 0) then
 alu(jw) = alu(jw) - tl*alu(jj)
 else
 s = s + tl*alu(jj)
 endif
 end do
 end do
! ----------------------- invert and store diagonal element.
 alu(ii) = alu(ii)-s
 if (alu(ii) == 0.0d0) goto 600
 alu(ii) = 1.0d0/alu(ii)
! ----------------------- reset pointer iw to zero
 iw(ii) = 0
 do i = js, jf
 201 iw(jlu(i)) = 0
 end do
 end do
 ierr = 0
 return
!      zero pivot :
 600 ierr = ii
 return
! ------- end-of-milu0 --------------------------------------------------
! -----------------------------------------------------------------------
 end subroutine milu0
 subroutine pgmres(n,im,rhs,sol,vv,eps,maxits,iout,aa,ja,ia,alu,jlu,ju,ierr)
! -----------------------------------------------------------------------
 implicit none
! BEGIN new declarations
 real :: eps1
 integer :: its
 integer :: j
 real :: ro
 real :: dnrm2
 integer :: i
 integer :: i1
 real :: ddot
 integer :: k
 integer :: k1
 real :: gam
 real :: epsmac
 integer :: ii
 integer :: jj
! END new declarations
 integer, intent(In) :: n
 integer, intent(In) :: im
 integer, intent(In) :: maxits
 integer :: iout
 integer, intent(Out) :: ierr
 integer, dimension(1:*), intent(In) :: ja
 integer, dimension(1:n + 1), intent(In) :: ia
 integer, dimension(1:*), intent(In) :: jlu
 integer, dimension(1:n), intent(In) :: ju
 real(kind=8), dimension(1:n,1:*), intent(InOut) :: vv
 real(kind=8), dimension(1:n), intent(InOut) :: rhs
 real(kind=8), dimension(1:n), intent(InOut) :: sol
 real(kind=8), dimension(1:*), intent(In) :: aa
 real(kind=8), dimension(1:*), intent(In) :: alu
 real(kind=8), intent(In) :: eps
! ----------------------------------------------------------------------*
!                                                                       *
!                  *** ILUT - Preconditioned GMRES ***                  *
!                                                                       *
! ----------------------------------------------------------------------*
!  This is a simple version of the ILUT preconditioned GMRES algorithm. *
!  The ILUT preconditioner uses a dual strategy for dropping elements   *
!  instead  of the usual level of-fill-in approach. See details in ILUT *
!  subroutine documentation. PGMRES uses the L and U matrices generated *
!  from the subroutine ILUT to precondition the GMRES algorithm.        *
!  The preconditioning is applied to the right. The stopping criterion  *
!  utilized is based simply on reducing the residual norm by epsilon.   *
!  This preconditioning is more reliable than ilu0 but requires more    *
!  storage. It seems to be much less prone to difficulties related to   *
!  strong nonsymmetries in the matrix. We recommend using a nonzero tol *
!  (tol=.005 or .001 usually give good results) in ILUT. Use a large    *
!  lfil whenever possible (e.g. lfil = 5 to 10). The higher lfil the    *
!  more reliable the code is. Efficiency may also be much improved.     *
!  Note that lfil=n and tol=0.0 in ILUT  will yield the same factors as *
!  Gaussian elimination without pivoting.                               *
!                                                                       *
!  ILU(0) and MILU(0) are also provided for comparison purposes         *
!  USAGE: first call ILUT or ILU0 or MILU0 to set up preconditioner and *
!  then call pgmres.                                                    *
! ----------------------------------------------------------------------*
!  Coded by Y. Saad - This version dated May, 7, 1990.                  *
! ----------------------------------------------------------------------*
!  parameters                                                           *
! -----------                                                           *
!  on entry:                                                            *
! ==========                                                            *
!                                                                       *
!  n     == integer. The dimension of the matrix.                       *
!  im    == size of krylov subspace:  should not exceed 50 in this      *
!           version (can be reset by changing parameter command for     *
!           kmax below)                                                 *
!  rhs   == real vector of length n containing the right hand side.     *
!           Destroyed on return.                                        *
!  sol   == real vector of length n containing an initial guess to the  *
!           solution on input. approximate solution on output           *
!  eps   == tolerance for stopping criterion. process is stopped        *
!           as soon as ( ||.|| is the euclidean norm):                  *
!           || current residual||/||initial residual|| <= eps           *
!  maxits== maximum number of iterations allowed                        *
!  iout  == output unit number number for printing intermediate results *
!           if (iout .le. 0) nothing is printed out.                    *
!                                                                       *
!  aa, ja,                                                              *
!  ia    == the input matrix in compressed sparse row format:           *
!           aa(1:nnz)  = nonzero elements of A stored row-wise in order *
!           ja(1:nnz) = corresponding column indices.                   *
!           ia(1:n+1) = pointer to beginning of each row in aa and ja.  *
!           here nnz = number of nonzero elements in A = ia(n+1)-ia(1)  *
!                                                                       *
!  alu,jlu== A matrix stored in Modified Sparse Row format containing   *
!            the L and U factors, as computed by subroutine ilut.       *
!                                                                       *
!  ju     == integer array of length n containing the pointers to       *
!            the beginning of each row of U in alu, jlu as computed     *
!            by subroutine ILUT.                                        *
!                                                                       *
!  on return:                                                           *
! ==========                                                            *
!  sol   == contains an approximate solution (upon successful return).  *
!  ierr  == integer. Error message with the following meaning.          *
!           ierr = 0 --> successful return.                             *
!           ierr = 1 --> convergence not achieved in itmax iterations.  *
!           ierr =-1 --> the initial guess seems to be the exact        *
!                        solution (initial residual computed was zero)  *
!                                                                       *
! ----------------------------------------------------------------------*
!                                                                       *
!  work arrays:                                                         *
! =============                                                         *
!  vv    == work array of length  n x (im+1) (used to store the Arnoli  *
!           basis)                                                      *
! ----------------------------------------------------------------------*
!  subroutines called :                                                 *
!  amux   : SPARSKIT routine to do the matrix by vector multiplication  *
!           delivers y=Ax, given x  -- see SPARSKIT/BLASSM/amux         *
!  lusol : combined forward and backward solves (Preconditioning ope.) *
!  BLAS1  routines.                                                     *
! ----------------------------------------------------------------------*
 integer, parameter :: kmax=50
 real(kind=8), dimension(1:kmax + 1,1:kmax) :: hh
 real(kind=8), dimension(1:kmax) :: c
 real(kind=8), dimension(1:kmax) :: s
 real(kind=8), dimension(1:kmax + 1) :: rs
 real(kind=8) :: t
! -------------------------------------------------------------
!  arnoldi size should not exceed kmax=50 in this version..
!  to reset modify paramter kmax accordingly.
! -------------------------------------------------------------
 data epsmac / 1.d-16 / 
 eps1 = 0.0
 its = 0
! -------------------------------------------------------------
!  outer loop starts here..
! -------------- compute initial residual vector --------------
 call amux (n, sol, vv, aa, ja, ia)
 do j=1,n
 vv(j,1) = rhs(j) - vv(j,1)
 end do
! -------------------------------------------------------------
 20 ro = dnrm2(n, vv, 1)
       if (iout  >  0 .and. its  ==  0)      write(iout, 199) its, ro
 if (ro == 0.0d0) goto 999
 t = 1.0d0/ ro
 do j=1, n
 vv(j,1) = vv(j,1)*t
 end do
 if (its == 0) eps1=eps*ro
!      ** initialize 1-st term  of rhs of hessenberg system..
 rs(1) = ro
 i = 0
 4 i=i+1
 its = its + 1
 i1 = i + 1
 call lusol(n,vv(1,i),rhs,alu,jlu,ju) 
 call amux (n, rhs, vv(1,i1), aa, ja, ia)
! -----------------------------------------
!      modified gram - schmidt...
! -----------------------------------------
 do j=1, i
 t = ddot(n, vv(1,j),1,vv(1,i1),1)
 hh(j,i) = t
 call daxpy(n, -t, vv(1,j), 1, vv(1,i1), 1)
 end do
 t = dnrm2(n, vv(1,i1), 1)
 hh(i1,i) = t
 if ( t == 0.0d0) goto 58
 t = 1.0d0/t
 do k=1,n
 vv(k,i1) = vv(k,i1)*t
 end do
! 
!      done with modified gram schimd and arnoldi step..
!      now  update factorization of hh
! 
 58 if (i == 1) goto 121
! --------perfrom previous transformations  on i-th column of h
 do k=2,i
 k1 = k-1
 t = hh(k1,i)
 hh(k1,i) = c(k1)*t + s(k1)*hh(k,i)
 hh(k,i) = -s(k1)*t + c(k1)*hh(k,i)
 end do
 121 gam = sqrt(hh(i,i)**2 + hh(i1,i)**2)
! 
!      if gamma is zero then any small value will do...
!      will affect only residual estimate
! 
 if (gam == 0.0d0) gam = epsmac
! 
!      get  next plane rotation
! 
 c(i) = hh(i,i)/gam
 s(i) = hh(i1,i)/gam
 rs(i1) = -s(i)*rs(i)
 rs(i) = c(i)*rs(i)
! 
!      detrermine residual norm and test for convergence-
! 
 hh(i,i) = c(i)*hh(i,i) + s(i)*hh(i1,i)
 ro = abs(rs(i1))
!  131   format(1h ,2e14.4)
       if (iout  >  0)      write(iout, 199) its, ro
 if (i < im .and. (ro > eps1)) goto 4
! 
!      now compute solution. first solve upper triangular system.
! 
 rs(i) = rs(i)/hh(i,i)
 do ii=2,i
 k=i-ii+1
 k1 = k+1
 t=rs(k)
 do j=k1,i
 t = t-hh(k,j)*rs(j)
 end do
 rs(k) = t/hh(k,k)
 end do
! 
!      form linear combination of v(*,i)'s to get solution
! 
 t = rs(1)
 do k=1, n
 rhs(k) = vv(k,1)*t
 end do
 do j=2, i
 t = rs(j)
 do k=1, n
 rhs(k) = rhs(k)+t*vv(k,j)
 end do
 end do
! 
!      call preconditioner.
! 
 call lusol(n,rhs,rhs,alu,jlu,ju) 
 do k=1, n
 sol(k) = sol(k) + rhs(k)
 end do
! 
!      restart outer loop  when necessary
! 
 if (ro <= eps1) goto 990
 if (its >= maxits) goto 991
! 
!      else compute residual vector and continue..
! 
 do j=1,i
 jj = i1-j+1
 rs(jj-1) = -s(jj-1)*rs(jj)
 rs(jj) = c(jj-1)*rs(jj)
 end do
 do j=1,i1
 t = rs(j)
 if (j == 1) t = t-1.0d0
 call daxpy (n, t, vv(1,j), 1, vv, 1)
 end do
 199   format('   its =', i4, ' res. norm =', d20.6)
!      restart outer loop.
 goto 20
 990 ierr = 0
 return
 991 ierr = 1
 return
 999 continue
 ierr = -1
 return
! -----------------end of pgmres ---------------------------------------
! -----------------------------------------------------------------------
 end subroutine pgmres
 subroutine lusol(n,y,x,alu,jlu,ju)
! BEGIN new declarations
! END new declarations
 implicit none
 real(kind=8), dimension(1:n), intent(Out) :: x
 real(kind=8), dimension(1:n), intent(In) :: y
 real(kind=8), dimension(1:*), intent(In) :: alu
 integer, intent(In) :: n
 integer, dimension(1:*), intent(In) :: jlu
 integer, dimension(1:*), intent(In) :: ju
! -----------------------------------------------------------------------
! 
!  This routine solves the system (LU) x = y, 
!  given an LU decomposition of a matrix stored in (alu, jlu, ju) 
!  modified sparse row format 
! 
! -----------------------------------------------------------------------
!  on entry:
!  n   = dimension of system 
!  y   = the right-hand-side vector
!  alu, jlu, ju 
!      = the LU matrix as provided from the ILU routines. 
! 
!  on return
!  x   = solution of LU x = y.     
! -----------------------------------------------------------------------
!  
!  Note: routine is in place: call lusol (n, x, x, alu, jlu, ju) 
!        will solve the system with rhs x and overwrite the result on x . 
! 
! -----------------------------------------------------------------------
!  local variables
! 
 integer :: i
 integer :: k
! 
!  forward solve
! 
 do i = 1, n
 x(i) = y(i)
 do k=jlu(i),ju(i)-1
 x(i) = x(i) - alu(k)* x(jlu(k))
 end do
 end do
! 
!      backward solve.
! 
 do i = n, 1, -1
 do k=ju(i),jlu(i+1)-1
 x(i) = x(i) - alu(k)*x(jlu(k))
 end do
 x(i) = alu(i)*x(i)
 end do
! 
 return
! ----------------end of lusol ------------------------------------------
! -----------------------------------------------------------------------
 end subroutine lusol
 subroutine lutsol(n,y,x,alu,jlu,ju)
! BEGIN new declarations
! END new declarations
 implicit none
 real(kind=8), dimension(1:n), intent(Out) :: x
 real(kind=8), dimension(1:n), intent(In) :: y
 real(kind=8), dimension(1:*), intent(In) :: alu
 integer, intent(In) :: n
 integer, dimension(1:*), intent(In) :: jlu
 integer, dimension(1:*), intent(In) :: ju
! -----------------------------------------------------------------------
! 
!  This routine solves the system  Transp(LU) x = y,
!  given an LU decomposition of a matrix stored in (alu, jlu, ju) 
!  modified sparse row format. Transp(M) is the transpose of M. 
! ----------------------------------------------------------------------- 
!  on entry:
!  n   = dimension of system 
!  y   = the right-hand-side vector
!  alu, jlu, ju 
!      = the LU matrix as provided from the ILU routines. 
! 
!  on return
!  x   = solution of transp(LU) x = y.   
! -----------------------------------------------------------------------
! 
!  Note: routine is in place: call lutsol (n, x, x, alu, jlu, ju) 
!        will solve the system with rhs x and overwrite the result on x . 
!  
! -----------------------------------------------------------------------
!  local variables
! 
 integer :: i
 integer :: k
! 
 do i = 1, n
 x(i) = y(i)
 end do
! 
!  forward solve (with U^T)
! 
 do i = 1, n
 x(i) = x(i) * alu(i)
 do k=ju(i),jlu(i+1)-1
 x(jlu(k)) = x(jlu(k)) - alu(k)* x(i)
 end do
 end do
!      
!      backward solve (with L^T)
!      
 do i = n, 1, -1 
 do k=jlu(i),ju(i)-1
 x(jlu(k)) = x(jlu(k)) - alu(k)*x(i)
 end do
 end do
! 
 return
! ----------------end of lutsol -----------------------------------------
! -----------------------------------------------------------------------
 end subroutine lutsol
 subroutine qsplit(a,ind,n,ncut)
! BEGIN new declarations
 implicit none
 integer :: mid
 integer :: j
! END new declarations
 real(kind=8), dimension(1:n), intent(InOut) :: a
 integer, dimension(1:n), intent(InOut) :: ind
 integer, intent(In) :: n
 integer, intent(In) :: ncut
! -----------------------------------------------------------------------
!      does a quick-sort split of a real array.
!      on input a(1:n). is a real array
!      on output a(1:n) is permuted such that its elements satisfy:
! 
!      abs(a(i)) .ge. abs(a(ncut)) for i .lt. ncut and
!      abs(a(i)) .le. abs(a(ncut)) for i .gt. ncut
! 
!      ind(1:n) is an integer array which permuted in the same way as a(*).
! -----------------------------------------------------------------------
 real(kind=8) :: tmp
 real(kind=8) :: abskey
 integer :: itmp
 integer :: first
 integer :: last
! -----
 first = 1
 last = n
 if (ncut < first .or. ncut > last) return
! 
!      outer loop -- while mid .ne. ncut do
! 
 1 mid = first
 abskey = abs(a(mid))
 do j=first+1, last
 if (abs(a(j)) > abskey) then
 mid = mid+1
!      interchange
 tmp = a(mid)
 itmp = ind(mid)
 a(mid) = a(j)
 ind(mid) = ind(j)
 a(j) = tmp
 ind(j) = itmp
 endif
 end do
! 
!      interchange
! 
 tmp = a(mid)
 a(mid) = a(first)
 a(first) = tmp
! 
 itmp = ind(mid)
 ind(mid) = ind(first)
 ind(first) = itmp
! 
!      test for while loop
! 
 if (mid == ncut) return
 if (mid > ncut) then
 last = mid-1
 else
 first = mid+1
 endif
 goto 1
! ----------------end-of-qsplit------------------------------------------
! -----------------------------------------------------------------------
 end subroutine qsplit
