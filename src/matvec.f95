subroutine amux(n,x,y,a,ja,ia)
! BEGIN declarations

    real(kind=8), dimension(*), intent(In) :: x, a
    real(kind=8), dimension(*), intent(Out) :: y
    integer, intent(In) :: n
    integer, dimension(*) :: ja
    integer, dimension(*), intent(In) :: ia
    real(kind=8) :: t
    integer :: i
    integer :: k
! END declarations
! -----------------------------------------------------------------------
!          A times a vector
! ----------------------------------------------------------------------- 
!  multiplies a matrix by a vector using the dot product form
!  Matrix A is stored in compressed sparse row storage.
! 
!  on entry:
! ----------
!  n     = row dimension of A
!  x     = real array of length equal to the column dimension of
!          the A matrix.
!  a, ja,
!     ia = input matrix in compressed sparse row format.
! 
!  on return:
! -----------
!  y     = real array of length n, containing the product y=Ax
! 
! -----------------------------------------------------------------------
    do i = 1,n
! 
!      compute the inner product of row i with vector x
!  
        t = 0.0d0
        do k=ia(i), ia(i+1)-1 
            t = t + a(k)*x(ja(k))
        end do
! 
!      store result in y(i) 
! 
        y(i) = t
    end do
! 
    return
! ---------end-of-amux---------------------------------------------------
! -----------------------------------------------------------------------
end subroutine amux

subroutine amuxms(n,x,y,a,ja)
! BEGIN declarations

    real(kind=8), dimension(*), intent(In) :: x, a
    real(kind=8), dimension(*), intent(Out) :: y
    integer, intent(In) :: n
    integer, dimension(*), intent(In) :: ja
    integer :: i,k
! END declarations
! -----------------------------------------------------------------------
!          A times a vector in MSR format
! -----------------------------------------------------------------------
!  multiplies a matrix by a vector using the dot product form
!  Matrix A is stored in Modified Sparse Row storage.
! 
!  on entry:
! ----------
!  n     = row dimension of A
!  x     = real array of length equal to the column dimension of
!          the A matrix.
!  a, ja,= input matrix in modified compressed sparse row format.
! 
!  on return:
! -----------
!  y     = real array of length n, containing the product y=Ax
! 
! -----------------------------------------------------------------------
!  local variables
! 

! -----------------------------------------------------------------------
    do i=1, n
        y(i) = a(i)*x(i)
    end do
 
    do i = 1,n
! 
!      compute the inner product of row i with vector x
! 
        do k=ja(i), ja(i+1)-1
            y(i) = y(i) + a(k) *x(ja(k))
        end do
    end do
! 
    return
! ---------end-of-amuxm--------------------------------------------------
! -----------------------------------------------------------------------
end subroutine amuxms

subroutine atmux(n,x,y,a,ja,ia)
! BEGIN new declarations
    real(kind=8), dimension(*), intent(In) :: x, a
    real(kind=8), dimension(*), intent(Out) :: y
    integer, intent(In) :: n
    integer, dimension(*), intent(In) :: ia
    integer, dimension(*) :: ja
    integer :: i
    integer :: k
! END new declarations
! -----------------------------------------------------------------------
!          transp( A ) times a vector
! ----------------------------------------------------------------------- 
!  multiplies the transpose of a matrix by a vector when the original
!  matrix is stored in compressed sparse row storage. Can also be
!  viewed as the product of a matrix by a vector when the original
!  matrix is stored in the compressed sparse column format.
! -----------------------------------------------------------------------
! 
!  on entry:
! ----------
!  n     = row dimension of A
!  x     = real array of length equal to the column dimension of
!          the A matrix.
!  a, ja,
!     ia = input matrix in compressed sparse row format.
! 
!  on return:
! -----------
!  y     = real array of length n, containing the product y=transp(A)*x
! 
! -----------------------------------------------------------------------

! -----------------------------------------------------------------------
! 
!      zero out output vector
!  
    do i=1,n
        y(i) = 0.0
    end do
! 
!  loop over the rows
! 
    do i = 1,n
        do k=ia(i), ia(i+1)-1 
            y(ja(k)) = y(ja(k)) + x(i)*a(k)
        end do
    end do
! 
    return
! -------------end-of-atmux---------------------------------------------- 
! -----------------------------------------------------------------------
end subroutine atmux

subroutine atmuxr(m,n,x,y,a,ja,ia)
! BEGIN declarations

    real(kind=8), dimension(*), intent(In) :: x, a
    real(kind=8), dimension(*), intent(Out) :: y
    integer, intent(In) :: m, n
    integer, dimension(*), intent(In) :: ia
    integer, dimension(*) :: ja
!      local variables 
! 
    integer :: i,k
! END declarations
! -----------------------------------------------------------------------
!          transp( A ) times a vector, A can be rectangular
! ----------------------------------------------------------------------- 
!  See also atmux.  The essential difference is how the solution vector
!  is initially zeroed.  If using this to multiply rectangular CSC 
!  matrices by a vector, m number of rows, n is number of columns.
! -----------------------------------------------------------------------
! 
!  on entry:
! ----------
!  m     = column dimension of A
!  n     = row dimension of A
!  x     = real array of length equal to the column dimension of
!          the A matrix.
!  a, ja,
!     ia = input matrix in compressed sparse row format.
! 
!  on return:
! -----------
!  y     = real array of length n, containing the product y=transp(A)*x
! 
! -----------------------------------------------------------------------

! -----------------------------------------------------------------------
! 
!      zero out output vector
!  
    do i=1,m
        y(i) = 0.0
    end do
! 
!  loop over the rows
! 
    do i = 1,n
        do k=ia(i), ia(i+1)-1 
            y(ja(k)) = y(ja(k)) + x(i)*a(k)
        end do
    end do
! 
    return
! -------------end-of-atmuxr--------------------------------------------- 
! -----------------------------------------------------------------------
end subroutine atmuxr

subroutine amuxe(n,x,y,na,ncol,a,ja)
! BEGIN declarations

    real(kind=8), dimension(n), intent(In) :: x, a
    real(kind=8), dimension(n), intent(Out) :: y
    real(kind=8), dimension(na,*), intent(In) :: a
    integer, intent(In) :: n, ncol
    integer :: na
    integer, dimension(na,*) :: ja
!  local variables
! 
    integer :: i,j
! -----------------------------------------------------------------------
! END declarations
! -----------------------------------------------------------------------
!         A times a vector in Ellpack Itpack format (ELL)               
! ----------------------------------------------------------------------- 
!  multiplies a matrix by a vector when the original matrix is stored 
!  in the ellpack-itpack sparse format.
! -----------------------------------------------------------------------
! 
!  on entry:
! ----------
!  n     = row dimension of A
!  x     = real array of length equal to the column dimension of
!          the A matrix.
!  na    = integer. The first dimension of arrays a and ja
!          as declared by the calling program.
!  ncol  = integer. The number of active columns in array a.
!          (i.e., the number of generalized diagonals in matrix.)
!  a, ja = the real and integer arrays of the itpack format
!          (a(i,k),k=1,ncol contains the elements of row i in matrix
!           ja(i,k),k=1,ncol contains their column numbers) 
! 
!  on return:
! -----------
!  y     = real array of length n, containing the product y=y=A*x
! 
! -----------------------------------------------------------------------
    y = 0.0
    do j=1,ncol
        do i = 1,n
            y(i) = y(i)+a(i,j)*x(ja(i,j))
        end do
    end do
! 
    return
! --------end-of-amuxe--------------------------------------------------- 
! -----------------------------------------------------------------------
end subroutine amuxe

subroutine amuxd(n,x,y,diag,ndiag,idiag,ioff)
! BEGIN new declarations

    integer, intent(In) :: n, idiag
    integer :: ndiag, j, k, io, i1,i2
    integer, dimension(idiag), intent(In) :: ioff
    real(kind=8), dimension(n), intent(In) :: x
    real(kind=8), dimension(n), intent(Out) :: y
    real(kind=8), dimension(ndiag,idiag), intent(In) :: diag
! -----------------------------------------------------------------------
! END new declarations
! -----------------------------------------------------------------------
!         A times a vector in Diagonal storage format (DIA) 
! ----------------------------------------------------------------------- 
!  multiplies a matrix by a vector when the original matrix is stored 
!  in the diagonal storage format.
! -----------------------------------------------------------------------
! 
!  on entry:
! ----------
!  n     = row dimension of A
!  x     = real array of length equal to the column dimension of
!          the A matrix.
!  ndiag  = integer. The first dimension of array adiag as declared in
!          the calling program.
!  idiag  = integer. The number of diagonals in the matrix.
!  diag   = real array containing the diagonals stored of A.
!  idiag  = number of diagonals in matrix.
!  diag   = real array of size (ndiag x idiag) containing the diagonals
!           
!  ioff   = integer array of length idiag, containing the offsets of the
!    	   diagonals of the matrix:
!           diag(i,k) contains the element a(i,i+ioff(k)) of the matrix.
! 
!  on return:
! -----------
!  y     = real array of length n, containing the product y=A*x
! 
! -----------------------------------------------------------------------
    y = 0.0
    do j=1, idiag
        io = ioff(j)
        i1 = max0(1,1-io)
        i2 = min0(n,n-io)
        do k=i1, i2
            y(k) = y(k)+diag(k,j)*x(k+io)
        end do
    end do
!  
    return
! ----------end-of-amuxd-------------------------------------------------
! -----------------------------------------------------------------------
end subroutine amuxd

subroutine amuxj(n,x,y,jdiag,a,ja,ia)
! BEGIN declarations

    integer, intent(In) :: n, jdiag
    integer, dimension(*) :: ja
    integer, dimension(*), intent(In) :: ia
    real(kind=8), dimension(n), intent(In) :: x
    real(kind=8), dimension(n), intent(Out) :: y
    real(kind=8), dimension(*), intent(In) :: a

!  local variables 
! 
    integer :: i, ii, k1, len, j
! END declarations
! -----------------------------------------------------------------------
!         A times a vector in Jagged-Diagonal storage format (JAD) 
! ----------------------------------------------------------------------- 
!  multiplies a matrix by a vector when the original matrix is stored 
!  in the jagged diagonal storage format.
! -----------------------------------------------------------------------
! 
!  on entry:
! ----------
!  n      = row dimension of A
!  x      = real array of length equal to the column dimension of
!          the A matrix.
!  jdiag  = integer. The number of jadded-diagonals in the data-structure.
!  a      = real array containing the jadded diagonals of A stored
!           in succession (in decreasing lengths) 
!  j      = integer array containing the colum indices of the 
!           corresponding elements in a.
!  ia     = integer array containing the lengths of the  jagged diagonals
! 
!  on return:
! -----------
!  y      = real array of length n, containing the product y=A*x
! 
!  Note:
! ------- 
!  Permutation related to the JAD format is not performed.
!  this can be done by:
!      call permvec (n,y,y,iperm) 
!  after the call to amuxj, where iperm is the permutation produced
!  by csrjad.
! -----------------------------------------------------------------------

    y = 0.0
 
    do ii=1, jdiag
        k1 = ia(ii)-1
        len = ia(ii+1)-k1-1
        do j=1,len
            y(j)= y(j)+a(k1+j)*x(ja(k1+j)) 
        end do
    end do
! 
    return
! ----------end-of-amuxj------------------------------------------------- 
! -----------------------------------------------------------------------
end subroutine amuxj

subroutine vbrmv(nr,nc,ia,ja,a,kvstr,kvstc,x,b)
! -----------------------------------------------------------------------
! BEGIN declarations

    integer, intent(In) :: nr
    integer :: nc, n, i, j, ii, jj, k, istart, istart
    real(kind = 8) :: xjj
    integer, dimension(nr + 1), intent(In) :: ia, kvstr
    integer, dimension(*), intent(In) :: ja, kvstc
    real(kind=8), dimension(*), intent(In) :: a, x
    real(kind=8), dimension(*), intent(Out) :: b
! END new declarations
! -----------------------------------------------------------------------
!      Sparse matrix-full vector product, in VBR format.
! -----------------------------------------------------------------------
!      On entry:
! --------------
!      nr, nc  = number of block rows and columns in matrix A
!      ia,ja,(),a,kvstr,kvstc = matrix A in variable block row format
!      x       = multiplier vector in full format
! 
!      On return:
! ---------------
!      b = product of matrix A times vector x in full format
! 
!      Algorithm:
! ---------------
!      Perform multiplication by traversing a in order.
! 
! -----------------------------------------------------------------------

    n = kvstc(nc+1)-1
    do i = 1, n
        b(i) = 0.d0
    enddo
! ---------------------------------
    k = 1
    do i = 1, nr
        istart = kvstr(i)
        istop = kvstr(i+1)-1
        do j = ia(i), ia(i+1)-1
            do jj = kvstc(ja(j)), kvstc(ja(j)+1)-1
                xjj = x(jj)
                do ii = istart, istop
                    b(ii) = b(ii) + xjj*a(k)
                    k = k + 1
                enddo
            enddo
        enddo
    enddo
! ---------------------------------
    return
end subroutine vbrmv

subroutine lsol(n,x,y,al,jal,ial)

    integer, intent(In) :: n
    integer, dimension(*) :: jal
    integer, dimension(n + 1), intent(In) :: ial
    real(kind=8), dimension(n), intent(Out) :: x
    real(kind=8), dimension(n), intent(In) :: y
    real(kind=8), dimension(*), intent(In) :: al
!  local variables 
! 
    integer :: k, j
    real(kind=8) :: t
! -----------------------------------------------------------------------

! -----------------------------------------------------------------------
!    solves    L x = y ; L = lower unit triang. /  CSR format
! ----------------------------------------------------------------------- 
!  solves a unit lower triangular system by standard (sequential )
!  forward elimination - matrix stored in CSR format. 
! -----------------------------------------------------------------------
! 
!  On entry:
! ---------- 
!  n      = integer. dimension of problem.
!  y      = real array containg the right side.
! 
!  al,
!  jal,
!  ial,    = Lower triangular matrix stored in compressed sparse row
!           format. 
! 
!  On return:
! ----------- 
!      x  = The solution of  L x  = y.
! --------------------------------------------------------------------
    x(1) = y(1) 
    do k = 2, n
        t = y(k) 
        do j = ial(k), ial(k+1)-1
            t = t-al(j)*x(jal(j))
        end do
        x(k) = t 
    end do
! 
    return
! ----------end-of-lsol-------------------------------------------------- 
! -----------------------------------------------------------------------
end subroutine lsol

subroutine ldsol(n,x,y,al,jal)

    integer, intent(In) :: n
    integer, dimension(*), intent(In) :: jal
    real(kind=8), dimension(n), intent(Out) :: x
    real(kind=8), dimension(n), intent(In) :: y, al
    real(kind=8), dimension(1:*), intent(In) :: al
!  local variables 
! 
    integer :: k, j
    real(kind=8) :: t
! -----------------------------------------------------------------------
! ----------------------------------------------------------------------- 
!      Solves L x = y    L = triangular. MSR format 
! -----------------------------------------------------------------------
!  solves a (non-unit) lower triangular system by standard (sequential) 
!  forward elimination - matrix stored in MSR format 
!  with diagonal elements already inverted (otherwise do inversion,
!  al(1:n) = 1.0/al(1:n),  before calling ldsol).
! -----------------------------------------------------------------------
! 
!  On entry:
! ---------- 
!  n      = integer. dimension of problem.
!  y      = real array containg the right hand side.
! 
!  al,
!  jal,   = Lower triangular matrix stored in Modified Sparse Row 
!           format. 
! 
!  On return:
! ----------- 
!      x = The solution of  L x = y .
! --------------------------------------------------------------------

    x(1) = y(1)*al(1) 
    do k = 2, n
        t = y(k) 
        do j = jal(k), jal(k+1)-1
            t = t - al(j)*x(jal(j))
        end do
        x(k) = al(k)*t 
    end do
    return
! ----------end-of-ldsol-------------------------------------------------
! -----------------------------------------------------------------------
end subroutine ldsol

subroutine lsolc(n,x,y,al,jal,ial)
    integer, intent(In) :: n
    integer, dimension(*) :: jal
    integer, dimension(*), intent(In) :: ial
    real(kind=8), dimension(n), intent(Out) :: x
    real(kind=8), dimension(n), intent(In) :: y
    real(kind=8), dimension(*), intent(In) :: al
!  local variables 
! 
    integer :: k, j
    real(kind=8) :: t
! -----------------------------------------------------------------------
! -----------------------------------------------------------------------
!        SOLVES     L x = y ;    where L = unit lower trang. CSC format
! -----------------------------------------------------------------------
!  solves a unit lower triangular system by standard (sequential )
!  forward elimination - matrix stored in CSC format. 
! -----------------------------------------------------------------------
! 
!  On entry:
! ---------- 
!  n      = integer. dimension of problem.
!  y      = real*8 array containg the right side.
! 
!  al,
!  jal,
!  ial,    = Lower triangular matrix stored in compressed sparse column 
!           format. 
! 
!  On return:
! ----------- 
!      x  = The solution of  L x  = y.
! -----------------------------------------------------------------------

    x = y
    do k = 1, n-1
        t = x(k) 
        do j = ial(k), ial(k+1)-1
            x(jal(j)) = x(jal(j)) - t*al(j) 
        end do
    end do
! 
    return
! ----------end-of-lsolc------------------------------------------------- 
! -----------------------------------------------------------------------
end subroutine lsolc

subroutine ldsolc(n,x,y,al,jal)
    integer, intent(In) :: n
    integer, dimension(*), intent(In) :: jal
    real(kind=8), dimension(n), intent(Out) :: x
    real(kind=8), dimension(n), intent(In) :: y
    real(kind=8), dimension(*), intent(In) :: al
!  local variables
! 
    integer :: k, j
    real(kind=8) :: t
! -----------------------------------------------------------------------

! -----------------------------------------------------------------------
!     Solves     L x = y ;    L = nonunit Low. Triang. MSC format 
! ----------------------------------------------------------------------- 
!  solves a (non-unit) lower triangular system by standard (sequential) 
!  forward elimination - matrix stored in Modified Sparse Column format 
!  with diagonal elements already inverted (otherwise do inversion,
!  al(1:n) = 1.0/al(1:n),  before calling ldsol).
! -----------------------------------------------------------------------
! 
!  On entry:
! ---------- 
!  n      = integer. dimension of problem.
!  y      = real array containg the right hand side.
! 
!  al,
!  jal,
!  ial,    = Lower triangular matrix stored in Modified Sparse Column
!            format.
! 
!  On return:
! ----------- 
!      x = The solution of  L x = y .
! --------------------------------------------------------------------
    x = y 
    do k = 1, n 
        x(k) = x(k)*al(k) 
        t = x(k) 
        do j = jal(k), jal(k+1)-1
            x(jal(j)) = x(jal(j)) - t*al(j) 
        end do
    end do
! 
    return
! ----------end-of-lsolc------------------------------------------------ 
! -----------------------------------------------------------------------
end subroutine ldsolc

subroutine ldsoll(n,x,y,al,jal,nlev,lev,ilev)

    integer :: k, n, ii, jrow, i
    integer, intent(In) :: nlev
    integer, dimension(*), intent(In) :: jal
    integer, dimension(nlev + 1), intent(In) :: ilev
    integer, dimension(n), intent(In) :: lev
    real(kind=8), dimension(n), intent(InOut) :: x
    real(kind=8), dimension(n), intent(In) :: y
    real(kind=8), dimension(*), intent(In) :: al
    real(kind=8) :: t
! -----------------------------------------------------------------------
!     Solves L x = y    L = triangular. Uses LEVEL SCHEDULING/MSR format 
! -----------------------------------------------------------------------
! 
!  On entry:
! ---------- 
!  n      = integer. dimension of problem.
!  y      = real array containg the right hand side.
! 
!  al,
!  jal,   = Lower triangular matrix stored in Modified Sparse Row 
!           format. 
!  nlev   = number of levels in matrix
!  lev    = integer array of length n, containing the permutation
!           that defines the levels in the level scheduling ordering.
!  ilev   = pointer to beginning of levels in lev.
!           the numbers lev(i) to lev(i+1)-1 contain the row numbers
!           that belong to level number i, in the level shcheduling
!           ordering.
! 
!  On return:
! ----------- 
!      x = The solution of  L x = y .
! --------------------------------------------------------------------

 
!      
!      outer loop goes through the levels. (SEQUENTIAL loop)
!      
    do ii=1, nlev
!      
!      next loop executes within the same level. PARALLEL loop
!      
        do i=ilev(ii), ilev(ii+1)-1 
            jrow = lev(i)
! 
!  compute inner product of row jrow with x
!  
            t = y(jrow) 
            do k=jal(jrow), jal(jrow+1)-1 
                t = t - al(k)*x(jal(k))
            end do
            x(jrow) = t*al(jrow) 
        end do
    end do
    return
! -----------------------------------------------------------------------
end subroutine ldsoll

subroutine usol(n,x,y,au,jau,iau)
    integer, intent(In) :: n
    integer, dimension(*) :: jau
    integer, dimension(n + 1), intent(In) :: iau
    real(kind=8), dimension(n), intent(InOut) :: x
    real(kind=8), dimension(n), intent(In) :: y
    real(kind=8), dimension(*), intent(In) :: au
!  local variables 
! 
    integer :: k,j
    real(kind=8) :: t
! -----------------------------------------------------------------------

! ----------------------------------------------------------------------- 
!              Solves   U x = y    U = unit upper triangular. 
! -----------------------------------------------------------------------
!  solves a unit upper triangular system by standard (sequential )
!  backward elimination - matrix stored in CSR format. 
! -----------------------------------------------------------------------
! 
!  On entry:
! ---------- 
!  n      = integer. dimension of problem.
!  y      = real array containg the right side.
! 
!  au,
!  jau,
!  iau,    = Lower triangular matrix stored in compressed sparse row
!           format. 
! 
!  On return:
! ----------- 
!      x = The solution of  U x = y . 
! -------------------------------------------------------------------- 


    x(n) = y(n) 
    do k = n-1,1,-1 
        t = y(k) 
        do j = iau(k), iau(k+1)-1
            t = t - au(j)*x(jau(j))
        end do
        x(k) = t 
    end do
! 
    return
! ----------end-of-usol-------------------------------------------------- 
! -----------------------------------------------------------------------
end subroutine usol

subroutine udsol(n,x,y,au,jau)
 
    integer, intent(In) :: n
    integer, dimension(*), intent(In) :: jau
    real(kind=8), dimension(n), intent(InOut) :: x
    real(kind=8), dimension(n), intent(In) :: y
    real(kind=8), dimension(*), intent(In) :: au
!  local variables 
! 
    integer :: k,j
    real(kind=8) :: t

! ----------------------------------------------------------------------- 
!              Solves   U x = y  ;   U = upper triangular in MSR format
! -----------------------------------------------------------------------
!  solves a non-unit upper triangular matrix by standard (sequential )
!  backward elimination - matrix stored in MSR format. 
!  with diagonal elements already inverted (otherwise do inversion,
!  au(1:n) = 1.0/au(1:n),  before calling).
! -----------------------------------------------------------------------
! 
!  On entry:
! ---------- 
!  n      = integer. dimension of problem.
!  y      = real array containg the right side.
! 
!  au,
!  jau,    = Lower triangular matrix stored in modified sparse row
!           format. 
! 
!  On return:
! ----------- 
!      x = The solution of  U x = y .
! --------------------------------------------------------------------
    x(n) = y(n)*au(n)
    do k = n-1,1,-1
        t = y(k) 
        do j = jau(k), jau(k+1)-1
            t = t - au(j)*x(jau(j))
        end do
        x(k) = au(k)*t 
    end do
! 
    return
! ----------end-of-udsol-------------------------------------------------
! -----------------------------------------------------------------------
end subroutine udsol

subroutine usolc(n,x,y,au,jau,iau)

    real(kind=8), dimension(*), intent(Out) :: x
    real(kind=8), dimension(*), intent(In) :: y, au
    integer, intent(In) :: n
    integer, dimension(*) :: jau
    integer, dimension(*), intent(In) :: iau
!  local variables 
!      
    integer :: k, j
    real(kind=8) :: t
! -----------------------------------------------------------------------
! -----------------------------------------------------------------------
!        SOUVES     U x = y ;    where U = unit upper trang. CSC format
! -----------------------------------------------------------------------
!  solves a unit upper triangular system by standard (sequential )
!  forward elimination - matrix stored in CSC format. 
! -----------------------------------------------------------------------
! 
!  On entry:
! ---------- 
!  n      = integer. dimension of problem.
!  y      = real*8 array containg the right side.
! 
!  au,
!  jau,
!  iau,    = Uower triangular matrix stored in compressed sparse column 
!           format. 
! 
!  On return:
! ----------- 
!      x  = The solution of  U x  = y.
! -----------------------------------------------------------------------

    do k=1,n
        x(k) = y(k) 
    end do
    do k = n,1,-1
        t = x(k) 
        do j = iau(k), iau(k+1)-1
            x(jau(j)) = x(jau(j)) - t*au(j) 
        end do
    end do
! 
    return
! ----------end-of-usolc------------------------------------------------- 
! -----------------------------------------------------------------------
end subroutine usolc
