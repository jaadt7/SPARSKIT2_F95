!*==dinfo1.f90 processed by SPAG 8.04RA 23:24 22 Feb 2024
!!SPAG Open source Personal, Educational or Academic User Clemson University  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
SUBROUTINE dinfo1(N,Iout,A,Ja,Ia,Valued,Title,Key,Type,Ao,Jao,Iao)
   USE ISO_FORTRAN_ENV                 
   IMPLICIT NONE
!
! Dummy argument declarations rewritten by SPAG
!
   INTEGER :: N
   INTEGER , INTENT(IN) :: Iout
   REAL(REAL64) , DIMENSION(*) :: A
   INTEGER , DIMENSION(*) :: Ja
   INTEGER , DIMENSION(N+1) :: Ia
   LOGICAL :: Valued
   CHARACTER(72) , INTENT(IN) :: Title
   CHARACTER(8) , INTENT(IN) :: Key
   CHARACTER(3) , INTENT(IN) :: Type
   REAL(REAL64) , DIMENSION(*) :: Ao
   INTEGER , INTENT(INOUT) , DIMENSION(*) :: Jao
   INTEGER , DIMENSION(N+1) :: Iao
!
! Local variable declarations rewritten by SPAG
!
   REAL(REAL64) :: amx
   REAL(REAL64) :: av , bndav , ddomc , ddomr , dist , fan , fas , fnorm , st , std
   REAL(REAL64) , DIMENSION(20) :: dcount
   INTEGER :: i , iband , imatch , ipos , j , jb1 , jb2 , job , k , k1 , k2 , ml , mu , n2 , nblk , ncol , ndiag , nlower , nnz ,  &
            & nrow , nsky , nupper , nzcol , nzdiag , nzmaxc , nzmaxr , nzminc , nzminr , nzrow
   INTEGER , DIMENSION(20) :: ioff
   INTEGER , SAVE :: ipar1
   LOGICAL :: sym
   CHARACTER(61) :: tmpst
   EXTERNAL ansym , avnz_col , bandpart , bandwidth , blkfnd , csrcsc , diag_domi , distaij , distdiag , frobnorm , nonz ,         &
          & nonz_lud , n_imp_diag , skyline
!
! End of declarations rewritten by SPAG
!
!----------------------------------------------------------------------c
!  SPARSKIT:  ELEMENTARY INFORMATION ROUTINE.                          c
!----------------------------------------------------------------------c
! info1 obtains a number of statistics on a sparse matrix and writes   c
! it into the output unit iout. The matrix is assumed                  c
! to be stored in the compressed sparse COLUMN format sparse a, ja, ia c
!----------------------------------------------------------------------c
! Modified Nov 1, 1989. 1) Assumes A is stored in column
! format. 2) Takes symmetry into account, i.e., handles Harwell-Boeing
!            matrices correctly.
!          ***  (Because of the recent modification the words row and
!            column may be mixed-up at occasions... to be checked...
!
! bug-fix July 25: 'upper' 'lower' mixed up in formats 108-107.
!
! On entry :
!-----------
! n	= integer. column dimension of matrix
! iout  = integer. unit number where the information it to be output.
! a	= real array containing the nonzero elements of the matrix
!	  the elements are stored by columns in order
!	  (i.e. column i comes before column i+1, but the elements
!         within each column can be disordered).
! ja	= integer array containing the row indices of elements in a
! ia	= integer array containing of length n+1 containing the
!         pointers to the beginning of the columns in arrays a and ja.
!	  It is assumed that ia(*) = 1 and ia(n+1) = nzz+1.
!
! valued= logical equal to .true. if values are provided and .false.
!         if only the pattern of the matrix is provided. (in that
!         case a(*) and ao(*) are dummy arrays.
!
! title = a 72-character title describing the matrix
!         NOTE: The first character in title is ignored (it is often
!         a one).
!
! key   = an 8-character key for the matrix
! type  = a 3-character string to describe the type of the matrix.
!         see harwell/Boeing documentation for more details on the
!         above three parameters.
!
! on return
!----------
! 1) elementary statistics on the matrix is written on output unit
!    iout. See below for detailed explanation of typical output.
! 2) the entries of a, ja, ia are sorted.
!
!----------
!
! ao	= real*8 array of length nnz used as work array.
! jao	= integer work array of length max(2*n+1,nnz)
! iao   = integer work array of length n+1
!
! Note  : title, key, type are the same paramaters as those
!         used for Harwell-Bowing matrices.
!
!-----------------------------------------------------------------------
! Output description:
!--------------------
! *** The following info needs to be updated.
!
! + A header containing the Title, key, type of the matrix and, if values
!   are not provided a message to that effect.
!    * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!    * SYMMETRIC STRUCTURE MEDIEVAL RUSSIAN TOWNS
!    *                    Key = RUSSIANT , Type = SSA
!    * No values provided - Information of pattern only
!    * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!
! +  dimension n, number of nonzero elements nnz, average number of
!    nonzero elements per column, standard deviation for this average.
! +  if the matrix is upper or lower triangular a message to that effect
!    is printed. Also the number of nonzeros in the strict upper
!    (lower) parts and the main diagonal are printed.
! +  weight of longest column. This is the largest number of nonzero
!    elements in a column encountered. Similarly for weight of
!    largest/smallest row.
! +  lower dandwidth as defined by
!          ml = max ( i-j, / all  a(i,j).ne. 0 )
! +  upper bandwidth as defined by
!          mu = max ( j-i, / all  a(i,j).ne. 0 )
!    NOTE that ml or mu can be negative. ml .lt. 0 would mean
!    that A is confined to the strict upper part above the diagonal
!    number -ml. Similarly for mu.
!
! +  maximun bandwidth as defined by
!    Max (  Max [ j ; a(i,j) .ne. 0 ] - Min [ j ; a(i,j) .ne. 0 ] )
!     i
! +  average bandwidth = average over all columns of the widths each column.
!
! +  If there are zero columns /or rows a message is printed
!    giving the number of such columns/rows.
!
! +  matching elements in A and transp(A) :this counts the number of
!    positions (i,j) such that if a(i,j) .ne. 0 then a(j,i) .ne. 0.
!    if this number is equal to nnz then the matrix is symmetric.
! +  Relative symmetry match : this is the ratio of the previous integer
!    over nnz. If this ratio is equal to one then the matrix has a
!    symmetric structure.
!
! +  average distance of a given element from the diagonal, standard dev.
!    the distance of a(i,j) is defined as iabs(j-i).
!
! +  Frobenius norm of A
!    Frobenius norm of 0.5*(A + transp(A))
!    Frobenius norm of 0.5*(A - transp(A))
!    these numbers provide information on the degree of symmetry
!    of the matrix. If the norm of the nonsymmetric part is
!    zero then the matrix is symmetric.
!
! + 90% of matrix is in the band of width k, means that
!   by moving away and in a symmetric manner from the main
!   diagonal you would have to include exactly k diagonals
!   (k is always odd), in order to include 90% of the nonzero
!   elements of A.  The same thing is then for 80%.
!
! + The total number of nonvoid diagonals, i.e., among the
!   2n-1 diagonals of the matrix which have at least one nonxero
!   element.
!
! +  Most important diagonals. The code selects a number of k
!    (k .le. 10) diagonals that are the most important ones, i.e.
!    that have the largest number of nonzero elements. Any diagonal
!    that has fewer than 1% of the nonzero elements of A is dropped.
!    the numbers printed are the offsets with respect to the
!    main diagonal, going from left tp right.
!    Thus 0 means the main diagonal -1 means the subdiagonal, and
!    +10 means the 10th upper diagonal.
! +  The accumulated percentages in the next line represent the
!    percentage of the nonzero elements represented by the diagonals
!    up the current one put together.
!    Thus:
!    *  The 10 most important diagonals are (offsets)    :             *
!    *     0     1     2    24    21     4    23    22    20    19     *
!    *  The accumulated percentages they represent are   :             *
!    *  40.4  68.1  77.7  80.9  84.0  86.2  87.2  88.3  89.4  90.4     *
!    *-----------------------------------------------------------------*
!    shows the offsets of the most important  diagonals and
!    40.4 represent ratio of the number of nonzero elements in the
!    diagonal zero (main diagonal) over the total number of nonzero
!    elements. the second number indicates that the diagonal 0 and the
!    diagonal 1 together hold 68.1% of the matrix, etc..
!
! +  Block structure:
!    if the matrix has a block structure then the block size is found
!    and printed. Otherwise the info1 will say that the matrix
!    does not have a block structure. Note that block struture has
!    a very specific meaning here. the matrix has a block structure
!    if it consists of square blocks that are dense. even if there
!    are zero elements in the blocks  they should be represented
!    otherwise it would be possible to determine the block size.
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
   DATA ipar1/1/
   WRITE (Iout,99003)
   WRITE (Iout,99001) Title(2:72) , Key , Type
99001 FORMAT (2x,' * ',a71,' *'/,2x,' *',20x,'Key = ',a8,' , Type = ',a3,25x,' *')
   IF ( .NOT.Valued ) WRITE (Iout,99002)
99002 FORMAT (2x,' * No values provided - Information on pattern only',23x,' *')
!---------------------------------------------------------------------
   nnz = Ia(N+1) - Ia(1)
   sym = ((Type(2:2)=='S') .OR. (Type(2:2)=='Z') .OR. (Type(2:2)=='s') .OR. (Type(2:2)=='z'))
!
   WRITE (Iout,99003)
   WRITE (Iout,99005) N , nnz
!-----------------------------------------------------------------------
99005 FORMAT (6x,' *  Dimension N                                      = ',i10,'  *'/6x,                                           &
             &' *  Number of nonzero elements                       = ',i10,'  *')
   job = 0
   IF ( Valued ) job = 1
   ipos = 1
   CALL csrcsc(N,job,ipos,A,Ja,Ia,Ao,Jao,Iao)
   CALL csrcsc(N,job,ipos,Ao,Jao,Iao,A,Ja,Ia)
!-------------------------------------------------------------------
! computing max bandwith, max number of nonzero elements per column
! min nonzero elements per column/row, row/column diagonal dominance
! occurences, average distance of an element from diagonal, number of
! elemnts in lower and upper parts, ...
!------------------------------------------------------------------
!    jao will be modified later, so we call skyline here
   CALL skyline(N,sym,Ja,Ia,Jao,Iao,nsky)
   CALL nonz_lud(N,Ja,Ia,nlower,nupper,ndiag)
   CALL avnz_col(N,Ja,Ia,av,st)
!------ write out info ----------------------------------------------
   IF ( sym ) nupper = nlower
   WRITE (Iout,99006) av , st
99006 FORMAT (6x,' *  Average number of nonzero elements/Column        = ',f10.4,'  *'/6x,                                         &
             &' *  Standard deviation for above average             = ',f10.4,'  *')
   IF ( nlower==0 ) WRITE (Iout,99013)
99013 FORMAT (6x,' *  The matrix is lower triangular ...       ',21x,' *')
   IF ( nupper==0 ) WRITE (Iout,99014)
99014 FORMAT (6x,' *  The matrix is upper triangular ...       ',21x,' *')
   WRITE (Iout,99015) nlower
99015 FORMAT (6x,' *  Nonzero elements in strict lower part            = ',i10,'  *')
   WRITE (Iout,99016) nupper
99016 FORMAT (6x,' *  Nonzero elements in strict upper part            = ',i10,'  *')
   WRITE (Iout,99017) ndiag
99017 FORMAT (6x,' *  Nonzero elements in main diagonal                = ',i10,'  *')
!
   CALL nonz(N,sym,Ja,Ia,Iao,nzmaxc,nzminc,nzmaxr,nzminr,nzcol,nzrow)
   WRITE (Iout,99007) nzmaxc , nzminc
!-----------------------------------------------------------------------
99007 FORMAT (6x,' *  Weight of longest column                         = ',i10,'  *'/6x,                                           &
             &' *  Weight of shortest column                        = ',i10,'  *')
!
   IF ( .NOT.sym ) WRITE (Iout,99008) nzmaxr , nzminr
99008 FORMAT (6x,' *  Weight of longest row                            = ',i10,'  *'/6x,                                           &
             &' *  Weight of shortest row                           = ',i10,'  *')
!
   IF ( nzcol/=0 ) WRITE (Iout,99024) nzcol
99024 FORMAT (6x,' *  There are zero columns. Number of such columns   = ',i10,'  *')
   IF ( nzrow/=0 ) WRITE (Iout,99023) nzrow
99023 FORMAT (6x,' *  There are zero rows. Number of such rows         = ',i10,'  *')
!
   CALL diag_domi(N,sym,Valued,A,Ja,Ia,Ao,Jao,Iao,ddomc,ddomr)
!-----------------------------------------------------------------------
! symmetry and near symmetry - Frobenius  norms
!-----------------------------------------------------------------------
   CALL frobnorm(N,sym,A,Ja,Ia,fnorm)
   CALL ansym(N,sym,A,Ja,Ia,Ao,Jao,Iao,imatch,av,fas,fan)
   CALL distaij(N,nnz,Ja,Ia,dist,std)
   amx = 0.0D0
   DO k = 1 , nnz
      amx = max(amx,abs(A(k)))
   ENDDO
   WRITE (Iout,99011) imatch , av , dist , std
99011 FORMAT (6x,' *  Matching elements in symmetry                    = ',i10,'  *'/6x,                                           &
             &' *  Relative Symmetry Match (symmetry=1)             = ',f10.4,'  *'/6x,                                            &
             &' *  Average distance of a(i,j)  from diag.           = ',e10.3,'  *'/6x,                                            &
             &' *  Standard deviation for above average             = ',e10.3,'  *')
   WRITE (Iout,99004)
   IF ( Valued ) THEN
      WRITE (Iout,99012) fnorm , fas , fan , amx , ddomr , ddomc
99012 FORMAT (6x,' *  Frobenius norm of A                              = ',e10.3,'  *'/6x,                                         &
             &' *  Frobenius norm of symmetric part                 = ',e10.3,'  *'/6x,                                            &
             &' *  Frobenius norm of nonsymmetric part              = ',e10.3,'  *'/6x,                                            &
             &' *  Maximum element in A                             = ',e10.3,'  *'/6x,                                            &
             &' *  Percentage of weakly diagonally dominant rows    = ',e10.3,'  *'/6x,                                            &
             &' *  Percentage of weakly diagonally dominant columns = ',e10.3,'  *')
      WRITE (Iout,99004)
   ENDIF
!-----------------------------------------------------------------------
!--------------------bandedness- main diagonals ----------------------- -
!-----------------------------------------------------------------------
   n2 = N + N - 1
   DO i = 1 , n2
      Jao(i) = 0
   ENDDO
   DO i = 1 , N
      k1 = Ia(i)
      k2 = Ia(i+1) - 1
      DO k = k1 , k2
         j = Ja(k)
         Jao(N+i-j) = Jao(N+i-j) + 1
      ENDDO
   ENDDO
!

   CALL bandwidth(N,Ja,Ia,ml,mu,iband,bndav)

!
!     write bandwidth information .
!
   WRITE (Iout,99009) ml , mu , iband , bndav
99009 FORMAT (6x,' *  Lower bandwidth  (max: i-j, a(i,j) .ne. 0)       = ',i10,'  *'/6x,                                           &
             &' *  Upper bandwidth  (max: j-i, a(i,j) .ne. 0)       = ',i10,'  *'/6x,                                              &
             &' *  Maximum Bandwidth                                = ',i10,'  *'/6x,                                              &
             &' *  Average Bandwidth                                = ',e10.3,'  *')
!
   WRITE (Iout,99010) nsky
99010 FORMAT (6x,' *  Number of nonzeros in skyline storage            = ',i10,'  *')
!
!         call percentage_matrix(n,nnz,ja,ia,jao,90,jb2)
!         call percentage_matrix(n,nnz,ja,ia,jao,80,jb1)
   nrow = N
   ncol = N
   CALL distdiag(nrow,ncol,Ja,Ia,Jao)
   CALL bandpart(N,Ia,Jao,90,jb2)
   CALL bandpart(N,Ia,Jao,80,jb1)
   WRITE (Iout,99020) 2*jb2 + 1 , 2*jb1 + 1
! 111	format(
!     * 6x,' *  They constitute the following % of A             = ',
!     * f8.1,' %  *')
99020 FORMAT (6x,' *  90% of matrix is in the band of width            = ',i10,'  *',/,6x,                                         &
             &' *  80% of matrix is in the band of width            = ',i10,'  *')
!-----------------------------------------------------------------
   nzdiag = 0
   n2 = N + N - 1
   DO i = 1 , n2
      IF ( Jao(i)/=0 ) nzdiag = nzdiag + 1
   ENDDO
   CALL n_imp_diag(N,nnz,Jao,ipar1,ndiag,ioff,dcount)
   WRITE (Iout,99025) nzdiag
99025 FORMAT (6x,' *  The total number of nonvoid diagonals is         = ',i10,'  *')
   WRITE (tmpst,'(10i6)') (ioff(j),j=1,ndiag)
   WRITE (Iout,99018) ndiag , tmpst
99018 FORMAT (6x,' *  The ',i2,' most important',' diagonals are (offsets)    : ',10x,'  *',/,6x,' *',a61,3x,' *')
   WRITE (tmpst,'(10f6.1)') (dcount(j),j=1,ndiag)
   WRITE (Iout,99019) tmpst
99019 FORMAT (6x,' *  The accumulated percentages they represent are ','  : ',10x,'  *',/,6x,' *',a61,3x,' *')
   WRITE (Iout,99004)
!     jump to next page -- optional //
!     write (iout,'(1h1)')
!-----------------------------------------------------------------------
!     determine block size if matrix is a block matrix..
!-----------------------------------------------------------------------
   CALL blkfnd(N,Ja,Ia,nblk)
   IF ( nblk<=1 ) THEN
      WRITE (Iout,99021)
99021 FORMAT (6x,' *  The matrix does not have a block structure ',19x,' *')
   ELSE
      WRITE (Iout,99022) nblk
99022 FORMAT (6x,' *  Block structure found with block size            = ',i10,'  *')
   ENDIF
   WRITE (Iout,99004)
!
!---------- done. Next define all the formats --------------------------
!
99003 FORMAT (2x,38(' *'))
99004 FORMAT (6x,' *',65('-'),'*')
!-------------------------- end of dinfo --------------------------
END SUBROUTINE dinfo1
