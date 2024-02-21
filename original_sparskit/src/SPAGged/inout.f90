!*==readmt.f90 processed by SPAG 8.04RA 15:47 21 Feb 2024
!!SPAG Open source Personal, Educational or Academic User Clemson University  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
!-----------------------------------------------------------------------
!-----------------------end-of-vbrinfo----------------------------------
!----------------------------------------------------------------------c
!                          S P A R S K I T                             c
!----------------------------------------------------------------------c
!                        INPUT-OUTPUT MODULE                           c
!----------------------------------------------------------------------c
! contents:                                                            c
!----------                                                            c
!  readmt : reads matrices in the Boeing/Harwell format.               c
!  prtmt  : prints matrices in the Boeing/Harwell format.              c
!  dump   : outputs matrix rows in a simple format (debugging purposes)c
!  pspltm : generates a post-script plot of the non-zero pattern of A  c
!  pltmt  : produces a 'pic' file for plotting a sparse matrix         c
!  smms   : write the matrx in a format used in SMMS package           c
!  readsm : reads matrics in coordinate format (as in SMMS package)    c
!  readsk : reads matrices in CSR format (simplified H/B formate).     c
!  skit   : writes matrics to a file, format same as above.            c
!  prtunf : writes matrics (in CSR format) unformatted                 c
!  readunf: reads unformatted data of matrics (in CSR format)          c
!----------------------------------------------------------------------c
SUBROUTINE readmt(Nmax,Nzmax,Job,Iounit,A,Ja,Ia,Rhs,Nrhs,Guesol,Nrow,Ncol,Nnz,Title,Key,Type,Ierr)
   USE ISO_FORTRAN_ENV                 
   IMPLICIT NONE
!
! Dummy argument declarations rewritten by SPAG
!
   INTEGER , INTENT(IN) :: Nmax
   INTEGER , INTENT(IN) :: Nzmax
   INTEGER , INTENT(INOUT) :: Job
   INTEGER , INTENT(IN) :: Iounit
   REAL(REAL64) , INTENT(OUT) , DIMENSION(Nzmax) :: A
   INTEGER , INTENT(OUT) , DIMENSION(Nzmax) :: Ja
   INTEGER , INTENT(OUT) , DIMENSION(Nmax+1) :: Ia
   REAL(REAL64) , INTENT(OUT) , DIMENSION(*) :: Rhs
   INTEGER , INTENT(INOUT) :: Nrhs
   CHARACTER(2) , INTENT(INOUT) :: Guesol
   INTEGER , INTENT(INOUT) :: Nrow
   INTEGER , INTENT(INOUT) :: Ncol
   INTEGER , INTENT(INOUT) :: Nnz
   CHARACTER(72) , INTENT(OUT) :: Title
   CHARACTER(8) , INTENT(OUT) :: Key
   CHARACTER(3) , INTENT(OUT) :: Type
   INTEGER , INTENT(INOUT) :: Ierr
!
! Local variable declarations rewritten by SPAG
!
   INTEGER :: i , iend , indcrd , len , lenrhs , n , neltvl , next , nrwindx , nvec , ptrcrd , rhscrd , totcrd , valcrd
   CHARACTER(16) :: indfmt , ptrfmt
   CHARACTER(20) :: rhsfmt , valfmt
   CHARACTER(3) :: rhstyp
!
! End of declarations rewritten by SPAG
!
!-----------------------------------------------------------------------
! this subroutine reads  a boeing/harwell matrix. handles right hand
! sides in full format only (no sparse right hand sides).
! Also the matrix must be in assembled forms.
! Author: Youcef Saad - Date: Sept., 1989
!         updated Oct 31, 1989.
!-----------------------------------------------------------------------
! on entry:
!---------
! nmax 	 =  max column dimension  allowed for matrix. The array ia should
!	    be of length at least ncol+1 (see below) if job.gt.0
! nzmax	 = max number of nonzeros elements allowed. the arrays a,
!          and ja should be of length equal to nnz (see below) if these
!          arrays are to be read (see job).
!
! job	 = integer to indicate what is to be read. (note: job is an
!          input and output parameter, it can be modified on return)
!          job = 0    read the values of ncol, nrow, nnz, title, key,
!                     type and return. matrix is not read and arrays
!                     a, ja, ia, rhs are not touched.
!          job = 1    read srtucture only, i.e., the arrays ja and ia.
!          job = 2    read matrix including values, i.e., a, ja, ia
!          job = 3    read matrix and right hand sides: a,ja,ia,rhs.
!		      rhs may contain initial guesses and exact
!                     solutions appended to the actual right hand sides.
!		      this will be indicated by the output parameter
!                     guesol [see below].
!
! nrhs   = integer. nrhs is an input as well as ouput parameter.
!          at input nrhs contains the total length of the array rhs.
!          See also ierr and nrhs in output parameters.
!
! iounit = logical unit number where to read the matrix from.
!
! on return:
!----------
! job    = on return job may be modified to the highest job it could
!          do: if job=2 on entry but no matrix values are available it
!          is reset to job=1 on return. Similarly of job=3 but no rhs
!          is provided then it is rest to job=2 or job=1 depending on
!          whether or not matrix values are provided.
!          Note that no error message is triggered (i.e. ierr = 0
!          on return in these cases. It is therefore important to
!          compare the values of job on entry and return ).
!
! a	 = the a matrix in the a, ia, ja (column) storage format
! ja 	 = row number of element a(i,j) in array a.
! ia     = pointer  array. ia(i) points to the beginning of column i.
!
! rhs    = real array of size nrow + 1 if available (see job)
!
! nrhs   = integer containing the number of right-hand sides found
!          each right hand side may be accompanied with an intial guess
!          and also the exact solution.
!
! guesol = a 2-character string indicating whether an initial guess
!          (1-st character) and / or the exact solution (2-nd
!          character) is provided with the right hand side.
!	   if the first character of guesol is 'G' it means that an
!          an intial guess is provided for each right-hand side.
!          These are appended to the right hand-sides in the array rhs.
!	   if the second character of guesol is 'X' it means that an
!          exact solution is provided for each right-hand side.
!          These are  appended to the right hand-sides
!          and the initial guesses (if any) in the array rhs.
!
! nrow   = number of rows in matrix
! ncol	 = number of columns in matrix
! nnz	 = number of nonzero elements in A. This info is returned
!          even if there is not enough space in a, ja, ia, in order
!          to determine the minimum storage needed.
!
! title  = character*72 = title of matrix test ( character a*72).
! key    = character*8  = key of matrix
! type   = charatcer*3  = type of matrix.
!          for meaning of title, key and type refer to documentation
!          Harwell/Boeing matrices.
!
! ierr   = integer used for error messages
!         * ierr  =  0 means that  the matrix has been read normally.
!         * ierr  =  1 means that  the array matrix could not be read
!         because ncol+1 .gt. nmax
!         * ierr  =  2 means that  the array matrix could not be read
!         because nnz .gt. nzmax
!         * ierr  =  3 means that  the array matrix could not be read
!         because both (ncol+1 .gt. nmax) and  (nnz .gt. nzmax )
!         * ierr  =  4 means that  the right hand side (s) initial
!         guesse (s) and exact solution (s)   could  not be
!         read because they are stored in sparse format (not handled
!         by this routine ...)
!         * ierr  =  5 means that the right-hand-sides, initial guesses
!         and exact solutions could not be read because the length of
!         rhs as specified by the input value of nrhs is not
!         sufficient to store them. The rest of the matrix may have
!         been read normally.
!
! Notes:
!-------
! 1) The file inout must be open (and possibly rewound if necessary)
!    prior to calling readmt.
! 2) Refer to the documentation on the Harwell-Boeing formats
!    for details on the format assumed by readmt.
!    We summarize the format here for convenience.
!
!    a) all lines in inout are assumed to be 80 character long.
!    b) the file consists of a header followed by the block of the
!       column start pointers followed by the block of the
!       row indices, followed by the block of the real values and
!       finally the numerical values of the right-hand-side if a
!       right hand side is supplied.
!    c) the file starts by a header which contains four lines if no
!       right hand side is supplied and five lines otherwise.
!       * first line contains the title (72 characters long) followed by
!         the 8-character identifier (name of the matrix, called key)
!        [ A72,A8 ]
!       * second line contains the number of lines for each
!         of the following data blocks (4 of them) and the total number
!         of lines excluding the header.
!        [5i4]
!       * the third line contains a three character string identifying
!         the type of matrices as they are referenced in the Harwell
!         Boeing documentation [e.g., rua, rsa,..] and the number of
!         rows, columns, nonzero entries.
!         [A3,11X,4I14]
!       * The fourth line contains the variable fortran format
!         for the following data blocks.
!         [2A16,2A20]
!       * The fifth line is only present if right-hand-sides are
!         supplied. It consists of three one character-strings containing
!         the storage format for the right-hand-sides
!         ('F'= full,'M'=sparse=same as matrix), an initial guess
!         indicator ('G' for yes), an exact solution indicator
!         ('X' for yes), followed by the number of right-hand-sides
!         and then the number of row indices.
!         [A3,11X,2I14]
!     d) The three following blocks follow the header as described
!        above.
!     e) In case the right hand-side are in sparse formats then
!        the fourth block uses the same storage format as for the matrix
!        to describe the NRHS right hand sides provided, with a column
!        being replaced by a right hand side.
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
   Ierr = 0
   lenrhs = Nrhs
!
   READ (Iounit,99001) Title , Key , totcrd , ptrcrd , indcrd , valcrd , rhscrd , Type , Nrow , Ncol , Nnz , neltvl , ptrfmt ,     &
                     & indfmt , valfmt , rhsfmt
99001 FORMAT (a72,a8/5I14/a3,11x,4I14/2A16,2A20)
!
   IF ( rhscrd>0 ) READ (Iounit,99002) rhstyp , Nrhs , nrwindx
99002 FORMAT (a3,11x,i14,i14)
!
! anything else to read ?
!
   IF ( Job<=0 ) RETURN
!     ---- check whether matrix is readable ------
   n = Ncol
   IF ( Ncol>Nmax ) Ierr = 1
   IF ( Nnz>Nzmax ) Ierr = Ierr + 2
   IF ( Ierr/=0 ) RETURN
!     ---- read pointer and row numbers ----------
   READ (Iounit,ptrfmt) (Ia(i),i=1,n+1)
   READ (Iounit,indfmt) (Ja(i),i=1,Nnz)
!     --- reading values of matrix if required....
   IF ( Job<=1 ) RETURN
!     --- and if available -----------------------
   IF ( valcrd<=0 ) THEN
      Job = 1
      RETURN
   ENDIF
   READ (Iounit,valfmt) (A(i),i=1,Nnz)
!     --- reading rhs if required ----------------
   IF ( Job<=2 ) RETURN
!     --- and if available -----------------------
   IF ( rhscrd<=0 ) THEN
      Job = 2
      Nrhs = 0
      RETURN
   ENDIF
!
!     --- read right-hand-side.--------------------
!
   IF ( rhstyp(1:1)=='M' ) THEN
      Ierr = 4
      RETURN
   ENDIF
!
   Guesol = rhstyp(2:3)
!
   nvec = 1
   IF ( Guesol(1:1)=='G' .OR. Guesol(1:1)=='g' ) nvec = nvec + 1
   IF ( Guesol(2:2)=='X' .OR. Guesol(2:2)=='x' ) nvec = nvec + 1
!
   len = Nrhs*Nrow
!
   IF ( len*nvec>lenrhs ) THEN
      Ierr = 5
      RETURN
   ENDIF
!
! read right-hand-sides
!
   next = 1
   iend = len
   READ (Iounit,rhsfmt) (Rhs(i),i=next,iend)
!
! read initial guesses if available
!
   IF ( Guesol(1:1)=='G' .OR. Guesol(1:1)=='g' ) THEN
      next = next + len
      iend = iend + len
      READ (Iounit,valfmt) (Rhs(i),i=next,iend)
   ENDIF
!
! read exact solutions if available
!
   IF ( Guesol(2:2)=='X' .OR. Guesol(2:2)=='x' ) THEN
      next = next + len
      iend = iend + len
      READ (Iounit,valfmt) (Rhs(i),i=next,iend)
   ENDIF
!
!--------- end of readmt -----------------------------------------------
!-----------------------------------------------------------------------
END SUBROUTINE readmt
!*==prtmt.f90 processed by SPAG 8.04RA 15:47 21 Feb 2024
!!SPAG Open source Personal, Educational or Academic User Clemson University  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
!-----------------------------------------------------------------------
SUBROUTINE prtmt(Nrow,Ncol,A,Ja,Ia,Rhs,Guesol,Title,Key,Type,Ifmt,Job,Iounit)
   USE ISO_FORTRAN_ENV                 
   IMPLICIT NONE
!
! Dummy argument declarations rewritten by SPAG
!
   INTEGER , INTENT(IN) :: Nrow
   INTEGER , INTENT(IN) :: Ncol
   REAL(REAL64) , INTENT(IN) , DIMENSION(*) :: A
   INTEGER , INTENT(IN) , DIMENSION(*) :: Ja
   INTEGER , INTENT(IN) , DIMENSION(*) :: Ia
   REAL(REAL64) , INTENT(IN) , DIMENSION(*) :: Rhs
   CHARACTER(2) , INTENT(IN) :: Guesol
   CHARACTER(72) , INTENT(IN) :: Title
   CHARACTER(8) , INTENT(IN) :: Key
   CHARACTER(3) , INTENT(IN) :: Type
   INTEGER , INTENT(INOUT) :: Ifmt
   INTEGER , INTENT(IN) :: Job
   INTEGER , INTENT(IN) :: Iounit
!
! Local variable declarations rewritten by SPAG
!
   INTEGER :: i , iend , ihead , indcrd , len , next , nnz , nperli , nrhs , nrwindx , ptrcrd , rhscrd , totcrd , valcrd
   CHARACTER(16) :: indfmt , ptrfmt
   CHARACTER(28) :: ix
   CHARACTER(3) :: rhstyp
   CHARACTER(20) :: valfmt
!
! End of declarations rewritten by SPAG
!
!-----------------------------------------------------------------------
! writes a matrix in Harwell-Boeing format into a file.
! assumes that the matrix is stored in COMPRESSED SPARSE COLUMN FORMAT.
! some limited functionality for right hand sides.
! Author: Youcef Saad - Date: Sept., 1989 - updated Oct. 31, 1989 to
! cope with new format.
!-----------------------------------------------------------------------
! on entry:
!---------
! nrow   = number of rows in matrix
! ncol	 = number of columns in matrix
! a	 = real*8 array containing the values of the matrix stored
!          columnwise
! ja 	 = integer array of the same length as a containing the column
!          indices of the corresponding matrix elements of array a.
! ia     = integer array of containing the pointers to the beginning of
!	   the row in arrays a and ja.
! rhs    = real array  containing the right-hand-side (s) and optionally
!          the associated initial guesses and/or exact solutions
!          in this order. See also guesol for details. the vector rhs will
!          be used only if job .gt. 2 (see below). Only full storage for
!          the right hand sides is supported.
!
! guesol = a 2-character string indicating whether an initial guess
!          (1-st character) and / or the exact solution (2-nd)
!          character) is provided with the right hand side.
!	   if the first character of guesol is 'G' it means that an
!          an intial guess is provided for each right-hand sides.
!          These are assumed to be appended to the right hand-sides in
!          the array rhs.
!	   if the second character of guesol is 'X' it means that an
!          exact solution is provided for each right-hand side.
!          These are assumed to be appended to the right hand-sides
!          and the initial guesses (if any) in the array rhs.
!
! title  = character*72 = title of matrix test ( character a*72 ).
! key    = character*8  = key of matrix
! type   = charatcer*3  = type of matrix.
!
! ifmt	 = integer specifying the format chosen for the real values
!	   to be output (i.e., for a, and for rhs-guess-sol if
!          applicable). The meaning of ifmt is as follows.
!	  * if (ifmt .lt. 100) then the D descriptor is used,
!           format Dd.m, in which the length (m) of the mantissa is
!           precisely the integer ifmt (and d = ifmt+6)
!	  * if (ifmt .gt. 100) then prtmt will use the
!           F- descriptor (format Fd.m) in which the length of the
!           mantissa (m) is the integer mod(ifmt,100) and the length
!           of the integer part is k=ifmt/100 (and d = k+m+2)
!	    Thus  ifmt= 4   means  D10.4  +.xxxxD+ee    while
!	          ifmt=104  means  F7.4   +x.xxxx
!	          ifmt=205  means  F9.5   +xx.xxxxx
!	    Note: formats for ja, and ia are internally computed.
!
! job	 = integer to indicate whether matrix values and
!	   a right-hand-side is available to be written
!          job = 1   write srtucture only, i.e., the arrays ja and ia.
!          job = 2   write matrix including values, i.e., a, ja, ia
!          job = 3   write matrix and one right hand side: a,ja,ia,rhs.
!	   job = nrhs+2 write matrix and nrhs successive right hand sides
!	   Note that there cannot be any right-hand-side if the matrix
!	   has no values. Also the initial guess and exact solutions when
!          provided are for each right hand side. For example if nrhs=2
!          and guesol='GX' there are 6 vectors to write.
!
!
! iounit = logical unit number where to write the matrix into.
!
! on return:
!----------
! the matrix a, ja, ia will be written in output unit iounit
! in the Harwell-Boeing format. None of the inputs is modofied.
!
! Notes: 1) This code attempts to pack as many elements as possible per
!        80-character line.
!        2) this code attempts to avoid as much as possible to put
!        blanks in the formats that are written in the 4-line header
!         (This is done for purely esthetical reasons since blanks
!        are ignored in format descriptors.)
!        3) sparse formats for right hand sides and guesses are not
!        supported.
!-----------------------------------------------------------------------
!--------------
!     compute pointer format
!--------------
   nnz = Ia(Ncol+1) - 1
   IF ( nnz==0 ) RETURN
   len = int(alog10(0.1+real(nnz+1))) + 1
   nperli = 80/len
   ptrcrd = Ncol/nperli + 1
   IF ( len>9 ) THEN
!         assign 101 to ix
      ix = '(1h(,i2,1HI,i2,1h) )'
   ELSE
      ix = '(1h(,i2,1HI,i1,1h) )'
   ENDIF
   WRITE (ptrfmt,ix) nperli , len
! 100  format(1h(,i2,1HI,i1,1h) )
! 101  format(1h(,i2,1HI,i2,1h) )
!----------------------------
! compute ROW index format
!----------------------------
   len = int(alog10(0.1+real(Nrow))) + 1
   nperli = min0(80/len,nnz)
   indcrd = (nnz-1)/nperli + 1
   WRITE (indfmt,'(1h(,i2,1HI,i1,1h))') nperli , len
!---------------
! compute values and rhs format (using the same for both)
!---------------
   valcrd = 0
   rhscrd = 0
! quit this part if no values provided.
   IF ( Job>1 ) THEN
!
      IF ( Ifmt>=100 ) THEN
         ihead = Ifmt/100
         Ifmt = Ifmt - 100*ihead
         len = ihead + Ifmt + 2
         nperli = 80/len
!
         IF ( len<=9 ) THEN
            ix = '(1h(,i2,1hF,i1,1h.,i1,1h) )'
         ELSEIF ( Ifmt<=9 ) THEN
            ix = '(1h(,i2,1hF,i2,1h.,i1,1h) )'
         ELSE
            ix = '(1h(,i2,1hF,i2,1h.,i2,1h) )'
         ENDIF
!
         WRITE (valfmt,ix) nperli , len , Ifmt
      ELSE
         len = Ifmt + 6
         nperli = 80/len
!     try to minimize the blanks in the format strings.
         IF ( nperli<=9 ) THEN
            IF ( len<=9 ) THEN
               ix = '(1h(,i1,1hD,i1,1h.,i1,1h) )'
            ELSEIF ( Ifmt<=9 ) THEN
               ix = '(1h(,i1,1hD,i2,1h.,i1,1h) )'
            ELSE
               ix = '(1h(,i1,1hD,i2,1h.,i2,1h) )'
            ENDIF
         ELSEIF ( len<=9 ) THEN
            ix = '(1h(,i2,1hD,i1,1h.,i1,1h) )'
         ELSEIF ( Ifmt<=9 ) THEN
            ix = '(1h(,i2,1hD,i2,1h.,i1,1h) )'
         ELSE
            ix = '(1h(,i2,1hD,i2,1h.,i2,1h) )'
         ENDIF
!-----------
         WRITE (valfmt,ix) nperli , len , Ifmt
!105     format(1h(,i1,1hD,i1,1h.,i1,1h) )
!106     format(1h(,i1,1hD,i2,1h.,i1,1h) )
!107     format(1h(,i1,1hD,i2,1h.,i2,1h) )
!108     format(1h(,i2,1hD,i1,1h.,i1,1h) )
!109     format(1h(,i2,1hD,i2,1h.,i1,1h) )
!110     format(1h(,i2,1hD,i2,1h.,i2,1h) )
!
      ENDIF
      valcrd = (nnz-1)/nperli + 1
      nrhs = Job - 2
      IF ( nrhs>=1 ) THEN
         i = (nrhs*Nrow-1)/nperli + 1
         rhscrd = i
         IF ( Guesol(1:1)=='G' .OR. Guesol(1:1)=='g' ) rhscrd = rhscrd + i
         IF ( Guesol(2:2)=='X' .OR. Guesol(2:2)=='x' ) rhscrd = rhscrd + i
         rhstyp = 'F'//Guesol
      ENDIF
   ENDIF
!
   totcrd = ptrcrd + indcrd + valcrd + rhscrd
!     write 4-line or five line header
   WRITE (Iounit,99001) Title , Key , totcrd , ptrcrd , indcrd , valcrd , rhscrd , Type , Nrow , Ncol , nnz , nrhs , ptrfmt ,      &
                      & indfmt , valfmt , valfmt
99001 FORMAT (a72,a8/5I14/a3,11x,4I14/2A16,2A20)
!-----------------------------------------------------------------------
   nrwindx = 0
   IF ( nrhs>=1 ) WRITE (Iounit,99002) rhstyp , nrhs , nrwindx
99002 FORMAT (A3,11x,i14,i14)
!
   WRITE (Iounit,ptrfmt) (Ia(i),i=1,Ncol+1)
   WRITE (Iounit,indfmt) (Ja(i),i=1,nnz)
   IF ( Job<=1 ) RETURN
   WRITE (Iounit,valfmt) (A(i),i=1,nnz)
   IF ( Job<=2 ) RETURN
   len = Nrow*nrhs
   next = 1
   iend = len
   WRITE (Iounit,valfmt) (Rhs(i),i=next,iend)
!
!     write initial guesses if available
!
   IF ( Guesol(1:1)=='G' .OR. Guesol(1:1)=='g' ) THEN
      next = next + len
      iend = iend + len
      WRITE (Iounit,valfmt) (Rhs(i),i=next,iend)
   ENDIF
!
!     write exact solutions if available
!
   IF ( Guesol(2:2)=='X' .OR. Guesol(2:2)=='x' ) THEN
      next = next + len
      iend = iend + len
      WRITE (Iounit,valfmt) (Rhs(i),i=next,iend)
   ENDIF
!
!----------end of prtmt ------------------------------------------------
!-----------------------------------------------------------------------
END SUBROUTINE prtmt
!*==dump.f90 processed by SPAG 8.04RA 15:47 21 Feb 2024
!!SPAG Open source Personal, Educational or Academic User Clemson University  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
!-----------------------------------------------------------------------
SUBROUTINE dump(I1,I2,Values,A,Ja,Ia,Iout)
   USE ISO_FORTRAN_ENV                 
   IMPLICIT NONE
!
! Dummy argument declarations rewritten by SPAG
!
   INTEGER , INTENT(IN) :: I1
   INTEGER , INTENT(IN) :: I2
   LOGICAL , INTENT(IN) :: Values
   REAL(REAL64) , INTENT(IN) , DIMENSION(*) :: A
   INTEGER , INTENT(IN) , DIMENSION(*) :: Ja
   INTEGER , INTENT(IN) , DIMENSION(*) :: Ia
   INTEGER , INTENT(IN) :: Iout
!
! Local variable declarations rewritten by SPAG
!
   INTEGER :: i , k , k1 , k2 , maxr
!
! End of declarations rewritten by SPAG
!
!-----------------------------------------------------------------------
! outputs rows i1 through i2 of a sparse matrix stored in CSR format
! (or columns i1 through i2 of a matrix stored in CSC format) in a file,
! one (column) row at a time in a nice readable format.
! This is a simple routine which is useful for debugging.
!-----------------------------------------------------------------------
! on entry:
!---------
! i1    = first row (column) to print out
! i2    = last row (column) to print out
! values= logical. indicates whether or not to print real values.
!         if value = .false. only the pattern will be output.
! a,
! ja,
! ia    =  matrix in CSR format (or CSC format)
! iout  = logical unit number for output.
!----------
! the output file iout will have written in it the rows or columns
! of the matrix in one of two possible formats (depending on the max
! number of elements per row. The values are output with only
! two digits of accuracy (D9.2). )
!-----------------------------------------------------------------------
!     local variables
!
! select mode horizontal or vertical
!
   maxr = 0
   DO i = I1 , I2
      maxr = max0(maxr,Ia(i+1)-Ia(i))
   ENDDO
 
   IF ( maxr<=8 ) THEN
!
! able to do one row acros line
!
      DO i = I1 , I2
         WRITE (Iout,99001) i
!
! formats :
!
99001    FORMAT (' ',34('-'),' row',i6,1x,34('-'))
         k1 = Ia(i)
         k2 = Ia(i+1) - 1
         WRITE (Iout,99002) (Ja(k),k=k1,k2)
99002    FORMAT (' col:',8(i5,'     :'))
         IF ( Values ) WRITE (Iout,99003) (A(k),k=k1,k2)
99003    FORMAT (' val:',8(D9.2,' :'))
      ENDDO
   ELSE
!
! unable to one row acros line. do three items at a time
! across a line
      DO i = I1 , I2
         IF ( Values ) THEN
            WRITE (Iout,99004) i
99004       FORMAT (' ',30('-'),' row',i3,1x,30('-'),/3('  columns :    values  * '))
         ELSE
            WRITE (Iout,99007) i
99007       FORMAT (' ',30('-'),' row',i3,1x,30('-'),/3('  column  :  column   *'))
         ENDIF
         k1 = Ia(i)
         k2 = Ia(i+1) - 1
         IF ( Values ) THEN
            WRITE (Iout,99005) (Ja(k),A(k),k=k1,k2)
!-------------xiiiiiihhhhhhddddddddd-*-
99005       FORMAT (3(' ',i6,'   :  ',D9.2,' * '))
         ELSE
            WRITE (Iout,99006) (Ja(k),k=k1,k2)
99006       FORMAT (6(' ',i5,'  *   '))
         ENDIF
      ENDDO
   ENDIF
!----end-of-dump--------------------------------------------------------
!-----------------------------------------------------------------------
END SUBROUTINE dump
!*==pspltm.f90 processed by SPAG 8.04RA 15:47 21 Feb 2024
!!SPAG Open source Personal, Educational or Academic User Clemson University  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
!-----------------------------------------------------------------------
SUBROUTINE pspltm(Nrow,Ncol,Mode,Ja,Ia,Title,Ptitle,Size,Munt,Nlines,Lines,Iunt)
   IMPLICIT NONE
!
! Dummy argument declarations rewritten by SPAG
!
   INTEGER , INTENT(IN) :: Nlines
   INTEGER , INTENT(IN) :: Nrow
   INTEGER , INTENT(IN) :: Ncol
   INTEGER , INTENT(IN) :: Mode
   INTEGER , INTENT(IN) , DIMENSION(*) :: Ja
   INTEGER , INTENT(IN) , DIMENSION(*) :: Ia
   CHARACTER(*) :: Title
   INTEGER , INTENT(IN) :: Ptitle
   REAL , INTENT(IN) :: Size
   CHARACTER(2) , INTENT(IN) :: Munt
   INTEGER , INTENT(IN) , DIMENSION(Nlines) :: Lines
   INTEGER , INTENT(IN) :: Iunt
!
! Local variable declarations rewritten by SPAG
!
   REAL :: botmrgn , delt , fnstit , frlw , lrmrgn , paperx , scfct , siz , u2dot , xl , xr , xtit , xx , yb , yt , ytit , ytitof ,&
         & yy
   REAL , SAVE :: conv , haf , zero
   INTEGER :: ii , ilast , isep , istart , k , kol , ltit , m , maxdim , n , nc , nr
   INTEGER , EXTERNAL :: lenstr
   LOGICAL , SAVE :: square
!
! End of declarations rewritten by SPAG
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
! PSPLTM - PostScript PLoTer of a (sparse) Matrix
! This version by loris renggli (renggli@masg1.epfl.ch), Dec 1991
! and Youcef Saad
!------
! Loris RENGGLI, Swiss Federal Institute of Technology, Math. Dept
! CH-1015 Lausanne (Switzerland)  -- e-mail:  renggli@masg1.epfl.ch
! Modified by Youcef Saad -- June 24, 1992 to add a few features:
! separation lines + acceptance of MSR format.
!-----------------------------------------------------------------------
! input arguments description :
!
! nrow   = number of rows in matrix
!
! ncol   = number of columns in matrix
!
! mode   = integer indicating whether the matrix is stored in
!           CSR mode (mode=0) or CSC mode (mode=1) or MSR mode (mode=2)
!
! ja     = column indices of nonzero elements when matrix is
!          stored rowise. Row indices if stores column-wise.
! ia     = integer array of containing the pointers to the
!          beginning of the columns in arrays a, ja.
!
! title  = character*(*). a title of arbitrary length to be printed
!          as a caption to the figure. Can be a blank character if no
!          caption is desired.
!
! ptitle = position of title; 0 under the drawing, else above
!
! size   = size of the drawing
!
! munt   = units used for size : 'cm' or 'in'
!
! nlines = number of separation lines to draw for showing a partionning
!          of the matrix. enter zero if no partition lines are wanted.
!
! lines  = integer array of length nlines containing the coordinates of
!          the desired partition lines . The partitioning is symmetric:
!          a horizontal line across the matrix will be drawn in
!          between rows lines(i) and lines(i)+1 for i=1, 2, ..., nlines
!          an a vertical line will be similarly drawn between columns
!          lines(i) and lines(i)+1 for i=1,2,...,nlines
!
! iunt   = logical unit number where to write the matrix into.
!-----------------------------------------------------------------------
! additional note: use of 'cm' assumes european format for paper size
! (21cm wide) and use of 'in' assumes american format (8.5in wide).
! The correct centering of the figure depends on the proper choice. Y.S.
!-----------------------------------------------------------------------
! external
! local variables ---------------------------------------------------
! change square to .true. if you prefer a square frame around
! a rectangular matrix
   DATA haf/0.5/ , zero/0.0/ , conv/2.54/ , square/.FALSE./
!-----------------------------------------------------------------------
   siz = Size
   nr = Nrow
   nc = Ncol
   n = nc
   IF ( Mode==0 ) n = nr
!      nnz = ia(n+1) - ia(1)
   maxdim = max(Nrow,Ncol)
   m = 1 + maxdim
   nc = nc + 1
   nr = nr + 1
!
! units (cm or in) to dot conversion factor and paper size
!
   IF ( Munt=='cm' .OR. Munt=='CM' ) THEN
      u2dot = 72.0/conv
      paperx = 21.0
   ELSE
      u2dot = 72.0
      paperx = 8.5*conv
      siz = siz*conv
   ENDIF
!
! left and right margins (drawing is centered)
!
   lrmrgn = (paperx-siz)/2.0
!
! bottom margin : 2 cm
!
   botmrgn = 2.0
! scaling factor
   scfct = siz*u2dot/m
! matrix frame line witdh
   frlw = 0.25
! font size for title (cm)
   fnstit = 0.5
   ltit = lenstr(Title)
! position of title : centered horizontally
!                     at 1.0 cm vertically over the drawing
   ytitof = 1.0
   xtit = paperx/2.0
   ytit = botmrgn + siz*nr/m + ytitof
! almost exact bounding box
   xl = lrmrgn*u2dot - scfct*frlw/2
   xr = (lrmrgn+siz)*u2dot + scfct*frlw/2
   yb = botmrgn*u2dot - scfct*frlw/2
   yt = (botmrgn+siz*nr/m)*u2dot + scfct*frlw/2
   IF ( ltit>0 ) yt = yt + (ytitof+fnstit*0.70)*u2dot
! add some room to bounding box
   delt = 10.0
   xl = xl - delt
   xr = xr + delt
   yb = yb - delt
   yt = yt + delt
!
! correction for title under the drawing
   IF ( Ptitle==0 .AND. ltit>0 ) THEN
      ytit = botmrgn + fnstit*0.3
      botmrgn = botmrgn + ytitof + fnstit*0.7
   ENDIF
! begin of output
!
   WRITE (Iunt,99001) '%!'
   WRITE (Iunt,99001) '%%Creator: PSPLTM routine'
   WRITE (Iunt,99003) '%%BoundingBox:' , xl , yb , xr , yt
99003 FORMAT (A,4(1x,F9.2))
   WRITE (Iunt,99001) '%%EndComments'
   WRITE (Iunt,99001) '/cm {72 mul 2.54 div} def'
   WRITE (Iunt,99001) '/mc {72 div 2.54 mul} def'
   WRITE (Iunt,99001) '/pnum { 72 div 2.54 mul 20 string'
   WRITE (Iunt,99001) 'cvs print ( ) print} def'
   WRITE (Iunt,99001) '/Cshow {dup stringwidth pop -2 div 0 rmoveto show} def'
!
! we leave margins etc. in cm so it is easy to modify them if
! needed by editing the output file
   WRITE (Iunt,99001) 'gsave'
   IF ( ltit>0 ) THEN
      WRITE (Iunt,*) '/Helvetica findfont ' , fnstit , ' cm scalefont setfont '
      WRITE (Iunt,*) xtit , ' cm ' , ytit , ' cm moveto '
      WRITE (Iunt,'(3A)') '(' , Title(1:ltit) , ') Cshow'
   ENDIF
   WRITE (Iunt,*) lrmrgn , ' cm ' , botmrgn , ' cm translate'
   WRITE (Iunt,*) siz , ' cm ' , m , ' div dup scale '
!-------
! draw a frame around the matrix
   WRITE (Iunt,*) frlw , ' setlinewidth'
   WRITE (Iunt,99001) 'newpath'
   WRITE (Iunt,99002) 0 , 0 , ' moveto'
   IF ( square ) THEN
      WRITE (Iunt,99002) m , 0 , ' lineto'
      WRITE (Iunt,99002) m , m , ' lineto'
      WRITE (Iunt,99002) 0 , m , ' lineto'
   ELSE
      WRITE (Iunt,99002) nc , 0 , ' lineto'
      WRITE (Iunt,99002) nc , nr , ' lineto'
      WRITE (Iunt,99002) 0 , nr , ' lineto'
   ENDIF
   WRITE (Iunt,99001) 'closepath stroke'
!
!     drawing the separation lines
!
   WRITE (Iunt,*) ' 0.2 setlinewidth'
   DO kol = 1 , Nlines
      isep = Lines(kol)
!
!     horizontal lines
!
      yy = real(Nrow-isep) + haf
      xx = real(Ncol+1)
      WRITE (Iunt,99004) zero , yy , ' moveto '
      WRITE (Iunt,99004) xx , yy , ' lineto stroke '
!
! vertical lines
!
      xx = real(isep) + haf
      yy = real(Nrow+1)
      WRITE (Iunt,99004) xx , zero , ' moveto '
      WRITE (Iunt,99004) xx , yy , ' lineto stroke '
   ENDDO
!
!----------- plotting loop ---------------------------------------------
!
   WRITE (Iunt,99001) '1 1 translate'
   WRITE (Iunt,99001) '0.8 setlinewidth'
   WRITE (Iunt,99001) '/p {moveto 0 -.40 rmoveto '
   WRITE (Iunt,99001) '           0  .80 rlineto stroke} def'
!
   DO ii = 1 , n
      istart = Ia(ii)
      ilast = Ia(ii+1) - 1
      IF ( Mode==1 ) THEN
         DO k = istart , ilast
            WRITE (Iunt,99002) ii - 1 , Nrow - Ja(k) , ' p'
         ENDDO
      ELSE
         DO k = istart , ilast
            WRITE (Iunt,99002) Ja(k) - 1 , Nrow - ii , ' p'
         ENDDO
! add diagonal element if MSR mode.
         IF ( Mode==2 ) WRITE (Iunt,99002) ii - 1 , Nrow - ii , ' p'
!
      ENDIF
   ENDDO
!-----------------------------------------------------------------------
   WRITE (Iunt,99001) 'showpage'
   RETURN
!
99001 FORMAT (A)
99002 FORMAT (2(I6,1x),A)
99004 FORMAT (2(F9.2,1x),A)
!-----------------------------------------------------------------------
END SUBROUTINE pspltm
!*==lenstr.f90 processed by SPAG 8.04RA 15:47 21 Feb 2024
!!SPAG Open source Personal, Educational or Academic User Clemson University  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
!
FUNCTION lenstr(S)
   IMPLICIT NONE
!
! Function and Dummy argument declarations rewritten by SPAG
!
   INTEGER :: lenstr
   CHARACTER(*) , INTENT(IN) :: S
!
! Local variable declarations rewritten by SPAG
!
   INTEGER :: n
!
! End of declarations rewritten by SPAG
!
!-----------------------------------------------------------------------
! return length of the string S
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
   n = len(S)
   DO WHILE ( S(n:n)==' ' )
      n = n - 1
      IF ( n<=0 ) THEN
         CALL spag_block_1
         RETURN
      ENDIF
   ENDDO
   CALL spag_block_1
CONTAINS
   SUBROUTINE spag_block_1
      lenstr = n
   END SUBROUTINE spag_block_1
!
!--------end-of-pspltm--------------------------------------------------
!-----------------------------------------------------------------------
END FUNCTION lenstr
!*==pltmt.f90 processed by SPAG 8.04RA 15:47 21 Feb 2024
!!SPAG Open source Personal, Educational or Academic User Clemson University  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
!-----------------------------------------------------------------------
SUBROUTINE pltmt(Nrow,Ncol,Mode,Ja,Ia,Title,Key,Type,Job,Iounit)
   IMPLICIT NONE
!
! Dummy argument declarations rewritten by SPAG
!
   INTEGER , INTENT(IN) :: Nrow
   INTEGER , INTENT(IN) :: Ncol
   INTEGER , INTENT(IN) :: Mode
   INTEGER , INTENT(IN) , DIMENSION(*) :: Ja
   INTEGER , INTENT(IN) , DIMENSION(*) :: Ia
   CHARACTER(72) , INTENT(IN) :: Title
   CHARACTER(8) , INTENT(IN) :: Key
   CHARACTER(3) , INTENT(IN) :: Type
   INTEGER , INTENT(IN) :: Job
   INTEGER , INTENT(IN) :: Iounit
!
! Local variable declarations rewritten by SPAG
!
   REAL :: hscale , ptsize , tiny , vscale , x , xht , xnrow , xshift , xwid , y , yshift
   INTEGER :: ii , ilast , ips , istart , k , maxdim , n , nnz
!
! End of declarations rewritten by SPAG
!
!-----------------------------------------------------------------------
! this subroutine creates a 'pic' file for plotting the pattern of
! a sparse matrix stored in general sparse format. it is not intended
! to be a means of plotting large matrices (it is very inefficient).
! It is however useful for small matrices and can be used for example
! for inserting matrix plots in a text. The size of the plot can be
! 7in x 7in or 5 in x 5in .. There is also an option for writing a
! 3-line header in troff (see description of parameter job).
! Author: Youcef Saad - Date: Sept., 1989
! See SPARSKIT/UNSUPP/ for a version of this to produce a post-script
! file.
!-----------------------------------------------------------------------
! nrow   = number of rows in matrix
!
! ncol	 = number of columns in matrix
!
! mode   = integer indicating whether the matrix is stored
!          row-wise (mode = 0) or column-wise (mode=1)
!
! ja     = column indices of nonzero elements when matrix is
!	   stored rowise. Row indices if stores column-wise.
! ia     = integer array of containing the pointers to the
!	   beginning of the columns in arrays a, ja.
!
! title  = character*71 = title of matrix test ( character a*71 ).
! key    = character*8  = key of matrix
! type   = character*3  = type of matrix.
!
! job    = this integer parameter allows to set a few minor
!          options. First it tells pltmt whether or not to
!          reduce the plot. The standard size of 7in is then
!          replaced by a 5in plot. It also tells pltmt whether or
!          not to append to the pic file a few 'troff' lines that
!          produce a centered caption includingg the title, key and
!          types as well as the size and number of nonzero elements.
!          job =  0 : do not reduce and do not make caption.
!          job =  1 : reduce and do not make caption.
!          job = 10 : do not reduce and make caption
!          job = 11 : reduce and make caption.
!          (i.e. trailing digit for reduction, leading digit for caption)
!
! iounit = logical unit number where to write the matrix into.
!
!-----------------------------------------------------------------------
! example of usage .
!-----------------
! In the fortran code:
!  a) read a Harwell/Boeing matrix
!          call readmt (.....)
!	   iout = 13
!  b) generate pic file:
!          call  pltmt (nrow,ncol,mode,ja,ia,title,key,type,iout)
!	   stop
! ---------
! Then in a unix environment plot the matrix by the command
!
!	pic FOR013.DAT | troff -me | lpr -Ppsx
!
!-----------------------------------------------------------------------
! notes: 1) Plots square as well as rectangular matrices.
!            (however not as much tested with rectangular matrices.)
!	  2) the dot-size is adapted according to the size of the
!            matrix.
!	  3) This is not meant at all as a way of plotting large
!            matrices. The pic file generaled will have one line for
!            each nonzero element. It is  only meant for use in
!	     such things as document poreparations etc..
!         4) The caption written will print the 71 character long
!            title. This may not be centered correctly if the
!            title has trailing blanks (a problem with Troff).
!            if you want the title centered then you can center
!            the string in title before calling pltmt.
!
!-----------------------------------------------------------------------
!-------
   n = Ncol
   IF ( Mode==0 ) n = Nrow
   nnz = Ia(n+1) - Ia(1)
   maxdim = max0(Nrow,Ncol)
   xnrow = real(Nrow)
   ptsize = 0.08
   hscale = (7.0-2.0*ptsize)/real(maxdim-1)
   vscale = hscale
   xwid = ptsize + real(Ncol-1)*hscale + ptsize
   xht = ptsize + real(Nrow-1)*vscale + ptsize
   xshift = (7.0-xwid)/2.0
   yshift = (7.0-xht)/2.0
!------
   IF ( mod(Job,10)==1 ) THEN
      WRITE (Iounit,99001)
99001 FORMAT ('.PS 5in',/,'.po 1.8i')
   ELSE
      WRITE (Iounit,99002)
99002 FORMAT ('.PS',/,'.po 0.7i')
   ENDIF
   WRITE (Iounit,99003)
99003 FORMAT ('box invisible wid 7.0 ht 7.0 with .sw at (0.0,0.0) ')
   WRITE (Iounit,99004) xwid , xht , xshift , yshift
99004 FORMAT ('box wid ',f5.2,' ht ',f5.2,' with .sw at (',f5.2,',',f5.2,')')
!
!     shift points slightly to account for size of dot , etc..
!
   tiny = 0.03
   IF ( mod(Job,10)==1 ) tiny = 0.05
   xshift = xshift + ptsize - tiny
   yshift = yshift + ptsize + tiny
!
!-----------------------------------------------------------------------
!
   ips = 8
   IF ( maxdim<=500 ) ips = 10
   IF ( maxdim<=300 ) ips = 12
   IF ( maxdim<=100 ) ips = 16
   IF ( maxdim<50 ) ips = 24
   WRITE (Iounit,99005) ips
99005 FORMAT ('.ps ',i2)
!
!-----------plottingloop ---------------------------------------------
!
   DO ii = 1 , n
      istart = Ia(ii)
      ilast = Ia(ii+1) - 1
      IF ( Mode/=0 ) THEN
         x = real(ii-1)
         DO k = istart , ilast
            y = xnrow - real(Ja(k))
            WRITE (Iounit,99006) xshift + x*hscale , yshift + y*vscale
         ENDDO
      ELSE
         y = xnrow - real(ii)
         DO k = istart , ilast
            x = real(Ja(k)-1)
            WRITE (Iounit,99006) xshift + x*hscale , yshift + y*vscale
         ENDDO
      ENDIF
   ENDDO
   WRITE (Iounit,99007)
99007 FORMAT ('.PE')
!     quit if caption not desired.
   IF ( (Job/10)/=1 ) RETURN
!
   WRITE (Iounit,99008) Key , Type , Title
99008 FORMAT ('.sp 4'/'.ll 7i'/'.ps 12'/'.po 0.7i'/'.ce 3'/,'Matrix:  ',a8,',  Type:  ',a3,/,a72)
   WRITE (Iounit,99009) Nrow , Ncol , nnz
99009 FORMAT ('Dimension: ',i4,' x ',i4,',  Nonzero elements: ',i5)
!-----------------------------------------------------------------------
99006 FORMAT ('"." at ',f6.3,',',f6.3,' ljust  ')
!----------------end-of-pltmt ------------------------------------------
!-----------------------------------------------------------------------
END SUBROUTINE pltmt
!*==smms.f90 processed by SPAG 8.04RA 15:47 21 Feb 2024
!!SPAG Open source Personal, Educational or Academic User Clemson University  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
!-----------------------------------------------------------------------
SUBROUTINE smms(N,First,Last,Mode,A,Ja,Ia,Iout)
   USE ISO_FORTRAN_ENV                 
   IMPLICIT NONE
!
! Dummy argument declarations rewritten by SPAG
!
   INTEGER , INTENT(IN) :: N
   INTEGER , INTENT(IN) :: First
   INTEGER , INTENT(IN) :: Last
   INTEGER , INTENT(IN) :: Mode
   REAL(REAL64) , INTENT(IN) , DIMENSION(*) :: A
   INTEGER , INTENT(IN) , DIMENSION(*) :: Ja
   INTEGER , INTENT(IN) , DIMENSION(*) :: Ia
   INTEGER , INTENT(IN) :: Iout
!
! Local variable declarations rewritten by SPAG
!
   LOGICAL :: csc , msr
   INTEGER :: i , k , k1 , k2
!
! End of declarations rewritten by SPAG
!
!-----------------------------------------------------------------------
! writes a matrix in Coordinate (SMMS) format --
!-----------------------------------------------------------------------
! on entry:
!---------
! n     = integer = size of matrix -- number of rows (columns if matrix
!         is stored columnwise)
! first  = first row (column) to be output. This routine will output
!          rows (colums) first to last.
! last   = last row (column) to be output.
! mode   = integer giving some information about the storage of the
!          matrix. A 3-digit decimal number. 'htu'
!         * u = 0 means that matrix is stored row-wise
!         * u = 1 means that matrix is stored column-wise
!         * t = 0 indicates that the matrix is stored in CSR format
!         * t = 1 indicates that the matrix is stored in MSR format.
!         * h = ... to be added.
! a,
! ja,
! ia    =  matrix in CSR or MSR format (see mode)
! iout  = output unit number.
!
! on return:
!----------
! the output file iout will have written in it the matrix in smms
! (coordinate format)
!
!-----------------------------------------------------------------------
!
! determine mode ( msr or csr )
!
   msr = .FALSE.
   csc = .FALSE.
   IF ( mod(Mode,10)==1 ) csc = .TRUE.
   IF ( (Mode/10)==1 ) msr = .TRUE.
 
   WRITE (Iout,*) N
   DO i = First , Last
      k1 = Ia(i)
      k2 = Ia(i+1) - 1
!     write (iout,*) ' row ', i
      IF ( msr ) WRITE (Iout,'(2i6,e22.14)') i , i , A(i)
      DO k = k1 , k2
         IF ( csc ) THEN
            WRITE (Iout,'(2i6,e22.14)') Ja(k) , i , A(k)
         ELSE
            WRITE (Iout,'(2i6,e22.14)') i , Ja(k) , A(k)
         ENDIF
      ENDDO
   ENDDO
!----end-of-smms--------------------------------------------------------
!-----------------------------------------------------------------------
END SUBROUTINE smms
!*==readsm.f90 processed by SPAG 8.04RA 15:47 21 Feb 2024
!!SPAG Open source Personal, Educational or Academic User Clemson University  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
!-----------------------------------------------------------------------
SUBROUTINE readsm(Nmax,Nzmax,N,Nnz,Ia,Ja,A,Iout,Ierr)
   USE ISO_FORTRAN_ENV                 
   IMPLICIT NONE
!
! Dummy argument declarations rewritten by SPAG
!
   INTEGER , INTENT(IN) :: Nmax
   INTEGER , INTENT(IN) :: Nzmax
   INTEGER , INTENT(INOUT) :: N
   INTEGER , INTENT(INOUT) :: Nnz
   INTEGER , INTENT(OUT) , DIMENSION(Nmax+1) :: Ia
   INTEGER , INTENT(OUT) , DIMENSION(Nzmax) :: Ja
   REAL(REAL64) , INTENT(OUT) , DIMENSION(Nzmax) :: A
   INTEGER , INTENT(IN) :: Iout
   INTEGER , INTENT(OUT) :: Ierr
!
! Local variable declarations rewritten by SPAG
!
   INTEGER :: i , j , k , row
   REAL(REAL64) :: x
!
! End of declarations rewritten by SPAG
!
   INTEGER :: spag_nextblock_1
   spag_nextblock_1 = 1
   SPAG_DispatchLoop_1: DO
      SELECT CASE (spag_nextblock_1)
      CASE (1)
!-----------------------------------------------------------------------
!     read a matrix in coordinate format as is used in the SMMS
!     package (F. Alvarado), i.e. the row is in ascending order.
!     Outputs the matrix in CSR format.
!-----------------------------------------------------------------------
! coded by Kesheng Wu on Oct 21, 1991 with the supervision of Y. Saad
!-----------------------------------------------------------------------
! on entry:
!---------
! nmax  = the maximum size of array
! nzmax = the maximum number of nonzeros
! iout  = the I/O unit that has the data file
!
! on return:
!----------
! n     = integer = size of matrix
! nnz   = number of non-zero entries in the matrix
! a,
! ja,
! ia    = matrix in CSR format
! ierr  = error code,
!         0 -- subroutine end with intended job done
!         1 -- error in I/O unit iout
!         2 -- end-of-file reached while reading n, i.e. a empty data file
!         3 -- n non-positive or too large
!         4 -- nnz is zero or larger than nzmax
!         5 -- data file is not orgnized in the order of ascending
!              row indices
!
! in case of errors:
!   n will be set to zero (0). In case the data file has more than nzmax
!   number of entries, the first nzmax entries will be read, and are not
!   cleared on return. The total number of entry is determined.
!   Ierr is set.
!-----------------------------------------------------------------------
!
         REWIND (Iout)
         Nnz = 0
         Ia(1) = 1
         row = 1
!
         READ (Iout,*,ERR=40,END=60) N
         IF ( (N<=0) .OR. (N>Nmax) ) THEN
!
!     problem with n
!
            Ierr = 3
            GOTO 80
         ELSE
            DO
!
               Nnz = Nnz + 1
               READ (Iout,*,ERR=40,END=20) i , j , x
 
!     set the pointers when needed
               IF ( i>row ) THEN
                  DO k = row + 1 , i
                     Ia(k) = Nnz
                  ENDDO
                  row = i
               ELSEIF ( i<row ) THEN
!
!     data entries not ordered
!
                  Ierr = 5
                  GOTO 80
               ENDIF
 
               Ja(Nnz) = j
               A(Nnz) = x
 
               IF ( Nnz>=Nzmax ) THEN
                  spag_nextblock_1 = 2
                  CYCLE SPAG_DispatchLoop_1
               ENDIF
            ENDDO
         ENDIF
 
!     normal return -- end of file reached
 20      Ia(row+1) = Nnz
         Nnz = Nnz - 1
         IF ( Nnz==0 ) THEN
            spag_nextblock_1 = 2
            CYCLE SPAG_DispatchLoop_1
         ENDIF
!
!     everything seems to be OK.
!
         Ierr = 0
         RETURN
!
!     error handling code
!
!     error in reading data entries
!
 40      Ierr = 1
         GOTO 80
!
!     empty file
!
 60      Ierr = 2
         GOTO 80
      CASE (2)
!
!     problem with nnz
!
         Ierr = 4
!
!     try to determine the real number of entries, in case needed
!
         IF ( Nnz>=Nzmax ) THEN
            DO
               READ (Iout,*,ERR=80,END=80) i , j , x
               Nnz = Nnz + 1
            ENDDO
         ENDIF
 80      N = 0
         EXIT SPAG_DispatchLoop_1
      END SELECT
   ENDDO SPAG_DispatchLoop_1
!----end-of-readsm------------------------------------------------------
!-----------------------------------------------------------------------
END SUBROUTINE readsm
!*==readsk.f90 processed by SPAG 8.04RA 15:47 21 Feb 2024
!!SPAG Open source Personal, Educational or Academic User Clemson University  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
!-----------------------------------------------------------------------
SUBROUTINE readsk(Nmax,Nzmax,N,Nnz,A,Ja,Ia,Iounit,Ierr)
   USE ISO_FORTRAN_ENV                 
   IMPLICIT NONE
!
! Dummy argument declarations rewritten by SPAG
!
   INTEGER , INTENT(IN) :: Nmax
   INTEGER , INTENT(IN) :: Nzmax
   INTEGER , INTENT(INOUT) :: N
   INTEGER , INTENT(INOUT) :: Nnz
   REAL(REAL64) , INTENT(OUT) , DIMENSION(Nzmax) :: A
   INTEGER , INTENT(OUT) , DIMENSION(Nzmax) :: Ja
   INTEGER , INTENT(INOUT) , DIMENSION(Nmax+1) :: Ia
   INTEGER , INTENT(IN) :: Iounit
   INTEGER , INTENT(OUT) :: Ierr
!
! Local variable declarations rewritten by SPAG
!
   INTEGER :: i
!
! End of declarations rewritten by SPAG
!
!-----------------------------------------------------------------------
! Reads matrix in Compressed Saprse Row format. The data is supposed to
! appear in the following order -- n, ia, ja, a
! Only square matrices accepted. Format has following features
! (1) each number is separated by at least one space (or end-of-line),
! (2) each array starts with a new line.
!-----------------------------------------------------------------------
! coded by Kesheng Wu on Oct 21, 1991 with supervision of Y. Saad
!-----------------------------------------------------------------------
! on entry:
!---------
! nmax 	 = max column dimension  allowed for matrix.
! nzmax	 = max number of nonzeros elements allowed. the arrays a,
!          and ja should be of length equal to nnz (see below).
! iounit = logical unit number where to read the matrix from.
!
! on return:
!----------
! ia,
! ja,
! a      = matrx in CSR format
! n      = number of rows(columns) in matrix
! nnz	 = number of nonzero elements in A. This info is returned
!          even if there is not enough space in a, ja, ia, in order
!          to determine the minimum storage needed.
! ierr   = error code,
!          0 : OK;
!          1 : error when try to read the specified I/O unit.
!          2 : end-of-file reached during reading of data file.
!          3 : array size in data file is negtive or larger than nmax;
!          4 : nunmer of nonzeros in data file is negtive or larger than nzmax
! in case of errors:
!---------
!     n is set to 0 (zero), at the same time ierr is set.
!-----------------------------------------------------------------------
!
!     read the size of the matrix
!
   REWIND (Iounit)
   READ (Iounit,*,ERR=100,END=200) N
   IF ( (N<=0) .OR. (N>Nmax) ) THEN
!
!     n non-positive or too large
      Ierr = 3
      N = 0
!
!     the real return statement
!
      N = 0
      RETURN
   ELSE
!
!     read the pointer array ia(*)
!
      READ (Iounit,*,ERR=100,END=200) (Ia(i),i=1,N+1)
!
!     Number of None-Zeros
!
      Nnz = Ia(N+1) - 1
      IF ( (Nnz<=0) .OR. (Nnz>Nzmax) ) THEN
!
!     NNZ non-positive or too large
         Ierr = 4
         N = 0
         RETURN
      ELSE
!
!     read the column indices array
!
         READ (Iounit,*,ERR=100,END=200) (Ja(i),i=1,Nnz)
!
!     read the matrix elements
!
         READ (Iounit,*,ERR=100,END=200) (A(i),i=1,Nnz)
!
!     normal return
!
         Ierr = 0
         RETURN
      ENDIF
   ENDIF
!
!     error handling code
!
!     error in reading I/O unit
 100  Ierr = 1
   N = 0
   RETURN
!
!     EOF reached in reading
 200  Ierr = 2
   N = 0
!---------end of readsk ------------------------------------------------
!-----------------------------------------------------------------------
END SUBROUTINE readsk
!*==skit.f90 processed by SPAG 8.04RA 15:47 21 Feb 2024
!!SPAG Open source Personal, Educational or Academic User Clemson University  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
!-----------------------------------------------------------------------
SUBROUTINE skit(N,A,Ja,Ia,Ifmt,Iounit,Ierr)
   USE ISO_FORTRAN_ENV                 
   IMPLICIT NONE
!
! Dummy argument declarations rewritten by SPAG
!
   INTEGER , INTENT(IN) :: N
   REAL(REAL64) , INTENT(IN) , DIMENSION(*) :: A
   INTEGER , INTENT(IN) , DIMENSION(*) :: Ja
   INTEGER , INTENT(IN) , DIMENSION(*) :: Ia
   INTEGER , INTENT(INOUT) :: Ifmt
   INTEGER , INTENT(IN) :: Iounit
   INTEGER , INTENT(OUT) :: Ierr
!
! Local variable declarations rewritten by SPAG
!
   INTEGER :: i , ihead , len , nnz , nperli
   CHARACTER(16) :: indfmt , ptrfmt
   CHARACTER(28) :: ix
   CHARACTER(20) :: valfmt
!
! End of declarations rewritten by SPAG
!
!-----------------------------------------------------------------------
!     Writes a matrix in Compressed Sparse Row format to an I/O unit.
!     It tryes to pack as many number as possible into lines of less than
!     80 characters. Space is inserted in between numbers for separation
!     to avoid carrying a header in the data file. This can be viewed
!     as a simplified Harwell-Boeing format.
!-----------------------------------------------------------------------
! Modified from subroutine prtmt written by Y. Saad
!-----------------------------------------------------------------------
! on entry:
!---------
! n      = number of rows(columns) in matrix
! a      = real*8 array containing the values of the matrix stored
!          columnwise
! ja     = integer array of the same length as a containing the column
!          indices of the corresponding matrix elements of array a.
! ia     = integer array of containing the pointers to the beginning of
!          the row in arrays a and ja.
! ifmt   = integer specifying the format chosen for the real values
!          to be output (i.e., for a, and for rhs-guess-sol if
!          applicable). The meaning of ifmt is as follows.
!          * if (ifmt .lt. 100) then the D descriptor is used,
!          format Dd.m, in which the length (m) of the mantissa is
!          precisely the integer ifmt (and d = ifmt+6)
!          * if (ifmt .gt. 100) then prtmt will use the
!          F- descriptor (format Fd.m) in which the length of the
!          mantissa (m) is the integer mod(ifmt,100) and the length
!          of the integer part is k=ifmt/100 (and d = k+m+2)
!          Thus  ifmt= 4   means  D10.4  +.xxxxD+ee    while
!          ifmt=104  means  F7.4   +x.xxxx
!          ifmt=205  means  F9.5   +xx.xxxxx
!          Note: formats for ja, and ia are internally computed.
!
! iounit = logical unit number where to write the matrix into.
!
! on return:
!----------
! ierr   = error code, 0 for normal 1 for error in writing to iounit.
!
! on error:
!--------
!     If error is encontacted when writing the matrix, the whole matrix
!     is written to the standard output.
!     ierr is set to 1.
!-----------------------------------------------------------------------
!--------------
!     compute pointer format
!--------------
   nnz = Ia(N+1)
   len = int(alog10(0.1+real(nnz))) + 2
   nnz = nnz - 1
   nperli = 80/len
 
   PRINT * , ' skit entries:' , N , nnz , len , nperli
 
   IF ( len>9 ) THEN
      ix = '(1h(,i2,1HI,i2,1h) )'
   ELSE
      ix = '(1h(,i2,1HI,i1,1h) )'
   ENDIF
   WRITE (ptrfmt,ix) nperli , len
 
!----------------------------
!     compute ROW index format
!----------------------------
   len = int(alog10(0.1+real(N))) + 2
   nperli = min0(80/len,nnz)
   WRITE (indfmt,'(1h(,i2,1HI,i1,1h) )') nperli , len
!---------------------------
!     compute value format
!---------------------------
   IF ( Ifmt>=100 ) THEN
      ihead = Ifmt/100
      Ifmt = Ifmt - 100*ihead
      len = ihead + Ifmt + 3
      nperli = 80/len
!
      IF ( len<=9 ) THEN
         ix = '(1h(,i2,1hF,i1,1h.,i1,1h) )'
      ELSEIF ( Ifmt<=9 ) THEN
         ix = '(1h(,i2,1hF,i2,1h.,i1,1h) )'
      ELSE
         ix = '(1h(,i2,1hF,i2,1h.,i2,1h) )'
      ENDIF
!
      WRITE (valfmt,ix) nperli , len , Ifmt
   ELSE
      len = Ifmt + 7
      nperli = 80/len
!     try to minimize the blanks in the format strings.
      IF ( nperli<=9 ) THEN
         IF ( len<=9 ) THEN
            ix = '(1h(,i1,1hD,i1,1h.,i1,1h) )'
         ELSEIF ( Ifmt<=9 ) THEN
            ix = '(1h(,i1,1hD,i2,1h.,i1,1h) )'
         ELSE
            ix = '(1h(,i1,1hD,i2,1h.,i2,1h) )'
         ENDIF
      ELSEIF ( len<=9 ) THEN
         ix = '(1h(,i2,1hD,i1,1h.,i1,1h) )'
      ELSEIF ( Ifmt<=9 ) THEN
         ix = '(1h(,i2,1hD,i2,1h.,i1,1h) )'
      ELSE
         ix = '(1h(,i2,1hD,i2,1h.,i2,1h) )'
      ENDIF
!-----------
      WRITE (valfmt,ix) nperli , len , Ifmt
   ENDIF
!
!     output the data
!
   WRITE (Iounit,*) N
   WRITE (Iounit,ptrfmt,ERR=100) (Ia(i),i=1,N+1)
   WRITE (Iounit,indfmt,ERR=100) (Ja(i),i=1,nnz)
   WRITE (Iounit,valfmt,ERR=100) (A(i),i=1,nnz)
!
!     done, if no trouble is encounted in writing data
!
   Ierr = 0
   RETURN
!
!     if can't write the data to the I/O unit specified, should be able to
!     write everything to standard output (unit 6)
!
 100  WRITE (0,*) 'Error, Can''t write data to sepcified unit' , Iounit
   WRITE (0,*) 'Write the matrix into standard output instead!'
   Ierr = 1
   WRITE (6,*) N
   WRITE (6,ptrfmt) (Ia(i),i=1,N+1)
   WRITE (6,indfmt) (Ja(i),i=1,nnz)
   WRITE (6,valfmt) (A(i),i=1,nnz)
!----------end of skit -------------------------------------------------
!-----------------------------------------------------------------------
END SUBROUTINE skit
!*==prtunf.f90 processed by SPAG 8.04RA 15:47 21 Feb 2024
!!SPAG Open source Personal, Educational or Academic User Clemson University  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
!-----------------------------------------------------------------------
SUBROUTINE prtunf(N,A,Ja,Ia,Iout,Ierr)
   USE ISO_FORTRAN_ENV                 
   IMPLICIT NONE
!
! Dummy argument declarations rewritten by SPAG
!
   INTEGER , INTENT(IN) :: N
   REAL(REAL64) , INTENT(IN) , DIMENSION(*) :: A
   INTEGER , INTENT(IN) , DIMENSION(*) :: Ja
   INTEGER , INTENT(IN) , DIMENSION(*) :: Ia
   INTEGER , INTENT(IN) :: Iout
   INTEGER , INTENT(OUT) :: Ierr
!
! Local variable declarations rewritten by SPAG
!
   INTEGER :: k , nnz
!
! End of declarations rewritten by SPAG
!
!-----------------------------------------------------------------------
! This subroutine dumps the arrays used for storing sparse compressed row
! format in machine code, i.e. unformatted using standard FORTRAN term.
!-----------------------------------------------------------------------
! First coded by Kesheng Wu on Oct 21, 1991 under the instruction of
! Prof. Y. Saad
!-----------------------------------------------------------------------
! On entry:
!     n: the size of the matrix (matrix is n X n)
!    ia: integer array stores the stariting position of each row.
!    ja: integer array stores the column indices of each entry.
!     a: the non-zero entries of the matrix.
!  iout: the unit number opened for storing the matrix.
! On return:
!  ierr: a error, 0 if everything's OK, else 1 if error in writing data.
! On error:
!  set ierr to 1.
!  No redirection is made, since direct the machine code to the standard
! output may cause unpridictable consequences.
!-----------------------------------------------------------------------
   nnz = Ia(N+1) - Ia(1)
!
   WRITE (UNIT=Iout,ERR=100) N
   WRITE (UNIT=Iout,ERR=100) (Ia(k),k=1,N+1)
   IF ( nnz>0 ) THEN
      WRITE (UNIT=Iout,ERR=100) (Ja(k),k=1,nnz)
      WRITE (UNIT=Iout,ERR=100) (A(k),k=1,nnz)
   ENDIF
!
   Ierr = 0
   RETURN
!
 100  Ierr = 1
END SUBROUTINE prtunf
!*==readunf.f90 processed by SPAG 8.04RA 15:47 21 Feb 2024
!!SPAG Open source Personal, Educational or Academic User Clemson University  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
!---------end of prtunf ------------------------------------------------
!
!-----------------------------------------------------------------------
SUBROUTINE readunf(Nmax,Nzmax,N,Nnz,A,Ja,Ia,Iounit,Ierr)
   USE ISO_FORTRAN_ENV                 
   IMPLICIT NONE
!
! Dummy argument declarations rewritten by SPAG
!
   INTEGER , INTENT(IN) :: Nmax
   INTEGER , INTENT(IN) :: Nzmax
   INTEGER , INTENT(INOUT) :: N
   INTEGER , INTENT(INOUT) :: Nnz
   REAL(REAL64) , INTENT(OUT) , DIMENSION(Nzmax) :: A
   INTEGER , INTENT(OUT) , DIMENSION(Nzmax) :: Ja
   INTEGER , INTENT(INOUT) , DIMENSION(Nmax+1) :: Ia
   INTEGER , INTENT(IN) :: Iounit
   INTEGER , INTENT(OUT) :: Ierr
!
! Local variable declarations rewritten by SPAG
!
   INTEGER :: k
!
! End of declarations rewritten by SPAG
!
!-----------------------------------------------------------------------
! This subroutine reads a matix store in machine code (FORTRAN
! unformatted form). The matrix is in CSR format.
!-----------------------------------------------------------------------
! First coded by Kesheng Wu on Oct 21, 1991 under the instruction of
! Prof. Y. Saad
!-----------------------------------------------------------------------
! On entry:
!    nmax: the maximum value of matrix size.
!   nzmax: the maximum number of non-zero entries.
!  iounit: the I/O unit that opened for reading.
! On return:
!       n: the actual size of array.
!     nnz: the actual number of non-zero entries.
! ia,ja,a: the matrix in CSR format.
!    ierr: a error code, it's same as that used in reaadsk
!          0 -- OK
!          1 -- error in reading iounit
!          2 -- end-of-file reached while reading data file
!          3 -- n is non-positive or too large
!          4 -- nnz is non-positive or too large
! On error:
!     return with n set to 0 (zero). nnz is kept if it's set already,
!     in case one want to use it to determine the size of array needed
!     to hold the data.
!-----------------------------------------------------------------------
!
!
   REWIND Iounit
!
   READ (UNIT=Iounit,ERR=100,END=200) N
   IF ( (N<=0) .OR. (N>Nmax) ) THEN
      Ierr = 3
      N = 0
      RETURN
   ELSE
!
      READ (UNIT=Iounit,ERR=100,END=200) (Ia(k),k=1,N+1)
!
      Nnz = Ia(N+1) - 1
      IF ( (Nnz<=0) .OR. (Nnz>Nzmax) ) THEN
         Ierr = 4
         N = 0
         RETURN
      ELSE
!
         READ (UNIT=Iounit,ERR=100,END=200) (Ja(k),k=1,Nnz)
         READ (UNIT=Iounit,ERR=100,END=200) (A(k),k=1,Nnz)
!
!     everything seems to be OK.
!
         Ierr = 0
         RETURN
      ENDIF
   ENDIF
!
!     error handling
!
 100  Ierr = 1
   N = 0
   RETURN
 200  Ierr = 2
   N = 0
END SUBROUTINE readunf
