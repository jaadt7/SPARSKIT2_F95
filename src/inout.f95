subroutine readmt(nmax,nzmax,job,iounit,a,ja,ia,rhs,nrhs,guesol,nrow,ncol,nnz,title,key,type,ierr)

     integer, intent(InOut) :: job, nrhs
     integer, intent(in) :: iounit, nmax,nzmax
     integer, intent(Out) :: ierr,nrow,ncol,nnz
     character(len=72), intent(out) :: title
     character(len=8), intent(out) :: key
     character(len=3), intent(out) :: type
     integer, dimension(nmax + 1), intent(out) :: ia
     integer, dimension(nzmax), intent(out) :: ja
     real(kind=8), dimension(nzmax), intent(out) :: a
     real(kind=8), dimension(*), intent(out) :: rhs
     character(len=2), intent(Out) :: guesol
 
     integer :: lenrhs,n,i,nvec,len,next,iend,totcrd,ptrcrd,indcrd,valcrd,rhscrd,neltvl,nrwindx
 
     character(len=16) :: ptrfmt, indfmt
     character(len=20) :: valfmt, rhsfmt
     character(len=3) :: rhstyp

! -----------------------------------------------------------------------
! -----------------------------------------------------------------------
!  this subroutine reads  a boeing/harwell matrix. handles right hand 
!  sides in full format only (no sparse right hand sides).
!  Also the matrix must be in assembled forms.
!  Author: Youcef Saad - Date: Sept., 1989
!          updated Oct 31, 1989.
! -----------------------------------------------------------------------
!  on entry:
! ---------
!  nmax 	 =  max column dimension  allowed for matrix. The array ia should 
!          be of length at least ncol+1 (see below) if job.gt.0
!  nzmax	 = max number of nonzeros elements allowed. the arrays a, 
!           and ja should be of length equal to nnz (see below) if these
!           arrays are to be read (see job).
!           
!  job	 = integer to indicate what is to be read. (note: job is an
!           input and output parameter, it can be modified on return)
!           job = 0    read the values of ncol, nrow, nnz, title, key,
!                      type and return. matrix is not read and arrays
!                      a, ja, ia, rhs are not touched.
!           job = 1    read srtucture only, i.e., the arrays ja and ia.
!           job = 2    read matrix including values, i.e., a, ja, ia
!           job = 3    read matrix and right hand sides: a,ja,ia,rhs.
!      	      rhs may contain initial guesses and exact 
!                      solutions appended to the actual right hand sides.
!      	      this will be indicated by the output parameter
!                      guesol [see below]. 
!                      
!  nrhs   = integer. nrhs is an input as well as ouput parameter.
!           at input nrhs contains the total length of the array rhs.
!           See also ierr and nrhs in output parameters.
! 
!  iounit = logical unit number where to read the matrix from.
! 
!  on return:
! ---------- 
!  job    = on return job may be modified to the highest job it could
!           do: if job=2 on entry but no matrix values are available it
!           is reset to job=1 on return. Similarly of job=3 but no rhs 
!           is provided then it is rest to job=2 or job=1 depending on 
!           whether or not matrix values are provided.
!           Note that no error message is triggered (i.e. ierr = 0 
!           on return in these cases. It is therefore important to
!           compare the values of job on entry and return ).
! 
!  a	 = the a matrix in the a, ia, ja (column) storage format
!  ja 	 = row number of element a(i,j) in array a.
!  ia     = pointer  array. ia(i) points to the beginning of column i.
! 
!  rhs    = real array of size nrow + 1 if available (see job)
! 
!  nrhs   = integer containing the number of right-hand sides found
!           each right hand side may be accompanied with an intial guess
!           and also the exact solution.
! 
!  guesol = a 2-character string indicating whether an initial guess 
!           (1-st character) and / or the exact solution (2-nd
!           character) is provided with the right hand side.
!         if the first character of guesol is 'G' it means that an
!           an intial guess is provided for each right-hand side.
!           These are appended to the right hand-sides in the array rhs.
!         if the second character of guesol is 'X' it means that an
!           exact solution is provided for each right-hand side.
!           These are  appended to the right hand-sides 
!           and the initial guesses (if any) in the array rhs.
! 
!  nrow   = number of rows in matrix
!  ncol	 = number of columns in matrix 
!  nnz	 = number of nonzero elements in A. This info is returned
!           even if there is not enough space in a, ja, ia, in order
!           to determine the minimum storage needed. 
! 
!  title  = character*72 = title of matrix test ( character a*72). 
!  key    = character*8  = key of matrix 
!  type   = charatcer*3  = type of matrix.
!           for meaning of title, key and type refer to documentation 
!           Harwell/Boeing matrices.
! 
!  ierr   = integer used for error messages 
!          * ierr  =  0 means that  the matrix has been read normally. 
!          * ierr  =  1 means that  the array matrix could not be read 
!          because ncol+1 .gt. nmax
!          * ierr  =  2 means that  the array matrix could not be read 
!          because nnz .gt. nzmax 
!          * ierr  =  3 means that  the array matrix could not be read 
!          because both (ncol+1 .gt. nmax) and  (nnz .gt. nzmax )
!          * ierr  =  4 means that  the right hand side (s) initial 
!          guesse (s) and exact solution (s)   could  not be
!          read because they are stored in sparse format (not handled
!          by this routine ...) 
!          * ierr  =  5 means that the right-hand-sides, initial guesses
!          and exact solutions could not be read because the length of 
!          rhs as specified by the input value of nrhs is not 
!          sufficient to store them. The rest of the matrix may have
!          been read normally.
!  
!  Notes:
! -------
!  1) The file inout must be open (and possibly rewound if necessary)
!     prior to calling readmt.
!  2) Refer to the documentation on the Harwell-Boeing formats
!     for details on the format assumed by readmt.
!     We summarize the format here for convenience.
!   
!     a) all lines in inout are assumed to be 80 character long.
!     b) the file consists of a header followed by the block of the 
!        column start pointers followed by the block of the
!        row indices, followed by the block of the real values and
!        finally the numerical values of the right-hand-side if a 
!        right hand side is supplied. 
!     c) the file starts by a header which contains four lines if no
!        right hand side is supplied and five lines otherwise.
!        * first line contains the title (72 characters long) followed by
!          the 8-character identifier (name of the matrix, called key)
!         [ A72,A8 ]
!        * second line contains the number of lines for each
!          of the following data blocks (4 of them) and the total number 
!          of lines excluding the header.
!         [5i4]
!        * the third line contains a three character string identifying
!          the type of matrices as they are referenced in the Harwell
!          Boeing documentation [e.g., rua, rsa,..] and the number of
!          rows, columns, nonzero entries.
!          [A3,11X,4I14]
!        * The fourth line contains the variable fortran format
!          for the following data blocks.
!          [2A16,2A20] 
!        * The fifth line is only present if right-hand-sides are 
!          supplied. It consists of three one character-strings containing
!          the storage format for the right-hand-sides 
!          ('F'= full,'M'=sparse=same as matrix), an initial guess 
!          indicator ('G' for yes), an exact solution indicator 
!          ('X' for yes), followed by the number of right-hand-sides
!          and then the number of row indices. 
!          [A3,11X,2I14] 
!      d) The three following blocks follow the header as described 
!         above.
!      e) In case the right hand-side are in sparse formats then 
!         the fourth block uses the same storage format as for the matrix
!         to describe the NRHS right hand sides provided, with a column
!         being replaced by a right hand side.
! -----------------------------------------------------------------------

     ierr = 0
     lenrhs = nrhs
! 
     read (iounit,10) title, key, totcrd, ptrcrd, indcrd, valcrd, rhscrd, type, nrow, ncol, nnz, neltvl, ptrfmt, indfmt, valfmt, rhsfmt
! 
     if (rhscrd  >  0) read (iounit,11) rhstyp, nrhs, nrwindx
! 
!  anything else to read ?
! 
     if (job <= 0) return
!      ---- check whether matrix is readable ------ 
     n = ncol
     if (ncol > nmax) ierr = 1
     if (nnz > nzmax) ierr = ierr + 2
     if (ierr /= 0) return
!      ---- read pointer and row numbers ---------- 
     read (iounit,ptrfmt) (ia (i), i = 1, n+1)
     read (iounit,indfmt) (ja (i), i = 1, nnz)
!      --- reading values of matrix if required....
     if (job <= 1) return
!      --- and if available ----------------------- 
     if (valcrd <= 0) then
          job = 1
          return
     endif
     
     read (iounit,valfmt) (a(i), i = 1, nnz)
!      --- reading rhs if required ---------------- 
     if (job <= 2) return
!      --- and if available ----------------------- 
     if (rhscrd <= 0) then
          job = 2
          nrhs = 0
          return
     endif
!      
!      --- read right-hand-side.-------------------- 
!      
     if (rhstyp(1:1)  ==  'M') then 
          ierr = 4
          return
     endif
! 
     guesol = rhstyp(2:3) 
!      
     nvec = 1 
     if (guesol(1:1)  ==  'G' .or. guesol(1:1)  ==  'g') nvec=nvec+1
     if (guesol(2:2)  ==  'X' .or. guesol(2:2)  ==  'x') nvec=nvec+1
!      
     len = nrhs*nrow 
!      
     if (len*nvec > lenrhs) then
          ierr = 5
          return
     endif
! 
!  read right-hand-sides
! 
     next = 1
     iend = len
     read(iounit,rhsfmt) (rhs(i), i = next, iend)
! 
!  read initial guesses if available
! 
     if (guesol(1:1)  ==  'G' .or. guesol(1:1)  ==  'g') then
          next = next+len
          iend = iend+ len
          read(iounit,valfmt) (rhs(i), i = next, iend)
     endif
!      
!  read exact solutions if available
! 
     if (guesol(2:2)  ==  'X' .or. guesol(2:2)  ==  'x') then
          next = next+len
          iend = iend+ len
          read(iounit,valfmt) (rhs(i), i = next, iend)
     endif

10 format (a72, a8 / 5i14 / a3, 11x, 4i14 / 2a16, 2a20)
11 format (a3,11x,i14,i14)
!      
     return
! --------- end of readmt -----------------------------------------------
! ----------------------------------------------------------------------- 
end subroutine readmt

subroutine prtmt(nrow,ncol,a,ja,ia,rhs,guesol,title,key,type,ifmt,job,iounit)
 
   
     integer, intent(In) :: job,iounit,nrow,ncol
     character(len=72), intent(In) :: title
     character(len=8), intent(In) :: key
     character(len=3), intent(In) :: type
     character(len=2), intent(In) :: guesol
     integer, dimension(*), intent(In) :: ja,ia
     real(kind=8), dimension(*), intent(In) :: a,rhs
 
     integer :: ihead, i, next, iend, totcrd,ptrcrd,indcrd,valcrd,rhscrd, &
               nnz,nrhs,len,nperli,nrwindx,ifmt

     character(len=16) :: ptrfmt
     character(len=16) :: indfmt
     character(len=20) :: valfmt
 
     character(len=3) :: rhstyp
     character(len=28) :: ix

! -----------------------------------------------------------------------
!  writes a matrix in Harwell-Boeing format into a file.
!  assumes that the matrix is stored in COMPRESSED SPARSE COLUMN FORMAT.
!  some limited functionality for right hand sides. 
!  Author: Youcef Saad - Date: Sept., 1989 - updated Oct. 31, 1989 to
!  cope with new format. 
! -----------------------------------------------------------------------
!  on entry:
! ---------
!  nrow   = number of rows in matrix
!  ncol	 = number of columns in matrix 
!  a	 = real*8 array containing the values of the matrix stored 
!           columnwise
!  ja 	 = integer array of the same length as a containing the column
!           indices of the corresponding matrix elements of array a.
!  ia     = integer array of containing the pointers to the beginning of 
!         the row in arrays a and ja.
!  rhs    = real array  containing the right-hand-side (s) and optionally
!           the associated initial guesses and/or exact solutions
!           in this order. See also guesol for details. the vector rhs will
!           be used only if job .gt. 2 (see below). Only full storage for
!           the right hand sides is supported. 
! 
!  guesol = a 2-character string indicating whether an initial guess 
!           (1-st character) and / or the exact solution (2-nd)
!           character) is provided with the right hand side.
!         if the first character of guesol is 'G' it means that an
!           an intial guess is provided for each right-hand sides. 
!           These are assumed to be appended to the right hand-sides in 
!           the array rhs.
!         if the second character of guesol is 'X' it means that an
!           exact solution is provided for each right-hand side.
!           These are assumed to be appended to the right hand-sides 
!           and the initial guesses (if any) in the array rhs.
! 
!  title  = character*72 = title of matrix test ( character a*72 ).
!  key    = character*8  = key of matrix 
!  type   = charatcer*3  = type of matrix.
! 
!  ifmt	 = integer specifying the format chosen for the real values
!         to be output (i.e., for a, and for rhs-guess-sol if 
!           applicable). The meaning of ifmt is as follows.
!        * if (ifmt .lt. 100) then the D descriptor is used,
!            format Dd.m, in which the length (m) of the mantissa is 
!            precisely the integer ifmt (and d = ifmt+6)
!        * if (ifmt .gt. 100) then prtmt will use the 
!            F- descriptor (format Fd.m) in which the length of the 
!            mantissa (m) is the integer mod(ifmt,100) and the length 
!            of the integer part is k=ifmt/100 (and d = k+m+2)
!          Thus  ifmt= 4   means  D10.4  +.xxxxD+ee    while
!                ifmt=104  means  F7.4   +x.xxxx
!                ifmt=205  means  F9.5   +xx.xxxxx
!          Note: formats for ja, and ia are internally computed.
! 
!  job	 = integer to indicate whether matrix values and
!         a right-hand-side is available to be written
!           job = 1   write srtucture only, i.e., the arrays ja and ia.
!           job = 2   write matrix including values, i.e., a, ja, ia
!           job = 3   write matrix and one right hand side: a,ja,ia,rhs.
!         job = nrhs+2 write matrix and nrhs successive right hand sides
!         Note that there cannot be any right-hand-side if the matrix
!         has no values. Also the initial guess and exact solutions when 
!           provided are for each right hand side. For example if nrhs=2 
!           and guesol='GX' there are 6 vectors to write.
!           
! 
!  iounit = logical unit number where to write the matrix into.
! 
!  on return:
! ---------- 
!  the matrix a, ja, ia will be written in output unit iounit
!  in the Harwell-Boeing format. None of the inputs is modofied.
!   
!  Notes: 1) This code attempts to pack as many elements as possible per
!         80-character line. 
!         2) this code attempts to avoid as much as possible to put
!         blanks in the formats that are written in the 4-line header
!          (This is done for purely esthetical reasons since blanks
!         are ignored in format descriptors.)
!         3) sparse formats for right hand sides and guesses are not
!         supported.
! -----------------------------------------------------------------------

! --------------
!      compute pointer format
! --------------
     nnz = ia(ncol+1) -1
     
     if (nnz == 0) return
 
     len = int ( alog10(0.1+real(nnz+1))) + 1
     nperli = 80/len
     ptrcrd = ncol/nperli + 1
     if (len > 9) then
!          assign 101 to ix
         ix = '(1h(,i2,1HI,i2,1h) )'
     else
         ix = '(1h(,i2,1HI,i1,1h) )'
     endif
     write (ptrfmt,ix) nperli,len
! ----------------------------
!  compute ROW index format
! ----------------------------
     len = int ( alog10(0.1+real(nrow) )) + 1
     nperli = min0(80/len,nnz)
     indcrd = (nnz-1)/nperli+1
     write (indfmt,'(1h(,i2,1HI,i1,1h))') nperli,len
! ---------------
!  compute values and rhs format (using the same for both)
! --------------- 
     valcrd = 0
     rhscrd = 0
!  quit this part if no values provided.
     block 
          if (job <= 1) exit !leaves block      
          if (ifmt >= 100) then
               ihead = ifmt/100
               ifmt = ifmt-100*ihead
               len = ihead+ifmt+2
               nperli = 80/len
!      
               if (len <= 9 ) then
                    ix = '(1h(,i2,1hF,i1,1h.,i1,1h) )'
               elseif (ifmt <= 9) then
                    ix = '(1h(,i2,1hF,i2,1h.,i1,1h) )'
               else 
                    ix = '(1h(,i2,1hF,i2,1h.,i2,1h) )'
               endif
!      
               write(valfmt,ix) nperli,len,ifmt
          else
               len = ifmt + 6
               nperli = 80/len
!      try to minimize the blanks in the format strings.
               if (nperli <= 9) then
                    if (len <= 9 ) then
                         ix = '(1h(,i1,1hD,i1,1h.,i1,1h) )'
                    elseif (ifmt <= 9) then
                         ix = '(1h(,i1,1hD,i2,1h.,i1,1h) )'
                    else 
                         ix = '(1h(,i1,1hD,i2,1h.,i2,1h) )'
                    endif
               else 
                    if (len <= 9 ) then
                         ix = '(1h(,i2,1hD,i1,1h.,i1,1h) )'
                    elseif (ifmt <= 9) then
                         ix = '(1h(,i2,1hD,i2,1h.,i1,1h) )'
                    else 
                         ix = '(1h(,i2,1hD,i2,1h.,i2,1h) )'
                    endif
               endif
! -----------
               write(valfmt,ix) nperli,len,ifmt
     
          endif 
          valcrd = (nnz-1)/nperli+1
          nrhs = job -2
 
          if (nrhs >= 1) then
               i = (nrhs*nrow-1)/nperli+1
               rhscrd = i
               if (guesol(1:1)  ==  'G' .or. guesol(1:1)  ==  'g')      rhscrd = rhscrd+i
               if (guesol(2:2)  ==  'X' .or. guesol(2:2)  ==  'x')      rhscrd = rhscrd+i
               rhstyp = 'F'//guesol
          endif 
     end block
!      
     totcrd = ptrcrd + indcrd + valcrd + rhscrd
!      write 4-line or five line header
     write(iounit,10) title,key,totcrd,ptrcrd,indcrd,valcrd,     rhscrd,type,nrow,ncol,nnz,nrhs,ptrfmt,indfmt,valfmt,valfmt
! -----------------------------------------------------------------------
     nrwindx = 0
     if (nrhs  >=  1) write (iounit,11) rhstyp, nrhs, nrwindx

 10 format (a72, a8 / 5i14 / a3, 11x, 4i14 / 2a16, 2a20)
 11 format(a3,11x,i14,i14)
!      
     write(iounit,ptrfmt) (ia (i), i = 1, ncol+1)
     write(iounit,indfmt) (ja (i), i = 1, nnz)
     if (job <= 1) return
     
     write(iounit,valfmt) (a(i), i = 1, nnz)
     
     if (job <= 2) return 
 
     len = nrow*nrhs 
     next = 1
     iend = len
     write(iounit,valfmt) (rhs(i), i = next, iend)
!      
!      write initial guesses if available
!      
     if (guesol(1:1)  ==  'G' .or. guesol(1:1)  ==  'g') then
          next = next+len
          iend = iend+ len
          write(iounit,valfmt) (rhs(i), i = next, iend)
     endif
!      
!      write exact solutions if available
!      
     if (guesol(2:2)  ==  'X' .or. guesol(2:2)  ==  'x') then
          next = next+len
          iend = iend+ len
          write(iounit,valfmt) (rhs(i), i = next, iend)
     endif
!      
     return
! ----------end of prtmt ------------------------------------------------ 
! -----------------------------------------------------------------------
end subroutine prtmt

subroutine dump(i1,i2,values,a,ja,ia,iout)

     integer :: k, iout, maxr, i, k1, k2
     integer, intent(In) :: i1, i2
     integer, dimension(*), intent(In) :: ia, ja 
     real(kind=8), dimension(*), intent(In) :: a
     logical, intent(In) :: values
! -----------------------------------------------------------------------
!  outputs rows i1 through i2 of a sparse matrix stored in CSR format 
!  (or columns i1 through i2 of a matrix stored in CSC format) in a file, 
!  one (column) row at a time in a nice readable format. 
!  This is a simple routine which is useful for debugging. 
! -----------------------------------------------------------------------
!  on entry:
! ---------
!  i1    = first row (column) to print out
!  i2    = last row (column) to print out 
!  values= logical. indicates whether or not to print real values.
!          if value = .false. only the pattern will be output.
!  a,
!  ja, 
!  ia    =  matrix in CSR format (or CSC format) 
!  iout  = logical unit number for output.
! ---------- 
!  the output file iout will have written in it the rows or columns 
!  of the matrix in one of two possible formats (depending on the max 
!  number of elements per row. The values are output with only 
!  two digits of accuracy (D9.2). )
! -----------------------------------------------------------------------
!      local variables
 
! 
!  select mode horizontal or vertical 
! 
     maxr = 0
     do i=i1, i2
          maxr = max0(maxr,ia(i+1)-ia(i))
     end do
     if (maxr <= 8) then
! 
!  able to do one row acros line
! 
          do i=i1, i2
               write(iout,100) i
               k1=ia(i)
               k2 = ia(i+1)-1
               write (iout,101) (ja(k),k=k1,k2)
               if (values) write (iout,102) (a(k),k=k1,k2)
          end do
     else 
! 
!  unable to one row acros line. do three items at a time
!  across a line 
          do i=i1, i2
               if (values) then
                    write(iout,200) i
               else
                    write(iout,203) i               
               endif
               
               k1=ia(i)
               k2 = ia(i+1)-1
 
               if (values) then
                    write (iout,201) (ja(k),a(k),k=k1,k2)
               else
                    write (iout,202) (ja(k),k=k1,k2)
               endif
          end do
     endif 
! 
!  formats :
! 
 100  format (1h ,34(1h-),' row',i6,1x,34(1h-) )
 101  format(' col:',8(i5,6h     : ))
 102  format(' val:',8(d9.2,2h :) )
 200  format (1h ,30(1h-),' row',i3,1x,30(1h-),/     3('  columns :    values  * ') )
! -------------xiiiiiihhhhhhddddddddd-*-
201 format(3(1h ,i6,6h : ,d9.2,3h * ) )
202 format(6(1h ,i5,6h * ) ) 
203  format (1h ,30(1h-),' row',i3,1x,30(1h-),/     3('  column  :  column   *') )
 return
! ----end-of-dump--------------------------------------------------------
! -----------------------------------------------------------------------
end subroutine dump

subroutine pspltm(nrow,ncol,mode,ja,ia,title,ptitle,size,munt,nlines,lines,iunt)
! -----------------------------------------------------------------------

     integer, intent(In) :: nlines, nrow, ncol, ptitle, mode
     integer :: m, kol, isep, iunt, lenstr, n, nr, nc, maxdim, istart, ilast, ii, k, ltit
     real :: haf, zero, lrmrgn, botmrgn,xtit,ytit,ytitof,fnstit,siz,xl,xr,yb,yt, &
               scfct,u2dot,frlw,delt,paperx,conv,xx,yy
     integer, dimension(*), intent(In) :: ja, ia
     integer, dimension(nlines), intent(In) :: lines
     real, intent(In) :: size
     character, intent(In) :: title
     character(len=2), intent(In) :: munt
     logical :: square
    
! ----------------------------------------------------------------------- 
!  PSPLTM - PostScript PLoTer of a (sparse) Matrix
!  This version by loris renggli (renggli@masg1.epfl.ch), Dec 1991
!  and Youcef Saad 
! ------
!  Loris RENGGLI, Swiss Federal Institute of Technology, Math. Dept
!  CH-1015 Lausanne (Switzerland)  -- e-mail:  renggli@masg1.epfl.ch
!  Modified by Youcef Saad -- June 24, 1992 to add a few features:
!  separation lines + acceptance of MSR format.
! -----------------------------------------------------------------------
!  input arguments description :
! 
!  nrow   = number of rows in matrix
! 
!  ncol   = number of columns in matrix 
! 
!  mode   = integer indicating whether the matrix is stored in 
!            CSR mode (mode=0) or CSC mode (mode=1) or MSR mode (mode=2) 
! 
!  ja     = column indices of nonzero elements when matrix is
!           stored rowise. Row indices if stores column-wise.
!  ia     = integer array of containing the pointers to the 
!           beginning of the columns in arrays a, ja.
! 
!  title  = character*(*). a title of arbitrary length to be printed 
!           as a caption to the figure. Can be a blank character if no
!           caption is desired.
! 
!  ptitle = position of title; 0 under the drawing, else above
! 
!  size   = size of the drawing  
! 
!  munt   = units used for size : 'cm' or 'in'
! 
!  nlines = number of separation lines to draw for showing a partionning
!           of the matrix. enter zero if no partition lines are wanted.
! 
!  lines  = integer array of length nlines containing the coordinates of 
!           the desired partition lines . The partitioning is symmetric: 
!           a horizontal line across the matrix will be drawn in 
!           between rows lines(i) and lines(i)+1 for i=1, 2, ..., nlines
!           an a vertical line will be similarly drawn between columns
!           lines(i) and lines(i)+1 for i=1,2,...,nlines 
! 
!  iunt   = logical unit number where to write the matrix into.
! ----------------------------------------------------------------------- 
!  additional note: use of 'cm' assumes european format for paper size
!  (21cm wide) and use of 'in' assumes american format (8.5in wide).
!  The correct centering of the figure depends on the proper choice. Y.S.
! -----------------------------------------------------------------------
 
!  change square to .true. if you prefer a square frame around
!  a rectangular matrix
     data haf / 0.5 / ,zero / 0.0 / ,conv / 2.54 / ,square / .false. / 
! -----------------------------------------------------------------------
     siz = size
     nr = nrow
     nc = ncol
     n = nc
     if (mode == 0) n = nr
!       nnz = ia(n+1) - ia(1) 
     maxdim = max(nrow, ncol)
     m = 1 + maxdim
     nc = nc+1
     nr = nr+1
! 
!  units (cm or in) to dot conversion factor and paper size
!  
     if (munt == 'cm' .or. munt == 'CM') then
          u2dot = 72.0/conv
          paperx = 21.0
     else
          u2dot = 72.0
          paperx = 8.5*conv
          siz = siz*conv
     end if
! 
!  left and right margins (drawing is centered)
!  
     lrmrgn = (paperx-siz)/2.0
! 
!  bottom margin : 2 cm
! 
     botmrgn = 2.0
!  scaling factor
     scfct = siz*u2dot/m
!  matrix frame line witdh
     frlw = 0.25
!  font size for title (cm)
     fnstit = 0.5
     ltit = lenstr(title)
!  position of title : centered horizontally
!                      at 1.0 cm vertically over the drawing
     ytitof = 1.0
     xtit = paperx/2.0
     ytit = botmrgn+siz*nr/m + ytitof
!  almost exact bounding box
     xl = lrmrgn*u2dot - scfct*frlw/2
     xr = (lrmrgn+siz)*u2dot + scfct*frlw/2
     yb = botmrgn*u2dot - scfct*frlw/2
     yt = (botmrgn+siz*nr/m)*u2dot + scfct*frlw/2
     if (ltit > 0) then
          yt = yt + (ytitof+fnstit*0.70)*u2dot
     end if
!  add some room to bounding box
     delt = 10.0
     xl = xl-delt
     xr = xr+delt
     yb = yb-delt
     yt = yt+delt
! 
!  correction for title under the drawing
     if (ptitle == 0 .and. ltit > 0) then
          ytit = botmrgn + fnstit*0.3
          botmrgn = botmrgn + ytitof + fnstit*0.7
     end if
!  begin of output
! 
     write(iunt,10) '%!'
     write(iunt,10) '%%Creator: PSPLTM routine'
     write(iunt,12) '%%BoundingBox:',xl,yb,xr,yt
     write(iunt,10) '%%EndComments'
     write(iunt,10) '/cm {72 mul 2.54 div} def'
     write(iunt,10) '/mc {72 div 2.54 mul} def'
     write(iunt,10) '/pnum { 72 div 2.54 mul 20 string'
     write(iunt,10) 'cvs print ( ) print} def'
     write(iunt,10)  '/Cshow {dup stringwidth pop -2 div 0 rmoveto show} def'
!  needed by editing the output file
! 
!  we leave margins etc. in cm so it is easy to modify them if
     write(iunt,10) 'gsave'
     
     if (ltit > 0) then
          write(iunt,*) '/Helvetica findfont ',fnstit,             ' cm scalefont setfont '
          write(iunt,*) xtit,' cm ',ytit,' cm moveto '
          write(iunt,'(3A)') '(',title(1:ltit),') Cshow'
     end if
     
     write(iunt,*) lrmrgn,' cm ',botmrgn,' cm translate'
     write(iunt,*) siz,' cm ',m,' div dup scale '
! ------- 
!  draw a frame around the matrix
     write(iunt,*) frlw,' setlinewidth'
     write(iunt,10) 'newpath'
     write(iunt,11) 0, 0, ' moveto'
 
     if (square) then
          write(iunt,11) m,0,' lineto'
          write(iunt,11) m, m, ' lineto'
          write(iunt,11) 0,m,' lineto'
     else
          write(iunt,11) nc,0,' lineto'
          write(iunt,11) nc,nr,' lineto'
          write(iunt,11) 0,nr,' lineto'
     end if
          
     write(iunt,10) 'closepath stroke'
! 
!      drawing the separation lines 
!  
     write(iunt,*)  ' 0.2 setlinewidth'
     
     do kol=1, nlines 
          isep = lines(kol) 
! 
!      horizontal lines 
! 
          yy = real(nrow-isep) + haf 
          xx = real(ncol+1) 
          write(iunt,13) zero, yy, ' moveto '
          write(iunt,13)  xx, yy, ' lineto stroke '
! 
!  vertical lines 
! 
          xx = real(isep) + haf 
          yy = real(nrow+1) 
          write(iunt,13) xx, zero,' moveto '
          write(iunt,13) xx, yy, ' lineto stroke '             
     end do
!  
! ----------- plotting loop ---------------------------------------------
! 
     write(iunt,10) '1 1 translate'
     write(iunt,10) '0.8 setlinewidth'
     write(iunt,10) '/p {moveto 0 -.40 rmoveto '
     write(iunt,10) '           0  .80 rlineto stroke} def'
!      
     do ii=1, n
          istart = ia(ii)
          ilast = ia(ii+1)-1 
          if (mode == 1) then
               do k=istart, ilast
                    write(iunt,11) ii-1, nrow-ja(k), ' p'
               end do
          else
               do k=istart, ilast
                    write(iunt,11) ja(k)-1, nrow-ii, ' p'
               end do
!  add diagonal element if MSR mode.
               if (mode  ==  2)          write(iunt,11) ii-1, nrow-ii, ' p' 
! 
          endif
     end do
! -----------------------------------------------------------------------
     write(iunt,10) 'showpage'
     return
! 
 10 format (a)
 11 format (2(i6,1x),a)
 12 format (a,4(1x,f9.2))
 13 format (2(f9.2,1x),a)
! -----------------------------------------------------------------------
end subroutine pspltm

integer function lenstr(s)
! -----------------------------------------------------------------------
!  return length of the string S
! -----------------------------------------------------------------------
     character(len=*) :: s
     integer :: len, n
! ----------------------------------------------------------------------- 
     n = len(s)

     do while (n > 0)
          if (s(n:n) == ' ') n = n - 1
     end do
     lenstr = n
! 
     return
! --------end-of-pspltm--------------------------------------------------
! -----------------------------------------------------------------------
end function lenstr

subroutine pltmt(nrow,ncol,mode,ja,ia,title,key,type,job,iounit)
     integer, intent(In) :: nrow, ncol, mode, job
     integer :: iounit, n, nnz, maxdim, ips, ii, istart,ilast,k
     real :: x,y, xnrow, ptsize, hscale, vscale, xwid, xht, xshift,yshift,tiny
     integer, dimension(*), intent(In) :: ja, ia
     character(len=8), intent(In) :: key
     character(len=72), intent(In) :: title
     character(len=3), intent(In) :: type
! -----------------------------------------------------------------------
!  this subroutine creates a 'pic' file for plotting the pattern of 
!  a sparse matrix stored in general sparse format. it is not intended
!  to be a means of plotting large matrices (it is very inefficient).
!  It is however useful for small matrices and can be used for example
!  for inserting matrix plots in a text. The size of the plot can be
!  7in x 7in or 5 in x 5in .. There is also an option for writing a 
!  3-line header in troff (see description of parameter job).
!  Author: Youcef Saad - Date: Sept., 1989
!  See SPARSKIT/UNSUPP/ for a version of this to produce a post-script
!  file. 
! -----------------------------------------------------------------------
!  nrow   = number of rows in matrix
! 
!  ncol	 = number of columns in matrix 
! 
!  mode   = integer indicating whether the matrix is stored
!           row-wise (mode = 0) or column-wise (mode=1)
! 
!  ja     = column indices of nonzero elements when matrix is
!         stored rowise. Row indices if stores column-wise.
!  ia     = integer array of containing the pointers to the 
!         beginning of the columns in arrays a, ja.
! 
!  title  = character*71 = title of matrix test ( character a*71 ).
!  key    = character*8  = key of matrix 
!  type   = character*3  = type of matrix. 
! 
!  job    = this integer parameter allows to set a few minor 
!           options. First it tells pltmt whether or not to
!           reduce the plot. The standard size of 7in is then
!           replaced by a 5in plot. It also tells pltmt whether or
!           not to append to the pic file a few 'troff' lines that 
!           produce a centered caption includingg the title, key and 
!           types as well as the size and number of nonzero elements.
!           job =  0 : do not reduce and do not make caption.
!           job =  1 : reduce and do not make caption.
!           job = 10 : do not reduce and make caption
!           job = 11 : reduce and make caption.
!           (i.e. trailing digit for reduction, leading digit for caption)
! 
!  iounit = logical unit number where to write the matrix into.
! 
! -----------------------------------------------------------------------
!  example of usage . 
! -----------------
!  In the fortran code:
!   a) read a Harwell/Boeing matrix
!           call readmt (.....)
!         iout = 13
!   b) generate pic file:
!           call  pltmt (nrow,ncol,mode,ja,ia,title,key,type,iout)
!         stop
!  ---------
!  Then in a unix environment plot the matrix by the command
! 
!      pic FOR013.DAT | troff -me | lpr -Ppsx
! 
! -----------------------------------------------------------------------
!  notes: 1) Plots square as well as rectangular matrices.
!             (however not as much tested with rectangular matrices.)
!        2) the dot-size is adapted according to the size of the
!             matrix.
!        3) This is not meant at all as a way of plotting large
!             matrices. The pic file generaled will have one line for
!             each nonzero element. It is  only meant for use in
!           such things as document poreparations etc..
!          4) The caption written will print the 71 character long
!             title. This may not be centered correctly if the
!             title has trailing blanks (a problem with Troff).
!             if you want the title centered then you can center
!             the string in title before calling pltmt. 
!        
! -----------------------------------------------------------------------
 

     n = ncol
     if (mode == 0) n = nrow
     nnz = ia(n+1) - ia(1) 
     maxdim = max0 (nrow, ncol)
     xnrow = real(nrow)
     ptsize = 0.08
     hscale = (7.0 -2.0*ptsize)/real(maxdim-1) 
     vscale = hscale 
     xwid = ptsize + real(ncol-1)*hscale + ptsize
     xht = ptsize + real(nrow-1)*vscale + ptsize
     xshift = (7.0-xwid)/2.0
     yshift = (7.0-xht)/2.0 
! ------
     if (mod(job,10) == 1) then
          write (iounit,88)
     else
          write (iounit,89)
     endif
     
     write(iounit,90) 
     write(iounit,91) xwid, xht, xshift, yshift

!      
!      
!      shift points slightly to account for size of dot , etc..
     tiny = 0.03
     if (mod(job,10) == 1) tiny = 0.05
     xshift = xshift + ptsize - tiny
     yshift = yshift + ptsize + tiny
!      
! -----------------------------------------------------------------------
!      
     ips = 8
     if (maxdim <= 500) ips = 10
     if (maxdim <= 300) ips = 12
     if (maxdim <= 100) ips = 16
     if (maxdim < 50) ips = 24
     write(iounit,92) ips
 
!      
! -----------plottingloop --------------------------------------------- 
!      
     do ii=1, n
          istart = ia(ii)
          ilast = ia(ii+1)-1 
          if (mode /= 0) then
               x = real(ii-1)
               do k=istart, ilast
                    y = xnrow-real(ja(k))
                    write(iounit,128) xshift+x*hscale, yshift+y*vscale
               end do
          else
               y = xnrow - real(ii)
               do k=istart, ilast
                    x = real(ja(k)-1)
                    write(iounit,128) xshift+x*hscale, yshift+y*vscale
               end do
          endif
     end do
     write (iounit, 129)
     !      quit if caption not desired. 
     if ( (job/10) /= 1) return
 !      
     write(iounit,127) key, type, title
     write(iounit,130) nrow,ncol,nnz
 
! -----------------------------------------------------------------------
 88   format('.PS 5in',/,'.po 1.8i')
 89   format('.PS',/,'.po 0.7i')
 90   format('box invisible wid 7.0 ht 7.0 with .sw at (0.0,0.0) ') 
 91   format('box wid ',f5.2,' ht ',f5.2,     ' with .sw at (',f5.2,',',f5.2,')' )
 92   format ('.ps ',i2)
 128  format(7h"." at ,f6.3,1h,,f6.3,8h ljust  )
     
 129  format('.PE')
 127  format('.sp 4'/'.ll 7i'/'.ps 12'/'.po 0.7i'/'.ce 3'/,     'Matrix:  ',a8,',  Type:  ',a3,/, &
      a72)
 130  format('Dimension: ',i4,' x ',i4,',  Nonzero elements: ',i5)
     return
! ----------------end-of-pltmt ------------------------------------------
! ----------------------------------------------------------------------- 
end subroutine pltmt

subroutine smms(n,first,last,mode,a,ja,ia,iout)

     integer :: i, k1, k2, k, iout
     integer, dimension(*), intent(In) :: ia, ja
     integer, intent(In) :: n, first, last, mode
     real(kind=8), dimension(*), intent(In) :: a
     logical :: msr, csc
! -----------------------------------------------------------------------
!  writes a matrix in Coordinate (SMMS) format -- 
! -----------------------------------------------------------------------
!  on entry:
! ---------
!  n     = integer = size of matrix -- number of rows (columns if matrix
!          is stored columnwise) 
!  first  = first row (column) to be output. This routine will output 
!           rows (colums) first to last. 
!  last   = last row (column) to be output. 
!  mode   = integer giving some information about the storage of the 
!           matrix. A 3-digit decimal number. 'htu' 
!          * u = 0 means that matrix is stored row-wise 
!          * u = 1 means that matrix is stored column-wise 
!          * t = 0 indicates that the matrix is stored in CSR format 
!          * t = 1 indicates that the matrix is stored in MSR format. 
!          * h = ... to be added. 
!  a,
!  ja,
!  ia    =  matrix in CSR or MSR format (see mode) 
!  iout  = output unit number.
! 
!  on return:
! ----------
!  the output file iout will have written in it the matrix in smms
!  (coordinate format)
! 
! -----------------------------------------------------------------------
 
! 
!  determine mode ( msr or csr )
! 
     msr = .false.
     csc = .false. 
     if (mod(mode,10) == 1) csc = .true.
     if ( (mode/10) == 1) msr = .true. 
     
     write (iout,*) n
 
     do i=first, last 
          k1=ia(i)
          k2 = ia(i+1)-1
!      write (iout,*) ' row ', i 
          if (msr) write(iout,'(2i6,e22.14)')  i, i, a(i) 
          do k=k1, k2
               if (csc) then
                    write(iout,'(2i6,e22.14)')  ja(k), i, a(k)
               else
                    write(iout,'(2i6,e22.14)')  i, ja(k), a(k)
               endif 
          end do
     end do
! ----end-of-smms--------------------------------------------------------
! -----------------------------------------------------------------------
end subroutine smms

subroutine readsm(nmax,nzmax,n,nnz,ia,ja,a,iout,ierr)
     integer, intent(Out) :: nnz, ierr
     integer, intent(In) :: nmax, nzmax
     integer :: row, n, iout, i, j, k
     integer, dimension(nmax + 1), intent(Out) :: ia
     integer, dimension(nzmax), intent(Out) :: ja
     real(kind=8), dimension(nzmax), intent(Out) :: a
     real(kind=8) :: x
! -----------------------------------------------------------------------
!      read a matrix in coordinate format as is used in the SMMS
!      package (F. Alvarado), i.e. the row is in ascending order.
!      Outputs the matrix in CSR format.
! -----------------------------------------------------------------------
!  coded by Kesheng Wu on Oct 21, 1991 with the supervision of Y. Saad
! -----------------------------------------------------------------------
!  on entry:
! ---------
!  nmax  = the maximum size of array
!  nzmax = the maximum number of nonzeros
!  iout  = the I/O unit that has the data file
! 
!  on return:
! ---------- 
!  n     = integer = size of matrix
!  nnz   = number of non-zero entries in the matrix
!  a,
!  ja, 
!  ia    = matrix in CSR format
!  ierr  = error code,
!          0 -- subroutine end with intended job done
!          1 -- error in I/O unit iout
!          2 -- end-of-file reached while reading n, i.e. a empty data file
!          3 -- n non-positive or too large
!          4 -- nnz is zero or larger than nzmax
!          5 -- data file is not orgnized in the order of ascending
!               row indices
! 
!  in case of errors:
!    n will be set to zero (0). In case the data file has more than nzmax
!    number of entries, the first nzmax entries will be read, and are not 
!    cleared on return. The total number of entry is determined.
!    Ierr is set.
! -----------------------------------------------------------------------
! 
     rewind(iout)
     nnz = 0
     ia(1) = 1
     row = 1
! 
     read (iout,*, err=1000, end=1010) n
     
     if ((n <= 0) .or. (n > nmax)) then
          ierr  = 3
          n = 0
          return
     end if

!    
     do while (nnz < nzmax)
          nnz = nnz + 1
          read (iout, *, err=1000, end=100) i, j, x
          if (i > row) then
               do k = row+1, i
                    ia(k) = nnz
               end do
               row = i
          else if (i < row) then
               ierr = 5
               n = 0
               return
          endif
          ja(nnz) = j
          a (nnz) = x
     end do

     if (nnz >= nzmax) then
          ierr = 4
! 
!      try to determine the real number of entries, in case needed
! 
          do
               read(iout,*,err = 210, end = 210) i,j,x
               nnz = nnz + 1
          end do
210       continue
          n = 0
          return
     endif

!      normal return -- end of file reached
100  ia(row+1) = nnz
     nnz = nnz - 1
     if (nnz == 0) then
          ierr = 4
          do
               read(iout,*,err = 220, end = 220) i,j,x
               nnz = nnz + 1
          end do
220       continue
          n = 0
          return
     end if

! 
!      everything seems to be OK.
! 
     ierr = 0
     return
! 
!      error handling code
! 
!      error in reading data entries
! 
1000 ierr = 1
     if (ierr == 1) then
          n = 0
          return
     end if
! 
!      empty file
! 
1010 ierr = 2
     if (ierr == 2) then
          n = 0
          return
     end if
! ----end-of-readsm------------------------------------------------------
! -----------------------------------------------------------------------
end subroutine readsm

subroutine readsk(nmax,nzmax,n,nnz,a,ja,ia,iounit,ierr)

     integer, intent(In) :: nmax, nzmax
     integer :: iounit, n, nnz, i
     integer, intent(Out) :: ierr
     integer, dimension(nmax + 1) :: ia
     integer, dimension(nzmax) :: ja
     real(kind=8), dimension(nzmax) :: a
! -----------------------------------------------------------------------
!  Reads matrix in Compressed Saprse Row format. The data is supposed to
!  appear in the following order -- n, ia, ja, a
!  Only square matrices accepted. Format has following features
!  (1) each number is separated by at least one space (or end-of-line), 
!  (2) each array starts with a new line.
! -----------------------------------------------------------------------
!  coded by Kesheng Wu on Oct 21, 1991 with supervision of Y. Saad
! -----------------------------------------------------------------------
!  on entry:
! ---------
!  nmax 	 = max column dimension  allowed for matrix.
!  nzmax	 = max number of nonzeros elements allowed. the arrays a, 
!           and ja should be of length equal to nnz (see below).
!  iounit = logical unit number where to read the matrix from.
! 
!  on return:
! ---------- 
!  ia,
!  ja,
!  a      = matrx in CSR format
!  n      = number of rows(columns) in matrix
!  nnz	 = number of nonzero elements in A. This info is returned
!           even if there is not enough space in a, ja, ia, in order
!           to determine the minimum storage needed.
!  ierr   = error code,
!           0 : OK;
!           1 : error when try to read the specified I/O unit.
!           2 : end-of-file reached during reading of data file.
!           3 : array size in data file is negtive or larger than nmax;
!           4 : nunmer of nonzeros in data file is negtive or larger than nzmax
!  in case of errors:
! ---------
!      n is set to 0 (zero), at the same time ierr is set.
! -----------------------------------------------------------------------
!      
!      read the size of the matrix
! 
     rewind(iounit)
     read (iounit, *, err=1000, end=1010) n
     if ((n <= 0).or.(n > nmax)) then ! n non-positive or too large
          ierr = 3
          n = 0
          return
     end if
!      
!      read the pointer array ia(*)
!      
     read (iounit, *, err=1000, end=1010) (ia(i), i=1, n+1)
!      
!      Number of None-Zeros
!      
     nnz = ia(n+1) - 1
     if ((nnz <= 0).or.(nnz > nzmax)) then !NNZ non-positive or too large
          ierr = 4
          n = 0
          return
     end if
!      
!      read the column indices array
!      
     read (iounit, *, err=1000, end=1010) (ja(i), i=1, nnz)
!      
!      read the matrix elements
!      
     read (iounit, *, err=1000, end=1010) (a(i), i=1, nnz)
!      
!      normal return
! 
     ierr = 0
     return
!      
!      error handling code
!      
!      error in reading I/O unit
1000      ierr = 1
          if (ierr == 1) then
               n = 0
               return
          end if

! 
!      EOF reached in reading
1010      ierr = 2
          if(ierr == 2) Then
               n = 0
               return
          end if
! 
! 
!      
! ---------end of readsk ------------------------------------------------
! -----------------------------------------------------------------------
end subroutine readsk

subroutine skit(n,a,ja,ia,ifmt,iounit,ierr)
! -----------------------------------------------------------------------
!      Writes a matrix in Compressed Sparse Row format to an I/O unit.
!      It tryes to pack as many number as possible into lines of less than
!      80 characters. Space is inserted in between numbers for separation
!      to avoid carrying a header in the data file. This can be viewed
!      as a simplified Harwell-Boeing format. 
! -----------------------------------------------------------------------
!  Modified from subroutine prtmt written by Y. Saad
! -----------------------------------------------------------------------
!  on entry:
! ---------
!  n      = number of rows(columns) in matrix
!  a      = real*8 array containing the values of the matrix stored 
!           columnwise
!  ja     = integer array of the same length as a containing the column
!           indices of the corresponding matrix elements of array a.
!  ia     = integer array of containing the pointers to the beginning of 
!           the row in arrays a and ja.
!  ifmt   = integer specifying the format chosen for the real values
!           to be output (i.e., for a, and for rhs-guess-sol if 
!           applicable). The meaning of ifmt is as follows.
!           * if (ifmt .lt. 100) then the D descriptor is used,
!           format Dd.m, in which the length (m) of the mantissa is 
!           precisely the integer ifmt (and d = ifmt+6)
!           * if (ifmt .gt. 100) then prtmt will use the 
!           F- descriptor (format Fd.m) in which the length of the 
!           mantissa (m) is the integer mod(ifmt,100) and the length 
!           of the integer part is k=ifmt/100 (and d = k+m+2)
!           Thus  ifmt= 4   means  D10.4  +.xxxxD+ee    while
!           ifmt=104  means  F7.4   +x.xxxx
!           ifmt=205  means  F9.5   +xx.xxxxx
!           Note: formats for ja, and ia are internally computed.
!      
!  iounit = logical unit number where to write the matrix into.
!      
!  on return:
! ----------
!  ierr   = error code, 0 for normal 1 for error in writing to iounit.
! 
!  on error:
! --------
!      If error is encontacted when writing the matrix, the whole matrix
!      is written to the standard output.
!      ierr is set to 1.
! -----------------------------------------------------------------------

     character(len=16) :: ptrfmt, indfmt
     character(len=20) :: valfmt
     integer :: iounit, n, len, nperli, nnz, i, ihead
     integer, intent(InOut) :: ifmt
     integer, dimension(*), intent(In) :: ja, ia
     integer, intent(Out) :: ierr
     real(kind=8), dimension(*), intent(In) :: a
     character(len=28) :: ix
! --------------
!      compute pointer format
! --------------
     nnz = ia(n+1)
     len = int ( alog10(0.1+real(nnz))) + 2
     nnz = nnz - 1
     nperli = 80/len
     print *, ' skit entries:', n, nnz, len, nperli
     if (len > 9) then
          ix = '(1h(,i2,1HI,i2,1h) )'
     else
          ix = '(1h(,i2,1HI,i1,1h) )'
     endif
     write (ptrfmt,ix) nperli,len
! ----------------------------
!      compute ROW index format
! ----------------------------
     len = int ( alog10(0.1+real(n) )) + 2
     nperli = min0(80/len,nnz)
     write (indfmt,'(1h(,i2,1HI,i1,1h) )') nperli,len
! ---------------------------
!      compute value format
! ---------------------------
     if (ifmt >= 100) then
          ihead = ifmt/100
          ifmt = ifmt-100*ihead
          len = ihead+ifmt+3
          nperli = 80/len
!      
          if (len <= 9 ) then
               ix = '(1h(,i2,1hF,i1,1h.,i1,1h) )'
          elseif (ifmt <= 9) then
               ix = '(1h(,i2,1hF,i2,1h.,i1,1h) )'
          else 
               ix = '(1h(,i2,1hF,i2,1h.,i2,1h) )'
          endif
!      
          write(valfmt,ix) nperli,len,ifmt
     else
          len = ifmt + 7
          nperli = 80/len
!      try to minimize the blanks in the format strings.
          if (nperli <= 9) then
               if (len <= 9 ) then
                    ix = '(1h(,i1,1hD,i1,1h.,i1,1h) )'
               elseif (ifmt <= 9) then
                    ix = '(1h(,i1,1hD,i2,1h.,i1,1h) )'
          else 
               ix = '(1h(,i1,1hD,i2,1h.,i2,1h) )'
          endif
     else 
          if (len <= 9 ) then
               ix = '(1h(,i2,1hD,i1,1h.,i1,1h) )'
          elseif (ifmt <= 9) then
               ix = '(1h(,i2,1hD,i2,1h.,i1,1h) )'
          else 
               ix = '(1h(,i2,1hD,i2,1h.,i2,1h) )'
          endif
     endif
! -----------
          write(valfmt,ix) nperli,len,ifmt
     endif 
!      
!      output the data
!      
     write(iounit, *) n
     write(iounit,ptrfmt,err=1000) (ia(i), i = 1, n+1)
     write(iounit,indfmt,err=1000) (ja(i), i = 1, nnz)
     write(iounit,valfmt,err=1000) ( a(i), i = 1, nnz)
! 
!      done, if no trouble is encounted in writing data
! 
     ierr = 0
     return
!      
!      if can't write the data to the I/O unit specified, should be able to
!      write everything to standard output (unit 6)
!      
1000 write(0, *) 'Error, Can''t write data to sepcified unit',iounit
     write(0, *) 'Write the matrix into standard output instead!'
     ierr = 1
     write(6,*) n
     write(6,ptrfmt) (ia(i), i=1, n+1)
     write(6,indfmt) (ja(i), i=1, nnz)
     write(6,valfmt) ( a(i), i=1, nnz)
     return
! ----------end of skit ------------------------------------------------- 
! -----------------------------------------------------------------------
end subroutine skit

subroutine prtunf(n,a,ja,ia,iout,ierr)
! -----------------------------------------------------------------------
!  This subroutine dumps the arrays used for storing sparse compressed row
!  format in machine code, i.e. unformatted using standard FORTRAN term.
! -----------------------------------------------------------------------
!  First coded by Kesheng Wu on Oct 21, 1991 under the instruction of
!  Prof. Y. Saad
! -----------------------------------------------------------------------
!  On entry:
!      n: the size of the matrix (matrix is n X n)
!     ia: integer array stores the stariting position of each row.
!     ja: integer array stores the column indices of each entry.
!      a: the non-zero entries of the matrix.
!   iout: the unit number opened for storing the matrix.
!  On return:
!   ierr: a error, 0 if everything's OK, else 1 if error in writing data.
!  On error:
!   set ierr to 1.
!   No redirection is made, since direct the machine code to the standard
!  output may cause unpridictable consequences.
! -----------------------------------------------------------------------

     integer :: k, iout, nnz
     integer, intent(In) :: n
     integer, intent(Out) :: ierr
     integer, dimension(*), intent(In) :: ia, ja
     real(kind=8), dimension(*), intent(In) :: a
     
     nnz = ia(n+1)-ia(1) 
! 
     write(unit=iout, err=1000)  n
     write(unit=iout, err=1000) (ia(k),k=1,n+1) 
     if (nnz > 0) then
          write(unit=iout, err=1000) (ja(k),k=1,nnz)
          write(unit=iout, err=1000) ( a(k),k=1,nnz)
     endif 
! 
     ierr = 0
     return
! 
1000 ierr = 1
     return
end subroutine prtunf

subroutine readunf(nmax,nzmax,n,nnz,a,ja,ia,iounit,ierr)
! -----------------------------------------------------------------------
!  This subroutine reads a matix store in machine code (FORTRAN
!  unformatted form). The matrix is in CSR format.
! -----------------------------------------------------------------------
!  First coded by Kesheng Wu on Oct 21, 1991 under the instruction of
!  Prof. Y. Saad
! -----------------------------------------------------------------------
!  On entry:
!     nmax: the maximum value of matrix size.
!    nzmax: the maximum number of non-zero entries.
!   iounit: the I/O unit that opened for reading.
!  On return:
!        n: the actual size of array.
!      nnz: the actual number of non-zero entries.
!  ia,ja,a: the matrix in CSR format.
!     ierr: a error code, it's same as that used in reaadsk
!           0 -- OK
!           1 -- error in reading iounit
!           2 -- end-of-file reached while reading data file
!           3 -- n is non-positive or too large
!           4 -- nnz is non-positive or too large
!  On error:
!      return with n set to 0 (zero). nnz is kept if it's set already,
!      in case one want to use it to determine the size of array needed
!      to hold the data.
! -----------------------------------------------------------------------
! 
     integer, intent(Out) :: ierr
     integer, intent(In) :: nmax, nzmax
     integer :: n, ioutit, nnz, k
     integer, dimension(nmax + 1) :: ia, ja
     real(kind=8), dimension(nzmax) :: a
! 
     rewind(iounit)
!       
     read (unit=iounit, err=1000, end=1010) n
     
     if ((n <= 0) .or. (n > nmax)) then
          ierr = 3
          n = 0
          return
     end if
! 
     read(unit=iounit, err=1000, end=1010) (ia(k),k=1,n+1)
! 
     nnz = ia(n+1) - 1
     if ((nnz <= 0) .or. (nnz > nzmax)) then
          ierr = 4
          n = 0
          return
     end if
! 
     read(unit=iounit, err=1000, end=1010) (ja(k),k=1,nnz)
     read(unit=iounit, err=1000, end=1010) (a(k),k=1,nnz)
! 
!      everything seems to be OK.
! 
     ierr = 0
     return
! 
!      error handling
! 
1000 ierr = 1
     if (ierr == 1) then
          n = 0
          return
     end if
1010 ierr = 2
     if (ierr == 2) then
          n = 0
          return
     end if
 end subroutine readunf