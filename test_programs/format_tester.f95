!*==chkfmt.f90 processed by SPAG 8.04RA 23:02 22 Feb 2024
!!SPAG Open source Personal, Educational or Academic User Clemson University  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
PROGRAM chkfmt
!----------------------------------------------------------------------c
!                          S P A R S K I T                             c
!----------------------------------------------------------------------c
! test suite for the Formats routines.                                 c
! tests all of the routines in the module formats.                     c
!----------------------------------------------------------------------c
! Note: the comments may not have been updated.
!
! Here is the sequence of what is done by this program.
! 1) call gen57bl togenerate a block matrix associated with a simple
!     5-point matrix on a 4 x 2 grid (2-D) with 2 degrees of freedom
!     per grid point.  Thus N = 16. This is produced in BSR format.
!    the pattern of the reduced  matrix is written in csr.mat.
!
! 2) the block format is translated into a compressed sparse row
!    format by bsrcsr. The result is dumped in file csr.mat
!    matrix is converted back to bsr format. block pattern shown again
!
! 3) the matrix is translated in dense format by csrdns.
!    result a 16 x 16 matrix is written in unit dns.mat.
!    This is a good file to look at to see what the matrix is
!    and to compare results of other formats with.
! 4) the dense matrix obtained in 3) is reconverted back to
!    csr format using dnscsr. Result appended to file csr.mat
! 5) The matrix obtained in 4) is converted in coordinate format
!    and the resulting matrix is written in file coo.mat
! 6) the result is converted back to csr format. matrix
!    appended to csr.mat.
! 7) result of 6) is converted to symmetric sparse row storage
!    (ssr) and the result is appended to csr.mat
! 8) result of 7) converted back to csr format and result is
!    appended to csr.mat
! 9) matrix resulting from 8) is converted to modified sparse
!    row format using csrmsr and result is written in msr.mat.
!10) the resulting matrix is converted back to csrformat and
!    result is appended to csr.mat
!11) result is converted to ellpack-itpack format with
!    csrell and result is printed in itp.mat
!12) result is converted back to csr format and appended to csr.mat
!12) result converted to csc format (transposition) using csrcsc
!    which should produce the same matrix here. result appended
!    to csr.mat. A second call to csrcsc is made on resulting
!    matrix.
!13) the subroutine csrdia is used to extract two diagonals
!    (offsets -1 and 0) and then all the diagonals of matrix.
!    results in dia.mat
!14) diacsr is then called to convert the diagonally stored matrix
!    back to csr format. result appended to csr.mat
!15) result is converted to band format (bnd) by calling
!    csrbnd. result dumped to bnd.mat
!16) result is converted back to csr format and appended to csr.mat
!17) result sorted by a call to csrcsc and then converted to
!    block format (csrbsr) and then back to csr format again.
!    result appedned to csr.mat.
!18) matrix converted to symmetric skyline format. result appended
!    to file band.mat
!19) matrix converted back to csr format and result appended to
!    csr.mat.
!20) result converted to jad format. result output in jad.mat
!21) result concverted back to csr fromat. appended to csr.mat
!------------------------------------------------------------------
   USE ISO_FORTRAN_ENV                 
   IMPLICIT NONE
!
! PARAMETER definitions rewritten by SPAG
!
   INTEGER , PARAMETER :: NXMAX = 10 , NMX = NXMAX*NXMAX , NNZMAX = 10*NMX
!
! Local variable declarations rewritten by SPAG
!
   REAL(REAL64) , DIMENSION(NNZMAX) :: a , a1
   REAL(REAL64) , DIMENSION(20,20) :: dns
   INTEGER :: i , idiag , idiag0 , ierr , imod , iout , j , job , k , k1 , k2 , kdiag , kend , kstart , len , lowd , maxcol , ml , &
            & mu , n , na , ndiag , nfree , nnz , nr , nx , ny , nz
   INTEGER , DIMENSION(NMX+1) :: ia , iwk2
   INTEGER , DIMENSION(NNZMAX) :: ia1 , ja , ja1
   INTEGER , DIMENSION(20) :: ioff
   INTEGER , DIMENSION(NMX*2+1) :: iwk
   INTEGER , SAVE :: ndns
   REAL(REAL64) , DIMENSION(7,100) :: stencil
   REAL(REAL64) , DIMENSION(NMX) :: wk
   EXTERNAL bndcsr , bsrcsr , coocsr , csrbnd , csrbsr , csrcoo , csrcsc , csrdia , csrdns , csrell , csrjad , csrlnk , csrmsr ,   &
          & csrssk , csrssr , diacsr , dnscsr , dump , ellcsr , gen57bl , infdia , jadcsr , lnkcsr , msrcsr , sskssr , ssrcsr
!
! End of declarations rewritten by SPAG
!
!
!-----------------------------------------------------------------------
   DATA ndns/20/
!----- open statements ----------------
   OPEN (UNIT=7,FILE='csr.mat')
   OPEN (UNIT=8,FILE='dns.mat')
   OPEN (UNIT=9,FILE='coo.mat')
   OPEN (UNIT=10,FILE='msr.mat')
   OPEN (UNIT=11,FILE='itp.mat')
   OPEN (UNIT=12,FILE='dia.mat')
   OPEN (UNIT=13,FILE='bnd.mat')
   OPEN (UNIT=14,FILE='jad.mat')
!
!---- dimension of grid
!
   nx = 4
   ny = 2
   nz = 1
   nfree = 2
!
!---- generate grid problem.
!
!      na = nx*ny*nz*5
   na = nfree*nfree
   CALL gen57bl(nx,ny,nz,nfree,na,nr,a1,ja1,ia1,iwk,stencil)
!      nr = n / nfree
   n = nr*nfree
!
!---- dump the reduced matrix
!
   iout = 7
   nnz = ia(n+1) - 1
   WRITE (iout,*) '-----------------------------------------'
   WRITE (iout,*) '  +++ Pattern of block matrix  (CSR) +++ '
   WRITE (iout,*) '-----------------------------------------'
!
   CALL dump(1,nr,.FALSE.,a1,ja1,ia1,7)
!
!---- convert to CSR format and dump result
!
   CALL bsrcsr(1,nr,nfree,na,a1,ja1,ia1,a,ja,ia)
   iout = 7
   nnz = ia(n+1) - 1
   WRITE (iout,*) '-----------------------------------------'
   WRITE (iout,*) '  +++  initial matrix in CSR format +++ '
   WRITE (iout,*) '-----------------------------------------'
   CALL dump(1,n,.TRUE.,a,ja,ia,7)
!
!---- convert back to BSR format and dump pattern again.
!
   CALL csrbsr(1,n,nfree,na,a,ja,ia,a1,ja1,ia1,iwk2,ierr)
   WRITE (iout,*) '-----------------------------------------'
   WRITE (iout,*) '  +++ Pattern of block matrix  (CSR) +++ '
   WRITE (iout,*) '-----------------------------------------'
   CALL dump(1,nr,.FALSE.,a1,ja1,ia1,7)
!
!----- convert to BSR format.
!
   CALL bsrcsr(1,nr,nfree,na,a1,ja1,ia1,a,ja,ia)
   iout = 7
   nnz = ia(n+1) - 1
   WRITE (iout,*) '-----------------------------------------'
   WRITE (iout,*) ' +++  matrix after BSRCSR conversion +++ '
   WRITE (iout,*) '-----------------------------------------'
   CALL dump(1,n,.TRUE.,a,ja,ia,7)
!
!----- convert to dense format.
!
   CALL csrdns(n,n,a,ja,ia,dns,ndns,ierr)
!
   iout = iout + 1
   WRITE (iout,*) '-----------------------------------------'
   WRITE (iout,*) '  +++  initial matrix in DENSE format+++ '
   WRITE (iout,*) '-----------------------------------------'
   WRITE (iout,'(4x,16i4)') (j,j=1,n)
   WRITE (iout,'(3x,65(1h-))')
   DO i = 1 , n
      WRITE (8,99001) i , (dns(i,j),j=1,n)
   ENDDO
!
!----- convert back to sparse format.
!
   CALL dnscsr(n,n,NNZMAX,dns,ndns,a1,ja1,ia1,ierr)
   WRITE (7,*) '-----------------------------------------'
   WRITE (7,*) '  +++ matrix after conversion from dnscsr +++ '
   WRITE (7,*) '-----------------------------------------'
   IF ( ierr/=0 ) WRITE (7,*) ' ***** ERROR FROM DNSCSR'
   IF ( ierr/=0 ) WRITE (7,*) '     IERR = ' , ierr
   CALL dump(1,n,.TRUE.,a1,ja1,ia1,7)
!
!     convert it to coordinate format.
!
   CALL csrcoo(n,3,NNZMAX,a,ja,ia,nnz,a1,ia1,ja1,ierr)
   iout = iout + 1
   IF ( ierr/=0 ) WRITE (iout,*) ' ***** ERROR IN CSRCOO'
   IF ( ierr/=0 ) WRITE (iout,*) '     IERR = ' , ierr
   WRITE (iout,*) '-----------------------------------------'
   WRITE (iout,*) ' +++ Matrix in coordinate format +++ '
   WRITE (iout,*) '-----------------------------------------'
   WRITE (iout,99002) (ia1(j),ja1(j),a1(j),j=1,nnz)
99002 FORMAT (' i =',i3,'    j = ',i3,'     a(i,j) = ',f4.1)
!
!     convert it back again to csr format
!
   CALL coocsr(n,nnz,a1,ia1,ja1,a,ja,ia)
 
   WRITE (7,*) '-----------------------------------------'
   WRITE (7,*) '  +++ matrix after conversion from coocsr +++ '
   WRITE (7,*) '-----------------------------------------'
   CALL dump(1,n,.TRUE.,a,ja,ia,7)
!
!     going to srs format
!
   CALL csrssr(n,a,ja,ia,NNZMAX,a1,ja1,ia1,ierr)
   WRITE (7,*) '-----------------------------------------'
   WRITE (7,*) '  +++ matrix after conversion to ssr format +++ '
   WRITE (7,*) '      (lower part only stored in csr format)    '
   WRITE (7,*) '-----------------------------------------'
   CALL dump(1,n,.TRUE.,a1,ja1,ia1,7)
!     back to csr
   CALL ssrcsr(3,1,n,a1,ja1,ia1,NNZMAX,a,ja,ia,iwk,iwk2,ierr)
   IF ( ierr/=0 ) WRITE (7,*) ' error in ssrcsr-IERR=' , ierr
   WRITE (7,*) '-----------------------------------------'
   WRITE (7,*) '  +++ matrix after conversion from ssrcsr +++ '
   WRITE (7,*) '-----------------------------------------'
   CALL dump(1,n,.TRUE.,a,ja,ia,7)
!---- msr format
   iout = iout + 1
   CALL csrmsr(n,a,ja,ia,a1,ja1,a1,ja1)
   WRITE (iout,*) '-----------------------------------------'
   WRITE (iout,*) '  +++ matrix in modified sparse row format +++'
   WRITE (iout,*) '-----------------------------------------'
   WRITE (iout,*) ' ** MAIN DIAGONAL '
   WRITE (iout,'(16f4.1)') (a1(k),k=1,n)
   WRITE (iout,*) ' ** POINTERS: '
   WRITE (iout,'(17i4)') (ja1(k),k=1,n+1)
   WRITE (iout,*) ' ** REMAINDER :'
   CALL dump(1,n,.TRUE.,a1,ja1,ja1,iout)
!-------
   CALL msrcsr(n,a1,ja1,a,ja,ia,wk,iwk2)
   WRITE (7,*) '-----------------------------------------'
   WRITE (7,*) ' +++ matrix after conversion from msrcsr +++'
   WRITE (7,*) '-----------------------------------------'
!
   CALL dump(1,n,.TRUE.,a,ja,ia,7)
!
   maxcol = 13
!
   CALL csrell(n,a,ja,ia,maxcol,a1,ja1,n,ndiag,ierr)
   iout = iout + 1
   IF ( ierr/=0 ) WRITE (iout,*) ' ***** ERROR IN CSRELL'
   IF ( ierr/=0 ) WRITE (iout,*) '     IERR = ' , ierr
   WRITE (iout,*) '-----------------------------------------'
   WRITE (iout,*) '  +++ matrix in ELLPACK-ITPACK format +++ '
   WRITE (iout,*) '-----------------------------------------'
   DO i = 1 , ndiag
      WRITE (iout,*) ' Column number: ' , i
      WRITE (iout,99003) (a1(n*(i-1)+k),k=1,n)
      WRITE (iout,99004) (ja1(n*(i-1)+k),k=1,n)
   ENDDO
   CALL ellcsr(n,a1,ja1,n,ndiag,a,ja,ia,NNZMAX,ierr)
   IF ( ierr/=0 ) WRITE (7,*) ' ***** ERROR IN ELLCSR'
   IF ( ierr/=0 ) WRITE (7,*) '     IERR = ' , ierr
   WRITE (7,*) '-----------------------------------------'
   WRITE (7,*) '  +++ matrix after conversion from ellcsr +++'
   WRITE (7,*) '-----------------------------------------'
   CALL dump(1,n,.TRUE.,a,ja,ia,7)
!
   CALL csrcsc(n,1,1,a,ja,ia,a1,ja1,ia1)
   WRITE (7,*) '-----------------------------------------'
   WRITE (7,*) '  +++ matrix after conversion from csrcsc  +++ '
   WRITE (7,*) '-----------------------------------------'
   CALL dump(1,n,.TRUE.,a1,ja1,ia1,7)
   CALL csrcsc(n,1,1,a1,ja1,ia1,a,ja,ia)
!
!--------test 1 : get main diagonal and subdiagonal
!     get some info on diagonals
   CALL infdia(n,ja,ia,iwk,idiag0)
   job = 0
   ioff(1) = 0
   ioff(2) = -1
   idiag = 2
   CALL csrdia(n,idiag,job,a,ja,ia,ndns,dns,ioff,a1,ja1,ia1,iwk)
   iout = iout + 1
   WRITE (iout,*) '-----------------------------------------'
   WRITE (iout,*) '  +++  diagonal format +++ '
   WRITE (iout,*) '-----------------------------------------'
   WRITE (iout,*) '  diagonals ioff = 0 and ioff = -1 '
   WRITE (iout,*) ' number of diag.s returned from csrdia=' , idiag
   DO kdiag = 1 , idiag
      WRITE (iout,*) ' diagonal offset = ' , ioff(kdiag)
      WRITE (iout,'(16f4.1)') (dns(k,kdiag),k=1,n)
   ENDDO
!     reverse conversion
   ndiag = ndns
   idiag = idiag0
   job = 10
   CALL csrdia(n,idiag,job,a,ja,ia,ndns,dns,ioff,a1,ja1,ia1,iwk)
   WRITE (iout,*) '-----------------------------------------'
   WRITE (iout,*) '  +++  second test diagonal format +++ '
   WRITE (iout,*) '         ** all diagonals of A  ** '
   WRITE (iout,*) '-----------------------------------------'
   WRITE (iout,*) ' number of diagonals on return from csrdia=' , idiag
   DO kdiag = 1 , idiag
      WRITE (iout,*) ' diagonal offset = ' , ioff(kdiag)
      WRITE (iout,'(16f4.1)') (dns(k,kdiag),k=1,n)
   ENDDO
!
!     reverse conversion
!
   job = 0
   CALL diacsr(n,job,idiag,dns,ndns,ioff,a,ja,ia)
!--------
   WRITE (7,*) '-----------------------------------------'
   WRITE (7,*) '  +++ matrix after conversion from diacsr  +++ '
   WRITE (7,*) '-----------------------------------------'
   CALL dump(1,n,.TRUE.,a,ja,ia,7)
!
!     checking the banded format
!
   lowd = 0
   job = 1
   CALL csrbnd(n,a,ja,ia,job,dns,ndns,lowd,ml,mu,ierr)
   iout = iout + 1
   IF ( ierr/=0 ) WRITE (iout,*) ' ***** ERROR IN CSRBND'
   IF ( ierr/=0 ) WRITE (iout,*) '     IERR = ' , ierr
   WRITE (iout,*) '-----------------------------------------'
   WRITE (iout,*) '       +++  banded  format +++ '
   WRITE (iout,*) ' bandwidth values found ml=' , ml , '  mu=' , mu
   WRITE (iout,*) '-----------------------------------------'
   WRITE (iout,'(4x,16i4)') (j,j=1,n)
   WRITE (iout,'(3x,65(1h-))')
   DO i = 1 , lowd
      WRITE (iout,99001) i , (dns(i,j),j=1,n)
   ENDDO
!
!     convert back to a, ja, ia format.
!
   len = NNZMAX
!--------
   CALL bndcsr(n,dns,ndns,lowd,ml,mu,a,ja,ia,len,ierr)
   WRITE (7,*) ' IERR IN BNDCSR = ' , ierr
   WRITE (7,*) '-----------------------------------------'
   WRITE (7,*) '  +++ matrix after conversion from bndcsr +++'
   WRITE (7,*) '-----------------------------------------'
   CALL dump(1,n,.TRUE.,a,ja,ia,7)
!
!     make sure it is sorted
!
   CALL csrcsc(n,1,1,a,ja,ia,a1,ja1,ia1)
!
!     checking skyline format.
!
   imod = 1
   CALL csrssk(n,imod,a1,ja1,ia1,a,ia,NNZMAX,ierr)
!
   IF ( ierr/=0 ) WRITE (iout,*) '     IERR = ' , ierr
   WRITE (iout,*) '-----------------------------------------'
   WRITE (iout,*) '    +++ Sym. Skyline format +++ '
   WRITE (iout,*) '-----------------------------------------'
   WRITE (iout,'(3x,65(1h-))')
!---------------------
!     create column values.
!---------------------
   DO i = 1 , n
      kend = ia(i+1) - 1
      kstart = ia(i)
      DO k = kstart , kend
         ja(k) = i - (kend-k)
      ENDDO
   ENDDO
!
   CALL dump(1,n,.TRUE.,a,ja,ia,iout)
!
!     back to ssr format..
!
   CALL sskssr(n,imod,a,ia,a1,ja1,ia1,NNZMAX,ierr)
   WRITE (7,*) '-----------------------------------------'
   WRITE (7,*) '  +++ matrix after conversion from sskcsr +++'
   WRITE (7,*) '-----------------------------------------'
   CALL dump(1,n,.TRUE.,a1,ja1,ia1,7)
!
!     checking jad format -----
!
!     first go back to the csr format ----
!
   CALL ssrcsr(3,1,n,a1,ja1,ia1,NNZMAX,a,ja,ia,iwk,iwk2,ierr)
!
   CALL csrjad(n,a,ja,ia,ndiag,iwk,a1,ja1,ia1)
!
   iout = iout + 1
   WRITE (iout,*) '-----------------------------------------'
   WRITE (iout,*) '   +++   matrix in JAD format +++ '
   WRITE (iout,*) '-----------------------------------------'
!
!     permutation array
!
   WRITE (iout,*) ' ** PERMUTATION ARRAY '
   WRITE (iout,'(17i4)') (iwk(k),k=1,n)
!     ------ diagonals
   DO i = 1 , ndiag
      WRITE (iout,*) ' J-diagonal number: ' , i
      k1 = ia1(i)
      k2 = ia1(i+1) - 1
      WRITE (iout,99003) (a1(k),k=k1,k2)
      WRITE (iout,99004) (ja1(k),k=k1,k2)
   ENDDO
!
!     back to csr format..
!
   CALL jadcsr(n,ndiag,a1,ja1,ia1,iwk,a,ja,ia)
!
   WRITE (7,*) '-----------------------------------------'
   WRITE (7,*) '  +++ matrix after conversion from jadcsr +++'
   WRITE (7,*) '-----------------------------------------'
   CALL dump(1,n,.TRUE.,a,ja,ia,7)
!-----------------------------------------------------------------------
!     checking the linked list format
!-----------------------------------------------------------------------
!
   nnz = ia(n+1) - ia(1)
   CALL csrlnk(n,ia,iwk)
!
!     print links in file 7 (no need for another file)
!
   iout = 7
   WRITE (iout,*) '-----------------------------------------'
   WRITE (iout,*) '   +++   matrix in LNK format +++ '
   WRITE (iout,*) '-----------------------------------------'
!
!     permutation array
!
   WRITE (iout,*) ' LINK ARRAY '
   WRITE (iout,*) ' ---------- '
   WRITE (iout,'(17i4)') (iwk(k),k=1,nnz)
!
!     back to csr format..
!
   CALL lnkcsr(n,a,ja,ia,iwk,a1,ja1,ia1)
!
   WRITE (7,*) '-----------------------------------------'
   WRITE (7,*) '  +++ matrix after conversion from lnkcsr +++'
   WRITE (7,*) '-----------------------------------------'
   CALL dump(1,n,.TRUE.,a,ja,ia,7)
99001 FORMAT (' ',i2,'|',16F4.1)
99003 FORMAT (' COEF  = ',16F4.0)
99004 FORMAT (' JCOEF = ',16I4)
!------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------end-of-chkfmt1----------------------------------------------
END PROGRAM chkfmt
