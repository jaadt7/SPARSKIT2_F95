!*==hb2ps.f90 processed by SPAG 8.04RA 23:57 22 Feb 2024
!!SPAG Open source Personal, Educational or Academic User Clemson University  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
PROGRAM hb2ps
   USE ISO_FORTRAN_ENV                 
   IMPLICIT NONE
!
! PARAMETER definitions rewritten by SPAG
!
   INTEGER , PARAMETER :: NMAX = 10000 , NZMAX = 100000
!
! Local variable declarations rewritten by SPAG
!
   REAL(REAL64) , DIMENSION(1) :: a , rhs
   CHARACTER(2) :: guesol
   INTEGER , DIMENSION(NMAX+1) :: ia
   INTEGER , DIMENSION(1) :: idummy
   INTEGER :: ierr , job , ncol , nnz , nrhs , nrow
   INTEGER , SAVE :: iin , iout , mode , nlines , ptitle
   INTEGER , DIMENSION(NZMAX) :: ja
   CHARACTER(8) :: key
   CHARACTER(2) , SAVE :: munt
   REAL , SAVE :: size
   CHARACTER(72) :: title
   REAL :: type
   EXTERNAL pspltm , readmt
!
! End of declarations rewritten by SPAG
!
!----------------------------------------------------------------------
! translates a harwell - boeing file into a post-script file. Usage:
!                   hb2ps < HB_file > Postscript_file
! where hb2ps is the executable generated from this program,
! HB_file is a file containing a matrix stored in Harwell-Boeing
! format and Postscript_file is a file to contain the post-script file.
!----------------------------------------------------------------------
   DATA iin/5/ , iout/6/ , size/5.0/ , nlines/0/ , ptitle/0/ , mode/0/
   DATA munt/'in'/
!-----------------------------------------------------------------------
   job = 1
   nrhs = 0
!
!     read matrix in Harwell-Boeing format
!
   CALL readmt(NMAX,NZMAX,job,iin,a,ja,ia,rhs,nrhs,guesol,nrow,ncol,nnz,title,key,type,ierr)
!
!     if   not readable return
!
   IF ( ierr/=0 ) THEN
      WRITE (iout,99001) ierr
!
99001 FORMAT (' **ERROR: Unable to read matrix',/,' Message returned fom readmt was ierr =',i3)
      STOP
   ENDIF
!
!     call post script generator
!
   CALL pspltm(nrow,ncol,mode,ja,ia,title,ptitle,size,munt,nlines,idummy,iout)
!-----------------------------------------------------------------------
END PROGRAM hb2ps
