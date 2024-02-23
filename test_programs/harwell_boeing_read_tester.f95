!*==hb2pic.f90 processed by SPAG 8.04RA 23:57 22 Feb 2024
!!SPAG Open source Personal, Educational or Academic User Clemson University  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
 
PROGRAM hb2pic
!------------------------------------------------------------------c
!
! reads a harwell-Boeing matrix and creates a pic file for pattern.
!
!------------------------------------------------------------------c
   USE ISO_FORTRAN_ENV                 
   IMPLICIT NONE
!
! PARAMETER definitions rewritten by SPAG
!
   INTEGER , PARAMETER :: NMAX = 5000 , NZMAX = 70000
!
! Local variable declarations rewritten by SPAG
!
   REAL(REAL64) , DIMENSION(NZMAX) :: a
   CHARACTER(2) :: guesol
   INTEGER , DIMENSION(NMAX+1) :: ia
   INTEGER :: ierr , iounit , job , mode , ncol , nnz , nrhs , nrow
   INTEGER , SAVE :: iin , iout
   INTEGER , DIMENSION(NZMAX) :: ja
   CHARACTER(8) :: key
   REAL(REAL64) , DIMENSION(1) :: rhs
   CHARACTER(72) :: title
   REAL(REAL64) :: type
   LOGICAL :: valued
   EXTERNAL pltmt , readmt
!
! End of declarations rewritten by SPAG
!
!--------------
   DATA iin/5/ , iout/6/
!--------------
   job = 2
   nrhs = 0
   CALL readmt(NMAX,NZMAX,job,iin,a,ja,ia,rhs,nrhs,guesol,nrow,ncol,nnz,title,key,type,ierr)
!---- if not readable return
   IF ( ierr/=0 ) THEN
      WRITE (iout,99001) ierr
99001 FORMAT (' **ERROR: Unable to read matrix',/,' Message returned fom readmt was ierr =',i3)
      STOP
   ENDIF
   valued = (job>=2)
!-------
   mode = 1
   iounit = 6
   job = 11
   CALL pltmt(nrow,ncol,mode,ja,ia,title,key,type,job,iout)
!-----------------------------------------------------------------------
END PROGRAM hb2pic
