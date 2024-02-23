!*==chkio.f90 processed by SPAG 8.04RA 23:57 22 Feb 2024
!!SPAG Open source Personal, Educational or Academic User Clemson University  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
PROGRAM chkio
!------------------------------------------------------------------c
! test suite for Part I : I/O routines.                            c
! tests the following : gen5pt.f, prtmt, readmt, amd pltmt.        c
! 1) generates a 100 x 100 5pt matrix,                             c
! 2) prints it with a given format in file 'first.mat'             c
! 3) reads the matrix from 'first.mat' using readmat               c
! 4) prints it again in file 'second.mat' in  a different format   c
! 5) makes 4 pic files to show the different options of pltmt.     c
!    these are in job0.pic, job01.pic, job10.pic, job11.pic        c
!                          coded by Y. Saad, RIACS, 08/31/1989.    c
!------------------------------------------------------------------c
   USE ISO_FORTRAN_ENV                 
   IMPLICIT NONE
!
! PARAMETER definitions rewritten by SPAG
!
   INTEGER , PARAMETER :: NXMAX = 20 , NMX = NXMAX*NXMAX
!
! Local variable declarations rewritten by SPAG
!
   REAL(REAL64) , DIMENSION(7*NMX) :: a
   REAL(REAL64) , DIMENSION(6) :: al
   CHARACTER(2) :: guesol
   INTEGER :: i , ierr , ifmt , iout , j , job , k , mode , n , ncol , nmax , nnz , nrhs , nrow , nx , ny , nz , nzmax
   INTEGER , DIMENSION(NMX) :: ia , iau
   INTEGER , DIMENSION(7*NMX) :: ja
   CHARACTER(8) :: key
   REAL(REAL64) , DIMENSION(3*NMX) :: rhs
   CHARACTER(72) :: title
   CHARACTER(3) :: type
   EXTERNAL gen57pt , pltmt , prtmt , readmt
!
! End of declarations rewritten by SPAG
!
!----- open statements ----------------
   OPEN (UNIT=7,FILE='first.mat')
   OPEN (UNIT=8,FILE='second.mat')
   OPEN (UNIT=20,FILE='job00.pic')
   OPEN (UNIT=21,FILE='job01.pic')
   OPEN (UNIT=22,FILE='job10.pic')
   OPEN (UNIT=23,FILE='job11.pic')
!
!---- dimension of grid
!
   nx = 10
   ny = 10
   nz = 1
   al(1) = 1.0D0
   al(2) = 0.0D0
   al(3) = 2.3D1
   al(4) = 0.4D0
   al(5) = 0.0D0
   al(6) = 8.2D-2
!
!---- generate grid problem.
!
   CALL gen57pt(nx,ny,nz,al,0,n,a,ja,ia,iau,rhs)
!
!---- create the Harwell-Boeing matrix. Start by defining title,
!     and type. them define format and print it.
!
   WRITE (title,99001) nx , ny
99001 FORMAT ('Five-point matrix on a square region',' using a ',I2,' by ',I2,' grid *SPARSKIT*')
   key = 'Fivept10'
   type = 'RSA'
   ifmt = 5
   job = 3
   guesol = 'GX'
!
! define a right hand side of ones, an initial guess of two's
! and an exact solution of three's.
!
   DO k = 1 , 3*n
      rhs(k) = real(1+(k-1)/n)
   ENDDO
!
   CALL prtmt(n,n,a,ja,ia,rhs,guesol,title,key,type,ifmt,job,7)
!---- read it again in same matrix a, ja, ia
   nmax = NMX
   nzmax = 7*NMX
   DO k = 1 , 3*n
      rhs(k) = 0.0
   ENDDO
   job = 3
!
   REWIND 7
!
   nrhs = 3*n
!
   CALL readmt(nmax,nzmax,job,7,a,ja,ia,rhs,nrhs,guesol,nrow,ncol,nnz,title,key,type,ierr)
   PRINT * , ' ierr = ' , ierr , ' nrhs ' , nrhs
!
! matrix read.  print it again in a different format
!
   ifmt = 102
   ncol = nrow
   job = 3
!
   CALL prtmt(nrow,ncol,a,ja,ia,rhs,guesol,title,key,type,ifmt,job,8)
!
!---- print four pic files
!
   mode = 0
   DO i = 1 , 2
      DO j = 1 , 2
         job = (i-1)*10 + j - 1
         iout = 20 + (i-1)*2 + j - 1
         CALL pltmt(nrow,ncol,mode,ja,ia,title,key,type,job,iout)
      ENDDO
   ENDDO
!--------
END PROGRAM chkio
