!*==zlatev.f90 processed by SPAG 8.04RA 11:58 23 Feb 2024
!!SPAG Open source Personal, Educational or Academic User Clemson University  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
 
PROGRAM zlatev
!-----------------------------------------------------------------------
!
! test suite for zlatev matrices. generates three matrices and
! writes them in three different files in Harwell-Boeing format.
!      zlatev1.mat  produced from matrf2
!      zlatev2.mat  produced from  dcn
!      zlatev3.mat  produced from  ecn
!
!-----------------------------------------------------------------------
   USE ISO_FORTRAN_ENV                 
   IMPLICIT NONE
!
! PARAMETER definitions rewritten by SPAG
!
   INTEGER , PARAMETER :: NMAX = 1000 , NZMAX = 20*NMAX
!
! Local variable declarations rewritten by SPAG
!
   REAL(REAL64) , DIMENSION(NZMAX) :: a
   REAL(REAL64) :: alpha , rhs
   CHARACTER(2) :: guesol
   INTEGER , DIMENSION(NZMAX) :: ia , ja
   INTEGER :: ic , ierr , ifmt , index , iout , job , m , n , ne , nn , nz
   INTEGER , DIMENSION(NMAX) :: iwk
   CHARACTER(3) :: key
   CHARACTER(72) :: title
   CHARACTER(8) :: type
   EXTERNAL coicsr , dcn , ecn , matrf2 , prtmt
!
! End of declarations rewritten by SPAG
!
!
   OPEN (UNIT=7,FILE='zlatev1.mat')
   OPEN (UNIT=8,FILE='zlatev2.mat')
   OPEN (UNIT=9,FILE='zlatev3.mat')
!
   m = 100
   n = m
   ic = n/2
   index = 10
   alpha = 5.0
   nn = NZMAX
!
! call matrf2
!
   CALL matrf2(m,n,ic,index,alpha,nn,nz,a,ia,ja,ierr)
   job = 1
!      do 110 i = 1, nz
!         print *, ia(i), ja(i), a(i)
! 110  continue
   CALL coicsr(n,nz,job,a,ja,ia,iwk)
!-----
   title = ' 1st matrix from zlatev examples                '
   type = 'RUA'
   key = ' ZLATEV1'
   iout = 7
   guesol = 'NN'
!
   ifmt = 3
   job = 2
!
! write result in H-B format.
!
! Replaces prtmt with smms in order to print matrix in format for
! SMMS instead.
!      call smms (n,1,n,0,a,ja,ia,iout)
   CALL prtmt(n,n,a,ja,ia,rhs,guesol,title,type,key,ifmt,job,iout)
 
!-------- second type of matrices dcn matrices ---------------
   n = 200
   nn = NZMAX
   ic = 20
!-------------------------------------------------------
! matrix of the type e(c,n)
!-------------------------------------------------------
   CALL dcn(a,ia,ja,n,ne,ic,nn,ierr)
!---------------------------------------------------
   CALL coicsr(n,ne,job,a,ja,ia,iwk)
   title = ' 2nd matrix from zlatev examples                '
   iout = iout + 1
   guesol = 'NN'
   type = 'RUA'
   key = ' ZLATEV2'
!
   ifmt = 3
   job = 2
!
! write result in second file
!
! Replaced prtmt with smms in order to print matrix in format for
! SMMS instead.
!      call smms (n,1,n,0,a,ja,ia,iout)
   CALL prtmt(n,n,a,ja,ia,rhs,guesol,title,type,key,ifmt,job,iout)
!-------------------------------------------------------
! matrix of the type e(c,n)
!-------------------------------------------------------
   n = 200
   ic = 20
   nn = NZMAX
!
! call ecn
!
   CALL ecn(n,ic,ne,ia,ja,a,nn,ierr)
   CALL coicsr(n,ne,job,a,ja,ia,iwk)
   title = ' 3nd matrix from zlatev examples                '
   guesol = 'NN'
   type = 'RUA'
   key = ' ZLATEV3'
   iout = iout + 1
!
   ifmt = 3
   job = 2
!
! write resulting matrix in third file
!
! Replaced prtmt with smms in order to print matrix in format for
! SMMS instead.
!      call smms (n,1,n,0,a,ja,ia,iout)
   CALL prtmt(n,n,a,ja,ia,rhs,guesol,title,type,key,ifmt,job,iout)
END PROGRAM zlatev
