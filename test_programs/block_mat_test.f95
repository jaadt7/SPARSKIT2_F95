!*==bfivept.f90 processed by SPAG 8.04RA 10:34 23 Feb 2024
!!SPAG Open source Personal, Educational or Academic User Clemson University  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
 
PROGRAM bfivept
   USE ISO_FORTRAN_ENV                 
   IMPLICIT NONE
!
! PARAMETER definitions rewritten by SPAG
!
   INTEGER , PARAMETER :: NXMAX = 20 , NMX = NXMAX*NXMAX*NXMAX , NTOT = NMX*25
!
! Local variable declarations rewritten by SPAG
!
   REAL(REAL64) , DIMENSION(NTOT) :: a , ao
   CHARACTER(2) :: guesol
   INTEGER , DIMENSION(NTOT) :: ia , iao , iau , ja , jao
   INTEGER :: ifmt , iout , job , n , na , nfree , nx , ny , nz
   CHARACTER(8) :: key
   CHARACTER(50) :: matfile
   REAL :: rhs
   REAL(REAL64) , DIMENSION(7,100) :: stencil
   CHARACTER(72) :: title
   CHARACTER(3) :: type
   EXTERNAL bsrcsr , gen57bl , prtmt
!
! End of declarations rewritten by SPAG
!
!-----------------------------------------------------------------------
! main program for generating BLOCK 5 point and 7-point matrices in the
! Harwell-Boeing format.  Creates a file with containing a
! harwell-boeing matrix.
!
! max block size = 5
! max number of grid points = 8000  = ( nx * ny * nz .le. 8000)
! matrix dimension =  (nx*ny*nz* Block-size**2) .le. 8000 * 25= 200,000
!
! typical session:
! Enter nx, ny, nz : 10 10 1
! enter block-size : 4
! enter filename for matrix: test.mat
! output matrix in data file : test.mat
!
! nz =1 will create a 2-D problem
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
   WRITE (6,*) '  '
   WRITE (6,'(22hEnter  nx, ny, nz   : ,$)')
   READ (5,*) nx , ny , nz
   WRITE (6,'(22hnfree (Block size)  : ,$)')
   READ (5,*) nfree
 
   WRITE (6,'(22hFilename for matrix : ,$)')
 
   READ (5,'(a50)') matfile
   OPEN (UNIT=7,FILE=matfile)
!
   WRITE (6,*) ' output in data file : ' , matfile
 
!------------------------------------------------------
   na = nfree*nfree
!
   CALL gen57bl(nx,ny,nz,nfree,na,n,a,ja,ia,iau,stencil)
!------------------------------------------------------
 
   PRINT * , ' n=' , n , ' nfree ' , nfree , ' na =' , na
 
   CALL bsrcsr(1,n,nfree,na,a,ja,ia,ao,jao,iao)
   n = n*nfree      ! Apr. 21, 1995
 
   guesol = 'NN'
 
   title = ' BLOCK 5-POINT TEST MATRIX FROM SPARSKIT               '
   type = 'RUA'
   key = 'BLOCK5PT'
!              12345678
   ifmt = 15
   job = 2
   iout = 7
   CALL prtmt(n,n,ao,jao,iao,rhs,guesol,title,key,type,ifmt,job,iout)
   PRINT * , ' output in data file : ' , matfile
!
END PROGRAM bfivept
