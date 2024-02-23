!*==fivept.f90 processed by SPAG 8.04RA 10:34 23 Feb 2024
!!SPAG Open source Personal, Educational or Academic User Clemson University  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
PROGRAM fivept
   USE ISO_FORTRAN_ENV                 
   IMPLICIT NONE
!
! PARAMETER definitions rewritten by SPAG
!
   INTEGER , PARAMETER :: NXMAX = 50 , NMX = NXMAX*NXMAX
!
! Local variable declarations rewritten by SPAG
!
   REAL(REAL64) , DIMENSION(7*NMX) :: a
   REAL(REAL64) , DIMENSION(6) :: al
   CHARACTER(2) :: guesol
   INTEGER , DIMENSION(NMX) :: ia , iau
   INTEGER :: ifmt , iout , job , n , nx , ny , nz
   INTEGER , DIMENSION(7*NMX) :: ja
   CHARACTER(8) :: key
   CHARACTER(50) :: matfile
   REAL(REAL64) , DIMENSION(NMX) :: rhs
   CHARACTER(72) :: title
   CHARACTER(3) :: type
   EXTERNAL gen57pt , prtmt
!
! End of declarations rewritten by SPAG
!
!-----------------------------------------------------------------------
! main program for generating 5 point and 7-point matrices in the
! Harwell-Boeing format.  Creates a file with containing a
! harwell-boeing matrix. typical session:
! user answer are after the colon
! Enter nx, ny, nz  : 10 10 1
! Filename for matrix: test.mat
! output matrix in data file : test.mat
!
! nz = 1 will create a 2-D problem
!
!-----------------------------------------------------------------------
!      implicit none
!-----------------------------------------------------------------------
   WRITE (6,*) '  '
   WRITE (6,'(22hEnter  nx, ny, nz   : ,$)')
   READ (5,*) nx , ny , nz
   WRITE (6,'(22hFilename for matrix : ,$)')
   READ (5,'(a50)') matfile
   OPEN (UNIT=7,FILE=matfile)
!
!     boundary condition is partly specified here
!
!      al(1) = 1.0D0
!      al(2) = 0.0D0
!      al(3) = 2.3D1
!      al(4) = 0.4D0
!      al(5) = 0.0D0
!      al(6) = 8.2D-2
   al(1) = 0.0D0
   al(2) = 0.0D0
   al(3) = 0.0D1
   al(4) = 0.0D0
   al(5) = 0.0D0
   al(6) = 0.0D0
!
   CALL gen57pt(nx,ny,nz,al,0,n,a,ja,ia,iau,rhs)
   iout = 7
!
!     write out the matrix
!
   guesol = 'NN'
   title = ' 5-POINT TEST MATRIX FROM SPARSKIT                    '
!          '123456789012345678901234567890123456789012345678901234567890
   type = 'RUA'
   key = 'SC5POINT'
!           12345678
   ifmt = 15
   job = 2
!   upper part only??
!      call getu (n, a, ja, ia, a, ja, ia)
   CALL prtmt(n,n,a,ja,ia,rhs,guesol,title,key,type,ifmt,job,iout)
   WRITE (6,*) ' output matrix in data file : ' , matfile
!
END PROGRAM fivept
