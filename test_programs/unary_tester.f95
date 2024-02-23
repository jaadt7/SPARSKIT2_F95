!*==chkfmt.f90 processed by SPAG 8.04RA 23:02 22 Feb 2024
!!SPAG Open source Personal, Educational or Academic User Clemson University  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
PROGRAM chkfmt
!----------------------------------------------------------------------c
!                          S P A R S K I T                             c
!----------------------------------------------------------------------c
! test suite for the unary routines.                                   c
! tests some of the routines in the module unary. Still needs to tests c
! many other routines.                                                 c
! Last update: May 2, 1994.
!----------------------------------------------------------------------c
   USE ISO_FORTRAN_ENV                 
   IMPLICIT NONE
!
! PARAMETER definitions rewritten by SPAG
!
   INTEGER , PARAMETER :: NXMAX = 10 , NMX = NXMAX*NXMAX , NNZMAX = 10*NMX
!
! Local variable declarations rewritten by SPAG
!
   REAL(REAL64) , DIMENSION(NNZMAX) :: a , a1 , rhs
   REAL(REAL64) , DIMENSION(6) :: al
   REAL(REAL64) , DIMENSION(20,20) :: dns
   INTEGER :: i , ierr , iout , job , n , nnz , nx , ny , nz
   INTEGER , DIMENSION(NMX+1) :: ia
   INTEGER , DIMENSION(NNZMAX) :: ia1 , ja , ja1
   INTEGER , DIMENSION(NMX*2+1) :: iwk
   INTEGER , SAVE :: ndns
   INTEGER , DIMENSION(16) :: perm
   INTEGER , DIMENSION(16) , SAVE :: qperm
   EXTERNAL csort , csrdns , dmpdns , dperm , dump , gen57pt
!
! End of declarations rewritten by SPAG
!
   DATA ndns/20/
   DATA qperm/1 , 3 , 6 , 8 , 9 , 11 , 14 , 16 , 2 , 4 , 5 , 7 , 10 , 12 , 13 , 15/
!
! define correct permutation
!
   DO i = 1 , 16
      perm(qperm(i)) = i
   ENDDO
!----- open statements ----------------
   OPEN (UNIT=7,FILE='unary.mat')
!
!---- dimension of grid
!
!     generate a 16 by 16 matrix (after eliminating boundary)
   nx = 6
   ny = 6
   nz = 1
   al(1) = 0.D0
   al(2) = 0.D0
   al(3) = 0.D0
   al(4) = 0.D0
   al(5) = 0.D0
   al(6) = 0.D0
!
!---- generate grid problem.
!
   CALL gen57pt(nx,ny,nz,al,0,n,a,ja,ia,iwk,rhs)
!
!---- write out the matrix
!
   iout = 7
   nnz = ia(n+1) - 1
   WRITE (iout,*) '-----------------------------------------'
   WRITE (iout,*) '  +++  initial matrix in CSR format +++ '
   WRITE (iout,*) '-----------------------------------------'
   CALL dump(1,n,.TRUE.,a,ja,ia,iout)
!
! call csrdns
!
   CALL csrdns(n,n,a,ja,ia,dns,ndns,ierr)
!
! write it out as a dense matrix.
!
   WRITE (iout,*) '-----------------------------------------'
   WRITE (iout,*) '  +++ initial matrix in DENSE format+++ '
   WRITE (iout,*) '-----------------------------------------'
   CALL dmpdns(n,n,ndns,dns,iout)
!
! red black ordering
!
   job = 1
!
   CALL dperm(n,a,ja,ia,a1,ja1,ia1,perm,perm,job)
!
   nnz = ia(n+1) - 1
   WRITE (iout,*) '-----------------------------------------'
   WRITE (iout,*) '  +++ red-black matrix in CSR format +++ '
   WRITE (iout,*) '-----------------------------------------'
   CALL dump(1,n,.TRUE.,a1,ja1,ia1,iout)
!
! sort matrix
!
   CALL csort(n,a1,ja1,ia1,.TRUE.)
   nnz = ia(n+1) - 1
   WRITE (iout,*) '-----------------------------------------'
   WRITE (iout,*) '  +++     matrix after sorting    +++ '
   WRITE (iout,*) '-----------------------------------------'
   CALL dump(1,n,.TRUE.,a1,ja1,ia1,iout)
!
!
! convert into dense format
!
   CALL csrdns(n,n,a1,ja1,ia1,dns,ndns,ierr)
   WRITE (iout,*) '-----------------------------------------'
   WRITE (iout,*) '  +++ red-black matrix in DENSE format+++ '
   WRITE (iout,*) '-----------------------------------------'
   CALL dmpdns(n,n,ndns,dns,iout)
END PROGRAM chkfmt
!*==dmpdns.f90 processed by SPAG 8.04RA 23:02 22 Feb 2024
!!SPAG Open source Personal, Educational or Academic User Clemson University  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
!-----------------------------------------------------------------------
SUBROUTINE dmpdns(Nrow,Ncol,Ndns,Dns,Iout)
   USE ISO_FORTRAN_ENV                 
   IMPLICIT NONE
!
! Dummy argument declarations rewritten by SPAG
!
   INTEGER , INTENT(IN) :: Ndns
   INTEGER , INTENT(IN) :: Nrow
   INTEGER , INTENT(IN) :: Ncol
   REAL(REAL64) , INTENT(IN) , DIMENSION(Ndns,*) :: Dns
   INTEGER , INTENT(IN) :: Iout
!
! Local variable declarations rewritten by SPAG
!
   CHARACTER(80) :: fmt
   INTEGER :: i , j , j1 , j2 , last
!
! End of declarations rewritten by SPAG
!
!-----------------------------------------------------------------------
! this subroutine prints out a dense matrix in a simple format.
! the zero elements of the matrix are omitted. The format for the
! nonzero elements is f4.1, i.e., very little precision is provided.
!-----------------------------------------------------------------------
! on entry
! --------
! nrow = row dimension of matrix
! ncol = column dimension of matrix
! ndns = first dimension of array dns.
! dns  = double dimensional array of size n x n containing the matrix
! iout = logical unit where to write matrix
!
! on return
! ---------
! matrix will be printed out on unit output iout.
!------------------------------------------------------------------------
!         local variables
!
! prints out a dense matrix -- without the zeros.
!
   WRITE (Iout,'(4x,16i4)') (j,j=1,Ncol)
   fmt(1:5) = '    |'
   j1 = 6
   DO j = 1 , Ncol
      j2 = j1 + 4
      fmt(j1:j2) = '----'
      j1 = j2
   ENDDO
   last = j1
   fmt(last:last) = '|'
   WRITE (Iout,*) fmt
!
! undo loop 1 ---
!
   j1 = 6
   DO j = 1 , Ncol
      j2 = j1 + 4
      fmt(j1:j2) = '   '
      j1 = j2
   ENDDO
!
   DO i = 1 , Nrow
      j1 = 6
      WRITE (fmt,99001) i
99001 FORMAT (' ',i2,' |')
      DO j = 1 , Ncol
         j2 = j1 + 4
         IF ( Dns(i,j)/=0.0 ) THEN
            WRITE (fmt(j1:j2),99002) Dns(i,j)
99002       FORMAT (f4.1)
         ENDIF
         j1 = j2
      ENDDO
      fmt(last:last) = '|'
      WRITE (Iout,*) fmt
   ENDDO
   fmt(1:5) = '    |'
   j1 = 6
   DO j = 1 , Ncol
      j2 = j1 + 4
      fmt(j1:j2) = '----'
      j1 = j2
   ENDDO
   fmt(last:last) = '|'
   WRITE (Iout,*) fmt
END SUBROUTINE dmpdns
