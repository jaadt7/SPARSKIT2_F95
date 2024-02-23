!*==riters.f90 processed by SPAG 8.04RA 00:36 23 Feb 2024
!!SPAG Open source Personal, Educational or Academic User Clemson University  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
PROGRAM riters
   USE ISO_FORTRAN_ENV                 
   IMPLICIT NONE
!
! PARAMETER definitions rewritten by SPAG
!
   INTEGER , PARAMETER :: NMAX = 5000 , NZMAX = 100000 , MAXITS = 60 , LWK = NMAX*40
!
! Local variable declarations rewritten by SPAG
!
   REAL(REAL64) , DIMENSION(NZMAX) :: a , au
   REAL(REAL64) , DIMENSION(16) :: fpar
   CHARACTER(2) :: guesol
   INTEGER :: i , ierr , iounit , job , lfil , ncol , nnz , nrhs , nrow , nwk
   INTEGER , DIMENSION(NMAX) :: ia
   INTEGER , DIMENSION(16) :: ipar
   INTEGER , DIMENSION(NMAX*3) :: iw
   INTEGER , DIMENSION(NZMAX) :: ja , jau , ju
   CHARACTER(8) :: key
   REAL(REAL64) , DIMENSION(NMAX) :: rhs , sol , xran
   CHARACTER(72) :: title
   REAL(REAL64) :: tol
   CHARACTER(3) :: type
   REAL(REAL64) , DIMENSION(NMAX*40) :: wk
   EXTERNAL amux , bcg , bcgstab , cg , cgnr , dbcg , dqgmres , fgmres , fom , gmres , ilut , readmt , runrc , tfqmr
!
! End of declarations rewritten by SPAG
!
!-----------------------------------------------------------------------
! test program for iters -- the basic iterative solvers
!
!     this program reads a Harwell/Boeing matrix from standard input
!     and solves the linear system with an artifical right-hand side
!     (the solution is a vector of (1,1,...,1)^T)
!-----------------------------------------------------------------------
!      implicit none
!      implicit real*8 (a-h,o-z)
!
!     set the parameters for the iterative solvers
!
   ipar(2) = 2
   ipar(3) = 1
   ipar(4) = LWK
   ipar(5) = 16
   ipar(6) = MAXITS
   fpar(1) = 1.0D-5
   fpar(2) = 1.0D-10
!--------------------------------------------------------------
!     read in a matrix from standard input
!--------------------------------------------------------------
   iounit = 5
   job = 2
   nrhs = 0
   CALL readmt(NMAX,NZMAX,job,iounit,a,ja,ia,a,nrhs,guesol,nrow,ncol,nnz,title,key,type,ierr)
   PRINT * , 'READ the matrix ' , key , type
   PRINT * , title
   PRINT *
!
!     set-up the preconditioner ILUT(15, 1E-4) ! new definition of lfil
!
   lfil = 15
   tol = 1.0D-4    ! this is too high for ilut for saylr1
   tol = 1.0D-7
   nwk = NZMAX
   CALL ilut(nrow,a,ja,ia,lfil,tol,au,jau,ju,nwk,wk,iw,ierr)
   ipar(2) = 2
!
!     generate a linear system with known solution
!
   DO i = 1 , nrow
      sol(i) = 1.0D0
      xran(i) = 0.D0
   ENDDO
   CALL amux(nrow,sol,rhs,a,ja,ia)
   PRINT * , ' '
   PRINT * , '	*** CG ***'
   CALL runrc(nrow,rhs,sol,ipar,fpar,wk,xran,a,ja,ia,au,jau,ju,cg)
   PRINT * , ' '
   PRINT * , '	*** BCG ***'
   CALL runrc(nrow,rhs,sol,ipar,fpar,wk,xran,a,ja,ia,au,jau,ju,bcg)
   PRINT * , ' '
   PRINT * , '	*** DBCG ***'
   CALL runrc(nrow,rhs,sol,ipar,fpar,wk,xran,a,ja,ia,au,jau,ju,dbcg)
   PRINT * , ' '
   PRINT * , '	*** CGNR ***'
   CALL runrc(nrow,rhs,sol,ipar,fpar,wk,xran,a,ja,ia,au,jau,ju,cgnr)
   PRINT * , ' '
   PRINT * , '	*** BCGSTAB ***'
   CALL runrc(nrow,rhs,sol,ipar,fpar,wk,xran,a,ja,ia,au,jau,ju,bcgstab)
   PRINT * , ' '
   PRINT * , '	*** TFQMR ***'
   CALL runrc(nrow,rhs,sol,ipar,fpar,wk,xran,a,ja,ia,au,jau,ju,tfqmr)
   PRINT * , ' '
   PRINT * , '	*** FOM ***'
   CALL runrc(nrow,rhs,sol,ipar,fpar,wk,xran,a,ja,ia,au,jau,ju,fom)
   PRINT * , ' '
   PRINT * , '	*** GMRES ***'
   CALL runrc(nrow,rhs,sol,ipar,fpar,wk,xran,a,ja,ia,au,jau,ju,gmres)
   PRINT * , ' '
   PRINT * , '	*** FGMRES ***'
   CALL runrc(nrow,rhs,sol,ipar,fpar,wk,xran,a,ja,ia,au,jau,ju,fgmres)
   PRINT * , ' '
   PRINT * , '	*** DQGMRES ***'
   CALL runrc(nrow,rhs,sol,ipar,fpar,wk,xran,a,ja,ia,au,jau,ju,dqgmres)
END PROGRAM riters
!-----end-of-main
!-----------------------------------------------------------------------
