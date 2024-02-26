!*==riters.f90 processed by SPAG 8.04RA 00:18 23 Feb 2024
!!SPAG Open source Personal, Educational or Academic User Clemson University  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
 
PROGRAM riters
   USE ISO_FORTRAN_ENV                 
   IMPLICIT NONE
!
! PARAMETER definitions rewritten by SPAG
!
   INTEGER , PARAMETER :: NMAX = 5000 , NZMAX = 100000 , MAXITS = 60 , LWK = NMAX*40
!
! COMMON variable declarations rewritten by SPAG
!
   REAL(REAL64) :: Alpha , Gammax , Gammay
   COMMON /func  / Gammax , Gammay , Alpha
!
! Local variable declarations rewritten by SPAG
!
   REAL(REAL64) , DIMENSION(NZMAX) :: a , au
   REAL(REAL64) , DIMENSION(NMAX) :: al , rhs , sol , xran
   REAL(REAL64) , DIMENSION(16) :: fpar
   INTEGER :: i , ierr , lfil , nrow , nwk , nx , ny , nz
   INTEGER , DIMENSION(NMAX) :: ia
   INTEGER , DIMENSION(16) :: ipar
   INTEGER , DIMENSION(NMAX*3) :: iw
   INTEGER , DIMENSION(NZMAX) :: ja , jau , ju
   REAL(REAL64) :: tol
   REAL(REAL64) , DIMENSION(NMAX*40) :: wk
   EXTERNAL amux , bcg , bcgstab , cg , cgnr , dbcg , dqgmres , fgmres , fom , gen57pt , gmres , ilut , runrc , tfqmr
!
! End of declarations rewritten by SPAG
!
!-----------------------------------------------------------------------
! test program for iters -- the basic iterative solvers
!
!     this program generates a sparse matrix using
!     GEN57PT and then solves a linear system with an
!     artificial rhs (the solution is a vector of (1,1,...,1)^T).
!-----------------------------------------------------------------------
!      implicit none
!      implicit real*8 (a-h,o-z)
!
!-----------------------------------------------------------------------
! pde to be discretized is :
!---------------------------
!
! -Lap u + gammax exp (xy)delx u + gammay exp (-xy) dely u +alpha u
!
! where Lap = 2-D laplacean, delx = part. der. wrt x,
! dely = part. der. wrt y.
! gammax, gammay, and alpha are passed via the commun func.
!
!-----------------------------------------------------------------------
!
! data for PDE:
!
   nx = 6
   ny = 6
   nz = 1
   Alpha = 0.0
   Gammax = 0.0
   Gammay = 0.0
!
!     set the parameters for the iterative solvers
!
   ipar(2) = 2
   ipar(3) = 1
   ipar(4) = LWK
   ipar(5) = 10
   ipar(6) = MAXITS
   fpar(1) = 1.0D-5
   fpar(2) = 1.0D-10
!--------------------------------------------------------------
! call GEN57PT to generate matrix in compressed sparse row format
!
!     al(1:6) are used to store part of the boundary conditions
!     (see documentation on GEN57PT.)
!--------------------------------------------------------------
   al(1) = 0.0
   al(2) = 0.0
   al(3) = 0.0
   al(4) = 0.0
   al(5) = 0.0
   al(6) = 0.0
   nrow = nx*ny*nz
   CALL gen57pt(nx,ny,nz,al,0,nrow,a,ja,ia,ju,rhs)
   PRINT * , 'RITERS: generated a finite difference matrix'
   PRINT * , '        grid size = ' , nx , ' X ' , ny , ' X ' , nz
   PRINT * , '        matrix size = ' , nrow
!
!     set-up the preconditioner ILUT(15, 1E-4)  ! new definition of lfil
!
   lfil = 3
   tol = 1.0D-4
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
