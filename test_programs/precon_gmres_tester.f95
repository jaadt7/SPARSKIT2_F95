!*==rilut.f90 processed by SPAG 8.04RA 00:32 23 Feb 2024
!!SPAG Open source Personal, Educational or Academic User Clemson University  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
PROGRAM rilut
!-----------------------------------------------------------------------
!     test program for ilut preconditioned gmres.
!     this program generates a sparse matrix using
!     matgen and then solves a linear system with an
!     artificial rhs.
!-----------------------------------------------------------------------
   USE ISO_FORTRAN_ENV                 
   IMPLICIT NONE
!
! PARAMETER definitions rewritten by SPAG
!
   INTEGER , PARAMETER :: NMAX = 5000 , NZMAX = 100000
!
! COMMON variable declarations rewritten by SPAG
!
   REAL(REAL64) :: Alpha , Gammax , Gammay
   COMMON /func  / Gammax , Gammay , Alpha
!
! Local variable declarations rewritten by SPAG
!
   REAL(REAL64) , DIMENSION(NZMAX) :: a , au
   REAL(REAL64) , DIMENSION(NMAX) :: al , rhs , x , xran , y
   REAL(REAL64) :: alph , eps , permtol , tol
   REAL(REAL64) , DIMENSION(16) :: fpar
   INTEGER , DIMENSION(NMAX) :: ia
   INTEGER :: ierr , im , iout , j , k , lfil , maxits , meth , n , nwk , nx , ny , nz
   INTEGER , DIMENSION(16) :: ipar
   INTEGER , DIMENSION(NMAX*2) :: iperm
   INTEGER , DIMENSION(NMAX*3) :: iw
   INTEGER , DIMENSION(NZMAX) :: ja , jau , ju , levs
   REAL(REAL64) , DIMENSION(NMAX,20) :: vv
   EXTERNAL amux , gen57pt , gmres , ilu0 , ilud , iluk , ilut , ilutp , milu0 , runrc
!
! End of declarations rewritten by SPAG
!
!
!     real t(2), t1, etime
!
!
!-----------------------------------------------------------------------
!     pde to be discretized is :
!---------------------------
!
!     -Lap u + gammax exp (xy)delx u + gammay exp (-xy) dely u +alpha u
!
!     where Lap = 2-D laplacean, delx = part. der. wrt x,
!     dely = part. der. wrt y.
!     gammax, gammay, and alpha are passed via the commun func.
!
!-----------------------------------------------------------------------
!
!     data for PDE:
!
   nx = 30
   ny = 30
   nz = 1
   Alpha = -50.0
   Gammax = 10.0
   Gammay = 10.0
!
!     data for preconditioner
!
   nwk = NZMAX
!
!     data for GMRES
!
   im = 10
   eps = 1.0D-07
   maxits = 100
   iout = 6
   permtol = 1.0
   ipar(2) = 2
   ipar(3) = 2
   ipar(4) = 20*NMAX
   ipar(5) = im
   ipar(6) = maxits
   fpar(1) = eps
   fpar(2) = 2.22D-16
!
!     same initial guess for gmres
!
!--------------------------------------------------------------
!     call gen57 to generate matrix in compressed sparse row format
!--------------------------------------------------------------
!
!     define part of the boundary condition here
!
   al(1) = 0.0
   al(2) = 1.0
   al(3) = 0.0
   al(4) = 0.0
   al(5) = 0.0
   al(6) = 0.0
   CALL gen57pt(nx,ny,nz,al,0,n,a,ja,ia,ju,rhs)
!
!     zero initial guess to the iterative solvers
!
   DO j = 1 , n
      xran(j) = 0.D0
   ENDDO
   PRINT * , 'RILUT:  generated a finite difference matrix'
   PRINT * , '        grid size = ' , nx , ' X ' , ny , ' X ' , nz
   PRINT * , '        matrix size = ' , n
!--------------------------------------------------------------
!     gnerate right han side = A * (1,1,1,...,1)**T
!--------------------------------------------------------------
   DO k = 1 , n
      x(k) = 1.0
   ENDDO
   CALL amux(n,x,y,a,ja,ia)
!--------------------------------------------------------------
!     test all different methods available:
!     ILU0, MILU0, ILUT and with different values of tol and lfil
!     ( from cheaper to more expensive preconditioners)
!     The more accurate the preconditioner the fewer iterations
!     are required in pgmres, in general.
!
!--------------------------------------------------------------
   DO meth = 1 , 15
      IF ( meth==2 ) THEN
         WRITE (iout,*) ' +++++ MILU(0) Preconditioner ++++ '
!        t1 = etime(t)
!        t1 = etime(t) - t1
         CALL milu0(n,a,ja,ia,au,jau,ju,iw,ierr)
      ELSEIF ( meth==3 ) THEN
         WRITE (iout,*) ' +++++ ILUT Preconditioner ++++ '
         WRITE (iout,*) ' +++++ tol = 0.0001, lfil=5 ++++ '
         tol = 0.0001
         lfil = 5
!
!        t1 = etime(t)
!         t1 = etime(t) - t1
         CALL ilut(n,a,ja,ia,lfil,tol,au,jau,ju,nwk,vv,iw,ierr)
      ELSEIF ( meth==4 ) THEN
         WRITE (iout,*) ' +++++ ILUT Preconditioner ++++ '
         WRITE (iout,*) ' +++++ tol = 0.0001, lfil=10 ++++ '
         tol = 0.0001
         lfil = 10
!
!        t1 = etime(t)
 
!        t1 = etime(t) - t1
         CALL ilut(n,a,ja,ia,lfil,tol,au,jau,ju,nwk,vv,iw,ierr)
      ELSEIF ( meth==5 ) THEN
         WRITE (iout,*) ' +++++ ILUT Preconditioner ++++ '
         WRITE (iout,*) ' +++++ tol = .0001, lfil=15 ++++ '
         tol = 0.0001
         lfil = 15
!
!        t1 = etime(t)
!        t1 = etime(t) - t1
         CALL ilut(n,a,ja,ia,lfil,tol,au,jau,ju,nwk,vv,iw,ierr)
      ELSEIF ( meth==6 ) THEN
         WRITE (iout,*) ' +++++ ILUTP Preconditioner ++++ '
         WRITE (iout,*) ' +++++ tol = 0.0001, lfil=5 ++++ '
         tol = 0.0001
         lfil = 5
!
!        t1 = etime(t)
!        t1 = etime(t) - t1
         CALL ilutp(n,a,ja,ia,lfil,tol,permtol,n,au,jau,ju,nwk,vv,iw,iperm,ierr)
      ELSEIF ( meth==7 ) THEN
         WRITE (iout,*) ' +++++ ILUTP Preconditioner ++++ '
         WRITE (iout,*) ' +++++ tol = 0.0001, lfil=10 ++++ '
         tol = 0.0001
         lfil = 10
!
!        t1 = etime(t)
!        t1 = etime(t) - t1
         CALL ilutp(n,a,ja,ia,lfil,tol,permtol,n,au,jau,ju,nwk,vv,iw,iperm,ierr)
      ELSEIF ( meth==8 ) THEN
!-----------------------------------------------------------------------------
         WRITE (iout,*) ' +++++ ILUTP Preconditioner ++++ '
         WRITE (iout,*) ' +++++ tol = .0001, lfil=15 ++++ '
         tol = 0.0001
         lfil = 15
!
!        t1 = etime(t)
!        t1 = etime(t) - t1
         CALL ilutp(n,a,ja,ia,lfil,tol,permtol,n,au,jau,ju,nwk,vv,iw,iperm,ierr)
      ELSEIF ( meth==9 ) THEN
         WRITE (iout,*) ' +++++ ILUK Preconditioner ++++ '
         WRITE (iout,*) ' +++++       lfil=0        ++++ '
         lfil = 0
!        t1 = etime(t)
         CALL iluk(n,a,ja,ia,lfil,au,jau,ju,levs,nwk,vv,iw,ierr)
         PRINT * , ' nnz for a =' , ia(n+1) - ia(1)
!        t1 = etime(t) - t1
         PRINT * , ' nnz for ilu =' , jau(n+1) - jau(1) + n
      ELSEIF ( meth==10 ) THEN
!
         WRITE (iout,*) ' +++++ ILUK Preconditioner ++++ '
         WRITE (iout,*) ' +++++       lfil=1        ++++ '
         lfil = 1
!        t1 = etime(t)
         CALL iluk(n,a,ja,ia,lfil,au,jau,ju,levs,nwk,vv,iw,ierr)
         PRINT * , ' nnz for a =' , ia(n+1) - ia(1)
!        t1 = etime(t) - t1
         PRINT * , ' nnz for ilu =' , jau(n+1) - jau(1) + n
      ELSEIF ( meth==11 ) THEN
!
         WRITE (iout,*) ' +++++ ILUK Preconditioner ++++ '
         WRITE (iout,*) ' +++++       lfil=3        ++++ '
         lfil = 3
!        t1 = etime(t)
         CALL iluk(n,a,ja,ia,lfil,au,jau,ju,levs,nwk,vv,iw,ierr)
         PRINT * , ' nnz for a =' , ia(n+1) - ia(1)
!        t1 = etime(t) - t1
         PRINT * , ' nnz for ilu =' , jau(n+1) - jau(1) + n
      ELSEIF ( meth==12 ) THEN
!
         WRITE (iout,*) ' +++++ ILUK Preconditioner ++++ '
         WRITE (iout,*) ' +++++       lfil=6        ++++ '
         lfil = 6
!        t1 = etime(t)
         CALL iluk(n,a,ja,ia,lfil,au,jau,ju,levs,nwk,vv,iw,ierr)
         PRINT * , ' nnz for a =' , ia(n+1) - ia(1)
!        t1 = etime(t) - t1
         PRINT * , ' nnz for ilu =' , jau(n+1) - jau(1) + n
      ELSEIF ( meth==13 ) THEN
!
!-----------------------------------------------------------------------
         WRITE (iout,*) ' +++++ ILUD Preconditioner ++++ '
         WRITE (iout,*) ' +++++ tol=0.075, alpha=0.0 ++++ '
         tol = 0.075
         alph = 0.0
!        t1 = etime(t)
         CALL ilud(n,a,ja,ia,alph,tol,au,jau,ju,nwk,vv,iw,ierr)
!
         PRINT * , ' nnz for a =' , ia(n+1) - ia(1)
!        t1 = etime(t) - t1
         PRINT * , ' nnz for ilu =' , jau(n+1) - jau(1) + n
      ELSEIF ( meth==14 ) THEN
!
         WRITE (iout,*) ' +++++ ILUD Preconditioner ++++ '
         WRITE (iout,*) ' +++++ tol=0.075, alpha=1.0 ++++ '
         tol = 0.075
         alph = 1.0
!        t1 = etime(t)
         CALL ilud(n,a,ja,ia,alph,tol,au,jau,ju,nwk,vv,iw,ierr)
         PRINT * , ' nnz for a =' , ia(n+1) - ia(1)
!        t1 = etime(t) - t1
         PRINT * , ' nnz for ilu =' , jau(n+1) - jau(1) + n
      ELSEIF ( meth==15 ) THEN
 
         WRITE (iout,*) ' +++++ ILUD Preconditioner ++++ '
         WRITE (iout,*) ' +++++ tol=0.01, alpha=1.0 ++++ '
         tol = 0.01
!        t1 = etime(t)
         CALL ilud(n,a,ja,ia,alph,tol,au,jau,ju,nwk,vv,iw,ierr)
         PRINT * , ' nnz for a =' , ia(n+1) - ia(1)
         PRINT * , ' nnz for ilu =' , jau(n+1) - jau(1) + n
      ELSE
         WRITE (iout,*) ' +++++ ILU(0) Preconditioner ++++ '
!        t1 = etime(t)
!        t1 = etime(t) - t1
         CALL ilu0(n,a,ja,ia,au,jau,ju,iw,ierr)
      ENDIF
!        t1 = etime(t) - t1
!        goto 100
!
!
!     check that return was succesful
!
!        print *, ' ILU factorization time ', t1
      PRINT * , ' Precon set-up returned with ierr ' , ierr
      IF ( ierr==0 ) THEN
!--------------------------------------------------------------
!     call GMRES
!--------------------------------------------------------------
         CALL runrc(n,y,x,ipar,fpar,vv,xran,a,ja,ia,au,jau,ju,gmres)
         PRINT * , 'GMRES return status = ' , ipar(1)
         WRITE (iout,*) '     '
      ENDIF
   ENDDO
!
!----------------------------------------------------------
!
   WRITE (iout,*) ' **** SOLUTION ****  '
   WRITE (iout,99001) (x(k),k=1,n)
99001 FORMAT (5D15.5)
END PROGRAM rilut
 
