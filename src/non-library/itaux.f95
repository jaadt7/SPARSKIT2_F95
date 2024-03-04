!*==runrc.f90 processed by SPAG 8.04RA 14:49 22 Feb 2024
!!SPAG Open source Personal, Educational or Academic User Clemson University  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
SUBROUTINE runrc(N,Rhs,Sol,Ipar,Fpar,Wk,Guess,A,Ja,Ia,Au,Jau,Ju,solver)
   USE ISO_FORTRAN_ENV                 
   IMPLICIT NONE
!
! Dummy argument declarations rewritten by SPAG
!
   INTEGER :: N
   REAL(REAL64) , DIMENSION(N) :: Rhs
   REAL(REAL64) , INTENT(INOUT) , DIMENSION(N) :: Sol
   INTEGER , INTENT(INOUT) , DIMENSION(16) :: Ipar
   REAL(REAL64) , DIMENSION(16) :: Fpar
   REAL(REAL64) , INTENT(INOUT) , DIMENSION(*) :: Wk
   REAL(REAL64) , INTENT(IN) , DIMENSION(N) :: Guess
   REAL(REAL64) , DIMENSION(*) :: A
   INTEGER , DIMENSION(*) :: Ja
   INTEGER , DIMENSION(N+1) :: Ia
   REAL(REAL64) , DIMENSION(*) :: Au
   INTEGER , DIMENSION(*) :: Jau
   INTEGER , DIMENSION(*) :: Ju
   EXTERNAL solver
!
! Local variable declarations rewritten by SPAG
!
   REAL(REAL64) , EXTERNAL :: dnrm2
   INTEGER :: i , iou
   INTEGER , SAVE :: its
   REAL(REAL64) , SAVE :: res
   EXTERNAL amux , atmux , lusol , lutsol
!
! End of declarations rewritten by SPAG
!
!-----------------------------------------------------------------------
!     the actual tester. It starts the iterative linear system solvers
!     with a initial guess suppied by the user.
!
!     The structure {au, jau, ju} is assumed to have the output from
!     the ILU* routines in ilut.f.
!
!-----------------------------------------------------------------------
!     local variables
!
!     real dtime, dt(2), time
!     external dtime
!
!     ipar(2) can be 0, 1, 2, please don't use 3
!
   IF ( Ipar(2)>2 ) THEN
      PRINT * , 'I can not do both left and right preconditioning.'
      RETURN
   ENDIF
!
!     normal execution
!
   its = 0
   res = 0.0D0
!
   DO i = 1 , N
      Sol(i) = Guess(i)
   ENDDO
!
   iou = 6
   Ipar(1) = 0
   SPAG_Loop_1_1: DO
!     time = dtime(dt)
      CALL solver(N,Rhs,Sol,Ipar,Fpar,Wk)
!
!     output the residuals
!
      IF ( Ipar(7)/=its ) THEN
         WRITE (iou,*) its , real(res)
         its = Ipar(7)
      ENDIF
      res = Fpar(5)
!
      IF ( Ipar(1)==1 ) THEN
         CALL amux(N,Wk(Ipar(8)),Wk(Ipar(9)),A,Ja,Ia)
         CYCLE
      ELSEIF ( Ipar(1)==2 ) THEN
         CALL atmux(N,Wk(Ipar(8)),Wk(Ipar(9)),A,Ja,Ia)
         CYCLE
      ELSEIF ( Ipar(1)==3 .OR. Ipar(1)==5 ) THEN
         CALL lusol(N,Wk(Ipar(8)),Wk(Ipar(9)),Au,Jau,Ju)
         CYCLE
      ELSEIF ( Ipar(1)==4 .OR. Ipar(1)==6 ) THEN
         CALL lutsol(N,Wk(Ipar(8)),Wk(Ipar(9)),Au,Jau,Ju)
         CYCLE
      ELSEIF ( Ipar(1)<=0 ) THEN
         IF ( Ipar(1)==0 ) THEN
            PRINT * , 'Iterative solver has satisfied convergence test.'
         ELSEIF ( Ipar(1)==-1 ) THEN
            PRINT * , 'Iterative solver has iterated too many times.'
         ELSEIF ( Ipar(1)==-2 ) THEN
            PRINT * , 'Iterative solver was not given enough work space.'
            PRINT * , 'The work space should at least have ' , Ipar(4) , ' elements.'
         ELSEIF ( Ipar(1)==-3 ) THEN
            PRINT * , 'Iterative solver is facing a break-down.'
         ELSE
            PRINT * , 'Iterative solver terminated. code =' , Ipar(1)
         ENDIF
      ENDIF
!     time = dtime(dt)
      WRITE (iou,*) Ipar(7) , real(Fpar(6))
      WRITE (iou,*) '# retrun code =' , Ipar(1) , '	convergence rate =' , Fpar(7)
!     write (iou, *) '# total execution time (sec)', time
!
!     check the error
!
      CALL amux(N,Sol,Wk,A,Ja,Ia)
      DO i = 1 , N
         Wk(N+i) = Sol(i) - 1.0D0
         Wk(i) = Wk(i) - Rhs(i)
      ENDDO
      WRITE (iou,*) '# the actual residual norm is' , dnrm2(N,Wk,1)
      WRITE (iou,*) '# the error norm is' , dnrm2(N,Wk(1+N),1)
!
      IF ( iou/=6 ) CLOSE (iou)
      EXIT SPAG_Loop_1_1
   ENDDO SPAG_Loop_1_1
END SUBROUTINE runrc
!*==distdot.f90 processed by SPAG 8.04RA 14:49 22 Feb 2024
!!SPAG Open source Personal, Educational or Academic User Clemson University  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
!-----end-of-runrc
!-----------------------------------------------------------------------
FUNCTION distdot(N,X,Ix,Y,Iy)
   USE ISO_FORTRAN_ENV                 
   IMPLICIT NONE
!
! Function and Dummy argument declarations rewritten by SPAG
!
   REAL(REAL64) :: distdot
   INTEGER :: N
   REAL(REAL64) , DIMENSION(*) :: X
   INTEGER :: Ix
   REAL(REAL64) , DIMENSION(*) :: Y
   INTEGER :: Iy
!
! Local variable declarations rewritten by SPAG
!
   REAL(REAL64) , EXTERNAL :: ddot
!
! End of declarations rewritten by SPAG
!
   distdot = ddot(N,X,Ix,Y,Iy)
END FUNCTION distdot
!*==afun.f90 processed by SPAG 8.04RA 14:49 22 Feb 2024
!!SPAG Open source Personal, Educational or Academic User Clemson University  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
!-----end-of-distdot
!-----------------------------------------------------------------------
!
FUNCTION afun()
   USE ISO_FORTRAN_ENV                 
   IMPLICIT NONE
!
! Function and Dummy argument declarations rewritten by SPAG
!
   REAL(REAL64) :: afun
!
! End of declarations rewritten by SPAG
!
   afun = -1.0D0
END FUNCTION afun
!*==bfun.f90 processed by SPAG 8.04RA 14:49 22 Feb 2024
!!SPAG Open source Personal, Educational or Academic User Clemson University  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
 
FUNCTION bfun()
   USE ISO_FORTRAN_ENV                 
   IMPLICIT NONE
!
! Function and Dummy argument declarations rewritten by SPAG
!
   REAL(REAL64) :: bfun
!
! End of declarations rewritten by SPAG
!
   bfun = -1.0D0
END FUNCTION bfun
!*==cfun.f90 processed by SPAG 8.04RA 14:49 22 Feb 2024
!!SPAG Open source Personal, Educational or Academic User Clemson University  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
 
FUNCTION cfun()
   USE ISO_FORTRAN_ENV                 
   IMPLICIT NONE
!
! Function and Dummy argument declarations rewritten by SPAG
!
   REAL(REAL64) :: cfun
!
! End of declarations rewritten by SPAG
!
   cfun = -1.0D0
END FUNCTION cfun
!*==dfun.f90 processed by SPAG 8.04RA 14:49 22 Feb 2024
!!SPAG Open source Personal, Educational or Academic User Clemson University  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
 
FUNCTION dfun(X,Y)
   USE ISO_FORTRAN_ENV                 
   IMPLICIT NONE
!
! COMMON variable declarations rewritten by SPAG
!
   REAL(REAL64) :: Alpha , Gammax , Gammay
   COMMON /func  / Gammax , Gammay , Alpha
!
! Function and Dummy argument declarations rewritten by SPAG
!
   REAL(REAL64) :: dfun
   REAL(REAL64) , INTENT(IN) :: X
   REAL(REAL64) , INTENT(IN) :: Y
 
!
! End of declarations rewritten by SPAG
!
   dfun = Gammax*exp(X*Y)
END FUNCTION dfun
!*==efun.f90 processed by SPAG 8.04RA 14:49 22 Feb 2024
!!SPAG Open source Personal, Educational or Academic User Clemson University  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
 
FUNCTION efun(X,Y)
   USE ISO_FORTRAN_ENV                 
   IMPLICIT NONE
!
! COMMON variable declarations rewritten by SPAG
!
   REAL(REAL64) :: Alpha , Gammax , Gammay
   COMMON /func  / Gammax , Gammay , Alpha
!
! Function and Dummy argument declarations rewritten by SPAG
!
   REAL(REAL64) :: efun
   REAL(REAL64) , INTENT(IN) :: X
   REAL(REAL64) , INTENT(IN) :: Y
!
! End of declarations rewritten by SPAG
!
   efun = Gammay*exp(-X*Y)
END FUNCTION efun
!*==ffun.f90 processed by SPAG 8.04RA 14:49 22 Feb 2024
!!SPAG Open source Personal, Educational or Academic User Clemson University  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
 
FUNCTION ffun()
   USE ISO_FORTRAN_ENV                 
   IMPLICIT NONE
!
! Function and Dummy argument declarations rewritten by SPAG
!
   REAL(REAL64) :: ffun
!
! End of declarations rewritten by SPAG
!
   ffun = 0.0D0
END FUNCTION ffun
!*==gfun.f90 processed by SPAG 8.04RA 14:49 22 Feb 2024
!!SPAG Open source Personal, Educational or Academic User Clemson University  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
 
FUNCTION gfun()
   USE ISO_FORTRAN_ENV                 
   IMPLICIT NONE
!
! COMMON variable declarations rewritten by SPAG
!
   REAL(REAL64) :: Alpha , Gammax , Gammay
   COMMON /func  / Gammax , Gammay , Alpha
!
! Function and Dummy argument declarations rewritten by SPAG
!
   REAL(REAL64) :: gfun
!
! End of declarations rewritten by SPAG
!
   gfun = Alpha
END FUNCTION gfun
!*==hfun.f90 processed by SPAG 8.04RA 14:49 22 Feb 2024
!!SPAG Open source Personal, Educational or Academic User Clemson University  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
 
FUNCTION hfun(X,Y,Z)
   USE ISO_FORTRAN_ENV                 
   IMPLICIT NONE
!
! COMMON variable declarations rewritten by SPAG
!
   REAL(REAL64) :: Alpha , Gammax , Gammay
   COMMON /func  / Gammax , Gammay , Alpha
!
! Function and Dummy argument declarations rewritten by SPAG
!
   REAL(REAL64) :: hfun
   REAL(REAL64) , INTENT(IN) :: X
   REAL(REAL64) , INTENT(IN) :: Y
   REAL(REAL64) , INTENT(IN) :: Z
!
! End of declarations rewritten by SPAG
!
   hfun = Alpha*sin(Gammax*X+Gammay*Y-Z)
END FUNCTION hfun
!*==betfun.f90 processed by SPAG 8.04RA 14:49 22 Feb 2024
!!SPAG Open source Personal, Educational or Academic User Clemson University  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
 
 
FUNCTION betfun()
   USE ISO_FORTRAN_ENV                 
   IMPLICIT NONE
!
! Function and Dummy argument declarations rewritten by SPAG
!
   REAL(REAL64) :: betfun

!
! End of declarations rewritten by SPAG
!
   betfun = 1.0
END FUNCTION betfun
!*==gamfun.f90 processed by SPAG 8.04RA 14:49 22 Feb 2024
!!SPAG Open source Personal, Educational or Academic User Clemson University  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
 
FUNCTION gamfun(Side)
   USE ISO_FORTRAN_ENV                 
   IMPLICIT NONE
!
! Function and Dummy argument declarations rewritten by SPAG
!
   REAL(REAL64) :: gamfun
   CHARACTER(2) , INTENT(IN) :: Side
!
! End of declarations rewritten by SPAG
!
   IF ( Side=='x2' ) THEN
      gamfun = 5.0
   ELSEIF ( Side=='y1' ) THEN
      gamfun = 2.0
   ELSEIF ( Side=='y2' ) THEN
      gamfun = 7.0
   ELSE
      gamfun = 0.0
   ENDIF
END FUNCTION gamfun
