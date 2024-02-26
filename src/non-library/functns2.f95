!*==afun.f90 processed by SPAG 8.04RA 14:49 22 Feb 2024
!!SPAG Open source Personal, Educational or Academic User Clemson University  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
!-----------------------------------------------------------------------
!     contains the functions needed for defining the PDE poroblems.
!
!     first for the scalar 5-point and 7-point PDE
!-----------------------------------------------------------------------
FUNCTION afun(X,Y,Z)
   USE ISO_FORTRAN_ENV                 
   IMPLICIT NONE
!
! Function and Dummy argument declarations rewritten by SPAG
!
   REAL(REAL64) :: afun
   REAL(REAL64) :: X
   REAL(REAL64) :: Y
   REAL(REAL64) :: Z
!
! End of declarations rewritten by SPAG
!
   afun = -1.0
END FUNCTION afun
!*==bfun.f90 processed by SPAG 8.04RA 14:49 22 Feb 2024
!!SPAG Open source Personal, Educational or Academic User Clemson University  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
 
FUNCTION bfun(X,Y,Z)
   USE ISO_FORTRAN_ENV                 
   IMPLICIT NONE
!
! Function and Dummy argument declarations rewritten by SPAG
!
   REAL(REAL64) :: bfun
   REAL(REAL64) :: X
   REAL(REAL64) :: Y
   REAL(REAL64) :: Z
!
! End of declarations rewritten by SPAG
!
   bfun = -1.0
END FUNCTION bfun
!*==cfun.f90 processed by SPAG 8.04RA 14:49 22 Feb 2024
!!SPAG Open source Personal, Educational or Academic User Clemson University  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
 
FUNCTION cfun(X,Y,Z)
   USE ISO_FORTRAN_ENV                 
   IMPLICIT NONE
!
! Function and Dummy argument declarations rewritten by SPAG
!
   REAL(REAL64) :: cfun
   REAL(REAL64) :: X
   REAL(REAL64) :: Y
   REAL(REAL64) :: Z
!
! End of declarations rewritten by SPAG
!
   cfun = -1.0D0
END FUNCTION cfun
!*==dfun.f90 processed by SPAG 8.04RA 14:49 22 Feb 2024
!!SPAG Open source Personal, Educational or Academic User Clemson University  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
 
FUNCTION dfun(X,Y,Z)
   USE ISO_FORTRAN_ENV                 
   IMPLICIT NONE
!
! Function and Dummy argument declarations rewritten by SPAG
!
   REAL(REAL64) :: dfun
   REAL(REAL64) :: X
   REAL(REAL64) :: Y
   REAL(REAL64) :: Z
!
! End of declarations rewritten by SPAG
!
   dfun = 10.D0
END FUNCTION dfun
!*==efun.f90 processed by SPAG 8.04RA 14:49 22 Feb 2024
!!SPAG Open source Personal, Educational or Academic User Clemson University  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
 
FUNCTION efun(X,Y,Z)
   USE ISO_FORTRAN_ENV                 
   IMPLICIT NONE
!
! Function and Dummy argument declarations rewritten by SPAG
!
   REAL(REAL64) :: efun
   REAL(REAL64) :: X
   REAL(REAL64) :: Y
   REAL(REAL64) :: Z
!
! End of declarations rewritten by SPAG
!
   efun = 0.0D0
END FUNCTION efun
!*==ffun.f90 processed by SPAG 8.04RA 14:49 22 Feb 2024
!!SPAG Open source Personal, Educational or Academic User Clemson University  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
 
FUNCTION ffun(X,Y,Z)
   USE ISO_FORTRAN_ENV                 
   IMPLICIT NONE
!
! Function and Dummy argument declarations rewritten by SPAG
!
   REAL(REAL64) :: ffun
   REAL(REAL64) :: X
   REAL(REAL64) :: Y
   REAL(REAL64) :: Z
!
! End of declarations rewritten by SPAG
!
   ffun = 0.0
END FUNCTION ffun
!*==gfun.f90 processed by SPAG 8.04RA 14:49 22 Feb 2024
!!SPAG Open source Personal, Educational or Academic User Clemson University  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
 
FUNCTION gfun(X,Y,Z)
   USE ISO_FORTRAN_ENV                 
   IMPLICIT NONE
!
! Function and Dummy argument declarations rewritten by SPAG
!
   REAL(REAL64) :: gfun
   REAL(REAL64) :: X
   REAL(REAL64) :: Y
   REAL(REAL64) :: Z
!
! End of declarations rewritten by SPAG
!
   gfun = 0.0
END FUNCTION gfun
!*==hfun.f90 processed by SPAG 8.04RA 14:49 22 Feb 2024
!!SPAG Open source Personal, Educational or Academic User Clemson University  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
 
FUNCTION hfun(X,Y,Z)
   USE ISO_FORTRAN_ENV                 
   IMPLICIT NONE
!
! Function and Dummy argument declarations rewritten by SPAG
!
   REAL(REAL64) :: hfun
   REAL(REAL64) :: X
   REAL(REAL64) :: Y
   REAL(REAL64) :: Z
!
! End of declarations rewritten by SPAG
!
   hfun = 0.0
END FUNCTION hfun
!*==betfun.f90 processed by SPAG 8.04RA 14:49 22 Feb 2024
!!SPAG Open source Personal, Educational or Academic User Clemson University  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
 
FUNCTION betfun(Side,X,Y,Z)
   USE ISO_FORTRAN_ENV                 
   IMPLICIT NONE
!
! Function and Dummy argument declarations rewritten by SPAG
!
   REAL(REAL64) :: betfun
   CHARACTER(2) :: Side
   REAL(REAL64) :: X
   REAL(REAL64) :: Y
   REAL(REAL64) :: Z
!
! End of declarations rewritten by SPAG
!
   betfun = 1.0
END FUNCTION betfun
!*==gamfun.f90 processed by SPAG 8.04RA 14:49 22 Feb 2024
!!SPAG Open source Personal, Educational or Academic User Clemson University  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
 
FUNCTION gamfun(Side,X,Y,Z)
   USE ISO_FORTRAN_ENV                 
   IMPLICIT NONE
!
! Function and Dummy argument declarations rewritten by SPAG
!
   REAL(REAL64) :: gamfun
   CHARACTER(2) , INTENT(IN) :: Side
   REAL(REAL64) :: X
   REAL(REAL64) :: Y
   REAL(REAL64) :: Z
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
!*==afunbl.f90 processed by SPAG 8.04RA 14:49 22 Feb 2024
!!SPAG Open source Personal, Educational or Academic User Clemson University  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
 
!-----------------------------------------------------------------------
!     functions for the block PDE's
!-----------------------------------------------------------------------
SUBROUTINE afunbl(Nfree,X,Y,Z,Coeff)
   USE ISO_FORTRAN_ENV                 
   IMPLICIT NONE
!
! Dummy argument declarations rewritten by SPAG
!
   INTEGER , INTENT(IN) :: Nfree
   REAL(REAL64) :: X
   REAL(REAL64) :: Y
   REAL(REAL64) :: Z
   REAL(REAL64) , INTENT(OUT) , DIMENSION(100) :: Coeff
!
! Local variable declarations rewritten by SPAG
!
   INTEGER :: i , j
!
! End of declarations rewritten by SPAG
!
   DO j = 1 , Nfree
      DO i = 1 , Nfree
         Coeff((j-1)*Nfree+i) = 0.0D0
      ENDDO
      Coeff((j-1)*Nfree+j) = -1.0D0
   ENDDO
END SUBROUTINE afunbl
!*==bfunbl.f90 processed by SPAG 8.04RA 14:49 22 Feb 2024
!!SPAG Open source Personal, Educational or Academic User Clemson University  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
 
SUBROUTINE bfunbl(Nfree,X,Y,Z,Coeff)
   USE ISO_FORTRAN_ENV                 
   IMPLICIT NONE
!
! Dummy argument declarations rewritten by SPAG
!
   INTEGER , INTENT(IN) :: Nfree
   REAL(REAL64) :: X
   REAL(REAL64) :: Y
   REAL(REAL64) :: Z
   REAL(REAL64) , INTENT(OUT) , DIMENSION(100) :: Coeff
!
! Local variable declarations rewritten by SPAG
!
   INTEGER :: i , j
!
! End of declarations rewritten by SPAG
!
   DO j = 1 , Nfree
      DO i = 1 , Nfree
         Coeff((j-1)*Nfree+i) = 0.0D0
      ENDDO
      Coeff((j-1)*Nfree+j) = -1.0D0
   ENDDO
END SUBROUTINE bfunbl
!*==cfunbl.f90 processed by SPAG 8.04RA 14:49 22 Feb 2024
!!SPAG Open source Personal, Educational or Academic User Clemson University  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
 
SUBROUTINE cfunbl(Nfree,X,Y,Z,Coeff)
   USE ISO_FORTRAN_ENV                 
   IMPLICIT NONE
!
! Dummy argument declarations rewritten by SPAG
!
   INTEGER , INTENT(IN) :: Nfree
   REAL(REAL64) :: X
   REAL(REAL64) :: Y
   REAL(REAL64) :: Z
   REAL(REAL64) , INTENT(OUT) , DIMENSION(100) :: Coeff
!
! Local variable declarations rewritten by SPAG
!
   INTEGER :: i , j
!
! End of declarations rewritten by SPAG
!
   DO j = 1 , Nfree
      DO i = 1 , Nfree
         Coeff((j-1)*Nfree+i) = 0.0D0
      ENDDO
      Coeff((j-1)*Nfree+j) = -1.0D0
   ENDDO
END SUBROUTINE cfunbl
!*==dfunbl.f90 processed by SPAG 8.04RA 14:49 22 Feb 2024
!!SPAG Open source Personal, Educational or Academic User Clemson University  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
 
SUBROUTINE dfunbl(Nfree,X,Y,Z,Coeff)
   USE ISO_FORTRAN_ENV                 
   IMPLICIT NONE
!
! Dummy argument declarations rewritten by SPAG
!
   INTEGER , INTENT(IN) :: Nfree
   REAL(REAL64) :: X
   REAL(REAL64) :: Y
   REAL(REAL64) :: Z
   REAL(REAL64) , INTENT(OUT) , DIMENSION(100) :: Coeff
!
! Local variable declarations rewritten by SPAG
!
   INTEGER :: i , j
!
! End of declarations rewritten by SPAG
!
   DO j = 1 , Nfree
      DO i = 1 , Nfree
         Coeff((j-1)*Nfree+i) = 0.0D0
      ENDDO
   ENDDO
END SUBROUTINE dfunbl
!*==efunbl.f90 processed by SPAG 8.04RA 14:49 22 Feb 2024
!!SPAG Open source Personal, Educational or Academic User Clemson University  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
 
SUBROUTINE efunbl(Nfree,X,Y,Z,Coeff)
   USE ISO_FORTRAN_ENV                 
   IMPLICIT NONE
!
! Dummy argument declarations rewritten by SPAG
!
   INTEGER , INTENT(IN) :: Nfree
   REAL(REAL64) :: X
   REAL(REAL64) :: Y
   REAL(REAL64) :: Z
   REAL(REAL64) , INTENT(OUT) , DIMENSION(100) :: Coeff
!
! Local variable declarations rewritten by SPAG
!
   INTEGER :: i , j
!
! End of declarations rewritten by SPAG
!
   DO j = 1 , Nfree
      DO i = 1 , Nfree
         Coeff((j-1)*Nfree+i) = 0.0D0
      ENDDO
   ENDDO
END SUBROUTINE efunbl
!*==ffunbl.f90 processed by SPAG 8.04RA 14:49 22 Feb 2024
!!SPAG Open source Personal, Educational or Academic User Clemson University  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
 
SUBROUTINE ffunbl(Nfree,X,Y,Z,Coeff)
   USE ISO_FORTRAN_ENV                 
   IMPLICIT NONE
!
! Dummy argument declarations rewritten by SPAG
!
   INTEGER , INTENT(IN) :: Nfree
   REAL(REAL64) :: X
   REAL(REAL64) :: Y
   REAL(REAL64) :: Z
   REAL(REAL64) , INTENT(OUT) , DIMENSION(100) :: Coeff
!
! Local variable declarations rewritten by SPAG
!
   INTEGER :: i , j
!
! End of declarations rewritten by SPAG
!
   DO j = 1 , Nfree
      DO i = 1 , Nfree
         Coeff((j-1)*Nfree+i) = 0.0D0
      ENDDO
   ENDDO
END SUBROUTINE ffunbl
!*==gfunbl.f90 processed by SPAG 8.04RA 14:49 22 Feb 2024
!!SPAG Open source Personal, Educational or Academic User Clemson University  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
 
SUBROUTINE gfunbl(Nfree,X,Y,Z,Coeff)
   USE ISO_FORTRAN_ENV                 
   IMPLICIT NONE
!
! Dummy argument declarations rewritten by SPAG
!
   INTEGER , INTENT(IN) :: Nfree
   REAL(REAL64) :: X
   REAL(REAL64) :: Y
   REAL(REAL64) :: Z
   REAL(REAL64) , INTENT(OUT) , DIMENSION(100) :: Coeff
!
! Local variable declarations rewritten by SPAG
!
   INTEGER :: i , j
!
! End of declarations rewritten by SPAG
!
   DO j = 1 , Nfree
      DO i = 1 , Nfree
         Coeff((j-1)*Nfree+i) = 0.0D0
      ENDDO
   ENDDO
END SUBROUTINE gfunbl
!*==xyk.f90 processed by SPAG 8.04RA 14:49 22 Feb 2024
!!SPAG Open source Personal, Educational or Academic User Clemson University  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
!-----------------------------------------------------------------------
!     The material property function xyk for the
!     finite element problem
!-----------------------------------------------------------------------
!        subroutine xyk(nel,xyke,x,y,ijk,node)
!        implicit real*8 (a-h,o-z)
!        dimension xyke(2,2), x(*), y(*), ijk(node,*)
!c
!c this is the identity matrix.
!c
!        xyke(1,1) = 1.0d0
!        xyke(2,2) = 1.0d0
!        xyke(1,2) = 0.0d0
!        xyke(2,1) = 0.0d0
!c
!        return
!        end
!
SUBROUTINE xyk(Xyke,X,Y)
   USE ISO_FORTRAN_ENV                 
   IMPLICIT NONE
!
! Dummy argument declarations rewritten by SPAG
!
   REAL(REAL64) , INTENT(OUT) , DIMENSION(2,2) :: Xyke
   REAL(REAL64) , INTENT(IN) :: X
   REAL(REAL64) , INTENT(IN) :: Y
!
! End of declarations rewritten by SPAG
!
 
   Xyke(1,1) = 1.
   Xyke(1,1) = exp(X+Y)
   Xyke(1,2) = 0.
   Xyke(2,1) = 0.
   Xyke(2,2) = 1.
   Xyke(2,2) = exp(X+Y)
END SUBROUTINE xyk
!*==funb.f90 processed by SPAG 8.04RA 14:49 22 Feb 2024
!!SPAG Open source Personal, Educational or Academic User Clemson University  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
 
FUNCTION funb(X,Y)
   USE ISO_FORTRAN_ENV                 
   IMPLICIT NONE
!
! Function and Dummy argument declarations rewritten by SPAG
!
   REAL(REAL64) :: funb
   REAL(REAL64) , INTENT(IN) :: X
   REAL(REAL64) :: Y
!
! End of declarations rewritten by SPAG
!
 
   funb = 0.
   funb = 2.5
   funb = 2*X
END FUNCTION funb
!*==func.f90 processed by SPAG 8.04RA 14:49 22 Feb 2024
!!SPAG Open source Personal, Educational or Academic User Clemson University  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
 
FUNCTION func(X,Y)
   USE ISO_FORTRAN_ENV                 
   IMPLICIT NONE
!
! Function and Dummy argument declarations rewritten by SPAG
!
   REAL(REAL64) :: func
   REAL(REAL64) :: X
   REAL(REAL64) , INTENT(IN) :: Y
!
! End of declarations rewritten by SPAG
!
 
   func = 0.
   func = -1.5
   func = -5*Y
END FUNCTION func
!*==fung.f90 processed by SPAG 8.04RA 14:49 22 Feb 2024
!!SPAG Open source Personal, Educational or Academic User Clemson University  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
 
FUNCTION fung(X,Y)
!  Right hand side corresponding to the exact solution of
!   u = exp(x+y)*x*(1.-x)*y*(1.-y)
!   (That exact solution is defined in the function exact)
   USE ISO_FORTRAN_ENV                 
   IMPLICIT NONE
!
! Function and Dummy argument declarations rewritten by SPAG
!
   REAL(REAL64) :: fung
   REAL(REAL64) , INTENT(IN) :: X
   REAL(REAL64) , INTENT(IN) :: Y
!
! Local variable declarations rewritten by SPAG
!
   REAL(REAL64) :: r
!
! End of declarations rewritten by SPAG
!
 
!        fung = 1. + x
!        fung = 2*y*(1.-y) + 2*x*(1.-x)
!        fung = 2*y*(1.-y) + 2*x*(1.-x) +2.5*(1-2*x)*y*(1-y)
!     1        - 1.5*(1.-2*y)*x*(1-x)
   r = exp(X+Y)
   fung = r*r*((X*X+3.*X)*Y*(1.-Y)+(Y*Y+3.*Y)*X*(1.-X)) + r*(r-2.*X)*(X*X+X-1.)*Y*(1.-Y) + r*(r+5.*Y)*(Y*Y+Y-1.)*X*(1.-X)
END FUNCTION fung
!*==exact.f90 processed by SPAG 8.04RA 14:49 22 Feb 2024
!!SPAG Open source Personal, Educational or Academic User Clemson University  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
 
FUNCTION exact(X,Y)
!  Exact Solution.
   USE ISO_FORTRAN_ENV                 
   IMPLICIT NONE
!
! Function and Dummy argument declarations rewritten by SPAG
!
   REAL(REAL64) :: exact
   REAL(REAL64) , INTENT(IN) :: X
   REAL(REAL64) , INTENT(IN) :: Y
!
! End of declarations rewritten by SPAG
!
 
   exact = exp(X+Y)*X*(1.-X)*Y*(1.-Y)
 
END FUNCTION exact
