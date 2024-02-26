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
   afun = -1.0D0
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
   bfun = -1.0D0
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
! Local variable declarations rewritten by SPAG
!
   REAL , SAVE :: gamma
!
! End of declarations rewritten by SPAG
!
   DATA gamma/100.0/
!     dfun = gamma * exp( x * y )
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
! Local variable declarations rewritten by SPAG
!
   REAL , SAVE :: gamma
!
! End of declarations rewritten by SPAG
!
   DATA gamma/100.0/
!     efun = gamma * exp( (- x) * y )
   efun = 0.D0
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
SUBROUTINE xyk(Nel,Xyke,X,Y,Ijk,Node)
   USE ISO_FORTRAN_ENV                 
   IMPLICIT NONE
!
! Dummy argument declarations rewritten by SPAG
!
   INTEGER , INTENT(IN) :: Node
   INTEGER :: Nel
   REAL(REAL64) , INTENT(OUT) , DIMENSION(2,2) :: Xyke
   REAL(REAL64) , DIMENSION(*) :: X
   REAL(REAL64) , DIMENSION(*) :: Y
   INTEGER , DIMENSION(Node,*) :: Ijk
!
! End of declarations rewritten by SPAG
!
!
!     this is the identity matrix.
!
   Xyke(1,1) = 1.0D0
   Xyke(2,2) = 1.0D0
   Xyke(1,2) = 0.0D0
   Xyke(2,1) = 0.0D0
 
END SUBROUTINE xyk
