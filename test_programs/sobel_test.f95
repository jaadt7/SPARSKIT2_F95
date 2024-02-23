!*==rsobel.f90 processed by SPAG 8.04RA 11:58 23 Feb 2024
!!SPAG Open source Personal, Educational or Academic User Clemson University  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
PROGRAM rsobel
   USE ISO_FORTRAN_ENV                 
   IMPLICIT NONE
!
! Local variable declarations rewritten by SPAG
!
   REAL(REAL64) , DIMENSION(1:1000) :: a , b , c
   INTEGER , DIMENSION(1:200) :: ia , ib , ic
   INTEGER :: ierr , n , ncolc , nrowc
   INTEGER , DIMENSION(1:1000) :: ja , jb , jc
   EXTERNAL dump , sobel
!
! End of declarations rewritten by SPAG
!
 
   WRITE (*,'(1x, 9hInput n:  ,$)')
   READ * , n
   CALL sobel(n,nrowc,ncolc,c,jc,ic,a,ja,ia,b,jb,ib,1000,ierr)
   PRINT * , 'ierr =' , ierr
   PRINT * , 'Nrow =' , nrowc , '	Ncol =' , ncolc
   CALL dump(1,nrowc,.TRUE.,c,jc,ic,6)
END PROGRAM rsobel
