!*==matprod.f90 processed by SPAG 8.04RA 18:15 22 Feb 2024
!!SPAG Open source Personal, Educational or Academic User Clemson University  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
PROGRAM matprod
!-----------------------------------------------------------------------
!      test program for some routines in BLASSM.f
!-----------------------------------------------------------------------
!      Last update: May 2, 1994
!-----------------------------------------------------------------------
   USE ISO_FORTRAN_ENV                 
   IMPLICIT NONE
!
! PARAMETER definitions rewritten by SPAG
!
   INTEGER , PARAMETER :: NXMAX = 30 , NMX = NXMAX*NXMAX , NZMAX = 7*NMX
!
! Local variable declarations rewritten by SPAG
!
   REAL(REAL64) , DIMENSION(NZMAX) :: a , b , c
   REAL(REAL64) , DIMENSION(6) :: al
   INTEGER , DIMENSION(NMX+1) :: ia , ib
   INTEGER , DIMENSION(NZMAX) :: ic , ja , jb , jc
   INTEGER :: ierr , ifmt , iout , j , jj , job , k , n , nx , ny , nz
   INTEGER , DIMENSION(NMX) :: iw
   REAL(REAL64) , DIMENSION(NMX) :: rhs , x , y , y1
   REAL(REAL64) :: s
   EXTERNAL amub , aplsb , aplsbt , apmbt , dump , gen57pt , ope , opet , ydfnorm
!
! End of declarations rewritten by SPAG
!
!-----------------------------------------------------------------------
 
   nx = 20
   ny = 20
   nz = 1
   al(1) = 1.0D0
   al(2) = 0.0D0
   al(3) = 2.3D1
   al(4) = 0.4D0
   al(5) = 0.0D0
   al(6) = 8.2D-2
   iout = 8
!-----------------------------------------------------------------------
   CALL gen57pt(nx,ny,nz,al,0,n,a,ja,ia,iw,rhs)
   CALL gen57pt(ny,nx,nz,al,0,n,b,jb,ib,iw,rhs)
!
   s = 3.812
!
!      call aplsb1(n,n,a,ja,ia,s,b,jb,ib,c,jc,ic,nzmax,ierr)
   CALL aplsb(n,n,a,ja,ia,s,b,jb,ib,c,jc,ic,NZMAX,iw,ierr)
   IF ( ierr/=0 ) PRINT * , ' ierr = ' , ierr
!
!       call dump (1,n,.true.,c,jc,ic,9)
!
   DO k = 1 , n
      x(k) = real(k)/real(n)
   ENDDO
!
   CALL ope(n,x,y1,a,ja,ia)
   CALL ope(n,x,y,b,jb,ib)
   DO j = 1 , n
      y1(j) = s*y(j) + y1(j)
   ENDDO
!
   CALL ope(n,x,y,c,jc,ic)
!------------------------------------------------------
   WRITE (6,*) ' ------------ checking APLSB --------------'
   CALL ydfnorm(n,y1,y,6)
!------------------------------------------------------
   ifmt = 103
!
   job = -1
!--------
   DO jj = 1 , 2
      WRITE (9,*) 'DUMP A____________________________'
      CALL dump(1,n,.TRUE.,a,ja,ia,9)
      WRITE (9,*) 'DUMP B____________________________'
      CALL dump(1,n,.TRUE.,b,jb,ib,9)
      CALL apmbt(n,n,job,a,ja,ia,b,jb,ib,c,jc,ic,NZMAX,iw,ierr)
      WRITE (9,*) 'DUMP C____________________________'
      CALL dump(1,n,.TRUE.,c,jc,ic,9)
      IF ( ierr/=0 ) PRINT * , ' ierr = ' , ierr
      CALL ope(n,x,y1,a,ja,ia)
      CALL opet(n,x,y,b,jb,ib)
      s = real(job)
      DO j = 1 , n
         y1(j) = y1(j) + s*y(j)
      ENDDO
!
      CALL ope(n,x,y,c,jc,ic)
!------------xs------------------------------------------
      WRITE (6,*) '  '
      WRITE (6,*) ' ------------ checking APMBT---------------'
      WRITE (6,*) ' ------------ with JOB = ' , job , ' -------------'
      CALL ydfnorm(n,y1,y,6)
!------------------------------------------------------
      job = job + 2
   ENDDO
!
   s = 0.1232445
   CALL aplsbt(n,n,a,ja,ia,s,b,jb,ib,c,jc,ic,NZMAX,iw,ierr)
!
   IF ( ierr/=0 ) PRINT * , ' ierr = ' , ierr
   CALL ope(n,x,y1,a,ja,ia)
   CALL opet(n,x,y,b,jb,ib)
   DO j = 1 , n
      y1(j) = y1(j) + s*y(j)
   ENDDO
!
   CALL ope(n,x,y,c,jc,ic)
!------------------------------------------------------
!------------------------------------------------------
   WRITE (6,*) '  '
   WRITE (6,*) ' ------------ checking APLSBT---------------'
   CALL ydfnorm(n,y1,y,6)
!-----------------------------------------------------------------------
! testing products
!-----------------------------------------------------------------------
   job = 1
   CALL amub(n,n,job,a,ja,ia,b,jb,ib,c,jc,ic,NZMAX,iw,ierr)
!
   IF ( ierr/=0 ) PRINT * , ' ierr = ' , ierr
   CALL ope(n,x,y,b,jb,ib)
   CALL ope(n,y,y1,a,ja,ia)
!
   CALL ope(n,x,y,c,jc,ic)
!-----------------------------------------------------------------------
   WRITE (6,*) '  '
   WRITE (6,*) ' ------------ checking AMUB  ---------------'
   CALL ydfnorm(n,y1,y,6)
!
END PROGRAM matprod
!*==ope.f90 processed by SPAG 8.04RA 18:15 22 Feb 2024
!!SPAG Open source Personal, Educational or Academic User Clemson University  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
!
!
SUBROUTINE ope(N,X,Y,A,Ja,Ia)
   USE ISO_FORTRAN_ENV                 
   IMPLICIT NONE
!
! Dummy argument declarations rewritten by SPAG
!
   INTEGER , INTENT(IN) :: N
   REAL(REAL64) , INTENT(IN) , DIMENSION(1) :: X
   REAL(REAL64) , INTENT(INOUT) , DIMENSION(1) :: Y
   REAL(REAL64) , INTENT(IN) , DIMENSION(1) :: A
   INTEGER , INTENT(IN) , DIMENSION(*) :: Ja
   INTEGER , INTENT(IN) , DIMENSION(*) :: Ia
!
! Local variable declarations rewritten by SPAG
!
   INTEGER :: i , k , k1 , k2
!
! End of declarations rewritten by SPAG
!
! sparse matrix * vector multiplication
!
   DO i = 1 , N
      k1 = Ia(i)
      k2 = Ia(i+1) - 1
      Y(i) = 0.0
      DO k = k1 , k2
         Y(i) = Y(i) + A(k)*X(Ja(k))
      ENDDO
   ENDDO
END SUBROUTINE ope
!*==opet.f90 processed by SPAG 8.04RA 18:15 22 Feb 2024
!!SPAG Open source Personal, Educational or Academic User Clemson University  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
!
SUBROUTINE opet(N,X,Y,A,Ja,Ia)
   USE ISO_FORTRAN_ENV                 
   IMPLICIT NONE
!
! Dummy argument declarations rewritten by SPAG
!
   INTEGER , INTENT(IN) :: N
   REAL(REAL64) , INTENT(IN) , DIMENSION(1) :: X
   REAL(REAL64) , INTENT(INOUT) , DIMENSION(1) :: Y
   REAL(REAL64) , INTENT(IN) , DIMENSION(1) :: A
   INTEGER , INTENT(IN) , DIMENSION(*) :: Ja
   INTEGER , INTENT(IN) , DIMENSION(*) :: Ia
!
! Local variable declarations rewritten by SPAG
!
   INTEGER :: i , j , k
!
! End of declarations rewritten by SPAG
!
! sparse matrix * vector multiplication
!
   DO j = 1 , N
      Y(j) = 0.0D0
   ENDDO
!
   DO i = 1 , N
      DO k = Ia(i) , Ia(i+1) - 1
         Y(Ja(k)) = Y(Ja(k)) + X(i)*A(k)
      ENDDO
   ENDDO
END SUBROUTINE opet
!*==ydfnorm.f90 processed by SPAG 8.04RA 18:15 22 Feb 2024
!!SPAG Open source Personal, Educational or Academic User Clemson University  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
!
SUBROUTINE ydfnorm(N,Y1,Y,Iout)
   USE ISO_FORTRAN_ENV                 
   IMPLICIT NONE
!
! Dummy argument declarations rewritten by SPAG
!
   INTEGER , INTENT(IN) :: N
   REAL(REAL64) , INTENT(IN) , DIMENSION(*) :: Y1
   REAL(REAL64) , INTENT(IN) , DIMENSION(*) :: Y
   INTEGER , INTENT(IN) :: Iout
!
! Local variable declarations rewritten by SPAG
!
   INTEGER :: k
   REAL(REAL64) :: t
!
! End of declarations rewritten by SPAG
!
!
   t = 0.0D0
   DO k = 1 , N
      t = t + (Y(k)-Y1(k))**2
   ENDDO
   t = sqrt(t)
   WRITE (Iout,*) '2-norm of error (exact answer-tested answer)=' , t
!-----------------------------------------------------------------------
END SUBROUTINE ydfnorm
