!*==dcopy.f90 processed by SPAG 8.04RA 00:18 23 Feb 2024
!!SPAG Open source Personal, Educational or Academic User Clemson University  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
SUBROUTINE dcopy(N,Dx,Incx,Dy,Incy)
   USE ISO_FORTRAN_ENV                 
   IMPLICIT NONE
!
! Dummy argument declarations rewritten by SPAG
!
   INTEGER , INTENT(IN) :: N
   REAL(REAL64) , INTENT(IN) , DIMENSION(1) :: Dx
   INTEGER , INTENT(IN) :: Incx
   REAL(REAL64) , INTENT(OUT) , DIMENSION(1) :: Dy
   INTEGER , INTENT(IN) :: Incy
!
! Local variable declarations rewritten by SPAG
!
   INTEGER :: i , ix , iy , m , mp1
!
! End of declarations rewritten by SPAG
!
!
!     copies a vector, x, to a vector, y.
!     uses unrolled loops for increments equal to one.
!     jack dongarra, linpack, 3/11/78.
!
!
   IF ( N<=0 ) RETURN
   IF ( Incx==1 .AND. Incy==1 ) THEN
!
!        code for both increments equal to 1
!
!
!        clean-up loop
!
      m = mod(N,7)
      IF ( m/=0 ) THEN
         DO i = 1 , m
            Dy(i) = Dx(i)
         ENDDO
         IF ( N<7 ) RETURN
      ENDIF
      mp1 = m + 1
      DO i = mp1 , N , 7
         Dy(i) = Dx(i)
         Dy(i+1) = Dx(i+1)
         Dy(i+2) = Dx(i+2)
         Dy(i+3) = Dx(i+3)
         Dy(i+4) = Dx(i+4)
         Dy(i+5) = Dx(i+5)
         Dy(i+6) = Dx(i+6)
      ENDDO
   ELSE
!
!        code for unequal increments or equal increments
!          not equal to 1
!
      ix = 1
      iy = 1
      IF ( Incx<0 ) ix = (-N+1)*Incx + 1
      IF ( Incy<0 ) iy = (-N+1)*Incy + 1
      DO i = 1 , N
         Dy(iy) = Dx(ix)
         ix = ix + Incx
         iy = iy + Incy
      ENDDO
      RETURN
   ENDIF
END SUBROUTINE dcopy
!*==ddot.f90 processed by SPAG 8.04RA 00:18 23 Feb 2024
!!SPAG Open source Personal, Educational or Academic User Clemson University  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
 
FUNCTION ddot(N,Dx,Incx,Dy,Incy)
   USE ISO_FORTRAN_ENV                 
   IMPLICIT NONE
!
! Function and Dummy argument declarations rewritten by SPAG
!
   REAL(REAL64) :: ddot
   INTEGER , INTENT(IN) :: N
   REAL(REAL64) , INTENT(IN) , DIMENSION(1) :: Dx
   INTEGER , INTENT(IN) :: Incx
   REAL(REAL64) , INTENT(IN) , DIMENSION(1) :: Dy
   INTEGER , INTENT(IN) :: Incy
!
! Local variable declarations rewritten by SPAG
!
   REAL(REAL64) :: dtemp
   INTEGER :: i , ix , iy , m , mp1
!
! End of declarations rewritten by SPAG
!
!
!     forms the dot product of two vectors.
!     uses unrolled loops for increments equal to one.
!     jack dongarra, linpack, 3/11/78.
!
!
   ddot = 0.0D0
   dtemp = 0.0D0
   IF ( N<=0 ) RETURN
   IF ( Incx==1 .AND. Incy==1 ) THEN
!
!        code for both increments equal to 1
!
!
!        clean-up loop
!
      m = mod(N,5)
      IF ( m/=0 ) THEN
         DO i = 1 , m
            dtemp = dtemp + Dx(i)*Dy(i)
         ENDDO
         IF ( N<5 ) THEN
            ddot = dtemp
            RETURN
         ENDIF
      ENDIF
      mp1 = m + 1
      DO i = mp1 , N , 5
         dtemp = dtemp + Dx(i)*Dy(i) + Dx(i+1)*Dy(i+1) + Dx(i+2)*Dy(i+2) + Dx(i+3)*Dy(i+3) + Dx(i+4)*Dy(i+4)
      ENDDO
      ddot = dtemp
   ELSE
!
!        code for unequal increments or equal increments
!          not equal to 1
!
      ix = 1
      iy = 1
      IF ( Incx<0 ) ix = (-N+1)*Incx + 1
      IF ( Incy<0 ) iy = (-N+1)*Incy + 1
      DO i = 1 , N
         dtemp = dtemp + Dx(ix)*Dy(iy)
         ix = ix + Incx
         iy = iy + Incy
      ENDDO
      ddot = dtemp
      RETURN
   ENDIF
END FUNCTION ddot
!*==dasum.f90 processed by SPAG 8.04RA 00:18 23 Feb 2024
!!SPAG Open source Personal, Educational or Academic User Clemson University  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
!
FUNCTION dasum(N,Dx,Incx)
   USE ISO_FORTRAN_ENV                 
   IMPLICIT NONE
!
! Function and Dummy argument declarations rewritten by SPAG
!
   REAL(REAL64) :: dasum
   INTEGER , INTENT(IN) :: N
   REAL(REAL64) , INTENT(IN) , DIMENSION(1) :: Dx
   INTEGER , INTENT(IN) :: Incx
!
! Local variable declarations rewritten by SPAG
!
   REAL(REAL64) :: dtemp
   INTEGER :: i , m , mp1 , nincx
!
! End of declarations rewritten by SPAG
!
!
!     takes the sum of the absolute values.
!     jack dongarra, linpack, 3/11/78.
!
!
   dasum = 0.0D0
   dtemp = 0.0D0
   IF ( N<=0 ) RETURN
   IF ( Incx==1 ) THEN
!
!        code for increment equal to 1
!
!
!        clean-up loop
!
      m = mod(N,6)
      IF ( m/=0 ) THEN
         DO i = 1 , m
            dtemp = dtemp + dabs(Dx(i))
         ENDDO
         IF ( N<6 ) THEN
            dasum = dtemp
            RETURN
         ENDIF
      ENDIF
      mp1 = m + 1
      DO i = mp1 , N , 6
         dtemp = dtemp + dabs(Dx(i)) + dabs(Dx(i+1)) + dabs(Dx(i+2)) + dabs(Dx(i+3)) + dabs(Dx(i+4)) + dabs(Dx(i+5))
      ENDDO
      dasum = dtemp
   ELSE
!
!        code for increment not equal to 1
!
      nincx = N*Incx
      DO i = 1 , nincx , Incx
         dtemp = dtemp + dabs(Dx(i))
      ENDDO
      dasum = dtemp
      RETURN
   ENDIF
END FUNCTION dasum
!*==daxpy.f90 processed by SPAG 8.04RA 00:18 23 Feb 2024
!!SPAG Open source Personal, Educational or Academic User Clemson University  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
 
SUBROUTINE daxpy(N,Da,Dx,Incx,Dy,Incy)
   USE ISO_FORTRAN_ENV                 
   IMPLICIT NONE
!
! Dummy argument declarations rewritten by SPAG
!
   INTEGER , INTENT(IN) :: N
   REAL(REAL64) , INTENT(IN) :: Da
   REAL(REAL64) , INTENT(IN) , DIMENSION(1) :: Dx
   INTEGER , INTENT(IN) :: Incx
   REAL(REAL64) , INTENT(INOUT) , DIMENSION(1) :: Dy
   INTEGER , INTENT(IN) :: Incy
!
! Local variable declarations rewritten by SPAG
!
   INTEGER :: i , ix , iy , m , mp1
!
! End of declarations rewritten by SPAG
!
!
!     constant times a vector plus a vector.
!     uses unrolled loops for increments equal to one.
!     jack dongarra, linpack, 3/11/78.
!
!
   IF ( N<=0 ) RETURN
   IF ( Da==0.0D0 ) RETURN
   IF ( Incx==1 .AND. Incy==1 ) THEN
!
!        code for both increments equal to 1
!
!
!        clean-up loop
!
      m = mod(N,4)
      IF ( m/=0 ) THEN
         DO i = 1 , m
            Dy(i) = Dy(i) + Da*Dx(i)
         ENDDO
         IF ( N<4 ) RETURN
      ENDIF
      mp1 = m + 1
      DO i = mp1 , N , 4
         Dy(i) = Dy(i) + Da*Dx(i)
         Dy(i+1) = Dy(i+1) + Da*Dx(i+1)
         Dy(i+2) = Dy(i+2) + Da*Dx(i+2)
         Dy(i+3) = Dy(i+3) + Da*Dx(i+3)
      ENDDO
   ELSE
!
!        code for unequal increments or equal increments
!          not equal to 1
!
      ix = 1
      iy = 1
      IF ( Incx<0 ) ix = (-N+1)*Incx + 1
      IF ( Incy<0 ) iy = (-N+1)*Incy + 1
      DO i = 1 , N
         Dy(iy) = Dy(iy) + Da*Dx(ix)
         ix = ix + Incx
         iy = iy + Incy
      ENDDO
      RETURN
   ENDIF
END SUBROUTINE daxpy
!*==dnrm2.f90 processed by SPAG 8.04RA 00:18 23 Feb 2024
!!SPAG Open source Personal, Educational or Academic User Clemson University  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
FUNCTION dnrm2(N,Dx,Incx)
   USE ISO_FORTRAN_ENV                 
   IMPLICIT NONE
!
! Function and Dummy argument declarations rewritten by SPAG
!
   REAL(REAL64) :: dnrm2
   INTEGER , INTENT(IN) :: N
   REAL(REAL64) , INTENT(IN) , DIMENSION(1) :: Dx
   INTEGER , INTENT(IN) :: Incx
!
! Local variable declarations rewritten by SPAG
!
   REAL(REAL64) , SAVE :: cuthi , cutlo , one , zero
   REAL(REAL64) :: hitest , sum , xmax
   INTEGER :: i , j , next , nn
!
! End of declarations rewritten by SPAG
!
   INTEGER :: spag_nextblock_1
   DATA zero , one/0.0D0 , 1.0D0/
!
!     euclidean norm of the n-vector stored in dx() with storage
!     increment incx .
!     if    n .le. 0 return with result = 0.
!     if n .ge. 1 then incx must be .ge. 1
!
!           c.l.lawson, 1978 jan 08
!
!     four phase method     using two built-in constants that are
!     hopefully applicable to all machines.
!         cutlo = maximum of  dsqrt(u/eps)  over all known machines.
!         cuthi = minimum of  dsqrt(v)      over all known machines.
!     where
!         eps = smallest no. such that eps + 1. .gt. 1.
!         u   = smallest positive no.   (underflow limit)
!         v   = largest  no.            (overflow  limit)
!
!     brief outline of algorithm..
!
!     phase 1    scans zero components.
!     move to phase 2 when a component is nonzero and .le. cutlo
!     move to phase 3 when a component is .gt. cutlo
!     move to phase 4 when a component is .ge. cuthi/m
!     where m = n for x() real and m = 2*n for complex.
!
!     values for cutlo and cuthi..
!     from the environmental parameters listed in the imsl converter
!     document the limiting values are as follows..
!     cutlo, s.p.   u/eps = 2**(-102) for  honeywell.  close seconds are
!                   univac and dec at 2**(-103)
!                   thus cutlo = 2**(-51) = 4.44089e-16
!     cuthi, s.p.   v = 2**127 for univac, honeywell, and dec.
!                   thus cuthi = 2**(63.5) = 1.30438e19
!     cutlo, d.p.   u/eps = 2**(-67) for honeywell and dec.
!                   thus cutlo = 2**(-33.5) = 8.23181d-11
!     cuthi, d.p.   same as s.p.  cuthi = 1.30438d19
!     data cutlo, cuthi / 8.232d-11,  1.304d19 /
!     data cutlo, cuthi / 4.441e-16,  1.304e19 /
   DATA cutlo , cuthi/8.232D-11 , 1.304D19/
   spag_nextblock_1 = 1
   SPAG_DispatchLoop_1: DO
      SELECT CASE (spag_nextblock_1)
      CASE (1)
!
         IF ( N>0 ) THEN
!
            ASSIGN 20 TO next
            sum = zero
            nn = N*Incx
!                                                 begin main loop
            i = 1
         ELSE
            dnrm2 = zero
            EXIT SPAG_DispatchLoop_1
         ENDIF
         spag_nextblock_1 = 2
      CASE (2)
         GOTO next
 20      IF ( dabs(Dx(i))>cutlo ) THEN
            spag_nextblock_1 = 5
            CYCLE SPAG_DispatchLoop_1
         ENDIF
         ASSIGN 40 TO next
         xmax = zero
!
!                        phase 1.  sum is zero
!
 40      IF ( Dx(i)==zero ) THEN
            spag_nextblock_1 = 6
            CYCLE SPAG_DispatchLoop_1
         ENDIF
         IF ( dabs(Dx(i))>cutlo ) THEN
            spag_nextblock_1 = 5
            CYCLE SPAG_DispatchLoop_1
         ENDIF
!
!                                prepare for phase 2.
         ASSIGN 60 TO next
         spag_nextblock_1 = 4
         CYCLE SPAG_DispatchLoop_1
      CASE (3)
!
!                                prepare for phase 4.
!
         i = j
         ASSIGN 80 TO next
         sum = (sum/Dx(i))/Dx(i)
         spag_nextblock_1 = 4
      CASE (4)
         xmax = dabs(Dx(i))
!
         sum = sum + (Dx(i)/xmax)**2
         spag_nextblock_1 = 6
         CYCLE SPAG_DispatchLoop_1
!
!                   phase 2.  sum is small.
!                             scale to avoid destructive underflow.
!
 60      IF ( dabs(Dx(i))>cutlo ) THEN
!
!
!                  prepare for phase 3.
!
            sum = (sum*xmax)*xmax
            spag_nextblock_1 = 5
            CYCLE SPAG_DispatchLoop_1
         ENDIF
!
!                     common code for phases 2 and 4.
!                     in phase 4 sum is large.  scale to avoid overflow.
!
 80      IF ( dabs(Dx(i))<=xmax ) THEN
            sum = sum + (Dx(i)/xmax)**2
         ELSE
            sum = one + sum*(xmax/Dx(i))**2
            xmax = dabs(Dx(i))
         ENDIF
         spag_nextblock_1 = 6
         CYCLE SPAG_DispatchLoop_1
      CASE (5)
!
!
!     for real or d.p. set hitest = cuthi/n
!     for complex      set hitest = cuthi/(2*n)
!
         hitest = cuthi/float(N)
!
!                   phase 3.  sum is mid-range.  no scaling.
!
         DO j = i , nn , Incx
            IF ( dabs(Dx(j))>=hitest ) THEN
               spag_nextblock_1 = 3
               CYCLE SPAG_DispatchLoop_1
            ENDIF
            sum = sum + Dx(j)**2
         ENDDO
         dnrm2 = dsqrt(sum)
         EXIT SPAG_DispatchLoop_1
      CASE (6)
!
         i = i + Incx
         IF ( i<=nn ) THEN
            spag_nextblock_1 = 2
            CYCLE SPAG_DispatchLoop_1
         ENDIF
!
!              end of main loop.
!
!              compute square root and adjust for scaling.
!
         dnrm2 = xmax*dsqrt(sum)
         EXIT SPAG_DispatchLoop_1
      END SELECT
   ENDDO SPAG_DispatchLoop_1
END FUNCTION dnrm2
!*==dscal.f90 processed by SPAG 8.04RA 00:18 23 Feb 2024
!!SPAG Open source Personal, Educational or Academic User Clemson University  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
 
SUBROUTINE dscal(N,Da,Dx,Incx)
   USE ISO_FORTRAN_ENV                 
   IMPLICIT NONE
!
! Dummy argument declarations rewritten by SPAG
!
   INTEGER , INTENT(IN) :: N
   REAL(REAL64) , INTENT(IN) :: Da
   REAL(REAL64) , INTENT(INOUT) , DIMENSION(1) :: Dx
   INTEGER , INTENT(IN) :: Incx
!
! Local variable declarations rewritten by SPAG
!
   INTEGER :: i , m , mp1 , nincx
!
! End of declarations rewritten by SPAG
!
!     scales a vector by a constant.
!     uses unrolled loops for increment equal to one.
!     jack dongarra, linpack, 3/11/78.
!
!
   IF ( N<=0 ) RETURN
   IF ( Incx==1 ) THEN
!
!        code for increment equal to 1
!
!
!        clean-up loop
!
      m = mod(N,5)
      IF ( m/=0 ) THEN
         DO i = 1 , m
            Dx(i) = Da*Dx(i)
         ENDDO
         IF ( N<5 ) RETURN
      ENDIF
      mp1 = m + 1
      DO i = mp1 , N , 5
         Dx(i) = Da*Dx(i)
         Dx(i+1) = Da*Dx(i+1)
         Dx(i+2) = Da*Dx(i+2)
         Dx(i+3) = Da*Dx(i+3)
         Dx(i+4) = Da*Dx(i+4)
      ENDDO
   ELSE
!
!        code for increment not equal to 1
!
      nincx = N*Incx
      DO i = 1 , nincx , Incx
         Dx(i) = Da*Dx(i)
      ENDDO
      RETURN
   ENDIF
END SUBROUTINE dscal
!*==dswap.f90 processed by SPAG 8.04RA 00:18 23 Feb 2024
!!SPAG Open source Personal, Educational or Academic User Clemson University  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
 
SUBROUTINE dswap(N,Dx,Incx,Dy,Incy)
   USE ISO_FORTRAN_ENV                 
   IMPLICIT NONE
!
! Dummy argument declarations rewritten by SPAG
!
   INTEGER , INTENT(IN) :: N
   REAL(REAL64) , INTENT(INOUT) , DIMENSION(1) :: Dx
   INTEGER , INTENT(IN) :: Incx
   REAL(REAL64) , INTENT(INOUT) , DIMENSION(1) :: Dy
   INTEGER , INTENT(IN) :: Incy
!
! Local variable declarations rewritten by SPAG
!
   REAL(REAL64) :: dtemp
   INTEGER :: i , ix , iy , m , mp1
!
! End of declarations rewritten by SPAG
!
!
!     interchanges two vectors.
!     uses unrolled loops for increments equal one.
!     jack dongarra, linpack, 3/11/78.
!
!
   IF ( N<=0 ) RETURN
   IF ( Incx==1 .AND. Incy==1 ) THEN
!
!       code for both increments equal to 1
!
!
!       clean-up loop
!
      m = mod(N,3)
      IF ( m/=0 ) THEN
         DO i = 1 , m
            dtemp = Dx(i)
            Dx(i) = Dy(i)
            Dy(i) = dtemp
         ENDDO
         IF ( N<3 ) RETURN
      ENDIF
      mp1 = m + 1
      DO i = mp1 , N , 3
         dtemp = Dx(i)
         Dx(i) = Dy(i)
         Dy(i) = dtemp
         dtemp = Dx(i+1)
         Dx(i+1) = Dy(i+1)
         Dy(i+1) = dtemp
         dtemp = Dx(i+2)
         Dx(i+2) = Dy(i+2)
         Dy(i+2) = dtemp
      ENDDO
   ELSE
!
!       code for unequal increments or equal increments not equal
!         to 1
!
      ix = 1
      iy = 1
      IF ( Incx<0 ) ix = (-N+1)*Incx + 1
      IF ( Incy<0 ) iy = (-N+1)*Incy + 1
      DO i = 1 , N
         dtemp = Dx(ix)
         Dx(ix) = Dy(iy)
         Dy(iy) = dtemp
         ix = ix + Incx
         iy = iy + Incy
      ENDDO
      RETURN
   ENDIF
END SUBROUTINE dswap
!*==idamax.f90 processed by SPAG 8.04RA 00:18 23 Feb 2024
!!SPAG Open source Personal, Educational or Academic User Clemson University  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
 
FUNCTION idamax(N,Dx,Incx)
   USE ISO_FORTRAN_ENV                 
   IMPLICIT NONE
!
! Function and Dummy argument declarations rewritten by SPAG
!
   INTEGER :: idamax
   INTEGER , INTENT(IN) :: N
   REAL(REAL64) , INTENT(IN) , DIMENSION(1) :: Dx
   INTEGER , INTENT(IN) :: Incx
!
! Local variable declarations rewritten by SPAG
!
   REAL(REAL64) :: dmax
   INTEGER :: i , ix
!
! End of declarations rewritten by SPAG
!
!
!     finds the index of element having max. absolute value.
!     jack dongarra, linpack, 3/11/78.
!
!
   idamax = 0
   IF ( N<1 ) RETURN
   idamax = 1
   IF ( N==1 ) RETURN
   IF ( Incx==1 ) THEN
!
!        code for increment equal to 1
!
      dmax = dabs(Dx(1))
      DO i = 2 , N
         IF ( dabs(Dx(i))>dmax ) THEN
            idamax = i
            dmax = dabs(Dx(i))
         ENDIF
      ENDDO
      RETURN
   ENDIF
!
!        code for increment not equal to 1
!
   ix = 1
   dmax = dabs(Dx(1))
   ix = ix + Incx
   DO i = 2 , N
      IF ( dabs(Dx(ix))>dmax ) THEN
         idamax = i
         dmax = dabs(Dx(ix))
      ENDIF
      ix = ix + Incx
   ENDDO
   RETURN
END FUNCTION idamax
!*==drot.f90 processed by SPAG 8.04RA 00:18 23 Feb 2024
!!SPAG Open source Personal, Educational or Academic User Clemson University  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
!
SUBROUTINE drot(N,Dx,Incx,Dy,Incy,C,S)
   USE ISO_FORTRAN_ENV                 
   IMPLICIT NONE
!
! Dummy argument declarations rewritten by SPAG
!
   INTEGER , INTENT(IN) :: N
   REAL(REAL64) , INTENT(INOUT) , DIMENSION(1) :: Dx
   INTEGER , INTENT(IN) :: Incx
   REAL(REAL64) , INTENT(INOUT) , DIMENSION(1) :: Dy
   INTEGER , INTENT(IN) :: Incy
   REAL(REAL64) , INTENT(IN) :: C
   REAL(REAL64) , INTENT(IN) :: S
!
! Local variable declarations rewritten by SPAG
!
   REAL(REAL64) :: dtemp
   INTEGER :: i , ix , iy
!
! End of declarations rewritten by SPAG
!
!
!     applies a plane rotation.
!     jack dongarra, linpack, 3/11/78.
!
!
   IF ( N<=0 ) RETURN
   IF ( Incx==1 .AND. Incy==1 ) THEN
!
!       code for both increments equal to 1
!
      DO i = 1 , N
         dtemp = C*Dx(i) + S*Dy(i)
         Dy(i) = C*Dy(i) - S*Dx(i)
         Dx(i) = dtemp
      ENDDO
      RETURN
   ENDIF
!
!       code for unequal increments or equal increments not equal
!         to 1
!
   ix = 1
   iy = 1
   IF ( Incx<0 ) ix = (-N+1)*Incx + 1
   IF ( Incy<0 ) iy = (-N+1)*Incy + 1
   DO i = 1 , N
      dtemp = C*Dx(ix) + S*Dy(iy)
      Dy(iy) = C*Dy(iy) - S*Dx(ix)
      Dx(ix) = dtemp
      ix = ix + Incx
      iy = iy + Incy
   ENDDO
   RETURN
END SUBROUTINE drot
!*==drotg.f90 processed by SPAG 8.04RA 00:18 23 Feb 2024
!!SPAG Open source Personal, Educational or Academic User Clemson University  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
!
SUBROUTINE drotg(Da,Db,C,S)
   USE ISO_FORTRAN_ENV                 
   IMPLICIT NONE
!
! Dummy argument declarations rewritten by SPAG
!
   REAL(REAL64) , INTENT(INOUT) :: Da
   REAL(REAL64) , INTENT(INOUT) :: Db
   REAL(REAL64) , INTENT(INOUT) :: C
   REAL(REAL64) , INTENT(INOUT) :: S
!
! Local variable declarations rewritten by SPAG
!
   REAL(REAL64) :: r , roe , scale , z
!
! End of declarations rewritten by SPAG
!
!
!     construct givens plane rotation.
!     jack dongarra, linpack, 3/11/78.
!
!
   roe = Db
   IF ( dabs(Da)>dabs(Db) ) roe = Da
   scale = dabs(Da) + dabs(Db)
   IF ( scale/=0.0D0 ) THEN
      r = scale*dsqrt((Da/scale)**2+(Db/scale)**2)
      r = dsign(1.0D0,roe)*r
      C = Da/r
      S = Db/r
   ELSE
      C = 1.0D0
      S = 0.0D0
      r = 0.0D0
   ENDIF
   z = 1.0D0
   IF ( dabs(Da)>dabs(Db) ) z = S
   IF ( dabs(Db)>=dabs(Da) .AND. C/=0.0D0 ) z = 1.0D0/C
   Da = r
   Db = z
END SUBROUTINE drotg
!*==ccopy.f90 processed by SPAG 8.04RA 00:18 23 Feb 2024
!!SPAG Open source Personal, Educational or Academic User Clemson University  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
!
SUBROUTINE ccopy(N,Cx,Incx,Cy,Incy)
   IMPLICIT NONE
!
! Dummy argument declarations rewritten by SPAG
!
   INTEGER , INTENT(IN) :: N
   COMPLEX , INTENT(IN) , DIMENSION(1) :: Cx
   INTEGER , INTENT(IN) :: Incx
   COMPLEX , INTENT(OUT) , DIMENSION(1) :: Cy
   INTEGER , INTENT(IN) :: Incy
!
! Local variable declarations rewritten by SPAG
!
   INTEGER :: i , ix , iy
!
! End of declarations rewritten by SPAG
!
!
!     copies a vector, x, to a vector, y.
!     jack dongarra, linpack, 3/11/78.
!
!
   IF ( N<=0 ) RETURN
   IF ( Incx==1 .AND. Incy==1 ) THEN
!
!        code for both increments equal to 1
!
      DO i = 1 , N
         Cy(i) = Cx(i)
      ENDDO
      RETURN
   ENDIF
!
!        code for unequal increments or equal increments
!          not equal to 1
!
   ix = 1
   iy = 1
   IF ( Incx<0 ) ix = (-N+1)*Incx + 1
   IF ( Incy<0 ) iy = (-N+1)*Incy + 1
   DO i = 1 , N
      Cy(iy) = Cx(ix)
      ix = ix + Incx
      iy = iy + Incy
   ENDDO
   RETURN
END SUBROUTINE ccopy
!*==cscal.f90 processed by SPAG 8.04RA 00:18 23 Feb 2024
!!SPAG Open source Personal, Educational or Academic User Clemson University  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
SUBROUTINE cscal(N,Ca,Cx,Incx)
   IMPLICIT NONE
!
! Dummy argument declarations rewritten by SPAG
!
   INTEGER , INTENT(IN) :: N
   COMPLEX , INTENT(IN) :: Ca
   COMPLEX , INTENT(INOUT) , DIMENSION(1) :: Cx
   INTEGER , INTENT(IN) :: Incx
!
! Local variable declarations rewritten by SPAG
!
   INTEGER :: i , nincx
!
! End of declarations rewritten by SPAG
!
!
!     scales a vector by a constant.
!     jack dongarra, linpack,  3/11/78.
!
!
   IF ( N<=0 ) RETURN
   IF ( Incx==1 ) THEN
!
!        code for increment equal to 1
!
      DO i = 1 , N
         Cx(i) = Ca*Cx(i)
      ENDDO
      RETURN
   ENDIF
!
!        code for increment not equal to 1
!
   nincx = N*Incx
   DO i = 1 , nincx , Incx
      Cx(i) = Ca*Cx(i)
   ENDDO
   RETURN
END SUBROUTINE cscal
!*==csrot.f90 processed by SPAG 8.04RA 00:18 23 Feb 2024
!!SPAG Open source Personal, Educational or Academic User Clemson University  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
!
SUBROUTINE csrot(N,Cx,Incx,Cy,Incy,C,S)
   IMPLICIT NONE
!
! Dummy argument declarations rewritten by SPAG
!
   INTEGER , INTENT(IN) :: N
   COMPLEX , INTENT(INOUT) , DIMENSION(1) :: Cx
   INTEGER , INTENT(IN) :: Incx
   COMPLEX , INTENT(INOUT) , DIMENSION(1) :: Cy
   INTEGER , INTENT(IN) :: Incy
   REAL , INTENT(IN) :: C
   REAL , INTENT(IN) :: S
!
! Local variable declarations rewritten by SPAG
!
   COMPLEX :: ctemp
   INTEGER :: i , ix , iy
!
! End of declarations rewritten by SPAG
!
!
!     applies a plane rotation, where the cos and sin (c and s) are real
!     and the vectors cx and cy are complex.
!     jack dongarra, linpack, 3/11/78.
!
!
   IF ( N<=0 ) RETURN
   IF ( Incx==1 .AND. Incy==1 ) THEN
!
!       code for both increments equal to 1
!
      DO i = 1 , N
         ctemp = C*Cx(i) + S*Cy(i)
         Cy(i) = C*Cy(i) - S*Cx(i)
         Cx(i) = ctemp
      ENDDO
      RETURN
   ENDIF
!
!       code for unequal increments or equal increments not equal
!         to 1
!
   ix = 1
   iy = 1
   IF ( Incx<0 ) ix = (-N+1)*Incx + 1
   IF ( Incy<0 ) iy = (-N+1)*Incy + 1
   DO i = 1 , N
      ctemp = C*Cx(ix) + S*Cy(iy)
      Cy(iy) = C*Cy(iy) - S*Cx(ix)
      Cx(ix) = ctemp
      ix = ix + Incx
      iy = iy + Incy
   ENDDO
   RETURN
END SUBROUTINE csrot
!*==cswap.f90 processed by SPAG 8.04RA 00:18 23 Feb 2024
!!SPAG Open source Personal, Educational or Academic User Clemson University  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
SUBROUTINE cswap(N,Cx,Incx,Cy,Incy)
   IMPLICIT NONE
!
! Dummy argument declarations rewritten by SPAG
!
   INTEGER , INTENT(IN) :: N
   COMPLEX , INTENT(INOUT) , DIMENSION(1) :: Cx
   INTEGER , INTENT(IN) :: Incx
   COMPLEX , INTENT(INOUT) , DIMENSION(1) :: Cy
   INTEGER , INTENT(IN) :: Incy
!
! Local variable declarations rewritten by SPAG
!
   COMPLEX :: ctemp
   INTEGER :: i , ix , iy
!
! End of declarations rewritten by SPAG
!
!
!     interchanges two vectors.
!     jack dongarra, linpack, 3/11/78.
!
!
   IF ( N<=0 ) RETURN
   IF ( Incx==1 .AND. Incy==1 ) THEN
!
!       code for both increments equal to 1
      DO i = 1 , N
         ctemp = Cx(i)
         Cx(i) = Cy(i)
         Cy(i) = ctemp
      ENDDO
      RETURN
   ENDIF
!
!       code for unequal increments or equal increments not equal
!         to 1
!
   ix = 1
   iy = 1
   IF ( Incx<0 ) ix = (-N+1)*Incx + 1
   IF ( Incy<0 ) iy = (-N+1)*Incy + 1
   DO i = 1 , N
      ctemp = Cx(ix)
      Cx(ix) = Cy(iy)
      Cy(iy) = ctemp
      ix = ix + Incx
      iy = iy + Incy
   ENDDO
   RETURN
END SUBROUTINE cswap
!*==csscal.f90 processed by SPAG 8.04RA 00:18 23 Feb 2024
!!SPAG Open source Personal, Educational or Academic User Clemson University  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
SUBROUTINE csscal(N,Sa,Cx,Incx)
   IMPLICIT NONE
!
! Dummy argument declarations rewritten by SPAG
!
   INTEGER , INTENT(IN) :: N
   REAL , INTENT(IN) :: Sa
   COMPLEX , INTENT(INOUT) , DIMENSION(1) :: Cx
   INTEGER , INTENT(IN) :: Incx
!
! Local variable declarations rewritten by SPAG
!
   INTEGER :: i , nincx
!
! End of declarations rewritten by SPAG
!
!
!     scales a complex vector by a real constant.
!     jack dongarra, linpack, 3/11/78.
!
!
   IF ( N<=0 ) RETURN
   IF ( Incx==1 ) THEN
!
!        code for increment equal to 1
!
      DO i = 1 , N
         Cx(i) = cmplx(Sa*real(Cx(i)),Sa*aimag(Cx(i)))
      ENDDO
      RETURN
   ENDIF
!
!        code for increment not equal to 1
!
   nincx = N*Incx
   DO i = 1 , nincx , Incx
      Cx(i) = cmplx(Sa*real(Cx(i)),Sa*aimag(Cx(i)))
   ENDDO
   RETURN
END SUBROUTINE csscal
