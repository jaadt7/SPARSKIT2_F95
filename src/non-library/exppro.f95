!*==expprod.f90 processed by SPAG 8.04RA 12:07 23 Feb 2024
!!SPAG Open source Personal, Educational or Academic User Clemson University  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
SUBROUTINE expprod(N,M,Eps,Tn,U,W,X,Y,A,Ioff,Ndiag,Iverboz)
   USE ISO_FORTRAN_ENV                 
   IMPLICIT NONE
!
! Dummy argument declarations rewritten by SPAG
!
   INTEGER :: N
   INTEGER :: M
   INTEGER :: Ndiag
   REAL(REAL64) :: Eps
   REAL(REAL64) :: Tn
   REAL(REAL64) , DIMENSION(N,M+1) :: U
   REAL(REAL64) , DIMENSION(N) :: W
   REAL(REAL64) , DIMENSION(N) :: X
   REAL(REAL64) , DIMENSION(N) :: Y
   REAL(REAL64) , DIMENSION(N,Ndiag) :: A
   INTEGER , DIMENSION(Ndiag) :: Ioff
   INTEGER , INTENT(IN) :: Iverboz
!
! Local variable declarations rewritten by SPAG
!
   INTEGER :: ierr , indic
   EXTERNAL exppro , oped
!
! End of declarations rewritten by SPAG
!
 
!-----------------------------------------------------------------------
! this subroutine computes an approximation to the vector
!
!        	w :=  exp( - A * tn ) * w
!
! for matrices stored in diagonal (DIA) format.
!
! this routine constitutes an interface for the routine exppro for
! matrices stored in diagonal (DIA) format.
!-----------------------------------------------------------------------
! ARGUMENTS
!----------
! see exppro for meaning of parameters n, m, eps, tn, u, w, x, y.
!
! a, ioff, and ndiag are the arguments of the matrix:
!
! a(n,ndiag) = a rectangular array with a(*,k) containing the diagonal
!              offset by ioff(k) (negative or positive or zero), i.e.,
!              a(i,jdiag) contains the element A(i,i+ioff(jdiag)) in
!              the usual dense storage scheme.
!
! ioff	     = integer array containing the offsets  of the ndiag diagonals
! ndiag      = integer. the number of diagonals.
! Iverboz    = integer. flag to print out diagnositics (1) or not (0)
!
!-----------------------------------------------------------------------
! local variables
!
   indic = 0
   SPAG_Loop_1_1: DO
      CALL exppro(N,M,Eps,Tn,U,W,X,Y,indic,ierr,Iverboz)
      IF ( indic==1 ) EXIT SPAG_Loop_1_1
!
!     matrix vector-product for diagonal storage --
!
      CALL oped(N,X,Y,A,Ioff,Ndiag)
   ENDDO SPAG_Loop_1_1
END SUBROUTINE expprod
!*==exppro.f90 processed by SPAG 8.04RA 12:07 23 Feb 2024
!!SPAG Open source Personal, Educational or Academic User Clemson University  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
!----------end-of-expprod-----------------------------------------------
!-----------------------------------------------------------------------
SUBROUTINE exppro(N,M,Eps,Tn,U,W,X,Y,Indic,Ierr,Iverboz)
   USE ISO_FORTRAN_ENV                 
   IMPLICIT NONE
!
! PARAMETER definitions rewritten by SPAG
!
   INTEGER , PARAMETER :: MMAX = 60
!
! Dummy argument declarations rewritten by SPAG
!
   INTEGER :: N
   INTEGER , INTENT(INOUT) :: M
   REAL(REAL64) :: Eps
   REAL(REAL64) , INTENT(IN) :: Tn
   REAL(REAL64) , DIMENSION(N,M+1) :: U
   REAL(REAL64) , DIMENSION(N) :: W
   REAL(REAL64) , DIMENSION(N) :: X
   REAL(REAL64) , DIMENSION(N) :: Y
   INTEGER , INTENT(INOUT) :: Indic
   INTEGER , INTENT(INOUT) :: Ierr
   INTEGER , INTENT(IN) :: Iverboz
!
! Local variable declarations rewritten by SPAG
!
   REAL(REAL64) , SAVE :: beta , dtl , errst , red , tcur , told
   REAL(REAL64) , DIMENSION(MMAX+2,MMAX+1) , SAVE :: hh
   INTEGER , SAVE :: ih , job
   LOGICAL :: verboz
   COMPLEX(REAL64) , DIMENSION(MMAX+1) , SAVE :: wkc
   REAL(REAL64) , DIMENSION(MMAX+1) , SAVE :: z
   EXTERNAL exphes , exppro_project
!
! End of declarations rewritten by SPAG
!
!     implicit  real*8 (a-h,o-z)
!-----------------------------------------------------------------------
!
! this subroutine computes an approximation to the vector
!
!        	w :=  exp( - A * tn ) * w
!
! where A is an arbitary matrix and w is a given input vector
! uses a dynamic estimation of internal time advancement (dt)
!-----------------------------------------------------------------------
! THIS IS A REVERSE COMMUNICATION IMPLEMENTATION.
!-------------------------------------------------
! USAGE: (see also comments on indic below).
!------
!
!      indic = 0
! 1    continue
!      call exppro (n, m, eps, tn, u, w, x, y, indic)
!      if (indic .eq. 1) goto 2 <-- indic .eq. 1 means job is finished
!      call matvec(n, x, y)     <--- user's matrix-vec. product
!                                    with x = input vector, and
!                                     y = result = A * x.
!      goto 1
! 2    continue
!      .....
!
!-----------------------------------------------------------------------
!
! en entry:
!----------
! n	= dimension of matrix
!
! m	= dimension of Krylov subspace (= degree of polynomial
!         approximation to the exponential used. )
!
! eps   = scalar indicating the relative error tolerated for the result.
!         the code will try to compute an answer such that
!         norm2(exactanswer-approximation) / norm2(w) .le. eps
!
! tn	= scalar by which to multiply matrix. (may be .lt. 0)
!         the code will compute an approximation to exp(- tn * A) w
!         and overwrite the result onto w.
!
! u	= work array of size n*(m+1) (used to hold the Arnoldi basis )
!
! w	= real array of length n = input vector to  which exp(-A) is
!         to be applied. this is also an output argument
!
! x, y  = two real work vectors of length at least  n each.
!         see indic for usage.
!
! Indic = integer used as indicator for the reverse communication.
!         in the first call enter indic = 0. See below for more.
!
! Iverboz = integer flag to print out diagnostics (1) or not (0)
!
! on return:
!-----------
!
! w     = contains the resulting vector exp(-A * tn ) * w when
!         exppro has finished (see indic)
!
! indic = indicator for the reverse communication protocole.
!       * INDIC .eq. 1  means that exppro has finished and w contains the
!         result.
!       * INDIC .gt. 1 ,  means that exppro has not finished and that
!         it is requesting another matrix vector product before
!         continuing. The user must compute Ax where A is the matrix
!         and x is the vector provided by exppro, and return the
!         result in y. Then exppro must be called again without
!         changing any other argument. typically this must be
!         implemented in a loop with exppro being called as long
!         indic is returned with a value .ne. 1.
!
! ierr  = error indicator.
!         ierr = 1 means phipro was called with indic=1 (not allowed)
!         ierr = -1 means that the input is zero the solution has been
!         unchanged.
!
! NOTES:  m should not exceed 60 in this version  (see mmax below)
!-----------------------------------------------------------------------
! written by Y. Saad -- version feb, 1991.
!-----------------------------------------------------------------------
! For reference see following papers :
! (1) E. Gallopoulos and Y. Saad: Efficient solution of parabolic
!     equations by Krylov approximation methods. RIACS technical
!     report 90-14.
! (2) Y.Saad: Analysis of some Krylov subspace approximations to the
!     matrix exponential operator. RIACS Tech report. 90-14
!-----------------------------------------------------------------------

   INTEGER :: spag_nextblock_1
   spag_nextblock_1 = 1

   if (Iverboz .eq. 0) then
      verboz = .false.
   else
      verboz = .true.
   endif

   SPAG_DispatchLoop_1: DO
      SELECT CASE (spag_nextblock_1)
      CASE (1)
!-----------------------------------------------------------------------
! indic = 3  means  passing through only with result of y= Ax to exphes
! indic = 2  means exphes has finished its job
! indic = 1  means exppro has finished its job (real end)/
!-----------------------------------------------------------------------
         Ierr = 0
         IF ( Indic==3 ) THEN
            spag_nextblock_1 = 3
            CYCLE SPAG_DispatchLoop_1
         ENDIF
         IF ( Indic==1 ) THEN
            Ierr = 1
            RETURN
         ENDIF
!-----
         ih = MMAX
         M = min0(M,MMAX)
         tcur = 0.0D0
         dtl = Tn - tcur
         job = -1
         spag_nextblock_1 = 2
      CASE (2)
!-------------------- outer loop -----------------------------
         IF ( verboz ) PRINT * , 'In EXPPRO, current time = ' , tcur , '---------'
!-------------------------------------------------------------
! ---- call exponential propagator ---------------------------
!-------------------------------------------------------------
         told = tcur
         spag_nextblock_1 = 3
      CASE (3)
         DO
!     if (told + dtl .gt. tn) dtl = tn-told
            CALL exphes(N,M,dtl,Eps,U,W,job,z,wkc,beta,errst,hh,ih,X,Y,Indic,Ierr,verboz)
!-----------------------------------------------------------------------
            IF ( Ierr/=0 ) RETURN
            IF ( Indic>=3 ) RETURN
            tcur = told + dtl
            IF ( verboz ) PRINT * , ' tcur now = ' , tcur , ' dtl = ' , dtl
!
!     relative error
!      if(verboz) print *, ' beta', beta
            errst = errst/beta
!---------
            IF ( (errst<=Eps) .AND. ((errst>Eps/100.0) .OR. (tcur==Tn)) ) THEN
!-------
!
               CALL exppro_project(N,M,U,z,W)
! never go beyond tcur
               job = 0
               dtl = dmin1(dtl,Tn-tcur)
               IF ( dabs(tcur+dtl)>dabs(Tn) ) dtl = Tn - tcur
               IF ( dabs(tcur)<dabs(Tn) ) THEN
                  spag_nextblock_1 = 2
                  CYCLE SPAG_DispatchLoop_1
               ENDIF
               Indic = 1
               EXIT SPAG_DispatchLoop_1
            ELSE
!
!     use approximation :  [ new err ] = fact**m  * [cur. error]
!
               red = (0.5*Eps/errst)**(1.0D0/dble(M))
               dtl = dtl*red
               IF ( dabs(told+dtl)>dabs(Tn) ) dtl = Tn - told
               IF ( verboz ) PRINT * , ' red =' , red , ' , reducing dt to: ' , dtl
!-------
               job = 1
            ENDIF
         ENDDO
         EXIT SPAG_DispatchLoop_1
      END SELECT
   ENDDO SPAG_DispatchLoop_1
!
END SUBROUTINE exppro
!*==exphes.f90 processed by SPAG 8.04RA 12:07 23 Feb 2024
!!SPAG Open source Personal, Educational or Academic User Clemson University  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
!----------end-of-expro-------------------------------------------------
!-----------------------------------------------------------------------
SUBROUTINE exphes(N,M,Dt,Eps,U,W,Job,Z,Wkc,Beta,Errst,Hh,Ih,X,Y,Indic,Ierr,verboz)
   USE ISO_FORTRAN_ENV                 
   IMPLICIT NONE
!
! PARAMETER definitions rewritten by SPAG
!
   INTEGER , PARAMETER :: NDMAX = 20
!
! Dummy argument declarations rewritten by SPAG
!
   INTEGER :: N
   INTEGER , INTENT(INOUT) :: M
   INTEGER :: Ih
   REAL(REAL64) , INTENT(INOUT) :: Dt
   REAL(REAL64) , INTENT(IN) :: Eps
   REAL(REAL64) , INTENT(INOUT) , DIMENSION(N,M+1) :: U
   REAL(REAL64) , DIMENSION(N) :: W
   INTEGER , INTENT(IN) :: Job
   REAL(REAL64) , INTENT(INOUT) , DIMENSION(M+1) :: Z
   COMPLEX(REAL64) , DIMENSION(M+1) :: Wkc
   REAL(REAL64) , INTENT(INOUT) :: Beta
   REAL(REAL64) , INTENT(INOUT) :: Errst
   REAL(REAL64) , INTENT(INOUT) , DIMENSION(Ih+2,M+1) :: Hh
   REAL(REAL64) , INTENT(OUT) , DIMENSION(N) :: X
   REAL(REAL64) , INTENT(IN) , DIMENSION(N) :: Y
   INTEGER , INTENT(INOUT) :: Indic
   INTEGER , INTENT(OUT) :: Ierr
!
! Local variable declarations rewritten by SPAG
!
   COMPLEX(REAL64) , DIMENSION(NDMAX+1) , SAVE :: alp , rd
   REAL(REAL64) , SAVE :: alp0 , fnorm , rm , t
   REAL(REAL64) , EXTERNAL :: ddot
   INTEGER , SAVE :: i , i0 , i1 , j , k , ldg , m1
   LOGICAL , intent(in):: verboz
   EXTERNAL hes , exppro_mgsr
!
! End of declarations rewritten by SPAG
!
!     implicit  real*8 (a-h,o-z)
!-----------------------------------------------------------------------
! this subroutine computes the Arnoldi basis and the corresponding
! coeffcient vector in the approximation
!
!        	w  ::= beta  Vm  ym
!               where ym = exp(- Hm *dt) * e1
!
! to the vector exp(-A dt) w where A is an arbitary matrix and
! w is a given input vector. In case job = 0 the arnoldi basis
! is recomputed. Otherwise the
! code assumes assumes that  u(*) contains an already computed
! arnoldi basis and computes only the y-vector (which is stored in v(*))
!-----------------------------------------------------------------------
! on entry:
!----------
! n	= dimension of matrix
!
! m	= dimension of Krylov subspace (= degree of polynomial
!         approximation to the exponential used. )
!
! dt	= scalar by which to multiply matrix. Can be viewed
!         as a time step. dt must be positive [to be fixed].
!
! eps   = scalar indicating the relative error tolerated for the result.
!         the code will try to compute an answer such that
!         norm2(exactanswer-approximation) / norm2(w) .le. eps
!
! u	= work array of size n*(m+1) to contain the Arnoldi basis
!
! w	= real array of length n = input vector to  which exp(-A) is
!         to be applied.
!
! job	= integer. job indicator. If job .lt.  0 then the Arnoldi
!         basis is recomputed. If job .gt. 0 then it is assumed
!         that the user wants to use a previously computed Krylov
!         subspace but a different dt. Thus the Arnoldi basis and
!         the Hessenberg matrix Hm are not recomputed.
!	  In that case the user should not modify the values of beta
!         and the matrices hh and u(n,*) when recalling phipro.
!         job = -1 : recompute basis and get an initial estimate for
!                    time step dt to be used.
!         job = 0  : recompute basis and do not alter dt.
!         job = 1  : do not recompute arnoldi basis.
!
! z     = real work array of  size (m+1)
! wkc   = complex*16 work array of size (m+1)
!
! hh    = work array of size size at least (m+2)*(m+1)
!
! ih+2	= first dimension of hh as declared in the calling program.
!         ih must be .ge. m.
!
!-----------------------------------------------------------------------
! on return:
!-----------
! w2	= resulting vector w2 = exp(-A *dt) * w
! beta  = real equal to the 2-norm of w. Needed if exppro will
!         be recalled with the same Krylov subspace and a different dt.
! errst = rough estimates of the 2-norm of the error.
! hh	= work array of dimension at least (m+2) x (m+1)
!
!-----------------------------------------------------------------------

!------use degree 14 chebyshev all the time --------------------------
   IF ( Indic>=3 ) THEN
      DO k = 1 , N
         U(k,i1) = Y(k)
      ENDDO
      i0 = 1
!
! switch  for Lanczos version
!     i0 = max0(1, i-1)
      CALL exppro_mgsr(N,i0,i1,U,Hh(1,i))
      fnorm = fnorm + ddot(i1,Hh(1,i),1,Hh(1,i),1)
      IF ( Hh(i1,i)==0.0 ) M = i
      IF ( i>=M ) THEN
!--------------done with arnoldi loop ---------------------------------
         rm = dble(M)
         fnorm = dsqrt(fnorm/rm)
!-------get  : beta*e1 into z
         m1 = M + 1
         DO i = 1 , m1
            Hh(i,m1) = 0.0
         ENDDO
!
!     compute initial dt when  job .lt. 1
!
         IF ( Job<0 ) THEN
!
!     t = eps / beta
!
            t = Eps
            DO k = 1 , M - 1
               t = t*(1.0D0-dble(M-k)/rm)
            ENDDO
!
            t = 2.0D0*rm*(t**(1.0D0/rm))/fnorm
            IF ( verboz ) PRINT * , ' t, dt = ' , t , Dt
            t = dmin1(dabs(Dt),t)
            Dt = dsign(t,Dt)
         ENDIF
         CALL spag_block_1
         RETURN
      ENDIF
   ELSE
!
!------input fraction expansion of rational function ------------------
! chebyshev (14,14)
      ldg = 7
      alp0 = 0.183216998528140087E-11
      alp(1) = (0.557503973136501826E+02,-0.204295038779771857E+03)
      rd(1) = (-0.562314417475317895E+01,0.119406921611247440E+01)
      alp(2) = (-0.938666838877006739E+02,0.912874896775456363E+02)
      rd(2) = (-0.508934679728216110E+01,0.358882439228376881E+01)
      alp(3) = (0.469965415550370835E+02,-0.116167609985818103E+02)
      rd(3) = (-0.399337136365302569E+01,0.600483209099604664E+01)
      alp(4) = (-0.961424200626061065E+01,-0.264195613880262669E+01)
      rd(4) = (-0.226978543095856366E+01,0.846173881758693369E+01)
      alp(5) = (0.752722063978321642E+00,0.670367365566377770E+00)
      rd(5) = (0.208756929753827868E+00,0.109912615662209418E+02)
      alp(6) = (-0.188781253158648576E-01,-0.343696176445802414E-01)
      rd(6) = (0.370327340957595652E+01,0.136563731924991884E+02)
      alp(7) = (0.143086431411801849E-03,0.287221133228814096E-03)
      rd(7) = (0.889777151877331107E+01,0.166309842834712071E+02)
!-----------------------------------------------------------------------
!
!     if job .gt. 0 skip arnoldi process:
!
      IF ( Job>0 ) THEN
         CALL spag_block_1
         RETURN
      ENDIF
!------normalize vector w and put in first column of u --
      Beta = dsqrt(ddot(N,W,1,W,1))
!-----------------------------------------------------------------------
      IF ( verboz ) PRINT * , ' In EXPHES, beta ' , Beta
      IF ( Beta==0.0D0 ) THEN
         Ierr = -1
         Indic = 1
         RETURN
      ENDIF
!
      t = 1.0D0/Beta
      DO j = 1 , N
         U(j,1) = W(j)*t
      ENDDO
!------------------Arnoldi loop -------------------------------------
!      fnorm = 0.0d0
      i1 = 1
   ENDIF
   i = i1
   i1 = i + 1
   DO k = 1 , N
      X(k) = U(k,i)
   ENDDO
   Indic = 3
   RETURN
CONTAINS
   SUBROUTINE spag_block_1
!
      Z(1) = Beta
      DO k = 2 , m1
         Z(k) = 0.0D0
      ENDDO
!-------get  : exp(H) * beta*e1
      CALL hes0(ldg,m1,Hh,Ih,Dt,Z,rd,alp,alp0,Wkc)
!-------error estimate
      Errst = dabs(Z(m1))
      IF ( verboz ) PRINT * , ' error estimate =' , Errst
!-----------------------------------------------------------------------
      Indic = 2
   END SUBROUTINE spag_block_1
END SUBROUTINE exphes
!*==mgsr.f90 processed by SPAG 8.04RA 12:07 23 Feb 2024
!!SPAG Open source Personal, Educational or Academic User Clemson University  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
!-----------------------------------------------------------------------
SUBROUTINE exppro_mgsr(N,I0,I1,Ss,R)
   USE ISO_FORTRAN_ENV                 
   IMPLICIT NONE
!
! Dummy argument declarations rewritten by SPAG
!
   INTEGER :: N
   INTEGER , INTENT(IN) :: I1
   INTEGER , INTENT(IN) :: I0
   REAL(REAL64) , INTENT(INOUT) , DIMENSION(N,I1) :: Ss
   REAL(REAL64) , INTENT(INOUT) , DIMENSION(I1) :: R
!
! Local variable declarations rewritten by SPAG
!
   REAL(REAL64) , EXTERNAL :: ddot
   REAL(REAL64) :: hinorm , t
   INTEGER :: i , it , j , k
   REAL(REAL64) , SAVE :: tet
   EXTERNAL daxpy
!
! End of declarations rewritten by SPAG
!
!     implicit  real*8 (a-h,o-z)
!-----------------------------------------------------------------------
! modified gram - schmidt  with  partial  reortho. the vector ss(*,i1) is
! orthogonalized against the first i vectors  of ss  (which  are  already
! orthogonal).  the coefficients of the orthogonalization are returned in
! the array r
!------------------------------------------------------------------------
! local variables
!
   DATA tet/10.0D0/
 
   DO j = 1 , I1
      R(j) = 0.0D0
   ENDDO
   i = I1 - 1
   it = 0
   SPAG_Loop_1_1: DO
      hinorm = 0.0D0
      it = it + 1
      IF ( i/=0 ) THEN
!
         DO j = I0 , i
            t = ddot(N,Ss(1,j),1,Ss(1,I1),1)
            hinorm = hinorm + t**2
            R(j) = R(j) + t
            CALL daxpy(N,-t,Ss(1,j),1,Ss(1,I1),1)
         ENDDO
         t = ddot(N,Ss(1,I1),1,Ss(1,I1),1)
      ENDIF
!
!     test for reorthogonalization see daniel et. al.
!     two reorthogonalization allowed ---
!
      IF ( t*tet>hinorm .OR. it>=2 ) THEN
         t = dsqrt(t)
         R(I1) = t
         IF ( t==0.0D0 ) RETURN
         t = 1.0D0/t
         DO k = 1 , N
            Ss(k,I1) = Ss(k,I1)*t
         ENDDO
         EXIT SPAG_Loop_1_1
      ENDIF
   ENDDO SPAG_Loop_1_1
END SUBROUTINE exppro_mgsr
!*==project.f90 processed by SPAG 8.04RA 12:07 23 Feb 2024
!!SPAG Open source Personal, Educational or Academic User Clemson University  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
!----------end-of-mgsr--------------------------------------------------
!-----------------------------------------------------------------------
SUBROUTINE exppro_project(N,M,U,V,W)
   USE ISO_FORTRAN_ENV                 
   IMPLICIT NONE
!
! Dummy argument declarations rewritten by SPAG
!
   INTEGER , INTENT(IN) :: N
   INTEGER , INTENT(IN) :: M
   REAL(REAL64) , INTENT(IN) , DIMENSION(N,M) :: U
   REAL(REAL64) , INTENT(IN) , DIMENSION(M) :: V
   REAL(REAL64) , INTENT(INOUT) , DIMENSION(N) :: W
!
! Local variable declarations rewritten by SPAG
!
   INTEGER :: j , k
!
! End of declarations rewritten by SPAG
!
!
!     computes the vector w = u * v
!
! local variables
!
 
   DO k = 1 , N
      W(k) = 0.D0
   ENDDO
 
   DO j = 1 , M
      DO k = 1 , N
         W(k) = W(k) + V(j)*U(k,j)
      ENDDO
   ENDDO
END SUBROUTINE exppro_project

!*==hes.f90 processed by SPAG 8.04RA 12:07 23 Feb 2024
!!SPAG Open source Personal, Educational or Academic User Clemson University  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
!-----------------------------------------------------------------------
SUBROUTINE hes0(Ndg,M1,Hh,Ih,Dt,Y,Root,Coef,Coef0,W2)
   USE ISO_FORTRAN_ENV                 
   IMPLICIT NONE
!
! PARAMETER definitions rewritten by SPAG
!
   INTEGER , PARAMETER :: M1MAX = 61
!
! Dummy argument declarations rewritten by SPAG
!
   INTEGER , INTENT(IN) :: Ndg
   INTEGER , INTENT(IN) :: M1
   INTEGER , INTENT(IN) :: Ih
   REAL(REAL64) , INTENT(IN) , DIMENSION(Ih+2,M1) :: Hh
   REAL(REAL64) , INTENT(IN) :: Dt
   REAL(REAL64) , INTENT(INOUT) , DIMENSION(M1) :: Y
   COMPLEX(REAL64) , INTENT(IN) , DIMENSION(Ndg) :: Root
   COMPLEX(REAL64) , INTENT(IN) , DIMENSION(Ndg) :: Coef
   REAL(REAL64) , INTENT(IN) :: Coef0
   COMPLEX(REAL64) , INTENT(INOUT) , DIMENSION(M1) :: W2
!
! Local variable declarations rewritten by SPAG
!
   COMPLEX(REAL64) , DIMENSION(M1MAX+1,M1MAX) :: hloc
   INTEGER :: i , ii , j
   COMPLEX(REAL64) :: t , zpiv
   REAL(REAL64) , DIMENSION(M1MAX) :: yloc
!
! End of declarations rewritten by SPAG
!
!     implicit  real*8 (a-h,o-z)
!--------------------------------------------------------------------
! computes exp ( H dt) * y    (1)
! where H = Hessenberg matrix (hh)
! y	  = arbitrary vector.
! ----------------------------
! ndg	= number of poles as determined by getrat
! m1    = dimension of hessenberg matrix
! hh	= hessenberg matrix (real)
! ih+2	= first dimension of hh
! dt	= scaling factor used for hh (see (1))
! y	= real vector. on return exp(H dt ) y is computed
!         and overwritten on y.
! root  = poles of the rational approximation to exp as
!         computed by getrat
! coef,
! coef0 = coefficients of partial fraction expansion
!
!  exp(t) ~ coef0 +  sum     Real [   coef(i) / (t - root(i)  ]
!                  i=1,ndg
!
! valid for real t.
! coef0 is real, coef(*) is a complex array.
!
!--------------------------------------------------------------------
! local variables
!
!
!      if (m1 .gt. m1max) print *, ' *** ERROR : In HES, M+1 TOO LARGE'
!
!     loop associated with the poles.
!
   DO j = 1 , M1
      yloc(j) = Y(j)
      Y(j) = Y(j)*Coef0
   ENDDO
!
   DO ii = 1 , Ndg
!
!     copy Hessenberg matrix into temporary
!
      DO j = 1 , M1
         DO i = 1 , j + 1
            hloc(i,j) = dcmplx(Dt*Hh(i,j))
         ENDDO
         hloc(j,j) = hloc(j,j) - Root(ii)
         W2(j) = dcmplx(yloc(j))
      ENDDO
!
! forward solve
!
      DO i = 2 , M1
         zpiv = hloc(i,i-1)/hloc(i-1,i-1)
         DO j = i , M1
            hloc(i,j) = hloc(i,j) - zpiv*hloc(i-1,j)
         ENDDO
         W2(i) = W2(i) - zpiv*W2(i-1)
      ENDDO
!
!     backward solve
!
      DO i = M1 , 1 , -1
         t = W2(i)
         DO j = i + 1 , M1
            t = t - hloc(i,j)*W2(j)
         ENDDO
         W2(i) = t/hloc(i,i)
      ENDDO
!
!     accumulate result in y.
!
      DO i = 1 , M1
         Y(i) = Y(i) + dble(Coef(ii)*W2(i))
      ENDDO
   ENDDO
END SUBROUTINE hes0
