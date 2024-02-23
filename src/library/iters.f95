!*==cg.f90 processed by SPAG 8.04RA 15:47 21 Feb 2024
!!SPAG Open source Personal, Educational or Academic User Clemson University  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
!---------end of readunf ----------------------------------------------
!----------------------------------------------------------------------c
!                          S P A R S K I T                             c
!----------------------------------------------------------------------c
!         Basic Iterative Solvers with Reverse Communication           c
!----------------------------------------------------------------------c
!     This file currently has several basic iterative linear system    c
!     solvers. They are:                                               c
!     CG       -- Conjugate Gradient Method                            c
!     CGNR     -- Conjugate Gradient Method (Normal Residual equation) c
!     BCG      -- Bi-Conjugate Gradient Method                         c
!     DBCG     -- BCG with partial pivoting                            c
!     BCGSTAB  -- BCG stabilized                                       c
!     TFQMR    -- Transpose-Free Quasi-Minimum Residual method         c
!     FOM      -- Full Orthogonalization Method                        c
!     GMRES    -- Generalized Minimum RESidual method                  c
!     FGMRES   -- Flexible version of Generalized Minimum              c
!                 RESidual method                                      c
!     DQGMRES  -- Direct versions of Quasi Generalize Minimum          c
!                 Residual method                                      c
!----------------------------------------------------------------------c
!     They all have the following calling sequence:
!      subroutine solver(n, rhs, sol, ipar, fpar, w)
!      integer n, ipar(16)
!      real*8 rhs(n), sol(n), fpar(16), w(*)
!     Where
!     (1) 'n' is the size of the linear system,
!     (2) 'rhs' is the right-hand side of the linear system,
!     (3) 'sol' is the solution to the linear system,
!     (4) 'ipar' is an integer parameter array for the reverse
!     communication protocol,
!     (5) 'fpar' is an floating-point parameter array storing
!     information to and from the iterative solvers.
!     (6) 'w' is the work space (size is specified in ipar)
!
!     They are preconditioned iterative solvers with reverse
!     communication. The preconditioners can be applied from either
!     from left or right or both (specified by ipar(2), see below).
!
!     Author: Kesheng John Wu (kewu@mail.cs.umn.edu) 1993
!
!     NOTES:
!
!     (1) Work space required by each of the iterative solver
!     routines is as follows:
!       CG      == 5 * n
!       CGNR    == 5 * n
!       BCG     == 7 * n
!       DBCG    == 11 * n
!       BCGSTAB == 8 * n
!       TFQMR   == 11 * n
!       FOM     == (n+3)*(m+2) + (m+1)*m/2 (m = ipar(5), default m=15)
!       GMRES   == (n+3)*(m+2) + (m+1)*m/2 (m = ipar(5), default m=15)
!       FGMRES  == 2*n*(m+1) + (m+1)*m/2 + 3*m + 2 (m = ipar(5),
!                  default m=15)
!       DQGMRES == n + lb * (2*n+4) (lb=ipar(5)+1, default lb = 16)
!
!     (2) ALL iterative solvers require a user-supplied DOT-product
!     routine named DISTDOT. The prototype of DISTDOT is
!
!     real*8 function distdot(n,x,ix,y,iy)
!     integer n, ix, iy
!     real*8 x(1+(n-1)*ix), y(1+(n-1)*iy)
!
!     This interface of DISTDOT is exactly the same as that of
!     DDOT (or SDOT if real == real*8) from BLAS-1. It should have
!     same functionality as DDOT on a single processor machine. On a
!     parallel/distributed environment, each processor can perform
!     DDOT on the data it has, then perform a summation on all the
!     partial results.
!
!     (3) To use this set of routines under SPMD/MIMD program paradigm,
!     several things are to be noted: (a) 'n' should be the number of
!     vector elements of 'rhs' that is present on the local processor.
!     (b) if RHS(i) is on processor j, it is expected that SOL(i)
!     will be on the same processor, i.e. the vectors are distributed
!     to each processor in the same way. (c) the preconditioning and
!     stopping criteria specifications have to be the same on all
!     processor involved, ipar and fpar have to be the same on each
!     processor. (d) DISTDOT should be replaced by a distributed
!     dot-product function.
!
!     ..................................................................
!     Reverse Communication Protocols
!
!     When a reverse-communication routine returns, it could be either
!     that the routine has terminated or it simply requires the caller
!     to perform one matrix-vector multiplication. The possible matrices
!     that involve in the matrix-vector multiplications are:
!     A       (the matrix of the linear system),
!     A^T     (A transposed),
!     Ml^{-1} (inverse of the left preconditioner),
!     Ml^{-T} (inverse of the left preconditioner transposed),
!     Mr^{-1} (inverse of the right preconditioner),
!     Mr^{-T} (inverse of the right preconditioner transposed).
!     For all the matrix vector multiplication, v = A u. The input and
!     output vectors are supposed to be part of the work space 'w', and
!     the starting positions of them are stored in ipar(8:9), see below.
!
!     The array 'ipar' is used to store the information about the solver.
!     Here is the list of what each element represents:
!
!     ipar(1) -- status of the call/return.
!     A call to the solver with ipar(1) == 0 will initialize the
!     iterative solver. On return from the iterative solver, ipar(1)
!     carries the status flag which indicates the condition of the
!     return. The status information is divided into two categories,
!     (1) a positive value indicates the solver requires a matrix-vector
!     multiplication,
!     (2) a non-positive value indicates termination of the solver.
!     Here is the current definition:
!       1 == request a matvec with A,
!       2 == request a matvec with A^T,
!       3 == request a left preconditioner solve (Ml^{-1}),
!       4 == request a left preconditioner transposed solve (Ml^{-T}),
!       5 == request a right preconditioner solve (Mr^{-1}),
!       6 == request a right preconditioner transposed solve (Mr^{-T}),
!      10 == request the caller to perform stopping test,
!       0 == normal termination of the solver, satisfied the stopping
!            criteria,
!      -1 == termination because iteration number is greater than the
!            preset limit,
!      -2 == return due to insufficient work space,
!      -3 == return due to anticipated break-down / divide by zero,
!            in the case where Arnoldi procedure is used, additional
!            error code can be found in ipar(12), where ipar(12) is
!            the error code of orthogonalization procedure MGSRO:
!               -1: zero input vector
!               -2: input vector contains abnormal numbers
!               -3: input vector is a linear combination of others
!               -4: trianguler system in GMRES/FOM/etc. has nul rank
!      -4 == the values of fpar(1) and fpar(2) are both <= 0, the valid
!            ranges are 0 <= fpar(1) < 1, 0 <= fpar(2), and they can
!            not be zero at the same time
!      -9 == while trying to detect a break-down, an abnormal number is
!            detected.
!     -10 == return due to some non-numerical reasons, e.g. invalid
!            floating-point numbers etc.
!
!     ipar(2) -- status of the preconditioning:
!       0 == no preconditioning
!       1 == left preconditioning only
!       2 == right preconditioning only
!       3 == both left and right preconditioning
!
!     ipar(3) -- stopping criteria (details of this will be
!     discussed later).
!
!     ipar(4) -- number of elements in the array 'w'. if this is less
!     than the desired size, it will be over-written with the minimum
!     requirement. In which case the status flag ipar(1) = -2.
!
!     ipar(5) -- size of the Krylov subspace (used by GMRES and its
!     variants), e.g. GMRES(ipar(5)), FGMRES(ipar(5)),
!     DQGMRES(ipar(5)).
!
!     ipar(6) -- maximum number of matrix-vector multiplies, if not a
!     positive number the iterative solver will run till convergence
!     test is satisfied.
!
!     ipar(7) -- current number of matrix-vector multiplies. It is
!     incremented after each matrix-vector multiplication. If there
!     is preconditioning, the counter is incremented after the
!     preconditioning associated with each matrix-vector multiplication.
!
!     ipar(8) -- pointer to the input vector to the requested matrix-
!     vector multiplication.
!
!     ipar(9) -- pointer to the output vector of the requested matrix-
!     vector multiplication.
!
!     To perform v = A * u, it is assumed that u is w(ipar(8):ipar(8)+n-1)
!     and v is stored as w(ipar(9):ipar(9)+n-1).
!
!     ipar(10) -- the return address (used to determine where to go to
!     inside the iterative solvers after the caller has performed the
!     requested services).
!
!     ipar(11) -- the result of the external convergence test
!     On final return from the iterative solvers, this value
!     will be reflected by ipar(1) = 0 (details discussed later)
!
!     ipar(12) -- error code of MGSRO, it is
!                  1 if the input vector to MGSRO is linear combination
!                    of others,
!                  0 if MGSRO was successful,
!                 -1 if the input vector to MGSRO is zero,
!                 -2 if the input vector contains invalid number.
!
!     ipar(13) -- number of initializations. During each initilization
!                 residual norm is computed directly from M_l(b - A x).
!
!     ipar(14) to ipar(16) are NOT defined, they are NOT USED by
!     any iterative solver at this time.
!
!     Information about the error and tolerance are stored in the array
!     FPAR. So are some internal variables that need to be saved from
!     one iteration to the next one. Since the internal variables are
!     not the same for each routine, we only define the common ones.
!
!     The first two are input parameters:
!     fpar(1) -- the relative tolerance,
!     fpar(2) -- the absolute tolerance (details discussed later),
!
!     When the iterative solver terminates,
!     fpar(3) -- initial residual/error norm,
!     fpar(4) -- target residual/error norm,
!     fpar(5) -- current residual norm (if available),
!     fpar(6) -- current residual/error norm,
!     fpar(7) -- convergence rate,
!
!     fpar(8:10) are used by some of the iterative solvers to save some
!     internal information.
!
!     fpar(11) -- number of floating-point operations. The iterative
!     solvers will add the number of FLOPS they used to this variable,
!     but they do NOT initialize it, nor add the number of FLOPS due to
!     matrix-vector multiplications (since matvec is outside of the
!     iterative solvers). To insure the correct FLOPS count, the
!     caller should set fpar(11) = 0 before invoking the iterative
!     solvers and account for the number of FLOPS from matrix-vector
!     multiplications and preconditioners.
!
!     fpar(12:16) are not used in current implementation.
!
!     Whether the content of fpar(3), fpar(4) and fpar(6) are residual
!     norms or error norms depends on ipar(3). If the requested
!     convergence test is based on the residual norm, they will be
!     residual norms. If the caller want to test convergence based the
!     error norms (estimated by the norm of the modifications applied
!     to the approximate solution), they will be error norms.
!     Convergence rate is defined by (Fortran 77 statement)
!     fpar(7) = log10(fpar(3) / fpar(6)) / (ipar(7)-ipar(13))
!     If fpar(7) = 0.5, it means that approximately every 2 (= 1/0.5)
!     steps the residual/error norm decrease by a factor of 10.
!
!     ..................................................................
!     Stopping criteria,
!
!     An iterative solver may be terminated due to (1) satisfying
!     convergence test; (2) exceeding iteration limit; (3) insufficient
!     work space; (4) break-down. Checking of the work space is
!     only done in the initialization stage, i.e. when it is called with
!     ipar(1) == 0. A complete convergence test is done after each
!     update of the solutions. Other conditions are monitored
!     continuously.
!
!     With regard to the number of iteration, when ipar(6) is positive,
!     the current iteration number will be checked against it. If
!     current iteration number is greater the ipar(6) than the solver
!     will return with status -1. If ipar(6) is not positive, the
!     iteration will continue until convergence test is satisfied.
!
!     Two things may be used in the convergence tests, one is the
!     residual 2-norm, the other one is 2-norm of the change in the
!     approximate solution. The residual and the change in approximate
!     solution are from the preconditioned system (if preconditioning
!     is applied). The DQGMRES and TFQMR use two estimates for the
!     residual norms. The estimates are not accurate, but they are
!     acceptable in most of the cases. Generally speaking, the error
!     of the TFQMR's estimate is less accurate.
!
!     The convergence test type is indicated by ipar(3). There are four
!     type convergence tests: (1) tests based on the residual norm;
!     (2) tests based on change in approximate solution; (3) caller
!     does not care, the solver choose one from above two on its own;
!     (4) caller will perform the test, the solver should simply continue.
!     Here is the complete definition:
!      -2 == || dx(i) || <= rtol * || rhs || + atol
!      -1 == || dx(i) || <= rtol * || dx(1) || + atol
!       0 == solver will choose test 1 (next)
!       1 == || residual || <= rtol * || initial residual || + atol
!       2 == || residual || <= rtol * || rhs || + atol
!     999 == caller will perform the test
!     where dx(i) denote the change in the solution at the ith update.
!     ||.|| denotes 2-norm. rtol = fpar(1) and atol = fpar(2).
!
!     If the caller is to perform the convergence test, the outcome
!     should be stored in ipar(11).
!     ipar(11) = 0 -- failed the convergence test, iterative solver
!     should continue
!     ipar(11) = 1 -- satisfied convergence test, iterative solver
!     should perform the clean up job and stop.
!
!     Upon return with ipar(1) = 10,
!     ipar(8)  points to the starting position of the change in
!              solution Sx, where the actual solution of the step is
!              x_j = x_0 + M_r^{-1} Sx.
!              Exception: ipar(8) < 0, Sx = 0. It is mostly used by
!              GMRES and variants to indicate (1) Sx was not necessary,
!              (2) intermediate result of Sx is not computed.
!     ipar(9)  points to the starting position of a work vector that
!              can be used by the caller.
!
!     NOTE: the caller should allow the iterative solver to perform
!     clean up job after the external convergence test is satisfied,
!     since some of the iterative solvers do not directly
!     update the 'sol' array. A typical clean-up stage includes
!     performing the final update of the approximate solution and
!     computing the convergence information (e.g. values of fpar(3:7)).
!
!     NOTE: fpar(4) and fpar(6) are not set by the accelerators (the
!     routines implemented here) if ipar(3) = 999.
!
!     ..................................................................
!     Usage:
!
!     To start solving a linear system, the user needs to specify
!     first 6 elements of the ipar, and first 2 elements of fpar.
!     The user may optionally set fpar(11) = 0 if one wants to count
!     the number of floating-point operations. (Note: the iterative
!     solvers will only add the floating-point operations inside
!     themselves, the caller will have to add the FLOPS from the
!     matrix-vector multiplication routines and the preconditioning
!     routines in order to account for all the arithmetic operations.)
!
!     Here is an example:
!     ipar(1) = 0	! always 0 to start an iterative solver
!     ipar(2) = 2	! right preconditioning
!     ipar(3) = 1	! use convergence test scheme 1
!     ipar(4) = 10000	! the 'w' has 10,000 elements
!     ipar(5) = 10	! use *GMRES(10) (e.g. FGMRES(10))
!     ipar(6) = 100	! use at most 100 matvec's
!     fpar(1) = 1.0E-6	! relative tolerance 1.0E-6
!     fpar(2) = 1.0E-10 ! absolute tolerance 1.0E-10
!     fpar(11) = 0.0	! clearing the FLOPS counter
!
!     After the above specifications, one can start to call an iterative
!     solver, say BCG. Here is a piece of pseudo-code showing how it can
!     be done,
!
! 10   call bcg(n,rhs,sol,ipar,fpar,w)
!      if (ipar(1).eq.1) then
!         call amux(n,w(ipar(8)),w(ipar(9)),a,ja,ia)
!         goto 10
!      else if (ipar(1).eq.2) then
!         call atmux(n,w(ipar(8)),w(ipar(9)),a,ja,ia)
!         goto 10
!      else if (ipar(1).eq.3) then
!         left preconditioner solver
!         goto 10
!      else if (ipar(1).eq.4) then
!         left preconditioner transposed solve
!         goto 10
!      else if (ipar(1).eq.5) then
!         right preconditioner solve
!         goto 10
!      else if (ipar(1).eq.6) then
!         right preconditioner transposed solve
!         goto 10
!      else if (ipar(1).eq.10) then
!         call my own stopping test routine
!         goto 10
!      else if (ipar(1).gt.0) then
!         ipar(1) is an unspecified code
!      else
!         the iterative solver terminated with code = ipar(1)
!      endif
!
!     This segment of pseudo-code assumes the matrix is in CSR format,
!     AMUX and ATMUX are two routines from the SPARSKIT MATVEC module.
!     They perform matrix-vector multiplications for CSR matrices,
!     where w(ipar(8)) is the first element of the input vectors to the
!     two routines, and w(ipar(9)) is the first element of the output
!     vectors from them. For simplicity, we did not show the name of
!     the routine that performs the preconditioning operations or the
!     convergence tests.
!-----------------------------------------------------------------------
SUBROUTINE cg(N,Rhs,Sol,Ipar,Fpar,W)
   USE ISO_FORTRAN_ENV                 
   IMPLICIT NONE
!
! Dummy argument declarations rewritten by SPAG
!
   INTEGER :: N
   REAL(REAL64) , DIMENSION(N) :: Rhs
   REAL(REAL64) , DIMENSION(N) :: Sol
   INTEGER , INTENT(INOUT) , DIMENSION(16) :: Ipar
   REAL(REAL64) , INTENT(INOUT) , DIMENSION(16) :: Fpar
   REAL(REAL64) , INTENT(INOUT) , DIMENSION(N,*) :: W
!
! Local variable declarations rewritten by SPAG
!
   REAL(REAL64) , SAVE :: alpha
   LOGICAL , EXTERNAL :: brkdn , stopbis
   REAL(REAL64) , EXTERNAL :: distdot
   INTEGER , SAVE :: i
   LOGICAL , SAVE :: lp , rp
   EXTERNAL bisinit , tidycg
!
! End of declarations rewritten by SPAG
!
!-----------------------------------------------------------------------
!     This is a implementation of the Conjugate Gradient (CG) method
!     for solving linear system.
!
!     NOTE: This is not the PCG algorithm. It is a regular CG algorithm.
!     To be consistent with the other solvers, the preconditioners are
!     applied by performing Ml^{-1} A Mr^{-1} P in place of A P in the
!     CG algorithm.  PCG uses the preconditioner differently.
!
!     fpar(7) is used here internally to store <r, r>.
!     w(:,1) -- residual vector
!     w(:,2) -- P, the conjugate direction
!     w(:,3) -- A P, matrix multiply the conjugate direction
!     w(:,4) -- temporary storage for results of preconditioning
!     w(:,5) -- change in the solution (sol) is stored here until
!               termination of this solver
!-----------------------------------------------------------------------
!     external functions used
!
!
!     local variables
!
!
!     check the status of the call
!
   IF ( Ipar(1)<=0 ) Ipar(10) = 0
   IF ( Ipar(10)==1 ) THEN
      Ipar(7) = Ipar(7) + 1
      Ipar(13) = 1
      DO i = 1 , N
         W(i,2) = Rhs(i) - W(i,3)
      ENDDO
      Fpar(11) = Fpar(11) + N
!
!     if left preconditioned
!
      IF ( lp ) THEN
         Ipar(1) = 3
         Ipar(9) = 1
         Ipar(10) = 2
         RETURN
      ENDIF
   ELSEIF ( Ipar(10)==2 ) THEN
   ELSEIF ( Ipar(10)==3 ) THEN
      CALL spag_block_2
      RETURN
   ELSEIF ( Ipar(10)==4 ) THEN
!
      IF ( lp ) THEN
         Ipar(1) = 3
         Ipar(8) = Ipar(9)
         Ipar(9) = N + N + 1
         Ipar(10) = 5
         RETURN
      ENDIF
      CALL spag_block_3
      RETURN
   ELSEIF ( Ipar(10)==5 ) THEN
      CALL spag_block_3
      RETURN
   ELSEIF ( Ipar(10)==6 ) THEN
      CALL spag_block_4
      RETURN
   ELSEIF ( Ipar(10)==7 ) THEN
      CALL spag_block_6
      RETURN
   ELSE
!
!     initialization
!
      CALL bisinit(Ipar,Fpar,5*N,1,lp,rp,W)
      IF ( Ipar(1)<0 ) RETURN
!
!     request for matrix vector multiplication A*x in the initialization
!
      Ipar(1) = 1
      Ipar(8) = N + 1
      Ipar(9) = Ipar(8) + N
      Ipar(10) = 1
      DO i = 1 , N
         W(i,2) = Sol(i)
      ENDDO
      RETURN
   ENDIF
!
   IF ( lp ) THEN
      DO i = 1 , N
         W(i,2) = W(i,1)
      ENDDO
   ELSE
      DO i = 1 , N
         W(i,1) = W(i,2)
      ENDDO
   ENDIF
!
   Fpar(7) = distdot(N,W,1,W,1)
   Fpar(11) = Fpar(11) + 2*N
   Fpar(3) = sqrt(Fpar(7))
   Fpar(5) = Fpar(3)
   IF ( abs(Ipar(3))==2 ) THEN
      Fpar(4) = Fpar(1)*sqrt(distdot(N,Rhs,1,Rhs,1)) + Fpar(2)
      Fpar(11) = Fpar(11) + 2*N
   ELSEIF ( Ipar(3)/=999 ) THEN
      Fpar(4) = Fpar(1)*Fpar(3) + Fpar(2)
   ENDIF
   CALL spag_block_1
CONTAINS
   SUBROUTINE spag_block_1
!
!     before iteration can continue, we need to compute A * p, which
!     includes the preconditioning operations
!
      IF ( rp ) THEN
         Ipar(1) = 5
         Ipar(8) = N + 1
         IF ( lp ) THEN
            Ipar(9) = Ipar(8) + N
         ELSE
            Ipar(9) = 3*N + 1
         ENDIF
         Ipar(10) = 3
         RETURN
      ENDIF
      CALL spag_block_2
   END SUBROUTINE spag_block_1
   SUBROUTINE spag_block_2
!
      Ipar(1) = 1
      IF ( rp ) THEN
         Ipar(8) = Ipar(9)
      ELSE
         Ipar(8) = N + 1
      ENDIF
      IF ( lp ) THEN
         Ipar(9) = 3*N + 1
      ELSE
         Ipar(9) = N + N + 1
      ENDIF
      Ipar(10) = 4
      RETURN
   END SUBROUTINE spag_block_2
   SUBROUTINE spag_block_3
!
!     continuing with the iterations
!
      Ipar(7) = Ipar(7) + 1
      alpha = distdot(N,W(1,2),1,W(1,3),1)
      Fpar(11) = Fpar(11) + 2*N
      IF ( brkdn(alpha,Ipar) ) THEN
         CALL spag_block_5
         RETURN
      ENDIF
      alpha = Fpar(7)/alpha
      DO i = 1 , N
         W(i,5) = W(i,5) + alpha*W(i,2)
         W(i,1) = W(i,1) - alpha*W(i,3)
      ENDDO
      Fpar(11) = Fpar(11) + 4*N
!
!     are we ready to terminate ?
!
      IF ( Ipar(3)==999 ) THEN
         Ipar(1) = 10
         Ipar(8) = 4*N + 1
         Ipar(9) = 3*N + 1
         Ipar(10) = 6
         RETURN
      ENDIF
      CALL spag_block_4
   END SUBROUTINE spag_block_3
   SUBROUTINE spag_block_4
      IF ( Ipar(3)==999 ) THEN
         IF ( Ipar(11)==1 ) THEN
            CALL spag_block_5
            RETURN
         ENDIF
      ELSEIF ( stopbis(N,Ipar,1,Fpar,W,W(1,2),alpha) ) THEN
         CALL spag_block_5
         RETURN
      ENDIF
!
!     continue the iterations
!
      alpha = Fpar(5)*Fpar(5)/Fpar(7)
      Fpar(7) = Fpar(5)*Fpar(5)
      DO i = 1 , N
         W(i,2) = W(i,1) + alpha*W(i,2)
      ENDDO
      Fpar(11) = Fpar(11) + 2*N
      CALL spag_block_1
      RETURN
   END SUBROUTINE spag_block_4
   SUBROUTINE spag_block_5
!
!     clean up -- necessary to accommodate the right-preconditioning
!
      IF ( rp ) THEN
         IF ( Ipar(1)<0 ) Ipar(12) = Ipar(1)
         Ipar(1) = 5
         Ipar(8) = 4*N + 1
         Ipar(9) = Ipar(8) - N
         Ipar(10) = 7
         RETURN
      ENDIF
      CALL spag_block_6
   END SUBROUTINE spag_block_5
   SUBROUTINE spag_block_6
      IF ( rp ) THEN
         CALL tidycg(N,Ipar,Fpar,Sol,W(1,4))
      ELSE
         CALL tidycg(N,Ipar,Fpar,Sol,W(1,5))
      ENDIF
   END SUBROUTINE spag_block_6
!
END SUBROUTINE cg
!*==cgnr.f90 processed by SPAG 8.04RA 15:47 21 Feb 2024
!!SPAG Open source Personal, Educational or Academic User Clemson University  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
!-----end-of-cg
!-----------------------------------------------------------------------
SUBROUTINE cgnr(N,Rhs,Sol,Ipar,Fpar,Wk)
   USE ISO_FORTRAN_ENV                 
   IMPLICIT NONE
!
! Dummy argument declarations rewritten by SPAG
!
   INTEGER :: N
   REAL(REAL64) , DIMENSION(N) :: Rhs
   REAL(REAL64) , DIMENSION(N) :: Sol
   INTEGER , INTENT(INOUT) , DIMENSION(16) :: Ipar
   REAL(REAL64) , INTENT(INOUT) , DIMENSION(16) :: Fpar
   REAL(REAL64) , INTENT(INOUT) , DIMENSION(N,*) :: Wk
!
! Local variable declarations rewritten by SPAG
!
   REAL(REAL64) , SAVE :: alpha , zz , zzm1
   LOGICAL , EXTERNAL :: brkdn , stopbis
   REAL(REAL64) , EXTERNAL :: distdot
   INTEGER , SAVE :: i
   LOGICAL , SAVE :: lp , rp
   EXTERNAL bisinit , tidycg
!
! End of declarations rewritten by SPAG
!
!-----------------------------------------------------------------------
!     CGNR -- Using CG algorithm solving A x = b by solving
!     Normal Residual equation: A^T A x = A^T b
!     As long as the matrix is not singular, A^T A is symmetric
!     positive definite, therefore CG (CGNR) will converge.
!
!     Usage of the work space:
!     wk(:,1) == residual vector R
!     wk(:,2) == the conjugate direction vector P
!     wk(:,3) == a scratch vector holds A P, or A^T R
!     wk(:,4) == a scratch vector holds intermediate results of the
!                preconditioning
!     wk(:,5) == a place to hold the modification to SOL
!
!     size of the work space WK is required = 5*n
!-----------------------------------------------------------------------
!     external functions used
!
!
!     local variables
!
!
!     check the status of the call
!
   IF ( Ipar(1)<=0 ) Ipar(10) = 0
   IF ( Ipar(10)==1 ) THEN
      Ipar(7) = Ipar(7) + 1
      Ipar(13) = Ipar(13) + 1
      DO i = 1 , N
         Wk(i,1) = Rhs(i) - Wk(i,2)
      ENDDO
      Fpar(11) = Fpar(11) + N
!
!     if left preconditioned, precondition the initial residual
!
      IF ( lp ) THEN
         Ipar(1) = 3
         Ipar(10) = 2
         RETURN
      ENDIF
   ELSEIF ( Ipar(10)==2 ) THEN
   ELSEIF ( Ipar(10)==3 ) THEN
      CALL spag_block_2
      RETURN
   ELSEIF ( Ipar(10)==4 ) THEN
!
      IF ( rp ) THEN
         Ipar(1) = 6
         Ipar(8) = Ipar(9)
         Ipar(9) = N + N + 1
         Ipar(10) = 5
         RETURN
      ENDIF
      CALL spag_block_3
      RETURN
   ELSEIF ( Ipar(10)==5 ) THEN
      CALL spag_block_3
      RETURN
   ELSEIF ( Ipar(10)==6 ) THEN
      CALL spag_block_4
      RETURN
   ELSEIF ( Ipar(10)==7 ) THEN
!
      IF ( lp ) THEN
         Ipar(1) = 3
         Ipar(8) = Ipar(9)
         Ipar(9) = N + N + 1
         Ipar(10) = 8
         RETURN
      ENDIF
      CALL spag_block_5
      RETURN
   ELSEIF ( Ipar(10)==8 ) THEN
      CALL spag_block_5
      RETURN
   ELSEIF ( Ipar(10)==9 ) THEN
      CALL spag_block_6
      RETURN
   ELSEIF ( Ipar(10)==10 ) THEN
      CALL spag_block_8
      RETURN
   ELSE
!
!     initialization
!
      CALL bisinit(Ipar,Fpar,5*N,1,lp,rp,Wk)
      IF ( Ipar(1)<0 ) RETURN
!
!     request for matrix vector multiplication A*x in the initialization
!
      Ipar(1) = 1
      Ipar(8) = 1
      Ipar(9) = 1 + N
      Ipar(10) = 1
      DO i = 1 , N
         Wk(i,1) = Sol(i)
      ENDDO
      RETURN
   ENDIF
!
   IF ( lp ) THEN
      DO i = 1 , N
         Wk(i,1) = Wk(i,2)
      ENDDO
   ENDIF
!
   zz = distdot(N,Wk,1,Wk,1)
   Fpar(11) = Fpar(11) + 2*N
   Fpar(3) = sqrt(zz)
   Fpar(5) = Fpar(3)
   IF ( abs(Ipar(3))==2 ) THEN
      Fpar(4) = Fpar(1)*sqrt(distdot(N,Rhs,1,Rhs,1)) + Fpar(2)
      Fpar(11) = Fpar(11) + 2*N
   ELSEIF ( Ipar(3)/=999 ) THEN
      Fpar(4) = Fpar(1)*Fpar(3) + Fpar(2)
   ENDIF
   CALL spag_block_1
CONTAINS
   SUBROUTINE spag_block_1
!
!     normal iteration begins here, first half of the iteration
!     computes the conjugate direction
!
!
!     request the caller to perform a A^T r --> wk(:,3)
!
      IF ( lp ) THEN
         Ipar(1) = 4
         Ipar(8) = 1
         IF ( rp ) THEN
            Ipar(9) = N + N + 1
         ELSE
            Ipar(9) = 3*N + 1
         ENDIF
         Ipar(10) = 3
         RETURN
      ENDIF
      CALL spag_block_2
   END SUBROUTINE spag_block_1
   SUBROUTINE spag_block_2
!
      Ipar(1) = 2
      IF ( lp ) THEN
         Ipar(8) = Ipar(9)
      ELSE
         Ipar(8) = 1
      ENDIF
      IF ( rp ) THEN
         Ipar(9) = 3*N + 1
      ELSE
         Ipar(9) = N + N + 1
      ENDIF
      Ipar(10) = 4
      RETURN
   END SUBROUTINE spag_block_2
   SUBROUTINE spag_block_3
!
      Ipar(7) = Ipar(7) + 1
      zzm1 = zz
      zz = distdot(N,Wk(1,3),1,Wk(1,3),1)
      Fpar(11) = Fpar(11) + 2*N
      IF ( brkdn(zz,Ipar) ) THEN
         CALL spag_block_7
         RETURN
      ENDIF
      IF ( Ipar(7)>3 ) THEN
         alpha = zz/zzm1
         DO i = 1 , N
            Wk(i,2) = Wk(i,3) + alpha*Wk(i,2)
         ENDDO
         Fpar(11) = Fpar(11) + 2*N
      ELSE
         DO i = 1 , N
            Wk(i,2) = Wk(i,3)
         ENDDO
      ENDIF
!
!     before iteration can continue, we need to compute A * p
!
      IF ( rp ) THEN
         Ipar(1) = 5
         Ipar(8) = N + 1
         IF ( lp ) THEN
            Ipar(9) = Ipar(8) + N
         ELSE
            Ipar(9) = 3*N + 1
         ENDIF
         Ipar(10) = 6
         RETURN
      ENDIF
      CALL spag_block_4
   END SUBROUTINE spag_block_3
   SUBROUTINE spag_block_4
!
      Ipar(1) = 1
      IF ( rp ) THEN
         Ipar(8) = Ipar(9)
      ELSE
         Ipar(8) = N + 1
      ENDIF
      IF ( lp ) THEN
         Ipar(9) = 3*N + 1
      ELSE
         Ipar(9) = N + N + 1
      ENDIF
      Ipar(10) = 7
      RETURN
   END SUBROUTINE spag_block_4
   SUBROUTINE spag_block_5
!
!     update the solution -- accumulate the changes in w(:,5)
!
      Ipar(7) = Ipar(7) + 1
      alpha = distdot(N,Wk(1,3),1,Wk(1,3),1)
      Fpar(11) = Fpar(11) + 2*N
      IF ( brkdn(alpha,Ipar) ) THEN
         CALL spag_block_7
         RETURN
      ENDIF
      alpha = zz/alpha
      DO i = 1 , N
         Wk(i,5) = Wk(i,5) + alpha*Wk(i,2)
         Wk(i,1) = Wk(i,1) - alpha*Wk(i,3)
      ENDDO
      Fpar(11) = Fpar(11) + 4*N
!
!     are we ready to terminate ?
!
      IF ( Ipar(3)==999 ) THEN
         Ipar(1) = 10
         Ipar(8) = 4*N + 1
         Ipar(9) = 3*N + 1
         Ipar(10) = 9
         RETURN
      ENDIF
      CALL spag_block_6
   END SUBROUTINE spag_block_5
   SUBROUTINE spag_block_6
      IF ( Ipar(3)==999 ) THEN
         IF ( Ipar(11)==1 ) THEN
            CALL spag_block_7
            RETURN
         ENDIF
      ELSEIF ( stopbis(N,Ipar,1,Fpar,Wk,Wk(1,2),alpha) ) THEN
         CALL spag_block_7
         RETURN
!
!     continue the iterations
!
      ENDIF
      CALL spag_block_1
      RETURN
   END SUBROUTINE spag_block_6
   SUBROUTINE spag_block_7
!
!     clean up -- necessary to accommodate the right-preconditioning
!
      IF ( rp ) THEN
         IF ( Ipar(1)<0 ) Ipar(12) = Ipar(1)
         Ipar(1) = 5
         Ipar(8) = 4*N + 1
         Ipar(9) = Ipar(8) - N
         Ipar(10) = 10
         RETURN
      ENDIF
      CALL spag_block_8
   END SUBROUTINE spag_block_7
   SUBROUTINE spag_block_8
      IF ( rp ) THEN
         CALL tidycg(N,Ipar,Fpar,Sol,Wk(1,4))
      ELSE
         CALL tidycg(N,Ipar,Fpar,Sol,Wk(1,5))
      ENDIF
   END SUBROUTINE spag_block_8
END SUBROUTINE cgnr
!*==bcg.f90 processed by SPAG 8.04RA 15:47 21 Feb 2024
!!SPAG Open source Personal, Educational or Academic User Clemson University  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
!-----end-of-cgnr
!-----------------------------------------------------------------------
SUBROUTINE bcg(N,Rhs,Sol,Ipar,Fpar,W)
   USE ISO_FORTRAN_ENV                 
   IMPLICIT NONE
!
! PARAMETER definitions rewritten by SPAG
!
   REAL(REAL64) , PARAMETER :: ONE = 1.0D0
!
! Dummy argument declarations rewritten by SPAG
!
   INTEGER :: N
   REAL(REAL64) , DIMENSION(N) :: Rhs
   REAL(REAL64) , DIMENSION(N) :: Sol
   INTEGER , INTENT(INOUT) , DIMENSION(16) :: Ipar
   REAL(REAL64) , INTENT(INOUT) , DIMENSION(16) :: Fpar
   REAL(REAL64) , INTENT(INOUT) , DIMENSION(N,*) :: W
!
! Local variable declarations rewritten by SPAG
!
   REAL(REAL64) , SAVE :: alpha
   LOGICAL , EXTERNAL :: brkdn , stopbis
   REAL(REAL64) , EXTERNAL :: distdot
   INTEGER , SAVE :: i
   LOGICAL , SAVE :: lp , rp
   EXTERNAL bisinit , tidycg
!
! End of declarations rewritten by SPAG
!
!-----------------------------------------------------------------------
!     BCG: Bi Conjugate Gradient method. Programmed with reverse
!     communication, see the header for detailed specifications
!     of the protocol.
!
!     in this routine, before successful return, the fpar's are
!     fpar(3) == initial residual norm
!     fpar(4) == target residual norm
!     fpar(5) == current residual norm
!     fpar(7) == current rho (rhok = <r, s>)
!     fpar(8) == previous rho (rhokm1)
!
!     w(:,1) -- r, the residual
!     w(:,2) -- s, the dual of the 'r'
!     w(:,3) -- p, the projection direction
!     w(:,4) -- q, the dual of the 'p'
!     w(:,5) -- v, a scratch vector to store A*p, or A*q.
!     w(:,6) -- a scratch vector to store intermediate results
!     w(:,7) -- changes in the solution
!-----------------------------------------------------------------------
!     external routines used
!
!
!
!     local variables
!
!
!     status of the program
!
   IF ( Ipar(1)<=0 ) Ipar(10) = 0
   IF ( Ipar(10)==1 ) THEN
      Ipar(7) = Ipar(7) + 1
      Ipar(13) = Ipar(13) + 1
      DO i = 1 , N
         W(i,1) = Rhs(i) - W(i,5)
      ENDDO
      Fpar(11) = Fpar(11) + N
      IF ( lp ) THEN
         Ipar(1) = 3
         Ipar(8) = 1
         Ipar(9) = N + 1
         Ipar(10) = 2
         RETURN
      ENDIF
   ELSEIF ( Ipar(10)==2 ) THEN
   ELSEIF ( Ipar(10)==3 ) THEN
      CALL spag_block_2
      RETURN
   ELSEIF ( Ipar(10)==4 ) THEN
!
      IF ( lp ) THEN
         Ipar(1) = 3
         Ipar(8) = Ipar(9)
         Ipar(9) = 4*N + 1
         Ipar(10) = 5
         RETURN
      ENDIF
      CALL spag_block_3
      RETURN
   ELSEIF ( Ipar(10)==5 ) THEN
      CALL spag_block_3
      RETURN
   ELSEIF ( Ipar(10)==6 ) THEN
      CALL spag_block_4
      RETURN
   ELSEIF ( Ipar(10)==7 ) THEN
      CALL spag_block_5
      RETURN
   ELSEIF ( Ipar(10)==8 ) THEN
!
      IF ( rp ) THEN
         Ipar(1) = 6
         Ipar(8) = Ipar(9)
         Ipar(9) = 4*N + 1
         Ipar(10) = 9
         RETURN
      ENDIF
      CALL spag_block_6
      RETURN
   ELSEIF ( Ipar(10)==9 ) THEN
      CALL spag_block_6
      RETURN
   ELSEIF ( Ipar(10)==10 ) THEN
      CALL spag_block_8
      RETURN
   ELSE
!
!     initialization, initial residual
!
      CALL bisinit(Ipar,Fpar,7*N,1,lp,rp,W)
      IF ( Ipar(1)<0 ) RETURN
!
!     compute initial residual, request a matvecc
!
      Ipar(1) = 1
      Ipar(8) = 3*N + 1
      Ipar(9) = Ipar(8) + N
      DO i = 1 , N
         W(i,4) = Sol(i)
      ENDDO
      Ipar(10) = 1
      RETURN
   ENDIF
!
   IF ( lp ) THEN
      DO i = 1 , N
         W(i,1) = W(i,2)
         W(i,3) = W(i,2)
         W(i,4) = W(i,2)
      ENDDO
   ELSE
      DO i = 1 , N
         W(i,2) = W(i,1)
         W(i,3) = W(i,1)
         W(i,4) = W(i,1)
      ENDDO
   ENDIF
!
   Fpar(7) = distdot(N,W,1,W,1)
   Fpar(11) = Fpar(11) + 2*N
   Fpar(3) = sqrt(Fpar(7))
   Fpar(5) = Fpar(3)
   Fpar(8) = ONE
   IF ( abs(Ipar(3))==2 ) THEN
      Fpar(4) = Fpar(1)*sqrt(distdot(N,Rhs,1,Rhs,1)) + Fpar(2)
      Fpar(11) = Fpar(11) + 2*N
   ELSEIF ( Ipar(3)/=999 ) THEN
      Fpar(4) = Fpar(1)*Fpar(3) + Fpar(2)
   ENDIF
   IF ( Ipar(3)>=0 .AND. Fpar(5)<=Fpar(4) ) THEN
      Fpar(6) = Fpar(5)
      CALL spag_block_7
      RETURN
   ENDIF
   CALL spag_block_1
CONTAINS
   SUBROUTINE spag_block_1
!
!     end of initialization, begin iteration, v = A p
!
      IF ( rp ) THEN
         Ipar(1) = 5
         Ipar(8) = N + N + 1
         IF ( lp ) THEN
            Ipar(9) = 4*N + 1
         ELSE
            Ipar(9) = 5*N + 1
         ENDIF
         Ipar(10) = 3
         RETURN
      ENDIF
      CALL spag_block_2
   END SUBROUTINE spag_block_1
   SUBROUTINE spag_block_2
!
      Ipar(1) = 1
      IF ( rp ) THEN
         Ipar(8) = Ipar(9)
      ELSE
         Ipar(8) = N + N + 1
      ENDIF
      IF ( lp ) THEN
         Ipar(9) = 5*N + 1
      ELSE
         Ipar(9) = 4*N + 1
      ENDIF
      Ipar(10) = 4
      RETURN
   END SUBROUTINE spag_block_2
   SUBROUTINE spag_block_3
!
      Ipar(7) = Ipar(7) + 1
      alpha = distdot(N,W(1,4),1,W(1,5),1)
      Fpar(11) = Fpar(11) + 2*N
      IF ( brkdn(alpha,Ipar) ) THEN
         CALL spag_block_7
         RETURN
      ENDIF
      alpha = Fpar(7)/alpha
      DO i = 1 , N
         W(i,7) = W(i,7) + alpha*W(i,3)
         W(i,1) = W(i,1) - alpha*W(i,5)
      ENDDO
      Fpar(11) = Fpar(11) + 4*N
      IF ( Ipar(3)==999 ) THEN
         Ipar(1) = 10
         Ipar(8) = 6*N + 1
         Ipar(9) = 5*N + 1
         Ipar(10) = 6
         RETURN
      ENDIF
      CALL spag_block_4
   END SUBROUTINE spag_block_3
   SUBROUTINE spag_block_4
      IF ( Ipar(3)==999 ) THEN
         IF ( Ipar(11)==1 ) THEN
            CALL spag_block_7
            RETURN
         ENDIF
      ELSEIF ( stopbis(N,Ipar,1,Fpar,W,W(1,3),alpha) ) THEN
         CALL spag_block_7
         RETURN
      ENDIF
!
!     A^t * x
!
      IF ( lp ) THEN
         Ipar(1) = 4
         Ipar(8) = 3*N + 1
         IF ( rp ) THEN
            Ipar(9) = 4*N + 1
         ELSE
            Ipar(9) = 5*N + 1
         ENDIF
         Ipar(10) = 7
         RETURN
      ENDIF
      CALL spag_block_5
   END SUBROUTINE spag_block_4
   SUBROUTINE spag_block_5
!
      Ipar(1) = 2
      IF ( lp ) THEN
         Ipar(8) = Ipar(9)
      ELSE
         Ipar(8) = 3*N + 1
      ENDIF
      IF ( rp ) THEN
         Ipar(9) = 5*N + 1
      ELSE
         Ipar(9) = 4*N + 1
      ENDIF
      Ipar(10) = 8
      RETURN
   END SUBROUTINE spag_block_5
   SUBROUTINE spag_block_6
!
      Ipar(7) = Ipar(7) + 1
      DO i = 1 , N
         W(i,2) = W(i,2) - alpha*W(i,5)
      ENDDO
      Fpar(8) = Fpar(7)
      Fpar(7) = distdot(N,W,1,W(1,2),1)
      Fpar(11) = Fpar(11) + 4*N
      IF ( brkdn(Fpar(7),Ipar) ) RETURN
      alpha = Fpar(7)/Fpar(8)
      DO i = 1 , N
         W(i,3) = W(i,1) + alpha*W(i,3)
         W(i,4) = W(i,2) + alpha*W(i,4)
      ENDDO
!
!     end of the iterations
!
      Fpar(11) = Fpar(11) + 4*N
      CALL spag_block_1
      RETURN
   END SUBROUTINE spag_block_6
   SUBROUTINE spag_block_7
!
!     some clean up job to do
!
      IF ( rp ) THEN
         IF ( Ipar(1)<0 ) Ipar(12) = Ipar(1)
         Ipar(1) = 5
         Ipar(8) = 6*N + 1
         Ipar(9) = Ipar(8) - N
         Ipar(10) = 10
         RETURN
      ENDIF
      CALL spag_block_8
   END SUBROUTINE spag_block_7
   SUBROUTINE spag_block_8
      IF ( rp ) THEN
         CALL tidycg(N,Ipar,Fpar,Sol,W(1,6))
      ELSE
         CALL tidycg(N,Ipar,Fpar,Sol,W(1,7))
      ENDIF
   END SUBROUTINE spag_block_8
!-----end-of-bcg
END SUBROUTINE bcg
!*==bcgstab.f90 processed by SPAG 8.04RA 15:47 21 Feb 2024
!!SPAG Open source Personal, Educational or Academic User Clemson University  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
!-----------------------------------------------------------------------
SUBROUTINE bcgstab(N,Rhs,Sol,Ipar,Fpar,W)
   USE ISO_FORTRAN_ENV                 
   IMPLICIT NONE
!
! PARAMETER definitions rewritten by SPAG
!
   REAL(REAL64) , PARAMETER :: ONE = 1.0D0
!
! Dummy argument declarations rewritten by SPAG
!
   INTEGER :: N
   REAL(REAL64) , DIMENSION(N) :: Rhs
   REAL(REAL64) , DIMENSION(N) :: Sol
   INTEGER , INTENT(INOUT) , DIMENSION(16) :: Ipar
   REAL(REAL64) , INTENT(INOUT) , DIMENSION(16) :: Fpar
   REAL(REAL64) , INTENT(INOUT) , DIMENSION(N,8) :: W
!
! Local variable declarations rewritten by SPAG
!
   REAL(REAL64) :: alpha , beta , omega , rho
   LOGICAL , EXTERNAL :: brkdn , stopbis
   REAL(REAL64) , EXTERNAL :: distdot
   INTEGER :: i
   LOGICAL , SAVE :: lp , rp
   EXTERNAL bisinit , tidycg
!
! End of declarations rewritten by SPAG
!
!-----------------------------------------------------------------------
!     BCGSTAB --- Bi Conjugate Gradient stabilized (BCGSTAB)
!     This is an improved BCG routine. (1) no matrix transpose is
!     involved. (2) the convergence is smoother.
!
!
!     Algorithm:
!     Initialization - r = b - A x, r0 = r, p = r, rho = (r0, r),
!     Iterate -
!     (1) v = A p
!     (2) alpha = rho / (r0, v)
!     (3) s = r - alpha v
!     (4) t = A s
!     (5) omega = (t, s) / (t, t)
!     (6) x = x + alpha * p + omega * s
!     (7) r = s - omega * t
!     convergence test goes here
!     (8) beta = rho, rho = (r0, r), beta = rho * alpha / (beta * omega)
!         p = r + beta * (p - omega * v)
!
!     in this routine, before successful return, the fpar's are
!     fpar(3) == initial (preconditionied-)residual norm
!     fpar(4) == target (preconditionied-)residual norm
!     fpar(5) == current (preconditionied-)residual norm
!     fpar(6) == current residual norm or error
!     fpar(7) == current rho (rhok = <r, r0>)
!     fpar(8) == alpha
!     fpar(9) == omega
!
!     Usage of the work space W
!     w(:, 1) = r0, the initial residual vector
!     w(:, 2) = r, current residual vector
!     w(:, 3) = s
!     w(:, 4) = t
!     w(:, 5) = v
!     w(:, 6) = p
!     w(:, 7) = tmp, used in preconditioning, etc.
!     w(:, 8) = delta x, the correction to the answer is accumulated
!               here, so that the right-preconditioning may be applied
!               at the end
!-----------------------------------------------------------------------
!     external routines used
!
!
!
!     local variables
!
!
!     where to go
!
   IF ( Ipar(1)>0 ) THEN
      IF ( Ipar(10)==1 ) THEN
         Ipar(7) = Ipar(7) + 1
         Ipar(13) = Ipar(13) + 1
         DO i = 1 , N
            W(i,1) = Rhs(i) - W(i,2)
         ENDDO
         Fpar(11) = Fpar(11) + N
         IF ( lp ) THEN
            Ipar(1) = 3
            Ipar(10) = 2
            RETURN
         ENDIF
         CALL spag_block_2
         RETURN
      ELSEIF ( Ipar(10)==2 ) THEN
         CALL spag_block_2
         RETURN
      ELSEIF ( Ipar(10)==3 ) THEN
         CALL spag_block_4
         RETURN
      ELSEIF ( Ipar(10)==4 ) THEN
         IF ( lp ) THEN
            Ipar(1) = 3
            Ipar(8) = Ipar(9)
            Ipar(9) = 4*N + 1
            Ipar(10) = 5
            RETURN
         ENDIF
         CALL spag_block_5
         RETURN
      ELSEIF ( Ipar(10)==5 ) THEN
         CALL spag_block_5
         RETURN
      ELSEIF ( Ipar(10)==6 ) THEN
         CALL spag_block_6
         RETURN
      ELSEIF ( Ipar(10)==7 ) THEN
         IF ( lp ) THEN
            Ipar(1) = 3
            Ipar(8) = Ipar(9)
            Ipar(9) = 3*N + 1
            Ipar(10) = 8
            RETURN
         ENDIF
      ELSEIF ( Ipar(10)==8 ) THEN
      ELSEIF ( Ipar(10)==9 ) THEN
         CALL spag_block_7
         RETURN
      ELSEIF ( Ipar(10)==10 ) THEN
         CALL spag_block_9
         RETURN
      ELSE
         CALL spag_block_1
         RETURN
      ENDIF
      Ipar(7) = Ipar(7) + 1
!
!     step (5)
      omega = distdot(N,W(1,4),1,W(1,4),1)
      Fpar(11) = Fpar(11) + N + N
      IF ( brkdn(omega,Ipar) ) THEN
         CALL spag_block_8
         RETURN
      ENDIF
      omega = distdot(N,W(1,4),1,W(1,3),1)/omega
      Fpar(11) = Fpar(11) + N + N
      IF ( brkdn(omega,Ipar) ) THEN
         CALL spag_block_8
         RETURN
      ENDIF
      Fpar(9) = omega
      alpha = Fpar(8)
!
!     step (6) and (7)
      DO i = 1 , N
         W(i,7) = alpha*W(i,6) + omega*W(i,3)
         W(i,8) = W(i,8) + W(i,7)
         W(i,2) = W(i,3) - omega*W(i,4)
      ENDDO
      Fpar(11) = Fpar(11) + 6*N + 1
!
!     convergence test
      IF ( Ipar(3)==999 ) THEN
         Ipar(1) = 10
         Ipar(8) = 7*N + 1
         Ipar(9) = 6*N + 1
         Ipar(10) = 9
         RETURN
      ENDIF
      IF ( .NOT.(stopbis(N,Ipar,2,Fpar,W(1,2),W(1,7),ONE)) ) THEN
         CALL spag_block_7
         RETURN
      ENDIF
      CALL spag_block_8
      RETURN
   ELSEIF ( Ipar(1)<0 ) THEN
      CALL spag_block_8
      RETURN
   ENDIF
   CALL spag_block_1
CONTAINS
   SUBROUTINE spag_block_1
!
!     call the initialization routine
!
      CALL bisinit(Ipar,Fpar,8*N,1,lp,rp,W)
      IF ( Ipar(1)<0 ) RETURN
!
!     perform a matvec to compute the initial residual
!
      Ipar(1) = 1
      Ipar(8) = 1
      Ipar(9) = 1 + N
      DO i = 1 , N
         W(i,1) = Sol(i)
      ENDDO
      Ipar(10) = 1
      RETURN
   END SUBROUTINE spag_block_1
   SUBROUTINE spag_block_2
!
      IF ( lp ) THEN
         DO i = 1 , N
            W(i,1) = W(i,2)
            W(i,6) = W(i,2)
         ENDDO
      ELSE
         DO i = 1 , N
            W(i,2) = W(i,1)
            W(i,6) = W(i,1)
         ENDDO
      ENDIF
!
      Fpar(7) = distdot(N,W,1,W,1)
      Fpar(11) = Fpar(11) + 2*N
      Fpar(5) = sqrt(Fpar(7))
      Fpar(3) = Fpar(5)
      IF ( abs(Ipar(3))==2 ) THEN
         Fpar(4) = Fpar(1)*sqrt(distdot(N,Rhs,1,Rhs,1)) + Fpar(2)
         Fpar(11) = Fpar(11) + 2*N
      ELSEIF ( Ipar(3)/=999 ) THEN
         Fpar(4) = Fpar(1)*Fpar(3) + Fpar(2)
      ENDIF
      IF ( Ipar(3)>=0 ) Fpar(6) = Fpar(5)
      IF ( Ipar(3)>=0 .AND. Fpar(5)<=Fpar(4) .AND. Ipar(3)/=999 ) THEN
         CALL spag_block_8
         RETURN
      ENDIF
      CALL spag_block_3
   END SUBROUTINE spag_block_2
   SUBROUTINE spag_block_3
!
!     beginning of the iterations
!
!     Step (1), v = A p
      IF ( rp ) THEN
         Ipar(1) = 5
         Ipar(8) = 5*N + 1
         IF ( lp ) THEN
            Ipar(9) = 4*N + 1
         ELSE
            Ipar(9) = 6*N + 1
         ENDIF
         Ipar(10) = 3
         RETURN
      ENDIF
      CALL spag_block_4
   END SUBROUTINE spag_block_3
   SUBROUTINE spag_block_4
!
      Ipar(1) = 1
      IF ( rp ) THEN
         Ipar(8) = Ipar(9)
      ELSE
         Ipar(8) = 5*N + 1
      ENDIF
      IF ( lp ) THEN
         Ipar(9) = 6*N + 1
      ELSE
         Ipar(9) = 4*N + 1
      ENDIF
      Ipar(10) = 4
      RETURN
   END SUBROUTINE spag_block_4
   SUBROUTINE spag_block_5
!
      Ipar(7) = Ipar(7) + 1
!
!     step (2)
      alpha = distdot(N,W(1,1),1,W(1,5),1)
      Fpar(11) = Fpar(11) + 2*N
      IF ( brkdn(alpha,Ipar) ) THEN
         CALL spag_block_8
         RETURN
      ENDIF
      alpha = Fpar(7)/alpha
      Fpar(8) = alpha
!
!     step (3)
      DO i = 1 , N
         W(i,3) = W(i,2) - alpha*W(i,5)
      ENDDO
      Fpar(11) = Fpar(11) + 2*N
!
!     Step (4): the second matvec -- t = A s
!
      IF ( rp ) THEN
         Ipar(1) = 5
         Ipar(8) = N + N + 1
         IF ( lp ) THEN
            Ipar(9) = Ipar(8) + N
         ELSE
            Ipar(9) = 6*N + 1
         ENDIF
         Ipar(10) = 6
         RETURN
      ENDIF
      CALL spag_block_6
   END SUBROUTINE spag_block_5
   SUBROUTINE spag_block_6
!
      Ipar(1) = 1
      IF ( rp ) THEN
         Ipar(8) = Ipar(9)
      ELSE
         Ipar(8) = N + N + 1
      ENDIF
      IF ( lp ) THEN
         Ipar(9) = 6*N + 1
      ELSE
         Ipar(9) = 3*N + 1
      ENDIF
      Ipar(10) = 7
      RETURN
   END SUBROUTINE spag_block_6
   SUBROUTINE spag_block_7
      IF ( Ipar(3)/=999 .OR. Ipar(11)/=1 ) THEN
!
!     step (8): computing new p and rho
         rho = Fpar(7)
         Fpar(7) = distdot(N,W(1,2),1,W(1,1),1)
         omega = Fpar(9)
         beta = Fpar(7)*Fpar(8)/(Fpar(9)*rho)
         DO i = 1 , N
            W(i,6) = W(i,2) + beta*(W(i,6)-omega*W(i,5))
         ENDDO
         Fpar(11) = Fpar(11) + 6*N + 3
         IF ( .NOT.(brkdn(Fpar(7),Ipar)) ) THEN
            CALL spag_block_3
            RETURN
!
!     end of an iteration
!
         ENDIF
      ENDIF
      CALL spag_block_8
   END SUBROUTINE spag_block_7
   SUBROUTINE spag_block_8
!
!     some clean up job to do
!
      IF ( rp ) THEN
         IF ( Ipar(1)<0 ) Ipar(12) = Ipar(1)
         Ipar(1) = 5
         Ipar(8) = 7*N + 1
         Ipar(9) = Ipar(8) - N
         Ipar(10) = 10
         RETURN
      ENDIF
      CALL spag_block_9
   END SUBROUTINE spag_block_8
   SUBROUTINE spag_block_9
      IF ( rp ) THEN
         CALL tidycg(N,Ipar,Fpar,Sol,W(1,7))
      ELSE
         CALL tidycg(N,Ipar,Fpar,Sol,W(1,8))
      ENDIF
   END SUBROUTINE spag_block_9
!
!-----end-of-bcgstab
END SUBROUTINE bcgstab
!*==tfqmr.f90 processed by SPAG 8.04RA 15:47 21 Feb 2024
!!SPAG Open source Personal, Educational or Academic User Clemson University  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
!-----------------------------------------------------------------------
SUBROUTINE tfqmr(N,Rhs,Sol,Ipar,Fpar,W)
   USE ISO_FORTRAN_ENV                 
   IMPLICIT NONE
!
! PARAMETER definitions rewritten by SPAG
!
   REAL(REAL64) , PARAMETER :: ONE = 1.0D0 , ZERO = 0.0D0
!
! Dummy argument declarations rewritten by SPAG
!
   INTEGER :: N
   REAL(REAL64) , DIMENSION(N) :: Rhs
   REAL(REAL64) , DIMENSION(N) :: Sol
   INTEGER , INTENT(INOUT) , DIMENSION(16) :: Ipar
   REAL(REAL64) , INTENT(INOUT) , DIMENSION(16) :: Fpar
   REAL(REAL64) , INTENT(INOUT) , DIMENSION(N,*) :: W
!
! Local variable declarations rewritten by SPAG
!
   REAL(REAL64) , SAVE :: alpha , eta , rho , sigma , tao , te , theta
   LOGICAL , EXTERNAL :: brkdn
   REAL(REAL64) , EXTERNAL :: distdot
   INTEGER , SAVE :: i
   LOGICAL , SAVE :: lp , rp
   EXTERNAL bisinit , tidycg
!
! End of declarations rewritten by SPAG
!
!-----------------------------------------------------------------------
!     TFQMR --- transpose-free Quasi-Minimum Residual method
!     This is developed from BCG based on the principle of Quasi-Minimum
!     Residual, and it is transpose-free.
!
!     It uses approximate residual norm.
!
!     Internally, the fpar's are used as following:
!     fpar(3) --- initial residual norm squared
!     fpar(4) --- target residual norm squared
!     fpar(5) --- current residual norm squared
!
!     w(:,1) -- R, residual
!     w(:,2) -- R0, the initial residual
!     w(:,3) -- W
!     w(:,4) -- Y
!     w(:,5) -- Z
!     w(:,6) -- A * Y
!     w(:,7) -- A * Z
!     w(:,8) -- V
!     w(:,9) -- D
!     w(:,10) -- intermediate results of preconditioning
!     w(:,11) -- changes in the solution
!-----------------------------------------------------------------------
!     external functions
!
!
!
!     local variables
!
!
!     status of the call (where to go)
!
   IF ( Ipar(1)<=0 ) Ipar(10) = 0
   IF ( Ipar(10)==1 ) THEN
      Ipar(7) = Ipar(7) + 1
      Ipar(13) = Ipar(13) + 1
      DO i = 1 , N
         W(i,1) = Rhs(i) - W(i,7)
         W(i,9) = ZERO
      ENDDO
      Fpar(11) = Fpar(11) + N
!
      IF ( lp ) THEN
         Ipar(1) = 3
         Ipar(9) = N + 1
         Ipar(10) = 2
         RETURN
      ENDIF
   ELSEIF ( Ipar(10)==2 ) THEN
   ELSEIF ( Ipar(10)==3 ) THEN
      CALL spag_block_2
      RETURN
   ELSEIF ( Ipar(10)==4 ) THEN
!
      IF ( lp ) THEN
         Ipar(1) = 3
         Ipar(8) = Ipar(9)
         Ipar(9) = 5*N + 1
         Ipar(10) = 5
         RETURN
      ENDIF
      CALL spag_block_3
      RETURN
   ELSEIF ( Ipar(10)==5 ) THEN
      CALL spag_block_3
      RETURN
   ELSEIF ( Ipar(10)==6 ) THEN
      CALL spag_block_4
      RETURN
   ELSEIF ( Ipar(10)==7 ) THEN
!
      IF ( lp ) THEN
         Ipar(1) = 3
         Ipar(8) = Ipar(9)
         Ipar(9) = 6*N + 1
         Ipar(10) = 8
         RETURN
      ENDIF
      CALL spag_block_5
      RETURN
   ELSEIF ( Ipar(10)==8 ) THEN
      CALL spag_block_5
      RETURN
   ELSEIF ( Ipar(10)==9 ) THEN
      CALL spag_block_6
      RETURN
   ELSEIF ( Ipar(10)==10 ) THEN
      CALL spag_block_8
      RETURN
   ELSE
!
!     initializations
!
      CALL bisinit(Ipar,Fpar,11*N,2,lp,rp,W)
      IF ( Ipar(1)<0 ) RETURN
      Ipar(1) = 1
      Ipar(8) = 1
      Ipar(9) = 1 + 6*N
      DO i = 1 , N
         W(i,1) = Sol(i)
      ENDDO
      Ipar(10) = 1
      RETURN
   ENDIF
   IF ( lp ) THEN
      DO i = 1 , N
         W(i,1) = W(i,2)
         W(i,3) = W(i,2)
      ENDDO
   ELSE
      DO i = 1 , N
         W(i,2) = W(i,1)
         W(i,3) = W(i,1)
      ENDDO
   ENDIF
!
   Fpar(5) = sqrt(distdot(N,W,1,W,1))
   Fpar(3) = Fpar(5)
   tao = Fpar(5)
   Fpar(11) = Fpar(11) + N + N
   IF ( abs(Ipar(3))==2 ) THEN
      Fpar(4) = Fpar(1)*sqrt(distdot(N,Rhs,1,Rhs,1)) + Fpar(2)
      Fpar(11) = Fpar(11) + N + N
   ELSEIF ( Ipar(3)/=999 ) THEN
      Fpar(4) = Fpar(1)*tao + Fpar(2)
   ENDIF
   te = ZERO
   rho = ZERO
   CALL spag_block_1
CONTAINS
   SUBROUTINE spag_block_1
!
!     begin iteration
!
      sigma = rho
      rho = distdot(N,W(1,2),1,W(1,3),1)
      Fpar(11) = Fpar(11) + N + N
      IF ( brkdn(rho,Ipar) ) THEN
         CALL spag_block_7
         RETURN
      ENDIF
      IF ( Ipar(7)==1 ) THEN
         alpha = ZERO
      ELSE
         alpha = rho/sigma
      ENDIF
      DO i = 1 , N
         W(i,4) = W(i,3) + alpha*W(i,5)
      ENDDO
      Fpar(11) = Fpar(11) + N + N
!
!     A * x -- with preconditioning
!
      IF ( rp ) THEN
         Ipar(1) = 5
         Ipar(8) = 3*N + 1
         IF ( lp ) THEN
            Ipar(9) = 5*N + 1
         ELSE
            Ipar(9) = 9*N + 1
         ENDIF
         Ipar(10) = 3
         RETURN
      ENDIF
      CALL spag_block_2
   END SUBROUTINE spag_block_1
   SUBROUTINE spag_block_2
!
      Ipar(1) = 1
      IF ( rp ) THEN
         Ipar(8) = Ipar(9)
      ELSE
         Ipar(8) = 3*N + 1
      ENDIF
      IF ( lp ) THEN
         Ipar(9) = 9*N + 1
      ELSE
         Ipar(9) = 5*N + 1
      ENDIF
      Ipar(10) = 4
      RETURN
   END SUBROUTINE spag_block_2
   SUBROUTINE spag_block_3
      Ipar(7) = Ipar(7) + 1
      DO i = 1 , N
         W(i,8) = W(i,6) + alpha*(W(i,7)+alpha*W(i,8))
      ENDDO
      sigma = distdot(N,W(1,2),1,W(1,8),1)
      Fpar(11) = Fpar(11) + 6*N
      IF ( brkdn(sigma,Ipar) ) THEN
         CALL spag_block_7
         RETURN
      ENDIF
      alpha = rho/sigma
      DO i = 1 , N
         W(i,5) = W(i,4) - alpha*W(i,8)
      ENDDO
      Fpar(11) = Fpar(11) + 2*N
!
!     the second A * x
!
      IF ( rp ) THEN
         Ipar(1) = 5
         Ipar(8) = 4*N + 1
         IF ( lp ) THEN
            Ipar(9) = 6*N + 1
         ELSE
            Ipar(9) = 9*N + 1
         ENDIF
         Ipar(10) = 6
         RETURN
      ENDIF
      CALL spag_block_4
   END SUBROUTINE spag_block_3
   SUBROUTINE spag_block_4
!
      Ipar(1) = 1
      IF ( rp ) THEN
         Ipar(8) = Ipar(9)
      ELSE
         Ipar(8) = 4*N + 1
      ENDIF
      IF ( lp ) THEN
         Ipar(9) = 9*N + 1
      ELSE
         Ipar(9) = 6*N + 1
      ENDIF
      Ipar(10) = 7
      RETURN
   END SUBROUTINE spag_block_4
   SUBROUTINE spag_block_5
      Ipar(7) = Ipar(7) + 1
      DO i = 1 , N
         W(i,3) = W(i,3) - alpha*W(i,6)
      ENDDO
!
!     update I
!
      theta = distdot(N,W(1,3),1,W(1,3),1)/(tao*tao)
      sigma = ONE/(ONE+theta)
      tao = tao*sqrt(sigma*theta)
      Fpar(11) = Fpar(11) + 4*N + 6
      IF ( brkdn(tao,Ipar) ) THEN
         CALL spag_block_7
         RETURN
      ENDIF
      eta = sigma*alpha
      sigma = te/alpha
      te = theta*eta
      DO i = 1 , N
         W(i,9) = W(i,4) + sigma*W(i,9)
         W(i,11) = W(i,11) + eta*W(i,9)
         W(i,3) = W(i,3) - alpha*W(i,7)
      ENDDO
      Fpar(11) = Fpar(11) + 6*N + 6
      IF ( Ipar(7)==1 ) THEN
         IF ( Ipar(3)==-1 ) THEN
            Fpar(3) = eta*sqrt(distdot(N,W(1,9),1,W(1,9),1))
            Fpar(4) = Fpar(1)*Fpar(3) + Fpar(2)
            Fpar(11) = Fpar(11) + N + N + 4
         ENDIF
      ENDIF
!
!     update II
!
      theta = distdot(N,W(1,3),1,W(1,3),1)/(tao*tao)
      sigma = ONE/(ONE+theta)
      tao = tao*sqrt(sigma*theta)
      Fpar(11) = Fpar(11) + 8 + 2*N
      IF ( brkdn(tao,Ipar) ) THEN
         CALL spag_block_7
         RETURN
      ENDIF
      eta = sigma*alpha
      sigma = te/alpha
      te = theta*eta
      DO i = 1 , N
         W(i,9) = W(i,5) + sigma*W(i,9)
         W(i,11) = W(i,11) + eta*W(i,9)
      ENDDO
      Fpar(11) = Fpar(11) + 4*N + 3
!
!     this is the correct over-estimate
!      fpar(5) = sqrt(real(ipar(7)+1)) * tao
!     this is an approximation
      Fpar(5) = tao
      IF ( Ipar(3)==999 ) THEN
         Ipar(1) = 10
         Ipar(8) = 10*N + 1
         Ipar(9) = 9*N + 1
         Ipar(10) = 9
         RETURN
      ELSEIF ( Ipar(3)<0 ) THEN
         Fpar(6) = eta*sqrt(distdot(N,W(1,9),1,W(1,9),1))
         Fpar(11) = Fpar(11) + N + N + 2
      ELSE
         Fpar(6) = Fpar(5)
      ENDIF
      IF ( Fpar(6)>Fpar(4) .AND. (Ipar(7)<Ipar(6) .OR. Ipar(6)<=0) ) THEN
         CALL spag_block_1
         RETURN
      ENDIF
      CALL spag_block_6
   END SUBROUTINE spag_block_5
   SUBROUTINE spag_block_6
      IF ( Ipar(3)==999 .AND. Ipar(11)==0 ) THEN
         CALL spag_block_1
         RETURN
      ENDIF
      CALL spag_block_7
   END SUBROUTINE spag_block_6
   SUBROUTINE spag_block_7
!
!     clean up
!
      IF ( rp ) THEN
         IF ( Ipar(1)<0 ) Ipar(12) = Ipar(1)
         Ipar(1) = 5
         Ipar(8) = 10*N + 1
         Ipar(9) = Ipar(8) - N
         Ipar(10) = 10
         RETURN
      ENDIF
      CALL spag_block_8
   END SUBROUTINE spag_block_7
   SUBROUTINE spag_block_8
      IF ( rp ) THEN
         CALL tidycg(N,Ipar,Fpar,Sol,W(1,10))
      ELSE
         CALL tidycg(N,Ipar,Fpar,Sol,W(1,11))
      ENDIF
   END SUBROUTINE spag_block_8
!
END SUBROUTINE tfqmr
!*==fom.f90 processed by SPAG 8.04RA 15:47 21 Feb 2024
!!SPAG Open source Personal, Educational or Academic User Clemson University  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
!-----end-of-tfqmr
!-----------------------------------------------------------------------
SUBROUTINE fom(N,Rhs,Sol,Ipar,Fpar,W)
   USE ISO_FORTRAN_ENV                 
   IMPLICIT NONE
!
! PARAMETER definitions rewritten by SPAG
!
   REAL(REAL64) , PARAMETER :: ONE = 1.0D0 , ZERO = 0.0D0
!
! Dummy argument declarations rewritten by SPAG
!
   INTEGER :: N
   REAL(REAL64) , DIMENSION(N) :: Rhs
   REAL(REAL64) , INTENT(INOUT) , DIMENSION(N) :: Sol
   INTEGER , INTENT(INOUT) , DIMENSION(16) :: Ipar
   REAL(REAL64) , INTENT(INOUT) , DIMENSION(16) :: Fpar
   REAL(REAL64) , INTENT(INOUT) , DIMENSION(*) :: W
!
! Local variable declarations rewritten by SPAG
!
   REAL(REAL64) , SAVE :: alpha , c , s
   REAL(REAL64) , EXTERNAL :: distdot
   INTEGER , SAVE :: hes , i , idx , ii , k , m , p2 , prs , ptr , vc , vrn , vs
   LOGICAL , SAVE :: lp , rp
   EXTERNAL bisinit , givens , mgsro
!
! End of declarations rewritten by SPAG
!
!-----------------------------------------------------------------------
!     This a version of The Full Orthogonalization Method (FOM)
!     implemented with reverse communication. It is a simple restart
!     version of the FOM algorithm and is implemented with plane
!     rotations similarly to GMRES.
!
!  parameters:
!  -----------
!     ipar(5) == the dimension of the Krylov subspace
!     after every ipar(5) iterations, the FOM will restart with
!     the updated solution and recomputed residual vector.
!
!     the work space in `w' is used as follows:
!     (1) the basis for the Krylov subspace, size n*(m+1);
!     (2) the Hessenberg matrix, only the upper triangular
!     portion of the matrix is stored, size (m+1)*m/2 + 1
!     (3) three vectors, all are of size m, they are
!     the cosine and sine of the Givens rotations, the third one holds
!     the residuals, it is of size m+1.
!
!     TOTAL SIZE REQUIRED == (n+3)*(m+2) + (m+1)*m/2
!     Note: m == ipar(5). The default value for this is 15 if
!     ipar(5) <= 1.
!-----------------------------------------------------------------------
!     external functions used
!
!
!
!     local variables, ptr and p2 are temporary pointers,
!     hes points to the Hessenberg matrix,
!     vc, vs point to the cosines and sines of the Givens rotations
!     vrn points to the vectors of residual norms, more precisely
!     the right hand side of the least square problem solved.
!
   INTEGER :: spag_nextblock_1
   spag_nextblock_1 = 1
   SPAG_DispatchLoop_1: DO
      SELECT CASE (spag_nextblock_1)
      CASE (1)
!
!     check the status of the call
!
         IF ( Ipar(1)<=0 ) Ipar(10) = 0
         IF ( Ipar(10)==1 ) THEN
            Ipar(7) = Ipar(7) + 1
            Ipar(13) = Ipar(13) + 1
            IF ( lp ) THEN
               DO i = 1 , N
                  W(N+i) = Rhs(i) - W(i)
               ENDDO
               Ipar(1) = 3
               Ipar(10) = 2
               RETURN
            ELSE
               DO i = 1 , N
                  W(i) = Rhs(i) - W(i)
               ENDDO
            ENDIF
            Fpar(11) = Fpar(11) + N
            spag_nextblock_1 = 3
            CYCLE SPAG_DispatchLoop_1
         ELSEIF ( Ipar(10)==2 ) THEN
            spag_nextblock_1 = 3
            CYCLE SPAG_DispatchLoop_1
         ELSEIF ( Ipar(10)==3 ) THEN
            spag_nextblock_1 = 5
            CYCLE SPAG_DispatchLoop_1
         ELSEIF ( Ipar(10)==4 ) THEN
!
            IF ( lp ) THEN
               Ipar(1) = 3
               Ipar(8) = Ipar(9)
               Ipar(9) = k*N + 1
               Ipar(10) = 5
               RETURN
            ENDIF
         ELSEIF ( Ipar(10)==5 ) THEN
         ELSEIF ( Ipar(10)==6 ) THEN
            spag_nextblock_1 = 7
            CYCLE SPAG_DispatchLoop_1
         ELSEIF ( Ipar(10)==7 ) THEN
            spag_nextblock_1 = 8
            CYCLE SPAG_DispatchLoop_1
         ELSE
!
!     initialization
!
            IF ( Ipar(5)<=1 ) THEN
               m = 15
            ELSE
               m = Ipar(5)
            ENDIF
            idx = N*(m+1)
            hes = idx + N
            vc = hes + (m+1)*m/2 + 1
            vs = vc + m
            vrn = vs + m
            i = vrn + m + 1
            CALL bisinit(Ipar,Fpar,i,1,lp,rp,W)
            IF ( Ipar(1)<0 ) RETURN
            spag_nextblock_1 = 2
            CYCLE SPAG_DispatchLoop_1
         ENDIF
!
!     Modified Gram-Schmidt orthogonalization procedure
!     temporary pointer 'ptr' is pointing to the current column of the
!     Hessenberg matrix. 'p2' points to the new basis vector
!
         Ipar(7) = Ipar(7) + 1
         ptr = k*(k-1)/2 + hes
         p2 = Ipar(9)
         CALL mgsro(.FALSE.,N,N,k+1,k+1,Fpar(11),W,W(ptr+1),Ipar(12))
         IF ( Ipar(12)<0 ) THEN
            spag_nextblock_1 = 6
            CYCLE SPAG_DispatchLoop_1
         ENDIF
!
!     apply previous Givens rotations to column.
!
         p2 = ptr + 1
         DO i = 1 , k - 1
            ptr = p2
            p2 = p2 + 1
            alpha = W(ptr)
            c = W(vc+i)
            s = W(vs+i)
            W(ptr) = c*alpha + s*W(p2)
            W(p2) = c*W(p2) - s*alpha
         ENDDO
!
!     end of one Arnoldi iteration, alpha will store the estimated
!     residual norm at current stage
!
         Fpar(11) = Fpar(11) + 6*k
 
         prs = vrn + k
         alpha = Fpar(5)
         IF ( W(p2)/=ZERO ) alpha = abs(W(p2+1)*W(prs)/W(p2))
         Fpar(5) = alpha
!
         IF ( (k>=m) .OR. (Ipar(3)>=0 .AND. alpha<=Fpar(4)) .OR. (Ipar(6)>0 .AND. Ipar(7)>=Ipar(6)) ) THEN
            spag_nextblock_1 = 6
            CYCLE SPAG_DispatchLoop_1
         ENDIF
!
         CALL givens(W(p2),W(p2+1),c,s)
         W(vc+k) = c
         W(vs+k) = s
         alpha = -s*W(prs)
         W(prs) = c*W(prs)
         W(prs+1) = alpha
!
         IF ( W(p2)==ZERO ) THEN
            spag_nextblock_1 = 6
            CYCLE SPAG_DispatchLoop_1
         ENDIF
         spag_nextblock_1 = 4
         CYCLE SPAG_DispatchLoop_1
      CASE (2)
!
!     request for matrix vector multiplication A*x in the initialization
!
         Ipar(1) = 1
         Ipar(8) = N + 1
         Ipar(9) = 1
         Ipar(10) = 1
         k = 0
         DO i = 1 , N
            W(N+i) = Sol(i)
         ENDDO
         RETURN
      CASE (3)
!
         alpha = sqrt(distdot(N,W,1,W,1))
         Fpar(11) = Fpar(11) + 2*N + 1
         IF ( Ipar(7)==1 .AND. Ipar(3)/=999 ) THEN
            IF ( abs(Ipar(3))==2 ) THEN
               Fpar(4) = Fpar(1)*sqrt(distdot(N,Rhs,1,Rhs,1)) + Fpar(2)
               Fpar(11) = Fpar(11) + 2*N
            ELSE
               Fpar(4) = Fpar(1)*alpha + Fpar(2)
            ENDIF
            Fpar(3) = alpha
         ENDIF
         Fpar(5) = alpha
         W(vrn+1) = alpha
         IF ( alpha<=Fpar(4) .AND. Ipar(3)>=0 .AND. Ipar(3)/=999 ) THEN
            Ipar(1) = 0
            Fpar(6) = alpha
            spag_nextblock_1 = 9
            CYCLE SPAG_DispatchLoop_1
         ENDIF
         alpha = ONE/alpha
         DO ii = 1 , N
            W(ii) = alpha*W(ii)
         ENDDO
         Fpar(11) = Fpar(11) + N
         spag_nextblock_1 = 4
      CASE (4)
!
!     request for (1) right preconditioning
!     (2) matrix vector multiplication
!     (3) left preconditioning
!
         k = k + 1
         IF ( rp ) THEN
            Ipar(1) = 5
            Ipar(8) = k*N - N + 1
            IF ( lp ) THEN
               Ipar(9) = k*N + 1
            ELSE
               Ipar(9) = idx + 1
            ENDIF
            Ipar(10) = 3
            RETURN
         ENDIF
         spag_nextblock_1 = 5
      CASE (5)
!
         Ipar(1) = 1
         IF ( rp ) THEN
            Ipar(8) = Ipar(9)
         ELSE
            Ipar(8) = (k-1)*N + 1
         ENDIF
         IF ( lp ) THEN
            Ipar(9) = idx + 1
         ELSE
            Ipar(9) = 1 + k*N
         ENDIF
         Ipar(10) = 4
         RETURN
      CASE (6)
         DO
!
!     update the approximate solution, first solve the upper triangular
!     system, temporary pointer ptr points to the Hessenberg matrix,
!     prs points to the right-hand-side (also the solution) of the system.
!
            ptr = hes + k*(k+1)/2
            prs = vrn + k
            IF ( W(ptr)==ZERO ) THEN
!
!     if the diagonal elements of the last column is zero, reduce k by 1
!     so that a smaller trianguler system is solved
!
               k = k - 1
               IF ( k>0 ) CYCLE
               Ipar(1) = -3
               Ipar(12) = -4
               spag_nextblock_1 = 9
               CYCLE SPAG_DispatchLoop_1
            ENDIF
            W(prs) = W(prs)/W(ptr)
            DO i = k - 1 , 1 , -1
               ptr = ptr - i - 1
               DO ii = 1 , i
                  W(vrn+ii) = W(vrn+ii) - W(prs)*W(ptr+ii)
               ENDDO
               prs = prs - 1
               W(prs) = W(prs)/W(ptr)
            ENDDO
!
            DO ii = 1 , N
               W(ii) = W(ii)*W(prs)
            ENDDO
            DO i = 1 , k - 1
               prs = prs + 1
               ptr = i*N
               DO ii = 1 , N
                  W(ii) = W(ii) + W(prs)*W(ptr+ii)
               ENDDO
            ENDDO
            Fpar(11) = Fpar(11) + 2*(k-1)*N + N + k*(k+1)
!
            IF ( rp ) THEN
               Ipar(1) = 5
               Ipar(8) = 1
               Ipar(9) = idx + 1
               Ipar(10) = 6
               RETURN
            ENDIF
            spag_nextblock_1 = 7
            CYCLE SPAG_DispatchLoop_1
         ENDDO
         spag_nextblock_1 = 7
      CASE (7)
!
         IF ( rp ) THEN
            DO i = 1 , N
               Sol(i) = Sol(i) + W(idx+i)
            ENDDO
         ELSE
            DO i = 1 , N
               Sol(i) = Sol(i) + W(i)
            ENDDO
         ENDIF
         Fpar(11) = Fpar(11) + N
!
!     process the complete stopping criteria
!
         IF ( Ipar(3)==999 ) THEN
            Ipar(1) = 10
            Ipar(8) = -1
            Ipar(9) = idx + 1
            Ipar(10) = 7
            RETURN
         ELSEIF ( Ipar(3)<0 ) THEN
            IF ( Ipar(7)<=m+1 ) THEN
               Fpar(3) = abs(W(vrn+1))
               IF ( Ipar(3)==-1 ) Fpar(4) = Fpar(1)*Fpar(3) + Fpar(2)
            ENDIF
            alpha = abs(W(vrn+k))
         ENDIF
         Fpar(6) = alpha
         spag_nextblock_1 = 8
      CASE (8)
!
!     do we need to restart ?
!
         IF ( Ipar(12)/=0 ) THEN
            Ipar(1) = -3
            spag_nextblock_1 = 9
            CYCLE SPAG_DispatchLoop_1
         ENDIF
         IF ( Ipar(7)<Ipar(6) .OR. Ipar(6)<=0 ) THEN
            IF ( Ipar(3)/=999 ) THEN
               IF ( Fpar(6)>Fpar(4) ) THEN
                  spag_nextblock_1 = 2
                  CYCLE SPAG_DispatchLoop_1
               ENDIF
            ELSEIF ( Ipar(11)==0 ) THEN
               spag_nextblock_1 = 2
               CYCLE SPAG_DispatchLoop_1
            ENDIF
         ENDIF
!
!     termination, set error code, compute convergence rate
!
         IF ( Ipar(1)>0 ) THEN
            IF ( Ipar(3)==999 .AND. Ipar(11)==1 ) THEN
               Ipar(1) = 0
            ELSEIF ( Ipar(3)/=999 .AND. Fpar(6)<=Fpar(4) ) THEN
               Ipar(1) = 0
            ELSEIF ( Ipar(7)>=Ipar(6) .AND. Ipar(6)>0 ) THEN
               Ipar(1) = -1
            ELSE
               Ipar(1) = -10
            ENDIF
         ENDIF
         spag_nextblock_1 = 9
      CASE (9)
         IF ( Fpar(3)/=ZERO .AND. Fpar(6)/=ZERO .AND. Ipar(7)>Ipar(13) ) THEN
            Fpar(7) = log10(Fpar(3)/Fpar(6))/dble(Ipar(7)-Ipar(13))
         ELSE
            Fpar(7) = ZERO
         ENDIF
         EXIT SPAG_DispatchLoop_1
      END SELECT
   ENDDO SPAG_DispatchLoop_1
END SUBROUTINE fom
!*==gmres.f90 processed by SPAG 8.04RA 15:47 21 Feb 2024
!!SPAG Open source Personal, Educational or Academic User Clemson University  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
!-----end-of-fom--------------------------------------------------------
!-----------------------------------------------------------------------
SUBROUTINE gmres(N,Rhs,Sol,Ipar,Fpar,W)
   USE ISO_FORTRAN_ENV                 
   IMPLICIT NONE
!
! PARAMETER definitions rewritten by SPAG
!
   REAL(REAL64) , PARAMETER :: ONE = 1.0D0 , ZERO = 0.0D0
!
! Dummy argument declarations rewritten by SPAG
!
   INTEGER :: N
   REAL(REAL64) , DIMENSION(N) :: Rhs
   REAL(REAL64) , INTENT(INOUT) , DIMENSION(N) :: Sol
   INTEGER , INTENT(INOUT) , DIMENSION(16) :: Ipar
   REAL(REAL64) , INTENT(INOUT) , DIMENSION(16) :: Fpar
   REAL(REAL64) , INTENT(INOUT) , DIMENSION(*) :: W
!
! Local variable declarations rewritten by SPAG
!
   REAL(REAL64) , SAVE :: alpha , c , s
   REAL(REAL64) , EXTERNAL :: distdot
   INTEGER , SAVE :: hess , i , idx , ii , k , m , p2 , ptr , vc , vrn , vs
   LOGICAL , SAVE :: lp , rp
   EXTERNAL bisinit , givens , mgsro
!
! End of declarations rewritten by SPAG
!
!-----------------------------------------------------------------------
!     This a version of GMRES implemented with reverse communication.
!     It is a simple restart version of the GMRES algorithm.
!
!     ipar(5) == the dimension of the Krylov subspace
!     after every ipar(5) iterations, the GMRES will restart with
!     the updated solution and recomputed residual vector.
!
!     the space of the `w' is used as follows:
!     (1) the basis for the Krylov subspace, size n*(m+1);
!     (2) the Hessenberg matrix, only the upper triangular
!     portion of the matrix is stored, size (m+1)*m/2 + 1
!     (3) three vectors, all are of size m, they are
!     the cosine and sine of the Givens rotations, the third one holds
!     the residuals, it is of size m+1.
!
!     TOTAL SIZE REQUIRED == (n+3)*(m+2) + (m+1)*m/2
!     Note: m == ipar(5). The default value for this is 15 if
!     ipar(5) <= 1.
!-----------------------------------------------------------------------
!     external functions used
!
!
!
!     local variables, ptr and p2 are temporary pointers,
!     hess points to the Hessenberg matrix,
!     vc, vs point to the cosines and sines of the Givens rotations
!     vrn points to the vectors of residual norms, more precisely
!     the right hand side of the least square problem solved.
!
   INTEGER :: spag_nextblock_1
   spag_nextblock_1 = 1
   SPAG_DispatchLoop_1: DO
      SELECT CASE (spag_nextblock_1)
      CASE (1)
!
!     check the status of the call
!
         IF ( Ipar(1)<=0 ) Ipar(10) = 0
         IF ( Ipar(10)==1 ) THEN
            Ipar(7) = Ipar(7) + 1
            Ipar(13) = Ipar(13) + 1
            IF ( lp ) THEN
               DO i = 1 , N
                  W(N+i) = Rhs(i) - W(i)
               ENDDO
               Ipar(1) = 3
               Ipar(10) = 2
               RETURN
            ELSE
               DO i = 1 , N
                  W(i) = Rhs(i) - W(i)
               ENDDO
            ENDIF
            Fpar(11) = Fpar(11) + N
            spag_nextblock_1 = 3
            CYCLE SPAG_DispatchLoop_1
         ELSEIF ( Ipar(10)==2 ) THEN
            spag_nextblock_1 = 3
            CYCLE SPAG_DispatchLoop_1
         ELSEIF ( Ipar(10)==3 ) THEN
            spag_nextblock_1 = 5
            CYCLE SPAG_DispatchLoop_1
         ELSEIF ( Ipar(10)==4 ) THEN
!
            IF ( lp ) THEN
               Ipar(1) = 3
               Ipar(8) = Ipar(9)
               Ipar(9) = k*N + 1
               Ipar(10) = 5
               RETURN
            ENDIF
         ELSEIF ( Ipar(10)==5 ) THEN
         ELSEIF ( Ipar(10)==6 ) THEN
            spag_nextblock_1 = 7
            CYCLE SPAG_DispatchLoop_1
         ELSEIF ( Ipar(10)==7 ) THEN
            spag_nextblock_1 = 8
            CYCLE SPAG_DispatchLoop_1
         ELSE
!
!     initialization
!
            IF ( Ipar(5)<=1 ) THEN
               m = 15
            ELSE
               m = Ipar(5)
            ENDIF
            idx = N*(m+1)
            hess = idx + N
            vc = hess + (m+1)*m/2 + 1
            vs = vc + m
            vrn = vs + m
            i = vrn + m + 1
            CALL bisinit(Ipar,Fpar,i,1,lp,rp,W)
            IF ( Ipar(1)<0 ) RETURN
            spag_nextblock_1 = 2
            CYCLE SPAG_DispatchLoop_1
         ENDIF
!
!     Modified Gram-Schmidt orthogonalization procedure
!     temporary pointer 'ptr' is pointing to the current column of the
!     Hessenberg matrix. 'p2' points to the new basis vector
!
         Ipar(7) = Ipar(7) + 1
         ptr = k*(k-1)/2 + hess
         p2 = Ipar(9)
         CALL mgsro(.FALSE.,N,N,k+1,k+1,Fpar(11),W,W(ptr+1),Ipar(12))
         IF ( Ipar(12)<0 ) THEN
            spag_nextblock_1 = 6
            CYCLE SPAG_DispatchLoop_1
         ENDIF
!
!     apply previous Givens rotations and generate a new one to eliminate
!     the subdiagonal element.
!
         p2 = ptr + 1
         DO i = 1 , k - 1
            ptr = p2
            p2 = p2 + 1
            alpha = W(ptr)
            c = W(vc+i)
            s = W(vs+i)
            W(ptr) = c*alpha + s*W(p2)
            W(p2) = c*W(p2) - s*alpha
         ENDDO
         CALL givens(W(p2),W(p2+1),c,s)
         W(vc+k) = c
         W(vs+k) = s
         p2 = vrn + k
         alpha = -s*W(p2)
         W(p2) = c*W(p2)
         W(p2+1) = alpha
!
!     end of one Arnoldi iteration, alpha will store the estimated
!     residual norm at current stage
!
         Fpar(11) = Fpar(11) + 6*k + 2
         alpha = abs(alpha)
         Fpar(5) = alpha
         IF ( .NOT.(k<m .AND. .NOT.(Ipar(3)>=0 .AND. alpha<=Fpar(4)) .AND. (Ipar(6)<=0 .OR. Ipar(7)<Ipar(6))) ) THEN
            spag_nextblock_1 = 6
            CYCLE SPAG_DispatchLoop_1
         ENDIF
         spag_nextblock_1 = 4
         CYCLE SPAG_DispatchLoop_1
      CASE (2)
!
!     request for matrix vector multiplication A*x in the initialization
!
         Ipar(1) = 1
         Ipar(8) = N + 1
         Ipar(9) = 1
         Ipar(10) = 1
         k = 0
         DO i = 1 , N
            W(N+i) = Sol(i)
         ENDDO
         RETURN
      CASE (3)
!
         alpha = sqrt(distdot(N,W,1,W,1))
         Fpar(11) = Fpar(11) + 2*N
         IF ( Ipar(7)==1 .AND. Ipar(3)/=999 ) THEN
            IF ( abs(Ipar(3))==2 ) THEN
               Fpar(4) = Fpar(1)*sqrt(distdot(N,Rhs,1,Rhs,1)) + Fpar(2)
               Fpar(11) = Fpar(11) + 2*N
            ELSE
               Fpar(4) = Fpar(1)*alpha + Fpar(2)
            ENDIF
            Fpar(3) = alpha
         ENDIF
         Fpar(5) = alpha
         W(vrn+1) = alpha
         IF ( alpha<=Fpar(4) .AND. Ipar(3)>=0 .AND. Ipar(3)/=999 ) THEN
            Ipar(1) = 0
            Fpar(6) = alpha
            spag_nextblock_1 = 9
            CYCLE SPAG_DispatchLoop_1
         ENDIF
         alpha = ONE/alpha
         DO ii = 1 , N
            W(ii) = alpha*W(ii)
         ENDDO
         Fpar(11) = Fpar(11) + N
         spag_nextblock_1 = 4
      CASE (4)
!
!     request for (1) right preconditioning
!     (2) matrix vector multiplication
!     (3) left preconditioning
!
         k = k + 1
         IF ( rp ) THEN
            Ipar(1) = 5
            Ipar(8) = k*N - N + 1
            IF ( lp ) THEN
               Ipar(9) = k*N + 1
            ELSE
               Ipar(9) = idx + 1
            ENDIF
            Ipar(10) = 3
            RETURN
         ENDIF
         spag_nextblock_1 = 5
      CASE (5)
!
         Ipar(1) = 1
         IF ( rp ) THEN
            Ipar(8) = Ipar(9)
         ELSE
            Ipar(8) = (k-1)*N + 1
         ENDIF
         IF ( lp ) THEN
            Ipar(9) = idx + 1
         ELSE
            Ipar(9) = 1 + k*N
         ENDIF
         Ipar(10) = 4
         RETURN
      CASE (6)
         DO
!
!     update the approximate solution, first solve the upper triangular
!     system, temporary pointer ptr points to the Hessenberg matrix,
!     p2 points to the right-hand-side (also the solution) of the system.
!
            ptr = hess + k*(k+1)/2
            p2 = vrn + k
            IF ( W(ptr)==ZERO ) THEN
!
!     if the diagonal elements of the last column is zero, reduce k by 1
!     so that a smaller trianguler system is solved [It should only
!     happen when the matrix is singular, and at most once!]
!
               k = k - 1
               IF ( k>0 ) CYCLE
               Ipar(1) = -3
               Ipar(12) = -4
               spag_nextblock_1 = 9
               CYCLE SPAG_DispatchLoop_1
            ENDIF
            W(p2) = W(p2)/W(ptr)
            DO i = k - 1 , 1 , -1
               ptr = ptr - i - 1
               DO ii = 1 , i
                  W(vrn+ii) = W(vrn+ii) - W(p2)*W(ptr+ii)
               ENDDO
               p2 = p2 - 1
               W(p2) = W(p2)/W(ptr)
            ENDDO
!
            DO ii = 1 , N
               W(ii) = W(ii)*W(p2)
            ENDDO
            DO i = 1 , k - 1
               ptr = i*N
               p2 = p2 + 1
               DO ii = 1 , N
                  W(ii) = W(ii) + W(p2)*W(ptr+ii)
               ENDDO
            ENDDO
            Fpar(11) = Fpar(11) + 2*k*N - N + k*(k+1)
!
            IF ( rp ) THEN
               Ipar(1) = 5
               Ipar(8) = 1
               Ipar(9) = idx + 1
               Ipar(10) = 6
               RETURN
            ENDIF
            spag_nextblock_1 = 7
            CYCLE SPAG_DispatchLoop_1
         ENDDO
         spag_nextblock_1 = 7
      CASE (7)
!
         IF ( rp ) THEN
            DO i = 1 , N
               Sol(i) = Sol(i) + W(idx+i)
            ENDDO
         ELSE
            DO i = 1 , N
               Sol(i) = Sol(i) + W(i)
            ENDDO
         ENDIF
         Fpar(11) = Fpar(11) + N
!
!     process the complete stopping criteria
!
         IF ( Ipar(3)==999 ) THEN
            Ipar(1) = 10
            Ipar(8) = -1
            Ipar(9) = idx + 1
            Ipar(10) = 7
            RETURN
         ELSEIF ( Ipar(3)<0 ) THEN
            IF ( Ipar(7)<=m+1 ) THEN
               Fpar(3) = abs(W(vrn+1))
               IF ( Ipar(3)==-1 ) Fpar(4) = Fpar(1)*Fpar(3) + Fpar(2)
            ENDIF
            Fpar(6) = abs(W(vrn+k))
         ELSE
            Fpar(6) = Fpar(5)
         ENDIF
         spag_nextblock_1 = 8
      CASE (8)
!
!     do we need to restart ?
!
         IF ( Ipar(12)/=0 ) THEN
            Ipar(1) = -3
            spag_nextblock_1 = 9
            CYCLE SPAG_DispatchLoop_1
         ENDIF
         IF ( (Ipar(7)<Ipar(6) .OR. Ipar(6)<=0) .AND. ((Ipar(3)==999 .AND. Ipar(11)==0) .OR. (Ipar(3)/=999 .AND. Fpar(6)>Fpar(4))) &
            & ) THEN
            spag_nextblock_1 = 2
            CYCLE SPAG_DispatchLoop_1
         ENDIF
!
!     termination, set error code, compute convergence rate
!
         IF ( Ipar(1)>0 ) THEN
            IF ( Ipar(3)==999 .AND. Ipar(11)==1 ) THEN
               Ipar(1) = 0
            ELSEIF ( Ipar(3)/=999 .AND. Fpar(6)<=Fpar(4) ) THEN
               Ipar(1) = 0
            ELSEIF ( Ipar(7)>=Ipar(6) .AND. Ipar(6)>0 ) THEN
               Ipar(1) = -1
            ELSE
               Ipar(1) = -10
            ENDIF
         ENDIF
         spag_nextblock_1 = 9
      CASE (9)
         IF ( Fpar(3)/=ZERO .AND. Fpar(6)/=ZERO .AND. Ipar(7)>Ipar(13) ) THEN
            Fpar(7) = log10(Fpar(3)/Fpar(6))/dble(Ipar(7)-Ipar(13))
         ELSE
            Fpar(7) = ZERO
         ENDIF
         EXIT SPAG_DispatchLoop_1
      END SELECT
   ENDDO SPAG_DispatchLoop_1
END SUBROUTINE gmres
!*==dqgmres.f90 processed by SPAG 8.04RA 15:47 21 Feb 2024
!!SPAG Open source Personal, Educational or Academic User Clemson University  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
!-----end-of-gmres
!-----------------------------------------------------------------------
SUBROUTINE dqgmres(N,Rhs,Sol,Ipar,Fpar,W)
   USE ISO_FORTRAN_ENV                 
   IMPLICIT NONE
!
! PARAMETER definitions rewritten by SPAG
!
   REAL(REAL64) , PARAMETER :: ONE = 1.0D0 , ZERO = 0.0D0 , DEPS = 1.0D-33
!
! Dummy argument declarations rewritten by SPAG
!
   INTEGER :: N
   REAL(REAL64) , DIMENSION(N) :: Rhs
   REAL(REAL64) , INTENT(INOUT) , DIMENSION(N) :: Sol
   INTEGER , INTENT(INOUT) , DIMENSION(16) :: Ipar
   REAL(REAL64) , INTENT(INOUT) , DIMENSION(16) :: Fpar
   REAL(REAL64) , INTENT(INOUT) , DIMENSION(*) :: W
!
! Local variable declarations rewritten by SPAG
!
   REAL(REAL64) , SAVE :: alpha , beta , c , psi , s
   REAL(REAL64) , EXTERNAL :: distdot
   LOGICAL , SAVE :: full , lp , rp
   INTEGER , SAVE :: i , ic , ihd , ihm , ii , is , iv , iw , j , j0 , jp1 , k , lb , ptr , ptrv , ptrw
   EXTERNAL bisinit , givens , mgsro
!
! End of declarations rewritten by SPAG
!
!-----------------------------------------------------------------------
!     DQGMRES -- Flexible Direct version of Quasi-General Minimum
!     Residual method. The right preconditioning can be varied from
!     step to step.
!
!     Work space used = n + lb * (2*n+4)
!     where lb = ipar(5) + 1 (default 16 if ipar(5) <= 1)
!-----------------------------------------------------------------------
!     local variables
!
!
!
!     where to go
!
   IF ( Ipar(1)<=0 ) Ipar(10) = 0
   IF ( Ipar(10)==1 ) THEN
!
      Ipar(7) = Ipar(7) + 1
      Ipar(13) = Ipar(13) + 1
      IF ( lp ) THEN
         DO i = 1 , N
            W(i) = Rhs(i) - W(i)
         ENDDO
         Ipar(1) = 3
         Ipar(8) = 1
         Ipar(9) = iv + 1
         Ipar(10) = 2
         RETURN
      ELSE
         DO i = 1 , N
            W(iv+i) = Rhs(i) - W(iv+i)
         ENDDO
      ENDIF
      Fpar(11) = Fpar(11) + N
   ELSEIF ( Ipar(10)==2 ) THEN
   ELSEIF ( Ipar(10)==3 ) THEN
      CALL spag_block_2
      RETURN
   ELSEIF ( Ipar(10)==4 ) THEN
!
      IF ( lp ) THEN
         Ipar(1) = 3
         Ipar(8) = Ipar(9)
         Ipar(9) = ptrw
         Ipar(10) = 5
         RETURN
      ENDIF
      CALL spag_block_3
      RETURN
   ELSEIF ( Ipar(10)==5 ) THEN
      CALL spag_block_3
      RETURN
   ELSEIF ( Ipar(10)==6 ) THEN
      CALL spag_block_4
      RETURN
   ELSE
!
!     locations of the work arrays. The arrangement is as follows:
!     w(1:n) -- temporary storage for the results of the preconditioning
!     w(iv+1:iw) -- the V's
!     w(iw+1:ic) -- the W's
!     w(ic+1:is) -- the COSINEs of the Givens rotations
!     w(is+1:ihm) -- the SINEs of the Givens rotations
!     w(ihm+1:ihd) -- the last column of the Hessenberg matrix
!     w(ihd+1:i) -- the inverse of the diagonals of the Hessenberg matrix
!
      IF ( Ipar(5)<=1 ) THEN
         lb = 16
      ELSE
         lb = Ipar(5) + 1
      ENDIF
      iv = N
      iw = iv + lb*N
      ic = iw + lb*N
      is = ic + lb
      ihm = is + lb
      ihd = ihm + lb
      i = ihd + lb
!
!     parameter check, initializations
!
      full = .FALSE.
      CALL bisinit(Ipar,Fpar,i,1,lp,rp,W)
      IF ( Ipar(1)<0 ) RETURN
      Ipar(1) = 1
      IF ( lp ) THEN
         DO ii = 1 , N
            W(iv+ii) = Sol(ii)
         ENDDO
         Ipar(8) = iv + 1
         Ipar(9) = 1
      ELSE
         DO ii = 1 , N
            W(ii) = Sol(ii)
         ENDDO
         Ipar(8) = 1
         Ipar(9) = iv + 1
      ENDIF
      Ipar(10) = 1
      RETURN
   ENDIF
!
   alpha = sqrt(distdot(N,W(iv+1),1,W(iv+1),1))
   Fpar(11) = Fpar(11) + (N+N)
   IF ( abs(Ipar(3))==2 ) THEN
      Fpar(4) = Fpar(1)*sqrt(distdot(N,Rhs,1,Rhs,1)) + Fpar(2)
      Fpar(11) = Fpar(11) + 2*N
   ELSEIF ( Ipar(3)/=999 ) THEN
      Fpar(4) = Fpar(1)*alpha + Fpar(2)
   ENDIF
   Fpar(3) = alpha
   Fpar(5) = alpha
   psi = alpha
   IF ( alpha<=Fpar(4) ) THEN
      Ipar(1) = 0
      Fpar(6) = alpha
      CALL spag_block_5
      RETURN
   ENDIF
   alpha = ONE/alpha
   DO i = 1 , N
      W(iv+i) = W(iv+i)*alpha
   ENDDO
   Fpar(11) = Fpar(11) + N
   j = 0
   CALL spag_block_1
CONTAINS
   SUBROUTINE spag_block_1
!
!     iterations start here
!
      j = j + 1
      IF ( j>lb ) j = j - lb
      jp1 = j + 1
      IF ( jp1>lb ) jp1 = jp1 - lb
      ptrv = iv + (j-1)*N + 1
      ptrw = iv + (jp1-1)*N + 1
      IF ( .NOT.full ) THEN
         IF ( j>jp1 ) full = .TRUE.
      ENDIF
      IF ( full ) THEN
         j0 = jp1 + 1
         IF ( j0>lb ) j0 = j0 - lb
      ELSE
         j0 = 1
      ENDIF
!
!     request the caller to perform matrix-vector multiplication and
!     preconditioning
!
      IF ( rp ) THEN
         Ipar(1) = 5
         Ipar(8) = ptrv
         Ipar(9) = ptrv + iw - iv
         Ipar(10) = 3
         RETURN
      ELSE
         DO i = 0 , N - 1
            W(ptrv+iw-iv+i) = W(ptrv+i)
         ENDDO
      ENDIF
      CALL spag_block_2
   END SUBROUTINE spag_block_1
   SUBROUTINE spag_block_2
!
      Ipar(1) = 1
      IF ( rp ) THEN
         Ipar(8) = Ipar(9)
      ELSE
         Ipar(8) = ptrv
      ENDIF
      IF ( lp ) THEN
         Ipar(9) = 1
      ELSE
         Ipar(9) = ptrw
      ENDIF
      Ipar(10) = 4
      RETURN
   END SUBROUTINE spag_block_2
   SUBROUTINE spag_block_3
!
!     compute the last column of the Hessenberg matrix
!     modified Gram-schmidt procedure, orthogonalize against (lb-1)
!     previous vectors
!
      CALL mgsro(full,N,N,lb,jp1,Fpar(11),W(iv+1),W(ihm+1),Ipar(12))
      IF ( Ipar(12)<0 ) THEN
         Ipar(1) = -3
         CALL spag_block_5
         RETURN
      ENDIF
      beta = W(ihm+jp1)
!
!     incomplete factorization (QR factorization through Givens rotations)
!     (1) apply previous rotations [(lb-1) of them]
!     (2) generate a new rotation
!
      IF ( full ) THEN
         W(ihm+jp1) = W(ihm+j0)*W(is+jp1)
         W(ihm+j0) = W(ihm+j0)*W(ic+jp1)
      ENDIF
      i = j0
      DO WHILE ( i/=j )
         k = i + 1
         IF ( k>lb ) k = k - lb
         c = W(ic+i)
         s = W(is+i)
         alpha = W(ihm+i)
         W(ihm+i) = c*alpha + s*W(ihm+k)
         W(ihm+k) = c*W(ihm+k) - s*alpha
         i = k
      ENDDO
      CALL givens(W(ihm+j),beta,c,s)
      IF ( full ) THEN
         Fpar(11) = Fpar(11) + 6*lb
      ELSE
         Fpar(11) = Fpar(11) + 6*j
      ENDIF
!
!     detect whether diagonal element of this column is zero
!
      IF ( abs(W(ihm+j))<DEPS ) THEN
         Ipar(1) = -3
         CALL spag_block_5
         RETURN
      ENDIF
      W(ihd+j) = ONE/W(ihm+j)
      W(ic+j) = c
      W(is+j) = s
!
!     update the W's (the conjugate directions) -- essentially this is one
!     step of triangular solve.
!
      ptrw = iw + (j-1)*N + 1
      IF ( full ) THEN
         DO i = j + 1 , lb
            alpha = -W(ihm+i)*W(ihd+i)
            ptr = iw + (i-1)*N + 1
            DO ii = 0 , N - 1
               W(ptrw+ii) = W(ptrw+ii) + alpha*W(ptr+ii)
            ENDDO
         ENDDO
      ENDIF
      DO i = 1 , j - 1
         alpha = -W(ihm+i)*W(ihd+i)
         ptr = iw + (i-1)*N + 1
         DO ii = 0 , N - 1
            W(ptrw+ii) = W(ptrw+ii) + alpha*W(ptr+ii)
         ENDDO
      ENDDO
!
!     update the solution to the linear system
!
      alpha = psi*c*W(ihd+j)
      psi = -s*psi
      DO i = 1 , N
         Sol(i) = Sol(i) + alpha*W(ptrw-1+i)
      ENDDO
      IF ( full ) THEN
         Fpar(11) = Fpar(11) + lb*(N+N)
      ELSE
         Fpar(11) = Fpar(11) + j*(N+N)
      ENDIF
!
!     determine whether to continue,
!     compute the desired error/residual norm
!
      Ipar(7) = Ipar(7) + 1
      Fpar(5) = abs(psi)
      IF ( Ipar(3)==999 ) THEN
         Ipar(1) = 10
         Ipar(8) = -1
         Ipar(9) = 1
         Ipar(10) = 6
         RETURN
      ENDIF
      IF ( Ipar(3)<0 ) THEN
         alpha = abs(alpha)
         IF ( Ipar(7)==2 .AND. Ipar(3)==-1 ) THEN
            Fpar(3) = alpha*sqrt(distdot(N,W(ptrw),1,W(ptrw),1))
            Fpar(4) = Fpar(1)*Fpar(3) + Fpar(2)
            Fpar(6) = Fpar(3)
         ELSE
            Fpar(6) = alpha*sqrt(distdot(N,W(ptrw),1,W(ptrw),1))
         ENDIF
         Fpar(11) = Fpar(11) + 2*N
      ELSE
         Fpar(6) = Fpar(5)
      ENDIF
      IF ( Ipar(1)>=0 .AND. Fpar(6)>Fpar(4) .AND. (Ipar(6)<=0 .OR. Ipar(7)<Ipar(6)) ) THEN
         CALL spag_block_1
         RETURN
      ENDIF
      CALL spag_block_4
   END SUBROUTINE spag_block_3
   SUBROUTINE spag_block_4
      IF ( Ipar(3)==999 .AND. Ipar(11)==0 ) THEN
         CALL spag_block_1
         RETURN
      ENDIF
      CALL spag_block_5
   END SUBROUTINE spag_block_4
   SUBROUTINE spag_block_5
!
!     clean up the iterative solver
!
      Fpar(7) = ZERO
      IF ( Fpar(3)/=ZERO .AND. Fpar(6)/=ZERO .AND. Ipar(7)>Ipar(13) ) Fpar(7) = log10(Fpar(3)/Fpar(6))/dble(Ipar(7)-Ipar(13))
      IF ( Ipar(1)>0 ) THEN
         IF ( Ipar(3)==999 .AND. Ipar(11)/=0 ) THEN
            Ipar(1) = 0
         ELSEIF ( Fpar(6)<=Fpar(4) ) THEN
            Ipar(1) = 0
         ELSEIF ( Ipar(6)>0 .AND. Ipar(7)>=Ipar(6) ) THEN
            Ipar(1) = -1
         ELSE
            Ipar(1) = -10
         ENDIF
      ENDIF
   END SUBROUTINE spag_block_5
END SUBROUTINE dqgmres
!*==fgmres.f90 processed by SPAG 8.04RA 15:47 21 Feb 2024
!!SPAG Open source Personal, Educational or Academic User Clemson University  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
!-----end-of-dqgmres
!-----------------------------------------------------------------------
SUBROUTINE fgmres(N,Rhs,Sol,Ipar,Fpar,W)
   USE ISO_FORTRAN_ENV                 
   IMPLICIT NONE
!
! PARAMETER definitions rewritten by SPAG
!
   REAL(REAL64) , PARAMETER :: ONE = 1.0D0 , ZERO = 0.0D0
!
! Dummy argument declarations rewritten by SPAG
!
   INTEGER :: N
   REAL(REAL64) , DIMENSION(N) :: Rhs
   REAL(REAL64) , INTENT(INOUT) , DIMENSION(N) :: Sol
   INTEGER , INTENT(INOUT) , DIMENSION(16) :: Ipar
   REAL(REAL64) , INTENT(INOUT) , DIMENSION(16) :: Fpar
   REAL(REAL64) , INTENT(INOUT) , DIMENSION(*) :: W
!
! Local variable declarations rewritten by SPAG
!
   REAL(REAL64) , SAVE :: alpha , c , s
   REAL(REAL64) , EXTERNAL :: distdot
   INTEGER , SAVE :: hess , i , idx , ii , iz , k , m , p2 , ptr , vc , vrn , vs
   LOGICAL , SAVE :: lp , rp
   EXTERNAL bisinit , givens , mgsro
!
! End of declarations rewritten by SPAG
!
!-----------------------------------------------------------------------
!     This a version of FGMRES implemented with reverse communication.
!
!     ipar(5) == the dimension of the Krylov subspace
!
!     the space of the `w' is used as follows:
!     >> V: the bases for the Krylov subspace, size n*(m+1);
!     >> W: the above bases after (left-)multiplying with the
!     right-preconditioner inverse, size m*n;
!     >> a temporary vector of size n;
!     >> the Hessenberg matrix, only the upper triangular portion
!     of the matrix is stored, size (m+1)*m/2 + 1
!     >> three vectors, first two are of size m, they are the cosine
!     and sine of the Givens rotations, the third one holds the
!     residuals, it is of size m+1.
!
!     TOTAL SIZE REQUIRED == n*(2m+1) + (m+1)*m/2 + 3*m + 2
!     Note: m == ipar(5). The default value for this is 15 if
!     ipar(5) <= 1.
!-----------------------------------------------------------------------
!     external functions used
!
!
!
!     local variables, ptr and p2 are temporary pointers,
!     hess points to the Hessenberg matrix,
!     vc, vs point to the cosines and sines of the Givens rotations
!     vrn points to the vectors of residual norms, more precisely
!     the right hand side of the least square problem solved.
!
   INTEGER :: spag_nextblock_1
   spag_nextblock_1 = 1
   SPAG_DispatchLoop_1: DO
      SELECT CASE (spag_nextblock_1)
      CASE (1)
!
!     check the status of the call
!
         IF ( Ipar(1)<=0 ) Ipar(10) = 0
         IF ( Ipar(10)==1 ) THEN
            Ipar(7) = Ipar(7) + 1
            Ipar(13) = Ipar(13) + 1
            Fpar(11) = Fpar(11) + N
            IF ( lp ) THEN
               DO i = 1 , N
                  W(N+i) = Rhs(i) - W(i)
               ENDDO
               Ipar(1) = 3
               Ipar(10) = 2
               RETURN
            ELSE
               DO i = 1 , N
                  W(i) = Rhs(i) - W(i)
               ENDDO
            ENDIF
            spag_nextblock_1 = 3
            CYCLE SPAG_DispatchLoop_1
         ELSEIF ( Ipar(10)==2 ) THEN
            spag_nextblock_1 = 3
            CYCLE SPAG_DispatchLoop_1
         ELSEIF ( Ipar(10)==3 ) THEN
            spag_nextblock_1 = 5
            CYCLE SPAG_DispatchLoop_1
         ELSEIF ( Ipar(10)==4 ) THEN
!
            IF ( lp ) THEN
               Ipar(1) = 3
               Ipar(8) = Ipar(9)
               Ipar(9) = k*N + 1
               Ipar(10) = 5
               RETURN
            ENDIF
         ELSEIF ( Ipar(10)==5 ) THEN
         ELSEIF ( Ipar(10)==6 ) THEN
            spag_nextblock_1 = 7
            CYCLE SPAG_DispatchLoop_1
         ELSE
!
!     initialization
!
            IF ( Ipar(5)<=1 ) THEN
               m = 15
            ELSE
               m = Ipar(5)
            ENDIF
            idx = N*(m+1)
            iz = idx + N
            hess = iz + N*m
            vc = hess + (m+1)*m/2 + 1
            vs = vc + m
            vrn = vs + m
            i = vrn + m + 1
            CALL bisinit(Ipar,Fpar,i,1,lp,rp,W)
            IF ( Ipar(1)<0 ) RETURN
            spag_nextblock_1 = 2
            CYCLE SPAG_DispatchLoop_1
         ENDIF
!
!     Modified Gram-Schmidt orthogonalization procedure
!     temporary pointer 'ptr' is pointing to the current column of the
!     Hessenberg matrix. 'p2' points to the new basis vector
!
         ptr = k*(k-1)/2 + hess
         p2 = Ipar(9)
         Ipar(7) = Ipar(7) + 1
         CALL mgsro(.FALSE.,N,N,k+1,k+1,Fpar(11),W,W(ptr+1),Ipar(12))
         IF ( Ipar(12)<0 ) THEN
            spag_nextblock_1 = 6
            CYCLE SPAG_DispatchLoop_1
         ENDIF
!
!     apply previous Givens rotations and generate a new one to eliminate
!     the subdiagonal element.
!
         p2 = ptr + 1
         DO i = 1 , k - 1
            ptr = p2
            p2 = p2 + 1
            alpha = W(ptr)
            c = W(vc+i)
            s = W(vs+i)
            W(ptr) = c*alpha + s*W(p2)
            W(p2) = c*W(p2) - s*alpha
         ENDDO
         CALL givens(W(p2),W(p2+1),c,s)
         W(vc+k) = c
         W(vs+k) = s
         p2 = vrn + k
         alpha = -s*W(p2)
         W(p2) = c*W(p2)
         W(p2+1) = alpha
         Fpar(11) = Fpar(11) + 6*k
!
!     end of one Arnoldi iteration, alpha will store the estimated
!     residual norm at current stage
!
         alpha = abs(alpha)
         Fpar(5) = alpha
         IF ( .NOT.(k<m .AND. .NOT.(Ipar(3)>=0 .AND. alpha<=Fpar(4)) .AND. (Ipar(6)<=0 .OR. Ipar(7)<Ipar(6))) ) THEN
            spag_nextblock_1 = 6
            CYCLE SPAG_DispatchLoop_1
         ENDIF
         spag_nextblock_1 = 4
         CYCLE SPAG_DispatchLoop_1
      CASE (2)
!
!     request for matrix vector multiplication A*x in the initialization
!
         Ipar(1) = 1
         Ipar(8) = N + 1
         Ipar(9) = 1
         Ipar(10) = 1
         k = 0
         DO ii = 1 , N
            W(ii+N) = Sol(ii)
         ENDDO
         RETURN
      CASE (3)
!
         alpha = sqrt(distdot(N,W,1,W,1))
         Fpar(11) = Fpar(11) + N + N
         IF ( Ipar(7)==1 .AND. Ipar(3)/=999 ) THEN
            IF ( abs(Ipar(3))==2 ) THEN
               Fpar(4) = Fpar(1)*sqrt(distdot(N,Rhs,1,Rhs,1)) + Fpar(2)
               Fpar(11) = Fpar(11) + 2*N
            ELSE
               Fpar(4) = Fpar(1)*alpha + Fpar(2)
            ENDIF
            Fpar(3) = alpha
         ENDIF
         Fpar(5) = alpha
         W(vrn+1) = alpha
         IF ( alpha<=Fpar(4) .AND. Ipar(3)>=0 .AND. Ipar(3)/=999 ) THEN
            Ipar(1) = 0
            Fpar(6) = alpha
            spag_nextblock_1 = 8
            CYCLE SPAG_DispatchLoop_1
         ENDIF
         alpha = ONE/alpha
         DO ii = 1 , N
            W(ii) = W(ii)*alpha
         ENDDO
         Fpar(11) = Fpar(11) + N
         spag_nextblock_1 = 4
      CASE (4)
!
!     request for (1) right preconditioning
!     (2) matrix vector multiplication
!     (3) left preconditioning
!
         k = k + 1
         IF ( rp ) THEN
            Ipar(1) = 5
            Ipar(8) = k*N - N + 1
            Ipar(9) = iz + Ipar(8)
            Ipar(10) = 3
            RETURN
         ELSE
            DO ii = 0 , N - 1
               W(iz+k*N-ii) = W(k*N-ii)
            ENDDO
         ENDIF
         spag_nextblock_1 = 5
      CASE (5)
!
         Ipar(1) = 1
         IF ( rp ) THEN
            Ipar(8) = Ipar(9)
         ELSE
            Ipar(8) = (k-1)*N + 1
         ENDIF
         IF ( lp ) THEN
            Ipar(9) = idx + 1
         ELSE
            Ipar(9) = 1 + k*N
         ENDIF
         Ipar(10) = 4
         RETURN
      CASE (6)
         DO
!
!     update the approximate solution, first solve the upper triangular
!     system, temporary pointer ptr points to the Hessenberg matrix,
!     p2 points to the right-hand-side (also the solution) of the system.
!
            ptr = hess + k*(k+1)/2
            p2 = vrn + k
            IF ( W(ptr)==ZERO ) THEN
!
!     if the diagonal elements of the last column is zero, reduce k by 1
!     so that a smaller trianguler system is solved [It should only
!     happen when the matrix is singular!]
!
               k = k - 1
               IF ( k>0 ) CYCLE
               Ipar(1) = -3
               Ipar(12) = -4
               spag_nextblock_1 = 8
               CYCLE SPAG_DispatchLoop_1
            ENDIF
            W(p2) = W(p2)/W(ptr)
            DO i = k - 1 , 1 , -1
               ptr = ptr - i - 1
               DO ii = 1 , i
                  W(vrn+ii) = W(vrn+ii) - W(p2)*W(ptr+ii)
               ENDDO
               p2 = p2 - 1
               W(p2) = W(p2)/W(ptr)
            ENDDO
!
            DO i = 0 , k - 1
               ptr = iz + i*N
               DO ii = 1 , N
                  Sol(ii) = Sol(ii) + W(p2)*W(ptr+ii)
               ENDDO
               p2 = p2 + 1
            ENDDO
            Fpar(11) = Fpar(11) + 2*k*N + k*(k+1)
!
!     process the complete stopping criteria
!
            IF ( Ipar(3)==999 ) THEN
               Ipar(1) = 10
               Ipar(8) = -1
               Ipar(9) = idx + 1
               Ipar(10) = 6
               RETURN
            ELSEIF ( Ipar(3)<0 ) THEN
               IF ( Ipar(7)<=m+1 ) THEN
                  Fpar(3) = abs(W(vrn+1))
                  IF ( Ipar(3)==-1 ) Fpar(4) = Fpar(1)*Fpar(3) + Fpar(2)
               ENDIF
               Fpar(6) = abs(W(vrn+k))
            ELSEIF ( Ipar(3)/=999 ) THEN
               Fpar(6) = Fpar(5)
            ENDIF
            spag_nextblock_1 = 7
            CYCLE SPAG_DispatchLoop_1
         ENDDO
         spag_nextblock_1 = 7
      CASE (7)
!
!     do we need to restart ?
!
         IF ( Ipar(12)/=0 ) THEN
            Ipar(1) = -3
            spag_nextblock_1 = 8
            CYCLE SPAG_DispatchLoop_1
         ENDIF
         IF ( (Ipar(7)<Ipar(6) .OR. Ipar(6)<=0) .AND. ((Ipar(3)==999 .AND. Ipar(11)==0) .OR. (Ipar(3)/=999 .AND. Fpar(6)>Fpar(4))) &
            & ) THEN
            spag_nextblock_1 = 2
            CYCLE SPAG_DispatchLoop_1
         ENDIF
!
!     termination, set error code, compute convergence rate
!
         IF ( Ipar(1)>0 ) THEN
            IF ( Ipar(3)==999 .AND. Ipar(11)==1 ) THEN
               Ipar(1) = 0
            ELSEIF ( Ipar(3)/=999 .AND. Fpar(6)<=Fpar(4) ) THEN
               Ipar(1) = 0
            ELSEIF ( Ipar(7)>=Ipar(6) .AND. Ipar(6)>0 ) THEN
               Ipar(1) = -1
            ELSE
               Ipar(1) = -10
            ENDIF
         ENDIF
         spag_nextblock_1 = 8
      CASE (8)
         IF ( Fpar(3)/=ZERO .AND. Fpar(6)/=ZERO .AND. Ipar(7)>Ipar(13) ) THEN
            Fpar(7) = log10(Fpar(3)/Fpar(6))/dble(Ipar(7)-Ipar(13))
         ELSE
            Fpar(7) = ZERO
         ENDIF
         EXIT SPAG_DispatchLoop_1
      END SELECT
   ENDDO SPAG_DispatchLoop_1
END SUBROUTINE fgmres
!*==dbcg.f90 processed by SPAG 8.04RA 15:47 21 Feb 2024
!!SPAG Open source Personal, Educational or Academic User Clemson University  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
!-----end-of-fgmres
!-----------------------------------------------------------------------
SUBROUTINE dbcg(N,Rhs,Sol,Ipar,Fpar,W)
   USE ISO_FORTRAN_ENV                 
   IMPLICIT NONE
!
! PARAMETER definitions rewritten by SPAG
!
   REAL(REAL64) , PARAMETER :: ONE = 1.0D0 , ZERO = 0.0D0
!
! Dummy argument declarations rewritten by SPAG
!
   INTEGER :: N
   REAL(REAL64) , DIMENSION(N) :: Rhs
   REAL(REAL64) , DIMENSION(N) :: Sol
   INTEGER , INTENT(INOUT) , DIMENSION(16) :: Ipar
   REAL(REAL64) , INTENT(INOUT) , DIMENSION(16) :: Fpar
   REAL(REAL64) , INTENT(INOUT) , DIMENSION(N,*) :: W
!
! Local variable declarations rewritten by SPAG
!
   REAL(REAL64) , SAVE :: beta , delta , res , ss , ss1 , t , umm , x , zeta
   REAL(REAL64) , EXTERNAL :: distdot
   LOGICAL , SAVE :: full , lp , rp
   INTEGER , SAVE :: i , i2 , indp , ip2 , j , ju , k , lb , lbm1 , np
   LOGICAL , DIMENSION(3) , SAVE :: perm
   REAL(REAL64) , DIMENSION(3) , SAVE :: u , usav , ypiv
   EXTERNAL bisinit , implu , tidycg , uppdir
!
! End of declarations rewritten by SPAG
!
!-----------------------------------------------------------------------
! Quasi GMRES method for solving a linear
! system of equations a * sol = y.  double precision version.
! this version is without restarting and without preconditioning.
! parameters :
! -----------
! n     = dimension of the problem
!
! y     = w(:,1) a temporary storage used for various operations
! z     = w(:,2) a work vector of length n.
! v     = w(:,3:4) size n x 2
! w     = w(:,5:6) size n x 2
! p     = w(:,7:9) work array of dimension n x 3
! del x = w(:,10)  accumulation of the changes in solution
! tmp   = w(:,11)  a temporary vector used to hold intermediate result of
!                  preconditioning, etc.
!
! sol   = the solution of the problem . at input sol must contain an
!         initial guess to the solution.
!    ***  note:   y is destroyed on return.
!
!-----------------------------------------------------------------------
! subroutines and functions called:
! 1) matrix vector multiplication and preconditioning through reverse
!     communication
!
! 2) implu, uppdir, distdot (blas)
!-----------------------------------------------------------------------
! aug. 1983  version.    author youcef saad. yale university computer
! science dept. some  changes made july 3, 1986.
! references: siam j. sci. stat. comp., vol. 5, pp. 203-228 (1984)
!-----------------------------------------------------------------------
!     local variables
!
!
!
!     where to go
!
   IF ( Ipar(1)<=0 ) Ipar(10) = 0
   IF ( Ipar(10)==1 ) THEN
      Ipar(7) = Ipar(7) + 1
      Ipar(13) = Ipar(13) + 1
      IF ( lp ) THEN
         DO i = 1 , N
            W(i,1) = Rhs(i) - W(i,2)
         ENDDO
         Ipar(1) = 3
         Ipar(8) = 1
         Ipar(9) = N + N + 1
         Ipar(10) = 2
         RETURN
      ELSE
         DO i = 1 , N
            W(i,3) = Rhs(i) - W(i,2)
         ENDDO
      ENDIF
      Fpar(11) = Fpar(11) + N
   ELSEIF ( Ipar(10)==2 ) THEN
   ELSEIF ( Ipar(10)==3 ) THEN
      CALL spag_block_2
      RETURN
   ELSEIF ( Ipar(10)==4 ) THEN
!
      IF ( lp ) THEN
         Ipar(1) = 3
         Ipar(8) = Ipar(9)
         Ipar(9) = 1
         Ipar(10) = 5
         RETURN
      ENDIF
      CALL spag_block_3
      RETURN
   ELSEIF ( Ipar(10)==5 ) THEN
      CALL spag_block_3
      RETURN
   ELSEIF ( Ipar(10)==6 ) THEN
      CALL spag_block_4
      RETURN
   ELSEIF ( Ipar(10)==7 ) THEN
!
      IF ( rp ) THEN
         Ipar(1) = 6
         Ipar(8) = Ipar(9)
         Ipar(9) = N + 1
         Ipar(10) = 8
         RETURN
      ENDIF
      CALL spag_block_5
      RETURN
   ELSEIF ( Ipar(10)==8 ) THEN
      CALL spag_block_5
      RETURN
   ELSEIF ( Ipar(10)==9 ) THEN
      CALL spag_block_6
      RETURN
   ELSEIF ( Ipar(10)==10 ) THEN
      CALL spag_block_8
      RETURN
   ELSE
!
!     initialization, parameter checking, clear the work arrays
!
      CALL bisinit(Ipar,Fpar,11*N,1,lp,rp,W)
      IF ( Ipar(1)<0 ) RETURN
      perm(1) = .FALSE.
      perm(2) = .FALSE.
      perm(3) = .FALSE.
      usav(1) = ZERO
      usav(2) = ZERO
      usav(3) = ZERO
      ypiv(1) = ZERO
      ypiv(2) = ZERO
      ypiv(3) = ZERO
!-----------------------------------------------------------------------
!     initialize constants for outer loop :
!-----------------------------------------------------------------------
      lb = 3
      lbm1 = 2
!
!     get initial residual vector and norm
!
      Ipar(1) = 1
      Ipar(8) = 1
      Ipar(9) = 1 + N
      DO i = 1 , N
         W(i,1) = Sol(i)
      ENDDO
      Ipar(10) = 1
      RETURN
   ENDIF
!
   Fpar(3) = sqrt(distdot(N,W(1,3),1,W(1,3),1))
   Fpar(11) = Fpar(11) + N + N
   Fpar(5) = Fpar(3)
   Fpar(7) = Fpar(3)
   zeta = Fpar(3)
   IF ( abs(Ipar(3))==2 ) THEN
      Fpar(4) = Fpar(1)*sqrt(distdot(N,Rhs,1,Rhs,1)) + Fpar(2)
      Fpar(11) = Fpar(11) + 2*N
   ELSEIF ( Ipar(3)/=999 ) THEN
      Fpar(4) = Fpar(1)*zeta + Fpar(2)
   ENDIF
   IF ( Ipar(3)>=0 .AND. Fpar(5)<=Fpar(4) ) THEN
      Fpar(6) = Fpar(5)
      CALL spag_block_7
      RETURN
   ENDIF
!
!     normalize first arnoldi vector
!
   t = ONE/zeta
   DO k = 1 , N
      W(k,3) = W(k,3)*t
      W(k,5) = W(k,3)
   ENDDO
   Fpar(11) = Fpar(11) + N
!
!     initialize constants for main loop
!
   beta = ZERO
   delta = ZERO
   i2 = 1
   indp = 0
   i = 0
   CALL spag_block_1
CONTAINS
   SUBROUTINE spag_block_1
!
!     main loop: i = index of the loop.
!
!-----------------------------------------------------------------------
      i = i + 1
!
      IF ( rp ) THEN
         Ipar(1) = 5
         Ipar(8) = (1+i2)*N + 1
         IF ( lp ) THEN
            Ipar(9) = 1
         ELSE
            Ipar(9) = 10*N + 1
         ENDIF
         Ipar(10) = 3
         RETURN
      ENDIF
      CALL spag_block_2
   END SUBROUTINE spag_block_1
   SUBROUTINE spag_block_2
!
      Ipar(1) = 1
      IF ( rp ) THEN
         Ipar(8) = Ipar(9)
      ELSE
         Ipar(8) = (1+i2)*N + 1
      ENDIF
      IF ( lp ) THEN
         Ipar(9) = 10*N + 1
      ELSE
         Ipar(9) = 1
      ENDIF
      Ipar(10) = 4
      RETURN
   END SUBROUTINE spag_block_2
   SUBROUTINE spag_block_3
!
!     A^t * x
!
      Ipar(7) = Ipar(7) + 1
      IF ( lp ) THEN
         Ipar(1) = 4
         Ipar(8) = (3+i2)*N + 1
         IF ( rp ) THEN
            Ipar(9) = N + 1
         ELSE
            Ipar(9) = 10*N + 1
         ENDIF
         Ipar(10) = 6
         RETURN
      ENDIF
      CALL spag_block_4
   END SUBROUTINE spag_block_3
   SUBROUTINE spag_block_4
!
      Ipar(1) = 2
      IF ( lp ) THEN
         Ipar(8) = Ipar(9)
      ELSE
         Ipar(8) = (3+i2)*N + 1
      ENDIF
      IF ( rp ) THEN
         Ipar(9) = 10*N + 1
      ELSE
         Ipar(9) = N + 1
      ENDIF
      Ipar(10) = 7
      RETURN
   END SUBROUTINE spag_block_4
   SUBROUTINE spag_block_5
!-----------------------------------------------------------------------
!     orthogonalize current v against previous v's and
!     determine relevant part of i-th column of u(.,.) the
!     upper triangular matrix --
!-----------------------------------------------------------------------
      Ipar(7) = Ipar(7) + 1
      u(1) = ZERO
      ju = 1
      k = i2
      IF ( i<=lbm1 ) ju = 0
      IF ( i<lb ) k = 0
      DO
         IF ( k==lbm1 ) k = 0
         k = k + 1
!
         IF ( k/=i2 ) THEN
            ss = delta
            ss1 = beta
            ju = ju + 1
            u(ju) = ss
         ELSE
            ss = distdot(N,W(1,1),1,W(1,4+k),1)
            Fpar(11) = Fpar(11) + 2*N
            ss1 = ss
            ju = ju + 1
            u(ju) = ss
         ENDIF
!
         DO j = 1 , N
            W(j,1) = W(j,1) - ss*W(j,k+2)
            W(j,2) = W(j,2) - ss1*W(j,k+4)
         ENDDO
         Fpar(11) = Fpar(11) + 4*N
!
         IF ( k==i2 ) THEN
!
!     end of Mod. Gram. Schmidt loop
!
            t = distdot(N,W(1,2),1,W(1,1),1)
!
            beta = sqrt(abs(t))
            delta = t/beta
!
            ss = ONE/beta
            ss1 = ONE/delta
!
!     normalize and insert new vectors
!
            ip2 = i2
            IF ( i2==lbm1 ) i2 = 0
            i2 = i2 + 1
!
            DO j = 1 , N
               W(j,i2+2) = W(j,1)*ss
               W(j,i2+4) = W(j,2)*ss1
            ENDDO
            Fpar(11) = Fpar(11) + 4*N
!-----------------------------------------------------------------------
!     end of orthogonalization.
!     now compute the coefficients u(k) of the last
!     column of the  l . u  factorization of h .
!-----------------------------------------------------------------------
            np = min0(i,lb)
            full = (i>=lb)
            CALL implu(np,umm,beta,ypiv,u,perm,full)
!-----------------------------------------------------------------------
!     update conjugate directions and solution
!-----------------------------------------------------------------------
            DO k = 1 , N
               W(k,1) = W(k,ip2+2)
            ENDDO
            CALL uppdir(N,W(1,7),np,lb,indp,W,u,usav,Fpar(11))
!-----------------------------------------------------------------------
            IF ( i/=1 ) THEN
               j = np - 1
               IF ( full ) j = j - 1
               IF ( .NOT.perm(j) ) zeta = -zeta*ypiv(j)
            ENDIF
            x = zeta/u(np)
            IF ( .NOT.(perm(np)) ) THEN
               DO k = 1 , N
                  W(k,10) = W(k,10) + x*W(k,1)
               ENDDO
               Fpar(11) = Fpar(11) + 2*N
            ENDIF
!-----------------------------------------------------------------------
            IF ( Ipar(3)==999 ) THEN
               Ipar(1) = 10
               Ipar(8) = 9*N + 1
               Ipar(9) = 10*N + 1
               Ipar(10) = 9
               RETURN
            ENDIF
            res = abs(beta*zeta/umm)
            Fpar(5) = res*sqrt(distdot(N,W(1,i2+2),1,W(1,i2+2),1))
            Fpar(11) = Fpar(11) + 2*N
            IF ( Ipar(3)<0 ) THEN
               Fpar(6) = x*sqrt(distdot(N,W,1,W,1))
               Fpar(11) = Fpar(11) + 2*N
               IF ( Ipar(7)<=3 ) THEN
                  Fpar(3) = Fpar(6)
                  IF ( Ipar(3)==-1 ) Fpar(4) = Fpar(1)*sqrt(Fpar(3)) + Fpar(2)
               ENDIF
            ELSE
               Fpar(6) = Fpar(5)
            ENDIF
            CALL spag_block_6
            RETURN
         ENDIF
      ENDDO
      CALL spag_block_6
   END SUBROUTINE spag_block_5
   SUBROUTINE spag_block_6
!---- convergence test -----------------------------------------------
      IF ( Ipar(3)==999 .AND. Ipar(11)==0 ) THEN
         CALL spag_block_1
         RETURN
      ENDIF
      IF ( Fpar(6)>Fpar(4) .AND. (Ipar(6)>Ipar(7) .OR. Ipar(6)<=0) ) THEN
         CALL spag_block_1
         RETURN
      ENDIF
!-----------------------------------------------------------------------
!     here the fact that the last step is different is accounted for.
!-----------------------------------------------------------------------
      IF ( perm(np) ) THEN
         x = zeta/umm
         DO k = 1 , N
            W(k,10) = W(k,10) + x*W(k,1)
         ENDDO
         Fpar(11) = Fpar(11) + 2*N
      ENDIF
      CALL spag_block_7
   END SUBROUTINE spag_block_6
   SUBROUTINE spag_block_7
!
!     right preconditioning and clean-up jobs
!
      IF ( rp ) THEN
         IF ( Ipar(1)<0 ) Ipar(12) = Ipar(1)
         Ipar(1) = 5
         Ipar(8) = 9*N + 1
         Ipar(9) = Ipar(8) + N
         Ipar(10) = 10
         RETURN
      ENDIF
      CALL spag_block_8
   END SUBROUTINE spag_block_7
   SUBROUTINE spag_block_8
      IF ( rp ) THEN
         CALL tidycg(N,Ipar,Fpar,Sol,W(1,11))
      ELSE
         CALL tidycg(N,Ipar,Fpar,Sol,W(1,10))
      ENDIF
   END SUBROUTINE spag_block_8
END SUBROUTINE dbcg
!*==implu.f90 processed by SPAG 8.04RA 15:47 21 Feb 2024
!!SPAG Open source Personal, Educational or Academic User Clemson University  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
!-----end-of-dbcg-------------------------------------------------------
!-----------------------------------------------------------------------
SUBROUTINE implu(Np,Umm,Beta,Ypiv,U,Permut,Full)
   USE ISO_FORTRAN_ENV                 
   IMPLICIT NONE
!
! Dummy argument declarations rewritten by SPAG
!
   INTEGER , INTENT(IN) :: Np
   REAL(REAL64) , INTENT(INOUT) :: Umm
   REAL(REAL64) , INTENT(IN) :: Beta
   REAL(REAL64) , INTENT(INOUT) , DIMENSION(*) :: Ypiv
   REAL(REAL64) , INTENT(INOUT) , DIMENSION(*) :: U
   LOGICAL , INTENT(INOUT) , DIMENSION(*) :: Permut
   LOGICAL , INTENT(IN) :: Full
!
! Local variable declarations rewritten by SPAG
!
   INTEGER :: k , npm1
   LOGICAL :: perm
   REAL(REAL64) :: x , xpiv
!
! End of declarations rewritten by SPAG
!
!-----------------------------------------------------------------------
!     performs implicitly one step of the lu factorization of a
!     banded hessenberg matrix.
!-----------------------------------------------------------------------
   npm1 = 0
   IF ( Np>1 ) THEN
      npm1 = Np - 1
!
!     -- perform  previous step of the factorization-
!
      DO k = 1 , npm1
         IF ( Permut(k) ) THEN
            x = U(k)
            U(k) = U(k+1)
            U(k+1) = x
         ENDIF
         U(k+1) = U(k+1) - Ypiv(k)*U(k)
      ENDDO
   ENDIF
!-----------------------------------------------------------------------
!     now determine pivotal information to be used in the next call
!-----------------------------------------------------------------------
   Umm = U(Np)
   perm = (Beta>abs(Umm))
   IF ( .NOT.perm ) THEN
      xpiv = Beta/Umm
   ELSE
      xpiv = Umm/Beta
      U(Np) = Beta
   ENDIF
   Permut(Np) = perm
   Ypiv(Np) = xpiv
   IF ( .NOT.Full ) RETURN
!     shift everything up if full...
   DO k = 1 , npm1
      Ypiv(k) = Ypiv(k+1)
      Permut(k) = Permut(k+1)
   ENDDO
!-----end-of-implu
END SUBROUTINE implu
!*==uppdir.f90 processed by SPAG 8.04RA 15:47 21 Feb 2024
!!SPAG Open source Personal, Educational or Academic User Clemson University  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
!-----------------------------------------------------------------------
SUBROUTINE uppdir(N,P,Np,Lbp,Indp,Y,U,Usav,Flops)
   USE ISO_FORTRAN_ENV                 
   IMPLICIT NONE
!
! PARAMETER definitions rewritten by SPAG
!
   REAL(REAL64) , PARAMETER :: ZERO = 0.0D0
!
! Dummy argument declarations rewritten by SPAG
!
   INTEGER , INTENT(IN) :: N
   INTEGER , INTENT(IN) :: Lbp
   REAL(REAL64) , INTENT(INOUT) , DIMENSION(N,Lbp) :: P
   INTEGER , INTENT(IN) :: Np
   INTEGER , INTENT(INOUT) :: Indp
   REAL(REAL64) , INTENT(INOUT) , DIMENSION(*) :: Y
   REAL(REAL64) , INTENT(IN) , DIMENSION(*) :: U
   REAL(REAL64) , INTENT(INOUT) , DIMENSION(*) :: Usav
   REAL(REAL64) , INTENT(INOUT) :: Flops
!
! Local variable declarations rewritten by SPAG
!
   INTEGER :: j , ju , k , npm1
   REAL(REAL64) :: x
!
! End of declarations rewritten by SPAG
!
!-----------------------------------------------------------------------
!     updates the conjugate directions p given the upper part of the
!     banded upper triangular matrix u.  u contains the non zero
!     elements of the column of the triangular matrix..
!-----------------------------------------------------------------------
!
   npm1 = Np - 1
   IF ( Np>1 ) THEN
      j = Indp
      ju = npm1
      DO
         IF ( j<=0 ) j = Lbp
         x = U(ju)/Usav(j)
         IF ( x/=ZERO ) THEN
            DO k = 1 , N
               Y(k) = Y(k) - x*P(k,j)
            ENDDO
            Flops = Flops + 2*N
         ENDIF
         j = j - 1
         ju = ju - 1
         IF ( ju<1 ) THEN
            CALL spag_block_1
            RETURN
         ENDIF
      ENDDO
   ENDIF
   CALL spag_block_1
CONTAINS
   SUBROUTINE spag_block_1
      Indp = Indp + 1
      IF ( Indp>Lbp ) Indp = 1
      Usav(Indp) = U(Np)
      DO k = 1 , N
         P(k,Indp) = Y(k)
      ENDDO
   END SUBROUTINE spag_block_1
!-----------------------------------------------------------------------
!-------end-of-uppdir---------------------------------------------------
END SUBROUTINE uppdir
!*==givens.f90 processed by SPAG 8.04RA 15:47 21 Feb 2024
!!SPAG Open source Personal, Educational or Academic User Clemson University  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
SUBROUTINE givens(X,Y,C,S)
   USE ISO_FORTRAN_ENV                 
   IMPLICIT NONE
!
! PARAMETER definitions rewritten by SPAG
!
   REAL(REAL64) , PARAMETER :: ZERO = 0.0D0 , ONE = 1.0D0
!
! Dummy argument declarations rewritten by SPAG
!
   REAL(REAL64) , INTENT(INOUT) :: X
   REAL(REAL64) , INTENT(INOUT) :: Y
   REAL(REAL64) , INTENT(INOUT) :: C
   REAL(REAL64) , INTENT(INOUT) :: S
!
! Local variable declarations rewritten by SPAG
!
   REAL(REAL64) :: t
!
! End of declarations rewritten by SPAG
!
!-----------------------------------------------------------------------
!     Given x and y, this subroutine generates a Givens' rotation c, s.
!     And apply the rotation on (x,y) ==> (sqrt(x**2 + y**2), 0).
!     (See P 202 of "matrix computation" by Golub and van Loan.)
!-----------------------------------------------------------------------
!
   IF ( X==ZERO .AND. Y==ZERO ) THEN
      C = ONE
      S = ZERO
   ELSEIF ( abs(Y)>abs(X) ) THEN
      t = X/Y
      X = sqrt(ONE+t*t)
      S = sign(ONE/X,Y)
      C = t*S
   ELSEIF ( abs(Y)<=abs(X) ) THEN
      t = Y/X
      Y = sqrt(ONE+t*t)
      C = sign(ONE/Y,X)
      S = t*C
   ELSE
!
!     X or Y must be an invalid floating-point number, set both to zero
!
      X = ZERO
      Y = ZERO
      C = ONE
      S = ZERO
   ENDIF
   X = abs(X*Y)
!
!     end of givens
!
END SUBROUTINE givens
!*==stopbis.f90 processed by SPAG 8.04RA 15:47 21 Feb 2024
!!SPAG Open source Personal, Educational or Academic User Clemson University  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
!-----end-of-givens
!-----------------------------------------------------------------------
FUNCTION stopbis(N,Ipar,Mvpi,Fpar,R,Delx,Sx)
   USE ISO_FORTRAN_ENV                 
   IMPLICIT NONE
!
! Function and Dummy argument declarations rewritten by SPAG
!
   INTEGER :: N
   LOGICAL :: stopbis
   INTEGER , INTENT(INOUT) , DIMENSION(16) :: Ipar
   INTEGER , INTENT(IN) :: Mvpi
   REAL(REAL64) , INTENT(INOUT) , DIMENSION(16) :: Fpar
   REAL(REAL64) , DIMENSION(N) :: R
   REAL(REAL64) , DIMENSION(N) :: Delx
   REAL(REAL64) , INTENT(IN) :: Sx
!
! Local variable declarations rewritten by SPAG
!
   REAL(REAL64) , EXTERNAL :: distdot
!
! End of declarations rewritten by SPAG
!
!-----------------------------------------------------------------------
!     function for determining the stopping criteria. return value of
!     true if the stopbis criteria is satisfied.
!-----------------------------------------------------------------------
   IF ( Ipar(11)==1 ) THEN
      stopbis = .TRUE.
   ELSE
      stopbis = .FALSE.
   ENDIF
   IF ( Ipar(6)>0 .AND. Ipar(7)>=Ipar(6) ) THEN
      Ipar(1) = -1
      stopbis = .TRUE.
   ENDIF
   IF ( stopbis ) RETURN
!
!     computes errors
!
   Fpar(5) = sqrt(distdot(N,R,1,R,1))
   Fpar(11) = Fpar(11) + 2*N
   IF ( Ipar(3)<0 ) THEN
!
!     compute the change in the solution vector
!
      Fpar(6) = Sx*sqrt(distdot(N,Delx,1,Delx,1))
      Fpar(11) = Fpar(11) + 2*N
      IF ( Ipar(7)<Mvpi+Mvpi+1 ) THEN
!
!     if this is the end of the first iteration, set fpar(3:4)
!
         Fpar(3) = Fpar(6)
         IF ( Ipar(3)==-1 ) Fpar(4) = Fpar(1)*Fpar(3) + Fpar(2)
      ENDIF
   ELSE
      Fpar(6) = Fpar(5)
   ENDIF
!
!     .. the test is struct this way so that when the value in fpar(6)
!       is not a valid number, STOPBIS is set to .true.
!
   IF ( Fpar(6)>Fpar(4) ) THEN
      stopbis = .FALSE.
      Ipar(11) = 0
   ELSE
      stopbis = .TRUE.
      Ipar(11) = 1
   ENDIF
!
END FUNCTION stopbis
!*==tidycg.f90 processed by SPAG 8.04RA 15:47 21 Feb 2024
!!SPAG Open source Personal, Educational or Academic User Clemson University  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
!-----end-of-stopbis
!-----------------------------------------------------------------------
SUBROUTINE tidycg(N,Ipar,Fpar,Sol,Delx)
   USE ISO_FORTRAN_ENV                 
   IMPLICIT NONE
!
! PARAMETER definitions rewritten by SPAG
!
   REAL(REAL64) , PARAMETER :: ZERO = 0.0D0
!
! Dummy argument declarations rewritten by SPAG
!
   INTEGER , INTENT(IN) :: N
   INTEGER , INTENT(INOUT) , DIMENSION(16) :: Ipar
   REAL(REAL64) , INTENT(INOUT) , DIMENSION(16) :: Fpar
   REAL(REAL64) , INTENT(INOUT) , DIMENSION(N) :: Sol
   REAL(REAL64) , INTENT(IN) , DIMENSION(N) :: Delx
!
! Local variable declarations rewritten by SPAG
!
   INTEGER :: i
!
! End of declarations rewritten by SPAG
!
!-----------------------------------------------------------------------
!     Some common operations required before terminating the CG routines
!-----------------------------------------------------------------------
!
   IF ( Ipar(12)/=0 ) THEN
      Ipar(1) = Ipar(12)
   ELSEIF ( Ipar(1)>0 ) THEN
      IF ( (Ipar(3)==999 .AND. Ipar(11)==1) .OR. Fpar(6)<=Fpar(4) ) THEN
         Ipar(1) = 0
      ELSEIF ( Ipar(7)>=Ipar(6) .AND. Ipar(6)>0 ) THEN
         Ipar(1) = -1
      ELSE
         Ipar(1) = -10
      ENDIF
   ENDIF
   IF ( Fpar(3)>ZERO .AND. Fpar(6)>ZERO .AND. Ipar(7)>Ipar(13) ) THEN
      Fpar(7) = log10(Fpar(3)/Fpar(6))/dble(Ipar(7)-Ipar(13))
   ELSE
      Fpar(7) = ZERO
   ENDIF
   DO i = 1 , N
      Sol(i) = Sol(i) + Delx(i)
   ENDDO
END SUBROUTINE tidycg
!*==brkdn.f90 processed by SPAG 8.04RA 15:47 21 Feb 2024
!!SPAG Open source Personal, Educational or Academic User Clemson University  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
!-----end-of-tidycg
!-----------------------------------------------------------------------
FUNCTION brkdn(Alpha,Ipar)
   USE ISO_FORTRAN_ENV                 
   IMPLICIT NONE
!
! PARAMETER definitions rewritten by SPAG
!
   REAL(REAL64) , PARAMETER :: ZERO = 0.0D0 , ONE = 1.0D0
!
! Function and Dummy argument declarations rewritten by SPAG
!
   LOGICAL :: brkdn
   REAL(REAL64) , INTENT(IN) :: Alpha
   INTEGER , INTENT(OUT) , DIMENSION(16) :: Ipar
!
! Local variable declarations rewritten by SPAG
!
   REAL(REAL64) :: beta
!
! End of declarations rewritten by SPAG
!
!-----------------------------------------------------------------------
!     test whether alpha is zero or an abnormal number, if yes,
!     this routine will return .true.
!
!     If alpha == 0, ipar(1) = -3,
!     if alpha is an abnormal number, ipar(1) = -9.
!-----------------------------------------------------------------------
   brkdn = .FALSE.
   IF ( Alpha>ZERO ) THEN
      beta = ONE/Alpha
      IF ( .NOT.beta>ZERO ) THEN
         brkdn = .TRUE.
         Ipar(1) = -9
      ENDIF
   ELSEIF ( Alpha<ZERO ) THEN
      beta = ONE/Alpha
      IF ( .NOT.beta<ZERO ) THEN
         brkdn = .TRUE.
         Ipar(1) = -9
      ENDIF
   ELSEIF ( Alpha==ZERO ) THEN
      brkdn = .TRUE.
      Ipar(1) = -3
   ELSE
      brkdn = .TRUE.
      Ipar(1) = -9
   ENDIF
END FUNCTION brkdn
!*==bisinit.f90 processed by SPAG 8.04RA 15:47 21 Feb 2024
!!SPAG Open source Personal, Educational or Academic User Clemson University  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
!-----end-of-brkdn
!-----------------------------------------------------------------------
SUBROUTINE bisinit(Ipar,Fpar,Wksize,Dsc,Lp,Rp,Wk)
   USE ISO_FORTRAN_ENV                 
   IMPLICIT NONE
!
! PARAMETER definitions rewritten by SPAG
!
   REAL(REAL64) , PARAMETER :: ZERO = 0.0D0 , ONE = 1.0D0
!
! Dummy argument declarations rewritten by SPAG
!
   INTEGER , INTENT(INOUT) , DIMENSION(16) :: Ipar
   REAL(REAL64) , INTENT(INOUT) , DIMENSION(16) :: Fpar
   INTEGER , INTENT(IN) :: Wksize
   INTEGER , INTENT(IN) :: Dsc
   LOGICAL , INTENT(OUT) :: Lp
   LOGICAL , INTENT(OUT) :: Rp
   REAL(REAL64) , INTENT(OUT) , DIMENSION(*) :: Wk
!
! Local variable declarations rewritten by SPAG
!
   INTEGER :: i
!
! End of declarations rewritten by SPAG
!
!-----------------------------------------------------------------------
!     some common initializations for the iterative solvers
!-----------------------------------------------------------------------
!
!     ipar(1) = -2 inidcate that there are not enough space in the work
!     array
!
   IF ( Ipar(4)<Wksize ) THEN
      Ipar(1) = -2
      Ipar(4) = Wksize
      RETURN
   ENDIF
!
   IF ( Ipar(2)>2 ) THEN
      Lp = .TRUE.
      Rp = .TRUE.
   ELSEIF ( Ipar(2)==2 ) THEN
      Lp = .FALSE.
      Rp = .TRUE.
   ELSEIF ( Ipar(2)==1 ) THEN
      Lp = .TRUE.
      Rp = .FALSE.
   ELSE
      Lp = .FALSE.
      Rp = .FALSE.
   ENDIF
   IF ( Ipar(3)==0 ) Ipar(3) = Dsc
!     .. clear the ipar elements used
   Ipar(7) = 0
   Ipar(8) = 0
   Ipar(9) = 0
   Ipar(10) = 0
   Ipar(11) = 0
   Ipar(12) = 0
   Ipar(13) = 0
!
!     fpar(1) must be between (0, 1), fpar(2) must be positive,
!     fpar(1) and fpar(2) can NOT both be zero
!     Normally returns ipar(1) = -4 to indicate any of above error
!
   IF ( Fpar(1)<ZERO .OR. Fpar(1)>=ONE .OR. Fpar(2)<ZERO .OR. (Fpar(1)==ZERO .AND. Fpar(2)==ZERO) ) THEN
      IF ( Ipar(1)==0 ) THEN
         Ipar(1) = -4
         RETURN
      ELSE
         Fpar(1) = 1.0D-6
         Fpar(2) = 1.0D-16
      ENDIF
   ENDIF
!     .. clear the fpar elements
   DO i = 3 , 10
      Fpar(i) = ZERO
   ENDDO
   IF ( Fpar(11)<ZERO ) Fpar(11) = ZERO
!     .. clear the used portion of the work array to zero
   DO i = 1 , Wksize
      Wk(i) = ZERO
   ENDDO
!
!-----end-of-bisinit
END SUBROUTINE bisinit
!*==mgsro.f90 processed by SPAG 8.04RA 15:47 21 Feb 2024
!!SPAG Open source Personal, Educational or Academic User Clemson University  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
!-----------------------------------------------------------------------
SUBROUTINE mgsro(Full,Lda,N,M,Ind,Ops,Vec,Hh,Ierr)
   USE ISO_FORTRAN_ENV                 
   IMPLICIT NONE
!
! PARAMETER definitions rewritten by SPAG
!
   REAL(REAL64) , PARAMETER :: ZERO = 0.0D0 , ONE = 1.0D0 , REORTH = 0.98D0
!
! Dummy argument declarations rewritten by SPAG
!
   INTEGER , INTENT(IN) :: Lda
   INTEGER , INTENT(IN) :: M
   LOGICAL , INTENT(IN) :: Full
   INTEGER :: N
   INTEGER , INTENT(IN) :: Ind
   REAL(REAL64) , INTENT(INOUT) :: Ops
   REAL(REAL64) , INTENT(INOUT) , DIMENSION(Lda,M) :: Vec
   REAL(REAL64) , INTENT(INOUT) , DIMENSION(M) :: Hh
   INTEGER , INTENT(OUT) :: Ierr
!
! Local variable declarations rewritten by SPAG
!
   REAL(REAL64) , EXTERNAL :: distdot
   REAL(REAL64) :: fct , nrm0 , nrm1 , thr
   INTEGER :: i , k
!
! End of declarations rewritten by SPAG
!
!-----------------------------------------------------------------------
!     MGSRO  -- Modified Gram-Schmidt procedure with Selective Re-
!               Orthogonalization
!     The ind'th vector of VEC is orthogonalized against the rest of
!     the vectors.
!
!     The test for performing re-orthogonalization is performed for
!     each indivadual vectors. If the cosine between the two vectors
!     is greater than 0.99 (REORTH = 0.99**2), re-orthogonalization is
!     performed. The norm of the 'new' vector is kept in variable NRM0,
!     and updated after operating with each vector.
!
!     full   -- .ture. if it is necessary to orthogonalize the ind'th
!               against all the vectors vec(:,1:ind-1), vec(:,ind+2:m)
!               .false. only orthogonalize againt vec(:,1:ind-1)
!     lda    -- the leading dimension of VEC
!     n      -- length of the vector in VEC
!     m      -- number of vectors can be stored in VEC
!     ind    -- index to the vector to be changed
!     ops    -- operation counts
!     vec    -- vector of LDA X M storing the vectors
!     hh     -- coefficient of the orthogonalization
!     ierr   -- error code
!               0 : successful return
!               -1: zero input vector
!               -2: input vector contains abnormal numbers
!               -3: input vector is a linear combination of others
!
!     External routines used: real*8 distdot
!-----------------------------------------------------------------------
!
!     compute the norm of the input vector
!
   nrm0 = distdot(N,Vec(1,Ind),1,Vec(1,Ind),1)
   Ops = Ops + N + N
   thr = nrm0*REORTH
   IF ( nrm0<=ZERO ) THEN
      Ierr = -1
      RETURN
   ELSEIF ( nrm0>ZERO .AND. ONE/nrm0>ZERO ) THEN
      Ierr = 0
   ELSE
      Ierr = -2
      RETURN
   ENDIF
!
!     Modified Gram-Schmidt loop
!
   IF ( Full ) THEN
      DO i = Ind + 1 , M
         fct = distdot(N,Vec(1,Ind),1,Vec(1,i),1)
         Hh(i) = fct
         DO k = 1 , N
            Vec(k,Ind) = Vec(k,Ind) - fct*Vec(k,i)
         ENDDO
         Ops = Ops + 4*N + 2
         IF ( fct*fct>thr ) THEN
            fct = distdot(N,Vec(1,Ind),1,Vec(1,i),1)
            Hh(i) = Hh(i) + fct
            DO k = 1 , N
               Vec(k,Ind) = Vec(k,Ind) - fct*Vec(k,i)
            ENDDO
            Ops = Ops + 4*N + 1
         ENDIF
         nrm0 = nrm0 - Hh(i)*Hh(i)
         IF ( nrm0<ZERO ) nrm0 = ZERO
         thr = nrm0*REORTH
      ENDDO
   ENDIF
!
   DO i = 1 , Ind - 1
      fct = distdot(N,Vec(1,Ind),1,Vec(1,i),1)
      Hh(i) = fct
      DO k = 1 , N
         Vec(k,Ind) = Vec(k,Ind) - fct*Vec(k,i)
      ENDDO
      Ops = Ops + 4*N + 2
      IF ( fct*fct>thr ) THEN
         fct = distdot(N,Vec(1,Ind),1,Vec(1,i),1)
         Hh(i) = Hh(i) + fct
         DO k = 1 , N
            Vec(k,Ind) = Vec(k,Ind) - fct*Vec(k,i)
         ENDDO
         Ops = Ops + 4*N + 1
      ENDIF
      nrm0 = nrm0 - Hh(i)*Hh(i)
      IF ( nrm0<ZERO ) nrm0 = ZERO
      thr = nrm0*REORTH
   ENDDO
!
!     test the resulting vector
!
   nrm1 = sqrt(distdot(N,Vec(1,Ind),1,Vec(1,Ind),1))
   Ops = Ops + N + N
   Hh(Ind) = nrm1
   IF ( nrm1<=ZERO ) THEN
      Ierr = -3
      RETURN
   ENDIF
!
!     scale the resulting vector
!
   fct = ONE/nrm1
   DO k = 1 , N
      Vec(k,Ind) = Vec(k,Ind)*fct
   ENDDO
   Ops = Ops + N + 1
!
!     normal return
!
   Ierr = 0
!     end surbotine mgsro
END SUBROUTINE mgsro
