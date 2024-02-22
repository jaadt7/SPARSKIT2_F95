!*==gen57pt.f90 processed by SPAG 8.04RA 15:47 21 Feb 2024
!!SPAG Open source Personal, Educational or Academic User Clemson University  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
!-----------------------------------------------------------------------
!------------------------end-of-csorted---------------------------------
!----------------------------------------------------------------------c
!                          S P A R S K I T                             c
!----------------------------------------------------------------------c
!    MATRIX GENERATION ROUTINES  -- FINITE DIFFERENCE MATRICES         c
!----------------------------------------------------------------------c
! contents:                                                            c
!----------                                                            c
! gen57pt  : generates 5-point and 7-point matrices.                   c
! gen57bl  : generates block 5-point and 7-point matrices.             c
!                                                                      c
! supporting routines:                                                 c
!---------                                                             c
! gensten  : generate the stencil (point version)                      c
! bsten    : generate the stencil (block version)                      c
! fdaddbc  : finite difference add boundary conditions                 c
! fdreduce : reduce the system to eliminate node with known values     c
! clrow    : clear a row of a CSR matrix                               c
! lctcsr   : locate the position of A(i,j) in CSR format               c
!----------------------------------------------------------------------c
SUBROUTINE gen57pt(Nx,Ny,Nz,Al,Mode,N,A,Ja,Ia,Iau,Rhs)
   USE ISO_FORTRAN_ENV                 
   IMPLICIT NONE
!
! PARAMETER definitions rewritten by SPAG
!
   REAL(REAL64) , PARAMETER :: ONE = 1.0D0
!
! Dummy argument declarations rewritten by SPAG
!
   INTEGER :: Nx
   INTEGER :: Ny
   INTEGER :: Nz
   REAL(REAL64) , DIMENSION(6) :: Al
   INTEGER :: Mode
   INTEGER :: N
   REAL(REAL64) , DIMENSION(*) :: A
   INTEGER , DIMENSION(*) :: Ja
   INTEGER , DIMENSION(*) :: Ia
   INTEGER , DIMENSION(*) :: Iau
   REAL(REAL64) , DIMENSION(*) :: Rhs
!
! Local variable declarations rewritten by SPAG
!
   LOGICAL :: genrhs , value
   REAL(REAL64) :: h , r
   INTEGER :: iedge , ix , iy , iz , kx , ky , kz , node
   REAL(REAL64) , DIMENSION(7) :: stencil
   EXTERNAL fdaddbc , fdreduce , getsten
!
! End of declarations rewritten by SPAG
!
!-----------------------------------------------------------------------
! On entry:
!
! nx      = number of grid points in x direction
! ny	  = number of grid points in y direction
! nz	  = number of grid points in z direction
! al      = array of size 6, carries the coefficient alpha of the
!           boundary conditions
! mode    = what to generate:
!           < 0 : generate the graph only,
!           = 0 : generate the matrix,
!           > 0 : generate the matrix and the right-hand side.
!
! On exit:
!
! n       = number of nodes with unknown values, ie number of rows
!           in the matrix
!
! a,ja,ia = resulting matrix in row-sparse format
!
! iau     = integer*n, containing the poisition of the diagonal element
!           in the a, ja, ia structure
!
! rhs     = the right-hand side
!
! External functions needed (must be supplied by caller)
!     afun, bfun, cfun, dfun, efun, ffun, gfun, hfun
!     betfun, gamfun
! They have the following prototype:
!     real*8 function xfun(x, y, z)
!     real*8 x, y, z
!-----------------------------------------------------------------------
! This subroutine computes the sparse matrix in compressed sparse row
! format for the elliptic equation:
!       d    du    d    du    d    du      du     du     du
! L u = --(A --) + --(B --) + --(C --) + D -- + E -- + F -- + G u = H u
!       dx   dx    dy   dy    dz   dz      dx     dy     dz
!
! with general Mixed Boundary conditions, on a rectangular 1-D,
! 2-D or 3-D grid using 2nd order centered difference schemes.
!
! The functions a, b, ..., g, h are known through the
! as afun, bfun, ..., gfun, hfun in this subroutine.
! NOTE: To obtain the correct matrix, any function that is not
! needed should be set to zero.  For example for two-dimensional
! problems, nz should be set to 1 and the functions cfun and ffun
! should be zero functions.
!
! The Boundary condition is specified in the following form:
!           du
!     alpha -- + beta u = gamma
!           dn
! Where alpha is constant at each side of the boundary surfaces.  Alpha
! is represented by parameter al.  It is expected to an array that
! contains enough elements to specify the boundaries for the problem,
! 1-D case needs two elements, 2-D needs 4 and 3-D needs 6.  The order
! of the boundaries in the array is left(west), right(east),
! bottom(south), top(north), front, rear.  Beta and gamma are functions
! of type real with three arguments x, y, z.  These two functions are
! known subroutine 'addbc' as betfun and gamfun.  They should following
! the same notion as afun ... hfun.  For more restriction on afun ...
! hfun, please read the documentation follows the subroutine 'getsten',
! and, for more on betfun and gamfun, please refer to the documentation
! under subroutine 'fdaddbc'.
!
! The nodes are ordered using natural ordering, first x direction, then
! y, then z.  The mesh size h is uniform and determined by grid points
! in the x-direction.
!
! The domain specified for the problem is [0 .ge. x .ge. 1],
! [0 .ge. y .ge. (ny-1)*h] and [0 .ge. z .ge. (nz-1)*h], where h is
! 1 / (nx-1).  Thus if non-Dirichlet boundary condition is specified,
! the mesh will have nx points along the x direction, ny along y and
! nz along z.  For 1-D case, both y and z value are assumed to zero
! when calling relavent functions that have three parameters.
! Similarly, for 2-D case, z is assumed to be zero.
!
! About the expectation of nx, ny and nz:
! nx is required to be .gt. 1 always;
! if the second dimension is present in the problem, then ny should be
! .gt. 1, else 1;
! if the third dimension is present in the problem, nz .gt. 1, else 1.
! when ny is 1, nz must be 1.
!-----------------------------------------------------------------------
!
!     stencil [1:7] has the following meaning:
!
!     center point = stencil(1)
!     west point = stencil(2)
!     east point = stencil(3)
!     south point = stencil(4)
!     north point = stencil(5)
!     front point = stencil(6)
!     back point = stencil(7)
!
!     al[1:6] carry the coefficient alpha in the similar order
!
!     west  side = al(1)
!     east  side = al(2)
!     south side = al(3)
!     north side = al(4)
!     front side = al(5)
!     back  side = al(6)
!
!                           al(4)
!                           st(5)
!                            |
!                            |
!                            |           al(6)
!                            |          .st(7)
!                            |     .
!         al(1)              | .             al(2)
!         st(2) ----------- st(1) ---------- st(3)
!                       .    |
!                   .        |
!               .            |
!            st(6)           |
!            al(5)           |
!                            |
!                           st(4)
!                           al(3)
!
!-------------------------------------------------------------------
!     some constants
!
!
!     local variables
!
!
!     nx has to be larger than 1
!
   IF ( Nx<=1 ) RETURN
   h = ONE/dble(Nx-1)
!
!     the mode
!
   value = (Mode>=0)
   genrhs = (Mode>0)
!
!     first generate the whole matrix as if the boundary condition does
!     not exist
!
   kx = 1
   ky = Nx
   kz = Nx*Ny
   iedge = 1
   node = 1
   DO iz = 1 , Nz
      DO iy = 1 , Ny
         DO ix = 1 , Nx
            Ia(node) = iedge
!
!     compute the stencil at the current node
!
            IF ( value ) CALL getsten(Ny,Nz,Mode,ix-1,iy-1,iz-1,stencil,h,r)
!     west
            IF ( ix>1 ) THEN
               Ja(iedge) = node - kx
               IF ( value ) A(iedge) = stencil(2)
               iedge = iedge + 1
            ENDIF
!     south
            IF ( iy>1 ) THEN
               Ja(iedge) = node - ky
               IF ( value ) A(iedge) = stencil(4)
               iedge = iedge + 1
            ENDIF
!     front plane
            IF ( iz>1 ) THEN
               Ja(iedge) = node - kz
               IF ( value ) A(iedge) = stencil(6)
               iedge = iedge + 1
            ENDIF
!     center node
            Ja(iedge) = node
            Iau(node) = iedge
            IF ( value ) A(iedge) = stencil(1)
            iedge = iedge + 1
!     east
            IF ( ix<Nx ) THEN
               Ja(iedge) = node + kx
               IF ( value ) A(iedge) = stencil(3)
               iedge = iedge + 1
            ENDIF
!     north
            IF ( iy<Ny ) THEN
               Ja(iedge) = node + ky
               IF ( value ) A(iedge) = stencil(5)
               iedge = iedge + 1
            ENDIF
!     back plane
            IF ( iz<Nz ) THEN
               Ja(iedge) = node + kz
               IF ( value ) A(iedge) = stencil(7)
               iedge = iedge + 1
            ENDIF
!     the right-hand side
            IF ( genrhs ) Rhs(node) = r
            node = node + 1
         ENDDO
      ENDDO
   ENDDO
   Ia(node) = iedge
!
!     Add in the boundary conditions
!
   CALL fdaddbc(Nx,Ny,Nz,A,Ja,Ia,Iau,Rhs,Al,h)
!
!     eliminate the boudary nodes from the matrix
!
   CALL fdreduce(Nx,Ny,Nz,Al,N,A,Ja,Ia,Iau,Rhs,stencil)
!
!     done
!
!-----end-of-gen57pt----------------------------------------------------
!-----------------------------------------------------------------------
END SUBROUTINE gen57pt
!*==getsten.f90 processed by SPAG 8.04RA 15:47 21 Feb 2024
!!SPAG Open source Personal, Educational or Academic User Clemson University  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
!-----------------------------------------------------------------------
SUBROUTINE getsten(Ny,Nz,Mode,Kx,Ky,Kz,Stencil,H,Rhs)
   USE ISO_FORTRAN_ENV                 
   IMPLICIT NONE
!
! PARAMETER definitions rewritten by SPAG
!
   REAL(REAL64) , PARAMETER :: ZERO = 0.0D0 , HALF = 0.5D0
!
! Dummy argument declarations rewritten by SPAG
!
   INTEGER , INTENT(IN) :: Ny
   INTEGER , INTENT(IN) :: Nz
   INTEGER , INTENT(IN) :: Mode
   INTEGER , INTENT(IN) :: Kx
   INTEGER , INTENT(IN) :: Ky
   INTEGER , INTENT(IN) :: Kz
   REAL(REAL64) , INTENT(INOUT) , DIMENSION(*) :: Stencil
   REAL(REAL64) , INTENT(IN) :: H
   REAL(REAL64) , INTENT(OUT) :: Rhs
!
! Local variable declarations rewritten by SPAG
!
   REAL(REAL64) , EXTERNAL :: afun , bfun , cfun , dfun , efun , ffun , gfun , hfun
   REAL(REAL64) :: cntr , coeff , hhalf , x , y , z
   INTEGER :: k
!
! End of declarations rewritten by SPAG
!
!-----------------------------------------------------------------------
!     This subroutine calculates the correct stencil values for
!     centered difference discretization of the elliptic operator
!     and the right-hand side
!
! L u = delx( A delx u ) + dely ( B dely u) + delz ( C delz u ) +
!	delx ( D u ) + dely (E u) + delz( F u ) + G u = H
!
!   For 2-D problems the discretization formula that is used is:
!
! h**2 * Lu == A(i+1/2,j)*{u(i+1,j) - u(i,j)} +
!	       A(i-1/2,j)*{u(i-1,j) - u(i,j)} +
!              B(i,j+1/2)*{u(i,j+1) - u(i,j)} +
!              B(i,j-1/2)*{u(i,j-1) - u(i,j)} +
!              (h/2)*D(i,j)*{u(i+1,j) - u(i-1,j)} +
!              (h/2)*E(i,j)*{u(i,j+1) - u(i,j-1)} +
!              (h/2)*E(i,j)*{u(i,j+1) - u(i,j-1)} +
!              (h**2)*G(i,j)*u(i,j)
!-----------------------------------------------------------------------
!     some constants
!
!
!     local variables
!
!
!     if mode < 0, we shouldn't have come here
!
   IF ( Mode<0 ) RETURN
!
   DO k = 1 , 7
      Stencil(k) = ZERO
   ENDDO
!
   hhalf = H*HALF
   x = H*dble(Kx)
   y = H*dble(Ky)
   z = H*dble(Kz)
   cntr = ZERO
!     differentiation wrt x:
   coeff = afun(x+hhalf,y,z)
   Stencil(3) = Stencil(3) + coeff
   cntr = cntr + coeff
!
   coeff = afun(x-hhalf,y,z)
   Stencil(2) = Stencil(2) + coeff
   cntr = cntr + coeff
!
   coeff = dfun(x,y,z)*hhalf
   Stencil(3) = Stencil(3) + coeff
   Stencil(2) = Stencil(2) - coeff
   IF ( Ny>1 ) THEN
!
!     differentiation wrt y:
!
      coeff = bfun(x,y+hhalf,z)
      Stencil(5) = Stencil(5) + coeff
      cntr = cntr + coeff
!
      coeff = bfun(x,y-hhalf,z)
      Stencil(4) = Stencil(4) + coeff
      cntr = cntr + coeff
!
      coeff = efun(x,y,z)*hhalf
      Stencil(5) = Stencil(5) + coeff
      Stencil(4) = Stencil(4) - coeff
      IF ( Nz>1 ) THEN
!
! differentiation wrt z:
!
         coeff = cfun(x,y,z+hhalf)
         Stencil(7) = Stencil(7) + coeff
         cntr = cntr + coeff
!
         coeff = cfun(x,y,z-hhalf)
         Stencil(6) = Stencil(6) + coeff
         cntr = cntr + coeff
!
         coeff = ffun(x,y,z)*hhalf
         Stencil(7) = Stencil(7) + coeff
         Stencil(6) = Stencil(6) - coeff
      ENDIF
   ENDIF
!
! contribution from function G:
!
   coeff = gfun(x,y,z)
   Stencil(1) = H*H*coeff - cntr
!
!     the right-hand side
!
   IF ( Mode>0 ) Rhs = H*H*hfun(x,y,z)
!
!------end-of-getsten---------------------------------------------------
!-----------------------------------------------------------------------
END SUBROUTINE getsten
!*==gen57bl.f90 processed by SPAG 8.04RA 15:47 21 Feb 2024
!!SPAG Open source Personal, Educational or Academic User Clemson University  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
!-----------------------------------------------------------------------
SUBROUTINE gen57bl(Nx,Ny,Nz,Nfree,Na,N,A,Ja,Ia,Iau,Stencil)
   USE ISO_FORTRAN_ENV                 
   IMPLICIT NONE
!
! PARAMETER definitions rewritten by SPAG
!
   REAL(REAL64) , PARAMETER :: ONE = 1.0D0
!
! Dummy argument declarations rewritten by SPAG
!
   INTEGER , INTENT(IN) :: Na
   INTEGER , INTENT(IN) :: Nx
   INTEGER :: Ny
   INTEGER :: Nz
   INTEGER :: Nfree
   INTEGER , INTENT(OUT) :: N
   REAL(REAL64) , INTENT(OUT) , DIMENSION(Na,1) :: A
   INTEGER , INTENT(OUT) , DIMENSION(*) :: Ja
   INTEGER , INTENT(OUT) , DIMENSION(*) :: Ia
   INTEGER , INTENT(OUT) , DIMENSION(*) :: Iau
   REAL(REAL64) , DIMENSION(7,1) :: Stencil
!
! Local variable declarations rewritten by SPAG
!
   REAL(REAL64) :: h
   INTEGER :: iedge , ix , iy , iz , k , kx , ky , kz , nfree2 , node
   EXTERNAL bsten
!
! End of declarations rewritten by SPAG
!
!     implicit real*8 (a-h,o-z)
!--------------------------------------------------------------------
! This subroutine computes the sparse matrix in compressed
! format for the elliptic operator
!
! L u = delx( a . delx u ) + dely ( b . dely u) + delz ( c . delz u ) +
!	delx ( d . u ) + dely (e . u) + delz( f . u ) + g . u
!
! Here u is a vector of nfree componebts and each of the functions
! a, b, c, d, e, f, g   is an (nfree x nfree) matrix depending of
! the coordinate (x,y,z).
! with Dirichlet Boundary conditions, on a rectangular 1-D,
! 2-D or 3-D grid using centered difference schemes.
!
! The functions a, b, ..., g are known through the
! subroutines  afunbl, bfunbl, ..., gfunbl. (user supplied) .
!
! uses natural ordering, first x direction, then y, then z
! mesh size h is uniform and determined by grid points
! in the x-direction.
!
! The output matrix is in Block -- Sparse Row format.
!
!--------------------------------------------------------------------
! parameters:
!-------------
! Input:
! ------
! nx      = number of points in x direction
! ny	  = number of points in y direction
! nz	  = number of points in z direction
! nfree   = number of degrees of freedom per point
! na	  = first dimension of array a as declared in calling
!           program. Must be .ge. nfree**2
!
! Output:
! ------
! n	  = dimension of matrix (output)
!
! a, ja, ia = resulting matrix in  Block Sparse Row format
!           a(1:nfree**2, j ) contains a nonzero block and ja(j)
!           contains the (block) column number of this block.
!           the block dimension of the matrix is n (output) and
!           therefore the total number of (scalar) rows is n x nfree.
!
! iau     = integer*n containing the position of the diagonal element
!           in the a, ja, ia structure
!
! Work space:
!------------
! stencil =  work array of size (7,nfree**2) [stores local stencils]
!
!--------------------------------------------------------------------
!
!     stencil (1:7,*) has the following meaning:
!
!     center point = stencil(1)
!     west point   = stencil(2)
!     east point   = stencil(3)
!     south point  = stencil(4)
!     north point  = stencil(5)
!     front point  = stencil(6)
!     back point   = stencil(7)
!
!
!                           st(5)
!                            |
!                            |
!                            |
!                            |          .st(7)
!                            |     .
!                            | .
!         st(2) ----------- st(1) ---------- st(3)
!                       .    |
!                   .        |
!               .            |
!            st(6)           |
!                            |
!                            |
!                           st(4)
!
!-------------------------------------------------------------------
!     some constants
!
!
!     local variables
!
!
   h = ONE/dble(Nx+1)
   kx = 1
   ky = Nx
   kz = Nx*Ny
   nfree2 = Nfree*Nfree
   iedge = 1
   node = 1
   DO iz = 1 , Nz
      DO iy = 1 , Ny
         DO ix = 1 , Nx
            Ia(node) = iedge
            CALL bsten(Ny,Nz,ix,iy,iz,Nfree,Stencil,h)
!     west
            IF ( ix>1 ) THEN
               Ja(iedge) = node - kx
               DO k = 1 , nfree2
                  A(k,iedge) = Stencil(2,k)
               ENDDO
               iedge = iedge + 1
            ENDIF
!     south
            IF ( iy>1 ) THEN
               Ja(iedge) = node - ky
               DO k = 1 , nfree2
                  A(k,iedge) = Stencil(4,k)
               ENDDO
               iedge = iedge + 1
            ENDIF
!     front plane
            IF ( iz>1 ) THEN
               Ja(iedge) = node - kz
               DO k = 1 , nfree2
                  A(k,iedge) = Stencil(6,k)
               ENDDO
               iedge = iedge + 1
            ENDIF
!     center node
            Ja(iedge) = node
            Iau(node) = iedge
            DO k = 1 , nfree2
               A(k,iedge) = Stencil(1,k)
            ENDDO
            iedge = iedge + 1
!     -- upper part
!     east
            IF ( ix<Nx ) THEN
               Ja(iedge) = node + kx
               DO k = 1 , nfree2
                  A(k,iedge) = Stencil(3,k)
               ENDDO
               iedge = iedge + 1
            ENDIF
!     north
            IF ( iy<Ny ) THEN
               Ja(iedge) = node + ky
               DO k = 1 , nfree2
                  A(k,iedge) = Stencil(5,k)
               ENDDO
               iedge = iedge + 1
            ENDIF
!     back plane
            IF ( iz<Nz ) THEN
               Ja(iedge) = node + kz
               DO k = 1 , nfree2
                  A(k,iedge) = Stencil(7,k)
               ENDDO
               iedge = iedge + 1
            ENDIF
!------next node -------------------------
            node = node + 1
         ENDDO
      ENDDO
   ENDDO
!
!     -- new version of BSR -- renumbering removed.
!     change numbering of nodes so that each ja(k) will contain the
!     actual column number in the original matrix of entry (1,1) of each
!     block (k).
!      do 101 k=1,iedge-1
!         ja(k) = (ja(k)-1)*nfree+1
! 101  continue
!
!      n = (node-1)*nfree
   N = node - 1
   Ia(node) = iedge
!--------------end-of-gen57bl-------------------------------------------
!-----------------------------------------------------------------------
END SUBROUTINE gen57bl
!*==bsten.f90 processed by SPAG 8.04RA 15:47 21 Feb 2024
!!SPAG Open source Personal, Educational or Academic User Clemson University  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
!-----------------------------------------------------------------------
SUBROUTINE bsten(Ny,Nz,Kx,Ky,Kz,Nfree,Stencil,H)
   USE ISO_FORTRAN_ENV                 
   IMPLICIT NONE
!
! PARAMETER definitions rewritten by SPAG
!
   REAL(REAL64) , PARAMETER :: ZERO = 0.0D0 , HALF = 0.5D0
!
! Dummy argument declarations rewritten by SPAG
!
   INTEGER , INTENT(IN) :: Ny
   INTEGER , INTENT(IN) :: Nz
   INTEGER , INTENT(IN) :: Kx
   INTEGER , INTENT(IN) :: Ky
   INTEGER , INTENT(IN) :: Kz
   INTEGER :: Nfree
   REAL(REAL64) , INTENT(INOUT) , DIMENSION(7,*) :: Stencil
   REAL(REAL64) , INTENT(IN) :: H
!
! Local variable declarations rewritten by SPAG
!
   REAL(REAL64) , DIMENSION(225) :: cntr , coeff
   REAL(REAL64) :: h2 , hhalf , x , y , z
   INTEGER :: i , k , nfree2
   EXTERNAL afunbl , bfunbl , cfunbl , dfunbl , efunbl , ffunbl , gfunbl
!
! End of declarations rewritten by SPAG
!
!-----------------------------------------------------------------------
!     This subroutine calcultes the correct block-stencil values for
!     centered difference discretization of the elliptic operator
!     (block version of stencil)
!
! L u = delx( a delx u ) + dely ( b dely u) + delz ( c delz u ) +
!       d delx ( u ) + e dely (u) + f delz( u ) + g u
!
!   For 2-D problems the discretization formula that is used is:
!
! h**2 * Lu == a(i+1/2,j)*{u(i+1,j) - u(i,j)} +
!              a(i-1/2,j)*{u(i-1,j) - u(i,j)} +
!              b(i,j+1/2)*{u(i,j+1) - u(i,j)} +
!              b(i,j-1/2)*{u(i,j-1) - u(i,j)} +
!              (h/2)*d(i,j)*{u(i+1,j) - u(i-1,j)} +
!              (h/2)*e(i,j)*{u(i,j+1) - u(i,j-1)} +
!              (h/2)*e(i,j)*{u(i,j+1) - u(i,j-1)} +
!              (h**2)*g(i,j)*u(i,j)
!-----------------------------------------------------------------------
!     some constants
!
!
!     local variables
!
!------------
   IF ( Nfree>15 ) THEN
      PRINT * , ' ERROR ** nfree too large '
      STOP
   ENDIF
!
   nfree2 = Nfree*Nfree
   DO k = 1 , nfree2
      cntr(k) = ZERO
      DO i = 1 , 7
         Stencil(i,k) = ZERO
      ENDDO
   ENDDO
!------------
   hhalf = H*HALF
   h2 = H*H
   x = H*dble(Kx)
   y = H*dble(Ky)
   z = H*dble(Kz)
! differentiation wrt x:
   CALL afunbl(Nfree,x+hhalf,y,z,coeff)
   DO k = 1 , nfree2
      Stencil(3,k) = Stencil(3,k) + coeff(k)
      cntr(k) = cntr(k) + coeff(k)
   ENDDO
!
   CALL afunbl(Nfree,x-hhalf,y,z,coeff)
   DO k = 1 , nfree2
      Stencil(2,k) = Stencil(2,k) + coeff(k)
      cntr(k) = cntr(k) + coeff(k)
   ENDDO
!
   CALL dfunbl(Nfree,x,y,z,coeff)
   DO k = 1 , nfree2
      Stencil(3,k) = Stencil(3,k) + coeff(k)*hhalf
      Stencil(2,k) = Stencil(2,k) - coeff(k)*hhalf
   ENDDO
   IF ( Ny>1 ) THEN
!
! differentiation wrt y:
!
      CALL bfunbl(Nfree,x,y+hhalf,z,coeff)
      DO k = 1 , nfree2
         Stencil(5,k) = Stencil(5,k) + coeff(k)
         cntr(k) = cntr(k) + coeff(k)
      ENDDO
!
      CALL bfunbl(Nfree,x,y-hhalf,z,coeff)
      DO k = 1 , nfree2
         Stencil(4,k) = Stencil(4,k) + coeff(k)
         cntr(k) = cntr(k) + coeff(k)
      ENDDO
!
      CALL efunbl(Nfree,x,y,z,coeff)
      DO k = 1 , nfree2
         Stencil(5,k) = Stencil(5,k) + coeff(k)*hhalf
         Stencil(4,k) = Stencil(4,k) - coeff(k)*hhalf
      ENDDO
      IF ( Nz>1 ) THEN
!
! differentiation wrt z:
!
         CALL cfunbl(Nfree,x,y,z+hhalf,coeff)
         DO k = 1 , nfree2
            Stencil(7,k) = Stencil(7,k) + coeff(k)
            cntr(k) = cntr(k) + coeff(k)
         ENDDO
!
         CALL cfunbl(Nfree,x,y,z-hhalf,coeff)
         DO k = 1 , nfree2
            Stencil(6,k) = Stencil(6,k) + coeff(k)
            cntr(k) = cntr(k) + coeff(k)
         ENDDO
!
         CALL ffunbl(Nfree,x,y,z,coeff)
         DO k = 1 , nfree2
            Stencil(7,k) = Stencil(7,k) + coeff(k)*hhalf
            Stencil(6,k) = Stencil(6,k) - coeff(k)*hhalf
         ENDDO
      ENDIF
   ENDIF
!
! discretization of  product by g:
!
   CALL gfunbl(Nfree,x,y,z,coeff)
   DO k = 1 , nfree2
      Stencil(1,k) = h2*coeff(k) - cntr(k)
   ENDDO
!
!------------end of bsten-----------------------------------------------
!-----------------------------------------------------------------------
END SUBROUTINE bsten
!*==fdreduce.f90 processed by SPAG 8.04RA 15:47 21 Feb 2024
!!SPAG Open source Personal, Educational or Academic User Clemson University  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
SUBROUTINE fdreduce(Nx,Ny,Nz,Alpha,N,A,Ja,Ia,Iau,Rhs,Stencil)
   USE ISO_FORTRAN_ENV                 
   IMPLICIT NONE
!
! PARAMETER definitions rewritten by SPAG
!
   REAL(REAL64) , PARAMETER :: ZERO = 0.0D0
!
! Dummy argument declarations rewritten by SPAG
!
   INTEGER , INTENT(IN) :: Nx
   INTEGER , INTENT(IN) :: Ny
   INTEGER , INTENT(IN) :: Nz
   REAL(REAL64) , INTENT(IN) , DIMENSION(*) :: Alpha
   INTEGER , INTENT(OUT) :: N
   REAL(REAL64) , INTENT(INOUT) , DIMENSION(*) :: A
   INTEGER , DIMENSION(*) :: Ja
   INTEGER , INTENT(INOUT) , DIMENSION(*) :: Ia
   INTEGER , INTENT(INOUT) , DIMENSION(*) :: Iau
   REAL(REAL64) , INTENT(INOUT) , DIMENSION(*) :: Rhs
   REAL(REAL64) , INTENT(INOUT) , DIMENSION(*) :: Stencil
!
! Local variable declarations rewritten by SPAG
!
   INTEGER :: i , iedge , j , k , kx , ky , kz , ld , lk , lx , ly , lz , nbnode , node , ux , uy , uz
   INTEGER , EXTERNAL :: lctcsr
   REAL(REAL64) :: val
!
! End of declarations rewritten by SPAG
!
!-----------------------------------------------------------------------
! This subroutine tries to reduce the size of the matrix by looking
! for Dirichlet boundary conditions at each surface and solve the boundary
! value and modify the right-hand side of related nodes, then clapse all
! the boundary nodes.
!-----------------------------------------------------------------------
!     parameters
!
!
!     local variables
!
!
!     The first half of this subroutine will try to change the right-hand
!     side of all the nodes that has a neighbor with Dirichlet boundary
!     condition, since in this case the value of the boundary point is
!     known.
!     Then in the second half, we will try to eliminate the boundary
!     points with known values (with Dirichlet boundary condition).
!
   kx = 1
   ky = Nx
   kz = Nx*Ny
   lx = 1
   ux = Nx
   ly = 1
   uy = Ny
   lz = 1
   uz = Nz
!
!     Here goes the first part. ----------------------------------------
!
!     the left (west) side
!
   IF ( Alpha(1)==ZERO ) THEN
      lx = 2
      DO k = 1 , Nz
         DO j = 1 , Ny
            node = (k-1)*kz + (j-1)*ky + 1
            nbnode = node + kx
            lk = lctcsr(nbnode,node,Ja,Ia)
            ld = Iau(node)
            val = Rhs(node)/A(ld)
!     modify the rhs
            Rhs(nbnode) = Rhs(nbnode) - A(lk)*val
         ENDDO
      ENDDO
   ENDIF
!
!     right (east) side
!
   IF ( Alpha(2)==ZERO ) THEN
      ux = Nx - 1
      DO k = 1 , Nz
         DO j = 1 , Ny
            node = (k-1)*kz + (j-1)*ky + Nx
            nbnode = node - kx
            lk = lctcsr(nbnode,node,Ja,Ia)
            ld = Iau(node)
            val = Rhs(node)/A(ld)
!     modify the rhs
            Rhs(nbnode) = Rhs(nbnode) - A(lk)*val
         ENDDO
      ENDDO
   ENDIF
!
!     if it's only 1-D, skip the following part
!
   IF ( Ny>1 ) THEN
!
!     the bottom (south) side
!
      IF ( Alpha(3)==ZERO ) THEN
         ly = 2
         DO k = 1 , Nz
            DO i = lx , ux
               node = (k-1)*kz + i
               nbnode = node + ky
               lk = lctcsr(nbnode,node,Ja,Ia)
               ld = Iau(node)
               val = Rhs(node)/A(ld)
!     modify the rhs
               Rhs(nbnode) = Rhs(nbnode) - A(lk)*val
            ENDDO
         ENDDO
      ENDIF
!
!     top (north) side
!
      IF ( Alpha(4)==ZERO ) THEN
         uy = Ny - 1
         DO k = 1 , Nz
            DO i = lx , ux
               node = (k-1)*kz + i + (Ny-1)*ky
               nbnode = node - ky
               lk = lctcsr(nbnode,node,Ja,Ia)
               ld = Iau(node)
               val = Rhs(node)/A(ld)
!     modify the rhs
               Rhs(nbnode) = Rhs(nbnode) - A(lk)*val
            ENDDO
         ENDDO
      ENDIF
!
!     if only 2-D skip the following section on z
!
      IF ( Nz>1 ) THEN
!
!     the front surface
!
         IF ( Alpha(5)==ZERO ) THEN
            lz = 2
            DO j = ly , uy
               DO i = lx , ux
                  node = (j-1)*ky + i
                  nbnode = node + kz
                  lk = lctcsr(nbnode,node,Ja,Ia)
                  ld = Iau(node)
                  val = Rhs(node)/A(ld)
!     modify the rhs
                  Rhs(nbnode) = Rhs(nbnode) - A(lk)*val
               ENDDO
            ENDDO
         ENDIF
!
!     rear surface
!
         IF ( Alpha(6)==ZERO ) THEN
            uz = Nz - 1
            DO j = ly , uy
               DO i = lx , ux
                  node = (Nz-1)*kz + (j-1)*ky + i
                  nbnode = node - kz
                  lk = lctcsr(nbnode,node,Ja,Ia)
                  ld = Iau(node)
                  val = Rhs(node)/A(ld)
!     modify the rhs
                  Rhs(nbnode) = Rhs(nbnode) - A(lk)*val
               ENDDO
            ENDDO
         ENDIF
      ENDIF
   ENDIF
!
!     now the second part ----------------------------------------------
!
!     go through all the actual nodes with unknown values, collect all
!     of them to form a new matrix in compressed sparse row format.
!
   kx = 1
   ky = ux - lx + 1
   kz = (uy-ly+1)*ky
   node = 1
   iedge = 1
   DO k = lz , uz
      DO j = ly , uy
         DO i = lx , ux
!
!     the corresponding old node number
            nbnode = ((k-1)*Ny+j-1)*Nx + i
!
!     copy the row into local stencil, copy is done is the exact
!     same order as the stencil is written into array a
            lk = Ia(nbnode)
            IF ( i>1 ) THEN
               Stencil(2) = A(lk)
               lk = lk + 1
            ENDIF
            IF ( j>1 ) THEN
               Stencil(4) = A(lk)
               lk = lk + 1
            ENDIF
            IF ( k>1 ) THEN
               Stencil(6) = A(lk)
               lk = lk + 1
            ENDIF
            Stencil(1) = A(lk)
            lk = lk + 1
            IF ( i<Nx ) THEN
               Stencil(3) = A(lk)
               lk = lk + 1
            ENDIF
            IF ( j<Ny ) THEN
               Stencil(5) = A(lk)
               lk = lk + 1
            ENDIF
            IF ( k<Nz ) Stencil(7) = A(lk)
!
!     first the ia pointer -- points to the beginning of each row
            Ia(node) = iedge
!
!     move the values from the local stencil to the new matrix
!
!     the neighbor on the left (west)
            IF ( i>lx ) THEN
               Ja(iedge) = node - kx
               A(iedge) = Stencil(2)
               iedge = iedge + 1
            ENDIF
!     the neighbor below (south)
            IF ( j>ly ) THEN
               Ja(iedge) = node - ky
               A(iedge) = Stencil(4)
               iedge = iedge + 1
            ENDIF
!     the neighbor in the front
            IF ( k>lz ) THEN
               Ja(iedge) = node - kz
               A(iedge) = Stencil(6)
               iedge = iedge + 1
            ENDIF
!     center node (itself)
            Ja(iedge) = node
            Iau(node) = iedge
            A(iedge) = Stencil(1)
            iedge = iedge + 1
!     the neighbor to the right (east)
            IF ( i<ux ) THEN
               Ja(iedge) = node + kx
               A(iedge) = Stencil(3)
               iedge = iedge + 1
            ENDIF
!     the neighbor above (north)
            IF ( j<uy ) THEN
               Ja(iedge) = node + ky
               A(iedge) = Stencil(5)
               iedge = iedge + 1
            ENDIF
!     the neighbor at the back
            IF ( k<uz ) THEN
               Ja(iedge) = node + kz
               A(iedge) = Stencil(7)
               iedge = iedge + 1
            ENDIF
!     the right-hand side
            Rhs(node) = Rhs(nbnode)
!------next node -------------------------
            node = node + 1
!
         ENDDO
      ENDDO
   ENDDO
!
   Ia(node) = iedge
!
!     the number of nodes in the final matrix is stored in n
!
   N = node - 1
!-----------------------------------------------------------------------
END SUBROUTINE fdreduce
!*==fdaddbc.f90 processed by SPAG 8.04RA 15:47 21 Feb 2024
!!SPAG Open source Personal, Educational or Academic User Clemson University  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
!-----end of fdreduce-----------------------------------------------------
!-----------------------------------------------------------------------
SUBROUTINE fdaddbc(Nx,Ny,Nz,A,Ja,Ia,Iau,Rhs,Al,H)
   USE ISO_FORTRAN_ENV                 
   IMPLICIT NONE
!
! PARAMETER definitions rewritten by SPAG
!
   REAL(REAL64) , PARAMETER :: HALF = 0.5D0 , ZERO = 0.0D0 , ONE = 1.0D0 , TWO = 2.0D0
!
! Dummy argument declarations rewritten by SPAG
!
   INTEGER , INTENT(IN) :: Nx
   INTEGER , INTENT(IN) :: Ny
   INTEGER , INTENT(IN) :: Nz
   REAL(REAL64) , INTENT(INOUT) , DIMENSION(7*Nx*Ny*Nz) :: A
   INTEGER , DIMENSION(7*Nx*Ny*Nz) :: Ja
   INTEGER , DIMENSION(Nx*Ny*Nz) :: Ia
   INTEGER , INTENT(IN) , DIMENSION(Nx*Ny*Nz) :: Iau
   REAL(REAL64) , INTENT(INOUT) , DIMENSION(Nx*Ny*Nz) :: Rhs
   REAL(REAL64) , INTENT(IN) , DIMENSION(6) :: Al
   REAL(REAL64) , INTENT(IN) :: H
!
! Local variable declarations rewritten by SPAG
!
   REAL(REAL64) , EXTERNAL :: afun , betfun , bfun , cfun , dfun , efun , ffun , gamfun
   REAL(REAL64) :: coeff , ctr , hhalf , x , y , z
   INTEGER :: i , j , k , kx , ky , kz , lx , ly , nbr , node , ux , uy
   INTEGER , EXTERNAL :: lctcsr
   CHARACTER(2) :: side
   EXTERNAL clrow
!
! End of declarations rewritten by SPAG
!
!-----------------------------------------------------------------------
! This subroutine will add the boundary condition to the linear system
! consutructed without considering the boundary conditions
!
! The Boundary condition is specified in the following form:
!           du
!     alpha -- + beta u = gamma
!           dn
! Alpha is stored in array AL.  The six side of the boundary appares
! in AL in the following order: left(west), right(east), bottom(south),
! top(north), front, back(rear). (see also the illustration in gen57pt)
! Beta and gamma appears as the functions, betfun and gamfun.
! They have the following prototype
!
! real*8 function xxxfun(x, y, z)
! real*8 x, y, z
!
! where x, y, z are vales in the range of [0, 1][0, (ny-1)*h]
! [0, (nz-1)*h]
!
! At the corners or boundary lines, the boundary conditions are applied
! in the follow order:
! 1) if one side is Dirichlet boundary condition, the Dirichlet boundary
!    condition is used;
! 2) if more than one sides are Dirichlet, the Direichlet condition
!    specified for X direction boundary will overwrite the one specified
!    for Y direction boundary which in turn has priority over Z
!     direction boundaries.
! 3) when all sides are non-Dirichlet, the average values are used.
!-----------------------------------------------------------------------
!     some constants
!
!
!     local variables
!
!         , gfun, hfun
!     external  gfun, hfun
!
   hhalf = HALF*H
   kx = 1
   ky = Nx
   kz = Nx*Ny
!
!     In 3-D case, we need to go through all 6 faces one by one. If
!     the actual dimension is lower, test on ny is performed first.
!     If ny is less or equals to 1, then the value of nz is not
!     checked.
!-----
!     the surface on the left (west) side
!     Concentrate on the contribution from the derivatives related to x,
!     The terms with derivative of x was assumed to be:
!
!     a(3/2,j,k)*[u(2,j,k)-u(1,j,k)] + a(1/2,j,k)*[u(0,j,k)-u(1,j,k)] +
!     h*d(1,j,k)*[u(2,j,k)-u(0,j,k)]/2
!
!     But they actually are:
!
!     2*{a(3/2,j,k)*[u(2,j,k)-u(1,j,k)] -
!     h*a(1,j,k)*[beta*u(1,j,k)-gamma]/alpha]} +
!     h*h*d(1,j,k)*[beta*u(1,j,k)-gamma]/alpha
!
!     Therefore, in terms of local stencil the right neighbor of a node
!     should be changed to 2*a(3/2,j,k),
!     The matrix never contains the left neighbor on this border, nothing
!     needs to be done about it.
!     The following terms should be added to the center stencil:
!     -a(3/2,j,k) + a(1/2,j,k) + [h*d(1,j,k)-2*a(1,j,k)]*h*beta/alpha
!
!     And these terms should be added to the corresponding right-hand side
!     [h*d(1,j,k)-2*a(1,j,k)]*h*gamma/alpha
!
!     Obviously, the formula do not apply for the Dirichlet Boundary
!     Condition, where alpha will be zero. In that case, we simply set
!     all the elements in the corresponding row to zero(0), then let
!     the diagonal element be beta, and the right-hand side be gamma.
!     Thus the value of u at that point will be set. Later on point
!     like this will be removed from the matrix, since they are of
!     know value before solving the system.(not done in this subroutine)
!
   x = ZERO
   side = 'x1'
   DO k = 1 , Nz
      z = (k-1)*H
      DO j = 1 , Ny
         y = (j-1)*H
         node = 1 + (j-1)*ky + (k-1)*kz
!
!     check to see if it's Dirichlet Boundary condition here
!
         IF ( Al(1)==ZERO ) THEN
            CALL clrow(node,A,Ia)
            A(Iau(node)) = betfun(side,x,y,z)
            Rhs(node) = gamfun(side,x,y,z)
         ELSE
!
!     compute the terms formulated above to modify the matrix.
!
!     the right neighbor is stroed in nbr'th posiiton in the a
            nbr = lctcsr(node,node+kx,Ja,Ia)
!
            coeff = TWO*afun(x,y,z)
            ctr = (H*dfun(x,y,z)-coeff)*H/Al(1)
            Rhs(node) = Rhs(node) + ctr*gamfun(side,x,y,z)
            ctr = afun(x-hhalf,y,z) + ctr*betfun(side,x,y,z)
            coeff = afun(x+hhalf,y,z)
            A(Iau(node)) = A(Iau(node)) - coeff + ctr
            A(nbr) = TWO*coeff
         ENDIF
      ENDDO
   ENDDO
!
!     the right (east) side boudary, similarly, the contirbution from
!     the terms containing the derivatives of x were assumed to be
!
!     a(nx+1/2,j,k)*[u(nx+1,j,k)-u(nx,j,k)] +
!     a(nx-1/2,j,k)*[u(nx-1,j,k)-u(nx,j,k)] +
!     d(nx,j,k)*[u(nx+1,j,k)-u(nx-1,j,k)]*h/2
!
!     Actualy they are:
!
!     2*{h*a(nx,j,k)*[gamma-beta*u(nx,j,k)]/alpha +
!     a(nx-1/2,j,k)*[u(nx-1,j,k)-u(nx,j,k)]} +
!     h*h*d(nx,j,k)*[gamma-beta*u(nx,j,k)]/alpha
!
!     The left stencil has to be set to 2*a(nx-1/2,j,k)
!
!     The following terms have to be added to the center stencil:
!
!     -a(nx-1/2,j,k)+a(nx+1/2,j,k)-[2*a(nx,j,k)+h*d(nx,j,k)]*beta/alpha
!
!     The following terms have to be added to the right-hand side:
!
!     -[2*a(nx,j,k)+h*d(nx,j,k)]*h*gamma/alpha
!
   x = ONE
   side = 'x2'
   DO k = 1 , Nz
      z = (k-1)*H
      DO j = 1 , Ny
         y = (j-1)*H
         node = (k-1)*kz + j*ky
!
         IF ( Al(2)==ZERO ) THEN
            CALL clrow(node,A,Ia)
            A(Iau(node)) = betfun(side,x,y,z)
            Rhs(node) = gamfun(side,x,y,z)
         ELSE
            nbr = lctcsr(node,node-kx,Ja,Ia)
!
            coeff = TWO*afun(x,y,z)
            ctr = (coeff+H*dfun(x,y,z))*H/Al(2)
            Rhs(node) = Rhs(node) - ctr*gamfun(side,x,y,z)
            ctr = afun(x+hhalf,y,z) - ctr*betfun(side,x,y,z)
            coeff = afun(x-hhalf,y,z)
            A(Iau(node)) = A(Iau(node)) - coeff + ctr
            A(nbr) = TWO*coeff
         ENDIF
      ENDDO
   ENDDO
!
!     If only one dimension, return now
!
   IF ( Ny<=1 ) RETURN
!
!     the bottom (south) side suface, This similar to the situation
!     with the left side, except all the function and realted variation
!     should be on the y.
!
!     These two block if statment here is to resolve the possible conflict
!     of assign the boundary value differently by different side of the
!     Dirichlet Boundary Conditions. They ensure that the edges that have
!     be assigned a specific value will not be reassigned.
!
   IF ( Al(1)==ZERO ) THEN
      lx = 2
   ELSE
      lx = 1
   ENDIF
   IF ( Al(2)==ZERO ) THEN
      ux = Nx - 1
   ELSE
      ux = Nx
   ENDIF
   y = ZERO
   side = 'y1'
   DO k = 1 , Nz
      z = (k-1)*H
      DO i = lx , ux
         x = (i-1)*H
         node = i + (k-1)*kz
!
         IF ( Al(3)==ZERO ) THEN
            CALL clrow(node,A,Ia)
            A(Iau(node)) = betfun(side,x,y,z)
            Rhs(node) = gamfun(side,x,y,z)
         ELSE
            nbr = lctcsr(node,node+ky,Ja,Ia)
!
            coeff = TWO*bfun(x,y,z)
            ctr = (H*efun(x,y,z)-coeff)*H/Al(3)
            Rhs(node) = Rhs(node) + ctr*gamfun(side,x,y,z)
            ctr = bfun(x,y-hhalf,z) + ctr*betfun(side,x,y,z)
            coeff = bfun(x,y+hhalf,z)
            A(Iau(node)) = A(Iau(node)) - coeff + ctr
            A(nbr) = TWO*coeff
         ENDIF
      ENDDO
   ENDDO
!
!     The top (north) side, similar to the right side
!
   y = (Ny-1)*H
   side = 'y2'
   DO k = 1 , Nz
      z = (k-1)*H
      DO i = lx , ux
         x = (i-1)*H
         node = (k-1)*kz + (Ny-1)*ky + i
!
         IF ( Al(4)==ZERO ) THEN
            CALL clrow(node,A,Ia)
            A(Iau(node)) = betfun(side,x,y,z)
            Rhs(node) = gamfun(side,x,y,z)
         ELSE
            nbr = lctcsr(node,node-ky,Ja,Ia)
!
            coeff = TWO*bfun(x,y,z)
            ctr = (coeff+H*efun(x,y,z))*H/Al(4)
            Rhs(node) = Rhs(node) - ctr*gamfun(side,x,y,z)
            ctr = bfun(x,y+hhalf,z) - ctr*betfun(side,x,y,z)
            coeff = bfun(x,y-hhalf,z)
            A(Iau(node)) = A(Iau(node)) - coeff + ctr
            A(nbr) = TWO*coeff
         ENDIF
      ENDDO
   ENDDO
!
!     If only has two dimesion to work on, return now
!
   IF ( Nz<=1 ) RETURN
!
!     The front side boundary
!
!     If the edges of the surface has been decided by Dirichlet Boundary
!     Condition, then leave them alone.
!
   IF ( Al(3)==ZERO ) THEN
      ly = 2
   ELSE
      ly = 1
   ENDIF
   IF ( Al(4)==ZERO ) THEN
      uy = Ny - 1
   ELSE
      uy = Ny
   ENDIF
!
   z = ZERO
   side = 'z1'
   DO j = ly , uy
      y = (j-1)*H
      DO i = lx , ux
         x = (i-1)*H
         node = i + (j-1)*ky
!
         IF ( Al(5)==ZERO ) THEN
            CALL clrow(node,A,Ia)
            A(Iau(node)) = betfun(side,x,y,z)
            Rhs(node) = gamfun(side,x,y,z)
         ELSE
            nbr = lctcsr(node,node+kz,Ja,Ia)
!
            coeff = TWO*cfun(x,y,z)
            ctr = (H*ffun(x,y,z)-coeff)*H/Al(5)
            Rhs(node) = Rhs(node) + ctr*gamfun(side,x,y,z)
            ctr = cfun(x,y,z-hhalf) + ctr*betfun(side,x,y,z)
            coeff = cfun(x,y,z+hhalf)
            A(Iau(node)) = A(Iau(node)) - coeff + ctr
            A(nbr) = TWO*coeff
         ENDIF
      ENDDO
   ENDDO
!
!     Similiarly for the top side of the boundary suface
!
   z = (Nz-1)*H
   side = 'z2'
   DO j = ly , uy
      y = (j-1)*H
      DO i = lx , ux
         x = (i-1)*H
         node = (Nz-1)*kz + (j-1)*ky + i
!
         IF ( Al(6)==ZERO ) THEN
            CALL clrow(node,A,Ia)
            A(Iau(node)) = betfun(side,x,y,z)
            Rhs(node) = gamfun(side,x,y,z)
         ELSE
            nbr = lctcsr(node,node-kz,Ja,Ia)
!
            coeff = TWO*cfun(x,y,z)
            ctr = (coeff+H*ffun(x,y,z))*H/Al(6)
            Rhs(node) = Rhs(node) - ctr*gamfun(side,x,y,z)
            ctr = cfun(x,y,z+hhalf) - ctr*betfun(side,x,y,z)
            coeff = cfun(x,y,z-hhalf)
            A(Iau(node)) = A(Iau(node)) - coeff + ctr
            A(nbr) = TWO*coeff
         ENDIF
      ENDDO
   ENDDO
!
!     all set
!
!-----------------------------------------------------------------------
END SUBROUTINE fdaddbc
!*==clrow.f90 processed by SPAG 8.04RA 15:47 21 Feb 2024
!!SPAG Open source Personal, Educational or Academic User Clemson University  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
!-----end of fdaddbc----------------------------------------------------
!-----------------------------------------------------------------------
SUBROUTINE clrow(I,A,Ia)
   USE ISO_FORTRAN_ENV                 
   IMPLICIT NONE
!
! Dummy argument declarations rewritten by SPAG
!
   INTEGER , INTENT(IN) :: I
   REAL(REAL64) , INTENT(OUT) , DIMENSION(*) :: A
   INTEGER , INTENT(IN) , DIMENSION(*) :: Ia
!
! Local variable declarations rewritten by SPAG
!
   INTEGER :: k
!
! End of declarations rewritten by SPAG
!
!-----------------------------------------------------------------------
!     clear the row i to all zero, but still keep the structure of the
!     CSR matrix
!-----------------------------------------------------------------------
   DO k = Ia(I) , Ia(I+1) - 1
      A(k) = 0.0D0
   ENDDO
!
!-----end of clrow------------------------------------------------------
END SUBROUTINE clrow
!*==lctcsr.f90 processed by SPAG 8.04RA 15:47 21 Feb 2024
!!SPAG Open source Personal, Educational or Academic User Clemson University  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
!-----------------------------------------------------------------------
FUNCTION lctcsr(I,J,Ja,Ia)
   IMPLICIT NONE
!
! Function and Dummy argument declarations rewritten by SPAG
!
   INTEGER :: lctcsr
   INTEGER , INTENT(IN) :: I
   INTEGER , INTENT(IN) :: J
   INTEGER , INTENT(IN) , DIMENSION(*) :: Ja
   INTEGER , INTENT(IN) , DIMENSION(*) :: Ia
!
! Local variable declarations rewritten by SPAG
!
   INTEGER :: k
!
! End of declarations rewritten by SPAG
!
!-----------------------------------------------------------------------
!     locate the position of a matrix element in a CSR format
!     returns -1 if the desired element is zero
!-----------------------------------------------------------------------
   lctcsr = -1
   k = Ia(I)
   DO WHILE ( k<Ia(I+1) .AND. (lctcsr==-1) )
      IF ( Ja(k)==J ) lctcsr = k
      k = k + 1
   ENDDO
!
!-----------------------------------------------------------------------
END FUNCTION lctcsr
