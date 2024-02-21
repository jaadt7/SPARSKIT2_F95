!*==lstif3.f90 processed by SPAG 8.04RA 15:47 21 Feb 2024
!!SPAG Open source Personal, Educational or Academic User Clemson University  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
!-----------------------------------------------------------------------
SUBROUTINE lstif3(Ske,Fe,Xe,Ye,xyk,funb,func,fung)
!---------------------------------------------------------------------------
!
!  This subroutine computes the local stiffness matrix for the
!   Diffusion-Convection Equation with the
!    variable cofficients, 'K(x,y), B(x,y), C(x,y) '
!
!  -Div( K(x,y) T(x,y)) + B(x,y) Tx + C(x,y) Ty = G
!
!   Here K(x,y) is a 2x2 Matrix, where each entry is a function of x and y.
!
!   K, B, C and G need to be supplied by user.
!    They need to be defined as externals in the calling routines.
!
!   PSI(i,x,y) : i-th shape fucntions on the standard triangle N, i=1, 2, 3
!      where N is the following.
!
!             (-1,1)
!                .
!                . .
!		 .   .
!		 .       .
!		 . . . . . . (1,-1)
!	       (-1,-1)
!
!   Local stiffness matrix is obtained by integral on the current
!   element. To do so, change the current coordinates to N
!    by Affine mapping, sending
!
!    (xe(1),ye(1))   ---> (-1,-1)
!    (xe(2),ye(2))   ---> (1,-1)
!    (xe(3),ye(3))   ---> (-1,1) .
!
!    Then we perform the integration on N
!     by Gaussian Quadrature with 9 points.
!
!---------------------------------------------------------------------------
!
!  on entry
!  ---------
!
!   xe       = x coordinates of the nodes in the current element.
!   ye       = y coordinates of the nodes in the current element.
!   xyk      = subroutine defining the function K(x,y).
!   funb     = function defining the function b(x,y).
!   func     = function defining the function c(x,y).
!   fung     = function defining the function g(x,y).
!
!---------------------------------------------------------------------------
!
!  on return
!  ---------
!
!   ske : Local Stiffness Matrix.( 3x3 in this subroutine.)
!   fe  : Local Load Vector.
!
!---------------------------------------------------------------------------
   USE ISO_FORTRAN_ENV                 
   IMPLICIT NONE
!
! Dummy argument declarations rewritten by SPAG
!
   REAL(REAL64) , INTENT(OUT) , DIMENSION(3,3) :: Ske
   REAL(REAL64) , INTENT(OUT) , DIMENSION(3) :: Fe
   REAL(REAL64) , INTENT(IN) , DIMENSION(3) :: Xe
   REAL(REAL64) , INTENT(IN) , DIMENSION(3) :: Ye
   EXTERNAL xyk
   REAL(REAL64) , EXTERNAL :: funb
   REAL(REAL64) , EXTERNAL :: func
   REAL(REAL64) , EXTERNAL :: fung
!
! Local variable declarations rewritten by SPAG
!
   REAL(REAL64) :: a11 , a12 , a21 , a22 , b1 , b2 , derv2 , r , rj1 , rj2 , rja , s , t1 , t2 , t3 , t4 , u11 , u12 , u21 , u22 , &
                 & v1 , v2 , w , x , x1 , x2 , x3 , y , y1 , y2 , y3
   REAL(REAL64) , DIMENSION(3,2) :: dn
   REAL(REAL64) , DIMENSION(9) , SAVE :: gau1 , gau2 , wei
   INTEGER :: i , j , k , npt
   REAL(REAL64) , EXTERNAL :: psi
   REAL(REAL64) , DIMENSION(2,2) :: xyke
!
! End of declarations rewritten by SPAG
!
 
!  Gau1 and Gau2 are the Gaussian Quadrature Points for the Traingle N,
!   and Wei, are the corresponding weights.
!
!   They are derived from the 1-D case by Reiterated integrals.
!
   DATA gau1/ - 0.8 , -0.1127016654 , 0.5745966692 , -0.8872983346 , -0.5 , -0.1127016654 , -0.9745966692 , -0.8872983346 , -0.8/
   DATA gau2/3* - 0.7745966692 , 3*0. , 3*0.7745966692/
   DATA wei/0.2738575107 , 0.4381720172 , 0.2738551072 , 0.2469135803 , 0.3950617284 , 0.2469135803 , 0.03478446464 ,              &
      & 0.05565514341 , 0.03478446464/
 
   npt = 9
 
!
!    Compute the Affine mappings from the current triangle to the
!     standard triangle N. Integration will be performed on that
!     triangle by Gaussian quadrature.
!
!    T = A X + B
!
!    A11, A12, A21, A22, B1, B2 will denote the entries of
!     A & B.
!
   x1 = Xe(1)
   x2 = Xe(2)
   x3 = Xe(3)
   y1 = Ye(1)
   y2 = Ye(2)
   y3 = Ye(3)
 
   rj1 = (x3-x1)*(y2-y3) - (x2-x3)*(y3-y1)
   rj2 = (x3-x1)*(y1-y2) - (x1-x2)*(y3-y1)
   a11 = 2*(y1-y3)/rj1
   a12 = 2*(x3-x1)/rj1
   a21 = 2*(y1-y2)/rj2
   a22 = 2*(x2-x1)/rj2
   b1 = 1. - a11*x2 - a12*y2
   b2 = -1. - a21*x2 - a22*y2
!
!  Compute the first order partial derivatives of the shape functions.
!   dn(i,1) and dn(i,2) are the first order partial derivativ of i-th shape function
!    with respect to x and y, respectively.
!
   dn(1,1) = -0.5*(a11+a21)
   dn(1,2) = -0.5*(a12+a22)
   dn(2,1) = 0.5*a11
   dn(2,2) = 0.5*a12
   dn(3,1) = 0.5*a21
   dn(3,2) = 0.5*a22
!  Compute the Jacobian associated with T.
   rja = a11*a22 - a12*a21
!
!  Find the inverse mapping of T
!
 
   u11 = a22/rja
   u12 = -a12/rja
   u21 = -a21/rja
   u22 = a11/rja
   v1 = -u11*b1 - u12*b2
   v2 = -u21*b1 - u22*b2
 
   DO i = 1 , 3
      t4 = 0.
      DO j = 1 , 3
         t1 = 0.
         t2 = 0.
         t3 = 0.
         DO k = 1 , npt
            r = gau1(k)
            s = gau2(k)
            w = wei(k)
 
            x = u11*r + u12*s + v1
            y = u21*r + u22*s + v2
 
            CALL xyk(xyke,x,y)
 
            derv2 = dn(i,1)*dn(j,1)*xyke(1,1) + dn(i,2)*dn(j,2)*xyke(2,2) + dn(i,1)*dn(j,2)*xyke(1,2) + dn(i,2)*dn(j,1)*xyke(2,1)
            IF ( j==1 ) t4 = t4 + w*fung(x,y)*psi(i,r,s)
 
            t1 = t1 + w*derv2
            t2 = t2 + w*funb(x,y)*psi(i,r,s)
            t3 = t3 + w*func(x,y)*psi(i,r,s)
         ENDDO
 
         Ske(i,j) = (t1+t2*dn(j,1)+t3*dn(j,2))/rja
      ENDDO
      Fe(i) = t4/rja
   ENDDO
 
END SUBROUTINE lstif3
!*==refall.f90 processed by SPAG 8.04RA 15:47 21 Feb 2024
!!SPAG Open source Personal, Educational or Academic User Clemson University  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
!--- end of lstif3 -----------------------------------------------------
!-----------------------------------------------------------------------
SUBROUTINE refall(Nx,Nelx,Ijk,Node,Ndeg,X,Y,Ichild,Iparnts,Nodcode,Nxmax,Nelmax,Ierr)
   USE ISO_FORTRAN_ENV                 
   IMPLICIT NONE
!
! Dummy argument declarations rewritten by SPAG
!
   INTEGER , INTENT(INOUT) :: Nx
   INTEGER , INTENT(IN) :: Node
   INTEGER , INTENT(IN) :: Ndeg
   INTEGER , INTENT(INOUT) :: Nelx
   INTEGER , INTENT(INOUT) , DIMENSION(Node,*) :: Ijk
   REAL(REAL64) , INTENT(INOUT) , DIMENSION(*) :: X
   REAL(REAL64) , INTENT(INOUT) , DIMENSION(*) :: Y
   INTEGER , INTENT(INOUT) , DIMENSION(Ndeg,1) :: Ichild
   INTEGER , INTENT(INOUT) , DIMENSION(2,Nx) :: Iparnts
   INTEGER , INTENT(INOUT) , DIMENSION(Nx) :: Nodcode
   INTEGER , INTENT(IN) :: Nxmax
   INTEGER , INTENT(IN) :: Nelmax
   INTEGER , INTENT(OUT) :: Ierr
!
! Local variable declarations rewritten by SPAG
!
   INTEGER :: i , ii , ipar1 , ipar2 , j , jchild , jj , jnod , k , k1 , k2 , last , nel , nelxnew , nxnew
   INTEGER , DIMENSION(20) :: inod , midnode
!
! End of declarations rewritten by SPAG
!
   INTEGER :: spag_nextblock_1
   spag_nextblock_1 = 1
   SPAG_DispatchLoop_1: DO
      SELECT CASE (spag_nextblock_1)
      CASE (1)
!-------------------------------------------------------------
! refines a finite element grid using triangular elements.
! uses mid points to refine all the elements of the grid.
!
! nx	= number of nodes at input
! nelx	= number of elements at input
! ijk	= connectivity matrix: for node k, ijk(*,k) point to the
!         nodes of element k.
! node  = first dimension of array ijk [should be >=3]
! ndeg	= first dimension of array ichild which is at least as large
!         as the max degree of each node
! x,y   = real*8 arrays containing the x(*) and y(*) coordinates
!	  resp. of the nodes.
! ichild= list of the children of a node: ichild(1,k) stores
!         the position in ichild(*,k)  of the last child so far.
!         (local use)
! iparnts= list of the 2 parents of each node.
!         (local use)
! nodcode= boundary information list for each node with the
!	   following meaning:
!	nodcode(i) = 0 -->  node i is internal
!	nodcode(i) = 1 -->  node i is a boundary but not a corner point
!	nodcode(i) = 2 -->  node i is a corner point.
! corner elements are used only to generate the grid by refinement
! since they do not  correspond to real elements.
! nxmax  = maximum number of nodes allowed. If during the algorithm
!          the number of nodes being created exceeds nxmax then
!	   refall  quits without modifying the (x,y) xoordinates
!	   and nx, nelx. ijk is modified. Also ierr is set to 1.
! nelmax = same as above for number of elements allowed. See ierr..
! ierr	 = error message:
!	   0 --> normal return
!	   1 --> refall quit because nxmax  was exceeded.
!	   2 --> refall quit because nelmax was exceeded.
!--------------------------------------------------------------
!---------------------------------------------------------------
! inilitialize lists of children and parents --
! data structure is as follows
! ichild(1,k) stores the position of last child of node k so far in list
! ichild(j,k) , j .ge. 2 = list of children of node k.
! iparnts(1,k) and iparnts(2,k) are the two parents of node k.
!---------------------------------------------------------------
!------ do a first check :
         IF ( Nx<Nxmax ) THEN
            IF ( Nelx>=Nelmax ) THEN
               spag_nextblock_1 = 3
               CYCLE SPAG_DispatchLoop_1
            ENDIF
!------ initialize
            DO k = 1 , Nx
               DO j = 2 , Ndeg
                  Ichild(j,k) = 0
               ENDDO
               Ichild(1,k) = 1
               Iparnts(1,k) = 0
               Iparnts(2,k) = 0
            ENDDO
!------- initialize nelxnew and nxnew
            nelxnew = Nelx
            nxnew = Nx
            Ierr = 0
!--------------------------------------------------------------
! main loop: scan all elements
!--------------------------------------------------------------
!     do 100 nel = nelx,1,-1
            DO nel = 1 , Nelx
! note : interesting question which order is best for parallelism?
! alternative order: do 100 nel = nelx, 1, -1
!
!------ unpack nodes of element
               DO i = 1 , Node
                  inod(i) = Ijk(i,nel)
! convention: node after last node = first node.
                  inod(Node+i) = inod(i)
                  midnode(i) = 0
               ENDDO
!--------------------------------------------------------------
! for each new potential node determine if it has already been
! numbered. a potential node is the middle of any two nodes ..
!--------------------------------------------------------------
               SPAG_Loop_3_1: DO ii = 1 , Node
                  k1 = inod(ii)
                  k2 = inod(ii+1)
!------- test for current pair :
                  last = Ichild(1,k1)
                  DO k = 2 , last
                     jchild = Ichild(k,k1)
                     ipar1 = Iparnts(1,jchild)
                     ipar2 = Iparnts(2,jchild)
                     IF ( ((ipar1==k1 .AND. ipar2==k2)) .OR. ((ipar2==k1 .AND. ipar1==k2)) ) THEN
! node has already been created and numbered ....
                        midnode(ii) = jchild
!... therefore it must be an internal node
                        Nodcode(jchild) = 0
!... and no new node to create.
                        CYCLE SPAG_Loop_3_1
                     ENDIF
!-----------------------------------------------------
                  ENDDO
!
! else  create a new node
!
                  nxnew = nxnew + 1
                  IF ( nxnew>Nxmax ) THEN
                     spag_nextblock_1 = 2
                     CYCLE SPAG_DispatchLoop_1
                  ENDIF
!-------
                  X(nxnew) = (X(k1)+X(k2))*0.5
                  Y(nxnew) = (Y(k1)+Y(k2))*0.5
                  midnode(ii) = nxnew
!
! update nodcode information -- normally min0(nodcode(k1),nodcode(k2))
!
                  Nodcode(nxnew) = min0(1,Nodcode(k1),Nodcode(k2))
!
! update parents and children's lists
!
                  Iparnts(1,nxnew) = k1
                  Iparnts(2,nxnew) = k2
!
                  last = last + 1
                  Ichild(last,k1) = nxnew
                  Ichild(1,k1) = last
!
                  last = Ichild(1,k2) + 1
                  Ichild(last,k2) = nxnew
                  Ichild(1,k2) = last
!
               ENDDO SPAG_Loop_3_1
!
!------- replace current element by new one
!
               DO i = 1 , Node
                  jnod = midnode(i)
                  Ijk(i,nel) = jnod
               ENDDO
!-------create new elements
               DO ii = 1 , Node
                  nelxnew = nelxnew + 1
                  IF ( nelxnew>Nelmax ) THEN
                     spag_nextblock_1 = 3
                     CYCLE SPAG_DispatchLoop_1
                  ENDIF
                  Ijk(1,nelxnew) = inod(ii)
                  k = ii
                  DO jj = 2 , Node
                     Ijk(jj,nelxnew) = midnode(k)
                     k = k + 2
                     IF ( k>Node ) k = k - Node
                  ENDDO
               ENDDO
!------ done !
            ENDDO
            Nx = nxnew
            Nelx = nelxnew
            RETURN
         ENDIF
         spag_nextblock_1 = 2
      CASE (2)
         Ierr = 1
         RETURN
      CASE (3)
         Ierr = 2
         EXIT SPAG_DispatchLoop_1
      END SELECT
   ENDDO SPAG_DispatchLoop_1
END SUBROUTINE refall
!*==checkref.f90 processed by SPAG 8.04RA 15:47 21 Feb 2024
!!SPAG Open source Personal, Educational or Academic User Clemson University  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
!
SUBROUTINE checkref(Nx,Nelx,Nodcode,Nbound,Nxnew,Nelxnew)
   IMPLICIT NONE
!
! Dummy argument declarations rewritten by SPAG
!
   INTEGER , INTENT(IN) :: Nx
   INTEGER , INTENT(IN) :: Nelx
   INTEGER , INTENT(IN) , DIMENSION(Nx) :: Nodcode
   INTEGER , INTENT(INOUT) :: Nbound
   INTEGER , INTENT(OUT) :: Nxnew
   INTEGER , INTENT(OUT) :: Nelxnew
!
! Local variable declarations rewritten by SPAG
!
   INTEGER :: j
!
! End of declarations rewritten by SPAG
!
!-------------------------------------------------------------
! returns the expected the new number of nodes and
! elemnts of refall is applied to current grid once.
!
! nx	= number of nodes at input
! nelx	= number of elements at input
! [ijk	= connectivity matrix: for node k, ijk(*,k) point to the
!         nodes of element k.] UNUSED here.
! nbound  = number of boundary points on entry - enter zero if
!	     unknown
!
! nodcode= boundary information list for each node with the
!	   following meaning:
!	nodcode(i) = 0 -->  node i is internal
!	nodcode(i) = 1 -->  node i is a boundary but not a corner point
!	nodcode(i) = 2 -->  node i is a corner point.
!
! nxnew  = new number of nodes if refall were to be applied
! nelxnew = same for nelx.
!--------------------------------------------------------------
!
   Nelxnew = Nelx*4
!
! count the number of boundary nodes
!
   IF ( Nbound==0 ) THEN
      DO j = 1 , Nx
         IF ( Nodcode(j)>=1 ) Nbound = Nbound + 1
      ENDDO
   ENDIF
! number of edges=[3*(number of elmts) + number of bound nodes ]/ 2
   Nxnew = Nx + (3*Nelx+Nbound)/2
   Nbound = 2*Nbound
END SUBROUTINE checkref
!*==unassbl.f90 processed by SPAG 8.04RA 15:47 21 Feb 2024
!!SPAG Open source Personal, Educational or Academic User Clemson University  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
!-----------------------------------------------------------------------
SUBROUTINE unassbl(A,Na,F,Nx,Nelx,Ijk,Node,X,Y,Ierr,xyk)
!-----------------------------------------------------------------------
! a      = un-assembled matrix on output
! na	 = 1-st dimension of a.  a(na,node,node)
!
! f      = right hand side (global load vector) in un-assembled form
! nx     = number of nodes at input
! nelx	 = number of elements at input
! ijk	 = connectivity matrix: for node k, ijk(*,k) point to the
!          nodes of element k.
! node	 = total number of nodal points in each element
!	   also second dimension of a.
!
! x,y   = real*8 arrays containing the $x$ and $y$ coordinates
!	  resp. of the nodes.
!         K11, K22, and K12 at that element.
! ierr	= error message integer .
!	  ierr = 0 --> normal return
!	  ierr = 1 --> negative area encountered (due to bad
!	           numbering of nodes of an element)
!
! xyk	= subroutine defining the material properties at each
!         element. Form:
! 	call xyk(nel,xyke,x,y,ijk,node)
!--------------------------------------------------------------
   USE ISO_FORTRAN_ENV                 
   IMPLICIT NONE
!
! Dummy argument declarations rewritten by SPAG
!
   INTEGER , INTENT(IN) :: Na
   INTEGER , INTENT(IN) :: Node
   REAL(REAL64) , INTENT(OUT) , DIMENSION(Na,Node,Node) :: A
   REAL(REAL64) , INTENT(OUT) , DIMENSION(Node,1) :: F
   INTEGER , INTENT(IN) :: Nx
   INTEGER , INTENT(IN) :: Nelx
   INTEGER , INTENT(IN) , DIMENSION(Node,1) :: Ijk
   REAL(REAL64) , DIMENSION(1) :: X
   REAL(REAL64) , DIMENSION(1) :: Y
   INTEGER , INTENT(INOUT) :: Ierr
   EXTERNAL xyk
!
! Local variable declarations rewritten by SPAG
!
   REAL(REAL64) :: det
   REAL(REAL64) , DIMENSION(3) :: fe , xe , ye
   INTEGER :: i , j , ka , kb , nel
   REAL(REAL64) , DIMENSION(3,3) :: ske
   REAL(REAL64) , DIMENSION(2,2) :: xyke
   EXTERNAL estif3
!
! End of declarations rewritten by SPAG
!
 
!--------------------------------------------------------------
!   initialize
!--------------------------------------------------------------
   DO i = 1 , Node
      DO j = 1 , Nx
         F(i,j) = 0.0D0
      ENDDO
   ENDDO
!---------------------------------------------------
! main loop
!---------------------------------------------------
   DO nel = 1 , Nelx
!
! get coordinetes of nodal points
!
      DO i = 1 , Node
         j = Ijk(i,nel)
         xe(i) = X(j)
         ye(i) = Y(j)
      ENDDO
!
! compute determinant
!
      det = xe(2)*(ye(3)-ye(1)) + xe(3)*(ye(1)-ye(2)) + xe(1)*(ye(2)-ye(3))
      IF ( det<=0. ) THEN
         PRINT * , 'nel' , nel , ' det = ' , det
         PRINT * , xe(1) , xe(2) , xe(3)
         PRINT * , ye(1) , ye(2) , ye(3)
      ENDIF
!
! set material properties
!
      CALL xyk(xyke,X,Y)
!
! construct element stiffness matrix
!
      Ierr = 0
      CALL estif3(nel,ske,det,xe,ye,xyke,Ierr)
      IF ( Ierr/=0 ) THEN
         WRITE (*,*) 'ERROR: estif3 gave an error' , Ierr
         RETURN
      ENDIF
!        write (8,'(9f8.4)') ((ske(i,j),j=1,3),i=1,3)
! assemble: add element stiffness matrix to global matrix
!
      DO ka = 1 , Node
         F(ka,nel) = fe(ka)
         DO kb = 1 , Node
            A(nel,ka,kb) = ske(ka,kb)
         ENDDO
      ENDDO
   ENDDO
END SUBROUTINE unassbl
!*==unassbl_lstif.f90 processed by SPAG 8.04RA 15:47 21 Feb 2024
!!SPAG Open source Personal, Educational or Academic User Clemson University  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
!-----------------------------------------------------------------------
SUBROUTINE unassbl_lstif(A,Na,F,Nx,Nelx,Ijk,Node,X,Y,Ierr,Xyk,Funb,Func,Fung)
!-----------------------------------------------------------------------
! a      = un-assembled matrix on output
!
! na	 = 1-st dimension of a.  a(na,node,node)
!
! f      = right hand side (global load vector) in un-assembled form
!
! nx     = number of nodes at input
!
! nelx	 = number of elements at input
!
! ijk	 = connectivity matrix: for node k, ijk(*,k) point to the
!          nodes of element k.
!
! node	 = total number of nodal points in each element
!	   also second dimension of a.
!
! x,y   = real*8 arrays containing the $x$ and $y$ coordinates
!	  resp. of the nodes.
!         K11, K22, and K12 at that element.
!
! ierr	= error message integer .
!	  ierr = 0 --> normal return
!	  ierr = 1 --> negative area encountered (due to bad
!	           numbering of nodes of an element)
!
! xyk	= subroutine defining the material properties at each
!         element. Form:  	call xyk(xyke,x,y)
!
! funb, = functions needed for the definition of lstif3 problem
! func,
! fung
!--------------------------------------------------------------
! moulitsa@cs.umn.edu : It uses lstif3 problem
!--------------------------------------------------------------
   USE ISO_FORTRAN_ENV                 
   IMPLICIT NONE
!
! Dummy argument declarations rewritten by SPAG
!
   INTEGER , INTENT(IN) :: Na
   INTEGER , INTENT(IN) :: Node
   REAL(REAL64) , INTENT(OUT) , DIMENSION(Na,Node,Node) :: A
   REAL(REAL64) , INTENT(OUT) , DIMENSION(Node,*) :: F
   INTEGER , INTENT(IN) :: Nx
   INTEGER , INTENT(IN) :: Nelx
   INTEGER , INTENT(IN) , DIMENSION(Node,*) :: Ijk
   REAL(REAL64) , INTENT(IN) , DIMENSION(*) :: X
   REAL(REAL64) , INTENT(IN) , DIMENSION(*) :: Y
   INTEGER , INTENT(OUT) :: Ierr
   REAL(REAL64) , EXTERNAL :: Xyk
   REAL(REAL64) , EXTERNAL :: Funb
   REAL(REAL64) , EXTERNAL :: Func
   REAL(REAL64) , EXTERNAL :: Fung
!
! Local variable declarations rewritten by SPAG
!
   REAL(REAL64) , DIMENSION(3) :: fe , xe , ye
   INTEGER :: i , j , ka , kb , nel
   REAL(REAL64) , DIMENSION(3,3) :: ske
   EXTERNAL lstif3
!
! End of declarations rewritten by SPAG
!
!--------------------------------------------------------------
!   initialize
!--------------------------------------------------------------
   DO i = 1 , Node
      DO j = 1 , Nx
         F(i,j) = 0.0D0
      ENDDO
   ENDDO
 
!---------------------------------------------------
! main loop
!---------------------------------------------------
   DO nel = 1 , Nelx
!
! get coordinetes of nodal points
!
      DO i = 1 , Node
         j = Ijk(i,nel)
         xe(i) = X(j)
         ye(i) = Y(j)
      ENDDO
!
! compute determinant
!
!	det=xe(2)*(ye(3)-ye(1))+xe(3)*(ye(1)-ye(2))+xe(1)*(ye(2)-ye(3))
!       if ( det .le. 0.) then
!         print *, 'nel', nel, ' det = ' , det
!         print *, xe(1), xe(2), xe(3)
!         print *, ye(1), ye(2), ye(3)
!       end if
!
! construct element stiffness matrix
!
      Ierr = 0
      CALL lstif3(ske,fe,xe,ye,Xyk,Funb,Func,Fung)
!        write (8,'(9f8.4)') ((ske(i,j),j=1,3),i=1,3)
!
! assemble: add element stiffness matrix to global matrix
!
      DO ka = 1 , Node
         F(ka,nel) = fe(ka)
         DO kb = 1 , Node
            A(nel,ka,kb) = ske(ka,kb)
         ENDDO
      ENDDO
 
   ENDDO
 
END SUBROUTINE unassbl_lstif
!*==assmbo.f90 processed by SPAG 8.04RA 15:47 21 Feb 2024
!!SPAG Open source Personal, Educational or Academic User Clemson University  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
!-----------------------------------------------------------------------
SUBROUTINE assmbo(Nx,Nelx,Node,Ijk,Nodcode,X,Y,A,Ja,Ia,F,Iwk,Jwk,Ierr,Xyk,Funb,Func,Fung)
!-----------------------------------------------------------------------
! nx     = number of nodes at input
!
! nelx	 = number of elements at input
!
! node	 = total number of nodal points in each element
!
! ijk	 = connectivity matrix: for node k, ijk(*,k) point to the
!          nodes of element k.
!
! nodcode= boundary information list for each node with the
!	   following meaning:
!	nodcode(i) = 0 -->  node i is internal
!	nodcode(i) = 1 -->  node i is a boundary but not a corner point
!	nodcode(i) = 2 -->  node i is a corner point (corner points
!
! x,y   = real arrays containing the $x$ and $y$ coordinates
!	  resp. of the nodes.
!
! a,ja,ia= assembled matrix on output
!
! f      = right hand side (global load vector)
!
! iwk,jwk = two integer work arrays.
!
! ierr	= error message integer .
!	  ierr = 0 --> normal return
!	  ierr = 1 --> negative area encountered (due to bad
!	           numbering of nodes of an element)
!
! xyk	= subroutine defining the material properties at each
!         element. Form:
! 	call xyk(nel,xyke,x,y,ijk,node) with on return
!         xyke =  material constant matrices.
!         for each element nel, xyke(1,nel),xyke(2,nel)
!         and xyke(3,nel) represent the constants
!         K11, K22, and K12 at that element.
!--------------------------------------------------------------
! moulitsa@cs.umn.edu : It has been modified so as to handle
!      more types of domains/meshes i.e.  |\ /|
!                                         | X |
!                                         |/ \|
!--------------------------------------------------------------
   USE ISO_FORTRAN_ENV                 
   IMPLICIT NONE
!
! Dummy argument declarations rewritten by SPAG
!
   INTEGER , INTENT(IN) :: Node
   INTEGER , INTENT(IN) :: Nx
   INTEGER , INTENT(IN) :: Nelx
   INTEGER , INTENT(IN) , DIMENSION(Node,1) :: Ijk
   INTEGER , INTENT(IN) , DIMENSION(1) :: Nodcode
   REAL(REAL64) , INTENT(IN) , DIMENSION(1) :: X
   REAL(REAL64) , INTENT(IN) , DIMENSION(1) :: Y
   REAL(REAL64) , INTENT(INOUT) , DIMENSION(*) :: A
   INTEGER , INTENT(INOUT) , DIMENSION(*) :: Ja
   INTEGER , INTENT(INOUT) , DIMENSION(1) :: Ia
   REAL(REAL64) , INTENT(INOUT) , DIMENSION(1) :: F
   INTEGER , INTENT(INOUT) , DIMENSION(1) :: Iwk
   INTEGER , INTENT(INOUT) , DIMENSION(1) :: Jwk
   INTEGER , INTENT(INOUT) :: Ierr
   REAL(REAL64) , EXTERNAL :: Xyk
   REAL(REAL64) , EXTERNAL :: Funb
   REAL(REAL64) , EXTERNAL :: Func
   REAL(REAL64) , EXTERNAL :: Fung
!
! Local variable declarations rewritten by SPAG
!
   REAL(REAL64) , DIMENSION(3) :: fe , xe , ye
   INTEGER :: i , ii , ilast , irowst , ista , isto , j , jj , k , ka , kb , knod , ksav , ksavn , nel
   REAL(REAL64) , DIMENSION(3,3) :: ske
   EXTERNAL lstif3
!
! End of declarations rewritten by SPAG
!
 
!--------------------------------------------------------------
!   initialize
!--------------------------------------------------------------
   DO i = 1 , Nx
      F(i) = 0.0
   ENDDO
! initialize  pointer arrays.
   DO k = 1 , Nx + 1
      Ia(k) = 1
      Jwk(k) = 0
   ENDDO
   DO k = 1 , Nelx
      DO j = 1 , Node
         knod = Ijk(j,k)
         Ia(knod) = Ia(knod) + 2
      ENDDO
   ENDDO
!---------------------------------------------------
   DO k = 1 , Nx
      IF ( Nodcode(k)>=1 ) Ia(k) = Ia(k) + 1
   ENDDO
!
   ksav = Ia(1)
   Ia(1) = 1
   DO j = 2 , Nx + 1
      ksavn = Ia(j)
      Ia(j) = Ia(j-1) + ksav
      Iwk(j-1) = Ia(j-1) - 1
      ksav = ksavn
   ENDDO
 
!-----------------
! main loop
!-----------------
   DO nel = 1 , Nelx
!
! get coordinates of nodal points
!
      DO i = 1 , Node
         j = Ijk(i,nel)
         xe(i) = X(j)
         ye(i) = Y(j)
      ENDDO
!
! compute determinant
!
!         det=xe(2)*(ye(3)-ye(1))+xe(3)*(ye(1)-ye(2))+xe(1)*(ye(2)-ye(3))
!
! set material properties
!
!        call xyk(nel,xyke,x,y,ijk,node)
!
! construct element stiffness matrix
!
      Ierr = 0
!
!        call evalg(nel, fe, xe, ye, fung, ierr)
!        call estif3(nel,ske,det,xe,ye,xyke,ierr)
      CALL lstif3(ske,fe,xe,ye,Xyk,Funb,Func,Fung)
      IF ( Ierr/=0 ) RETURN
!
! assemble: add element stiffness matrix to global matrix
!
      DO ka = 1 , Node
         ii = Ijk(ka,nel)
         F(ii) = F(ii) + fe(ka)
!
! unpack row into jwk1
!
         irowst = Ia(ii)
         ilast = Iwk(ii)
         DO k = irowst , ilast
            Jwk(Ja(k)) = k
         ENDDO
!
         DO kb = 1 , Node
!
! column number = jj
!
            jj = Ijk(kb,nel)
            k = Jwk(jj)
            IF ( k==0 ) THEN
               ilast = ilast + 1
               Jwk(jj) = ilast
               Ja(ilast) = jj
               A(ilast) = ske(ka,kb)
            ELSE
               A(k) = A(k) + ske(ka,kb)
            ENDIF
         ENDDO
! refresh jwk
         DO k = irowst , ilast
            Jwk(Ja(k)) = 0
         ENDDO
         Iwk(ii) = ilast
      ENDDO
!
   ENDDO
 
! squeeze away the zero entries
! added so as to handle more type of domains/meshes
   DO i = 1 , Nx
      ista = Ia(i)
      isto = Ia(i+1) - 1
      SPAG_Loop_2_1: DO j = ista , isto
         IF ( Ja(j)==0 ) THEN
            Iwk(i) = j - ista
            EXIT SPAG_Loop_2_1
         ENDIF
      ENDDO SPAG_Loop_2_1
   ENDDO
 
   DO i = 2 , Nx
      ksav = Ia(i)
      Ia(i) = Ia(i-1) + Iwk(i-1)
      ksavn = Ia(i)
      DO j = 0 , Iwk(i) - 1
         Ja(ksavn+j) = Ja(ksav+j)
         A(ksavn+j) = A(ksav+j)
      ENDDO
   ENDDO
   Ia(Nx+1) = Ia(Nx) + Iwk(Nx)
 
END SUBROUTINE assmbo
!*==assmbo2.f90 processed by SPAG 8.04RA 15:47 21 Feb 2024
!!SPAG Open source Personal, Educational or Academic User Clemson University  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
!-----------------------------------------------------------------------
SUBROUTINE assmbo2(Nx,Nelx,Node,Ijk,X,Y,A,Ja,Ia,F,Iwk,Jwk,Ierr,Xyk,Funb,Func,Fung)
!-----------------------------------------------------------------------
! nx     = number of nodes at input
!
! nelx	 = number of elements at input
!
! node	 = total number of nodal points in each element
!
! ijk	 = connectivity matrix: for node k, ijk(*,k) point to the
!          nodes of element k.
!
! x,y   = real arrays containing the $x$ and $y$ coordinates
!	  resp. of the nodes.
!
! a,ja,ia= assembled matrix on output
!
! f      = right hand side (global load vector)
!
! iwk,jwk = two integer work arrays.
!
! ierr	= error message integer .
!	  ierr = 0 --> normal return
!	  ierr = 1 --> negative area encountered (due to bad
!	           numbering of nodes of an element)
!
! xyk	= subroutine defining the material properties at each
!         element. Form:
! 	call xyk(nel,xyke,x,y,ijk,node) with on return
!         xyke =  material constant matrices.
!         for each element nel, xyke(1,nel),xyke(2,nel)
!         and xyke(3,nel) represent the constants
!         K11, K22, and K12 at that element.
!--------------------------------------------------------------
!
! moulitsa@cs.umn.edu : This routine yields the same results
!   as assmbo. It differs in that it constructs the ia array
!   by creating a list with the adjacent nodes for each node
!
!--------------------------------------------------------------
   USE ISO_FORTRAN_ENV                 
   IMPLICIT NONE
!
! Dummy argument declarations rewritten by SPAG
!
   INTEGER , INTENT(IN) :: Node
   INTEGER , INTENT(IN) :: Nx
   INTEGER , INTENT(IN) :: Nelx
   INTEGER , INTENT(IN) , DIMENSION(Node,*) :: Ijk
   REAL(REAL64) , INTENT(IN) , DIMENSION(*) :: X
   REAL(REAL64) , INTENT(IN) , DIMENSION(*) :: Y
   REAL(REAL64) , INTENT(INOUT) , DIMENSION(*) :: A
   INTEGER , INTENT(INOUT) , DIMENSION(*) :: Ja
   INTEGER , INTENT(INOUT) , DIMENSION(*) :: Ia
   REAL(REAL64) , INTENT(INOUT) , DIMENSION(*) :: F
   INTEGER , INTENT(INOUT) , DIMENSION(*) :: Iwk
   INTEGER , INTENT(INOUT) , DIMENSION(*) :: Jwk
   INTEGER , INTENT(INOUT) :: Ierr
   REAL(REAL64) , EXTERNAL :: Xyk
   REAL(REAL64) , EXTERNAL :: Funb
   REAL(REAL64) , EXTERNAL :: Func
   REAL(REAL64) , EXTERNAL :: Fung
!
! Local variable declarations rewritten by SPAG
!
   REAL(REAL64) , DIMENSION(3) :: fe , xe , ye
   INTEGER :: i , ii , ilast , irowst , j , jj , k , ka , kb , knod , ksav , ksavn , nedges , nel
   INTEGER , DIMENSION(500) :: kwk
   REAL(REAL64) , DIMENSION(3,3) :: ske
   EXTERNAL lstif3
!
! End of declarations rewritten by SPAG
!
!--------------------------------------------------------------
!   initialize
!--------------------------------------------------------------
   DO i = 1 , Nx
      F(i) = 0.0
      Iwk(i) = 0
      kwk(i) = 0
   ENDDO
 
! iwk : how many elements a node belongs to
   DO k = 1 , Nelx
      DO j = 1 , Node
         knod = Ijk(j,k)
         Iwk(knod) = Iwk(knod) + 1
      ENDDO
   ENDDO
!
! iwk : prepare for csr like format
   ksav = Iwk(1)
   Iwk(1) = 1
   DO j = 2 , Nx + 1
      ksavn = Iwk(j)
      Iwk(j) = Iwk(j-1) + ksav
      ksav = ksavn
   ENDDO
!
! jwk : list of elements a node belongs to
   k = 1
   DO i = 1 , Nelx
      DO j = 1 , Node
         knod = Ijk(j,i)
         k = Iwk(knod)
         Jwk(k) = i
         Iwk(knod) = Iwk(knod) + 1
      ENDDO
   ENDDO
 
! iwk : transform iwk back to what it was
   DO i = Nx + 1 , 2 , -1
      Iwk(i) = Iwk(i-1)
   ENDDO
   Iwk(1) = 1
 
! kwk : mark edges that a node is associated with
   nedges = 1
   Ia(1) = 1
   DO i = 1 , Nx
      kwk(i) = i
      DO j = Iwk(i) , Iwk(i+1) - 1
         DO k = 1 , Node
            knod = Ijk(k,Jwk(j))
            IF ( kwk(knod)/=i ) THEN
               kwk(knod) = i
               nedges = nedges + 1
            ENDIF
         ENDDO
      ENDDO
      Ia(i+1) = nedges
   ENDDO
   DO i = 2 , Nx + 1
      Ia(i) = Ia(i) + i - 1
      Iwk(i-1) = Ia(i-1) - 1
      Jwk(i) = 0
   ENDDO
   Jwk(1) = 0
 
!-----------------
! main loop
!-----------------
   DO nel = 1 , Nelx
!
! get coordinates of nodal points
!
      DO i = 1 , Node
         j = Ijk(i,nel)
         xe(i) = X(j)
         ye(i) = Y(j)
      ENDDO
!
! compute determinant
!
!         det=xe(2)*(ye(3)-ye(1))+xe(3)*(ye(1)-ye(2))+xe(1)*(ye(2)-ye(3))
!
! set material properties
!
!        call xyk(nel,xyke,x,y,ijk,node)
!
! construct element stiffness matrix
!
      Ierr = 0
!
!        call evalg(nel, fe, xe, ye, fung, ierr)
!        call estif3(nel,ske,fe,det,xe,ye,xyke,ierr)
      CALL lstif3(ske,fe,xe,ye,Xyk,Funb,Func,Fung)
      IF ( Ierr/=0 ) RETURN
!
! assemble: add element stiffness matrix to global matrix
!
      DO ka = 1 , Node
         ii = Ijk(ka,nel)
         F(ii) = F(ii) + fe(ka)
!
! unpack row into jwk1
!
         irowst = Ia(ii)
         ilast = Iwk(ii)
         DO k = irowst , ilast
            Jwk(Ja(k)) = k
         ENDDO
!
         DO kb = 1 , Node
!
! column number = jj
!
            jj = Ijk(kb,nel)
            k = Jwk(jj)
            IF ( k==0 ) THEN
               ilast = ilast + 1
               Jwk(jj) = ilast
               Ja(ilast) = jj
               A(ilast) = ske(ka,kb)
            ELSE
               A(k) = A(k) + ske(ka,kb)
            ENDIF
         ENDDO
! refresh jwk
         DO k = irowst , ilast
            Jwk(Ja(k)) = 0
         ENDDO
         Iwk(ii) = ilast
      ENDDO
!
   ENDDO
 
END SUBROUTINE assmbo2
!*==chkelmt.f90 processed by SPAG 8.04RA 15:47 21 Feb 2024
!!SPAG Open source Personal, Educational or Academic User Clemson University  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
!-----------------------------------------------------------------------
SUBROUTINE chkelmt(X,Y,Nelx,Ijk,Node)
   USE ISO_FORTRAN_ENV                 
   IMPLICIT NONE
!
! Dummy argument declarations rewritten by SPAG
!
   INTEGER , INTENT(IN) :: Node
   REAL(REAL64) , INTENT(IN) , DIMENSION(*) :: X
   REAL(REAL64) , INTENT(IN) , DIMENSION(*) :: Y
   INTEGER , INTENT(IN) :: Nelx
   INTEGER , INTENT(INOUT) , DIMENSION(Node,*) :: Ijk
!
! Local variable declarations rewritten by SPAG
!
   REAL(REAL64) :: det
   INTEGER :: j , nel
!
! End of declarations rewritten by SPAG
!
!-----------------------------------------------------------------------
! this subsourine checks the labeling within each elment and reorders
! the nodes in they ar not correctly ordered.
!-----------------------------------------------------------------------
   DO nel = 1 , Nelx
      det = X(Ijk(2,nel))*(Y(Ijk(3,nel))-Y(Ijk(1,nel))) + X(Ijk(3,nel))*(Y(Ijk(1,nel))-Y(Ijk(2,nel))) + X(Ijk(1,nel))              &
          & *(Y(Ijk(2,nel))-Y(Ijk(3,nel)))
!
! if determinant negative exchange last two nodes of elements.
!
      IF ( det<0.0D0 ) THEN
         j = Ijk(2,nel)
         Ijk(2,nel) = Ijk(3,nel)
         Ijk(3,nel) = j
      ENDIF
   ENDDO
!
END SUBROUTINE chkelmt
!*==dlauny.f90 processed by SPAG 8.04RA 15:47 21 Feb 2024
!!SPAG Open source Personal, Educational or Academic User Clemson University  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
!-----------------------------------------------------------------------
SUBROUTINE dlauny(X,Y,Nodes,Elmnts,Nemax,Nelmnt)
   USE ISO_FORTRAN_ENV                 
   IMPLICIT NONE
!
! Dummy argument declarations rewritten by SPAG
!
   INTEGER , INTENT(IN) :: Nodes
   INTEGER , INTENT(IN) :: Nemax
   REAL(REAL64) , INTENT(INOUT) , DIMENSION(Nodes) :: X
   REAL(REAL64) , INTENT(INOUT) , DIMENSION(Nodes) :: Y
   INTEGER , INTENT(INOUT) , DIMENSION(Nemax,3) :: Elmnts
   INTEGER , INTENT(INOUT) :: Nelmnt
!
! Local variable declarations rewritten by SPAG
!
   REAL(REAL64) :: cx , cy , dx , dy , r2 , rn2 , x2 , x3 , xl , xmax , xmin , xr , y2 , y3 , yl , ymax , ymin , yr , z
   INTEGER :: i , i1 , i2 , i3 , ie , in , j , je , k , l , match , nart , ndel , newel , nn
!
! End of declarations rewritten by SPAG
!
!
! code written by P.K. Sweby
! simple delauney triangulation routine (non optimal)
!
!     ******************************************************************
!     *                                                                *
!     * Performs a Delaunay triangularisation of a region given a set  *
!     * of mesh points.                                                *
!     *   X,Y    :- 1D arrays holding coordinates of mesh points.      *
!     *             dimensioned AT LEAST NODES+3.                      *
!     *   NODES  :- number of mesh points.                             *
!     *   ELMNTS :- INTEGER array, dimensioned NEMAX x 3, which on exit*
!     *             contains the index of global nodes associated with *
!     *             each element.                                      *
!     *   NELMNT :- on exit contains the number of elements in the     *
!     *             triangularisation.                                 *
!     *                                                                *
!     *                                   P.K.Sweby                    *
!     *                                                                *
!     ******************************************************************
!
   nn = 0
!
!      pi=4.0*atan(1.0)
!
!     Calculate artificial nodes NODES+i i=1,2,3,4 and construct first
!     two (artificial) elements.
!
   xmin = X(1)
   xmax = X(1)
   ymin = Y(1)
   ymax = Y(1)
   DO i = 2 , Nodes
      xmin = min(xmin,X(i))
      xmax = max(xmax,X(i))
      ymin = min(ymin,Y(i))
      ymax = max(ymax,Y(i))
   ENDDO
   dx = xmax - xmin
   dy = ymax - ymin
   xl = xmin - 4.0*dx
   xr = xmax + 4.0*dx
   yl = ymin - 4.0*dy
   yr = ymax + 4.0*dy
   X(Nodes+1) = xl
   Y(Nodes+1) = yl
   X(Nodes+2) = xl
   Y(Nodes+2) = yr
   X(Nodes+3) = xr
   Y(Nodes+3) = yr
   X(Nodes+4) = xr
   Y(Nodes+4) = yl
   Elmnts(1,1) = Nodes + 1
   Elmnts(1,2) = Nodes + 2
   Elmnts(1,3) = Nodes + 3
   Elmnts(2,1) = Nodes + 3
   Elmnts(2,2) = Nodes + 4
   Elmnts(2,3) = Nodes + 1
   Nelmnt = 2
   DO in = 1 , Nodes
!
!     add one mesh point at a time and remesh locally if necessary
!
      ndel = 0
      newel = 0
      DO ie = 1 , Nelmnt
!
!     is point in insided circumcircle of element IE ?
!
         i1 = Elmnts(ie,1)
         i2 = Elmnts(ie,2)
         i3 = Elmnts(ie,3)
         x2 = X(i2) - X(i1)
         x3 = X(i3) - X(i1)
         y2 = Y(i2) - Y(i1)
         y3 = Y(i3) - Y(i1)
         z = (x2*(x2-x3)+y2*(y2-y3))/(y2*x3-y3*x2)
         cx = 0.5*(x3-z*y3)
         cy = 0.5*(y3+z*x3)
         r2 = cx**2 + cy**2
         rn2 = ((X(in)-X(i1)-cx)**2+(Y(in)-Y(i1)-cy)**2)
         IF ( rn2<=r2 ) THEN
!
!     yes it is inside,create new elements and mark old for deletion.
!
            DO j = 1 , 3
               DO k = 1 , 3
                  Elmnts(Nelmnt+newel+j,k) = Elmnts(ie,k)
               ENDDO
               Elmnts(Nelmnt+newel+j,j) = in
            ENDDO
            newel = newel + 3
            Elmnts(ie,1) = 0
            ndel = ndel + 1
         ENDIF
!
      ENDDO
!
!     If IN was inside circumcircle of more than 1 element then will
!     have created 2 identical new elements: delete them both.
!
      IF ( ndel>1 ) THEN
         DO ie = Nelmnt + 1 , Nelmnt + newel - 1
            DO je = ie + 1 , Nelmnt + newel
               match = 0
               DO k = 1 , 3
                  DO l = 1 , 3
                     IF ( Elmnts(ie,k)==Elmnts(je,l) ) match = match + 1
                  ENDDO
               ENDDO
               IF ( match==3 ) THEN
                  Elmnts(ie,1) = 0
                  Elmnts(je,1) = 0
                  ndel = ndel + 2
               ENDIF
            ENDDO
         ENDDO
      ENDIF
!
!     delete any elements
!
      nn = Nelmnt + newel
      ie = 1
      SPAG_Loop_2_1: DO
         IF ( Elmnts(ie,1)==0 ) THEN
            DO j = ie , nn - 1
               DO k = 1 , 3
                  Elmnts(j,k) = Elmnts(j+1,k)
               ENDDO
            ENDDO
            nn = nn - 1
            ie = ie - 1
         ENDIF
         ie = ie + 1
         IF ( ie>nn ) THEN
            Nelmnt = nn
            EXIT SPAG_Loop_2_1
         ENDIF
      ENDDO SPAG_Loop_2_1
   ENDDO
!
!     finally remove elements containing artificial nodes
!
   ie = 1
   SPAG_Loop_1_2: DO
      nart = 0
      DO l = 1 , 3
         IF ( Elmnts(ie,l)>Nodes ) nart = nart + 1
      ENDDO
      IF ( nart>0 ) THEN
         DO j = ie , nn - 1
            DO k = 1 , 3
               Elmnts(j,k) = Elmnts(j+1,k)
            ENDDO
         ENDDO
         Nelmnt = Nelmnt - 1
         ie = ie - 1
      ENDIF
      ie = ie + 1
      IF ( ie>Nelmnt ) EXIT SPAG_Loop_1_2
   ENDDO SPAG_Loop_1_2
END SUBROUTINE dlauny
!*==estif3.f90 processed by SPAG 8.04RA 15:47 21 Feb 2024
!!SPAG Open source Personal, Educational or Academic User Clemson University  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
!-----------------------------------------------------------------------
SUBROUTINE estif3(Nel,Ske,Det,Xe,Ye,Xyke,Ierr)
!-----------------------------------------------------------------------
! this subroutine constructs the element stiffness matrix for heat
! condution problem
!
!                  - Div ( K(x,y) Grad u ) = f
!                    u = 0 on boundary
!
! using 3-node triangular elements arguments:
! nel        = element number
! ske        = element stiffness matrix
! [fe        = element load vector] unused here.
! det        = 2*area of the triangle
! xy, ye= coordinates of the three nodal points in an element.
! xyke  = material constants (kxx, kxy, kyx, kyy)
!
!------------------------------------------------------------------------
   USE ISO_FORTRAN_ENV                 
   IMPLICIT NONE
!
! Dummy argument declarations rewritten by SPAG
!
   INTEGER :: Nel
   REAL(REAL64) , INTENT(OUT) , DIMENSION(3,3) :: Ske
   REAL(REAL64) :: Det
   REAL(REAL64) , DIMENSION(3) :: Xe
   REAL(REAL64) , DIMENSION(3) :: Ye
   REAL(REAL64) , INTENT(IN) , DIMENSION(2,2) :: Xyke
   INTEGER :: Ierr
!
! Local variable declarations rewritten by SPAG
!
   REAL(REAL64) :: area , t
   REAL(REAL64) , DIMENSION(3,2) :: dn
   INTEGER :: i , j , k , l
   EXTERNAL gradi3
!
! End of declarations rewritten by SPAG
!
!
! initialize
!
   area = 0.5*Det
!
   DO i = 1 , 3
      DO j = 1 , 3
         Ske(i,j) = 0.0D0
      ENDDO
   ENDDO
!
! get first gradient of shape function
!
   CALL gradi3(Nel,Xe,Ye,dn,Det,Ierr)
   IF ( Ierr/=0 ) RETURN
!
   DO i = 1 , 3
      DO j = 1 , 3
         t = 0.0D0
         DO k = 1 , 2
            DO l = 1 , 2
               t = t + Xyke(k,l)*dn(i,k)*dn(j,l)
            ENDDO
         ENDDO
         Ske(i,j) = t*area
      ENDDO
   ENDDO
!
END SUBROUTINE estif3
!*==gradi3.f90 processed by SPAG 8.04RA 15:47 21 Feb 2024
!!SPAG Open source Personal, Educational or Academic User Clemson University  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
!-------------------------------------------------------
SUBROUTINE gradi3(Nel,Xe,Ye,Dn,Det,Ierr)
!-------------------------------------------------------
! constructs the first derivative of the shape functions.
! arguments:
! nel       = element nuumber
! xy, ye    = coordinates of the three nodal points in an element.
! dn        = gradients (1-st derivatives) of the shape functions.
! area      = area of the triangle
!
!-------------------------------------------------------
   USE ISO_FORTRAN_ENV                 
   IMPLICIT NONE
!
! Dummy argument declarations rewritten by SPAG
!
   INTEGER , INTENT(IN) :: Nel
   REAL(REAL64) , INTENT(IN) , DIMENSION(3) :: Xe
   REAL(REAL64) , INTENT(IN) , DIMENSION(3) :: Ye
   REAL(REAL64) , INTENT(OUT) , DIMENSION(3,2) :: Dn
   REAL(REAL64) , INTENT(IN) :: Det
   INTEGER , INTENT(OUT) :: Ierr
!
! Local variable declarations rewritten by SPAG
!
   REAL(REAL64) , SAVE :: eps
!
! End of declarations rewritten by SPAG
!
   DATA eps/1.D-17/
! compute area
   Ierr = 0
   IF ( Det<=eps ) THEN
!
      Ierr = 3
      PRINT * , 'ERROR:negative area encountered at elmt: ' , Nel
      RETURN
   ENDIF
!
   Dn(1,1) = (Ye(2)-Ye(3))/Det
   Dn(2,1) = (Ye(3)-Ye(1))/Det
   Dn(3,1) = (Ye(1)-Ye(2))/Det
   Dn(1,2) = (Xe(3)-Xe(2))/Det
   Dn(2,2) = (Xe(1)-Xe(3))/Det
   Dn(3,2) = (Xe(2)-Xe(1))/Det
!
   RETURN
!        write(iout,*) det,(xe(i),ye(i),i=1,3)
END SUBROUTINE gradi3
!*==hsourc.f90 processed by SPAG 8.04RA 15:47 21 Feb 2024
!!SPAG Open source Personal, Educational or Academic User Clemson University  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
!-----------------------------------------------------------------------
SUBROUTINE hsourc(Indic,Nelx,Node,X,Y,Ijk,Fs,F)
   USE ISO_FORTRAN_ENV                 
   IMPLICIT NONE
!
! Dummy argument declarations rewritten by SPAG
!
   INTEGER , INTENT(IN) :: Node
   INTEGER , INTENT(IN) :: Indic
   INTEGER , INTENT(IN) :: Nelx
   REAL(REAL64) , INTENT(IN) , DIMENSION(*) :: X
   REAL(REAL64) , INTENT(IN) , DIMENSION(*) :: Y
   INTEGER , INTENT(IN) , DIMENSION(Node,*) :: Ijk
   REAL(REAL64) , INTENT(IN) , DIMENSION(*) :: Fs
   REAL(REAL64) , INTENT(INOUT) , DIMENSION(*) :: F
!
! Local variable declarations rewritten by SPAG
!
   REAL(REAL64) :: areao3 , det
   INTEGER :: i , ii , j , jnod , ka , nel
   REAL(REAL64) , DIMENSION(3) :: xe , ye
!
! End of declarations rewritten by SPAG
!
!
! generates the load vector f in assembled/unassembled form from the
! the element contributions fs.
! indic = indicates if f is to be assembled (1) or not (zero)
! note: f(*) not initilazed. because might use values from boundary
! conditions.
!
   jnod = 0
   DO nel = 1 , Nelx
!
! get coordinates of nodal points
!
      DO i = 1 , Node
         j = Ijk(i,nel)
         xe(i) = X(j)
         ye(i) = Y(j)
      ENDDO
!
! compute determinant
!
      det = xe(2)*(ye(3)-ye(1)) + xe(3)*(ye(1)-ye(2)) + xe(1)*(ye(2)-ye(3))
! area3 = area/3
      areao3 = det/6.0
!
! contributions to nodes in the element
!
      IF ( Indic==0 ) THEN
         DO ka = 1 , Node
            jnod = jnod + 1
            F(jnod) = Fs(nel)*areao3
         ENDDO
      ELSE
         DO ka = 1 , Node
            ii = Ijk(ka,nel)
            F(ii) = F(ii) + Fs(nel)*areao3
         ENDDO
      ENDIF
!
   ENDDO
END SUBROUTINE hsourc
!*==bound.f90 processed by SPAG 8.04RA 15:47 21 Feb 2024
!!SPAG Open source Personal, Educational or Academic User Clemson University  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
!----- end of hsourc ---------------------------------------------------
!-----------------------------------------------------------------------
SUBROUTINE bound(Nx,Nelx,Ijk,Nodcode,Node,Nint,Iperm,X,Y,Wk,Iwk)
!-----------------------------------------------------------------------
! this routine counts the number of boundary points and
! reorders the points in such a way that the boundary nodes
! are last.
!
! nx, nelx, ijk, nodcode, node: see other subroutines
! iperm = permutation array from old orderin to new ordering,
! iwk   = reverse permutation array or return.
! wk	= real work array
! On return
! x, y, nodecode, are permuted
! ijk  is updated according to new oerdering.
! nint = number of interior points.
!
!-----------------------------------------------------------------------
   USE ISO_FORTRAN_ENV                 
   IMPLICIT NONE
!
! Dummy argument declarations rewritten by SPAG
!
   INTEGER , INTENT(IN) :: Node
   INTEGER , INTENT(IN) :: Nx
   INTEGER , INTENT(IN) :: Nelx
   INTEGER , INTENT(INOUT) , DIMENSION(Node,1) :: Ijk
   INTEGER , INTENT(INOUT) , DIMENSION(1) :: Nodcode
   INTEGER , INTENT(INOUT) :: Nint
   INTEGER , INTENT(INOUT) , DIMENSION(1) :: Iperm
   REAL(REAL64) , INTENT(INOUT) , DIMENSION(1) :: X
   REAL(REAL64) , INTENT(INOUT) , DIMENSION(1) :: Y
   REAL(REAL64) , INTENT(INOUT) , DIMENSION(1) :: Wk
   INTEGER , INTENT(INOUT) , DIMENSION(1) :: Iwk
!
! Local variable declarations rewritten by SPAG
!
   INTEGER :: j , k , knod , nbound , nel
!
! End of declarations rewritten by SPAG
!
 
!     put all boundary points at the end, backwards
   Nint = 1
   nbound = Nx
   DO j = 1 , Nx
      IF ( Nodcode(j)==0 ) THEN
         Iperm(Nint) = j
         Nint = Nint + 1
      ELSE
         Iperm(nbound) = j
         nbound = nbound - 1
      ENDIF
   ENDDO
!-------------------------------------------------------------------
   Nint = Nint - 1
!
! permute x's
!
   DO k = 1 , Nx
      Wk(k) = X(k)
   ENDDO
   DO k = 1 , Nx
      X(k) = Wk(Iperm(k))
   ENDDO
!
!     permute the y's
!
   DO k = 1 , Nx
      Wk(k) = Y(k)
   ENDDO
   DO k = 1 , Nx
      Y(k) = Wk(Iperm(k))
   ENDDO
!
!     permute the boundary information
!
   DO k = 1 , Nx
      Iwk(k) = Nodcode(k)
   ENDDO
   DO k = 1 , Nx
      Nodcode(k) = Iwk(Iperm(k))
   ENDDO
!
!     get reverse permutation
!
   DO k = 1 , Nx
      Iwk(Iperm(k)) = k
   ENDDO
!
!     update the elements connectivity matrix
!
   DO nel = 1 , Nelx
      DO j = 1 , Node
         knod = Ijk(j,nel)
         Ijk(j,nel) = Iwk(knod)
      ENDDO
   ENDDO
END SUBROUTINE bound
!*==symbound.f90 processed by SPAG 8.04RA 15:47 21 Feb 2024
!!SPAG Open source Personal, Educational or Academic User Clemson University  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
!-----------------------------------------------------------------------
SUBROUTINE symbound(Nx,Nelx,Ijk,Nodcode,Node,Nint,Iperm,Iwk)
!-----------------------------------------------------------------------
!     this routine is a symbolic version of routine bound.
!
!   nx, nelx, ijk, nodcode, node: see other subroutines
!   iperm = permutation array from old orderin to new ordering,
!   iwk   = reverse permutation array or return.
!   [wk	  = real work array] unused
!   On return
!   ijk   = is updated according to new oerdering.
!   nint  = number of interior points.
!
!-----------------------------------------------------------------------
   USE ISO_FORTRAN_ENV                 
   IMPLICIT NONE
!
! Dummy argument declarations rewritten by SPAG
!
   INTEGER , INTENT(IN) :: Node
   INTEGER , INTENT(IN) :: Nx
   INTEGER , INTENT(IN) :: Nelx
   INTEGER , INTENT(INOUT) , DIMENSION(Node,*) :: Ijk
   INTEGER , INTENT(INOUT) , DIMENSION(*) :: Nodcode
   INTEGER , INTENT(INOUT) :: Nint
   INTEGER , INTENT(INOUT) , DIMENSION(*) :: Iperm
   INTEGER , INTENT(INOUT) , DIMENSION(*) :: Iwk
!
! Local variable declarations rewritten by SPAG
!
   INTEGER :: j , k , knod , nbound , nel
!
! End of declarations rewritten by SPAG
!
 
! put all boundary points at the end, backwards
   Nint = 1
   nbound = Nx
   DO j = 1 , Nx
      IF ( Nodcode(j)==0 ) THEN
         Iperm(Nint) = j
         Nint = Nint + 1
      ELSE
         Iperm(nbound) = j
         nbound = nbound - 1
      ENDIF
   ENDDO
!-------------------------------------------------------------------
   Nint = Nint - 1
!
! permute the boundary information
!
   DO k = 1 , Nx
      Iwk(k) = Nodcode(k)
   ENDDO
   DO k = 1 , Nx
      Nodcode(k) = Iwk(Iperm(k))
   ENDDO
!
! get reverse permutation
!
   DO k = 1 , Nx
      Iwk(Iperm(k)) = k
   ENDDO
!
! update the elements connectivity matrix
!
   DO nel = 1 , Nelx
      DO j = 1 , Node
         knod = Ijk(j,nel)
         Ijk(j,nel) = Iwk(knod)
      ENDDO
   ENDDO
END SUBROUTINE symbound
!*==diric.f90 processed by SPAG 8.04RA 15:47 21 Feb 2024
!!SPAG Open source Personal, Educational or Academic User Clemson University  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
!-----------------------------------------------------------------------
SUBROUTINE diric(Nint,A,Ja,Ia)
!--------------------------------------------------------------
! this routine takes into account the boundary conditions
! and removes the unnecessary boundary points.
!--------------------------------------------------------------
   USE ISO_FORTRAN_ENV                 
   IMPLICIT NONE
!
! Dummy argument declarations rewritten by SPAG
!
   INTEGER :: Nint
   REAL(REAL64) , DIMENSION(*) :: A
   INTEGER , DIMENSION(*) :: Ja
   INTEGER , DIMENSION(*) :: Ia
!
! Local variable declarations rewritten by SPAG
!
   INTEGER :: nc , nr
   EXTERNAL submat
!
! End of declarations rewritten by SPAG
!
! call extract from UNARY
   CALL submat(1,1,Nint,1,Nint,A,Ja,Ia,nr,nc,A,Ja,Ia)
   WRITE (*,*) 'nr=' , nr , 'nc=' , nc
!----------- end of diric -------------------------------------
END SUBROUTINE diric
!*==symdiric.f90 processed by SPAG 8.04RA 15:47 21 Feb 2024
!!SPAG Open source Personal, Educational or Academic User Clemson University  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
!-----------------------------------------------------------------------
SUBROUTINE symdiric(Nint,A,Ja,Ia)
!--------------------------------------------------------------
! this routine takes into account the boundary conditions
! and removes the unnecessary boundary points.
!--------------------------------------------------------------
   USE ISO_FORTRAN_ENV                 
   IMPLICIT NONE
!
! Dummy argument declarations rewritten by SPAG
!
   INTEGER :: Nint
   REAL(REAL64) , DIMENSION(*) :: A
   INTEGER , DIMENSION(*) :: Ja
   INTEGER , DIMENSION(*) :: Ia
!
! Local variable declarations rewritten by SPAG
!
   INTEGER :: nc , nr
   EXTERNAL submat
!
! End of declarations rewritten by SPAG
!
!     call submat from UNARY, with job = 0,
!     meaning no movement of real values.
   CALL submat(0,1,Nint,1,Nint,A,Ja,Ia,nr,nc,A,Ja,Ia)
!----------- end of symdiric -------------------------------------
END SUBROUTINE symdiric
!*==cleannods.f90 processed by SPAG 8.04RA 15:47 21 Feb 2024
!!SPAG Open source Personal, Educational or Academic User Clemson University  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
!-----------------------------------------------------------------------
SUBROUTINE cleannods(Nx,X,Y,Nelx,Ijk,Node,Nodcode,Iperm)
   USE ISO_FORTRAN_ENV                 
   IMPLICIT NONE
!
! Dummy argument declarations rewritten by SPAG
!
   INTEGER , INTENT(INOUT) :: Nx
   INTEGER , INTENT(IN) :: Nelx
   INTEGER , INTENT(IN) :: Node
   REAL(REAL64) , INTENT(INOUT) , DIMENSION(Nx) :: X
   REAL(REAL64) , INTENT(INOUT) , DIMENSION(Nx) :: Y
   INTEGER , INTENT(INOUT) , DIMENSION(Node,Nelx) :: Ijk
   INTEGER , INTENT(INOUT) , DIMENSION(*) :: Nodcode
   INTEGER , INTENT(INOUT) , DIMENSION(Nx) :: Iperm
!
! Local variable declarations rewritten by SPAG
!
   INTEGER :: i , indx , j , k , nel
!
! End of declarations rewritten by SPAG
!
!      implicit none
!-----------------------------------------------------------------------
!     this routine removes the nodes that do not belong to any element
!     (spurious points) and relabels the ijk array accordingly.
!-----------------------------------------------------------------------
!
   DO j = 1 , Nx
      Iperm(j) = 0
   ENDDO
!
   DO nel = 1 , Nelx
      DO i = 1 , Node
         k = Ijk(i,nel)
         Iperm(k) = nel
      ENDDO
   ENDDO
!
   indx = 0
   DO j = 1 , Nx
      IF ( Iperm(j)/=0 ) THEN
         indx = indx + 1
         Iperm(indx) = j
         X(indx) = X(j)
         Y(indx) = Y(j)
         Nodcode(indx) = Nodcode(j)
      ENDIF
   ENDDO
!
!     update nx
!
   Nx = indx
!
!     old number to new numbers
!
   DO j = 1 , Nx
      Iperm(Nx+Iperm(j)) = j
   ENDDO
!
!
!     change all node numbers in ijk
!
   DO nel = 1 , Nelx
      DO i = 1 , Node
         k = Ijk(i,nel)
         k = Iperm(Nx+k)
         Ijk(i,nel) = k
      ENDDO
   ENDDO
!-----------------------------------------------------------------------
!-----end-of-cleannod---------------------------------------------------
END SUBROUTINE cleannods
!*==cleanel.f90 processed by SPAG 8.04RA 15:47 21 Feb 2024
!!SPAG Open source Personal, Educational or Academic User Clemson University  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
!-----------------------------------------------------------------------
SUBROUTINE cleanel(Nelx,Ijk,Node,Nodcode,Nodexc)
   IMPLICIT NONE
!
! Dummy argument declarations rewritten by SPAG
!
   INTEGER , INTENT(INOUT) :: Nelx
   INTEGER , INTENT(IN) :: Node
   INTEGER , INTENT(INOUT) , DIMENSION(Node,Nelx) :: Ijk
   INTEGER , INTENT(IN) , DIMENSION(*) :: Nodcode
   INTEGER , INTENT(IN) :: Nodexc
!
! Local variable declarations rewritten by SPAG
!
   LOGICAL :: exclude
   INTEGER :: i , k , nel
!
! End of declarations rewritten by SPAG
!
!      implicit none
!-----------------------------------------------------------------------
!     this routine remove certain types of elements from the mesh
!     An element whose nodes are all labelled by the same label
!     nodexc are removed. nelx is changed accordingly on return.
!-----------------------------------------------------------------------
   nel = 1
   SPAG_Loop_1_1: DO
      exclude = .TRUE.
      DO i = 1 , Node
         k = Ijk(i,nel)
         exclude = (exclude .AND. Nodcode(k)==Nodexc)
      ENDDO
!
      IF ( exclude ) THEN
         DO i = 1 , Node
            Ijk(i,nel) = Ijk(i,Nelx)
         ENDDO
         Nelx = Nelx - 1
      ELSE
         nel = nel + 1
      ENDIF
      IF ( nel>Nelx ) EXIT SPAG_Loop_1_1
   ENDDO SPAG_Loop_1_1
!-----------------------------------------------------------------------
!-----end-of-cleanel----------------------------------------------------
END SUBROUTINE cleanel
!*==psi.f90 processed by SPAG 8.04RA 15:47 21 Feb 2024
!!SPAG Open source Personal, Educational or Academic User Clemson University  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
!---------------------------------------------------------------------------
!  Piecewise linear fucntions on triangle.
FUNCTION psi(I,R,S)
   USE ISO_FORTRAN_ENV                 
   IMPLICIT NONE
!
! Function and Dummy argument declarations rewritten by SPAG
!
   REAL(REAL64) :: psi
   INTEGER , INTENT(IN) :: I
   REAL(REAL64) , INTENT(IN) :: R
   REAL(REAL64) , INTENT(IN) :: S
!
! End of declarations rewritten by SPAG
!
 
   IF ( I==2 ) THEN
      psi = (R+1.)/2.
      RETURN
   ELSEIF ( I==3 ) THEN
      psi = (S+1.)/2.
      RETURN
   ENDIF
   psi = -(R+S)/2.
   RETURN
END FUNCTION psi
