!*==genfea.f90 processed by SPAG 8.04RA 15:47 21 Feb 2024
!!SPAG Open source Personal, Educational or Academic User Clemson University  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
!----------------------------------------------------------------------c
!                          S P A R S K I T                             c
!----------------------------------------------------------------------c
!      MATRIX GENERATION ROUTINES - FINITE ELEMENT MATRICES            c
!----------------------------------------------------------------------c
! contents:                                                            c
!----------                                                            c
! genfea       : generates finite element matrices in assembled form   c
! genfea_wbc   : generates finite element matrices in assembled form   c
!                without applying the boundary conditions              c
! genfeu       : generates finite element matrices in unassembled form c
! genfeu_wbc   : generates finite element matrices in unassembled form c
!                without applying the boundary conditions              c
! genfeu_lstif : generates finite element matrices in unassembled form c
!                using the lstif problem appearing in elmtlib2.f       c
! assmb1       : assembles an unassembled matrix (produced by genfeu)  c
!----------------------------------------------------------------------c
SUBROUTINE genfea(Nx,Nelx,Node,Job,X,Y,Ijk,Nodcode,Fs,Nint,A,Ja,Ia,F,Iwk,Jwk,Ierr,Xyk)
   USE ISO_FORTRAN_ENV                 
   IMPLICIT NONE
!
! Dummy argument declarations rewritten by SPAG
!
   INTEGER :: Node
   INTEGER :: Nx
   INTEGER :: Nelx
   INTEGER , INTENT(IN) :: Job
   REAL(REAL64) , DIMENSION(*) :: X
   REAL(REAL64) , DIMENSION(*) :: Y
   INTEGER , DIMENSION(Node,*) :: Ijk
   INTEGER , DIMENSION(*) :: Nodcode
   REAL(REAL64) , DIMENSION(*) :: Fs
   INTEGER :: Nint
   REAL(REAL64) , DIMENSION(*) :: A
   INTEGER , DIMENSION(*) :: Ja
   INTEGER , DIMENSION(*) :: Ia
   REAL(REAL64) , DIMENSION(*) :: F
   INTEGER , DIMENSION(*) :: Iwk
   INTEGER , DIMENSION(*) :: Jwk
   INTEGER :: Ierr
   REAL , EXTERNAL :: Xyk
!
! Local variable declarations rewritten by SPAG
!
   INTEGER :: indic
   EXTERNAL assmbo , bound , diric , funb , func , fung , hsourc
!
! End of declarations rewritten by SPAG
!
!-----------------------------------------------------------------------
! this subroutine generates a finite element matrix in assembled form.
! the matrix is assembled in compressed sparse row format. See genfeu
! for matrices in unassembled form. The user must provide the grid,
! (coordinates x, y and connectivity matrix ijk) as well as some
! information on the nodes (nodcode) and the material properties
! (the function K(x,y) above) in the form of a subroutine xyk.
!----------------------------------------------------------------------
!
! on entry:
! ---------
!
! nx	    = integer . the number of nodes in the grid .
! nelx	    = integer . the number of elements in the grid.
! node      = integer = the number of nodes per element (should be
!             set to three in this version). also the first dimension
!             of ijk
! job	    = integer. If job=0, it is assumed that there is no heat
!             source (i.e. fs = 0) and the right hand side
!             produced will therefore be a zero vector.
!             If job = 1 on entry then the contributions from the
!             heat source in each element are taken into account.
!
! x, y      = two real arrays containing the coordinates of the nodes.
!
! ijk       =  an integer array containing the connectivity matrix.
!              ijk(i,nel), i=1,2,..node, is the list of the nodes
!              constituting the element nel, ans listed in
!              counter clockwise order.
!
! nodcode   = an integer array containing the boundary information for
!             each node with the following meaning.
!	nodcode(i) = 0 -->  node i is internal
!	nodcode(i) = 1 -->  node i is a boundary but not a corner point
!	nodcode(i) = 2 -->  node i is a corner node. [This node and the
!             corresponmding element are discarded.]
!
! fs	    = real array of length nelx on entry containing the heat
!             source for each element (job = 1 only)
!
! xyk	    = subroutine defining the material properties at each
!	      element. Form:
! 	      call xyk(nel,xyke,x,y,ijk,node) with on return
!             xyke =  material constant matrices.
!	      for each element nel, xyke(1,nel),xyke(2,nel)
!             and xyke(3,nel) represent the constants
!             K11, K22, and K12 at that element.
!
! on return
! ---------
! nint	    = integer. The number of active (nonboundary) nodes. Also
!             equal to the dimension of the assembled matrix.
!
! a, ja, ia = assembled matrix in compressed sparse row format.
!
! f	    = real array containing the right hand for the linears
!             system to solve.
!
! ierr	    = integer. Error message. If (ierr .ne. 0) on return
!             it means that one of the elements has a negative or zero
!             area probably because of a bad ordering of the nodes
!             (see ijk above). Use the subroutine chkelmt to reorder
!             the nodes properly if necessary.
! iwk, jwk  = two integer work arrays of length nx each.
!
!-----------------------------------------------------------------------
!
   Ierr = 0
!
!     take into boundary conditions to remove boundary nodes.
!
   CALL bound(Nx,Nelx,Ijk,Nodcode,Node,Nint,Jwk,X,Y,F,Iwk)
!
!     assemble the matrix
!
   CALL assmbo(Nx,Nelx,Node,Ijk,Nodcode,X,Y,A,Ja,Ia,F,Iwk,Jwk,Ierr,Xyk,funb,func,fung)
!
!     if applicable (job .eq. 1) get heat source function
!
   indic = 1
   IF ( Job==1 ) CALL hsourc(indic,Nelx,Node,X,Y,Ijk,Fs,F)
!
!     call diric for Dirichlet conditions
!
   CALL diric(Nint,A,Ja,Ia)
!     done
!------end of genfea ---------------------------------------------------
!-----------------------------------------------------------------------
END SUBROUTINE genfea
!*==genfea_wbc.f90 processed by SPAG 8.04RA 15:47 21 Feb 2024
!!SPAG Open source Personal, Educational or Academic User Clemson University  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
SUBROUTINE genfea_wbc(Nx,Nelx,Node,Job,X,Y,Ijk,Nodcode,Fs,A,Ja,Ia,F,Iwk,Jwk,Ierr,Xyk)
   USE ISO_FORTRAN_ENV                 
   IMPLICIT NONE
!
! Dummy argument declarations rewritten by SPAG
!
   INTEGER :: Node
   INTEGER :: Nx
   INTEGER :: Nelx
   INTEGER , INTENT(IN) :: Job
   REAL(REAL64) , DIMENSION(*) :: X
   REAL(REAL64) , DIMENSION(*) :: Y
   INTEGER , DIMENSION(Node,*) :: Ijk
   INTEGER , DIMENSION(*) :: Nodcode
   REAL(REAL64) , DIMENSION(*) :: Fs
   REAL(REAL64) , DIMENSION(*) :: A
   INTEGER , DIMENSION(*) :: Ja
   INTEGER , DIMENSION(*) :: Ia
   REAL(REAL64) , DIMENSION(*) :: F
   INTEGER , DIMENSION(*) :: Iwk
   INTEGER , DIMENSION(*) :: Jwk
   INTEGER :: Ierr
   REAL , EXTERNAL :: Xyk
!
! Local variable declarations rewritten by SPAG
!
   INTEGER :: indic
   EXTERNAL assmbo , funb , func , fung , hsourc
!
! End of declarations rewritten by SPAG
!
!-----------------------------------------------------------------------
! this subroutine generates a finite element matrix in assembled form.
! the matrix is assembled in compressed sparse row format. See genfeu
! for matrices in unassembled form. The user must provide the grid,
! (coordinates x, y and connectivity matrix ijk) as well as some
! information on the nodes (nodcode) and the material properties
! (the function K(x,y) above) in the form of a subroutine xyk.
!----------------------------------------------------------------------
! Irene Moulitsas, moulitsa@cs.umn.edu  : It does not apply boundary
!                 conditions; variable nint is eliminated
!----------------------------------------------------------------------
!
! on entry:
! ---------
!
! nx	    = integer . the number of nodes in the grid .
! nelx	    = integer . the number of elements in the grid.
! node      = integer = the number of nodes per element (should be
!             set to three in this version). also the first dimension
!             of ijk
! job	    = integer. If job=0, it is assumed that there is no heat
!             source (i.e. fs = 0) and the right hand side
!             produced will therefore be a zero vector.
!             If job = 1 on entry then the contributions from the
!             heat source in each element are taken into account.
!
! x, y      = two real arrays containing the coordinates of the nodes.
!
! ijk       =  an integer array containing the connectivity matrix.
!              ijk(i,nel), i=1,2,..node, is the list of the nodes
!              constituting the element nel, ans listed in
!              counter clockwise order.
!
! nodcode   = an integer array containing the boundary information for
!             each node with the following meaning.
!	nodcode(i) = 0 -->  node i is internal
!	nodcode(i) = 1 -->  node i is a boundary but not a corner point
!	nodcode(i) = 2 -->  node i is a corner node. [This node and the
!             corresponmding element are discarded.]
!
! fs	    = real array of length nelx on entry containing the heat
!             source for each element (job = 1 only)
!
! xyk	    = subroutine defining the material properties at each
!	      element. Form:
! 	      call xyk(nel,xyke,x,y,ijk,node) with on return
!             xyke =  material constant matrices.
!	      for each element nel, xyke(1,nel),xyke(2,nel)
!             and xyke(3,nel) represent the constants
!             K11, K22, and K12 at that element.
!
! on return
! ---------
! a, ja, ia = assembled matrix in compressed sparse row format.
!
! f	    = real array containing the right hand for the linears
!             system to solve.
!
! ierr	    = integer. Error message. If (ierr .ne. 0) on return
!             it means that one of the elements has a negative or zero
!             area probably because of a bad ordering of the nodes
!             (see ijk above). Use the subroutine chkelmt to reorder
!             the nodes properly if necessary.
! iwk, jwk  = two integer work arrays of length nx each.
!
!-----------------------------------------------------------------------
!
   Ierr = 0
!
!     assemble the matrix
!
   CALL assmbo(Nx,Nelx,Node,Ijk,Nodcode,X,Y,A,Ja,Ia,F,Iwk,Jwk,Ierr,Xyk,funb,func,fung)
!
!     if applicable (job .eq. 1) get heat source function
!
   indic = 1
   IF ( Job==1 ) CALL hsourc(indic,Nelx,Node,X,Y,Ijk,Fs,F)
!
!     done
!------end of genfea_wbc -----------------------------------------------
!-----------------------------------------------------------------------
END SUBROUTINE genfea_wbc
!*==genfeu.f90 processed by SPAG 8.04RA 15:47 21 Feb 2024
!!SPAG Open source Personal, Educational or Academic User Clemson University  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
!-----------------------------------------------------------------------
SUBROUTINE genfeu(Nx,Nelx,Node,Job,X,Y,Ijk,Nodcode,Fs,Nint,A,Na,F,Iwk,Jwk,Ierr,Xyk)
   USE ISO_FORTRAN_ENV                 
   IMPLICIT NONE
!
! Dummy argument declarations rewritten by SPAG
!
   INTEGER :: Node
   INTEGER :: Na
   INTEGER :: Nx
   INTEGER :: Nelx
   INTEGER , INTENT(IN) :: Job
   REAL(REAL64) , DIMENSION(*) :: X
   REAL(REAL64) , DIMENSION(*) :: Y
   INTEGER , DIMENSION(Node,*) :: Ijk
   INTEGER , DIMENSION(*) :: Nodcode
   REAL(REAL64) , DIMENSION(*) :: Fs
   INTEGER :: Nint
   REAL(REAL64) , DIMENSION(Na,Node,Node) :: A
   REAL(REAL64) , DIMENSION(*) :: F
   INTEGER , DIMENSION(*) :: Iwk
   INTEGER , DIMENSION(*) :: Jwk
   INTEGER :: Ierr
   REAL , EXTERNAL :: Xyk
!
! Local variable declarations rewritten by SPAG
!
   INTEGER :: indic
   EXTERNAL bound , hsourc , unassbl
!
! End of declarations rewritten by SPAG
!
!-----------------------------------------------------------------------
! this subroutine generates finite element matrices for heat
! condution problem
!
!                  - Div ( K(x,y) Grad u ) = f
!                    u = 0 on boundary
!
! (with Dirichlet boundary conditions). The matrix is returned
! in unassembled form. The user must provide the grid,
! (coordinates x, y and connectivity matrix ijk) as well as some
! information on the nodes (nodcode) and the material properties
! (the function K(x,y) above) in the form of a subroutine xyk.
!
!----------------------------------------------------------------------
!
! on entry:
! ---------
!
! nx	    = integer . the number of nodes in the grid .
! nelx	    = integer . the number of elements in the grid.
! node      = integer = the number of nodes per element (should be
!             set to three in this version). also the first dimension
!             of ijk
! job	    = integer. If job=0, it is assumed that there is no heat
!             source (i.e. fs = 0) and the right hand side
!             produced will therefore be a zero vector.
!             If job = 1 on entry then the contributions from the
!             heat source in each element are taken into account.
!
! na	    = integer. The first dimension of the array a.
!             a is declared as an array of dimension a(na,node,node).
!
! x, y      = two real arrays containing the coordinates of the nodes.
!
! ijk       =  an integer array containing the connectivity matrix.
!              ijk(i,nel), i=1,2,..node, is the list of the nodes
!              constituting the element nel, ans listed in
!              counter clockwise order.
!
! xyk	    = subroutine defining the material properties at each
!	      element. Form:
! 	      call xyk(nel,xyke,x,y,ijk,node) with on return
!             xyke =  material constant matrices.
!	      for each element nel, xyke(1,nel),xyke(2,nel)
!             and xyke(3,nel) represent the constants
!             K11, K22, and K12 at that element.
!
! nodcode   = an integer array containing the boundary information for
!             each node with the following meaning.
!	nodcode(i) = 0 -->  node i is internal
!	nodcode(i) = 1 -->  node i is a boundary but not a corner point
!	nodcode(i) = 2 -->  node i is a corner node. [This node and the
!             corresponmding element are discarded.]
!
! fs	    = real array of length nelx on entry containing the heat
!             source for each element (job = 1 only)
!
! on return
! ---------
! nint	    = integer. The number of active (nonboundary) nodes. Also
!             equal to the dimension of the assembled matrix.
!
! a         = matrix in unassembled form. a(nel,*,*) contains the
!             element matrix for element nel.
!
! f	    = real array containing the right hand for the linears
!             system to solve, in assembled form.
!
! ierr	    = integer. Error message. If (ierr .ne. 0) on return
!             it means that one of the elements has a negative or zero
!             area probably because of a bad ordering of the nodes
!             (see ijk above). Use the subroutine chkelmt to reorder
!             the nodes properly if necessary.
! iwk, jwk  = two integer work arrays of length nx each.
!
!-----------------------------------------------------------------------
!
   Ierr = 0
!
!     take boundary conditions into account to move boundary nodes to
!     the end..
!
   CALL bound(Nx,Nelx,Ijk,Nodcode,Node,Nint,Jwk,X,Y,F,Iwk)
!
!     assemble the matrix
!
   CALL unassbl(A,Na,F,Nx,Nelx,Ijk,Node,X,Y,Ierr,Xyk)
!
!     if applicable (job .eq. 1) get heat source function
!
   indic = 0
   IF ( Job==1 ) CALL hsourc(indic,Nelx,Node,X,Y,Ijk,Fs,F)
!
!     done
!
END SUBROUTINE genfeu
!*==genfeu_wbc.f90 processed by SPAG 8.04RA 15:47 21 Feb 2024
!!SPAG Open source Personal, Educational or Academic User Clemson University  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
!----- end of genfeu ----------------------------------------------------
SUBROUTINE genfeu_wbc(Nx,Nelx,Node,Job,X,Y,Ijk,Fs,A,Na,F,Ierr,Xyk)
   USE ISO_FORTRAN_ENV                 
   IMPLICIT NONE
!
! Dummy argument declarations rewritten by SPAG
!
   INTEGER :: Node
   INTEGER :: Na
   INTEGER :: Nx
   INTEGER :: Nelx
   INTEGER , INTENT(IN) :: Job
   REAL(REAL64) , DIMENSION(*) :: X
   REAL(REAL64) , DIMENSION(*) :: Y
   INTEGER , DIMENSION(Node,*) :: Ijk
   REAL(REAL64) , DIMENSION(*) :: Fs
   REAL(REAL64) , DIMENSION(Na,Node,Node) :: A
   REAL(REAL64) , DIMENSION(*) :: F
   INTEGER :: Ierr
   REAL , EXTERNAL :: Xyk
!
! Local variable declarations rewritten by SPAG
!
   INTEGER :: indic
   EXTERNAL hsourc , unassbl
!
! End of declarations rewritten by SPAG
!
!-----------------------------------------------------------------------
! this subroutine generates finite element matrices for heat
! condution problem
!
!                  - Div ( K(x,y) Grad u ) = f
!                    u = 0 on boundary
!
! (with Dirichlet boundary conditions). The matrix is returned
! in unassembled form. The user must provide the grid,
! (coordinates x, y and connectivity matrix ijk) as well as some
! information on the nodes (nodcode) and the material properties
! (the function K(x,y) above) in the form of a subroutine xyk.
!
!----------------------------------------------------------------------
! moulitsa@cs   : It does not apply boundary conditions
!                 variable nint is eliminated
!----------------------------------------------------------------------
!
! on entry:
! ---------
!
! nx	    = integer . the number of nodes in the grid .
! nelx	    = integer . the number of elements in the grid.
! node      = integer = the number of nodes per element (should be
!             set to three in this version). also the first dimension
!             of ijk
! job	    = integer. If job=0, it is assumed that there is no heat
!             source (i.e. fs = 0) and the right hand side
!             produced will therefore be a zero vector.
!             If job = 1 on entry then the contributions from the
!             heat source in each element are taken into account.
!
! na	    = integer. The first dimension of the array a.
!             a is declared as an array of dimension a(na,node,node).
!
! x, y      = two real arrays containing the coordinates of the nodes.
!
! ijk       =  an integer array containing the connectivity matrix.
!              ijk(i,nel), i=1,2,..node, is the list of the nodes
!              constituting the element nel, ans listed in
!              counter clockwise order.
!
! xyk	    = subroutine defining the material properties at each
!	      element. Form:
! 	      call xyk(nel,xyke,x,y,ijk,node) with on return
!             xyke =  material constant matrices.
!	      for each element nel, xyke(1,nel),xyke(2,nel)
!             and xyke(3,nel) represent the constants
!             K11, K22, and K12 at that element.
!
! nodcode   = an integer array containing the boundary information for
!             each node with the following meaning.
!	nodcode(i) = 0 -->  node i is internal
!	nodcode(i) = 1 -->  node i is a boundary but not a corner point
!	nodcode(i) = 2 -->  node i is a corner node. [This node and the
!             corresponmding element are discarded.]
!
! fs	    = real array of length nelx on entry containing the heat
!             source for each element (job = 1 only)
!
! on return
! ---------
! a         = matrix in unassembled form. a(nel,*,*) contains the
!             element matrix for element nel.
!
! f	    = real array containing the right hand for the linears
!             system to solve, in assembled form.
!
! ierr	    = integer. Error message. If (ierr .ne. 0) on return
!             it means that one of the elements has a negative or zero
!             area probably because of a bad ordering of the nodes
!             (see ijk above). Use the subroutine chkelmt to reorder
!             the nodes properly if necessary.
!
!-----------------------------------------------------------------------
!
   Ierr = 0
!
!     assemble the matrix
!
   CALL unassbl(A,Na,F,Nx,Nelx,Ijk,Node,X,Y,Ierr,Xyk)
!
!     if applicable (job .eq. 1) get heat source function
!
   indic = 0
   IF ( Job==1 ) CALL hsourc(indic,Nelx,Node,X,Y,Ijk,Fs,F)
!
!     done
!
END SUBROUTINE genfeu_wbc
!*==genfeu_lstif.f90 processed by SPAG 8.04RA 15:47 21 Feb 2024
!!SPAG Open source Personal, Educational or Academic User Clemson University  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
!----- end of genfeu_wbc -----------------------------------------------
SUBROUTINE genfeu_lstif(Nx,Nelx,Node,Job,X,Y,Ijk,Fs,A,Na,F,Ierr,Xyk)
   USE ISO_FORTRAN_ENV                 
   IMPLICIT NONE
!
! Dummy argument declarations rewritten by SPAG
!
   INTEGER :: Node
   INTEGER :: Na
   INTEGER :: Nx
   INTEGER :: Nelx
   INTEGER , INTENT(IN) :: Job
   REAL(REAL64) , DIMENSION(*) :: X
   REAL(REAL64) , DIMENSION(*) :: Y
   INTEGER , DIMENSION(Node,*) :: Ijk
   REAL(REAL64) , DIMENSION(*) :: Fs
   REAL(REAL64) , DIMENSION(Na,Node,Node) :: A
   REAL(REAL64) , DIMENSION(*) :: F
   INTEGER :: Ierr
   REAL , EXTERNAL :: Xyk
!
! Local variable declarations rewritten by SPAG
!
   INTEGER :: indic
   EXTERNAL funb , func , fung , hsourc , unassbl_lstif
!
! End of declarations rewritten by SPAG
!
!-----------------------------------------------------------------------
! this subroutine generates finite element matrices using unassmbl_lstif.
! The matrix is returned in unassembled form.
! The user must provide the grid, coordinates x, y and connectivity matrix
! ijk) as well as some information on the nodes (nodcode) and the material
! properties (the function K(x,y) above) in the form of a subroutine xyk.
!
!----------------------------------------------------------------------
! moulitsa@cs.umn.edu   : It does not apply boundary conditions
!                         variable nint is eliminated
!----------------------------------------------------------------------
!
! on entry:
! ---------
!
! nx	    = integer . the number of nodes in the grid .
! nelx	    = integer . the number of elements in the grid.
! node      = integer = the number of nodes per element (should be
!             set to three in this version). also the first dimension
!             of ijk
! job	    = integer. If job=0, it is assumed that there is no heat
!             source (i.e. fs = 0) and the right hand side
!             produced will therefore be a zero vector.
!             If job = 1 on entry then the contributions from the
!             heat source in each element are taken into account.
!
! na	    = integer. The first dimension of the array a.
!             a is declared as an array of dimension a(na,node,node).
!
! x, y      = two real arrays containing the coordinates of the nodes.
!
! ijk       =  an integer array containing the connectivity matrix.
!              ijk(i,nel), i=1,2,..node, is the list of the nodes
!              constituting the element nel, ans listed in
!              counter clockwise order.
!
! xyk	    = subroutine defining the material properties at each
!	      element. Form:
! 	      call xyk(nel,xyke,x,y,ijk,node) with on return
!             xyke =  material constant matrices.
!	      for each element nel, xyke(1,nel),xyke(2,nel)
!             and xyke(3,nel) represent the constants
!             K11, K22, and K12 at that element.
!
! nodcode   = an integer array containing the boundary information for
!             each node with the following meaning.
!
! fs	    = real array of length nelx on entry containing the heat
!             source for each element (job = 1 only)
!
! on return
! ---------
! a         = matrix in unassembled form. a(nel,*,*) contains the
!             element matrix for element nel.
!
! f	    = real array containing the right hand for the linears
!             system to solve, in assembled form.
!
! ierr	    = integer. Error message. If (ierr .ne. 0) on return
!             it means that one of the elements has a negative or zero
!             area probably because of a bad ordering of the nodes
!             (see ijk above). Use the subroutine chkelmt to reorder
!             the nodes properly if necessary.
! jwk  = integer work array of length nx
!
!-----------------------------------------------------------------------
!
   Ierr = 0
!
!     assemble the matrix
!
   CALL unassbl_lstif(A,Na,F,Nx,Nelx,Ijk,Node,X,Y,Ierr,Xyk,funb,func,fung)
!
!     if applicable (job .eq. 1) get heat source function
!
   indic = 0
   IF ( Job==1 ) CALL hsourc(indic,Nelx,Node,X,Y,Ijk,Fs,F)
!
!     done
!
END SUBROUTINE genfeu_lstif
!*==assmb1.f90 processed by SPAG 8.04RA 15:47 21 Feb 2024
!!SPAG Open source Personal, Educational or Academic User Clemson University  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
!----- end of genfeu_lstif ---------------------------------------------
!-----------------------------------------------------------------------
SUBROUTINE assmb1(U,Nu,A,Ja,Ia,Fu,F,Nx,Nelx,Ijk,Nodcode,Node,Iwk,Jwk)
!--------------------------------------------------------------
! u	 = unassembled matrix u(na,node,node)
! nu	 = 1-st dimension of u
! a,ja,ia= assembled matrix on output
! fu	 = unassembled right hand side
! f      = right hand side (global load vector) assembled
! nx     = number of nodes at input
! nelx	 = number of elements at input
! ijk	 = connectivity matrix: for node k, ijk(*,k) point to the
!          nodes of element k.
! node	 = total number of nodal points in each element
!
! nodcode= boundary information list for each node with the
!	   following meaning:
!	nodcode(i) = 0 -->  node i is internal
!	nodcode(i) = 1 -->  node i is a boundary but not a corner point
!	nodcode(i) = 2 -->  node i is a corner point (corner points
!
! x,y   = real*8 arrays containing the $x$ and $y$ coordinates
!	  resp. of the nodes.
!         K11, K22, and K12 at that element.
! iwk,jwk = two integer work arrays.
! ierr	= error message integer .
!	  ierr = 0 --> normal return
!	  ierr = 1 --> negative area encountered (due to bad
!	           numbering of nodes of an element- see
!		   message printed in unit iout). not used..
! iout	= output unit (not used here).
!--------------------------------------------------------------
   USE ISO_FORTRAN_ENV                 
   IMPLICIT NONE
!
! Dummy argument declarations rewritten by SPAG
!
   INTEGER , INTENT(IN) :: Nu
   INTEGER , INTENT(IN) :: Node
   REAL(REAL64) , INTENT(IN) , DIMENSION(Nu,Node,Node) :: U
   REAL(REAL64) , INTENT(INOUT) , DIMENSION(*) :: A
   INTEGER , INTENT(INOUT) , DIMENSION(*) :: Ja
   INTEGER , INTENT(INOUT) , DIMENSION(*) :: Ia
   REAL(REAL64) , INTENT(IN) , DIMENSION(Node,*) :: Fu
   REAL(REAL64) , INTENT(INOUT) , DIMENSION(*) :: F
   INTEGER , INTENT(IN) :: Nx
   INTEGER , INTENT(IN) :: Nelx
   INTEGER , INTENT(IN) , DIMENSION(Node,*) :: Ijk
   INTEGER , INTENT(IN) , DIMENSION(*) :: Nodcode
   INTEGER , INTENT(INOUT) , DIMENSION(*) :: Iwk
   INTEGER , INTENT(INOUT) , DIMENSION(*) :: Jwk
!
! Local variable declarations rewritten by SPAG
!
   INTEGER :: i , ii , ilast , irowst , j , jj , k , ka , kb , knod , ksav , ksavn , nel
!
! End of declarations rewritten by SPAG
!
!     max number of nonzeros per row allowed  = 200
!--------------------------------------------------------------
!     initialize
!--------------------------------------------------------------
   DO i = 1 , Nx
      F(i) = 0.0D0
   ENDDO
!
!     initialize  pointer arrays.
!
   DO k = 1 , Nx + 1
      Ia(k) = 1
      Jwk(k) = 0
   ENDDO
   DO k = 1 , Nelx
      DO j = 1 , Node
         knod = Ijk(j,k)
         Ia(knod) = Ia(knod) + 1
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
!     main loop
!-----------------
   DO nel = 1 , Nelx
!
!     get nodal points
!
      DO ka = 1 , Node
         ii = Ijk(ka,nel)
         F(ii) = F(ii) + Fu(ka,nel)
!
!     unpack row into jwk1
!
         irowst = Ia(ii)
         ilast = Iwk(ii)
         DO k = irowst , ilast
            Jwk(Ja(k)) = k
         ENDDO
!
         DO kb = 1 , Node
!
!     column number = jj
!
            jj = Ijk(kb,nel)
            k = Jwk(jj)
            IF ( k==0 ) THEN
               ilast = ilast + 1
               Jwk(jj) = ilast
               Ja(ilast) = jj
               A(ilast) = U(nel,ka,kb)
            ELSE
               A(k) = A(k) + U(nel,ka,kb)
            ENDIF
         ENDDO
!     refresh jwk
         DO k = irowst , ilast
            Jwk(Ja(k)) = 0
         ENDDO
         Iwk(ii) = ilast
      ENDDO
!
   ENDDO
!---------end-of-assmb1----------------------------------------------
END SUBROUTINE assmb1
