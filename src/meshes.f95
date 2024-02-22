!*==inmesh.f90 processed by SPAG 8.04RA 15:47 21 Feb 2024
!!SPAG Open source Personal, Educational or Academic User Clemson University  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
SUBROUTINE inmesh(Nmesh,Iin,Nx,Nelx,Node,X,Y,Nodcode,Ijk,Iperm)
   USE ISO_FORTRAN_ENV                 
   IMPLICIT NONE
!
! Dummy argument declarations rewritten by SPAG
!
   INTEGER , INTENT(INOUT) :: Nx
   INTEGER , INTENT(INOUT) :: Nelx
   INTEGER :: Node
   INTEGER , INTENT(IN) :: Nmesh
   INTEGER , INTENT(IN) :: Iin
   REAL(REAL64) , DIMENSION(*) :: X
   REAL(REAL64) , DIMENSION(*) :: Y
   INTEGER , DIMENSION(Nx) :: Nodcode
   INTEGER , DIMENSION(Node,Nelx) :: Ijk
   INTEGER , DIMENSION(Nx) :: Iperm
!
! Local variable declarations rewritten by SPAG
!
   INTEGER :: i , ii , j
   EXTERNAL fmesh1 , fmesh2 , fmesh3 , fmesh4 , fmesh5 , fmesh6 , fmesh7 , fmesh8 , fmesh9
!
! End of declarations rewritten by SPAG
!
!-----------------------------------------------------------------------
! this subroutine selects and initializes a mesh among a few
! choices. So far there are 9 initial meshes provided and the user can
! also enter his own mesh as a 10th option.
!
! on entry:
!---------
! nmesh	    = integer indicating the mesh chosen. nmesh=1,...,9
!             corresponds to one of the 9 examples supplied by
!             SPARSKIT. nmesh = 0 is a user supplied initial mesh.
!             see below for additional information for the format.
! iin       = integer containing the I/O unit number where to read
!             the data from in case nmesh = 1. A dummy integer
!             otherwise.
! node      = integer = the number of nodes per element (should be
!             set to node=3 in this version). node is also the first
!             dimension of the array ijk.
!
! on return
! ---------
! nx	    = integer . the number of nodes
! nelx	    = integer . the number of elements
! x, y      = two real arrays containing the coordinates of the nodes.
! nodcode   = an integer array containing the boundary information for
!             each node with the following meaning.
!	nodcode(i) = 0 -->  node i is internal
!	nodcode(i) = 1 -->  node i is a boundary but not a corner point
!	nodcode(i) = 2 -->  node i is a corner node.
!
! ijk(node,*)= an integer array containing the connectivity matrix.
!
!-----------------------------------------------------------------------
! format for user supplied mesh (when nmesh = 7)
!
! option nmesh = 0, is a user definied initial mesh.
!---------
! format is as follows:
! line 1: two integers, the first containing the number of nodes
!         the second the number of elements.
! line 2: to line nx+1:  node information.
!        enter the following, one line per node:
!     * the number of the node in the numbering chosen (integer
!          taking the values 1 to nx), followed by,
!        * the coordinates of the nodes (2 reals)  followed by
!	   the boundary information, an integer taking one of the
!          values 0, 1, or 2,  with the meaning explained above.
!
! line nx+2 to nx+nelx+1: connectivity matrix
!       enter the following one line per element:
!       * number of the element in the numbering chosen, followed by
!       * The three numbers of the nodes (according to the above numbering
!       of the nodes) that constitute the element, in a counter clock-wise
!       order (this is in fact not important since it is checked by the
!       subroutine chkelemt).
!
! AN EXAMPLE: consisting on one single element (a triangle)
!------------
!    3    1
!    1    0.0000    0.0000    2
!    2    4.0000    0.0000    2
!    3    0.0000    4.0000    2
!    1    1    2    3
!
!-----------------------------------------------------------------------
!     local variables
!
!     print *, ' ----- nmesh = ', nmesh
   IF ( Nmesh+1==1 ) THEN
!
!-------option 0 : reading mesh from IO unit iin.
!
      READ (Iin,*) Nx , Nelx
!
      DO i = 1 , Nx
         READ (Iin,*) ii , X(ii) , Y(ii) , Nodcode(ii)
      ENDDO
      DO i = 1 , Nelx
         READ (Iin,*) ii , (Ijk(j,ii),j=1,Node)
         IF ( ii>Nelx ) Nelx = ii
      ENDDO
   ELSEIF ( Nmesh+1==3 ) THEN
      CALL fmesh2(Nx,Nelx,Node,X,Y,Nodcode,Ijk)
   ELSEIF ( Nmesh+1==4 ) THEN
      CALL fmesh3(Nx,Nelx,Node,X,Y,Nodcode,Ijk)
   ELSEIF ( Nmesh+1==5 ) THEN
      CALL fmesh4(Nx,Nelx,Node,X,Y,Nodcode,Ijk)
   ELSEIF ( Nmesh+1==6 ) THEN
      CALL fmesh5(Nx,Nelx,Node,X,Y,Nodcode,Ijk)
   ELSEIF ( Nmesh+1==7 ) THEN
      CALL fmesh6(Nx,Nelx,Node,X,Y,Nodcode,Ijk)
   ELSEIF ( Nmesh+1==8 ) THEN
      CALL fmesh7(Nx,Nelx,Node,X,Y,Nodcode,Ijk,Iperm)
   ELSEIF ( Nmesh+1==9 ) THEN
      CALL fmesh8(Nx,Nelx,Node,X,Y,Nodcode,Ijk,Iperm)
   ELSEIF ( Nmesh+1==10 ) THEN
      CALL fmesh9(Nx,Nelx,Node,X,Y,Nodcode,Ijk,Iperm)
   ELSE
      CALL fmesh1(Nx,Nelx,Node,X,Y,Nodcode,Ijk)
   ENDIF
!-----------------------------------------------------------------------
!     and return
!-----------------------------------------------------------------------
END SUBROUTINE inmesh
!*==fmesh1.f90 processed by SPAG 8.04RA 15:47 21 Feb 2024
!!SPAG Open source Personal, Educational or Academic User Clemson University  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
!-----------------------------------------------------------------------
SUBROUTINE fmesh1(Nx,Nelx,Node,X,Y,Nodcode,Ijk)
!--------------------------------------------------------------
!
! initial mesh for a simple square with two elemnts
!      3             4
!       --------------
!       |          . |
!       |   2    .   |
!       |      .     |
!       |   .    1   |
!       | .          |
!       --------------
!      1              2
!--------------------------------------------------------------
! input parameters: node = first dimensoin of ijk (must be .ge. 3)
! output parameters:
!    nx    = number of nodes
!    nelx = number of elemnts
!    (x(1:nx), y(1:nx)) = coordinates of nodes
!    nodcode(1:nx) = integer code for each node with the
!	    following meening:
!	nodcode(i) = 0 -->  node i is internal
!	nodcode(i) = 1 -->  node i is a boundary but not a corner point
!	nodcode(i) = 2 -->  node i is a corner point.
!   ijk(1:3,1:nelx) = connectivity matrix. for a given element
!	    number nel, ijk(k,nel), k=1,2,3 represent the nodes
!	    composing the element nel.
!--------------------------------------------------------------
   USE ISO_FORTRAN_ENV                 
   IMPLICIT NONE
!
! Dummy argument declarations rewritten by SPAG
!
   INTEGER , INTENT(IN) :: Node
   INTEGER , INTENT(INOUT) :: Nx
   INTEGER , INTENT(INOUT) :: Nelx
   REAL(REAL64) , INTENT(OUT) , DIMENSION(*) :: X
   REAL(REAL64) , INTENT(OUT) , DIMENSION(*) :: Y
   INTEGER , INTENT(OUT) , DIMENSION(*) :: Nodcode
   INTEGER , INTENT(OUT) , DIMENSION(Node,*) :: Ijk
!
! Local variable declarations rewritten by SPAG
!
   INTEGER , DIMENSION(2) , SAVE :: ijk1 , ijk2 , ijk3
   INTEGER :: k
   REAL(REAL64) , DIMENSION(4) , SAVE :: x1 , y1
!
! End of declarations rewritten by SPAG
!
!--------------------------------------------------------------
! coordinates of nodal points
!--------------------------------------------------------------
   DATA x1/0.0 , 1.0 , 0.0 , 1.0/
   DATA y1/0.0 , 0.0 , 1.0 , 1.0/
!
!------------------|--|
! elements         1  2
!------------------|--|
   DATA ijk1/1 , 1/
   DATA ijk2/2 , 4/
   DATA ijk3/4 , 3/
!
   Nx = 4
!
   DO k = 1 , Nx
      X(k) = x1(k)
      Y(k) = y1(k)
      Nodcode(k) = 1
   ENDDO
!
   Nodcode(2) = 2
   Nodcode(3) = 2
!
   Nelx = 2
!
   DO k = 1 , Nelx
      Ijk(1,k) = ijk1(k)
      Ijk(2,k) = ijk2(k)
      Ijk(3,k) = ijk3(k)
   ENDDO
!
END SUBROUTINE fmesh1
!*==fmesh2.f90 processed by SPAG 8.04RA 15:47 21 Feb 2024
!!SPAG Open source Personal, Educational or Academic User Clemson University  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
!-----------------------------------------------------------------------
SUBROUTINE fmesh2(Nx,Nelx,Node,X,Y,Nodcode,Ijk)
!---------------------------------------------------------------
! initial mesh for a simple D-shaped region with 4 elemnts
!       6
!       | .
!       |    .
!       |      .
!       |   4     .
!       |           .
!     4 -------------- 5
!       |          . |
!       |   3    .   |
!       |      .     |
!       |   .    2   |
!       | .          |
!       --------------
!       | 2         . 3
!       |         .
!       |   1   .
!       |     .
!       |  .
!       |.
!       1
!--------------------------------------------------------------
! input parameters: node = first dimensoin of ijk (must be .ge. 3)
! output parameters:
!    nx    = number of nodes
!    nelx = number of elemnts
!    (x(1:nx), y(1:nx)) = coordinates of nodes
!    nodcode(1:nx) = integer code for each node with the
!	    following meening:
!	nodcode(i) = 0 -->  node i is internal
!	nodcode(i) = 1 -->  node i is a boundary but not a corner point
!	nodcode(i) = 2 -->  node i is a corner point.
!   ijk(1:3,1:nelx) = connectivity matrix. for a given element
!	    number nel, ijk(k,nel), k=1,2,3 represent the nodes
!	    composing the element nel.
!
!--------------------------------------------------------------
   USE ISO_FORTRAN_ENV                 
   IMPLICIT NONE
!
! Dummy argument declarations rewritten by SPAG
!
   INTEGER , INTENT(IN) :: Node
   INTEGER , INTENT(INOUT) :: Nx
   INTEGER , INTENT(INOUT) :: Nelx
   REAL(REAL64) , INTENT(OUT) , DIMENSION(*) :: X
   REAL(REAL64) , INTENT(OUT) , DIMENSION(*) :: Y
   INTEGER , INTENT(OUT) , DIMENSION(*) :: Nodcode
   INTEGER , INTENT(OUT) , DIMENSION(Node,*) :: Ijk
!
! Local variable declarations rewritten by SPAG
!
   INTEGER , DIMENSION(4) , SAVE :: ijk1 , ijk2 , ijk3
   INTEGER :: k
   REAL(REAL64) , DIMENSION(6) , SAVE :: x1 , y1
!
! End of declarations rewritten by SPAG
!
!--------------------------------------------------------------
! coordinates of nodal points
!--------------------------------------------------------------
   DATA x1/0.0 , 0.0 , 1.0 , 0.0 , 1.0 , 0.0/
   DATA y1/0.0 , 1.0 , 1.0 , 2.0 , 2.0 , 3.0/
!
!------------------|--|--|--|
! elements         1  2  3  4
!------------------|--|--|--|
   DATA ijk1/1 , 2 , 2 , 4/
   DATA ijk2/3 , 3 , 5 , 5/
   DATA ijk3/2 , 5 , 4 , 6/
!
   Nx = 6
!
   DO k = 1 , Nx
      X(k) = x1(k)
      Y(k) = y1(k)
      Nodcode(k) = 1
   ENDDO
!
   Nelx = 4
!
   DO k = 1 , Nelx
      Ijk(1,k) = ijk1(k)
      Ijk(2,k) = ijk2(k)
      Ijk(3,k) = ijk3(k)
   ENDDO
!
END SUBROUTINE fmesh2
!*==fmesh3.f90 processed by SPAG 8.04RA 15:47 21 Feb 2024
!!SPAG Open source Personal, Educational or Academic User Clemson University  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
!-----------------------------------------------------------------------
SUBROUTINE fmesh3(Nx,Nelx,Node,X,Y,Nodcode,Ijk)
!---------------------------------------------------------------
! initial mesh for a C-shaped region composed of 10 elements --
!
!
!      10           11            12
!       ---------------------------
!       |          . |          . |
!       |  7     .   |   9    .   |
!       |      .     |      .     |
!       |   .    8   |   .   10   |
!       | .          | .          |
!     7 ---------------------------
!       |          . |8           9
!       |   5    .   |
!       |      .     |
!       |   .    6   |
!     4 | .          |5           6
!       ---------------------------
!       |          . |          . |
!       |   1    .   |  3     .   |
!       |      .     |      .     |
!       |   .    2   |   .   4    |
!       | .          | .          |
!       ---------------------------
!      1             2            3
!
!--------------------------------------------------------------
! input parameters: node = first dimensoin of ijk (must be .ge. 3)
!    nx    = number of nodes
!    nelx = number of elemnts
!    (x(1:nx), y(1:nx)) = coordinates of nodes
!    nodcode(1:nx) = integer code for each node with the
!	    following meening:
!	nodcode(i) = 0 -->  node i is internal
!	nodcode(i) = 1 -->  node i is a boundary but not a corner point
!	nodcode(i) = 2 -->  node i is a corner point.
!   ijk(1:3,1:nelx) = connectivity matrix. for a given element
!	    number nel, ijk(k,nel), k=1,2,3 represent the nodes
!	    composing the element nel.
!
!--------------------------------------------------------------
   USE ISO_FORTRAN_ENV                 
   IMPLICIT NONE
!
! Dummy argument declarations rewritten by SPAG
!
   INTEGER , INTENT(IN) :: Node
   INTEGER , INTENT(INOUT) :: Nx
   INTEGER , INTENT(INOUT) :: Nelx
   REAL(REAL64) , INTENT(OUT) , DIMENSION(*) :: X
   REAL(REAL64) , INTENT(OUT) , DIMENSION(*) :: Y
   INTEGER , INTENT(OUT) , DIMENSION(*) :: Nodcode
   INTEGER , INTENT(OUT) , DIMENSION(Node,*) :: Ijk
!
! Local variable declarations rewritten by SPAG
!
   INTEGER , DIMENSION(10) , SAVE :: ijk1 , ijk2 , ijk3
   INTEGER :: k
   REAL(REAL64) , DIMENSION(12) , SAVE :: x1 , y1
!
! End of declarations rewritten by SPAG
!
!--------------------------------------------------------------
! coordinates of nodal points
!--------------------------------------------------------------
   DATA x1/0.0 , 1.0 , 2.0 , 0.0 , 1.0 , 2.0 , 0.0 , 1.0 , 2.0 , 0.0 , 1.0 , 2.0/
   DATA y1/0.0 , 0.0 , 0.0 , 1.0 , 1.0 , 1.0 , 2.0 , 2.0 , 2.0 , 3.0 , 3.0 , 3.0/
!
!------------------|--|--|--|--|--|--|---|---|---|
! elements         1  2  3  4  5  6  7   8   9  10
!------------------|--|--|--|--|--|--|---|---|---|
   DATA ijk1/1 , 1 , 2 , 2 , 4 , 4 , 7 , 7 , 8 , 8/
   DATA ijk2/5 , 2 , 6 , 3 , 8 , 5 , 11 , 8 , 12 , 9/
   DATA ijk3/4 , 5 , 5 , 6 , 7 , 8 , 10 , 11 , 11 , 12/
!
   Nx = 12
!
   DO k = 1 , Nx
      X(k) = x1(k)
      Y(k) = y1(k)
      Nodcode(k) = 1
   ENDDO
!
   Nodcode(3) = 2
   Nodcode(10) = 2
   Nodcode(9) = 2
!
   Nelx = 10
!
   DO k = 1 , Nelx
      Ijk(1,k) = ijk1(k)
      Ijk(2,k) = ijk2(k)
      Ijk(3,k) = ijk3(k)
   ENDDO
!
END SUBROUTINE fmesh3
!*==fmesh4.f90 processed by SPAG 8.04RA 15:47 21 Feb 2024
!!SPAG Open source Personal, Educational or Academic User Clemson University  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
!-----------------------------------------------------------------------
SUBROUTINE fmesh4(Nx,Nelx,Node,X,Y,Nodcode,Ijk)
!-----------------------------------------------------------------------
! initial mesh for a C-shaped region composed of 10 elements --
!      10                   11
!       +------------------+ .
!       | .                |    .
!       |    .       8     |       . 12
!       |        .         |  9   . |
!       |     7      .     |   .    |
!     7 |                . | .   10 |
!       -------------------+--------+ 9
!       |                 .| 8
!       |     5       .    |
!       |         .        |
!       |    .       6     |
!       |.                 | 5      6
!    4  +------------------+--------+
!       |               .  | .   4  |
!       |    1       .     |    .   |
!       |        .         |  3    .| 3
!       |    .        2    |    .
!       | .                | .
!       --------------------
!       1                  2
!--------------------------------------------------------------
! input parameters: node = first dimensoin of ijk (must be .ge. 3)
!    nx    = number of nodes
!    nelx = number of elemnts
!    (x(1:nx), y(1:nx)) = coordinates of nodes
!    nodcode(1:nx) = integer code for each node with the
!	    following meening:
!	nodcode(i) = 0 -->  node i is internal
!	nodcode(i) = 1 -->  node i is a boundary but not a corner point
!	nodcode(i) = 2 -->  node i is a corner point.
!   ijk(1:3,1:nelx) = connectivity matrix. for a given element
!	    number nel, ijk(k,nel), k=1,2,3 represent the nodes
!	    composing the element nel.
!
!--------------------------------------------------------------
   USE ISO_FORTRAN_ENV                 
   IMPLICIT NONE
!
! Dummy argument declarations rewritten by SPAG
!
   INTEGER , INTENT(IN) :: Node
   INTEGER , INTENT(INOUT) :: Nx
   INTEGER , INTENT(INOUT) :: Nelx
   REAL(REAL64) , INTENT(OUT) , DIMENSION(*) :: X
   REAL(REAL64) , INTENT(OUT) , DIMENSION(*) :: Y
   INTEGER , INTENT(OUT) , DIMENSION(*) :: Nodcode
   INTEGER , INTENT(OUT) , DIMENSION(Node,*) :: Ijk
!
! Local variable declarations rewritten by SPAG
!
   INTEGER , DIMENSION(10) , SAVE :: ijk1 , ijk2 , ijk3
   INTEGER :: k
   REAL(REAL64) , DIMENSION(12) , SAVE :: x1 , y1
!
! End of declarations rewritten by SPAG
!
!--------------------------------------------------------------
! coordinates of nodal points
!--------------------------------------------------------------
   DATA x1/0.0 , 1.0 , 1.5 , 0.0 , 1.0 , 1.5 , 0.0 , 1.0 , 1.5 , 0.0 , 1.0 , 1.5/
   DATA y1/0.0 , 0.0 , 0.5 , 1.0 , 1.0 , 1.0 , 2.0 , 2.0 , 2.0 , 3.0 , 3.0 , 2.5/
!
!------------------|--|--|--|--|--|--|---|---|---|
! elements         1  2  3  4  5  6  7   8   9  10
!------------------|--|--|--|--|--|--|---|---|---|
   DATA ijk1/1 , 1 , 2 , 5 , 4 , 4 , 7 , 10 , 8 , 8/
   DATA ijk2/5 , 2 , 3 , 3 , 8 , 5 , 8 , 8 , 12 , 9/
   DATA ijk3/4 , 5 , 5 , 6 , 7 , 8 , 10 , 11 , 11 , 12/
!
   Nx = 12
!
   DO k = 1 , Nx
      X(k) = x1(k)
      Y(k) = y1(k)
      Nodcode(k) = 1
   ENDDO
!
   Nodcode(6) = 2
   Nodcode(9) = 2
!
   Nelx = 10
!
   DO k = 1 , Nelx
      Ijk(1,k) = ijk1(k)
      Ijk(2,k) = ijk2(k)
      Ijk(3,k) = ijk3(k)
   ENDDO
!
END SUBROUTINE fmesh4
!*==fmesh5.f90 processed by SPAG 8.04RA 15:47 21 Feb 2024
!!SPAG Open source Personal, Educational or Academic User Clemson University  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
!-----------------------------------------------------------------------
SUBROUTINE fmesh5(Nx,Nelx,Node,X,Y,Nodcode,Ijk)
!---------------------------------------------------------------
!     initial mesh for a whrench shaped region composed of 14 elements --
!
!                                      13            15
!                                        . ----------.           |-3
!                                      .   .   13  .   .         |
!                                   .   12   .   .  14    .      |
! 9        10        11       12  .            . 14        . 16  |
! ----------------------------------------------------------     |-2
! |       . |       . |       . |            . |                 |
! | 1   .   |  3  .   |  5  .   |    7   .     |                 |
! |   .  2  |   .  4  |   .  6  |     .    8   |                 |
! |.        |.        |.        | .            |                 |
! -----------------------------------------------------------    |-1
! 1         2         3       4  .           6 .           . 8   |
!                                   .   9    .   .   11   .      |
!                                      .   .  10    .   .        |
!                                        .___________.           |-0
!                                       5             7
!
! 0---------1--------2----------3--------------4-------------5
!--------------------------------------------------------------
! input parameters: node = first dimensoin of ijk (must be .ge. 3)
!    nx    = number of nodes
!    nelx = number of elemnts
!    (x(1:nx), y(1:nx)) = coordinates of nodes
!    nodcode(1:nx) = integer code for each node with the
!	    following meening:
!	nodcode(i) = 0 -->  node i is internal
!	nodcode(i) = 1 -->  node i is a boundary but not a corner point
!	nodcode(i) = 2 -->  node i is a corner point.
!   ijk(1:3,1:nelx) = connectivity matrix. for a given element
!	    number nel, ijk(k,nel), k=1,2,3 represent the nodes
!	    composing the element nel.
!
!--------------------------------------------------------------
   USE ISO_FORTRAN_ENV                 
   IMPLICIT NONE
!
! Dummy argument declarations rewritten by SPAG
!
   INTEGER , INTENT(IN) :: Node
   INTEGER , INTENT(INOUT) :: Nx
   INTEGER , INTENT(INOUT) :: Nelx
   REAL(REAL64) , INTENT(OUT) , DIMENSION(*) :: X
   REAL(REAL64) , INTENT(OUT) , DIMENSION(*) :: Y
   INTEGER , INTENT(OUT) , DIMENSION(*) :: Nodcode
   INTEGER , INTENT(OUT) , DIMENSION(Node,*) :: Ijk
!
! Local variable declarations rewritten by SPAG
!
   INTEGER , DIMENSION(14) , SAVE :: ijk1 , ijk2 , ijk3
   INTEGER :: k
   REAL(REAL64) , DIMENSION(16) , SAVE :: x1 , y1
!
! End of declarations rewritten by SPAG
!
!--------------------------------------------------------------
!     coordinates of nodal points
!--------------------------------------------------------------
   DATA x1/0. , 1. , 2. , 3. , 3.5 , 4. , 4.5 , 5. , 0. , 1. , 2. , 3. , 3.5 , 4. , 4.5 , 5./
   DATA y1/1. , 1. , 1. , 1. , 0. , 1. , 0. , 1. , 2. , 2. , 2. , 2. , 3. , 2. , 3. , 2./
!
!------------------|--|--|--|--|--|--|---|---|---|--|---|---|---|
! elements         1  2  3  4  5  6  7   8   9  10  11  12  13  14
!------------------|--|--|--|--|--|--|---|---|---|--|---|---|---|
   DATA ijk1/1 , 1 , 2 , 2 , 3 , 3 , 4 , 4 , 4 , 5 , 6 , 12 , 14 , 14/
   DATA ijk2/10 , 2 , 11 , 3 , 12 , 4 , 14 , 6 , 5 , 7 , 7 , 14 , 15 , 16/
   DATA ijk3/9 , 10 , 10 , 11 , 11 , 12 , 12 , 14 , 6 , 6 , 8 , 13 , 13 , 15/
!
   Nx = 16
!
   DO k = 1 , Nx
      X(k) = x1(k)
      Y(k) = y1(k)
      Nodcode(k) = 1
   ENDDO
!
   Nodcode(9) = 2
   Nodcode(8) = 2
   Nodcode(16) = 2
!
   Nelx = 14
!
   DO k = 1 , Nelx
      Ijk(1,k) = ijk1(k)
      Ijk(2,k) = ijk2(k)
      Ijk(3,k) = ijk3(k)
   ENDDO
!
END SUBROUTINE fmesh5
!*==fmesh6.f90 processed by SPAG 8.04RA 15:47 21 Feb 2024
!!SPAG Open source Personal, Educational or Academic User Clemson University  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
!-----------------------------------------------------------------------
SUBROUTINE fmesh6(Nx,Nelx,Node,X,Y,Nodcode,Ijk)
!---------------------------------------------------------------
! this generates a finite element mesh for an ellipse-shaped
! domain.
!---------------------------------------------------------------
   USE ISO_FORTRAN_ENV                 
   IMPLICIT NONE
!
! Dummy argument declarations rewritten by SPAG
!
   INTEGER :: Node
   INTEGER , INTENT(INOUT) :: Nx
   INTEGER :: Nelx
   REAL(REAL64) , DIMENSION(*) :: X
   REAL(REAL64) , DIMENSION(*) :: Y
   INTEGER , INTENT(OUT) , DIMENSION(*) :: Nodcode
   INTEGER , DIMENSION(Node,*) :: Ijk
!
! Local variable declarations rewritten by SPAG
!
   REAL(REAL64) :: a , ar , b , br , delr , pi , theta
   INTEGER :: i , j , k , nd , nemax , nr
   INTEGER , DIMENSION(200,3) :: ijktr
   INTEGER , DIMENSION(200) :: nel
   EXTERNAL chkelmt , dlauny
!
! End of declarations rewritten by SPAG
!
!--------------------------------------------------------------
!     coordinates of nodal points
!--------------------------------------------------------------
   nd = 8
   nr = 3
!
!     define axes of ellipse
!
   a = 2.0
   b = 1.30
!
   Nx = 1
   pi = 4.0*atan(1.0)
   theta = 2.0*pi/real(nd)
   X(1) = 0.0
   Y(1) = 0.0
   delr = a/real(nr)
   Nx = 0
   DO i = 1 , nr
      ar = real(i)*delr
      br = ar*b/a
      DO j = 1 , nd
         Nx = Nx + 1
         X(Nx) = a + ar*cos(real(j)*theta)
         Y(Nx) = b + br*sin(real(j)*theta)
!            write (13,*) ' nod ', nx, ' x,y', x(nx), y(nx)
         Nodcode(Nx) = 0
         IF ( i==nr ) Nodcode(Nx) = 1
      ENDDO
   ENDDO
!
   nemax = 200
   CALL dlauny(X,Y,Nx,ijktr,nemax,Nelx)
!
!     print *, ' delauny -- nx, nelx ', nx, nelx
   DO j = 1 , Nx
      nel(j) = 0
   ENDDO
!     transpose ijktr into ijk and count the number of
!     elemnts to which each node belongs
!
   DO j = 1 , Nelx
      DO k = 1 , Node
         i = ijktr(j,k)
         Ijk(k,j) = i
         nel(i) = nel(i) + 1
      ENDDO
   ENDDO
!
!     take care of ordering within each element
!
   CALL chkelmt(X,Y,Nelx,Ijk,Node)
!
END SUBROUTINE fmesh6
!*==fmesh7.f90 processed by SPAG 8.04RA 15:47 21 Feb 2024
!!SPAG Open source Personal, Educational or Academic User Clemson University  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
!--------------------------------------------------------
SUBROUTINE fmesh7(Nx,Nelx,Node,X,Y,Nodcode,Ijk,Iperm)
   USE ISO_FORTRAN_ENV                 
   IMPLICIT NONE
!
! Dummy argument declarations rewritten by SPAG
!
   INTEGER , INTENT(INOUT) :: Nx
   INTEGER :: Nelx
   INTEGER :: Node
   REAL(REAL64) , INTENT(INOUT) , DIMENSION(*) :: X
   REAL(REAL64) , INTENT(INOUT) , DIMENSION(*) :: Y
   INTEGER , INTENT(INOUT) , DIMENSION(Nx) :: Nodcode
   INTEGER , DIMENSION(Node,Nelx) :: Ijk
   INTEGER , DIMENSION(Nx) :: Iperm
!
! Local variable declarations rewritten by SPAG
!
   REAL(REAL64) :: a , arx , ary , b , delr , excl , pi , rad , theta , xcntr , xnew , ycntr , ynew
   INTEGER :: i , j , k , nemax , nr , nsec
   INTEGER , DIMENSION(200,3) :: ijktr
   INTEGER , DIMENSION(200) :: nel
   INTEGER , SAVE :: nodexc
   REAL(REAL64) , SAVE :: x1 , x2 , y1 , y2
   EXTERNAL chkelmt , cleanel , cleannods , clos2bdr , dlauny
!
! End of declarations rewritten by SPAG
!
!---------------------------------------------------------------
!     this generates a U-shaped domain with an elliptic inside.
!     then a Delauney triangulation is used to generate the mesh.
!     mesh needs to be post-processed -- see inmesh --
!---------------------------------------------------------------
!--------------------------------------------------------------
!     coordinates of nodal points
!--------------------------------------------------------------
   DATA x1/0.0/ , y1/0.0/ , x2/6.0/ , y2/6.0/ , nodexc/1/
   xcntr = x2/3.0
   ycntr = y2/2.0
   rad = 1.8
   nsec = 20
   nr = 3
!
!     exclusion zone near the boundary
!
   excl = 0.02*x2
!-----------------------------------------------------------------------
!     enter the four corner points.
!-----------------------------------------------------------------------
   Nx = 1
   X(Nx) = x1
   Y(Nx) = y1
   Nodcode(Nx) = 1
   Nx = Nx + 1
   X(Nx) = x2
   Y(Nx) = y1
   Nodcode(Nx) = 1
   Nx = Nx + 1
   X(Nx) = x2
   Y(Nx) = y2
   Nodcode(Nx) = 1
   Nx = Nx + 1
   X(Nx) = x1
   Y(Nx) = y2
   Nodcode(Nx) = 1
!
!     define axes of ellipse
!
   a = 2.0
   b = 1.30
!-----------------------------------------------------------------------
   pi = 4.0*atan(1.0)
   delr = a/real(nr)
   DO i = 1 , nsec
      theta = 2.0*real(i-1)*pi/real(nsec)
      xnew = xcntr + rad*cos(theta)
      ynew = ycntr + rad*b*sin(theta)/a
      IF ( (xnew<x2) .AND. (xnew>x1) .AND. (ynew<y2) .AND. (ynew>y1) ) THEN
         Nx = Nx + 1
         X(Nx) = xnew
         Y(Nx) = ynew
         Nodcode(Nx) = nodexc
         arx = delr*cos(theta)
         ary = delr*b*sin(theta)/a
         SPAG_Loop_2_1: DO
!
!     while inside domain do:
!
            xnew = X(Nx) + arx
            ynew = Y(Nx) + ary
            IF ( xnew>=x2 ) THEN
               X(Nx) = x2
               Nodcode(Nx) = 1
            ELSEIF ( xnew<=x1 ) THEN
               X(Nx) = x1
               Nodcode(Nx) = 1
            ELSEIF ( ynew>=y2 ) THEN
               Y(Nx) = y2
               Nodcode(Nx) = 1
            ELSEIF ( ynew<=y1 ) THEN
               Y(Nx) = y1
               Nodcode(Nx) = 1
            ELSE
               Nx = Nx + 1
               X(Nx) = xnew
               Y(Nx) = ynew
               Nodcode(Nx) = 0
               CALL clos2bdr(Nx,xnew,ynew,X,Y,x1,x2,y1,y2,excl,Nodcode)
            ENDIF
!        write (13,*) ' nod ', nx, ' x,y', x(nx), y(nx)
!     *         ,' arx--ary ', arx, ary
            arx = arx*1.2
            ary = ary*1.2
            IF ( Nodcode(Nx)>0 ) EXIT SPAG_Loop_2_1
         ENDDO SPAG_Loop_2_1
      ENDIF
   ENDDO
!
   nemax = 200
   CALL dlauny(X,Y,Nx,ijktr,nemax,Nelx)
!
!     print *, ' delauny -- nx, nelx ', nx, nelx
   DO j = 1 , Nx
      nel(j) = 0
   ENDDO
!
!     transpose ijktr into ijk and count the number of
!     elemnts to which each node belongs
!
   DO j = 1 , Nelx
      DO k = 1 , Node
         i = ijktr(j,k)
         Ijk(k,j) = i
         nel(i) = nel(i) + 1
      ENDDO
   ENDDO
!
!     this mesh needs cleaning up --
!
   CALL cleanel(Nelx,Ijk,Node,Nodcode,nodexc)
   CALL cleannods(Nx,X,Y,Nelx,Ijk,Node,Nodcode,Iperm)
!
!     take care of ordering within each element
!
   CALL chkelmt(X,Y,Nelx,Ijk,Node)
END SUBROUTINE fmesh7
!*==fmesh8.f90 processed by SPAG 8.04RA 15:47 21 Feb 2024
!!SPAG Open source Personal, Educational or Academic User Clemson University  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
!-----------------------------------------------------------------------
SUBROUTINE fmesh8(Nx,Nelx,Node,X,Y,Nodcode,Ijk,Iperm)
   USE ISO_FORTRAN_ENV                 
   IMPLICIT NONE
!
! Dummy argument declarations rewritten by SPAG
!
   INTEGER , INTENT(INOUT) :: Nx
   INTEGER :: Nelx
   INTEGER :: Node
   REAL(REAL64) , DIMENSION(*) :: X
   REAL(REAL64) , DIMENSION(*) :: Y
   INTEGER , INTENT(INOUT) , DIMENSION(Nx) :: Nodcode
   INTEGER , DIMENSION(Node,Nelx) :: Ijk
   INTEGER , DIMENSION(Nx) :: Iperm
!
! Local variable declarations rewritten by SPAG
!
   REAL(REAL64) :: a , arx , ary , b , excl , pi , rad , radi , theta , xcntr , xnew , ycntr , ynew
   INTEGER :: i , j , k , nemax , nr , nsec
   INTEGER , DIMENSION(1500,3) :: ijktr
   INTEGER , DIMENSION(1500) :: nel
   INTEGER , SAVE :: nodexc
   REAL(REAL64) , SAVE :: x1 , x2 , y1 , y2
   EXTERNAL chkelmt , cleanel , cleannods , clos2bdr , dlauny
!
! End of declarations rewritten by SPAG
!
!---------------------------------------------------------------
!     this generates a small rocket type shape inside a rectangle
!     then a Delauney triangulation is used to generate the mesh.
!     mesh needs to be post-processed -- see inmesh --
!---------------------------------------------------------------
!--------------------------------------------------------------
!     coordinates of corners + some additional data
!--------------------------------------------------------------
   DATA x1/0.0/ , y1/0.0/ , x2/6.0/ , y2/6.0/ , nodexc/3/
   xcntr = 4.0
   ycntr = y2/2.0
   rad = 0.6
!
!     exclusion zone near the boundary.
!
   excl = 0.02*x2
   nsec = 30
   nr = 4
!-----------------------------------------------------------------------
!     enter the four corner points.
!-----------------------------------------------------------------------
   Nx = 1
   X(Nx) = x1
   Y(Nx) = y1
   Nodcode(Nx) = 1
   Nx = Nx + 1
   X(Nx) = x2
   Y(Nx) = y1
   Nodcode(Nx) = 1
   Nx = Nx + 1
   X(Nx) = x2
   Y(Nx) = y2
   Nodcode(Nx) = 1
   Nx = Nx + 1
   X(Nx) = x1
   Y(Nx) = y2
   Nodcode(Nx) = 1
!
!     define axes of ellipse /circle / object
!
   a = 2.0
   b = 1.0
!-----------------------------------------------------------------------
   pi = 4.0*atan(1.0)
   DO i = 1 , nsec
      theta = 2.0*real(i-1)*pi/real(nsec)
      IF ( theta>pi ) theta = theta - 2.0*pi
      radi = rad*(1.0+0.05*((pi/2.0)**2-theta**2)**2)/(1.0+0.05*(pi/2.0)**4)
      arx = radi*cos(theta)/real(nr)
      ary = radi*sin(theta)/real(nr)
!      a hack!
!         arx = (abs(theta)+0.25)*cos(theta)/real(nr)
!         ary = (abs(theta)+0.25)*sin(theta)/real(nr)
!
      xnew = xcntr + radi*cos(theta)
      ynew = ycntr + radi*b*sin(theta)/a
      IF ( (xnew<x2) .AND. (xnew>x1) .AND. (ynew<y2) .AND. (ynew>y1) ) THEN
         Nx = Nx + 1
         X(Nx) = xnew
         Y(Nx) = ynew
         Nodcode(Nx) = nodexc
         SPAG_Loop_2_1: DO
!
!     while inside domain do:
!
            xnew = xnew + arx
            ynew = ynew + ary
            IF ( xnew>=x2 ) THEN
               X(Nx) = x2
               Nodcode(Nx) = 1
            ELSEIF ( xnew<=x1 ) THEN
               X(Nx) = x1
               Nodcode(Nx) = 1
            ELSEIF ( ynew>=y2 ) THEN
               Y(Nx) = y2
               Nodcode(Nx) = 1
            ELSEIF ( ynew<=y1 ) THEN
               Y(Nx) = y1
               Nodcode(Nx) = 1
!
!     else we can add this as interior point
!
            ELSE
               Nx = Nx + 1
               X(Nx) = xnew
               Y(Nx) = ynew
               Nodcode(Nx) = 0
!
!     do something if point is too close to boundary
!
               CALL clos2bdr(Nx,xnew,ynew,X,Y,x1,x2,y1,y2,excl,Nodcode)
            ENDIF
!
            arx = arx*1.1
            ary = ary*1.1
            IF ( Nodcode(Nx)/=0 ) EXIT SPAG_Loop_2_1
         ENDDO SPAG_Loop_2_1
      ENDIF
   ENDDO
!
   nemax = 1500
   CALL dlauny(X,Y,Nx,ijktr,nemax,Nelx)
!
   PRINT * , ' delauney -- nx, nelx ' , Nx , Nelx
   DO j = 1 , Nx
      nel(j) = 0
   ENDDO
!-----------------------------------------------------------------------
!     transpose ijktr into ijk and count the number of
!     elemnts to which each node belongs
!-----------------------------------------------------------------------
   DO j = 1 , Nelx
      DO k = 1 , Node
         i = ijktr(j,k)
         Ijk(k,j) = i
         nel(i) = nel(i) + 1
      ENDDO
   ENDDO
!
!     this mesh needs cleaning up --
!
   CALL cleanel(Nelx,Ijk,Node,Nodcode,nodexc)
   CALL cleannods(Nx,X,Y,Nelx,Ijk,Node,Nodcode,Iperm)
!
!     take care of ordering within each element
!
   CALL chkelmt(X,Y,Nelx,Ijk,Node)
END SUBROUTINE fmesh8
!*==fmesh9.f90 processed by SPAG 8.04RA 15:47 21 Feb 2024
!!SPAG Open source Personal, Educational or Academic User Clemson University  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
!-----------------------------------------------------------------------
SUBROUTINE fmesh9(Nx,Nelx,Node,X,Y,Nodcode,Ijk,Iperm)
   USE ISO_FORTRAN_ENV                 
   IMPLICIT NONE
!
! Dummy argument declarations rewritten by SPAG
!
   INTEGER , INTENT(INOUT) :: Nx
   INTEGER :: Nelx
   INTEGER :: Node
   REAL(REAL64) , INTENT(INOUT) , DIMENSION(*) :: X
   REAL(REAL64) , INTENT(INOUT) , DIMENSION(*) :: Y
   INTEGER , INTENT(INOUT) , DIMENSION(Nx) :: Nodcode
   INTEGER , DIMENSION(Node,Nelx) :: Ijk
   INTEGER , DIMENSION(Nx) :: Iperm
!
! Local variable declarations rewritten by SPAG
!
   REAL(REAL64) :: arx , ary , delr , excl , pi , rad , theta , xcntr , xnew , ycntr , ynew
   INTEGER :: i , j , k , nemax , nr , nsec
   INTEGER , DIMENSION(1500,3) :: ijktr
   INTEGER , DIMENSION(1500) :: nel
   INTEGER , SAVE :: nodexc
   REAL(REAL64) , SAVE :: x1 , x2 , y1 , y2
   EXTERNAL chkelmt , cleanel , cleannods , clos2bdr , dlauny
!
! End of declarations rewritten by SPAG
!
!---------------------------------------------------------------
!     this generates a U-shaped domain with an elliptic inside.
!     then a Delauney triangulation is used to generate the mesh.
!     mesh needs to be post-processed -- see inmesh --
!---------------------------------------------------------------
!--------------------------------------------------------------
!     coordinates of nodal points
!--------------------------------------------------------------
   DATA x1/0.0/ , y1/0.0/ , x2/11.0/ , y2/5.5/ , nodexc/3/
   xcntr = 1.50
   ycntr = y2/2.0
   rad = 0.6
   nsec = 30
   nr = 3
!
!-----------------------------------------------------------------------
!     enter the four corner points.
!-----------------------------------------------------------------------
   Nx = 1
   X(Nx) = x1
   Y(Nx) = y1
   Nodcode(Nx) = 1
   Nx = Nx + 1
   X(Nx) = x2
   Y(Nx) = y1
   Nodcode(Nx) = 1
   Nx = Nx + 1
   X(Nx) = x2
   Y(Nx) = y2
   Nodcode(Nx) = 1
   Nx = Nx + 1
   X(Nx) = x1
   Y(Nx) = y2
   Nodcode(Nx) = 1
!
!     define axes of ellipse
!
!-----------------------------------------------------------------------
   pi = 4.0*atan(1.0)
   delr = rad/real(nr)
   DO i = 1 , nsec
      theta = 2.0*real(i-1)*pi/real(nsec)
      xnew = xcntr + rad*cos(theta)
      ynew = ycntr + rad*sin(theta)
      IF ( (xnew<x2) .AND. (xnew>x1) .AND. (ynew<y2) .AND. (ynew>y1) ) THEN
         Nx = Nx + 1
         X(Nx) = xnew
         Y(Nx) = ynew
         Nodcode(Nx) = nodexc
         arx = delr*cos(theta)
         ary = delr*sin(theta)
!
!     exclusion zone near the boundary
!
!         excl = 0.1*delr
         excl = 0.15*delr
         SPAG_Loop_2_1: DO
!
!     while inside domain do:
!
            xnew = X(Nx) + arx
            ynew = Y(Nx) + ary
            IF ( xnew>=x2 ) THEN
               X(Nx) = x2
               Nodcode(Nx) = 1
            ELSEIF ( xnew<=x1 ) THEN
               X(Nx) = x1
               Nodcode(Nx) = 1
            ELSEIF ( ynew>=y2 ) THEN
               Y(Nx) = y2
               Nodcode(Nx) = 1
            ELSEIF ( ynew<=y1 ) THEN
               Y(Nx) = y1
               Nodcode(Nx) = 1
            ELSE
               Nx = Nx + 1
               X(Nx) = xnew
               Y(Nx) = ynew
               Nodcode(Nx) = 0
               CALL clos2bdr(Nx,xnew,ynew,X,Y,x1,x2,y1,y2,excl,Nodcode)
            ENDIF
            arx = arx*1.1
            ary = ary*1.1
            excl = excl*1.1
            IF ( Nodcode(Nx)>0 ) EXIT SPAG_Loop_2_1
         ENDDO SPAG_Loop_2_1
      ENDIF
   ENDDO
!
   nemax = 1500
   CALL dlauny(X,Y,Nx,ijktr,nemax,Nelx)
!
!     print *, ' delauny -- nx, nelx ', nx, nelx
   DO j = 1 , Nx
      nel(j) = 0
   ENDDO
!
!     transpose ijktr into ijk and count the number of
!     elemnts to which each node belongs
!
   DO j = 1 , Nelx
      DO k = 1 , Node
         i = ijktr(j,k)
         Ijk(k,j) = i
         nel(i) = nel(i) + 1
      ENDDO
   ENDDO
!
!     this mesh needs cleaning up --
!
   CALL cleanel(Nelx,Ijk,Node,Nodcode,nodexc)
   CALL cleannods(Nx,X,Y,Nelx,Ijk,Node,Nodcode,Iperm)
!
!     take care of ordering within each element
!
   CALL chkelmt(X,Y,Nelx,Ijk,Node)
END SUBROUTINE fmesh9
!*==clos2bdr.f90 processed by SPAG 8.04RA 15:47 21 Feb 2024
!!SPAG Open source Personal, Educational or Academic User Clemson University  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
!-----------------------------------------------------------------------
SUBROUTINE clos2bdr(Nx,Xnew,Ynew,X,Y,X1,X2,Y1,Y2,Excl,Nodcode)
   USE ISO_FORTRAN_ENV                 
   IMPLICIT NONE
!
! Dummy argument declarations rewritten by SPAG
!
   INTEGER , INTENT(IN) :: Nx
   REAL(REAL64) , INTENT(IN) :: Xnew
   REAL(REAL64) , INTENT(IN) :: Ynew
   REAL(REAL64) , INTENT(INOUT) , DIMENSION(Nx) :: X
   REAL(REAL64) , INTENT(INOUT) , DIMENSION(Nx) :: Y
   REAL(REAL64) , INTENT(IN) :: X1
   REAL(REAL64) , INTENT(IN) :: X2
   REAL(REAL64) , INTENT(IN) :: Y1
   REAL(REAL64) , INTENT(IN) :: Y2
   REAL(REAL64) , INTENT(IN) :: Excl
   INTEGER , INTENT(OUT) , DIMENSION(Nx) :: Nodcode
!
! End of declarations rewritten by SPAG
!
!-----------------------------------------------------------------------
!     takes care of case where a point generated is too close to the
!     boundary  -- in this case projects the previous point to the
!     rectangle boundary == that makes some exclusion criterion
!     violated... does a simple job.
!-----------------------------------------------------------------------
   IF ( Xnew>=X2-Excl ) THEN
      X(Nx) = X2
      Y(Nx) = Y(Nx-1)
      Nodcode(Nx) = 1
   ENDIF
   IF ( Xnew<=X1+Excl ) THEN
      X(Nx) = X1
      Y(Nx) = Y(Nx-1)
      Nodcode(Nx) = 1
   ENDIF
   IF ( Ynew>=Y2-Excl ) THEN
      Y(Nx) = Y2
      X(Nx) = X(Nx-1)
      Nodcode(Nx) = 1
   ENDIF
   IF ( Ynew<=Y1+Excl ) THEN
      Y(Nx) = Y1
      X(Nx) = X(Nx-1)
      Nodcode(Nx) = 1
   ENDIF
!
END SUBROUTINE clos2bdr
