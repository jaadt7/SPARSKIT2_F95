!*==dblstr.f90 processed by SPAG 8.04RA 15:47 21 Feb 2024
!!SPAG Open source Personal, Educational or Academic User Clemson University  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
!----------------------------------------------------------------------c
!                          S P A R S K I T                             c
!----------------------------------------------------------------------c
!               REORDERING ROUTINES -- LEVEL SET BASED ROUTINES        c
!----------------------------------------------------------------------c
! dblstr   : doubled stripe partitioner
! rdis     : recursive dissection partitioner
! dse2way  : distributed site expansion usuing sites from dblstr
! dse      : distributed site expansion usuing sites from rdis
!------------- utility routines -----------------------------------------
! BFS      : Breadth-First search traversal algorithm
! add_lvst : routine to add a level -- used by BFS
! stripes  : finds the level set structure
! stripes0 : finds a trivial one-way partitioning from level-sets
! perphn   : finds a pseudo-peripheral node and performs a BFS from it.
! mapper4  : routine used by dse and dse2way to do center expansion
! get_domns: routine to find subdomaine from linked lists found by
!            mapper4.
! add_lk   : routine to add entry to linked list -- used by mapper4.
! find_ctr : routine to locate an approximate center of a subgraph.
! rversp   : routine to reverse a given permutation (e.g., for RCMK)
! maskdeg  : integer function to compute the `masked' of a node
!-----------------------------------------------------------------------
SUBROUTINE dblstr(N,Ja,Ia,Ip1,Ip2,Nfirst,Riord,Ndom,Map,Mapptr,Mask,Levels,Iwk)
   IMPLICIT NONE
!
! Dummy argument declarations rewritten by SPAG
!
   INTEGER :: N
   INTEGER , DIMENSION(*) :: Ja
   INTEGER , DIMENSION(*) :: Ia
   INTEGER :: Ip1
   INTEGER :: Ip2
   INTEGER :: Nfirst
   INTEGER , DIMENSION(*) :: Riord
   INTEGER , INTENT(INOUT) :: Ndom
   INTEGER , DIMENSION(*) :: Map
   INTEGER , INTENT(INOUT) , DIMENSION(*) :: Mapptr
   INTEGER , DIMENSION(*) :: Mask
   INTEGER , DIMENSION(*) :: Levels
   INTEGER , INTENT(INOUT) , DIMENSION(*) :: Iwk
!
! Local variable declarations rewritten by SPAG
!
   INTEGER :: idom , init , j , jdom , k , kdom , maskval , ndp1 , nextdom , nlev , numnod
   EXTERNAL bfs , perphn , stripes
!
! End of declarations rewritten by SPAG
!
!-----------------------------------------------------------------------
!     this routine does a two-way partitioning of a graph using
!     level sets recursively. First a coarse set is found by a
!     simple cuthill-mc Kee type algorithm. Them each of the large
!     domains is further partitioned into subsets using the same
!     technique. The ip1 and ip2 parameters indicate the desired number
!     number of partitions 'in each direction'. So the total number of
!     partitions on return ought to be equal (or close) to ip1*ip2
!----------------------parameters----------------------------------------
! on entry:
!---------
! n      = row dimension of matrix == number of vertices in graph
! ja, ia = pattern of matrix in CSR format (the ja,ia arrays of csr data
!          structure)
! ip1    = integer indicating the number of large partitions ('number of
!          paritions in first direction')
! ip2    = integer indicating the number of smaller partitions, per
!          large partition, ('number of partitions in second direction')
! nfirst = number of nodes in the first level that is input in riord
! riord  = (also an ouput argument). on entry riord contains the labels
!          of the nfirst nodes that constitute the first level.
! on return:
!-----------
! ndom   = total number of partitions found
! map    = list of nodes listed partition by partition from partition 1
!          to paritition ndom.
! mapptr = pointer array for map. All nodes from position
!          k1=mapptr(idom),to position k2=mapptr(idom+1)-1 in map belong
!          to partition idom.
! work arrays:
!-------------
! mask   = array of length n, used to hold the partition number of each
!          node for the first (large) partitioning.
!          mask is also used as a marker of  visited nodes.
! levels = integer array of length .le. n used to hold the pointer
!          arrays for the various level structures obtained from BFS.
!
!-----------------------------------------------------------------------
   maskval = 1
   DO j = 1 , N
      Mask(j) = maskval
   ENDDO
   Iwk(1) = 0
   CALL bfs(N,Ja,Ia,Nfirst,Iwk,Mask,maskval,Riord,Levels,nlev)
!
!     init = riord(1)
!     call perphn (ja,ia,mask,maskval,init,nlev,riord,levels)
   CALL stripes(nlev,Riord,Levels,Ip1,Map,Mapptr,Ndom)
!-----------------------------------------------------------------------
   IF ( Ip2==1 ) RETURN
   ndp1 = Ndom + 1
!
!     pack info into array iwk
!
   DO j = 1 , Ndom + 1
      Iwk(j) = ndp1 + Mapptr(j)
   ENDDO
   DO j = 1 , Mapptr(Ndom+1) - 1
      Iwk(ndp1+j) = Map(j)
   ENDDO
   DO idom = 1 , Ndom
      j = Iwk(idom)
      numnod = Iwk(idom+1) - Iwk(idom)
      init = Iwk(j)
      DO k = j , Iwk(idom+1) - 1
      ENDDO
   ENDDO
 
   DO idom = 1 , Ndom
      DO k = Mapptr(idom) , Mapptr(idom+1) - 1
         Mask(Map(k)) = idom
      ENDDO
   ENDDO
   nextdom = 1
!
!     jdom = counter for total number of (small) subdomains
!
   jdom = 1
   Mapptr(jdom) = 1
!-----------------------------------------------------------------------
   DO idom = 1 , Ndom
      maskval = idom
      Nfirst = 1
      numnod = Iwk(idom+1) - Iwk(idom)
      j = Iwk(idom)
      init = Iwk(j)
      nextdom = Mapptr(jdom)
!  note:    old version uses iperm array
      CALL perphn(numnod,Ja,Ia,init,Mask,maskval,nlev,Riord,Levels)
!
      CALL stripes(nlev,Riord,Levels,Ip2,Map(nextdom),Mapptr(jdom),kdom)
!
      Mapptr(jdom) = nextdom
      DO j = jdom , jdom + kdom - 1
         Mapptr(j+1) = nextdom + Mapptr(j+1) - 1
      ENDDO
      jdom = jdom + kdom
   ENDDO
!
   Ndom = jdom - 1
END SUBROUTINE dblstr
!*==rdis.f90 processed by SPAG 8.04RA 15:47 21 Feb 2024
!!SPAG Open source Personal, Educational or Academic User Clemson University  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
!-----------------------------------------------------------------------
SUBROUTINE rdis(N,Ja,Ia,Ndom,Map,Mapptr,Mask,Levels,Size,Iptr)
   IMPLICIT NONE
!
! Dummy argument declarations rewritten by SPAG
!
   INTEGER , INTENT(IN) :: Ndom
   INTEGER :: N
   INTEGER , DIMENSION(*) :: Ja
   INTEGER , DIMENSION(*) :: Ia
   INTEGER , INTENT(INOUT) , DIMENSION(*) :: Map
   INTEGER , INTENT(INOUT) , DIMENSION(*) :: Mapptr
   INTEGER , INTENT(INOUT) , DIMENSION(*) :: Mask
   INTEGER , DIMENSION(*) :: Levels
   INTEGER , INTENT(INOUT) , DIMENSION(Ndom) :: Size
   INTEGER , INTENT(INOUT) , DIMENSION(Ndom) :: Iptr
!
! Local variable declarations rewritten by SPAG
!
   INTEGER :: idom , init , j , k , ko , lev , maskval , maxsiz , nextdom , nextsiz , nlev , wantsiz
   EXTERNAL perphn
!
! End of declarations rewritten by SPAG
!
!-----------------------------------------------------------------------
!     recursive dissection algorithm for partitioning.
!     initial graph is cut in two - then each time, the largest set
!     is cut in two until we reach desired number of domains.
!-----------------------------------------------------------------------
!     input
!     n, ja, ia = graph
!     ndom      = desired number of subgraphs
!     output
!     ------
!     map, mapptr  = pointer array data structure for domains.
!             if k1 = mapptr(i), k2=mapptr(i+1)-1 then
!             map(k1:k2) = points in domain number i
!    work arrays:
!    -------------
!    mask(1:n)    integer
!    levels(1:n)  integer
!    size(1:ndom) integer
!    iptr(1:ndom) integer
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
   idom = 1
   nextdom = 0
!-----------------------------------------------------------------------
!     size(i) = size of domnain  i
!     iptr(i)  = index of first element of domain i
!-----------------------------------------------------------------------
   Size(idom) = N
   Iptr(idom) = 1
   DO j = 1 , N
      Mask(j) = 1
   ENDDO
   SPAG_Loop_1_1: DO
!
!     domain loop
!
!
!     select domain with largest size
!
      maxsiz = 0
      DO j = 1 , idom
         IF ( Size(j)>maxsiz ) THEN
            maxsiz = Size(j)
            nextdom = j
         ENDIF
      ENDDO
!
!     do a Prphn/ BFS on nextdom
!
      maskval = nextdom
      init = Iptr(nextdom)
      CALL perphn(N,Ja,Ia,init,Mask,maskval,nlev,Map,Levels)
!
!     determine next subdomain
!
      nextsiz = 0
      wantsiz = maxsiz/2
      idom = idom + 1
      lev = nlev
      DO WHILE ( nextsiz<wantsiz )
         DO k = Levels(lev) , Levels(lev+1) - 1
            Mask(Map(k)) = idom
         ENDDO
         nextsiz = nextsiz + Levels(lev+1) - Levels(lev)
         lev = lev - 1
      ENDDO
!
      Size(nextdom) = Size(nextdom) - nextsiz
      Size(idom) = nextsiz
!
!     new initial point = last point of previous domain
!
      Iptr(idom) = Map(Levels(nlev+1)-1)
!       iptr(idom) = map(levels(lev)+1)
!      iptr(idom) = 1
!
! alternative
!      lev = 1
!      do while (nextsiz .lt. wantsiz)
!         do k = levels(lev), levels(lev+1)-1
!            mask(map(k)) = idom
!         enddo
!         nextsiz = nextsiz + levels(lev+1) - levels(lev)
!         lev = lev+1
!      enddo
!
!     set size of new domain and adjust previous one
!
!      size(idom) = nextsiz
!      size(nextdom) = size(nextdom) - nextsiz
!      iptr(idom) = iptr(nextdom)
!      iptr(nextdom) = map(levels(lev))
 
      IF ( idom>=Ndom ) THEN
!
!     domains found -- build data structure
!
         Mapptr(1) = 1
         DO idom = 1 , Ndom
            Mapptr(idom+1) = Mapptr(idom) + Size(idom)
         ENDDO
         DO k = 1 , N
            idom = Mask(k)
            ko = Mapptr(idom)
            Map(ko) = k
            Mapptr(idom) = ko + 1
         ENDDO
!
!     reset pointers
!
         DO j = Ndom , 1 , -1
            Mapptr(j+1) = Mapptr(j)
         ENDDO
         Mapptr(1) = 1
         EXIT SPAG_Loop_1_1
      ENDIF
   ENDDO SPAG_Loop_1_1
!
END SUBROUTINE rdis
!*==dse2way.f90 processed by SPAG 8.04RA 15:47 21 Feb 2024
!!SPAG Open source Personal, Educational or Academic User Clemson University  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
!-----------------------------------------------------------------------
SUBROUTINE dse2way(N,Ja,Ia,Ip1,Ip2,Nfirst,Riord,Ndom,Dom,Idom,Mask,Jwk,Link)
!-----------------------------------------------------------------------
!     uses centers obtained from dblstr partition to get new partition
!-----------------------------------------------------------------------
!     input: n, ja, ia   = matrix
!     nfirst = number of first points
!     riord  = riord(1:nfirst) initial points
!     output
!     ndom   = number of domains
!     dom, idom = pointer array structure for domains.
!     mask , jwk, link = work arrays,
!-----------------------------------------------------------------------
   IMPLICIT NONE
!
! Dummy argument declarations rewritten by SPAG
!
   INTEGER :: N
   INTEGER , DIMENSION(*) :: Ja
   INTEGER , DIMENSION(*) :: Ia
   INTEGER :: Ip1
   INTEGER :: Ip2
   INTEGER :: Nfirst
   INTEGER , DIMENSION(*) :: Riord
   INTEGER :: Ndom
   INTEGER , DIMENSION(*) :: Dom
   INTEGER , DIMENSION(*) :: Idom
   INTEGER , DIMENSION(*) :: Mask
   INTEGER , DIMENSION(*) :: Jwk
   INTEGER , DIMENSION(*) :: Link
!
! Local variable declarations rewritten by SPAG
!
   INTEGER :: i , init , k , maskval , mid , nouter , nsiz , outer
   EXTERNAL dblstr , find_ctr , get_domns2 , mapper4
!
! End of declarations rewritten by SPAG
!
!
!-----------------------------------------------------------------------
!     local variables
   CALL dblstr(N,Ja,Ia,Ip1,Ip2,Nfirst,Riord,Ndom,Dom,Idom,Mask,Link,Jwk)
!
   nouter = 3
!-----------------------------------------------------------------------
 
   DO outer = 1 , nouter
!
!     set masks
!
      DO i = 1 , Ndom
         DO k = Idom(i) , Idom(i+1) - 1
            Mask(Dom(k)) = i
         ENDDO
      ENDDO
!
!     get centers
!
      DO i = 1 , Ndom
         nsiz = Idom(i+1) - Idom(i)
         init = Dom(Idom(i))
         maskval = i
!
!         use link for local riord -- jwk for other arrays --
!
         CALL find_ctr(N,nsiz,Ja,Ia,init,Mask,maskval,Link,Jwk,mid,Jwk(nsiz+1))
         Riord(i) = mid
      ENDDO
!
!     do level-set expansion from centers -- save previous diameter
!
      CALL mapper4(N,Ja,Ia,Ndom,Riord,Jwk,Mask,Link)
      CALL get_domns2(Ndom,Riord,Link,Jwk,Dom,Idom)
!-----------------------------------------------------------------------
   ENDDO
END SUBROUTINE dse2way
!*==dse.f90 processed by SPAG 8.04RA 15:47 21 Feb 2024
!!SPAG Open source Personal, Educational or Academic User Clemson University  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
!-----------------------------------------------------------------------
SUBROUTINE dse(N,Ja,Ia,Ndom,Riord,Dom,Idom,Mask,Jwk,Link)
   IMPLICIT NONE
!
! Dummy argument declarations rewritten by SPAG
!
   INTEGER :: N
   INTEGER , DIMENSION(*) :: Ja
   INTEGER , DIMENSION(*) :: Ia
   INTEGER :: Ndom
   INTEGER , DIMENSION(*) :: Riord
   INTEGER , DIMENSION(*) :: Dom
   INTEGER , DIMENSION(*) :: Idom
   INTEGER , DIMENSION(*) :: Mask
   INTEGER , DIMENSION(*) :: Jwk
   INTEGER , DIMENSION(*) :: Link
!
! Local variable declarations rewritten by SPAG
!
   INTEGER :: i , init , k , maskval , mid , nouter , nsiz , outer
   EXTERNAL find_ctr , get_domns2 , mapper4 , rdis
!
! End of declarations rewritten by SPAG
!
!-----------------------------------------------------------------------
!     uses centers produced from rdis to get a new partitioning --
!     see calling sequence in rdis..
!-----------------------------------------------------------------------
!     local variables
!-----------------------------------------------------------------------
   nouter = 3
!
   CALL rdis(N,Ja,Ia,Ndom,Dom,Idom,Mask,Link,Jwk,Jwk(Ndom+1))
!
!     initial points =
!
   DO outer = 1 , nouter
!
!     set masks
!
      DO i = 1 , Ndom
         DO k = Idom(i) , Idom(i+1) - 1
            Mask(Dom(k)) = i
         ENDDO
      ENDDO
!
!     get centers
!
      DO i = 1 , Ndom
         nsiz = Idom(i+1) - Idom(i)
         init = Dom(Idom(i))
         maskval = i
!
!         use link for local riord -- jwk for other arrays --
!
 
         CALL find_ctr(N,nsiz,Ja,Ia,init,Mask,maskval,Link,Jwk,mid,Jwk(nsiz+1))
         Riord(i) = mid
      ENDDO
!
!     do level-set expansion from centers -- save previous diameter
!
      CALL mapper4(N,Ja,Ia,Ndom,Riord,Jwk,Mask,Link)
      CALL get_domns2(Ndom,Riord,Link,Jwk,Dom,Idom)
!-----------------------------------------------------------------------
   ENDDO
END SUBROUTINE dse
!*==bfs.f90 processed by SPAG 8.04RA 15:47 21 Feb 2024
!!SPAG Open source Personal, Educational or Academic User Clemson University  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
!-----------------------------------------------------------------------
SUBROUTINE bfs(N,Ja,Ia,Nfirst,Iperm,Mask,Maskval,Riord,Levels,Nlev)
   IMPLICIT NONE
!
! Dummy argument declarations rewritten by SPAG
!
   INTEGER , INTENT(IN) :: N
   INTEGER , DIMENSION(*) :: Ja
   INTEGER , DIMENSION(*) :: Ia
   INTEGER , INTENT(IN) :: Nfirst
   INTEGER , INTENT(IN) , DIMENSION(N) :: Iperm
   INTEGER , INTENT(INOUT) , DIMENSION(N) :: Mask
   INTEGER :: Maskval
   INTEGER , INTENT(INOUT) , DIMENSION(*) :: Riord
   INTEGER , INTENT(OUT) , DIMENSION(*) :: Levels
   INTEGER , INTENT(INOUT) :: Nlev
!
! Local variable declarations rewritten by SPAG
!
   INTEGER :: iend , ii , istart , j , nod
   LOGICAL :: permut
   EXTERNAL add_lvst
!
! End of declarations rewritten by SPAG
!
!-----------------------------------------------------------------------
! finds the level-structure (breadth-first-search or CMK) ordering for a
! given sparse matrix. Uses add_lvst. Allows an set of nodes to be
! the initial level (instead of just one node).
!-------------------------parameters------------------------------------
! on entry:
!---------
!     n      = number of nodes in the graph
!     ja, ia = pattern of matrix in CSR format (the ja,ia arrays of csr data
!     structure)
!     nfirst = number of nodes in the first level that is input in riord
!     iperm  = integer array indicating in which order to  traverse the graph
!     in order to generate all connected components.
!     if iperm(1) .eq. 0 on entry then BFS will traverse the nodes
!     in the  order 1,2,...,n.
!
!     riord  = (also an ouput argument). On entry riord contains the labels
!     of the nfirst nodes that constitute the first level.
!
!     mask   = array used to indicate whether or not a node should be
!     condidered in the graph. see maskval.
!     mask is also used as a marker of  visited nodes.
!
!     maskval= consider node i only when:  mask(i) .eq. maskval
!     maskval must be .gt. 0.
!     thus, to consider all nodes, take mask(1:n) = 1.
!     maskval=1 (for example)
!
!     on return
!     ---------
!     mask   = on return mask is restored to its initial state.
!     riord  = `reverse permutation array'. Contains the labels of the nodes
!     constituting all the levels found, from the first level to
!     the last.
!     levels = pointer array for the level structure. If lev is a level
!     number, and k1=levels(lev),k2=levels(lev+1)-1, then
!     all the nodes of level number lev are:
!     riord(k1),riord(k1+1),...,riord(k2)
!     nlev   = number of levels found
!-----------------------------------------------------------------------
!
   permut = (Iperm(1)/=0)
!
!     start pointer structure to levels
!
   Nlev = 0
!
!     previous end
!
   istart = 0
   ii = 0
!
!     current end
!
   iend = Nfirst
!
!     intialize masks to zero -- except nodes of first level --
!
   DO j = 1 , Nfirst
      Mask(Riord(j)) = 0
   ENDDO
   SPAG_Loop_1_2: DO
!-----------------------------------------------------------------------
!
      Nlev = Nlev + 1
      Levels(Nlev) = istart + 1
      CALL add_lvst(istart,iend,Riord,Ja,Ia,Mask,Maskval)
      IF ( istart>=iend ) THEN
         SPAG_Loop_2_1: DO
            ii = ii + 1
            IF ( ii<=N ) THEN
               nod = ii
               IF ( permut ) nod = Iperm(nod)
               IF ( Mask(nod)/=Maskval ) CYCLE
!
!     start a new level
!
               istart = iend
               iend = iend + 1
               Riord(iend) = nod
               Mask(nod) = 0
               EXIT SPAG_Loop_2_1
            ENDIF
!-----------------------------------------------------------------------
            Levels(Nlev+1) = iend + 1
            DO j = 1 , iend
               Mask(Riord(j)) = Maskval
            ENDDO
            EXIT SPAG_Loop_1_2
         ENDDO SPAG_Loop_2_1
      ENDIF
   ENDDO SPAG_Loop_1_2
!-----------------------------------------------------------------------
END SUBROUTINE bfs
!*==add_lvst.f90 processed by SPAG 8.04RA 15:47 21 Feb 2024
!!SPAG Open source Personal, Educational or Academic User Clemson University  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
!-----------------------------------------------------------------------
SUBROUTINE add_lvst(Istart,Iend,Riord,Ja,Ia,Mask,Maskval)
   IMPLICIT NONE
!
! Dummy argument declarations rewritten by SPAG
!
   INTEGER , INTENT(INOUT) :: Istart
   INTEGER , INTENT(INOUT) :: Iend
   INTEGER , INTENT(INOUT) , DIMENSION(*) :: Riord
   INTEGER , INTENT(IN) , DIMENSION(*) :: Ja
   INTEGER , INTENT(IN) , DIMENSION(*) :: Ia
   INTEGER , INTENT(INOUT) , DIMENSION(*) :: Mask
   INTEGER , INTENT(IN) :: Maskval
!
! Local variable declarations rewritten by SPAG
!
   INTEGER :: i , ir , j , k , nod
!
! End of declarations rewritten by SPAG
!
!-------------------------------------------------------------
!     adds one level set to the previous sets..
!     span all nodes of previous mask
!-------------------------------------------------------------
   nod = Iend
   DO ir = Istart + 1 , Iend
      i = Riord(ir)
      DO k = Ia(i) , Ia(i+1) - 1
         j = Ja(k)
         IF ( Mask(j)==Maskval ) THEN
            nod = nod + 1
            Mask(j) = 0
            Riord(nod) = j
         ENDIF
      ENDDO
   ENDDO
   Istart = Iend
   Iend = nod
END SUBROUTINE add_lvst
!*==stripes.f90 processed by SPAG 8.04RA 15:47 21 Feb 2024
!!SPAG Open source Personal, Educational or Academic User Clemson University  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
!-----------------------------------------------------------------------
SUBROUTINE stripes(Nlev,Riord,Levels,Ip,Map,Mapptr,Ndom)
   IMPLICIT NONE
!
! Dummy argument declarations rewritten by SPAG
!
   INTEGER , INTENT(IN) :: Nlev
   INTEGER , INTENT(IN) , DIMENSION(*) :: Riord
   INTEGER , INTENT(IN) , DIMENSION(Nlev+1) :: Levels
   INTEGER , INTENT(IN) :: Ip
   INTEGER , INTENT(OUT) , DIMENSION(*) :: Map
   INTEGER , INTENT(OUT) , DIMENSION(*) :: Mapptr
   INTEGER , INTENT(INOUT) :: Ndom
!
! Local variable declarations rewritten by SPAG
!
   INTEGER :: ib , ilev , k , ktr , nsiz , psiz
!
! End of declarations rewritten by SPAG
!
!-----------------------------------------------------------------------
!    this is a post processor to BFS. stripes uses the output of BFS to
!    find a decomposition of the adjacency graph by stripes. It fills
!    the stripes level by level until a number of nodes .gt. ip is
!    is reached.
!---------------------------parameters-----------------------------------
! on entry:
! --------
! nlev   = number of levels as found by BFS
! riord  = reverse permutation array produced by BFS --
! levels = pointer array for the level structure as computed by BFS. If
!          lev is a level number, and k1=levels(lev),k2=levels(lev+1)-1,
!          then all the nodes of level number lev are:
!                      riord(k1),riord(k1+1),...,riord(k2)
!  ip    = number of desired partitions (subdomains) of about equal size.
!
! on return
! ---------
! ndom     = number of subgraphs (subdomains) found
! map      = node per processor list. The nodes are listed contiguously
!            from proc 1 to nproc = mpx*mpy.
! mapptr   = pointer array for array map. list for proc. i starts at
!            mapptr(i) and ends at mapptr(i+1)-1 in array map.
!-----------------------------------------------------------------------
! local variables.
!
   Ndom = 1
   ib = 1
! to add: if (ip .le. 1) then ...
   nsiz = Levels(Nlev+1) - Levels(1)
   psiz = (nsiz-ib)/max(1,(Ip-Ndom+1)) + 1
   Mapptr(Ndom) = ib
   ktr = 0
   DO ilev = 1 , Nlev
!
!     add all nodes of this level to domain
!
      DO k = Levels(ilev) , Levels(ilev+1) - 1
         Map(ib) = Riord(k)
         ib = ib + 1
         ktr = ktr + 1
         IF ( ktr>=psiz .OR. k>=nsiz ) THEN
            Ndom = Ndom + 1
            Mapptr(Ndom) = ib
            psiz = (nsiz-ib)/max(1,(Ip-Ndom+1)) + 1
            ktr = 0
         ENDIF
!
      ENDDO
   ENDDO
   Ndom = Ndom - 1
END SUBROUTINE stripes
!*==stripes0.f90 processed by SPAG 8.04RA 15:47 21 Feb 2024
!!SPAG Open source Personal, Educational or Academic User Clemson University  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
!-----------------------------------------------------------------------
SUBROUTINE stripes0(Ip,Nlev,Il,Ndom,Iptr)
   IMPLICIT NONE
!
! Dummy argument declarations rewritten by SPAG
!
   INTEGER , INTENT(IN) :: Ip
   INTEGER , INTENT(IN) :: Nlev
   INTEGER , INTENT(IN) , DIMENSION(*) :: Il
   INTEGER , INTENT(OUT) :: Ndom
   INTEGER , INTENT(OUT) , DIMENSION(*) :: Iptr
!
! Local variable declarations rewritten by SPAG
!
   INTEGER :: iband , ilev , ktr
!
! End of declarations rewritten by SPAG
!
!-----------------------------------------------------------------------
!     This routine is a simple level-set partitioner. It scans
!     the level-sets as produced by BFS from one to nlev.
!     each time the number of nodes in the accumulated set of
!     levels traversed exceeds the parameter ip, this set defines
!     a new subgraph.
!-------------------------parameter-list---------------------------------
! on entry:
! --------
! ip     = desired number of nodes per subgraph.
! nlev   = number of levels found  as output by BFS
! il     = integer array containing the pointer array for
!          the level data structure as output by BFS.
!          thus il(lev+1) - il(lev) = the number of
!          nodes that constitute the level numbe lev.
! on return
! ---------
! ndom   = number of sungraphs found
! iptr   = pointer array for the sugraph data structure.
!          thus, iptr(idom) points to the first level that
!          consistutes the subgraph number idom, in the
!          level data structure.
!-----------------------------------------------------------------------
   ktr = 0
   iband = 1
   Iptr(iband) = 1
!-----------------------------------------------------------------------
 
   DO ilev = 1 , Nlev
      ktr = ktr + Il(ilev+1) - Il(ilev)
      IF ( ktr>Ip ) THEN
         iband = iband + 1
         Iptr(iband) = ilev + 1
         ktr = 0
      ENDIF
!
   ENDDO
!-----------returning --------------------
   Iptr(iband) = Nlev + 1
   Ndom = iband - 1
!-----------------------------------------------------------------------
!-----end-of-stripes0---------------------------------------------------
END SUBROUTINE stripes0
!*==maskdeg.f90 processed by SPAG 8.04RA 15:47 21 Feb 2024
!!SPAG Open source Personal, Educational or Academic User Clemson University  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
!-----------------------------------------------------------------------
FUNCTION maskdeg(Ja,Ia,Nod,Mask,Maskval)
   IMPLICIT NONE
!
! Function and Dummy argument declarations rewritten by SPAG
!
   INTEGER :: maskdeg
   INTEGER , INTENT(IN) , DIMENSION(*) :: Ja
   INTEGER , INTENT(IN) , DIMENSION(*) :: Ia
   INTEGER , INTENT(IN) :: Nod
   INTEGER , INTENT(IN) , DIMENSION(*) :: Mask
   INTEGER , INTENT(IN) :: Maskval
!
! Local variable declarations rewritten by SPAG
!
   INTEGER :: deg , k
!
! End of declarations rewritten by SPAG
!
!-----------------------------------------------------------------------
   deg = 0
   DO k = Ia(Nod) , Ia(Nod+1) - 1
      IF ( Mask(Ja(k))==Maskval ) deg = deg + 1
   ENDDO
   maskdeg = deg
END FUNCTION maskdeg
!*==perphn.f90 processed by SPAG 8.04RA 15:47 21 Feb 2024
!!SPAG Open source Personal, Educational or Academic User Clemson University  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
!-----------------------------------------------------------------------
SUBROUTINE perphn(N,Ja,Ia,Init,Mask,Maskval,Nlev,Riord,Levels)
   IMPLICIT NONE
!
! Dummy argument declarations rewritten by SPAG
!
   INTEGER :: N
   INTEGER , DIMENSION(*) :: Ja
   INTEGER , DIMENSION(*) :: Ia
   INTEGER , INTENT(INOUT) :: Init
   INTEGER , DIMENSION(*) :: Mask
   INTEGER :: Maskval
   INTEGER :: Nlev
   INTEGER , INTENT(INOUT) , DIMENSION(*) :: Riord
   INTEGER , DIMENSION(*) :: Levels
!
! Local variable declarations rewritten by SPAG
!
   INTEGER :: deg , j , mindeg , nfirst , nlevp , nod
   INTEGER , DIMENSION(1) :: iperm
   INTEGER , EXTERNAL :: maskdeg
   EXTERNAL bfs
!
! End of declarations rewritten by SPAG
!
!-----------------------------------------------------------------------
!     finds a peripheral node and does a BFS search from it.
!-----------------------------------------------------------------------
!     see routine  dblstr for description of parameters
! input:
!-------
! ja, ia  = list pointer array for the adjacency graph
! mask    = array used for masking nodes -- see maskval
! maskval = value to be checked against for determing whether or
!           not a node is masked. If mask(k) .ne. maskval then
!           node k is not considered.
! init    = init node in the pseudo-peripheral node algorithm.
!
! output:
!-------
! init    = actual pseudo-peripherial node found.
! nlev    = number of levels in the final BFS traversal.
! riord   =
! levels  =
!-----------------------------------------------------------------------
   nlevp = 0
   SPAG_Loop_1_1: DO
      Riord(1) = Init
      nfirst = 1
      iperm(1) = 0
!
      CALL bfs(N,Ja,Ia,nfirst,iperm,Mask,Maskval,Riord,Levels,Nlev)
      IF ( Nlev>nlevp ) THEN
         mindeg = N + 1
         DO j = Levels(Nlev) , Levels(Nlev+1) - 1
            nod = Riord(j)
            deg = maskdeg(Ja,Ia,nod,Mask,Maskval)
            IF ( deg<mindeg ) THEN
               Init = nod
               mindeg = deg
            ENDIF
         ENDDO
         nlevp = Nlev
         CYCLE
      ENDIF
      EXIT SPAG_Loop_1_1
   ENDDO SPAG_Loop_1_1
END SUBROUTINE perphn
!*==mapper4.f90 processed by SPAG 8.04RA 15:47 21 Feb 2024
!!SPAG Open source Personal, Educational or Academic User Clemson University  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
!-----------------------------------------------------------------------
SUBROUTINE mapper4(N,Ja,Ia,Ndom,Nodes,Levst,Marker,Link)
   IMPLICIT NONE
!
! Dummy argument declarations rewritten by SPAG
!
   INTEGER :: N
   INTEGER :: Ndom
   INTEGER , INTENT(IN) , DIMENSION(*) :: Ja
   INTEGER , INTENT(IN) , DIMENSION(*) :: Ia
   INTEGER , DIMENSION(*) :: Nodes
   INTEGER , INTENT(INOUT) , DIMENSION(2*Ndom) :: Levst
   INTEGER , INTENT(INOUT) , DIMENSION(N) :: Marker
   INTEGER , INTENT(INOUT) , DIMENSION(*) :: Link
!
! Local variable declarations rewritten by SPAG
!
   INTEGER :: i , idom , ii , ilast , isiz , j , kk , lkend , next , nod , nodprev , nsize , nstuck
   INTEGER , EXTERNAL :: mindom
   EXTERNAL add_lk
!
! End of declarations rewritten by SPAG
!
!-----------------------------------------------------------------------
!     finds domains given ndom centers -- by doing a level set expansion
!-----------------------------------------------------------------------
!     on entry:
!     ---------
!     n      = dimension of matrix
!     ja, ia = adajacency list of matrix (CSR format without values) --
!     ndom   = number of subdomains (nr output by coarsen)
!     nodes  = array of size at least n. On input the first ndom entries
!              of nodes should contain the labels of the centers of the
!              ndom domains from which to do the expansions.
!
!     on return
!     ---------
!     link  = linked list array for the ndom domains.
!     nodes = contains the list of nodes of the domain corresponding to
!             link. (nodes(i) and link(i) are related to the same node).
!
!     levst = levst(j) points to beginning of subdomain j in link.
!
!     work arrays:
!     -----------
!     levst : work array of length 2*ndom -- contains the beginning and
!     end of  current level in link.
!     beginning of last level in link for each processor.
!     also ends in levst(ndom+i)
!     marker : work array of length n.
!
!     Notes on implementation:
!     -----------------------
!     for j .le. ndom link(j) is <0  and indicates the end of the
!     linked list. The most recent element added to the linked
!     list is added at the end of the list (traversal=backward)
!     For  j .le. ndom, the value of -link(j) is the size of
!     subdomain j.
!
!-----------------------------------------------------------------------
!     local variables
   INTEGER :: spag_nextblock_1
   ii = 0
!
   lkend = Ndom
   nod = Ndom
   nstuck = 0
   spag_nextblock_1 = 1
   SPAG_DispatchLoop_1: DO
      SELECT CASE (spag_nextblock_1)
      CASE (1)
!
!     initilaize nodes and link arrays
!
         DO j = 1 , N
            Marker(j) = 0
         ENDDO
!
         DO j = 1 , Ndom
            Link(j) = -1
            Marker(Nodes(j)) = j
            Levst(j) = j
            Levst(Ndom+j) = j
         ENDDO
!
!     ii = next untouched node for restarting new connected component.
!
         
      CASE (2)
!-----------------------------------------------------------------------
         idom = mindom(N,Ndom,Link)
         spag_nextblock_1 = 3
      CASE (3)
!-----------------------------------------------------------------------
!     begin level-set loop
!-----------------------------------------------------------------------
         nodprev = nod
         ilast = Levst(Ndom+idom)
         Levst(Ndom+idom) = lkend
         next = Levst(idom)
!
!     linked list traversal loop
!
         isiz = 0
         nsize = Link(idom)
         DO
            i = Nodes(next)
            isiz = isiz + 1
!
!     adjacency list traversal loop
!
            DO kk = Ia(i) , Ia(i+1) - 1
               j = Ja(kk)
               IF ( Marker(j)==0 ) CALL add_lk(j,nod,idom,lkend,Levst,Link,Nodes,Marker)
            ENDDO
!
!     if last element of the previous level not reached continue
!
            IF ( next>ilast ) THEN
               next = Link(next)
               IF ( next>0 ) CYCLE
            ENDIF
            spag_nextblock_1 = 4
            CYCLE SPAG_DispatchLoop_1
         ENDDO
         spag_nextblock_1 = 4
      CASE (4)
!-----------------------------------------------------------------------
!     end level-set traversal --
!-----------------------------------------------------------------------
         IF ( nodprev==nod ) THEN
!
!     link(idom) >0 indicates that set is stuck  --
!
            Link(idom) = -Link(idom)
            nstuck = nstuck + 1
         ENDIF
!
         IF ( nstuck<Ndom ) THEN
            spag_nextblock_1 = 2
            CYCLE SPAG_DispatchLoop_1
         ENDIF
!
!     reset sizes --
!
         DO j = 1 , Ndom
            IF ( Link(j)>0 ) Link(j) = -Link(j)
         ENDDO
!
         IF ( nod==N ) RETURN
         DO
!
!     stuck. add first non assigned point to smallest domain
!
            ii = ii + 1
            IF ( ii<=N ) THEN
               IF ( Marker(ii)/=0 ) CYCLE
               idom = 0
               isiz = N + 1
               DO kk = Ia(ii) , Ia(ii+1) - 1
                  i = Marker(Ja(kk))
                  IF ( i/=0 ) THEN
                     nsize = abs(Link(i))
                     IF ( nsize<isiz ) THEN
                        isiz = nsize
                        idom = i
                     ENDIF
                  ENDIF
               ENDDO
!
!     if no neighboring domain select smallest one
!
               IF ( idom==0 ) idom = mindom(N,Ndom,Link)
!
!     add ii to sudomain idom at end of linked list
!
               CALL add_lk(ii,nod,idom,lkend,Levst,Link,Nodes,Marker)
               spag_nextblock_1 = 3
               CYCLE SPAG_DispatchLoop_1
            ENDIF
            EXIT SPAG_DispatchLoop_1
         ENDDO
         EXIT SPAG_DispatchLoop_1
      END SELECT
   ENDDO SPAG_DispatchLoop_1
END SUBROUTINE mapper4
!*==get_domns2.f90 processed by SPAG 8.04RA 15:47 21 Feb 2024
!!SPAG Open source Personal, Educational or Academic User Clemson University  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
!-----------------------------------------------------------------------
SUBROUTINE get_domns2(Ndom,Nodes,Link,Levst,Riord,Iptr)
   IMPLICIT NONE
!
! Dummy argument declarations rewritten by SPAG
!
   INTEGER , INTENT(IN) :: Ndom
   INTEGER , INTENT(IN) , DIMENSION(*) :: Nodes
   INTEGER , INTENT(IN) , DIMENSION(*) :: Link
   INTEGER , INTENT(IN) , DIMENSION(*) :: Levst
   INTEGER , INTENT(OUT) , DIMENSION(*) :: Riord
   INTEGER , INTENT(OUT) , DIMENSION(*) :: Iptr
!
! Local variable declarations rewritten by SPAG
!
   INTEGER :: ii , j , next , nod
!
! End of declarations rewritten by SPAG
!
!-----------------------------------------------------------------------
!     constructs the subdomains from its linked list data structure
!-----------------------------------------------------------------------
!     input:
!     ndom  = number of subdomains
!     nodes = sequence of nodes are produced by mapper4.
!     link  = link list array as produced by mapper4.
!     on return:
!----------
!     riord = contains the nodes in each subdomain in succession.
!     iptr  = pointer in riord for beginnning of each subdomain.
!     Thus subdomain number i consists of nodes
!     riord(k1),riord(k1)+1,...,riord(k2)
!     where k1 = iptr(i), k2= iptr(i+1)-1
!
!-----------------------------------------------------------------------
!     local variables
   nod = 1
   Iptr(1) = nod
   DO j = 1 , Ndom
      next = Levst(j)
      SPAG_Loop_2_1: DO
         ii = Nodes(next)
         Riord(nod) = ii
         nod = nod + 1
         next = Link(next)
         IF ( next<=0 ) THEN
            Iptr(j+1) = nod
            EXIT SPAG_Loop_2_1
         ENDIF
      ENDDO SPAG_Loop_2_1
   ENDDO
!
!-----------------------------------------------------------------------
END SUBROUTINE get_domns2
!*==mindom.f90 processed by SPAG 8.04RA 15:47 21 Feb 2024
!!SPAG Open source Personal, Educational or Academic User Clemson University  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
!-----------------------------------------------------------------------
FUNCTION mindom(N,Ndom,Link)
   IMPLICIT NONE
!
! Function and Dummy argument declarations rewritten by SPAG
!
   INTEGER , INTENT(IN) :: N
   INTEGER :: mindom
   INTEGER , INTENT(IN) :: Ndom
   INTEGER , INTENT(IN) , DIMENSION(N) :: Link
!
! Local variable declarations rewritten by SPAG
!
   INTEGER :: i , isiz , nsize
!
! End of declarations rewritten by SPAG
!
!-----------------------------------------------------------------------
!     returns  the domain with smallest size
!-----------------------------------------------------------------------
!      locals
!
!
   mindom = 0
   isiz = N + 1
   DO i = 1 , Ndom
      nsize = -Link(i)
      IF ( nsize>=0 ) THEN
         IF ( nsize<isiz ) THEN
            isiz = nsize
            mindom = i
         ENDIF
      ENDIF
   ENDDO
END FUNCTION mindom
!*==add_lk.f90 processed by SPAG 8.04RA 15:47 21 Feb 2024
!!SPAG Open source Personal, Educational or Academic User Clemson University  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
!-----------------------------------------------------------------------
SUBROUTINE add_lk(New,Nod,Idom,Lkend,Levst,Link,Nodes,Marker)
   IMPLICIT NONE
!
! Dummy argument declarations rewritten by SPAG
!
   INTEGER , INTENT(IN) :: New
   INTEGER , INTENT(INOUT) :: Nod
   INTEGER , INTENT(IN) :: Idom
   INTEGER , INTENT(INOUT) :: Lkend
   INTEGER , INTENT(INOUT) , DIMENSION(*) :: Levst
   INTEGER , INTENT(INOUT) , DIMENSION(*) :: Link
   INTEGER , INTENT(OUT) , DIMENSION(*) :: Nodes
   INTEGER , INTENT(OUT) , DIMENSION(*) :: Marker
!
! Local variable declarations rewritten by SPAG
!
   INTEGER :: ktop
!
! End of declarations rewritten by SPAG
!
!-----------------------------------------------------------------------
!     inserts new element to linked list from the tail.
!-----------------------------------------------------------------------
!     adds one entry (new) to linked list and ipdates everything.
!     new  = node to be added
!     nod  = current number of marked nodes
!     idom = domain to which new is to be added
!     ndom = total number of domains
!     lkend= location of end of structure (link and nodes)
!     levst= pointer array for link, nodes
!     link = link array
!     nodes= nodes array --
!     marker = marker array == if marker(k) =0 then node k is not
!              assigned yet.
!-----------------------------------------------------------------------
!      locals
!
   Lkend = Lkend + 1
   Nodes(Lkend) = New
   Nod = Nod + 1
   Marker(New) = Idom
   ktop = Levst(Idom)
   Link(Lkend) = ktop
   Link(Idom) = Link(Idom) - 1
   Levst(Idom) = Lkend
!-----------------------------------------------------------------------
!-------end-of-add_lk---------------------------------------------------
END SUBROUTINE add_lk
!*==find_ctr.f90 processed by SPAG 8.04RA 15:47 21 Feb 2024
!!SPAG Open source Personal, Educational or Academic User Clemson University  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
!-----------------------------------------------------------------------
SUBROUTINE find_ctr(N,Nsiz,Ja,Ia,Init,Mask,Maskval,Riord,Levels,Center,Iwk)
   IMPLICIT NONE
!
! Dummy argument declarations rewritten by SPAG
!
   INTEGER :: N
   INTEGER , INTENT(IN) :: Nsiz
   INTEGER , DIMENSION(*) :: Ja
   INTEGER , DIMENSION(*) :: Ia
   INTEGER :: Init
   INTEGER , DIMENSION(*) :: Mask
   INTEGER :: Maskval
   INTEGER , DIMENSION(*) :: Riord
   INTEGER , DIMENSION(*) :: Levels
   INTEGER , INTENT(OUT) :: Center
   INTEGER , DIMENSION(*) :: Iwk
!
! Local variable declarations rewritten by SPAG
!
   INTEGER :: init0 , k , kl , kr , midlev , newmask , nlev , nlev0
   EXTERNAL perphn
!
! End of declarations rewritten by SPAG
!
!-----------------------------------------------------------------------
!     finds a center point of a subgraph --
!-----------------------------------------------------------------------
!     n, ja, ia = graph
!     nsiz = size of current domain.
!     init = initial node in search
!     mask
!     maskval
!-----------------------------------------------------------------------
!     local variables
   CALL perphn(N,Ja,Ia,Init,Mask,Maskval,nlev,Riord,Levels)
!-----------------------------------------------------------------------
!     midlevel = level which cuts domain into 2 roughly equal-size
!     regions
!
   midlev = 1
   k = 0
   SPAG_Loop_1_1: DO
      k = k + Levels(midlev+1) - Levels(midlev)
      IF ( k*2<Nsiz ) THEN
         midlev = midlev + 1
         CYCLE
      ENDIF
!-----------------------------------------------------------------------
      newmask = N + Maskval
!
!     assign temporary masks to mid-level elements
!
      DO k = Levels(midlev) , Levels(midlev+1) - 1
         Mask(Riord(k)) = newmask
      ENDDO
!
!     find pseudo-periph node for mid-level `line'
!
      kr = 1
      kl = kr + Nsiz
      init0 = Riord(Levels(midlev))
      CALL perphn(N,Ja,Ia,init0,Mask,newmask,nlev0,Iwk(kr),Iwk(kl))
!-----------------------------------------------------------------------
!     restore  mask to initial state
!-----------------------------------------------------------------------
      DO k = Levels(midlev) , Levels(midlev+1) - 1
         Mask(Riord(k)) = Maskval
      ENDDO
!-----------------------------------------------------------------------
!     define center
!-----------------------------------------------------------------------
      midlev = 1 + (nlev0-1)/2
      k = Iwk(kl+midlev-1)
      Center = Iwk(k)
      EXIT SPAG_Loop_1_1
   ENDDO SPAG_Loop_1_1
!-----------------------------------------------------------------------
END SUBROUTINE find_ctr
!*==rversp.f90 processed by SPAG 8.04RA 15:47 21 Feb 2024
!!SPAG Open source Personal, Educational or Academic User Clemson University  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
!-----------------------------------------------------------------------
SUBROUTINE rversp(N,Riord)
   IMPLICIT NONE
!
! Dummy argument declarations rewritten by SPAG
!
   INTEGER , INTENT(IN) :: N
   INTEGER , INTENT(INOUT) , DIMENSION(N) :: Riord
!
! Local variable declarations rewritten by SPAG
!
   INTEGER :: j , k
!
! End of declarations rewritten by SPAG
!
!-----------------------------------------------------------------------
!     this routine does an in-place reversing of the permutation array
!     riord --
!-----------------------------------------------------------------------
   DO j = 1 , N/2
      k = Riord(j)
      Riord(j) = Riord(N-j+1)
      Riord(N-j+1) = k
   ENDDO
END SUBROUTINE rversp
