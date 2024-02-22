!*==blccnx.f90 processed by SPAG 8.04RA 15:47 21 Feb 2024
!!SPAG Open source Personal, Educational or Academic User Clemson University  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
!----------------------------------------------------------------------c
!                          S P A R S K I T                             c
!----------------------------------------------------------------------c
!          REORDERING ROUTINES -- STRONGLY CONNECTED COMPONENTS        c
!----------------------------------------------------------------------c
!     Contributed by:
!     Laura C. Dutto - email: dutto@cerca.umontreal.ca
!                      July 1992 - Update: March 1994
!-----------------------------------------------------------------------
!     CONTENTS:
!     --------
!     blccnx : Driver routine to reduce the structure of a  matrix
!              to its strongly connected components.
!     cconex : Main routine to compute the strongly connected components
!              of a (block diagonal) matrix.
!     anccnx : We put in ICCNEX the vertices marked in the component MCCNEX.
!     newcnx : We put in ICCNEX the vertices marked in the component
!              MCCNEX. We modify also the vector KPW.
!     blccn1 : Parallel computation of the connected components of a
!              matrix. The parallel loop is performed only if the matrix
!              has a block diagonal structure.
!     ccnicopy:We copy an integer vector into anothoer.
!     compos : We calculate the composition between two permutation
!              vectors.
!     invlpw : We calculate the inverse of a permutation vector.
!     numini : We initialize a vector to the identity.
!     tbzero : We initialize to ZERO an integer vector.
!     iplusa : Given two integers IALPHA and IBETA, for an integer vector
!              IA we calculate IA(i) = ialpha + ibeta * ia(i)
!
!----------------------------------------------------------------------c
SUBROUTINE blccnx(N,Nbloc,Nblcmx,Nsbloc,Job,Lpw,Amat,Ja,Ia,Iout,Ier,Izs,Nw)
   USE ISO_FORTRAN_ENV                 
   IMPLICIT NONE
!
! Dummy argument declarations rewritten by SPAG
!
   INTEGER :: N
   INTEGER :: Nblcmx
   INTEGER , INTENT(IN) :: Nw
   INTEGER :: Nbloc
   INTEGER , DIMENSION(0:Nblcmx) :: Nsbloc
   INTEGER :: Job
   INTEGER , DIMENSION(N) :: Lpw
   REAL(REAL64) , DIMENSION(*) :: Amat
   INTEGER , DIMENSION(*) :: Ja
   INTEGER , DIMENSION(N+1) :: Ia
   INTEGER :: Iout
   INTEGER , INTENT(INOUT) :: Ier
   INTEGER , DIMENSION(Nw) :: Izs
!
! Local variable declarations rewritten by SPAG
!
   CHARACTER(6) :: chsubr
   INTEGER :: iamat , ibloc , iend , iiat , iiend , ijat , ikpw , ilccnx , ilpw , imark , ipos , ireal , long1 , long2 , mxccex ,  &
            & mxptbl , nbloc0 , nfree , ntb
   LOGICAL :: impr
   EXTERNAL blccn1 , ccnicopy , compos , csrcsc , dcopy , dperm , tbzero
!
! End of declarations rewritten by SPAG
!
!-----------------------------------------------------------------------
!
!     This routine determines if the matrix given by the structure
!     IA et JA is irreductible. If not, it orders the unknowns such
!     that all the consecutive unknowns in KPW between NSBLOC(i-1)+1
!     and NSBLOC(i) belong to the ith component of the matrix.
!     The numerical values of the matrix are in AMAT. They are modified
!     only if JOB = 1 and if we have more than one connected component.
!
!     On entry:
!     --------
!     n      = row and column dimension of the matrix
!     nblcmx = maximum number of connected components allowed. The size
!              of NSBLOC is nblcmx + 1 (in fact, it starts at 0).
!     job    = integer indicating the work to be done:
!              job = 1  if the permutation LPW is modified, we
!                       permute not only the structure of the matrix
!                       but also its numerical values.
!              job.ne.1 if the permutation LPW is modified, we permute
!                       the structure of the matrix ignoring real values.
!     iout   = impression parameter. If 0 < iout < 100, we print
!              comments and error messages on unit IOUT.
!     nw     = length of the work vector IZS.
!
!     Input / output:
!     --------------
!     nbloc  = number of connected components of the matrix. If the
!              matrix is not irreductible, nbloc > 1. We allow
!              nbloc > 1 on entry; in this case we calculate the
!              number of connected components in each previous one.
!     nsbloc = integer array of length NBLOC + 1 containing the pointers
!              to the first node of each component on the old (input)
!              and on the new (output) ordering.
!     lpw    = integer array of length N corresponding to the
!              permutation of the unknowns. We allow LPW to be a vector
!              different from the identity on input.
!     amat   = real*8 values of the matrix given by the structure IA, JA.
!     ja     = integer array of length NNZERO (= IA(N+1)-IA(1)) corresponding
!              to the column indices of nonzero elements of the matrix, stored
!              rowwise. It is modified only if the matrix has more
!              than one connected component.
!     ia     = integer array of length N+1 corresponding to the
!              pointer to the beginning of each row in JA (compressed
!              sparse row storage). It is modified only if
!              the matrix has more than one connected component.
!
!     On return:
!     ----------
!     ier    = integer. Error message. Normal return ier = 0.
!
!     Work space:
!     ----------
!     izs    = integer vector of length NW
!
!-----------------------------------------------------------------------
!     Laura C. Dutto - email: dutto@cerca.umontreal.ca
!                      July 1992 - Update: March 1994
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
   Ier = 0
   impr = Iout>0 .AND. Iout<=99
   ntb = Ia(N+1) - 1
   mxccex = max(Nblcmx,20)
!.....The matrix AMAT is a real*8 vector
   ireal = 2
!
!.....MXPTBL: maximal number of vertices by block
   mxptbl = 0
   DO ibloc = 1 , Nbloc
      mxptbl = max(mxptbl,Nsbloc(ibloc)-Nsbloc(ibloc-1))
   ENDDO
!
   long1 = Nbloc*mxptbl
   long2 = Nbloc*(mxccex+1)
!.....Dynamic allocation of memory
   iend = 1
   iiend = iend
   ilpw = iiend
   ikpw = ilpw + N
   ilccnx = ikpw + long1
   imark = ilccnx + long2
   iend = imark + N
   IF ( iend>Nw ) THEN
      CALL spag_block_2
      RETURN
   ENDIF
!
   nbloc0 = Nbloc
   chsubr = 'BLCCN1'
!.....We determine if the matrix has more than NBLOC0 connected components.
   CALL blccn1(N,Nbloc,Nblcmx,Nsbloc,Izs(ilpw),Izs(ikpw),Ia,Ja,Izs(imark),mxccex,Izs(ilccnx),mxptbl,Iout,Ier)
   IF ( Ier/=0 ) THEN
!
      IF ( impr ) WRITE (Iout,99001) chsubr , Ier
!
99001 FORMAT (' ***BLCCNX*** ERROR IN ',a6,'. IER = ',i8)
   ELSE
!
      IF ( Nbloc>nbloc0 ) THEN
!..........The matrix has more than NBLOC0 conneted components. So, we
!..........modify the vectors IA and JA to take account of the new permutation.
         nfree = iend - ikpw
         CALL tbzero(Izs(ikpw),nfree)
         iiat = ikpw
         ijat = iiat + N + 1
         iamat = ijat + ntb
         iend = iamat
         IF ( Job==1 ) iend = iamat + ireal*ntb
         IF ( iend>Nw ) THEN
            CALL spag_block_2
            RETURN
         ENDIF
!
!..........We copy IA and JA on IAT and JAT respectively
         CALL ccnicopy(N+1,Ia,Izs(iiat))
         CALL ccnicopy(ntb,Ja,Izs(ijat))
         IF ( Job==1 ) CALL dcopy(ntb,Amat,1,Izs(iamat),1)
         CALL dperm(N,Izs(iamat),Izs(ijat),Izs(iiat),Amat,Ja,Ia,Izs(ilpw),Izs(ilpw),Job)
         ipos = 1
!..........We sort columns inside JA.
         CALL csrcsc(N,Job,ipos,Amat,Ja,Ia,Izs(iamat),Izs(ijat),Izs(iiat))
!         CALL csrcsc(N,Job,ipos,Izs(iamat),Izs(ijat),Izs(iiat),Amat,Ja,Ia)
      ENDIF
!.....We modify the ordering of unknowns in LPW
      CALL compos(N,Lpw,Izs(ilpw))
   ENDIF
   CALL spag_block_1
CONTAINS
   SUBROUTINE spag_block_1
!
      nfree = iend - iiend
      CALL tbzero(Izs(iiend),nfree)
      iend = iiend
      RETURN
   END SUBROUTINE spag_block_1
   SUBROUTINE spag_block_2
      IF ( impr ) WRITE (Iout,99002) Nw , iend
99002 FORMAT (' ***BLCCNX*** THERE IS NOT ENOUGH MEMORY IN THE WORK',' VECTOR.'/13X,' ALLOWED MEMORY = ',I10,'  - NEEDED',         &
             &' MEMORY = ',I10)
      IF ( Ier==0 ) Ier = -1
      CALL spag_block_1
      RETURN
   END SUBROUTINE spag_block_2
END SUBROUTINE blccnx
!*==cconex.f90 processed by SPAG 8.04RA 15:47 21 Feb 2024
!!SPAG Open source Personal, Educational or Academic User Clemson University  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
! **********************************************************************
SUBROUTINE cconex(N,Icol0,Mxccnx,Lccnex,Kpw,Ia,Ja,Mark,Iout,Ier)
   IMPLICIT NONE
!
! Dummy argument declarations rewritten by SPAG
!
   INTEGER :: N
   INTEGER , INTENT(INOUT) :: Mxccnx
   INTEGER , INTENT(IN) :: Icol0
   INTEGER , INTENT(INOUT) , DIMENSION(0:Mxccnx) :: Lccnex
   INTEGER , INTENT(INOUT) , DIMENSION(N) :: Kpw
   INTEGER , INTENT(IN) , DIMENSION(N+1) :: Ia
   INTEGER , INTENT(IN) , DIMENSION(*) :: Ja
   INTEGER , INTENT(INOUT) , DIMENSION(N) :: Mark
   INTEGER , INTENT(IN) :: Iout
   INTEGER , INTENT(OUT) :: Ier
!
! Local variable declarations rewritten by SPAG
!
   INTEGER :: i , iccnex , ideb , idiag , index , ipos , ir , j , jref , mccnex , nccnex , newind , noeicc , nwindx
   LOGICAL :: impr
   EXTERNAL anccnx , newcnx , numini , tbzero
!
! End of declarations rewritten by SPAG
!
!-----------------------------------------------------------------------
!
!     This routine determines if the matrix given by the structure
!     IA and JA is irreductible. If not, it orders the unknowns such
!     that all the consecutive unknowns in KPW between LCCNEX(i-1)+1
!     and LCCNEX(i) belong to the ith component of the matrix.
!     The structure of the matrix could be nonsymmetric.
!     The diagonal vertices (if any) will belong to the last connected
!     component (convention).
!
!     On entry:
!     --------
!     n      = row and column dimension of the matrix
!     icol0  = the columns of the matrix are between ICOL0+1 and ICOL0+N
!     iout   = impression parameter. If 0 < IOUT < 100, we print
!              comments and error messages on unit IOUT.
!     ia     = integer array of length N+1 corresponding to the
!              pointer to the beginning of each row in JA (compressed
!              sparse row storage).
!     ja     = integer array of length NNZERO (= IA(N+1)-IA(1))
!              corresponding to the column indices of nonzero elements
!              of the matrix, stored rowwise.
!
!     Input/Output:
!     ------------
!     mxccnx = maximum number of connected components allowed on input,
!              and number of connected components of the matrix, on output.
!
!     On return:
!     ----------
!     lccnex = integer array of length MXCCNX + 1 containing the pointers
!              to the first node of each component, in the vector KPW.
!     kpw    = integer array of length N corresponding to the
!              inverse of permutation vector.
!     ier    = integer. Error message. Normal return ier = 0.
!
!     Work space:
!     ----------
!     mark   = integer vector of length N
!
!-----------------------------------------------------------------------
!     Laura C. Dutto - email: dutto@cerca.umontreal.ca
!                      July 1992 - Update: March 1994
!-----------------------------------------------------------------------
   INTEGER :: spag_nextblock_1
   Ier = 0
   ipos = Ia(1) - 1
   impr = Iout>0 .AND. Iout<=99
   spag_nextblock_1 = 1
   SPAG_DispatchLoop_1: DO
      SELECT CASE (spag_nextblock_1)
      CASE (1)
!-----------------------------------------------------------------------
         Ier = 0
         ipos = Ia(1) - 1
         impr = Iout>0 .AND. Iout<=99
!
         nccnex = 0
!.....We initialize MARK to zero. At the end of the algorithm, it would
!.....indicate the number of connected component associated with the vertex.
!.....The number (-1) indicates that the row associated with this vertex
!.....is a diagonal row. This value could be modified because we accept
!.....a non symmetric matrix. All the diagonal vertices will be put in
!.....the same connected component.
         CALL tbzero(Mark,N)
         spag_nextblock_1 = 2
      CASE (2)
!
         DO i = 1 , N
            IF ( Mark(i)==0 ) THEN
               ideb = i
               spag_nextblock_1 = 3
               CYCLE SPAG_DispatchLoop_1
            ENDIF
         ENDDO
!.......................................................................
!
!     We have partitioned the graph in its connected components!
!
!.......................................................................
!
!.....All the vertices have been already marked. Before modifying KPW
!.....(if necessary), we put the diagonal vertex (if any) in the last
!.....connected component.
         CALL tbzero(Lccnex(1),nccnex)
!
         idiag = 0
         DO i = 1 , N
            iccnex = Mark(i)
            IF ( iccnex==-1 ) THEN
               idiag = idiag + 1
               IF ( idiag==1 ) THEN
                  nccnex = nccnex + 1
                  IF ( nccnex>Mxccnx ) THEN
                     spag_nextblock_1 = 5
                     CYCLE SPAG_DispatchLoop_1
                  ENDIF
                  IF ( impr ) WRITE (Iout,99003)
! 323  format(' ***CCONEX*** ERROR IN ',A6,'. IER = ',I8)
99003             FORMAT (/' ***CCONEX*** THE LAST CONNECTED COMPONENT WILL',' HAVE THE DIAGONAL VERTICES.')
               ENDIF
               Mark(i) = nccnex
            ELSE
               Lccnex(iccnex) = Lccnex(iccnex) + 1
            ENDIF
         ENDDO
         IF ( idiag>=1 ) Lccnex(nccnex) = idiag
!
         IF ( nccnex==1 ) THEN
            Lccnex(nccnex) = N
            spag_nextblock_1 = 4
            CYCLE SPAG_DispatchLoop_1
         ENDIF
!
         iccnex = 1
         DO WHILE ( iccnex<=nccnex )
            IF ( Lccnex(iccnex)<=0 ) THEN
               DO i = 1 , N
                  IF ( Mark(i)>=iccnex ) Mark(i) = Mark(i) - 1
               ENDDO
               nccnex = nccnex - 1
               DO mccnex = iccnex , nccnex
                  Lccnex(mccnex) = Lccnex(mccnex+1)
               ENDDO
            ELSE
               iccnex = iccnex + 1
            ENDIF
         ENDDO
!
         index = 0
         DO iccnex = 1 , nccnex
            noeicc = Lccnex(iccnex)
            Lccnex(iccnex) = index
            index = index + noeicc
         ENDDO
         IF ( index/=N ) THEN
!
            IF ( impr ) WRITE (Iout,99001) index , N
!
99001       FORMAT (' ***CCONEX*** ERROR TRYING TO DETERMINE THE NUMBER',' OF CONNECTED COMPONENTS.'/13X,' NUMBER OF MARKED',      &
                   &' VERTICES =',i7,3x,'TOTAL NUMBER OF VERTICES =',I7)
            spag_nextblock_1 = 6
            CYCLE SPAG_DispatchLoop_1
         ELSE
!
!.....We define correctly KPW
            DO i = 1 , N
               iccnex = Mark(i)
               index = Lccnex(iccnex) + 1
               Kpw(index) = i
               Lccnex(iccnex) = index
            ENDDO
            spag_nextblock_1 = 4
            CYCLE SPAG_DispatchLoop_1
         ENDIF
      CASE (3)
!
         IF ( Ia(ideb+1)-Ia(ideb)==1 ) THEN
!..........The row is a diagonal row.
            Mark(ideb) = -1
            spag_nextblock_1 = 2
            CYCLE SPAG_DispatchLoop_1
         ENDIF
         iccnex = nccnex + 1
         IF ( iccnex>Mxccnx ) THEN
            spag_nextblock_1 = 5
            CYCLE SPAG_DispatchLoop_1
         ENDIF
         index = 0
         newind = 0
         jref = 0
         Mark(ideb) = iccnex
         index = index + 1
         Kpw(index) = ideb
         DO
!
            jref = jref + 1
            ideb = Kpw(jref)
 
            DO ir = Ia(ideb) - ipos , Ia(ideb+1) - ipos - 1
               j = Ja(ir) - Icol0
               mccnex = Mark(j)
               IF ( mccnex<=0 ) THEN
                  index = index + 1
                  Kpw(index) = j
                  Mark(j) = iccnex
               ELSEIF ( mccnex==iccnex ) THEN
               ELSEIF ( mccnex>iccnex ) THEN
!.............We realize that the connected component MCCNX is,
!.............in fact, included in this one. We modify MARK and KPW.
                  CALL newcnx(N,mccnex,iccnex,index,Kpw,Mark)
                  IF ( mccnex==nccnex ) nccnex = nccnex - 1
               ELSE
!.............We realize that the previously marked vertices belong,
!.............in fact, to the connected component ICCNX. We modify MARK.
                  CALL anccnx(N,iccnex,mccnex,Mark,nwindx)
                  iccnex = mccnex
                  newind = newind + nwindx
               ENDIF
            ENDDO
            IF ( jref>=index ) THEN
!
!.....We have finished with this connected component.
               index = index + newind
               IF ( iccnex==nccnex+1 ) nccnex = nccnex + 1
               spag_nextblock_1 = 2
               CYCLE SPAG_DispatchLoop_1
            ENDIF
         ENDDO
         spag_nextblock_1 = 4
      CASE (4)
!
         Mxccnx = nccnex
         Lccnex(0) = nccnex
         IF ( nccnex==1 ) CALL numini(N,Kpw)
         RETURN
      CASE (5)
         IF ( impr ) WRITE (Iout,99002) nccnex , Mxccnx
99002    FORMAT (' ***CCONEX*** THE ALLOWED NUMBER OF CONNECTED COMPONENTS',' IS NOT ENOUGH.'/13X,' NECESSARY NUMBER = ',I4,5x,    &
                &' ALLOWED NUMBER = ',I4)
         spag_nextblock_1 = 6
      CASE (6)
         Ier = -1
         RETURN
      END SELECT
   ENDDO SPAG_DispatchLoop_1
END SUBROUTINE cconex
!*==anccnx.f90 processed by SPAG 8.04RA 15:47 21 Feb 2024
!!SPAG Open source Personal, Educational or Academic User Clemson University  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
! **********************************************************************
SUBROUTINE anccnx(N,Mccnex,Iccnex,Mark,Ncount)
   IMPLICIT NONE
!
! Dummy argument declarations rewritten by SPAG
!
   INTEGER , INTENT(IN) :: N
   INTEGER , INTENT(IN) :: Mccnex
   INTEGER , INTENT(IN) :: Iccnex
   INTEGER , INTENT(INOUT) , DIMENSION(N) :: Mark
   INTEGER , INTENT(INOUT) :: Ncount
!
! Local variable declarations rewritten by SPAG
!
   INTEGER :: i
!
! End of declarations rewritten by SPAG
!
!-----------------------------------------------------------------------
!
!     We put in ICCNEX the vertices marked in the component MCCNEX.
!
!-----------------------------------------------------------------------
!     include "NSIMPLIC"
!-----------------------------------------------------------------------
!     Laura C. Dutto - email: dutto@cerca.umontreal.ca - December 1993
!-----------------------------------------------------------------------
   Ncount = 0
   DO i = 1 , N
      IF ( Mark(i)==Mccnex ) THEN
         Mark(i) = Iccnex
         Ncount = Ncount + 1
      ENDIF
   ENDDO
!
END SUBROUTINE anccnx
!*==newcnx.f90 processed by SPAG 8.04RA 15:47 21 Feb 2024
!!SPAG Open source Personal, Educational or Academic User Clemson University  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
! **********************************************************************
SUBROUTINE newcnx(N,Mccnex,Iccnex,Index,Kpw,Mark)
   IMPLICIT NONE
!
! Dummy argument declarations rewritten by SPAG
!
   INTEGER , INTENT(IN) :: N
   INTEGER , INTENT(IN) :: Mccnex
   INTEGER , INTENT(IN) :: Iccnex
   INTEGER , INTENT(INOUT) :: Index
   INTEGER , INTENT(OUT) , DIMENSION(*) :: Kpw
   INTEGER , INTENT(INOUT) , DIMENSION(N) :: Mark
!
! Local variable declarations rewritten by SPAG
!
   INTEGER :: i
!
! End of declarations rewritten by SPAG
!
!-----------------------------------------------------------------------
!
!     We put in ICCNEX the vertices marked in the component MCCNEX. We
!     modify also the vector KPW.
!
!-----------------------------------------------------------------------
!     include "NSIMPLIC"
!-----------------------------------------------------------------------
!     Laura C. Dutto - email: dutto@cerca.umontreal.ca - December 1993
!-----------------------------------------------------------------------
   DO i = 1 , N
      IF ( Mark(i)==Mccnex ) THEN
         Mark(i) = Iccnex
         Index = Index + 1
         Kpw(Index) = i
      ENDIF
   ENDDO
!
END SUBROUTINE newcnx
!*==blccn1.f90 processed by SPAG 8.04RA 15:47 21 Feb 2024
!!SPAG Open source Personal, Educational or Academic User Clemson University  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
! **********************************************************************
SUBROUTINE blccn1(N,Nbloc,Nblcmx,Nsbloc,Lpw,Kpw,Ia,Ja,Mark,Mxccex,Lccnex,Mxptbl,Iout,Ier)
   IMPLICIT NONE
!
! Dummy argument declarations rewritten by SPAG
!
   INTEGER , INTENT(IN) :: N
   INTEGER , INTENT(INOUT) :: Nbloc
   INTEGER , INTENT(IN) :: Mxccex
   INTEGER , INTENT(IN) :: Mxptbl
   INTEGER , INTENT(IN) :: Nblcmx
   INTEGER , INTENT(INOUT) , DIMENSION(0:Nbloc) :: Nsbloc
   INTEGER , DIMENSION(N) :: Lpw
   INTEGER , DIMENSION(Mxptbl*Nbloc) :: Kpw
   INTEGER , DIMENSION(N+1) :: Ia
   INTEGER , DIMENSION(*) :: Ja
   INTEGER , DIMENSION(N) :: Mark
   INTEGER , DIMENSION((Mxccex+1)*Nbloc) :: Lccnex
   INTEGER , INTENT(IN) :: Iout
   INTEGER , INTENT(INOUT) :: Ier
!
! Local variable declarations rewritten by SPAG
!
   CHARACTER(6) :: chsubr
   INTEGER :: ibloc , icc , ik0 , ik1 , ilccnx , info , ins0 , isor , kpibl , lcc0 , nccnex , newblc , nsb , nsfin , ntb0
   LOGICAL :: impr
   EXTERNAL cconex , invlpw , iplusa , numini
!
! End of declarations rewritten by SPAG
!
!-----------------------------------------------------------------------
!
!     This routine determines if the matrix given by the structure
!     IA et JA is irreductible. If not, it orders the unknowns such
!     that all the consecutive unknowns in KPW between NSBLOC(i-1)+1
!     and NSBLOC(i) belong to the ith component of the matrix.
!
!     On entry:
!     --------
!     n      = row and column dimension of the matrix
!     nblcmx = The size of NSBLOC is nblcmx + 1 (in fact, it starts at 0).
!     ia     = integer array of length N+1 corresponding to the
!              pointer to the beginning of each row in JA (compressed
!              sparse row storage).
!     ja     = integer array of length NNZERO (= IA(N+1)-IA(1)) corresponding
!              to the column indices of nonzero elements of the matrix,
!              stored rowwise.
!     mxccex = maximum number of connected components allowed by block.
!     mxptbl = maximum number of points (or unknowns) in each connected
!              component (mxptbl .le. n).
!     iout   = impression parameter. If 0 < iout < 100, we print
!              comments and error messages on unit IOUT.
!
!     Input/Output:
!     ------------
!     nbloc  = number of connected components of the matrix. If the
!              matrix is not irreductible, nbloc > 1. We allow
!              nbloc > 1 on entry; in this case we calculate the
!              number of connected components in each previous one.
!     nsbloc = integer array of length NBLOC + 1 containing the pointers
!              to the first node of each component on the new ordering.
!              Normally, on entry you put: NBLOC = 1, NSBLOC(0) = 0,
!              NSBLOC(NBLOC) = N.
!
!     On return:
!     ----------
!     lpw    = integer array of length N corresponding to the
!              permutation vector (the row i goes to lpw(i)).
!     ier    = integer. Error message. Normal return ier = 0.
!
!     Work space:
!     ----------
!     kpw    = integer vector of length MXPTBL*NBLOC necessary for parallel
!              computation.
!     mark   = integer vector of length N
!     lccnex = integer vector of length (MXCCEX+1)*NBLOC necessary for parallel
!              computation.
!
!-----------------------------------------------------------------------
!     Laura C. Dutto - e-mail: dutto@cerca.umontreal.ca
!                      Juillet 1992. Update: March 1994
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
   Ier = 0
   nsb = 0
   impr = Iout>0 .AND. Iout<=99
   isor = 0
!
   chsubr = 'CCONEX'
   newblc = 0
!$DOACROSS if(nbloc.gt.1), LOCAL(ibloc, ik0, ik1, ins0, ntb0,
!$&   nccnex, ilccnx, info, kpibl), REDUCTION(ier, newblc)
   DO ibloc = 1 , Nbloc
      ik0 = Nsbloc(ibloc-1)
      ik1 = Nsbloc(ibloc)
      ntb0 = Ia(ik0+1)
      IF ( Ia(ik1+1)-ntb0>1 ) THEN
         ntb0 = ntb0 - 1
         ins0 = ik1 - ik0
!........We need more memory place for KPW1 because of parallel computation
         kpibl = (ibloc-1)*Mxptbl
         CALL numini(ins0,Kpw(kpibl+1))
         nccnex = Mxccex
         ilccnx = (Mxccex+1)*(ibloc-1) + 1
!.......................................................................
!
!        Call to the main routine: CCONEX
!
!.......................................................................
         CALL cconex(ins0,ik0,nccnex,Lccnex(ilccnx),Kpw(kpibl+1),Ia(ik0+1),Ja(ntb0+1),Mark(ik0+1),isor,info)
         Ier = Ier + info
         IF ( info==0 .AND. nccnex>=1 ) THEN
!
!........We add the new connected components on NEWBLC
            newblc = newblc + nccnex
!........We define LPW different from the identity only if there are more
!........than one connected component in this block
            IF ( nccnex==1 ) THEN
               CALL numini(ins0,Lpw(ik0+1))
            ELSE
               CALL invlpw(ins0,Kpw(kpibl+1),Lpw(ik0+1))
            ENDIF
            CALL iplusa(ins0,ik0,1,Lpw(ik0+1))
         ENDIF
      ENDIF
   ENDDO
!
   IF ( Ier/=0 ) THEN
!
      IF ( impr ) WRITE (Iout,99001) chsubr , Ier
!
99001 FORMAT (' ***BLCCN1*** ERROR IN ',a6,'. IER = ',i8)
   ELSEIF ( newblc/=Nbloc ) THEN
      IF ( newblc>Nblcmx ) THEN
         IF ( impr ) WRITE (Iout,99002) newblc , Nblcmx
99002    FORMAT (' ***BLCCN1*** THE MEMORY SPACE ALLOWED FOR NSBLOC IS',' NOT ENOUGH.'/13X,' NUMBER (NECESSARY) OF CONNECTED',     &
                &' COMPONENTS = ',I5/13X,' MAXIMAL NUMBER OF BLOCKS',14x,'= ',i5)
         IF ( Ier==0 ) Ier = -1
      ELSE
!
!.....We modify the number of blocks to indicate the number of connected
!.....components in the matrix.
         newblc = 0
         nsfin = 0
!DIR$ NEXT SCALAR
         DO ibloc = 1 , Nbloc
            ilccnx = (Mxccex+1)*(ibloc-1) + 1
            nccnex = Lccnex(ilccnx)
            IF ( nccnex>1 .AND. impr ) WRITE (Iout,99003) ibloc , nccnex
99003       FORMAT (' *** The block ',i3,' has ',i3,' strongly connected',' components. The number of vertices by component is:')
            lcc0 = 0
!DIR$ NEXT SCALAR
            DO icc = 1 , nccnex
               newblc = newblc + 1
               nsb = Lccnex(ilccnx+icc)
!...........Be careful! In LCCNEX we have the cumulated number of vertices
               Nsbloc(newblc) = nsfin + nsb
               IF ( nccnex>1 .AND. impr ) WRITE (Iout,99004) icc , nsb - lcc0
99004          FORMAT (5x,'Component No.',i3,' - Number of vertices = ',i6)
               lcc0 = nsb
            ENDDO
            nsfin = nsfin + nsb
         ENDDO
         Nbloc = newblc
      ENDIF
   ENDIF
!
   RETURN
END SUBROUTINE blccn1
!*==ccnicopy.f90 processed by SPAG 8.04RA 15:47 21 Feb 2024
!!SPAG Open source Personal, Educational or Academic User Clemson University  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
!***********************************************************************
SUBROUTINE ccnicopy(N,Ix,Iy)
   IMPLICIT NONE
!
! Dummy argument declarations rewritten by SPAG
!
   INTEGER , INTENT(IN) :: N
   INTEGER , INTENT(IN) , DIMENSION(N) :: Ix
   INTEGER , INTENT(OUT) , DIMENSION(N) :: Iy
!
! Local variable declarations rewritten by SPAG
!
   INTEGER :: i
!
! End of declarations rewritten by SPAG
!
!.......................................................................
!     We copy the vector IX on the vector IY
!.......................................................................
!.......................................................................
   IF ( N<=0 ) RETURN
!$DOACROSS if(n .gt. 250), local(i)
   DO i = 1 , N
      Iy(i) = Ix(i)
   ENDDO
!
END SUBROUTINE ccnicopy
!*==compos.f90 processed by SPAG 8.04RA 15:47 21 Feb 2024
!!SPAG Open source Personal, Educational or Academic User Clemson University  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
!***********************************************************************
SUBROUTINE compos(N,Lpw0,Lpw1)
   IMPLICIT NONE
!
! Dummy argument declarations rewritten by SPAG
!
   INTEGER , INTENT(IN) :: N
   INTEGER , INTENT(INOUT) , DIMENSION(N) :: Lpw0
   INTEGER , INTENT(IN) , DIMENSION(N) :: Lpw1
!
! Local variable declarations rewritten by SPAG
!
   INTEGER :: i0
!
! End of declarations rewritten by SPAG
!
!-----------------------------------------------------------------------
!
!     We take account of the original order of unknowns. We put the
!     final result on LPW0.
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!     Laura C. Dutto - Mars 1994
!-----------------------------------------------------------------------
!$DOACROSS if(n .gt. 250), local(i0)
   DO i0 = 1 , N
      Lpw0(i0) = Lpw1(Lpw0(i0))
   ENDDO
!
END SUBROUTINE compos
!*==invlpw.f90 processed by SPAG 8.04RA 15:47 21 Feb 2024
!!SPAG Open source Personal, Educational or Academic User Clemson University  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
! **********************************************************************
SUBROUTINE invlpw(N,Lpw,Kpw)
   IMPLICIT NONE
!
! Dummy argument declarations rewritten by SPAG
!
   INTEGER , INTENT(IN) :: N
   INTEGER , INTENT(IN) , DIMENSION(N) :: Lpw
   INTEGER , INTENT(OUT) , DIMENSION(N) :: Kpw
!
! Local variable declarations rewritten by SPAG
!
   INTEGER :: i0 , i1
!
! End of declarations rewritten by SPAG
!
!.......................................................................
!
!     KPW is the inverse of LPW
!
!.......................................................................
!.......................................................................
!     Laura C. Dutto - Novembre 1993
!.......................................................................
!$DOACROSS if(n .gt. 200), local(i0, i1)
   DO i0 = 1 , N
      i1 = Lpw(i0)
      Kpw(i1) = i0
   ENDDO
!
END SUBROUTINE invlpw
!*==numini.f90 processed by SPAG 8.04RA 15:47 21 Feb 2024
!!SPAG Open source Personal, Educational or Academic User Clemson University  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
! **********************************************************************
SUBROUTINE numini(N,Lpw)
   IMPLICIT NONE
!
! Dummy argument declarations rewritten by SPAG
!
   INTEGER , INTENT(IN) :: N
   INTEGER , INTENT(OUT) , DIMENSION(N) :: Lpw
!
! Local variable declarations rewritten by SPAG
!
   INTEGER :: i
!
! End of declarations rewritten by SPAG
!
!.......................................................................
!.......................................................................
!
!     The vector LPW is initialized as the identity.
!
!.......................................................................
!     Laura C. Dutto - Novembre 1993
!.......................................................................
!$DOACROSS if(n .gt. 250), local(i)
   DO i = 1 , N
      Lpw(i) = i
   ENDDO
!
END SUBROUTINE numini
!*==tbzero.f90 processed by SPAG 8.04RA 15:47 21 Feb 2024
!!SPAG Open source Personal, Educational or Academic User Clemson University  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
!***********************************************************************
SUBROUTINE tbzero(M,Nmot)
   IMPLICIT NONE
!
! Dummy argument declarations rewritten by SPAG
!
   INTEGER , INTENT(IN) :: Nmot
   INTEGER , INTENT(OUT) , DIMENSION(Nmot) :: M
!
! Local variable declarations rewritten by SPAG
!
   INTEGER :: i
!
! End of declarations rewritten by SPAG
!
!.......................................................................
!     We initialize to ZERO an integer vector of length NMOT.
!.......................................................................
!.......................................................................
   IF ( Nmot<=0 ) RETURN
!$DOACROSS if(nmot.gt.500), LOCAL(i)
   DO i = 1 , Nmot
      M(i) = 0
   ENDDO
END SUBROUTINE tbzero
!*==iplusa.f90 processed by SPAG 8.04RA 15:47 21 Feb 2024
!!SPAG Open source Personal, Educational or Academic User Clemson University  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
! **********************************************************************
SUBROUTINE iplusa(N,Nalpha,Nbeta,Ia)
   IMPLICIT NONE
!
! Dummy argument declarations rewritten by SPAG
!
   INTEGER , INTENT(IN) :: N
   INTEGER , INTENT(IN) :: Nalpha
   INTEGER , INTENT(IN) :: Nbeta
   INTEGER , INTENT(INOUT) , DIMENSION(N) :: Ia
!
! Local variable declarations rewritten by SPAG
!
   INTEGER :: i , nmax
!
! End of declarations rewritten by SPAG
!
!.......................................................................
!
!     We add NALPHA to each element of NBETA * IA:
!
!            ia(i) = nalpha + nbeta * ia(i)
!
!.......................................................................
!.......................................................................
!     Laura C. Dutto - February 1994
!.......................................................................
   IF ( N<=0 ) RETURN
!
   nmax = 500
   IF ( Nalpha==0 ) THEN
      IF ( Nbeta==1 ) RETURN
      IF ( Nbeta==-1 ) THEN
!$DOACROSS if(n .gt. nmax), local (i)
         DO i = 1 , N
            Ia(i) = -Ia(i)
         ENDDO
      ELSE
!$DOACROSS if(n .gt. nmax/2), local (i)
         DO i = 1 , N
            Ia(i) = Nbeta*Ia(i)
         ENDDO
      ENDIF
      RETURN
   ENDIF
   IF ( Nbeta==0 ) THEN
!$DOACROSS if(n .gt. nmax), local (i)
      DO i = 1 , N
         Ia(i) = Nalpha
      ENDDO
      RETURN
   ENDIF
   IF ( Nbeta==-1 ) THEN
!$DOACROSS if(n .gt. nmax/2), local (i)
      DO i = 1 , N
         Ia(i) = Nalpha - Ia(i)
      ENDDO
   ELSEIF ( Nbeta==1 ) THEN
!$DOACROSS if(n .gt. nmax/2), local (i)
      DO i = 1 , N
         Ia(i) = Nalpha + Ia(i)
      ENDDO
   ELSE
!$DOACROSS if(n .gt. nmax/3), local (i)
      DO i = 1 , N
         Ia(i) = Nalpha + Nbeta*Ia(i)
      ENDDO
   ENDIF
!
END SUBROUTINE iplusa
