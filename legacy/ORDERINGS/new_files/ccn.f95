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

subroutine blccnx(n, nbloc, nblcmx, nsbloc, job, lpw, amat, ja, ia, iout, ier, izs, nw)
   integer, intent(In) :: n, nblcmx, job, nw
   integer, intent(InOut) :: nbloc
   integer, intent(Out) :: ier
   integer :: ntb, iout, mxccex, ireal, mxptbl, ibloc, long1, long2, iend, iiend, ilpw, ikpw, ilccnx, &
              imark, nbloc0, nfree, iiat, iamat, ipos
   integer, dimension(nw), intent(InOut) :: izs
   integer, dimension(n), intent(InOut) :: lpw
   integer, dimension(nblcmx), intent(InOut) :: nsbloc
   integer, dimension(n + 1), intent(In) :: ia
   integer, dimension(*), intent(In) :: ja
   real(kind=8), dimension(1:*), intent(In) :: amat
   logical :: impr
   character(len=6) :: chsubr
! -----------------------------------------------------------------------
!
!      This routine determines if the matrix given by the structure
!      IA et JA is irreductible. If not, it orders the unknowns such
!      that all the consecutive unknowns in KPW between NSBLOC(i-1)+1
!      and NSBLOC(i) belong to the ith component of the matrix.
!      The numerical values of the matrix are in AMAT. They are modified
!      only if JOB = 1 and if we have more than one connected component.
!
!      On entry:
!      --------
!      n      = row and column dimension of the matrix
!      nblcmx = maximum number of connected components allowed. The size
!               of NSBLOC is nblcmx + 1 (in fact, it starts at 0).
!      job    = integer indicating the work to be done:
!               job = 1  if the permutation LPW is modified, we
!                        permute not only the structure of the matrix
!                        but also its numerical values.
!               job.ne.1 if the permutation LPW is modified, we permute
!                        the structure of the matrix ignoring real values.
!      iout   = impression parameter. If 0 < iout < 100, we print
!               comments and error messages on unit IOUT.
!      nw     = length of the work vector IZS.
!
!      Input / output:
!      --------------
!      nbloc  = number of connected components of the matrix. If the
!               matrix is not irreductible, nbloc > 1. We allow
!               nbloc > 1 on entry; in this case we calculate the
!               number of connected components in each previous one.
!      nsbloc = integer array of length NBLOC + 1 containing the pointers
!               to the first node of each component on the old (input)
!               and on the new (output) ordering.
!      lpw    = integer array of length N corresponding to the
!               permutation of the unknowns. We allow LPW to be a vector
!               different from the identity on input.
!      amat   = real*8 values of the matrix given by the structure IA, JA.
!      ja     = integer array of length NNZERO (= IA(N+1)-IA(1)) corresponding
!               to the column indices of nonzero elements of the matrix, stored
!               rowwise. It is modified only if the matrix has more
!               than one connected component.
!      ia     = integer array of length N+1 corresponding to the
!               pointer to the beginning of each row in JA (compressed
!               sparse row storage). It is modified only if
!               the matrix has more than one connected component.
!
!      On return:
!      ----------
!      ier    = integer. Error message. Normal return ier = 0.
!
!      Work space:
!      ----------
!      izs    = integer vector of length NW
!
! -----------------------------------------------------------------------
!      Laura C. Dutto - email: dutto@cerca.umontreal.ca
!                       July 1992 - Update: March 1994
! -----------------------------------------------------------------------

! -----------------------------------------------------------------------

   ier = 0
   impr = iout > 0 .and. iout <= 99
   ntb = ia(n + 1) - 1
   mxccex = max(nblcmx, 20)
! .....The matrix AMAT is a real*8 vector
   ireal = 2
!
! .....MXPTBL: maximal number of vertices by block
   mxptbl = 0
   do ibloc = 1, nbloc
      mxptbl = max(mxptbl, nsbloc(ibloc) - nsbloc(ibloc - 1))
   end do
!
   long1 = nbloc*mxptbl
   long2 = nbloc*(mxccex + 1)
! .....Dynamic allocation of memory
   iend = 1
   iiend = iend
   ilpw = iiend
   ikpw = ilpw + n
   ilccnx = ikpw + long1
   imark = ilccnx + long2
   iend = imark + n
   if (iend > nw) then
      if (impr) write (iout, 320) nw, iend
      if (ier == 0) ier = -1
      nfree = iend - iiend
      call tbzero(izs(iiend), nfree)
      iend = iiend
      return
   end if

!
   nbloc0 = nbloc
   chsubr = 'BLCCN1'
! .....We determine if the matrix has more than NBLOC0 connected components.
   call blccn1(n, nbloc, nblcmx, nsbloc, izs(ilpw), izs(ikpw), ia, ja, izs(imark), mxccex, izs(ilccnx), mxptbl, &
               iout, ier)
   if (ier /= 0) then
      if (impr) write (iout, 310) chsubr, ier
      nfree = iend - iiend
      call tbzero(izs(iiend), nfree)
      iend = iiend
      return
   end if
!
   if (nbloc > nbloc0) then
! ..........The matrix has more than NBLOC0 conneted components. So, we
! ..........modify the vectors IA and JA to take account of the new permutation.
      nfree = iend - ikpw
      call tbzero(izs(ikpw), nfree)
      iiat = ikpw
      ijat = iiat + n + 1
      iamat = ijat + ntb
      iend = iamat
      if (job == 1) iend = iamat + ireal*ntb
      if (iend > nw) then
         IF (IMPR) WRITE (IOUT, 320) nw, iend
         if (ier .eq. 0) ier = -1
         nfree = iend - iiend
         call tbzero(izs(iiend), nfree)
         iend = iiend
         return

      end if
!
! ..........We copy IA and JA on IAT and JAT respectively
      call ccnicopy(n + 1, ia, izs(iiat))
      call ccnicopy(ntb, ja, izs(ijat))
      if (job == 1) call dcopy(ntb, amat, 1, izs(iamat), 1)
      call dperm(n, izs(iamat), izs(ijat), izs(iiat), amat, ja, ia, izs(ilpw), izs(ilpw), job)
      ipos = 1
! ..........We sort columns inside JA.
      call csrcsc(n, job, ipos, amat, ja, ia, izs(iamat), izs(ijat), izs(iiat))
      call csrcsc(n, job, ipos, izs(iamat), izs(ijat), izs(iiat), amat, ja, ia)
   end if
! .....We modify the ordering of unknowns in LPW
   call compos(n, lpw, izs(ilpw))
!
   nfree = iend - iiend
   call tbzero(izs(iiend), nfree)
   iend = iiend
   return
!
310 format(' ***BLCCNX*** ERROR IN ', a6, '. IER = ', i8)
320 format(' ***BLCCNX*** THERE IS NOT ENOUGH MEMORY IN THE WORK', ' VECTOR.'/13x, &
          ' ALLOWED MEMORY = ', i10, '  - NEEDED', ' MEMORY = ', i10)
end subroutine blccnx

subroutine cconex(n, icol0, mxccnx, lccnex, kpw, ia, ja, mark, iout, ier)

   integer, intent(In) :: n, icol0
   integer, intent(InOut) :: mxccnx
   integer, intent(Out) :: ier
   integer :: ipos, iout, nccnex, i, ideb, iccnex, index, newind, jref, ir, ipo, j, mccnex, &
              nwindx, idiag, noeicc
   integer, dimension(n + 1), intent(In) :: ia
   integer, dimension(0:mxccnx), intent(InOut) :: lccnex
   integer, dimension(n), intent(Out) :: kpw
   integer, dimension(*), intent(In) :: ja
   integer, dimension(n), intent(InOut) :: mark
   logical :: impr
! -----------------------------------------------------------------------
!
!      This routine determines if the matrix given by the structure
!      IA and JA is irreductible. If not, it orders the unknowns such
!      that all the consecutive unknowns in KPW between LCCNEX(i-1)+1
!      and LCCNEX(i) belong to the ith component of the matrix.
!      The structure of the matrix could be nonsymmetric.
!      The diagonal vertices (if any) will belong to the last connected
!      component (convention).
!
!      On entry:
!      --------
!      n      = row and column dimension of the matrix
!      icol0  = the columns of the matrix are between ICOL0+1 and ICOL0+N
!      iout   = impression parameter. If 0 < IOUT < 100, we print
!               comments and error messages on unit IOUT.
!      ia     = integer array of length N+1 corresponding to the
!               pointer to the beginning of each row in JA (compressed
!               sparse row storage).
!      ja     = integer array of length NNZERO (= IA(N+1)-IA(1))
!               corresponding to the column indices of nonzero elements
!               of the matrix, stored rowwise.
!
!      Input/Output:
!      ------------
!      mxccnx = maximum number of connected components allowed on input,
!               and number of connected components of the matrix, on output.
!
!      On return:
!      ----------
!      lccnex = integer array of length MXCCNX + 1 containing the pointers
!               to the first node of each component, in the vector KPW.
!      kpw    = integer array of length N corresponding to the
!               inverse of permutation vector.
!      ier    = integer. Error message. Normal return ier = 0.
!
!      Work space:
!      ----------
!      mark   = integer vector of length N
!
! -----------------------------------------------------------------------
!      Laura C. Dutto - email: dutto@cerca.umontreal.ca
!                       July 1992 - Update: March 1994
! -----------------------------------------------------------------------

! -----------------------------------------------------------------------
   ier = 0
   ipos = ia(1) - 1
   impr = iout > 0 .and. iout <= 99
!
   nccnex = 0
! .....We initialize MARK to zero. At the end of the algorithm, it would
! .....indicate the number of connected component associated with the vertex.
! .....The number (-1) indicates that the row associated with this vertex
! .....is a diagonal row. This value could be modified because we accept
! .....a non symmetric matrix. All the diagonal vertices will be put in
! .....the same connected component.
   call tbzero(mark, n)
!
5  do i = 1, n
      if (mark(i) == 0) then
         ideb = i
         goto 15
      end if
   end do
   goto 35
!
15 if (ia(ideb + 1) - ia(ideb) == 1) then
! ..........The row is a diagonal row.
      mark(ideb) = -1
      goto 5
   end if
   iccnex = nccnex + 1
   if (iccnex > mxccnx) goto 220
   index = 0
   newind = 0
   jref = 0
   mark(ideb) = iccnex
   index = index + 1
   kpw(index) = ideb
!
20 jref = jref + 1
   ideb = kpw(jref)
   do ir = ia(ideb) - ipos, ia(ideb + 1) - ipos - 1
      j = ja(ir) - icol0
      mccnex = mark(j)
      if (mccnex <= 0) then
         index = index + 1
         kpw(index) = j
         mark(j) = iccnex
      else if (mccnex == iccnex) then
         goto 30
      else if (mccnex > iccnex) then
! .............We realize that the connected component MCCNX is,
! .............in fact, included in this one. We modify MARK and KPW.
         call newcnx(n, mccnex, iccnex, index, kpw, mark)
         if (mccnex == nccnex) nccnex = nccnex - 1
      else
! .............We realize that the previously marked vertices belong,
! .............in fact, to the connected component ICCNX. We modify MARK.
         call anccnx(n, iccnex, mccnex, mark, nwindx)
         iccnex = mccnex
         newind = newind + nwindx
      end if
30    continue
   end do
   if (jref < index) goto 20
!
! .....We have finished with this connected component.
   index = index + newind
   if (iccnex == nccnex + 1) nccnex = nccnex + 1
   goto 5
! .......................................................................
!
!      We have partitioned the graph in its connected components!
!
! .......................................................................
35 continue
!
! .....All the vertices have been already marked. Before modifying KPW
! .....(if necessary), we put the diagonal vertex (if any) in the last
! .....connected component.
   call tbzero(lccnex(1), nccnex)
!
   idiag = 0
   do i = 1, n
      iccnex = mark(i)
      if (iccnex == -1) then
         idiag = idiag + 1
         if (idiag == 1) then
            nccnex = nccnex + 1
            if (nccnex > mxccnx) goto 220
            if (impr) write (iout, 340)
         end if
         mark(i) = nccnex
      else
         lccnex(iccnex) = lccnex(iccnex) + 1
      end if
   end do
   if (idiag >= 1) lccnex(nccnex) = idiag
!
   if (nccnex == 1) then
      lccnex(nccnex) = n
      goto 40
   end if
!
   iccnex = 1
8  if (iccnex > nccnex) goto 12
   if (lccnex(iccnex) <= 0) then
   do i = 1, n
      if (mark(i) >= iccnex) mark(i) = mark(i) - 1
   end do
   nccnex = nccnex - 1
   do mccnex = iccnex, nccnex
      lccnex(mccnex) = lccnex(mccnex + 1)
   end do
   else
   iccnex = iccnex + 1
   end if
   goto 8
!
12 index = 0
   do iccnex = 1, nccnex
      noeicc = lccnex(iccnex)
      lccnex(iccnex) = index
      index = index + noeicc
   end do
   if (index /= n) goto 210
!
! .....We define correctly KPW
   do i = 1, n
      iccnex = mark(i)
      index = lccnex(iccnex) + 1
      kpw(index) = i
      lccnex(iccnex) = index
   end do
!
40 mxccnx = nccnex
   lccnex(0) = nccnex
   call numini(n, kpw)
   return
!
210 if (impr) write (iout, 310) index, n
   goto 235
220 if (impr) write (iout, 320) nccnex, mxccnx
   goto 235
235 ier = -1
   return
!
310 format(' ***CCONEX*** ERROR TRYING TO DETERMINE THE NUMBER', &
          ' OF CONNECTED COMPONENTS.'/13x, ' NUMBER OF MARKED', ' VERTICES =', i7, 3x, &
          'TOTAL NUMBER OF VERTICES =', i7)
320 format(' ***CCONEX*** THE ALLOWED NUMBER OF CONNECTED COMPONENTS', &
          ' IS NOT ENOUGH.'/13x, ' NECESSARY NUMBER = ', i4, 5x, ' ALLOWED NUMBER = ', i4)
!  323  format(' ***CCONEX*** ERROR IN ',A6,'. IER = ',I8)
340 format(/' ***CCONEX*** THE LAST CONNECTED COMPONENT WILL', &
           ' HAVE THE DIAGONAL VERTICES.')
end subroutine cconex
subroutine anccnx(n, mccnex, iccnex, mark, ncount)
! -----------------------------------------------------------------------
!
!      We put in ICCNEX the vertices marked in the component MCCNEX.
!
! -----------------------------------------------------------------------
!      include "NSIMPLIC"
! BEGIN new declarations
   implicit none
   integer, intent(In) :: n
   integer, intent(In) :: mccnex
   integer, intent(In) :: iccnex
   integer, intent(Out) :: ncount
   integer :: i
! END new declarations
   integer, dimension(1:n), intent(InOut) :: mark
! -----------------------------------------------------------------------
!      Laura C. Dutto - email: dutto@cerca.umontreal.ca - December 1993
! -----------------------------------------------------------------------
   ncount = 0
   do i = 1, n
   if (mark(i) == mccnex) then
      mark(i) = iccnex
      ncount = ncount + 1
   end if
   end do
!
   return
end subroutine anccnx
subroutine newcnx(n, mccnex, iccnex, index, kpw, mark)
! -----------------------------------------------------------------------
!
!      We put in ICCNEX the vertices marked in the component MCCNEX. We
!      modify also the vector KPW.
!
! -----------------------------------------------------------------------
!      include "NSIMPLIC"
! BEGIN new declarations
   implicit none
   integer, intent(In) :: n
   integer, intent(In) :: mccnex
   integer, intent(In) :: iccnex
   integer, intent(InOut) :: index
   integer :: i
! END new declarations
   integer, dimension(1:*), intent(Out) :: kpw
   integer, dimension(1:n), intent(InOut) :: mark
! -----------------------------------------------------------------------
!      Laura C. Dutto - email: dutto@cerca.umontreal.ca - December 1993
! -----------------------------------------------------------------------
   do i = 1, n
   if (mark(i) == mccnex) then
      mark(i) = iccnex
      index = index + 1
      kpw(index) = i
   end if
   end do
!
   return
end subroutine newcnx
subroutine blccn1(n, nbloc, nblcmx, nsbloc, lpw, kpw, ia, ja, mark, mxccex, lccnex, mxptbl, iout, ier)
! BEGIN new declarations
   implicit none
   integer :: n
   integer, intent(InOut) :: nbloc
   integer, intent(In) :: nblcmx
   integer, intent(In) :: mxccex
   integer, intent(In) :: mxptbl
   integer :: iout
   integer, intent(Out) :: ier
   integer :: nsb
   integer :: isor
   integer :: newblc
   integer :: ibloc
   integer :: ik0
   integer :: ik1
   integer :: ntb0
   integer :: ins0
   integer :: kpibl
   integer :: nccnex
   integer :: ilccnx
   integer :: info
   integer :: nsfin
   integer :: lcc0
   integer :: icc
! END new declarations
   integer, dimension(1:n), intent(InOut) :: lpw
   integer, dimension(1:mxptbl*nbloc), intent(InOut) :: kpw
   integer, dimension(1:n + 1), intent(In) :: ia
   integer, dimension(1:*), intent(In) :: ja
   integer, dimension(1:(mxccex + 1)*nbloc), intent(InOut) :: lccnex
   integer, dimension(0:nbloc), intent(InOut) :: nsbloc
   integer, dimension(1:n), intent(InOut) :: mark
   logical :: impr
! -----------------------------------------------------------------------
!
!      This routine determines if the matrix given by the structure
!      IA et JA is irreductible. If not, it orders the unknowns such
!      that all the consecutive unknowns in KPW between NSBLOC(i-1)+1
!      and NSBLOC(i) belong to the ith component of the matrix.
!
!      On entry:
!      --------
!      n      = row and column dimension of the matrix
!      nblcmx = The size of NSBLOC is nblcmx + 1 (in fact, it starts at 0).
!      ia     = integer array of length N+1 corresponding to the
!               pointer to the beginning of each row in JA (compressed
!               sparse row storage).
!      ja     = integer array of length NNZERO (= IA(N+1)-IA(1)) corresponding
!               to the column indices of nonzero elements of the matrix,
!               stored rowwise.
!      mxccex = maximum number of connected components allowed by block.
!      mxptbl = maximum number of points (or unknowns) in each connected
!               component (mxptbl .le. n).
!      iout   = impression parameter. If 0 < iout < 100, we print
!               comments and error messages on unit IOUT.
!
!      Input/Output:
!      ------------
!      nbloc  = number of connected components of the matrix. If the
!               matrix is not irreductible, nbloc > 1. We allow
!               nbloc > 1 on entry; in this case we calculate the
!               number of connected components in each previous one.
!      nsbloc = integer array of length NBLOC + 1 containing the pointers
!               to the first node of each component on the new ordering.
!               Normally, on entry you put: NBLOC = 1, NSBLOC(0) = 0,
!               NSBLOC(NBLOC) = N.
!
!      On return:
!      ----------
!      lpw    = integer array of length N corresponding to the
!               permutation vector (the row i goes to lpw(i)).
!      ier    = integer. Error message. Normal return ier = 0.
!
!      Work space:
!      ----------
!      kpw    = integer vector of length MXPTBL*NBLOC necessary for parallel
!               computation.
!      mark   = integer vector of length N
!      lccnex = integer vector of length (MXCCEX+1)*NBLOC necessary for parallel
!               computation.
!
! -----------------------------------------------------------------------
!      Laura C. Dutto - e-mail: dutto@cerca.umontreal.ca
!                       Juillet 1992. Update: March 1994
! -----------------------------------------------------------------------
   character(len=6) :: chsubr
! -----------------------------------------------------------------------
   ier = 0
   nsb = 0
   impr = iout .gt. 0 .and. iout .le. 99
   isor = 0
!
   chsubr = 'CCONEX'
   newblc = 0
! $DOACROSS if(nbloc.gt.1), LOCAL(ibloc, ik0, ik1, ins0, ntb0,
! $&   nccnex, ilccnx, info, kpibl), REDUCTION(ier, newblc)
   do ibloc = 1, nbloc
      ik0 = nsbloc(ibloc - 1)
      ik1 = nsbloc(ibloc)
      ntb0 = ia(ik0 + 1)
      if (ia(ik1 + 1) - ntb0 <= 1) goto 100
      ntb0 = ntb0 - 1
      ins0 = ik1 - ik0
! ........We need more memory place for KPW1 because of parallel computation
      kpibl = (ibloc - 1)*mxptbl
      call numini(ins0, kpw(kpibl + 1))
      nccnex = mxccex
      ilccnx = (mxccex + 1)*(ibloc - 1) + 1
! .......................................................................
!
!         Call to the main routine: CCONEX
!
! .......................................................................
      call cconex(ins0, ik0, nccnex, lccnex(ilccnx), kpw(kpibl + 1), ia(ik0 + 1), ja(ntb0 + 1), mark(ik0 + 1), &
                  isor, info)
      ier = ier + info
      if (info /= 0 .or. nccnex < 1) goto 100
!
! ........We add the new connected components on NEWBLC
      newblc = newblc + nccnex
! ........We define LPW different from the identity only if there are more
! ........than one connected component in this block
      if (nccnex == 1) then
         call numini(ins0, lpw(ik0 + 1))
      else
         call invlpw(ins0, kpw(kpibl + 1), lpw(ik0 + 1))
      end if
      call iplusa(ins0, ik0, 1, lpw(ik0 + 1))
100   continue
   end do
!
   if (ier /= 0) goto 218
   if (newblc == nbloc) goto 120
   if (newblc > nblcmx) goto 230
!
! .....We modify the number of blocks to indicate the number of connected
! .....components in the matrix.
   newblc = 0
   nsfin = 0
! DIR$ NEXT SCALAR
   do ibloc = 1, nbloc
      ilccnx = (mxccex + 1)*(ibloc - 1) + 1
      nccnex = lccnex(ilccnx)
      if (nccnex > 1 .and. impr) write (iout, 420) ibloc, nccnex
      lcc0 = 0
! DIR$ NEXT SCALAR
      do icc = 1, nccnex
         newblc = newblc + 1
         nsb = lccnex(ilccnx + icc)
! ...........Be careful! In LCCNEX we have the cumulated number of vertices
         nsbloc(newblc) = nsfin + nsb
         if (nccnex > 1 .and. impr) write (iout, 425) icc, nsb - lcc0
         lcc0 = nsb
      end do
      nsfin = nsfin + nsb
   end do
   nbloc = newblc
!
120 return
!
218 if (impr) write (iout, 318) chsubr, ier
   goto 120
230 if (impr) write (iout, 330) newblc, nblcmx
   if (ier == 0) ier = -1
   goto 120
!
318 format(' ***BLCCN1*** ERROR IN ', a6, '. IER = ', i8)
330 format(' ***BLCCN1*** THE MEMORY SPACE ALLOWED FOR NSBLOC IS', ' NOT ENOUGH.'/13x, &
          ' NUMBER (NECESSARY) OF CONNECTED', ' COMPONENTS = ', i5/13x, &
          ' MAXIMAL NUMBER OF BLOCKS', 14x, '= ', i5)
420 format(' *** The block ', i3, ' has ', i3, ' strongly connected', &
          ' components. The number of vertices by component is:')
425 format(5x, 'Component No.', i3, ' - Number of vertices = ', i6)
end subroutine blccn1
subroutine ccnicopy(n, ix, iy)
! .......................................................................
!      We copy the vector IX on the vector IY
! .......................................................................
! BEGIN new declarations
   implicit none
   integer, intent(In) :: n
   integer :: i
! END new declarations
   integer, dimension(1:n), intent(In) :: ix
   integer, dimension(1:n), intent(Out) :: iy
! .......................................................................
   if (n <= 0) return
! $DOACROSS if(n .gt. 250), local(i)
   do i = 1, n
      iy(i) = ix(i)
   end do
!
   return
end subroutine ccnicopy
subroutine compos(n, lpw0, lpw1)
! -----------------------------------------------------------------------
!
!      We take account of the original order of unknowns. We put the
!      final result on LPW0.
!
! -----------------------------------------------------------------------
! BEGIN new declarations
   implicit none
   integer, intent(In) :: n
   integer :: i0
! END new declarations
   integer, dimension(1:n), intent(Out) :: lpw0
   integer, dimension(1:n), intent(In) :: lpw1
! -----------------------------------------------------------------------
!      Laura C. Dutto - Mars 1994
! -----------------------------------------------------------------------
! $DOACROSS if(n .gt. 250), local(i0)
   do i0 = 1, n
      lpw0(i0) = lpw1(lpw0(i0))
   end do
!
   return
end subroutine compos
subroutine invlpw(n, lpw, kpw)
! .......................................................................
!
!      KPW is the inverse of LPW
!
! .......................................................................
! BEGIN new declarations
   implicit none
   integer, intent(In) :: n
   integer :: i0
   integer :: i1
! END new declarations
   integer, dimension(1:n), intent(In) :: lpw
   integer, dimension(1:n), intent(Out) :: kpw
! .......................................................................
!      Laura C. Dutto - Novembre 1993
! .......................................................................
! $DOACROSS if(n .gt. 200), local(i0, i1)
   do i0 = 1, n
      i1 = lpw(i0)
      kpw(i1) = i0
   end do
!
   return
end subroutine invlpw
subroutine numini(n, lpw)
! .......................................................................
! BEGIN new declarations
   implicit none
   integer, intent(In) :: n
   integer :: i
! END new declarations
   integer, dimension(1:n), intent(Out) :: lpw
! .......................................................................
!
!      The vector LPW is initialized as the identity.
!
! .......................................................................
!      Laura C. Dutto - Novembre 1993
! .......................................................................
! $DOACROSS if(n .gt. 250), local(i)
   do i = 1, n
      lpw(i) = i
   end do
!
   return
end subroutine numini
subroutine tbzero(m, nmot)
! .......................................................................
!      We initialize to ZERO an integer vector of length NMOT.
! .......................................................................
! BEGIN new declarations
   implicit none
   integer, intent(In) :: nmot
   integer :: i
! END new declarations
   integer, dimension(1:nmot), intent(Out) :: m
! .......................................................................
   if (nmot <= 0) return
! $DOACROSS if(nmot.gt.500), LOCAL(i)
   do i = 1, nmot
      m(i) = 0
   end do
   return
end subroutine tbzero
subroutine iplusa(n, nalpha, nbeta, ia)
! .......................................................................
!
!      We add NALPHA to each element of NBETA * IA:
!
!             ia(i) = nalpha + nbeta * ia(i)
!
! .......................................................................
! BEGIN new declarations
   implicit none
   integer, intent(In) :: n
   integer, intent(In) :: nalpha
   integer, intent(In) :: nbeta
   integer :: nmax
   integer :: i
! END new declarations
   integer, dimension(1:n), intent(InOut) :: ia
! .......................................................................
!      Laura C. Dutto - February 1994
! .......................................................................
   if (n <= 0) return
!
   nmax = 500
   if (nalpha == 0) then
      if (nbeta == 1) return
      if (nbeta == -1) then
! $DOACROSS if(n .gt. nmax), local (i)
         do i = 1, n
            ia(i) = -ia(i)
         end do
      else
! $DOACROSS if(n .gt. nmax/2), local (i)
         do i = 1, n
            ia(i) = nbeta*ia(i)
         end do
      end if
      return
   end if
   if (nbeta == 0) then
! $DOACROSS if(n .gt. nmax), local (i)
      do i = 1, n
         ia(i) = nalpha
      end do
      return
   end if
   if (nbeta == -1) then
! $DOACROSS if(n .gt. nmax/2), local (i)
      do i = 1, n
         ia(i) = nalpha - ia(i)
      end do
   else if (nbeta == 1) then
! $DOACROSS if(n .gt. nmax/2), local (i)
      do i = 1, n
         ia(i) = nalpha + ia(i)
      end do
   else
! $DOACROSS if(n .gt. nmax/3), local (i)
      do i = 1, n
         ia(i) = nalpha + nbeta*ia(i)
      end do
   end if
!
   return
end subroutine iplusa
