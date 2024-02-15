! ----------------------------------------------------------------------c
!                           S P A R S K I T                             c
! ----------------------------------------------------------------------c
!                   INFORMATION ROUTINES. INFO MODULE                   c
! ----------------------------------------------------------------------c
!  bandwidth :  Computes  ml     = lower_bandwidth(A)                   c
!                         mu     = upper_bandwidth(A)                   c
!                         iband  = max_bandwidth(A)                     c
!                         bndav  = average_bandwidth(A)                 c
!  nonz      :  Computes  nzmaxc = max_column_length(A)                 c
!                         nzminc = min_column_length(A)                 c
!                         nzmaxr = max_row_length(A)                    c
!                         nzminr = min_row_length(A)                    c
!                         nzcol  = zero_column_number(A)                c
!                         nzrow  = zero_row_number(A)                   c
!  diag_domi :  Computes  ddomc  = diag_domi_column_percentage(A)       c
!                         ddomr  = diag_domi_row_percentage(A)          c
!  frobnorm  :  Computes  Fnorm  = Frobenius_norm(A)                    c
!  ansym     :  Computes  fas    = sym_part_Frobenius_norm(A)           c
!                         fan    = nonsym_part_Frobenius_norm(A)        c
!                         imatch = matching_elements_number(A)          c
!                         av     = relative_sym_match(A)                c
!  distaij   :  Computes  dist   = average_dist_of_a(i,j)(A)            c
!                         std    = standard_deviation(A)                c
!  skyline   :  Computes  nsky   = nonzero_number_in_skyline(A)         c
!  distdiag  :  Computes  dist   = element_number_in_eachdiag(A)        c
!  bandpart  :  Computes  band   = bandwidth_width(A)                   c
!  n_imp_diag:  Computes  ndiag  = important_diag_number(A)             c
!  nonz_lud  :  Computes  nlower = nonzero_number_of_lower_part(A)      c
!                         nupper = nonzero_number_of_upper_part(A)      c
!                         ndiag  = nonzero_number_of_maindiag(A)        c
!  avnz_col  :  Computes  av     = average_nonzero_number_in_column(A)  c
!                         st     = standard_deviation(A)                c
!  vbrinfo   :  Print info on matrix in variable block row format       c
! ----------------------------------------------------------------------c
subroutine bandwidth(n,ja,ia,ml,mu,iband,bndav)

     integer, intent(In) :: n
     integer, intent(Out) :: ml,mu,iband
     integer, dimension(*), intent(In) :: ja
     integer, dimension(n + 1), intent(In) :: ia
     real(kind=8), intent(Out) :: bndav
! -----------------------------------------------------------------------
!  this routine computes the lower, upper, maximum, and average 
!  bandwidths.     revised -- July 12, 2001  -- bug fix -- YS. 
! -----------------------------------------------------------------------
!  On Entry:
! ----------
!  n     = integer. column dimension of matrix
!  a     = real array containing the nonzero elements of the matrix
!          the elements are stored by columns in order
!          (i.e. column i comes before column i+1, but the elements
!          within each column can be disordered).
!  ja    = integer array containing the row indices of elements in a
!  ia    = integer array containing of length n+1 containing the
!          pointers to the beginning of the columns in arrays a and ja.
!          It is assumed that ia(*) = 1 and ia(n+1) = nzz+1.
! 
!  on return
! ----------
!  ml    = lower bandwidth as defined by
!         ml = max(i-j | all  a(i,j).ne. 0)
!  mu    = upper bandwidth as defined by
!         mu = max ( j-i | all  a(i,j).ne. 0 )
!  iband =  maximum bandwidth as defined by
!          iband = Max (  Max [ j | a(i,j) .ne. 0 ] - 
!                         Min [ j | a(i,j) .ne. 0 ] )
!  bndav = Average Bandwidth          
! -----------------------------------------------------------------------
!      locals
     integer :: max, j0, j1,jminc,jmaxc,i
! -----------------------------------------------------------------------
     ml = -n
     mu = -n
     bndav = 0.0d0
     iband = 0 
     do i=1,n
          j0 = ia(i)
          j1 = ia(i+1) - 1
          jminc = ja(j0)
          jmaxc = ja(j1)
          ml = max(ml,i-jminc)
          mu = max(mu,jmaxc-i)
          iband = max(iband,jmaxc-jminc+1)
          bndav = bndav+real( jmaxc-jminc+1)
     end do
     bndav = bndav/real(n)
     return
! -----end-of-bandwidth--------------------------------------------------
! -----------------------------------------------------------------------
end subroutine bandwidth

subroutine nonz(n,sym,ja,ia,iao,nzmaxc,nzminc,nzmaxr,nzminr,nzcol,nzrow)

     integer, intent(In) :: n
     integer, intent(Out) :: nzmaxc,nzminc,nzmaxr,nzminr,nzcol,nzrow
     integer, dimension(*), intent(In) :: ja
     integer, dimension(n + 1), intent(In) :: ia,iao
     logical, intent(In) :: sym
! ----------------------------------------------------------------------
!      this routine computes maximum numbers of nonzero elements 
!      per column/row, minimum numbers of nonzero elements per column/row, 
!      and  numbers of zero columns/rows.
! ----------------------------------------------------------------------
!      On Entry:
! ----------
!      n     = integer column dimension of matrix
!      ja    = integer array containing the row indices of elements in a
!      ia    = integer array containing of length n+1 containing the
!      pointers to the beginning of the columns in arrays a and ja.
!      It is assumed that ia(*) = 1 and ia(n+1) = nzz+1.
!      iao   = similar array for the transpose of the matrix. 
!      sym   = logical variable indicating whether or not the matrix is 
!      stored in symmetric mode.
!      on return
! ----------
!      nzmaxc = max length of columns
!      nzminc = min length of columns 
!      nzmaxr = max length of rows
!  nzminr = min length of rows
!      nzcol  = number of zero columns
!      nzrow = number of zero rows
! -----------------------------------------------------------------------
     integer :: i,j0,j0r,j1r,indiag,k,j1,lenc,lenr
! 

     nzmaxc = 0
     nzminc = n
     nzmaxr = 0
     nzminr = n
     nzcol = 0
     nzrow = 0
! -----------------------------------------------------------------------
     do i = 1, n
          j0 = ia(i)
          j1 = ia(i+1) 
          j0r = iao(i)
          j1r = iao(i+1)
          indiag = 0
          do k=j0, j1-1 
               if (ja(k) == i) indiag = 1
          end do
!          
          lenc = j1-j0
          lenr = j1r-j0r
!          
          if (sym) lenc = lenc + lenr - indiag
          if (lenc <= 0) nzcol = nzcol +1
          nzmaxc = max0(nzmaxc,lenc)
          nzminc = min0(nzminc,lenc)
          if (lenr <= 0) nzrow = nzrow+1
          nzmaxr = max0(nzmaxr,lenr)
          nzminr = min0(nzminr,lenr)
     end do
     return
end subroutine nonz

subroutine diag_domi(n,sym,valued,a,ja,ia,ao,jao,iao,ddomc,ddomr)
     
     real(kind=8), dimension(*), intent(In) :: a, ao
     real(kind=8), intent(Out) :: ddomc, ddomr
     integer, intent(In) :: n
     integer, dimension(*), intent(In) :: ja, jao
     integer, dimension(n + 1), intent(In) :: ia, iao
     logical, intent(In) :: sym, valued
! -----------------------------------------------------------------
!      this routine computes the percentage of weakly diagonally 
!      dominant rows/columns
! -----------------------------------------------------------------
!      on entry:
!      ---------
!      n     = integer column dimension of matrix
!  a     = real array containing the nonzero elements of the matrix
!      the elements are stored by columns in order
!      (i.e. column i comes before column i+1, but the elements
!      within each column can be disordered).
!      ja    = integer array containing the row indices of elements in a
!  ia    = integer array containing of length n+1 containing the
!      pointers to the beginning of the columns in arrays a and ja.
!      It is assumed that ia(*) = 1 and ia(n+1) = nzz+1.
!      ao    = real array containing the nonzero elements of the matrix
!      the elements are stored by rows in order
!      (i.e. row i comes before row i+1, but the elements
!      within each row can be disordered).
!  ao,jao, iao, 
!      structure for transpose of a 
!      sym   = logical variable indicating whether or not the matrix is
!      symmetric.
!      valued= logical equal to .true. if values are provided and .false.
!          if only the pattern of the matrix is provided. (in that
!      case a(*) and ao(*) are dummy arrays.
!      
!      ON RETURN
! ----------
!      ddomc = percentage of weakly diagonally dominant columns
!      ddomr = percentage of weakly diagonally dominant rows
! -------------------------------------------------------------------
!      locals
     integer :: i, j0, j1, k, j
     real(kind=8) :: aii, dsumr, dsumc

!      number of diagonally dominant columns and number of diagonally dominant rows
!      real arithmetic used to avoid problems.. YS. 03/27/01 
    
     ddomc = 0.0
     ddomr = 0.0
     do i = 1, n
          j0 = ia(i)
          j1 = ia(i+1) - 1
          if (valued) then
               aii = 0.0d0
               dsumc = 0.0d0
               do k=j0,j1
                    j = ja(k) 
                    if (j == i) then
                         aii = abs(a(k))
                    else
                         dsumc = dsumc + abs(a(k))
                    endif
               end do
               dsumr = 0.0d0
               if (.not. sym) then
                    do k=iao(i), iao(i+1)-1
                         if (jao(k) /= i) dsumr = dsumr+abs(ao(k))
                    end do
               else
                    dsumr = dsumc
               endif
               if (dsumc <= aii) ddomc = ddomc + 1.0
               if (dsumr <= aii) ddomr = ddomr + 1.0
          endif
     end do
     ddomr = ddomr / real(n)
     ddomc = ddomc / real(n)
     return
! -----------------------------------------------------------------------
! --------end-of-diag_moni-----------------------------------------------
end subroutine diag_domi

subroutine frobnorm(n,sym,a,ja,ia,fnorm)

     integer, intent(In) :: n
     real(kind=8), dimension(*), intent(In) :: a
     real(kind=8), intent(Out) :: fnorm
     integer, dimension(*), intent(In) :: ja
     integer, dimension(n + 1), intent(In) :: ia
     logical, intent(In) :: sym
! --------------------------------------------------------------------------
!      this routine computes the Frobenius norm of A.
! --------------------------------------------------------------------------
!      on entry:
! -----------
!  n      = integer colum dimension of matrix
!  a      = real array containing the nonzero elements of the matrix
!           the elements are stored by columns in order
!           (i.e. column i comes before column i+1, but the elements
!           within each column can be disordered).
!  ja     = integer array containing the row indices of elements in a.
!  ia     = integer array containing of length n+1 containing the
!           pointers to the beginning of the columns in arrays a and 
!           ja. It is assumed that ia(*)= 1 and ia(n+1)  = nnz +1.
!  sym    = logical variable indicating whether or not the matrix is
!           symmetric.
! 
!  on return
! -----------
!  Fnorm  = Frobenius norm of A.
! --------------------------------------------------------------------------
     real(kind=8) :: fdiag
     integer :: i, k

     fdiag = 0.0
     fnorm = 0.0
     do i =1,n
          do k = ia(i), ia(i+1)-1
               if (ja(k) == i) then
                    fdiag = fdiag + a(k)**2
               else
                    fnorm = fnorm + a(k)**2
               endif
          enddo 
     enddo 
     if (sym) then
          fnorm = 2*fnorm +fdiag
     else
          fnorm = fnorm + fdiag
     endif
     fnorm = sqrt(fnorm)
     return
end subroutine frobnorm

subroutine ansym(n,sym,a,ja,ia,ao,jao,iao,imatch,av,fas,fan)
! -----------------------------------------------------------------------
! ---------------------------------------------------------------------
!      this routine computes the Frobenius norm of the symmetric and
!      non-symmetric parts of A, computes number of matching elements
!      in symmetry and relative symmetry match. 
! ---------------------------------------------------------------------
!  on entry:
! ----------
!  n   = integer column dimension of matrix
!  a   = real array containing the nonzero elements of the matrix
!        the elements are stored by columns in order
!        (i.e. column i comes before column i+1, but the elements
!        within each column can be disordered).
!  ja  = integer array containing the row indices of elements in a
!  ia  = integer array containing of length n+1 containing the
!        pointers to the beginning of the columns in arrays a and ja.
!        It is assumed that ia(*) = 1 and ia(n+1) = nzz+1.
!  sym = logical variable indicating whether or not the matrix is
!        symmetric.
!  on return
! ----------
!  fas   = Frobenius norm of symmetric part
!  fan   = Frobenius norm of non-symmetric part
!  imatch = number of matching elements in symmetry
!  av     = relative symmetry match (symmetry = 1)
!  ao,jao,iao = transpose of A just as a, ja, ia contains 
!               information of A.
     integer :: nnz, i, k1, k2, k1max, k2max, j1, j2
     real(kind=8), dimension(*), intent(In) :: a, ao
     real(kind=8), intent(Out) :: fas, fan, av
     real(kind=8) :: fnorm, st
     integer, intent(In) :: n
     integer, dimension(*), intent(In) :: ja, jao
     integer, dimension(n + 1), intent(In) :: ia, iao
     integer, intent(Out) :: imatch
     logical, intent(In) :: sym
! -----------------------------------------------------------------------
     st = 0.0
     nnz = ia(n+1)-ia(1)
     call csrcsc(n,1,1,a,ja,ia,ao,jao,iao)
     select case (sym)
          case(.true.)
               call frobnorm(n, sym,a,ja,ia,fnorm)
               if (sym) then
                    imatch = nnz
                    fas = fnorm
                    fan = 0.0
               else
                    if (imatch == nnz) then
                         st = 0.0
                    else 
                         st = 0.5 * (fnorm**2 - st)
                         if (st < 0.0) st = 0.0
                    end if
                    fas = sqrt(fas + st)
                    fan = sqrt(fan + st)
               end if
               av = real(imatch)/real(nnz)
          case(.false.)
          
               st = 0.0
               fas = 0.0
               fan = 0.0
               imatch = 0
               do i = 1,n
                    k1 = ia(i)
                    k2 = iao(i)
                    k1max = ia(i+1) - 1
                    k2max = iao(i+1) - 1
                    do while(k1 > k1max .or. k2 > k2max)
                         j1 = ja(k1)
                         j2 = jao(k2)
                         if (j1 /= j2) then
                              k1 = k1 + 1
                              k2 = k2 + 1
                              if (j1 < j2) k2 = k2 - 1
                              if (j1 > j2) k1 = k1 - 1
                         else
                              fas = fas + (a(k1) + ao(k2))**2
                              fan = fan + (a(k1) - ao(k2))**2
                              st = st + a(k1)**2
                              imatch = imatch + 1
                              k1 = k1 + 1
                              k2 = k2 + 1
                              if (j1 < j2) k2 = k2 - 1
                              if (j1 > j2) k1 = k1 - 1
                         end if
                    end do 
               end do
     end select
     av = real(imatch)/real(nnz)
     return
end subroutine ansym

subroutine distaij(n,nnz,ja,ia,dist,std)
 
     integer, intent(In) :: n, nnz
     integer :: i, j0, j1,k ,j
     real(kind=8), intent(Out) :: dist, std
     integer, dimension(*), intent(In) :: ja
     integer, dimension(n + 1), intent(In) :: ia
! -----------------------------------------------------------------------
!      this routine computes the average distance of a(i,j) from diag and
!      standard deviation  for this average.
! -----------------------------------------------------------------------
!  On entry :
! -----------
!  n     = integer. column dimension of matrix
!  nnz   = number of nonzero elements of matrix
!  ja    = integer array containing the row indices of elements in a
!  ia    = integer array containing of length n+1 containing the
!          pointers to the beginning of the columns in arrays a and ja.
!          It is assumed that ia(*) = 1 and ia(n+1) = nzz+1.
!  sym   = logical variable indicating whether or not the matrix is
!          symmetric.
!  on return
! ----------
!  dist  = average distance of a(i,j) from diag.
!  std   = standard deviation for above average.
! -----------------------------------------------------------------------
! 
!  distance of an element from diagonal.
! 
     dist = 0.0
     std = 0.0
     do i=1,n
          j0 = ia(i)
          j1 = ia(i+1) - 1
          do k=j0, j1
               j=ja(k)
               dist = dist + real(iabs(j-i) )
          end do
     end do
     dist = dist/real(nnz)
     do i = 1, n 
          do k=ia(i), ia(i+1) - 1
               std=std+(dist-real(iabs(ja(k)-i)))**2
          end do
     end do
     std = sqrt(std/ real(nnz))
     return
end subroutine distaij

subroutine skyline(n,sym,ja,ia,jao,iao,nsky)
     integer :: i, j0, j0r, jminc,jminr, nskyl, nskyu
     integer, intent(In) :: n
     integer, dimension(*), intent(In) :: ja, jao
     integer, dimension(n + 1), intent(In) :: ia, iao
     integer, intent(Out) :: nsky
     logical, intent(In) :: sym
! -------------------------------------------------------------------
!  this routine computes the number of nonzeros in the skyline storage.
! -------------------------------------------------------------------
! 
!  On entry :
! -----------
!  n     = integer. column dimension of matrix
!  ja    = integer array containing the row indices of elements in a
!  ia    = integer array containing of length n+1 containing the
!          pointers to the beginning of the columns in arrays a and ja.
!          It is assumed that ia(*) = 1 and ia(n+1) = nzz+1.
!  iao   = integer array containing of length n+1 containing the
!          pointers to the beginning of the rows in arrays ao and jao.
!          It is assumed that iao(*) = 1 and iao(n+1) = nzz+1.
!  jao   = integer array containing the column indices of elements in ao.
!  sym   = logical variable indicating whether or not the matrix is
!          symmetric.
!  on return
! ----------
!  nsky  = number of nonzeros in skyline storage
! -------------------------------------------------------------------
! 
!  nskyu = skyline storage for upper part

!  nskyl = skyline storage for lower part

     data nkyu, nskyl /0,0/
     do i=1,n
          j0 = ia(i)
          j0r = iao(i)
          jminc = ja(j0)
          jminr = jao(j0r)
          if (sym) jminc = jminr
          nskyl = nskyl + i-jminr + 1
          nskyu = nskyu + i-jminc + 1
     end do
     nsky = nskyl+nskyu-n
     if (sym) nsky = nskyl
     return
end subroutine skyline

subroutine distdiag(nrow,ncol,ja,ia,dist)

     integer :: n2, i,k1,k2,k,j
     integer, intent(In) :: nrow, ncol
     integer, dimension(*), intent(In) :: ja
     integer, dimension(nrow + 1), intent(In) :: ia
     integer, dimension(*), intent(Out) :: dist
! ----------------------------------------------------------------------
!  this routine computes the numbers of elements in each diagonal. 
! ----------------------------------------------------------------------
!  On entry :
! -----------
!  nrow  = integer. row dimension of matrix
!  ncol  = integer. column dimension of matrix
!  ja    = integer array containing the row indices of elements in a
!  ia    = integer array containing of length n+1 containing the
!          pointers to the beginning of the columns in arrays a and ja.
!          It is assumed that ia(*) = 1 and ia(n+1) = nzz+1.
!  on return
! ----------
!  dist  = integer array containing the numbers of elements in each of 
!          the nrow+ncol-1 diagonals of A. dist(k) contains the 
!          number of elements in diagonal '-nrow+k'.  k ranges from 
!          1 to (nrow+ncol-1).
! ----------------------------------------------------------------------
! 
     n2 = nrow + ncol - 1
     do i=1, n2
          dist(i) = 0
     end do
     
     do i=1, nrow
          k1 = ia(i)
          k2 = ia(i+1) -1
          do k=k1, k2
               j = ja(k)
               dist(nrow+j-i) = dist(nrow+j-i) +1
          end do
     end do
     return
end subroutine distdiag

subroutine bandpart(n,ia,dist,nper,band)

     integer :: nnz, iacc, j, n
     integer, dimension(n + 1), intent(In) :: ia
     integer, dimension(*), intent(In) :: dist
     integer, intent(In) :: nper
     integer, intent(Out) :: band
! -------------------------------------------------------------------------
!  this routine computes the bandwidth of the banded matrix, which contains
!  'nper' percent of the original matrix.
! -------------------------------------------------------------------------
!  On entry :
! -----------
!  n     = integer. column dimension of matrix
!  a, ja, ia = matrix in CSR format. only ia used/
!  ia    = integer array containing of length n+1 containing the
!          pointers to the beginning of the columns in arrays a and ja.
!          It is assumed that ia(*) = 1 and ia(n+1) = nzz+1.
!  dist  = integer array containing the numbers of elements in the 
!          matrix with different distance of row indices and column 
!          indices.
!  nper  = percentage of matrix  within the bandwidth
!  on return
! ----------
!  band  = the width of the bandwidth
! ----------------------------------------------------------------------
     nnz = ia(n+1)-ia(1)
     iacc = dist(n)
     band = 0
     j = 1
     iacc = iacc + dist(n+j) + dist(n-j)
     do while (iacc*100 <= nnz*nper)
          j = j + 1
          iacc = iacc + dist(n + j) + dist (n - j)
          band = band + 1
     end do
     return 
end subroutine bandpart

subroutine n_imp_diag(n,nnz,dist,ipar1,ndiag,ioff,dcount)
     integer :: n2, itot, ii, jmax, i, k, j
     real(kind=8), dimension(*), intent(Out) :: dcount
     integer, intent(In) :: n, nnz, ipar1
     integer, dimension(*), intent(InOut) :: dist
     integer, intent(Out) :: ndiag
     integer, dimension(*), intent(Out) :: ioff
! -----------------------------------------------------------------------
!      this routine computes the most important diagonals.
! -----------------------------------------------------------------------
! 
!  On entry :
! -----------
!  n     = integer. column dimension of matrix
!  nnz   = number of nonzero elements of matrix
!  dist  = integer array containing the numbers of elements in the 
!          matrix with different distance of row indices and column 
!          indices. ipar1 = percentage of nonzero  elements of A that 
!          a diagonal should have in order to be an important diagonal 
!  on return
! ----------
!  ndiag = number of the most important diagonals
!  ioff  = the offsets with respect to the main diagonal
!  dcount= the accumulated percentages
! -----------------------------------------------------------------------
     n2 = 2*n - 1
     ndiag = min0(n2,10)
     itot = 0
     ii = 0
!      sort diagonals by decreasing order of weights.
     do while (ii < ndiag)
          jmax = 0
          i = 1
          do k = 1, n2
               j = dist(k)
               if (j >= jmax) then
                    i = k
                    jmax = j
               end if
          end do
!      permute ----
!      save offsets and accumulated count if diagonal is acceptable
!      (if it has at least ipar1*nnz/100 nonzero elements)
!      quite if no more acceptable diagonals --
!      
          if  (100*jmax < ipar1*nnz) exit
          ii = ii + 1
          ioff(ii) = i - n
          dist(i) = -jmax
          itot = itot + jmax
          dcount(ii) = real(100*itot)/real(nnz) 
     end do
     ndiag = ii
     return
! -----------------------------------------------------------------------
end subroutine n_imp_diag

subroutine nonz_lud(n,ja,ia,nlower,nupper,ndiag)
     integer :: i,j,k
     integer, intent(In) :: n
     integer, dimension(*), intent(In) :: ja
     integer, dimension(n + 1), intent(In) :: ia
     integer, intent(Out) :: nlower, nupper, ndiag
! -----------------------------------------------------------------------
!  this routine computes the number of nonzero elements in strict lower
!  part, strict upper part, and main diagonal.
! -----------------------------------------------------------------------
! 
!  On entry :
! -----------
!  n     = integer. column dimension of matrix
!  ja    = integer array containing the row indices of elements in a
!  ia    = integer array containing of length n+1 containing the
!          pointers to the beginning of the columns in arrays a and ja.
!          It is assumed that ia(*) = 1 and ia(n+1) = nzz+1.
!  on return
! ----------
!  nlower= number of nonzero elements in strict lower part
!  nupper= number of nonzero elements in strict upper part
!  ndiag = number of nonzero elements in main diagonal
! ------------------------------------------------------------------- 
! 
!  number of nonzero elements in upper part
! 
     nupper = 0
     ndiag = 0
     do i=1,n
!      indiag = nonzero diagonal element indicator
          do k=ia(i), ia(i+1)-1
               j=ja(k)
               if (j < i) nupper = nupper+1
               if (j == i) ndiag = ndiag + 1 
          end do
     end do
     nlower = ia(n+1)-1-nupper-ndiag
     return
end subroutine nonz_lud

subroutine avnz_col(n,ia,av,st)
      integer :: i, j0, j1, lenc
     real(kind=8), intent(Out) :: av, st
     integer, intent(In) :: n
     integer, dimension(n + 1), intent(In) :: ia
! ---------------------------------------------------------------------
!      this routine computes average number of nonzero elements/column and
!      standard deviation for this average
! ---------------------------------------------------------------------
! 
!  On entry :
! -----------
!  n     = integer. column dimension of matrix
!  ja    = integer array containing the row indices of elements in a
!  ia    = integer array containing of length n+1 containing the
!          pointers to the beginning of the columns in arrays a and ja.
!          It is assumed that ia(*) = 1 and ia(n+1) = nzz+1.
!  ndiag = number of the most important diagonals
!  On return
! ----------
!  av    = average number of nonzero elements/column
!  st    = standard deviation for this average
!  Notes
! ---------
!  standard deviation will not be correct for symmetric storage. 
! ----------------------------------------------------------------------
!      standard deviatioan for the average
     st = 0.0d0
!      average and standard deviation
!      
     av = real(ia(n+1)-1)/real(n)
!      
!      will be corrected later.
!      
     do i=1,n
          j0 = ia(i)
          j1 = ia(i+1) - 1
!      indiag = nonzero diagonal element indicator
          lenc = j1+1-j0
          st = st + (real(lenc) - av)**2
     end do
!      
     st = sqrt( st / real(n) )
     return
end subroutine avnz_col

subroutine vbrinfo(nr,nc,kvstr,kvstc,ia,ka,iwk)

     integer, intent(In) :: nr, nc
     integer, dimension(nr + 1), intent(In) :: kvstr
     integer, dimension(nc + 1), intent(In) :: kvstc
     integer, dimension(nr + 1), intent(In) :: ia
     integer, dimension(*), intent(In) :: ka
     integer, dimension(2*nr + nc + 2), intent(Out) :: iwk
! -----------------------------------------------------------------------
!      Print info on matrix in variable block row format.
! -----------------------------------------------------------------------
!      On entry:
! --------------
!      nr,nc   = matrix block row and block column dimension
!      kvstr   = first row number for each block row
!      kvstc   = first column number for each block column
!      ia,ja,ka,a = input matrix in VBR format
!      iout    = unit number for printed output
! 
!      On return:
! ---------------
!      Printed output to unit number specified in iout.  If a non-square
!      matrix is provided, the analysis will be performed on the block
!      rows, otherwise a row/column conformal partitioning will be used.
! 
!      Work space:
! ----------------
!      iwk(1:nb+1) = conformal block partitioning
!         (nb is unknown at start but is no more than nr+nc)
!      iwk(nb+2:nb+2+nr) = frequency of each blocksize
!      The workspace is not assumed to be initialized to zero, nor is it
!      left that way.
! 
! -----------------------------------------------------------------------
! -----local variables
     integer :: n, nb, nnz, nnzb,i,j,neq,max,num
     character(len=101) :: tmpst
     integer, dimension(10) :: bsiz, freq
! -----------------------------------------------------------------------
     n = kvstr(nr+1)-1
     nnz = ka(ia(nr+1)) - ka(1)
     nnzb = ia(nr+1) - ia(1)
     print '(6x,A,65(1h-),A)','*','*'
     
     print '(6x,A,i10,A,6x,i10,A,6x,f10.4,A)', '* Number of rows = ',n,' * Number of nonzero elements = ',nnz, &
                         ' * Average number of nonzero elements/row = ',real(nnz)/real(n), ' *'
     
     print '(6x,A,i10,6x,A,i10,6x,A,f10.4,A)','* Number of block rows = ',nr,' * Number of nonzero blocks = ',nnzb, ' Average number of &
                         nonzero blocks/ block rows = ',real(nnzb)/real(nr), ' * ' 
     
! -----if non-square matrix, do analysis on block rows,
!      else do analysis on conformal partitioning
     if (kvstr(nr+1) /= kvstc(nc+1)) then
          print '(6x,A,6x)', '* Non-square matrix * Performing analysis on block rows *'
          
          do i = 1, nr+1
               iwk(i) = kvstr(i)
          enddo
          nb = nr
     else
          call kvstmerge(nr, kvstr, nc, kvstc, nb, iwk)
          if ((nr  /=  nc) .or. (nc  /=  nb)) print '(6x,A,i10,A)', '* Non row-column conformal partitioning supplied; &
                                                            using conformal partitioning. Number of bl rows = ', nb,' *'
     endif
! -----accumulate frequencies of each blocksize

     max = 1
     iwk(1+nb+2) = 0
     
     do i = 1, nb
          neq = iwk(i+1) - iwk(i)
          if (neq > max) then
               do j = max+1, neq
                    iwk(j+nb+2) = 0
               enddo
               max = neq
          endif
          iwk(neq+nb+2) = iwk(neq+nb+2) + 1
     enddo
! -----store largest 10 of these blocksizes
     num = 0
     do i = max, 1, -1
          if ((iwk(i+nb+2) /= 0) .and. (num < 10)) then
               num = num + 1
               bsiz(num) = i
               freq(num) = iwk(i+nb+2)
          endif
     enddo
! -----print information about blocksizes
     print '(6x,A,i10,A)', '* Number of different blocksizes = ',num, ' *'
     write (tmpst,'(10i6)') (bsiz(j),j=1,num)
     write (iout,110) num,tmpst
     write (tmpst,'(10i6)') (freq(j),j=1,num)
     write (iout,111) tmpst
     print '(6x,A,65(1h-),A)','*','*'
! -----------------------------------------------------------------------
 110   format(6x,' *  The ', i2, ' largest dimension nodes',     ' have dimension    : ',10x,'  *', &
      /, 6x,' *',a61,3x,' *')
 111   format(6x,' *  The frequency of nodes these ',     'dimensions are      : ',10x,'  *',/, 6x, &
      ' *',a61,3x,' *')
! ---------------------------------
     return
end subroutine vbrinfo