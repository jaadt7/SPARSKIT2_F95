!*==bandwidth.f90 processed by SPAG 8.04RA 15:47 21 Feb 2024
!!SPAG Open source Personal, Educational or Academic User Clemson University  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
!----------------------------------------------------------------------c
!                          S P A R S K I T                             c
!----------------------------------------------------------------------c
!                  INFORMATION ROUTINES. INFO MODULE                   c
!----------------------------------------------------------------------c
! bandwidth :  Computes  ml     = lower_bandwidth(A)                   c
!                        mu     = upper_bandwidth(A)                   c
!                        iband  = max_bandwidth(A)                     c
!                        bndav  = average_bandwidth(A)                 c
! nonz      :  Computes  nzmaxc = max_column_length(A)                 c
!                        nzminc = min_column_length(A)                 c
!                        nzmaxr = max_row_length(A)                    c
!                        nzminr = min_row_length(A)                    c
!                        nzcol  = zero_column_number(A)                c
!                        nzrow  = zero_row_number(A)                   c
! diag_domi :  Computes  ddomc  = diag_domi_column_percentage(A)       c
!                        ddomr  = diag_domi_row_percentage(A)          c
! frobnorm  :  Computes  Fnorm  = Frobenius_norm(A)                    c
! ansym     :  Computes  fas    = sym_part_Frobenius_norm(A)           c
!                        fan    = nonsym_part_Frobenius_norm(A)        c
!                        imatch = matching_elements_number(A)          c
!                        av     = relative_sym_match(A)                c
! distaij   :  Computes  dist   = average_dist_of_a(i,j)(A)            c
!                        std    = standard_deviation(A)                c
! skyline   :  Computes  nsky   = nonzero_number_in_skyline(A)         c
! distdiag  :  Computes  dist   = element_number_in_eachdiag(A)        c
! bandpart  :  Computes  band   = bandwidth_width(A)                   c
! n_imp_diag:  Computes  ndiag  = important_diag_number(A)             c
! nonz_lud  :  Computes  nlower = nonzero_number_of_lower_part(A)      c
!                        nupper = nonzero_number_of_upper_part(A)      c
!                        ndiag  = nonzero_number_of_maindiag(A)        c
! avnz_col  :  Computes  av     = average_nonzero_number_in_column(A)  c
!                        st     = standard_deviation(A)                c
! vbrinfo   :  Print info on matrix in variable block row format       c
!----------------------------------------------------------------------c
SUBROUTINE bandwidth(N,Ja,Ia,Ml,Mu,Iband,Bndav)
   USE ISO_FORTRAN_ENV                 
   IMPLICIT NONE
!
! Dummy argument declarations rewritten by SPAG
!
   INTEGER , INTENT(IN) :: N
   INTEGER , INTENT(IN) , DIMENSION(*) :: Ja
   INTEGER , INTENT(IN) , DIMENSION(N+1) :: Ia
   INTEGER , INTENT(INOUT) :: Ml
   INTEGER , INTENT(INOUT) :: Mu
   INTEGER , INTENT(INOUT) :: Iband
   REAL(REAL64) , INTENT(INOUT) :: Bndav
!
! Local variable declarations rewritten by SPAG
!
   INTEGER :: i , j0 , j1 , jmaxc , jminc
!
! End of declarations rewritten by SPAG
!
!-----------------------------------------------------------------------
! this routine computes the lower, upper, maximum, and average
! bandwidths.     revised -- July 12, 2001  -- bug fix -- YS.
!-----------------------------------------------------------------------
! On Entry:
!----------
! n     = integer. column dimension of matrix
! a     = real array containing the nonzero elements of the matrix
!         the elements are stored by columns in order
!         (i.e. column i comes before column i+1, but the elements
!         within each column can be disordered).
! ja    = integer array containing the row indices of elements in a
! ia    = integer array containing of length n+1 containing the
!         pointers to the beginning of the columns in arrays a and ja.
!         It is assumed that ia(*) = 1 and ia(n+1) = nzz+1.
!
! on return
!----------
! ml    = lower bandwidth as defined by
!        ml = max(i-j | all  a(i,j).ne. 0)
! mu    = upper bandwidth as defined by
!        mu = max ( j-i | all  a(i,j).ne. 0 )
! iband =  maximum bandwidth as defined by
!         iband = Max (  Max [ j | a(i,j) .ne. 0 ] -
!                        Min [ j | a(i,j) .ne. 0 ] )
! bndav = Average Bandwidth
!-----------------------------------------------------------------------
!     locals
!-----------------------------------------------------------------------
   Ml = -N
   Mu = -N
   Bndav = 0.0D0
   Iband = 0
   DO i = 1 , N
      j0 = Ia(i)
      j1 = Ia(i+1) - 1
      jminc = Ja(j0)
      jmaxc = Ja(j1)
      Ml = max(Ml,i-jminc)
      Mu = max(Mu,jmaxc-i)
      Iband = max(Iband,jmaxc-jminc+1)
      Bndav = Bndav + real(jmaxc-jminc+1)
   ENDDO
   Bndav = Bndav/real(N)
!-----end-of-bandwidth--------------------------------------------------
!-----------------------------------------------------------------------
END SUBROUTINE bandwidth
!*==nonz.f90 processed by SPAG 8.04RA 15:47 21 Feb 2024
!!SPAG Open source Personal, Educational or Academic User Clemson University  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
!-----------------------------------------------------------------------
SUBROUTINE nonz(N,Sym,Ja,Ia,Iao,Nzmaxc,Nzminc,Nzmaxr,Nzminr,Nzcol,Nzrow)
   IMPLICIT NONE
!
! Dummy argument declarations rewritten by SPAG
!
   INTEGER , INTENT(IN) :: N
   LOGICAL , INTENT(IN) :: Sym
   INTEGER , INTENT(IN) , DIMENSION(*) :: Ja
   INTEGER , INTENT(IN) , DIMENSION(N+1) :: Ia
   INTEGER , INTENT(IN) , DIMENSION(N+1) :: Iao
   INTEGER , INTENT(INOUT) :: Nzmaxc
   INTEGER , INTENT(INOUT) :: Nzminc
   INTEGER , INTENT(INOUT) :: Nzmaxr
   INTEGER , INTENT(INOUT) :: Nzminr
   INTEGER , INTENT(INOUT) :: Nzcol
   INTEGER , INTENT(INOUT) :: Nzrow
!
! Local variable declarations rewritten by SPAG
!
   INTEGER :: i , indiag , j0 , j0r , j1 , j1r , k , lenc , lenr
!
! End of declarations rewritten by SPAG
!
!----------------------------------------------------------------------
!     this routine computes maximum numbers of nonzero elements
!     per column/row, minimum numbers of nonzero elements per column/row,
!     and  numbers of zero columns/rows.
!----------------------------------------------------------------------
!     On Entry:
!----------
!     n     = integer column dimension of matrix
!     ja    = integer array containing the row indices of elements in a
!     ia    = integer array containing of length n+1 containing the
!     pointers to the beginning of the columns in arrays a and ja.
!     It is assumed that ia(*) = 1 and ia(n+1) = nzz+1.
!     iao   = similar array for the transpose of the matrix.
!     sym   = logical variable indicating whether or not the matrix is
!     stored in symmetric mode.
!     on return
!----------
!     nzmaxc = max length of columns
!     nzminc = min length of columns
!     nzmaxr = max length of rows
! nzminr = min length of rows
!     nzcol  = number of zero columns
!     nzrow = number of zero rows
!-----------------------------------------------------------------------
!
   Nzmaxc = 0
   Nzminc = N
   Nzmaxr = 0
   Nzminr = N
   Nzcol = 0
   Nzrow = 0
!-----------------------------------------------------------------------
   DO i = 1 , N
      j0 = Ia(i)
      j1 = Ia(i+1)
      j0r = Iao(i)
      j1r = Iao(i+1)
      indiag = 0
      DO k = j0 , j1 - 1
         IF ( Ja(k)==i ) indiag = 1
      ENDDO
!
      lenc = j1 - j0
      lenr = j1r - j0r
!
      IF ( Sym ) lenc = lenc + lenr - indiag
      IF ( lenc<=0 ) Nzcol = Nzcol + 1
      Nzmaxc = max0(Nzmaxc,lenc)
      Nzminc = min0(Nzminc,lenc)
      IF ( lenr<=0 ) Nzrow = Nzrow + 1
      Nzmaxr = max0(Nzmaxr,lenr)
      Nzminr = min0(Nzminr,lenr)
   ENDDO
END SUBROUTINE nonz
!*==diag_domi.f90 processed by SPAG 8.04RA 15:47 21 Feb 2024
!!SPAG Open source Personal, Educational or Academic User Clemson University  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
!-----------------------------------------------------------------------
SUBROUTINE diag_domi(N,Sym,Valued,A,Ja,Ia,Ao,Jao,Iao,Ddomc,Ddomr)
   USE ISO_FORTRAN_ENV                 
   IMPLICIT NONE
!
! Dummy argument declarations rewritten by SPAG
!
   INTEGER , INTENT(IN) :: N
   LOGICAL , INTENT(IN) :: Sym
   LOGICAL , INTENT(IN) :: Valued
   REAL(REAL64) , INTENT(IN) , DIMENSION(*) :: A
   INTEGER , INTENT(IN) , DIMENSION(*) :: Ja
   INTEGER , INTENT(IN) , DIMENSION(N+1) :: Ia
   REAL(REAL64) , INTENT(IN) , DIMENSION(*) :: Ao
   INTEGER , INTENT(IN) , DIMENSION(*) :: Jao
   INTEGER , INTENT(IN) , DIMENSION(N+1) :: Iao
   REAL(REAL64) , INTENT(INOUT) :: Ddomc
   REAL(REAL64) , INTENT(INOUT) :: Ddomr
!
! Local variable declarations rewritten by SPAG
!
   REAL(REAL64) :: aii , dsumc , dsumr
   INTEGER :: i , j , j0 , j1 , k
!
! End of declarations rewritten by SPAG
!
!-----------------------------------------------------------------
!     this routine computes the percentage of weakly diagonally
!     dominant rows/columns
!-----------------------------------------------------------------
!     on entry:
!     ---------
!     n     = integer column dimension of matrix
! a     = real array containing the nonzero elements of the matrix
!     the elements are stored by columns in order
!     (i.e. column i comes before column i+1, but the elements
!     within each column can be disordered).
!     ja    = integer array containing the row indices of elements in a
! ia    = integer array containing of length n+1 containing the
!     pointers to the beginning of the columns in arrays a and ja.
!     It is assumed that ia(*) = 1 and ia(n+1) = nzz+1.
!     ao    = real array containing the nonzero elements of the matrix
!     the elements are stored by rows in order
!     (i.e. row i comes before row i+1, but the elements
!     within each row can be disordered).
! ao,jao, iao,
!     structure for transpose of a
!     sym   = logical variable indicating whether or not the matrix is
!     symmetric.
!     valued= logical equal to .true. if values are provided and .false.
!         if only the pattern of the matrix is provided. (in that
!     case a(*) and ao(*) are dummy arrays.
!
!     ON RETURN
!----------
!     ddomc = percentage of weakly diagonally dominant columns
!     ddomr = percentage of weakly diagonally dominant rows
!-------------------------------------------------------------------
!     locals
!     number of diagonally dominant columns
!     real arithmetic used to avoid problems.. YS. 03/27/01
   Ddomc = 0.0
!     number of diagonally dominant rows
   Ddomr = 0.0
   DO i = 1 , N
      j0 = Ia(i)
      j1 = Ia(i+1) - 1
      IF ( Valued ) THEN
         aii = 0.0D0
         dsumc = 0.0D0
         DO k = j0 , j1
            j = Ja(k)
            IF ( j==i ) THEN
               aii = abs(A(k))
            ELSE
               dsumc = dsumc + abs(A(k))
            ENDIF
         ENDDO
         dsumr = 0.0D0
         IF ( .NOT.Sym ) THEN
            DO k = Iao(i) , Iao(i+1) - 1
               IF ( Jao(k)/=i ) dsumr = dsumr + abs(Ao(k))
            ENDDO
         ELSE
            dsumr = dsumc
         ENDIF
         IF ( dsumc<=aii ) Ddomc = Ddomc + 1.0
         IF ( dsumr<=aii ) Ddomr = Ddomr + 1.0
      ENDIF
   ENDDO
   Ddomr = Ddomr/real(N)
   Ddomc = Ddomc/real(N)
!-----------------------------------------------------------------------
!--------end-of-diag_moni-----------------------------------------------
END SUBROUTINE diag_domi
!*==frobnorm.f90 processed by SPAG 8.04RA 15:47 21 Feb 2024
!!SPAG Open source Personal, Educational or Academic User Clemson University  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
!-----------------------------------------------------------------------
SUBROUTINE frobnorm(N,Sym,A,Ja,Ia,Fnorm)
   USE ISO_FORTRAN_ENV                 
   IMPLICIT NONE
!
! Dummy argument declarations rewritten by SPAG
!
   INTEGER , INTENT(IN) :: N
   LOGICAL , INTENT(IN) :: Sym
   REAL(REAL64) , INTENT(IN) , DIMENSION(*) :: A
   INTEGER , INTENT(IN) , DIMENSION(*) :: Ja
   INTEGER , INTENT(IN) , DIMENSION(N+1) :: Ia
   REAL(REAL64) , INTENT(INOUT) :: Fnorm
!
! Local variable declarations rewritten by SPAG
!
   REAL(REAL64) :: fdiag
   INTEGER :: i , k
!
! End of declarations rewritten by SPAG
!
!--------------------------------------------------------------------------
!     this routine computes the Frobenius norm of A.
!--------------------------------------------------------------------------
!     on entry:
!-----------
! n      = integer colum dimension of matrix
! a      = real array containing the nonzero elements of the matrix
!          the elements are stored by columns in order
!          (i.e. column i comes before column i+1, but the elements
!          within each column can be disordered).
! ja     = integer array containing the row indices of elements in a.
! ia     = integer array containing of length n+1 containing the
!          pointers to the beginning of the columns in arrays a and
!          ja. It is assumed that ia(*)= 1 and ia(n+1)  = nnz +1.
! sym    = logical variable indicating whether or not the matrix is
!          symmetric.
!
! on return
!-----------
! Fnorm  = Frobenius norm of A.
!--------------------------------------------------------------------------
   fdiag = 0.0
   Fnorm = 0.0
   DO i = 1 , N
      DO k = Ia(i) , Ia(i+1) - 1
         IF ( Ja(k)==i ) THEN
            fdiag = fdiag + A(k)**2
         ELSE
            Fnorm = Fnorm + A(k)**2
         ENDIF
      ENDDO
   ENDDO
   IF ( Sym ) THEN
      Fnorm = 2*Fnorm + fdiag
   ELSE
      Fnorm = Fnorm + fdiag
   ENDIF
   Fnorm = sqrt(Fnorm)
END SUBROUTINE frobnorm
!*==ansym.f90 processed by SPAG 8.04RA 15:47 21 Feb 2024
!!SPAG Open source Personal, Educational or Academic User Clemson University  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
!----------------------------------------------------------------------
SUBROUTINE ansym(N,Sym,A,Ja,Ia,Ao,Jao,Iao,Imatch,Av,Fas,Fan)
!---------------------------------------------------------------------
!     this routine computes the Frobenius norm of the symmetric and
!     non-symmetric parts of A, computes number of matching elements
!     in symmetry and relative symmetry match.
!---------------------------------------------------------------------
! on entry:
!----------
! n   = integer column dimension of matrix
! a   = real array containing the nonzero elements of the matrix
!       the elements are stored by columns in order
!       (i.e. column i comes before column i+1, but the elements
!       within each column can be disordered).
! ja  = integer array containing the row indices of elements in a
! ia  = integer array containing of length n+1 containing the
!       pointers to the beginning of the columns in arrays a and ja.
!       It is assumed that ia(*) = 1 and ia(n+1) = nzz+1.
! sym = logical variable indicating whether or not the matrix is
!       symmetric.
! on return
!----------
! fas   = Frobenius norm of symmetric part
! fan   = Frobenius norm of non-symmetric part
! imatch = number of matching elements in symmetry
! av     = relative symmetry match (symmetry = 1)
! ao,jao,iao = transpose of A just as a, ja, ia contains
!              information of A.
!-----------------------------------------------------------------------
   USE ISO_FORTRAN_ENV                 
   IMPLICIT NONE
!
! Dummy argument declarations rewritten by SPAG
!
   INTEGER :: N
   LOGICAL :: Sym
   REAL(REAL64) , DIMENSION(*) :: A
   INTEGER , DIMENSION(*) :: Ja
   INTEGER , DIMENSION(N+1) :: Ia
   REAL(REAL64) , DIMENSION(*) :: Ao
   INTEGER , DIMENSION(*) :: Jao
   INTEGER , DIMENSION(N+1) :: Iao
   INTEGER , INTENT(INOUT) :: Imatch
   REAL(REAL64) , INTENT(OUT) :: Av
   REAL(REAL64) , INTENT(INOUT) :: Fas
   REAL(REAL64) , INTENT(INOUT) :: Fan
!
! Local variable declarations rewritten by SPAG
!
   REAL(REAL64) :: fnorm , st
   INTEGER :: i , j1 , j2 , k1 , k1max , k2 , k2max , nnz
   EXTERNAL csrcsc , frobnorm
!
! End of declarations rewritten by SPAG
!
!-----------------------------------------------------------------------
   st = 0.0
   nnz = Ia(N+1) - Ia(1)
   CALL csrcsc(N,1,1,A,Ja,Ia,Ao,Jao,Iao)
   IF ( .NOT.(Sym) ) THEN
      st = 0.0D0
      Fas = 0.0D0
      Fan = 0.0D0
      Imatch = 0
      DO i = 1 , N
         k1 = Ia(i)
         k2 = Iao(i)
         k1max = Ia(i+1) - 1
         k2max = Iao(i+1) - 1
!
         DO WHILE ( k1<=k1max .AND. k2<=k2max )
!
            j1 = Ja(k1)
            j2 = Jao(k2)
            IF ( j1==j2 ) THEN
               Fas = Fas + (A(k1)+Ao(k2))**2
               Fan = Fan + (A(k1)-Ao(k2))**2
               st = st + A(k1)**2
               Imatch = Imatch + 1
            ENDIF
            k1 = k1 + 1
            k2 = k2 + 1
            IF ( j1<j2 ) k2 = k2 - 1
            IF ( j1>j2 ) k1 = k1 - 1
         ENDDO
      ENDDO
      Fas = 0.25D0*Fas
      Fan = 0.25D0*Fan
   ENDIF
   CALL frobnorm(N,Sym,Ao,Jao,Iao,fnorm)
   IF ( Sym ) THEN
      Imatch = nnz
      Fas = fnorm
      Fan = 0.0D0
   ELSE
      IF ( Imatch==nnz ) THEN
         st = 0.0D0
      ELSE
         st = 0.5D0*(fnorm**2-st)
         IF ( st<0.0D0 ) st = 0.0D0
      ENDIF
      Fas = sqrt(Fas+st)
      Fan = sqrt(Fan+st)
   ENDIF
   Av = real(Imatch)/real(nnz)
END SUBROUTINE ansym
!*==distaij.f90 processed by SPAG 8.04RA 15:47 21 Feb 2024
!!SPAG Open source Personal, Educational or Academic User Clemson University  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
!------end-of-ansym-----------------------------------------------------
!-----------------------------------------------------------------------
SUBROUTINE distaij(N,Nnz,Ja,Ia,Dist,Std)
   USE ISO_FORTRAN_ENV                 
   IMPLICIT NONE
!
! Dummy argument declarations rewritten by SPAG
!
   INTEGER , INTENT(IN) :: N
   INTEGER , INTENT(IN) :: Nnz
   INTEGER , INTENT(IN) , DIMENSION(*) :: Ja
   INTEGER , INTENT(IN) , DIMENSION(N+1) :: Ia
   REAL(REAL64) , INTENT(INOUT) :: Dist
   REAL(REAL64) , INTENT(INOUT) :: Std
!
! Local variable declarations rewritten by SPAG
!
   INTEGER :: i , j , j0 , j1 , k
!
! End of declarations rewritten by SPAG
!
!-----------------------------------------------------------------------
!     this routine computes the average distance of a(i,j) from diag and
!     standard deviation  for this average.
!-----------------------------------------------------------------------
! On entry :
!-----------
! n     = integer. column dimension of matrix
! nnz   = number of nonzero elements of matrix
! ja    = integer array containing the row indices of elements in a
! ia    = integer array containing of length n+1 containing the
!         pointers to the beginning of the columns in arrays a and ja.
!         It is assumed that ia(*) = 1 and ia(n+1) = nzz+1.
! sym   = logical variable indicating whether or not the matrix is
!         symmetric.
! on return
!----------
! dist  = average distance of a(i,j) from diag.
! std   = standard deviation for above average.
!-----------------------------------------------------------------------
!
! distance of an element from diagonal.
!
   Dist = 0.0
   Std = 0.0
   DO i = 1 , N
      j0 = Ia(i)
      j1 = Ia(i+1) - 1
      DO k = j0 , j1
         j = Ja(k)
         Dist = Dist + real(iabs(j-i))
      ENDDO
   ENDDO
   Dist = Dist/real(Nnz)
   DO i = 1 , N
      DO k = Ia(i) , Ia(i+1) - 1
         Std = Std + (Dist-real(iabs(Ja(k)-i)))**2
      ENDDO
   ENDDO
   Std = sqrt(Std/real(Nnz))
END SUBROUTINE distaij
!*==skyline.f90 processed by SPAG 8.04RA 15:47 21 Feb 2024
!!SPAG Open source Personal, Educational or Academic User Clemson University  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
!-----------------------------------------------------------------------
SUBROUTINE skyline(N,Sym,Ja,Ia,Jao,Iao,Nsky)
   USE ISO_FORTRAN_ENV                 
   IMPLICIT NONE
!
! Dummy argument declarations rewritten by SPAG
!
   INTEGER , INTENT(IN) :: N
   LOGICAL , INTENT(IN) :: Sym
   INTEGER , INTENT(IN) , DIMENSION(*) :: Ja
   INTEGER , INTENT(IN) , DIMENSION(N+1) :: Ia
   INTEGER , INTENT(IN) , DIMENSION(*) :: Jao
   INTEGER , INTENT(IN) , DIMENSION(N+1) :: Iao
   INTEGER , INTENT(OUT) :: Nsky
!
! Local variable declarations rewritten by SPAG
!
   INTEGER :: i , j0 , j0r , jminc , jminr , nskyl , nskyu
!
! End of declarations rewritten by SPAG
!
!-------------------------------------------------------------------
! this routine computes the number of nonzeros in the skyline storage.
!-------------------------------------------------------------------
!
! On entry :
!-----------
! n     = integer. column dimension of matrix
! ja    = integer array containing the row indices of elements in a
! ia    = integer array containing of length n+1 containing the
!         pointers to the beginning of the columns in arrays a and ja.
!         It is assumed that ia(*) = 1 and ia(n+1) = nzz+1.
! iao   = integer array containing of length n+1 containing the
!         pointers to the beginning of the rows in arrays ao and jao.
!         It is assumed that iao(*) = 1 and iao(n+1) = nzz+1.
! jao   = integer array containing the column indices of elements in ao.
! sym   = logical variable indicating whether or not the matrix is
!         symmetric.
! on return
!----------
! nsky  = number of nonzeros in skyline storage
!-------------------------------------------------------------------
!
! nskyu = skyline storage for upper part
   nskyu = 0
! nskyl = skyline storage for lower part
   nskyl = 0
   DO i = 1 , N
      j0 = Ia(i)
      j0r = Iao(i)
 
      jminc = Ja(j0)
      jminr = Jao(j0r)
      IF ( Sym ) jminc = jminr
 
      nskyl = nskyl + i - jminr + 1
      nskyu = nskyu + i - jminc + 1
 
   ENDDO
   Nsky = nskyl + nskyu - N
   IF ( Sym ) Nsky = nskyl
END SUBROUTINE skyline
!*==distdiag.f90 processed by SPAG 8.04RA 15:47 21 Feb 2024
!!SPAG Open source Personal, Educational or Academic User Clemson University  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
!-----------------------------------------------------------------------
SUBROUTINE distdiag(Nrow,Ncol,Ja,Ia,Dist)
   USE ISO_FORTRAN_ENV                 
   IMPLICIT NONE
!
! Dummy argument declarations rewritten by SPAG
!
   INTEGER , INTENT(IN) :: Nrow
   INTEGER , INTENT(IN) :: Ncol
   INTEGER , INTENT(IN) , DIMENSION(*) :: Ja
   INTEGER , INTENT(IN) , DIMENSION(Nrow+1) :: Ia
   INTEGER , INTENT(INOUT) , DIMENSION(*) :: Dist
!
! Local variable declarations rewritten by SPAG
!
   INTEGER :: i , j , k , k1 , k2 , n2
!
! End of declarations rewritten by SPAG
!
!----------------------------------------------------------------------
! this routine computes the numbers of elements in each diagonal.
!----------------------------------------------------------------------
! On entry :
!-----------
! nrow  = integer. row dimension of matrix
! ncol  = integer. column dimension of matrix
! ja    = integer array containing the row indices of elements in a
! ia    = integer array containing of length n+1 containing the
!         pointers to the beginning of the columns in arrays a and ja.
!         It is assumed that ia(*) = 1 and ia(n+1) = nzz+1.
! on return
!----------
! dist  = integer array containing the numbers of elements in each of
!         the nrow+ncol-1 diagonals of A. dist(k) contains the
!         number of elements in diagonal '-nrow+k'.  k ranges from
!         1 to (nrow+ncol-1).
!----------------------------------------------------------------------
!
   n2 = Nrow + Ncol - 1
   DO i = 1 , n2
      Dist(i) = 0
   ENDDO
   DO i = 1 , Nrow
      k1 = Ia(i)
      k2 = Ia(i+1) - 1
      DO k = k1 , k2
         j = Ja(k)
         Dist(Nrow+j-i) = Dist(Nrow+j-i) + 1
      ENDDO
   ENDDO
END SUBROUTINE distdiag
!*==bandpart.f90 processed by SPAG 8.04RA 15:47 21 Feb 2024
!!SPAG Open source Personal, Educational or Academic User Clemson University  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
!-----------------------------------------------------------------------
SUBROUTINE bandpart(N,Ia,Dist,Nper,Band)
   USE ISO_FORTRAN_ENV                 
   IMPLICIT NONE
!
! Dummy argument declarations rewritten by SPAG
!
   INTEGER , INTENT(IN) :: N
   INTEGER , INTENT(IN) , DIMENSION(N+1) :: Ia
   INTEGER , INTENT(IN) , DIMENSION(*) :: Dist
   INTEGER , INTENT(IN) :: Nper
   INTEGER , INTENT(INOUT) :: Band
!
! Local variable declarations rewritten by SPAG
!
   INTEGER :: iacc , j , nnz
!
! End of declarations rewritten by SPAG
!
!-------------------------------------------------------------------------
! this routine computes the bandwidth of the banded matrix, which contains
! 'nper' percent of the original matrix.
!-------------------------------------------------------------------------
! On entry :
!-----------
! n     = integer. column dimension of matrix
! a, ja, ia = matrix in CSR format. only ia used/
! ia    = integer array containing of length n+1 containing the
!         pointers to the beginning of the columns in arrays a and ja.
!         It is assumed that ia(*) = 1 and ia(n+1) = nzz+1.
! dist  = integer array containing the numbers of elements in the
!         matrix with different distance of row indices and column
!         indices.
! nper  = percentage of matrix  within the bandwidth
! on return
!----------
! band  = the width of the bandwidth
!----------------------------------------------------------------------
   nnz = Ia(N+1) - Ia(1)
   iacc = Dist(N)
   Band = 0
   j = 0
   SPAG_Loop_1_1: DO
      j = j + 1
      iacc = iacc + Dist(N+j) + Dist(N-j)
      IF ( iacc*100<=nnz*Nper ) THEN
         Band = Band + 1
         CYCLE
      ENDIF
      EXIT SPAG_Loop_1_1
   ENDDO SPAG_Loop_1_1
END SUBROUTINE bandpart
!*==n_imp_diag.f90 processed by SPAG 8.04RA 15:47 21 Feb 2024
!!SPAG Open source Personal, Educational or Academic User Clemson University  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
!-----------------------------------------------------------------------
SUBROUTINE n_imp_diag(N,Nnz,Dist,Ipar1,Ndiag,Ioff,Dcount)
   USE ISO_FORTRAN_ENV                 
   IMPLICIT NONE
!
! Dummy argument declarations rewritten by SPAG
!
   INTEGER , INTENT(IN) :: N
   INTEGER , INTENT(IN) :: Nnz
   INTEGER , INTENT(INOUT) , DIMENSION(*) :: Dist
   INTEGER , INTENT(IN) :: Ipar1
   INTEGER , INTENT(INOUT) :: Ndiag
   INTEGER , INTENT(OUT) , DIMENSION(*) :: Ioff
   REAL(REAL64) , INTENT(OUT) , DIMENSION(*) :: Dcount
!
! Local variable declarations rewritten by SPAG
!
   INTEGER :: i , ii , itot , j , jmax , k , n2
!
! End of declarations rewritten by SPAG
!
!-----------------------------------------------------------------------
!     this routine computes the most important diagonals.
!-----------------------------------------------------------------------
!
! On entry :
!-----------
! n     = integer. column dimension of matrix
! nnz   = number of nonzero elements of matrix
! dist  = integer array containing the numbers of elements in the
!         matrix with different distance of row indices and column
!         indices. ipar1 = percentage of nonzero  elements of A that
!         a diagonal should have in order to be an important diagonal
! on return
!----------
! ndiag = number of the most important diagonals
! ioff  = the offsets with respect to the main diagonal
! dcount= the accumulated percentages
!-----------------------------------------------------------------------
   n2 = N + N - 1
   Ndiag = 10
   Ndiag = min0(n2,Ndiag)
   itot = 0
   ii = 0
   SPAG_Loop_1_1: DO
!     sort diagonals by decreasing order of weights.
      jmax = 0
      i = 1
      DO k = 1 , n2
         j = Dist(k)
         IF ( j>=jmax ) THEN
            i = k
            jmax = j
         ENDIF
      ENDDO
!     permute ----
!     save offsets and accumulated count if diagonal is acceptable
!     (if it has at least ipar1*nnz/100 nonzero elements)
!     quite if no more acceptable diagonals --
!
      IF ( jmax*100<Ipar1*Nnz ) THEN
         Ndiag = ii
         EXIT SPAG_Loop_1_1
      ELSE
         ii = ii + 1
         Ioff(ii) = i - N
         Dist(i) = -jmax
         itot = itot + jmax
         Dcount(ii) = real(100*itot)/real(Nnz)
         IF ( ii>=Ndiag ) THEN
            Ndiag = ii
            EXIT SPAG_Loop_1_1
         ENDIF
      ENDIF
   ENDDO SPAG_Loop_1_1
!-----------------------------------------------------------------------
END SUBROUTINE n_imp_diag
!*==nonz_lud.f90 processed by SPAG 8.04RA 15:47 21 Feb 2024
!!SPAG Open source Personal, Educational or Academic User Clemson University  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
!-----------------------------------------------------------------------
SUBROUTINE nonz_lud(N,Ja,Ia,Nlower,Nupper,Ndiag)
   USE ISO_FORTRAN_ENV                 
   IMPLICIT NONE
!
! Dummy argument declarations rewritten by SPAG
!
   INTEGER , INTENT(IN) :: N
   INTEGER , INTENT(IN) , DIMENSION(*) :: Ja
   INTEGER , INTENT(IN) , DIMENSION(N+1) :: Ia
   INTEGER , INTENT(OUT) :: Nlower
   INTEGER , INTENT(INOUT) :: Nupper
   INTEGER , INTENT(INOUT) :: Ndiag
!
! Local variable declarations rewritten by SPAG
!
   INTEGER :: i , j , k
!
! End of declarations rewritten by SPAG
!
!-----------------------------------------------------------------------
! this routine computes the number of nonzero elements in strict lower
! part, strict upper part, and main diagonal.
!-----------------------------------------------------------------------
!
! On entry :
!-----------
! n     = integer. column dimension of matrix
! ja    = integer array containing the row indices of elements in a
! ia    = integer array containing of length n+1 containing the
!         pointers to the beginning of the columns in arrays a and ja.
!         It is assumed that ia(*) = 1 and ia(n+1) = nzz+1.
! on return
!----------
! nlower= number of nonzero elements in strict lower part
! nupper= number of nonzero elements in strict upper part
! ndiag = number of nonzero elements in main diagonal
!-------------------------------------------------------------------
!
! number of nonzero elements in upper part
!
   Nupper = 0
   Ndiag = 0
 
   DO i = 1 , N
!     indiag = nonzero diagonal element indicator
      DO k = Ia(i) , Ia(i+1) - 1
         j = Ja(k)
         IF ( j<i ) Nupper = Nupper + 1
         IF ( j==i ) Ndiag = Ndiag + 1
      ENDDO
   ENDDO
   Nlower = Ia(N+1) - 1 - Nupper - Ndiag
END SUBROUTINE nonz_lud
!*==avnz_col.f90 processed by SPAG 8.04RA 15:47 21 Feb 2024
!!SPAG Open source Personal, Educational or Academic User Clemson University  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
!-----------------------------------------------------------------------
SUBROUTINE avnz_col(N,Ia,Av,St)
   USE ISO_FORTRAN_ENV                 
   IMPLICIT NONE
!
! Dummy argument declarations rewritten by SPAG
!
   INTEGER , INTENT(IN) :: N
   INTEGER , INTENT(IN) , DIMENSION(N+1) :: Ia
   REAL(REAL64) , INTENT(INOUT) :: Av
   REAL(REAL64) , INTENT(INOUT) :: St
!
! Local variable declarations rewritten by SPAG
!
   INTEGER :: i , j0 , j1 , lenc
!
! End of declarations rewritten by SPAG
!
!---------------------------------------------------------------------
!     this routine computes average number of nonzero elements/column and
!     standard deviation for this average
!---------------------------------------------------------------------
!
! On entry :
!-----------
! n     = integer. column dimension of matrix
! ja    = integer array containing the row indices of elements in a
! ia    = integer array containing of length n+1 containing the
!         pointers to the beginning of the columns in arrays a and ja.
!         It is assumed that ia(*) = 1 and ia(n+1) = nzz+1.
! ndiag = number of the most important diagonals
! On return
!----------
! av    = average number of nonzero elements/column
! st    = standard deviation for this average
! Notes
!---------
! standard deviation will not be correct for symmetric storage.
!----------------------------------------------------------------------
!     standard deviatioan for the average
   St = 0.0D0
!     average and standard deviation
!
   Av = real(Ia(N+1)-1)/real(N)
!
!     will be corrected later.
!
   DO i = 1 , N
      j0 = Ia(i)
      j1 = Ia(i+1) - 1
!     indiag = nonzero diagonal element indicator
      lenc = j1 + 1 - j0
      St = St + (real(lenc)-Av)**2
   ENDDO
!
   St = sqrt(St/real(N))
END SUBROUTINE avnz_col
!*==vbrinfo.f90 processed by SPAG 8.04RA 15:47 21 Feb 2024
!!SPAG Open source Personal, Educational or Academic User Clemson University  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
SUBROUTINE vbrinfo(Nr,Nc,Kvstr,Kvstc,Ia,Ka,Iwk,Iout)
   IMPLICIT NONE
!
! Dummy argument declarations rewritten by SPAG
!
   INTEGER :: Nr
   INTEGER :: Nc
   INTEGER , DIMENSION(Nr+1) :: Kvstr
   INTEGER , DIMENSION(Nc+1) :: Kvstc
   INTEGER , INTENT(IN) , DIMENSION(Nr+1) :: Ia
   INTEGER , INTENT(IN) , DIMENSION(*) :: Ka
   INTEGER , INTENT(INOUT) , DIMENSION(Nr+Nc+2+Nr) :: Iwk
   INTEGER , INTENT(IN) :: Iout
!
! Local variable declarations rewritten by SPAG
!
   INTEGER , DIMENSION(10) :: bsiz , freq
   INTEGER :: i , j , max , n , nb , neq , nnz , nnzb , num
   CHARACTER(101) :: tmpst
   EXTERNAL kvstmerge
!
! End of declarations rewritten by SPAG
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!     Print info on matrix in variable block row format.
!-----------------------------------------------------------------------
!     On entry:
!--------------
!     nr,nc   = matrix block row and block column dimension
!     kvstr   = first row number for each block row
!     kvstc   = first column number for each block column
!     ia,ja,ka,a = input matrix in VBR format
!     iout    = unit number for printed output
!
!     On return:
!---------------
!     Printed output to unit number specified in iout.  If a non-square
!     matrix is provided, the analysis will be performed on the block
!     rows, otherwise a row/column conformal partitioning will be used.
!
!     Work space:
!----------------
!     iwk(1:nb+1) = conformal block partitioning
!        (nb is unknown at start but is no more than nr+nc)
!     iwk(nb+2:nb+2+nr) = frequency of each blocksize
!     The workspace is not assumed to be initialized to zero, nor is it
!     left that way.
!
!-----------------------------------------------------------------------
!-----local variables
!-----------------------------------------------------------------------
   n = Kvstr(Nr+1) - 1
   nnz = Ka(Ia(Nr+1)) - Ka(1)
   nnzb = Ia(Nr+1) - Ia(1)
   WRITE (Iout,99001)
   WRITE (Iout,99002) n , nnz , real(nnz)/real(n)
99002 FORMAT (6x,' *  Number of rows                                   = ',i10,'  *'/6x,                                           &
             &' *  Number of nonzero elements                       = ',i10,'  *'/6x,                                              &
             &' *  Average number of nonzero elements/Row           = ',f10.4,'  *')
   WRITE (Iout,99003) Nr , nnzb , real(nnzb)/real(Nr)
99003 FORMAT (6x,' *  Number of block rows                             = ',i10,'  *'/6x,                                           &
             &' *  Number of nonzero blocks                         = ',i10,'  *'/6x,                                              &
             &' *  Average number of nonzero blocks/Block row       = ',f10.4,'  *')
!-----if non-square matrix, do analysis on block rows,
!     else do analysis on conformal partitioning
   IF ( Kvstr(Nr+1)/=Kvstc(Nc+1) ) THEN
      WRITE (Iout,99004)
99004 FORMAT (6x,' *  Non-square matrix.                                 ','            *'/6x,                                     &
             &' *  Performing analysis on block rows.                 ','            *')
      DO i = 1 , Nr + 1
         Iwk(i) = Kvstr(i)
      ENDDO
      nb = Nr
   ELSE
      CALL kvstmerge(Nr,Kvstr,Nc,Kvstc,nb,Iwk)
      IF ( (Nr/=Nc) .OR. (Nc/=nb) ) WRITE (Iout,99005) nb
99005 FORMAT (6x,' *  Non row-column conformal partitioning supplied.    ','            *'/6x,                                     &
             &' *  Using conformal partitioning.  Number of bl rows = ',i10,'  *')
   ENDIF
!-----accumulate frequencies of each blocksize
   max = 1
   Iwk(1+nb+2) = 0
   DO i = 1 , nb
      neq = Iwk(i+1) - Iwk(i)
      IF ( neq>max ) THEN
         DO j = max + 1 , neq
            Iwk(j+nb+2) = 0
         ENDDO
         max = neq
      ENDIF
      Iwk(neq+nb+2) = Iwk(neq+nb+2) + 1
   ENDDO
!-----store largest 10 of these blocksizes
   num = 0
   DO i = max , 1 , -1
      IF ( (Iwk(i+nb+2)/=0) .AND. (num<10) ) THEN
         num = num + 1
         bsiz(num) = i
         freq(num) = Iwk(i+nb+2)
      ENDIF
   ENDDO
!-----print information about blocksizes
   WRITE (Iout,99006) num
99006 FORMAT (6x,' *  Number of different blocksizes                   = ',i10,'  *')
   WRITE (tmpst,'(10i6)') (bsiz(j),j=1,num)
   WRITE (Iout,99007) num , tmpst
99007 FORMAT (6x,' *  The ',i2,' largest dimension nodes',' have dimension    : ',10x,'  *',/,6x,' *',a61,3x,' *')
   WRITE (tmpst,'(10i6)') (freq(j),j=1,num)
   WRITE (Iout,99008) tmpst
99008 FORMAT (6x,' *  The frequency of nodes these ','dimensions are      : ',10x,'  *',/,6x,' *',a61,3x,' *')
   WRITE (Iout,99001)
!-----------------------------------------------------------------------
99001 FORMAT (6x,' *',65('-'),'*')
!---------------------------------
END SUBROUTINE vbrinfo
