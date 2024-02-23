!*==matrf2.f90 processed by SPAG 8.04RA 15:47 21 Feb 2024
!!SPAG Open source Personal, Educational or Academic User Clemson University  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
SUBROUTINE matrf2(M,N,C,Index,Alpha,Nn,Nz,A,Snr,Rnr,Fejlm)
   USE ISO_FORTRAN_ENV                 
   IMPLICIT NONE
!
! Dummy argument declarations rewritten by SPAG
!
   INTEGER , INTENT(IN) :: Nn
   INTEGER , INTENT(IN) :: M
   INTEGER , INTENT(IN) :: N
   INTEGER , INTENT(IN) :: C
   INTEGER , INTENT(IN) :: Index
   REAL(REAL64) , INTENT(INOUT) :: Alpha
   INTEGER , INTENT(INOUT) :: Nz
   REAL(REAL64) , INTENT(OUT) , DIMENSION(Nn) :: A
   INTEGER , INTENT(OUT) , DIMENSION(Nn) :: Snr
   INTEGER , INTENT(OUT) , DIMENSION(Nn) :: Rnr
   INTEGER , INTENT(OUT) :: Fejlm
!
! Local variable declarations rewritten by SPAG
!
   REAL(REAL64) :: alpha1
   INTEGER :: i , index1 , j , j1 , k , m1 , m2 , n2 , nz1 , rr1 , rr2 , rr3
!
! End of declarations rewritten by SPAG
!
!--------------------------------------------------------------------
!
!   PURPOSE
!   -------
!   The subroutine generates sparse (rectangular or square) matrices.
!   The dimensions of the matrix and the average number of nonzero
!   elements per row can be specified by the user. Moreover, the user
!   can also change the sparsity pattern and the condition number of the
!   matrix. The non-zero elements of the desired matrix will be
!   accumulated (in an arbitrary order) in the first NZ positions of
!   array A. The column and the row numbers of the non-zero element
!   stored in A(I), I=1,...,NZ, will be found in SNR(I) and RNR(I),
!   respectively. The matrix generated by this subroutine is of the
!   class F(M,N,C,R,ALPHA) (see reference).
!
!   Note: If A is the sparse matrix of type F(M,N,C,R,ALPHA), then
!
!           min|A(i,j)| = 1/ALPHA,
!
!           max|A(i,j)| = max(INDEX*N - N,10*ALPHA).
!
!
!   CONTRIBUTOR: Ernest E. Rothman
!                Cornell Theory Center/Cornell National Supercomputer
!                Facility.
!                e-mail address: BITNET:   eer@cornellf
!                                INTERNET: eer@cornellf.tn.cornell.edu
!
!   minor modifications by Y. Saad. April 26, 1990.
!
!   Note: This subroutine has been copied from the following reference.
!         The allowable array sizes have been changed.
!
!   REFERENCE: Zlatev, Zahari; Schaumburg, Kjeld; Wasniewski, Jerzy;
!      "A testing Scheme for Subroutines Solving Large Linear Problems",
!      Computers and Chemistry, Vol. 5, No. 2-3, pp. 91-100, 1981.
!
!
!   INPUT PARAMETERS
!   ----------------
!   M    - Integer. The number of rows in the desired matrix.
!          N < M+1 < 9000001 must be specified.
!
!   N    - Integer. The number of columns in the desired matrix.
!          21 < N < 9000001 must be specified.
!
!   C    - Integer. The sparsity pattern can be changed by means of this
!          parameter.  10 < C < N-10  must be specified.
!
!   INDEX - Integer.  The average number of non-zero elements per row in
!           the matrix will be equal to INDEX.
!           1 < INDEX < N-C-8 must be specified.
!
!   ALPHA - Real. The condition number of the matrix can be changed
!           BY THIS PARAMETER. ALPHA > 0.0 MUST BE SPECIFIED.
!           If ALPHA is approximately equal to 1.0 then the generated
!           matrix is well-conditioned. Large values of ALPHA will
!           usually produce ill-conditioned matrices. Note that no
!           round-off errors during the computations in this subroutine
!           are made if ALPHA = 2**I (where I is an arbitrary integer
!           which produces numbers in the machine range).
!
!   NN    - Integer. The length of arrays A, RNR, and SNR (see below).
!           INDEX*M+109 < NN < 9000001 must be specified.
!
!
!   OUTPUT PARAMETERS
!   -----------------
!   NZ    - Integer. The number of non-zero elements in the matrix.
!
!   A(NN) - Real array. The non-zero elements of the matrix generated
!           are accumulated in the first NZ locations of array A.
!
!   SNR(NN) - INTEGER array. The column number of the non-zero element
!           kept in A(I), I=1,...NZ, is stored in SNR(I).
!
!   RNR(NN) - Integer array. The row number of the non-zero element
!           kept in A(I), I=1,...NZ, is stored in RNR(I).
!
!   FEJLM - Integer. FEJLM=0 indicates that the call is successful.
!           Error diagnostics are given by means of positive values of
!           this parameter as follows:
!             FEJLM = 1    -  N       is out of range.
!             FEJLM = 2    -  M       is out of range.
!             FEJLM = 3    -  C       is out of range.
!             FEJLM = 4    -  INDEX   is out of range.
!             FEJLM = 5    -  NN      is out of range.
!             FEJLM = 7    -  ALPHA   is out of range.
!
!
!
!
   INTEGER :: spag_nextblock_1
   spag_nextblock_1 = 1
   SPAG_DispatchLoop_1: DO
      SELECT CASE (spag_nextblock_1)
      CASE (1)
         m1 = M
         Fejlm = 0
         nz1 = Index*M + 110
         k = 1
         alpha1 = Alpha
         index1 = Index - 1
!
!  Check the parameters.
!
         IF ( N>=22 ) THEN
            IF ( N<=9000000 ) THEN
               IF ( M<N ) THEN
                  spag_nextblock_1 = 2
                  CYCLE SPAG_DispatchLoop_1
               ENDIF
               IF ( M>9000000 ) THEN
                  spag_nextblock_1 = 2
                  CYCLE SPAG_DispatchLoop_1
               ENDIF
               IF ( C<11 ) THEN
                  spag_nextblock_1 = 3
                  CYCLE SPAG_DispatchLoop_1
               ENDIF
               IF ( N-C<11 ) THEN
                  spag_nextblock_1 = 3
                  CYCLE SPAG_DispatchLoop_1
               ENDIF
               IF ( Index<1 ) THEN
                  Fejlm = 4
               ELSEIF ( N-C-Index<9 ) THEN
                  Fejlm = 4
               ENDIF
               IF ( Nn<nz1 ) THEN
                  spag_nextblock_1 = 4
                  CYCLE SPAG_DispatchLoop_1
               ENDIF
               IF ( Nn>9000000 ) THEN
                  spag_nextblock_1 = 4
                  CYCLE SPAG_DispatchLoop_1
               ENDIF
               IF ( Alpha>0.0 ) THEN
!
!  End of the error check. Begin to generate the non-zero elements of
!  the required matrix.
!
                  DO i = 1 , N
                     A(i) = 1.0D0
                     Snr(i) = i
                     Rnr(i) = i
                  ENDDO
                  Nz = N
                  j1 = 1
                  IF ( index1/=0 ) THEN
                     DO j = 1 , index1
                        j1 = -j1
                        DO i = 1 , N
                           A(Nz+i) = dfloat(j1*j*i)
                           IF ( i+C+j-1<=N ) Snr(Nz+i) = i + C + j - 1
                           IF ( i+C+j-1>N ) Snr(Nz+i) = C + i + j - 1 - N
                           Rnr(Nz+i) = i
                        ENDDO
                        Nz = Nz + N
                     ENDDO
                  ENDIF
                  rr1 = 10
                  rr2 = Nz
                  rr3 = 1
                  DO
                     DO i = 1 , rr1
                        A(rr2+i) = Alpha*dfloat(i)
                        Snr(rr2+i) = N - rr1 + i
                        Rnr(rr2+i) = rr3
                     ENDDO
                     IF ( rr1==1 ) THEN
                        Nz = Nz + 55
                        DO
                           m1 = m1 - N
                           Alpha = 1.0D0/Alpha
                           IF ( m1<=0 ) THEN
                              Alpha = 1.0D0/alpha1
                              rr1 = 1
                              rr2 = Nz
                              DO
                                 DO i = 1 , rr1
                                    A(rr2+i) = Alpha*dfloat(rr1+1-i)
                                    Snr(rr2+i) = i
                                    Rnr(rr2+i) = N - 10 + rr1
                                 ENDDO
                                 IF ( rr1==10 ) THEN
                                    Nz = Nz + 55
                                    Alpha = alpha1
                                    EXIT SPAG_DispatchLoop_1
                                 ELSE
                                    rr2 = rr2 + rr1
                                    rr1 = rr1 + 1
                                 ENDIF
                              ENDDO
                           ELSE
                              n2 = k*N
                              IF ( m1>=N ) m2 = N
                              IF ( m1<N ) m2 = m1
                              DO i = 1 , m2
                                 A(Nz+i) = Alpha*dfloat(k+1)
                                 Snr(Nz+i) = i
                                 Rnr(Nz+i) = n2 + i
                              ENDDO
                              Nz = Nz + m2
                              IF ( index1/=0 ) THEN
                                 j1 = 1
                                 DO j = 1 , index1
                                    j1 = -j1
                                    DO i = 1 , m2
                                       A(Nz+i) = Alpha*dfloat(j*j1)*(dfloat((k+1)*i)+1.0D0)
                                       IF ( i+C+j-1<=N ) Snr(Nz+i) = i + C + j - 1
                                       IF ( i+C+j-1>N ) Snr(Nz+i) = C + i + j - 1 - N
                                       Rnr(Nz+i) = n2 + i
                                    ENDDO
                                    Nz = Nz + m2
                                 ENDDO
                              ENDIF
                              k = k + 1
                           ENDIF
                        ENDDO
                     ELSE
                        rr2 = rr2 + rr1
                        rr1 = rr1 - 1
                        rr3 = rr3 + 1
                     ENDIF
                  ENDDO
               ELSE
                  Fejlm = 6
                  RETURN
               ENDIF
            ENDIF
         ENDIF
         Fejlm = 1
         RETURN
      CASE (2)
         Fejlm = 2
         RETURN
      CASE (3)
         Fejlm = 3
         RETURN
      CASE (4)
         Fejlm = 5
         RETURN
      END SELECT
   ENDDO SPAG_DispatchLoop_1
END SUBROUTINE matrf2
!*==dcn.f90 processed by SPAG 8.04RA 15:47 21 Feb 2024
!!SPAG Open source Personal, Educational or Academic User Clemson University  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
SUBROUTINE dcn(Ar,Ia,Ja,N,Ne,Ic,Nn,Ierr)
   USE ISO_FORTRAN_ENV                 
   IMPLICIT NONE
!
! Dummy argument declarations rewritten by SPAG
!
   INTEGER , INTENT(IN) :: Nn
   REAL(REAL64) , INTENT(OUT) , DIMENSION(Nn) :: Ar
   INTEGER , INTENT(OUT) , DIMENSION(Nn) :: Ia
   INTEGER , INTENT(OUT) , DIMENSION(Nn) :: Ja
   INTEGER , INTENT(IN) :: N
   INTEGER , INTENT(INOUT) :: Ne
   INTEGER , INTENT(IN) :: Ic
   INTEGER , INTENT(OUT) :: Ierr
!
! Local variable declarations rewritten by SPAG
!
   INTEGER :: i , icount , ilast , it , j
!
! End of declarations rewritten by SPAG
!
!-----------------------------------------------------------------------
!
!   PURPOSE
!   -------
!   The subroutine generates sparse (square) matrices of the type
!   D(N,C).  This type of matrix has the following characteristics:
!   1's in the diagonal, three bands at the distance C above the
!   diagonal (and reappearing cyclicly under it), and a 10 x 10
!   triangle of elements in the upper right-hand corner.
!   Different software libraries require different storage schemes.
!   This subroutine generates the matrix in the storage  by
!   indices mode.
!
!
!   Note: If A is the sparse matrix of type D(N,C), then
!
!       min|A(i,j)| = 1,     max|A(i,j)| = max(1000,N + 1)
!
!
!
!   CONTRIBUTOR: Ernest E. Rothman
!                Cornell Theory Center/Cornell National Supercomputer
!                Facility.
!                e-mail address: BITNET:   eer@cornellf
!                                INTERNET: eer@cornellf.tn.cornell.edu
!
!
!   REFERENCE
!   ---------
!   1) Zlatev, Zahari; Schaumburg, Kjeld; Wasniewski, Jerzy;
!      "A Testing Scheme for Subroutines Solving Large Linear Problems",
!       Computers and Chemistry, Vol. 5, No. 2-3, pp. 91-100, 1981.
!   2) Osterby, Ole and Zletev, Zahari;
!      "Direct Methods for Sparse Matrices";
!       Springer-Verlag 1983.
!
!
!
!   INPUT PARAMETERS
!   ----------------
!   N    - Integer. The size of the square matrix.
!          N > 13 must be specified.
!
!   NN   - Integer. The dimension of integer arrays IA and JA and
!          real array AR. Must be at least NE.
!
!   IC   - Integer. The sparsity pattern can be changed by means of this
!          parameter.  0 < IC < N-12  must be specified.
!
!
!   OUTPUT PARAMETERS
!   -----------------
!   NE   - Integer. The number of nonzero elements in the sparse matrix
!          of the type D(N,C). NE = 4*N + 55.
!
!   AR(NN) - Real array. (Double precision)
!            Stored entries of a sparse matrix to be generated by this
!            subroutine.
!            NN is greater then or equal to, NE, the number of
!            nonzeros including a mandatory diagonal entry for
!            each row. Entries are stored by indices.
!
!   IA(NN) - Integer array.
!            Pointers to specify rows for the stored nonzero entries
!            in AR.
!
!   JA(NN) - Integer array.
!            Pointers to specify columns for the stored nonzero entries
!            in AR.
!
!   IERR   - Error parameter is returned as zero on successful
!             execution of the subroutine.
!             Error diagnostics are given by means of positive values
!             of this parameter as follows:
!             IERR = 1    -  N       is out of range.
!             IERR = 2    -  IC      is out of range.
!             IERR = 3    -  NN      is out of range.
!
!----------------------------------------------------------------------
!
   Ierr = 0
!
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!heck the input parameters:
!
   IF ( N<=13 ) THEN
      Ierr = 1
      RETURN
   ENDIF
   IF ( Ic<=0 .OR. Ic>=N-12 ) THEN
      Ierr = 2
      RETURN
   ENDIF
   Ne = 4*N + 55
   IF ( Nn<Ne ) THEN
      Ierr = 3
      RETURN
   ENDIF
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
! Begin to generate the nonzero elements as well as the row and column
! pointers:
!
   DO i = 1 , N
      Ar(i) = 1.0D0
      Ia(i) = i
      Ja(i) = i
   ENDDO
   ilast = N
   DO i = 1 , N - Ic
      it = ilast + i
      Ar(it) = 1.0 + dfloat(i)
      Ia(it) = i
      Ja(it) = i + Ic
   ENDDO
   ilast = ilast + N - Ic
   DO i = 1 , N - Ic - 1
      it = ilast + i
      Ar(it) = -dfloat(i)
      Ia(it) = i
      Ja(it) = i + Ic + 1
   ENDDO
   ilast = ilast + N - Ic - 1
   DO i = 1 , N - Ic - 2
      it = ilast + i
      Ar(it) = 16.0D0
      Ia(it) = i
      Ja(it) = i + Ic + 2
   ENDDO
   ilast = ilast + N - Ic - 2
   icount = 0
   DO j = 1 , 10
      DO i = 1 , 11 - j
         icount = icount + 1
         it = ilast + icount
         Ar(it) = 100.0D0*dfloat(j)
         Ia(it) = i
         Ja(it) = N - 11 + i + j
      ENDDO
   ENDDO
   icount = 0
   ilast = 55 + ilast
   DO i = N - Ic + 1 , N
      icount = icount + 1
      it = ilast + icount
      Ar(it) = 1.0D0 + dfloat(i)
      Ia(it) = i
      Ja(it) = i - N + Ic
   ENDDO
   ilast = ilast + Ic
   icount = 0
   DO i = N - Ic , N
      icount = icount + 1
      it = ilast + icount
      Ar(it) = -dfloat(i)
      Ia(it) = i
      Ja(it) = i - N + Ic + 1
   ENDDO
   ilast = ilast + Ic + 1
   icount = 0
   DO i = N - Ic - 1 , N
      icount = icount + 1
      it = ilast + icount
      Ar(it) = 16.0D0
      Ia(it) = i
      Ja(it) = i - N + Ic + 2
   ENDDO
!     ilast = ilast + ic + 2
!     if(ilast.ne.4*n+55) then
!     write(*,*)' ilast equal to ', ilast
!     write(*,*)' ILAST, the number of nonzeros, should = ', 4*n + 55
!     stop
!     end if
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
END SUBROUTINE dcn
!*==ecn.f90 processed by SPAG 8.04RA 15:47 21 Feb 2024
!!SPAG Open source Personal, Educational or Academic User Clemson University  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
SUBROUTINE ecn(N,Ic,Ne,Ia,Ja,Ar,Nn,Ierr)
   USE ISO_FORTRAN_ENV                 
   IMPLICIT NONE
!
! Dummy argument declarations rewritten by SPAG
!
   INTEGER , INTENT(IN) :: Nn
   INTEGER , INTENT(IN) :: N
   INTEGER , INTENT(IN) :: Ic
   INTEGER , INTENT(INOUT) :: Ne
   INTEGER , INTENT(OUT) , DIMENSION(Nn) :: Ia
   INTEGER , INTENT(OUT) , DIMENSION(Nn) :: Ja
   REAL(REAL64) , INTENT(OUT) , DIMENSION(Nn) :: Ar
   INTEGER , INTENT(OUT) :: Ierr
!
! Local variable declarations rewritten by SPAG
!
   INTEGER :: i , ilast , it
!
! End of declarations rewritten by SPAG
!
!----------------------------------------------------------------------
!
!   PURPOSE
!   -------
!   The subroutine generates sparse (square) matrices of the type
!   E(N,C).  This type of matrix has the following characteristics:
!   Symmetric, positive-definite, N x N matrices with 4 in the diagonal
!   and -1 in the two sidediagonal and in the two bands at the distance
!   C from the diagonal. These matrices are similar to matrices obtained
!   from using the five-point formula in the discretization of the
!   elliptic PDE.
!
!
!   Note: If A is the sparse matrix of type E(N,C), then
!
!       min|A(i,j)| = 1,     max|A(i,j)| = 4
!
!
!
!   CONTRIBUTOR: Ernest E. Rothman
!                Cornell Theory Center/Cornell National Supercomputer
!                Facility.
!                e-mail address: BITNET:   eer@cornellf
!                                INTERNET: eer@cornellf.tn.cornell.edu
!
!
!   REFERENCE
!   ---------
!   1) Zlatev, Zahari; Schaumburg, Kjeld; Wasniewski, Jerzy;
!      "A Testing Scheme for Subroutines Solving Large Linear Problems",
!       Computers and Chemistry, Vol. 5, No. 2-3, pp. 91-100, 1981.
!   2) Osterby, Ole and Zletev, Zahari;
!      "Direct Methods for Sparse Matrices";
!       Springer-Verlag 1983.
!
!
!
!   INPUT PARAMETERS
!   ----------------
!   N    - Integer. The size of the square matrix.
!          N > 2 must be specified.
!
!   NN   - Integer. The dimension of integer arrays IA and JA and
!          real array AR. Must be at least NE.
!
!   NN  - Integer. The dimension of integer array JA. Must be at least
!          NE.
!
!   IC   - Integer. The sparsity pattern can be changed by means of this
!          parameter.  1 < IC < N   must be specified.
!
!
!
!   OUTPUT PARAMETERS
!   -----------------
!   NE   - Integer. The number of nonzero elements in the sparse matrix
!          of the type E(N,C). NE = 5*N - 2*IC - 2 .
!
!   AR(NN)  - Real array.
!             Stored entries of the sparse matrix A.
!             NE is the number of nonzeros including a mandatory
!             diagonal entry for each row.
!
!   IA(NN)  - Integer array.(Double precision)
!             Pointers to specify rows for the stored nonzero entries
!             in AR.
!
!   JA(NN) - Integer array.
!             Pointers to specify columns for the stored nonzero entries
!             in AR.
!
!   IERR    - Error parameter is returned as zero on successful
!             execution of the subroutine.
!             Error diagnostics are given by means of positive values
!             of this parameter as follows:
!             IERR = 1    -  N       is out of range.
!             IERR = 2    -  IC      is out of range.
!             IERR = 3    -  NN      is out of range.
!
!---------------------------------------------------------------------
!
!
   Ierr = 0
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
!heck the input parameters:
!
   IF ( N<=2 ) THEN
      Ierr = 1
      RETURN
   ENDIF
   IF ( Ic<=1 .OR. Ic>=N ) THEN
      Ierr = 2
      RETURN
   ENDIF
!
   Ne = 5*N - 2*Ic - 2
   IF ( Nn<Ne ) THEN
      Ierr = 3
      RETURN
   ENDIF
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
! Begin to generate the nonzero elements as well as the row and column
! pointers:
!
   DO i = 1 , N
      Ar(i) = 4.0D0
      Ia(i) = i
      Ja(i) = i
   ENDDO
   ilast = N
   DO i = 1 , N - 1
      it = ilast + i
      Ar(it) = -1.0D0
      Ia(it) = i + 1
      Ja(it) = i
   ENDDO
   ilast = ilast + N - 1
   DO i = 1 , N - 1
      it = ilast + i
      Ar(it) = -1.0D0
      Ia(it) = i
      Ja(it) = i + 1
   ENDDO
   ilast = ilast + N - 1
   DO i = 1 , N - Ic
      it = ilast + i
      Ar(it) = -1.0D0
      Ia(it) = i + Ic
      Ja(it) = i
   ENDDO
   ilast = ilast + N - Ic
   DO i = 1 , N - Ic
      it = ilast + i
      Ar(it) = -1.0D0
      Ia(it) = i
      Ja(it) = i + Ic
   ENDDO
!      ilast = ilast + n-ic
!      if(ilast.ne.5*n-2*ic-2) then
!      write(*,*)' ilast equal to ', ilast
!      write(*,*)' ILAST, the no. of nonzeros, should = ', 5*n-2*ic-2
!      stop
!      end if
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
END SUBROUTINE ecn
 