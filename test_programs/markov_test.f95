!*==markov.f90 processed by SPAG 8.04RA 11:58 23 Feb 2024
!!SPAG Open source Personal, Educational or Academic User Clemson University  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
PROGRAM markov
   USE ISO_FORTRAN_ENV                 
   IMPLICIT NONE
!
! PARAMETER definitions rewritten by SPAG
!
   INTEGER , PARAMETER :: NMAX = 5000 , NZMAX = 4*NMAX
!
! Local variable declarations rewritten by SPAG
!
   REAL(REAL64) , DIMENSION(NZMAX) :: a
   INTEGER , DIMENSION(NMAX+1) :: ia
   INTEGER :: ifmt , iout , job , m , n
   INTEGER , DIMENSION(NZMAX) :: ja
   CHARACTER(8) :: key
   CHARACTER(72) :: title
   CHARACTER(3) :: type
   REAL :: x
   EXTERNAL markgen , prtmt
!
! End of declarations rewritten by SPAG
!
!-----------------------------------------------------------------------
!
! program to generate a Markov chain matrix (to test eigenvalue routines
! or algorithms for singular systems (in which case use I-A ))
! the matrix models  simple random walk on a triangular grid.
! see additional comments in subroutine.
! -----
! just compile this segment and link to the rest of sparskit
! (uses subroutine prtmt from MATGEN)
! will create a matrix in the HARWELL/BOEING format and put it in
! the file markov.mat
!
!-----------------------------------------------------------------------
!
   OPEN (UNIT=11,FILE='markov.mat')
!
! read - in grid size - will not accept too large grids.
!
   WRITE (6,'(17hEnter grid-size: ,$)')
   READ * , m
   IF ( m*(m+1)>2*NMAX ) THEN
      PRINT * , ' m too large - unable to produce matrix '
      STOP
   ENDIF
!
! call generator.
!
   CALL markgen(m,n,a,ja,ia)
!-----------------------------------------------------------------------
   title = ' Test matrix from SPARSKIT - markov chain model           '
   key = 'randwk01'
   type = 'rua'
   iout = 11
   job = 2
   ifmt = 10
   CALL prtmt(n,n,a,ja,ia,x,'NN',title,key,type,ifmt,job,iout)
END PROGRAM markov
!*==markgen.f90 processed by SPAG 8.04RA 11:58 23 Feb 2024
!!SPAG Open source Personal, Educational or Academic User Clemson University  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
!
SUBROUTINE markgen(M,N,A,Ja,Ia)
   USE ISO_FORTRAN_ENV                 
   IMPLICIT NONE
!
! Dummy argument declarations rewritten by SPAG
!
   INTEGER , INTENT(IN) :: M
   INTEGER , INTENT(OUT) :: N
   REAL(REAL64) , INTENT(INOUT) , DIMENSION(*) :: A
   INTEGER , INTENT(OUT) , DIMENSION(*) :: Ja
   INTEGER , INTENT(OUT) , DIMENSION(*) :: Ia
!
! Local variable declarations rewritten by SPAG
!
   REAL(REAL64) :: cst , pd , pu
   REAL(REAL64) , SAVE :: half
   INTEGER :: i , ix , j , jax , jmax
!
! End of declarations rewritten by SPAG
!
!-----------------------------------------------------------------------
! matrix generator for a markov model of a random walk on a triang. grid
!-----------------------------------------------------------------------
! this subroutine generates a test matrix that models a random
! walk on a triangular grid. This test example was used by
! G. W. Stewart ["{SRRIT} - a FORTRAN subroutine to calculate the
! dominant invariant subspaces of a real matrix",
! Tech. report. TR-514, University of Maryland (1978).] and in a few
! papers on eigenvalue problems by Y. Saad [see e.g. LAA, vol. 34,
! pp. 269-295 (1980) ]. These matrices provide reasonably easy
! test problems for eigenvalue algorithms. The transpose of the
! matrix  is stochastic and so it is known that one is an exact
! eigenvalue. One seeks the eigenvector of the transpose associated
! with the eigenvalue unity. The problem is to calculate the
! steady state probability distribution of the system, which is
! the eigevector associated with the eigenvalue one and scaled in
! such a way that the sum all the components is equal to one.
!-----------------------------------------------------------------------
! parameters
!------------
! on entry :
!----------
! m     = integer. number of points in each direction.
!
! on return:
!----------
! n     = integer. The dimension of the matrix. (In fact n is known
!         to be equal to (m(m+1))/2      )
! a,
! ja,
! ia    = the matrix stored in CSR format.
!
!-----------------------------------------------------------------------
! Notes: 1) the code will actually compute the transpose of the
! stochastic matrix that contains the transition probibilities.
!        2) It should also be possible to have a matrix generator
! with an additional parameter (basically redefining `half' below
! to be another parameter and changing the rest accordingly, but
! this is not as simple as it sounds). This is not likely to provide
! any more interesting matrices.
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
   DATA half/0.5D0/
!
   cst = half/real(M-1)
!
!     --- ix counts the grid point (natural ordering used), i.e.,
!     --- the row number of the matrix.
!
   ix = 0
   jax = 1
   Ia(1) = jax
!
!     sweep y coordinates
!
   DO i = 1 , M
      jmax = M - i + 1
!
!     sweep x coordinates
!
      DO j = 1 , jmax
         ix = ix + 1
         IF ( j/=jmax ) THEN
            pd = cst*real(i+j-1)
!
!     north
!
            A(jax) = pd
            IF ( i==1 ) A(jax) = A(jax) + pd
            Ja(jax) = ix + 1
            jax = jax + 1
!     east
            A(jax) = pd
            IF ( j==1 ) A(jax) = A(jax) + pd
            Ja(jax) = ix + jmax
            jax = jax + 1
         ENDIF
!     south
         pu = half - cst*real(i+j-3)
         IF ( j>1 ) THEN
            A(jax) = pu
            Ja(jax) = ix - 1
            jax = jax + 1
         ENDIF
!     west
         IF ( i>1 ) THEN
            A(jax) = pu
            Ja(jax) = ix - jmax - 1
            jax = jax + 1
         ENDIF
         Ia(ix+1) = jax
      ENDDO
   ENDDO
   N = ix
END SUBROUTINE markgen
