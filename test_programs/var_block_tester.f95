!*==rvbr.f90 processed by SPAG 8.04RA 23:02 22 Feb 2024
!!SPAG Open source Personal, Educational or Academic User Clemson University  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
PROGRAM rvbr
   USE ISO_FORTRAN_ENV                 
   IMPLICIT NONE
!
! PARAMETER definitions rewritten by SPAG
!
   INTEGER , PARAMETER :: NXMAX = 10 , NMX = NXMAX*NXMAX , NNZMAX = 10*NMX
!
! Local variable declarations rewritten by SPAG
!
   REAL(REAL64) , DIMENSION(NNZMAX) :: a , a1 , b
   REAL(REAL64) , DIMENSION(NMX) :: ans , rhs , x
   INTEGER :: i , ierr , job , maxblock , n , na , nc , nfree , nr , nx , ny , nz
   INTEGER , DIMENSION(NMX+1) :: ia , ib , kvstc , kvstr
   INTEGER , DIMENSION(NNZMAX) :: ia1 , ja , ja1
   INTEGER , DIMENSION(NMX*2+1) :: iwk
   INTEGER , DIMENSION(NMX*10) :: jb , kb
   REAL(REAL64) , EXTERNAL :: rnd
   REAL(REAL64) , DIMENSION(7,100) :: stencil
   EXTERNAL amux , bsrcsr , csrvbr , gen57bl , vbrcsr , vbrinfo , vbrmv
!
! End of declarations rewritten by SPAG
!
!-----------------------------------------------------------------------
!     SPARSKIT test program for Variable Block Matrix Support
!-----------------------------------------------------------------------
!     This program tests all three conversion routines of csrvbr.
!     For each conversion to VBR, the format is converted back to CSR
!     with vbrcsr.  The subroutines csrkvstr, csrkvstc, and kvstmerge
!     become tested in the process.  The subroutines vbrinfo and vbrmv
!     are also tested.
!-----------------------------------------------------------------------
!-----dimension of grid
   nx = 4
   ny = 2
   nz = 1
   nfree = 2
!-----generate grid problem.
   na = nfree*nfree
   CALL gen57bl(nx,ny,nz,nfree,na,n,a1,ja1,ia1,iwk,stencil)
!-----convert matrix to CSR
   CALL bsrcsr(1,n,nfree,na,a1,ja1,ia1,a,ja,ia)
   n = n*nfree
!     call dump(1, n, .true., a, ja, ia, 6)
!-----generate random x vector for testing matrix-vector product
   DO i = 1 , n
      x(i) = rnd()
   ENDDO
!-----generate correct solution for matrix-vector product
   CALL amux(n,x,ans,a,ja,ia)
   DO job = 0 , 2
      print *, 'Testing job = ' , job
      IF ( job==0 ) THEN
!-----------maximum blocksize for random block partitioning
         maxblock = n/4
!-----------generate random block partitioning for rows
         nr = 1
         kvstr(1) = 1
         SPAG_Loop_2_1: DO
            nr = nr + 1
            kvstr(nr) = kvstr(nr-1) + int(rnd()*maxblock) + 1
            IF ( kvstr(nr)>=n+1 ) THEN
               kvstr(nr) = n + 1
               nr = nr - 1
!-----------generate random block partitioning for columns
               nc = 1
               kvstc(1) = 1
               DO
                  nc = nc + 1
                  kvstc(nc) = kvstc(nc-1) + int(rnd()*maxblock) + 1
                  IF ( kvstc(nc)>=n+1 ) THEN
                     kvstc(nc) = n + 1
                     nc = nc - 1
                     EXIT SPAG_Loop_2_1
                  ENDIF
               ENDDO
            ENDIF
         ENDDO SPAG_Loop_2_1
      ENDIF
!--------convert to VBR format------------------------------------------
      CALL csrvbr(n,ia,ja,a,nr,nc,kvstr,kvstc,ib,jb,kb,b,job,iwk,NMX*10,NNZMAX,ierr)
!--------convert back to CSR format-------------------------------------
      CALL vbrcsr(ia1,ja1,a1,nr,kvstr,kvstc,ib,jb,kb,b,NNZMAX,ierr)
!--------compare original and converted CSR structures if job not 0
      print *, 'Checking conversions....'
      IF ( job/=0 ) THEN
         DO i = 1 , n
            IF ( ia(i)/=ia1(i) ) THEN
               print *, 'csrvbr or vbrcsr conversion mismatch'
               STOP
            ENDIF
         ENDDO
         DO i = 1 , ia(n+1) - 1
            IF ( (ja(i)/=ja1(i)) .OR. (a(i)/=a1(i)) ) THEN
               print *, 'csrvbr or vbrcsr conversion mismatch'
               STOP
            ENDIF
         ENDDO
      ENDIF
!--------test vbrinfo---------------------------------------------------
      CALL vbrinfo(nr,nc,kvstr,kvstc,ib,jb,kb,iwk,6)
!--------test vbrmv-----------------------------------------------------
      CALL vbrmv(nr,nc,ib,jb,b,kvstr,kvstc,x,rhs)
!--------compare answer with answer computed with CSR format
      DO i = 1 , n
         IF ( abs(ans(i)-rhs(i))>abs(0.001D0*ans(i)) ) THEN
            print *, 'VBR matrix-vector product is erroneous ' , i
            STOP
         ENDIF
      ENDDO
!--------fill CSR structure with garbage
      DO i = 1 , ia1(n+1) - 1
         ja1(i) = -1
         a1(i) = -1.D0
      ENDDO
      DO i = 1 , n + 1
         ia1(i) = -1
      ENDDO
!--------fill VBR structure with garbage
      DO i = 1 , kb(ib(nr+1)) - 1
         b(i) = -1.D0
      ENDDO
      DO i = 1 , ib(nr+1)
         jb(i) = -1
         kb(i) = -1
      ENDDO
      DO i = 1 , nr + 1
         ib(i) = -1
      ENDDO
!--------fill kvstr and kvstc with garbage
      DO i = 1 , nr + 1
         kvstr(i) = -1
      ENDDO
      DO i = 1 , nc + 1
         kvstc(i) = -1
      ENDDO
!--------fill rhs with garbage
      DO i = 1 , n
         rhs(i) = -1.D0
      ENDDO
!-----endloop on job
   ENDDO
END PROGRAM rvbr
!*==rnd.f90 processed by SPAG 8.04RA 23:02 22 Feb 2024
!!SPAG Open source Personal, Educational or Academic User Clemson University  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
!-----------------------------------------------------------------------
FUNCTION rnd()
   USE ISO_FORTRAN_ENV                 
   IMPLICIT NONE
   REAL(REAL64) :: rnd
!
! Local variable declarations rewritten by SPAG
!
   INTEGER , SAVE :: ia , ic , im , jran
!
! End of declarations rewritten by SPAG
!
   DATA im/6075/ , ia/106/ , ic/1283/ , jran/1/
   jran = mod(jran*ia+ic,im)
   rnd = dble(jran)/dble(im)
END FUNCTION rnd
!-----------------------------------------------------------------------
!     Coded by Edmond Chow, chow@cs.umn.edu
!-----------------------------------------------------------------------
