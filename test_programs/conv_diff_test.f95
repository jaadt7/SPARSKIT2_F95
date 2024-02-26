!*==convdif.f90 processed by SPAG 8.04RA 10:54 23 Feb 2024
!!SPAG Open source Personal, Educational or Academic User Clemson University  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
PROGRAM convdif
!-----------------------------------------------------------------------
! this driver will generate a finite element matrix for the
! convection-diffusion problem
!
!                  - Div ( K(x,y) Grad u ) + C grad u = f
!                    u = 0 on boundary
!
! (Dirichlet boundary conditions).
!-----------------------------------------------------------------------
! this code will prompt for desired mesh (from 0 to 9, with one being
! a user input one) and will general the matrix in harewell-Boeing format
! in the file mat.hb. It also generates two post script files, one
! showing the pattern of the matrix (mat.ps) and the other showing the
! corresponding mesh (msh.ps).
!-----------------------------------------------------------------------
! the structure and organization of the fem codes follow very closely
! that of the book by Noborou Kikuchi (Finite element methods in
! mechanics,   Cambridge Univ. press, 1986).
!-----------------------------------------------------------------------
! coded Y. Saad and S. Ma -- this version dated August 11, 1993
!-----------------------------------------------------------------------
   USE ISO_FORTRAN_ENV                 
   IMPLICIT NONE
!
! PARAMETER definitions rewritten by SPAG
!
   INTEGER , PARAMETER :: MAXNX = 8000 , MAXNEL = 15000
!
! Local variable declarations rewritten by SPAG
!
   REAL(REAL64) , DIMENSION(7*MAXNX) :: a
   REAL(REAL64) , DIMENSION(3*MAXNX) :: f
   INTEGER , DIMENSION(MAXNX) :: ia , iperm , iwk , jwk , nodcode
   INTEGER , DIMENSION(12,MAXNEL) :: ichild
   INTEGER :: ierr , ifmt , ii , job , mode , n , nb , ncol , nelmax , nelxnew , nmesh , nref , nxmax , nxnew , ptitle
   INTEGER , SAVE :: iin , iout , na , ndeg , nelx , node , nx
   INTEGER , DIMENSION(3,MAXNEL) :: ijk
   INTEGER , DIMENSION(2,MAXNX) :: iparnts
   INTEGER , DIMENSION(7*MAXNX) :: ja
   CHARACTER(8) :: key
   CHARACTER(20) :: matfile
   CHARACTER(2) :: munt
   REAL :: size
   CHARACTER(72) :: title
   CHARACTER(3) :: type
   REAL(REAL64) , DIMENSION(MAXNX) :: x , y
   EXTERNAL assmbo , checkref , funb , func , fung , inmesh , prtmt , psgrid , pspltm , refall , xyk
!
! End of declarations rewritten by SPAG
!
!
!
   DATA iin/7/ , node/3/ , nx/0/ , nelx/0/ , iout/8/ , ndeg/12/ , na/3000/
!--------------------------------------------------------------
! choose starting mesh   ---
!--------------------------------------------------------------
!     files for output
!
   OPEN (UNIT=10,FILE='mat.hb')
   OPEN (UNIT=11,FILE='msh.ps')
   OPEN (UNIT=12,FILE='mat.ps')
!-----------------------------------------------------------------------
   PRINT * , ' enter chosen mesh '
   READ (*,*) nmesh
   IF ( nmesh==0 ) THEN
      PRINT * , 'enter input file for initial mesh '
      READ (*,'(a20)') matfile
      OPEN (UNIT=7,FILE=matfile)
   ENDIF
!-----------------------------------------------------------------------
   PRINT * , ' Enter the number of refinements desired '
   READ (*,*) nref
   CALL inmesh(nmesh,iin,nx,nelx,node,x,y,nodcode,ijk,iperm)
!
!     ...REFINE THE GRID
!
   nxmax = MAXNX
   nelmax = MAXNEL
   nb = 0
   DO ii = 1 , nref
!
!     estimate the number nx and nelx at next refinement level.
!
      CALL checkref(nx,nelx,nodcode,nb,nxnew,nelxnew)
      IF ( nxnew>nxmax .OR. nelxnew>nelmax ) THEN
         PRINT * , ' Was able to do only ' , ii - 1 , '  refinements'
         CALL spag_block_1
         STOP
      ENDIF
!
!     ...if OK refine  all elements
!
      CALL refall(nx,nelx,ijk,node,ndeg,x,y,ichild,iparnts,nodcode,nxmax,nelmax,ierr)
      IF ( ierr/=0 ) PRINT * , '** ERROR IN REFALL : ierr =' , ierr
   ENDDO
   CALL spag_block_1
CONTAINS
   SUBROUTINE spag_block_1
!-----------------------------------------------------------------------
      job = 0
!-----------------------------------------------------------------------
!     assemble the matrix in CSR format
!-----------------------------------------------------------------------
      CALL assmbo(nx,nelx,node,ijk,nodcode,x,y,a,ja,ia,f,iwk,jwk,ierr,xyk,funb,func,fung)
      n = nx
!----------------------Harewell-Boeing matrix---------------------------
!---------------------1---------2---------3---------5---------6
!            12345678901234567890123456789012345678901234567890
      title = '1Sample matrix from SPARSKIT  '
      key = 'SPARSKIT'
      type = 'rua'
      ifmt = 6
      job = 2
      CALL prtmt(n,n,a,ja,ia,f,'NN',title,key,type,ifmt,job,10)
!----------------------Plot of mesh-------------------------------------
!-----------------------------------------------------------------------
      size = 6.0
      munt = 'in'
      mode = 0
      title = 'Finite element mesh '
      ptitle = 1
      CALL psgrid(nx,ja,ia,x,y,title,ptitle,size,munt,11)
!      hsize  = 5.6
!      vsize = 5.6
!      xleft = 0.0
!      bot = 0.0
!      job = 30
!      call texgrd(nx,ja,ia,x,y,munt,size,vsize,hsize,
!     *     xleft,bot,job,title,ptitle,ijk,node,nelx,11)
!----------------------Plot of matrix-pattern---------------------------
!-----------------------------------------------------------------------
      size = 5.5
      mode = 0
      title = 'Assembled Matrix'
      ptitle = 1
      ncol = 0
      iout = 12
      CALL pspltm(n,n,mode,ja,ia,title,ptitle,size,munt,ncol,iwk,12)
   END SUBROUTINE spag_block_1
!      xleft = 0.00
!      bot = 0.70
!      job =  0
!      call texplt(nx,nx,mode,ja,ia,munt,size,vsize,hsize,xleft,bot,
!     *     job,title,ptitle,ncol,iwk,12)
!-----end-of-program-convdif--------------------------------------------
!-----------------------------------------------------------------------
END PROGRAM convdif
