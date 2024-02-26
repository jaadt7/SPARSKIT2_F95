!*==rmatvec.f90 processed by SPAG 8.04RA 18:18 22 Feb 2024
!!SPAG Open source Personal, Educational or Academic User Clemson University  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
PROGRAM rmatvec
   USE ISO_FORTRAN_ENV                 
   IMPLICIT NONE
!
! PARAMETER definitions rewritten by SPAG
!
   INTEGER , PARAMETER :: NMAX = 10000 , NZMAX = 80000
!
! Local variable declarations rewritten by SPAG
!
   REAL(REAL64) , DIMENSION(NZMAX) :: a1 , a2
   INTEGER , DIMENSION(NMAX) :: ia1 , ia2 , iwk1 , iwk2
   INTEGER :: idiag , ierr , ii , j , jdiag , jj , k , n , ndiag , nfree , nlev , nx , ny , nz
   INTEGER , DIMENSION(11) , SAVE :: idim
   INTEGER , DIMENSION(10) :: ioff
   INTEGER , SAVE :: iout
   INTEGER , DIMENSION(NZMAX) :: ja1 , ja2 , jad
   REAL(REAL64) :: scal
   REAL(REAL64) , DIMENSION(100) :: stencil
   REAL(REAL64) , DIMENSION(NMAX) :: x , y , y0 , y1
   EXTERNAL amux , amuxd , amuxe , amuxj , atmux , csrcsc , csrdia , csrell , csrjad , csrmsr , dvperm , errpr , gen57bl , getl ,  &
          & getu , jadcsr , ldsol , ldsolc , ldsoll , levels , lsol , udsol , udsolc , usol
!
! End of declarations rewritten by SPAG
!
!-----------------------------------------------------------------------
! This test program tests all the subroutines in matvec.
! it generates matrices and transforms them in appropriate formats
! and then call the appropriate routines.
!-----------------------------------------------------------------------
! common used only to generate nonsymmetric matrices
!      common /gam/ gamma, gamma1, cvar
   DATA idim/4 , 10 , 15 , 40 , 50 , 60 , 70 , 80 , 90 , 100 , 200/
   DATA iout/6/
!
!     initialize common gam
!
!      gamma = 0.5
!      gamma1 = 1.0
!      cvar = 1.0
!-----------------------------------------------------------------------
!     ii loop corresponds to size of problem
!-----------------------------------------------------------------------
   DO ii = 1 , 3
      WRITE (iout,*) '---------------- ii ' , ii , '--------------------'
      nfree = 1
      nx = idim(ii)
      ny = nx
!-----------------------------------------------------------------------
!     jj loop corresponds to 2-D and 3-D problems.
!-----------------------------------------------------------------------
      DO jj = 1 , 2
         WRITE (iout,*) '     ----------- jj ' , jj , ' -------------'
         nz = 1
         IF ( jj==2 ) nz = 10
!
!     call matrix generation routine --
!     (strange to use block version to generate 1 x 1 blocks...)
!
         CALL gen57bl(nx,ny,nz,1,1,n,a1,ja1,ia1,ia2,stencil)
!
!     initialize x
!
         DO j = 1 , n
            x(j) = real(j)
         ENDDO
!
! initial call to get `` exact '' answer in y0
!
         CALL amux(n,x,y0,a1,ja1,ia1)
!-----------------------------------------------------------------------
! TESTING AMUXE
!-----------------------------------------------------------------------
!
!     convert to itpack format -----
!
         CALL csrell(n,a1,ja1,ia1,7,a2,jad,n,ndiag,ierr)
         CALL amuxe(n,x,y,n,ndiag,a2,jad)
         CALL errpr(n,y,y0,iout,'amuxe ')
!-----------------------------------------------------------------------
! TESTING AMUXD
!-----------------------------------------------------------------------
!
!     convert to diagonal format
!
         idiag = 7
         CALL csrdia(n,idiag,10,a1,ja1,ia1,NMAX,a2,ioff,a2,ja2,ia2,jad)
         CALL amuxd(n,x,y,a2,NMAX,idiag,ioff)
         CALL errpr(n,y,y0,iout,'amuxd ')
!-----------------------------------------------------------------------
! TESTING ATMUX
!-----------------------------------------------------------------------
!
!    convert to csc format (transpose)
!
         CALL csrcsc(n,1,1,a1,ja1,ia1,a2,ja2,ia2)
         CALL atmux(n,x,y,a2,ja2,ia2)
         CALL errpr(n,y,y0,iout,'atmux ')
!-----------------------------------------------------------------------
! TESTING AMUXJ
!-----------------------------------------------------------------------
!
!     convert to jagged diagonal format
!
         CALL csrjad(n,a1,ja1,ia1,jdiag,jad,a2,ja2,ia2)
         CALL amuxj(n,x,y,jdiag,a2,ja2,ia2)
         CALL dvperm(n,y,jad)
         CALL errpr(n,y,y0,iout,'amuxj ')
!
! convert back
 
         CALL jadcsr(n,jdiag,a2,ja2,ia2,jad,a1,ja1,ia1)
         CALL amux(n,x,y,a1,ja1,ia1)
         CALL errpr(n,y,y0,iout,'jadcsr')
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!       triangular systems solutions
!-----------------------------------------------------------------------
! TESTING LDSOL
!-----------------------------------------------------------------------
         CALL getl(n,a1,ja1,ia1,a2,ja2,ia2)
         CALL amux(n,x,y0,a2,ja2,ia2)
         CALL atmux(n,x,y1,a2,ja2,ia2)
         CALL csrmsr(n,a2,ja2,ia2,a2,ja2,y,iwk2)
         DO k = 1 , n
            a2(k) = 1.0D0/a2(k)
         ENDDO
         CALL ldsol(n,y,y0,a2,ja2)
         CALL errpr(n,x,y,iout,'ldsol ')
!-----------------------------------------------------------------------
! TESTING LDSOLL
!-----------------------------------------------------------------------
         CALL levels(n,ja2,ja2,nlev,jad,iwk1,iwk2)
         CALL ldsoll(n,y,y0,a2,ja2,nlev,jad,iwk1)
         CALL errpr(n,x,y,iout,'ldsoll')
!-----------------------------------------------------------------------
! TESTING UDSOLC
!-----------------------------------------------------------------------
! here we take advantage of the fact that the MSR format for U
! is the MSC format for L
!
         CALL udsolc(n,y,y1,a2,ja2)
         CALL errpr(n,x,y,iout,'udsolc')
!-----------------------------------------------------------------------
! TESTING LSOL
!-----------------------------------------------------------------------
! here we exploit the fact that with MSR format a, ja, ja is actually
! the correct data structure for the strict lower triangular part of
! the CSR format. First rescale matrix.
!
         scal = 0.1
         DO k = ja2(1) , ja2(n+1) - 1
            a2(k) = a2(k)*scal
         ENDDO
         CALL amux(n,x,y0,a2,ja2,ja2)
         DO j = 1 , n
            y0(j) = x(j) + y0(j)
         ENDDO
         CALL lsol(n,y,y0,a2,ja2,ja2)
         CALL errpr(n,x,y,iout,'lsol  ')
!-----------------------------------------------------------------------
! TESTING UDSOL
!-----------------------------------------------------------------------
         CALL getu(n,a1,ja1,ia1,a2,ja2,ia2)
         CALL amux(n,x,y0,a2,ja2,ia2)
         CALL atmux(n,x,y1,a2,ja2,ia2)
         CALL csrmsr(n,a2,ja2,ia2,a2,ja2,y,jad)
         DO k = 1 , n
            a2(k) = 1.0D0/a2(k)
         ENDDO
         CALL udsol(n,y,y0,a2,ja2)
         CALL errpr(n,x,y,iout,'udsol ')
!-----------------------------------------------------------------------
! TESTING LDSOLC
!-----------------------------------------------------------------------
! here we take advantage of the fact that the MSR format for L
! is the MSC format for U
!
         CALL ldsolc(n,y,y1,a2,ja2)
         CALL errpr(n,x,y,iout,'ldsolc')
!-----------------------------------------------------------------------
! TESTING USOL
!-----------------------------------------------------------------------
! here we exploit the fact that with MSR format a, ja, ja is actually
! the correct data structure for the strict lower triangular part of
! the CSR format. First rescale matrix.
!
         scal = 0.1
         DO k = ja2(1) , ja2(n+1) - 1
            a2(k) = a2(k)*scal
         ENDDO
         CALL amux(n,x,y1,a2,ja2,ja2)
         DO j = 1 , n
            y1(j) = x(j) + y1(j)
         ENDDO
         CALL usol(n,y,y1,a2,ja2,ja2)
         CALL errpr(n,x,y,iout,'usol  ')
!-----------------------------------------------------------------------
!   --  END --
!-----------------------------------------------------------------------
      ENDDO
   ENDDO
!---------------end-of-main---------------------------------------------
!-----------------------------------------------------------------------
END PROGRAM rmatvec
!*==errpr.f90 processed by SPAG 8.04RA 18:18 22 Feb 2024
!!SPAG Open source Personal, Educational or Academic User Clemson University  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
!-----------------------------------------------------------------------
SUBROUTINE errpr(N,Y,Y1,Iout,Msg)
   USE ISO_FORTRAN_ENV                 
   IMPLICIT NONE
!
! Dummy argument declarations rewritten by SPAG
!
   INTEGER , INTENT(IN) :: N
   REAL(REAL64) , INTENT(IN) , DIMENSION(*) :: Y
   REAL(REAL64) , INTENT(IN) , DIMENSION(*) :: Y1
   INTEGER , INTENT(IN) :: Iout
   CHARACTER(6) , INTENT(IN) :: Msg
!
! Local variable declarations rewritten by SPAG
!
   INTEGER :: k
   REAL(REAL64) :: t
!
! End of declarations rewritten by SPAG
!
   t = 0.0D0
   DO k = 1 , N
      t = t + (Y(k)-Y1(k))**2
   ENDDO
   t = sqrt(t)
   WRITE (Iout,*) ' 2-norm of difference in ' , Msg , ' =' , t
END SUBROUTINE errpr
 
