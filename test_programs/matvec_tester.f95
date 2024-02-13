module parameters
    integer, parameter :: nmax=10000, nzmax = 80000
end module

program rmatvec

    use parameters

! -----------------------------------------------------------------------
!  This test program tests all the subroutines in matvec.
!  it generates matrices and transforms them in appropriate formats
!  and then call the appropriate routines. 
! -----------------------------------------------------------------------
! BEGIN declarations
    integer :: ii, nfree, nx, ny, jj, nz, n, j, ndiag, ierr, idiag, jdiag, k, nlev
    real :: scal
    integer, dimension(nmax) :: ia1, ia2, ja1, ja2, iwk1, iwk2
    integer, dimension(1:nzmax) :: ja1, ja2, jad
    integer, dimension(11) :: idim = (/ 4,10,15,40,50,60,70,80,90,100,200 / )
    integer, dimension(10) :: ioff
    real(kind=8), dimension(nzmax) :: a1, a2
    real(kind=8), dimension(nmax) :: x, y, y0, y1
    real(kind=8), dimension(100) :: stencil
! END declarations

! -----------------------------------------------------------------------
!      ii loop corresponds to size of problem 
! -----------------------------------------------------------------------
    do ii = 1, 3 
        print *, '---------------- ii ',ii,'--------------------' 
        nfree = 1
        nx = idim(ii) 
        ny = nx
! -----------------------------------------------------------------------
!      jj loop corresponds to 2-D and 3-D problems.
! -----------------------------------------------------------------------
        do jj=1, 2
            print *, '     ----------- jj ',jj,' -------------' 
            nz = 1
            if (jj == 2) nz = 10
!      
!      call matrix generation routine --
!      (strange to use block version to generate 1 x 1 blocks...)
!      
            call gen57bl (nx,ny,nz,1,1,n,a1,ja1,ia1,ia2,stencil)
!      
!      initialize x
!      
            do j=1, n
                x(j) = real(j)
            end do
!      
!  initial call to get `` exact '' answer in y0
! 
            call amux(n,x,y0, a1, ja1, ia1) 
! -----------------------------------------------------------------------
!  TESTING AMUXE
! -----------------------------------------------------------------------
!      
!      convert to itpack format -----
!      
            call csrell (n,a1,ja1,ia1,7,a2,jad,n,ndiag,ierr)
            call amuxe (n, x, y, n, ndiag, a2,jad)
            call errpr(n,y,y0,'amuxe ')

! -----------------------------------------------------------------------
!  TESTING AMUXD
! -----------------------------------------------------------------------
!      
!      convert to diagonal format
! 
            idiag = 7
            call csrdia (n, idiag,10,a1, ja1, ia1, nmax, a2, ioff, a2, ja2, ia2, jad) 
            call amuxd (n,x,y,a2,nmax,idiag,ioff) 
            call errpr(n,y,y0,'amuxd ')
! -----------------------------------------------------------------------
!  TESTING ATMUX
! -----------------------------------------------------------------------
! 
!     convert to csc format (transpose)
!     
            call csrcsc (n,1,1,a1,ja1,ia1,a2,ja2,ia2)
            call atmux (n, x, y, a2, ja2, ia2)
            call errpr(n,y,y0,'atmux ')
! -----------------------------------------------------------------------
!  TESTING AMUXJ
! -----------------------------------------------------------------------
!      
!      convert to jagged diagonal format
!      
            call csrjad (n,a1,ja1,ia1, jdiag, jad, a2, ja2, ia2) 
            call amuxj (n, x, y, jdiag, a2, ja2, ia2) 
            call dvperm (n, y, jad) 
            call errpr(n,y,y0,'amuxj ')
! 
!  convert back
            call jadcsr (n, jdiag, a2, ja2, ia2, jad, a1, ja1, ia1)
            call amux (n, x, y, a1, ja1, ia1)
            call errpr(n,y, y0,'jadcsr')
! -----------------------------------------------------------------------
! -----------------------------------------------------------------------
!        triangular systems solutions
! -----------------------------------------------------------------------
!  TESTING LDSOL
! -----------------------------------------------------------------------
            call getl (n, a1, ja1, ia1, a2, ja2, ia2)
            call amux (n,x,y0, a2, ja2, ia2)
            call atmux(n,x,y1, a2, ja2, ia2)
            call csrmsr (n, a2, ja2, ia2, a2, ja2, y, iwk2)
            
            do k = 1,n
                a2(k) = 1.d0/a2(k)
            end do
            
            call ldsol (n, y, y0, a2, ja2)
            call errpr(n,x,y,'ldsol ')
! -----------------------------------------------------------------------
!  TESTING LDSOLL
! ----------------------------------------------------------------------- 
            call levels (n, ja2, ja2, nlev, jad, iwk1, iwk2) 
            call ldsoll (n, y, y0, a2, ja2, nlev, jad, iwk1) 
            call errpr(n,x,y,'ldsoll')
! -----------------------------------------------------------------------
!  TESTING UDSOLC
! -----------------------------------------------------------------------
!  here we take advantage of the fact that the MSR format for U
!  is the MSC format for L
!  
            call udsolc (n, y, y1, a2, ja2)
            call errpr(n,x,y,'udsolc')
! -----------------------------------------------------------------------
!  TESTING LSOL 
! -----------------------------------------------------------------------
!  here we exploit the fact that with MSR format a, ja, ja is actually
!  the correct data structure for the strict lower triangular part of
!  the CSR format. First rescale matrix. 
! 
            scal = 0.1
            do k=ja2(1), ja2(n+1)-1
                a2(k)=a2(k)*scal
            end do
            
            call amux(n, x, y0, a2, ja2, ja2) 
            
            do j=1,n
                y0(j) = x(j) + y0(j) 
            end do

            call lsol (n, y, y0, a2, ja2, ja2) 
            call errpr(n,x,y,'lsol  ')
! -----------------------------------------------------------------------
!  TESTING UDSOL
! -----------------------------------------------------------------------
            call getu (n, a1, ja1, ia1, a2, ja2, ia2)
            call amux (n,x,y0, a2, ja2, ia2)
            call atmux(n,x,y1, a2, ja2, ia2)
            call csrmsr (n, a2, ja2, ia2, a2, ja2, y, jad) 
 
            do k=1,n
                a2(k) = 1.0d0/ a2(k)
            end do
            
            call udsol (n, y, y0, a2, ja2)
            call errpr(n,x,y,'udsol ')

! -----------------------------------------------------------------------
!  TESTING LDSOLC
! -----------------------------------------------------------------------
!  here we take advantage of the fact that the MSR format for L
!  is the MSC format for U 
!  
            call ldsolc (n, y, y1, a2, ja2)
            call errpr(n,x,y,'ldsolc')

! -----------------------------------------------------------------------
!  TESTING USOL 
! -----------------------------------------------------------------------
!  here we exploit the fact that with MSR format a, ja, ja is actually
!  the correct data structure for the strict lower triangular part of
!  the CSR format. First rescale matrix. 
! 
            scal = 0.1
            
            do k=ja2(1), ja2(n+1)-1
                a2(k)=a2(k)*scal
            end do
            
            call amux(n, x, y1, a2, ja2, ja2) 

            do j=1,n
                y1(j) = x(j) + y1(j) 
            end do
 
            call usol (n, y, y1, a2, ja2, ja2) 
            call errpr(n,x,y,'usol  ')
! -----------------------------------------------------------------------
!    --  END --
! ----------------------------------------------------------------------- 
        end do
    end do
! ---------------end-of-main---------------------------------------------
! -----------------------------------------------------------------------
end program rmatvec

subroutine errpr(n,y,y1,msg)
    use parameters
! BEGIN new declarations
 integer, intent(In) :: n
 integer :: k
 real(kind=8), dimension(nmax), intent(In) :: y, y1
 real(kind=8) :: t, sqrt
 character(len=6), intent(In) :: msg
! END new declarations
    t = 0.0d0
    do k=1,n
        t = t+(y(k)-y1(k))**2
    end do
 
    t = sqrt(t)
        
    print *, ' 2-norm of difference in ',msg,' =', t
    return

end subroutine errpr
