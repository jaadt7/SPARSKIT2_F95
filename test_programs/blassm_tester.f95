module parameters
        integer, parameter :: nxmax = 30, nmx = 900, nzmax = 6300
end module parameters

program matprod
! -----------------------------------------------------------------------
!       test program for some routines in BLASSM.f95
! ----------------------------------------------------------------------- 
!       Last update: May 2, 1994
! ----------------------------------------------------------------------- 
        use parameters
! BEGIN new declarations
        integer :: nx,ny,nz,n,ierr,k,j,ifmt,job,jj
        real :: s
! 
        integer, dimension(nmx + 1) :: ia, ib
        integer, dimension(nzmax) :: ic, ja, jb, jc
        integer, dimension(nmx) :: iw
!       
        real(kind=8), dimension(nzmax) :: a, b, c
        real(kind=8), dimension(nmx) :: x, y, y1, rhs
        real(kind=8), dimension(6) :: al
!END new declarations 

!initializing
        nx = 20
        ny = 20
        nz = 1
        al(1) = 1.0d0
        al(2) = 0.0d0
        al(3) = 2.3d1
        al(4) = 0.4d0
        al(5) = 0.0d0
        al(6) = 8.2d-2
       
! -----------------------------------------------------------------------
        call gen57pt (nx,ny,nz,al,0,n,a,ja,ia,iw,rhs)
        call gen57pt (ny,nx,nz,al,0,n,b,jb,ib,iw,rhs)
! 
        s = 3.812
!         
        call aplsb(n,n,a,ja,ia,s,b,jb,ib,c,jc,ic,nzmax,iw,ierr)
        if (ierr  /=  0) print *,' ierr = ',ierr
! 
 
        do k=1,n
           x(k) = real(k)/real(n)
        end do
! 
        call ope(n,x,y1,a,ja,ia) 
        call ope(n,x,y,b,jb,ib)

        y1 = s*y+y1 
!      
        call ope(n,x,y,c,jc,ic) 
! ------------------------------------------------------
        print *, ' ------------ checking APLSB --------------'
        call ydfnorm(n,y1,y,6) 
! ------------------------------------------------------
        ifmt = 103
!      
        job = -1
! --------             
        do jj=1,2
          print *, 'DUMP A____________________________'
          call dump (1,n,.true.,a,ja,ia,9)
          print *, 'DUMP B____________________________'
          call dump (1,n,.true.,b,jb,ib,9) 
          call apmbt(n,n,job,a,ja,ia,b,jb,ib,c,jc,ic,nzmax,iw,ierr)
          print *, 'DUMP C____________________________'
          call dump (1,n,.true.,c,jc,ic,9) 
          if (ierr  /=  0) print *,' ierr = ',ierr
          call ope(n,x,y1,a,ja,ia) 
          call opet(n,x,y,b,jb,ib) 
          s = real(job) 
          
          y1 = y1 + s*y
          
!      
          call ope(n,x,y,c,jc,ic) 
! ------------xs------------------------------------------
          print *, '  '
          print *, ' ------------ checking APMBT---------------'
          print *, ' ------------ with JOB = ',job,' -------------'
          call ydfnorm(n,y1,y,6) 
! ------------------------------------------------------
          job = job + 2 
        end do
! 
        s = 0.1232445
        call aplsbt(n,n,a,ja,ia,s,b,jb,ib,c,jc,ic,nzmax,iw,ierr)
! 
        if (ierr  /=  0) print *,' ierr = ',ierr
        call ope(n,x,y1,a,ja,ia) 
        call opet(n,x,y,b,jb,ib) 
        
        y1 = y1 + s*y
! 
        call ope(n,x,y,c,jc,ic) 
! ------------------------------------------------------
! ------------------------------------------------------
        print *,  '  '
        print *, ' ------------ checking APLSBT---------------'
        call ydfnorm(n,y1,y,6) 
! ----------------------------------------------------------------------- 
!  testing products
! -----------------------------------------------------------------------
        job = 1 
        call amub (n,n,job,a,ja,ia,b,jb,ib,c,jc,ic,nzmax,iw,ierr)
! 
        if (ierr  /=  0) print *,' ierr = ',ierr
        call ope(n,x,y,b,jb,ib) 
        call ope(n,y,y1,a,ja,ia) 
! 
        call ope(n,x,y,c,jc,icr) 
! ----------------------------------------------------------------------- 
        print *, '  '
        print *, ' ------------ checking AMUB  ---------------'
 
        call ydfnorm(n,y1,y) 
! 

end program matprod

subroutine ope(n,x,y,a,ja,ia)
        use parameters
! BEGIN declarations
        integer, intent(In) :: n
        integer :: i,k1,k2,k
        real(kind = 8), dimension(nmx), intent(In) :: x
        real(kind = 8), dimension(nzmax),intent(In) :: a
        integer, dimension(nmx+1), intent(In) :: ia
        integer, dimension(nzmax), intent(in) :: ja
        real(kind = 8), dimension(nmx), intent(out) :: y
! END declarations

!  sparse matrix * vector multiplication
!         
        do i=1,n
            k1 = ia(i)
            k2 = ia(i+1) -1
            y(i) = 0.0
           
            do k=k1,k2
                y(i) = y(i) + a(k)*x(ja(k)) 
            end do
        end do
 
        return
end subroutine ope

subroutine opet(n,x,y,a,ja,ia)
        use parameters

! BEGIN declarations
        integer, intent(In) :: n
        integer :: j,i,k
        real(kind=8), dimension(nmx), intent(In) :: x
        real(kind=8), dimension(nmx), intent(Out) :: y
        real(kind=8), dimension(nzmax), intent(In) :: a
        integer, dimension(nmx + 1), intent(In) :: ia
        integer, dimension(nzmax), intent(in) :: ja

! END declarations
 
!  sparse matrix * vector multiplication
!         
        y = 0.d0
! 
        do i=1,n
           do k=ia(i), ia(i+1)-1 
                y(ja(k)) = y(ja(k)) + x(i)*a(k)
           end do
        end do
        
        return
end subroutine opet
 
subroutine ydfnorm(n,y1,y)
        use parameters

! BEGIN declarations
        integer, intent(In) :: n, iout
        real :: t
        integer :: k
        real(kind=8), dimension(nmx), intent(In) :: y, y1
! END declarations
! 
        t = 0.0d0
        do k=1,n
           t = t+(y(k)-y1(k))**2
        end do
        t = sqrt(t) 
        print *, '2-norm of error (exact answer-tested answer)=',t
! -----------------------------------------------------------------------
 return
 end subroutine ydfnorm
