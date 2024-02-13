subroutine amub(nrow,ncol,job,a,ja,ia,b,jb,ib,c,jc,ic,nzmax,iw,ierr)
! BEGIN declarations
        integer, intent(In) :: nrow, nrow, ncol, job, nzmax
        integer, intent(Out) :: ierr
        integer :: len,ii,ka,jj,kb,jcol,jpos,k
        real(kind=8), dimension(:), intent(In) :: a,b
        real(kind=8), dimension(:), intent(Out) :: c
        integer, dimension(:), intent(In) :: ja, jb
        integer, dimension(:), intent(Out) :: jc
        integer, dimension(nrow + 1), intent(In) :: ia
        integer, dimension(:), intent(In) :: ib
        integer, dimension(:), intent(Out) :: ic
        integer, dimension(ncol), intent(inOut) :: iw

! END declarations
! -----------------------------------------------------------------------
!  performs the matrix by matrix product C = A B 
! -----------------------------------------------------------------------
!  on entry:
!  ---------
!  nrow  = integer. The row dimension of A = row dimension of C
!  ncol  = integer. The column dimension of B = column dimension of C
!  job   = integer. Job indicator. When job = 0, only the structure
!                   (i.e. the arrays jc, ic) is computed and the
!                   real values are ignored
! 
!  a,
!  ja,
!  ia   = Matrix A in compressed sparse row format.
!  
!  b, 
!  jb, 
!  ib    =  Matrix B in compressed sparse row format.
! 
!  nzmax = integer. The  length of the arrays c and jc.
!          amub will stop if the result matrix C  has a number 
!          of elements that exceeds exceeds nzmax. See ierr.
!  
!  on return:
! ----------
!  c, 
!  jc, 
!  ic    = resulting matrix C in compressed sparse row sparse format.
!            
!  ierr  = integer. serving as error message. 
!          ierr = 0 means normal return,
!          ierr .gt. 0 means that amub stopped while computing the
!          i-th row  of C with i=ierr, because the number 
!          of elements in C exceeds nzmax.
! 
!  work arrays:
! ------------
!  iw    = integer work array of length equal to the number of
!          columns in A.
!  Note: 
! -------
!    The row dimension of B is not needed. However there is no checking 
!    on the condition that ncol(A) = nrow(B). 
! 
! ----------------------------------------------------------------------- 
        real(kind=8) :: scal
        logical :: values
        values = (job \= 0) 
        len = 0
        ic(1) = 1 
        ierr = 0
!       initialize array iw.
        iw = 0
! 
        do ii=1, nrow 
!      row i 
           do ka=ia(ii), ia(ii+1)-1 
                if (values) scal = a(ka)
                jj = ja(ka)
                do kb=ib(jj),ib(jj+1)-1
                     jcol = jb(kb)
                     jpos = iw(jcol)
                     if (jpos == 0) then
                        len = len + 1
                        if (len > nzmax) then
                                ierr = ii
                                return
                        endif
                        
                        jc(len) = jcol
                        iw(jcol)= len
                        if (values) c(len) = scal*b(kb)
                     else
                        if (values) c(jpos) = c(jpos) + scal*b(kb)
                     endif
                end do
           end do
           do k=ic(ii), len
               iw(jc(k)) = 0
           end do
           ic(ii+1) = len+1
        end do
        return
! -------------end-of-amub-----------------------------------------------
! -----------------------------------------------------------------------
end subroutine amub

subroutine aplb(nrow,ncol,job,a,ja,ia,b,jb,ib,c,jc,ic,nzmax,iw,ierr)
! BEGIN declarations

        integer, intent(In) :: nrow,ncol,job,nzmax
        integer, intent(Out) :: ierr
        integer :: len,ii,ka,jcol,kb,jpos,k
        real(kind=8), dimension(:), intent(In) :: a, b
        real(kind=8), dimension(:), intent(Out) :: c
        integer, dimension(:), intent(In) :: ja, jb
        integer, dimension(:), intent(Out) :: jc
        integer, dimension(nrow + 1), intent(In) :: ia, ib
        integer, dimension(nrow + 1), intent(Out) :: ic
        integer, dimension(ncol), intent(InOut) :: iw
! END  declarations

! -----------------------------------------------------------------------
! -----------------------------------------------------------------------
!  performs the matrix sum  C = A+B. 
! -----------------------------------------------------------------------
!  on entry:
!  ---------
!  nrow     = integer. The row dimension of A and B
!  ncol  = integer. The column dimension of A and B.
!  job   = integer. Job indicator. When job = 0, only the structure
!                   (i.e. the arrays jc, ic) is computed and the
!                   real values are ignored.
! 
!  a,
!  ja,
!  ia   = Matrix A in compressed sparse row format.
!  
!  b, 
!  jb, 
!  ib     =  Matrix B in compressed sparse row format.
! 
!  nzmax	= integer. The  length of the arrays c and jc.
!          amub will stop if the result matrix C  has a number 
!          of elements that exceeds exceeds nzmax. See ierr.
!  
!  on return:
! ----------
!  c, 
!  jc, 
!  ic	= resulting matrix C in compressed sparse row sparse format.
!          
!  ierr	= integer. serving as error message. 
!          ierr = 0 means normal return,
!          ierr .gt. 0 means that amub stopped while computing the
!          i-th row  of C with i=ierr, because the number 
!          of elements in C exceeds nzmax.
! 
!  work arrays:
! ------------
!  iw	= integer work array of length equal to the number of
!          columns in A.
! ----------------------------------------------------------------------
        logical :: values
        values = (job \= 0) 
        ierr = 0
        len = 0
        ic(1) = 1 
        iw = 0
!      
        do ii=1, nrow
!      row i 
             do ka=ia(ii), ia(ii+1)-1 
                  len = len+1
                  jcol = ja(ka)
                  if (len > nzmax) then
                          ierr = ii
                          return
                  end if
                  jc(len) = jcol 
                  if (values) c(len) = a(ka) 
                  iw(jcol)= len
             end do
!      
             do kb=ib(ii),ib(ii+1)-1
                  jcol = jb(kb)
                  jpos = iw(jcol)
                  if (jpos == 0) then
                           len = len+1
                           if (len > nzmax) then
                                   ierr = ii
                                   return
                           end if
                           jc(len) = jcol
                           if (values) c(len) = b(kb)
                           iw(jcol)= len
                  else
                           if (values) c(jpos) = c(jpos) + b(kb)
                  endif
             end do
 
             do k=ic(ii), len
                 iw(jc(k)) = 0
             end do
             ic(ii+1) = len+1
        end do

return
! ------------end of aplb ----------------------------------------------- 
! -----------------------------------------------------------------------
end subroutine aplb
 
subroutine aplb1(nrow,ncol,job,a,ja,ia,b,jb,ib,c,jc,ic,nzmax,ierr)
! BEGIN declarations
        integer, intent(In) :: nrow, ncol, job, nzmax
        integer, intent(Out) :: ierr
        integer :: kc, i, ka,kb,kamax,kbmax,j1,j2
        real(kind=8), dimension(:), intent(In) :: a,b
        real(kind=8), dimension(:), intent(Out) :: c
        integer, dimension(:), intent(In) :: ja, jb
        integer, dimension(:), intent(Out) :: jc
        integer, dimension(nrow + 1), intent(In) :: ia, ib
        integer, dimension(nrow + 1), intent(Out) :: ic
! END declarations
! -----------------------------------------------------------------------
!  performs the matrix sum  C = A+B for matrices in sorted CSR format.
!  the difference with aplb  is that the resulting matrix is such that
!  the elements of each row are sorted with increasing column indices in
!  each row, provided the original matrices are sorted in the same way. 
! -----------------------------------------------------------------------
!  on entry:
!  ---------
!  nrow	= integer. The row dimension of A and B
!  ncol  = integer. The column dimension of A and B.
!  job   = integer. Job indicator. When job = 0, only the structure
!                   (i.e. the arrays jc, ic) is computed and the
!                   real values are ignored.
! 
!  a,
!  ja,
!  ia   = Matrix A in compressed sparse row format with entries sorted
!  
!  b, 
!  jb, 
!  ib	=  Matrix B in compressed sparse row format with entries sorted
!         ascendly in each row   
! 
!  nzmax	= integer. The  length of the arrays c and jc.
!          amub will stop if the result matrix C  has a number 
!          of elements that exceeds exceeds nzmax. See ierr.
!  
!  on return:
! ----------
!  c, 
!  jc, 
!  ic	= resulting matrix C in compressed sparse row sparse format
!          with entries sorted ascendly in each row. 
!          
!  ierr	= integer. serving as error message. 
!          ierr = 0 means normal return,
!          ierr .gt. 0 means that amub stopped while computing the
!          i-th row  of C with i=ierr, because the number 
!          of elements in C exceeds nzmax.
! 
!  Notes: 
! -------
!      this will not work if any of the two input matrices is not sorted
! -----------------------------------------------------------------------
        logical :: values
        values = (job \= 0) 
        ierr = 0
        kc = 1
        ic(1) = kc 
! 
        do i=1, nrow
        ka = ia(i)
        kb = ib(i)
        kamax = ia(i+1)-1
        kbmax = ib(i+1)-1
        do while (ka <= kamax .or. kb <= kb)
                 
                if (ka <= kamax) then
                   j1 = ja(ka)
                else
                   j1 = ncol+1
                endif
 
                if (kb <= kbmax) then 
                        j2 = jb(kb) 
                else 
                        j2 = ncol+1
                endif

                if (kc >= nzmax) then
                        ierr = i
                        return
                end if

! 
!      three cases
!      
                select case (.true.)
                        case (j1 == j2)
                                if (values) c(kc) = a(ka) + b(kb)
                                jc(kc) = j1
                                ka = ka + 1
                                kb = kb + 1 
                                kc = kc + 1
                        case (j1 < j2)
                                jc(kc) = j1
                                if (values) c(kc) = a (ka)
                                ka = ka + 1
                                kc = kc + 1
                        case (j1 > j2)
                                jc(kc) = j2
                                if (values) c(kc) = b(kb)
                                kb = kb + 1
                                kc = kc + 1
                end select

        end do
 
        return
! ------------end-of-aplb1----------------------------------------------- 
! -----------------------------------------------------------------------
end subroutine aplb1

subroutine aplsb(nrow,ncol,a,ja,ia,s,b,jb,ib,c,jc,ic,nzmax,iw,ierr)
! BEGIN new declarations

        integer, intent(In) :: nrow, ncol, nzmax
        integer, intent(Out) :: ierr
        integer :: len, ii, ka, jcol, kb, jpos, k
        real(kind=8), dimension(:), intent(In) :: a, b
        real(kind=8), dimension(:), intent(Out) :: c
        real(kind=8), intent(In) :: s
        integer, dimension(:), intent(In) :: ja,jb
        integer, dimension(*), intent(Out) :: jc
        integer, dimension(nrow + 1), intent(In) :: ia, ib
        integer, dimension(nrow + 1), intent(Out) :: ic
        integer, dimension(ncol), intent(inout) :: iw
! END new declarations
! -----------------------------------------------------------------------
! -----------------------------------------------------------------------
!  performs the matrix sum  C = A+s*B. 
! -----------------------------------------------------------------------
!  on entry:
!  ---------
!  nrow     = integer. The row dimension of A and B
!  ncol  = integer. The column dimension of A and B.
!  job   = integer. Job indicator. When job = 0, only the structure
!                   (i.e. the arrays jc, ic) is computed and the
!                   real values are ignored.
! 
!  a,
!  ja,
!  ia   = Matrix A in compressed sparse row format.
!  
!  s    = real*8 - coefficient that multiplies B.
!  b, 
!  jb, 
!  ib     =  Matrix B in compressed sparse row format.
! 
!  nzmax	= integer. The  length of the arrays c and jc.
!          amub will stop if the result matrix C  has a number 
!          of elements that exceeds exceeds nzmax. See ierr.
!  
!  on return:
! ----------
!  c, 
!  jc, 
!  ic	= resulting matrix C in compressed sparse row sparse format.
!          
!  ierr	= integer. serving as error message. 
!          ierr = 0 means normal return,
!          ierr .gt. 0 means that aplsb1 stopped while computing the
!          i-th row  of C with i=ierr, because the number 
!          of elements in C exceeds nzmax.
! 
!  work arrays:
! ------------
!  iw	= integer work array of length equal to the number of
!          columns in A.
!  note: expanded  row implementation. Does not require column indices to
!        be sorted.
! ----------------------------------------------------------------------
        
        ierr = 0
        len = 0
        ic(1) = 1
        iw = 0
!      
        do ii=1, nrow
!      copy row ii  to C
             do ka=ia(ii), ia(ii+1)-1 
                len = len+1
                jcol = ja(ka)
                if (len > nzmax) then
                        ierr = ii
                        return
                end if

                jc(len) = jcol 
                c(len) = a(ka) 
                iw(jcol)= len
             end do
!      
             do kb=ib(ii),ib(ii+1)-1
                jcol = jb(kb)
                jpos = iw(jcol)
                if (jpos == 0) then
                   len = len+1
                   if (len > nzmax) then
                          ierr = ii
                          return
                   end if

                    jc(len) = jcol
                    c(len) = s*b(kb)
                    iw(jcol)= len
                else
                    c(jpos) = c(jpos) + s*b(kb)
                endif

             end do
 
             do k=ic(ii), len
                iw(jc(k)) = 0
             end do
             ic(ii+1) = len+1
        end do
 
        return
! ------------end of aplsb ---------------------------------------------- 
! -----------------------------------------------------------------------
end subroutine aplsb

subroutine aplsb1(nrow,ncol,a,ja,ia,s,b,jb,ib,c,jc,ic,nzmax,ierr)
! BEGIN new declarations

        integer, intent(In) :: nrow, ncol, nzmax
        integer, intent(Out) :: ierr
        integer :: kc, i, ka, kb, kamax, kbmax, j1, j1
        real(kind=8), dimension(:), intent(In) :: a, b
        real(kind=8), dimension(:), intent(Out) :: c
        real(kind=8), intent(In) :: s
        integer, dimension(:), intent(In) :: ja, jb
        integer, dimension(:), intent(Out) :: jc
        integer, dimension(nrow + 1), intent(In) :: ia, ib
        integer, dimension(nrow + 1), intent(Out) :: ic
! END new declarations
! -----------------------------------------------------------------------
!  performs the operation C = A+s B for matrices in sorted CSR format.
!  the difference with aplsb is that the resulting matrix is such that
!  the elements of each row are sorted with increasing column indices in
!  each row, provided the original matrices are sorted in the same way. 
! -----------------------------------------------------------------------
!  on entry:
!  ---------
!  nrow	= integer. The row dimension of A and B
!  ncol  = integer. The column dimension of A and B.
! 
!  a,
!  ja,
!  ia   = Matrix A in compressed sparse row format with entries sorted
! 
!  s	= real. scalar factor for B.
!  
!  b, 
!  jb, 
!  ib	=  Matrix B in compressed sparse row format with entries sorted
!         ascendly in each row   
! 
!  nzmax	= integer. The  length of the arrays c and jc.
!          amub will stop if the result matrix C  has a number 
!          of elements that exceeds exceeds nzmax. See ierr.
!  
!  on return:
! ----------
!  c, 
!  jc, 
!  ic	= resulting matrix C in compressed sparse row sparse format
!          with entries sorted ascendly in each row. 
!          
!  ierr	= integer. serving as error message. 
!          ierr = 0 means normal return,
!          ierr .gt. 0 means that amub stopped while computing the
!          i-th row  of C with i=ierr, because the number 
!          of elements in C exceeds nzmax.
! 
!  Notes: 
! -------
!      this will not work if any of the two input matrices is not sorted
! -----------------------------------------------------------------------
        ierr = 0
        kc = 1
        ic(1) = kc 
! 
!      the following loop does a merge of two sparse rows + adds  them.
!  
        do i=1, nrow
             ka = ia(i)
             kb = ib(i)
             kamax = ia(i+1)-1
             kbmax = ib(i+1)-1 
             
             do while (ka <= kamax .or. kb <= kbmax)

!      
                if (ka <= kamax) then
                    j1 = ja(ka)
                else
!      take j1 large enough  that always j2 .lt. j1
                    j1 = ncol+1
                endif
                if (kb <= kbmax) then 
                   j2 = jb(kb) 
                else 
!      similarly take j2 large enough  that always j1 .lt. j2 
                    j2 = ncol+1
                endif
!      
!      three cases
!      
                select case(.true.)
                case (j1 == j2)
                        c(kc)  = a(ka) + s*b(kb)
                        jc(kc) = j1
                        
                        ka = ka + 1
                        kb = kb + 1
                        kc = kc + 1
                case (j1 < j2)
                        jc(kc) = j1
                        c(kc)  =  a(ka)

                        ka = ka + 1
                        kc = kc + 1
                case (j1 > j2)
                        jc(kc) = j2
                        c(kc)  = s*b(kb)
                        
                        kb = kb + 1
                        kc = kc + 1
                end select  

                if (kc > nzmax) then
                        ierr = i
                        return
                end if

             end do 
! 
!      end while loop
! 
 
             ic(i+1) = kc
        end do
        return
! ------------end-of-aplsb1 --------------------------------------------- 
! -----------------------------------------------------------------------
end subroutine aplsb1

subroutine apmbt(nrow,ncol,job,a,ja,ia,b,jb,ib,c,jc,ic,nzmax,iw,ierr)
! BEGIN declarations
 
        integer, intent(In) :: nrow, ncol, job,nzmax
        integer, intent(Out) :: ierr
        integer :: nnza, nnzb, len, ljob,ipos, k, ii, ka, jcol, jpos, i
        real(kind=8), dimension(:), intent(In) :: a,b
        real(kind=8), dimension(:), intent(InOut) :: c
        integer, dimension(:), intent(In) :: ja,jb
        integer, dimension(:), intent(InOut) :: jc
        integer, dimension(nrow + 1), intent(In) :: ia, ib
        integer, dimension(:), intent(InOut) :: ic
        integer, dimension(:), intent(InOut) :: iw
! END declarations
! -----------------------------------------------------------------------
!  performs the matrix sum  C = A + transp(B) or C = A - transp(B) 
! -----------------------------------------------------------------------
!  on entry:
!  ---------
!  nrow	= integer. The row dimension of A and transp(B)
!  ncol  = integer. The column dimension of A. Also the row 
!                   dimension of B. 
! 
!  job	= integer. if job = -1, apmbt will compute C= A - transp(B)
!          (structure + values) 
!          if (job .eq. 1)  it will compute C=A+transp(A) 
!          (structure+ values) 
!          if (job .eq. 0) it will compute the structure of
!          C= A+/-transp(B) only (ignoring all real values).
!          any other value of job will be treated as  job=1
!  a,
!  ja,
!  ia    = Matrix A in compressed sparse row format.
! 
!  b, 
!  jb, 
!  ib	=  Matrix B in compressed sparse row format.
! 
!  nzmax	= integer. The  length of the arrays c, jc, and ic.
!          amub will stop if the result matrix C  has a number 
!          of elements that exceeds exceeds nzmax. See ierr.
!  
!  on return:
! ----------
!  c, 
!  jc, 
!  ic	= resulting matrix C in compressed sparse row format.
!          
!  ierr	= integer. serving as error message. 
!          ierr = 0 means normal return.
!          ierr = -1 means that nzmax was .lt. either the number of
!          nonzero elements of A or the number of nonzero elements in B.
!          ierr .gt. 0 means that amub stopped while computing the
!          i-th row  of C with i=ierr, because the number 
!          of elements in C exceeds nzmax.
! 
!  work arrays:
! ------------
!  iw	= integer work array of length at least max(ncol,nrow) 
! 
!  Notes:
! ------- It is important to note that here all of three arrays c, ic, 
!         and jc are assumed to be of length nnz(c). This is because 
!         the matrix is internally converted in coordinate format.
!         
! -----------------------------------------------------------------------
        logical :: values
        values = (job \= 0) 
! 
        ierr = 0
        iw = 0
!      
        nnza = ia(nrow+1)-1
        nnzb = ib(ncol+1)-1
        len = nnzb

        if (nzmax < nnzb .or. nzmax < nnza) then
                ierr = -1
                return
        endif
!      
!  trasnpose matrix b into c
! 
        ljob = 0
        if (values) ljob = 1
        ipos = 1
        call csrcsc (ncol,ljob,ipos,b,jb,ib,c,jc,ic) 
! ----------------------------------------------------------------------- 
        if (job == -1) then
                do k=1,len
                        c(k) = -c(k)
                end do
        endif
! 
! --------------- main loop --------------------------------------------
! 
        do ii=1, nrow
                do k = ic(ii),ic(ii+1)-1
                        iw(jc(k)) = k
                end do
! -----------------------------------------------------------------------     
                do ka = ia(ii), ia(ii+1)-1 
                        jcol = ja(ka)
                        jpos = iw(jcol)
                        if (jpos == 0) then
! 
!      if fill-in append in coordinate format to matrix.
!  
                                len = len+1
                                if (len > nzmax) then
                                        ierr = ii
                                        return
                                end if
                        
                                jc(len) = jcol
                                ic(len) = ii
                                if (values) c(len) = a(ka)
                        else
!      else do addition.
                                if (values) c(jpos) = c(jpos) + a(ka)
                        endif
                end do
                
                do k=ic(ii), ic(ii+1)-1
                        iw(jc(k)) = 0
                end do
        end do
!      
!      convert first part of matrix (without fill-ins) into coo format
!      
        ljob = 2
        if (values) ljob = 3
        do i=1, nrow+1
                iw(i) = ic(i) 
        end do
        call csrcoo (nrow,ljob,nnzb,c,jc,iw,nnzb,c,ic,jc,ierr)
! 
!      convert the whole thing back to csr format. 
!  
        ljob = 0
        if (values) ljob = 1
        call coicsr (nrow,len,ljob,c,jc,ic,iw)
 
        return
! --------end-of-apmbt---------------------------------------------------
! -----------------------------------------------------------------------
end subroutine apmbt

subroutine aplsbt(nrow,ncol,a,ja,ia,s,b,jb,ib,c,jc,ic,nzmax,iw,ierr)
! BEGIN new declarations
        integer, intent(In) :: nrow, ncol, nzmax
        integer, intent(Out) :: ierr
        integer :: j, nnza, nnzb, len, ljob, ipos, k, ii, ka, jcol, jpos, i
        real(kind=8), dimension(:), intent(In) :: a, b
        real(kind=8), dimension(:), intent(InOut) :: c
        real(kind=8), intent(In) :: s
        integer, dimension(:), intent(In) :: ja, jb
        integer, dimension(:), intent(InOut) :: jc
        integer, dimension(nrow + 1), intent(In) :: ia, ib
        integer, dimension(:), intent(Out) :: ic
        integer, dimension(:), intent(inout) :: iw
! END new declarations
! -----------------------------------------------------------------------
!  performs the matrix sum  C = A + transp(B).
! -----------------------------------------------------------------------
!  on entry:
!  ---------
!  nrow	= integer. The row dimension of A and transp(B)
!  ncol  = integer. The column dimension of A. Also the row 
!                   dimension of B. 
! 
!  a,
!  ja,
!  ia    = Matrix A in compressed sparse row format.
! 
!  s	= real. scalar factor for B.
! 
!  
!  b, 
!  jb, 
!  ib	=  Matrix B in compressed sparse row format.
! 
!  nzmax	= integer. The  length of the arrays c, jc, and ic.
!          amub will stop if the result matrix C  has a number 
!          of elements that exceeds exceeds nzmax. See ierr.
!  
!  on return:
! ----------
!  c, 
!  jc, 
!  ic	= resulting matrix C in compressed sparse row format.
!          
!  ierr	= integer. serving as error message. 
!          ierr = 0 means normal return.
!          ierr = -1 means that nzmax was .lt. either the number of
!          nonzero elements of A or the number of nonzero elements in B.
!          ierr .gt. 0 means that amub stopped while computing the
!          i-th row  of C with i=ierr, because the number 
!          of elements in C exceeds nzmax.
! 
!  work arrays:
! ------------
!  iw	= integer work array of length at least max(nrow,ncol) 
! 
!  Notes:
! ------- It is important to note that here all of three arrays c, ic, 
!         and jc are assumed to be of length nnz(c). This is because 
!         the matrix is internally converted in coordinate format.
!         
! -----------------------------------------------------------------------
        ierr = 0
        do j=1, ncol
                iw(j) = 0
        end do
!      
        nnza = ia(nrow+1)-1
        nnzb = ib(ncol+1)-1
        len = nnzb
        if (nzmax < nnzb .or. nzmax < nnza) then
                ierr = -1
                return
        endif
!      
!      transpose matrix b into c
! 
        ljob = 1
        ipos = 1
        call csrcsc (ncol,ljob,ipos,b,jb,ib,c,jc,ic) 
        do k=1,len
                c(k) = c(k)*s
        end do
!      
!      main loop. add rows from ii = 1 to nrow.
!      
        do ii=1, nrow
!      iw is used as a system to recognize whether there
!      was a nonzero element in c. 
                do k = ic(ii),ic(ii+1)-1
                        iw(jc(k)) = k
                end do
!      
                do ka = ia(ii), ia(ii+1)-1 
                        jcol = ja(ka)
                        jpos = iw(jcol)
                        if (jpos == 0) then
!      
!      if fill-in append in coordinate format to matrix.
!      
                                len = len+1
                                if (len > nzmax) then
                                        ierr = ii
                                        return
                                end if
                                jc(len) = jcol 
                                ic(len) = ii
                                c(len) = a(ka)
                        else
!      else do addition.
                                c(jpos) = c(jpos) + a(ka)
                        endif
                end do
 
                do k=ic(ii), ic(ii+1)-1
                        iw(jc(k)) = 0
                end do
        end do
!      
!      convert first part of matrix (without fill-ins) into coo format
!      
        ljob = 3
        do i=1, nrow+1
                iw(i) = ic(i) 
        end do
        call csrcoo (nrow,ljob,nnzb,c,jc,iw,nnzb,c,ic,jc,ierr)
! 
!      convert the whole thing back to csr format. 
!  
        ljob = 1
        call coicsr (nrow,len,ljob,c,jc,ic,iw)
        return
! --------end-of-aplsbt--------------------------------------------------
! -----------------------------------------------------------------------
end subroutine aplsbt

subroutine diamua(nrow,job,a,ja,ia,diag,b,jb,ib)
! BEGIN new declarations
        integer, intent(In) :: nrow, job
        integer :: ii, k1, k2, k
        real(kind=8), dimension(:), intent(In) :: a
        real(kind=8), dimension(:), intent(Out) :: b
        real(kind=8), dimension(nrow), intent(In) :: diag
        real(kind=8) :: scal
        integer, dimension(*), intent(In) :: ja
        integer, dimension(*), intent(Out) :: jb
        integer, dimension(nrow + 1), intent(In) :: ia
        integer, dimension(nrow + 1), intent(Out) :: ib
! END new declarations
! -----------------------------------------------------------------------
!  performs the matrix by matrix product B = Diag * A  (in place) 
! -----------------------------------------------------------------------
!  on entry:
!  ---------
!  nrow	= integer. The row dimension of A
! 
!  job   = integer. job indicator. Job=0 means get array b only
!          job = 1 means get b, and the integer arrays ib, jb.
! 
!  a,
!  ja,
!  ia   = Matrix A in compressed sparse row format.
!  
!  diag = diagonal matrix stored as a vector dig(1:n)
! 
!  on return:
! ----------
! 
!  b, 
!  jb, 
!  ib	= resulting matrix B in compressed sparse row sparse format.
!          
!  Notes:
! -------
!  1)        The column dimension of A is not needed. 
!  2)        algorithm in place (B can take the place of A).
!            in this case use job=0.
! -----------------------------------------------------------------
        do ii=1,nrow
!      
!      normalize each row 
!      
                k1 = ia(ii)
                k2 = ia(ii+1)-1
                scal = diag(ii) 
                
                do k=k1, k2
                        b(k) = a(k)*scal
                end do
        end do
!      
        if (job == 0) return
!      
        do ii=1, nrow+1
                ib(ii) = ia(ii)
        end do
 
        do k=ia(1), ia(nrow+1) -1 
                jb(k) = ja(k)
        end do
        return
! ----------end-of-diamua------------------------------------------------
! -----------------------------------------------------------------------
end subroutine diamua

subroutine amudia(nrow,job,a,ja,ia,diag,b,jb,ib)
! BEGIN new declarations
        integer, intent(In) :: nrow, job
        integer :: ii, k1, k2, k
        real(kind=8), dimension(*), intent(In) :: a
        real(kind=8), dimension(*), intent(Out) :: b
        real(kind=8), dimension(nrow), intent(In) :: diag
        integer, dimension(:), intent(In) :: ja
        integer, dimension(:), intent(Out) :: jb
        integer, dimension(nrow + 1), intent(In) :: ia
        integer, dimension(nrow + 1), intent(Out) :: ib
! END new declarations
! -----------------------------------------------------------------------
!  performs the matrix by matrix product B = A * Diag  (in place) 
! -----------------------------------------------------------------------
!  on entry:
!  ---------
!  nrow	= integer. The row dimension of A
! 
!  job   = integer. job indicator. Job=0 means get array b only
!          job = 1 means get b, and the integer arrays ib, jb.
! 
!  a,
!  ja,
!  ia   = Matrix A in compressed sparse row format.
!  
!  diag = diagonal matrix stored as a vector dig(1:n)
! 
!  on return:
! ----------
! 
!  b, 
!  jb, 
!  ib	= resulting matrix B in compressed sparse row sparse format.
!          
!  Notes:
! -------
!  1)        The column dimension of A is not needed. 
!  2)        algorithm in place (B can take the place of A).
! -----------------------------------------------------------------
        do ii=1,nrow
!      
!      scale each element 
!      
                k1 = ia(ii)
                k2 = ia(ii+1)-1
 
                do k=k1, k2
                        b(k) = a(k)*diag(ja(k)) 
                end do
        end do
!      
        if (job == 0) return
!      
        do ii=1, nrow+1
                ib(ii) = ia(ii)
        end do
 
        do k=ia(1), ia(nrow+1) -1 
                jb(k) = ja(k)
        end do
        return
! -----------------------------------------------------------------------
! -----------end-of-amudiag----------------------------------------------
end subroutine amudia

subroutine aplsca(nrow,a,ja,ia,scal,iw)
! BEGIN new declarations
        integer, intent(In) :: nrow
        integer :: icount, j, ko, ii, k1, k2, k
        real(kind=8), dimension(:), intent(InOut) :: a
        real(kind=8), intent(In) :: scal
        integer, dimension(:), intent(InOut) :: ja
        integer, dimension(nrow + 1), intent(InOut) :: ia
        integer, dimension(:), intent(InOut) :: iw
! END new declarations
! -----------------------------------------------------------------------
!  Adds a scalar to the diagonal entries of a sparse matrix A :=A + s I 
! -----------------------------------------------------------------------
!  on entry:
!  ---------
!  nrow	= integer. The row dimension of A
! 
!  a,
!  ja,
!  ia    = Matrix A in compressed sparse row format.
!  
!  scal  = real. scalar to add to the diagonal entries. 
! 
!  on return:
! ----------
! 
!  a, 
!  ja, 
!  ia	= matrix A with diagonal elements shifted (or created).
!          
!  iw    = integer work array of length n. On return iw will
!          contain  the positions of the diagonal entries in the 
!          output matrix. (i.e., a(iw(k)), ja(iw(k)), k=1,...n,
!          are the values/column indices of the diagonal elements 
!          of the output matrix. ). 
! 
!  Notes:
! -------
!      The column dimension of A is not needed. 
!      important: the matrix a may be expanded slightly to allow for
!      additions of nonzero elements to previously nonexisting diagonals.
!      The is no checking as to whether there is enough space appended
!      to the arrays a and ja. if not sure allow for n additional 
!      elemnts. 
!      coded by Y. Saad. Latest version July, 19, 1990
! -----------------------------------------------------------------------
        logical :: test
! 
        call diapos (nrow,ja,ia,iw)
        icount = 0
        do j=1, nrow
                if (iw(j) == 0) then
                        icount = icount+1
                else
                        a(iw(j)) = a(iw(j)) + scal 
                endif
        end do
! 
!      if no diagonal elements to insert in data structure return.
! 
        if (icount == 0) return
! 
!  shift the nonzero elements if needed, to allow for created 
!  diagonal elements. 
! 
        ko = ia(nrow+1)+icount
! 
!      copy rows backward
! 
        do ii=nrow, 1, -1 
!      
!      go through  row ii
!      
                k1 = ia(ii)
                k2 = ia(ii+1)-1 
                ia(ii+1) = ko
                test = (iw(ii) .eq. 0) 
                do k = k2,k1,-1 
                        j = ja(k)
                        if (test .and. (j < ii)) then 
                                test = .false. 
                                ko = ko - 1
                                a(ko) = scal 
                                ja(ko) = ii
                                iw(ii) = ko
                        endif
                        ko = ko-1
                        a(ko) = a(k) 
                        ja(ko) = j
                end do
!      diagonal element has not been added yet.
                if (test) then
                        ko = ko-1
                        a(ko) = scal 
                        ja(ko) = ii
                        iw(ii) = ko
                endif
        end do
        ia(1) = ko 
        return
! -----------------------------------------------------------------------
! ----------end-of-aplsca------------------------------------------------ 
end subroutine aplsca

subroutine apldia(nrow,job,a,ja,ia,diag,b,jb,ib,iw)
! BEGIN new declarations
        integer, intent(In) :: nrow, job
        integer :: nnz, k ,icount, j, ko, ii, k1, k2
        real(kind=8), dimension(:), intent(In) :: a
        real(kind=8), dimension(:), intent(Out) :: b
        real(kind=8), dimension(nrow), intent(In) :: diag
        integer, dimension(:), intent(In) :: ja
        integer, dimension(:), intent(Out) :: jb
        integer, dimension(nrow + 1), intent(In) :: ia
        integer, dimension(nrow + 1), intent(Out) :: ib
        integer, dimension(:), intent(InOut) :: iw
! END new declarations
! -----------------------------------------------------------------------
!  Adds a diagonal matrix to a general sparse matrix:  B = A + Diag 
! -----------------------------------------------------------------------
!  on entry:
!  ---------
!  nrow	= integer. The row dimension of A
! 
!  job   = integer. job indicator. Job=0 means get array b only
!          (i.e. assume that a has already been copied into array b,
!          or that algorithm is used in place. ) For all practical 
!          purposes enter job=0 for an in-place call and job=1 otherwise
!  
!          Note: in case there are missing diagonal elements in A, 
!          then the option job =0 will be ignored, since the algorithm 
!          must modify the data structure (i.e. jb, ib) in this 
!          situation.
! 
!  a,
!  ja,
!  ia   = Matrix A in compressed sparse row format.
!      
!  diag = diagonal matrix stored as a vector dig(1:n)
! 
!  on return:
! ----------
! 
!  b, 
!  jb, 
!  ib	= resulting matrix B in compressed sparse row sparse format.
! 
! 
!  iw    = integer work array of length n. On return iw will
!          contain  the positions of the diagonal entries in the 
!          output matrix. (i.e., a(iw(k)), ja(iw(k)), k=1,...n,
!          are the values/column indices of the diagonal elements 
!          of the output matrix. ). 
! 
!  Notes:
! -------
!  1)        The column dimension of A is not needed. 
!  2)        algorithm in place (b, jb, ib, can be the same as
!            a, ja, ia, on entry). See comments for parameter job.
! 
!  coded by Y. Saad. Latest version July, 19, 1990
! -----------------------------------------------------------------
        logical :: test
! 
!      copy integer arrays into b's data structure if required
! 
        if (job /= 0) then 
                nnz = ia(nrow+1)-1
                do k=1, nnz
                        jb(k) = ja(k)
                        b(k) = a(k) 
                end do
 
                do k=1, nrow+1
                        ib(k) = ia(k)
                end do
        endif 
! 
!      get positions of diagonal elements in data structure.
!      
        call diapos (nrow,ja,ia,iw)
!      
!      count number of holes in diagonal and add diag(*) elements to
!      valid diagonal entries.
!      
        icount = 0
 
        do j=1, nrow
                if (iw(j) == 0) then
                        icount = icount+1
                else
                        b(iw(j)) = a(iw(j)) + diag(j) 
                endif
        end do
!      
!      if no diagonal elements to insert return
!      
        if (icount == 0) return
!      
!      shift the nonzero elements if needed, to allow for created 
!      diagonal elements. 
!      
        ko = ib(nrow+1)+icount
!      
!      copy rows backward
!      
        do ii=nrow, 1, -1 
!      
!      go through  row ii
!      
                k1 = ib(ii)
                k2 = ib(ii+1)-1 
                ib(ii+1) = ko
                test = (iw(ii) .eq. 0) 
                
                do k = k2,k1,-1 
                        j = jb(k)
                        if (test .and. (j < ii)) then 
                                test = .false. 
                                ko = ko - 1
                                b(ko) = diag(ii) 
                                jb(ko) = ii
                                iw(ii) = ko
                        endif
                        
                        ko = ko-1
                        b(ko) = a(k) 
                        jb(ko) = j
                end do
!      diagonal element has not been added yet.
                if (test) then
                        ko = ko-1
                        b(ko) = diag(ii) 
                        jb(ko) = ii
                        iw(ii) = ko
                endif
        end do
        ib(1) = ko 
        return
! -----------------------------------------------------------------------
! ------------end-of-apldiag---------------------------------------------
 end subroutine apldia
