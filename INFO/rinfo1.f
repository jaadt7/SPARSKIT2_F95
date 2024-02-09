      program info1
c---------------------------------------------------------------------- 
c usage    info1.ex < HB_file 
c 
c where info1 is the executable generated by makefile, HB_file is a 
c file containing a matrix stored in Harwell-Boeing matrices.
c Info1 will then dump the information into the standard output.
c
c To use with larger matrices, increase nmax and nzmax.
c---------------------------------------------------------------------- 
      implicit none 
      integer nmax, nzmax
      parameter (nmax = 30000, nzmax = 800000)
      integer ia(nmax+1),ia1(nmax+1),ja(nzmax),ja1(nzmax) 
      real*8  a(nzmax),a1(nzmax),rhs(1)
      character title*72, type*3, key*8, guesol*2 
      logical valued
c
      integer job, iin, nrow,ncol,nnz,ierr, nrhs, iout
c--------------
      data iin /5/, iout/6/
c--------------
      job = 2
      nrhs = 0
      call readmt (nmax,nzmax,job,iin,a,ja,ia, rhs, nrhs,
     +     guesol,nrow,ncol,nnz,title,key,type,ierr)
c---- if not readable return 
      if (ierr .ne. 0) then
         write (iout,100) ierr
 100     format(' **ERROR: Unable to read matrix',/,
     +        ' Message returned fom readmt was ierr =',i3)
         stop
      endif
      valued = (job .ge. 2)
c-------
      call dinfo1(ncol,iout,a,ja,ia,valued,title,key,type,a1,ja1,ia1)
c--------------------end------------------------------------------------ 
c-----------------------------------------------------------------------
      end
