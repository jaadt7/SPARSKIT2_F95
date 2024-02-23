!*==psgrid.f90 processed by SPAG 8.04RA 10:54 23 Feb 2024
!!SPAG Open source Personal, Educational or Academic User Clemson University  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
!-----------------------------------------------------------------------
SUBROUTINE psgrid(Npts,Ja,Ia,Xx,Yy,Title,Ptitle,Size,Munt,Iunt)
!-----------------------------------------------------------------------
!     plots a symmetric graph defined by ja,ia and the coordinates
!     xx(*),yy(*)
!-----------------------------------------------------------------------
! npts   = number of points in mesh
! ja, ia = connectivity of mesh -- as given by pattern of sparse
!            matrix.
! xx, yy = cordinates of the points.
!
! title  = character*(*). a title of arbitrary length to be printed
!          as a caption to the figure. Can be a blank character if no
!          caption is desired.
!
! ptitle = position of title; 0 under the drawing, else above
!
! size   = size of the drawing
!
! munt   = units used for size : 'cm' or 'in'
!
!     iunt   = logical unit number where to write the matrix into.
!-----------------------------------------------------------------------
   USE ISO_FORTRAN_ENV                 
   IMPLICIT NONE
!
! Dummy argument declarations rewritten by SPAG
!
   INTEGER , INTENT(IN) :: Npts
   INTEGER , INTENT(IN) , DIMENSION(*) :: Ja
   INTEGER , INTENT(IN) , DIMENSION(*) :: Ia
   REAL(REAL64) , INTENT(IN) , DIMENSION(Npts) :: Xx
   REAL(REAL64) , INTENT(IN) , DIMENSION(Npts) :: Yy
   CHARACTER(*) :: Title
   INTEGER , INTENT(IN) :: Ptitle
   REAL , INTENT(IN) :: Size
   CHARACTER(2) , INTENT(IN) :: Munt
   INTEGER , INTENT(IN) :: Iunt
!
! Local variable declarations rewritten by SPAG
!
   REAL :: botmrgn , delt , dimen , fnstit , frlw , lrmrgn , paperx , scfct , siz , u2dot , xdim , xi , xj , xl , xr , xtit , yb , &
         & ydim , yi , yj , yt , ytit , ytitof
   REAL , SAVE :: conv , haf , zero
   INTEGER :: ii , j , k , ltit , nr
   INTEGER , EXTERNAL :: lenstr
   REAL(REAL64) :: xmax , xmin , ymax , ymin
!
! End of declarations rewritten by SPAG
!
!-----------------------------------------------------------------------
!     local variables --------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
   DATA haf/0.5/ , zero/0.0/ , conv/2.54/
   siz = Size
!
!     get max and min dimensions
!
   xmin = Xx(1)
   xmax = xmin
   ymin = Yy(1)
   ymax = ymin
   DO j = 2 , Npts
      xmax = max(xmax,Xx(j))
      xmin = min(xmin,Xx(j))
      ymax = max(ymax,Yy(j))
      ymin = min(ymin,Yy(j))
   ENDDO
!-----------------------------------------------------------------------
!      n = npts
   nr = Npts
   xdim = xmax - xmin
   ydim = ymax - ymin
   dimen = max(xdim,ydim)
!-----------------------------------------------------------------------
   PRINT * , ' xmin' , xmin , ' xmax' , xmax
   PRINT * , ' ymin' , ymin , ' ymax' , ymax , ' dimen ' , dimen
!-----------------------------------------------------------------------
!
! units (cm or in) to dot conversion factor and paper size
!
   IF ( Munt=='cm' .OR. Munt=='CM' ) THEN
      u2dot = 72.0/conv
      paperx = 21.0
   ELSE
      u2dot = 72.0
      paperx = 8.5*conv
      siz = siz*conv
   ENDIF
!
! left and right margins (drawing is centered)
!
   lrmrgn = (paperx-siz)/2.0
!
! bottom margin : 2 cm
!
   botmrgn = 2.0/dimen
!     scaling factor
   scfct = siz*u2dot/dimen
!     frame line witdh
   frlw = 0.25/dimen
!     font siz for title (cm)
   fnstit = 0.5/dimen
   ltit = lenstr(Title)
!
!     position of title : centered horizontally
!                     at 1.0 cm vertically over the drawing
   ytitof = 1.0/dimen
   xtit = paperx/2.0
   ytit = botmrgn + siz*nr/dimen + ytitof
! almost exact bounding box
   xl = lrmrgn*u2dot - scfct*frlw/2
   xr = (lrmrgn+siz)*u2dot + scfct*frlw/2
   yb = botmrgn*u2dot - scfct*frlw/2
   yt = (botmrgn+siz*ydim/dimen)*u2dot + scfct*frlw/2
   IF ( ltit>0 ) yt = yt + (ytitof+fnstit*0.70)*u2dot
! add some room to bounding box
   delt = 10.0
   xl = xl - delt
   xr = xr + delt
   yb = yb - delt
   yt = yt + delt
!
! correction for title under the drawing
!
   IF ( Ptitle==0 .AND. ltit>0 ) THEN
      ytit = botmrgn + fnstit*0.3
      botmrgn = botmrgn + ytitof + fnstit*0.7
   ENDIF
!
!     begin output
!
   WRITE (Iunt,99001) '%!'
   WRITE (Iunt,99001) '%%Creator: PSPLTM routine'
   WRITE (Iunt,99003) '%%BoundingBox:' , xl , yb , xr , yt
99003 FORMAT (A,4F9.2)
   WRITE (Iunt,99001) '%%EndComments'
   WRITE (Iunt,99001) '/cm {72 mul 2.54 div} def'
   WRITE (Iunt,99001) '/mc {72 div 2.54 mul} def'
   WRITE (Iunt,99001) '/pnum { 72 div 2.54 mul 20 string'
   WRITE (Iunt,99001) 'cvs print ( ) print} def'
   WRITE (Iunt,99001) '/Cshow {dup stringwidth pop -2 div 0 rmoveto show} def'
!
!     we leave margins etc. in cm so it is easy to modify them if
!     needed by editing the output file
!
   WRITE (Iunt,99001) 'gsave'
   IF ( ltit>0 ) THEN
      WRITE (Iunt,*) '/Helvetica findfont' , fnstit , ' cm scalefont setfont'
      WRITE (Iunt,*) xtit , ' cm' , ytit , ' cm moveto'
      WRITE (Iunt,'(3A)') '(' , Title(1:ltit) , ') Cshow'
   ENDIF
!
   WRITE (Iunt,*) lrmrgn , ' cm ' , botmrgn , ' cm translate'
   WRITE (Iunt,*) siz , ' cm ' , dimen , ' div dup scale '
!
!     draw a frame around the matrix  // REMOVED
!
!       del = 0.005
!       del2 = del*2.0
!       write(iunt,*) del, ' setlinewidth'
!       write(iunt,10) 'newpath'
!       write(iunt,11) -del2, -del2, ' moveto'
!       write(iunt,11) dimen+del2,-del2,' lineto'
!       write(iunt,11) dimen+del2, dimen+del2, ' lineto'
!       write(iunt,11) -del2,dimen+del2,' lineto'
!       write(iunt,10) 'closepath stroke'
!
!----------- plotting loop ---------------------------------------------
   WRITE (Iunt,*) ' 0.01 setlinewidth'
!
   DO ii = 1 , Npts
 
!     if (mask(ii) .eq. 0) goto 1
      xi = Xx(ii) - xmin
      yi = Yy(ii) - ymin
!     write (iout+1,*) ' ******** ii pt', xi, yi, xmin, ymin
      DO k = Ia(ii) , Ia(ii+1) - 1
         j = Ja(k)
         IF ( j>ii ) THEN
            xj = Xx(j) - xmin
            yj = Yy(j) - ymin
!     write (iout+1,*) ' j pt -- j= ',j, 'pt=', xj, yj
!
!     draw a line from ii to j
!
            WRITE (Iunt,99002) xi , yi , ' moveto '
            WRITE (Iunt,99002) xj , yj , ' lineto'
            WRITE (Iunt,99001) 'closepath stroke'
         ENDIF
      ENDDO
   ENDDO
!-----------------------------------------------------------------------
   WRITE (Iunt,99001) 'showpage'
   RETURN
!
99001 FORMAT (A)
99002 FORMAT (2F9.2,A)
99004 FORMAT (2F9.2,A)
!-----------------------------------------------------------------------
END SUBROUTINE psgrid
