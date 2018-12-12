c
c     POSTPAK    - low level plotting routines
c                 to support Postcript output
c                 using plotpak calls
c                 together with some typset calls
c                 implemented to exploit the font 
c                 capabilities of the Laserwriter II
c
c		  This version is an attempt to combine the
c	 	  modifications of Brian and Phil into a 
c		  single library. It includes the font
c		  changes of Brian's version as well as
c		  Phil's minor alterations to stop
c		  padding strings out to 80 characters.
c
c						MS 1/12/92
c     Current Version - Dec 1992     
c
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
C*    TYPSET SIMULATION
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~bk~
      SUBROUTINE TYPSET( XB0 , YB0 )
C
c++   This routine resets the origin in plotter units to
c++   (xbo ,yb0) for an interface to the R.S.E.S. plotter
c++   libraries based on the PLOTPAK calls
c     also set macros for Laserwriter fonts
c
c     1      LucidaSans
c     2      LucidaSans-Bold
c     3      Palatino-Roman  
c     4      Times-Roman
c     5      LucidaSans-Italic
c     6      Palatino-Italic
c     7      Times-Italic
c     8      Palatino-Bold
c     9      Palatino-BoldItalic
c     10     Symbol
c
C
      COMMON/P00000/LPLOT,IROT,IL34,A,B,C,D,ASP,THET
      COMMON/L00000/PSCA,xo,yo
      COMMON/T00000/zxb0,zyb0,fchar,ffac
      character*4 fchar
C
      ZXB0 = XB0
      ZYB0 = YB0
      write(LPLOT,*) '/f0 {/Helvetica findfont '
      write(LPLOT,*) ' exch scalefont setfont} def'
      write(LPLOT,*) '/f1 {/LucidaSans findfont '
      write(LPLOT,*) ' exch scalefont setfont} def'
      write(LPLOT,*) '/f2 {/LucidaSans-Bold findfont '
      write(LPLOT,*) ' exch scalefont setfont} def'
      write(LPLOT,*) '/f3 {/Palatino-Roman findfont '
      write(LPLOT,*) ' exch scalefont setfont} def'
      write(LPLOT,*) '/f4 {/Times-Roman findfont '
      write(LPLOT,*) ' exch scalefont setfont} def'
      write(LPLOT,*) '/f5 {/LucidaSans-Italic findfont '
      write(LPLOT,*) ' exch scalefont setfont} def'
      write(LPLOT,*) '/f6 {/Palatino-Italic findfont '
      write(LPLOT,*) ' exch scalefont setfont} def'
      write(LPLOT,*) '/f7 {/Times-Italic findfont '
      write(LPLOT,*) ' exch scalefont setfont} def'
      write(LPLOT,*) '/f8 {/Helvetica-Bold findfont '
      write(LPLOT,*) ' exch scalefont setfont} def'
      write(LPLOT,*) '/f9 {/Palatino-BoldItalic findfont '
      write(LPLOT,*) ' exch scalefont setfont} def'
      write(LPLOT,*) '/fA {/Symbol findfont '
      write(LPLOT,*) ' exch scalefont setfont} def'
      RETURN
      END
c------------------------------------------------------------bk
      subroutine zpick(ifont,n,is)
c     
c     choose a Laserwriter font
c
c     0      Helvetica
c     1      LucidaSans
c     2      LucidaSans-Bold
c     3      Palatino-Roman
c     4      Times-Roman
c     5      LucidaSans-Italic
c     6      Palatino-Italic
c     7      Times-Italic
c     8      Palatino-Bold
c     9      Palatino-BoldItalic
c     10     Symbol
c
      COMMON/P00000/LPLOT,IROT,IL34,A,B,C,D,ASP,THET
      COMMON/L00000/PSCA,xo,yo
      COMMON/T00000/zxb0,zyb0,fchar,ffac
      
      character*4 fchar
      if (ifont.lt.10) then
        write(fchar,51) ifont
  51    format(' f',i1,' ')
      elseif (ifont.eq.10) then
        write(fchar,52) 
  52    format(' fA ')
      endif
      return
      end
c-------------------------------------------------------------bk
      SUBROUTINE TYPNUM(X,Y,SIZE1,FNUM1,ANGLE,NDEC1)
C
C-- THIS ROUTINE CONVERTS FNUM1 TO THE APPROPRIATE FIXED DECIMAL
C-- EQUIVALENT AND PLOTS IT TO ANY DEGREE OF ACCURACY WITH ROUNDING.
C
      COMMON/P00000/LPLOT,IROT,IL34,A,B,C,D,ASP,THET
      COMMON/L00000/PSCA,xo,yo
      COMMON/T00000/zxb0,zyb0,fchar,ffac
      character*4 fchar
      CHARACTER*8 IFORM
      CHARACTER*40 CWORD
      integer iword(10)
      equivalence (iword(1),cword)
C
C    create the format expression in IFORM
C
      NDEC = NDEC1
      SIZE=ABS(SIZE1)
      RN=FNUM1
      iadd = 1
      if(ndec.ge.0) iadd=2
      if(rn.lt.0) iadd=iadd+1
      if(abs(rn).lt.1) iadd=iadd+1
      if(rn.eq.0) then
        if(ndec.lt.0) np = 1
        if(ndec.ge.0) np = 2
      else
        pow = alog10(abs(rn))
        np = int(pow+0.01)+iadd
      endif
      if(ndec.lt.0) nsf = -np
      if(ndec.ge.0) nsf = 10*(np+ndec)+ndec
c       write(6,*) nsf,rn
c 
      IF(NSF.LT.0)GO TO 20
      ITOT=NSF/10
      IDPL=MOD(NSF,10)
      Write(IFORM,55)ITOT,IDPL
   55 FORMAT('(F',I2,'.',I1,')')
      Write(CWORD,IFORM)RN
      GO TO 30
C
C    for integer format
C
   20 ITOT=-NSF
      Write(IFORM,65)ITOT
   65 FORMAT('(I',I2,')   ')
      IR=INT(RN)
      Write(CWORD,IFORM)IR
c      write(6,fmt='(a)') iword
c     send to plotter
   30 CALL TYPSTR(X,Y,SIZE,IWORD,ANGLE,ITOT)
      RETURN
      END
c---------------------------------------------------------bk
c     SUBROUTINE TYPSTR(X,Y,SIZE,KWORD,ANGL,NCHAR) old version
c
      SUBROUTINE TYPSTR(X,Y,SIZE,IWORD,ANGL,NCHAR)
C
C     writes a Hollerith string on the plot--plotter units
C
c** Modified to write only a string of length min(nchar,80), rather
c** than always padding with zeros to produce 80 chars - Phil Cummins 1/92
c
c					Modified 19/1/93 use of integer strings
c					removed to fix spurious bug. 
c				
      COMMON/P00000/LPLOT,IROT,IL34,A,B,C,D,ASP,THET
      COMMON/L00000/PSCA,xo,yo
      COMMON/T00000/zxb0,zyb0,fchar,ffac
      character*4 fchar
c     integer kword(1),jword(20)
      CHARACTER*80 IWORD
c					Modified by MS 19/1/93 to make
c					iword the same length as cword
      CHARACTER*80 CWORD
c
c     equivalence (jword(1),iword)
c				
c     do 10 j=1,20
c      jword(j) = kword(j)
c10   continue
c
c     write(6,*)' typstr: angl',angl,' nchar',nchar
c     write(6,*)' typstr: x',x,'y',y,'size',size
c     write(6,*)' typstr: string'
c     write(6,fmt='(a80)')iword
      nch = nchar
      IF(NCHAR.GT.77)NCH=77
C
C     select character orientation
C
      IF(ANGL.EQ.0.0)GO TO 40
C
C     select character size
C
c  40 SZ=ABS(SIZE)*PSCA*1.5
   40 SZ=ABS(SIZE)*PSCA
c     write(6,*)' sz ',sz,' psca ',psca,' size ',size
      do 50 k=1,nch
        cword(k:k) = iword(k:k)
   50 continue
c      do 51 k=nch+1,60
c        cword(k:k) = ' '
c   51 continue
C
C      move pen to symbol location
C
      IP=3
      IF(SIZE.LT.0.0)IP=-3
      CALL PLOT(X+zxb0,Y+zyb0,IP)
      write(LPLOT,fmt='(f9.3,a)') angl,' rotate'
C
c     write character string
C
      write(LPLOT,fmt='(f9.3,x,a)') sz,fchar
c     write(LPLOT,*) 'The value of sz = ',sz
c     write(LPLOT,*) 'The value of size = ',size
c     write(LPLOT,*) 'The value of psca = ',psca
c     write(LPLOT,*) 'The value of x = ',x
c     write(LPLOT,*) 'The value of y = ',y
      write(LPLOT,fmt='(x,a,a,a,/,a)') '(',cword(1:nch),')',' show'
C
C     reset character orientation if necessary
C
      IF(ANGL.EQ.0.0)RETURN
      bngl = -angl 
      write(LPLOT,fmt='(f9.3,a)') bngl,' rotate'
      RETURN
      END
c
c     LASPAK    - low level plotting routines
c                 to support Postcript output
c                 using plotpak calls
c     Current Version - Sept. 1988     blnk
c     Last modified -  03/04/89          mb    
c     pimask added 21/07/89            blnk
C-----------------------------------------------------------------bk
      subroutine ASPECT(RASP)
C
C     sets width to height ratio for characters to RASP
C
      COMMON/P00000/LPLOT,IROT,IL34,A,B,C,D,ASP,THET
      COMMON/L00000/PSCA,xo,yo
      ASP = RASP
      return
      end
C-----------------------------------------------------------------
      subroutine CIRCLE(RADIUS,NSIDES)
C
C    draws circle centred at current pen location, with
C    circumference divided into NSIDES straight segments
C
      COMMON/P00000/LPLOT,IROT,IL34,A,B,C,D,ASP,THET
      COMMON/L00000/PSCA,xo,yo
c     CHARACTER*40 IWORD
C
      cx = xo
      cy = yo
      nsid = NSIDES
      if(nsid.EQ.0) nsid = 72
      rpa = radius
      ANG = 6.283185308/float(nsid)
      xv = rpa+cx
      yv = cy
      call plot(xv,yv,3)
      sta=0.0
      do 30 i=1,nsid
        sta = sta+ang
        xv = rpa*cos(sta)+cx
        yv = rpa*sin(sta)+cy
        call plot(xv,yv,2)
 30   continue
      xo = cx
      yo = cy
      RETURN
      END
C----------------------------------------------------------------------
      subroutine CSYMBL(X,Y,IP,SIZE,INT)
C
C      writes a centered symbol at location (X,Y). The symbol is
C      is selected from the list below by INT for 1<INT<10
C      and is circle,triangle,square,pentagon,hexagon,heptagon,
C      octagon for 11<INT<17
C      if INT lies in the range 21<INT<27 then
C      a filled symbol is plotted
C
c
      COMMON/P00000/LPLOT,IROT,IL34,A,B,C,D,ASP,THET
      COMMON/L00000/PSCA,xo,yo
      CHARACTER*1 ISYM
c     CHARACTER*40 IWORD
      DIMENSION ISYM(10)
      DIMENSION ICIR(10)
      DATA ISYM/'O','X','*','+','#','$','@','8','H','Z'/
      DATA ICIR/20,3,4,5,6,7,8,9,10,11/
      iny = INT
      ifill = 0
      if(iny.GT.20) then
        iny = iny-10
        ifill=1
      endif
C     write (6,*) iny,ifill
      IF(iny.GE.11)GO TO 20
C
C    select character size
C
      call plotu(x,y,ip)
      call symbol(x-0.5*size,y-0.5*size,size,isym(iny),0.0,1)
      return
c
C     move pen to symbol location, symbol is written after move
C
   20 CALL PLOTU(X,Y,IP)
      CALL CIRCLE(SIZE*0.75,ICIR(iny-10))
      if(ifill.eq.1) then
        write(LPLOT,*) 'closepath fill'
      endif
      RETURN
C
       END
C---------------------------------------------------------------bk-
      subroutine DASHLN(LDASH,LPAT)
C
C     defines  the style for line drawing
C      ldash < 0   -  reset to solid line
C         or = 12
C
C      ldash = 0   -  dots at calling points
C            = 1   -  dots
C            = 2   -  half dash
C            = 3   -  long dash
C            = 4   -  chain dotted
C            = 5   -  long and short
C            = LPL,*   -  long and two short
C
C       lpat - percentage of diagonal of paper used for
C              a pattern
C              if (lpat.eq.0 ) lpat = 2
C
      COMMON/P00000/LPLOT,IROT,IL34,A,B,C,D,ASP,THET
      COMMON/L00000/PSCA,xo,yo
c       CHARACTER*40 IWORD
C
        integer*4 ldash,lpat
        write(LPLOT,*) ' stroke' 
        if(ldash.lt.0)  write(LPLOT,*) '[] 0 setdash'
        if(ldash.eq.12) write(LPLOT,*) '[] 0 setdash'
c
        if(ldash.eq.0) write(LPLOT,*) '[2 8] 0 setdash'
        if(ldash.eq.1) write(LPLOT,*) '[2 8] 0 setdash'
        if(ldash.eq.2) write(LPLOT,*) '[4 4] 0 setdash'
        if(ldash.eq.3) write(LPLOT,*) '[8 8] 0 setdash'
        if(ldash.eq.4) write(LPLOT,*) '[6 2 2 2] 0 setdash'
        if(ldash.eq.5) write(LPLOT,*) '[8 4 4 4] 0 setdash'
        if(ldash.eq.6) write(LPLOT,*) '[6 4 4 4 4 4] 0 setdash'
        xv=psca*xo
        yv=psca*yo
        write(LPLOT,fmt='(f9.3,1x,f9.3,3a)') xv,yv,' pM'
C
        return
        end
C---------------------------------------------------------------bk-
      subroutine ITALIC(THETA)
C
C     defines angle of slant for labels in degrees
c     null for laspak
C
      return
      end
C------------------------------------------------------------------
      subroutine NUMBER(X,Y,SIZE,RN,ANGL,NSF)
C
C     writes a number on the plot: if NSF=klm, format is Fkl.m
C      if NSF=-lm, RN is fixed to an integer and format is Ilm.
C
      COMMON/P00000/LPLOT,IROT,IL34,A,B,C,D,ASP,THET
      COMMON/L00000/PSCA,xo,yo
      CHARACTER*8 IFORM
      CHARACTER*40 IWORD
C
C    create the format expression in IFORM
C
      IF(NSF.LT.0)GO TO 20
      ITOT=NSF/10
      IDPL=MOD(NSF,10)
      Write(IFORM,55)ITOT,IDPL
   55 FORMAT('(F',I2,'.',I1,')')
      Write(IWORD,IFORM)RN
      GO TO 30
C
C    for integer format
C
   20 ITOT=-NSF
      Write(IFORM,65)ITOT
   65 FORMAT('(I',I2,')   ')
      IR=IFIX(RN)
      Write(IWORD,IFORM)IR
C
C     encode number and send to plotter
C
   30 CALL SYMBOL(X,Y,SIZE,IWORD,ANGL,ITOT)
      RETURN
      END
C----------------------------------------------------------------
      subroutine PLOT(X,Y,I)
C
C     Raises (I=3) or lowers (I=2) pen and moves to coordinates
C      (X,Y) if I>0 or to current position plus (X,Y) if I<0
C
      COMMON/P00000/LPLOT,IROT,IL34,A,B,C,D,ASP,THET
      COMMON/L00000/PSCA,xo,yo
c     CHARACTER*40 IWORD
      DATA IUP/1/
      II=IABS(I)
C
C     Rotate plot by 90 degrees if necessary
C
      XP=X
      YP=Y
      IF(IROT.EQ.0)GO TO 30
      YP=X
      XP=-Y
      if(il34.eq.0) then
       if(I.GT.0) XP = 27.2-Y
      else if (il34.eq.1) then
       if(I.GT.0) XP = 40.1-Y
      end if
C
C    convert to points
C
   30 XV= PSCA*XP
      YV= PSCA*YP
C
C     plot
C
      if(I.eq.2)  write(LPLOT,fmt='(f11.3,1x,f11.3,a)') xv,yv,' pL'
      if(I.eq.3)  write(LPLOT,fmt='(f11.3,1x,f11.3,a)') xv,yv,' pM'
      if(I.eq.-2) write(LPLOT,fmt='(f11.3,1x,f11.3,a)') xv,yv,' rL'
      if(I.eq.-3) write(LPLOT,fmt='(f11.3,1x,f11.3,a)') xv,yv,' rM'
c 
c      if(I.gt.0) then
       xo = x   
       yo = y
c      endif
      RETURN
      END
C------------------------------------------------------------------
      subroutine PLOTU(X,Y,II)
C
C     scales user coordinates to plotter coordinates
C
      COMMON/P00000/LPLOT,IROT,IL34,A,B,C,D,ASP,THET
      COMMON/L00000/PSCA,xo,yo
      XP=A*X
      YP=C*Y
      IF(II.LT.0)GO TO 10
      XP=XP+B
      YP=YP+D
   10 CALL PLOT(XP,YP,II)
      RETURN
      END
C---------------------------------------------------------gh--bk-
      subroutine SCALE(XMIN,XMAX,PX1,PX2,YMIN,YMAX,PY1,PY2)
C
C      sets up scale factors used in PLOTU and other routines
C
      COMMON/P00000/LPLOT,IROT,IL34,A,B,C,D,ASP,THET
      COMMON/L00000/PSCA,xo,yo
      common/p00001/a1,a2,b1,b2,c1,c2,d1,d2
c
      A=(PX2-PX1)/(XMAX-XMIN)
      B=PX1-A*XMIN
      C=(PY2-PY1)/(YMAX-YMIN)
      D=PY1-C*YMIN
c
      a1 = xmin
      a2 = xmax
      b1 = px1
      b2 = px2
      c1 = ymin
      c2 = ymax
      d1 = py1
      d2 = py2
c
      RETURN
      END
C-----------------------------------------------------------------
      subroutine SYMBOL(X,Y,SIZE,IWORD,ANGL,NCHAR)
C
C     writes a Hollerith string on the plot--plotter units
C
c** Modified to write only a string of length min(nchar,80), rather
c** than always padding with zeros to produce 80 chars - Phil Cummins 1/92
      COMMON/P00000/LPLOT,IROT,IL34,A,B,C,D,ASP,THET
      COMMON/L00000/PSCA,xo,yo
      CHARACTER*80 IWORD
      CHARACTER*80 CWORD
      nch = nchar
      IF(NCHAR.GT.80)NCH=80
C
C     select character orientation
C
      ab=0.0
      if(irot.ne.0) ab=90.0
C
C     select character size
C
   40 SZ=ABS(SIZE)*PSCA*1.5
      do 50 k=1,nch
        cword(k:k) = iword(k:k)
   50 continue
c      do 51 k=nch+1,80
c        cword(k:k) = ' '
c   51 continue
C
C      move pen to symbol location
C
      IP=3
      IF(SIZE.LT.0.0)IP=-3
      CALL PLOT(X,Y,IP)
      ang=angl+ab
      write(LPLOT,fmt='(f9.3,a)') ang,' rotate'
C
c     write character string
C
      write(LPLOT,fmt='(f9.3,x,a)') sz,' hV'
      write(LPLOT,fmt='(x,a,a,a,/,a)') '(',cword(1:nch),')',' show'
C
C     reset character orientation if necessary
C
      bng = -ang 
      write(LPLOT,fmt='(f9.3,a)') bng,' rotate'
10    RETURN
      END
c-----------------------------------------------------------------bk
      subroutine SYMBU(X,Y,SIZE,IWORD,ANGL,NCHAR)
C
C     writes a Hollerith string on the plot--user units
C
      COMMON/P00000/LPLOT,IROT,IL34,A,B,C,D,ASP,THET
      COMMON/L00000/PSCA,xo,yo
      CHARACTER*80 IWORD
      CHARACTER*80 CWORD
      nch = nchar
      IF(NCHAR.GT.80)NCH=80
C
C     select character orientation
C
      IF(ANGL.EQ.0.0)GO TO 40
C
C     select character size
C
   40 SZ=ABS(SIZE)*PSCA*1.5
      do 50 k=1,80
        cword(k:k) = iword(k:k)
   50 continue
      do 51 k=nch+1,80
        cword(k:k) = ' '
   51 continue
C
C      move pen to symbol location
C
      IP=3
      IF(SIZE.LT.0.0)IP=-3
      CALL PLOTU(X,Y,IP)
      write(LPLOT,fmt='(f9.3,a)') angl,' rotate'
C
c     write character string
C
      write(LPLOT,fmt='(f9.3,x,a)') sz,' hV'
      write(LPLOT,fmt='(x,a,a,a,/,a)') '(',cword,')',' show'
C
C     reset character orientation if necessary
C
      IF(ANGL.EQ.0.0)RETURN
      bngl = -angl 
      write(LPLOT,fmt='(f9.3,a)') bngl,' rotate'
      RETURN
      END
c-----------------------------------------------------------------bk
      subroutine FILLTYP(it,spac,ian)
C
C     specifies fill type for shading
C
C     it    - 
C     spac  - line spacing in cm converted to gray tone
C     ian   - 
C
      COMMON/P00000/LPLOT,IROT,IL34,A,B,C,D,ASP,THET
      COMMON/L00000/PSCA,xo,yo
c     CHARACTER*40 IWORD
C
      spec = spac
c     if(spec.gt.0.9) spec=0.9
      write(LPLOT,fmt='(f9.3,a)') spec,' setgray'
C
      return
      end
C---------------------------------------------------------------bk-
      subroutine SHADRT(XI,YI)
C
C     Shades rectangle defined by coordinate increments xinc,yinc
C
        COMMON/P00000/LPLOT,IROT,IL34,A,B,C,D,ASP,THET
        COMMON/L00000/PSCA,xo,yo
c       character*40 iword
C
        xin = a*xi
        yin = c*yi
        call edgert(xi,yi)
        write(LPLOT,*) 'fill'
C
        return
        end
C---------------------------------------------------------------bk-
      subroutine EDGERT(XI,YI)
C
C     Edges rectangle defined by coordinate increments xinc,yinc
C
        COMMON/P00000/LPLOT,IROT,IL34,A,B,C,D,ASP,THET
        COMMON/L00000/PSCA,xo,yo
c       character*40 iword
C
        xin = a*xi
        yin = c*yi
        call plot(xin,0.0,-2)
        call plot(0.0,yin,-2)
        call plot(-xin,0.0,-2)
        call plot(0.0,-yin,-2)
        write(LPLOT,*) 'closepath'
C
        return
        end
C-------------------------------------------------------------------
      subroutine PIMAG4 (xori,yori,xxl,yyl,nsxx,nsyy,arr,
     &                   nth11,nth22,sth11,sth22)
c------
c         pimag4         bitmap image plotting in postscript 
c                        (4 bit representation)
c                        xori - x origin for bit image
c                        yori - y origin for bit image
c                        xxl  - x dimension of block
c                        yyl  - y dimension of block
c                        arr(nsxx,nsyy) - data array
c                        nth11,nth22  - range of levels (1,16)
c                        sth11,sth22  - corresponding array values
c
c     n.b.   grey scale inverted i.e. 0 - white 16 - black
c------
      COMMON/P00000/LPLOT,IROT,IL34,A,B,C,D,ASP,THET
      COMMON/L00000/PSCA,xo,yo
      real*4 xl,yl,sth1,sth2
      integer*4 nsx,nsy,nth1,nth2
      character*1 car,chex(16)
      dimension car(256,256)
      dimension arr(nsxx,nsyy)
      data chex/'0','1','2','3','4','5','6','7','8','9',
     &          'A','B','C','D','E','F'/
c    
      do 100 j=1,256
        do 110 k=1,256
          car(k,j) = ' '
 110    continue
 100  continue
c                                     create Hex map
      xl=xxl
      yl=yyl
      nsx=nsxx
      nsy=nsyy
      nth1=nth11
      nth2=nth22
      sth1=sth11
      sth2=sth22
      nsx1 = nsx
      nr = MOD(nsx,2)
      if(nr.eq.1) nsx1=nsx1+1
      nsx2 = nsx1/2
      do 10 j=1,nsy
        do 20 i=1,nsx1
          npix=int(float(nth1)+float(nth2-nth1)/(sth2-sth1)*
     &         (arr(i,j)-sth1))
          if(npix.lt.nth1) npix=nth1
          if(npix.gt.nth2) npix=nth2
          car(i,j) = chex(17-npix)
c          write(6,*) i,j, arr(i,j),npix,car(i,j)
20      continue
10    continue     
c100   format(Z1)
c                                     write Postscript commands 
      xso = xori*psca
      yso = yori*psca
      xscal = xl*psca
      yscal = yl*psca
      write(LPLOT,*) 'gsave'
      write(LPLOT,81) xso,yso
 81    format(1x,2f8.3,' translate')
      write(LPLOT,82) xscal,yscal
 82    format(1x,2f8.3,' scale')
      write(LPLOT,83) nsx,nsy,nsx,nsy
 83    format(' /imPr {',2i4,' 4 [',i4,' 0 0',i4,' 0 0 ]',/
     ^       ' { currentfile iLn readhexstring pop} image } def')
      write(LPLOT,84) nsx2
 84    format(' /iLn ',i4,' string def',/,' imPr')
      write(LPLOT,85) ((car(ip,jp),ip=1,nsx1),jp=1,nsy)
 85    format(1x,64a1)
      write(LPLOT,*) 'grestore'
c
      return 
      end
C-------------------------------------------------------------------
      subroutine PIMAG8(xori,yori,xxl,yyl,nsxx,nsyy,arr,
     &                   nth11,nth22,sth11,sth22)
c------
c         pimag8         bitmap image plotting in postscript 
c                        (8 bit representation)
c                        xori - x origin for bit image
c                        yori - y origin for bit image
c                        xxl  - x dimension of block
c                        yyl  - y dimension of block
c                        arr(nsxx,nsyy) - data array
c                        nth11,nth22  - range of levels (1,256)
c                        sth11,sth22  - corresponding array values
c
c     n.b.   grey scale inverted i.e. 1 - white 256 - black
c------
      COMMON/P00000/LPLOT,IROT,IL34,A,B,C,D,ASP,THET
      COMMON/L00000/PSCA,xo,yo
      real*4 xl,yl,sth1,sth2
      integer*4 nsx,nsy,nth1,nth2
      character*1 chex(16)
      character*2 car,ccar
      dimension car(256,256)
      dimension arr(nsxx,nsyy)
      data chex/'0','1','2','3','4','5','6','7','8','9',
     &          'A','B','C','D','E','F'/
c    
      do 100 j=1,256
        do 110 k=1,256
          car(k,j) = '  '
 110    continue
 100  continue
c                                     create Hex map
      xl=xxl
      yl=yyl
      nsx=nsxx
      nsy=nsyy
      nth1=nth11
      nth2=nth22
      nmax = max(nth11,nth22)
      nmin = min(nth11,nth22)
      sth1=sth11
      sth2=sth22
      do 10 j=1,nsy
        do 20 i=1,nsx
          npix=int(float(nth1)+float(nth2-nth1)/(sth2-sth1)*
     &         (arr(i,j)-sth1))
c					map array values outside
c					sth1 and sth2 to pen boundaries
          if(npix.lt.nmin) npix=nmin
          if(npix.gt.nmax) npix=nmax
c         if(npix.lt.nth1) npix=nth1 
c         if(npix.gt.nth2) npix=nth2
          npix=256-npix
          np1=npix/16+1
          np2=MOD(npix,16)+1
          ccar(1:1)=chex(np1)
          ccar(2:2)=chex(np2)
          car(i,j) = ccar 
c         write(6,*) i,j, arr(i,j),npix,car(i,j)
20      continue
10    continue     
c                                     write Postscript commands 
      xso = xori*psca
      yso = yori*psca
      xscal = xl*psca
      yscal = yl*psca
      write(LPLOT,*) 'gsave'
      write(LPLOT,81) xso,yso
 81    format(1x,2f8.3,' translate')
      write(LPLOT,82) xscal,yscal
 82    format(1x,2f8.3,' scale')
      write(LPLOT,83) nsx,nsy,nsx,nsy
 83    format(' /imPr {',2i4,' 8 [',i4,' 0 0',i4,' 0 0 ]',/
     ^       ' { currentfile iLn readhexstring pop} image } def')
      write(LPLOT,84) nsx
 84    format(' /iLn ',i4,' string def',/,' imPr')
      write(LPLOT,85) ((car(ip,jp),ip=1,nsx),jp=1,nsy)
 85    format(1x,32a2)
      write(LPLOT,*) 'grestore'
c
      return 
      end

c----------------------------------------------------------------------

      Subroutine pcirclef (x,y,radius)

      COMMON/P00000/LPLOT,IROT,IL34,A,B,C,D,ASP,THET
      COMMON/L00000/PSCA,xo,yo
      DATA IUP/1/
	i=3
      II=IABS(I)
C
C     Rotate plot by 90 degrees if necessary
C
      XP=X
      YP=Y
      IF(IROT.EQ.0)GO TO 30
      YP=X
      XP=-Y
      if(il34.eq.0) then
       if(I.GT.0) XP = 27.2-Y
      else if (il34.eq.1) then
       if(I.GT.0) XP = 40.1-Y
      end if
C
C    convert to points
C
   30 XV= PSCA*XP
      YV= PSCA*YP
C
C     plot
C
c     if(I.eq.3)  write(LPLOT,fmt='(f11.3,1x,f11.3,a)') xv,yv,' pM'
c 
       xo = x   
       yo = y

      write (lplot,*) xv,yv,radius,0,360,' arc fill stroke'

      return
      end

