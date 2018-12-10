c---------------------------------------------------------------------c
c                                                                     c
c      COMPUTER PROGRAMS IN SEISMOLOGY                                c
c      VOLUME I                                                       c
c                                                                     c
c      PROGRAM: NUMBER                                                c
c                                                                     c
c      COPYRIGHT (C)  1997 R. B. Herrmann                             c
c                                                                     c
c      Department of Earth and Atmospheric Sciences                   c
c      Saint Louis University                                         c
c      221 North Grand Boulevard                                      c
c      St. Louis, Missouri 63103                                      c
c      U. S. A.                                                       c
c                                                                     c
c---------------------------------------------------------------------c
      subroutine number (xpage,ypage,height,fpn,angle,ndec)
c.....     xpage,ypage coordinates of lower left corner of number.
c.....     height   height of plotted number.
c.....     fpn      floating point number to be plotted.
c.....     angle    angle at which number is plotted, in degrees.
c.....     ndec     number of decimal places to be drawn.
c.....     this version of number requires the symbol version with
c.....     999. x, y feature, and  nc = 0 feature.
c.....
c.....     ndec .ge. 1000 .le. 1009 invokes an exponential format
c.....     of the Fortran Format Ew.d, where d = ndec - 1000
c.....     and w = 7 + d. The output will look like
c.....     sd.dddddEsdd, where s is the sign and d is a digit
c.....     In all cases, no more than 9 significant figures are
c.....     permitted to the right of the decimal place.
c.....
c.....     ndec .ge. 2000 .le. 2009 invokes exponential but
c.....     plotted in scientific notation, e.g., 1.0 10 ^ exp
      character num*20
      character cexp*4
        common/Scplot/x0,y0
        common/Xcplot/xold,yold,xcur,ycur

        if(xpage .lt. 999.0) then
            x0 = xpage
        endif
        if(ypage .lt. 999.0) then
            y0 = ypage
        endif
c-----
c    to avoid dregs in memory null first character
c-----
      num(1:1)= ' '
      ii=0
c-----
      n = 0
        fp = fpn
        if(fpn.lt.0.0)then
                n=n+1
                num(n:n)='-'
                fp=abs(fp)
        endif
      if(ndec.lt.1000)then
                nzflag = 0
                mdec = -ndec
c               if(ndec.lt.0)then
c                       mdec = mdec - 1
c               endif
                if(ndec.eq.-1)then
                        nzflag = 1
                        mdec = 1
                endif
                call fpack(num,n,fp,mdec,19,nzflag)
                nman = 0
                nexp = 0
        else if(ndec.ge.1000)then
                       mdec = 1
                if(ndec .ge.1000 .and. ndec .le. 1009)then
                       mdec = ndec -1000
                else if(ndec .ge.2000 .and. ndec .le. 2009)then
                       mdec = ndec -2000
                else 
                       mdec = 9
                endif
                mdec = - mdec
                if(fp.eq.0.0)then
                        iex=0
                else
                        aex = alog10(fp)
                        iex = aex
c----------careful check of integer exponent
                        if(iex.lt.0.and.float(iex).gt.aex)iex=iex-1
                endif
                fp = fp * 10.** (- iex)
                nzflag = 0
                call fpack(num,n,fp,mdec,14,nzflag)
c---put in exponent assuming iex < 100
                nman=n
                n=n+1
                num(n:n)='E'
                if(iex.lt.0)then
                        n=n+1
                        num(n:n)='-'
                        iex = -iex
                else
                        n=n+1
                        num(n:n)='+'
                endif
                nexp = n
                n=n+1
                jex = iex
                jex = mod(jex/10,10)
                num(n:n) = char(ichar('0' )+ jex)
                jex = mod(iex,10)
                n = n + 1
                num(n:n) = char(ichar('0' )+ jex)
        endif
        if(ndec .le. 1009)then
                call symbol(x0,y0,height,num,angle,n)
        else if(ndec .ge. 2000 .and. ndec.le.2009)then
c-----
c       save the exponent, replace the E-EX with 'x10'
c-----
                cexp(1:3) = num(nexp:n)
                num(nman+1:nman+3)='x10'
            xx0 = x0
            yy0 = y0
                call symbol(x0,y0,height,num(1:nman+3),
     1               angle,nman+3)
c-----
c       get the proper position because of rotation
c-----
                ang = angle*3.1415927/180.0
                ca = cos(ang)
                sa = sin(ang)
                xl = height*(nman+3)
                yl = 0.6*height
            xx0 = xx0 + xl*ca
            yy0 = yy0 + xl*sa
                xxx0 = xx0  - yl*sa
                yyy0 = yy0  + yl*ca
                ht = 0.7*height
                call symbol(xxx0,yyy0,ht,cexp(1:3),angle,3)
c-----
c           now position at end of proper string
c-----
            xl =  3.0*ht
            xcur = xx0 + xl*ca
            ycur = yy0 + xl*sa
            x0 = xcur
            y0 = ycur
            
c-----
        endif
        return
        end
        subroutine fpack(num,n,fpv,mdec,mwid,nzflag)
        character num*(*)
        fp = fpv
        nzflg = nzflag
c-----since we have a maximum field width of mwid and since
c-----the . takes a position, we can only have mwid -1 figures
      m = mdec 
        maxm = 9
        if(m.gt.maxm)m=maxm
        if(m.lt. -maxm)m= -maxm
        mm = m
        if(m.gt.0) mm = mm - 1
c----do a simple rounding up
        fp = fp + 0.5 * 10.** (mm)
c---- get position of first significant figure
c
c     5 4 3 2 1 0 -1 -2 -3
c
c     d d d d d .  d  d  d
c
c----

        iex = alog10(fp) 
        if(fp.ge.1.0)then
                iex = iex + 1
        else
                iex = iex - 1
        endif
c----
c----   procede from most significant down
c----   but never less than 10**-9
c----
        
        if(fp.le.1.0e-9)then
                fp = 0.0
                iex = -9
        endif
        ilw = mdec
        if(fp.ge.1.0)then
                iup=iex
        else
                iup = 1
        endif
        if(m.le.0)then
                ilw= m+1
        else
                ilw = m 
        endif
        if(iex.gt.9)then
                ilw = 1
                nzflg = 1
        endif
        do 100 ipos=iup,ilw,-1
                k = fp * 10.**( -ipos +1 )
                if(n.lt.mwid)then
                        n = n + 1
                        num(n:n) = char(ichar('0')+k)
                endif
                fp = fp - (float(k) * 10.**(ipos-1))
        if(ipos.eq.1.and.nzflg.eq.0.and.n.lt.mwid)then
                        n = n + 1
                        num(n:n) = '.'
                endif
 100    continue
        return
        end
