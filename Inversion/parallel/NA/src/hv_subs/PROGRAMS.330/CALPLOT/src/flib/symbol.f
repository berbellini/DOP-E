c---------------------------------------------------------------------c
c                                                                     c
c      COMPUTER PROGRAMS IN SEISMOLOGY                                c
c      VOLUME I                                                       c
c                                                                     c
c      PROGRAM: SYMBOL                                                c
c                                                                     c
c      COPYRIGHT (C)  1997 R. B. Herrmann                             c
c                                                                     c
c      Department of Earth and Atmospheric Sciences                   c
c      Saint Louis University                                         c
c      3507 Laclede Avenue                                            c
c      St. Louis, Missouri 63103                                      c
c      U. S. A.                                                       c
c                                                                     c
c---------------------------------------------------------------------c
        subroutine symbol (xloc,yloc,height,inbuf,angle,nocar) 
c-----
c       produces either 1) one symbol of choice or 2) character string
c
c       xloc    x-coordinate lower left corner of symbol
c               if nocar < 0 and symbol is centered, center x point
c       yloc    y-coordinate lower left corner of symbol or string
c       height  height in inches of symbol
c       inbuf   nocar > 0 inbuf is ascii string of characters
c               nocar long
c               nocar < 0 inbuf is a single character of an
c                         integer equivalent, e.g., char(inteq)
c
c                         of the symbol to be plotted
c       angle   angle of rotation of symbol in degrees
c       nocar   > 0 number of characters in inbuf string
c               -1 pen is up during move after which single symbol
c                  is plotted
c               -2 pen is down, a line is drawn from present position
c                  to point where symbol is drawn
c-----
        real*4 xloc, yloc, height, angle
        integer*4 nocar
        integer*4 ix, iy, iht, iang, nchar
        character inbuf*(*) 
        common/Scplot/x0,y0
        common/Xcplot/xold,yold,xcur,ycur
        common/Ycplot/xstp,ystp
        common/Zcplot/iunit,ifont

            if(xloc .lt. 999.0) x0 = xloc
            if(yloc .lt. 999.0) y0 = yloc
            if(nocar .gt. 0 .and. len(inbuf).lt.nocar)then
                nchar = len(inbuf)
            else
                nchar = nocar
            endif
            xcur=x0
            ycur=y0
c-----
c       convert into inches if necessary
c-----
            xn=xcur*xstp
            yn=ycur*ystp
c-----
c       current position is in inches, ala calcomp
c       but to use unix plot filters we need integer*2
c       the factor 1000 is really 11000 counts / 11.0 inches
c-----
            xx = 1000. * (xn + xold)
            yy = 1000. * (yn + yold)
            if(xx.gt.1000000000.0)then
                ix = 1000000000
            else if(xx.lt.-1000000000.0)then
                ix = -1000000000
            else
                ix = xx
            endif
            if(yy.gt.1000000000.0)then
                iy = 1000000000
            else if(yy.lt.-1000000000.0)then
                iy = -1000000000
            else
                iy = yy
            endif
c-----
c       convert height to proper units
c-----
                xht = (1000.0 * height * xstp)
            if(xht.lt.0.0)then
                iht = 0
            else if(xht.gt.1000000000.0)then
                iht = 1000000000
            else
                iht = xht
            endif
            iang = angle
            call dfgsym(ix,iy,iht,iang,nchar,inbuf)
c-----
c       adjust current position
c-----
            if(nocar .gt. -1)then
                ang = angle * 3.1415927/180.0
                x0 = x0 + height * cos(ang) * nocar
                y0 = y0 + height * sin(ang) * nocar
                xcur=x0
                ycur=y0
            endif
        return
        end
