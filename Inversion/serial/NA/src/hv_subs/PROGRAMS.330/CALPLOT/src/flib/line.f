c---------------------------------------------------------------------c
c                                                                     c
c      COMPUTER PROGRAMS IN SEISMOLOGY                                c
c      VOLUME I                                                       c
c                                                                     c
c      PROGRAM: LINE                                                  c
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
        subroutine line(x,y,n,inc,lintyp,inteq)
c-----
c       x       array of abscissa
c       y       array of ordinates
c       n       number of points to be plotted
c       inc     plot every inc'th point
c               (there are n*inc + inc + 1 points in array)
c       lintyp  != 0 plot symbol inteq at each abs(lintyp) point
c               =0   plot line only
c               <0   plot only symbols, no connecting line
c               >0   line connects symbol
c       inteq   if lintyp != 0, inteq symbol plotted at point
c-----
        real SYMSIZ
       parameter (SYMSIZ = 0.07)
        dimension x(1),y(1)
        integer dj
c-----
c       start plotting from end of line closest to current pen position
c-----
        l1 = n*inc + 1
        l2 = l1 + inc
        l3 = l1 - inc
        x1 = x(l1)
        dx = x(l2)
        y1 = y(l1)
        dy = y(l2)
c-----
c       determine current pen position
c-----
        call where(cx,cy,cf)
c-----
c       determine which end of plot array is closest to current pen
c       position
c-----
        dl = amax1(abs((x(1)-x1)/dx-cx),
     1             abs((y(1)-y1)/dy-cy))
        dr = amax1(abs((x(l3)-x1)/dx-cx),
     1             abs((y(l3)-y1)/dy-cy))
        if(dr.lt.dl)then
                dj = -inc
                j1 = l3
                j2 = 1
        else
                dj = inc
                j1 = 1
                j2 = l3
        endif
c-----
c       do the plotting
c-----
        ia=iabs(lintyp)
        ipen = 3
        icode = -1
        do 100 m=j1,j2,dj
                xn = (x(m)-x1)/dx
                yn = (y(m)-y1)/dy
                call plot(xn,yn,ipen)
                if(ia.gt.0)then
                if(mod(m,ia).eq.1.or.ia.eq.1)call symbol(xn,
     1              yn,SYMSIZ,char(inteq),0.0,icode)
                endif
                if(lintyp.ge.0)ipen=2
  100   continue
c-----
c       penup
c-----
        call plot(xn,yn,3)
        return
        end
