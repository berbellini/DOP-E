c---------------------------------------------------------------------c
c                                                                     c
c      COMPUTER PROGRAMS IN SEISMOLOGY                                c
c      VOLUME I                                                       c
c                                                                     c
c      PROGRAM: PLTLGD                                                c
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
        subroutine pltlgd(x,y,n,x1,y1,deltax,deltay,lintyp,inteq,
     1      ht,nocx,nocy,ipat,xlen)
        real*4 x(1),y(1)
c-----
c       plot y(x) in linear, log or semilog plots. Check
c       for bounds of plot
c       ipat    integer providing bit pattern for dashed lines
c       xlen    length of each on bit in pattern ipat
c-----
        if(ht.le.0.0)then
            hht = 0.07
        else
            hht = ht
        endif
        if(nocx.gt.0)then
            xmin = 10.**x1
            xmax = 10.**(x1+nocx)
        endif
        if(nocy.gt.0)then
            ymin = 10.**y1
            ymax = 10.**(y1+nocy)
        endif
c-----
c       plot, checking that bounds are not exceeded in log plots
c-----
        ipen = 3
        do 100 i=1,n
            if(nocx.gt.0)then
                if(x(i).ge.xmin .and. x(i).le.xmax)then
                    xx = (alog10(x(i))-x1)*deltax
                elseif(x(i).lt.xmin)then
                    xx = 0.0
                else
                    xx =  nocx * deltax
                endif
            else
                xx = (x(i)-x1)/deltax
            endif
            if(nocy.gt.0)then
                if(y(i).ge.ymin .and. y(i).le.ymax)then
                    yy = (alog10(y(i))-y1)*deltay
                else if(y(i).lt.ymin)then
                    yy = 0.0
                else
                    yy =   deltay *nocy
                endif
            else
                yy = (y(i)-y1)/deltay
            endif
            if(ipat.eq.0 .or. ipen.eq.3)then
                call plot(xx,yy,ipen)
            else
                call plotd(xx,yy,ipat,xlen)
            endif
            if(lintyp.ne.0)call symbol(xx,yy,hht,char(inteq),0.0,-1)
            if(lintyp.ge.0)ipen=2
  100   continue
        return
        end
