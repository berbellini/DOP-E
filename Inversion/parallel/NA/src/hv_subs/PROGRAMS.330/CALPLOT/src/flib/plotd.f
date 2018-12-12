c---------------------------------------------------------------------c
c                                                                     c
c      COMPUTER PROGRAMS IN SEISMOLOGY                                c
c      VOLUME I                                                       c
c                                                                     c
c      PROGRAM: PLOTD                                                 c
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
        subroutine plotd(xx,yy,jpat,xlen)
c-----
c       draw a line between the current point
c       and the coordinates (xx,yy) using
c       the pattern ipat. ipat is considered
c       to be a binary pattern with 5 bits
c       set to either 0 or 1. The length of
c       each bit is xlen. The 0 bit indicates
c       a dark vector and the 1 bit a plotted vector
c-----
        save is,xsv,ipatld
        data is,xsv,ipatld/0,0.0,-1/
c-----
c       is  - index to bit in pattern ipat
c       xsv - offset within current bit in xlen units
c       ipatld  - previous pattern. If pattern changes reinitialize
c-----
        ipat = jpat
        if(ipat.ne.ipatld)then
            is = 0
            xsv = 0.0
            if(ipat.gt.31 .or. ipat.lt.0)ipat = 0
            ipatld = ipat
        endif
c-----
c       test for solid line
c-----
        if(ipat.eq.0)call plot(xx,yy,2)
        if(ipat.eq.0)return
c-----
c       get current pen position
c-----
        call where(x0,y0,fct)
c-----
c       determine length of segment to be plotted, is zero return
c-----
        s = sqrt( (xx-x0)**2 + (yy-y0)**2 )
        if(s.le.0.001)return
c-----
c       determine direction cosines for parametric representation 
c       of line also remember that the ultimate precision 
c       is 0.001 units in plot space
c-----
        if(abs(xx-x0).lt.0.001)then
            sx = 0.0
        else
            sx = (xx-x0)/s
        endif
        if(abs(yy-y0).lt.0.001)then
            sy = 0.0
        else
            sy = (yy-y0)/s
        endif
        xs = x0
        ys = y0
c-----
c       plot segments
c-----
 1000   continue
        call segplt(s,sx,sy,xlen,xsv,is,xs,ys,ipat,iret)
        if(iret.eq.1)goto 1000
c-----
c       exit subroutine
c-----
        return
        end

        subroutine segplt(s,sx,sy,xlen,xsv,is,xs,ys,ipat,iret)
c-----
c       plot bit pattern segments, either partial or complete
c-----
c       s   - length of segment yet to be plotted
c       sx  - direction cosine in x-direction
c       sy  - direction cosine in y-direction
c       xlen    - length of bit segment
c       xsv - fraction of current segment drawn
c       is  - index to bit in pattern ipat
c       xs  - current x position
c       ys  - current y position
c       ipat    - bit pattern for plot
c       iret    - -1 partial segment drawn
c       iret    - +1 complete segment drawn
c-----
        dsmx = (1.0-xsv)*xlen
        dsx = dsmx*sx
        dsy = dsmx*sy
        if(s .lt. dsmx)then
            xs = xs + s*sx
            ys = ys + s*sy
            xsv = xsv + s/xlen
            if(xsv.ge.1.0)then
                xsv=0.0
                is = is + 1
                if(is.gt.5)is = 0
            endif
            iret = -1
        elseif(s .ge. dsmx) then
            xs = xs + dsx
            ys = ys + dsy
            xsv = 0.0
            iret = +1
            if(dsmx.gt.0.0)then
                s = s -dsmx
            else
                s = s - xlen
            endif
            if(s.lt.0.0)s = 0.0
        endif
            if(mod(ipat/2**is , 2).eq.1)then
                ipen = 2
            else
                ipen = 3
            endif
            call plot(xs,ys,ipen)
            if(iret.eq.+1)then
                is = is + 1
                if(is.gt.5)is=0
            endif
        return
        end
