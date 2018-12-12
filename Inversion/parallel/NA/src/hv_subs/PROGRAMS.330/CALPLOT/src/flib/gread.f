c---------------------------------------------------------------------c
c                                                                     c
c      COMPUTER PROGRAMS IN SEISMOLOGY                                c
c      VOLUME I                                                       c
c                                                                     c
c      ROUTINE: GREAD                                               c
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
c
c       Read the PageNum'th page from the file fname.
c       Clip the plot page about the limits (Xlow,Ylow) -> (Xhigh,Yhigh)
c       The origin is redefined to be the lower left corner of the
c       clip region. This corner is then shifted to
c       position (X0, Y0). Then apply a scale factor sclfac
c
c
        subroutine gread(fname, X0, Y0, Xlow, Ylow, Xhigh, Yhigh, 
     1      PageNum, tsclx, tscly)
        character fname*(*)
        real X0, Y0, Xlow, Ylow, Xhigh, Yhigh
        integer PageNum
        real tsclx, tscly
        integer*4 lstr
        integer*4 NumX, NumY, LowX, LowY, HighX, HighY, Num
        integer*4 Sclx, Scly
        Num = PageNum
        Sclx = (1000.0 * tsclx)
        Scly = (1000.0 * tscly)
        call doconv(X0-Xlow*tsclx, Y0-Ylow*tscly, NumX, NumY)
        call doconv(Xlow, Ylow, LowX, LowY)
        call doconv(Xhigh, Yhigh, HighX, HighY)
        lstr = len(fname)
        call dfread(lstr,fname, NumX, NumY, 
     1      LowX, LowY, HighX, HighY, Num, Sclx, Scly)
        return
        end

        subroutine  doconv(xn, yn, ix, iy)
        real xn, yn
        integer*4 ix, iy
        real xx, yy
        common/Xcplot/xold,yold,xcur,ycur
            xx = 1000. * (xn + xold)
            yy = 1000. * (yn + yold)
            if(xx .gt. 1000000000.0)then
                ix = 1000000000
            else if(xx .lt. -1000000000.0)then
                ix = -1000000000
            else
                ix = xx 
            endif
            if(yy .gt. 1000000000.0)then
                iy = 1000000000
            else if(yy .lt. -1000000000.0)then
                iy = -1000000000
            else
                iy = yy 
            endif
        return
        end
