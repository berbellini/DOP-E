c---------------------------------------------------------------------c
c                                                                     c
c      COMPUTER PROGRAMS IN SEISMOLOGY                                c
c      VOLUME I                                                       c
c                                                                     c
c      PROGRAM: GCONT                                                 c
c                                                                     c
c      COPYRIGHT (C)  2004 R. B. Herrmann                             c
c                                                                     c
c      Department of Earth and Atmospheric Sciences                   c
c      Saint Louis University                                         c
c      3507 Laclede Avenue                                            c
c      St. Louis, Missouri 63103                                      c
c      U. S. A.                                                       c
c                                                                     c
c---------------------------------------------------------------------c
        subroutine gcont(x,y)
c-----
c       basic routine to move point on plotter
c
c       x       abscissa in inches
c       y       ordinate in inches
c-----
        common/Gcplot/owidth
        common/Xcplot/xold,yold,xcur,ycur
        common/Ycplot/xstp,ystp
        common/Zcplot/iunit,ifont
c-----
c      xold and yold are absolute plotter coordinates
c      xcur and ycur are coordinates of point relative
c            to current plotter position in current
c            software units
c-----
        integer*4 ix,iy
c-----
c       x and y are page coordinates in units given by
c       last call gunit()
c-----
        xcur=x
        ycur=y
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
        call dfcont(ix,iy)
      return
      end
