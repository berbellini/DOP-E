c---------------------------------------------------------------------c
c                                                                     c
c      COMPUTER PROGRAMS IN SEISMOLOGY                                c
c      VOLUME I                                                       c
c                                                                     c
c      PROGRAM: SHADEP                                                c
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
        subroutine shadep(n, xarr, yarr)
c-----
c       shade a polygon
c-----
c       n       - number of x,y coordinate pairs
c       xarr    - x-array of abscissa
c       yarr    - y-array of ordinates
c-----
        common/Xcplot/xold,yold,xcur,ycur
        common/Ycplot/xstp,ystp
c-----
c      xold and yold are absolute plotter coordinates
c      xcur and ycur are coordinates of point relative
c            to current plotter position in current
c            software units
c-----
        integer*4 n
        real*4 xarr(n), yarr(n)
c-----
c       since we must convert to integers for calplot and
c       since we need to alloc() space for the integer arrays
c       we let the C routine do all the work, and 
c       pass all necessary information to the C routine
c-----
        call dffilp(n,xarr,yarr,xold,yold,xcur,ycur,xstp,ystp)
        return
        end
