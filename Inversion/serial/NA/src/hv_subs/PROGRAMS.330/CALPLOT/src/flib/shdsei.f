c---------------------------------------------------------------------c
c                                                                     c
c      COMPUTER PROGRAMS IN SEISMOLOGY                                c
c      VOLUME I                                                       c
c                                                                     c
c      PROGRAM: SHDSEI                                                c
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
      subroutine shdsei(x1,y1,ixy,istnd,iplmn)
c-----
c       permit shading of seismic traces
c-----
c       x1  - x coordinate of first point local zero
c       y1  - y coordinate of first point local zero
c       ixy - time axis is parallel to x-axis (0)
c                                      y-axis (1)
c       istnd   - 0 turn off shading  
c             1 turn on  shading
c       iplmn   - 0 shade positive amplitudes
c           - 1 shade negative amplitudes
c-----
c       For example, if time is parallel to x-axis
c       and iplmn = 0, then values of y > y1 will be shaded from y to y1
c-----
        common/Xcplot/xold,yold,xcur,ycur
        common/Ycplot/xstp,ystp
c-----
c      xold and yold are absolute plotter coordinates
c      xcur and ycur are coordinates of point relative
c            to current plotter position in current
c            software units
c-----
        integer*4 jx1, jy1, jxy, jstnd, jplmn
        xcur = x1
        ycur = y1
        xn=xcur*xstp
        yn=ycur*ystp
c-----
c current position is in inches, ala calcomp
c but to use unix plot filters we need integer*2
c the factor 1000 is really 11000 counts / 11.0 inches
c-----
        xx = 1000. * (xn + xold)
        yy = 1000. * (yn + yold)
        if(xx.gt. 1000000000.0)then
            ix = 1000000000
        else if(xx.lt.-1000000000.0)then
            ix = -1000000000
        else
            ix = xx
        endif
        if(yy.gt. 1000000000.0)then
            iy = 1000000000
        else if(yy.lt.-1000000000.0)then
            iy = -1000000000
        else
            iy = yy
        endif
        jx1 = ix
        jy1 = iy
        jxy = ixy
        jstnd = istnd
        jplmn = iplmn
        call dffils(jx1,jy1,jxy,jstnd,jplmn)
        return
        end
