c---------------------------------------------------------------------c
c                                                                     c
c      COMPUTER PROGRAMS IN SEISMOLOGY                                c
c      VOLUME I                                                       c
c                                                                     c
c      PROGRAM: SHADET                                                c
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
        subroutine shadet(x1,y1,x2,y2,x3,y3,ipatx,ipaty,xlen,ylen)
c-----
c       shade a triangular region
c-----
c       x1,y1   - coordinate of one corner of triangle
c       x2,y2   - coordinate of second corner of triangle
c       x3,y3   - coordinate of third corner of trianlge
c             (one or more corners may be identical, e.g.,
c              degenerates to a line or a point )
c       ipatx   - pattern for shading used in drawing lines 
c             parallel to x-axis
c       ipaty   - pattern for shading used in drawing lines
c             parallel to y-axis
c             (patterns are integer from 0 to 31)
c       xlen    - length in inches of one bit of pattern
c       ylen    - length in inches of one bit of pattern for y-axis
c             (see discussion of plotd, etc)
c-----
c-----
      common/Xcplot/xold,yold,xcur,ycur
        common/Ycplot/xstp,ystp
c-----
c      xold and yold are absolute plotter coordinates
c      xcur and ycur are coordinates of point relative
c            to current plotter position in current
c            software units
c-----
        integer*4 jx1,jy1,jx2,jy2,jx3,jy3,jpatx,jpaty,jlenx,jleny
        real*4 arr(6)
        integer*4 iarr(6)
        integer*2 jarr(6)
        arr(1) = x1
        arr(2) = y1
        arr(3) = x2
        arr(4) = y2
        arr(5) = x3
        arr(6) = y3
c-----
c current position is in inches, ala calcomp
c but to use unix plot filters we need integer*2
c the factor 1000 is really 11000 counts / 11.0 inches
c-----
        do 100 i=1,6
            if(i.eq.1 .or. i.eq.3 .or. i.eq.5)then
                arr(i) = arr(i) * xstp
                arr(i) = 1000.0*(arr(i) + xold)
                iarr(i) = arr(i)
            else
                arr(i) = arr(i) * ystp
                arr(i) = 1000.0*(arr(i) + yold)
                iarr(i) = arr(i)
            endif
            if(iarr(i) .gt. 1000000000)iarr(i)=1000000000
            if(iarr(i) .lt.-1000000000)iarr(i)=-1000000000
                jarr(i) = iarr(i)
  100   continue
        jx1 = jarr(1)
        jy1 = jarr(2)
        jx2 = jarr(3)
        jy2 = jarr(4)
        jx3 = jarr(5)
        jy3 = jarr(6)
        jpatx = ipatx
        jpaty = ipaty
        jlenx = ( 1000.0 * xstp * xlen)
        jleny = ( 1000.0 * ystp * ylen)
c-----
c       protect against later divide errors
c-----
        if(jlenx .le. 0)jlenx = 1
        if(jleny .le. 0)jleny = 1
        call dffilt(jx1,jy1,jx2,jy2,jx3,jy3,jpatx,jpaty,jlenx,jleny)
        return
        end
