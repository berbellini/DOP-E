c---------------------------------------------------------------------c
c                                                                     c
c      COMPUTER PROGRAMS IN SEISMOLOGY                                c
c      VOLUME I                                                       c
c                                                                     c
c      PROGRAM: PLTSCL                                                c
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
        subroutine pltscl(x,axlen,n,x1,deltax,nocx)
c-----
c       PLTSCL
c           GIVEN ARRAY, ESTABLISHES SCALING PARAMETERS
c           FOR LOG scale, sets according to largest value
c       
c           x   -   array of values
c           axlen   -   length of axis in inches
c           n   -   number of points
c           x1  -   value of first point of axis
c                   returned value
c           deltax  -   inches per cycle if log
c                   units per inch is linear
c           nocx    -   positive - number of cycles on
c                   log scale
c
c                   zero or negative - linear scale
c-----
        real*4 x(1)
        if(nocx.gt.0)then
            xmax = -1.0e+38
            do 100 i=1,n
                if(x(i).gt.xmax)xmax = x(i)
  100       continue
            if(xmax.gt.0.0)then
                xmax = alog10(xmax)
                if(xmax.gt.0.0)then
                    lmax = xmax +1
                else
                    lmax = xmax
                endif
                x1 = lmax
            else
                x1 = 0.0
            endif
            deltax = axlen/nocx
            x1 = x1 - nocx
        else
            call gscale(x,axlen,n,1)
            x1 = x(n+1)
            deltax = x(n+2)
        endif
        return
        end
