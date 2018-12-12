c---------------------------------------------------------------------c
c                                                                     c
c      COMPUTER PROGRAMS IN SEISMOLOGY                                c
c      VOLUME I                                                       c
c                                                                     c
c      PROGRAM: GWIDTH                                                c
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
        subroutine gwidth(width) 
        common/Gcplot/owidth
        common/Ycplot/xstp,ystp
        integer*4 ix
            owidth = width
            xx = 1000.0*xstp*owidth
            if(xx.gt.1000000000.0)then
                ix = 1000000000
            else
                ix = xx
            endif
            if(ix .lt.0)ix = 0
            call dfgwid(ix)
        return
        end
