c---------------------------------------------------------------------c
c                                                                     c
c      COMPUTER PROGRAMS IN SEISMOLOGY                                c
c      VOLUME I                                                       c
c                                                                     c
c      PROGRAM: FRAME                                                 c
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
        subroutine frame()
c-----
c       start new plot on new page
c       move current position to lower left corner
c       define new origin
c       do not redefine factor scaling, thougn
c-----
        common/Ycplot/xstp,ystp
        common/Xcplot/xold,yold,xcur,ycur
        integer*4 mode
        call plot(-xold/xstp,-yold/ystp,-3)
        xcur=0.0
        ycur=0.0
        mode = 0
        call dferas(mode)
        return
        end

        subroutine gframe(mode)
c-----
c       start new plot on new page
c       move current position to lower left corner
c       define new origin
c       do not redefine factor scaling, thougn
c
c       mode    1   automatically go to next page without 
c               user interaction under XVIG
c       mode    0   wait for user Next from window
c-----
        common/Ycplot/xstp,ystp
        common/Xcplot/xold,yold,xcur,ycur
        integer*4 pode
        call plot(-xold/xstp,-yold/ystp,-3)
        xcur=0.0
        ycur=0.0
        pode = mode
        call dferas(pode)
        return
        end
