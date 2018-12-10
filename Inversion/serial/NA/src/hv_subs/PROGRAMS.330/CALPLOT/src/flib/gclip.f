c---------------------------------------------------------------------c
c                                                                     c
c      COMPUTER PROGRAMS IN SEISMOLOGY                                c
c      VOLUME I                                                       c
c                                                                     c
c      PROGRAM: CLIP                                                  c
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
        subroutine  gclip(cmd,  xlow, ylow, xhigh, yhigh )
        character cmd*(*)
        real*4 xlow, ylow, xhigh, yhigh
c-----
c       Set clip region
c       cmd Char*(*)    "ON" or "OFF"
c       xlow    R*4     Coordinates of bounding box
c       ylow    R*4     
c       xhigh   R*4     
c       yhigh   R*4     
c-----
        common/Xcplot/xold,yold,xcur,ycur
        common/Ycplot/xstp,ystp

        integer*4 ilw, jlw, iup, jup
        integer*4 icmd
        xx = 1000. * (xlow*xstp + xold)
        yy = 1000. * (ylow*xstp + yold)
        if(xx .gt. 1000000000.0)then
            ilw = 1000000000
        else if(xx .lt. -1000000000.0)then
            ilw = -1000000000
        else
            ilw = xx 
        endif
        if(yy .gt. 1000000000.0)then
            jlw = 1000000000
        else if(yy .lt. -1000000000.0)then
            jlw = -1000000000
        else
            jlw = yy 
        endif
        xx = 1000. * (xhigh*xstp + xold)
        yy = 1000. * (yhigh*xstp + yold)
        if(xx .gt. 1000000000.0)then
            iup = 1000000000
        else if(xx .lt. -1000000000.0)then
            iup = -1000000000
        else
            iup = xx 
        endif
        if(yy .gt. 1000000000.0)then
            jup = 1000000000
        else if(yy .lt. -1000000000.0)then
            jup = -1000000000
        else
            jup = yy 
        endif
        if(cmd(1:1) .eq. 'o' .or. cmd(1:1) .eq. 'O')then
            if(cmd(2:2) .eq. 'n' .or. cmd(2:2) .eq. 'N')then
                icmd = 1
            else
                icmd = 0
            endif
        else
            icmd = 0
        endif
        call dfclip(icmd, ilw, jlw, iup, jup)
        return
        end
