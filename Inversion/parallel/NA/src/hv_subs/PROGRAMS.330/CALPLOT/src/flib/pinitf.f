c---------------------------------------------------------------------c
c                                                                     c
c      COMPUTER PROGRAMS IN SEISMOLOGY                                c
c      VOLUME I                                                       c
c                                                                     c
c      PROGRAM: PINITF                                                c
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
        subroutine  pinitf(string)
        character string*(*)
        character constr*2
            constr = ' '
            call  ginitf(string, constr)
        return
        end

        subroutine  ginitf(string, constr)
        character string*(*), constr*(*)
        common/Xcplot/xold,yold,xcur,ycur
        common/Scplot/xs,ys
        common/Zcplot/iunit,ifont
        integer*4 x0,y0,x1,y1
        integer*4 zero, mone
        integer*4 lstr, lcon
c       data x0/0/,y0/0/,x1/7620/,y1/7620/
        x0 = 0
        y0 = 0
        x1 = 7620
        y1 = 7620
        xs = 0.0
        ys = 0.0
        xold=0.0
        yold=0.0
        xcur=0.0
        ycur=0.0
        iunit = 0
        ifont = 0
        mone = -1
        zero = 0
        lstr = len(string)
        lcon = len(constr)
        call dfopen(string,constr,lstr,lcon)
c-----
c map 11" to max tek screen
c this is really done in plot where
c we let 11" = 11000 points
c-----
        call dfspce(x0,y0,x1,y1)
        call dfclip(zero, mone, mone, mone, mone)
        call gwidth(0.0)
        call factor(1.0)
        call plot(0.0,0.0,-3)
        return
        end
