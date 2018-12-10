c---------------------------------------------------------------------c
c                                                                     c
c      COMPUTER PROGRAMS IN SEISMOLOGY                                c
c      VOLUME I                                                       c
c                                                                     c
c      PROGRAM: ALGAXE                                                c
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
      subroutine algaxe(xaxlen,yaxlen,nocx,nocy,ttlx,ttly,mtx,mty,x1,
     1  y1,deltax,deltay)
c-----
c       ALGAXE
c           PLOTS AXES
c           xaxlen  -   length of x axis in inches
c                   (if <= 0.0 do not plot )
c           yaxlen  -   length of y axis in inches
c                   (if <= 0.0 do not plot )
c           nocx    -   number of logarithmic cycles
c                   on x axis if log plot
c                   = 0 if linear plot
c           nocy    -   number of logarithmic cycles
c                   = 0 if linear plot
c           ttlx    -   x axis title a character string
c           ttly    -   y axis title a character string
c           mtx -   number of characters in x axis
c                   title
c           mty -   number of characters in y ayis
c                   title
c           x1  -   first value in x-axis
c           y1  -   first value in y-axis
c                   if log plot these are the powers of 10
c           deltax  -   if linear number of x-units per inch
c                   if log length of cycle
c           deltay  -   yaxis
c                   on y axis if log plot
c-----
      character*(*) ttlx,ttly
      slt = 0.02*yaxlen
      sst = 0.01 * yaxlen
      sp = -0.06*yaxlen
      ss = 0.035*yaxlen
      ssp = sp + ss - 0.06
      ttlp = -0.11*yaxlen - 0.1
      sttl = 0.035*yaxlen
      xnum = 1
      yl = y1
      yu = y1 + nocy
      if(abs(yl).ge.10. .or. abs(yu).ge.10. )xnum = xnum + 1.
      if(abs(yl).ge.100..or. abs(yu).ge.100.)xnum = xnum + 1.
      if(y1.lt.0) xnum = xnum + 1.0
      xpo = x1
      ypo = y1
c-----
c       draw x-axis
c-----
        if(xaxlen .gt.0.0)then
            if(nocx.eq.0) then
                call axis(0.0,0.0,ttlx,-mtx,xaxlen,0.0,x1,deltax)
            else
                call plot(0.0,-slt,3)
                call plot(0.0,0.0,2)
                anocx = nocx
                factx = xaxlen/anocx
                call symbol(-.6*ss,sp,ss,'10',0.0,2)
                call number(999.,ssp,0.5*ss,x1,0.0,-1)
                call plot(0.0,0.0,3)
                do 3 j = 1,nocx
                    do 2 i = 1,10
                        x = i
                        x = alog10(x) *factx + (j-1)*factx
                        if(i.ne.1)then
                            call plot(x,0.0,2)
                            call plot(x,-sst,2)
                        endif
                        call plot(x,0.0,3)
    2               continue
                    call plot(x,-slt,2)
                    call symbol(x-.6*ss,sp,ss,'10',0.0,2)
                    xpo = xpo + 1.0
                    call number(999.,ssp,0.5*ss,xpo,0.0,-1)
                    call plot(x,0.0,3)
    3           continue
                xtl = mtx
                xtl = (xaxlen-xtl*sttl)/2.0
                call symbol(xtl,ttlp,sttl,ttlx,0.0,mtx)
            endif
            call plot(0.0,0.0,3)
        endif
c-----
c       draw y-axis
c-----
        if(yaxlen .gt. 0.0)then
            if(nocy.eq.0) then
                call axis(0.0,0.0,ttly,mty,yaxlen,90.,y1,deltay)
            else
                call plot(-slt,0.0,3)
                call plot(0.0,0.0,2)
                anocy = nocy
                sp = sp - (xnum - 1.5) * 0.5 * ss
                ttlp = ttlp - (xnum-1.)*0.5*ss
                facty = yaxlen/anocy
                call symbol(sp-0.4,-0.5*ss,ss,'10',0.0,2)
                call number(999.,.5*ss-.06,.5*ss,y1,0.0,-1)
                call plot(0.0,0.0,3)
                do 9 j = 1,nocy
                    do 8 i = 1,10
                        y = i
                        y = alog10(y) * facty + (j-1)*facty
                        if(i.ne.1)then
                            call plot(0.0,y,2)
                            call plot(-sst,y,2)
                        endif
                        call plot(0.0,y,3)
    8               continue
                    call plot(-slt,y,2)
                    call symbol(sp-.4,y-.5*ss,ss,'10',0.0,2)
                    ypo = ypo + 1
                    call number(999.,y+.5*ss-.06,.5*ss,ypo,0.0,-1)
                    call plot(0.0,y,3)
    9           continue
                ytl=mty
                ytl = (yaxlen-ytl*sttl)/2.0
                call symbol(ttlp-.2,ytl,sttl,ttly,90.,mty)
            endif
        endif
        return
        end

