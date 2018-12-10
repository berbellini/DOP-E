c-----
c       changes
c
c       01 JAN 2004 adjusted position of Y-axis title in doliny dolny 
c           for case when numbers are negative
c----
        subroutine gbox(xl,yl,xh,yh)
        call plot(xl,yl,3)
        call plot(xh,yl,2)
        call plot(xh,yh,2)
        call plot(xl,yh,2)
        call plot(xl,yl,2)
        return
        end

        subroutine gcent(xx,yy,ht,string,angle)
            character string*(*)
            il = lgstr(string)
            t = angle*3.1415927/180.0
            ct = cos(t)
            st = sin(t)
            rl = 0.5*il*ht
            
            xxx = xx - rl*ct
            yyy = yy - rl*st
            call symbol(xxx,yyy,ht,string,angle,il)
        return
        end

        subroutine gright(xx,yy,ht,string,angle)
            character string*(*)
            il = lgstr(string)
            t = angle*3.1415927/180.0
            ct = cos(t)
            st = sin(t)
            rl = 1.0*il*ht
            
            xxx = xx - rl*ct
            yyy = yy - rl*st
            call symbol(xxx,yyy,ht,string,angle,il)
        return
        end

        subroutine gleft(xx,yy,ht,string,angle)
            character string*(*)
            il = lgstr(string)
            t = angle*3.1415927/180.0
            ct = cos(t)
            st = sin(t)
            rl = 0.0*il*ht
            
            xxx = xx - rl*ct
            yyy = yy - rl*st
            call symbol(xxx,yyy,ht,string,angle,il)
        return
        end

        subroutine dology(x0,y0,yleng,symax,symin,
     1      sizey,ticlft,lablft,dopow,ly,titley)
c-----
c       put in tic marks along the y-axis 
c
c       x0  R*4 - position of bottom side of axis
c       y0  R*4 - position of bottom side of axis
c       yleng   R*4 - length of Y axis
c       ymax    R*4 - maximum value of number corresponding to
c                   top
c       ymin    R*4 - minimum value of number corresponding to
c                   bottom
c       nocy    I*4 - number of cycles along the axis
c       sizey   R*4 - size of numbers, e.g., 10
c       ticlft  L   - .true. tics go left
c                 .false. tics go right
c       lablft  L   - .true. 10**N goes to left of axis
c                 .false. goes to right
c       dopow   L   - .true. put in 10**N labels
c       ly  I*4 length of y-axis title string
c       titley  Ch  y-axis title string
c-----
        real*4 x0, y0, yleng, ymax
        integer*4 nocy
        logical ticlft, lablft, dopow
        integer ly
        character titley*(*)
        character title*80

        real*4 xnum(9)
        data xnum/1.0 ,2.0 ,3.0 ,4.0 ,5.0 ,6.0 ,7.0 ,8.0 ,9.0/

        ymin = symin
        ymax = symax
c-----
c       put in the axis
c-----
        ycen = y0 + 0.5*yleng
        call plot(x0,y0,3)
        call plot(x0,y0+yleng,2)
c-----
c       set up positioning for 10**N labels
c-----
c-----
c       compute limits for scale
c-----
        ymxlog = alog10(ymax)
        ymmin = alog10(ymin)
        nocy = ymxlog - ymmin + 1
        iy = ymmin
        if(ymmin .lt. 0)iy = iy - 1

        ymlog = alog10(ymax)
        iiy = ymmin
        if(ymlog .lt. 0)iiy = iiy - 1

        ja = max(abs(iy),abs(iiy))
        if(ja.lt.10)then
            jpow = 1
        else if(ja.ge.10 .and. ja.lt.100)then
            jpow = 2
        else if(ja.ge.100 .and. ja.lt.1000)then
            jpow = 3
        endif
        ii = min(iy,iiy)
        if(ii.lt.0)jpow = jpow + 1
        jpow = jpow + 2
        xshft = (2 + jpow*0.7)*sizey

        if(lablft)then
            if(ticlft)then
                xpos = x0 - 0.50*sizey - xshft
                xpos2 = x0 - 0.50*sizey - xshft -1.2*sizey
                xpos1 = xpos + 0.7*sizey
            else
                xpos = x0 + 0.50*sizey - xshft
                xpos2 = x0 + 0.50*sizey - xshft -1.2*sizey
                xpos1 = xpos + 0.7*sizey
            endif
        else
            if(ticlft)then
                xpos = x0 + 1.00*sizey
                xpos2 = x0 + 0.00*sizey + 0.2*sizey + xshft
                xpos1 = xpos + 0.7*sizey
            else
                xpos = x0 + 2.0*sizey
                xpos2 = x0 + 1.0*sizey + 0.2*sizey + xshft
                xpos1 = xpos + 0.7*sizey
            endif
        endif

        ifirst = 1
        ylow = y0
        yscal = yleng/alog10(ymax/ymin)
        do 100 ii=iy,iy+nocy+2
            tenpow = 10.0**ii
            do 200 jj=1,9
                yval = xnum(jj)*tenpow
c-----
c               put in tics where appropriate
c-----
                if(yval .ge. ymin .and. yval.le.ymax)then
                    yv = alog10(yval/ymin)
                    yy = y0 + yv*yscal
                    if(jj.eq.1)then
                        ticlen = 1.5*sizey
c-----
c                       put in 10**N
c-----
                        if(ifirst.eq.1)then
                            ypos = yy
                        else
                            ypos = yy - 0.5*sizey
                    if(yval.eq.ymax)ypos=yy-1.4*sizey
                            if(ypos.lt.y0)ypos = y0
                        endif
                        if(dopow)then
                        ypos1 = ypos + 0.7*sizey
            call symbol(xpos,ypos,sizey,'10',0.0,2)
            call number(999.,ypos1,0.7*sizey,real(ii),0.0,-1)
                        endif
                    else
                        ticlen = sizey
                    endif
                        ifirst = ifirst + 1
                    call plot(x0,yy,3)
                    if(ticlft)then  
                        call plot(x0-ticlen,yy,2)
                    else
                        call plot(x0+ticlen,yy,2)
                    endif
                endif
  200       continue
            call plot(x0,yy,3)
  100   continue
c-----
c       put in the title if we put in the numbers
c-----
        if(dopow .and. ly.gt.0)then
            title = ' '
            title(1:ly) = titley(1:ly)
            if(lablft)then
                call gcent(xpos2,ycen,1.2*sizey,title, 90.0)
            else
                call gcent(xpos2,ycen,1.2*sizey,title,-90.0)
            endif
        endif
        return
        end

        subroutine dologx(x0,y0,xleng,sxmax,sxmin,
     1      sizex,ticup,labtop,dopow,lx,titlex)
c-----
c       put in tic marks along the x-axis 
c
c       x0  R*4 - position of left side of axis
c       y0  R*4 - position of left side of axis
c       xleng   R*4 - length of X axis
c       xmax    R*4 - maximum value of number corresponding to
c                   far right
c       xmin    R*4 - maximum value of number corresponding to
c                   far left
c       sizex   R*4 - size of numbers, e.g., 10
c       ticup   L   - .true. tics go up
c                 .false. tics go down
c       labtop  L   - .true. 10**N goes above axis
c                 .false. goes below
c       dopow   L   - .true. put in 10**N labels
c       lx  I*4 length of x-axis title string
c       titlex  Ch  x-axis title string
c-----
        real*4 x0, y0, xleng, xmax
        integer*4 nocx
        logical ticup, labtop, dopow
        integer lx
        character titlex*(*)
        character title*80

        real*4 xnum(9)
        data xnum/1.0 ,2.0 ,3.0 ,4.0 ,5.0 ,6.0 ,7.0 ,8.0 ,9.0/

        xmin = sxmin
        xmax = sxmax
c-----
c       put in the axis
c-----
        xcen = x0 + 0.5*xleng
        call plot(x0,y0,3)
        call plot(x0+xleng,y0,2)
c-----
c       set up positioning for 10**N labels
c-----
        if(labtop)then
            if(ticup)then
                ypos  = y0   + 2.0 * sizex
                ypos2 = y0   + 3.7 * sizex
                ypos1 = ypos + 0.7*sizex
            else
                ypos  = y0   + sizex
                ypos2 = y0   + 2.7*sizex
                ypos1 = ypos + 0.7*sizex
            endif
        else
            if(ticup)then
                ypos  = y0   - 2.0*sizex
                ypos2 = y0   - 3.7*sizex
                ypos1 = ypos + 0.7*sizex
            else
                ypos  = y0   - 3.0*sizex
                ypos2 = y0   - 4.7*sizex
                ypos1 = ypos + 0.7*sizex
            endif
        endif
c-----
c       compute limits for scale
c-----
C       xmin = xmax / 10.0**nocx
        xmxlog = alog10(xmax)
        xmmin = alog10(xmin)
        nocx = xmxlog - xmmin + 1
        ix = xmmin
        if(xmmin .lt. 0)ix = ix - 1
        ifirst = 1
        do 100 ii=ix,ix+nocx+2
            tenpow = 10.0**ii
            ja = abs(ii)
            if(ja.lt.10)then
                jpow = 1
            else if(ja.ge.10 .and. ja.lt.100)then
                jpow = 2
            else if(ja.ge.100 .and. ja.lt.1000)then
                jpow = 3
            endif
            if(ii.lt.0)jpow = jpow + 1
            xscal = xleng/alog10(xmax/xmin)
            do 200 jj=1,9
                xval = xnum(jj)*tenpow
c-----
c               put in tics where appropriate
c-----
                if(xval .ge. xmin .and. xval.le.xmax)then
                    xv = alog10(xval/xmin)
                    xx = x0 + xv*xscal
                    if(jj.eq.1)then
                        ticlen = 1.5*sizex
c-----
c                       put in 10**N
c-----
                    if(ifirst.eq.1)then
                            xpos = xx
                    else
                        xshft = (2 + jpow*0.7)*sizex
                        if(xx+xshft .gt. x0+xleng)then
                            xpos = x0+xleng - xshft
                        else
                            xpos = xx - xshft/2.0
                        endif
                        if(xpos.lt.x0)xpos = x0
                    endif
                    if(dopow)then
                    call symbol(xpos,ypos,sizex,'10',0.0,2)
            call number(999.,ypos1,0.7*sizex,real(ii),0.0,-1)
                        endif
                    else
                        ticlen = sizex
                    endif
                        ifirst = ifirst + 1
                    call plot(xx,y0,3)
                    if(ticup)then   
                        call plot(xx,y0+ticlen,2)
                    else
                        call plot(xx,y0-ticlen,2)
                    endif
                endif
  200       continue
            call plot(xx,y0,3)
  100   continue
        if(dopow .and. lx.gt.0)then
            title = ' '
            title(1:lx) = titlex(1:lx)
            call gcent(xcen,ypos2,1.2*sizex,title,0.0)
        endif
        return
        end

        subroutine dnlinx(x0,y0,xleng,xmax,xmin,
     1      sizex,ticup,labtop,dopow,llx,titlex)
c-----
c       put in tic marks along the x-axis 
c
c       x0  R*4 - position of left side of axis
c       y0  R*4 - position of left side of axis
c       xleng   R*4 - length of X axis
c       xmax    R*4 - maximum value of number corresponding to
c                   far right
c       xmin    R*4 - maximum value of number corresponding to
c                   far left
c       sizex   R*4 - size of numbers, e.g., 10
c       ticup   L   - .true. tics go up
c                 .false. tics go down
c       labtop  L   - .true. number goes above axis
c                 .false. goes below
c       dopow   L   - .true. put in number labels
c       llx I*4 length of x-axis title string
c       titlex  Ch  x-axis title string
c-----
        real*4 x0, y0, xleng, xmax, xmin, sizex
        logical ticup, labtop, dopow
        integer llx
        character titlex*(*)
        call dolnx(x0,y0,xleng,xmax,xmin,
     1      sizex,ticup,labtop,dopow,llx,titlex,.false.)
        return
        end

        subroutine dolinx(x0,y0,xleng,xmax,xmin,
     1      sizex,ticup,labtop,dopow,llx,titlex)
c-----
c       put in tic marks along the x-axis 
c
c       x0  R*4 - position of left side of axis
c       y0  R*4 - position of left side of axis
c       xleng   R*4 - length of X axis
c       xmax    R*4 - maximum value of number corresponding to
c                   far right
c       xmin    R*4 - maximum value of number corresponding to
c                   far left
c       sizex   R*4 - size of numbers, e.g., 10
c       ticup   L   - .true. tics go up
c                 .false. tics go down
c       labtop  L   - .true. number goes above axis
c                 .false. goes below
c       dopow   L   - .true. put in number labels
c       llx I*4 length of x-axis title string
c       titlex  Ch  x-axis title string
c-----
        real*4 x0, y0, xleng, xmax, xmin, sizex
        logical ticup, labtop, dopow
        integer llx
        character titlex*(*)
        call dolnx(x0,y0,xleng,xmax,xmin,
     1      sizex,ticup,labtop,dopow,llx,titlex,.true.)
        return
        end

        subroutine dolnx(x0,y0,xleng,xmax,xmin,
     1      sizex,ticup,labtop,dopow,llx,titlex,doasis)
c-----
c       put in tic marks along the x-axis 
c
c       x0  R*4 - position of left side of axis
c       y0  R*4 - position of left side of axis
c       xleng   R*4 - length of X axis
c       xmax    R*4 - maximum value of number corresponding to
c                   far right
c       xmin    R*4 - maximum value of number corresponding to
c                   far left
c       sizex   R*4 - size of numbers, e.g., 10
c       ticup   L   - .true. tics go up
c                 .false. tics go down
c       labtop  L   - .true. number goes above axis
c                 .false. goes below
c       dopow   L   - .true. put in number labels
c       llx I*4 length of x-axis title string
c       titlex  Ch  x-axis title string
c       doasis  L   .true. plot values
c               .false. annotate with negative (useful for depth plot)
c-----
        real*4 x0, y0, xleng, xmax, xmin, sizex
        logical ticup, labtop, dopow
        integer llx
        character titlex*(*)
        logical doasis
        character title*80

        logical dosci


        lx = llx
        if(lx.lt.0)lx = 0
c-----
c       put in the axis
c-----
        xl = x0 + xleng
        xcen = x0 + 0.5*xleng
        call plot(x0,y0,3)
        call plot(xl,y0,2)
c-----
c       set up positioning for tics
c-----
        tlen1 = 0.6*sizex
        tlen2 = 1.2*sizex
        if(labtop)then
            if(ticup)then
                ypos = y0 + 2.0 * sizex
                ypos1 = ypos + 0.7*sizex
                ypos2 = y0 + 4.0 * sizex
            else
                ypos = y0 + sizex
                ypos1 = ypos + 0.7*sizex
                ypos2 = y0 + 3.0 * sizex
            endif
        else
            if(ticup)then
                ypos = y0 - 2.0*sizex
                ypos1 = ypos + 0.7*sizex
                ypos2 = y0 - 4.5*sizex
            else
                ypos = y0 - 3.0*sizex
                ypos1 = ypos + 0.7*sizex
                ypos2 = y0 - 5.5*sizex
            endif
        endif
c-----
c       compute limits for scale
c-----
        xmn = min(xmin,xmax)
        xmx = max(xmin,xmax)
        xxv = xmn
c-----
c       safety
c-----
        if(xmn.eq.xmx)xmx = xmn + 1.0
c-----
c       get some idea of the size of the number to be plotted
c-----
        xval = max(abs(xmn),abs(xmx))
c-----
c       we set the maximum number of decimal digits as 3, 
c       if the number is greater than 1000 or less than 0.01, use
c       scientific notation
c       get the maximum number of digits for the number
c-----
        if(xval.le.0.01 .and.xval.gt.0.0 )then
            dosci = .true.
            xnorm = alog10(1.0/xval)
            nxnorm = xnorm+1
            xnorm = 10.0**( nxnorm)
            nscalx = -nxnorm
        else if( xval.ge.10000.0)then
            dosci = .true.
            xnorm = alog10(xval)
            nxnorm = xnorm
            nscalx =  nxnorm
            xnorm = 10.0**(-nxnorm)
        else
            dosci = .false.
            xnorm = 1.0
            nxnorm = 0.0
            nscalx = 0
        endif
c-----
c       choose the maximum number of tics somehow on
c       xleng, and size of numbers and necessary decimals
c-----
        dxx = xnorm*(XMX - XMN)/5.0
        if(dxx.eq.0.0)dxx = 1.0
c-----
c       get start limits for search - we will place at most 10 tics
c       here
c-----
        n = alog10(dxx)
        if(dxx.lt.1.0)n = n -1
        ifac = (10.0**(-n)*dxx)
        
c----- 
c       now get a new minimum value that is some multiple 
c       of ifac 10**n 
c----- 
        dxx = ifac * 10.0**(n)
        ilow = xnorm*xmn/dxx -1
        iup  = xnorm*xmx/dxx +1
        if(ifac.eq.1)then
            ndxxx = 5
        else if(ifac.eq.2)then
            ndxxx = 5
        else
            ndxxx = ifac
        endif
        dxxx = dxx/ndxxx

c-----
c       loop for labels
c-----
c       xe is the end of the previous label
c-----
        xc = xmn
        call plot(x0,y0,3)
c-----
c       here
c       xx is the position of the plot corresponding to a value of xxv
c       xb is the position of the beginning of the number ot be plotted
c       xe is the position of the end of the last number 
c           plotted + one space
c       xl is the end of the line
c-----
c-----
c       consider numbers between ilow*dxx and iup*dxx
c-----
        xe = x0
        do 1000 i=ilow,iup
            do 1001 j=0,ndxxx-1
            
            xxv = i*dxx + j*dxxx
            xxv = xxv/xnorm
            if(xxv.ge.xmn .and. xxv .le.xmx)then
c-----
c               We can plot this number/tic
c-----
                xx = x0 + (xxv - xmn)*xleng/(xmx - xmn)
                call plot(xx,y0,3)
                if(j.eq.0 )then
c-----
c       estimate the number of digits for proper labeling
c       only do this once for the plot
c           dosci   = .false.    +0.000
c               = .true.     +0.000 10^+00
c       also Never overlay plot labels
c----
c       get width of number string
c-----
                    if(doasis)then
                        xv = xxv*xnorm
                    else
                        xv = - xxv*xnorm
                    endif
                    nshft = 1
                    if(xv .lt. 0.0)nshft=nshft+1
                    if(dosci)then
                        nshft = nshft + 1
                        ndec = 2
                        nshft = nshft +  ndec
                    else 
                        if(abs(dxx) .gt. 1.0)then
                            ndec = -1
                        else
                            if(ABS(dxx).lt.0.1)then
                                ndec = 3
                            else
                                ndec = 2
                            endif
                            nshft = nshft + ndec
                        endif
                    if(abs(xv) .ge. 10.0) nshft=nshft+1
                    if(abs(xv) .ge. 100.0)nshft=nshft+1
                    if(abs(xv) .ge. 1000.0)nshft=nshft+1
                    endif
                    xxx = 0.5*(nshft -0.5)*sizex
                    xb = xx - xxx
                    if((xb - xe) .gt. 2.*sizex)then
                    if(xb.ge.x0.and.(xx+xxx).le.xl)then
                    if(dopow)then
                    call number(xb,ypos,sizex,xv,0.0,ndec)
                    endif
                    xe = xb + xxx + sizex
                    endif
                        ticlen = tlen2
                    else
                        ticlen = 0.8*tlen2
                    endif
                else
                    ticlen = tlen1
                endif
                    call plot(xx,y0,3)
                    if(ticup)then   
                        call plot(xx,y0+ticlen,2)
                    else
                        call plot(xx,y0-ticlen,2)
                    endif
                    call plot(xx,y0,3)
            endif
 1001       continue
 1000   continue
c-----
c       put in the title if we put in the numbers
c-----
        if(dopow )then
            sizexx = 1.2*sizex
            title = ' '
            if(lx.gt.0)then
                title(1:lx) = titlex(1:lx)
            else
                lx = 0
            endif
            if(dosci)then
                title(lx +1:lx+5)=' *10 '
                call gcent(xcen,ypos2,sizexx,title,0.0)
                xe = xcen + 0.5*sizexx*(lx+4)
                call number(xe,ypos2+0.7*sizexx,0.7*sizexx,
     1              real(nscalx),0.0,-1)
            else
                call gcent(xcen,ypos2,sizexx,title,0.0)
            endif
        endif
        return
        end

        subroutine dnliny(x0,y0,yleng,ymax,ymin,
     1      sizey,ticlft,lablft,dopow,lly,titley)
c-----
c       put in tic marks along the y-axis 
c
c       x0  R*4 - position of bottom side of axis
c       y0  R*4 - position of bottom side of axis
c       yleng   R*4 - length of Y axis
c       ymax    R*4 - maximum value of number corresponding to
c                   top
c       ymax    R*4 - maximum value of number corresponding to
c                   bottom
c       sizey   R*4 - size of numbers, e.g., 10
c       ticlft  L   - .true. tics go left
c                 .false. tics go right
c       lablft  L   - .true. numbers goes to left of axis
c                 .false. goes to right
c       dopow   L   - .true. put in numbers labels
c       ly  I*4 length of y-axis title string
c       titley  Ch  y-axis title string
c-----
        real*4 x0, y0, yleng, ymax, ymin, sizey
        logical ticlft, lablft, dopow
        integer lly
        character titley*(*)
        call dolny(x0,y0,yleng,ymax,ymin,
     1      sizey,ticlft,lablft,dopow,lly,titley,.false.)
        return
        end

        subroutine doliny(x0,y0,yleng,ymax,ymin,
     1      sizey,ticlft,lablft,dopow,lly,titley)
c-----
c       put in tic marks along the y-axis 
c
c       x0  R*4 - position of bottom side of axis
c       y0  R*4 - position of bottom side of axis
c       yleng   R*4 - length of Y axis
c       ymax    R*4 - maximum value of number corresponding to
c                   top
c       ymax    R*4 - maximum value of number corresponding to
c                   bottom
c       sizey   R*4 - size of numbers, e.g., 10
c       ticlft  L   - .true. tics go left
c                 .false. tics go right
c       lablft  L   - .true. numbers goes to left of axis
c                 .false. goes to right
c       dopow   L   - .true. put in numbers labels
c       ly  I*4 length of y-axis title string
c       titley  Ch  y-axis title string
c-----
        real*4 x0, y0, yleng, ymax, ymin, sizey
        logical ticlft, lablft, dopow
        integer lly
        character titley*(*)
        call dolny(x0,y0,yleng,ymax,ymin,
     1      sizey,ticlft,lablft,dopow,lly,titley,.true.)
        return
        end

        subroutine dolny(x0,y0,yleng,ymax,ymin,
     1      sizey,ticlft,lablft,dopow,lly,titley,doasis)
c-----
c       put in tic marks along the y-axis 
c
c       x0  R*4 - position of bottom side of axis
c       y0  R*4 - position of bottom side of axis
c       yleng   R*4 - length of Y axis
c       ymax    R*4 - maximum value of number corresponding to
c                   top
c       ymax    R*4 - maximum value of number corresponding to
c                   bottom
c       sizey   R*4 - size of numbers, e.g., 10
c       ticlft  L   - .true. tics go left
c                 .false. tics go right
c       lablft  L   - .true. numbers goes to left of axis
c                 .false. goes to right
c       dopow   L   - .true. put in numbers labels
c       ly  I*4 length of y-axis title string
c       titley  Ch  y-axis title string
c       doasis  L   .true. plot values
c               .false. annotate with negative (useful for depth plot)
c-----
        real*4 x0, y0, yleng, ymax, ymin, sizey
        logical ticlft, lablft, dopow
        integer lly
        character titley*(*)
        logical doasis

        character title*80
        logical dosci

        character cmd*10

        ly = lly
        if(ly.lt.0)ly = 0
c-----
c       put in the axis
c-----
        yl = y0 + yleng
        ycen = y0 + 0.5*yleng
        call plot(x0,y0,3) 
        call plot(x0,yl,2) 
c-----
c       set up positioning for tics
c-----
        tlen1 = 0.6*sizey
        tlen2 = 1.2*sizey
        if(lablft)then
            if(ticlft)then
                xpos   = x0 - 2*sizey
            else
                xpos   = x0 - sizey
            endif
            cmd = 'RIGHT'
        else
            if(ticlft)then
                xpos   = x0 + 1.00*sizey
            else
                xpos   = x0 + 2.0*sizey
            endif
            cmd = 'RIGHT'
        endif
c-----
c       compute limits for scale
c-----
        ymn = min(ymin,ymax)
        ymx = max(ymin,ymax)
        yyv = ymn
c-----
c       safety
c-----
        if(ymn.eq.ymx)ymx = ymn + 1.0
c-----
c       get some idea of the size of the number to be plotted
c-----
        yval = max(abs(ymn),abs(ymx))
c-----
c       we set the maximum number of decimal digits as 3, 
c       if the number is greater than 1000 or less than 0.01, use
c       scientific notation
c       get the maximum number of digits for the number
c-----
        if(yval.le.0.01 .and.yval.gt.0.0 )then
            dosci = .true.
            ynorm = alog10(1.0/yval)
            nynorm = ynorm+1
            ynorm = 10.0**( nynorm)
            nscaly = -nynorm
            maxdig = 4
        else if( yval.ge.10000.0)then
            dosci = .true.
            ynorm = alog10(yval)
            nynorm = ynorm
            ynorm = 10.0**(-nynorm)
            nscaly =  nynorm
            maxdig = 4
        else
            dosci = .false.
            ynorm = 1.0
            nynorm = 0.0
            nscaly = 0
            if(yval .ge. 1. .and. yval .lt. 10.)then
                        maxdig = 4
                else if(yval .ge. 10.0 .and. yval .lt. 100.)then
                        maxdig = 2
                else if(yval .ge. 100.0 .and. yval .lt. 1000.)then
                        maxdig = 2
                else
                        maxdig = 4
            endif
        endif
        if(ymn .lt. 0.0 .or. ymx.lt.0.0)then
            maxdig = maxdig + 1
        endif
        if(.not.lablft)then
            xpos = xpos + maxdig*sizey
        endif
c-----
c       choose the maximum number of tics somehow on
c       yleng, and size of numbers and necessary decimals
c-----
        dyy = ynorm*(ymx - ymn)/5.0
        if(dyy.eq.0.0)dyy = 1.0
c-----
c       get start limits for search - we will place at most 10 tics
c       here
c-----
        n = alog10(dyy)
        if(dyy.lt.1.0)n = n -1
        ifac = (10.0**(-n)*dyy)
        
c----- 
c       now get a new minimum value that is some multiple 
c       of ifac 10**n 
c----- 
        dyy = ifac * 10.0**(n)
        ilow = ynorm*ymn/dyy -1
        iup  = ynorm*ymx/dyy +1
        if(ifac.eq.1)then
            ndyyy = 5
        else if(ifac.eq.2)then
            ndyyy = 5
        else
            ndyyy = ifac
        endif
        dyyy = dyy/ndyyy

c-----
c       loop for labels
c-----
c       ye is the end of the previous label
c-----
        yc = ymn
        call plot(x0,y0,3)
c-----
c       here
c       yy is the position of the plot corresponding to a value of yyv
c       yend is the position of the end of the last number plotted
c-----
        ye = y0 
        nshftmx = 0
c-----
c       consider numbers between ilow*dyy and iup*dyy
c-----
        do 1000 i=ilow,iup
            do 1001 j=0,ndyyy-1
            
            yyv = i*dyy + j*dyyy
            
            yyv = yyv/ynorm
            if(yyv.ge.ymn .and. yyv .le.ymx)then
                    yy = y0 + (yyv - ymn)*yleng/(ymx - ymn)
                    call plot(x0,yy,3)
                    if(j.eq.0)then
                    if(doasis )then
                        yv = yyv * ynorm
                    else
                        yv = - yyv * ynorm
                    endif
                    nshft = 1
                    if(yv .lt. 0.0)nshft=nshft+1
                    if(dosci)then
                        nshft = nshft +1
                        ndec = 2
                        nshft = nshft + 2
                    else 
                        if(abs(dyy) .gt. 1)then
                            ndec = -1
                        else 
                            if(ABS(dyy).lt.0.1)then
                                ndec = 3
                            else
                                ndec = 2
                            endif
                            nshft = nshft + ndec
                        endif
                    if(abs(yv) .ge. 10.0) nshft=nshft+1
                    if(abs(yv) .ge. 100.0) nshft=nshft+1
                    if(abs(yv) .ge. 1000.0) nshft=nshft+1

                    endif
                    if(nshft .gt. nshftmx)nshftmx = nshft
                    yp = yy - 0.5*sizey
                    yb = yy
                    if((yb - ye ) .gt. 2.*sizey)then
                if(yy.ge.min(ye,yl) .and. yy.le.max(ye,yl))then
                    if(dopow)then
                call mynum(xpos,yp,sizey,yv,0.0,ndec,cmd)
                    endif
                        ye = yb 
                    endif
                        ticlen = tlen2
                    else
                        ticlen = 0.8*tlen2
                    endif
                else 
                    ticlen = tlen1
                endif
                call plot(x0,yy,3)
                if(ticlft)then  
                    call plot(x0-ticlen,yy,2)
                else
                    call plot(x0+ticlen,yy,2)
                endif
            endif
 1001       continue
 1000   continue
c-----
c       put in the title if we put in the numbers
c-----
        if(lablft)then
            if(ticlft)then
                xpos   = x0 - 2*sizey
                naxdig = -(nshftmx+1+2)
                xpos2  = x0 + naxdig* sizey -1.2*sizey
            else
                xpos   = x0 - sizey
                naxdig = -(nshftmx+0+2)
                xpos2  = x0 + naxdig* sizey -1.2*sizey
            endif
        else
            if(ticlft)then
                xpos   = x0 + 1.00*sizey
                naxdig =  (nshftmx+0+2)
                xpos2  = x0 + naxdig* sizey +1.0*sizey
            else
                xpos   = x0 + 2.0*sizey
                naxdig =  (nshftmx+1+2)
                xpos2  = x0 + naxdig* sizey +1.0*sizey
            endif
        endif
        if(dopow .and. ly .gt.0 )then
            sizeyy = 1.2*sizey
            title = ' '
            if(ly.gt.0)then
                title = titley(1:ly)
            else
                ly = 0
            endif
            if(dosci)then
                title(ly +1:ly+5)=' *10 '
                if(lablft)then
                call gcent(xpos2,ycen,sizeyy,title, 90.0)
                    ye = ycen + 0.5*sizeyy*(ly+4)
                call number(xpos2-0.7*sizey,ye,0.7*sizeyy,
     1              real(nscaly), 90.0,-1)
                else
                call gcent(xpos2,ycen,sizeyy,title,-90.0)
                    ye = ycen - 0.5*sizeyy*(ly+4)
                call number(xpos2+0.7*sizeyy,ye,0.7*sizeyy,
     1              real(nscaly),-90.0,-1)
                endif
            else
                if(lablft)then
                call gcent(xpos2,ycen,sizeyy,title, 90.0)
                else
                call gcent(xpos2,ycen,sizeyy,title,-90.0)
                endif
            endif
        endif
        return
        end


        subroutine fillit(cmd,rad,x0,y0)
        character cmd*(*)
        real xval(370), yval(370)
c-----
c       fill in a solid symbol
c-----
        ipatx = 0
        ipaty = 0
        xlen = 0.01
        ylen = 0.01
        r2 = rad 
        if(cmd(1:2).eq.'SQ')then
            jj = 0
            do 400 i=45,405,90
                ang = 3.1415927*i/180.0
                sa = sin(ang)
                ca = cos(ang)
                xnew = x0 + r2*ca
                ynew = y0 + r2*sa
                jj = jj + 1
                xval(jj) = xnew
                yval(jj) = ynew
  400       continue
            call shadep(4,xval,yval)
        else if(cmd(1:2).eq.'TR')then
            jj = 0
            do 300 i=90,450,120
                ang = 3.1415927*i/180.0
                sa = sin(ang)
                ca = cos(ang)
                xnew = x0 + r2*ca
                ynew = y0 + r2*sa
                jj = jj + 1
                xval(jj) = xnew
                yval(jj) = ynew
  300       continue
            call shadep(3,xval,yval)
        else if(cmd(1:2).eq.'HX')then
            x1 = x0 - 0.866*r2
            x2 = x0 + 0.866*r2
            x3 = x0
            y1 = y0 - 0.500*r2
            y2 = y0 - 0.500*r2
            y3 = y0 + r2
            call shadet(x1,y1,x2,y2,x3,y3,
     1          ipatx,ipaty,xlen,ylen)
            y1 = y0 + 0.500*r2
            y2 = y0 + 0.500*r2
            y3 = y0 - r2
            call shadet(x1,y1,x2,y2,x3,y3,
     1          ipatx,ipaty,xlen,ylen)
        else if(cmd(1:2).eq.'DI')then
            jj = 0
            do 200 i=0,360,90
                ang = 3.1415927*i/180.0
                sa = sin(ang)
                ca = cos(ang)
                xnew = x0 + r2*ca
                ynew = y0 + r2*sa
                jj = jj + 1
                xval(jj) = xnew
                yval(jj) = ynew
  200       continue
            call shadep(4,xval,yval)
        else if(cmd(1:2).eq.'CI')then
            jj = 0
            do 100 i=0,360,10
                ang = 3.1415927*i/180.0
                sa = sin(ang)
                ca = cos(ang)
                xnew = x0 + r2*ca
                ynew = y0 + r2*sa
                jj = jj + 1
                xval(jj) = xnew
                yval(jj) = ynew
  100       continue
            call shadep(jj,xval,yval)
        endif
        return
        end

        subroutine curvit(cmd,rad,x0,y0)
        character cmd*(*)
        real xval(370), yval(370)
c-----
c       plot an outline
c-----
        ipatx = 0
        ipaty = 0
        xlen = 0.01
        ylen = 0.01
        r2 = rad 
        if(cmd(1:2).eq.'SQ')then
            jj = 0
            do 400 i=45,405,90
                ang = 3.1415927*i/180.0
                sa = sin(ang)
                ca = cos(ang)
                xnew = x0 + r2*ca
                ynew = y0 + r2*sa
                jj = jj + 1
                xval(jj) = xnew
                yval(jj) = ynew
  400       continue
            call drawcv(jj,xval,yval)
        else if(cmd(1:2).eq.'TR')then
            jj = 0
            do 300 i=90,450,120
                ang = 3.1415927*i/180.0
                sa = sin(ang)
                ca = cos(ang)
                xnew = x0 + r2*ca
                ynew = y0 + r2*sa
                jj = jj + 1
                xval(jj) = xnew
                yval(jj) = ynew
  300       continue
            call drawcv(jj,xval,yval)
C       else if(cmd(1:2).eq.'HX')then
C           x1 = x0 - 0.866*r2
C           x2 = x0 + 0.866*r2
C           x3 = x0
C           y1 = y0 - 0.500*r2
C           y2 = y0 - 0.500*r2
C           y3 = y0 + r2
C           call shadet(x1,y1,x2,y2,x3,y3,
C     1         ipatx,ipaty,xlen,ylen)
C           y1 = y0 + 0.500*r2
C           y2 = y0 + 0.500*r2
C           y3 = y0 - r2
C           call shadet(x1,y1,x2,y2,x3,y3,
C     1         ipatx,ipaty,xlen,ylen)
        else if(cmd(1:2).eq.'DI')then
            jj = 0
            do 200 i=0,360,90
                ang = 3.1415927*i/180.0
                sa = sin(ang)
                ca = cos(ang)
                xnew = x0 + r2*ca
                ynew = y0 + r2*sa
                jj = jj + 1
                xval(jj) = xnew
                yval(jj) = ynew
  200       continue
            call drawcv(jj,xval,yval)
        else if(cmd(1:2).eq.'CI')then
            jj = 0
            do 100 i=0,360,10
                ang = 3.1415927*i/180.0
                sa = sin(ang)
                ca = cos(ang)
                xnew = x0 + r2*ca
                ynew = y0 + r2*sa
                jj = jj + 1
                xval(jj) = xnew
                yval(jj) = ynew
  100       continue
            call drawcv(jj,xval,yval)
        endif
        return
        end

        subroutine drawcv(jj,xval,yval)
        real xval(jj),yval(jj)
        call plot(xval(1), yval(1), 3)
        do 1000 i=1,jj
            call plot(xval(i), yval(i), 2)
 1000   continue
        call plot(xval(jj), yval(jj), 3)
        return
        end
        

        subroutine rclip(xmin,xmax,ymin,ymax,xc1,yc1,xc2,yc2,
     1          x0,yy0,x1,yy1,iplt)
c-----
c       perform a line clipping for the plot region
c-----
c       xmin    R*4 - minimum value of X
c       xmax    R*4 - maximum value of X
c       ymin    R*4 - minimum value of Y
c       ymax    R*4 - maximum value of Y
c       xc1 R*4 - coordinate #1 of plotted line segment
c       yc1 R*4
c       xc2 R*4 - coordinate #2 of plotted line segment
c       yc2 R*4
c       x0  R*4 - first coordinate
c       yy0 R*4
c       x1  R*4 - second coordinate
c       yy1 R*4
c       iplt    I*4 - > 0 plot the segment, otherwise do not
c-----
        real*4 xmin, xmax, ymin,ymax
        real*4 xc1, yc1, xc2, yc2
        real*4 xx0, yz0, xx1, yz1
c-----
c       Fortran implementation of Cohen and Sutherland Line 
c           Clipping Routine
c-----
        integer c0, c1, c

        iplt = 0
        xx0 = x0
        yz0 = yy0
        xx1 = x1
        yz1 = yy1
        call linecode(xx0,yz0,ymax,xmin,ymin,xmax,c0)
        call linecode(xx1,yz1,ymax,xmin,ymin,xmax,c1)
        if( xx0 .ne. xx1)slope = (yz1-yz0)/(xx1-xx0)
        if( yz0 .ne. yz1) slopeinv = (xx1-xx0)/(yz1-yz0)
 1000   continue
        if(c0 .eq. 0 .and. c1 .eq.0)go to 1001
            if( (mod(c0,2).eq.1 .and. mod(c1,2).eq.1).or.
     1      (mod(c0/2,2).eq.1 .and. mod(c1/2,2).eq.1).or.
     2      (mod(c0/4,2).eq.1 .and. mod(c1/4,2).eq.1).or.
     3      (mod(c0/8,2).eq.1 .and. mod(c1/8,2).eq.1))return

                if(c0 .eq. c1) return 
                if(c0 .eq. 0 )then
                        c = c1
                else 
                        c = c0
            endif
                if(mod(c,2).eq.1)then
                        y = yz0 + slope*(xmin-xx0)
                        x = xmin
                endif
                if(mod(c/2,2).eq.1)then 
                y = yz0 + slope*(xmax-xx0)
                        x = xmax
                endif
                if(mod(c/8,2).eq.1)then
                        x = slopeinv*(ymax-yz0) + xx0
                        y = ymax
                endif
                if(mod(c/4,2).eq.1)then
                        x = slopeinv*(ymin-yz0) + xx0
                        y = ymin
                endif
                if(c .eq. c0)then
                        xx0 = x
                        yz0 = y
                call linecode(xx0,yz0,ymax,xmin,ymin,xmax,c0)
                else 
                        xx1 = x
                        yz1 = y
                call linecode(xx1,yz1,ymax,xmin,ymin,xmax,c1)
                endif
        go to 1000
 1001   continue
        xc1 = xx0
        xc2 = xx1
        yc1 = yz0
        yc2 = yz1
        iplt = 1
        return
        end
c       
        subroutine linecode(x,y,ymax,xmin,ymin,xmax,c)
        real*4 x, y, ymax, xmin, ymin, xmax
        integer*4 c
        c = 0

        if(x .lt. xmin) then
            c = 1
        else if(x .gt. xmax)then
            c = 2
        endif
        if(y .lt. ymin) then
            c = c + 4
        else if(y.gt.ymax)then
            c = c + 8
        endif
        return
        end
c---------------------------------------------------------------------c
c                                                                     c
c      COMPUTER PROGRAMS IN SEISMOLOGY                                c
c      VOLUME I                                                       c
c                                                                     c
c      PROGRAM: NUMBER                                                c
c                                                                     c
c      COPYRIGHT (C)  1986, 1994 R. B. Herrmann                       c
c                                                                     c
c      Department of Earth and Atmospheric Sciences                   c
c      Saint Louis University                                         c
c      221 North Grand Boulevard                                      c
c      St. Louis, Missouri 63103                                      c
c      U. S. A.                                                       c
c                                                                     c
c---------------------------------------------------------------------c
      subroutine mynum (xpage,ypage,height,fpn,angle,ndec,cmd)
c.....     xpage,ypage coordinates of lower left corner of number.
c.....     height   height of plotted number.
c.....     fpn      floating point number to be plotted.
c.....     angle    angle at which number is plotted, in degrees.
c.....     ndec     number of decimal places to be drawn.
c.....     this version of number requires the symbol version with
c.....     999. x, y feature, and  nc = 0 feature.
c.....
c.....     ndec .ge. 1000 .le. 1009 invokes an exponential format
c.....     of the Fortran Format Ew.d, where d = ndec - 1000
c.....     and w = 7 + d. The output will look like
c.....     sd.dddddEsdd, where s is the sign and d is a digit
c.....     In all cases, no more than 9 significant figures are
c.....     permitted to the right of the decimal place.
c.....
c.....     ndec .ge. 2000 .le. 2009 invokes exponential but
c.....     plotted in scientific notation, e.g., 1.0 10 ^ exp
        character cmd*(*)
      character num*20
      character cexp*4
        common/Scplot/x0,y0
        common/Xcplot/xold,yold,xcur,ycur

        if(xpage .lt. 999.0) then
            x0 = xpage
        endif
        if(ypage .lt. 999.0) then
            y0 = ypage
        endif
c-----
c    to avoid dregs in memory null first character
c-----
      num(1:1)= ' '
      ii=0
c-----
      n = 0
        fp = fpn
        if(fpn.lt.0.0)then
                n=n+1
                num(n:n)='-'
                fp=abs(fp)
        endif
      if(ndec.lt.1000)then
                nzflag = 0
                mdec = -ndec
c               if(ndec.lt.0)then
c                       mdec = mdec - 1
c               endif
                if(ndec.eq.-1)then
                        nzflag = 1
                        mdec = 1
                endif
                call ffpack(num,n,fp,mdec,19,nzflag)
        else if(ndec.ge.1000)then
                       mdec = 1
                if(ndec .ge.1000 .and. ndec .le. 1009)then
                       mdec = ndec -1000
                else if(ndec .ge.2000 .and. ndec .le. 2009)then
                       mdec = ndec -2000
                else 
                       mdec = 9
                endif
                mdec = - mdec
                if(fp.eq.0.0)then
                        iex=0
                else
                        aex = alog10(fp)
                        iex = aex
c----------careful check of integer exponent
                        if(iex.lt.0.and.float(iex).gt.aex)iex=iex-1
                endif
                fp = fp * 10.** (- iex)
                nzflag = 0
                call ffpack(num,n,fp,mdec,14,nzflag)
c---put in exponent assuming iex < 100
                nman=n
                n=n+1
                num(n:n)='E'
                if(iex.lt.0)then
                        n=n+1
                        num(n:n)='-'
                        iex = -iex
                else
                        n=n+1
                        num(n:n)='+'
                endif
                nexp = n
                n=n+1
                jex = iex
                jex = mod(jex/10,10)
                num(n:n) = char(ichar('0' )+ jex)
                jex = mod(iex,10)
                n = n + 1
                num(n:n) = char(ichar('0' )+ jex)
        endif
        if(ndec .le. 1009)then
        if(cmd.eq.'LEFT')then
                call symbol(x0,y0,height,num,angle,n)
        else if(cmd.eq.'CENTER')then
                call symbol(x0-0.5*n*height,y0,height,num,
     1               angle,n)
        else if(cmd.eq.'RIGHT')then
                call symbol(x0-n*height,y0,height,num,angle,n)
        endif
        else if(ndec .ge. 2000 .and. ndec.le.2009)then
c-----
c       save the exponent, replace the E-EX with 'x10'
c-----
                cexp(1:3) = num(nexp:n)
                num(nman+1:nman+3)='x10'
            xx0 = x0
            yy0 = y0
                call symbol(x0,y0,height,num(1:nman+3),
     1               angle,nman+3)
c-----
c       get the proper position because of rotation
c-----
                ang = angle*3.1415927/180.0
                ca = cos(ang)
                sa = sin(ang)
                xl = height*(nman+3)
                yl = 0.6*height
            xx0 = xx0 + xl*ca
            yy0 = yy0 + xl*sa
                xxx0 = xx0  - yl*sa
                yyy0 = yy0  + yl*ca
                ht = 0.7*height
                call symbol(xxx0,yyy0,ht,cexp(1:3),angle,3)
c-----
c           now position at end of proper string
c-----
            xl =  3.0*ht
            xcur = xx0 + xl*ca
            ycur = yy0 + xl*sa
            x0 = xcur
            y0 = ycur
            
c-----
        endif
        return
        end
        subroutine ffpack(num,n,fpv,mdec,mwid,nzflag)
        character num*(*)
        fp = fpv
        nzflg = nzflag
c-----since we have a maximum field width of mwid and since
c-----the . takes a position, we can only have mwid -1 figures
      m = mdec 
        maxm = 9
        if(m.gt.maxm)m=maxm
        if(m.lt. -maxm)m= -maxm
        mm = m
        if(m.gt.0) mm = mm - 1
c----do a simple rounding up
        fp = fp + 0.5 * 10.** (mm)
c---- get position of first significant figure
c
c     5 4 3 2 1 0 -1 -2 -3
c
c     d d d d d .  d  d  d
c
c----

        iex = alog10(fp) 
        if(fp.ge.1.0)then
                iex = iex + 1
        else
                iex = iex - 1
        endif
c----
c----   procede from most significant down
c----   but never less than 10**-9
c----
        
        if(fp.le.1.0e-9)then
                fp = 0.0
                iex = -9
        endif
        ilw = mdec
        if(fp.ge.1.0)then
                iup=iex
        else
                iup = 1
        endif
        if(m.le.0)then
                ilw= m+1
        else
                ilw = m 
        endif
        if(iex.gt.9)then
                ilw = 1
                nzflg = 1
        endif
        do 100 ipos=iup,ilw,-1
                k = fp * 10.**( -ipos +1 )
                if(n.lt.mwid)then
                        n = n + 1
                        num(n:n) = char(ichar('0')+k)
                endif
                fp = fp - (float(k) * 10.**(ipos-1))
        if(ipos.eq.1.and.nzflg.eq.0.and.n.lt.mwid)then
                        n = n + 1
                        num(n:n) = '.'
                endif
 100    continue
        return
        end

      subroutine msupsc (xpage,ypage,height,str,npow,angle)
c.....     xpage,ypage coordinates of lower left corner of number.
c.....     height   height of plotted number.
c.....     str  Character string
c.....     npow power
c.....     fpn      floating point number to be plotted.
c.....     angle    angle at which number is plotted, in degrees.
        character str*(*)
        common/Scplot/x0,y0
        common/Xcplot/xold,yold,xcur,ycur

        if(xpage .lt. 999.0) then
            x0 = xpage
        endif
        if(ypage .lt. 999.0) then
            y0 = ypage
        endif
        xx0 = x0
        yy0 = y0
        il = lgstr(str)
        call symbol(x0,y0,height,str(1:il), angle,il)
c-----
c       get the proper position because of rotation
c-----
                ang = angle*3.1415927/180.0
                ca = cos(ang)
                sa = sin(ang)
                xl = height*(il)
                yl = 0.6*height
            xx0 = xx0 + xl*ca
            yy0 = yy0 + xl*sa
                xxx0 = xx0  - yl*sa
                yyy0 = yy0  + yl*ca
                ht = 0.7*height
                call number(xxx0,yyy0,ht,real(npow),angle,-1)
c-----
c           now position at end of proper string
c-----
            xl =  3.0*ht
            xcur = xx0 + xl*ca
            ycur = yy0 + xl*sa
            x0 = xcur
            y0 = ycur
            
c-----
        return
        end

        subroutine gsubsc(x,y,ht,s1,n1,s2,n2)
        character s1*(*), s2*(*)
        common/savfon/infont
        if(n1.lt.0)then
            call gfont(4)
            nn = -n1
        else
            nn = n1
            call gfont(infont)
        endif
        if(nn.gt.0)call symbol(x,y,ht,s1,0.0,nn)
        if(n2.lt.0)then
            call gfont(4)
            n = -n2
        else
            n = n2
            call gfont(infont)
        endif
        if(n.gt.0)call symbol(x+nn*ht+0.4*ht,y-0.4*ht,0.7*ht,s2,0.0,n)
        if(n1.lt.0 .or. n2.lt.0)call gfont(infont)
        return
        end
        
        subroutine gsupsc(x,y,ht,s1,n1,s2,n2)
        character s1*(*), s2*(*)
        common/savfon/infont
        if(n1.lt.0)then
            call gfont(4)
            nn = -n1
        else
            nn = n1
            call gfont(infont)
        endif
        if(nn.gt.0)call symbol(x,y,ht,s1,0.0,nn)
        if(n2.lt.0)then
            call gfont(4)
            n = -n2
        else
            n = n2
            call gfont(infont)
        endif
        if(n.gt.0)call symbol(x+nn*ht+0.4*ht,y+0.4*ht,0.7*ht,s2,0.0,n)
        if(n1.lt.0 .or. n2.lt.0)call gfont(infont)
        return
        end

        subroutine gsubsup(x,y,ht,s1,n1,s2,n2,s3,n3)
        character s1*(*), s2*(*), s3*(*)
        common/savfon/infont
        if(n1.lt.0)then
            call gfont(4)
            nn = -n1
        else
            nn = n1
            call gfont(infont)
        endif
        if(nn.gt.0)call symbol(x,y,ht,s1,0.0,nn)
        if(n2.lt.0)then
            call gfont(4)
            n = -n2
        else
            n = n2
            call gfont(infont)
        endif
        if(n.gt.0)call symbol(x+nn*ht+0.4*ht,y-0.4*ht,0.7*ht,s2,0.0,n)
        if(n3.lt.0)then
            call gfont(4)
            n = -n3
        else
            n = n3
            call gfont(infont)
        endif
        if(n.gt.0)call symbol(x+nn*ht+0.4*ht,y+0.4*ht,0.7*ht,s3,0.0,n)
        if(n1.lt.0 .or. n2.lt.0 .or. n3.lt.0)call gfont(infont)
        return
        end
