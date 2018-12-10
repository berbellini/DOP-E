        program calplt
c-----
c       changes
c       26 JAN 2010
c       add curvit as a command to be able to place a border
c       around symbol defineed using the sfill command
c-----
c       program to to script based CALPLOT plots
c-----
        common/savfon/infont
        common/savpen/ipen,ipencolor
        character cmd*20, string*120
        character t1*80,s1*80,t2*80,s2*80,t3*80,s3*80
        character u1*80,u2*80,u3*80
        logical verby
        logical bold
        integer lstr
        character lq*2, rq*2

        lq = ' '''
        rq = ''' '

c----
c       parse command line
c-----
        call gcmdln(verby)
        if(verby)call help()
c----
c       initialize graphics
c-----
        call pinitf('CALPLT.PLT')
        infont = 0
        ipencolor = 1
        widsav = 0.001
        
c-----
c       open file for saving all commands
c-----
        open(1,file='CALPLT.cmd',status='unknown',form=
     1      'formatted',access='sequential')
        rewind 1
c-----
c       begin processing
c-----
 1000   continue
        if(verby)write(6,*)'Enter command'
        read(5,'(a)',end=9999,err=9999)cmd
        if(cmd(1:4).ne.'HELP')write(1,'(a)')cmd
        if(cmd(1:6).eq.'SYMBOL')then
            if(verby)then
                write(6,*)'SYMBOL:x,y,ht,string,angle,nchar'
            endif
            read(5,*)xx,yy,ht,string,angle,nchar
            ls = lstr(string)
            write(1,*)xx,yy,ht,lq,string(1:ls),rq,angle,nchar
            if(nchar.gt.0)then
                nchar = lstr(string)
                call symbol(xx,yy,ht,string,angle,nchar)
            else
                read(string,'(i10)')ival
                call symbol(xx,yy,ht,char(ival),angle,-1)
            endif
        else if(cmd(1:5).eq.'ARRAY')then
            if(verby)then
                write(6,*)'ARRAY: Enter file of x,y pairs in quotes'
            endif
            read(5,*)string
            ls = lstr(string)
            write(1,*)lq,string(1:ls),rq
            call doxy(string)
        else if(cmd(1:6).eq.'CIRCLE')then
            if(verby)then
                write(6,*)'CIRCLE:rad,x0,y0'
            endif
            read(5,*)rad,x0,y0
            write(1,*)rad,x0,y0
            call circle(rad,x0,y0)
        else if(cmd(1:6).eq.'SQUARE')then
            if(verby)then
                write(6,*)'SQUARE:rad,x0,y0'
            endif
            read(5,*)rad,x0,y0
            write(1,*)rad,x0,y0
            call square(rad,x0,y0)
        else if(cmd(1:7).eq.'HEXAGON')then
            if(verby)then
                write(6,*)'HEXAGON:rad,x0,y0'
            endif
            read(5,*)rad,x0,y0
            write(1,*)rad,x0,y0
            call hexagon(rad,x0,y0)
        else if(cmd(1:6).eq.'SFILL')then
            if(verby)then
                write(6,*)
     1          'SFILL:string_type,height,x0,y0 where', 
     1              ' type = SQuare, CIrcle',
     2              ' DIamond, TRiangle, HXagon' 
            endif
            read(5,*)string,rad,x0,y0
            ls = lstr(string)
            write(1,*)lq,string(1:ls),rq,rad,x0,y0
            call fillit(string(1:ls),rad,x0,y0,.true.)
        else if(cmd(1:4).eq.'CURV')then
            if(verby)then
                write(6,*)
     1          'CURVIT:string_type,height,x0,y0 where', 
     1              ' type = SQuare, CIrcle',
     2              ' DIamond, TRiangle, HXagon' 
            endif
            read(5,*)string,rad,x0,y0
            ls = lstr(string)
            write(1,*)lq,string(1:ls),rq,rad,x0,y0
            call fillit(string(1:ls),rad,x0,y0,.false.)
        else if(cmd(1:3).eq.'ARC')then
            if(verby)then
                write(6,*)'ARC:rad,x0,y0,ang1,ang2'
            endif
            read(5,*)rad,x0,y0,ang1,ang2
            write(1,*)rad,x0,y0,ang1,ang2
            call arc(rad,x0,y0,ang1,ang2)
        else if(cmd(1:4).eq.'HELP')then
            call help()
        else if(cmd(1:6).eq.'FACTOR')then
            if(verby)then
                write(6,*)
     1              'FACTOR:scaling_factor'
            endif
            read(5,*)fct
            write(1,*)fct
            call factor(fct)
        else if(cmd(1:7).eq.'COMMENT')then
            if(verby)then
                write(6,*)'Enter comment (quotes)'
            endif
            read(5,*)string
            ls = lstr(string)
            write(1,*)lq,string(1:ls),rq
        elseif(cmd(1:6).eq.'CENTER')then
            if(verby)then
                write(6,*)'CENTER:x,y,ht,string,angle'
            endif
            read(5,*)xx,yy,ht,string,angle
            ls = lstr(string)
            write(1,*)xx,yy,ht,lq,string(1:ls),rq,angle
            call center(xx,yy,ht,string,nchar,angle)
            call symbol(xx,yy,ht,string,angle,nchar)
        elseif(cmd(1:5).eq.'RIGHT')then
            if(verby)then
                write(6,*)'RIGHT:x,y,ht,string,angle'
            endif
            read(5,*)xx,yy,ht,string,angle
            ls = lstr(string)
            write(1,*)xx,yy,ht,lq,string(1:ls),rq,angle
            call right(xx,yy,ht,string,nchar,angle)
            call symbol(xx,yy,ht,string,angle,nchar)
        elseif(cmd(1:4).eq.'LEFT')then
            if(verby)then
                write(6,*)'LEFT:x,y,ht,string,angle'
            endif
            read(5,*)xx,yy,ht,string,angle
            ls = lstr(string)
            write(1,*)xx,yy,ht,lq,string(1:ls),rq,angle
            call symbol(xx,yy,ht,string,angle,ls)
        elseif(cmd(1:6).eq.'NUMBER')then
            if(verby)then
                write(6,*)'NUMBER:x,y,ht,number,angle,nchar'
            endif
            read(5,*)xx,yy,ht,fpn,angle,ndec
            write(1,*)xx,yy,ht,fpn,angle,ndec
            call number(xx,yy,ht,fpn,angle,ndec)
        elseif(cmd(1:6).eq.'NEWPEN')then
            if(verby)then
                write(6,*)'NEWPEN:ipen'
            endif
            read(5,*)ipencolor
            write(1,*)ipencolor
            call newpen(ipencolor)
        elseif(cmd(1:5).eq.'PLOTD')then
            if(verby)then
                write(6,*)'PLOTD:x,y,ipat,xlen'
            endif
            read(5,*)xx,yy,ipat,xlen
            write(1,*)xx,yy,ipat,xlen
            call plotd(xx,yy,ipat,xlen)
        elseif(cmd(1:6).eq.'SHADER')then
            if(verby)then
                write(6,*)'SHADER:x1,y1,x2,y2,ipatx,ipaty,xlen,ylen'
            endif
            read(5,*)x1,y1,x2,y2,ipatx,ipaty,xlen,ylen
            write(1,*)x1,y1,x2,y2,ipatx,ipaty,xlen,ylen
            call shader(x1,y1,x2,y2,ipatx,ipaty,xlen,ylen)
        else if(cmd(1:4) .eq. 'PLOT')then
            if(verby)then
                write(6,*)'PLOT:x,y,ipen'
            endif
            read(5,*)xx,yy,ipen
            write(1,*)xx,yy,ipen
            call plot(xx,yy,ipen)
        else if(cmd(1:5) .eq. 'FRAME')then
            call frame()
        else if(cmd(1:6) .eq. 'GWIDTH')then
            if(verby)then
                write(6,*)'GWIDTH:width'
            endif
            read(5,*)width
            write(1,*)width
            widsav = width
            call gwidth(width)
        else if(cmd(1:5) .eq. 'GFONT')then
            if(verby)then
                write(6,*)'GFONT:ifont'
            endif
            read(5,*)ifont
            write(1,*)ifont
            call gfont(ifont)
            infont = ifont
        else if(cmd(1:5) .eq. 'GUNIT')then
            if(verby)then
                write(6,*)'GUNIT:in or cm (string)'
            endif
            read(5,*)string
            ls = lstr(string)
            write(1,*)lq,string(1:ls),rq
            if(string(1:2) .eq. 'in' .or. string(1:2) .eq. 'IN')then
                call gunit('in')
            else if(string(1:2).eq. 'cm' .or. string(1:2).eq.'CM')then
                call gunit('cm')
            endif
        else if(cmd(1:4).eq.'PEND')then
            go to 9999
        else if(cmd(1:4).eq.'LINE')then
            if(verby)then
                write(6,*)'LINE:xl,yl,xu,yu'
            endif
            read(5,*)xl,yl,xu,yu
            write(1,*)xl,yl,xu,yu
            call plot(xl,yl,3)
            call plot(xu,yu,2)
        else if(cmd(1:5).eq.'ARROW')then
            if(verby)then
                write(6,*)'ARROW:xs,ys,xe,ye,bold(TF)'
            endif
            read(5,*)xs,ys,xe,ye,bold
            write(1,*)xs,ys,xe,ye,bold
            call arrow(xs,ys,xe,ye,bold)
        else if(cmd(1:5).eq.'SUBSC')then
            if(verby)then
                write(6,*)'SUBSC:x,y,ht,s1,n1,s2,n2',
     1              ' - n1,n2  neg = font 4'
            endif
            read(5,*)x,y,ht,s1,n1,s2,n2
            ls1 = lstr(s1)
            ls2 = lstr(s2)
            write(1,*)x,y,ht,lq,s1(1:ls1),rq,n1,lq,s2(1:ls2),rq,n2
            call subsc(x,y,ht,s1,n1,s2,n2)
        else if(cmd(1:5).eq.'SUPSC')then
            if(verby)then
                write(6,*)'SUPSC:x,y,ht,s1,n1,s2,n2',
     1              ' - n1,n2  neg = font 4'
            endif
            read(5,*)x,y,ht,s1,n1,s2,n2
            ls1 = lstr(s1)
            ls2 = lstr(s2)
            write(1,*)x,y,ht,lq,s1(1:ls1),rq,n1,lq,s2(1:ls2),rq,n2
            call supsc(x,y,ht,s1,n1,s2,n2)
        else if(cmd(1:6).eq.'SUBSUP')then
            if(verby)then
                write(6,*)'SUBSUP,y,ht,s1,n1,s2,n2,s3,n3',
     1              ' - n1,n2,n3  neg = font 4'
            endif
            read(5,*)x,y,ht,s1,n1,s2,n2,s3,n3
            ls1 = lstr(s1)
            ls2 = lstr(s2)
            ls3 = lstr(s3)
            write(1,*)x,y,ht,lq,s1(1:ls1),rq,n1,lq,s2(1:ls2),rq,
     1          n2,lq,s3(1:ls3),rq,n3
            call subsup(x,y,ht,s1,n1,s2,n2,s3,n3)
        else if(cmd(1:5).eq.'AXES3' .and. cmd.ne.'AXES3SU')then
            if(verby)then
                write(6,*)'AXES3:x0,y0,xl,yl,zl,',
     1              'ang1,ang2,ang3,',
     2              't1,s1,t2,s2,t3,s3,cmd'
            endif
            read(5,*) x0,y0,xl,yl,zl,ang1,
     1          ang2,ang3,t1,s1,t2,s2,t3,s3,string
            ls1 = lstr(s1)
            lt1 = lstr(t1)
            ls2 = lstr(s2)
            lt2 = lstr(t2)
            ls3 = lstr(s3)
            lt3 = lstr(t3)
            ls = lstr(string)
            write(1,*) x0,y0,xl,yl,zl,ang1,
     1          ang2,ang3,
     2          lq,t1(1:lt1),rq,lq,s1(1:ls1),rq,
     3          lq,t2(1:lt2),rq,lq,s2(1:ls2),rq,
     4          lq,t3(1:lt3),rq,lq,s3(1:ls3),rq,
     5          lq,string(1:ls),rq
            call axes3su(x0,y0,xl,yl,zl,ang1,
     1          ang2,ang3,t1,s1,' ',t2,s2,' ',t3,s3,' ',string)
        else if(cmd(1:7).eq.'AXES3SU')then
            if(verby)then
                write(6,*)
     1          'AXES3SU:x0,y0,xl,yl,zl,ang1,ang2,ang3,t1,s1,u1',
     2           't2,s2,u2,t3,s3,u3,cmd'
            endif
            read(5,*) x0,y0,xl,yl,zl,ang1,
     1          ang2,ang3,t1,s1,u1,t2,s2,u2,t3,s3,u3,string
            ls1 = lstr(s1)
            lt1 = lstr(t1)
            lu1 = lstr(u1)
            ls2 = lstr(s2)
            lt2 = lstr(t2)
            lu2 = lstr(u2)
            ls3 = lstr(s3)
            lt3 = lstr(t3)
            lu3 = lstr(u3)
            ls = lstr(string)
            write(1,*) x0,y0,xl,yl,zl,ang1,ang2,ang3,
     1          lq,t1(1:lt1),rq,lq,s1(1:ls1),rq,lq,u1(1:lu1),rq,
     2          lq,t2(1:lt2),rq,lq,s2(1:ls2),rq,lq,u2(1:lu2),rq,
     3          lq,t3(1:lt3),rq,lq,s3(1:ls3),rq,lq,u3(1:lu3),rq,
     4          lq,string(1:ls),rq
            call axes3su(x0,y0,xl,yl,zl,ang1,
     1          ang2,ang3,t1,s1,u1,t2,s2,u2,t3,s3,u3,string)
        else if(cmd(1:3).eq.'BOX')then
            if(verby)then
                write(6,*)'BOX:xl,yl,xu,yu,string,ht'
            endif
            read(5,*)xl,yl,xu,yu,string,ht
            ls = lstr(string)
            write(1,*)xl,yl,xu,yu,lq,string(1:ls),rq,ht
            call newpen(0)
            call shader(xl,yl,xu,yu,0,0,0.1,0.1)
            call newpen(ipencolor)
            call plot(xl,yl,3)
            call plot(xl,yu,2)
            call plot(xu,yu,2)
            call plot(xu,yl,2)
            call plot(xl,yl,2)
            l = lstr(string)
            if(l.gt.0)then
                if(ht.eq.0.0)then
                    ht1 = abs(xl-xu)/real(l)
                    ht2 = abs(yl-yu)
                    if(ht2 .lt. ht1)then
                        ht = 0.6*ht2
                    else
                        ht = 0.6*ht1
                    endif
                endif
                xx = 0.5*(xu+xl)
                yy = 0.5*(yu+yl) - 0.5*ht
                angle = 0.0
                call center(xx,yy,ht,string,nchar,angle)
                call gwidth(0.0)
                call symbol(xx,yy,ht,string,angle,nchar)
                call gwidth(widsav)
            endif
            
        endif
        goto 1000
 9999   continue
        call pend()
        close(1)
        end

        function lstr (str)
c-----
c     this routine determines the length of a character string.
c-----
        character*(*) str
        integer*4 lstr
        lstr=len(str)
        do 10 i=lstr,1,-1
            if (str(i:i).ne.' ') go to 20
  10    continue
  20    lstr=i
        return
        end

        subroutine center(xx,yy,ht,string,il,angle)
            character string*(*)
            il = lstr(string)
            t = angle*3.1415927/180.0
            ct = cos(t)
            st = sin(t)
            rl = 0.5*il*ht
            
            xx = xx - rl*ct
            yy = yy - rl*st
        return
        end

        subroutine right(xx,yy,ht,string,il,angle)
            character string*(*)
            il = lstr(string)
            t = angle*3.1415927/180.0
            ct = cos(t)
            st = sin(t)
            rl = 1.0*il*ht
            
            xx = xx - rl*ct
            yy = yy - rl*st
        return
        end

        subroutine gcmdln(verby)
c-----
c       parse command line arguments and return control
c       parameters
c
c       requires subroutine mgtarg(i,name) to return
c         the i'th argument in the string name
c
c       and the function mnmarg() to return the number
c         of arguments excluding program name
c         The first argument is i = 1
c
c-----
c       verby   L    - verbose output
c-----
        character name*20
        integer *4 mnmarg
        logical verby
c-----
c       set up defaults for poor usage test
c-----
        verby = .false.
        nmarg = mnmarg()
        i = 0
   11   i = i + 1
        if(i.gt.nmarg)goto 13
        call mgtarg(i,name)
        if(name(1:2).eq.'-V')then
            verby = .true.
        endif
        go to 11
   13   continue
        return
        end

        subroutine circle(rad,x0,y0)
c-----
c       draw a circle of radius rad centered at (x0,y0)
c-----
        degrad = 3.1415927/180.0
        ipen = 3
        do 100 i=0,360
            xx = x0 + rad*cos(degrad*i)
            yy = y0 + rad*sin(degrad*i)
            call plot(xx,yy,ipen)
            ipen = 2
  100   continue
        call plot(xx,yy,3)
        return
        end

        subroutine square(rad,x0,y0)
c-----
c       draw a square of leg rad centered at (x0,y0)
c-----
        x1 = x0 - 0.5 * rad
        y1 = y0 - 0.5 * rad
        x2 = x0 + 0.5 * rad
        y2 = y0 + 0.5 * rad
        call plot(x1,y1,3)
        call plot(x2,y1,2)
        call plot(x2,y2,2)
        call plot(x1,y2,2)
        call plot(x1,y1,2)
        call plot(x1,y1,3)

        return
        end

        subroutine hexagon(rad,x0,y0)
        real*4 x0, y0, rad
        real xv(7), yv(7)
c-----
c       draw a octagon of side rad  centered at (x0,y0)
c-----
        degrad = 3.1415927/180.0
c-----
c       if necessary, reset the region behind the focal circle to the
c       background color
c-----
c-----
c       generate corners of octagon
c-----  
        degrad = 3.1415927/180.0
        j = 1
        do 100 i=0,360,60
            ang = degrad * (-30.0 + i )
            c = cos(ang)
            s = sin(ang)
            xv(j) = x0 + 0.5*rad * c
            yv(j) = y0 + 0.5*rad * s
            j = j + 1
  100   continue
        do 200 i=1,7
            if(i.eq.1)then
                call plot(xv(i),yv(i),3)
            else 
                call plot(xv(i),yv(i),2)
            endif
 200    continue
        call plot(xv(1),yv(1),3)
        return
        end

        subroutine fillit(cmd,rad,x0,y0,dofill)
        character cmd*(*)
        real xval(370), yval(370)
        logical dofill
c-----
c       fill in a solid symbol
c-----
        ipatx = 0
        ipaty = 0
        xlen = 0.01
        ylen = 0.01
        r2 = rad / 2.0
        if(cmd(1:2).eq.'SQ')then
            dx = r2
            dy = r2
            x1 = x0 - dx
            x2 = x0 + dx
            y1 = y0 - dy
            y2 = y0 + dy
            xval(1) = x1
            yval(1) = y1
            xval(2) = x2
            yval(2) = y1
            xval(3) = x2
            yval(3) = y2
            xval(4) = x1
            yval(4) = y2
            xval(5) = x1
            yval(5) = y1
            if(dofill)then
                  call shadep(4,xval,yval)
            else
                  call drawp(5,xval,yval)
            endif
C           call shader(x1,y1,x2,y2,ipatx,ipaty,xlen,ylen)
        else if(cmd(1:2).eq.'TR')then
            x1 = x0 - 0.866*r2
            x2 = x0 + 0.866*r2
            x3 = x0
            y1 = y0 - 0.500*r2
            y2 = y0 - 0.500*r2
            y3 = y0 + r2
            xval(1) = x1
            yval(1) = y1
            xval(2) = x2
            yval(2) = y2
            xval(3) = x3
            yval(3) = y3
            xval(4) = x1
            yval(4) = y1
            if(dofill)then
                  call shadep(3,xval,yval)
            else
                  call drawp(4,xval,yval)
            endif
C           call shadet(x1,y1,x2,y2,x3,y3,
C     1         ipatx,ipaty,xlen,ylen)
        else if(cmd(1:2).eq.'HX')then
            x1 = x0 - 0.866*r2
            x2 = x0 + 0.866*r2
            x3 = x0
            y1 = y0 - 0.500*r2
            y2 = y0 - 0.500*r2
            y3 = y0 + r2
            xval(1) = x1
            yval(1) = y1
            xval(2) = x2
            yval(2) = y2
            xval(3) = x3
            yval(3) = y3
            xval(4) = x1
            yval(4) = y1
C            call shadet(x1,y1,x2,y2,x3,y3,
C     1          ipatx,ipaty,xlen,ylen)
            if(dofill)then
                  call shadep(3,xval,yval)
            else
                  call drawp(4,xval,yval)
            endif
            y1 = y0 + 0.500*r2
            y2 = y0 + 0.500*r2
            y3 = y0 - r2
            xval(1) = x1
            yval(1) = y1
            xval(2) = x2
            yval(2) = y2
            xval(3) = x3
            yval(3) = y3
            xval(4) = x1
            yval(4) = y1
            if(dofill)then
                  call shadep(3,xval,yval)
            else
                  call drawp(4,xval,yval)
            endif
C            call shadet(x1,y1,x2,y2,x3,y3,
C     1          ipatx,ipaty,xlen,ylen)
        else if(cmd(1:2).eq.'DI')then
            dx = r2
            dy = r2
C           call shadet(x0,y0+dy,x0-dx,y0,x0+dx,y0,
C     1         ipatx,ipaty,xlen,ylen)
C           call shadet(x0,y0-dy,x0-dx,y0,x0+dx,y0,
C     1         ipatx,ipaty,xlen,ylen)
            xval(1) = x0 - dx
            yval(1) = y0
            xval(2) = x0 
            yval(2) = y0 - dy
            xval(3) = x0 + dx 
            yval(3) = y0 
            xval(4) = x0  
            yval(4) = y0 + dy
            xval(5) = x0 - dx
            yval(5) = y0
            if(dofill)then
                  call shadep(4,xval,yval)
            else
                  call drawp(5,xval,yval)
            endif
        else if(cmd(1:2).eq.'CI')then
            xold = x0 + r2
            yold = y0
            jj = 0
            do 100 i=0,360,10
                ang = 3.1415927*i/180.0
                sa = sin(ang)
                ca = cos(ang)
                xnew = x0 + r2*ca
                ynew = y0 + r2*sa
C               call shadet(x0,y0,xold,yold,xnew,ynew,
C     1             ipatx,ipaty,xlen,ylen)
                xold = xnew
                yold = ynew
                jj = jj + 1
                xval(jj) = xnew
                yval(jj) = ynew
  100       continue
            xval(jj+1) = xval(1)
            yval(jj+1) = yval(1)
            if(dofill)then
                  call shadep(jj,xval,yval)
            else
                  call drawp(jj+1,xval,yval)
            endif
        endif
        return
        end

        subroutine drawp(jj,xval,yval)
        real xval(jj),yval(jj)
        call plot(xval(1), yval(1), 3)
        do 1000 i=1,jj
            call plot(xval(i), yval(i), 2)
 1000   continue
        call plot(xval(jj), yval(jj), 3)
        return
        end

            

        subroutine arc(rad,x0,y0,ang1,ang2)
c-----
c       draw a arc of radius rad centered at (x0,y0)
c-----
        degrad = 3.1415927/180.0
        i1 = ang1
        i2 = ang2
        if(i2.lt.i1)then
            inc = -1
        else
            inc = 1
        endif
        ipen = 3
        do 100 i=i1,i2,inc
            xx = x0 + rad*cos(degrad*i)
            yy = y0 + rad*sin(degrad*i)
            call plot(xx,yy,ipen)
            ipen = 2
  100   continue
        call plot(xx,yy,3)
        return
        end

        subroutine arrow(xx,yy,xxx,yyy,bold)
        logical bold
        if(bold)then
            call gwidth(0.06)
        else
            call gwidth(0.02)
        endif
        call plot(xx,yy,3)
        call plot(xxx,yyy,2)
        call plot(xx,yy,2)
c-----
c       get direction cosines for arrow tip
c-----
        dx = xxx - xx
        dy = yyy - yy
        r = sqrt(dx**2 + dy**2)
        if(r .le. 1.0e-5)return
        c = dx/r
        s = dy/r
        dx = -0.2*c - 0.1*s
        dy = -0.2*s + 0.1*c
        call plot(xxx+dx,yyy+dy,3)
        dx = -0.2*c + 0.1*s
        dy = -0.2*s - 0.1*c
        call plot(xxx, yyy, 2)
        call plot(xxx+dx,yyy+dy,3)
        call plot(xxx, yyy, 2)
        call plot(xxx, yyy, 3)
        if(bold)call gwidth(0.02)

        return
        end
        
        subroutine subsc(x,y,ht,s1,n1,s2,n2)
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
        if(n.gt.0)call symbol(x+nn*ht,y-0.4*ht,0.7*ht,s2,0.0,n)
        if(n1.lt.0 .or. n2.lt.0)call gfont(infont)
        return
        end
        
        subroutine subsup(x,y,ht,s1,n1,s2,n2,s3,n3)
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
        if(n.gt.0)call symbol(x+nn*ht,y-0.4*ht,0.7*ht,s2,0.0,n)
        if(n3.lt.0)then
            call gfont(4)
            n = -n3
        else
            n = n3
            call gfont(infont)
        endif
        if(n.gt.0)call symbol(x+nn*ht,y+0.4*ht,0.7*ht,s3,0.0,n)
        if(n1.lt.0 .or. n2.lt.0 .or. n3.lt.0)call gfont(infont)
        return
        end
        
        subroutine supsc(x,y,ht,s1,n1,s2,n2)
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
        if(n.gt.0)call symbol(x+nn*ht,y+0.4*ht,0.7*ht,s2,0.0,n)
        if(n1.lt.0 .or. n2.lt.0)call gfont(infont)
        return
        end

        subroutine axes3su(x0,y0,x1,x2,x3,a1,a2,a3
     1      ,t1,s1,u1,t2,s2,u2,t3,s3,u3,cmd)
        character t1*(*),s1*(*),t2*(*),s2*(*),t3*(*),s3*(*)
        character u1*(*),u2*(*),u3*(*)
        character tit*80, sub*80, sup*80
        character cmd*(*)

        call gwidth(0.02)

c-----
c       process each axis
c-----
        do 1000 i=1,3
            if(i.eq.1)then
                tit = t1
                sub = s1
                sup = u1
                xl = x1
                ang = a1
            else if(i.eq.2)then
                tit = t2
                sub = s2
                sup = u2
                xl = x2
                ang = a2
            else if(i.eq.3)then
                tit = t3
                sub = s3
                sup = u3
                xl = x3
                ang = a3
            endif
            lt = lstr(tit)
            ls = lstr(sub)
            if(ls.eq.1 .and. sub(1:1).eq.' ')lu = 0
            lu = lstr(sup)
            if(lu.eq.1 .and. sup(1:1).eq.' ')lu = 0
            ht = 0.3
            hs = 0.7*ht
            ag = ang*3.1415927/180.0
            s = sin(ag)
            c = cos(ag)
            xn = xl*c
            yn = xl*s
            xo = 0.5*ht*c
            yo = 0.5*ht*s
            if(yo.lt.0.0)then
                xo = 3.0*xo
                yo = 3.0*yo
            endif
            if(cmd.eq.'F')then
                call arrow(x0-xn,y0-yn,x0+xn,y0+yn,.false.)
            else
                call arrow(x0,y0,x0+xn,y0+yn,.false.)
            endif
            call symbol(x0+xn+xo,y0+yn+yo,ht,tit,0.0,lt)
            if(ls .gt. 0)then
            call symbol(x0+xn+xo+lt*ht,y0+yn+yo-0.6*hs,hs,sub,0.0,ls)
            endif
            if(lu .gt. 0)then
            call symbol(x0+xn+xo+lt*ht,y0+yn+yo+0.6*hs,hs,sup,0.0,lu)
            endif
 1000   continue
        return
        end

        subroutine help()
c-----
c       provide assistance
c-----
        write(6,*)
     1  'FACTOR:scaling_factor'
        write(6,*)
     1  'SYMBOL:x,y,ht,string,angle,nchar'
        write(6,*)
     1  'CIRCLE:rad,x0,y0'
        write(6,*)
     1  'SFILL:string_type,height,x0,y0 where type = SQuare, CIrcle',
     2       'DIamond, TRiangle, HXagon' 
        write(6,*)
     1  'CURVIT:string_type,height,x0,y0 where type = SQuare, CIrcle',
     2       'DIamond, TRiangle, HXagon' 
        write(6,*)
     1  'ARC:rad,x0,y0,ang1,ang2'
        write(6,*)
     1  'HELP: (show this menu)'
        write(6,*)
     1  'CENTER:x,y,ht,string,angle'
        write(6,*)
     1  'LEFT:x,y,ht,string,angle'
        write(6,*)
     1  'RIGHT:x,y,ht,string,angle'
        write(6,*)
     1  'NUMBER:x,y,ht,string,angle,nchar'
        write(6,*)
     1  'NEWPEN:ipen'
        write(6,*)
     1  'PLOTD:x,y,ipat,xlen'
        write(6,*)
     1  'SHADER:x1,y1,x2,y2,ipatx,ipaty,xlen,ylen'
        write(6,*)
     1  'PLOT:x,y,ipen'
        write(6,*)
     1  'FRAME: (start new plot frame)'
        write(6,*)
     1  'GWIDTH:width'
        write(6,*)
     1  'GFONT:ifont'
        write(6,*)
     1  'GUNIT:in or cm (string)'
        write(6,*)
     1  'PEND:'
        write(6,*)
     1  'LINE:xl,yl,xu,yu'
        write(6,*)
     1  'BOX:xl,yl,xu,yu,string,ht'
        write(6,*)
     1  'COMMENT: comment for CALPLT.PLT in quotes'
        write(6,*)
     1  'ARROW:xs,ys,xe,ye,bold(TF)'
        write(6,*)
     1  'SUBSC:x,y,ht,s1,n1,s2,n2: - n1,n2 or neg = font 4'
        write(6,*)
     1  'SUPSC:x,y,ht,s1,n1,s2,n2: - n1,n2  neg = font 4'
        write(6,*)
     1  'SUBSUP:x,y,ht,s1,n1,s2,n2,s3,n3: - n1,n2,n3  neg = font 4'
        write(6,*)
     1  'AXES3:x0,y0,xl,yl,zl,ang1,ang2,ang3,t1,s1,t2,s2,',
     2  't3,s3,cmd'
        write(6,*)
     1  'AXES3SU:x0,y0,xl,yl,zl,ang1,ang2,ang3,t1,s1,u1,t2,s2,u2,',
     2  't3,s3,u3,cmd'
        write(6,*)
     1  'ARRAY: enter quoted file name containing x,y pairs '
        
        return
        end

        subroutine doxy(fname)
        character fname*(*)
        open(2,file=fname,access='sequential',form='formatted',
     1      status='unknown')
        ipen = 3
 1000   continue
            read(2,*,end=9999,err=9999)x,y
            call plot(x,y,ipen)
            ipen = 2
        go to 1000
 9999   continue
        call where(xx,yy,fct)
c-----
c       clean up
c-----
        call plot(xx,yy,3)
        rewind 2
        close(2)
        return
        end
