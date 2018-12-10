        program fmmfit
c----------------------------------------------------------------------c
c                                                                      c
c      COMPUTER PROGRAMS IN SEISMOLOGY                                 c
c      VOLUME V                                                        c
c                                                                      c
c      PROGRAM: fmdfit                                                 c
c                                                                      c
c      COPYRIGHT 2002                                                  c
c      R. B. Herrmann                                                  c
c      Department of Earth and Atmospheric Sciences                    c
c      Saint Louis University                                          c
c      221 North Grand Boulevard                                       c
c      St. Louis, Missouri 63103                                       c
c      U. S. A.                                                        c
c                                                                      c
c----------------------------------------------------------------------c
c Revision history:
c       31 JUL 2002 - programs created
c       21 APR 2004 - corrected type in subroutine usage
c
c       Plot focal mechanism grid search results from 
c       srdisp96 and ..... for a given depth
c----------------------------------------------------------------------c
        implicit none
        integer LER, LIN, LOT
        parameter (LIN=5,LOT=6,LER=0)                          

        real x0, xlen, y0, ylen
        
        real dipmn,dipmx, ddip, stkmn,stkmx, dstk
        integer jid
        character instr*170
        integer ls, lsep
        real h, stk, dip, rak, rr, rl, xmwr, xmwl, best
        logical doit
        integer lgstr
        real x6, x7
        character ostr*10

        call gcmdln(dipmn,dipmx,ddip,
     1      stkmn,stkmx,dstk,x0,y0,xlen,ylen,jid)

        if(jid.lt.0)then
            call pinitf('FMMFIT.PLT')
        else
            write(ostr,'(i2.2)')jid
            ls = lgstr(ostr)
            call pinitf('FMMFIT'//ostr(1:ls)//'.PLT')
        endif
c-----
c       draw a box
c-----
        call gbox(x0,y0,x0+xlen,y0+ylen)
c-----
c       put in the axes
c-----
        call dolinx(x0,y0,xlen,stkmn,stkmx,
     1      0.10,.true.,.false.,.true., 6, 'Strike')
        call dolinx(x0,y0+ylen,xlen,stkmn,stkmx,
     1      0.10,.false.,.false.,.false., 6, 'Strike')
        call doliny(x0,y0,ylen,dipmn,dipmx,
     1      0.10,.false.,.true.,.true.,3,'Dip')
        call doliny(x0+xlen,y0,ylen,dipmn,dipmx,
     1      0.10,.true.,.true.,.false.,3,'Dip')
c-----
c       read data and plot only within the search bounds
c-----
 1000   continue
        read(5,'(a)',end=9000)instr
c-----
c       check format
c-----
            ls = lgstr(instr)
            lsep = index(instr,'SRFGRD96')
            doit = .false.
            if(lsep.gt.0 .and.ls.gt.8)then
                read(instr(lsep+8:ls),*)h,stk,dip,rak,rr,rl,
     1              xmwr,xmwl,best
                doit = .true.
            endif
            lsep = index(instr,'WVFGRD96')
            if(lsep.gt.0 .and.ls.gt.8)then
                read(instr(lsep+8:ls),*)h,stk,dip,rak,xmwr,best
                doit = .true.
            endif
            lsep = index(instr,'WVFMTD96')
            if(lsep.gt.0 .and.ls.gt.8)then
                read(instr(lsep+8:ls),*)h,stk,dip,rak,xmwr,x6,x7,best
                doit = .true.
            endif
            lsep = index(instr,'WVFMT96')
            if(lsep.gt.0 .and.ls.gt.7)then
                read(instr(lsep+7:ls),*)h,stk,dip,rak,xmwr,x6,x7,best
                doit = .true.
            endif
c-----
c       if have valid entry, plot it
c-----
            if(doit .and. best.gt.0.0)then
                call plotit(stk,dip,rak,best,dipmn,dipmx,ddip,
     1              stkmn,stkmx,dstk,x0,y0,xlen,ylen)
            endif
        go to 1000
 9000   continue
        call pend()
        end

        subroutine plotit(stk,dip,rak,best,dipmn,dipmx,ddip,
     1      stkmn,stkmx,dstk,x0,y0,xlen,ylen)
c-----
c       stk R   - strike
c       dip R   - dip
c       rak R   - rake
c       best    R   - best value [-1,1]
c       dipmn   R   - minimum dip value for y-axis
c       dipmx   R   - maximum dip value for y-axis
c       stkmn   R   - minimum strike value for x-axis
c       stkmx   R   - maximum strike value for x-axis
c       x0  R   - x-coordinate of lower left corner
c       y0  R   - y-coordinate of lower left corner
c       xlen    R   - length of x-axis
c       ylen    R   - length of y-axis
c-----
        implicit none
        real stk,dip,rak,best,dipmn,dipmx,ddip,dstk
        real stkmn,stkmx,x0,y0,xlen,ylen
        real xx, yy
        integer jpen
        real rakrad, c, s
        real arlen, t, u
        if(dip.lt.dipmn .or. dip.gt.dipmx .or. stk.lt.stkmn
     1      .or. stk.gt.stkmx)return
c-----
c       map best [0,1] to colors [BLUE,RED] = [1100, 1000]
c-----
        jpen = 1100 - 100.0*best
        xx = x0 + xlen * (stk - stkmn)/(stkmx - stkmn)
        yy = y0 + ylen * (dip - dipmn)/(dipmx - dipmn)
        rakrad = rak*3.1415927/180.0
        c = cos(rakrad)
        s = sin(rakrad)
        call newpen(jpen)
        call plot(xx,yy,3)
c-----
c       compute the arrow length
c-----
        t = ylen*ddip/(dipmx-dipmn)
        u = xlen*dstk/(stkmx-stkmn)
        if(t.gt.u)then
            arlen = u
        else
            arlen = t
        endif

        xx = xx + best*arlen*s
        yy = yy + best*arlen*c
        call plot(xx,yy,2)
        return
        end



        subroutine gcmdln(dipmn,dipmx,ddip,
     1      stkmn,stkmx,dstk,x0,y0,xlen,ylen,jid)
c-----
c       parse command line arguments
c
c       requires subroutine mgtarg() and function mnmarg()
c
c-----
c       dipmn,dipmx,ddip    R*4 - dip range search parameters
c       stkmn,stkmx,dstk    R*4 - strike range search parameters
c       x0,y0           R*4 - lower left corner of plot
c       xlen,ylen       R*4 - length of X and Y axes
c       jid         I*4 - file ID for PLOT
c-----
        implicit none
        integer LER, LIN, LOT
        parameter (LER=0, LIN=5, LOT=6)
c-----
c       parameters for grid search
c-----
        real dipmn, dipmx, ddip
        real stkmn, stkmx, dstk
        real x0, y0, xlen, ylen
        integer jid

        character*25 names
        integer*4 mnmarg
        integer nmarg
        integer i
c-----
c       initialize variables
c-----
        dipmn = 30  
        dipmx = 90
        ddip = 15
        stkmn = 0
        stkmx = 360
        dstk = 10
        x0 = 1.0
        y0 = 1.0
        xlen = 8.0
        ylen = 6.0
        jid = -1
c-----
c       process command line arguments
c-----
        nmarg=mnmarg()
        i = 0
   11   i = i + 1
        if(i.gt.nmarg)go to 13
            call mgtarg(i,names)
            if(names(1:4).eq.'-DMN')then
                i=i+1
                call mgtarg(i,names)
                call chtofp(names,dipmn)
            else if(names(1:4).eq.'-DMX')then
                i=i+1
                call mgtarg(i,names)
                call chtofp(names,dipmx)
            else if(names(1:3).eq.'-DD')then
                i=i+1
                call mgtarg(i,names)
                call chtofp(names,ddip)
            else if(names(1:4).eq.'-SMN')then
                i=i+1
                call mgtarg(i,names)
                call chtofp(names,stkmn)
            else if(names(1:4).eq.'-SMX')then
                i=i+1
                call mgtarg(i,names)
                call chtofp(names,stkmx)
            else if(names(1:3).eq.'-DS')then
                i=i+1
                call mgtarg(i,names)
                call chtofp(names,dstk)
            else if(names(1:3).eq.'-ID')then
                i=i+1
                call mgtarg(i,names)
                read(names,'(bn,i10)')jid
            else if(names(1:2).eq.'-?')then
                call usage(' ')
            else if(names(1:2).eq.'-h')then
                call usage(' ')
            endif
        goto 11
   13   continue
        return
        end
        
        subroutine usage(ostr)
        implicit none
        character ostr*(*)
        integer LER, LIN, LOT
        parameter (LER=0, LIN=5, LOT=6)
        integer ls
        integer lgstr
        ls = lgstr(ostr)
        if(ostr.ne.' ') write(LER,*)ostr(1:ls)
        write(LER,*)'fmmfit -DMN dipmn -DMX dipmx -DD ddip'
        write(LER,*)'         -SMN stkmn -SMX stkmx -DS dstk -ID jid'
        write(LER,*)
     1  ' -DMN dipmn  (default   30)  Minimum dip'
        write(LER,*)
     1  ' -DMX dipmx  (default   90)  Maximum dip'
        write(LER,*)
     1  ' -DD ddip    (default   15)  dip increment'
        write(LER,*)
     1  ' -SMN stkmn  (default    0)  Minimum strike'
        write(LER,*)
     1  ' -SMX stkmx  (default  350)  Maximum strike'
        write(LER,*)
     1  ' -DS dstk    (default   10)  strike increment'
        write(LER,*)
     1  ' -ID jid     (default   0)  Integer ID for naming plot'
        write(LER,*)
     1  ' -?                   This online help'
        write(LER,*)
     1  ' -h                   This online help'
        stop
        end

        subroutine chtofp(str,fout)
c------
c       routine to convert string to floating point
c       The E format is accepted as well as free form
c       input
c
c       This assumes that the string is 20 characters long
c-----
        implicit none
        character*(*) str
        real*4 fout
        integer*4 lgstr
        integer i, l
        logical hase

        l = lgstr(str)
c------
c       If the string str contains an E or e, then
c       we parse the result using an E format.
c------
        hase = .false.
        do 100 i=1,l
            if(str(i:i) .eq.'e' .or. str(i:i).eq.'E')then
                hase = .true.
            endif
  100   continue
c-----
c       read floating point number
c-----
        if(hase)then
            read(str,'(bn,e20.13)')fout
        else
            read(str,'(bn,f20.0)') fout
        endif
        return
        end
