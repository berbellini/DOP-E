c-----
c  INTERNAL PARAMETERS
c
c  rhdr    = SAC waveform real header fields (70).
c  ihdr    = SAC waveform integer header fields (40).
c  chdr    = SAC waveform character header fields (48).
c  x(i)    = array containing the waveform data.
c  delta   = waveform sampling interval in seconds.
c  btime   = waveform begin time.
c-----
        common/sachdr/rhdr,ihdr,chdr
        real*4 rhdr(70)
        integer*4 ihdr(40)
        character*8 chdr(24)
c-----
c      define storage space for trace data
c-----
        parameter (LN=131072)
        common/sacx/x
        real*4 x(LN)
        logical sacbin
c-----
c      name for SAC binary data file
c-----
        character name*80
        parameter (LOT=6)
c-----
c      parse command line arguments
c-----
        call gcmdln(sacbin,name)
        if(name .eq. ' ')stop
c-----
c  Read input waveform data in SAC binary format.
c-----
            if(sacbin)then
                call brsac (1,LN,name,x,nerr)
            else
                call arsac (1,LN,name,x,nerr)
            endif
            if(nerr.ne.-1)then
                call pinitf('SHWSAC.PLT')
                call parsesac(LN,x)
                call pend()
            endif
        end

        subroutine gcmdln(sacbin,fname)
        logical sacbin
        character fname*(*)

        character name*80
        nmarg = mnmarg()
        sacbin = .true.
        fname = ' '
        i = 0
   10   continue
            i = i + 1
            if(i.gt.nmarg)go to 11
            call mgtarg(i,name)
            if(name(1:2).eq.'-A')then
                sacbin = .false.
            else if(name(1:2).eq.'-B')then
                sacbin = .true.
            else if(name(1:2).eq.'-?')then
                call usage()
            else if(name(1:2).eq.'-h')then
                call usage()
            else 
                fname = name
            endif
        go to 10
   11   continue
        return
        end

        subroutine usage()
        parameter (LER=6)
        write(LER,*)'Usage: shwsac [-B] [-A] [-?] file_name'
        write(LER,*)
     1      ' -B    file is a binary SAC file (default)'
        write(LER,*)
     1      ' -A    file is an ASCII SAC file (default)'
        write(LER,*)
     1      ' -?    Command help - do not run'
        write(LER,*)
     1      ' -h    Command help - do not run'
        write(LER,*)
     1      ' file_name Name of SAC file (required)'
        stop
        end

        subroutine parsesac(LN,x)

        parameter(LIN=5, LOT=6, LER=0)

        real*4 rval
        integer*4 ival
        character*8 cval
        logical lval

        common/sachdr/rhdr,ihdr,chdr
        real*4 rhdr(70)
        integer*4 ihdr(40)
        character*8 chdr(24)
        real*4 x(LN)
        character*8 rstr(70), istr(40), cstr(24)
        data (rstr(i),i=1,45)/
     1  'DELTA'  , 'DEPMIN'  , 'DEPMAX'   , 'SCALE'   , 'ODELTA'   , 
     1  'B'       , 'E'      , 'O'        , 'A'       , 'FMT'      , 
     1  'T0'      , 'T1'      , 'T2'      , 'T3'      , 'T4'      , 
     1  'T5'      , 'T6'      , 'T7'      , 'T8'      , 'T9'      , 
     1  'F'       , 'RESP0'   , 'RESP1'   , 'RESP2'   , 'RESP3'   , 
     1  'RESP4'   , 'RESP5'   , 'RESP6'   , 'RESP7'   , 'RESP8'   , 
     1  'RESP9'   , 'STLA'    , 'STLO'    , 'STEL'    , 'STDP'    , 
     1  'EVLA'    , 'EVLO'    , 'EVEL'    , 'EVDP'    , 'MAG'  , 
     1  'USER0'   , 'USER1'   , 'USER2'   , 'USER3'   , 'USER4'    /
        data (rstr(i),i=46,70)/
     1  'USER5'   ,'USER6'   , 'USER7'   , 'USER8'   , 'USER9'   , 
     1  'DIST'    ,'AZ'      , 'BAZ'     , 'GCARC'   , 'SB'      , 
     1  'SDELTA'  ,'DEPMEN'  , 'CMPAZ'   , 'CMPINC'  , 'XMINIMUM', 
     1  'XMAXIMUM','YMINIMUM', 'YMAXIMUM', 'ADJTM'   , 'TIMMAX'  , 
     1  'TIMMIN'  ,'FHDR67'  , 'FHDR68'  , 'FHDR69'  , 'FHDR70'   /
        data istr/
     1  'NZYEAR'  ,'NZJDAY'  ,'NZHOUR'  ,'NZMIN'   ,'NZSEC'   , 
     1  'NZMSEC'  ,'NVHDR'   ,'NINF'    ,'NHST'    ,'NPTS'    , 
     1  'NSNPTS'  ,'NSN'     ,'NXSIZE'  ,'NYSIZE'  ,'NHDR15'  , 
     1  'IFTYPE'  ,'IDEP'    ,'IZTYPE'  ,'IHDR4'   ,'IINST'   , 
     1  'ISTREG'  ,'IEVREG'  ,'IEVTYP'  ,'IQUAL'   ,'ISYNTH'  , 
     1  'IHDR11'  ,'IHDR12'  ,'IHDR13'  ,'IHDR14'  ,'IHDR15'  , 
     1  'IHDR16'  ,'IHDR17'  ,'IHDR18'  ,'IHDR19'  ,'IHDR20'  , 
     1  'LEVEN'   ,'LPSPOL'  ,'LOVROK'  ,'LCALDA'  ,'LHDR5'    /
        data cstr/
     1  'KSTNM'   ,'KEVNM'   ,'KEVNMC'  ,'KHOLE'   , 
     1  'KO'      ,'KA'      ,'KT0'     ,'KT1'     , 
     1  'KT2'     ,'KT3'     ,'KT4'     ,'KT5'     , 
     1  'KT6'     ,'KT7'     ,'KT8'     ,'KT9'     , 
     1  'KF'      ,'KUSER0'  ,'KUSER1'  ,'KUSER2'  , 
     1  'KCMPNM'  ,'KNETWK'  ,'KDATRD'  ,'KINST'   
     1  /
c-----
c      formats
c-----
   10   format(' REAL        INDEX  NAME',
     1      '        INT VALUE   REAL VALUE')
   11   format(13x,i5,2x,a8,e13.5,e13.5)
   20   format('INTEGER     INDEX  NAME',
     1      '        INT VALUE    INT VALUE')
   21   format(13x,i5,2x,a8,i13,i13)
   30   format(' ENNUMERATED INDEX  NAME',
     1      '        INT VALUE    ENU VALUE')
   31   format(13x,i5,2x,a8,i13,5x,a8)
   40   format(' LOGICAL     INDEX  NAME',
     1      '        INT VALUE    LOG VALUE')
   41   format(13x,i5,2x,a8,i13,5x,l8)
   50   format(' CHARACTER   INDEX  NAME',
     1      '        INT VALUE   CHAR VALUE')
   51   format(13x,i5,2x,a8,5x,a8,5x,a8)
c-----
c      output real header
c-----
        write(LOT,10)
        do 1000 i=1,70
            call getfhv(rstr(i),rval,nerr)
            write(LOT,11)i,rstr(i),rhdr(i), rval
 1000   continue
c-----
c      output integer header
c-----
        write(LOT,20)
C       do 2000 i=1,15
        do 2000 i=1,40
            call getnhv(istr(i),ival,nerr)
            write(LOT,21)i,istr(i),ihdr(i), ival
 2000   continue
c-----
c      output enumerated header
c-----
        write(LOT,30)
        do 3000 i=16,25
            call getihv(istr(i),cval,nerr)
            write(LOT,31)i,istr(i),ihdr(i), cval
 3000   continue
c-----
c      output logical header
c-----
        write(LOT,40)
        do 4000 i=36,40
            call getlhv(istr(i),lval,nerr)
            write(LOT,41)i,istr(i),ihdr(i), lval
 4000   continue
c-----
c      output character header
c-----
        write(LOT,50)
        do 5000 i=1,24
            call getkhv(cstr(i),cval,nett)
            write(LOT,51)i,cstr(i),chdr(i), cval
 5000   continue
        call getnhv('NPTS    ',npts,nerr)
        inc = npts/10
        do 9000 i=1,npts,inc
            write(LOT,*)i,x(i)
 9000   continue
        write(LOT,*)npts,x(npts)
c-----
c      plot it
c-----
        npts = ihdr(10)
        call plotit(x,npts,rhdr(1),rhdr(2),rhdr(3))
        return
        end

        subroutine plotit(x,npts,rhdr1,rhdr2,rhdr3)
        real*4 x(npts)
C       ymax = abs(rhdr2)
C       if(ymax.lt.abs(rhdr3))ymax = abs(rhdr3)
        ymx = x(1)
        ymn = x(1)
        do 1001 i=2,npts
            if(x(i) .gt.ymx)ymx = x(i)
            if(x(i) .lt.ymn)ymn = x(i)
 1001   continue
        
        if(ymn .eq.0.0 .or.ymn.eq.ymx)then
            ymnv = -1.0
            ymxv =  1.0
        else
            ymnv = ymn
            ymxv = ymx  
        endif 
        if(ymx.eq.ymn)then
            ymn = -ymx
        endif
        if(npts.le.0)return
        dx = 8.0/npts
        ipen = 3
        xx = 1.0
        do 1000 i=1,npts
            xx = xx + dx
            yy = 3.5 + 2.0*(x(i) - ymnv)/(ymxv - ymnv)
            call plot(xx,yy,ipen)
            ipen= 2
 1000   continue
        call symbol(  1.0,0.6,0.10,'NPTS=',0.0,5)
        call number(999.0,0.6,0.10,real(npts),0.0,-1)
        call symbol(  3.0,0.6,0.10,' MAX=',0.0,5)
        call number(999.0,0.6,0.10,ymx,0.0,1003)
        call symbol(  1.0,0.4,0.10,'DT=',0.0,3)
        call number(999.0,0.4,0.10,rhdr1,0.0,1003)
        call symbol(  3.0,0.4,0.10,' MIN=',0.0,5)
        call number(999.0,0.4,0.10,ymn,0.0,1003)
        return
        end
