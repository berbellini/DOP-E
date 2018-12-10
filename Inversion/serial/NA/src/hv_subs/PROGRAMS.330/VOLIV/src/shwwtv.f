        subroutine shwwtv(mmax,dd,wc,nf10)
c-----
c       show the layer weighting for Velocity or Thickness 
c
c       mmax    I   - number of layers
c       dd  R   - array of weights
c       wc  L   - array of smoothing controls
c                 .true.  smoothing
c                 .false. no smoothing
c       nf10    I   - array indicating if VP/VS or VP is fixed
c-----
        integer LER, LIN, LOT
        parameter(LER=0,LIN=5,LOT=6)

        integer NL, NL2, NLAY
        parameter(NL=200,NLAY=200,NL2=NL+NL)

        real dd(NL2)
        logical wc(NL2)
        integer nf10(NL)

        integer mmax

        common/isomod/dl(NLAY),va(NLAY),vb(NLAY),rho(NLAY),
     1      qa(NLAY),qb(NLAY),etap(NLAY),etas(NLAY), 
     2      frefp(NLAY), frefs(NLAY)
        common/depref/refdep
        common/modtit/title
        character*80 title

        real dtop(NLAY), dbot(NLAY)
        character*4 cstr(2)
        character*4 ostr(NLAY)
        character str*80
        cstr(1)=' Lyr'
        cstr(2)=' Bdy'
c-----
c       formats
c-----
c       formats for layer weighting
c-----
   21   format( i3,1x,f7.2,a4,1x,i1,f10.2,'-',f10.2)
c-----
c       formats for boundary weighting
c-----
   31   format( i3,1x,f7.2,a4,1x,i1,f10.2,1x,10x)
c-----
c       title formats
c-----
   19   format( 'V  Model Weighting Parameters: Large value',
     1  ' forces change at boundary or layer   ')
   20   format( '  I   DD(I) Inv S     depth/range        ',
     1          '  I   DD(I) Inv S     depth/range')
   23   format('----------------------------------------'   ,
     1  '----------------------------------------'   )
   

        depth = 0.0
        dtop(1) = 0.0
        do 1000 i=1,mmax-1
            dbot(i)   = dtop(i) + dl(i)
            dtop(i+1) = dbot(i)
            if(wc(i))then
                ostr(i) = cstr(2)
            else
                ostr(i) = cstr(1)
            endif
 1000   continue
        dbot(mmax) = 9999.
        if(wc(mmax))then
            ostr(mmax) = cstr(2)
        else
            ostr(mmax) = cstr(1)
        endif

            write(LOT,19)
            write(LOT,20)
            write(LOT,23)
            if(mod(mmax,2).eq.0)then
                mup = mmax/2
            else
                mup = mmax/2 + 1
            endif

            do 2000 ilw=1,mup
                str = ' '
                iup = ilw + mup
                j=ilw
                k=iup
                jq = j 
                kq = k 
                if(wc(jq))then
                    write(str(1:38),31)jq,dd(jq),ostr(jq),
     1                  nf10(jq),dbot(j)
                else
                    write(str(1:38),21)jq,dd(jq),ostr(jq),
     1                  nf10(jq),dtop(j),dbot(j)
                endif
                str(39:43) = '  |  '
                if(kq.le.mmax)then
                if(wc(kq))then
                    write(str(42:79),31)kq,dd(kq),ostr(kq),
     1                  nf10(kq),dbot(k)
                else
                    write(str(42:79),21)kq,dd(kq),ostr(kq),
     1                  nf10(kq),dtop(k),dbot(k)
                endif
                endif
                write(LOT,'(a)')str
 2000       continue
            write(LOT,23)
            WRITE(LOT,*)'Bdy - get velocity change at boundary'
            WRITE(LOT,*)'Lyr - get velocity in layer'
            WRITE(LOT,*)'S =0 Vp fixed, S=1 Vp/Vs fixed in layer'
            WRITE(LOT,*)'Use option 30 to change how Vp obtained'
            WRITE(LOT,*)'Use option 31 to change layer weight'
            WRITE(LOT,*)'Use option 48 to change layer smoothing'
            WRITE(LOT,*)'Use option 45 to redisplay this menu'
        return
        end
