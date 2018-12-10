
        subroutine fquery(str,fval,lquiet,i)
c-----
c       output string, pronpt for query
c-----
c       str Ch* output string
c       fval    R   float value
c       lquiet  L   .true. silent mode
c-----
        character str*(*)
        real fval
        logical lquiet
        character name*80
        integer LIN, LER, LOT
        parameter (LIN=5, LOT=6, LER=0)
        
        if(lquiet)then
            i = i + 1
            call mgtarg(i,name)
            read(name,'(bn,f10.0)')fval
        else
            write(LOT,*)str,fval
 1000       continue
            write(LOT,*)'Enter new float value'
            read(LIN,*,err=1000,end=1000)fval
        endif
        return
        end
        
        subroutine iquery(str,ival,lquiet,i)
c-----
c       output string, pronpt for query
c-----
c       str Ch* output string
c       ival    I   integer value
c       lquiet  L   .true. silent mode
c-----
        character str*(*)
        logical lquiet
        integer ival, i
        character name*80
        integer LIN, LER, LOT
        parameter (LIN=5, LOT=6, LER=0)
        
        if(lquiet)then
            i = i + 1
            call mgtarg(i,name)
            read(name,'(bn,i10)')ival
        else
            write(LOT,*)str,ival
 1000       continue
            write(LOT,*)'Enter new integer value'
            read(LIN,*,err=1000,end=1000)ival
        endif
        return
        end
        
        subroutine ibound(ilw,iup,ival,idef)
        implicit none
        integer ilw,iup,ival,idef
        integer LIN, LER, LOT
        parameter (LIN=5, LOT=6, LER=0)
c-----
c       Perform a bounds check on the integer ival
c
c       ival must be     ilw <= ival <= iup
c
c       otherwise  reset ival to the default value idef
c-----
        if(ival.lt.ilw .or. ival.gt.iup)then
            ival = idef
            WRITE(LOT,*)'Outside bounds. Setting to default:',
     1          idef
        endif
        return
        end

        subroutine fbound(flw,fup,fval,fdef)
        implicit none
        real flw, fup, fval, fdef
        integer LIN, LER, LOT
        parameter (LIN=5, LOT=6, LER=0)
c-----
c       Perform a bounds check on the real fval
c
c       ival must be     flw <= fval <= fup
c
c       otherwise  reset ival to the default value fdef
c-----
        if(fval.lt.flw .or. fval.gt.fup)then
            fval = fdef
            WRITE(LOT,*)'Outside bounds. Setting to default:',
     1          fdef
        endif
        return
        end

