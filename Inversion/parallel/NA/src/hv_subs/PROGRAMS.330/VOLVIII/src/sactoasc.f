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
        real*4 x(LN)
c-----
c      name for SAC binary data file
c-----
        character name1*80
        character name2*80
        nmarg = mnmarg()
        if(nmarg .ne. 2)call usage()
        call mgtarg(1,name1)
        call mgtarg(2,name2)
c-----
c  Read input waveform data in SAC binary format.
c-----
            call brsac (1,LN,name1,x,nerr)
            if(nerr.ne.-1)call awsac (1,LN,name2,x)
        end


        subroutine usage()
        parameter (LER=0)
        write(LER,*)'Usage: sactoasc  SAC_BINARY_FILE SAC_ASCII'
        write(LER,*)' Convert SAC BINARY FILE TO SAC ASCII'
        stop 
        end

