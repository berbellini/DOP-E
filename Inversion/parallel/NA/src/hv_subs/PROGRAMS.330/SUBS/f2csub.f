        real*8 function cdabs(Z)
        complex*16 Z
        cdabs = zabs(Z)
        return
        end

        complex*16 function cdexp(Z)
        complex*16 Z
        cdexp = zexp(Z)
        return
        end

        complex*16 function cdlog(Z)
        complex*16 Z
        cdlog = zlog(Z)
        return
        end

        complex*16 function cdsqrt(Z)
        complex*16 Z
        cdsqrt = zsqrt(Z)
        return
        end

        real*8 function dreal(Z)
        complex*16 Z
        dreal = real(Z)
        return
        end
        
