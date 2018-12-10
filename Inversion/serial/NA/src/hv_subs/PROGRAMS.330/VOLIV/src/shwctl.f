
        subroutine shwctl()
c-----
c       display all processing parameters
c       in a manner that facilitates assigning new values
c
c       For simplicity just read it from the file
c       instead of passing all arguments through the
c       command line
c-----
        integer LER, LIN, LOT
        parameter(LER=0,LIN=5,LOT=6)

        integer NL
        parameter (NL=200)
        integer nf10(NL)
        data nf10/NL*1/

        call gttmp0(iprog,itot,nf1,nf2,nf34,nf5,nf67,nf10,nfilt,
     1      nup,dlam,qaqb,wref,ivcs,
     2      lsin,twmn,twmx,iter,nurf,indp,pval,
     4              sigv,sigr,sigg,
     3      id2 ,idum2,idum3,idum4,idum5,
     4      rdum1,rdum2,rdum3,rdum4,rdum5)
        if(iprog.eq.1)then
            write(LOT,'(a)')'Inversion controls for surf96'
        else if(iprog.eq.2)then
            write(LOT,'(a)')'Inversion controls for rftn'
        else if(iprog.eq.3)then
            write(LOT,'(a)')'Inversion controls for joint96'
        endif
    1   format('Cmd   Value   Description')
    2   format(i3,3x,i3,6x,a)
    3   format(i3,f8.3,4x,a)
    4   format(3x,3x,i3,6x,a)
    5   format(i3,1x,f8.4,2x,a)
    6   format(15x,a)
    7   format(i3,g11.3,1x,a)
        write(LOT,1)
        WRITE(LOT,4)nf1,' 1 Variance based residual of fit'
        WRITE(LOT,6)    ' 0 Variance based on observed std observation'
        if(iprog.eq.1 .or. iprog.eq.3)then
        WRITE(LOT,4)nf34,' Maximum number of Love     modes to use'
        WRITE(LOT,4)nf67,' Maximum number of Rayleigh modes to use'
        endif
        nf8 = mod(nfilt,2)
        nf11= nfilt/2
        WRITE(LOT,4)iter,' Current iteration'
        WRITE(LOT,4)nurf,' Number of receiver functions to be inverted'
        WRITE(LOT,4)lsin,' 2 last inversion for Vs'
        WRITE(LOT,6)     ' 3 last inversion for Q inverse'
        WRITE(LOT,6)     ' 4 last inversion for Vs-Q inverse'
        WRITE(LOT,2) 5,indp,' 0 Layer thickness  inversion'
        WRITE(LOT,6)        ' 1 Layer velocity/Q inversion'
        WRITE(LOT,3)32,dlam,' Damping value    (default value 1.0)'
        if(iprog.eq.2 .or. iprog.eq.3)then
        WRITE(LOT,3)33,twmn,' Minimum window for RFTN (default -5.0 s)'
        WRITE(LOT,3)34,twmx,' Maximum window for RFTN (default 20.0 s)'
        endif
        WRITE(LOT,2)35,ivcs,' 0 non-causal Vs - Q relation (default)'
        WRITE(LOT,6)        ' 1 Decoupled causal'
        WRITE(LOT,6)        ' 2 Fully coupled causal'
        WRITE(LOT,2)36,nf11,' 0 No smoothing constraint'
        WRITE(LOT,6)        ' 1 Differential smoothing constraint'
        if(iprog.eq.1 .or. iprog.eq.3)then
        WRITE(LOT,5)40,sigv,' Std error of fit floor for velocity disp'
        endif
        if(iprog.eq.1 )then
        WRITE(LOT,7)41,sigg,' Std error of fit floor for gamma disp'
        endif
        if(iprog.eq.2 .or. iprog.eq.3)then
        WRITE(LOT,5)42,sigr,' Std error of fit floor for RFTN '
        endif
        if(iprog.eq.3)then
        WRITE(LOT,3)43,pval,' Joint inversion influence parameter'
        WRITE(LOT,6)        ' RFTN = 0 <= p <= 1 Surface Wave '
        endif
        if(iprog.eq.2 .or. iprog.eq.3)then
        WRITE(LOT,2)44,id2 ,' 0 Match observed RFTN window'
        WRITE(LOT,6)        ' 1 Use 2x times series to compute RFTN'
        endif

   23   format('----------------------------------------'   ,
     1  '----------------------------------------'   )
        WRITE(LOT,23)
        WRITE(LOT,*)'Use menu command cmd to change value'


        return
        end
