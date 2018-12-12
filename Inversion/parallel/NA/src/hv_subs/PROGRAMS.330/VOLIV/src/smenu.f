        subroutine menu()
        integer LOT, LER, LOUT
        parameter (LOT=6, LER=0, LIN=5)
        integer NMENU
        parameter (NMENU=52)
        integer doprnt(NMENU)
c-----
c       LIST MENU OPTIONS
c       Output only non-blank entries. This would be
c       much easier to manage in C
c-----
        integer i,j,nprnt,nout,ll
        integer lnstr
        character*40 out(52)
        character*80 outl
        character outstr*12
        data outstr/'SURF96 MENU'/
        data (out(i),i=1,10)/
     1 ' 0- Display menu                        '
     1,' 1- Get  Dispersion                     '
     1,' 2- Invert Dispersion Data              '
     1,' 3- Invert Gamma Data                   '
     1,' 4- Simultaneous Gamma/Disp  Inversion  '
     1,' 5- Set Thick(0)/Velocity(1) Inversion  '
     1,' 6- Update Model (need lam)             '
     1,' 7- Plot RFTN/Dispersion/Velocity Model '
     1,' 8- Plot Gamma / Qb inverse  Model '
     1,' 9- Plot Resolution Kernel         '/
        data (out(i),i=11,20)/
     1 '10- List Singular Values                '
     1,'11- Gamma Data Partial Derivatives      '
     1,'12- Gamma Data Dispersion               '
     1,'13- List Q Model (need lam)             '
     1,'14- Q(beta) Resolving Kernel (need lam) '
     1,'                                        '
     1,'16- Velocity Data Partial Derivatives   '
     1,'17- Velocity Data Dispersion            '
     1,'18- List Velocity Model (need lam)      '
     1,'19- Velocity Resolving Kernels(need lam)'/
        data (out(i),i=21,30)/
     1 '22- ASCII Gamma Dispersion (file name)  '
     1,'23- ASCII Q(beta) File (file name,lam)  '
     1,'24- ASCII Q Resolving  (file name,lam)  '
     1,'                                        '
     1,'27- ASCII Veloc Dispersion (file name)  '
     1,'28- ASCII Model File  (file name, lam)  '
     1,'29- ASCII Vel Resolving(file name,lam)  '
     1,'30- (0) Fix Vp,(1) Fix Vp/Vs            ' 
     1,'31- Change dd(i), enter i,dd(i)         '
     1,'32- Enter Damping Factor (lam)          '/
        data (out(i),i=31,40)/
     1 '                                        '
     1,'                                        ' 
     1,'35- Inversion: (0) Non-Causal (default) '
     1,'               (1) Decoupled Causal     '
     1,'               (2) Coupled Causal       '
     1,'36- Smoothing: (0) Global reset none    '
     1,'               (1) Global reset diff    '
     1,'                                        '
     1,'37- Reset Number of Iterations          '
     1,'38- Temporary End                       '/
        data (out(i),i=41,50)/
     1 '39- Permanent End                       '
     1,'40- Enter Sigv minimum                  '
     1,'41- Enter Sigg minimum                  '
     1,'                                        '
     1,'                                        ' 
     1,'                                        '
     1,'45- Show Velocity Weights               ' 
     1,'46- Show Qinv     Weights               ' 
     1,'47- Show Inversion Controls             ' 
     1,'48- Modify Individual Layer Smoothing   '/
        data (out(i),i=51,52)/
     1 '                                        '
     1,'                                        '/
c-----
c       determine the number of non blank entries
c-----
        nprnt = 0
        do 50 j=1,NMENU
            if(out(j)(1:20).ne.'                    ')then
                nprnt = nprnt + 1
                doprnt(nprnt) = j
            endif
   50   continue
c-----
c       determine the halfway point in the output list
c-----
        if(mod(nprnt,2).eq.1)then
            nout = nprnt/2 + 1
        else
            nout = nprnt/2
        endif
c-----
c       output the menu
c-----
        write(LOT,1)outstr
    1   format(36x,a)
        do 100 i=1,nout
            outl = ' '
            outl(1:40)=out(doprnt(i))
            if(i.eq.nout)then
                if(mod(nprnt,2).eq.0)then
                outl(41:80)=out(doprnt(i+nout))
                endif
            else
                outl(41:80)=out(doprnt(i+nout))
            endif
            ll = lnstr(outl)
            write(LOT,'(1x,a)')outl(1:ll)
  100   continue
        write(LOT,*)'Enter Command at READY Prompt'
        return
        end
