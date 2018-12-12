      program ffilt96 
c---------------------------------------------------------------------c 
c                                                                     c 
c      COMPUTER PROGRAMS IN SEISMOLOGY                                c 
c      VOLUME V                                                       c 
c                                                                     c
c      PROGRAM: FFILT96                                               c
c                                                                     c
c      COPYRIGHT 1996                                                 c
c      R. B. Herrmann                                                 c
c      Department of Earth and Atmospheric Sciences                   c
c      Saint Louis University                                         c
c      221 North Grand Boulevard                                      c
c      St. Louis, Missouri 63103                                      c
c      U. S. A.                                                       c
c                                                                     c
c---------------------------------------------------------------------c
c       Changes
c
c       17 JUL 2000 modified Seoul korea to permit robustness in reading 
c           the SAC response file. A * in column 1 to indicate a
c           comment works, a blank or empty line is not assumed to
c           be a 0+0i entry, number of zeros do not have to be
c           specified if they are zero
c       11 SEP 2000 add TP, TSV TSH to FILE96 format
c       14 JAN 2001 add A, C, F, L, N and density information
c       23 JUL 2009 - ensured consistent call to FFT routine with
c               complex arguments. We now call zfour
c-----
c
c     ffilt96 -DEMEAN -TAPER -PCT pct -FL fl -FH fh -A -R 
c           -PZ pole_zero_file -W water
c
c     This program applies or removes a filter response for
c     a SAC pole/zero format. The input and output
c     are  file96(V) time series files
c
C       HOWEVER the SAC pole_zero file must be complete, SAC permits
c       repeated zeros to be ignored
c-----
c-----
c-----
c-----
c       variables in common blocks
c
c       iftype  I*4 File type
c               1 - single trace
c               3 - three component
c               16 - Green's function
c               21 - Green's function
c
c       iobsyn  I*4 1 - observed
c               2 - synthetic
c
c       itmfrq  I*4 1 - time series
c               2 - Fourier spectra (not implemented)
c
c       iunit   I*4 1   - counts 
c               2   - cm
c               3   - cm/sec
c               4   - cm/sec/sec
c               5   - m
c               6   - m/sec
c               7   - m/sec/sec
c               8   - microns
c               9   - microns/sec
c               10  - microns/sec/sec
c       junit   I*4
c               11  - Pa  (nt/m^2)
c               12  - MPa  (mega Pascal)
c
c       cfilt   C*80    comment on filtering operations
c
c       keyear  I*4 event year
c       kemon   I*4 event mon
c       keday   I*4 event day
c       kehour  I*4 event hour
c       kemin   I*4 event minute
c       esec    R*4 event second
c       evlat   R*4 event latitude
c       evlon   R*4 event longitude
c       evdep   R*4 event depth
c       
c       stname  C*8 station name
c               NOTE: AH(6), SEED(5), CSS(6), SAC(8)
c       
c       stlat   R*4 station latitude
c       stlon   R*4 station longitude
c       stelev  R*4 station elevation
c
c
c       distkm  R*4 epicentral distance in km
c       distdg  R*4 epicentral distance in degrees (111.195 km/degree)
c       evstaz  R*4 event -> station azimuth, degrees east of north
c       stevaz  R*4 station -> event azimuth, degrees east of north
c       
c       cpulse  C*80    pulse description
c
c       ccomnt  C*80    comment
c
c       jsrc    I*4 Array of indices to indicate if traces are 
c                   present
c               iftype =  1  jsrc(1) indicates if trace is 
c                   present
c               iftype =  3  jsrc(i),i=1,5 indicates if trace 
c                   is present
c                   i = Z, 2=N, 3=E, 4=R, 5=T
c                   (we could carry all five, but the real 
c                   purpose  is to permit no more than 
c                   3 traces, but 2 may all that are 
c                   present in a real data set)
c
c-----

        common/ihdr96/ ftype, iftype, iobsyn, itmfrq, iunit, idecon,
     1      keyear, kemon, keday, kehour, kemin,
     2      ksyear, ksmon, ksday, kshour, ksmin,
     3      jsrc, junit
        integer*4 ftype, iftype, iobsyn, itmfrq, iunit, idecon
        integer*4 keyear, kemon, keday, kehour, kemin
        integer*4 ksyear, ksmon, ksday, kshour, ksmin
        integer*4 jsrc(21), junit

        common/rhdr96/ esec, distkm, distdg, evstaz, stevaz,
     1      evlat, evlon, evdep, stlat, stlon, stelev,
     2      TP, TSV, TSH, SA, SC, SF, SL, SN, SR
        real*4 esec, distkm, distdg, evstaz, stevaz
        real*4 evlat, evlon, evdep 
        real*4 stlat, stlon, stelev
        REAL*4 TP, TSV, TSH
        REAL*4 SA, SC, SF, SL, SN, SR

        common/chdr96/ stname, cfilt, cpulse, ccomnt
        character*8 stname*8, cfilt*80, cpulse*80, ccomnt*80    

        character stcomp*8
        real*4 cmpaz, cmpinc, cmpdt
        real*4 ssec
c-----
c       Pole zero description of complex transfer function
c
c       H(s) = const Product [ s - zero(i) ]/ Product [s - pole(i) ]
c                  
c-----
c       zpoles      Cmp*8   - complex poles
c       zzeros      Cmp*8   - complex zeros
c       npoles      I*4 - number of poles
c       nzeros      I*4 - number of zeros
c       const       R*4 - constant
c       npordr      I*4 - 0 do not use, 1 - >
c                       order in increasing angular frequency
c                       for single or double poles, separately
c       nptyp       I*4 - 1 = single pole, 2=double
c                     2 = part of complex conjucgate pair
c                     0 = do not use, part of complex conjugate pair
c       npardr      I*4 - 0 do not use, 1 - >
c                       order in increasing angular frequency
c                       for all poles, separately
c       pwn     R*4 - natural frequency
c       phn     R*4 - damping for second order
c       pconst      R*4 - modified constant accounting for low pass
c-----
        parameter (NPZ=100)
        common/pzresp/zpoles(NPZ), zzeros(NPZ), npoles, nzeros, 
     1      const, npordr(NPZ), nptyp(NPZ), npardr(NPZ),
     2      pwn(NPZ), phn(NPZ), pconst
        complex zpoles, zzeros
        integer*4 npoles, nzeros, npordr, nptyp, npardr
        real*4 const, pwn, phn, pconst
c-----
c       command line arguments
c-----
        logical domean, dotape, doappl
        real*4 fl, fh, pct
        character fname*80
c-----
c       variables local to the program
c-----
        parameter (LER=6, LIN=5, LOT=6)
        integer NSAMP
        parameter (NSAMP=16384)
        integer*4 npts
        real*4 x(NSAMP)
c-----
c       call machine dependent initialization
c-----
        call mchdep()
c-----
c       parse command line arguments
c-----
        call gcmdln(domean, dotape, doappl, fl, fh, pct, fname,water)
c-----
c       get the sac response
c-----
        call gtsrsp(fname,ierr)
        if(ierr.lt.0)go to 9999
c-----
c       process
c-----
 1000   continue
            call rdhd96(LIN,nerr)
            if(nerr .lt. 0 )go to 9999
c-----
c           modify the comment string
c-----
            call wrhd96(LOT,nerr)
c-----
c           process the traces
c-----
            do 200 jgrn=1,21
c-----
c           initialize input array  
c-----
                if(jsrc(jgrn).ne.0)then
                    call rdtr96(LIN,stcomp,cmpinc,cmpaz,
     1                  cmpdt, npts, ksyear, ksmon, 
     2                  ksday, kshour, ksmin, ssec, 
     3                  x,nerr,NSAMP)
                    if(nerr .lt. 0 )go to 9999
c-----
c               the same, do not recompute the pulse
c-----
c-----
c               set up filter coefficients which 
c               depend upon dt
c-----
                    dt = cmpdt
                    fnyq = 0.5/dt

                    if(domean)call dmean(x,npts)
                    if(dotape)call dtape(x,npts,pct)
                    call filt(doappl,x,npts,dt,
     1                  water,fl,fh)

c-----
c                   output the filtered trace
c-----
                    call wrtr96(LOT,stcomp,cmpinc,cmpaz,
     1                  cmpdt, npts, ksyear, ksmon, 
     2                  ksday, kshour, ksmin, ssec, 
     3                  x,nerr,NSAMP)
                endif
  200       continue
        go to 1000
 9999   continue
        end

        subroutine gcmdln(domean, dotape, doappl, fl, fh, 
     1      pct, fname,water)
c-----
c       
c     ffilt96 -DEMEAN -TAPER -PCT pct -FL fl -FH fh -A -R 
c           -PZ pole_zero_file -W water
c
c       parse command line arguments and return control
c       parameters
c
c       domean  L   - .true. remove mean from time seeries
c       dotape  L   - .true. cosine taper ends
c       doappl  L   - .true. apply filter
c                 .false. remove filter
c       fl  R*4 - low cut corner of filter operation
c       fh  R*4 - high cut corner of filter operation
c       pct R*4 - percent of taper.
c       fname   C*80    - name of SAC response file
c       water   R*4 - water level for inversion, fraction of maximum
c
c-----
        logical domean, dotape, doappl
        real*4 fl, fh, pct, water
        character fname*80
        
        character*20 name
        integer*4 nmarg
        nmarg = mnmarg()
        domean = .false.
        dotape = .false.
        doappl = .true.
        fl = -1.0
        fh = -1.0
        pct = 5.0
        water = 0.01
        fname = ' '
        i = 0
   11   i = i + 1
            if(i.gt.nmarg)goto 13
            call mgtarg(i,name)
            if(name(1:3).eq.'-DE')then
                domean = .true.
            else if(name(1:3).eq.'-TA')then
                dotape = .true.
            else if(name(1:2).eq.'-A')then
                doappl = .true.
            else if(name(1:2).eq.'-R')then
                doappl = .false.
            else if(name(1:3).eq.'-FL')then
                i = i + 1
                call mgtarg(i,name)
                read(name,'(bn,f10.0)')fl
            else if(name(1:3).eq.'-FH')then
                i = i + 1
                call mgtarg(i,name)
                read(name,'(bn,f10.0)')fh
            else if(name(1:3).eq.'-PC')then
                i = i + 1
                call mgtarg(i,name)
                read(name,'(bn,f10.0)')pct
            else if(name(1:3).eq.'-PZ')then
                i = i + 1
                call mgtarg(i,fname)
            else if(name(1:2).eq.'-W')then
                i = i + 1
                call mgtarg(i,name)
                read(name,'(bn,f10.0)')water
                if(water.gt.1.0)water = 1.0
                if(water.le.0.0)water= 0.00001
            else if(name(1:2).eq.'-?')then
                call usage(' ')
            else if(name(1:2).eq.'-h')then
                call usage(' ')
            endif
        go to 11
   13   continue
        return
        end

        subroutine usage(str)
        parameter (LER=6, LIN=5, LOT=6)
        character str*(*)
        write(LER,*)'ffilt96:',str
        write(LER,*)'USAGE: ',
     1  'ffilt96  [-?] [-h]',
     2  '-DEMEAN -TAPER -PCT pct ',     
     3  ' -A -R -PZ pole_zero_file -W water_level'
        write(LER,*)
     1  ' -DEMEAN  (default false) remove mean before filter'
        write(LER,*)
     1  ' -TAPER   (default false) apply taper before filter'
        write(LER,*)
     1  ' -PCT pct (default 5.0 percent) taper percentage'
        write(LER,*)
     1  ' -A       (default true) apply filter'
        write(LER,*)
     1  ' -R       (default false) remove filter'
        write(LER,*)
     1  ' -PZ pole_zero_file (none)  SAC response file'
        write(LER,*)
     1  ' -W water_level (0.01)  for removal to not divide by zero'
        write(LER,*)
     1  ' -?           Write this help message'
        write(LER,*)
     1  ' -h           Write this help message'
        stop
        end

        subroutine gtsrsp(fname,ierr)
c-----
c       GeT Sac ReSPonse
c       read the SAC instrument response file, and fill the
c       common block with the instrument response
c
c       fname   C*(*)   - name of pole zero transfer function in 
c                   SAC format
c       ierr    I*4 - < 0 means im proper SAC file
c-----
        character fname*(*)
c-----
c       Pole zero description of complex transfer function
c
c       H(s) = const Product [ s - zero(i) ]/ Product [s - pole(i) ]
c                  
c-----
c       zpoles      Cmp*8   - complex poles
c       zzeros      Cmp*8   - complex zeros
c       npoles      I*4 - number of poles
c       nzeros      I*4 - number of zeros
c       const       R*4 - constant
c       npordr      I*4 - 0 do not use, 1 - >
c                       order in increasing angular frequency
c                       for single or double poles, separately
c       nptyp       I*4 - 1 = single pole, 2=double
c                     2 = part of complex conjucgate pair
c                     0 = do not use, part of complex conjugate pair
c       npardr      I*4 - 0 do not use, 1 - >
c                       order in increasing angular frequency
c                       for all poles, separately
c       pwn     R*4 - natural frequency
c       phn     R*4 - damping for second order
c       pconst      R*4 - modified constant accounting for low pass
c-----
        parameter (NPZ=100)
        common/pzresp/zpoles(NPZ), zzeros(NPZ), npoles, nzeros, 
     1      const, npordr(NPZ), nptyp(NPZ), npardr(NPZ),
     2      pwn(NPZ), phn(NPZ), pconst
        complex zpoles, zzeros
        integer*4 npoles, nzeros, npordr, nptyp, npardr
        real*4 const, pwn, phn, pconst
        
        real*4 tmp(NPZ)
        integer*4 ktmp(NPZ)

        logical ext, inpole, inzero
        logical hascon, haspol, haszer
        parameter (LSTR=80)
        character str*(LSTR)

c-----
c       first initialize the common black
c-----
        do 500 i=1,NPZ
            zpoles(i) = cmplx(0.0,0.0)
            zzeros(i) = cmplx(0.0,0.0)
  500   continue

        npoles = 0
        nzeros = 0

        inquire(file=fname,exist=ext)
        if (.not. ext) return

        open(1,file=fname,status='old',form='formatted',
     1      access='sequential')
        rewind 1
c-----
c       read the input line by line, parsing it
c-----
        inpole = .false.
        inzero = .false.
        hascon = .false.
        haspol = .false.
        haszer = .false.
 1000   continue
        read(1,'(a)',end=9999,err=9999)str
        call detab(str,LSTR)
            l = LSTR
            if(str(1:3).eq.'CON' .or. str(1:3).eq.'con')then
                call gtblnk(str,ipos,1,l)
                read(str(ipos:l),'(e30.3)')const
                hascon = .true.
            else if(str(1:1).eq.'P' .or. str(1:1).eq.'p')then
                call gtblnk(str,ipos,1,l)
                read(str(ipos:l),'(i30)')npoles
                inpole = .true.
                np = 0
                haspol = .true.
            else if(str(1:1).eq.'Z' .or. str(1:1).eq.'z')then
                call gtblnk(str,ipos,1,l)
                read(str(ipos:l),'(i30)')nzeros
                inzero = .true.
                nz = 0
                haszer = .true.
            else if(str(1:1).eq.'*')then
c-----
c               in comment
c-----
                continue
            else
c-----
c               get pole/zero
c-----
                lstrt = 1
                call getnb(str,lstrt,l,l1,l2)
                if(l2.gt.0)then
                    read(str(l1:l2),'(f30.0)')rx
                    lstrt=l2+1
                    call getnb(str,lstrt,l,l1,l2)
                    read(str(l1:l2),'(f30.0)')ry
                    if(inpole)then
                        np = np + 1
                        zpoles(np) = cmplx(rx,ry)
                    else if(inzero)then
                        nz = nz + 1
                        zzeros(nz) = cmplx(rx,ry)
                    endif
                endif
            endif
        go to 1000
 9999   continue
        close (1)
        if(.not. hascon .and. .not. haspol .and. .not. haszer)then
            ierr = -1
            return
        else
            ierr = 1
        endif
c-----
c       now classify the poles according to order of natural
c       frequency, whether simple or double
c       Once done, revise the constant
c-----
        do 2000 i=1,npoles
            x = real(zpoles(i))
            y = aimag(zpoles(i))
            nptyp(i) = 0
            pwn(i) = sqrt(x*x + y*y)
c-----
c       get single pole, not part of a complex conjugate pair
c-----  
            if(y.eq.0.0 .or. abs(y).lt.1.0e-6*abs(x))then
                nptyp(i) = 1
c-----
c       complex conjugate pair
c-----
            else if(y .gt.0.0)then
                nptyp(i) = 2
                phn(i) = -x/pwn(i)
            endif
 2000   continue
        do 2500 i=1,npoles
            npordr(i) = 0
            npardr(i) = 0
 2500   continue
c-----
c       now order the 1st and 2nd order poles by 
c           increasing angular frequency
c-----
        do 3200 korder=1,2
            j = 0
            do 3000 i=1,npoles
                if(nptyp(i).eq.korder)then
                    j = j + 1
                    ktmp(j) = j
                    tmp(j) = pwn(i)
                endif
 3000       continue
            call sort(xtmp,ktmp,j)
c----
c           now indicate the order by the npordr flag
c-----
            j = 0
            do 3100 i=1,npoles
                if(nptyp(i).eq.korder)then
                    j = j + 1
                    npordr(i) = ktmp(j)
                endif
 3100       continue
 3200   continue
c-----
c       now order all  poles by increasing angular frequency
c-----
        j = 0
        do 4000 i=1,npoles
            if(nptyp(i).gt.0)then
                j = j + 1
                ktmp(j) = j
                tmp(j) = pwn(i)
            endif
 4000   continue
        call sort(xtmp,ktmp,j)
c----
c       now indicate the order by the npordr flag
c-----
        j = 0
        do 4100 i=1,npoles
            if(nptyp(i).gt.0)then
                j = j + 1
                npardr(i) = ktmp(j)
            endif
 4100   continue
        return
        end

        subroutine getnb(instr,lstrt,ls,l1,l2)
c-----
c       determine non-black string range
c       
c       instr   Ch  - examined character string
c       lstrt   I   - start position for search
c       ls  I   - length of string
c       l1  I   - index of first non-blank
c       l2  I   - index of last  non-blank
c                   -1 means that no non-black found
c               This program does not guarantee that the
c               resultant string is numberical
c-----
        character instr*(*)
        integer lstrt,ls,l1,l2
        character tab*1
        tab=char(9)
        l1 = lstrt
        l2 = -1
        igotit = 0
        do 1000 i=lstrt,ls
            if(igotit.eq.0)then
            if(instr(i:i).ne.' ' .and. instr(i:i).ne.tab)then
                l1 = i
                l2 = i
                igotit = 1
            endif
            else if(igotit.eq.1)then
            if(instr(i:i).ne.' ' .and. instr(i:i).ne.tab)then
                l2 = i
            else 
                igotit = 2
            endif
            endif
 1000   continue
        return
        end

        subroutine gtblnk(str,ipos,k,l)
        character str*(*)
c-----
c       find the first blank in the stream
c-----
        ipos = k
        do 1000 i=k,l
            if(str(i:i).eq.' ')then
                ipos = i
                return
            endif
 1000   continue
        return
        end

       subroutine sort(x,key,n)
c-----
c     Starting with x(1) ,,, x(n)
c     return   the xarray sorted in increasing order
c     also return the pointers key to the initial array. 
c     For example given x = [ 3, 1, 2 ]
c     the returned values are
c                       x = [ 1, 2, 3 ]        
c                     key = [ 2, 3, 1 ]
c-----
c        Reference: http://en.wikipedia.org/wiki/Bubble_sort
c-----
       integer n
       real x(n)
       integer key(n)

       do i=1,n
           key(i) = i
       enddo
       do i = n, 1, -1
           do j = 1 , i -1
               if(x(j) .gt. x(j+1))then
                   tmp = x(j)
                   x(j) = x(j+1)
                   x(j+1) = tmp
                   ktmp = key(j)
                   key(j) = key(j+1)
                   key(j+1) = ktmp
                endif
           enddo
       enddo
       return
       end

        subroutine detab(str,l)
        character str*(*)
        integer l
        do 1000 i=1,l
            if(str(i:i).eq.char(9))str(i:i) = ' '
 1000   continue    
        return
        end

        subroutine dmean(x,npts)
c-----
c       remove the mean from the time series.
c       Note this will usually remove the DC offset of the
c       trace, but for a one sided pulse this may seriously
c       affect the trace and its spectrum
c-----
c       x   R*4 - time series array
c       npts    I*4 - numbe of samples
c-----
        integer*4 npts
        real*4 x(npts)
        real*8 dsum
c-----
c       remove the mean
c-----
        if(npts.gt.1)then
            dsum  = 0.0d+00
            do 1000 i=1,npts
                dsum = dsum + dble(x(i))
 1000       continue
            dsum = dsum / npts
            xmean = sngl (dsum)
            do 2000 i= 1, npts
                x(i) = x(i) - xmean
 2000       continue
        endif
        return
        end

        subroutine dtape(x,npts,pct)
c-----
c       taper the edges of the trace
c       
c       x   R*4 - input time series
c       npts    I*4 - number of points
c       pct R*4 - percentage of taper
c                 The taper covers pct percent of the total
c                 length of the time series from each end
c           NOTE CHECK WITH SAC CONCEPT
c-----
        integer*4 npts
        real*4 x(npts), pct
        n1 = (pct * npts) / 100
        if(n1.le.1)return
c-----
c       perfrom the taper
c           x(1) = 0.0
c           x(n1) = 1.0
c           ...
c           x(npts - n1 + 1) = 1.0
c           x(npts) = 0.0
c-----
        x(1) = 0.0
        x(npts) = 0.0
        do 1000 i=2,n1
            call costap(i,n1,xx)
            x(i) = xx * x(i)
            x(npts-n1+i) = xx * x(npts-n1+i)
 1000   continue
        return
        end
        
        subroutine costap(i,n,xx)
c-----
c       provide a cosine taper
c-----
c       i   I*4 - pointer 1 <= i <= n
c       n   I*4 - 
c       xx  R*4 - taper function
c                   xx = 0.0 for i = 1
c                   xx = 1.0 for i = n
c                   xx = 1.0 - cos( (i-1)* pi / 2(n-1))
c-----
        pi2 = 3.1415927 / 2.0
c-----
c       safety check
c-----
        if(n.le.1)return
c-----
c       general taper function
c-----
        if(i .le. 1)then
            xx = 0.0
        else if(i.gt.1 .and. i.lt.n)then
            xx = 1.0 - cos( (i-1)*pi2/(n-1))
        else if(i.ge.n)then
            xx = 1.0
        endif
        return
        end



        subroutine filt(doappl,x,npts,dt,water,fl,fh)
c-----
c       doappl  L   - .true. apply filter
c                 .false. remove filter with water level percentage
c       x   R*4 - time series array
c       npts    I*4 - number of points in time series
c       dt  R*4 - sampling interval
c       water   R*4 - a number 0 <= water <=1 used to avoid
c                 divide by zero
c       fl  R*4 - An alternate way to define water level.   
c       fh  R*4 - we use the minimum of the response at fl and fh   
c-----
        logical doappl
        integer*4 npts
        real*4 x(npts), water

        integer*4 NSAMP
        parameter (NSAMP=16384)
        common/zz/z(NSAMP)
        complex z
        common/zr/r(NSAMP)
        complex r
        real*8 za, zmax, zb

        complex zresp
        logical doflfh

c-----
c       get a power of 2
c-----
        n = npts
        call npow2(n)
        n21 = n /2 + 1
c-----
c       zero fill and get FFT
c-----
        do 1000 i=1,n
            if(i.le.npts)then
                z(i) = cmplx(x(i),0.0)
            else
                z(i) = cmplx(0.0,0.0)
            endif
 1000   continue
        call zfour(z,n,-1,dt,df)
c-----
c       check water level
c-----
        doflfh = .false.
        zb = 1.0e+38
        if(fl .gt. 0.0)then
            call zfilt(zresp,fl)
            zb = cabs(zresp)
            doflfh = .true.
        endif
        if(fh .gt. 0.0) then
            call zfilt(zresp,fh)
            if(cabs(zresp) .lt. zb)zb = cabs(zresp)
            doflfh = .true.
        endif
c------
c       get the instrument response
c-----
        if(doappl)then
            do 2000 i=1,n21
                freq = (i-1)*df
                call zfilt(zresp,freq)
                z(i) = z(i) * zresp
                if(i.gt.1)then
                    z(n+2 -i) = conjg(z(i))
                endif
 2000       continue
        else
c-----
c       deconvolve the filter response, first by
c       getting the maximum response value over all 
c           required frequencies
c-----
            zmax = 0.0
            do 3000 i=1,n21
                freq = (i-1)*df
                call zfilt(zresp,freq)
                r(i) = zresp
                if(.not. doflfh)then
                    za = cabs(zresp)
                    if(za .gt. zmax)zmax = za
                endif
 3000       continue
c-----
c       apply the deconvolution using a water level
c       The model is O=output, I=input, H=filter
c
c       O = H I
c       I ~=  conjg(H) O / ( H conjg(H) + sigma**2)
c-----
c-----
            
        write(0,*)'water=',water
        write(0,*)'fl,fh',fl,fh
            do 4000 i=1,n21
                if(doflfh)then
                    denom = r(i)*conjg(r(i)) + 
     1              zb*zb
                else
                    denom = r(i)*conjg(r(i)) + 
     1              water*water*zmax*zmax
                endif
                z(i) = z(i)*conjg(r(i))/denom
                if(i.gt.1)then
                    z(n+2 -i) = conjg(z(i))
                endif
 4000       continue
        endif
        z(n21) = cmplx(real(z(n21)),0.0)
        call zfour(z,n,+1,dt,df)
c-----
c       obtain the inverse FFT
c-----
        do 5000 i=1,npts
            x(i) = real(z(i))
 5000   continue
        return
        end

        subroutine npow2(npts)
c-----
c       Given npts, determine the N=2**m such that N >= npts
c       return the new ntps
c-----
        integer*4 nsamp, npts
        nsamp = npts
        npts = 1
 1000   continue
        npts = 2*npts
        if(npts.lt.nsamp)goto 1000
        return
        end

        subroutine zfour(zarr,nn,isign,dt,df) 
c-----
c     THE input is a complex array
c     which has numbers stored in memory as
c     R1, I1, R2, I2, ..., Rnn, Inn
c     where nn must be a power of 2 R and I are the real and imaginary
c     parts of the complex number
c
c     For isign -1 this is a complex time series
c     For isign +1 this is a complex frequency series with
c        index 1 (in fortran corresponding to f=0
c              2                              f=df
c            nn/2 + 1                         f = 1/2dt = Nyquist
c            nn - 1                           f = -2df
c             nn                              f = -df

c-----
c     the cooley-tookey fast fourier transform in usasi basic fortran
c     transform(j) = sum(datc(i)*w**((i-1)(j-1)), where i and j run
c     from 1 to nn and w = exp(isign*2*pi*sqrt(-1)/nn).  datc is a one-
c     dimensional complex array (i.e., the real and imaginary parts of
c     datc are located immediately adjacent in storage, such as fortran
c     places them) whose length nn is a power of two.  isign
c     is +1 or -1, giving the sign of the transform.  transform values
c     are returned in array datc, replacing the input datc.  the time is
c     proportional to n*log2(n), rather than the usual n**2
c     rms resolution error being bounded by 6*sqrt(i)*log2(nn)*2**(-b),
c     b is the number of bits in the floating point fraction.
c
c     the program computes df from dt, dt from df and checks to see
c     if they are consistent. In addition, the transforms are multiplied
c     by dt or df to make the results dimensionally correct
c
c     This is a slightly modified version of the original Brenner routine
c     The changes were to add physical dimensions to the transform
c     and to make it all complex
c-----
        complex zarr(*) 
        integer nn, isign
        real dt, df

        complex ztemp
c-----
c       ensure that the dt and df are defined and
c       consistent
c-----
        if(dt.eq.0.0) dt = 1./(nn*df) 
        if(df.eq.0.0) df = 1./(nn*dt) 
        if(dt.ne.(nn*df)) df = 1./(nn*dt) 
c-----
c       now begin the transform
c-----
        jj = 1
        do 5 ii=1,nn 
        if(ii .lt. jj) then
              ztemp = zarr(jj)
              zarr(jj) = zarr(ii)
              zarr(ii) = ztemp
        endif
        mm = nn/2
    3   continue
        if(jj.le.mm) then
            go to 55
        else 
              jj = jj - mm
              mm = mm /2
              if(mm.lt.1)then
                  go to 55
              else
                  go to 3
              endif
        endif
   55   continue
        jj = jj + mm
    5   continue
        mmax = 1 
c-----
    6 continue
        if(mmax .lt. nn)then
            go to 7
        else if(mmax .ge. nn)then
            go to 10
        endif
    7   continue
        istep= 2*mmax 
        theta = 6.283185307/(isign*2.0*mmax) 
        sinth=sin(theta/2.) 
        wstpr=-2.*sinth*sinth 
        wstpi=sin(theta) 
        wr=1.0 
        wi=0.0 
        do 9 mm=1,mmax
              do 8 ii=mm,nn,istep
                    jj=ii+mmax
                    ztemp=cmplx(wr,wi)*zarr(jj)
                    zarr(jj) = zarr(ii) - ztemp
                    zarr(ii) = zarr(ii) + ztemp
    8         continue
c-----
c       use trig relations to compute the next sin/cos
c       without actually calling sin() or cos()
c-----
              tempr = wr 
              wr = wr*wstpr-wi*wstpi + wr 
              wi = wi*wstpr+tempr*wstpi + wi 
    9   continue
        mmax = istep 
        go to 6 
c-----
c       transform is done
c-----
   10   continue 
c-----
c     give the arrays the proper physical dimensions
c-----
        if(isign.lt.0)then
c-----
c             time to frequency domain
c-----
              do  ii = 1,nn
                    zarr(ii) = zarr(ii) * dt
              enddo
        else
c-----
c             frequency to time domain
c-----
              do ii = 1,nn
                    zarr(ii) = zarr(ii) * df
              enddo
        endif
        return
        end

        subroutine zfilt(zresp,freq)
c-----
c       evaluate the SAC response at frequency freq
c-----
c       zresp   Cplx    - filter fresponse
c       freq    R*4 - frequency
c-----
        complex zresp
        real*4 freq
        complex*16 s
        complex*16 zzprod
c-----
c       common blocks
c-----
c       Pole zero description of complex transfer function
c
c       H(s) = const Product [ s - zero(i) ]/ Product [s - pole(i) ]
c                  
c-----
c       zpoles      Cmp*8   - complex poles
c       zzeros      Cmp*8   - complex zeros
c       npoles      I*4 - number of poles
c       nzeros      I*4 - number of zeros
c       const       R*4 - constant
c       npordr      I*4 - 0 do not use, 1 - >
c                       order in increasing angular frequency
c                       for single or double poles, separately
c       nptyp       I*4 - 1 = single pole, 2=double
c                     2 = part of complex conjucgate pair
c                     0 = do not use, part of complex conjugate pair
c       npardr      I*4 - 0 do not use, 1 - >
c                       order in increasing angular frequency
c                       for all poles, separately
c       pwn     R*4 - natural frequency
c       phn     R*4 - damping for second order
c       pconst      R*4 - modified constant accounting for low pass
c-----
        parameter (NPZ=100)
        common/pzresp/zpoles(NPZ), zzeros(NPZ), npoles, nzeros, 
     1      const, npordr(NPZ), nptyp(NPZ), npardr(NPZ),
     2      pwn(NPZ), phn(NPZ), pconst
        complex zpoles, zzeros
        integer*4 npoles, nzeros, npordr, nptyp, npardr
        real*4 const, pwn, phn, pconst
        
        s = dcmplx(0.0d+00,dble(6.2831853*freq))
        if(npoles.lt.nzeros)then
            np = nzeros
        else
            np = npoles
        endif
        zzprod = dcmplx(dble(const),0.0d+00)
        do 1000 i=1,np
            if(i.le.nzeros)then
                zzprod = zzprod * (s - zzeros(i))
            endif
            if(i.le.npoles)then
                zzprod = zzprod / (s - zpoles(i))
            endif
 1000   continue

c-----
c       convert to single precision
c-----
        if(cdabs(zzprod).lt.1.0e-37)then
            zresp = cmplx(0.0,0.0)
        else
            zresp = zzprod
        endif
        return
        end
