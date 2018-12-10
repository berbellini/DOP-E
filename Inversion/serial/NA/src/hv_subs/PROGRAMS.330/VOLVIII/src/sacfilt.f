      program sacfilt 
c---------------------------------------------------------------------c 
c                                                                    c 
c     COMPUTER PROGRAMS IN SEISMOLOGY                                c 
c     VOLUME V                                                       c 
c                                                                    c
c     PROGRAM: SACFILT                                               c
c                                                                    c
c     COPYRIGHT 2000                                                 c
c     R. B. Herrmann                                                 c
c     Department of Earth and Atmospheric Sciences                   c
c     Saint Louis University                                         c
c     221 North Grand Boulevard                                      c
c     St. Louis, Missouri 63103                                      c
c     U. S. A.                                                       c
c                                                                    c
c---------------------------------------------------------------------c
c      CHANGES
c      21 MAY 2002 When response is removed, set the MINPER MAXPER
c              on the basis of the water level
c      30 AUG 2002 More work on reading SAC response files in gtsrsp
c      23 MAY 2004 Added FREQLIMITS from GSAC and also permit NPOLES
c              with no poles listed in the manner of SAC
c      09 JAN 2005     LER not defined in subroutine gcmdln 
c          baker@usgs.gov
c       23 JUL 2009 - ensured consistent call to FFT routine with
c               complex arguments. We now call zfour
c-----
c
c    sacfilt -DEMEAN -TAPER -FREQLIMITS f1 f2 f3 f4 -A -R 
c      -AMP amp_file -PHA phase_file -SACIN binary_sac_in
c      -SACOUT binary_sac_out 
c
c    This program applies or removes a complex instrument response
c      given by evalresp(IRIS). The input and output files
c      are SAC binary files.  
c
c      The program searches for keywords POLES, ZEROS and CONSTANT
c      starting in column 1
c      The program ignores blank lines
c      and also lines starting with a *
c      if the number of zeros defined is not matched by
c      entries, then those zeros are assumed to be 0+0i
c-----
c-----
c-----
c      command line arguments
c-----
        logical domean, dotape, doappl
        real*4 f1, f2, f3, f4
        character fname*80, sacin*80, sacout*80
        logical doftaper
c-----
c      variables local to the program
c-----
        integer LER, LIN, LOT
        parameter (LER=0, LIN=5, LOT=6)
        parameter (LN=131072)
        integer*4 npts
        real*4 x(LN)
c-----
c      SAC header files
c-----
        common/sachdr/rhdr,ihdr,chdr
        real*4 rhdr(70)
        integer*4 ihdr(40)
        character*8 chdr(24)
c-----
c      Pole zero description of complex transfer function
c
c      H(s) = const Product [ s - zero(i) ]/ Product [s - pole(i) ]
c                 
c-----
c      zpoles      Cmp*8   - complex poles
c      zzeros      Cmp*8   - complex zeros
c      npoles      I*4 - number of poles
c      nzeros      I*4 - number of zeros
c      const       R*4 - constant
c      npordr      I*4 - 0 do not use, 1 - >
c                      order in increasing angular frequency
c                      for single or double poles, separately
c      nptyp       I*4 - 1 = single pole, 2=double
c                    2 = part of complex conjucgate pair
c                    0 = do not use, part of complex conjugate pair
c      npardr      I*4 - 0 do not use, 1 - >
c                      order in increasing angular frequency
c                      for all poles, separately
c      pwn     R*4 - natural frequency
c      phn     R*4 - damping for second order
c      pconst      R*4 - modified constant accounting for low pass
c-----
        parameter (NPZ=100)
        common/pzresp/zpoles(NPZ), zzeros(NPZ), npoles, nzeros, 
     1      const, npordr(NPZ), nptyp(NPZ), npardr(NPZ),
     2      pwn(NPZ), phn(NPZ), pconst
        complex zpoles, zzeros
        integer*4 npoles, nzeros, npordr, nptyp, npardr
        real*4 const, pwn, phn, pconst
        logical verbos
c-----
c      call machine dependent initialization
c-----
        call mchdep()
c-----
c      parse command line arguments
c-----
        call gcmdln(domean, dotape, doappl, f1, f2, f3, f4,
     1      sacin, sacout, fname, doftaper, verbos)
        write(LER,*)'domean :',domean
        write(LER,*)'dotape :',dotape
        write(LER,*)'doappl :',doappl
        if(doftaper)then
            write(LER,*)'f1     :',f1    
            write(LER,*)'f2     :',f2    
            write(LER,*)'f3     :',f3    
            write(LER,*)'f4     :',f4    
        endif
        write(LER,*)'sacin  :',sacin 
        write(LER,*)'sacout :',sacout
        write(LER,*)'fpolezero:',fname  
c-----
c      get the sac response
c-----
        call gtsrsp(fname,ierr,verbos)
        if(ierr.lt.0)then
            WRITE(LER,*)'Error in reading response file'
            go to 9999
        endif
c-----
c      process
c-----
        IRU = 1
        call brsac (IRU,LN,sacin,x,nerr)
        if(nerr.ge.0)then
            call getnhv('NPTS',npts,nerr)
            call getfhv('DELTA',dt,nerr)
            fnyq = 0.5/dt

            if(domean)call dmean(x,npts)
            if(dotape)call dtape(x,npts,pct)
            call getfhv('USER1', opermin, ierr)
            call getfhv('USER2', opermax, ierr)
            call filt(doappl,x,npts,dt,f1,f2,f3,f4,
     1          doftaper,permin,permax)
            call scmxmn(x,npts,depmax,depmin,depmen,indmax,indmin)
            call getfhv('B     ', btime , ierr)
            call setfhv('TIMMAX', btime + indmax*dt, ierr)
            call setfhv('TIMMIN', btime + indmin*dt, ierr)
            call setfhv('DEPMAX', depmax, ierr)
            call setfhv('DEPMIN', depmin, ierr)
            call setfhv('DEPMEN', depmen, ierr)
c-----
c          permin and permax are set only if the instrument
c          is removed
c-----
            if(opermin .ge. 0.0 .and. opermax .ge.0.0)then
                if(permin.lt.opermin)permin = opermin
                if(permax.gt.opermax)permin = opermax
                call setfhv('USER1', permin, ierr)
                call setfhv('USER2', permax, ierr)
                call setkhv('KUSER1', 'PER_MIN  ', ierr)
                call setkhv('KUSER2', 'PER_MAX  ', ierr)
            endif
c-----
c          output the filtered trace
c-----
            IWU = 2
            call bwsac(IWU,LN,sacout,x)
        else
            WRITE(LER,*)'Trace not processed: nerr=',nerr
            if(nerr.eq.-1)WRITE(LER,*)'Sac file does not exist'
            if(nerr.eq.-2)then
                WRITE(LER,*)'NPTS in file',
     1      ' exceeds dimension',LN
            endif
        endif
            
 9999   continue
        end

        subroutine gcmdln(domean, dotape, doappl, f1, f2, f3, f4,
     1      sacin, sacout, fname, doftaper, verbos)
c-----
c      sacfilt -DEMEAN -TAPER -FREQLIMITS f1 f2 f3 f4 -A -R 
c          -PZ pole_zero_file  -SACIN binary_sac_in
c          -SACOUT binary_sac_out 
c
c      parse command line arguments and return control
c      parameters
c
c      domean  L   - .true. remove mean from time seeries
c      dotape  L   - .true. cosine taper ends
c      doappl  L   - .true. apply filter
c                .false. remove filter
c      doftaper L  If true do frequency taper
c      f1  R*4 Apply zero phase band bass [f1,f2,f3,f4]
c      f2  R*4 where response = 1 for f2 <= f <= f3
c      f3  R*4
c      f4  R*4
c      sacin   Ch*80   - binary sac trace file for input
c      sacout  Ch*80   - binary sac trace file for out
c      fname   Ch*80   - name of polezero file
c
c-----
        integer LER, LIN, LOT
        parameter (LER=0, LIN=5, LOT=6)
        logical domean, dotape, doappl
        real*4 f1, f2, f3, f4
        character sacin*80, sacout*80, fname*80
        logical verbos
        logical doftaper
        
        character*80 name
        integer*4 nmarg
        nmarg = mnmarg()
        domean = .false.
        dotape = .false.
        doappl = .true.
        doftaper = .false.
        f1 = -10.
        f2 = -5.
        f3 = 1.0e+6
        f4 = 1.0e+7
        fname = ' '
        sacin  = ' '
        sacout = ' '
        verbos = .false.
        i = 0
   11   i = i + 1
            if(i.gt.nmarg)goto 13
            call mgtarg(i,name)
            if(name(1:3).eq.'-DE')then
                domean = .true.
            else if(name(1:3).eq.'-TA')then
                dotape = .true.
            else if(name(1:2).eq.'-A'.and.name(1:3).ne.'-AM')then
                doappl = .true.
            else if(name(1:2).eq.'-R')then
                doappl = .false.
            else if(name(1:3).eq.'-FR' .or. name(1:3).eq.'-fr')then
                i = i + 1
                call mgtarg(i,name)
                read(name,'(bn,f10.0)')f1
                i = i + 1
                call mgtarg(i,name)
                read(name,'(bn,f10.0)')f2
                i = i + 1
                call mgtarg(i,name)
                read(name,'(bn,f10.0)')f3
                i = i + 1
                call mgtarg(i,name)
                read(name,'(bn,f10.0)')f4
                if(f4.gt.f3 .and. f3.gt.f2 .and. f2.gt.f1
     1              .and. f1.gt.0.0)then
                    doftaper = .true.
                else
                    doftaper = .false.
                    write(LER,*)'error in FREQLIMITS'
                endif
            else if(name(1:3).eq.'-PZ')then
                i = i + 1
                call mgtarg(i,fname)
            else if(name(1:5).eq.'-SACI')then
                i = i + 1
                call mgtarg(i,sacin)
            else if(name(1:5).eq.'-SACO')then
                i = i + 1
                call mgtarg(i,sacout)
            else if(name(1:2).eq.'-?')then
                call usage(' ')
            else if(name(1:2).eq.'-h')then
                call usage(' ')
            else if(name(1:2).eq.'-V')then
                verbos = .true.
            endif
        go to 11
   13   continue
        return
        end

        subroutine usage(str)
        parameter (LER=0, LIN=5, LOT=6)
        character str*(*)
        write(LER,*)'sacfilt:',str
        write(LER,*)'USAGE: ',
     1       'sacfilt -DEMEAN -TAPER -FREQLIMITS f1 f2 f3 f4 -A -R'
        write(LER,*)
     1  '       -PZ polezero_file -SACIN binary_sac_in'
        write(LER,*)
     1  '       -SACOUT binary_sac_out '
        write(LER,*)
     1  ' -DEMEAN  (default false) remove mean before filter'
        write(LER,*)
     1  ' -TAPER   (default false) apply taper before filter'
        write(LER,*)
     1  ' -FREQLIMITS f1 f2 f3 f4 (default all pass)',
     2  '          Mute for f < f1, f > f4, taper between f1-f2',
     3  '          f3-f4, pass between f2-f3'
        write(LER,*)
     1  ' -A       (default true) apply filter'
        write(LER,*)
     1  ' -R       (default false) remove filter'
        write(LER,*)
     1  ' -PZ pole_zero_file (none)  SAC response file'
        write(LER,*)
     1  ' -SACIN  binary_sac_input  file (none)'
        write(LER,*)
     1  ' -SACOUT binary_sac_output file (none)'
        write(LER,*)
     1  ' -?           Write this help message'
        write(LER,*)
     1  ' -h           Write this help message'
        write(LER,*)
     1  ' SAC header values set'
        write(LER,*)
     1  '  USER1 :  per_min       KUSER1:  PER_MIN'
        write(LER,*)
     1  '  USER2 :  per_max       KUSER2:  PER_MAX'
        stop
        end

        subroutine gtsrsp(fname,ierr,verbos)
c-----
c      GeT Sac ReSPonse
c      read the SAC instrument response file, and fill the
c      common block with the instrument response
c
c      fname   C*(*)   - name of pole zero transfer function in 
c                  SAC format
c      ierr    I*4 - < 0 means im proper SAC file
c-----
        character fname*(*)
        logical verbos
c-----
c      Pole zero description of complex transfer function
c
c      H(s) = const Product [ s - zero(i) ]/ Product [s - pole(i) ]
c                 
c-----
c      zpoles      Cmp*8   - complex poles
c      zzeros      Cmp*8   - complex zeros
c      npoles      I*4 - number of poles
c      nzeros      I*4 - number of zeros
c      const       R*4 - constant
c      npordr      I*4 - 0 do not use, 1 - >
c                      order in increasing angular frequency
c                      for single or double poles, separately
c      nptyp       I*4 - 1 = single pole, 2=double
c                    2 = part of complex conjucgate pair
c                    0 = do not use, part of complex conjugate pair
c      npardr      I*4 - 0 do not use, 1 - >
c                      order in increasing angular frequency
c                      for all poles, separately
c      pwn     R*4 - natural frequency
c      phn     R*4 - damping for second order
c      pconst      R*4 - modified constant accounting for low pass
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
        parameter (LSTR=180)
        character str*(LSTR)
        integer ls

        integer LER, LIN, LOT
        parameter (LER=0, LIN=5, LOT=6)
c-----
c      first initialize the common black
c-----
        do 500 i=1,NPZ
            zpoles(i) = cmplx(0.0,0.0)
            zzeros(i) = cmplx(0.0,0.0)
  500   continue

        npoles = 0
        nzeros = 0

        ierr = -1
        inquire(file=fname,exist=ext)
        if (.not. ext) return

        open(1,file=fname,status='old',form='formatted',
     1      access='sequential')
        rewind 1
c-----
c      read the input line by line, parsing it
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
                ls = lgstr(str(1:l))
                call chtofp(str(ipos:ls), const)
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
c              in comment
c-----
                continue
            else
c-----
c              get pole/zero
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
c      now classify the poles according to order of natural
c      frequency, whether simple or double
c      Once done, revise the constant
c-----
        do 2000 i=1,npoles
            x = real(zpoles(i))
            y = aimag(zpoles(i))
            nptyp(i) = 0
            pwn(i) = sqrt(x*x + y*y)
c-----
c      get single pole, not part of a complex conjugate pair
c-----  
            if(y.eq.0.0 .or. abs(y).lt.1.0e-6*abs(x))then
                nptyp(i) = 1
c-----
c      complex conjugate pair
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
c      now order the 1st and 2nd order poles by increasing 
c      angular frequency
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
c          now indicate the order by the npordr flag
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
c      now order all  poles by increasing angular frequency
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
c      now indicate the order by the npordr flag
c-----
        j = 0
        do 4100 i=1,npoles
            if(nptyp(i).gt.0)then
                j = j + 1
                npardr(i) = ktmp(j)
            endif
 4100   continue
c-----
c      verbose output
c-----
        if(verbos)then
            write(LOT,*)'CONSTANT',const
            write(LOT,*)'ZEROS   ',nzeros
            do 5000 nz=1,nzeros
                write(LOT,*)zzeros(nz)
 5000       continue
            write(LOT,*)'POLES   ',npoles
            do 5100 np=1,npoles
                write(LOT,*)zpoles(np)
 5100       continue
        endif
        return
        end

        subroutine getnb(instr,lstrt,ls,l1,l2)
c-----
c      determine non-black string range
c      
c      instr   Ch  - examined character string
c      lstrt   I   - start position for search
c      ls  I   - length of string
c      l1  I   - index of first non-blank
c      l2  I   - index of last  non-blank
c                  -1 means that no non-black found
c              This program does not guarantee that the
c              resultant string is numberical
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
c      find the first blank in the stream
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

        subroutine detab(str,L)
c-----
c      remove tabs from input string
c-----
        character str*(*)
        integer L
        do 1000 i=1,L
            if(str(i:i).eq.char(9))str(i:i) = ' '
 1000   continue    
        return
        end

        subroutine dmean(x,npts)
c-----
c      remove the mean from the time series.
c      Note this will usually remove the DC offset of the
c      trace, but for a one sided pulse this may seriously
c      affect the trace and its spectrum
c-----
c      x   R*4 - time series array
c      npts    I*4 - numbe of samples
c-----
        integer*4 npts
        real*4 x(npts)
        real*8 dsum
c-----
c      remove the mean
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
c      taper the edges of the trace
c      
c      x   R*4 - input time series
c      npts    I*4 - number of points
c      pct R*4 - percentage of taper
c                The taper covers pct percent of the total
c                length of the time series from each end
c          NOTE CHECK WITH SAC CONCEPT
c-----
        integer*4 npts
        real*4 x(npts), pct
        n1 = (pct * npts) / 100
        if(n1.le.1)return
c-----
c      perfrom the taper
c          x(1) = 0.0
c          x(n1) = 1.0
c          ...
c          x(npts - n1 + 1) = 1.0
c          x(npts) = 0.0
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
c      provide a cosine taper
c-----
c      i   I*4 - pointer 1 <= i <= n
c      n   I*4 - 
c      xx  R*4 - taper function
c                  xx = 0.0 for i = 1
c                  xx = 1.0 for i = n
c                  xx = 1.0 - cos( (i-1)* pi / 2(n-1))
c-----
        pi2 = 3.1415927 / 2.0
c-----
c      safety check
c-----
        if(n.le.1)return
c-----
c      general taper function
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

        subroutine filt(doappl,x,npts,dt,f1,f2,f3,f4,
     1          doftaper,permin,permax)
c-----
c      doappl  L   - .true. apply filter
c                .false. remove filter with water level percentage
c      x   R*4 - time series array
c      npts    I*4 - number of points in time series
c      dt  R*4 - sampling interval
c      water   R*4 - a number 0 <= water <=1 used to avoid
c                divide by zero
c      fl  R*4 - An alternate way to define water level.   
c      fh  R*4 - we use the minimum of the response at fl and fh   
c      permin  R*4 minimum period 
c      permax  R*4 maximum period of passband defined approximately as
c              frequencies where response is > 2 water*zmax or
c              2 zb
c-----
        logical doappl
        integer*4 npts
        real*4 x(npts)
        real dt
        real f1, f2, f3, f4
        logical doftaper
        real permin, permax
        real taper3
        

        integer*4 NSAMP
        parameter (NSAMP=131072)
        common/zz/z(NSAMP)
        complex z
        common/zr/r(NSAMP)
        complex r

        complex zresp

c-----
c      get a power of 2
c-----
        n = npts
        call npow2(n)
        n21 = n /2 + 1
c-----
c      zero fill and get FFT
c-----
        do 1000 i=1,n
            if(i.le.npts)then
                z(i) = cmplx(x(i),0.0)
            else
                z(i) = cmplx(0.0,0.0)
            endif
 1000   continue
        call zfour(z,n,-1,dt,df)
c------
c      get the instrument response
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
c      deconvolve the filter response
c-----
            do 3000 i=1,n21
                if(i.eq.1)then
                    freq = 0.01*df
                else
                    freq = (i-1)*df
                endif
                call zfilt(zresp,freq)
                if(doftaper)then
                    if(freq.lt.f2)then
                        taper= taper3(f1,f2,freq)
                    else if(freq.gt.f3)then
                        taper= taper3(f4,f3,freq)
                    else 
                        taper = 1.0
                    endif
                else
                    taper = 1.0
                endif
                z(i) = z(i)*taper/zresp
                if(i.gt.1)then
                    z(n+2 -i) = conjg(z(i))
                endif
 3000       continue
c-----
c          set frequency limits for the header
c-----
            if(doftaper)then
                    permax = 1.0/f2
                    permin = 1.0/f3
            else
                    permax = n*dt
                    permin = 1.0/(2.0*dt)
            endif
        endif
        z(n21) = cmplx(real(z(n21)),0.0)
        call zfour(z,n,+1,dt,df)
c-----
c      obtain the inverse FFT
c-----
        do 5000 i=1,npts
            x(i) = real(z(i))
 5000   continue
        return
        end

        subroutine npow2(npts)
c-----
c      Given npts, determine the N=2**m such that N >= npts
c      return the new ntps
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
c      evaluate the SAC response at frequency freq
c-----
c      zresp   Cplx    - filter fresponse
c      freq    R*4 - frequency
c-----
        complex zresp
        real*4 freq
        complex*16 s
        complex*16 zzprod
c-----
c      common blocks
c-----
c      Pole zero description of complex transfer function
c
c      H(s) = const Product [ s - zero(i) ]/ Product [s - pole(i) ]
c                 
c-----
c      zpoles      Cmp*8   - complex poles
c      zzeros      Cmp*8   - complex zeros
c      npoles      I*4 - number of poles
c      nzeros      I*4 - number of zeros
c      const       R*4 - constant
c      npordr      I*4 - 0 do not use, 1 - >
c                      order in increasing angular frequency
c                      for single or double poles, separately
c      nptyp       I*4 - 1 = single pole, 2=double
c                    2 = part of complex conjucgate pair
c                    0 = do not use, part of complex conjugate pair
c      npardr      I*4 - 0 do not use, 1 - >
c                      order in increasing angular frequency
c                      for all poles, separately
c      pwn     R*4 - natural frequency
c      phn     R*4 - damping for second order
c      pconst      R*4 - modified constant accounting for low pass
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
c      convert to single precision
c-----
        if(cdabs(zzprod).lt.1.0e-37)then
            zresp = cmplx(0.0,0.0)
        else
            zresp = zzprod
        endif
        return
        end

        subroutine chtofp(str,fout)
c------
c      routine to convert string to floating point
c      The E format is accepted as well as free form
c      input
c
c      This assumes that the string is 20 characters long
c-----
        character*(*) str
        real*4 fout
        integer*4 lgstr
        logical hase
        l = lgstr(str)
c------
c      If the string str contains an E or e, then
c      we parse the result using an E format.
c------
        hase = .false.
        do 100 i=1,l
            if(str(i:i) .eq.'e' .or. str(i:i).eq.'E')then
                hase = .true.
            endif
  100   continue
c-----
c      read floating point number
c-----
        if(hase)then
            read(str,'(bn,e20.0)')fout
        else
            read(str,'(bn,f20.0)') fout
        endif
        return
        end

        real function taper3(xl,xh,x)
        real xl, xh, x
        real p
c-----
c      cubic taper between 0 and 1 for xl <= x <= xh,
c      with = 1 for x = xh and 0 for x = xl,
c-----
        if(xl .eq. xh)then
            taper3 = 1.0
        else
            p = (x - xl)/(xh - xl)
            p = 2.0*p - 1.0
            if(p .le. -1.0)then
                taper3 = 0.0
            else if(p .ge. 1.0)then
                taper3 = 1.0
            else 
                taper3 = 0.5 + 0.75*p*(1.0 - p*p/3.0)
            endif
        endif
        return
        end
