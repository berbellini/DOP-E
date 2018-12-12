c----------------------------------------------------------------------c
c                                                                      c
c      COMPUTER PROGRAMS IN SEISMOLOGY                                 c
c      VOLUME III                                                      c
c                                                                      c
c      PROGRAM:   SLAT2D96                                             c
c                                                                      c
c      COPYRIGHT 1996                                                  c
c      R. B. Herrmann                                                  c
c      Department of Earth and Atmospheric Sciences                    c
c      Saint Louis University                                          c
c      221 North Grand Boulevard                                       c
c      St. Louis, Missouri 63103                                       c
c      U. S. A.                                                        c
c                                                                      c
c----------------------------------------------------------------------c
c     Read in control file to make pseudo- s(rl)egn96(III) output      c
c     file for use by spulse96(III) to make seismograms for laterally  c
c     varying earth model.                                             c
c                                                                      c
c     NOTE: spulse96(III) can only be run for ONE distance using this  c
c     new eigenfunction file                                           c
c                                                                      c
c     26 MAY 2007 - careful reworking to make the slat2d work
c     26 AUG 2011 - note that the model output for spulse96 will not 
c           give correct travel time since spulse96 only permits a 1-D
c           model                                                                      c
c----------------------------------------------------------------------c
        parameter (LIN=5, LOT=6, LER=0)
        parameter (NL=1000)
        dimension d(NL),a(NL),b(NL),rho(NL),qa(NL),qb(NL)
        integer*4 mmax,nper,ifunc,kmode

        dimension dd(NL),aa(NL),bb(NL),rhorho(NL),qaqa(NL),qbqb(NL)
        integer*4 nnper,iifunc,kkmode

        character*120 names
        character cmdfil*120

        character mname*120
        integer ipar(10)
        real*4 fpar(10)

        integer ipardb(10)
        real*4 fpardb(10)
c-----
c-----
c       machine dependent initialization
c-----
        call mchdep()
c-----
c       get command line arguments
c-----
        call gcmdln(dist,cmdfil)
c-----
c       open command file which tells width of each region and
c       the data file name of each.
c
c       The command file is copied into an internal work array
c-----
        open(1,access='sequential',form='unformatted',status=
     1      'scratch')
        rewind 1
c-----
c       set up command file
c-----
        call getcmd(distmx,cmdfil,nmod,dist,nwhich)
c-----
c       set up double buffering temporary storage file
c-----
        open(3,access='sequential',form='unformatted',
     1      status='scratch')
        rewind 3
        open(4,access='sequential',form='unformatted',
     1      status='scratch')
        rewind 3
        rewind 4
        rewind 1
c-----
c       start the work, and we will double buffer with a temporary array
c-----
        iout = 3
        iin  = 4
c-----
c       curdst is the distance from the source to the left of
c       the current section
c-----
        curdst = 0.0
        kfunc = 0
        do  2000 nsec=1,nwhich
            rewind iin
            rewind iout
            read(1,end=2000,err=2000)x,names
            open(2,file=names,access='sequential',status='old',
     1          form='unformatted')
            if(nsec.eq.1)then
                call getmdl(2,mmax,d,a,b,rho,qa,qb,nper,
     1              dphs,dphr,
     2              mname,ipar,fpar)
                call putmdl(iout,mmax,d,a,b,rho,qa,qb,
     1              nper,dphs,dphr,
     2              mname,ipar,fpar)
            else if(nsec.gt.1)then
                call getmdl(2,mmax,d,a,b,rho,qa,qb,
     1              nnper,dphs,dphr,
     2              mname,ipar,fpar)
                call getmdl(iin,mmax,dd,aa,bb,rhorho,
     1              qaqa,qbqb,nnper,ddphs,ddphr,
     2              mname,ipardb,fpardb)
                call putmdl(iout,mmax,dd,aa,bb,rhorho,
     1              qaqa,qbqb,nnper,ddphs,ddphr,
     2              mname,ipardb,fpardb)
            endif
            if(nsec.eq.1)then
c-----
c                  check on whether source is in fluid
c-----
                   ipar2sv = ipar(2)
            else if(nsec.eq.nwhich)then
c-----
c                  check on whether receiver is in a fluid
c-----
                   ipar3sv = ipar(3)
            endif
 2100       continue
            if(nsec.eq.1)then
                call gethed(2,ifunc,kmode,t0,ierr)
                kkmode = kmode
                tt0 = t0
                iifunc = ifunc
                if(kfunc.eq.0)kfunc = ifunc
            else
                call gethed(2,iifunc,kkmode,tt0,ierr)
                call gethed(iin,ifunc,kmode,t0,ierr)
            endif
            if(kmode.gt.kkmode)then
                lmode = kkmode
                mmode = kmode
            else
                lmode = kmode
                mmode = kkmode
            endif
c-----
c           kmode is the minimum number of common modes in both sets
c           mmode is the maximum number of modes in either set
c           mmode >= lmode
c           (this means later that it may be necessary to skip some)
c-----
            if(ifunc.gt.0.and.iifunc.ne.kfunc)then
                call usage('inconsistent type for file: '// names)
            endif
            call puthed(iout,ifunc,lmode,t0)
            lorr = ifunc
            if(ifunc.le.0) go to 2200
            if(kmode.le.0) go to 2100
            do 2300 i=1,mmode
                sumkr = 0.0
                sumgr = 0.0
                sumgv = 0.0
                if(nsec.eq.1)then
                    call getegn(2,lorr,1,wvno,u,gamma,
     1              sur,sdur,suz,sduz,sare,wvnsrc,sur0,
     2              rur,rdur,ruz,rduz,rare,wvnrec,rur0,
     3              sumkr,sumgr,sumgv,ierr)
                    wvnor = wvno
                    gvr = u
                    gammar = gamma
                else
c-----
c       here unit 2 has the information on the receiver, and the
c       current iin has the information on the source and also
c       the summation of phase, group and attenuation delays
c-----
                    if(i.le.kkmode)then
                        call getegn(2,lorr,1,wvnor,gvr,gammar,
     1                      x1,x2,x3,x4,x5,x6,x7,
     2                      rur,rdur,ruz,rduz,rare,wvnrec,rur0,
     3                       x10,x11,x12,ierr)
                    endif
                    if(i.le.kmode)then
                        call getegn(iin,lorr,0,wvno,u,gamma,
     1                      sur,sdur,suz,sduz,sare,wvnsrc,sur0,
     2                      x1,x2,x3,x4,x5,x6,x7,
     3                      sumkr,sumgr,sumgv,ierr)
                    endif
                endif
                if(i.le.lmode)then
                    if(nsec.eq.nwhich)then
                        sumkr = sumkr + (dist-curdst)*wvnor
                        sumgr = sumgr + (dist-curdst)*gammar
                        sumgv = sumgv + (dist-curdst)/gvr
                    else
                        sumkr = sumkr + x*wvnor
                        sumgr = sumgr + x*gammar
                        sumgv = sumgv + x/gvr
                    endif
                    call putegn(iout,lorr,0,wvno,u,gamma,
     1              sur,sdur,suz,sduz,sare,wvnsrc,sur0,
     2              rur,rdur,ruz,rduz,rare,wvnrec,rur0,
     3              sumkr,sumgr,sumgv)
                endif

 2300       continue
            goto 2100
 2200   continue
c-----
c       do double buffering
c-----
            itmp = iin
            iin = iout
            iout = itmp
            curdst = curdst + x
            close (2)
 2000   continue
        close (1)
c-----
c       we must be careful here
c       The velocity model must be that of the source medium
c       However, we must output the correctly the conditions 
c       whether the source is in fluid and the receiver is in fluid
c
c       we also override ipar(3) concerning whether the receiver is
c       in the fluid
c-----
c-----
c       set up outfil and generate output by perturbing source
c       terms and computing true phase delay and attenuation average
c-----
        if(kfunc.eq.1)then
            open(2,file='slatl96.egn',form='unformatted',
     1          access='sequential',status='unknown')
        else if(kfunc.eq.2)then
            open(2,file='slatr96.egn',form='unformatted',
     1          access='sequential',status='unknown')
        endif
        rewind 2
            rewind iin
            call getmdl(iin ,mmax,d,a,b,rho,qa,qb,nper,dphs,dphr,
     2              mname,ipar,fpar)
c-----
c           override the fluid definition of source or receiver
c           so that the source is the first model, and the receiver 
c           is the nwhich model
c-----
            ipar(2) = ipar2sv
            ipar(3) = ipar3sv

            call putmdl(2,mmax,d,a,b,rho,qa,qb,nper,dphs,dphr,
     2              mname,ipar,fpar)
 3100       continue
                call gethed(iin,ifunc,kmode,t0,ierr)
c-----
c       proper termination required
c-----
                if(ifunc.lt.0)then
                    call puthed(2,ifunc,kmode,t0)
                else
                    call puthed(2,ifunc+2,kmode,t0)
                endif
            lorr = ifunc
            if(ifunc.le.0) go to 3200
            if(kmode.le.0) go to 3100
            do 3300 i=1,kmode
                call getegn(iin,lorr,0,wvno,u,gamma,
     1              sur,sdur,suz,sduz,sare,wvnsrc,sur0,
     2              rur,rdur,ruz,rduz,rare,wvnrec,rur0,
     3              sumkr,sumgr,sumgv,ierr)
                avg = sumgr/dist
                avk = sumkr/dist
                avv = dist/sumgv
                wvno = avk
                gamma = avg
                u = avv
                call putegn(2,lorr+2,1,wvno,u,gamma,
     1              sur,sdur,suz,sduz,sare,wvnsrc,sur0,
     2              rur,rdur,ruz,rduz,rare,wvnrec,rur0,
     3              sumkr,sumgr,sumgv)

 3300       continue
        goto 3100
 3200   continue
            close (2)
        end

        subroutine getcmd(distmx,cmdfil,nmod,dist,nwhich)
c-----
c       read ascii command file, and save as binary for
c       rapid use. Also use the opportunity to get the
c       total width of the structure
c-----
c       distmx  - R*4   - width of the entire section
c       cmdfil  - C*(*) - name of command file
c       nmod    - I*4   - number of sections in lateral model
c       dist    - R*4   - desired distance for the receiver
c       nwhich  - I*4   - which section that receiver is in
c-----
c       cmdfil consists a one or more lines of entries consisting of
c       width_of_region 'file_name'
c       The quotes are required for list directed input of strings
c-----
        real*4 distmx

        character cmdfil*(*)
        character names*120
        integer*4 nmod
        open(2,file=cmdfil,access='sequential',form='formatted',
     1      status='old')
        rewind 2
c-----
c       distmx is the maximum width of the structure, starting at 0
c       For dist > distmx, the last structure continues infinitely
c-----
        distmx = 0.0
        nmod = 0
 1000   continue
            read(2,*,err=1001,end=1001)x,names
            nmod = nmod + 1
            if(distmx .lt. dist)then
                nwhich = nmod
            endif
            distmx = distmx + x
            write(1)x,names
            goto 1000
 1001   continue
        close (2)
        return
        end

        subroutine gcmdln(dist,cmdfil)
c-----
c     parse command line arguments and return control
c     parameters
c
c     requires subroutine mgtarg(i,name) to return
c           the i'th argument in the string name
c
c     and the function mnmarg() to return the number
c           of arguments excluding program name
c           The first argument is i = 1
c-----
c       dist    - R*4   desired distance
c       cmdfil  - C*(*) name of command file
c-----
      character*120 name
      character cmdfil*(*)
      real*4 dist
      integer*4 mnmarg
      nmarg = mnmarg()
      dist = -1.0
      cmdfil = ' '
      i = 0
   11 i = i + 1
      if(i.gt.nmarg)goto 13
            call mgtarg(i,name)
            if(name(1:2).eq.'-D')then
                  i=i+1
                  call mgtarg(i,name)
                  read(name,'(f20.0)')dist
            else if(name(1:2).eq.'-C')then
                  i=i+1
                  call mgtarg(i,cmdfil)
            else if(name(1:2).eq.'-?')then
            call usage(' ')
            else if(name(1:2).eq.'-h')then
            call usage(' ')
            endif
            go to 11
   13 continue
      if(dist .le. 0.0 .or. cmdfil .eq. ' ')call usage(' ')
      return
      end

        subroutine usage(str)
        character str*(*)
        integer LER
        parameter (LER=0)
        write(LER,*)'Usage: ',str
        write(LER,*)'USAGE: slat2d96 -D dist -C cmdfil -? -h'
        write(LER,*)
     1  ' -D dist   (no default) desired distance'
        write(LER,*)
     1  ' -C cmdfil (no default) file with distance-eigenfunction',
     2      'information'
        write(LER,*)
     1  ' -?         On line help'
        write(LER,*)
     1  ' -h         On line help'

        stop
        end

