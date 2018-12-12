        program wvfmch96
c---------------------------------------------------------------------c
c                                                                   c
c    COMPUTER PROGRAMS IN SEISMOLOGY                                c
c    VOLUME VIII                                                    c
c                                                                   c
c    PROGRAM: WVFMCH96                                              c
c                                                                   c
c    COPYRIGHT 1994 Charles J. Ammon                                c
c                                                                   c
c    Department of Geosciences                                      c
c    Penn State University                                          c
c    Deike Bldg                                                     c
c    University Park, PA 16802                                      c
c    U. S. A.                                                       c
c                                                                   c
c---------------------------------------------------------------------c
c     FIXES
c     21 OCT 2002 - corrections to SAC header for prediction
c     07 JAN 2003 - replaced Numerical Recipes routines
c     16 SEP 2004 - the predicted have the correct 
c                   year day hr min sec msec
c     16 FEB 2006 - output the observed file as an .obs
c     13 MAR 2006 - restructure - build in dynamic 
c                   weighting and time shift
c         following ideas of Harley Benz
c     04 OCT 2006 - In an effort to both comply with and promulgate 
c          international standards, the NEIC Policy and Procedures 
c          Working Group has made the decision to adopt the 
c          international standard for converting scalar moment 
c          to moment magnitude proposed by the IASPEI Commission on
c          Seismological Observation and Interpretation (CoSOI) 
c          Working Group on Magnitude --
c     
c         old     Mw = 2/3 log Mo - 10.7    Mo in dyne-cm
c         new     Mw = 2/3(log Mo - 16.1)   Mo in dyne-cm
c                 Mw = 2/3(log Mo -  9.1)   Mo in N-m
c       15 MAR 2007 - put in error checking for time series with
c         less than two points - so that not used
c       19 JUL 2007 - Whoa - removed the cross-correlation time 
c         shift stuff for the synthetic put preserved the absolute time
c       11 AUG 2007 - cleaned up code so that it works correctly when the
c                     moment tensor is specified. In addition added a new
c                     source -CRACK which requires -DIP dip -STK dipdir 
c                     MW/M0 to define a expanding/closing crack
c       09 NOV 2007 - modified the code to remove any limit on the maximum
c                     number of waveforms. This is done by processing each
c                     input line separately, instead of reading all into
c                     memory.  The code change required the introduction of
c                     a logical flag to indicate termination of the
c                     input sequence.  This means that the front-end
c                     is now slightly different than that of wvfgrd96 
c	05 JAN 2009 - ensure that KNETWK is set in the .obs and .pre files
c-----
c
c    deviatoric moment tensor inversion
c
c    Time function is fixed, so problem is linear
c
c    Author:  Chuck Ammon, UC Santa Cruz, Saint Louis University
c             parts by George Randall, U South Carolina
c
c    Version: 2.1 September, 1994
c
c    Notation from Langston (1981)
c      x is north, y is east, z is down
c
c-----
c     maxpts  Maximum number of points in time series
c     nprm    Dimension of moment tensor matrix
c     mxwvs   maximum number of observed traces considered
c     nmax    dimension of array
c-----
        integer maxpts, nprm, nmax, mxwvs
        parameter(maxpts=2048, nprm=6, mxwvs=1000, nmax =mxwvs*maxpts)
c-----
c     common blocks for trace and Green functions for station/component
c-----
        common /gftns/ vss(maxpts),vds(maxpts),vdd(maxpts),vex(maxpts),
     1                response(maxpts)
        real vss, vds, vdd, vex, response
        common /observations/ obs(maxpts)
        real obs
c-----
c     common block for matrix inversion and solution
c-----
        common/mat/a,asave,datav,datasv,covmatrix,mtensor
        real a(nmax,nprm), asave(nmax,nprm)
        real datav(nmax), datasv(nmax)
        real covmatrix(nprm,nprm)
        real mtensor(6)
c-----
c     common blocks for traces and trace attributes
c-----
        common/rtrace/az,wt,user1,user2,
     1      disto, depthg, dt, distg,
     2      tpg, tsg, bg, og,
     3      tpo, tso, bo, oo,
     4      sum_sq_amp, fr_errors, obswt, misfit
        real az(mxwvs), wt(mxwvs), user1(mxwvs), user2(mxwvs)
        real disto(mxwvs), depthg(mxwvs), dt(mxwvs), distg(mxwvs)
        real tpg(mxwvs), tsg(mxwvs), bg(mxwvs), og(mxwvs)
        real tpo(mxwvs), tso(mxwvs), bo(mxwvs), oo(mxwvs)
        real sum_sq_amp(mxwvs), fr_errors(mxwvs), misfit(mxwvs)
        real obswt(mxwvs)
        common/itrace/wvtype, mshft, npts, douse, itotal_pts,
     1      nzyear,nzjday, nzhour, nzmin, nzsec, nzmsec
        integer wvtype(mxwvs), mshft(mxwvs), npts(mxwvs), douse(mxwvs)
        integer itotal_pts
        integer nzyear(mxwvs), nzjday(mxwvs), nzhour(mxwvs), 
     1      nzmin(mxwvs), nzsec(mxwvs), nzmsec(mxwvs)
        common/ctrace/obs_name, gftn_name, kstnm, kcmpnm, knetwk
        character*256 obs_name(mxwvs), gftn_name(mxwvs)
        character*8 kstnm(mxwvs), kcmpnm(mxwvs), knetwk(mxwvs)

        common/event/evla, evlo
        real evla, evlo
c-----
c     common for best double couple and for solution vector
c-----
        real pred(nmax)
        real stk0,dip0,rak0,stk1,dip1,rak1,xmw,pclvd

        real btime,  m0
        real total_sum_sq, total_sum_sq_wt

        real dip, stk, rake, xmom
        real xmt(3,3)
        integer isds
c-----
c     variables for this routine
c-----

        integer ibegin
        integer stdin,stdout,outunit,inunit
        integer ipass
        integer nwaves
        logical ldoexit

        character*256 infile
        character*256 outfile

        stdin = 5
        stdout = 6
        outunit = 4
        inunit = 11


c-----
c     parse command line arguments
c-----
        call gcmdln(dip,rake,strike,xmom,isds,xmt)

c-----
c     initialize everything
c-----
        ldoexit = .false.
c-----
c       read the control from the standard input
c       and open files reading waveforms and output
c-----
        call getcontrol(stdin,stdout,inunit,outunit,
     1       infile,outfile)

c-----
c       now loop over the input
c       we only require ONE trace in memory. Thus we never exceed the
c       mxwvs variable
c-----
 1000   continue
        if(ldoexit )then
                go to 9999
        endif
              call initialize(inunit, nwaves,ldoexit,
     1            m0,douse,obswt,
     2            outfile, wvtype, obs_name, gftn_name, wt, mshft)
c-----
c                Output some event and inversion 
c                     information (using the first data file)
c-----
              call getevent_data(obs_name(1),outunit,'wvfmch96')
c-----

c-----
c   END of SETUP
c-----
c     begin the three passes
c     pass 1: do the inversion, determine time shift, set weights
c     pass 2: do the inversion with shifted time series
c     pass 3: final run output obs and pre traces
c-----
              do ipass=2,2

              ibegin = 1
              total_sum_sq = 0
              total_sum_sq_wt = 0
              itotal_pts = 0
c-----
c    Now read in the data, GFs, and set up the matrix equations
c-----
              call getobsgf(nwaves,ibegin,
     1             btime,outunit,stdout,
     2             total_sum_sq,total_sum_sq_wt,ipass)

c-----
c            MAKE THE SYNTHETIC
c-----
             call domtfwd(ibegin,itotal_pts,outunit,m0,dip,rake,strike,
     1            xmom,isds,xmt)
c-----
c     output waveforms and compute fit to individual traces
c-----
                 call dosumup(nwaves,outunit,m0,
     1           total_sum_sq, total_sum_sq_wt, 
     4           ipass)
        enddo
        go to 1000
 9999   continue
c-----
c     be sure to close all open files
c-----
        close( inunit)
        close(outunit)

        end
c-----
c    END of MAIN PROGRAM
c-----
        subroutine getcontrol(stdin,stdout,inunit,outunit,
     1       infile,outfile)
        character infile*(*), outfile*(*)
        integer stdin, stdout,inunit,outunit
        logical lexists
c-----
c       read the standard input to get the file names
c-----
c-----
        write(stdout,*)'What is the input file name?'
        read(stdin,*) infile
c-----
        write(stdout,*)'What is the output file name?'
        read(stdin,*) outfile
c-----
        INQUIRE(file=infile, EXIST = lexists)
        if(.not. lexists) then
           write(stdout,*)'Input file not found - quitting.'
           stop
        endif
c-----
c       open the input file
c-----
        open(unit = inunit, file = infile)
c-----
c                Output some event and inversion 
c                     information (using the first data file)
c-----
              open(unit = outunit, file = outfile)
        return
        end

        subroutine gcmdln(dip,rake,strike,xmom,isds,xmt)
        parameter (LER=0, LIN=5, LOT=6)
        real dip, rake, strike, xmom, xmt(3,3)
        integer isds
c-----
c     parse command line arguments
c
c     requires subroutine mgtarg() and function mnmarg()
c
c-----
c     dip R*4 - dip of fault
c     rake    R*4 - rake angle of motion on fault
c     strike  R*4 - strike of fault
c     isds    I*4 - indicator of couple source description
c               -1 none given
c                0 strike, dip, rake
c                1 moment tensor
c                2 explosion
c     xmom    R*4 - seismic moment of strike, dip, rake form
c               1.0 == 1.0 dyne-cm for km,gm,km/s system
c     isds    I*4 - indicator of couple source description
c                   -1 none given
c                    0 strike, dip, rake
c                    1 moment tensor
c                    2 explosion
c                      crack - requires dip and strike and Mw/Mo
c                     actually this is the moment tensor solution
c-----
        character*25 name
        integer*4 mnmarg
        logical docrack

        dip = 0.0
        rake = 0.0
        strike = 0.0
        isds = -1
        xmom=1.0
        docrack = .false.
        do 120 i=1,3
            do 121 j=1,3
                xmt(j,i) = 0.0
  121       continue
  120   continue

        nmarg=mnmarg()
        if(nmarg.le.0)then
            call usage()
        endif
        i = 0
   11   i = i + 1
        if(i.gt.nmarg)go to 13
            call mgtarg(i,name)
            if(name(1:2).eq.'-D')then
                i=i+1
                call mgtarg(i,name)
                call chtofp(name,dip)
                isds = 0
            else if(name(1:2).eq.'-S')then
                i=i+1
                call mgtarg(i,name)
                call chtofp(name,strike)
                isds = 0
            else if(name(1:2).eq.'-R'.and.name(1:3).ne.'-RO')then
                i=i+1
                call mgtarg(i,name)
                call chtofp(name,rake)
                isds = 0
            else if(name(1:3).eq.'-M0' .or. name(1:3).eq.'-MO')then
                i=i+1
                call mgtarg(i,name)
                call chtofp(name,xmom)
            else if(name(1:3).eq.'-MW' .or. name(1:3).eq.'-Mw')then
                i=i+1
                call mgtarg(i,name)
                call chtofp(name,xmom)
                xmom = 10.**(1.5*xmom + 16.10)
            else if(name(1:3).eq.'-xx' .or. name(1:3).eq.'-XX')then
                i=i+1
                call mgtarg(i,name)
                call chtofp(name,xmt(1,1))
                isds = 1
            else if(name(1:3).eq.'-yy' .or. name(1:3).eq.'-YY')then
                i=i+1
                call mgtarg(i,name)
                call chtofp(name,xmt(2,2))
                isds = 1
            else if(name(1:3).eq.'-zz' .or. name(1:3).eq.'-ZZ')then
                i=i+1
                call mgtarg(i,name)
                call chtofp(name,xmt(3,3))
                isds = 1
            else if(name(1:3).eq.'-xy' .or. name(1:3).eq.'-XY')then
                i=i+1
                call mgtarg(i,name)
                call chtofp(name,xmt(1,2))
                xmt(2,1) = xmt(1,2)
                isds = 1
            else if(name(1:3).eq.'-xz' .or. name(1:3).eq.'-XZ')then
                i=i+1
                call mgtarg(i,name)
                call chtofp(name,xmt(1,3))
                xmt(3,1) = xmt(1,3)
                isds = 1
            else if(name(1:3).eq.'-yz' .or. name(1:3).eq.'-YZ')then
                i=i+1
                call mgtarg(i,name)
                call chtofp(name,xmt(2,3))
                xmt(3,2) = xmt(2,3)
                isds = 1
            else if(name(1:2).eq.'-?')then
                call usage()
            else if(name(1:2).eq.'-h')then
                call usage()
            endif
        goto 11
   13   continue
c-----
c       special handling for the crack source. convert it to a moment
c       tensor
c       Nakano, M., and H. Kumagi (2005). Waveform inversion of volcano-seismic 
c            signals assuming possible source geometries,
c            Geophys. Res. Letters 32, L12302, doi:10.1029/2005GL022666.
c       the conversion to moment tensor must be deferred until later in the code
c       since the lambda and mu are required at the source depth
c-----
        if(docrack)then
              isds = 1
              call makecrack(xmt,strike,dip,rake,xmom)
        endif
c-----
c       convert everything to the units required for the synthetics.
c       The synthetics are generated using KM for distance and
c       layer thickness, KM/SEC for velocity, and GM/CC for density.
c       The conversions below cause accoount for the disjoint between
c       these mixed units and the CM CM/SEC and GM/CC required
c-----

        if(isds.eq.1)then
                xmt(1,1) = xmt(1,1) / 1.0E+20
                xmt(1,2) = xmt(1,2) / 1.0E+20
                xmt(1,3) = xmt(1,3) / 1.0E+20
                xmt(2,3) = xmt(2,3) / 1.0E+20
                xmt(3,3) = xmt(3,3) / 1.0E+20
                xmt(2,2) = xmt(2,2) / 1.0E+20
                xmt(3,2) = xmt(2,3)
                xmt(3,1) = xmt(1,3)
                xmt(2,1) = xmt(1,2)
        else if(isds.ne.1)then
                 xmom   = abs(xmom)   / 1.0E+20
        endif
c-----
c     convert everything to the units required for the synthetics.
c     The synthetics are generated using KM for distance and
c     layer thickness, KM/SEC for velocity, and GM/CC for density.
c-----
        return
        end
        
        subroutine usage()
        parameter (LER=0, LIN=5, LOT=6)
        write(LER,*)'wvfmch96 -D Dip -S Stk -R Rake',
     1      ' -M0 Mom -MW Mw -CRACK'
        write(LER,*)' -XX Mxx -YY Myy -ZZ Mzz -XY -Mxy -XZ Mxz -YZ Myz'
        write(LER,*)
     1  ' '
        write(LER,*)
     1  ' -D Dip               dip of fault plane'
        write(LER,*)
     1  ' -S Strike            strike of fault plane'
        write(LER,*)
     1  ' -R Rake              slip angle on fault plane'
        write(LER,*)
     1  ' -M0 Moment (def=1.0) Seismic moment in units of dyne-cm'
        write(LER,*)
     1  ' -MW mw               Moment Magnitude  '
        write(LER,*)
     1  ' -CRACK       (default no ) use crack model: Dip here is dip ',
     2  '  of crack, Strike is dip direction of crack',
     3  '  rake and Mw or M0, where M0 > 0 and'
        write(LER,*)
     1  '               Rake > 0 for expanding crack',
     2  '                    < 0 for closing crack. '
        write(LER,*)
     1  '               M0= sgn(rake) mu DELTA Volume'
        write(LER,*)
     1  '               Poisson ratio = 0.25'
        write(LER,*)
     1  ' -XX Mxx -YY Myy -ZZ Mzz  Moment tensor elements in units of'
        write(LER,*)
     1  ' -XY Mxy -XZ Mxz -YZ Myz    dyne-cm'
        write(LER,*)
     1  ' -?                   This online help'
        write(LER,*)
     1  ' -h                   This online help'
        stop
        end

        subroutine chtofp(str,fout)
c------
c     routine to convert string to floating point
c     The E format is accepted as well as free form
c     input
c
c     This assumes that the string is 20 characters long
c-----
        character*(*) str
        real*4 fout
        integer*4 lgstr
        logical hase
        l = lgstr(str)
c------
c     If the string str contains an E or e, then
c     we parse the result using an E format.
c------
        hase = .false.
        do 100 i=1,l
            if(str(i:i) .eq.'e' .or. str(i:i).eq.'E')then
                hase = .true.
            endif
  100   continue
c-----
c     read floating point number
c-----
        if(hase)then
            read(str,'(bn,e20.0)')fout
        else
            read(str,'(bn,f20.0)') fout
        endif
        return
        end

        subroutine initialize(inunit, nwaves, ldoexit,
     1      m0,douse,obswt,
     2      outfile, wvtype, obs_name, gftn_name, wt, mshft)
        implicit none
        integer maxpts, nprm, nmax, mxwvs
        parameter(maxpts=2048, nprm=6, mxwvs=1000, nmax =mxwvs*maxpts)
        integer inunit, nwaves
        integer wvtype(mxwvs), mshft(mxwvs), douse(mxwvs)
        integer mpts
        real m0, wt(mxwvs), obswt(mxwvs)
        real depmax, depmin
        character*256 obs_name(mxwvs), gftn_name(mxwvs), outfile, infile
        integer i, nerr
        logical ldoexit

        inunit = 11
c-----
c     m0 is the moment scaling factor for the Green s functions
c     For Computer programs in seismology with earth model in
c     km/sec and gm/cc, m0 == 1.0+20
c-----
        m0 = 1.0e+20

        i = 1
        nwaves = 0
    5   continue
            read(inunit,*,end=6)wvtype(i),obs_name(i),gftn_name(i),
     1          wt(i)
            nwaves = i
            mshft(i) = 0
            douse(i) = 1
            obswt(i) = 1
c-----
c     protection against exceeding array dimension
c-----
            if(nwaves.eq.mxwvs)go to 6
c-----
c       ensure that there are at least two points per time series
c       this is to avoid the problem with 0 points
c-----
            call brsach(1,obs_name(i),nerr)
            call getnhv('NPTS    ',mpts,nerr)
            call getfhv('DEPMAX  ',depmax,nerr)
            call getfhv('DEPMIN  ',depmin,nerr)
            if(mpts.le.4)then
                  WRITE(0,*)'Rejecting: ',i,' MPTS=',mpts,obs_name(i)
                  go to 5
            endif
            return
    6   continue
        ldoexit = .true.
        return
        end

        subroutine getobsgf(nwaves,ibegin,
     1                  btime,ounit,stdout,
     2          total_sum_sq,total_sum_sq_wt,ipass)
        implicit none
        integer nwaves
        integer maxpts, nprm, nmax, mxwvs
        parameter(maxpts=2048, nprm=6, mxwvs=1000, nmax =mxwvs*maxpts)
c-----
c     common blocks for trace and Green functions for station/component
c-----
        common /gftns/ vss(maxpts),vds(maxpts),vdd(maxpts),vex(maxpts),
     1                response(maxpts)
        real vss, vds, vdd, vex, response
        common /observations/ obs(maxpts)
        real obs
c-----
c     common blocks for traces and trace attributes
c-----
        common/rtrace/az,wt,user1,user2,
     1      disto, depthg, dt, distg,
     2      tpg, tsg, bg, og,
     3      tpo, tso, bo, oo,
     4      sum_sq_amp, fr_errors, obswt, misfit
        real az(mxwvs), wt(mxwvs), user1(mxwvs), user2(mxwvs)
        real disto(mxwvs), depthg(mxwvs), dt(mxwvs), distg(mxwvs)
        real tpg(mxwvs), tsg(mxwvs), bg(mxwvs), og(mxwvs)
        real tpo(mxwvs), tso(mxwvs), bo(mxwvs), oo(mxwvs)
        real sum_sq_amp(mxwvs), fr_errors(mxwvs), misfit(mxwvs)
        real obswt(mxwvs)
        common/itrace/wvtype, mshft, npts, douse, itotal_pts,
     1      nzyear,nzjday, nzhour, nzmin, nzsec, nzmsec
        integer wvtype(mxwvs), mshft(mxwvs), npts(mxwvs), douse(mxwvs)
        integer itotal_pts
        integer nzyear(mxwvs), nzjday(mxwvs), nzhour(mxwvs), 
     1      nzmin(mxwvs), nzsec(mxwvs), nzmsec(mxwvs)
        common/ctrace/obs_name, gftn_name, kstnm, kcmpnm, knetwk
        character*256 obs_name(mxwvs), gftn_name(mxwvs)
        character*8 kstnm(mxwvs), kcmpnm(mxwvs), knetwk(mxwvs)


        integer etype,  ibegin
        integer ounit, stdout

        real total_sum_sq, total_sum_sq_wt

c-----
c     common block for matrix inversion and solution
c-----
        common/mat/a,asave,datav,datasv,covmatrix,mtensor
        real a(nmax,nprm), asave(nmax,nprm)
        real datav(nmax), datasv(nmax)
        real covmatrix(nprm,nprm)
        real mtensor(6)

        real btime
        real deg_to_rad, rad_to_deg
        parameter(deg_to_rad = 0.017453292, rad_to_deg = 57.29577951)
        integer i, j, m
        real st_az
        integer ipass

        real varnoise, varsignal, deltao, varwt

        do 125 i = 1, nwaves
            call getwaveforms(gftn_name(i),obs_name(i),wvtype(i),
     1          npts(i),dt(i),btime,ounit,
     2          kstnm(i),kcmpnm(i),knetwk(i),disto(i),distg(i),
     3          depthg(i),az(i),user1(i),user2(i),
     4          tpg(i),tsg(i),bg(i),og(i),
     5          tpo(i),tso(i),bo(i),oo(i),
     6          deltao,
     6          nzyear(i),nzjday(i),nzhour(i),nzmin(i),
     7          nzsec(i), nzmsec(i))      
        call shiftall(mshft(i),npts(i),bg(i),bo(i),dt(i))
C       WRITE(0,*)kstnm(i),' ',mshft(i), npts(i)
c-----
c      compute the norm of the data
c-----
          etype = 2
          call getnorm(obs,npts(i),etype,sum_sq_amp(i))
c-----
c     NEW 06 APR 2006
c     get variance of window before the P and after the P
c     NOTE this assumes that there is a P time
c     the variance weighting is in the range of 
c          1  < S/N   < inf
c         0.05 < varwt < 1 
c     NOTE THIS IS COMPLETELY ADHOC AND NEETS MORE WORK
c     basically is the waveform is good visually 
c                   it should carry more weight
c     use sqrt variance since we really want this to sort of reject
c     S/N on order of 2 - hence the 0.25 term
c     also check for bad numbers
c-----
        if(deltao .gt.0.0)then
            j = 1 + (tpo(i) - bo(i))/deltao
            if(j.gt.1 .and. j.lt. npts(i)/2)then
                call getnorm(obs,j,etype,varnoise)
                call getnorm(obs(j),j+j,etype,varsignal)
                varwt = 1.0 - 1.95/(sqrt(0.25*varsignal/varnoise) + 1)
c-----
c     safety
c-----
                if(varwt.le.0.0)varwt = 0.001
            else
                varwt = 1.0
            endif
        else
            varwt = 1.0
        endif
        

c-----
c      SET UP WEIGHTING
c
c      if the input wt >= 0, 
c         then the weight is simply abs(wt(i)) (no sum_sq_amp
c         factor is included
c-----
        if(douse(i) .eq. 1)then
            obswt(i) = 1.0 * varwt
        else
            obswt(i) = 0.1 * varwt
        endif
c
C         write(ounit,*) ' Sum of Square Amplitudes: ',sum_sq_amp(i)
C         write(ounit,*) ' Weight:  ',obswt(i)
C         write(ounit,*) ' '
c-----
c
c      keep track of the total sum of square amplitudes to provide 
c         "normalized"
c       error measurements
c-----
          total_sum_sq = total_sum_sq + sum_sq_amp(i)
          total_sum_sq_wt = total_sum_sq_wt + (sum_sq_amp(i) * obswt(i))
c-----
c
c      SET-UP THE INVERSION
c
c      compute the response for individual mt elements
c      and store in the system matrix a
c-----
          if(ibegin .gt. nmax) then
            write(stdout,*) 'ERROR: Overflowing the a matrix'
            write(stdout,*) 'ibegin = ',ibegin
            write(stdout,*) 'max allowed =', nmax
            stop
          endif
c     
c-----
c      a copy of the "a matrix" is saved  (before weighting)
c         for later computation of errors and predicted fits
c         a quick alternative is to just read in the seismograms
c         and build the A matrix again later
c
c-----
c     in what follows below the unweighted Green and obs are saved in
c     asave and datasv respectively, the weighted are saved in
c     a and datav respectively
c-----
          st_az = az(i) * deg_to_rad
          do 100 m = 1, nprm
             call getmtresp(m,wvtype(i),st_az,npts(i))

             call appendvector(asave(1,m),nmax,ibegin,
     1                        response,npts(i))
             call appendvector(a(1,m),nmax,ibegin,
     1                        response,npts(i))
100     continue
c-----
c     FORCE DEVIATORIC
c     THIS IS DONE BY APPLYING THE CONSTRAINT THAT
c     m11 + m22 + m33 = 0 , or
c     mtensor(1) + mtensor(2) + mtensor(6) = 0
c
c     If this were not done, we would have a 
c                   general moment tensor inversion
c-----
        do 141 j=1,npts(i)
            a(j,1) = a(j,1) - a(j,6)
            a(j,2) = a(j,2) - a(j,6)
            a(j,6) = 0.0
  141   continue
            
c-----
c      next set up the data vector (save a copy here as well)
c----- 
          call appendvector(datasv,nmax,ibegin,obs,npts(i))
          call appendvector(datav ,nmax,ibegin,obs,npts(i))
c-----
c      apply the weighting
c-----
          if(obswt(i) .ne. 1.0) then
            call apply_weights(a,nmax,ibegin,npts(i),1,nprm,
     1          datav,obswt(i))
          endif
c-----
c      keep track of the number of rows in a matrix
c-----
          ibegin = ibegin + npts(i)
c-----
c      ibegin points to the next unused row of the a matrix
c-----
  125   continue
c-----
c   end loop over waveforms
c-----
        do 127 i = 1, nwaves
          itotal_pts = itotal_pts + npts(i)
127   continue
        write(ounit,*) ' '
        write(ounit,*)
     1    'Using a total of ',itotal_pts, ' points in inversion.'
        write(ounit,*) ' '
        return
        end

        subroutine domtfwd(ibegin,itotal_pts,ounit,m0,
     1     dip,rake,strike,xmom,isds,xmt)
        implicit none
        integer maxpts, nprm, nmax, mxwvs, isds
        parameter(maxpts=2048, nprm=6, mxwvs=1000, nmax =mxwvs*maxpts)
        integer ibegin
        real dip,rake,strike,xmom, xmt(3,3)
c-----
c     common block for matrix inversion and solution
c-----
        common/mat/a,asave,datav,datasv,covmatrix,mtensor
        real a(nmax,nprm), asave(nmax,nprm)
        real datav(nmax), datasv(nmax)
        real covmatrix(nprm,nprm)
        real mtensor(6)
c-----
c     common for best double couple and for solution vector
c-----
        common/rbstdc/stk0,dip0,rak0,stk1,dip1,rak1,xmw,pclvd,pred
        real pred(nmax)
        real stk0,dip0,rak0,stk1,dip1,rak1,xmw,pclvd

        real m0
        integer itotal_pts, ounit
        double precision amat(nmax,nprm), sval(nprm)
        double precision u(nmax,nprm), v(nmax,nprm)
        double precision tmp(nprm)
        double precision tsum
        integer nr, nc
        integer j, n, m, i, k
        integer ierr
        integer istk, idip, irake
c-----
c       create the moment tensor
c-----
        if(isds.eq.0)then
               istk = strike
               idip = dip
               irake = rake
               call dislocation_to_mij(istk,idip,irake,mtensor)
C       WRITE(6,*)'istk,idip,irake:',istk,idip,irake
C       WRITE(6,*)'mtensor:',mtensor
               do 2195 m=1,nprm
                   mtensor(m) = mtensor(m) * xmom
 2195          continue
        else if (isds.eq.1)then
c-----
c              from out_mte in mtsubs.f
c-----
               mtensor(1) = xmt(1,1)
               mtensor(2) = xmt(2,2)
               mtensor(3) = xmt(1,2)
               mtensor(4) = xmt(1,3)
               mtensor(5) = xmt(2,3)
               mtensor(6) = xmt(3,3)
        endif
c-----
c       compute the major/minor double couple from the moment tensor
c-----
        call mtdec(mtensor,m0,ounit,stk0,dip0,rak0,
     1      stk1,dip1,rak1,xmw,pclvd)
C       WRITE(6,*)'stk ,dip ,rake:',strike ,dip ,rake
C       WRITE(6,*)'xmom          :',xmom
C       WRITE(6,*)'stk0,dip0,rak0:',stk0,dip0,rak0
C       WRITE(6,*)'stk1,dip1,rak1:',stk1,dip1,rak1
C       WRITE(6,*)'xmw,pclvd:',xmw,pclvd
c-----
c    COMPUTE THE PREDICTED WAVEFORMS
c-----
c-----
c   compute the best-fitting predictions
c-----
        do 1280 j = 1, itotal_pts
            pred(j) = 0
 1280   continue



        do 2300 m = 1, nprm
            do 2290 j = 1, itotal_pts
                pred(j) = pred(j) + mtensor(m)*asave(j,m)
 2290       continue
 2300   continue
        return
        end

        subroutine dosumup(nwaves,ounit,m0,
     1      total_sum_sq, total_sum_sq_wt, 
     4      ipass)
        implicit none
        integer maxpts, nprm, nmax, mxwvs
        parameter(maxpts=2048, nprm=6, mxwvs=1000, nmax =mxwvs*maxpts)

c-----
c     common block for matrix inversion and solution
c-----
        common/mat/a,asave,datav,datasv,covmatrix,mtensor
        real a(nmax,nprm), asave(nmax,nprm)
        real datav(nmax), datasv(nmax)
        real covmatrix(nprm,nprm)
        real mtensor(6)
c-----
c     common blocks for traces and trace attributes
c-----
        common/rtrace/az,wt,user1,user2,
     1      disto, depthg, dt, distg,
     2      tpg, tsg, bg, og,
     3      tpo, tso, bo, oo,
     4      sum_sq_amp, fr_errors, obswt, misfit
        real az(mxwvs), wt(mxwvs), user1(mxwvs), user2(mxwvs)
        real disto(mxwvs), depthg(mxwvs), dt(mxwvs), distg(mxwvs)
        real tpg(mxwvs), tsg(mxwvs), bg(mxwvs), og(mxwvs)
        real tpo(mxwvs), tso(mxwvs), bo(mxwvs), oo(mxwvs)
        real sum_sq_amp(mxwvs), fr_errors(mxwvs), misfit(mxwvs)
        real obswt(mxwvs)
        common/itrace/wvtype, mshft, npts, douse, itotal_pts,
     1      nzyear,nzjday, nzhour, nzmin, nzsec, nzmsec
        integer wvtype(mxwvs), mshft(mxwvs), npts(mxwvs), douse(mxwvs)
        integer itotal_pts
        integer nzyear(mxwvs), nzjday(mxwvs), nzhour(mxwvs), 
     1      nzmin(mxwvs), nzsec(mxwvs), nzmsec(mxwvs)
        common/ctrace/obs_name, gftn_name, kstnm, kcmpnm, knetwk
        character*256 obs_name(mxwvs), gftn_name(mxwvs)
        character*8 kstnm(mxwvs), kcmpnm(mxwvs), knetwk(mxwvs)
c-----
c     common for best double couple and for solution vector
c-----
        common/rbstdc/stk0,dip0,rak0,stk1,dip1,rak1,xmw,pclvd,pred
        real pred(nmax)
        real stk0,dip0,rak0,stk1,dip1,rak1,xmw,pclvd

c

        real total_sum_sq, total_sum_sq_wt
        real m0

        integer ipass


        integer nwaves, ounit
        real total_err, total_wtd
        real err, factor
        real xx, xy, yy
        real tmpsig
        integer ibegin, i, n
        integer slen

        integer iPz, iPr, iSH
        real xmis
c-----
c     function prototype
c-----
        real vdot
        integer lgstr
c-----
c     get the fit between observed and predicted
c-----
        xx = vdot(datasv, datasv, itotal_pts)
        xy = vdot(datasv, pred  , itotal_pts)
        yy = vdot(pred  , pred  , itotal_pts)
c-----
c     write the observed .obs and predicted .pre waveforms
c-----
        iPz = 1
        iPr = 2
        iSH = 3
        total_err = 0.0
        total_wtd = 0.0
        ibegin = 1
        n = 0
        write(ounit,1)
    1   format('Station  Comp      Wt   Sum_Sq_Err ',
     1      'Frac  Factor    Misfit   TimeShift')
        do 350 i = 1, nwaves
c-----
c         write out the prediction as a sac file
c-----
            if(ipass.eq.2 )then
                call preobssv(i,ibegin)
            endif
c-----
c      COMPUTE THE ERRORS (scale factor = 1.0, l2 norm)
c
c
c WARNING!!!!!  There are several ways to compute 
c the error or misfit between
c  the observed and synthetic seismograms.  
c This part of the program has been
c  changes from the original code of Ammon, so please beware.  
c I have modified
c  the definition used to compute the misfit.  
c The formulation that I use is
c  misfit = ( data - syn )**2 / data**2, so if data = syn, misfit = 0;
c  if data = -syn, misfit = 4; if syn = 0, misfit = 1.  
c This is similar to
c  the error function of Pasyanos et al. (1996).  
c This is the error that is
c  passed through to the subroutine the prints the ASCII solution.
c
c-----

            err = 0.0

    2   format(a8,1x,a8,1x,f6.3,e10.3,f6.3,e10.3,e10.3,f10.2)
C            write(ounit,2)kstnm(i), kcmpnm(i),obswt(i),
C     1          err,err/sum_sq_amp(i),factor,misfit(i)
C     2          ,mshft(i)*dt(i)
c 
            ibegin = ibegin + npts(i)
c
            fr_errors(i) = err/sum_sq_amp(i)
c
  350   continue
c    end loop over waveforms
c
c
c*******************************************************************
c    Graphical error summary
c
C       call bar_chart(nwaves,obs_name,fr_errors,12,50,ounit)
c
c***********************************************************************
c
c    OUTPUT the error for all waveforms
c
c    output the sum of the square errors divided by the number of traces
c      minus the number of unknowns (free parameters)
c
        tmpsig = 1.0 / float(n-6)
c
        write(ounit,*) ' '
        write(ounit,*)
     1 'Unweighted SSE,Weighted SSE,Unweighted RMS,Weighted RMS'
        write(ounit,*) total_err*tmpsig, total_wtd*tmpsig,
     1      sqrt(total_err*tmpsig), sqrt(total_wtd*tmpsig)
c
c Output sum of square errors divided by the 
c sum of the square amplitudes.
c
        write(ounit,*)
     1 'Total Error divided by sum of all seismogram amplitudes'
        write(ounit,*)
     1 'Unweighted, weighted'
        write(ounit,*) total_err/total_sum_sq, total_wtd/total_sum_sq_wt
c-----
c     Now output the Moment Tensor Solution Again but with
c     Error estimates
c-----
        call out_mte(mtensor,total_err*tmpsig,covmatrix,nprm,m0,ounit)
        close (3)
        return
        end

        subroutine shiftall(mshft,npts,bg,bo,dt)
        integer mshft, npts
c-----
c     apply the shift to the data stream - also change 
c     the number of points
c     This would be simpler if we were to use new arrays but to save
c     space we will overwrite arrays
c
c     mshft - shift in samples: =0 no shift
c                   >0 shift observed to left, change
c                     npts and bo
c                   <0 shift green to left, change
c                     npts and bg
c                  |=0 shorten length of time series
c     npts - number of time samples in series - returned as
c             npts - abs(mshft)
c     bg   - SAC header to time stamp of B for first sample for Green
c     bo   - SAC header to time stamp of B for first sample for obs
c     dt   - sample interval
c-----
        integer maxpts, nprm, nmax, mxwvs
        parameter(maxpts=2048, nprm=6, mxwvs=1000, nmax =mxwvs*maxpts)

        real vss, vds, vdd, vex, response
        common/gftns/vss(maxpts),vds(maxpts),vdd(maxpts),vex(maxpts),
     1                response(maxpts)
        common /observations/ obs(maxpts)
        real obs
        integer i, j, n
        n = npts
        if(mshft.eq.0)then
            return
        else if(mshft.lt.0)then
            npts = npts + mshft
            bg = bg - mshft*dt
            do i=1,npts
                j = i - mshft
                if(j.ge.1 .and. j.le.n)then
                    vss(i) = vss(j)
                    vds(i) = vds(j)
                    vdd(i) = vdd(j)
                    vex(i) = vex(j)
                else
                    vss(i) = vss(npts)
                    vds(i) = vds(npts)
                    vdd(i) = vdd(npts)
                    vex(i) = vex(npts)
                endif   
            enddo
        else if(mshft.gt.0)then
            npts = npts - mshft
            bo = bo + mshft*dt
            do i=1,npts
                j = i + mshft
                if(j.ge.1 .and. j.le.n)then
                    obs(i) = obs(j)
                else
                    obs(i) = obs(npts)
                endif   
            enddo
        endif
c-----
c     safety
c-----
        if(npts.lt.0)npts = 0
        return
        end

        subroutine preobssv(i,ibegin)
        implicit none
        integer i

        integer maxpts, nprm, nmax, mxwvs
        parameter(maxpts=2048, nprm=6, mxwvs=1000, nmax =mxwvs*maxpts)
c-----
c     common block for matrix inversion and solution
c-----
        common/mat/a,asave,datav,datasv,covmatrix,mtensor
        real a(nmax,nprm), asave(nmax,nprm)
        real datav(nmax), datasv(nmax)
        real covmatrix(nprm,nprm)
        real mtensor(6)
c-----
c     common blocks for traces and trace attributes
c-----
        common/rtrace/az,wt,user1,user2,
     1      disto, depthg, dt, distg,
     2      tpg, tsg, bg, og,
     3      tpo, tso, bo, oo,
     4      sum_sq_amp, fr_errors, obswt, misfit
        real az(mxwvs), wt(mxwvs), user1(mxwvs), user2(mxwvs)
        real disto(mxwvs), depthg(mxwvs), dt(mxwvs), distg(mxwvs)
        real tpg(mxwvs), tsg(mxwvs), bg(mxwvs), og(mxwvs)
        real tpo(mxwvs), tso(mxwvs), bo(mxwvs), oo(mxwvs)
        real sum_sq_amp(mxwvs), fr_errors(mxwvs), misfit(mxwvs)
        real obswt(mxwvs)
        common/itrace/wvtype, mshft, npts, douse, itotal_pts,
     1      nzyear,nzjday, nzhour, nzmin, nzsec, nzmsec
        integer wvtype(mxwvs), mshft(mxwvs), npts(mxwvs), douse(mxwvs)
        integer itotal_pts
        integer nzyear(mxwvs), nzjday(mxwvs), nzhour(mxwvs), 
     1      nzmin(mxwvs), nzsec(mxwvs), nzmsec(mxwvs)
        common/ctrace/obs_name, gftn_name, kstnm, kcmpnm, knetwk
        character*256 obs_name(mxwvs), gftn_name(mxwvs)
        character*8 kstnm(mxwvs), kcmpnm(mxwvs), knetwk(mxwvs)
c-----
c     common for best double couple and for solution vector
c-----
        common/rbstdc/stk0,dip0,rak0,stk1,dip1,rak1,xmw,pclvd,pred
        real pred(nmax)
        real stk0,dip0,rak0,stk1,dip1,rak1,xmw,pclvd
c
        integer ibegin, nerr
        integer slen

        integer iSH

        character wvfile*256
c-----
c     function prototype
c-----
        integer lgstr
c-----
c     save the predicted abd observed traces
c-----
c         set up predicted trace for this observation
c         but design it so that it permits an overlay of
c         observed and predicted
c
c         convert from dist(km) to degrees using conversion
c         constant build within GSAC
c-----
            call newhdr()
            call setnhv('NZYEAR  ',nzyear(i),nerr)
            call setnhv('NZJDAY  ',nzjday(i),nerr)
            call setnhv('NZHOUR  ',nzhour(i),nerr)
            call setnhv('NZMIN   ',nzmin(i) ,nerr)
            call setnhv('NZSEC   ',nzsec(i) ,nerr)
            call setnhv('NZMSEC  ',nzmsec(i),nerr)
            call setkhv('KCMPNM  ',kcmpnm(i),nerr)
            call setkhv('KSTNM   ',kstnm (i),nerr)
            call setkhv('KNETWK  ',knetwk(i),nerr)
            call setfhv('EVDP    ',depthg (i),nerr)
            call setfhv('DIST    ',disto (i),nerr)
            call setfhv('GCARC   ',disto (i)/111.1946,nerr)
            call setfhv('AZ      ',az    (i),nerr)
            call setfhv('BAZ     ',-12345.  ,nerr)
            call setfhv('USER1   ',user1 (i),nerr)
            call setfhv('USER2   ',user2 (i),nerr)
            call setfhv('O       ',og    (i)-og(i),nerr)
            call setfhv('USER9   ', mshft(i)*dt(i),nerr)
            if(wvtype(i).eq.iSH)then
                call setfhv('T1      ',tsg   (i),nerr)
            else
                call setfhv('T0      ',tsg   (i),nerr)
            endif
        
            wvfile = ' '
            slen = lgstr(obs_name(i))
            wvfile(1:slen+4) = obs_name(i)(1:slen)//'.pre'
            if(douse(i).eq.1)then
c-----
c     if we used the trace, we will output the observed and predicted
c     we also take time to ensure that the P and S picks are
c     preserved for comparison
c-----
            if(mshft(i).gt.0)then
            call setfhv('A       ',tpg(i)-og(i)+mshft(i)*dt(i),nerr)
            if(tso(i).gt. -12345.0 .and. tsg(i).gt.-12345.)then
            if(wvtype(i).eq.iSH)then
            call setfhv('T1      ',tsg(i)-og(i)+mshft(i)*dt(i),nerr)
            else
            call setfhv('T0      ',tsg(i)-og(i)+mshft(i)*dt(i),nerr)
            endif
            endif
            call wsac1(wvfile,pred(ibegin),npts(i),bo(i),dt(i),nerr)
            else

            call setfhv('A       ',bo(i)+tpg(i)-bg(i),nerr)
            if(tso(i).gt. -12345.0 .and. tsg(i).gt.-12345.)then
            if(wvtype(i).eq.iSH)then
            call setfhv('T1      ',bo(i)+tsg(i)-bg(i),nerr)
            else
            call setfhv('T0      ',bo(i)+tsg(i)-bg(i),nerr)
            endif
            endif
            call wsac1(wvfile,pred(ibegin),npts(i),bo(i),dt(i),nerr)
            endif
            endif
c------
c     output the observed
c-----
            call setfhv('A       ',tpo(i),nerr)
            call setfhv('B       ',bo(i),nerr)
            call setfhv('O       ',oo(i),nerr)
            call setfhv('USER9   ', -12345.,nerr)
            wvfile = ' '
            slen = lgstr(obs_name(i))
            wvfile(1:slen+4) = obs_name(i)(1:slen)//'.obs'
            if(douse(i).eq.1)then
            call wsac1(wvfile,datasv(ibegin),npts(i),bo(i),dt(i),nerr)
            endif
        return
        end

        subroutine  makecrack(x,strike,dip,rake,xmom)
        real strike, dip, rake, x(3,3), xmom
c-----
c       dip     R*4      - dip of crack with respect to the horizontal
c       strike  R*4      - strike of crack plane - the plane dips
c                          down to the right
c       rake    R*4      - > 0 opening crack, moment > 0
c                          < 0 closing crack, moment < 0
c        x      R*4      - moment tensor
c-----
c       NOTE THIS ASSUMES THAT lambda = mu, isotropic medium with
c           Poisson ratio = 0.25
c-----
c-----
c       special handling for the crack source. convert it to a moment tensor
c       Nakano, M., and H. Kumagi (2005). Waveform inversion of volcano-seismic 
c            signals assuming possible source geometries,
c            Geophys. Res. Letters 32, L12302, doi:10.1029/2005GL022666.
c       the conversion to moment tensor must be deferred until later in the code
c       since the lambda and mu are required at the source depth
c
c       Note they define a theta as the angle that the plane normal makes with
c       respect to the vertical. By definition, theta = dip
c-----
       real degrad, cs,ss, xm, lamdmu
       real ct, st, s2t, s2s
       degrad = 3.1415927/180.0
       cs = cos(degrad*strike)
       ss = sin(degrad*strike)
       ct = cos(degrad*dip   )
       st = sin(degrad*dip   )

       s2s = 2.*ss*cs
       s2t = 2.*st*ct
c-----
c      opening crack
c-----
       if(rake.gt.0.0)then
              xm = xmom
c-----
c      closing crack
c-----
       else if(rake .le. 0.0)then
              xm = - xmom
       endif
c-----
c      lambda/mu
c-----
       lamdmu = 1.0
       x(1,1) = xm * ( lamdmu + 2.* st*st*cs*cs)
       x(1,2) = xm * ( st*st*s2s )
       x(1,3) = xm * ( s2t*cs )
       x(2,2) = xm * ( lamdmu + 2.*st*st*ss*ss )
       x(2,3) = xm * ( s2t * ss)
       x(3,3) = xm * ( lamdmu + 2.*ct*ct)
       x(2,1) = x(1,2)
       x(3,1) = x(1,3)
       x(3,2) = x(2,3)
       return 
       end
