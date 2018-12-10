       program wvfgrd96
c---------------------------------------------------------------------c
c                                                                     c
c      COMPUTER PROGRAMS IN SEISMOLOGY                                c
c      VOLUME VIII                                                    c
c                                                                     c
c      PROGRAM: WVFGRD96                                              c
c                                                                     c
c      COPYRIGHT 1994 Charles J. Ammon                                c
c                                                                     c
c      Department of Geosciences                                      c
c      Penn State University                                          c
c      Deike Bldg                                                     c
c      University Park, PA 16802                                      c
c      U. S. A.                                                       c
c                                                                     c
c---------------------------------------------------------------------c
c       FIXES
c       21 OCT 2002 - corrections to SAC header for prediction
c      COMPUTER PROGRAMS IN SEISMOLOGY                                c
c       16 SEP 2004 - the predicted have the correct year day 
c           hr min sec msec
c       02 NOV 2004 - implement time shift - also implement 
c           precomputation
c               of dot products to speed up computations a low
c       16 FEB 2006 - output the observed file as an .obs
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
c       19 MAR 2007 - implemented a two pass - to permit a fast 
c         crude followed by a
c         refined search to save time
c       30 AUG 2007 - the observed and predicted traces now output
c         have the time shift set in the header in USER9 and, more
c         importantly, have the traces shifted for alignment
c	05 JAN 2009 - ensure that KNETWK is set in the .obs and .pre files
c       05 SEP 2010 - the predicted has USER5 set which is the percentage
c                     of set, e.g., 100 * correlation coefficient-squared     
c-----
c
c      dislocation source grid search
c
c      Time function is fixed, so problem is linear
c
c      Author:  Chuck Ammon, UC Santa Cruz
c               parts by George Randall, U South Carolina
c
c      Version: 1.1 September, 1994
c
c      Notation from Langston (1981)
c        x is north, y is east, z is down
c        all signs are opposite Aki and Richards due to a different
c        choice of the direction for the fault normal
c-----
        implicit none
        integer maxpts, nprm, NMAX, mxwvs
        parameter(maxpts=2048, nprm=6, mxwvs=1000, NMAX =mxwvs*maxpts)
        real DEG_TO_RAD, RAD_TO_DEG
        parameter(DEG_TO_RAD = 0.017453292, RAD_TO_DEG = 57.29577951)

        common/gftns/vss(maxpts),vds(maxpts),vdd(maxpts),vex(maxpts),
     1      response(maxpts)
        real vss, vds, vdd, vex, response
        common /observations/ obs(maxpts)
        real obs

        common/event/evla, evlo
        real evla, evlo

        real pred(maxpts*mxwvs)
        real a(NMAX,nprm),  asave(NMAX,nprm)
        real datav(NMAX), datasv(NMAX)
        real mtensor(nprm), btime, az(mxwvs), wt(mxwvs)
        real user1(mxwvs), user2(mxwvs)
        real sum_sq_amp(mxwvs)
        real obswt(mxwvs), err_l2, min_err_l2
        real rvar_l2
        real vwt, rwt, twt
        real total_sum_sq, total_sum_sq_wt
        real st_az
        real factor
        real total_err, total_wtd
        real best_l2_factor
        real disto(mxwvs), depth(mxwvs), dt(mxwvs), distg(mxwvs)
        real tpg(mxwvs), tsg(mxwvs), bg(mxwvs), og(mxwvs)
        real tpo(mxwvs), tso(mxwvs), bo(mxwvs), oo(mxwvs)
        real xmwl2
        real best_l2_mw
        real rr(3)
        real depmax, depmin
        integer ixx,iyy,ixy,ixz,iyz,izz, wvtype(mxwvs), iPz, iPr, iSH
        integer i, nerr, etype, ibegin, npts(mxwvs)
        integer nzyear(mxwvs), nzjday(mxwvs), nzhour(mxwvs), 
     1      nzmin(mxwvs), nzsec(mxwvs), nzmsec(mxwvs)
        integer slen, lgstr
        integer itotal_pts, best_l2(3)
        integer stdin,stdout,ounit,inunit
        integer nwaves
        integer j, m, n
        integer istrike, idip, irake
        integer jrake
        integer mpts
        integer ipass
        integer istkmn, istkmx, istkdd
        integer idipmn, idipmx, idipdd
        integer irakmn, irakmx, irakdd

       
        character*256 obs_name(mxwvs),gftn_name(mxwvs),infile,ofile
        character*259 wvfile
        character*8 kstnm(mxwvs), kcmpnm(mxwvs), knetwk(mxwvs)
c-----
c       New Additions R.B. Herrmann, Saint Louis University
c           28 Oct 2004 - first modifications to permit time shift
c               initially just precompute crosscorrelations
c           
c-----
c       yy = autocorrelation of each trace with trace weighted by w
c       yg = cross-correlation between Green and trace
c       gg = cross-correlation among the Green
c-----
c       We will do a time shift for a total of nshft points
c       these will be defined by klw ,..., kup
c       e.g., nsfht = 3, (odd number always) klw = -1 kup = 1
c           k=1,nshft  1  2  3  4  5  6  7  8  9
c           shift     -4 -3 -2 -1  0  1  2  3  4
c           shift =   k - 5  where 5 = nshft/2 +1 = koff
c-----
        integer KSHFT
        parameter (KSHFT=800)
        real gg(mxwvs,nprm,nprm)
        real vdot, vdots

        real yy(mxwvs)
        real yg(mxwvs,nprm,KSHFT)

        real mijg(maxpts,nprm)
        real oodot, opdot, ppdot
        real too, top, tpp
        real oor(3), opr(3), ppr(3)
        real ttop
        integer k
        integer kshift(mxwvs)
        integer mshft(mxwvs)
        integer nshft, nsft, koff
        real deltao
        double  precision  rxx, rxy, ryy, xxyymx, pctfit
        

        
c
        stdin = 5
        stdout = 6
        ounit = 10
        inunit = 11
c
        ixx = 1
        iyy = 2
        ixy = 3
        ixz = 4
        iyz = 5
        izz = 6
c
        vwt = 1
        rwt = 1
        twt = 1
c
        iPz = 1
        iPr = 2
        iSH = 3
c
        ibegin = 1
c
        total_sum_sq = 0
        total_sum_sq_wt = 0
c
        itotal_pts = 0

        min_err_l2 = 1e9
c-----
c       parse the command line
c-----
        call gcmdln(nsft)
c-----
c       safety check 
c       We must not exceed the dimension of KSHFT. We also
c       want an odd number to search.
c-----
        if(nsft .lt. 0)nsft = 0
        if (nsft .gt. KSHFT/2)then
            nsft = KSHFT/2
        endif
        koff = nsft + 1
        nshft = nsft + nsft + 1
        
c
        write(stdout,*)'What is the input file name?'
        read(stdin,*) infile
c
        write(stdout,*)'What is the output file name?'
        read(stdin,*) ofile
        open(unit = ounit, file = ofile)
        rewind(ounit)
        write(ounit,*)'******* Dislocation Grid Search ************'
c
        open(unit = inunit, file = infile)
        rewind(inunit)
        i = 1
        nwaves = 0
    5   continue
            read(inunit,*,end=6)wvtype(i),obs_name(i),gftn_name(i),
     1          wt(i)
            nwaves = i 
            if(abs(wt(i)) .eq.0.0)wt(i) = 1.0e-8
c-----
c       protection against exceeding array dimension
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
            if(depmax.eq.depmin)then
                  WRITE(0,*)'Rejecting: ',i,' DEPMAX=',depmax,
     1            ' DEPMIN=',depmin,obs_name(i)
                  go to 5
            endif
c-----
c           OK this is a valid trace set
c-----
            i = i + 1
        go to 5
    6   continue
        close(inunit)
c
c
c***********************************************************************
c      Output some event and inversion 
c           information (using the first data file)
c    
       call getevent_data(obs_name(1),ounit,'wvfgrd96')
c
c***********************************************************************
c
c      Now read in the data, GFs, and set up the matrix equations
c       
        do 125 i = 1, nwaves

            call getwaveforms(gftn_name(i),obs_name(i),wvtype(i),
     1          npts(i),dt(i),btime,ounit,
     2          kstnm(i),kcmpnm(i),knetwk(i),disto(i),distg(i),
     3          depth(i),az(i),user1(i),user2(i),
     4          tpg(i),tsg(i),bg(i),og(i),
     6          tpo(i),tso(i),bo(i),oo(i),
     6          deltao,
     4          nzyear(i),nzjday(i),nzhour(i),nzmin(i),
     5          nzsec(i), nzmsec(i))      
            st_az = az(i) * DEG_TO_RAD
c-----
c           compute the norm of the data
c-----
            etype = 2
            call getnorm(obs,npts(i),etype,sum_sq_amp(i))
c-----
c       form weighted dot products
c-----
            do 101 m = 1, nprm
                call getmtresp(m,wvtype(i),st_az,npts(i))
                call appendvector(mijg(1,m),maxpts,1,
     1              response,npts(i))
  101       continue
            yy(i) = vdot(obs,obs,npts(i))*abs(wt(i))
            do 1002 m=1,nprm
                do 1004 k=1,nshft
                yg(i,m,k) = vdots(obs,mijg(1,m),npts(i),k-koff)
     1              *abs(wt(i))
 1004           continue
C               yg(i,m,1) = vdot(obs,mijg(1,m),npts(i))*abs(wt(i))
                do 1003 n=1,nprm
                    gg(i,m,n)= vdot(mijg(1,m),mijg(1,n),
     1                  npts(i))*abs(wt(i))
 1003           continue
 1002       continue
            
c-----
c           SET UP WEIGHTING
c
c           if the input wt >= 0, 
c           then the weight is simply abs(wt(i)) (no sum_sq_amp
c           factor is included
c-----
            if(wt(i) .ge. 0.0) then
                obswt(i) = abs(wt(i))
            endif

            write(ounit,*) ' Sum of Square Amplitudes: ',sum_sq_amp(i)
            write(ounit,*) ' Weight:  ',obswt(i)
            write(ounit,*) ' '
c-----
c
c           keep track of the total sum of square amplitudes to provide 
c           "normalized"
c           error measurements
c-----
            total_sum_sq = total_sum_sq + sum_sq_amp(i)
            total_sum_sq_wt = total_sum_sq_wt + 
     1              (sum_sq_amp(i) * obswt(i))
c-----
c
c           SET-UP THE INVERSION
c
c           compute the response for individual mt elements
c           and store in the system matrix a
c-----
            if(ibegin .gt. NMAX) then
                write(stdout,*) 'ERROR: Overflowing the a matrix'
                write(stdout,*) 'ibegin = ',ibegin
                write(stdout,*) 'max allowed =', NMAX
                stop
            endif
c-----
c           a copy of the "a matrix" is saved  (before weighting)
c           for later computation of errors and predicted fits
c               a quick alternative is to just read in the seismograms
c           and build the A matrix again later
c-----
            do 100 m = 1, nprm
                call getmtresp(m,wvtype(i),st_az,npts(i))
                call appendvector(a(1,m),NMAX,ibegin,
     1              response,npts(i))
                call appendvector(asave(1,m),NMAX,ibegin,
     1              response,npts(i))
  100       continue
c-----
c           next set up the data vector (save a copy here as well)
c-----
            call appendvector(datav ,NMAX,ibegin,obs,npts(i))
            call appendvector(datasv,NMAX,ibegin,obs,npts(i))
c-----
c           apply the weighting
c-----
            if(obswt(i) .ne. 1.0) then
                call apply_weights(a,NMAX,ibegin,npts(i),
     1              1,nprm,datav,sqrt(obswt(i)))
            endif
c-----
c           keep track of the number of rows in a matrix
c-----
            ibegin = ibegin + npts(i)
c-----
c           ibegin points to the next unused row of the a matrix
c-----
  125   continue
c-----
c       end loop over waveforms
c-----
        do 127 i = 1, nwaves
            itotal_pts = itotal_pts + npts(i)
  127   continue

        write(ounit,*) ' '
        write(ounit,*)
     1      'Using a total of ',itotal_pts, ' points in inversion.'
        write(ounit,*) ' '
c-----
c       END of SETUP 
c-----
       write(ounit,*) 
     1  'Strike Dip Rake Mw Reduction_Variance'

c-----
c       LOOP over strike dip and rake
c       for each set of angles, compute the corresponding moment tensor
c       and then compute the predicted responses for that Mij
c-----

        rvar_l2 = -2.
        do 1001 ipass=1,2
        if(ipass.eq.1)then
            istkmn = 0
            istkmx = 355
            istkdd = 10
            idipmn = 0
            idipmx = 90
            idipdd = 10
            irakmn = -90
            irakmx = 90
            irakdd = 10
        else if(ipass.eq.2)then
            istkdd = 5
            istkmn = best_l2(1) - 4.*istkdd
            istkmx = best_l2(1) + 4.*istkdd
            idipdd = 5
            idipmn = best_l2(2) - 4.*idipdd
            idipmx = best_l2(2) + 4.*idipdd
            if(idipmx.gt.90)then
                 idipmx=90
                 idipmn = idipmn - 4.*idipdd
            endif
            if(idipmn.lt. 0)then
                 idipmx= 0
                 idipmx = idipmx + 4.*idipdd
            endif
            irakdd = 5
            irakmn = best_l2(3) - 4.*irakdd
            irakmx = best_l2(3) + 4.*irakdd
        endif

        do 1000 istrike = istkmn, istkmx, istkdd
            do  900 idip    = idipmn, idipmx, idipdd
                do  800 irake   = irakmn, irakmx, irakdd

                call dislocation_to_mij(istrike,idip,
     1              irake,mtensor)
C
c-----
c       apply new methodology
c-----
c       we need to compute the o.o p.p o.p from which we get the
c       moment and the correlation coefficient
c-----
c       eventually we will check the cross correlation to 
c           get the lag time
c-----
            oodot = 0.0
            opdot = 0.0
            ppdot = 0.0
            oor(1) = 0.0
            oor(2) = 0.0
            oor(3) = 0.0
            opr(1) = 0.0
            opr(2) = 0.0
            opr(3) = 0.0
            ppr(1) = 0.0
            ppr(2) = 0.0
            ppr(3) = 0.0
        do 3125 i = 1, nwaves
                too = yy(i)
                tpp = 0.0
                do 3300 m = 1,nprm
                do 3301 n = 1,nprm
                    tpp = tpp + gg(i,m,n)*mtensor(m)*
     1                  mtensor(n)
 3301           continue
 3300           continue
c-----
c       here we search to find the time shift corresponding to the
c       best fit
c       We just search the "o dot p" for the maximum value
c-----
                top = -1.0e+38
                do 3302 k = 1, nshft
                    ttop = 0.0
                    do 3303 m = 1, nprm
                    ttop = ttop + yg(i,m,k)*mtensor(m)
 3303               continue
                    if(ttop.gt.top)then
                        top = ttop
                        kshift(i) = k - koff
                    endif
 3302           continue
                oodot = oodot + too
                opdot = opdot + top
                ppdot = ppdot + tpp
                j = wvtype(i)
                oor(j) = oor(j) + too
                opr(j) = opr(j) + top
                ppr(j) = ppr(j) + tpp
                    
 3125   continue
        do 3126 i=1,3
        if(oor(i).eq.0.0 .or. ppr(i).eq.0.0)then
            rr(i) = -2.0
        else
            rr(i) = opr(i)/sqrt(oor(i)*ppr(i))
        endif
 3126   continue
        factor = (opdot/ppdot)
        if(factor.lt.0.0)then
            rr(1) = - rr(1)
            rr(2) = - rr(2)
            rr(3) = - rr(3)
            factor = abs(factor)
            if(irake.lt.0)then
                jrake = irake + 180
            else
                jrake = irake - 180
            endif
        else
            jrake = irake
        endif
        xmwl2= (20.0+ alog10(factor) - 16.10)/1.5
        write(ounit,3)depth(1),istrike,idip,jrake,xmwl2,
     1      opdot*opdot/(oodot*ppdot),rr
    3   format('WVFGRD96 ',f6.1,3i6,f7.2,f7.4,3f8.4)
c
c-----
c       define the best mechanism and also save the corresponding time
c       shift
c-----
                if(opdot*opdot/(oodot*ppdot) .gt. rvar_l2)then
                    rvar_l2 = opdot*opdot/(oodot*ppdot)
                    min_err_l2 = total_sum_sq_wt*(1. - rvar_l2)
                    best_l2_factor = factor
                    best_l2(1) = istrike
                    best_l2(2) = idip
                    best_l2(3) = jrake
                    best_l2_mw = xmwl2
                    do 801 i=1,nwaves
                        mshft(i) = kshift(i)
  801               continue
                endif
  800           continue
  900       continue
 1000   continue
 1001   continue
c-----
c      FINISHED WITH THE GRID SEARCH
c      
c      output info about the best fitting L2 dislocation
c-----
        write(ounit,1500) 'best l2: ',
     1      best_l2(1),best_l2(2),best_l2(3),best_l2_factor
c-----
c       compute the best fitting moment tensor
c-----
        call dislocation_to_mij(best_l2(1),best_l2(2),
     1          best_l2(3),mtensor)
        do 1270 i = 1,nprm
            mtensor(i) = mtensor(i) * best_l2_factor
 1270   continue
c-----
c       compute the best-fitting predictions
c-----
        do 1280 j = 1, itotal_pts
            pred(j) = 0
 1280   continue
        do 2300 m = 1, nprm
            do 2290 j = 1, itotal_pts
                pred(j) = pred(j) + mtensor(m)*asave(j,m)
 2290       continue
 2300   continue
c-----
c     write out the prediction as a sac file
c-----
        ibegin = 1
        do 1400 i = 1, nwaves
c-----
c              apply the shift by padding with the last 
c              point if necessary
c-----
            call doshift(mshft(i),pred(ibegin),bg(i),datasv(ibegin),
     1            bo(i),npts(i),dt(i))
c-----
c     get correlation coefficient of fit between
c     observed and predictedy
c-----
            rxy = vdot(datasv(ibegin),pred(ibegin)  ,npts(i))
            ryy = vdot(pred(ibegin)  ,pred(ibegin)  ,npts(i))
            rxx = vdot(datasv(ibegin),datasv(ibegin),npts(i))
            xxyymx = amax1(rxx,ryy)
            pctfit=100.0*(rxy*rxy)/(rxx*ryy + 0.0001*xxyymx*xxyymx)
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
            call setfhv('DIST    ',distg (i),nerr)
            call setfhv('EVDP    ',depth (i),nerr)
            call setfhv('AZ      ',az    (i),nerr)
            call setfhv('BAZ     ',-12345.  ,nerr)
            call setfhv('USER1   ',user1 (i),nerr)
            call setfhv('USER2   ',user2 (i),nerr)
            call setfhv('A       ',tpg   (i)-og(i),nerr)
            call setfhv('B       ',bg    (i)-og(i),nerr)
            call setfhv('O       ',og    (i)-og(i),nerr)
            call setfhv('USER5   ',sngl(pctfit),nerr)
            call setfhv('USER9   ',mshft(i)*dt(i),nerr)
            if(wvtype(i).eq.iSH)then
                call setfhv('T1      ',tsg   (i),nerr)
            else
                call setfhv('T0      ',tsg   (i),nerr)
            endif
            wvfile = ' '
            slen = lgstr(obs_name(i))
            wvfile(1:slen+4) = obs_name(i)(1:slen)//'.pre'     
            call wsac1(wvfile,pred(ibegin),npts(i),
     1          bg(i)-og(i),dt(i),nerr)
c------
c       output the observed
c-----
            call setfhv('A       ',tpo   (i),nerr)
            call setfhv('B       ',bo    (i),nerr)
            call setfhv('O       ',oo    (i),nerr)
            call setfhv('USER5   ', -12345.,nerr)
            call setfhv('USER9   ', -12345.,nerr)
            wvfile = ' '
            slen = lgstr(obs_name(i))
            wvfile(1:slen+4) = obs_name(i)(1:slen)//'.obs'
            call wsac1(wvfile,datasv(ibegin),npts(i),bo(i),dt(i),nerr)
            ibegin = ibegin + npts(i)
 1400   continue

            close(ounit)
 1500   format(a9,1x,i4,1x,i2,1x,i4,f10.3)
c
c-----
c       output L2 summary file
c-----
        open(3,file='fmdfit.dat',access='sequential',
     1      form='formatted',status='unknown')
        rewind 3
c-----
c       ASSUME GREEN s FUNCTIONS ARE IN SAME PHYSICAL UNITS AS
c       OBSERVED - but that Greens are for a Moment of 1.0e+20 dyne cm
c-----
        
        write(3,3)depth(1), best_l2(1),best_l2(2),best_l2(3),
     1      best_l2_mw,rvar_l2
        close (3)
        end
c-----
c       END of MAIN PROGRAM
c-----

        subroutine gcmdln(nsft)
        parameter (LER=0, LIN=5, LOT=6)
c-----
c       parse command line arguments
c
c       requires subroutine mgtarg() and function mnmarg()
c
c-----
c       dip R*4 - dip of fault
c       rake    R*4 - rake angle of motion on fault
c       strike  R*4 - strike of fault
c       isds    I*4 - indicator of couple source description
c                 -1 none given
c                  0 strike, dip, rake
c                  1 moment tensor
c                  2 explosion
c       xmom    R*4 - seismic moment of strike, dip, rake form
c                 1.0 == 1.0 dyne-cm for km,gm,km/s system
c-----
        character*25 name
        integer*4 mnmarg
        integer nsft

        nsft = 0
        nmarg=mnmarg()
        i = 0
   11   i = i + 1
        if(i.gt.nmarg)go to 13
            call mgtarg(i,name)
            if(name(1:2).eq.'-N')then
                i=i+1
                call mgtarg(i,name)
                read(name,'(bn,i10)')nsft
            else if(name(1:2).eq.'-?')then
                call usage()
            else if(name(1:2).eq.'-h')then
                call usage()
            endif
        goto 11
   13   continue
        return
        end
        
        subroutine usage()
        parameter (LER=0, LIN=5, LOT=6)
        write(LER,*)'wvfgrd96 -N nshft ',
     1      ' -h -?'
        write(LER,*)
     1  ' '
        write(LER,*)
     1  ' -N nshift  (default 0 ) number of +- time shifts'
        write(LER,*)
     1  ' -?                   This online help'
        write(LER,*)
     1  ' -h                   This online help'
        stop
        end

        subroutine doshift(mshft,pred,bg,obs,bo,npts,dt)
        implicit none
        integer mshft, npts
        real pred(npts), bg, obs(npts),bo,dt
c-----
c       apply the shift to the data stream - also change 
c       the number of points
c       This would be simpler if we were to use new arrays but to save
c       space we will overwrite arrays
c
c       mshft - shift in samples: =0 no shift
c                     >0 shift observed to left, change
c                       npts and bo
c                     <0 shift green to left, change
c                       npts and bg
c                    |=0 shorten length of time series
c       npts - number of time samples in series - returned as
c               npts - abs(mshft)
c       bg   - SAC header to time stamp of B for first sample for Green
c       bo   - SAC header to time stamp of B for first sample for obs
c       dt   - sample interval
c-----
        integer maxpts, nprm, nmax, mxwvs
        parameter(maxpts=2048, nprm=6, mxwvs=1000, nmax =mxwvs*maxpts)

        integer i, j, n, mpts
        n = npts
        if(mshft.eq.0)then
            return
        else if(mshft.lt.0)then
            mpts = npts + mshft
            bg = bg - mshft*dt
            do i=1,mpts
                j = i - mshft
                if(j.ge.1 .and. j.le.n)then
                    pred(i) = pred(j)
                else
                    pred(i) = pred(npts)
                endif   
            enddo
        else if(mshft.gt.0)then
            mpts = npts - mshft
            bo = bo + mshft*dt
            do i=1,mpts
                j = i + mshft
                if(j.ge.1 .and. j.le.n)then
                    obs(i) = obs(j)
                else
                    obs(i) = obs(npts)
                endif   
            enddo
        endif
        return
        end
