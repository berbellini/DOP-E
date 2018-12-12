        program wvfmtd96
c---------------------------------------------------------------------c
c                                                                     c
c      COMPUTER PROGRAMS IN SEISMOLOGY                                c
c      VOLUME VIII                                                    c
c                                                                     c
c      PROGRAM: WVFMTD96                                              c
c                                                                     c
c      COPYRIGHT 2002 Robert B. Herrmann  and C.  J.  Ammon           c 
c                                                                     c
c      Department of Earth and Atmospheric Sciences                   c
c      Saint Louis University                                         c
c      221 North Grand Boulevard                                      c
c      St. Louis, Missouri 63103                                      c
c      U. S. A.                                                       c
c                                                                     c
c      Department of Geosciences                                      c
c      Penn State University                                          c
c      Deike Bldg                                                     c
c      University Park, PA 16802                                      c
c                                                                     c
c---------------------------------------------------------------------c
c       FIXES
c       21 OCT 2002 - corrections to SAC header for prediction
c       07 JAN 2003 - replaced Numerical Recipes routines
c       16 SEP 2004 - the predicted have the correct 
c           year day hr min sec msec
c       16 FEB 2006 - output the observed file as an .obs
c       13 MAR 2006 - restructure - build in dynamic 
c           weighting and time shift
c           following ideas of Harley Benz
c       04 OCT 2006 - In an effort to both comply with and promulgate 
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
c	05 JAN 2009 - ensure that KNETWK is set in the .obs and .pre files
c       21 MAY 2010 - Carol Bryan noted a big problem with the
c                     deviatoric constraint at line 450
c       03 APR 2011 - the predicted has USER5 set which is the percentage
c                     of set, e.g., 100 * correlation coefficient-squared
c-----
c
c      deviatoric moment tensor inversion
c
c      Time function is fixed, so problem is linear
c
c      Author:  Chuck Ammon, UC Santa Cruz, Saint Louis University
c               parts by George Randall, U South Carolina
c
c      Version: 2.1 September, 1994
c
c      Notation from Langston (1981)
c        x is north, y is east, z is down
c
c-----
c       maxpts  Maximum number of points in time series
c       nprm    Dimension of momenttensor matrix
c       mxwvs   maximum number of observed traces considered
c       nmax    dimension of array
c-----
        integer maxpts, nprm, nmax, mxwvs
        parameter(maxpts=2048, nprm=6, mxwvs=1000, nmax=mxwvs*maxpts)
c-----
c       common blocks for trace and Green functions for 
c           station/component
c-----
        common /gftns/ vss(maxpts),vds(maxpts),vdd(maxpts),vex(maxpts),
     1                response(maxpts)
        real vss, vds, vdd, vex, response
        common /observations/ obs(maxpts)
        real obs
c-----
c       common block for matrix inversion and solution
c-----
        common/mat/a,asave,datav,datasv,covmatrix,mtensor
        real a(nmax,nprm), asave(nmax,nprm)
        real datav(nmax), datasv(nmax)
        real covmatrix(nprm,nprm)
        real mtensor(6)
c-----
c       common blocks for traces and trace attributes
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
c       common for best double couple and for solution vector
c-----
        real pred(nmax)
        real stk0,dip0,rak0,stk1,dip1,rak1,xmw,pclvd

        real btime,  m0
        real total_sum_sq, total_sum_sq_wt
c-----
c       variables for this routine
c-----
        integer ibegin
        integer stdin,stdout,ounit,inunit
        integer ipass
        integer nwaves

        character*256 ofile


c-----
c       initialize everything
c-----
        call initialize(stdin, stdout, ounit, inunit, nwaves,
     1      m0,douse,obswt,
     2      ofile, wvtype, obs_name, gftn_name, wt, mshft)
c-----
c      Output some event and inversion 
c           information (using the first data file)
c-----
        open(unit = ounit, file = ofile)
        write(ounit,*)'Number of traces:',nwaves
        write(ounit,*)'*********** L2 INVERSION **********************'
        call getevent_data(obs_name(1),ounit,'wvfmtd96')
c-----
c       begin the three passes
c       pass 1: do the inversion, determine time shift, set weights
c       pass 2: do the inversion with shifted time series
c       pass 3: final run output obs and pre traces
c-----
        do ipass=1,2

        ibegin = 1
        total_sum_sq = 0
        total_sum_sq_wt = 0
        itotal_pts = 0
c-----
c      Now read in the data, GFs, and set up the matrix equations
c-----
            call getobsgf(nwaves,ibegin,
     1                  btime,ounit,stdout,
     2          total_sum_sq,total_sum_sq_wt,ipass)

c-----
c     END of SETUP
c-----
c     SOLVE THE SYSTEM OF EQUATION
c-----
        call domtinv(ibegin,itotal_pts,ounit,m0)
c-----
c       output waveforms and compute fit to individual traces
c-----
            call dosumup(nwaves,ounit,m0,
     1      total_sum_sq, total_sum_sq_wt, 
     4      ipass)
c
        enddo
        close(ounit)
c-----
c       be sure to close all open files
c-----
        end
c-----
c      END of MAIN PROGRAM
c-----

        subroutine initialize(stdin, stdout, ounit, inunit, nwaves,
     1      m0,douse,obswt,
     2      ofile, wvtype, obs_name, gftn_name, wt, mshft)
        implicit none
        integer maxpts, nprm, nmax, mxwvs
        parameter(maxpts=2048, nprm=6, mxwvs=1000, nmax=mxwvs*maxpts)
        integer stdin, stdout, ounit, inunit, nwaves
        integer wvtype(mxwvs), mshft(mxwvs), douse(mxwvs)
        integer mpts
        real m0, wt(mxwvs), obswt(mxwvs)
        real depmax, depmin
        character*256 obs_name(mxwvs), gftn_name(mxwvs), ofile, infile
        integer i, nerr
        logical lexists

        stdin = 5
        stdout = 6
        ounit = 4
        inunit = 11
c-----
        write(stdout,*)'What is the input file name?'
        read(stdin,*) infile
c-----
        write(stdout,*)'What is the output file name?'
        read(stdin,*) ofile
c-----
        INQUIRE(file=infile, EXIST = lexists)
        if(.not. lexists) then
           write(stdout,*)'Input file not found - quitting.'
           stop
        endif
c-----
c       ASSUME GREEN s FUNCTIONS ARE IN SAME PHYSICAL UNITS AS
c       OBSERVED - but that Greens are for a Moment of 1.0e+20 dyne cm
c-----
c-----
c       m0 is the moment scaling factor for the Green s functions
c       For Computer programs in seismology with earth model in
c       km/sec and gm/cc, m0 == 1.0+20
c-----
        m0 = 1.0e+20

        open(unit = inunit, file = infile)
        i = 1
        nwaves = 0
    5   continue
            read(inunit,*,end=6)wvtype(i),obs_name(i),gftn_name(i),
     1          wt(i)
            nwaves = i
            mshft(i) = 0
            douse(i) = 1
            obswt(i) = 1
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
        return
        end

        subroutine getobsgf(nwaves,ibegin,
     1                  btime,ounit,stdout,
     2          total_sum_sq,total_sum_sq_wt,ipass)
        implicit none
        integer nwaves
        integer maxpts, nprm, nmax, mxwvs
        parameter(maxpts=2048, nprm=6, mxwvs=1000, nmax=mxwvs*maxpts)
c-----
c       common blocks for trace and Green functions for 
c           station/component
c-----
        common /gftns/ vss(maxpts),vds(maxpts),vdd(maxpts),vex(maxpts),
     1                response(maxpts)
        real vss, vds, vdd, vex, response
        common /observations/ obs(maxpts)
        real obs
c-----
c       common blocks for traces and trace attributes
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
c       common block for matrix inversion and solution
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
c        compute the norm of the data
c-----
          etype = 2
          call getnorm(obs,npts(i),etype,sum_sq_amp(i))
c-----
c       NEW 06 APR 2006
c       get variance of window before the P and after the P
c       NOTE this assumes that there is a P time
c       the variance weighting is in the range of 
c            1  < S/N   < inf
c           0.05 < varwt < 1 
c       NOTE THIS IS COMPLETELY ADHOC AND NEETS MORE WORK
c       basically is the waveform is good visually it should 
c           carry more weight
c       use sqrt variance since we really want this to sort of reject
c       S/N on order of 2 - hence the 0.25 term
c       also check for bad numbers
c-----
        if(deltao .gt.0.0)then
            j = 1 + (tpo(i) - bo(i))/deltao
            if(j.gt.1 .and. j.lt. npts(i)/2)then
                call getnorm(obs,j,etype,varnoise)
                call getnorm(obs(j),j+j,etype,varsignal)
                varwt = 1.0 - 1.95/(sqrt(0.25*varsignal/varnoise) + 1)
c-----
c       safety
c-----
                if(varwt.le.0.0)varwt = 0.001
            else
                varwt = 1.0
            endif
        else
            varwt = 1.0
        endif
        

c-----
c        SET UP WEIGHTING
c
c        if the input wt >= 0, 
c           then the weight is simply abs(wt(i)) (no sum_sq_amp
c           factor is included
c-----
        if(douse(i) .eq. 1)then
            obswt(i) = 1.0 * varwt * abs(wt(i))
        else
            obswt(i) = 0.1 * varwt * abs(wt(i))
        endif
c
          write(ounit,*) ' Sum of Square Amplitudes: ',sum_sq_amp(i)
          write(ounit,*) ' Weight:  ',obswt(i)
          write(ounit,*) ' '
c-----
c
c        keep track of the total sum of square amplitudes to provide 
c           "normalized"
c         error measurements
c-----
          total_sum_sq = total_sum_sq + sum_sq_amp(i)
          total_sum_sq_wt = total_sum_sq_wt + (sum_sq_amp(i) * obswt(i))
c-----
c
c        SET-UP THE INVERSION
c
c        compute the response for individual mt elements
c        and store in the system matrix a
c-----
          if(ibegin .gt. nmax) then
            write(stdout,*) 'ERROR: Overflowing the a matrix'
            write(stdout,*) 'ibegin = ',ibegin
            write(stdout,*) 'max allowed =', nmax
            stop
          endif
c       
c-----
c        a copy of the "a matrix" is saved  (before weighting)
c           for later computation of errors and predicted fits
c           a quick alternative is to just read in the seismograms
c           and build the A matrix again later
c
c-----
c       in what follows below the unweighted Green and obs are saved in
c       asave and datasv respectively, the weighted are saved in
c       a and datav respectively
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
c       FORCE DEVIATORIC
c       THIS IS DONE BY APPLYING THE CONSTRAINT THAT
c       m11 + m22 + m33 = 0 , or
c       mtensor(1) + mtensor(2) + mtensor(6) = 0
c
c       If this were not done, we would have a general 
c           moment tensor inversion
c-----
C        do 141 j=1,npts(i)
C           a(j,1) = a(j,1) - a(j,6)
C           a(j,2) = a(j,2) - a(j,6)
C           a(j,6) = 0.0
C 141   continue
c-----
c        error fixed noted by Carol Bryan 20 MAY 2010
c-----
        do 141 j=1,npts(i)
            a(j + ibegin -1,1) = a(j + ibegin -1,1) - a(j + ibegin -1,6)
            a(j + ibegin -1,2) = a(j + ibegin -1,2) - a(j + ibegin -1,6)
            a(j + ibegin -1,6) = 0.0
  141   continue
            
c-----
c        next set up the data vector (save a copy here as well)
c----- 
          call appendvector(datasv,nmax,ibegin,obs,npts(i))
          call appendvector(datav ,nmax,ibegin,obs,npts(i))
c-----
c        apply the weighting
c-----
          if(obswt(i) .ne. 1.0) then
            call apply_weights(a,nmax,ibegin,npts(i),1,nprm,
     1          datav,obswt(i))
          endif
c-----
c        keep track of the number of rows in a matrix
c-----
          ibegin = ibegin + npts(i)
c-----
c        ibegin points to the next unused row of the a matrix
c-----
  125   continue
c-----
c     end loop over waveforms
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

        subroutine domtinv(ibegin,itotal_pts,ounit,m0)
        implicit none
        integer maxpts, nprm, nmax, mxwvs
        parameter(maxpts=2048, nprm=6, mxwvs=1000, nmax=mxwvs*maxpts)
        integer ibegin
c-----
c       common block for matrix inversion and solution
c-----
        common/mat/a,asave,datav,datasv,covmatrix,mtensor
        real a(nmax,nprm), asave(nmax,nprm)
        real datav(nmax), datasv(nmax)
        real covmatrix(nprm,nprm)
        real mtensor(6)
c-----
c       common for best double couple and for solution vector
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
c-----
c     set up the number of rows in the a matrix
c-----
        nr = ibegin -1
c----------
c       Estimate 6 components of the moment tensor
c       using singular value decomposition
c
c       A X = Y
c                T
c       A = U L V 
c              -1 T
c       X = V L  U Y
c-----
        nc = nprm
c-----
c       solve the normal equations using an SVD
c-----
c-----
c       SVD
c           nmax - number of rows in declaration of A
c           nr   - actual number of rows in A for this data set
c           nc   - number of columns declared in A
c           .true. indicates we want the U and V matrices
c           sval - singular values
c           ierr 0 for normal return
c-----
        do 142 j=1,nc
            do 143 i=1,nr
                amat(i,j) = a(i,j)
  143       continue
  142   continue
        call svd(nmax,nr,nc,amat,sval,.true.,u,.true.,v,ierr,tmp)
c
        write(ounit,*)' '     
        write(ounit,*)' Singular Values of A:'     
        do 150 i = 1,nprm
            write(ounit,*) sval(i)
  150   continue
c-----
c       solve the Ax = b system
c       A X = Y
c                T
c       A = U L V 
c              -1 T
c       X = V L  U Y
c-----
        do 160 i=1+ierr,nprm
            tsum = 0.0
            if(sval(i).ne.0.0)then
            do 161 j=1,nr
                tsum = tsum + u(j,i)*datav(j)/sval(i)
  161       continue
            endif
            tmp(i) = tsum
  160   continue
        do 163 i=1+ierr,nprm
            tsum = 0.0
            do 162 j=1,nprm
                tsum = tsum + v(i,j)*tmp(j)
  162       continue
            mtensor(i) = tsum
  163   continue
c-----
c       define the 6 th element of the moment tensor from the
c       deviatoric constraint
c-----
        mtensor(6) = -(mtensor(1) + mtensor(2))

        call mtdec(mtensor,m0,ounit,stk0,dip0,rak0,
     1      stk1,dip1,rak1,xmw,pclvd)
c-----
c   Compute the Covariance Matrix
c       
c       COV(I,J) = SUM V(i,k)V(j,k)/sval(k)^2
c
c-----
        do 170 i=1,nprm
            do 171 j=1,nprm
                tsum = 0.0
                do 172 k=1,nprm
                    if(sval(k).ne.0.0)then
                    tsum=tsum+v(i,k)*v(j,k)/sval(k)**2
                    endif
  172           continue
                covmatrix(i,j) = tsum
  171       continue
  170   continue
        call print_cm(covmatrix,nprm,ounit)
c-----
c      COMPUTE THE PREDICTED WAVEFORMS
c----- 
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
        return
        end

        subroutine dosumup(nwaves,ounit,m0,
     1      total_sum_sq, total_sum_sq_wt, 
     4      ipass)
        implicit none
        integer maxpts, nprm, nmax, mxwvs
        parameter(maxpts=2048, nprm=6, mxwvs=1000, nmax=mxwvs*maxpts)

c-----
c       common block for matrix inversion and solution
c-----
        common/mat/a,asave,datav,datasv,covmatrix,mtensor
        real a(nmax,nprm), asave(nmax,nprm)
        real datav(nmax), datasv(nmax)
        real covmatrix(nprm,nprm)
        real mtensor(6)
c-----
c       common blocks for traces and trace attributes
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
c       common for best double couple and for solution vector
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
c       function prototype
c-----
        real vdot
        integer lgstr
c-----
c       get the fit between observed and predicted
c-----
        xx = vdot(datasv, datasv, itotal_pts)
        xy = vdot(datasv, pred  , itotal_pts)
        yy = vdot(pred  , pred  , itotal_pts)
c-----
c       write the observed .obs and predicted .pre waveforms
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
c           write out the prediction as a sac file
c-----
            if(ipass.eq.2 )then
                call preobssv(i,ibegin)
            endif
c-----
c        COMPUTE THE ERRORS (scale factor = 1.0, l2 norm)
c
c    WARNING!!!!!  There are several ways to compute 
c    the error or misfit between the observed and synthetic 
c    seismograms.  This part of the program has been changes 
c    from the original code of Ammon, so please beware.  
c    I have modified the definition used to compute the misfit.  
c    The formulation that I use is misfit = ( data - syn )**2 / data**2,
c     so if data = syn, misfit = 0; if data = -syn, misfit = 4; 
c    if syn = 0, misfit = 1.  This is similar to the error 
c    function of Pasyanos et al. (1996).  This is the error 
c    that is passed through to the subroutine the 
c    prints the ASCII solution.  
c-----

            err = 0.0
            call docross(datasv(ibegin),pred(ibegin),npts(i),
     1          err,xmis,mshft(i),ipass)
            misfit(i) = xmis
C           if(misfit(i).gt.0.9)douse(i) = 0
            if(misfit(i).gt.1.5)douse(i) = 0
            if(obswt(i) .ne. 0.0) then
                    total_err = total_err + err
                    total_wtd = total_wtd + err * obswt(i)
                    n = n + npts(i)
            endif
c-----
c        compute a scale factor for the seismogram fit
c          not used, just a statistic
c-----
            call scale_fit(datasv(ibegin),pred(ibegin),npts(i),factor)
c-----
c        the value of sum_sq_amp(i) != 0, 
c        this is checked when first computed
c-----

    2   format(a8,1x,a8,1x,f6.3,e10.3,f6.3,e10.3,e10.3,f10.2)
            write(ounit,2)kstnm(i), kcmpnm(i),obswt(i),
     1          err,err/sum_sq_amp(i),factor,misfit(i)
     2          ,mshft(i)*dt(i)
c 
            ibegin = ibegin + npts(i)
c
            fr_errors(i) = err/sum_sq_amp(i)
c
  350   continue
c      end loop over waveforms
c
c
c**********************************************************************
c      Graphical error summary
c
C       call bar_chart(nwaves,obs_name,fr_errors,12,50,ounit)
c
c***********************************************************************
c
c      OUTPUT the error for all waveforms
c
c      output the sum of the square errors divided by the 
c        number of traces
c        minus the number of unknowns (free parameters)
c
        tmpsig = 1.0 / float(n-6)
c
        write(ounit,*) ' '
        write(ounit,*)
     1 'Unweighted SSE,Weighted SSE,Unweighted RMS,Weighted RMS'
        write(ounit,*) total_err*tmpsig, total_wtd*tmpsig,
     1      sqrt(total_err*tmpsig), sqrt(total_wtd*tmpsig)
c
c      Output sum of square errors divided by the sum of the 
c        square amplitudes.
c
        write(ounit,*)
     1 'Total Error divided by sum of all seismogram amplitudes'
        write(ounit,*)
     1 'Unweighted, weighted'
        write(ounit,*) total_err/total_sum_sq, 
     1     total_wtd/total_sum_sq_wt
c-----
c       Now output the Moment Tensor Solution Again but with
c       Error estimates
c-----
        call out_mte(mtensor,total_err*tmpsig,covmatrix,nprm,m0,ounit)
c-----
c       output L2 summary file
c       depthg  - depthg of the Greens functions
c       stk - strike of best double couple
c       dip - dip    of best double couple
c       rake    - rake   of best double couple
c       xmw - moment magnitude of best double couple
c       1.0-total_err/total_sum_sq - raw data 0 zero is best
c       sqrt(total_err*tmpsig)  - std err
c       xy/sqrt(xx*yy) - cosine of angle between observed and predicted
c       1.0-total_wtd/total_sum_sq_wt - weighted data zero is best
c       sqrt(total_wtd*tmpsig)  - weighted std err
c       pcvld - percent CLVD - 0 means pure double couple
c-----
        if(ipass.eq.2)then
        open(3,file='fmdfit.dat',access='sequential',
     1      form='formatted',status='unknown')
        rewind 3
        
        write(3,3)depthg(1), stk0,dip0,rak0,xmw, 
     1     1.0-total_err/total_sum_sq,
     1      sqrt(total_err*tmpsig), 
     2      1.0-total_wtd/total_sum_sq_wt, 
     3      xy/sqrt(xx*yy), 
     4      sqrt(total_wtd*tmpsig) , pclvd
    3   format('WVFMTD96 ',f6.1,3f6.0,f7.2,f10.3,e10.3,
     1     2f10.3,e10.3,f6.1)
        call outmsg(nwaves)
        close (3)
        endif
        return
        end

        subroutine docross(o,p,npts,sumerr,xmis,mshft,ipass)
        implicit none
        integer maxpts, nprm, nmax, mxwvs
        parameter(maxpts=2048, nprm=6, mxwvs=1000, nmax=mxwvs*maxpts)
c-----
c       create an envelope
c       cross-correlate the envelope to get the shift
c       use the shift to generate the misfit post-shift
c       ipass = 1 - compute the shift
c             = 2 - do not compute shift just get sum err
c-----
        integer npts
        real o(npts), p(npts)
        real xmis, sumerr
        integer mshft, ipass

        real oe(maxpts), pe(maxpts)
        integer i,is,j
        real dot, dotsv
        real xnum, xden
        real xx, xy, yy

c-----
c       calculate the envelope
c-----
        call do_envelope(o,oe,npts)
        call do_envelope(p,pe,npts)
c-----
c       determine the best time shift on the basin of a dot product
c-----
        dotsv = -1.0
        do is=-20,20
            xx = 0.0
            xy = 0.0
            yy = 0.0
            do i=1,npts
                j = i - is
                if(j.ge.1 .and. j.le.npts)then
                    xx = xx + (o(i) * o(i))
                    xy = xy + (o(i) * p(j))
                    yy = yy + (p(j) * p(j))
                endif
            enddo
            dot = xy/sqrt(xx*yy)
            if(dot.gt.dotsv)then
                dotsv = dot
                mshft = is
            endif
        enddo
c-----
c       we will determine the goodness of fit for a +- 20 point shift
c       For speed I will do a dot product - then I will compute the
c       misfit = SUM(o - p)**2 / SUM (o)**2
c-----
c-----
c       calculate the misfit by applying the shift
c-----
        xnum = 0.0
        xden = 0.0
        if(ipass.eq.2)mshft=0
        do i=1,npts
            j = i - mshft
            if(j.ge.1 .and. j.le.npts)then
                xnum = xnum + (o(i) - p(j))**2
                xden = xden + (o(i)       )**2
            endif
        enddo
        xmis = xnum/ xden
        sumerr = xnum
        
        
        return
        end

        subroutine shiftall(mshft,npts,bg,bo,dt)
        integer mshft, npts
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
        parameter(maxpts=2048, nprm=6, mxwvs=1000, nmax=mxwvs*maxpts)

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
c       safety
c-----
        if(npts.lt.0)npts = 0
        return
        end

c------
c       Routines to obtain an envelope for cross-correlation
c       to determine time shifts
c-----
c-----
c       implementation of an envelope function using 
c       IIR filters
C       Reference: http://www.biochem.oulu.fi/~oniemita/dsp/hilbert/
c       July 2003, Olli Niemitalo, o@iki.fi
c       ---
c       create an all pass filter and one roughly 
c           90 degrees out of phase
c       note - the IIR is causal so there will be a slight delay
c-----
        subroutine do_envelope(x,e,n)
        integer n
        real x(n), e(n)
        integer maxpts
        parameter (maxpts=2048)
        real h1(maxpts), h2(maxpts), t(maxpts)
        real a1(4), a2(4)
        data a1/0.6923878,0.9360654322959,
     1      0.9882295226860,0.9987488452737/
        data a2/0.4021921162426,0.8561710882420,
     1      0.9722909545651,0.9952884791278/
c-----
c       state 1
c-----
        call h_copy(x,h1,n)
        do i=1,4
            call h_sect(h1,t,n,a1(i))
            call h_copy(t,h1,n)
        enddo
        call h_delay(h1,t,n)
        call h_advance(t,h1,n,2)
c-----
c       state 2
c-----
        call h_copy(x,h2,n)
        do i=1,4
            call h_sect(h2,t,n,a2(i))
            call h_copy(t,h2,n)
        enddo
        call h_advance(t,h2,n,2)
c-----
c       compute the envelope
c-----
        call do_env(h1,h2,e,N)
        return
        end

        subroutine do_env(x,y,z,n)
        integer n
        real x(n), y(n), z(n)
        integer i
        do i=1,n
            z(i) = sqrt(x(i)*x(i) + y(i)*y(i))
        enddo
        return
        end

        subroutine h_copy(x,y,n)
        integer n
        real x(n), y(n)
        integer i
        do i=1,n
            y(i) = x(i)
        enddo
        return
        end

        subroutine h_sect(x,y,n,a)
        integer n
        real x(n), y(n)
        real a
        integer i
        real a2, xm1, xm2, ym1, ym2
        xm1 = x(1)
        xm2 = xm1
        ym1 = x(1)
        ym2 = ym1
        a2 = a * a
        do i=1,n
            y(i) = a2*(x(i) + ym2) - xm2
            xm2 = xm1
            xm1 = x(i)
            ym2 = ym1
            ym1 = y(i)
        enddo
        return
        end

        subroutine h_delay(x,y,n)
        integer n
        real x(n), y(n)
        integer i
        y(1) = x(1)
        do i=2,n
            y(i) = x(i-1)
        enddo
        return
        end

        subroutine h_advance(x,y,n,nad)
        integer n, nad
        real x(n), y(n)
        integer i
        do i=1,n-nad
            if(i.le.n-nad)then
                y(i) = x(i+nad)
            else
                y(i) = x(n)
            endif
        enddo
        return
        end

        subroutine h_reverse(x,y,n)
        integer n
        real x(n), y(n)
        integer i
        do i=1,n
            y(i) = x(n+1-i)
        enddo
        return
        end

        subroutine preobssv(i,ibegin)
        implicit none
        integer i

        integer maxpts, nprm, nmax, mxwvs
        parameter(maxpts=2048, nprm=6, mxwvs=1000, nmax=mxwvs*maxpts)
c-----
c       common block for matrix inversion and solution
c-----
        common/mat/a,asave,datav,datasv,covmatrix,mtensor
        real a(nmax,nprm), asave(nmax,nprm)
        real datav(nmax), datasv(nmax)
        real covmatrix(nprm,nprm)
        real mtensor(6)
c-----
c       common blocks for traces and trace attributes
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
c       common for best double couple and for solution vector
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
c       function prototype
c-----
        integer lgstr
        real vdot
        double  precision  rxx, rxy, ryy, xxyymx, pctfit
c-----
c       save the predicted abd observed traces
c-----
c           set up predicted trace for this observation
c           but design it so that it permits an overlay of
c           observed and predicted
c
c           convert from dist(km) to degrees using conversion
c           constant build within GSAC
c-----
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
            call setfhv('EVDP    ',depthg (i),nerr)
            call setfhv('DIST    ',disto (i),nerr)
            call setfhv('GCARC   ',disto (i)/111.1946,nerr)
            call setfhv('AZ      ',az    (i),nerr)
            call setfhv('BAZ     ',-12345.  ,nerr)
            call setfhv('USER1   ',user1 (i),nerr)
            call setfhv('USER2   ',user2 (i),nerr)
            call setfhv('O       ',og    (i)-og(i),nerr)
            call setfhv('USER5   ',sngl(pctfit),nerr)
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
c       if we used the trace, we will output the observed and predicted
c       we also take time to ensure that the P and S picks are
c       preserved for comparison
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
c       output the observed
c-----
            call setfhv('A       ',tpo(i),nerr)
            call setfhv('B       ',bo(i),nerr)
            call setfhv('O       ',oo(i),nerr)
            call setfhv('USER5   ', -12345.,nerr)
            call setfhv('USER9   ', -12345.,nerr)
            wvfile = ' '
            slen = lgstr(obs_name(i))
            wvfile(1:slen+4) = obs_name(i)(1:slen)//'.obs'
            if(douse(i).eq.1)then
            call wsac1(wvfile,datasv(ibegin),npts(i),bo(i),dt(i),nerr)
            endif
        return
        end

        subroutine outmsg(nwaves)
        implicit none
        integer nwaves 

        integer maxpts, nprm, nmax, mxwvs
        parameter(maxpts=2048, nprm=6, mxwvs=1000, nmax=mxwvs*maxpts)
c-----
c       common blocks for traces and trace attributes
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
c       common for best double couple and for solution vector
c-----
        common/rbstdc/stk0,dip0,rak0,stk1,dip1,rak1,xmw,pclvd,pred
        real pred(nmax)
        real stk0,dip0,rak0,stk1,dip1,rak1,xmw,pclvd
c-----
c       internal variables
c-----
        integer kyear, kdoy, kmonth, kday, khour, kminute, ksecond
        integer kdate
        real rsecond
        real relat, relon, rm0
        character orgstr*18
        character kstr*32
        real*8 timeo
        character clat*1, clon*1
        real evlat, evlon
        real scalar_mom
        real evms, evmb
        integer nchar, proj
        integer i, nwave
        real maxgap
c-----
c       function prototype
c-----
        integer lgstr
c-----
c       We follow Harley Benz s email format which gives 
c       a mechanism plot 
c       and also tabulates the misfit for the stations used
c-----
c       get origin time from first observed trace
c       we first get the epoch time of the reference time
c-----
c-----
c       define the origin time time as a string for use with the
c       the location summary file YYYYMMDDHHMMSS.msg
c       we only do this for the first trace
c       get month and day, the get reference time, then get origin
c       time, then convert back to human
c-----
        call mnthdy(nzyear(1),nzjday(1),kmonth,kday)
        call htoe(nzyear(1),kmonth, kday, nzhour(1),
     1      nzmin(1), float(nzsec(1))+0.001* nzmsec(1),timeo)
        timeo = timeo + oo(1)
        call etoh(timeo,kdate,kstr,kdoy,kyear,kmonth,kday,
     1      khour,kminute,rsecond)
        ksecond = rsecond
        write(orgstr,1)kyear,kmonth,kday,khour,kminute,ksecond
    1   format(i4.4,i2.2,i2.2,i2.2,i2.2,i2.2,'.msg')
c-----
c       determine the number of traces used as well as the
c       azimuthal gap
c-----
        nwave = 0
        do i=1,nwaves
            if(douse(i).eq.1)nwave = nwave + 1
        enddo
        call azgap(az,nwaves,maxgap)
c-----
c       now open the message file
c-----
        open(2,file=orgstr,access='sequential',form='formatted',
     1      status='unknown')
        rewind 2
            if( evla .lt. 0.0 ) then
                evlat = -evla
                clat = 'S'
            else
                evlat =  evla
                clat = 'N'
            endif
            if( evlo .lt. 0.0 ) then
                evlon = -evlo
                clon = 'W'
            else
                evlon =  evlo
                clon = 'E'
            endif
            
        write(2,102)xmw,' '
        write(2,*)
        write(2,103)kyear,kmonth,kday,khour,kminute,ksecond,
     1                  ' '
        write(2,104)evlat,clat,evlon,clon,int(depthg(1))
c-----
c       DEBUG
c----
        evmb= 0
        evms = 0
c-----
c       end debug
c-----  
        write(2,106)xmw,evmb,evms

        scalar_mom = 10.0** (1.5*xmw + 16.10)
        write(2,107)scalar_mom
        write(2,*)' '
        write(2,105)nwaves,nint(maxgap)
        write(2,115)pclvd
        write(2,116)nint(user1(1)+0.49), nint(user2(1)+0.49) 
c
c Write out the beachball in an ASCII format
c
        write(2,109)
        write(2,110) nint(stk0),nint(dip0),nint(rak0)
        write(2,111) nint(stk1),nint(dip1),nint(rak1)
        write(2,*)'      '
        nchar = 40
        proj = 1
        call rtextbeach(stk0,dip0,rak0,nchar,proj,2)
        write(2,112)
        write(2,113) (kstnm(i),int(disto(i)/111.1946),int(az(i)),
     1                    misfit(i),i=1,nwaves)
        write(2,*)'Misfit > 1 heavily downweighted'

        close(2)
c-----
c       format statements
c-----
  102   format('Mw=',f3.1,1x,a)
  103   format(i4.4,'/',i2.2,'/',i2.2,1x,i2.2,':',i2.2,':',i2.2,4x,
     1        'EVID: ',a)
  104   format('LAT: ',f5.2,a1,', LON: ',f6.2,a1,
     1         ', DEPTH: ',i3)
  106   format('Mw: ',f3.1,', mb: ',f3.1,', Ms: ',f3.1)
  107   format('Mo=',e11.3,' (dyne-cm)')
  105   format('No. of Stations: ',i3,', GAP: ',i3)
  115   format('Error (clvd/dc)*100= ',f6.2)
  116   format('Bandpass: ',i3, ' - ', i3,' sec')
  109   format('Nodal plane parameters',/,
     1      '      strike',5x,'dip',4x,'rake')
  110   format('NP1: ',3x,3(i4,4x))
  111   format('NP2: ',3x,3(i4,4x))
  112   format(/,'Stat Dist Az Misfit   Stat Dist Az Misfit')
  113   format(a4,1x,i3,1x,i3,1x,f5.2,4x,a4,1x,i3,1x,i3,1x,f5.2)


        return 
        end

        subroutine rtextbeach(strike,dip,rake,nchar,proj,lun)
c-----
c       Text based focal mechanism plot with account for
c       terminal/printer aspect ratio
c
c       this also plots the P and T axes
c
c       This is based on Chuck Ammon s lrpmech which is more
c       general in that it will ploit the entire moment tensor.
c       However that code does not plot the P-T axes
c
c       strike  R   - fault strike
c       dip R   - fault dip
c       rake    R   - fault rake
c       nchar   I   - diameter in horizontal direction
c       proj    I   - 0 stereographic (Wulff)
c                 1 equal area    (schmidt
c       lun I   - file logical unit number for output
c-----
        real strike,dip,rake
        integer nchar,proj,lun
        real pi, pi2, ar
        parameter (pi =3.1415926535898)
        parameter (pi2=1.5707963267998)
        parameter (ar=18./33.)
        integer NNC
        parameter (NNC=101)
        integer KK(NNC,NNC)
        integer lx, ly, lc
        character ostr*101
c-----
c       rbh
c-----
        real a11, a12, a13
        real a21, a22, a23
        real a31, a32, a33
        real p1, p2, p3
        real t1, t2, t3
        real xc, yc
        real x, y, z
        integer ifm
        real rad
        real xn
        integer ncx, ncy

        real degrad
        real sr, cr, sd, cd, ss, cs
        real tt, tp, pt, pp
        real tnorm, pnorm
        integer ita, jta, ipa, jpa

        nc = nchar
        nc = ( 2 * float(nc/2) )
        nr = nint(float(nchar)*ar + 0.5)
        nr = ( 2 * (nr/2) )

        x0 = (nc+1)/2
        y0 = (nr+1)/2
        r0 = 0.99*x0

        degrad = pi/180.0
        sr = sin(degrad*rake)
        cr = cos(degrad*rake)
        sd = sin(degrad*dip )
        cd = cos(degrad*dip )
        ss = sin(degrad*strike )
        cs = cos(degrad*strike )

        
        a11 = cr*cs + sr*cd*ss
        a12 = cr*ss - sr*cd*cs
        a13 = - sr*sd
        a21 = - ss*sd
        a22 =   cs*sd
        a23 =  -cd
        a31 = cs*sr - cd*cr*ss
        a32 = ss*sr + cd*cr*cs 
        a33 = cr*sd

        t1 = 0.707*(a11 + a21)
        t2 = 0.707*(a12 + a22)
        t3 = 0.707*(a13 + a23)
        tnorm = sqrt(t1*t1 + t2*t2 + t3*t3)
        t1 = t1 / tnorm
        t2 = t2 / tnorm
        t3 = t3 / tnorm
        if(t3.lt.0.0)then
            t1 = - t1
            t2 = - t2
            t3 = - t3
        endif
        

        p1 = 0.707*(a11 - a21)
        p2 = 0.707*(a12 - a22)
        p3 = 0.707*(a13 - a23)
        pnorm = sqrt(p1*p1 + p2*p2 + p3*p3)
        p1 = p1 / pnorm
        p2 = p2 / pnorm
        p3 = p3 / pnorm
        if(p3.lt.0.0)then
            p1 = - p1
            p2 = - p2
            p3 = - p3
        endif


c-----
c       get trend and plunge of the T axis
c-----
        tt = atan2(t2,t1)
        tp = asin(t3)
c-----
c       get trend and plunge of the P axis
c-----
        pt = atan2(p2,p1)
        pp = asin(p3)
c-----
c       now compute the location of the P and T axes in the grid
c       we actually know where it intersects the focal spere from the
c       t11, t12 coordinates - so just convert to screen coordinates
c-----
c-----
c       now draw initialize the matrix, determine the actual size of
c       the region
c       Notation:    0 blank
c               +1 positive (#)
c               -1 negative (-)
c               +2 T-axis   (T)
c               -2 P-axis   (P)
c-----

        do lx=1,NNC
            do ly=1,NNC
                KK(lx,ly) = 0
            enddo
        enddo
        do i=1,nr
            yy = (y0 - float(i))/ar
            do j=1,nc
                xx = float(j) - x0
                r = sqrt(xx*xx + yy*yy) / r0
                if(r .le. 1.0)then
                    phi = atan2(xx,yy)
                    if(proj.eq.1)then
                        theta = 2*asin(r/1.414213562)
                    else
                        theta = 2*atan(r) 
                    endif
        
c-----
c       now fill in the focal sphere using the grid to define the
c       trend and plunge
c-----
c-----
c       now fill the focal sphere with the double couple first motion
c       value to avoid inverse mapping we will just do an overkill of
c       very fine sampling of the trend and plunge angles
c-----
c-----
c       define the vector that intersects the focal sphere
c-----
                x = sin(theta)*cos(phi)
                y = sin(theta)*sin(phi)
                z = cos(theta)
c-----
c       re-define this vector in terms of the focal mechanism
c       coordinte system
c-----
                xc = a11*x + a12*y + a13*z
                yc = a21*x + a22*y + a23*z
c-----
c       get the first motion
c-----
                if(xc*yc .gt.0.0)then
                    ifm = +1
                else
                    ifm = -1
                endif
c-----
c       project the (x,y,z) point on the unit sphere 
c       according to the transformation - now xc, yc represents
c       coordinates of the transformed point 
c       rad is the radius in the plane
c-----
                KK(i,j) = ifm
                endif
            enddo
        enddo
c-----
c       to do the P and T axes we must convert the trend and plunge to a
c       point on the new coordinates - we first use the projection to
c       map the axis onto the display plane - do the mapping, then
c       check the bounds
c-----
c       T - axis
c-----
        call plotpt(KK,NNC,nr,nc,tt,tp,'T',proj,x0,y0,r0,ar)
        call plotpt(KK,NNC,nr,nc,pt,pp,'P',proj,x0,y0,r0,ar)
                
c-----
c-----
c       output the mechanism
c-----
        do ly=1,nr
            do lx=1,nc
                lc = KK(ly,lx)
                if(lc.eq.-2)then
                    ostr(lx:lx) = 'P'
                else if(lc.eq.-1)then
                    ostr(lx:lx) = '-'
                else if(lc.eq. 0)then
                    ostr(lx:lx) = ' '
                else if(lc.eq.+1)then
                    ostr(lx:lx) = '#'
                else if(lc.eq.+2)then
                    ostr(lx:lx) = 'T'
                else
                    ostr(lx:lx) = '+'
                endif
                    
            enddo
            write(lun,1)ostr(1:nc)
        enddo
    1   format(1x,a)
c-----
c       end rbh
c-----
        return
        end
        subroutine plotpt(KK,NNC,nr,nc,t,p,caxis,proj,x0,y0,r0,ar)
        implicit none
        integer NNC
        integer KK(NNC,NNC)
        integer nr, nc, proj
        real t, p, x0, y0, r0, ar
        character caxis*(*)

        real pi2
        parameter (pi2=1.5707963267998)
        real r, xx, yy
        integer i, j, lx, ly
        integer ifm
        
        if(proj.eq.1)then
            r = 1.414213562*sin((p-pi2)/2.)
        else
            r = tan((p-pi2)/2.)
        endif
        xx = - r0*r*sin(t)
        yy = - r0*r*cos(t)
        i = ( y0 - ar*yy + 0.49)
        j = ( x0 + xx + 0.49)
        if(caxis(1:1).eq.'T')then
            ifm =  2
        else if(caxis(1:1).eq.'P')then
            ifm = -2
        else
            ifm =  0
        endif
c-----
c       safety checks
c-----
        if(j.gt.nc)j=nc
        if(j.lt.1 )j=1
        if(i.gt.nr)i=nr
        if(i.lt.1 )i=1
        do ly=i-1,i+1
           do lx=j-1,j+1
              if(ly.ge.1 .and. ly.le.nr)then
                 if(lx.ge.1 .and. lx.le.nc)then
                    KK(ly,lx) = 0
                    if(lx.eq.j .and. ly .eq.i)then
                       KK(ly,lx) = ifm
                    endif
                 endif
              endif
           enddo
        enddo
        return
        end

        subroutine azgap(az,npts,maxgap)
        parameter (mxwvs=1000)
        real az(mxwvs), azdif(mxwvs), aztmp(mxwvs), maxgap, vmin
        integer npts
c
c Subroutine to compute the maximum azimuth gap between
c   stations
c
c Copy azimuths to a temporary array
c
        do 4 i = 1, npts
            aztmp(i) = az(i)
4       continue
c
c Sort station in increasing order of azimuth
c
        do 5  j = 1, npts
        vmin = aztmp(j)
        do 10 i = j+1, npts
           if(aztmp(i).lt.vmin) then
                vmin = aztmp(i)
                aztmp(i) = aztmp(j)
                aztmp(j) = vmin
           endif
10      continue
5       continue
c
c Compute the difference between the stations
c
        do 30 i = 1, npts-1
           azdif(i) = aztmp(i+1) - aztmp(i)
30      continue
        azdif(npts) = 360 - aztmp(npts) + aztmp(1)
c
c Find the maximum azimuthal gap
c
        maxgap = azdif(1)
        do 40 i = 2, npts
            if( azdif(i) .gt. maxgap) maxgap = azdif(i)
40      continue
        return
        end
