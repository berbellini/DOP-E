C-----
C       Changes
C       22 SEP 2002 - Error in getmtresp
C           do 10 and do 11 loops have vdd(i) vex(i) instead of
C               vdd(j), vex(j)
C           also force response = 0 for SH for izz to avoid dregs
c     27 NOV 2002   relaxed check on equal DT s by using 0.01 relative
C           error e.g., replaced dt.ne.dt2 by (abs(dt-dt2).gt.0.01*dt)
C           
C-----
        subroutine scale_fit(obs,pred,npts,scl)
c-----
c     Given the model Y = a X, the least sqaures estimate of
c     a is     T     T
c         a = Y Y / Y X
c    
c     obs R   - array of observations
c     pred    R   - array of predictions
c     npts    I   - number of obseervations
c     scl R   - scale factor
c-----
        implicit none
        integer npts
        real obs(npts), pred(npts)
        real scl
        real numerator, denominator
        integer i
        numerator = 0.0
        denominator = 0.0
              
        do 10 i = 1, npts
            numerator = numerator + obs(i)*pred(i)
            denominator = denominator + pred(i)*pred(i)
   10   continue
        scl = numerator / denominator
       
        return
        end
       

       
        subroutine sum_err(obs,pred,s,npts,etype,err)
c-----
c     Obtain the Norm divided by number of observations
c     obs R   - array of observations
c     pred    R   - array of predictions
c     s   R   - scale factor for model waveforms
c     npts    I   - number of data points
c     etype   I   - norm type 1: L1 norm, 2: L2 norm
c     err R   - output errors
c-----
        implicit none
        integer npts, etype
        real obs(npts), pred(npts)
        real s, res
        real err

        integer l1, l2
        integer i
       
        l2 = 2
        l1 = 1
       
        err = 0
       
c-----
c     Simple L2 norm
c-----
        if(etype .eq. l2)then
            do 10 i = 1, npts
                res = obs(i) - s*pred(i)
                err = err + res*res
   10       continue
            err = err/float(npts)
        endif
c-----
c     Simple L1 norm
c-----
        if(etype .eq. l1)then
            do 20 i = 1, npts
                res = obs(i) - s*pred(i)
                err = err + abs(res)
   20       continue
            err = err/float(npts)
        endif
        return
        end

         
        subroutine getnorm(x,n,etype,norm)
c------
c     determine the norm of the vector X
c     X   R   - vector
c     N   I   - number of points in vector
c     ETYPE   I   - 1 return L1 norm
c               2 return L2 norm
c-----
        implicit none
        integer n, etype
        real x(n), norm

        integer l1,l2,i
      
        l1 = 1
        l2 = 2
        norm = 0.0
      
        if(etype .eq. l1) then
            do 1 i = 1, n
                norm = norm + abs(x(i))
    1       continue
        else if(etype .eq. l2)then
            do 2 i = 1, n
                norm = norm + x(i)*x(i)
    2       continue
        endif
       
        return
        end


        subroutine getwaveforms(gftn_name,obs_name,
     1      iwvtype,npts,dt,btime,ounit,kstnm,kcmpnm,knetwk,
     2      disto,distg,depthg,az,user1,user2,
     3      tpg,tsg,bg,og,
     3      tpo,tso,bo,oo,
     3      deltao,
     3      nzyear,nzjday, nzhour, nzmin, nzsec, nzmsec)            
c-----
c
c routine that reads in and stores the green functions and observations
c         into the common blocks
c
c-----
c     Note the nzyear,nzjday, nzhour, nzmin, nzsec, nzmsec
c     refers to the origin time of the observed
c-----
        implicit none
        character*256 obs_name, gftn_name
        integer iwvtype, npts, ounit
        real dt, btime, disto, distg, depthg, az, user1, user2
        real tpo, tso, bo, oo
        real tpg, tsg, bg, og
        real deltao
        integer nzyear,nzjday, nzhour, nzmin, nzsec, nzmsec
        character*8 kstnm, kcmpnm, knetwk

        
        integer maxpts, nprm, nmax, mxwvs 
        parameter(maxpts=2048, nprm=6, mxwvs=1000, nmax=mxwvs*maxpts)
        real deg_to_rad, rad_to_deg
        parameter(deg_to_rad = 0.017453292, rad_to_deg = 57.29577951)

        integer stdin,stdout
        integer iPz,iPr,iSH
        integer lgstr
        common/gftns/vss(maxpts),vds(maxpts),vdd(maxpts),vex(maxpts),
     &                response(maxpts)
        real vss, vds, vdd, vex, response

        common /observations/ obs(maxpts)
        real obs

        integer nptsdat, nptsgrn
        character*256 infile
        integer nerr, sleng, sleno
        real dt2
        real otime, second
        character*32 ascfmt
        real*8 epoch
        integer month, day
        integer*4 date
        character str*32
      
        iPz = 1
        iPr = 2
        iSH = 3
       
        stdin = 5
        stdout = 6
c-----
c
c     read in the fundamental fault responses
c     created using f96tosac -G in Computer Programs in Seismology
c
c-----       
        sleng = lgstr(gftn_name)
        sleno = lgstr(obs_name)
        infile = ' '
        if(iwvtype .eq. iPz) then
            infile(1:sleng+4) = gftn_name(1:sleng)//'.ZSS'     
            call rsac1(infile,vss,nptsgrn,btime,dt,maxpts,nerr)
            infile(1:sleng+4) = gftn_name(1:sleng)//'.ZDS'     
            call rsac1(infile,vds,nptsgrn,btime,dt,maxpts,nerr)
            infile(1:sleng+4) = gftn_name(1:sleng)//'.ZDD'     
            call rsac1(infile,vdd,nptsgrn,btime,dt,maxpts,nerr)
            infile(1:sleng+4) = gftn_name(1:sleng)//'.ZEX'     
            call rsac1(infile,vex,nptsgrn,btime,dt,maxpts,nerr)
c-----
c         convert Greens from cm to meters
c             or cm/s to m/sec - displacement / velocity and
c             acceleration comes from hpulse96 and deconvolution
c             of real data and anything done in the
c             script DOSTA
c-----
            call vmul(vss,nptsgrn,0.01)
            call vmul(vds,nptsgrn,0.01)
            call vmul(vdd,nptsgrn,0.01)
            call vmul(vex,nptsgrn,0.01)
            call getfhv('O       ',og    , nerr)
            call getfhv('A       ',tpg   , nerr)
            call getfhv('T0      ',tsg   , nerr)
        else if(iwvtype .eq. iPr) then
            infile(1:sleng+4) = gftn_name(1:sleng)//'.RSS'     
            call rsac1(infile,vss,nptsgrn,btime,dt,maxpts,nerr)
            infile(1:sleng+4) = gftn_name(1:sleng)//'.RDS'     
            call rsac1(infile,vds,nptsgrn,btime,dt,maxpts,nerr)
            infile(1:sleng+4) = gftn_name(1:sleng)//'.RDD'     
            call rsac1(infile,vdd,nptsgrn,btime,dt,maxpts,nerr)
            infile(1:sleng+4) = gftn_name(1:sleng)//'.REX'     
            call rsac1(infile,vex,nptsgrn,btime,dt,maxpts,nerr)
c-----
c         convert Greens from cm to meters
c-----
            call vmul(vss,nptsgrn,0.01)
            call vmul(vds,nptsgrn,0.01)
            call vmul(vdd,nptsgrn,0.01)
            call vmul(vex,nptsgrn,0.01)
            call getfhv('O       ',og    , nerr)
            call getfhv('A       ',tpg   , nerr)
            call getfhv('T0      ',tsg   , nerr)
        else if(iwvtype .eq. iSH) then
            infile(1:sleng+4) = gftn_name(1:sleng)//'.TSS'     
            call rsac1(infile,vss,nptsgrn,btime,dt,maxpts,nerr)
            infile(1:sleng+4) = gftn_name(1:sleng)//'.TDS'     
            call rsac1(infile,vds,nptsgrn,btime,dt,maxpts,nerr)
c-----
c         convert Greens from cm to meters
c-----
            call vmul(vss,nptsgrn,0.01)
            call vmul(vds,nptsgrn,0.01)
            call getfhv('O       ',og    , nerr)
            call getfhv('A       ',tpg   , nerr)
            call getfhv('T1      ',tsg   , nerr)
        endif
 
        call getfhv('EVDP    ',depthg,nerr)
        call getfhv('DIST    ',distg, nerr)
        call getfhv('B       ',bg    , nerr)
 
c-----
c     read in the observed
c-----             
        sleno = lgstr(obs_name)
        infile = ' '
        infile(1:sleno) = obs_name(1:sleno)     
        call rsac1(infile,obs,nptsdat,btime,dt2,maxpts,nerr) 
      
        call getnhv('NZYEAR',nzyear,nerr)
        call getnhv('NZJDAY',nzjday,nerr)
        call getnhv('NZHOUR',nzhour,nerr)
        call getnhv('NZMIN' ,nzmin ,nerr)
        call getnhv('NZSEC' ,nzsec ,nerr)
        call getnhv('NZMSEC',nzmsec,nerr)
        call getfhv('O'     ,otime, nerr)
        call getfhv('A'     ,tpo, nerr)
        call getfhv('T1'    ,tso, nerr)
        call getfhv('B'     ,bo, nerr)
        call getfhv('O'     ,oo, nerr)
        call getfhv('DELTA' ,deltao, nerr)

c-----
c     at this point, the date refers to the reference time. We want to
c     convert this to the origin time
c-----
        call mnthdy(nzyear,nzjday,month,day)
        second = nzsec + 0.001*nzmsec
        call htoe(nzyear,month,day,nzhour,nzmin,second,epoch)
        epoch = epoch + otime
        call etoh(epoch,date,str,nzjday,
     1      nzyear,month,day,nzhour,nzmin,second)
        nzsec = second
        nzmsec = 1000.0*(second - nzsec)
        tpo = tpo - oo
        tso = tso - oo
        bo  = bo  - oo
        oo = 0.0

        call getkhv('KSTNM   ',kstnm ,nerr)
        call getkhv('KCMPNM  ',kcmpnm,nerr)
        call getkhv('KNETWK  ',knetwk,nerr)
        call getfhv('DIST    ',disto ,nerr)
        call getfhv('AZ      ',az    ,nerr)
        call getfhv('USER1   ',user1 ,nerr)
        call getfhv('USER2   ',user2 ,nerr)
 
c-----
c     error check
c-----
        if(abs(dt - dt2).gt. 0.01*dt)then
            write(stdout,*) 'ERROR: dt not equal for ',
     1          obs_name
            write(6,*)dt,dt2
            stop
        endif
c-----
c     set shortest npts for length
c-----
        if(nptsdat .ne. nptsgrn) then
            write(ascfmt,'(a6,i2,a1)')'(a30,a',sleno,')'
            write(ounit,ascfmt) '  WARNING: npts not equal for ',
     1          obs_name
            npts = min0( nptsgrn, nptsdat )
            write(ounit,*) ' Truncating to ',npts,' points.'
        else
            npts = nptsdat
        endif
        return
        end
       

        subroutine getmtresp(i,iwvtype,az,npts)
c-----
c
c   routine that combines the fundamental fault responses of 
c       Herrmann into moment tensor responses
c
c     NOTE AZ is in RADIANS
c-----
        implicit none
        integer i, iwvtype, npts
        real az

        integer maxpts, nprm, nmax, mxwvs
        real deg_to_rad, rad_to_deg
        parameter(maxpts=2048, nprm=6, mxwvs=1000, nmax=mxwvs*maxpts)
        parameter(deg_to_rad = 0.017453292, rad_to_deg = 57.29577951)
       
        integer ixx,iyy,ixy,ixz,izz,iyz
        integer iPz, iPr, iSH
        real c2a, s2a, ca, sa

        real vss, vds, vdd, vex, response
        common/gftns/vss(maxpts),vds(maxpts),vdd(maxpts),vex(maxpts),
     &                response(maxpts)
           
        integer j
              
        ixx = 1
        iyy = 2
        ixy = 3
        ixz = 4
        iyz = 5
        izz = 6

        iPz = 1
        iPr = 2
        iSH = 3

        c2a = cos(2.0 * az)
        s2a = sin(2.0 * az)
        ca = cos(az)
        sa = sin(az)

        if(iwvtype .eq. iPz .or. iwvtype .eq. iPr)then
      
            if(i .eq. ixx)then
                do 10 j = 1, npts
                    response(j)=  0.5*c2a*vss(j)
     1              -vdd(j)/6.0 + vex(j)/3.0
   10           continue
            else if(i .eq. iyy)then
                do 11 j = 1, npts
                    response(j)= -0.5*c2a*vss(j)
     1              -vdd(j)/6.0 + vex(j)/3.0
   11           continue
            else if(i .eq. ixy)then
                do 12 j = 1, npts
                    response(j) = s2a * vss(j)
   12           continue
            else if(i .eq. ixz)then
                do 13 j = 1, npts
                    response(j) = ca * vds(j)
   13           continue
            else if(i .eq. iyz)then
                do 14 j = 1, npts
                    response(j) = sa * vds(j)
   14           continue
            else if(i .eq. izz)then
                do 15 j = 1, npts
                    response(j) = (vdd(j) + vex(j))/3.0
   15           continue
            endif

        else if(iwvtype .eq. iSH) then
      
            if(i .eq. ixx)then
                do 20 j = 1, npts
                    response(j) = 0.5*s2a*vss(j)
   20       continue
            else if(i .eq. iyy)then
                do 21 j = 1, npts
                    response(j) = -0.5*s2a*vss(j)
   21           continue
            else if(i .eq. ixy)then
                do 22 j = 1, npts
                    response(j) = -c2a*vss(j)
   22           continue
            else if(i .eq. ixz)then
                do 23 j = 1, npts
                    response(j) =   sa * vds(j)
   23           continue
            else if(i .eq. iyz)then
                do 24 j = 1, npts
                    response(j) = - ca * vds(j)
   24           continue
            else if(i .eq. izz)then
                do 25 j = 1, npts
                    response(j) = 0.0
   25           continue
            endif

        endif

        return
        end

c******************************************************************
       
        subroutine appendvector(outv,mxpts,ibegin,inv,npts)
c-----
c     append a vector to an existing vector
c-----
        implicit none
        integer mxpts, ibegin, npts
        real inv(npts), outv(mxpts)
        integer  i, k
      
        do 1 i = 1, npts
            k = i + ibegin - 1
            outv(k) = inv(i)
    1   continue

        return
        end
      
c******************************************************************

        subroutine atrans_a(a,nr,nc,nmax,ata)
c-----
c     given the matrix A form the sqaure matrix
c      T
c     A A
c-----
c     A   R   initial matrix
c             of dimention (nmax,nc) only nr rows
c     nr  I   number of rows in matrix A
c     nc  I   number of columns in matrix A
c     nmax    I   Maximum number of rows in A
c     ata R   Resultant matrix
c-----
        implicit none
        integer nr, nc, nmax
        real a(nmax,nc), ata(nc,nc)
        integer r, i, j
      
        do 3 i = 1, nc
            do 2 j = 1, nc
                ata(i,j) = 0.0
                do 1 r = 1, nr
                    ata(i,j) = ata(i,j) + a(r,i)*a(r,j)
    1           continue
    2       continue
    3   continue

        return
        end
      
c******************************************************************    
  
        subroutine atrans_b(a,nr,nc,nmax,b,atb)
c-----
c     given the matrices A and vector B, form
c             T
c     vector A B
c-----
c     A   R   initial matrix
c             of dimention (nmax,nc) only nr rows
c     nr  I   number of rows in matrix A
c     nc  I   number of columns in matrix A
c     nmax    I   Maximum number of rows in A
c     b   R   b - vector with nr or nmax rows
c     atb R   Resultant vector
c-----
        implicit none
        integer nr, nc, nmax
        real a(nmax,nc), atb(nc), b(nr)
        integer  i, r
      
        do 2 i = 1, nc
            atb(i) = 0.0
            do 1 r = 1, nr
                atb(i) = atb(i) + a(r,i)*b(r)
    1       continue
    2   continue
        return
        end
      
c******************************************************************


        subroutine print_cm(cm,n,ounit)
c-----
c     print the COVARIANCE MATRIX
c-----
        implicit none
        integer LIN, LOT, LER
        parameter (LIN=5, LOT=6, LER=0)
        integer n, ounit
        real cm(n,n)
      
        integer i,j
      
      
        write(ounit,90) 
        write(ounit,91)
            
        do 10 i = 1,n
            write(ounit,100) (cm(i,j), j = 1,i)
   10   continue

   90   format(/,'Covariance Matrix' )
   91   format('  CXX         CYY         CXY        CXZ',
     1      '        CYZ         CZZ')

  100   format(6(e10.3,1x))

        write(ounit,*)

        return
        end

       
        subroutine apply_weights(a,nrmax,ir,nr,ic,nc,b,wt)
c-----
c
c routine to scale a set of equations Ax = b by scaling the elements of
c        a and b
c
c   nr = number of rows to scale, nc = ncols
c   ir = first row to scale, ic = first col to scale
c
c   nrmax = maximum dimension of a and b in main program
c
c-----
        implicit none
        integer nrmax,nr,nc,ic,ir
        real a(nrmax,nc),b(nrmax),wt
      
        integer i,j,n
      
        n = ir + nr
      
        do 2 i = ir, n
            b(i) = b(i) * wt
            do 1 j = ic, nc
                a(i,j) = a(i,j) * wt
    1       continue
    2   continue
        return
        end

c-----
c
c   routine to compute the deviatoric moment tensor for an 
c      input strike,dip, and rake
c
c   set up to return the mij in lagnston (1981) notation
c      which differ from aki and richard notation by a sign on all
c      components
c     THIS IS THE NOTATION IN UDIAS p325
c                             LAY and WALLACE p343
c
c-----
        subroutine dislocation_to_mij(istrike,idip,irake,m)
c-----
c     convert dislocation specification in terms of
c     strike, dip and rake to a moment tensor
c
c     These definitions agree with Udias (17.24) and Lay and Wallace
c
c     Note for a dislocation source, M11 + M22 + M33 = 0
c-----
        implicit none
        integer istrike,idip,irake
        real m(6)
c
        real d2r
        parameter(d2r = 0.017453292)

        integer ixx, iyy, ixy, ixz, iyz, izz
        real sind, cosd, sin2d, cos2d
        real sins, coss, sin2s, cos2s
        real sinr, cosr

        ixx = 1
        iyy = 2
        ixy = 3
        ixz = 4
        iyz = 5
        izz = 6

        sind   =  sin(d2r * idip)
        cosd   =  cos(d2r * idip)
        sin2d  =  sin(d2r * 2 * idip)
        cos2d  =  cos(d2r * 2 * idip)
        
        sins   =  sin(d2r * istrike)
        coss   =  cos(d2r * istrike)
        sin2s  =  sin(d2r * 2 * istrike)
        cos2s  =  cos(d2r * 2 * istrike)

        sinr   =  sin(d2r * irake)
        cosr   =  cos(d2r * irake)
        
        m(ixx) = -sind*cosr*sin2s - sin2d*sinr*sins*sins
        m(iyy) =  sind*cosr*sin2s - sin2d*sinr*coss*coss
        m(ixy) =  sind*cosr*cos2s + 0.5*sin2d*sinr*sin2s
        m(ixz) = -cosd*cosr*coss  - cos2d*sinr*sins

C---- 
C       changes 22 SEP 02 to make agree with fmplot 
C       m(iyz) =  cosd*cosr*sins  - cos2d*sinr*coss
c-----
        m(iyz) = -cosd*cosr*sins  + cos2d*sinr*coss
        m(izz) =  sin2d*sinr

C       WRITE(0,*)'STK:',istrike
C       WRITE(0,*)'DIP:',idip
C       WRITE(0,*)'RAK:',irake
C       WRITE(0,*)'Mxx:',m(ixx)
C       WRITE(0,*)'Myy:',m(iyy)
C       WRITE(0,*)'Mxy:',m(ixy)
C       WRITE(0,*)'Myz:',m(iyz)
C       WRITE(0,*)'Mxz:',m(ixz)
C       WRITE(0,*)'Mzz:',m(izz)
      
        return
        end

        subroutine getevent_data(fname,ounit,pname)
c-----
c     routine to output some event data using first file
c-----
        implicit none
        character*256 fname
        integer ounit
        character pname*(*)
      
        real origin_time,obs
        integer slen1,slen2
        character*32 ref_date, ref_time,ascfmt
        integer nz(6)

        integer ihr, imin, isec, imsec
        integer nerr
        real btime, dt2
        integer nptsdat, nptsgrn
        integer ls
        integer lgstr

        common/event/evla, evlo
        real evla, evlo

c-----
c   read in the data header
c-----
        call rsac1(fname,obs,nptsdat,btime,dt2,1,nerr) 
 
        call getfhv('EVLA',evla,nerr)
        call getfhv('EVLO',evlo,nerr)
        call getfhv('O',origin_time,nerr)
      
        call getnhv('NZYEAR',nz(1),nerr)
        call getnhv('NZJDAY',nz(2),nerr)
        call getnhv('NZHOUR',nz(3),nerr)
        call getnhv('NZMIN',nz(4),nerr)
        call getnhv('NZSEC',nz(5),nerr)
        call getnhv('NZMSEC',nz(6),nerr)
      
        ihr = origin_time / 3600
        imin = (origin_time - 3600*ihr) / 60
        isec = (origin_time - 3600*ihr - 60*imin)
        imsec = 1000* (origin_time - 3600*ihr - 60*imin - isec)
      
        call kadate(nz(1),nz(2),32,ref_date,nerr)
        call katime(nz(3),nz(4),nz(5),nz(6),32,ref_time,nerr)
      
        ls = lgstr(pname)
        write(ounit,*) 'Moment Tensor Inversion using ', pname(1:ls)
        write(ounit,*)' '

        write(ounit,*)' '
        write(ounit,*)'Information from header of first data file: '
       
        slen1 = lgstr(ref_date)
        slen2 = lgstr(ref_time)
       
        write(ascfmt,'(a6,i2,a5,i2,a1)')
     &      '(a17,a'  , slen1 ,  ',4x,a' ,  slen2,   ')'
    
        write(ounit,ascfmt)'Reference Time: ',ref_date, ref_time
        write(ounit,*)'Origin time (s) relative to reference: ',
     1      origin_time
        write(ounit,'(a16,f7.3,2x,f8.3)')' Event Lat,Lon: ',evla, evlo
        write(ounit,*)' '
      
        return
        end

        subroutine getrr(rr, obs, pre, npts, wvtype, nwaves)
        implicit none
c-----
c     compute correlation coefficient between observed and
c     predicted for Z, R and T separately
c-----
        integer maxpts, nprm, nmax, mxwvs
        parameter(maxpts=2048, nprm=6, mxwvs=1000, nmax=mxwvs*maxpts)

        real rr(3)
        real obs(nmax), pre(nmax)
        
        integer npts(mxwvs), wvtype(mxwvs), nwaves

        real sum0(3),  sumxx(3), sumxy(3), sumyy(3)
        integer i, ibegin, j, k, m
c-----
c     initialize
c-----
        do 1000 i=1,3
            sum0(i)  = 0.0
            sumxx(i) = 0.0
            sumxy(i) = 0.0
            sumyy(i) = 0.0
 1000   continue

        ibegin = 1
        do 2000 i=1,nwaves
            j = wvtype(i)
            do 2100 k=1,npts(i)
                m = ibegin + k - 1
                sum0(j) = sum0(j) + 1
                sumxx(j) = sumxx(j) + obs(m)*obs(m)
                sumxy(j) = sumxy(j) + obs(m)*pre(m)
                sumyy(j) = sumyy(j) + pre(m)*pre(m)
 2100       continue
            ibegin = ibegin +npts(i)
 2000   continue
c-----
c     compute vector dot product
c-----
        do 3000 i=1,3
            if(sumxx(i).gt.0.0 .and. sumyy(i).gt.0.0)then
                rr(i) = sumxy(i)/(sqrt(sumxx(i))*sqrt(sumyy(i)))
            else
                rr(i) = -2.0
            endif
 3000   continue
        return
        end
                
        subroutine vmul(x,n,fac)
c-----
c     multiply a vector by a scalar
c-----
        implicit none
        integer n
        real x(n)
        real fac
        integer i
        do 1000 i=1,n
            x(i) = x(i) * fac
 1000   continue
        return
        end
                
        subroutine out_mte(tensor,var,covmatrix,nprm,m0,ounit)
c-----
c     output moment tensor with error bounds
c
c     tensor      R*4 moment tensor (in units of m0)
c     var     R*4 variance of fit = SUM err^2 / ndf
c     covmatrix   R*4 variance/covariance matrix
c     nprm        I   array dimension - will be 6
c     m0      R*4 Scale factor The CPS synthetics
c                 are for a moment of 1.0E+20 dyne-cm
c     ounit       I   output file unit
c-----
        implicit none
        real tensor(6)
        real var
        integer nprm
        real covmatrix(nprm,nprm)
        real m0
        integer ounit
        character*4 mijstr(6)
        integer i
        real fnorm, fnormn, fnormd
        data mijstr/' xx ',' yy ',' xy ',' xz ',' yz ',' zz '/
        write(ounit,*)' '
        write(ounit,'(a,a,a)')' ij ',
     1      '      Mij   ',
     2      '  StdErr Mij'
        fnormn = 0
        fnormd = 0
        do 1000 i=1,6
            write(ounit,'(a,2e12.4)')mijstr(i),
     1          tensor(i)*m0,
     2          sqrt(abs(covmatrix(i,i))*var)*m0    
            fnormn = fnormn + abs(covmatrix(i,i)*var)
            fnormd = fnormd + tensor(i)*tensor(i)
c-----
c         count the off diagonal elements twice
c-----
            if(i.eq.3)then
                fnormn = fnormn + abs(covmatrix(i,i)*var)
                fnormd = fnormd + tensor(i)*tensor(i)
            endif
            if(i.eq.4)then
                fnormn = fnormn + abs(covmatrix(i,i)*var)
                fnormd = fnormd + tensor(i)*tensor(i)
            endif
            if(i.eq.5)then
                fnormn = fnormn + abs(covmatrix(i,i)*var)
                fnormd = fnormd + tensor(i)*tensor(i)
            endif
 1000   continue
        fnorm = sqrt(fnormn)/ sqrt(fnormd)
        write(ounit,'(a,1x,e12.4)')' Fnorm = ',fnorm
        return
        end

        function vdot(x,y,n)
        integer n
        real x(n), y(n)
        vdot = 0.0
        do i=1,n
            vdot = vdot + x(i)*y(i)
        enddo
        return
        end

        function vdots(x,y,n,k)
c-----
c     compute a dot product but with a time shift of
c     k samples, e.g., we shift x
c
c     y:  01 02 03 04 05 06 07
c     x:  01 02 03 04 05 06 07  if k= 0
c         02 03 04 05 06 07     if k= 1
c            01 02 03 04 05 06  if k=-1
c     this is really a cross-correlation
c-----
        integer n
        real x(n), y(n)
        if(k.le.0)then
            ilw = 1 - k
            iup = n
        else
            ilw = 1
            iup = n - k
        endif
        vdots = 0.0
        do i=ilw, iup
            vdots = vdots + y(i)*x(i+k)
        enddo
        return
        end

