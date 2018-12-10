        program fdelaz
c-----
c  return AZ, BAZ, DELKM or DELDEG from command line of
c  
c  delaz -ELAT elat -ELON elon -SLAT slat -SLON slon 
c     [-AZ -BAZ -DELKM -DELDEG]
c-----
        integer LER, LOT, LIN
        parameter (LER=0, LIN=5, LOT=6)
        logical dodelkm, dodeldg, doaz, dobaz
        real*4 elat, elon, slat, slon, del, az, baz, delkm, saz, caz
        call gcmdln(elat,elon,slat,slon,dodelkm,dodeldg,doaz,dobaz)
        call delaz(elat,elon,slat,slon,del,az,baz,delkm,saz,caz)
        write(LOT,*)elat,elon,slat,slon,dodelkm,dodeldg,doaz,dobaz
        write(LOT,*)del,delkm,az,baz
        end

        subroutine gcmdln(elat,elon,slat,slon,dodelkm,
     1      dodeldg,doaz,dobaz)
        logical dodelkm, dodeldg, doaz, dobaz
        real*4 elat, elon, slat, slon
        character names*20
        integer mnmarg
        
        dodelkm = .false.
        dodeldg = .false.
        doaz = .false.
        dobaz = .false.
        elat = 0.0
        elon = 0.0
        slat = 0.0
        slon = 0.0

        nmarg = mnmarg()
        i = 0
 1000   continue
            i = i + 1
            if(i.gt.nmarg)go to 2000
            call mgtarg(i,names)
            if(names(1:5).eq.'-ELAT' 
     1          .or. names(1:5).eq.'-elat')then
                i = i + 1
                call mgtarg(i,names)
                read(names,'(bn,f20.0)')elat
            else if(names(1:5).eq.'-ELON'
     1          .or.names(1:5).eq.'-elon')then
                i = i + 1
                call mgtarg(i,names)
                read(names,'(bn,f20.0)')elon
            else if(names(1:5).eq.'-SLAT'
     1          .or.names(1:5).eq.'-slat')then
                i = i + 1
                call mgtarg(i,names)
                read(names,'(bn,f20.0)')slat
            else if(names(1:5).eq.'-SLON'
     1          .or.names(1:5).eq.'-slon')then
                i = i + 1
                call mgtarg(i,names)
                read(names,'(bn,f20.0)')slon
            else if(names(1:3).eq.'-AZ' 
     1          .or. names(1:3).eq.'-az')then
                doaz = .true.
            else if(names(1:4).eq.'-BAZ' 
     1          .or.names(1:4).eq.'-baz')then
                dobaz = .true.
            else if(names(1:6).eq.'-DELKM'
     1          .or.names(1:6).eq.'-delkm')then
                dodelkm = .true.
            else if(names(1:7).eq.'-DELDEG'
     1          .or.names(1:7).eq.'-deldeg')then
                dodeldg = .true.
            else if(names(1:2).eq.'-h' 
     1          .or. names(1:2).eq.'-?')then
                call usage()
            endif
        go to 1000
 2000   continue
        return
        end
                
        subroutine usage()
        stop
        end


        subroutine delaz(elat,elon,slat,slon,del,az,baz,delkm,saz,caz)
c-----
c     elat    R*4 Epicenter latitude (degrees)
c     elon    R*4 Epicenter longitude (degrees)
c     slat    R*4 Station latitude (degrees)
c     slon    R*4 Station longitude (degrees)
c     del R*4 Epicentral Distance (degrees)
c     az  R*4 Epicenter to Station Azimuth
c     baz R*4 Station to Epicenter Backazimuth
c     delkm   R*4 Epicentral Distance (km)
c     saz R*4 Sine of Azimuth
c     caz R*4 Cosine of Azimuth
c-----
        real*4 rad, e2, re
        data rad/0.017453292/,e2/0.9932315/, re/6371.003/
        call dircos(elat,elon,selat,celat,selon,celon,ea,eb,ec)
        call dircos(slat,slon,sslat,cslat,sslon,cslon,sa,sb,sc)

c-----
c     compute distance
c     Choose correct formula for short and large distances
c-----
        cdel = (ea*sa + eb*sb + ec*sc)
c-----
c     if DEL = [0,20)
c-----
        if(cdel .gt. 0.9396)then
            fac = (ea-sa)**2 + (eb-sb)**2 + (ec-sc)**2
            fac = sqrt(fac)/2.0
            del = 2.0*asin(fac)
c-----
c     if DEL = [20,160]
c-----
        else if(cdel. le. 0.9396 .and. cdel. ge. -0.9396)then
            del = acos(cdel)
c-----
c     if DEL = (160,180]
c-----
        else
            fac = (ea+sa)**2 + (eb+sb)**2 + (ec+sc)**2
            fac = sqrt(fac)/2.0
            del = 2.0*acos(fac)
        endif
        
c-----
c     check for station or epicenter at pole
c-----
        if(elat.eq.90.0 .or. slat.eq.-90.0)then
            az = 180.0
            baz = 0.0
            saz =  0.0
            caz = -1.0
        else if(elat.eq.-90.0 .or. slat.eq.90.0)then
            az = 0.0
            baz = 180.0
            saz = 0.0
            caz = 1.0
        else
            saz = celat*(cslat * sin(rad*(slon - elon)))
            caz = (sslat - cdel*selat)
            fac = sqrt(saz**2 + caz**2)
            if(fac.gt.0.0)then
                saz = saz / fac
                caz = caz / fac
                az = atan2(saz,caz)
            
                sbz = - cslat*(celat * sin(rad*(slon - elon)))
                cbz = (selat - cdel*sslat)
                baz = atan2(sbz,cbz)
            else
                az = 0.0
                caz = 1.0
                saz = 0.0
                baz = 180.0
            endif
            az = az / rad
            baz = baz / rad
        endif
        delkm = del * re
        del = del / rad

c-----
c     put az and baz in the range [0,360)
c-----
        if(az .lt. 0.0)az = az + 360.0
        if(baz .lt. 0.0)baz = baz + 360.0
        return
        end
        
        subroutine dircos(lat,lon,slat,clat,slon,clon,aa,bb,cc)
        real*4 lat,lon,slat,slon,clat,clon, aa, bb, cc
        real*4 rad, e2
        data rad/0.017453292/,e2/0.993305615/
c-----
c     convert geographic latitude to geocentric
c     Use flattening of Chovitz (1981) f= 1/298.257 adopted by IUGG 1980
c     
c     The relation between geocentric and geographic latitude is
c     tan phi c = ( 1 - f)^2 tan phi g
c
c     To avoid problems at poles, define sin phi c and cos phi c
c     so that the relation holds ans also that s^2 + c^2 = 1
c
c     For geographic to geocentric use e2 = (1-f)^2 = 0.993395615
c     For geographic to equidistant use e2 =(1-f)^1.5=0.99776354
c     Brown, R. J. (1984). On the determination of source-receiver
c         distances using a new equidistant latitude,
c         Geophys. J. R. astr. Soc. 76, 445-459.
c     
c-----
        e2 = 0.99776354
        c = cos(rad*lat)
        s = sin(rad*lat)
        e4 = e2**2
        fac = sqrt(e4 + (1.0-e4)*c*c)
        slat = e2 * s /fac
        clat =      c /fac
        slon = sin(rad*lon)
        clon = cos(rad*lon)

        aa = clat * clon
        bb = clat * slon
        cc = slat
        
        return
        end
        function mnmarg()
c---------------------------------------------------------------------c
c                                                                   c
c    COMPUTER PROGRAMS IN SEISMOLOGY                                c
c    VOLUME V                                                       c
c                                                                   c
c    SUBROUTINE: MNMARG                                             c
c                                                                   c
c    COPYRIGHT 1996 R. B. Herrmann                                  c
c                                                                   c
c    Department of Earth and Atmospheric Sciences                   c
c    Saint Louis University                                         c
c    221 North Grand Boulevard                                      c
c    St. Louis, Missouri 63103                                      c
c    U. S. A.                                                       c
c                                                                   c
c---------------------------------------------------------------------c
        integer mnmarg
            mnmarg = iargc()
        return
        end
        subroutine mgtarg(i,name)
c---------------------------------------------------------------------c
c                                                                   c
c    COMPUTER PROGRAMS IN SEISMOLOGY                                c
c    VOLUME V                                                       c
c                                                                   c
c    SUBROUTINE: MGTARG                                             c
c                                                                   c
c    COPYRIGHT 1996 R. B. Herrmann                                  c
c                                                                   c
c    Department of Earth and Atmospheric Sciences                   c
c    Saint Louis University                                         c
c    221 North Grand Boulevard                                      c
c    St. Louis, Missouri 63103                                      c
c    U. S. A.                                                       c
c                                                                   c
c---------------------------------------------------------------------c
c-----
c     return the i'th command line argument
c
c     This version works on SUN, IBM RS6000
c-----
        character name*(*)
            call getarg(i,name)
        return
        end
