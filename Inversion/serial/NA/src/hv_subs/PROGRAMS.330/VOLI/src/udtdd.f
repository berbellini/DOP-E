c----
c     return P-wave dt/dD in sec/km as function of 
c         epicentral distance and source depth
c
c     udtdd   -D DELDEG -H DEPTH_MK
c     udtdd   -S binary_sacfile_with_headers_set
c
c
c----
        call gcmdln(delta, depth)
c----
c     
c----
        ideg = delta
        if(ideg .lt. 1 .or. ideg.gt.105)then
            dtdd = -12345.0
            write(6,*)dtdd
        else
                call jbtime(delta,1.0,0.0,time,dtdlat,dtdlon,dtdel2,
     1          dtdz,depth,'P       ',anin)
            dtdd = dtdel2/111.195
            write(6,*)dtdd
        endif
        end

        subroutine gcmdln(delta,depth)
c----
c     delta   R   - epicentral distance in degrees
c     depth   R   - source depth in km
c----
        real delta, depth
        character names*80
        integer mnmarg
        integer i, nmarg

        i = 0
        nmarg = mnmarg()

        delta = 10.0
        depth = 10.0
 1000   continue
            i = i + 1
            if(i.gt.nmarg)go to 2000
            call mgtarg(i,names)
            if(names(1:6).eq.'-GCARC')then
                i = i + 1
                call mgtarg(i,names)
                read(names,'(bn,f10.0)')delta
            else if(names(1:5).eq.'-EVDP')then
                i = i + 1
                call mgtarg(i,names)
                read(names,'(bn,f10.0)')depth
            else if(names(1:2).eq.'-h' .or. names(1:2).eq.'-?')then
                call usage()
            endif
        go to 1000
 2000   continue
        return
        end

        subroutine usage()
        integer LER, LOT, LIN
        parameter (LER=0, LIN=5, LOT=6)
        write(LER,*)'Usage: udtdd -GCARC gcarc -EVDP evdp [-h] [-?]'
        write(LER,*)'   Compute P-wave ray parameter in sec/km'
        write(LER,*)'   Example: udtdd -GCARC 31.0 -EVDP 170.0'
        write(LER,*)
     1      ' -GCARC gcarc   (default 10.0) epicentral dist deg'
        write(LER,*)
     1      ' -EVDP evdp     (default 10.0) source depth  km'
        write(LER,*)
     1      ' -h                            this command help'
        write(LER,*)
     1      ' -?                            this command help'
        
        stop
        end

        subroutine jbtime(del,caz,saz,time,dtdlat,dtdlon,dtdel,
     1      dtdz,h,phase,anin)
c----
c     del R*4 - distance in degrees
c     caz R*4 - cos(azimuth)
c     saz R*4 - sin(azimuth)  
c     time    R*4 - travel time in seconds (returned)
c     dtdlat  R*4 - change in travel time for one km NS position change 
c     dtdlon  R*4 - change in travel time for one km EW position change
c     dtdel   R*4 - change in time for one degree in distance
c     h   R*4 - focal depth in km
c     phase   C*8 - Phase identification, e.g., 'P      '
c     anin    R*8 - sin of takeoff angle
c----
        real*4 del, caz, saz, time, dtdlat, dtdlon, dtdel, h
        character phase*8

      common /tabl/ tp,ts,bofp,bofs,cp,cs,pv,sv
      real*4 tp(14,21),ts(14,24),bofp(17),bofs(20),cp(21),cs(24)
      real*4 pv(14),sv(14)
c----
c     get velocity at source depth
c----
      if(h.le.33.) then
        j=1
        vp=pv(j)+(pv(j+1)-pv(j))*h/33.
        vs=sv(j)+(sv(j+1)-sv(j))*h/33.
      else
        do 310 kl=1,12
           ak=kl
           q=33.+ak*63.38
           if(h.le.q)go to 315
  310      continue
  315      j=kl+1
           r=h-q
           r = abs(r)
           vp=pv(j+1)+(pv(j)-pv(j+1))*r/63.38
           vs=sv(j+1)+(sv(j)-sv(j+1))*r/63.38
      endif
c------------------------------------------------------------------
c   calculate the p and s arrival times 
c------------------------------------------------------------------
        call pstime(h,del,0,pt,st)
        if(h.lt.60.0)then
            call pstime(h     ,del,0,ptm,stm)
            call pstime(h+30.0,del,0,ptp,stp)
            dh = 30.0
        else if(h.ge.60.0 .and. h.le.730.0)then
            call pstime(h-63.38,del,0,ptm,stm)
            call pstime(h+63.38,del,0,ptp,stp)
            dh = 2.0*63.38
        else if(h.gt.793.0)then
            call pstime(h     ,del,0,ptp,stp)
            call pstime(h-30.0,del,0,ptm,stm)
            dh = 30.0
        endif
            
        if(phase.eq.'P')then
            time = pt
            dtdz = (ptp - ptm)/(dh)
        else if(phase.eq.'S')then
            time = st
            dtdz = (stp - stm)/(dh)
        endif
c------------------------------------------------------------------
c   The application of the J-B table is limited. If the depth is
c   too large, it is applied to smaller distance.
c------------------------------------------------------------------
c     S does not exist
c----
        if(h.gt.97..and.del.gt.107.) go to 147 
        if(h.gt.287..and.del.gt.106.) go to 147 
        if(h.gt.540..and.del.gt.105.) go to 147 
        if(h.gt.667..and.del.gt.104.) go to 147 
        if(h.gt.800..and.del.gt.103.) go to 147 
        goto 148
  147   continue
            st = -1.0
  148   continue
c----
c     P does not exist
c----

        if(h.gt.33..and.del.gt.105.) go to 145 
        if(h.gt.287..and.del.gt.104.) go to 145 
        if(h.gt.540..and.del.gt.103.) go to 145 
        if(h.gt.800..and.del.gt.102.) go to 145 
        go to 146
  145   continue
            pt = -1.0
  146   continue
      call psang(del,san,pan,h,vp,vs,dtddp,dtdds)
c----
c     get information on d2d/dDelta2
c     If positive ray goes up and anin > 0
c     If negative ray goes down and anin < 0
c----
      call psang(del+1.0,san1,pan1,h,vp,vs,dtddp1,dtdds1)
        if(phase.eq.'P')then
            dtdlat = -dtddp * caz / 111.195
            dtdlon = -dtddp * saz / 111.195
            dtdel = dtddp 
            anin = sin(3.1415927*pan/180.0)
            if(dtddp1 .lt. dtddp)anin = - anin
        else if(phase.eq.'S')then
            dtdlat = -dtdds * caz / 111.195
            dtdlon = -dtdds * saz / 111.195
            dtdel = dtdds 
            anin = sin(3.1415927*san/180.0)
            if(dtdds1 .lt. dtdds)anin = - anin
        endif
        return
      end

      subroutine pstime(h,angi,nfd,pt,st)
c-----------------------------------------------------------------------
c   this subroutine calculate the p & s travel time 
c-----------------------------------------------------------------------
c    Cubic spline interpolation is employed to determine the travel time 
c-----------------------------------------------------------------------
            call bfit(h,angi,nfd,pt,st)
      return
      end

      subroutine psang(dis,san,pan,h,vp,vs,dtddp,dtdds)
c----------------------------------------------------------------------
c   employ bfit subroutine to calculate the 1st derivative
c----------------------------------------------------------------------
      call bfit(h,dis,1,yp,ys)
      degrad = 57.29577951
        dtdds = ys
      u = ys*vs*degrad/(6371.-h)
      u = abs(u)
      if(u.lt.1.) then
         san = asin(u)
      else
         san = 1.5708
      endif
      san=san*57.2958

        dtddp = yp
      c=yp*vp*degrad/(6371.-h)
      c = abs(c)
      if(c.lt.1.) then
         pan = asin(c)
      else
         pan = 1.5708
      endif
      pan=pan*57.2958
      return
      end

      subroutine bfit(h,dis,nfd,yp,ys)
      dimension ut(4),ss(4)
      dimension tp(14,21),ts(14,24),bofp(17),bofs(20),cp(21),cs(24)
      dimension pv(14),sv(14)
      common /tabl/ tp,ts,bofp,bofs,cp,cs,pv,sv
            p1(t)=.25*t**3
            p2(t)=-(1.-t)**2*(1.+t)*.75+1.
            cp1(t)=t**4/16.
            cp2(t)=t*(1./4.+t*(3./8.+t*(1./4.-3./16.*t)))
            dp1(t)=.75*t**2
            dp2(t)=.75 + t*(1.5-2.25*t)
            ddp1(t)=1.5*t
            ddp2(t)=-4.5*t+1.5
c-----------------------------------------------------------------------
c     calculate the coefficients at particular depth
c     if the depth is not exactly at 0.01R, 0.02R ...
c     we use the linear interporation for those coefficients
c-----------------------------------------------------------------------
         if (h .le. 0.) then
              do 11 j=1,21
              cp(j)=tp(1,j)/10.
  11         continue
              do 12 j=1,24
              cs(j)=ts(1,j)/10.
  12         continue
         elseif ( h .gt. 0. .and. h .le. 33.) then
             do 13 j=1,21
             ptemp=(h - 0.)*(tp(2,j)-tp(1,j))/33.
             cp(j) = (tp(1,j) + ptemp)/10.
  13         continue
             do 14 j=1,24
             stemp=(h-0.)*(ts(2,j)-ts(1,j))/33.
             cs(j)=(ts(1,j)+stemp)/10.
  14         continue
         elseif (h.gt.33..and.h.le.797.52) then
            m=0
            r1=63.71
  111       r2=33.+m*r1
            r3=33.+(m+1)*r1
             if (h.gt.r2.and. h.le.r3) then
                do 20 j=1,21
                ptemp=(h-r2)*(tp(3+m,j)-tp(2+m,j))/r1
                cp(j)=(tp(2+m,j)+ptemp)/10.
   20           continue
                do 21 j=1,24
                stemp=(h-r2)*(ts(3+m,j)-ts(2+m,j))/r1
                cs(j)=(ts(2+m,j)+stemp)/10.
   21           continue
             else 
                m=m+1
                if (m .lt. 12) go to 111
             endif
         elseif (h .gt. 797.52) then
              do 25 j=1,21
               cp(j)=(tp(14,j))/10.
   25         continue
              do 26 j=1,24
               cs(j)=(ts(14,j))/10.
   26         continue
         endif
c-----------------------------------------------------------------------
c     calculate the break points
c-----------------------------------------------------------------------
          h1p=(35.0-0.0)/9.
          bofp(1)=0.
          do 30 n=1,9
          bofp(n+1)=bofp(n)+h1p
  30     continue
          h2p=(107.0-35.0)/6.
          bofp(11)=35.0
          do 31 n=11,16
          bofp(n+1)=bofp(n)+h2p
  31    continue
          h1s=(70.0-0.0)/14.
          bofs(1)=0.
          do 32 n=1,14
          bofs(n+1)=bofs(n)+h1s
  32     continue
          h2s=(109.0-70.0)/4.
          bofs(16)=70.0
          do 33 n=16,19
          bofs(n+1)=bofs(n)+h2s
  33    continue
c-----------------------------------------------------------------------
c     determine the distance between particular break points
c-----------------------------------------------------------------------
        if ( dis .le. 35.0 ) then
              jn=1
 100          if( dis .le. bofp(jn+1)) go to 10
               jn=jn+1
               go to 100
 10            u=(dis-bofp(jn))/h1p
        else
              jn=11
 200          if( dis .le. bofp(jn+1)) go to 210
               jn=jn+1
               go to 200
 210           u=(dis-bofp(jn))/h2p
        endif
        if ( dis .le. 70.0 ) then
              jns=1
 300          if( dis .le. bofs(jns+1)) go to 310
               jns=jns+1
               go to 300
 310            us=(dis-bofs(jns))/h1s
        else
              jns=16
 400          if( dis.le.bofs(jns+1)) go to 410
               jns=jns+1
               go to 400
 410           us=(dis-bofs(jns))/h2s
        endif
c-----------------------------------------------------------------------
c  calculate 
c     i. nfd=0 -- interpolation-- travel time (sec)
c    ii. nfd=1 -- 1st derivative --slowness (sec/degree)
c   iii. nfd=2 -- 2nd derivatives --(sec/degree)**2
c-----------------------------------------------------------------------
        v=1.-u
        vv=1.-us
          if(nfd.eq.0) then
             ut(1)=p1(v)
             ut(2)=p2(v)
             ut(3)=p2(u)
             ut(4)=p1(u)
             ss(1)=p1(vv)
             ss(2)=p2(vv)
             ss(3)=p2(us)
             ss(4)=p1(us)
          elseif (nfd.eq.1) then
            if(dis .le. 35.0) then
              ut(1)=-dp1(v)/h1p
              ut(2)=-dp2(v)/h1p
              ut(3)=dp2(u)/h1p
              ut(4)=dp1(u)/h1p
            elseif(dis .gt. 35.0) then
              ut(1)=-dp1(v)/h2p
              ut(2)=-dp2(v)/h2p
              ut(3)=dp2(u)/h2p
              ut(4)=dp1(u)/h2p
            endif
            if(dis .le. 70.0) then
              ss(1)=-dp1(vv)/h1s
              ss(2)=-dp2(vv)/h1s
              ss(3)=dp2(us)/h1s
              ss(4)=dp1(us)/h1s
            elseif(dis .gt. 70.0) then
              ss(1)=-dp1(vv)/h2s
              ss(2)=-dp2(vv)/h2s
              ss(3)=dp2(us)/h2s
              ss(4)=dp1(us)/h2s
            endif
          elseif (nfd.eq.2) then
            if (dis .le. 35.0) then
             ut(1)=ddp1(v)/h1p**2
             ut(2)=ddp2(v)/h1p**2
             ut(3)=ddp2(u)/h1p**2
             ut(4)=ddp1(u)/h1p**2
            else 
             ut(1)=ddp1(v)/h2p**2
             ut(2)=ddp2(v)/h2p**2
             ut(3)=ddp2(u)/h2p**2
             ut(4)=ddp1(u)/h2p**2
            endif
            if (dis .le. 70.0) then
             ss(1)=ddp1(vv)/h1s**2
             ss(2)=ddp2(vv)/h1s**2
             ss(3)=ddp2(us)/h1s**2
             ss(4)=ddp1(us)/h1s**2
            else 
             ss(1)=ddp1(vv)/h2s**2
             ss(2)=ddp2(vv)/h2s**2
             ss(3)=ddp2(us)/h2s**2
             ss(4)=ddp1(us)/h2s**2
            endif
          elseif (nfd.eq.3) then
            if (dis .le. 35.0) then
             ut(1)=-cp1(v)*h1p
             ut(2)=-cp2(v)*h1p
             ut(3)=cp2(u)*h1p
             ut(4)=cp1(u)*h1p
            else 
             ut(1)=-cp1(v)*h2p
             ut(2)=-cp2(v)*h2p
             ut(3)=cp2(u)*h2p
             ut(4)=cp1(u)*h2p
            endif
            if (dis .le. 70.0) then
             ss(1)=-cp1(vv)*h1s
             ss(2)=-cp2(vv)*h1s
             ss(3)=cp2(us)*h1s
             ss(4)=cp1(us)*h1s
            else 
             ss(1)=-cp1(vv)*h2s
             ss(2)=-cp2(vv)*h2s
             ss(3)=cp2(us)*h2s
             ss(4)=cp1(us)*h2s
            endif
          endif       
        if (dis .gt. 35.0) jn=jn+2
        if (dis .gt. 70.0) jns=jns+2
      yp=cp(jn)*ut(1)+cp(jn+1)*ut(2)+cp(jn+2)*ut(3)+cp(jn+3)*ut(4)
      ys=cs(jns)*ss(1)+cs(jns+1)*ss(2)+cs(jns+2)*ss(3)+cs(jns+3)*ss(4)
        return
        end


      blockdata
      common /tabl/ tp,ts,bofp,bofs,cp,cs,pv,sv
      real*4 tp(14,21),ts(14,24),bofp(17),bofs(20),cp(21),cs(24)
      real*4 pv(14),sv(14)
        data ((tp(i,j),j=1,21),i=1,3)/
     1  -617.,   46.,  417.,  781., 1142., 1481., 1816., 2072., 2319.,
     2  2548., 2775., 2991., 2085., 2779., 3442., 4038., 4554., 5002.,
     3  5373., 5732., 6071.,
     4  -246.,   11.,  397.,  756., 1118., 1456., 1788., 2040., 2289.,
     5  2516., 2746., 2950., 2052., 2747., 3410., 4005., 4519., 4968.,
     6  5337., 5698., 6033.,
     7    65.,   10.,  399.,  745., 1100., 1435., 1754., 1999., 2249.,
     8  2474., 2703., 2906., 2017., 2703., 3366., 3959., 4471., 4918.,
     9  5286., 5647., 5978./
        data ((tp(i,j),j=1,21),i=4,6)/
     1   246.,   47.,  400.,  741., 1083., 1418., 1718., 1963., 2209.,
     2  2436., 2661., 2867., 1969., 2664., 3322., 3913., 4425., 4869.,
     3  5237., 5597., 5929.,
     4   342.,  100.,  407.,  741., 1072., 1402., 1681., 1929., 2170.,
     5  2397., 2621., 2833., 1930., 2626., 3281., 3867., 4379., 4821.,
     6  5188., 5550., 5877.,
     7   394.,  160.,  421.,  742., 1064., 1385., 1646., 1898., 2132.,
     8  2362., 2581., 2800., 1894., 2587., 3241., 3823., 4333., 4775.,
     9  5141., 5504., 5829./
        data ((tp(i,j),j=1,21),i=7,9)/
     1   445.,  216.,  443.,  745., 1062., 1360., 1617., 1867., 2098.,
     2  2327., 2544., 2768., 1869., 2549., 3205., 3782., 4290., 4731.,
     3  5096., 5458., 5779.,
     4   478.,  274.,  467.,  753., 1059., 1333., 1591., 1836., 2066.,
     5  2293., 2512., 2730., 1840., 2516., 3170., 3745., 4249., 4687.,
     6  5052., 5413., 5737.,
     7   508.,  329.,  492.,  763., 1045., 1313., 1569., 1808., 2039.,
     8  2265., 2481., 2701., 1806., 2487., 3136., 3708., 4211., 4645.,
     9  5010., 5370., 5699./
        data ((tp(i,j),j=1,21),i=10,12)/
     1   529.,  382.,  520.,  771., 1038., 1299., 1549., 1785., 2015.,
     2  2238., 2455., 2669., 1790., 2459., 3106., 3674., 4175., 4606.,
     3  4971., 5331., 5659.,
     4   569.,  428.,  551.,  780., 1035., 1290., 1532., 1767., 1994.,
     5  2215., 2431., 2644., 1750., 2439., 3077., 3643., 4139., 4573.,
     6  4928., 5310., 5459.,
     7   596.,  474.,  580.,  794., 1035., 1283., 1518., 1752., 1975.,
     8  2195., 2411., 2619., 1732., 2419., 3051., 3613., 4107., 4538.,
     9  4893., 5276., 5419./
        data ((tp(i,j),j=1,21),i=13,14)/
     1   632.,  517.,  613.,  809., 1038., 1277., 1510., 1738., 1960.,
     2  2178., 2392., 2600., 1723., 2399., 3028., 3584., 4076., 4503.,
     3  4858., 5242., 5381.,
     4   662.,  562.,  644.,  826., 1043., 1274., 1503., 1728., 1948.,
     5  2164., 2376., 2582., 1714., 2383., 3005., 3558., 4046., 4470.,
     6  4823., 5209., 5345./
        data ((ts(i,j),j=1,24),i=1,3)/
     1 -1328.,  100.,  909., 1758., 2550., 3344., 3936., 4469., 4990.,
     2  5500., 5988., 6460., 6914., 7353., 7773., 8173., 8555., 7388.,
     3  8182., 8906., 9548.,10094.,10637.,11164.,
     4  -656.,   20.,  886., 1714., 2514., 3297., 3884., 4415., 4936.,
     5  5446., 5934., 6404., 6859., 7297., 7716., 8115., 8501., 7324.,
     6  8126., 8847., 9490.,10034.,10578.,11104.,
     7    -2.,   -6.,  891., 1684., 2487., 3238., 3807., 4344., 4860.,
     8  5368., 5855., 6324., 6777., 7215., 7633., 8032., 8412., 7263.,
     9  8038., 8761., 9399., 9943.,10487.,11023./
        data ((ts(i,j),j=1,24),i=4,6)/
     1   370.,   45.,  893., 1666., 2460., 3181., 3734., 4275., 4787.,
     2  5295., 5778., 6247., 6699., 7135., 7551., 7949., 8324., 7176.,
     3  7956., 8675., 9311., 9855.,10398.,10927.,
     4   630.,  122.,  901., 1655., 2437., 3119., 3667., 4205., 4718.,
     5  5224., 5705., 6174., 6622., 7059., 7472., 7870., 8243., 7101.,
     6  7876., 8593., 9225., 9769.,10312.,10841.,
     7   795.,  218.,  916., 1651., 2414., 3056., 3607., 4137., 4651.,
     8  5155., 5635., 6101., 6549., 6984., 7396., 7793., 8165., 7026.,
     9  7799., 8513., 9141., 9685.,10228.,10755./
        data ((ts(i,j),j=1,24),i=7,9)/
     1   877.,  327.,  935., 1655., 2386., 2997., 3550., 4073., 4589.,
     2  5088., 5567., 6032., 6479., 6911., 7323., 7717., 8091., 6954.,
     3  7723., 8436., 9060., 9603.,10147.,10673.,
     4   928.,  437.,  960., 1665., 2349., 2947., 3492., 4014., 4528.,
     5  5025., 5503., 5964., 6411., 6840., 7251., 7643., 8016., 6888.,
     6  7649., 8360., 8984., 9524.,10068.,10593.,
     7  1001.,  535., 1000., 1662., 2314., 2901., 3441., 3961., 4475.,
     8  4969., 5443., 5904., 6348., 6776., 7184., 7575., 7943., 6820.,
     9  7580., 8289., 8908., 9449., 9993.,10517./
        data ((ts(i,j),j=1,24),i=10,12)/
     1  1055.,  631., 1039., 1662., 2287., 2862., 3395., 3914., 4427.,
     2  4916., 5389., 5847., 6290., 6716., 7122., 7512., 7877., 6765.,
     3  7515., 8222., 8837., 9377., 9921.,10446.,
     4  1097.,  725., 1080., 1669., 2265., 2829., 3356., 3874., 4382.,
     5  4870., 5340., 5797., 6237., 6661., 7066., 7451., 7817., 6673.,
     6  7463., 8155., 8776., 9297., 9891.,10019.,
     7  1131.,  815., 1124., 1675., 2250., 2800., 3324., 3841., 4344.,
     8  4828., 5296., 5751., 6187., 6610., 7011., 7394., 7763., 6619.,
     9  7408., 8095., 8712., 9232., 9827., 9951./
        data ((ts(i,j),j=1,24),i=13,14)/
     1  1183.,  897., 1173., 1685., 2240., 2778., 3299., 3813., 4310.,
     2  4792., 5257., 5708., 6142., 6560., 6960., 7342., 7702., 6562.,
     3  7357., 8035., 8656., 9158., 9805., 9394.,
     4  1233.,  976., 1224., 1698., 2236., 2761., 3280., 3789., 4282.,
     5  4759., 5221., 5669., 6101., 6515., 6913., 7291., 7649., 1399.,
     6   890., 1399., 2186., 2943., 3692., 4402./
      data (pv(i),sv(i),i=1,14)/
     1   6.75,  3.90,  7.75,  4.35,  7.94,  4.44,  8.13,  4.54,
     2   8.33,  4.64,  8.54,  4.74,  8.75,  4.85,  8.97,  4.96,
     3   9.50,  5.23,  9.91,  5.46, 10.26,  5.67, 10.55,  5.85,
     4  10.77,  6.00, 10.99,  6.12/
      end

