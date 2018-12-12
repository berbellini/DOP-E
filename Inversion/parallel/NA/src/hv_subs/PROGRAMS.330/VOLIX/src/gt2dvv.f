
        subroutine gt2dvv(lun,mmax,listmd,xmin,xmax,zmin,zmax)
        logical listmd
        character string*80

        integer LIN, LOT, LER
        parameter(LIN=5, LOT=6, LER=0)
c-----
c       common for laterally varying model
c-----
        integer NLAYER, NNODE, NDIST
        parameter (NLAYER=20, NNODE=20, NDIST=100)
        common/latbnd/xc, zc, nc, ii
        real*4 xc(NLAYER,NNODE), zc(NLAYER,NNODE)
        integer*4 nc(NLAYER), ii(NLAYER,NNODE)

        common/latval/vp,vs,rho,qp,etap,qs,etas
        real*4 vp(2,NLAYER), vs(2,NLAYER), rho(2,NLAYER), 
     1      qp(2,NLAYER), etap(2,NLAYER), 
     2      qs(2,NLAYER), etas(2,NLAYER)
c-----
c       xc(i,j), yc(i,j) j=1,nc(i)
c           coordinate pairs of layer boundary
c       ii(i,j) a parameter for interpolation used by Cerveny Program
c           set to -1 to avoid spline interpolation
c       vp(i)   P veloctiy in Layer
c       vs(i)   S velocity in layer
c       rho(i)  density velocity in layer
c       qp(i)   Quality factor in layer
c       qs(i)   Quality factor in layer
c       etap(i)     frequency dependence in layer
c       etas(i)     frequency dependence in layer
c-----

        common/depref/refdep

c-----
c       get model specifically for 2-D isotropic variable velocity
c-----
c-----
c            MLYR (NC(i), i=1,MLYR)
c            ((XC(I,J),ZC(I,J),J=1,NC(i)),I=1,MLYR)
c-----
c-----
c       get comment line
c-----
        read(lun,'(a)')string
        read(lun,'(a)')string
        read(lun,*)mmax, (nc(i),i=1,mmax)
        deplow = 0.0
        dephih = 0.0
        zmin = 0.0
        do 900 i=1,mmax
            read(lun,*)(xc(i,j),zc(i,j),j=1,nc(i))
  900   continue
        do 901 i=1,mmax
            do 902 j=1,nc(i)
                ii(i,j) = -1
                if((deplow+zc(i,j)).lt.dephih)then
                    dephih = deplow+zc(i,j)
                endif
  902       continue
            deplow = dephih
  901   continue
        refdep = 0.0
        irefdp = 0
        read(lun,'(a)')string
        do 1000 j=1,mmax
                read(lun,*,err=9000,end=9000)vd,vvp,vvs,vrho,
     1              vqp,vqs
                vp(1,j) = vvp
                vs(1,j) = vvs
                rho(1,j) = vrho
                qp(1,j) = vqp
                qs(1,j) = vqs

                read(lun,*,err=9000,end=9000)vd,vvp,vvs,vrho,
     1              vqp,vqs
                vp(2,j) = vvp
                vs(2,j) = vvs
                rho(2,j) = vrho
                qp(2,j) = vqp
                qs(2,j) = vqs
                if(vd.lt.0.0)then
                    vd = - vd
                    refdep = refdep + vd
                    irefdp = j
                endif
 1000   continue
        zmax = dephih
    1   format(' LAYER             H      P-VEL     S-VEL   DENSITY  ')
    2   format(' ',i5,5x,4f10.3)
    3   format(' ','-SURFACE ','- - - - - ','- - - - - ',
     1      '- - - - - ','- - - - - -')
        if(mmax.gt.0)then
            if(listmd)then
            ierr = 0
            write(LOT,1)
        write(LOT,*)'MMAX=',mmax
            do 2000 i=1,mmax
                if(i.lt.mmax)then
                    d = zc(i+1,1) - zc(i,1)
                else
                    d = 0.0
                endif

                write(LOT,2)
     1              i,d,vp(1,i),vs(1,i),rho(1,i)
                write(LOT,2)
     1              i,d,vp(2,i),vs(2,i),rho(2,i)
                if(i.eq.irefdp)write(LOT,3)
 2000       continue
            endif
        else 
            ierr = -3
            write(LER,*)'Error in model file'
        endif
        return
 9000   continue
            ierr = -3
            write(LER,*)'Error in model file'
        return
        end

