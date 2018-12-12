        subroutine putmod(wlun,mname,mmax,title,iunit,iiso,iflsph,
     1      idimen,icnvel,lverby)
c-----
c       purpose model input - output for Transverse Isotropic models
c-----
c       CHANGES
c       19 APR 2002 - created
c       03 MAY 2002 - permit write to standard output
c       This model specification is designed to be as 
c           general as possible
c
c       Input lines
c       Line 01: MODEL
c       Line 02: Model Name
c       Line 03: ISOTROPIC or ANISOTROPIC or 
c           TRANSVERSELY ANISOTROPIC
c       Line 04: Model Units, First character is length (k for kilometer
c           second is mass (g for gm/cc), third is time (s for time)
c       Line 05: FLAT EARTH or SPHERICAL EARTH
c       Line 06: 1-D, 2-D or 3-D
c       Line 07: CONSTANT VELOCITY
c       Line 08: open for future use
c       Line 09: open for future use
c       Line 10: open for future use
c       Line 11: open for future use
c       Lines 12-end:   These are specific to the model
c           For ISOTROPIC the entries are
c           Layer Thickness, P-velocity, S-velocity, Density, Qp, Qs,
c           Eta-P, Eta S (Eta is frequency dependence), 
c           FreqRefP, FreqRefP
c-----
cMODEL
cTEST MODEL.01
cISOTROPIC
cKGS
cFLAT EARTH
c1-D
cCONSTANT VELOCITY
cLINE08
cLINE09
cLINE10
cLINE11
c H    VP    VS    RHO    QP    QS    ETAP    ETAS  REFP   REFS
c1.0    5.0 3.0 2.5 0.0 0.0 0.0 0.0 1.0 1.0
c2.0    5.1 3.1 2.6 0.0 0.0 0.0 0.0 1.0 1.0
c7.0    6.0 3.5 2.8 0.0 0.0 0.0 0.0 1.0 1.0
c10.0   6.5 3.8 2.9 0.0 0.0 0.0 0.0 1.0 1.0
c20.0   7.0 4.0 3.0 0.0 0.0 0.0 0.0 1.0 1.0
c40.0   8.0 4.7 3.3 0.0 0.0 0.0 0.0 1.0 1.0
c-----
c-----
c       wlun    I*4 - logical unit for writing model file. This
c                 unit is released after the use of this routine
c       mname   C*(*)   - model name
c       mmax    I*4 - number of layers in the model, last layer is
c                    halfspace
c       title   C*(*)   - title of the model file
c       iunit   I*4 - 0 Kilometer, Gram, Sec
c       iiso    I*4 - 0 isotropic 
c                 1 transversely anisotropic 
c                 2 general anisotropic 
c       iflsph  I*4 - 0 flat earth model
c                 1 spherical earth model
c       idimen  I*4 - 1 1-D
c               - 2 2-D
c               - 3 3-D
c       icnvel  I*4 - 0 constant velocity
c                 1 variable velocity
c       lverby  L   - .false. quiet output
c------
        implicit none

        character mname*(*), title*(*)
        integer wlun
        integer*4 mmax, iunit, iiso, iflsph, idimen, icnvel
        logical lverby
c-----
c       LIN I*4 - logical unit for standard input
c       LOT I*4 - logical unit for standard output
c-----
        integer lun
        integer LIN, LOT, LER
        parameter (LIN=5,LOT=6,LER=0)

        integer NL
        parameter (NL=200)
        common/timodel/d(NL),TA(NL),TC(NL),TL(NL),TN(NL),TF(NL),
     1      TRho(NL),
     2      qa(NL),qb(NL),etap(NL),etas(NL), 
     3      frefp(NL), frefs(NL)
        real d,TA,TC,TN,TL,TF,TRho,qa,qb,etap,etas,frefp,frefs
        real vpv, vph, vsv, vsh, vpf

        common/depref/refdep
        real refdep

        integer lt, lgstr
        real curdep
        integer j

        logical ext
        character cmnt1*110
        character cmnt2*110
        cmnt1(  1: 11) = '      H(KM)'
        cmnt1( 12: 22) = '  VPV(KM/S)'
        cmnt1( 23: 33) = '  VSV(KM/S)'
        cmnt1( 34: 44) = ' RHO(GM/CC)'
        cmnt1( 45: 55) = '         QP'
        cmnt1( 56: 66) = '         QS'
        cmnt1( 67: 77) = '       ETAP'
        cmnt1( 78: 88) = '       ETAS'
        cmnt1( 89: 99) = '      FREFP'
        cmnt1(100:110) = '      FREFS'

        cmnt2(  1: 11) = '           '
        cmnt2( 12: 22) = '  VPH(KM/S)'
        cmnt2( 23: 33) = '  VSH(KM/S)'
        cmnt2( 34: 44) = '  VPF(KM/S)'
        cmnt2( 45: 55) = '           '
        cmnt2( 56: 66) = '           '
        cmnt2( 67: 77) = '           '
        cmnt2( 78: 88) = '           '
        cmnt2( 89: 99) = '           '
        cmnt2(100:110) = '           '

        lt = lgstr(title)
c-----
c       test to see if the file exists
c-----
        if(MNAME(1:6).eq.'stdout' .or. mname(1:6).eq.'STDOUT')then
c-----
c           do not open anything, use standard output
c-----
            lun = LOT
        else
            inquire(file=mname,exist=ext)
            if(ext .and.  lverby)then
                write(LER,*)'Overwriting Existing model File'
            endif
            lun = wlun
c-----
c           open the file
c-----
            open(lun,file=mname,status='unknown',form='formatted',
     1          access='sequential')
            rewind lun
        endif
c-----
c       verify the file type
c-----
c-----
c       LINE 01
c-----
        write(lun,'(a)')'MODEL.01'
c-----
c       LINE 02
c-----
        write(lun,'(a)')title(1:lt)
c-----
c       LINE 03
c-----
        if(iiso.eq.0)then
            write(lun,'(a)')'ISOTROPIC'
        else if(iiso.eq.1)then
            write(lun,'(a)')'TRANSVERSE ISOTROPIC'
        else if(iiso.eq.2)then
            write(lun,'(a)')'ANISOTROPIC'
        endif
c-----
c       LINE 04
c-----
        write(lun,'(a)')'KGS'
c-----
c       LINE 05
c-----
        if(iflsph.eq.0)then
            write(lun,'(a)')'FLAT EARTH'
        else if(iflsph.eq.1)then
            write(lun,'(a)')'SPHERICAL EARTH'
        endif
c-----
c       LINE 06
c-----
        if(idimen.eq.1)then
            write(lun,'(a)')'1-D'
        else if(idimen.eq.2)then
            write(lun,'(a)')'2-D'
        else if(idimen.eq.3)then
            write(lun,'(a)')'3-D'
        endif
c-----
c       LINE 07
c-----
        if(icnvel.eq.0)then
            write(lun,'(a)')'CONSTANT VELOCITY'
        else if(icnvel.eq.1)then
            write(lun,'(a)')'VARIABLE VELOCITY'
        endif
c-----
c       put lines 8 through 11
c-----
        write(lun,'(a)')'LINE08'
        write(lun,'(a)')'LINE09'
        write(lun,'(a)')'LINE10'
        write(lun,'(a)')'LINE11'
c-----
c       put model specifically for 1-D flat isotropic
c-----
c-----
c       put comment line
c-----
        write(lun,'(a)')cmnt1(1:110)
        write(lun,'(a)')cmnt2(1:110)
        curdep = 0.0
        
        do 1000 j=1,mmax
            curdep = curdep + abs(d(j))
            if(curdep .le. refdep)d(j) = - d(j)
            vpv = sqrt(TC(j)/TRho(j))
            vph = sqrt(TA(j)/TRho(j))
            vsv = sqrt(TL(j)/TRho(j))
            vsh = sqrt(TN(j)/TRho(j))
            vpf = sqrt(TF(j)/TRho(j))
            write(lun,'(f11.4,3f9.4,6g11.3)')d(j),vpv,vsv,
     1          TRho(j),qa(j),qb(j),etap(j),etas(j),
     2          frefp(j),frefs(j)
            write(lun,'(11x,3f9.4)')vph,vsh,vpf
 1000   continue
        if(lun.ne.LOT)close (lun)
        return
        end
