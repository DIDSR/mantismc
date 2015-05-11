CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C    PPPPP   EEEEEE  N    N  EEEEEE  L        OOOO   PPPPP   EEEEEE    C
C    P    P  E       NN   N  E       L       O    O  P    P  E         C
C    P    P  E       N N  N  E       L       O    O  P    P  E         C
C    PPPPP   EEEE    N  N N  EEEE    L       O    O  PPPPP   EEEE      C
C    P       E       N   NN  E       L       O    O  P       E         C
C    P       EEEEEE  N    N  EEEEEE  LLLLLL   OOOO   P       EEEEEE    C
C                                                                      C
C                                                   (version 2005).    C
C                                                                      C
C  Subroutine package for Monte Carlo simulation of coupled electron-  C
C  photon transport in homogeneous media.                              C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C  PENELOPE/PENGEOM (version 2005)                                     C
C  Copyright (c) 2001-2005                                             C
C  Universitat de Barcelona                                            C
C                                                                      C
C  Permission to use, copy, modify, distribute and sell this software  C
C  and its documentation for any purpose is hereby granted without     C
C  fee, provided that the above copyright notice appears in all        C
C  copies and that both that copyright notice and this permission      C
C  notice appear in all supporting documentation. The Universitat de   C
C  Barcelona makes no representations about the suitability of this    C
C  software for any purpose. It is provided "as is" without express    C
C  or implied warranty.                                                C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C  The convention used to name the interaction simulation subroutines is
C  the following:
C  - The first letter indicates the particle (E for electrons, P for
C    positrons, G for photons).
C  - The second and third letters denote the interaction mechanism
C    (EL for elastic, IN for inelastic, SI for inner shell ionization,
C    BR for bremsstrahlung, AN for annihilation, RA for Rayleigh, CO for
C    Compton, PH for photoelectric and PP for pair-production).
C  - The random sampling routines have three-letter names. Auxiliary
C    routines, which perform specific calculations, have longer names,
C    with the fourth and subsequent letters and/or numbers indicating
C    the kind of calculation (TX for total x-section, DX for differen-
C    tial x-section) or action (W for write data on a file, R for read
C    data from a file, I for initialization of simulation algorithm).
C
C  The present subroutines may print warning and error messages in
C  the I/O unit 26. This is the standard output unit in the example
C  main programs PENSLAB, PENCYL and PENMAIN.
C
C  Subroutine PEMATW connects files to the I/O units 3 (input) and 7
C  (oputput). However, this does not conflict with the main program,
C  because PEMATW is not invoked during simulation. It is only used by
C  the program MATERIAL, to generate the material data files.
C
C  *********************************************************************
C                       SUBROUTINE PEINIT
C  *********************************************************************
      SUBROUTINE PEINIT(EMAX,NMAT,IRD,IWR,INFO)
C
C  Input of material data and initialization of simulation routines.
C
C  Each material is defined through the input file (unit=IRD), which is
C  created by the program 'material' using information contained in the
C  database. This file can be modified by the user if more accurate in-
C  teraction data are available. Data files for different materials must
C  be concatenated in a single input file, the M-th material in this
C  file is identified by the index M.
C
C  Input arguments:
C    EMAX ... maximum particle energy (kinetic energy for electrons and
C             positrons) used in the simulation. Note: Positrons with
C             energy E may produce photons with energy E+1.022E6.
C    NMAT ... number of materials in the geometry.
C    IRD .... input unit.
C    IWD .... output unit.
C    INFO ... determines the amount of information that is written on
C             the output file,
C               INFO=1, minimal (composition data only).
C               INFO=2, medium (same information as in the material
C                 definition data file, useful to check that the struc-
C                 ture of the latter is correct).
C               INFO=3 or larger, full information, including tables of
C                 interaction properties used in the simulation.
C
C  For the preliminary computations, PEINIT needs to know the absorption
C  energies EABS(KPAR,M) and the simulation parameters C1(M), C2(M),
C  WCC(M) and WCR(M). This information is introduced through the named
C  common block /CSIMPA/ that has to be loaded before invoking the
C  PEINIT subroutine.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N)
      CHARACTER*3 LIT
C  ****  Main-PENELOPE common.
      COMMON/TRACK/E,X,Y,Z,U,V,W,WGHT,KPAR,IBODY,M,ILB(5)
C  ****  Simulation parameters.
      PARAMETER (MAXMAT=10)
      COMMON/CSIMPA/EABS(3,MAXMAT),C1(MAXMAT),C2(MAXMAT),WCC(MAXMAT),
     1  WCR(MAXMAT)
      COMMON/CECUTR/ECUTR(MAXMAT)
C
      COMMON/CERSEC/IERSEC
      IERSEC=0
C
C  ****  Lower limit of the energy grid.
C
      EMIN=1.0D35
      DO M=1,NMAT
        EMIN=MIN(EMIN,EABS(1,M),EABS(2,M),EABS(3,M))
      ENDDO
      IF(EMIN.LT.50.0D0) EMIN=50.0D0
C
      WRITE(IWR,2000)
 2000 FORMAT(/1X,34('*'),/1X,'**   PENELOPE  (version 2005)   **',
     1  /1X,34('*'))
      WRITE(IWR,2001) EMIN,EMAX
 2001 FORMAT(/1X,'EMIN =',1P,E11.4,' eV,  EMAX =',E11.4,' eV')
      IF(EMAX.LT.EMIN+10.0D0) STOP 'The energy interval is too narrow.'
      IF(NMAT.GT.MAXMAT) THEN
        WRITE(IWR,2002) NMAT,MAXMAT,NMAT
 2002   FORMAT(/1X,'*** PENELOPE cannot handle ',I2,' different mater',
     1  'ials.'/5X,'Edit the source file and change the parameter ',
     2  'MAXMAT = ',I2,' to MAXMAT = ',I2)
        STOP 'PEINIT. Too many materials.'
      ENDIF
      IF(INFO.GT.2) WRITE(IWR,2102)
 2102 FORMAT(/1X,'NOTE: 1 mtu = 1 g/cm**2')
C
      CALL EGRID(EMIN,EMAX)  ! Defines the simulation energy grid.
      CALL ESI0  ! Initializes electron impact ionization routines.
      CALL PSI0  ! Initializes positron impact ionization routines.
      CALL GPH0  ! Initializes photoelectric routines.
      CALL RELAX0  ! Initializes atomic relaxation routines.
C
      DO M=1,NMAT
        IF(M.EQ.1) LIT='st'
        IF(M.EQ.2) LIT='nd'
        IF(M.EQ.3) LIT='rd'
        IF(M.GT.3) LIT='th'
        WRITE(IWR,2003) M,LIT
 2003   FORMAT(//1X,22('*')/1X,'**  ',I2,A2,' material   **',
     1    /1X,22('*'))
C
C  ****  Energy limits and thresholds.
C
        WRITE(IWR,2004)
 2004   FORMAT(/1X,'*** Simulation parameters:')
        IF(EABS(1,M).LT.50.0D0) THEN
          EABS(1,M)=50.0D0
          WRITE(IWR,2005)
 2005     FORMAT(1X,'*** Warning: electron absorption energy has ',
     1      'been set to 50 eV')
        ENDIF
        WRITE(IWR,2006) EABS(1,M)
 2006   FORMAT(5X,'Electron absorption energy =',1P,E11.4,' eV')
C
        IF(EABS(2,M).LT.50.0D0) THEN
          EABS(2,M)=50.0D0
          WRITE(IWR,2007)
 2007     FORMAT(1X,'*** Warning: photon absorption energy has ',
     1    'been set to 50 eV')
        ENDIF
        WRITE(IWR,2008) EABS(2,M)
 2008   FORMAT(7X,'Photon absorption energy =',1P,E11.4,' eV')
C
        IF(EABS(3,M).LT.50.0D0) THEN
          EABS(3,M)=50.0D0
          WRITE(IWR,2009)
 2009     FORMAT(1X,'*** Warning: positron absorption energy has ',
     1      'been set to 50 eV')
        ENDIF
        WRITE(IWR,2010) EABS(3,M)
 2010   FORMAT(5X,'Positron absorption energy =',1P,E11.4,' eV')
C
        C1(M)=MIN(0.2D0,ABS(C1(M)))
        C2(M)=MIN(0.2D0,ABS(C2(M)))
        WCC(M)=MIN(ABS(WCC(M)),EMAX)
        IF(WCR(M).LT.0.0D0) WRITE(IWR,2011)
 2011   FORMAT(1X,'*** Warning: soft radiative losses are switched off')
        WRITE(IWR,2012) C1(M),C2(M),WCC(M),MAX(WCR(M),10.0D0)
 2012   FORMAT(6X,'C1 =',1P,E11.4,',       C2 =',E11.4,/5X,'WCC =',
     2    E11.4,' eV,   WCR =',E11.4,' eV',/)
C
        ECUTR(M)=MIN(EABS(1,M),EABS(2,M))
        CALL PEMATR(M,IRD,IWR,INFO)
      ENDDO
      RETURN
      END
C  *********************************************************************
C                       SUBROUTINE EGRID
C  *********************************************************************
      SUBROUTINE EGRID(EMIN,EMAX)
C
C  This subroutine sets the energy grid where transport functions are
C  tabulated. The grid is logarithmically spaced and we assume that it
C  is dense enough to permit accurate linear log-log interpolation of
C  the tabulated functions.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N)
C  ****  Energy grid and interpolation constants for the current energy.
      PARAMETER (NEGP=200)
      COMMON/CEGRID/EL,EU,ET(NEGP),DLEMP(NEGP),DLEMP1,DLFC,
     1  XEL,XE,XEK,KE
C
C  ****  Consistency of the interval end-points.
C
      IF(EMIN.LT.50.0D0) EMIN=50.0D0
      IF(EMIN.GT.EMAX-1.0D0) THEN
        WRITE(26,2100) EMIN,EMAX
 2100   FORMAT(/3X,'EMIN =',1P,E11.4,' eV,  EMAX =',E11.4,' eV')
        STOP 'EGRID. The energy interval is too narrow.'
      ENDIF
C
C  ****  Energy grid points.
C
      EL=0.99999D0*EMIN
      EU=1.00001D0*EMAX
      DLFC=LOG(EU/EL)/DBLE(NEGP-1)
      DLEMP1=LOG(EL)
      DLEMP(1)=DLEMP1
      ET(1)=EL
      DO I=2,NEGP
        DLEMP(I)=DLEMP(I-1)+DLFC
        ET(I)=EXP(DLEMP(I))
      ENDDO
      DLFC=1.0D0/DLFC
C
C  NOTE: To determine the interval KE where the energy E is located, we
C  do the following,
C     XEL=LOG(E)
C     XE=1.0D0+(XEL-DLEMP1)*DLFC
C     KE=XE
C     XEK=XE-KE  ! 'fractional' part of XE (used for interpolation).
C
      RETURN
      END
C  *********************************************************************
C                       SUBROUTINE PEMATR
C  *********************************************************************
      SUBROUTINE PEMATR(M,IRD,IWR,INFO)
C
C  This subroutine reads the definition file of material M (unit IRD)
C  and initializes the simulation routines for this material. Informa-
C  tion is written on unit IWR, the amount of written information is
C  determined by the value of INFO.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N)
      CHARACTER*2 LASYMB
      CHARACTER*60 NAME,LNAME
      PARAMETER (A0B=5.291772083D-9)  ! Bohr radius (cm)
      PARAMETER (HREV=27.2113834D0)  ! Hartree energy (eV)
      PARAMETER (AVOG=6.02214199D23)  ! Avogadro's number
      PARAMETER (REV=5.10998902D5)  ! Electron rest energy (eV)
      PARAMETER (SL=137.03599976D0)  ! Speed of light (1/alpha)
      PARAMETER (PI=3.1415926535897932D0, FOURPI=4.0D0*PI)
C  ****  Composition data.
      PARAMETER (MAXMAT=10)
      COMMON/COMPOS/STF(MAXMAT,30),ZT(MAXMAT),AT(MAXMAT),RHO(MAXMAT),
     1  VMOL(MAXMAT),IZ(MAXMAT,30),NELEM(MAXMAT)
C  ****  Element data.
      COMMON/CADATA/ATW(99),EPX(99),RA1(99),RA2(99),RA3(99),RA4(99),
     1  RA5(99),RSCR(99),ETA(99),EB(99,30),IFI(99,30),IKS(99,30),
     2  NSHT(99),LASYMB(99)
C  ****  Simulation parameters.
      COMMON/CSIMPA/EABS(3,MAXMAT),C1(MAXMAT),C2(MAXMAT),WCC(MAXMAT),
     1  WCR(MAXMAT)
      COMMON/CECUTR/ECUTR(MAXMAT)
C  ****  Energy grid and interpolation constants for the current energy.
      PARAMETER (NEGP=200)
      COMMON/CEGRID/EL,EU,ET(NEGP),DLEMP(NEGP),DLEMP1,DLFC,
     1  XEL,XE,XEK,KE
      COMMON/CRANGE/RANGE(3,MAXMAT,NEGP),RANGEL(3,MAXMAT,NEGP)
C  ****  E/P inelastic collisions.
      PARAMETER (NO=64)
      COMMON/CEIN/EXPOT(MAXMAT),OP2(MAXMAT),F(MAXMAT,NO),UI(MAXMAT,NO),
     1  WRI(MAXMAT,NO),KZ(MAXMAT,NO),KS(MAXMAT,NO),NOSC(MAXMAT)
      COMMON/CEINAC/EINAC(MAXMAT,NEGP,NO),PINAC(MAXMAT,NEGP,NO)
      COMMON/CEINTF/T1EI(NEGP),T2EI(NEGP),T1PI(NEGP),T2PI(NEGP)
C  ****  Partial cross sections of individual shells/oscillators.
      COMMON/CEIN00/SXH0(NO),SXH1(NO),SXH2(NO),SXS0(NO),SXS1(NO),
     1              SXS2(NO),SXT0(NO),SXT1(NO),SXT2(NO)
      COMMON/CPIN00/SYH0(NO),SYH1(NO),SYH2(NO),SYS0(NO),SYS1(NO),
     1              SYS2(NO),SYT0(NO),SYT1(NO),SYT2(NO)
C  ****  Compton scattering.
      PARAMETER (NOCO=64)
      COMMON/CGCO/FCO(MAXMAT,NOCO),UICO(MAXMAT,NOCO),FJ0(MAXMAT,NOCO),
     2  KZCO(MAXMAT,NOCO),KSCO(MAXMAT,NOCO),NOSCCO(MAXMAT)
C  ****  Electron simulation tables.
      COMMON/CEIMFP/SEHEL(MAXMAT,NEGP),SEHIN(MAXMAT,NEGP),
     1  SEISI(MAXMAT,NEGP),SEHBR(MAXMAT,NEGP),SEAUX(MAXMAT,NEGP),
     2  SETOT(MAXMAT,NEGP),CSTPE(MAXMAT,NEGP),RSTPE(MAXMAT,NEGP),
     3  DEL(MAXMAT,NEGP),W1E(MAXMAT,NEGP),W2E(MAXMAT,NEGP),
     4  RNDCE(MAXMAT,NEGP),AE(MAXMAT,NEGP),BE(MAXMAT,NEGP),
     5  T1E(MAXMAT,NEGP),T2E(MAXMAT,NEGP)
C  ****  Positron simulation tables.
      COMMON/CPIMFP/SPHEL(MAXMAT,NEGP),SPHIN(MAXMAT,NEGP),
     1  SPISI(MAXMAT,NEGP),SPHBR(MAXMAT,NEGP),SPAN(MAXMAT,NEGP),
     2  SPAUX(MAXMAT,NEGP),SPTOT(MAXMAT,NEGP),CSTPP(MAXMAT,NEGP),
     3  RSTPP(MAXMAT,NEGP),W1P(MAXMAT,NEGP),W2P(MAXMAT,NEGP),
     4  RNDCP(MAXMAT,NEGP),AP(MAXMAT,NEGP),BP(MAXMAT,NEGP),
     5  T1P(MAXMAT,NEGP),T2P(MAXMAT,NEGP)
C  ****  Photon simulation tables.
      COMMON/CGIMFP/SGRA(MAXMAT,NEGP),SGCO(MAXMAT,NEGP),
     1  SGPH(MAXMAT,NEGP),SGPP(MAXMAT,NEGP),SGAUX(MAXMAT,NEGP)
      PARAMETER (NDIM=1500)
      COMMON/CGPH01/ER(NDIM),XSR(NDIM),NPHD
C  ****  Auxiliary arrays.
      DIMENSION EIT(NEGP),EITL(NEGP),FL(NEGP),F1(NEGP),F2(NEGP),
     1  F3(NEGP),F4(NEGP),A(NEGP),B(NEGP),C(NEGP),D(NEGP),RADY(NEGP)
C  ****  Inner shell ionization by electron and positron impact.
      PARAMETER (NRP=6000)
      COMMON/CESI0/XESI(NRP,9),IESIF(99),IESIL(99),NSESI(99),NCURE
      COMMON/CPSI0/XPSI(NRP,9),IPSIF(99),IPSIL(99),NSPSI(99),NCURP
C
C  ****  Rescaling of calculated MFPs.
C  When ISCALE is set equal to 1, the program rescales the calculated
C  total cross sections and MFPs to reproduce the cross sections and
C  stopping powers read from the input material data file. To avoid this
C  rescaling, set ISCALE=0.
C
      ISCALE=1
C
C
C  ************  Material characteristics.
C
      IF(M.GT.MAXMAT) STOP 'PEMATR. Too many materials.'
C
      LNAME=' PENELOPE (v. 2005)  Material data file ...............'
      READ(IRD,'(A55)') NAME
      IF(NAME.NE.LNAME) THEN
        WRITE(IWR,'(/1X,''I/O error. Corrupt material data file.'')')
        WRITE(IWR,'(''     The first line is: '',A55)') NAME
        WRITE(IWR,'(''     ... and should be: '',A55)') LNAME
        STOP 'PEMATR. Corrupt material data file.'
      ENDIF
      WRITE(IWR,'(A55)') NAME
C
      READ(IRD,5001) LNAME
 5001 FORMAT(11X,A60)
      WRITE(IWR,1001) LNAME
 1001 FORMAT(' Material: ',A60)
      READ(IRD,5002) RHO(M)
 5002 FORMAT(15X,E15.8)
      WRITE(IWR,1002) RHO(M)
 1002 FORMAT(' Mass density =',1P,E15.8,' g/cm**3')
      READ(IRD,5003) NELEM(M)
 5003 FORMAT(37X,I3)
      WRITE(IWR,1003) NELEM(M)
 1003 FORMAT(' Number of elements in the molecule = ',I2)
      IF(NELEM(M).GT.30) STOP 'PEMATR. Too many elements.'
      ZT(M)=0.0D0
      AT(M)=0.0D0
      DO I=1,NELEM(M)
        READ(IRD,5004) IZ(M,I),STF(M,I)
 5004   FORMAT(18X,I3,19X,E15.8)
        WRITE(IWR,1004) LASYMB(IZ(M,I)),IZ(M,I),STF(M,I)
 1004   FORMAT(3X,' Element: ',A2,' (Z=',I2,'), atoms/molecule =',1P,
     1    E15.8)
        ZT(M)=ZT(M)+STF(M,I)*IZ(M,I)
        IZZ=IZ(M,I)
        AT(M)=AT(M)+ATW(IZZ)*STF(M,I)
      ENDDO
      VMOL(M)=AVOG*RHO(M)/AT(M)
      IF(INFO.GE.2) WRITE(IWR,
     1  '(/1X,''Molecular density = '',1P,E15.8,'' 1/cm**3'')') VMOL(M)
      OP2(M)=FOURPI*ZT(M)*VMOL(M)*A0B**3*HREV**2
      OMEGA=SQRT(OP2(M))
C
      IF(INFO.GE.2) THEN
        WRITE(IWR,'(/1X,''*** Electron/positron inelastic'',
     1    '' scattering.'')')
        WRITE(IWR,
     1  '(1X,''Plasma energy = '',1P,E15.8,'' eV'')') OMEGA
      ENDIF
C
      READ(IRD,5005) EXPOT(M)
 5005 FORMAT(25X,E15.8)
      WRITE(IWR,1005) EXPOT(M)
 1005 FORMAT(' Mean excitation energy =',1P,E15.8,' eV')
C
C  ****  E/P inelastic collisions.
C
      READ(IRD,5006) NOSC(M)
 5006 FORMAT(24X,I3)
      IF(INFO.GE.2.OR.NOSC(M).GT.NO) WRITE(IWR,1006) NOSC(M)
 1006 FORMAT(' Number of oscillators =',I3)
      IF(NOSC(M).GT.NO) STOP 'PEMATR. Too many oscillators.'
      IF(INFO.GE.2) WRITE(IWR,5507)
 5507 FORMAT(/11X,'Fi',12X,'Ui (eV)',9X,'Wi (eV)',6X,
     1  'KZ  KS',/1X,60('-'))
      EXPT=0.0D0
      DO I=1,NOSC(M)
        READ(IRD,5007) F(M,I),UI(M,I),WRI(M,I),KZ(M,I),KS(M,I)
 5007   FORMAT(4X,1P,3E16.8,2I4)
        IF(INFO.GE.2) WRITE(IWR,1007) I,F(M,I),UI(M,I),WRI(M,I),
     1   KZ(M,I),KS(M,I)
 1007   FORMAT(I4,1P,3E16.8,2I4)
        EXPT=EXPT+F(M,I)*LOG(WRI(M,I))
      ENDDO
      EXPT=EXP(EXPT/ZT(M))
C
      IF(ABS(EXPT-EXPOT(M)).GT.1.0D-8*EXPOT(M)) THEN
        WRITE(26,*) 'EXPT      =',EXPT
        WRITE(26,*) 'EXPOT (M) =',EXPOT(M)
        WRITE(26,*) 'Inconsistent oscillator data.'
        STOP 'PEMATR. Inconsistent oscillator data.'
      ENDIF
C
C  ****  Compton scattering.
C
      READ(IRD,5106) NOSCCO(M)
 5106 FORMAT(19X,I3)
      IF(INFO.GE.2.OR.NOSCCO(M).GT.NOCO) THEN
        WRITE(IWR,'(/1X,''*** Compton scattering '',
     1    ''(Impulse Approximation).'')')
        WRITE(IWR,1106) NOSCCO(M)
 1106   FORMAT(1X,'Number of shells =',I3)
        IF(NOSCCO(M).GT.NOCO) STOP 'PEMATR. Too many shells.'
      ENDIF
      IF(INFO.GE.2) WRITE(IWR,5607)
 5607 FORMAT(/11X,'Fi',12X,'Ui (eV)',9X,'Ji(0)',8X,'KZ  KS',
     1       /1X,60('-'))
      DO I=1,NOSCCO(M)
        READ(IRD,5107) FCO(M,I),UICO(M,I),FJ0(M,I),KZCO(M,I),KSCO(M,I)
 5107   FORMAT(4X,1P,3E16.8,2I4)
        IF(INFO.GE.2) WRITE(IWR,1207) I,FCO(M,I),UICO(M,I),FJ0(M,I),
     1    KZCO(M,I),KSCO(M,I)
 1207   FORMAT(I4,1P,3E16.8,2I4)
        FJ0(M,I)=FJ0(M,I)*SL
      ENDDO
C
C  ************  Atomic relaxation data.
C
      DO I=1,NELEM(M)
        CALL RELAXR(IRD,IWR,INFO)
      ENDDO
C
C  ****  Electron and positron interaction properties.
C
C  ****  Scaled bremss x-section.
      IF(WCR(M).GE.0.0D0) THEN
        WCR(M)=MAX(WCR(M),10.0D0)
        WCRM=WCR(M)
        IBREMS=1
      ELSE
        WCRM=10.0D0
        WCR(M)=10.0D0
        IBREMS=0
      ENDIF
      CALL EBRR(WCRM,M,IRD,IWR,INFO)
C  ****  Bremss angular distribution.
      CALL BRAR(M,IRD,IWR,INFO)
C
C  ****  Stopping powers.
C
      READ(IRD,5008) NDATA
 5008 FORMAT(58X,I4)
      IF(INFO.GE.2) WRITE(IWR,1008) NDATA
 1008 FORMAT(/1X,'*** Stopping powers for electrons and positrons',
     1  ',  NDATA =',I4)
      IF(NDATA.GT.NEGP) STOP 'PEMATR. Too many data points (1).'
      IF(INFO.GE.2) WRITE(IWR,1108)
 1108 FORMAT(/2X,'Energy',5X,'Scol,e-',5X,'Srad,e-',5X,'Scol,e+',
     1  5X,'Srad,e+',/3X,'(eV)',5X,'(MeV/mtu)',3X,'(MeV/mtu)',3X,
     2  '(MeV/mtu)',3X,'(MeV/mtu)',/1X,58('-'))
      DO I=1,NDATA
        READ(IRD,*) EIT(I),F1(I),F2(I),F3(I),F4(I)
        IF(INFO.GE.2) WRITE(IWR,'(1P,E10.3,5E12.5)')
     1    EIT(I),F1(I),F2(I),F3(I),F4(I)
        EITL(I)=LOG(EIT(I))
      ENDDO
C
C  ****  Inelastic x-section for electrons.
C
      WCCM=WCC(M)
      DIFMAX=0.0D0
      DO I=1,NDATA
        IF(EIT(I).GE.EL.AND.EIT(I).LE.EU) THEN
          STPI=F1(I)*RHO(M)*1.0D6  ! Collision stopping power.
          CALL EINTX(EIT(I),WCCM,XH0,XH1,XH2,XS0,XS1,XS2,XT1,XT2,
     1     DELTA,M)
          STPC=(XS1+XH1)*VMOL(M)
          DIFMAX=MAX(DIFMAX,ABS(STPC-STPI)/(0.01D0*STPI))
        ENDIF
      ENDDO
C
      IF(DIFMAX.GT.1.0D-3.AND.ISCALE.EQ.1) THEN
        ICOR=1
        DO I=1,NDATA
          FL(I)=LOG(F1(I)*RHO(M)*1.0D6)
        ENDDO
        CALL SPLINE(EITL,FL,A,B,C,D,0.0D0,0.0D0,NDATA)
      ELSE
        ICOR=0
      ENDIF
C
      DO I=1,NEGP
        CALL EINTX(ET(I),WCCM,XH0,XH1,XH2,XS0,XS1,XS2,XT1,XT2,DELTA,M)
        STPC=(XS1+XH1)*VMOL(M)
        IF(ICOR.EQ.1) THEN
          EC=DLEMP(I)
          CALL FINDI(EITL,EC,NDATA,J)
          STPI=EXP(A(J)+EC*(B(J)+EC*(C(J)+EC*D(J))))
          FACT=STPI/STPC
        ELSE
          FACT=1.0D0
        ENDIF
        CSTPE(M,I)=STPC*FACT
        SEHIN(M,I)=LOG(MAX(XH0*VMOL(M)*FACT,1.0D-35))
        W1E(M,I)=XS1*VMOL(M)*FACT
        W2E(M,I)=XS2*VMOL(M)*FACT
        T1EI(I)=XT1*VMOL(M)*FACT
        T2EI(I)=XT2*VMOL(M)*FACT
        DEL(M,I)=DELTA
C
        SXHT=0.0D0
        DO KO=1,NOSC(M)
          EINAC(M,I,KO)=SXHT
          SXHT=SXHT+SXH0(KO)
        ENDDO
        IF(SXHT.GT.1.0D-35) THEN
          FNORM=1.0D0/SXHT
          DO KO=1,NOSC(M)
            EINAC(M,I,KO)=EINAC(M,I,KO)*FNORM
          ENDDO
        ELSE
          DO KO=1,NOSC(M)
            EINAC(M,I,KO)=1.0D0
          ENDDO
        ENDIF
      ENDDO
C
C  ****  Bremss x-section for electrons.
C
      DIFMAX=0.0D0
      DO I=1,NDATA
        IF(EIT(I).GE.EL.AND.EIT(I).LT.EU) THEN
          STPI=F2(I)*RHO(M)*1.0D6  ! Radiative stopping power.
          CALL EBRTX(EIT(I),WCRM,XH0,XH1,XH2,XS1,XS2)
          STPR=(XS1+XH1)*VMOL(M)
          DIFMAX=MAX(DIFMAX,ABS(STPR-STPI)/(0.01D0*STPI))
        ENDIF
      ENDDO
C
      IF(DIFMAX.GT.1.0D-3.AND.ISCALE.EQ.1) THEN
        ICOR=1
        DO I=1,NDATA
          FL(I)=LOG(F2(I)*RHO(M)*1.0D6)
        ENDDO
        CALL SPLINE(EITL,FL,A,B,C,D,0.0D0,0.0D0,NDATA)
      ELSE
        ICOR=0
      ENDIF
C
      DO I=1,NEGP
        CALL EBRTX(ET(I),WCRM,XH0,XH1,XH2,XS1,XS2)
        STPR=(XS1+XH1)*VMOL(M)
        IF(ICOR.EQ.1) THEN
          EC=DLEMP(I)
          CALL FINDI(EITL,EC,NDATA,J)
          STPI=EXP(A(J)+EC*(B(J)+EC*(C(J)+EC*D(J))))
          FACT=STPI/STPR
        ELSE
          FACT=1.0D0
        ENDIF
        RSTPE(M,I)=STPR*FACT
        SEHBR(M,I)=LOG(MAX(XH0*VMOL(M)*FACT,1.0D-35))
        IF(IBREMS.EQ.1) THEN
          W1E(M,I)=W1E(M,I)+XS1*VMOL(M)*FACT
          W2E(M,I)=W2E(M,I)+XS2*VMOL(M)*FACT
        ENDIF
      ENDDO
C
C  ****  Electron range as a function of energy.
C
      DO I=1,NEGP
        F1(I)=1.0D0/(CSTPE(M,I)+RSTPE(M,I))
      ENDDO
      CALL SPLINE(ET,F1,A,B,C,D,0.0D0,0.0D0,NEGP)
      RANGE(1,M,1)=1.0D-8
      RANGEL(1,M,1)=LOG(RANGE(1,M,1))
      DO I=2,NEGP
        XL=ET(I-1)
        XU=ET(I)
        CALL SINTEG(ET,A,B,C,D,XL,XU,DR,NEGP)
        RANGE(1,M,I)=RANGE(1,M,I-1)+DR
        RANGEL(1,M,I)=LOG(RANGE(1,M,I))
      ENDDO
C
C  ****  Electron radiation yield as a function of energy.
C
      DO I=1,NEGP
        F1(I)=RSTPE(M,I)/(CSTPE(M,I)+RSTPE(M,I))
      ENDDO
      CALL SPLINE(ET,F1,A,B,C,D,0.0D0,0.0D0,NEGP)
      RADY(1)=0.0D0
      DO I=2,NEGP
        XL=ET(I-1)
        XU=ET(I)
        CALL SINTEG(ET,A,B,C,D,XL,XU,DR,NEGP)
        RADY(I)=RADY(I-1)+DR
      ENDDO
C
C  ****  Print electron stopping power tables.
C
      IF(INFO.GE.3) THEN
        WRITE(IWR,1009)
 1009   FORMAT(/1X,'PENELOPE >>>  Stopping powers for electrons')
        WRITE(IWR,1010)
 1010   FORMAT(/3X,'Energy',8X,'Scol',9X,'Srad',9X,'range',5X,
     1    'Rad. Yield',6X,'delta',/4X,'(eV)',7X,'(eV/mtu)',5X,
     2    '(eV/mtu)',7X,'(mtu)',/1X,77('-'))
        DO I=1,NEGP
          WRITE(IWR,'(1P,6(E12.5,1X))') ET(I),CSTPE(M,I)/RHO(M),
     1      RSTPE(M,I)/RHO(M),RANGE(1,M,I)*RHO(M),RADY(I)/ET(I),
     2      DEL(M,I)
        ENDDO
      ENDIF
C
C  ****  Inelastic x-section for positrons.
C
      DIFMAX=0.0D0
      DO I=1,NDATA
        IF(EIT(I).GE.EL.AND.EIT(I).LE.EU) THEN
          STPI=F3(I)*RHO(M)*1.0D6  ! Collision stopping power.
          CALL PINTX(EIT(I),WCCM,XH0,XH1,XH2,XS0,XS1,XS2,XT1,XT2,
     1      DELTA,M)
          STPC=(XS1+XH1)*VMOL(M)
          DIFMAX=MAX(DIFMAX,ABS(STPC-STPI)/(0.01D0*STPI))
        ENDIF
      ENDDO
C
      IF(DIFMAX.GT.1.0D-3.AND.ISCALE.EQ.1) THEN
        ICOR=1
        DO I=1,NDATA
          FL(I)=LOG(F3(I)*RHO(M)*1.0D6)
        ENDDO
        CALL SPLINE(EITL,FL,A,B,C,D,0.0D0,0.0D0,NDATA)
      ELSE
        ICOR=0
      ENDIF
C
      DO I=1,NEGP
        CALL PINTX(ET(I),WCCM,XH0,XH1,XH2,XS0,XS1,XS2,XT1,XT2,DELTA,M)
        STPC=(XS1+XH1)*VMOL(M)
        IF(ICOR.EQ.1) THEN
          EC=DLEMP(I)
          CALL FINDI(EITL,EC,NDATA,J)
          STPI=EXP(A(J)+EC*(B(J)+EC*(C(J)+EC*D(J))))
          FACT=STPI/STPC
        ELSE
          FACT=1.0D0
        ENDIF
        CSTPP(M,I)=STPC*FACT
        SPHIN(M,I)=LOG(MAX(XH0*VMOL(M)*FACT,1.0D-35))
        W1P(M,I)=XS1*VMOL(M)*FACT
        W2P(M,I)=XS2*VMOL(M)*FACT
        T1PI(I)=XT1*VMOL(M)*FACT
        T2PI(I)=XT2*VMOL(M)*FACT
C
        SXHT=0.0D0
        DO KO=1,NOSC(M)
          PINAC(M,I,KO)=SXHT
          SXHT=SXHT+SYH0(KO)
        ENDDO
        IF(SXHT.GT.1.0D-35) THEN
          FNORM=1.0D0/SXHT
          DO KO=1,NOSC(M)
            PINAC(M,I,KO)=PINAC(M,I,KO)*FNORM
          ENDDO
        ELSE
          DO KO=1,NOSC(M)
            PINAC(M,I,KO)=1.0D0
          ENDDO
        ENDIF
      ENDDO
C
C  ****  Bremss x-section for positrons.
C
      DIFMAX=0.0D0
      DO I=1,NDATA
        IF(EIT(I).GE.EL.AND.EIT(I).LT.EU) THEN
          STPI=F4(I)*RHO(M)*1.0D6  ! Radiative stopping power.
          CALL PBRTX(EIT(I),WCRM,XH0,XH1,XH2,XS1,XS2)
          STPR=(XS1+XH1)*VMOL(M)
          DIFMAX=MAX(DIFMAX,ABS(STPR-STPI)/(0.01D0*STPI))
        ENDIF
      ENDDO
C
      IF(DIFMAX.GT.1.0D-3.AND.ISCALE.EQ.1) THEN
        ICOR=1
        DO I=1,NDATA
          FL(I)=LOG(F4(I)*RHO(M)*1.0D6)
        ENDDO
        CALL SPLINE(EITL,FL,A,B,C,D,0.0D0,0.0D0,NDATA)
      ELSE
        ICOR=0
      ENDIF
C
      DO I=1,NEGP
        CALL PBRTX(ET(I),WCRM,XH0,XH1,XH2,XS1,XS2)
        STPR=(XS1+XH1)*VMOL(M)
        IF(ICOR.EQ.1) THEN
          EC=DLEMP(I)
          CALL FINDI(EITL,EC,NDATA,J)
          STPI=EXP(A(J)+EC*(B(J)+EC*(C(J)+EC*D(J))))
          FACT=STPI/STPR
        ELSE
          FACT=1.0D0
        ENDIF
        RSTPP(M,I)=STPR*FACT
        SPHBR(M,I)=LOG(MAX(XH0*VMOL(M)*FACT,1.0D-35))
        IF(IBREMS.EQ.1) THEN
          W1P(M,I)=W1P(M,I)+XS1*VMOL(M)*FACT
          W2P(M,I)=W2P(M,I)+XS2*VMOL(M)*FACT
        ENDIF
      ENDDO
C
C  ****  Positron annihilation.
C
      DO I=1,NEGP
        CALL PANTX(ET(I),TXAN)
        SPAN(M,I)=LOG(TXAN*ZT(M)*VMOL(M))
      ENDDO
C
C  ****  Positron range as a function of energy.
C
      DO I=1,NEGP
        F3(I)=1.0D0/(CSTPP(M,I)+RSTPP(M,I))
      ENDDO
      CALL SPLINE(ET,F3,A,B,C,D,0.0D0,0.0D0,NEGP)
      RANGE(3,M,1)=1.0D-8
      RANGEL(3,M,1)=LOG(RANGE(3,M,1))
      DO I=2,NEGP
        XL=ET(I-1)
        XU=ET(I)
        CALL SINTEG(ET,A,B,C,D,XL,XU,DR,NEGP)
        RANGE(3,M,I)=RANGE(3,M,I-1)+DR
        RANGEL(3,M,I)=LOG(RANGE(3,M,I))
      ENDDO
C
C  ****  Positron radiation yield as a function of energy.
C
      DO I=1,NEGP
        F3(I)=RSTPP(M,I)/(CSTPP(M,I)+RSTPP(M,I))
      ENDDO
      CALL SPLINE(ET,F3,A,B,C,D,0.0D0,0.0D0,NEGP)
      RADY(1)=0.0D0
      DO I=2,NEGP
        XL=ET(I-1)
        XU=ET(I)
        CALL SINTEG(ET,A,B,C,D,XL,XU,DR,NEGP)
        RADY(I)=RADY(I-1)+DR
      ENDDO
C
C  ****  Print positron stopping power tables.
C
      IF(INFO.GE.3) THEN
        WRITE(IWR,1011)
 1011   FORMAT(/1X,'PENELOPE >>>  Stopping powers for positrons')
        WRITE(IWR,1012)
 1012   FORMAT(/3X,'Energy',8X,'Scol',9X,'Srad',9X,'range',5X,
     1    'Rad. Yield',3X,'annih. mfp.',/4X,'(eV)',7X,'(eV/mtu)',5X,
     2    '(eV/mtu)',7X,'(mtu)',21X,'(mtu)',/1X,77('-'))
        DO I=1,NEGP
          WRITE(IWR,'(1P,6(E12.5,1X))') ET(I),CSTPP(M,I)/RHO(M),
     1     RSTPP(M,I)/RHO(M),RANGE(3,M,I)*RHO(M),RADY(I)/ET(I),
     2     RHO(M)/EXP(SPAN(M,I))
        ENDDO
      ENDIF
C
      IF(INFO.GE.3) WRITE(IWR,'(/1X,''PENELOPE >>>  Soft stopping '',
     1  ''power and energy straggling'')')
      IF(INFO.GE.3) WRITE(IWR,1111)
 1111 FORMAT(/3X,'Energy',7X,'SStp,e-',6X,'SStr,e-',5X,'STP,e-',
     1  5X,'SStp,e+',6X,'SStr,e+',5X,'Stp,e+',/4X,'(eV)',7X,
     2  '(eV/mtu)',4X,'(eV**2/mtu)',2X,'soft/tot',3X,'(eV/mtu)',4X,
     3  '(eV**2/mtu)',2X,'soft/tot',/1X,84('-'))
      FMTU=1.0D0/RHO(M)
      DO I=1,NEGP
C  ****  Soft energy-loss events are switched off when W1 is too small.
        FSOFTE=W1E(M,I)/(CSTPE(M,I)+RSTPE(M,I))
        FSOFTP=W1P(M,I)/(CSTPP(M,I)+RSTPP(M,I))
        IF(W1E(M,I).GT.1.0D-5*(CSTPE(M,I)+RSTPE(M,I))) THEN
          W1EW=W1E(M,I)
          W2EW=W2E(M,I)
          W1E(M,I)=LOG(MAX(W1E(M,I),1.0D-35))
          W2E(M,I)=LOG(MAX(W2E(M,I),1.0D-35))
        ELSE
          W1EW=0.0D0
          W2EW=0.0D0
          W1E(M,I)=-100.0D0
          W2E(M,I)=-100.0D0
        ENDIF
        IF(W1P(M,I).GT.1.0D-5*(CSTPP(M,I)+RSTPP(M,I))) THEN
          W1PW=W1P(M,I)
          W2PW=W2P(M,I)
          W1P(M,I)=LOG(MAX(W1P(M,I),1.0D-35))
          W2P(M,I)=LOG(MAX(W2P(M,I),1.0D-35))
        ELSE
          W1PW=0.0D0
          W2PW=0.0D0
          W1P(M,I)=-100.0D0
          W2P(M,I)=-100.0D0
        ENDIF
        IF(INFO.GE.3) WRITE(IWR,'(1P,E12.5,2(2E13.5,E10.2))')
     1    ET(I),W1EW*FMTU,W2EW*FMTU,FSOFTE,W1PW*FMTU,W2PW*FMTU,FSOFTP
      ENDDO
C
C  ****  Elastic scattering of electrons and positrons.
C
      CALL EELR(M,IRD,IWR,INFO)
      CALL EELdR(M,IRD,IWR,INFO)  ! Uses the ELSEP database.
C
      IF(INFO.GE.3) THEN
        WRITE(IWR,'(/1X,''PENELOPE >>>  Soft angular transport coef'',
     1    ''ficients'')')
        WRITE(IWR,1113)
 1113   FORMAT(/3X,'Energy',6X,'SITMFP1,e-',3X,'SITMFP2,e-',3X,
     1    'SITMFP1,e+',3X,'SITMFP2,e+',/4X,'(eV)',8X,'(1/mtu)',6X,
     2    '(1/mtu)',6X,'(1/mtu)',6X,'(1/mtu)',/1X,64('-'))
        DO I=1,NEGP
          WRITE(IWR,'(1P,E12.5,4E13.5)')
     1    ET(I),EXP(T1E(M,I))*FMTU,EXP(T2E(M,I))*FMTU,
     2    EXP(T1P(M,I))*FMTU,EXP(T2P(M,I))*FMTU
        ENDDO
      ENDIF
C
C  ****  Inner shell ionization by electron and positron impact.
C
      CALL ESIR(M,IRD,IWR,INFO)
      CALL PSIR(M,IRD,IWR,INFO)
C
      DO I=1,NEGP
        SEIST=0.0D0
        DO J=1,NELEM(M)
          IZZ=IZ(M,J)
          INDC=IESIF(IZZ)-1
          DO ISH=1,NSESI(IZZ)
            WCUT=EB(IZZ,ISH)
            IF(WCUT.GT.ECUTR(M).AND.WCUT.LT.ET(I)) THEN
              PCSI=EXP(XESI(INDC+I,ISH))
              IF(PCSI.GT.1.1D-35) SEIST=SEIST+PCSI*STF(M,J)
            ENDIF
          ENDDO
        ENDDO
        SEISI(M,I)=LOG(MAX(SEIST*VMOL(M),1.0D-35))
      ENDDO
C
      DO I=1,NEGP
        SPIST=0.0D0
        DO J=1,NELEM(M)
          IZZ=IZ(M,J)
          INDC=IPSIF(IZZ)-1
          DO ISH=1,NSPSI(IZZ)
            WCUT=EB(IZZ,ISH)
            IF(WCUT.GT.ECUTR(M).AND.WCUT.LT.ET(I)) THEN
              PCSI=EXP(XPSI(INDC+I,ISH))
              IF(PCSI.GT.1.1D-35) SPIST=SPIST+PCSI*STF(M,J)
            ENDIF
          ENDDO
        ENDDO
        SPISI(M,I)=LOG(MAX(SPIST*VMOL(M),1.0D-35))
      ENDDO
C
C  ****  Print electron mean free path tables.
C
      IF(INFO.GE.3) THEN
        WRITE(IWR,1109)
 1109   FORMAT(/1X,'PENELOPE >>>  Electron mean free paths (hard',
     1    ' events)',/1X,'*** NOTE: The MFP for inner shell ionization',
     2    ' (isi) is listed only for',/11X,'completeness. The MFP ',
     3    'for inelastic collisions (in) has been',/11X,'calculated by',
     4    ' considering all inelastic events, including isi.')
        WRITE(IWR,1110)
 1110   FORMAT(/3X,'Energy',8X,'MFPel',8X,'MFPin',8X,'MFPbr',7X,
     1   'MFPtot',7X,'MFPisi',/4X,'(eV)',9X,'(mtu)',8X,'(mtu)',8X,
     2   '(mtu)',8X,'(mtu)',8X,'(mtu)',/1X,77('-'))
      ENDIF
      DO I=1,NEGP
        FPEL=EXP(SEHEL(M,I))
        FPIN=EXP(SEHIN(M,I))
        FPSI=EXP(SEISI(M,I))
        FPBR=EXP(SEHBR(M,I))
        FPTOT=FPEL+FPIN+FPBR
        IF(INFO.GE.3) WRITE(IWR,'(1P,6(E12.5,1X))') ET(I),RHO(M)/FPEL,
     1    RHO(M)/FPIN,RHO(M)/FPBR,RHO(M)/FPTOT,RHO(M)/FPSI
        SEAUX(M,I)=LOG(1.0D-35)
        SETOT(M,I)=LOG(FPTOT+FPSI)
      ENDDO
C
      DO I=2,NEGP-1
        IF(EXP(SETOT(M,I)).GT.1.005D0*EXP(SETOT(M,I-1)).AND.
     1     EXP(SETOT(M,I)).GT.1.005D0*EXP(SETOT(M,I+1)).AND.
     2     ET(I).GT.EABS(1,M).AND.ET(I).LT.1.0D6) THEN
          WRITE(IWR,1112) ET(I)
          WRITE(26,1112) ET(I)
 1112     FORMAT(/1X,'WARNING: The electron hard IMFP has a maximum',
     1      ' at E =',1P,E13.6,' eV')
        ENDIF
      ENDDO
C
C  ****  Print positron mean free path tables.
C
      IF(INFO.GE.3) THEN
        WRITE(IWR,1209)
 1209   FORMAT(/1X,'PENELOPE >>>  Positron mean free paths (hard',
     1    ' events)',/1X,'*** NOTE: The MFP for inner shell ionization',
     2    ' (isi) is listed only for',/11X,'completeness. The MFP ',
     3    'for inelastic collisions (in) has been',/11X,'calculated by',
     4    ' considering all inelastic events, including isi.')
        WRITE(IWR,1210)
 1210   FORMAT(/3X,'Energy',8X,'MFPel',8X,'MFPin',8X,'MFPbr',7X,
     1   'MFPan',8X,'MFPtot',7X,'MFPisi',/4X,'(eV)',9X,'(mtu)',8X,
     2   '(mtu)',8X,'(mtu)',7X,'(mtu)',9X,'(mtu)',8X,'(mtu)',
     3   /1X,90('-'))
      ENDIF
      DO I=1,NEGP
        FPEL=EXP(SPHEL(M,I))
        FPIN=EXP(SPHIN(M,I))
        FPSI=EXP(SPISI(M,I))
        FPBR=EXP(SPHBR(M,I))
        FPAN=EXP(SPAN(M,I))
        FPTOT=FPEL+FPIN+FPBR+FPAN
        IF(INFO.GE.3) WRITE(IWR,'(1P,7(E12.5,1X))') ET(I),RHO(M)/FPEL,
     1    RHO(M)/FPIN,RHO(M)/FPBR,RHO(M)/FPAN,RHO(M)/FPTOT,RHO(M)/FPSI
        SPAUX(M,I)=LOG(1.0D-35)
        SPTOT(M,I)=LOG(FPTOT+FPSI)
      ENDDO
C
      DO I=2,NEGP-1
        IF(EXP(SPTOT(M,I)).GT.1.005D0*EXP(SPTOT(M,I-1)).AND.
     1     EXP(SPTOT(M,I)).GT.1.005D0*EXP(SPTOT(M,I+1)).AND.
     2     ET(I).GT.EABS(3,M).AND.ET(I).LT.1.0D6) THEN
          WRITE(IWR,1212) ET(I)
          WRITE(26,1212) ET(I)
 1212     FORMAT(/1X,'WARNING: The positron hard IMFP has a maximum',
     1      ' at E =',1P,E13.6,' eV')
        ENDIF
      ENDDO
C
C  ************  Photon interaction properties.
C
      READ(IRD,5009) NDATA
 5009 FORMAT(67X,I4)
      IF(INFO.GE.2) WRITE(IWR,1013) NDATA
 1013 FORMAT(/1X,'*** Rayleigh, Compton and pair-production cross',
     1  ' sections,  NDATA =',I4)
      IF(NDATA.GT.NEGP) STOP 'PEMATR. Too many data points (2).'
      IF(INFO.GE.2) WRITE(IWR,1114)
 1114 FORMAT(/2X,'Energy',5X,'CS_Rayl',5X,'CS_Comp',5X,'CS_pair',/3X,
     1  '(eV)',6X,'(cm**2)',5X,'(cm**2)',5X,'(cm**2)',/1X,46('-'))
      DO I=1,NDATA
        READ(IRD,*) EIT(I),F1(I),F2(I),F3(I)
        IF(INFO.GE.2) WRITE(IWR,'(1P,E10.3,5E12.5)')
     1    EIT(I),F1(I),F2(I),F3(I)
        EITL(I)=LOG(EIT(I))
      ENDDO
C
C  ****  Rayleigh scattering.
C
      VMOLL=LOG(VMOL(M))
      DO I=1,NDATA
        FL(I)=LOG(F1(I))
      ENDDO
      CALL SPLINE(EITL,FL,A,B,C,D,0.0D0,0.0D0,NDATA)
      DO I=1,NEGP
        EC=DLEMP(I)
        CALL FINDI(EITL,EC,NDATA,J)
        SGRA(M,I)=A(J)+EC*(B(J)+EC*(C(J)+EC*D(J)))+VMOLL
      ENDDO
C
      CALL GRAI(M)
C
C  ****  Compton scattering.
C
      DO I=1,NDATA
        FL(I)=LOG(MAX(F2(I),1.0D-35))
      ENDDO
      CALL SPLINE(EITL,FL,A,B,C,D,0.0D0,0.0D0,NDATA)
      DO I=1,NEGP
        EC=DLEMP(I)
        CALL FINDI(EITL,EC,NDATA,J)
        SGCO(M,I)=A(J)+EC*(B(J)+EC*(C(J)+EC*D(J)))+VMOLL
      ENDDO
C
C  ****  Electron-positron pair production.
C
      NP=0
      DO 1 I=1,NDATA
        IF(EIT(I).LT.1.023D6) GO TO 1
        NP=NP+1
        F4(NP)=EITL(I)
        FL(NP)=LOG(MAX(F3(I),1.0D-35))
    1 CONTINUE
      CALL SPLINE(F4,FL,A,B,C,D,0.0D0,0.0D0,NP)
      DO I=1,NEGP
        IF(ET(I).LT.1.023D6) THEN
          SGPP(M,I)=-80.6D0
        ELSE
          EC=DLEMP(I)
          CALL FINDI(F4,EC,NP,J)
          SGPP(M,I)=A(J)+EC*(B(J)+EC*(C(J)+EC*D(J)))+VMOLL
        ENDIF
      ENDDO
C
C  ****  Photoelectric absorption.
C
      CALL GPHR(M,IRD,IWR,INFO)
C
C  ****  Photon 'range' (= mean free path, slightly underestimated).
C
      DO KE=1,NEGP
        PRA=EXP(SGRA(M,KE))
        PCO=EXP(SGCO(M,KE))
        IF(ET(KE).LT.1.023D6) THEN
          PPP=1.0D-35
        ELSE
          PPP=EXP(SGPP(M,KE))
        ENDIF
        PPH=SGPH(M,KE)
        PT=(PRA+PCO+PPP+PPH)
        RANGE(2,M,KE)=1.0D0/PT
        RANGEL(2,M,KE)=LOG(RANGE(2,M,KE))
      ENDDO
C
      IF(INFO.GE.3) THEN
        WRITE(IWR,1014)
 1014   FORMAT(/1X,'PENELOPE >>>  Photon mass attenuation coefficients')
        WRITE(IWR,1015)
 1015   FORMAT(/3X,'Energy',6X,'Rayleigh',6X,'Compton',4X,'Photoelect.',
     1    5X,'Pair',9X,'Total',/4X,'(eV)',8X,'(1/mtu)',6X,'(1/mtu)',6X,
     2    '(1/mtu)',6X,'(1/mtu)',6X,'(1/mtu)',/1X,77('-'))
      ENDIF
      DO I=1,NPHD
        IF(ER(I).GE.EL.AND.ER(I).LE.EU) THEN
          XEL=LOG(ER(I))
          XE=1.0D0+(XEL-DLEMP1)*DLFC
          KE=XE
          XEK=XE-KE
          PRA=EXP(SGRA(M,KE)+(SGRA(M,KE+1)-SGRA(M,KE))*XEK)/RHO(M)
          PCO=EXP(SGCO(M,KE)+(SGCO(M,KE+1)-SGCO(M,KE))*XEK)/RHO(M)
          IF(ER(I).LT.1.022D6) THEN
            PPP=1.0D-35
          ELSE
            PPP=EXP(SGPP(M,KE)+(SGPP(M,KE+1)-SGPP(M,KE))*XEK)/RHO(M)
          ENDIF
          PPH=XSR(I)*VMOL(M)/RHO(M)
          PT=PRA+PCO+PPP+PPH
          IF(INFO.GE.3) THEN
            WRITE(IWR,'(1P,E12.5,5E13.5)') ER(I),PRA,PCO,PPH,PPP,PT
          ENDIF
        ENDIF
      ENDDO
C
C  ****  Pair production. Initialization routine.
C
      CALL GPP0(M)
C
      DO IE=1,NEGP
        SGAUX(M,IE)=LOG(1.0D-35)
      ENDDO
C
      LNAME=' PENELOPE (v. 2005)  End of material data file ........'
      READ(IRD,'(A55)') NAME
      IF(NAME.NE.LNAME) THEN
        WRITE(IWR,'(/1X,''I/O error. Corrupt material data file.'')')
        WRITE(IWR,'(''      The last line is: '',A55)') NAME
        WRITE(IWR,'(''     ... and should be: '',A55)') LNAME
        STOP 'PEMATR. Corrupt material data file.'
      ENDIF
      WRITE(IWR,'(A55)') NAME
C
      RETURN
      END
C  *********************************************************************
C                       SUBROUTINE PEMATW
C  *********************************************************************
      SUBROUTINE PEMATW
C
C  This subroutine generates the material definition file, a part of the
C  input data file of PENELOPE.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N)
      CHARACTER*2 LASYMB
      CHARACTER*5 CH5
      CHARACTER*80 PFILE
      CHARACTER*60 NAME
      PARAMETER (A0B=5.291772083D-9)  ! Bohr radius (cm)
      PARAMETER (HREV=27.2113834D0)  ! Hartree energy (eV)
      PARAMETER (AVOG=6.02214199D23)  ! Avogadro's number
      PARAMETER (REV=5.10998902D5)  ! Electron rest energy (eV)
      PARAMETER (SL=137.03599976D0)  ! Speed of light (1/alpha)
      PARAMETER (TREV=2.0D0*REV)
      PARAMETER (PI=3.1415926535897932D0, FOURPI=4.0D0*PI)
C  ****  Composition data.
      PARAMETER (MAXMAT=10)
      COMMON/COMPOS/STF(MAXMAT,30),ZT(MAXMAT),AT(MAXMAT),RHO(MAXMAT),
     1  VMOL(MAXMAT),IZ(MAXMAT,30),NELEM(MAXMAT)
      DIMENSION FBW(30)
C  ****  Element data.
      COMMON/CADATA/ATW(99),EPX(99),RA1(99),RA2(99),RA3(99),RA4(99),
     1  RA5(99),RSCR(99),ETA(99),EB(99,30),IFI(99,30),IKS(99,30),
     2  NSHT(99),LASYMB(99)
      DIMENSION CP(99,30)
C  ****  Simulation parameters.
      COMMON/CSIMPA/EABS(3,MAXMAT),C1(MAXMAT),C2(MAXMAT),WCC(MAXMAT),
     1  WCR(MAXMAT)
C  ****  Energy grid and interpolation constants for the current energy.
      PARAMETER (NEGP=200)
      COMMON/CEGRID/EL,EU,ET(NEGP),DLEMP(NEGP),DLEMP1,DLFC,
     1  XEL,XE,XEK,KE
C  ****  E/P inelastic collisions.
      PARAMETER (NO=64)
      COMMON/CEIN/EXPOT(MAXMAT),OP2(MAXMAT),F(MAXMAT,NO),UI(MAXMAT,NO),
     1  WRI(MAXMAT,NO),KZ(MAXMAT,NO),KS(MAXMAT,NO),NOSC(MAXMAT)
C  ****  Compton scattering.
      PARAMETER (NOCO=64)
      COMMON/CGCO/FCO(MAXMAT,NOCO),UICO(MAXMAT,NOCO),FJ0(MAXMAT,NOCO),
     2  KZCO(MAXMAT,NOCO),KSCO(MAXMAT,NOCO),NOSCCO(MAXMAT)
      PARAMETER (NOM=400)
      DIMENSION FF(NOM),UUI(NOM),FFJ0(NOM),WWRI(NOM),KKZ(NOM),KKS(NOM)
      DIMENSION FFT(NOM),UIT(NOM),WRIT(NOM),KZT(NOM),KST(NOM)
      DIMENSION FC(NOM),UIC(NOM),FJ0C(NOM),KZC(NOM),KSC(NOM)
C  ****  Bremsstrahlung emission.
      PARAMETER (NBE=57, NBW=32)
      COMMON/CEBR/WB(NBW),PBCUT(MAXMAT,NEGP),WBCUT(MAXMAT,NEGP),
     1  PDFB(MAXMAT,NEGP,NBW),PACB(MAXMAT,NEGP,NBW),ZBR2(MAXMAT)
C
      DIMENSION EGRT(16)
      DATA EGRT/1.0D0,1.25D0,1.50D0,1.75D0,2.00D0,2.50D0,3.00D0,
     1 3.50D0,4.00D0,4.50D0,5.00D0,6.00D0,7.00D0,8.00D0,9.00D0,
     2 1.00D1/
C  ****  'Standard' energy grid.
      DIMENSION ES(NEGP),ESS(NEGP)
      DIMENSION EIT(NEGP),XG0(NEGP),PPE(NEGP),PPT(NEGP)
C
      M=1
C
      WRITE(6,*) '      '
      WRITE(6,*) ' Select one option (1 or 2):'
      WRITE(6,*) '   1: Enter composition data from the keyboard'
      WRITE(6,*) '   2: Read them from the file pdcompos.p05'
      READ(5,*) IREAD
      IF(IREAD.EQ.2) GO TO 901
C
C  ************  Data entered from keyboard.
C
      WRITE(6,*) ' Enter material name, for your information ',
     1  '(no more than 60 characters) ...'
      READ(5,'(A60)') NAME
      WRITE(6,'(1X,'' Material: '',A60)') NAME
C
      ZT(M)=0.0D0
      AT(M)=0.0D0
      EXPOT(M)=0.0D0
C  ****  Chemical formula or fractions by weight.
  800 CONTINUE
      WRITE(6,*) '      '
      WRITE(6,*) ' Chemical formula:'
      WRITE(6,*) ' Number of elements in the molecule...'
      READ(5,*) NELEM(M)
      IF(NELEM(M).LT.1.OR.NELEM(M).GT.30) THEN
        WRITE(6,*) ' NELEM must be positive and less than 31.'
        WRITE(6,*) ' Please, enter an allowed value.'
        GO TO 800
      ENDIF
      WRITE(6,'(1X,'' Number of elements = '',I2)') NELEM(M)
      WRITE(6,*) '      '
      WRITE(6,*) ' Select one option (1 or 2):'
      WRITE(6,*) '   1: Enter chemical (stoichiometric) formula'
      WRITE(6,*) '   2: Enter fraction by weight of each element'
      READ(5,*) IREAD2
C
      IF(IREAD2.EQ.2) THEN
        WRITE(6,*) '      '
        WRITE(6,*) ' Weight fractions...'
        DO I=1,NELEM(M)
 801      CONTINUE
          IF(I.EQ.1) THEN
            WRITE(6,231)
 231        FORMAT(/1X,' Enter atomic number and fraction by weight',
     1        ' of the first element ...')
          ELSE IF(I.EQ.2) THEN
            WRITE(6,232)
 232        FORMAT(/1X,' Enter atomic number and fraction by weight',
     1        ' of the second element ...')
          ELSE IF(I.EQ.3) THEN
            WRITE(6,233)
 233        FORMAT(/1X,' Enter atomic number and fraction by weight',
     1        ' of the third element ...')
          ELSE
            WRITE(6,234) I
 234        FORMAT(/1X,' Enter atomic number and fraction by weight',
     1        ' of the',I3,'-th element ...')
          ENDIF
          READ(5,*) IZZ,FBW(I)
          IF(IZZ.LT.1.OR.IZZ.GT.99) THEN
            WRITE(6,*) ' The atomic number must be in the range 1-99.'
            WRITE(6,*) ' Please, enter an allowed value.'
            GO TO 801
          ENDIF
          IF(FBW(I).LE.0.0D0) THEN
            WRITE(6,*) ' The fraction by weight must be positive.'
            WRITE(6,*) ' Please, enter a positive value.'
            GO TO 801
          ENDIF
          WRITE(6,*) '      '
          IZ(M,I)=IZZ
          IF(I.GT.1) THEN
            DO K=1,I-1
              IF(IZZ.EQ.IZ(M,K)) THEN
              WRITE(6,*) ' STOP. This element has been declared twice.'
                WRITE(6,*) ' Type ''enter'' to continue...'
                READ(*,*)
                STOP
              ENDIF
            ENDDO
          ENDIF
          WRITE(6,'(1X,'' Element: '',A2,'' (Z='',I2,''), fraction '',
     1    ''by weight ='',1P,E12.5)') LASYMB(IZZ),IZZ,FBW(I)
          STF(M,I)=FBW(I)/ATW(IZZ)
        ENDDO
C
        STFM=0.0D0
        DO I=1,NELEM(M)
          IF(STF(M,I).GT.STFM) THEN
            STFM=STF(M,I)
          ENDIF
        ENDDO
        IF(STFM.LT.1.0D-16) THEN
          WRITE(6,*) ' STOP. Fractions by weight are too small.'
          WRITE(6,*) ' Type ''enter'' to continue...'
          READ(*,*)
          STOP
        ENDIF
        DO I=1,NELEM(M)
          STF(M,I)=STF(M,I)/STFM
        ENDDO
C
        WRITE(6,*) '      '
        DO I=1,NELEM(M)
          IZZ=IZ(M,I)
          WRITE(6,'(1X,'' Element: '',A2,'' (Z='',I2,''), atoms/'',
     1      ''molecule ='',1P,E12.5)') LASYMB(IZZ),IZZ,STF(M,I)
          ZT(M)=ZT(M)+IZZ*STF(M,I)
          AT(M)=AT(M)+ATW(IZZ)*STF(M,I)
          EXPOT(M)=EXPOT(M)+IZZ*LOG(EPX(IZZ))*STF(M,I)
        ENDDO
      ELSE
C
        WRITE(6,*) '      '
        WRITE(6,*) ' Stoichiometric indices...'
        DO I=1,NELEM(M)
 802      CONTINUE
          IF(I.EQ.1) THEN
            WRITE(6,2001)
 2001     FORMAT(/1X,' Enter atomic number and number of atoms/molec',
     1      'ule of the first element ...')
          ELSE IF(I.EQ.2) THEN
            WRITE(6,2002)
 2002     FORMAT(/1X,' Enter atomic number and number of atoms/molec',
     1      'ule of the second element ...')
          ELSE IF(I.EQ.3) THEN
            WRITE(6,2003)
 2003     FORMAT(/1X,' Enter atomic number and number of atoms/molec',
     1      'ule of the third element ...')
          ELSE
            WRITE(6,2004) I
 2004     FORMAT(/1X,' Enter atomic number and number of atoms/molec',
     1      'ule of the',I3,'-th element ...')
          ENDIF
          READ(5,*) IZZ,STF(M,I)
          IF(IZZ.LT.1.OR.IZZ.GT.99) THEN
            WRITE(6,*) ' The atomic number must be in the range 1-99.'
            WRITE(6,*) ' Please, enter an allowed value.'
            GO TO 802
          ENDIF
          IF(STF(M,I).LE.0.0D0) THEN
            WRITE(6,*) ' The stoichiometric fraction must be positive.'
            WRITE(6,*) ' Please, enter a positive value.'
            GO TO 802
          ENDIF
          IZ(M,I)=IZZ
          WRITE(6,'(1X,'' Element: '',A2,'' (Z='',I2,''), atoms/'',
     1      ''molecule ='',1P,E12.5)') LASYMB(IZZ),IZZ,STF(M,I)
          IF(I.GT.1) THEN
            DO K=1,I-1
              IF(IZZ.EQ.IZ(M,K)) THEN
              WRITE(6,*) ' STOP. This element has been declared twice.'
                WRITE(6,*) ' Type ''enter'' to continue...'
                READ(*,*)
                STOP
              ENDIF
            ENDDO
          ENDIF
          ZT(M)=ZT(M)+IZZ*STF(M,I)
          AT(M)=AT(M)+ATW(IZZ)*STF(M,I)
          EXPOT(M)=EXPOT(M)+IZZ*LOG(EPX(IZZ))*STF(M,I)
        ENDDO
      ENDIF
C
      EXPOT(M)=EXP(EXPOT(M)/ZT(M))
      WRITE(6,*) '      '
      WRITE(6,'(1X,'' The calculated mean excitation energy I'',
     1  '' is '',1P,E12.5,'' eV'')') EXPOT(M)
      WRITE(6,*) ' Do you want to change it?   (1=yes,2=no)'
      READ(5,*) IYESNO
      IF(IYESNO.EQ.1) THEN
 803    CONTINUE
        WRITE(6,*) '      '
        WRITE(6,*) ' Enter mean excitation energy (eV) ...'
        READ(5,*) EXPOT(M)
        WRITE(6,'(1X,'' Mean excitation energy ='',1P,E12.5,
     1    '' eV'')') EXPOT(M)
        IF(EXPOT(M).LT.1.0D0) THEN
          WRITE(6,*) ' The mean exc. energy must be larger than 1 eV.'
          WRITE(6,*) ' Please, enter a valid value.'
          GO TO 803
        ENDIF
      ENDIF
C
 804  CONTINUE
      WRITE(6,*) '      '
      WRITE(6,*) ' Enter mass density (g/cm**3) ...'
      READ(5,*) RHO(M)
      WRITE(6,'(1X,'' Density = '',1P,E12.5,'' g/cm**3'')') RHO(M)
      IF(RHO(M).LE.0.0D0) THEN
        WRITE(6,*) ' The mass density must be positive.'
        WRITE(6,*) ' Please, enter a positive value.'
        GO TO 804
      ENDIF
      VMOL(M)=AVOG*RHO(M)/AT(M)
      GO TO 904
C
C  ************  Material data read from file compdata.p05.
C
  901 CONTINUE
      WRITE(6,*) ' Enter material identification number ...'
      READ(5,*) IDNUM
      IF(IDNUM.LT.1.OR.IDNUM.GT.280) THEN
        WRITE(6,*) ' STOP. The allowed material ID numbers are 1-280.'
        WRITE(6,*) ' Type ''enter'' to continue...'
        READ(*,*)
        STOP
      ENDIF
C
      OPEN(3,FILE='pdcompos.p05')
      DO I=1,15
        READ(3,'(A60)') NAME
      ENDDO
      DO K1=1,300
        READ(3,'(I3,1X,A60)',END=902) IORD,NAME
        READ(3,*) NELEM(M),ZOA,EXPOT(M),RHO(M)
        IF(NELEM(M).GT.30) THEN
          WRITE(6,*) ' STOP. NELEM cannot be larger than 30.'
          WRITE(6,*) ' Type ''enter'' to continue...'
          READ(*,*)
          STOP
        ENDIF
        IF(NELEM(M).LT.1) THEN
          WRITE(6,*) ' STOP. NELEM must be positive.'
          WRITE(6,*) ' Type ''enter'' to continue...'
          READ(*,*)
          STOP
        ENDIF
        DO I=1,NELEM(M)
          READ(3,*) IZ(M,I),XNOTH,STF(M,I)
        ENDDO
        IF(IORD.EQ.IDNUM) GO TO 903
      ENDDO
  902 CONTINUE
      WRITE(6,*) ' STOP. Abnormal termination of file ''pdcompos.p05''.'
      WRITE(6,*) ' Type ''enter'' to continue...'
      READ(*,*)
      STOP
  903 CONTINUE
      CLOSE(UNIT=3)
      WRITE(6,'(/2X,I3,1X,A60)') IORD,NAME
C
      ZT(M)=0.0D0
      AT(M)=0.0D0
      DO I=1,NELEM(M)
        IZZ=IZ(M,I)
        IF(IZZ.LT.1.OR.IZZ.GT.99) THEN
          WRITE(6,'(1X,'' Element:    (Z='',I2,''), atoms/'',
     1    ''molecule ='',1P,E12.5)') LASYMB(IZZ),IZZ,STF(M,I)
          WRITE(6,*) ' STOP. Wrong atomic number.'
          WRITE(6,*) ' Type ''enter'' to continue...'
          READ(*,*)
          STOP
        ENDIF
        WRITE(6,'(1X,'' Element: '',A2,'' (Z='',I2,''), atoms/'',
     1    ''molecule ='',1P,E12.5)') LASYMB(IZZ),IZZ,STF(M,I)
        IF(STF(M,I).LE.0.0D0) THEN
          WRITE(6,*) ' STOP. STF must be positive.'
          WRITE(6,*) ' Type ''enter'' to continue...'
          READ(*,*)
          STOP
        ENDIF
        IF(I.GT.1) THEN
          DO K=1,I-1
            IF(IZZ.EQ.IZ(M,K)) THEN
              WRITE(6,*) ' STOP. Element has been declared twice.'
              WRITE(6,*) ' Type ''enter'' to continue...'
              STOP
            ENDIF
          ENDDO
        ENDIF
        ZT(M)=ZT(M)+IZZ*STF(M,I)
        AT(M)=AT(M)+ATW(IZZ)*STF(M,I)
      ENDDO
      WRITE(6,'(/1X,'' Density = '',1P,E12.5,'' g/cm**3'')') RHO(M)
      WRITE(6,'(/1X,'' Number of electrons per molecule = '',
     1  1P,E12.5)') ZT(M)
      IF(RHO(M).LE.0.0D0) THEN
        WRITE(6,*) ' STOP. The density must be positive.'
        WRITE(6,*) ' Type ''enter'' to continue...'
        STOP
      ENDIF
      WRITE(6,'(1X,'' Mean excitation energy ='',1P,E12.5,
     1  '' eV'')') EXPOT(M)
      VMOL(M)=AVOG*RHO(M)/AT(M)
  904 CONTINUE
C
C  ************  Atomic configuration.
C
      DO I=1,99
        NSHT(I)=0
        DO J=1,30
          EB(I,J)=0.0D0
          CP(I,J)=0.0D0
          IFI(I,J)=0
          IKS(I,J)=0
        ENDDO
      ENDDO
C
      DO I=1,NELEM(M)
        IZZ=IZ(M,I)
C  ****  Loads element data only once. NSHT(IZZ) is used as a status
C        indicator.
        IF(NSHT(IZZ).EQ.0) THEN
          OPEN(3,FILE='pdatconf.p05')
          DO J=1,19
            READ(3,'(A5)') CH5
          ENDDO
          NS=0
          IZZT=0
          DO J=1,150000
            READ(3,2005,END=905) IIZ,IS,CH5,IIF,IE,CCP
 2005       FORMAT(I3,1X,I2,1X,A5,1X,I1,1X,I6,E9.2)
            IF(IIZ.EQ.IZZ) THEN
              NS=NS+1
              IF(NS.GT.30) THEN
                WRITE(6,'(/1X,''NS ='',I4)') NS
                WRITE(6,*) ' STOP. Too many shells.'
                WRITE(6,*) ' Type ''enter'' to continue...'
                READ(*,*)
                STOP
              ENDIF
              IF(IS.LT.1.OR.IS.GT.30) THEN
                WRITE(6,'(/1X,''IS ='',I4)') IS
                WRITE(6,*) ' STOP. Wrong shell number.'
                WRITE(6,*) ' Type ''enter'' to continue...'
                READ(*,*)
                STOP
              ENDIF
              IZZT=IZZT+IIF
              EB(IZZ,IS)=IE
              CP(IZZ,IS)=CCP
              IFI(IZZ,IS)=IIF
              IKS(IZZ,NS)=IS
            ENDIF
          ENDDO
  905     CONTINUE
          NSHT(IZZ)=NS
          IF(IZZ.NE.IZZT) THEN
            WRITE(6,*) ' STOP. Unbalanced charges (element).'
            WRITE(6,*) ' Type ''enter'' to continue...'
            READ(*,*)
            STOP
          ENDIF
          CLOSE(3)
        ENDIF
      ENDDO
C
C  ************  E/P inelastic scattering model parameters.
C
C  ****  Set the oscillator table (analogous to ICRU37).
C
      DO I=1,NO
        F(M,I)=0.0D0
        UI(M,I)=0.0D0
        WRI(M,I)=0.0D0
        KZ(M,I)=0
        KS(M,I)=0
      ENDDO
C
      DO I=1,NOM
         FF(I)=0.0D0
        UUI(I)=0.0D0
        WWRI(I)=0.0D0
        FFJ0(I)=0.0D0
        KKZ(I)=0
        KKS(I)=0
      ENDDO
      FT=0.0D0
C  ****  The 1st oscillator corresponds to the conduction band, which is
C  tentatively assumed to consist of outer shells with binding energies
C  less than 13 eV.
      NOS=1
      FF(1)=0.0D0
      UUI(1)=0.0D0
      FFJ0(1)=0.0D0
      KKZ(1)=0
      KKS(1)=30
      DO I=1,NELEM(M)
        IZZ=IZ(M,I)
        DO K=1,30
          JS=IKS(IZZ,K)
          IF(JS.GT.0) THEN
            IF(IFI(IZZ,JS).GT.0) THEN
              NOS=NOS+1
              IF(NOS.GT.NOM) THEN
                WRITE(6,*) ' STOP. Too many oscillators. The parameter',
     1            ' NOM should be enlarged.'
                WRITE(6,*) ' Type ''enter'' to continue...'
                READ(*,*)
                STOP
              ENDIF
              FF(NOS)=IFI(IZZ,JS)*STF(M,I)
              UUI(NOS)=EB(IZZ,JS)
              FFJ0(NOS)=CP(IZZ,JS)
              KKZ(NOS)=IZZ
              IF(IZZ.GT.6.AND.JS.LT.10) THEN
                KKS(NOS)=JS
              ELSE
                KKS(NOS)=30
              ENDIF
              FT=FT+FF(NOS)
              IF(EB(IZZ,K).LT.13.0D0) FF(1)=FF(1)+IFI(IZZ,JS)*STF(M,I)
            ENDIF
          ENDIF
        ENDDO
      ENDDO
C
      IF(ABS(FT-ZT(M)).GT.1.0D-10*ZT(M)) THEN
        WRITE(6,*) ' STOP. Unbalanced charges (compound).'
        WRITE(6,*) ' Type ''enter'' to continue...'
        READ(*,*)
        STOP
      ENDIF
C  ****  Oscillators are sorted by increasing energies.
      DO I=1,NOS-1
        DO J=I+1,NOS
          IF(UUI(I).GE.UUI(J)) THEN
            SAVE=UUI(I)
            UUI(I)=UUI(J)
            UUI(J)=SAVE
            SAVE=FF(I)
            FF(I)=FF(J)
            FF(J)=SAVE
            SAVE=FFJ0(I)
            FFJ0(I)=FFJ0(J)
            FFJ0(J)=SAVE
            ISAVE=KKZ(I)
            KKZ(I)=KKZ(J)
            KKZ(J)=ISAVE
            ISAVE=KKS(I)
            KKS(I)=KKS(J)
            KKS(J)=ISAVE
          ENDIF
        ENDDO
      ENDDO
C
C  ************  Plasma energy and conduction band excitations.
C
      OP2(M)=FOURPI*ZT(M)*VMOL(M)*A0B**3*HREV**2
      OMEGA=SQRT(OP2(M))
      EPP=OMEGA*SQRT(FF(1)/ZT(M))
      WRITE(6,'(/1X,'' Estimated oscillator strength and '',
     1  ''energy of the plasmon:''/2X,''Fcb ='',1P,E12.5,'',  Wcb ='',
     2  E12.5,'' eV'',/2X,''(for insulators, these quantiti'',
     1  ''es should be set equal to zero)'')') FF(1),EPP
C
      WRITE(6,'(/1X,'' Do you wish to change the Fcb and Wcb values? '',
     1  ''  (1=yes,2=no)''/2X,''(type 2 if you are not sure...)'')')
      READ(5,*) IPLOSP
C
      IF(IPLOSP.EQ.1) THEN
        WRITE(6,*) '      '
        WRITE(6,2007)
 2007   FORMAT(1X,' Enter the oscillator strength Fcb and ',
     1    'energy Wcb (in eV) of the plasmon ...')
        READ(5,*) FP,EP
        IF(FP.LT.0.5D0) THEN
          FP=0.0D0
          EP=0.0D0
          GO TO 975
        ENDIF
        IF(EP.LT.0.1D0) THEN
          EP=OMEGA*SQRT(FP/ZT(M))
        ENDIF
      ELSE
        EP=EPP
        FP=FF(1)
      ENDIF
  975 CONTINUE
      WRITE(6,'(/1X,'' Fcb ='',1P,E12.5,'',  Wcb ='',E12.5,'' eV'')')
     1   FP,EP
      IF(FP.GT.ZT(M)+1.0D-13) THEN
        WRITE(6,*) ' STOP. FP is too large.'
        WRITE(6,*) ' Type ''enter'' to continue...'
        READ(*,*)
        STOP
      ENDIF
      IF(EP.LT.1.0D0.OR.FP.LT.0.5D0) THEN
C  ****  Insulator. There is no conduction band.
        IFCB=0
        DO J=1,NOS-1
          FF(J)=FF(J+1)
          UUI(J)=UUI(J+1)
          FFJ0(J)=FFJ0(J+1)
          KKZ(J)=KKZ(J+1)
          KKS(J)=KKS(J+1)
        ENDDO
        NOS=NOS-1
      ELSE
C  ****  Conductor. Outer shells are 'moved' to the c.b.
        IFCB=1
        IDEAD=0
        FPP=FP
        I=1
  907   I=I+1
        IF(FF(I).LT.FPP) THEN
          FPP=FPP-FF(I)
          FF(I)=0.0D0
          IDEAD=IDEAD+1
          GO TO 907
        ELSE
          FF(I)=FF(I)-FPP
          IF(ABS(FF(I)).LT.1.0D-12) THEN
            FP=FP+FF(I)
            FF(I)=0.0D0
            IDEAD=IDEAD+1
          ENDIF
        ENDIF
        FF(1)=FP
        UUI(1)=0.0D0
        WWRI(1)=EP
        FFJ0(1)=0.75D0/SQRT(3.0D0*PI*PI*VMOL(M)*A0B**3*FP)
        KKZ(1)=0
        KKS(1)=30
        IF(IDEAD.GT.0) THEN
          DO J=2,NOS-IDEAD
            FF(J)=FF(J+IDEAD)
            UUI(J)=UUI(J+IDEAD)
            FFJ0(J)=FFJ0(J+IDEAD)
            KKZ(J)=KKZ(J+IDEAD)
            KKS(J)=KKS(J+IDEAD)
          ENDDO
          NOS=NOS-IDEAD
        ENDIF
      ENDIF
C  ****  Check f-sum rule.
      SUM=0.0D0
      DO J=1,NOS
        SUM=SUM+FF(J)
      ENDDO
      IF(ABS(SUM-ZT(M)).GT.1.0D-6*ZT(M)) THEN
        WRITE(6,*) ' STOP. Inconsistent oscillator strength data.'
        WRITE(6,*) ' Type ''enter'' to continue...'
        READ(*,*)
        STOP
      ENDIF
      IF(ABS(SUM-ZT(M)).GT.1.0D-12*ZT(M)) THEN
        FACT=ZT(M)/SUM
        DO J=1,NOS
          FF(J)=FACT*FF(J)
        ENDDO
      ENDIF
C
C  ****  Initial parameters for Compton scattering (before grouping).
C
      NOSTC=NOS
      CSUMT=0.0D0
      DO I=1,NOSTC
        FC(I)=FF(I)
        UIC(I)=UUI(I)
        FJ0C(I)=FFJ0(I)
        KZC(I)=KKZ(I)
        KSC(I)=KKS(I)
        CSUMT=CSUMT+FC(I)*FJ0C(I)
      ENDDO
C
C  ************  Sternheimer's adjustment factor.
C
      IF(NOS.GT.1) THEN
        TST=ZT(M)*DLOG(EXPOT(M))
        AAL=0.5D0
        AAU=10.0D0
  908   AA=0.5D0*(AAL+AAU)
        SUM=0.0D0
        DO I=1,NOS
          IF(I.EQ.1.AND.IFCB.EQ.1) THEN
            SUM=SUM+FF(1)*DLOG(WWRI(1))
          ELSE
            WI2=(AA*UUI(I))**2
     1         +0.666666666666666D0*(FF(I)/ZT(M))*OMEGA**2
            WWRI(I)=DSQRT(WI2)
            SUM=SUM+FF(I)*DLOG(WWRI(I))
          ENDIF
        ENDDO
        IF(SUM.LT.TST) THEN
          AAL=AA
        ELSE
          AAU=AA
        ENDIF
        IF(AAU-AAL.GT.1.0D-14*AA) GO TO 908
      ELSE
        UUI(1)=DABS(UUI(1))
        WWRI(1)=EXPOT(M)
      ENDIF
      WRITE(6,'(1X,'' Sternheimer adjustment factor = '',1P,E12.5)') AA
C  ****  Verification.
      EXPT=FF(1)*LOG(WWRI(1))
      TST=FF(1)
      IF(NOS.GT.1) THEN
        DO I=2,NOS
          EXPT=EXPT+FF(I)*LOG(WWRI(I))
          TST=TST+FF(I)
        ENDDO
      ENDIF
C
      IF(ABS(TST-ZT(M)).GT.1.0D-8*ZT(M)) THEN
        WRITE(6,*) ' STOP. Inconsistent oscillator strength data.'
        WRITE(6,*) ' TST-ZT(M) =',TST-ZT(M)
        WRITE(6,*) ' Type ''enter'' to continue...'
        READ(*,*)
        STOP
      ENDIF
      EXPT=EXP(EXPT/ZT(M))
      IF(ABS(EXPT-EXPOT(M)).GT.1.0D-8*EXPOT(M)) THEN
        WRITE(6,*) 'EXPT-EXPOT(M) =',EXPT-EXPOT(M)
        WRITE(6,*) 'Error in the calculation of the Sternheimer factor.'
        WRITE(6,*) '      '
        WRITE(6,'(1X,'' Number of oscillators  = '',I3)') NOS
        DO I=1,NOS
          WRITE(6,'(I4,1P,4E13.5,2I4)')
     1      I,FF(I),UUI(I),WWRI(I),FFJ0(I),KKZ(I),KKS(I)
        ENDDO
        WRITE(6,*) ' STOP. Inconsistent oscillator strength data (2).'
        WRITE(6,*) ' Type ''enter'' to continue...'
        READ(*,*)
        STOP
      ENDIF
C
C  ****  Oscillators with resonance energies differing by a factor less
C        than RGROUP are grouped as a single oscillator.
C
      RGROUP=1.05D0
      NOST=NOS
      DO I=1,NOST
        FFT(I)=FF(I)
        UIT(I)=UUI(I)
        WRIT(I)=WWRI(I)
        KZT(I)=KKZ(I)
        KST(I)=KKS(I)
      ENDDO
C  ****  Oscillators are sorted by increasing energies.
      IF(NOST.GT.1) THEN
        DO I=1,NOST-1
          DO J=I+1,NOST
            IF(WRIT(I).GT.WRIT(J)) THEN
              SAVE=FFT(I)
              FFT(I)=FFT(J)
              FFT(J)=SAVE
              SAVE=UIT(I)
              UIT(I)=UIT(J)
              UIT(J)=SAVE
              SAVE=WRIT(I)
              WRIT(I)=WRIT(J)
              WRIT(J)=SAVE
              ISAVE=KZT(I)
              KZT(I)=KZT(J)
              KZT(J)=ISAVE
              ISAVE=KST(I)
              KST(I)=KST(J)
              KST(J)=ISAVE
            ENDIF
          ENDDO
        ENDDO
      ENDIF
C
  910 CONTINUE
      IELIM=0
      IF(NOST.GT.2) THEN
        DO 911 I=1,NOST-1
        IF(WRIT(I).LT.1.0D0.OR.WRIT(I+1).LT.1.0D0) GO TO 911
        IF(WRIT(I+1).GT.RGROUP*WRIT(I)) GO TO 911
        WRIT(I)=EXP((FFT(I)*LOG(WRIT(I))+FFT(I+1)*LOG(WRIT(I+1)))
     1         /(FFT(I)+FFT(I+1)))
        UIT(I)=(FFT(I)*UIT(I)+FFT(I+1)*UIT(I+1))/(FFT(I)+FFT(I+1))
        FFT(I)=FFT(I)+FFT(I+1)
        IF(KZT(I).NE.KZT(I+1)) KZT(I)=0
        IF(KST(I).LT.10.OR.KST(I+1).LT.10) THEN
          KST(I)=MIN(KST(I),KST(I+1))
        ELSE
          KST(I)=30
        ENDIF
        IF(I.LT.NOST-1) THEN
          DO J=I+1,NOST-1
            FFT(J)=FFT(J+1)
            UIT(J)=UIT(J+1)
            WRIT(J)=WRIT(J+1)
            KZT(J)=KZT(J+1)
            KST(J)=KST(J+1)
          ENDDO
        ENDIF
        IELIM=IELIM+1
        FFT(NOST)=0.0D0
        UIT(NOST)=0.0D0
        WRIT(NOST)=0.0D0
        KZT(NOST)=0
        KST(NOST)=0
  911   CONTINUE
      ENDIF
      IF(IELIM.GT.0) THEN
        NOST=NOST-IELIM
        GO TO 910
      ENDIF
C  ****  E/P inelastic model parameters transferred to the final
C        arrays.
      IF(NOST.LT.NO) THEN
        NOSC(M)=NOST
        DO I=1,NOSC(M)
          F(M,I)=FFT(I)
          UI(M,I)=UIT(I)
          WRI(M,I)=WRIT(I)
          KZ(M,I)=KZT(I)
          KS(M,I)=KST(I)
        ENDDO
      ELSE
        RGROUP=RGROUP**2
        GO TO 910
      ENDIF
      WRITE(6,'(1X,'' E/P in. grouping factor = '',1P,E12.5)') RGROUP
C
C  ************  Compton (impulse approximation) parameters.
C
C  ****  Selection of the lowest ionization energy for inner shells.
C  Only the K, L and M shells with ionization energy less than that of
C  the N1 shell of the heaviest element in the material are considered
C  as inner shells. As a result, the inner/outer character of an atomic
C  shell depends on the composition of the material.
C
      IZMAX=0
      DO I=1,NELEM(M)
        IZMAX=MAX(IZMAX,IZ(M,I))
      ENDDO
      WISCUT=MAX(200.0D0,EB(IZMAX,10))
C
C  ****  Outer shells with ionization energies differing by a factor
C  less than RGROUP are grouped as a single shell.
C
      RGROUP=1.50D0
C  ****  Shells are sorted by increasing ionization energies.
      IF(NOSTC.GT.1) THEN
        DO I=1,NOSTC-1
          DO J=I+1,NOSTC
            IF(UIC(I).GT.UIC(J)) THEN
              SAVE=FC(I)
              FC(I)=FC(J)
              FC(J)=SAVE
              SAVE=UIC(I)
              UIC(I)=UIC(J)
              UIC(J)=SAVE
              SAVE=FJ0C(I)
              FJ0C(I)=FJ0C(J)
              FJ0C(J)=SAVE
              ISAVE=KZC(I)
              KZC(I)=KZC(J)
              KZC(J)=ISAVE
              ISAVE=KSC(I)
              KSC(I)=KSC(J)
              KSC(J)=ISAVE
            ENDIF
          ENDDO
        ENDDO
      ENDIF
C
  920 CONTINUE
      IELIM=0
      IF(NOSTC.GT.2) THEN
        DO 921 I=1,NOST-1
        IF(UIC(I).GT.WISCUT) GO TO 921
        IF(UIC(I).LT.1.0D0.OR.UIC(I+1).LT.1.0D0) GO TO 921
        IF(UIC(I+1).GT.RGROUP*UIC(I)) GO TO 921
        UIC(I)=(FC(I)*UIC(I)+FC(I+1)*UIC(I+1))/(FC(I)+FC(I+1))
        FJ0C(I)=(FC(I)*FJ0C(I)+FC(I+1)*FJ0C(I+1))/(FC(I)+FC(I+1))
        FC(I)=FC(I)+FC(I+1)
        IF(KZC(I).NE.KZC(I+1)) KZC(I)=0
        KSC(I)=30
        IF(I.LT.NOSTC-1) THEN
          DO J=I+1,NOSTC-1
            FC(J)=FC(J+1)
            UIC(J)=UIC(J+1)
            FJ0C(J)=FJ0C(J+1)
            KZC(J)=KZC(J+1)
            KSC(J)=KSC(J+1)
          ENDDO
        ENDIF
        IELIM=IELIM+1
        FC(NOSTC)=0.0D0
        UIC(NOSTC)=0.0D0
        FJ0C(NOSTC)=0.0D0
        KZC(NOSTC)=0
        KSC(NOSTC)=0
  921   CONTINUE
      ENDIF
      IF(IELIM.GT.0) THEN
        NOSTC=NOSTC-IELIM
        GO TO 920
      ENDIF
C  ****  Compton scattering model parameters transferred to the final
C        arrays.
      IF(NOSTC.LT.NOCO) THEN
        NOSCCO(M)=NOSTC
        DO I=1,NOSCCO(M)
          FCO(M,I)=FC(I)
          UICO(M,I)=UIC(I)
          FJ0(M,I)=FJ0C(I)
          KZCO(M,I)=KZC(I)
          KSCO(M,I)=KSC(I)
          CSUMT=CSUMT-FCO(M,I)*FJ0(M,I)
        ENDDO
        IF(ABS(CSUMT).GT.1.0D-9) THEN
          WRITE(6,*) ' STOP. Error in grouping the Compton profiles.'
          WRITE(6,'(''  Residual sum ='',1p,e12.5)') ABS(CSUMT)
          WRITE(6,*) ' Type ''enter'' to continue...'
          READ(*,*)
          STOP
        ENDIF
      ELSE
        RGROUP=RGROUP**2
        GO TO 920
      ENDIF
      WRITE(6,'(1X,'' Compton grouping factor = '',1P,E12.5)') RGROUP
C
C  ************  PENELOPE's input file.
C
      WRITE(6,*) '      '
      WRITE(6,*) ' PENELOPE''s material data file is being created.'
      WRITE(6,*) ' Enter path+name for this file (up to 80 characte',
     1   'rs) ...'
      READ(5,'(A80)') PFILE
      OPEN(7,FILE=PFILE)
C
      WRITE(7,1000)
 1000 FORMAT(1X,
     1  'PENELOPE (v. 2005)  Material data file ...............')
      WRITE(7,1001) NAME
 1001 FORMAT(' Material: ',A60)
      WRITE(7,1002) RHO(M)
 1002 FORMAT(' Mass density =',1P,E15.8,' g/cm**3')
      WRITE(7,1003) NELEM(M)
 1003 FORMAT(' Number of elements in the molecule = ',I2)
      DO I=1,NELEM(M)
        WRITE(7,1004) IZ(M,I),STF(M,I)
 1004   FORMAT('   atomic number =',I3,',  atoms/molecule =',1P,E15.8)
      ENDDO
      WRITE(7,1005) EXPOT(M)
 1005 FORMAT(' Mean excitation energy =',1P,E15.8,' eV')
C
      WRITE(7,1006) NOSC(M)
 1006 FORMAT(' Number of oscillators =',I3,' (E/P inelastic model)')
      DO I=1,NOSC(M)
        WRITE(7,1007) I,F(M,I),UI(M,I),WRI(M,I),KZ(M,I),KS(M,I)
 1007   FORMAT(I4,1P,3E16.8,2I4)
      ENDDO
C
      WRITE(7,1106) NOSCCO(M)
 1106 FORMAT(' Number of shells =',I3,' (Compton IA model)')
      DO I=1,NOSCCO(M)
        WRITE(7,1107) I,FCO(M,I),UICO(M,I),FJ0(M,I),KZCO(M,I),
     1    KSCO(M,I)
 1107   FORMAT(I4,1P,3E16.8,2I4)
        FJ0(M,I)=FJ0(M,I)*SL
      ENDDO
C
C  ****  Atomic relaxation data.
C
      DO I=1,NELEM(M)
        IZZ=IZ(M,I)
        CALL RELAXW(IZZ,7)
      ENDDO
C
C  ****  Energy grid (standard).
C
      NES=0
      IGRID=0
      FGRID=1.0D0
  101 IGRID=IGRID+1
      EV=EGRT(IGRID)*FGRID
      IF(IGRID.EQ.16) THEN
        IGRID=1
        FGRID=10.0D0*FGRID
      ENDIF
      IF(EV.LT.49.0D0) GO TO 101
      NES=NES+1
      ES(NES)=EV
      IF(EV.LT.1.0D9) GO TO 101
C
      EMIN=50.0D0
      EMAX=1.0D9
      CALL EGRID(EMIN,EMAX)
C
      WCRM=10.0D0
      WCCM=0.0D0
C
C  ****  Bremsstrahlung emission,
C
      CALL EBRW(WCRM,M,7)
      ZEQ=SQRT(ZBR2(M))
      CALL BRAW(ZEQ,7)
C
      WRITE(7,1008) NES
 1008 FORMAT(1X,'*** Stopping powers for electrons and positrons',
     1  ',  NDATA =',I4)
      DO IE=1,NES
        EE=ES(IE)
        CALL EINTX(EE,WCCM,XH0,XH1,XH2,XS0,XS1,XS2,XT1,XT2,DELTA,M)
        ESTP=(XS1+XH1)*VMOL(M)*1.0D-6/RHO(M)
        CALL EBRTX(EE,WCRM,XH0,XH1,XH2,XS1,XS2)
        ERSTP=(XS1+XH1)*VMOL(M)*1.0D-6/RHO(M)
        CALL PINTX(EE,WCCM,XH0,XH1,XH2,XS0,XS1,XS2,XT1,XT2,DELTA,M)
        PSTP=(XS1+XH1)*VMOL(M)*1.0D-6/RHO(M)
        CALL PBRTX(EE,WCRM,XH0,XH1,XH2,XS1,XS2)
        PRSTP=(XS1+XH1)*VMOL(M)*1.0D-6/RHO(M)
        WRITE(7,'(1P,E10.3,5E12.5)') EE,ESTP,ERSTP,PSTP,PRSTP
      ENDDO
C
C  **** Electron and positron elastic x-sections.
C
      CALL EELW(M,7)
      CALL EELdW(M,7)  ! Uses the ELSEP database.
C
C  **** Electron and positron inner shell ionization x-sections.
C
      CALL ESIW(M,7)
      CALL PSIW(M,7)
C
C  ****  Photon x-sections.
C
      CALL GPPW(EIT,XG0,NPTAB,M)
      DO I=1,NES
        PPE(I)=0.0D0
      ENDDO
      CALL MERGE2(ES,PPE,EIT,XG0,ESS,PPT,NES,NPTAB,NESS)
C
      WRITE(7,1009) NESS
 1009 FORMAT(1X,'*** Rayleigh, Compton and pair-production cross',
     1  ' sections,  NDATA =',I4)
      DO IE=1,NESS
        EE=ESS(IE)
        CALL GCOTX(EE,CSC,M)
        IF(CSC.LT.1.0D-35) CSC=0.0D0
        CALL GRATX(EE,CSR,M)
        IF(EE.LT.TREV+3.0D0) PPT(IE)=0.0D0
        WRITE(7,'(1P,E10.3,5E12.5)') EE,CSR,CSC,PPT(IE)
      ENDDO
C
      CALL GPHW(M,7)
C
      WRITE(7,1099)
 1099 FORMAT(1X,
     1  'PENELOPE (v. 2005)  End of material data file ........')
      CLOSE(7)
      RETURN
      END
C  *********************************************************************
C                  FUNCTION PRANGE
C  *********************************************************************
      FUNCTION PRANGE(E,KPAR,M)
C
C  This function computes the range (in cm) of particles of type KPAR
C  and energy E in material M. For electrons and positrons, the output
C  is the CSDA range. For photons, the delivered value is the mean free
C  path (=inverse attenuation coefficient).
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N)
C  ****  Energy grid and interpolation constants for the current energy.
      PARAMETER (MAXMAT=10)
      PARAMETER (NEGP=200)
      COMMON/CEGRID/EL,EU,ET(NEGP),DLEMP(NEGP),DLEMP1,DLFC,
     1  XEL,XE,XEK,KE
      COMMON/CRANGE/RANGE(3,MAXMAT,NEGP),RANGEL(3,MAXMAT,NEGP)
C
      IF(E.LT.EL) THEN
        EE=EL
      ELSE IF (E.GT.EU) THEN
        EE=EU
      ELSE
        EE=E
      ENDIF
      XEL=LOG(EE)
      XE=1.0D0+(XEL-DLEMP1)*DLFC
      KE=XE
      IF(KE.LT.1) KE=1
      IF(KE.GE.NEGP) KE=NEGP-1
      XEK=XE-KE
C
      PRANGE=EXP(RANGEL(KPAR,M,KE)
     1      +(RANGEL(KPAR,M,KE+1)-RANGEL(KPAR,M,KE))*XEK)
      RETURN
      END
C  *********************************************************************
C                  FUNCTION PHMFP
C  *********************************************************************
      FUNCTION PHMFP(E,KPAR,M,ICOL)
C
C  This function computes the mean free path (in cm) of particles of
C  type KPAR and energy E between hard interactions of kind ICOL in
C  material M. If ICOL does not correspond to a hard interaction type,
C  the result is set equal to 1.0D16.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N)
C  ****  Energy grid and interpolation constants for the current energy.
      PARAMETER (NEGP=200)
      COMMON/CEGRID/EL,EU,ET(NEGP),DLEMP(NEGP),DLEMP1,DLFC,
     1  XEL,XE,XEK,KE
C  ****  Composition data.
      PARAMETER (MAXMAT=10)
      COMMON/COMPOS/STF(MAXMAT,30),ZT(MAXMAT),AT(MAXMAT),RHO(MAXMAT),
     1  VMOL(MAXMAT),IZ(MAXMAT,30),NELEM(MAXMAT)
C  ****  Electron simulation tables.
      COMMON/CEIMFP/SEHEL(MAXMAT,NEGP),SEHIN(MAXMAT,NEGP),
     1  SEISI(MAXMAT,NEGP),SEHBR(MAXMAT,NEGP),SEAUX(MAXMAT,NEGP),
     2  SETOT(MAXMAT,NEGP),CSTPE(MAXMAT,NEGP),RSTPE(MAXMAT,NEGP),
     3  DEL(MAXMAT,NEGP),W1E(MAXMAT,NEGP),W2E(MAXMAT,NEGP),
     4  RNDCE(MAXMAT,NEGP),AE(MAXMAT,NEGP),BE(MAXMAT,NEGP),
     5  T1E(MAXMAT,NEGP),T2E(MAXMAT,NEGP)
C  ****  Photon simulation tables.
      COMMON/CGIMFP/SGRA(MAXMAT,NEGP),SGCO(MAXMAT,NEGP),
     1  SGPH(MAXMAT,NEGP),SGPP(MAXMAT,NEGP),SGAUX(MAXMAT,NEGP)
C  ****  Photoelectric cross sections.
      PARAMETER (NTP=8000)
      COMMON/CGPH00/EPH(NTP),XPH(NTP,10),IPHF(99),IPHL(99),NPHS(99),NCUR
C  ****  Positron simulation tables.
      COMMON/CPIMFP/SPHEL(MAXMAT,NEGP),SPHIN(MAXMAT,NEGP),
     1  SPISI(MAXMAT,NEGP),SPHBR(MAXMAT,NEGP),SPAN(MAXMAT,NEGP),
     2  SPAUX(MAXMAT,NEGP),SPTOT(MAXMAT,NEGP),CSTPP(MAXMAT,NEGP),
     3  RSTPP(MAXMAT,NEGP),W1P(MAXMAT,NEGP),W2P(MAXMAT,NEGP),
     4  RNDCP(MAXMAT,NEGP),AP(MAXMAT,NEGP),BP(MAXMAT,NEGP),
     5  T1P(MAXMAT,NEGP),T2P(MAXMAT,NEGP)
C
      IF(E.LE.EL.OR.E.GE.EU) THEN
        PHMFP=1.0D50
        RETURN
      ENDIF
      XEL=LOG(E)
      XE=1.0D0+(XEL-DLEMP1)*DLFC
      KE=XE
      IF(KE.LT.1) KE=1
      IF(KE.GE.NEGP) KE=NEGP-1
      XEK=XE-KE
C
      HMFP=1.0D-35
      IF(KPAR.EQ.1) THEN
        IF(ICOL.EQ.2) THEN
          HMFP=EXP(SEHEL(M,KE)+(SEHEL(M,KE+1)-SEHEL(M,KE))*XEK)
        ELSE IF(ICOL.EQ.3) THEN
          HMFP=EXP(SEHIN(M,KE)+(SEHIN(M,KE+1)-SEHIN(M,KE))*XEK)
        ELSE IF(ICOL.EQ.4) THEN
          HMFP=EXP(SEHBR(M,KE)+(SEHBR(M,KE+1)-SEHBR(M,KE))*XEK)
        ELSE IF(ICOL.EQ.5) THEN
          HMFP=EXP(SEISI(M,KE)+(SEISI(M,KE+1)-SEISI(M,KE))*XEK)
        ELSE IF(ICOL.EQ.8) THEN
          HMFP=EXP(SEAUX(M,KE)+(SEAUX(M,KE+1)-SEAUX(M,KE))*XEK)
        ENDIF
      ELSE IF(KPAR.EQ.2) THEN
        IF(ICOL.EQ.1) THEN
          HMFP=EXP(SGRA(M,KE)+(SGRA(M,KE+1)-SGRA(M,KE))*XEK)
        ELSE IF(ICOL.EQ.2) THEN
          HMFP=EXP(SGCO(M,KE)+(SGCO(M,KE+1)-SGCO(M,KE))*XEK)
        ELSE IF(ICOL.EQ.3) THEN
          PTOT=0.0D0
          DO IEL=1,NELEM(M)
            IZZ=IZ(M,IEL)
C  ****  Binary search.
            I=IPHF(IZZ)
            IU=IPHL(IZZ)
    1       IT=(I+IU)/2
            IF(XEL.GT.EPH(IT)) THEN
              I=IT
            ELSE
              IU=IT
            ENDIF
            IF(IU-I.GT.1) GO TO 1
C
            DEE=EPH(I+1)-EPH(I)
            IF(DEE.GT.1.0D-15) THEN
              PCSL=XPH(I,1)+(XPH(I+1,1)-XPH(I,1))*(XEL-EPH(I))/DEE
            ELSE
              PCSL=XPH(I,1)
            ENDIF
            PTOT=PTOT+STF(M,IEL)*EXP(PCSL)
          ENDDO
          HMFP=PTOT*VMOL(M)
        ELSE IF(ICOL.EQ.4.AND.E.GT.1.022D6) THEN
          HMFP=EXP(SGPP(M,KE)+(SGPP(M,KE+1)-SGPP(M,KE))*XEK)
        ELSE IF(ICOL.EQ.6) THEN
          HMFP=EXP(SGAUX(M,KE)+(SGAUX(M,KE+1)-SGAUX(M,KE))*XEK)
        ENDIF
      ELSE IF(KPAR.EQ.3) THEN
        IF(ICOL.EQ.2) THEN
          HMFP=EXP(SPHEL(M,KE)+(SPHEL(M,KE+1)-SPHEL(M,KE))*XEK)
        ELSE IF(ICOL.EQ.3) THEN
          HMFP=EXP(SPHIN(M,KE)+(SPHIN(M,KE+1)-SPHIN(M,KE))*XEK)
        ELSE IF(ICOL.EQ.4) THEN
          HMFP=EXP(SPHBR(M,KE)+(SPHBR(M,KE+1)-SPHBR(M,KE))*XEK)
        ELSE IF(ICOL.EQ.5) THEN
          HMFP=EXP(SPISI(M,KE)+(SPISI(M,KE+1)-SPISI(M,KE))*XEK)
        ELSE IF(ICOL.EQ.6) THEN
          HMFP=EXP(SPAN(M,KE)+(SPAN(M,KE+1)-SPAN(M,KE))*XEK)
        ELSE IF(ICOL.EQ.8) THEN
          HMFP=EXP(SPAUX(M,KE)+(SPAUX(M,KE+1)-SPAUX(M,KE))*XEK)
        ENDIF
      ENDIF
C
      PHMFP=1.0D0/MAX(HMFP,1.0D-35)
      RETURN
      END
C  *********************************************************************
C                       SUBROUTINE EEL
C  *********************************************************************
      SUBROUTINE EEL(A,B,RNDC,RMU)
C
C  Simulation of hard elastic events. MW model.
C
C  Input arguments:
C    A, B ... angular distribution parameters.
C    RNDC ... cutoff probability.
C  Output values:
C    RMU .... angular deflection, =(1-CDT)/2.
C
      use mtpar
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N)
!     ! EXTERNALRAND
C
      A1=A+1.0D0
      IF(B.GE.0.0D0) THEN
C
C  ****  Case I.
C
        RMUAV=A*A1*LOG(A1/A)-A
        B1=1.0D0-B
        RND0=B1*A1*RMUAV/(A+RMUAV)
        status=rmt(r)
        RND=RNDC+r*(1.0D0-RNDC)
        IF(RND.LT.RND0) THEN
          RMU=RND*A/(B1*A1-RND)
        ELSE IF(RND.GT.RND0+B) THEN
          RNDMB=RND-B
          RMU=RNDMB*A/(B1*A1-RNDMB)
        ELSE
          RMU=RMUAV
        ENDIF
      ELSE
C
C  ****  Case II.
C
        BB=-B
        B1=1.0D0-BB
        RMUC=RNDC*A/(B1*A1-RNDC)
        PW=B1*A*(1.0D0-RMUC)/(A+RMUC)
        status=rmt(r)
        IF(r*(BB+PW).LT.BB) THEN
          status=rmt(r)
          RMU=0.5D0*(1.0D0+SQRT(r))
        ELSE
          status=rmt(r)
          RNDRC=r*(1.0D0-RMUC)
          RMU=(A*RNDRC+A1*RMUC)/(A1-RNDRC)
        ENDIF
      ENDIF
C
      RETURN
      END
C  *********************************************************************
C                       SUBROUTINE EELR
C  *********************************************************************
      SUBROUTINE EELR(M,IRD,IWR,INFO)
C
C  This subroutine reads elastic cross sections for electrons and posi-
C  trons in material M from the material data file. It also initializes
C  the algorithm for simulation of elastic scattering of electrons and
C  positrons.
C
      use mtpar
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N)
C  ****  Simulation parameters.
      PARAMETER (MAXMAT=10)
      COMMON/CSIMPA/EABS(3,MAXMAT),C1(MAXMAT),C2(MAXMAT),WCC(MAXMAT),
     1  WCR(MAXMAT)
C  ****  Composition data.
      COMMON/COMPOS/STF(MAXMAT,30),ZT(MAXMAT),AT(MAXMAT),RHO(MAXMAT),
     1  VMOL(MAXMAT),IZ(MAXMAT,30),NELEM(MAXMAT)
C  ****  Energy grid and interpolation constants for the current energy.
      PARAMETER (NEGP=200)
      COMMON/CEGRID/EL,EU,ET(NEGP),DLEMP(NEGP),DLEMP1,DLFC,
     1  XEL,XE,XEK,KE
C  ****  Electron simulation tables.
      COMMON/CEIMFP/SEHEL(MAXMAT,NEGP),SEHIN(MAXMAT,NEGP),
     1  SEISI(MAXMAT,NEGP),SEHBR(MAXMAT,NEGP),SEAUX(MAXMAT,NEGP),
     2  SETOT(MAXMAT,NEGP),CSTPE(MAXMAT,NEGP),RSTPE(MAXMAT,NEGP),
     3  DEL(MAXMAT,NEGP),W1E(MAXMAT,NEGP),W2E(MAXMAT,NEGP),
     4  RNDCE(MAXMAT,NEGP),AE(MAXMAT,NEGP),BE(MAXMAT,NEGP),
     5  T1E(MAXMAT,NEGP),T2E(MAXMAT,NEGP)
C  ****  Positron simulation tables.
      COMMON/CPIMFP/SPHEL(MAXMAT,NEGP),SPHIN(MAXMAT,NEGP),
     1  SPISI(MAXMAT,NEGP),SPHBR(MAXMAT,NEGP),SPAN(MAXMAT,NEGP),
     2  SPAUX(MAXMAT,NEGP),SPTOT(MAXMAT,NEGP),CSTPP(MAXMAT,NEGP),
     3  RSTPP(MAXMAT,NEGP),W1P(MAXMAT,NEGP),W2P(MAXMAT,NEGP),
     4  RNDCP(MAXMAT,NEGP),AP(MAXMAT,NEGP),BP(MAXMAT,NEGP),
     5  T1P(MAXMAT,NEGP),T2P(MAXMAT,NEGP)
C  ****  Total and transport cross sections.
      COMMON/CEEL00/EIT(NEGP),XE0(NEGP),XE1(NEGP),XE2(NEGP),XP0(NEGP),
     1  XP1(NEGP),XP2(NEGP),T1E0(NEGP),T2E0(NEGP),T1P0(NEGP),T2P0(NEGP),
     2  EITL(NEGP),FL(NEGP),A(NEGP),B(NEGP),C(NEGP),D(NEGP)
      COMMON/CEINTF/T1EI(NEGP),T2EI(NEGP),T1PI(NEGP),T2PI(NEGP)
C
C  ****  Reading the input cross section table.
C
      READ(IRD,2001) NDATA
 2001 FORMAT(60X,I4)
      IF(INFO.GE.2) WRITE(IWR,1001) NDATA
 1001 FORMAT(/1X,'***  Electron and positron elastic cross sections',
     1  ',  NDATA =',I4)
      IF(INFO.GE.2) WRITE(IWR,1101)
 1101 FORMAT(/2X,'Energy',7X,'CS0,e-',6X,'CS1,e-',6X,'CS2,e-',6X,
     1  'CS0,e+',6X,'CS1,e+',6X,'CS2,e+',/3X,'(eV)',8X,'(cm**2)',
     2  5X,'(cm**2)',5X,'(cm**2)',5X,'(cm**2)',5X,'(cm**2)',5X,
     3  '(cm**2)',/1X,84('-'))
      DO I=1,NDATA
        READ(IRD,*) EIT(I),XE0(I),XE1(I),XE2(I),XP0(I),XP1(I),XP2(I)
        IF(INFO.GE.2) WRITE(IWR,'(1P,7E12.5)')
     1    EIT(I),XE0(I),XE1(I),XE2(I),XP0(I),XP1(I),XP2(I)
        EITL(I)=LOG(EIT(I))
      ENDDO
C
C  ****  Elastic scattering of electrons.
C
      DO I=1,NDATA
        FL(I)=LOG(XE0(I))
      ENDDO
      CALL SPLINE(EITL,FL,A,B,C,D,0.0D0,0.0D0,NDATA)
      DO I=1,NEGP
        EC=DLEMP(I)
        CALL FINDI(EITL,EC,NDATA,J)
        XE0(I)=EXP(A(J)+EC*(B(J)+EC*(C(J)+EC*D(J))))
      ENDDO
C
      DO I=1,NDATA
        FL(I)=LOG(XE1(I))
      ENDDO
      CALL SPLINE(EITL,FL,A,B,C,D,0.0D0,0.0D0,NDATA)
      DO I=1,NEGP
        EC=DLEMP(I)
        CALL FINDI(EITL,EC,NDATA,J)
        XE1(I)=EXP(A(J)+EC*(B(J)+EC*(C(J)+EC*D(J))))
      ENDDO
C
      DO I=1,NDATA
        FL(I)=LOG(XE2(I))
      ENDDO
      CALL SPLINE(EITL,FL,A,B,C,D,0.0D0,0.0D0,NDATA)
      DO I=1,NEGP
        EC=DLEMP(I)
        CALL FINDI(EITL,EC,NDATA,J)
        XE2(I)=EXP(A(J)+EC*(B(J)+EC*(C(J)+EC*D(J))))
      ENDDO
C
      DO I=1,NEGP
        XS0=XE0(I)
        XS1=XE1(I)
        XS2=XE2(I)
        FPEL=1.0D0/(XS0*VMOL(M))
        FPT1=1.0D0/(XS1*VMOL(M))
        FPST=ET(I)/(CSTPE(M,I)+RSTPE(M,I))
        XS0H=1.0D0/(VMOL(M)*MAX(FPEL,MIN(C1(M)*FPT1,C2(M)*FPST)))
        CALL EEL0(XS0,XS1,XS2,XS0H,AAA,BBB,RNDC,XS1S,XS2S)
        SEHEL(M,I)=XS0H*VMOL(M)
        RNDCE(M,I)=RNDC
        AE(M,I)=AAA
        BE(M,I)=BBB
        T1E0(I)=XS1S
        T1E(M,I)=T1EI(I)+XS1S*VMOL(M)
        T2E0(I)=XS2S
        T2E(M,I)=T2EI(I)+XS2S*VMOL(M)
      ENDDO
C
C  ****  Print electron elastic scattering tables.
C
      IF(INFO.GE.3) WRITE(IWR,1002)
 1002 FORMAT(/1X,'PENELOPE >>>  Elastic scattering of electrons')
      IF(INFO.GE.3) WRITE(IWR,1003)
 1003 FORMAT(/3X,'E (eV)',6X,'MFP (mtu)',3X,'TMFP1 (mtu)',2X,
     1  'MFPh (mtu)',8X,'A',11X,'B',11X,'RNDC',/1X,90('-'))
      DO I=1,NEGP
        FP0=RHO(M)/(XE0(I)*VMOL(M))
        FP1=RHO(M)/(XE1(I)*VMOL(M))
        HMFP=RHO(M)/SEHEL(M,I)
        IF(INFO.GE.3) WRITE(IWR,'(1P,7(E12.5,1X))') ET(I),FP0,FP1,
     1    HMFP,AE(M,I),BE(M,I),RNDCE(M,I)
        SEHEL(M,I)=LOG(SEHEL(M,I))
        AE(M,I)=LOG(AE(M,I))
C  ****  Soft scattering events are switched off when T1E is too small.
        IF(T1E(M,I).GT.1.0D-6*XE1(I)*VMOL(M)) THEN
          T1E(M,I)=LOG(MAX(T1E(M,I),1.0D-35))
          T2E(M,I)=LOG(MAX(T2E(M,I),1.0D-35))
        ELSE
          T1E(M,I)=-100.0D0
          T2E(M,I)=-100.0D0
        ENDIF
      ENDDO
C
C  ****  Elastic scattering of positrons.
C
      DO I=1,NDATA
        FL(I)=LOG(XP0(I))
      ENDDO
      CALL SPLINE(EITL,FL,A,B,C,D,0.0D0,0.0D0,NDATA)
      DO I=1,NEGP
        EC=DLEMP(I)
        CALL FINDI(EITL,EC,NDATA,J)
        XP0(I)=EXP(A(J)+EC*(B(J)+EC*(C(J)+EC*D(J))))
      ENDDO
C
      DO I=1,NDATA
        FL(I)=LOG(XP1(I))
      ENDDO
      CALL SPLINE(EITL,FL,A,B,C,D,0.0D0,0.0D0,NDATA)
      DO I=1,NEGP
        EC=DLEMP(I)
        CALL FINDI(EITL,EC,NDATA,J)
        XP1(I)=EXP(A(J)+EC*(B(J)+EC*(C(J)+EC*D(J))))
      ENDDO
C
      DO I=1,NDATA
        FL(I)=LOG(XP2(I))
      ENDDO
      CALL SPLINE(EITL,FL,A,B,C,D,0.0D0,0.0D0,NDATA)
      DO I=1,NEGP
        EC=DLEMP(I)
        CALL FINDI(EITL,EC,NDATA,J)
        XP2(I)=EXP(A(J)+EC*(B(J)+EC*(C(J)+EC*D(J))))
      ENDDO
C
      DO I=1,NEGP
        XS0=XP0(I)
        XS1=XP1(I)
        XS2=XP2(I)
        FPEL=1.0D0/(XS0*VMOL(M))
        FPT1=1.0D0/(XS1*VMOL(M))
        FPST=ET(I)/(CSTPP(M,I)+RSTPP(M,I))
        XS0H=1.0D0/(VMOL(M)*MAX(FPEL,MIN(C1(M)*FPT1,C2(M)*FPST)))
        CALL EEL0(XS0,XS1,XS2,XS0H,AAA,BBB,RNDC,XS1S,XS2S)
        SPHEL(M,I)=XS0H*VMOL(M)
        RNDCP(M,I)=RNDC
        AP(M,I)=AAA
        BP(M,I)=BBB
        T1P0(I)=XS1S
        T1P(M,I)=T1PI(I)+XS1S*VMOL(M)
        T2P0(I)=XS2S
        T2P(M,I)=T2PI(I)+XS2S*VMOL(M)
      ENDDO
C
C  ****  Print positron elastic scattering tables.
C
      IF(INFO.GE.3) WRITE(IWR,1004)
 1004 FORMAT(/1X,'PENELOPE >>>  Elastic scattering of positrons')
      IF(INFO.GE.3) WRITE(IWR,1005)
 1005 FORMAT(/3X,'E (eV)',6X,'MFP (mtu)',3X,'TMFP1 (mtu)',2X,
     1  'MFPh (mtu)',8X,'A',11X,'B',11X,'RNDC',/1X,90('-'))
      DO I=1,NEGP
        FP0=RHO(M)/(XP0(I)*VMOL(M))
        FP1=RHO(M)/(XP1(I)*VMOL(M))
        HMFP=RHO(M)/SPHEL(M,I)
        IF(INFO.GE.3) WRITE(IWR,'(1P,7(E12.5,1X))') ET(I),FP0,FP1,
     1    HMFP,AP(M,I),BP(M,I),RNDCP(M,I)
        SPHEL(M,I)=LOG(SPHEL(M,I))
        AP(M,I)=LOG(AP(M,I))
C  ****  Soft scattering events are switched off when T1P is too small.
        IF(T1P(M,I).GT.1.0D-6*XP1(I)*VMOL(M)) THEN
          T1P(M,I)=LOG(MAX(T1P(M,I),1.0D-35))
          T2P(M,I)=LOG(MAX(T2P(M,I),1.0D-35))
        ELSE
          T1P(M,I)=-100.0D0
          T2P(M,I)=-100.0D0
        ENDIF
      ENDDO
      RETURN
      END
C  *********************************************************************
C                       SUBROUTINE EEL0
C  *********************************************************************
      SUBROUTINE EEL0(XS0,XS1,XS2,XS0H,A,B,RNDC,XS1S,XS2S)
C
C     This subroutine determines the parameters of the MW model for
C  elastic scattering of electrons and positrons and initializes the
C  mixed simulation algorithm (for particles with a given energy).
C
C  Input arguments:
C    XS0 .... total x-section (cm**2).
C    XS1 .... 1st transport x-section (cm**2).
C    XS2 .... 2nd transport x-section (cm**2).
C    XS0H ... suggested x-section for hard events (cm**2).
C
C  Output values:
C    A, B ... angular distribution parameters.
C    RNDC ... cutoff probability.
C    XS0H ... adopted x-section for hard events (cm**2).
C    XS1S ... 1st transport x-section for soft events (cm**2).
C    XS2S ... 2nd transport x-section for soft events (cm**2).
C
      use mtpar
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N)
C
      IF(XS0.LT.0.0D0) STOP 'EEL0. Negative total cross section.'
      RMU1=XS1/(2.0D0*XS0)
      RMU2=(3.0D0*XS1-XS2)/(6.0D0*XS0)
      IF(RMU1.LT.0.0D0.OR.RMU1.GT.1.0D0) THEN
        WRITE(26,1001) XS0,XS1
 1001   FORMAT(/3X,'*** The arguments in subroutine EEL0 are inconsi',
     1    'stent.'/7X,'XS0 = ',1P,E14.7,', XS1 = ',E14.7)
        STOP 'EEL0. Inconsistent arguments.'
      ENDIF
C
C  ****  Wentzel screening parameter.
C
      A=1.0D0
   10 A=A+A
      TST=A*(A+1.0D0)*LOG((1.0D0+A)/A)-A-RMU1
      IF(TST.LT.0.0D0) GO TO 10
      AU=A
      AL=0.0D0
    1 A=0.5D0*(AL+AU)
      TST=A*(A+1.0D0)*LOG((1.0D0+A)/A)-A-RMU1
      IF(TST.GT.0.0D0) THEN
        AU=A
      ELSE
        AL=A
      ENDIF
      IF(ABS(TST).GT.1.0D-15) GO TO 1
C  ****  At high energies, when truncation errors in the input tables
C  are significant, we use delta scattering.
      IF(RMU2-RMU1**2.LT.1.0D-12.OR.A.LT.1.0D-9) THEN
        B=1.0D0
        RNDC=1.0D0-XS0H/XS0
        IF(RNDC.LT.1.0D-14) RNDC=0.0D0
        XS1S=XS1*RNDC
        XS2S=XS2*RNDC
        RETURN
      ENDIF
C
      RMU1W=A*(A+1.0D0)*LOG((1.0D0+A)/A)-A
      RMU2W=A*(1.0D0-2.0D0*RMU1W)
      B=(RMU2W-RMU2)/(RMU2W-RMU1W*RMU1W)
C
C  ****  Case I.
C
      IF(B.GT.0.0D0) THEN
        RNDC=1.0D0-XS0H/XS0
        IF(RNDC.LT.1.0D-6) THEN
          RNDC=0.0D0
          XS0H=XS0
          XS1S=0.0D0
          XS2S=0.0D0
          RETURN
        ENDIF
C
        A1=A+1.0D0
        B1=1.0D0-B
        RND0=B1*A1*RMU1/(A+RMU1)
        RNDC=1.0D0-XS0H/XS0
        IF(RNDC.LT.RND0) THEN
          RMUC=RNDC*A/(B1*A1-RNDC)
          XS1S=B1*A*A1*(LOG((A+RMUC)/A)-(RMUC/(A+RMUC)))
          XS2S=B1*(A*A1*RMUC**2/(A+RMUC))-2.0D0*A*XS1S
        ELSE IF(RNDC.GT.RND0+B) THEN
          RNDMB=RNDC-B
          RMUC=RNDMB*A/(B1*A1-RNDMB)
          XS1S=B1*A*A1*(LOG((A+RMUC)/A)-(RMUC/(A+RMUC)))
          XS2S=B1*(A*A1*RMUC**2/(A+RMUC))-2.0D0*A*XS1S
          XS1S=XS1S+B*RMU1
          XS2S=XS2S+B*RMU1**2
        ELSE
          RMUC=RMU1
          WB=RNDC-RND0
          XS1S=B1*A*A1*(LOG((A+RMUC)/A)-(RMUC/(A+RMUC)))
          XS2S=B1*(A*A1*RMUC**2/(A+RMUC))-2.0D0*A*XS1S
          XS1S=XS1S+WB*RMU1
          XS2S=XS2S+WB*RMU1**2
        ENDIF
        XS2S=6.0D0*XS0*(XS1S-XS2S)
        XS1S=2.0D0*XS0*XS1S
        RETURN
      ENDIF
      IF(B.GT.-1.0D-12) THEN
        B=0.0D0
        RNDC=1.0D0-XS0H/XS0
        A1=A+1.0D0
        RMUC=RNDC*A/(A1-RNDC)
        XS1S=A*A1*(LOG((A+RMUC)/A)-(RMUC/(A+RMUC)))
        XS2S=(A*A1*RMUC**2/(A+RMUC))-2.0D0*A*XS1S
        XS2S=6.0D0*XS0*(XS1S-XS2S)
        XS1S=2.0D0*XS0*XS1S
        RETURN
      ENDIF
C
C  ****  Case II.
C
      C1=8.333333333333333D-1
      C2=7.083333333333333D-1
      D1=C2-RMU2
      D2=C1-RMU1
      D3=C2*RMU1-C1*RMU2
      AL=1.0D-24
      AU=A
    2 A=0.5D0*(AL+AU)
      RMU1W=A*(A+1.0D0)*LOG((1.0D0+A)/A)-A
      RMU2W=A*(1.0D0-2.0D0*RMU1W)
      F=D1*RMU1W-D2*RMU2W-D3
      IF(F.LT.0.0D0) THEN
        AL=A
      ELSE
        AU=A
      ENDIF
      IF(AU-AL.GT.1.0D-14*A) GO TO 2
      B=(RMU1W-RMU1)/(C1-RMU1W)
C
      RNDC=1.0D0-XS0H/XS0
      IF(RNDC.LT.1.0D-10) THEN
        RNDC=0.0D0
        XS0H=XS0
        XS1S=0.0D0
        XS2S=0.0D0
        RETURN
      ENDIF
      A1=A+1.0D0
      B1=1.0D0+B
      RNDCM=B1*A1*0.5D0/(A+0.5D0)
      IF(RNDC.GT.RNDCM) THEN
        RNDC=RNDCM
        XS0H=XS0*(1.0D0-RNDC)
      ENDIF
      RMUC=RNDC*A/(B1*A1-RNDC)
      XS1S=B1*A*A1*(LOG((A+RMUC)/A)-(RMUC/(A+RMUC)))
      XS2S=B1*(A*A1*RMUC**2/(A+RMUC))-2.0D0*A*XS1S
      XS2S=6.0D0*XS0*(XS1S-XS2S)
      XS1S=2.0D0*XS0*XS1S
      RETURN
      END
C  *********************************************************************
C                       SUBROUTINE EELW
C  *********************************************************************
      SUBROUTINE EELW(M,IWR)
C
C  This subroutine generates a table of integrated cross sections for
C  elastic scattering of electrons and positrons in material M, and
C  writes it on the material definition file. Data are read from the
C  files 'pdeelZZ.p05'.
C
      use mtpar
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N)
C
      CHARACTER*12 FILEN
      CHARACTER*1 LDIG(10),LDIG1,LDIG2
      DATA LDIG/'0','1','2','3','4','5','6','7','8','9'/
C  ****  Composition data.
      PARAMETER (MAXMAT=10)
      COMMON/COMPOS/STF(MAXMAT,30),ZT(MAXMAT),AT(MAXMAT),RHO(MAXMAT),
     1  VMOL(MAXMAT),IZ(MAXMAT,30),NELEM(MAXMAT)
C  ****  Total and transport cross sections.
      PARAMETER (NEGP=200)
C  ****  Total and transport cross sections.
      COMMON/CEEL00/EIT(NEGP),XE0(NEGP),XE1(NEGP),XE2(NEGP),XP0(NEGP),
     1  XP1(NEGP),XP2(NEGP),T1E0(NEGP),T2E0(NEGP),T1P0(NEGP),T2P0(NEGP),
     2  EITL(NEGP),FL(NEGP),A(NEGP),B(NEGP),C(NEGP),D(NEGP)
C
C  ****  Building the cross section table.
C
      DO I=1,NEGP
        XE0(I)=0.0D0
        XE1(I)=0.0D0
        XE2(I)=0.0D0
        XP0(I)=0.0D0
        XP1(I)=0.0D0
        XP2(I)=0.0D0
      ENDDO
C
      DO IEL=1,NELEM(M)
        IZZ=IZ(M,IEL)
        WGHT=STF(M,IEL)
        NLD=IZZ
        NLD1=NLD-10*(NLD/10)
        NLD2=(NLD-NLD1)/10
        LDIG1=LDIG(NLD1+1)
        LDIG2=LDIG(NLD2+1)
        FILEN='pdeel'//LDIG2//LDIG1//'.p05'
        OPEN(3,FILE=FILEN)
        READ(3,*) IZZZ
        IF(IZZZ.NE.IZZ) STOP 'EELW. Wrong file.'
        DO I=1,NEGP
          READ(3,*,END=1) EIT(I),XE0P,XE1P,XE2P,XP0P,XP1P,XP2P
          XE0(I)=XE0(I)+WGHT*XE0P
          XE1(I)=XE1(I)+WGHT*XE1P
          XE2(I)=XE2(I)+WGHT*XE2P
          XP0(I)=XP0(I)+WGHT*XP0P
          XP1(I)=XP1(I)+WGHT*XP1P
          XP2(I)=XP2(I)+WGHT*XP2P
          NPTAB=I
        ENDDO
    1   CONTINUE
        CLOSE(3)
      ENDDO
C
C  ****  Write final x-section table.
C
      WRITE(IWR,2001) NPTAB
 2001 FORMAT(1X,'***  Electron and positron elastic cross sections',
     1  ',  NDATA =',I4)
      DO I=1,NPTAB
        WRITE(IWR,2002) EIT(I),XE0(I),XE1(I),XE2(I),XP0(I),XP1(I),
     1    XP2(I)
      ENDDO
 2002 FORMAT(1P,E10.3,6E12.5)
      RETURN
      END
C  *********************************************************************
C                       SUBROUTINE EIN
C  *********************************************************************
      SUBROUTINE EIN(E,DELTA,DE,EP,CDT,ES,CDTS,M,IOSC)
C
C  Random sampling of hard inelastic collisions of electrons.
C
C  Sternheimer-Liljequist GOS model
C
C  Input arguments:
C    E ....... electron energy (eV).
C    M ....... material where electrons propagate.
C    DELTA ... Fermi's density effect correction.
C  Output arguments:
C    DE ...... energy loss (eV)
C    EP ...... energy of the scattered electron (eV)
C    CDT ..... cosine of the polar scattering angle.
C    ES ...... energy of the emitted secondary electron (eV)
C    CDTS .... polar cosine of direction of the secondary electron.
C    IOSC .... index of the oscillator that has been 'ionized'.
C
      use mtpar
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N)
      PARAMETER (REV=5.10998902D5)  ! Electron rest energy (eV)
      PARAMETER (TREV=2.0D0*REV)
C  ****  Simulation parameters.
      PARAMETER (MAXMAT=10)
      COMMON/CSIMPA/EABS(3,MAXMAT),C1(MAXMAT),C2(MAXMAT),WCC(MAXMAT),
     1  WCR(MAXMAT)
      COMMON/CECUTR/ECUTR(MAXMAT)
C  ****  Energy grid and interpolation constants for the current energy.
      PARAMETER (NEGP=200)
      COMMON/CEGRID/EL,EU,ET(NEGP),DLEMP(NEGP),DLEMP1,DLFC,
     1  XEL,XE,XEK,KE
C  ****  E/P inelastic collisions.
      PARAMETER (NO=64)
      COMMON/CEIN/EXPOT(MAXMAT),OP2(MAXMAT),F(MAXMAT,NO),UI(MAXMAT,NO),
     1  WRI(MAXMAT,NO),KZ(MAXMAT,NO),KS(MAXMAT,NO),NOSC(MAXMAT)
      COMMON/CEINAC/EINAC(MAXMAT,NEGP,NO),PINAC(MAXMAT,NEGP,NO)
C
!     ! EXTERNALRAND
C
      WCCM=WCC(M)
      IF(WCCM.GT.E) THEN
        DE=0.0D0
        EP=E
        CDT=1.0D0
        ES=0.0D0
        CDTS=0.0D0
        IOSC=NO
        RETURN
      ENDIF
C  ****  Energy grid point.
      PK=(XEL-DLEMP(KE))*DLFC
      status=rmt(r)
      IF(r.LT.PK) THEN
        JE=KE+1
      ELSE
        JE=KE
      ENDIF
C
C  ************  Selection of the active oscillator.
C
      status=rmt(r)
      TST=r
C  ****  Bipartition search.
      IOSC=1
      JO=NOSC(M)+1
    1 IT=(IOSC+JO)/2
      IF(TST.GT.EINAC(M,JE,IT)) THEN
        IOSC=IT
      ELSE
        JO=IT
      ENDIF
      IF(JO-IOSC.GT.1) GO TO 1
C
C  ****  Constants.
C
      RB=E+TREV
      GAM=1.0D0+E/REV
      GAM2=GAM*GAM
      BETA2=(GAM2-1.0D0)/GAM2
      AMOL=((GAM-1.0D0)/GAM)**2
C
C  ************  Partial cross sections of the active oscillator.
C
      WI=WRI(M,IOSC)
C  ****  Distant excitations.
      IF(WI.GT.WCCM.AND.WI.LT.E) THEN
        CPS=E*RB
        CP=SQRT(CPS)
        XHDT0=MAX(LOG(GAM2)-BETA2-DELTA,0.0D0)
        IF(WI.GT.1.0D-6*E) THEN
          CPP=SQRT((E-WI)*(E-WI+TREV))
          QM=SQRT((CP-CPP)**2+REV*REV)-REV
        ELSE
          QM=WI**2/(BETA2*TREV)
          QM=QM*(1.0D0-QM/TREV)
        ENDIF
        IF(QM.LT.WI) THEN
          XHDL=LOG(WI*(QM+TREV)/(QM*(WI+TREV)))/WI
          XHDT=XHDT0/WI
        ELSE
          QM=WI
          XHDL=0.0D0
          XHDT=0.0D0
        ENDIF
      ELSE
        QM=WI
        CPS=0.0D0
        CP=0.0D0
        XHDL=0.0D0
        XHDT=0.0D0
      ENDIF
C  ****  Close collisions.
      WMAXC=0.5D0*E
      WCL=MAX(WCCM,WI)
      RCL=WCL/E
      IF(WCL.LT.WMAXC) THEN
        RL1=1.0D0-RCL
        XHC=(AMOL*(0.5D0-RCL)+1.0D0/RCL-1.0D0/RL1
     1       +(1.0D0-AMOL)*LOG(RCL/RL1))/E
      ELSE
        XHC=0.0D0
      ENDIF
C
      XHTOT=XHC+XHDL+XHDT
      IF(XHTOT.LT.1.0D-35) THEN
        DE=0.0D0
        EP=E
        CDT=1.0D0
        ES=0.0D0
        CDTS=0.0D0
        IOSC=NO
        RETURN
      ENDIF
C
      status=rmt(r)
      TST=r*XHTOT
C
C  ************  Hard close collision.
C
      TS1=XHC
      IF(TST.LT.TS1) THEN
        A=5.0D0*AMOL
        ARCL=A*0.5D0*RCL
    2   CONTINUE
        status=rmt(r)
        FB=(1.0D0+ARCL)*r
        IF(FB.LT.1.0D0) THEN
          RK=RCL/(1.0D0-FB*(1.0D0-(RCL+RCL)))
        ELSE
          RK=RCL+(FB-1.0D0)*(0.5D0-RCL)/ARCL
        ENDIF
        RK2=RK*RK
        RKF=RK/(1.0D0-RK)
        PHI=1.0D0+RKF**2-RKF+AMOL*(RK2+RKF)
        status=rmt(r)
        IF(r*(1.0D0+A*RK2).GT.PHI) GO TO 2
C  ****  Energy and scattering angle (primary electron).
        DE=RK*E
        EP=E-DE
        CDT=SQRT(EP*RB/(E*(RB-DE)))
C  ****  Energy and emission angle of the delta ray.
        IF(KS(M,IOSC).LT.10) THEN
          IF(UI(M,IOSC).GT.ECUTR(M)) THEN
            ES=DE-UI(M,IOSC)
          ELSE
            ES=DE
          ENDIF
        ELSE
          ES=DE
        ENDIF
        CDTS=SQRT(DE*RB/(E*(DE+TREV)))
        RETURN
      ENDIF
C
C  ************  Hard distant longitudinal collision.
C
      TS1=TS1+XHDL
      DE=WI
      EP=E-DE
      IF(TST.LT.TS1) THEN
        QS=QM/(1.0D0+QM/TREV)
        status=rmt(r)
        Q=QS/(((QS/DE)*(1.0D0+DE/TREV))**r-(QS/TREV))
        QTREV=Q*(Q+TREV)
        CPPS=EP*(EP+TREV)
        CDT=(CPPS+CPS-QTREV)/(2.0D0*CP*SQRT(CPPS))
        IF(CDT.GT.1.0D0) CDT=1.0D0
C  ****  Energy and emission angle of the delta ray.
        IF(KS(M,IOSC).LT.10) THEN
          ES=DE-UI(M,IOSC)
        ELSE
          ES=DE
        ENDIF
        CDTS=0.5D0*(DE*(E+RB-DE)+QTREV)/SQRT(CPS*QTREV)
        IF(CDTS.GT.1.0D0) CDTS=1.0D0
        RETURN
      ENDIF
C
C  ****  Hard distant transverse collision.
C
      CDT=1.0D0
C  ****  Energy and emission angle of the delta ray.
      IF(KS(M,IOSC).LT.10) THEN
        IF(UI(M,IOSC).GT.ECUTR(M)) THEN
          ES=DE-UI(M,IOSC)
        ELSE
          ES=DE
        ENDIF
      ELSE
        ES=DE
      ENDIF
      CDTS=0.5D0
C
      RETURN
      END
C  *********************************************************************
C                       SUBROUTINE EINTX
C  *********************************************************************
      SUBROUTINE EINTX(E,WCC,XH0,XH1,XH2,XS0,XS1,XS2,XT1,XT2,DELTA,M)
C
C  Integrated cross sections for inelastic collisions of electrons of
C  energy E in material M, restricted to energy losses larger than and
C  less than the cutoff energy WCC.
C
C  Sternheimer-Liljequist GOS model.
C
C  Output arguments:
C    XH0 ... total cross section for hard colls. (cm**2).
C    XH1 ... stopping cross section for hard colls. (eV*cm**2)
C    XH2 ... straggling cross section for hard colls. (eV**2*cm**2).
C    XS0 ... total cross section for soft colls. (cm**2).
C    XS1 ... stopping cross section for soft colls. (eV*cm**2)
C    XS2 ... straggling cross section for soft colls. (eV**2*cm**2).
C    XT1 ... 1st transport cross section for soft colls. (cm**2).
C    XT2 ... 2nd transport cross section for soft colls. (cm**2).
C    DELTA ... Fermi's density effect correction.
C
      use mtpar
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N)
      PARAMETER (REV=5.10998902D5)  ! Electron rest energy (eV)
      PARAMETER (ELRAD=2.817940285D-13)  ! Class. electron radius (cm)
      PARAMETER (PI=3.1415926535897932D0, TREV=2.0D0*REV,
     1  PIELR2=PI*ELRAD*ELRAD)
      PARAMETER (MAXMAT=10)
C  ****  Composition data.
      COMMON/COMPOS/STF(MAXMAT,30),ZT(MAXMAT),AT(MAXMAT),RHO(MAXMAT),
     1  VMOL(MAXMAT),IZ(MAXMAT,30),NELEM(MAXMAT)
C  ****  E/P inelastic collisions.
      PARAMETER (NO=64)
      COMMON/CEIN/EXPOT(MAXMAT),OP2(MAXMAT),F(MAXMAT,NO),UI(MAXMAT,NO),
     1  WRI(MAXMAT,NO),KZ(MAXMAT,NO),KS(MAXMAT,NO),NOSC(MAXMAT)
C  ****  Partial cross sections of individual shells/oscillators.
      COMMON/CEIN00/SXH0(NO),SXH1(NO),SXH2(NO),SXS0(NO),SXS1(NO),
     1              SXS2(NO),SXT0(NO),SXT1(NO),SXT2(NO)
C
      COMMON/CEIN01/EE,BETA2,CP,CPS,XHDT0,AMOL,MOM
C
C  ****  Constants.
C
      GAM=1.0D0+E/REV
      GAM2=GAM*GAM
      BETA2=(GAM2-1.0D0)/GAM2
      CONST=PIELR2*TREV/BETA2
      XHDT0=LOG(GAM2)-BETA2
C
      CPS=E*(E+TREV)
      CP=SQRT(CPS)
      AMOL=(E/(E+REV))**2
C
C  ************  Density effect.
C
C  ****  Sternheimer's resonance energy (WL2=L**2).
      TST=ZT(M)/(GAM2*OP2(M))
      WL2=0.0D0
      FDEL=0.0D0
      DO I=1,NOSC(M)
        FDEL=FDEL+F(M,I)/(WRI(M,I)**2+WL2)
      ENDDO
      IF(FDEL.LT.TST) THEN
        DELTA=0.0D0
        GO TO 3
      ENDIF
      WL2=WRI(M,NOSC(M))**2
    1 WL2=WL2+WL2
      FDEL=0.0D0
      DO I=1,NOSC(M)
        FDEL=FDEL+F(M,I)/(WRI(M,I)**2+WL2)
      ENDDO
      IF(FDEL.GT.TST) GO TO 1
      WL2L=0.0D0
      WL2U=WL2
    2 WL2=0.5D0*(WL2L+WL2U)
      FDEL=0.0D0
      DO I=1,NOSC(M)
        FDEL=FDEL+F(M,I)/(WRI(M,I)**2+WL2)
      ENDDO
      IF(FDEL.GT.TST) THEN
        WL2L=WL2
      ELSE
        WL2U=WL2
      ENDIF
      IF(WL2U-WL2L.GT.1.0D-12*WL2) GO TO 2
C  ****  Density effect correction (delta).
      DELTA=0.0D0
      DO I=1,NOSC(M)
        DELTA=DELTA+F(M,I)*LOG(1.0D0+WL2/WRI(M,I)**2)
      ENDDO
      DELTA=(DELTA/ZT(M))-WL2/(GAM2*OP2(M))
    3 CONTINUE
C
C  ****  Shell-oscillator cross sections.
C
      DO I=1,NOSC(M)
        SXH0(I)=0.0D0
        SXH1(I)=0.0D0
        SXH2(I)=0.0D0
        SXS0(I)=0.0D0
        SXS1(I)=0.0D0
        SXS2(I)=0.0D0
        SXT0(I)=0.0D0
        SXT1(I)=0.0D0
        SXT2(I)=0.0D0
      ENDDO
      XH0=0.0D0
      XH1=0.0D0
      XH2=0.0D0
      XS0=0.0D0
      XS1=0.0D0
      XS2=0.0D0
      XT0=0.0D0
      XT1=0.0D0
      XT2=0.0D0
C
      DO I=1,NOSC(M)
        W=WRI(M,I)
        CALL EINTX1(E,W,DELTA,WCC,H0,H1,H2,S0,S1,S2,R0,R1,R2)
        SXH0(I)=F(M,I)*CONST*H0
        SXH1(I)=F(M,I)*CONST*H1
        SXH2(I)=F(M,I)*CONST*H2
        SXS0(I)=F(M,I)*CONST*S0
        SXS1(I)=F(M,I)*CONST*S1
        SXS2(I)=F(M,I)*CONST*S2
        SXT0(I)=F(M,I)*CONST*R0
        SXT1(I)=F(M,I)*CONST*2.0D0*R1
        SXT2(I)=F(M,I)*CONST*6.0D0*(R1-R2)
        XH0=XH0+SXH0(I)
        XH1=XH1+SXH1(I)
        XH2=XH2+SXH2(I)
        XS0=XS0+SXS0(I)
        XS1=XS1+SXS1(I)
        XS2=XS2+SXS2(I)
        XT0=XT0+SXT0(I)
        XT1=XT1+SXT1(I)
        XT2=XT2+SXT2(I)
      ENDDO
C
      RETURN
      END
C  *********************************************************************
C                       SUBROUTINE EINTX1
C  *********************************************************************
      SUBROUTINE EINTX1(E,WI,DELTA,WCC,H0,H1,H2,S0,S1,S2,R0,R1,R2)
C
C  Integrated cross sections for inelastic collisions of electrons with
C  a single-shell oscillator, restricted to energy losses larger than
C  WCC. The present subroutine is not self-consistent; some of the par-
C  ameters are precalculated by subroutine EINTX and entered through
C  the common block CEIN01.
C
      use mtpar
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N)
C  ****  Fundamental constants and related quantities.
      PARAMETER (REV=5.10998902D5)  ! Electron rest energy (eV)
      PARAMETER (TREV=2.0D0*REV)
C  ****  Coefficients transferred from subroutine EINTX.
      COMMON/CEIN01/EE,BETA2,CP,CPS,XHDT0,AMOL,MOM
      EXTERNAL EINSDX
C
      EE=E
      H0=0.0D0
      H1=0.0D0
      H2=0.0D0
      S0=0.0D0
      S1=0.0D0
      S2=0.0D0
      R0=0.0D0
      R1=0.0D0
      R2=0.0D0
      IF(E.LT.WI) RETURN
C
C  ****  Distant interactions.
C
      CP1S=(E-WI)*(E-WI+TREV)
      CP1=SQRT(CP1S)
      A=4.0D0*CP*CP1
      B=(CP-CP1)**2
      RMU1=(WI*(WI+TREV)-B)/A
C  ****  Distant longitudinal interactions.
      IF(WI.GT.1.0D-6*E) THEN
        QM=SQRT((CP-CP1)**2+REV**2)-REV
      ELSE
        QM=WI*WI/(BETA2*TREV)
        QM=QM*(1.0D0-QM/TREV)
      ENDIF
      IF(QM.LT.WI) THEN
        SDL1=LOG(WI*(QM+TREV)/(QM*(WI+TREV)))
      ELSE
        SDL1=0.0D0
      ENDIF
C  ****  Distant transverse interactions.
      IF(SDL1.GT.0.0D0) THEN
        SDT1=MAX(XHDT0-DELTA,0.0D0)
        SD1=SDL1+SDT1
        IF(WCC.GT.WI) THEN
          S1=SD1
          S0=SD1/WI
          S2=SD1*WI
C  ****  Soft distant transport moments of orders 0-2.
          BA=B/A
          R0=LOG((RMU1+BA)/BA)
          R1=RMU1-BA*R0
          R2=BA**2*R0+0.5D0*RMU1*(RMU1-2.0D0*BA)
          R0=R0/WI
          R1=R1/WI
          R2=R2/WI
        ELSE
          H1=SD1
          H0=SD1/WI
          H2=SD1*WI
        ENDIF
      ENDIF
C
C  ****  Close collisions (Moller's cross section).
C
      WL=MAX(WCC,WI)
      WU=0.5D0*E
      IF(WL.LT.WU-1.0D-5) THEN
        H0=H0+(1.0D0/(E-WU))-(1.0D0/(E-WL))
     1    -(1.0D0/WU)+(1.0D0/WL)
     2    +(1.0D0-AMOL)*LOG(((E-WU)*WL)/((E-WL)*WU))/E
     3    +AMOL*(WU-WL)/E**2
        H1=H1+LOG(WU/WL)+(E/(E-WU))-(E/(E-WL))
     1    +(2.0D0-AMOL)*LOG((E-WU)/(E-WL))
     2    +AMOL*(WU**2-WL**2)/(2.0D0*E**2)
        H2=H2+(2.0D0-AMOL)*(WU-WL)+(WU*(2.0D0*E-WU)/(E-WU))
     1    -(WL*(2.0D0*E-WL)/(E-WL))
     2    +(3.0D0-AMOL)*E*LOG((E-WU)/(E-WL))
     2    +AMOL*(WU**3-WL**3)/(3.0D0*E**2)
        WU=WL
      ENDIF
C
      WL=WI
      IF(WL.GT.WU-1.0D-5) RETURN
      S0=S0+(1.0D0/(E-WU))-(1.0D0/(E-WL))
     1  -(1.0D0/WU)+(1.0D0/WL)
     2  +(1.0D0-AMOL)*LOG(((E-WU)*WL)/((E-WL)*WU))/E
     3  +AMOL*(WU-WL)/E**2
      S1=S1+LOG(WU/WL)+(E/(E-WU))-(E/(E-WL))
     1  +(2.0D0-AMOL)*LOG((E-WU)/(E-WL))
     2  +AMOL*(WU**2-WL**2)/(2.0D0*E**2)
      S2=S2+(2.0D0-AMOL)*(WU-WL)+(WU*(2.0D0*E-WU)/(E-WU))
     1  -(WL*(2.0D0*E-WL)/(E-WL))
     2  +(3.0D0-AMOL)*E*LOG((E-WU)/(E-WL))
     2  +AMOL*(WU**3-WL**3)/(3.0D0*E**2)
C  ****  Soft close transport moments of orders 0-2.
      CP2S=(E-WU)*(E-WU+TREV)
      CP2=SQRT(CP2S)
      RMU2=(WU*(WU+TREV)-(CP-CP2)**2)/(4.0D0*CP*CP2)
      MOM=0
      R0=R0+SUMGA(EINSDX,RMU1,RMU2,1.0D-7)
      MOM=1
      R1=R1+SUMGA(EINSDX,RMU1,RMU2,1.0D-7)
      MOM=2
      R2=R2+SUMGA(EINSDX,RMU1,RMU2,1.0D-7)
C
      RETURN
      END
C  *********************************************************************
C                       FUNCTION EINSDX
C  *********************************************************************
      FUNCTION EINSDX(RMU)
C
C  Angular differential cross section for soft close inelastic colli-
C  sions of electrons.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N)
      PARAMETER (REV=5.10998902D5)  ! Electron rest energy (eV)
C  ****  Coefficients transferred from subroutine EINTX.
      COMMON/CEIN01/E,BETA2,CP,CPS,XHDT0,AMOL,MOM
C
      AUX=2.0D0*RMU*(1.0D0-RMU)
      DEN=E*AUX+REV
      W=CPS*AUX/DEN
      DWDMU=CPS*REV*(2.0D0-4.0D0*RMU)/DEN**2
      EINSDX=(1.0D0+(W/(E-W))**2-(1.0D0-AMOL)*(W/(E-W))+AMOL*(W/E)**2)
     1      *DWDMU*RMU**MOM/W**2
      RETURN
      END
C  *********************************************************************
C                       SUBROUTINE PIN
C  *********************************************************************
      SUBROUTINE PIN(E,DELTA,DE,EP,CDT,ES,CDTS,M,IOSC)
C
C  Random sampling of hard inelastic collisions of positrons.
C
C  Sternheimer-Liljequist GOS model
C
C  Input arguments:
C    E ....... positron energy (eV).
C    M ....... material where positrons propagate.
C    DELTA ... Fermi's density effect correction.
C  Output arguments:
C    DE ...... energy loss (eV)
C    EP ...... energy of the scattered positron (eV)
C    CDT ..... cosine of the polar scattering angle.
C    ES ...... energy of the emitted secondary electron (eV)
C    CDTS .... polar cosine of direction of the secondary electron.
C    IOSC .... index of the oscillator that has been 'ionized'.
C
      use mtpar
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N)
      PARAMETER (REV=5.10998902D5)  ! Electron rest energy (eV)
      PARAMETER (TREV=2.0D0*REV)
C  ****  Simulation parameters.
      PARAMETER (MAXMAT=10)
      COMMON/CSIMPA/EABS(3,MAXMAT),C1(MAXMAT),C2(MAXMAT),WCC(MAXMAT),
     1  WCR(MAXMAT)
      COMMON/CECUTR/ECUTR(MAXMAT)
C  ****  Energy grid and interpolation constants for the current energy.
      PARAMETER (NEGP=200)
      COMMON/CEGRID/EL,EU,ET(NEGP),DLEMP(NEGP),DLEMP1,DLFC,
     1  XEL,XE,XEK,KE
C  ****  E/P inelastic collisions.
      PARAMETER (NO=64)
      COMMON/CEIN/EXPOT(MAXMAT),OP2(MAXMAT),F(MAXMAT,NO),UI(MAXMAT,NO),
     1  WRI(MAXMAT,NO),KZ(MAXMAT,NO),KS(MAXMAT,NO),NOSC(MAXMAT)
      COMMON/CEINAC/EINAC(MAXMAT,NEGP,NO),PINAC(MAXMAT,NEGP,NO)
C
!     ! EXTERNALRAND
C
      WCCM=WCC(M)
      IF(WCCM.GT.E) THEN
        DE=0.0D0
        EP=E
        CDT=1.0D0
        ES=0.0D0
        CDTS=0.0D0
        IOSC=NO
        RETURN
      ENDIF
C  ****  Energy grid point.
      PK=(XEL-DLEMP(KE))*DLFC
      status=rmt(r)
      IF(r.LT.PK) THEN
        JE=KE+1
      ELSE
        JE=KE
      ENDIF
C
C  ************  Selection of the active oscillator.
C
      status=rmt(r)
      TST=r
C  ****  Bipartition search.
      IOSC=1
      JO=NOSC(M)+1
    1 IT=(IOSC+JO)/2
      IF(TST.GT.PINAC(M,JE,IT)) THEN
        IOSC=IT
      ELSE
        JO=IT
      ENDIF
      IF(JO-IOSC.GT.1) GO TO 1
C
C  ****  Constants.
C
      RB=E+TREV
      GAM=1.0D0+E/REV
      GAM2=GAM*GAM
      BETA2=(GAM2-1.0D0)/GAM2
      G12=(GAM+1.0D0)**2
      AMOL=(E/(E+REV))**2
      BHA1=AMOL*(2.0D0*G12-1.0D0)/(GAM2-1.0D0)
      BHA2=AMOL*(3.0D0+1.0D0/G12)
      BHA3=AMOL*2.0D0*GAM*(GAM-1.0D0)/G12
      BHA4=AMOL*(GAM-1.0D0)**2/G12
C
C  ************  Partial cross sections of the active oscillator.
C
      WI=WRI(M,IOSC)
C  ****  Distant excitations.
      IF(WI.GT.WCCM.AND.WI.LT.E) THEN
        CPS=E*RB
        CP=SQRT(CPS)
        XHDT0=MAX(LOG(GAM2)-BETA2-DELTA,0.0D0)
        IF(WI.GT.1.0D-6*E) THEN
          CPP=SQRT((E-WI)*(E-WI+TREV))
          QM=SQRT((CP-CPP)**2+REV*REV)-REV
        ELSE
          QM=WI**2/(BETA2*TREV)
          QM=QM*(1.0D0-QM/TREV)
        ENDIF
        IF(QM.LT.WI) THEN
          XHDL=LOG(WI*(QM+TREV)/(QM*(WI+TREV)))/WI
          XHDT=XHDT0/WI
        ELSE
          QM=WI
          XHDL=0.0D0
          XHDT=0.0D0
        ENDIF
      ELSE
        QM=WI
        CPS=0.0D0
        CP=0.0D0
        XHDL=0.0D0
        XHDT=0.0D0
      ENDIF
C  ****  Close collisions.
      WMAXC=E
      WCL=MAX(WCCM,WI)
      RCL=WCL/E
      IF(WCL.LT.WMAXC) THEN
        RL1=1.0D0-RCL
        XHC=((1.0D0/RCL-1.0D0)+BHA1*LOG(RCL)+BHA2*RL1
     1     +(BHA3/2.0D0)*(RCL**2-1.0D0)+(BHA4/3.0D0)*(1.0D0-RCL**3))/E
      ELSE
        XHC=0.0D0
      ENDIF
C
      XHTOT=XHC+XHDL+XHDT
      IF(XHTOT.LT.1.0D-35) THEN
        DE=0.0D0
        EP=E
        CDT=1.0D0
        ES=0.0D0
        CDTS=0.0D0
        IOSC=NO
        RETURN
      ENDIF
C
      status=rmt(r)
      TST=r*XHTOT
C
C  ************  Hard close collision.
C
      TS1=XHC
      IF(TST.LT.TS1) THEN
    2   CONTINUE
        status=rmt(r)
        RK=RCL/(1.0D0-r*(1.0D0-RCL))
        PHI=1.0D0-RK*(BHA1-RK*(BHA2-RK*(BHA3-BHA4*RK)))
      status=rmt(r)
        IF(r.GT.PHI) GO TO 2
C  ****  Energy and scattering angle (primary positron).
        DE=RK*E
        EP=E-DE
        CDT=SQRT(EP*RB/(E*(RB-DE)))
C  ****  Energy and emission angle of the delta ray.
        IF(KS(M,IOSC).LT.10) THEN
          IF(UI(M,IOSC).GT.ECUTR(M)) THEN
            ES=DE-UI(M,IOSC)
          ELSE
            ES=DE
          ENDIF
        ELSE
          ES=DE
        ENDIF
        CDTS=SQRT(DE*RB/(E*(DE+TREV)))
        RETURN
      ENDIF
C
C  ************  Hard distant longitudinal collision.
C
      TS1=TS1+XHDL
      DE=WI
      EP=E-DE
      IF(TST.LT.TS1) THEN
        QS=QM/(1.0D0+QM/TREV)
      status=rmt(r)
        Q=QS/(((QS/DE)*(1.0D0+DE/TREV))**r-(QS/TREV))
        QTREV=Q*(Q+TREV)
        CPPS=EP*(EP+TREV)
        CDT=(CPPS+CPS-QTREV)/(2.0D0*CP*SQRT(CPPS))
        IF(CDT.GT.1.0D0) CDT=1.0D0
C  ****  Energy and emission angle of the delta ray.
        IF(KS(M,IOSC).LT.10) THEN
          ES=DE-UI(M,IOSC)
        ELSE
          ES=DE
        ENDIF
        CDTS=0.5D0*(DE*(E+RB-DE)+QTREV)/SQRT(CPS*QTREV)
        IF(CDTS.GT.1.0D0) CDTS=1.0D0
        RETURN
      ENDIF
C
C  ****  Hard distant transverse collision.
C
      CDT=1.0D0
C  ****  Energy and emission angle of the delta ray.
      IF(KS(M,IOSC).LT.10) THEN
        IF(UI(M,IOSC).GT.ECUTR(M)) THEN
          ES=DE-UI(M,IOSC)
        ELSE
          ES=DE
        ENDIF
      ELSE
        ES=DE
      ENDIF
      CDTS=0.5D0
C
      RETURN
      END
C  *********************************************************************
C                       SUBROUTINE PINTX
C  *********************************************************************
      SUBROUTINE PINTX(E,WCC,XH0,XH1,XH2,XS0,XS1,XS2,XT1,XT2,DELTA,M)
C
C  Integrated cross sections for inelastic collisions of positrons of
C  energy E in material M, restricted to energy losses larger than and
C  less than the cutoff energy WCC.
C
C  Sternheimer-Liljequist GOS model.
C
C  Output arguments:
C    XH0 ... total cross section for hard colls. (cm**2).
C    XH1 ... stopping cross section for hard colls. (eV*cm**2)
C    XH2 ... straggling cross section for hard colls. (eV**2*cm**2).
C    XS0 ... total cross section for soft colls. (cm**2).
C    XS1 ... stopping cross section for soft colls. (eV*cm**2)
C    XS2 ... straggling cross section for soft colls. (eV**2*cm**2).
C    XT1 ... 1st transport cross section for soft colls. (cm**2).
C    XT2 ... 2nd transport cross section for soft colls. (cm**2).
C    DELTA ... Fermi's density effect correction.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N)
      PARAMETER (REV=5.10998902D5)  ! Electron rest energy (eV)
      PARAMETER (ELRAD=2.817940285D-13)  ! Class. electron radius (cm)
      PARAMETER (PI=3.1415926535897932D0, TREV=2.0D0*REV,
     1  PIELR2=PI*ELRAD*ELRAD)
      PARAMETER (MAXMAT=10)
C  ****  Composition data.
      COMMON/COMPOS/STF(MAXMAT,30),ZT(MAXMAT),AT(MAXMAT),RHO(MAXMAT),
     1  VMOL(MAXMAT),IZ(MAXMAT,30),NELEM(MAXMAT)
C  ****  E/P inelastic collisions.
      PARAMETER (NO=64)
      COMMON/CEIN/EXPOT(MAXMAT),OP2(MAXMAT),F(MAXMAT,NO),UI(MAXMAT,NO),
     1  WRI(MAXMAT,NO),KZ(MAXMAT,NO),KS(MAXMAT,NO),NOSC(MAXMAT)
C  ****  Partial cross sections of individual shells/oscillators.
      COMMON/CPIN00/SYH0(NO),SYH1(NO),SYH2(NO),SYS0(NO),SYS1(NO),
     1              SYS2(NO),SYT0(NO),SYT1(NO),SYT2(NO)
C
      COMMON/CPIN01/EE,BETA2,CP,CPS,XHDT0,BHA1,BHA2,BHA3,BHA4,MOM
C
C  ****  Constants.
C
      GAM=1.0D0+E/REV
      GAM2=GAM*GAM
      BETA2=(GAM2-1.0D0)/GAM2
      CONST=PIELR2*TREV/BETA2
      XHDT0=LOG(GAM2)-BETA2
C
      CPS=E*(E+TREV)
      CP=SQRT(CPS)
      G12=(GAM+1.0D0)**2
      AMOL=(E/(E+REV))**2
      BHA1=AMOL*(2.0D0*G12-1.0D0)/(GAM2-1.0D0)
      BHA2=AMOL*(3.0D0+1.0D0/G12)
      BHA3=AMOL*2.0D0*GAM*(GAM-1.0D0)/G12
      BHA4=AMOL*(GAM-1.0D0)**2/G12
C
C  ************  Density effect.
C
C  ****  Sternheimer's resonance energy (WL2=L**2).
      TST=ZT(M)/(GAM2*OP2(M))
      WL2=0.0D0
      FDEL=0.0D0
      DO I=1,NOSC(M)
        FDEL=FDEL+F(M,I)/(WRI(M,I)**2+WL2)
      ENDDO
      IF(FDEL.LT.TST) THEN
        DELTA=0.0D0
        GO TO 3
      ENDIF
      WL2=WRI(M,NOSC(M))**2
    1 WL2=WL2+WL2
      FDEL=0.0D0
      DO I=1,NOSC(M)
        FDEL=FDEL+F(M,I)/(WRI(M,I)**2+WL2)
      ENDDO
      IF(FDEL.GT.TST) GO TO 1
      WL2L=0.0D0
      WL2U=WL2
    2 WL2=0.5D0*(WL2L+WL2U)
      FDEL=0.0D0
      DO I=1,NOSC(M)
        FDEL=FDEL+F(M,I)/(WRI(M,I)**2+WL2)
      ENDDO
      IF(FDEL.GT.TST) THEN
        WL2L=WL2
      ELSE
        WL2U=WL2
      ENDIF
      IF(WL2U-WL2L.GT.1.0D-12*WL2) GO TO 2
C  ****  Density effect correction (delta).
      DELTA=0.0D0
      DO I=1,NOSC(M)
        DELTA=DELTA+F(M,I)*LOG(1.0D0+WL2/WRI(M,I)**2)
      ENDDO
      DELTA=(DELTA/ZT(M))-WL2/(GAM2*OP2(M))
    3 CONTINUE
C
C  ****  Shell-oscillator cross sections.
C
      DO I=1,NOSC(M)
        SYH0(I)=0.0D0
        SYH1(I)=0.0D0
        SYH2(I)=0.0D0
        SYS0(I)=0.0D0
        SYS1(I)=0.0D0
        SYS2(I)=0.0D0
        SYT0(I)=0.0D0
        SYT1(I)=0.0D0
        SYT2(I)=0.0D0
      ENDDO
      XH0=0.0D0
      XH1=0.0D0
      XH2=0.0D0
      XS0=0.0D0
      XS1=0.0D0
      XS2=0.0D0
      XT0=0.0D0
      XT1=0.0D0
      XT2=0.0D0
C
      DO I=1,NOSC(M)
        W=WRI(M,I)
        CALL PINTX1(E,W,DELTA,WCC,H0,H1,H2,S0,S1,S2,R0,R1,R2)
        SYH0(I)=F(M,I)*CONST*H0
        SYH1(I)=F(M,I)*CONST*H1
        SYH2(I)=F(M,I)*CONST*H2
        SYS0(I)=F(M,I)*CONST*S0
        SYS1(I)=F(M,I)*CONST*S1
        SYS2(I)=F(M,I)*CONST*S2
        SYT0(I)=F(M,I)*CONST*R0
        SYT1(I)=F(M,I)*CONST*2.0D0*R1
        SYT2(I)=F(M,I)*CONST*6.0D0*(R1-R2)
        XH0=XH0+SYH0(I)
        XH1=XH1+SYH1(I)
        XH2=XH2+SYH2(I)
        XS0=XS0+SYS0(I)
        XS1=XS1+SYS1(I)
        XS2=XS2+SYS2(I)
        XT0=XT0+SYT0(I)
        XT1=XT1+SYT1(I)
        XT2=XT2+SYT2(I)
      ENDDO
C
      RETURN
      END
C  *********************************************************************
C                       SUBROUTINE PINTX1
C  *********************************************************************
      SUBROUTINE PINTX1(E,WI,DELTA,WCC,H0,H1,H2,S0,S1,S2,R0,R1,R2)
C
C  Integrated cross sections for inelastic collisions of positrons with
C  a single-shell oscillator, restricted to energy losses larger than
C  WCC. The present subroutine is not self-consistent; some of the par-
C  ameters are precalculated by subroutine PINTX and entered through
C  the common block CPIN01.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N)
C  ****  Fundamental constants and related quantities.
      PARAMETER (REV=5.10998902D5)  ! Electron rest energy (eV)
      PARAMETER (TREV=2.0D0*REV)
C  ****  Coefficients transferred from subroutine EINTX.
      COMMON/CPIN01/EE,BETA2,CP,CPS,XHDT0,BHA1,BHA2,BHA3,BHA4,MOM
      EXTERNAL PINSDX
C
      EE=E
      H0=0.0D0
      H1=0.0D0
      H2=0.0D0
      S0=0.0D0
      S1=0.0D0
      S2=0.0D0
      R0=0.0D0
      R1=0.0D0
      R2=0.0D0
      IF(E.LT.WI) RETURN
C
C  ****  Distant interactions.
C
      CP1S=(E-WI)*(E-WI+TREV)
      CP1=SQRT(CP1S)
      A=4.0D0*CP*CP1
      B=(CP-CP1)**2
      RMU1=(WI*(WI+TREV)-B)/A
C  ****  Distant longitudinal interactions.
      IF(WI.GT.1.0D-6*E) THEN
        QM=SQRT((CP-CP1)**2+REV**2)-REV
      ELSE
        QM=WI*WI/(BETA2*TREV)
        QM=QM*(1.0D0-QM/TREV)
      ENDIF
      IF(QM.LT.WI) THEN
        SDL1=LOG(WI*(QM+TREV)/(QM*(WI+TREV)))
      ELSE
        SDL1=0.0D0
      ENDIF
C  ****  Distant transverse interactions.
      IF(SDL1.GT.0.0D0) THEN
        SDT1=MAX(XHDT0-DELTA,0.0D0)
        SD1=SDL1+SDT1
        IF(WCC.GT.WI) THEN
          S1=SD1
          S0=SD1/WI
          S2=SD1*WI
C  ****  Soft distant transport moments of orders 0-2.
          BA=B/A
          R0=LOG((RMU1+BA)/BA)
          R1=RMU1-BA*R0
          R2=BA**2*R0+0.5D0*RMU1*(RMU1-2.0D0*BA)
          R0=R0/WI
          R1=R1/WI
          R2=R2/WI
        ELSE
          H1=SD1
          H0=SD1/WI
          H2=SD1*WI
        ENDIF
      ENDIF
C
C  ****  Close collisions (Bhabha's cross section).
C
      WL=MAX(WCC,WI)
      WU=E
      IF(WL.LT.WU-1.0D-5) THEN
        H0=H0+(1.0D0/WL)-(1.0D0/WU)-BHA1*LOG(WU/WL)/E
     1    +BHA2*(WU-WL)/E**2-BHA3*(WU**2-WL**2)/(2.0D0*E**3)
     2    +BHA4*(WU**3-WL**3)/(3.0D0*E**4)
        H1=H1+LOG(WU/WL)-BHA1*(WU-WL)/E
     1    +BHA2*(WU**2-WL**2)/(2.0D0*E**2)
     2    -BHA3*(WU**3-WL**3)/(3.0D0*E**3)
     3    +BHA4*(WU**4-WL**4)/(4.0D0*E**4)
        H2=H2+WU-WL-BHA1*(WU**2-WL**2)/(2.0D0*E)
     1    +BHA2*(WU**3-WL**3)/(3.0D0*E**2)
     2    -BHA3*(WU**4-WL**4)/(4.0D0*E**3)
     3    +BHA4*(WU**5-WL**5)/(5.0D0*E**4)
        WU=WL
      ENDIF
C
      WL=WI
      IF(WL.GT.WU-1.0D-5) RETURN
        S0=S0+(1.0D0/WL)-(1.0D0/WU)-BHA1*LOG(WU/WL)/E
     1    +BHA2*(WU-WL)/E**2-BHA3*(WU**2-WL**2)/(2.0D0*E**3)
     2    +BHA4*(WU**3-WL**3)/(3.0D0*E**4)
        S1=S1+LOG(WU/WL)-BHA1*(WU-WL)/E
     1    +BHA2*(WU**2-WL**2)/(2.0D0*E**2)
     2    -BHA3*(WU**3-WL**3)/(3.0D0*E**3)
     3    +BHA4*(WU**4-WL**4)/(4.0D0*E**4)
        S2=S2+WU-WL-BHA1*(WU**2-WL**2)/(2.0D0*E)
     1    +BHA2*(WU**3-WL**3)/(3.0D0*E**2)
     2    -BHA3*(WU**4-WL**4)/(4.0D0*E**3)
     3    +BHA4*(WU**5-WL**5)/(5.0D0*E**4)
C  ****  Soft close transport moments of orders 0-2.
      IF(WU.LT.E-1.0D0) THEN
        CP2S=(E-WU)*(E-WU+TREV)
        CP2=SQRT(CP2S)
        RMU2=(WU*(WU+TREV)-(CP-CP2)**2)/(4.0D0*CP*CP2)
      ELSE
        RMU2=0.5D0
      ENDIF
      MOM=0
      R0=R0+SUMGA(PINSDX,RMU1,RMU2,1.0D-7)
      MOM=1
      R1=R1+SUMGA(PINSDX,RMU1,RMU2,1.0D-7)
      MOM=2
      R2=R2+SUMGA(PINSDX,RMU1,RMU2,1.0D-7)
C
      RETURN
      END
C  *********************************************************************
C                       FUNCTION PINSDX
C  *********************************************************************
      FUNCTION PINSDX(RMU)
C
C  Angular differential cross section for soft close inelastic colli-
C  sions of positrons.
C
      use mtpar
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N)
      PARAMETER (REV=5.10998902D5)  ! Electron rest energy (eV)
C  ****  Coefficients transferred from subroutine EINTX.
      COMMON/CPIN01/E,BETA2,CP,CPS,XHDT0,BHA1,BHA2,BHA3,BHA4,MOM
C
      AUX=2.0D0*RMU*(1.0D0-RMU)
      DEN=E*AUX+REV
      W=CPS*AUX/DEN
      DWDMU=CPS*REV*(2.0D0-4.0D0*RMU)/DEN**2
      WE=W/E
      PINSDX=(1.0D0-WE*(BHA1-WE*(BHA2-WE*(BHA3-WE*BHA4))))
     1      *DWDMU*RMU**MOM/W**2
      RETURN
      END
C  *********************************************************************
C                       SUBROUTINE ESI
C  *********************************************************************
      SUBROUTINE ESI(IZZ,ISH)
C
C  Ionization of inner shells by impact of electrons.
C
C  Output arguments:
C    IZZ .... atomic number of the element where ionization has ocurred.
C    ISH .... atomic electron shell that has been ionized.
C
      use mtpar
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N)
      CHARACTER*2 LASYMB
C  ****  Main-PENELOPE common.
      COMMON/TRACK/E,X,Y,Z,U,V,W,WGHT,KPAR,IBODY,M,ILB(5)
C  ****  Simulation parameters.
      PARAMETER (MAXMAT=10)
      COMMON/CSIMPA/EABS(3,MAXMAT),C1(MAXMAT),C2(MAXMAT),WCC(MAXMAT),
     1  WCR(MAXMAT)
      COMMON/CECUTR/ECUTR(MAXMAT)
C  ****  Energy grid and interpolation constants for the current energy.
      PARAMETER (NEGP=200)
      COMMON/CEGRID/EL,EU,ET(NEGP),DLEMP(NEGP),DLEMP1,DLFC,
     1  XEL,XE,XEK,KE
C  ****  Element data.
      COMMON/CADATA/ATW(99),EPX(99),RA1(99),RA2(99),RA3(99),RA4(99),
     1  RA5(99),RSCR(99),ETA(99),EB(99,30),IFI(99,30),IKS(99,30),
     2  NSHT(99),LASYMB(99)
C  ****  Composition data.
      COMMON/COMPOS/STF(MAXMAT,30),ZT(MAXMAT),AT(MAXMAT),RHO(MAXMAT),
     1  VMOL(MAXMAT),IZ(MAXMAT,30),NELEM(MAXMAT)
C  ****  Inner shell ionization by electron impact.
      PARAMETER (NRP=6000)
      COMMON/CESI0/XESI(NRP,9),IESIF(99),IESIL(99),NSESI(99),NCURE
      DIMENSION PACSI(120),IZSI(120),ISHSI(120)
C
      ! EXTERNALRAND
C
      KT=1
      SEIST=0.0D0
      PACSI(1)=0.0D0
C
      DO J=1,NELEM(M)
        IZZ=IZ(M,J)
        INDC=IESIF(IZZ)-1
        DO ISH=1,NSESI(IZZ)
          WCUT=EB(IZZ,ISH)
          IF(WCUT.GT.ECUTR(M).AND.WCUT.LT.E) THEN
            PCSI=EXP(XESI(INDC+KE,ISH)
     1          +(XESI(INDC+KE+1,ISH)-XESI(INDC+KE,ISH))*XEK)
            IF(PCSI.GT.1.1D-35) THEN
              SEIST=SEIST+PCSI*STF(M,J)
              IZSI(KT)=IZZ
              ISHSI(KT)=ISH
              KT=KT+1
              PACSI(KT)=SEIST
            ENDIF
          ENDIF
        ENDDO
      ENDDO
C
      IF(KT.EQ.1) THEN
        IZZ=0
        ISH=0
        RETURN
      ENDIF
C
      status=rmt(r)
      TST=PACSI(KT)*r
C  ****  Bipartition search.
      IS=1
      JS=KT
    1 IT=(IS+JS)/2
      IF(TST.GT.PACSI(IT)) IS=IT
      IF(TST.LE.PACSI(IT)) JS=IT
      IF(JS-IS.GT.1) GO TO 1
C
      IZZ=IZSI(IS)
      ISH=ISHSI(IS)
      RETURN
      END
C  *********************************************************************
C                       SUBROUTINE ESIR
C  *********************************************************************
      SUBROUTINE ESIR(M,IRD,IWR,INFO)
C
C  This subroutine reads cross sections for inner-shell ionization by
C  electron impact of the elements in material M and prepares simulation
C  tables.
C
      use mtpar
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N)
      CHARACTER*2 LASYMB
C  ****  Energy grid and interpolation constants for the current energy.
      PARAMETER (NEGP=200)
      COMMON/CEGRID/EL,EU,ET(NEGP),DLEMP(NEGP),DLEMP1,DLFC,
     1  XEL,XE,XEK,KE
C  ****  Composition data.
      PARAMETER (MAXMAT=10)
      COMMON/COMPOS/STF(MAXMAT,30),ZT(MAXMAT),AT(MAXMAT),RHO(MAXMAT),
     1  VMOL(MAXMAT),IZ(MAXMAT,30),NELEM(MAXMAT)
      COMMON/CECUTR/ECUTR(MAXMAT)
C  ****  Element data.
      COMMON/CADATA/ATW(99),EPX(99),RA1(99),RA2(99),RA3(99),RA4(99),
     1  RA5(99),RSCR(99),ETA(99),EB(99,30),IFI(99,30),IKS(99,30),
     2  NSHT(99),LASYMB(99)
C  ****  Electron impact ionization cross section tables.
      PARAMETER (NRP=6000)
      COMMON/CESI0/XESI(NRP,9),IESIF(99),IESIL(99),NSESI(99),NCURE
      PARAMETER (NDIM=1000)
      DIMENSION E(NDIM),XESIR(NDIM,9),X(NDIM),Y(NDIM)
C
C  ************  Read element x-section tables
C
      DO IEL=1,NELEM(M)
        READ(IRD,1001) IZZ,NSHR,NDATA
 1001   FORMAT(47X,I3,11X,I3,10X,I4)
        IF(INFO.GE.2) WRITE(IWR,2001) IZZ,NSHR,NDATA
 2001   FORMAT(/1X,'***  Electron impact ionization cross sections, ',
     1    ' IZ =',I3,',  NSHELL =',I3,',  NDATA =',I4)
        IF(IZZ.NE.IZ(M,IEL)) STOP 'ESIR. Corrupt material data file.'
        IF(NDATA.GT.NDIM) STOP 'ESIR. Too many data points.'
        IF(NSHR.GT.9) STOP 'ESIR. Too many shells.'
        DO IE=1,NDATA
          READ(IRD,*) E(IE),(XESIR(IE,IS),IS=1,NSHR)
        ENDDO
C
C  ****  Remove shells with ionization energies less than 100 eV.
C
        IF(NSHR.GT.1) THEN
          NSHA=NSHR
          DO IS=NSHA,1,-1
            IF(EB(IZZ,IS).LT.100.0D0) THEN
              NSHR=NSHR-1
            ELSE
              GO TO 1
            ENDIF
          ENDDO
          IF(NSHR.LT.1) NSHR=1
        ENDIF
    1   CONTINUE
C
        IF(INFO.GE.2) THEN
          IF(NSHR.EQ.1) THEN
            WRITE(IWR,2102)
          ELSE IF(NSHR.EQ.2) THEN
            WRITE(IWR,2202)
          ELSE IF(NSHR.EQ.3) THEN
            WRITE(IWR,2302)
          ELSE IF(NSHR.EQ.4) THEN
            WRITE(IWR,2402)
          ELSE IF(NSHR.EQ.5) THEN
            WRITE(IWR,2502)
          ELSE IF(NSHR.EQ.6) THEN
            WRITE(IWR,2602)
          ELSE IF(NSHR.EQ.7) THEN
            WRITE(IWR,2702)
          ELSE IF(NSHR.EQ.8) THEN
            WRITE(IWR,2802)
          ELSE
            WRITE(IWR,2902)
          ENDIF
          DO IE=1,NDATA
            WRITE(IWR,'(1P,12E12.5)') E(IE),(XESIR(IE,IS),IS=1,NSHR)
          ENDDO
        ENDIF
 2102   FORMAT(/3X,'Energy',7X,'CS-K',
     1    /4X,'(eV)',2X,5X,'(cm**2)',/1X,24('-'))
 2202   FORMAT(/3X,'Energy',7X,'CS-K',8X,'CS-L1',
     1    /4X,'(eV)',2X,2(5X,'(cm**2)'),/1X,36('-'))
 2302   FORMAT(/3X,'Energy',7X,'CS-K',8X,'CS-L1',7X,
     1    'CS-L2',/4X,'(eV)',2X,3(5X,'(cm**2)'),/1X,48('-'))
 2402   FORMAT(/3X,'Energy',7X,'CS-K',8X,'CS-L1',7X,
     1    'CS-L2',7X,'CS-L3',
     2    /4X,'(eV)',2X,4(5X,'(cm**2)'),/1X,60('-'))
 2502   FORMAT(/3X,'Energy',7X,'CS-K',8X,'CS-L1',7X,
     1    'CS-L2',7X,'CS-L3',7X,'CS-M1',
     2    /4X,'(eV)',2X,5(5X,'(cm**2)'),/1X,72('-'))
 2602   FORMAT(/3X,'Energy',7X,'CS-K',8X,'CS-L1',7X,
     1    'CS-L2',7X,'CS-L3',7X,'CS-M1',7X,'CS-M2',
     2    /4X,'(eV)',2X,6(5X,'(cm**2)'),/1X,84('-'))
 2702   FORMAT(/3X,'Energy',7X,'CS-K',8X,'CS-L1',7X,
     1    'CS-L2',7X,'CS-L3',7X,'CS-M1',7X,'CS-M2',7X,'CS-M3',
     2    /4X,'(eV)',2X,7(5X,'(cm**2)'),/1X,96('-'))
 2802   FORMAT(/3X,'Energy',7X,'CS-K',8X,'CS-L1',7X,
     1    'CS-L2',7X,'CS-L3',7X,'CS-M1',7X,'CS-M2',7X,'CS-M3',7X,
     2    'CS-M4',/4X,'(eV)',2X,8(5X,'(cm**2)'),
     3    /1X,108('-'))
 2902   FORMAT(/3X,'Energy',7X,'CS-K',8X,'CS-L1',7X,
     1    'CS-L2',7X,'CS-L3',7X,'CS-M1',7X,'CS-M2',7X,'CS-M3',7X,
     2    'CS-M4',7X,'CS-M5',/4X,'(eV)',2X,9(5X,'(cm**2)'),
     3    /1X,120('-'))
C
        NSESI(IZZ)=NSHR
        IF(IESIF(IZZ).EQ.0) THEN
          IESIF(IZZ)=NCURE+1
          IF(NCURE+NEGP.GT.NRP) THEN
            WRITE(IWR,*) 'Insufficient memory storage in ESIR.'
            WRITE(IWR,*) 'Increase the value of the parameter NRP to',
     1        NCURE+NEGP
            WRITE(26,*) 'Insufficient memory storage in ESIR.'
            WRITE(26,*) 'Increase the value of the parameter NRP to',
     1        NCURE+NEGP
            STOP 'ESIR. Insufficient memory storage.'
          ENDIF
          DO IS=1,NSHR
            N=0
            DO I=1,NDATA
              IF(XESIR(I,IS).GT.1.0D-35) THEN
                N=N+1
                X(N)=LOG(E(I))
                IF(N.GT.1) X(N)=MAX(X(N),X(N-1)+1.0D-6)
                Y(N)=LOG(XESIR(I,IS))
              ENDIF
            ENDDO
            IF(N.GT.4) THEN
              DO I=1,NEGP
                IC=NCURE+I
                XC=DLEMP(I)
                IF(XC.GT.X(1)) THEN
                  CALL FINDI(X,XC,N,J)
                  IF(J.EQ.N) J=N-1
                  DX=X(J+1)-X(J)
                  IF(DX.GT.1.0D-6) THEN
                    XESI(IC,IS)=Y(J)+(XC-X(J))*(Y(J+1)-Y(J))/DX
                  ELSE
                    XESI(IC,IS)=(Y(J+1)+Y(J))/2.0D0
                  ENDIF
                ELSE
                  XESI(IC,IS)=-80.6D0
                ENDIF
              ENDDO
            ELSE
              DO I=1,NEGP
                IC=NCURE+I
                XESI(IC,IS)=-80.6D0
              ENDDO
            ENDIF
          ENDDO
          NCURE=NCURE+NEGP
          IESIL(IZZ)=NCURE
        ENDIF
      ENDDO
C
      RETURN
      END
C  *********************************************************************
C                       SUBROUTINE ESI0
C  *********************************************************************
      SUBROUTINE ESI0
C
C  This subroutine sets all variables in common /CESI0/ to zero.
C
C  It has to be invoked before reading the first material definition
C  file.
C
      use mtpar
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N)
C  ****  Ionization cross section tables.
      PARAMETER (NRP=6000)
      COMMON/CESI0/XESI(NRP,9),IESIF(99),IESIL(99),NSESI(99),NCURE
C
      DO I=1,99
        IESIF(I)=0
        IESIL(I)=0
        NSESI(I)=0
      ENDDO
C
      DO I=1,NRP
        DO J=1,9
          XESI(I,J)=-80.6D0
        ENDDO
      ENDDO
      NCURE=0
      RETURN
      END
C  *********************************************************************
C                       SUBROUTINE ESIW
C  *********************************************************************
      SUBROUTINE ESIW(M,IWR)
C
C  This subroutine generates tables of cross sections for inner-shell
C  ionization by electron impact for the elements in material M and
C  writes them on the material data file.
C
C  Data are read from the files 'pdesiZZ.p05'.
C
      use mtpar
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N)
C
      CHARACTER*12 FILEN
      CHARACTER*1 LDIG(10),LDIG1,LDIG2
      DATA LDIG/'0','1','2','3','4','5','6','7','8','9'/
C  ****  Composition data.
      PARAMETER (MAXMAT=10)
      COMMON/COMPOS/STF(MAXMAT,30),ZT(MAXMAT),AT(MAXMAT),RHO(MAXMAT),
     1  VMOL(MAXMAT),IZ(MAXMAT,30),NELEM(MAXMAT)
C
      PARAMETER (NES=400)
      DIMENSION E(NES),XESIR(NES,9)
C
      DO IEL=1,NELEM(M)
        IZZ=IZ(M,IEL)
        NLD=IZZ
        NLD1=NLD-10*(NLD/10)
        NLD2=(NLD-NLD1)/10
        LDIG1=LDIG(NLD1+1)
        LDIG2=LDIG(NLD2+1)
        FILEN='pdesi'//LDIG2//LDIG1//'.p05'
        OPEN(3,FILE=FILEN)
        READ(3,*) IZZZ,NSHR
        IF(IZZZ.NE.IZZ) STOP 'ESIW. Corrupt data file.'
        IF(NSHR.GT.9) STOP 'ESIW. Too many shells.'
        DO IE=1,NES
          READ(3,*,END=1) E(IE),(XESIR(IE,IS),IS=1,NSHR)
          NPTAB=IE
          IF(E(IE).GT.0.999D9) GO TO 1
        ENDDO
    1   CONTINUE
        CLOSE(3)
        WRITE(IWR,2001) IZZ,NSHR,NPTAB
 2001 FORMAT(1X,'***  Electron ionization cross sections,  IZ =',I3,
     1  ',  NSHELL =',I3,',  NDATA =',I4)
        DO IE=1,NPTAB
          WRITE(IWR,2002) E(IE),(XESIR(IE,IS),IS=1,NSHR)
        ENDDO
 2002 FORMAT(1P,10E12.5)
      ENDDO
      RETURN
      END
C  *********************************************************************
C                       SUBROUTINE PSI
C  *********************************************************************
      SUBROUTINE PSI(IZZ,ISH)
C
C  Ionization of inner shells by impact of positrons.
C
C  Output arguments:
C    IZZ .... atomic number of the element where ionization has ocurred.
C    ISH .... atomic electron shell that has been ionized.
C
      use mtpar
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N)
      CHARACTER*2 LASYMB
C  ****  Main-PENELOPE common.
      COMMON/TRACK/E,X,Y,Z,U,V,W,WGHT,KPAR,IBODY,M,ILB(5)
C  ****  Simulation parameters.
      PARAMETER (MAXMAT=10)
      COMMON/CSIMPA/EABS(3,MAXMAT),C1(MAXMAT),C2(MAXMAT),WCC(MAXMAT),
     1  WCR(MAXMAT)
      COMMON/CECUTR/ECUTR(MAXMAT)
C  ****  Energy grid and interpolation constants for the current energy.
      PARAMETER (NEGP=200)
      COMMON/CEGRID/EL,EU,ET(NEGP),DLEMP(NEGP),DLEMP1,DLFC,
     1  XEL,XE,XEK,KE
C  ****  Element data.
      COMMON/CADATA/ATW(99),EPX(99),RA1(99),RA2(99),RA3(99),RA4(99),
     1  RA5(99),RSCR(99),ETA(99),EB(99,30),IFI(99,30),IKS(99,30),
     2  NSHT(99),LASYMB(99)
C  ****  Composition data.
      COMMON/COMPOS/STF(MAXMAT,30),ZT(MAXMAT),AT(MAXMAT),RHO(MAXMAT),
     1  VMOL(MAXMAT),IZ(MAXMAT,30),NELEM(MAXMAT)
C  ****  Inner shell ionization by positron impact.
      PARAMETER (NRP=6000)
      COMMON/CPSI0/XPSI(NRP,9),IPSIF(99),IPSIL(99),NSPSI(99),NCURP
      DIMENSION PACSI(120),IZSI(120),ISHSI(120)
C
      ! EXTERNALRAND
C
      KT=1
      SPIST=0.0D0
      PACSI(1)=0.0D0
C
      DO J=1,NELEM(M)
        IZZ=IZ(M,J)
        INDC=IPSIF(IZZ)-1
        DO ISH=1,NSPSI(IZZ)
          WCUT=EB(IZZ,ISH)
          IF(WCUT.GT.ECUTR(M).AND.WCUT.LT.E) THEN
            PCSI=EXP(XPSI(INDC+KE,ISH)
     1          +(XPSI(INDC+KE+1,ISH)-XPSI(INDC+KE,ISH))*XEK)
            IF(PCSI.GT.1.1D-35) THEN
              SPIST=SPIST+PCSI*STF(M,J)
              IZSI(KT)=IZZ
              ISHSI(KT)=ISH
              KT=KT+1
              PACSI(KT)=SPIST
            ENDIF
          ENDIF
        ENDDO
      ENDDO
C
      IF(KT.EQ.1) THEN
        IZZ=0
        ISH=0
        RETURN
      ENDIF
C
        status=rmt(r)
      TST=PACSI(KT)*r
C  ****  Bipartition search.
      IS=1
      JS=KT
    1 IT=(IS+JS)/2
      IF(TST.GT.PACSI(IT)) IS=IT
      IF(TST.LE.PACSI(IT)) JS=IT
      IF(JS-IS.GT.1) GO TO 1
C
      IZZ=IZSI(IS)
      ISH=ISHSI(IS)
      RETURN
      END
C  *********************************************************************
C                       SUBROUTINE PSIR
C  *********************************************************************
      SUBROUTINE PSIR(M,IRD,IWR,INFO)
C
C  This subroutine reads cross sections for inner-shell ionization by
C  positron impact of the elements in material M and prepares simulation
C  tables.
C
      use mtpar 
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N)
      CHARACTER*2 LASYMB
C  ****  Energy grid and interpolation constants for the current energy.
      PARAMETER (NEGP=200)
      COMMON/CEGRID/EL,EU,ET(NEGP),DLEMP(NEGP),DLEMP1,DLFC,
     1  XEL,XE,XEK,KE
C  ****  Composition data.
      PARAMETER (MAXMAT=10)
      COMMON/COMPOS/STF(MAXMAT,30),ZT(MAXMAT),AT(MAXMAT),RHO(MAXMAT),
     1  VMOL(MAXMAT),IZ(MAXMAT,30),NELEM(MAXMAT)
      COMMON/CECUTR/ECUTR(MAXMAT)
C  ****  Element data.
      COMMON/CADATA/ATW(99),EPX(99),RA1(99),RA2(99),RA3(99),RA4(99),
     1  RA5(99),RSCR(99),ETA(99),EB(99,30),IFI(99,30),IKS(99,30),
     2  NSHT(99),LASYMB(99)
C  ****  Positron impact ionization cross section tables.
      PARAMETER (NRP=6000)
      COMMON/CPSI0/XPSI(NRP,9),IPSIF(99),IPSIL(99),NSPSI(99),NCURP
      PARAMETER (NDIM=1000)
      DIMENSION E(NDIM),XPSIR(NDIM,9),X(NDIM),Y(NDIM)
C
C  ************  Read element x-section tables
C
      DO IEL=1,NELEM(M)
        READ(IRD,1001) IZZ,NSHR,NDATA
 1001   FORMAT(47X,I3,11X,I3,10X,I4)
        IF(INFO.GE.2) WRITE(IWR,2001) IZZ,NSHR,NDATA
 2001   FORMAT(/1X,'***  Positron impact ionization cross sections, ',
     1    ' IZ =',I3,',  NSHELL =',I3,',  NDATA =',I4)
        IF(IZZ.NE.IZ(M,IEL)) STOP 'PSIR. Corrupt material data file.'
        IF(NDATA.GT.NDIM) STOP 'PSIR. Too many data points.'
        IF(NSHR.GT.9) STOP 'PSIR. Too many shells.'
        DO IE=1,NDATA
          READ(IRD,*) E(IE),(XPSIR(IE,IS),IS=1,NSHR)
        ENDDO
C
C  ****  Remove shells with ionization energies less than 100 eV.
C
        IF(NSHR.GT.1) THEN
          NSHA=NSHR
          DO IS=NSHA,1,-1
            IF(EB(IZZ,IS).LT.100.0D0) THEN
              NSHR=NSHR-1
            ELSE
              GO TO 1
            ENDIF
          ENDDO
          IF(NSHR.LT.1) NSHR=1
        ENDIF
    1   CONTINUE
C
        IF(INFO.GE.2) THEN
          IF(NSHR.EQ.1) THEN
            WRITE(IWR,2102)
          ELSE IF(NSHR.EQ.2) THEN
            WRITE(IWR,2202)
          ELSE IF(NSHR.EQ.3) THEN
            WRITE(IWR,2302)
          ELSE IF(NSHR.EQ.4) THEN
            WRITE(IWR,2402)
          ELSE IF(NSHR.EQ.5) THEN
            WRITE(IWR,2502)
          ELSE IF(NSHR.EQ.6) THEN
            WRITE(IWR,2602)
          ELSE IF(NSHR.EQ.7) THEN
            WRITE(IWR,2702)
          ELSE IF(NSHR.EQ.8) THEN
            WRITE(IWR,2802)
          ELSE
            WRITE(IWR,2902)
          ENDIF
          DO IE=1,NDATA
            WRITE(IWR,'(1P,12E12.5)') E(IE),(XPSIR(IE,IS),IS=1,NSHR)
          ENDDO
        ENDIF
 2102   FORMAT(/3X,'Energy',7X,'CS-K',
     1    /4X,'(eV)',2X,5X,'(cm**2)',/1X,24('-'))
 2202   FORMAT(/3X,'Energy',7X,'CS-K',8X,'CS-L1',
     1    /4X,'(eV)',2X,2(5X,'(cm**2)'),/1X,36('-'))
 2302   FORMAT(/3X,'Energy',7X,'CS-K',8X,'CS-L1',7X,
     1    'CS-L2',/4X,'(eV)',2X,3(5X,'(cm**2)'),/1X,48('-'))
 2402   FORMAT(/3X,'Energy',7X,'CS-K',8X,'CS-L1',7X,
     1    'CS-L2',7X,'CS-L3',
     2    /4X,'(eV)',2X,4(5X,'(cm**2)'),/1X,60('-'))
 2502   FORMAT(/3X,'Energy',7X,'CS-K',8X,'CS-L1',7X,
     1    'CS-L2',7X,'CS-L3',7X,'CS-M1',
     2    /4X,'(eV)',2X,5(5X,'(cm**2)'),/1X,72('-'))
 2602   FORMAT(/3X,'Energy',7X,'CS-K',8X,'CS-L1',7X,
     1    'CS-L2',7X,'CS-L3',7X,'CS-M1',7X,'CS-M2',
     2    /4X,'(eV)',2X,6(5X,'(cm**2)'),/1X,84('-'))
 2702   FORMAT(/3X,'Energy',7X,'CS-K',8X,'CS-L1',7X,
     1    'CS-L2',7X,'CS-L3',7X,'CS-M1',7X,'CS-M2',7X,'CS-M3',
     2    /4X,'(eV)',2X,7(5X,'(cm**2)'),/1X,96('-'))
 2802   FORMAT(/3X,'Energy',7X,'CS-K',8X,'CS-L1',7X,
     1    'CS-L2',7X,'CS-L3',7X,'CS-M1',7X,'CS-M2',7X,'CS-M3',7X,
     2    'CS-M4',/4X,'(eV)',2X,8(5X,'(cm**2)'),
     3    /1X,108('-'))
 2902   FORMAT(/3X,'Energy',7X,'CS-K',8X,'CS-L1',7X,
     1    'CS-L2',7X,'CS-L3',7X,'CS-M1',7X,'CS-M2',7X,'CS-M3',7X,
     2    'CS-M4',7X,'CS-M5',/4X,'(eV)',2X,9(5X,'(cm**2)'),
     3    /1X,120('-'))
C
        NSPSI(IZZ)=NSHR
        IF(IPSIF(IZZ).EQ.0) THEN
          IPSIF(IZZ)=NCURP+1
          IF(NCURP+NEGP.GT.NRP) THEN
            WRITE(IWR,*) 'Insufficient memory storage in PSIR.'
            WRITE(IWR,*) 'Increase the value of the parameter NRP to',
     1        NCURP+NEGP
            WRITE(26,*) 'Insufficient memory storage in PSIR.'
            WRITE(26,*) 'Increase the value of the parameter NRP to',
     1        NCURP+NEGP
            STOP 'PSIR. Insufficient memory storage.'
          ENDIF
          DO IS=1,NSHR
            N=0
            DO I=1,NDATA
              IF(XPSIR(I,IS).GT.1.0D-35) THEN
                N=N+1
                X(N)=LOG(E(I))
                IF(N.GT.1) X(N)=MAX(X(N),X(N-1)+1.0D-6)
                Y(N)=LOG(XPSIR(I,IS))
              ENDIF
            ENDDO
            IF(N.GT.4) THEN
              DO I=1,NEGP
                IC=NCURP+I
                XC=DLEMP(I)
                IF(XC.GT.X(1)) THEN
                  CALL FINDI(X,XC,N,J)
                  IF(J.EQ.N) J=N-1
                  DX=X(J+1)-X(J)
                  IF(DX.GT.1.0D-6) THEN
                    XPSI(IC,IS)=Y(J)+(XC-X(J))*(Y(J+1)-Y(J))/DX
                  ELSE
                    XPSI(IC,IS)=(Y(J+1)+Y(J))/2.0D0
                  ENDIF
                ELSE
                  XPSI(IC,IS)=-80.6D0
                ENDIF
              ENDDO
            ELSE
              DO I=1,NEGP
                IC=NCURP+I
                XPSI(IC,IS)=-80.6D0
              ENDDO
            ENDIF
          ENDDO
          NCURP=NCURP+NEGP
          IPSIL(IZZ)=NCURP
        ENDIF
      ENDDO
C
      RETURN
      END
C  *********************************************************************
C                       SUBROUTINE PSI0
C  *********************************************************************
      SUBROUTINE PSI0
C
C  This subroutine sets all variables in common /CPSI0/ to zero.
C
C  It has to be invoked before reading the first material definition
C  file.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N)
C  ****  Ionization cross section tables.
      PARAMETER (NRP=6000)
      COMMON/CPSI0/XPSI(NRP,9),IPSIF(99),IPSIL(99),NSPSI(99),NCURP
C
      DO I=1,99
        IPSIF(I)=0
        IPSIL(I)=0
        NSPSI(I)=0
      ENDDO
C
      DO I=1,NRP
        DO J=1,9
          XPSI(I,J)=-80.6D0
        ENDDO
      ENDDO
      NCURP=0
      RETURN
      END
C  *********************************************************************
C                       SUBROUTINE PSIW
C  *********************************************************************
      SUBROUTINE PSIW(M,IWR)
C
C  This subroutine generates tables of cross sections for inner-shell
C  ionization by positron impact for the elements in material M and
C  writes them on the material data file.
C
C  Data are read from the files 'pdpsiZZ.p05'.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N)
C
      CHARACTER*12 FILEN
      CHARACTER*1 LDIG(10),LDIG1,LDIG2
      DATA LDIG/'0','1','2','3','4','5','6','7','8','9'/
C  ****  Composition data.
      PARAMETER (MAXMAT=10)
      COMMON/COMPOS/STF(MAXMAT,30),ZT(MAXMAT),AT(MAXMAT),RHO(MAXMAT),
     1  VMOL(MAXMAT),IZ(MAXMAT,30),NELEM(MAXMAT)
C
      PARAMETER (NES=400)
      DIMENSION E(NES),XPSIR(NES,9)
C
      DO IEL=1,NELEM(M)
        IZZ=IZ(M,IEL)
        NLD=IZZ
        NLD1=NLD-10*(NLD/10)
        NLD2=(NLD-NLD1)/10
        LDIG1=LDIG(NLD1+1)
        LDIG2=LDIG(NLD2+1)
        FILEN='pdpsi'//LDIG2//LDIG1//'.p05'
        OPEN(3,FILE=FILEN)
        READ(3,*) IZZZ,NSHR
        IF(IZZZ.NE.IZZ) STOP 'ESIW. Corrupt data file.'
        IF(NSHR.GT.9) STOP 'ESIW. Too many shells.'
        DO IE=1,NES
          READ(3,*,END=1) E(IE),(XPSIR(IE,IS),IS=1,NSHR)
          NPTAB=IE
          IF(E(IE).GT.0.999D9) GO TO 1
        ENDDO
    1   CONTINUE
        CLOSE(3)
        WRITE(IWR,2001) IZZ,NSHR,NPTAB
 2001 FORMAT(1X,'***  Positron ionization cross sections,  IZ =',I3,
     1  ',  NSHELL =',I3,',  NDATA =',I4)
        DO IE=1,NPTAB
          WRITE(IWR,2002) E(IE),(XPSIR(IE,IS),IS=1,NSHR)
        ENDDO
 2002 FORMAT(1P,10E12.5)
      ENDDO
      RETURN
      END
C  *********************************************************************
C                       SUBROUTINE EBR
C  *********************************************************************
      SUBROUTINE EBR(E,W,M)
C
C  Simulation of bremsstrahlung emission by electrons or positrons in
C  the material M.
C
      use mtpar
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N)
C  ****  Simulation parameters.
      PARAMETER (MAXMAT=10)
      COMMON/CSIMPA/EABS(3,MAXMAT),C1(MAXMAT),C2(MAXMAT),WCC(MAXMAT),
     1  WCR(MAXMAT)
C  ****  Energy grid and interpolation constants for the current energy.
      PARAMETER (NEGP=200)
      COMMON/CEGRID/EL,EU,ET(NEGP),DLEMP(NEGP),DLEMP1,DLFC,
     1  XEL,XE,XEK,KE
C  ****  Bremsstrahlung emission.
      PARAMETER (NBW=32)
      COMMON/CEBR/WB(NBW),PBCUT(MAXMAT,NEGP),WBCUT(MAXMAT,NEGP),
     1  PDFB(MAXMAT,NEGP,NBW),PACB(MAXMAT,NEGP,NBW),ZBR2(MAXMAT)
C
      ! EXTERNALRAND
C
      IF(WCR(M).GT.E) THEN
        W=0.0D0
        RETURN
      ENDIF
C
C  ****  Selection of the energy grid point.
C
        status=rmt(r)
      IF(r.LT.XEK) THEN
        IE=KE+1
      ELSE
        IE=KE
      ENDIF
C  ****  Pointer.
    1 CONTINUE
        status=rmt(r)
      PT=PBCUT(M,IE)+r*(PACB(M,IE,NBW)-PBCUT(M,IE))
C  ****  Binary search of the W-interval.
      I=1
      J=NBW
    2 K=(I+J)/2
      IF(PT.GT.PACB(M,IE,K)) THEN
        I=K
      ELSE
        J=K
      ENDIF
      IF(J-I.GT.1) GO TO 2
C  ****  Sampling the photon energy (rejection method).
      W1=WB(I)
      W2=WB(I+1)
      DW=W2-W1
      DP=PDFB(M,IE,I+1)-PDFB(M,IE,I)
      B=DP/DW
      A=PDFB(M,IE,I)-B*W1
      IF(W1.LT.WBCUT(M,IE)) W1=WBCUT(M,IE)
      IF(W2.LT.W1) THEN
        WRITE(26,*) ' **** WARNING: EBR. Conflicting end-point values.'
        W=W1
        RETURN
      ENDIF
      PMAX=MAX(A+B*W1,A+B*W2)
    3 CONTINUE
        status=rmt(r)
      W=W1*(W2/W1)**r
        status=rmt(r)
      IF(r*PMAX.GT.A+B*W) GO TO 3
      W=W*E
      IF(W.LT.WCR(M)) GO TO 1
      RETURN
      END
C  *********************************************************************
C                       SUBROUTINE EBRR
C  *********************************************************************
      SUBROUTINE EBRR(WCR,M,IRD,IWR,INFO)
C
C  This subroutine reads the bremss scaled cross section for electrons
C  in material M from the material data file. It computes restricted
C  integrated cross sections and initializes the algorithm for simula-
C  tion of bremss emission by electrons and positrons.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N)
C  ****  Physical constants (values adopted by Seltzer and Berger).
      PARAMETER (REV=5.10998902D5)  ! Electron rest energy (eV)
      PARAMETER (ELRAD=2.817940285D-13)  ! Class. electron radius (cm)
      PARAMETER (TREV=2.0D0*REV)
C  ****  Energy grid and interpolation constants for the current energy.
      PARAMETER (NEGP=200)
      COMMON/CEGRID/EL,EU,ET(NEGP),DLEMP(NEGP),DLEMP1,DLFC,
     1  XEL,XE,XEK,KE
      DIMENSION A(NEGP),B(NEGP),C(NEGP),D(NEGP),PAC(NEGP),PDF(NEGP)
C  ****  Bremsstrahlung emission.
      PARAMETER (MAXMAT=10)
      PARAMETER (NBE=57, NBW=32)
      COMMON/CEBR/WB(NBW),PBCUT(MAXMAT,NEGP),WBCUT(MAXMAT,NEGP),
     1  PDFB(MAXMAT,NEGP,NBW),PACB(MAXMAT,NEGP,NBW),ZBR2(MAXMAT)
      COMMON/CEBR01/EBT(NBE),XS(NBE,NBW),TXS(NBE),P0(NEGP,NBW),X(NBE),
     1  Y(NBE),MLAST
      DIMENSION WB0(NBW)
      DATA WB0/1.0D-12,0.025D0,0.05D0,0.075D0,0.1D0,0.15D0,0.2D0,0.25D0,
     1   0.3D0,0.35D0,0.4D0,0.45D0,0.5D0,0.55D0,0.6D0,0.65D0,0.7D0,
     2   0.75D0,0.8D0,0.85D0,0.9D0,0.925D0,0.95D0,0.97D0,0.99D0,
     3   0.995D0,0.999D0,0.9995D0,0.9999D0,0.99995D0,0.99999D0,1.0D0/
C
C  ****  Reading the scaled cross section table.
C
      READ(IRD,5001) ZBR,NBER
 5001 FORMAT(45X,E12.5,10X,I4)
      IF(INFO.GE.2) WRITE(IWR,2001) ZBR,NBER
 2001 FORMAT(/1X,'*** Electron scaled bremss x-section,  ZEQ =',
     1  1P,E12.5,',  NDATA =',I4)
      IF(NBER.NE.NBE) STOP 'EBRR. Inconsistent format.'
      ZBR2(M)=ZBR*ZBR
C
      DO IE=1,NBE
        READ(IRD,5002) EBT(IE),(XS(IE,IW),IW=1,NBW),TXS(IE)
        IF(INFO.GE.2) WRITE(IWR,2002)
     1    EBT(IE),(XS(IE,IW),IW=1,NBW),TXS(IE)
        X(IE)=LOG(EBT(IE))
      ENDDO
 5002 FORMAT(E9.2,5E12.5,/9X,5E12.5,/9X,5E12.5,/9X,5E12.5,
     1  /9X,5E12.5,/9X,5E12.5,/9X,2E12.5,36X,E10.3)
 2002 FORMAT(1P,E9.2,5E12.5,/9X,5E12.5,/9X,5E12.5,/9X,5E12.5,
     1  /9X,5E12.5,/9X,5E12.5,/9X,2E12.5,36X,E10.3)
      MLAST=M
C
C  ************  Initialization of the calculation routines.
C
      DO I=1,NBW
        WB(I)=WB0(I)
      ENDDO
C
C  ****  Compute the scaled energy loss distribution and sampling
C        parameters for the energies in the simulation grid.
C
C  ****  Interpolation in E.
C
      DO IW=1,NBW
        DO IE=1,NBE
          Y(IE)=LOG(XS(IE,IW))
        ENDDO
        CALL SPLINE(X,Y,A,B,C,D,0.0D0,0.0D0,NBE)
        DO I=1,NEGP
          ELL=DLEMP(I)
          CALL FINDI(X,ELL,NBE,J)
          IF(ELL.GT.X(1)) THEN
            P0(I,IW)=EXP(A(J)+ELL*(B(J)+ELL*(C(J)+ELL*D(J))))
          ELSE
            F1=A(1)+X(1)*(B(1)+X(1)*(C(1)+X(1)*D(1)))
            FP1=B(1)+X(1)*(2.0D0*C(1)+X(1)*3.0D0*D(1))
            P0(I,IW)=EXP(F1+FP1*(ELL-X(1)))
          ENDIF
        ENDDO
      ENDDO
C
      DO IE=1,NEGP
        DO IW=1,NBW
          PDF(IW)=P0(IE,IW)
        ENDDO
C
        CALL RLPAC(WB,PDF,PAC,NBW)
        DO IW=1,NBW
          PDFB(M,IE,IW)=PDF(IW)
          PACB(M,IE,IW)=PAC(IW)
        ENDDO
C  ****  The cutoff scaled energy loss is slightly modified to ensure
C that the sampling routine EBR covers the allowed energy loss interval.
        IF(IE.LT.NEGP) THEN
          XC=WCR/ET(IE+1)
        ELSE
          XC=WCR/ET(NEGP)
        ENDIF
        IF(XC.LT.1.0D0) THEN
          PBCUT(M,IE)=RLMOM(WB,PDF,XC,NBW,-1)
          WBCUT(M,IE)=XC
        ELSE
          PBCUT(M,IE)=RLMOM(WB,PDF,1.0D0,NBW,-1)
          WBCUT(M,IE)=1.0D0
        ENDIF
      ENDDO
C
      RETURN
      END
C  *********************************************************************
C                       SUBROUTINE EBRW
C  *********************************************************************
      SUBROUTINE EBRW(WCR,M,IWR)
C
C  This subroutine generates a table of the scaled energy-loss cross
C  section for bremsstrahlung emission by electrons in material M. Data
C  are read from the files 'pdebrZZ.p05'.
C
C  The cutoff energy loss WCR is entered as an argument for consistency
C  with subroutine EBRR.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N)
C  ****  Physical constants (values adopted by Seltzer and Berger).
      PARAMETER (REV=5.10998902D5)  ! Electron rest energy (eV)
      PARAMETER (ELRAD=2.817940285D-13)  ! Class. electron radius (cm)
      PARAMETER (TREV=2.0D0*REV)
C
      CHARACTER*12 FILEN
      CHARACTER*1 LDIG(10),LDIG1,LDIG2
      DATA LDIG/'0','1','2','3','4','5','6','7','8','9'/
C  ****  Composition data.
      PARAMETER (MAXMAT=10)
      COMMON/COMPOS/STF(MAXMAT,30),ZT(MAXMAT),AT(MAXMAT),RHO(MAXMAT),
     1  VMOL(MAXMAT),IZ(MAXMAT,30),NELEM(MAXMAT)
C  ****  Energy grid and interpolation constants for the current energy.
      PARAMETER (NEGP=200)
      COMMON/CEGRID/EL,EU,ET(NEGP),DLEMP(NEGP),DLEMP1,DLFC,
     1  XEL,XE,XEK,KE
      DIMENSION A(NEGP),B(NEGP),C(NEGP),D(NEGP)
C  ****  Bremsstrahlung emission.
      PARAMETER (NBE=57, NBW=32)
      COMMON/CEBR/WB(NBW),PBCUT(MAXMAT,NEGP),WBCUT(MAXMAT,NEGP),
     1  PDFB(MAXMAT,NEGP,NBW),PACB(MAXMAT,NEGP,NBW),ZBR2(MAXMAT)
      COMMON/CEBR01/EBT(NBE),XS(NBE,NBW),TXS(NBE),P0(NEGP,NBW),X(NBE),
     1  Y(NBE),MLAST
      DIMENSION WB0(NBW),PDF(NBE)
      DATA WB0/1.0D-12,0.025D0,0.05D0,0.075D0,0.1D0,0.15D0,0.2D0,0.25D0,
     1   0.3D0,0.35D0,0.4D0,0.45D0,0.5D0,0.55D0,0.6D0,0.65D0,0.7D0,
     2   0.75D0,0.8D0,0.85D0,0.9D0,0.925D0,0.95D0,0.97D0,0.99D0,
     3   0.995D0,0.999D0,0.9995D0,0.9999D0,0.99995D0,0.99999D0,1.0D0/
C
C  ****  'Equivalent' atomic number.
C
      SUMZ2=0.0D0
      SUMS=0.0D0
      DO IEL=1,NELEM(M)
        SUMZ2=SUMZ2+STF(M,IEL)*IZ(M,IEL)**2
        SUMS=SUMS+STF(M,IEL)
      ENDDO
      ZBR2(M)=SUMZ2/SUMS
C
C  ****  Building the scaled cross section table.
C
      DO IE=1,NBE
        TXS(IE)=0.0D0
        DO IW=1,NBW
          XS(IE,IW)=0.0D0
        ENDDO
      ENDDO
C
      DO IEL=1,NELEM(M)
        IZZ=IZ(M,IEL)
        WGHT=STF(M,IEL)*IZZ*IZZ/ZBR2(M)
        NLD=IZZ
        NLD1=NLD-10*(NLD/10)
        NLD2=(NLD-NLD1)/10
        LDIG1=LDIG(NLD1+1)
        LDIG2=LDIG(NLD2+1)
        FILEN='pdebr'//LDIG2//LDIG1//'.p05'
        OPEN(3,FILE=FILEN)
        READ(3,*) IZZZ
        IF(IZZZ.NE.IZZ) STOP 'EBRW. Corrupt file.'
        DO IE=1,NBE
          READ(3,1001) EBT(IE),(PDF(IW),IW=1,NBW),TXSP
 1001     FORMAT(E9.2,5E12.5,/9X,5E12.5,/9X,5E12.5,/9X,5E12.5,
     1      /9X,5E12.5,/9X,5E12.5,/9X,2E12.5,36X,E10.3)
          TXS(IE)=TXS(IE)+WGHT*TXSP
          DO IW=1,NBW
            XS(IE,IW)=XS(IE,IW)+WGHT*PDF(IW)
          ENDDO
        ENDDO
        CLOSE(3)
      ENDDO
C
C  ****  The energy loss spectrum is re-normalized to reproduce the
C        total scaled cross section of Berger and Seltzer.
C
      DO IE=1,NBE
        DO IW=1,NBW
          X(IW)=WB0(IW)
          Y(IW)=XS(IE,IW)
        ENDDO
        RSUM=RLMOM(X,Y,1.0D0,NBW,0)
        FACT=(EBT(IE)+REV)*1.0D-27*137.03604D0/(ELRAD**2*(EBT(IE)+TREV))
        FNORM=TXS(IE)/(RSUM*FACT)
        TST=100.0D0*ABS(FNORM-1.0D0)
        IF(TST.GT.1.0D0) STOP 'EBRW. Check the bremss database file.'
        DO IW=1,NBW
          XS(IE,IW)=XS(IE,IW)*FNORM
        ENDDO
      ENDDO
C
C  ****  Write output scaled x-section table.
C
      WRITE(IWR,2001) SQRT(ZBR2(M)),NBE
 2001 FORMAT(1X,'*** Electron scaled bremss x-section,  ZEQ =',1P,E12.5,
     1  ',  NDATA =',I4)
      DO IE=1,NBE
        WRITE(IWR,2002) EBT(IE),(XS(IE,IW),IW=1,NBW),TXS(IE)
      ENDDO
 2002 FORMAT(1P,E9.2,5E12.5,/9X,5E12.5,/9X,5E12.5,/9X,5E12.5,
     1  /9X,5E12.5,/9X,5E12.5,/9X,2E12.5,36X,E10.3)
      MLAST=M
C
C  ************  Initialization of the calculation routines.
C
      DO I=1,NBW
        WB(I)=WB0(I)
      ENDDO
C
C  ****  Compute the scaled energy loss distribution and sampling
C        parameters for the energies in the simulation grid.
C
C  ****  Interpolation in E.

      DO IE=1,NBE
        X(IE)=LOG(EBT(IE))
      ENDDO
      DO IW=1,NBW
        DO IE=1,NBE
          Y(IE)=LOG(XS(IE,IW))
        ENDDO
        CALL SPLINE(X,Y,A,B,C,D,0.0D0,0.0D0,NBE)
        DO I=1,NEGP
          ELL=DLEMP(I)
          IF(ELL.GT.X(1)) THEN
            CALL FINDI(X,ELL,NBE,J)
            P0(I,IW)=EXP(A(J)+ELL*(B(J)+ELL*(C(J)+ELL*D(J))))
          ELSE
            F1=A(1)+X(1)*(B(1)+X(1)*(C(1)+X(1)*D(1)))
            FP1=B(1)+X(1)*(2.0D0*C(1)+X(1)*3.0D0*D(1))
            P0(I,IW)=EXP(F1+FP1*(ELL-X(1)))
          ENDIF
        ENDDO
      ENDDO
C
      RETURN
      END
C  *********************************************************************
C                       SUBROUTINE EBRTX
C  *********************************************************************
      SUBROUTINE EBRTX(E,WCR,XH0,XH1,XH2,XS1,XS2)
C
C  Integrated cross sections for bremss emission by electrons of energy
C  E, restricted to energy losses larger than and less than the cutoff
C  energy WCR. This subroutine uses data that are overwritten by EBRR,
C  so that the active material is the one loaded last.
C
C  Output arguments:
C    XH0 ... total cross section for hard emission (cm**2).
C    XH1 ... stopping cross section for hard emission (eV*cm**2)
C    XH2 ... straggling cross section for hard emission (eV**2*cm**2).
C    XS1 ... stopping cross section for soft emission (eV*cm**2)
C    XS2 ... straggling cross section for soft emission (eV**2*cm**2).
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N)
      PARAMETER (REV=5.10998902D5)  ! Electron rest energy (eV)
      PARAMETER (TREV=2.0D0*REV)
C  ****  Energy grid and interpolation constants for the current energy.
      PARAMETER (NEGP=200)
      COMMON/CEGRID/EL,EU,ET(NEGP),DLEMP(NEGP),DLEMP1,DLFC,
     1  XEL,XE,XEK,KE
C  ****  Bremsstrahlung emission.
      PARAMETER (MAXMAT=10)
      PARAMETER (NBE=57, NBW=32)
      COMMON/CEBR/WB(NBW),PBCUT(MAXMAT,NEGP),WBCUT(MAXMAT,NEGP),
     1  PDFB(MAXMAT,NEGP,NBW),PACB(MAXMAT,NEGP,NBW),ZBR2(MAXMAT)
      COMMON/CEBR01/EBT(NBE),XS(NBE,NBW),TXS(NBE),P0(NEGP,NBW),X(NBE),
     1  Y(NBE),MLAST
C
      XEL=MAX(LOG(E),DLEMP1)
      XE=1.0D0+(XEL-DLEMP1)*DLFC
      KE=XE
      XEK=XE-KE
C  ****  Global x-section factor.
      FACT=ZBR2(MLAST)*((E+REV)**2/(E*(E+TREV)))*1.0D-27
C
C  ****  Moments of the scaled bremss x-section.
C
      WCRE=WCR/E
      DO IW=1,NBW
        X(IW)=WB(IW)
        Y(IW)=P0(KE,IW)
      ENDDO
      XH0A=RLMOM(X,Y,X(NBW),NBW,-1)-RLMOM(X,Y,WCRE,NBW,-1)
      XS1A=RLMOM(X,Y,WCRE,NBW,0)
      XS2A=RLMOM(X,Y,WCRE,NBW,1)
      XH1A=RLMOM(X,Y,X(NBW),NBW,0)-XS1A
      XH2A=RLMOM(X,Y,X(NBW),NBW,1)-XS2A
      DO IW=1,NBW
        Y(IW)=P0(MIN(KE+1,NEGP),IW)
      ENDDO
      XH0B=RLMOM(X,Y,X(NBW),NBW,-1)-RLMOM(X,Y,WCRE,NBW,-1)
      XS1B=RLMOM(X,Y,WCRE,NBW,0)
      XS2B=RLMOM(X,Y,WCRE,NBW,1)
      XH1B=RLMOM(X,Y,X(NBW),NBW,0)-XS1B
      XH2B=RLMOM(X,Y,X(NBW),NBW,1)-XS2B
C
      XH0=((1.0D0-XEK)*XH0A+XEK*XH0B)*FACT
      XS1=((1.0D0-XEK)*XS1A+XEK*XS1B)*FACT*E
      XH1=((1.0D0-XEK)*XH1A+XEK*XH1B)*FACT*E
      XS2=((1.0D0-XEK)*XS2A+XEK*XS2B)*FACT*E*E
      XH2=((1.0D0-XEK)*XH2A+XEK*XH2B)*FACT*E*E
      RETURN
      END
C  *********************************************************************
C                       SUBROUTINE PBRTX
C  *********************************************************************
      SUBROUTINE PBRTX(E,WCR,XH0,XH1,XH2,XS1,XS2)
C
C  Integrated cross sections for bremss emission by positrons of energy
C  E, restricted to energy losses larger than and less than the cutoff
C  energy WCR. This subroutine uses data that are overwritten by EBRR,
C  so that the active material is the one loaded last.
C
C  Output arguments:
C    XH0 ... total cross section for hard emission (cm**2).
C    XH1 ... stopping cross section for hard emission (eV*cm**2)
C    XH2 ... straggling cross section for hard emission (eV**2*cm**2).
C    XS1 ... stopping cross section for soft emission (eV*cm**2)
C    XS2 ... straggling cross section for soft emission (eV**2*cm**2).
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N)
      PARAMETER (REV=5.10998902D5)  ! Electron rest energy (eV)
      PARAMETER (TREV=2.0D0*REV)
C  ****  Energy grid and interpolation constants for the current energy.
      PARAMETER (NEGP=200)
      COMMON/CEGRID/EL,EU,ET(NEGP),DLEMP(NEGP),DLEMP1,DLFC,
     1  XEL,XE,XEK,KE
C  ****  Bremsstrahlung emission.
      PARAMETER (MAXMAT=10)
      PARAMETER (NBE=57, NBW=32)
      COMMON/CEBR/WB(NBW),PBCUT(MAXMAT,NEGP),WBCUT(MAXMAT,NEGP),
     1  PDFB(MAXMAT,NEGP,NBW),PACB(MAXMAT,NEGP,NBW),ZBR2(MAXMAT)
      COMMON/CEBR01/EBT(NBE),XS(NBE,NBW),TXS(NBE),P0(NEGP,NBW),X(NBE),
     1  Y(NBE),MLAST
C
      XEL=MAX(LOG(E),DLEMP1)
      XE=1.0D0+(XEL-DLEMP1)*DLFC
      KE=XE
      XEK=XE-KE
C  ****  Global x-section factor.
      FACT=ZBR2(MLAST)*((E+REV)**2/(E*(E+TREV)))*1.0D-27
C  ****  Positron correction factor.
      T=LOG(1.0D0+1.0D6*E/(REV*ZBR2(MLAST)))
      FPOS=1.0D0-EXP(-T*(1.2359D-1-T*(6.1274D-2-T*(3.1516D-2-T
     1    *(7.7446D-3-T*(1.0595D-3-T*(7.0568D-5-T*1.8080D-6)))))))
      FACT=FACT*FPOS
C
C  ****  Moments of the scaled bremss x-section.
C
      WCRE=WCR/E
      DO IW=1,NBW
        X(IW)=WB(IW)
        Y(IW)=P0(KE,IW)
      ENDDO
      XH0A=RLMOM(X,Y,X(NBW),NBW,-1)-RLMOM(X,Y,WCRE,NBW,-1)
      XS1A=RLMOM(X,Y,WCRE,NBW,0)
      XS2A=RLMOM(X,Y,WCRE,NBW,1)
      XH1A=RLMOM(X,Y,X(NBW),NBW,0)-XS1A
      XH2A=RLMOM(X,Y,X(NBW),NBW,1)-XS2A
      DO IW=1,NBW
        Y(IW)=P0(MIN(KE+1,NEGP),IW)
      ENDDO
      XH0B=RLMOM(X,Y,X(NBW),NBW,-1)-RLMOM(X,Y,WCRE,NBW,-1)
      XS1B=RLMOM(X,Y,WCRE,NBW,0)
      XS2B=RLMOM(X,Y,WCRE,NBW,1)
      XH1B=RLMOM(X,Y,X(NBW),NBW,0)-XS1B
      XH2B=RLMOM(X,Y,X(NBW),NBW,1)-XS2B
C
      XH0=((1.0D0-XEK)*XH0A+XEK*XH0B)*FACT
      XS1=((1.0D0-XEK)*XS1A+XEK*XS1B)*FACT*E
      XH1=((1.0D0-XEK)*XH1A+XEK*XH1B)*FACT*E
      XS2=((1.0D0-XEK)*XS2A+XEK*XS2B)*FACT*E*E
      XH2=((1.0D0-XEK)*XH2A+XEK*XH2B)*FACT*E*E
      RETURN
      END
C  *********************************************************************
C                       FUNCTION RLMOM
C  *********************************************************************
      FUNCTION RLMOM(X,FCT,XC,NP,MOM)
C
C  Calculation of the integral of (X**MOM)*FCT(X) over the interval from
C  X(1) to XC, obtained by linear interpolation on a table of FCT.
C  The independent variable X is assumed to take only positive values.
C
C    X ....... array of values of the variable (in increasing order).
C    FCT ..... corresponding FCT values.
C    NP ...... number of points in the table.
C    XC ...... upper limit of the integral, X(1).LE.XC.LE.X(NP).
C    MOM ..... moment order (GE.-1).
C    RLMOM ... integral of (X**MOM)*FCT(X) over the interval from X(1)
C              to XC.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N)
      PARAMETER (EPS=1.0D-35)
      DIMENSION X(NP),FCT(NP)
C
      IF(MOM.LT.-1) STOP 'RLMOM. Error code 0.'
      IF(NP.LT.2) STOP 'RLMOM. Error code 1.'
      IF(X(1).LT.0.0D0) STOP 'RLMOM. Error code 2'
      DO I=2,NP
        IF(X(I).LT.0.0D0)
     1    STOP 'RLMOM. Error code 3'
        IF(X(I).LT.X(I-1)) STOP 'RLMOM. Error code 4.'
      ENDDO
C
      RLMOM=0.0D0
      IF(XC.LT.X(1)) RETURN
      IEND=0
      XT=MIN(XC,X(NP))
      DO I=1,NP-1
        X1=MAX(X(I),EPS)
        Y1=FCT(I)
        X2=MAX(X(I+1),EPS)
        Y2=FCT(I+1)
        IF(XT.LT.X2) THEN
          XTC=XT
          IEND=1
        ELSE
          XTC=X2
        ENDIF
        DX=X2-X1
        DY=Y2-Y1
        IF(ABS(DX).GT.1.0D-14*ABS(DY)) THEN
          B=DY/DX
          A=Y1-B*X1
          IF(MOM.EQ.-1) THEN
            DS=A*LOG(XTC/X1)+B*(XTC-X1)
          ELSE
            DS=A*(XTC**(MOM+1)-X1**(MOM+1))/DBLE(MOM+1)
     1        +B*(XTC**(MOM+2)-X1**(MOM+2))/DBLE(MOM+2)
          ENDIF
        ELSE
          DS=0.5D0*(Y1+Y2)*(XTC-X1)*XTC**MOM
        ENDIF
        RLMOM=RLMOM+DS
        IF(IEND.NE.0) RETURN
      ENDDO
      RETURN
      END
C  *********************************************************************
C                       SUBROUTINE RLPAC
C  *********************************************************************
      SUBROUTINE RLPAC(X,PDF,PAC,NP)
C
C  Cumulative distribution function of PDF(X)/X, obtained from linear
C  interpolation on a table of PDF.
C  The independent variable X is assumed to take only positive values.
C
C    X ..... array of values of the variable (in increasing order).
C    PDF ... corresponding PDF values.
C    PAC ... cumulative probability function.
C    NP .... number of points in the table.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N)
      PARAMETER (EPS=1.0D-35)
      DIMENSION X(NP),PDF(NP),PAC(NP)
C
      PAC(1)=0.0D0
      DO I=2,NP
        X1=MAX(X(I-1),EPS)
        Y1=PDF(I-1)
        X2=MAX(X(I),EPS)
        Y2=PDF(I)
        DX=X2-X1
        DY=Y2-Y1
        B=DY/DX
        A=Y1-B*X1
        DS=A*LOG(X2/X1)+B*(X2-X1)
        PAC(I)=PAC(I-1)+DS
      ENDDO
      RETURN
      END
C  *********************************************************************
C                       SUBROUTINE EBRA
C  *********************************************************************
      SUBROUTINE EBRA(E,DE,CDT,M)
C
C  Random sampling of the initial direction of bremss photons, relative
C  to the direction of the projectile.
C  Numerical fit/interpolation of partial-wave shape functions given by
C  Kissel, Quarles and Pratt; ANDT 28(1993)381.
C
C  Input parameters:
C    M ..... material where the projectile moves.
C    E ..... kinetic energy of the projectile.
C    DE .... energy of the emitted photon.
C  Output parameter:
C    CDT ... cosine of the polar emission angle.
C
      use mtpar
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N)
      PARAMETER (REV=5.10998902D5)  ! Electron rest energy (eV)
      PARAMETER (TREV=2.0D0*REV)
C  ****  Bremsstrahlung angular distributions.
      PARAMETER (MAXMAT=10)
      COMMON/CBRANG/BET(6),BK(21),BP1(MAXMAT,6,21,4),BP2(MAXMAT,6,21,4),
     1              ZBEQ(MAXMAT)
C
      ! EXTERNALRAND
C
C  ****  Distribution parameters.
C
      BETA=SQRT(E*(E+TREV))/(E+REV)
C
C  A pure dipole distribution is used for E>500 keV.
      IF(E.GT.500.0D3) THEN
        status=rmt(r)
        CDT=2.0D0*r-1.0D0
        status=rmt(r)
        IF(r.GT.0.75D0) THEN
          IF(CDT.LT.0.0D0) THEN
            CDT=-(-CDT)**0.333333333333333D0
          ELSE
            CDT=CDT**0.333333333333333D0
          ENDIF
        ENDIF
        CDT=(CDT+BETA)/(1.0D0+BETA*CDT)
        RETURN
      ENDIF
C
      IF(BETA.GT.BET(6)) THEN
        IE=6
        GO TO 20
      ENDIF
      IF(BETA.LT.BET(1)) THEN
        IE=1
        GO TO 20
      ENDIF
      IE=1
      IE1=6
   10 IET=(IE+IE1)/2
      IF(BETA.GT.BET(IET)) THEN
        IE=IET
      ELSE
        IE1=IET
      ENDIF
      IF(IE1-IE.GT.1) GO TO 10
   20 CONTINUE
C
      RK=1.0D0+20.0D0*DE/E
      IK=MIN(INT(RK),20)
C
      P10=BP1(M,IE,IK,1)+BETA*(BP1(M,IE,IK,2)
     1   +BETA*(BP1(M,IE,IK,3)+BETA*BP1(M,IE,IK,4)))
      P11=BP1(M,IE,IK+1,1)+BETA*(BP1(M,IE,IK+1,2)
     1   +BETA*(BP1(M,IE,IK+1,3)+BETA*BP1(M,IE,IK+1,4)))
      P1=P10+(RK-IK)*(P11-P10)
C
      P20=BP2(M,IE,IK,1)+BETA*(BP2(M,IE,IK,2)
     1   +BETA*(BP2(M,IE,IK,3)+BETA*BP2(M,IE,IK,4)))
      P21=BP2(M,IE,IK+1,1)+BETA*(BP2(M,IE,IK+1,2)
     1   +BETA*(BP2(M,IE,IK+1,3)+BETA*BP2(M,IE,IK+1,4)))
      P2=P20+(RK-IK)*(P21-P20)
C
C  ****  Sampling from the Lorentz-transformed dipole distributions.
C
      P1=MIN(EXP(P1)/BETA,1.0D0)
      BETAP=MIN(MAX(BETA*(1.0D0+P2/BETA),0.0D0),0.999999999D0)
C
        status=rmt(r)
      IF(r.LT.P1) THEN
        status=rmt(r)
    1   CDT=2.0D0*r-1.0D0
        status=rmt(r)
        IF(2.0D0*r.GT.1.0D0+CDT*CDT) GO TO 1
      ELSE
        status=rmt(r)
    2   CDT=2.0D0*r-1.0D0
        status=rmt(r)
        IF(r.GT.1.0D0-CDT*CDT) GO TO 2
      ENDIF
      CDT=(CDT+BETAP)/(1.0D0+BETAP*CDT)
C
      RETURN
      END
C  *********************************************************************
C                       SUBROUTINE BRAR
C  *********************************************************************
      SUBROUTINE BRAR(M,IRD,IWR,INFO)
C
C  This subroutine reads bremsstrahlung angular distribution parameters
C  of material M from the material data file. It also initializes the
C  algorithm for generation of the initial direction of bremsstrahlung
C  photons.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N)
      PARAMETER (REV=5.10998902D5)  ! Electron rest energy (eV)
      PARAMETER (TREV=2.0D0*REV)
C  ****  Bremsstrahlung emission.
      PARAMETER (MAXMAT=10, NEGP=200, NBW=32)
      COMMON/CEBR/WB(NBW),PBCUT(MAXMAT,NEGP),WBCUT(MAXMAT,NEGP),
     1  PDFB(MAXMAT,NEGP,NBW),PACB(MAXMAT,NEGP,NBW),ZBR2(MAXMAT)
C  ****  Bremsstrahlung angular distributions.
      COMMON/CBRANG/BET(6),BK(21),BP1(MAXMAT,6,21,4),BP2(MAXMAT,6,21,4),
     1              ZBEQ(MAXMAT)
C
      DIMENSION E(6),RK(4),Q1R(6,4),Q2R(6,4),Q1(6,21),Q2(6,21)
      DIMENSION X(6),A(6),B(6),C(6),D(6)
C
      E(1)=1.0D3
      E(2)=5.0D3
      E(3)=1.0D4
      E(4)=5.0D4
      E(5)=1.0D5
      E(6)=5.0D5
C
      RK(1)=0.0D0
      RK(2)=0.6D0
      RK(3)=0.8D0
      RK(4)=0.95D0
C
      DO IE=1,6
        BET(IE)=SQRT(E(IE)*(E(IE)+TREV))/(E(IE)+REV)
      ENDDO
C
C  ****  Grid of reduced photon energies.
C
      DO IK=1,21
        BK(IK)=(IK-1)*0.05D0
      ENDDO
C
C  ****  Read angular distribution parameters from file (unit IRD).
C
      READ(IRD,5001) ZEQ,NDATA
 5001 FORMAT(40X,E12.5,10X,I4)
      IF(INFO.GE.2) WRITE(IWR,2001) ZEQ,NDATA
 2001 FORMAT(/1X,'*** Bremss angular distribution,  ZEQ =',
     1  1P,E12.5,',  NDATA =',I4)
      IF(NDATA.NE.24) STOP 'BRAR. Inconsistent data.'
      ZBEQ(M)=MIN(MAX(ZEQ,2.0D0),92.0D0)
C
      DO IE1=1,6
        DO IK1=1,4
          READ(IRD,*) IE,IK,ER,RKR,Q1RR,Q2RR
          IF((ABS(ER-E(IE)).LT.1.0D-6).AND.
     1       (ABS(RKR-RK(IK)).LT.1.0D-6)) THEN
            Q1R(IE,IK)=Q1RR/ZEQ
            Q2R(IE,IK)=Q2RR
          ELSE
            WRITE(26,*) 'Corrupt data file (pdbrang.p05).'
            STOP 'BRAR. Corrupt data file (pdbrang.p05).'
          ENDIF
        ENDDO
      ENDDO
C
      IF(INFO.GE.2) THEN
        DO IE=1,6
          DO IK=1,4
            WRITE(IWR,2002) E(IE),RK(IK),Q1R(IE,IK)*ZEQ,Q2R(IE,IK)
          ENDDO
        ENDDO
 2002   FORMAT(1P,E10.3,E11.3,2E15.7)
      ENDIF
C
C  ****  Expanded table of distribution parameters.
C
      DO IE=1,6
        DO IK=1,4
          X(IK)=LOG(Q1R(IE,IK))
        ENDDO
        CALL SPLINE(RK,X,A,B,C,D,0.0D0,0.0D0,4)
        DO IK=1,21
          CALL FINDI(RK,BK(IK),4,J)
          Q1(IE,IK)=A(J)+BK(IK)*(B(J)+BK(IK)*(C(J)+BK(IK)*D(J)))
        ENDDO
        DO IK=1,4
          X(IK)=Q2R(IE,IK)
        ENDDO
        CALL SPLINE(RK,X,A,B,C,D,0.0D0,0.0D0,4)
        DO IK=1,21
          CALL FINDI(RK,BK(IK),4,J)
          Q2(IE,IK)=A(J)+BK(IK)*(B(J)+BK(IK)*(C(J)+BK(IK)*D(J)))
        ENDDO
      ENDDO
C
C  ****  ... and natural cubic spline interpolations.
C
      DO IK=1,21
        DO IE=1,6
          X(IE)=Q1(IE,IK)
        ENDDO
        CALL SPLINE(BET,X,A,B,C,D,0.0D0,0.0D0,6)
        DO IE=1,6
          BP1(M,IE,IK,1)=A(IE)
          BP1(M,IE,IK,2)=B(IE)
          BP1(M,IE,IK,3)=C(IE)
          BP1(M,IE,IK,4)=D(IE)
        ENDDO
        DO IE=1,6
          X(IE)=Q2(IE,IK)
        ENDDO
        CALL SPLINE(BET,X,A,B,C,D,0.0D0,0.0D0,6)
        DO IE=1,6
          BP2(M,IE,IK,1)=A(IE)
          BP2(M,IE,IK,2)=B(IE)
          BP2(M,IE,IK,3)=C(IE)
          BP2(M,IE,IK,4)=D(IE)
        ENDDO
      ENDDO
C
      RETURN
      END
C  *********************************************************************
C                       SUBROUTINE BRAW
C  *********************************************************************
      SUBROUTINE BRAW(ZEQ,IWR)
C
C  This subroutine generates the parameters of the angular distribution
C  of bremsstrahlung photons for the element of atomic number ZEQ. In
C  the case of compounds (and mixtures) ZEQ is the average atomic number
C  of the elements in the molecule. The evaluated parameters are written
C  on the material definition file. Data are read from the database file
C  'pdbrang.p05'.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N)
      PARAMETER (REV=5.10998902D5)  ! Electron rest energy (eV)
      PARAMETER (TREV=2.0D0*REV)
C
      PARAMETER (MAXMAT=10)
C  ****  Bremsstrahlung angular distributions.
      COMMON/CBRANG/BET(6),BK(21),BP1(MAXMAT,6,21,4),BP2(MAXMAT,6,21,4),
     1              ZBEQ(MAXMAT)
C
      DIMENSION Z(6),E(6),RK(4),P1(6,6,4),P2(6,6,4),Q1(6,4),Q2(6,4)
      DIMENSION X(6),Y(6),A(6),B(6),C(6),D(6)
C
      Z(1)=2.0D0
      Z(2)=8.0D0
      Z(3)=13.0D0
      Z(4)=47.0D0
      Z(5)=79.0D0
      Z(6)=92.0D0
C
      E(1)=1.0D3
      E(2)=5.0D3
      E(3)=1.0D4
      E(4)=5.0D4
      E(5)=1.0D5
      E(6)=5.0D5
C
      RK(1)=0.0D0
      RK(2)=0.6D0
      RK(3)=0.8D0
      RK(4)=0.95D0
C
C  ****  Read database file.
C
      OPEN(3,FILE='pdbrang.p05')
      DO IK1=1,4
        DO IZ1=1,6
          DO IE1=1,6
            READ(3,*) IZ,IE,IK,ZR,ER,RKR,P1R,P2R
             IF((ABS(ZR-Z(IZ)).LT.1.0D-6).AND.
     1          (ABS(ER-E(IE)).LT.1.0D-6).AND.
     2          (ABS(RKR-RK(IK)).LT.1.0D-6)) THEN
               P1(IZ,IE,IK)=P1R
               P2(IZ,IE,IK)=P2R
             ELSE
               WRITE(26,*) 'Corrupt data file (pdbrang.p05).'
               STOP 'BRAW. Corrupt data file (pdbrang.p05).'
             ENDIF
          ENDDO
        ENDDO
      ENDDO
      CLOSE(3)
C
C  ****  Interpolation in Z.
C
      DO IE=1,6
        DO IK=1,4
          DO IZ=1,6
            X(IZ)=LOG(P1(IZ,IE,IK))
            Y(IZ)=P2(IZ,IE,IK)
          ENDDO
          CALL SPLINE(Z,X,A,B,C,D,0.0D0,0.0D0,6)
          CALL FINDI(Z,ZEQ,6,I)
          Q1(IE,IK)=EXP(A(I)+ZEQ*(B(I)+ZEQ*(C(I)+ZEQ*D(I))))
          CALL SPLINE(Z,Y,A,B,C,D,0.0D0,0.0D0,6)
          CALL FINDI(Z,ZEQ,6,I)
          Q2(IE,IK)=A(I)+ZEQ*(B(I)+ZEQ*(C(I)+ZEQ*D(I)))
        ENDDO
      ENDDO
C
C  ****  Write final table of parameters.
C
      NDATA=24
      WRITE(IWR,2001) ZEQ,NDATA
 2001 FORMAT(1X,'*** Bremss angular distribution,  ZEQ =',
     1  1P,E12.5,',  NDATA =',I4)
      DO IE=1,6
        DO IK=1,4
          WRITE(IWR,2002) IE,IK,E(IE),RK(IK),Q1(IE,IK),Q2(IE,IK)
        ENDDO
      ENDDO
 2002 FORMAT(2I2,1P,2E11.3,2E15.7)
C
      RETURN
      END
C  *********************************************************************
C                       SUBROUTINE PANR
C  *********************************************************************
      SUBROUTINE PANR
C
C  Simulation of annihilation of positrons at rest.
C
      use mtpar
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N)
      PARAMETER (REV=5.10998902D5)  ! Electron rest energy (eV)
      PARAMETER (PI=3.1415926535897932D0, TWOPI=PI+PI, TREV=2.0D0*REV)
C  ****  Main-PENELOPE common.
      COMMON/TRACK/E,X,Y,Z,U,V,W,WGHT,KPAR,IBODY,M,ILB(5)
      DIMENSION ILBA(5)
C  ****  Simulation parameters.
      PARAMETER (MAXMAT=10)
      COMMON/CSIMPA/EABS(3,MAXMAT),C1(MAXMAT),C2(MAXMAT),WCC(MAXMAT),
     1  WCR(MAXMAT)
C
      ! EXTERNALRAND
C
      IF(REV.LT.EABS(2,M)) RETURN
C
      US=U
      VS=V
      WS=W
        status=rmt(r)
      CDT1=-1.0D0+2.0D0*r
        status=rmt(r)
      DF=TWOPI*r
      CALL DIRECT(CDT1,DF,US,VS,WS)
      ILBA(1)=ILB(1)+1
      ILBA(2)=3
      ILBA(3)=6
      ILBA(4)=0
      ILBA(5)=ILB(5)
      CALL STORES(REV,X,Y,Z,US,VS,WS,WGHT,2,ILBA)
      CALL STORES(REV,X,Y,Z,-US,-VS,-WS,WGHT,2,ILBA)
C
      RETURN
      END
C  *********************************************************************
C                       SUBROUTINE PAN
C  *********************************************************************
      SUBROUTINE PAN(E1,CDT1,E2,CDT2)
C
C  Simulation of positron annihilation (either at rest or in flight).
C  Ei and CDTi are the energies and polar direction cosines of the two
C  annihilation photons.
C
      use mtpar
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N)
      PARAMETER (REV=5.10998902D5)  ! Electron rest energy (eV)
      PARAMETER (PI=3.1415926535897932D0, TWOPI=PI+PI, TREV=2.0D0*REV)
C  ****  Main-PENELOPE common.
      COMMON/TRACK/E,X,Y,Z,U,V,W,WGHT,KPAR,IBODY,M,ILB(5)
C  ****  Simulation parameters.
      PARAMETER (MAXMAT=10)
      COMMON/CSIMPA/EABS(3,MAXMAT),C1(MAXMAT),C2(MAXMAT),WCC(MAXMAT),
     1  WCR(MAXMAT)
C
      ! EXTERNALRAND
C
C  ****  Slow positrons (assumed at rest).
C
      IF(E.LT.EABS(3,M)) THEN
        E1=0.5D0*(E+TREV)
        E2=E1
        status=rmt(r)
        CDT1=-1.0D0+2.0D0*r
        CDT2=-CDT1
      ELSE
C  ****  Annihilation in flight (two photons with energy and directions
C        determined from the dcs and energy-momentum conservation).
        GAM=1.0D0+MAX(E,1.0D0)/REV
        GAM21=SQRT(GAM*GAM-1.0D0)
        ANI=1.0D0+GAM
        CHIMIN=1.0D0/(ANI+GAM21)
        RCHI=(1.0D0-CHIMIN)/CHIMIN
        GT0=ANI*ANI-2.0D0
    1   CONTINUE
        status=rmt(r)
        CHI=CHIMIN*RCHI**r
        GREJ=ANI*ANI*(1.0D0-CHI)+GAM+GAM-1.0D0/CHI
        status=rmt(r)
        IF(r*GT0.GT.GREJ) GO TO 1
C
        DET=E+TREV
        E1=CHI*DET
        CDT1=(ANI-1.0D0/CHI)/GAM21
        CHIP=1.0D0-CHI
        E2=DET-E1
        CDT2=(ANI-1.0D0/CHIP)/GAM21
      ENDIF
C
      RETURN
      END
C  *********************************************************************
C                       SUBROUTINE PANTX
C  *********************************************************************
      SUBROUTINE PANTX(E,TXS)
C
C  Total cross section (per electron) for annihilation of positrons with
C  kinetic energy E. Computed from Heitler's dcs formula for annihila-
C  tion with free electrons at rest.
C
C  Output argument:
C    XST ... total annihilation cross section (cm**2).
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N)
      PARAMETER (REV=5.10998902D5)  ! Electron rest energy (eV)
      PARAMETER (ELRAD=2.817940285D-13)  ! Class. electron radius (cm)
      PARAMETER (PI=3.1415926535897932D0, PIELR2=PI*ELRAD*ELRAD)
C
      GAM=1.0D0+MAX(E,1.0D0)/REV
      GAM2=GAM*GAM
      F2=GAM2-1.0D0
      F1=SQRT(F2)
      TXS=PIELR2*((GAM2+4.0D0*GAM+1.0D0)*LOG(GAM+F1)/F2
     1   -(GAM+3.0D0)/F1)/(GAM+1.0D0)
      RETURN
      END
C  *********************************************************************
C                       SUBROUTINE GRA
C  *********************************************************************
      SUBROUTINE GRA(E,CDT,M)
C
C  Random sampling of coherent (Rayleigh) scattering.
C
      use mtpar
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N)
      CHARACTER*2 LASYMB
      PARAMETER (REV=5.10998902D5)  ! Electron rest energy (eV)
C  ****  Composition data.
      PARAMETER (MAXMAT=10)
      COMMON/COMPOS/STF(MAXMAT,30),ZT(MAXMAT),AT(MAXMAT),RHO(MAXMAT),
     1  VMOL(MAXMAT),IZ(MAXMAT,30),NELEM(MAXMAT)
C  ****  Element data.
      COMMON/CADATA/ATW(99),EPX(99),RA1(99),RA2(99),RA3(99),RA4(99),
     1  RA5(99),RSCR(99),ETA(99),EB(99,30),IFI(99,30),IKS(99,30),
     2  NSHT(99),LASYMB(99)
C
      COMMON/CGRA/X2COH(MAXMAT,241),PDCOH(MAXMAT,241),FLCOH
C
      ! EXTERNALRAND
C
      X2MAX=2.0D0*LOG(41.2148D0*E/REV)
      IF(X2MAX.LT.X2COH(M,2)) THEN
        JM=1
      ELSE IF(X2MAX.GE.X2COH(M,240)) THEN
        JM=240
      ELSE
        JM=1+(X2MAX-X2COH(M,1))/FLCOH
      ENDIF
      RUMAX=PDCOH(M,JM)+(PDCOH(M,JM+1)-PDCOH(M,JM))
     1     *(X2MAX-X2COH(M,JM))/(X2COH(M,JM+1)-X2COH(M,JM))
C
    1 CONTINUE
        status=rmt(r)
      RU=RUMAX+LOG(r)
C  ****  Binary search.
      J=1
      JU=JM+1
    2 JT=(J+JU)/2
      IF(RU.GT.PDCOH(M,JT)) THEN
        J=JT
      ELSE
        JU=JT
      ENDIF
      IF(JU-J.GT.1) GO TO 2
C
      DENOM=PDCOH(M,J+1)-PDCOH(M,J)
      IF(DENOM.GT.1.0D-12) THEN
        X2RAT=X2COH(M,J)+((X2COH(M,J+1)-X2COH(M,J))
     1       *(RU-PDCOH(M,J))/DENOM)-X2MAX
      ELSE
        X2RAT=X2COH(M,J)-X2MAX
      ENDIF
      CDT=1.0D0-2.0D0*EXP(X2RAT)
C  ****  Rejection.
      G=0.5D0*(1.0D0+CDT*CDT)
        status=rmt(r)
      IF(r.GT.G) GO TO 1
      RETURN
      END
C  *********************************************************************
C                       SUBROUTINE GRAI
C  *********************************************************************
      SUBROUTINE GRAI(M)
C
C  Initialization of random sampling for Rayleigh (coherent) scattering
C  of photons.
C
      use mtpar
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N)
      CHARACTER*2 LASYMB
C  ****  Composition data.
      PARAMETER (MAXMAT=10)
      COMMON/COMPOS/STF(MAXMAT,30),ZT(MAXMAT),AT(MAXMAT),RHO(MAXMAT),
     1  VMOL(MAXMAT),IZ(MAXMAT,30),NELEM(MAXMAT)
C  ****  Element data.
      COMMON/CADATA/ATW(99),EPX(99),RA1(99),RA2(99),RA3(99),RA4(99),
     1  RA5(99),RSCR(99),ETA(99),EB(99,30),IFI(99,30),IKS(99,30),
     2  NSHT(99),LASYMB(99)
C
      COMMON/CGRA/X2COH(MAXMAT,241),PDCOH(MAXMAT,241),FLCOH
      COMMON/CGRA00/EE,MM,MOM
      EXTERNAL GRADX1
C
      MM=M
      SUM=0.0D0
      XL=0.0D0
      XU=1.0D-4
      FACT=(1.0D6/XU)**(1.0D0/240.0D0)
      SUM=SUMGA(GRADX1,XL,XU,1.0D-10)
      X2COH(M,1)=XU
      PDCOH(M,1)=SUM
      DO I=2,241
        XL=XU
        XU=XL*FACT
        SUM=SUMGA(GRADX1,XL,XU,1.0D-10)
        X2COH(M,I)=XU
        PDCOH(M,I)=PDCOH(M,I-1)+SUM
      ENDDO
      DO I=1,241
        X2COH(M,I)=LOG(X2COH(M,I))
        PDCOH(M,I)=LOG(PDCOH(M,I))
      ENDDO
      FLCOH=LOG(FACT)
      RETURN
      END
C  *********************************************************************
C                       SUBROUTINE GRATX
C  *********************************************************************
      SUBROUTINE GRATX(E,CS,M)
C
C  Total cross section for Rayleigh (coherent) photon scattering. Born
C  approximation with analytical atomic form factors.
C
C  Input arguments:
C    E ........ photon energy (eV).
C    M ........ material where photons propagate.
C  Output argument:
C    CS ....... coherent total cross section (cm**2/molecule).
C
      use mtpar
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N)
      PARAMETER (REV=5.10998902D5)  ! Electron rest energy (eV)
      PARAMETER (ELRAD=2.817940285D-13)  ! Class. electron radius (cm)
      PARAMETER (PI=3.1415926535897932D0, PIELR2=PI*ELRAD*ELRAD)
C  ****  Composition data.
      PARAMETER (MAXMAT=10)
      COMMON/COMPOS/STF(MAXMAT,30),ZT(MAXMAT),AT(MAXMAT),RHO(MAXMAT),
     1  VMOL(MAXMAT),IZ(MAXMAT,30),NELEM(MAXMAT)
      COMMON/CGRA00/FACTE,MM,MOM
C
      EXTERNAL GRADX
C
      MM=M
      MOM=0
      IZZ=0
      DO I=1,NELEM(M)
        IZZ=MAX(IZZ,IZ(M,I))
      ENDDO
      EC=MIN(E,5.0D5*IZZ)
      FACTE=849.3315D0*(EC/REV)**2
      CS=SUMGA(GRADX,-1.0D0,1.0D0,1.0D-6)
      CS=PIELR2*CS*(EC/E)**2
      RETURN
      END
C  *********************************************************************
C                       FUNCTION GRADX
C  *********************************************************************
      FUNCTION GRADX(CDT)
C
C  Differential x-section for Rayleigh scattering.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N)
      PARAMETER (REV=5.10998902D5)  ! Electron rest energy (eV)
C
      COMMON/CGRA00/FACTE,M,MOM
C
      X2=FACTE*(1.0D0-CDT)
      GRADX=(1.0D0+CDT*CDT)*GRADX1(X2)
      IF(MOM.GT.0) GRADX=GRADX*((1.0D0-CDT)*0.5D0)**MOM
      RETURN
      END
C  *********************************************************************
C                       FUNCTION GRADX1
C  *********************************************************************
      FUNCTION GRADX1(X2)
C
C  Squared molecular form factor (additivity rule).
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N)
      CHARACTER*2 LASYMB
      PARAMETER (SL=137.03599976D0)  ! Speed of light (1/alpha)
C  ****  Composition data.
      PARAMETER (MAXMAT=10)
      COMMON/COMPOS/STF(MAXMAT,30),ZT(MAXMAT),AT(MAXMAT),RHO(MAXMAT),
     1  VMOL(MAXMAT),IZ(MAXMAT,30),NELEM(MAXMAT)
C  ****  Element data.
      COMMON/CADATA/ATW(99),EPX(99),RA1(99),RA2(99),RA3(99),RA4(99),
     1  RA5(99),RSCR(99),ETA(99),EB(99,30),IFI(99,30),IKS(99,30),
     2  NSHT(99),LASYMB(99)
C
      COMMON/CGRA00/FACTE,M,MOM
C
      X=SQRT(X2)
      GRADX1=0.0D0
      DO I=1,NELEM(M)
C  ****  Atomic form factors.
        IZZ=IZ(M,I)
        FA=IZZ*(1.0D0+X2*(RA1(IZZ)+X*(RA2(IZZ)+X*RA3(IZZ))))
     1    /(1.0D0+X2*(RA4(IZZ)+X2*RA5(IZZ)))**2
        IF(IZZ.GT.10.AND.FA.LT.2.0D0) THEN
          PA=(IZZ-0.3125D0)/SL
          PG=SQRT(1.0D0-PA*PA)
          PQ=2.426311D-2*X/PA
          FB=SIN(2.0D0*PG*ATAN2(PQ,1.0D0))/(PG*PQ*(1.0D0+PQ*PQ)**PG)
          FA=MAX(FA,FB)
        ENDIF
        GRADX1=GRADX1+STF(M,I)*FA**2
      ENDDO
      RETURN
      END
C  *********************************************************************
C                       SUBROUTINE GCO
C  *********************************************************************
      SUBROUTINE GCO(E,DE,EP,CDT,ES,CDTS,M,ISHELL)
C
C  Random sampling of incoherent (Compton) scattering of photons. Relat-
C  ivistic impulse approximation with analytical one-electron Compton
C  profiles.
C
C  Input arguments:
C    E ........ incident photon energy (eV).
C    M ........ material where photons propagate.
C  Output argument:
C    DE ....... energy loss (eV)
C    EP ....... energy of the scattered photon (eV)
C    CDT ...... cosine of the polar scattering angle.
C    ES ....... energy of the emitted electron (eV)
C    CDTS ..... polar cosine of direction of the electron.
C    ISHELL ... index of the shell that has been 'ionized'.
C
      use mtpar
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N)
      PARAMETER (REV=5.10998902D5)  ! Electron rest energy (eV)
      PARAMETER (D2=1.4142135623731D0, D1=1.0D0/D2, D12=0.5D0)
C  ****  Composition data.
      PARAMETER (MAXMAT=10)
      COMMON/COMPOS/STF(MAXMAT,30),ZT(MAXMAT),AT(MAXMAT),RHO(MAXMAT),
     1  VMOL(MAXMAT),IZ(MAXMAT,30),NELEM(MAXMAT)
      COMMON/CECUTR/ECUTR(MAXMAT)
C  ****  Compton scattering.
      PARAMETER (NOCO=64)
      COMMON/CGCO/FCO(MAXMAT,NOCO),UICO(MAXMAT,NOCO),FJ0(MAXMAT,NOCO),
     2  KZCO(MAXMAT,NOCO),KSCO(MAXMAT,NOCO),NOSCCO(MAXMAT)
C
      DIMENSION RN(NOCO),PAC(NOCO)
      ! EXTERNALRAND
C
      EK=E/REV
      EK2=EK+EK+1.0D0
      EK3=EK*EK
      EK1=EK3-EK2-1.0D0
      TAUMIN=1.0D0/EK2
      TAUM2=TAUMIN*TAUMIN
      A1=LOG(EK2)
      A2=A1+2.0D0*EK*(1.0D0+EK)*TAUM2
      IF(E.GT.5.0D6) GO TO 4
C
C  ****  Incoherent scattering function for theta=PI.
C
      S0=0.0D0
      DO I=1,NOSCCO(M)
        IF(UICO(M,I).LT.E) THEN
          AUX=E*(E-UICO(M,I))*2.0D0
          PZOMC=FJ0(M,I)*(AUX-REV*UICO(M,I))
     1         /(REV*SQRT(AUX+AUX+UICO(M,I)**2))
          IF(PZOMC.GT.0.0D0) THEN
            RNI=1.0D0-0.5D0*EXP(D12-(D1+D2*PZOMC)**2)
          ELSE
            RNI=0.5D0*EXP(D12-(D1-D2*PZOMC)**2)
          ENDIF
          S0=S0+FCO(M,I)*RNI
        ENDIF
      ENDDO
C
C  ****  Sampling tau.
C
    1 CONTINUE
        status=rmt(r)
      IF(r*A2.LT.A1) THEN
        status=rmt(r)
        TAU=TAUMIN**r
      ELSE
        status=rmt(r)
        TAU=SQRT(1.0D0+r*(TAUM2-1.0D0))
      ENDIF
      CDT1=(1.0D0-TAU)/(EK*TAU)
C  ****  Incoherent scattering function.
      S=0.0D0
      DO I=1,NOSCCO(M)
        IF(UICO(M,I).LT.E) THEN
          AUX=E*(E-UICO(M,I))*CDT1
          PZOMC=FJ0(M,I)*(AUX-REV*UICO(M,I))
     1         /(REV*SQRT(AUX+AUX+UICO(M,I)**2))
          IF(PZOMC.GT.0.0D0) THEN
            RN(I)=1.0D0-0.5D0*EXP(D12-(D1+D2*PZOMC)**2)
          ELSE
            RN(I)=0.5D0*EXP(D12-(D1-D2*PZOMC)**2)
          ENDIF
          S=S+FCO(M,I)*RN(I)
          PAC(I)=S
        ELSE
          PAC(I)=S-1.0D-6
        ENDIF
      ENDDO
C  ****  Rejection function.
      TST=S*(1.0D0+TAU*(EK1+TAU*(EK2+TAU*EK3)))
     1   /(EK3*TAU*(1.0D0+TAU*TAU))
        status=rmt(r)
      IF(r*S0.GT.TST) GO TO 1
      CDT=1.0D0-CDT1
C
C  ****  Target electron shell.
C
    2 CONTINUE
        status=rmt(r)
      TST=S*r
      DO I=1,NOSCCO(M)
        IF(PAC(I).GT.TST) THEN
          ISHELL=I
          GO TO 3
        ENDIF
      ENDDO
      ISHELL=NOSCCO(M)
    3 CONTINUE
C
C  ****  Projected momentum of the target electron.
C
        status=rmt(r)
      A=r*RN(ISHELL)
      IF(A.LT.0.5D0) THEN
        PZOMC=(D1-SQRT(D12-LOG(A+A)))/(D2*FJ0(M,ISHELL))
      ELSE
        PZOMC=(SQRT(D12-LOG(2.0D0-A-A))-D1)/(D2*FJ0(M,ISHELL))
      ENDIF
      IF(PZOMC.LT.-1.0D0) GO TO 2
C
C  ****  F(EP) rejection.
C
      XQC=1.0D0+TAU*(TAU-2.0D0*CDT)
      AF=SQRT(XQC)*(1.0D0+TAU*(TAU-CDT)/XQC)
      IF(AF.GT.0.0D0) THEN
        FPZMAX=1.0D0+AF*0.2D0
      ELSE
        FPZMAX=1.0D0-AF*0.2D0
      ENDIF
      FPZ=1.0D0+AF*MAX(MIN(PZOMC,0.2D0),-0.2D0)
        status=rmt(r)
      IF(r*FPZMAX.GT.FPZ) GO TO 2
C
C  ****  Energy of the scattered photon.
C
      T=PZOMC**2
      B1=1.0D0-T*TAU*TAU
      B2=1.0D0-T*TAU*CDT
      IF(PZOMC.GT.0.0D0) THEN
        EP=E*(TAU/B1)*(B2+SQRT(ABS(B2*B2-B1*(1.0D0-T))))
      ELSE
        EP=E*(TAU/B1)*(B2-SQRT(ABS(B2*B2-B1*(1.0D0-T))))
      ENDIF
      GO TO 6
C
C  ****  No Doppler broadening for E greater than 5 MeV.
C
    4 CONTINUE
        status=rmt(r)
      IF(r*A2.LT.A1) THEN
        status=rmt(r)
        TAU=TAUMIN**r
      ELSE
        status=rmt(r)
        TAU=SQRT(1.0D0+r*(TAUM2-1.0D0))
      ENDIF
C  ****  Rejection function.
      TST=(1.0D0+TAU*(EK1+TAU*(EK2+TAU*EK3)))
     1   /(EK3*TAU*(1.0D0+TAU*TAU))
        status=rmt(r)
      IF(r.GT.TST) GO TO 4
      EP=TAU*E
      CDT=1.0D0-(1.0D0-TAU)/(EK*TAU)
C
C  ****  Target electron shell.
C
        status=rmt(r)
      TST=ZT(M)*r
      S=0.0D0
      DO I=1,NOSCCO(M)
        S=S+FCO(M,I)
        IF(S.GT.TST) THEN
          ISHELL=I
          GO TO 5
        ENDIF
      ENDDO
      ISHELL=NOSCCO(M)
    5 CONTINUE
      IF(EP.GT.E-UICO(M,ISHELL)) GO TO 4
C
    6 CONTINUE
      DE=E-EP
      IF(KSCO(M,ISHELL).LT.10) THEN
        IF(UICO(M,ISHELL).GT.ECUTR(M)) THEN
          ES=DE-UICO(M,ISHELL)
        ELSE
          ES=DE
        ENDIF
      ELSE
        ES=DE
      ENDIF
C
      Q2=E*E+EP*(EP-2.0D0*E*CDT)
      IF(Q2.GT.1.0D-12) THEN
        CDTS=(E-EP*CDT)/SQRT(Q2)
      ELSE
        CDTS=1.0D0
      ENDIF
C
      RETURN
      END
C  *********************************************************************
C                       SUBROUTINE GCOTX
C  *********************************************************************
      SUBROUTINE GCOTX(E,CS,M)
C
C  Total cross section for incoherent (Compton) scattering. Relativistic
C  Impulse approximation with analytical Compton profiles.
C
C  Input arguments:
C    E ........ photon energy (eV).
C    M ........ material where photons propagate.
C  Output argument:
C    CS ....... incoherent total cross section (cm**2/molecule).
C
      use mtpar
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N)
      PARAMETER (MAXMAT=10)
      PARAMETER (REV=5.10998902D5)  ! Electron rest energy (eV)
      PARAMETER (ELRAD=2.817940285D-13)  ! Class. electron radius (cm)
      PARAMETER (PI=3.1415926535897932D0, PIELR2=PI*ELRAD*ELRAD)
C  ****  Compton scattering.
      PARAMETER (NOCO=64)
      COMMON/CGCO/FCO(MAXMAT,NOCO),UICO(MAXMAT,NOCO),FJ0(MAXMAT,NOCO),
     2  KZCO(MAXMAT,NOCO),KSCO(MAXMAT,NOCO),NOSCCO(MAXMAT)
C
      COMMON/CGCO00/EE,EP,MM
      EXTERNAL GCODX
C
      IF(E.LT.5.0D6) THEN
        EE=E
        MM=M
        CS=SUMGA(GCODX,-1.0D0,1.0D0,1.0D-5)
      ELSE
C  ****  Klein-Nishina total cross section.
        EK=E/REV
        EK3=EK*EK
        EK2=1.0D0+EK+EK
        EK1=EK3-EK2-1.0D0
        T0=1.0D0/(1.0D0+EK+EK)
        CSL=0.5D0*EK3*T0*T0+EK2*T0+EK1*LOG(T0)-1.0D0/T0
        CS=0.0D0
        DO 1 I=1,NOSCCO(M)
          TAU=(E-UICO(M,I))/E
          IF(TAU.LT.T0) GO TO 1
          CSU=0.5D0*EK3*TAU*TAU+EK2*TAU+EK1*LOG(TAU)-1.0D0/TAU
          CS=CS+FCO(M,I)*(CSU-CSL)
    1   CONTINUE
        CS=PIELR2*CS/(EK*EK3)
      ENDIF
      RETURN
      END
C  *********************************************************************
C                        FUNCTION GCODX
C  *********************************************************************
      FUNCTION GCODX(CDT)
C
C  Single differential cross section for photon Compton scattering, dif-
C  ferential in the direction of the scattered photon only. Evaluated
C  from the incoherent scattering fucntion.
C
C  The energy E of the primary photon is entered through common CGCO00.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N)
      PARAMETER (MAXMAT=10)
      PARAMETER (REV=5.10998902D5)  ! Electron rest energy (eV)
      PARAMETER (ELRAD=2.817940285D-13)  ! Class. electron radius (cm)
      PARAMETER (PI=3.1415926535897932D0, PIELR2=PI*ELRAD*ELRAD)
      PARAMETER (D2=1.4142135623731D0, D1=1.0D0/D2, D12=0.5D0)
C  ****  Compton scattering.
      PARAMETER (NOCO=64)
      COMMON/CGCO/FCO(MAXMAT,NOCO),UICO(MAXMAT,NOCO),FJ0(MAXMAT,NOCO),
     2  KZCO(MAXMAT,NOCO),KSCO(MAXMAT,NOCO),NOSCCO(MAXMAT)
C
      COMMON/CGCO00/E,EP,M
C
      CDT1=1.0D0-CDT
C  ****  Energy of the Compton line.
      EOEC=1.0D0+(E/REV)*CDT1
      ECOE=1.0D0/EOEC
C  ****  Incoherent scattering function (analytical profile).
      SIA=0.0D0
      DO 1 I=1,NOSCCO(M)
        IF(E.LT.UICO(M,I)) GO TO 1
        AUX=E*(E-UICO(M,I))*CDT1
        PZIMAX=(AUX-REV*UICO(M,I))/(REV*SQRT(AUX+AUX+UICO(M,I)**2))
        X=FJ0(M,I)*PZIMAX
        IF(X.GT.0.0D0) THEN
          SIAP=1.0D0-0.5D0*EXP(D12-(D1+D2*X)**2)
        ELSE
          SIAP=0.5D0*EXP(D12-(D1-D2*X)**2)
        ENDIF
        SIA=SIA+FCO(M,I)*SIAP
    1 CONTINUE
C  ****  Klein-Nishina X-factor.
      XKN=EOEC+ECOE-1.0D0+CDT*CDT
C  ****  Differential cross section.
      GCODX=PIELR2*ECOE**2*XKN*SIA
      RETURN
      END
C  *********************************************************************
C                       SUBROUTINE GPH
C  *********************************************************************
      SUBROUTINE GPH(ES,IZZ,ISH)
C
C  Simulation of photoelectric absorption in material M.
C
C  Output arguments:
C    ES .... kinetic energy of the photoelectron.
C    IZZ ... atomic number of the element where absorption has ocurred.
C    ISH ... atomic electron shell that has been ionized.
C
C  NOTE: JUMP uses a photoelectric cross section that is slightly larger
C  than its 'true' value. To correct for this, the photon is allowed to
C  'survive' a photoelectric event. Survival of the photon is flagged by
C  setting IZZ=0, ISH=0, ES=0.0D0 (the energy E of the photon is kept
C  unaltered.
C
      use mtpar
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N)
      CHARACTER*2 LASYMB
C  ****  Main-PENELOPE common.
      COMMON/TRACK/E,X,Y,Z,U,V,W,WGHT,KPAR,IBODY,M,ILB(5)
C  ****  Simulation parameters.
      PARAMETER (MAXMAT=10)
      COMMON/CSIMPA/EABS(3,MAXMAT),C1(MAXMAT),C2(MAXMAT),WCC(MAXMAT),
     1  WCR(MAXMAT)
      COMMON/CECUTR/ECUTR(MAXMAT)
C  ****  Energy grid and interpolation constants for the current energy.
      PARAMETER (NEGP=200)
      COMMON/CEGRID/EL,EU,ET(NEGP),DLEMP(NEGP),DLEMP1,DLFC,
     1  XEL,XE,XEK,KE
C  ****  Element data.
      COMMON/CADATA/ATW(99),EPX(99),RA1(99),RA2(99),RA3(99),RA4(99),
     1  RA5(99),RSCR(99),ETA(99),EB(99,30),IFI(99,30),IKS(99,30),
     2  NSHT(99),LASYMB(99)
C  ****  Composition data.
      COMMON/COMPOS/STF(MAXMAT,30),ZT(MAXMAT),AT(MAXMAT),RHO(MAXMAT),
     1  VMOL(MAXMAT),IZ(MAXMAT,30),NELEM(MAXMAT)
C  ****  Photon simulation tables.
      COMMON/CGIMFP/SGRA(MAXMAT,NEGP),SGCO(MAXMAT,NEGP),
     1  SGPH(MAXMAT,NEGP),SGPP(MAXMAT,NEGP),SGAUX(MAXMAT,NEGP)
C  ****  Photoelectric cross sections.
      PARAMETER (NTP=8000)
      COMMON/CGPH00/EPH(NTP),XPH(NTP,10),IPHF(99),IPHL(99),NPHS(99),NCUR
      DIMENSION ACP(35),IP(35)
C
      PARAMETER (PI=3.1415926535897932D0, TWOPI=PI+PI)
C
      ! EXTERNALRAND
C
C  ****  Partial attenuation coefficients.
C
      PTOT=0.0D0
      DO IEL=1,NELEM(M)
        IZZ=IZ(M,IEL)
C  ****  Binary search.
        I=IPHF(IZZ)
        IU=IPHL(IZZ)
    1   IT=(I+IU)/2
        IF(XEL.GT.EPH(IT)) THEN
          I=IT
        ELSE
          IU=IT
        ENDIF
        IF(IU-I.GT.1) GO TO 1
C
        IP(IEL)=I
        DEE=EPH(I+1)-EPH(I)
        IF(DEE.GT.1.0D-15) THEN
          PCSL=XPH(I,1)+(XPH(I+1,1)-XPH(I,1))*(XEL-EPH(I))/DEE
        ELSE
          PCSL=XPH(I,1)
        ENDIF
        PTOT=PTOT+STF(M,IEL)*EXP(PCSL)
        ACP(IEL)=PTOT
      ENDDO
      IF(PTOT*VMOL(M).GT.SGPH(M,KE)) THEN
        WRITE(26,*) 'WARNING: SGPH is less than the actual mac.'
      ENDIF
C
C  ****  Sample the active element.
C
        status=rmt(r)
      TST=r*SGPH(M,KE)/VMOL(M)
      DO IEL=1,NELEM(M)
        IF(ACP(IEL).GT.TST) THEN
          IELAC=IEL
          IZZ=IZ(M,IEL)
          GO TO 2
        ENDIF
      ENDDO
C  ****  Delta interaction. Introduced to correct for the use an upper
C        bound of the photoelectric attenuation coefficient.
      IZZ=0  ! Flags delta interactions.
      ISH=0
      ES=0.0D0
      RETURN
C
    2 CONTINUE
      IF(NPHS(IZZ).EQ.1) THEN
        ISH=1
        GO TO 3
      ELSE
C  ****  Selection of the active shell.
        I=IP(IELAC)
        DEE=EPH(I+1)-EPH(I)
        PIS=0.0D0
        IF(DEE.GT.1.0D-15) THEN
          PTOT=EXP(XPH(I,1)+(XPH(I+1,1)-XPH(I,1))*(XEL-EPH(I))/DEE)
          DO IS=1,NPHS(IZZ)
            J=IS+1
            PCSL=XPH(I,J)+(XPH(I+1,J)-XPH(I,J))*(XEL-EPH(I))/DEE
            PIS=PIS+EXP(PCSL)
            ACP(IS)=PIS
          ENDDO
        ELSE
          PTOT=EXP(XPH(I,1))
          DO IS=1,NPHS(IZZ)
            PIS=PIS+EXP(XPH(I,IS+1))
            ACP(IS)=PIS
          ENDDO
        ENDIF
        status=rmt(r)
        TST=r*PTOT
        DO IS=1,NPHS(IZZ)
          IF(ACP(IS).GT.TST) THEN
            ISH=IS
            GO TO 3
          ENDIF
        ENDDO
        ISH=10
      ENDIF
C
C  ****  Photoelectron emission.
C
    3 CONTINUE
      IF(ISH.LT.10) THEN
        EBB=EB(IZZ,ISH)
        IF(EBB.GT.ECUTR(M)) THEN
          ES=E-EBB
        ELSE
          ES=E
          ISH=10
        ENDIF
      ELSE
        ES=E
      ENDIF
C
      RETURN
      END
C  *********************************************************************
C                       SUBROUTINE SAUTER
C  *********************************************************************
      SUBROUTINE SAUTER(ES,CDTS)
C
C  Random sampling of the initial direction of photoelectrons from the
C  Sauter distribution.
C
      use mtpar
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N)
      PARAMETER (REV=5.10998902D5)  ! Electron rest energy (eV)
      ! EXTERNALRAND
C
      IF(ES.GT.1.0D9) THEN
        CDTS=1.0D0
        RETURN
      ENDIF
      GAM=1.0D0+ES/REV
      GAM2=GAM*GAM
      BETA=SQRT((GAM2-1.0D0)/GAM2)
      AC=1.0D0/BETA-1.0D0
      A1=0.5D0*BETA*GAM*(GAM-1.0D0)*(GAM-2.0D0)
      A2=AC+2.0D0
      GTMAX=2.0D0*(A1+1.0D0/AC)
    1 CONTINUE
        status=rmt(r)
      RU=r
      TSAM=2.0D0*AC*(2.0D0*RU+A2*SQRT(RU))/(A2*A2-4.0D0*RU)
      GTR=(2.0D0-TSAM)*(A1+1.0D0/(AC+TSAM))
        status=rmt(r)
      IF(r*GTMAX.GT.GTR) GO TO 1
      CDTS=1.0D0-TSAM
      RETURN
      END
C  *********************************************************************
C                       SUBROUTINE GPHR
C  *********************************************************************
      SUBROUTINE GPHR(M,IRD,IWR,INFO)
C
C  This subroutine reads photoelectric cross sections of the elements in
C  material M and prepares simulation tables.
C
C  NOTE: The array SGPH(M,IE) defines a piecewise constant function that
C  is larger than the actual photoelectric cross section. SGPH(M,IE) is
C  defined as the largest value of the photoelectric x-section in the
C  energy interval from ET(IE) to ET(IE+1). The photon mean free path is
C  sampled by using this 'augmented' cross section and, to compensate
C  for this, the photon survives (i.e., it is not absorbed) with a prob-
C  ability such that the 'exact' photoelectric attenuation coefficient
C  is reproduced. This trick allows subroutine JUMP to disregard the
C  existence of absorption edges in the photoelectric x-section and to
C  perform the simulation of photons faster.
C
      use mtpar
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N)
      CHARACTER*2 LASYMB
C  ****  Energy grid and interpolation constants for the current energy.
      PARAMETER (NEGP=200)
      COMMON/CEGRID/EL,EU,ET(NEGP),DLEMP(NEGP),DLEMP1,DLFC,
     1  XEL,XE,XEK,KE
C  ****  Composition data.
      PARAMETER (MAXMAT=10)
      COMMON/COMPOS/STF(MAXMAT,30),ZT(MAXMAT),AT(MAXMAT),RHO(MAXMAT),
     1  VMOL(MAXMAT),IZ(MAXMAT,30),NELEM(MAXMAT)
      COMMON/CECUTR/ECUTR(MAXMAT)
C  ****  Element data.
      COMMON/CADATA/ATW(99),EPX(99),RA1(99),RA2(99),RA3(99),RA4(99),
     1  RA5(99),RSCR(99),ETA(99),EB(99,30),IFI(99,30),IKS(99,30),
     2  NSHT(99),LASYMB(99)
C  ****  Photon simulation tables.
      COMMON/CGIMFP/SGRA(MAXMAT,NEGP),SGCO(MAXMAT,NEGP),
     1  SGPH(MAXMAT,NEGP),SGPP(MAXMAT,NEGP),SGAUX(MAXMAT,NEGP)
C  ****  Photoelectric cross sections.
      PARAMETER (NTP=8000)
      COMMON/CGPH00/EPH(NTP),XPH(NTP,10),IPHF(99),IPHL(99),NPHS(99),NCUR
      PARAMETER (NDIM=1500)
      COMMON/CGPH01/ER(NDIM),XSR(NDIM),NPHD
      DIMENSION XGPHR(NDIM,10),X1(NDIM),Y1(NDIM),X2(NDIM),Y2(NDIM)
C
C  ************  Read element x-section tables
C
      DO IEL=1,NELEM(M)
        READ(IRD,1001) IZZ,NSHR,NDATA
 1001   FORMAT(41X,I3,11X,I3,10X,I4)
        IF(INFO.GE.2) WRITE(IWR,2001) IZZ,NSHR,NDATA
 2001   FORMAT(/1X,'***  Photoelectric cross sections,  IZ =',I3,
     1    ',  NSHELL =',I3,',  NDATA =',I4)
        IF(IZZ.NE.IZ(M,IEL)) STOP 'GPHR. Corrupt material data file.'
        IF(NDATA.GT.NDIM) STOP 'GPHR. Too many data points.'
        IF(NSHR.GT.9) STOP 'GPHR. Too many shells.'
        DO IE=1,NDATA
          READ(IRD,*) ER(IE),(XGPHR(IE,IS),IS=1,NSHR+1)
        ENDDO
C
C  ****  Remove shells with ionization energies less than 100 eV.
C
        IF(NSHR.GT.1) THEN
          NSHA=NSHR
          DO IS=NSHA,1,-1
            IF(EB(IZZ,IS).LT.100.0D0) THEN
              NSHR=NSHR-1
            ELSE
              GO TO 1
            ENDIF
          ENDDO
          IF(NSHR.LT.1) NSHR=1
        ENDIF
    1   CONTINUE
C
        IF(INFO.GE.2) THEN
          IF(NSHR.EQ.1) THEN
            WRITE(IWR,2102)
          ELSE IF(NSHR.EQ.2) THEN
            WRITE(IWR,2202)
          ELSE IF(NSHR.EQ.3) THEN
            WRITE(IWR,2302)
          ELSE IF(NSHR.EQ.4) THEN
            WRITE(IWR,2402)
          ELSE IF(NSHR.EQ.5) THEN
            WRITE(IWR,2502)
          ELSE IF(NSHR.EQ.6) THEN
            WRITE(IWR,2602)
          ELSE IF(NSHR.EQ.7) THEN
            WRITE(IWR,2702)
          ELSE IF(NSHR.EQ.8) THEN
            WRITE(IWR,2802)
          ELSE
            WRITE(IWR,2902)
          ENDIF
          DO IE=1,NDATA
            WRITE(IWR,'(1P,26E12.5)') ER(IE),(XGPHR(IE,IS),IS=1,NSHR+1)
          ENDDO
        ENDIF
 2102   FORMAT(/3X,'Energy',5X,'CS-total',6X,'CS-K',
     1    /4X,'(eV)',2X,2(5X,'(cm**2)'),/1X,36('-'))
 2202   FORMAT(/3X,'Energy',5X,'CS-total',6X,'CS-K',8X,'CS-L1',7X,
     1    /4X,'(eV)',2X,3(5X,'(cm**2)'),/1X,48('-'))
 2302   FORMAT(/3X,'Energy',5X,'CS-total',6X,'CS-K',8X,'CS-L1',7X,
     1    'CS-L2',/4X,'(eV)',2X,4(5X,'(cm**2)'),/1X,60('-'))
 2402   FORMAT(/3X,'Energy',5X,'CS-total',6X,'CS-K',8X,'CS-L1',7X,
     1    'CS-L2',7X,'CS-L3',
     2    /4X,'(eV)',2X,5(5X,'(cm**2)'),/1X,72('-'))
 2502   FORMAT(/3X,'Energy',5X,'CS-total',6X,'CS-K',8X,'CS-L1',7X,
     1    'CS-L2',7X,'CS-L3',7X,'CS-M1',
     2    /4X,'(eV)',2X,6(5X,'(cm**2)'),/1X,84('-'))
 2602   FORMAT(/3X,'Energy',5X,'CS-total',6X,'CS-K',8X,'CS-L1',7X,
     1    'CS-L2',7X,'CS-L3',7X,'CS-M1',7X,'CS-M2',
     2    /4X,'(eV)',2X,7(5X,'(cm**2)'),/1X,96('-'))
 2702   FORMAT(/3X,'Energy',5X,'CS-total',6X,'CS-K',8X,'CS-L1',7X,
     1    'CS-L2',7X,'CS-L3',7X,'CS-M1',7X,'CS-M2',7X,'CS-M3',
     2    /4X,'(eV)',2X,8(5X,'(cm**2)'),/1X,108('-'))
 2802   FORMAT(/3X,'Energy',5X,'CS-total',6X,'CS-K',8X,'CS-L1',7X,
     1    'CS-L2',7X,'CS-L3',7X,'CS-M1',7X,'CS-M2',7X,'CS-M3',7X,
     2    'CS-M4',/4X,'(eV)',2X,9(5X,'(cm**2)'),
     3    /1X,120('-'))
 2902   FORMAT(/3X,'Energy',5X,'CS-total',6X,'CS-K',8X,'CS-L1',7X,
     1    'CS-L2',7X,'CS-L3',7X,'CS-M1',7X,'CS-M2',7X,'CS-M3',7X,
     2    'CS-M4',7X,'CS-M5',/4X,'(eV)',2X,10(5X,'(cm**2)'),
     3    /1X,132('-'))
C
        IF(NPHS(IZZ).EQ.0) THEN
          IPHF(IZZ)=NCUR+1
          IF(NCUR+NDATA.GT.NTP) THEN
            WRITE(IWR,*) 'Insufficient memory storage in GPHR.'
            WRITE(IWR,*) 'Increase the value of the parameter NTP to',
     1        NCUR+NDATA
            WRITE(26,*) 'Insufficient memory storage in GPHR.'
            WRITE(26,*) 'Increase the value of the parameter NTP to',
     1        NCUR+NDATA
            STOP 'GPHR. Insufficient memory storage.'
          ENDIF
          DO IE=1,NDATA
            IC=NCUR+IE
            EPH(IC)=LOG(ER(IE))
            DO IS=1,NSHR+1
              XPH(IC,IS)=LOG(MAX(XGPHR(IE,IS),1.0D-35))
            ENDDO
          ENDDO
          NCUR=NCUR+NDATA
          IPHL(IZZ)=NCUR
          NPHS(IZZ)=NSHR
        ENDIF
      ENDDO
C
C  ****  Total photoelectric attenuation coefficient.
C
      IZZ=IZ(M,1)
      IST=IPHF(IZZ)
      LST=IPHL(IZZ)
      NPHD=0
      DO I=IST,LST
        NPHD=NPHD+1
        ER(NPHD)=EXP(EPH(I))
        XSR(NPHD)=STF(M,1)*EXP(XPH(I,1))
      ENDDO
      IF(NELEM(M).GT.1) THEN
        DO IEL=2,NELEM(M)
          N1=NPHD
          DO I=1,N1
            X1(I)=ER(I)
            Y1(I)=XSR(I)
          ENDDO
          IZZ=IZ(M,IEL)
          IST=IPHF(IZZ)
          LST=IPHL(IZZ)
          N2=0
          DO I=IST,LST
            N2=N2+1
            X2(N2)=EXP(EPH(I))
            Y2(N2)=STF(M,IEL)*EXP(XPH(I,1))
          ENDDO
          CALL MERGE2(X1,Y1,X2,Y2,ER,XSR,N1,N2,NPHD)
        ENDDO
      ENDIF
C
C  ****  Total photoelectric cross section at the simulation grid points
C        (slightly increased to simplify the interpolation).
C
      DO I=1,NPHD
        X1(I)=LOG(ER(I))
        Y1(I)=LOG(XSR(I)*VMOL(M))
      ENDDO
C
      DO IE=1,NEGP-1
        EG1=DLEMP(IE)
        CALL FINDI(X1,EG1,NPHD,I1)
        IF(I1.EQ.NPHD) I1=NPHD-1
        DX=X1(I1+1)-X1(I1)
        IF(DX.GT.1.0D-15) THEN
          F1=Y1(I1)+(Y1(I1+1)-Y1(I1))*(EG1-X1(I1))/DX
        ELSE
          F1=Y1(I1)
        ENDIF
        EG2=DLEMP(IE+1)
        CALL FINDI(X1,EG2,NPHD,I2)
        IF(I2.EQ.NPHD) I2=NPHD-1
        DX=X1(I2+1)-X1(I2)
        IF(DX.GT.1.0D-15) THEN
          F2=Y1(I2)+(Y1(I2+1)-Y1(I2))*(EG2-X1(I2))/DX
        ELSE
          F2=Y1(I2)
        ENDIF
C  ****  To avoid interpolating the attenuation coefficient tables, we
C        replace the photoelectric inverse mean free path (imfp) in each
C        energy interval by its upper bound. The increase of the imfp
C        is interpreted as the imfp of delta interactions.
        FM=MAX(F1,F2)
        IF(I1+1.LE.I2) THEN
          DO I=I1+1,I2
           FM=MAX(FM,Y1(I))
          ENDDO
        ENDIF
        SGPH(M,IE)=EXP(FM)
      ENDDO
      SGPH(M,NEGP)=SGPH(M,NEGP-1)
C
      RETURN
      END
C  *********************************************************************
C                       SUBROUTINE GPH0
C  *********************************************************************
      SUBROUTINE GPH0
C
C  This subroutine sets all variables in common /CGPH00/ to zero.
C  It has to be invoked before reading the first material definition
C  file.
C
      use mtpar
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N)
      CHARACTER*2 LASYMB
C  ****  Element data.
      COMMON/CADATA/ATW(99),EPX(99),RA1(99),RA2(99),RA3(99),RA4(99),
     1  RA5(99),RSCR(99),ETA(99),EB(99,30),IFI(99,30),IKS(99,30),
     2  NSHT(99),LASYMB(99)
C  ****  Photoelectric cross sections.
      PARAMETER (NTP=8000)
      COMMON/CGPH00/EPH(NTP),XPH(NTP,10),IPHF(99),IPHL(99),NPHS(99),NCUR
C
      DO I=1,99
        NPHS(I)=0
        IPHF(I)=0
        IPHL(I)=0
      ENDDO
C
      DO I=1,NTP
        EPH(I)=0.0D0
        DO J=1,10
          XPH(I,J)=1.0D-35
        ENDDO
      ENDDO
      NCUR=0
      RETURN
      END
C  *********************************************************************
C                       SUBROUTINE GPHW
C  *********************************************************************
      SUBROUTINE GPHW(M,IWR)
C
C  This subroutine generates the table of photoelectric cross sections
C  for photons in material M and writes it on the material data file.
C  Data are read from the files 'pdgphZZ.p05'.
C
      use mtpar
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N)
C
      CHARACTER*12 FILEN
      CHARACTER*1 LDIG(10),LDIG1,LDIG2
      DATA LDIG/'0','1','2','3','4','5','6','7','8','9'/
C  ****  Composition data.
      PARAMETER (MAXMAT=10)
      COMMON/COMPOS/STF(MAXMAT,30),ZT(MAXMAT),AT(MAXMAT),RHO(MAXMAT),
     1  VMOL(MAXMAT),IZ(MAXMAT,30),NELEM(MAXMAT)
C
      PARAMETER (NPHM=400)
      DIMENSION XS(10),E0(NPHM),XS0(NPHM,10)
C
      DO IEL=1,NELEM(M)
        IZZ=IZ(M,IEL)
        NLD=IZZ
        NLD1=NLD-10*(NLD/10)
        NLD2=(NLD-NLD1)/10
        LDIG1=LDIG(NLD1+1)
        LDIG2=LDIG(NLD2+1)
        FILEN='pdgph'//LDIG2//LDIG1//'.p05'
        OPEN(3,FILE=FILEN)
        READ(3,*) IZZZ,NSHR
        IF(IZZZ.NE.IZZ) STOP 'GPHW. Corrupt file.'
        IF(NSHR.GT.9) STOP 'GPHW. Too many shells.'
        NPTAB=0
        DO IE=1,1000
          READ(3,*,END=1) ER,(XS(IS),IS=1,NSHR+1)
          IF(ER.GT.49.9D0.AND.ER.LT.1.01D9) THEN
            NPTAB=NPTAB+1
            E0(NPTAB)=ER
            DO IS=1,NSHR+1
              XS0(NPTAB,IS)=XS(IS)*1.0D-24
            ENDDO
          ENDIF
        ENDDO
    1   CONTINUE
        CLOSE(3)
        WRITE(IWR,2001) IZZ,NSHR,NPTAB
 2001 FORMAT(1X,'***  Photoelectric cross sections,  IZ =',I3,
     1  ',  NSHELL =',I3,',  NDATA =',I4)
        DO IE=1,NPTAB
          WRITE(IWR,2002) E0(IE),(XS0(IE,IS),IS=1,NSHR+1)
        ENDDO
 2002 FORMAT(1P,26E12.5)
      ENDDO
      RETURN
      END
C  *********************************************************************
C                       SUBROUTINE GPP
C  *********************************************************************
      SUBROUTINE GPP(EE,CDTE,EP,CDTP)
C
C  Random sampling of electron-positron pair production by photons.
C  Bethe-Heitler differential cross section.
C
C  Output values:
C    EE .....  kinetic energy of the electron.
C    CDTE ...  polar direction cosine of the electron.
C    EP .....  kinetic energy of the positron.
C    CDTP ...  polar direction cosine of the positron.
C
      use mtpar
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N)
      PARAMETER (REV=5.10998902D5)  ! Electron rest energy (eV)
      PARAMETER (SL=137.03599976D0)  ! Speed of light (1/alpha)
      PARAMETER (PI=3.1415926535897932D0, TWOPI=PI+PI, TREV=2.0D0*REV)
C  ****  Main-PENELOPE common.
      COMMON/TRACK/E,X,Y,Z,U,V,W,WGHT,KPAR,IBODY,M,ILB(5)
C  ****  Simulation parameters.
      PARAMETER (MAXMAT=10)
      COMMON/CSIMPA/EABS(3,MAXMAT),C1(MAXMAT),C2(MAXMAT),WCC(MAXMAT),
     1  WCR(MAXMAT)
C  ****  Pair-production cross section parameters.
      COMMON/CGPP00/ZEQPP(MAXMAT),F0(MAXMAT,2),BCB(MAXMAT)
C
      ! EXTERNALRAND
C
      EKI=REV/E
      IF(E.LT.1.1D6) THEN
        status=rmt(r)
        EPS=EKI+(1.0D0-2.0D0*EKI)*r
        GO TO 3
      ENDIF
C  ****  Low-energy and Coulomb corrections.
      ALZ=ZEQPP(M)/SL
      T=SQRT(2.0D0*EKI)
      F00=(-1.774D0-1.210D1*ALZ+1.118D1*ALZ*ALZ)*T
     1    +(8.523D0+7.326D1*ALZ-4.441D1*ALZ*ALZ)*T**2
     2    -(1.352D1+1.211D2*ALZ-9.641D1*ALZ*ALZ)*T**3
     3    +(8.946D0+6.205D1*ALZ-6.341D1*ALZ*ALZ)*T**4
      G0=F0(M,2)+F00
      BMIN=4.0D0*EKI/BCB(M)
      CALL SCHIFF(BMIN,G1,G2)
      G1MIN=G1+G0
      G2MIN=G2+G0
      XR=0.5D0-EKI
      A1=6.666666666666666D-1*G1MIN*XR**2
      P1=A1/(A1+G2MIN)
C  ****  Random sampling of EPS.
    1 CONTINUE
        status=rmt(r)
      IF(r.GT.P1) GO TO 2
        status=rmt(r)
      RU2M1=2.0D0*r-1.0D0
      IF(RU2M1.LT.0.0D0) THEN
        EPS=0.5D0-XR*ABS(RU2M1)**3.333333333333333D-1
      ELSE
        EPS=0.5D0+XR*RU2M1**3.333333333333333D-1
      ENDIF
      B=EKI/(BCB(M)*EPS*(1.0D0-EPS))
      CALL SCHIFF(B,G1,G2)
      G1=MAX(G1+G0,0.0D0)
        status=rmt(r)
      IF(r*G1MIN.GT.G1) GO TO 1
      GO TO 3
    2 CONTINUE
        status=rmt(r)
      EPS=EKI+2.0D0*XR*r
      B=EKI/(BCB(M)*EPS*(1.0D0-EPS))
      CALL SCHIFF(B,G1,G2)
      G2=MAX(G2+G0,0.0D0)
        status=rmt(r)
      IF(r*G2MIN.GT.G2) GO TO 1
    3 CONTINUE
C  ****  Electron.
      EE=EPS*E-REV
      IF(EE.GT.EABS(1,M)) THEN
        status=rmt(r)
        CDTE=2.0D0*r-1.0D0
        A1=EE+REV
        A2=SQRT(EE*(EE+TREV))
        CDTE=(CDTE*A1+A2)/(A1+CDTE*A2)
      ELSE
        CDTE=1.0D0
      ENDIF
C  ****  Positron.
      EP=(1.0D0-EPS)*E-REV
      IF(EP.GT.EABS(3,M)) THEN
        status=rmt(r)
        CDTP=2.0D0*r-1.0D0
        A1=EP+REV
        A2=SQRT(EP*(EP+TREV))
        CDTP=(CDTP*A1+A2)/(A1+CDTP*A2)
      ELSE
        CDTP=1.0D0
      ENDIF
      RETURN
      END
C  *********************************************************************
C                       SUBROUTINE SCHIFF
C  *********************************************************************
      SUBROUTINE SCHIFF(B,G1,G2)
C
C  Screening functions F1(B) and F2(B) in the Bethe-Heitler differential
C  cross section for pair production.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N)
      PARAMETER (PI=3.1415926535897932D0, TWOPI=PI+PI)
      B2=B*B
      F1=2.0D0-2.0D0*LOG(1.0D0+B2)
      F2=F1-6.666666666666666D-1
      IF(B.LT.1.0D-10) THEN
        F1=F1-TWOPI*B
      ELSE
        A0=4.0D0*B*ATAN2(1.0D0,B)
        F1=F1-A0
        F2=F2+2.0D0*B2*(4.0D0-A0-3.0D0*LOG((1.0D0+B2)/B2))
      ENDIF
      G1=0.5D0*(3.0D0*F1-F2)
      G2=0.25D0*(3.0D0*F1+F2)
      RETURN
      END
C  *********************************************************************
C                       SUBROUTINE GPP0
C  *********************************************************************
      SUBROUTINE GPP0(M)
C
C  Initialization of the sampling algorithm for electron-positron pair
C  production by photons in material M. Bethe-Heitler differential cross
C  section.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N)
      PARAMETER (REV=5.10998902D5)  ! Electron rest energy (eV)
      PARAMETER (SL=137.03599976D0)  ! Speed of light (1/alpha)
C  ****  Element data.
      CHARACTER*2 LASYMB
      COMMON/CADATA/ATW(99),EPX(99),RA1(99),RA2(99),RA3(99),RA4(99),
     1  RA5(99),RSCR(99),ETA(99),EB(99,30),IFI(99,30),IKS(99,30),
     2  NSHT(99),LASYMB(99)
C  ****  Composition data.
      PARAMETER (MAXMAT=10)
      COMMON/COMPOS/STF(MAXMAT,30),ZT(MAXMAT),AT(MAXMAT),RHO(MAXMAT),
     1  VMOL(MAXMAT),IZ(MAXMAT,30),NELEM(MAXMAT)
C  ****  Pair-production cross section parameters.
      COMMON/CGPP00/ZEQPP(MAXMAT),F0(MAXMAT,2),BCB(MAXMAT)
C
C  ***  Effective atomic number.
C
      FACT=0.0D0
      DO I=1,NELEM(M)
        IZZ=IZ(M,I)
        FACT=FACT+IZZ*ATW(IZZ)*STF(M,I)
      ENDDO
      ZEQPP(M)=FACT/AT(M)
      IZZ=ZEQPP(M)+0.25D0
      IF(IZZ.LE.0) IZZ=1
      IF(IZZ.GT.99) IZZ=99
C  ****  DBM Coulomb correction.
      ALZ=ZEQPP(M)/SL
      A=ALZ*ALZ
      FC=A*(0.202059D0-A*(0.03693D0-A*(0.00835D0-A*(0.00201D0-A*
     1 (0.00049D0-A*(0.00012D0-A*0.00003D0)))))+1.0D0/(A+1.0D0))
C  ****  Screening functions and low-energy correction.
      BCB(M)=2.0D0/RSCR(IZZ)
      F0(M,1)=4.0D0*LOG(RSCR(IZZ))
      F0(M,2)=F0(M,1)-4.0D0*FC
      RETURN
      END
C  *********************************************************************
C                       SUBROUTINE GPPW
C  *********************************************************************
      SUBROUTINE GPPW(EIT,XG0,NPTAB,M)
C
C  This subroutine generates a table of electron-positron pair produc-
C  tion cross sections for photons in material M. Data are read from the
C  files 'pdgppZZ.p05'.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N)
C
      CHARACTER*12 FILEN
      CHARACTER*1 LDIG(10),LDIG1,LDIG2
      DATA LDIG/'0','1','2','3','4','5','6','7','8','9'/
C  ****  Composition data.
      PARAMETER (MAXMAT=10)
      COMMON/COMPOS/STF(MAXMAT,30),ZT(MAXMAT),AT(MAXMAT),RHO(MAXMAT),
     1  VMOL(MAXMAT),IZ(MAXMAT,30),NELEM(MAXMAT)
C  ****
      PARAMETER (NEGP=200)
      DIMENSION EIT(NEGP),XG0(NEGP)
C
C  ****  Building the cross section table.
C
      DO I=1,NEGP
        XG0(I)=0.0D0
      ENDDO
C
      DO IEL=1,NELEM(M)
        IZZ=IZ(M,IEL)
        WGHT=STF(M,IEL)*1.0D-24
        NLD=IZZ
        NLD1=NLD-10*(NLD/10)
        NLD2=(NLD-NLD1)/10
        LDIG1=LDIG(NLD1+1)
        LDIG2=LDIG(NLD2+1)
        FILEN='pdgpp'//LDIG2//LDIG1//'.p05'
        OPEN(3,FILE=FILEN)
        READ(3,*) IZZZ
        IF(IZZZ.NE.IZZ) STOP 'GPPW. Corrupt file.'
        DO I=1,NEGP
          READ(3,*,END=1) EIT(I),XG0P
          XG0(I)=XG0(I)+WGHT*XG0P
          NPTAB=I
          IF(EIT(I).GT.0.999D9) GO TO 1
        ENDDO
    1   CONTINUE
        CLOSE(3)
      ENDDO
      RETURN
      END
C  *********************************************************************
C                       SUBROUTINE RELAX
C  *********************************************************************
      SUBROUTINE RELAX(IZ,IS)
C
C  This subroutine simulates the relaxation of a singly ionized atom of
C  the element IZ with a vacancy in the IS shell (the K shell or an L
C  or M subshell). This initial vacancy is filled by electrons from
C  outer shells through radiative and non-radiative transitions, which
C  may produce additional vacancies.
C
C  We use the following notation to designate the possible transitions:
C  *  Radiative: IS0-IS1 (an electron from the IS1 shell fills the
C     vacancy in the IS0 shell, leaving a hole in the IS1 shell).
C  *  Non-radiative: IS0-IS1-IS2 (an electron from the IS1 shell fills
C     the vacancy in the IS0 shell, and the released energy is taken
C     away by an electron in the IS2 shell; this process leaves two
C     holes, in the IS1 and IS2 shells).
C  The de-excitation cascade (i.e. the set of transitions that occur for
C  a given initial vacancy) is sampled from the transition probabilities
C  contained in the Livermore Evaluated Atomic Data Library (EADL). The
C  energy of the radiation emitted in each transition is read from the
C  PENELOPE database.
C
C  The simulation of the de-excitation cascade is discontinued either
C  when the K, L and M shells have been filled up or when there is not
C  enough energy to produce 'active' radiation (with energy larger than
C  EABS). The excitation energy of the residual ion is assumed to be
C  deposited locally. We disregard the emission and transport of soft
C  x-rays and slow electrons, whose energies are less than the binding
C  energy of the N1 shell of the heavier element in the medium. This
C  sets a lower limit for the energy interval that can be covered by the
C  simulation program in a consistent way.
C
C  De-excitation data for the loaded elements are stored in the common
C  block /CRELAX/, in a form designed to minimize the amount of memory
C  and to facilitate the random sampling. The quantities in the common
C  block are the following:
C  IFIRST(99,9) ... de-excitation data for a vacancy in the shell IS of
C     the element IZ start at the position K=IFIRST(IZ,IS) in the
C     storage arrays. The allowed values for IS are 1 to 9 (K shell and
C     L and M subshells).
C  ILAST(99,9) ... the de-excitation data for a vacancy in the shell IS
C     of the element IZ end at the position K=ILAST(IZ,IS) in the
C     storage arrays.
C  IS1(K), IS2(K) ... shells that are active in the transition (see the
C     shell label code below). For radiative transitions, IS2(K)=0.
C  P(K) ... relative probability for the transition IS-IS1(K)-IS2(K).
C  ET(K) ... energy of the secondary particle emitted in the transition.
C  F(K), IAL(K) ... cutoff and alias values (Walker's sampling method).
C
C  ---------------------------------------------------------------------
C  Label code IS for electron shells:
C      1 = K  (1s1/2),     11 = N2 (4p1/2),     21 = O5 (5d5/2),
C      2 = L1 (2s1/2),     12 = N3 (4p3/2),     22 = O6 (5f5/2),
C      3 = L2 (2p1/2),     13 = N4 (4d3/2),     23 = O7 (5f7/2),
C      4 = L3 (2p3/2),     14 = N5 (4d5/2),     24 = P1 (6s1/2),
C      5 = M1 (3s1/2),     15 = N6 (4f5/2),     25 = P2 (6p1/2),
C      6 = M2 (3p1/2),     16 = N7 (4f7/2),     26 = P3 (6p3/2),
C      7 = M3 (3p3/2),     17 = O1 (5s1/2),     27 = P4 (6d3/2),
C      8 = M4 (3d3/2),     18 = O2 (5p1/2),     28 = P5 (6d5/2),
C      9 = M5 (3d5/2),     19 = O3 (5p3/2),     29 = Q1 (7s1/2),
C     10 = N1 (4s1/2),     20 = O4 (5d3/2),     30 = outer shells.
C  ---------------------------------------------------------------------
C
      use mtpar
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N)
      CHARACTER*2 LASYMB
      PARAMETER (PI=3.1415926535897932D0, TWOPI=PI+PI)
      DIMENSION ISV(256)
C  ****  Element data.
      COMMON/CADATA/ATW(99),EPX(99),RA1(99),RA2(99),RA3(99),RA4(99),
     1  RA5(99),RSCR(99),ETA(99),EB(99,30),IFI(99,30),IKS(99,30),
     2  NSHT(99),LASYMB(99)
C  ****  Main-PENELOPE commons.
      PARAMETER (MAXMAT=10)
      COMMON/TRACK/E,X,Y,Z,U,V,W,WGHT,KPAR,IBODY,M,ILB(5)
      COMMON/CSIMPA/EABS(3,MAXMAT),C1(MAXMAT),C2(MAXMAT),WCC(MAXMAT),
     1  WCR(MAXMAT)
      COMMON/CECUTR/ECUTR(MAXMAT)
C  ****  Atomic relaxation data.
      PARAMETER (NRX=15000)
      COMMON/CRELAX/P(NRX),ET(NRX),F(NRX),IAL(NRX),IS1(NRX),IS2(NRX),
     1              IFIRST(99,9),ILAST(99,9),NCUR,KS,MODER
      COMMON/CHIST/ILBA(5)
C
      ! EXTERNALRAND
C
C  ****  Initialization.
C
      IF(IZ.LT.6.OR.IS.GT.9) RETURN
C  ****  If the shell ionization energy is less than ECUTR, the cascade
C        is not followed.
      IF(EB(IZ,IS).LT.ECUTR(M)) RETURN
C
      NV=1
      ISV(1)=IS
C
C  ****  Next transition.
C
    1 CONTINUE
      ISP=ISV(NV)
      KF=IFIRST(IZ,ISP)
      KL=ILAST(IZ,ISP)
      NV=NV-1
      IF(KL.GT.KF) THEN
C  ****  Walker's sampling algorithm.
        status=rmt(r)
        RN=r*(KL-KF+1)
        K1=INT(RN)
        TST=RN-K1
        IF(TST.GT.F(KF+K1)) THEN
          KS=IAL(KF+K1)
        ELSE
          KS=KF+K1
        ENDIF
      ELSE
        KS=KF
      ENDIF
C  ****  If MODER=0, control is returned to the calling program after
C  determining the first transition, KS. Useful for testing the random
C  sampling. For normal operation, we can comment out the following
C  statement.
      IF(MODER.EQ.0) RETURN
C
C  ****  Fluorescence radiation.
C
      IS1K=IS1(KS)
      IS2K=IS2(KS)
      IF(IS2K.EQ.0) THEN
        KPARS=2
        IF(IS1K.LT.10) THEN
          IF(EB(IZ,IS1K).GT.ECUTR(M)) THEN
            NV=NV+1
            ISV(NV)=IS1K
          ENDIF
        ENDIF
      ELSE
        KPARS=1
        IF(IS1K.LT.10) THEN
          IF(EB(IZ,IS1K).GT.ECUTR(M)) THEN
            NV=NV+1
            ISV(NV)=IS1K
          ENDIF
        ENDIF
        IF(IS2K.LT.10) THEN
          IF(EB(IZ,IS2K).GT.ECUTR(M)) THEN
            NV=NV+1
            ISV(NV)=IS2K
          ENDIF
        ENDIF
      ENDIF
C
C  ****  The emitted particle is stored in the secondary stack when
C        its energy ET(K) is greater than EABS.
C
      IF(ET(KS).GT.EABS(KPARS,M)) THEN
C  ****  Initial direction (isotropic).
        status=rmt(r)
        WS=-1.0D0+2.0D0*r
        SDTS=SQRT(1.0D0-WS*WS)
        status=rmt(r)
        DF=TWOPI*r
        US=COS(DF)*SDTS
        VS=SIN(DF)*SDTS
        ILBA(1)=ILB(1)+1
        ILBA(2)=KPAR
        ILBA(4)=IZ*1000000+ISP*10000+IS1K*100+IS2K
        ILBA(5)=ILB(5)
        CALL STORES(ET(KS),X,Y,Z,US,VS,WS,WGHT,KPARS,ILBA)
      ENDIF
C
C  ****  Is there any K-, L- or M-vacancy unfilled?
C
      IF(NV.GT.0) GO TO 1
      RETURN
      END
C  *********************************************************************
C                       SUBROUTINE RELAX0
C  *********************************************************************
      SUBROUTINE RELAX0
C
C  This subroutine sets all variables in common /CRELAX/ to zero.
C
      use mtpar
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N)
      CHARACTER*2 LASYMB
C  ****  Element data.
      COMMON/CADATA/ATW(99),EPX(99),RA1(99),RA2(99),RA3(99),RA4(99),
     1  RA5(99),RSCR(99),ETA(99),EB(99,30),IFI(99,30),IKS(99,30),
     2  NSHT(99),LASYMB(99)
C  ****  Atomic relaxation data.
      PARAMETER (NRX=15000)
      COMMON/CRELAX/P(NRX),ET(NRX),F(NRX),IS0(NRX),IS1(NRX),IS2(NRX),
     1              IFIRST(99,9),ILAST(99,9),NCUR,KS,MODER
C
      DO I=1,99
        NSHT(I)=0
        DO J=1,9
          IFIRST(I,J)=0
          ILAST(I,J)=0
        ENDDO
        DO J=1,30
          IFI(I,J)=0
          EB(I,J)=0.0D0
          IKS(I,J)=0
        ENDDO
      ENDDO
C
      DO I=1,NRX
        IS0(I)=0
        IS1(I)=0
        IS2(I)=0
        ET(I)=0.0D0
        P(I)=0.0D0
        F(I)=0.0D0
      ENDDO
      NCUR=0
      RETURN
      END
C  *********************************************************************
C                       SUBROUTINE RELAXR
C  *********************************************************************
      SUBROUTINE RELAXR(IRD,IWR,INFO)
C
C  This subroutine reads atomic relaxation data, for a single element,
C  from the material definition file (unit IRD) and initializes the
C  algorithm for random sampling of de-excitation cascades of this
C  element. (See heading comments in subroutine RELAX)
C
      use mtpar
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N)
      CHARACTER*2 LASYMB
      CHARACTER*5 CH5
C  ****  Element data.
      COMMON/CADATA/ATW(99),EPX(99),RA1(99),RA2(99),RA3(99),RA4(99),
     1  RA5(99),RSCR(99),ETA(99),EB(99,30),IFI(99,30),IKS(99,30),
     2  NSHT(99),LASYMB(99)
C  ****  Atomic relaxation data.
      PARAMETER (NRX=15000, NTRAN=1000)
      COMMON/CRELAX/P(NRX),ET(NRX),F(NRX),IS0(NRX),IS1(NRX),IS2(NRX),
     1              IFIRST(99,9),ILAST(99,9),NCUR,KS,MODER
      DIMENSION JS0(NTRAN),JS1(NTRAN),JS2(NTRAN),PR(NTRAN),ER(NTRAN),
     1          WW(NTRAN),FF(NTRAN),KK(NTRAN),IORD(NTRAN),ISR(NTRAN),
     2          IQQ(30),EE(30)
C
      CHARACTER*2 LSHELL(0:30)
      DATA LSHELL/'  ','K ','L1','L2','L3','M1','M2','M3','M4','M5',
     1   'N1','N2','N3','N4','N5','N6','N7','O1','O2','O3','O4',
     1   'O5','O6','O7','P1','P2','P3','P4','P5','Q1','X '/
C
      MODER=1  ! RELAX normal operation mode.
C
C  ****  Input transition data.
C
      READ(IRD,5001) IZ,NSHR,NT
 5001 FORMAT(16X,I3,18X,I3,23X,I4)
      IF(INFO.GE.2) WRITE(IWR,2001) IZ,NSHR,NT
 2001 FORMAT(/1X,'*** RELAX:  Z =',I3,',  NO. OF SHELLS =',I3,
     1          ',  NO. OF TRANSITIONS =',I4)
      IF(NT.GT.NTRAN) STOP 'RELAXR. NTRAN needs to be increased.'
      IF(NCUR+NT.GT.NRX) THEN
        WRITE(IWR,*) 'Insufficient memory storage in RELAXR.'
        WRITE(IWR,*) 'Increase the value of the parameter NRX to',
     1    NCUR+NT
        WRITE(26,*) 'Insufficient memory storage in RELAXR.'
        WRITE(26,*) 'Increase the value of the parameter NRX to',
     1    NCUR+NT
        STOP 'RELAXR. Insufficient memory storage.'
      ENDIF
C
      IF(INFO.GE.2) WRITE(IWR,2002)
 2002 FORMAT(/2X,'i',2X,' Shell    f',5X,'Ui (eV)',/1X,28('-'))
      DO IS=1,NSHR
        READ(IRD,5002) ISR(IS),CH5,IQQ(IS),EE(IS)
 5002   FORMAT(1X,I3,1X,A5,1X,I1,E16.8)
        IF(INFO.GE.2) WRITE(IWR,2003) ISR(IS),CH5,LSHELL(ISR(IS)),
     1    IQQ(IS),EE(IS)
 2003   FORMAT(1X,I3,1X,A5,1X,A2,2X,I1,1P,E16.8)
      ENDDO
C
      IF(NT.GT.0) THEN
        IF(INFO.GE.2) WRITE(IWR,2004)
 2004   FORMAT(/2X,'S0 S1 S2',3X,'Probability',5X,'Energy (eV)',
     1    /1X,42('-'))
        DO I=1,NT
          READ(IRD,*) JS0(I),JS1(I),JS2(I),PR(I),ER(I)
          IF(INFO.GE.2) WRITE(IWR,2005) LSHELL(JS0(I)),LSHELL(JS1(I)),
     1      LSHELL(JS2(I)),PR(I),ER(I)
 2005     FORMAT(1X,3(1X,A2),1P,2E16.8)
          IF(PR(I).LT.1.0D-35) THEN
            WRITE(IWR,*) 'Negative transition probability?'
            STOP 'RELAXR. Negative transition probability?'
          ENDIF
          IF(JS1(I).GT.16) JS1(I)=30
          IF(JS2(I).GT.16) JS2(I)=30
        ENDDO
      ENDIF
C
C  ****  Check if this element's data have already been loaded.
C
      IF(IFIRST(IZ,1).NE.0) RETURN
C
      NSHT(IZ)=NSHR
      DO IS=1,NSHR
        IFI(IZ,ISR(IS))=IQQ(IS)
        EB(IZ,ISR(IS))=EE(IS)
      ENDDO
      IF(NT.EQ.0) THEN
        DO IS=1,9
          IFIRST(IZ,IS)=NCUR+1
          ILAST(IZ,IS)=NCUR+1
C  ****  The array IS0 contains the alias values.
          IS0(NCUR+1)=NCUR+1
          P(NCUR+1)=1.0D0
          F(NCUR+1)=1.0D0
          ET(NCUR+1)=0.0D0
          IS1(NCUR+1)=1
          IS2(NCUR+1)=1
          NCUR=NCUR+1
        ENDDO
        RETURN
      ENDIF
C
C  ****  Walker's aliasing.
C
      DO IS=1,9
        N=0
        DO J=1,NT
          IF(JS0(J).EQ.IS) THEN
            N=N+1
            IORD(N)=J
            WW(N)=PR(J)
          ENDIF
        ENDDO
        IF(N.GT.1) THEN
          CALL IRND0(WW,FF,KK,N)
          IFIRST(IZ,IS)=NCUR+1
          ILAST(IZ,IS)=NCUR+N
          DO L=1,N
            P(NCUR+L)=WW(L)
            F(NCUR+L)=FF(L)
            ET(NCUR+L)=ER(IORD(L))
C  ****  The array IS0 contains the alias values.
            IS0(NCUR+L)=IFIRST(IZ,IS)+KK(L)-1
            IS1(NCUR+L)=JS1(IORD(L))
            IS2(NCUR+L)=JS2(IORD(L))
          ENDDO
          NCUR=NCUR+N
        ELSE
          NCUR=NCUR+1
          IFIRST(IZ,IS)=NCUR
          ILAST(IZ,IS)=NCUR
          IS0(NCUR+1)=NCUR
          P(NCUR)=1.0D0
          F(NCUR)=1.0D0
          ET(NCUR)=ER(1)
          IS1(NCUR)=JS1(1)
          IS2(NCUR)=JS2(1)
        ENDIF
      ENDDO
C
C  ****  Verify that transition probabilities are correctly reproduced.
C
      TST=0.0D0
      DO IS=1,9
        I0=IFIRST(IZ,IS)
        IN=ILAST(IZ,IS)
        PT=0.0D0
        DO I=I0,IN
          PT=PT+P(I)
        ENDDO
        DO I=I0,IN
          PPI=0.0D0
          DO J=I0,IN
            IF(IS0(J).EQ.I) PPI=PPI+(1.0D0-F(J))
          ENDDO
          PPI=(PPI+F(I))/DBLE(IN-I0+1)
          TST=MAX(TST,ABS(1.0D0-PPI*PT/P(I)))
        ENDDO
      ENDDO
      IF(TST.GT.1.0D-12) STOP 'RELAXR. Rounding error is too large.'
      RETURN
      END
C  *********************************************************************
C                       SUBROUTINE RELAXW
C  *********************************************************************
      SUBROUTINE RELAXW(IZ,IWR)
C
C  This subroutine produces a table of atomic relaxation data for the
C  element IZ, and prints it on unit IWR. The output table is part of
C  PENELOPE's material definition file.
C
C  Data are read from file 'pdrelax.p05', which contains data pertaining
C  to singly ionized atoms with the initial vacancy in one of the K, L
C  and M shells. This file was prepared from the Livermore Evaluated
C  Atomic Data Library (EADL).
C
C  NOTE: The transition probabilities and emission energies can be
C  modified by editing the material data file. For each initial vacancy,
C  the sum of transition probabilities _must_ be equal to unity.
C
      use mtpar
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N)
      CHARACTER*2 LASYMB
      CHARACTER*5 CSH5(30)
C  ****  ELEMENT DATA.
      COMMON/CADATA/ATW(99),EPX(99),RA1(99),RA2(99),RA3(99),RA4(99),
     1  RA5(99),RSCR(99),ETA(99),EB(99,30),IFI(99,30),IKS(99,30),
     2  NSHT(99),LASYMB(99)
C
      PARAMETER (NM=1000)
      DIMENSION IS0(NM),IS1(NM),IS2(NM),P(NM),EI(NM),EE(99)
C
      DATA CSH5/'1s1/2','2s1/2','2p1/2','2p3/2','3s1/2','3p1/2',
     1          '3p3/2','3d3/2','3d5/2','4s1/2','4p1/2','4p3/2',
     2          '4d3/2','4d5/2','4f5/2','4f7/2','5s1/2','5p1/2',
     3          '5p3/2','5d3/2','5d5/2','5f5/2','5f7/2','6s1/2',
     4          '6p1/2','6p3/2','6d3/2','6d5/2','7s1/2',' free'/
C
      IF(NSHT(IZ).LE.0) STOP 'RELAXW. The element is not loaded.'
C
      OPEN(3,FILE='pdrelax.p05')
      READ(3,*,END=1) IZR,IS0R
      NT=0
      DO I=1,150000
        READ(3,*,END=1) IZR,IS0R,IS1R,IS2R,PR,EIN
        IF(IZR.EQ.IZ) THEN
          NT=NT+1
          IS0(NT)=IS0R
          IS1(NT)=IS1R
          IS2(NT)=IS2R
          P(NT)=PR
          EI(NT)=EIN
        ENDIF
      ENDDO
    1 CONTINUE
      CLOSE(3)
C
      WRITE(IWR,2001) IZ,NSHT(IZ),NT
 2001 FORMAT(1X,'*** RELAX:  Z =',I3,',  NO. OF SHELLS =',I3,
     1          ',  NO. OF TRANSITIONS =',I4)
C
      DO I=1,99
        EE(I)=0.0D0
      ENDDO
      DO J=1,30
        KS=IKS(IZ,J)
        IF(KS.GT.0) THEN
          IF(IFI(IZ,KS).NE.0) THEN
            WRITE(IWR,2002) KS,CSH5(KS),IFI(IZ,KS),EB(IZ,KS)
 2002       FORMAT(1X,I3,1X,A5,1X,I1,1P,E16.8)
            EE(KS)=EB(IZ,KS)
          ENDIF
        ENDIF
      ENDDO
C
      IF(NT.GT.0) THEN
        DO I=1,NT
          IF(IS2(I).EQ.0) THEN
            IF(EI(I).LT.1.0D0) THEN
              ET=EE(IS0(I))-EE(IS1(I))
            ELSE
              ET=EI(I)
            ENDIF
          ELSE
            IF(EI(I).LT.1.0D0) THEN
              ET=EE(IS0(I))-EE(IS1(I))-EE(IS2(I))
            ELSE
              ET=EI(I)
            ENDIF
          ENDIF
          ET=MAX(ET,1.0D0)
          WRITE(IWR,2003) IS0(I),IS1(I),IS2(I),P(I),ET
        ENDDO
 2003   FORMAT(1X,3I3,1P,2E16.8)
      ENDIF
      RETURN
      END
C  *********************************************************************
C                       BLOCK DATA PENDAT
C  *********************************************************************
      BLOCK DATA PENDAT
C
C  Physical data for the elements Z=1-99.
C
      use mtpar
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N)
      CHARACTER*2 LASYMB
C
      COMMON/CADATA/ATW(99),EPX(99),RA1(99),RA2(99),RA3(99),RA4(99),
     1  RA5(99),RSCR(99),ETA(99),EB(99,30),IFI(99,30),IKS(99,30),
     2  NSHT(99),LASYMB(99)
C
C  ************  Chemical symbols of the elements.
C
      DATA LASYMB/'H ','He','Li','Be','B ','C ','N ','O ','F ',
     1     'Ne','Na','Mg','Al','Si','P ','S ','Cl','Ar','K ',
     2     'Ca','Sc','Ti','V ','Cr','Mn','Fe','Co','Ni','Cu',
     3     'Zn','Ga','Ge','As','Se','Br','Kr','Rb','Sr','Y ',
     4     'Zr','Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd','In',
     5     'Sn','Sb','Te','I ','Xe','Cs','Ba','La','Ce','Pr',
     6     'Nd','Pm','Sm','Eu','Gd','Tb','Dy','Ho','Er','Tm',
     7     'Yb','Lu','Hf','Ta','W ','Re','Os','Ir','Pt','Au',
     8     'Hg','Tl','Pb','Bi','Po','At','Rn','Fr','Ra','Ac',
     9     'Th','Pa','U ','Np','Pu','Am','Cm','Bk','Cf','Es'/
C
C  ************  Atomic weights (g/mol).
C
      DATA ATW  /1.0079D0,4.0026D0,6.9410D0,9.0122D0,1.0811D1,
     1  1.2011D1,1.4007D1,1.5999D1,1.8998D1,2.0179D1,2.2990D1,
     2  2.4305D1,2.6982D1,2.8086D1,3.0974D1,3.2066D1,3.5453D1,
     3  3.9948D1,3.9098D1,4.0078D1,4.4956D1,4.7880D1,5.0942D1,
     4  5.1996D1,5.4938D1,5.5847D1,5.8933D1,5.8690D1,6.3546D1,
     5  6.5390D1,6.9723D1,7.2610D1,7.4922D1,7.8960D1,7.9904D1,
     6  8.3800D1,8.5468D1,8.7620D1,8.8906D1,9.1224D1,9.2906D1,
     7  9.5940D1,9.7907D1,1.0107D2,1.0291D2,1.0642D2,1.0787D2,
     8  1.1241D2,1.1482D2,1.1871D2,1.2175D2,1.2760D2,1.2690D2,
     9  1.3129D2,1.3291D2,1.3733D2,1.3891D2,1.4012D2,1.4091D2,
     A  1.4424D2,1.4491D2,1.5036D2,1.5196D2,1.5725D2,1.5893D2,
     B  1.6250D2,1.6493D2,1.6726D2,1.6893D2,1.7304D2,1.7497D2,
     C  1.7849D2,1.8095D2,1.8385D2,1.8621D2,1.9020D2,1.9222D2,
     D  1.9508D2,1.9697D2,2.0059D2,2.0438D2,2.0720D2,2.0898D2,
     E  2.0898D2,2.0999D2,2.2202D2,2.2302D2,2.2603D2,2.2703D2,
     F  2.3204D2,2.3104D2,2.3803D2,2.3705D2,2.3905D2,2.4306D2,
     G  2.4707D2,2.4707D2,2.5108D2,2.5208D2/
C
C  ************  Mean excitation energies of the elements (eV).
C
      DATA EPX / 19.2D0, 41.8D0, 40.0D0, 63.7D0, 76.0D0, 81.0D0,
     1   82.0D0, 95.0D0,115.0D0,137.0D0,149.0D0,156.0D0,166.0D0,
     2  173.0D0,173.0D0,180.0D0,174.0D0,188.0D0,190.0D0,191.0D0,
     3  216.0D0,233.0D0,245.0D0,257.0D0,272.0D0,286.0D0,297.0D0,
     4  311.0D0,322.0D0,330.0D0,334.0D0,350.0D0,347.0D0,348.0D0,
     5  343.0D0,352.0D0,363.0D0,366.0D0,379.0D0,393.0D0,417.0D0,
     6  424.0D0,428.0D0,441.0D0,449.0D0,470.0D0,470.0D0,469.0D0,
     7  488.0D0,488.0D0,487.0D0,485.0D0,491.0D0,482.0D0,488.0D0,
     8  491.0D0,501.0D0,523.0D0,535.0D0,546.0D0,560.0D0,574.0D0,
     9  580.0D0,591.0D0,614.0D0,628.0D0,650.0D0,658.0D0,674.0D0,
     A  684.0D0,694.0D0,705.0D0,718.0D0,727.0D0,736.0D0,746.0D0,
     B  757.0D0,790.0D0,790.0D0,800.0D0,810.0D0,823.0D0,823.0D0,
     C  830.0D0,825.0D0,794.0D0,827.0D0,826.0D0,841.0D0,847.0D0,
     D  878.0D0,890.0D0,902.0D0,921.0D0,934.0D0,939.0D0,952.0D0,
     E  966.0D0,980.0D0/
C
C  ************  Atomic form factor parameters.
C
      DATA RA1/0.0D0, 3.9265D+0, 4.3100D+1, 5.2757D+1, 2.5021D+1,
     1     1.2211D+1, 9.3229D+0, 3.2455D+0, 2.4197D+0, 1.5985D+0,
     2     3.0926D+1, 1.5315D+1, 7.7061D+0, 3.9493D+0, 2.2042D+0,
     3     1.9453D+1, 1.9354D+1, 8.0374D+0, 8.3779D+1, 5.7370D+1,
     4     5.2310D+1, 4.7514D+1, 4.3785D+1, 4.2048D+1, 3.6972D+1,
     5     3.3849D+1, 3.1609D+1, 2.8763D+1, 2.7217D+1, 2.4263D+1,
     6     2.2403D+1, 1.8606D+1, 1.5143D+1, 1.4226D+1, 1.1792D+1,
     7     9.7574D+0, 1.2796D+1, 1.2854D+1, 1.2368D+1, 1.0208D+1,
     8     8.2823D+0, 7.4677D+0, 7.6028D+0, 6.1090D+0, 5.5346D+0,
     9     4.2340D+0, 4.0444D+0, 4.2905D+0, 4.7950D+0, 5.1112D+0,
     A     5.2407D+0, 5.2153D+0, 5.1639D+0, 4.8814D+0, 5.8054D+0,
     B     6.6724D+0, 6.5104D+0, 6.3364D+0, 6.2889D+0, 6.3028D+0,
     C     6.3853D+0, 6.3475D+0, 6.5779D+0, 6.8486D+0, 7.0993D+0,
     D     7.6122D+0, 7.9681D+0, 8.3481D+0, 6.3875D+0, 8.0042D+0,
     E     8.0820D+0, 7.6940D+0, 7.1927D+0, 6.6751D+0, 6.1623D+0,
     F     5.8335D+0, 5.5599D+0, 4.6551D+0, 4.4327D+0, 4.7601D+0,
     G     5.2872D+0, 5.6084D+0, 5.7680D+0, 5.8041D+0, 5.7566D+0,
     H     5.6541D+0, 6.3932D+0, 6.9313D+0, 7.0027D+0, 6.8796D+0,
     I     6.4739D+0, 6.2405D+0, 6.0081D+0, 5.5708D+0, 5.3680D+0,
     J     5.8660D+0, 5.6375D+0, 5.1719D+0, 4.8989D+0/
      DATA RA2/0.0D0, 1.3426D-1, 9.4875D+1,-1.0896D+2,-4.5494D+1,
     1    -1.9572D+1,-1.2382D+1,-3.6827D+0,-2.4542D+0,-1.4453D+0,
     2     1.3401D+2, 7.9717D+1, 6.2164D+1, 4.0300D+1, 3.1682D+1,
     3    -1.3639D+1,-1.5950D+1,-5.1523D+0, 1.8351D+2, 1.2205D+2,
     4     1.0007D+2, 8.5632D+1, 7.9145D+1, 6.3675D+1, 6.2954D+1,
     5     5.6601D+1, 5.4171D+1, 4.8752D+1, 3.8062D+1, 3.9933D+1,
     6     4.8343D+1, 4.2137D+1, 3.4617D+1, 2.9430D+1, 2.4010D+1,
     7     1.9744D+1, 4.0009D+1, 5.1614D+1, 5.0456D+1, 3.9088D+1,
     8     2.6824D+1, 2.2953D+1, 2.4773D+1, 1.6893D+1, 1.4548D+1,
     9     9.7226D+0, 1.0192D+1, 1.1153D+1, 1.3188D+1, 1.4733D+1,
     A     1.5644D+1, 1.5939D+1, 1.5923D+1, 1.5254D+1, 2.0748D+1,
     B     2.6901D+1, 2.7032D+1, 2.4938D+1, 2.1528D+1, 2.0362D+1,
     C     1.9474D+1, 1.8238D+1, 1.7898D+1, 1.9174D+1, 1.9023D+1,
     D     1.8194D+1, 1.8504D+1, 1.8955D+1, 1.4276D+1, 1.7558D+1,
     E     1.8651D+1, 1.7984D+1, 1.6793D+1, 1.5469D+1, 1.4143D+1,
     F     1.3149D+1, 1.2255D+1, 9.2352D+0, 8.6067D+0, 9.7460D+0,
     G     1.1749D+1, 1.3281D+1, 1.4326D+1, 1.4920D+1, 1.5157D+1,
     H     1.5131D+1, 1.9489D+1, 2.3649D+1, 2.4686D+1, 2.4760D+1,
     I     2.1519D+1, 2.0099D+1, 1.8746D+1, 1.5943D+1, 1.4880D+1,
     J     1.6345D+1, 1.5283D+1, 1.3061D+1, 1.1553D+1/
      DATA RA3/0.0D0, 2.2648D+0, 1.0579D+3, 8.6177D+2, 2.4422D+2,
     1     7.8788D+1, 3.8293D+1, 1.2564D+1, 6.9091D+0, 3.7926D+0,
     2     0.0000D+0, 0.0000D+0, 1.6759D-9, 1.3026D+1, 3.0569D+0,
     3     1.5521D+2, 1.2815D+2, 4.7378D+1, 9.2802D+2, 4.7508D+2,
     4     3.6612D+2, 2.7582D+2, 2.1008D+2, 1.5903D+2, 1.2322D+2,
     5     9.2898D+1, 7.1345D+1, 5.1651D+1, 3.8474D+1, 2.7410D+1,
     6     1.9126D+1, 1.0889D+1, 5.3479D+0, 8.2223D+0, 5.0837D+0,
     7     2.8905D+0, 2.7457D+0, 6.7082D-1, 0.0000D+0, 0.0000D+0,
     8     0.0000D+0, 0.0000D+0, 0.0000D+0, 0.0000D+0, 0.0000D+0,
     9     0.0000D+0, 0.0000D+0, 0.0000D+0, 0.0000D+0, 0.0000D+0,
     A     0.0000D+0, 0.0000D+0, 0.0000D+0, 0.0000D+0, 0.0000D+0,
     B     0.0000D+0, 0.0000D+0, 0.0000D+0, 1.7264D-1, 2.7322D-1,
     C     3.9444D-1, 4.5648D-1, 6.2286D-1, 7.2468D-1, 8.4296D-1,
     D     1.1698D+0, 1.2994D+0, 1.4295D+0, 0.0000D+0, 8.1570D-1,
     E     6.9349D-1, 4.9536D-1, 3.1211D-1, 1.5931D-1, 2.9512D-2,
     F     0.0000D+0, 0.0000D+0, 0.0000D+0, 0.0000D+0, 0.0000D+0,
     G     0.0000D+0, 0.0000D+0, 0.0000D+0, 0.0000D+0, 0.0000D+0,
     H     0.0000D+0, 0.0000D+0, 0.0000D+0, 0.0000D+0, 0.0000D+0,
     I     0.0000D+0, 0.0000D+0, 0.0000D+0, 0.0000D+0, 0.0000D+0,
     J     0.0000D+0, 0.0000D+0, 0.0000D+0, 0.0000D+0/
      DATA RA4/1.1055D1,6.3519D0,4.7367D+1, 3.9402D+1, 2.2896D+1,
     1     1.3979D+1, 1.0766D+1, 6.5252D+0, 5.1631D+0, 4.0524D+0,
     2     2.7145D+1, 1.8724D+1, 1.4782D+1, 1.1608D+1, 9.7750D+0,
     3     1.6170D+1, 1.5249D+1, 9.1916D+0, 5.4499D+1, 4.1381D+1,
     4     3.7395D+1, 3.3815D+1, 3.1135D+1, 2.8273D+1, 2.6140D+1,
     5     2.3948D+1, 2.2406D+1, 2.0484D+1, 1.8453D+1, 1.7386D+1,
     6     1.7301D+1, 1.5388D+1, 1.3411D+1, 1.2668D+1, 1.1133D+1,
     7     9.8081D+0, 1.3031D+1, 1.4143D+1, 1.3815D+1, 1.2077D+1,
     8     1.0033D+1, 9.2549D+0, 9.5338D+0, 7.9076D+0, 7.3263D+0,
     9     5.9996D+0, 6.0087D+0, 6.2660D+0, 6.7914D+0, 7.1501D+0,
     A     7.3367D+0, 7.3729D+0, 7.3508D+0, 7.1465D+0, 8.2731D+0,
     B     9.3745D+0, 9.3508D+0, 8.9897D+0, 8.4566D+0, 8.2690D+0,
     C     8.1398D+0, 7.9183D+0, 7.9123D+0, 8.1677D+0, 8.1871D+0,
     D     8.1766D+0, 8.2881D+0, 8.4227D+0, 7.0273D+0, 8.0002D+0,
     E     8.1440D+0, 7.9104D+0, 7.5685D+0, 7.1970D+0, 6.8184D+0,
     F     6.5469D+0, 6.3056D+0, 5.4844D+0, 5.2832D+0, 5.5889D+0,
     G     6.0919D+0, 6.4340D+0, 6.6426D+0, 6.7428D+0, 6.7636D+0,
     H     6.7281D+0, 7.5729D+0, 8.2808D+0, 8.4400D+0, 8.4220D+0,
     I     7.8662D+0, 7.5993D+0, 7.3353D+0, 6.7829D+0, 6.5520D+0,
     J     6.9181D+0, 6.6794D+0, 6.1735D+0, 5.8332D+0/
      DATA RA5/0.0D0, 4.9828D+0, 5.5674D+1, 3.0902D+1, 1.1496D+1,
     1     4.8936D+0, 2.5506D+0, 1.2236D+0, 7.4698D-1, 4.7042D-1,
     2     4.7809D+0, 4.6315D+0, 4.3677D+0, 4.9269D+0, 2.6033D+0,
     3     9.6229D+0, 7.2592D+0, 4.1634D+0, 1.3999D+1, 8.6975D+0,
     4     6.9630D+0, 5.4681D+0, 4.2653D+0, 3.2848D+0, 2.7354D+0,
     5     2.1617D+0, 1.7030D+0, 1.2826D+0, 9.7080D-1, 7.2227D-1,
     6     5.0874D-1, 3.1402D-1, 1.6360D-1, 3.2918D-1, 2.3570D-1,
     7     1.5868D-1, 1.5146D-1, 9.7662D-2, 7.3151D-2, 6.4206D-2,
     8     4.8945D-2, 4.3189D-2, 4.4368D-2, 3.3976D-2, 3.0466D-2,
     9     2.4477D-2, 3.7202D-2, 3.7093D-2, 3.8161D-2, 3.8576D-2,
     A     3.8403D-2, 3.7806D-2, 3.4958D-2, 3.6029D-2, 4.3087D-2,
     B     4.7069D-2, 4.6452D-2, 4.2486D-2, 4.1517D-2, 4.1691D-2,
     C     4.2813D-2, 4.2294D-2, 4.5287D-2, 4.8462D-2, 4.9726D-2,
     D     5.5097D-2, 5.6568D-2, 5.8069D-2, 1.2270D-2, 3.8006D-2,
     E     3.5048D-2, 3.0050D-2, 2.5069D-2, 2.0485D-2, 1.6151D-2,
     F     1.4631D-2, 1.4034D-2, 1.1978D-2, 1.1522D-2, 1.2375D-2,
     G     1.3805D-2, 1.4954D-2, 1.5832D-2, 1.6467D-2, 1.6896D-2,
     H     1.7166D-2, 1.9954D-2, 2.2497D-2, 2.1942D-2, 2.1965D-2,
     I     2.0005D-2, 1.8927D-2, 1.8167D-2, 1.6314D-2, 1.5522D-2,
     J     1.3141D-2, 1.2578D-2, 1.2238D-2, 1.0427D-2/
C
C  ************  Pair production cross section parameters.
C
C  ****  Screening parameter (R mc/hbar).
      DATA RSCR  /1.2281D2,7.3167D1,6.9228D1,6.7301D1,6.4696D1,
     1   6.1228D1,5.7524D1,5.4033D1,5.0787D1,4.7851D1,4.6373D1,
     2   4.5401D1,4.4503D1,4.3815D1,4.3074D1,4.2321D1,4.1586D1,
     3   4.0953D1,4.0524D1,4.0256D1,3.9756D1,3.9144D1,3.8462D1,
     4   3.7778D1,3.7174D1,3.6663D1,3.5986D1,3.5317D1,3.4688D1,
     5   3.4197D1,3.3786D1,3.3422D1,3.3068D1,3.2740D1,3.2438D1,
     6   3.2143D1,3.1884D1,3.1622D1,3.1438D1,3.1142D1,3.0950D1,
     7   3.0758D1,3.0561D1,3.0285D1,3.0097D1,2.9832D1,2.9581D1,
     8   2.9411D1,2.9247D1,2.9085D1,2.8930D1,2.8721D1,2.8580D1,
     9   2.8442D1,2.8312D1,2.8139D1,2.7973D1,2.7819D1,2.7675D1,
     A   2.7496D1,2.7285D1,2.7093D1,2.6911D1,2.6705D1,2.6516D1,
     B   2.6304D1,2.6108D1,2.5929D1,2.5730D1,2.5577D1,2.5403D1,
     C   2.5245D1,2.5100D1,2.4941D1,2.4790D1,2.4655D1,2.4506D1,
     D   2.4391D1,2.4262D1,2.4145D1,2.4039D1,2.3922D1,2.3813D1,
     E   2.3712D1,2.3621D1,2.3523D1,2.3430D1,2.3331D1,2.3238D1,
     F   2.3139D1,2.3048D1,2.2967D1,2.2833D1,2.2694D1,2.2624D1,
     G   2.2545D1,2.2446D1,2.2358D1,2.2264D1/
C  ****  Asymptotic triplet contribution (eta).
      DATA ETA   /1.1570D0,1.1690D0,1.2190D0,1.2010D0,1.1890D0,
     1   1.1740D0,1.1760D0,1.1690D0,1.1630D0,1.1570D0,1.1740D0,
     2   1.1830D0,1.1860D0,1.1840D0,1.1800D0,1.1780D0,1.1750D0,
     3   1.1700D0,1.1800D0,1.1870D0,1.1840D0,1.1800D0,1.1770D0,
     4   1.1660D0,1.1690D0,1.1660D0,1.1640D0,1.1620D0,1.1540D0,
     5   1.1560D0,1.1570D0,1.1580D0,1.1570D0,1.1580D0,1.1580D0,
     6   1.1580D0,1.1660D0,1.1730D0,1.1740D0,1.1750D0,1.1700D0,
     7   1.1690D0,1.1720D0,1.1690D0,1.1680D0,1.1640D0,1.1670D0,
     8   1.1700D0,1.1720D0,1.1740D0,1.1750D0,1.1780D0,1.1790D0,
     9   1.1800D0,1.1870D0,1.1940D0,1.1970D0,1.1960D0,1.1940D0,
     A   1.1940D0,1.1940D0,1.1940D0,1.1940D0,1.1960D0,1.1970D0,
     B   1.1960D0,1.1970D0,1.1970D0,1.1980D0,1.1980D0,1.2000D0,
     C   1.2010D0,1.2020D0,1.2040D0,1.2050D0,1.2060D0,1.2080D0,
     D   1.2070D0,1.2080D0,1.2120D0,1.2150D0,1.2180D0,1.2210D0,
     E   1.2240D0,1.2270D0,1.2300D0,1.2370D0,1.2430D0,1.2470D0,
     F   1.2500D0,1.2510D0,1.2520D0,1.2550D0,1.2560D0,1.2570D0,
     G   1.2590D0,1.2620D0,1.2620D0,1.2650D0/
C
      END

CXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
C                                                                      C
C  NUMERICAL TOOLS     (Francesc Salvat. Barcelona. February, 2001.)   C
C                                                                      C
CXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

C  *********************************************************************
C                       SUBROUTINE MERGE2
C  *********************************************************************
      SUBROUTINE MERGE2(X1,Y1,X2,Y2,X,Y12,N1,N2,N)
C
C  This subroutine merges two tables (X1,Y1), (X2,Y2) of two functions
C  to produce a table (X,Y12) of the sum of these functions, with abs-
C  cissas in increasing order. N1, N2 and N are the numbers of grid
C  points in the input and merged tables. A discontinuity in the func-
C  tion is described by giving twice the abscissa. Log-log linear int-
C  erpolation is used to interpolate the input tables.
C
      use mtpar
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N)
      PARAMETER (EPS=1.0D-10,NP=1500,NP2=NP+NP)
      DIMENSION X1(NP),Y1(NP),X2(NP),Y2(NP),X(NP2),Y12(NP2)
C
      CALL SORT2(X1,Y1,N1)
      CALL SORT2(X2,Y2,N2)
C
      DO I1=1,N1
        X(I1)=X1(I1)
      ENDDO
      N=N1
      DO I2=1,N2
        XC=X2(I2)
        CALL FINDI(X1,XC,N1,I1)
        IF(I1.EQ.N1) I1=N1-1
        TST1=ABS(XC-X1(I1))
        TST2=ABS(XC-X1(I1+1))
        TST12=MIN(TST1,TST2)
        IF(I2.GT.1) THEN
          TST3=ABS(XC-X2(I2-1))
        ELSE
          TST3=1.0D0
        ENDIF
        IF(I2.LT.N2) THEN
          TST4=ABS(XC-X2(I2+1))
        ELSE
          TST4=1.0D0
        ENDIF
        TST34=MIN(TST3,TST4)
        TST=EPS*XC
        IF(TST34.GT.TST) THEN
          IF(TST12.GT.TST) THEN
            N=N+1
            X(N)=XC
          ENDIF
        ELSE
          N=N+1
          X(N)=XC
        ENDIF
      ENDDO
C
C  ****  Sort and clean the merged grid.
C
    1 CONTINUE
      DO I=1,N-1
        IMIN=I
        XMIN=X(I)
        DO J=I+1,N
          IF(X(J).LT.XMIN) THEN
            IMIN=J
            XMIN=X(J)
          ENDIF
        ENDDO
        SAVE=X(I)
        X(I)=X(IMIN)
        X(IMIN)=SAVE
      ENDDO
C
      DO I=1,N-2
        IF(X(I).GT.X(I+2)*(1.0D0-EPS)) THEN
          X(I+1)=X(N)
          N=N-1
          GO TO 1
        ENDIF
      ENDDO
C
      DO I=1,N
        XC=X(I)
        IF(I.LT.N) THEN
          IF(X(I).GT.X(I+1)*(1.0D0-EPS)) XC=X(I)*(1.0D0-EPS)
        ENDIF
        IF(I.GT.1) THEN
          IF(X(I).LT.X(I-1)*(1.0D0+EPS)) XC=X(I)*(1.0D0+EPS)
        ENDIF
        CALL FINDI(X1,XC,N1,J)
        IF(J.EQ.N1) J=N1-1
        IF(X1(J+1).GT.X1(J)+EPS) THEN
          YI1=EXP(LOG(Y1(J))+LOG(XC/X1(J))*LOG(Y1(J+1)/Y1(J))
     1     /LOG(X1(J+1)/X1(J)))
        ELSE
          YI1=Y1(J)
        ENDIF
        CALL FINDI(X2,XC,N2,J)
        IF(J.EQ.N2) J=N2-1
        IF(X2(J+1).GT.X2(J)+EPS) THEN
          YI2=EXP(LOG(Y2(J))+LOG(XC/X2(J))*LOG(Y2(J+1)/Y2(J))
     1     /LOG(X2(J+1)/X2(J)))
        ELSE
          YI2=Y2(J)
        ENDIF
        Y12(I)=YI1+YI2
        IF(Y12(I).LT.1.0D-35) Y12(I)=1.0D-35
      ENDDO
C
      RETURN
      END
C  *********************************************************************
C                       SUBROUTINE SORT2
C  *********************************************************************
      SUBROUTINE SORT2(X,Y,N)
C
C  This subroutine sorts a table (X,Y) of a function with n data points.
C  A discontinuity of the function is described by giving twice the abs-
C  cissa. It is assumed that the function is strictly positive (negative
C  values of Y are set to zero).
C
      use mtpar
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N)
      PARAMETER (EPS=1.0D-10,NP=1000)
      DIMENSION X(N),Y(N),IORDER(NP)
C
      IF(N.EQ.1) RETURN
      DO I=1,N
        IORDER(I)=I
        IF(Y(I).LT.1.0D-35) Y(I)=1.0D-35
      ENDDO
C
      DO 1 I=1,N-1
        IMIN=I
        XMIN=X(I)
        DO J=I+1,N
          IF(X(J).LT.XMIN) THEN
            IMIN=J
            XMIN=X(J)
          ENDIF
        ENDDO
        SAVE=X(I)
        X(I)=X(IMIN)
        X(IMIN)=SAVE
        SAVE=Y(I)
        Y(I)=Y(IMIN)
        Y(IMIN)=SAVE
        ISAVE=IORDER(I)
        IORDER(I)=IORDER(IMIN)
        IORDER(IMIN)=ISAVE
        IF(I.EQ.1) GO TO 1
        IF(IORDER(I).LT.IORDER(I-1).AND.ABS(X(I)-X(I-1)).LT.1.0D-15)
     1    THEN
          SAVE=X(I-1)
          X(I-1)=X(I)
          X(I)=SAVE
          SAVE=Y(I-1)
          Y(I-1)=Y(I)
          Y(I)=SAVE
          ISAVE=IORDER(I-1)
          IORDER(I-1)=IORDER(I)
          IORDER(I)=ISAVE
        ENDIF
    1 CONTINUE
      I=N
      IF(IORDER(I).LT.IORDER(I-1).AND.ABS(X(I)-X(I-1)).LT.1.0D-15) THEN
        SAVE=X(I-1)
        X(I-1)=X(I)
        X(I)=SAVE
        SAVE=Y(I-1)
        Y(I-1)=Y(I)
        Y(I)=SAVE
      ENDIF
      RETURN
      END
C  *********************************************************************
C                       SUBROUTINE SPLINE
C  *********************************************************************
      SUBROUTINE SPLINE(X,Y,A,B,C,D,S1,SN,N)
C
C  Cubic spline interpolation of tabulated data.
C
C  Input:
C     X(I) (I=1:N) ... grid points (the X values must be in increasing
C                      order).
C     Y(I) (I=1:N) ... corresponding function values.
C     S1,SN .......... second derivatives at X(1) and X(N). The natural
C                      spline corresponds to taking S1=SN=0.
C     N .............. number of grid points.
C  Output:
C     A(1:N),B(1:N),C(1:N),D(1:N) ... spline coefficients.
C
C  The interpolating cubic polynomial in the I-th interval, from X(I) to
C  X(I+1), is
C               P(x) = A(I)+x*(B(I)+x*(C(I)+x*D(I)))
C
C  Reference: M.J. Maron, 'Numerical Analysis: a Practical Approach',
C             MacMillan Publ. Co., New York, 1982.
C
      use mtpar
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N)
      DIMENSION X(N),Y(N),A(N),B(N),C(N),D(N)
C
      IF(N.LT.4) THEN
        WRITE(26,10) N
   10   FORMAT(5X,'Spline interpolation cannot be performed with',
     1    I4,' points. Stop.')
        STOP 'SPLINE. N is less than 4.'
      ENDIF
      N1=N-1
      N2=N-2
C  ****  Auxiliary arrays H(=A) and DELTA(=D).
      DO I=1,N1
        IF(X(I+1)-X(I).LT.1.0D-13) THEN
          WRITE(26,11)
   11     FORMAT(5X,'Spline X values not in increasing order. Stop.')
          STOP 'SPLINE. X values not in increasing order.'
        ENDIF
        A(I)=X(I+1)-X(I)
        D(I)=(Y(I+1)-Y(I))/A(I)
      ENDDO
C  ****  Symmetric coefficient matrix (augmented).
      DO I=1,N2
        B(I)=2.0D0*(A(I)+A(I+1))
        K=N1-I+1
        D(K)=6.0D0*(D(K)-D(K-1))
      ENDDO
      D(2)=D(2)-A(1)*S1
      D(N1)=D(N1)-A(N1)*SN
C  ****  Gauss solution of the tridiagonal system.
      DO I=2,N2
        R=A(I)/B(I-1)
        B(I)=B(I)-R*A(I)
        D(I+1)=D(I+1)-R*D(I)
      ENDDO
C  ****  The SIGMA coefficients are stored in array D.
      D(N1)=D(N1)/B(N2)
      DO I=2,N2
        K=N1-I+1
        D(K)=(D(K)-A(K)*D(K+1))/B(K-1)
      ENDDO
      D(N)=SN
C  ****  Spline coefficients.
      SI1=S1
      DO I=1,N1
        SI=SI1
        SI1=D(I+1)
        H=A(I)
        HI=1.0D0/H
        A(I)=(HI/6.0D0)*(SI*X(I+1)**3-SI1*X(I)**3)
     1      +HI*(Y(I)*X(I+1)-Y(I+1)*X(I))
     2      +(H/6.0D0)*(SI1*X(I)-SI*X(I+1))
        B(I)=(HI/2.0D0)*(SI1*X(I)**2-SI*X(I+1)**2)
     1      +HI*(Y(I+1)-Y(I))+(H/6.0D0)*(SI-SI1)
        C(I)=(HI/2.0D0)*(SI*X(I+1)-SI1*X(I))
        D(I)=(HI/6.0D0)*(SI1-SI)
      ENDDO
C  ****  Natural cubic spline for X.GT.X(N).
      FN=Y(N)
      FNP=B(N1)+X(N)*(2.0D0*C(N1)+X(N)*3.0D0*D(N1))
      A(N)=FN-X(N)*FNP
      B(N)=FNP
      C(N)=0.0D0
      D(N)=0.0D0
C
      RETURN
      END
C  *********************************************************************
C                       SUBROUTINE FINDI
C  *********************************************************************
      SUBROUTINE FINDI(X,XC,N,I)
C
C  Finds the interval (X(I),X(I+1)) that contains the value XC.
C
C  Input:
C     X(I) (I=1:N) ... grid points (the X values must be in increasing
C                      order).
C     XC ............. point to be located.
C     N  ............. number of grid points.
C  Output:
C     I .............. interval index.
C
      use mtpar
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N)
      DIMENSION X(N)
C
      IF(XC.GT.X(N)) THEN
        I=N
        RETURN
      ENDIF
      IF(XC.LT.X(1)) THEN
        I=1
        RETURN
      ENDIF
      I=1
      I1=N
    1 IT=(I+I1)/2
      IF(XC.GT.X(IT)) THEN
        I=IT
      ELSE
        I1=IT
      ENDIF
      IF(I1-I.GT.1) GO TO 1
      RETURN
      END
C  *********************************************************************
C                       SUBROUTINE SINTEG
C  *********************************************************************
      SUBROUTINE SINTEG(X,A,B,C,D,XL,XU,SUM,N)
C
C  Computes the integral of a cubic spline function.
C
C  Input:
C     X(I) (I=1:N) ... grid points (the X values must be in increasing
C                      order).
C     A(1:N),B(1:N),C(1:N),D(1:N) ... spline coefficients.
C     N  ............. number of grid points.
C     XL ............. lower limit of the integral.
C     XU ............. upper limit of the integral.
C  Output:
C     SUM ............ value of the integral.
C
      use mtpar
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N)
      DIMENSION X(N),A(N),B(N),C(N),D(N)
C  ****  Set integration limits in increasing order.
      IF(XU.GT.XL) THEN
        XLL=XL
        XUU=XU
        SIGN=1.0D0
      ELSE
        XLL=XU
        XUU=XL
        SIGN=-1.0D0
      ENDIF
C  ****  Check integral limits.
      IF(XLL.LT.X(1).OR.XUU.GT.X(N)) THEN
        WRITE(26,10)
   10   FORMAT(5X,'Integral limits out of range. Stop.')
        STOP 'SINTEG. Integral limits out of range.'
      ENDIF
C  ****  Find involved intervals.
      SUM=0.0D0
      CALL FINDI(X,XLL,N,IL)
      CALL FINDI(X,XUU,N,IU)
C
      IF(IL.EQ.IU) THEN
C  ****  Only a single interval involved.
        X1=XLL
        X2=XUU
        SUM=X2*(A(IL)+X2*((B(IL)/2)+X2*((C(IL)/3)+X2*D(IL)/4)))
     1     -X1*(A(IL)+X1*((B(IL)/2)+X1*((C(IL)/3)+X1*D(IL)/4)))
      ELSE
C  ****  Contributions from several intervals.
        X1=XLL
        X2=X(IL+1)
        SUM=X2*(A(IL)+X2*((B(IL)/2)+X2*((C(IL)/3)+X2*D(IL)/4)))
     1     -X1*(A(IL)+X1*((B(IL)/2)+X1*((C(IL)/3)+X1*D(IL)/4)))
        IL=IL+1
        DO I=IL,IU
          X1=X(I)
          X2=X(I+1)
          IF(I.EQ.IU) X2=XUU
          SUMP=X2*(A(I)+X2*((B(I)/2)+X2*((C(I)/3)+X2*D(I)/4)))
     1        -X1*(A(I)+X1*((B(I)/2)+X1*((C(I)/3)+X1*D(I)/4)))
          SUM=SUM+SUMP
        ENDDO
      ENDIF
      SUM=SIGN*SUM
      RETURN
      END
C  *********************************************************************
C                       FUNCTION SUMGA
C  *********************************************************************
      FUNCTION SUMGA(FCT,XL,XU,TOL)
C
C  This function calculates the value SUMGA of the integral of the
C  (external) function FCT over the interval (XL,XU) using the 20-point
C  Gauss quadrature method with an adaptive bipartition scheme.
C
C  TOL is the tolerance, i.e. maximum allowed relative error; it should
C  not exceed 1.0D-13. A warning message is written in unit 6 when the
C  required accuracy is not attained.
C
      use mtpar
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N)
      PARAMETER (NP=10, NST=256, NCALLS=20000)
      DIMENSION X(NP),W(NP),S(NST),SN(NST),XR(NST),XRN(NST)
C  ****  Gauss 20-point integration formula.
C  Abscissas.
      DATA X/7.6526521133497334D-02,2.2778585114164508D-01,
     1       3.7370608871541956D-01,5.1086700195082710D-01,
     2       6.3605368072651503D-01,7.4633190646015079D-01,
     3       8.3911697182221882D-01,9.1223442825132591D-01,
     4       9.6397192727791379D-01,9.9312859918509492D-01/
C  Weights.
      DATA W/1.5275338713072585D-01,1.4917298647260375D-01,
     1       1.4209610931838205D-01,1.3168863844917663D-01,
     2       1.1819453196151842D-01,1.0193011981724044D-01,
     3       8.3276741576704749D-02,6.2672048334109064D-02,
     4       4.0601429800386941D-02,1.7614007139152118D-02/
C  ****  Error control.
      CTOL=MIN(MAX(TOL,1.0D-13),1.0D-2)
      PTOL=0.01D0*CTOL
      ERR=1.0D35
C  ****  Gauss integration from XL to XU.
      H=XU-XL
      SUMGA=0.0D0
      A=0.5D0*(XU-XL)
      B=0.5D0*(XL+XU)
      C=A*X(1)
      D=W(1)*(FCT(B+C)+FCT(B-C))
      DO I1=2,NP
        C=A*X(I1)
        D=D+W(I1)*(FCT(B+C)+FCT(B-C))
      ENDDO
      ICALL=NP+NP
      LH=1
      S(1)=D*A
      XR(1)=XL
C  ****  Adaptive bipartition scheme.
    1 CONTINUE
      HO=H
      H=0.5D0*H
      SUMR=0.0D0
      LHN=0
      DO I=1,LH
        SI=S(I)
        XA=XR(I)
        XB=XA+H
        XC=XA+HO
        A=0.5D0*(XB-XA)
        B=0.5D0*(XB+XA)
        C=A*X(1)
        D=W(1)*(FCT(B+C)+FCT(B-C))
        DO I2=2,NP
          C=A*X(I2)
          D=D+W(I2)*(FCT(B+C)+FCT(B-C))
        ENDDO
        S1=D*A
        A=0.5D0*(XC-XB)
        B=0.5D0*(XC+XB)
        C=A*X(1)
        D=W(1)*(FCT(B+C)+FCT(B-C))
        DO I3=2,NP
          C=A*X(I3)
          D=D+W(I3)*(FCT(B+C)+FCT(B-C))
        ENDDO
        S2=D*A
        ICALL=ICALL+4*NP
        S12=S1+S2
        IF(ABS(S12-SI).LE.MAX(PTOL*ABS(S12),1.0D-35)) THEN
          SUMGA=SUMGA+S12
        ELSE
          SUMR=SUMR+S12
          LHN=LHN+2
          IF(LHN.GE.NST) GO TO 2
          SN(LHN)=S2
          XRN(LHN)=XB
          SN(LHN-1)=S1
          XRN(LHN-1)=XA
        ENDIF
        IF(ICALL.GT.NCALLS) GO TO 2
      ENDDO
      ERR=ABS(SUMR)/MAX(ABS(SUMR+SUMGA),1.0D-35)
      IF(ERR.LT.CTOL.OR.LHN.EQ.0) RETURN
      LH=LHN
      DO I=1,LH
        S(I)=SN(I)
        XR(I)=XRN(I)
      ENDDO
      GO TO 1
C  ****  Warning (low accuracy) message.
    2 CONTINUE
      WRITE(26,11)
   11 FORMAT(/2X,'>>> SUMGA. Gauss adaptive-bipartition quadrature.')
      WRITE(26,12) XL,XU,TOL
   12 FORMAT(2X,'XL =',1P,E19.12,',  XU =',E19.12,',  TOL =',E8.1)
      WRITE(26,13) ICALL,SUMGA,ERR,LHN
   13 FORMAT(2X,'NCALLS = ',I5,',  SUMGA =',1P,E20.13,',  ERR =',E8.1,
     1      /2X,'Number of open subintervals = ',I3)
      WRITE(26,14)
   14 FORMAT(2X,'WARNING: the required accuracy has not been ',
     1  'attained.'/)
      RETURN
      END
C  *********************************************************************
C                       SUBROUTINE IRND0
C  *********************************************************************
      SUBROUTINE IRND0(W,F,K,N)
C
C  Initialization of Walker's aliasing algorithm for random sampling
C  from discrete probability distributions.
C
C  Input arguments:
C    N ........ number of different values of the random variable.
C    W(1:N) ... corresponding point probabilities (not necessarily
C               normalized to unity).
C  Output arguments:
C    F(1:N) ... cutoff values.
C    K(1:N) ... alias values.
C
      use mtpar
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N)
      DIMENSION W(N),F(N),K(N)
C  ****  Renormalization.
      WS=0.0D0
      DO I=1,N
        IF(W(I).LT.0.0D0) STOP 'IRND0. Negative point probability.'
        WS=WS+W(I)
      ENDDO
      WS=DBLE(N)/WS
      DO I=1,N
        K(I)=I
        F(I)=W(I)*WS
      ENDDO
      IF(N.EQ.1) RETURN
C  ****  Cutoff and alias values.
      DO I=1,N-1
        HLOW=1.0D0
        HIGH=1.0D0
        ILOW=0
        IHIGH=0
        DO J=1,N
          IF(K(J).EQ.J) THEN
            IF(F(J).LT.HLOW) THEN
              HLOW=F(J)
              ILOW=J
            ELSE IF(F(J).GT.HIGH) THEN
              HIGH=F(J)
              IHIGH=J
            ENDIF
          ENDIF
        ENDDO
        IF(ILOW.EQ.0.OR.IHIGH.EQ.0) RETURN
        K(ILOW)=IHIGH
        F(IHIGH)=HIGH+HLOW-1.0D0
      ENDDO
      RETURN
      END
C  *********************************************************************
C                       SUBROUTINE DIRECT
C  *********************************************************************
      SUBROUTINE DIRECT(CDT,DF,U,V,W)
C
C  This subroutine computes the new direction cosines of the particle
C  velocity after a collision with given polar and azimuthal scattering
C  angles.
C
C  Input:  U,V,W ... initial direction cosines,
C          CDT ..... cosine of the polar scattering angle,
C          DF ...... azimuthal scattering angle (rad).
C
C  Output: U,V,W ... new direction cosines.
C          CDT and DF remain unchanged.
C
      use mtpar
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N)
C
      SDF=SIN(DF)
      CDF=COS(DF)
C  ****  Ensure normalization.
      DXY=U*U+V*V
      DXYZ=DXY+W*W
      IF(ABS(DXYZ-1.0D0).GT.1.0D-14) THEN
        FNORM=1.0D0/SQRT(DXYZ)
        U=FNORM*U
        V=FNORM*V
        W=FNORM*W
        DXY=U*U+V*V
      ENDIF
C
      IF(DXY.GT.1.0D-28) THEN
        SDT=SQRT((1.0D0-CDT*CDT)/DXY)
        UP=U
        U=U*CDT+SDT*(UP*W*CDF-V*SDF)
        V=V*CDT+SDT*(V*W*CDF+UP*SDF)
        W=W*CDT-DXY*SDT*CDF
      ELSE
        SDT=SQRT(1.0D0-CDT*CDT)
        V=SDT*SDF
        IF(W.GT.0.0D0) THEN
          U=SDT*CDF
          W=CDT
        ELSE
          U=-SDT*CDF
          W=-CDT
        ENDIF
      ENDIF
C
      RETURN
      END
C  *********************************************************************
C                         FUNCTION RAND
C  *********************************************************************
!     FUNCTION r
C
C  This is an adapted version of subroutine RANECU written by F. James
C  (Comput. Phys. Commun. 60 (1990) 329-344), which has been modified to
C  give a single random number at each call.
C
C  The 'seeds' ISEED1 and ISEED2 must be initialized in the main program
C  and transferred through the named common block /RSEED/.
C
C  Some compilers incorporate an intrinsic random number generator with
C  the same name (but with different argument lists). To avoid conflict,
C  it is advisable to declare RAND as an external function in all sub-
C  programs that call it.
C
!     IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N)
!     PARAMETER (USCALE=1.0D0/2.147483563D9)
!     COMMON/RSEED/ISEED1,ISEED2
!
!     I1=ISEED1/53668
!     ISEED1=40014*(ISEED1-I1*53668)-I1*12211
!     IF(ISEED1.LT.0) ISEED1=ISEED1+2147483563
!
!     I2=ISEED2/52774
!     ISEED2=40692*(ISEED2-I2*52774)-I2*3791
!     IF(ISEED2.LT.0) ISEED2=ISEED2+2147483399
!
!     IZ=ISEED1-ISEED2
!     IF(IZ.LT.1) IZ=IZ+2147483562
!     RAND=IZ*USCALE
!
!     RETURN
!     END

CXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
C                                                                      X
C    TRANSPORT ROUTINES (Francesc Salvat. Barcelona. April, 2002.)     X
C                                                                      X
CXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

C  *********************************************************************
C                       SUBROUTINE CLEANS
C  *********************************************************************
      SUBROUTINE CLEANS
C
C  This subroutine initializes the secondary stack. It must be called
C  before starting the simulation of each primary track.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N)
C  ****  Secondary stack.
      PARAMETER (NMS=1000)
      COMMON/SECST/ES(NMS),XS(NMS),YS(NMS),ZS(NMS),US(NMS),
     1   VS(NMS),WS(NMS),WGHTS(NMS),KS(NMS),IBODYS(NMS),MS(NMS),
     2   ILBS(5,NMS),NSEC
C
      NSEC=0
      RETURN
      END
C  *********************************************************************
C                       SUBROUTINE START
C  *********************************************************************
      SUBROUTINE START
C
C  This subroutine forces the next event to be an artificial soft event.
C  It must be called when a new (primary or secondary) electron or posi-
C  tron track is started and when it crosses an interface.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N)
C  ****  Main-PENELOPE common.
      COMMON/TRACK/E,X,Y,Z,U,V,W,WGHT,KPAR,IBODY,M,ILB(5)
C  ****  Energy grid and interpolation constants for the current energy.
      PARAMETER (NEGP=200)
      COMMON/CEGRID/EL,EU,ET(NEGP),DLEMP(NEGP),DLEMP1,DLFC,
     1  XEL,XE,XEK,KE
C  ****  Current simulation state (used only in class II simulation).
      COMMON/CJUMP1/MODE,KSOFTE,KSOFTI,KDELTA
C
      IF(E.LT.1.00000001D0*EL.OR.E.GT.0.99999999D0*EU) THEN
        WRITE(26,1000) KPAR,E,EL,EU
 1000   FORMAT(/3X,'*** Energy out of range. KPAR = ',I2,',  E = ',
     1  1P,E12.5,' eV',/7X,'EMIN = ',E12.5,' eV,  EMAX = ',E12.5,
     2  ' eV.'/7X,'Check the values of EABS(KPAR,M) and EMAX.')
        STOP 'START. E out of range.'
      ENDIF
      MODE=0
      RETURN
      END
C  *********************************************************************
C                       SUBROUTINE JUMP
C  *********************************************************************
      SUBROUTINE JUMP(DSMAX,DS)
C
C  Calculation of the free path from the starting point to the position
C  of the next event and of the probabilities of occurrence of different
C  events.
C
C  Arguments:
C    DSMAX ... maximum allowed step length (input),
C    DS ...... step length (output).
C
      use mtpar
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N)
C  ****  Main-PENELOPE common.
      COMMON/TRACK/E,X,Y,Z,U,V,W,WGHT,KPAR,IBODY,M,ILB(5)
C  ****  Simulation parameters.
      PARAMETER (MAXMAT=10)
      COMMON/CSIMPA/EABS(3,MAXMAT),C1(MAXMAT),C2(MAXMAT),WCC(MAXMAT),
     1  WCR(MAXMAT)
C  ****  Energy grid and interpolation constants for the current energy.
      PARAMETER (NEGP=200)
      COMMON/CEGRID/EL,EU,ET(NEGP),DLEMP(NEGP),DLEMP1,DLFC,
     1  XEL,XE,XEK,KE
C  ****  Electron simulation tables.
      COMMON/CEIMFP/SEHEL(MAXMAT,NEGP),SEHIN(MAXMAT,NEGP),
     1  SEISI(MAXMAT,NEGP),SEHBR(MAXMAT,NEGP),SEAUX(MAXMAT,NEGP),
     2  SETOT(MAXMAT,NEGP),CSTPE(MAXMAT,NEGP),RSTPE(MAXMAT,NEGP),
     3  DEL(MAXMAT,NEGP),W1E(MAXMAT,NEGP),W2E(MAXMAT,NEGP),
     4  RNDCE(MAXMAT,NEGP),AE(MAXMAT,NEGP),BE(MAXMAT,NEGP),
     5  T1E(MAXMAT,NEGP),T2E(MAXMAT,NEGP)
C  ****  Positron simulation tables.
      COMMON/CPIMFP/SPHEL(MAXMAT,NEGP),SPHIN(MAXMAT,NEGP),
     1  SPISI(MAXMAT,NEGP),SPHBR(MAXMAT,NEGP),SPAN(MAXMAT,NEGP),
     2  SPAUX(MAXMAT,NEGP),SPTOT(MAXMAT,NEGP),CSTPP(MAXMAT,NEGP),
     3  RSTPP(MAXMAT,NEGP),W1P(MAXMAT,NEGP),W2P(MAXMAT,NEGP),
     4  RNDCP(MAXMAT,NEGP),AP(MAXMAT,NEGP),BP(MAXMAT,NEGP),
     5  T1P(MAXMAT,NEGP),T2P(MAXMAT,NEGP)
C  ****  Current IMFPs.
      COMMON/CJUMP0/P(8),ST,DST,DS1,W1,W2,T1,T2
      COMMON/CJUMP1/MODE,KSOFTE,KSOFTI,KDELTA
C
      ! EXTERNALRAND
C
      IF(KPAR.EQ.1) THEN
C
C  ************  Electrons (KPAR=1).
C
        XEL=LOG(E)
        XE=1.0D0+(XEL-DLEMP1)*DLFC
        KE=XE
        XEK=XE-KE
        IF(MODE.EQ.1) THEN
          CALL EIMFP(1)
          DS=DS1
          RETURN
        ENDIF
        CALL EIMFP(2)
C
C  ****  Inverse hard mean free path (interaction probability per unit
C        path length).
C
        ST=P(2)+P(3)+P(4)+P(5)+P(8)
        DSMAXP=DSMAX
C
C  ****  Soft stopping interactions.
C        KSOFTI=1, soft stopping is active,
C        KSOFTI=0, soft stopping is not active.
C
        IF(W1.GT.1.0D-20) THEN
          KSOFTI=1
C  ****  The maximum step length, DSMAXP, is determined in terms of the
C        input DSMAX value (which is specified by the user) and the mean
C        free path for hard interactions (1/ST).
          DSMC=4.0D0/ST
          IF(DSMAXP.GT.DSMC) THEN
            DSMAXP=DSMC
          ELSE IF(DSMAXP.LT.1.0D-8) THEN
            DSMAXP=DSMC
          ENDIF
C  ****  The value of DSMAXP is randomized to eliminate dose artifacts
C        at the end of the first step.
        status=rmt(r)
          DSMAXP=(0.5D0+r*0.5D0)*DSMAXP
C
C  ****  Upper bound for the interaction probability along the step
C        (including artificial energy straggling).
C
          EDE0=W1*DSMAXP
          VDE0=W2*DSMAXP
          SPRIME=(W1E(M,KE+1)-W1E(M,KE))/(ET(KE+1)-ET(KE))
          EDEM=EDE0*MAX(1.0D0-0.5D0*SPRIME*EDE0,0.75D0)
          VDEM=VDE0*MAX(1.0D0-(0.5D0*((W2E(M,KE+1)-W2E(M,KE))
     1        /(ET(KE+1)-ET(KE)))+SPRIME)*EDE0,0.75D0)
          W21=VDEM/EDEM
          IF(EDEM.GT.9.0D0*W21) THEN
            ELOWER=MAX(E-(EDEM+3.0D0*SQRT(VDEM)),EABS(1,M))
          ELSE IF(EDEM.GT.3.0D0*W21) THEN
            ELOWER=MAX(E-(EDEM+SQRT(3.0D0*VDEM)),EABS(1,M))
          ELSE
            ELOWER=MAX(E-1.5D0*(EDEM+W21),EABS(1,M))
          ENDIF
          XE1=1.0D0+(LOG(ELOWER)-DLEMP1)*DLFC
          KE1=XE1
          XEK1=XE1-KE1
          STLWR=EXP(SETOT(M,KE1)+(SETOT(M,KE1+1)-SETOT(M,KE1))*XEK1)
          ST=MAX(ST,STLWR)
        ELSE
          KSOFTI=0
        ENDIF
C
C  ****  Soft elastic scattering.
C        KSOFTE=1, soft scattering is active,
C        KSOFTE=0, soft scattering is not active.
C
        IF(T1.GT.1.0D-20) THEN
          KSOFTE=1
        ELSE
          KSOFTE=0
        ENDIF
C
C  ****  Delta interactions.
C        KDELTA=0, a hard interaction follows,
C        KDELTA=1, a delta interaction follows.
C
        status=rmt(r)
        DST=-LOG(r)/ST
        IF(DST.LT.DSMAXP) THEN
          KDELTA=0
        ELSE
          DST=DSMAXP
          KDELTA=1
        ENDIF
C
        IF(KSOFTE+KSOFTI.EQ.0) THEN
          DS=DST
          DS1=0.0D0
          MODE=1
        ELSE
        status=rmt(r)
          DS=DST*r
          DS1=DST-DS
        ENDIF
        RETURN
      ELSE IF(KPAR.EQ.3) THEN
C
C  ************  Positrons (KPAR=3).
C
        XEL=LOG(E)
        XE=1.0D0+(XEL-DLEMP1)*DLFC
        KE=XE
        XEK=XE-KE
        IF(MODE.EQ.1) THEN
          CALL PIMFP(1)
          DS=DS1
          RETURN
        ENDIF
        CALL PIMFP(2)
C
C  ****  Inverse hard mean free path (interaction probability per unit
C        path length).
C
        ST=P(2)+P(3)+P(4)+P(5)+P(6)+P(8)
        DSMAXP=DSMAX
C
C  ****  Soft stopping interactions.
C        KSOFTI=1, soft stopping is active,
C        KSOFTI=0, soft stopping is not active.
C
        IF(W1.GT.1.0D-20) THEN
          KSOFTI=1
C  ****  The maximum step length, DSMAXP, is determined in terms of the
C        input DSMAX value (which is specified by the user) and the mean
C        free path for hard interactions (1/ST).
          DSMC=4.0D0/ST
          IF(DSMAXP.GT.DSMC) THEN
            DSMAXP=DSMC
          ELSE IF(DSMAXP.LT.1.0D-8) THEN
            DSMAXP=DSMC
          ENDIF
C  ****  The value of DSMAXP is randomized to eliminate dose artifacts
C        at the end of the first step.
        status=rmt(r)
          DSMAXP=(0.5D0+r*0.5D0)*DSMAXP
C
C  ****  Upper bound for the interaction probability along the step
C        (including artificial energy straggling).
C
          EDE0=W1*DSMAXP
          VDE0=W2*DSMAXP
          SPRIME=(W1P(M,KE+1)-W1P(M,KE))/(ET(KE+1)-ET(KE))
          EDEM=EDE0*MAX(1.0D0-0.5D0*SPRIME*EDE0,0.75D0)
          VDEM=VDE0*MAX(1.0D0-(0.5D0*((W2P(M,KE+1)-W2P(M,KE))
     1        /(ET(KE+1)-ET(KE)))+SPRIME)*EDE0,0.75D0)
          W21=VDEM/EDEM
          IF(EDEM.GT.9.0D0*W21) THEN
            ELOWER=MAX(E-(EDEM+3.0D0*SQRT(VDEM)),EABS(3,M))
          ELSE IF(EDEM.GT.3.0D0*W21) THEN
            ELOWER=MAX(E-(EDEM+SQRT(3.0D0*VDEM)),EABS(3,M))
          ELSE
            ELOWER=MAX(E-1.5D0*(EDEM+W21),EABS(3,M))
          ENDIF
          XE1=1.0D0+(LOG(ELOWER)-DLEMP1)*DLFC
          KE1=XE1
          XEK1=XE1-KE1
          STLWR=EXP(SPTOT(M,KE1)+(SPTOT(M,KE1+1)-SPTOT(M,KE1))*XEK1)
          ST=MAX(ST,STLWR)
        ELSE
          KSOFTI=0
        ENDIF
C
C  ****  Soft elastic scattering.
C        KSOFTE=1, soft scattering is active,
C        KSOFTE=0, soft scattering is not active.
C
        IF(T1.GT.1.0D-20) THEN
          KSOFTE=1
        ELSE
          KSOFTE=0
        ENDIF
C
C  ****  Delta interactions.
C        KDELTA=0, a hard interaction follows,
C        KDELTA=1, a delta interaction follows.
C
        status=rmt(r)
        DST=-LOG(r)/ST
        IF(DST.LT.DSMAXP) THEN
          KDELTA=0
        ELSE
          DST=DSMAXP
          KDELTA=1
        ENDIF
C
        IF(KSOFTE+KSOFTI.EQ.0) THEN
          DS=DST
          DS1=0.0D0
          MODE=1
        ELSE
        status=rmt(r)
          DS=DST*r
          DS1=DST-DS
        ENDIF
        RETURN
      ELSE
C
C  ************  Photons (KPAR=2).
C
        XEL=LOG(E)
        XE=1.0D0+(XEL-DLEMP1)*DLFC
        KE=XE
        XEK=XE-KE
C
        CALL GIMFP
        ST=P(1)+P(2)+P(3)+P(4)+P(8)
        status=rmt(r)
        DS=-LOG(r)/ST
      ENDIF
      RETURN
      END
C  *********************************************************************
C                       SUBROUTINE KNOCK
C  *********************************************************************
      SUBROUTINE KNOCK(DE,ICOL)
C
C     Simulation of random hinges and hard interaction events.
C
C  Output arguments:
C    DE ..... energy deposited by the particle in the material. It is
C             usually equal to the difference between the energies
C             before and after the interaction.
C    ICOL ... kind of interaction suffered by the particle.
C
      use mtpar
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N)
      CHARACTER*2 LASYMB
      PARAMETER (PI=3.1415926535897932D0, TWOPI=PI+PI)
      PARAMETER (REV=5.10998902D5)  ! Electron rest energy (eV)
      PARAMETER (TREV=2.0D0*REV)
      PARAMETER (TRUNC=1.01538698D0)
C  ****  Main-PENELOPE common.
      COMMON/TRACK/E,X,Y,Z,U,V,W,WGHT,KPAR,IBODY,M,ILB(5)
      COMMON/CHIST/ILBA(5)
C  ****  Simulation parameters.
      PARAMETER (MAXMAT=10)
      COMMON/CSIMPA/EABS(3,MAXMAT),C1(MAXMAT),C2(MAXMAT),WCC(MAXMAT),
     1  WCR(MAXMAT)
C  ****  Composition data.
      COMMON/COMPOS/STF(MAXMAT,30),ZT(MAXMAT),AT(MAXMAT),RHO(MAXMAT),
     1  VMOL(MAXMAT),IZ(MAXMAT,30),NELEM(MAXMAT)
C  ****  Energy grid and interpolation constants for the current energy.
      PARAMETER (NEGP=200)
      COMMON/CEGRID/EL,EU,ET(NEGP),DLEMP(NEGP),DLEMP1,DLFC,
     1  XEL,XE,XEK,KE
C  ****  Element data.
      COMMON/CADATA/ATW(99),EPX(99),RA1(99),RA2(99),RA3(99),RA4(99),
     1  RA5(99),RSCR(99),ETA(99),EB(99,30),IFI(99,30),IKS(99,30),
     2  NSHT(99),LASYMB(99)
C  ****  E/P inelastic collisions.
      PARAMETER (NO=64)
      COMMON/CEIN/EXPOT(MAXMAT),OP2(MAXMAT),F(MAXMAT,NO),UI(MAXMAT,NO),
     1  WRI(MAXMAT,NO),KZ(MAXMAT,NO),KS(MAXMAT,NO),NOSC(MAXMAT)
C  ****  Compton scattering.
      PARAMETER (NOCO=64)
      COMMON/CGCO/FCO(MAXMAT,NOCO),UICO(MAXMAT,NOCO),FJ0(MAXMAT,NOCO),
     2  KZCO(MAXMAT,NOCO),KSCO(MAXMAT,NOCO),NOSCCO(MAXMAT)
C  ****  Bremsstrahlung emission.
      PARAMETER (NBW=32)
      COMMON/CEBR/WB(NBW),PBCUT(MAXMAT,NEGP),WBCUT(MAXMAT,NEGP),
     1  PDFB(MAXMAT,NEGP,NBW),PACB(MAXMAT,NEGP,NBW),ZBR2(MAXMAT)
C  ****  Electron simulation tables.
      COMMON/CEIMFP/SEHEL(MAXMAT,NEGP),SEHIN(MAXMAT,NEGP),
     1  SEISI(MAXMAT,NEGP),SEHBR(MAXMAT,NEGP),SEAUX(MAXMAT,NEGP),
     2  SETOT(MAXMAT,NEGP),CSTPE(MAXMAT,NEGP),RSTPE(MAXMAT,NEGP),
     3  DEL(MAXMAT,NEGP),W1E(MAXMAT,NEGP),W2E(MAXMAT,NEGP),
     4  RNDCE(MAXMAT,NEGP),AE(MAXMAT,NEGP),BE(MAXMAT,NEGP),
     5  T1E(MAXMAT,NEGP),T2E(MAXMAT,NEGP)
C  ****  Positron simulation tables.
      COMMON/CPIMFP/SPHEL(MAXMAT,NEGP),SPHIN(MAXMAT,NEGP),
     1  SPISI(MAXMAT,NEGP),SPHBR(MAXMAT,NEGP),SPAN(MAXMAT,NEGP),
     2  SPAUX(MAXMAT,NEGP),SPTOT(MAXMAT,NEGP),CSTPP(MAXMAT,NEGP),
     3  RSTPP(MAXMAT,NEGP),W1P(MAXMAT,NEGP),W2P(MAXMAT,NEGP),
     4  RNDCP(MAXMAT,NEGP),AP(MAXMAT,NEGP),BP(MAXMAT,NEGP),
     5  T1P(MAXMAT,NEGP),T2P(MAXMAT,NEGP)
C  ****  Current IMFPs.
      COMMON/CJUMP0/P(8),ST,DST,DS1,W1,W2,T1,T2
      COMMON/CJUMP1/MODE,KSOFTE,KSOFTI,KDELTA
C
      COMMON/CELSEP/EELMAX(MAXMAT),PELMAX(MAXMAT)
C
      ! EXTERNALRAND
C
      IF(KPAR.EQ.2) GO TO 2000
      IF(KPAR.EQ.3) GO TO 3000
C
C  ************  Electrons (KPAR=1).
C
 1000 CONTINUE
      IF(MODE.EQ.1) GO TO 1100
C
C  ****  Artificial soft event (ICOL=1).
C
      ICOL=1
      MODE=1
C
      IF(KSOFTI.EQ.0) THEN
        DE=0.0D0
      ELSE
        EDE0=W1*DST
        VDE0=W2*DST
        SPRIME=(W1E(M,KE+1)-W1E(M,KE))/(ET(KE+1)-ET(KE))
        EDE=EDE0*MAX(1.0D0-0.5D0*SPRIME*EDE0,0.75D0)
        VDE=VDE0*MAX(1.0D0-(0.5D0*((W2E(M,KE+1)-W2E(M,KE))
     1       /(ET(KE+1)-ET(KE)))+SPRIME)*EDE0,0.75D0)
C  ****  Generation of random values DE distributed from 0 to infinity
C        with mean EDE and variance VDE.
        SIGMA=SQRT(VDE)
        IF(SIGMA.LT.0.333333333D0*EDE) THEN
C  ****  Sampling from a gaussian distribution (Box-Muller method).
 1001     CONTINUE
        status=rmt(r)
          RN=SQRT(-2.0D0*LOG(r))
        status=rmt(r)
          RN=RN*SIN(TWOPI*r)*TRUNC
          IF(ABS(RN).GE.3.0D0) GO TO 1001
          DE=EDE+RN*SIGMA
        ELSE
        status=rmt(r)
          RU=r
          EDE2=EDE*EDE
          VDE3=3.0D0*VDE
          IF(EDE2.LT.VDE3) THEN
            PNULL=(VDE3-EDE2)/(VDE3+3.0D0*EDE2)
            IF(RU.LT.PNULL) THEN
              DE=0.0D0
            ELSE
              DE=1.5D0*(EDE+VDE/EDE)*(RU-PNULL)/(1.0D0-PNULL)
            ENDIF
          ELSE
            DE=EDE+(2.0D0*RU-1.0D0)*SQRT(VDE3)
          ENDIF
        ENDIF
      ENDIF
C
      E=E-DE
      IF(E.LT.EABS(1,M)) THEN
        DE=E+DE
        E=0.0D0
        RETURN
      ENDIF
      IF(KSOFTE.EQ.0) RETURN
C
C  ****  Bielajew's randomly alternate hinge.
C
        status=rmt(r)
      IF(r.GT.0.5D0.AND.DE.GT.1.0D-3) THEN
        XEL=LOG(E)
        XE=1.0D0+(XEL-DLEMP1)*DLFC
        KE=XE
        XEK=XE-KE
        IF(T1E(M,KE+1).GT.-78.3D0) THEN
          T1=EXP(T1E(M,KE)+(T1E(M,KE+1)-T1E(M,KE))*XEK)
          T2=EXP(T2E(M,KE)+(T2E(M,KE+1)-T2E(M,KE))*XEK)
        ELSE
          T1=0.0D0
          T2=0.0D0
        ENDIF
        IF(T1.LT.1.0D-35) RETURN
      ENDIF
C  ****  1st and 2nd moments of the angular distribution.
      EMU1=0.5D0*(1.0D0-EXP(-DST*T1))
      EMU2=EMU1-(1.0D0-EXP(-DST*T2))/6.0D0
C  ****  Sampling from a two-bar histogram with these moments.
      PNUM=2.0D0*EMU1-3.0D0*EMU2
      PDEN=1.0D0-2.0D0*EMU1
      PB=PNUM/PDEN
      PA=PDEN+PB
        status=rmt(r)
      RND=r
      IF(RND.LT.PA) THEN
        CDT=1.0D0-2.0D0*PB*(RND/PA)
      ELSE
        CDT=1.0D0-2.0D0*(PB+(1.0D0-PB)*((RND-PA)/(1.0D0-PA)))
      ENDIF
        status=rmt(r)
      DF=TWOPI*r
      CALL DIRECT(CDT,DF,U,V,W)
      RETURN
C
C  ************  Hard event.
C
 1100 CONTINUE
      MODE=0
C  ****  A delta interaction (ICOL=7) occurs when the maximum
C        allowed step length is exceeded.
      IF(KDELTA.EQ.1) THEN
        ICOL=7
        DE=0.0D0
        RETURN
      ENDIF
C  ****  Interaction probabilities.
      STNOW=P(2)+P(3)+P(4)+P(5)
C  ****  Random sampling of the interaction type.
        status=rmt(r)
      STS=MAX(STNOW,ST)*r
      SS=P(2)
      IF(SS.GT.STS) GO TO 1200
      SS=SS+P(3)
      IF(SS.GT.STS) GO TO 1300
      SS=SS+P(4)
      IF(SS.GT.STS) GO TO 1400
      SS=SS+P(5)
      IF(SS.GT.STS) GO TO 1500
      SS=SS+P(8)
      IF(SS.GT.STS) GO TO 1800
C  ****  A delta interaction (ICOL=7) may occur when the total
C        interaction probability per unit pathlength, ST, is
C        larger than STNOW.
      ICOL=7
      DE=0.0D0
      RETURN
C
C  ****  Hard elastic collision (ICOL=2).
C
 1200 ICOL=2
      TRNDC=RNDCE(M,KE)+(RNDCE(M,KE+1)-RNDCE(M,KE))*XEK
      IF(E.GE.EELMAX(M)) THEN
        TA=EXP(AE(M,KE)+(AE(M,KE+1)-AE(M,KE))*XEK)
        TB=BE(M,KE)+(BE(M,KE+1)-BE(M,KE))*XEK
        CALL EEL(TA,TB,TRNDC,RMU)
      ELSE
        CALL EELd(TRNDC,RMU)  ! Uses the ELSEP database.
      ENDIF
      CDT=1.0D0-(RMU+RMU)
        status=rmt(r)
      DF=TWOPI*r
      CALL DIRECT(CDT,DF,U,V,W)
      DE=0.0D0
      RETURN
C
C  ****  Hard inelastic collision (ICOL=3).
C
 1300 ICOL=3
      DELTA=DEL(M,KE)+(DEL(M,KE+1)-DEL(M,KE))*XEK
      CALL EIN(E,DELTA,DE,EP,CDT,ES,CDTS,M,IOSC)
C  ****  Scattering angles (primary electron).
        status=rmt(r)
      DF=TWOPI*r
C  ****  Delta ray.
      IF(ES.GT.EABS(1,M)) THEN
        DFS=DF+PI
        US=U
        VS=V
        WS=W
        CALL DIRECT(CDTS,DFS,US,VS,WS)
        ILBA(1)=ILB(1)+1
        ILBA(2)=KPAR
        ILBA(3)=ICOL
        ILBA(4)=0
        ILBA(5)=ILB(5)
        CALL STORES(ES,X,Y,Z,US,VS,WS,WGHT,1,ILBA)
      ENDIF
C  ****  New energy and direction.
      IF(EP.GT.EABS(1,M)) THEN
        E=EP
        CALL DIRECT(CDT,DF,U,V,W)
        RETURN
      ENDIF
      DE=E
      E=0.0D0
      RETURN
C
C  ****  Hard bremsstrahlung emission (ICOL=4).
C
 1400 ICOL=4
      CALL EBR(E,DE,M)
C  ****  Bremsstrahlung photon.
      IF(DE.GT.EABS(2,M)) THEN
        CALL EBRA(E,DE,CDTS,M)
        status=rmt(r)
        DFS=TWOPI*r
        US=U
        VS=V
        WS=W
        CALL DIRECT(CDTS,DFS,US,VS,WS)
        ILBA(1)=ILB(1)+1
        ILBA(2)=KPAR
        ILBA(3)=ICOL
        ILBA(4)=0
        ILBA(5)=ILB(5)
        CALL STORES(DE,X,Y,Z,US,VS,WS,WGHT,2,ILBA)
      ENDIF
C  ****  New energy.
      E=E-DE
      IF(E.GT.EABS(1,M)) RETURN
      DE=E+DE
      E=0.0D0
      RETURN
C
C  ****  Ionization of an inner shell (ICOL=5) -does not affect the
C        primary electron.
C
 1500 ICOL=5
      DE=0.0D0
      CALL ESI(IZA,ISA)
      IF(IZA.LT.1) RETURN
      ILBA(3)=ICOL
      CALL RELAX(IZA,ISA)
      RETURN
C
C  ****  Auxiliary fictitious mechanism (ICOL=8).
C
 1800 ICOL=8
      DE=0.0D0
      CALL EAUX
      RETURN
C
C  ************  Photons (KPAR=2).
C
 2000 CONTINUE
C
        status=rmt(r)
      STS=ST*r
      SS=P(1)
      IF(SS.GT.STS) GO TO 2100
      SS=SS+P(2)
      IF(SS.GT.STS) GO TO 2200
      SS=SS+P(3)
      IF(SS.GT.STS) GO TO 2300
      SS=SS+P(4)
      IF(SS.GT.STS) GO TO 2400
      SS=SS+P(8)
      IF(SS.GT.STS) GO TO 2800
C
C  ****  Rayleigh scattering (ICOL=1).
C
 2100 ICOL=1
      DE=0.0D0
      CALL GRA(E,CDT,M)
        status=rmt(r)
      DF=TWOPI*r
      CALL DIRECT(CDT,DF,U,V,W)
      RETURN
C
C  ****  Compton scattering (ICOL=2).
C
 2200 ICOL=2
      CALL GCO(E,DE,EP,CDT,ES,CDTS,M,ISHELL)
      IZA=KZCO(M,ISHELL)
      ISA=KSCO(M,ISHELL)
        status=rmt(r)
      DF=TWOPI*r
      IF(IZA.GT.0.AND.ISA.LT.10) THEN
        ILBA(3)=ICOL
        CALL RELAX(IZA,ISA)
      ENDIF
C  ****  Compton electron.
      IF(ES.GT.EABS(1,M)) THEN
        DFS=DF+PI
        US=U
        VS=V
        WS=W
        CALL DIRECT(CDTS,DFS,US,VS,WS)
        ILBA(1)=ILB(1)+1
        ILBA(2)=KPAR
        ILBA(3)=ICOL
        ILBA(4)=IZA*1000000+ISA
        ILBA(5)=ILB(5)
        CALL STORES(ES,X,Y,Z,US,VS,WS,WGHT,1,ILBA)
      ENDIF
C  ****  New direction and energy.
      IF(EP.LT.EABS(2,M)) THEN
        DE=E
        E=0.0D0
      ELSE
        CALL DIRECT(CDT,DF,U,V,W)
        E=EP
      ENDIF
      RETURN
C
C  ****  Photoelectric absorption (ICOL=3).
C
 2300 ICOL=3
      CALL GPH(ES,IZA,ISA)
C  ****  Delta interaction. Introduced to correct for the use of an
C        upper bound of the photoelectric attenuation coefficient.
      IF(IZA.EQ.0) THEN
        ICOL=7
        DE=0.0D0
        RETURN
      ENDIF
C
      IF(ES.GT.EABS(1,M)) THEN
        CALL SAUTER(ES,CDTS)
        status=rmt(r)
        DFS=TWOPI*r
        US=U
        VS=V
        WS=W
        CALL DIRECT(CDTS,DFS,US,VS,WS)
        ILBA(1)=ILB(1)+1
        ILBA(2)=KPAR
        ILBA(3)=ICOL
        ILBA(4)=IZA*1000000+ISA
        ILBA(5)=ILB(5)
        CALL STORES(ES,X,Y,Z,US,VS,WS,WGHT,1,ILBA)
      ENDIF
      IF(ISA.LT.10) THEN
        ILBA(3)=ICOL
        CALL RELAX(IZA,ISA)
      ENDIF
      DE=E
      E=0.0D0
      RETURN
C
C  ****  Electron-positron pair production (ICOL=4).
C
 2400 ICOL=4
      CALL GPP(EE,CDTE,EP,CDTP)
      DE=E
C  ****  Electron.
      IF(EE.GT.EABS(1,M)) THEN
        status=rmt(r)
        DF=TWOPI*r
        US=U
        VS=V
        WS=W
        CALL DIRECT(CDTE,DF,US,VS,WS)
        ILBA(1)=ILB(1)+1
        ILBA(2)=KPAR
        ILBA(3)=ICOL
        ILBA(4)=0
        ILBA(5)=ILB(5)
        CALL STORES(EE,X,Y,Z,US,VS,WS,WGHT,1,ILBA)
      ENDIF
C  ****  Positron.
      IF(EP.GT.EABS(3,M)) THEN
        status=rmt(r)
        DF=TWOPI*r
        CALL DIRECT(CDTP,DF,U,V,W)
        ILBA(1)=ILB(1)+1
        ILBA(2)=KPAR
        ILBA(3)=ICOL
        ILBA(4)=0
        ILBA(5)=ILB(5)
        CALL STORES(EP,X,Y,Z,U,V,W,WGHT,3,ILBA)
C  ****  The positron carries a 'latent' energy of 1022 keV.
        DE=DE-TREV
      ELSE
        CALL PANR
      ENDIF
      E=0.0D0
      RETURN
C
C  ****  Auxiliary fictitious mechanism (ICOL=8).
C
 2800 ICOL=8
      DE=0.0D0
      CALL GAUX
      RETURN
C
C  ************  Positrons (KPAR=3).
C
 3000 CONTINUE
      IF(MODE.EQ.1) GO TO 3100
C
C  ****  Artificial soft event (ICOL=1).
C
      ICOL=1
      MODE=1
C
      IF(KSOFTI.EQ.0) THEN
        DE=0.0D0
      ELSE
        EDE0=W1*DST
        VDE0=W2*DST
        SPRIME=(W1P(M,KE+1)-W1P(M,KE))/(ET(KE+1)-ET(KE))
        EDE=EDE0*MAX(1.0D0-0.5D0*SPRIME*EDE0,0.75D0)
        VDE=VDE0*MAX(1.0D0-(0.5D0*((W2P(M,KE+1)-W2P(M,KE))
     1       /(ET(KE+1)-ET(KE)))+SPRIME)*EDE0,0.75D0)
C  ****  Generation of random values DE distributed from 0 to infinity
C        with mean EDE and variance VDE.
        SIGMA=SQRT(VDE)
        IF(SIGMA.LT.0.333333333D0*EDE) THEN
C  ****  Sampling from a gaussian distribution (Box-Muller method).
 3001     CONTINUE
        status=rmt(r)
          RN=SQRT(-2.0D0*LOG(r))
        status=rmt(r)
     1      *SIN(TWOPI*r)*TRUNC
          IF(ABS(RN).GE.3.0D0) GO TO 3001
          DE=EDE+RN*SIGMA
        ELSE
        status=rmt(r)
          RU=r
          EDE2=EDE*EDE
          VDE3=3.0D0*VDE
          IF(EDE2.LT.VDE3) THEN
            PNULL=(VDE3-EDE2)/(VDE3+3.0D0*EDE2)
            IF(RU.LT.PNULL) THEN
              DE=0.0D0
            ELSE
              DE=1.5D0*(EDE+VDE/EDE)*(RU-PNULL)/(1.0D0-PNULL)
            ENDIF
          ELSE
            DE=EDE+(2.0D0*RU-1.0D0)*SQRT(VDE3)
          ENDIF
        ENDIF
      ENDIF
C
      E=E-DE
C  ****  Annihilation at rest.
      IF(E.LT.EABS(3,M)) THEN
        DE=E+DE+TREV
        E=0.0D0
        CALL PANR
        RETURN
      ENDIF
      IF(KSOFTE.EQ.0) RETURN
C
C  ****  Bielajew's randomly alternate hinge.
C
        status=rmt(r)
      IF(r.GT.0.5D0.AND.DE.GT.1.0D-3) THEN
        XEL=LOG(E)
        XE=1.0D0+(XEL-DLEMP1)*DLFC
        KE=XE
        XEK=XE-KE
        IF(T1E(M,KE+1).GT.-78.3D0) THEN
          T1=EXP(T1P(M,KE)+(T1P(M,KE+1)-T1P(M,KE))*XEK)
          T2=EXP(T2P(M,KE)+(T2P(M,KE+1)-T2P(M,KE))*XEK)
        ELSE
          T1=0.0D0
          T2=0.0D0
        ENDIF
        IF(T1.LT.1.0D-35) RETURN
      ENDIF
C  ****  1st and 2nd moments of the angular distribution.
      EMU1=0.5D0*(1.0D0-EXP(-DST*T1))
      EMU2=EMU1-(1.0D0-EXP(-DST*T2))/6.0D0
C  ****  Sampling from a two-bar histogram with these moments.
      PNUM=2.0D0*EMU1-3.0D0*EMU2
      PDEN=1.0D0-2.0D0*EMU1
      PB=PNUM/PDEN
      PA=PDEN+PB
        status=rmt(r)
      RND=r
      IF(RND.LT.PA) THEN
        CDT=1.0D0-2.0D0*PB*(RND/PA)
      ELSE
        CDT=1.0D0-2.0D0*(PB+(1.0D0-PB)*((RND-PA)/(1.0D0-PA)))
      ENDIF
        status=rmt(r)
      DF=TWOPI*r
      CALL DIRECT(CDT,DF,U,V,W)
      RETURN
C
C  ************  Hard event.
C
 3100 CONTINUE
      MODE=0
C  ****  A delta interaction (ICOL=7) occurs when the maximum
C        allowed step length is exceeded.
      IF(KDELTA.EQ.1) THEN
        ICOL=7
        DE=0.0D0
        RETURN
      ENDIF
C  ****  Interaction probabilities.
      STNOW=P(2)+P(3)+P(4)+P(5)+P(6)
C  ****  Random sampling of the interaction type.
        status=rmt(r)
      STS=MAX(STNOW,ST)*r
      SS=P(2)
      IF(SS.GT.STS) GO TO 3200
      SS=SS+P(3)
      IF(SS.GT.STS) GO TO 3300
      SS=SS+P(4)
      IF(SS.GT.STS) GO TO 3400
      SS=SS+P(5)
      IF(SS.GT.STS) GO TO 3500
      SS=SS+P(6)
      IF(SS.GT.STS) GO TO 3600
      SS=SS+P(8)
      IF(SS.GT.STS) GO TO 3800
C  ****  A delta interaction (ICOL=7) may occur when the total
C        interaction probability per unit pathlength, ST, is
C        larger than STNOW.
      ICOL=7
      DE=0.0D0
      RETURN
C
C  ****  Hard elastic collision (ICOL=2).
C
 3200 ICOL=2
      TRNDC=RNDCP(M,KE)+(RNDCP(M,KE+1)-RNDCP(M,KE))*XEK
      IF(E.GE.PELMAX(M)) THEN
        TA=EXP(AP(M,KE)+(AP(M,KE+1)-AP(M,KE))*XEK)
        TB=BP(M,KE)+(BP(M,KE+1)-BP(M,KE))*XEK
        CALL EEL(TA,TB,TRNDC,RMU)
      ELSE
        CALL PELd(TRNDC,RMU)  ! Uses the ELSEP database.
      ENDIF
      CDT=1.0D0-2.0D0*RMU
        status=rmt(r)
      DF=TWOPI*r
      CALL DIRECT(CDT,DF,U,V,W)
      DE=0.0D0
      RETURN
C
C  ****  Hard inelastic collision (ICOL=3).
C
 3300 ICOL=3
      DELTA=DEL(M,KE)+(DEL(M,KE+1)-DEL(M,KE))*XEK
      CALL PIN(E,DELTA,DE,EP,CDT,ES,CDTS,M,IOSC)
C  ****  Scattering angles (primary particle).
        status=rmt(r)
      DF=TWOPI*r
C  ****  Delta ray.
      IF(ES.GT.EABS(1,M)) THEN
        DFS=DF+PI
        US=U
        VS=V
        WS=W
        CALL DIRECT(CDTS,DFS,US,VS,WS)
        ILBA(1)=ILB(1)+1
        ILBA(2)=KPAR
        ILBA(3)=ICOL
        ILBA(4)=0
        ILBA(5)=ILB(5)
        CALL STORES(ES,X,Y,Z,US,VS,WS,WGHT,1,ILBA)
      ENDIF
C  ****  New energy and direction.
      IF(EP.GT.EABS(3,M)) THEN
        E=EP
        CALL DIRECT(CDT,DF,U,V,W)
        RETURN
      ENDIF
      DE=E+TREV
      E=0.0D0
C  ****  Annihilation at rest.
      CALL PANR
      RETURN
C
C  ****  Hard bremsstrahlung emission (ICOL=4).
C
 3400 ICOL=4
      CALL EBR(E,DE,M)
C  ****  Bremsstrahlung photon.
      IF(DE.GT.EABS(2,M)) THEN
        CALL EBRA(E,DE,CDTS,M)
        status=rmt(r)
        DFS=TWOPI*r
        US=U
        VS=V
        WS=W
        CALL DIRECT(CDTS,DFS,US,VS,WS)
        ILBA(1)=ILB(1)+1
        ILBA(2)=KPAR
        ILBA(3)=ICOL
        ILBA(4)=0
        ILBA(5)=ILB(5)
        CALL STORES(DE,X,Y,Z,US,VS,WS,WGHT,2,ILBA)
      ENDIF
C  ****  New energy.
      E=E-DE
      IF(E.GT.EABS(3,M)) RETURN
      DE=E+DE+TREV
      E=0.0D0
C  ****  Annihilation at rest.
      CALL PANR
      RETURN
C
C  ****  Ionization of an inner shell (ICOL=5) -does not affect the
C        primary positron.
C
 3500 ICOL=5
      DE=0.0D0
      CALL PSI(IZA,ISA)
      IF(IZA.LT.1) RETURN
      ILBA(3)=ICOL
      CALL RELAX(IZA,ISA)
      RETURN
C
C  ****  Positron annihilation in flight (ICOL=6).
C
 3600 ICOL=6
      CALL PAN(E1,CDT1,E2,CDT2)
        status=rmt(r)
      DF=TWOPI*r
      IF(E1.GT.EABS(2,M)) THEN
        US=U
        VS=V
        WS=W
        CALL DIRECT(CDT1,DF,US,VS,WS)
        ILBA(1)=ILB(1)+1
        ILBA(2)=KPAR
        ILBA(3)=ICOL
        ILBA(4)=0
        ILBA(5)=ILB(5)
        CALL STORES(E1,X,Y,Z,US,VS,WS,WGHT,2,ILBA)
      ENDIF
      IF(E2.GT.EABS(2,M)) THEN
        DF=DF+PI
        US=U
        VS=V
        WS=W
        CALL DIRECT(CDT2,DF,US,VS,WS)
        ILBA(1)=ILB(1)+1
        ILBA(2)=KPAR
        ILBA(3)=ICOL
        ILBA(4)=0
        ILBA(5)=ILB(5)
        CALL STORES(E2,X,Y,Z,US,VS,WS,WGHT,2,ILBA)
      ENDIF
      DE=E+TREV
      E=0.0D0
      RETURN
C
C  ****  Auxiliary fictitious mechanism (ICOL=8).
C
 3800 ICOL=8
      DE=0.0D0
      CALL PAUX
      RETURN
      END
C  *********************************************************************
C                       SUBROUTINE STORES
C  *********************************************************************
      SUBROUTINE STORES(EI,XI,YI,ZI,UI,VI,WI,WGHTI,KPARI,ILBI)
C
C  This subroutine stores the initial state of a new secondary particle
C  in the secondary stack. The input values are:
C     EI ........... initial energy,
C     XI, YI, ZI ... initial position coordinates.
C     UI, VI, WI ... initial direction cosines,
C     WGHTI ........ weight (=1 in analogue simulation),
C     KPARI ........ kind of particle (1: electron, 2: photon,
C                    3: positron).
C     ILBI(5) ...... particle labels.
C
C  The parameter NMS fixes the size of the secondary stack (i.e. the
C  maximum number of particles that can be stored). If this number is
C  exceeded, a warning message is printed on unit 26. When the memory
C  storage is exhausted, each new secondary particle is stored on the
C  position of the less energetic secondary electron or photon already
C  produced, which is thus discarded.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N)
      PARAMETER (NMS=1000)
      DIMENSION ILBI(5)
C  ****  Main-PENELOPE common.
      COMMON/TRACK/E,X,Y,Z,U,V,W,WGHT,KPAR,IBODY,M,ILB(5)
C  ****  Secondary stack.
      COMMON/SECST/ES(NMS),XS(NMS),YS(NMS),ZS(NMS),US(NMS),
     1   VS(NMS),WS(NMS),WGHTS(NMS),KS(NMS),IBODYS(NMS),MS(NMS),
     2   ILBS(5,NMS),NSEC
      COMMON/CERSEC/IERSEC
C
      IF(NSEC.LT.NMS) THEN
        NSEC=NSEC+1
        IS=NSEC
      ELSE
        IF(IERSEC.EQ.0) THEN
          WRITE(26,1001)
          WRITE(6,1001)
 1001     FORMAT(/3X,'*** WARNING: (STORES) not enough storage for ',
     1      'secondaries.',/16X,'EABS(KPAR,m) or the parameter NMS ',
     2      'should be enlarged.')
          IERSEC=1
        ENDIF
        NSEC=NMS
        EME=1.0D35
        EMG=1.0D35
        IE=0
        IG=0
        DO 1 I=1,NMS
        IF(KS(I).EQ.1) THEN
          IF(ES(I).LT.EME) THEN
            EME=ES(I)
            IE=I
          ENDIF
        ELSE IF(KS(I).EQ.2) THEN
          IF (ES(I).LT.EMG) THEN
            EMG=ES(I)
            IG=I
          ENDIF
        ENDIF
    1   CONTINUE
C
        IF(IE.GT.0) THEN
          IS=IE
        ELSE IF(IG.GT.0) THEN
          IS=IG
        ELSE
          WRITE(26,1002)
 1002     FORMAT(/3X,'*** Not enough storage for secondary positrons.',
     1      /7X,'JOB ABORTED.')
          STOP 'STORES. Not enough storage for secondary positrons.'
        ENDIF
C
        DO I=IS,NMS-1
          ES(I)=ES(I+1)
          XS(I)=XS(I+1)
          YS(I)=YS(I+1)
          ZS(I)=ZS(I+1)
          US(I)=US(I+1)
          VS(I)=VS(I+1)
          WS(I)=WS(I+1)
          WGHTS(I)=WGHTS(I+1)
          KS(I)=KS(I+1)
          IBODYS(I)=IBODYS(I+1)
          MS(I)=MS(I+1)
          DO IB=1,5
            ILBS(IB,I)=ILBS(IB,I+1)
          ENDDO
        ENDDO
        IS=NMS
      ENDIF
C
      ES(IS)=EI
      XS(IS)=XI
      YS(IS)=YI
      ZS(IS)=ZI
      US(IS)=UI
      VS(IS)=VI
      WS(IS)=WI
      WGHTS(IS)=WGHTI
      KS(IS)=KPARI
      IBODYS(IS)=IBODY
      MS(IS)=M
      ILBS(1,IS)=ILBI(1)
      ILBS(2,IS)=ILBI(2)
      ILBS(3,IS)=ILBI(3)
      ILBS(4,IS)=ILBI(4)
      ILBS(5,IS)=ILBI(5)
      RETURN
      END
C  *********************************************************************
C                       SUBROUTINE SECPAR
C  *********************************************************************
      SUBROUTINE SECPAR(LEFT)
C
C  This subroutine delivers the initial state of a secondary particle
C  produced during the previous simulation of the shower. This particle
C  is removed from the secondary stack, so that it will be lost if a new
C  call to SECPAR is performed before simulating its trajectory up to
C  the end.
C
C  LEFT is the number of particles in the secondary stack at the calling
C  time. When LEFT=0, the simulation of the shower has been completed.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N)
C  ****  Main-PENELOPE common.
      COMMON/TRACK/E,X,Y,Z,U,V,W,WGHT,KPAR,IBODY,M,ILB(5)
C  ****  Secondary stack.
      PARAMETER (NMS=1000)
      COMMON/SECST/ES(NMS),XS(NMS),YS(NMS),ZS(NMS),US(NMS),
     1   VS(NMS),WS(NMS),WGHTS(NMS),KS(NMS),IBODYS(NMS),MS(NMS),
     2   ILBS(5,NMS),NSEC
C
      IF(NSEC.GT.0) THEN
        LEFT=NSEC
        E=ES(NSEC)
        X=XS(NSEC)
        Y=YS(NSEC)
        Z=ZS(NSEC)
        U=US(NSEC)
        V=VS(NSEC)
        W=WS(NSEC)
        WGHT=WGHTS(NSEC)
        KPAR=KS(NSEC)
        IBODY=IBODYS(NSEC)
        M=MS(NSEC)
        DO I=1,5
          ILB(I)=ILBS(I,NSEC)
        ENDDO
        NSEC=NSEC-1
      ELSE
        LEFT=0
      ENDIF
      RETURN
      END
C  *********************************************************************
C                       SUBROUTINE EIMFP
C  *********************************************************************
      SUBROUTINE EIMFP(IEND)
C
C  This subroutine computes the inverse mean free paths for hard inter-
C  actions of electrons with the current energy in material M.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N)
C  ****  Main-PENELOPE common.
      COMMON/TRACK/E,X,Y,Z,U,V,W,WGHT,KPAR,IBODY,M,ILB(5)
C  ****  Energy grid and interpolation constants for the current energy.
      PARAMETER (NEGP=200)
      COMMON/CEGRID/EL,EU,ET(NEGP),DLEMP(NEGP),DLEMP1,DLFC,
     1  XEL,XE,XEK,KE
C  ****  Electron simulation tables.
      PARAMETER (MAXMAT=10)
      COMMON/CEIMFP/SEHEL(MAXMAT,NEGP),SEHIN(MAXMAT,NEGP),
     1  SEISI(MAXMAT,NEGP),SEHBR(MAXMAT,NEGP),SEAUX(MAXMAT,NEGP),
     2  SETOT(MAXMAT,NEGP),CSTPE(MAXMAT,NEGP),RSTPE(MAXMAT,NEGP),
     3  DEL(MAXMAT,NEGP),W1E(MAXMAT,NEGP),W2E(MAXMAT,NEGP),
     4  RNDCE(MAXMAT,NEGP),AE(MAXMAT,NEGP),BE(MAXMAT,NEGP),
     5  T1E(MAXMAT,NEGP),T2E(MAXMAT,NEGP)
C  ****  Current IMFPs.
      COMMON/CJUMP0/P(8),ST,DST,DS1,W1,W2,T1,T2
C
      P(2)=EXP(SEHEL(M,KE)+(SEHEL(M,KE+1)-SEHEL(M,KE))*XEK)
      P(3)=EXP(SEHIN(M,KE)+(SEHIN(M,KE+1)-SEHIN(M,KE))*XEK)
      P(4)=EXP(SEHBR(M,KE)+(SEHBR(M,KE+1)-SEHBR(M,KE))*XEK)
      P(5)=EXP(SEISI(M,KE)+(SEISI(M,KE+1)-SEISI(M,KE))*XEK)
      P(8)=EXP(SEAUX(M,KE)+(SEAUX(M,KE+1)-SEAUX(M,KE))*XEK)
      IF(IEND.EQ.1) RETURN
      IF(W1E(M,KE+1).GT.-78.3D0) THEN
        W1=EXP(W1E(M,KE)+(W1E(M,KE+1)-W1E(M,KE))*XEK)
        W2=EXP(W2E(M,KE)+(W2E(M,KE+1)-W2E(M,KE))*XEK)
      ELSE
        W1=0.0D0
        W2=0.0D0
      ENDIF
      IF(T1E(M,KE+1).GT.-78.3D0) THEN
        T1=EXP(T1E(M,KE)+(T1E(M,KE+1)-T1E(M,KE))*XEK)
        T2=EXP(T2E(M,KE)+(T2E(M,KE+1)-T2E(M,KE))*XEK)
      ELSE
        T1=0.0D0
        T2=0.0D0
      ENDIF
      RETURN
      END
C  *********************************************************************
C                       SUBROUTINE PIMFP
C  *********************************************************************
      SUBROUTINE PIMFP(IEND)
C
C  This subroutine computes the inverse mean free paths for hard inter-
C  actions of positrons with the current energy in material M.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N)
C  ****  Main-PENELOPE common.
      COMMON/TRACK/E,X,Y,Z,U,V,W,WGHT,KPAR,IBODY,M,ILB(5)
C  ****  Energy grid and interpolation constants for the current energy.
      PARAMETER (NEGP=200)
      COMMON/CEGRID/EL,EU,ET(NEGP),DLEMP(NEGP),DLEMP1,DLFC,
     1  XEL,XE,XEK,KE
C  ****  Positron simulation tables.
      PARAMETER (MAXMAT=10)
      COMMON/CPIMFP/SPHEL(MAXMAT,NEGP),SPHIN(MAXMAT,NEGP),
     1  SPISI(MAXMAT,NEGP),SPHBR(MAXMAT,NEGP),SPAN(MAXMAT,NEGP),
     2  SPAUX(MAXMAT,NEGP),SPTOT(MAXMAT,NEGP),CSTPP(MAXMAT,NEGP),
     3  RSTPP(MAXMAT,NEGP),W1P(MAXMAT,NEGP),W2P(MAXMAT,NEGP),
     4  RNDCP(MAXMAT,NEGP),AP(MAXMAT,NEGP),BP(MAXMAT,NEGP),
     5  T1P(MAXMAT,NEGP),T2P(MAXMAT,NEGP)
C  ****  Current IMFPs.
      COMMON/CJUMP0/P(8),ST,DST,DS1,W1,W2,T1,T2
C
      P(2)=EXP(SPHEL(M,KE)+(SPHEL(M,KE+1)-SPHEL(M,KE))*XEK)
      P(3)=EXP(SPHIN(M,KE)+(SPHIN(M,KE+1)-SPHIN(M,KE))*XEK)
      P(4)=EXP(SPHBR(M,KE)+(SPHBR(M,KE+1)-SPHBR(M,KE))*XEK)
      P(5)=EXP(SPISI(M,KE)+(SPISI(M,KE+1)-SPISI(M,KE))*XEK)
      P(6)=EXP(SPAN(M,KE)+(SPAN(M,KE+1)-SPAN(M,KE))*XEK)
      P(8)=EXP(SPAUX(M,KE)+(SPAUX(M,KE+1)-SPAUX(M,KE))*XEK)
      IF(IEND.EQ.1) RETURN
      IF(W1P(M,KE+1).GT.-78.3D0) THEN
        W1=EXP(W1P(M,KE)+(W1P(M,KE+1)-W1P(M,KE))*XEK)
        W2=EXP(W2P(M,KE)+(W2P(M,KE+1)-W2P(M,KE))*XEK)
      ELSE
        W1=0.0D0
        W2=0.0D0
      ENDIF
      IF(T1P(M,KE+1).GT.-78.3D0) THEN
        T1=EXP(T1P(M,KE)+(T1P(M,KE+1)-T1P(M,KE))*XEK)
        T2=EXP(T2P(M,KE)+(T2P(M,KE+1)-T2P(M,KE))*XEK)
      ELSE
        T1=0.0D0
        T2=0.0D0
      ENDIF
      RETURN
      END
C  *********************************************************************
C                       SUBROUTINE GIMF
C  *********************************************************************
      SUBROUTINE GIMFP
C
C  This subroutine computes the inverse mean free paths for interactions
C  of photons with the current energy in material M.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N)
C  ****  Main-PENELOPE common.
      COMMON/TRACK/E,X,Y,Z,U,V,W,WGHT,KPAR,IBODY,M,ILB(5)
C  ****  Energy grid and interpolation constants for the current energy.
      PARAMETER (NEGP=200)
      COMMON/CEGRID/EL,EU,ET(NEGP),DLEMP(NEGP),DLEMP1,DLFC,
     1  XEL,XE,XEK,KE
C  ****  Photon simulation tables.
      PARAMETER (MAXMAT=10)
      COMMON/CGIMFP/SGRA(MAXMAT,NEGP),SGCO(MAXMAT,NEGP),
     1  SGPH(MAXMAT,NEGP),SGPP(MAXMAT,NEGP),SGAUX(MAXMAT,NEGP)
C  ****  Current IMFPs.
      COMMON/CJUMP0/P(8),ST,DST,DS1,W1,W2,T1,T2
C
      P(1)=EXP(SGRA(M,KE)+(SGRA(M,KE+1)-SGRA(M,KE))*XEK)
      P(2)=EXP(SGCO(M,KE)+(SGCO(M,KE+1)-SGCO(M,KE))*XEK)
      P(3)=SGPH(M,KE)
      IF(E.LT.1.023D6) THEN
        P(4)=0.0D0
      ELSE
        P(4)=EXP(SGPP(M,KE)+(SGPP(M,KE+1)-SGPP(M,KE))*XEK)
      ENDIF
      P(8)=EXP(SGAUX(M,KE)+(SGAUX(M,KE+1)-SGAUX(M,KE))*XEK)
      RETURN
      END
C  *********************************************************************
C                       SUBROUTINE EAUX
C  *********************************************************************
      SUBROUTINE EAUX
C
C  Auxiliary interaction mechanism for electrons, definable by the user.
C  Usually it is not active.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N)
C  ****  Main-PENELOPE common.
      COMMON/TRACK/E,X,Y,Z,U,V,W,WGHT,KPAR,IBODY,M,ILB(5)
C  ****  Simulation parameters.
      PARAMETER (MAXMAT=10)
      COMMON/CSIMPA/EABS(3,MAXMAT),C1(MAXMAT),C2(MAXMAT),WCC(MAXMAT),
     1  WCR(MAXMAT)
C
      WRITE(26,1000)
 1000 FORMAT(1X,'Warning: Subroutine EAUX has been entered.')
      RETURN
      END
C  *********************************************************************
C                       SUBROUTINE PAUX
C  *********************************************************************
      SUBROUTINE PAUX
C
C  Auxiliary interaction mechanism for positrons, definable by the user.
C  Usually it is not active.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N)
C  ****  Main-PENELOPE common.
      COMMON/TRACK/E,X,Y,Z,U,V,W,WGHT,KPAR,IBODY,M,ILB(5)
C  ****  Simulation parameters.
      PARAMETER (MAXMAT=10)
      COMMON/CSIMPA/EABS(3,MAXMAT),C1(MAXMAT),C2(MAXMAT),WCC(MAXMAT),
     1  WCR(MAXMAT)
C
      WRITE(26,1000)
 1000 FORMAT(1X,'Warning: Subroutine PAUX has been entered.')
      RETURN
      END
C  *********************************************************************
C                       SUBROUTINE GAUX
C  *********************************************************************
      SUBROUTINE GAUX
C
C  Auxiliary interaction mechanism for photons, definable by the user.
C  Usually it is not active.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N)
C  ****  Main-PENELOPE common.
      COMMON/TRACK/E,X,Y,Z,U,V,W,WGHT,KPAR,IBODY,M,ILB(5)
C  ****  Simulation parameters.
      PARAMETER (MAXMAT=10)
      COMMON/CSIMPA/EABS(3,MAXMAT),C1(MAXMAT),C2(MAXMAT),WCC(MAXMAT),
     1  WCR(MAXMAT)
C
      WRITE(26,1000)
 1000 FORMAT(1X,'Warning: Subroutine GAUX has been entered.')
      RETURN
      END


CXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
CXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX


C
C     *********************************
C     *  SUBROUTINE PACKAGE PENELAST  *
C     *********************************
C
C     The following subroutines perform class II simulation of elastic
C  scattering of electrons and positrons using the numerical cross
C  sections of the ELSEP database, which covers the energy range from
C  50 eV to 100 MeV.
C
C   ****  DO NOT USE THESE SUBROUTINES FOR ENERGIES ABOVE 100 MeV  ****
C
C                                         Francesc Salvat. October 2004.
C
C  *********************************************************************
C                       SUBROUTINE EELdW
C  *********************************************************************
      SUBROUTINE EELdW(M,IWR)
C
C  This subroutine generates a table of differential cross sections for
C  elastic scattering of electrons and positrons in material M, and
C  writes it on the material definition file. Data are read from the
C  ELSEP elastic database files.
C 
      use mtpar
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N)
      PARAMETER (PI=3.1415926535897932D0)
C  ****  Composition data.
      PARAMETER (MAXMAT=10)
      COMMON/COMPOS/STF(MAXMAT,30),ZT(MAXMAT),AT(MAXMAT),RHO(MAXMAT),
     1  VMOL(MAXMAT),IZ(MAXMAT,30),NELEM(MAXMAT)
      DIMENSION IZM(30),STFM(30)
C
C  ****  Elastic scattering simulation tables.
C
      PARAMETER (NE=96,NA=606)
      COMMON/CDCSEP/ETS(NE),ETL(NE),TH(NA),THR(NA),XMU(NA),XMUL(NA),
     1       ECS(NE),ETCS1(NE),ETCS2(NE),EDCS(NE,NA),
     2       PCS(NE),PTCS1(NE),PTCS2(NE),PDCS(NE,NA),
     3       DCSI(NA),DCSIL(NA),CSI,TCS1I,TCS2I
C
      DO I=1,NELEM(M)
        IZM(I)=IZ(M,I)
        STFM(I)=STF(M,I)
      ENDDO
      CALL ELINIT(IZM,STFM,NELEM(M))
C
C  ************  Write final DCS tables.
C
C  ****  Electrons.
C
      IELEC=-1
      WRITE(IWR,2001)
 2001 FORMAT(1X,'***  Electron elastic differential cross sections')
      DO IE=1,NE
        DO  K=1,NA
          DCSI(K)=EDCS(IE,K)
        ENDDO
        ECS0=4.0D0*PI*RMOMX(XMU,DCSI,0.0D0,1.0D0,NA,0)
        ECS1=4.0D0*PI*RMOMX(XMU,DCSI,0.0D0,1.0D0,NA,1)
        ECS2=4.0D0*PI*RMOMX(XMU,DCSI,0.0D0,1.0D0,NA,2)
        TCS1=2.0D0*ECS1
        TCS2=6.0D0*(ECS1-ECS2)
        WRITE(IWR,'(I3,1P,E10.3,5E12.5)')
     1    IELEC,ETS(IE),ECS0,TCS1,TCS2
        WRITE(IWR,'(1P,10E12.5)') (EDCS(IE,K),K=1,NA)
C  ****  Consistency test.
        TS0=(ECS0-ECS(IE))/ECS(IE)
        TS1=(TCS1-ETCS1(IE))/ETCS1(IE)
        TS2=(TCS2-ETCS2(IE))/ETCS2(IE)
        TSTE=MAX(ABS(TS0),ABS(TS1),ABS(TS2))
        IF(TSTE.GT.1.0D-2) THEN
          WRITE(IWR,'('' E='',1P,E12.5)') ETS(IE)
          WRITE(IWR,'(3X,3E12.5)') ECS0,TCS1,TCS2
          WRITE(IWR,'(3X,3E12.5)') ECS(IE),ETCS1(IE),ETCS2(IE)
          WRITE(IWR,*) ' Electron cross section data are corrupt.'
          WRITE(6,'('' E='',1P,E12.5)') ETS(IE)
          STOP ' Electron cross section data are corrupt.'
        ENDIF
      ENDDO
C
C  ****  Positrons.
C
      IELEC=+1
      WRITE(IWR,2002)
 2002 FORMAT(1X,'***  Positron elastic differential cross sections')
      DO IE=1,NE
        DO  K=1,NA
          DCSI(K)=PDCS(IE,K)
        ENDDO
        ECS0=4.0D0*PI*RMOMX(XMU,DCSI,0.0D0,1.0D0,NA,0)
        ECS1=4.0D0*PI*RMOMX(XMU,DCSI,0.0D0,1.0D0,NA,1)
        ECS2=4.0D0*PI*RMOMX(XMU,DCSI,0.0D0,1.0D0,NA,2)
        TCS1=2.0D0*ECS1
        TCS2=6.0D0*(ECS1-ECS2)
        WRITE(IWR,'(I3,1P,E10.3,5E12.5)')
     1    IELEC,ETS(IE),ECS0,TCS1,TCS2
        WRITE(IWR,'(1P,10E12.5)') (PDCS(IE,K),K=1,NA)
C  ****  Consistency test.
        TS0=(ECS0-PCS(IE))/PCS(IE)
        TS1=(TCS1-PTCS1(IE))/PTCS1(IE)
        TS2=(TCS2-PTCS2(IE))/PTCS2(IE)
        TSTE=MAX(ABS(TS0),ABS(TS1),ABS(TS2))
        IF(TSTE.GT.1.0D-2) THEN
          WRITE(IWR,'('' E='',1P,E12.5)') ETS(IE)
          WRITE(IWR,'(3X,3E12.5)') ECS0,TCS1,TCS2
          WRITE(IWR,'(3X,3E12.5)') PCS(IE),PTCS1(IE),PTCS2(IE)
          WRITE(IWR,*) ' Positron cross section data are corrupt.'
          WRITE(6,'('' E='',1P,E12.5)') ETS(IE)
          STOP ' Positron cross section data are corrupt.'
        ENDIF
      ENDDO
      RETURN
      END
C  *********************************************************************
C                       SUBROUTINE EELdR
C  *********************************************************************
      SUBROUTINE EELdR(M,IRD,IWR,INFO)
C
C     This subroutine reads elastic cross sections for electrons and
c  positrons in material M from the elastic scattering database. It also
C  initializes the algorithm for simulation of elastic collisions of
C  electrons and positrons.
C
      use mtpar
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N)
      CHARACTER*50 CTEXT
      PARAMETER (PI=3.1415926535897932D0)
C  ****  Composition data.
      PARAMETER (MAXMAT=10)
      COMMON/COMPOS/STF(MAXMAT,30),ZT(MAXMAT),AT(MAXMAT),RHO(MAXMAT),
     1  VMOL(MAXMAT),IZ(MAXMAT,30),NELEM(MAXMAT)
C  ****  Energy grid and interpolation constants for the current energy.
      PARAMETER (NEGP=200)
      COMMON/CEGRID/EL,EU,ET(NEGP),DLEMP(NEGP),DLEMP1,DLFC,
     1  XEL,XE,XEK,KE
C  ****  Simulation parameters.
      COMMON/CSIMPA/EABS(3,MAXMAT),C1(MAXMAT),C2(MAXMAT),WCC(MAXMAT),
     1  WCR(MAXMAT)
      COMMON/CEINTF/T1EI(NEGP),T2EI(NEGP),T1PI(NEGP),T2PI(NEGP)
C  ****  Electron simulation tables.
      COMMON/CEIMFP/SEHEL(MAXMAT,NEGP),SEHIN(MAXMAT,NEGP),
     1  SEISI(MAXMAT,NEGP),SEHBR(MAXMAT,NEGP),SEAUX(MAXMAT,NEGP),
     2  SETOT(MAXMAT,NEGP),CSTPE(MAXMAT,NEGP),RSTPE(MAXMAT,NEGP),
     3  DEL(MAXMAT,NEGP),W1E(MAXMAT,NEGP),W2E(MAXMAT,NEGP),
     4  RNDCE(MAXMAT,NEGP),AE(MAXMAT,NEGP),BE(MAXMAT,NEGP),
     5  T1E(MAXMAT,NEGP),T2E(MAXMAT,NEGP)
C  ****  Positron simulation tables.
      COMMON/CPIMFP/SPHEL(MAXMAT,NEGP),SPHIN(MAXMAT,NEGP),
     1  SPISI(MAXMAT,NEGP),SPHBR(MAXMAT,NEGP),SPAN(MAXMAT,NEGP),
     2  SPAUX(MAXMAT,NEGP),SPTOT(MAXMAT,NEGP),CSTPP(MAXMAT,NEGP),
     3  RSTPP(MAXMAT,NEGP),W1P(MAXMAT,NEGP),W2P(MAXMAT,NEGP),
     4  RNDCP(MAXMAT,NEGP),AP(MAXMAT,NEGP),BP(MAXMAT,NEGP),
     5  T1P(MAXMAT,NEGP),T2P(MAXMAT,NEGP)
C  ****  Total and transport cross sections.
      COMMON/CEEL00/EIT(NEGP),XE0(NEGP),XE1(NEGP),XE2(NEGP),XP0(NEGP),
     1  XP1(NEGP),XP2(NEGP),T1E0(NEGP),T2E0(NEGP),T1P0(NEGP),T2P0(NEGP),
     2  EITL(NEGP),FL(NEGP),A(NEGP),B(NEGP),C(NEGP),D(NEGP)
C
C  ****  Elastic scattering simulation tables.
C
      PARAMETER (NE=96,NA=606)
      COMMON/CDCSEP/ETS(NE),ETL(NE),TH(NA),THR(NA),XMU(NA),XMUL(NA),
     1       ECS(NE),ETCS1(NE),ETCS2(NE),EDCS(NE,NA),
     2       PCS(NE),PTCS1(NE),PTCS2(NE),PDCS(NE,NA),
     3       DCSI(NA),DCSIL(NA),CSI,TCS1I,TCS2I
C
      PARAMETER (NM=512)
      COMMON/CQRND/XTI(NM),PACI(NM),AI(NM),BI(NM),NPI
      COMMON/CQRND1/ITTLI(NM),ITTUI(NM),NPM1I
C
      PARAMETER (NP=128)
      COMMON/CEELDB/XSE(NP,NEGP,MAXMAT),PSE(NP,NEGP,MAXMAT),
     1              ASE(NP,NEGP,MAXMAT),BSE(NP,NEGP,MAXMAT),
     2              ITLE(NP,NEGP,MAXMAT),ITUE(NP,NEGP,MAXMAT)
      COMMON/CPELDB/XSP(NP,NEGP,MAXMAT),PSP(NP,NEGP,MAXMAT),
     1              ASP(NP,NEGP,MAXMAT),BSP(NP,NEGP,MAXMAT),
     2              ITLP(NP,NEGP,MAXMAT),ITUP(NP,NEGP,MAXMAT)
      COMMON/CELSEP/EELMAX(MAXMAT),PELMAX(MAXMAT)
C
      DIMENSION EGRD(16)
      DATA EGRD/1.0D0,1.25D0,1.50D0,1.75D0,2.00D0,2.50D0,3.00D0,
     1  3.50D0,4.00D0,4.50D0,5.00D0,6.00D0,7.00D0,8.00D0,9.00D0,
     2  1.00D1/
C
      EXTERNAL DCSEL
C
C  ****  Energy mesh points (in eV).
C
      IE=0
      IGRID=10
      FGRID=10.0D0
   10 IGRID=IGRID+1
      EV=EGRD(IGRID)*FGRID
      IF(IGRID.EQ.16) THEN
        IGRID=1
        FGRID=10.0D0*FGRID
      ENDIF
      IE=IE+1
      ETS(IE)=EV
      ETL(IE)=LOG(ETS(IE))
C     WRITE(6,'(I5,1P,3E15.7)') IE,ETS(IE)
      IF(IE.LT.NE) GO TO 10
C
C  ****  Angular grid (TH in deg, XMU=(1.0D0-COS(TH))/2).
C
      I=1
      TH(I)=0.0D0
      THR(I)=TH(I)*PI/180.0D0
      XMU(I)=(1.0D0-COS(THR(I)))/2.0D0
      XMUL(I)=LOG(1.0D-35)
      I=2
      TH(I)=1.0D-4
      THR(I)=TH(I)*PI/180.0D0
      XMU(I)=(1.0D0-COS(THR(I)))/2.0D0
      XMUL(I)=LOG(MAX(XMU(I),1.0D-35))
   20 CONTINUE
      I=I+1
      IF(TH(I-1).LT.0.9999D-3) THEN
        TH(I)=TH(I-1)+2.5D-5
      ELSE IF(TH(I-1).LT.0.9999D-2) THEN
        TH(I)=TH(I-1)+2.5D-4
      ELSE IF(TH(I-1).LT.0.9999D-1) THEN
        TH(I)=TH(I-1)+2.5D-3
      ELSE IF(TH(I-1).LT.0.9999D+0) THEN
        TH(I)=TH(I-1)+2.5D-2
      ELSE IF(TH(I-1).LT.0.9999D+1) THEN
        TH(I)=TH(I-1)+1.0D-1
      ELSE IF(TH(I-1).LT.2.4999D+1) THEN
        TH(I)=TH(I-1)+2.5D-1
      ELSE
        TH(I)=TH(I-1)+5.0D-1
      ENDIF
      THR(I)=TH(I)*PI/180.0D0
      XMU(I)=MAX((1.0D0-COS(THR(I)))/2.0D0,1.0D-35)
      XMUL(I)=LOG(MAX(XMU(I),1.0D-35))
      IF(I.LT.NA) GO TO 20
C     DO I=1,NA
C       WRITE(6,'(I5,1P,3E15.7)') I,TH(I),THR(I),XMU(I)
C     ENDDO
C
C  ****  Read elastic DCS tables.
C
      IF(INFO.GE.3) WRITE(IWR,2001)
 2001 FORMAT(1X,'***  Electron elastic differential cross sections')
      IELEC=-1
      READ(IRD,'(A50)') CTEXT
      DO IE=1,NE
        READ(IRD,'(I3,1P,E10.3,5E12.5)')
     1    IELEC,ETSIE,ECS(IE),ETCS1(IE),ETCS2(IE)
C       WRITE(6,'(I5,1P,3E15.7)') IE,ETSIE,ECS(IE)
        IF(INFO.GE.3) WRITE(IWR,'(I3,1P,E10.3,5E12.5)')
     1    IELEC,ETS(IE),ECS(IE),ETCS1(IE),ETCS2(IE)
        IF(IELEC.NE.-1.OR.ABS(ETSIE-ETS(IE)).GT.0.1D0) THEN
          WRITE(IWR,*) ' Error reading electron elastic DCS data.'
          STOP 'Error reading electron elastic DCS data.'
        ENDIF
        READ(IRD,'(1P,10E12.5)') (EDCS(IE,K),K=1,NA)
        IF(INFO.GE.3) WRITE(IWR,'(1P,10E12.5)') (EDCS(IE,K),K=1,NA)
C  ****  Consistency test.
        DO K=1,NA
          DCSI(K)=EDCS(IE,K)
        ENDDO
        ECS0=4.0D0*PI*RMOMX(XMU,DCSI,0.0D0,1.0D0,NA,0)
        ECS1=4.0D0*PI*RMOMX(XMU,DCSI,0.0D0,1.0D0,NA,1)
        ECS2=4.0D0*PI*RMOMX(XMU,DCSI,0.0D0,1.0D0,NA,2)
        TCS1=2.0D0*ECS1
        TCS2=6.0D0*(ECS1-ECS2)
        TS0=(ECS0-ECS(IE))/ECS(IE)
        TS1=(TCS1-ETCS1(IE))/ETCS1(IE)
        TS2=(TCS2-ETCS2(IE))/ETCS2(IE)
        TSTE=MAX(ABS(TS0),ABS(TS1),ABS(TS2))
        IF(TSTE.GT.1.0D-4) THEN
          WRITE(IWR,'('' E='',1P,E12.5)') ETS(IE)
          WRITE(6,'('' E='',1P,E12.5)') ETS(IE)
          WRITE(IWR,*) ' Electron cross section data are corrupt.'
          STOP ' Electron cross section data are corrupt.'
        ENDIF
      ENDDO
C
      IF(INFO.GE.3) WRITE(IWR,2002)
 2002 FORMAT(1X,'***  Positron elastic differential cross sections')
      IELEC=+1
      READ(IRD,'(A50)') CTEXT
      DO IE=1,NE
        READ(IRD,'(I3,1P,E10.3,5E12.5)')
     1    IELEC,ETSIE,PCS(IE),PTCS1(IE),PTCS2(IE)
C       WRITE(6,'(I5,1P,3E15.7)') IE,ETSIE,PCS(IE)
        IF(INFO.GE.3) WRITE(IWR,'(I3,1P,E10.3,5E12.5)')
     1    IELEC,ETS(IE),PCS(IE),PTCS1(IE),PTCS2(IE)
        IF(IELEC.NE.+1.OR.ABS(ETSIE-ETS(IE)).GT.0.1D0) THEN
          WRITE(IWR,*) ' Error reading positron elastic DCS data.'
          STOP 'Error reading positron elastic DCS data.'
        ENDIF
        READ(IRD,'(1P,10E12.5)') (PDCS(IE,K),K=1,NA)
        IF(INFO.GE.3) WRITE(IWR,'(1P,10E12.5)') (PDCS(IE,K),K=1,NA)
C  ****  Consistency test.
        DO K=1,NA
          DCSI(K)=PDCS(IE,K)
        ENDDO
        ECS0=4.0D0*PI*RMOMX(XMU,DCSI,0.0D0,1.0D0,NA,0)
        ECS1=4.0D0*PI*RMOMX(XMU,DCSI,0.0D0,1.0D0,NA,1)
        ECS2=4.0D0*PI*RMOMX(XMU,DCSI,0.0D0,1.0D0,NA,2)
        TCS1=2.0D0*ECS1
        TCS2=6.0D0*(ECS1-ECS2)
        TS0=(ECS0-PCS(IE))/PCS(IE)
        TS1=(TCS1-PTCS1(IE))/PTCS1(IE)
        TS2=(TCS2-PTCS2(IE))/PTCS2(IE)
        TSTE=MAX(ABS(TS0),ABS(TS1),ABS(TS2))
        IF(TSTE.GT.1.0D-4) THEN
          WRITE(IWR,'('' E='',1P,E12.5)') ETS(IE)
          WRITE(6,'('' E='',1P,E12.5)') ETS(IE)
          WRITE(IWR,*) ' Positron cross section data are corrupt.'
          STOP ' Positron cross section data are corrupt.'
        ENDIF
      ENDDO
C
C  ************  Electrons.
C
      IEME=0
      DO KE=1,NEGP
        IF(ET(KE).GT.0.999999D8) GO TO 100
        CALL DCSEL0(ET(KE),-1)
        CALL QRNDI0(DCSEL,0.0D0,1.0D0,NP,ERRM)
C       WRITE(6,'(''  E, ERRM ='',1P,2E13.6)') ET(KE),ERRM
        DO I=1,NP
          XSE(I,KE,M)=XTI(I)
          PSE(I,KE,M)=PACI(I)
          ASE(I,KE,M)=AI(I)
          BSE(I,KE,M)=BI(I)
          ITLE(I,KE,M)=ITTLI(I)
          ITUE(I,KE,M)=ITTUI(I)
        ENDDO
        CALL QRNDMS(0.0D0,1.0D0,XM0A,XM1,XM2)
        ECS0=CSI
        ECS1=CSI*XM1/XM0A
        ECS2=CSI*XM2/XM0A
        XE0(KE)=ECS0
        XE1(KE)=2.0D0*ECS1
        XE2(KE)=6.0D0*(ECS1-ECS2)
C
        FPEL=1.0D0/(XE0(KE)*VMOL(M))
        FPT1=1.0D0/(XE1(KE)*VMOL(M))
        FPST=ET(KE)/(CSTPE(M,KE)+RSTPE(M,KE))
        XS0H=1.0D0/(VMOL(M)*MAX(FPEL,MIN(C1(M)*FPT1,C2(M)*FPST)))
        RNDC=MAX(1.0D0-XS0H/XE0(KE),1.0D-10)
        IF(RNDC.LT.1.0D-6) RNDC=0.0D0
        RNDCE(M,KE)=RNDC
C
        RU=RNDC
        I=1
        J=NP
    1   K=(I+J)/2
        IF(RU.GT.PSE(K,KE,M)) THEN
          I=K
        ELSE
          J=K
        ENDIF
        IF(J-I.GT.1) GO TO 1
C
        RR=RU-PSE(I,KE,M)
        DPRO=PSE(I+1,KE,M)-PSE(I,KE,M)
        IF(DPRO.LT.1.0D-10) THEN
          RMUC=XSE(I,KE,M)
        ELSE
          CI=(1.0D0+ASE(I,KE,M)+BSE(I,KE,M))*DPRO
          RMUC=XSE(I,KE,M)+(CI*RR/(DPRO**2+(DPRO*ASE(I,KE,M)
     1        +BSE(I,KE,M)*RR)*RR))*(XSE(I+1,KE,M)-XSE(I,KE,M))
        ENDIF
C
C  ****  Moments of the PDF on the restricted interval (0,RMUC).
C        Total and transport cross sections for soft interactions.
C
        CALL QRNDMS(0.0D0,RMUC,XM0,XM1,XM2)
        ECS1=CSI*XM1/XM0A
        ECS2=CSI*XM2/XM0A
        TCS1=2.0D0*ECS1
        TCS2=6.0D0*(ECS1-ECS2)
        SEHEL(M,KE)=XS0H*VMOL(M)
        T1E0(KE)=TCS1
        T1E(M,KE)=T1EI(KE)+TCS1*VMOL(M)
        T2E0(KE)=TCS2
        T2E(M,KE)=T2EI(KE)+TCS2*VMOL(M)
        IEME=KE
      ENDDO
  100 CONTINUE
      EELMAX(M)=MIN(ET(IEME)-1.0D0,0.999999D8)
C
C  ****  Print electron elastic scattering tables.
C
      IF(INFO.GE.3) WRITE(IWR,1002)
 1002 FORMAT(/1X,'PENELOPE >>>  Elastic scattering of electrons',
     1  ' (ELSEP database)')
      IF(INFO.GE.3) WRITE(IWR,1003)
 1003 FORMAT(/3X,'E (eV)',6X,'MFP (mtu)',3X,'TMFP1 (mtu)',2X,
     1  'MFPh (mtu)',/1X,50('-'))
      DO I=1,IEME
        FP0=RHO(M)/(XE0(I)*VMOL(M))
        FP1=RHO(M)/(XE1(I)*VMOL(M))
        HMFP=RHO(M)/SEHEL(M,I)
        IF(INFO.GE.3) WRITE(IWR,'(1P,7(E12.5,1X))') ET(I),FP0,FP1,
     1    HMFP
        SEHEL(M,I)=LOG(SEHEL(M,I))
C  ****  Soft scattering events are switched off when T1E is too small.
        IF(T1E(M,I).GT.1.0D-6*XE1(I)*VMOL(M)) THEN
          T1E(M,I)=LOG(MAX(T1E(M,I),1.0D-35))
          T2E(M,I)=LOG(MAX(T2E(M,I),1.0D-35))
        ELSE
          T1E(M,I)=-100.0D0
          T2E(M,I)=-100.0D0
        ENDIF
      ENDDO
C
C  ************  Positrons.
C
      IEMP=0
      DO KE=1,NEGP
        IF(ET(KE).GT.0.999999D8) GO TO 200
        CALL DCSEL0(ET(KE),+1)
        CALL QRNDI0(DCSEL,0.0D0,1.0D0,NP,ERRM)
        DO I=1,NP
          XSP(I,KE,M)=XTI(I)
          PSP(I,KE,M)=PACI(I)
          ASP(I,KE,M)=AI(I)
          BSP(I,KE,M)=BI(I)
          ITLP(I,KE,M)=ITTLI(I)
          ITUP(I,KE,M)=ITTUI(I)
        ENDDO
        CALL QRNDMS(0.0D0,1.0D0,XM0A,XM1,XM2)
        ECS0=CSI
        ECS1=CSI*XM1/XM0A
        ECS2=CSI*XM2/XM0A
        XP0(KE)=ECS0
        XP1(KE)=2.0D0*ECS1
        XP2(KE)=6.0D0*(ECS1-ECS2)
C
        FPEL=1.0D0/(XP0(KE)*VMOL(M))
        FPT1=1.0D0/(XP1(KE)*VMOL(M))
        FPST=ET(KE)/(CSTPP(M,KE)+RSTPP(M,KE))
        XS0H=1.0D0/(VMOL(M)*MAX(FPEL,MIN(C1(M)*FPT1,C2(M)*FPST)))
        RNDC=MAX(1.0D0-XS0H/XP0(KE),1.0D-10)
        IF(RNDC.LT.1.0D-6) RNDC=0.0D0
        RNDCP(M,KE)=RNDC
C
        RU=RNDC
        I=1
        J=NP
    2   K=(I+J)/2
        IF(RU.GT.PSP(K,KE,M)) THEN
          I=K
        ELSE
          J=K
        ENDIF
        IF(J-I.GT.1) GO TO 2
C
        RR=RU-PSP(I,KE,M)
        DPRO=PSP(I+1,KE,M)-PSP(I,KE,M)
        IF(DPRO.LT.1.0D-10) THEN
          RMUC=XSP(I,KE,M)
        ELSE
          CI=(1.0D0+ASP(I,KE,M)+BSP(I,KE,M))*DPRO
          RMUC=XSP(I,KE,M)+(CI*RR/(DPRO**2+(DPRO*ASP(I,KE,M)
     1        +BSP(I,KE,M)*RR)*RR))*(XSP(I+1,KE,M)-XSP(I,KE,M))
        ENDIF
C
C  ****  Moments of the PDF on the restricted interval (0,RMUC).
C        Total and transport cross sections for soft interactions.
C
        CALL QRNDMS(0.0D0,RMUC,XM0,XM1,XM2)
        ECS1=CSI*XM1/XM0A
        ECS2=CSI*XM2/XM0A
        TCS1=2.0D0*ECS1
        TCS2=6.0D0*(ECS1-ECS2)
        SPHEL(M,KE)=XS0H*VMOL(M)
        T1P0(KE)=TCS1
        T1P(M,KE)=T1PI(KE)+TCS1*VMOL(M)
        T2P0(KE)=TCS2
        T2P(M,KE)=T2PI(KE)+TCS2*VMOL(M)
        IEMP=KE
      ENDDO
  200 CONTINUE
      PELMAX(M)=MIN(ET(IEMP)-1.0D0,0.999999D8)
C
C  ****  Print positron elastic scattering tables.
C
      IF(INFO.GE.3) WRITE(IWR,1004)
 1004 FORMAT(/1X,'PENELOPE >>>  Elastic scattering of positrons',
     1  ' (ELSEP database)')
      IF(INFO.GE.3) WRITE(IWR,1005)
 1005 FORMAT(/3X,'E (eV)',6X,'MFP (mtu)',3X,'TMFP1 (mtu)',2X,
     1  'MFPh (mtu)',/1X,50('-'))
      DO I=1,IEMP
        FP0=RHO(M)/(XP0(I)*VMOL(M))
        FP1=RHO(M)/(XP1(I)*VMOL(M))
        HMFP=RHO(M)/SPHEL(M,I)
        IF(INFO.GE.3) WRITE(IWR,'(1P,7(E12.5,1X))') ET(I),FP0,FP1,
     1    HMFP
        SPHEL(M,I)=LOG(SPHEL(M,I))
C  ****  Soft scattering events are switched off when T1P is too small.
        IF(T1P(M,I).GT.1.0D-6*XP1(I)*VMOL(M)) THEN
          T1P(M,I)=LOG(MAX(T1P(M,I),1.0D-35))
          T2P(M,I)=LOG(MAX(T2P(M,I),1.0D-35))
        ELSE
          T1P(M,I)=-100.0D0
          T2P(M,I)=-100.0D0
        ENDIF
      ENDDO
      RETURN
      END
C  *********************************************************************
C                       SUBROUTINE EELd
C  *********************************************************************
      SUBROUTINE EELd(RNDC,RMU)
C
C     Simulation of electron hard elastic events. Cross sections from
C  the ELSEP numerical database.
C
C  Argument value:
C    RNDC ... cutoff value of the uniform random number
C             (only hard events are simulated).
C    RMU .... sampled angular deflection, =(1-CDT)/2.
C
      use mtpar
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N)
C  ****  Main-PENELOPE common.
      COMMON/TRACK/E,X,Y,Z,U,V,W,WGHT,KPAR,IBODY,M,ILB(5)
C  ****  Energy grid and interpolation constants for the current energy.
      PARAMETER (NEGP=200)
      COMMON/CEGRID/EL,EU,ET(NEGP),DLEMP(NEGP),DLEMP1,DLFC,
     1  XEL,XE,XEK,KE
C  ****  Electron simulation tables.
      PARAMETER (MAXMAT=10)
      PARAMETER (NP=128,NPM1=NP-1)
      COMMON/CEELDB/XSE(NP,NEGP,MAXMAT),PSE(NP,NEGP,MAXMAT),
     1              ASE(NP,NEGP,MAXMAT),BSE(NP,NEGP,MAXMAT),
     2              ITLE(NP,NEGP,MAXMAT),ITUE(NP,NEGP,MAXMAT)
C
      ! EXTERNALRAND
C  ****  Energy grid point.
      PK=(XEL-DLEMP(KE))*DLFC
        status=rmt(r)
      IF(r.LT.PK) THEN
        JE=KE+1
      ELSE
        JE=KE
      ENDIF
C  ****  Pointer.
        status=rmt(r)
      RU=RNDC+r*(1.0D0-RNDC)
C  ****  Selection of the interval (binary search in a restricted
C        interval).
      ITN=RU*NPM1+1
      I=ITLE(ITN,JE,M)
      J=ITUE(ITN,JE,M)
      IF(J-I.LT.2) GO TO 2
    1 K=(I+J)/2
      IF(RU.GT.PSE(K,JE,M)) THEN
        I=K
      ELSE
        J=K
      ENDIF
      IF(J-I.GT.1) GO TO 1
C  ****  Sampling from the rational inverse cumulative distribution.
    2 CONTINUE
      PP=PSE(I,JE,M)
      RR=RU-PP
      IF(RR.GT.1.0D-16) THEN
        XX=XSE(I,JE,M)
        AA=ASE(I,JE,M)
        BB=BSE(I,JE,M)
        D=PSE(I+1,JE,M)-PP
        CD=(1.0D0+AA+BB)*D
        RMU=XX+(CD*RR/(D*D+(AA*D+BB*RR)*RR))*(XSE(I+1,JE,M)-XX)
      ELSE
        RMU=XSE(I,JE,M)
      ENDIF
      RETURN
      END
C  *********************************************************************
C                       SUBROUTINE PELd
C  *********************************************************************
      SUBROUTINE PELd(RNDC,RMU)
C
C     Simulation of positron hard elastic events. Cross sections from
C  the ELSEP numerical database.
C
C  Argument value:
C    RNDC ... cutoff value of the uniform random number
C             (only hard events are simulated).
C    RMU .... sampled angular deflection, =(1-CDT)/2.
C
      use mtpar
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N)
C  ****  Main-PENELOPE common.
      COMMON/TRACK/E,X,Y,Z,U,V,W,WGHT,KPAR,IBODY,M,ILB(5)
C  ****  Energy grid and interpolation constants for the current energy.
      PARAMETER (NEGP=200)
      COMMON/CEGRID/EL,EU,ET(NEGP),DLEMP(NEGP),DLEMP1,DLFC,
     1  XEL,XE,XEK,KE
C  ****  Positron simulation tables.
      PARAMETER (MAXMAT=10)
      PARAMETER (NP=128,NPM1=NP-1)
      COMMON/CPELDB/XSP(NP,NEGP,MAXMAT),PSP(NP,NEGP,MAXMAT),
     1              ASP(NP,NEGP,MAXMAT),BSP(NP,NEGP,MAXMAT),
     2              ITLP(NP,NEGP,MAXMAT),ITUP(NP,NEGP,MAXMAT)
C
      ! EXTERNALRAND
C  ****  Energy grid point.
      PK=(XEL-DLEMP(KE))*DLFC
        status=rmt(r)
      IF(r.LT.PK) THEN
        JE=KE+1
      ELSE
        JE=KE
      ENDIF
C  ****  Pointer.
        status=rmt(r)
      RU=RNDC+r*(1.0D0-RNDC)
C  ****  Selection of the interval (binary search in a restricted
C        interval).
      ITN=RU*NPM1+1
      I=ITLP(ITN,JE,M)
      J=ITUP(ITN,JE,M)
      IF(J-I.LT.2) GO TO 2
    1 K=(I+J)/2
      IF(RU.GT.PSP(K,JE,M)) THEN
        I=K
      ELSE
        J=K
      ENDIF
      IF(J-I.GT.1) GO TO 1
C  ****  Sampling from the rational inverse cumulative distribution.
    2 CONTINUE
      PP=PSP(I,JE,M)
      RR=RU-PP
      IF(RR.GT.1.0D-16) THEN
        XX=XSP(I,JE,M)
        AA=ASP(I,JE,M)
        BB=BSP(I,JE,M)
        D=PSP(I+1,JE,M)-PP
        CD=(1.0D0+AA+BB)*D
        RMU=XX+(CD*RR/(D*D+(AA*D+BB*RR)*RR))*(XSP(I+1,JE,M)-XX)
      ELSE
        RMU=XSP(I,JE,M)
      ENDIF
      RETURN
      END
C  *********************************************************************
C                 SUBROUTINES ELINIT, DCSEL0 AND DCSEL
C  *********************************************************************
C
C     These subroutines read the ELSEP database for elastic scattering
C  of electrons and positrons by neutral atoms, and generate the molecu-
C  lar DCS of a compound for arbitrary values of the energy (from 50 eV
C  to 100 MeV) and the angular deflection RMU=(1-COS(THETA))/2.
C
C  Other subprograms needed: subroutine SPLINE,
C                            function FINDI.
C
C  --> Subroutine ELINIT reads atomic elastic DCSs for electrons and
C  positrons from the database files and determines a table of the mole-
C  cular DCS, for the standard grid of energies and angular deflections,
C  as the incoherent sum of atomic DCSs. It is assumed that the database
C  files are in the same directory as the binary executable module. If
C  you wish to keep the database files on a separate directory, you can
C  edit the present source file and change the string FILE1, which con-
C  tains the names of the database files of the element, to include the
C  directory path.
C
C  --> Subroutine DCSEL0(E,IELEC) initializes the calculation of DCSs
C  for electrons (IELEC=-1) or positrons (IELEC=+1) with kinetic energy
C  E (eV). It builds a table of DCS values for the standard grid of
C  angular deflections RMU using log-log cubic spline interpolation in E
C  of the tables prepared by subroutine ELINIT.
C
C  --> Function DCSEL(RMU) computes the DCS in (cm**2/sr) at RMU by
C  linear log-log interpolation of the values tabulated by subroutine
C  DCSEL0. Notice that the delivered DCSEL value corresponds to the
C  particle and energy defined in the last call to subroutine DCSEL0.
C
C  EXAMPLE: To calculate cross sections for water, our main program must
C  contain the following definitions and calls:
C
C     CALL ELINIT(IZ,STF,NELEM) with NELEM=2
C                                    IZ(1)=1, STF(1)=2   <-- 2 H atoms
C                                    IZ(2)=8, STF(2)=1   <-- 1 O atom
C  (This sets the standard tabulation of cross sections for water.)
C
C  Now, to obtain the DCS for electrons with a kinetic energy of 10 keV
C  at RMU=0.0D0 (forward scattering) we simply insert the following two
C  statements in the main program,
C
C     CALL DCSEL0(1.0D4,-1)
C     DCS=DCSEL(0.0D0)
C
C  *********************************************************************
C                       SUBROUTINE ELINIT
C  *********************************************************************
      SUBROUTINE ELINIT(IZ,STF,NELEM)
C
C  This subroutine reads atomic elastic cross sections for electrons and
C  positrons from the database files and determines the molecular cross
C  section as the incoherent sum of atomic cross sections.
C
C  Input arguments:
C    IZ (1:NELEM) .... atomic numbers of the elements in the compound.
C    STF (1:NELEM) ... stoichiometric indices.
C    NELEM ........... number of different elements.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N)
      PARAMETER (PI=3.1415926535897932D0)
C
      DIMENSION IZ(NELEM),STF(NELEM)
C
      PARAMETER (NE=96,NA=606)
      COMMON/CDCSEP/ET(NE),ETL(NE),TH(NA),THR(NA),XMU(NA),XMUL(NA),
     1       ECS(NE),ETCS1(NE),ETCS2(NE),EDCS(NE,NA),
     2       PCS(NE),PTCS1(NE),PTCS2(NE),PDCS(NE,NA),
     3       DCSI(NA),DCSIL(NA),CSI,TCS1I,TCS2I
C
      DIMENSION EGRD(16)
      DATA EGRD/1.0D0,1.25D0,1.50D0,1.75D0,2.00D0,2.50D0,3.00D0,
     1  3.50D0,4.00D0,4.50D0,5.00D0,6.00D0,7.00D0,8.00D0,9.00D0,
     2  1.00D1/
C
      CHARACTER*1 LIT10(10),LIT1,LIT2,LIT3
      DATA LIT10/'0','1','2','3','4','5','6','7','8','9'/
      CHARACTER*12 FILE1
C
C  ****  Energy mesh points (in eV).
C
      IE=0
      IGRID=10
      FGRID=10.0D0
   10 IGRID=IGRID+1
      EV=EGRD(IGRID)*FGRID
      IF(IGRID.EQ.16) THEN
        IGRID=1
        FGRID=10.0D0*FGRID
      ENDIF
      IE=IE+1
      ET(IE)=EV
      ETL(IE)=LOG(ET(IE))
C     WRITE(6,'(I5,1P,3E15.7)') IE,ET(IE)
      IF(IE.LT.NE) GO TO 10
C
C  ****  Angular grid (TH in deg, XMU=(1.0D0-COS(TH))/2).
C
      I=1
      TH(I)=0.0D0
      THR(I)=TH(I)*PI/180.0D0
      XMU(I)=(1.0D0-COS(THR(I)))/2.0D0
      XMUL(I)=LOG(1.0D-35)
      I=2
      TH(I)=1.0D-4
      THR(I)=TH(I)*PI/180.0D0
      XMU(I)=(1.0D0-COS(THR(I)))/2.0D0
      XMUL(I)=LOG(MAX(XMU(I),1.0D-35))
   20 CONTINUE
      I=I+1
      IF(TH(I-1).LT.0.9999D-3) THEN
        TH(I)=TH(I-1)+2.5D-5
      ELSE IF(TH(I-1).LT.0.9999D-2) THEN
        TH(I)=TH(I-1)+2.5D-4
      ELSE IF(TH(I-1).LT.0.9999D-1) THEN
        TH(I)=TH(I-1)+2.5D-3
      ELSE IF(TH(I-1).LT.0.9999D+0) THEN
        TH(I)=TH(I-1)+2.5D-2
      ELSE IF(TH(I-1).LT.0.9999D+1) THEN
        TH(I)=TH(I-1)+1.0D-1
      ELSE IF(TH(I-1).LT.2.4999D+1) THEN
        TH(I)=TH(I-1)+2.5D-1
      ELSE
        TH(I)=TH(I-1)+5.0D-1
      ENDIF
      THR(I)=TH(I)*PI/180.0D0
      XMU(I)=MAX((1.0D0-COS(THR(I)))/2.0D0,1.0D-35)
      XMUL(I)=LOG(MAX(XMU(I),1.0D-35))
      IF(I.LT.NA) GO TO 20
C
C     DO I=1,NA
C       WRITE(6,'(I5,1P,3E15.7)') I,TH(I),THR(I),XMU(I)
C     ENDDO
C
      DO IE=1,NE
        ECS(IE)=0.0D0
        ETCS1(IE)=0.0D0
        ETCS2(IE)=0.0D0
        PCS(IE)=0.0D0
        PTCS1(IE)=0.0D0
        PTCS2(IE)=0.0D0
        DO IA=1,NA
          EDCS(IE,IA)=0.0D0
          PDCS(IE,IA)=0.0D0
        ENDDO
      ENDDO
C
C  ****  Read atomic DCS tables and compute the molecular DCS as the
C        incoherent sum of atomic DCSs.
C
      DO IEL=1,NELEM
        IZZ=IZ(IEL)
        STFF=STF(IEL)
        NS=IZ(IEL)
        IF(NS.GT.89999) NS=89999
        NS1=NS-10*(NS/10)
        NS=(NS-NS1)/10
        NS2=NS-10*(NS/10)
        NS=(NS-NS2)/10
        NS3=NS-10*(NS/10)
        LIT1=LIT10(NS1+1)
        LIT2=LIT10(NS2+1)
        LIT3=LIT10(NS3+1)
C
        FILE1='eeldx'//LIT3//LIT2//LIT1//'.p05'
        OPEN(3,FILE=FILE1)
        DO IE=1,NE
          READ(3,'(I3,I4,1P,E10.3,5E12.5)')
     1      IELEC,IZR,ENR,CSE,TCS1E,TCS2E
          IF(IELEC.NE.-1.OR.IZR.NE.IZZ.OR.ABS(ENR-ET(IE)).GT.1.0D-3)
     1      STOP 'Corrupt data file.'
          READ(3,'(1P,10E12.5)') (DCSI(IA),IA=1,NA)
          ECS(IE)=ECS(IE)+STFF*CSE
          ETCS1(IE)=ETCS1(IE)+STFF*TCS1E
          ETCS2(IE)=ETCS2(IE)+STFF*TCS2E
          DO IA=1,NA
            EDCS(IE,IA)=EDCS(IE,IA)+STFF*DCSI(IA)
          ENDDO
        ENDDO
        CLOSE(3)
C
        FILE1='peldx'//LIT3//LIT2//LIT1//'.p05'
        OPEN(3,FILE=FILE1)
        DO IE=1,NE
          READ(3,'(I3,I4,1P,E10.3,5E12.5)')
     1      IELEC,IZR,ENR,CSP,TCS1P,TCS2P
          IF(IELEC.NE.+1.OR.IZR.NE.IZZ.OR.ABS(ENR-ET(IE)).GT.1.0D-3)
     1      STOP 'Corrupt data file.'
          READ(3,'(1P,10E12.5)') (DCSI(IA),IA=1,NA)
          PCS(IE)=PCS(IE)+STFF*CSP
          PTCS1(IE)=PTCS1(IE)+STFF*TCS1P
          PTCS2(IE)=PTCS2(IE)+STFF*TCS2P
          DO IA=1,NA
            PDCS(IE,IA)=PDCS(IE,IA)+STFF*DCSI(IA)
          ENDDO
        ENDDO
        CLOSE(3)
      ENDDO
      RETURN
      END
C  *********************************************************************
C                       SUBROUTINE DCSEL0
C  *********************************************************************
      SUBROUTINE DCSEL0(E,IELEC)
C
C     This subroutine computes a table of the molecular elastic DCSs for
C  electrons (IELEC=-1) or positrons (IELEC=+1) with kinetic energy  E
C  (in eV) by log-log cubic spline interpolation in E of the DCS table
C  prepared by subroutine ELINIT.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N)
C
      PARAMETER (NE=96,NA=606)
      COMMON/CDCSEP/ET(NE),ETL(NE),TH(NA),THR(NA),XMU(NA),XMUL(NA),
     1       ECS(NE),ETCS1(NE),ETCS2(NE),EDCS(NE,NA),
     2       PCS(NE),PTCS1(NE),PTCS2(NE),PDCS(NE,NA),
     3       DCSI(NA),DCSIL(NA),CSI,TCS1I,TCS2I
C
      DIMENSION Y(NE),A(NE),B(NE),C(NE),D(NE)
C
      IF(E.LT.49.999D0.OR.E.GT.1.0001D8) THEN
        WRITE(6,'('' Error in DCSEL0: Energy out of range.'')')
        STOP 'Energy out of range.'
      ENDIF
C
      EL=LOG(E)
      CALL FINDI(ETL,EL,NE,JE)
C
      DO IA=1,NA
        DO IE=1,NE
          IF(IELEC.EQ.-1) THEN
            Y(IE)=LOG(EDCS(IE,IA))
          ELSE
            Y(IE)=LOG(PDCS(IE,IA))
          ENDIF
        ENDDO
        CALL SPLINE(ETL,Y,A,B,C,D,0.0D0,0.0D0,NE)
        DCSIL(IA)=A(JE)+EL*(B(JE)+EL*(C(JE)+EL*D(JE)))
        DCSI(IA)=EXP(DCSIL(IA))
      ENDDO
C
      DO IE=1,NE
        IF(IELEC.EQ.-1) THEN
          Y(IE)=LOG(ECS(IE))
        ELSE
          Y(IE)=LOG(PCS(IE))
        ENDIF
      ENDDO
      CALL SPLINE(ETL,Y,A,B,C,D,0.0D0,0.0D0,NE)
      CSI=EXP(A(JE)+EL*(B(JE)+EL*(C(JE)+EL*D(JE))))
C
      DO IE=1,NE
        IF(IELEC.EQ.-1) THEN
          Y(IE)=LOG(ETCS1(IE))
        ELSE
          Y(IE)=LOG(PTCS1(IE))
        ENDIF
      ENDDO
      CALL SPLINE(ETL,Y,A,B,C,D,0.0D0,0.0D0,NE)
      TCS1I=EXP(A(JE)+EL*(B(JE)+EL*(C(JE)+EL*D(JE))))
C
      DO IE=1,NE
        IF(IELEC.EQ.-1) THEN
          Y(IE)=LOG(ETCS2(IE))
        ELSE
          Y(IE)=LOG(PTCS2(IE))
        ENDIF
      ENDDO
      CALL SPLINE(ETL,Y,A,B,C,D,0.0D0,0.0D0,NE)
      TCS2I=EXP(A(JE)+EL*(B(JE)+EL*(C(JE)+EL*D(JE))))
C
      RETURN
      END
C  *********************************************************************
C                       FUNCTION DCSEL
C  *********************************************************************
      FUNCTION DCSEL(RMU)
C
C  This function computes the DCS in (cm**2/sr) by linear log-log inter-
C  polation in RMU=(1-cos(theta))/2. It is initialized by subroutine
C  DCSEL0, which must be called before using DCSEL.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N)
C
      PARAMETER (NE=96,NA=606)
      COMMON/CDCSEP/ET(NE),ETL(NE),TH(NA),THR(NA),XMU(NA),XMUL(NA),
     1       ECS(NE),ETCS1(NE),ETCS2(NE),EDCS(NE,NA),
     2       PCS(NE),PTCS1(NE),PTCS2(NE),PDCS(NE,NA),
     3       DCSI(NA),DCSIL(NA),CSI,TCS1I,TCS2I
C
      RMUL=LOG(MAX(RMU,1.0D-35))
      CALL FINDI(XMUL,RMUL,NA,I)
      DCSEL=EXP(DCSIL(I)+(DCSIL(I+1)-DCSIL(I))
     1     *((RMUL-XMUL(I))/(XMUL(I+1)-XMUL(I))))
      RETURN
      END
C  *********************************************************************
C                       FUNCTION RMOMX
C  *********************************************************************
      FUNCTION RMOMX(X,PDF,XD,XU,NP,MOM)
C
C     Calculation of momenta of a pdf, PDF(X), obtained from linear
C  log-log interpolation of the input table. The independent variable X
C  is assumed to take only positive values.
C
C     X ........ array of variable values (in increasing order).
C     PDF ...... corresponding PDF values (must be non-negative).
C     NP ....... number of points in the table.
C     XD, XU ... limits of the integration interval.
C     MOM ...... moment order.
C     RMOM = INTEGRAL (X**N)*PDF(X) dX over the interval (XD,XU)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N)
      PARAMETER (EPS=1.0D-12,ZERO=1.0D-35)
      DIMENSION X(NP),PDF(NP)
C
      IF(NP.LT.2) STOP 'RMOMX. Error code 1.'
      IF(X(1).LT.0.0D0.OR.PDF(1).LT.0.0D0) THEN
        WRITE(6,*) 'X(1),PDF(1) =',X(1),PDF(1)
        STOP 'RMOMX. Error code 2'
      ENDIF
      DO I=2,NP
        IF(X(I).LT.0.0D0.OR.PDF(I).LT.0.0D0) THEN
          WRITE(6,*) 'I,X(I),PDF(I) =',I,X(I),PDF(I)
          STOP 'RMOMX. Error code 3'
        ENDIF
        IF(X(I).LT.X(I-1)) STOP 'RMOMX. Error code 4.'
      ENDDO
C
      XLOW=MAX(X(1),XD)
      IF(XLOW.LT.ZERO) XLOW=ZERO
      XUP=MIN(X(NP),XU)
C
      IF(XLOW.GE.XUP) THEN
        WRITE(6,*) ' WARNING: XLOW is greater than XUP in RMOMX.'
        WRITE(6,*) ' XLOW =',XLOW,',   XUP =',XUP
        RMOMX=0.0D0
        RETURN
      ENDIF
C
      IL=1
      IU=NP-1
      DO I=1,NP-1
        IF(X(I).LT.XLOW) IL=I
        IF(X(I).LT.XUP) IU=I
      ENDDO
C
C  ****  A single interval.
C
      IF(IU.EQ.IL) THEN
        XIL=LOG(MAX(X(IL),ZERO))
        XFL=LOG(X(IL+1))
        YIL=LOG(MAX(PDF(IL),ZERO))
        YFL=LOG(MAX(PDF(IL+1),ZERO))
        X1=XLOW
        X2=XUP
        DEN=XFL-XIL
        IF(ABS(DEN).GT.ZERO) THEN
          Y1=EXP(YIL+(YFL-YIL)*(LOG(X1)-XIL)/DEN)*X1**MOM
          Y2=EXP(YIL+(YFL-YIL)*(LOG(X2)-XIL)/DEN)*X2**MOM
        ELSE
          Y1=EXP(YIL)*X1**MOM
          Y2=EXP(YIL)*X2**MOM
        ENDIF
        DXL=LOG(X2)-LOG(X1)
        DYL=LOG(MAX(Y2,ZERO))-LOG(MAX(Y1,ZERO))
        IF(ABS(DXL).GT.EPS*ABS(DYL)) THEN
          AP1=1.0D0+(DYL/DXL)
          IF(ABS(AP1).GT.EPS) THEN
            DSUM=(Y2*X2-Y1*X1)/AP1
          ELSE
            DSUM=Y1*X1*DXL
          ENDIF
        ELSE
          DSUM=0.5D0*(Y1+Y2)*(X2-X1)
        ENDIF
        RMOMX=DSUM
        RETURN
      ENDIF
C
C  ****  Multiple intervals.
C
      XIL=LOG(MAX(X(IL),ZERO))
      XFL=LOG(X(IL+1))
      YIL=LOG(MAX(PDF(IL),ZERO))
      YFL=LOG(MAX(PDF(IL+1),ZERO))
      X1=XLOW
      DEN=XFL-XIL
      IF(ABS(DEN).GT.ZERO) THEN
        Y1=EXP(YIL+(YFL-YIL)*(LOG(X1)-XIL)/DEN)*X1**MOM
      ELSE
        Y1=EXP(YIL)*X1**MOM
      ENDIF
      X2=X(IL+1)
      Y2=MAX(PDF(IL+1),ZERO)*X2**MOM
      DXL=LOG(X2)-LOG(X1)
      DYL=LOG(MAX(Y2,ZERO))-LOG(MAX(Y1,ZERO))
      IF(ABS(DXL).GT.EPS*ABS(DYL)) THEN
        AP1=1.0D0+(DYL/DXL)
        IF(ABS(AP1).GT.EPS) THEN
          DSUM=(Y2*X2-Y1*X1)/AP1
        ELSE
          DSUM=Y1*X1*DXL
        ENDIF
      ELSE
        DSUM=0.5D0*(Y1+Y2)*(X2-X1)
      ENDIF
      RMOMX=DSUM
C
      IF(IU.GT.IL+1) THEN
        DO I=IL+1,IU-1
          X1=X(I)
          Y1=MAX(PDF(I),ZERO)*X1**MOM
          X2=X(I+1)
          Y2=MAX(PDF(I+1),ZERO)*X2**MOM
          DXL=LOG(X2)-LOG(X1)
          DYL=LOG(MAX(Y2,ZERO))-LOG(MAX(Y1,ZERO))
          IF(ABS(DXL).GT.EPS*ABS(DYL)) THEN
            AP1=1.0D0+(DYL/DXL)
            IF(ABS(AP1).GT.EPS) THEN
              DSUM=(Y2*X2-Y1*X1)/AP1
            ELSE
              DSUM=Y1*X1*DXL
            ENDIF
          ELSE
            DSUM=0.5D0*(Y1+Y2)*(X2-X1)
          ENDIF
          RMOMX=RMOMX+DSUM
        ENDDO
      ENDIF
C
      X1=X(IU)
      Y1=MAX(PDF(IU),ZERO)*X1**MOM
      XIL=LOG(X(IU))
      XFL=LOG(X(IU+1))
      YIL=LOG(MAX(PDF(IU),ZERO))
      YFL=LOG(MAX(PDF(IU+1),ZERO))
      X2=XUP
      DEN=XFL-XIL
      IF(ABS(DEN).GT.ZERO) THEN
        Y2=EXP(YIL+(YFL-YIL)*(LOG(X2)-XIL)/DEN)*X2**MOM
      ELSE
        Y2=EXP(YIL)*X2**MOM
      ENDIF
      DXL=LOG(X2)-LOG(X1)
      DYL=LOG(MAX(Y2,ZERO))-LOG(MAX(Y1,ZERO))
      IF(ABS(DXL).GT.EPS*ABS(DYL)) THEN
        AP1=1.0D0+(DYL/DXL)
        IF(ABS(AP1).GT.EPS) THEN
          DSUM=(Y2*X2-Y1*X1)/AP1
        ELSE
          DSUM=Y1*X1*DXL
        ENDIF
      ELSE
        DSUM=0.5D0*(Y1+Y2)*(X2-X1)
      ENDIF
      RMOMX=RMOMX+DSUM
C
      RETURN
      END


C  *********************************************************************
C                        SUBROUTINE QRNDI0
C  *********************************************************************
      SUBROUTINE QRNDI0(PDF,XMIN,XMAX,N,ERRM)
C
C     Initialization of the QRNDI algorithm for random sampling of a
C  continuous random variable X from a probability distribution function
C  PDF(X) defined in the domain (XMIN,XMAX).
C
C  Other subprograms needed: subroutine SIMPSU.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N)
      PARAMETER (EPS=1.0D-10, ZERO=1.0D-34, ZEROT=0.1D0*ZERO)
      PARAMETER (NM=512)
C
C     The information used by the sampling function QRND is exported
C  through the common block /CQRND/,
      COMMON/CQRND/X(NM),PAC(NM),A(NM),B(NM),NP
C  where
C    X(I) ...... grid points, in increasing order.
C    PAC(I) .... value of the cumulative pdf at X(I).
C    A(I), B(I) ... rational inverse cumulative distribution parameters.
C    NP ........ number of grid points (.LE.NM).
C
C  The sampling function QRNDI finds the X-interval using binary search
C  within precalculated limits. It uses the parameters in the common
C  blocks /CQRND/ and /CQRND1/,
      COMMON/CQRND1/ITTL(NM),ITTU(NM),NPM1
C  where
C    ITTL(I) ... largest J for which PAC(J)<(I-1)/(NP-1),
C    ITTU(I) ... smallest K for which PAC(K)>I/(NP-1),
C    NPM1 ...... =NP-1.
C
      PARAMETER (NIP=51)
      DIMENSION AREA(NM),ERR(NM),C(NM),XI(NIP),PDFI(NIP),SUMI(NIP)
      EXTERNAL PDF
C
      NUNIF=8
      IF(N.LE.NUNIF) THEN
        WRITE(6,'('' Error in QRNDI0: N must be larger than NUNIF ='',
     1    I3)') NUNIF
        STOP 'QRNDI0: N must be larger than NUNIF.'
      ENDIF
      IF(N.GT.NM) THEN
        WRITE(6,'('' Error in QRNDI0: N must be less than NM=512.'')')
        STOP 'QRNDI0: N must be less than NM=512.'
      ENDIF
      IF(XMIN.GE.XMAX-EPS) THEN
        WRITE(6,'('' Error in QRNDI0: XMIN must be larger than XMAX.'')
     1    ')
        STOP 'QRNDI0: XMIN must be larger than XMAX.'
      ENDIF
C
C  ****  We start with a grid of NUNIF points uniformly spaced in the
C        interval (XMIN,XMAX).
C
      NP=NUNIF
      DX=(XMAX-XMIN)/DBLE(NP-1)
      X(1)=XMIN
      DO I=1,NP-1
        X(I+1)=XMIN+I*DX
        DXI=(X(I+1)-X(I))/DBLE(NIP-1)
        PDFMAX=0.0D0
        DO K=1,NIP
          XI(K)=X(I)+DBLE(K-1)*DXI
          PDFI(K)=MAX(PDF(XI(K)),ZEROT)
          PDFMAX=MAX(PDFMAX,PDFI(K))
        ENDDO
        CALL SIMPSU(DXI,PDFI,SUMI,NIP)
        AREA(I)=SUMI(NIP)
        FACT=1.0D0/AREA(I)
        DO K=1,NIP
          SUMI(K)=FACT*SUMI(K)
        ENDDO
C  ****  When the PDF vanishes at one of the interval end points, it is
C        slightly modified.
        IF(PDFI(1).LT.ZERO) PDFI(1)=1.0D-5*PDFMAX
        IF(PDFI(NIP).LT.ZERO) PDFI(NIP)=1.0D-5*PDFMAX
C
        PLI=PDFI(1)*FACT
        PUI=PDFI(NIP)*FACT
        B(I)=1.0D0-1.0D0/(PLI*PUI*DX*DX)
        A(I)=(1.0D0/(PLI*DX))-1.0D0-B(I)
        C(I)=1.0D0+A(I)+B(I)
C
C  ****  ERR(I) is defined as the integral of the absolute difference
C        between the rational interpolation and the true PDF, extended
C        over the interval (X(I),X(I+1)).
C
        ICASE=1
  100   CONTINUE
        ERR(I)=0.0D0
        DO K=1,NIP
          RR=SUMI(K)
          PAP=AREA(I)*(1.0D0+(A(I)+B(I)*RR)*RR)**2/
     1       ((1.0D0-B(I)*RR*RR)*C(I)*(X(I+1)-X(I)))
          IF(K.EQ.1.OR.K.EQ.NIP) THEN
            ERR(I)=ERR(I)+0.5D0*ABS(PAP-PDFI(K))
          ELSE
            ERR(I)=ERR(I)+ABS(PAP-PDFI(K))
          ENDIF
        ENDDO
        ERR(I)=ERR(I)*DXI
C  ****  If ERR(I) is too large, the PDF is approximated by a uniform
C        distribution.
        IF(ERR(I).GT.0.10D0*AREA(I).AND.ICASE.EQ.1) THEN
          B(I)=0.0D0
          A(I)=0.0D0
          C(I)=1.0D0
          ICASE=2
          GO TO 100
        ENDIF
      ENDDO
      X(NP)=XMAX
C
C  ****  New grid points are added by halving the subinterval with the
C        largest absolute error.
C
  200 CONTINUE
      ERRM=0.0D0
      LMAX=1
      DO I=1,NP-1
C  ****  ERRM is the largest of the interval errors ERR(I).
        IF(ERR(I).GT.ERRM) THEN
          ERRM=ERR(I)
          LMAX=I
        ENDIF
      ENDDO
C
      NP=NP+1
      DO I=NP,LMAX+1,-1
        X(I)=X(I-1)
        A(I)=A(I-1)
        B(I)=B(I-1)
        C(I)=C(I-1)
        ERR(I)=ERR(I-1)
        AREA(I)=AREA(I-1)
      ENDDO
      X(LMAX+1)=0.5D0*(X(LMAX)+X(LMAX+2))
      DO I=LMAX,LMAX+1
        DX=X(I+1)-X(I)
        DXI=(X(I+1)-X(I))/DBLE(NIP-1)
        PDFMAX=0.0D0
        DO K=1,NIP
          XI(K)=X(I)+DBLE(K-1)*DXI
          PDFI(K)=MAX(PDF(XI(K)),ZEROT)
          PDFMAX=MAX(PDFMAX,PDFI(K))
        ENDDO
        CALL SIMPSU(DXI,PDFI,SUMI,NIP)
        AREA(I)=SUMI(NIP)
        FACT=1.0D0/AREA(I)
        DO K=1,NIP
          SUMI(K)=FACT*SUMI(K)
        ENDDO
C
        IF(PDFI(1).LT.ZERO) PDFI(1)=1.0D-5*PDFMAX
        IF(PDFI(NIP).LT.ZERO) PDFI(NIP)=1.0D-5*PDFMAX
        PLI=PDFI(1)*FACT
        PUI=PDFI(NIP)*FACT
        B(I)=1.0D0-1.0D0/(PLI*PUI*DX*DX)
        A(I)=(1.0D0/(PLI*DX))-1.0D0-B(I)
        C(I)=1.0D0+A(I)+B(I)
C
        ICASE=1
  300   CONTINUE
        ERR(I)=0.0D0
        DO K=1,NIP
          RR=SUMI(K)
          PAP=AREA(I)*(1.0D0+(A(I)+B(I)*RR)*RR)**2/
     1       ((1.0D0-B(I)*RR*RR)*C(I)*(X(I+1)-X(I)))
          IF(K.EQ.1.OR.K.EQ.NIP) THEN
            ERR(I)=ERR(I)+0.5D0*ABS(PAP-PDFI(K))
          ELSE
            ERR(I)=ERR(I)+ABS(PAP-PDFI(K))
          ENDIF
        ENDDO
        ERR(I)=ERR(I)*DXI
C
        IF(ERR(I).GT.0.10D0*AREA(I).AND.ICASE.EQ.1) THEN
          B(I)=0.0D0
          A(I)=0.0D0
          C(I)=1.0D0
          ICASE=2
          GO TO 300
        ENDIF
      ENDDO
C
      IF(NP.LT.N) GO TO 200
      NPM1=NP-1
C
C  ****  Renormalization.
C
      WS=0.0D0
      DO I=1,NPM1
        WS=WS+AREA(I)
      ENDDO
      WS=1.0D0/WS
      ERRM=0.0D0
      DO I=1,NPM1
        AREA(I)=AREA(I)*WS
        ERR(I)=ERR(I)*WS
        ERRM=MAX(ERRM,ERR(I))
      ENDDO
C
      PAC(1)=0.0D0
      DO I=1,NPM1
        PAC(I+1)=PAC(I)+AREA(I)
      ENDDO
      PAC(NP)=1.0D0
C
C  ****  Precalculated limits for the initial binary search in
C        subroutine QRNDI.
C
      BIN=1.0D0/DBLE(NPM1)
      ITTL(1)=1
      DO I=2,NPM1
        PTST=(I-1)*BIN
        DO J=ITTL(I-1),NP
          IF(PAC(J).GT.PTST) THEN
            ITTL(I)=J-1
            ITTU(I-1)=J
            GO TO 400
          ENDIF
        ENDDO
  400   CONTINUE
      ENDDO
      ITTU(NPM1)=NP
      ITTL(NP)=NP-1
      ITTU(NP)=NP
C
      RETURN
      END
C  *********************************************************************
C                       SUBROUTINE SIMPSU
C  *********************************************************************
      SUBROUTINE SIMPSU(H,F,SUMF,N)
C
C     Simpson's integration of a uniformly tabulated function.
C
C     H ...... step length,
C     F ...... array of function values (ordered abscissas),
C     SUMF ... array of integral values defined as
C              SUMF(I)=INTEGRAL(F) from X(1) to X(I)=X(1)+(I-1)*H,
C     N ...... number of points in the table.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N)
      PARAMETER (F1O12=1.0D0/12.0D0)
      DIMENSION F(N),SUMF(N)
      CONS=H*F1O12
      HCONS=0.5D0*CONS
      IF(N.LT.4) STOP 'Not enough data points.'
      SUMF(1)=0.0D0
      SUMF(2)=SUMF(1)+(5.0D0*F(1)+8.0D0*F(2)-F(3))*CONS
      DO I=3,N-1
        SUMF(I)=SUMF(I-1)+(13.0D0*(F(I-1)+F(I))-F(I+1)-F(I-2))*HCONS
      ENDDO
      SUMF(N)=SUMF(N-1)+(5.0D0*F(N)+8.0D0*F(N-1)-F(N-2))*CONS
      RETURN
      END
C  *********************************************************************
C                       SUBROUTINE QRNDMS
C  *********************************************************************
      SUBROUTINE QRNDMS(XD,XU,XM0,XM1,XM2)
C
C     Calculation of (restricted) momenta of a pdf, PDF(X), obtained
C  from its QRND rational approximation.
C
C     XD, XU ... limits of the integration interval.
C     XM0 ...... momentum of 0th order.
C     XM1 ...... momentum of 1st order.
C     XM2 ...... momentum of 2nd order.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N)
C
      PARAMETER (NM=512)
      COMMON/CQRND/XT(NM),PAC(NM),A(NM),B(NM),NP
      COMMON/CQRNDS/IGR
      PARAMETER (NIP=75)
      DIMENSION XI(NIP),FUN(NIP),SUM(NIP)
C
      XM0=0.0D0
      XM1=0.0D0
      XM2=0.0D0
      DO I=1,NP-1
        IF(XT(I+1).GE.XD.AND.XT(I).LE.XU) THEN
          X1=MAX(XT(I),XD)
          X2=MIN(XT(I+1),XU)
          IGR=I
C
          DX=(X2-X1)/DBLE(NIP-1)
          DO K=1,NIP
            XI(K)=X1+DBLE(K-1)*DX
            FUN(K)=QRNDM(XI(K))
          ENDDO
          CALL SIMPSU(DX,FUN,SUM,NIP)
          XM0=XM0+SUM(NIP)
C
          DO K=1,NIP
            FUN(K)=FUN(K)*XI(K)
          ENDDO
          CALL SIMPSU(DX,FUN,SUM,NIP)
          XM1=XM1+SUM(NIP)
C
          DO K=1,NIP
            FUN(K)=FUN(K)*XI(K)
          ENDDO
          CALL SIMPSU(DX,FUN,SUM,NIP)
          XM2=XM2+SUM(NIP)
        ENDIF
      ENDDO
      RETURN
      END
C  *********************************************************************
      FUNCTION QRNDM(X)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N)
      PARAMETER (NM=512)
      COMMON/CQRND/XT(NM),PAC(NM),A(NM),B(NM),NP
      COMMON/CQRNDS/IGR
C
      TAU=(X-XT(IGR))/(XT(IGR+1)-XT(IGR))
      CON1=2.0D0*B(IGR)*TAU
      CI=1.0D0+A(IGR)+B(IGR)
      CON2=CI-A(IGR)*TAU
      IF(ABS(CON1).GT.1.0D-10*ABS(CON2)) THEN
        ETA=CON2*(1.0D0-SQRT(1.0D0-2.0D0*TAU*CON1/CON2**2))/CON1
      ELSE
        ETA=TAU/CON2
      ENDIF
      QRNDM=(PAC(IGR+1)-PAC(IGR))*(1.0D0+(A(IGR)+B(IGR)*ETA)*ETA)**2
     1  /((1.0D0-B(IGR)*ETA*ETA)*CI*(XT(IGR+1)-XT(IGR)))
C
      RETURN
      END
