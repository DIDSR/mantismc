CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C        PPPPP   EEEEEE  N    N   GGGG   EEEEEE   OOOO   M    M        C
C        P    P  E       NN   N  G    G  E       O    O  MM  MM        C
C        P    P  E       N N  N  G       E       O    O  M MM M        C
C        PPPPP   EEEE    N  N N  G  GGG  EEEE    O    O  M    M        C
C        P       E       N   NN  G    G  E       O    O  M    M        C
C        P       EEEEEE  N    N   GGGG   EEEEEE   OOOO   M    M        C
C                                                                      C
C                                                   (version 2005).    C
C                                                                      C
C  Modular quadric geometry subroutines for Monte Carlo simulation     C
C  (using PENELOPE or detailed procedures) of radiation transport in   C
C  material systems.                                                   C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C  Changes introduced after version 2003 (with dates of introduction):
C
C  May 29, 2004. Quadric surfaces defined in implicit form can now be
C       rotated and translated.
C
C  May 29, 2004. The format of the input file has been modified to
C       allow 4-digit labels (needed for complex geometries). Geometry
C       files written for the 2003 version of PENGEOM must be adapted
C       to the new format (this requires small editing work).
C
C  June 6, 2004. Geometry definition files may include comment lines
C       (starting with the string 'C ' or 'c ') anywhere in the file.
C       Blank lines are not permitted.
C
C  June 6, 2004. A 'CLONE' operation for modules has been defined. This
C       simplifies the definition of repeated material structures.
C
C  June 7, 2004. Redundant surfaces (within round-off tolerance) are
C       removed. This speeds up the tracking of particles and reduces
C       the risk of having overlaps of adjacent bodies separated by
C       (nearly) coinciding surfaces.
C
C  June 7, 2004. A geometry report containing the genealogical modular
C       tree is issued when the parameter IREP (subroutine GEOMIN) is
C       given a positive value. This report is written on the output
C       file 'pengeom_tree.rep'. IREP is the unit number of this file;
C       you must make sure that it does not coincide with any of the
C       currently open units.
C
C  July 31, 2004. Subroutine STEP has been rewritten to implement a
C       more straightforward tracking strategy. Part of the iterative
C       operations has been separated in two new subroutines, STEPSI
C       (intersections with surfaces) and STEPLB (identification of
C       bodies from their surface side pointers).
C       Bugs found by Andreu Badal and Camille Guillerminet have been
C       corrected.
C
C  October 1, 2004. Each body can be made part of an impact detector.
C       The definition of impact detectors is through the common block
C         COMMON/QKDET/KDET(NB)
C       After the call to subroutine GEOMIN, a body KB is made part of
C       detector number IDET by setting KDET(KB)=IDET. Then, subroutine
C       STEP will stop particles when they enter body KB from vacuum or
C       from another body that is not part of detector IDET.
C
C       This new feature is used in the geometry viewer GVIEW2D to
C       display the various bodies intersected by the screen plane in
C       different colours. Notice that the body labels shown on the
C       screen are the ones used internally by PENELOPE. These internal
C       labels, which may be different from the user labels, are used
C       for specifying the simulation parameter DSMAX, for activating
C       variance reduction methods locally, and to define detectors.
C
C  April 22, 2005. Bodies can be limited not only by other bodies, but
C  also by predefined modules.
C
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
C  *********************************************************************
C                       SUBROUTINE GEOMIN
C  *********************************************************************
      SUBROUTINE GEOMIN(PARINP,NPINP,NMAT,NBOD,IRD,IWR)
C
C     This subroutine initialises the geometry package for Monte Carlo
C  simulation of particle transport.
C
C  Input arguments:
C     PARINP .... array containing optional parameters, which may
C                 replace the ones entered from the input file.
C     NPINP ..... number of parameters defined in PARINP (.ge.0).
C     IRD ....... input file unit (opened in the main program).
C     IWR ....... output file unit (opened in the main program).
C  Output arguments:
C     NMAT ...... number of different materials in full bodies (exclud-
C                 ing void regions).
C     NBOD ...... Number of defined bodies.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), CHARACTER*8 (L),
     1    INTEGER*4 (I-K,M-N)
      CHARACTER*12 GFILE
      CHARACTER*72 BLINE,DEFB,DEF
      PARAMETER (NS=100000,NB=40000,NX=1000)
      PARAMETER (PI=3.1415926535897932D0)
      PARAMETER (LNUL='00000000',LSUR='SURFACE ',LIND='INDICES=',
     1           LBOD='BODY    ',LMAT='MATERIAL',LMOD='MODULE  ',
     2          LOPEN=')       ',LEND='END     ',LONE='11111111',
     3           LXSC='X-SCALE=',LYSC='Y-SCALE=',LZSC='Z-SCALE=',
     4           LTHE='  THETA=',LPHI='    PHI=',LOME='  OMEGA=',
     5           LXSH='X-SHIFT=',LYSH='Y-SHIFT=',LZSH='Z-SHIFT=',
     6           LAXX='    AXX=',LAXY='    AXY=',LAXZ='    AXZ=',
     7           LAYY='    AYY=',LAYZ='    AYZ=',LAZZ='    AZZ=',
     8            LAX='     AX=', LAY='     AY=', LAZ='     AZ=',
     9            LA0='     A0=',LDEG=') DEG   ',LRAD=') RAD   ',
     A          LCLON='CLONE   ',LINC='INCLUDE ',LFIL='FILE=   ')
      DIMENSION PARINP(NPINP),LARRAY(8),KQ(5),KM(NS)
      DIMENSION IDESC(NB),IDONE(NB),ISCL(NS),ISOR(NS),IBCL(NB),IBOR(NB)
      COMMON/QSURF/AXX(NS),AXY(NS),AXZ(NS),AYY(NS),AYZ(NS),AZZ(NS),
     1    AX(NS),AY(NS),AZ(NS),A0(NS),NSURF,NS2,KPLANE(NS)
      COMMON/QBODY/KBODY(NB,NX),KALIAB(NB),KBOMO(NB)
      COMMON/QKDET/KDET(NB)
      COMMON/QTREE/NBODY,MATER(NB),KMOTH(NB),KDGHT(NB,NX),
     1    KSURF(NB,NX),KFLAG(NB,NX),KALIAS(NS),KSLAST
      PARAMETER (NLM=20)
      COMMON/QINFL/DEF(NLM),DEFB(NB,NLM),NLB(NB)
C
C  **** When the parameter IREP is set equal to a positive value, a
C  detailed report with the structure of the genealogical tree of
C  modules and with information on redundant surfaces is written on the
C  file 'pengeom_tree.rep' (UNIT=IREP).
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! FORCING NOT TO SAVE A FILE PENGEOM_TREE.REP FOR LINUX RUNS
      PARAMETER (IREP=00)  ! A detailed geometry report is written.
      ! PARAMETER (IREP=99)  ! A detailed geometry report is written.
C
C  ************  Initialise parameters.
C
C  ****  Surface coefficients, alias and KPLANE.
C
      DO KS=1,NS
        KALIAS(KS)=-10000  ! Surface aliases (user labels).
        AXX(KS)=0.0D0
        AXY(KS)=0.0D0
        AXZ(KS)=0.0D0
        AYY(KS)=0.0D0
        AYZ(KS)=0.0D0
        AZZ(KS)=0.0D0
        AX(KS)=0.0D0
        AY(KS)=0.0D0
        AZ(KS)=0.0D0
        A0(KS)=0.0D0
        KPLANE(KS)=0
C  KPLANE=1 if the surface is a plane.
      ENDDO
C  **** Alias, material, mother, daughters, surfaces and side pointers
C  of bodies. The last second component of the double arrays is the
C  number of used components. For example, KSURF(KB,NX) is the number
C  of surfaces that limit the body or module KB. The number of values
C  stored in array KFLAG(KB,--), however, is given by KSURF(KB,NX).
      DO KB=1,NB
C  **** The array KDET(KB) is used to label bodies that are parts of
C  impact detectors. Detectors must be defined in the main program,
C  after the call to subroutine GEOMIN.
        KDET(KB)=0  ! KDET(KB).ne.0 if body KB is part of a detector.
        KALIAB(KB)=-10000  ! Body aliases (user labels).
        MATER(KB)=0  ! Material in body KB.
        KBOMO(KB)=0  ! 0 for bodies, 1 for modules.
        KMOTH(KB)=0  ! Mother of body KB (unique).
        DO K2=1,NX
          KBODY(KB,K2)=0  ! Limiting bodies of body KB (keep it small).
C  KBODY(KB,K2)=1, if K2 is a body that limits KB.
C  KBODY(KB,K2)=0, if body K2 does not limit body KB.
          KDGHT(KB,K2)=0  ! Daughters of module KB (should be small).
C  KDGHT(KB,K2), is the K2-th daughter (body or module) of module KB.
          KSURF(KB,K2)=0  ! Limiting surfaces of body KB.
C  KSURF(KB,K2), is the K2-th limiting surface of body KB.
          KFLAG(KB,K2)=5  ! Surface side pointers of material bodies.
C  KFLAG(KB,KS)=1, if KS is a limiting surface of KB and KB is inside KS
C                  (i.e. side pointer =-1).
C  KFLAG(KB,KS)=2, if KS is a limiting surface of KB and KB is outside
C                  KS (i.e. side pointer =+1).
C  KFLAG(KB,KS)=3, if KB is a body and KS does not directly limit KB,
C                  but KS appears in the definition of a body that
C                  limits KB.
C  KFLAG(KB,KS)=4, if KB is a module and KS limits one of its daughters
C                  (bodies and submodules), but does not appear in the
C                  definition of KB.
C  KFLAG(KB,KS)=5, otherwise.
        ENDDO
      ENDDO
      NMAT=0
      NSURF=0
      NBODY=0
      ICLONE=0
C
C  ************  Reading the geometry input file.
C
      IR=IRD
      IW=IWR
      IF(IW.EQ.IR) THEN
        WRITE(IW,'(''SUBROUTINE GEOMIN. Input arguments.'')')
        WRITE(IW,'(''IRD ='',I4,'',  IWR ='',I4)') IRD,IWR
        WRITE(IW,'(''*** The input and output units must be '',
     1     ''different.'')')
        STOP
      ENDIF
      IRI=MAX(IR,IW)+1
C
 1    CONTINUE
      READ(IR,'(A72)') BLINE
      READ(BLINE,'(A8)') LKEYW
      IF(LKEYW.NE.LNUL) THEN
        IF(IR.EQ.IRD) THEN
          WRITE(IW,'(A72)') BLINE
        ELSE
          WRITE(IW,'(A2,A72)') 'C ',BLINE
        ENDIF
        GOTO 1
      ELSE
        IF(IR.EQ.IRD) THEN
          WRITE(IW,'(64(''0''))')
        ENDIF
      ENDIF
C
 2    CONTINUE
      READ(IR,'(A72)') BLINE
      IF(BLINE(1:2).EQ.'C '.OR.BLINE(1:2).EQ.'c ') THEN
        WRITE(IW,'(A72)') BLINE
        GOTO 2
      ENDIF
      READ(BLINE,'(A8,1X,I4,7A8))') LKEYW,IT,(LARRAY(I),I=1,7)
      IF(LKEYW.EQ.LSUR) THEN
        GOTO 100
      ELSE IF(LKEYW.EQ.LBOD) THEN
        GOTO 200
      ELSE IF(LKEYW.EQ.LMOD) THEN
        GOTO 300
      ELSE IF(LKEYW.EQ.LCLON) THEN
        GOTO 400
      ELSE IF(LKEYW.EQ.LINC) THEN
        GOTO 500
      ELSE IF(LKEYW.EQ.LEND) THEN
        GOTO 600
      ELSE
        WRITE(IW,'(A72)') BLINE
        WRITE(IW,'(''*** What do you mean?'')')
        STOP
      ENDIF
C
C  ************  Surfaces.
C
 100  CONTINUE
      READ(BLINE,'(9X,I4,7A8))') IT,(LARRAY(I),I=1,7)
      IF(NSURF.GT.0) THEN
        DO KS0=1,NSURF
          IF(IT.EQ.KALIAS(KS0)) THEN
            WRITE(IW,'(A72)') BLINE
            WRITE(IW,'(''*** Same label for two surfaces.'')')
            STOP
          ENDIF
        ENDDO
      ENDIF
      NSURF=NSURF+1
      KS=NSURF
      WRITE(IW,'(A8,''('',I4,7A8)') LKEYW,KS,(LARRAY(I),I=1,7)
      IF(KS.GT.NS) THEN
        WRITE(IW,'(''*** The parameter NS must be increased.'')')
        STOP
      ENDIF
      KALIAS(KS)=IT
C  ****  Indices.
 190  CONTINUE
      READ(IR,'(A72)') BLINE
      IF(BLINE(1:2).EQ.'C '.OR.BLINE(1:2).EQ.'c ') THEN
        WRITE(IW,'(A72)') BLINE
        GOTO 190
      ENDIF
      READ(BLINE,'(A8,1X,I2,1X,I2,1X,I2,1X,I2,1X,I2)')
     1  LKEYW,(KQ(I),I=1,5)
      WRITE(IW,'(A8,''('',4(I2,'',''),I2,'')'')')
     1  LKEYW,(KQ(I),I=1,5)
      IMODE=MAX(ABS(KQ(1)),ABS(KQ(2)),ABS(KQ(3)),ABS(KQ(4)),ABS(KQ(5)))
      IF(IMODE.GT.1.OR.LKEYW.NE.LIND) THEN
        WRITE(IW,'(''*** Incorrect surface indices.'')')
        STOP
      ENDIF
      IF(IMODE.EQ.0) GOTO 103
C
C  ****  Reduced form.
C
      XSCALE=1.0D0
      YSCALE=1.0D0
      ZSCALE=1.0D0
      OMEGA=0.0D0
      THETA=0.0D0
      PHI=0.0D0
      XSHIFT=0.0D0
      YSHIFT=0.0D0
      ZSHIFT=0.0D0
C  ****  Parameters of the quadric.
 101  CONTINUE
      READ(IR,'(A72)') BLINE
      IF(BLINE(1:2).EQ.'C '.OR.BLINE(1:2).EQ.'c ') THEN
        WRITE(IW,'(A72)') BLINE
        GOTO 101
      ENDIF
      READ(BLINE,'(A8,1X,E22.15,1X,I3,A8)') LKEYW,VALUE,ICHPAR,LANGLE
      IF(LKEYW.EQ.LNUL) GOTO 102
      IF(ICHPAR.LE.NPINP) THEN
        IF(ICHPAR.GT.0) THEN
          VALUE=PARINP(ICHPAR)
          ICHPAR=-ICHPAR  ! Switches off the parameter changing option.
        ENDIF
      ELSE
        WRITE(IW,'(A8,''('',1P,E22.15,'','',I3,A8)')
     1    LKEYW,VALUE,ICHPAR,LANGLE
        WRITE(IW,'(''*** NPINP is too small (check PARINP).'')')
        STOP
      ENDIF
C  ****  Transformation parameters.
      IF(LKEYW.EQ.LXSC) THEN
        WRITE(IW,'(A8,''('',1P,E22.15,'','',I3,A8)')
     1    LKEYW,VALUE,ICHPAR,LOPEN
        IF(VALUE.LT.1.0D-15) THEN
          WRITE(IW,'(''*** Scale factor less than 1.0D-15.'')')
          STOP
        ENDIF
        XSCALE=VALUE
      ELSE IF(LKEYW.EQ.LYSC) THEN
        WRITE(IW,'(A8,''('',1P,E22.15,'','',I3,A8)')
     1    LKEYW,VALUE,ICHPAR,LOPEN
        IF(VALUE.LT.1.0D-15) THEN
          WRITE(IW,'(''*** Scale factor less than 1.0D-15.'')')
          STOP
        ENDIF
        YSCALE=VALUE
      ELSE IF(LKEYW.EQ.LZSC) THEN
        WRITE(IW,'(A8,''('',1P,E22.15,'','',I3,A8)')
     1    LKEYW,VALUE,ICHPAR,LOPEN
        IF(VALUE.LT.1.0D-15) THEN
          WRITE(IW,'(''*** Scale factor less than 1.0D-15.'')')
          STOP
        ENDIF
        ZSCALE=VALUE
      ELSE IF(LKEYW.EQ.LOME) THEN
        IF(LANGLE.EQ.LRAD) THEN
          WRITE(IW,'(A8,''('',1P,E22.15,'','',I3,A8)')
     1      LKEYW,VALUE,ICHPAR,LRAD
          OMEGA=VALUE
        ELSE
          WRITE(IW,'(A8,''('',1P,E22.15,'','',I3,A8)')
     1      LKEYW,VALUE,ICHPAR,LDEG
          OMEGA=VALUE*PI/180.0D0
        ENDIF
      ELSE IF(LKEYW.EQ.LTHE) THEN
        IF(LANGLE.EQ.LRAD) THEN
          WRITE(IW,'(A8,''('',1P,E22.15,'','',I3,A8)')
     1      LKEYW,VALUE,ICHPAR,LRAD
          THETA=VALUE
        ELSE
          WRITE(IW,'(A8,''('',1P,E22.15,'','',I3,A8)')
     1      LKEYW,VALUE,ICHPAR,LDEG
          THETA=VALUE*PI/180.0D0
        ENDIF
      ELSE IF(LKEYW.EQ.LPHI) THEN
        IF(LANGLE.EQ.LRAD) THEN
          WRITE(IW,'(A8,''('',1P,E22.15,'','',I3,A8)')
     1      LKEYW,VALUE,ICHPAR,LRAD
          PHI=VALUE
        ELSE
          WRITE(IW,'(A8,''('',1P,E22.15,'','',I3,A8)')
     1      LKEYW,VALUE,ICHPAR,LDEG
          PHI=VALUE*PI/180.0D0
        ENDIF
      ELSE IF(LKEYW.EQ.LXSH) THEN
        WRITE(IW,'(A8,''('',1P,E22.15,'','',I3,A8)')
     1    LKEYW,VALUE,ICHPAR,LOPEN
        XSHIFT=VALUE
      ELSE IF(LKEYW.EQ.LYSH) THEN
        WRITE(IW,'(A8,''('',1P,E22.15,'','',I3,A8)')
     1    LKEYW,VALUE,ICHPAR,LOPEN
        YSHIFT=VALUE
      ELSE IF(LKEYW.EQ.LZSH) THEN
        WRITE(IW,'(A8,''('',1P,E22.15,'','',I3,A8)')
     1    LKEYW,VALUE,ICHPAR,LOPEN
        ZSHIFT=VALUE
      ELSE
        WRITE(IW,'(A72)') BLINE
        WRITE(IW,'(''*** What do you mean?'')')
        STOP
      ENDIF
      GOTO 101
C  ****  Expanded quadric.
 102  CONTINUE
      QXX=KQ(1)/XSCALE**2
      QXY=0.0D0
      QXZ=0.0D0
      QYY=KQ(2)/YSCALE**2
      QYZ=0.0D0
      QZZ=KQ(3)/ZSCALE**2
      QX=0.0D0
      QY=0.0D0
      QZ=KQ(4)/ZSCALE
      Q0=KQ(5)
C  ****  Rotated and shifted quadric.
      CALL ROTSHF(OMEGA,THETA,PHI,XSHIFT,YSHIFT,ZSHIFT,
     1            QXX,QXY,QXZ,QYY,QYZ,QZZ,QX,QY,QZ,Q0)
      IF(MAX(ABS(MIN(QXX,QXY,QXZ,QYY,QYZ,QZZ)),
     1   ABS(MAX(QXX,QXY,QXZ,QYY,QYZ,QZZ))).LT.1.0D-30) KPLANE(KS)=1
      AXX(KS)=QXX
      AXY(KS)=QXY
      AXZ(KS)=QXZ
      AYY(KS)=QYY
      AYZ(KS)=QYZ
      AZZ(KS)=QZZ
      AX(KS)=QX
      AY(KS)=QY
      AZ(KS)=QZ
      A0(KS)=Q0
      WRITE(IW,'(64(''0''))')
      GOTO 2
C
C  ****  Implicit form.
C
 103  CONTINUE
      QXX=0.0D0
      QXY=0.0D0
      QXZ=0.0D0
      QYY=0.0D0
      QYZ=0.0D0
      QZZ=0.0D0
      QX=0.0D0
      QY=0.0D0
      QZ=0.0D0
      Q0=0.0D0
C
 193  CONTINUE
      READ(IR,'(A72)') BLINE
      IF(BLINE(1:2).EQ.'C '.OR.BLINE(1:2).EQ.'c ') THEN
        WRITE(IW,'(A72)') BLINE
        GOTO 193
      ENDIF
      READ(BLINE,'(A8,1X,E22.15,1X,I3,A8)') LKEYW,VALUE,ICHPAR,LANGLE
      IF(LKEYW.EQ.LNUL) GOTO 107
      IF(LKEYW.EQ.LONE) GOTO 104
      IF(ICHPAR.LE.NPINP) THEN
        IF(ICHPAR.GT.0) THEN
          VALUE=PARINP(ICHPAR)
          ICHPAR=-ICHPAR
        ENDIF
      ELSE
        WRITE(IW,'(A8,''('',1P,E22.15,'','',I3,A8)')
     1    LKEYW,VALUE,ICHPAR,LANGLE
        WRITE(IW,'(''*** NPINP is too small (check PARINP).'')')
        STOP
      ENDIF
      IF(LKEYW.EQ.LAXX) THEN
        QXX=VALUE
        WRITE(IW,'(A8,''('',1P,E22.15,'','',I3,A8)')
     1    LAXX,QXX,ICHPAR,LOPEN
      ELSE IF(LKEYW.EQ.LAXY) THEN
        QXY=VALUE
        WRITE(IW,'(A8,''('',1P,E22.15,'','',I3,A8)')
     1    LAXY,QXY,ICHPAR,LOPEN
      ELSE IF(LKEYW.EQ.LAXZ) THEN
        QXZ=VALUE
        WRITE(IW,'(A8,''('',1P,E22.15,'','',I3,A8)')
     1    LAXZ,QXZ,ICHPAR,LOPEN
      ELSE IF(LKEYW.EQ.LAYY) THEN
        QYY=VALUE
        WRITE(IW,'(A8,''('',1P,E22.15,'','',I3,A8)')
     1    LAYY,QYY,ICHPAR,LOPEN
      ELSE IF(LKEYW.EQ.LAYZ) THEN
        QYZ=VALUE
        WRITE(IW,'(A8,''('',1P,E22.15,'','',I3,A8)')
     1    LAYZ,QYZ,ICHPAR,LOPEN
      ELSE IF(LKEYW.EQ.LAZZ) THEN
        QZZ=VALUE
        WRITE(IW,'(A8,''('',1P,E22.15,'','',I3,A8)')
     1    LAZZ,QZZ,ICHPAR,LOPEN
      ELSE IF(LKEYW.EQ.LAX) THEN
        QX=VALUE
        WRITE(IW,'(A8,''('',1P,E22.15,'','',I3,A8)')
     1    LAX,QX,ICHPAR,LOPEN
      ELSE IF(LKEYW.EQ.LAY) THEN
        QY=VALUE
        WRITE(IW,'(A8,''('',1P,E22.15,'','',I3,A8)')
     1    LAY,QY,ICHPAR,LOPEN
      ELSE IF(LKEYW.EQ.LAZ) THEN
        QZ=VALUE
        WRITE(IW,'(A8,''('',1P,E22.15,'','',I3,A8)')
     1    LAZ,QZ,ICHPAR,LOPEN
      ELSE IF(LKEYW.EQ.LA0) THEN
        Q0=VALUE
        WRITE(IW,'(A8,''('',1P,E22.15,'','',I3,A8)')
     1    LA0,Q0,ICHPAR,LOPEN
      ELSE
        WRITE(IW,'(A72)') BLINE
        WRITE(IW,'(''*** What do you mean?'')')
        STOP
      ENDIF
      GOTO 193
C  ****  Transformation parameters.
 104  CONTINUE
      WRITE(IW,'(64(''1''))')
      OMEGA=0.0D0
      THETA=0.0D0
      PHI=0.0D0
      XSHIFT=0.0D0
      YSHIFT=0.0D0
      ZSHIFT=0.0D0
C
 105  CONTINUE
      READ(IR,'(A72)') BLINE
      IF(BLINE(1:2).EQ.'C '.OR.BLINE(1:2).EQ.'c ') THEN
        WRITE(IW,'(A72)') BLINE
        GOTO 105
      ENDIF
      READ(BLINE,'(A8,1X,E22.15,1X,I3,A8)') LKEYW,VALUE,ICHPAR,LANGLE
      IF(LKEYW.EQ.LNUL) GOTO 106
      IF(ICHPAR.LE.NPINP) THEN
        IF(ICHPAR.GT.0) THEN
          VALUE=PARINP(ICHPAR)
          ICHPAR=-ICHPAR  ! Switches off the parameter changing option.
        ENDIF
      ELSE
        WRITE(IW,'(A8,''('',1P,E22.15,'','',I3,A8)')
     1    LKEYW,VALUE,ICHPAR,LANGLE
        WRITE(IW,'(''*** NPINP is too small (check PARINP).'')')
        STOP
      ENDIF
C
      IF(LKEYW.EQ.LOME) THEN
        IF(LANGLE.EQ.LRAD) THEN
          WRITE(IW,'(A8,''('',1P,E22.15,'','',I3,A8)')
     1      LKEYW,VALUE,ICHPAR,LRAD
          OMEGA=VALUE
        ELSE
          WRITE(IW,'(A8,''('',1P,E22.15,'','',I3,A8)')
     1      LKEYW,VALUE,ICHPAR,LDEG
          OMEGA=VALUE*PI/180.0D0
        ENDIF
      ELSE IF(LKEYW.EQ.LTHE) THEN
        IF(LANGLE.EQ.LRAD) THEN
          WRITE(IW,'(A8,''('',1P,E22.15,'','',I3,A8)')
     1      LKEYW,VALUE,ICHPAR,LRAD
          THETA=VALUE
        ELSE
          WRITE(IW,'(A8,''('',1P,E22.15,'','',I3,A8)')
     1      LKEYW,VALUE,ICHPAR,LDEG
          THETA=VALUE*PI/180.0D0
        ENDIF
      ELSE IF(LKEYW.EQ.LPHI) THEN
        IF(LANGLE.EQ.LRAD) THEN
          WRITE(IW,'(A8,''('',1P,E22.15,'','',I3,A8)')
     1      LKEYW,VALUE,ICHPAR,LRAD
          PHI=VALUE
        ELSE
          WRITE(IW,'(A8,''('',1P,E22.15,'','',I3,A8)')
     1      LKEYW,VALUE,ICHPAR,LDEG
          PHI=VALUE*PI/180.0D0
        ENDIF
      ELSE IF(LKEYW.EQ.LXSH) THEN
        WRITE(IW,'(A8,''('',1P,E22.15,'','',I3,A8)')
     1    LKEYW,VALUE,ICHPAR,LOPEN
        XSHIFT=VALUE
      ELSE IF(LKEYW.EQ.LYSH) THEN
        WRITE(IW,'(A8,''('',1P,E22.15,'','',I3,A8)')
     1    LKEYW,VALUE,ICHPAR,LOPEN
        YSHIFT=VALUE
      ELSE IF(LKEYW.EQ.LZSH) THEN
        WRITE(IW,'(A8,''('',1P,E22.15,'','',I3,A8)')
     1    LKEYW,VALUE,ICHPAR,LOPEN
        ZSHIFT=VALUE
      ELSE
        WRITE(IW,'(A72)') BLINE
        WRITE(IW,'(''*** What do you mean?'')')
        STOP
      ENDIF
      GOTO 105
C  ****  Rotation and translation of the surface.
 106  CONTINUE
      CALL ROTSHF(OMEGA,THETA,PHI,XSHIFT,YSHIFT,ZSHIFT,
     1            QXX,QXY,QXZ,QYY,QYZ,QZZ,QX,QY,QZ,Q0)
C
 107  CONTINUE
      TSTL=MIN(QXX,QXY,QXZ,QYY,QYZ,QZZ)
      TSTU=MAX(QXX,QXY,QXZ,QYY,QYZ,QZZ)
      IF(MAX(ABS(TSTL),ABS(TSTU)).LT.1.0D-30) KPLANE(KS)=1
      AXX(KS)=QXX
      AXY(KS)=QXY
      AXZ(KS)=QXZ
      AYY(KS)=QYY
      AYZ(KS)=QYZ
      AZZ(KS)=QZZ
      AX(KS)=QX
      AY(KS)=QY
      AZ(KS)=QZ
      A0(KS)=Q0
      WRITE(IW,'(64(''0''))')
      GOTO 2
C
C  ************  Bodies.
C
 200  CONTINUE
      READ(BLINE,'(9X,I4,7A8))') IT,(LARRAY(I),I=1,7)
      IF(NBODY.GT.0) THEN
        DO KB0=1,NBODY
          IF(IT.EQ.KALIAB(KB0)) THEN
            WRITE(IW,'(A72)') BLINE
            WRITE(IW,
     1        '(''*** Same label for two bodies (or modules).'')')
            STOP
          ENDIF
        ENDDO
      ENDIF
      NBODY=NBODY+1
      WRITE(IW,'(A8,''('',I4,7A8)') LKEYW,NBODY,(LARRAY(I),I=1,7)
      IF(NBODY.GT.NB) THEN
        WRITE(IW,'(''*** The parameter NB must be increased.'')')
        STOP
      ENDIF
      NL=1
      WRITE(DEF(NL),'(A8,''('',I4,7A8)') LKEYW,NBODY,(LARRAY(I),I=1,7)
      KALIAB(NBODY)=IT
C
 295  CONTINUE
      READ(IR,'(A72)') BLINE
      IF(BLINE(1:2).EQ.'C '.OR.BLINE(1:2).EQ.'c ') THEN
        WRITE(IW,'(A72)') BLINE
        GOTO 295
      ENDIF
      READ(BLINE,'(A8,1X,I4)') LKEYW,IMAT
      IF(LKEYW.NE.LMAT) THEN
        WRITE(IW,'(A72)') BLINE
        WRITE(IW,'(''*** Incorrect material definition line.'')')
        STOP
      ENDIF
      WRITE(IW,'(A8,''('',I4,'')'')') LKEYW,IMAT
      NL=MIN(NL+1,NLM)
      WRITE(DEF(NL),'(A8,''('',I4,'')'')') LKEYW,IMAT
      IF(IMAT.LT.0) IMAT=0
      MATER(NBODY)=IMAT
      NMAT=MAX(NMAT,IMAT)
C
 201  CONTINUE
      READ(IR,'(A72)') BLINE
      IF(BLINE(1:2).EQ.'C '.OR.BLINE(1:2).EQ.'c ') THEN
        WRITE(IW,'(A72)') BLINE
        GOTO 201
      ENDIF
      READ(BLINE,'(A8,1X,I4,17X,I2)') LKEYW,KIN,INDS
      IF(LKEYW.EQ.LNUL) GOTO 208
      IF(LKEYW.EQ.LSUR) THEN
C  ****  Surface.
        DO KS0=1,NSURF
          IF(KIN.EQ.KALIAS(KS0)) THEN
            KS=KS0
            GOTO 202
          ENDIF
        ENDDO
        WRITE(IW,'(A72)') BLINE
        WRITE(IW,'(''*** Undefined surface label.'')')
        STOP
 202    CONTINUE
        WRITE(IW,'(A8,''('',I4,''), SIDE POINTER=('',I2,'')'')')
     1    LKEYW,KS,INDS
        NL=MIN(NL+1,NLM)
        WRITE(DEF(NL),'(A8,''('',I4,''), SIDE POINTER=('',I2,'')'')')
     1    LKEYW,KS,INDS
C
        KST=KSURF(NBODY,NX)
        IF(KST.GT.0) THEN
          KST0=KST
          DO K=1,KST0
            IF(KS.EQ.KSURF(NBODY,K)) THEN
              IF(KFLAG(NBODY,K).LT.3) THEN
                WRITE(IW,'(''*** The last limiting surface has '',
     1            ''been defined twice.'')')
                STOP
              ELSE
                KST=K  ! KS limits a limiting body.
                GOTO 203
              ENDIF
            ENDIF
          ENDDO
        ENDIF
        KST=KST+1
        IF(KST.GE.NX) THEN
          WRITE(IW,'(''*** The parameter NX is too small.'')')
          STOP
        ENDIF
        KSURF(NBODY,NX)=KST
        KSURF(NBODY,KST)=KS
 203    CONTINUE
        IF(INDS.EQ.-1) THEN
          KFLAG(NBODY,KST)=1
        ELSE IF(INDS.EQ.1) THEN
          KFLAG(NBODY,KST)=2
        ELSE
          WRITE(IW,'(A72)') BLINE
          WRITE(IW,'(''*** Check side pointer value.'')')
          STOP
        ENDIF
      ELSE IF(LKEYW.EQ.LBOD) THEN
C  ****  Body.
        IF(NBODY.GT.1) THEN
          DO KB0=1,NBODY-1
            IF(KIN.EQ.KALIAB(KB0)) THEN
              KB=KB0
              GOTO 204
            ENDIF
          ENDDO
          WRITE(IW,'(A72)') BLINE
          WRITE(IW,'(''*** Undefined body label.'')')
          STOP
        ENDIF
 204    CONTINUE
        WRITE(IW,'(A8,''('',I4,'')'')') LKEYW,KB
        IF(KBOMO(KB).NE.0) THEN
          WRITE(IW,'(''*** This body is a module.'')')
          STOP
        ENDIF
        NL=MIN(NL+1,NLM)
        WRITE(DEF(NL),'(A8,''('',I4,'')'')') LKEYW,KB
C
        KN1=KSURF(KB,NX)
        DO 205 KS1=1,KN1
          KSURF1=KSURF(KB,KS1)
          KN2=KSURF(NBODY,NX)
          DO KS2=1,KN2
            KSURF2=KSURF(NBODY,KS2)
            IF(KSURF2.EQ.KSURF1) GOTO 205
          ENDDO
          KST=KN2+1
          IF(KST.GE.NX) THEN
            WRITE(IW,'(A72)') BLINE
            WRITE(IW,'(''*** The parameter NX is too small.'')')
            STOP
          ENDIF
          KSURF(NBODY,NX)=KST
          KSURF(NBODY,KST)=KSURF1
          KFLAG(NBODY,KST)=3
 205    CONTINUE
        KBODY(NBODY,NX)=KBODY(NBODY,NX)+1
        KBODY(NBODY,KBODY(NBODY,NX))=KB
      ELSE IF(LKEYW.EQ.LMOD) THEN
C  ****  Module.
        IF(NBODY.GT.1) THEN
          DO KB0=1,NBODY-1
            IF(KIN.EQ.KALIAB(KB0)) THEN
              KB=KB0
              GOTO 206
            ENDIF
          ENDDO
          WRITE(IW,'(A72)') BLINE
          WRITE(IW,'(''*** Undefined body label.'')')
          STOP
        ENDIF
 206    CONTINUE
        WRITE(IW,'(A8,''('',I4,'')'')') LKEYW,KB
        IF(KBOMO(KB).NE.1) THEN
          WRITE(IW,'(''*** This module is a body.'')')
          STOP
        ENDIF
        NL=MIN(NL+1,NLM)
        WRITE(DEF(NL),'(A8,''('',I4,'')'')') LKEYW,KB
C
        KN1=KSURF(KB,NX)
        DO 207 KS1=1,KN1
          KSURF1=KSURF(KB,KS1)
          IF(KFLAG(KB,KS1).GT.2) GOTO 207
          KN2=KSURF(NBODY,NX)
          DO KS2=1,KN2
            KSURF2=KSURF(NBODY,KS2)
            IF(KSURF2.EQ.KSURF1) GOTO 207
          ENDDO
          KST=KN2+1
          IF(KST.GE.NX) THEN
            WRITE(IW,'(A72)') BLINE
            WRITE(IW,'(''*** The parameter NX is too small.'')')
            STOP
          ENDIF
          KSURF(NBODY,NX)=KST
          KSURF(NBODY,KST)=KSURF1
          KFLAG(NBODY,KST)=3
 207    CONTINUE
        KBODY(NBODY,NX)=KBODY(NBODY,NX)+1
        KBODY(NBODY,KBODY(NBODY,NX))=KB
      ELSE
        WRITE(IW,'(A72)') BLINE
        WRITE(IW,'(''*** What do you mean?'')')
        STOP
      ENDIF
      GOTO 201
 208  CONTINUE
      WRITE(IW,'(64(''0''))')
      NL=MIN(NL+1,NLM)
      WRITE(DEF(NL),'(64(''0''))')
C
      NLB(NBODY)=NL
      DO I=1,NL
        DEFB(NBODY,I)=DEF(I)
C       WRITE(IW,'(''C '',A72)') DEFB(NBODY,I)
      ENDDO
      GOTO 2
C
C  ************  Modules.
C
 300  CONTINUE
      READ(BLINE,'(9X,I4,7A8))') IT,(LARRAY(I),I=1,7)
      IF(NBODY.GT.0) THEN
        DO KB0=1,NBODY
        IF(IT.EQ.KALIAB(KB0)) THEN
          WRITE(IW,'(A72)') BLINE
          WRITE(IW,
     1      '(''*** Same label for two bodies (or modules).'')')
          STOP
        ENDIF
        ENDDO
      ENDIF
      NBODY=NBODY+1
      WRITE(IW,'(A8,''('',I4,7A8)') LKEYW,NBODY,(LARRAY(I),I=1,7)
      IF(NBODY.GT.NB) THEN
        WRITE(IW,'(''*** The parameter NB must be increased.'')')
        STOP
      ENDIF
      NL=1
      WRITE(DEF(NL),'(A8,''('',I4,7A8)') LKEYW,NBODY,(LARRAY(I),I=1,7)
      KALIAB(NBODY)=IT
C
 391  CONTINUE
      READ(IR,'(A72)') BLINE
      IF(BLINE(1:2).EQ.'C '.OR.BLINE(1:2).EQ.'c ') THEN
        WRITE(IW,'(A72)') BLINE
        GOTO 391
      ENDIF
      READ(BLINE,'(A8,1X,I4)') LKEYW,IMAT
      WRITE(IW,'(A8,''('',I4,'')'')') LKEYW,IMAT
      IF(LKEYW.NE.LMAT) THEN
        WRITE(IW,'(''*** Incorrect material definition line.'')')
        STOP
      ENDIF
      NL=MIN(NL+1,NLM)
      WRITE(DEF(NL),'(A8,''('',I4,'')'')') LKEYW,IMAT
      IF(IMAT.LT.0) IMAT=0
      MATER(NBODY)=IMAT
      NMAT=MAX(NMAT,IMAT)
C
      KDGHT(NBODY,NX)=1
      KDGHT(NBODY,1)=NBODY
C
C  ****  Limiting surfaces and components.
C
 301  CONTINUE
      READ(IR,'(A72)') BLINE
      IF(BLINE(1:2).EQ.'C '.OR.BLINE(1:2).EQ.'c ') THEN
        WRITE(IW,'(A72)') BLINE
        GOTO 301
      ENDIF
      READ(BLINE,'(A8,1X,I4,17X,I2)') LKEYW,KIN,INDS
      IF(LKEYW.EQ.LNUL.OR.LKEYW.EQ.LONE) THEN
        KDT=KDGHT(NBODY,NX)
C  ****  Sort daughters in increasing order.
        IF(KDT.GT.1) THEN
          DO KI=1,KDT-1
            KBMIN=KDGHT(NBODY,KI)
            KMIN=KI
            DO KJ=KI+1,KDT
              IF(KDGHT(NBODY,KJ).LT.KBMIN) THEN
                KBMIN=KDGHT(NBODY,KJ)
                KMIN=KJ
              ENDIF
            ENDDO
            IF(KMIN.NE.KI) THEN
              KSAVE=KDGHT(NBODY,KI)
              KDGHT(NBODY,KI)=KDGHT(NBODY,KMIN)
              KDGHT(NBODY,KMIN)=KSAVE
            ENDIF
          ENDDO
        ENDIF
        IF(LKEYW.EQ.LONE) GOTO 308
        IF(LKEYW.EQ.LNUL) GOTO 312
      ENDIF
C
C  ****  Limiting surface.
C
      IF(LKEYW.EQ.LSUR) THEN
        DO KS0=1,NSURF
          IF(KIN.EQ.KALIAS(KS0)) THEN
            KS=KS0
            GOTO 302
          ENDIF
        ENDDO
        WRITE(IW,'(A8,''('',I4,''), SIDE POINTER=('',I2,'')'')')
     1    LKEYW,KIN,INDS
        WRITE(IW,'(''*** Undefined surface label.'')')
        STOP
 302    CONTINUE
        WRITE(IW,'(A8,''('',I4,''), SIDE POINTER=('',I2,'')'')')
     1    LKEYW,KS,INDS
        NL=MIN(NL+1,NLM)
        WRITE(DEF(NL),'(A8,''('',I4,''), SIDE POINTER=('',I2,'')'')')
     1    LKEYW,KS,INDS
C
        KST=KSURF(NBODY,NX)
        IF(KST.GT.0) THEN  ! Check whether KS is in the list.
          KST0=KST
          DO K=1,KST0
            IF(KS.EQ.KSURF(NBODY,K)) THEN
              IF(KFLAG(NBODY,K).LT.3) THEN
                WRITE(IW,'(''*** The last limiting surface has '',
     1            ''been defined twice.'')')
                STOP
              ELSE
                KST=K
                GOTO 303
              ENDIF
            ENDIF
          ENDDO
        ENDIF
        KST=KST+1
        IF(KST.GE.NX) THEN
          WRITE(IW,'(''*** The parameter NX is too small.'')')
          STOP
        ENDIF
        KSURF(NBODY,NX)=KST
        KSURF(NBODY,KST)=KS
 303    CONTINUE
        IF(INDS.EQ.-1) THEN
          KFLAG(NBODY,KST)=1
        ELSE IF(INDS.EQ.1) THEN
          KFLAG(NBODY,KST)=2
        ELSE
          WRITE(IW,'(''*** Check side pointer value.'')')
          STOP
        ENDIF
C
C  ****  Body.
C
      ELSE IF(LKEYW.EQ.LBOD) THEN
        IF(NBODY.GT.1) THEN  ! Check whether KB is in the list.
          DO KB0=1,NBODY-1
            IF(KIN.EQ.KALIAB(KB0)) THEN
              KB=KB0
              GOTO 304
            ENDIF
          ENDDO
          WRITE(IW,'(A8,''('',I4,'')'')') LKEYW,KIN
          WRITE(IW,'(''*** Undefined body label.'')')
          STOP
        ENDIF
 304    CONTINUE
        WRITE(IW,'(A8,''('',I4,'')'')') LKEYW,KB
        IF(KBOMO(KB).NE.0) THEN
          WRITE(IW,'(''*** This body is a module.'')')
          STOP
        ENDIF
        NL=MIN(NL+1,NLM)
        WRITE(DEF(NL),'(A8,''('',I4,'')'')') LKEYW,KB
        IF(KMOTH(KB).GT.0.AND.KMOTH(KB).NE.NBODY) THEN
          WRITE(IW,'(''*** You are trying to assign two mothers to ''
     1      ''the last body.'')')
          STOP
        ENDIF
        KMOTH(KB)=NBODY
        KDT=KDGHT(NBODY,NX)+1
        KDGHT(NBODY,NX)=KDT
        KDGHT(NBODY,KDT)=KB
        DO K2=1,KBODY(KB,NX)
          IF(KMOTH(KBODY(KB,K2)).EQ.0) KMOTH(KBODY(KB,K2))=NBODY
        ENDDO
C  ****  Assign genealogical attributes to the sisters of body KB.
        K2M=KBODY(KB,NX)
        DO K2=1,K2M
          KBD=KBODY(KB,K2)
          IF(KMOTH(KBD).EQ.0) KMOTH(KBD)=NBODY
          IDGHT=0
          DO K=1,KDT
            IF(KDGHT(NBODY,K).EQ.KBD) IDGHT=K
          ENDDO
          IF(IDGHT.EQ.0) THEN
            KDGHT(NBODY,NX)=KDGHT(NBODY,NX)+1
            KDGHT(NBODY,KDGHT(NBODY,NX))=KB
          ENDIF
        ENDDO
C  ****  Surfaces of the sister bodies.
        KN1=KSURF(KB,NX)
        DO 317 KS1=1,KN1
          IF(KFLAG(KB,KS1).GT.3) GOTO 317
          KSURF1=KSURF(KB,KS1)
          KN2=KSURF(NBODY,NX)
          DO KS2=1,KN2
            KSURF2=KSURF(NBODY,KS2)
            IF(KSURF2.EQ.KSURF1) GOTO 317
          ENDDO
          KN2=KN2+1
          IF(KN2.GE.NX) THEN
            WRITE(IW,'(''*** The parameter NX is too small.'')')
            STOP
          ENDIF
          KSURF(NBODY,NX)=KN2
          KSURF(NBODY,KN2)=KSURF1
          KFLAG(NBODY,KN2)=4
 317    CONTINUE
C
C  ****  Module.
C
      ELSE IF(LKEYW.EQ.LMOD) THEN  ! Check whether KB is in the list.
        IF(NBODY.GT.0) THEN
          DO KB0=1,NBODY-1
            IF(KIN.EQ.KALIAB(KB0)) THEN
              KB=KB0
              GOTO 305
            ENDIF
          ENDDO
          WRITE(IW,'(A8,''('',I4,'')'')') LKEYW,KIN
          WRITE(IW,'(''*** Undefined body label.'')')
          STOP
        ENDIF
 305    CONTINUE
        WRITE(IW,'(A8,''('',I4,'')'')') LKEYW,KB
        IF(KBOMO(KB).NE.1) THEN
          WRITE(IW,'(''*** This module is a body.'')')
          STOP
        ENDIF
        NL=MIN(NL+1,NLM)
        WRITE(DEF(NL),'(A8,''('',I4,'')'')') LKEYW,KB
        IF(KMOTH(KB).GT.0.AND.KMOTH(KB).NE.NBODY) THEN
          WRITE(IW,'(''*** You are trying to assign two mothers to '',
     1      ''the last module.'')')
          STOP
        ENDIF
        KMOTH(KB)=NBODY
        KDT=KDGHT(NBODY,NX)+1
        KDGHT(NBODY,NX)=KDT
        KDGHT(NBODY,KDT)=KB
C  ****  Limiting surfaces.
        KN1=KSURF(KB,NX)
        DO 307 KS1=1,KN1
          IF(KFLAG(KB,KS1).GT.2) GOTO 307
          KSURF1=KSURF(KB,KS1)
          KN2=KSURF(NBODY,NX)
          DO KS2=1,KN2
            KSURF2=KSURF(NBODY,KS2)
            IF(KSURF2.EQ.KSURF1) GOTO 307
          ENDDO
          KN2=KN2+1
          IF(KN2.GE.NX) THEN
            WRITE(IW,'(''*** The parameter NX is too small.'')')
            STOP
          ENDIF
          KSURF(NBODY,NX)=KN2
          KSURF(NBODY,KN2)=KSURF1
          KFLAG(NBODY,KN2)=4
 307    CONTINUE
C
      ELSE
        WRITE(IW,'(A72)') BLINE
        WRITE(IW,'(''*** What do you mean?'')')
        STOP
      ENDIF
      GOTO 301
C
C  ****  Transformation parameters.
C
 308  CONTINUE
      WRITE(IW,'(64(''1''))')
      NL=MIN(NL+1,NLM)
      WRITE(DEF(NL),'(64(''1''))')
      OMEGA=0.0D0
      THETA=0.0D0
      PHI=0.0D0
      XSHIFT=0.0D0
      YSHIFT=0.0D0
      ZSHIFT=0.0D0
C
 309  CONTINUE
      READ(IR,'(A72)') BLINE
      IF(BLINE(1:2).EQ.'C '.OR.BLINE(1:2).EQ.'c ') THEN
        WRITE(IW,'(A72)') BLINE
        GOTO 309
      ENDIF
      READ(BLINE,'(A8,1X,E22.15,1X,I3,A8)') LKEYW,VALUE,ICHPAR,LANGLE
      IF(LKEYW.EQ.LNUL) GOTO 310
      IF(ICHPAR.LE.NPINP) THEN
        IF(ICHPAR.GT.0) THEN
          VALUE=PARINP(ICHPAR)
          ICHPAR=-ICHPAR
        ENDIF
      ELSE
        WRITE(IW,'(A72)') BLINE
        WRITE(IW,'(''*** NPINP is too small (check PARINP).'')')
        STOP
      ENDIF
      IF(LKEYW.EQ.LOME) THEN
        IF(LANGLE.EQ.LRAD) THEN
          WRITE(IW,'(A8,''('',1P,E22.15,'','',I3,A8)')
     1      LKEYW,VALUE,ICHPAR,LRAD
          NL=MIN(NL+1,NLM)
          WRITE(DEF(NL),'(A8,''('',1P,E22.15,'','',I3,A8)')
     1      LKEYW,VALUE,ICHPAR,LRAD
          OMEGA=VALUE
        ELSE
          WRITE(IW,'(A8,''('',1P,E22.15,'','',I3,A8)')
     1      LKEYW,VALUE,ICHPAR,LDEG
          NL=MIN(NL+1,NLM)
          WRITE(DEF(NL),'(A8,''('',1P,E22.15,'','',I3,A8)')
     1      LKEYW,VALUE,ICHPAR,LDEG
          OMEGA=VALUE*PI/180.0D0
        ENDIF
      ELSE IF(LKEYW.EQ.LTHE) THEN
        IF(LANGLE.EQ.LRAD) THEN
          WRITE(IW,'(A8,''('',1P,E22.15,'','',I3,A8)')
     1      LKEYW,VALUE,ICHPAR,LRAD
          NL=MIN(NL+1,NLM)
          WRITE(DEF(NL),'(A8,''('',1P,E22.15,'','',I3,A8)')
     1      LKEYW,VALUE,ICHPAR,LRAD
          THETA=VALUE
        ELSE
          WRITE(IW,'(A8,''('',1P,E22.15,'','',I3,A8)')
     1      LKEYW,VALUE,ICHPAR,LDEG
          NL=MIN(NL+1,NLM)
          WRITE(DEF(NL),'(A8,''('',1P,E22.15,'','',I3,A8)')
     1      LKEYW,VALUE,ICHPAR,LDEG
          THETA=VALUE*PI/180.0D0
        ENDIF
      ELSE IF(LKEYW.EQ.LPHI) THEN
        IF(LANGLE.EQ.LRAD) THEN
          WRITE(IW,'(A8,''('',1P,E22.15,'','',I3,A8)')
     1      LKEYW,VALUE,ICHPAR,LRAD
          NL=MIN(NL+1,NLM)
          WRITE(DEF(NL),'(A8,''('',1P,E22.15,'','',I3,A8)')
     1      LKEYW,VALUE,ICHPAR,LRAD
          PHI=VALUE
        ELSE
          WRITE(IW,'(A8,''('',1P,E22.15,'','',I3,A8)')
     1      LKEYW,VALUE,ICHPAR,LDEG
          NL=MIN(NL+1,NLM)
          WRITE(DEF(NL),'(A8,''('',1P,E22.15,'','',I3,A8)')
     1      LKEYW,VALUE,ICHPAR,LDEG
          PHI=VALUE*PI/180.0D0
        ENDIF
      ELSE IF(LKEYW.EQ.LXSH) THEN
        WRITE(IW,'(A8,''('',1P,E22.15,'','',I3,A8)')
     1    LKEYW,VALUE,ICHPAR,LOPEN
        NL=MIN(NL+1,NLM)
        WRITE(DEF(NL),'(A8,''('',1P,E22.15,'','',I3,A8)')
     1    LKEYW,VALUE,ICHPAR,LOPEN
        XSHIFT=VALUE
      ELSE IF(LKEYW.EQ.LYSH) THEN
        WRITE(IW,'(A8,''('',1P,E22.15,'','',I3,A8)')
     1    LKEYW,VALUE,ICHPAR,LOPEN
        NL=MIN(NL+1,NLM)
        WRITE(DEF(NL),'(A8,''('',1P,E22.15,'','',I3,A8)')
     1    LKEYW,VALUE,ICHPAR,LOPEN
        YSHIFT=VALUE
      ELSE IF(LKEYW.EQ.LZSH) THEN
        WRITE(IW,'(A8,''('',1P,E22.15,'','',I3,A8)')
     1    LKEYW,VALUE,ICHPAR,LOPEN
        NL=MIN(NL+1,NLM)
        WRITE(DEF(NL),'(A8,''('',1P,E22.15,'','',I3,A8)')
     1    LKEYW,VALUE,ICHPAR,LOPEN
        ZSHIFT=VALUE
      ELSE
        WRITE(IW,'(A72)') BLINE
        WRITE(IW,'(''*** What do you mean?'')')
        STOP
      ENDIF
      GOTO 309
C
C  ****  Rotation and translation of the module contents (its surfaces).
C
 310  CONTINUE
      DO KS=1,NSURF
        KM(KS)=0
      ENDDO
      DO KB=1,NBODY
        KBMOTH=KB  ! We need to transform all the descendants.
 311    CONTINUE
        IF(KBMOTH.EQ.NBODY) THEN
          KNS=KSURF(KB,NX)
          DO KSS=1,KNS
            KS=KSURF(KB,KSS)
            IF(KFLAG(KB,KSS).LT.5.AND.KM(KS).EQ.0) THEN
              QXX=AXX(KS)
              QXY=AXY(KS)
              QXZ=AXZ(KS)
              QYY=AYY(KS)
              QYZ=AYZ(KS)
              QZZ=AZZ(KS)
              QX=AX(KS)
              QY=AY(KS)
              QZ=AZ(KS)
              Q0=A0(KS)
              CALL ROTSHF(OMEGA,THETA,PHI,XSHIFT,YSHIFT,ZSHIFT,
     1                    QXX,QXY,QXZ,QYY,QYZ,QZZ,QX,QY,QZ,Q0)
              AXX(KS)=QXX
              AXY(KS)=QXY
              AXZ(KS)=QXZ
              AYY(KS)=QYY
              AYZ(KS)=QYZ
              AZZ(KS)=QZZ
              AX(KS)=QX
              AY(KS)=QY
              AZ(KS)=QZ
              A0(KS)=Q0
              KM(KS)=1
            ENDIF
          ENDDO
        ELSE
         KBMOTH=KMOTH(KBMOTH)  ! Moves a generation up (grandmother).
         IF(KBMOTH.GT.0) GOTO 311
        ENDIF
      ENDDO
C
 312  CONTINUE
      WRITE(IW,'(64(''0''))')
      NL=MIN(NL+1,NLM)
      WRITE(DEF(NL),'(64(''0''))')
C
      NLB(NBODY)=NL
      KBOMO(NBODY)=1
      DO I=1,NL
        DEFB(NBODY,I)=DEF(I)
C       WRITE(IW,'(''C '',A72)') DEFB(NBODY,I)
      ENDDO
      GOTO 2
C
C  ************  Clone a module.
C
 400  CONTINUE
      ICLONE=ICLONE+1
      READ(BLINE,'(9X,I4,7A8))') KALIAC,(LARRAY(I),I=1,7)
      WRITE(IW,'(''C '',A72)') BLINE
      IF(NBODY.GT.0) THEN
        DO KB0=1,NBODY
        IF(KALIAC.EQ.KALIAB(KB0)) THEN
          WRITE(IW,'(A72)') BLINE
          WRITE(IW,
     1      '(''*** Same label for two bodies or modules.'')')
          STOP
        ENDIF
        ENDDO
      ENDIF
C
 401  CONTINUE
      READ(IR,'(A72)') BLINE
      IF(BLINE(1:2).EQ.'C '.OR.BLINE(1:2).EQ.'c ') THEN
        WRITE(IW,'(A72)') BLINE
        GOTO 401
      ENDIF
      READ(BLINE,'(A8,1X,I4)') LKEYW,KKK
      WRITE(IW,'(''C '',A8,''('',I4,'')'')') LKEYW,KKK
      IF(LKEYW.NE.LMOD) THEN
        WRITE(IW,'(A72)') BLINE
        WRITE(IW,'(''*** The cloned object must be a module.'')')
        STOP
      ENDIF
      IF(NBODY.EQ.0) THEN
        WRITE(IW,'(A72)') BLINE
        WRITE(IW,'(''*** This module is not defined.'')')
        STOP
      ENDIF
      DO KB0=1,NBODY
        IF(KKK.EQ.KALIAB(KB0)) THEN
          KORIG=KB0
          IF(KBOMO(KORIG).NE.1) THEN
            WRITE(IW,'(''*** The cloned object must be a module.'')')
            WRITE(IW,'(''*** The selected object is a body.'')')
            STOP
          ENDIF
          GOTO 402
        ENDIF
      ENDDO
      WRITE(IW,'(''*** The label does not correspond to a module.'')')
      STOP
C
 402  CONTINUE
      READ(IR,'(A72)') BLINE
      IF(BLINE(1:2).EQ.'C '.OR.BLINE(1:2).EQ.'c ') THEN
        WRITE(IW,'(''C '',A72)') BLINE
        GOTO 402
      ENDIF
      READ(BLINE,'(A8,1X,I4)') LKEYW
      IF(LKEYW.EQ.LONE.OR.LKEYW.EQ.LNUL) THEN
        GOTO 403
      ELSE
        WRITE(IW,'(A72)') BLINE
        WRITE(IW,'(''*** What do you mean?'')')
        STOP
      ENDIF
C
C  ****  Transformation parameters.
C
 403  CONTINUE
      OMEGA=0.0D0
      THETA=0.0D0
      PHI=0.0D0
      XSHIFT=0.0D0
      YSHIFT=0.0D0
      ZSHIFT=0.0D0
      IF(LKEYW.EQ.LNUL) GOTO 405
C
 404  CONTINUE
      READ(IR,'(A72)') BLINE
      IF(BLINE(1:2).EQ.'C '.OR.BLINE(1:2).EQ.'c ') THEN
        WRITE(IW,'(''C '',A72)') BLINE
        GOTO 404
      ENDIF
      READ(BLINE,'(A8,1X,E22.15,1X,I3,A8)') LKEYW,VALUE,ICHPAR,LANGLE
      IF(LKEYW.EQ.LNUL) GOTO 405
      IF(ICHPAR.LE.NPINP) THEN
        IF(ICHPAR.GT.0) THEN
          VALUE=PARINP(ICHPAR)
          ICHPAR=-ICHPAR
        ENDIF
      ELSE
        WRITE(IW,'(A72)') BLINE
        WRITE(IW,'(''*** NPINP is too small (check PARINP).'')')
        STOP
      ENDIF
      IF(LKEYW.EQ.LOME) THEN
        IF(LANGLE.EQ.LRAD) THEN
          WRITE(IW,'(''C '',A8,''('',1P,E22.15,'','',I3,A8)')
     1      LKEYW,VALUE,ICHPAR,LRAD
          OMEGA=VALUE
        ELSE
          WRITE(IW,'(''C '',A8,''('',1P,E22.15,'','',I3,A8)')
     1      LKEYW,VALUE,ICHPAR,LDEG
          OMEGA=VALUE*PI/180.0D0
        ENDIF
      ELSE IF(LKEYW.EQ.LTHE) THEN
        IF(LANGLE.EQ.LRAD) THEN
          WRITE(IW,'(''C '',A8,''('',1P,E22.15,'','',I3,A8)')
     1      LKEYW,VALUE,ICHPAR,LRAD
          THETA=VALUE
        ELSE
          WRITE(IW,'(''C '',A8,''('',1P,E22.15,'','',I3,A8)')
     1      LKEYW,VALUE,ICHPAR,LDEG
          THETA=VALUE*PI/180.0D0
        ENDIF
      ELSE IF(LKEYW.EQ.LPHI) THEN
        IF(LANGLE.EQ.LRAD) THEN
          WRITE(IW,'(''C '',A8,''('',1P,E22.15,'','',I3,A8)')
     1      LKEYW,VALUE,ICHPAR,LRAD
          PHI=VALUE
        ELSE
          WRITE(IW,'(''C '',A8,''('',1P,E22.15,'','',I3,A8)')
     1      LKEYW,VALUE,ICHPAR,LDEG
          PHI=VALUE*PI/180.0D0
        ENDIF
      ELSE IF(LKEYW.EQ.LXSH) THEN
        WRITE(IW,'(''C '',A8,''('',1P,E22.15,'','',I3,A8)')
     1    LKEYW,VALUE,ICHPAR,LOPEN
        XSHIFT=VALUE
      ELSE IF(LKEYW.EQ.LYSH) THEN
        WRITE(IW,'(''C '',A8,''('',1P,E22.15,'','',I3,A8)')
     1    LKEYW,VALUE,ICHPAR,LOPEN
        YSHIFT=VALUE
      ELSE IF(LKEYW.EQ.LZSH) THEN
        WRITE(IW,'(''C '',A8,''('',1P,E22.15,'','',I3,A8)')
     1    LKEYW,VALUE,ICHPAR,LOPEN
        ZSHIFT=VALUE
      ELSE
        WRITE(IW,'(A72)') BLINE
        WRITE(IW,'(''*** What do you mean?'')')
        STOP
      ENDIF
      GOTO 404
 405  CONTINUE
C  ****  Determine all the descendants of module KORIG.
      ND=1
      IDESC(1)=KORIG  ! Descendants of the cloned module.
      IDONE(1)=0      ! The descendants have not yet been identified.
 406  CONTINUE
      NDC=ND
      KDG=0
      DO I=1,NDC
        IF(IDONE(I).EQ.0) THEN
          KB=IDESC(I)
          IF(KBODY(KB,NX).GT.0) THEN
            DO J=1,KBODY(KB,NX)
              ND=ND+1
              IDESC(ND)=KBODY(KB,J)  ! New descendant.
              IDONE(ND)=0
              KDG=KDG+1
            ENDDO
          ELSE IF(KDGHT(KB,NX).GT.0) THEN
            DO J=1,KDGHT(KB,NX)
              IF(KDGHT(KB,J).NE.KB) THEN
                ND=ND+1
                IDESC(ND)=KDGHT(KB,J)  ! New descendant.
                IDONE(ND)=0
                KDG=KDG+1
              ENDIF
            ENDDO
          ENDIF
          IDONE(I)=1  ! The descendants of KB=IDEC(I) have been listed.
        ENDIF
      ENDDO
      IF(KDG.GT.0) GOTO 406
C
      IN=0
      DO I=1,NB
        IBCL(I)=0  ! Label of a cloned body or module.
        IBOR(I)=0  ! Label of the original body or module.
      ENDDO
      KSD=NSURF
      DO I=1,NS
        ISCL(I)=0  ! Label of a cloned surface.
        ISOR(I)=0  ! Label of the original surface.
      ENDDO
      KBD=NBODY
      DO 407 KBB=1,NBODY
        DO ID=ND,1,-1
          KB=IDESC(ID)
          IF(KBB.EQ.KB) THEN
            KBD=KBD+1
            IBCL(KB)=KBD
            IBOR(KBD)=KB
            MATER(KBD)=MATER(KB)
C  ****  Clone the surfaces of the original module and its descendants.
            KSURF(KBD,NX)=KSURF(KB,NX)
            DO KSS=KSURF(KB,NX),1,-1
              KS=KSURF(KB,KSS)
              IF(KFLAG(KB,KSS).LT.3.AND.ISCL(KS).EQ.0) THEN
                KSD=KSD+1
                IF(KSD.GT.NS) THEN
                  WRITE(IW,
     1              '(''*** The parameter NS must be increased.'')')
                  STOP
                ENDIF
                ISCL(KS)=KSD
                ISOR(KSD)=KS
                QXX=AXX(KS)
                QXY=AXY(KS)
                QXZ=AXZ(KS)
                QYY=AYY(KS)
                QYZ=AYZ(KS)
                QZZ=AZZ(KS)
                QX=AX(KS)
                QY=AY(KS)
                QZ=AZ(KS)
                Q0=A0(KS)
                CALL ROTSHF(OMEGA,THETA,PHI,XSHIFT,YSHIFT,ZSHIFT,
     1                      QXX,QXY,QXZ,QYY,QYZ,QZZ,QX,QY,QZ,Q0)
                AXX(KSD)=QXX
                AXY(KSD)=QXY
                AXZ(KSD)=QXZ
                AYY(KSD)=QYY
                AYZ(KSD)=QYZ
                AZZ(KSD)=QZZ
                AX(KSD)=QX
                AY(KSD)=QY
                AZ(KSD)=QZ
                A0(KSD)=Q0
                WRITE(IW,'(A8,''('',I4,'')'')') LSUR,KSD
                WRITE(IW,'(A8,''('',4(I2,'',''),I2,'')'')')
     1            LIND,IN,IN,IN,IN,IN
                IF(ABS(AXX(KSD)).GT.1.0D-20)
     1            WRITE(IW,'(A8,''('',1P,E22.15,'','',I3,A8)')
     2            LAXX,AXX(KSD),IN,LOPEN
                IF(ABS(AXY(KSD)).GT.1.0D-20)
     1            WRITE(IW,'(A8,''('',1P,E22.15,'','',I3,A8)')
     2            LAXY,AXY(KSD),IN,LOPEN
                IF(ABS(AXZ(KSD)).GT.1.0D-20)
     1            WRITE(IW,'(A8,''('',1P,E22.15,'','',I3,A8)')
     2            LAXZ,AXZ(KSD),IN,LOPEN
                IF(ABS(AYY(KSD)).GT.1.0D-20)
     1            WRITE(IW,'(A8,''('',1P,E22.15,'','',I3,A8)')
     2            LAYY,AYY(KSD),IN,LOPEN
                IF(ABS(AYZ(KSD)).GT.1.0D-20)
     1            WRITE(IW,'(A8,''('',1P,E22.15,'','',I3,A8)')
     2            LAYZ,AYZ(KSD),IN,LOPEN
                IF(ABS(AZZ(KSD)).GT.1.0D-20)
     1            WRITE(IW,'(A8,''('',1P,E22.15,'','',I3,A8)')
     2            LAZZ,AZZ(KSD),IN,LOPEN
                IF(ABS(AX(KSD)).GT.1.0D-20)
     1            WRITE(IW,'(A8,''('',1P,E22.15,'','',I3,A8)')
     2            LAX,AX(KSD),IN,LOPEN
                IF(ABS(AY(KSD)).GT.1.0D-20)
     1            WRITE(IW,'(A8,''('',1P,E22.15,'','',I3,A8)')
     2            LAY,AY(KSD),IN,LOPEN
                IF(ABS(AZ(KSD)).GT.1.0D-20)
     1            WRITE(IW,'(A8,''('',1P,E22.15,'','',I3,A8)')
     2            LAZ,AZ(KSD),IN,LOPEN
                IF(ABS(A0(KSD)).GT.1.0D-20)
     1            WRITE(IW,'(A8,''('',1P,E22.15,'','',I3,A8)')
     2            LA0,A0(KSD),IN,LOPEN
                  WRITE(IW,'(64(''0''))')
                  KALIAS(KSD)=KALIAS(KS)+ICLONE*10000
              ENDIF
            ENDDO
            GOTO 407
          ENDIF
        ENDDO
 407  CONTINUE
C  ****  Clone the original module and its descendants.
      DO KB=NBODY+1,KBD
        IF(KB.GT.NB) THEN
          WRITE(IW,'(''*** The parameter NB must be increased.'')')
          STOP
        ENDIF
        KBO=IBOR(KB)
C  ****  The new element is assigned the alias of the 'father' plus
C  the ICLONE value times 10000 (this allows keeping track of the
C  cloning generation).
        KALIAB(KB)=KALIAB(KBO)+ICLONE*10000
        IF(KMOTH(KBO).GT.0) THEN
          KMOTH(KB)=IBCL(KMOTH(KBO))
        ELSE
          KMOTH(KB)=0
        ENDIF
        KBOMO(KB)=KBOMO(KBO)
        DO I=1,NLM
          DEFB(KB,I)=DEFB(KBO,I)
        ENDDO
        BLINE=DEFB(KBO,1)
        IF(KBOMO(KB).EQ.0) THEN
          LKEYW=LBOD
          WRITE(IW,'(A8,''('',I4,A51)') LBOD,KB,BLINE(14:64)
        ELSE IF(KBOMO(KB).EQ.1) THEN
          LKEYW=LMOD
          WRITE(IW,'(A8,''('',I4,A51)') LMOD,KB,BLINE(14:64)
        ELSE
          WRITE(IW,'(''KBOMO('',I4,'') ='',I4)') KB,KBOMO(KB)
          WRITE(IW,'(''*** Something wrong...'')')
          STOP
        ENDIF
        WRITE(IW,'(A8,''('',I4,'')'')') LMAT,MATER(KB)
        DO KS=1,KSURF(KB,NX)
          KSURF(KB,KS)=ISCL(KSURF(KBO,KS))
          KFLAG(KB,KS)=KFLAG(KBO,KS)
          IF(KFLAG(KB,KS).LT.3) THEN
            IF(KFLAG(KB,KS).EQ.1) THEN
              INDS=-1
            ELSE
              INDS=+1
            ENDIF
            WRITE(IW,'(A8,''('',I4,''), SIDE POINTER=('',I2,'')'')')
     1        LSUR,KSURF(KB,KS),INDS
          ENDIF
        ENDDO
        IF(KBOMO(KB).EQ.0) THEN
          KBODY(KB,NX)=KBODY(KBO,NX)
          DO I=1,KBODY(KB,NX)
            KBODY(KB,I)=IBCL(KBODY(KBO,I))
            KBB=KBODY(KB,I)
            IF(KBOMO(KB).EQ.0) THEN
              LKEYW=LBOD
              WRITE(IW,'(A8,''('',I4,'')'')') LBOD,KBB
            ELSE IF(KBB.NE.KB) THEN
              LKEYW=LMOD
              WRITE(IW,'(A8,''('',I4,'')'')') LMOD,KBB
            ENDIF
            IF(KBB.GE.KB) THEN
              WRITE(IW,'(''*** The limiting body or module is not '',
     1          ''yet defined'')')
              STOP
            ENDIF
          ENDDO
        ELSE
          KDGHT(KB,NX)=KDGHT(KBO,NX)
          DO I=1,KDGHT(KB,NX)
            KDGHT(KB,I)=IBCL(KDGHT(KBO,I))
            KBB=KDGHT(KB,I)
            IF(KBB.NE.KB) THEN
              IF(KBOMO(KBB).EQ.0) THEN
                LKEYW=LBOD
                WRITE(IW,'(A8,''('',I4,'')'')') LBOD,KBB
              ELSE
                LKEYW=LMOD
                WRITE(IW,'(A8,''('',I4,'')'')') LMOD,KBB
              ENDIF
              IF(KBB.GE.KB) THEN
                WRITE(IW,'(''*** The limiting body or module is not '',
     1            ''yet defined'')')
                STOP
              ENDIF
            ENDIF
          ENDDO
        ENDIF
        WRITE(IW,'(64(''0''))')
      ENDDO
      KALIAB(KBD)=KALIAC
      KMOTH(KBD)=0
      NBODY=KBD
      NSURF=KSD
      GOTO 2
C
C  ************  Included geometry file.
C
 500  CONTINUE
      READ(IR,'(A72)') BLINE
      IF(BLINE(1:2).EQ.'C '.OR.BLINE(1:2).EQ.'c ') THEN
        WRITE(IW,'(A72)') BLINE
        GOTO 500
      ENDIF
      READ(BLINE,'(A8,1X,A12,6A8)') LKEYW,GFILE,(LARRAY(I),I=1,6)
      IF(LKEYW.NE.LFIL) THEN
        WRITE(IW,'(A72)') BLINE
        WRITE(IW,'(''*** What do you mean?'')')
        STOP
      ENDIF
 501  CONTINUE
      READ(IR,'(A72)') BLINE
      IF(BLINE(1:2).EQ.'C '.OR.BLINE(1:2).EQ.'c ') THEN
        WRITE(IW,'(A72)') BLINE
        GOTO 501
      ENDIF
      READ(BLINE,'(8A8)') LKEYW,(LARRAY(I),I=1,6)
      IF(LKEYW.EQ.LNUL) THEN
        IF(IR.EQ.IRI) THEN
          WRITE(IW,'(A72)') BLINE
          WRITE(IW,'(''*** Too many include levels.'')')
          STOP
        ELSE
          IR=IRI
          OPEN(IR,FILE=GFILE)
        ENDIF
        GOTO 1
      ENDIF
C
C  ************  End-line in the input file.
C
 600  CONTINUE
      IF(IR.EQ.IRI) THEN
        CLOSE(IR)
        IR=IRD
        GOTO 2
      ENDIF
      IF(NBODY.EQ.1) THEN
        KDGHT(1,NX)=1
        KDGHT(1,1)=1
      ENDIF
C
C  ****  Define the enclosure, if it has not been defined yet.
C
      NORPH=0
      DO KB=1,NBODY
        IF(KMOTH(KB).EQ.0) NORPH=NORPH+1
      ENDDO
C  ****  The next line serves only to avoid a warning message issued
C        by certain compilers.
      NT=0
      IF(NORPH.GT.1) THEN
        KS=NSURF+1
        AXX(KS)=1.0D0
        AYY(KS)=1.0D0
        AZZ(KS)=1.0D0
        A0(KS)=-1.0D30
        KB=NBODY+1
        KBOMO(KB)=1
        KSURF(KB,NX)=1
        KSURF(KB,1)=KS
        KFLAG(KB,1)=1
        KDGHT(KB,NX)=1
        KDGHT(KB,1)=KB
        DO KBB=1,NBODY
          IF(KMOTH(KBB).EQ.0) THEN
            NT=KDGHT(KB,NX)+1
            KDGHT(KB,NX)=NT
            KDGHT(KB,NT)=KBB
            KN1=KSURF(KBB,NX)
            DO 601 KS1=1,KN1
              IF(KFLAG(KBB,KS1).GT.3) GOTO 601
              KSURF1=KSURF(KBB,KS1)
              KN2=KSURF(KB,NX)
              DO KS2=1,KN2
                IF(KSURF(KB,KS2).EQ.KSURF1) GOTO 601
              ENDDO
              KN2=KN2+1
              IF(KN2.GE.NX) THEN
                WRITE(IW,'(''*** The parameter NX is too small.'')')
                STOP
              ENDIF
              KSURF(KB,NX)=KN2
              KSURF(KB,KN2)=KSURF1
              KFLAG(KB,KN2)=4
 601        CONTINUE
            KMOTH(KBB)=KB
          ENDIF
        ENDDO
        NSURF=KS
        NBODY=KB
        IF(NBODY.GT.NB) THEN
          WRITE(IW,'(''*** The parameter NB must be increased.'')')
          STOP
        ENDIF
        IF(NSURF.GT.NS) THEN
          WRITE(IW,'(''*** The parameter NS must be increased.'')')
          STOP
        ENDIF
C  ****  Sort daughters in increasing order.
        IF(NT.GT.1) THEN
          DO KI=1,NT-1
            KBMIN=KDGHT(NBODY,KI)
            KMIN=KI
            DO KJ=KI+1,NT
              IF(KDGHT(NBODY,KJ).LT.KBMIN) THEN
                KBMIN=KDGHT(NBODY,KJ)
                KMIN=KJ
              ENDIF
            ENDDO
            IF(KMIN.NE.KI) THEN
              KSAVE=KDGHT(NBODY,KI)
              KDGHT(NBODY,KI)=KDGHT(NBODY,KMIN)
              KDGHT(NBODY,KMIN)=KSAVE
            ENDIF
          ENDDO
        ENDIF
      ENDIF
C
      WRITE(IW,'(A8,1X,55(''0''))') LKEYW
      NBOD=NBODY
      NS2=2*NSURF
C
C  ****  Duplicated surfaces (within round-off tolerance) are removed.
C
      IF(IREP.GT.0) THEN
        OPEN(IREP,FILE='pengeom_tree.rep')
        IF(NSURF.LT.3) RETURN
        WRITE(IREP,'(''C'')')
        WRITE(IREP,'(''C *****************************************'')')
        WRITE(IREP,'(''C ****     PENGEOM (version 2005)      ****'')')
        WRITE(IREP,'(''C ****  Constructive Quadric Geometry  ****'')')
        WRITE(IREP,'(''C *****************************************'')')
        WRITE(IREP,'(''C'')')
      ENDIF
      IWRITE=0
      TOL=1.0D-12
      DO KS=1,NSURF
        KM(KS)=0
      ENDDO
      DO 703 KS=1,NSURF-1
        IF(KM(KS).NE.0) GOTO 703
        F=MAX(AXX(KS),AXY(KS),AXZ(KS),AYY(KS),AYZ(KS),AZZ(KS),
     1        AX(KS),AY(KS),AZ(KS),A0(KS))
        FM=MIN(AXX(KS),AXY(KS),AXZ(KS),AYY(KS),AYZ(KS),AZZ(KS),
     1         AX(KS),AY(KS),AZ(KS),A0(KS))
        IF(ABS(FM).GT.ABS(F)) F=FM
        DO 702 KSP=KS+1,NSURF
          IF(KM(KSP).NE.0) GOTO 702
          FP=MAX(AXX(KSP),AXY(KSP),AXZ(KSP),AYY(KSP),AYZ(KSP),AZZ(KSP),
     1           AX(KSP),AY(KSP),AZ(KSP),A0(KSP))
          FM=MIN(AXX(KSP),AXY(KSP),AXZ(KSP),AYY(KSP),AYZ(KSP),AZZ(KSP),
     1           AX(KSP),AY(KSP),AZ(KSP),A0(KSP))
          IF(ABS(FM).GT.ABS(FP)) FP=FM
C
          IF(ABS(AXX(KS)*FP-AXX(KSP)*F).LE.TOL*ABS(AXX(KS)*FP).AND.
     1       ABS(AXY(KS)*FP-AXY(KSP)*F).LE.TOL*ABS(AXY(KS)*FP).AND.
     2       ABS(AXZ(KS)*FP-AXZ(KSP)*F).LE.TOL*ABS(AXZ(KS)*FP).AND.
     3       ABS(AYY(KS)*FP-AYY(KSP)*F).LE.TOL*ABS(AYY(KS)*FP).AND.
     4       ABS(AYZ(KS)*FP-AYZ(KSP)*F).LE.TOL*ABS(AYZ(KS)*FP).AND.
     5       ABS(AZZ(KS)*FP-AZZ(KSP)*F).LE.TOL*ABS(AZZ(KS)*FP).AND.
     6       ABS(AX(KS)*FP-AX(KSP)*F).LE.TOL*ABS(AX(KS)*FP).AND.
     7       ABS(AY(KS)*FP-AY(KSP)*F).LE.TOL*ABS(AY(KS)*FP).AND.
     8       ABS(AZ(KS)*FP-AZ(KSP)*F).LE.TOL*ABS(AZ(KS)*FP).AND.
     9       ABS(A0(KS)*FP-A0(KSP)*F).LE.TOL*ABS(A0(KS)*FP)) THEN
            IF(IREP.GT.0) THEN
              IF(IWRITE.EQ.0) THEN
                WRITE(IREP,'(''C'')')
                WRITE(IREP,'(''C'')')
                WRITE(IREP,'(''C ************  Modifications of '',
     1            ''user-defined surfaces.'')')
                WRITE(IREP,'(''C'')')
                IWRITE=1
              ENDIF
              WRITE(IREP,'(''C SURFACE ('',I4,'') has been replaced ''
     1          ''by SURFACE ('',I4,'')'' )') KSP,KS
            ENDIF
C  ****  Check whether the two surface functions have the same global
C        sign (ISPF=0) or not (ISPF=1).
            ISPF=0
            IF(ABS(AXX(KS)).GT.TOL.AND.AXX(KS)*AXX(KSP).LT.0.0D0) ISPF=1
            IF(ABS(AXY(KS)).GT.TOL.AND.AXY(KS)*AXY(KSP).LT.0.0D0) ISPF=1
            IF(ABS(AXZ(KS)).GT.TOL.AND.AXZ(KS)*AXZ(KSP).LT.0.0D0) ISPF=1
            IF(ABS(AYY(KS)).GT.TOL.AND.AYY(KS)*AYY(KSP).LT.0.0D0) ISPF=1
            IF(ABS(AYZ(KS)).GT.TOL.AND.AYZ(KS)*AYZ(KSP).LT.0.0D0) ISPF=1
            IF(ABS(AZZ(KS)).GT.TOL.AND.AZZ(KS)*AZZ(KSP).LT.0.0D0) ISPF=1
            IF(ABS(AX(KS)).GT.TOL.AND.AX(KS)*AX(KSP).LT.0.0D0) ISPF=1
            IF(ABS(AY(KS)).GT.TOL.AND.AY(KS)*AY(KSP).LT.0.0D0) ISPF=1
            IF(ABS(AZ(KS)).GT.TOL.AND.AZ(KS)*AZ(KSP).LT.0.0D0) ISPF=1
            IF(ABS(A0(KS)).GT.TOL.AND.A0(KS)*A0(KSP).LT.0.0D0) ISPF=1
C
            AXX(KSP)=AXX(KS)
            AXY(KSP)=AXY(KS)
            AXZ(KSP)=AXZ(KS)
            AYY(KSP)=AYY(KS)
            AYZ(KSP)=AYZ(KS)
            AZZ(KSP)=AZZ(KS)
            AX(KSP)=AX(KS)
            AY(KSP)=AY(KS)
            AZ(KSP)=AZ(KS)
            A0(KSP)=A0(KS)
            KPLANE(KSP)=KPLANE(KS)
            KM(KSP)=KS
C
            KBI=1
 701        CONTINUE
            DO KB=KBI,NBODY
              KSL=0
              KSLP=0
              NSB=KSURF(KB,NX)
              DO K=1,NSB
                IF(KSURF(KB,K).EQ.KS) KSL=K
                IF(KSURF(KB,K).EQ.KSP) KSLP=K
              ENDDO
              IF(KSLP.GT.0) THEN
                KSURF(KB,KSLP)=KS
C  ****  If the implicit equations of surfaces KS and KSP differ by a
C        global sign, the side pointer of KSP must be reversed.
                IF(ISPF.EQ.1) THEN
                  IF(KFLAG(KB,KSLP).EQ.1) THEN
                    KFLAG(KB,KSLP)=2
                  ELSE IF(KFLAG(KB,KSLP).EQ.2) THEN
                    KFLAG(KB,KSLP)=1
                  ENDIF
                ENDIF
                IF(KSL.GT.0) THEN
                  KFL=KFLAG(KB,KSL)
                  KFLP=KFLAG(KB,KSLP)
                  IF(MIN(KFL,KFLP).LT.3) THEN
                     KFLAG(KB,KSL)=MIN(KFL,KFLP)
                  ELSE IF(KBOMO(KB).EQ.0) THEN
                     KFLAG(KB,KSL)=3
                  ELSE IF(KBOMO(KB).EQ.1) THEN
                     KFLAG(KB,KSL)=4
                  ELSE
                     KFLAG(KB,KSL)=5
                  ENDIF
                  DO K=KSLP,NSB-1
                    KSURF(KB,K)=KSURF(KB,K+1)
                    KFLAG(KB,K)=KFLAG(KB,K+1)
                  ENDDO
                  KSURF(KB,NSB)=0
                  KFLAG(KB,NSB)=5
                  KSURF(KB,NX)=NSB-1
                  IF(KB.LT.NBODY) THEN
                    KBI=KB
                    GOTO 701
                  ENDIF
                ENDIF
              ENDIF
            ENDDO
          ENDIF
 702    CONTINUE
 703  CONTINUE
C
      IF(IREP.LE.0) RETURN
      WRITE(IREP,'(''C'')')
      WRITE(IREP,'(''C'')')
      WRITE(IREP,'(''C ************  Genealogical tree. '')')
      WRITE(IREP,'(''C'')')
      DO KB=1,NBODY
        IF(KBOMO(KB).EQ.0) THEN
          WRITE(IREP,'(''C '',/''C *** BODY   ='',I5,'',  KMOTH ='',
     1      I5,'',  MAT ='',I3)') KB,KMOTH(KB),MATER(KB)
          IF(KBODY(KB,NX).GT.0)
     1    WRITE(IREP,'(''C KBODY ='',15I5,100(/''C '',7X,15I5))')
     2      (KBODY(KB,K2),K2=1,KBODY(KB,NX))
        ELSE IF(KBOMO(KB).EQ.1) THEN
          WRITE(IREP,'(''C '',/''C *** MODULE ='',I5,'',  KMOTH ='',
     1      I5,'',  MAT ='',I3)') KB,KMOTH(KB),MATER(KB)
          WRITE(IREP,'(''C KDGHT ='',15I5,100(/''C '',7X,15I5))')
     1      (KDGHT(KB,K2),K2=1,KDGHT(KB,NX))
        ELSE
          WRITE(IREP,'(''C '',/''C ****** ERROR: the label '',I5,
     1      '' does not correspond to a body.'')') KB
          WRITE(IW,'(''C '',/''C ****** ERROR: the label '',I5,
     1      '' does not correspond to a body.'')') KB
          STOP
        ENDIF
        WRITE(IREP,'(''C KSURF ='',15I5,100(/''C '',7X,15I5))')
     1    (KSURF(KB,KS),KS=1,KSURF(KB,NX))
        WRITE(IREP,'(''C KFLAG ='',15I5,100(/''C '',7X,15I5))')
     1    (KFLAG(KB,KS),KS=1,KSURF(KB,NX))
      ENDDO
C
C  ****  Easiness test.
C
      NBU=0
      NSU=0
      DO KB=1,NBODY
        IF(KBOMO(KB).EQ.0) THEN
          NBU=MAX(NBU,KBODY(KB,NX))
        ELSE IF(KBOMO(KB).EQ.1) THEN
          NBU=MAX(NBU,KDGHT(KB,NX))
        ENDIF
        NSE=0
        DO K=1,KSURF(KB,NX)
          IF(KFLAG(KB,K).LT.5) NSE=NSE+1
        ENDDO
        NSU=MAX(NSU,NSE)
      ENDDO
      WRITE(IREP,'(''C'')')
      WRITE(IREP,'(''C'')')
      WRITE(IREP,'(''C ************  Adequacy of the geometry defin'',
     1  ''ition.'')')
      WRITE(IREP,'(''C'')')
      WRITE(IREP,'(''C The largest number of bodies in a module or'')')
      WRITE(IREP,'(''C     bodies limiting a single body is ........'',
     1  ''.... '',I4)') NBU
      WRITE(IREP,'(''C The largest number of limiting surfaces for'')')
      WRITE(IREP,'(''C     a single body or module is ...............'',
     1  ''... '',I4)') NSU
      WRITE(IREP,'(''C'')')
      IF(NBODY.LT.15.AND.NSURF.LT.15) THEN
        WRITE(IREP,'(''C The simulation of this geometry will be rela'',
     1    ''tively fast,'')')
        WRITE(IREP,'(''C     no further optimization seems to be re'',
     1    ''quired.'')')
      ELSE IF(NBU.LT.10.AND.NSU.LT.10) THEN
        WRITE(IREP,'(''C The simulation of this geometry will be rela'',
     1    ''tively fast,'')')
        WRITE(IREP,'(''C     no further optimization seems to be re'',
     1    ''quired.'')')
      ELSE IF(NBU.LT.15.AND.NSU.LT.20) THEN
        WRITE(IREP,'(''C The simulation of this geometry is expected '',
     1    ''to be slow,'')')
        WRITE(IREP,'(''C     try to split complex bodies into seve'',
     1    ''ral modules.'')')
      ELSE IF(NBU.LT.25.AND.NSU.LT.30) THEN
        WRITE(IREP,'(''C The simulation of this geometry will be ve'',
     1    ''ry slow,'')')
        WRITE(IREP,'(''C     you should optimize the structure of '',
     1    ''the module tree.'')')
      ELSE
        WRITE(IREP,'(''C Simulating this geometry will be extremely '',
     1    ''slow.'')')
      ENDIF
      WRITE(IREP,'(''C'')')
      WRITE(IREP,'(''C ************  The end.'')')
      CLOSE(IREP)
C
      RETURN
      END
C  *********************************************************************
C                       SUBROUTINE ROTSHF
C  *********************************************************************
      SUBROUTINE ROTSHF(OMEGA,THETA,PHI,DX,DY,DZ,
     1                  AXX,AXY,AXZ,AYY,AYZ,AZZ,AX,AY,AZ,A0)
C
C     This subroutine rotates and shifts a quadric surface.
C
C  Input parameters:
C     OMEGA, THETA, PHI ... Euler rotation angles,
C     DX, DY, DZ .......... components of the displacement vector,
C     AXX, ..., A0 ........ coefficients of the initial quadric.
C
C  Output parameters:
C     AXX, ..., A0 ........ coefficients of the transformed quadric.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N)
      DIMENSION R(3,3),A2(3,3),B2(3,3),A1(3),B1(3),D1(3)
C
C  ****  Initial quadric.
C
      B2(1,1)=AXX
      B2(1,2)=0.5D0*AXY
      B2(1,3)=0.5D0*AXZ
      B2(2,1)=B2(1,2)
      B2(2,2)=AYY
      B2(2,3)=0.5D0*AYZ
      B2(3,1)=B2(1,3)
      B2(3,2)=B2(2,3)
      B2(3,3)=AZZ
      B1(1)=AX
      B1(2)=AY
      B1(3)=AZ
      B0=A0
      D1(1)=DX
      D1(2)=DY
      D1(3)=DZ
C
C  ****  Rotation matrix.
C
      STHETA=SIN(THETA)
      CTHETA=COS(THETA)
      SPHI=SIN(PHI)
      CPHI=COS(PHI)
      SOMEGA=SIN(OMEGA)
      COMEGA=COS(OMEGA)
C
      R(1,1)=CPHI*CTHETA*COMEGA-SPHI*SOMEGA
      R(1,2)=-CPHI*CTHETA*SOMEGA-SPHI*COMEGA
      R(1,3)=CPHI*STHETA
      R(2,1)=SPHI*CTHETA*COMEGA+CPHI*SOMEGA
      R(2,2)=-SPHI*CTHETA*SOMEGA+CPHI*COMEGA
      R(2,3)=SPHI*STHETA
      R(3,1)=-STHETA*COMEGA
      R(3,2)=STHETA*SOMEGA
      R(3,3)=CTHETA
C
C  ****  Rotated quadric.
C
      DO I=1,3
        A1(I)=0.0D0
        DO J=1,3
          A1(I)=A1(I)+R(I,J)*B1(J)
          A2(I,J)=0.0D0
          DO M=1,3
            DO K=1,3
              A2(I,J)=A2(I,J)+R(I,K)*B2(K,M)*R(J,M)
            ENDDO
          ENDDO
        ENDDO
      ENDDO
C
C  ****  Shifted-rotated quadric.
C
      DO I=1,3
        A2D=0.0D0
        DO J=1,3
          A2D=A2D+A2(I,J)*D1(J)
        ENDDO
        B1(I)=A1(I)-2.0D0*A2D
        B0=B0+D1(I)*(A2D-A1(I))
      ENDDO
C
      AXX=A2(1,1)
      IF(ABS(AXX).LT.1.0D-15) AXX=0.0D0
      AXY=A2(1,2)+A2(2,1)
      IF(ABS(AXY).LT.1.0D-15) AXY=0.0D0
      AXZ=A2(1,3)+A2(3,1)
      IF(ABS(AXZ).LT.1.0D-15) AXZ=0.0D0
      AYY=A2(2,2)
      IF(ABS(AYY).LT.1.0D-15) AYY=0.0D0
      AYZ=A2(2,3)+A2(3,2)
      IF(ABS(AYZ).LT.1.0D-15) AYZ=0.0D0
      AZZ=A2(3,3)
      IF(ABS(AZZ).LT.1.0D-15) AZZ=0.0D0
      AX=B1(1)
      IF(ABS(AX).LT.1.0D-15) AX=0.0D0
      AY=B1(2)
      IF(ABS(AY).LT.1.0D-15) AY=0.0D0
      AZ=B1(3)
      IF(ABS(AZ).LT.1.0D-15) AZ=0.0D0
      A0=B0
      IF(ABS(A0).LT.1.0D-15) A0=0.0D0
      RETURN
      END
C  *********************************************************************
C                       SUBROUTINE LOCATE
C  *********************************************************************
      SUBROUTINE LOCATE
C
C     This subroutine determines the body that contains the point with
C  coordinates (X,Y,Z).
C
C  Input values (common /TRACK/):
C     X, Y, Z ... coordinates of the particle,
C     U, V, W ... direction of movement.
C
C  Output values (common /TRACK/):
C     IBODY ..... body where the particle moves,
C     MAT  ...... material in IBODY,
C                    MAT=0, indicates a void region.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N)
      PARAMETER (NS=100000,NB=40000,NX=1000)
      PARAMETER (EPS=1.0D-8)
      COMMON/QSURF/AXX(NS),AXY(NS),AXZ(NS),AYY(NS),AYZ(NS),AZZ(NS),
     1    AX(NS),AY(NS),AZ(NS),A0(NS),NSURF,NS2,KPLANE(NS)
      COMMON/QTREE/NBODY,MATER(NB),KMOTH(NB),KDGHT(NB,NX),
     1    KSURF(NB,NX),KFLAG(NB,NX),KSP(NS),KSLAST
C
      COMMON/TRACK/E,X,Y,Z,U,V,W,WGHT,KPAR,IBODY,MAT,ILB(5)
C
      DO I=1,NSURF
        KSP(I)=0
      ENDDO
      KB0=NBODY
 100  CONTINUE
      DO 101 KSS=1,KSURF(KB0,NX)
        KS=KSURF(KB0,KSS)
        IF(KSP(KS).NE.0.AND.KFLAG(KB0,KSS).GT.4) GOTO 101
        IF(KPLANE(KS).EQ.0) THEN
          F=X*(AXX(KS)*X+AXY(KS)*Y+AXZ(KS)*Z+AX(KS))
     1     +Y*(AYY(KS)*Y+AYZ(KS)*Z+AY(KS))
     2     +Z*(AZZ(KS)*Z+AZ(KS))+A0(KS)
        ELSE
          F=X*AX(KS)+Y*AY(KS)+Z*AZ(KS)+A0(KS)
        ENDIF
C
C  ****  Point close to a surface.
C
        IF(ABS(F).LT.1.0D-12) THEN
C  ****  Quadric surface gradient.
          DFX=2.0D0*AXX(KS)*X+AXY(KS)*Y+AXZ(KS)*Z+AX(KS)
          DFY=AXY(KS)*X+2.0D0*AYY(KS)*Y+AYZ(KS)*Z+AY(KS)
          DFZ=AXZ(KS)*X+AYZ(KS)*Y+2.0D0*AZZ(KS)*Z+AZ(KS)
          DF=SQRT(DFX*DFX+DFY*DFY+DFZ*DFZ)
          DOTP=U*DFX+V*DFY+W*DFZ
          IF(DOTP.GT.0.0D0) THEN
            ISN=1
          ELSE
            ISN=-1
          ENDIF
C  ****  Shift the point a small distance along the gradient direction.
          IF(DF.GT.1.0D-12) THEN
            X=X+ISN*EPS*DFX/DF
            Y=Y+ISN*EPS*DFY/DF
            Z=Z+ISN*EPS*DFZ/DF
          ELSE
            X=X+ISN*EPS
            Y=Y+ISN*EPS
            Z=Z+ISN*EPS
          ENDIF
          DO I=1,NSURF
            KSP(I)=0
          ENDDO
          KB0=NBODY
          GOTO 100
        ENDIF
C
        IF(F.LT.0.0D0) THEN
          KSP(KS)=1
        ELSE
          KSP(KS)=2
        ENDIF
 101  CONTINUE
C
C  ****  Determine the module or body that contains the point.
C
      DO 102 KBB=1,KDGHT(KB0,NX)
        KB=KDGHT(KB0,KBB)
        DO KSS=1,KSURF(KB,NX)
          KS=KSURF(KB,KSS)
          KF=KFLAG(KB,KSS)
          IF(KF.LT.3.AND.KSP(KS).NE.KF) GOTO 102
        ENDDO
        IF(KB.EQ.KB0) THEN
          IBODY=KB
          MAT=MATER(KB)
          RETURN
        ELSE IF(KDGHT(KB,NX).GT.1) THEN
          KB0=KB
          GOTO 100
        ELSE
          IBODY=KB
          MAT=MATER(KB)
          RETURN
        ENDIF
 102  CONTINUE
      IBODY=NBODY+1
      MAT=0
      RETURN
      END
C  *********************************************************************
C                       SUBROUTINE STEPSI
C  *********************************************************************
      SUBROUTINE STEPSI(KB,S,IS,NSC)
C
C     This subroutine determines the intersections of the trajectory
C  with the limiting surfaces of body KB. The intersections are added to
C  the list and sorted in decreasing order. STEPSI should be invoked
C  only from within subroutine STEP.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N)
      PARAMETER (NS=100000,NB=40000,NX=1000)
      PARAMETER (EPS=1.0D-8,ONE=1.0D0+EPS)
      PARAMETER (NS2M=2*NS)
      DIMENSION S(NS2M),IS(NS2M)
      COMMON/QSURF/AXX(NS),AXY(NS),AXZ(NS),AYY(NS),AYZ(NS),AZZ(NS),
     1    AX(NS),AY(NS),AZ(NS),A0(NS),NSURF,NS2,KPLANE(NS)
      COMMON/QTREE/NBODY,MATER(NB),KMOTH(NB),KDGHT(NB,NX),
     1    KSURF(NB,NX),KFLAG(NB,NX),KSP(NS),KSLAST
C
      COMMON/TRACK/E,X,Y,Z,U,V,W,WGHT,KPAR,IBODY,MAT,ILB(5)
C
C  ************  Determine surface crossings.
C
      DO 100 KSS=1,KSURF(KB,NX)
C  ****  Notice that intersections with a given surface are calculated
C        only once. The side pointer of a surface must be changed each
C        time the surface is crossed.
        KFL=KFLAG(KB,KSS)
        IF(KFL.GT.4) GOTO 100
        KS=KSURF(KB,KSS)
        IF(KSP(KS).NE.0) GOTO 100
        IF(KPLANE(KS).EQ.0) THEN
          A=U*(AXX(KS)*U+AXY(KS)*V+AXZ(KS)*W)
     1     +V*(AYY(KS)*V+AYZ(KS)*W)+W*AZZ(KS)*W
          XXX=AXX(KS)*X+AXY(KS)*Y+AXZ(KS)*Z+AX(KS)
          YYY=AYY(KS)*Y+AYZ(KS)*Z+AY(KS)
          ZZZ=AZZ(KS)*Z+AZ(KS)
          B=U*(AXX(KS)*X+XXX)+V*(AXY(KS)*X+AYY(KS)*Y+YYY)
     1     +W*(AXZ(KS)*X+AYZ(KS)*Y+AZZ(KS)*Z+ZZZ)
          C=X*XXX+Y*YYY+Z*ZZZ+A0(KS)
        ELSE
          A=0.0D0
          B=U*AX(KS)+V*AY(KS)+W*AZ(KS)
          C=X*AX(KS)+Y*AY(KS)+Z*AZ(KS)+A0(KS)
        ENDIF
C
        IF(C.LT.0.0D0) THEN
          KSP(KS)=1
        ELSE
          KSP(KS)=2
        ENDIF
C
        IF(ABS(A).LT.1.0D-16) THEN
          IF(ABS(B).GT.1.0D-16*ABS(C)) THEN
            S1=-C/B
            IF(S1.GT.0.0D0) THEN
              NSC=NSC+1
              IS(NSC)=KS
              S(NSC)=S1
            ENDIF
          ENDIF
        ELSE
          DISCR=B*B-4.0D0*A*C
          IF(DISCR.GT.1.0D-35) THEN
            R2A=0.5D0/A
            D=SQRT(DISCR)
            S1=(-B-D)*R2A
            IF(S1.GT.0.0D0) THEN
              NSC=NSC+1
              IS(NSC)=KS
              S(NSC)=S1
            ENDIF
            S2=(-B+D)*R2A
            IF(S2.GT.0.0D0) THEN
              NSC=NSC+1
              IS(NSC)=KS
              S(NSC)=S2
            ENDIF
          ENDIF
        ENDIF
 100  CONTINUE
C
C  ****  Sort surface distances in decreasing order.
C
      IF(NSC.LT.2) RETURN
      DO KI=1,NSC-1
        SMAX=S(KI)
        KMAX=KI
        DO KJ=KI+1,NSC
          IF(S(KJ).GT.SMAX) THEN
            SMAX=S(KJ)
            KMAX=KJ
          ENDIF
        ENDDO
        IF(KMAX.NE.KI) THEN
          SMAX=S(KI)
          S(KI)=S(KMAX)
          S(KMAX)=SMAX
          KKMAX=IS(KI)
          IS(KI)=IS(KMAX)
          IS(KMAX)=KKMAX
        ENDIF
      ENDDO
C
      RETURN
      END
C  *********************************************************************
C                       SUBROUTINE STEPLB
C  *********************************************************************
      SUBROUTINE STEPLB(KB,S,IS,NSC,IERR)
C
C     This subroutine finds the body or module that has the prescribed
C  side pointers for the analysed surfaces. It should be invoked only
C  from within subroutine STEP.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N)
      PARAMETER (NS=100000,NB=40000,NX=1000)
      PARAMETER (EPS=1.0D-8,ONE=1.0D0+EPS)
      PARAMETER (NS2M=2*NS)
      DIMENSION S(NS2M),IS(NS2M)
      COMMON/QSURF/AXX(NS),AXY(NS),AXZ(NS),AYY(NS),AYZ(NS),AZZ(NS),
     1    AX(NS),AY(NS),AZ(NS),A0(NS),NSURF,NS2,KPLANE(NS)
      COMMON/QTREE/NBODY,MATER(NB),KMOTH(NB),KDGHT(NB,NX),
     1    KSURF(NB,NX),KFLAG(NB,NX),KSP(NS),KSLAST
C
      COMMON/TRACK/E,X,Y,Z,U,V,W,WGHT,KPAR,IBODY,MAT,ILB(5)
C
C  ****  Locate the particle at the sifted position..
C
      DO 100 KBB=1,KDGHT(KB,NX)
        KB1=KDGHT(KB,KBB)
        NST=KSURF(KB1,NX)
        DO KSS=1,NST
          KS=KSURF(KB1,KSS)
          KF=KFLAG(KB1,KSS)
          IF(KF.LT.3.AND.KSP(KS).NE.KF) GOTO 100
        ENDDO
        IBODY=KB1
        IF(KB1.EQ.KB) THEN
          IERR=0
          MAT=MATER(KB1)
          RETURN
        ELSE IF(KDGHT(KB1,NX).GT.1) THEN
          IERR=-1  ! The point is inside a submodule.
          RETURN
        ELSE
          IERR=0
          MAT=MATER(KB1)
          RETURN
        ENDIF
 100  CONTINUE
      IERR=1  ! The point is outside the module or body.
      IBODY=KMOTH(KB)
      IF(IBODY.EQ.0) THEN
        IBODY=NBODY+1
        MAT=0
      ENDIF
C
      RETURN
      END
C  *********************************************************************
C                       SUBROUTINE STEP
C  *********************************************************************
      SUBROUTINE STEP(DS,DSEF,NCROSS)
C
C     This subroutine handles the geometrical part of the track simula-
C  tion. The particle starts from the point (X,Y,Z) and travels a length
C  DS in the direction (U,V,W) within the material where it moves. When
C  the track leaves the initial material, the particle is stopped just
C  after entering the next material body (void regions with MAT=0 are
C  crossed automatically). Furthermore, when the particle arrives from
C  a void region, it is stopped just after entering the first material
C  body.
C
C  Input values (common /TRACK/):
C     X, Y, Z ... coordinates of the initial point,
C     U, V, W ... direction cosines of the displacement,
C     IBODY ..... body where the initial point is located,
C     MAT ....... material in body IBODY.
C
C  Input argument:
C     DS ........ path length to travel.
C
C  Output arguments:
C     DSEF....... travelled path length before leaving the initial
C                 material or completing the jump (less than DS if the
C                 track crosses an interface),
C     NCROSS .... = 0 if the whole step is contained in the initial
C                   material,
C                 .gt.0 if the particle has crossed an interface, i.e.
C                   if it has entered a new material.
C
C  Output values (common /TRACK/):
C     X, Y, Z ... coordinates of the final position,
C     IBODY ..... body where the final point is located,
C     MAT ....... material in IBODY. The value MAT=0 indicates that the
C                 particle has escaped from the system.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N)
      PARAMETER (NS=100000,NB=40000,NX=1000)
      PARAMETER (EPS=1.0D-8,EPSS=0.01D0*EPS,ONE=1.0D0+1.0D-4*EPS)
      PARAMETER (NS2M=2*NS)
      DIMENSION S(NS2M),IS(NS2M)
      COMMON/QSURF/AXX(NS),AXY(NS),AXZ(NS),AYY(NS),AYZ(NS),AZZ(NS),
     1    AX(NS),AY(NS),AZ(NS),A0(NS),NSURF,NS2,KPLANE(NS)
C
C  ****  When NCROSS.ne.0, the output value of the variable KSLAST is
C        the label of the last surface crossed by the particle before
C        entering a material body. KSLAST is used for rendering in 3D.
      COMMON/QTREE/NBODY,MATER(NB),KMOTH(NB),KDGHT(NB,NX),
     1    KSURF(NB,NX),KFLAG(NB,NX),KSP(NS),KSLAST
      COMMON/QKDET/KDET(NB)
C
      COMMON/TRACK/E,X,Y,Z,U,V,W,WGHT,KPAR,IBODY,MAT,ILB(5)
C
      DSEF=0.0D0
      NCROSS=0
      KSLAST=0
C
      NSC=0        ! Number of surface crossings ahead of the particle.
      DO I=1,NSURF
        KSP(I)=0   ! Side pointers of the evaluated surfaces.
      ENDDO
      MAT0=MAT     ! Initial material.
      IF(MAT.EQ.0) THEN
        DSRES=1.0D34
      ELSE
        DSRES=DS   ! Residual path length.
      ENDIF
C
C  ************  The particle enters from outside the enclosure.
C
      IF(IBODY.GT.NBODY) THEN
        KB1=NBODY
        CALL STEPSI(KB1,S,IS,NSC)
        IF(NSC.EQ.0) GOTO 400
        NSCT=NSC
        NST=KSURF(KB1,NX)
        DO KI=NSCT,1,-1
C  ***   The particle crosses a surface.
          KSLAST=IS(KI)
          IF(KSP(KSLAST).EQ.1) THEN
            KSP(KSLAST)=2
          ELSE
            KSP(KSLAST)=1
          ENDIF
          DSEF=DSEF+S(KI)
          DSP=ONE*S(KI)+EPSS
          X=X+DSP*U
          Y=Y+DSP*V
          Z=Z+DSP*W
          NSC=NSC-1
          IF(NSC.GT.0) THEN
            DO I=1,NSC
              S(I)=S(I)-DSP
            ENDDO
          ENDIF
C
          DO KSS=1,NST
            KS1=KSURF(KB1,KSS)
            KF=KFLAG(KB1,KSS)
            IF(KF.LT.3.AND.KSP(KS1).NE.KF) GOTO 101
          ENDDO
C  ****  The particle enters the enclosure.
 100      CONTINUE
          CALL STEPLB(KB1,S,IS,NSC,IERR)
          IF(IERR.EQ.-1) THEN
C  ****  The particle enters a submodule.
            KB1=IBODY
            CALL STEPSI(KB1,S,IS,NSC)
            GOTO 100
          ENDIF
C  ****  The particle enters a material body.
          IF(MAT.NE.0) THEN
            NCROSS=1
            RETURN
          ENDIF
          KB1=NBODY
          GOTO 300
 101      CONTINUE
        ENDDO
        GOTO 400
      ENDIF
C
C  ************  Surface crossings.
C
 200  CONTINUE
      SMIN=1.0D35  ! Distance to the nearest interface.
      KB0=IBODY
      DO 201 KSS=1,KSURF(KB0,NX)
        KFL=KFLAG(KB0,KSS)
        IF(KFL.GT.4) GOTO 201
        KS=KSURF(KB0,KSS)
        IF(KPLANE(KS).EQ.0) THEN
          A=U*(AXX(KS)*U+AXY(KS)*V+AXZ(KS)*W)
     1     +V*(AYY(KS)*V+AYZ(KS)*W)+W*AZZ(KS)*W
          XXX=AXX(KS)*X+AXY(KS)*Y+AXZ(KS)*Z+AX(KS)
          YYY=AYY(KS)*Y+AYZ(KS)*Z+AY(KS)
          ZZZ=AZZ(KS)*Z+AZ(KS)
          B=U*(AXX(KS)*X+XXX)+V*(AXY(KS)*X+AYY(KS)*Y+YYY)
     1     +W*(AXZ(KS)*X+AYZ(KS)*Y+AZZ(KS)*Z+ZZZ)
          C=X*XXX+Y*YYY+Z*ZZZ+A0(KS)
        ELSE
          A=0.0D0
          B=U*AX(KS)+V*AY(KS)+W*AZ(KS)
          C=X*AX(KS)+Y*AY(KS)+Z*AZ(KS)+A0(KS)
        ENDIF
C
        IF(C.LT.0.0D0) THEN
          KSP(KS)=1
        ELSE
          KSP(KS)=2
        ENDIF
C  ****  Consistency test, and position correction if not passed.
        IF(KFL.LT.3.AND.KFL.NE.KSP(KS)) THEN
          CALL LOCATE
          IF(MAT.EQ.MAT0) THEN
            DSEF=DSEF+EPS
            DO I=1,NSURF
              KSP(I)=0
            ENDDO
            NSC=0
            GOTO 200
          ENDIF
          NCROSS=1
          KSLAST=KS
          RETURN
        ENDIF
C
        IF(ABS(A).LT.1.0D-16) THEN
          IF(ABS(B).GT.1.0D-16*ABS(C)) THEN
            S1=-C/B
            IF(S1.GT.0.0D0) THEN
              IF(S1.LT.SMIN) THEN
                SMIN=S1
              ENDIF
              NSC=NSC+1
              IS(NSC)=KS  ! Surface label.
              S(NSC)=S1   ! Path length to the surface KS.
            ENDIF
          ENDIF
        ELSE
          DISCR=B*B-4.0D0*A*C
          IF(DISCR.GT.1.0D-35) THEN
            R2A=0.5D0/A
            D=SQRT(DISCR)
            S1=(-B-D)*R2A
            IF(S1.GT.0.0D0) THEN
              IF(S1.LT.SMIN) THEN
                SMIN=S1
              ENDIF
              NSC=NSC+1
              IS(NSC)=KS
              S(NSC)=S1
            ENDIF
            S2=(-B+D)*R2A
            IF(S2.GT.0.0D0) THEN
              IF(S2.LT.SMIN) THEN
                SMIN=S2
              ENDIF
              NSC=NSC+1
              IS(NSC)=KS
              S(NSC)=S2
            ENDIF
          ENDIF
        ENDIF
 201  CONTINUE
C
C  ****  The particle remains in the same material body.
C
      IF(MAT.NE.0.AND.DSRES.LT.SMIN) THEN
        IF(MAT.EQ.MAT0) DSEF=DSEF+DSRES
        X=X+DSRES*U
        Y=Y+DSRES*V
        Z=Z+DSRES*W
        RETURN
      ENDIF
C
C  ************  Surface crossings (sister bodies and modules).
C
      KB1=KMOTH(KB0)
      IF(KB1.EQ.0) KB1=KB0
      CALL STEPSI(KB1,S,IS,NSC)
C
C  ************  New position.
C
 300  CONTINUE
      NSCT=NSC
      DO KI=NSCT,1,-1
        MATL=MAT
C  ****  The step ends within the body.
        IF(S(KI).GT.DSRES) THEN
          IF(MAT.EQ.MAT0) THEN
            X=X+DSRES*U
            Y=Y+DSRES*V
            Z=Z+DSRES*W
            DSEF=DSEF+DSRES
            RETURN
          ENDIF
        ENDIF
C  ***   The particle crosses a surface.
        KSLAST=IS(KI)
        IF(KSP(KSLAST).EQ.1) THEN
          KSP(KSLAST)=2
        ELSE
          KSP(KSLAST)=1
        ENDIF
        DSP=S(KI)
C  ****  When the particle crosses a surface, it is shifted forward a
C        small additional length to make sure that it is placed behind
C        the surface.
        DSPT=ONE*DSP+EPSS
        X=X+DSPT*U
        Y=Y+DSPT*V
        Z=Z+DSPT*W
        IF(MAT.EQ.MAT0) THEN
          DSEF=DSEF+DSP
          DSRES=DSRES-DSP
        ENDIF
        NSC=NSC-1
        IF(NSC.GT.0) THEN
          DO I=1,NSC
            S(I)=S(I)-DSP
          ENDDO
        ENDIF
C
        IBODYL=IBODY
 301    CONTINUE
        CALL STEPLB(KB1,S,IS,NSC,IERR)
        IF(IERR.EQ.-1) THEN
C  ****  The particle enters a submodule.
          KB1=IBODY
          CALL STEPSI(KB1,S,IS,NSC)
          GOTO 301
        ELSE IF(IERR.EQ.1) THEN
C  ****  The particle leaves the body or module.
          IF(IBODY.LE.NBODY) THEN
            KB1=IBODY
            CALL STEPSI(KB1,S,IS,NSC)
            GOTO 301
          ELSE
C  ****  The particle leaves the enclosure.
            IF(MAT.NE.MATL) NCROSS=NCROSS+1
            DSP=1.0D30
            X=X+DSP*U
            Y=Y+DSP*V
            Z=Z+DSP*W
            RETURN
          ENDIF
        ENDIF
        IF(IBODY.EQ.IBODYL) GOTO 302
        KB1=IBODY
C  ****  The particle continues flying when it enters a new body of the
C        same material which is not a new detector...
        IF(MAT.EQ.MATL) THEN
          IF(KDET(IBODY).EQ.0) GOTO 300
          IF(KDET(IBODY).EQ.KDET(IBODYL)) GOTO 300
        ENDIF
C  ****  The particle continues flying when it enters a void region...
        IF(MAT.EQ.0) THEN
          IF(MATL.EQ.MAT0) THEN
            NCROSS=NCROSS+1
          ENDIF
          DSRES=1.0D30
          GOTO 300
C  ****  ... and stops when it penetrates a new material body or a
C        detector.
        ELSE IF(MAT.NE.MATL.OR.KDET(IBODY).NE.0) THEN
          NCROSS=NCROSS+1
          X=X+EPS*U
          Y=Y+EPS*V
          Z=Z+EPS*W
          RETURN
        ENDIF
 302    CONTINUE
      ENDDO
C  ****  The particle leaves the enclosure.
 400  CONTINUE
      DSP=1.0D30
      IBODY=NBODY+1
      X=X+DSP*U
      Y=Y+DSP*V
      Z=Z+DSP*W
C
      RETURN
      END
