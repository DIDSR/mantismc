!*******************************************************************
!*                           PENAUX                                *
!*                                                                 *
!* Short description:                                              *
!*   Miscelaneous routines that serve to initialize and monitor a  *
!*   PENELOPE main program.                                        *
!*                                                                 *
!* Dependencies:                                                   *
!*   from PENELOPE:                                                *
!*   -> common /TRACK/                                             *
!*   -> common /RSEED/                                             *
!*   -> common /CSIMPA/ (includes MAXMAT, max number of materials) *
!*   -> routines PEINIT, PANR, PHMFP                               *
!*   from PENGEOM:                                                 *
!*   -> routine GEOMIN                                             *
!*   -> common /QSURF/                                             *
!*   -> common /QTREE/                                             *
!*   from PENVARED                                                 *
!*   -> common /CFORCE/                                            *
!*   from LIBTIM:                                                  *
!*   -> routines CPUTIME and REALTIME                              *
!*   from other penEasy files:                                     *
!*   -> routines TALLY and REPORT                                  *
!*                                                                 *
!* Compatible with PENELOPE versions:                              *
!*   2005                                                          *
!*                                                                 *
!* Last update:                                                    *
!*   2005-08-01                                                    *
!*     Adapted to PENELOPE'2005 and MANTIS                         *
!*******************************************************************


      subroutine iniconfig(funit)
!*******************************************************************
!*    Reads simulation configuration.                              *
!*    To be called during initialization.                          *
!*                                                                 *
!*    Input:                                                       *
!*      funit -> input unit (must be already open)                 *
!*    Output:                                                      *
!*      /ctrsim/ /rseed/                                           *
!*******************************************************************
      implicit none
      integer funit

      integer*4 seed1,seed2
      common/rseed/seed1,seed2
      integer*4 nhist,atime,refres,ncalls,lastim
      real*8 time0,accura
      common /ctrsim/ time0,accura,atime,refres,ncalls,lastim,nhist
      character*80 buffer
      character*(*) secid,eos
      parameter (secid='[SECTION CONFIG v.2003-04-21]')
      parameter (eos='[END OF SECTION]')
      integer in,finduf,idpresent

      write(*,*) ' '
      write(*,'(a)')
     & '>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'
      read(funit,'(a80)') buffer
      idpresent = index(buffer,secid)
      if (idpresent.eq.0) then
        write(*,'(a)') 'iniconfig:ERROR: incorrect section;'
        write(*,'(a,a)') '  expecting to find: ',secid
        write(*,'(a,a)') '  found instead:     ',buffer
        stop
      endif
      write(*,'(a)') secid

      write(*,'(a)') 'No. of histories requested: '
      read(funit,*) nhist
      write(*,'(i11)') nhist
      write(*,'(a)') 'Allotted time (s):'
      read(funit,*) atime
      write(*,'(i11)') atime
      write(*,'(a)') 'Accuracy limit (%):'
      read(funit,*) accura
      write(*,'(1pe12.5)') accura
      write(*,'(a)') 'Refresh frequency:'
      read(funit,*) refres
      write(*,'(i11)') refres
      if (refres.gt.50000) then
        write(*,*) 'iniconfig:ERROR: refresh time must be < 50000'
        stop
      endif
      write(*,'(a)') 'Random seeds:'
      read(funit,*) seed1,seed2
      if (seed1.eq.0.and.seed2.eq.0) then
        ! Reads seeds from an external file:
        write(*,'(a)') '  (reading seeds for external file):'
        read(funit,'(1x,a30)') buffer
        write(*,'(2x,a)') buffer
        write(*,*) ' '
        in = finduf()
        open(in,file=buffer,status='old')
        read(in,*) seed1,seed2
        close(in)
      endif
      write(*,'(2(1x,i11))') seed1,seed2
      ncalls = 0

      read(funit,'(a80)') buffer
      if (index(buffer,eos).eq.0) then
        write(*,*) 'iniconfig:ERROR: End-Of-Section mark not found'
        write(*,'(a,a)') '  expecting to find: ',eos
        write(*,'(a,a)') '  found instead:     ',buffer
        stop
      endif
      write(*,'(a)') '>>>> Config initialization finished >>>>'
      end


      subroutine inipen(funit,e0)
!*******************************************************************
!*    Reads simulation parameters and initializes PENELOPE and     *
!*    PENGEOM. To be called during initialization.                 *
!*                                                                 *
!*    Input:                                                       *
!*      funit -> input unit (must be already open)                 *
!*      e0 -> source energy (eV)                                   *
!*******************************************************************
      implicit none
      integer funit
      real*8 e0

      integer*4 maxmat
      parameter (maxmat=10)
      real*8 eabs,c1,c2,wcc,wcr
      common/csimpa/eabs(3,maxmat),c1(maxmat),c2(maxmat),
     &              wcc(maxmat),wcr(maxmat)
      real*8 dsmaxval
      common /ctrsi1/ dsmaxval(maxmat)
      character*80 buffer,mfilen
      character*(*) secid,eos
      parameter (secid=
     &'[SECTION PENELOPE+PENGEOM v.2003-04-21]')
      parameter (eos='[END OF SECTION]')
      integer in,finduf,imat,i,idpresent
      integer*4 nmat,nmatg,nbody
      real*8 parinp,infty,emax
      parameter (infty=1.0d30)

      write(*,*) ' '
      write(*,'(a)')
     & '>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'
      read(funit,'(a80)') buffer
      idpresent = index(buffer,secid)
      if (idpresent.eq.0) then
        write(*,'(a)') 'inipen:ERROR: incorrect section;'
        write(*,'(a,a)') '  expecting to find: ',secid
        write(*,'(a,a)') '  found instead:     ',buffer
        stop
      endif
      write(*,'(a)') secid

!     *** Find a valid file unit:
      in = finduf()

!     **** Init geometry:
      write(*,'(a)') 'Opening geometry data file:'
      read(funit,'(1x,a30)') buffer
      write(*,'(a)') buffer
      open(in,file=buffer)
      write(*,'(a)') 'Now calling GEOMIN...'
      write(*,'(a)') ' '
      call geomin(parinp,0,nmatg,nbody,in,6)
      close(in)
      write(*,'(a)') ' '
      write(*,'(a)') '>>>> Geometry initialization finished >>>>'

      write(*,'(a)') ' '
      write(*,'(a)')
     & '>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'
      read(funit,'(1x,a30)') mfilen  ! Material file name
!     *** Read simulation parameters:
      write(*,'(a)') ' '
      write(*,'(a)') 'SIMULATION PARAMETERS'
      write(*,'(a)') 'No. materials expected in MAT file (0=auto):'
      read(funit,*) nmat
      write(*,'(i3)') nmat
      if (nmat.gt.maxmat) then
        write(*,*) 'inipen:ERROR: Too many materials; enlarge MAXMAT'
        stop
      endif
      write(*,'(a86)') 'Mat EABS(e-)   EABS(ph)   EABS(e+)   C1'//
     &           '         C2        WCC         WCR        DSMAX'
      if (nmat.eq.0) then
        write(*,'(a)') '(default values requested)'
        nmat=nmatg
        do i=1,nmat
          eabs(1,i) = min(max(e0*1.0d-2,1.0d2),1.0d6)
          eabs(2,i) = min(max(e0*1.0d-3,1.0d2),1.0d6)
          eabs(3,i) = eabs(1,i)
          c1(i)     = 0.1d0
          c2(i)     = 0.1d0
          wcc(i)    = min(e0*1.0d-2,eabs(1,i))
          wcr(i)    = min(e0*1.0d-3,eabs(2,i))
          dsmaxval(i)  = infty
          write(*,'(i2,8(1x,1pe10.3))') i,eabs(1,i),eabs(2,i),
     &      eabs(3,i),c1(i),c2(i),wcc(i),wcr(i),dsmaxval(i)
        enddo
      else
        ! User has given simulation parameters explicitly:
        if (nmatg.gt.nmat) then
          write(*,*)
     &      'inipen:ERROR: There are more materials in the geometry'
          write(*,*)
     &      '   file than expected in the materials file'
          stop
        endif
        read(funit,'(a80)') buffer
        do i=1,nmat
          read(funit,*) imat,eabs(1,i),eabs(2,i),eabs(3,i),
     &              c1(i),c2(i),wcc(i),wcr(i),dsmaxval(i)
          write(*,'(i2,8(1x,1pe10.3))') imat,eabs(1,i),eabs(2,i),
     &              eabs(3,i),c1(i),c2(i),wcc(i),wcr(i),dsmaxval(i)
          if (imat.ne.i) then
            write(*,*)
     &      'inipen:ERROR: Materials must be ordered sequentially:'
            write(*,'(i3)') imat
            stop
          endif
        enddo
      endif

      ! Init PENELOPE kernel:
      emax = 1.0d9
      write(*,'(a)') 'Opening material data file:'
      write(*,'(a)') mfilen
      open(in,file=mfilen)
      write(*,'(a)') 'Now calling PEINIT...'
      write(*,'(a)') ' '
      call peinit(emax,nmat,in,6,1)
      close(in)

      read(funit,'(a80)') buffer
      if (index(buffer,eos).eq.0) then
        write(*,*) 'inipen:ERROR: End-Of-Section mark not found'
        write(*,'(a,a)') '  expecting to find: ',eos
        write(*,'(a,a)') '  found instead:     ',buffer
        stop
      endif
      write(*,'(a)') ' '
      write(*,'(a)') '>>>> PENELOPE initialization finished >>>>'
      end


      logical function endsim(n)
!*******************************************************************
!*    Determines whether the simulation is done or not and writes  *
!*    progress reports.                                            *
!*                                                                 *
!*    Input:                                                       *
!*      n -> current history no.                                   *
!*    Output:                                                      *
!*      .true. if done                                             *
!*******************************************************************
      implicit none
      integer*4 n

      integer*4 nhist,atime,refres,ncalls,lastim
      real*8 time0,accura
      common /ctrsim/ time0,accura,atime,refres,ncalls,lastim,nhist
      logical needfresh
      integer*4 now
      real*8 realtime,cputime,uncert

      endsim = .false.

!     *** Check no. of histories:
      if (n.ge.nhist) endsim = .true.

!     *** Check time:
      if (atime.gt.0) then
        if (realtime().gt.atime) endsim = .true.  ! Real time
      else
        if (cputime().gt.-atime) endsim = .true.  ! User (CPU) time
      endif

!     *** Progress report:
      if (endsim) then
        call report(n,cputime(),0,uncert)
        return
      endif

!     *** Check if a refresh is needed:
      needfresh = .false.
      now = realtime() ! Ensures one call per day at least
      if (refres.gt.0) then
        if (now-lastim.gt.refres) then
          lastim = now
          needfresh = .true.
        endif
      else
        ncalls = ncalls+1
        if (ncalls.ge.-refres) then
          ncalls = 0
          needfresh = .true.
        endif
      endif

!     *** Refresh:
      if (needfresh) then
        call comand(n)
        call report(n,cputime(),0,uncert)    ! Writes progress report
        if (uncert.lt.accura) endsim = .true.  ! Required accuracy attained
      endif
      end


      subroutine dostop
!*******************************************************************
!*    Forces the simulation to stop after completion of current    *
!*    history by re-setting the no. of histories to simulate to    *
!*    zero.                                                        *
!*******************************************************************
      implicit none
      integer*4 nhist,atime,refres,ncalls,lastim
      real*8 time0,accura
      common /ctrsim/ time0,accura,atime,refres,ncalls,lastim,nhist

      nhist = 0
      end


      real*8 function dsmax()
!*******************************************************************
!*    Maximum step length as passed to JUMP; this function is an   *
!*    interface to common /CTRSI1/                                 *
!*                                                                 *
!*    Input:                                                       *
!*      /TRACK/                                                    *
!*******************************************************************
      implicit none
      integer*4 kpar,ibody,mat,ilb
      real*8 e,x,y,z,u,v,w,wght
      common/track/e,x,y,z,u,v,w,wght,kpar,ibody,mat,ilb(5)
      integer*4 maxmat
      parameter (maxmat=10)
      real*8 dsmaxval
      common /ctrsi1/ dsmaxval(maxmat)

      dsmax = dsmaxval(mat)
      end


      logical function absorb()
!*******************************************************************
!*    Checks whether or not a particle is absorbed given its       *
!*    current energy.                                              *
!*                                                                 *
!*    Input:                                                       *
!*      /TRACK/                                                    *
!*    Output:                                                      *
!*      .true. if absorbed, .false. else                           *
!*    Comments:                                                    *
!*      -> Tallies deposited energy when necessary; to do so it    *
!*         calls TALLY.                                            *
!*      -> Absorbs particles above 1 GeV; these can be created by  *
!*         SOURCE or be the result of the annihilation of a e+.    *
!*      -> This function is well suited to implement range         *
!*         rejection by making the returned value dependent        *
!*         on the particle's position, material, energy, etc.      *
!*         Notice that in order to preserve the calculation of     *
!*         the fluence unbiased an electron or positron should     *
!*         _never_ be rejected _inside_the_detector_               *
!*******************************************************************
      implicit none
      integer*4 kpar,ibody,mat,ilb
      real*8 e,x,y,z,u,v,w,wght
      common/track/e,x,y,z,u,v,w,wght,kpar,ibody,mat,ilb(5)
      integer*4 maxmat
      parameter (maxmat=10)
      real*8 eabs,c1,c2,wcc,wcr
      common/csimpa/eabs(3,maxmat),c1(maxmat),c2(maxmat),
     &              wcc(maxmat),wcr(maxmat)
      real*8 mc2,twomc2
      parameter (mc2=5.10998902d5,twomc2=2.0d0*mc2)  ! Electron rest energy

      if (mat.eq.0) then  ! vacuum
        absorb = .true.
        return
      endif
      if (e.gt.eabs(kpar,mat).and.e.lt.1.0d9) then
        absorb = .false.
        return
      endif

      absorb = .true.
      call tally(-98,e)         ! Tallies remaining kinetic energy
      if (kpar.eq.3.and.e.gt.0.0d0) then  ! Precaution in case rejected positron
        call panr               ! Puts two annihilation photons in the stack
        call tally(-97,twomc2)  ! Tally the extra 2mc^2
      endif
      if (e.gt.1.0d9) then
        write(*,'(a)') ' '
        write(*,'(a)') '***************'
        write(*,'(a,1pe12.5,a)')
     &  'absorb:WARNING: particle with energy out of range: ',e,' eV'
        write(*,'(a)') '  (it has been locally absorbed)'
        write(*,'(a)') '***************'
        write(*,'(a)') ' '
      endif
      end


      subroutine comand(n)
!*******************************************************************
!*    Reads and executes commands from an external file, allowing  *
!*    in-flight steering of the simulation.                        *
!*                                                                 *
!*    Input:                                                       *
!*      n -> current no. of histories already simulated.           *
!*      -> From file command.in.  Possible codes are:              *
!*          0     -> Do nothing                                    *
!*          1     -> Obsolete, unused                              *
!*          2 <N> -> Reset No of histories to <N>                  *
!*          3 <t> -> Reset simulation time to <t>, in s            *
!*          4 <A> -> Reset accuracy limit                          *
!*          5 <R> -> Reset refresh frequency to <R>                *
!*          6 <string> -> write <string> to stdout                 *
!*    Output:                                                      *
!*      -> Command is executed and file command.in reset to default*
!*    Comments:                                                    *
!*      -> All command arguments (such as <N>) must be preceded    *
!*         by a line feed (as in '2 <return> 200000').             *
!*      -> <string> must not exceed 80 chars.                      *
!*      -> Note that the simulation may be stopped immediately by  *
!*         resetting <N> to 0.                                     *
!*      -> Be careful before saving file command.in ;              *
!*         a syntax error may cause the program to abort.          *
!*    Dependencies:                                                *
!*      None.                                                      *
!*******************************************************************
      implicit none
      integer*4 n

      integer*4 nhist,atime,refres,ncalls,lastim
      real*8 time0,accura
      common /ctrsim/ time0,accura,atime,refres,ncalls,lastim,nhist
      character*80 buffer
      logical ffound
      integer com,in,finduf,error
      integer*4 xatime,xnp,xrefre
      real*8 xaccur

      inquire(file='command.in',exist=ffound)
      if (.not.ffound) return
      com = 0
      in = finduf()

!     *** Read command file:
      open(in,file='command.in',iostat=error)
      read(in,*,iostat=error) com
      if (error.ne.0.or.com.eq.0) goto 20
      write(*,*) ' '
      write(*,*) 'comand: Command received when Nhist was:'
      write(*,'(2x,i10)') n
      write(*,*) '        Command description:'
      if (com.eq.1) then
        write(*,*) '  code 1, obsolete; it has been ignored.'
      else if (com.eq.2) then
        read(in,*,iostat=error) xnp
        if (error.ne.0) goto 20
        nhist = xnp
        write(*,*) ' max Nhist reset to:'
        write(*,'(2x,i10)') nhist
      else if (com.eq.3) then
        read(in,*,iostat=error) xatime
        if (error.ne.0) goto 20
        atime = xatime
        write(*,*) ' max time reset to:'
        write(*,'(2x,i11)') atime
      else if (com.eq.4) then
        read(in,*,iostat=error) xaccur
        if (error.ne.0) goto 20
        accura = xaccur
        write(*,*) ' accuracy limit reset to:'
        write(*,'(2x,1pe12.5)') accura
      else if (com.eq.5) then
        read(in,*,iostat=error) xrefre
        if (error.ne.0) goto 20
        refres = xrefre
        write(*,*) ' refresh frequency reset to:'
        write(*,'(2x,i10)') refres
      else if (com.eq.6) then
        read(in,'(a80)',iostat=error) buffer
        if (error.ne.0) goto 20
        write(*,*) '*** Message from command.in  follows:'
        write(*,'(1x,a80)') buffer
      endif
 20   continue
      close(in,status='delete',iostat=error)

!     *** Reset command file:
      open(in,file='command.in',iostat=error)
      if (error.ne.0) return
      write(in,*,iostat=error) ' 0'
      if (error.ne.0) return    ! Write protected file
      write(in,*) '  '
      write(in,*) '*** Codes:  '
      write(in,*) '0     -> Do nothing'
      write(in,*) '1     -> Obsolete, unused'
      write(in,*) '2 <N> -> Reset No of histories'
      write(in,*) '3 <t> -> Reset simulation time'
      write(in,*) '4 <A> -> Reset accuracy limit'
      write(in,*) '5 <R> -> Reset refresh frequency'
      write(in,*)
     & '6 <string> -> Write <string> (up to 80 chars) to stdout'
      write(in,*) ' '
      close(in,iostat=error)
      end


      integer function finduf()
!*******************************************************************
!*    Finds a valid (i.e. unused) file unit.                       *
!*******************************************************************
      implicit none
      logical used
      integer maxuf
      parameter (maxuf=17)

      finduf = 6
      do
        finduf = finduf+1
        if (finduf.gt.maxuf) then
          write(*,*) ' '
          write(*,*) 'finduf:ERROR: Unable to find a valid unit'
          stop
        endif
        inquire(finduf,opened=used)
        if (.not.used) return
      enddo
      end


      subroutine iniforce(funit,e0)
!*******************************************************************
!*    Initializes the interaction forcing routines                 *
!*                                                                 *
!*    Input:                                                       *
!*      funit -> input unit (must be already open)                 *
!*      e0 -> source energy (eV)                                   *
!*******************************************************************
      implicit none
      integer funit

      real*8 e0
      integer ns,nb,nx
      PARAMETER (NS=100000,NB=40000,NX=1000)
      integer*4 nbody,mater,kmoth,kdght,ksurf,kflag,kalias,kslast
      common/qtree/nbody,mater(nb),kmoth(nb),kdght(nb,nx),
     1    ksurf(nb,nx),kflag(nb,nx),kalias(ns),kslast
      real*8 force
      common/cforce/force(nb,3,8)
      integer*4 maxmat
      parameter (maxmat=10)
      logical analog
      real*8 minwght,maxwght
      common /frc001/ analog(maxmat,3),minwght,maxwght
      logical active
      character*80 buffer
      character*(*) secid,eos
      parameter (secid=
     &'[SECTION INTERACTION FORCING v.2005-06-02]')
      parameter (eos='[END OF SECTION]')
      integer i,j,idpresent,mat,kpar,icol,error
      real*8 fforce,phmfp,hmfp

      write(*,*) ' '
      write(*,'(a)')
     & '>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'
      read(funit,'(a80)') buffer
      idpresent = index(buffer,secid)
      if (idpresent.eq.0) then
        write(*,'(a)') 'iniconfig:ERROR: incorrect section;'
        write(*,'(a,a)') '  expecting to find: ',secid
        write(*,'(a,a)') '  found instead:     ',buffer
        stop
      endif
      write(*,'(a)') secid

      ! Reset all:
      do icol=1,8
        do kpar=1,3
          do i=1,nb
            force(i,kpar,icol) = 1.0d0
          enddo
        enddo
      enddo
      do kpar=1,3
        do mat=1,maxmat
          analog(mat,kpar) = .true.
        enddo
      enddo

      write(*,'(a)') 'MAT -- KPAR -- ICOL -- Forcing -- Wmin -- Wmax:'
      active = .false.
      do
        read(funit,*,iostat=error) mat,kpar,icol,fforce,minwght,maxwght
        if (error.ne.0) then
          write(*,'(a)')
     &     'iniforce:ERROR: unable to read line containing:'
          write(*,'(a)')
     &     '  MAT,KPAR,ICOL,FORCING,MINWGHT,MAXWGHT'
          write(*,'(a)')
     &     '  Recall to use a line with MAT=0 to end the list.'
          stop
        endif
        if (mat.eq.0) exit  ! End of list
        write(*,'(3(1x,i5),3(1x,1pe12.5))')
     &    mat,kpar,icol,fforce,minwght,maxwght
        active = .true.
        if (mat.lt.0) then
          write(*,*) 'iniforce:ERROR: MAT must be positive'
          stop
        endif
        if (kpar.lt.1.or.kpar.gt.3) then
          write(*,*) 'iniforce:ERROR: KPAR must be in [1,3]'
          stop
        endif
        if (icol.ne.0.and.(icol.lt.2.or.icol.gt.6)) then
          write(*,*) 'iniforce:ERROR: ICOL must be 0 or in [2,6]'
          stop
        endif
        if (fforce.le.1.0d0) then
          write(*,*) 'iniforce:ERROR: FORCING must be in > 1'
          stop
        endif
        do i=1,nb  ! Sets forcing for all bodies of the selected material
          if (mater(i).eq.mat) then
            if (icol.ne.0) then
              force(i,kpar,icol) = fforce
            else  ! Force all interactions
              do j=1,8
                force(i,kpar,j) = fforce
              enddo
            endif
          endif
        enddo
        analog(mat,kpar) = .false.
      enddo

      if (active) then
        write(*,'(a)') ' '
        write(*,'(a,1pe12.5)')
     &    'INFO: Hard mean free paths at E(eV) =',e0
        write(*,*) 'MAT -- KPAR -- ICOL -- HMFP(cm)'
        do mat=1,maxmat
          do kpar=1,3
            if (analog(mat,kpar)) cycle
            do icol=1,8
              if (kpar.eq.1.and.(icol.lt.2.or.icol.gt.5)) cycle
              if (kpar.eq.2.and.icol.gt.4) cycle
              if (kpar.eq.3.and.(icol.lt.2.or.icol.gt.6)) cycle
              hmfp = phmfp(e0,kpar,mat,icol)
              write(*,'(i3,1x,i1,1x,i1,1x,1pe12.5)')
     &          mat,kpar,icol,hmfp
            enddo
          enddo
        enddo
        write(*,'(a)') ' '
        write(*,'(a)') '**** iniforce: WARNING ****'
        write(*,'(a)')
     &    '  Interaction forcing may bias pulse height spectra.'
      endif

      read(funit,'(a80)') buffer
      if (index(buffer,eos).eq.0) then
        write(*,*) 'iniforce:ERROR: End-Of-Section mark not found'
        write(*,'(a,a)') '  expecting to find: ',eos
        write(*,'(a,a)') '  found instead:     ',buffer
        stop
      endif

      if (active) then
        write(*,'(a)') '>>>> Forcing initialization finished >>>>'
      else
        write(*,'(a)') '>>>> Interaction forcing is INACTIVE >>>>'
      endif
      end


      subroutine isforcing
!*******************************************************************
!*    Determines whether interaction forcing is to be applied      *
!*                                                                 *
!*    Output:                                                      *
!*      -> logical forcing (in a common block)                     *
!*******************************************************************
      integer*4 kpar,ibody,mat,ilb
      real*8 e,x,y,z,u,v,w,wght
      common/track/e,x,y,z,u,v,w,wght,kpar,ibody,mat,ilb(5)
      integer*4 maxmat
      parameter (maxmat=10)
      logical analog
      real*8 minwght,maxwght
      common /frc001/ analog(maxmat,3),minwght,maxwght
      logical forcing
      common /frc002/ forcing

      forcing = .false.
      if (mat.eq.0) return
      if (analog(mat,kpar).or.wght.lt.minwght.or.wght.gt.maxwght)
     &  return
      forcing = .true.
      end


      subroutine snormal(gradx,grady,gradz)
!*******************************************************************
!*    When called immediately after a surface intersection (with   *
!*    STEP) it provides the normal to that surface                 *
!*                                                                 *
!*    Output:                                                      *
!*      grad{xyz} -> gradient to the surface (normalized)          *
!*******************************************************************
      implicit none
      real*8 gradx,grady,gradz

      integer*4 kpar,ibody,mat,ilb
      real*8 e,x,y,z,u,v,w,wght
      common/track/e,x,y,z,u,v,w,wght,kpar,ibody,mat,ilb(5)
      integer ns,nb,nx
      PARAMETER (NS=100000,NB=40000,NX=1000)
      integer*4 nsurf,ns2,kplane
      real*8 axx,axy,axz,ayy,ayz,azz,ax,ay,az,a0
      common/qsurf/axx(ns),axy(ns),axz(ns),ayy(ns),ayz(ns),azz(ns),
     1    ax(ns),ay(ns),az(ns),a0(ns),nsurf,ns2,kplane(ns)
      integer nbody,mater,kmoth,kdght,ksurf,kflag,kalias,kslast
      common/qtree/nbody,mater(nb),kmoth(nb),kdght(nb,nx),
     1    ksurf(nb,nx),kflag(nb,nx),kalias(ns),kslast
      real*8 norm

      ! Quadric surface gradient:
      gradx = 2.0d0*axx(kslast)*x+axy(kslast)*y+axz(kslast)*z+
     &        ax(kslast)
      grady = axy(kslast)*x+2.0d0*ayy(kslast)*y+ayz(kslast)*z+
     &        ay(kslast)
      gradz = axz(kslast)*x+ayz(kslast)*y+2.0d0*azz(kslast)*z+
     &        az(kslast)
      norm = sqrt(gradx*gradx+grady*grady+gradz*gradz)
      if (norm.gt.0.0d0) then
        norm = 1.0d0/norm
        gradx = gradx*norm
        grady = grady*norm
        gradz = gradz*norm
      endif
      end


      integer function strlen(string)
c*******************************************************************
c*    Returns length of string excluding trailing blanks           *
c*                                                                 *
c*    Input:                                                       *
c*      string -> character string                                 *
c*******************************************************************
      implicit none
      character*(*) string

      do strlen=len(string),1,-1
        if (string(strlen:strlen).ne.' ') return
      enddo
      strlen = 0
      end


! **** end of file *************************************************

