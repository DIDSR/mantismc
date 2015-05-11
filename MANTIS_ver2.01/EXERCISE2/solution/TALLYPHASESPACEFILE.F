!*******************************************************************
!*                          TALLY                                  *
!*                    PHASE SPACE FILE                             *
!*                                                                 *
!* Short description:                                              *
!*   Tally routines for radiation transport calculations with      *
!*   PENELOPE.                                                     *
!*                                                                 *
!*   Writes to a file the state of all particles that reach        *
!*   a given material, considered as a particle sink.              *
!*                                                                 *
!* Dependencies:                                                   *
!*   from PENELOPE:                                                *
!*   -> common /TRACK/                                             *
!*   -> common /RSEED/                                             *
!*   from other penEasy files:                                     *
!*   -> routine FINDUF                                             *
!*   -> routine STRLEN                                             *
!*                                                                 *
!* Compatible with PENELOPE versions:                              *
!*   2005                                                          *
!*                                                                 *
!* Last update:                                                    *
!*   2005-08-01                                                    *
!*     Adapted to MANTIS                                           *
!*******************************************************************


      subroutine PSFtally(mode,nhnow)
!*******************************************************************
!*    Input:                                                       *
!*      mode -> Identifies the state of the calling procedure      *
!*      nhnow -> current history number (only when mode=1)         *
!*******************************************************************
      implicit none
      integer mode
      real*8 nhnow

      integer*4 kpar,ibody,mat,ilb
      real*8 e,x,y,z,u,v,w,wght
      common/track/e,x,y,z,u,v,w,wght,kpar,ibody,mat,ilb(5)
      logical active
      integer detmat,psfunit
      integer*4 nhist,nhlast,npar
      common /scopsf/ nhist,nhlast,npar,detmat,psfunit,active
      character*80 fmtstr
      integer*4 dn

      if (.not.active) return

      if (mode.eq.4.or.mode.eq.-99) then  ! New material
        if (mat.ne.detmat) return
        dn = nhist-nhlast
        fmtstr = '(i1,8(1x,1pe12.5),1x,i2,1x,i10)'
        write(psfunit,fmtstr) kpar,e,x,y,z,u,v,w,wght,ilb(5),dn
        nhlast = nhist
        npar = npar+1
        mat = 0            ! Forces main to stop this particle's simulation

      else if (mode.eq.1) then  ! Generation of a new history
        nhist = nhnow+0.1d0     ! Safe real to integer conversion

      endif
      end


      subroutine PSFreport(n,cputim,screen,uncert)
!*******************************************************************
!*    Input:                                                       *
!*      n -> no. of histories simulated                            *
!*      cputim -> elapsed CPU time                                 *
!*      screen -> 1 if output goes to screen, 0 else               *
!*    Output:                                                      *
!*      uncert -> overall uncertainty in % (at 2.sigma)            *
!*    Comments:                                                    *
!*      -> 'cputim' should not include initialization procedures;  *
!*         enter 0 or neg. if not available.                       *
!*******************************************************************
      implicit none
      integer screen
      integer*4 n
      real*8 cputim,uncert

      integer*4 seed1,seed2
      common/rseed/seed1,seed2
      logical active
      integer detmat,psfunit
      integer*4 nhist,nhlast,npar
      common /scopsf/ nhist,nhlast,npar,detmat,psfunit,active
      character*100 jobname
      common /cjob/ jobname
      character*100 outfilen
      integer out,finduf,error,strlen

      uncert = 0.0d0    ! Needed in case of inactivity
      if (.not.active) return
      uncert = 100.0d0  ! Never stop because of uncertainty

!     *** Prepare output files:
      out = 6    ! Default is to write to screen
      if (screen.eq.0) then
        out = finduf()
        outfilen =
     &    jobname(:strlen(jobname))//'tallyPhaseSpaceFile.dat'
        open(out,file=outfilen,iostat=error)
        if (error.eq.0) write(out,*,iostat=error) ' '
        if (error.ne.0) then
          write(*,*) ' '
          write(*,*) '********************************************'
          write(*,*)
     &      'PSFreport:ERROR: cannot write output data file;'
          write(*,*) '  writing to screen instead'
          write(*,*) '********************************************'
          write(*,*) ' '
          out = 6
        endif
      endif

      write(out,*) ' '
      write(out,'(a)')
     &'#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'
      write(out,'(a)') '# [SECTION REPORT PSF]'
      write(out,'(a)') '# No. of particles written to PSF:'
      write(out,'(i10)') npar

!     *** Generic report:
      write(out,'(a)') ' '
      write(out,'(a)') '# Performance report'
      write(out,'(a)') '#   Random seeds:'
      write(out,'(a,i10)') '#   ',seed1
      write(out,'(a,i10)') '#   ',seed2
      write(out,'(a)') '#   No. of histories simulated [N]:'
      write(out,'(a,i10)') '#   ',n
      write(out,'(a)') '#   CPU time [t] (s):'
      write(out,'(a,1pe12.5)') '#   ',cputim
      if (cputim.gt.0.0d0) then
        write(out,'(a)') '#   Speed (histories/s):'
        write(out,'(a,1pe12.5)') '#   ',n/cputim
      endif
      write(out,'(a)') '#'
      write(out,'(a)') '# Have a nice day.'
      write(out,'(a)') ' '

      if (out.ne.6) close(out)
      end


      subroutine PSFinitally(funit)
!*******************************************************************
!*    Initializes. To be called before TALLY                       *
!*                                                                 *
!*    Input:                                                       *
!*      funit -> input unit (must be already open)                 *
!*******************************************************************
      implicit none
      integer funit

      logical active
      integer detmat,psfunit
      integer*4 nhist,nhlast,npar
      common /scopsf/ nhist,nhlast,npar,detmat,psfunit,active
      character*80 buffer
      character*(*) secid,eos
      parameter (secid=
     &'[SECTION TALLY PHASE SPACE FILE v.2004-04-06]')
      parameter (eos='[END OF SECTION]')
      integer finduf,idpresent

      active = .true.
      write(*,*) ' '
      write(*,'(a)')
     & '>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'
      read(funit,'(a80)') buffer
      idpresent = index(buffer,secid)
      if (buffer(1:1).eq.'#'.and.idpresent.ne.0) then
        active = .false.
        write(*, '(a)')
     &    '>>>> Tally Phase Space File is INACTIVE >>>>'
        do
          read(funit,'(a80)',end=10) buffer  ! Reads dummy lines of inactive tally
          if (index(buffer,eos).ne.0) return  ! End-Of Section found
        enddo
 10     return  ! End-Of-File, something went wrong
      else if (idpresent.eq.0) then
        write(*,'(a)') 'PSFinitally:ERROR: incorrect section;'
        write(*,'(a,a)') '  expecting to find: ',secid
        write(*,'(a,a)') '  found instead:     ',buffer
        stop
      endif
      write(*,'(a)') secid

      write(*,'(a)') 'Detection material set to:'
      read(funit,*) detmat
      write(*,'(i3)') detmat
      if (detmat.le.0) then
        write(*,*) 'PSFinitally:ERROR: detection material must be >0'
        stop
      endif
      write(*,'(a)') 'PSF filename:'
      read(funit,'(1x,a30)') buffer
      write(*,'(a)') buffer
      psfunit = finduf()
      write(*,'(a)') 'Opening PSF as unit:'
      write(*,'(i2)') psfunit
      open(psfunit,file=buffer)
      write(psfunit,'(a)')
     &  '# kpar : e : x : y : z : u : v : w : wght : ilb(5) : DeltaN'

      !*** Init vars:
      nhist  = 0
      nhlast = 0
      npar   = 0

      read(funit,'(a80)') buffer
      if (index(buffer,eos).eq.0) then
        write(*,*) 'PSFinitally:ERROR: End-Of-Section mark not found'
        write(*,'(a,a)') '  expecting to find: ',eos
        write(*,'(a,a)') '  found instead:     ',buffer
        stop
      endif
      write(*,'(a)') '>>>> PSF tally initialization finished >>>>'
      end


! **** end of file *************************************************

