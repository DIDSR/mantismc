!*******************************************************************
!*                          TALLY                                  *
!*                 PARTICLE TRACK STRUCTURE                        *
!*                                                                 *
!* Short description:                                              *
!*   Tally routines for radiation transport calculations with      *
!*   PENELOPE.                                                     *
!*                                                                 *
!*   Writes to a file the position and value of energy loss events *
!*   so that particle tracks can be visualized.                    *
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
!*   2005-08-01 Created                                            *
!*                                                                 *
!*   Josep Sempau, email: josep.sempau@upc.es                      *
!*   Universitat Politecnica de Catalunya, Barcelona, Spain        *
!*******************************************************************


      subroutine PTStally(mode,arg)
!*******************************************************************
!*    Input:                                                       *
!*      mode -> Identifies the state of the calling routine        *
!*      arg -> energy loss (mode<0) or history no. (mode=1)        *
!*******************************************************************
      implicit none
      integer mode
      real*8 arg

      integer*4 kpar,ibody,mat,ilb
      real*8 e,x,y,z,u,v,w,wght
      common/track/e,x,y,z,u,v,w,wght,kpar,ibody,mat,ilb(5)
      logical active,report
      integer ptsunit
      integer*4 nhlast,ntrack,period
      real*8 xlast,ylast,zlast
      common /scopts/ xlast,ylast,zlast,
     &   nhlast,ntrack,period,ptsunit,report,active
      integer*4 n
      real*8 de

      if (.not.active) return

      if (mode.eq.1) then  ! Start of a new history
        n = arg+0.1d0
        if (n-nhlast.ge.period) then
          report = .true.  ! This history will be reported
          nhlast = n
          ntrack = ntrack+1
        else
          report = .false.
        endif
        return
      endif

      if (.not.report) return

      if (mode.lt.0) then
        de = arg  ! Energy scored (negative if particle retrieved from stack)
        if (mode.eq.-99) then  ! A new particle starts
          if (ilb(1).eq.1) de = 0.0d0 ! The particle is a primary
          write(ptsunit,*) ' '
        endif
        xlast = x
        ylast = y
        zlast = z

      else if (mode.eq.4) then
        if (mat.eq.0) then  ! The particle escaped
          de = 0.0d0
          x = xlast+arg*u  ! Move the particle up to the exiting surface
          y = ylast+arg*v
          z = zlast+arg*w
        else
          xlast = x
          ylast = y
          zlast = z
          return  ! No need to report
        endif

      else
        return  ! No need to report

      endif

      write(ptsunit,'(i1,2(1x,i6),6(1x,1pe12.5),1x,i10)')
     &  kpar,ibody,mat,x,y,z,e,de,wght,ilb(5)
      end


      subroutine PTSreport(n,cputim,screen,uncert)
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
      logical active,report
      integer ptsunit
      integer*4 nhlast,ntrack,period
      real*8 xlast,ylast,zlast
      common /scopts/ xlast,ylast,zlast,
     &   nhlast,ntrack,period,ptsunit,report,active
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
     &   jobname(:strlen(jobname))//'tallyParticleTrackStructure.dat'
        open(out,file=outfilen,iostat=error)
        if (error.eq.0) write(out,*,iostat=error) ' '
        if (error.ne.0) then
          write(*,*) ' '
          write(*,*) '********************************************'
          write(*,*)
     &      'PTSreport:ERROR: cannot write output data file;'
          write(*,*) '  writing to screen instead'
          write(*,*) '********************************************'
          write(*,*) ' '
          out = 6
        endif
      endif

      write(out,*) ' '
      write(out,'(a)')
     &'#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'
      write(out,'(a)') '# [SECTION REPORT PTS]'
      write(out,'(a)') '# No. of track histories written to PTS file:'
      write(out,'(i10)') ntrack

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


      subroutine PTSinitally(funit)
!*******************************************************************
!*    Initializes. To be called before TALLY                       *
!*                                                                 *
!*    Input:                                                       *
!*      funit -> input unit (must be already open)                 *
!*******************************************************************
      implicit none
      integer funit

      logical active,report
      integer ptsunit
      integer*4 nhlast,ntrack,period
      real*8 xlast,ylast,zlast
      common /scopts/ xlast,ylast,zlast,
     &   nhlast,ntrack,period,ptsunit,report,active
      character*80 buffer
      character*(*) secid,eos
      parameter (secid=
     &'[SECTION TALLY PARTICLE TRACK STRUCTURE v.2005-08-01]')
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
     &    '>>>> Tally Particle Track Structure is INACTIVE >>>>'
        do
          read(funit,'(a80)',end=10) buffer  ! Reads dummy lines of inactive tally
          if (index(buffer,eos).ne.0) return  ! End-Of Section found
        enddo
 10     return  ! End-Of-File, something went wrong
      else if (idpresent.eq.0) then
        write(*,'(a)') 'PTSinitally:ERROR: incorrect section;'
        write(*,'(a,a)') '  expecting to find: ',secid
        write(*,'(a,a)') '  found instead:     ',buffer
        stop
      endif
      write(*,'(a)') secid

      write(*,'(a)') 'Track data will be written to file:'
      read(funit,'(1x,a30)') buffer
      write(*,'(a)') buffer
      ptsunit = finduf()
      write(*,'(a)') 'Opening file as unit: (commented out by aldo) '
      write(*,'(i2)') ptsunit
      open(ptsunit,file=buffer)
      write(ptsunit,'(a)')
     &  '# kpar : body : mat : x : y : z : e : elost : wght : ilb(5)'
      write(*,'(a)') 'Writing one track in every:'
      read(funit,*) period
      write(*,'(i3)') period

      !*** Init vars:
      nhlast = 0
      ntrack = 0
      report = .false.

      read(funit,'(a80)') buffer
      if (index(buffer,eos).eq.0) then
        write(*,*) 'PTSinitally:ERROR: End-Of-Section mark not found'
        write(*,'(a,a)') '  expecting to find: ',eos
        write(*,'(a,a)') '  found instead:     ',buffer
        stop
      endif
      write(*,'(a)') '>>>> PTS tally initialization finished >>>>'
      end


! **** end of file *************************************************

