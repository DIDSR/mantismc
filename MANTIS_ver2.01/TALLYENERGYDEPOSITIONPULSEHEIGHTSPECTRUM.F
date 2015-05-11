!*******************************************************************
!*                          TALLY                                  *
!*             ENERGY DEPOSITION (COLLISION ESTIMATOR)             *
!*                AND PULSE HEIGHT SPECTRUM                        *
!*                                                                 *
!* Short description:                                              *
!*   Tally routines for radiation transport calculations with      *
!*   PENELOPE.                                                     *
!*                                                                 *
!*   Calculates the spectrum of the energy deposited in a given    *
!*   material and the total energy deposition with the collision   *
!*   estimator.                                                    *
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
!*											 *
!* SEE COPYRIGHT NOTICE IN FILE DISCLAIMER.PDF				 *
!*******************************************************************



      subroutine EPStally(mode,eloss)
!*******************************************************************
!*    Input:                                                       *
!*      mode -> Identifies the state of the calling procedure      *
!*      eloss -> energy deposition                                 *
!*******************************************************************
      implicit none
      integer mode
      real*8 eloss

      integer*4 kpar,ibody,mat,ilb
      real*8 e,x,y,z,u,v,w,wght
      common/track/e,x,y,z,u,v,w,wght,kpar,ibody,mat,ilb(5)
      logical active
      integer detmat,maxbin,nbin
      parameter (maxbin=1000)
      real*8 edptmp,count,emin,ebin,iebin,edep,edep2
      common /scoeps/ count(maxbin),edptmp,edep,edep2,emin,ebin,
     &                iebin,nbin,detmat,active
      integer bin

      if (.not.active) return

      if (mode.le.0) then
        if (mat.ne.detmat) return
        edptmp = edptmp+eloss*wght

      else if (mode.eq.6) then
        edep   = edep +edptmp
        edep2  = edep2+edptmp**2
        bin = (edptmp-emin)*iebin+1.0d0
        edptmp = 0.0d0
        if (bin.lt.1.or.bin.gt.nbin) return
        count(bin) = count(bin)+1.0d0

      endif
      end


      subroutine EPSreport(n,cputim,screen,maxuncert)
!*******************************************************************
!*    Input:                                                       *
!*      n -> no. of histories simulated                            *
!*      cputim -> elapsed CPU time                                 *
!*      screen -> 1 if output goes to screen, 0 else               *
!*    Output:                                                      *
!*      maxuncert -> max overall uncertainty in % (at 2.sigma)     *
!*    Comments:                                                    *
!*      -> 'cputim' should not include initialization procedures;  *
!*         enter 0 or neg. if not available.                       *
!*******************************************************************
      implicit none
      integer screen
      integer*4 n
      real*8 cputim,maxuncert

      integer*4 seed1,seed2
      common/rseed/seed1,seed2
      logical active
      integer detmat,maxbin,nbin
      parameter (maxbin=1000)
      real*8 edptmp,count,emin,ebin,iebin,edep,edep2
      common /scoeps/ count(maxbin),edptmp,edep,edep2,emin,ebin,
     &                iebin,nbin,detmat,active
      character*100 jobname
      common /cjob/ jobname
      character*100 outfilen
      integer i,nchan,out,finduf,error,strlen
      real*8 q,q2,sigma,eff,avesig,maxq,emiddle,elow,uncert

      maxuncert = 0.0d0  ! Needed in case of inactivity
      if (.not.active) return

!     *************************************
!     *** Total energy deposited report ***
!     *************************************

!     *** Prepare output files:
      out = 6    ! Default is to write to screen
      if (screen.eq.0) then
        out = finduf()
        outfilen =
     &    jobname(:strlen(jobname))//'tallyEnergyDeposition.dat'
        open(out,file=outfilen,iostat=error)
        if (error.ne.0) then
          write(*,*) ' '
          write(*,*) '*******************************************'
          write(*,*) 'EPSreport:ERROR: cannot open output data file;'
          write(*,*) '  writing to screen instead'
          write(*,*) '*******************************************'
          write(*,*) ' '
          out = 6
        endif
      endif

      write(out,*) ' '
      write(out,'(a)')
     &'#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'
      write(out,'(a)') '# [SECTION REPORT ENERGY DEPOSITION]'
      write(out,'(a)') '# Total Energy Deposition(eV) : +-2sigma'
      q = edep/n
      q2 = edep2/n
      sigma = (q2-q**2)/n
      sigma = sqrt(max(sigma,0.0d0))
      write(out,'(1pe12.5,1x,1pe8.1)') q,2.0d0*sigma
!     *** Evaluate rel. uncertainty:
      uncert = 100.0d0
      if (q.gt.0.0d0) uncert = 200.0d0*sigma/q
      maxuncert = max(maxuncert,uncert)

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
      write(out,'(a)')
     & '#   Average uncertainty (above 1/2 max score) in % [uncert]:'
      write(out,'(a,1pe12.5)') '#   ',uncert
      eff = n*uncert**2
      if (eff.gt.0.0d0) then
        write(out,'(a)') '#   Intrinsic efficiency [N*uncert^2]^-1:'
        write(out,'(a,1pe12.5)') '#   ',1.0d0/eff
      endif
      eff = cputim*uncert**2
      if (eff.gt.0.0d0) then
        write(out,'(a)') '#   Absolute efficiency [t*uncert^2]^-1:'
        write(out,'(a,1pe12.5)') '#   ',1.0d0/eff
      endif
      write(out,'(a)') '#'
      write(out,'(a)') '# Have a nice day.'
      write(out,'(a)') ' '

      if (out.ne.6) close(out)

!     ***********************
!     *** Spectrum report ***
!     ***********************

!     *** Prepare output files:
      out = 6    ! Default is to write to screen
      if (screen.eq.0) then
        out = finduf()
        outfilen =
     &    jobname(:strlen(jobname))//'tallyPulseHeightSpectrum.dat'
        open(out,file=outfilen,iostat=error)
        if (error.ne.0) then
          write(*,*) ' '
          write(*,*) '*******************************************'
          write(*,*) 'EPSreport:ERROR: cannot open output data file;'
          write(*,*) '  writing to screen instead'
          write(*,*) '*******************************************'
          write(*,*) ' '
          out = 6
        endif
      endif

!     *** Evaluate 1/2 of the max score:
      avesig = 0.0d0
      nchan = 0
      maxq = 0.0d0
      do i=1,nbin
        maxq = max(maxq,count(i))
      enddo
      maxq = 0.5d0*maxq

      write(out,'(a)') ' '
      write(out,'(a)')
     &'#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'
      write(out,'(a)') '# [SECTION REPORT PULSE HEIGHT SPECTRUM]'
      write(out,'(a)') '#'
      write(out,'(a)')
     &  '# Two energies per row are given, namely,'
      write(out,'(a)')
     &  '#   the low end and the middle point of each channel'
      write(out,'(a)') '#'
      write(out,'(a)')
     &  '# E_low(eV) : E_middle(eV) : counts(1/eV) : +-2sigma'
      do i=1,nbin
        q = count(i)/n
        sigma = q*(1.0d0-q)/n
        sigma = sqrt(max(sigma,0.0d0))
        q = q*iebin
        sigma = sigma*iebin
        elow = emin+ebin*(i-1)
        emiddle = elow+ebin*0.5d0
        write(out,'(3(1x,1pe12.5),1x,1pe8.1)')
     &    elow,emiddle,q,2.0d0*sigma
!       *** Evaluate average uncertainty for scores above 1/2 max score:
        if (count(i).gt.maxq) then
          avesig = avesig+(sigma/q)**2
          nchan = nchan+1
        endif
      enddo
      uncert = 100.0d0
      if (nchan.gt.0) uncert = 200.0d0*sqrt(avesig/nchan)
      maxuncert = max(maxuncert,uncert)

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
      write(out,'(a)')
     & '#   Average uncertainty (above 1/2 max score) in % [uncert]:'
      write(out,'(a,1pe12.5)') '#   ',uncert
      eff = n*uncert**2
      if (eff.gt.0.0d0) then
        write(out,'(a)') '#   Intrinsic efficiency [N*uncert^2]^-1:'
        write(out,'(a,1pe12.5)') '#   ',1.0d0/eff
      endif
      eff = cputim*uncert**2
      if (eff.gt.0.0d0) then
        write(out,'(a)') '#   Absolute efficiency [t*uncert^2]^-1:'
        write(out,'(a,1pe12.5)') '#   ',1.0d0/eff
      endif
      write(out,'(a)') '#'
      write(out,'(a)') '# Have a nice day.'
      write(out,'(a)') ' '

      if (out.ne.6) close(out)
      end


      subroutine EPSinitally(funit)
!*******************************************************************
!*    Initializes. To be called before TALLY.                      *
!*                                                                 *
!*    Input:                                                       *
!*      funit -> input unit (must be already open)                 *
!*******************************************************************
      implicit none
      integer funit

      logical active
      integer detmat,maxbin,nbin
      parameter (maxbin=1000)
      real*8 edptmp,count,emin,ebin,iebin,edep,edep2
      common /scoeps/ count(maxbin),edptmp,edep,edep2,emin,ebin,
     &                iebin,nbin,detmat,active
      character*80 buffer
      character*(*) secid,eos
      parameter (secid=
     &'[SECTION TALLY ENERGY DEPOSITION PULSE SPECTRUM v.2003-04-01]')
      parameter (eos='[END OF SECTION]')
      integer i,idpresent
      real*8 emax

      active = .true.
      write(*,*) ' '
      write(*,'(a)')
     & '>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'
      read(funit,'(a80)') buffer
      idpresent = index(buffer,secid)
      if (buffer(1:1).eq.'#'.and.idpresent.ne.0) then
        active = .false.
        write(*, '(a)')
     &'>>>> Tally Energy Deposition Pulse Height Spectrum is INACTIVE'
        do
          read(funit,'(a80)',end=10) buffer  ! Reads dummy lines of inactive tally
          if (index(buffer,eos).ne.0) return  ! End-Of Section found
        enddo
 10     return  ! End-Of-File, something went wrong
      else if (idpresent.eq.0) then
        write(*,'(a)') 'EPSinitally:ERROR: incorrect section;'
        write(*,'(a,a)') '  expecting to find: ',secid
        write(*,'(a,a)') '  found instead:     ',buffer
        stop
      endif
      write(*,'(a)') secid

      write(*,'(a)') 'Detection material set to:'
      read(funit,*) detmat
      write(*,'(i3)') detmat
      write(*,'(a)')
     &  'Spectrum limits [Emin,Emax] (eV) and no. channels:'
      read(funit,*) emin,emax,nbin
      write(*,'(2(1x,1pe12.5),1x,i5)') emin,emax,nbin
      ebin = (emax-emin)/nbin
      iebin = 1.0d0/ebin
!     *** Add one bin to store, separately, counts with E=Emax:
      nbin = nbin+1
      if (nbin.gt.maxbin) then
        write(*,*)
     &    'EPSinitally:ERROR: Too many channels; increase MAXBIN'
        stop
      endif

!     *** Clear counters:
      edptmp = 0.0d0
      edep   = 0.0d0
      edep2  = 0.0d0
      do i=1,nbin
        count(i) = 0.0d0
      enddo

      read(funit,'(a80)') buffer
      if (index(buffer,eos).eq.0) then
        write(*,*) 'EPSinitally:ERROR: End-Of-Section mark not found'
        write(*,'(a,a)') '  expecting to find: ',eos
        write(*,'(a,a)') '  found instead:     ',buffer
        stop
      endif
      write(*,'(a)') '>>>> EPS tally initialization finished >>>>'
      end


! **** end of file *************************************************

