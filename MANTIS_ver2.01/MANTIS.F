!*******************************************************************
!*                          MANTIS                                 *
!*                                                                 *
!* Short description:                                              *
!*   MANTIS is the coupling of penEasy (i.e. PENELOPE) + DETECT2.  *
!*                                                                 *
!* Dependencies:                                                   *
!*   from PENELOPE:                                                *
!*   -> common /TRACK/                                             *
!*   -> routines CLEANS,START,SECPAR                               *
!*   from PENGEOM:                                                 *
!*   -> routine STEP                                               *
!*   from PENVARED:                                                *
!*   -> routine JUMPF                                              *
!*   from other penEasy files:                                     *
!*   -> common /FRC002/                                            *
!*   -> routines in penaux.f, penpatch.f, sourceXX.f, tallyXX.f    *
!*      and timing.f                                               *
!*   from DETECT2                                                  *
!*   -> routines DETECT2REPORT,INIDETECT2                          *
!*        and, in tallyDETECT2.f,                                  *
!*   -> routines DETECT2,DETECT2EOH                                *
!*                                                                 *
!* Compatible with PENELOPE versions:                              *
!*   2005                                                          *
!*                                                                 *
!* Last update:                                                    *
!*   2005-08-01                                                    *
!*     Created.                                                    *
!*                                                                 *
!* Josep Sempau, email: josep.sempau@upc.es                        *
!*   Universitat Politecnica de Catalunya, Barcelona, Spain        *
!* Aldo Badano, email: agb@cdrh.fda.gov                            *
!*   Food and Drug Admin, Rockville, MD USA                        *
!*                                                                 *
!* SEE COPYRIGHT NOTICE IN FILE DISCLAIMER.PDF				 *
!*******************************************************************

!*******************************************************************
!*    Includes                                                     *
!*******************************************************************
      ! PENELOPE routines:
      include 'penelope.f'
      include 'pengeom.f'
      include 'penvared.f'

      ! Auxiliary routines:
      include 'penaux.f'
      ! include 'penpatch.f'
      include 'timing.f'

      ! Source models (see documentation for a detailed description):
      include 'sourceBoxIsotropicGaussSpectrum.f'
      include 'sourcePhaseSpaceFile.f'
      ! <you may add your own here>

      ! Tallies (see documentation for a detailed description):
      include 'tallySpatialDoseDistrib_c1.f'
      include 'tallySpatialDoseDistrib_c2.f'
      include 'tallySpatialDoseDistrib_c3.f'
      include 'tallyEnergyDepositionPulseHeightSpectrum.f'
      ! include 'tallyFluenceTrackLength.f'
      include 'tallyPhaseSpaceFile.f'
      include 'tallyParticleCurrentSpectrum.f'
      include 'tallyParticleTrackStructure.f'
      include 'tallyDetect2.f'
      ! <you may add your own here>


!*******************************************************************
!*    MAIN                                                         *
!*******************************************************************
      program main
      implicit none
      integer*4 kpar,ibody,mat,ilb
      real*8 e,x,y,z,u,v,w,wght
      common/track/e,x,y,z,u,v,w,wght,kpar,ibody,mat,ilb(5)
      logical forcing
      common /frc002/ forcing
      logical endsim,absorb
      integer*4 n,ncross,icol,left
      real*8 ds,dsef,de,uncert,cputime,dsmax

      write(*,*) ' '
      write(*,*)
     & '>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'
      write(*,*)
     & '>>>> This is MANTIS 2005 >>>>'
      write(*,*)
     & '>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'
      write(*,*) ' '
      call initime ! Writes date on the screen
      call init    ! Initializes the penEasy+PENELOPE+PENGEOM system
      call treset  ! Resets simulation timer
      n = 0        ! Resets history counter

      history: do              ! Each iteration simulates a new history
        n = n+1                ! Updates history counter
        call cleans            ! Empties the stack
        call source(n)         ! Puts primary particles (from the same history) in stack
        call tally(1,dble(n))  ! The simulation of this history begins

        particle: do                       ! Each iteration simulates a new particle
          call secpar(left)                ! Retrieves a particle from the stack
          if (left.eq.0) exit particle     ! Stack was empty
          call tally(-99,-e)               ! The simulation of this particle begins
          call start                       ! Resets transport mechanics
          call isforcing                   ! Determines if forcing will be applied (variance reduct)

          interact: do                     ! Each iteration simulates an interaction
            if (absorb()) exit interact    ! Checks particle absorption
            if (forcing) then
              call jumpf(dsmax(),ds)       ! Gets distance DS until next interaction
            else
              call jump(dsmax(),ds)        ! Interaction forcing (see PENELOPE manual)
            endif
            call step(ds,dsef,ncross)      ! Advance up to interaction point or interface
            if (ncross.eq.0) then
              call tally(3,ds)             ! Moved a distance DS, no interface crossed
            else
              call tally(4,dsef)           ! Moved a distance DSEF, interface found
              if (mat.eq.0) exit interact  ! New material is vacuum => gone
              call start                   ! New material => reset transport mechanics
              call isforcing               ! Determines if forcing will be applied (variance reduct)
              cycle interact
            endif
            if (forcing) then
              call knockf(de,icol)         ! Interaction forcing (see PENELOPE manual)
            else
              call knock(de,icol)          ! Simulates an interaction
            endif
            call tally(-int(icol),de)      ! Tallies kinetic energy released
          enddo interact
        enddo particle

        call tally(6,dble(n))            ! End-of-history bookkeeping
        if (endsim(n)) exit history      ! Simulation is finished
      enddo history

      call report(n,cputime(),1,uncert)  ! Writes report on the screen
      call endtime                       ! Writes timing on the screen
      stop
      end


      subroutine init
!*******************************************************************
!*    Initializes                                                  *
!*******************************************************************
      implicit none
      character*100 jobname
      common /cjob/ jobname

      character*100 filen
      integer funit,finduf,strlen
      real*8 emax

      write(*,*) ' '
      write(*,'(a)') '>>>> Init started >>>>'
      write(*,'(a)') 'init: Reading jobname from mantis.job'
      funit = finduf()

!    CHANGE TO READ IN COMMAND-LINE ARGUMENT
!      open(funit,file='mantis.job',status='old')
!      read(funit,'(1x,a30)') jobname
!      close(funit)
!      write(*,*) ' '
      call getarg(1,jobname)
!    END OF CHANGE
!     open(funit,file='mantis.job',status='old')
!     read(funit,'(1x,a30)') jobname
!     close(funit)
!     write(*,*) ' '
      write(*,'(a,a)') '  Jobname is: ',jobname

!    CHANGE TO READ PAR FILES 
      filen =
     &   jobname(:strlen(jobname))
	open(funit,file=trim(filen)//'.pen',status='old')
      call iniconfig(funit)      ! Simulation config
      call inisource(funit,emax) ! Source models
      call inipen(funit,emax)    ! Init PENELOPE+PENGEOM
      call initally(funit)       ! Tallies
      call iniforce(funit,emax)  ! Interaction forcing
      close(funit)
      write(*,'(a)') ' '
      write(*,'(a)') ' '
      write(*,'(a)')
     & '>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'
      write(*,'(a)')
     & '             *****************************'
      write(*,'(a)')
     & '             init: Initialization finished'
      write(*,'(a)')
     & '             *****************************'
      write(*,'(a)')
     & '>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'
      write(*,'(a)') ' '
      end


!*******************************************************************
!*******************************************************************
!*    Find below the source routines modifiable by users.          *
!*    Add your own models or delete the unwated ones following     *
!*    these rules:                                                 *
!*    -> All source models require: i) an initialization routine   *
!*       that must be called by INISOURCE (below); ii) a particle  *
!*       generation routine that must be called by SOURCE (see     *
!*       below).                                                   *
!*    -> Notice that the number and ordering of the source         *
!*       initialization routines has to be consistent with the     *
!*       sections included in the input file for penEasy.          *
!*******************************************************************
!*******************************************************************

      subroutine inisource(funit,e0)
!*******************************************************************
!*    Init routines for source models                              *
!*                                                                 *
!*    Input:                                                       *
!*      funit -> input unit (must be already open)                 *
!*    Output:                                                      *
!*      e0 -> max source energy (eV)                               *
!*******************************************************************
      implicit none
      integer funit
      real*8 e0

      integer active
      real*8 e

      e0 = 0.0d0  ! Will return the largest energy generated
      active = 0  ! Counts how many active sources are defined

      call BIGSinisrc(funit,e)
      if (e.ge.0.0d0) active = active+1
      if (e0.lt.e) e0 = e

      call PSFinisrc(funit,e)
      if (e.ge.0.0d0) active = active+1
      if (e0.lt.e) e0 = e

      if (active.lt.1) then
        write(*,*) 'inisrc:ERROR: No active sources defined'
        stop
      endif
      if (active.gt.1) then
        write(*,*)
     &    'inisrc:ERROR: There is more than one source active.'
        write(*,*)
     &    '  Please comment out the unused ones in the input file.'
        stop
      endif
      end


      subroutine source(n)
!*******************************************************************
!*    Source models                                                *
!*                                                                 *
!*    Input:                                                       *
!*      n -> top history counter                                   *
!*******************************************************************
      implicit none
      integer*4 n

      call BIGSsource(n)
      call PSFsource(n)
      end


!*******************************************************************
!*******************************************************************
!*    Find below the tally routines modifiable by users.           *
!*    Add your own tallies or delete the unwated ones following    *
!*    these rules:                                                 *
!*    -> All tally methods require: i) an initialization routine   *
!*       that must be called by INITALLY (below); ii) a tally      *
!*       routine that must be called by TALLY (below); iii) a      *
!*       reporting routine that must be called by REPORT (below).  *
!*    -> Notice that the number and ordering of the tally          *
!*       initialization routines has to be consistent with the     *
!*       sections included in the input file for penEasy.          *
!*******************************************************************
!*******************************************************************

      subroutine initally(funit)
!*******************************************************************
!*    Init routines for tallying                                   *
!*                                                                 *
!*    Input:                                                       *
!*      funit -> input unit (must be already open)                 *
!*******************************************************************
      implicit none
      integer funit

      call SDDinitally_c1(funit)
      call SDDinitally_c2(funit)
      call SDDinitally_c3(funit)
      call EPSinitally(funit)
      ! call FTLinitally  ! Not needed; removed for compilation simplicity
      call PSFinitally(funit)
      call PCSinitally(funit)
      call PTSinitally(funit)
      call iniDetect2		! Init Detect2 system
      end


      subroutine tally(mode,arg)
!*******************************************************************
!*    Tally routines                                               *
!*******************************************************************
      implicit none
      integer mode
      real*8 arg

      call SDDtally_c1(mode,arg)
      call SDDtally_c2(mode,arg)
      call SDDtally_c3(mode,arg)
      call EPStally(mode,arg)
      ! call FTLtally(mode,arg)
      call PSFtally(mode,arg)
      call PCStally(mode,arg)
      call PTStally(mode,arg)
      call detect2tally(mode,arg)
      end


      subroutine report(n,cputim,screen,uncert)
!*******************************************************************
!*    Reporting routines for each tally                            *
!*                                                                 *
!*    Input:                                                       *
!*      n -> no. of histories simulated                            *
!*      cputim -> elapsed CPU time                                 *
!*      screen -> 1 if output goes to screen, 0 else               *
!*    Output:                                                      *
!*      uncert -> overall uncertainty in % (at 2.sigma)            *
!*******************************************************************
      implicit none
      integer screen
      integer*4 n
      real*8 cputim,uncert
      real*8 unc

      uncert = 0.0d0  ! Will return max uncertainty
      call SDDreport_c1(n,cputim,screen,unc)
      ! if (uncert.lt.unc) uncert = unc
      call SDDreport_c2(n,cputim,screen,unc)
      ! if (uncert.lt.unc) uncert = unc
      call SDDreport_c3(n,cputim,screen,unc)
      ! if (uncert.lt.unc) uncert = unc
      call EPSreport(n,cputim,screen,unc)
      ! if (uncert.lt.unc) uncert = unc
      ! call FTLreport(n,cputim,screen,unc)
      ! if (uncert.lt.unc) uncert = unc
      call PSFreport(n,cputim,screen,unc)
      ! if (uncert.lt.unc) uncert = unc
      call PCSreport(n,cputim,screen,unc)
      ! if (uncert.lt.unc) uncert = unc
      call PTSreport(n,cputim,screen,unc)
      ! if (uncert.lt.unc) uncert = unc
      call detect2report(n,cputim,uncert)  ! This 'tally' determines the overall uncert
      end


!>>>> End Of File >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

