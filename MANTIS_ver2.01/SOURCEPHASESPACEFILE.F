!*******************************************************************
!*                         SOURCE                                  *
!*                   PHASE SPACE FILE                              *
!*                                                                 *
!* Short description:                                              *
!*   Generation of primary particle states for radiation transport *
!*   calculations with PENELOPE.                                   *
!*                                                                 *
!*   Initial particle states are read from an external phase space *
!*   file (PSF).                                                   *
!*                                                                 *
!* Dependencies:                                                   *
!*   from PENELOPE:                                                *
!*   -> common /TRACK/                                             *
!*   -> routine STORES                                             *
!*   from PENGEOM:                                                 *
!*   -> routines LOCATE,STEP                                       *
!*   from PENAUX.F                                                 *
!*   -> routine DOSTOP                                             *
!*   -> routine FINDUF                                             *
!*   from other penEasy libraries:                                 *
!*   -> routine TALLY                                              *
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


      subroutine PSFsource(n)
!*******************************************************************
!*    Input:                                                       *
!*      n -> top history counter                                   *
!*    Output:                                                      *
!*      through /track/ and sec stack                              *
!*      n -> top history counter                                   *
!*******************************************************************
      implicit none
      integer*4 n

      integer*4 kpar,ibody,mat,ilb
      real*8 e,x,y,z,u,v,w,wght
      common/track/e,x,y,z,u,v,w,wght,kpar,ibody,mat,ilb(5)
      integer*4 kpars,ilbs,dns
      real*8 es,xs,ys,zs,us,vs,ws,wghts
      common /srcpsf/ es,xs,ys,zs,us,vs,ws,wghts,kpars,ilbs(5),dns
      logical active
      integer in,split
      integer*4 nlin
      real*8 rot,xshift,yshift,zshift
      common /srcps1/ rot(3,3),xshift,yshift,zshift,nlin,split,in,
     &                active
      logical getpar
      integer i
      integer*4 ncross
      real*8 infty,dsef,norm
      parameter (infty=1.0d30)

      if (.not.active) return

      n = n-1+dns  ! Correct history counter from main program
      do
        ! Load previously stored particle into active common:
        kpar = kpars
        wght = wghts/split
        e = es
        ! Rotate and translate position:
        x = xshift+rot(1,1)*xs+rot(1,2)*ys+rot(1,3)*zs
        y = yshift+rot(2,1)*xs+rot(2,2)*ys+rot(2,3)*zs
        z = zshift+rot(3,1)*xs+rot(3,2)*ys+rot(3,3)*zs
        ! Rotate direction and renormalize to double precision:
        u = rot(1,1)*us+rot(1,2)*vs+rot(1,3)*ws
        v = rot(2,1)*us+rot(2,2)*vs+rot(2,3)*ws
        w = rot(3,1)*us+rot(3,2)*vs+rot(3,3)*ws
        norm = 1.0d0/sqrt(u**2+v**2+w**2)
        u = u*norm
        v = v*norm
        w = w*norm
        call locate  ! Sets ibody and mat from /track/ values
        if (mat.eq.0) call step(infty,dsef,ncross)  ! Where it is aiming at
        ilb(5) = ilbs(5)  ! Only ILB(5) is transferred through the PSF
        do i=1,split  ! particle splitting:
          call stores(e,x,y,z,u,v,w,wght,kpar,ilb)
          call tally(0,e)
        enddo

        ! Read a new particle and store for later calls:
        if (.not.getpar()) then
          call dostop  ! Forces simulation to stop
          write(*,*)  ' '
          write(*,'(a)')
     &      'PSFsource:WARNING: PSF exhausted; stopped forcefully'
          exit
        endif
        if (dns.ne.0) exit  ! While top primary is not exhausted
      enddo
      end


      subroutine PSFinisrc(funit,e0)
!*******************************************************************
!*    Initializes. To be called before SOURCE.                     *
!*                                                                 *
!*    Input:                                                       *
!*      funit -> input unit (must be already open)                 *
!*    Output:                                                      *
!*      e0 -> mean source energy (eV)                              *
!*******************************************************************
      implicit none
      integer funit
      real*8 e0

      integer*4 kpars,ilbs,dns
      real*8 es,xs,ys,zs,us,vs,ws,wghts
      common /srcpsf/ es,xs,ys,zs,us,vs,ws,wghts,kpars,ilbs(5),dns
      logical active
      integer in,split
      integer*4 nlin
      real*8 rot,xshift,yshift,zshift
      common /srcps1/ rot(3,3),xshift,yshift,zshift,nlin,split,in,
     &                active
      logical getpar
      character*80 psfnam,buffer
      character*(*) secid,eos
      parameter (secid=
     &'[SECTION SOURCE PHASE SPACE FILE v.2004-11-16]')
      parameter (eos='[END OF SECTION]')
      integer finduf,error,validate,idpresent
      integer*4 ntop,nele,npho,npos
      real*8 omega,theta,phi,comega,ctheta,cphi,somega,stheta,sphi
      real*8 dummy

      active = .true.
      write(*,*) ' '
      write(*,'(a)')
     & '>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'
      read(funit,'(a80)') buffer
      idpresent = index(buffer,secid)
      if (buffer(1:1).eq.'#'.and.idpresent.ne.0) then
        active = .false.
        e0 = -1.0d0
        write(*, '(a)')
     &    '>>>> Source Phase Space File is INACTIVE >>>>'
        do
          read(funit,'(a80)',end=10) buffer  ! Reads dummy lines of inactive source
          if (index(buffer,eos).ne.0) return  ! End-Of Section found
        enddo
 10     return  ! End-Of-File, something went wrong
      else if (idpresent.eq.0) then
        write(*,'(a)') 'PSFinisrc:ERROR: incorrect section;'
        write(*,'(a,a)') '  expecting to find: ',secid
        write(*,'(a,a)') '  found instead:     ',buffer
        stop
      endif
      write(*,'(a)') secid

      write(*,'(a)') 'PSF filename:'
      read(funit,'(1x,a30)') psfnam
      write(*,'(a)') psfnam

      write(*,'(a)') 'Splitting factor:'
      read(funit,*) split
      write(*,'(i10)') split
      if (split.lt.1) then
        write(*,*) 'PSFinisrc:ERROR: split < 1'
        stop
      endif

      write(*,'(a)')
     & 'Euler angles (deg) to rotate position and direction vectors:'
      read(funit,*) omega,theta,phi
      write(*,'(3(1x,1pe12.5))') omega,theta,phi
      ! Calculate rotation matrix:
      sphi   = sin(phi)
      cphi   = cos(phi)
      stheta = sin(theta)
      ctheta = cos(theta)
      somega = sin(omega)
      comega = cos(omega)
      rot(1,1) = cphi*ctheta*comega-sphi*somega
      rot(1,2) = -cphi*ctheta*somega-sphi*comega
      rot(1,3) = cphi*stheta
      rot(2,1) = sphi*ctheta*comega+cphi*somega
      rot(2,2) = -sphi*ctheta*somega+cphi*comega
      rot(2,3) = sphi*stheta
      rot(3,1) = -stheta*comega
      rot(3,2) = stheta*somega
      rot(3,3) = ctheta

      write(*,'(a)') 'Cartesian components of position shift (cm):'
      read(funit,*) xshift,yshift,zshift
      write(*,'(3(1x,1pe12.5))') xshift,yshift,zshift

      read(funit,*) validate
      if (validate.eq.1) then
        ! Pre-read PSF, validate and statistics:
        write(*,'(a)') 'Starting PSF validation'
        nlin = 0
        e0 = 0.0d0
        ntop = 0
        nele = 0
        npho = 0
        npos = 0
        in = finduf()
        open(in,file=psfnam,status='old',iostat=error)
        if (error.ne.0) then
          write(*,*) 'PSFinisrc:ERROR: cannot open the PSF'
          stop
        endif
        do
          if (.not.getpar()) exit  ! EOF reached
          ! Count particles in PSF:
          ntop = ntop+dns  ! Note that nlin is increased by GETPAR
          if (kpars.eq.1) then
            nele = nele+1
          else if (kpars.eq.2) then
            npho = npho+1
          else if (kpars.eq.3) then
            npos = npos+1
          else
            write(*,*) 'PSFinisrc:ERROR: invalid KPAR: in line:'
            write(*,'(2(1x,i10))') kpars,nlin
            stop
          endif
          if (es.lt.0.0d0.or.es.gt.1.0d9) then
            write(*,*)
     &        'PSFinisrc:ERROR: invalid energy(eV): in line:'
            write(*,'(1pe12.5,1x,i10)') es,nlin
            stop
          endif
          if (es.gt.e0) e0 = es
          if (us**2+vs**2+ws**2.lt.1.0d-30) then
            write(*,*)
     &        'PSFinisrc:ERROR: null vector direction found in line:'
            write(*,'(i10)') nlin
            stop
          endif
        enddo
        close(in)
        write(*,'(a)') 'PSF statistics:'
        write(*,'(a)') '  No. electrons:'
        write(*,'(2x,i10)') nele
        write(*,'(a)') '  No. photons:'
        write(*,'(2x,i10)') npho
        write(*,'(a)') '  No. positrons:'
        write(*,'(2x,i10)') npos
        write(*,'(a)') '  No. particles, total:'
        write(*,'(2x,i10)') nele+npho+npos
        write(*,'(a)') '  No. top primary histories:'
        write(*,'(2x,i10)') ntop
        write(*,'(a)') '  Max energy(eV):'
        write(*,'(2x,1pe12.5)') e0
        write(*,'(a)') ' (max energy declared in input file ignored)'
        read(funit,*) dummy
      else  ! Do not validate the PSF
        write(*,'(a)') '** User refused to pre-validate PSF **'
        write(*,'(a)') 'Max energy (eV) declared in input file:'
        read(funit,*) e0
        write(*,'(1pe12.5)') e0
      endif

      ! Prepare for 1st call to SOURCE:
      nlin = 0
      in = finduf()
      open(in,file=psfnam,status='old',iostat=error)
      if (error.ne.0) then
        write(*,*) 'PSFinisrc:ERROR: cannot open the PSF'
        stop
      endif
      if (.not.getpar()) then
        write(*,*) 'PSFinisrc:ERROR: PSF is empty'
        stop
      endif

      read(funit,'(a80)') buffer
      if (index(buffer,eos).eq.0) then
        write(*,*) 'PSFinisrc:ERROR: End-Of-Section mark not found'
        write(*,'(a,a)') '  expecting to find: ',eos
        write(*,'(a,a)') '  found instead:     ',buffer
        stop
      endif
      write(*,'(a)') '>>>> PSF source initialization finished >>>>'
      end


      logical function getpar()
!*******************************************************************
!*    Reads a new particle from the PSF.                           *
!*                                                                 *
!*    Output:                                                      *
!*      -> returns .false. if EOF has been reached, else .true.    *
!*      -> particle state in /srcpsf/                              *
!*******************************************************************
      implicit none
      integer*4 kpars,ilbs,dns
      real*8 es,xs,ys,zs,us,vs,ws,wghts
      common /srcpsf/ es,xs,ys,zs,us,vs,ws,wghts,kpars,ilbs(5),dns
      logical active
      integer in,split
      integer*4 nlin
      real*8 rot,xshift,yshift,zshift
      common /srcps1/ rot(3,3),xshift,yshift,zshift,nlin,split,in,
     &                active
      character*256 buffer
      integer error

      do
        read(in,'(a256)',end=10,iostat=error) buffer
        if (error.ne.0) then
          write(*,*)
     &      'getpar:ERROR: unable to read PSF line; last line read:'
          write(*,'(i10)') nlin
          stop
        endif
        nlin = nlin+1    ! lines read
        if (buffer(1:1).ne.'#') exit  ! A non-comment line found
      enddo

      read(buffer,*,iostat=error)
     &  kpars,es,xs,ys,zs,us,vs,ws,wghts,ilbs(5),dns
      if (error.ne.0) then
        write(*,*)
     &    'getpar:ERROR: invalid or missing datum in PSF line:'
        write(*,'(i10)') nlin
        write(*,*) '  line contents:'
        write(*,'(a)') buffer
        stop
      endif
      getpar = .true.
      return

 10   getpar = .false.  ! EOF code
      end


! **** end of file *************************************************

