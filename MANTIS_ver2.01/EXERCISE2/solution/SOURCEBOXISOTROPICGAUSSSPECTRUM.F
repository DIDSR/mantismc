!*******************************************************************
!*                         SOURCE                                  *
!*               BOX-ISOTROPIC-GAUSS-SPECTRUM                      *
!*                                                                 *
!* Short description:                                              *
!*   Generation of primary particle states for radiation transport *
!*   calculations with PENELOPE.                                   *
!*                                                                 *
!* Dependencies:                                                   *
!*   from PENELOPE:                                                *
!*   -> common /TRACK/                                             *
!*   -> routines DIRECT,STORES,RAND                                *
!*   from PENGEOM:                                                 *
!*   -> routines LOCATE,STEP                                       *
!*   from other penEasy libraries:                                 *
!*   -> routine TALLY                                              *
!*                                                                 *
!* Compatible with PENELOPE versions:                              *
!*   2005                                                          *
!*                                                                 *
!* Last update:                                                    *
!*   2004-11-14                                                    *
!*     Adapted to MANTIS                                           *
!*******************************************************************


      subroutine BIGSsource(n)
!*******************************************************************
!*    Input:                                                       *
!*      n -> History no.                                           *
!*    Output:                                                      *
!*      through /track/ and sec stack                              *
!*******************************************************************
      use mtpar
      implicit none
      integer*4 n
      integer status

      integer*4 kpar,ibody,mat,ilb
      real*8 e,x,y,z,u,v,w,wght
      common/track/e,x,y,z,u,v,w,wght,kpar,ibody,mat,ilb(5)
      logical warned,srcpoint,active
      integer parsrc,matsrc,nspc,dim
      parameter (dim=10000)
      real*8 shots,usrc,vsrc,wsrc,cossrc,espc,pspc,despc,rot
      real*8 xsrc,ysrc,zsrc,dxsrc,dysrc,dzsrc
      common /srcbig/ rot(3,3),espc(dim),pspc(dim),despc(dim),shots,
     &                cossrc,usrc,vsrc,wsrc,xsrc,ysrc,zsrc,dxsrc,
     &                dysrc,dzsrc,parsrc,matsrc,nspc,warned,srcpoint,
     &                active
      integer ntrial,i,seeki
      integer*4 ncross
      real*8 costhe,phi,rand,infty,dsef,xrot,yrot,zrot
      real*8 pi,dospi,one3rd,randno,r
      parameter (pi=3.1415926535897932d0,dospi=2.0d0*pi,infty=1.0d30)
      parameter (one3rd=1.0d0/3.0d0)

      if (.not.active) return

      ntrial = 0
      kpar = parsrc    ! Particle type

      direction: do
        u = usrc
        v = vsrc
        w = wsrc
        status=rmt(r)
        costhe = cossrc+(1.0d0-cossrc)*r
        status=rmt(r)
        phi = dospi*r
        call direct(costhe,phi,u,v,w)

        position: do
          ntrial = ntrial+1
          if (n.eq.1.and.ntrial.gt.100) then
            write(*,*)
     &    'BIGSsource:ERROR: could not generate a valid particle '//
     &    'after 100 trials;'
            write(*,*)
     &       '  box enclosure is not in declared '//
     &       'material or particle direction is not correctly aimed'
            stop
          endif
          ! Sample coordinates relative to the box center:
          if (srcpoint) then  ! No need to sample for a point source
            x = xsrc
            y = ysrc
            z = zsrc
          else
            status=rmt(r)
            x = dxsrc*(r-0.5d0)
            status=rmt(r)
            y = dysrc*(r-0.5d0)
            status=rmt(r)
            z = dzsrc*(r-0.5d0)
            ! Rotate:
            xrot = rot(1,1)*x+rot(1,2)*y+rot(1,3)*z
            yrot = rot(2,1)*x+rot(2,2)*y+rot(2,3)*z
            zrot = rot(3,1)*x+rot(3,2)*y+rot(3,3)*z
            ! Shift to the box center position:
            x = xrot+xsrc
            y = yrot+ysrc
            z = zrot+zsrc
          endif
          call locate  ! Finds body and material (needs U,V,W)

          ! Accept or reject?
          if (matsrc.eq.0) then    ! Any material will do, even vacuum
            if (mat.eq.0) call step(infty,dsef,ncross)
            exit direction         ! Always accept
          else if (srcpoint) then  ! Particle must aim at specified material
            if (mat.eq.matsrc) exit direction
            do
              call step(infty,dsef,ncross)
              if (mat.eq.0.or.ncross.eq.0) cycle direction  ! Not aiming at MATSRC
              if (mat.eq.matsrc) exit direction  ! Accepted
            enddo
          else if (mat.eq.matsrc) then  ! Materials match
            exit direction
          endif                    ! Re-sample position if materials do not match
        enddo position
      enddo direction

      ! Kinetic energy:
      if (nspc.eq.2.and.pspc(2).eq.0.0d0) then  ! Gaussian, sample using Box-Muller
        do
          status=rmt(r)
          randno=r
          status=rmt(r)
          e = espc(2)+espc(1)*  ! These two values were set to Emean and sigma respect.
     &         sqrt(-2.0d0*log(randno))*sin(dospi*r)
          if (e.gt.0.0d0) exit
        enddo
      else    ! User-defined spectrum
        status=rmt(r)
        randno = r
        i = seeki(pspc,randno,nspc)
        e = espc(i)+(randno-pspc(i))*despc(i)
      endif

      ! Average no. of shots per call:
      shots = shots+ntrial
      if (.not.warned.and.shots.gt.n*10.0d0) then
        write(*,*) ' '
        write(*,'(a)')
     &    '***************'
        write(*,'(a)')
     &    'BIGSsource:WARNING: source effectiveness is too low !'
        write(*,'(a)')
     &    '  Redefine source parameters appropiately.'
        write(*,'(a)')
     &    '  Histories and shots per source call so far:'
        write(*,'(2x,i10,1x,1pe12.5)') n,shots/n
        write(*,'(a)')
     &    '***************'
        write(*,*) ' '
        warned = .true.
      endif

      wght = 1.0d0
      ilb(1) = 1    ! Particle is 1st generation (a primary)
      ilb(5) = 0    ! Optional label (trasferred to descendants)

      call stores(e,x,y,z,u,v,w,wght,kpar,ilb)  ! Pushes particle to stack
      call tally(0,e)                           ! Deposits its kinetic E
      end


      subroutine BIGSinisrc(funit,e0)
!*******************************************************************
!*    Initializes the source.                                      *
!*                                                                 *
!*    Input:                                                       *
!*      funit -> input unit (must be already open)                 *
!*    Output:                                                      *
!*      e0 -> max source energy (eV) (Emean for Gaussian)          *
!*******************************************************************
      implicit none
      integer funit
      real*8 e0

      logical warned,srcpoint,active
      integer parsrc,matsrc,nspc,dim
      parameter (dim=10000)
      real*8 shots,usrc,vsrc,wsrc,cossrc,espc,pspc,despc,rot
      real*8 xsrc,ysrc,zsrc,dxsrc,dysrc,dzsrc
      common /srcbig/ rot(3,3),espc(dim),pspc(dim),despc(dim),shots,
     &                cossrc,usrc,vsrc,wsrc,xsrc,ysrc,zsrc,dxsrc,
     &                dysrc,dzsrc,parsrc,matsrc,nspc,warned,srcpoint,
     &                active
      character*80 buffer
      character*(*) secid,eos
      parameter (secid=
     &'[SECTION SOURCE BOX ISOTROPIC GAUSS SPECTRUM v.2004-11-13]')
      parameter (eos='[END OF SECTION]')
      integer j,idpresent
      real*8 omega,theta,phi,comega,somega,ctheta,stheta,cphi,sphi
      real*8 prob,pi,norm,fwhm2sig,emean,sigma,deg2rad
      parameter (pi=3.1415926535897932d0,deg2rad=pi/180.0d0)
      parameter (fwhm2sig=4.2466090014400952d-1)

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
     &    '>>>> Source Box Isotropic Gauss Spectrum is INACTIVE >>>>'
        do
          read(funit,'(a80)',end=10) buffer  ! Reads dummy lines of inactive source
          if (index(buffer,eos).ne.0) return  ! End-Of Section found
        enddo
 10     return  ! End-Of-File, something went wrong
      else if (idpresent.eq.0) then
        write(*,'(a)') 'BIGSinisrc:ERROR: incorrect section;'
        write(*,'(a,a)') '  expecting to find: ',secid
        write(*,'(a,a)') '  found instead:     ',buffer
        stop
      endif
      write(*,'(a)') secid

      ! Type of particle:
      write(*,'(a)') 'Particle type:'
      read(funit,*) parsrc
      write(*,'(i1)') parsrc
      if (parsrc.ne.1.and.parsrc.ne.2.and.parsrc.ne.3) then
        write(*,*) 'BIGSinisrc:ERROR: invalid particle type'
        stop
      endif

      ! Energy:
      write(*,'(a)') 'Reading energy spectrum'
      read(funit,'(a80)') buffer
      write(*,'(a)') '  Energy(eV)   Probability   Channel#'
      nspc = 0
      pspc(1) = 0.0d0
      do
        nspc = nspc+1
        read(funit,*) espc(nspc),prob
        write(*,'(2(1x,1pe12.5),1x,i5)') espc(nspc),prob,nspc
        if (espc(nspc).lt.espc(max(nspc-1,1))) then
          write(*,*) 'BIGSinisrc:ERROR: decreasing energy'
          stop
        endif
        if (prob.lt.0.0d0) exit    ! End of spectrum
        if (nspc.ge.dim) then
          write(*,*) 'BIGSinisrc:ERROR: too many bins in spectrum;'
          write(*,*) '              enlarge DIM'
          stop
        endif
        pspc(nspc+1) = pspc(nspc)+prob
      enddo
      write(*,'(a)') 'No. of channels read:'
      write(*,'(i5)') nspc-1
      if (nspc.lt.2) then
        write(*,*)
     &    'BIGSinisrc:ERROR: at least 1 channel must be defined'
        stop
      endif
      if (pspc(nspc).gt.0.0d0) then
        write(*,'(a)')
     &    'Sum of relative probabilities before normalization:'
        write(*,'(1pe12.5)') pspc(nspc)
        do j=1,nspc    ! Normalize to unity
          pspc(j) = pspc(j)/pspc(nspc)
        enddo
      else
        write(*,'(a)')
     &    'Sum of probabilities was nul; assuming Gaussian spectrum'
        if (nspc.gt.2) then
          write(*,*)
     &      'BIGSinisrc:ERROR: Gaussian requires only 1 channel'
          stop
        endif
        emean = 0.5d0*(espc(1)+espc(2))
        sigma = fwhm2sig*(espc(2)-espc(1))
        write(*,'(a)')
     &    'Mean energy, FWHM and sigma (all in eV) of the Gaussian:'
        write(*,'(3(1x,1pe12.5))') emean,espc(2)-espc(1),sigma
        if (emean.lt.0.0d0) then
          write(*,*) 'BIGSinisrc:ERROR: negative mean energy'
          stop
        endif
        espc(2) = emean
        espc(1) = sigma
      endif
      e0 = espc(nspc)  ! Set max (or Gaussian mean) energy
      do j=1,nspc-1    ! Prepare auxiliar array for sampling
        despc(j) = 0.0d0
        if (pspc(j+1).gt.pspc(j))
     &    despc(j) = (espc(j+1)-espc(j))/(pspc(j+1)-pspc(j))
      enddo

      ! Position:
      write(*,'(a)') 'Center coordinates (cm):'
      read(funit,*) xsrc,ysrc,zsrc
      write(*,'(3(1x,1pe12.5))') xsrc,ysrc,zsrc
      write(*,'(a)') 'Box sides (cm):'
      read(funit,*) dxsrc,dysrc,dzsrc
      write(*,'(3(1x,1pe12.5))') dxsrc,dysrc,dzsrc
      write(*,'(a)') 'Euler angles (deg):'
      read(funit,*) omega,theta,phi
      write(*,'(3(1x,1pe12.5))') omega,theta,phi
      ! Rotation matrix:
      somega = sin(omega*deg2rad)
      comega = cos(omega*deg2rad)
      stheta = sin(theta*deg2rad)
      ctheta = cos(theta*deg2rad)
      sphi = sin(phi*deg2rad)
      cphi = cos(phi*deg2rad)
      rot(1,1) = cphi*ctheta*comega-sphi*somega
      rot(1,2) = -cphi*ctheta*somega-sphi*comega
      rot(1,3) = cphi*stheta
      rot(2,1) = sphi*ctheta*comega+cphi*somega
      rot(2,2) = -sphi*ctheta*somega+cphi*comega
      rot(2,3) = sphi*stheta
      rot(3,1) = -stheta*comega
      rot(3,2) = stheta*somega
      rot(3,3) = ctheta
      ! Material:
      if (dxsrc.eq.0.0d0.and.dysrc.eq.0.0d0.and.dzsrc.eq.0.0d0)
     &  srcpoint = .true.  ! Needed to interpret source material
      write(*,'(a)') 'Source material:'
      read(funit,*) matsrc
      write(*,'(i2)') matsrc
      if (srcpoint.and.matsrc.ne.0) write(*,'(a)')
     &  '  (interpreted as the material particles must be aiming at)'
      if (matsrc.lt.0) then
        write(*,*) 'BIGSinisrc:ERROR: negative materials are invalid'
        stop
      endif

      ! Direction:
      write(*,'(a)') 'Direction vector (u,v,w):'
      read(funit,*) usrc,vsrc,wsrc
      write(*,'(3(1x,1pe12.5))') usrc,vsrc,wsrc
      write(*,'(a)') 'Aperture (deg) and cos():'
      read(funit,*) theta
      cossrc = cos(theta*pi/180.0d0)
      write(*,'(2(1x,1pe12.5))') theta,cossrc
      if (theta.lt.0.0d0.or.theta.gt.180.0d0) then
        write(*,*) 'BIGSinisrc:ERROR: aperture must be in [0,180].'
        stop
      endif
      norm = sqrt(usrc**2+vsrc**2+wsrc**2)
      if (norm.lt.1.0d-30) then
        if (theta.ne.180.0d0) then
          write(*,*)
     &  'BIGSinisrc:ERROR: null direction only valid if aperture=180'
          stop
        endif
        wsrc = 1.0d0
      else
        usrc = usrc/norm
        vsrc = vsrc/norm
        wsrc = wsrc/norm
      endif

      ! Init performance vars:
      shots = 0.0d0
      warned = .false.

      read(funit,'(a80)') buffer
      if (index(buffer,eos).eq.0) then
        write(*,*) 'BIGSinisrc:ERROR: End-Of-Section mark not found'
        write(*,'(a,a)') '  expecting to find: ',eos
        write(*,'(a,a)') '  found instead:     ',buffer
        stop
      endif
      write(*,'(a)') '>>>> BIGS source initialization finished >>>>'
      end


      integer function seeki(x,xc,n)
!*******************************************************************
!*    Finds the interval (x(i),x(i+1)] containing the value xc.    *
!*                                                                 *
!*    Input:                                                       *
!*      x(1..n) -> data array                                      *
!*      xc -> point to be located                                  *
!*      n -> no. of data points                                    *
!*    Output:                                                      *
!*      index i of the semiopen interval where xc lies             *
!*    Comments:                                                    *
!*      -> If xc=x(1) then i=1 is returned.                        *
!*      -> If xc is outside the closed interval [x(1),x(n)]  the   *
!*         execution is aborted.                                   *
!*******************************************************************
      implicit none
      integer n
      real*8 xc,x(n)

      integer itop,imid

      if(xc.gt.x(n)) then
        write(*,*) 'seeki error: value outside range, xc>x(n):'
        write(*,*) xc,x(n)
        stop
      endif
      if(xc.lt.x(1)) then
        write(*,*) 'seeki error: value outside range, xc<x(1):'
        write(*,*) xc,x(1)
        stop
      endif

      seeki = 1
      itop = n
      do
        imid = (seeki+itop)/2
        if(xc.gt.x(imid)) then
          seeki = imid
        else
          itop = imid
        endif
        if(itop-seeki.le.1) exit
      enddo
      end


! **** end of file *************************************************

