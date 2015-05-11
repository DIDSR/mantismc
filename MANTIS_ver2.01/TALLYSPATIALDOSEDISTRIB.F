!*******************************************************************
!*                          TALLY                                  *
!*             SPATIAL (3D) DOSE DISTRIBUTION                      *
!*                                                                 *
!* Short description:                                              *
!*   Tally routines for radiation transport calculations with      *
!*   PENELOPE.                                                     *
!*                                                                 *
!*   Dose distribution in volume elements (voxels) superimposed    *
!*   on a solid-body PENGEOM geometry.                             *
!*                                                                 *
!* Dependencies:                                                   *
!*   from PENELOPE:                                                *
!*   -> common /TRACK/                                             *
!*   -> common /RSEED/                                             *
!*   -> common /COMPOS/                                            *
!*   from PENGEOM:                                                 *
!*   -> routine LOCATE                                             *
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



      subroutine SDDtally(mode,arg)
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
      logical active
      integer prtvox
      integer*4 maxvox
      parameter(maxvox=1000000)
      integer*4 nx,ny,nz,nynz,nlast,nhist
      real*8 edptmp,edep,edep2,idens,xmin,ymin,zmin
      real*8 dx,dy,dz,idx,idy,idz
      common /scosdd/ edptmp(maxvox),edep(maxvox),edep2(maxvox),
     &                idens(maxvox),nlast(maxvox),xmin,ymin,zmin,
     &                dx,dy,dz,idx,idy,idz,nx,ny,nz,nynz,nhist,
     &                prtvox,active
      integer*4 vox,getabs,i,j,k
      real*8 eps
      parameter(eps=1.0d-6)

      if (.not.active) return

      if (mode.le.0) then
        if (arg.eq.0.0d0) return    ! Nothing to deposit
!       *** Check if particle is inside tally region:
        i = (x-xmin)*idx+1.0d0
        if (i.lt.1.or.i.gt.nx) return
        j = (y-ymin)*idy+1.0d0
        if (j.lt.1.or.j.gt.ny) return
        k = (z-zmin)*idz+1.0d0
        if (k.lt.1.or.k.gt.nz) return
        vox = getabs(i,j,k)
!       *** Transfer partial tally to totals only when a new history visits:
        if (nhist.ne.nlast(vox)) then
          edep(vox)  = edep(vox) +edptmp(vox)
          edep2(vox) = edep2(vox)+edptmp(vox)**2
          edptmp(vox)= arg*wght
          nlast(vox) = nhist
        else
          edptmp(vox) = edptmp(vox)+arg*wght
        endif

      else if (mode.eq.1) then
        nhist = arg+eps

      endif
      end


      subroutine SDDreport(n,cputim,screen,uncert)
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
      integer prtvox
      integer*4 maxvox
      parameter(maxvox=1000000)
      integer*4 nx,ny,nz,nynz,nlast,nhist
      real*8 edptmp,edep,edep2,idens,xmin,ymin,zmin
      real*8 dx,dy,dz,idx,idy,idz
      common /scosdd/ edptmp(maxvox),edep(maxvox),edep2(maxvox),
     &                idens(maxvox),nlast(maxvox),xmin,ymin,zmin,
     &                dx,dy,dz,idx,idy,idz,nx,ny,nz,nynz,nhist,
     &                prtvox,active
      character*100 jobname
      common /cjob/ jobname
      character*100 outfilen
      integer nchan,out,finduf,error,strlen
      integer*4 vox,getabs,i,j,k,nxef,nyef,nzef
      real*8 q,q2,sigma,eff,avesig,maxq,fact,x,y,z
      real*8 xminef,yminef,zminef,xmiddle,ymiddle,zmiddle

      uncert = 0.0d0  ! Needed in case of inactivity
      if (.not.active) return

!     *** Prepare output files:
      out = 6    ! Default is to write to screen
      if (screen.eq.0) then
        out = finduf()
        outfilen =
     &    jobname(:strlen(jobname))//'tallySpatialDoseDistrib.dat'
        open(out,file=outfilen,iostat=error)
        if (error.ne.0) then
          write(*,*) ' '
          write(*,*) '*******************************************'
          write(*,*) 'SDDreport:ERROR: cannot open output data file;'
          write(*,*) '  writing to screen instead'
          write(*,*) '*******************************************'
          write(*,*) ' '
          out = 6
        endif
      endif

!     *** Dump counters and obtain max score:
      avesig = 0.0d0
      nchan = 0
      maxq = 0.0d0
      do vox=1,nx*ny*nz
        if (nlast(vox).ne.0) then
          edep(vox)  = edep(vox) +edptmp(vox)
          edep2(vox) = edep2(vox)+edptmp(vox)**2
          edptmp(vox)= 0.0d0
          nlast(vox) = 0
        endif
        maxq = max(maxq,edep(vox))    ! 1/2 of the max score
      enddo
      maxq = 0.5d0*maxq

!     *** Write info:
      write(out,'(a)') ' '
      write(out,'(a)')
     &'#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'
      write(out,'(a)') '# [SECTION REPORT SPATIAL DOSE DISTRIB]'
      write(out,'(a)') '# Dose is in eV*cm3/(g*[dx]*[dy]*[dz])'
      write(out,'(a)')
     &'#   [dx] (in cm) is present only if nxBin is not null;'
      write(out,'(a)')
     &'#   the same applies to [dy],[dz] with respect to nyBin,nzBin'

      nxef = nx
      xminef = xmin
      if (dx.eq.0.0d0) then
        nxef = 0
        xminef = 0.0d0
      endif
      nyef = ny
      yminef = ymin
      if (dy.eq.0.0d0) then
        nyef = 0
        yminef = 0.0d0
      endif
      nzef = nz
      zminef = zmin
      if (dz.eq.0.0d0) then
        nzef = 0
        zminef = 0.0d0
      endif
      write(out,'(a)') '#   nxBin,nyBin,nzBin:'
      write(out,'(a,3(1x,i5))') '#   ',nxef,nyef,nzef
      write(out,'(a)') '#   xmin,dx,ymin,dy,zmin,dz (cm):'
      write(out,'(a,3(2x,1pe12.5,1x,1pe12.5))')
     &  '#   ',xminef,dx,yminef,dy,zminef,dz

!     *** Write header:
      write(out,'(a)') '#'
      if (prtvox.eq.1) then
        write(out,'(a)')
     &    '# Two values per bin are given, namely,'
        write(out,'(a)')
     &    '#   the low end and the middle point of each bin'
      endif
      write(out,'(a,$)') '# '
      if (prtvox.eq.1) then
        if (dx.gt.0.0d0)
     &    write(out,'(a,$)') 'xBinIndex : xLow(cm) : xMiddle(cm) : '
        if (dy.gt.0.0d0)
     &    write(out,'(a,$)') 'yBinIndex : yLow(cm) : yMiddle(cm) : '
        if (dz.gt.0.0d0)
     &    write(out,'(a,$)') 'zBinIndex : zLow(cm) : zMiddle(cm) : '
      endif
      write(out,'(a)') 'dose : +-2sigma'

!     *** Write data:
      do k=1,nz
        do j=1,ny
          do i=1,nx
            vox = getabs(i,j,k)
            fact = idens(vox)
            if (dx.gt.0.0d0) then
              fact = fact*idx
              if (prtvox.eq.1) then
                x = xmin+dx*(i-1)
                xmiddle = x+dx*0.5d0
                write(out,'(1x,i5,2(1x,1pe12.5),$)') i,x,xmiddle
              endif
            endif
            if (dy.gt.0.0d0) then
              fact = fact*idy
              if (prtvox.eq.1) then
                y = ymin+dy*(j-1)
                ymiddle = y+dy*0.5d0
                write(out,'(1x,i5,2(1x,1pe12.5),$)') j,y,ymiddle
              endif
            endif
            if (dz.gt.0.0d0) then
              fact = fact*idz
              if (prtvox.eq.1) then
                z = zmin+dz*(k-1)
                zmiddle = z+dz*0.5d0
                write(out,'(1x,i5,2(1x,1pe12.5),$)') k,z,zmiddle
              endif
            endif
            q = edep(vox)/n
            q2 = edep2(vox)/n
            sigma = (q2-q**2)/n
            sigma = sqrt(max(sigma,0.0d0))
            q = q*fact
            sigma = sigma*fact
            write(out,'(1x,1pe12.5,1x,1pe8.1)') q,2.0d0*sigma
!           *** Evaluate average uncertainty for scores above 1/2 max score:
            if (edep(vox).gt.maxq) then
              avesig = avesig+(sigma/q)**2
              nchan = nchan+1
            endif
          enddo
        enddo
      enddo
      uncert = 100.0d0
      if (nchan.gt.0) uncert = 200.0d0*sqrt(avesig/nchan)

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


      subroutine SDDinitally(funit)
!*******************************************************************
!*    Initializes. To be called before TALLY                       *
!*******************************************************************
      implicit none
      integer funit

      integer*4 kpar,ibody,mat,ilb
      real*8 e,x,y,z,u,v,w,wght
      common/track/e,x,y,z,u,v,w,wght,kpar,ibody,mat,ilb(5)
      integer maxmat
      parameter (maxmat=10)
      integer*4 iz,nelem
      real*8 stf,zt,at,rho,vmol
      common/compos/stf(maxmat,30),zt(maxmat),at(maxmat),rho(maxmat),
     1  vmol(maxmat),iz(maxmat,30),nelem(maxmat)
      logical active
      integer prtvox
      integer*4 maxvox
      parameter(maxvox=1000000)
      integer*4 nx,ny,nz,nynz,nlast,nhist
      real*8 edptmp,edep,edep2,idens,xmin,ymin,zmin
      real*8 dx,dy,dz,idx,idy,idz
      common /scosdd/ edptmp(maxvox),edep(maxvox),edep2(maxvox),
     &                idens(maxvox),nlast(maxvox),xmin,ymin,zmin,
     &                dx,dy,dz,idx,idy,idz,nx,ny,nz,nynz,nhist,
     &                prtvox,active
      character*80 buffer
      character*(*) secid,eos
      parameter (secid=
     &'[SECTION TALLY SPATIAL DOSE DISTRIB v.2003-04-01]')
      parameter (eos='[END OF SECTION]')
      integer*4 i,j,k,vox,getabs,idpresent
      real*8 xmax,ymax,zmax,infty,dens,novoxs
      parameter (infty=1.0d30)

      active = .true.
      write(*,*) ' '
      write(*,'(a)')
     & '>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'
      read(funit,'(a80)') buffer
      idpresent = index(buffer,secid)
      if (buffer(1:1).eq.'#'.and.idpresent.ne.0) then
        active = .false.
        write(*, '(a)')
     &    '>>>> Tally Spatial Dose Distrib is INACTIVE >>>>'
        do
          read(funit,'(a80)',end=10) buffer  ! Reads dummy lines of inactive tally
          if (index(buffer,eos).ne.0) return  ! End-Of Section found
        enddo
 10     return  ! End-Of-File, something went wrong
      else if (idpresent.eq.0) then
        write(*,'(a)') 'SDDinitally:ERROR: incorrect section;'
        write(*,'(a,a)') '  expecting to find: ',secid
        write(*,'(a,a)') '  found instead:     ',buffer
        stop
      endif
      write(*,'(a)') secid

!     *** Read geometry parameters:
      write(*,'(a)') 'xmin,xmax,nx:'
      read(funit,*) xmin,xmax,nx
      write(*,'(2(1x,1pe12.5),1x,i5)') xmin,xmax,nx
      if (xmin.gt.xmax.or.nx.lt.0) then
        write(*,*) 'SDDinitally:ERROR: xmin > xmax  or  nx < 0'
        stop
      endif
      if (nx.eq.0) then
!       *** Redefine parameters when dose is integrated along x-axis:
        xmin = -0.5d0*infty
        dx = 0.0d0
        idx = 1.0d0/infty
        nx = 1
      else
        dx = (xmax-xmin)/nx
        idx = 1.0d0/dx
      endif
      write(*,'(a)') 'ymin,ymax,ny:'
      read(funit,*) ymin,ymax,ny
      write(*,'(2(1x,1pe12.5),1x,i5)') ymin,ymax,ny
      if (ymin.gt.ymax.or.ny.lt.0) then
        write(*,*) 'SDDinitally:ERROR: ymin > ymax  or  ny < 0'
        stop
      endif
      if (ny.eq.0) then
        ymin = -0.5d0*infty
        dy = 0.0d0
        idy = 1.0d0/infty
        ny = 1
      else
        dy = (ymax-ymin)/ny
        idy = 1.0d0/dy
      endif
      write(*,'(a)') 'zmin,zmax,nz:'
      read(funit,*) zmin,zmax,nz
      write(*,'(2(1x,1pe12.5),1x,i5)') zmin,zmax,nz
      if (zmin.gt.zmax.or.nz.lt.0) then
        write(*,*) 'SDDinitally:ERROR: zmin > zmax  or  nz < 0'
        stop
      endif
      if (nz.eq.0) then
        zmin = -0.5d0*infty
        dz = 0.0d0
        idz = 1.0d0/infty
        nz = 1
      else
        dz = (zmax-zmin)/nz
        idz = 1.0d0/dz
      endif
      novoxs = dble(nx)*dble(ny)*dble(nz)
      if (novoxs.gt.maxvox) then
        write(*,*)
     &    'SDDinitally:ERROR: Too many voxels; increase MAXVOX.'
        stop
      endif
      nynz = ny*nz

!     *** Print option, decides whether or not REPORT prints coordinates:
      write(*,'(a)') 'Print option:'
      read(funit,*) prtvox
      if (prtvox.eq.1) then
        write(*,'(a)') 'yes'
      else
        write(*,'(a)') 'no'
      endif

!     *** Init arrays:
      do i=1,nx
        x = xmin+dx*(i-0.5d0)
        if (dx.eq.0.0d0) x = xmax
        do j=1,ny
          y = ymin+dy*(j-0.5d0)
          if (dy.eq.0.0d0) y = ymax
          do k=1,nz
            vox = getabs(i,j,k)
            edptmp(vox) = 0.0d0
            edep(vox)   = 0.0d0
            edep2(vox)  = 0.0d0
            nlast(vox) = 0
            z = zmin+dz*(k-0.5d0)
            if (dz.eq.0.0d0) z = zmax
            u = 0.0d0
            v = 0.0d0
            w = 1.0d0
            call locate
            dens = 0.0d0
            idens(vox) = 0.0d0
            if (mat.ne.0) dens = rho(mat)
            if (dens.gt.0.0d0) idens(vox) = 1.0d0/dens
          enddo
        enddo
      enddo

      read(funit,'(a80)') buffer
      if (index(buffer,eos).eq.0) then
        write(*,*) 'SDDinitally:ERROR: End-Of-Section mark not found'
        write(*,'(a,a)') '  expecting to find: ',eos
        write(*,'(a,a)') '  found instead:     ',buffer
        stop
      endif
      write(*,'(a)') '>>>> SDD tally initialization finished >>>>'
      end


      integer*4 function getabs(i,j,k)
!*******************************************************************
!*    Gets the absolute voxel index from the cartesian ones.       *
!*                                                                 *
!*    Input:                                                       *
!*      {i,j,k} -> cartesian voxel indices                         *
!*******************************************************************
      implicit none
      integer*4 i,j,k

      logical active
      integer prtvox
      integer*4 maxvox
      parameter(maxvox=1000000)
      integer*4 nx,ny,nz,nynz,nlast,nhist
      real*8 edptmp,edep,edep2,idens,xmin,ymin,zmin
      real*8 dx,dy,dz,idx,idy,idz
      common /scosdd/ edptmp(maxvox),edep(maxvox),edep2(maxvox),
     &                idens(maxvox),nlast(maxvox),xmin,ymin,zmin,
     &                dx,dy,dz,idx,idy,idz,nx,ny,nz,nynz,nhist,
     &                prtvox,active

      getabs = k+(j-1)*nz+(i-1)*nynz
      end


! **** end of file *************************************************
