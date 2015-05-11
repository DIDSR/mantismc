!*******************************************************************
!*                          TALLY                                  *
!*                         DETECT2                                 *
!*       (generation and transport of optical photons)             *
!*                                                                 *
!* Short description:                                              *
!*   These routines link penEasy to DETECT-2. The latter simulates *
!*   the transport of optical photons in matter.                   *
!*                                                                 *
!* Dependencies:                                                   *
!*   from PENELOPE:                                                *
!*   -> common /TRACK/                                             *
!*   from DETECT2                                                  *
!*   -> routines DETECT2,DETECT2EOH                                *
!*                                                                 *
!* Compatible with PENELOPE versions:                              *
!*   2005                                                          *
!*                                                                 *
!* Comments:                                                       *
!*   -> Not adapted to variance reduction; do not use with WGHT!=1.*
!*   -> This tally does not require the init and report routines;  *
!*      they are defined inside the DETECT2 code.                  *
!*                                                                 *
!* Last update:                                                    *
!*   2005-08-01                                                    *
!*     Created.                                                    *
!*                                                                 *
!* Josep Sempau, email: josep.sempau@upc.es                        *
!*   Universitat Politecnica de Catalunya, Barcelona, Spain        *
!* Aldo Badano, email: agb@cdrh.fda.gov                            *
!*   Food and Drug Admin, Rockville, MD USA                        *
!*******************************************************************


      subroutine detect2tally(mode,arg)
!*******************************************************************
!*    Input:                                                       *
!*      mode -> Identifies the state of the calling procedure      *
!*      arg -> energy loss (modified on output) or hist no.        *
!*******************************************************************
      implicit none
      integer mode
      real*8 arg

      integer*4 kpar,ibody,mat,ilb
      real*8 e,x,y,z,u,v,w,wght
      common/track/e,x,y,z,u,v,w,wght,kpar,ibody,mat,ilb(5)
      ! Secondary stack:
      real*8 es,xs,ys,zs,us,vs,ws,wghts
      integer*4 ks,ibodys,ms,ilbs,nsec,nms
      parameter (nms=1000)
      common/secst/es(nms),xs(nms),ys(nms),zs(nms),us(nms),
     1   vs(nms),ws(nms),wghts(nms),ks(nms),ibodys(nms),ms(nms),
     2   ilbs(5,nms),nsec
      integer nsecold
      integer*4 nhist
      common /scodt2/ nhist,nsecold
      integer crystalMat,i,left
      parameter (crystalMat=1)
      real*8 de

      if (mode.lt.0.and.mat.eq.crystalMat) then
        if (nsec.ge.nms) then
          write(*,*) ' '
          write(*,*)
     &      'main:WARNING: stack overflow; try to increase'
          write(*,*)
     &      '  absorption energies or parameter NMS.'
          write(*,*)
     &      '  The current history has been discarded.'
          call cleans  ! Clear the secondary stack
          mat = 0      ! Force main to stop this particle's simulation
          return
        endif
        de = arg
        do i=nsecold+1,nsec    ! For all newly created particles
          de = de-es(i)        ! Edep=W-Estack
        enddo
        if (de.gt.0.0d0) then  ! Protects from ISI-like events
          call stores(e,x,y,z,u,v,w,wght,kpar,ilb)  ! Put particle in stack
          e = de  ! Informs DETECT2 about the amount of energy locally deposited
          call detect2(nhist)
          call secpar(left)    ! Recover particle previously stored
        endif

      else if (mode.eq.3.or.mode.eq.4.or.mode.eq.-99) then
        ! No. of secs before a KNOCK, change of MAT or a new particle happens
        nsecold = nsec

      else if (mode.eq.6) then
        call detect2EOH  ! EOH tasks

      else if (mode.eq.1) then
        nhist = arg+0.1d0  ! Safe real-to-integer conversion

      endif
      end


! **** end of file *************************************************

