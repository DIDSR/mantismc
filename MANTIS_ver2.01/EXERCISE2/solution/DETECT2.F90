!-------------------------------------------------------------------!
!                            DETECT-II
!       (subroutines for MANTIS: PENELOPE/DETECT-II combination)
!
! Aldo Badano (aldo.badano@fda.hhs.gov) and Josep Sempau
! Original development by: Aldo Badano, Michael J. Flynn
!
! current  version dII_vFDA1.0                        2002
! beta     version dII_vFDA2.0                        2003
! fda2     -
! fda3     addition of detector reflectivity
! fda4     addition of normal deviate subroutine for gain variance
! fda5     serial adjustment in beowolf and parallelization starts
! fda6     new thin-film description (thanks to Shu-Jen Lee)
! fda6p    for wrapper using penelope
! fda6pen  wrapper for penelope (Josep Sempau's visit to FDA)  2004
! DETECT2.F90 subroutine for MANTIS                            2004
! detect2.f90 subroutine for MANTIS (using penelope2005)
!            (Josep Sempau's 2nd visit to FDA)                 2005
! detect2.f90 parallelizable random number generator       Sep/2005 
! detect2.f90 Intel intrinsic MT random number generator   Jan/2006 
! detect2.f90 parallelizable MT random number generator    Jan/2007 
!--------------------------------------------------------------------!
! modules and type declarations to run MT parallel number generation
    ! module mtpar
    ! use mkl_vsl_type
    ! use mkl_vsl
    ! type (vsl_stream_state) :: stream
! 
    ! contains
! 
    ! function irmt(seed,seq)
        ! integer :: brng,seq
        ! integer seed
        ! brng=vsl_brng_mt2203+seq
        ! status = vslnewstream( stream, brng,  seed )
    ! end function irmt 
! 
    ! function rmt(r)
        ! integer :: status,n
        ! real*8 :: a, b, o(1), r
        ! integer :: method, brng
        ! a=0.D0
        ! b=1.D0
        ! method=0
        ! n=1
        ! status = vdrnguniform( method, stream, n, o, a, b)
        ! r=o(1)
    ! end function rmt
! 
    ! end module mtpar
!--------------------------------------------------------------------!
    module global_variables

    integer :: status
    integer, parameter :: fine=selected_real_kind(9,9)
    integer, parameter :: long=selected_int_kind(9)

! global_variables new in penelope version detect2
    integer :: lost,emerged                     ! checking lost photons
    integer :: matn                             ! # materials
    integer :: pmean                            ! mean opticals/keV
    integer :: ehgen                               ! # eh-pairs per nm
    logical :: countf                           ! signal or count
    real*8  :: imgdx,imgdy                   ! image pixel size
    integer*4 dii_ncross
    real*8 dii_dseff
    integer, allocatable, dimension(:,:) :: surftype ! surf defs
    integer :: vtype                                ! type of variance calc
    integer :: vstep                                ! frequency of variance calc
    integer :: model                                ! model of gain variations
    integer :: emodel                                ! model of scint. energy dep
    real*8 :: verror                        ! variance (error)
    real*8 :: mom1,mom2,mom0,sf,sf_old
    real*8 :: errorsf,totalsf,squarredsf        ! moments and Swank factor
    integer, allocatable, dimension(:) :: dethistvar
    logical :: reportplot                       ! flag for saving data
    integer*4 xm,xmp,xmp_old                    ! counters for histories per dii
    character infile*100                        ! code for jobname

! global_variables contains variable, namelist, and function definitions
    logical :: flagdie                           ! flag for perfect absorbers
    logical :: flaga                             ! flag for all-use
    logical :: inside_flag,flag_pp               ! flag for pressed powder model
    real*8 ::  pmax_mu,pmax_si2
    integer*4, dimension(:), allocatable :: produced ! produced photons
    integer*4 :: absorbed               ! count of absorbed photons
    integer*4 :: scatters               ! count of scatterings
    integer*4:: p,pzero ! loop counter for histories
    integer*4 :: c      ! loop counter for events
    integer :: cmax              ! maximum number of events allowed
    integer :: ab                ! flag for collision type
    integer :: matnumber         ! # of mats in this run
    integer, dimension(:), allocatable :: wpoints ! points w-dept properties
    integer, dimension(:,:), allocatable :: wavepoints ! w-dept properties
    integer :: k,l,i,j,jj,n,m,iw   ! general loop counters
    integer :: bins,surftemp1,surftemp
    integer :: nzcell            ! number of z-slabs
    integer :: nofre            ! flag for no-fresnel analysis
    integer :: ptype            ! flag for source shapes
    integer :: scat_type        ! flag for scattering types
    integer :: icell            ! cell counter
    integer :: expslab            ! slab for exponential source
    integer :: ptemp,temp,ttemp        ! temp
    integer :: bwave,bsensorwave    ! # of bins for source spectrum
    integer :: iplane            ! flag for surface planar sources
    integer :: sg            ! flag for source geometry
    integer :: sd            ! flag for source directionallity
    integer :: se            ! flag for source wavelength
    integer :: sp            ! flag for source polarization
    integer :: periodic_elements    ! periodic number
    integer :: surfdefault        ! surface type default
    integer*4, dimension(:), allocatable :: detec
    integer*4 :: pmax ! count of photons per x-ray interaction
    integer*4 :: counted    ! count of counted photons
    integer*4 :: sumscat    ! count of scatterings
    real*8 :: detected        ! count of detected photons
    integer, dimension(:), allocatable :: periodic_mat
    integer*4 :: died    ! count of photons died in perfect absorber
    integer :: mat1, mat2    ! flags defining material of cell
    integer :: bt        ! flag for display application binning type
    integer :: be        ! flag for scattering event statistics
    integer :: bw        ! flag for wavelength distribution binning
    integer :: ba        ! flag for angular distribution binning type
    integer :: bz        ! flag for z distribution binning type
    integer :: b3d        ! flag for 3d paths
    integer :: seed        ! initial seed for rngs
    integer :: e1,e2        ! range for wavelength interpolation
    integer :: rint        ! integer random number [0:999]
    integer :: xsize,ysize    ! binning array size
    integer, dimension(:,:), allocatable :: d3  ! psf for 3D plots (x,y,z)
    integer, dimension(:), allocatable :: detlsf  ! output psf for detector runs
    integer, dimension(:,:), allocatable :: d3d  ! 3D plots for detector planes
    integer, dimension(:,:,:), allocatable :: matpen   ! material flag
    integer, dimension(:), allocatable :: pwave  ! exit wavelength distribution
    integer, dimension(:), allocatable :: binrad      ! radius binning array
    integer, dimension(:), allocatable :: scatev   ! # of scatterings variable
    integer, dimension(:), allocatable  :: anglepen ! angular binning array
    integer, dimension(:), allocatable :: nxcell,nycell   ! # of x,y cells
    integer :: fileunit,fileunit2,fileunit3     ! finding unit numbers (from penelope)

    real*8 :: dcoll        ! distance to collision
    real*8 :: dcoll_pp    ! distance to grain boundary in pp model
    real*8 :: dbound    ! distance to boundary
    real*8 :: t,ttest
    real*8 :: maxtemp,maxpec
    real*8 :: beta        ! angle of restriction in roughness
    real*8 :: iangle        ! angle for specular reflectance calculations
    real*8 :: radius        ! radius of exiting photon
    real*8 :: rado        ! radius for circular source
    real*8 :: atest        ! flag for mirror preformance
    real*8 :: expl        ! exponential absorption for expl sources
    real*8 :: r,rr           ! random number
    real*8 :: rrr            ! auxilliary random number
    real*8 :: tthick
    real*8 :: radmaxtemp        ! radmax temporary for loop usage
    real*8 :: maxall          ! max possible distance in object
    real*8 :: para        ! parallel component
    real*8 :: posx,posy        ! storage variables for position
    real*8 :: phimax,phimaxr
    real*8 :: time1,time2     ! time tracking
    real*8 :: grain_size        ! grain size in microns for pp model
    real*8 :: fresnel_pa,fresnel_pe    ! reflection coefficients
    real*8 :: ain,aout              ! incidence/transmitted angles and cos(sin)
    real*8 :: cosain,sinaout  ! incidence/transmitted angles and cos(sin)
    real*8 :: r1, r2      ! radii for directional emission binning
    real*8 :: dd          ! inside radius of test pattern
    real*8 :: ddo         ! outside radius of test pattern
    real*8 :: cone        ! angle for restricted sourcing
    real*8 :: absdiffrac     ! absorbed fraction for absdifuser
    real*8 :: absroughfrac     ! absorbed fraction for rough
    real*8 :: abssmoothfrac    ! absorbed fraction for rough
    real*8 :: detrefl        ! reflectivity of detector plane
    real*8 :: absthinfrac        ! absorbed fraction for thin film
    real*8 :: absmirfrac        ! absorbed fraction for mirror
    real*8 :: lsfraction         ! fraction of lamb vs mirror reflector
    real*8 :: phipen,theta   ! coordinate system angles
    real*8 :: sinphi      ! sin(phi)
    real*8 :: fresnel     ! probability of Fresnel reflection
    real*8 :: meandist    ! averaged distance
    real*8 :: radmax      ! maximum radius for binning
    real*8 :: psensor_actual
    real*8, dimension(:), allocatable :: periodic
    real*8, dimension(:), allocatable :: wspect,wsensor ! w/source spectrum
    real*8, dimension(:), allocatable :: pspect,pspec,psensor ! spectrum
    real*8, dimension(:), allocatable :: fabs   ! absorption fractions
    real*8, dimension(:), allocatable :: thickpen ! z-slabs thickness
    real*8, dimension(3) :: ipos ! launch site or plane
    real*8, dimension(4) :: ilimits ! limit planar sources in 2 directions
    real*8, dimension(3) :: maxpos ! max position x,y,z (cm)
    real*8, dimension(3) :: pb   ! pencil beam location in xy plane
    real*8, dimension(3) :: normal,normal_old ! normal to surfaces
    real*8, dimension(:,:), allocatable :: dycell,dxcell ! dims of x/y cells
    real*8, dimension(:,:), allocatable :: pmat ! material properties
    real*8, dimension(:,:,:), allocatable :: edprop ! w-dependent props
    real*4, parameter :: pipen=3.141592654, twopipen=6.283185308
    real*8, parameter :: adjust=1.0D-8

!!!!! MODIFIED FOR DEPTH OF INTERACTION    
real*8  :: depth_x,depth_y,depth       
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     
! this variables introduced by the use of James/Marsaglia 1990 rng
    real*8, dimension(97) :: uu
    real*8 :: ran_c,ran_cd,ran_cm,ran_s,ran_t,uni
    integer :: i97,j97,ij,kl,ran_i,ran_j,ran_k,ran_l,ran_m,ran_ii,ran_jj,ivec

! this variables introduced by thin film routine
    logical :: tf_flag,sens_flag
    logical, dimension(:), allocatable :: emissive_flag
    integer :: tf_number
    real*8, dimension(:), allocatable :: tf_mat,tf_th

! type definition
    type photon
        integer :: events                     ! # events along track
        integer :: slab,x,y                   ! slab and cell #
        real*8, dimension(3) :: pos       ! position x,y,z (cm)
        real*8, dimension(3) :: pos_old   ! old position
        real*8, dimension(3) :: pos_tf    ! tf position
        real*8, dimension(3) :: dcos      ! directional cosines
        real*8, dimension(3) :: dcos_old  ! old dcos
        real*8, dimension(3) :: dcos_tf   ! tf dcos
        real*8, dimension(3) :: polar     ! polarization
        real*8 :: polar_perp        ! perp component
        real*8 :: wavel             ! wavelenght (nm)
        real*8 :: wavel_old        ! old wavelength
    end type photon
    type(photon) pho

! namelist definitions
    namelist/mcparameters/ pmax_mu,cmax,model,emodel
    namelist/sourcetype/ sg,sd,se,sp,nofre
    namelist/scattering_type/ scat_type
    namelist/binangle/ phimax,beta,cone
    namelist/pressed_powder/ grain_size
    namelist/surftypes/ abssmoothfrac,absdiffrac,absroughfrac,absmirfrac, &
                        & absthinfrac,lsfraction,detrefl
    namelist/startposition/ ipos
    namelist/planarposition/ iplane,ipos
    namelist/planarposition_type/ ptype
    namelist/planarposition_limits/ ilimits
    namelist/surf/ surfdefault
    namelist/tf/ tf_number
    namelist/pencilbeam/ pb
    namelist/incidence_angle/ iangle
    namelist/binwave/ bwave
    namelist/binsensorwave/ bsensorwave
    namelist/binningtype/ reportplot
    namelist/binarraysize/ xsize,ysize,imgdx,imgdy
    namelist/variance/ vtype,verror,vstep

! functions definitions
    contains
        function nor(v1,v2)
            real*8,dimension(3) :: v1,v2,nor
            nor(1)=v1(2)*v2(3)-v2(2)*v1(3)
            nor(2)=v1(3)*v2(1)-v2(3)*v1(1)
            nor(3)=v1(1)*v2(2)-v2(1)*v1(2)
        end function nor
        function normalize(v)
            real*8 :: length
            real*8,dimension(3) :: v,normalize
            length=sqrt(v(1)*v(1)+v(2)*v(2)+v(3)*v(3))
            normalize=v/length
        end function normalize
        function gainfactor_csi_1(vv)
            real*8 :: vv,v
			real*8 :: gainfactor_csi_1
            v=vv*1.e-3
            if (v .gt. 60) then
                gainfactor_csi_1=1.
            else
                gainfactor_csi_1= 1.D-11*v**6-5.D-9*v**5+7.D-7*v**4  &
                & -5.D-5*v**3+0.0019D0*v**2-0.0361D0*v+1.2546D0
            end if
        end function gainfactor_csi_1
        function gainfactor_csi_2(vv)
            real*8 :: vv,v
			real*8 :: gainfactor_csi_2
            v=vv*1.e-3
            if (v .gt. 60) then
                gainfactor_csi_2=1.
            else
                gainfactor_csi_2= -1.D-12*v**6+2.D-10*v**5-3.D-6*v**3  &
                & -.0003D0*v**2+0.0145D0*v+0.7234D0
            end if
        end function gainfactor_csi_2

    end module global_variables
!----------------------------------------------------------------------!
! The twister mersenne prng from 
! http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/mt.html
! Fortran translation by Hiroshi Takano.  Jan. 13, 1999.
!   genrand()      -&gt; double precision function grnd()
!   sgenrand(seed) -&gt; subroutine sgrnd(seed)
! Fortran version rewritten as an F90 module and mt state saving and getting
! subroutines added by Richard Woloshyn. (rwww@triumf.ca). June 30, 1999
 module mtpar
! Default seed
   integer, parameter :: defaultsd = 4357
! Period parameters
   integer, parameter :: mtN = 624, mtN1 = mtN + 1
! the array for the state vector
   integer, save, dimension(0:mtN-1) :: mt
   integer, save                   :: mti = mtN1
!
    contains
!
!Initialization subroutine
   subroutine sgrnd(seed)
     implicit none
!      setting initial seeds to mt[N] using
!      the generator Line 25 of Table 1 in
!      [KNUTH 1981, The Art of Computer Programming
!         Vol. 2 (2nd Ed.), pp102]
     integer, intent(in) :: seed
     mt(0) = iand(seed,-1)
     do mti=1,mtN-1
       mt(mti) = iand(69069 * mt(mti-1),-1)
     enddo
     return
   end subroutine sgrnd
!
!!Random number generator
    !!real(8) function grnd()
    function rmt(r)
      implicit integer(a-z)
      real*8 :: r
!! Period parameters
      integer, parameter :: mtM = 397, MATA  = -1727483681
!!                                    constant vector a
      integer, parameter :: LMASK =  2147483647
!!                                    least significant r bits
      integer, parameter :: UMASK = -LMASK - 1
!!                                    most significant w-r bits
!! Tempering parameters
      integer, parameter :: TMASKB= -1658038656, TMASKC= -272236544
      dimension mag01(0:1)
      data mag01/0, MATA/
      save mag01
!!                        mag01(x) = x * MATA for x=0,1
      TSHFTU(y)=ishft(y,-11)
      TSHFTS(y)=ishft(y,7)
      TSHFTT(y)=ishft(y,15)
      TSHFTL(y)=ishft(y,-18)
      if(mti.ge.mtN) then
!!                       generate N words at one time
        if(mti.eq.mtN+1) then
!!                            if sgrnd() has not been called,
          call sgrnd( defaultsd )
!!                              a default initial seed is used
        endif
!
        do kk=0,mtN-mtM-1
            y=ior(iand(mt(kk),UMASK),iand(mt(kk+1),LMASK))
            mt(kk)=ieor(ieor(mt(kk+mtM),ishft(y,-1)),mag01(iand(y,1)))
        enddo
        do kk=mtN-mtM,mtN-2
            y=ior(iand(mt(kk),UMASK),iand(mt(kk+1),LMASK))
            mt(kk)=ieor(ieor(mt(kk+(mtM-mtN)),ishft(y,-1)),mag01(iand(y,1)))
        enddo
        y=ior(iand(mt(mtN-1),UMASK),iand(mt(0),LMASK))
        mt(mtN-1)=ieor(ieor(mt(mtM-1),ishft(y,-1)),mag01(iand(y,1)))
        mti = 0
      endif

      y=mt(mti)
      mti = mti + 1
      y=ieor(y,TSHFTU(y))
      y=ieor(y,iand(TSHFTS(y),TMASKB))
      y=ieor(y,iand(TSHFTT(y),TMASKC))
      y=ieor(y,TSHFTL(y))
  
      if(y .lt. 0) then
        r=(dble(y)+2.0d0**32)/(2.0d0**32-1.0d0)
      else
        r=dble(y)/(2.0d0**32-1.0d0)
      endif
      return
      end function rmt 
!
    end module mtpar
!----------------------------------------------------------------------!
    module nr

        IMPLICIT NONE
        INTEGER, PARAMETER :: NPAR_ARTH=16,NPAR2_ARTH=8

        CONTAINS

        FUNCTION arth(first,increment,n)
            REAL, INTENT(IN) :: first,increment
        INTEGER, INTENT(IN) :: n
            REAL, DIMENSION(n) :: arth
        INTEGER :: k,k2
            REAL :: temp
        if (n > 0) arth(1)=first
        if (n <= NPAR_ARTH) then
        do k=2,n
        arth(k)=arth(k-1)+increment
        end do
        else
        do k=2,NPAR2_ARTH
        arth(k)=arth(k-1)+increment
        end do
        temp=increment*NPAR2_ARTH
        k=NPAR2_ARTH
        do
        if (k >= n) exit
        k2=k+k
        arth(k+1:min(k2,n))=temp+arth(1:min(k,n-k))
        temp=temp+temp
        k=k2
        end do
        end if
        END FUNCTION arth

        FUNCTION gammln(xx)
        IMPLICIT NONE
            REAL, INTENT(IN) :: xx
            REAL :: gammln
            REAL :: tmp,x
            REAL :: stp = 2.5066282746310005
            REAL, DIMENSION(6) :: coef = (/76.18009172947146,&
        -86.50532032941677,24.01409824083091,&
        -1.231739572450155,0.1208650973866179e-2,&
        -0.5395239384953e-5/)
        x=xx
        tmp=x+5.5
        tmp=(x+0.5)*log(tmp)-tmp
        gammln=tmp+log(stp*(1.000000000190015+ &
        & sum(coef(:)/arth(x+1.0,1.0,size(coef))))/x)

        END FUNCTION gammln

    end module nr
!-----------------------------------------------------------------------
    module poisson

        contains

        function poidev(xpm)
        use nr
        use global_variables
        use mtpar
        IMPLICIT NONE
        REAL, INTENT(IN) :: xpm
        REAL :: poidev
        REAL :: em,tpoi,ypoi
        REAL, SAVE :: alxm,g,oldm=-1.0,sq

        if (xpm < 12.0) then
        if (xpm /= oldm) then
        oldm=xpm
        g=exp(-xpm)
        end if
        em=-1
        tpoi=1.0
        do
        em=em+1.0
        status=rmt(r)
        tpoi=tpoi*real(r)
        if (tpoi <= g) exit
        end do
        else
        if (xpm /= oldm) then
            oldm=xpm
            sq=sqrt(2.0*abs(xpm))
            alxm=log(xpm)
            g=xpm*alxm-gammln(xpm+1.0)
        end if
        do
        do
        status=rmt(r)
        ypoi=tan(pipen*real(r))
        em=sq*ypoi+xpm
        if (em >= 0.0) exit
        end do
        em=int(em)
        tpoi=0.9*(1.0+ypoi**2)*exp(em*alxm-gammln(em+1.0)-g)
        status=rmt(r)
        if (r <= tpoi) exit
        end do
        end if

        poidev=em
        END FUNCTION poidev

    end module poisson
!----------------------------------------------------------------------!
subroutine detect2(xmpl)
! DETECT-II subroutine for PENELOPE runs
! - xmpl refers to the history according to penelope

    use global_variables
    use poisson
    use mtpar
    implicit none
    external dii_source,distcoll,scatter,binevscat,collana
    external boundana,tallyd2

    integer*4 kpar,ibody,mat,ilb
    integer*4 ibody_old,mat_old
    real*8 e,x,y,z,u,v,w,wght
    common/track/e,x,y,z,u,v,w,wght,kpar,ibody,mat,ilb(5)

    integer*4, INTENT(IN) ::  xmpl

! count histories according to penelope and detect2
    xmp=xmpl
    if (xmp .gt. xmp_old) then
        xm=xm+1         ! xm is the histories in phosphor
        xmp_old=xmp
    end if

! energy-dependence of the gain not used in version MANTIS 1.0
! model for variance due to energy dependence of the gain
!    select case (emodel)
!    case (0)
         pmean=nint(pmax_mu*e)
!    case (1)
!        pmean=nint(pmax_mu*e * gainfactor_csi_1(e))
!    case (2)
!        pmean=nint(pmax_mu*e * gainfactor_csi_2(e))
!    end select
!    if (pmean .le. 0) return

! record deposition of energy into file TEMPORARY
    if (reportplot .eqv. .true.) then
       if (xm >= 1 .and. xm <= 5) then
          write(100+xm,'(4(1x,es15.5),1x,i1,1x,i3,1x,i3,es12.2)') &
             & x*1.e4,y*1.e4,z*1.e4,e,mat,xm,pmean
       else if (xm > 5) then
          close(101); close(102)
          close(103); close(104); close(105)
       end if
    end if

    select case (model)
    case (0)
       pmax=pmean
    case (1)
       call gasdev_s(r)
        pmax=pmean+nint(r*pmean)
    case (2)
       pmax=nint(poidev(real(pmean)))
    end select
    if (pmax .le. 0.) return
    produced(xm)=produced(xm)+pmax
    ipos(1)=x
    ipos(2)=y
    ipos(3)=z
    ibody_old=ibody
    mat_old=mat

! start histories; pmax from normal (mu=pmax_mu and si^2=pmax_si2)
    histories: do p=1,pmax
       call dii_source
       x=pho%pos(1)
       y=pho%pos(2)
       z=pho%pos(3)

!!!!!! MODIFICATION FOR DEPTH OF INTERACTION
       x=depth_x;y=depth_y;z=depth   
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

       u=pho%dcos(1)
       v=pho%dcos(2)
       w=pho%dcos(3)
       ibody=ibody_old
       mat=mat_old

! loop over all events
       event: do c=1,cmax

! counting of photons lost
             if (c==cmax) then
                lost=lost+1
                if (be==2) call binevscat
                pho%events=pho%events+1
                exit event
             end if

! record optical paths into file (TEMPORARY)
             if (reportplot .eqv. .true.) then
                if (xm>= 1.and.xm<=5.and.c.lt.20.and.p==1) then
                   write(110+xm,'(3(1x,es14.4),1x,i3,1x,i2,i3,1x,i3)') &
                   & x*1.e4,y*1.e4,z*1.e4,mat,p,c
                else if (xm > 5) then
                   close(111);close(112);close(113);close(114);close(115)
                end if
             end if

! counting photons into diode material (mat #2)
             if (mat==2) then
                   call tallyd2
                   exit event
             end if

! move optical photon
             call distcoll
             mat1=mat
             call step(dcoll,dii_dseff,dii_ncross)

! counting photons emerging out thru vacuum
            if (mat == 0 .or. mat==3) then
                emerged=emerged+1
                exit event
             end if

             mat2=mat
             meandist=meandist+dii_dseff

! test for collision in cell
             collision: if (dii_ncross == 0) then
                 meandist=meandist+dcoll
                 ab=0
                 call collana
                 select case (ab)
                    case (1)
                        absorbed=absorbed+1
                        if (be==2) call binevscat
                        exit event
                     case (2)
                        scatters=scatters+1
                        if (scat_type .eq. 1) then
                           call scatter
                        else if (scat_type .eq. 2) then
                           call isoscatter
                        end if
                        u=pho%dcos(1)
                        v=pho%dcos(2)
                        w=pho%dcos(3)
                        if (be==2) call binevscat
                     end select

! if no collision, analyze boundary crossing
              else 

! boundary crossing analisys
                 flagdie=.false.
                 call snormal(normal(1),normal(2),normal(3))
                 call boundana
! counting of photons that died
                 if (flagdie .eqv. .true.) then
                     died=died+1
                     flagdie=.false.
                     if (be==2) call binevscat
                     exit event
                 end if
              end if collision
           end do event
    end do histories
    return

end subroutine detect2
!------------------------------------------------------------
subroutine detect2EOH
! update counters for variance calculation

    use global_variables
    use mtpar
    implicit none

end subroutine detect2EOH
!------------------------------------------------------------
subroutine boundana

! subroutine boundana returns new photon after a boundary intersection
    use global_variables
    use mtpar
    implicit none

    integer*4 kpar,ibody,mat,ilb
    real*8 e,x,y,z,u,v,w,wght
    common/track/e,x,y,z,u,v,w,wght,kpar,ibody,mat,ilb(5)


! surfaces defined from 2 materials
    surftemp=surftype(mat1,mat2)

! check surface type with flag
    select case (surftemp)
    case(1)
        call smooth
    case(11)
        call smoothabs
    case(2)
        call rough
        call smooth
    case(3)
        call rough
        call mirror
    case(4)
        call mirror
        flaga=.true.
    case(5)
        call killer
    case(6)
        call absdifuser
        flaga=.true.
    case(8)
        call transparent
    case(9)
    ! determine if reflection occurs with random number
         status=rmt(r)
         if (r<=detrefl) then
             call difuser
             flaga=.true.
         else
        ! flaga=.false.     ! indicates refraction but no change in direction
        ! flagdie=.false.
             call smooth
         end if
    case(25)
        status=rmt(r)
        if (r > lsfraction) then
            call absdifuser
        else
            call mirror
        end if
        flaga=.true.
    end select
    pho%dcos=normalize(pho%dcos)
    u=pho%dcos(1)
    v=pho%dcos(2)
    w=pho%dcos(3)

! if reflection, un-do additional displacement into next material
! by moving the photon with xstep into "old" material
    if (flaga .eqv. .true.) then
        call step(1.D30,dii_dseff,dii_ncross)
        pho%pos(1)=x
        pho%pos(2)=y
        pho%pos(3)=z
    end if
    return

end subroutine boundana
!---------------------------------------------------------------------
subroutine rough

! subroutine rough model a surface whose normal is randomly tilted
! makes use of the smooth routine
    use global_variables
    use mtpar
    implicit none
    real*8, dimension(3) :: normal_pert

! procedure is adding small vector to normal and check for angle
! n'=n+beta*pert where pert is isotropic and beta is constant (10%)
    status=rmt(r)
    normal_pert(3) = 2.*r-1.
    rr=sqrt(1.-r*r)
    status=rmt(r)
    theta=r*twopipen
    normal_pert(1)=rr*cos(theta)
    normal_pert(2)=rr*sin(theta)
    normal=normalize(normal+beta*normal_pert)
    return

end subroutine rough
!---------------------------------------------------------------------
subroutine smooth

! for smooth (flat) surfaces
    use global_variables
    use mtpar
    implicit none
    real*8 :: cosdef,cosaout
    external locate,fresnels

! when indexes are equal, return same dcos and same polarization
    if (pmat(mat1,1)==pmat(mat2,1)) then
        flaga=.false.
        return
    end if

! store current normal to plane of incidence and dcos
    normal_old=normal
    pho%dcos_old=pho%dcos

! compute cosine of incidence cosain (normal is surface normal)
    cosain=dot_product(normal,pho%dcos)
    if (cosain .lt. -1.0) cosain=-1.0
    if (cosain .gt. 1.0) cosain=1.0
    cosdef=max(0.D0,1.D0-cosain**2)

! compute reflected vector
    pho%dcos=2.*dot_product(-pho%dcos,normal)*normal + pho%dcos  ! specular ray

! case of total internal reflection (no change in polarization)
! initialize flaga and then count for cases:
! flaga= .F.refraction .true.reflection

    flaga=.false.
    if (abs(sqrt(cosdef)*(pmat(mat1,1)/pmat(mat2,1)))>1.) then
        flaga=.true.
        return
    else
        sinaout=pmat(mat1,1)/pmat(mat2,1)*sqrt(cosdef)
        cosaout=sqrt(abs(1.0-sinaout**2))

! compute probability of reflection using Fresnel's law
        call fresnels

! if reflection, no action needed; if refraction, update dcos
        if (flaga .eqv. .true.) then
               return
        else
               pho%dcos = pmat(mat1,1)/pmat(mat2,1)*pho%dcos_old -  &
               & (cosaout -  pmat(mat1,1)/pmat(mat2,1) * cosain)*normal
               pho%dcos=normalize(pho%dcos)

! compute new perpendicular component using probability of
! transmission for perpendicular is 1-reflection given by fresnel_pe
! if reflection occurs, no change in perpendicular component pho%polar_perp
               if (sp/=1) then
                      if (pho%polar_perp < adjust) pho%polar_perp = adjust
                             pho%polar_perp=pho%polar_perp*sqrt(1.-fresnel_pe)
                             where (normal .eq. 1.) pho%polar=pho%polar_perp
                             pho%polar=normalize(pho%polar)
                      end if
        end if
     end if

    return

end subroutine smooth
!---------------------------------------------------------------------
subroutine smoothabs

! for smooth (flat) surfaces with absorption
    use global_variables
    use mtpar
    implicit none
    status=rmt(r)
    if (r<abssmoothfrac) then
        flagdie=.true.
        return
    end if
    call smooth

end subroutine smoothabs
!---------------------------------------------------------------------
subroutine transparent
! for surfaces that act as no change occured
    use global_variables
    implicit none
    flaga=.false.     ! indicates refraction but no change in direction
    flagdie=.false.
    return

end subroutine transparent
!---------------------------------------------------------------------
subroutine killer
! for surfaces that act as perfect absorbers
    use global_variables
    implicit none
    flagdie=.true.
    return

end subroutine killer
!----------------------------------------------------------------------
subroutine absdifuser

! for surfaces that act as partial absorbers with an absorbed fraction
! of absdiffrac and randomly lambertian reflection of survivors

    use global_variables
    use mtpar
    implicit none
    real*8, dimension(3) :: dcos_temp,dcos_pa,dcos_nor,cross
    real*8 :: dot

! determine if absorption occurs with random number
    status=rmt(r)
    if (r<absdiffrac) then
        flagdie=.true.
        return
    end if

! reflect photon with lambertian random directional cosine
    status=rmt(r)
    phipen=acos(1.-1.99999*r)*.5
!   dcos_temp(1)=-cos(phipen)
    status=rmt(r)
    dcos_temp(1)=cos(phipen) * sign(1.0d0,r-0.5d0)
    status=rmt(r)
    theta=r*twopipen
    sinphi=sqrt(1.D0-dcos_temp(1)**2)
    status=rmt(r)
    dcos_temp(2)=sinphi*cos(theta) * sign(1.0d0,r-0.5d0)
    dcos_temp(3)=abs(sinphi*sin(theta)) * sign(normal(3),-1.0d0)

    pho%dcos(1)=dcos_temp(1)
    pho%dcos(2)=dcos_temp(2)
    pho%dcos(3)=dcos_temp(3)

! ROTATION along an axis (i.e., the normal)
! compute the dot product of the vector and the rotation axis.
    dot=dot_product(dcos_temp,normal)

! compute the parallel component of the vector.
    dcos_pa= dot * normal

! compute the normal component of the vector.
    dcos_nor(1)= dcos_temp(1)-dcos_pa(1)
    dcos_nor(2)= dcos_temp(2)-dcos_pa(2)
    dcos_nor(3)= dcos_temp(3)-dcos_pa(3)
    dcos_nor=normalize(dcos_nor)

! compute a second vector, lying in the plane, perpendicular
    cross(1)=normal(2)*dcos_nor(3)-normal(3)*dcos_nor(2)
    cross(2)=normal(3)*dcos_nor(1)-normal(1)*dcos_nor(3)
    cross(3)=normal(1)*dcos_nor(2)-normal(2)*dcos_nor(1)
    cross=normalize(cross)

! rotated vector is the parallel component plus the rotated component.
    pho%dcos(1)=dcos_pa(1) + cos(dot)* dcos_nor(1) + sin(dot)*cross(1)
    pho%dcos(2)=dcos_pa(2) + cos(dot)* dcos_nor(2) + sin(dot)*cross(2)
    pho%dcos(3)=dcos_pa(3) + cos(dot)* dcos_nor(3) + sin(dot)*cross(3)
    pho%dcos=normalize(pho%dcos)

! END ROTATION
  return

end subroutine absdifuser
!---------------------------------------------------------------------
subroutine difuser

! for randomly lambertian reflection

    use global_variables
    use mtpar
    implicit none
    real*8, dimension(3) :: dcos_temp,dcos_pa,dcos_nor,cross
    real*8 :: dot

! reflect photon with lambertian random directional cosine
    status=rmt(r)
    phipen=acos(1.-1.99999*r)*.5
!   dcos_temp(1)=-cos(phipen)
    status=rmt(r)
    dcos_temp(1)=cos(phipen) * sign(1.0d0,r-0.5d0)
    status=rmt(r)
    theta=r*twopipen
    sinphi=sqrt(1.D0-dcos_temp(1)**2)
    status=rmt(r)
    dcos_temp(2)=sinphi*cos(theta) * sign(1.0d0,r-0.5d0)
    dcos_temp(3)=abs(sinphi*sin(theta)) * sign(normal(3),-1.0d0)

    pho%dcos(1)=dcos_temp(1)
    pho%dcos(2)=dcos_temp(2)
    pho%dcos(3)=dcos_temp(3)

! ROTATION along an axis (i.e., the normal)
! compute the dot product of the vector and the rotation axis.
    dot=dot_product(dcos_temp,normal)

! compute the parallel component of the vector.
    dcos_pa= dot * normal

! compute the normal component of the vector.
    dcos_nor(1)= dcos_temp(1)-dcos_pa(1)
    dcos_nor(2)= dcos_temp(2)-dcos_pa(2)
    dcos_nor(3)= dcos_temp(3)-dcos_pa(3)
    dcos_nor=normalize(dcos_nor)

! compute a second vector, lying in the plane, perpendicular
    cross(1)=normal(2)*dcos_nor(3)-normal(3)*dcos_nor(2)
    cross(2)=normal(3)*dcos_nor(1)-normal(1)*dcos_nor(3)
    cross(3)=normal(1)*dcos_nor(2)-normal(2)*dcos_nor(1)
    cross=normalize(cross)

! rotated vector is the parallel component plus the rotated component.
    pho%dcos(1)=dcos_pa(1) + cos(dot)* dcos_nor(1) + sin(dot)*cross(1)
    pho%dcos(2)=dcos_pa(2) + cos(dot)* dcos_nor(2) + sin(dot)*cross(2)
    pho%dcos(3)=dcos_pa(3) + cos(dot)* dcos_nor(3) + sin(dot)*cross(3)
    pho%dcos=normalize(pho%dcos)

! END ROTATION
  return

end subroutine difuser
!---------------------------------------------------------------------
subroutine mirror

! model of specular reflector. It may absorb a fraction absmirfrac.
    use global_variables
    use mtpar
    implicit none

! determine if absorption occurs with random number
    status=rmt(r)
    if (r<absmirfrac) then
        flagdie=.true.
        return
    end if

! compute reflected vector
    pho%dcos=2.*dot_product(-pho%dcos,normal)*normal + pho%dcos  ! specular ray
    return

end subroutine mirror
!---------------------------------------------------------------------
subroutine collana

! subroutine collana determines collision type
    use global_variables
    use mtpar
    integer*4 kpar,ibody,mat,ilb
    real*8 e,x,y,z,u,v,w,wght
    common/track/e,x,y,z,u,v,w,wght,kpar,ibody,mat,ilb(5)

! determine type of collision with random number
! note that fabs is the fractional probability of
! absorption (abs/abs+scat)
    if (fabs(mat)<=2e-12) then
        ab=2
    else
        status=rmt(r)
        if (r<=fabs(mat)) then
            ab=1        ! ABSORPTION EVENT
        else
            ab=2        ! SCATTERING EVENT
        end if
    end if
    return

end subroutine collana
!----------------------------------------------------------------------
subroutine distcoll

! subroutine distcoll computes the distance to collision based
! on a total cross section
    use global_variables
    use mtpar
    implicit none

    integer*4 kpar,ibody,mat,ilb
    real*8 e,x,y,z,u,v,w,wght
    common/track/e,x,y,z,u,v,w,wght,kpar,ibody,mat,ilb(5)

    status=rmt(r)

! with precision (9,38), log function must be restricted
! test to control overflow in the log function
! this might be machine-dependent. Problem may be solved by using
! qlog instead of log (a kind=16 function), or dlog (a kind=8
! type function). If no such precision is needed or wanted, an
! if test is faster overall than using higher precision.
! However, by modifying the rng to include 1 and exclude 0
! the normal log function can be used. This is done in version2

    if (pmat(mat,5)<=1e-12) then
        dcoll=1e12
    else
        dcoll=-1./pmat(mat,5)*log(r)
    end if
    return

end subroutine distcoll
!---------------------------------------------------------------------
subroutine fresnels

! subroutine fresnels computes the probabilities for fresnel
! reflection according to incidence angle and polarization
! and returns a flag to indicate the outcome when compared to a
! random number
    use global_variables
    use mtpar
    implicit none
    real*8 :: sindiff,sinsum,tandiff,tansum
    real*8 :: adiff,asum  ! diff and sum of ain and aout
    real*8 :: fsum

! compute probability of reflection using Fresnel's law
! first, remove problem of ain=aout=0
    ain=acos(cosain)
    aout=asin(sinaout)
    if (ain==aout) then
        if (ain==0.) then
            fresnel=pmat(mat2,1)/pmat(mat1,1)-1.
            fresnel=fresnel/(pmat(mat2,1)/pmat(mat1,1)+1.)
            fresnel=fresnel**2
        else
            fresnel=0.
        end if
    else

! compute polar_perp as the component of light polarized along the
! perpendicular direction with respect to the incidence plane.
        adiff=ain-aout
        asum=ain+aout
        sindiff=sin(adiff)
        sinsum=sin(asum)
        tandiff=tan(adiff)
        tansum=tan(asum)

! compute fresnel probabilities of reflection for parallel & perpendicular
        if (sp==1) then
            fresnel=sindiff*sindiff/sinsum/sinsum
            fresnel=fresnel+tandiff*tandiff/tansum/tansum
            fresnel=fresnel*.5
        else
            call perp_fraction
            fresnel_pa=tandiff/tansum
            fresnel_pa=fresnel_pa*fresnel_pa
            fresnel_pe=sindiff/sinsum
            fresnel_pe=fresnel_pe**2

! normalize fresnel to unity
            fsum=sqrt(fresnel_pe**2+fresnel_pa**2)
            fresnel_pa=fresnel_pa/fsum
            fresnel_pe=fresnel_pe/fsum
            if (fresnel_pa > 1.) fresnel_pa=1.
            if (fresnel_pe > 1.) fresnel_pe=1.

! add them according to weight of the components
! note: to find the parallel component, use para+perp^2=1
            para=1.-pho%polar_perp*pho%polar_perp
            fresnel=para*fresnel_pa
            fresnel=fresnel+pho%polar_perp*pho%polar_perp*fresnel_pe
            fresnel=fresnel/(pho%polar_perp*pho%polar_perp+para)
        end if
    end if
    flaga=.false.
    status=rmt(r)
    if (fresnel>=r) flaga=.true. ! reflection occurs, set flaga
    return

end subroutine fresnels
!---------------------------------------------------------------------
subroutine perp_fraction

! computes the perpendicular component with respect to current plane
    use global_variables
    implicit none

! compute the perp component using the normal
    pho%polar_perp=dot_product(pho%polar,normal)

! to compute parallel component using para^2+perp^2=1
    para=sqrt(abs(1.-pho%polar_perp*pho%polar_perp))

! Note: pho%polar_perp has to be positive
! a negative value in polar is equivalent to rotate the
! polarization vector by pipen, which does not make any difference
! for computing the reflection coefficient for both components.

    if (pho%polar_perp<0.) pho%polar_perp=-pho%polar_perp
    return

end subroutine perp_fraction
!---------------------------------------------------------------------
subroutine material
! subroutine material read material properties

    use global_variables
    implicit none
    integer :: wavesca              ! loop wavelength
    real*8 :: musca               ! scatt. coeff. at wavelength wavesca
    character(len=30), dimension(:), allocatable :: mater  ! material file
    character(len=30), dimension(:), allocatable :: mati   ! path-corr file
    integer :: finduf

! open file
    
    fileunit2=finduf()
    open(fileunit2,file=trim(infile)//'.d2m',status='old')     ! material descriptor
    read(fileunit2,*);read(fileunit2,*);read(fileunit2,*);read(fileunit2,*)
    read(fileunit2,*) matnumber

! allocate pmat, mati and mater
    allocate(wpoints(0:matnumber))
    allocate(pmat(0:matnumber,5))
    allocate(wavepoints(0:matnumber,30))
    allocate(edprop(0:matnumber,5,30))
    allocate(mater(matnumber))
    allocate(mati(matnumber))
    allocate(emissive_flag(matnumber))
    allocate(fabs(0:matnumber))

! initialization
    wpoints=0;pmat=0.;wavepoints=0;edprop=0.;fabs=0.
    read(fileunit2,*)
    read(fileunit2,*) (mater(j),j=1,matnumber)
    close(fileunit2)
    
! re-direct path
    mati=mater

! open material files needed and read properties
    print*, 'materials for optical transport'
    do j=1,matnumber
        k=20+j
        print*,j,mati(j)
        open(k,file=mati(j),status='old')
! advance 23 lines to the first optical properties line and read
        read(20+j,*);read(20+j,*);read(20+j,*);read(20+j,*);read(20+j,*)
        read(20+j,*);read(20+j,*);read(20+j,*);read(20+j,*);read(20+j,*)
        read(20+j,*);read(20+j,*);read(20+j,*);read(20+j,*);read(20+j,*)
        read(20+j,*);read(20+j,*);read(20+j,*);read(20+j,*);read(20+j,*)
        read(20+j,*);read(20+j,*)
        read(20+j,*) wpoints(j)

! check number of points is less or equal than 30
    if (wpoints(j)>30) then
        print*,'Use no more than 30 wavelengths in material files'
        stop
    end if
    do k=1,wpoints(j)
        read(20+j,*) wavepoints(j,k),edprop(j,1,k)
    end do
    read(20+j,*)

! if index is more than 9, then divide by 10!
    do k=1,wpoints(j)
        if (edprop(j,1,k) >  99.) then
            edprop(j,1,k)=edprop(j,1,k)/100.
        else if (edprop(j,1,k) >  9.)  then
            edprop(j,1,k)=edprop(j,1,k)/10.
        end if
    end do
    do k=1,wpoints(j)
          read(20+j,*) wavepoints(j,k),edprop(j,3,k)
    end do
        read(20+j,*)
        read(20+j,*) wavesca
        read(20+j,*) musca
        read(20+j,*) emissive_flag(j)

! build the scattering coefficient table with the lambda^-4 dependency
    do k=1,wpoints(j)
        edprop(j,4,k)=musca/wavepoints(j,k)/wavepoints(j,k)
        edprop(j,4,k)=edprop(j,4,k)*wavesca*wavesca
        edprop(j,4,k)=edprop(j,4,k)/wavepoints(j,k)/wavepoints(j,k)
        edprop(j,4,k)=edprop(j,4,k)*wavesca*wavesca
    end do
        close(20+j)
    end do

! define material zero as vacuum/air with constant properties
    edprop(0,3:5,:)=1e-12
    edprop(0,1,:)=1.
    edprop(0,2,:)=1.
    pmat(0,1)=1.;pmat(0,3)=0.;pmat(0,4)=0.;pmat(0,5)=0.

! build material cross section tables
    edprop(:,5,:)=edprop(:,3,:)+edprop(:,4,:)
    deallocate(mater,mati)

end subroutine material
!------------------------------------------------------------
SUBROUTINE gasdev_s(harvest)

    use global_variables
    use mtpar
    IMPLICIT NONE
        real*8, INTENT(OUT) :: harvest
!    Returns in harvest a normally distributed deviate with zero mean and
!    unit variance, using source of uniform deviates.
    real*8 :: rsq
    real*8, DIMENSION(1) :: v1,v2
    real*8, SAVE :: g
    LOGICAL, SAVE :: gaus_stored=.false.

    if (gaus_stored) then                ! We have an extra deviate handy,
        harvest=g                        ! so return it,
        gaus_stored=.false.              ! and unset the flag
    else                                 ! We don't have extra deviate, so do
        do
            status=rmt(r)
            v1(1)=r                 ! two uniform numbers in -1 to +1
            status=rmt(r)
            v2(1)=r 
            v1(1)=2.0*v1(1)-1.0
            v2(1)=2.0*v2(1)-1.0
            rsq=v1(1)**2+v2(1)**2        ! see if they are in the unit circle,
            if (rsq > 0.0 .and. rsq < 1.0) exit
        end do                           ! otherwise try again.
        rsq=sqrt(-2.D0*log(rsq)/rsq)      ! make Box-Muller transformation to
                                         ! get two normal deviates.
        harvest=v1(1)*rsq
        g=v2(1)*rsq
        gaus_stored=.true.                ! Set flag.
    end if

END SUBROUTINE gasdev_s
!---------------------------------------------------------------------
subroutine detect2report(dii_n,dii_cputim,dii_uncert)
! report generates file with final results

    use global_variables
    implicit none
    integer*4 dii_n
    real*8 dii_cputim, dii_uncert
    real*4 :: fwhm
    real*4, dimension(:,:), allocatable :: nd3d  ! normalized d3d
    real*4, dimension(:), allocatable :: nd1d  ! normalized d3d
    real*8, dimension(bins) :: binsum        ! binrad
    integer, allocatable, dimension(:) :: dethist
    real*4, dimension(:), allocatable :: scasum  ! normalized scatterings/h
    integer finduf
    integer maxcollect

    allocate(scasum(cmax))

! preliminary report calculations, if no photons counted, warning and stop
    if (sum(detec) == 0) then
        print*,'No photons counted (no d2 output at ',xm,' histories).'
        return
    end if

! calculate Swank factor
    allocate(dethistvar(0:maxval(detec)))
    dethistvar=0
    do j=1,ubound(detec,1)
        dethistvar(detec(j))=dethistvar(detec(j))+1
    end do
    mom0=0.;mom1=0.;mom2=0.
    do j=1,maxval(detec)
        mom0=mom0+ dethistvar(j)
        mom1=mom1+ j*dethistvar(j)
        mom2=mom2+ j*j*dethistvar(j)
    end do
    deallocate(dethistvar)
    sf=mom1*mom1/mom0/mom2
    if ( sf  /= sf_old) then
        errorsf=(sf-sf_old)/sf
        totalsf=(totalsf+sf)/2.
        squarredsf=(squarredsf+sf**2)/2.
        sf_old=sf
        dii_uncert=sf
     end if

! write output
    fileunit=finduf()
    open(fileunit,file=trim(infile)//'case.out')
    write(fileunit,*)'Case options:'
    write(fileunit,*)sg,sd,se,sp
    write(fileunit,*)'dii statistics:'
    write(fileunit,*)'# of histories                                 ',xmp
    write(fileunit,*)'# of primaries in phosphor             ',xm
    write(fileunit,*)'optical photons (on average)       ',&
        & real(sum(produced))/(real(xm))
    write(fileunit,*)'optical photons (total)                  ',sum(produced)
    write(fileunit,*)'absorbed                            ',absorbed
    write(fileunit,*)'killed in surface                   ',died
    write(fileunit,*)'counted                             ',counted
    write(fileunit,*)'photons emerged outside screen      ',emerged
    write(fileunit,*)'total number of scattering events   ',scatters
    write(fileunit,*)'mean # scatterings per histories    ', &
        & real(scatters)/sum(produced)
    write(fileunit,*)'mean # absorptions per histories    ', &
        & real(absorbed)/sum(produced)
    write(fileunit,*)'mean distance per history (cm)      ', &
        & meandist/sum(produced)
    write(fileunit,*)'time/track (seconds)[real-time]     ',dii_cputim/sum(produced)
    write(fileunit,*)'collection efficiency [%]           ', &
      & real(counted)/real(sum(produced))*100.
    write(fileunit,*)'lost somewhere (just checking)      ',lost
    write(fileunit,*)'Swank factor          ',sf,'+/-',errorsf
    write(fileunit,*)'m0,m1,m2  ',nint(mom0),nint(mom1),mom2
    write(fileunit,*)'m1/histories (light output per xray) ',mom1/xm
    write(fileunit,*)'xm/total histories (interaction eff) ',real(xm)/real(dii_n)
    
    ! fileunit2=finduf()
    ! open(fileunit2,file=trim(infile)//'raw.out')    ! results
    ! do j=1,bins
        ! write(fileunit2,*) (j-1)*radmax/real(bins),binrad(j)
    ! end do
    ! close(fileunit2)

    ! open(fileunit2,file=trim(infile)//'area.out')    ! results
    ! radmaxtemp=radmax/bins
    ! do j=1,bins
       ! binsum(j)=real(binrad(j))/(pipen*((j*radmaxtemp)* &
    ! & (j*radmaxtemp)-((j-1)*radmaxtemp)*((j-1)*radmaxtemp)))
    ! end do
    ! binsum=binsum/maxval(binsum)
    ! do j=1,bins
       ! write(fileunit2,*) (2*j-2)*radmax*5/bins*1.e4,binsum(j)
       ! write(fileunit2,*) j*radmax/bins*1.e4,binsum(j)
    ! end do
    ! close(fileunit2)

! estimate FWHM
    fw: do j=bins,1,-1
    if (binsum(j) >= .5) then
        fwhm=real(j*radmax)/real(bins)*1.e4 * 2.
    exit fw
    end if
    end do fw
    write(fileunit,*)'estimated FWHM                          ',fwhm
    close(fileunit)
    open(fileunit2,file=trim(infile)//'3det.out')    ! results
    allocate(nd3d(xsize,ysize))
    ! UNCOMMENT FOR NORMALIZED OUTPUT
    nd3d=real(d3d)/real(maxval(d3d))
    if (ysize==1) then
       do j=1,xsize
          write(fileunit2,*) nd3d(j,1)
       end do
    else if (ysize>1 .and. xsize>1) then
        do k=1,xsize
           do j=1,ysize
               write(fileunit2,*) nd3d(j,k)
           end do
           write(fileunit2,*)
        end do
    else if (xsize==1) then
        do j=1,ysize
          write(fileunit2,*) nd3d(1,j)
        end do
    end if
    close(fileunit2)

    ! open(fileunit2,file=trim(infile)//'1det.out')    ! results
    allocate(nd1d(xsize))
    ! do j=1,xsize
        ! nd1d(j)=sum(nd3d(j,:))
    ! end do
    ! nd1d=nd1d/max(1.,maxval(nd1d))
    ! do j=1,xsize
       ! write(fileunit2,*) j,nd1d(j)
    ! end do
    deallocate(nd3d,nd1d)
    ! close(fileunit2)

! collection statistics
! INCREASE IF NEEDED (E.G., MEV IMAGING)
    maxcollect=10000
    allocate(dethist(0:maxcollect))
    ! EQUALIZE LENGHTS allocate(dethist(0:maxval(detec)))
    dethist=0
    do j=1,min(ubound(detec,1),maxcollect)
        dethist(detec(j))=dethist(detec(j))+1
    end do
    ! UNCOMMENT FOR NORMALIZED OUTPUT 
	dethist=dethist/sum(dethist)
    open(fileunit2,file=trim(infile)//'collect.out')
    write(fileunit2,*) "# collection statistics"
    do j=1,min(maxval(detec),maxcollect)
       write(fileunit2,*) j,dethist(j)
    end do
    write(fileunit2,*) '# total ',sum(detec),sum(dethist)
    close(fileunit2)
    open(fileunit2,file=trim(infile)//'collect10.out')
    write(fileunit2,*) "# collection statistics"
    do j=1,9990,10
       write(fileunit2,*) j+5,sum(dethist(j:j+10))
    end do
    write(fileunit2,*) '# total ',sum(detec),sum(dethist)
    close(fileunit2)
    deallocate(dethist)
    return

end subroutine detect2report
!---------------------------------------------------------------------
subroutine scatter

! subroutine scatter returns new photon pos and dcos after scattering
    use global_variables
    use mtpar
    implicit none
    real*8 :: sphi,stheta           ! scattering angles
    real*8 :: rtemp            ! temporary variable
    real*8, dimension(3,3) :: rotsph     ! rotation matrix about z-axis
    real*8, dimension(3,3) :: rotsth     ! rotation matrix about xy-plane
    real*8, dimension(3) :: scadcos      ! scattering dcos

! store initial direction
    pho%dcos_old=pho%dcos
! sampling for Rayleigh scattering distribution function
! p(x,y)=1-cos(x)**2sin(y)**2
! where x is the scattering angle theta and
! y is the azymuthal angle phi as defined in Kerker, p.33
! use conditional probability pdf(x|y)

! first, sample y from p(y) resulting from integrating out
! the x dependencies in p(x,y)
    do
        status=rmt(r)
        r=r*pipen
        sphi=r
        rtemp=sin(r)**2
        rrr=(2.-rtemp)/2.
        status=rmt(r)
        if (r<rrr) exit
    end do

    do
        status=rmt(r)
        r=r*pipen
        stheta=r
        rrr=(1.-cos(r)**2*rtemp)
        status=rmt(r)
        if (r<rrr) exit
    end do

! then define vector with scattering angles
    rrr=sin(sphi)
    scadcos(1)=rrr*cos(stheta)
    scadcos(2)=rrr*sin(stheta)
    scadcos(3)=cos(sphi)

! computation of phi and theta
    phipen=acos(pho%dcos(3))
    sinphi=sin(phipen)
    if (sinphi/=0.) then
        if (abs(pho%dcos(2)/sinphi)>=1) then
            if (abs(pho%dcos(1)/sinphi)>=1) then
        status=rmt(r)
        theta=r*twopipen
            else
                theta=acos(pho%dcos(1)/sinphi)
            end if
        else
            theta=asin(pho%dcos(2)/sinphi)
        end if
    else
    status=rmt(r)
    theta=r*twopipen
    end if

! then rotate to original direction using initial angles
    rotsph=0; rotsph(2,2)=1; rotsph(1,1)=pho%dcos(3)
    rotsph(1,3)=sin(phipen); rotsph(3,1)=-rotsph(1,3)
    rotsph(3,3)=rotsph(1,1)
    rotsth=0; rotsth(3,3)=1; rotsth(1,1)=cos(theta)
    rotsth(2,2)=rotsth(1,1); rotsth(1,2)=-sin(theta)
    rotsth(2,1)=rotsth(1,2)

! rotation according to direction of propagation
    pho%dcos=matmul(rotsph,scadcos)
    pho%dcos=matmul(rotsth,pho%dcos)

! computate phi and theta for polarization vector
    phipen=acos(pho%polar(3))
    sinphi=sin(phipen)
    if (sinphi/=0.) then
        if (abs(pho%polar(2)/sinphi)>=1) then
            if (abs(pho%polar(1)/sinphi)>=1) then
                status=rmt(r)
                theta=r*twopipen
            else
                theta=acos(pho%polar(1)/sinphi)
            end if
        else
            theta=asin(pho%polar(2)/sinphi)
        end if
    else
        status=rmt(r)
        theta=r*twopipen
    end if

! define rotation based on polarization vector
    rotsph=0; rotsph(2,2)=1; rotsph(1,1)=pho%polar(3)
    rotsph(1,3)=sin(phipen); rotsph(3,1)=-rotsph(1,3)
    rotsph(3,3)=rotsph(1,1)
    rotsth=0; rotsth(3,3)=1; rotsth(1,1)=cos(theta)
    rotsth(2,2)=rotsth(1,1); rotsth(1,2)=-sin(theta)
    rotsth(2,1)=rotsth(1,2)

! rotation of directional cosines according to polarization
    pho%dcos=matmul(rotsph,scadcos)
    pho%dcos=matmul(rotsth,pho%dcos)

! normalize directional cosines
    pho%dcos=normalize(pho%dcos)

! finally, if needed, recompute the polarization vector
    if (sp/=1) then

! define normal to the plane of scattering
        pho%polar=nor(pho%dcos,pho%dcos_old)
        pho%polar=normalize(pho%polar)
    else
        print*,'Handling of Rayleigh scattering without polarization is not valid '
        print*,'(see Van de Hulst, Multiple light scattering)'
        print*,'Please, re-state your problem in input'
        stop
    end if

    return

end subroutine scatter
!---------------------------------------------------------------------
subroutine isoscatter

! subroutine scatter returns new photon pos and dcos after scattering
    use global_variables
    use mtpar
    implicit none
    pho%dcos_old=pho%dcos
    status=rmt(r)
    r=r*2.-1.
    pho%dcos(3)=r
    rr=sqrt(1.-r*r)
    status=rmt(r)
    theta=r*twopipen
    pho%dcos(1)=rr*cos(theta)
    pho%dcos(2)=rr*sin(theta)
    return

end subroutine isoscatter
!---------------------------------------------------------------------
subroutine inidetect2
! old subroutine setup
! -jobname contains the *-character jobname for IO

    use global_variables
    use mtpar
    implicit none
    external material,surface
    integer*4 nhist,atime,refres,ncalls,lastim
    real*8 time0,accura
    integer :: seedpen1,seedpen2
	common /ctrsim/ time0,accura,atime,refres,ncalls,lastim,nhist
	
	integer*4 seed1,seed2
    common/rseed/seed1,seed2

!!!!! MODIFICATION FOR DEPTH OF INTERACTION
    character*4 arg1 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
    integer :: finduf
    character*100 jobname
    common /cjob/ jobname

    infile=jobname

!!!!! MODIFICATION FOR DEPTH OF INTERACTION
    call getarg(2,arg1)
    read(arg1,*) depth_x
    call getarg(3,arg1)
    read(arg1,*) depth_y
    call getarg(4,arg1)
    read(arg1,*) depth
    depth_x=depth_x*1.0000e-4
    depth_y=depth_y*1.0000e-4
    depth=depth*1.0000e-4
    print*,depth_x,depth_y,depth
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ehgen=1.;bins=100;lost=0;emerged=0;xmp_old=0;xmp=0;xm=0
    allocate(binrad(bins))


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! initialize rng  USE PENELOPE's SEEDS? 
    print*,'using PENELOPE seeds #1 for the MT initialization'
    print*,seed1
    print*,'and #2 for the sub-sequence (this should reflect sub-job unit ID'
    print*,seed2
	call sgrnd(seed1)    
	! Intel intrinsic MT status =  irmt(seed1,seed2)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! open file
    fileunit=finduf()
    open(fileunit,file=trim(infile)//'.d2c',status='old')  ! running parameters

! read running MC parameters pmax,cmax from namelist
    rewind(fileunit);read(fileunit,nml=mcparameters)
    rewind(fileunit);read(fileunit,nml=variance)

! if no optical photons just return
    if (pmax_mu == 0.) return

! read source parameters sg,sd,se,sp from namelist
    rewind(fileunit);read(fileunit,nml=sourcetype)
    if (sd==41) then
        rewind(fileunit);read(fileunit,nml=incidence_angle)
        iangle=iangle*pipen/180.
    end if

! read scattering flag
    rewind(fileunit);read(fileunit,nml=scattering_type)

! read binning angle phimax and angle of restriction for rough surfaces
    rewind(fileunit);read(fileunit,nml=binangle)

! confirm beta is <= 60 degrees
    if (beta>60.) then
        print*,'Please, use beta  <= 60 degrees !'
        stop
    end if

! read parameters related to surface type and sources
    rewind(fileunit);read(fileunit,nml=surftypes)

! correct for percentage due to script handling
    absdiffrac=absdiffrac/100.
    abssmoothfrac=abssmoothfrac/100.
    absroughfrac=absroughfrac/100.
    absthinfrac=absthinfrac/100.
    absmirfrac=absmirfrac/100.
    detrefl=detrefl/100.
    close(fileunit)

! read material descriptor
    call material

! read surface non-default descriptor
    call surface

! read source multiple line wavelength distribution if se=3
    select case (se)
    case (1)
        fileunit2=finduf()
        open(fileunit2,file=trim(infile)//'.d2d',status='old')
        read(fileunit2,nml=binwave); rewind(fileunit2)
! allocate wspect and pspect with bwave
        allocate(wspect(bwave),pspec(bwave),pspect(bwave),pwave(bwave))

! read spectrum distribution
        read(fileunit2,*);read(fileunit2,*);read(fileunit2,*)
        read(fileunit2,*);read(fileunit2,*)
        read(fileunit2,*);read(fileunit2,*);read(fileunit2,*)
        read(fileunit2,*);read(fileunit2,*)
        do j=1,bwave
            read(fileunit2,*) wspect(j),pspect(j)
        end do
        close (fileunit2)
        if (sum(pspect)>1.1) then
            do iw=1,bwave
                pspec(iw)=sum(pspect(1:iw))
            end do
            maxpec=maxval(pspec)
            pspec=pspec/maxpec
        else
            pspec=pspect/maxval(pspect)
        end if
        deallocate(pspect)

        fileunit2=finduf()
        open(fileunit2,file=trim(infile)//'.d2x',status='old')
        read(fileunit2,nml=binsensorwave); rewind(fileunit2)

! allocate wspect and pspect with bwave
        allocate(wsensor(bsensorwave),psensor(bsensorwave))

! read spectrum distribution
        read(fileunit2,*);read(fileunit2,*);read(fileunit2,*)
        read(fileunit2,*);read(fileunit2,*)
        read(fileunit2,*);read(fileunit2,*);read(fileunit2,*)
        read(fileunit2,*);read(fileunit2,*)
        do j=1,bsensorwave
            read(fileunit2,*) wsensor(j),psensor(j)
        end do
        close (fileunit2)

    end select

    fileunit=finduf()
    open(fileunit,file=trim(infile)//'.d2c',status='old')  ! running parameters

! read binning types wanted
    rewind(fileunit);read(fileunit,nml=binningtype)
    rewind(fileunit);read(fileunit,nml=binarraysize)
    imgdx=imgdx*1.e-4
    imgdy=imgdy*1.e-4
    rewind(fileunit);read(fileunit,nml=binangle)
    close(fileunit)

! should it be larger due to secondaries?
    allocate(detec(1:nhist))
    allocate(produced(1:nhist))
    detec=0;produced=0
    allocate(d3d(xsize,ysize)); d3d=0
    allocate(detlsf(xsize));    detlsf=0    

! initialize counters and flags
    flaga=.true.;flagdie=.false.;absorbed=0;ab=0
    counted=0;died=0;binrad=0;detected=0
    scatters=0
    meandist=0.;d3d=0;pwave=0
    radmax=max(xsize*imgdx/2.,ysize*imgdy/2.)

! convert phimax from degrees to radians
    if (phimax==0.) phimax=90.
    phimaxr=phimax*pipen/180.
    scatters=0

! OPTION reportplot
    if (reportplot .eqv. .true.) then
        open(101,file=trim(infile)//'1xray_events.out')
        open(102,file=trim(infile)//'2xray_events.out')
        open(103,file=trim(infile)//'3xray_events.out')
        open(104,file=trim(infile)//'4xray_events.out')
        open(105,file=trim(infile)//'5xray_events.out')
        open(111,file=trim(infile)//'1optical_events.out')
        open(112,file=trim(infile)//'2optical_events.out')
        open(113,file=trim(infile)//'3optical_events.out')
        open(114,file=trim(infile)//'4optical_events.out')
        open(115,file=trim(infile)//'5optical_events.out')
    end if
    print*,'iniDetect2 done...'
    return

end subroutine iniDetect2
!---------------------------------------------------------------------
subroutine dii_source
! subroutine source returns all the properties
! of the initial photon. It samples the angles and computes
! the directional cosines, samples the wavelength and the polarization.
! This is done according to flags for different source types.

    use global_variables
    use mtpar
    implicit none

! select source directionallity using flag variable sd defined in input.dat
! select source geometry using flag variable sg defined in input.dat

    pho%pos=ipos
    select case (sd)
       case (1)
         call isotropic
       case (3)
         call isoall
    end select
    pho%dcos=normalize(pho%dcos)

! sample the polarization component polar (perpendicular)

    select case (sp)
    case (1)
! use equal polarization components, which will determine the average
! Fresnel's equation. As it will not be used, assing 0. to eventually check it
    case (2)
! polarization vector has to be perpendicular to the ray direction.
! Thus, this condition implies that sampling 2 out of 3 coordinates
! will be enough to define it. Then, define randomly the
! perpendicular component perp
       status=rmt(r)
       pho%polar(1)=2*r-1
       status=rmt(r)
       pho%polar(2)=2*r-1
       pho%polar(3)=pho%polar(1)*pho%dcos(1)+pho%polar(2)*pho%dcos(2)
       pho%polar(3)=pho%polar(3)/max(1.0d-12,pho%dcos(3))
       pho%polar=normalize(pho%polar)
       status=rmt(r)
       pho%polar_perp=r
    end select

! sample the photon wavelength
    select case (se)
    case (1)
      pho%wavel_old=pho%wavel
      status=rmt(r)
      sample_wave: do j=1,bwave
        if (r<=pspec(j)) then
            pho%wavel=wspect(j)
            exit sample_wave
        end if
      end do sample_wave
! using the selected wavelength, assign to pmat the properties of
! all materials involved in the case, from the energy dependent properties
! if it's different than the previous
      if (pho%wavel/=0. .and. pho%wavel_old/=pho%wavel) call mat_properties

    case (2)
      pho%wavel_old=pho%wavel
      status=rmt(r)
      pho%wavel=400.+r*300.
! using the selected wavelength, assign to pmat the properties of
! all materials involved in the case, from the energy dependent properties
      if (pho%wavel/=0. .and. pho%wavel_old/=pho%wavel) call mat_properties

    case (3)
      pho%wavel_old=pho%wavel
      status=rmt(r)
      pho%wavel=400.+r*1600.
! using the selected wavelength, assign to pmat the properties of
! all materials involved in the case, from the energy dependent properties
      if (pho%wavel/=0. .and. pho%wavel_old/=pho%wavel) call mat_properties

    end select

! initialize # of events
    pho%events=0
    return

end subroutine dii_source
!---------------------------------------------------------------------
subroutine mat_properties

! given the photon wavelength assigns the material properties
    use global_variables
    implicit none

    select case (se)
    case (1)

! loop over all materials and get material properties for this wavelength
        do j=1,matnumber

! locate range for interpolation in [400,800] (nm)
        locatew: do k=1,wpoints(j)
            if (pho%wavel<=wavepoints(j,k)) then
               e2=k
               exit locatew
            end if
        end do locatew
        if (k==1) then
            e1=e2
        else
            e1=k-1
        end if
        if (e2==0) e2=k-1

! perform linear interpolation, noting the singular case of equal wavelengths
        if (wavepoints(j,e2)<=wavepoints(j,e1)) then
            pmat(j,:)=edprop(j,:,e1)
            if (pmat(j,5)/=0.) then
                fabs(j)=pmat(j,3)/pmat(j,5)
            else
                fabs(j)=1e-12
            end if
        else
            pmat(j,:)=(edprop(j,:,e2)-edprop(j,:,e1))
            pmat(j,:)=pmat(j,:)/(wavepoints(j,e2)-wavepoints(j,e1))
            pmat(j,:)=pmat(j,:)*(pho%wavel-wavepoints(j,e1))
            pmat(j,:)=edprop(j,:,e1)+pmat(j,:)

! also compute the absorption fraction as (abs/(abs+sca))
            if (pmat(j,5)/=0.) then
                fabs(j)=pmat(j,3)/pmat(j,5)
            else
                fabs(j)=1e-12
            end if
        end if
    end do

    case (2,3)

! loop over all materials and get material properties for this wavelength

    do j=1,matnumber
! assumes that index of refraction is not dependent on wavelength
! and that only change is into the scattering coefficient musca
! according to lambda^-4
           pmat(j,1)=edprop(j,1,1)
           pmat(j,3)=edprop(j,3,1)
           pmat(j,4)=edprop(j,4,1)*(wavepoints(j,1)/pho%wavel)**4
           pmat(j,5)=pmat(j,3)+pmat(j,4)
! also compute the absorption fraction as (abs/(abs+sca))
           if (pmat(j,5)/=0.) then
              fabs(j)=pmat(j,3)/pmat(j,5)
           else
              fabs(j)=1e-12
           end if
       end do
    end select
    return

end subroutine mat_properties
!---------------------------------------------------------------------
subroutine isotropic

! sample the angles phi and theta and compute directional cosines
    use global_variables
    use mtpar
    implicit none
    status=rmt(r)
    pho%dcos(3)=r
    rr=sqrt(1.-r*r)
    status=rmt(r)
    theta=r*twopipen
    pho%dcos(1)=rr*cos(theta)
    pho%dcos(2)=rr*sin(theta)
    return

end subroutine isotropic
!---------------------------------------------------------------------
subroutine isoall

! sample angles phi and theta and compute directional cosines for isotropic
! source in the complete 4PI space
    use global_variables
    use mtpar
    implicit none
    status=rmt(r)
    r=2.*r-1.
    pho%dcos(3)=r
    rr=sqrt(1.-r*r)
    status=rmt(r)
    theta=r*twopipen
    pho%dcos(1)=rr*cos(theta)
    pho%dcos(2)=rr*sin(theta)
    return

end subroutine isoall
!---------------------------------------------------------------------
subroutine surface

! subroutine surface reads non-default surface types from surface.in
    use global_variables
    implicit none
    integer :: finduf

! open file
    fileunit2=finduf()
    open(fileunit2,file=trim(infile)//'.d2s',status='old') ! surface descriptor
! FORMAT OF SURFACE.IN
! MATN MAT11         MAT2                 ...         MATN
! MAT1 0                 SURF1-2         ...         SURF1-N
! MAT2 SURF2-1         0                 ...         SURF2-N
! ...
! MATN SURFN-1        SURFN-2         ...        0
    read(fileunit2,*);read(fileunit2,*);read(fileunit2,*);read(fileunit2,*)
    read(fileunit2,*);read(fileunit2,*);read(fileunit2,*)
    read(fileunit2,*) matn
    read(fileunit2,*)
    allocate(surftype(1:matn,1:matn))
    surftype=1
    do i=1,matn
      read(fileunit2,*) temp,(surftype(j,i) , j=1,matn)
    end do
    close(fileunit2)
    return

end subroutine surface
!---------------------------------------------------------------------
subroutine tallyd2

! subroutine tallyd2 decides if the photon has been succesfull
! and if so convert location to absolute and bin
    use global_variables
    use poisson
    implicit none

    integer*4 kpar,ibody,mat,ilb
    real*8 e,x,y,z,u,v,w,wght
    common/track/e,x,y,z,u,v,w,wght,kpar,ibody,mat,ilb(5)

! to add variance reduction with interaction forcing multiply
! all quanitites with the weight "wght"

    countf = .FALSE.
    call sensor
    if (countf .eqv. .FALSE.) then
        absorbed=absorbed+1
        if (be==2) call binevscat
        flaga=.TRUE.
    else if (countf .eqv. .TRUE.) then
        counted=counted+1
        detec(xm)=detec(xm)+1
        flaga=.TRUE.

! get ehgen from Poisson deviate as P(lambda)
! with E(eV)=1240/lambda(nm) and 3.6 eV/ehp
        ehgen=nint(poidev(real(pho%wavel/344.44)))
        k= nint(y/real(imgdy) + ysize/2. )
        j= nint(x/real(imgdx) + xsize/2. )
        if (j <= ubound(d3d,1) .and. k <= ubound(d3d,2) .and. &
           & j >= lbound(d3d,1) .and. k >= lbound(d3d,2)) then
             d3d(j,k)=d3d(j,k)+ehgen     ! *pho%wavel
        end if
        radius=sqrt(x**2+y**2)
        j=ceiling(radius/radmax*bins)
        if (j <= ubound(binrad,1) .and.  &
           & j >= lbound(binrad,1)) binrad(j)=binrad(j)+1
    end if
    return

end subroutine tallyd2
!----------------------------------------------------------------------
subroutine sensor

! subroutine sensor decides if the optical photon has been succesfull
! generating a count
    use global_variables
    use mtpar
    implicit none
    if (bsensorwave == 1) then
        countf = .TRUE.
    else
        snsr: do i=1,bsensorwave
            if (wsensor(i) <= pho%wavel) then
               status=rmt(r)
               if (r <= psensor(i)) countf = .true.
               exit snsr
            end if
        end do snsr
    end if
    return

end subroutine sensor
!----------------------------------------------------------------------
subroutine binevscat

! subroutine binevscat bins the number of scattering events
    use global_variables
    use mtpar
    implicit none
    j=max(1,floor(real(scatters)/cmax))
    scatev(j)=scatev(j)+1
    sumscat=sumscat+scatters
    return

end subroutine binevscat
!----------------------------------------------------------------------
