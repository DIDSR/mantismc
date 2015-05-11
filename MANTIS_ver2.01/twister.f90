!----------------------------------------------------------------------!
! MANTIS RNG: FUNCTION FOR MULTIPLE RANDOM NUMBER GENERATORS IN MANTIS !
!                                                                      !
! 1. Replace FUNCTION RAND in PENELOPE with the one below              !
! 2. Select the choice rng and uncomment that "function rmt"           !
!   a-penelope's rngs                                                  !
!   b-mersenne twister (coded)                                         !
!   c-mersenne twister (Intel intrinsic function)                      !
!                                                                      !
!   Note 1: option a and b require initialization call in mantis.f     !
!   Note 2: option b requires additional code in mt_intel.f90          !
!----------------------------------------------------------------------!
C  *********************************************************************
C                         FUNCTION RAND
C  *********************************************************************
      REAL*8 FUNCTION RAND(DUMMY)
	  real*8 rand,dummy
	  integer status
      status=rmt(rand)
	  RETURN
      END
!----------------------------------------------------------------------!
!   a-penelope's rngs                                                  !
!----------------------------------------------------------------------!
C  *********************************************************************
C                         FUNCTION RAND
C  *********************************************************************
!     FUNCTION RAND(DUMMY)
C
C  This is an adapted version of subroutine RANECU written by F. James
C  (Comput. Phys. Commun. 60 (1990) 329-344), which has been modified to
C  give a single random number at each call.
C
C  The 'seeds' ISEED1 and ISEED2 must be initialized in the main program
C  and transferred through the named common block /RSEED/.
C
C  Some compilers incorporate an intrinsic random number generator with
C  the same name (but with different argument lists). To avoid conflict,
C  it is advisable to declare RAND as an external function in all sub-
C  programs that call it.
C
!     IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N)
!     PARAMETER (USCALE=1.0D0/2.147483563D9)
!     COMMON/RSEED/ISEED1,ISEED2
C
!     I1=ISEED1/53668
!     ISEED1=40014*(ISEED1-I1*53668)-I1*12211
!     IF(ISEED1.LT.0) ISEED1=ISEED1+2147483563
C
!     I2=ISEED2/52774
!     ISEED2=40692*(ISEED2-I2*52774)-I2*3791
!     IF(ISEED2.LT.0) ISEED2=ISEED2+2147483399
C
!     IZ=ISEED1-ISEED2
!     IF(IZ.LT.1) IZ=IZ+2147483562
!     RAND=IZ*USCALE
C
!     RETURN
!     END
!----------------------------------------------------------------------!
!   b-mersenne twister (coded)                                         !
!----------------------------------------------------------------------!
! The twister mersenne prng from 
! http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/mt.html
! Fortran translation by Hiroshi Takano.  Jan. 13, 1999.
!   genrand()      -&gt; double precision function grnd()
!   sgenrand(seed) -&gt; subroutine sgrnd(seed)
! Fortran version rewritten as an F90 module and mt state saving and getting
! subroutines added by Richard Woloshyn. (rwww@triumf.ca). June 30, 1999

! Initialization function
     function sgrnd(seed)
     implicit none
!      setting initial seeds to mt[N] using
!      the generator Line 25 of Table 1 in
!      [KNUTH 1981, The Art of Computer Programming
!       Vol. 2 (2nd Ed.), pp102]
! Default seed
     integer, parameter :: defaultsd = 4357
! Period parameters
     integer, parameter :: mtN = 624, mtN1 = mtN + 1
! the array for the state vector
     integer, save, dimension(0:mtN-1) :: mt
     integer, save                   :: mti = mtN1
     integer, intent(in) :: seed
     mt(0) = iand(seed,-1)
     do mti=1,mtN-1
       mt(mti) = iand(69069 * mt(mti-1),-1)
     enddo
     return
     end function sgrnd

! Random number generator function
      REAL*8 FUNCTION rmt(r)
      real*8 :: rmt,r
! Default seed
      integer, parameter :: defaultsd = 4357
! Period parameters
      integer, parameter :: mtN = 624, mtN1 = mtN + 1
      integer, save, dimension(0:mtN-1) :: mt
      integer, save                   :: mti = mtN1
      integer, parameter :: mtM = 397, MATA  = -1727483681
      integer, parameter :: LMASK =  2147483647
      integer, parameter :: UMASK = -LMASK - 1
      integer, parameter :: TMASKB= -1658038656, TMASKC= -272236544
      dimension mag01(0:1)
      data mag01/0, MATA/
      save mag01

      TSHFTU(y)=ishft(y,-11)
      TSHFTS(y)=ishft(y,7)
      TSHFTT(y)=ishft(y,15)
      TSHFTL(y)=ishft(y,-18)
      if(mti.ge.mtN) then
        if(mti.eq.mtN+1) then
          call sgrnd( defaultsd )
        endif
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

!----------------------------------------------------------------------!
!   c-mersenne twister (Intel intrinsic function)                      !
!----------------------------------------------------------------------!
      function irmt(seed,seq)
         use mkl_vsl_type
         use mkl_vsl
         type (vsl_stream_state) :: stream
         integer :: brng,seq
         integer seed      
		 
         brng=vsl_brng_mt2203+seq
         status = vslnewstream( stream, brng,  seed )
      end function irmt 
! 
      function rmt(r)
         use mkl_vsl_type
         use mkl_vsl
         type (vsl_stream_state) :: stream
         integer :: status,n
         real*8 :: a, b, o(1), r
         integer :: method, brng
         a=0.D0
         b=1.D0
         method=0
         n=1
         status = vdrnguniform( method, stream, n, o, a, b)
         r=o(1)
      end function rmt

!! file: mkl_vsl.fi
!!
!!                             INTEL CONFIDENTIAL
!!  Copyright(C) 2005 Intel Corporation. All Rights Reserved.
!!  The source code contained  or  described herein and all documents related to
!!  the source code ("Material") are owned by Intel Corporation or its suppliers
!!  or licensors.  Title to the  Material remains with  Intel Corporation or its
!!  suppliers and licensors. The Material contains trade secrets and proprietary
!!  and  confidential  information of  Intel or its suppliers and licensors. The
!!  Material  is  protected  by  worldwide  copyright  and trade secret laws and
!!  treaty  provisions. No part of the Material may be used, copied, reproduced,
!!  modified, published, uploaded, posted, transmitted, distributed or disclosed
!!  in any way without Intel's prior express written permission.
!!  No license  under any  patent, copyright, trade secret or other intellectual
!!  property right is granted to or conferred upon you by disclosure or delivery
!!  of the Materials,  either expressly, by implication, inducement, estoppel or
!!  otherwise.  Any  license  under  such  intellectual property  rights must be
!!  express and approved by Intel in writing.
!!++
!!  Fortran 90 VSL interface.
!!--

      MODULE MKL_VSL_TYPE


!!++
!!  Definitions for VSL functions return values (errors, warnings)
!!--

!!  "No error" status
      INTEGER VSL_STATUS_OK
      INTEGER VSL_ERROR_OK
      PARAMETER (VSL_STATUS_OK = 0)
      PARAMETER (VSL_ERROR_OK  = 0)

!!  Common errors (-1..-999)
      INTEGER VSL_ERROR_FEATURE_NOT_IMPLEMENTED
      INTEGER VSL_ERROR_UNKNOWN
      INTEGER VSL_ERROR_BADARGS
      INTEGER VSL_ERROR_MEM_FAILURE
      INTEGER VSL_ERROR_NULL_PTR
      PARAMETER (VSL_ERROR_FEATURE_NOT_IMPLEMENTED = -1)
      PARAMETER (VSL_ERROR_UNKNOWN                 = -2)
      PARAMETER (VSL_ERROR_BADARGS                 = -3)
      PARAMETER (VSL_ERROR_MEM_FAILURE             = -4)
      PARAMETER (VSL_ERROR_NULL_PTR                = -5)

!!  RNG errors (-1000..-1999)
!!  brng errors
      INTEGER VSL_ERROR_INVALID_BRNG_INDEX
      INTEGER VSL_ERROR_LEAPFROG_UNSUPPORTED
      INTEGER VSL_ERROR_SKIPAHEAD_UNSUPPORTED
      INTEGER VSL_ERROR_BRNGS_INCOMPATIBLE
      INTEGER VSL_ERROR_BAD_STREAM
      INTEGER VSL_ERROR_BRNG_TABLE_FULL
      INTEGER VSL_ERROR_BAD_STREAM_STATE_SIZE
      INTEGER VSL_ERROR_BAD_WORD_SIZE
      INTEGER VSL_ERROR_BAD_NSEEDS
      INTEGER VSL_ERROR_BAD_NBITS
      PARAMETER (VSL_ERROR_INVALID_BRNG_INDEX      = -1000)
      PARAMETER (VSL_ERROR_LEAPFROG_UNSUPPORTED    = -1002)
      PARAMETER (VSL_ERROR_SKIPAHEAD_UNSUPPORTED   = -1003)
      PARAMETER (VSL_ERROR_BRNGS_INCOMPATIBLE      = -1005)
      PARAMETER (VSL_ERROR_BAD_STREAM              = -1006)
      PARAMETER (VSL_ERROR_BRNG_TABLE_FULL         = -1007)
      PARAMETER (VSL_ERROR_BAD_STREAM_STATE_SIZE   = -1008)
      PARAMETER (VSL_ERROR_BAD_WORD_SIZE           = -1009)
      PARAMETER (VSL_ERROR_BAD_NSEEDS              = -1010)
      PARAMETER (VSL_ERROR_BAD_NBITS               = -1011)

!! abstract stream related errors
      INTEGER VSL_ERROR_BAD_UPDATE
      INTEGER VSL_ERROR_NO_NUMBERS
      INTEGER VSL_ERROR_INVALID_ABSTRACT_STREAM
      PARAMETER (VSL_ERROR_BAD_UPDATE              = -1120)
      PARAMETER (VSL_ERROR_NO_NUMBERS              = -1121)
      PARAMETER (VSL_ERROR_INVALID_ABSTRACT_STREAM = -1122)

!! read/write stream to file errors
      INTEGER VSL_ERROR_FILE_CLOSE
      INTEGER VSL_ERROR_FILE_OPEN
      INTEGER VSL_ERROR_FILE_WRITE
      INTEGER VSL_ERROR_FILE_READ

      INTEGER VSL_ERROR_BAD_FILE_FORMAT
      INTEGER VSL_ERROR_UNSUPPORTED_FILE_VER
      PARAMETER (VSL_ERROR_FILE_CLOSE           = -1100)
      PARAMETER (VSL_ERROR_FILE_OPEN            = -1101)
      PARAMETER (VSL_ERROR_FILE_WRITE           = -1102)
      PARAMETER (VSL_ERROR_FILE_READ            = -1103)

      PARAMETER (VSL_ERROR_BAD_FILE_FORMAT      = -1110)
      PARAMETER (VSL_ERROR_UNSUPPORTED_FILE_VER = -1111)


!!++
!!  CONV/CORR RELATED MACRO DEFINITIONS
!!--
      INTEGER VSL_CONV_MODE_AUTO
      INTEGER VSL_CORR_MODE_AUTO
      INTEGER VSL_CONV_MODE_DIRECT
      INTEGER VSL_CORR_MODE_DIRECT
      INTEGER VSL_CONV_MODE_FFT
      INTEGER VSL_CORR_MODE_FFT
      INTEGER VSL_CONV_PRECISION_SINGLE
      INTEGER VSL_CORR_PRECISION_SINGLE
      INTEGER VSL_CONV_PRECISION_DOUBLE
      INTEGER VSL_CORR_PRECISION_DOUBLE
      PARAMETER (VSL_CONV_MODE_AUTO        = 0)
      PARAMETER (VSL_CORR_MODE_AUTO        = 0)
      PARAMETER (VSL_CONV_MODE_DIRECT      = 1)
      PARAMETER (VSL_CORR_MODE_DIRECT      = 1)
      PARAMETER (VSL_CONV_MODE_FFT         = 2)
      PARAMETER (VSL_CORR_MODE_FFT         = 2)
      PARAMETER (VSL_CONV_PRECISION_SINGLE = 1)
      PARAMETER (VSL_CORR_PRECISION_SINGLE = 1)
      PARAMETER (VSL_CONV_PRECISION_DOUBLE = 2)
      PARAMETER (VSL_CORR_PRECISION_DOUBLE = 2)


!!++
!!  BASIC RANDOM NUMBER GENERATOR (BRNG) RELATED MACRO DEFINITIONS
!!--


!!  MAX NUMBER OF BRNGS CAN BE REGISTERED IN VSL
!!  No more than VSL_MAX_REG_BRNGS basic generators can be registered in VSL
!!  (including predefined basic generators).
!!
!!  Change this number to increase/decrease number of BRNGs can be registered.
      INTEGER VSL_MAX_REG_BRNGS
      PARAMETER (VSL_MAX_REG_BRNGS = 512)

!!  PREDEFINED BRNG NAMES
      INTEGER VSL_BRNG_SHIFT
      INTEGER VSL_BRNG_INC

      INTEGER VSL_BRNG_MCG31
      INTEGER VSL_BRNG_R250
      INTEGER VSL_BRNG_MRG32K3A
      INTEGER VSL_BRNG_MCG59
      INTEGER VSL_BRNG_WH
      INTEGER VSL_BRNG_SOBOL
      INTEGER VSL_BRNG_NIEDERR
      INTEGER VSL_BRNG_MT19937
      INTEGER VSL_BRNG_MT2203
      INTEGER VSL_BRNG_IABSTRACT
      INTEGER VSL_BRNG_DABSTRACT
      INTEGER VSL_BRNG_SABSTRACT

      PARAMETER (VSL_BRNG_SHIFT=20)
      PARAMETER (VSL_BRNG_INC=ISHFT(1, VSL_BRNG_SHIFT))

      PARAMETER (VSL_BRNG_MCG31    =VSL_BRNG_INC)
      PARAMETER (VSL_BRNG_R250     =VSL_BRNG_MCG31    +VSL_BRNG_INC)
      PARAMETER (VSL_BRNG_MRG32K3A =VSL_BRNG_R250     +VSL_BRNG_INC)
      PARAMETER (VSL_BRNG_MCG59    =VSL_BRNG_MRG32K3A +VSL_BRNG_INC)
      PARAMETER (VSL_BRNG_WH       =VSL_BRNG_MCG59    +VSL_BRNG_INC)
      PARAMETER (VSL_BRNG_SOBOL    =VSL_BRNG_WH       +VSL_BRNG_INC)
      PARAMETER (VSL_BRNG_NIEDERR  =VSL_BRNG_SOBOL    +VSL_BRNG_INC)
      PARAMETER (VSL_BRNG_MT19937  =VSL_BRNG_NIEDERR  +VSL_BRNG_INC)
      PARAMETER (VSL_BRNG_MT2203   =VSL_BRNG_MT19937  +VSL_BRNG_INC)
      PARAMETER (VSL_BRNG_IABSTRACT=VSL_BRNG_MT2203   +VSL_BRNG_INC)
      PARAMETER (VSL_BRNG_DABSTRACT=VSL_BRNG_IABSTRACT+VSL_BRNG_INC)
      PARAMETER (VSL_BRNG_SABSTRACT=VSL_BRNG_DABSTRACT+VSL_BRNG_INC)

!!  LEAPFROG METHOD FOR GRAY-CODE BASED QUASI-RANDOM NUMBER BASIC GENERATORS
!!  VSL_BRNG_SOBOL and VSL_BRNG_NIEDERR are Gray-code based quasi-random number
!!  basic generators. In contrast to pseudorandom number basic generators,
!!  quasi-random ones take the dimension as initialization parameter.
!!
!!  Suppose that quasi-random number generator (QRNG) dimension is S. QRNG
!!  sequence is a sequence of S-dimensional vectors:
!!
!!     x0=(x0[0],x0[1],...,x0[S-1]),x1=(x1[0],x1[1],...,x1[S-1]),...
!!
!!  VSL treats the output of any basic generator as 1-dimensional, however:
!!
!!     x0[0],x0[1],...,x0[S-1],x1[0],x1[1],...,x1[S-1],...
!!
!!  Because of nature of VSL_BRNG_SOBOL and VSL_BRNG_NIEDERR QRNGs,
!!  the only S-stride Leapfrog method is supported for them. In other words,
!!  user can generate subsequences, which consist of fixed elements of
!!  vectors x0,x1,... For example, if 0 element is fixed, the following
!!  subsequence is generated:
!!
!!     x0[1],x1[1],x2[1],...
!!
!!  To use the s-stride Leapfrog method with given QRNG, user should call
!!  vslLeapfrogStream function with parameter k equal to element to be fixed
!!  (0<=k<S) and parameter nstreams equal to VSL_QRNG_LEAPFROG_COMPONENTS.
      INTEGER VSL_QRNG_LEAPFROG_COMPONENTS
      PARAMETER (VSL_QRNG_LEAPFROG_COMPONENTS = Z"7FFFFFFF")


!!  INITIALIZATION METHODS FOR USER-DESIGNED BASIC RANDOM NUMBER GENERATORS.
!!  Each BRNG must support at least VSL_INIT_METHOD_STANDARD initialization
!!  method. In addition, VSL_INIT_METHOD_LEAPFROG and VSL_INIT_METHOD_SKIPAHEAD
!!  initialization methods can be supported.
!!
!!  If VSL_INIT_METHOD_LEAPFROG is not supported then initialization routine
!!  must return VSL_ERROR_LEAPFROG_UNSUPPORTED error code.
!!
!!  If VSL_INIT_METHOD_SKIPAHEAD is not supported then initialization routine
!!  must return VSL_ERROR_SKIPAHEAD_UNSUPPORTED error code.
!!
!!  If there is no error during initialization, the initialization routine must
!!  return VSL_ERROR_OK code.
      INTEGER VSL_INIT_METHOD_STANDARD
      INTEGER VSL_INIT_METHOD_LEAPFROG
      INTEGER VSL_INIT_METHOD_SKIPAHEAD
      PARAMETER (VSL_INIT_METHOD_STANDARD  = 0)
      PARAMETER (VSL_INIT_METHOD_LEAPFROG  = 1)
      PARAMETER (VSL_INIT_METHOD_SKIPAHEAD = 2)

!!++
!!  TRANSFORMATION METHOD NAMES FOR DISTRIBUTION RANDOM NUMBER GENERATORS
!!  VSL interface allows more than one generation method in a distribution
!!  transformation subroutine. Following macro definitions are used to
!!  specify generation method for given distribution generator.
!!
!!  Method name macro is constructed as
!!
!!     VSL_METHOD_<Precision><Distribution>_<Method>
!!
!!  where
!!
!!     <Precision> - S (single precision) or D (double precision)
!!     <Distribution> - probability distribution
!!     <Method> - method name
!!
!!  VSL_METHOD_<Precision><Distribution>_<Method> should be used with
!!  vsl<precision>Rng<Distribution> function only, where
!!
!!     <precision> - s (single) or d (double)
!!     <Distribution> - probability distribution
!!--

!! Uniform
!!
!! <Method>   <Short Description>
!! STD        standard method. Currently there is only one method for this
!!            distribution generator
      INTEGER VSL_METHOD_SUNIFORM_STD
      INTEGER VSL_METHOD_DUNIFORM_STD
      INTEGER VSL_METHOD_IUNIFORM_STD
      PARAMETER (VSL_METHOD_SUNIFORM_STD = 0)
      PARAMETER (VSL_METHOD_DUNIFORM_STD = 0)
      PARAMETER (VSL_METHOD_IUNIFORM_STD = 0)

!! Uniform Bits
!!
!! <Method>   <Short Description>
!! STD        standard method. Currently there is only one method for this
!!            distribution generator
      INTEGER VSL_METHOD_IUNIFORMBITS_STD
      PARAMETER (VSL_METHOD_IUNIFORMBITS_STD = 0)

!! Gaussian
!!
!! <Method>   <Short Description>
!! BOXMULLER  generates normally distributed random number x thru the pair of
!!            uniformly distributed numbers u1 and u2 according to the formula:
!!
!!               x=sqrt(-ln(u1))*sin(2*Pi*u2)
!!
!! BOXMULLER2 generates pair of normally distributed random numbers x1 and x2
!!            thru the pair of uniformly dustributed numbers u1 and u2
!!            according to the formula
!!
!!               x1=sqrt(-ln(u1))*sin(2*Pi*u2)
!!               x2=sqrt(-ln(u1))*cos(2*Pi*u2)
!!
!!            NOTE: implementation correctly works with odd vector lengths
      INTEGER VSL_METHOD_SGAUSSIAN_BOXMULLER
      INTEGER VSL_METHOD_SGAUSSIAN_BOXMULLER2
      INTEGER VSL_METHOD_DGAUSSIAN_BOXMULLER
      INTEGER VSL_METHOD_DGAUSSIAN_BOXMULLER2
      PARAMETER (VSL_METHOD_SGAUSSIAN_BOXMULLER  = 0)
      PARAMETER (VSL_METHOD_SGAUSSIAN_BOXMULLER2 = 1)
      PARAMETER (VSL_METHOD_DGAUSSIAN_BOXMULLER  = 0)
      PARAMETER (VSL_METHOD_DGAUSSIAN_BOXMULLER2 = 1)

!! GaussianMV - multivariate (correlated) normal
!! Multivariate (correlated) normal random number generator is based on
!! uncorrelated Gaussian random number generator (see vslsRngGaussian and
!! vsldRngGaussian functions):
!!
!! <Method>   <Short Description>
!! BOXMULLER  generates normally distributed random number x thru the pair of
!!            uniformly distributed numbers u1 and u2 according to the formula:
!!
!!               x=sqrt(-ln(u1))*sin(2*Pi*u2)
!!
!! BOXMULLER2 generates pair of normally distributed random numbers x1 and x2
!!            thru the pair of uniformly dustributed numbers u1 and u2
!!            according to the formula
!!
!!               x1=sqrt(-ln(u1))*sin(2*Pi*u2)
!!               x2=sqrt(-ln(u1))*cos(2*Pi*u2)
!!
!!            NOTE: implementation correctly works with odd vector lengths
      INTEGER VSL_METHOD_SGAUSSIANMV_BOXMULLER
      INTEGER VSL_METHOD_SGAUSSIANMV_BOXMULLER2
      INTEGER VSL_METHOD_DGAUSSIANMV_BOXMULLER
      INTEGER VSL_METHOD_DGAUSSIANMV_BOXMULLER2
      PARAMETER (VSL_METHOD_SGAUSSIANMV_BOXMULLER  = 0)
      PARAMETER (VSL_METHOD_SGAUSSIANMV_BOXMULLER2 = 1)
      PARAMETER (VSL_METHOD_DGAUSSIANMV_BOXMULLER  = 0)
      PARAMETER (VSL_METHOD_DGAUSSIANMV_BOXMULLER2 = 1)

!! Exponential
!!
!! <Method>   <Short Description>
!! ICDF       inverse cumulative distribution function method
      INTEGER VSL_METHOD_SEXPONENTIAL_ICDF
      INTEGER VSL_METHOD_DEXPONENTIAL_ICDF
      PARAMETER (VSL_METHOD_SEXPONENTIAL_ICDF = 0)
      PARAMETER (VSL_METHOD_DEXPONENTIAL_ICDF = 0)

!! Laplace
!!
!! <Method>   <Short Description>
!! ICDF       inverse cumulative distribution function method
!!
!! ICDF - inverse cumulative distribution function method:
!!
!!           x=+/-ln(u) with probability 1/2,
!!
!!        where
!!
!!           x - random number with Laplace distribution,
!!           u - uniformly distributed random number
      INTEGER VSL_METHOD_SLAPLACE_ICDF
      INTEGER VSL_METHOD_DLAPLACE_ICDF
      PARAMETER (VSL_METHOD_SLAPLACE_ICDF = 0)
      PARAMETER (VSL_METHOD_DLAPLACE_ICDF = 0)

!! Weibull
!!
!! <Method>   <Short Description>
!! ICDF       inverse cumulative distribution function method
      INTEGER VSL_METHOD_SWEIBULL_ICDF
      INTEGER VSL_METHOD_DWEIBULL_ICDF
      PARAMETER (VSL_METHOD_SWEIBULL_ICDF = 0)
      PARAMETER (VSL_METHOD_DWEIBULL_ICDF = 0)

!! Cauchy
!!
!! <Method>   <Short Description>
!! ICDF       inverse cumulative distribution function method
      INTEGER VSL_METHOD_SCAUCHY_ICDF
      INTEGER VSL_METHOD_DCAUCHY_ICDF
      PARAMETER (VSL_METHOD_SCAUCHY_ICDF = 0)
      PARAMETER (VSL_METHOD_DCAUCHY_ICDF = 0)

!! Rayleigh
!!
!! <Method>   <Short Description>
!! ICDF       inverse cumulative distribution function method
      INTEGER VSL_METHOD_SRAYLEIGH_ICDF
      INTEGER VSL_METHOD_DRAYLEIGH_ICDF
      PARAMETER (VSL_METHOD_SRAYLEIGH_ICDF = 0)
      PARAMETER (VSL_METHOD_DRAYLEIGH_ICDF = 0)

!! Lognormal
!!
!! <Method>   <Short Description>
!! ICDF       inverse cumulative distribution function method
      INTEGER VSL_METHOD_SLOGNORMAL_ICDF
      INTEGER VSL_METHOD_DLOGNORMAL_ICDF
      PARAMETER (VSL_METHOD_SLOGNORMAL_ICDF = 0)
      PARAMETER (VSL_METHOD_DLOGNORMAL_ICDF = 0)

!! Gumbel
!!
!! <Method>   <Short Description>
!! ICDF       inverse cumulative distribution function method
      INTEGER VSL_METHOD_SGUMBEL_ICDF
      INTEGER VSL_METHOD_DGUMBEL_ICDF
      PARAMETER (VSL_METHOD_SGUMBEL_ICDF = 0)
      PARAMETER (VSL_METHOD_DGUMBEL_ICDF = 0)

!! Gamma
!!
!! <Method>     <Short Description>
!! GNORM        nonlinear transformation of gaussian numbers
!! alpha>1,     based on acceptance/rejection method with
!!              squeezes
!!
!! alpha>=0.6,  rejection from the Weibull distribution
!! alpha<1
!!
!! alpha<0.6,   transformation of exponential power distribution
!!              (EPD), EPD random numbers are generated using
!!              by means of acceptance/rejection technique
      INTEGER(KIND=4) VSL_METHOD_SGAMMA_GNORM
      INTEGER(KIND=4) VSL_METHOD_DGAMMA_GNORM
      PARAMETER (VSL_METHOD_SGAMMA_GNORM = 0)
      PARAMETER (VSL_METHOD_DGAMMA_GNORM = 0)

!! Beta
!!
!! <Method>     <Short Description>
!! CJA - stands for first letters of Cheng, Johnk, and Atkinson
!! Cheng      - generation of beta random numbers of the second kind
!! min(p,q)>1   based on acceptance/rejection technique and its
!!              transformation to beta random numbers of the first kind;
!!
!! Johnk,     - if q + K*p^2+C<=0, K=0.852..., C=-0.956...
!! Atkinson,    algorithm of Johnk: beta distributed random number
!! max(p,q)<1   is generated as u1^(1/p) / (u1^(1/p)+u2^(1/q)),
!!              if u1^(1/p)+u2^(1/q)<=1;
!!              otherwise switching algorithm of Atkinson:
!!              interval (0,1) is divided into two domains (0,t) and (t,1),
!!              on each interval acceptance/rejection technique with
!!              convenient majorizing function is used;
!!
!! Atkinson   - switching algorithm of Atkinson is used
!! min(p,q)<1   (with another point t, see short description above);
!! max(p,q)>1
!!
!! ICDF       - inverse cumulative distribution function method according
!!              to formulas x=1-u^(1/q) for p = 1, and x = u^(1/p) for q=1,
!!              where x is beta distributed random number,
!!              u - uniformly distributed random number.
!!              for p=q=1 beta distribution reduces to uniform distribution.

      INTEGER(KIND=4) VSL_METHOD_SBETA_CJA
      INTEGER(KIND=4) VSL_METHOD_DBETA_CJA
      PARAMETER (VSL_METHOD_SBETA_CJA = 0)
      PARAMETER (VSL_METHOD_DBETA_CJA = 0)

!! Bernoulli
!!
!! <Method>   <Short Description>
!! ICDF       inverse cumulative distribution function method
      INTEGER VSL_METHOD_IBERNOULLI_ICDF
      PARAMETER (VSL_METHOD_IBERNOULLI_ICDF = 0)

!! Geometric
!!
!! <Method>   <Short Description>
!! ICDF       inverse cumulative distribution function method
      INTEGER VSL_METHOD_IGEOMETRIC_ICDF
      PARAMETER (VSL_METHOD_IGEOMETRIC_ICDF = 0)

!! Binomial
!!
!! <Method>   <Short Description>
!! BTPE       for ntrial*min(p,1-p)>30 acceptance/rejection method with
!!            decomposition onto 4 regions:
!!
!!               * 2 parallelograms;
!!               * triangle;
!!               * left exponential tail;
!!               * right exponential tail.
!!
!!            othewise table lookup method is used
      INTEGER VSL_METHOD_IBINOMIAL_BTPE
      PARAMETER (VSL_METHOD_IBINOMIAL_BTPE = 0)

!! Hypergeometric
!!
!! <Method>   <Short Description>
!! H2PE       if mode of distribution is large, acceptance/rejection method is
!!            used with decomposition onto 3 regions:
!!
!!               * rectangular;
!!               * left exponential tail;
!!               * right exponential tail.
!!
!!            othewise table lookup method is used
      INTEGER VSL_METHOD_IHYPERGEOMETRIC_H2PE
      PARAMETER (VSL_METHOD_IHYPERGEOMETRIC_H2PE = 0)

!! Poisson
!!
!! <Method>   <Short Description>
!! PTPE       if lambda>=27, acceptance/rejection method is used with
!!            decomposition onto 4 regions:
!!
!!               * 2 parallelograms;
!!               * triangle;
!!               * left exponential tail;
!!               * right exponential tail.
!!
!!            othewise table lookup method is used
!!
!! POISNORM   for lambda>=1 method is based on Poisson inverse CDF
!!            approximation by Gaussian inverse CDF; for lambda<1
!!            table lookup method is used.
      INTEGER VSL_METHOD_IPOISSON_PTPE
      INTEGER VSL_METHOD_IPOISSON_POISNORM
      PARAMETER (VSL_METHOD_IPOISSON_PTPE     = 0)
      PARAMETER (VSL_METHOD_IPOISSON_POISNORM = 1)

!! Poisson
!!
!! <Method>   <Short Description>
!! POISNORM   for lambda>=1 method is based on Poisson inverse CDF
!!            approximation by Gaussian inverse CDF; for lambda<1
!!            ICDF method is used.
      INTEGER VSL_METHOD_IPOISSONV_POISNORM
      PARAMETER (VSL_METHOD_IPOISSONV_POISNORM = 0)

!! Negbinomial
!!
!! <Method>   <Short Description>
!! NBAR       if (a-1)*(1-p)/p>=100, acceptance/rejection method is used with
!!            decomposition onto 5 regions:
!!
!!               * rectangular;
!!               * 2 trapezoid;
!!               * left exponential tail;
!!               * right exponential tail.
!!
!!            othewise table lookup method is used.
      INTEGER VSL_METHOD_INEGBINOMIAL_NBAR
      PARAMETER (VSL_METHOD_INEGBINOMIAL_NBAR = 0)

!!++
!!  MATRIX STORAGE SCHEMES
!!--

!! Some multivariate random number generators, e.g. GaussianMV, operate
!! with matrix parameters. To optimize matrix parameters usage VSL offers
!! following matrix storage schemes. (See VSL documentation for more details).
!!
!! FULL     - whole matrix is stored
!! PACKED   - lower/higher triangular matrix is packed in 1-dimensional array
!! DIAGONAL - diagonal elements are packed in 1-dimensional array
      INTEGER VSL_MATRIX_STORAGE_FULL
      INTEGER VSL_MATRIX_STORAGE_PACKED
      INTEGER VSL_MATRIX_STORAGE_DIAGONAL
      PARAMETER (VSL_MATRIX_STORAGE_FULL     = 0)
      PARAMETER (VSL_MATRIX_STORAGE_PACKED   = 1)
      PARAMETER (VSL_MATRIX_STORAGE_DIAGONAL = 2)

!!++
!!  TYPEDEFS
!!--

!!  VSL STREAM STATE POINTER
!!  This structure is to store VSL stream state address allocated by
!!  VSLNEWSTREAM subroutine.
      TYPE VSL_STREAM_STATE
          INTEGER*4 descriptor1
          INTEGER*4 descriptor2
      END TYPE VSL_STREAM_STATE

      TYPE VSL_CONV_TASK
          INTEGER*4 descriptor1
          INTEGER*4 descriptor2
      END TYPE VSL_CONV_TASK

      TYPE VSL_CORR_TASK
          INTEGER*4 descriptor1
          INTEGER*4 descriptor2
      END TYPE VSL_CORR_TASK

!!  BASIC RANDOM NUMBER GENERATOR PROPERTIES STRUCTURE
!!  The structure describes the properties of given basic generator, e.g. size
!!  of the stream state structure, pointers to function implementations, etc.
!!
!!  BRNG properties structure fields:
!!  StreamStateSize - size of the stream state structure (in bytes)
!!  WordSize        - size of base word (in bytes). Typically this is 4 bytes.
!!  NSeeds          - number of words necessary to describe generator's state
!!  NBits           - number of bits actually used in base word. For example,
!!                    only 31 least significant bits are actually used in
!!                    basic random number generator MCG31m1 with 4-byte base
!!                    word. NBits field is useful while interpreting random
!!                    words as a sequence of random bits.
!!  IncludesZero    - FALSE if 0 cannot be generated in integer-valued
!!                    implementation; TRUE if 0 can be potentially generated in
!!                    integer-valued implementation.
!!  InitStream      - pointer to stream state initialization function
!!  sBRng           - pointer to single precision implementation
!!  dBRng           - pointer to double precision implementation
!!  iBRng           - pointer to integer-value implementation
      TYPE VSL_BRNG_PROPERTIES
          INTEGER streamstatesize
          INTEGER nseeds
          INTEGER includeszero
          INTEGER wordsize
          INTEGER nbits
          INTEGER initstream
          INTEGER sbrng
          INTEGER dbrng
          INTEGER ibrng
      END TYPE VSL_BRNG_PROPERTIES

      END MODULE MKL_VSL_TYPE

      MODULE MKL_VSL

      USE MKL_VSL_TYPE


!!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!!==============================================================================
!!------------------------------------------------------------------------------

      INTERFACE
        INTEGER FUNCTION vsldconvnewtask( task, mode, dims, xshape,
     &                                    yshape, zshape )
            USE MKL_VSL_TYPE
          TYPE(VSL_CONV_TASK) :: task
          INTEGER             :: mode,dims,xshape(dims),yshape(dims)
          INTEGER             :: zshape(dims)
        END FUNCTION
      END INTERFACE

      INTERFACE
        INTEGER FUNCTION vslsconvnewtask( task, mode, dims, xshape,
     &                                    yshape, zshape )
            USE MKL_VSL_TYPE
          TYPE(VSL_CONV_TASK) :: task
          INTEGER             :: mode,dims,xshape(dims),yshape(dims)
          INTEGER             :: zshape(dims)
        END FUNCTION
      END INTERFACE

      INTERFACE
        INTEGER FUNCTION vsldcorrnewtask( task, mode, dims, xshape,
     &                                    yshape, zshape )
            USE MKL_VSL_TYPE
          TYPE(VSL_CORR_TASK) :: task
          INTEGER             :: mode,dims,xshape(dims),yshape(dims)
          INTEGER             :: zshape(dims)
        END FUNCTION
      END INTERFACE

      INTERFACE
        INTEGER FUNCTION vslscorrnewtask( task, mode, dims, xshape,
     &                                    yshape, zshape )
            USE MKL_VSL_TYPE
          TYPE(VSL_CORR_TASK) :: task
          INTEGER             :: mode,dims,xshape(dims),yshape(dims)
          INTEGER             :: zshape(dims)
        END FUNCTION
      END INTERFACE

      INTERFACE
        INTEGER FUNCTION vsldconvnewtask1d( task, mode, xshape, yshape,
     &                                      zshape )
            USE MKL_VSL_TYPE
          TYPE(VSL_CONV_TASK) :: task
          INTEGER             :: mode,xshape,yshape,zshape
        END FUNCTION
      END INTERFACE

      INTERFACE
        INTEGER FUNCTION vslsconvnewtask1d( task, mode, xshape, yshape,
     &                                      zshape )
            USE MKL_VSL_TYPE
          TYPE(VSL_CONV_TASK) :: task
          INTEGER             :: mode,xshape,yshape,zshape
        END FUNCTION
      END INTERFACE

      INTERFACE
        INTEGER FUNCTION vsldcorrnewtask1d( task, mode, xshape, yshape,
     &                                      zshape )
            USE MKL_VSL_TYPE
          TYPE(VSL_CORR_TASK) :: task
          INTEGER             :: mode,xshape,yshape,zshape
        END FUNCTION
      END INTERFACE

      INTERFACE
        INTEGER FUNCTION vslscorrnewtask1d( task, mode, xshape, yshape,
     &                                      zshape )
            USE MKL_VSL_TYPE
          TYPE(VSL_CORR_TASK) :: task
          INTEGER             :: mode,xshape,yshape,zshape
        END FUNCTION
      END INTERFACE

      INTERFACE
        INTEGER FUNCTION vsldconvnewtaskx( task, mode, dims, xshape,
     &                                     yshape, zshape, x, xstride )
            USE MKL_VSL_TYPE
          TYPE(VSL_CONV_TASK) :: task
          INTEGER             :: mode,dims
          INTEGER             :: xshape(dims),yshape(dims)
          INTEGER             :: zshape(dims),xstride(dims)
          REAL(8),DIMENSION(*):: x
        END FUNCTION
      END INTERFACE

      INTERFACE
        INTEGER FUNCTION vslsconvnewtaskx( task, mode, dims, xshape,
     &                                     yshape, zshape, x, xstride )
            USE MKL_VSL_TYPE
          TYPE(VSL_CONV_TASK) :: task
          INTEGER             :: mode,dims
          INTEGER             :: xshape(dims),yshape(dims)
          INTEGER             :: zshape(dims),xstride(dims)
          REAL(4),DIMENSION(*):: x
        END FUNCTION
      END INTERFACE

      INTERFACE
        INTEGER FUNCTION vsldcorrnewtaskx( task, mode, dims, xshape,
     &                                     yshape, zshape, x, xstride )
            USE MKL_VSL_TYPE
          TYPE(VSL_CORR_TASK) :: task
          INTEGER             :: mode,dims
          INTEGER             :: xshape(dims),yshape(dims)
          INTEGER             :: zshape(dims),xstride(dims)
          REAL(8),DIMENSION(*):: x
        END FUNCTION
      END INTERFACE

      INTERFACE
        INTEGER FUNCTION vslscorrnewtaskx( task, mode, dims, xshape,
     &                                     yshape, zshape, x, xstride )
            USE MKL_VSL_TYPE
          TYPE(VSL_CORR_TASK) :: task
          INTEGER             :: mode,dims
          INTEGER             :: xshape(dims),yshape(dims)
          INTEGER             :: zshape(dims),xstride(dims)
          REAL(4),DIMENSION(*):: x
        END FUNCTION
      END INTERFACE

      INTERFACE
        INTEGER FUNCTION vsldconvnewtaskx1d( task, mode, xshape,
     &                                       yshape, zshape, x,xstride)
            USE MKL_VSL_TYPE
          TYPE(VSL_CONV_TASK) :: task
          INTEGER             :: mode,xshape,yshape,zshape,xstride
          REAL(8),DIMENSION(*):: x
        END FUNCTION
      END INTERFACE

      INTERFACE
        INTEGER FUNCTION vslsconvnewtaskx1d( task, mode, xshape,
     &                                       yshape, zshape, x,xstride)
            USE MKL_VSL_TYPE
          TYPE(VSL_CONV_TASK) :: task
          INTEGER             :: mode,xshape,yshape,zshape,xstride
          REAL(4),DIMENSION(*):: x
        END FUNCTION
      END INTERFACE

      INTERFACE
        INTEGER FUNCTION vsldcorrnewtaskx1d( task, mode, xshape,
     &                                       yshape, zshape, x,xstride)
            USE MKL_VSL_TYPE
          TYPE(VSL_CORR_TASK) :: task
          INTEGER             :: mode,xshape,yshape,zshape,xstride
          REAL(8),DIMENSION(*):: x
        END FUNCTION
      END INTERFACE

      INTERFACE
        INTEGER FUNCTION vslscorrnewtaskx1d( task, mode, xshape,
     &                                       yshape, zshape, x,xstride)
            USE MKL_VSL_TYPE
          TYPE(VSL_CORR_TASK) :: task
          INTEGER             :: mode,xshape,yshape,zshape,xstride
          REAL(4),DIMENSION(*):: x
        END FUNCTION
      END INTERFACE

      INTERFACE
        INTEGER FUNCTION vslconvdeletetask( task )
            USE MKL_VSL_TYPE
          TYPE(VSL_CONV_TASK) :: task
        END FUNCTION
      END INTERFACE

      INTERFACE
        INTEGER FUNCTION vslcorrdeletetask( task )
            USE MKL_VSL_TYPE
          TYPE(VSL_CORR_TASK) :: task
        END FUNCTION
      END INTERFACE

      INTERFACE
        INTEGER FUNCTION vslconvcopytask( desttask, srctask )
            USE MKL_VSL_TYPE
          TYPE(VSL_CONV_TASK) :: desttask
          TYPE(VSL_CONV_TASK) :: srctask
        END FUNCTION
      END INTERFACE

      INTERFACE
        INTEGER FUNCTION vslcorrcopytask( desttask, srctask )
            USE MKL_VSL_TYPE
          TYPE(VSL_CORR_TASK) :: desttask
          TYPE(VSL_CORR_TASK) :: srctask
        END FUNCTION
      END INTERFACE

      INTERFACE
        INTEGER FUNCTION vslConvSetMode( task, newmode )
            USE MKL_VSL_TYPE
          TYPE(VSL_CONV_TASK) :: task
          INTEGER             :: newmode
        END FUNCTION
      END INTERFACE

      INTERFACE
        INTEGER FUNCTION vslCorrSetMode( task, newmode )
            USE MKL_VSL_TYPE
          TYPE(VSL_CORR_TASK) :: task
          INTEGER             :: newmode
        END FUNCTION
      END INTERFACE

      INTERFACE
        INTEGER FUNCTION vslConvSetInternalPrecision( task, precision )
            USE MKL_VSL_TYPE
          TYPE(VSL_CONV_TASK) :: task
          INTEGER             :: precision
        END FUNCTION
      END INTERFACE

      INTERFACE
        INTEGER FUNCTION vslCorrSetInternalPrecision( task, precision )
            USE MKL_VSL_TYPE
          TYPE(VSL_CORR_TASK) :: task
          INTEGER             :: precision
        END FUNCTION
      END INTERFACE

      INTERFACE
        INTEGER FUNCTION vslConvSetStart( task, start )
            USE MKL_VSL_TYPE
          TYPE(VSL_CONV_TASK) :: task
          INTEGER,DIMENSION(*):: start
        END FUNCTION
      END INTERFACE

      INTERFACE
        INTEGER FUNCTION vslCorrSetStart( task, start )
            USE MKL_VSL_TYPE
          TYPE(VSL_CORR_TASK) :: task
          INTEGER,DIMENSION(*):: start
        END FUNCTION
      END INTERFACE

      INTERFACE
        INTEGER FUNCTION vslConvSetDecimation( task, decimation )
            USE MKL_VSL_TYPE
          TYPE(VSL_CONV_TASK) :: task
          INTEGER,DIMENSION(*):: decimation
        END FUNCTION
      END INTERFACE

      INTERFACE
        INTEGER FUNCTION vslCorrSetDecimation( task, decimation )
            USE MKL_VSL_TYPE
          TYPE(VSL_CORR_TASK) :: task
          INTEGER,DIMENSION(*):: decimation
        END FUNCTION
      END INTERFACE


      INTERFACE
        INTEGER FUNCTION vsldconvexec( task, x, xstride, y, ystride, z,
     &                                 zstride )
            USE MKL_VSL_TYPE
          TYPE(VSL_CONV_TASK)      :: task
          INTEGER,DIMENSION(*)     :: xstride,ystride,zstride
          REAL(KIND=8),DIMENSION(*):: x,y,z
        END FUNCTION
      END INTERFACE

      INTERFACE
        INTEGER FUNCTION vslsconvexec( task, x, xstride, y, ystride, z,
     &                                 zstride )
            USE MKL_VSL_TYPE
          TYPE(VSL_CONV_TASK)      :: task
          INTEGER,DIMENSION(*)     :: xstride,ystride,zstride
          REAL(KIND=4),DIMENSION(*):: x,y,z
        END FUNCTION
      END INTERFACE

      INTERFACE
        INTEGER FUNCTION vsldcorrexec( task, x, xstride, y, ystride, z,
     &                                 zstride )
            USE MKL_VSL_TYPE
          TYPE(VSL_CORR_TASK)      :: task
          INTEGER,DIMENSION(*)     :: xstride,ystride,zstride
          REAL(KIND=8),DIMENSION(*):: x,y,z
        END FUNCTION
      END INTERFACE

      INTERFACE
        INTEGER FUNCTION vslscorrexec( task, x, xstride, y, ystride, z,
     &                                 zstride )
            USE MKL_VSL_TYPE
          TYPE(VSL_CORR_TASK)      :: task
          INTEGER,DIMENSION(*)     :: xstride,ystride,zstride
          REAL(KIND=4),DIMENSION(*):: x,y,z
        END FUNCTION
      END INTERFACE

      INTERFACE
        INTEGER FUNCTION vsldconvexec1d( task, x, xstride, y, ystride,
     &                                   z, zstride )
            USE MKL_VSL_TYPE
          TYPE(VSL_CONV_TASK)      :: task
          INTEGER                  :: xstride,ystride,zstride
          REAL(KIND=8),DIMENSION(*):: x,y,z
        END FUNCTION
      END INTERFACE

      INTERFACE
        INTEGER FUNCTION vslsconvexec1d( task, x, xstride, y, ystride,
     &                                   z, zstride )
            USE MKL_VSL_TYPE
          TYPE(VSL_CONV_TASK)      :: task
          INTEGER                  :: xstride,ystride,zstride
          REAL(KIND=4),DIMENSION(*):: x,y,z
        END FUNCTION
      END INTERFACE

      INTERFACE
        INTEGER FUNCTION vsldcorrexec1d( task, x, xstride, y, ystride,
     &                                   z, zstride )
            USE MKL_VSL_TYPE
          TYPE(VSL_CORR_TASK)      :: task
          INTEGER                  :: xstride,ystride,zstride
          REAL(KIND=8),DIMENSION(*):: x,y,z
        END FUNCTION
      END INTERFACE

      INTERFACE
        INTEGER FUNCTION vslscorrexec1d( task, x, xstride, y, ystride,
     &                                   z, zstride )
            USE MKL_VSL_TYPE
          TYPE(VSL_CORR_TASK)      :: task
          INTEGER                  :: xstride,ystride,zstride
          REAL(KIND=4),DIMENSION(*):: x,y,z
        END FUNCTION
      END INTERFACE

      INTERFACE
        INTEGER FUNCTION vsldconvexecx( task, y, ystride, z, zstride )
            USE MKL_VSL_TYPE
          TYPE(VSL_CONV_TASK)      :: task
          INTEGER,DIMENSION(*)     :: ystride,zstride
          REAL(KIND=8),DIMENSION(*):: y,z
        END FUNCTION
      END INTERFACE

      INTERFACE
        INTEGER FUNCTION vslsconvexecx( task, y, ystride, z, zstride )
            USE MKL_VSL_TYPE
          TYPE(VSL_CONV_TASK)      :: task
          INTEGER,DIMENSION(*)     :: ystride,zstride
          REAL(KIND=4),DIMENSION(*):: y,z
        END FUNCTION
      END INTERFACE

      INTERFACE
        INTEGER FUNCTION vsldcorrexecx( task, y, ystride, z, zstride )
            USE MKL_VSL_TYPE
          TYPE(VSL_CORR_TASK)      :: task
          INTEGER,DIMENSION(*)     :: ystride,zstride
          REAL(KIND=8),DIMENSION(*):: y,z
        END FUNCTION
      END INTERFACE

      INTERFACE
        INTEGER FUNCTION vslscorrexecx( task, y, ystride, z, zstride )
            USE MKL_VSL_TYPE
          TYPE(VSL_CORR_TASK)      :: task
          INTEGER,DIMENSION(*)     :: ystride,zstride
          REAL(KIND=4),DIMENSION(*):: y,z
        END FUNCTION
      END INTERFACE

      INTERFACE
        INTEGER FUNCTION vsldconvexecx1d( task, y, ystride, z,z stride)
            USE MKL_VSL_TYPE
          TYPE(VSL_CONV_TASK)      :: task
          INTEGER                  :: ystride,zstride
          REAL(KIND=8),DIMENSION(*):: y,z
        END FUNCTION
      END INTERFACE

      INTERFACE
        INTEGER FUNCTION vslsconvexecx1d( task, y, ystride, z, zstride)
            USE MKL_VSL_TYPE
          TYPE(VSL_CONV_TASK)      :: task
          INTEGER                  :: ystride,zstride
          REAL(KIND=4),DIMENSION(*):: y,z
        END FUNCTION
      END INTERFACE

      INTERFACE
        INTEGER FUNCTION vsldcorrexecx1d( task, y, ystride, z, zstride)
            USE MKL_VSL_TYPE
          TYPE(VSL_CORR_TASK)      :: task
          INTEGER                  :: ystride,zstride
          REAL(KIND=8),DIMENSION(*):: y,z
        END FUNCTION
      END INTERFACE

      INTERFACE
        INTEGER FUNCTION vslscorrexecx1d( task, y, ystride, z, zstride)
            USE MKL_VSL_TYPE
          TYPE(VSL_CORR_TASK)      :: task
          INTEGER                  :: ystride,zstride
          REAL(KIND=4),DIMENSION(*):: y,z
        END FUNCTION
      END INTERFACE

!!++
!!  VSL CONTINUOUS DISTRIBUTION GENERATOR FUNCTION INTERFACES.
!!--

!!  Uniform distribution
      INTERFACE
        INTEGER FUNCTION vsrnguniform( method, stream, n, r, a, b )
            USE MKL_VSL_TYPE
          INTEGER,INTENT(IN)       :: method
          TYPE(VSL_STREAM_STATE)   :: stream
          INTEGER,INTENT(IN)       :: n
          REAL(KIND=4),INTENT(OUT) :: r(n)
          REAL(KIND=4),INTENT(IN)  :: a
          REAL(KIND=4),INTENT(IN)  :: b
        END FUNCTION
      END INTERFACE

      INTERFACE
        INTEGER FUNCTION vdrnguniform( method, stream, n, r, a, b )
            USE MKL_VSL_TYPE
          INTEGER,INTENT(IN)       :: method
          TYPE(VSL_STREAM_STATE)   :: stream
          INTEGER,INTENT(IN)       :: n
          REAL(KIND=8),INTENT(OUT) :: r(n)
          REAL(KIND=8),INTENT(IN)  :: a
          REAL(KIND=8),INTENT(IN)  :: b
        END FUNCTION
      END INTERFACE

!!  Gaussian distribution
      INTERFACE
        INTEGER FUNCTION vsrnggaussian( method, stream, n, r, a, sigma)
            USE MKL_VSL_TYPE
          INTEGER,INTENT(IN)       :: method
          TYPE(VSL_STREAM_STATE)   :: stream
          INTEGER,INTENT(IN)       :: n
          REAL(KIND=4),INTENT(OUT) :: r(n)
          REAL(KIND=4),INTENT(IN)  :: a
          REAL(KIND=4),INTENT(IN)  :: sigma
        END FUNCTION
      END INTERFACE

      INTERFACE
        INTEGER FUNCTION vdrnggaussian( method, stream, n, r, a, sigma)
            USE MKL_VSL_TYPE
          INTEGER,INTENT(IN)       :: method
          TYPE(VSL_STREAM_STATE)   :: stream
          INTEGER,INTENT(IN)       :: n
          REAL(KIND=8),INTENT(OUT) :: r(n)
          REAL(KIND=8),INTENT(IN)  :: a
          REAL(KIND=8),INTENT(IN)  :: sigma
        END FUNCTION
      END INTERFACE

!!  GaussianMV distribution
      INTERFACE
        INTEGER FUNCTION vsrnggaussianmv( method, stream, n, r, dimen,
     &                              mstorage, a, t )
            USE MKL_VSL_TYPE
          INTEGER,INTENT(IN)       :: method
          TYPE(VSL_STREAM_STATE)   :: stream
          INTEGER,INTENT(IN)       :: n
          INTEGER,INTENT(IN)       :: dimen
          REAL(KIND=4),INTENT(OUT) :: r(dimen,n)
          INTEGER,INTENT(IN)       :: mstorage
          REAL(KIND=4),INTENT(IN)  :: a(dimen)
          REAL(KIND=4),INTENT(IN)  :: t(dimen,dimen)
        END FUNCTION
      END INTERFACE

      INTERFACE
        INTEGER FUNCTION vdrnggaussianmv( method, stream, n, r, dimen,
     &                              mstorage, a, t )
            USE MKL_VSL_TYPE
          INTEGER,INTENT(IN)       :: method
          TYPE(VSL_STREAM_STATE)   :: stream
          INTEGER,INTENT(IN)       :: n
          REAL(KIND=8),INTENT(OUT) :: r(n)
          INTEGER,INTENT(IN)       :: dimen
          INTEGER,INTENT(IN)       :: mstorage
          REAL(KIND=8),INTENT(IN)  :: a(dimen)
          REAL(KIND=8),INTENT(IN)  :: t(dimen,dimen)
        END FUNCTION
      END INTERFACE

!!  Exponential distribution
      INTERFACE
        INTEGER FUNCTION vsrngexponential( method, stream, n, r, a,
     &                                     beta )
            USE MKL_VSL_TYPE
          INTEGER,INTENT(IN)       :: method
          TYPE(VSL_STREAM_STATE)   :: stream
          INTEGER,INTENT(IN)       :: n
          REAL(KIND=4),INTENT(OUT) :: r(n)
          REAL(KIND=4),INTENT(IN)  :: a
          REAL(KIND=4),INTENT(IN)  :: beta
        END FUNCTION
      END INTERFACE

      INTERFACE
        INTEGER FUNCTION vdrngexponential( method, stream, n, r, a,
     &                                     beta )
            USE MKL_VSL_TYPE
          INTEGER,INTENT(IN)       :: method
          TYPE(VSL_STREAM_STATE)   :: stream
          INTEGER,INTENT(IN)       :: n
          REAL(KIND=8),INTENT(OUT) :: r(n)
          REAL(KIND=8),INTENT(IN)  :: a
          REAL(KIND=8),INTENT(IN)  :: beta
        END FUNCTION
      END INTERFACE

!!  Laplace distribution
      INTERFACE
        INTEGER FUNCTION vsrnglaplace( method, stream, n, r, a, beta )
            USE MKL_VSL_TYPE
          INTEGER,INTENT(IN)       :: method
          TYPE(VSL_STREAM_STATE)   :: stream
          INTEGER,INTENT(IN)       :: n
          REAL(KIND=4),INTENT(OUT) :: r(n)
          REAL(KIND=4),INTENT(IN)  :: a
          REAL(KIND=4),INTENT(IN)  :: beta
        END FUNCTION
      END INTERFACE

      INTERFACE
        INTEGER FUNCTION vdrnglaplace( method, stream, n, r, a, beta )
            USE MKL_VSL_TYPE
          INTEGER,INTENT(IN)       :: method
          TYPE(VSL_STREAM_STATE)   :: stream
          INTEGER,INTENT(IN)       :: n
          REAL(KIND=8),INTENT(OUT) :: r(n)
          REAL(KIND=8),INTENT(IN)  :: a
          REAL(KIND=8),INTENT(IN)  :: beta
        END FUNCTION
      END INTERFACE

!!  Weibull distribution
      INTERFACE
        INTEGER FUNCTION vsrngweibull( method, stream, n, r, alpha, a,
     &                                 beta )
            USE MKL_VSL_TYPE
          INTEGER,INTENT(IN)       :: method
          TYPE(VSL_STREAM_STATE)   :: stream
          INTEGER,INTENT(IN)       :: n
          REAL(KIND=4),INTENT(OUT) :: r(n)
          REAL(KIND=4),INTENT(IN)  :: alpha
          REAL(KIND=4),INTENT(IN)  :: a
          REAL(KIND=4),INTENT(IN)  :: beta
        END FUNCTION
      END INTERFACE

      INTERFACE
        INTEGER FUNCTION vdrngweibull( method, stream, n, r, alpha, a,
     &                                 beta )
            USE MKL_VSL_TYPE
          INTEGER,INTENT(IN)       :: method
          TYPE(VSL_STREAM_STATE)   :: stream
          INTEGER,INTENT(IN)       :: n
          REAL(KIND=8),INTENT(OUT) :: r(n)
          REAL(KIND=8),INTENT(IN)  :: alpha
          REAL(KIND=8),INTENT(IN)  :: a
          REAL(KIND=8),INTENT(IN)  :: beta
        END FUNCTION
      END INTERFACE

!!  Cauchy distribution
      INTERFACE
        INTEGER FUNCTION vsrngcauchy( method, stream, n, r, a, beta )
            USE MKL_VSL_TYPE
          INTEGER,INTENT(IN)       :: method
          TYPE(VSL_STREAM_STATE)   :: stream
          INTEGER,INTENT(IN)       :: n
          REAL(KIND=4),INTENT(OUT) :: r(n)
          REAL(KIND=4),INTENT(IN)  :: a
          REAL(KIND=4),INTENT(IN)  :: beta
        END FUNCTION
      END INTERFACE

      INTERFACE
        INTEGER FUNCTION vdrngcauchy( method, stream, n, r, a, beta )
            USE MKL_VSL_TYPE
          INTEGER,INTENT(IN)       :: method
          TYPE(VSL_STREAM_STATE)   :: stream
          INTEGER,INTENT(IN)       :: n
          REAL(KIND=8),INTENT(OUT) :: r(n)
          REAL(KIND=8),INTENT(IN)  :: a
          REAL(KIND=8),INTENT(IN)  :: beta
        END FUNCTION
      END INTERFACE

!!  Rayleigh distribution
      INTERFACE
        INTEGER FUNCTION vsrngrayleigh( method, stream, n, r, a, beta )
            USE MKL_VSL_TYPE
          INTEGER,INTENT(IN)       :: method
          TYPE(VSL_STREAM_STATE)   :: stream
          INTEGER,INTENT(IN)       :: n
          REAL(KIND=4),INTENT(OUT) :: r(n)
          REAL(KIND=4),INTENT(IN)  :: a
          REAL(KIND=4),INTENT(IN)  :: beta
        END FUNCTION
      END INTERFACE

      INTERFACE
        INTEGER FUNCTION vdrngrayleigh( method, stream, n, r, a, beta )
            USE MKL_VSL_TYPE
          INTEGER,INTENT(IN)       :: method
          TYPE(VSL_STREAM_STATE)   :: stream
          INTEGER,INTENT(IN)       :: n
          REAL(KIND=8),INTENT(OUT) :: r(n)
          REAL(KIND=8),INTENT(IN)  :: a
          REAL(KIND=8),INTENT(IN)  :: beta
        END FUNCTION
      END INTERFACE

!!  Lognormal distribution
      INTERFACE
        INTEGER FUNCTION vsrnglognormal( method, stream, n, r, a,sigma,
     &                             b, beta )
            USE MKL_VSL_TYPE
          INTEGER,INTENT(IN)       :: method
          TYPE(VSL_STREAM_STATE)   :: stream
          INTEGER,INTENT(IN)       :: n
          REAL(KIND=4),INTENT(OUT) :: r(n)
          REAL(KIND=4),INTENT(IN)  :: a
          REAL(KIND=4),INTENT(IN)  :: sigma
          REAL(KIND=4),INTENT(IN)  :: b
          REAL(KIND=4),INTENT(IN)  :: beta
        END FUNCTION
      END INTERFACE

      INTERFACE
        INTEGER FUNCTION vdrnglognormal( method, stream, n, r, a,sigma,
     &                             b, beta )
            USE MKL_VSL_TYPE
          INTEGER,INTENT(IN)       :: method
          TYPE(VSL_STREAM_STATE)   :: stream
          INTEGER,INTENT(IN)       :: n
          REAL(KIND=8),INTENT(OUT) :: r(n)
          REAL(KIND=8),INTENT(IN)  :: a
          REAL(KIND=8),INTENT(IN)  :: sigma
          REAL(KIND=8),INTENT(IN)  :: b
          REAL(KIND=8),INTENT(IN)  :: beta
        END FUNCTION
      END INTERFACE

!!  Gumbel distribution
      INTERFACE
        INTEGER FUNCTION vsrnggumbel( method, stream, n, r, a, beta )
            USE MKL_VSL_TYPE
          INTEGER,INTENT(IN)       :: method
          TYPE(VSL_STREAM_STATE)   :: stream
          INTEGER,INTENT(IN)       :: n
          REAL(KIND=4),INTENT(OUT) :: r(n)
          REAL(KIND=4),INTENT(IN)  :: a
          REAL(KIND=4),INTENT(IN)  :: beta
        END FUNCTION
      END INTERFACE

      INTERFACE
        INTEGER FUNCTION vdrnggumbel( method, stream, n, r, a, beta )
            USE MKL_VSL_TYPE
          INTEGER,INTENT(IN)       :: method
          TYPE(VSL_STREAM_STATE)   :: stream
          INTEGER,INTENT(IN)       :: n
          REAL(KIND=8),INTENT(OUT) :: r(n)
          REAL(KIND=8),INTENT(IN)  :: a
          REAL(KIND=8),INTENT(IN)  :: beta
        END FUNCTION
      END INTERFACE

!!  Gamma distribution
      INTERFACE
        INTEGER FUNCTION vsrnggamma( method, stream, n, r, alpha, a,
     &                               beta )
          USE MKL_VSL_TYPE
          INTEGER,INTENT(IN)             :: method
          TYPE(VSL_STREAM_STATE)         :: stream
          INTEGER,INTENT(IN)             :: n
          REAL(KIND=4),INTENT(OUT)       :: r(n)
          REAL(KIND=4),INTENT(IN)        :: alpha
          REAL(KIND=4),INTENT(IN)        :: a
          REAL(KIND=4),INTENT(IN)        :: beta
        END FUNCTION
      END INTERFACE

      INTERFACE
        INTEGER FUNCTION vdrnggamma( method, stream, n, r, alpha, a,
     &                               beta )
          USE MKL_VSL_TYPE
          INTEGER(KIND=4),INTENT(IN)     :: method
          TYPE(VSL_STREAM_STATE)         :: stream
          INTEGER(KIND=4),INTENT(IN)     :: n
          REAL(KIND=8),INTENT(OUT)       :: r(n)
          REAL(KIND=8),INTENT(IN)        :: alpha
          REAL(KIND=8),INTENT(IN)        :: a
          REAL(KIND=8),INTENT(IN)        :: beta
        END FUNCTION
      END INTERFACE

!!  Beta distribution
      INTERFACE
        INTEGER FUNCTION vsrngbeta( method, stream, n, r, p, q, a,
     &                              beta )
          USE MKL_VSL_TYPE
          INTEGER,INTENT(IN)             :: method
          TYPE(VSL_STREAM_STATE)         :: stream
          INTEGER,INTENT(IN)             :: n
          REAL(KIND=4),INTENT(OUT)       :: r(n)
          REAL(KIND=4),INTENT(IN)        :: p
          REAL(KIND=4),INTENT(IN)        :: q
          REAL(KIND=4),INTENT(IN)        :: a
          REAL(KIND=4),INTENT(IN)        :: beta
        END FUNCTION
      END INTERFACE

      INTERFACE
        INTEGER FUNCTION vdrngbeta( method, stream, n, r, p, q, a,
     &                              beta )
          USE MKL_VSL_TYPE
          INTEGER(KIND=4),INTENT(IN)     :: method
          TYPE(VSL_STREAM_STATE)         :: stream
          INTEGER(KIND=4),INTENT(IN)     :: n
          REAL(KIND=8),INTENT(OUT)       :: r(n)
          REAL(KIND=8),INTENT(IN)        :: p
          REAL(KIND=8),INTENT(IN)        :: q
          REAL(KIND=8),INTENT(IN)        :: a
          REAL(KIND=8),INTENT(IN)        :: beta
        END FUNCTION
      END INTERFACE

!!++
!!  VSL DISCRETE DISTRIBUTION GENERATOR FUNCTION INTERFACES.
!!--

!!  Uniform distribution
      INTERFACE
        INTEGER FUNCTION virnguniform( method, stream, n, r, a, b )
            USE MKL_VSL_TYPE
          INTEGER,INTENT(IN)     :: method
          TYPE(VSL_STREAM_STATE) :: stream
          INTEGER,INTENT(IN)     :: n
          INTEGER,INTENT(OUT)    :: r(n)
          INTEGER,INTENT(IN)     :: a
          INTEGER,INTENT(IN)     :: b
        END FUNCTION
      END INTERFACE

!!  UniformBits distribution
      INTERFACE
        INTEGER FUNCTION virnguniformbits( method, stream, n, r )
            USE MKL_VSL_TYPE
          INTEGER,INTENT(IN)     :: method
          TYPE(VSL_STREAM_STATE) :: stream
          INTEGER,INTENT(IN)     :: n
          INTEGER,INTENT(OUT)    :: r(n)
        END FUNCTION
      END INTERFACE

!!  Bernoulli distribution
      INTERFACE
        INTEGER FUNCTION virngbernoulli( method, stream, n, r, p )
            USE MKL_VSL_TYPE
          INTEGER,INTENT(IN)     :: method
          TYPE(VSL_STREAM_STATE) :: stream
          INTEGER,INTENT(IN)     :: n
          INTEGER,INTENT(OUT)    :: r(n)
          REAL(KIND=8),INTENT(IN):: p
        END FUNCTION
      END INTERFACE

!!  Geometric distribution
      INTERFACE
        INTEGER FUNCTION virnggeometric( method, stream, n, r, p )
            USE MKL_VSL_TYPE
          INTEGER,INTENT(IN)     :: method
          TYPE(VSL_STREAM_STATE) :: stream
          INTEGER,INTENT(IN)     :: n
          INTEGER,INTENT(OUT)    :: r(n)
          REAL(KIND=8),INTENT(IN):: p
        END FUNCTION
      END INTERFACE

!!  Binomial distribution
      INTERFACE
        INTEGER FUNCTION virngbinomial(method, stream, n, r, ntrial, p)
            USE MKL_VSL_TYPE
          INTEGER,INTENT(IN)     :: method
          TYPE(VSL_STREAM_STATE) :: stream
          INTEGER,INTENT(IN)     :: n
          INTEGER,INTENT(OUT)    :: r(n)
          INTEGER,INTENT(IN)     :: ntrial
          REAL(KIND=8),INTENT(IN):: p
        END FUNCTION
      END INTERFACE

!!  Hypergeometric distribution
      INTERFACE
        INTEGER FUNCTION virnghypergeometric( method, stream, n, r, l,
     &                                        s, m )
            USE MKL_VSL_TYPE
          INTEGER,INTENT(IN)     :: method
          TYPE(VSL_STREAM_STATE) :: stream
          INTEGER,INTENT(IN)     :: n
          INTEGER,INTENT(OUT)    :: r(n)
          INTEGER,INTENT(IN)     :: l
          INTEGER,INTENT(IN)     :: s
          INTEGER,INTENT(IN)     :: m
        END FUNCTION
      END INTERFACE

!!  Poisson distribution
      INTERFACE
        INTEGER FUNCTION virngpoisson( method, stream, n, r, lambda )
            USE MKL_VSL_TYPE
          INTEGER,INTENT(IN)     :: method
          TYPE(VSL_STREAM_STATE) :: stream
          INTEGER,INTENT(IN)     :: n
          INTEGER,INTENT(OUT)    :: r(n)
          REAL(KIND=8),INTENT(IN):: lambda
        END FUNCTION
      END INTERFACE

!!  PoissonV distribution
      INTERFACE
        INTEGER FUNCTION virngpoissonv( method, stream, n, r, lambda )
            USE MKL_VSL_TYPE
          INTEGER,INTENT(IN)     :: method
          TYPE(VSL_STREAM_STATE) :: stream
          INTEGER,INTENT(IN)     :: n
          INTEGER,INTENT(OUT)    :: r(n)
          REAL(KIND=8),INTENT(IN):: lambda(n)
        END FUNCTION
      END INTERFACE

!!  Negbinomial distribution
      INTERFACE
        INTEGER FUNCTION virngnegbinomial( method, stream, n, r, a, p )
            USE MKL_VSL_TYPE
          INTEGER,INTENT(IN)     :: method
          TYPE(VSL_STREAM_STATE) :: stream
          INTEGER,INTENT(IN)     :: n
          INTEGER,INTENT(OUT)    :: r(n)
          REAL(KIND=8),INTENT(IN):: a
          REAL(KIND=8),INTENT(IN):: p
        END FUNCTION
      END INTERFACE

!!++
!!  VSL SERVICE FUNCTION INTERFACES.
!!--

!! NewStream - stream creation/initialization
      INTERFACE
        INTEGER FUNCTION vslnewstream( stream, brng, seed )
            USE MKL_VSL_TYPE
          TYPE(VSL_STREAM_STATE) :: stream
          INTEGER,INTENT(IN)     :: brng
          INTEGER,INTENT(IN)     :: seed
        END FUNCTION
      END INTERFACE

!! NewStreamEx - advanced stream creation/initialization
      INTERFACE
        INTEGER FUNCTION vslnewstreamex( stream, brng, n, params )
            USE MKL_VSL_TYPE
          TYPE(VSL_STREAM_STATE) :: stream
          INTEGER,INTENT(IN)     :: brng
          INTEGER,INTENT(IN)     :: n
          INTEGER,INTENT(IN)     :: params(n)
        END FUNCTION
      END INTERFACE

!!    INEWABSTRACTSTREAM
      INTERFACE
        INTEGER FUNCTION vslinewabstractstream( stream, n, ibuf, ifunc)
          USE MKL_VSL_TYPE
          TYPE(VSL_STREAM_STATE),INTENT(OUT) :: stream
          INTEGER(KIND=4),INTENT(IN)         :: n
          INTEGER(KIND=4),INTENT(IN)         :: ibuf(n)
          INTEGER(KIND=4),EXTERNAL           :: ifunc
        END FUNCTION
      END INTERFACE

!!    DNEWABSTRACTSTREAM
      INTERFACE
        INTEGER FUNCTION vsldnewabstractstream( stream, n, dbuf, a, b,
     &                                          dfunc )
          USE MKL_VSL_TYPE
          TYPE(VSL_STREAM_STATE),INTENT(OUT) :: stream
          INTEGER(KIND=4),INTENT(IN)         :: n
          REAL(KIND=8)   ,INTENT(IN)         :: dbuf(n)
          REAL(KIND=8)   ,INTENT(IN)         :: a
          REAL(KIND=8)   ,INTENT(IN)         :: b
          INTEGER(KIND=4),EXTERNAL           :: dfunc
        END FUNCTION
      END INTERFACE

!!    SNEWABSTRACTSTREAM
      INTERFACE
        INTEGER FUNCTION vslsnewabstractstream( stream, n, sbuf, a, b,
     &                                          sfunc )
          USE MKL_VSL_TYPE
          TYPE(VSL_STREAM_STATE),INTENT(OUT) :: stream
          INTEGER(KIND=4),INTENT(IN)         :: n
          REAL(KIND=4)   ,INTENT(IN)         :: sbuf(n)
          REAL(KIND=4)   ,INTENT(IN)         :: a
          REAL(KIND=4)   ,INTENT(IN)         :: b
          INTEGER(KIND=4),EXTERNAL           :: sfunc
        END FUNCTION
      END INTERFACE

!! DeleteStream - delete stream
      INTERFACE
        INTEGER FUNCTION vsldeletestream( stream )
            USE MKL_VSL_TYPE
          TYPE(VSL_STREAM_STATE) :: stream
        END FUNCTION
      END INTERFACE

!! CopyStream - copy all stream information
      INTERFACE
        INTEGER FUNCTION vslcopystream( newstream, srcstream )
            USE MKL_VSL_TYPE
          TYPE(VSL_STREAM_STATE) :: newstream
          TYPE(VSL_STREAM_STATE) :: srcstream
        END FUNCTION
      END INTERFACE

!! CopyStreamState - copy stream state only
      INTERFACE
        INTEGER FUNCTION vslcopystreamstate( deststream, srcstream )
            USE MKL_VSL_TYPE
          TYPE(VSL_STREAM_STATE) :: deststream
          TYPE(VSL_STREAM_STATE) :: srcstream
        END FUNCTION
      END INTERFACE

!! LeapfrogStream - leapfrog method
      INTERFACE
        INTEGER FUNCTION vslleapfrogstream( stream, k, nstreams )
            USE MKL_VSL_TYPE
          TYPE(VSL_STREAM_STATE) :: stream
          INTEGER,INTENT(IN)     :: k
          INTEGER,INTENT(IN)     :: nstreams
        END FUNCTION
      END INTERFACE

!! SkipAheadStream - skip-ahead method
      INTERFACE
        INTEGER FUNCTION vslskipaheadstream( stream, nskip )
            USE MKL_VSL_TYPE
          TYPE(VSL_STREAM_STATE) :: stream
          INTEGER,INTENT(IN)     :: nskip
        END FUNCTION
      END INTERFACE

!! GetBrngProperties - get BRNG properties
      INTERFACE
        INTEGER FUNCTION vslgetbrngproperties( brng, properties )
            USE MKL_VSL_TYPE
          INTEGER,INTENT(IN) :: brng
          TYPE(VSL_BRNG_PROPERTIES),INTENT(OUT) :: properties
        END FUNCTION
      END INTERFACE

!! GetNumRegBrngs - get number of registered BRNGs
      INTERFACE
        INTEGER FUNCTION vslgetnumregbrngs( )
        END FUNCTION
      END INTERFACE

!! GetStreamStateBrng - get BRNG associated with given stream
      INTERFACE
        INTEGER FUNCTION vslgetstreamstatebrng( stream )
            USE MKL_VSL_TYPE
          TYPE(VSL_STREAM_STATE) :: stream
        END FUNCTION
      END INTERFACE

!! RegisterBrng - register new BRNG
      INTERFACE
        INTEGER FUNCTION vslregisterbrng( properties )
            USE MKL_VSL_TYPE
          TYPE(VSL_BRNG_PROPERTIES)      :: properties
        END FUNCTION
      END INTERFACE

!! SaveStreamF - save stream to file
      INTERFACE
        INTEGER FUNCTION vslsavestreamf( stream, fname )
            USE MKL_VSL_TYPE
          CHARACTER(*)           :: fname
          TYPE(VSL_STREAM_STATE) :: stream
        END FUNCTION
      END INTERFACE

!! LoadStreamF - save stream to file
      INTERFACE
        INTEGER FUNCTION vslloadstreamf( stream, fname )
            USE MKL_VSL_TYPE
          CHARACTER(*)           :: fname
          TYPE(VSL_STREAM_STATE) :: stream
        END FUNCTION
      END INTERFACE

      END MODULE MKL_VSL

!--------------------------------------------------------------------!