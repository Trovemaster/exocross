module accuracy
  implicit none
  private
  public sik, ik, hik, rk, ark, out, inp, safe_max,safe_min,max_exp, pi, twopi, cl, wl, c2
  public accuracyInitialize
  public planck,avogno,vellgt,boltz,bohr,todebye
  public epsil,small_,sqrt2,sqrt3,rad,fititermax,aston,hartree,ev,my_fmt
  !
  integer, parameter :: sik         = selected_int_kind(4)       ! Small integers
  integer, parameter :: ik          = selected_int_kind(8)       ! "Normal" integers. This must map on
                                                                 ! C "int" type, or the DX interface won't
                                                                 ! work correctly.
  integer, parameter :: hik         = selected_int_kind(8)       ! "Pointer" integers - sufficient to store
                                                                 ! memory address
  integer, parameter :: drk         = selected_real_kind(12,25)  ! "Double" reals and complex (complexi? :-)
  integer, parameter :: rk          = selected_real_kind(12,25)  ! "Normal" reals and complex (complexi? :-)
  integer, parameter :: ark         = selected_real_kind(25,32)  ! "Accurate" reals and complex (complexi? :-)
  integer, parameter :: inp         = 5                          ! Output I/O channel
  integer, parameter :: out         = 6                          ! Output I/O channel
  integer, parameter :: nfilelegendre = 101                      ! Damp-output channel for eigenfunction 


  real(rk), parameter  :: safe_max   = exp(log(huge(1.0_rk))*0.25_rk)   ! Largest number we want to work with
  real(rk), parameter  :: safe_min   = exp(log(tiny(1.0_rk))*0.25_rk)   ! Smallest number we want to work with
                                                                        ! (The somewhat convoluted syntax is used to standard F95 
                                                                        !  which prohibits non-integer exponents in this context)
  real(rk), parameter  :: max_exp    =log(safe_max)                     ! Largest number OK for exponentiating
  real(rk), parameter  :: small_     =epsilon(1.0_rk)                   ! a positive model number that is almost 
                                                                        ! negligible compared to unity in the current model
                                                                 ! epsil - antisymmetric tensor (Levi-Civita symbol)
  real(rk), parameter :: epsil(3,3,3)=reshape( (/0,0,0,0,0,-1,0,1,0,0,0,1,0,0,0,& 
                                                -1,0,0,0,-1,0,1,0,0,0,0,0/), (/3,3,3/) )
  real(rk), parameter :: pi    = 4.0_rk * atan2(1.0_rk,1.0_rk)   ! PI=3.14...
  real(rk), parameter :: twopi = 2.0_rk * pi                     ! 2*PI=6.28...
  real(rk), parameter :: sqrt2 = sqrt(2._rk)                     ! \sqrt{2}
  real(rk), parameter :: sqrt3 = sqrt(3._rk)                     ! \sqrt{3}
  real(rk), parameter :: rad   = 180._rk/pi                      ! radian = 180/pi
  integer, parameter  :: cl          = 80                        ! Max character string length
  integer, parameter  :: wl          = 500                       ! Very large max character string length 
  character(len=cl)   :: my_fmt                                  !text variable containing formats for reads/writes
  integer, parameter  :: fititermax  = 200                       ! Max number of iterations in different fittings 


  ! physical constants -- All constants updated 21 March 2012 from the NIST
  ! website http://physics.nist.gov/cuu/Constants/index.html
  real(drk), parameter :: planck     =  6.62606957e-27_rk         ! Planck constant in (non-SI) erg*second
  real(drk), parameter :: avogno     =  6.02214129e+23_rk         ! Avogadro constant
  real(drk), parameter :: vellgt     =  2.99792458e+10_rk         ! Speed of light constant in (non-SI) cm/second
  real(drk), parameter :: boltz      =  1.3806488e-16_rk          ! Boltzmann constant in (non-SI) erg/Kelvin
  real(drk), parameter :: bohr       =  0.52917720859_rk          ! bohr constant in Angstroms
  real(drk), parameter :: hartree    =  219474.6313705_rk         ! hartree in cm-1
  real(drk), parameter :: ev         =  8065.54465_rk             ! ev in cm-1
  real(drk), parameter :: uma        =  1.660538921e-24_rk        ! unified atomic mass unit [=mass(C-12)/12 ] in grams
  real(drk), parameter :: aston      =  planck/(8._rk*PI**2*vellgt*uma*1e-16_rk)  !rotational factor in cm-1 amu Ang^2
  real(drk), parameter :: todebye    =  2.541765_rk               ! a.u. in debye
  real(drk), parameter :: c2         =  1.4387751601679204d0                 ! second radiative constant
  
  contains

    subroutine accuracyInitialize  ! This initialization routine is now superfluous but kept for compatibility
     !
    end subroutine accuracyInitialize

end module accuracy

