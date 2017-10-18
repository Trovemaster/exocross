module VoigtKampffCollection_module
    use accuracy
    implicit none

    public 	VoigtKampffCollection

    private

    integer(ik),parameter :: EXPAND_SIZE = 50

    real(rk),parameter    ::    REFERENCE_NU0 = 1.0
    real(rk),parameter    ::    DISTANCE_MAGIC_NUMBER = 2.0

    type ::        VoigtKampff
        private
        integer(ik)    ::    m_Npoints
        real(rk)    ::    m_res
        real(rk)    ::    m_lorentz_cutoff
        integer(ik)    ::    m_middle_point
        real(rk)     ::    m_gammaD
        real(rk)     ::    m_gammaL
        real(rk)    ::    m_mag
        logical        ::    normalize
        logical        ::    constructed = .false.
        real(rk),pointer    ::    m_voigt_grid(:)


    contains
        procedure,public    ::    construct    =>    construct_voigt
    
        procedure,public    ::    compute        =>     compute_voigt
    
        procedure,public        ::    destroy

        procedure,public    ::    get_gammaL

    end type    
    
    
    
    type ::	 VoigtKampffCollection

        type(VoigtKampff),allocatable :: fast_voigts(:)
        integer(ik) :: capacity
        integer(ik) :: n
        logical     :: normalize
        real(rk)    :: dpwcoeff
        real(rk)    :: offset
        real(rk)    :: dfreq

    contains
        procedure,public  :: construct
        procedure,public  :: expand_voigts
        procedure,public  :: generate_indices
        procedure,public  :: get_size
        procedure,private :: add_voigt
        procedure,private :: search_gamma
        procedure,public  :: compute => compute_fast_voigt
    end type

    
    

contains

    subroutine construct(this,dpwcoeff,offset,dfreq,normalize)
        class(VoigtKampffCollection)	::	this
        logical,intent(in)		::	normalize
        real(rk),intent(in)		::	dpwcoeff,offset,dfreq

        this%n = 0
        this%capacity = 0

        this%dpwcoeff =dpwcoeff
        this%offset = offset
        this%dfreq = dfreq
        this%normalize = normalize

        call this%expand_voigts()

    end subroutine

    subroutine expand_voigts(this)
        class(VoigtKampffCollection) :: this
        type(VoigtKampff),allocatable :: temp(:)
        integer(ik) :: old_capacity

        if(this%capacity == 0) then
            this%capacity = this%capacity + EXPAND_SIZE
            allocate(this%fast_voigts(this%capacity))
            return
        endif

        allocate(temp(this%capacity))
        temp(:) = this%fast_voigts(:)

        deallocate(this%fast_voigts)

        old_capacity=this%capacity

        this%capacity = this%capacity + EXPAND_SIZE

        allocate(this%fast_voigts(this%capacity))

        this%fast_voigts(1:old_capacity) = temp(1:old_capacity)

        deallocate(temp)

    end subroutine

    integer(ik) function get_size(this)
        class(VoigtKampffCollection),intent(in)	::	this

        get_size = this%n
    end function


    subroutine generate_indices(this,gammaL,gamma_idx,Jmax)
        class(VoigtKampffCollection) :: this
        real(rk),intent(in) :: gammaL
        integer(ik),intent(inout) :: gamma_idx
        integer(ik),intent(in) :: Jmax

        integer(ik) :: ido,jdo
        real(rk)    :: val
        integer     :: idx

        !Search of the gamma value, if its similar then resue the index otherwise we add antoher voigt

        val = gammaL
        idx = this%search_gamma(val)
              !If we didnt find it then add a voigt
        if(idx==0) then
           gamma_idx = this%add_voigt(val)
        else
           gamma_idx = idx
        endif

    end subroutine

    subroutine compute_fast_voigt(this,freq,intens,abscoef,ib,ie,start_nu,nu0,index)
        class(VoigtKampffCollection),intent(in)	::	this
        real(rk),intent(in)     :: freq(:),abscoef,start_nu,nu0
        integer(ik),intent(in)  :: ib,ie
        real(rk),intent(inout)  :: intens(:)
        integer,intent(in)      :: index
        if(index > this%n) stop "FastVoigt::Out of bounds"

        call this%fast_voigts(index)%compute(freq,intens,abscoef,ib,ie,start_nu,nu0)

    end subroutine

    integer(ik) function add_voigt(this,gammaL)
        class(VoigtKampffCollection) :: this
        real(rk),intent(in) :: gammaL

        this%n = this%n+1

        if(this%n > this%capacity) call this%expand_voigts

        add_voigt = this%n

        call this%fast_voigts(this%n)%construct(this%dpwcoeff,gammaL,this%dfreq,this%offset,this%normalize)
   
    end function

    integer(ik) function search_gamma(this,gammaL)
        class(VoigtKampffCollection) :: this
        real(rk),intent(in) :: gammaL

        search_gamma = 0

        if(this%n == 0) return
        search_gamma = binarySearch_R(this%fast_voigts(1:this%n),gammaL)


    end function

recursive function binarySearch_R (a, value) result (bsresult)
    type(VoigtKampff), intent(in) :: a(:)
    real(rk)	::	value
    integer          :: bsresult, mid
 
    mid = size(a)/2 + 1
    if (size(a) == 0) then
        bsresult = 0        ! not found
    else if(is_close(a(mid)%get_gammaL(),value,1e-5_rk)) then
        bsresult = mid      ! SUCCESS!!
    else if (a(mid)%get_gammaL() > value) then
        bsresult= binarySearch_R(a(:mid-1), value)
    else if (a(mid)%get_gammaL() < value) then
        bsresult = binarySearch_R(a(mid+1:), value)
        if (bsresult /= 0) then
            bsresult = mid + bsresult
        end if
   endif
end function binarySearch_R



    !Check if the lorentzian is the same
    logical function is_close(a, b, rel_tol )
      real(rk),intent(in)	::	a,b,rel_tol
      	is_close = abs(a-b) <= rel_tol * max(abs(a), abs(b))
    end function



  subroutine construct_voigt(this,pGammaD,pGammaL,pRes,pLorentzCutoff,pNormalize)
      class(VoigtKampff),intent(inout)    ::    this
      real(rk),intent(in)            ::    pGammaD,pGammaL
      real(rk),intent(in)            ::    pRes
      real(rk),intent(in)            ::    pLorentzCutoff
      logical,intent(in)            ::    pNormalize

      integer(ik)                ::    ierr,ido

      real(rk)                ::    nu,nu0

  
      this%m_res = pRes
      this%m_gammaD = pGammaD
      this%m_gammaL = pGammaL
      this%m_lorentz_cutoff = pLorentzCutoff
     this%normalize = pNormalize
      this%m_Npoints = int(2.0*this%m_lorentz_cutoff/this%m_res,ik);
      this%m_middle_point = (this%m_Npoints/2) + 1
      !Construct the voigt grid

      this%m_mag = 0.0

      allocate(this%m_voigt_grid(this%m_Npoints),stat=ierr)

      !call ArrayStart('Voigt-F-Grid',info,size(this%m_voigt_grid),kind(this%m_voigt_grid))
      nu0 = REFERENCE_NU0
      nu = nu0-this%m_lorentz_cutoff

      do ido=1,this%m_Npoints
      
          this%m_voigt_grid(ido) = voigt_humlicek(nu,nu0,this%m_gammaD,this%m_gammaL)
          this%m_mag = this%m_mag + this%m_voigt_grid(ido)
          nu = nu+ this%m_res
      enddo

      this%m_mag = this%m_mag*this%m_res

      this%constructed = .true.

  end subroutine

  subroutine destroy(this)
     class(VoigtKampff),intent(inout)    ::    this
     
     if(.not. this%constructed)return
     
     if(associated(this%m_voigt_grid)) then
         deallocate(this%m_voigt_grid)
     endif
  end subroutine

  
  subroutine compute_voigt(this,freq,intens,abscoef,ib,ie,start_nu,nu0)
      class(VoigtKampff),intent(in)    ::    this
      real(rk),intent(in)        ::    freq(:),abscoef,start_nu,nu0
      integer(ik),intent(in)        ::    ib,ie
      real(rk),intent(inout)        ::    intens(:)

  
      integer                ::    center_point,dist
      integer                ::    middle_shift
      integer                ::    start_dist,end_dist,left_start,left_end,right_start,right_end
      integer                ::    ib_rel,ie_rel,num_hum_points,ido
      integer                ::    Npoints
      real(rk)            ::    gammaD,hum_res,mag,nu
      real(rk),allocatable        ::    temp_humlicek(:)

  
      center_point = (nu0-start_nu)/this%m_res + 1
      middle_shift = (center_point - this%m_middle_point) + 1

  

      dist = max(DISTANCE_MAGIC_NUMBER/this%m_res,3.0)

      Npoints = ie - ib + 1

      ib_rel = ib - middle_shift
    ie_rel = ie - middle_shift

      gammaD = this%m_gammaD*nu0


    mag = this%m_mag

    start_dist = max(center_point - dist + 1, ib)
    end_dist = min(center_point + dist, ie)
    
    left_start = ib_rel
    left_end = min(this%m_middle_point - dist, ie_rel)
    right_start = max(this%m_middle_point + dist, ib_rel)
    right_end = ie_rel

    !print *,mag


    !If we have any humlicek points
    if(start_dist < end_dist) then


        if(this%normalize) then 
            num_hum_points = end_dist - start_dist + 1;
            allocate(temp_humlicek(num_hum_points))

            do ido=start_dist,end_dist
                nu=freq(ido)
                hum_res = voigt_humlicek(nu,nu0,gammaD,this%m_gammaL)
                temp_humlicek(ido-start_dist+1) = hum_res
                mag = mag - this%m_voigt_grid(ido -middle_shift+1)*this%m_res
                mag = mag + hum_res*this%m_res
            enddo
    
    
            intens(start_dist:end_dist) = intens(start_dist:end_dist) + temp_humlicek(:)*abscoef/mag
    
            !print *,mag
            deallocate(temp_humlicek)
        else
            do ido=start_dist,end_dist
                nu=freq(ido)
                intens(ido) = intens(ido) + voigt_humlicek(nu,nu0,gammaD,this%m_gammaL)*abscoef
            enddo    
        
            mag = 1.0
    
    
        endif
      endif
      if (left_start < left_end ) then
          call vectorized_voigt(intens(ib-ib_rel+left_start:ib-ib_rel+left_end),this%m_voigt_grid(left_start:left_end),(abscoef/mag))

      endif

      if (right_start < right_end ) then
          call vectorized_voigt(intens(ib-ib_rel+right_start:ib-ib_rel+right_end),this%m_voigt_grid(right_start:right_end),(abscoef/mag))
      endif

  end subroutine
          
  real(rk) function get_gammaL(this)
      class(VoigtKampff),intent(in)    ::    this

      get_gammaL = this%m_gammaL

      return
  end function




!-----------Humlicek methods----------------------!

 ! This should be heavily vectorized by the compiler
  subroutine vectorized_voigt(intens,voigt,abscoeff)
      real(rk),intent(inout)    ::    intens(:)
      real(rk),intent(in)    ::    voigt(:)
      real(rk),intent(in)    ::    abscoeff

      intens(:) = intens(:)+voigt(:)*abscoeff
  end subroutine


  function voigt_humlicek(nu,nu0,alpha,gamma) result(f)
   !
   real(rk),intent(in) :: nu,nu0,alpha,gamma
   real(rk) :: f,sigma,x,y

   real(rk),parameter ::  RT2LN2 = 1.1774100225154747_rk     ! sqrt(2ln2)
   real(rk),parameter ::  RTPI   = 1.7724538509055159_rk     ! sqrt(pi)
   real(rk),parameter ::  RT2PI  = 2.5066282746310002_rk     ! sqrt(2pi)
   real(rk),parameter ::  RT2    = 1.4142135623730951_rk     ! sqrt(pi)
   !
   sigma = alpha / RT2LN2
   x = (nu0 - nu) / sigma / RT2
   y = gamma / sigma / RT2
   !
   f = humlic(x, y) / sigma / RT2PI
   !
  end function voigt_humlicek

  !
  function humlic(x,y) result (f)
    real(rk),intent(in) :: x,y
    real(rk) :: f
    complex(rk) :: t,humlic1,u
    real(rk) :: s
    !
    T = cmplx(y,-x)
    S = abs(x) + y
    !
    if (S >= 15) then
        ! Region I
        humlic1 = T*0.5641896_rk/(0.5_rk+T*T)
        !fprintf(stdout, "I")
        !
        f = real(humlic1,rk)
        return
    endif
    !
    if (S >= 5.5) then
        ! Region II
        U = T * T;
        humlic1 = T * (1.410474_rk + U*0.5641896_rk)/(0.75_rk + U*(3.0_rk+U))
        !
        !fprintf(stdout, "II");
        !
        f = real(humlic1,rk)
        return
        !
    endif
    !
    if (y >= 0.195_rk * abs(x) - 0.176_rk) then
        ! Region III
        humlic1 = (16.4955_rk+T*(20.20933_rk+T*(11.96482_rk &
                +T*(3.778987_rk+T*.5642236_rk)))) / (16.4955_rk+T*(38.82363_rk &
                +T*(39.27121_rk+T*(21.69274_rk+T*(6.699398_rk+T)))))
        !fprintf(stdout, "III");
        f = real(humlic1,rk)
        return
        !
    endif
    !
    ! Region IV
    U = T * T
    humlic1 = exp(U)-T*(36183.31_rk-U*(3321.9905_rk-U*(1540.787_rk-U*(219.0313_rk-U* &
       (35.76683_rk-U*(1.320522_rk-U*.56419_rk))))))/(32066.6_rk-U*(24322.84_rk-U* &
       (9022.228_rk-U*(2186.181_rk-U*(364.2191_rk-U*(61.57037_rk-U*(1.841439_rk-U)))))))
    !fprintf(stdout, "IV");
    !
    f = real(humlic1,rk)
    !
    end function humlic




end module VoigtKampffCollection_module
