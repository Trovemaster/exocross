module VoigtKampff_module
    use accuracy
    implicit none

    public VoigtKampff

    private

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


contains
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
end module VoigtKampff_module

