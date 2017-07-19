module VoigtKampffCollection_module
	use accuracy
	use VoigtKampff_module
	implicit none
	
	public 	VoigtKampffCollection
	
	private
	
	integer(ik),parameter		::	EXPAND_SIZE = 50
	
	type	::	 VoigtKampffCollection
		
		type(VoigtKampff),allocatable	::	fast_voigts(:)
		
		integer(ik)				::	capacity
		integer(ik)				::	n
		logical					::	normalize
		real(rk)				::	dpwcoeff
		real(rk)				::	offset
		real(rk)				::	dfreq
		
		
		
	contains
		procedure,public		::	construct
		procedure,public		::	expand_voigts
		procedure,public		::	generate_indices
		procedure,public		::	get_size
		procedure,private		::	add_voigt
		procedure,private		::	search_gamma
		procedure,public		::	compute => compute_fast_voigt
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
		class(VoigtKampffCollection)	::	this
		type(VoigtKampff),allocatable	::	temp(:)
		integer(ik)			::	old_capacity
		
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
		class(VoigtKampffCollection)	::	this
		real(rk),intent(in)		::	gammaL(0:Jmax,-1:1)
		integer(ik),intent(out)		::	gamma_idx(0:Jmax,-1:1)
		integer(ik),intent(in)		::	Jmax
		
		integer(ik)			::	ido,jdo
		real(rk)			::	val
		integer				::	idx
		
		!Search of the gamma value, if its similar then resue the index otherwise we add antoher voigt
		do ido=0,Jmax
			do jdo=-1,1
				val = gammaL(ido,jdo)
				idx = this%search_gamma(val)
				!If we didnt find it then add a voigt
				if(idx==0) then
					gamma_idx(ido,jdo) = this%add_voigt(val)
				else
					gamma_idx(ido,jdo) = idx
				endif
			enddo
		enddo
				
		
		
	
	end subroutine
	
	subroutine compute_fast_voigt(this,freq,intens,abscoef,ib,ie,start_nu,nu0,index)
	  	class(VoigtKampffCollection),intent(in)	::	this
	  	real(rk),intent(in)		::	freq(:),abscoef,start_nu,nu0
	  	integer(ik),intent(in)		::	ib,ie
	  	real(rk),intent(inout)		::	intens(:)
	  	integer,intent(in)		::	index
	  	if(index > this%n) stop "FastVoigt::Out of bounds"
	  	
	  	call this%fast_voigts(index)%compute(freq,intens,abscoef,ib,ie,start_nu,nu0)
	  	
	 end subroutine
	
	integer(ik) function add_voigt(this,gammaL)
		class(VoigtKampffCollection)	::	this
		real(rk),intent(in)			::	gammaL
		
		this%n = this%n+1
		
		if(this%n > this%capacity) call this%expand_voigts
		
		add_voigt = this%n
		
		call this%fast_voigts(this%n)%construct(this%dpwcoeff,gammaL,this%dfreq,this%offset,this%normalize)
		
	end function
		
		
		
	
	
	integer(ik) function search_gamma(this,gammaL)
		class(VoigtKampffCollection)	::	this
		real(rk),intent(in)		::	gammaL	
		
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








end module VoigtKampffCollection_module
