module Algo916_module
	use accuracy
	implicit none
	
	real(rk),parameter :: TOL=1.43e-17
contains


	subroutine sigma_ab(x,y,sigma1,sigma2,sigma3,a,ex2)
		real(rk),intent(in)		::	a,ex2,y
		real(rk),intent(out)		::	sigma1,sigma2,sigma3
		real(rk)			::	x
		real(rk)			::	f,f3p,f3n
		real(rk)			::	an,an3p,an3n
		
		real(rk)			::	yy,e2axn,e2axp,ena2n2
		
		integer(ik)			::	n,n0,n3p,n3n
		
		
		sigma1 = 0.0
		sigma2 = 0.0
		sigma3 = 0.0
		
		yy = y*y
		
		if(x < 0.0) x = -x
		
		n0 = int(ceiling(x/a))

		
		do n=1,8
			n3p=n0+n-1
			n3n=n0-n
			an = a*n
			
			ena2n2 = exp(-a*a*n*n)
			
			
			an3p = a*n3p
			an3n=a*n3n
			
			f = 1.0/ (an * an + yy)
			f3p = 1.0/ (an3p * an3p + yy)
			f3n = 1.0 / (an3n * an3n + yy)
			
			sigma1 = sigma1 + f * ena2n2*ex2
			sigma2 = sigma2 + f * exp(-2.0*a*x*n)*ena2n2*ex2
			sigma3 = sigma3 + f3p * exp(-(an3p - x) * (an3p - x))
	
			if(n3n >= 1) sigma3 = sigma3 + f3n * exp(-(an3n - x) * (an3n - x));
		enddo
	end subroutine


	real(rk) function algo916(x,y)
		real(rk),intent(in)	::	y
		real(rk)		::	x
		real(rk)		::	s1,s2,s3
		real(rk)		::	ex2,xy,a2ipi,cos2xy,sinxy,t1
		real(rk),parameter	::	a  = 0.5
		
		
		ex2 = exp(-x * x);
		
		if(x /= 0.0 .and. y /= 0.0) call sigma_ab(x, y, s1, s2, s3, a, ex2);
		
		xy = x * y
		a2ipi = 2.0 * a / pi
		cos2xy = cos(2.0 * xy)
		sinxy = sin(xy)

		t1 = ex2 *  ERFC_SCALED(y) * cos2xy
		t1 = t1 +  a2ipi * x * sinxy * ex2 * sinxy / xy
		t1 = t1 + a2ipi * y * (-cos2xy * s1 + 0.5 * (s2 + s3))
	
		if(x == 0) t1 = ERFC_SCALED(y)
		if(y == 0) t1 = ex2;	
		
		algo916 = t1
	end function	
	
	  function voigt_algo916(nu,nu0,alpha,gamma) result(f)
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
	   f = algo916(x, y) / sigma / RT2PI
	   !
	  end function voigt_algo916

end module
