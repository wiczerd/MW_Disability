module nr_lib

  use nrtype
  use nrutil, only : arth_d, assert, assert_eq, nrerror
	
  implicit none

	INTERFACE gammln
		MODULE PROCEDURE gammln_s, gammln_v
	END INTERFACE

	INTERFACE gser
		MODULE PROCEDURE gser_s, gser_v
	END INTERFACE

	INTERFACE gcf
		MODULE PROCEDURE gcf_s, gcf_v
	END INTERFACE

	INTERFACE gammp
		MODULE PROCEDURE gammp_s, gammp_v
	END INTERFACE

	INTERFACE ran1
		MODULE PROCEDURE ran1_s, ran1_v
	END INTERFACE

	INTERFACE locate
		MODULE PROCEDURE locate_s, locate_d
	END INTERFACE

contains
	
	!-------------------------------------------------------------------------------------------------------------------------------------------
	function linspace(a,b,n)
		
		real(DP), intent(in) :: a,b
		integer, intent(in) :: n
		real(DP), dimension(n) :: linspace	
		real(DP) :: d
		integer :: i
			
		d=(b-a)/(real(n-1))
		linspace(1)=a
		do i=2,n
			linspace(i)=linspace(i-1)+d
		enddo
			
	end function linspace
	!--------------------------------------------------------------------------------------------------------------------------------------------

	!--------------------------------------------------------------------------------------------------------------------------------------------
	subroutine tauchen(sigma_eps,rho,m,N,states,trans_mat)

		use nrtype
		implicit none

		real(DP), intent(in) :: sigma_eps,rho
		integer, intent(in) :: m, N
		real(DP), dimension(N), intent(inout) :: states
		real(DP), dimension(N,N), intent(inout) :: trans_mat

		integer :: j,k
		real(DP) :: d

		states=linspace(-m*(sigma_eps/(1.0_dp-rho**2.0_dp))**(0.5_dp),m*(sigma_eps/(1.0_dp-rho**2.0_dp))**(0.5_dp),N)
		d=states(2)-states(1)

		!write(*,*) cdf(0.0_dp)
		!write(*,*) cdf(1.0_dp)
		!write(*,*) cdf(2.0_dp)

		do j=1,N
			trans_mat(j,1)=cdf((states(1)+d/2.0_dp-rho*states(j))/sqrt(sigma_eps))
			do k=2,N-1
				trans_mat(j,k)=cdf((states(k)+d/2.0_dp-rho*states(j))/sqrt(sigma_eps))-cdf((states(k)-d/2.0_dp-rho*states(j))/sqrt(sigma_eps))
			enddo
			trans_mat(j,N)=1.0_dp-cdf((states(N)-d/2.0_dp-rho*states(j))/sqrt(sigma_eps))
		enddo

		contains
			function cdf(eps)
				real(DP), intent(in) :: eps
				real(DP) :: cdf
				cdf=0.5_dp+0.5_dp*erf_s(eps/sqrt(2.0_dp))
			end function cdf

	end subroutine tauchen
	!--------------------------------------------------------------------------------------------------------------------------------------------
	subroutine rouwenhorst(sigma_eps, rho, N, states, trans_mat)
		use nrtype
		implicit none

		real(DP), intent(in) :: sigma_eps, rho
		integer, intent(in) :: N
		real(DP), dimension(N), intent(inout) :: states
		real(DP), dimension(N,N), intent(inout) :: trans_mat

		integer :: i, j
		real(DP) :: sigma_z, psi, p, q

		sigma_z=sqrt(sigma_eps**2.0_dp/(1-rho**2.0_dp))
		psi=sqrt(N-1.0_dp)*sigma_z
		do i=1,N
			states(i)=-psi+2.0_dp*psi/(N-1.0_dp)*(i-1.0_dp)
		enddo

		p=(1.0_dp+rho)/2.0_dp
		q=p

		trans_mat=pi_rouw(N,p,q)

		contains
 
			recursive function pi_rouw(N,p,q) result(res)

				implicit none

				integer, intent(in) :: N
				real(DP), intent(in) :: p,q
				real(DP), dimension(N,N) :: res

				real(DP), dimension(N-1,N-1) :: tmp
				integer :: i, j

				if (N==2) then
					res(1,1)=p
					res(1,2)=1-p
					res(2,1)=1-q
					res(2,2)=q
					return
				else
					tmp=pi_rouw(N-1,p,q)
					do j=1,N
						res(1,j)=factorial(N-1)/(factorial(j-1)*factorial((N-1)-(j-1)))*p**(N-j)*(1-p)**(j-1)
						res(N,j)=factorial(N-1)/(factorial(j-1)*factorial((N-1)-(j-1)))*(1-q)**(N-j)*q**(j-1)
					enddo
					do i=2,N-1
						res(i,1)=p*tmp(i,1)
						do j=2,N-1
							res(i,j)=p*tmp(i,j)+(1-p)*tmp(i,j-1)
						enddo
						res(i,N)=(1-p)*tmp(i,N-1)
					enddo
				endif

			end function pi_rouw

	end subroutine rouwenhorst
	!--------------------------------------------------------------------------------------------------------------------------------------------
	function factorial(k)
		implicit none
		integer, intent(in) :: k
		integer :: factorial
		integer :: i
		factorial = 1		
		do i=1,k
			factorial=factorial*i
		enddo		
	end function factorial
	!--------------------------------------------------------------------------------------------------------------------------------------------

	!--------------------------------------------------------------------------------------------------------------------------------------------
	FUNCTION locate_s(xx,x)
	
		USE nrtype
		
		IMPLICIT NONE
		
		REAL(SP), DIMENSION(:), INTENT(IN) :: xx
		REAL(SP), INTENT(IN) :: x
		INTEGER(I4B) :: locate_s
		
		!Given an array xx(1:N), and given a value x, returns a value j such that x is between
		!xx(j) and xx(j + 1). xx must be monotonic, either increasing or decreasing. j = 0 or
		!j = N is returned to indicate that x is out of range.
		
		INTEGER(I4B) :: n,jl,jm,ju
		LOGICAL :: ascnd
		
		n=size(xx)
		ascnd = (xx(n) >= xx(1))
		jl=0 
		ju=n+1 
		do
			if (ju-jl <= 1) exit
			jm=(ju+jl)/2
			if (ascnd .eqv. (x >= xx(jm))) then
				jl=jm 
			else
				ju=jm 
			end if
		end do
		if (x == xx(1)) then
			locate_s=1
		else if (x == xx(n)) then
			locate_s=n-1
		else
			locate_s=jl
		end if
	END FUNCTION locate_s
	!--------------------------------------------------------------------------------------------------------------------------------------------

	!--------------------------------------------------------------------------------------------------------------------------------------------
	FUNCTION locate_d(xx,x)
	
		USE nrtype
		
		IMPLICIT NONE
		
		REAL(DP), DIMENSION(:), INTENT(IN) :: xx
		REAL(DP), INTENT(IN) :: x
		INTEGER(I4B) :: locate_d
		
		!Given an array xx(1:N), and given a value x, returns a value j such that x is between
		!xx(j) and xx(j + 1). xx must be monotonic, either increasing or decreasing. j = 0 or
		!j = N is returned to indicate that x is out of range.
		
		INTEGER(I4B) :: n,jl,jm,ju
		LOGICAL :: ascnd
		
		n=size(xx)
		ascnd = (xx(n) >= xx(1))
		jl=0 
		ju=n+1 
		do
			if (ju-jl <= 1) exit
			jm=(ju+jl)/2
			if (ascnd .eqv. (x >= xx(jm))) then
				jl=jm 
			else
				ju=jm 
			end if
		end do
		if (x == xx(1)) then
			locate_d=1
		else if (x == xx(n)) then
			locate_d=n-1
		else
			locate_d=jl
		end if
	END FUNCTION locate_d
	!--------------------------------------------------------------------------------------------------------------------------------------------
	
	!--------------------------------------------------------------------------------------------------------------------------------------------
	SUBROUTINE tridag_ser(a,b,c,r,u)
	
		USE nrtype;
		USE nrutil, ONLY : assert_eq,nrerror
		
		IMPLICIT NONE
		REAL(DP), DIMENSION(:), INTENT(IN) :: a,b,c,r
		REAL(DP), DIMENSION(:), INTENT(OUT) :: u
	
		!Solves for a vector u of size N the tridiagonal linear set given by equation (2.4.1) using a
		!serial algorithm. Input vectors b (diagonal elements) and r (right-hand sides) have size N,
		!while a and c (off-diagonal elements) are size N - 1.

		REAL(DP), DIMENSION(size(b)) :: gam
		INTEGER(I4B) :: n,j
		REAL(DP) :: bet
		n=assert_eq((/size(a)+1,size(b),size(c)+1,size(r),size(u)/),"tridag_ser")
		bet=b(1)
		if (bet == 0.0) call nrerror("tridag_ser: Error at code stage 1")
		u(1)=r(1)/bet
		do j=2,n
			gam(j)=c(j-1)/bet
			bet=b(j)-a(j-1)*gam(j)
			if (bet == 0.0) &
				call nrerror("tridag_ser: Error at code stage 2")
			u(j)=(r(j)-a(j-1)*u(j-1))/bet
		end do
		do j=n-1,1,-1
			u(j)=u(j)-gam(j+1)*u(j+1)
		end do
	END SUBROUTINE tridag_ser
	!--------------------------------------------------------------------------------------------------------------------------------------------
	
	!--------------------------------------------------------------------------------------------------------------------------------------------
	! LINEAR INTERPOLATION
	!--------------------------------------------------------------------------------------------------------------------------------------------
	function lin_interp(xa,ya,x)
	
		use nrtype
		USE nrutil, ONLY : assert_eq
		
		real(DP), dimension(:), intent(in) :: xa, ya
		real(DP), intent(in) :: x
		real(DP) :: lin_interp
		
		!given arrays xa and ya of length n, the first containing the grid and the second containing the corresponding function values, calculates
		!linear_interpolation of the function at x
		
		integer(I4B) :: n, j
		real(DP) :: A, B
		
		n=assert_eq(size(xa),size(ya),"linest")
		j=locate(xa,x)
		A=(xa(j+1)-x)/(xa(j+1)-xa(j))
		B=1-A
		
		lin_interp=A*ya(j)+B*ya(j+1)		
		
	end function lin_interp
	!--------------------------------------------------------------------------------------------------------------------------------------------
	function lin_interp_d(xa,ya,x)
	
		use nrtype
		USE nrutil, ONLY : assert_eq
		
		real(DP), dimension(:), intent(in) :: xa, ya
		real(DP), intent(in) :: x
		real(DP) :: lin_interp_d
		
		!given arrays xa and ya of length n, the first containing the grid and the second containing the corresponding function values, calculates
		!linear_interpolation of the function at x
		
		integer(I4B) :: n, j
		real(DP) :: A_prime, B_prime
		
		n=assert_eq(size(xa),size(ya),"linest")
		j=locate(xa,x)

		!I ADD IN EXTRAPOLATION HERE
		IF (j <1) THEN
			lin_interp_d = ya(2)+(xa(2)-x)*(ya(2)-ya(1))/(xa(2)-xa(1))

		ELSEIF (j == n) THEN
			lin_interp_d = ya(n-1)+(x-xa(n-1))*(ya(n)-ya(n-1))/(xa(n)-xa(n-1))

		ELSE
			A_prime=-1/(xa(j+1)-xa(j))
			B_prime=-A_prime
		
			lin_interp_d=A_prime*ya(j)+B_prime*ya(j+1)		

		ENDIF
		
	end function lin_interp_d
	!------------------------------------------------------------------------------------------------------------------------------------------------
	
	!--------------------------------------------------------------------------------------------------------------------------------------------
	! CUBIC SPLINE INTERPOLATION ROUTINES
	!--------------------------------------------------------------------------------------------------------------------------------------------
	SUBROUTINE spline(x,y,yp1,ypn,y2)
	
		USE nrtype;
		USE nrutil, ONLY : assert_eq
		!USE nr, ONLY : tridag
		
		IMPLICIT NONE
		
		REAL(DP), DIMENSION(:), INTENT(IN) :: x,y
		REAL(DP), INTENT(IN) :: yp1,ypn
		REAL(DP), DIMENSION(:), INTENT(OUT) :: y2

		!Given arrays x and y of length N containing a tabulated function, i.e., yi = f(xi), with x1 <
		!x2 < ... < xN, and given values yp1 and ypn for the first derivative of the interpolating
		!function at points 1 and N, respectively, this routine returns an array y2 of length N
		!that contains the second derivatives of the interpolating function at the tabulated points
		!xi. If yp1 and/or ypn are equal to 1 × 1030 or larger, the routine is signaled to set the
		!corresponding boundary condition for a natural spline, with zero second derivative on that
		!boundary.
		
		INTEGER(I4B) :: n
		REAL(DP), DIMENSION(size(x)) :: a,b,c,r
		n=assert_eq(size(x),size(y),size(y2),"spline")
		c(1:n-1)=x(2:n)-x(1:n-1)
		r(1:n-1)=6.0_sp*((y(2:n)-y(1:n-1))/c(1:n-1))
		r(2:n-1)=r(2:n-1)-r(1:n-2)
		a(2:n-1)=c(1:n-2)
		b(2:n-1)=2.0_sp*(c(2:n-1)+a(2:n-1))
		b(1)=1.0
		b(n)=1.0
		if (yp1 > 0.99e30_dp) then 
			r(1)=0.0
			c(1)=0.0
		else
			r(1)=(3.0_dp/(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
			c(1)=0.5
		end if 
		if (ypn > 0.99e30_dp) then
			r(n)=0.0
			a(n)=0.0
		else
			r(n)=(-3.0_dp/(x(n)-x(n-1)))*((y(n)-y(n-1))/(x(n)-x(n-1))-ypn)
			a(n)=0.5
		end if
		
		call tridag_ser(a(2:n),b(1:n),c(1:n-1),r(1:n),y2(1:n))
	END SUBROUTINE spline	
	!------------------------------------------------------------------------------------------------
	FUNCTION splint(xa,ya,y2a,x)
	
		USE nrtype;	
		USE nrutil, ONLY : assert_eq,nrerror
		!USE nr, ONLY: locate
		
		IMPLICIT NONE

		REAL(DP), DIMENSION(:), INTENT(IN) :: xa,ya,y2a
		REAL(DP), INTENT(IN) :: x
		REAL(DP) :: splint

		!Given the arrays xa and ya, which tabulate a function (with the xai ’s in increasing or
		!decreasing order), and given the array y2a, which is the output from spline above, and
		!given a value of x, this routine returns a cubic-spline interpolated value. The arrays xa, ya
		!and y2a are all of the same size.

		INTEGER(I4B) :: khi,klo,n
		REAL(SP) :: a,b,h
		n=assert_eq(size(xa),size(ya),size(y2a),"splint")
		klo=max(min(locate(xa,x),n-1),1)
		
		!We will find the right place in the table by means of locate’s bisection algorithm. This is
		!optimal if sequential calls to this routine are at random values of x. If sequential calls are in
		!order, and closely spaced, one would do better to store previous values of klo and khi and
		!test if they remain appropriate on the next call.

		khi=klo+1
		h=xa(khi)-xa(klo)
		if (h == 0.0) call nrerror("bad xa input in splint")
		a=(xa(khi)-x)/h
		b=(x-xa(klo))/h
		splint=a*ya(klo)+b*ya(khi)+((a**3-a)*y2a(klo)+(b**3-b)*y2a(khi))*(h**2)/6.0_dp
		
	END FUNCTION splint
	!-------------------------------------------------------------------------------------------------------------------------------------------
	function dsplint(xa,ya,y2a,x)
		USE nrtype;	
		USE nrutil, ONLY : assert_eq,nrerror
		!USE nr, ONLY: locate
		
		IMPLICIT NONE

		REAL(DP), DIMENSION(:), INTENT(IN) :: xa,ya,y2a
		REAL(DP), INTENT(IN) :: x
		REAL(DP) :: dsplint
		
		!Given the arrays xa and ya, which tabulate a function (with the xai ’s in increasing or
		!decreasing order), and given the array y2a, which is the output from spline above, and
		!given a value of x, this routine returns a cubic-spline interpolated derivative value. The arrays xa, ya
		!and y2a are all of the same size.
		
		INTEGER(I4B) :: khi,klo,n
		REAL(SP) :: a,b,h
		n=assert_eq(size(xa),size(ya),size(y2a),"splint")
		klo=max(min(locate(xa,x),n-1),1)
		
		!We will find the right place in the table by means of locate’s bisection algorithm. This is
		!optimal if sequential calls to this routine are at random values of x. If sequential calls are in
		!order, and closely spaced, one would do better to store previous values of klo and khi and
		!test if they remain appropriate on the next call.

		khi=klo+1
		h=xa(khi)-xa(klo)
		if (h == 0.0) call nrerror("bad xa input in splint")
		a=(xa(khi)-x)/h
		b=(x-xa(klo))/h
		dsplint=(ya(khi)-ya(klo))/h-(3.0_dp*a**2.0_dp-1.0_dp)*h*y2a(klo)/6.0_dp+(3.0_dp*b**2.0_dp-1.0_dp)*h*y2a(khi)/6.0_dp
		
	end function dsplint	
	!-------------------------------------------------------------------------------------------------------------------------------------------
	! BRENT'S METHOD FOR MINIMIZATION IN ONE DIMENSION WITH DERIVATIVE INFORMATION
	!-------------------------------------------------------------------------------------------------------------------------------------------
	FUNCTION dbrent(ax,bx,cx,func,dfunc,tol,xmin)
	
		USE nrtype
		USE nrutil, ONLY : nrerror
		
		IMPLICIT NONE
		REAL(DP), INTENT(IN) :: ax,bx,cx,tol
		REAL(DP), INTENT(OUT) :: xmin
		REAL(DP) :: dbrent
		INTERFACE
			FUNCTION func(x)
				USE nrtype
				IMPLICIT NONE
				REAL(DP), INTENT(IN) :: x
				real(DP) :: func
			END FUNCTION func
			FUNCTION dfunc(x)
				USE nrtype
				IMPLICIT NONE
				REAL(DP), INTENT(IN) :: x
				REAL(DP) :: dfunc
			END FUNCTION dfunc
		END INTERFACE
		INTEGER(I4B), PARAMETER :: ITMAX=100
		REAL(DP), PARAMETER :: ZEPS=1.0e-3_dp*epsilon(ax)

		!Given a function func and its derivative function dfunc, and given a bracketing triplet of
		!abscissas ax, bx, cx [such that bx is between ax and cx, and func(bx) is less than both
		!func(ax) and func(cx)], this routine isolates the minimum to a fractional precision of
		!about tol using a modification of Brent’s method that uses derivatives. The abscissa of
		!the minimum is returned as xmin, and the minimum function value is returned as dbrent,
		!the returned function value.
		!Parameters: Maximum allowed number of iterations, and a small number that protects
		!against trying to achieve fractional accuracy for a minimum that happens to be exactly
		!zero.

		INTEGER(I4B) :: iter
		REAL(DP) :: a,b,d,d1,d2,du,dv,dw,dx,e,fu,fv,fw,fx,olde,tol1,tol2,&
		u,u1,u2,v,w,x,xm
		
		!Comments following will point out only differences from the routine brent. Read that
		!routine first.
		
		LOGICAL :: ok1,ok2 !Will be used as flags for whether proposed steps are acceptable or not.
		a=min(ax,cx) 
		b=max(ax,cx)
		v=bx
		w=v
		x=v
		e=0.0_dp
		fx=func(x)
		fv=fx
		fw=fx
		dx=dfunc(x) !All our housekeeping chores are doubled by the necessity of moving
								!derivative values around as well	as function values.
		dv=dx
		dw=dx
		do iter=1,ITMAX
			xm=0.5_dp*(a+b)
			tol1=tol*abs(x)+ZEPS
			tol2=2.0_dp*tol1
			if (abs(x-xm) <= (tol2-0.5_dp*(b-a))) exit
			if (abs(e) > tol1) then
				d1=2.0_dp*(b-a) !Initialize these d’s to an out-of-bracket value.
				d2=d1 
				if (dw /= dx) d1=(w-x)*dx/(dx-dw) !Secant method with each point.
				if (dv /= dx) d2=(v-x)*dx/(dx-dv)
				
				!Which of these two estimates of d shall we take? We will insist that they be within
				!the bracket, and on the side pointed to by the derivative at x:
				
				u1=x+d1
				u2=x+d2
				ok1=((a-u1)*(u1-b) > 0.0) .and. (dx*d1 <= 0.0)
				ok2=((a-u2)*(u2-b) > 0.0) .and. (dx*d2 <= 0.0)
				olde=e !Movement on the step before last.
				e=d
				if (ok1 .or. ok2) then !Take only an acceptable d, and if both are acceptable, then take
															 !the smallest one.
					if (ok1 .and. ok2) then
						d=merge(d1,d2, abs(d1) < abs(d2))
					else
						d=merge(d1,d2,ok1)
					end if
					if (abs(d) <= abs(0.5_dp*olde)) then
						u=x+d
						if (u-a < tol2 .or. b-u < tol2) &
							d=sign(tol1,xm-x)
					else
						e=merge(a,b, dx >= 0.0)-x
						!Decide which segment by the sign of the derivative.
						d=0.5_dp*e !Bisect, not golden section.
					end if
				else
					e=merge(a,b, dx >= 0.0)-x
					d=0.5_dp*e !Bisect, not golden section.
				end if
			else
				e=merge(a,b, dx >= 0.0)-x
				d=0.5_dp*e !Bisect, not golden section.
			end if
			if (abs(d) >= tol1) then
				u=x+d
				fu=func(u)
			else
				u=x+sign(tol1,d)
				fu=func(u)
				!If the minimum step in the downhill direction takes us uphill, then we are done.
				if (fu > fx) exit
			end if
			du=dfunc(u) !Now all the housekeeping, sigh.
			if (fu <= fx) then
				if (u >= x) then
					a=x
				else
					b=x
				end if
				call mov3(v,fv,dv,w,fw,dw)
				call mov3(w,fw,dw,x,fx,dx)
				call mov3(x,fx,dx,u,fu,du)
			else
				if (u < x) then
					a=u
				else
					b=u
				end if
				if (fu <= fw .or. w == x) then
					call mov3(v,fv,dv,w,fw,dw)
					call mov3(w,fw,dw,u,fu,du)
				else if (fu <= fv .or. v == x .or. v == w) then
					call mov3(v,fv,dv,u,fu,du)
				end if
			end if
		end do
		
		if (iter > ITMAX) call nrerror("dbrent:exceeded maximum iterations")
		xmin=x
		dbrent=fx

		CONTAINS
			SUBROUTINE mov3(a,b,c,d,e,f)
				REAL(DP), INTENT(IN) :: d,e,f
				REAL(DP), INTENT(OUT) :: a,b,c
				a=d
				b=e
				c=f
			END SUBROUTINE mov3
	
	END FUNCTION dbrent
	!-----------------------------------------------------------------------------------------------

	!-----------------------------------------------------------------------------------------------
	! ZBRENT FOR ROOT FINDING
	!-----------------------------------------------------------------------------------------------
	FUNCTION zbrent(func,x1,x2,tol)
		USE nrtype; USE nrutil, ONLY :nrerror
		IMPLICIT NONE
		REAL(DP), INTENT(IN) :: x1,x2,tol
		REAL(DP) :: zbrent
		INTERFACE
			FUNCTION func(x)
				USE nrtype
				IMPLICIT NONE
				REAL(DP), INTENT(IN) :: x
				REAL(DP) :: func
			END FUNCTION func
		END INTERFACE
		INTEGER(I4B), PARAMETER :: ITMAX=100
		REAL(DP), PARAMETER :: EPS=epsilon(x1)

		!Using Brent’s method, find the root of a function func known to lie between x1 and x2.
		!The root, returned as zbrent, will be refined until its accuracy is tol.
		!Parameters: Maximum allowed number of iterations, and machine floating-point precision.

		INTEGER(I4B) :: iter
		REAL(DP) :: a,b,c,d,e,fa,fb,fc,p,q,r,s,tol1,xm
		a=x1
		b=x2
		fa=func(a)
		fb=func(b)
		if ((fa > 0.0 .and. fb > 0.0) .or. (fa < 0.0 .and. fb < 0.0)) &
		call nrerror("root must be bracketed for zbrent")
		c=b
		fc=fb
		do iter=1,ITMAX
		if ((fb > 0.0 .and. fc > 0.0) .or. (fb < 0.0 .and. fc < 0.0)) then
			c=a !Rename a, b, c and adjust bounding interval d.
			fc=fa 
			d=b-a
			e=d
		end if
		if (abs(fc) < abs(fb)) then
			a=b
			b=c
			c=a
			fa=fb
			fb=fc
			fc=fa
		end if
		tol1=2.0_dp*EPS*abs(b)+0.5_dp*tol !Convergence check.
		xm=0.5_dp*(c-b)
		if (abs(xm) <= tol1 .or. fb == 0.0) then
			zbrent=b
			RETURN
		end if
		if (abs(e) >= tol1 .and. abs(fa) > abs(fb)) then
			s=fb/fa !Attempt inverse quadratic interpolation.
			if (a == c) then
				p=2.0_sp*xm*s
				q=1.0_sp-s
				else
				q=fa/fc
				r=fb/fc
				p=s*(2.0_dp*xm*q*(q-r)-(b-a)*(r-1.0_dp))
				q=(q-1.0_dp)*(r-1.0_dp)*(s-1.0_dp)
			end if
			if (p > 0.0) q=-q !Check whether in bounds.
			p=abs(p)
			if (2.0_dp*p < min(3.0_dp*xm*q-abs(tol1*q),abs(e*q))) then
				e=d !Accept interpolation.
				d=p/q
				else
				d=xm !Interpolation failed; use bisection.
				e=d
			end if
			else !Bounds decreasing too slowly; use bisection
			d=xm 
			e=d
		end if
		a=b !Move last best guess to a.
		fa=fb
		b=b+merge(d,sign(tol1,xm), abs(d) > tol1 ) !Evaluate new trial root.
		fb=func(b)
		end do
		call nrerror("zbrent:exceeded maximum iterations")
		zbrent=b
	END FUNCTION zbrent
	!-----------------------------------------------------------------------------------------------

	!-----------------------------------------------------------------------------------------------
	! BRENT'S METHOD (NO DERIVATIVES)
	!-----------------------------------------------------------------------------------------------
	FUNCTION brent(ax,bx,cx,func,tol,xmin)
		USE nrtype; USE nrutil, ONLY :nrerror
		IMPLICIT NONE
		REAL(DP), INTENT(IN) :: ax,bx,cx,tol
		REAL(DP), INTENT(OUT) :: xmin
		REAL(DP) :: brent
		INTERFACE
			FUNCTION func(x)
				USE nrtype
				IMPLICIT NONE
				REAL(DP), INTENT(IN) :: x
				REAL(DP) :: func
			END FUNCTION func
		END INTERFACE
		INTEGER(I4B), PARAMETER :: ITMAX=100
		REAL(DP), PARAMETER :: CGOLD=0.3819660_dp,ZEPS=1.0e-3_dp*epsilon(ax)

		!Given a function func, and given a bracketing triplet of abscissas ax, bx, cx (such that bx
		!is between ax and cx,a n d func(bx) is less than both func(ax) and func(cx)), this
		!routine isolates the minimum to a fractional precision of about tol using Brent’s method.
		!The abscissa of the minimum is returned as xmin, and the minimum function value is
		!returned as brent, the returned function value.
		!Parameters: Maximum allowed number of iterations;golden ratio;and a small number that
		!protects against trying to achieve fractional accuracy for a minimum that happens to be
		!exactly zero.
		
		INTEGER(I4B) :: iter
		REAL(DP) :: a,b,d,e,etemp,fu,fv,fw,fx,p,q,r,tol1,tol2,u,v,w,x,xm
		a=min(ax,cx) !a and b must be in ascending order, though the input abscissas need not be.
		b=max(ax,cx) 
		v=bx !Initializations...
		w=v
		x=v
		e=0.0 !This will be the distance moved on the step before last.
		fx=func(x) 
		fv=fx
		fw=fx
		do iter=1,ITMAX !Main program loop.
			xm=0.5_dp*(a+b)
			tol1=tol*abs(x)+ZEPS
			tol2=2.0_dp*tol1
			if (abs(x-xm) <= (tol2-0.5_dp*(b-a))) then !Test for done here.
				xmin=x !Arrive here ready to exit with best values.
				brent=fx
				RETURN
			end if
			if (abs(e) > tol1) then !Construct a trial parabolic ?t.
				r=(x-w)*(fx-fv)
				q=(x-v)*(fx-fw)
				p=(x-v)*q-(x-w)*r
				q=2.0_sp*(q-r)
				if (q > 0.0) p=-p
				q=abs(q)
				etemp=e
				e=d
				if (abs(p) >= abs(0.5_dp*q*etemp) .or. p <= q*(a-x) .or. p >= q*(b-x)) then
					!The above conditions determine the acceptability of the parabolic ?t. Here it is
					!not o.k., so we take the golden section step into the larger of the two segments.
					e=merge(a-x,b-x, x >= xm )
					d=CGOLD*e
				else !Take the parabolic step.
					d=p/q
					u=x+d
					if (u-a < tol2 .or. b-u < tol2) d=sign(tol1,xm-x)
				end if
			else !Take the golden section step into the larger of the two segments.
				e=merge(a-x,b-x, x >= xm ) 
				d=CGOLD*e
			end if
			u=merge(x+d,x+sign(tol1,d), abs(d) >= tol1 )
			!Arrive here with d computed either from parabolic fit, or else from golden section.
			fu=func(u)
			!This is the one function evaluation per iteration.
			if (fu <= fx) then !Now we have to decide what to do with our function evaluation. Housekeeping follows:
				if (u >= x) then 
					a=x
					else
					b=x
				end if
				call shft(v,w,x,u)
				call shft(fv,fw,fx,fu)
			else
				if (u < x) then
					a=u
					else
					b=u
				end if
				if (fu <= fw .or. w == x) then
					v=w
					fv=fw
					w=u
					fw=fu
				else if (fu <= fv .or. v == x .or. v == w) then
					v=u
					fv=fu
				end if
			end if
		end do !Done with housekeeping. Back for another iteration.
		call nrerror("brent:exceed maximum iterations") 
		CONTAINS
			SUBROUTINE shft(a,b,c,d)
				REAL(DP), INTENT(OUT) :: a
				REAL(DP), INTENT(INOUT) :: b,c
				REAL(DP), INTENT(IN) :: d
				a=b
				b=c
				c=d
			END SUBROUTINE shft
	END FUNCTION brent

FUNCTION golden(ax,bx,func,xmax)
USE nrtype
! Golden search
	IMPLICIT NONE
	REAL(DP), INTENT(IN) :: ax,bx
	REAL(DP), INTENT(OUT) :: xmax
	REAL(DP), PARAMETER :: TOL=epsilon(1.0_DP)**(1.0/2.0)
	REAL(DP) :: golden
	INTERFACE
		FUNCTION func(x)
		USE nrtype
		IMPLICIT NONE
		REAL(DP), INTENT(IN) :: x
		REAL(DP) :: func
		END FUNCTION func
	END INTERFACE
	REAL(DP), PARAMETER :: R=0.618033988749895_DP,C=1.0_DP-R
	REAL(DP) :: f1,f2,x0,x1,x2,x3

	x0=ax
	x3=bx
	x1=R*ax+C*bx
	x2=C*ax+R*bx
	f1=func(x1)
	f2=func(x2)
	DO
		IF (abs(x3-x0) <= TOL*(abs(x1)+abs(x2))) EXIT
		IF (f2 > f1) THEN
			x0=x1
			x1=x2
			x2=R*x2+C*x3
			f1=f2
			f2=func(x2)
		ELSE
			x3=x2
			x2=x1
			x1=R*x1+C*x0
			f2=f1
			f1=func(x1)
		END IF
	END DO
	IF (f1 > f2) THEN
		golden=f1
		xmax=x1
	ELSE
		golden=f2
		xmax=x2
	END IF
END FUNCTION golden
	!-----------------------------------------------------------------------------------------------
	! FUNCTIONS USED TO CALCULATE CDF OF NORMAL DISTRIBUTION
	!-----------------------------------------------------------------------------------------------
	FUNCTION erf_s(x)
		USE nrtype
		!USE nr, ONLY : gammp
		IMPLICIT NONE
		REAL(DP), INTENT(IN) :: x
		REAL(DP) :: erf_s
		!Returns the error function erf(x).
		erf_s=gammp(0.5_dp,x**2)
		if (x < 0.0) erf_s=-erf_s
	END FUNCTION erf_s
	!-------------------------------------------------------------------------------------------------
	FUNCTION gammp_s(a,x)
		USE nrtype; USE nrutil, ONLY : assert
		!USE nr, ONLY : gcf,gser
		IMPLICIT NONE
		REAL(DP), INTENT(IN) :: a,x
		REAL(DP) :: gammp_s
		!Returns the incomplete gamma function P(a,x).
		call assert( x >= 0.0, a > 0.0, "gammp_s args")
		if (x<a+1.0_dp) then !Use the series representation.
			gammp_s=gser(a,x)
		else !Use the continued fraction representation and take its complement.
			gammp_s=1.0_dp-gcf(a,x) 
		end if
	END FUNCTION gammp_s
	!------------------------------------------------------------------------------------------------
	FUNCTION gammp_v(a,x)
		USE nrtype; USE nrutil, ONLY : assert,assert_eq
		!USE nr, ONLY : gcf,gser
		IMPLICIT NONE
		REAL(DP), DIMENSION(:), INTENT(IN) :: a,x
		REAL(DP), DIMENSION(size(x)) :: gammp_v
		LOGICAL(LGT), DIMENSION(size(x)) :: mask
		INTEGER(I4B) :: ndum
		ndum=assert_eq(size(a),size(x),"gammp_v")
		call assert( all(x >= 0.0), all(a > 0.0), "gammp_vargs")
		mask = (x<a+1.0_dp)
		gammp_v=merge(gser(a,merge(x,0.0_dp,mask)), &
		1.0_dp-gcf(a,merge(x,0.0_dp,.not. mask)),mask)
	END FUNCTION gammp_v
	!-------------------------------------------------------------------------------------------------
	FUNCTION gcf_s(a,x,gln)
		USE nrtype; USE nrutil, ONLY : nrerror
		!USE nr, ONLY : gammln
		IMPLICIT NONE
		REAL(DP), INTENT(IN) :: a,x
		REAL(DP), OPTIONAL, INTENT(OUT) :: gln
		REAL(DP) :: gcf_s
		INTEGER(I4B), PARAMETER :: ITMAX=100
		REAL(DP), PARAMETER :: EPS=epsilon(x),FPMIN=tiny(x)/EPS
	
		!Returns the incomplete gamma function Q(a,x) evaluated by its continued fraction repre-
		!sentation as gammcf. Also optionally returns lnG(a) as gln.
		!Parameters: ITMAX is the maximum allowed number of iterations; EPS is the relative accu-
		!racy; FPMIN is a number near the smallest representable ?oating-point number.
		
		INTEGER(I4B) :: i
		REAL(DP) :: an,b,c,d,del,h
		if (x == 0.0) then
			gcf_s=1.0
			RETURN
		end if
		b=x+1.0_dp-a !Set up for evaluating continued fraction by modified Lentz’s method (§5.2) with b0 =0.
		c=1.0_dp/FPMIN 
		d=1.0_dp/b
		h=d
		do i=1,ITMAX !Iterate to convergence.
			an=-i*(i-a)
			b=b+2.0_dp
			d=an*d+b
			if (abs(d) < FPMIN) d=FPMIN
			c=b+an/c
			if (abs(c) < FPMIN) c=FPMIN
			d=1.0_dp/d
			del=d*c
			h=h*del
			if (abs(del-1.0_dp) <= EPS) exit
		end do
		if (i > ITMAX) call nrerror("a too large, ITMAX too small in gcf_s")
		if (present(gln)) then
			gln=gammln(a)
			gcf_s=exp(-x+a*log(x)-gln)*h !Put factors in front.
		else
			gcf_s=exp(-x+a*log(x)-gammln(a))*h
		end if
	END FUNCTION gcf_s
	!-------------------------------------------------------------------------------------------------
	FUNCTION gcf_v(a,x,gln)
		USE nrtype; USE nrutil, ONLY : assert_eq,nrerror
		!USE nr, ONLY : gammln
		IMPLICIT NONE
		REAL(DP), DIMENSION(:), INTENT(IN) :: a,x
		REAL(DP), DIMENSION(:), OPTIONAL, INTENT(OUT) :: gln
		REAL(DP), DIMENSION(size(a)) :: gcf_v
		INTEGER(I4B), PARAMETER :: ITMAX=100
		REAL(DP), PARAMETER :: EPS=epsilon(x),FPMIN=tiny(x)/EPS
		INTEGER(I4B) :: i
		REAL(DP), DIMENSION(size(a)) :: an,b,c,d,del,h
		LOGICAL(LGT), DIMENSION(size(a)) :: converged,zero
		i=assert_eq(size(a),size(x),"gcf_v")
		zero=(x == 0.0)
		where (zero)
			gcf_v=1.0
		elsewhere
			b=x+1.0_dp-a
			c=1.0_dp/FPMIN
			d=1.0_dp/b
			h=d
		end where
		converged=zero
		do i=1,ITMAX
			where (.not. converged)
				an=-i*(i-a)
				b=b+2.0_dp
				d=an*d+b
				d=merge(FPMIN,d, abs(d)<FPMIN )
				c=b+an/c
				c=merge(FPMIN,c, abs(c)<FPMIN )
				d=1.0_dp/d
				del=d*c
				h=h*del
				converged = (abs(del-1.0_dp)<=EPS)
			end where
			if (all(converged)) exit
		end do
		if (i > ITMAX) call nrerror("a too large, ITMAX too small in gcf_v")
		if (present(gln)) then
			if (size(gln) < size(a)) call nrerror("gser: Not enough space for gln")
			gln=gammln(a)
			where (.not. zero) gcf_v=exp(-x+a*log(x)-gln)*h
		else
			where (.not. zero) gcf_v=exp(-x+a*log(x)-gammln(a))*h
		end if
	END FUNCTION gcf_v
	!-------------------------------------------------------------------------------------------------
	FUNCTION gser_s(a,x,gln)
		USE nrtype; USE nrutil, ONLY : nrerror
		!USE nr, ONLY : gammln
		IMPLICIT NONE
		REAL(DP), INTENT(IN) :: a,x
		REAL(DP), OPTIONAL, INTENT(OUT) :: gln
		REAL(DP) :: gser_s
		INTEGER(I4B), PARAMETER :: ITMAX=100
		REAL(DP), PARAMETER :: EPS=epsilon(x)
		INTEGER(I4B) :: n
		REAL(DP) :: ap,del,summ
		
		if (x == 0.0) then
			gser_s=0.0
			RETURN
		end if
		ap=a
		summ=1.0_dp/a
		del=summ
		do n=1,ITMAX
			ap=ap+1.0_dp
			del=del*x/ap
			summ=summ+del
			if (abs(del) < abs(summ)*EPS) exit
		end do
		if (n > ITMAX) call nrerror("a too large, ITMAX too small in gser_s")
		if (present(gln)) then
			gln=gammln(a)
			gser_s=summ*exp(-x+a*log(x)-gln)
		else
			gser_s=summ*exp(-x+a*log(x)-gammln(a))
		end if
	END FUNCTION gser_s
	!------------------------------------------------------------------------------------------------
	FUNCTION gser_v(a,x,gln)
		USE nrtype; USE nrutil, ONLY : assert_eq,nrerror
		!USE nr, ONLY : gammln
		IMPLICIT NONE
		REAL(DP), DIMENSION(:), INTENT(IN) :: a,x
		REAL(DP), DIMENSION(:), OPTIONAL, INTENT(OUT) :: gln
		REAL(DP), DIMENSION(size(a)) :: gser_v
		INTEGER(I4B), PARAMETER :: ITMAX=100
		REAL(DP), PARAMETER :: EPS=epsilon(x)
		INTEGER(I4B) :: n
		REAL(DP), DIMENSION(size(a)) :: ap,del,summ
		LOGICAL(LGT), DIMENSION(size(a)) :: converged,zero
		n=assert_eq(size(a),size(x),"gser_v")
		zero=(x == 0.0)
		where (zero) gser_v=0.0
		ap=a
		summ=1.0_dp/a
		del=summ
		converged=zero
		do n=1,ITMAX
			where (.not. converged)
				ap=ap+1.0_dp
				del=del*x/ap
				summ=summ+del
				converged = (abs(del) < abs(summ)*EPS)
			end where
			if (all(converged)) exit
		end do
		if (n > ITMAX) call nrerror("a too large, ITMAX too small in gser_v")
		if (present(gln)) then
			if (size(gln) < size(a)) call nrerror("gser: Not enough space for gln")
			gln=gammln(a)
			where (.not. zero) gser_v=summ*exp(-x+a*log(x)-gln)
		else
			where (.not. zero) gser_v=summ*exp(-x+a*log(x)-gammln(a))
		end if
	END FUNCTION gser_v
	!-----------------------------------------------------------------------------------------------
	FUNCTION gammln_s(xx)
		USE nrtype;
		USE nrutil, ONLY : arth_d,assert
		IMPLICIT NONE
		REAL(DP), INTENT(IN) :: xx
		REAL(DP) :: gammln_s
		!Returns the value ln[G(xx)] for xx > 0.
		REAL(DP) :: tmp,x
		REAL(DP) :: stp = 2.5066282746310005_dp
		REAL(DP), DIMENSION(6) :: coef = (/76.18009172947146_dp,&
		-86.50532032941677_dp,24.01409824083091_dp,&
		-1.231739572450155_dp,0.1208650973866179e-2_dp,&
		-0.5395239384953e-5_dp/)

		call assert(xx > 0.0, "gammln_s arg")
		x=xx
		tmp=x+5.5_dp
		tmp=(x+0.5_dp)*log(tmp)-tmp
		gammln_s=tmp+log(stp*(1.000000000190015_dp+&
		sum(coef(:)/arth_d(x+1.0_dp,1.0_dp,size(coef))))/x)
	END FUNCTION gammln_s
	!-----------------------------------------------------------------------------------------------
	FUNCTION gammln_v(xx)
		USE nrtype;
		USE nrutil, ONLY: assert
		IMPLICIT NONE
		INTEGER(I4B) :: i
		REAL(DP), DIMENSION(:), INTENT(IN) :: xx
		REAL(DP), DIMENSION(size(xx)) :: gammln_v
		REAL(DP), DIMENSION(size(xx)) :: ser,tmp,x,y
		REAL(DP) :: stp = 2.5066282746310005_dp
		REAL(DP), DIMENSION(6) :: coef = (/76.18009172947146_dp,&
		-86.50532032941677_dp,24.01409824083091_dp,&
		-1.231739572450155_dp,0.1208650973866179e-2_dp,&
		-0.5395239384953e-5_dp/)
		if (size(xx) == 0) RETURN
		call assert(all(xx > 0.0), "gammln_varg")
		x=xx
		tmp=x+5.5_dp
		tmp=(x+0.5_dp)*log(tmp)-tmp
		ser=1.000000000190015_dp
		y=x
		do i=1,size(coef)
			y=y+1.0_dp
			ser=ser+coef(i)/y
		end do
		gammln_v=tmp+log(stp*ser/x)
	END FUNCTION gammln_v
	!-----------------------------------------------------------------------------------------------

	!-----------------------------------------------------------------------------------------------
	! RANDOM NUMBER GENERATION AND MARKOV CHAIN REALIZATION ROUTINES
	!-----------------------------------------------------------------------------------------------
	FUNCTION ran(idum)
		IMPLICIT NONE
		INTEGER, PARAMETER :: K4B=selected_int_kind(9)
		INTEGER(K4B), INTENT(INOUT) :: idum
		REAL :: ran

		!“Minimal” random number generator of Park and Miller combined with a Marsaglia shift
		!sequence. Returns a uniform random deviate between 0.0 and 1.0 (exclusive of the endpoint
		!values). This fully portable, scalar generator has the “traditional” (not Fortran 90) calling
		!sequence with a random deviate as the returned function value: call with idum a negative
		!integer to initialize; thereafter, do not alter idum except to reinitialize. The period of this
		!generator is about 3.1 × 10^18.

		INTEGER(K4B), PARAMETER :: IA=16807,IM=2147483647,IQ=127773,IR=2836
		REAL, SAVE :: am
		INTEGER(K4B), SAVE :: ix=-1,iy=-1,k
		if (idum <= 0 .or. iy < 0) then !Initialize.
			am=nearest(1.0,-1.0)/IM
			iy=ior(ieor(888889999,abs(idum)),1)
			ix=ieor(777755555,abs(idum))
			idum=abs(idum)+1 !Set idum positive.
		end if
		ix=ieor(ix,ishft(ix,13)) !Marsaglia shift sequence with period 2^32 - 1.
		ix=ieor(ix,ishft(ix,-17))
		ix=ieor(ix,ishft(ix,5))
		k=iy/IQ !Park-Miller sequence by Schrage’s method, period 2^31 - 2.
		iy=IA*(iy-k*IQ)-IR*k 
		if (iy < 0) iy=iy+IM
		ran=am*ior(iand(IM,ieor(ix,iy)),1) !Combine the two generators with masking to ensure nonzero value.
	END FUNCTION ran
	!--------------------------------------------------------------------------------------------------------
	SUBROUTINE ran1_s(harvest)
		USE nrtype
		USE ran_state, ONLY: K4B,amm,lenran,ran_init, &
		iran0,jran0,kran0,nran0,mran0,rans
		IMPLICIT NONE
		REAL(SP), INTENT(OUT) :: harvest
		!Lagged Fibonacci generator combined with two Marsaglia shift sequences. On output, re-
		!turns as harvest a uniform random deviate between 0.0 and 1.0 (exclusive of the endpoint
		!values). This generator has the same calling and initialization conventions as Fortran 90’s
		!random number routine. Use ran seed to initialize or reinitialize to a particular sequence.
		!The period of this generator is about 8.5×10^37 , and it fully vectorizes. Validity of the integer
		!model assumed by this generator is tested at initialization.

		if (lenran < 1) call ran_init(1) !Initialization routine in ran_state.
		rans=iran0-kran0 !Update Fibonacci generator, which
		if (rans < 0) rans=rans+2147483579_k4b !has period p^2+p+1, p=2^31-69
		iran0=jran0
		jran0=kran0
		kran0=rans
		nran0=ieor(nran0,ishft(nran0,13)) !Update Marsaglia shift sequence.
		nran0=ieor(nran0,ishft(nran0,-17))
		nran0=ieor(nran0,ishft(nran0,5))
		!Once only per cycle, advance sequence by 1, shortening its period to 2^32 - 2.
		if (nran0 == 1) nran0=270369_k4b
		mran0=ieor(mran0,ishft(mran0,5)) !Update Marsaglia shift sequence with period 2^32 - 1.
		mran0=ieor(mran0,ishft(mran0,-13)) 
		mran0=ieor(mran0,ishft(mran0,6))
		rans=ieor(nran0,rans)+mran0
		!Combine the generators. The above statement has wrap-around addition.
		harvest=amm*merge(rans,not(rans), rans<0 ) !Make the result positive definite (note that amm is negative).
	END SUBROUTINE ran1_s
	!--------------------------------------------------------------------------------------------------------
	SUBROUTINE ran1_v(harvest)
		USE nrtype
		USE ran_state, ONLY: K4B,amm,lenran,ran_init, &
		iran,jran,kran,nran,mran,ranv
		IMPLICIT NONE
		REAL(DP), DIMENSION(:), INTENT(OUT) :: harvest
		INTEGER(K4B) :: n
		n=size(harvest)
		if (lenran < n+1) call ran_init(n+1)
		ranv(1:n)=iran(1:n)-kran(1:n)
		where (ranv(1:n) < 0) ranv(1:n)=ranv(1:n)+2147483579_k4b
		iran(1:n)=jran(1:n)
		jran(1:n)=kran(1:n)
		kran(1:n)=ranv(1:n)
		nran(1:n)=ieor(nran(1:n),ishft(nran(1:n),13))
		nran(1:n)=ieor(nran(1:n),ishft(nran(1:n),-17))
		nran(1:n)=ieor(nran(1:n),ishft(nran(1:n),5))
		where (nran(1:n) == 1) nran(1:n)=270369_k4b
		mran(1:n)=ieor(mran(1:n),ishft(mran(1:n),5))
		mran(1:n)=ieor(mran(1:n),ishft(mran(1:n),-13))
		mran(1:n)=ieor(mran(1:n),ishft(mran(1:n),6))
		ranv(1:n)=ieor(nran(1:n),ranv(1:n))+mran(1:n)
		harvest=amm*merge(ranv(1:n),not(ranv(1:n)), ranv(1:n)<0 )
	END SUBROUTINE ran1_v
	!--------------------------------------------------------------------------------------------------------
	SUBROUTINE gasdev_s(harvest)
		USE nrtype
		!USE nr, ONLY : ran1
		IMPLICIT NONE
		REAL(SP), INTENT(OUT) :: harvest

		!Returns in harvest a normally distributed deviate with zero mean and unit variance, using
		!ran1 as the source of uniform deviates.

		REAL(SP) :: rsq,v1,v2
		REAL(SP), SAVE :: g
		LOGICAL, SAVE :: gaus_stored=.false.
		if (gaus_stored) then !We have an extra deviate handy,
		harvest=g !so return it,
		gaus_stored=.false. !and unset the flag.
		else !We don’t have an extra deviate handy, so
		do
		call ran1(v1) !pick two uniform numbers in the square ex-
		call ran1(v2) !tending from -1 to +1 in each direction,
		v1=2.0_sp*v1-1.0_sp
		v2=2.0_sp*v2-1.0_sp
		rsq=v1**2+v2**2 !see if they are in the unit circle,
		if (rsq > 0.0 .and. rsq < 1.0) exit
		end do !otherwise try again.
		rsq=sqrt(-2.0_sp*log(rsq)/rsq) !Now make the Box-Muller transformation to
		harvest=v1*rsq !get two normal deviates. Return one and
		g=v2*rsq !save the other for next time.
		gaus_stored=.true. !Set flag.
		end if
	END SUBROUTINE gasdev_s
	!--------------------------------------------------------------------------------------------------------
	SUBROUTINE gasdev_v(harvest)
		USE nrtype; USE nrutil, ONLY : array_copy
		!USE nr, ONLY : ran1
		IMPLICIT NONE
		REAL(DP), DIMENSION(:), INTENT(OUT) :: harvest
		REAL(DP), DIMENSION(size(harvest)) :: rsq,v1,v2
		REAL(DP), ALLOCATABLE, DIMENSION(:), SAVE :: g
		INTEGER(I4B) :: n,ng,nn,m
		INTEGER(I4B), SAVE :: last_allocated=0
		LOGICAL, SAVE :: gaus_stored=.false.
		LOGICAL, DIMENSION(size(harvest)) :: mask
		n=size(harvest)
		if (n /= last_allocated) then
		if (last_allocated /= 0) deallocate(g)
		allocate(g(n))
		last_allocated=n
		gaus_stored=.false.
		end if
		if (gaus_stored) then
		harvest=g
		gaus_stored=.false.
		else
		ng=1
		do
		if (ng > n) exit
		call ran1(v1(ng:n))
		call ran1(v2(ng:n))
		v1(ng:n)=2.0_sp*v1(ng:n)-1.0_dp
		v2(ng:n)=2.0_sp*v2(ng:n)-1.0_dp
		rsq(ng:n)=v1(ng:n)**2+v2(ng:n)**2
		mask(ng:n)=(rsq(ng:n)>0.0 .and. rsq(ng:n)<1.0)
		call array_copy(pack(v1(ng:n),mask(ng:n)),v1(ng:),nn,m)
		v2(ng:ng+nn-1)=pack(v2(ng:n),mask(ng:n))
		rsq(ng:ng+nn-1)=pack(rsq(ng:n),mask(ng:n))
		ng=ng+nn
		end do
		rsq=sqrt(-2.0_sp*log(rsq)/rsq)
		harvest=v1*rsq
		g=v2*rsq
		gaus_stored=.true.
		end if
	END SUBROUTINE gasdev_v
	!--------------------------------------------------------------------------------------------------------
	subroutine markovgen(transmat,out,istart)
		use nrtype
		use nrutil, only: assert_eq, nrerror
		implicit none

		real(DP), dimension(:,:), intent(in) :: transmat
		integer, dimension(:), intent(inout) :: out
		integer, intent(in) :: istart
		!integer, intent(inout) :: iseed

		! generates a vector of realizations of an m-state markov chain represented by transmat.
		! the vector out is filled with integers in the range 1-m.
		! the starting state istart must be in that range as well.
		! iseed is an integer seed for the random number generator
		! based on the routine of the same name in NR 3rd edition C++ code

		integer :: i, ilo, ihi, ii, j, m, n
		real(DP), dimension(size(transmat,dim=1),size(transmat,dim=2)) :: cum
		real(DP), dimension(size(out)) :: rr
		real(DP) :: r

		m=assert_eq(size(transmat,dim=1),size(transmat,dim=2),"markovgen")
		if (istart<0 .or. istart>m) call nrerror("markovgen: istart outside range of state indices!")
		n=size(out)

		cum=transmat
		do i=1,m
			do j=2,m
				cum(i,j)=cum(i,j)+cum(i,j-1)
			enddo
			if (abs(cum(i,m)-1)>0.01) call nrerror("markovgen: transmat not a transition matrix!")
		enddo

		call ran1_v(rr)

		j=istart
		out(1)=j
		do ii=2,n
			r=rr(ii)/cum(j,m)
			ilo=1
			ihi=m+1
			do
				if (ihi-ilo<=1) exit
				i=ishft(ihi+ilo,-1)
				if (r>cum(j,i-1)) then
					ilo=i
				else
					ihi=i
				endif
			enddo
			j=ilo
			out(ii)=j
		enddo

	end subroutine markovgen
	!-----------------------------------------------------------------------------------------------------------

	!-----------------------------------------------------------------------------------------------------------
	! STATISTICAL DESCRIPTION OF DATA
	!-----------------------------------------------------------------------------------------------------------
	SUBROUTINE moment(data,ave,adev,sdev,var,skew,curt)
		USE nrtype; USE nrutil, ONLY : nrerror
		IMPLICIT NONE
		REAL(DP), INTENT(OUT) :: ave,adev,sdev,var,skew,curt
		REAL(DP), DIMENSION(:), INTENT(IN) :: data

		!Given an array of data, this routine returns its mean ave, average deviation adev, standard
		!deviation sdev,variance var, skewness skew, and kurtosis curt.

		INTEGER(I4B) :: n
		REAL(DP) :: ep
		REAL(DP), DIMENSION(size(data)) :: p,s
		n=size(data)
		if (n <= 1) call nrerror("moment: n must be at least 2")
		ave=sum(data(:))/n !First pass to get the mean.
		s(:)=data(:)-ave !Second pass to get the ?rst (absolute), second, third, and
		ep=sum(s(:)) !fourth moments of the deviation from the mean.
		adev=sum(abs(s(:)))/n
		p(:)=s(:)*s(:)
		var=sum(p(:))
		p(:)=p(:)*s(:)
		skew=sum(p(:))
		p(:)=p(:)*s(:)
		curt=sum(p(:))
		var=(var-ep**2/n)/(n-1) !Corrected two-pass formula.
		sdev=sqrt(var)
		if (var /= 0.0) then
		skew=skew/(n*sdev**3)
		curt=curt/(n*var**2)-3.0_dp
		else
		call nrerror("moment: no skew or kurtosis when zero variance")
		end if
	END SUBROUTINE moment
	!----------------------------------------------------------------------------------------------------------
	SUBROUTINE avevar(data,ave,var)
		USE nrtype
		IMPLICIT NONE
		REAL(DP), DIMENSION(:), INTENT(IN) :: data
		REAL(DP), INTENT(OUT) :: ave,var

		!Given array data, returns its mean as ave and its variance as var.

		INTEGER(I4B) :: n
		REAL(DP), DIMENSION(size(data)) :: s
		n=size(data)
		ave=sum(data(:))/n
		s(:)=data(:)-ave
		var=dot_product(s,s)
		var=(var-sum(s)**2/n)/(n-1) !Corrected two-pass formula (14.1.8).
	END SUBROUTINE avevar
	!--------------------------------------------------------------------------------------------------------------

!******************************************************************************************!
!(3) CHEBYSHEV INTERPOLATION
!******************************************************************************************!

	!********************************************************************************************!
	!Function: cheb_nodes
	!	This Function calculates  m Chebyshev interpolation nodes for an interval [a,b] 
	!								
	!	Inputs:	
	!			m = # of nodes
	!			a = min interval
	!			b = max interval
	!	Outputs
	!			cheb_nodes(m) = Chebyshev nodes
	!********************************************************************************************!

Function cheb_nodes(mcheb,acheb,bcheb)

	!------------------------------------------------
		integer, intent(in) :: mcheb
		real(DP), intent(in) :: acheb, bcheb
		real(DP), dimension(mcheb) :: cheb_nodes
		integer :: i
		real(DP), dimension(mcheb) :: x
		real(DP) :: mm
	!------------------------------------------------

mm = mcheb*1.0_dp

	do i=1,mcheb
		x(i) =  cos((2.0*i -1.0)*3.14/(2.0*mcheb))	
		cheb_nodes(i) = acheb + (x(i) + 1.0)*(bcheb-acheb)/2.0	
	end do


end function cheb_nodes

	!********************************************************************************************

	!********************************************************************************************!
	!Function: cheb_weight(m)
	!	This Function calculate weights for Chebyshev integration of a normal variable with m nodes 
	!								
	!	Inputs:	
	!			m = number of nodes
	!
	!			
	!	Outputs
	!			cheb_weight= w(x)
	!********************************************************************************************!

Function cheb_weight(m_in)
!------------------------------------------------
		integer, intent(in) :: m_in
		real(DP),  dimension(m_in) :: cheb_weight
		integer :: i
		real(DP),  dimension(m_in) :: node
!------------------------------------------------
node = cheb_nodes(m_in,-1.0_dp,1.0_dp)

	do i =1,m_in
	    cheb_weight(i) = (1-node(i)**2)**(0.5)
	end do

end function cheb_weight

	!********************************************************************************************

	!********************************************************************************************!
	!Function: cheb_coeff
	!	This Function calculates n Chebyshev coefficients given m interpolation points and ygrid
	!								
	!	Inputs:	
	!			n = degree of polynomial
	!			m = # of nodes
	!			ygrid = known points
	!	Outputs
	!			cheb_coeff(n+1) = Chebyshev coefficients
	!********************************************************************************************!

Function cheb_coeff(ncheb,mcheb,yygrid)
	!------------------------------------------------
		integer, intent(in) :: ncheb, mcheb
		real(DP), dimension(mcheb), intent(in) :: yygrid
		real(DP),dimension(ncheb) :: cheb_coeff
		integer :: i, j, n, m
		real(DP), dimension(mcheb) :: x
		real(DP),dimension(mcheb,ncheb) :: t		!one vector of basis (dim(n+1)) for any m interpolation nodes
		real(DP) :: tt, anorm, bnorm
	!------------------------------------------------

! starting:
	anorm = -1.0
	bnorm = 1.0
	m = mcheb
	n = ncheb
	x = cheb_nodes(m,anorm,bnorm)		!note that I'm using the normalized interval [-1,1]
	
!calculate n coeff of the basis function for each node
	do i = 1,m
		t(i,:) = cheb_basis(x(i),ncheb)		
	end do

	do j = 1,n
		cheb_coeff(j) = 0		!initializing k(j)
		tt = 0					!initializing k(j)
		do i = 1,m
			cheb_coeff(j) = cheb_coeff(j) + yygrid(i)*t(i,j)
			tt = tt + t(i,j)**2
		end do
	
	cheb_coeff(j) = cheb_coeff(j)/tt

	end do


end function cheb_coeff
	!********************************************************************************************

	!********************************************************************************************!
	!Function: cheb_basis(x,n)
	!	This Function evaluates the basis function for Chebyshev polynomials of order n evaluated 
	!			at x on the interval [-1,1]
	!								
	!	Inputs:	
	!			n = degree of polynomial
	!			x = Eval @ F(x)
	!			
	!	Outputs
	!			cheb_basis = F(x)
	!********************************************************************************************!

Function cheb_basis(xin,nin)
!------------------------------------------------
		real(DP), intent(in) :: xin
		integer, intent(in) :: nin
		real(DP),  dimension(nin) :: cheb_basis
		integer :: i
!------------------------------------------------

cheb_basis(1) = 1;
cheb_basis(2) = xin;

IF(nin>2) then
	do i =3,nin
	    cheb_basis(i) = 2*xin*cheb_basis(i-1)-cheb_basis(i-2)
	end do
EndIF
end function cheb_basis

	!********************************************************************************************

	!********************************************************************************************!
	!Function: cheby(x,n,m,ygrid,min,max)
	!	This Function calculates F(x) by interpolation using chebyshev polynomials
	!								
	!	Inputs:	
	!			n = degree of polynomial
	!			m = # of nodes						**	NOTES: m>=n  **
	!			ygrid = known points
	!			x = Eval @ F(x)
	!			
	!	Outputs
	!			cheb_basis = F(x)
	!********************************************************************************************!

Function cheby(xin,nin,minn,yygrid,minin,maxin)
!------------------------------------------------
	real(DP), intent(in) :: xin, minin, maxin
	integer,intent(in) :: nin, minn
	real(DP),intent(in), dimension(minn) :: yygrid
	real(DP):: cheby, fx
	real,dimension(nin) :: k, b
	integer :: i, m, n
	real(DP) :: xnorm, x, min, max
!------------------------------------------------
m=minn
n=nin
x=xin
min=minin
max=maxin

k = cheb_coeff(n,m,yygrid)

xnorm = (2*(x - min)/(max - min))-1

b = cheb_basis(xnorm,n)

fx = 0
	do i=1,n
		fx = fx + k(i)*b(i)
	end do
cheby = fx
end function cheby

	!********************************************************************************************!

	!********************************************************************************************

	!********************************************************************************************!
	!Function: linspace(xmin,xmax,np)
	!	This Function creates a linearly spaced grid
	!								
	!	Inputs:	
	!			xmin = minimum
	!			xmax = maximum
	!			np = number of points
	!			
	!	Outputs
	!			linspace = 1 x np vector from xmin to xmax
	!********************************************************************************************!

Function lingrid(xmin,xmax,np)
!------------------------------------------------
	real(DP), intent(in) :: xmin, xmax
	integer,intent(in) :: np
	real(DP), dimension(np) :: lingrid
	real(DP)	::	d
	integer :: i
!------------------------------------------------

d = (xmax-xmin)/np
lingrid(1) = xmin

DO i = 2,np
	lingrid(i) = d + lingrid(i-1) 
endDO

end function lingrid

	!********************************************************************************************!


	!********************************************************************************************
	!********************************************************************************************!
	!Function: nrmlpdf(mu1,sig1,Z1)
	!	This Function calculates  normal pdf
	!********************************************************************************************!

Function mynrmlpdf(mu2,sig2,Z2)
!------------------------------------------------
	real(DP), intent(in) :: mu2, sig2, Z2
	real(DP)	::	mynrmlpdf
!------------------------------------------------


mynrmlpdf = (2.71828183**(-((Z2-mu2)**2.0)/(2.0*sig2)))*((2.0*3.14159265*sig2)**-0.5)


end function mynrmlpdf

	!********************************************************************************************!
	!Function: nrmlcdf(aa,bb,mu1,sig1,Z1)
	!	This Function calculates  normal cdf
	!********************************************************************************************!

Function nrmlcdf(mu2,sig2,Z2)
!------------------------------------------------
	real(DP), intent(in) :: mu2, sig2, Z2
	real(DP)	::	nrmlcdf
!------------------------------------------------


nrmlcdf = 0.5*(1+erf_s((Z2-mu2)/((2*(sig2))**(0.5))))


end function nrmlcdf

	!********************************************************************************************!


	!********************************************************************************************!
	!Function: tr_nrmlcdf(aa,bb,mu1,sig1,Z1)
	!	This Function calculates truncated normal cdf
	!********************************************************************************************!

Function tr_nrmlcdf(aa,bb,mu1,sig1,Z1)
!------------------------------------------------
	real(DP), intent(in) :: aa, bb, mu1, sig1, Z1
	real(DP)	::	tr_nrmlcdf, z1x, z1a, z1b
!------------------------------------------------


z1x = (Z1-mu1)/sig1
z1a = (aa-mu1)/sig1
z1b = (bb-mu1)/sig1

tr_nrmlcdf = (nrmlcdf(mu1,sig1,z1x)-nrmlcdf(mu1,sig1,z1a))/(nrmlcdf(mu1,sig1,z1b)-nrmlcdf(mu1,sig1,z1a))

end function tr_nrmlcdf

	!********************************************************************************************!

	!********************************************************************************************!
	!Gaussian-Hermite for expectation of normal dist
	! E[f(y)] = pi^(-0.5)\sum w_i f(sqrt(2)*sig*x+mu)
	!		Where y ~ N(mu,sig)
	!********************************************************************************************!

	SUBROUTINE gauher(x,w)
	USE nrtype; USE nrutil, ONLY : arth,assert_eq,nrerror
	IMPLICIT NONE
	REAL(DP), DIMENSION(:), INTENT(OUT) :: x,w
	REAL(DP), PARAMETER :: EPS=3.0e-13_dp,PIM4=0.7511255444649425_dp
	INTEGER(I4B) :: its,j,m,n
	INTEGER(I4B), PARAMETER :: MAXIT=10
	REAL(DP) :: anu
	REAL(DP), PARAMETER :: C1=9.084064e-01_dp,C2=5.214976e-02_dp,&
		C3=2.579930e-03_dp,C4=3.986126e-03_dp
	REAL(DP), DIMENSION((size(x)+1)/2) :: rhs,r2,r3,theta
	REAL(DP), DIMENSION((size(x)+1)/2) :: p1,p2,p3,pp,z,z1
	LOGICAL(LGT), DIMENSION((size(x)+1)/2) :: unfinished
	n=assert_eq(size(x),size(w),'gauher')
	m=(n+1)/2
	anu=2.0_dp*n+1.0_dp
	rhs=arth(3,4,m)*PI/anu
	r3=rhs**(1.0_sp/3.0_dp)
	r2=r3**2
	theta=r3*(C1+r2*(C2+r2*(C3+r2*C4)))
	z=sqrt(anu)*cos(theta)
	unfinished=.true.
	do its=1,MAXIT
		where (unfinished)
			p1=PIM4
			p2=0.0
		end where
		do j=1,n
			where (unfinished)
				p3=p2
				p2=p1
				p1=z*sqrt(2.0_dp/j)*p2-sqrt(real(j-1,dp)/real(j,dp))*p3
			end where
		end do
		where (unfinished)
			pp=sqrt(2.0_dp*n)*p2
			z1=z
			z=z1-p1/pp
			unfinished=(abs(z-z1) > EPS)
		end where
		if (.not. any(unfinished)) exit
	end do
	if (its == MAXIT+1) call nrerror('too many iterations in gauher')
	x(1:m)=z
	x(n:n-m+1:-1)=-z
	w(1:m)=2.0_dp/pp**2
	w(n:n-m+1:-1)=w(1:m)
	END SUBROUTINE gauher

		!********************************************************************************************!


	!********************************************************************************************!

	!********************************************************************************************!
	!Inverse Erf function
	! 
	!********************************************************************************************!
       function dinvnorm(p)
      real(DP), intent(in) :: p
      real(DP) :: a1,a2,a3,a4,a5,a6, p_low,p_high, b1,b2,b3,b4,b5, c1,c2,c3,c4,c5,c6, d1,d2,d3,d4, z,q,r, dinvnorm  

      a1=-39.6968302866538
      a2=220.946098424521
      a3=-275.928510446969
      a4=138.357751867269
      a5=-30.6647980661472
      a6=2.50662827745924
      b1=-54.4760987982241
      b2=161.585836858041
      b3=-155.698979859887
      b4=66.8013118877197
      b5=-13.2806815528857
      c1=-0.00778489400243029
      c2=-0.322396458041136
      c3=-2.40075827716184
      c4=-2.54973253934373
      c5=4.37466414146497
      c6=2.93816398269878
      d1=0.00778469570904146
      d2=0.32246712907004
      d3=2.445134137143
      d4=3.75440866190742
      p_low=0.02425
      p_high=1-p_low

      if(p < p_low) then 

		q=dsqrt(-2*dlog(p))
		z=(((((c1*q+c2)*q+c3)*q+c4)*q+c5)*q+c6)/((((d1*q+d2)*q+d3)*q+d4)*q+1)
	  else
		if ((p > p_low).and.(p < p_high)) then
		  q=p-0.5
		  r=q*q
		  z=(((((a1*r+a2)*r+a3)*r+a4)*r+a5)*r+a6)*q/(((((b1*r+b2)*r+b3)*r+b4)*r+b5)*r+1)
		else
			q=dsqrt(-2*dlog(1-p))
		    z=-(((((c1*q+c2)*q+c3)*q+c4)*q+c5)*q+c6)/((((d1*q+d2)*q+d3)*q+d4)*q+1)
		endif
	  endif
dinvnorm=z

endfunction



Subroutine D_inssor (XDONT)
!  Sorts XDONT into increasing order (Insertion sort)
! __________________________________________________________
!  This subroutine uses insertion sort. It does not use any
!  work array and is faster when XDONT is of very small size
!  (< 20), or already almost sorted, but worst case behavior
!  can happen fairly probably (initially inverse sorted).
!  In many cases, the quicksort or merge sort method is faster.
!  Michel Olagnon - Apr. 2000
! __________________________________________________________
! __________________________________________________________
! __________________________________________________________
      Real (DP), Dimension (:), Intent (InOut) :: XDONT
! __________________________________________________________
      Real (DP) :: XWRK, XMIN
!
! __________________________________________________________
!
      Integer :: ICRS, IDCR, NDON
!
      NDON = Size (XDONT)
!
! We first bring the minimum to the first location in the array.
! That way, we will have a "guard", and when looking for the
! right place to insert a value, no loop test is necessary.
!
      If (XDONT (1) < XDONT (NDON)) Then
          XMIN = XDONT (1)
      Else
          XMIN = XDONT (NDON)
          XDONT (NDON) = XDONT (1)
      Endif
      Do IDCR = NDON-1, 2, -1
         XWRK = XDONT(IDCR)
         IF (XWRK < XMIN) Then
            XDONT (IDCR) = XMIN
            XMIN = XWRK
         End If
      End Do
      XDONT (1) = XMIN
!
! The first value is now the minimum
! Loop over the array, and when a value is smaller than
! the previous one, loop down to insert it at its right place.
!
      Do ICRS = 3, NDON
         XWRK = XDONT (ICRS)
         IDCR = ICRS - 1
         If (XWRK < XDONT(IDCR)) Then
            XDONT (ICRS) = XDONT (IDCR)
            IDCR = IDCR - 1
            Do
               If (XWRK >= XDONT(IDCR)) Exit
               XDONT (IDCR+1) = XDONT (IDCR)
               IDCR = IDCR - 1
            End Do
            XDONT (IDCR+1) = XWRK
         End If
      End Do
!
      Return
!
End Subroutine D_inssor



	
end module nr_lib
