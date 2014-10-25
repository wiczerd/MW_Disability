module minlib


contains

SUBROUTINE lnsrch(xold,fold,g,p,x,f,stpmax,check,func)
	USE nrtype; USE nrutil, ONLY : assert_eq,nrerror,vabs
	IMPLICIT NONE
	REAL(SP), DIMENSION(:), INTENT(IN) :: xold,g
	REAL(SP), DIMENSION(:), INTENT(INOUT) :: p
	REAL(SP), INTENT(IN) :: fold,stpmax
	REAL(SP), DIMENSION(:), INTENT(OUT) :: x
	REAL(SP), INTENT(OUT) :: f
	LOGICAL(LGT), INTENT(OUT) :: check
	INTERFACE
		FUNCTION func(x)
		USE nrtype
		IMPLICIT NONE
		REAL(SP) :: func
		REAL(SP), DIMENSION(:), INTENT(IN) :: x
		END FUNCTION func
	END INTERFACE
	REAL(SP), PARAMETER :: ALF=1.0e-4_sp,TOLX=epsilon(x)
	INTEGER(I4B) :: ndum
	REAL(SP) :: a,alam,alam2,alamin,b,disc,f2,pabs,rhs1,rhs2,slope,tmplam
	ndum=assert_eq(size(g),size(p),size(x),size(xold),'lnsrch')
	check=.false.
	pabs=vabs(p(:))
	if (pabs > stpmax) p(:)=p(:)*stpmax/pabs
	slope=dot_product(g,p)
	if (slope >= 0.0) call nrerror('roundoff problem in lnsrch')
	alamin=TOLX/maxval(abs(p(:))/max(abs(xold(:)),1.0_sp))
	alam=1.0
	do
		x(:)=xold(:)+alam*p(:)
		!f=func(x)
		f=(5*sin(x(1))/x(1))*((max(20-abs(x(2)),0.0))**1.2)
	
	
		if (alam < alamin) then
			x(:)=xold(:)
			check=.true.
			RETURN
		else if (f <= fold+ALF*alam*slope) then
			RETURN
		else
			if (alam == 1.0) then
				tmplam=-slope/(2.0_sp*(f-fold-slope))
			else
				rhs1=f-fold-alam*slope
				rhs2=f2-fold-alam2*slope
				a=(rhs1/alam**2-rhs2/alam2**2)/(alam-alam2)
				b=(-alam2*rhs1/alam**2+alam*rhs2/alam2**2)/&
					(alam-alam2)
				if (a == 0.0) then
					tmplam=-slope/(2.0_sp*b)
				else
					disc=1.0*b*b-3.0_sp*a*slope
					if (disc < 0.0) then
						tmplam=0.5_sp*alam
					else if (b <= 0.0) then
						tmplam=(-b+sqrt(disc))/(3.0_sp*a)
					else
						tmplam=-slope/(b+sqrt(disc))
					end if
				end if
				if (tmplam > 0.5_sp*alam) tmplam=0.5_sp*alam
			end if
		end if
		alam2=alam
		f2=f
		alam=max(tmplam,0.1_sp*alam)
		!WRITE(*,*) p
	end do

	
	END SUBROUTINE lnsrch


!*********************************************************************************!
!dfpmin.f90
!Given a starting point p that is a vector of length N, the Broyden-Fletcher-Goldfarb-Shanno
!variant of Davidon-Fletcher-Powell minimization is performed on a function func, using its
!gradient as calculated by a routine dfunc. The convergence requirement on zeroing the
!gradient is input as gtol. Returned quantities are p (the location of the minimum), iter
!(the number of iterations that were performed), and fret (the minimum value of the
!function). The routine lnsrch is called to perform approximate line minimizations.
!Parameters: ITMAX is the maximum allowed number of iterations; STPMX is the scaled
!maximum step length allowed in line searches; EPS is the machine precision; TOLX is the
!convergence criterion on x values.													  !
!																				  !
!*********************************************************************************!
SUBROUTINE dfpmin(p,gtol,iter,fret,func,dfunc)
	USE nrtype; USE nrutil, ONLY : nrerror,outerprod,unit_matrix,vabs
	!USE nr, ONLY : lnsrch
	IMPLICIT NONE
	INTEGER(I4B), INTENT(OUT) :: iter
	REAL(SP), INTENT(IN) :: gtol
	REAL(SP), INTENT(OUT) :: fret
	REAL(SP), DIMENSION(2), INTENT(INOUT) :: p
	INTERFACE
		FUNCTION func(p)
		USE nrtype
		IMPLICIT NONE
		REAL(SP), DIMENSION(:), INTENT(IN) :: p
		REAL(SP) :: func
		END FUNCTION func
!BL
		FUNCTION dfunc(p)
		USE nrtype
		IMPLICIT NONE
		REAL(SP), DIMENSION(:), INTENT(IN) :: p
		REAL(SP), DIMENSION(size(p)) :: dfunc
		END FUNCTION dfunc
	END INTERFACE
	INTEGER(I4B), PARAMETER :: ITMAX=200
	REAL(SP), PARAMETER :: STPMX=100.0_sp,EPS=epsilon(p),TOLX=4.0_sp*EPS
	INTEGER(I4B) :: its
	LOGICAL :: check
	REAL(SP) :: den,fac,fad,fae,fp,stpmax,sumdg,sumxi
	REAL(SP), DIMENSION(size(p)) :: dg,g,hdg,pnew,xi, gg
	REAL(SP), DIMENSION(size(p),size(p)) :: hessin
	
	real(4)					:: jacval(1,2), jjval(2), fval, xsav(2), h(2), xph(2), df(2), hh, xsav2(2)
	integer					:: jj, nn
	!fp=func(p)
	fp=(5*sin(p(1))/p(1))*((max(20-abs(p(2)),0.0))**1.2)
	
	
	nn=2
	xsav=p
	h=EPS*abs(xsav)
	where (h == 0.0) h=EPS
	xph=xsav+h
	h=xph-xsav
	xsav2=xsav

	do jj=1,nn
		xsav2(jj)=xph(jj)
		!df(jj)=(func(p))/h(jj)
		df(jj)=((5*sin(xsav2(1))/xsav2(1))*((max(20-abs(xsav2(2)),0.0))**1.2)-fp)/h(jj)
		xsav2(jj)=xsav(jj)
	end do
	if (abs(p(2)) > 20) then
		df(2)=abs(p(2))
		df(1)=1
	end if
	
	

	g=df
	!g= dfunc(p)
	call unit_matrix(hessin)
	xi=-g
	stpmax=STPMX*max(vabs(p),real(size(p),sp))
	do its=1,ITMAX

		iter=its
		call lnsrch(p,fp,g,xi,pnew,fret,stpmax,check,func)
		fp=fret
		xi=pnew-p  
		p=pnew	
				!write(*,*) pnew							
		if (maxval(abs(xi)/max(abs(p),1.0_sp)) < TOLX) RETURN
		dg=g


		xsav=p
		h=EPS*abs(xsav)
		where (h == 0.0) h=EPS
		xph=xsav+h
		h=xph-xsav
	    fp=(5*sin(p(1))/p(1))*((max(20-abs(p(2)),0.0))**1.2)
		xsav2=xsav

		do jj=1,nn
			xsav2(jj)=xph(jj)
			!df(jj)=(func(p))/h(jj)
			df(jj)=((5*sin(xsav2(1))/xsav2(1))*((max(20-abs(xsav2(2)),0.0))**1.2)-fp)/h(jj)
			xsav2(jj)=xsav(jj)
			
		end do
		if (abs(p(2)) > 20) then
			df(2)=abs(p(2))
			df(1)=1
		end if

		!WRITE(*,*)  df



		g=df
		!g=dfunc(p)
		den=max(fret,1.0_sp)
		if (maxval(abs(g)*max(abs(p),1.0_sp)/den) < gtol) RETURN
		dg=g-dg
		hdg=matmul(hessin,dg)
		fac=dot_product(dg,xi)
		fae=dot_product(dg,hdg)
		sumdg=dot_product(dg,dg)
		sumxi=dot_product(xi,xi)
		if (fac > sqrt(EPS*sumdg*sumxi)) then
			fac=1.0_sp/fac
			fad=1.0_sp/fae
			dg=fac*xi-fad*hdg
			hessin=hessin+fac*outerprod(xi,xi)-&
				fad*outerprod(hdg,hdg)+fae*outerprod(dg,dg)
		end if
		xi=-matmul(hessin,g)
	end do
	call nrerror('dfpmin: too many iterations')
	END SUBROUTINE dfpmin

!*********************************************************************************!
!Sobol quasi random number generator package
!
!*********************************************************************************!
subroutine inhalt ( flag, dimen, atmost, quasi )

!*****************************************************************************80
!
!! INHALT initializes the Halton quasirandom number generator.
!
!  Discussion:
!
!    INHALT first checks whether the user-supplied dimension DIMEN of
!    the quasirandom vectors is acceptable (between 2 and 40).
!
!    INHALT then calculates a tolerance parameter E to make the program work
!    correctly in finite precision arithmetic and a parameter DELTA
!    to check that E works.  If the test is not passed, then ATMOST
!    is too big relative to the machine precision.
!
!    Otherwise, INHALT computes and returns the first vector QUASI.
!    For the following values of QUASI, it is necessary to call GOHALT.
!
!  Modified:
!
!    18 March 2003
!
!  Reference:
!
!    Bennett Fox,
!    Algorithm 647:
!    Implementation and Relative Efficiency of Quasirandom 
!    Sequence Generators,
!    ACM Transactions on Mathematical Software,
!    Volume 12, Number 4, pages 362-376, 1986.
!
!    J H Halton and G B Smith,
!    Algorithm 247: Radical-Inverse Quasi-Random Point Sequence,
!    Communications of the ACM,
!    Volume 7, 1964, pages 701-702.
!
!  Parameters:
!
!    Output, logical FLAG(2), error flags.
!    FLAG(1) is FALSE if the input value of DIMEN is unacceptable.
!    FLAG(2) is FALSE if the input value of ATMOST is unacceptable.
!
!    Input, integer DIMEN, the spatial dimension.  DIMEN should
!    satisfy: 2 <= DIMEN <= 40.
!
!    Input, integer ATMOST, the maximum number of quasirandom
!    vectors to be computed.
!
!    Output, double precision QUASI(DIMEN), the first element of
!    the Halton sequence.
!
!  Local Parameters:
!
!    Local, double precision PRIME(40), the first 40 primes.
!
!  Global Parameters:
!
!    Stored in common block /HALTON/:
!
!    Global, double precision E, a tolerance.
!
!    Global, double precision PRIME_INV(40), the reciprocals of the
!    first 40 primes.
!
!    Global, integer S, the spatial dimension.
!
  implicit none

  integer dimen
  integer, parameter :: dim_max = 40

  integer atmost
  double precision delta
  double precision e
  logical flag(2)
  integer prime(dim_max)
  double precision prime_inv(dim_max)
  double precision quasi(dimen)
  integer s
  double precision small

  common /halton/ e, prime_inv, s

  save /halton/
!
!  Check DIMEN.
!
  flag(1) = .true.
  flag(2) = .true.

  s = dimen

  if ( s < 2 .or. dim_max < s ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'INHALT - Fatal error!'
    write ( *, '(a)' ) '  The spatial dimension S should satisfy:'
    write ( *, '(a,i6)' ) '    2 <= S <= ', dim_max
    write ( *, '(a,i6)' ) '  But this input value is S = ', s
    flag(1) = .false.
    return
  end if
!
!  Set the primes.
!
  prime(1:dim_max) = (/ &
      2,   3,   5,   7,  11,  13,  17,  19,  23,  29, &
     31,  37,  41,  43,  47,  53,  59,  61,  67,  71, &
     73,  79,  83,  89,  97, 101, 103, 107, 109, 113, &
    127, 131, 137, 139, 149, 151, 157, 163, 167, 173 /)
!
!  Compute the tolerance and make the check.
!
  small = epsilon ( small )

  e = 0.9D+00 * ( 1.0D+00 / ( dble ( atmost * prime(s) ) ) - 10.0D+00 * small )

  delta = 100.0D+00 * small * dble ( atmost + 1 ) * log10 ( dble ( atmost ) )

  if ( 0.09D+00 * ( e - 10.0D+00 * small ) < delta ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'INHALT - Fatal error!'
    write ( *, '(a)' ) '  The value of ATMOST is too great.'
    flag(2) = .false.
    return
  end if
!
!  Set the inverse primes.
!
  prime_inv(1:dim_max) = 1.0D+00 / dble ( prime(1:dim_max) )
!
!  Compute the first vector.
!
  quasi(1:s) = prime_inv(1:s)

  return
end subroutine
subroutine gohalt ( quasi )

!*****************************************************************************80
!
!! GOHALT generates a new quasirandom Halton vector with each call.
!
!  Discussion:
!
!    The routine adapts key ideas from Halton and Smith.
!
!    The routine INHALT must be called once before using
!    this routine.
!
!  Reference:
!
!    Bennett Fox,
!    Algorithm 647:
!    Implementation and Relative Efficiency of Quasirandom 
!    Sequence Generators,
!    ACM Transactions on Mathematical Software,
!    Volume 12, Number 4, pages 362-376, 1986.
!
!    J H Halton and G B Smith,
!    Algorithm 247: Radical-Inverse Quasi-Random Point Sequence,
!    Communications of the ACM,
!    Volume 7, 1964, pages 701-702.
!
!  Parameters:
!
!    Input/output, real QUASI(DIMEN), on input, the previous
!    quasirandom vector; on output, the next quasirandom vector.
!    On the first call, the input value should be the output
!    value given by INHALT.
!
!  Global Parameters:
!
!    In labelled common block /HALTON/:
!
!    Global, double precision E, a tolerance.
!
!    Global, double precision PRIME_INV(40), the reciprocal of
!    the first 40 prime numbers.
!
!    Global, integer S, the spatial dimension.
!
  implicit none

  integer, parameter :: dim_max = 40

  double precision e
  double precision f
  double precision g
  double precision h
  integer i
  double precision prime_inv(dim_max)
  double precision quasi(*)
  integer s
  double precision t

  common /halton/ e, prime_inv, s

  save /halton/
!
!  Generate QUASI one component at a time, using radix 1/PRIME(K) for 
!  component K.
!
  do i = 1, s

    t = prime_inv(i)
    f = 1.0D+00 - quasi(i)
    g = 1.0D+00
    h = prime_inv(i)

    do

      if ( e <= f - h ) then
        exit
      end if
!
!  This checks whether Q + H > 1 - E.
!
      g = h
      h = h * t
!
!  If this is the (K-1)-st time this statement is reached, check whether
!  QUASI(I) + R**(-K) > 1-E.
!
    end do
!
!  For the appropriate I (depending on how many times the loop above
!  is executed), add H**(I+1) + H**(I) - 1
!  to the old QUASI(I) to get the next QUASI(I).
!
    quasi(i) = g + h - f

  end do

  return
end subroutine

!*********************************************************************************!
!Amoeba
!
!*********************************************************************************!
		SUBROUTINE amoeba(p,y,ftol,func,iter)
	USE nrtype; USE nrutil, ONLY : assert_eq,imaxloc,iminloc,nrerror,swap
	IMPLICIT NONE
	INTEGER(I4B), INTENT(OUT) :: iter
	REAL(SP), INTENT(IN) :: ftol
	REAL(SP), DIMENSION(:), INTENT(INOUT) :: y
	REAL(SP), DIMENSION(:,:), INTENT(INOUT) :: p
	INTERFACE
		FUNCTION func(x)
		USE nrtype
		IMPLICIT NONE
		REAL(SP), DIMENSION(:), INTENT(IN) :: x
		REAL(SP) :: func
		END FUNCTION func
	END INTERFACE
	INTEGER(I4B), PARAMETER :: ITMAX=5000
	REAL(SP), PARAMETER :: TINY=1.0e-10
	INTEGER(I4B) :: ihi,ndim
	REAL(SP), DIMENSION(size(p,2)) :: psum
	call amoeba_private
	CONTAINS
!BL
	SUBROUTINE amoeba_private
	IMPLICIT NONE
	INTEGER(I4B) :: i,ilo,inhi
	REAL(SP) :: rtol,ysave,ytry,ytmp
	ndim=assert_eq(size(p,2),size(p,1)-1,size(y)-1,'amoeba')
	iter=0
	psum(:)=sum(p(:,:),dim=1)
	do
		ilo=iminloc(y(:))
		ihi=imaxloc(y(:))
		ytmp=y(ihi)
		y(ihi)=y(ilo)
		inhi=imaxloc(y(:))
		y(ihi)=ytmp
		rtol=2.0_sp*abs(y(ihi)-y(ilo))/(abs(y(ihi))+abs(y(ilo))+TINY)
		if (rtol < ftol) then
			call swap(y(1),y(ilo))
			call swap(p(1,:),p(ilo,:))
			RETURN
		end if
		if (iter >= ITMAX) call nrerror('ITMAX exceeded in amoeba')
		ytry=amotry(-1.0_sp)
		iter=iter+1
		if (ytry <= y(ilo)) then
			ytry=amotry(2.0_sp)
			iter=iter+1
		else if (ytry >= y(inhi)) then
			ysave=y(ihi)
			ytry=amotry(0.5_sp)
			iter=iter+1
			if (ytry >= ysave) then
				p(:,:)=0.5_sp*(p(:,:)+spread(p(ilo,:),1,size(p,1)))
				do i=1,ndim+1
					!if (i /= ilo) y(i)=func(p(i,:))
					if (i /= ilo) y(i)=((5*sin(p(i,1)))/p(i,1))*((max(20-abs(p(i,2)),0.0))**1.2)
				end do
				iter=iter+ndim
				psum(:)=sum(p(:,:),dim=1)
			end if
		end if
		!WRITE(*,*) y
	end do
	END SUBROUTINE amoeba_private
!BL

	FUNCTION amotry(fac)
	IMPLICIT NONE
	REAL(SP), INTENT(IN) :: fac
	REAL(SP) :: amotry
	REAL(SP) :: fac1,fac2,ytry
	REAL(SP), DIMENSION(size(p,2)) :: ptry
	fac1=(1.0_sp-fac)/ndim
	fac2=fac1-fac
	ptry(:)=psum(:)*fac1-p(ihi,:)*fac2
	ytry=((5*sin(ptry(1)))/ptry(1))*((max(20-abs(ptry(2)),0.0))**1.2)
	!ytry=func(ptry)
	if (ytry < y(ihi)) then
		y(ihi)=ytry
		psum(:)=psum(:)-p(ihi,:)+ptry(:)
		p(ihi,:)=ptry(:)
	end if
	amotry=ytry
	END FUNCTION amotry
	END SUBROUTINE amoeba


!******************************************************
!Computes forward-difference approximation to Jacobian. On input, x is the point at which
!the Jacobian is to be evaluated, and fvec is the vector of function values at the point,
!both arrays of length N. df is the N × N output Jacobian. FUNCTION funcv(x) is a
!fixed-name, user-supplied routine that returns the vector of functions at x.
!Parameter: EPS is the approximate square root of the machine precision.
!***********************************************************
!	SUBROUTINE fdjac(x,fvec,df)
!	USE nrtype; USE nrutil, ONLY : assert_eq
!	IMPLICIT NONE
!	REAL(SP), DIMENSION(:), INTENT(IN) :: fvec
!	REAL(SP), DIMENSION(:), INTENT(INOUT) :: x
!	REAL(SP), DIMENSION(:,:), INTENT(OUT) :: df
!	INTERFACE
!		FUNCTION funcv(x)
!		USE nrtype
!		IMPLICIT NONE
!		REAL(SP), DIMENSION(:), INTENT(IN) :: x
!		REAL(SP), DIMENSION(size(x)) :: funcv
!		END FUNCTION funcv
!	END INTERFACE

!	REAL(SP), PARAMETER :: EPS=1.0e-4_sp
!	INTEGER(I4B) :: j,n
!	REAL(SP), DIMENSION(size(x)) :: xsav,xph,h
!	n=assert_eq(size(x),size(fvec),size(df,1),size(df,2),'fdjac')
!	xsav=x
!	h=EPS*abs(xsav)
!	where (h == 0.0) h=EPS
!	xph=xsav+h
!	h=xph-xsav
	
!	do j=1,n
!		x(j)=xph(j)
!		df(:,j)=(funcv(x)-fvec(:))/h(j)
!		x(j)=xsav(j)
!	end do


!END SUBROUTINE fdjac





END MODULE minlib
