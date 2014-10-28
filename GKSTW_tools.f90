module GKSTW_tools

use nrtype
implicit none
!    INTEGER, PARAMETER :: SP = KIND(1.0D0)
!    INTEGER, PARAMETER :: DP = KIND(1.0D0)
!    INTEGER, PARAMETER :: WP = DP

contains
subroutine mat2csv(A,fname,append)
	real(8), dimension(:,:), intent(in) :: A
	character(LEN=*), intent(in) :: fname
	integer, intent(in), optional :: append
	CHARACTER(LEN=*), PARAMETER  :: FMT = "(G20.12)"
	CHARACTER(LEN=20) :: FMT_1
	integer :: r,c,ri,ci
	r = size(A,1)
	c = size(A,2)
	if(present(append)) then
		if(append .eq. 1) then 
			open(1, file=fname,ACCESS='APPEND', POSITION='APPEND')
		else
			open(1, file=fname)
		endif
	else
		open(1, file=fname)
	endif
	write(FMT_1, "(A1,I2,A7)") "(", c, "G20.12)"
	do ri=1,r
		!write(1,FMT_1) (A(ri,ci), ci = 1,c)
		do ci = 1,c-1
			write(1,FMT, advance='no') A(ri,ci)
		enddo
		write(1,FMT) A(ri,c)
	enddo
	write(1,*) " "! trailing space
	close(1)

end subroutine mat2csv

subroutine vec2csv(A,fname,append)
	real(8), dimension(:), intent(in) :: A
	character(len=*), intent(in) :: fname
	integer, intent(in), optional :: append
	integer :: r,ri
	r = size(A,1)
	if(present(append)) then
		if(append .eq. 1) then 
			open(1, file=fname,ACCESS='APPEND', POSITION='APPEND')
		else
			open(1, file=fname)
		endif
	else
		open(1, file=fname) 
	endif
	do ri=1,r
		write(1,*) A(ri)
	end do
	close(1)

end subroutine vec2csv


subroutine rouwenhorst(N,mu,rho,sigma,grid,PP)

!Purpose:    Finds a Markov chain whose sample paths approximate those of
!            the AR(1) process
!                z(t+1) = (1-rho)*mu + rho * z(t) + eps(t+1)
!            where eps are normal with stddev sigma
!
!Format:     [Z, PI] = rouwenhorst(N,mu,rho,sigma)
!
!Input:      N       scalar, number of nodes for Z
!            mu      scalar, unconditional mean of process
!            rho     scalar
!            sigma   scalar, std. dev. of epsilons
!
!Output:     Z       N*1 vector, nodes for Z
!            PI      N*N matrix, transition probabilities
!Adapted from Martin Floden's Matlab code by David Wiczer

	integer, intent(in)	:: N
	real(8), intent(in)	:: mu,rho,sigma
	real(8), dimension(N,N)	, intent(out)	:: PP
	real(8), dimension(N)	, intent(out)	:: grid
	real(8), dimension(N,N)	:: PPZ, PZP, ZPP	
	real(8)	:: sigmaz, p, fi
	integer :: i
	PP	= 0.0
	PPZ	= 0.0
	PZP	= 0.0
	ZPP	= 0.0		
	sigmaz	= sigma / sqrt(1.0-rho*rho)
	p	= (1.0+rho)/2.0
	PP(1,1:2)	= (/ p	 , 1.0-p/)
	PP(2,1:2)	= (/1.0-p,  p 	/)
	if (N>=3) then
	do i= 3,N
		PPZ	= 0.0
		PZP	= 0.0
		ZPP	= 0.0
		PPZ(1:i-1,2:i)	= PP(1:i-1,1:i-1)
		PZP(2:i,1:i-1)	= PP(1:i-1,1:i-1)
		PZP(2:i,2:i)	= PP(1:i-1,1:i-1)
		PP(1:i,1:i) 	= p*PP(1:i,1:i) + (1.0 - p)*PPZ(1:i,1:i) + (1.0 - p)*PZP(1:i,1:i) + p*ZPP(1:i,1:i)
		PP(2:i-1,1:i)	= PP(2:i-1,1:i)/2
	enddo
	endif
	fi	= sqrt(real(N)-1.0)*sigmaz
	grid	= (/ (i, i=0,N-1) /)
	grid	= grid *(2.0*fi / (real(N) - 1.0)) - fi + mu
end subroutine rouwenhorst

SUBROUTINE grid(x,xmin,xmax,s)
! Purpose: Generate grid x on [xmin,xmax] using spacing parameter s set as follows:
! s=1		linear spacing
! s>1		left skewed grid spacing with power s
! 0<s<1		right skewed grid spacing with power s
! s<0		geometric spacing with distances changing by a factor -s^(1/(n-1)), (>1 grow, <1 shrink)
! s=-1		logarithmic spacing with distances changing by a factor (xmax-xmin+1)^(1/(n-1))
! s=0		logarithmic spacing with distances changing by a factor (xmax/xmin)^(1/(n-1)), only if xmax,xmin>0
	implicit none
	real(8), dimension(:), intent(out) :: x
	real(8), intent(in) :: xmin,xmax,s
	real(8) :: c ! growth rate of grid subintervals for logarithmic spacing
	integer :: n,i
	n=size(x)
	forall(i=1:n) x(i)=real(i-1,8)/real(n-1,8)
	if (s>0.0) then
		x=x**s*(xmax-xmin)+xmin
		if (s==1.0) then
!			PRINT '(a,i8,a,f6.3,a,f6.3,a)', 'Using ',n,' equally spaced grid points over domain [',xmin,',',xmax,']'
		else
!			PRINT '(a,i8,a,f6.3,a,f6.3,a,f6.3,a)', 'Using ',n,' skewed spaced grid points with power ',s,' over domain [',xmin,',',xmax,']'
		end if
	else
		if (s==-1.0) then
			c=xmax-xmin+1
!		elseif (s==0.0_WP) then
!			if (xmin>0.0_WP) then
!				c=xmax/xmin
!			else
!				STOP 'grid: can not use logarithmic spacing for nonpositive values'
!			end if
		else
			c=-s
		end if
!		PRINT '(a,i8,a,f6.3,a,f6.3,a,f6.3,a)', 'Using ',n,' logarithmically spaced grid points with growth rate ',c,' over domain [',xmin,',',xmax,']'
		x=((xmax-xmin)/(c-1))*(c**x)-((xmax-c*xmin)/(c-1))
	end if
end SUBROUTINE grid

subroutine kron(TransA,TransB,A,B,TransC,C)
	real(8), dimension(:,:),intent(in)	:: A
	real(8), dimension(:,:),intent(in)	:: B
	!integer, intent(in)			:: rC,cC
	real(8), dimension(:,:),intent(out) 	:: C
	character*1, intent(in)	:: TransA,TransB,TransC
	integer			:: rA,cA,rB,cB, ri1,ri2,ci1,ci2
	if (TransA .eq. 'N') then
		rA = size(A,1)
		cA = size(A,2)
	else
		rA = size(A,2)
		cA = size(A,1)
	endif

	if (TransB .eq. 'N') then
		rB = size(B,1)
		cB = size(B,2)
	else
		rB = size(B,2)
		cB = size(B,1)
	endif
	!if (rA*rB .ne. rC) then
	!	stop 'Sizes are incorrect in kron product'
	!endif
	if (TransC .eq. 'N') then
		do ri1 = 1,rA
		do ci1 = 1,cA
			do ri2=1,rB
				forall(ci2 = 1:cB) C(rB*(ri1-1)+ri2,cB*(ci1-1)+ci2) = A(ri1,ci1)*B(ri2,ci2)
			enddo
		enddo
		enddo
	else
		do ri1 = 1,rA
		do ci1 = 1,cA
			do ri2=1,rB
				forall(ci2 = 1:cB) C(cB*(ci1-1)+ci2,rB*(ri1-1)+ri2) = A(ri1,ci1)*B(ri2,ci2)
			enddo
		enddo
		enddo	
	endif
end subroutine kron

function max_arr(arr1,arr2,max_loc) ! this should have a flag for single crossing
	real(8), dimension(:), intent(in) 	:: arr1,arr2
	integer, optional, dimension(:), intent(out)	:: max_loc
	real(8), dimension(size(arr1,1))	:: max_arr
	integer	:: i, ret_max=0
	if(present(max_loc)) then
		ret_max =1
	endif
	do i = 1,size(arr1,1)
		if (arr1(i)>=arr2(i)) then
			max_arr(i) = arr1(i)
			if(ret_max==1) then
				max_loc(i) = 1
			endif
		else
			max_arr(i) = arr2(i)
			if(ret_max==1) then
				max_loc(i) = 2
			endif
		endif
	enddo
	return
end function max_arr

function ones_mat(r,c)
	integer, intent(in) :: r,c
	integer	:: ri,ci
	real(8), dimension(r,c) :: ones_mat
	do ri =1,r
		forall (ci=1:c) ones_mat(ri,ci) = 1.000
	enddo
end function ones_mat

function ones_vec(r)
	integer, intent(in) :: r
	integer	:: ri
	real(8), dimension(r) :: ones_vec
	do ri =1,r
		ones_vec(ri) = 1.000
	enddo
end function ones_vec


function zeros_mat(r,c)
	integer, intent(in) :: r,c
	integer	:: ri,ci
	real(8), dimension(r,c) :: zeros_mat
	do ri =1,r
		forall (ci=1:c) zeros_mat(ri,ci) = 0.000
	enddo
end function zeros_mat

function zeros_vec(r)
	integer, intent(in) :: r
	integer	:: ri
	real(8), dimension(r) :: zeros_vec
	do ri =1,r
		zeros_vec(ri) = 0.000
	enddo
end function zeros_vec

subroutine invmat(A, invA)
	real(8), dimension(:,:), intent(in) :: A
	real(8), dimension(size(A,1),size(A,2)), intent(out) :: invA
	real(8), dimension(size(A,1)*size(A,1)) :: wk
	integer, dimension(size(A,1) + 1) :: ipiv
	integer :: n, info
	external DGETRF
	external DGETRI
	
	invA = A
	n = size(A,1)
	
	call DGETRF(n,n,invA,n,ipiv,info) ! first computes the LU factorization
	
	if (info /= 0) then
		print *, 'Matrix is singular.  Make it less singular before proceeding'
	endif
	call DGETRI(n,invA,n,ipiv,wk,n,info)
	if(info /=0) then
		print *, 'Matrix inversion failed, though it is not singular'
	endif
end subroutine invmat

subroutine addeqmat(A,B) ! A = A+B
	real(8), dimension(:,:), intent(inout)	:: A
	real(8), dimension(:,:), intent(in)	:: B	
	integer :: r,c,n,m
	n = size(A,1)
	m = size(A,2)
	do r=1,n
		forall (c=1:m) A(r,c) = A(r,c) + B(r,c)
	end do	
end subroutine addeqmat

subroutine subeqmat(A,B) ! A = A-B
	real(8), dimension(:,:), intent(inout) 	:: A
	real(8), dimension(:,:), intent(in) 	:: B
	integer :: r,c,n,m
	n = size(A,1)
	m = size(A,2)
	do r=1,n
		forall (c=1:m)	A(r,c) = A(r,c) - B(r,c)
	end do	
end subroutine subeqmat

function factorial(X)
	integer, intent(in) :: X
	integer :: factorial
	integer :: i
	factorial = 1
	if (X .gt. 1) then
	do i=2,X
		factorial = i*factorial
	enddo
	endif
	
end function factorial

! Local functional approximation by cubic and linear splines
subroutine spline(x,y,y2,yp1,ypn)
! this is the NR version of the spline

	IMPLICIT NONE
	REAL(WP), DIMENSION(:), INTENT(IN) :: x,y
	REAL(WP), DIMENSION(:), INTENT(OUT) :: y2
	REAL(WP), INTENT(IN), OPTIONAL :: yp1,ypn
	REAL(WP), DIMENSION(size(x)-1) :: u
	real(wp) :: p,qn,si,un
	INTEGER :: n,i,k
	n=size(x)
	IF (size(y)/=n .or. size(y2)/=n) THEN
		PRINT *, 'spline: x,y and y2 must be of the same size'
		STOP 'program terminated by spline'
	END IF

	IF (present(yp1)) THEN
		y2(1)=-0.5_wp
		u(1) = (1.0_wp/(x(2)-x(1) ))*( (y(2)-y(1))/(x(2)-x(1))-yp1)
	ELSE
		y2(1)=0.0_WP
		u(1)=0.0_WP
	END IF
	do i =2,n-1
		si	= (x(i)-x(i-1))/(x(i+1)-x(i-1))
		p	= si*y2(i-1)+2.0_wp
		y2(i)	= (si-1.0_wp)/p
		u(i)	= (6.0_wp*((y(i+1)-y(i))/(x(i+1)-x(i))-(y(i)-y(i-1)) &
			& /(x(i)-x(i-1)))/(x(i+1)-x(i-1))-si*u(i-1))/p
	enddo

	IF (present(ypn)) THEN
		qn = 0.5_wp
		un = (3.0_wp/(x(n)-x(n-1)))*( ypn-(y(n)-y(n-1))/(x(n)-x(n-1)) )
	ELSE
		qn = 0.0_wp
		un = 0.0_wp
!		y2(n)=y2(n-1)
!		a(n-1)=0.0_WP
	END IF
	y2(n) = (un-qn*u(n-1))/(qn*y2(n-1)+1.0_wp)
	do k=n-1,1,-1
		y2(k) = y2(k)*y2(k+1)+u(k)
	enddo



end subroutine spline


FUNCTION splint(x,y,y2,xi)
! cubic interpolation of function y on grid x at interpolation point xi
	IMPLICIT NONE
	REAL(WP), DIMENSION(:), INTENT(IN) :: x,y,y2
	REAL(WP), INTENT(IN) :: xi
	REAL(WP) :: splint
	REAL(WP) :: a,b,d,xhr
	INTEGER :: n,i
	n=size(x)
	IF (size(y)/=n .or. size(y2)/=n) THEN
		PRINT *, 'splint: x,y and y2 must be of the same size'
		STOP 'program terminated by splint'
	END IF
	i=max(min(locate(x,xi),n-1),1)
	d=x(i+1)-x(i)
	xhr = xi
	! push it back if too much extrapolation
	if(xi-x(n) .ge. d) then
		xhr = x(n)+d
	elseif(-xi+x(1) .ge. d) then
		xhr = x(1)-d
	endif
	
	if (d == 0.0) STOP 'bad x input in splint'
	a=(x(i+1)-xhr)/d
	b=(xhr-x(i))/d
	if((xhr .ge. x(1)) .and. (xhr .le. x(n))) then
		splint=a*y(i)+b*y(i+1)+((a**3-a)*y2(i)+(b**3-b)*y2(i+1))*(d**2)/6.0_WP
	elseif( xhr .ge. x(n) ) then
		splint = (y(n)-y(n-1))/(x(n)-x(n-1))*(xhr -x(n)) + y(n)
	elseif( xhr .le. x(1) ) then
		splint = (y(2)-y(1))/(x(2)-x(1))*(xhr -x(1)) + y(1)
	endif
	
END FUNCTION splint

FUNCTION dsplint(x,y,y2,xi)
! derivative implied by cubic interpolation of function y on grid x at interpolation point xi
	IMPLICIT NONE
	REAL(WP), DIMENSION(:), INTENT(IN) :: x,y,y2
	REAL(WP), INTENT(IN) :: xi
	REAL(WP) :: dsplint
	REAL(WP) :: a,b,d, xhr
	INTEGER :: n,i
	n=size(x)
	IF (size(y)/=n .or. size(y2)/=n) THEN
		PRINT *, 'splint: x,y and y2 must be of the same size'
		STOP 'program terminated by splint'
	END IF
	i=max(min(locate(x,xi),n-1),1)
	d=x(i+1)-x(i)
	xhr = xi
	if(xi-x(n) .ge. d) xhr = x(n)+d
	if (d == 0.0) STOP 'bad x input in dsplint'
	a=(x(i+1)-xhr)/d
	b=(xhr-x(i))/d
	dsplint=(y(i+1)-y(i))/d+((3*b**2-1)*y2(i+1)-(3*a**2-1)*y2(i))*d/6.0_WP
END FUNCTION dsplint

FUNCTION linint(x,y,xi)
! linear interpolation of function y on grid x at interpolation point xi
	IMPLICIT NONE
	REAL(WP), DIMENSION(:), INTENT(IN) :: x,y
	REAL(WP), INTENT(IN) :: xi
	REAL(WP) :: linint
	REAL(WP) :: a,b,d,xhr
	INTEGER :: n,i
	n=size(x)
	IF (size(y)/=n) THEN
		PRINT *, 'linint: x and y must be of the same size'
		STOP 'program terminated by linint'
	END IF
	i=max(min(locate(x,xi),n-1),1)
	d=x(i+1)-x(i)
	xhr = xi
	if(xi-x(n) .ge. d) xhr = x(n)+d	
	IF (d == 0.0) STOP 'bad x input in splint'
	a=(x(i+1)-xhr)/d
	b=(xhr-x(i))/d
	linint=a*y(i)+b*y(i+1)
END FUNCTION linint


FUNCTION dlinint(x,y,xi)
! derivative implied by linear interpolation of function y on grid x at interpolation point xi
	IMPLICIT NONE
	REAL(WP), DIMENSION(:), INTENT(IN) :: x,y
	REAL(WP), INTENT(IN) :: xi
	REAL(WP) :: dlinint
	REAL(WP) :: dx,dy
	INTEGER :: n,i
	n=size(x)
	IF (size(y)/=n) THEN
		PRINT *, 'linint: x and y must be of the same size'
		STOP 'program terminated by linint'
	END IF
	i=max(min(locate(x,xi),n-1),1)
	dx=x(i+1)-x(i)
	IF (dx == 0.0) STOP 'bad x input in splint'
	dy= y(i+1) - y(i)
	dlinint=dy/dx
END FUNCTION dlinint


SUBROUTINE linintv(x,y,xi,yi)
! linear interpolation of function y on grid x at interpolation vector xi
    IMPLICIT NONE
    REAL(WP), DIMENSION(:), INTENT(IN)  :: x,y,xi
    REAL(WP), DIMENSION(:), INTENT(OUT) :: yi
    REAL(WP) :: a,b,d
    INTEGER :: m,n,i,j
    n=size(x)
    IF (size(y)/=n) THEN
        PRINT *, 'linintv: x and y must be of the same size'
        STOP 'program terminated by linintv'
    END IF
    m=size(xi)
    IF (size(yi)/=m) THEN
        PRINT *, 'linintv: xi and yi must be of the same size'
        STOP 'program terminated by linintv'
    END IF
    DO j=1,m
        i=max(min(locate(x,xi(j)),n-1),1)
        d=x(i+1)-x(i)
        IF (d == 0.0) THEN
            STOP 'bad x input in linintv'
        END IF
        a=(x(i+1)-xi(j))/d
        b=(xi(j)-x(i))/d
        yi(j)=a*y(i)+b*y(i+1)
    END DO
END SUBROUTINE linintv

function bilinint(x,y,f,xiyi)
	implicit none
	real(wp), dimension(:), intent(in) :: x,y
	real(wp), dimension(:,:), intent(in) :: f
	real(wp), dimension(:), intent(in) :: xiyi
	real(wp) :: fq11,fq21,fq12,fq22, dx,dy
	
	real(wp) :: bilinint
	integer  :: x1,x2,y1,y2

	x1=0
	y1=0
	if(size(x)/=size(f,1)) stop 'x,f grids not the same length in bilinear interpolation'
	if(size(y)/=size(f,2)) stop 'y,f grids not the same length in bilinear interpolation'
	
	x1 = locate(x,xiyi(1))
	if(x1 .ge. size(x)) then
		x2 = size(x)
		x1 = x2-1
	elseif(x1 .le. 1) then
		x2 = 2
		x1 = 1
	else
		x2 = x1+1
	endif
	y1 = locate(y,xiyi(2))
	if(y1 .ge. size(y)) then
		y2 = size(y)
		y1 = y2-1
	elseif(y1 .le. 1) then
		y2=2
		y1=1
	else
		y2 = y1+1
	endif
	
	dx = x(x2) - x(x1)
	dy = y(y2) - y(y1)

	fq11 = f(x1,y1)
	fq21 = f(x2,y1)
	fq12 = f(x1,y2)
	fq22 = f(x2,y2)
	
	bilinint = (fq11*(x(x2)-xiyi(1))*(y(y2)-xiyi(2)) &
		&+ fq21*(xiyi(1)-x(x1))*(y(y2)-xiyi(2))  &
		&+ fq12*(x(x2)-xiyi(1))*(xiyi(2)-y(y1))  &
		&+ fq22*(xiyi(1)-x(x1))*(xiyi(2)-y(y1)))/(dx*dy)
end function

subroutine dbilinint(x,y,f,xiyi,dxdy)
	implicit none
	real(wp), dimension(:), intent(in) :: x,y
	real(wp), dimension(:,:), intent(in) :: f
	real(wp), dimension(:), intent(in) :: xiyi
	real(wp) :: fq11,fq21,fq12,fq22, dx,dy
	
	real(wp), dimension(2) :: dxdy
	integer  :: x1,x2,y1,y2

	if(size(x)/=size(f,1)) stop 'x,f grids not the same length in bilinear interpolation'
	if(size(y)/=size(f,2)) stop 'y,f grids not the same length in bilinear interpolation'
	
	x1 = locate(x,xiyi(1))
	if(x1 .ge. size(x)) then
		x2 = size(x)
		x1 = x2 - 1
	elseif(x1 .le. 1) then
		x1 =1 
		x2 =2
	else 
		x2 = x1+1
	endif
	y1 = locate(y,xiyi(2))
	if(y1 .ge. size(y)) then
		y2 = size(y)
		y1 = y2-1
	elseif(y1 .le. 1) then
		y1 = 1
		y2 = 2
	else
		y2 = y1+1
	endif
	
	dx = x(x2) - x(x1)
	dy = y(y2) - y(y1)

	fq11 = f(x1,y1)
	fq21 = f(x2,y1)
	fq12 = f(x1,y2)
	fq22 = f(x2,y2)
	
	dxdy(1) = ( -fq11*(y(y2)-xiyi(2)) &
		&+ fq21*(y(y2)-xiyi(2))  &
		&- fq12*(xiyi(2)-y(y1))  &
		&+ fq22*(xiyi(2)-y(y1)))/(dx*dy)
		
	dxdy(2) = (-fq11*(x(x2)-xiyi(1)) &
		&- fq21*(xiyi(1)-x(x1))  &
		&+ fq12*(x(x2)-xiyi(1))  &
		&+ fq22*(xiyi(1)-x(x1)))/(dx*dy)		
	
end subroutine

function bilinint_v(x,y,f,xiyi)
	implicit none
	real(wp), dimension(:), intent(in) :: x,y
	real(wp), dimension(:), intent(in) :: f
	real(wp), dimension(:), intent(in) :: xiyi
	real(wp) :: fq11,fq21,fq12,fq22, dx,dy
	
	real(wp) :: bilinint_v
	integer  :: x1,x2,y1,y2,Nx,Ny

	Nx = size(x)
	Ny = size(y)
	if(Nx*Ny/=size(f,1)) stop 'x,f grids not the same length in bilinear interpolation'
	
		
	
	x1 = locate(x,xiyi(1))
	if(x1 .ge. size(x)) then
		x2 = x1
		x1 = x1 - 1 
	else 
		x2 = x1+1
	endif
	y1 = locate(y,xiyi(2))
	if(y1 .ge. size(y)) then
		y2 = y1
		y1 = y1-1
	else
		y2 = y1+1
	endif
	
	dx = x(x2) - x(x1)
	dy = y(y2) - y(y1)

	fq11 = f((x1-1)*Ny+y1)
	fq21 = f((x2-1)*Ny+y1)
	fq12 = f((x1-1)*Ny+y2)
	fq22 = f((x2-1)*Ny+y2)
	
	bilinint_v = (fq11*(x(x2)-xiyi(1))*(y(y2)-xiyi(2)) &
		&+ fq21*(xiyi(1)-x(x1))*(y(y2)-xiyi(2))  &
		&+ fq12*(x(x2)-xiyi(1))*(xiyi(2)-y(y1))  &
		&+ fq22*(xiyi(1)-x(x1))*(xiyi(2)-y(y1)))/(dx*dy)
end function


function bisplint(x,y,f,coefs,xiyi)
!This is just a place holder, just calls bilinear interpolation right now
	implicit none
	real(wp), dimension(:), intent(in) :: x,y
	real(wp), dimension(:,:), intent(in) :: f,coefs
	real(wp), dimension(:), intent(in) :: xiyi

	real(wp) :: bisplint
	
	bisplint = bilinint(x,y,f,xiyi)

end function bisplint


function nlinint(x,f,xi, dims)
	implicit none
	real(wp), dimension(:), intent(in) :: x ! will have each dimension flattened
	real(wp), dimension(:), intent(in) :: f ! f defined over the flattened dimensions
	real(wp), dimension(:), intent(in) :: dims ! array containing the sizes of each dimension
	real(wp), dimension(:), intent(in) :: xi! where the interpolation should happen, lenght is going to be the number of dimensions
	real(wp) :: fq1s,fq2s

	real(wp) :: nlinint
	integer  :: Nh,si,sii
	integer, dimension(size(dims))  :: x1,x2, xL,xH
	real(wp), dimension(size(dims))  :: dx1

	Nh = size(dims)
	! put in some checks on the size of x f
	if(size(x)/=sum(dims) ) stop 'x needs to contain all of the dimensions on which f is defined'
	if(size(f)/=product(dims)) stop 'f is not defined on all of the dimensions of x'

	do si=1,Nh
		xL(si) = 1
		xH(si) = dims(1)
		if(si .gt. 1) then
			xL(si) = xH(si-1)+1
			xH(si) = xH(si-1)+dims(si)
		endif
		x1(si) = locate(x(xL(si):xH(si)),xi(si))
		if(x1(si) .ge. dims(si)) then
			x2(si) = x1(si)
			x1(si) = x1(si) - 1
		else
			x2(si) = x1(si)+1
		endif

		dx1(si) = (xi(si) - x(xL(si) + x1(si)-1 ))/( x(x2(si)+xL(si)-1) - x(x1(si)+xL(si)-1) )

	enddo

	nlinint = 0.0_dp
	do si=1,Nh
		fq1s = 0._dp
		fq2s = 0._dp		
		do sii =1,Nh 
			if(si .ne. sii)then
!				fq1s = dx1(sii)*f() +(1._dp-dx1(sii))*f()+ fq1s
!				fq1s = dx1(sii)*f() +(1._dp-dx1(sii))*f()+ fq1s			
			endif
		enddo	
	
	enddo


end function

PURE FUNCTION locate(xx,x)
! locates x on array xx and returns the integer index
	IMPLICIT NONE
	REAL(WP), DIMENSION(:), INTENT(IN) :: xx
	REAL(WP), INTENT(IN) :: x
	INTEGER :: locate
	INTEGER :: n,il,im,iu
	n=size(xx)
	il=0
	iu=n+1
	do
		if (iu-il <= 1) exit
			im=(iu+il)/2
		if (x >= xx(im)) then
			il=im
		else
			iu=im
		end if
	end do
	if (x <= xx(1)+epsilon(xx(1))) then
		locate=1
	else if (x >= xx(n)-epsilon(xx(n))) then
		locate=n-1
	else
		locate=il
	end if
END FUNCTION locate

FUNCTION locate_retvals(xx,x,xhi,xlo)
! locates x on array xx and returns the integer index
! output arguments xhi and xlo are the x values above and below x on xx
	IMPLICIT NONE
	REAL(WP), DIMENSION(:), INTENT(IN) :: xx
	REAL(WP), INTENT(OUT) :: xhi,xlo
	REAL(WP), INTENT(IN) :: x
	INTEGER :: locate_retvals
	INTEGER :: n,il,im,iu
	n=size(xx)
	il=0
	iu=n+1
	do
		if (iu-il <= 1) exit
			im=(iu+il)/2
		if (x >= xx(im)) then
			il=im
		else
			iu=im
		end if
	end do
	if (x <= xx(1)+epsilon(xx(1))) then
		locate_retvals = 1
		xhi = xx(2)
		xlo = xx(1)
	else if (x >= xx(n)-epsilon(xx(n))) then
		locate_retvals = n-1
		xlo = xx(n-1)
		xhi = xx(n)
	else
		locate_retvals = il
		xhi = xx(iu)
		xlo = xx(il)

	end if

END FUNCTION locate_retvals


PURE FUNCTION locate_idx(xx,idx,x)
! locates x on array xx(idx) where index integers idx map the naturals to the sorting
	IMPLICIT NONE
	REAL(WP), DIMENSION(:), INTENT(IN) :: xx
	INTEGER , DIMENSION(:), INTENT(IN) :: idx
	REAL(WP), INTENT(IN) :: x
	INTEGER :: locate_idx
	INTEGER :: n,il,im,iu
	n=size(xx)
	if(n .ne. size(idx)) then
		locate_idx = -1
		return
	endif
	il=0
	iu=n+1
	do
		if (iu-il <= 1) exit
			im=(iu+il)/2
		if (x >= xx(idx(im))) then
			il=im
		else
			iu=im
		end if
	end do
	if (x <= xx(idx(1))+epsilon(xx(idx(1)))) then
		locate_idx=1
	else if (x >= xx( idx(n) )-epsilon(xx( idx(n) ))) then
		locate_idx=n-1
	else
		locate_idx=il
	end if
END FUNCTION locate_idx


FUNCTION brent(ax,bx,cx,func,xmin, funcp,info,tol_in,niter)
! inputs: 
! 	ax,bx,cx define the domain for the optimum
!	funcp is a vector of parameters
! 	func is a function that takes x and funcp
! outputs:
!	brent is the function value at the optimum
!	xmin is the arg min 
! 	info is the status, 0 for sucess and 1 for max iterations

	IMPLICIT NONE
	REAL(WP), INTENT(IN) :: ax,bx,cx
	REAL(WP), INTENT(IN), optional :: tol_in
	REAL(WP), INTENT(OUT) :: xmin
	integer , intent(out) :: info
	integer , intent(out), optional :: niter
	REAL(WP) :: brent
	real(wp), dimension(:), intent(in) :: funcp ! a vector of function parameters
	INTERFACE
		FUNCTION func(x, funcp)
		use nrtype
!		USE mkl95_precision, ONLY: WP => DP
		IMPLICIT NONE
		REAL(WP), INTENT(IN) :: x
		REAL(WP), INTENT(IN), dimension(:) :: funcp
		REAL(WP) :: func
		END FUNCTION func
	END INTERFACE
	INTEGER, PARAMETER :: ITMAX=100
	real(wp) :: TOL
	REAL(WP), PARAMETER :: CGOLD=0.381966011250105_WP,ZEPS=1.0e-3_WP*epsilon(ax)
	INTEGER :: iter
	REAL(WP) :: a,b,d,e,etemp,fu,fv,fw,fx,p,q,r,tol1,tol2,u,v,w,x,xm
	info = 0

	if(present(tol_in) .eqv. .true.) then 
		tol = tol_in
	else
		tol =sqrt(epsilon(ax))
	endif

	a=min(ax,cx)
	b=max(ax,cx)
	v=bx
	w=v
	x=v
	e=0.0
	!if(present(funcp)) then
	fx=func(x, funcp)
	!else
	!	fx=func(x)
	!endif
	fv=fx
	fw=fx
	do iter=1,ITMAX
		xm=0.5_WP*(a+b)
		tol1=tol*abs(x)+ZEPS
		tol2=2.0_WP*tol1
		if (abs(x-xm) <= (tol2-0.5_WP*(b-a))) then
			xmin=x
			brent=fx
			if( (present(niter).eqv. .true.) ) niter = iter-1
			RETURN
		end if
		if (abs(e) > tol1) then
			r=(x-w)*(fx-fv)
			q=(x-v)*(fx-fw)
			p=(x-v)*q-(x-w)*r
			q=2.0_WP*(q-r)
			if (q > 0.0) p=-p
			q=abs(q)
			etemp=e
			e=d
			if (abs(p) >= abs(0.5_WP*q*etemp) .or. &
				p <= q*(a-x) .or. p >= q*(b-x)) then
				e=merge(a-x,b-x, x >= xm )
				d=CGOLD*e
			else
				d=p/q
				u=x+d
				if (u-a < tol2 .or. b-u < tol2) d=sign(tol1,xm-x)
			end if
		else
			e=merge(a-x,b-x, x >= xm )
			d=CGOLD*e
		end if
		u=merge(x+d,x+sign(tol1,d), abs(d) >= tol1 )
		!if(present(funcp)) then
		fu=func(u, funcp)
		!else
		!	fu=func(u)
		!endif
		if (fu <= fx) then
			if (u >= x) then
				a=x
			else
				b=x
			end if
			v=w
			fv=fw
			w=x
			fw=fx
			x=u
			fx=fu
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
	end do
	info = 1
	brent = fx
END FUNCTION brent


FUNCTION zbrent(func,x1,x2,funcp,tol,flag)
	IMPLICIT NONE
	REAL(WP), INTENT(IN) :: x1,x2,tol
	REAL(WP) :: zbrent
	real(wp), dimension(:), intent(in) :: funcp ! a vector of function parameters
	integer , intent(out) :: flag
	INTERFACE
		FUNCTION func(x, funcp)
		use nrtype
!		USE mkl95_precision, ONLY: WP => DP
		IMPLICIT NONE
		REAL(WP), INTENT(IN) :: x
		REAL(WP), INTENT(IN), dimension(:) :: funcp
		REAL(WP) :: func
		END FUNCTION func
	END INTERFACE
	INTEGER, PARAMETER :: ITMAX=100
	REAL(WP), PARAMETER :: EPS=epsilon(x1)
	INTEGER :: iter
	REAL(WP) :: a,b,c,d,e,fa,fb,fc,p,q,r,s,tol1,xm
	a=x1
	b=x2
	fa=func(a, funcp)
	fb=func(b, funcp)
	if ((fa > 0.0 .and. fb > 0.0) .or. (fa < 0.0 .and. fb < 0.0)) then
		if(abs(fa) < abs(fb)) then
			zbrent = a
		else
			zbrent = b
		endif
		flag = -1
		return
		!STOP 'root must be bracketed for zbrent'
	endif
	c=b
	fc=fb
	do iter=1,ITMAX
		if ((fb > 0.0 .and. fc > 0.0) .or. (fb < 0.0 .and. fc < 0.0)) then
			c=a
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
		tol1=2.0_WP*EPS*abs(b)+0.5_WP*tol
		xm=0.5_WP*(c-b)
		if (abs(xm) <= tol1 .or. fb == 0.0) then
			zbrent=b
			flag = 0
			RETURN
		end if
		if (abs(e) >= tol1 .and. abs(fa) > abs(fb)) then
			s=fb/fa
			if (a == c) then
				p=2.0_WP*xm*s
				q=1.0_WP-s
			else
				q=fa/fc
				r=fb/fc
				p=s*(2.0_WP*xm*q*(q-r)-(b-a)*(r-1.0_WP))
				q=(q-1.0_WP)*(r-1.0_WP)*(s-1.0_WP)
			end if
			if (p > 0.0) q=-q
			p=abs(p)
			if (2.0_WP*p  <  min(3.0_WP*xm*q-abs(tol1*q),abs(e*q))) then
				e=d
				d=p/q
			else
				d=xm
				e=d
			end if
		else
			d=xm
			e=d
		end if
		a=b
		fa=fb
		b=b+merge(d,sign(tol1,xm), abs(d) > tol1 )
		fb=func(b, funcp)
	end do
	!STOP 'zbrent: exceeded maximum iterations'
	zbrent = b
	flag = -2
	zbrent=b
END FUNCTION zbrent



function trapezoid_refine ( a, b, m, f, q )

!    This routine is designed to be an efficient way to carry out a 
!    sequence of integral estimates, using the trapezoidal rule
!    and a nested sequence of evaluation points, in such a way that
!    the function is not re-evaluated unnecessarily.
!
!    The user calls first with M = 1 and Q = 0 to get a 2 point
!    integral estimate.  On the second call, the user sets M = 2,
!    and the input value of Q should be the integral estimate returned
!    on the previous call.  By incrementing M and passing the previous
!    estimate, the user gets a sequence of increasingly improved
!    integral estimates:
!
!    q = 0.0
!    m = 1
!    do
!      q_new = trapezoid_refine ( a, b, m, f, q )
!      if ( satisfied ) then
!        exit
!      end if
!      q = q_new
!      m = m + 1
!    end do
!
!    The number of points used on each step of the iteration is:
!
!    M   N 
!    1   2
!    2   3
!    3   5
!    4   9
!    5  17
!    6  33
!    m   2^(m-1)+1
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 February 2013
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A, B, the endpoints of the interval.
!
!    Input, integer ( kind = 4 ) M, the current level of refinement.
!    On the first call, M should be 1.  On each subsequent call,
!    the user should increment M by 1.
!
!    Input, external real ( kind = 8 ) F, the name of the function
!    which evaluates the integrand, and whose form is:
!      function f ( x )
!      real ( kind = 8 ) f
!      real ( kind = 8 ) x
!      f = ...
!
!    Input, real ( kind = 8 ) Q, the integral estimate return on
!    the previous call.  But on the first call, set Q to zero.
!
!    Output, real ( kind = 8 ) TRAPEZOID_REFINE, the improved 
!    integral estimate.
!
  implicit none

  real (wp)  :: a, b
  real ( kind = 8 ), external :: f
  integer :: i,k,m

  real (wp) :: q, trapezoid_refine,value,x

  if ( m < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TRAPEZOID_REFINE - Fatal error!'
    write ( *, '(a)' ) '  Illegal input value of M.'
    stop
  else if ( m == 1 ) then
    value = ( b - a ) / 2.0D+00 * ( f ( a ) + f ( b ) )
  else if ( 2 <= m ) then
    k = 2 ** ( m - 2 )
    value = 0.0D+00
    do i = 1, k
      x = ( real ( 2 * k - 2 * i + 1, kind = 8 ) * a   &
          + real (         2 * i - 1, kind = 8 ) * b ) &
          / real ( 2 * k,             kind = 8 )
      value = value + f ( x )
    end do

    value = 0.5D+00 * q + ( b - a ) * value / real ( 2 * k, kind = 8 )

  end if

  trapezoid_refine = value

  return
end function trapezoid_refine

function romberg_trap ( a, b, f, tol, m )

!*****************************************************************************80
!
!! ROMBERG_TRAP approximates an integral by extrapolating the trapezoidal rule.
!
!  Discussion:
!
!    This routine computes a sequence of integral estimates involving the 
!    trapezoid rule.  Extrapolation of successive trapezoidal estimates 
!    produces an estimate with a higher rate of convergence.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 February 2013
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A, B, the endpoints of the interval.
!
!    Input, external real ( kind = 8 ) F, the name of the function
!    which evaluates the integrand, and whose form is:
!      function f ( x )
!      real ( kind = 8 ) f
!      real ( kind = 8 ) x
!      f = ...
!
!    Input, real ( kind = 8 ) TOL, the tolerance.  When two successive integral
!    estimates differ by less than this value, the iteration will halt.
!
!    Output, integer ( kind = 4 ) M, the number of trapezoidal estimates
!    that were required.
!
!    Output, real ( kind = 8 ) ROMBERG_TRAP, the integral estimate.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ), external :: f
  integer ( kind = 4 ) m
  real ( kind = 8 ) q
  real ( kind = 8 ) q_new
  real ( kind = 8 ) r
  real ( kind = 8 ) r_new
  real ( kind = 8 ) romberg_trap
  real ( kind = 8 ) tol


  q_new = 0.0D+00
  r_new = 0.0D+00
  m = 1

  do

    q = q_new
    q_new = trapezoid_refine ( a, b, m, f, q )

    if ( m == 1 ) then
      r_new = q_new
    else
      r = r_new
      r_new = ( 4.0D+00 * q_new - q ) / 3.0D+00
      if ( abs ( r_new - r ) .lt. tol * ( 1.0 + abs ( r_new ) ) ) then
        exit
      end if
    end if

    if ( 20 <= m ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'ROMBERG_TRAP - Fatal error!'
      write ( *, '(a)' ) '  No convergence in 20 iterations.'
      write ( *, '(a)' ) '  The algorithm is halting.'
      stop
    end if

    m = m + 1

  end do

  romberg_trap = r_new

  return
end function romberg_trap


SUBROUTINE qsortd(x,ind,n)
	 
	! Code converted using TO_F90 by Alan Miller
	! Date: 2002-12-18  Time: 11:55:47
	! input is an array x, and output is an index array ind.  N should be the size of the two arrays

	IMPLICIT NONE
	INTEGER, PARAMETER  :: dp = SELECTED_REAL_KIND(12, 60)

	REAL (dp), INTENT(IN)  :: x(:)
	INTEGER, INTENT(OUT)   :: ind(:)
	INTEGER, INTENT(IN)    :: n


	! NOTE -- IU AND IL MUST BE DIMENSIONED >= LOG(N) WHERE LOG HAS BASE 2.

	INTEGER   :: iu(21), il(21)
	INTEGER   :: m, i, j, k, l, ij, it, itt, indx
	REAL      :: r
	REAL (dp) :: t

	! LOCAL PARAMETERS -

	! IU,IL =  TEMPORARY STORAGE FOR THE UPPER AND LOWER
	!            INDICES OF PORTIONS OF THE ARRAY X
	! M =      INDEX FOR IU AND IL
	! I,J =    LOWER AND UPPER INDICES OF A PORTION OF X
	! K,L =    INDICES IN THE RANGE I,...,J
	! IJ =     RANDOMLY CHOSEN INDEX BETWEEN I AND J
	! IT,ITT = TEMPORARY STORAGE FOR INTERCHANGES IN IND
	! INDX =   TEMPORARY INDEX FOR X
	! R =      PSEUDO RANDOM NUMBER FOR GENERATING IJ
	! T =      CENTRAL ELEMENT OF X

	IF (n <= 0) RETURN

	! INITIALIZE IND, M, I, J, AND R

	DO  i = 1, n
	  ind(i) = i
	END DO
	m = 1
	i = 1
	j = n
	r = .375

	! TOP OF LOOP

	20 IF (i >= j) GO TO 70
	IF (r <= .5898437) THEN
	  r = r + .0390625
	ELSE
	  r = r - .21875
	END IF

	! INITIALIZE K

	30 k = i

	! SELECT A CENTRAL ELEMENT OF X AND SAVE IT IN T

	ij = i + r*(j-i)
	it = ind(ij)
	t = x(it)

	! IF THE FIRST ELEMENT OF THE ARRAY IS GREATER THAN T,
	!   INTERCHANGE IT WITH T

	indx = ind(i)
	IF (x(indx) > t) THEN
	  ind(ij) = indx
	  ind(i) = it
	  it = indx
	  t = x(it)
	END IF

	! INITIALIZE L

	l = j

	! IF THE LAST ELEMENT OF THE ARRAY IS LESS THAN T,
	!   INTERCHANGE IT WITH T

	indx = ind(j)
	IF (x(indx) >= t) GO TO 50
	ind(ij) = indx
	ind(j) = it
	it = indx
	t = x(it)

	! IF THE FIRST ELEMENT OF THE ARRAY IS GREATER THAN T,
	!   INTERCHANGE IT WITH T

	indx = ind(i)
	IF (x(indx) <= t) GO TO 50
	ind(ij) = indx
	ind(i) = it
	it = indx
	t = x(it)
	GO TO 50

	! INTERCHANGE ELEMENTS K AND L

	40 itt = ind(l)
	ind(l) = ind(k)
	ind(k) = itt

	! FIND AN ELEMENT IN THE UPPER PART OF THE ARRAY WHICH IS
	!   NOT LARGER THAN T

	50 l = l - 1
	indx = ind(l)
	IF (x(indx) > t) GO TO 50

	! FIND AN ELEMENT IN THE LOWER PART OF THE ARRAY WHCIH IS NOT SMALLER THAN T

	60 k = k + 1
	indx = ind(k)
	IF (x(indx) < t) GO TO 60

	! IF K <= L, INTERCHANGE ELEMENTS K AND L

	IF (k <= l) GO TO 40
	! SAVE THE UPPER AND LOWER SUBSCRIPTS OF THE PORTION OF THE
	!   ARRAY YET TO BE SORTED

	IF (l-i > j-k) THEN
	  il(m) = i
	  iu(m) = l
	  i = k
	  m = m + 1
	  GO TO 80
	END IF

	il(m) = k
	iu(m) = j
	j = l
	m = m + 1
	GO TO 80

	! BEGIN AGAIN ON ANOTHER UNSORTED PORTION OF THE ARRAY

	70 m = m - 1
	IF (m == 0) RETURN
	i = il(m)
	j = iu(m)

	80 IF (j-i >= 11) GO TO 30
	IF (i == 1) GO TO 20
	i = i - 1

	! SORT ELEMENTS I+1,...,J.  NOTE THAT 1 <= I < J AND J-I < 11.

	90 i = i + 1
	IF (i == j) GO TO 70
	indx = ind(i+1)
	t = x(indx)
	it = indx
	indx = ind(i)
	IF (x(indx) <= t) GO TO 90
	k = i

	100 ind(k+1) = ind(k)
	k = k - 1
	indx = ind(k)
	IF (t < x(indx)) GO TO 100

	ind(k+1) = it
	GO TO 90
END SUBROUTINE qsortd


end module utils


MODULE Solve_NonLin
	use utils, only:mat2csv

! Corrections to FUNCTION Enorm - 28 November 2003

IMPLICIT NONE
INTEGER, PARAMETER  :: dp = SELECTED_REAL_KIND(14, 60)
PRIVATE
PUBLIC  :: hbrd, hybrd


CONTAINS


SUBROUTINE hbrd(fcn, n, x, fvec, epsfcn, tol, info, diag)
 
! Code converted using TO_F90 by Alan Miller
! Date: 2003-07-15  Time: 13:27:42

INTEGER, INTENT(IN)        :: n
REAL (dp), INTENT(IN OUT)  :: x(n)
REAL (dp), INTENT(IN OUT)  :: fvec(n)
REAL (dp), INTENT(IN)      :: epsfcn
REAL (dp), INTENT(IN)      :: tol
INTEGER, INTENT(OUT)       :: info
REAL (dp), INTENT(OUT)     :: diag(n)

! EXTERNAL fcn
INTERFACE
  SUBROUTINE FCN(N, X, FVEC, IFLAG)
    IMPLICIT NONE
    INTEGER, PARAMETER  :: dp = SELECTED_REAL_KIND(14, 60)
    INTEGER, INTENT(IN)      :: n
    REAL (dp), INTENT(IN)    :: x(n)
    REAL (dp), INTENT(OUT)   :: fvec(n)
    INTEGER, INTENT(IN OUT)  :: iflag
  END SUBROUTINE FCN
END INTERFACE

!   **********

!   SUBROUTINE HBRD

!     THE PURPOSE OF HBRD IS TO FIND A ZERO OF A SYSTEM OF N NONLINEAR
!   FUNCTIONS IN N VARIABLES BY A MODIFICATION OF THE POWELL HYBRID METHOD.
!   THIS IS DONE BY USING THE MORE GENERAL NONLINEAR EQUATION SOLVER HYBRD.
!   THE USER MUST PROVIDE A SUBROUTINE WHICH CALCULATES THE FUNCTIONS.
!   THE JACOBIAN IS THEN CALCULATED BY A FORWARD-DIFFERENCE APPROXIMATION.

!   THE SUBROUTINE STATEMENT IS

!     SUBROUTINE HBRD(N, X, FVEC, EPSFCN, TOL, INFO, WA, LWA)

!   WHERE

!     FCN IS THE NAME OF THE USER-SUPPLIED SUBROUTINE WHICH CALCULATES
!       THE FUNCTIONS.  FCN MUST BE DECLARED IN AN EXTERNAL STATEMENT
!       IN THE USER CALLING PROGRAM, AND SHOULD BE WRITTEN AS FOLLOWS.

!       SUBROUTINE FCN(N, X, FVEC, IFLAG)
!       INTEGER N,IFLAG
!       REAL X(N),FVEC(N)
!       ----------
!       CALCULATE THE FUNCTIONS AT X AND RETURN THIS VECTOR IN FVEC.
!       ---------
!       RETURN
!       END

!       THE VALUE OF IFLAG NOT BE CHANGED BY FCN UNLESS
!       THE USER WANTS TO TERMINATE THE EXECUTION OF HBRD.
!       IN THIS CASE SET IFLAG TO A NEGATIVE INTEGER.

!     N IS A POSITIVE INTEGER INPUT VARIABLE SET TO THE NUMBER
!       OF FUNCTIONS AND VARIABLES.

!     X IS AN ARRAY OF LENGTH N. ON INPUT X MUST CONTAIN AN INITIAL
!       ESTIMATE OF THE SOLUTION VECTOR.  ON OUTPUT X CONTAINS THE
!       FINAL ESTIMATE OF THE SOLUTION VECTOR.

!     FVEC IS AN OUTPUT ARRAY OF LENGTH N WHICH CONTAINS
!       THE FUNCTIONS EVALUATED AT THE OUTPUT X.

!     EPSFCN IS AN INPUT VARIABLE USED IN DETERMINING A SUITABLE STEP LENGTH
!       FOR THE FORWARD-DIFFERENCE APPROXIMATION.  THIS APPROXIMATION ASSUMES
!       THAT THE RELATIVE ERRORS IN THE FUNCTIONS ARE OF THE ORDER OF EPSFCN.
!       IF EPSFCN IS LESS THAN THE MACHINE PRECISION, IT IS ASSUMED THAT THE
!       RELATIVE ERRORS IN THE FUNCTIONS ARE OF THE ORDER OF THE MACHINE
!       PRECISION.

!     TOL IS A NONNEGATIVE INPUT VARIABLE.  TERMINATION OCCURS WHEN THE
!       ALGORITHM ESTIMATES THAT THE RELATIVE ERROR BETWEEN X AND THE SOLUTION
!       IS AT MOST TOL.

!     INFO IS AN INTEGER OUTPUT VARIABLE.  IF THE USER HAS TERMINATED
!       EXECUTION, INFO IS SET TO THE (NEGATIVE) VALUE OF IFLAG.
!       SEE DESCRIPTION OF FCN.  OTHERWISE, INFO IS SET AS FOLLOWS.

!       INFO = 0   IMPROPER INPUT PARAMETERS.

!       INFO = 1   ALGORITHM ESTIMATES THAT THE RELATIVE ERROR
!                  BETWEEN X AND THE SOLUTION IS AT MOST TOL.

!       INFO = 2   NUMBER OF CALLS TO FCN HAS REACHED OR EXCEEDED 200*(N+1).

!       INFO = 3   TOL IS TOO SMALL. NO FURTHER IMPROVEMENT IN
!                  THE APPROXIMATE SOLUTION X IS POSSIBLE.

!       INFO = 4   ITERATION IS NOT MAKING GOOD PROGRESS.

!   SUBPROGRAMS CALLED

!     USER-SUPPLIED ...... FCN

!     MINPACK-SUPPLIED ... HYBRD

!   ARGONNE NATIONAL LABORATORY. MINPACK PROJECT. MARCH 1980.
!   BURTON S. GARBOW, KENNETH E. HILLSTROM, JORGE J. MORE

! Reference:
! Powell, M.J.D. 'A hybrid method for nonlinear equations' in Numerical Methods
!      for Nonlinear Algebraic Equations', P.Rabinowitz (editor), Gordon and
!      Breach, London 1970.
!   **********
INTEGER    :: maxfev, ml, mode, mu, nfev, nprint
REAL (dp)  :: xtol
REAL (dp), PARAMETER  :: factor = 100.0_dp, zero = 0.0_dp

info = 0

!     CHECK THE INPUT PARAMETERS FOR ERRORS.


IF (n <= 0 .OR. epsfcn < zero .OR. tol < zero) GO TO 20

!     CALL HYBRD.

maxfev = 400*(n + 1)
xtol = tol
ml = n - 1
mu = n - 1
mode = 2
nprint = 0
CALL hybrd(fcn, n, x, fvec, xtol, maxfev, ml, mu, epsfcn, diag, mode,  &
           factor, nprint, info, nfev)
IF (info == 5) info = 4
20 RETURN

!     LAST CARD OF SUBROUTINE HBRD.

END SUBROUTINE hbrd



SUBROUTINE hybrd(fcn, n, x, fvec, xtol, maxfev, ml, mu, epsfcn, diag, mode,  &
                 factor, nprint, info, nfev)

INTEGER, INTENT(IN)        :: n
REAL (dp), INTENT(IN OUT)  :: x(n)
REAL (dp), INTENT(IN OUT)  :: fvec(n)
REAL (dp), INTENT(IN)      :: xtol
INTEGER, INTENT(IN OUT)    :: maxfev
INTEGER, INTENT(IN OUT)    :: ml
INTEGER, INTENT(IN)        :: mu
REAL (dp), INTENT(IN)      :: epsfcn
REAL (dp), INTENT(OUT)     :: diag(n)
INTEGER, INTENT(IN)        :: mode
REAL (dp), INTENT(IN)      :: factor
INTEGER, INTENT(IN OUT)    :: nprint
INTEGER, INTENT(OUT)       :: info
INTEGER, INTENT(OUT)       :: nfev

! EXTERNAL fcn
INTERFACE
  SUBROUTINE FCN(N, X, FVEC, IFLAG)
    IMPLICIT NONE
    INTEGER, PARAMETER  :: dp = SELECTED_REAL_KIND(14, 60)
    INTEGER, INTENT(IN)      :: n
    REAL (dp), INTENT(IN)    :: x(n)
    REAL (dp), INTENT(OUT)   :: fvec(n)
    INTEGER, INTENT(IN OUT)  :: iflag
  END SUBROUTINE FCN
END INTERFACE

!   **********

!   SUBROUTINE HYBRD

!   THE PURPOSE OF HYBRD IS TO FIND A ZERO OF A SYSTEM OF N NONLINEAR
!   FUNCTIONS IN N VARIABLES BY A MODIFICATION OF THE POWELL HYBRID METHOD.
!   THE USER MUST PROVIDE A SUBROUTINE WHICH CALCULATES THE FUNCTIONS.
!   THE JACOBIAN IS THEN CALCULATED BY A FORWARD-DIFFERENCE APPROXIMATION.

!   THE SUBROUTINE STATEMENT IS

!     SUBROUTINE HYBRD(FCN, N, X, FVEC, XTOL, MAXFEV, ML, MU, EPSFCN,
!                      DIAG, MODE, FACTOR, NPRINT, INFO, NFEV, FJAC,
!                      LDFJAC, R, LR, QTF, WA1, WA2, WA3, WA4)

!   WHERE

!     FCN IS THE NAME OF THE USER-SUPPLIED SUBROUTINE WHICH CALCULATES
!       THE FUNCTIONS.  FCN MUST BE DECLARED IN AN EXTERNAL STATEMENT IN
!       THE USER CALLING PROGRAM, AND SHOULD BE WRITTEN AS FOLLOWS.

!       SUBROUTINE FCN(N, X, FVEC, IFLAG)
!       INTEGER N, IFLAG
!       REAL X(N), FVEC(N)
!       ----------
!       CALCULATE THE FUNCTIONS AT X AND
!       RETURN THIS VECTOR IN FVEC.
!       ---------
!       RETURN
!       END

!       THE VALUE OF IFLAG SHOULD NOT BE CHANGED BY FCN UNLESS
!       THE USER WANTS TO TERMINATE EXECUTION OF HYBRD.
!       IN THIS CASE SET IFLAG TO A NEGATIVE INTEGER.

!     N IS A POSITIVE INTEGER INPUT VARIABLE SET TO THE NUMBER
!       OF FUNCTIONS AND VARIABLES.

!     X IS AN ARRAY OF LENGTH N.  ON INPUT X MUST CONTAIN AN INITIAL
!       ESTIMATE OF THE SOLUTION VECTOR.  ON OUTPUT X CONTAINS THE FINAL
!       ESTIMATE OF THE SOLUTION VECTOR.

!     FVEC IS AN OUTPUT ARRAY OF LENGTH N WHICH CONTAINS
!       THE FUNCTIONS EVALUATED AT THE OUTPUT X.

!     XTOL IS A NONNEGATIVE INPUT VARIABLE.  TERMINATION OCCURS WHEN THE
!       RELATIVE ERROR BETWEEN TWO CONSECUTIVE ITERATES IS AT MOST XTOL.

!     MAXFEV IS A POSITIVE INTEGER INPUT VARIABLE.  TERMINATION OCCURS WHEN
!       THE NUMBER OF CALLS TO FCN IS AT LEAST MAXFEV BY THE END OF AN
!       ITERATION.

!     ML IS A NONNEGATIVE INTEGER INPUT VARIABLE WHICH SPECIFIES THE
!       NUMBER OF SUBDIAGONALS WITHIN THE BAND OF THE JACOBIAN MATRIX.
!       IF THE JACOBIAN IS NOT BANDED, SET ML TO AT LEAST N - 1.

!     MU IS A NONNEGATIVE INTEGER INPUT VARIABLE WHICH SPECIFIES THE NUMBER
!       OF SUPERDIAGONALS WITHIN THE BAND OF THE JACOBIAN MATRIX.
!       IF THE JACOBIAN IS NOT BANDED, SET MU TO AT LEAST N - 1.

!     EPSFCN IS AN INPUT VARIABLE USED IN DETERMINING A SUITABLE STEP LENGTH
!       FOR THE FORWARD-DIFFERENCE APPROXIMATION.  THIS APPROXIMATION
!       ASSUMES THAT THE RELATIVE ERRORS IN THE FUNCTIONS ARE OF THE ORDER
!       OF EPSFCN. IF EPSFCN IS LESS THAN THE MACHINE PRECISION,
!       IT IS ASSUMED THAT THE RELATIVE ERRORS IN THE FUNCTIONS ARE OF THE
!       ORDER OF THE MACHINE PRECISION.

!     DIAG IS AN ARRAY OF LENGTH N. IF MODE = 1 (SEE BELOW),
!       DIAG IS INTERNALLY SET.  IF MODE = 2, DIAG MUST CONTAIN POSITIVE
!       ENTRIES THAT SERVE AS MULTIPLICATIVE SCALE FACTORS FOR THE
!       VARIABLES.

!     MODE IS AN INTEGER INPUT VARIABLE. IF MODE = 1, THE VARIABLES WILL BE
!       SCALED INTERNALLY.  IF MODE = 2, THE SCALING IS SPECIFIED BY THE
!       INPUT DIAG.  OTHER VALUES OF MODE ARE EQUIVALENT TO MODE = 1.

!     FACTOR IS A POSITIVE INPUT VARIABLE USED IN DETERMINING THE
!       INITIAL STEP BOUND. THIS BOUND IS SET TO THE PRODUCT OF
!       FACTOR AND THE EUCLIDEAN NORM OF DIAG*X IF NONZERO, OR ELSE
!       TO FACTOR ITSELF. IN MOST CASES FACTOR SHOULD LIE IN THE
!       INTERVAL (.1,100.). 100. IS A GENERALLY RECOMMENDED VALUE.

!     NPRINT IS AN INTEGER INPUT VARIABLE THAT ENABLES CONTROLLED
!       PRINTING OF ITERATES IF IT IS POSITIVE. IN THIS CASE,
!       FCN IS CALLED WITH IFLAG = 0 AT THE BEGINNING OF THE FIRST
!       ITERATION AND EVERY NPRINT ITERATIONS THEREAFTER AND
!       IMMEDIATELY PRIOR TO RETURN, WITH X AND FVEC AVAILABLE
!       FOR PRINTING. IF NPRINT IS NOT POSITIVE, NO SPECIAL CALLS
!       OF FCN WITH IFLAG = 0 ARE MADE.

!     INFO IS AN INTEGER OUTPUT VARIABLE. IF THE USER HAS
!       TERMINATED EXECUTION, INFO IS SET TO THE (NEGATIVE)
!       VALUE OF IFLAG. SEE DESCRIPTION OF FCN. OTHERWISE,
!       INFO IS SET AS FOLLOWS.

!       INFO = 0   IMPROPER INPUT PARAMETERS.

!       INFO = 1   RELATIVE ERROR BETWEEN TWO CONSECUTIVE ITERATES
!                  IS AT MOST XTOL.

!       INFO = 2   NUMBER OF CALLS TO FCN HAS REACHED OR EXCEEDED MAXFEV.

!       INFO = 3   XTOL IS TOO SMALL. NO FURTHER IMPROVEMENT IN
!                  THE APPROXIMATE SOLUTION X IS POSSIBLE.

!       INFO = 4   ITERATION IS NOT MAKING GOOD PROGRESS, AS
!                  MEASURED BY THE IMPROVEMENT FROM THE LAST
!                  FIVE JACOBIAN EVALUATIONS.

!       INFO = 5   ITERATION IS NOT MAKING GOOD PROGRESS, AS MEASURED BY
!                  THE IMPROVEMENT FROM THE LAST TEN ITERATIONS.

!     NFEV IS AN INTEGER OUTPUT VARIABLE SET TO THE NUMBER OF CALLS TO FCN.

!     FJAC IS AN OUTPUT N BY N ARRAY WHICH CONTAINS THE ORTHOGONAL MATRIX Q
!       PRODUCED BY THE QR FACTORIZATION OF THE FINAL APPROXIMATE JACOBIAN.

!     LDFJAC IS A POSITIVE INTEGER INPUT VARIABLE NOT LESS THAN N
!       WHICH SPECIFIES THE LEADING DIMENSION OF THE ARRAY FJAC.

!     R IS AN OUTPUT ARRAY OF LENGTH LR WHICH CONTAINS THE
!       UPPER TRIANGULAR MATRIX PRODUCED BY THE QR FACTORIZATION
!       OF THE FINAL APPROXIMATE JACOBIAN, STORED ROWWISE.

!     LR IS A POSITIVE INTEGER INPUT VARIABLE NOT LESS THAN (N*(N+1))/2.

!     QTF IS AN OUTPUT ARRAY OF LENGTH N WHICH CONTAINS
!       THE VECTOR (Q TRANSPOSE)*FVEC.

!     WA1, WA2, WA3, AND WA4 ARE WORK ARRAYS OF LENGTH N.

!   SUBPROGRAMS CALLED

!     USER-SUPPLIED ...... FCN

!     MINPACK-SUPPLIED ... DOGLEG,SPMPAR,ENORM,FDJAC1,
!                          QFORM,QRFAC,R1MPYQ,R1UPDT

!     FORTRAN-SUPPLIED ... ABS,MAX,MIN,MIN,MOD

!   ARGONNE NATIONAL LABORATORY. MINPACK PROJECT. MARCH 1980.
!   BURTON S. GARBOW, KENNETH E. HILLSTROM, JORGE J. MORE

!   **********

INTEGER    :: i, iflag, iter, j, jm1, l, lr, msum, ncfail, ncsuc, nslow1,  &
              nslow2
INTEGER    :: iwa(1),print_lev=4
LOGICAL    :: jeval, sing
REAL (dp)  :: actred, delta, epsmch, fnorm, fnorm1, pnorm, prered,   &
              ratio, sum, temp, xnorm
REAL (dp), PARAMETER  :: one = 1.0_dp, p1 = 0.1_dp, p5 = 0.5_dp,   &
                         p001 = 0.001_dp, p0001 = 0.0001_dp, zero = 0.0_dp

! The following were workspace arguments
REAL (dp)  :: fjac(n,n), r(n*(n+1)/2), qtf(n), wa1(n), wa2(n),   &
              wa3(n), wa4(n)

!     EPSMCH IS THE MACHINE PRECISION.

epsmch = EPSILON(1.0_dp)

info = 0
iflag = 0
nfev = 0
lr = n*(n+1)/2

!     CHECK THE INPUT PARAMETERS FOR ERRORS.

IF (n > 0 .AND. xtol >= zero .AND. maxfev > 0 .AND. ml >= 0 .AND. mu >=  &
    0 .AND. factor > zero ) THEN
IF (mode == 2) THEN
  diag(1:n) = one
END IF

!     EVALUATE THE FUNCTION AT THE STARTING POINT AND CALCULATE ITS NORM.

iflag = 1
CALL fcn(n, x, fvec, iflag)
nfev = 1
IF (iflag >= 0) THEN
  fnorm = enorm(n, fvec)
  
!   DETERMINE THE NUMBER OF CALLS TO FCN NEEDED TO COMPUTE THE JACOBIAN MATRIX.
  
  msum = MIN(ml+mu+1,n)
  
!     INITIALIZE ITERATION COUNTER AND MONITORS.
  
  iter = 1
  ncsuc = 0
  ncfail = 0
  nslow1 = 0
  nslow2 = 0
  
!     BEGINNING OF THE OUTER LOOP.
  
  20 jeval = .true.
  
!        CALCULATE THE JACOBIAN MATRIX.
  
  iflag = 2
  CALL fdjac1(fcn, n, x, fvec, fjac, n, iflag, ml, mu, epsfcn, wa1, wa2)
  if(print_lev>=4) call mat2csv(fjac,'fjac.csv')
  
  nfev = nfev + msum
  IF (iflag >= 0) THEN
    
!        COMPUTE THE QR FACTORIZATION OF THE JACOBIAN.
    
    CALL qrfac(n, n, fjac, n, .false., iwa, 1, wa1, wa2, wa3)
    
!        ON THE FIRST ITERATION AND IF MODE IS 1, SCALE ACCORDING
!        TO THE NORMS OF THE COLUMNS OF THE INITIAL JACOBIAN.
    
    IF (iter == 1) THEN
      IF (mode /= 2) THEN
        DO  j = 1, n
          diag(j) = wa2(j)
          IF (wa2(j) == zero) diag(j) = one
        END DO
      END IF
      
!        ON THE FIRST ITERATION, CALCULATE THE NORM OF THE SCALED X
!        AND INITIALIZE THE STEP BOUND DELTA.
      
      wa3(1:n) = diag(1:n) * x(1:n)
      xnorm = enorm(n, wa3)
      delta = factor * xnorm
      IF (delta == zero) delta = factor
    END IF
    
!        FORM (Q TRANSPOSE)*FVEC AND STORE IN QTF.
    
    qtf(1:n) = fvec(1:n)
    DO  j = 1, n
      IF (fjac(j,j) /= zero) THEN
        sum = zero
        DO  i = j, n
          sum = sum + fjac(i,j) * qtf(i)
        END DO
        temp = -sum / fjac(j,j)
        DO  i = j, n
          qtf(i) = qtf(i) + fjac(i,j) * temp
        END DO
      END IF
    END DO
    
!        COPY THE TRIANGULAR FACTOR OF THE QR FACTORIZATION INTO R.
    
    sing = .false.
    DO  j = 1, n
      l = j
      jm1 = j - 1
      IF (jm1 >= 1) THEN
        DO  i = 1, jm1
          r(l) = fjac(i,j)
          l = l + n - i
        END DO
      END IF
      r(l) = wa1(j)
      IF (wa1(j) == zero) sing = .true.
    END DO
    
!        ACCUMULATE THE ORTHOGONAL FACTOR IN FJAC.
    
    CALL qform(n, n, fjac, n, wa1)
    
!        RESCALE IF NECESSARY.
    
    IF (mode /= 2) THEN
      DO  j = 1, n
        diag(j) = MAX(diag(j), wa2(j))
      END DO
    END IF
    
!        BEGINNING OF THE INNER LOOP.
    
!           IF REQUESTED, CALL FCN TO ENABLE PRINTING OF ITERATES.
    
    120 IF (nprint > 0) THEN
      iflag = 0
      IF (MOD(iter-1, nprint) == 0) CALL fcn(n, x, fvec, iflag)
      IF (iflag < 0) GO TO 190
    END IF
    
!           DETERMINE THE DIRECTION P.
    
    CALL dogleg(n, r, lr, diag, qtf, delta, wa1, wa2, wa3)
    
!           STORE THE DIRECTION P AND X + P. CALCULATE THE NORM OF P.
    
    DO  j = 1, n
      wa1(j) = -wa1(j)
      wa2(j) = x(j) + wa1(j)
      wa3(j) = diag(j) * wa1(j)
    END DO
    pnorm = enorm(n, wa3)
    
!           ON THE FIRST ITERATION, ADJUST THE INITIAL STEP BOUND.
    
    IF (iter == 1) delta = MIN(delta, pnorm)
    
!           EVALUATE THE FUNCTION AT X + P AND CALCULATE ITS NORM.
    
    iflag = 1
    CALL fcn(n, wa2, wa4, iflag)
    nfev = nfev + 1
    IF (iflag >= 0) THEN
      fnorm1 = enorm(n, wa4)
      
!           COMPUTE THE SCALED ACTUAL REDUCTION.
      
      actred = -one
      IF (fnorm1 < fnorm) actred = one - (fnorm1/fnorm) ** 2
      
!           COMPUTE THE SCALED PREDICTED REDUCTION.
      
      l = 1
      DO  i = 1, n
        sum = zero
        DO  j = i, n
          sum = sum + r(l) * wa1(j)
          l = l + 1
        END DO
        wa3(i) = qtf(i) + sum
      END DO
      temp = enorm(n, wa3)
      prered = zero
      IF (temp < fnorm) prered = one - (temp/fnorm) ** 2
      
!           COMPUTE THE RATIO OF THE ACTUAL TO THE PREDICTED REDUCTION.
      
      ratio = zero
      IF (prered > zero) ratio = actred / prered
      
!           UPDATE THE STEP BOUND.
      
      IF (ratio < p1) THEN
        ncsuc = 0
        ncfail = ncfail + 1
        delta = p5 * delta
      ELSE
        ncfail = 0
        ncsuc = ncsuc + 1
        IF (ratio >= p5 .OR. ncsuc > 1) delta = MAX(delta,pnorm/p5)
        IF (ABS(ratio-one) <= p1) delta = pnorm / p5
      END IF
      
!           TEST FOR SUCCESSFUL ITERATION.
      
      IF (ratio >= p0001) THEN
        
!           SUCCESSFUL ITERATION. UPDATE X, FVEC, AND THEIR NORMS.
        
        DO  j = 1, n
          x(j) = wa2(j)
          wa2(j) = diag(j) * x(j)
          fvec(j) = wa4(j)
        END DO
        xnorm = enorm(n, wa2)
        fnorm = fnorm1
        iter = iter + 1
      END IF
      
!           DETERMINE THE PROGRESS OF THE ITERATION.
      
      nslow1 = nslow1 + 1
      IF (actred >= p001) nslow1 = 0
      IF (jeval) nslow2 = nslow2 + 1
      IF (actred >= p1) nslow2 = 0
      
!           TEST FOR CONVERGENCE.
      
      IF (delta <= xtol*xnorm .OR. fnorm == zero) info = 1
      IF (info == 0) THEN
        
!           TESTS FOR TERMINATION AND STRINGENT TOLERANCES.
        
        IF (nfev >= maxfev) info = 2
        IF (p1*MAX(p1*delta, pnorm) <= epsmch*xnorm) info = 3
        IF (nslow2 == 5) info = 4
        IF (nslow1 == 10) info = 5
        IF (info == 0) THEN
          
!           CRITERION FOR RECALCULATING JACOBIAN APPROXIMATION
!           BY FORWARD DIFFERENCES.
          
          IF (ncfail /= 2) THEN
            
!           CALCULATE THE RANK ONE MODIFICATION TO THE JACOBIAN
!           AND UPDATE QTF IF NECESSARY.
            
            DO  j = 1, n
              sum = zero
              DO  i = 1, n
                sum = sum + fjac(i,j) * wa4(i)
              END DO
              wa2(j) = (sum-wa3(j)) / pnorm
              wa1(j) = diag(j) * ((diag(j)*wa1(j))/pnorm)
              IF (ratio >= p0001) qtf(j) = sum
            END DO
            
!           COMPUTE THE QR FACTORIZATION OF THE UPDATED JACOBIAN.
            
            CALL r1updt(n, n, r, lr, wa1, wa2, wa3, sing)
            CALL r1mpyq(n, n, fjac, n, wa2, wa3)
            CALL r1mpyq(1, n, qtf, 1, wa2, wa3)
            
!           END OF THE INNER LOOP.
            
            jeval = .false.
            GO TO 120
          END IF
          
!        END OF THE OUTER LOOP.
          
          GO TO 20
        END IF
      END IF
    END IF
  END IF
END IF
END IF

!     TERMINATION, EITHER NORMAL OR USER IMPOSED.

190 IF (iflag < 0) info = iflag
iflag = 0
IF (nprint > 0) CALL fcn(n, x, fvec, iflag)
RETURN

!     LAST CARD OF SUBROUTINE HYBRD.

END SUBROUTINE hybrd



SUBROUTINE dogleg(n, r, lr, diag, qtb, delta, x, wa1, wa2)

INTEGER, INTENT(IN)        :: n
INTEGER, INTENT(IN)        :: lr
REAL (dp), INTENT(IN)      :: r(lr)
REAL (dp), INTENT(IN)      :: diag(n)
REAL (dp), INTENT(IN)      :: qtb(n)
REAL (dp), INTENT(IN)      :: delta
REAL (dp), INTENT(IN OUT)  :: x(n)
REAL (dp), INTENT(OUT)     :: wa1(n)
REAL (dp), INTENT(OUT)     :: wa2(n)


!     **********

!     SUBROUTINE DOGLEG

!     GIVEN AN M BY N MATRIX A, AN N BY N NONSINGULAR DIAGONAL
!     MATRIX D, AN M-VECTOR B, AND A POSITIVE NUMBER DELTA, THE
!     PROBLEM IS TO DETERMINE THE CONVEX COMBINATION X OF THE
!     GAUSS-NEWTON AND SCALED GRADIENT DIRECTIONS THAT MINIMIZES
!     (A*X - B) IN THE LEAST SQUARES SENSE, SUBJECT TO THE
!     RESTRICTION THAT THE EUCLIDEAN NORM OF D*X BE AT MOST DELTA.

!     THIS SUBROUTINE COMPLETES THE SOLUTION OF THE PROBLEM
!     IF IT IS PROVIDED WITH THE NECESSARY INFORMATION FROM THE
!     QR FACTORIZATION OF A. THAT IS, IF A = Q*R, WHERE Q HAS
!     ORTHOGONAL COLUMNS AND R IS AN UPPER TRIANGULAR MATRIX,
!     THEN DOGLEG EXPECTS THE FULL UPPER TRIANGLE OF R AND
!     THE FIRST N COMPONENTS OF (Q TRANSPOSE)*B.

!     THE SUBROUTINE STATEMENT IS

!       SUBROUTINE DOGLEG(N,R,LR,DIAG,QTB,DELTA,X,WA1,WA2)

!     WHERE

!       N IS A POSITIVE INTEGER INPUT VARIABLE SET TO THE ORDER OF R.

!       R IS AN INPUT ARRAY OF LENGTH LR WHICH MUST CONTAIN THE UPPER
!         TRIANGULAR MATRIX R STORED BY ROWS.

!       LR IS A POSITIVE INTEGER INPUT VARIABLE NOT LESS THAN
!         (N*(N+1))/2.

!       DIAG IS AN INPUT ARRAY OF LENGTH N WHICH MUST CONTAIN THE
!         DIAGONAL ELEMENTS OF THE MATRIX D.

!       QTB IS AN INPUT ARRAY OF LENGTH N WHICH MUST CONTAIN THE FIRST
!         N ELEMENTS OF THE VECTOR (Q TRANSPOSE)*B.

!       DELTA IS A POSITIVE INPUT VARIABLE WHICH SPECIFIES AN UPPER
!         BOUND ON THE EUCLIDEAN NORM OF D*X.

!       X IS AN OUTPUT ARRAY OF LENGTH N WHICH CONTAINS THE DESIRED
!         CONVEX COMBINATION OF THE GAUSS-NEWTON DIRECTION AND THE
!         SCALED GRADIENT DIRECTION.

!       WA1 AND WA2 ARE WORK ARRAYS OF LENGTH N.

!     SUBPROGRAMS CALLED

!       MINPACK-SUPPLIED ... SPMPAR,ENORM

!       FORTRAN-SUPPLIED ... ABS,MAX,MIN,SQRT

!     ARGONNE NATIONAL LABORATORY. MINPACK PROJECT. MARCH 1980.
!     BURTON S. GARBOW, KENNETH E. HILLSTROM, JORGE J. MORE

!     **********
INTEGER    :: i, j, jj, jp1, k, l
REAL (dp)  :: alpha, bnorm, epsmch, gnorm, qnorm, sgnorm, sum, temp

!     EPSMCH IS THE MACHINE PRECISION.

epsmch = EPSILON(1.0_dp)

!     FIRST, CALCULATE THE GAUSS-NEWTON DIRECTION.

jj = (n*(n+1)) / 2 + 1
DO  k = 1, n
  j = n - k + 1
  jp1 = j + 1
  jj = jj - k
  l = jj + 1
  sum = 0.0
  IF (n >= jp1) THEN
    DO  i = jp1, n
      sum = sum + r(l) * x(i)
      l = l + 1
    END DO
  END IF
  temp = r(jj)
  IF (temp == 0.0) THEN
    l = j
    DO  i = 1, j
      temp = MAX(temp,ABS(r(l)))
      l = l + n - i
    END DO
    temp = epsmch * temp
    IF (temp == 0.0) temp = epsmch
  END IF
  x(j) = (qtb(j)-sum) / temp
END DO

!     TEST WHETHER THE GAUSS-NEWTON DIRECTION IS ACCEPTABLE.

DO  j = 1, n
  wa1(j) = 0.0
  wa2(j) = diag(j) * x(j)
END DO
qnorm = enorm(n, wa2)
IF (qnorm > delta) THEN
  
!     THE GAUSS-NEWTON DIRECTION IS NOT ACCEPTABLE.
!     NEXT, CALCULATE THE SCALED GRADIENT DIRECTION.
  
  l = 1
  DO  j = 1, n
    temp = qtb(j)
    DO  i = j, n
      wa1(i) = wa1(i) + r(l) * temp
      l = l + 1
    END DO
    wa1(j) = wa1(j) / diag(j)
  END DO
  
!     CALCULATE THE NORM OF THE SCALED GRADIENT AND TEST FOR
!     THE SPECIAL CASE IN WHICH THE SCALED GRADIENT IS ZERO.
  
  gnorm = enorm(n, wa1)
  sgnorm = 0.0
  alpha = delta / qnorm
  IF (gnorm /= 0.0) THEN
    
!     CALCULATE THE POINT ALONG THE SCALED GRADIENT
!     AT WHICH THE QUADRATIC IS MINIMIZED.
    
    DO  j = 1, n
      wa1(j) = (wa1(j)/gnorm) / diag(j)
    END DO
    l = 1
    DO  j = 1, n
      sum = 0.0
      DO  i = j, n
        sum = sum + r(l) * wa1(i)
        l = l + 1
      END DO
      wa2(j) = sum
    END DO
    temp = enorm(n, wa2)
    sgnorm = (gnorm/temp) / temp
    
!     TEST WHETHER THE SCALED GRADIENT DIRECTION IS ACCEPTABLE.
    
    alpha = 0.0
    IF (sgnorm < delta) THEN
      
!     THE SCALED GRADIENT DIRECTION IS NOT ACCEPTABLE.
!     FINALLY, CALCULATE THE POINT ALONG THE DOGLEG
!     AT WHICH THE QUADRATIC IS MINIMIZED.
      
      bnorm = enorm(n, qtb)
      temp = (bnorm/gnorm) * (bnorm/qnorm) * (sgnorm/delta)
      temp = temp - (delta/qnorm) * (sgnorm/delta) ** 2 + SQRT((  &
          temp-(delta/qnorm))**2+(1.0-(delta/qnorm)**2)*(1.0-( sgnorm/delta)**2))
      alpha = ((delta/qnorm)*(1.0-(sgnorm/delta)**2)) / temp
    END IF
  END IF
  
!     FORM APPROPRIATE CONVEX COMBINATION OF THE GAUSS-NEWTON
!     DIRECTION AND THE SCALED GRADIENT DIRECTION.
  
  temp = (1.0-alpha) * MIN(sgnorm,delta)
  DO  j = 1, n
    x(j) = temp * wa1(j) + alpha * x(j)
  END DO
END IF
RETURN
END SUBROUTINE dogleg


SUBROUTINE fdjac1(fcn, n, x, fvec, fjac, ldfjac, iflag, ml, mu, epsfcn,   &
                  wa1, wa2)

INTEGER, INTENT(IN)        :: n
REAL (dp), INTENT(IN OUT)  :: x(n)
REAL (dp), INTENT(IN)      :: fvec(n)
INTEGER, INTENT(IN)        :: ldfjac
REAL (dp), INTENT(OUT)     :: fjac(ldfjac,n)
INTEGER, INTENT(IN OUT)    :: iflag
INTEGER, INTENT(IN)        :: ml
INTEGER, INTENT(IN)        :: mu
REAL (dp), INTENT(IN)      :: epsfcn
REAL (dp), INTENT(IN OUT)  :: wa1(n)
REAL (dp), INTENT(OUT)     :: wa2(n)

! EXTERNAL fcn
INTERFACE
  SUBROUTINE FCN(N, X, FVEC, IFLAG)
    IMPLICIT NONE
    INTEGER, PARAMETER  :: dp = SELECTED_REAL_KIND(14, 60)
    INTEGER, INTENT(IN)      :: n
    REAL (dp), INTENT(IN)    :: x(n)
    REAL (dp), INTENT(OUT)   :: fvec(n)
    INTEGER, INTENT(IN OUT)  :: iflag
  END SUBROUTINE FCN
END INTERFACE

!   **********

!   SUBROUTINE FDJAC1

!   THIS SUBROUTINE COMPUTES A FORWARD-DIFFERENCE APPROXIMATION TO THE N BY N
!   JACOBIAN MATRIX ASSOCIATED WITH A SPECIFIED PROBLEM OF N FUNCTIONS IN N
!   VARIABLES.  IF THE JACOBIAN HAS A BANDED FORM, THEN FUNCTION EVALUATIONS
!   ARE SAVED BY ONLY APPROXIMATING THE NONZERO TERMS.

!   THE SUBROUTINE STATEMENT IS

!     SUBROUTINE FDJAC1(FCN,N,X,FVEC,FJAC,LDFJAC,IFLAG,ML,MU,EPSFCN,
!                       WA1,WA2)

!   WHERE

!     FCN IS THE NAME OF THE USER-SUPPLIED SUBROUTINE WHICH CALCULATES
!       THE FUNCTIONS.  FCN MUST BE DECLARED IN AN EXTERNAL STATEMENT IN
!       THE USER CALLING PROGRAM, AND SHOULD BE WRITTEN AS FOLLOWS.

!       SUBROUTINE FCN(N,X,FVEC,IFLAG)
!       INTEGER N,IFLAG
!       REAL X(N),FVEC(N)
!       ----------
!       CALCULATE THE FUNCTIONS AT X AND
!       RETURN THIS VECTOR IN FVEC.
!       ----------
!       RETURN
!       END

!       THE VALUE OF IFLAG SHOULD NOT BE CHANGED BY FCN UNLESS
!       THE USER WANTS TO TERMINATE EXECUTION OF FDJAC1.
!       IN THIS CASE SET IFLAG TO A NEGATIVE INTEGER.

!     N IS A POSITIVE INTEGER INPUT VARIABLE SET TO THE NUMBER
!       OF FUNCTIONS AND VARIABLES.

!     X IS AN INPUT ARRAY OF LENGTH N.

!     FVEC IS AN INPUT ARRAY OF LENGTH N WHICH MUST CONTAIN THE
!       FUNCTIONS EVALUATED AT X.

!     FJAC IS AN OUTPUT N BY N ARRAY WHICH CONTAINS THE
!       APPROXIMATION TO THE JACOBIAN MATRIX EVALUATED AT X.

!     LDFJAC IS A POSITIVE INTEGER INPUT VARIABLE NOT LESS THAN N
!       WHICH SPECIFIES THE LEADING DIMENSION OF THE ARRAY FJAC.

!     IFLAG IS AN INTEGER VARIABLE WHICH CAN BE USED TO TERMINATE
!       THE EXECUTION OF FDJAC1.  SEE DESCRIPTION OF FCN.

!     ML IS A NONNEGATIVE INTEGER INPUT VARIABLE WHICH SPECIFIES
!       THE NUMBER OF SUBDIAGONALS WITHIN THE BAND OF THE
!       JACOBIAN MATRIX. IF THE JACOBIAN IS NOT BANDED, SET
!       ML TO AT LEAST N - 1.

!     EPSFCN IS AN INPUT VARIABLE USED IN DETERMINING A SUITABLE
!       STEP LENGTH FOR THE FORWARD-DIFFERENCE APPROXIMATION. THIS
!       APPROXIMATION ASSUMES THAT THE RELATIVE ERRORS IN THE
!       FUNCTIONS ARE OF THE ORDER OF EPSFCN. IF EPSFCN IS LESS
!       THAN THE MACHINE PRECISION, IT IS ASSUMED THAT THE RELATIVE
!       ERRORS IN THE FUNCTIONS ARE OF THE ORDER OF THE MACHINE PRECISION.

!     MU IS A NONNEGATIVE INTEGER INPUT VARIABLE WHICH SPECIFIES
!       THE NUMBER OF SUPERDIAGONALS WITHIN THE BAND OF THE
!       JACOBIAN MATRIX. IF THE JACOBIAN IS NOT BANDED, SET
!       MU TO AT LEAST N - 1.

!     WA1 AND WA2 ARE WORK ARRAYS OF LENGTH N.  IF ML + MU + 1 IS AT
!       LEAST N, THEN THE JACOBIAN IS CONSIDERED DENSE, AND WA2 IS
!       NOT REFERENCED.

!   SUBPROGRAMS CALLED

!     MINPACK-SUPPLIED ... SPMPAR

!     FORTRAN-SUPPLIED ... ABS,MAX,SQRT

!   ARGONNE NATIONAL LABORATORY. MINPACK PROJECT. MARCH 1980.
!   BURTON S. GARBOW, KENNETH E. HILLSTROM, JORGE J. MORE

!   **********
INTEGER    :: i, j, k, msum
REAL (dp)  :: eps, epsmch, h, temp
REAL (dp), PARAMETER  :: zero = 0.0_dp

!     EPSMCH IS THE MACHINE PRECISION.

epsmch = EPSILON(1.0_dp)

eps = SQRT(MAX(epsfcn, epsmch))
msum = ml + mu + 1
IF (msum >= n) THEN
  
!        COMPUTATION OF DENSE APPROXIMATE JACOBIAN.
  
  DO  j = 1, n
    temp = x(j)
    h = eps * ABS(temp)
    IF (h == zero) h = eps
    x(j) = temp + h
    CALL fcn(n, x, wa1, iflag)
    IF (iflag < 0) EXIT
    x(j) = temp
    DO  i = 1, n
      fjac(i,j) = (wa1(i)-fvec(i)) / h
    END DO
  END DO
ELSE
  
!        COMPUTATION OF BANDED APPROXIMATE JACOBIAN.
  
  DO  k = 1, msum
    DO  j = k, n, msum
      wa2(j) = x(j)
      h = eps * ABS(wa2(j))
      IF (h == zero) h = eps
      x(j) = wa2(j) + h
    END DO
    CALL fcn(n, x, wa1, iflag)
    IF (iflag < 0) EXIT
    DO  j = k, n, msum
      x(j) = wa2(j)
      h = eps * ABS(wa2(j))
      IF (h == zero) h = eps
      DO  i = 1, n
        fjac(i,j) = zero
        IF (i >= j-mu .AND. i <= j+ml) fjac(i,j) = (wa1(i)-fvec(i)) / h
      END DO
    END DO
  END DO
END IF
RETURN

!     LAST CARD OF SUBROUTINE FDJAC1.

END SUBROUTINE fdjac1



SUBROUTINE qform(m, n, q, ldq, wa)

INTEGER, INTENT(IN)     :: m
INTEGER, INTENT(IN)     :: n
INTEGER, INTENT(IN)     :: ldq
REAL (dp), INTENT(OUT)  :: q(ldq,m)
REAL (dp), INTENT(OUT)  :: wa(m)


!   **********

!   SUBROUTINE QFORM

!   THIS SUBROUTINE PROCEEDS FROM THE COMPUTED QR FACTORIZATION OF AN M BY N
!   MATRIX A TO ACCUMULATE THE M BY M ORTHOGONAL MATRIX Q FROM ITS FACTORED FORM.

!   THE SUBROUTINE STATEMENT IS

!     SUBROUTINE QFORM(M,N,Q,LDQ,WA)

!   WHERE

!     M IS A POSITIVE INTEGER INPUT VARIABLE SET TO THE NUMBER
!       OF ROWS OF A AND THE ORDER OF Q.

!     N IS A POSITIVE INTEGER INPUT VARIABLE SET TO THE NUMBER OF COLUMNS OF A.

!     Q IS AN M BY M ARRAY. ON INPUT THE FULL LOWER TRAPEZOID IN
!       THE FIRST MIN(M,N) COLUMNS OF Q CONTAINS THE FACTORED FORM.
!       ON OUTPUT Q HAS BEEN ACCUMULATED INTO A SQUARE MATRIX.

!     LDQ IS A POSITIVE INTEGER INPUT VARIABLE NOT LESS THAN M
!       WHICH SPECIFIES THE LEADING DIMENSION OF THE ARRAY Q.

!     WA IS A WORK ARRAY OF LENGTH M.

!   SUBPROGRAMS CALLED

!     FORTRAN-SUPPLIED ... MIN

!   ARGONNE NATIONAL LABORATORY. MINPACK PROJECT. MARCH 1980.
!   BURTON S. GARBOW, KENNETH E. HILLSTROM, JORGE J. MORE

!   **********
INTEGER    :: i, j, jm1, k, l, minmn, np1
REAL (dp)  :: sum, temp
REAL (dp), PARAMETER  :: one = 1.0_dp, zero = 0.0_dp

!     ZERO OUT UPPER TRIANGLE OF Q IN THE FIRST MIN(M,N) COLUMNS.

minmn = MIN(m,n)
IF (minmn >= 2) THEN
  DO  j = 2, minmn
    jm1 = j - 1
    DO  i = 1, jm1
      q(i,j) = zero
    END DO
  END DO
END IF

!     INITIALIZE REMAINING COLUMNS TO THOSE OF THE IDENTITY MATRIX.

np1 = n + 1
IF (m >= np1) THEN
  DO  j = np1, m
    DO  i = 1, m
      q(i,j) = zero
    END DO
    q(j,j) = one
  END DO
END IF

!     ACCUMULATE Q FROM ITS FACTORED FORM.

DO  l = 1, minmn
  k = minmn - l + 1
  DO  i = k, m
    wa(i) = q(i,k)
    q(i,k) = zero
  END DO
  q(k,k) = one
  IF (wa(k) /= zero) THEN
    DO  j = k, m
      sum = zero
      DO  i = k, m
        sum = sum + q(i,j) * wa(i)
      END DO
      temp = sum / wa(k)
      DO  i = k, m
        q(i,j) = q(i,j) - temp * wa(i)
      END DO
    END DO
  END IF
END DO
RETURN

!     LAST CARD OF SUBROUTINE QFORM.

END SUBROUTINE qform


SUBROUTINE qrfac(m, n, a, lda, pivot, ipvt, lipvt, rdiag, acnorm, wa)

INTEGER, INTENT(IN)        :: m
INTEGER, INTENT(IN)        :: n
INTEGER, INTENT(IN)        :: lda
REAL (dp), INTENT(IN OUT)  :: a(lda,n)
LOGICAL, INTENT(IN)        :: pivot
INTEGER, INTENT(IN)        :: lipvt
INTEGER, INTENT(OUT)       :: ipvt(lipvt)
REAL (dp), INTENT(OUT)     :: rdiag(n)
REAL (dp), INTENT(OUT)     :: acnorm(n)
REAL (dp), INTENT(OUT)     :: wa(n)


!   **********

!   SUBROUTINE QRFAC

!   THIS SUBROUTINE USES HOUSEHOLDER TRANSFORMATIONS WITH COLUMN PIVOTING
!   (OPTIONAL) TO COMPUTE A QR FACTORIZATION OF THE M BY N MATRIX A.
!   THAT IS, QRFAC DETERMINES AN ORTHOGONAL MATRIX Q, A PERMUTATION MATRIX P,
!   AND AN UPPER TRAPEZOIDAL MATRIX R WITH DIAGONAL ELEMENTS OF NONINCREASING
!   MAGNITUDE, SUCH THAT A*P = Q*R.  THE HOUSEHOLDER TRANSFORMATION FOR
!   COLUMN K, K = 1,2,...,MIN(M,N), IS OF THE FORM

!                         T
!         I - (1/U(K))*U*U

!   WHERE U HAS ZEROS IN THE FIRST K-1 POSITIONS.  THE FORM OF THIS
!   TRANSFORMATION AND THE METHOD OF PIVOTING FIRST APPEARED IN THE
!   CORRESPONDING LINPACK SUBROUTINE.

!   THE SUBROUTINE STATEMENT IS

!     SUBROUTINE QRFAC(M,N,A,LDA,PIVOT,IPVT,LIPVT,RDIAG,ACNORM,WA)

!   WHERE

!     M IS A POSITIVE INTEGER INPUT VARIABLE SET TO THE NUMBER OF ROWS OF A.

!     N IS A POSITIVE INTEGER INPUT VARIABLE SET TO THE NUMBER
!       OF COLUMNS OF A.

!     A IS AN M BY N ARRAY.  ON INPUT A CONTAINS THE MATRIX FOR WHICH THE
!       QR FACTORIZATION IS TO BE COMPUTED.  ON OUTPUT THE STRICT UPPER
!       TRAPEZOIDAL PART OF A CONTAINS THE STRICT UPPER TRAPEZOIDAL PART OF R,
!       AND THE LOWER TRAPEZOIDAL PART OF A CONTAINS A FACTORED FORM OF Q
!       (THE NON-TRIVIAL ELEMENTS OF THE U VECTORS DESCRIBED ABOVE).

!     LDA IS A POSITIVE INTEGER INPUT VARIABLE NOT LESS THAN M
!       WHICH SPECIFIES THE LEADING DIMENSION OF THE ARRAY A.

!     PIVOT IS A LOGICAL INPUT VARIABLE.  IF PIVOT IS SET TRUE,
!       THEN COLUMN PIVOTING IS ENFORCED.  IF PIVOT IS SET FALSE,
!       THEN NO COLUMN PIVOTING IS DONE.

!     IPVT IS AN INTEGER OUTPUT ARRAY OF LENGTH LIPVT.  IPVT DEFINES THE
!       PERMUTATION MATRIX P SUCH THAT A*P = Q*R.
!       COLUMN J OF P IS COLUMN IPVT(J) OF THE IDENTITY MATRIX.
!       IF PIVOT IS FALSE, IPVT IS NOT REFERENCED.

!     LIPVT IS A POSITIVE INTEGER INPUT VARIABLE.  IF PIVOT IS FALSE,
!       THEN LIPVT MAY BE AS SMALL AS 1.  IF PIVOT IS TRUE, THEN
!       LIPVT MUST BE AT LEAST N.

!     RDIAG IS AN OUTPUT ARRAY OF LENGTH N WHICH CONTAINS THE
!       DIAGONAL ELEMENTS OF R.

!     ACNORM IS AN OUTPUT ARRAY OF LENGTH N WHICH CONTAINS THE NORMS OF
!       THE CORRESPONDING COLUMNS OF THE INPUT MATRIX A.
!       IF THIS INFORMATION IS NOT NEEDED, THEN ACNORM CAN COINCIDE WITH RDIAG.

!     WA IS A WORK ARRAY OF LENGTH N. IF PIVOT IS FALSE, THEN WA
!       CAN COINCIDE WITH RDIAG.

!   SUBPROGRAMS CALLED

!     MINPACK-SUPPLIED ... SPMPAR,ENORM

!     FORTRAN-SUPPLIED ... MAX,SQRT,MIN

!   ARGONNE NATIONAL LABORATORY. MINPACK PROJECT. MARCH 1980.
!   BURTON S. GARBOW, KENNETH E. HILLSTROM, JORGE J. MORE

!   **********
INTEGER    :: i, j, jp1, k, kmax, minmn
REAL (dp)  :: ajnorm, epsmch, sum, temp
REAL (dp), PARAMETER  :: one = 1.0_dp, p05 = 0.05_dp, zero = 0.0_dp

!     EPSMCH IS THE MACHINE PRECISION.

epsmch = EPSILON(1.0_dp)

!     COMPUTE THE INITIAL COLUMN NORMS AND INITIALIZE SEVERAL ARRAYS.

DO  j = 1, n
  acnorm(j) = enorm(m, a(1:,j))
  rdiag(j) = acnorm(j)
  wa(j) = rdiag(j)
  IF (pivot) ipvt(j) = j
END DO

!     REDUCE A TO R WITH HOUSEHOLDER TRANSFORMATIONS.

minmn = MIN(m,n)
DO  j = 1, minmn
  IF (pivot) THEN
    
!        BRING THE COLUMN OF LARGEST NORM INTO THE PIVOT POSITION.
    
    kmax = j
    DO  k = j, n
      IF (rdiag(k) > rdiag(kmax)) kmax = k
    END DO
    IF (kmax /= j) THEN
      DO  i = 1, m
        temp = a(i,j)
        a(i,j) = a(i,kmax)
        a(i,kmax) = temp
      END DO
      rdiag(kmax) = rdiag(j)
      wa(kmax) = wa(j)
      k = ipvt(j)
      ipvt(j) = ipvt(kmax)
      ipvt(kmax) = k
    END IF
  END IF
  
!        COMPUTE THE HOUSEHOLDER TRANSFORMATION TO REDUCE THE
!        J-TH COLUMN OF A TO A MULTIPLE OF THE J-TH UNIT VECTOR.
  
  ajnorm = enorm(m-j+1, a(j:,j))
  IF (ajnorm /= zero) THEN
    IF (a(j,j) < zero) ajnorm = -ajnorm
    DO  i = j, m
      a(i,j) = a(i,j) / ajnorm
    END DO
    a(j,j) = a(j,j) + one
    
!        APPLY THE TRANSFORMATION TO THE REMAINING COLUMNS AND UPDATE THE NORMS.
    
    jp1 = j + 1
    IF (n >= jp1) THEN
      DO  k = jp1, n
        sum = zero
        DO  i = j, m
          sum = sum + a(i,j) * a(i,k)
        END DO
        temp = sum / a(j,j)
        DO  i = j, m
          a(i,k) = a(i,k) - temp * a(i,j)
        END DO
        IF (.NOT.(.NOT.pivot.OR.rdiag(k) == zero)) THEN
          temp = a(j,k) / rdiag(k)
          rdiag(k) = rdiag(k) * SQRT(MAX(zero,one-temp**2))
          IF (p05*(rdiag(k)/wa(k))**2 <= epsmch) THEN
            rdiag(k) = enorm(m-j, a(jp1:,k))
            wa(k) = rdiag(k)
          END IF
        END IF
      END DO
    END IF
  END IF
  rdiag(j) = -ajnorm
END DO
RETURN

!     LAST CARD OF SUBROUTINE QRFAC.

END SUBROUTINE qrfac



SUBROUTINE r1mpyq(m, n, a, lda, v, w)

INTEGER, INTENT(IN)        :: m
INTEGER, INTENT(IN)        :: n
INTEGER, INTENT(IN)        :: lda
REAL (dp), INTENT(IN OUT)  :: a(lda,n)
REAL (dp), INTENT(IN)      :: v(n)
REAL (dp), INTENT(IN)      :: w(n)


!   **********

!   SUBROUTINE R1MPYQ

!   GIVEN AN M BY N MATRIX A, THIS SUBROUTINE COMPUTES A*Q WHERE
!   Q IS THE PRODUCT OF 2*(N - 1) TRANSFORMATIONS

!         GV(N-1)*...*GV(1)*GW(1)*...*GW(N-1)

!   AND GV(I), GW(I) ARE GIVENS ROTATIONS IN THE (I,N) PLANE WHICH
!   ELIMINATE ELEMENTS IN THE I-TH AND N-TH PLANES, RESPECTIVELY.
!   Q ITSELF IS NOT GIVEN, RATHER THE INFORMATION TO RECOVER THE
!   GV, GW ROTATIONS IS SUPPLIED.

!   THE SUBROUTINE STATEMENT IS

!     SUBROUTINE R1MPYQ(M, N, A, LDA, V, W)

!   WHERE

!     M IS A POSITIVE INTEGER INPUT VARIABLE SET TO THE NUMBER OF ROWS OF A.

!     N IS A POSITIVE INTEGER INPUT VARIABLE SET TO THE NUMBER OF COLUMNS OF A.

!     A IS AN M BY N ARRAY.  ON INPUT A MUST CONTAIN THE MATRIX TO BE
!       POSTMULTIPLIED BY THE ORTHOGONAL MATRIX Q DESCRIBED ABOVE.
!       ON OUTPUT A*Q HAS REPLACED A.

!     LDA IS A POSITIVE INTEGER INPUT VARIABLE NOT LESS THAN M
!       WHICH SPECIFIES THE LEADING DIMENSION OF THE ARRAY A.

!     V IS AN INPUT ARRAY OF LENGTH N. V(I) MUST CONTAIN THE INFORMATION
!       NECESSARY TO RECOVER THE GIVENS ROTATION GV(I) DESCRIBED ABOVE.

!     W IS AN INPUT ARRAY OF LENGTH N. W(I) MUST CONTAIN THE INFORMATION
!       NECESSARY TO RECOVER THE GIVENS ROTATION GW(I) DESCRIBED ABOVE.

!   SUBROUTINES CALLED

!     FORTRAN-SUPPLIED ... ABS, SQRT

!   ARGONNE NATIONAL LABORATORY. MINPACK PROJECT. MARCH 1980.
!   BURTON S. GARBOW, KENNETH E. HILLSTROM, JORGE J. MORE

!   **********
INTEGER    :: i, j, nmj, nm1
REAL (dp)  :: COS, SIN, temp
REAL (dp), PARAMETER  :: one = 1.0_dp

!     APPLY THE FIRST SET OF GIVENS ROTATIONS TO A.

nm1 = n - 1
IF (nm1 >= 1) THEN
  DO  nmj = 1, nm1
    j = n - nmj
    IF (ABS(v(j)) > one) COS = one / v(j)
    IF (ABS(v(j)) > one) SIN = SQRT(one-COS**2)
    IF (ABS(v(j)) <= one) SIN = v(j)
    IF (ABS(v(j)) <= one) COS = SQRT(one-SIN**2)
    DO  i = 1, m
      temp = COS * a(i,j) - SIN * a(i,n)
      a(i,n) = SIN * a(i,j) + COS * a(i,n)
      a(i,j) = temp
    END DO
  END DO
  
!     APPLY THE SECOND SET OF GIVENS ROTATIONS TO A.
  
  DO  j = 1, nm1
    IF (ABS(w(j)) > one) COS = one / w(j)
    IF (ABS(w(j)) > one) SIN = SQRT(one-COS**2)
    IF (ABS(w(j)) <= one) SIN = w(j)
    IF (ABS(w(j)) <= one) COS = SQRT(one-SIN**2)
    DO  i = 1, m
      temp = COS * a(i,j) + SIN * a(i,n)
      a(i,n) = -SIN * a(i,j) + COS * a(i,n)
      a(i,j) = temp
    END DO
  END DO
END IF
RETURN

!     LAST CARD OF SUBROUTINE R1MPYQ.

END SUBROUTINE r1mpyq



SUBROUTINE r1updt(m, n, s, ls, u, v, w, sing)

INTEGER, INTENT(IN)        :: m
INTEGER, INTENT(IN)        :: n
INTEGER, INTENT(IN)        :: ls
REAL (dp), INTENT(IN OUT)  :: s(ls)
REAL (dp), INTENT(IN)      :: u(m)
REAL (dp), INTENT(IN OUT)  :: v(n)
REAL (dp), INTENT(OUT)     :: w(m)
LOGICAL, INTENT(OUT)       :: sing


!   **********

!   SUBROUTINE R1UPDT

!   GIVEN AN M BY N LOWER TRAPEZOIDAL MATRIX S, AN M-VECTOR U,
!   AND AN N-VECTOR V, THE PROBLEM IS TO DETERMINE AN
!   ORTHOGONAL MATRIX Q SUCH THAT

!                 T
!         (S + U*V )*Q

!   IS AGAIN LOWER TRAPEZOIDAL.

!   THIS SUBROUTINE DETERMINES Q AS THE PRODUCT OF 2*(N - 1) TRANSFORMATIONS

!         GV(N-1)*...*GV(1)*GW(1)*...*GW(N-1)

!   WHERE GV(I), GW(I) ARE GIVENS ROTATIONS IN THE (I,N) PLANE
!   WHICH ELIMINATE ELEMENTS IN THE I-TH AND N-TH PLANES, RESPECTIVELY.
!   Q ITSELF IS NOT ACCUMULATED, RATHER THE INFORMATION TO RECOVER THE GV,
!   GW ROTATIONS IS RETURNED.

!   THE SUBROUTINE STATEMENT IS

!     SUBROUTINE R1UPDT(M,N,S,LS,U,V,W,SING)

!   WHERE

!     M IS A POSITIVE INTEGER INPUT VARIABLE SET TO THE NUMBER OF ROWS OF S.

!     N IS A POSITIVE INTEGER INPUT VARIABLE SET TO THE NUMBER
!       OF COLUMNS OF S.  N MUST NOT EXCEED M.

!     S IS AN ARRAY OF LENGTH LS. ON INPUT S MUST CONTAIN THE LOWER
!       TRAPEZOIDAL MATRIX S STORED BY COLUMNS. ON OUTPUT S CONTAINS
!       THE LOWER TRAPEZOIDAL MATRIX PRODUCED AS DESCRIBED ABOVE.

!     LS IS A POSITIVE INTEGER INPUT VARIABLE NOT LESS THAN
!       (N*(2*M-N+1))/2.

!     U IS AN INPUT ARRAY OF LENGTH M WHICH MUST CONTAIN THE VECTOR U.

!     V IS AN ARRAY OF LENGTH N. ON INPUT V MUST CONTAIN THE VECTOR V.
!       ON OUTPUT V(I) CONTAINS THE INFORMATION NECESSARY TO
!       RECOVER THE GIVENS ROTATION GV(I) DESCRIBED ABOVE.

!     W IS AN OUTPUT ARRAY OF LENGTH M. W(I) CONTAINS INFORMATION
!       NECESSARY TO RECOVER THE GIVENS ROTATION GW(I) DESCRIBED ABOVE.

!     SING IS A LOGICAL OUTPUT VARIABLE.  SING IS SET TRUE IF ANY OF THE
!       DIAGONAL ELEMENTS OF THE OUTPUT S ARE ZERO.  OTHERWISE SING IS
!       SET FALSE.

!   SUBPROGRAMS CALLED

!     MINPACK-SUPPLIED ... SPMPAR

!     FORTRAN-SUPPLIED ... ABS,SQRT

!   ARGONNE NATIONAL LABORATORY. MINPACK PROJECT. MARCH 1980.
!   BURTON S. GARBOW, KENNETH E. HILLSTROM, JORGE J. MORE, JOHN L. NAZARETH

!   **********
INTEGER    :: i, j, jj, l, nmj, nm1
REAL (dp)  :: COS, cotan, giant, SIN, TAN, tau, temp
REAL (dp), PARAMETER  :: one = 1.0_dp, p5 = 0.5_dp, p25 = 0.25_dp, zero = 0.0_dp

!     GIANT IS THE LARGEST MAGNITUDE.

giant = HUGE(1.0_dp)

!     INITIALIZE THE DIAGONAL ELEMENT POINTER.

jj = (n*(2*m-n+1)) / 2 - (m-n)

!     MOVE THE NONTRIVIAL PART OF THE LAST COLUMN OF S INTO W.

l = jj
DO  i = n, m
  w(i) = s(l)
  l = l + 1
END DO

!     ROTATE THE VECTOR V INTO A MULTIPLE OF THE N-TH UNIT VECTOR
!     IN SUCH A WAY THAT A SPIKE IS INTRODUCED INTO W.

nm1 = n - 1
IF (nm1 >= 1) THEN
  DO  nmj = 1, nm1
    j = n - nmj
    jj = jj - (m-j+1)
    w(j) = zero
    IF (v(j) /= zero) THEN
      
!        DETERMINE A GIVENS ROTATION WHICH ELIMINATES THE J-TH ELEMENT OF V.
      
      IF (ABS(v(n)) < ABS(v(j))) THEN
        cotan = v(n) / v(j)
        SIN = p5 / SQRT(p25+p25*cotan**2)
        COS = SIN * cotan
        tau = one
        IF (ABS(COS)*giant > one) tau = one / COS
      ELSE
        TAN = v(j) / v(n)
        COS = p5 / SQRT(p25+p25*TAN**2)
        SIN = COS * TAN
        tau = SIN
      END IF
      
!        APPLY THE TRANSFORMATION TO V AND STORE THE INFORMATION
!        NECESSARY TO RECOVER THE GIVENS ROTATION.
      
      v(n) = SIN * v(j) + COS * v(n)
      v(j) = tau
      
!        APPLY THE TRANSFORMATION TO S AND EXTEND THE SPIKE IN W.
      
      l = jj
      DO  i = j, m
        temp = COS * s(l) - SIN * w(i)
        w(i) = SIN * s(l) + COS * w(i)
        s(l) = temp
        l = l + 1
      END DO
    END IF
  END DO
END IF

!     ADD THE SPIKE FROM THE RANK 1 UPDATE TO W.

DO  i = 1, m
  w(i) = w(i) + v(n) * u(i)
END DO

!     ELIMINATE THE SPIKE.

sing = .false.
IF (nm1 >= 1) THEN
  DO  j = 1, nm1
    IF (w(j) /= zero) THEN
      
!        DETERMINE A GIVENS ROTATION WHICH ELIMINATES THE
!        J-TH ELEMENT OF THE SPIKE.
      
      IF (ABS(s(jj)) < ABS(w(j))) THEN
        cotan = s(jj) / w(j)
        SIN = p5 / SQRT(p25 + p25*cotan**2)
        COS = SIN * cotan
        tau = one
        IF (ABS(COS)*giant > one) tau = one / COS
      ELSE
        TAN = w(j) / s(jj)
        COS = p5 / SQRT(p25+p25*TAN**2)
        SIN = COS * TAN
        tau = SIN
      END IF
      
!        APPLY THE TRANSFORMATION TO S AND REDUCE THE SPIKE IN W.
      
      l = jj
      DO  i = j, m
        temp = COS * s(l) + SIN * w(i)
        w(i) = -SIN * s(l) + COS * w(i)
        s(l) = temp
        l = l + 1
      END DO
      
!        STORE THE INFORMATION NECESSARY TO RECOVER THE GIVENS ROTATION.
      
      w(j) = tau
    END IF
    
!        TEST FOR ZERO DIAGONAL ELEMENTS IN THE OUTPUT S.
    
    IF (s(jj) == zero) sing = .true.
    jj = jj + (m-j+1)
  END DO
END IF

!     MOVE W BACK INTO THE LAST COLUMN OF THE OUTPUT S.

l = jj
DO  i = n, m
  s(l) = w(i)
  l = l + 1
END DO
IF (s(jj) == zero) sing = .true.
RETURN

!     LAST CARD OF SUBROUTINE R1UPDT.

END SUBROUTINE r1updt


FUNCTION enorm(n, x) RESULT(fn_val)
 
INTEGER, INTENT(IN)    :: n
REAL (dp), INTENT(IN)  :: x(n)
REAL (dp)              :: fn_val

!   **********

!   FUNCTION ENORM

!   GIVEN AN N-VECTOR X, THIS FUNCTION CALCULATES THE EUCLIDEAN NORM OF X.

!   THE EUCLIDEAN NORM IS COMPUTED BY ACCUMULATING THE SUM OF SQUARES IN THREE
!   DIFFERENT SUMS.  THE SUMS OF SQUARES FOR THE SMALL AND LARGE COMPONENTS
!   ARE SCALED SO THAT NO OVERFLOWS OCCUR.  NON-DESTRUCTIVE UNDERFLOWS ARE
!   PERMITTED.  UNDERFLOWS AND OVERFLOWS DO NOT OCCUR IN THE COMPUTATION OF THE UNSCALED
!   SUM OF SQUARES FOR THE INTERMEDIATE COMPONENTS.
!   THE DEFINITIONS OF SMALL, INTERMEDIATE AND LARGE COMPONENTS DEPEND ON
!   TWO CONSTANTS, RDWARF AND RGIANT.  THE MAIN RESTRICTIONS ON THESE CONSTANTS
!   ARE THAT RDWARF**2 NOT UNDERFLOW AND RGIANT**2 NOT OVERFLOW.
!   THE CONSTANTS GIVEN HERE ARE SUITABLE FOR EVERY KNOWN COMPUTER.

!   THE FUNCTION STATEMENT IS

!     REAL FUNCTION ENORM(N, X)

!   WHERE

!     N IS A POSITIVE INTEGER INPUT VARIABLE.

!     X IS AN INPUT ARRAY OF LENGTH N.

!   SUBPROGRAMS CALLED

!     FORTRAN-SUPPLIED ... ABS,SQRT

!   ARGONNE NATIONAL LABORATORY. MINPACK PROJECT. MARCH 1980.
!   BURTON S. GARBOW, KENNETH E. HILLSTROM, JORGE J. MORE

!   **********
INTEGER    :: i
REAL (dp)  :: agiant, floatn, s1, s2, s3, xabs, x1max, x3max
REAL (dp), PARAMETER  :: rdwarf = 1.0D-100, rgiant = 1.0D+100

s1 = 0.0_dp
s2 = 0.0_dp
s3 = 0.0_dp
x1max = 0.0_dp
x3max = 0.0_dp
floatn = n
agiant = rgiant / floatn
DO  i = 1, n
  xabs = ABS(x(i))
  IF (xabs <= rdwarf .OR. xabs >= agiant) THEN
    IF (xabs > rdwarf) THEN
      
!              SUM FOR LARGE COMPONENTS.
      
      IF (xabs > x1max) THEN
        s1 = 1.0_dp + s1 * (x1max/xabs) ** 2
        x1max = xabs
      ELSE
        s1 = s1 + (xabs/x1max) ** 2
      END IF
    ELSE
      
!              SUM FOR SMALL COMPONENTS.
      
      IF (xabs > x3max) THEN
        s3 = 1.0_dp + s3 * (x3max/xabs) ** 2
        x3max = xabs
      ELSE
        IF (xabs /= 0.0_dp) s3 = s3 + (xabs/x3max) ** 2
      END IF
    END IF
  ELSE
    
!           SUM FOR INTERMEDIATE COMPONENTS.
    
    s2 = s2 + xabs ** 2
  END IF
END DO

!     CALCULATION OF NORM.

IF (s1 /= 0.0_dp) THEN
  fn_val = x1max * SQRT(s1 + (s2/x1max)/x1max)
ELSE
  IF (s2 /= 0.0_dp) THEN
    IF (s2 >= x3max) fn_val = SQRT(s2*(1.0_dp + (x3max/s2)*(x3max*s3)))
    IF (s2 < x3max) fn_val = SQRT(x3max*((s2/x3max) + (x3max*s3)))
  ELSE
    fn_val = x3max * SQRT(s3)
  END IF
END IF
RETURN
END FUNCTION enorm

END MODULE Solve_NonLin


MODULE random
! A module for random number generation from the following distributions:
!
!     Distribution                    Function/subroutine name
!
!     Normal (Gaussian)               random_normal
!     Gamma                           random_gamma
!     Chi-squared                     random_chisq
!     Exponential                     random_exponential
!     Weibull                         random_Weibull
!     Beta                            random_beta
!     t                               random_t
!     Multivariate normal             random_mvnorm
!     Generalized inverse Gaussian    random_inv_gauss
!     Poisson                         random_Poisson
!     Binomial                        random_binomial1   *
!                                     random_binomial2   *
!     Negative binomial               random_neg_binomial
!     von Mises                       random_von_Mises
!     Cauchy                          random_Cauchy
!
!  Generate a random ordering of the integers 1 .. N
!                                     random_order
!     Initialize (seed) the uniform random number generator for ANY compiler
!                                     seed_random_number

!     Lognormal - see note below.

!  ** Two functions are provided for the binomial distribution.
!  If the parameter values remain constant, it is recommended that the
!  first function is used (random_binomial1).   If one or both of the
!  parameters change, use the second function (random_binomial2).

! The compilers own random number generator, SUBROUTINE RANDOM_NUMBER(r),
! is used to provide a source of uniformly distributed random numbers.

! N.B. At this stage, only one random number is generated at each call to
!      one of the functions above.

! The module uses the following functions which are included here:
! bin_prob to calculate a single binomial probability
! lngamma  to calculate the logarithm to base e of the gamma function

! Some of the code is adapted from Dagpunar's book:
!     Dagpunar, J. 'Principles of random variate generation'
!     Clarendon Press, Oxford, 1988.   ISBN 0-19-852202-9
!
! In most of Dagpunar's routines, there is a test to see whether the value
! of one or two floating-point parameters has changed since the last call.
! These tests have been replaced by using a logical variable FIRST.
! This should be set to .TRUE. on the first call using new values of the
! parameters, and .FALSE. if the parameter values are the same as for the
! previous call.

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Lognormal distribution
! If X has a lognormal distribution, then log(X) is normally distributed.
! Here the logarithm is the natural logarithm, that is to base e, sometimes
! denoted as ln.  To generate random variates from this distribution, generate
! a random deviate from the normal distribution with mean and variance equal
! to the mean and variance of the logarithms of X, then take its exponential.

! Relationship between the mean & variance of log(X) and the mean & variance
! of X, when X has a lognormal distribution.
! Let m = mean of log(X), and s^2 = variance of log(X)
! Then
! mean of X     = exp(m + 0.5s^2)
! variance of X = (mean(X))^2.[exp(s^2) - 1]

! In the reverse direction (rarely used)
! variance of log(X) = log[1 + var(X)/(mean(X))^2]
! mean of log(X)     = log(mean(X) - 0.5var(log(X))

! N.B. The above formulae relate to population parameters; they will only be
!      approximate if applied to sample values.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Version 1.13, 2 October 2000
! Changes from version 1.01
! 1. The random_order, random_Poisson & random_binomial routines have been
!    replaced with more efficient routines.
! 2. A routine, seed_random_number, has been added to seed the uniform random
!    number generator.   This requires input of the required number of seeds
!    for the particular compiler from a specified I/O unit such as a keyboard.
! 3. Made compatible with Lahey's ELF90.
! 4. Marsaglia & Tsang algorithm used for random_gamma when shape parameter > 1.
! 5. INTENT for array f corrected in random_mvnorm.

!     Author: Alan Miller
!             CSIRO Division of Mathematical & Information Sciences
!             Private Bag 10, Clayton South MDC
!             Clayton 3169, Victoria, Australia
!     Phone: (+61) 3 9545-8016      Fax: (+61) 3 9545-8080
!     e-mail: amiller @ bigpond.net.au

IMPLICIT NONE
REAL, PRIVATE      :: zero = 0.0, half = 0.5, one = 1.0, two = 2.0,   &
                      vsmall = TINY(1.0), vlarge = HUGE(1.0)
PRIVATE            :: integral
INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(12, 60)


CONTAINS


FUNCTION random_normal() RESULT(fn_val)

! Adapted from the following Fortran 77 code
!      ALGORITHM 712, COLLECTED ALGORITHMS FROM ACM.
!      THIS WORK PUBLISHED IN TRANSACTIONS ON MATHEMATICAL SOFTWARE,
!      VOL. 18, NO. 4, DECEMBER, 1992, PP. 434-435.

!  The function random_normal() returns a normally distributed pseudo-random
!  number with zero mean and unit variance.

!  The algorithm uses the ratio of uniforms method of A.J. Kinderman
!  and J.F. Monahan augmented with quadratic bounding curves.

REAL :: fn_val

!     Local variables
REAL     :: s = 0.449871, t = -0.386595, a = 0.19600, b = 0.25472,           &
            r1 = 0.27597, r2 = 0.27846, u, v, x, y, q

!     Generate P = (u,v) uniform in rectangle enclosing acceptance region

DO
  CALL RANDOM_NUMBER(u)
  CALL RANDOM_NUMBER(v)
  v = 1.7156 * (v - half)

!     Evaluate the quadratic form
  x = u - s
  y = ABS(v) - t
  q = x**2 + y*(a*y - b*x)

!     Accept P if inside inner ellipse
  IF (q < r1) EXIT
!     Reject P if outside outer ellipse
  IF (q > r2) CYCLE
!     Reject P if outside acceptance region
  IF (v**2 < -4.0*LOG(u)*u**2) EXIT
END DO

!     Return ratio of P's coordinates as the normal deviate
fn_val = v/u
RETURN

END FUNCTION random_normal



FUNCTION random_gamma(s, first) RESULT(fn_val)

! Adapted from Fortran 77 code from the book:
!     Dagpunar, J. 'Principles of random variate generation'
!     Clarendon Press, Oxford, 1988.   ISBN 0-19-852202-9

!     FUNCTION GENERATES A RANDOM GAMMA VARIATE.
!     CALLS EITHER random_gamma1 (S > 1.0)
!     OR random_exponential (S = 1.0)
!     OR random_gamma2 (S < 1.0).

!     S = SHAPE PARAMETER OF DISTRIBUTION (0 < REAL).

REAL, INTENT(IN)    :: s
LOGICAL, INTENT(IN) :: first
REAL                :: fn_val

IF (s <= zero) THEN
  WRITE(*, *) 'SHAPE PARAMETER VALUE MUST BE POSITIVE'
  STOP
END IF

IF (s > one) THEN
  fn_val = random_gamma1(s, first)
ELSE IF (s < one) THEN
  fn_val = random_gamma2(s, first)
ELSE
  fn_val = random_exponential()
END IF

RETURN
END FUNCTION random_gamma



FUNCTION random_gamma1(s, first) RESULT(fn_val)

! Uses the algorithm in
! Marsaglia, G. and Tsang, W.W. (2000) `A simple method for generating
! gamma variables', Trans. om Math. Software (TOMS), vol.26(3), pp.363-372.

! Generates a random gamma deviate for shape parameter s >= 1.

REAL, INTENT(IN)    :: s
LOGICAL, INTENT(IN) :: first
REAL                :: fn_val

! Local variables
REAL, SAVE  :: c, d
REAL        :: u, v, x

IF (first) THEN
  d = s - one/3.
  c = one/SQRT(9.0*d)
END IF

! Start of main loop
DO

! Generate v = (1+cx)^3 where x is random normal; repeat if v <= 0.

  DO
    x = random_normal()
    v = (one + c*x)**3
    IF (v > zero) EXIT
  END DO

! Generate uniform variable U

  CALL RANDOM_NUMBER(u)
  IF (u < one - 0.0331*x**4) THEN
    fn_val = d*v
    EXIT
  ELSE IF (LOG(u) < half*x**2 + d*(one - v + LOG(v))) THEN
    fn_val = d*v
    EXIT
  END IF
END DO

RETURN
END FUNCTION random_gamma1



FUNCTION random_gamma2(s, first) RESULT(fn_val)

! Adapted from Fortran 77 code from the book:
!     Dagpunar, J. 'Principles of random variate generation'
!     Clarendon Press, Oxford, 1988.   ISBN 0-19-852202-9

! FUNCTION GENERATES A RANDOM VARIATE IN [0,INFINITY) FROM
! A GAMMA DISTRIBUTION WITH DENSITY PROPORTIONAL TO
! GAMMA2**(S-1) * EXP(-GAMMA2),
! USING A SWITCHING METHOD.

!    S = SHAPE PARAMETER OF DISTRIBUTION
!          (REAL < 1.0)

REAL, INTENT(IN)    :: s
LOGICAL, INTENT(IN) :: first
REAL                :: fn_val

!     Local variables
REAL       :: r, x, w
REAL, SAVE :: a, p, c, uf, vr, d

IF (s <= zero .OR. s >= one) THEN
  WRITE(*, *) 'SHAPE PARAMETER VALUE OUTSIDE PERMITTED RANGE'
  STOP
END IF

IF (first) THEN                        ! Initialization, if necessary
  a = one - s
  p = a/(a + s*EXP(-a))
  IF (s < vsmall) THEN
    WRITE(*, *) 'SHAPE PARAMETER VALUE TOO SMALL'
    STOP
  END IF
  c = one/s
  uf = p*(vsmall/a)**s
  vr = one - vsmall
  d = a*LOG(a)
END IF

DO
  CALL RANDOM_NUMBER(r)
  IF (r >= vr) THEN
    CYCLE
  ELSE IF (r > p) THEN
    x = a - LOG((one - r)/(one - p))
    w = a*LOG(x)-d
  ELSE IF (r > uf) THEN
    x = a*(r/p)**c
    w = x
  ELSE
    fn_val = zero
    RETURN
  END IF

  CALL RANDOM_NUMBER(r)
  IF (one-r <= w .AND. r > zero) THEN
    IF (r*(w + one) >= one) CYCLE
    IF (-LOG(r) <= w) CYCLE
  END IF
  EXIT
END DO

fn_val = x
RETURN

END FUNCTION random_gamma2



FUNCTION random_chisq(ndf, first) RESULT(fn_val)

!     Generates a random variate from the chi-squared distribution with
!     ndf degrees of freedom

INTEGER, INTENT(IN) :: ndf
LOGICAL, INTENT(IN) :: first
REAL                :: fn_val

fn_val = two * random_gamma(half*ndf, first)
RETURN

END FUNCTION random_chisq



FUNCTION random_exponential() RESULT(fn_val)

! Adapted from Fortran 77 code from the book:
!     Dagpunar, J. 'Principles of random variate generation'
!     Clarendon Press, Oxford, 1988.   ISBN 0-19-852202-9

! FUNCTION GENERATES A RANDOM VARIATE IN [0,INFINITY) FROM
! A NEGATIVE EXPONENTIAL DlSTRIBUTION WlTH DENSITY PROPORTIONAL
! TO EXP(-random_exponential), USING INVERSION.

REAL  :: fn_val

!     Local variable
REAL  :: r

DO
  CALL RANDOM_NUMBER(r)
  IF (r > zero) EXIT
END DO

fn_val = -LOG(r)
RETURN

END FUNCTION random_exponential



FUNCTION random_Weibull(a) RESULT(fn_val)

!     Generates a random variate from the Weibull distribution with
!     probability density:
!                      a
!               a-1  -x
!     f(x) = a.x    e

REAL, INTENT(IN) :: a
REAL             :: fn_val

!     For speed, there is no checking that a is not zero or very small.

fn_val = random_exponential() ** (one/a)
RETURN

END FUNCTION random_Weibull



FUNCTION random_beta(aa, bb, first) RESULT(fn_val)

! Adapted from Fortran 77 code from the book:
!     Dagpunar, J. 'Principles of random variate generation'
!     Clarendon Press, Oxford, 1988.   ISBN 0-19-852202-9

! FUNCTION GENERATES A RANDOM VARIATE IN [0,1]
! FROM A BETA DISTRIBUTION WITH DENSITY
! PROPORTIONAL TO BETA**(AA-1) * (1-BETA)**(BB-1).
! USING CHENG'S LOG LOGISTIC METHOD.

!     AA = SHAPE PARAMETER FROM DISTRIBUTION (0 < REAL)
!     BB = SHAPE PARAMETER FROM DISTRIBUTION (0 < REAL)

REAL, INTENT(IN)    :: aa, bb
LOGICAL, INTENT(IN) :: first
REAL                :: fn_val

!     Local variables
REAL, PARAMETER  :: aln4 = 1.3862944
REAL             :: a, b, g, r, s, x, y, z
REAL, SAVE       :: d, f, h, t, c
LOGICAL, SAVE    :: swap

IF (aa <= zero .OR. bb <= zero) THEN
  WRITE(*, *) 'IMPERMISSIBLE SHAPE PARAMETER VALUE(S)'
  STOP
END IF

IF (first) THEN                        ! Initialization, if necessary
  a = aa
  b = bb
  swap = b > a
  IF (swap) THEN
    g = b
    b = a
    a = g
  END IF
  d = a/b
  f = a+b
  IF (b > one) THEN
    h = SQRT((two*a*b - f)/(f - two))
    t = one
  ELSE
    h = b
    t = one/(one + (a/(vlarge*b))**b)
  END IF
  c = a+h
END IF

DO
  CALL RANDOM_NUMBER(r)
  CALL RANDOM_NUMBER(x)
  s = r*r*x
  IF (r < vsmall .OR. s <= zero) CYCLE
  IF (r < t) THEN
    x = LOG(r/(one - r))/h
    y = d*EXP(x)
    z = c*x + f*LOG((one + d)/(one + y)) - aln4
    IF (s - one > z) THEN
      IF (s - s*z > one) CYCLE
      IF (LOG(s) > z) CYCLE
    END IF
    fn_val = y/(one + y)
  ELSE
    IF (4.0*s > (one + one/d)**f) CYCLE
    fn_val = one
  END IF
  EXIT
END DO

IF (swap) fn_val = one - fn_val
RETURN
END FUNCTION random_beta



FUNCTION random_t(m) RESULT(fn_val)

! Adapted from Fortran 77 code from the book:
!     Dagpunar, J. 'Principles of random variate generation'
!     Clarendon Press, Oxford, 1988.   ISBN 0-19-852202-9

! FUNCTION GENERATES A RANDOM VARIATE FROM A
! T DISTRIBUTION USING KINDERMAN AND MONAHAN'S RATIO METHOD.

!     M = DEGREES OF FREEDOM OF DISTRIBUTION
!           (1 <= 1NTEGER)

INTEGER, INTENT(IN) :: m
REAL                :: fn_val

!     Local variables
REAL, SAVE      :: s, c, a, f, g
REAL            :: r, x, v

REAL, PARAMETER :: three = 3.0, four = 4.0, quart = 0.25,   &
                   five = 5.0, sixteen = 16.0
INTEGER         :: mm = 0

IF (m < 1) THEN
  WRITE(*, *) 'IMPERMISSIBLE DEGREES OF FREEDOM'
  STOP
END IF

IF (m /= mm) THEN                    ! Initialization, if necessary
  s = m
  c = -quart*(s + one)
  a = four/(one + one/s)**c
  f = sixteen/a
  IF (m > 1) THEN
    g = s - one
    g = ((s + one)/g)**c*SQRT((s+s)/g)
  ELSE
    g = one
  END IF
  mm = m
END IF

DO
  CALL RANDOM_NUMBER(r)
  IF (r <= zero) CYCLE
  CALL RANDOM_NUMBER(v)
  x = (two*v - one)*g/r
  v = x*x
  IF (v > five - a*r) THEN
    IF (m >= 1 .AND. r*(v + three) > f) CYCLE
    IF (r > (one + v/s)**c) CYCLE
  END IF
  EXIT
END DO

fn_val = x
RETURN
END FUNCTION random_t



SUBROUTINE random_mvnorm(n, h, d, f, first, x, ier)

! Adapted from Fortran 77 code from the book:
!     Dagpunar, J. 'Principles of random variate generation'
!     Clarendon Press, Oxford, 1988.   ISBN 0-19-852202-9

! N.B. An extra argument, ier, has been added to Dagpunar's routine

!     SUBROUTINE GENERATES AN N VARIATE RANDOM NORMAL
!     VECTOR USING A CHOLESKY DECOMPOSITION.

! ARGUMENTS:
!        N = NUMBER OF VARIATES IN VECTOR
!           (INPUT,INTEGER >= 1)
!     H(J) = J'TH ELEMENT OF VECTOR OF MEANS
!           (INPUT,REAL)
!     X(J) = J'TH ELEMENT OF DELIVERED VECTOR
!           (OUTPUT,REAL)
!
!    D(J*(J-1)/2+I) = (I,J)'TH ELEMENT OF VARIANCE MATRIX (J> = I)
!            (INPUT,REAL)
!    F((J-1)*(2*N-J)/2+I) = (I,J)'TH ELEMENT OF LOWER TRIANGULAR
!           DECOMPOSITION OF VARIANCE MATRIX (J <= I)
!            (OUTPUT,REAL)

!    FIRST = .TRUE. IF THIS IS THE FIRST CALL OF THE ROUTINE
!    OR IF THE DISTRIBUTION HAS CHANGED SINCE THE LAST CALL OF THE ROUTINE.
!    OTHERWISE SET TO .FALSE.
!            (INPUT,LOGICAL)

!    ier = 1 if the input covariance matrix is not +ve definite
!        = 0 otherwise

INTEGER, INTENT(IN)   :: n
REAL(8), INTENT(IN)      :: h(:), d(:)   ! d(n*(n+1)/2)
REAL(8), INTENT(IN OUT)  :: f(:)         ! f(n*(n+1)/2)
REAL(8), INTENT(OUT)     :: x(:)
LOGICAL, INTENT(IN)   :: first
INTEGER, INTENT(OUT)  :: ier

!     Local variables
INTEGER       :: j, i, m
REAL(8)          :: y, v
INTEGER, SAVE :: n2

IF (n < 1) THEN
  WRITE(*, *) 'SIZE OF VECTOR IS NON POSITIVE'
  STOP
END IF

ier = 0
IF (first) THEN                        ! Initialization, if necessary
  n2 = 2*n
  IF (d(1) < zero) THEN
    ier = 1
    RETURN
  END IF

  f(1) = SQRT(d(1))
  y = one/f(1)
  DO j = 2,n
    f(j) = d(1+j*(j-1)/2) * y
  END DO

  DO i = 2,n
    v = d(i*(i-1)/2+i)
    DO m = 1,i-1
      v = v - f((m-1)*(n2-m)/2+i)**2
    END DO

    IF (v < zero) THEN
      ier = 1
      RETURN
    END IF

    v = SQRT(v)
    y = one/v
    f((i-1)*(n2-i)/2+i) = v
    DO j = i+1,n
      v = d(j*(j-1)/2+i)
      DO m = 1,i-1
        v = v - f((m-1)*(n2-m)/2+i)*f((m-1)*(n2-m)/2 + j)
      END DO ! m = 1,i-1
      f((i-1)*(n2-i)/2 + j) = v*y
    END DO ! j = i+1,n
  END DO ! i = 2,n
END IF

x(1:n) = h(1:n)
DO j = 1,n
  y = random_normal()
  DO i = j,n
    x(i) = x(i) + f((j-1)*(n2-j)/2 + i) * y
  END DO ! i = j,n
END DO ! j = 1,n

RETURN
END SUBROUTINE random_mvnorm



FUNCTION random_inv_gauss(h, b, first) RESULT(fn_val)

! Adapted from Fortran 77 code from the book:
!     Dagpunar, J. 'Principles of random variate generation'
!     Clarendon Press, Oxford, 1988.   ISBN 0-19-852202-9

! FUNCTION GENERATES A RANDOM VARIATE IN [0,INFINITY] FROM
! A REPARAMETERISED GENERALISED INVERSE GAUSSIAN (GIG) DISTRIBUTION
! WITH DENSITY PROPORTIONAL TO  GIG**(H-1) * EXP(-0.5*B*(GIG+1/GIG))
! USING A RATIO METHOD.

!     H = PARAMETER OF DISTRIBUTION (0 <= REAL)
!     B = PARAMETER OF DISTRIBUTION (0 < REAL)

REAL, INTENT(IN)    :: h, b
LOGICAL, INTENT(IN) :: first
REAL                :: fn_val

!     Local variables
REAL            :: ym, xm, r, w, r1, r2, x
REAL, SAVE      :: a, c, d, e
REAL, PARAMETER :: quart = 0.25

IF (h < zero .OR. b <= zero) THEN
  WRITE(*, *) 'IMPERMISSIBLE DISTRIBUTION PARAMETER VALUES'
  STOP
END IF

IF (first) THEN                        ! Initialization, if necessary
  IF (h > quart*b*SQRT(vlarge)) THEN
    WRITE(*, *) 'THE RATIO H:B IS TOO SMALL'
    STOP
  END IF
  e = b*b
  d = h + one
  ym = (-d + SQRT(d*d + e))/b
  IF (ym < vsmall) THEN
    WRITE(*, *) 'THE VALUE OF B IS TOO SMALL'
    STOP
  END IF

  d = h - one
  xm = (d + SQRT(d*d + e))/b
  d = half*d
  e = -quart*b
  r = xm + one/xm
  w = xm*ym
  a = w**(-half*h) * SQRT(xm/ym) * EXP(-e*(r - ym - one/ym))
  IF (a < vsmall) THEN
    WRITE(*, *) 'THE VALUE OF H IS TOO LARGE'
    STOP
  END IF
  c = -d*LOG(xm) - e*r
END IF

DO
  CALL RANDOM_NUMBER(r1)
  IF (r1 <= zero) CYCLE
  CALL RANDOM_NUMBER(r2)
  x = a*r2/r1
  IF (x <= zero) CYCLE
  IF (LOG(r1) < d*LOG(x) + e*(x + one/x) + c) EXIT
END DO

fn_val = x

RETURN
END FUNCTION random_inv_gauss



FUNCTION random_Poisson(mu, first) RESULT(ival)
!**********************************************************************
!     Translated to Fortran 90 by Alan Miller from:
!                           RANLIB
!
!     Library of Fortran Routines for Random Number Generation
!
!                    Compiled and Written by:
!
!                         Barry W. Brown
!                          James Lovato
!
!             Department of Biomathematics, Box 237
!             The University of Texas, M.D. Anderson Cancer Center
!             1515 Holcombe Boulevard
!             Houston, TX      77030
!
! This work was supported by grant CA-16672 from the National Cancer Institute.

!                    GENerate POIsson random deviate

!                            Function

! Generates a single random deviate from a Poisson distribution with mean mu.

!                            Arguments

!     mu --> The mean of the Poisson distribution from which
!            a random deviate is to be generated.
!                              REAL mu

!                              Method

!     For details see:

!               Ahrens, J.H. and Dieter, U.
!               Computer Generation of Poisson Deviates
!               From Modified Normal Distributions.
!               ACM Trans. Math. Software, 8, 2
!               (June 1982),163-179

!     TABLES: COEFFICIENTS A0-A7 FOR STEP F. FACTORIALS FACT
!     COEFFICIENTS A(K) - FOR PX = FK*V*V*SUM(A(K)*V**K)-DEL

!     SEPARATION OF CASES A AND B

!     .. Scalar Arguments ..
REAL, INTENT(IN)    :: mu
LOGICAL, INTENT(IN) :: first
INTEGER             :: ival
!     ..
!     .. Local Scalars ..
REAL          :: b1, b2, c, c0, c1, c2, c3, del, difmuk, e, fk, fx, fy, g,  &
                 omega, px, py, t, u, v, x, xx
REAL, SAVE    :: s, d, p, q, p0
INTEGER       :: j, k, kflag
LOGICAL, SAVE :: full_init
INTEGER, SAVE :: l, m
!     ..
!     .. Local Arrays ..
REAL, SAVE    :: pp(35)
!     ..
!     .. Data statements ..
REAL, PARAMETER :: a0 = -.5, a1 = .3333333, a2 = -.2500068, a3 = .2000118,  &
                   a4 = -.1661269, a5 = .1421878, a6 = -.1384794,   &
                   a7 = .1250060

REAL, PARAMETER :: fact(10) = (/ 1., 1., 2., 6., 24., 120., 720., 5040.,  &
                                 40320., 362880. /)

!     ..
!     .. Executable Statements ..
IF (mu > 10.0) THEN
!     C A S E  A. (RECALCULATION OF S, D, L IF MU HAS CHANGED)

  IF (first) THEN
    s = SQRT(mu)
    d = 6.0*mu*mu

!             THE POISSON PROBABILITIES PK EXCEED THE DISCRETE NORMAL
!             PROBABILITIES FK WHENEVER K >= M(MU). L=IFIX(MU-1.1484)
!             IS AN UPPER BOUND TO M(MU) FOR ALL MU >= 10 .

    l = mu - 1.1484
    full_init = .false.
  END IF


!     STEP N. NORMAL SAMPLE - random_normal() FOR STANDARD NORMAL DEVIATE

  g = mu + s*random_normal()
  IF (g > 0.0) THEN
    ival = g

!     STEP I. IMMEDIATE ACCEPTANCE IF ival IS LARGE ENOUGH

    IF (ival>=l) RETURN

!     STEP S. SQUEEZE ACCEPTANCE - SAMPLE U

    fk = ival
    difmuk = mu - fk
    CALL RANDOM_NUMBER(u)
    IF (d*u >= difmuk*difmuk*difmuk) RETURN
  END IF

!     STEP P. PREPARATIONS FOR STEPS Q AND H.
!             (RECALCULATIONS OF PARAMETERS IF NECESSARY)
!             .3989423=(2*PI)**(-.5)  .416667E-1=1./24.  .1428571=1./7.
!             THE QUANTITIES B1, B2, C3, C2, C1, C0 ARE FOR THE HERMITE
!             APPROXIMATIONS TO THE DISCRETE NORMAL PROBABILITIES FK.
!             C=.1069/MU GUARANTEES MAJORIZATION BY THE 'HAT'-FUNCTION.

  IF (.NOT. full_init) THEN
    omega = .3989423/s
    b1 = .4166667E-1/mu
    b2 = .3*b1*b1
    c3 = .1428571*b1*b2
    c2 = b2 - 15.*c3
    c1 = b1 - 6.*b2 + 45.*c3
    c0 = 1. - b1 + 3.*b2 - 15.*c3
    c = .1069/mu
    full_init = .true.
  END IF

  IF (g < 0.0) GO TO 50

!             'SUBROUTINE' F IS CALLED (KFLAG=0 FOR CORRECT RETURN)

  kflag = 0
  GO TO 70

!     STEP Q. QUOTIENT ACCEPTANCE (RARE CASE)

  40 IF (fy-u*fy <= py*EXP(px-fx)) RETURN

!     STEP E. EXPONENTIAL SAMPLE - random_exponential() FOR STANDARD EXPONENTIAL
!             DEVIATE E AND SAMPLE T FROM THE LAPLACE 'HAT'
!             (IF T <= -.6744 THEN PK < FK FOR ALL MU >= 10.)

  50 e = random_exponential()
  CALL RANDOM_NUMBER(u)
  u = u + u - one
  t = 1.8 + SIGN(e, u)
  IF (t <= (-.6744)) GO TO 50
  ival = mu + s*t
  fk = ival
  difmuk = mu - fk

!             'SUBROUTINE' F IS CALLED (KFLAG=1 FOR CORRECT RETURN)

  kflag = 1
  GO TO 70

!     STEP H. HAT ACCEPTANCE (E IS REPEATED ON REJECTION)

  60 IF (c*ABS(u) > py*EXP(px+e) - fy*EXP(fx+e)) GO TO 50
  RETURN

!     STEP F. 'SUBROUTINE' F. CALCULATION OF PX, PY, FX, FY.
!             CASE ival < 10 USES FACTORIALS FROM TABLE FACT

  70 IF (ival>=10) GO TO 80
  px = -mu
  py = mu**ival/fact(ival+1)
  GO TO 110

!             CASE ival >= 10 USES POLYNOMIAL APPROXIMATION
!             A0-A7 FOR ACCURACY WHEN ADVISABLE
!             .8333333E-1=1./12.  .3989423=(2*PI)**(-.5)

  80 del = .8333333E-1/fk
  del = del - 4.8*del*del*del
  v = difmuk/fk
  IF (ABS(v)>0.25) THEN
    px = fk*LOG(one + v) - difmuk - del
  ELSE
    px = fk*v*v* (((((((a7*v+a6)*v+a5)*v+a4)*v+a3)*v+a2)*v+a1)*v+a0) - del
  END IF
  py = .3989423/SQRT(fk)
  110 x = (half - difmuk)/s
  xx = x*x
  fx = -half*xx
  fy = omega* (((c3*xx + c2)*xx + c1)*xx + c0)
  IF (kflag <= 0) GO TO 40
  GO TO 60

!---------------------------------------------------------------------------
!     C A S E  B.    mu < 10
!     START NEW TABLE AND CALCULATE P0 IF NECESSARY

ELSE
  IF (first) THEN
    m = MAX(1, INT(mu))
    l = 0
    p = EXP(-mu)
    q = p
    p0 = p
  END IF

!     STEP U. UNIFORM SAMPLE FOR INVERSION METHOD

  DO
    CALL RANDOM_NUMBER(u)
    ival = 0
    IF (u <= p0) RETURN

!     STEP T. TABLE COMPARISON UNTIL THE END PP(L) OF THE
!             PP-TABLE OF CUMULATIVE POISSON PROBABILITIES
!             (0.458=PP(9) FOR MU=10)

    IF (l == 0) GO TO 150
    j = 1
    IF (u > 0.458) j = MIN(l, m)
    DO k = j, l
      IF (u <= pp(k)) GO TO 180
    END DO
    IF (l == 35) CYCLE

!     STEP C. CREATION OF NEW POISSON PROBABILITIES P
!             AND THEIR CUMULATIVES Q=PP(K)

    150 l = l + 1
    DO k = l, 35
      p = p*mu / k
      q = q + p
      pp(k) = q
      IF (u <= q) GO TO 170
    END DO
    l = 35
  END DO

  170 l = k
  180 ival = k
  RETURN
END IF

RETURN
END FUNCTION random_Poisson



FUNCTION random_binomial1(n, p, first) RESULT(ival)

! FUNCTION GENERATES A RANDOM BINOMIAL VARIATE USING C.D.Kemp's method.
! This algorithm is suitable when many random variates are required
! with the SAME parameter values for n & p.

!    P = BERNOULLI SUCCESS PROBABILITY
!           (0 <= REAL <= 1)
!    N = NUMBER OF BERNOULLI TRIALS
!           (1 <= INTEGER)
!    FIRST = .TRUE. for the first call using the current parameter values
!          = .FALSE. if the values of (n,p) are unchanged from last call

! Reference: Kemp, C.D. (1986). `A modal method for generating binomial
!            variables', Commun. Statist. - Theor. Meth. 15(3), 805-813.

INTEGER, INTENT(IN) :: n
REAL, INTENT(IN)    :: p
LOGICAL, INTENT(IN) :: first
INTEGER             :: ival

!     Local variables

INTEGER         :: ru, rd
INTEGER, SAVE   :: r0
REAL            :: u, pd, pu
REAL, SAVE      :: odds_ratio, p_r
REAL, PARAMETER :: zero = 0.0, one = 1.0

IF (first) THEN
  r0 = (n+1)*p
  p_r = bin_prob(n, p, r0)
  odds_ratio = p / (one - p)
END IF

CALL RANDOM_NUMBER(u)
u = u - p_r
IF (u < zero) THEN
  ival = r0
  RETURN
END IF

pu = p_r
ru = r0
pd = p_r
rd = r0
DO
  rd = rd - 1
  IF (rd >= 0) THEN
    pd = pd * (rd+1) / (odds_ratio * (n-rd))
    u = u - pd
    IF (u < zero) THEN
      ival = rd
      RETURN
    END IF
  END IF

  ru = ru + 1
  IF (ru <= n) THEN
    pu = pu * (n-ru+1) * odds_ratio / ru
    u = u - pu
    IF (u < zero) THEN
      ival = ru
      RETURN
    END IF
  END IF
END DO

!     This point should not be reached, but just in case:

ival = r0
RETURN

END FUNCTION random_binomial1



FUNCTION bin_prob(n, p, r) RESULT(fn_val)
!     Calculate a binomial probability

INTEGER, INTENT(IN) :: n, r
REAL, INTENT(IN)    :: p
REAL                :: fn_val

!     Local variable
REAL                :: one = 1.0

fn_val = EXP( lngamma(DBLE(n+1)) - lngamma(DBLE(r+1)) - lngamma(DBLE(n-r+1)) &
              + r*LOG(p) + (n-r)*LOG(one - p) )
RETURN

END FUNCTION bin_prob



FUNCTION lngamma(x) RESULT(fn_val)

! Logarithm to base e of the gamma function.
!
! Accurate to about 1.e-14.
! Programmer: Alan Miller

! Latest revision of Fortran 77 version - 28 February 1988

REAL (dp), INTENT(IN) :: x
REAL (dp)             :: fn_val

!       Local variables

REAL (dp) :: a1 = -4.166666666554424D-02, a2 = 2.430554511376954D-03,  &
             a3 = -7.685928044064347D-04, a4 = 5.660478426014386D-04,  &
             temp, arg, product, lnrt2pi = 9.189385332046727D-1,       &
             pi = 3.141592653589793D0
LOGICAL   :: reflect

!       lngamma is not defined if x = 0 or a negative integer.

IF (x > 0.d0) GO TO 10
IF (x /= INT(x)) GO TO 10
fn_val = 0.d0
RETURN

!       If x < 0, use the reflection formula:
!               gamma(x) * gamma(1-x) = pi * cosec(pi.x)

10 reflect = (x < 0.d0)
IF (reflect) THEN
  arg = 1.d0 - x
ELSE
  arg = x
END IF

!       Increase the argument, if necessary, to make it > 10.

product = 1.d0
20 IF (arg <= 10.d0) THEN
  product = product * arg
  arg = arg + 1.d0
  GO TO 20
END IF

!  Use a polynomial approximation to Stirling's formula.
!  N.B. The real Stirling's formula is used here, not the simpler, but less
!       accurate formula given by De Moivre in a letter to Stirling, which
!       is the one usually quoted.

arg = arg - 0.5D0
temp = 1.d0/arg**2
fn_val = lnrt2pi + arg * (LOG(arg) - 1.d0 + &
                  (((a4*temp + a3)*temp + a2)*temp + a1)*temp) - LOG(product)
IF (reflect) THEN
  temp = SIN(pi * x)
  fn_val = LOG(pi/temp) - fn_val
END IF
RETURN
END FUNCTION lngamma



FUNCTION random_binomial2(n, pp, first) RESULT(ival)
!**********************************************************************
!     Translated to Fortran 90 by Alan Miller from:
!                              RANLIB
!
!     Library of Fortran Routines for Random Number Generation
!
!                      Compiled and Written by:
!
!                           Barry W. Brown
!                            James Lovato
!
!               Department of Biomathematics, Box 237
!               The University of Texas, M.D. Anderson Cancer Center
!               1515 Holcombe Boulevard
!               Houston, TX      77030
!
! This work was supported by grant CA-16672 from the National Cancer Institute.

!                    GENerate BINomial random deviate

!                              Function

!     Generates a single random deviate from a binomial
!     distribution whose number of trials is N and whose
!     probability of an event in each trial is P.

!                              Arguments

!     N  --> The number of trials in the binomial distribution
!            from which a random deviate is to be generated.
!                              INTEGER N

!     P  --> The probability of an event in each trial of the
!            binomial distribution from which a random deviate
!            is to be generated.
!                              REAL P

!     FIRST --> Set FIRST = .TRUE. for the first call to perform initialization
!               the set FIRST = .FALSE. for further calls using the same pair
!               of parameter values (N, P).
!                              LOGICAL FIRST

!     random_binomial2 <-- A random deviate yielding the number of events
!                from N independent trials, each of which has
!                a probability of event P.
!                              INTEGER random_binomial

!                              Method

!     This is algorithm BTPE from:

!         Kachitvichyanukul, V. and Schmeiser, B. W.
!         Binomial Random Variate Generation.
!         Communications of the ACM, 31, 2 (February, 1988) 216.

!**********************************************************************

!*****DETERMINE APPROPRIATE ALGORITHM AND WHETHER SETUP IS NECESSARY

!     ..
!     .. Scalar Arguments ..
REAL, INTENT(IN)    :: pp
INTEGER, INTENT(IN) :: n
LOGICAL, INTENT(IN) :: first
INTEGER             :: ival
!     ..
!     .. Local Scalars ..
REAL            :: alv, amaxp, f, f1, f2, u, v, w, w2, x, x1, x2, ynorm, z, z2
REAL, PARAMETER :: zero = 0.0, half = 0.5, one = 1.0
INTEGER         :: i, ix, ix1, k, mp
INTEGER, SAVE   :: m
REAL, SAVE      :: p, q, xnp, ffm, fm, xnpq, p1, xm, xl, xr, c, al, xll,  &
                   xlr, p2, p3, p4, qn, r, g

!     ..
!     .. Executable Statements ..

!*****SETUP, PERFORM ONLY WHEN PARAMETERS CHANGE

IF (first) THEN
  p = MIN(pp, one-pp)
  q = one - p
  xnp = n * p
END IF

IF (xnp > 30.) THEN
  IF (first) THEN
    ffm = xnp + p
    m = ffm
    fm = m
    xnpq = xnp * q
    p1 = INT(2.195*SQRT(xnpq) - 4.6*q) + half
    xm = fm + half
    xl = xm - p1
    xr = xm + p1
    c = 0.134 + 20.5 / (15.3 + fm)
    al = (ffm-xl) / (ffm - xl*p)
    xll = al * (one + half*al)
    al = (xr - ffm) / (xr*q)
    xlr = al * (one + half*al)
    p2 = p1 * (one + c + c)
    p3 = p2 + c / xll
    p4 = p3 + c / xlr
  END IF

!*****GENERATE VARIATE, Binomial mean at least 30.

  20 CALL RANDOM_NUMBER(u)
  u = u * p4
  CALL RANDOM_NUMBER(v)

!     TRIANGULAR REGION

  IF (u <= p1) THEN
    ix = xm - p1 * v + u
    GO TO 110
  END IF

!     PARALLELOGRAM REGION

  IF (u <= p2) THEN
    x = xl + (u-p1) / c
    v = v * c + one - ABS(xm-x) / p1
    IF (v > one .OR. v <= zero) GO TO 20
    ix = x
  ELSE

!     LEFT TAIL

    IF (u <= p3) THEN
      ix = xl + LOG(v) / xll
      IF (ix < 0) GO TO 20
      v = v * (u-p2) * xll
    ELSE

!     RIGHT TAIL

      ix = xr - LOG(v) / xlr
      IF (ix > n) GO TO 20
      v = v * (u-p3) * xlr
    END IF
  END IF

!*****DETERMINE APPROPRIATE WAY TO PERFORM ACCEPT/REJECT TEST

  k = ABS(ix-m)
  IF (k <= 20 .OR. k >= xnpq/2-1) THEN

!     EXPLICIT EVALUATION

    f = one
    r = p / q
    g = (n+1) * r
    IF (m < ix) THEN
      mp = m + 1
      DO i = mp, ix
        f = f * (g/i-r)
      END DO

    ELSE IF (m > ix) THEN

      ix1 = ix + 1
      DO i = ix1, m
        f = f / (g/i-r)
      END DO
    END IF

    IF (v > f) THEN
      GO TO 20
    ELSE
      GO TO 110
    END IF
  END IF

!     SQUEEZING USING UPPER AND LOWER BOUNDS ON LOG(F(X))

  amaxp = (k/xnpq) * ((k*(k/3. + .625) + .1666666666666)/xnpq + half)
  ynorm = -k * k / (2.*xnpq)
  alv = LOG(v)
  IF (alv<ynorm - amaxp) GO TO 110
  IF (alv>ynorm + amaxp) GO TO 20

!     STIRLING'S (actually de Moivre's) FORMULA TO MACHINE ACCURACY FOR
!     THE FINAL ACCEPTANCE/REJECTION TEST

  x1 = ix + 1
  f1 = fm + one
  z = n + 1 - fm
  w = n - ix + one
  z2 = z * z
  x2 = x1 * x1
  f2 = f1 * f1
  w2 = w * w
  IF (alv - (xm*LOG(f1/x1) + (n-m+half)*LOG(z/w) + (ix-m)*LOG(w*p/(x1*q)) +    &
      (13860.-(462.-(132.-(99.-140./f2)/f2)/f2)/f2)/f1/166320. +               &
      (13860.-(462.-(132.-(99.-140./z2)/z2)/z2)/z2)/z/166320. +                &
      (13860.-(462.-(132.-(99.-140./x2)/x2)/x2)/x2)/x1/166320. +               &
      (13860.-(462.-(132.-(99.-140./w2)/w2)/w2)/w2)/w/166320.) > zero) THEN
    GO TO 20
  ELSE
    GO TO 110
  END IF

ELSE
!     INVERSE CDF LOGIC FOR MEAN LESS THAN 30
  IF (first) THEN
    qn = q ** n
    r = p / q
    g = r * (n+1)
  END IF

  90 ix = 0
  f = qn
  CALL RANDOM_NUMBER(u)
  100 IF (u >= f) THEN
    IF (ix > 110) GO TO 90
    u = u - f
    ix = ix + 1
    f = f * (g/ix - r)
    GO TO 100
  END IF
END IF

110 IF (pp > half) ix = n - ix
ival = ix
RETURN

END FUNCTION random_binomial2




FUNCTION random_neg_binomial(sk, p) RESULT(ival)

! Adapted from Fortran 77 code from the book:
!     Dagpunar, J. 'Principles of random variate generation'
!     Clarendon Press, Oxford, 1988.   ISBN 0-19-852202-9

! FUNCTION GENERATES A RANDOM NEGATIVE BINOMIAL VARIATE USING UNSTORED
! INVERSION AND/OR THE REPRODUCTIVE PROPERTY.

!    SK = NUMBER OF FAILURES REQUIRED (Dagpunar's words!)
!       = the `power' parameter of the negative binomial
!           (0 < REAL)
!    P = BERNOULLI SUCCESS PROBABILITY
!           (0 < REAL < 1)

! THE PARAMETER H IS SET SO THAT UNSTORED INVERSION ONLY IS USED WHEN P <= H,
! OTHERWISE A COMBINATION OF UNSTORED INVERSION AND
! THE REPRODUCTIVE PROPERTY IS USED.

REAL, INTENT(IN)   :: sk, p
INTEGER            :: ival

!     Local variables
! THE PARAMETER ULN = -LOG(MACHINE'S SMALLEST REAL NUMBER).

REAL, PARAMETER    :: h = 0.7
REAL               :: q, x, st, uln, v, r, s, y, g
INTEGER            :: k, i, n

IF (sk <= zero .OR. p <= zero .OR. p >= one) THEN
  WRITE(*, *) 'IMPERMISSIBLE DISTRIBUTION PARAMETER VALUES'
  STOP
END IF

q = one - p
x = zero
st = sk
IF (p > h) THEN
  v = one/LOG(p)
  k = st
  DO i = 1,k
    DO
      CALL RANDOM_NUMBER(r)
      IF (r > zero) EXIT
    END DO
    n = v*LOG(r)
    x = x + n
  END DO
  st = st - k
END IF

s = zero
uln = -LOG(vsmall)
IF (st > -uln/LOG(q)) THEN
  WRITE(*, *) ' P IS TOO LARGE FOR THIS VALUE OF SK'
  STOP
END IF

y = q**st
g = st
CALL RANDOM_NUMBER(r)
DO
  IF (y > r) EXIT
  r = r - y
  s = s + one
  y = y*p*g/s
  g = g + one
END DO

ival = x + s + half
RETURN
END FUNCTION random_neg_binomial



FUNCTION random_von_Mises(k, first) RESULT(fn_val)

!     Algorithm VMD from:
!     Dagpunar, J.S. (1990) `Sampling from the von Mises distribution via a
!     comparison of random numbers', J. of Appl. Statist., 17, 165-168.

!     Fortran 90 code by Alan Miller
!     CSIRO Division of Mathematical & Information Sciences

!     Arguments:
!     k (real)        parameter of the von Mises distribution.
!     first (logical) set to .TRUE. the first time that the function
!                     is called, or the first time with a new value
!                     for k.   When first = .TRUE., the function sets
!                     up starting values and may be very much slower.

REAL, INTENT(IN)     :: k
LOGICAL, INTENT(IN)  :: first
REAL                 :: fn_val

!     Local variables

INTEGER          :: j, n
INTEGER, SAVE    :: nk
REAL, PARAMETER  :: pi = 3.14159265
REAL, SAVE       :: p(20), theta(0:20)
REAL             :: sump, r, th, lambda, rlast
REAL (dp)        :: dk

IF (first) THEN                        ! Initialization, if necessary
  IF (k < zero) THEN
    WRITE(*, *) '** Error: argument k for random_von_Mises = ', k
    RETURN
  END IF

  nk = k + k + one
  IF (nk > 20) THEN
    WRITE(*, *) '** Error: argument k for random_von_Mises = ', k
    RETURN
  END IF

  dk = k
  theta(0) = zero
  IF (k > half) THEN

!     Set up array p of probabilities.

    sump = zero
    DO j = 1, nk
      IF (j < nk) THEN
        theta(j) = ACOS(one - j/k)
      ELSE
        theta(nk) = pi
      END IF

!     Numerical integration of e^[k.cos(x)] from theta(j-1) to theta(j)

      CALL integral(theta(j-1), theta(j), p(j), dk)
      sump = sump + p(j)
    END DO
    p(1:nk) = p(1:nk) / sump
  ELSE
    p(1) = one
    theta(1) = pi
  END IF                         ! if k > 0.5
END IF                           ! if first

CALL RANDOM_NUMBER(r)
DO j = 1, nk
  r = r - p(j)
  IF (r < zero) EXIT
END DO
r = -r/p(j)

DO
  th = theta(j-1) + r*(theta(j) - theta(j-1))
  lambda = k - j + one - k*COS(th)
  n = 1
  rlast = lambda

  DO
    CALL RANDOM_NUMBER(r)
    IF (r > rlast) EXIT
    n = n + 1
    rlast = r
  END DO

  IF (n .NE. 2*(n/2)) EXIT         ! is n even?
  CALL RANDOM_NUMBER(r)
END DO

fn_val = SIGN(th, (r - rlast)/(one - rlast) - half)
RETURN
END FUNCTION random_von_Mises



SUBROUTINE integral(a, b, result, dk)

!     Gaussian integration of exp(k.cosx) from a to b.

REAL (dp), INTENT(IN) :: dk
REAL, INTENT(IN)      :: a, b
REAL, INTENT(OUT)     :: result

!     Local variables

REAL (dp)  :: xmid, range, x1, x2,                                    &
  x(3) = (/0.238619186083197_dp, 0.661209386466265_dp, 0.932469514203152_dp/), &
  w(3) = (/0.467913934572691_dp, 0.360761573048139_dp, 0.171324492379170_dp/)
INTEGER    :: i

xmid = (a + b)/2._dp
range = (b - a)/2._dp

result = 0._dp
DO i = 1, 3
  x1 = xmid + x(i)*range
  x2 = xmid - x(i)*range
  result = result + w(i)*(EXP(dk*COS(x1)) + EXP(dk*COS(x2)))
END DO

result = result * range
RETURN
END SUBROUTINE integral



FUNCTION random_Cauchy() RESULT(fn_val)

!     Generate a random deviate from the standard Cauchy distribution

REAL     :: fn_val

!     Local variables
REAL     :: v(2)

DO
  CALL RANDOM_NUMBER(v)
  v = two*(v - half)
  IF (ABS(v(2)) < vsmall) CYCLE               ! Test for zero
  IF (v(1)**2 + v(2)**2 < one) EXIT
END DO
fn_val = v(1) / v(2)

RETURN
END FUNCTION random_Cauchy



SUBROUTINE random_order(order, n)

!     Generate a random ordering of the integers 1 ... n.

INTEGER, INTENT(IN)  :: n
INTEGER, INTENT(OUT) :: order(n)

!     Local variables

INTEGER :: i, j, k
REAL    :: wk

DO i = 1, n
  order(i) = i
END DO

!     Starting at the end, swap the current last indicator with one
!     randomly chosen from those preceeding it.

DO i = n, 2, -1
  CALL RANDOM_NUMBER(wk)
  j = 1 + i * wk
  IF (j < i) THEN
    k = order(i)
    order(i) = order(j)
    order(j) = k
  END IF
END DO

RETURN
END SUBROUTINE random_order



SUBROUTINE seed_random_number(iounit)

INTEGER, INTENT(IN)  :: iounit

! Local variables

INTEGER              :: k
INTEGER, ALLOCATABLE :: seed(:)

CALL RANDOM_SEED(SIZE=k)
ALLOCATE( seed(k) )

WRITE(*, '(a, i2, a)')' Enter ', k, ' integers for random no. seeds: '
READ(*, *) seed
WRITE(iounit, '(a, (7i10))') ' Random no. seeds: ', seed
CALL RANDOM_SEED(PUT=seed)

DEALLOCATE( seed )

RETURN
END SUBROUTINE seed_random_number


function alnorm ( x, upper )

!*****************************************************************************80
!
!! ALNORM computes the cumulative density of the standard normal distribution.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    13 January 2008
!
!  Author:
!
!    Original FORTRAN77 version by David Hill.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    David Hill,
!    Algorithm AS 66:
!    The Normal Integral,
!    Applied Statistics,
!    Volume 22, Number 3, 1973, pages 424-427.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, is one endpoint of the semi-infinite interval
!    over which the integration takes place.
!
!    Input, logical UPPER, determines whether the upper or lower
!    interval is to be integrated:
!    .TRUE.  => integrate from X to + Infinity;
!    .FALSE. => integrate from - Infinity to X.
!
!    Output, real ( kind = 8 ) ALNORM, the integral of the standard normal
!    distribution over the desired interval.
!
  implicit none

  real ( kind = 8 ), parameter :: a1 = 5.75885480458D+00
  real ( kind = 8 ), parameter :: a2 = 2.62433121679D+00
  real ( kind = 8 ), parameter :: a3 = 5.92885724438D+00
  real ( kind = 8 ) alnorm
  real ( kind = 8 ), parameter :: b1 = -29.8213557807D+00
  real ( kind = 8 ), parameter :: b2 = 48.6959930692D+00
  real ( kind = 8 ), parameter :: c1 = -0.000000038052D+00
  real ( kind = 8 ), parameter :: c2 = 0.000398064794D+00
  real ( kind = 8 ), parameter :: c3 = -0.151679116635D+00
  real ( kind = 8 ), parameter :: c4 = 4.8385912808D+00
  real ( kind = 8 ), parameter :: c5 = 0.742380924027D+00
  real ( kind = 8 ), parameter :: c6 = 3.99019417011D+00
  real ( kind = 8 ), parameter :: con = 1.28D+00
  real ( kind = 8 ), parameter :: d1 = 1.00000615302D+00
  real ( kind = 8 ), parameter :: d2 = 1.98615381364D+00
  real ( kind = 8 ), parameter :: d3 = 5.29330324926D+00
  real ( kind = 8 ), parameter :: d4 = -15.1508972451D+00
  real ( kind = 8 ), parameter :: d5 = 30.789933034D+00
  real ( kind = 8 ), parameter :: ltone = 7.0D+00
  real ( kind = 8 ), parameter :: p = 0.398942280444D+00
  real ( kind = 8 ), parameter :: q = 0.39990348504D+00
  real ( kind = 8 ), parameter :: r = 0.398942280385D+00
  logical up
  logical upper
  real ( kind = 8 ), parameter :: utzero = 18.66D+00
  real ( kind = 8 ) x
  real ( kind = 8 ) y
  real ( kind = 8 ) z

  up = upper
  z = x

  if ( z < 0.0D+00 ) then
    up = .not. up
    z = - z
  end if

  if ( ltone < z .and. ( ( .not. up ) .or. utzero < z ) ) then

    if ( up ) then
      alnorm = 0.0D+00
    else
      alnorm = 1.0D+00
    end if

    return

  end if

  y = 0.5D+00 * z * z

  if ( z <= con ) then

    alnorm = 0.5D+00 - z * ( p - q * y &
      / ( y + a1 + b1 &
      / ( y + a2 + b2 &
      / ( y + a3 ))))

  else

    alnorm = r * exp ( - y ) &
      / ( z + c1 + d1 &
      / ( z + c2 + d2 &
      / ( z + c3 + d3 &
      / ( z + c4 + d4 &
      / ( z + c5 + d5 &
      / ( z + c6 ))))))

  end if

  if ( .not. up ) then
    alnorm = 1.0D+00 - alnorm
  end if

  return
end function alnorm

END MODULE random

module pchip
! this is part of nms, numerical analysis library from : http://people.sc.fsu.edu/~jburkardt/f_src/nms/nms.html


implicit none

contains


function pchst ( arg1, arg2 )

!*****************************************************************************80
!
!! PCHST: PCHIP sign-testing routine.
!
!  Discussion:
!
!    This routine essentially computes the sign of ARG1 * ARG2.
!
!    The object is to do this without multiplying ARG1 * ARG2, to avoid
!    possible over/underflow problems.
!
!  Author:
!
!    Fred Fritsch,
!    Mathematics and Statistics Division,
!    Lawrence Livermore National Laboratory.
!
!    FORTRAN90 translation by John Burkardt.
!
!  Reference:
!
!    Fred Fritsch, Ralph Carlson,
!    Monotone Piecewise Cubic Interpolation,
!    SIAM Journal on Numerical Analysis,
!    Volume 17, Number 2, April 1980, pages 238-246.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) ARG1, ARG2, two values to check.
!
!    Output, real ( kind = 8 ) PCHST,
!    -1.0, if ARG1 and ARG2 are of opposite sign.
!     0.0, if either argument is zero.
!    +1.0, if ARG1 and ARG2 are of the same sign.
!


  real ( kind = 8 ) arg1
  real ( kind = 8 ) arg2
  real ( kind = 8 ) pchst

  pchst = sign ( 1.0D+00, arg1 ) * sign ( 1.0D+00, arg2 )

  if ( arg1 == 0.0D+00 .or. arg2 == 0.0D+00 ) then
    pchst = 0.0D+00
  end if

end function

function pchid ( n, x, f, d, incfd, skip, ia, ib, ierr )

!*****************************************************************************80
!
!! PCHID evaluates the definite integral of a piecewise cubic Hermite function.
!
!  Description:
!
!    PCHID evaluates the definite integral of a cubic Hermite function
!    over the interval [X(IA), X(IB)].  The endpoints of the integration
!    interval must be data points.
!
!    To provide compatibility with PCHIM and PCHIC, includes an
!    increment between successive values of the F and D arrays.
!
!  Author:
!
!    Fred Fritsch,
!    Mathematics and Statistics Division,
!    Lawrence Livermore National Laboratory.
!
!    FORTRAN90 translation by John Burkardt.
!
!  Reference:
!
!    Fred Fritsch, Ralph Carlson,
!    Monotone Piecewise Cubic Interpolation,
!    SIAM Journal on Numerical Analysis,
!    Volume 17, Number 2, April 1980, pages 238-246.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of data points.
!    N must be at least 2.
!
!    Input, real ( kind = 8 ) X(N), the strictly increasing independent
!    variable values.
!
!    Input, real ( kind = 8 ) F(INCFD,N), the function values.
!    F(1+(I-1)*INCFD) is the value corresponding to X(I).
!
!    Input, real ( kind = 8 ) F(INCFD,N), the derivative values.
!    D(1+(I-1)*INCFD) is the value corresponding to X(I).
!
!    Input, integer ( kind = 4 ) INCFD, increment between successive values
!    in F and D.
!
!    Input/output, logical SKIP, should be set to TRUE if the user wishes to
!    skip checks for validity of preceding parameters, or to FALSE otherwise.
!    This will save time in case these checks have already been performed
!    say, in PCHIM or PCHIC.  SKIP will be set to TRUE on return with
!    IERR = 0 or -4.
!
!    Input, integer ( kind = 4 ) IA, IB, the indices in the X array for the
!    limits of integration.  Both must be in the range [1,N].
!
!    Output, integer ( kind = 4 ) IERR, error flag.
!    0, no errors.
!    -1, if N < 2.
!    -2, if INCFD < 1.
!    -3, if the X array is not strictly increasing.
!    -4, if IA or IB is out of range.
!
!    Output, real ( kind = 8 ) PCHID, the value of the requested integral.
!


  integer ( kind = 4 ) incfd
  integer ( kind = 4 )  n

  real ( kind = 8 ) d(incfd,n)
  real ( kind = 8 ) f(incfd,n)
  real ( kind = 8 ) h
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ia
  integer ( kind = 4 ) ib
  integer ( kind = 4 ) ierr
  integer ( kind = 4 ) iup
  integer ( kind = 4 ) low
  real ( kind = 8 ) pchid
  logical skip
  real ( kind = 8 ) sum2
  real ( kind = 8 ) value
  real ( kind = 8 ) x(n)

  if ( .not. skip ) then

    if ( n < 2 ) then
      ierr = -1
      !call xerror ('pchid -- number of data points less than two', ierr, 1)
      return
    end if

    if ( incfd < 1 ) then
      ierr = -2
      !call xerror ('pchid -- increment less than one', ierr, 1)
      return
    end if

    do i = 2, n
      if ( x(i) <= x(i-1) ) then
        ierr = -3
        !call xerror ('pchid -- x-array not strictly increasing', ierr, 1)
        return
      end if
    end do

  end if

  skip = .true.

  if ( ia < 1 .or. n < ia ) then
    go to 5004
  end if

  if ( ib < 1 .or. n < ib ) then
    go to 5004
  end if

  ierr = 0
!
!  Compute integral value.
!
  if ( ia == ib ) then

    value = 0.0D+00

  else

    low = min ( ia, ib )
    iup = max ( ia, ib ) - 1
    sum2 = 0.0D+00

    do i = low, iup
      h = x(i+1) - x(i)
      sum2 = sum2 + h * &
        ( ( f(1,i) + f(1,i+1) ) + ( d(1,i) - d(1,i+1) ) * ( h / 6.0D+00 ) )
    end do

    value = 0.5D+00 * sum2

    if ( ib < ia ) then
      value = -value
    end if

  end if

  pchid = value

  return
!
!  error returns.
!
 5004 continue
!
!  ia or ib out of range return.
!
  ierr = -4
  !call xerror ('pchid -- ia or ib out of range', ierr, 1)

end function

function chfiv ( x1, x2, f1, f2, d1, d2, a, b, ierr )

!*****************************************************************************80
!
!! CHFIV evaluates the integral of a cubic polynomial in Hermite form.
!
!  Discussion:
!
!    CHFIV is called by PCHIA to evaluate the integral of a single cubic (in
!    Hermite form) over an arbitrary interval (A,B).
!
!  Author:
!
!    Fred Fritsch,
!    Mathematics and Statistics Division,
!    Lawrence Livermore National Laboratory.
!
!    FORTRAN90 translation by John Burkardt.
!
!  Reference:
!
!    Fred Fritsch, Ralph Carlson,
!    Monotone Piecewise Cubic Interpolation,
!    SIAM Journal on Numerical Analysis,
!    Volume 17, Number 2, April 1980, pages 238-246.
!
!    David Kahaner, Cleve Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1989,
!    ISBN: 0-13-627258-4,
!    LC: TA345.K34.
!
!  Parameters:
!
!    Output, real ( kind = 8 ) VALUE, the value of the requested integral.
!
!    Input, real ( kind = 8 ) X1, X2, the endpoints of the interval of
!    definition of the cubic.  X1 and X2 must be distinct.
!
!    Input, real ( kind = 8 ) F1, F2, the values of the function at X1
!    and X2, respectively.
!
!    Input, real ( kind = 8 ) D1, D2, the derivative values at the ends
!    of the interval.
!
!    Input, real ( kind = 8 ) A, B, the endpoints of interval of integration.
!
!    Output, integer ( kind = 4 ) IERR, error flag.
!    0, no errors.
!    -1, X1 == X2.
!


  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) chfiv
  real ( kind = 8 ) d1
  real ( kind = 8 ) d2
  real ( kind = 8 ) dterm
  real ( kind = 8 ) f1
  real ( kind = 8 ) f2
  real ( kind = 8 ) fterm
  real ( kind = 8 ) h
  integer ( kind = 4 ) ierr
  real ( kind = 8 ) phia1
  real ( kind = 8 ) phia2
  real ( kind = 8 ) phib1
  real ( kind = 8 ) phib2
  real ( kind = 8 ) psia1
  real ( kind = 8 ) psia2
  real ( kind = 8 ) psib1
  real ( kind = 8 ) psib2
  real ( kind = 8 ) ta1
  real ( kind = 8 ) ta2
  real ( kind = 8 ) tb1
  real ( kind = 8 ) tb2
  real ( kind = 8 ) ua1
  real ( kind = 8 ) ua2
  real ( kind = 8 ) ub1
  real ( kind = 8 ) ub2
  real ( kind = 8 ) x1
  real ( kind = 8 ) x2
!
!  Check input.
!
  if ( x1 == x2 ) then
    ierr = -1
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'CHFIV - Fatal error!'
    write ( *, '(a)' ) '  X1 = X2.'
    stop
  end if

  ierr = 0
!
!  Compute integral.
!
  h = x2 - x1
  ta1 = ( a - x1 ) / h
  ta2 = ( x2 - a ) / h
  tb1 = ( b - x1 ) / h
  tb2 = ( x2 - b ) / h

  ua1 = ta1 * ta1 * ta1
  phia1 = ua1 * ( 2.0D+00 - ta1 )
  psia1 = ua1 * ( 3.0D+00 * ta1 - 4.0D+00 )
  ua2 = ta2 * ta2 * ta2
  phia2 =  ua2 * ( 2.0D+00 - ta2)
  psia2 = -ua2 * ( 3.0D+00 * ta2 - 4.0D+00 )

  ub1 = tb1 * tb1 * tb1
  phib1 = ub1 * ( 2.0D+00 - tb1 )
  psib1 = ub1 * ( 3.0D+00 * tb1 - 4.0D+00 )
  ub2 = tb2 * tb2 * tb2
  phib2 =  ub2 * ( 2.0D+00 - tb2 )
  psib2 = -ub2 * ( 3.0D+00 * tb2 - 4.0D+00 )

  fterm =   f1 * ( phia2 - phib2 ) + f2 * ( phib1 - phia1 )
  dterm = ( d1 * ( psia2 - psib2 ) + d2 * ( psib1 - psia1 ) ) * ( h / 6.0D+00 )

  chfiv = 0.5D+00 * h * ( fterm + dterm )

end function

function chfmc ( d1, d2, delta )

!*****************************************************************************80
!
!! CHFMC determines the monotonicity properties of a cubic polynomial.
!
!  Discussion:
!
!    CHFMC is called by PCHMC to determine the monotonicity properties
!    of the cubic with boundary derivative values D1, D2 and chord
!    slope DELTA.
!
!  Author:
!
!    Fred Fritsch,
!    Mathematics and Statistics Division,
!    Lawrence Livermore National Laboratory.
!
!    FORTRAN90 translation by John Burkardt.
!
!  Reference:
!
!    Fred Fritsch, Ralph Carlson,
!    Monotone Piecewise Cubic Interpolation,
!    SIAM Journal on Numerical Analysis,
!    Volume 17, Number 2, April 1980, pages 238-246.
!
!    David Kahaner, Cleve Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1989,
!    ISBN: 0-13-627258-4,
!    LC: TA345.K34.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) D1, D2, the derivative values at the ends
!    of the interval.
!
!    Input, real ( kind = 8 ) DELTA, the data slope over that interval.
!
!    Output, integer ( kind = 4 ) CHFMC, indicates the monotonicity of the
!    cubic segment:
!    -1, if function is strictly decreasing;
!     0, if function is constant;
!     1, if function is strictly increasing;
!     2, if function is non-monotonic;
!     3, if unable to determine.
!


  real ( kind = 8 ) a
  real ( kind = 8 ) b
  integer ( kind = 4 ) chfmc
  real ( kind = 8 ) d1
  real ( kind = 8 ) d2
  real ( kind = 8 ) delta
  real ( kind = 8 ) eps
  integer ( kind = 4 ) ismon
  integer ( kind = 4 ) itrue
  real ( kind = 8 ) phi

  eps = 10.0D+00 * epsilon ( eps )
!
!  Make the check.
!
  if ( delta == 0.0D+00 ) then

    if ( d1 == 0.0D+00 .and. d2 == 0.0D+00 ) then
      ismon = 0
    else
      ismon = 2
    end if

  else

     itrue = sign ( 1.0D+00, delta)
     a = d1 / delta
     b = d2 / delta

     if ( a < 0.0D+00 .or. b < 0.0D+00 ) then
       ismon = 2
     else if ( a <= 3.0D+00 - eps  .and. b <= 3.0D+00 -eps ) then
!
!  Inside square (0,3)x(0,3) implies OK.
!
       ismon = itrue
     else if ( 4.0D+00 + eps < a .and. 4.0D+00 + eps < b ) then
!
!  Outside square (0,4)x(0,4) implies nonmonotonic.
!
       ismon = 2
     else
!
!  Must check against boundary of ellipse.
!
        a = a - 2.0D+00
        b = b - 2.0D+00
        phi = ( ( a * a + b * b ) + a * b ) - 3.0D+00

        if ( phi < -eps ) then
          ismon = itrue
        else if ( eps < phi ) then
          ismon = 2
        else
!
!  Too close to boundary to tell,
!  in the presence of round-off errors.
!
           ismon = 3
        end if
    end if
  end if

  chfmc = ismon

end function

function pchdf ( k, x, s, ierr )

!*****************************************************************************80
!
!! PCHDF approximates a derivative using divided differences of data.
!
!  Discussion:
!
!    The routine uses a divided difference formulation to compute a K-point
!    approximation to the derivative at X(K) based on the data in X and S.
!
!    It is called by PCHCE and PCHSP to compute 3 and 4 point boundary
!    derivative approximations.
!
!  Author:
!
!    Fred Fritsch,
!    Mathematics and Statistics Division,
!    Lawrence Livermore National Laboratory.
!
!    FORTRAN90 translation by John Burkardt.
!
!  Reference:
!
!    Carl deBoor,
!    A Practical Guide to Splines,
!    Pages 10-16,
!    Springer-Verlag, 1978.
!
!    Fred Fritsch, Ralph Carlson,
!    Monotone Piecewise Cubic Interpolation,
!    SIAM Journal on Numerical Analysis,
!    Volume 17, Number 2, April 1980, pages 238-246.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) K, is the order of the desired derivative
!    approximation.  K must be at least 3.
!
!    Input, real ( kind = 8 ) X(K), contains the K values of the independent
!    variable.  X need not be ordered, but the values must be distinct.
!
!    Input/output, real ( kind = 8 ) S(K-1).  On input, the associated slope
!    values:
!      S(I) = ( F(I+1)-F(I))/(X(I+1)-X(I))
!    On output, S is overwritten.
!
!    Output, integer ( kind = 4 ) IERR, error flag.
!    0, no error.
!    -1, if K < 2.
!
!    Output, real ( kind = 8 ) PCHDF, the desired derivative approximation if
!    IERR=0 or to zero if IERR=-1.
!


  integer ( kind = 4 ) k

  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierr
  integer ( kind = 4 ) j
  real ( kind = 8 ) pchdf
  real ( kind = 8 ) s(k-1)
  real ( kind = 8 ) value
  real ( kind = 8 ) x(k)
!
!  Check for legal value of K.
!
  if ( k < 3 ) then
    ierr = -1
    !call xerror ( 'pchdf -- k less than three', ierr, 1 )
    pchdf = 0.0D+00
    return
  end if
!
!  Compute coefficients of interpolating polynomial.
!
  do j = 2, k-1
    do i = 1, k-j
      s(i) = ( s(i+1) - s(i) ) / ( x(i+j) - x(i) )
    end do
  end do
!
!  Evaluate the derivative at X(K).
!
  value = s(1)

  do i = 2, k-1
    value = s(i) + value * ( x(k) - x(i) )
  end do

  ierr = 0
  pchdf = value

end function

function pchia ( n, x, f, d, incfd, skip, a, b, ierr )

!*****************************************************************************80
!
!! PCHIA evaluates the integral of a piecewise cubic Hermite function.
!
!  Description:
!
!    PCHIA evaluates the definite integral of a cubic Hermite function
!    over the interval [A, B].
!
!    To provide compatibility with PCHIM and PCHIC, includes an
!    increment between successive values of the F and D arrays.
!
!  Author:
!
!    Fred Fritsch,
!    Mathematics and Statistics Division,
!    Lawrence Livermore National Laboratory.
!
!    FORTRAN90 translation by John Burkardt.
!
!  Reference:
!
!    Fred Fritsch, Ralph Carlson,
!    Monotone Piecewise Cubic Interpolation,
!    SIAM Journal on Numerical Analysis,
!    Volume 17, Number 2, April 1980, pages 238-246.
!
!  Parameters:
!
!    Output, real ( kind = 8 ) VALUE, the value of the requested integral.
!
!    Input, integer ( kind = 4 ) N, the number of data points.
!    N must be at least 2.
!
!    Input, real ( kind = 8 ) X(N), the strictly increasing independent
!    variable values.
!
!    Input, real ( kind = 8 ) F(INCFD,N), the function values.
!    F(1+(I-1)*INCFD) is the value corresponding to X(I).
!
!    Input, real ( kind = 8 ) D(INCFD,N), the derivative values.
!    D(1+(I-1)*INCFD) is the value corresponding to X(I).
!
!    Input, integer ( kind = 4 ) INCFD, increment between successive values
!    in F and D.
!
!    Input/output, logical SKIP, controls whether data validity checks
!    should be made.  Setting the input value to FALSE will skip the checks.
!    On output with 0 <= IERR, SKIP will be set to TRUE.
!
!    Input, real ( kind = 8 ) A, B, the limits of integration.  The
!    integration interval is normally contained within [X(1),X(N)], but
!    this is not required.
!
!    Output, integer ( kind = 4 ) IERR, error flag.
!    0, no errors.
!    1, if A is outside the interval [X(1),X(N)].
!    2, if B is outside the interval [X(1),X(N)].
!    3, if both of the above are true.  This means that either [A,B] contains
!       the data interval or the intervals do not intersect at all.
!    -1, if N < 2.
!    -2, if INCFD < 1.
!    -3, if the X array is not strictly increasing.
!


  integer ( kind = 4 ) incfd
  integer ( kind = 4 ) n

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  !real ( kind = 8 ) chfiv
  real ( kind = 8 ) d(incfd,n)
  real ( kind = 8 ) f(incfd,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ia
  integer ( kind = 4 ) ib
  integer ( kind = 4 ) ierd
  integer ( kind = 4 ) ierr
  integer ( kind = 4 ) ierv
  integer ( kind = 4 ) il
  integer ( kind = 4 ) ir
  real ( kind = 8 ) pchia
  !real ( kind = 8 ) pchid
  logical skip
  real ( kind = 8 ) value
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) xa
  real ( kind = 8 ) xb
!
!  Check arguments.
!
  if ( .not. skip ) then

    if ( n < 2 ) then
      ierr = -1
      !call xerror ('pchia -- number of data points less than two', ierr, 1)
      return
    end if

    if ( incfd < 1 ) then
      ierr = -2
      !call xerror ('pchia -- increment less than one', ierr, 1)
      return
    end if

    do i = 2, n
      if ( x(i) <= x(i-1) ) then
        ierr = -3
        !call xerror ('pchia -- x-array not strictly increasing', ierr, 1)
        return
      end if
    end do

    skip = .true.

  end if

  ierr = 0

  if ( a < x(1) .or. x(n) < a ) then
    ierr = ierr + 1
  end if

  if ( b < x(1) .or. x(n) < b ) then
    ierr = ierr + 2
  end if
!
!  Compute integral value.
!
  if ( a == b ) then

    value = 0.0D+00

  else

    xa = min (a, b)
    xb = max (a, b)
!
!  Interval is to left of X(2), so use first cubic.
!
    if ( xb <= x(2) ) then

      value = chfiv ( x(1), x(2), f(1,1), f(1,2), &
        d(1,1), d(1,2), a, b, ierv )

      if ( ierv < 0 ) then
        ierr = -4
        !call xerror ('pchia -- trouble in chfiv', ierr, 1)
        return
      end if
!
!  Interval is to right of x(n-1), so use last cubic.
!
    else if ( x(n-1) <= xa ) then

      value = chfiv ( x(n-1), x(n), f(1,n-1), f(1,n), &
        d(1,n-1), d(1,n), a, b, ierv )

      if ( ierv < 0 ) then
        ierr = -4
        !call xerror ('pchia -- trouble in chfiv', ierr, 1)
        return
      end if
!
!  Normal case -- xa<xb, xa<x(n-1), x(2) < xb.
!
!  Locate ia and ib such that
!  x(ia-1) < xa <= x(ia) <= x(ib) <= xb <= x(ib+1)
!
    else

      ia = 1
      do i = 1, n-1
        if ( x(i) < xa ) then
          ia = i + 1
        end if
      end do
!
!  IA = 1 implies xa<x(1) .  Otherwise,
!  ia is largest index such that x(ia-1)<xa,.
!
      ib = n
      do i = n, ia, -1
        if ( xb < x(i) ) then
          ib = i - 1
        end if
      end do
!
!  IB = N implies X(N) < XB.  Otherwise,
!  ib is smallest index such that xb<x(ib+1) .
!
!  Compute the integral.
!
      ierv = 0
      if ( ib < ia ) then
!
!  This means IB = IA-1 and (A,B) is a subset of (x(ib),x(ia)).
!
           value = chfiv ( x(ib), x(ia), f(1,ib), f(1,ia), &
             d(1,ib), d(1,ia), a, b, ierv )

           if ( ierv < 0 ) then
             ierr = -4
             !call xerror ('pchia -- trouble in chfiv', ierr, 1)
             return
           end if

        else
!
!  First compute integral over (x(ia),x(ib)).
!
           if ( ib == ia ) then
              value = 0.0D+00
           else

              value = pchid ( n, x, f, d, incfd, skip, ia, ib, ierd )

              if ( ierd < 0 ) then
                ierr = -5
                !call xerror ('pchia -- trouble in pchid', ierr, 1)
                return
              end if

           end if
!
!  Then add on integral over ( XA, X(IA) ).
!
           if ( xa < x(ia) ) then
              il = max ( 1, ia-1 )
              ir = il + 1

              value = value + chfiv ( x(il), x(ir), f(1,il), f(1,ir), &
                d(1,il), d(1,ir), xa, x(ia), ierv )

              if ( ierv < 0 ) then
                ierr = -4
                !call xerror ('pchia -- trouble in chfiv', ierr, 1)
                return
              end if

           end if
!
!  Then add on integral over ( X(IB), XB ).
!
           if ( x(ib) < xb ) then
              ir = min ( ib+1, n )
              il = ir - 1

              value = value + chfiv ( x(il), x(ir), f(1,il), f(1,ir), &
                d(1,il), d(1,ir), x(ib), xb, ierv )

              if ( ierv < 0 ) then
                ierr = -4
                !call xerror ('pchia -- trouble in chfiv', ierr, 1)
                return
              end if

           end if
!
!  Adjust sign if necessary.
!
           if ( b < a ) then
             value = -value
           end if

        end if
     end if
  end if

  pchia = value

  return
end function

function pchqa ( n, x, f, d, a, b, ierr )

!*****************************************************************************80
!
!! PCHQA: easy to use cubic Hermite or spline integration.
!
!  Discussion:
!
!    PCHQA evaluates the definite integral of a cubic Hermite or spline
!    function over the interval [A, B].  This is an easy to use driver
!    for the routine PCHIA.
!
!  Author:
!
!    Fred Fritsch,
!    Mathematics and Statistics Division,
!    Lawrence Livermore National Laboratory.
!
!    FORTRAN90 translation by John Burkardt.
!
!  Reference:
!
!    David Kahaner, Cleve Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1989,
!    ISBN: 0-13-627258-4,
!    LC: TA345.K34.
!
!    Fred Fritsch, Ralph Carlson,
!    Monotone Piecewise Cubic Interpolation,
!    SIAM Journal on Numerical Analysis,
!    Volume 17, Number 2, April 1980, pages 238-246.
!
!    Fred Fritsch,
!    Piecewise Cubic Hermite Interpolation Package, Final Specifications,
!    Lawrence Livermore National Laboratory,
!    Computer Documentation UCID-30194, August 1982.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of data points.  N must be
!    at least 2.
!
!    Input, real ( kind = 8 ) X(N), the strictly increasing independent
!    variable values.
!
!    Input, real ( kind = 8 ) F(N), the function values.  F(I) is the value
!    corresponding to X(I).
!
!    Input, real ( kind = 8 ) D(N), the derivative values.  D(I) is the value
!    corresponding to X(I).
!
!    Input, real ( kind = 8 ) A, B, the limits of integration.  The
!    interval [A,B] is normally contained in [X(1),X(N)], but this is
!    not required.
!
!    Output, integer ( kind = 4 ) IERR, error flag.
!    0, no errors).
!    1, if A is outside the interval [X(1),X(N)].
!    2, if B is outside the interval [X(1),X(N)].
!    3, if both of the above are true.  This means that either [A,B] contains
!       the data interval or the intervals do not intersect at all.
!    -1, if N < 2 .
!    -3, if the X array is not strictly increasing.
!
!    Output, real ( kind = 8 ) PCHQA, the value of the requested integral.
!


  integer ( kind = 4 ) n

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) d(n)
  real ( kind = 8 ) f(n)
  integer ( kind = 4 ) ierr
  integer ( kind = 4 ), save :: incfd = 1
  !real ( kind = 8 ) pchia
  real ( kind = 8 ) pchqa
  logical, save :: skip = .true.
  real ( kind = 8 ) x(n)

  pchqa  =  pchia ( n, x, f, d, incfd, skip, a, b, ierr )


end function

subroutine pchsw ( dfmax, iextrm, d1, d2, h, slope, ierr )

!*****************************************************************************80
!
!! PCHSW: the PCHCS switch excursion limiter.
!
!  Discussion:
!
!    This routine is called by PCHCS to adjust D1 and D2 if necessary to
!    insure that the extremum on this interval is not further than DFMAX
!    from the extreme data value.
!
!  Author:
!
!    Fred Fritsch,
!    Mathematics and Statistics Division,
!    Lawrence Livermore National Laboratory.
!
!    FORTRAN90 translation by John Burkardt.
!
!  Reference:
!
!    Fred Fritsch, Ralph Carlson,
!    Monotone Piecewise Cubic Interpolation,
!    SIAM Journal on Numerical Analysis,
!    Volume 17, Number 2, April 1980, pages 238-246.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) DFMAX, the maximum allowed difference between
!    F(IEXTRM) and the cubic determined by the derivative values D1 and D2.
!    DFMAX should be nonnegative.
!
!    Input, integer ( kind = 4 ) IEXTRM, the index of the extreme data value,
!    which should be 1 or 2.
!
!    Input/output, real ( kind = 8 ) D1, D2, the derivative values at the
!    ends of the interval.  It is assumed that D1 * D2 <= 0.  On output,
!    the values may be modified if necessary to meet the restriction
!    imposed by DFMAX.
!
!    Input, real ( kind = 8 ) H, interval length.  H should be positive.
!
!    Input, real ( kind = 8 ) SLOPE, the data slope on the interval.
!
!    Output, integer ( kind = 4 ) IERR, error flag.
!    0, no errors.
!    -1, assumption on D1 and D2 is not satisfied.
!    -2, quadratic equation locating extremum has negative descriminant
!        (should never occur).
!
!  Local variables:
!
!    RHO is the ratio of the data slope to the derivative being tested.
!
!    LAMBDA is the ratio of D2 to D1.
!
!    THAT = T-hat(rho) is the normalized location of the extremum.
!
!    PHI is the normalized value of P(X)-f1 at X = xhat = x-hat(rho),
!      where  that = (xhat - x1)/h .
!      that is, p(xhat)-f1 = D * H * PHI,  where d=d1 or d2.
!      similarly,  p(xhat)-f2 = d*h*(phi-rho) .
!
!    SMALL should be a few orders of magnitude greater than macheps.
!


  real ( kind = 8 ) cp
  real ( kind = 8 ) d1
  real ( kind = 8 ) d2
  real ( kind = 8 ) dfmax
  real ( kind = 8 ) dmax
  real ( kind = 8 ), parameter :: fact = 100.0D+00
  real ( kind = 8 ) h
  integer ( kind = 4 ) ierr
  integer ( kind = 4 ) iextrm
  real ( kind = 8 ) lambda
  real ( kind = 8 ) nu
  real ( kind = 8 ) phi
  real ( kind = 8 ) radcal
  real ( kind = 8 ) rho
  real ( kind = 8 ) sigma
  real ( kind = 8 ) slope
  real ( kind = 8 ) small
  real ( kind = 8 ) that
  real ( kind = 8 ), parameter :: third = 0.33333D+00

  small = fact * epsilon ( 1.0D+00 )

  if ( d1 == 0.0D+00 ) then
!
!  Special case -- D1 == 0.0D+00 .
!
!  If D2 is also zero, this routine should not have been called.
!
     if ( d2 == 0.0D+00 ) then
       ierr = -1
      ! call xerror ('pchsw -- d1 and/or d2 invalid', ierr, 1)
       return
     end if

     rho = slope / d2
!
!  Extremum is outside interval when 1/3 <= RHO.
!
     if ( third <= rho ) then
       ierr = 0
       return
     end if

     that = ( 2.0D+00 * ( 3.0D+00 * rho - 1.0D+00 ) ) &
       / ( 3.0D+00 * ( 2.0D+00 * rho - 1.0D+00 ) )

     phi = that**2 * ( ( 3.0D+00 * rho - 1.0D+00 ) / 3.0D+00 )
!
!  Convert to distance from F2 if IEXTRM /= 1.
!
     if ( iextrm /= 1 ) then
       phi = phi - rho
     end if
!
!  Test for exceeding limit, and adjust accordingly.
!
     dmax = dfmax / ( h * abs ( phi ) )
     if ( dmax < abs ( d2 ) ) then
       d2 = sign ( dmax, d2 )
     end if

  else

     rho = slope / d1
     lambda = -d2 / d1
     if ( d2 == 0.0D+00 ) then
!
!  Special case -- D2 == 0.0D+00 .
!
!  Extremum is outside interval when 1/3 <= RHO.
!
        if ( third <= rho ) then
          ierr = 0
          return
        end if

        cp = 2.0D+00 - 3.0D+00 * rho
        nu = 1.0D+00 - 2.0D+00 * rho
        that = 1.0D+00 / ( 3.0D+00 * nu )

     else

        if ( lambda <= 0.0D+00 ) then
          ierr = -1
         ! call xerror ('pchsw -- d1 and/or d2 invalid', ierr, 1)
          return
        end if
!
!  Normal case, D1 and D2 both nonzero, opposite signs.
!
        nu = 1.0D+00 - lambda - 2.0D+00 * rho
        sigma = 1.0D+00 - rho
        cp = nu + sigma

        if ( small < abs ( nu ) ) then

          radcal = ( nu - ( 2.0D+00 * rho + 1.0D+00 ) ) * nu + sigma**2

          if ( radcal < 0.0D+00 ) then
            ierr = -2
           ! call xerror ( 'pchsw -- negative radical', ierr, 1)
            return
          end if

          that = ( cp - sqrt ( radcal ) ) / ( 3.0D+00 * nu )

        else

          that = 1.0D+00 / ( 2.0D+00 * sigma )

        end if

     end if

     phi = that * ( ( nu * that - cp ) * that + 1.0D+00 )
!
!  Convert to distance from F2 if IEXTRM /= 1.
!
     if ( iextrm /= 1 ) then
       phi = phi - rho
     end if
!
!  Test for exceeding limit, and adjust accordingly.
!
     dmax = dfmax / ( h * abs ( phi ) )

     if ( dmax < abs ( d1 ) ) then
        d1 = sign ( dmax, d1 )
        d2 = -lambda * d1
     end if

  end if

  ierr = 0

  return
end subroutine

subroutine chfdv ( x1, x2, f1, f2, d1, d2, ne, xe, fe, de, next, ierr )

!*****************************************************************************80
!
!! CHFDV evaluates a cubic polynomial and its derivative given in Hermite form.
!
!  Discussion:
!
!    CHFDV evaluates a cubic polynomial and its first derivative.
!    The cubic polynomial is given in Hermite form.  The evaluation
!    is carried out at an array of points.
!
!    This routine was designed for use by PCHFD, but it may also be
!    useful directly as an evaluator for a piecewise cubic Hermite
!    function in applications, such as graphing, where the interval
!    is known in advance.
!
!    If only function values are required, use CHFEV instead.
!
!  Author:
!
!    Fred Fritsch,
!    Mathematics and Statistics Division,
!    Lawrence Livermore National Laboratory.
!
!    FORTRAN90 translation by John Burkardt.
!
!  Reference:
!
!    Fred Fritsch, Ralph Carlson,
!    Monotone Piecewise Cubic Interpolation,
!    SIAM Journal on Numerical Analysis,
!    Volume 17, Number 2, April 1980, pages 238-246.
!
!    David Kahaner, Cleve Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1989,
!    ISBN: 0-13-627258-4,
!    LC: TA345.K34.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X1, X2, the endpoints of the interval of
!    definition of  the cubic.  X1 and X2 must be distinct.
!
!    Input, real ( kind = 8 ) F1, F2, the values of the function at X1 and
!    X2, respectively.
!
!    Input, real ( kind = 8 ) D1, D2, the derivative values at the ends
!     of the interval.
!
!    Input, integer ( kind = 4 ) NE, the number of evaluation points.
!
!    Input, real ( kind = 8 ) XE(NE), the points at which the functions are to
!    be evaluated.  If any of the XE are outside the interval
!    [X1,X2], a warning error is returned in next.
!
!    Output, real ( kind = 8 ) FE(NE), DE(NE), the values of the cubic
!    function and its derivative at the points XE(*).
!
!    Output, integer ( kind = 4 ) NEXT(2), indicates the number of
!    extrapolation points:
!    NEXT(1) = number of evaluation points to left of interval.
!    NEXT(2) = number of evaluation points to right of interval.
!
!    Output, integer ( kind = 4 ) IERR, error flag.
!    0, no errors.
!    -1, NE < 1.
!    -2, X1 == X2.
!


  integer ( kind = 4 ) ne

  real ( kind = 8 ) c2
  real ( kind = 8 ) c2t2
  real ( kind = 8 ) c3
  real ( kind = 8 ) c3t3
  real ( kind = 8 ) d1
  real ( kind = 8 ) d2
  real ( kind = 8 ) de(ne)
  real ( kind = 8 ) del1
  real ( kind = 8 ) del2
  real ( kind = 8 ) delta
  real ( kind = 8 ) f1
  real ( kind = 8 ) f2
  real ( kind = 8 ) fe(ne)
  real ( kind = 8 ) h
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierr
  integer ( kind = 4 ) next(2)
  real ( kind = 8 ) x
  real ( kind = 8 ) x1
  real ( kind = 8 ) x2
  real ( kind = 8 ) xe(ne)
  real ( kind = 8 ) xma
  real ( kind = 8 ) xmi
!
!  Check arguments.
!
  if ( ne < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'CHFDV - Fatal error!'
    write ( *, '(a)' ) '  The number of evaluation points was less than 1.'
    stop
  end if

  h = x2 - x1

  if ( h == 0.0D+00 ) then
    ierr = -2
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'CHFDV - Fatal error!'
    write ( *, '(a)' ) '  The interval endpoints are equal.'
    return
  end if
!
!  Initialize.
!
  ierr = 0
  next(1) = 0
  next(2) = 0
  xmi = min ( 0.0D+00, h )
  xma = max ( 0.0D+00, h )
!
!  Compute cubic coefficients expanded about X1.
!
  delta = ( f2 - f1 ) / h
  del1 = ( d1 - delta ) / h
  del2 = ( d2 - delta ) / h

  c2 = -( del1 + del1 + del2 )
  c2t2 = c2 + c2
  c3 = ( del1 + del2 ) / h
  c3t3 = c3 + c3 + c3
!
!  Evaluation loop.
!
  do i = 1, ne

    x = xe(i) - x1
    fe(i) = f1 + x * ( d1 + x * ( c2 + x * c3 ) )
    de(i) = d1 + x * ( c2t2 + x * c3t3 )
!
!  Count extrapolation points.
!
    if ( x < xmi ) then
      next(1) = next(1) + 1
    end if

    if ( xma < x ) then
      next(2) = next(2) + 1
    end if

  end do

  return
end subroutine

subroutine chfev ( x1, x2, f1, f2, d1, d2, ne, xe, fe, next, ierr )

!*****************************************************************************80
!
!! CHFEV evaluates a cubic polynomial given in Hermite form.
!
!  Discussion:
!
!    This routine evaluates a cubic polynomial given in Hermite form at an
!    array of points.  While designed for use by PCHFE, it may
!    be useful directly as an evaluator for a piecewise cubic
!    Hermite function in applications, such as graphing, where
!    the interval is known in advance.
!
!    The cubic polynomial is determined by function values
!    F1, F2 and derivatives D1, D2 on the interval [X1,X2].
!
!  Author:
!
!    Fred Fritsch,
!    Mathematics and Statistics Division,
!    Lawrence Livermore National Laboratory.
!
!    FORTRAN90 translation by John Burkardt.
!
!  Reference:
!
!    Fred Fritsch, Ralph Carlson,
!    Monotone Piecewise Cubic Interpolation,
!    SIAM Journal on Numerical Analysis,
!    Volume 17, Number 2, April 1980, pages 238-246.
!
!    David Kahaner, Cleve Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1989,
!    ISBN: 0-13-627258-4,
!    LC: TA345.K34.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X1, X2, the endpoints of the interval of
!    definition of the cubic.  X1 and X2 must be distinct.
!
!    Input, real ( kind = 8 ) F1, F2, the values of the function at X1 and
!    X2, respectively.
!
!    Input, real ( kind = 8 ) D1, D2, the derivative values at X1 and
!    X2, respectively.
!
!    Input, integer ( kind = 4 ) NE, the number of evaluation points.
!
!    Input, real ( kind = 8 ) XE(NE), the points at which the function is to
!    be evaluated.  If any of the XE are outside the interval
!    [X1,X2], a warning error is returned in NEXT.
!
!    Output, real ( kind = 8 ) FE(NE), the value of the cubic function
!    at the points XE.
!
!    Output, integer ( kind = 4 ) NEXT(2), the number of extrapolation points:
!    NEXT(1) = number of evaluation points to the left of interval.
!    NEXT(2) = number of evaluation points to the right of interval.
!
!    Output, integer ( kind = 4 ) IERR, error flag.
!    0, no errors.
!    -1, NE < 1.
!    -2, X1 == X2.
!


  integer ( kind = 4 ) ne

  real ( kind = 8 ) c2
  real ( kind = 8 ) c3
  real ( kind = 8 ) d1
  real ( kind = 8 ) d2
  real ( kind = 8 ) del1
  real ( kind = 8 ) del2
  real ( kind = 8 ) delta
  real ( kind = 8 ) f1
  real ( kind = 8 ) f2
  real ( kind = 8 ) fe(ne)
  real ( kind = 8 ) h
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierr
  integer ( kind = 4 ) next(2)
  real ( kind = 8 ) x
  real ( kind = 8 ) x1
  real ( kind = 8 ) x2
  real ( kind = 8 ) xe(ne)
  real ( kind = 8 ) xma
  real ( kind = 8 ) xmi

  if ( ne < 1 ) then
    ierr = -1
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'CHFEV - Fatal error!'
    write ( *, '(a)' ) '  Number of evaluation points is less than 1.'
    write ( *, '(a,i6)' ) '  NE = ', ne
    stop
  end if

  h = x2 - x1

  if ( h == 0.0D+00 ) then
    ierr = -2
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'CHFEV - Fatal error!'
    write ( *, '(a)' ) '  The interval [X1,X2] is of zero length.'
    stop
  end if
!
!  Initialize.
!
  ierr = 0
  next(1) = 0
  next(2) = 0
  xmi = min ( 0.0D+00, h )
  xma = max ( 0.0D+00, h )
!
!  Compute cubic coefficients expanded about X1.
!
  delta = ( f2 - f1 ) / h
  del1 = ( d1 - delta ) / h
  del2 = ( d2 - delta ) / h
  c2 = -( del1 + del1 + del2 )
  c3 = ( del1 + del2 ) / h
!
!  Evaluation loop.
!
  do i = 1, ne

    x = xe(i) - x1
    fe(i) = f1 + x * ( d1 + x * ( c2 + x * c3 ) )
!
!  Count the extrapolation points.
!
    if ( x < xmi ) then
      next(1) = next(1) + 1
    end if

    if ( xma < x ) then
      next(2) = next(2) + 1
    end if

  end do

  return
end subroutine



subroutine pchce ( ic, vc, n, x, h, slope, d, incfd, ierr )

!*****************************************************************************80
!
!! PCHCE is called by PCHIC to set end derivatives as requested by the user.
!
!  Discussion:
!
!    PCHCE must be called after interior derivative values have been set.
!
!    To facilitate two-dimensional applications, includes an increment
!    between successive values of the D array.
!
!    One could reduce the number of arguments and amount of local
!    storage, at the expense of reduced code clarity, by passing in
!    the array WK, rather than splitting it into H and SLOPE, and
!    increasing its length enough to incorporate STEMP and XTEMP.
!
!    The two monotonicity checks only use the sufficient conditions.
!    thus, it is possible (but unlikely) for a boundary condition to
!    be changed, even though the original interpolant was monotonic.
!    At least the result is a continuous function of the data.
!
!  Author:
!
!    FORTRAN77 original version by Fred Fritsch.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Fred Fritsch, Ralph Carlson,
!    Monotone Piecewise Cubic Interpolation,
!    SIAM Journal on Numerical Analysis,
!    Volume 17, Number 2, April 1980, pages 238-246.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) IC(2), specifies the desired boundary
!    conditions:
!    IC(1) = IBEG, desired condition at beginning of data.
!    IC(2) = IEND, desired condition at end of data.
!    See the prologue to PCHIC for details.
!
!    Input, real ( kind = 8 ) VC(2), specifies desired boundary values, as
!    indicated above.  VC(1) need be set only if IC(1) = 2 or 3.
!    VC(2) need be set only if IC(2) = 2 or 3.
!
!    Input, integer ( kind = 4 ) N, the number of data points.
!    N must be at least 2.
!
!    Input, real ( kind = 8 ) X(N), the strictly increasing independent
!    variable values.
!
!    Input, real ( kind = 8 ) H(N), interval lengths.  H(I) = X(I+1)-X(I),
!    for I = 1 to N-1.
!
!    Input, real ( kind = 8 ) SLOPE(N), the data slopes.
!    SLOPE(I) = ( Y(I+1) - Y(I) ) / H(I), for I = 1 to N-1.
!
!    Input/output, real ( kind = 8 ) D(INCFD,N), the derivative values at the
!    data points.  The value corresponding to X(I) must be stored in
!    D(1+(I-1)*INCFD).  On output, the value of D at X(1) and/or X(N) is
!    changed, if necessary, to produce the requested boundary conditions.
!
!    Input, integer ( kind = 4 ) INCFD, increment between successive values
!    in D.  This argument is provided primarily for 2-d applications.
!
!    Output, integer ( kind = 4 ) IERR, error flag.
!    0, no errors.
!    1, if IBEG < 0 and D(1) had to be adjusted for monotonicity.
!    2, if IEND < 0 and D(1+(N-1)*INCFD) had to be adjusted for monotonicity.
!    3, if both of the above are true.
!


  integer ( kind = 4 ) incfd
  integer ( kind = 4 ) n

  real ( kind = 8 ) d(incfd,n)
  real ( kind = 8 ) h(n)
  integer ( kind = 4 ) ibeg
  integer ( kind = 4 ) ic(2)
  integer ( kind = 4 ) iend
  integer ( kind = 4 ) ierf
  integer ( kind = 4 ) ierr
  integer ( kind = 4 ) index
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
 ! real ( kind = 8 ) pchdf
 ! real ( kind = 8 ) pchst
  real ( kind = 8 ) slope(n)
  real ( kind = 8 ) stemp(3)
  real ( kind = 8 ) vc(2)
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) xtemp(4)

  ibeg = ic(1)
  iend = ic(2)
  ierr = 0
!
!  Set to default boundary conditions if N is too small.
!
  if ( n < abs ( ibeg ) ) then
    ibeg = 0
  end if

  if ( n < abs ( iend ) ) then
    iend = 0
  end if
!
!  Treat beginning boundary condition.
!
  if ( ibeg == 0 ) then
    go to 2000
  end if

  k = abs ( ibeg )

  if ( k == 1 ) then
!
!  Boundary value provided.
!
     d(1,1) = vc(1)
  else if ( k == 2 ) then
!
!  Boundary second derivative provided.
!
     d(1,1) = 0.5D+00 * &
      ( ( 3.0D+00 * slope(1) - d(1,2) ) - 0.5D+00 * vc(1) * h(1) )

  else if ( k < 5 ) then
!
!  Use K-point derivative formula.
!  Pick up first K points, in reverse order.
!
     do j = 1, k
       index = k-j+1
       xtemp(j) = x(index)
       if ( j < k ) then
         stemp(j) = slope(index-1)
       end if
     end do

     d(1,1) = pchdf ( k, xtemp, stemp, ierf )

     if ( ierf /= 0 ) then
       ierr = -1
       !call xerror ('PCHCE -- error return from pchdf', ierr, 1)
       return
     end if

  else
!
!  Use 'not a knot' condition.
!
     d(1,1) = ( 3.0D+00 * ( h(1) * slope(2) + h(2) * slope(1) ) &
       - 2.0D+00 * ( h(1) + h(2) ) * d(1,2) - h(1) * d(1,3) ) / h(2)
  end if
!
!  Check d(1,1) for compatibility with monotonicity.
!
  if ( ibeg <= 0 ) then

    if ( slope(1) == 0.0D+00 ) then
      if ( d(1,1) /= 0.0D+00 ) then
        d(1,1) = 0.0D+00
        ierr = ierr + 1
      end if
    else if ( pchst ( d(1,1), slope(1) ) < 0.0D+00 ) then
      d(1,1) = 0.0D+00
      ierr = ierr + 1
    else if ( 3.0D+00 * abs ( slope(1) ) < abs ( d(1,1) ) ) then
      d(1,1) = 3.0D+00 * slope(1)
      ierr = ierr + 1
    end if

  end if

2000 continue
!
!  Treat end boundary condition.
!
  if ( iend == 0 ) then
    return
  end if

  k = abs ( iend )
  if ( k == 1 ) then
!
!  Boundary value provided.
!
     d(1,n) = vc(2)

  else if ( k == 2 ) then
!
!  Boundary second derivative provided.
!
     d(1,n) = 0.5D+00 * ( (3.0D+00 * slope(n-1) - d(1,n-1)) &
       + 0.5D+00 * vc(2) * h(n-1) )

  else if ( k < 5 ) then
!
!  Use K-point derivative formula.  Pick up last K points.
!
     do j = 1, k
       index = n - k + j
       xtemp(j) = x(index)
       if ( j < k ) then
         stemp(j) = slope(index)
       end if
     end do

     d(1,n) = pchdf ( k, xtemp, stemp, ierf )

     if ( ierf /= 0 ) then
       ierr = -1
       !call xerror ('pchce -- error return from pchdf', ierr, 1)
       return
     end if

  else
!
!  Use 'not a knot' condition.
!
     d(1,n) = ( 3.0D+00 * ( h(n-1) * slope(n-2) + h(n-2) * slope(n-1) ) &
       - 2.0D+00 * ( h(n-1) + h(n-2)) * d(1,n-1) - h(n-1) * d(1,n-2) ) / h(n-2)
  end if

  if ( 0 < iend ) then
    return
  end if
!
!  Check D(1,n) for compatibility with monotonicity.
!
  if ( slope(n-1) == 0.0D+00 ) then
    if ( d(1,n) /= 0.0D+00 ) then
       d(1,n) = 0.0D+00
       ierr = ierr + 2
    end if
  else if ( pchst ( d(1,n), slope(n-1) ) < 0.0D+00 ) then
    d(1,n) = 0.0D+00
    ierr = ierr + 2
  else if ( 3.0D+00 * abs ( slope(n-1) ) < abs ( d(1,n) ) ) then
    d(1,n) = 3.0D+00 * slope(n-1)
    ierr = ierr + 2
  end if

  return
end subroutine


subroutine pchci ( n, h, slope, d, incfd )

!*****************************************************************************80
!
!! PCHCI sets derivatives for a monotone piecewise cubic Hermite interpolant.
!
!  Discussion:
!
!    Default boundary conditions are provided which are compatible
!    with monotonicity.  If the data are only piecewise monotonic, the
!    interpolant will have an extremum at each point where monotonicity
!    switches direction.
!
!    To facilitate two-dimensional applications, includes an increment
!    between successive values of the D array.
!
!    The resulting piecewise cubic Hermite function should be identical
!    (within roundoff error) to that produced by PCHIM.
!
!  Author:
!
!    Fred Fritsch,
!    Mathematics and Statistics Division,
!    Lawrence Livermore National Laboratory.
!
!    FORTRAN90 translation by John Burkardt.
!
!  Reference:
!
!    Fred Fritsch, Ralph Carlson,
!    Monotone Piecewise Cubic Interpolation,
!    SIAM Journal on Numerical Analysis,
!    Volume 17, Number 2, April 1980, pages 238-246.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of data points.
!    N must be at least 2.
!
!    Input, real ( kind = 8 ) H(N), interval lengths.  H(I) = X(I+1)-X(I),
!    for I = 1 to N-1.
!
!    Input, real ( kind = 8 ) SLOPE(N), the data slopes.
!    SLOPE(I) = ( Y(I+1) - Y(I) ) / H(I), for I = 1 to N-1.
!
!    Output, real ( kind = 8 ) D(INCFD,N), the derivative values at the data
!    points.  If the data are monotonic, these values will determine a monotone
!    cubic Hermite function.  The value corresponding to X(I) is stored in
!    D(1+(I-1)*INCFD).
!
!    Input, integer ( kind = 4 ) INCFD, increment between successive values
!    in D.  This argument is provided primarily for 2D applications.
!


  integer ( kind = 4 ) incfd
  integer ( kind = 4 ) n

  real ( kind = 8 ) d(incfd,n)
  real ( kind = 8 ) del1
  real ( kind = 8 ) del2
  real ( kind = 8 ) dmax
  real ( kind = 8 ) dmin
  real ( kind = 8 ) drat1
  real ( kind = 8 ) drat2
  real ( kind = 8 ) h(n)
  real ( kind = 8 ) hsum
  real ( kind = 8 ) hsumt3
  integer ( kind = 4 ) i
  integer ( kind = 4 ) nless1
  !real ( kind = 8 ) pchst
  real ( kind = 8 ) slope(n)
  real ( kind = 8 ) w1
  real ( kind = 8 ) w2

  nless1 = n - 1
  del1 = slope(1)
!
!  Special case N=2 -- use linear interpolation.
!
  if ( nless1 <= 1 ) then
    d(1,1) = del1
    d(1,n) = del1
    return
  end if
!
!  Normal case, 3 <= N.
!
  del2 = slope(2)
!
!  Set D(1) via non-centered three point formula, adjusted to be
!  shape preserving.
!
  hsum = h(1) + h(2)
  w1 = ( h(1) + hsum ) / hsum
  w2 = -h(1) / hsum
  d(1,1) = w1 * del1 + w2 * del2

  if ( pchst ( d(1,1), del1 ) <= 0.0D+00 ) then
    d(1,1) = 0.0D+00
  else if ( pchst ( del1, del2 ) < 0.0D+00 ) then
!
!  Need do this check only if monotonicity switches.
!
    dmax = 3.0D+00 * del1

    if ( abs ( dmax ) < abs ( d(1,1) ) ) then
      d(1,1) = dmax
    end if

  end if
!
!  Loop through interior points.
!
  do i = 2, nless1

    if ( i /= 2 ) then
      hsum = h(i-1) + h(i)
      del1 = del2
      del2 = slope(i)
    end if
!
!  Set D(I)=0 unless data are strictly monotonic.
!
    d(1,i) = 0.0D+00
!
!  Use Brodlie modification of Butland formula.
!
    if ( 0.0D+00 < pchst ( del1, del2 ) ) then

      hsumt3 = hsum + hsum + hsum
      w1 = ( hsum + h(i-1)) / hsumt3
      w2 = ( hsum + h(i)  ) / hsumt3
      dmax = max ( abs ( del1 ), abs ( del2 ) )
      dmin = min ( abs ( del1 ), abs ( del2 ) )
      drat1 = del1 / dmax
      drat2 = del2 / dmax
      d(1,i) = dmin / ( w1 * drat1 + w2 * drat2 )

    end if

  end do
!
!  Set D(N) via non-centered three point formula, adjusted to
!  be shape preserving.
!
  w1 = -h(n-1) / hsum
  w2 = ( h(n-1) + hsum ) / hsum
  d(1,n) = w1 * del1 + w2 * del2

  if ( pchst ( d(1,n), del2 ) <= 0.0D+00 ) then
    d(1,n) = 0.0D+00
  else if ( pchst ( del1, del2 ) < 0.0D+00 ) then
    dmax = 3.0D+00 * del2
    if ( abs ( dmax ) < abs ( d(1,n) ) ) then
      d(1,n) = dmax
    end if
  end if

  return
end subroutine

subroutine pchcs ( switch, n, h, slope, d, incfd, ierr )

!*****************************************************************************80
!
!! PCHCS adjusts the curve produced by PCHIM so it is more "visually pleasing".
!
!  Discussion:
!
!    PCHCS is called by PCHIC to adjust the values of D in the vicinity of a
!    switch in direction of monotonicity, to produce a more "visually
!    pleasing" curve than that given by PCHIM.
!
!  Author:
!
!    Fred Fritsch,
!    Mathematics and Statistics Division,
!    Lawrence Livermore National Laboratory.
!
!    FORTRAN90 translation by John Burkardt.
!
!  Reference:
!
!    Fred Fritsch, Ralph Carlson,
!    Monotone Piecewise Cubic Interpolation,
!    SIAM Journal on Numerical Analysis,
!    Volume 17, Number 2, April 1980, pages 238-246.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) SWITCH, indicates the amount of control desired
!    over local excursions from data.
!
!    Input, integer ( kind = 4 ) N, the number of data points.
!    N must be at least 3.
!
!    Input, real ( kind = 8 ) H(N), interval lengths.  H(I) = X(I+1)-X(I),
!    for I = 1 to N-1.
!
!    Input, real ( kind = 8 ) SLOPE(N), the data slopes.
!    SLOPE(I) = ( Y(I+1) - Y(I) ) / H(I), for I = 1 to N-1.
!
!    Input/output, real ( kind = 8 ) D(N).  On input, the derivative values at
!    the data points, as determined by PCHIC.  On output, derivatives in the
!    vicinity of switches in direction of monotonicity may be adjusted to
!    produce a more "visually pleasing" curve.  The value corresponding to
!    X(I) is stored in D(1+(I-1)*INCFD).
!
!    Input, integer ( kind = 4 ) INCFD, increment between successive values
!    in D.  This argument is provided primarily for 2D applications.
!
!    Output, integer ( kind = 4 ) IERR, error flag.
!    0, no errors.
!    negative, trouble in PCHSW.
!


  integer ( kind = 4 ) incfd
  integer ( kind = 4 ) n

  real ( kind = 8 ) d(incfd,n)
  real ( kind = 8 ) del(3)
  real ( kind = 8 ) dext
  real ( kind = 8 ) dfloc
  real ( kind = 8 ) dfmx
  real ( kind = 8 ) fact
  real ( kind = 8 ) h(n)
  real ( kind = 8 ), parameter :: fudge = 4.0D+00
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierr
  integer ( kind = 4 ) indx
  integer ( kind = 4 ) k
  integer ( kind = 4 ) nless1
 ! real ( kind = 8 ) pchst
  real ( kind = 8 ) slmax
  real ( kind = 8 ) slope(n)
  real ( kind = 8 ) switch
  real ( kind = 8 ) wtave(2)
!
!  Define inline function for weighted average of slopes.
!
  real ( kind = 8 ) pchsd, s1, s2, h1, h2

  pchsd ( s1, s2, h1, h2 ) = ( h2 / ( h1 + h2 ) ) &
    * s1 + ( h1 / ( h1 + h2 ) ) * s2
!
!  Initialize.
!
  ierr = 0
  nless1 = n - 1
!
!  Loop over segments.
!
  do i = 2, nless1

     if ( pchst ( slope(i-1), slope(i) ) )  100, 300, 900

  100    continue
!
!  Slope switches monotonicity at i-th point
!
!  Do not change D if 'up-down-up'.
!
        if ( 2 < i ) then
          if ( 0.0D+00 < pchst ( slope(i-2), slope(i) ) ) then
            cycle
          end if
        end if

        if ( i < nless1 ) then
          if ( 0.0D+00 < pchst ( slope(i+1), slope(i-1) ) ) then
            cycle
          end if
        end if
!
!  Compute provisional value for D(1,i).
!
        dext = pchsd ( slope(i-1), slope(i), h(i-1), h(i) )
!
!  Determine which interval contains the extremum.
!
        if ( pchst ( dext, slope(i-1) ) )  200, 900, 250

  200       continue
!
!  DEXT and slope(i-1) have opposite signs.
!  extremum is in (x(i-1),x(i)).
!
           k = i - 1
!
!  Set up to compute new values for D(1,i-1) and D(1,i).
!
           wtave(2) = dext
           if ( 1 < k ) then
             wtave(1) = pchsd (slope(k-1), slope(k), h(k-1), h(k))
           end if
           go to 400

  250       continue
!
!  DEXT and SLOPE(I) have opposite signs.
!  The extremum is in (x(i),x(i+1)).
!
           k = i
!
!  Set up to compute new values for D(1,i) and D(1,i+1).
!
           wtave(1) = dext
           if ( k < nless1 ) then
             wtave(2) = pchsd ( slope(k), slope(k+1), h(k), h(k+1) )
           end if
           go to 400

  300    continue
!
!  At least one of SLOPE(I-1) and slope(i) is zero.
!  Check for flat-topped peak
!
        if ( i == nless1 ) then
          cycle
        end if

        if ( 0.0D+00 <= pchst ( slope(i-1), slope(i+1) ) ) then
          cycle
        end if
!
!  We have flat-topped peak on (x(i),x(i+1)).
!
        k = i
!
!  Set up to compute new values for d(1,i) and d(1,i+1).
!
        wtave(1) = pchsd ( slope(k-1), slope(k), h(k-1), h(k) )
        wtave(2) = pchsd ( slope(k), slope(k+1), h(k), h(k+1) )

  400    continue
!
!  At this point we have determined that there will be an extremum
!  on (x(k),x(k+1)), where k=i or i-1, and have set array WTAVE.
!  wtave(1) is a weighted average of slope(k-1) and slope(k), if 1 < K
!  wtave(2) is a weighted average of slope(k) and slope(k+1), if k<n-1
!
     slmax = abs ( slope(k) )
     if ( 1 < k ) then
       slmax = max ( slmax, abs ( slope(k-1) ) )
     end if

     if ( k < nless1 ) then
       slmax = max ( slmax, abs ( slope(k+1) ) )
     end if

     if ( 1 < k ) then
       del(1) = slope(k-1) / slmax
     end if

     del(2) = slope(k) / slmax

     if ( k < nless1 ) then
       del(3) = slope(k+1) / slmax
     end if

     if ( 1 < k .and. k < nless1 ) then
!
!  Normal case -- extremum is not in a boundary interval.
!
        fact = fudge * abs ( del(3) * ( del(1) - del(2) ) &
          * ( wtave(2) / slmax ) )

        d(1,k) = d(1,k) + min ( fact, 1.0D+00 ) * ( wtave(1) - d(1,k) )
        fact = fudge * abs ( del(1) * ( del(3) - del(2) ) &
          * ( wtave(1) / slmax ) )
        d(1,k+1) = d(1,k+1) + min ( fact, 1.0D+00 ) * ( wtave(2) - d(1,k+1) )
     else
!
!  Special case K=1 (which can occur only if I=2) or
!  k=nless1 (which can occur only if i=nless1).
!
        fact = fudge * abs ( del(2) )
        d(1,i) = min ( fact, 1.0D+00 ) * wtave(i-k+1)
!
!  Note that i-k+1 = 1 if k=i  (=nless1),
!            i-k+1 = 2 if k=i-1(=1).
!
     end if
!
!  Adjust if necessary to limit excursions from data.
!
     if ( switch <= 0.0D+00 ) then
       cycle
     end if

     dfloc = h(k) * abs ( slope(k) )

     if ( 1 < k ) then
       dfloc = max ( dfloc, h(k-1) * abs ( slope(k-1) ) )
     end if

     if ( k < nless1 ) then
       dfloc = max ( dfloc, h(k+1) * abs ( slope(k+1) ) )
     end if

     dfmx = switch * dfloc
     indx = i-k+1
!
!  INDX = 1 if K = I,
!  INDX = 2 if K = I-1.
!
     call pchsw ( dfmx, indx, d(1,k), d(1,k+1), h(k), slope(k), ierr )

     if ( ierr /= 0 ) then
       return
     end if

  900 continue

  end do

  return
end subroutine




subroutine pchev ( n, x, f, d, nval, xval, fval, dval, ierr )

!*****************************************************************************80
!
!! PCHEV evaluates a piecewise cubic Hermite or spline function.
!
!  Discussion:
!
!    PCHEV evaluates the function and first derivative of a piecewise
!    cubic Hermite or spline function at an array of points XVAL.
!
!    The evaluation will be most efficient if the elements of XVAL are
!    increasing relative to X; that is, for all J <= K,
!      X(I) <= XVAL(J)
!    implies
!      X(I) <= XVAL(K).
!
!    If any of the XVAL are outside the interval [X(1),X(N)],
!    values are extrapolated from the nearest extreme cubic,
!    and a warning error is returned.
!
!  Author:
!
!    Fred Fritsch,
!    Mathematics and Statistics Division,
!    Lawrence Livermore National Laboratory.
!
!    FORTRAN90 translation by John Burkardt.
!
!  Reference:
!
!    David Kahaner, Cleve Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1989,
!    ISBN: 0-13-627258-4,
!    LC: TA345.K34.
!
!    Fred Fritsch, Ralph Carlson,
!    Monotone Piecewise Cubic Interpolation,
!    SIAM Journal on Numerical Analysis,
!    Volume 17, Number 2, April 1980, pages 238-246.
!
!    Fred Fritsch,
!    Piecewise Cubic Hermite Interpolation Package, Final Specifications,
!    Lawrence Livermore National Laboratory,
!    Computer Documentation UCID-30194, August 1982.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of data points.
!    N must be at least 2.
!
!    Input, real ( kind = 8 ) X(N), the strictly increasing independent
!    variable values.
!
!    Input, real ( kind = 8 ) F(N), the function values.  F(I) is the value
!    corresponding to X(I).
!
!    Input, real ( kind = 8 ) D(N), the derivative values.  D(i) is the value
!    corresponding to X(I).
!
!    Input, integer ( kind = 4 ) NVAL, the number of points at which the
!    functions are to be evaluated.
!
!    Input, real ( kind = 8 ) XVAL(NVAL), the points at which the functions
!    are to be evaluated.
!
!    Output, real ( kind = 8 ) FVAL(NVAL), the values of the cubic Hermite
!    function at XVAL.
!
!    Output, real ( kind = 8 ) DVAL(NVAL), the derivatives of the cubic
!    Hermite function at XVAL.
!
!    Output, integer ( kind = 4 ) IERR, error flag.
!    0, no errors.
!    positive, means that extrapolation was performed at IERR points.
!    -1, if N < 2.
!    -3, if the X array is not strictly increasing.
!    -4, if NVAL < 1.
!    -5, if an error has occurred in CHFDV.
!


  integer ( kind = 4 ) n
  integer ( kind = 4 ) nval

  real ( kind = 8 ) d(n)
  real ( kind = 8 ) dval(nval)
  real ( kind = 8 ) f(n)
  real ( kind = 8 ) fval(nval)
  integer ( kind = 4 ) ierr
  integer ( kind = 4 ), save :: incfd = 1
  logical, save :: skip = .true.
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) xval(nval)

  call pchfd ( n, x, f, d, incfd, skip, nval, xval, fval, dval, ierr )

  return
end subroutine


subroutine pchez ( n, x, f, d, spline, wk, lwk, ierr )

!*****************************************************************************80
!
!! PCHEZ carries out easy to use spline or cubic Hermite interpolation.
!
!  Discussion:
!
!    This routine sets derivatives for spline (two continuous derivatives)
!    or Hermite cubic (one continuous derivative) interpolation.
!    Spline interpolation is smoother, but may not "look" right if the
!    data contains both "steep" and "flat" sections.  Hermite cubics
!    can produce a "visually pleasing" and monotone interpolant to
!    monotone data.
!
!    This routine is an easy to use driver for the PCHIP routines.
!    Various boundary conditions are set to default values by PCHEZ.
!    Many other choices are available in the subroutines PCHIC,
!    PCHIM and PCHSP.
!
!    Use PCHEV to evaluate the resulting function and its derivative.
!
!    If SPLINE is TRUE, the interpolating spline satisfies the default
!    "not-a-knot" boundary condition, with a continuous third derivative
!    at X(2) and X(N-1).
!
!    If SPLINE is FALSE, the interpolating Hermite cubic will be monotone
!    if the input data is monotone.  Boundary conditions are computed from
!    the derivative of a local quadratic unless this alters monotonicity.
!
!  Author:
!
!    Fred Fritsch,
!    Mathematics and Statistics Division,
!    Lawrence Livermore National Laboratory.
!
!    FORTRAN90 translation by John Burkardt.
!
!  Reference:
!
!    David Kahaner, Cleve Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1989,
!    ISBN: 0-13-627258-4,
!    LC: TA345.K34.
!
!    Fred Fritsch, Ralph Carlson,
!    Monotone Piecewise Cubic Interpolation,
!    SIAM Journal on Numerical Analysis,
!    Volume 17, Number 2, April 1980, pages 238-246.
!
!    Fred Fritsch, Judy Butland,
!    A Method for Constructing Local Monotone Piecewise
!    Cubic Interpolants,
!    SIAM Journal on Scientific and Statistical Computing,
!    Volume 5, Number 2, 1984, pages 300-304.
!
!    Carl deBoor,
!    A Practical Guide to Splines, Chapter IV,
!    Springer-Verlag,
!    New York, 1978.
!
!    Fred Fritsch,
!    Piecewise Cubic Hermite Interpolation Package, Final Specifications,
!    Lawrence Livermore National Laboratory,
!    Computer Documentation UCID-30194, August 1982.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of data points.
!    N must be at least 2.
!
!    Input, real ( kind = 8 ) X(N), the strictly increasing independent
!    variable values.
!
!    Input, real ( kind = 8 ) F(N), the function values.  F(I) is the value
!    corresponding to X(I).
!
!    Output, real ( kind = 8 ) D(N), the derivative values at the data points.
!
!    Input, logical SPLINE, specifies if the interpolant is to be a spline
!    with two continuous derivatives (SPLINE is TRUE), or a Hermite cubic
!    interpolant with one continuous derivative (SPLINE is FALSE).
!
!    Workspace, real ( kind = 8 ) WK(LWK), required only if SPLINE is TRUE.
!
!    Input, integer ( kind = 4 ) LWK, the length of the work array WK, which
!    must be at least 2*N.  However, WK is not needed if SPLINE is FALSE,
!    and in this case LWK is arbitrary.
!
!    Output, integer ( kind = 4 ) IERR, error flag.
!    0, no errors.
!    positive, can only occur when SPLINE is FALSE,  means that there were
!      IERR switches in the direction of monotonicity.  When SPLINE is
!      FALSE, PCHEZ guarantees that if the input data is monotone, the
!      interpolant will be too.  This warning is to alert you to the fact
!      that the input data was not monotone.
!    -1, if N < 2.
!    -3, if the X array is not strictly increasing.
!    -7, if LWK is less than 2*N and SPLINE is TRUE.
!


  integer ( kind = 4 ) lwk
  integer ( kind = 4 ) n

  real ( kind = 8 ) d(n)
  real ( kind = 8 ) f(n)
  integer ( kind = 4 ), save, dimension ( 2 ) :: ic = (/ 0, 0 /)
  integer ( kind = 4 ) ierr
  integer ( kind = 4 ), parameter :: incfd = 1
  logical spline
  real ( kind = 8 ) vc(2)
  real ( kind = 8 ) wk(lwk)
  real ( kind = 8 ) x(n)

  if ( spline ) then
    call pchsp ( ic, vc, n, x, f, d, incfd, wk, lwk, ierr )
  else
    call pchim ( n, x, f, d, incfd, ierr )
  end if

  return
end subroutine


subroutine pchfd ( n, x, f, d, incfd, skip, ne, xe, fe, de, ierr )

!*****************************************************************************80
!
!! PCHFD evaluates a piecewise cubic Hermite function.
!
!  Discsussion:
!
!    PCHFD evaluates a piecewise cubic Hermite function and its first
!    derivative at an array of points.  PCHFD may be used by itself
!    for Hermite interpolation, or as an evaluator for PCHIM
!    or PCHIC.
!
!    PCHFD evaluates the cubic Hermite function and its first derivative
!    at the points XE.
!
!    If only function values are required, use PCHFE instead.
!
!    To provide compatibility with PCHIM and PCHIC, includes an
!    increment between successive values of the F and D arrays.
!
!    Most of the coding between the call to CHFDV and the end of
!    the IR loop could be eliminated if it were permissible to
!    assume that XE is ordered relative to X.
!
!    CHFDV does not assume that X1 is less than X2.  Thus, it would
!    be possible to write a version of PCHFD that assumes a strictly
!    decreasing X array by simply running the IR loop backwards
!    and reversing the order of appropriate tests.
!
!    The present code has a minor bug, which I have decided is not
!    worth the effort that would be required to fix it.
!    If XE contains points in [X(N-1),X(N)], followed by points less than
!    X(N-1), followed by points greater than X(N), the extrapolation points
!    will be counted (at least) twice in the total returned in IERR.
!
!    The evaluation will be most efficient if the elements of XE are
!    increasing relative to X; that is, for all J <= K,
!      X(I) <= XE(J)
!    implies
!      X(I) <= XE(K).
!
!    If any of the XE are outside the interval [X(1),X(N)],
!    values are extrapolated from the nearest extreme cubic,
!    and a warning error is returned.
!
!  Modified:
!
!    13 August 2005
!
!  Author:
!
!    Fred Fritsch,
!    Mathematics and Statistics Division,
!    Lawrence Livermore National Laboratory.
!
!    FORTRAN90 translation by John Burkardt.
!
!  Reference:
!
!    Fred Fritsch, Ralph Carlson,
!    Monotone Piecewise Cubic Interpolation,
!    SIAM Journal on Numerical Analysis,
!    Volume 17, Number 2, April 1980, pages 238-246.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of data points.
!    N must be at least 2.
!
!    Input, real ( kind = 8 ) X(N), the strictly increasing independent
!    variable values.
!
!    Input, real ( kind = 8 ) F(INCFD,N), the function values.
!    F(1+(I-1)*INCFD) is the value corresponding to X(I).
!
!    Input, real ( kind = 8 ) D(INCFD,N), the derivative values.
!    D(1+(I-1)*INCFD) is the value corresponding to X(I).
!
!    Input, integer ( kind = 4 ) INCFD, increment between successive values in
!    F and D.
!
!    Input/output, logical SKIP, controls whether data validity checks
!    should be made.  Setting the input value to FALSE will skip the checks.
!    On output with 0 <= IERR, SKIP will be set to TRUE.
!
!    Input, integer ( kind = 4 ) NE, the number of evaluation points.
!
!    Input, real ( kind = 8 ) XE(NE), points at which the function is
!    to be evaluated.
!
!    Output, real ( kind = 8 ) FE(NE), the values of the cubic Hermite
!    function at XE.
!
!    Output, real ( kind = 8 ) DE(NE), the derivative of the cubic
!    Hermite function at XE.
!
!    Output, integer ( kind = 4 ) IERR, error flag.
!    0, no errors.
!    positive, means that extrapolation was performed at IERR points.
!    -1, if N < 2.
!    -2, if INCFD < 1.
!    -3, if the X array is not strictly increasing.
!    -4, if NE < 1.
!    -5, if an error has occurred in the lower-level routine CHFDV.
!


  integer ( kind = 4 ) incfd
  integer ( kind = 4 ) n
  integer ( kind = 4 ) ne

  real ( kind = 8 ) d(incfd,n)
  real ( kind = 8 ) de(ne)
  real ( kind = 8 ) f(incfd,n)
  real ( kind = 8 ) fe(ne)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierc
  integer ( kind = 4 ) ierr
  integer ( kind = 4 ) ir
  integer ( kind = 4 ) j
  integer ( kind = 4 ) j_first
  integer ( kind = 4 ) j_new
  integer ( kind = 4 ) j_save
  integer ( kind = 4 ) next(2)
  integer ( kind = 4 ) nj
  logical skip
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) xe(ne)
!
!  Check arguments.
!
  if ( .not. skip ) then

    if ( n < 2 ) then
      ierr = -1
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'PCHFD - Fatal error!'
      write ( *, '(a)' ) '  Number of data points less than 2.'
      return
    end if

    if ( incfd < 1 ) then
      ierr = -2
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'PCHFD - Fatal error!'
      write ( *, '(a)' ) '  Increment less than 1.'
      return
    end if

    do i = 2, n
      if ( x(i) <= x(i-1) ) then
        ierr = -3
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'PCHFD - Fatal error!'
        write ( *, '(a)' ) '  X array not strictly increasing.'
        return
      end if
    end do

  end if

  if ( ne < 1 ) then
    ierr = -4
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PCHFD - Fatal error!'
    write ( *, '(a)' ) '  Number of evaluation points less than 1.'
    return
  end if

  ierr = 0
  skip = .true.
!
!  Loop over intervals.
!  The interval index is IL = IR - 1.
!  The interval is X(IL) <= X < X(IR).
!
  j_first = 1
  ir = 2

  do
!
!  Skip out of loop if have processed all evaluation points.
!
    if ( ne < j_first ) then
      exit
    end if
!
!  Locate all points in interval.
!
    j_save = ne + 1

    do j = j_first, ne
      if ( x(ir) <= xe(j) ) then
        j_save = j
        if ( ir == n ) then
          j_save = ne + 1
        end if
        exit
      end if
    end do
!
!  Have located first point beyond interval.
!
    j = j_save

    nj = j - j_first
!
!  Skip evaluation if no points in interval.
!
    if ( nj /= 0 ) then
!
!  Evaluate cubic at XE(J_FIRST:J-1).
!
      call chfdv ( x(ir-1), x(ir), f(1,ir-1), f(1,ir), d(1,ir-1), d(1,ir), &
        nj, xe(j_first), fe(j_first), de(j_first), next, ierc )

      if ( ierc < 0 ) then
        ierr = -5
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'PCHFD - Fatal error!'
        write ( *, '(a)' ) '  Error return from CHFDV.'
        return
      end if
!
!  In the current set of XE points, there are NEXT(2) to the right of X(IR).
!
      if ( next(2) /= 0 ) then

        if ( ir < n ) then
          ierr = -5
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'PCHFD - Fatal error!'
          write ( *, '(a)' ) '  IR < N.'
          return
        end if
!
!  These are actually extrapolation points.
!
        ierr = ierr + next(2)

      end if
!
!  In the current set of XE points, there are NEXT(1) to the left of X(IR-1).
!
      if ( next(1) /= 0 ) then
!
!  These are actually extrapolation points.
!
        if ( ir <= 2 ) then
          ierr = ierr + next(1)
!
!  XE is not ordered relative to X, so must adjust evaluation interval.
!
!  First, locate first point to left of X(IR-1).
!
        else

          j_new = -1

          do i = j_first, j-1
            if ( xe(i) < x(ir-1) ) then
              j_new = i
              exit
            end if
          end do

          if ( j_new == -1 ) then
            ierr = -5
            write ( *, '(a)' ) ' '
            write ( *, '(a)' ) 'PCHFD - Fatal error!'
            write ( *, '(a)' ) '  Could not bracket the data point.'
            return
          end if
!
!  Reset J.  This will be the new J_FIRST.
!
          j = j_new
!
!  Now find out how far to back up in the X array.
!
          do i = 1, ir-1
            if ( xe(j) < x(i) ) then
              exit
            end if
          end do
!
!  At this point, either XE(J) < X(1) or X(i-1) <= XE(J) < X(I) .
!
!  Reset IR, recognizing that it will be incremented before cycling.
!
          ir = max ( 1, i-1 )

        end if

      end if

      j_first = j

    end if

    ir = ir + 1

    if ( n < ir ) then
      exit
    end if

  end do

  return
end subroutine


subroutine pchfe ( n, x, f, d, incfd, skip, ne, xe, fe, ierr )

!*****************************************************************************80
!
!! PCHFE evaluates a piecewise cubic Hermite function at an array of points.
!
!  Description:
!
!    PCHFE may be used by itself for Hermite interpolation, or as an
!    evaluator for PCHIM or PCHIC.
!
!    PCHFE evaluates the cubic Hermite function at the points XE.
!
!    To provide compatibility with PCHIM and PCHIC, the routine includes an
!    increment between successive values of the F and D arrays.
!
!    Most of the coding between the call to CHFEV and the end of
!    the IR loop could be eliminated if it were permissible to
!    assume that XE is ordered relative to X.
!
!    CHFEV does not assume that X1 is less than X2.  Thus, it would
!    be possible to write a version of pchfe that assumes a strictly
!    decreasing X array by simply running the IR loop backwards
!    and reversing the order of appropriate tests.
!
!    The present code has a minor bug, which I have decided is not
!    worth the effort that would be required to fix it.
!    If XE contains points in [X(N-1),X(N)], followed by points less than
!    X(N-1), followed by points greater than X(N), the extrapolation points
!    will be counted (at least) twice in the total returned in IERR.
!
!    The evaluation will be most efficient if the elements of XE are
!    increasing relative to X; that is, for all J <= K,
!      X(I) <= XE(J)
!    implies
!      X(I) <= XE(K).
!
!    If any of the XE are outside the interval [X(1),X(N)],
!    values are extrapolated from the nearest extreme cubic,
!    and a warning error is returned.
!
!  Modified:
!
!    12 August 2005
!
!  Author:
!
!    Fred Fritsch,
!    Mathematics and Statistics Division,
!    Lawrence Livermore National Laboratory.
!
!    FORTRAN90 translation by John Burkardt.
!
!  Reference:
!
!    Fred Fritsch, Ralph Carlson,
!    Monotone Piecewise Cubic Interpolation,
!    SIAM Journal on Numerical Analysis,
!    Volume 17, Number 2, April 1980, pages 238-246.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of data points.
!    N must be at least 2.
!
!    Input, real ( kind = 8 ) X(N), the strictly increasing independent
!    variable values.
!
!    Input, real ( kind = 8 ) F(INCFD,N), the function values.
!    F(1+(I-1)*INCFD) is the value corresponding to X(I).
!
!    Input, real ( kind = 8 ) D(INCFD,N), the derivative values.
!    D(1+(I-1)*INCFD) is the value corresponding to X(I).
!
!    Input, integer ( kind = 4 ) INCFD, increment between successive values
!    in F and D.
!
!    Input/output, logical SKIP, controls whether data validity checks
!    should be made.  Setting the input value to FALSE will skip the checks.
!    On output with 0 <= IERR, SKIP will be set to TRUE.
!
!    Input, integer ( kind = 4 ) NE, the number of evaluation points.
!
!    Input, real ( kind = 8 ) XE(NE), points at which the function is to
!    be evaluated.
!
!    Output, real ( kind = 8 ) FE(NE), the values of the cubic Hermite
!    function at XE.
!
!    Output, integer ( kind = 4 ) IERR, error flag.
!    0, no errors.
!    positive, means that extrapolation was performed at IERR points.
!    -1, if N < 2.
!    -2, if INCFD < 1.
!    -3, if the X array is not strictly increasing.
!    -4, if NE < 1.
!    -5, error in CHFEV.
!


  integer ( kind = 4 ) incfd
  integer ( kind = 4 ) n
  integer ( kind = 4 ) ne

  real ( kind = 8 ) d(incfd,n)
  real ( kind = 8 ) f(incfd,n)
  real ( kind = 8 ) fe(ne)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierc
  integer ( kind = 4 ) ierr
  integer ( kind = 4 ) ir
  integer ( kind = 4 ) j
  integer ( kind = 4 ) j_first
  integer ( kind = 4 ) j_new
  integer ( kind = 4 ) j_save
  integer ( kind = 4 ) next(2)
  integer ( kind = 4 ) nj
  logical skip
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) xe(ne)
!
!  Check arguments.
!
  if ( .not. skip ) then

    if ( n < 2 ) then
      ierr = -1
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'PCHFE - Fatal error!'
      write ( *, '(a)' ) '  Number of data points less than 2.'
      return
    end if

    if ( incfd < 1 ) then
      ierr = -2
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'PCHFE - Fatal error!'
      write ( *, '(a)' ) '  Increment less than 1.'
      return
    end if

    do i = 2, n
      if ( x(i) <= x(i-1) ) then
        ierr = -3
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'PCHFE - Fatal error!'
        write ( *, '(a)' ) '  X array not strictly increasing.'
        return
      end if
    end do

  end if

  if ( ne < 1 ) then
    ierr = -4
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PCHFE - Fatal error!'
    write ( *, '(a)' ) '  Number of evaluation points less than 1.'
    return
  end if

  ierr = 0
  skip = .true.
!
!  Loop over intervals.
!  The interval index is IL = IR-1.
!  The interval is X(IL) <= X < X(IR).
!
  j_first = 1
  ir = 2

  do
!
!  Skip out of the loop if have processed all evaluation points.
!
    if ( ne < j_first ) then
      exit
    end if
!
!  Locate all points in the interval.
!
    j_save = ne + 1

    do j = j_first, ne
      if ( x(ir) <= xe(j) ) then
        j_save = j
        if ( ir == n ) then
          j_save = ne + 1
        end if
        exit
      end if
    end do
!
!  Have located first point beyond interval.
!
    j = j_save

    nj = j - j_first
!
!  Skip evaluation if no points in interval.
!
    if ( nj /= 0 ) then
!
!  Evaluate cubic at XE(J_FIRST:J-1).
!
      call chfev ( x(ir-1), x(ir), f(1,ir-1), f(1,ir), d(1,ir-1), d(1,ir), &
        nj, xe(j_first), fe(j_first), next, ierc )

      if ( ierc < 0 ) then
        ierr = -5
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'PCHFE - Fatal error!'
        write ( *, '(a)' ) '  Error return from CHFEV.'
        return
      end if
!
!  In the current set of XE points, there are NEXT(2) to the right of X(IR).
!
      if ( next(2) /= 0 ) then

        if ( ir < n ) then
          ierr = -5
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'PCHFE - Fatal error!'
          write ( *, '(a)' ) '  IR < N.'
          return
        end if
!
!  These are actually extrapolation points.
!
        ierr = ierr + next(2)

      end if
!
!  In the current set of XE points, there are NEXT(1) to the left of X(IR-1).
!
      if ( next(1) /= 0 ) then
!
!  These are actually extrapolation points.
!
        if ( ir <= 2 ) then
          ierr = ierr + next(1)
        else

          j_new = -1

          do i = j_first, j-1
            if ( xe(i) < x(ir-1) ) then
              j_new = i
              exit
            end if
          end do

          if ( j_new == -1 ) then
            ierr = -5
            write ( *, '(a)' ) ' '
            write ( *, '(a)' ) 'PCHFE - Fatal error!'
            write ( *, '(a)' ) '  Could not bracket the data point.'
            return
          end if
!
!  Reset J.  This will be the new J_FIRST.
!
          j = j_new
!
!  Now find out how far to back up in the X array.
!
          do i = 1, ir-1
            if ( xe(j) < x(i) ) then
              exit
            end if
          end do
!
!  At this point, either XE(J) < X(1) or X(i-1) <= XE(J) < X(I) .
!
!  Reset IR, recognizing that it will be incremented before cycling.
!
          ir = max ( 1, i-1 )

        end if

      end if

      j_first = j

    end if

    ir = ir + 1

    if ( n < ir ) then
      exit
    end if

  end do

  return
end subroutine


subroutine pchic ( ic, vc, switch, n, x, f, d, incfd, wk, nwk, ierr )

!*****************************************************************************80
!
!! PCHIC sets derivatives for a piecewise monotone cubic Hermite interpolant.
!
!  Description:
!
!    PCHIC sets derivatives needed to determine a piecewise monotone
!    piecewise cubic Hermite interpolant to the data given in X and F
!    satisfying the boundary conditions specified by IC and VC.
!
!    The treatment of points where monotonicity switches direction is
!    controlled by argument switch.
!
!    To facilitate two-dimensional applications, includes an increment
!    between successive values of the F and D arrays.
!
!    The resulting piecewise cubic Hermite function may be evaluated
!    by PCHFE or PCHFD.
!
!    User control is available over boundary conditions and
!    treatment of points where monotonicity switches direction.
!
!  Author:
!
!    Fred Fritsch,
!    Mathematics and Statistics Division,
!    Lawrence Livermore National Laboratory.
!
!    FORTRAN90 translation by John Burkardt.
!
!  Reference:
!
!    Fred Fritsch, Ralph Carlson,
!    Monotone Piecewise Cubic Interpolation,
!    SIAM Journal on Numerical Analysis,
!    Volume 17, Number 2, April 1980, pages 238-246.
!
!    Fred Fritsch, Judy Butland,
!    A Method for Constructing Local Monotone Piecewise
!    Cubic Interpolants,
!    SIAM Journal on Scientific and Statistical Computing,
!    Volume 5, Number 2, 1984, pages 300-304.
!
!    Fred Fritsch,
!    Piecewise Cubic Hermite Interpolation Package, Final Specifications,
!    Lawrence Livermore National Laboratory,
!    Computer Documentation UCID-30194, August 1982.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) IC(2), specifies desired boundary conditions:
!    IC(1) = IBEG, desired condition at beginning of data.
!    IC(2) = IEND, desired condition at end of data.
!
!    IBEG = 0  for the default boundary condition (the same as used by PCHIM).
!    If IBEG/=0, then its sign indicates whether the boundary derivative is
!    to be adjusted, if necessary, to be compatible with monotonicity:
!    0 < IBEG, if no adjustment is to be performed.
!    IBEG < 0, if the derivative is to be adjusted for monotonicity.
!
!    Allowable values for the magnitude of IBEG are:
!    1, if first derivative at x(1) is given in VC(1).
!    2, if second derivative at x(1) is given in VC(1).
!    3, to use the 3-point difference formula for D(1).
!       This reverts to the default boundary condition if N < 3.
!    4, to use the 4-point difference formula for D(1).
!       This reverts to the default boundary condition if N < 4.
!    5, to set D(1) so that the second derivative is continuous at X(2).
!       This reverts to the default boundary condition if N < 4.
!       This option is somewhat analogous to the "not a knot"
!       boundary condition provided by PCHSP.
!
!    An error return is taken if 5 < |IBEG|.
!    Only in case IBEG <= 0 is it guaranteed that the interpolant will be
!    monotonic in the first interval.  If the returned value of D(1) lies
!    between zero and 3 * SLOPE(1), the interpolant will be monotonic.  This
!    is not checked if 0 < IBEG.
!    If IBEG < 0 and D(1) had to be changed to achieve monotonicity, a
!    warning error is returned.
!
!    IEND may take on the same values as IBEG, but applied to the derivative
!    at X(N).  In case IEND = 1 or 2, the value is given in VC(2).
!
!    An error return is taken if 5 < |IEND|.
!    Only in case IEND <= 0  is it guaranteed that the interpolant will be
!    monotonic in the last interval.  If the returned value of
!    D(1+(N-1)*INCFD) lies between zero and 3 * SLOPE(N-1), the interpolant
!    will be monotonic.  This is not checked if 0 < IEND.
!    If IEND < 0 and D(1+(N-1)*INCFD) had to be changed to achieve
!    monotonicity, a warning error is returned.
!
!    Input, real ( kind = 8 ) VC(2), specifies desired boundary values,
!    as indicated above.
!    VC(1) need be set only if IC(1) = 1 or 2.
!    VC(2) need be set only if IC(2) = 1 or 2.
!
!    Input, integer ( kind = 4 ) SWITCH, indicates the desired treatment of
!    points where the direction of monotonicity switches:
!    * set SWITCH to zero if the interpolant is required to be monotonic in
!      each interval, regardless of monotonicity of data.  This will cause D
!      to be set to zero at all switch points, thus forcing extrema there.
!      The result of using this option with the default boundary conditions
!      will be identical to using PCHIM, but will generally cost more
!      compute time.  This option is provided only to facilitate comparison
!      of different switch and/or boundary conditions.
!    * set SWITCH nonzero to use a formula based on the 3-point difference
!      formula in the vicinity of switch points.  If SWITCH is positive, the
!      interpolant on each interval containing an extremum is controlled to
!      not deviate from the data by more than SWITCH * DFLOC, where DFLOC is the
!      maximum of the change of F on this interval and its two immediate
!      neighbors.  If SWITCH is negative, no such control is to be imposed.
!
!    Input, integer ( kind = 4 ) N, the number of data points.
!    N must be at least 2.
!
!    Input, real ( kind = 8 ) X(N), the strictly increasing independent
!    variable values.
!
!    Input, real ( kind = 8 ) F(INCFD,N), the dependent values to be
!    interpolated.  F(1+(I-1)*INCFD) is the value corresponding to X(I).
!
!    Output, real ( kind = 8 ) D(INCFD,N), the derivative values at the data
!    points.  These values will determine a monotone cubic Hermite function
!    on each subinterval on which the data are monotonic, except possibly
!    adjacent to switches in monotonicity.  The value corresponding to X(I)
!    is stored in D(1+(I-1)*INCFD).
!
!    Input, integer ( kind = 4 ) INCFD, increment between successive values
!    in F and D.
!
!    Workspace, real ( kind = 8 ) WK(NWK).  The user may wish to know that
!    the returned values, for I = 1 to N-1, are:
!      WK(I)     = H(I)     = X(I+1) - X(I)
!      WK(N-1+I) = SLOPE(I) = ( F(1,I+1) - F(1,I)) / H(I)
!
!    Input, integer ( kind = 4 ) NWK, the length of the work array, which must
!    be at least 2 * ( N - 1 ).
!
!    Output, integer ( kind = 4 ) IERR, error flag.
!    0, no errors.
!    1, if IBEG < 0 and D(1) had to be adjusted for monotonicity.
!    2, if IEND < 0 and D(1+(N-1)*INCFD) had to be adjusted for monotonicity.
!    3, if both of the above are true.
!    -1, if N < 2.
!    -2, if INCFD < 1.
!    -3, if the X array is not strictly increasing.
!    -4, if 5 < abs ( IBEG ).
!    -5, if 5 < abs ( IEND ).
!    -6, if both of the above are true.
!    -7, if NWK < 2 * ( N - 1 ).
!


  integer ( kind = 4 ) incfd
  integer ( kind = 4 ) n
  integer ( kind = 4 ) nwk

  real ( kind = 8 ) d(incfd,n)
  real ( kind = 8 ) f(incfd,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ibeg
  integer ( kind = 4 ) ic(2)
  integer ( kind = 4 ) iend
  integer ( kind = 4 ) ierr
  integer ( kind = 4 ) nless1
  real ( kind = 8 ) switch
  real ( kind = 8 ) vc(2)
  real ( kind = 8 ) wk(nwk)
  real ( kind = 8 ) x(n)

  if ( n < 2 ) then
    ierr = -1
    !call xerror ('pchic -- number of data points less than two', ierr, 1)
    return
  end if

  if ( incfd < 1 ) then
    ierr = -2
    !call xerror ('pchic -- increment less than one', ierr, 1)
    return
  end if

  do i = 2, n
    if ( x(i) <= x(i-1) ) then
      ierr = -3
      !call xerror ('pchic -- x-array not strictly increasing', ierr, 1)
      return
    end if
  end do

  ibeg = ic(1)
  iend = ic(2)
  ierr = 0

  if ( 5 < abs ( ibeg ) ) then
    ierr = ierr - 1
  end if

  if ( 5 < abs ( iend ) ) then
    ierr = ierr - 2
  end if

  if ( ierr < 0 ) then
    ierr = ierr - 3
   ! call xerror ('pchic -- ic out of range', ierr, 1)
    return
  end if
!
!  Function definition is ok -- go on.
!
  nless1 = n - 1
  if ( nwk < 2 * nless1 ) then
    ierr = -7
   ! call xerror ('pchic -- work array too small', ierr, 1)
    return
  end if
!
!  Set up H and slope arrays.
!
  do i = 1, nless1
    wk(i) = x(i+1) - x(i)
    wk(nless1+i) = (f(1,i+1) - f(1,i)) / wk(i)
  end do
!
!  Special case n=2 -- use linear interpolation.
!
  if ( 1 < nless1 ) then
    go to 1000
  end if

  d(1,1) = wk(2)
  d(1,n) = wk(2)
  go to 3000
!
!  Normal case  (3 <= N) .
!
 1000 continue
!
!  Set interior derivatives and default end conditions.
!
  call pchci ( n, wk(1), wk(n), d, incfd )
!
!  Set derivatives at points where monotonicity switches direction.
!
  if ( switch /= 0.0D+00 ) then

    call pchcs (switch, n, wk(1), wk(n), d, incfd, ierr)

    if ( ierr /= 0 ) then
      ierr = -8
      !call xerror ('pchic -- error return from pchcs', ierr, 1)
      return
    end if

  end if
!
!  Set end conditions.
!
 3000 continue

  if ( ibeg == 0 .and. iend == 0 ) then
    return
  end if

  call pchce ( ic, vc, n, x, wk(1), wk(n), d, incfd, ierr )

  if ( ierr < 0 ) then
    ierr = -9
    !call xerror ('pchic -- error return from pchce', ierr, 1)
    return
  end if

  return
end subroutine


subroutine pchim ( n, x, f, d, incfd, ierr )

!*****************************************************************************80
!
!! PCHIM sets derivatives for a piecewise cubic Hermite interpolant.
!
!  Discussion:
!
!    The routine set derivatives needed to determine a monotone piecewise
!    cubic Hermite interpolant to given data.
!
!    The interpolant will have an extremum at each point where
!    monotonicity switches direction.  See PCHIC if user control is desired
!    over boundary or switch conditions.
!
!    To facilitate two-dimensional applications, includes an increment
!    between successive values of the F and D arrays.
!
!    The resulting piecewise cubic Hermite function may be evaluated
!    by PCHFE or PCHFD.
!
!  Author:
!
!    Fred Fritsch,
!    Mathematics and Statistics Division,
!    Lawrence Livermore National Laboratory.
!
!    FORTRAN90 translation by John Burkardt.
!
!  Reference:
!
!    Fred Fritsch, Ralph Carlson,
!    Monotone Piecewise Cubic Interpolation,
!    SIAM Journal on Numerical Analysis,
!    Volume 17, Number 2, April 1980, pages 238-246.
!
!    Fred Fritsch, Judy Butland,
!    A Method for Constructing Local Monotone Piecewise
!    Cubic Interpolants,
!    SIAM Journal on Scientific and Statistical Computing,
!    Volume 5, Number 2, 1984, pages 300-304.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of data points.
!    N must be at least 2.
!
!    Input, real ( kind = 8 ) X(N), the strictly increasing independent
!    variable values.
!
!    Input, real ( kind = 8 ) F(INCFD,N), dependent variable values to be
!    interpolated.  F(1+(I-1)*INCFD) is the value corresponding to X(I).
!    PCHIM is designed for monotonic data, but it will work for any F-array.
!    It will force extrema at points where monotonicity switches direction.
!    If some other treatment of switch points is desired, PCHIC should be
!    used instead.
!
!    Output, real ( kind = 8 ) D(INCFD,N), the derivative values at the
!    data points.  If the data are monotonic, these values will determine
!    a monotone cubic Hermite function.  The value corresponding to X(I)
!    is stored in D(1+(I-1)*INCFD).
!
!    Input, integer ( kind = 4 ) INCFD, the increment between successive
!    values in F and D.
!
!    Output, integer ( kind = 4 ) IERR, error flag.
!    0, no errors.
!    positive, means IERR switches in the direction of monotonicity
!    were detected.
!    -1, if N < 2.
!    -2, if INCFD < 1.
!    -3, if the X array is not strictly increasing.
!


  integer ( kind = 4 ) incfd
  integer ( kind = 4 ) n

  real ( kind = 8 ) d(incfd,n)
  real ( kind = 8 ) del1
  real ( kind = 8 ) del2
  real ( kind = 8 ) dmax
  real ( kind = 8 ) dmin
  real ( kind = 8 ) drat1
  real ( kind = 8 ) drat2
  real ( kind = 8 ) dsave
  real ( kind = 8 ) f(incfd,n)
  real ( kind = 8 ) h1
  real ( kind = 8 ) h2
  real ( kind = 8 ) hsum
  real ( kind = 8 ) hsumt3
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierr
  integer ( kind = 4 ) nless1
  !real ( kind = 8 ) pchst
  real ( kind = 8 ) temp
  real ( kind = 8 ) w1
  real ( kind = 8 ) w2
  real ( kind = 8 ) x(n)
!
!  Check the arguments.
!
  if ( n < 2 ) then
    ierr = -1
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PCHIM - Fatal error!'
    write ( *, '(a)' ) '  Number of data points less than 2.'
    return
  end if

  if ( incfd < 1 ) then
    ierr = -2
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PCHIM - Fatal error!'
    write ( *, '(a)' ) '  Increment less than 1.'
    return
  end if

  do i = 2, n
    if ( x(i) <= x(i-1) ) then
      ierr = -3
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'PCHIM - Fatal error!'
      write ( *, '(a)' ) '  X array not strictly increasing.'
      return
    end if
  end do

  ierr = 0
  nless1 = n - 1
  h1 = x(2) - x(1)
  del1 = ( f(1,2) - f(1,1) ) / h1
  dsave = del1
!
!  Special case N=2, use linear interpolation.
!
  if ( n == 2 ) then
    d(1,1) = del1
    d(1,n) = del1
    return
  end if
!
!  Normal case, 3 <= N.
!
  h2 = x(3) - x(2)
  del2 = ( f(1,3) - f(1,2) ) / h2
!
!  Set D(1) via non-centered three point formula, adjusted to be
!  shape preserving.
!
  hsum = h1 + h2
  w1 = ( h1 + hsum ) / hsum
  w2 = -h1 / hsum
  d(1,1) = w1 * del1 + w2 * del2

  if ( pchst ( d(1,1), del1 ) <= 0.0D+00 ) then

    d(1,1) = 0.0D+00
!
!  Need do this check only if monotonicity switches.
!
  else if ( pchst ( del1, del2 ) < 0.0D+00 ) then

     dmax = 3.0D+00 * del1

     if ( abs ( dmax ) < abs ( d(1,1) ) ) then
       d(1,1) = dmax
     end if

  end if
!
!  Loop through interior points.
!
  do i = 2, nless1

    if ( 2 < i ) then
      h1 = h2
      h2 = x(i+1) - x(i)
      hsum = h1 + h2
      del1 = del2
      del2 = ( f(1,i+1) - f(1,i) ) / h2
    end if
!
!  Set D(I)=0 unless data are strictly monotonic.
!
    d(1,i) = 0.0D+00

    temp = pchst ( del1, del2 )

    if ( temp < 0.0D+00 ) then

      ierr = ierr + 1
      dsave = del2
!
!  Count number of changes in direction of monotonicity.
!
    else if ( temp == 0.0D+00 ) then

      if ( del2 /= 0.0D+00 ) then
        if ( pchst ( dsave, del2 ) < 0.0D+00 ) then
          ierr = ierr + 1
        end if
        dsave = del2
      end if
!
!  Use Brodlie modification of Butland formula.
!
    else

      hsumt3 = 3.0D+00 * hsum
      w1 = ( hsum + h1 ) / hsumt3
      w2 = ( hsum + h2 ) / hsumt3
      dmax = max ( abs ( del1 ), abs ( del2 ) )
      dmin = min ( abs ( del1 ), abs ( del2 ) )
      drat1 = del1 / dmax
      drat2 = del2 / dmax
      d(1,i) = dmin / ( w1 * drat1 + w2 * drat2 )

    end if

  end do
!
!  Set D(N) via non-centered three point formula, adjusted to be
!  shape preserving.
!
  w1 = -h2 / hsum
  w2 = ( h2 + hsum ) / hsum
  d(1,n) = w1 * del1 + w2 * del2

  if ( pchst ( d(1,n), del2 ) <= 0.0D+00 ) then
    d(1,n) = 0.0D+00
  else if ( pchst ( del1, del2 ) < 0.0D+00 ) then
!
!  Need do this check only if monotonicity switches.
!
    dmax = 3.0D+00 * del2

    if ( abs ( dmax ) < abs ( d(1,n) ) ) then
      d(1,n) = dmax
    end if

  end if

  return
end subroutine


subroutine pchmc ( n, x, f, d, incfd, skip, ismon, ierr )

!*****************************************************************************80
!
!! PCHMC: piecewise cubic Hermite monotonicity checker.
!
!  Discussion:
!
!    PCHMC checks a cubic Hermite function for monotonicity.
!
!    To provide compatibility with PCHIM and PCHIC, the routine includes an
!    increment between successive values of the F and D arrays.
!
!  Modified:
!
!    31 August 2002
!
!  Author:
!
!    Fred Fritsch,
!    Mathematics and Statistics Division,
!    Lawrence Livermore National Laboratory.
!
!    FORTRAN90 translation by John Burkardt.
!
!  Reference:
!
!    Fred Fritsch, Ralph Carlson,
!    Monotone Piecewise Cubic Interpolation,
!    SIAM Journal on Numerical Analysis,
!    Volume 17, Number 2, April 1980, pages 238-246.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of data points.  N must be at
!    least 2.
!
!    Input, real ( kind = 8 ) X(N), the strictly increasing independent
!    variable values.
!
!    Input, real ( kind = 8 ) F(INCFD,N), the function values.
!    F(1+(I-1)*INCFD) is the value corresponding to X(I).
!
!    Input, real ( kind = 8 ) D(INCFD,N), the derivative values.
!    D(1+(I-1)*INCFD) is the value corresponding to X(I).
!
!    Input, integer ( kind = 4 ) INCFD, increment between successive values
!    in F and D.
!
!    Input/output, logical SKIP.  On input, should be set to TRUE if the
!    user wishes to skip checks for validity of preceding parameters, or
!    to FALSE otherwise.  This will save time in case these checks have
!    already been performed.  SKIP will be set to TRUE on normal return.
!
!    Output, integer ( kind = 4 ) ISMON(N), indicates the intervals on which the
!    piecewise cubic Hermite function is monotonic.
!    For data interval [X(I),X(I+1)], and 1 <= I <= N-1, ISMON(I) is
!    -1, if function is strictly decreasing;
!     0, if function is constant;
!     1, if function is strictly increasing;
!     2, if function is non-monotonic;
!     3, if unable to determine.  This means that the D values are near the
!       boundary of the monotonicity region.  A small increase produces
!       non-monotonicity; decrease, strict monotonicity.
!    ISMON(N) indicates whether the entire function is monotonic on [X(1),X(N)].
!
!    Output, integer ( kind = 4 ) IERR, error flag.
!    0, no errors.
!    -1, if N < 2.
!    -2, if INCFD < 1.
!    -3, if the X array is not strictly increasing.
!


  integer ( kind = 4 ) incfd
  integer ( kind = 4 )  n

  !integer ( kind = 4 ) chfmc
  real ( kind = 8 ) d(incfd,n)
  real ( kind = 8 ) delta
  real ( kind = 8 ) f(incfd,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierr
  integer ( kind = 4 ) ismon(n)
  integer ( kind = 4 ) nseg
  logical skip
  real ( kind = 8 ) x(n)

  if ( .not. skip ) then

    if ( n < 2 ) then
      ierr = -1
      !call xerror ('pchmc -- number of data points less than two', ierr, 1)
      return
    end if

    if ( incfd < 1 ) then
      ierr = -2
      !call xerror ('pchmc -- increment less than one', ierr, 1)
      return
    end if

    do i = 2, n
      if ( x(i) <= x(i-1) ) then
        ierr = -3
       ! call xerror ('pchmc -- x-array not strictly increasing', ierr, 1)
        return
      end if
    end do

    skip = .true.

  end if

  nseg = n - 1

  do i = 1, nseg

     delta = ( f(1,i+1) - f(1,i) ) / ( x(i+1) - x(i) )

     ismon(i) = chfmc ( d(1,i), d(1,i+1), delta )

     if ( i == 1 ) then
        ismon(n) = ismon(1)
     else
!
!  Need to figure out cumulative monotonicity from 'multiplication table'--
!
!                    *      i s m o n (i)
!                     *  -1   0   1   2   3
!               i      *--------------------*
!               s   -1 i -1  -1   2   2   3 i
!               m    0 i -1   0   1   2   3 i
!               o    1 i  2   1   1   2   3 i
!               n    2 i  2   2   2   2   2 i
!              (n)   3 i  3   3   3   2   3 i
!                      *--------------------*
!
!  If equal or already declared nonmonotonic, no change needed.
!
        if ( ismon(i) /= ismon(n) .and. ismon(n) /= 2 ) then
           if ( 1 < max ( ismon(i), ismon(n) ) ) then
!
!  At least one is either 'no' or 'maybe'.
!
              if ( ismon(i) == 2 ) then
                 ismon(n) = 2
              else
                 ismon(n) = 3
              end if
           else if ( ismon(i) * ismon(n) < 0 ) then
!
!  Both monotonic, but in opposite senses.
!
              ismon(n) = 2
           else
!
!  At this point, one is zero, the other is +-1.
!
              ismon(n) = ismon(n) + ismon(i)
           end if
        end if
     end if

  end do

  ierr = 0

  return
end subroutine




subroutine pchsp ( ic, vc, n, x, f, d, incfd, wk, nwk, ierr )

!*****************************************************************************80
!
!! PCHSP sets derivatives needed for Hermite cubic spline interpolant.
!
!  Description:
!
!    PCHSP sets derivatives needed to determine the Hermite representation
!    of the cubic spline interpolant to given data, with specified boundary
!    conditions.
!
!    To facilitate two-dimensional applications, includes an increment
!    between successive values of the F and D arrays.
!
!    The resulting piecewise cubic Hermite function may be evaluated
!    by PCHFE or PCHFD.
!
!    This is a modified version of Carl de Boor's cubic spline routine CUBSPL.
!
!  Author:
!
!    Fred Fritsch,
!    Mathematics and Statistics Division,
!    Lawrence Livermore National Laboratory.
!
!    FORTRAN90 translation by John Burkardt.
!
!  Reference:
!
!    Carl de Boor,
!    A Practical Guide to Splines,
!    Springer-Verlag (new york, 1978).
!
!    Fred Fritsch, Ralph Carlson,
!    Monotone Piecewise Cubic Interpolation,
!    SIAM Journal on Numerical Analysis,
!    Volume 17, Number 2, April 1980, pages 238-246.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) IC(2), specifies desired boundary conditions:
!    IC(1) = IBEG, desired condition at beginning of data.
!    0, to set D(1) so that the third derivative is continuous at X(2).
!      This is the "not a knot" condition provided by de Boor's cubic spline
!      routine CUBSPL, and is the default boundary condition here.
!    1, if first derivative at X(1) is given in VC(1).
!    2, if second derivative at X(1) is given in VC(1).
!    3, to use the 3-point difference formula for D(1).
!      Reverts to the default boundary condition if N < 3.
!    4, to use the 4-point difference formula for D(1).
!      Reverts to the default boundary condition if N < 4.
!    For the "natural" boundary condition, use ibeg=2 and vc(1)=0.
!    IC(2) = IEND, desired condition at end of data.
!    IEND may take on the same values as IBEG, but applied to derivative at
!    X(N).  In case IEND = 1 or 2, the value is given in VC(2).
!
!    Input, real ( kind = 8 ) VC(2), specifies desired boundary values,
!    as indicated above.  VC(1) need be set only if IC(1) = 1 or 2.
!    VC(2) need be set only if IC(2) = 1 or 2.
!
!    Input, integer ( kind = 4 ) N, the number of data points.  N must be
!    at least 2.
!
!    Input, real ( kind = 8 ) X(N), the strictly increasing independent
!    variable values.
!
!    Input, real ( kind = 8 ) F(INCFD,N), the dependent values to be
!    interpolated.  F(1+(I-1)*INCFD) is the value corresponding to X(I).
!
!    Output, real ( kind = 8 ) D(INCFD,N), the derivative values at the
!    data points.  These values will determine the cubic spline interpolant
!    with the requested boundary conditions.  The value corresponding to
!    X(I) is stored in D(1+(I-1)*INCFD).
!
!    Input, integer ( kind = 4 ) INCFD, increment between successive values
!    in F and D.
!
!    Workspace, real ( kind = 8 ) WK(NWK).
!
!    Input, integer ( kind = 4 ) NWK, the size of WK, which must be
!    at least 2 * N.
!
!    Output, integer ( kind = 4 ) IERR, error flag.
!    0, no errors.
!    -1, if N < 2.
!    -2, if INCFD < 1.
!    -3, if the X array is not strictly increasing.
!    -4, if IBEG < 0 or 4 < IBEG.
!    -5, if IEND < 0 or 4 < IEND.
!    -6, if both of the above are true.
!    -7, if NWK is too small.
!    -8, in case of trouble solving the linear system
!        for the interior derivative values.
!


  integer ( kind = 4 ) incfd
  integer ( kind = 4 ) n

  real ( kind = 8 ) d(incfd,n)
  real ( kind = 8 ) f(incfd,n)
  real ( kind = 8 ) g
  integer ( kind = 4 ) ibeg
  integer ( kind = 4 ) ic(2)
  integer ( kind = 4 ) iend
  integer ( kind = 4 ) ierr
  integer ( kind = 4 ) index
  integer ( kind = 4 ) j
  integer ( kind = 4 ) nwk
  !real ( kind = 8 ) pchdf
  real ( kind = 8 ) stemp(3)
  real ( kind = 8 ) vc(2)
  real ( kind = 8 ) wk(2,n)
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) xtemp(4)

  if ( n < 2 ) then
    ierr = -1
    !call xerror ('pchsp -- number of data points less than two', ierr, 1)
    return
  end if

  if ( incfd < 1 ) then
    ierr = -2
    !call xerror ('pchsp -- increment less than one', ierr, 1)
    return
  end if

  do j = 2, n
    if ( x(j) <= x(j-1) ) then
      ierr = -3
      !call xerror ('pchsp -- x-array not strictly increasing', ierr, 1)
      return
    end if
  end do

  ibeg = ic(1)
  iend = ic(2)
  ierr = 0

  if ( ibeg < 0 .or. 4 < ibeg ) then
    ierr = ierr - 1
  end if

  if ( iend < 0 .or. 4 < iend ) then
    ierr = ierr - 2
  end if

  if ( ierr < 0 ) then
    go to 5004
  end if
!
!  Function definition is ok -- go on.
!
  if ( nwk < 2 * n ) then
    go to 5007
  end if
!
!  Compute first differences of X sequence and store in wk(1,.). also,
!  compute first divided difference of data and store in wk(2,.).
!
  do j = 2, n
    wk(1,j) = x(j) - x(j-1)
    wk(2,j) = ( f(1,j) - f(1,j-1) ) / wk(1,j)
  end do
!
!  Set to default boundary conditions if N is too small.
!
  if ( n < ibeg ) then
    ibeg = 0
  end if

  if ( n < iend ) then
    iend = 0
  end if
!
!  Set up for boundary conditions.
!
  if ( ibeg == 1 .or. ibeg == 2 ) then
     d(1,1) = vc(1)
  else if ( 2 < ibeg ) then
!
!  Pick up first IBEG points, in reverse order.
!
     do j = 1, ibeg
       index = ibeg - j + 1
       xtemp(j) = x(index)
       if ( j < ibeg ) then
         stemp(j) = wk(2,index)
       end if
     end do

     d(1,1) = pchdf ( ibeg, xtemp, stemp, ierr )
     if ( ierr /= 0 ) then
       go to 5009
     end if

     ibeg = 1
  end if

  if ( iend == 1 .or. iend == 2 ) then
     d(1,n) = vc(2)
  else if ( 2 < iend ) then
!
!  Pick up last IEND points.
!
     do j = 1, iend
       index = n - iend + j
       xtemp(j) = x(index)
       if ( j < iend ) then
         stemp(j) = wk(2,index+1)
       end if
     end do

     d(1,n) = pchdf ( iend, xtemp, stemp, ierr )

     if ( ierr /= 0 ) then
       go to 5009
     end if

     iend = 1

  end if
!
!  Begin coding from cubspl
!
!  A tridiagonal linear system for the unknown slopes S(1:N) of
!  F at X(1:N) is generated and then solved by Gauss elimination,
!  with s(j) ending up in d(1,j), all j.
!  wk(1,.) and wk(2,.) are used for temporary storage.
!
!  Construct first equation from first boundary condition, of the form
!    wk(2,1) * s(1) + wk(1,1) * s(2) = D(1,1)
!
  if ( ibeg == 0 ) then
!
!  No condition at left end and N = 2.
!
     if ( n == 2 ) then
        wk(2,1) = 1.0D+00
        wk(1,1) = 1.0D+00
        d(1,1) = 2.0D+00 * wk(2,2)
!
!  Not-a-knot condition at left end and 2 < N.
!
     else
        wk(2,1) = wk(1,3)
        wk(1,1) = wk(1,2) + wk(1,3)
        d(1,1) =(( wk(1,2) + 2.0D+00 * wk(1,1) ) * wk(2,2) * wk(1,3) &
                             + wk(1,2)**2 * wk(2,3) ) / wk(1,1)
     end if
  else if ( ibeg == 1 ) then
!
!  Slope prescribed at left end.
!
     wk(2,1) = 1.0D+00
     wk(1,1) = 0.0D+00
  else
!
!  Second derivative prescribed at left end.
!
     wk(2,1) = 2.0D+00
     wk(1,1) = 1.0D+00
     d(1,1) = 3.0D+00 * wk(2,2) - 0.5D+00 * wk(1,2) * d(1,1)
  end if
!
!  If there are interior knots, generate the corresponding equations and
!  carry out the forward pass of Gauss elimination, after which the J-th
!  equation reads
!
!    wk(2,j) * s(j) + wk(1,j) * s(j+1) = d(1,j).
!
  if ( 1 < n-1 ) then
    do j = 2, n-1
        if ( wk(2,j-1) == 0.0D+00 ) then
          go to 5008
        end if
        g = -wk(1,j+1) / wk(2,j-1)
        d(1,j) = g * d(1,j-1) + 3.0D+00 &
          * ( wk(1,j) * wk(2,j+1) + wk(1,j+1) * wk(2,j) )
        wk(2,j) = g * wk(1,j-1) + 2.0D+00 * ( wk(1,j) + wk(1,j+1) )
    end do
  end if
!
!  Construct last equation from second boundary condition, of the form
!
!    (-g * wk(2,n-1)) * s(n-1) + wk(2,n) * s(n) = d(1,n)
!
!  If slope is prescribed at right end, one can go directly to back-
!  substitution, since arrays happen to be set up just right for it
!  at this point.
!
  if ( iend == 1 ) then
    go to 30
  end if

  if ( iend == 0 ) then
     if ( n == 2 .and. ibeg == 0 ) then
!
!  Not-a-knot at right endpoint and at left endpoint and N = 2.
!
        d(1,2) = wk(2,2)
        go to 30
     else if ( n == 2 .or. ( n == 3 .and. ibeg == 0 ) ) then
!
!  Either ( N = 3 and not-a-knot also at left) or (N=2 and *not*
!  not-a-knot at left end point).
!
        d(1,n) = 2.0D+00 * wk(2,n)
        wk(2,n) = 1.0D+00
        if ( wk(2,n-1) == 0.0D+00 ) then
          go to 5008
        end if
        g = -1.0D+00 / wk(2,n-1)
     else
!
!  Not-a-knot and 3 <= N, and either 3 < N or also not-a-
!  knot at left end point.
!
        g = wk(1,n-1) + wk(1,n)
!
!  Do not need to check following denominators (x-differences).
!
        d(1,n) = ( ( wk(1,n) + 2.0D+00 * g ) * wk(2,n) * wk(1,n-1) &
          + wk(1,n)**2 * ( f(1,n-1) - f(1,n-2) ) / wk(1,n-1) ) / g
        if ( wk(2,n-1) == 0.0D+00 ) then
          go to 5008
        end if
        g = -g / wk(2,n-1)
        wk(2,n) = wk(1,n-1)
     end if
  else
!
!  Second derivative prescribed at right endpoint.
!
     d(1,n) = 3.0D+00 *wk(2,n) + 0.5D+00 * wk(1,n) * d(1,n)
     wk(2,n) = 2.0D+00
     if ( wk(2,n-1) == 0.0D+00 ) then
       go to 5008
     end if
     g = -1.0D+00 / wk(2,n-1)
  end if
!
!  Complete forward pass of Gauss elimination.
!
  wk(2,n) = g * wk(1,n-1) + wk(2,n)

  if ( wk(2,n) == 0.0D+00 ) then
    go to 5008
  end if

  d(1,n) = ( g * d(1,n-1) + d(1,n) ) / wk(2,n)
!
!  Carry out back substitution.
!
   30 continue

  do j = n-1, 1, -1
    if ( wk(2,j) == 0.0D+00 ) then
      go to 5008
    end if
    d(1,j) = ( d(1,j) - wk(1,j) * d(1,j+1) ) / wk(2,j)
  end do

  return
!
!  error returns.
!
 5004 continue
!
!  ic out of range return.
!
  ierr = ierr - 3
  !call xerror ('pchsp -- ic out of range', ierr, 1)
  return

 5007 continue
!
!  nwk too small return.
!
  ierr = -7
  !call xerror ('pchsp -- work array too small', ierr, 1)
  return

 5008 continue
!
!  singular system.
!  theoretically, this can only occur if successive x-values
!  are equal, which should already have been caught (ierr=-3).
!
  ierr = -8
  !call xerror ('pchsp -- singular linear system', ierr, 1)
  return
!
 5009 continue
!
!  error return from pchdf.
!  this case should never occur.
!
  ierr = -9
  !call xerror ('pchsp -- error return from pchdf', ierr, 1)

  return
end subroutine

end module


