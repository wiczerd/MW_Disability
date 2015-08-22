module V0para
!************************************************************************************************!
!This module:
!	IN: Edit to define paramters for V0main
!	OUT: parameter values and grids
!
!************************************************************************************************!

!use nrtype
!use nr_lib, only:  cheb_nodes, cheb_weight, mynrmlpdf
!INCLUDE 'link_fnl_shared.h'
!use DNORDF_int
implicit none
save
!***Unit number***********************************!
character(LEN=10), parameter ::    sfile = 'one'	!Where to save things

!**Environmental Parameters**********************************************************************!
real(8), parameter ::	youngD = 20.0, &	!Length of initial young period
		oldD = 10.0, &		!Length of each old period
		tlen =12., &		!Number of periods per year (monthly)	
		Longev = 78.- 25., &		!Median longevity	
		xi0 = 0.16, &		!Probability of DI accept for d=0
		xi1 = 0.22, &		!Probability of DI accept for d=1
		xi2 = 0.50, &		!Probability of DI accept for d=2
		ageW = 0.02, &		!Coefficient on Age in Mincer
		ageW2 = -0.0001, &	!Coefficient on Age^2 in Mincer
		ageD = 0.1, &		!Coefficient on Age over 45 in Disability hazard (exponential)
		UIrr = 0.4, &		!Replacement Rate in UI
		DItest1 = 1.0, &	!Earnings Index threshold 1 (See Kopecky later)
		DItest2 = 1.5, &	!Earnings Index threshold 2
		DItest3 = 2.0, & 	!Earnings Index threshold 3	
		dRiskL = 0.5,&		!Lower bound on occupation-related extra disability risk (mult factor)
		dRiskH = 1.5, &		!Upper bound on occupation-related extra disability risk (mult factor)
		zRiskL = 0.5,&		!Lower bound on occupation-related extra economic risk (mult factor)
		zRiskH = 1.5,&		!Upper bound on occupation-related extra economic risk (mult factor)
		R = dexp(0.03/tlen)	!People can save

integer, parameter ::  	   oldN = 1, &!2!Number of old periods
		TT = oldN+2		!Total number of periods, oldN periods plus young and retired

!----------------------------------------------------------------------------!
	
 !Markov Chains----------------------------------------------------!
 !Baseline transition probabilities
 !Disability
 !	  1          2          3	
 !   1 (1-pid1)     pid1        0
 !   2    0       (1-pid2)     pid2
 !   3    0          0          1          
 !Technology
 !	  l          2          3	
 !   1 (1-piz1)     piz1        0
 !   2   pz2     (1-pz2-pz3)   pz3
 !   3    0         piz4      1-piz4 

real(8), parameter :: 	pid1 = 0.005, &	!Probability d0->d1
		pid2 = 0.05, &  !Probability d1->d2
		piz1 = 0.05, &  !Probability zl->z0
		piz2 = 0.01, &  !Probability z0->zl
		piz3 = 0.30, &  !Probability z0->zh
		piz4 = 0.20	!Probability zh->z0 (average duration of expansion)

!-------------------------------------------------------------------!			

!**Programming Parameters***********************!
integer, parameter ::	nal = 3, &!11	!Number of individual alpha types 
			nbi = 1,  &		!Number of indiVidual beta types
			ndi = 2,  &!3		!Number of individual disability risk types
			nj  = 1,  &		!Number of occupations (downward TFP risk variation)
			nd  = 3,  &		!Number of disability extents
			ne  = 3, &!10		!Points on earnings grid
			na  = 50, &!200	!Points on assets grid
			nz  = 3,  &		!Number of Occ TFP Shocks
			maxiter = 2000, &!	!Tolerance parameter	
			iaa_lowindow = 5,& 	!how far below to begin search
			iaa_hiwindow = 5, &	!how far above to keep searching
			Nsim = 200, &!		!how many agents to draw
			Ndat = 5000, & 		!size of data, for estimation
			Tsim = int(tlen*Longev)+1, &	!how many periods to solve
			Nk   = TT-1+(nd-1)*2+2*2!number of regressors - each period, each health and leading, occupation dynamics


! thse relate to how we compute it, i.e. what's continuous, what's endogenous, etc. 
logical, parameter ::	del_contin = .false., &	!make delta draws take continuous values or stay on the grid
			al_contin = .false.,&	!make alpha draws continuous
			j_rand = .false. 	!randomly assign j, or let choose.


real(8), parameter ::   Vtol     = 0.0001, & 	!Tolerance on V-dist
!		beti_mu  = 0.0,    & 	!Mean of beta_i wage parameter (Log Normal)
!		beti_sig = 0.0,    & 	!Var of beta_i wage parameter (Log Normal)
!		di_lambd = 1.0,    &	!Shape of disability dist. (Exponential)
		amax 	 = 10.0,   &	!Max on Asset Grid
		amin = 0.0	   	!Min on Asset Grid
								   	

!**To build***************************!
real(8) :: 	alfgrid(nal), &		!Alpha_i grid- individual wage type parameter
		beti(nbi), &		!Beta_i grid- individual wage type parameter
		occz(nj), &		!Occupation-specific z risk
		wtau(TT-1), &		!Age-specific wage parameter
		wd(nd),   &		!Disability-specific wage parameter
		ptau(TT), &		!Probability of aging
		dtau(TT-1), &		!Proportional age-related disability risk
		delgrid(ndi), &		!Individual specific disability risk
		zgrid(nz), &		!TFP shock grid
		xi(nd),&		!DI acceptance probability
		agrid(na),&		!Assets grid
		egrid(ne),&		!Earnings Index Grid
		pialf(nal,nal),&	!Alpha_i transition matrix
		piz(nz,nz,nj),&		!TFP transition matrix
		pid(nd,nd,ndi,TT-1),&	!Disability transition matrix
		Njdist(nj)		!Fraction in each occupation
		
integer :: 	dgrid(nd), &		! just enumerate the d states
		agegrid(TT)		! the mid points of the ages

!***preferences and technologies that may change
real(8) :: 	beta= dexp(-.03/tlen),&	!People are impatient (3% annual discount rate to start)
		nu = 2.50, &		!Psychic cost of applying for DI	
		alfrho = 0.95, &	!Peristance of Alpha_i type
		alfmu = 0.0,&		!Mean of Alpha_i type
		alfsig = 0.1,&		!Unconditional StdDev of Alpha_i type (Normal)
		b = 0.05,&		!Home production income
		rhho = 0.2,&		!Probability of finding a job when long-term unemployed (David)
		pphi = 0.2, &		!Probability moving to LTU (5 months)
		prob_t(TT), &		!Probability of being in each age group to start
		prborn_t(Tsim),&	!probability of being born at each point t
		amenityscale = 1.	!scale parameter of gumbel distribution for occ choice

!Preferences----------------------------------------------------------------!
! u(c,p,d) = 1/(1-gam)*(c*e^(theta*d)*e^(eta*p))^(1-gam)

real(8) :: 	gam	= 1.0, &	!IES
		eta 	= -0.2, &	!Util cost of participation
		theta 	= -0.22		!Util cost of disability	

integer :: print_lev, verbose
		
contains
subroutine setparams()

	logical, parameter :: lower= .FALSE. 
	integer:: i, j, k, t
	real(8):: summy, emin, emax, step, &
		  node, nodeH, nodeL, midL, midH
	real(8), parameter :: pival = 4.D0*datan(1.D0) !number pi

	!***For Now, Simple Grids for Comparative Statics***!
	!Occupation Specific Transistion Multipliers: 
	!Extra probability of bad z transitions
	if(nj>1) then
		do i=1,nj
			occz(i) = zRiskL +(zRiskH-zRiskL)*dble(i-1)/dble(nj-1)
		enddo
	else
		occz(1) = 1.
	endif

	!Individual- Specific Things
	!Individual exposure to TFP shocks (beta)
	beti(1) = 1.0
	!beti(2) = 1.2 

	!Individual Wage component (alpha)- grid= 2 std. deviations
	!Nodes chosen as Gauss-Chebychev
	emin = alfmu-2*alfsig
	emax = alfmu+2*alfsig
	summy = 0
	do i=1,nal
		k = nal-i+1
		node = cos(pival*(2.*k-1.)/(2.*nal))
		nodeL = cos(pival*dble(2*max(k-1,1)-1)/dble(2*nal))
		nodeH = cos(pival*dble(2*min(k+1,nal)-1)/dble(2*nal))
		alfgrid(i) = ((node+1.)/2.)*(emax-emin)+emin
		if (i == 1) then
			midH = ((nodeH-node)/2.)+node
			!pialf(:,i) = DNORDF(((midH+1)/2)*(2-2)-2)
			pialf(:,i) = alnorm((((midH+1)/2.)*(2-2)-2),lower)
			summy = summy + pialf(1,1)
		elseif (i == nal) then
			pialf(:,i) = 1-summy	
		else
			midL = node-((node-nodeL)/2.)
			midH = ((nodeH-node)/2.)+node
			pialf(:,i) = alnorm((((midH+1)/2)*(2-2)-2),lower)-alnorm((((midL+1)/2)*(2-2)-2),lower)	 
			summy = summy + pialf(1,i)
		endif
	enddo


	!Pdf of alpha- N(alfmu,alfsig)
	 !Probability of landing in bin centered at node
	do i=1,nal	!Current ai
		do j=1,nal	!ai'
			pialf(i,j) = pialf(i,j)*(1-alfrho)
		enddo
		pialf(i,i) = pialf(i,i) + alfrho	!Larger probability of staying
	enddo

	forall(i=1:nd) dgrid(i) = i

	!Extra disability risk (uniform distributed)
	do i=1,ndi
		delgrid(i) = dRiskL +dble(i-1)*(dRiskH-dRiskL)/dble(ndi-1)
	enddo

	!TFP 
	zgrid(1) = 0.5		!Structural Decline
	zgrid(2) = 0.9		!Recession
	zgrid(3) = 1		!Normal Times

	!Age-Specific Things

	! age grid building
	agegrid(1) = youngD/2
	do t = 2,TT-1
		agegrid(t) = youngD + oldD*(t-2) + oldD/2
	enddo
	agegrid(TT) = Longev/2 - (youngD + oldD*oldN) /2 + (youngD + oldD*oldN)

	!Wage Bonus
!	wtau(1) = ageW*agegrid(1)+ageW2*agegrid(1)**2			     !Young
!	do t=2,TT-1							
!		wtau(t) = ageW*agegrid(t)+ageW2*agegrid(t)**2 !Old
!	enddo
	wtau(1) = 0.
	wtau(2) = 0.15
	if(size(wtau)>=3) wtau(3) = -0.060
	if(size(wtau)>=4) wtau(4) = -0.1

	!Aging Probability (actually, probability of not aging)
	! Mean Duration = (pr(age))^(-1)-1 <--in 1/tlen units
	ptau(1) = 1-(tlen*youngD+1)**(-1)

	do t=2,TT-1
		ptau(t) = 1-(tlen*oldD+1)**(-1)
	enddo
	ptau(TT) = 1-((Longev-youngD+oldN*oldD)*tlen-1)**(-1)

	!initial age structure
	prob_t(1) = youngD/Longev
	do t=2,TT-1
		prob_t(t) = oldD/(Longev - 25.)
	enddo
	prob_t(TT) = 1.-sum(prob_t)
	!prob of getting born
	do t=2,Tsim
		prborn_t(t) = 0.01/tlen !1% population growth per year
	enddo
	prborn_t(1) = 1. - sum(prborn_t(2:Tsim))
	
	!Age-related disability risk
	dtau(1) = 0.5	!Young's Risk
	do t=2,TT-1
		dtau(t) = dexp(ageD*t*oldD)	!Old (exponential)
	enddo

	!Disability Extent-Specific Things
	!Wage Penalty 
	wd(1) = 0	!Healthy, no penalty
	wd(2) = -0.1	!Partially Disabled, small penalty	
	wd(3) = -0.2	!Full Disabled, large penalty

	!DI Acceptance probability, convert to monthly
	xi(1) = 1.-(1.-xi0)**(1./tlen)
	xi(2) = 1.-(1.-xi1)**(1./tlen)
	xi(3) = 1.-(1.-xi2)**(1./tlen)

	!Earnings Grid
	!Make linear from lowest possible wage (disabled entrant, lowest types)
	emin = dexp(beti(1)*minval(zgrid)+minval(alfgrid)+wtau(1)+wd(nd))
	!... to highest, maximizing over t
	!wtmax = int(min(floor(ageW/(2*ageW2)),TT-1))
	emax = dexp(beti(nbi)*maxval(zgrid)+maxval(alfgrid)+wtau(TT-1)+wd(1))
	step = (emax-emin)/dble(ne-1)
	do i=1,ne
		egrid(i) = emin+step*dble(i-1)
	enddo

	!Assets Grid
	do i=1,na
		agrid(i)=dble(i-1)/dble(na-1)
		agrid(i)=agrid(i)**2*(amax-amin)+amin
	enddo

	!Make Markov transition matrices with all wierd invidual stuff
	!Disability: pid(id,id';i,t) <---indv. type and age specific
	do i=1,ndi
	do t=1,TT-1

		pid(1,2,i,TT-t) = 1.-(1.-pid1*dtau(TT-t)*delgrid(i)) &	!Partial Disability 
					& **(1./tlen)	
		pid(1,1,i,TT-t) = 1.-pid(1,2,i,TT-t)			!Stay healthy
		pid(1,3,i,TT-t) = 0.					!Full Disability
		pid(2,1,i,TT-t) = 0.					!Monotone
		pid(2,3,i,TT-t) = 1.-(1.-pid2*dtau(TT-t)*delgrid(i)) &	!Full Disability
					& **(1./tlen)	
		pid(2,2,i,TT-t) = 1.-pid(2,3,i,TT-t)			!Stay Partial
		pid(3,1,i,TT-t) = 0.					!Full is absorbing State
		pid(3,2,i,TT-t) = 0.
		pid(3,3,i,TT-t) = 1.
	enddo
	enddo
		
	!Technology: piz(iz,iz';j) <--- occupations differ in downside risk       
	do j=1,nj
		piz(1,1,j) = (1-piz1)**(1./tlen)	!Stay in really bad shock
		piz(1,2,j) = 1-piz(1,1,j)		!Move to low shock
		piz(1,3,j) = 0.
		piz(2,2,j) = (1-piz2*occz(j)-piz3) &	!Stay in low shock
				&**(1./tlen)		! 
		piz(2,3,j) = 1-(1-piz3)**(1./tlen)	!Move to good state
		piz(2,1,j) = 1-piz(2,2,j)-piz(2,3,j)	!Move to really bad shock (occupations affect it)
		piz(3,1,j) = 0.				!Must go through low to get to really bad
		piz(3,3,j) = (1-piz4*occz(j))**(1./tlen)	  !Stay in high shock
		piz(3,2,j) = 1-piz(3,3,j)		!Move to low shock
		
	EndDO
	
	! distribution across occupations
	Njdist(1) = 0.5
!	Njdist(2) = 0.1
!	Njdist(3) = 0.2
!	Njdist(4) = 0.2
	Njdist = Njdist/sum(Njdist)
	

end subroutine setparams

!****************************************************************************************************************************************!
!		FUNCTIONS		
!
!****************************************************************************************************************************************!

		
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
  real ( kind = 8 ), parameter :: rr = 0.398942280385D+00
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

    alnorm = rr * exp ( - y ) &
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



subroutine random_normal(fn_val)

! Adapted from the following Fortran 77 code
!      ALGORITHM 712, COLLECTED ALGORITHMS FROM ACM.
!      THIS WORK PUBLISHED IN TRANSACTIONS ON MATHEMATICAL SOFTWARE,
!      VOL. 18, NO. 4, DECEMBER, 1992, PP. 434-435.

!  The function random_normal() returns a normally distributed pseudo-random
!  number with zero mean and unit variance.

!  The algorithm uses the ratio of uniforms method of A.J. Kinderman
!  and J.F. Monahan augmented with quadratic bounding curves.

	REAL(8),intent(out) :: fn_val

	!     Local variables
	REAL(8)     :: s = 0.449871, t = -0.386595, a = 0.19600, b = 0.25472,           &
		    r1 = 0.27597, r2 = 0.27846, u, v, x, y, q, half = 0.5

	!     Generate P = (u,v) uniform in rectangle enclosing acceptance region

	do
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
	END do

	!     Return ratio of P's coordinates as the normal deviate
	fn_val = v/u

end subroutine random_normal

subroutine random_gumbel(fn_val)
! written by D Wiczer, based on random_exponential credited to:

! Adapted from Fortran 77 code from the book:
!     Dagpunar, J. 'Principles of random variate generation'
!     Clarendon Press, Oxford, 1988.   ISBN 0-19-852202-9

! FUNCTION GENERATES A RANDOM VARIATE IN [0,INFINITY) FROM
! A NEGATIVE EXPONENTIAL DlSTRIBUTION WlTH DENSITY PROPORTIONAL
! TO EXP(-random_exponential), USING INVERSION.

REAL(8),intent(out)  :: fn_val

!     Local variable
REAL(8)  :: r

do
  call random_number(r)
  if (r > 0.) exit
end do

fn_val = -log(-log(r) )


END subroutine random_gumbel


end module V0para
