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
		oldD = 5.0, &		!Length of each old period
		tlen =12., &		!Number of periods per year (monthly)	
		Longev = 78.- 25., &		!Median longevity	
		ageW = 0.0373076, &	!Coefficient on Age in Mincer,
		ageW2 = -0.0007414, &	!Coefficient on Age^2 in Mincer
		ageD = 0.1, &		!Coefficient on Age over 45 in Disability hazard (exponential)
		UIrr = 0.4, &		!Replacement Rate in UI
		eligY  = 0.407,&	!Fraction young who are eligable
		R = dexp(0.03/tlen)	!People can save

integer, parameter ::  	   oldN = 4,&	!4!Number of old periods
		TT = oldN+2		!Total number of periods, oldN periods plus young and retired

!----------------------------------------------------------------------------!
	
 !Transition probabilities ----------------------------------------------------!
 !Baseline transition probabilities
 !Disability
 !	  1          2          3	
 !   1 (1-pid1)     pid1        0
 !   2    0       (1-pid2)     pid2
 !   3    0          0          1          
 !Technology
 !	Blocks (PiY1, PiY2) discretize AR process given (zrho,zmu,zsig) such as:
 !  PiYY =
 !	  z1YY	   z2YY
 ! z1YY	  piz11	  (1-piz11)
 ! z2YY	(1-piz22)  piz22
 ! Block transition is at rate 1/Tblock
 ! Pi =
 ! 	  zzY1	   zzY2
 ! zzY1   PiYY 	 1/Tblock*PiYY
 ! zzY2	   0	    PiYY
! annual probabilities
real(8) ::	pid1	= 0.074, &	!Probability d0->d1
		pid2	= 0.014, &  	!Probability d1->d2
		pid0	= 0.587, &	!Probability d1->d0
		dRiskL	= 0.5,&		!Lower bound on occupation-related extra disability risk (mult factor)
		dRiskH	= 1.5		!Upper bound on occupation-related extra disability risk (mult factor)

!**Programming Parameters***********************!
integer, parameter ::	nal = 4,  &!11		!Number of individual alpha types 
			nbi = 1,  &		!Number of indiVidual beta types
			ndi = 1,  &!3		!Number of individual disability risk types
			nj  = 2,  &		!Number of occupations (downward TFP risk variation)
			nd  = 3,  &		!Number of disability extents
			ne  = 4, &!10		!Points on earnings grid
			na  = 100, &!100		!Points on assets grid
			nz  = 6,  &		!Number of Occ TFP Shocks (MUST BE multiple of 2)
			maxiter = 2000, &!2000	!Tolerance parameter	
			iaa_lowindow = 5,& 	!how far below to begin search
			iaa_hiwindow = 5, &	!how far above to keep searching
			Nsim = 5000, &!		!how many agents to draw
			Ndat = 5000, & 		!size of data, for estimation
			Tsim = int(tlen*(2010-1980)), &	!how many periods to solve for simulation
			Nk   = TT+(nd-1)*2+2	!number of regressors - each age-1, each health and leading, occupation dynamics + 1 constant


! thse relate to how we compute it, i.e. what's continuous, what's endogenous, etc. 
logical, parameter ::	del_contin = .false., &	!make delta draws take continuous values or stay on the grid
			al_contin = .false.,&	!make alpha draws continuous
			j_rand = .false. 	!randomly assign j, or let choose.


real(8), parameter ::   Vtol = 1e-5, & 	!Tolerance on V-dist
!		beti_mu  = 0.0,    & 	!Mean of beta_i wage parameter (Log Normal)
!		beti_sig = 0.0,    & 	!Var of beta_i wage parameter (Log Normal)
!		di_lambd = 1.0,    &	!Shape of disability dist. (Exponential)
		amax 	 = 10.0,   &	!Max on Asset Grid
		amin = 0.0	   	!Min on Asset Grid


!**To build***************************!
real(8) :: 	alfgrid(nal), &		!Alpha_i grid- individual wage type parameter
		beti(nbi), &		!Beta_i grid- individual wage type parameter
		wtau(TT-1), &		!Age-specific wage parameter
		wd(nd),   &		!Disability-specific wage parameter
		ptau(TT), &		!Probability of aging
		dtau(TT-1), &		!Proportional age-related disability risk
		delgrid(ndi), &		!Individual specific disability risk
		zscale(nj),&		!scales occupation TFP in second period.  
		zgrid(nz,nj), &		!TFP shock grid
		xi(nd,TT-1), &		!DI acceptance probability
		agrid(na),&		!Assets grid
		egrid(ne),&		!Earnings Index Grid
		pialf(nal,nal),&	!Alpha_i transition matrix
		piz(nz,nz),&		!TFP transition matrix
		pid(nd,nd,ndi,TT-1),&	!Disability transition matrix
		prob_age(TT), &		!Probability of being in each age group to start
		prborn_t(Tsim),&	!probability of being born at each point t
		hazborn_t(Tsim), &	!hazard of being born at each point t
		Njdist(nj)		!Fraction in each occupation
		
integer :: 	dgrid(nd), &		! just enumerate the d states
		agegrid(TT)		! the mid points of the ages

!***preferences and technologies that may change
real(8) :: 	beta= 1./R,&	!People are impatient (3% annual discount rate to start)
		nu = 0.005, &		!Psychic cost of applying for DI - proportion of potential payout	
!	Idiosyncratic income risk
		alfrho = 0.988, &	!Peristance of Alpha_i type
		alfmu = 0.0,&		!Mean of Alpha_i type
		alfsig = 0.015**0.5,&	!Unconditional StdDev of Alpha_i type (Normal)
		b = 0.05,&		!Home production income
		lrho = 0.2,&		!Probability of finding a job when long-term unemployed (David)
		srho = 0.5, &		!Probability of finding a job when short-term unemployed
		pphi = 0.2, &		!Probability moving to LTU (5 months)
		xsep = 0.015, &		!Separation probability into unemployment
!	Agregate income risk
		Tblock	= 20.,	&	!Expected time before structural change (years)
		zrho	= 0.9,	&	!Persistence of the AR process
		zmu	= 0.,	&	!Drift of the AR process, should always be 0
		zsig	= 0.015**0.5,&	!Unconditional standard deviation of AR process
!		
		amenityscale = 1.,&	!scale parameter of gumbel distribution for occ choice
		xi0Y = 0.297, &		!Probability of DI accept for d=0, young
		xi1Y = 0.427, &		!Probability of DI accept for d=1, young
		xi2Y = 0.478, &		!Probability of DI accept for d=2, young
		xi0M = 0.315, &		!Probability of DI accept for d=0, middle age
		xi1M = 0.450, &		!Probability of DI accept for d=1, middle age
		xi2M = 0.503, &		!Probability of DI accept for d=2, middle age
		xi0O = 0.315, &		!Probability of DI accept for d=0, old
		xi1O = 0.450, &		!Probability of DI accept for d=1, old
		xi2O = 0.503, &		!Probability of DI accept for d=2, old
		xizcoef = -2., &		!This targets the change from 1-30% -> 1-16% acceptance rate with z deterioration
		DItest1 = 1.0, &	!Earnings Index threshold 1 (These change based on the average wage)
		DItest2 = 1.5, &	!Earnings Index threshold 2
		DItest3 = 2.0, & 	!Earnings Index threshold 3
		smthELPM = 1.		!Smoothing for the LPM


!Preferences----------------------------------------------------------------!
! u(c,p,d) = 1/(1-gam)*(c*e^(theta*d)*e^(eta*p))^(1-gam)

real(8) :: 	gam	= 1.5, &	!IES
		eta 	= -0.185, &	!Util cost of participation
		theta 	= -0.448		!Util cost of disability	

integer :: print_lev, verbose
		
contains
subroutine setparams()

	logical, parameter :: lower= .FALSE. 
	integer:: i, j, k, t
	real(8):: summy, emin, emax, step, &
		  node, nodei,nodek,nodemidH,nodemidL, &
		  alfcondsig,alfcondsigt,alfrhot,alfsigt
		  
	real(8), parameter :: pival = 4.D0*datan(1.D0) !number pi

	real(8) :: prob_age_tsim(TT,Tsim),pop_size(Tsim),cumprnborn_t(Tsim)

	!Individual- Specific Things
	!Individual exposure to TFP shocks (beta)
	beti(1) = 1.0
	!beti(2) = 1.2 

	!Individual Wage component (alpha)- grid= 2 std. deviations
	!Nodes >2 chosen as chebychev zeros, node 1 chosen to make unemployment
	
	alfrhot = alfrho**(1./tlen)
	alfcondsig = (alfsig**2*(1-alfrho**2))**0.5
	alfsigt = (alfsig**2/tlen)**0.5
	alfcondsigt = (alfsigt**2*(1-alfrhot**2))**0.5
	call rouwenhorst(nal-1,alfmu,alfrho,alfcondsig,alfgrid(2:nal),pialf(2:nal,2:nal))
	! this is the transformation if I use the annual-frequency shocks and experience them ever 1/tlen periods
	do i=2,nal
		pialf(i,i) = pialf(i,i)**(1./tlen)
		do k=2,nal
			if(k /= i) pialf(i,k) = 1. - (1.-pialf(i,k))**(1./tlen)
		enddo
		summy = sum(pialf(i,2:nal))
		if(summy /=1 ) pialf(i,2:nal)=pialf(i,2:nal)/summy !this is just numerical error
	enddo

	! value of alfgrid(1) chosen later
	pialf(1,1) = (1-srho)
	pialf(1,2:nal) = srho/dble(nal - 1)
	pialf(2:nal,1) = xsep
	forall(i=2:nal,k=2:nal) pialf(i,k) = pialf(i,k)*(1-xsep)


	forall(i=1:nd) dgrid(i) = i

	!Extra disability risk (uniform distributed)
	if(ndi>1) then
		do i=1,ndi
			delgrid(i) = dRiskL +dble(i-1)*(dRiskH-dRiskL)/dble(ndi-1)
		enddo
	else
		delgrid(1) = 0.5*(dRiskH + dRiskL)
	endif

	! TFP
	zscale = 0. 
	do j=1,nj
		zscale(j) = zsig*(dble(j)-dble(nj-1)/2.-1.)/dble(nj-1)
	enddo
	call settfp()

	
	!Age-Specific Things

	! age grid building
	agegrid(1) = youngD/2
	do i = 2,TT-1
		agegrid(i) = youngD + oldD*(i-2) + oldD/2
	enddo
	agegrid(TT) = Longev/2 - (youngD + oldD*oldN) /2 + (youngD + oldD*oldN)

	!Wage Bonus 0.0373076,-0.0007414
	do i=1,TT-1							
		wtau(i) = ageW*agegrid(i)+ageW2*agegrid(i)**2 !Old
	enddo

	!Aging Probability (actually, probability of not aging)
	! Mean Duration = (pr(age))^(-1)-1 <--in 1/tlen units
	ptau(1) = 1-(tlen*youngD)**(-1)

	do i=2,TT-1
		ptau(i) = 1-(tlen*oldD)**(-1)
	enddo
	! rate exit retirement (only one way to go.... down)
	ptau(TT) = 1-((Longev-(youngD+oldN*oldD))*tlen)**(-1)


	! BALANCED AGES:
	!initial age structure
	prob_age(1) = youngD/Longev
	do i=2,TT-1
		prob_age(i) = oldD/Longev
	enddo
	prob_age(TT) = (Longev- (youngD+oldN*oldD))/Longev
	prob_age = prob_age/sum(prob_age) !just to be sure
	pop_size(1) = sum(prob_age(1:TT-1))
	!evolution of age-structure
	prob_age_tsim(:,1) = prob_age
	cumprnborn_t = 1.
	!prob of getting born
	do t=2,Tsim
		i=1
		prborn_t(t) = prob_age_tsim(TT-1,t-1)*(1.-ptau(TT-1)) !keeps const labor force.  Otherwise =>0.01/tlen 1% population growth per year
		prob_age_tsim(i,t) = prob_age_tsim(i,t-1)*ptau(i) + prborn_t(t)
		do i = 2,TT
			prob_age_tsim(i,t) = prob_age_tsim(i,t-1)*ptau(i) + prob_age_tsim(i-1,t-1)*(1.-ptau(i-1))
		enddo
		pop_size(t) = sum(prob_age_tsim(:,t)) !should always be 1
		!prborn_t(t) = hazborn_t(t)*cumprnborn_t(t-1) !not actually hazard yet because not have initial prob
	enddo
	prborn_t(1) =  max(1.-sum(prborn_t(2:Tsim)),0.1) ! need to have some positive mass alive when the survey starts
	cumprnborn_t(1) = 1. - prborn_t(1)
	hazborn_t(1) = prborn_t(1)
	do t=2,Tsim
		hazborn_t(t) = prborn_t(t)/cumprnborn_t(t-1)
		cumprnborn_t(t) = (1.-prborn_t(t))*cumprnborn_t(t-1)
	enddo

	!DATA AGES: 


	
	!Age-related disability risk
	dtau(1) = 0.5	!Young's Risk
	do i=2,TT-1
		dtau(i) = dexp(ageD*dble(i)*oldD)	!Old (exponential)
	enddo

	!Disability Extent-Specific Things
	!Wage Penalty 
	wd(1) = 0		!Healthy, no penalty
	wd(2) = -0.2973112	!Partially Disabled, small penalty	
	wd(3) = -0.50111	!Full Disabled, large penalty

	!DI Acceptance probability for each z status

	xi(1,1) = xi0Y
	xi(2,1) = xi1Y
	xi(3,1) = xi2Y
	xi(1,2) = xi0M
	xi(2,2) = xi1M
	xi(3,2) = xi2M
	if(TT-1>=3) then
		xi(1,2) = xi0M
		xi(2,2) = xi1M
		xi(3,2) = xi2M
	endif
	if(TT-1>=4) then
		xi(1,4) = xi0O
		xi(2,4) = xi1O
		xi(3,4) = xi2O
	endif
	if(TT-1>=5) then
		xi(1,5) = xi0O
		xi(2,5) = xi1O
		xi(3,5) = xi2O
	endif
	
	!set alfgrid(1) now that everything else is known, -1 just for good measure, so they really don't want to work
	alfgrid(1) = (beti(1)*minval(zgrid)+(alfgrid(2))+wtau(1)+wd(nd))-1.

	!Earnings Grid
	!Make linear from lowest possible wage (disabled entrant, lowest types)
	emin = dexp(beti(1)*minval(zgrid)+alfgrid(2)+wtau(1)+wd(nd))
	!... to highest, maximizing over t
	!wtmax = int(min(floor(ageW/(2*ageW2)),TT-1))
	emax = dexp(beti(nbi)*maxval(zgrid)+maxval(alfgrid)+maxval(wtau)+wd(1))
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
	do j=1,ndi
	do i=1,TT-1

		pid(1,2,j,i) = 1.-(1.-pid1*dtau(i)*delgrid(j))&	!Partial Disability 
				& **(1./tlen)	
		pid(1,1,j,i) = 1.-pid(1,2,j,i)			!Stay healthy
		pid(1,3,j,i) = 0.				!Full Disability
		pid(2,1,j,i) = 1.-(1.-pid0)**(1./tlen)		!revert
		pid(2,3,j,i) = 1.-(1.-pid2*dtau(i)*delgrid(j))&	!Full Disability
					& **(1./tlen)	
		pid(2,2,j,i) = 1.-pid(2,3,j,i)-pid(2,1,j,i)	!Stay Partial
		pid(3,1,j,i) = 0.				!Full is absorbing State
		pid(3,2,j,i) = 0.
		pid(3,3,j,i) = 1.
	enddo
	enddo
		
	! Initial distribution (just for debugging) of people across occupations
	!Njdist(1) = 0.5
	!if(nj>1) then
	do j=1,nj
		Njdist(j) = 1./dble(nj)
	enddo
	!endif
	!Njdist = Njdist/sum(Njdist)

end subroutine setparams

subroutine settfp()

	integer :: i, k,j
	logical, parameter :: lower= .FALSE. 
	real(8) ::  summy, zcondsig
	real(8) :: zrhot,zsigt, zcondsigt

	zrhot = zrho**(1./tlen)
	zsigt = (zsig**2/tlen)**(0.5)
	zcondsig = ((zsig**2)*(1.-zrho**2))**(0.5)
	zcondsigt = ((zsigt**2)*(1.-zrhot**2))**(0.5)
	!first set transition probabilities at an annual basis
	call rouwenhorst(nz/2,zmu,zrho,zcondsig,zgrid(1:nz/2,1),piz(1:nz/2,1:nz/2))
	do j=2,nj
		zgrid(1:nz/2,j) = zgrid(1:nz/2,1)
	enddo
	!adjust the mean for after the transition
	do j=1,nj
		do i=(nz/2+1),nz
			zgrid(i,j) = zscale(j) + zgrid(i-nz/2,j)
		enddo
	enddo
	!adjust transition probabilities to come only ever tlen periods
	do i=1,nz/2
		piz(i,i) = piz(i,i)**(1./tlen)
		do k=1,nz/2
			if(k /= i) piz(i,k) = 1.-(1.-piz(i,k))**(1./tlen)
		enddo
		summy = sum(piz(i,1:nz/2))
		if(summy /= 1) piz(i,1:nz/2) = piz(i,1:nz/2)/summy !this is correcting numerical error
	enddo

	piz(nz/2+1:nz,nz/2+1:nz)= piz(1:nz/2,1:nz/2)
	piz(1:nz/2,nz/2+1:nz) 	= 1./(Tblock*tlen)*piz(nz/2+1:nz,nz/2+1:nz) ! struct change
	piz(1:nz/2,1:nz/2) = (1.-1./(Tblock*tlen))*piz(1:nz/2,1:nz/2) ! make it markov
	piz(nz/2+1:nz,1:nz/2) = 1./(Tblock*tlen)*piz(nz/2+1:nz,nz/2+1:nz) ! go back
	piz(nz/2+1:nz,1+nz/2:nz) = (1.-1./(Tblock*tlen))*piz(nz/2+1:nz,1+nz/2:nz) ! go back
	
	
end subroutine settfp

!****************************************************************************************************************************************!
!		FUNCTIONS		
!
!****************************************************************************************************************************************!


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
!based on Kopecky & Suen (2010) by David Wiczer

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
	sigmaz	= sigma / (1.0-rho*rho)**0.5
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
		ZPP(2:i,2:i)	= PP(1:i-1,1:i-1)
		PP(1:i,1:i) 	= p*PP(1:i,1:i) + (1.0 - p)*PPZ(1:i,1:i) + (1.0 - p)*PZP(1:i,1:i) + p*ZPP(1:i,1:i)
		PP(2:i-1,1:i)	= PP(2:i-1,1:i)/2
	enddo
	endif
	fi	= (dble(N)-1.0)**0.5*sigmaz
	grid	= (/ (i, i=0,N-1) /)
	grid	= grid *(2.0*fi / (dble(N) - 1.0)) - fi + mu
end subroutine rouwenhorst


subroutine rand_num_closed(harvest)
	!ensures we're drawing uniform [0,1] on a closed interval
	real(8), intent(inout) :: harvest
	integer :: i
	rngbound: do
		call random_number(harvest)
		harvest = harvest*1.1D0 - 0.05D0
		if(harvest<= 1.D0 .and. harvest>=0.D0) then
			exit rngbound
		endif
	enddo rngbound
	

end subroutine rand_num_closed

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
