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
real(8), parameter :: 	   beta= 0.996, & 	 !People are impatient
			   R = 1/beta, &	 !People can save
			   youngD = 20.0, &	 !Length of initial young period
			   oldD = 2.0, &	 !Length of each old period
			   tlength =12, &	 !Number of periods per year (monthly)	
			   Longev = 78, &	 !Median longevity	
			   xi0 = 0.16, &	 !Probability of DI accept for d=0
			   xi1 = 0.22, &	 !Probability of DI accept for d=1
	 		   xi2 = 0.50, &	 !Probability of DI accept for d=2
			   nu = 2.50, &		 !Psychic cost of applying for DI	
			   ageW = 0.02, &	 !Coefficient on Age in Mincer
			   ageW2 = -0.0001, &	 !Coefficient on Age^2 in Mincer
			   ageD = 0.1, &	 !Coefficient on Age over 45 in Disability hazard (exponential)
			   UIrr = 0.3, &	 !Replacement Rate in UI
			   DItest1 = 1.0, &	 !Earnings Index threshold 1 (See Kopecky later)
			   DItest2 = 1.5, &	 !Earnings Index threshold 2
			   DItest3 = 2.0, & 	 !Earnings Index threshold 3	
			   alfii = 0.95, &	 !Peristance of Alpha_i type
			   alfmu = 0.0,&	 !Mean of Alpha_i type
		  	   alfsig = 0.1,&	 !StdDev of Alpha_i type (Normal)
			   dRiskL = 0.5,&	 !Lower bound on occupation-related extra disability risk (mult factor)
			   dRiskH = 1.5, &	 !Upper bound on occupation-related extra disability risk (mult factor)
			   zRiskL = 0.5,&	 !Lower bound on occupation-related extra economic risk (mult factor)
			   zRiskH = 1.5,&	 !Upper bound on occupation-related extra economic risk (mult factor)
			   b = 0.05,&		 !Home production
			   rhho = 0.2,&		 !Probability of finding a job when long-term unemployed (David)
			   pphi = 0.2		 !Probability moving to LTU (5 months)

integer, parameter ::  	   oldN = 1,	 &	 !Number of old periods
			   TT = oldN+2		 !Total number of periods

 !Preferences----------------------------------------------------------------!
 ! u(c,p,d) = 1/(1-gam)*(c*e^(theta*d)*e^(eta*p))^(1-gam)

real(8), parameter :: 	   gam= 2.0, & 		 !IES
			   exval = 2.71828, &    !Value of constant e 	
			   eta = -0.20, &	 !Util cost of participation
			   theta = -0.22, &	 !Util cost of disability	
			   pival = 3.14159265	 !The number pi
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
			piz4 = 0.10	!Probability zh->z0

!-------------------------------------------------------------------!			

!**Programming Parameters***********************!
integer, parameter :: 	nai = 11, &		!Number of individual alpha types 
			nbi = 1,  &		!Number of indiVidual beta types
			ndi = 3,  &		!Number of individual disability types
			nj  = 1,  &		!Number of occupations (downward TFP risk variation)
			nd  = 3,  &		!Number of disability extents
			ne  = 10, &		!Points on earnings grid
			na  = 200, &		!Points on assets grid
			nz  = 3,  &		!Number of Occ TFP Shocks
			maxiter = 2, &	!Tolerance parameter	
			iaa_lowindow = 10,& 	!how far below to begin search
			iaa_hiwindow = 25	!how far above to keep searching


real(8), parameter ::   Vtol     = 0.0001, & 	!Tolerance on V-dist
!			alfi_mu  = 0.0,    & 	!Mean of alpha_i wage parameter (Log Normal)
!			alfi_sig = 0.001,    & 	!Var of alpha_i wage parameter (Log Normal)
!			beti_mu  = 0.0,    & 	!Mean of beta_i wage parameter (Log Normal)
!			beti_sig = 0.0,    & 	!Var of beta_i wage parameter (Log Normal)
!			di_lambd = 1.0,    &	!Shape of disability dist. (Exponential)
			amax 	 = 10.0,   &	!Max on Asset Grid
			amin = 0.0	   	!Min on Asset Grid
								   	

!**To build***************************!
real(8) :: 		alfi(nai), &		!Alpha_i grid- individual wage type parameter
			beti(nbi), &		!Beta_i grid- individual wage type parameter
			occz(nj), &		!Occupation-specific z risk
			wtau(TT-1), &		!Age-specific wage parameter
			wd(nd),   &		!Disability-specific wage parameter
			ptau(TT), &		!Probability of aging
			dtau(TT-1), &		!Age-related disability risk
			dtype(ndi), &		!Individual specific disability risk
			zgrid(nz), &		!TFP shock grid
			xi(nd),&		!DI acceptance probability
			agrid(na),&		!Assets grid
			egrid(ne),&		!Earnings Index Grid
			pialf(nai,nai),&	!Alpha_i transition matrix
			piz(nz,nz,nj),&		!TFP transition matrix
			pid(nd,nd,ndi,TT-1)	!Disability transition matrix

		
contains
subroutine setparams()

			logical, parameter :: lower= .FALSE. 
			integer:: i, j, k, t
			real(8):: summy, emin, emax, step, &
				  node, nodeH, nodeL, midL, midH
			real(8), parameter :: pival = 3.14159265	 !The number pi

			!***For Now, Simple Grids for Comparative Statics***!
			!Occupation Specific Transistion Multipliers: 
				!Extra probability of bad z transitions
				DO i=1,nj
					occz(i) = zRiskL +(i-1)*(zRiskH-zRiskL)/(nz-1)
				EndDO

	
			!Individual- Specific Things
				!Individual exposure to TFP shocks (beta)
				beti(1) = 1.0
				!beti(2) = 1.2 

				!Individual Wage component (alpha)- grid= 2 std. deviations
				!Nodes chosen as Gauss-Chebychev
				emin = alfmu-2*alfsig
				emax = alfmu+2*alfsig
				summy = 0
				DO i=1,nai
					k = nai-i+1
					node = COS(pival*(2.*k-1.)/(2.*nai))
					nodeL = COS(pival*(2.*max(k-1,1)-1.)/(2.*nai))
					nodeH = COS(pival*(2.*min(k+1,nai)-1.)/(2*nai))
					alfi(i) = ((node+1.)/2.)*(emax-emin)+emin
					IF (i .EQ. 1) THEN
						midH = ((nodeH-node)/2.)+node
						!pialf(:,i) = DNORDF(((midH+1)/2)*(2-2)-2)
						pialf(:,i) = alnorm((((midH+1)/2.)*(2-2)-2),lower)
						summy = summy + pialf(1,1)
					ELSEIF (i .EQ. nai) THEN
						pialf(:,i) = 1-summy	
					ELSE
						midL = node-((node-nodeL)/2.)
						midH = ((nodeH-node)/2.)+node
						pialf(:,i) = alnorm((((midH+1)/2)*(2-2)-2),lower)-alnorm((((midL+1)/2)*(2-2)-2),lower)	 
						summy = summy + pialf(1,i)
					EndIF
				EndDO


				!Pdf of alpha- N(alfmu,alfsig)
				 !Probability of landing in bin centered at node
				DO i=1,nai	!Current ai
					DO j=1,nai	!ai'
						pialf(i,j) = pialf(i,j)*(1-alfii)
					EndDO
					pialf(i,i) = pialf(i,i) + alfii	!Larger probability of staying
				EndDO

				!Extra disability risk
				DO i=1,ndi
					dtype(i) = dRiskL +dble(i-1)*(dRiskH-dRiskL)/dble(ndi)
				EndDO


			!TFP 
				zgrid(1) = 0.5		!Structural Decline
				zgrid(2) = 0.9		!Recession
				zgrid(3) = 1		!Normal Times

			!Age-Specific Things
				!Wage Bonus
				wtau(1) = ageW*(25.+youngD/2.)+ageW2*((25.+youngD/2.)**2)				     !Young
				DO t=2,TT-1							
					wtau(t) = ageW*(25.+youngD+t*oldD-oldD/2.)+ageW2*((25.+youngD+t*oldD-oldD/2.)**2) !Old
				ENDDO
				do t=1,TT-1
					wtau(t) = log(wtau(t))
				enddo

				!Aging Probability (actually, probability of not aging)
				! Mean Duration = (pr(age))^(-1)-1 <--in 1/tlength units
				ptau(1) = 1-(tlength*youngD+1)**(-1)

				DO t=2,TT-1
					ptau(t) = 1-(tlength*oldD+1)**(-1)
				ENDDO
				ptau(TT) = 1-((Longev-25+youngD+oldN*oldD)*tlength-1)**(-1)

				!Age-related disability risk
				dtau(1) = 0.5	!Young's Risk
				DO t=2,TT-1
					dtau(t) = dexp(ageD*t*oldD)	!Old (exponential)
				ENDDO		

			!Disability Extent-Specific Things
				!Wage Penalty 
				wd(1) = 0	!Healthy, no penalty
				wd(2) = -0.1	!Partially Disabled, small penalty	
				wd(3) = -0.2	!Full Disabled, large penalty

				!DI Acceptance probability
				xi(1) = xi0
				xi(2) = xi1
				xi(3) = xi2

			!Earnings Grid
				!Make linear from lowest possible wage (disabled entrant, lowest types)
				emin = dexp(beti(1)*zgrid(1)+minval(alfi)+wtau(1)+wd(nd))
				!... to highest, maximizing over t
				!wtmax = int(min(floor(ageW/(2*ageW2)),TT-1))
				emax = dexp(beti(nbi)*zgrid(nz)+maxval(alfi)+wtau(TT-1)+wd(1))
				step = (emax-emin)/dble(ne-1)
				DO i=1,ne
					egrid(i) = emin+step*dble(i-1)
				ENDdo

			!Assets Grid
				DO i=1,na
				 agrid(i)=dble(i-1)/dble(na-1)
				 agrid(i)=agrid(i)**2*(amax-amin)+amin
				ENDdo

		!Make Markov transition matrices with all wierd invidual stuff
		!Disability: pid(id,id';i,t) <---indv. type and age specific
		DO i=1,ndi
		DO t=1,TT-1
	           pid(1,1,i,TT-t) = 1-pid1*dtau(TT-t)*dtype(i)		!Stay healthy
		   pid(1,2,i,TT-t) = pid1*dtau(TT-t)*dtype(i)		!Partial Disability
		   pid(1,3,i,TT-t) = 0					!Full Disability
		   pid(2,1,i,TT-t) = 0					!Monotone
		   pid(2,2,i,TT-t) = 1-pid2*dtau(TT-t)*dtype(i)	!Stay Partial
		   pid(2,3,i,TT-t) = pid2*dtau(TT-t)*dtype(i)		!Full Disability
		   pid(3,1,i,TT-t) = 0		!Full is absorbing State
		   pid(3,2,i,TT-t) = 0
		   pid(3,3,i,TT-t) = 1
	        EndDO
		EndDO
		!Technology: piz(iz,iz';j) <--- occupations differ in downside risk       
		DO j=1,nj
		   piz(1,1,j) = 1-piz1   	!Stay in really bad shock
		   piz(1,2,j) = piz1		!Move to low shock
		   piz(1,3,j) = 0
		   piz(2,1,j) = piz2*occz(j)	  !Move to really bad shock (occupations affect it)
		   piz(2,2,j) = 1-piz2*occz(j)-piz3 !Stay in low shock
		   piz(2,3,j) = piz3		  !Move to high shock
		   piz(3,1,j) = 0		  !Must go through low to get to really bad
		   piz(3,2,j) = piz4*occz(j)	  !Move to low shock
		   piz(3,3,j) = 1-piz4*occz(j)	  !Stay in high shock
		EndDO

	

end subroutine setparams

!****************************************************************************************************************************************!
!8888888888888888888888888888888888888888		FUNCTIONS		888888888888888888888888888888888888888888888888888888888!
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

								  
end module V0para
