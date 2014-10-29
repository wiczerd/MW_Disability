module V0para
!************************************************************************************************!
!This module:
!	IN: Edit to define paramters for V0main
!	OUT: parameter values and grids
!
!************************************************************************************************!


!use nrtype
!use nr_lib, only:  cheb_nodes, cheb_weight, mynrmlpdf
implicit none
save
!***Unit number***********************************!
character(LEN=10), parameter ::    sfile = 'one'	!Where to save things

!**Environmental Parameters**********************************************************************!
real(8), parameter :: 	   beta= 0.9, & 	 !People are impatient
			   R = 1/beta, &	 !People can save
			   youngD = 20.0, &	 !Length of initial young period
			   oldD = 2.0, &	 !Length of each old period
			   tlength =9, &	 !Number of periods per year	
			   Longev = 78, &	 !Median longevity	
			   xi0 = 0.16, &	 !Probability of DI accept for d=0
			   xi1 = 0.22, &	 !Probability of DI accept for d=1
	 		   xi2 = 0.50, &	 !Probability of DI accept for d=2
			   nu = 0.10, &		 !Psychic cost of applying for DI	
			   ageW = 0.02, &	 !Coefficient on Age in Mincer
			   ageW2 = -0.0001, &	 !Coefficient on Age^2 in Mincer
			   ageD = 0.1, &	 !Coefficient on Age over 45 in Disability hazard (exponential)
			   UIrr = 0.3, &	 !Replacement Rate in UI
			   DItest1 = 1.0, &	 !Earnings Index threshold 1 (See Kopecky later)
			   DItest2 = 1.5, &	 !Earnings Index threshold 2
			   DItest3 = 2.0, & 	 !Earnings Index threshold 3	
			   alfii = 0.95		 !Peristance of Alpha_i type

integer, parameter ::  	   oldN = 11,	 &	 !Number of old periods
			   TT = oldN+2		 !Total number of periods

 !Preferences----------------------------------------------------------------!
 ! u(c,p,d) = 1/(1-gam)*(c*e^(theta*d)*e^(eta*p))^(1-gam)

real(8), parameter :: 	   gam= 0.5, & 		 !IES
			   exval = 2.71828, &    !Value of constant e 	
			   eta = -0.20, &	 !Util cost of participation
			   theta = -0.22	 !Util cose of disability	
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
integer, parameter ::  nai = 2,  &	!Number of individual alpha types 
		       nbi = 2,  &	!Number of indiidual beta types
		       ndi = 2,  &	!Number of individual disability types
		       nj  = 2,  &	!Number of occupations (for now 2 types, give TFP risk)
		       nd  = 3,  &	!Number of disability extents
		       ne  = 10, &	!Points on earnings grid
		       na  = 80, &	!Points on assets grid
		       nz  = 3,  &	!Number of Occ TFP Shocks
		       maxiter = 1000   !Tolerance parameter	
		       	

real(8), parameter ::  Vtol     = 0.0001, & 	!Tolerance on V-dist
			alfi_mu  = 0.0,    & 	!Mean of alpha_i wage parameter (Log Normal)
			alfi_sig = 0.0,    & 	!Var of alpha_i wage parameter (Log Normal)
			beti_mu  = 0.0,    & 	!Mean of beta_i wage parameter (Log Normal)
			beti_sig = 0.0,    & 	!Var of beta_i wage parameter (Log Normal)
			di_lambd = 1.0,    &	!Shape of disability dist. (Exponential)
			amax = 10.0, 	   &	!Max on Asset Grid
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

			character(len=24) :: param_name
			integer:: i, j, k, unitno, t, ii 
			real(8):: summy, meps2(ne), agrid2(na), emin, emax, wtmax, step

			!***For Now, Simple Grids for Comparative Statics***!
			!Occupation Specific Transistion Multipliers: 
				!Extra probability of bad z transitions
				!Extra probability of bad d transitions
				!Comp Stat --> 1=(l,l) 2=(l,h) 3=(h,l) 4=(l,l) 
				occz(1) = 1.0
				occz(2) = 1.2
			open (newunit=unitno,file ='occtfp.txt',status ='replace')
			write (unitno,*) occz
			close (unitno)

	
			!Individual- Specific Things
				!Wage components: alpha_i, beta_i
				!Full Model-->Joint Distributed log-normal.
				!CompStat  --> 1=(h,l) 2=(h,h) 3=(l,h) 4=(l,l)
				alfi(1) = 0.5
				alfi(2) = 0.2
				beti(1) = 1.0
				beti(2) = 1.2 

				!Transition on alpha
				DO i=1,nai	!Current ai
				DO j=1,nai	!ai'
					IF (i .EQ. j) THEN
					pialf(i,j) = alfii
					ELSE 
					pialf(i,j) = 1-alfii
					EndIF
				EndDO
				EndDO

				!Extra disability risk
				dtype(1) = 1
				dtype(2) = 1.2

			!TFP 
				zgrid(1) = 0.5		!Structural Decline
				zgrid(2) = 0.9		!Recession
				zgrid(3) = 1		!Normal Times

			!Age-Specific Things
				!Wage Bonus
				wtau(1) = ageW*(25+youngD/2)+ageW2*((25+youngD/2)**2)				     !Young
				DO t=2,TT-1							
					wtau(t) = ageW*(25+youngD+t*oldD-oldD/2)+ageW2*((25+youngD+t*oldD-oldD/2)**2) !Old
				ENDDO

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
					dtau(t) = (exval**(ageD*t*oldD))	!Old (exponential)
				ENDDO		

			!Disability Extent-Specific Things
				!Wage Penalty 
				wd(1) = 0	!Healthy, no penalty
				wd(2) = -0.1	!Partially Disabled, small penalty	
				wd(3) = -0.5	!Full Disabled, large penalty

				!DI Acceptance probability
				xi(1) = xi0
				xi(2) = xi1
				xi(3) = xi2

			!Earnings Grid
				!Make linear from lowest possible wage (disabled entrant, lowest types)
				emin = alfi(3)+beti(2)*zgrid(1)+wtau(1)+wd(3)
				!... to highest, maximizing over t
				wtmax = int(min(floor(ageW/(2*ageW2)),TT-1))
				emax = alfi(1)+beti(1)*zgrid(3)+wtau(wtmax)+wd(1)
				step = (emax-emin)/(ne-1)
				DO i=1,ne
					egrid(i) = emin+step*(i-1)
				ENDdo

			!Assets Grid
				step = (amax-amin)/(na-1)
				agrid(1) = 0.0
				DO i=2,na
					agrid(i) = amin+amax/((na-i+1)**(1-gam))-1.1
				ENDdo

		!Make Markov transition matrices with all wierd invidual stuff
		!Disability: pid(id,id';i,t) <---indv. type and age specific
		DO i=1,ndi
		DO t=1,TT-1
	           pid(1,1,i,TT-t) = 1-pid1*dtau(TT-t)*dtype(i)	!Stay healthy
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


			
		     !xxxJUNKxxx!Construct general weights for gauss-chebyshev quadratures
				!weps = cheb_weight(neps)
				!DO i = 1,na
				!	meps(i) = mynrmlpdf(abar,siga,agrid(i))
				!EndDO
				!agrid = cheb_nodes(na,acheb_min,acheb_max)
				!WRITE(*,*) hgrid			
			
end subroutine setparams


					  
end module V0para
