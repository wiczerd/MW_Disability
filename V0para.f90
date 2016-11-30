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
character(len=12) :: caselabel
character(len=10) :: callog = "callog.log"
integer           :: fcallog = 7

integer, parameter:: dp=kind(0.d0) ! double precision

!**Environmental Parameters**********************************************************************!
real(8), parameter ::	youngD = 20., &	!Length of initial young period
		oldD = 5., &		!Length of each old period
		tlen =12., &		!Number of periods per year (monthly)	
		Longev = 82.- 25.,&	!Median longevity	
		ageW = 0.0373076, &	!Coefficient on Age in Mincer,
		ageW2 = -.0007414,&	!Coefficient on Age^2 in Mincer
		ageD = 0.1, &		!Coefficient on Age over 45 in Disability hazard (exponential)
		UIrr = 0.4, &		!Replacement Rate in UI
		eligY  = 0.407,&	!Fraction young who are eligable
		R =1.,&!dexp(.02/tlen),&	!People can save at 3% (not quite the rate they'd like)
		avg_unrt = 0.055,&	!average rate of unemployment over the period.
		avg_undur = 3.,&	! average months of unemployment
		upd_zscl = 0.1,&		! rate at which to update zshift
		upd_wgtrnd = 0.05		! rate at which update wage_trend

integer, parameter :: oldN = 4,&	!4!Number of old periods
		TT = oldN+2, &		!Total number of periods, oldN periods plus young and retired
		itlen = 12		! just an integer version of tlen so I don't have to keep casting
!----------------------------------------------------------------------------!

!**Programming Parameters***********************!
integer, parameter ::	nal = 5,  &!7		!Number of individual alpha types 
			ntr = 10,  &		        !Number of occupation trend points
			ndi = 2,  &		    	!Number of individual disability risk types
			nj  = 16, &!16          !Number of occupations
			nd  = 3,  &		        !Number of disability extents
			ne  = 5, &!5	        !Points on earnings grid - should be 1 if hearnlw = .true.
			na  = 50, &!100	        !Points on assets grid
			nz  = 2,  &		        !Number of aggregate shock states
			maxiter = 2000, &		!Tolerance parameter	
			Nsim = 16000, &!1000*nj !how many agents to draw
			Ndat = 5000, &          !size of data, for estimation
			Tsim = itlen*(2010-1980), &	!how many periods to solve for simulation
			struc_brk = 20,&	    ! when does the structural break happen
			Nk   = TT+(nd-1)*2+2,&	!number of regressors - each age-1, each health and leading, occupation dynamics + 1 constant
			fread = 10
			


! thse relate to how we compute it
logical            :: al_contin  = .true.,&	!make alpha draws continuous or stay on the grid
					  z_flowrts	 = .true.,& !make zj just control flow rates and not productivity (makes the next irrelevant)
					  zj_contin	 = .false.,& !make zj draws continous
					  z_regimes	 = .false.,&!different z regimes?
					  ineligNoNu = .false.,&! do young ineligable also pay the nu cost when they are ineligable?
					  dieyoung   = .false.
					  
					  
! these relate to what's changing over the simulation/across occupation
logical           ::  del_by_occ = .true.,& !delta is fully determined by occupation, right now alternative is fully random
					  j_regimes  = .true.,& !different pref shifts
					  j_rand     = .true.,&! randomly assign j, or let choose.
					  w_strchng	 = .true.,& ! w gets fed a structural change sequence
					  demog_dat	 = .true.,& !do the demographics follow 
					  NBER_tseq  = .true.	!just feed in NBER recessions?
					  


real(8), parameter ::  amax 	 = 10.0,   &	!Max on Asset Grid
					   amin = 0.0	   	!Min on Asset Grid


!**To build***************************!
real(8) :: 	alfgrid(nal), &		!Alpha_i grid- individual wage type parameter
		trgrid(ntr), &		!Wage trend grid- individual wage type parameter
		wtau(TT-1), &		!Age-specific wage parameter
		wd(nd),   &		!Disability-specific wage parameter
		ptau(TT), &		!Probability of aging
		dtau(TT-1), &		!Proportional age-related disability risk
		delgrid(ndi), &		!Individual specific disability risk
		delwt(ndi,nj),&		!The occupation-specific probability of getting a particular delta
		delcumwt(ndi+1,nj),&!cumulative dist
		occdel(nj),&		!The occupation,specific mean delta
		zshift(nj),&		!shifts occupation TFP in second period.  
		zscale(nj),&		!scales occupation TFP relative to the aggregate shock.  
		zgrid(nz,nj), &		!TFP shock grid
!		xi(nd,TT-1), &		!DI acceptance probability
		xi_d(nd),&			!1st round DI acceptance probability
		agrid(na),&			!Assets grid
		egrid(ne),&			!Earnings Index Grid
		pialf(nal,nal),&	!Alpha_i transition matrix
		ergpialf(nal),&		!Alpha_i ergodic distribution
		piz(nz,nz),&		!TFP transition matrix
		ergpiz(nz),& !ergodic TFP distribution
		pid(nd,nd,ndi,TT-1),&	!Disability transition matrix
!
		prob_age(TT), &		!Probability of being in each age group to start
		prborn_t(Tsim),&	!probability of being born at each point t
		hazborn_t(Tsim), &	!hazard of being born at each point t
		PrD3age(TT-1), &	!Fraction of D2 at each age
		PrDeath(nd,TT),&	!Probability of death during work-life

		jshift(nj,Tsim),&!Preference shift to ensure proper proportions, 2 regimes
		wage_trend(Tsim,nj),&!trend in wages
		wage_lev(nj),&		!occupation-specific differences in wage-level

!		targets for occupations
		seprisk(nz,nj),&	!occupation-cycle specific job separation
		fndrate(nz,nj),&	!occupation-cycle specific job finding rates
		occwg_trend(Tsim,nj),& !trend in occupation wage
		occwg_lev(nj),&		!level of occupation wage
		
		occsz0(nj),&		!Fraction in each occupation
		occpr_trend(Tsim,nj)!trend in occupation choice
		
integer :: 	dgrid(nd), &		! just enumerate the d states
		agegrid(TT)		! the mid points of the ages

!***preferences and technologies that may change
real(8) :: 	beta= dexp(-.05/tlen),&	!People are impatient (5% annual discount rate to start)
		nu = 5, &		!Psychic cost of applying for DI - proportion of potential payout
		util_const = 0.,&	!Give life some value
!	Idiosyncratic income risk
		alfrho = 0.988, &	!Peristence of Alpha_i type
		alfmu = 0.0,&		!Mean of Alpha_i type
		alfsig = 0.015**0.5,&	!Unconditional StdDev of Alpha_i type (Normal)
		b = 0.05,&		!Home production income
		lrho = 0.5,&		!Discount in the probability of finding a job when long-term unemployed (David)
		srho = 0.5, &		!Probability of finding a job when short-term unemployed
		pphi = 0.2, &		!Probability moving to LTU (5 months)
		xsep = 0.015, &		!Separation probability into unemployment
!	Agregate income risk
		zrho	= 0.95,	&	!Persistence of the AR process
		zmu		= 0.,	&	!Drift of the AR process, should always be 0
		zsig	= 0.15**0.5,&	!Unconditional standard deviation of AR process
! 	Health risk grid
		dRiskL	= 0.95,&	!Lower bound on occupation-related extra disability risk (multiplicative)
		dRiskH	= 1.05,&		!Upper bound on occupation-related extra disability risk

		wmean	= 1.,&		! to set the average wage on which disability stuff is set
!		
		amenityscale = 1.,&	!scale parameter of gumbel distribution for occ choice
		vscale		 = 1.,&	!will adjust to scale the discrete choice.  

		proc_time1 = 3.6,&	!time to process an application 
		proc_time2 = 28.05,&!time to process an application in stage 2 (28.05-3.6)
		xizcoef = 0.2, &	!change in acceptance rate with z deterioration
		xiagecoef = 0.,&	!increase in vocational acceptance due to age
		voc_age	= 0.25,&	!target for increase in vocation due to age
		xi_d1shift = -0.,&	!worse hlth stage acceptance for d=1
		xi_d3shift = 0.,&	!better hlth stage acceptance for d=3
		maxwin,minwin,&		!frac limits for earnings for DI probs
		DItest1 = 1.0, &	!Earnings Index threshold 1 (These change based on the average wage)
		DItest2 = 1.5, &	!Earnings Index threshold 2
		DItest3 = 2.0, & 	!Earnings Index threshold 3
		smthELPM = 1.		!Smoothing for the LPM
!		
		
		
! some handy programming terms
integer :: 		Tblock_exp	= 2e4,	&	!Expected time before structural change (years)
			Tblock_sim = struc_brk,&		!The actual time before structural change (years)
			ialUn	= 1 ,&		! the index of alpha that signifies an exogenously displaced worker
			ialL	= 2 ,&
			tri0	= ntr/2	! the index of the trgrid that is 0

!**** calibration targets
real(8) :: emp_persist = 0.98 ,&
		emp_std = 0.01 ,&
		apprt_target = .01,&	!target for application rates (to be filled below)
		dirt_target = 0.018,&	!target for di rates
		voc_accept = 0.25,&		!fraction of admissions from vocational criteria, target 1985
		hlth_accept = 0.75		!fraction taken based on health criteria, target 1985
		


!Preferences----------------------------------------------------------------!
! u(c,p,d) = 1/(1-gam)*(c*e^(theta*d)*e^(eta*p))^(1-gam)

real(8) :: 	gam	= 1.5, &	!IES
		eta 	= -0.197, &	!Util cost of participation
		theta 	= -0.224		!Util cost of disability	

integer :: print_lev, verbose
logical :: simp_concav = .false.

real(8) ::  Vtol = 5e-6 	!Tolerance on V-dist
real(8) ::  simtol = 5.e-5_dp !tolerance on simulations

contains
subroutine setparams()

	logical, parameter :: lower= .FALSE. 
	integer:: i, j, k, t,tri
	real(8):: summy, emin, emax, step, &
		  node, nodei,nodek,nodemidH,nodemidL, &
		  alfcondsig,alfcondsigt,alfrhot,alfsigt, &
		  mean_uloss,sd_uloss
		  
	real(8), parameter :: pival = 4.D0*datan(1.D0) !number pi

	real(8) :: prob_age_tsim(TT,Tsim),pop_size(Tsim),cumprnborn_t(Tsim), age_occ_read(6,18), age_read(TT), maxADL_read(nj),avgADL, &
		& occbody_trend_read(Tsim,17), wage_trend_read(Tsim,17), wage_lev_read(nj), UE_occ_read(nz,16),EU_occ_read(nz,16),apprt_read(50,2),&
		& pid_tmp(nd,nd+1,TT-1),causal_phys_read(nj)
	

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	verbose = 3
	print_lev = 2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!	
	! Read things in:
	
	!Read in the occupation size among the young.
	open(unit=fread, file="hpShareO.csv")
	do t=1,Tsim
		read(fread,*) occbody_trend_read(t,:)
	enddo
	close(fread)
	
	!Read in the occupation finding and separation rates.
	open(unit=fread, file="UE_occ.csv")
	do t=1,nz
		read(fread,*) UE_occ_read(t,:)
	enddo
	close(fread)
	open(unit=fread, file="EU_occ.csv")
	do t=1,nz
		read(fread,*) EU_occ_read(t,:)
	enddo
	close(fread)
	

	!read initial distributions of age X occ
	open(unit= fread, file = "initial_AGE_OCC.csv")
	do i=1,6
		read(fread,*) age_occ_read(i,:)
	enddo
	close(fread)
	!read initial distribution of age 
	open(unit= fread, file = "AGES_80.csv")
	read(fread,*) age_read(:)
	close(fread)

	open(unit= fread, file = "wageTrend.csv")
	do t=1,Tsim
		read(fread,*) wage_trend_read(t,:)
	enddo
	close(fread)
	open(unit= fread, file = "wageLev.csv")
	do j=1,nj
		read(fread,*) wage_lev_read(j)
	enddo
	close(fread)
	
	!Read in the disability means by occuaption
	open(unit= fread, file="maxADL.csv")
	do j=1,nj
		read(fread, *,iostat=k) maxADL_read(j)
	enddo
	close(fread)
	!Read in the disability means by occuaption
	open(unit= fread, file="causal_phys.csv")
	do j=1,nj
		read(fread, *,iostat=k) causal_phys_read(j)
	enddo
	close(fread)
	
	!Read in the disability application rates
	open(unit= fread, file="Annual_apprt.csv")
	do t=1,50
		read(fread, *,iostat=k) apprt_read(t,:)
	enddo
	close(fread)
	
	
	
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!Earnings
	

	!Individual Wage component (alpha)- grid= 2 std. deviations
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
	alfgrid(1) = alfmu - 4*alfsig !aribtrary... meant to be small

	ergpialf = 0.
	ergpialf(1) = 0.
	ergpialf(2) = alnorm( ((alfgrid(3)+alfgrid(2))/2.-alfmu) /alfsig,.false.)
	if( nal > 3) then
		do i=3,(nal-1)
			ergpialf(i) = ( &
				&	alnorm( ((alfgrid(i+1)+alfgrid(i))/2.-alfmu) /alfsig,.false.)- &
				&	alnorm( ((alfgrid(i-1)+alfgrid(i))/2.-alfmu) /alfsig,.false.) )
		enddo
	endif
	ergpialf(nal) = 1.-sum(ergpialf(2:(nal-1)))

	! the probabilities associated with going into the alpha term that is unemployment go to zero.
	pialf(1,1) = 0.
	mean_uloss = -0.154996  
	!this comes from the SIPP, code in CVW: DTall[ lfstat_wave==1 & seam==T, mean(wavewage,na.rm = T)] - DTall[ lfstat_wave==1 & seam==T & shift(EU_wave,type="lag")==T, mean(wavewage,na.rm = T)]
	sd_uloss   = (0.5888828)**0.5
	
	pialf(1,2) = alnorm( ((alfgrid(3)+alfgrid(2))/2.-mean_uloss) /sd_uloss,.false.)
	do i= 3,(nal-1)
		pialf(1,i) = ( &
				&	alnorm( ((alfgrid(i+1)+alfgrid(i))/2.- mean_uloss ) /sd_uloss,.false.)- &
				&	alnorm( ((alfgrid(i-1)+alfgrid(i))/2.- mean_uloss ) /sd_uloss,.false.) )
	enddo
	pialf(1,nal) = 1.-sum(pialf(1,2:(nal-1) ))
	pialf(2:nal,1) = 0. !exogenous separation is not built into alpha transitions

	maxwin = exp(maxval(alfgrid)) !will overwrite in the main code
	minwin = exp(minval(alfgrid)) !will overwrite in the main code

	!~~~~~~~~~~~~~~~~~~~~~~~~~~~
	! Occupation wage component
	do i=1,nj
		occwg_lev(i) = wage_lev_read(i)
		do t=1,Tsim	
			occwg_trend(t,i) = wage_trend_read(t,i+1)
			!if( occwg_trend(t,i) <minval(alfgrid)) occwg_trend(t,i) = minval(alfgrid)
			!if( occwg_trend(t,i) >maxval(alfgrid)) occwg_trend(t,i) = maxval(alfgrid)
		enddo
	enddo
	!initialize the input to the observed
	wage_lev = 0.!occwg_lev
	wage_trend = occwg_trend

	!Wage-trend grid-setup
	trgrid(1) = minval(occwg_trend)*0.9_dp
	if(ntr>1) then
		trgrid(ntr) = maxval(occwg_trend)*1.1_dp
		do tri=2,(ntr-1)
			trgrid(tri) = dble(tri-1)*(trgrid(ntr)-trgrid(1))/dble(ntr-1) + trgrid(1)
		enddo
		! set one of the grid points to 0:
		tri0 = minloc( abs(trgrid) ,dim=1)
		trgrid(tri0) = 0._dp
	else
		trgrid(1) = 0.
		tri0 = 1
	endif



	!read these numberrs in already
	seprisk = 0.
	fndrate = 0.
	do i =1,nz
		do j=1,nj
			k = 3-i !flip, recession is 1 and expansion is 2 in markov chain
			fndrate(i,j) = UE_occ_read(k,j)
			seprisk(i,j) = EU_occ_read(k,j)
		enddo
	enddo

	
	! TFP
	zshift = 0. 
	zscale = 1. 

	call settfp()

	
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!	
	!Demographics

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

	!initial age structure
	prob_age = age_read
	!prob_age(TT) = 1. - (youngD+oldN*oldD)/Longev ! balanced fraction of retirees to begin
	!prob_age(1) = age_occ_read(1,nj+1)/100./(1.-prob_age(TT) ) ! balanced : youngD/Longev
	!do i=2,TT-1
	!	prob_age(i) = age_occ_read(i,nj+1)/100. /(1.-prob_age(TT) )
	!enddo
	
	prob_age = prob_age/sum(prob_age) !just to be sure it adds to 1
	pop_size(1) = sum(prob_age(1:TT-1))
	!evolution of age-structure
!	prob_age_tsim(:,1) = prob_age
	cumprnborn_t = 1.
	!prob of getting born
	do t=2,Tsim
		i=1
		prborn_t(t) = 1.-(1.-0.005)**(1./tlen) +  1.-ptau(TT) !1% per year
		!prborn_t(t) = prob_age_tsim(TT-1,t-1)*(1.-ptau(TT-1)) !keeps const labor force.  
		!prob_age_tsim(i,t) = prob_age_tsim(i,t-1)*ptau(i) + prborn_t(t)
		!do i = 2,TT
		!	prob_age_tsim(i,t) = prob_age_tsim(i,t-1)*ptau(i) + prob_age_tsim(i-1,t-1)*(1.-ptau(i-1))
		!enddo
		!pop_size(t) = sum(prob_age_tsim(:,t)) !should always be 1
		!prborn_t(t) = hazborn_t(t)*cumprnborn_t(t-1) !not actually hazard yet because not have initial prob
	enddo
	prborn_t(1) =  product(1.-prborn_t(2:Tsim)) ! need to have some positive mass alive when the survey starts
	cumprnborn_t(1) = 1. - prborn_t(1)
	hazborn_t(1) = prborn_t(1)
	do t=2,Tsim
		hazborn_t(t) = prborn_t(t)/cumprnborn_t(t-1)
		cumprnborn_t(t) = (1.-prborn_t(t))*cumprnborn_t(t-1)
	enddo
!~  	hazborn_t = prborn_t



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
!	occupation structure
	summy = 0.
	do j=1,nj
		occpr_trend(1,j) = occbody_trend_read( 1,j+1 )
		summy = occpr_trend(1,j) + summy
	enddo
	occpr_trend(1,:) = occpr_trend(1,:)/summy
	do t=2,Tsim
		summy = 0.
		do j=1,nj
			occpr_trend(t,j) = (occbody_trend_read( t,j+1 ) - occbody_trend_read(t-1,j+1 )*(1.-ptau(1)))/prborn_t(t)
			summy = occpr_trend(t,j) + summy
		enddo
		occpr_trend(t,:) = occpr_trend(t,:)/summy
	enddo



! Initial distribution of people across occupations
	!occsz0(1) = 0.5
	!if(nj>1) then
	do j=1,nj
		occsz0(j) = age_occ_read(6,j+1)/100.
		jshift(j,:) = 0.
	enddo
	!endif
	!ensure it sums to 1
	summy =0.
	do j=1,nj
		summy = occsz0(j) + summy
	enddo
	occsz0(j) = occsz0(j)/summy

! Disability stuff	
	forall(i=1:nd) dgrid(i) = i

	!occupation-specific factor
	avgADL = 0.
	do j=1,nj
		avgADL = maxADL_read(j)*occsz0(j) + avgADL
	enddo
	!make mean 1 for occdel:
	summy = 1.-(1.-avgADL)**(1./(tlen*youngD+tlen*dble(oldN)*oldD))
	
	do j=1,nj
		occdel(j) =  (1.-(1.- (avgADL+ causal_phys_read(j)) )**(1./(tlen*youngD+tlen*dble(oldN)*oldD))) /summy
	enddo
	!will set this in draw_del
	delwt	 = 1._dp/dble(ndi) ! initialize with equal weight
		

	!Extra disability risk (uniformly spaced)
	if( maxval(occdel) > dRiskH ) then
		dRiskH = maxval(occdel)
!		dRiskL = 1. - (dRiskH-1.)
	endif
	if(dRiskL > minval(occdel)) then
		dRiskL = minval(occdel)
!		dRiskH = 1+( 1.-dRiskL )
	endif
	if(ndi>1) then
		do i=1,ndi
			delgrid(i) = dRiskL +dble(i-1)*(dRiskH-dRiskL)/dble(ndi-1)
		enddo
	else
		delgrid(1) = 0.5*(dRiskH + dRiskL)
	endif
	if(del_by_occ .eqv. .false.) delgrid = 1.
		
	!Disability Extent-Specific Things
	!Wage Penalty 
	wd(1) = 0		!Healthy, no penalty
	wd(2) = -0.2973112	!Partially Disabled, small penalty	
	wd(3) = -0.50111	!Full Disabled, large penalty

	!DI Acceptance probability for each d,t status
	xiagecoef = voc_age
	xi_d1shift = (1491*.603+2211*0.546)/(1491+2211) - (347*.581+752*.655)/(347+752) !differences from Lahiri, Vaughn, Wixon : denial rate for hlth 1,2 vs. 3,4
	xi_d3shift = (1491*.603+2211*0.546)/(1491+2211) - .484
	
	! initialize a few starting values
	xi_d(1) = 0.001!xi_d(2)+xi_d1shift	!
	xi_d(2) = xi_d(1)-xi_d1shift !xi_d(3)-xi_d3shift 
	xi_d(3) = xi_d(2)+xi_d3shift !

	
	
	if(xizcoef == 0.) then
		xizcoef = (0.4_dp - xi_d(3))/0.5_dp !dlog((0.4_dp - xi_d(2))/(1._dp - xi_d(2)))/dlog(0.5_dp) !using d=2 for the average health of applicant and 0.5 for average wage between minwage maxwage
		xizcoef = 1.1_dp * xizcoef !just and arbitrary adjustment
	endif
	!DI applications target.  Average of 1980-1985
	apprt_target = 0.
	do t=1,50
		if(apprt_read(t,1)>=1980 .and. apprt_read(t,1)<=1985 )  apprt_target = apprt_read(t,2)/6. + apprt_target
	enddo


!!!!!!!!!!!!!!!!!!!!!!!!!!
! 	Asset/SSDI wealth things:

	!set alfgrid(1) now that everything else is known, -1 just for good measure, so they really don't want to work
	alfgrid(1) = (minval(zgrid)+(alfgrid(2))+wtau(1)+wd(nd))-1.

	!Earnings Grid
	!Make linear from lowest possible wage (disabled entrant, lowest types)
	emin = dexp(trgrid(1) +alfgrid(2)+wtau(1)+wd(nd))
	!... to highest, maximizing over t
	!wtmax = int(min(floor(ageW/(2*ageW2)),TT-1))
	emax = dexp(trgrid(ntr) +maxval(alfgrid)+maxval(wtau)+wd(1))
	step = (emax-emin)/dble(ne-1)
	do i=1,ne
		egrid(i) = emin+step*dble(i-1)
	enddo

	!Assets Grid
	do i=1,na
		agrid(i)=dble(i-1)/dble(na-1)
		agrid(i)=agrid(i)**2*(amax-amin)+amin
	enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!Make 2-year Markov transition matrices with all wierd invidual stuff
	!Disability: pid(id,id';i,t) <---indv. type and age specific
	!  					& 0 & 1 & 2 & dead\\
!    pid_tmp(1,:,1) = (/0.978 , 0.018 , 0.003 , 0.001/)
!    pid_tmp(2,:,1) = (/0.417 , 0.527 , 0.056 , 0.000/)
!    pid_tmp(3,:,1) = (/0.000 , 0.000 , 0.952 , 0.048/)
!    !\multicolumn{4}{c}{Age \quad 22-45}\\

!    pid_tmp(1,:,2) = (/0.960 , 0.024 , 0.014 , 0.002/)
!    pid_tmp(2,:,2) = (/0.305 , 0.603 , 0.071 , 0.022/)
!    pid_tmp(3,:,2) = (/0.000 , 0.000 , 0.950 , 0.050/)
    
!    pid_tmp(1,:,3) = (/0.960 , 0.024 , 0.014 , 0.002/)
!    pid_tmp(2,:,3) = (/0.305 , 0.603 , 0.071 , 0.022/)
!    pid_tmp(3,:,3) = (/0.000 , 0.000 , 0.950 , 0.050/)
!    !\multicolumn{4}{c}{Age \quad 46-55}\\

!    pid_tmp(1,:,4) =  (/0.945 , 0.038 , 0.010 , 0.008/)
!    pid_tmp(2,:,4) =  (/0.375 , 0.455 , 0.134 , 0.037/)
!    pid_tmp(3,:,4) =  (/0.000 , 0.000 , 0.985 , 0.015/)
!    !\multicolumn{4}{c}{Age \quad 56-60}\\

!    pid_tmp(1,:,5) =  (/0.896 , 0.072 , 0.017 , 0.016/)
!    pid_tmp(2,:,5) =  (/0.174 , 0.692 , 0.089 , 0.045/)
!    pid_tmp(3,:,5) =  (/0.000 , 0.000 , 0.986 , 0.014/)
    !\multicolumn{4}{c}{Age \quad 61-65}\\
	pid_tmp(1,:,1) = (/ .9562 ,  .0344  , .0084, 9.8e-04/)
	pid_tmp(2,:,1) = (/ .5528 ,  .3662  , .0801, 9.2e-04/)
	pid_tmp(3,:,1) = (/	  0.  ,    0.   , .9979,   .0021/)
	
	pid_tmp(1,:,2) = (/ .9378, .0468, .0137, .0018/)
	pid_tmp(2,:,2) = (/ .4429, .4325, .1198, .0048/)
	pid_tmp(3,:,2) = (/   0.0,   0.0, .9852, .0148/)

	pid_tmp(1,:,3) = (/ .9378, .0468, .0137, .0018/)
	pid_tmp(2,:,3) = (/ .4429, .4325, .1198, .0048/)
	pid_tmp(3,:,3) = (/   0.0,   0.0, .9852, .0148/)

	pid_tmp(1,:,4) = (/ .9293, .0518, .0159, .0029/)
	pid_tmp(2,:,4) = (/ .4545, .3929, .1338, .0188/)
	pid_tmp(3,:,4) = (/   0.0,   0.0, .9947, .0053/)
	
	pid_tmp(1,:,5) = (/ .8944, .0797, .0213, .0046/)
	pid_tmp(2,:,5) = (/ .3665, .5391, .0818, .0126/)
	pid_tmp(3,:,5) = (/   0.0,   0.0, .9701, .0299/)
	

	
	
	! make stochastic--minus that death probability
	do t =1,TT-1
		do i=1,3
			pid_tmp(i,1:3,t) = pid_tmp(i,1:3,t)/(1. - pid_tmp(i,4,t))
		enddo
	enddo
	
	! convert to monthly and multiply by delgrid (was a 2-year transition matrix)
	do i=1,TT-1
		do j=1,ndi
		pid(1,2,j,i) = (1. - ( 1.-pid_tmp(1,2,i) )**(0.50_dp/tlen)) *delgrid(j)
		pid(1,3,j,i) = (1. - ( 1.-pid_tmp(1,3,i) )**(0.50_dp/tlen)) *delgrid(j)
		pid(1,1,j,i) = 1.- pid(1,2,j,i) - pid(1,3,j,i)
		
		pid(2,1,j,i) = (1. - ( 1.-pid_tmp(2,1,i) )**(0.5_dp/tlen)) 
		pid(2,3,j,i) = (1. - ( 1.-pid_tmp(2,3,i) )**(0.5_dp/tlen)) *delgrid(j)
		pid(2,2,j,i) = 1. - pid(2,1,j,i) - pid(2,3,j,i)
		
		pid(3,3,j,i) = 1.
		pid(3,1:2,j,i) = 0.
		
		enddo
	enddo
	
	!was: PrD3age = (/0.1,0.17,0.21,0.27,0.34 /)
	PrD3age = (/0.0444_dp,0.0756_dp,0.0933_dp,0.1201_dp,0.1617_dp /)
	
	PrDeath(:,1:(TT-1)) = 1.-(1._dp-pid_tmp(:,4,:))**(0.5_dp/tlen)
	PrDeath(:,TT) = PrDeath(:,TT-1)
	
end subroutine setparams

subroutine settfp()

	integer :: i, k,j,nzblock
	logical, parameter :: lower= .FALSE. 
	real(8) ::  summy, zcondsig
	real(8) :: zrhot,zsigt, zcondsigt
	real(8), allocatable :: piblock(:,:),ergpi1(:,:),ergpi2(:,:)
	real(8) :: Zzgrid(nz)

	external dgemm

	if(z_regimes .eqv. .true.) then
		nzblock = nz/2
	else
		nzblock = nz
	endif
	
	allocate(piblock(nzblock,nzblock))
	allocate(ergpi1(nzblock,nzblock))
	allocate(ergpi2(nzblock,nzblock))
	piblock = 0.
	ergpi1  = 0.
	ergpi2  = 0.


	if( nz>2 ) then !just taking nber probabilities then

		zrhot = zrho**(1./tlen)
		zsigt = zsig**(1/tlen)
		zcondsig = ((zsig**2)*(1.-zrho**2))**(0.5)
		zcondsigt = ((zsigt**2)*(1.-zrhot**2))**(0.5)
		!first set transition probabilities at an annual basis
		call rouwenhorst(nzblock,zmu,zrho,zcondsig,zgrid(1:nzblock,1),piz(1:nzblock,1:nzblock))

		if(z_regimes .eqv. .true.) then
			!adjust the mean for after the transition
			do j=2,nj
				zgrid(1:nz/2,j) = zgrid(1:nz/2,1)
			enddo
			do j=1,nj
				do i=(nz/2+1),nz
					zgrid(i,j) = zshift(j) + zgrid(i-nz/2,j)
				enddo
			enddo
			!adjust transition probabilities to come only ever tlen periods
			do i=1,nz/2
				piz(i,i) = piz(i,i)**(1./tlen)
				do k=1,nz/2
					if(k /= i) piz(i,k) = 1.-(1.-piz(i,k))**(1./tlen)
				enddo
				summy = sum(piz(i,1:nzblock))
				if(summy /= 1) piz(i,1:nzblock) = piz(i,1:nzblock)/summy !this is correcting numerical error
			enddo

			piz(nz/2+1:nz,nz/2+1:nz)= piz(1:nz/2,1:nz/2)
			piz(1:nz/2,nz/2+1:nz) 	= 1./(dble(Tblock_exp)*tlen)*piz(nz/2+1:nz,nz/2+1:nz) ! struct change
			piz(1:nz/2,1:nz/2) = (1.-1./(dble(Tblock_exp)*tlen))*piz(1:nz/2,1:nz/2) ! make it markov
			piz(nz/2+1:nz,1:nz/2) = 1./(dble(Tblock_exp)*tlen)*piz(nz/2+1:nz,nz/2+1:nz) ! go back
			piz(nz/2+1:nz,1+nz/2:nz) = (1.-1./(dble(Tblock_exp)*tlen))*piz(nz/2+1:nz,1+nz/2:nz) ! go back

			piblock = piz(1:nz/2,1:nz/2)
			
		else
			piblock= piz
			Zzgrid = zgrid(:,1)
			do j=1,nj
				zgrid(:,j) = Zzgrid*zscale(j)
			enddo
		endif

	else
		! probability of entering recession  4/(58+12+92+120+73) = 0.011267606
		piz(2,1) = 4./(58.+12.+92.+120.+73.)
		! probability of exit from recession 4/(6+16+8+8+18) = 0.071428571
		piz(1,2) = 4./(6.+16.+8.+8.+18.)
		piz(1,1) = 1.-piz(1,2)
		piz(2,2) = 1.-piz(2,1)
		piblock = piz
		if( z_flowrts .eqv. .true.) then
			Zzgrid = 0.
			zgrid = 0.
		else
			Zzgrid  = (/ zmu-zsig, zmu+zsig/)
			do j=1,nj
					zgrid(:,j) = Zzgrid*zscale(j)
			enddo
		endif
	endif
		!DGEMM('N','N',  M,  N,    K, ALPHA,  A,     M,    B,       K,  BETA,     C,       M)
	call dgemm('n','n',nzblock,nzblock,nzblock,1._dp,piblock,nzblock, piblock, nzblock, 0._dp, ergpi1,nzblock)
	do i=1,10000
			call dgemm('n','n',nzblock,nzblock,nzblock,1._dp,piblock,nzblock, ergpi1, nzblock, 0._dp, ergpi2,nzblock)
			ergpi1 = ergpi2
	enddo
	if(z_regimes .eqv. .true.) then
		do i=1,nz/2
			ergpiz(i) = ergpi1(1,i)
			ergpiz(i+nz/2) = ergpi1(1,i)
		enddo
	else
		ergpiz = ergpi1(1,:)
	endif
	deallocate(ergpi1,ergpi2,piblock)
		
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

	real(8),intent(out) :: fn_val

	!     Local variables
	real(8)     :: s = 0.449871, t = -0.386595, a = 0.19600, b = 0.25472,           &
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

real(8),intent(out)  :: fn_val

!     Local variable
real(8)  :: r

do
  call random_number(r)
  if (r > 0.) exit
end do

fn_val = -log(-log(r) )


END subroutine random_gumbel

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
	real(8), INTENT(IN) :: ax,bx,cx
	real(8), INTENT(IN), optional :: tol_in
	real(8), INTENT(OUT) :: xmin
	integer , intent(out) :: info
	integer , intent(out), optional :: niter
	real(8) :: brent
	real(8), dimension(:), intent(in) :: funcp ! a vector of function parameters
	INTERFACE
		FUNCTION func(x, funcp)
!		use nrtype
!		USE mkl95_precision, ONLY: 8 => DP
		IMPLICIT NONE
		real(8), INTENT(IN) :: x
		real(8), INTENT(IN), dimension(:) :: funcp
		real(8) :: func
		END FUNCTION func
	END INTERFACE
	INTEGER, PARAMETER :: ITMAX=100
	real(8) :: TOL
	real(8), PARAMETER :: CGOLD=0.381966011250105_8,ZEPS=1.0e-3_8*epsilon(ax)
	INTEGER :: iter
	real(8) :: a,b,d,e,etemp,fu,fv,fw,fx,p,q,r,tol1,tol2,u,v,w,x,xm
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
		xm=0.5_8*(a+b)
		tol1=tol*abs(x)+ZEPS
		tol2=2.0_8*tol1
		if (abs(x-xm) <= (tol2-0.5_8*(b-a))) then
			xmin=x
			brent=fx
			if( (present(niter).eqv. .true.) ) niter = iter-1
			RETURN
		end if
		if (abs(e) > tol1) then
			r=(x-w)*(fx-fv)
			q=(x-v)*(fx-fw)
			p=(x-v)*q-(x-w)*r
			q=2.0_8*(q-r)
			if (q > 0.0) p=-p
			q=abs(q)
			etemp=e
			e=d
			if (abs(p) >= abs(0.5_8*q*etemp) .or. &
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
	real(8), INTENT(IN) :: x1,x2,tol
	real(8) :: zbrent
	real(8), dimension(:), intent(in) :: funcp ! a vector of function parameters
	integer , intent(out) :: flag
	INTERFACE
		FUNCTION func(x, funcp)
!		use nrtype
!		USE mkl95_precision, ONLY: 8 => DP
		IMPLICIT NONE
		real(8), INTENT(IN) :: x
		real(8), INTENT(IN), dimension(:) :: funcp
		real(8) :: func
		END FUNCTION func
	END INTERFACE
	INTEGER, PARAMETER :: ITMAX=100
	real(8), PARAMETER :: EPS=epsilon(x1)
	INTEGER :: iter
	real(8) :: a,b,c,d,e,fa,fb,fc,p,q,r,s,tol1,xm
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
		tol1=2.0_8*EPS*abs(b)+0.5_8*tol
		xm=0.5_8*(c-b)
		if (abs(xm) <= tol1 .or. fb == 0.0) then
			zbrent=b
			flag = 0
			RETURN
		end if
		if (abs(e) >= tol1 .and. abs(fa) > abs(fb)) then
			s=fb/fa
			if (a == c) then
				p=2.0_8*xm*s
				q=1.0_8-s
			else
				q=fa/fc
				r=fb/fc
				p=s*(2.0_8*xm*q*(q-r)-(b-a)*(r-1.0_8))
				q=(q-1.0_8)*(r-1.0_8)*(s-1.0_8)
			end if
			if (p > 0.0) q=-q
			p=abs(p)
			if (2.0_8*p  <  min(3.0_8*xm*q-abs(tol1*q),abs(e*q))) then
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


subroutine hptrend(t,y,phi, yt,yd)
!
!  Detrend time series with Hodrick-Prescott method:
!
!      sum  (y(t)-yt(t))^2 + phi sum [(yt(t+1)-yt(t))-(yt(t)-yt(t-1))]^2
!     t=1:T                     t=2:T-1     
!
!  for the data in y.   A larger phi results in a smoother trend series.  
!  For quarterly data Hodrick and Prescott (1980) use phi=1600.
!
!  Also returned are the series with no trend:  yd =y-yt
!
!          

!  Ellen R. McGrattan,  4-23-87
!  Revised, 5-16-06, ERM

!  References
!  ----------
!  [1] Hodrick, Robert J. and Edward C. Prescott, "Post-War U.S. Business 
!        Cycles: An Empirical Investigation," Working Paper, Carnegie-Mellon 
!        University, November, 1980.
!  [2] Prescott, Edward C., "Theory Ahead of Business Cycle Measurement," 
!        QUARTERLY REVIEW, Federal Reserve Bank of Minneapolis, Fall 1986.
!      


	implicit none
	integer, parameter             :: maxit=100,nm=40
	integer, intent(in)            :: t
	integer, dimension(t+1)        :: ip
	integer, dimension(t)          :: iwork,ju
	integer, dimension(5*t-6)      :: ia,ja,jlu
	integer                        :: i,j,info
	real(8), intent(in)               :: phi
	real(8)                           :: eps
	real(8), dimension(t),intent(in)  :: y
	real(8), dimension(t),intent(out) :: yd,yt
	real(8), dimension(t)             :: yn
	real(8), dimension(t,t)           :: a
	real(8), dimension(5*t-6)         :: s,alu
	real(8), dimension(t,nm+1)        :: v
	real(8), dimension(t-2)           :: x2
	real(8), dimension(t-3)           :: x3
	real(8), dimension(t-4)           :: x4

	external dgesv

	!if (t<301) then <- DGW comment: Ellen has a version for large and for small datasets.  Only implemented here is the small dataset version

	a              = 0.
	a(1,1:3)       = (/ 1.+phi, -2.*phi, phi /)
	a(2,1:4)       = (/ -2.*phi, 1.+5.*phi, -4.*phi, phi /)
	do i=3,t-2
		a(i,i-2:i+2) = (/ phi, -4.*phi, 1.+6.*phi, -4.*phi, phi /)
	enddo
	a(t-1,t-3:t)   = (/ phi, -4.*phi, 1.+5.*phi, -2.*phi /)
	a(t,t-2:t)     = (/ phi, -2.*phi, 1.+phi /)

	yt             = y
	call dgesv(t,1,a,t,iwork,yt,t,info)
	if (info /= 0) then
		write(*,*) 'Error: problem with dgesv'
		return
	endif

	!else

	!	s(1:3)         = (/ 1.+phi, -2.*phi, phi /)
	!	ia(1:3)        = (/ 1,1,1 /)
	!	ja(1:3)        = (/ 1,2,3 /)
	!	s(4:7)         = (/ -2.*phi, 1.+5.*phi, -4.*phi, phi /)
	!	ia(4:7)        = (/ 2,2,2,2 /)
	!	ja(4:7)        = (/ 1,2,3,4 /)
	!	j              = 8
	!	do i=3,t-2
	!		s(j:j+4)     = (/ phi, -4.*phi, 1.+6.*phi, -4.*phi, phi /)
	!		ia(j:j+4)    = (/ i,i,i,i,i /)
	!		ja(j:j+4)    = (/ i-2,i-1,i,i+1,i+2 /)
	!		j            = j+5
	!	enddo
	!	s(j:j+3)       = (/ phi, -4.*phi, 1.+5.*phi, -2.*phi /)
	!	ia(j:j+3)      = (/ t-1,t-1,t-1,t-1 /)
	!	ja(j:j+3)      = (/ t-3,t-2,t-1,t /)
	!	j              = j+4
	!	s(j:j+2)       = (/ phi, -2.*phi, 1.+phi /)
	!	ia(j:j+2)      = (/ t,t,t /)
	!	ja(j:j+2)      = (/ t-2,t-1,t /)
	!	j              = 1
	!	ip(1)          = 1
	!	do i=2,5*t-6
	!		if (ia(i) /= ia(i-1)) then
	!			j          = j+1
	!			ip(j)      = i
	!		endif
	!	enddo
	!	ip(t+1)        = 5*t-5

	!	call ilu0(t,s,ja,ip,alu,jlu,ju,iwork,info)
	!	if (info /= 0) then
	!		write(*,*) 'Error: problem with ilu0'
	!		return
	!	endif
	!	eps            = 1.e-10
	!	yn             = y
	!	call pgmres(t,nm,yn,yt,v,eps,maxit,0,s,ja,ip,alu,jlu,ju,info)
	!	if (info /= 0) then
	!		write(*,*) 'Error: problem with pgmres'
	!		return
	!	endif

	!endif
	yd  = y-yt

end subroutine hptrend


! qsort stuff

recursive subroutine quicksort(a, first, last)
! quicksort.f -*-f90-*-
! Author: t-nissie
! License: GPLv3
! Gist: https://gist.github.com/t-nissie/479f0f16966925fa29ea
	implicit none
	real*8  a(*), x, t
	integer first, last
	integer i, j

	x = a( (first+last) / 2 )
	i = first
	j = last
	do
		do while (a(i) < x)
			i=i+1
		end do
		do while (x < a(j))
			j=j-1
		end do
		if (i >= j) exit
		t = a(i);  a(i) = a(j);  a(j) = t
		i=i+1
		j=j-1
	end do
	if (first < i-1) call quicksort(a, first, i-1)
	if (j+1 < last)  call quicksort(a, j+1, last)
end subroutine quicksort


end module V0para



