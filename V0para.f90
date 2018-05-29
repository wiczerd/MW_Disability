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


logical :: dbg_skip = .false. !skip stuff for a minimal debug

!**Environmental Parameters**********************************************************************!
real(8), parameter ::	youngD = 15., &	!Length of initial young period
		oldD = 5., &		!Length of each old period
		tlen =12., &		!Number of periods per year (monthly)
		Longev = 82.- 30.,&	!Median longevity
		UIrr = 0.8, &		!Replacement Rate in UI
		eligY  = 0.834,&	!Fraction young who are eligable
		R = dexp(.02/tlen),&	!People can save in the backyard
		upd_zscl = 0.1,&		! rate at which to update zshift
		upd_wgtrnd = 0.01,&		! rate at which update wage_trend
		smth_diaward = 0.1		! number between 0,1 for how much to mix the smoothed awards

integer, parameter :: oldN = 4,&	!4!Number of old periods
		TT = oldN+2, &		!Total number of periods, oldN periods plus young and retired
		itlen = 12		! just an integer version of tlen so I don't have to keep casting
!----------------------------------------------------------------------------!

!**Programming Parameters***********************!
integer, parameter ::	nal = 6,  &!5		!Number of individual alpha types
			ntr = 7,    &!7	        !Number of occupation trend points
			ndi = 2,    &		    !Number of individual disability risk types
			nl	= 2,    &			!Number of finding/separation rates
			nd  = 3,    &		    !Number of disability extents
			ne  = 5,    &!5	        !Points on earnings grid - should be 1 if hearnlw = .true.
			na  = 30,   &!50	    !Points on assets grid
			nz  = 2,    &		    !Number of aggregate shock states
			nj  = 16,   &!16		!Number of occupations
			Nskill  = 2,&			!number of skills that define occupations. First is always physical
			NKpolyT = 1,&			!polynomial order for time trend for occupation or number of spline segments
			NTpolyT = 2,& 			!polynomial order for time trend overall or number of spline segments
			Nknots   = 4,& 			! Number of knots
			maxiter = 2000, &		!Tolerance parameter
			Nsim = 40000,&!5000*nj	!how many agents to draw
			year0 = 1984, &			!when simulation starts and stops
			yearT = 2013, &
			TossYears = 5, & 		!number of years to run and throwaway
			Tsim = itlen*(yearT - year0+1 + TossYears), &	!how many periods to solve for simulation
			init_yrs = 3,&			!how many years for calibration to initial state of things
			init0_yrs= TossYears,&	!how many years buffer before calibration to initial state of things
			struc_brk = 20,&	    ! when does the structural break happen
			Nk   = TT+(nd-1)*2+2,&	!number of regressors - each age-1, each health and leading, occupation dynamics + 1 constant
			fread = 10


! thse relate to how we compute it. Mostly for debugging purposes
logical            :: tr_spline  = .true.,& 	! use spline or global polynomials for fitting trend
					  al_contin  = .true.,&		!make alpha draws continuous or stay on the grid
					  zj_contin	 = .false.,&	!make zj draws continous
					  ineligNoNu = .false.,&	!do young ineligable also pay the nu cost when they are ineligable?
					  dieyoung   = .true.,&		!do the young die (rate associated with health state)
					  w_strchng	 = .true.,&     ! w gets fed a structural change sequence
					  wglev_0	 = .false. ,&  	!should the initial wage level be 0 for all occupations
					  readshocks = .false.		!readshocks from disk?


! these relate to what's changing over the simulation/across occupation
logical           ::  del_by_occ = .true.,& !delta is fully determined by occupation, right now alternative is fully random
					  j_regimes  = .true.,& !different pref shifts
					  j_rand     = .true.,&! randomly assign j, or let choose.
					  demog_dat	 = .true.,& !do the demographics follow
					  wtr_by_occ = .true.,& ! do we feed in occupation-specific trends for wages
					  occ_dat    = .true.,& ! do we
					  NBER_tseq  = .true.,&	!just feed in NBER recessions?
					  RAS_pid    = .true.,& !balance the health transition matrix
					  buscyc	 = .true.

logical			  ::  run_experiments = .false., &
					  run_cal = .false., &
					  refine_cal = .false., &
					  est_elasticity = .false.

real(8), parameter ::  amax 	 = 24.0,   &	!Max on Asset Grid
					   amin = 0.0	   	!Min on Asset Grid


!**To build***************************!
real(8) :: 	alfgrid(nal,nd), &	!Alpha_i grid- individual wage type parameter
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
		pialf(nal,nal,nd),&	!Alpha_i transition matrix
		ergpialf(nal),&		!Alpha_i ergodic distribution
		piz(nz,nz),&		!TFP transition matrix
		ergpiz(nz),& 		!ergodic TFP distribution
		pid(nd,nd,ndi,TT-1),&!Disability transition matrix
!
		prob_age(TT-1,Tsim), &!Probability of being in each age group over the whole history
		prborn_t(Tsim),&	!probability of being born at each point t
		hazborn_t(Tsim), &	!hazard of being born at each point t
		prborn_constpop(Tsim),&	!probability of being born at each point t if population structure stays constant
		hazborn_constpop(Tsim), &	!hazard of being born at each point t if population structure stays constant

		PrDage(nd,TT), &	!Fraction of each D at each age
		PrDageDel(nd,TT,ndi), &	!Ergodic distribution of each D at each Age and Delta (implied by pid)
		PrDeath(nd,TT),&	!Probability of death during work-life
!
		tr_decls(11),& !0,.1,.2,.3,.4,.5,.6,.7,.8,.9,1
		wg_decls(11),& !0,.1,.2,.3,.4,.5,.6,.7,.8,.9,1
!
		jshift(nj,Tsim),&!Preference shift to ensure proper proportions, 2 regimes
		wage_trend(Tsim,nj),&!trend in wages
		wage_lev(nj),&		!occupation-specific differences in wage-level
		!wage_coef(Nskill*2+NTpolyT+5),& !occupation-specific differences in wage-level
		!wage_coef(Nskill*NKpolyT+Nskill+NTpolyT+5),& !occupation-specific differences in wage-level
		sepgrid(nl,nz),&		!grid for separation rates
		fndgrid(nl,nz),&		!grid for finding rates
		sepwt(nl,nj,nz),&		!grid for separation rates
		fndwt(nl,nj,nz),&		!grid for finding rates
		seprt_mul=1.,fndrt_mul=1.,&   !multipliers for separation and finding rates
!		targets for occupations
		seprisk(nz,nj),&	!occupation-cycle specific job separation
		fndrate(nz,nj),&	!occupation-cycle specific job finding rates
		occ_onet(nj,Nskill),&!physical and 3 KSA
		occwg_datcoef_sqr(Nskill+1,NKpolyT+1),& !Only used with NKpolyT>=2. Then these are the data coefficients for wage regression. also includes 0-order and time-only
		occwg_datcoef(Nskill*2+NTpolyT+5),& !!coefficients for wage regression. First is cubic in time, then linear in skill dimension. Then 2 for age profile, 2 for health dummies, 1 const
		occwg_datspline(Nskill+Nknots-1+Nskill*(Nknots-1)+5),& !!coefficients for wage spline regression. First is levels for skills, then cubic spline in time. Then spline for each skill. Then 2 for age profile, 2 for health dummies, 1 const
		occwg_dattrend(Tsim,nj),& !trend in occupation wage
		occwg_datlev(nj),&		!level of occupation wage

		occsz0(nj),&		   !Fraction in each occupation
		occpr_trend(Tsim,nj) !trend in occupation choice

real(8), allocatable :: wage_coef(:), & !occupation-specific differences in wage-level
						tr_knots(:)  !will be the knot points

integer :: 	dgrid(nd)	! just enumerate the d states
real(8)	::	agegrid(TT)		! the mid points of the ages

!***preferences and technologies that may change
real(8) :: 	beta= dexp(-.05/tlen),&	!People are impatient (5% annual discount rate to start)
		nu = 1.e-3, &		!Psychic cost of applying for DI - proportion of potential wage
		util_const = 0.,&	!Give life some value
		Fd(nd) = 0.,&			! Fixed cost of participating in labor

!	Idiosyncratic income process
		alfrho(nd) = 0.988, &	!Peristence of Alpha_i type
		alfmu(nd) = 0.0,&		!Mean of Alpha_i type
		alfcondsig(nd) = 0.015**0.5,&	!Conditional StdDev of Alpha_i type (Normal)

		LTUrr = 0.4,&			!Home production income
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
		proc_time1 =   2.5,&!The average time to decision	(could be 2.5 for 'meets criteria' or 3.64 for initial decision)
		proc_time2 = 14.12,&!The average time to decision	(could be 28.05 for appeal that continues)
		xizcoef    = 0.1, &	!change in acceptance rate with z deterioration
		xizd1coef  = 0.0, &	!change in acceptance rate with z deterioration if d=1
		xizd23coef = 0.1, &	!change in acceptance rate with z deterioration if d=2 or 3
		xiagecoef  = 0.,&	!increase in vocational acceptance due to age
		voc_age	   = 0.25,&	!target for increase in vocation due to age
		xi_d1shift = -0.,&	!worse hlth stage acceptance for d=1
		xi_d3shift = 0.,&	!better hlth stage acceptance for d=3

		DItest1 = 0.3076, &	!Earnings Index threshold 1 (These change based on the average wage)
		DItest2 = 1.8570, &	!Earnings Index threshold 2
		DItest3 = 3.4093,&  !Earnings Index threshold 3
		smth_dicont = 1.	!Smoothing for the di application value
!


! some handy programming terms
integer :: 		Tblock_exp	= 2000,	&	!Expected time before structural change (years)
			Tblock_sim = struc_brk,&		!The actual time before structural change (years). For z_regime
			ialUn	= 1 ,&		! the index of alpha that signifies an exogenously displaced worker
			ialL	= 2 ,&
			tri0	= ntr/2,&	! the index of the trgrid that is 0
			Nnuisance = 3+(nd-1)	!addl parameters in the regression that we never target: constant, age quadratic, health dummies

logical  :: cal_on_iter_wgtrend = .true.
integer  :: cal_niter = 0
real(8)  :: cal_obj = 1.
real(8), allocatable  :: wc_guess_nolev(:), wc_guess_lev(:)
logical  :: cal_on_grad = .false.

!remove this!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
real(8) :: tbase_out(Tsim, Nknots-1)


!**** calibration targets
real(8) :: apprt_target = .01,&	!target for application rates (to be filled below)
		dirt_target = 0.018,&	!target for di rates
		diaward_target = 0.0038,& !target for new award rate
		d1_diawardfrac_target = 0.16,&
		voc_acc_target = 0.25,&		!fraction of admissions from vocational criteria, target 1985
		hlth_acc_target = 0.75,&		!fraction taken based on health criteria, target 1985
		avg_unrt = 0.055,&	!average rate of unemployment over the period.
		avg_undur = 3.,&	! average months of unemployment
		avg_frt   = 0.4,&	! average rate of long-term unemployment
		p1d2_target = -.197,&	! how much less d=2 participate
		p1d3_target = -.649  	! how much less d=3 participate



!Preferences----------------------------------------------------------------!
! u(c,p,d) = 1/(1-gam)*(c*e^(theta*d)*e^(eta*p))^(1-gam)

real(8) :: 	gam	= 1.5, &	!IES
		eta 	= -0.197, &	!Util cost of participation
		theta 	= -0.224		!Util cost of disability

integer :: print_lev, verbose
logical :: simp_concav = .false.

real(8) ::  Vtol = 5e-6 	!Tolerance on V-dist
real(8) ::  simtol =1e-6_dp !tolerance on simulations

contains
subroutine setparams()

	logical, parameter :: lower= .FALSE.
	integer:: i, j, k, t,d, tri,iter,it
	real(8):: summy, emin, emax, step, &
		   alfsig(nd),alfcondsigt(nd),alfrhot(nd),alfsigt(nd), &
		  mean_uloss,sd_uloss,tbase(Nknots-1)

	real(8), parameter :: pival = 4.D0*datan(1.D0) !number pi

	real(8) :: age_occ_read(6,18), age_read(38,TT), maxADL_read(16)
	real(8) :: occbody_trend_read(yearT-year0+1,17),occbody_trend_interp(Tsim,nj)
	real(8) :: UE_occ_read(2,16),EU_occ_read(2,16), apprt_read(50,2), ONET_read(16,4)
	real(8) :: pid_tmp(nd,nd,TT-1),causal_phys_read(16), PrDDp_Age_read(15,4), PrD_Age_read(6,4),pid_in_read(6,5),PrDeath_in_read(15)
	real(8) :: age_read_wkr(38),occpr_read_wkr(yearT-year0+1)
	real(8) :: wage_coef_O2_read(17),wage_coef_O3_read(21),wage_coef_O1_read(22),wage_coef_CS_read(25)


	real(8) :: pid1(nd,nd),r1(nd),s1(nd),PrDage_tp1(nd,TT-1)

	real(8) :: Hdist_read(5,nd+1),Hmat_read(7,9)

	real(8) :: bornin,p1,p2,p3,p1p,p2p,p3p,junk,pi21
	real(8) :: pNy,pNm,Ny,Nm,dy,dm, Nbar,totborn,prH,prL, bN(Tsim)

	real(8) :: wr(nd),wi(nd),vl(nd,nd),vr(nd,nd),abnrm,rcondv(nd),scl(nd),sdec(nd,nd),rconde(nd),wrk(nd*(nd+6))
	integer :: ilo,ihi,iwrk(nd*2-2),lwrk,status
	CHARACTER(80) :: headers


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	verbose = 1
	print_lev = 2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	! Read things in:

	!Read in the occupation size among the young.
	open(unit=fread, file="Male_SOC_Demog.csv")
	read(fread,*) headers !throw away the first row
	do t=1,(yearT-year0+1)
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
	open(unit= fread, file = "Male_Age_Demog.csv")
	read(fread,*) headers
	do i=1,38
		read(fread,*) age_read(i,:)
	enddo
	close(fread)

	if( tr_spline .eqv. .false. ) then
		open(unit= fread, file = "OLSWageTrend_O2.csv")
		do j=1,17
			read(fread,*) wage_coef_O2_read(j)
		enddo
		close(fread)

		open(unit= fread, file = "OLSWageTrend_O3.csv")
		do j=1,21
			read(fread,*) wage_coef_O3_read(j)
		enddo
		close(fread)

		open(unit= fread, file = "OLSWageTrend_O1.csv")
		do j=1,22
			read(fread,*) wage_coef_O1_read(j)
		enddo
		close(fread)
	else
		if(Nskill ==3) then
			open(unit= fread, file = "OLSWageTrend_CS2.csv")
			do j=1,25
				read(fread,*) wage_coef_CS_read(j)
			enddo
		else
			open(unit= fread, file = "OLSWageTrend_CS1.csv")
			do j=1,17
				read(fread,*) wage_coef_CS_read(j)
			enddo
		endif

		close(fread)

	endif

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
	open(unit= fread, file="ONETpca.csv")
	do j=1,nj
		read(fread, *,iostat=k) ONET_read(j,:) !first column is label, then Phys, then non-phys
	enddo
	close(fread)

	open(unit= fread, file="Hdist.csv")
	do j=1,size(Hdist_read,1)
		read(fread, *,iostat=k) Hdist_read(j,:) !first column is label, then Phys, then non-phys
	enddo
	close(fread)
	!Hdist rows:
	!state label
	!age<45
	!age>45 & age<56
	!age>55 & age<61
	!age>60 & age<66


	open(unit= fread, file="HmatIn.csv")
	do j=1,size(Hmat_read,1)
		read(fread, *,iostat=k) Hmat_read(j,:) !first column is label, then Phys, then non-phys
	enddo
	close(fread)
	!Hmat rows:
	!state transition label
	!coef on occ health
	!age>45 & age<56
	!age>55 & age<61
	!age>60 & age<66
	!age>65
	!base rate

	!!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!* OLD input files
	!read in health transitions estimated directly from PSID
	open(unit=fread, file="PrDDp_Age.csv")
	do j=1,15
		read(fread, *,iostat=k) PrDDp_Age_read(j,:)
	enddo
	close(fread)


	!read in the disability rates by age
	open(unit=fread, file="PrD_Age.csv")
	do j=1,6
		read(fread, *,iostat=k) PrD_Age_read(j,:)
	enddo
	close(fread)

	!read in the death rates by age X disability
	open(unit=fread, file="PrDeath_in.csv")
	do j=1,15
		read(fread, *,iostat=k) PrDeath_in_read(j)
	enddo
	close(fread)

	!read in the health transition rates by age --- computed in matlab to match ss rates
	open(unit=fread, file="pid_in.csv")
	do j=1,6
		read(fread, *,iostat=k) pid_in_read(j,:)
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

	!values from GMM with Heckit residuals
	alfrho = (/ 0.95260904 ,0.93609776 ,0.93609776 /)	!Peristence of Alpha_i type
	alfmu = 0.0	!Mean of Alpha_i type
	alfcondsig = (/0.01446509**0.5, 0.02344495**0.5, 0.02344495**0.5/) !Conditional StdDev of Alpha_i type (Normal)



	!Individual Wage component (alpha)- grid= 2 std. deviations
	do d = 1,nd
		alfrhot(d) = alfrho(d)**(1./tlen)
		alfsig(d) = (alfcondsig(d)**2/(1-alfrho(d)**2))**0.5
		alfsigt(d) = (alfsig(d)**2/tlen)**0.5
		alfcondsigt(d) = (alfsigt(d)**2*(1-alfrhot(d)**2))**0.5
		if(d==1) then
			call rouwenhorst(nal-1,alfmu(d),alfrho(d),alfcondsig(d),alfgrid(2:nal,d),pialf(2:nal,2:nal,d))
		else
			alfgrid(2:nal,d)=alfgrid(2:nal,1)
			call tauchen(nal-1,alfmu(d) , alfrho(d),alfcondsig(d),alfgrid(2:nal,d) ,pialf(2:nal,2:nal,d) ,readzvec=.true.)
		endif
		! this is the transformation if I use the annual-frequency shocks and experience them ever 1/tlen periods
		do i=2,nal
			pialf(i,i,d) = pialf(i,i,d)**(1./tlen)
			do k=2,nal
				if(k /= i) pialf(i,k,d) = 1. - (1.-pialf(i,k,d))**(1./tlen)
			enddo
			summy = sum(pialf(i,2:nal,d))
			if(summy /=1 ) pialf(i,2:nal,d)=pialf(i,2:nal,d)/summy !this is just numerical error
		enddo
		alfgrid(1,d) = min(log(LTUrr),alfgrid(2,d)) !aribtrary... meant to be small

		! ergpialf = 0.
		! ergpialf(1) = 0.
		! ergpialf(2) = alnorm( ((alfgrid(3)+alfgrid(2))/2.-alfmu) /alfsig,.false.)
		! if( nal > 3) then
		! 	do i=3,(nal-1)
		! 		ergpialf(i) = ( &
		! 			&	alnorm( ((alfgrid(i+1)+alfgrid(i))/2.-alfmu) /alfsig,.false.)- &
		! 			&	alnorm( ((alfgrid(i-1)+alfgrid(i))/2.-alfmu) /alfsig,.false.) )
		! 	enddo
		! endif
		! ergpialf(nal) = 1.-sum(ergpialf(2:(nal-1)))

		! the probabilities associated with going into the alpha term that is unemployment go to zero.
		pialf(1,1,d) = 0.
		mean_uloss = -0.154996
		!this comes from the SIPP, code in CVW: DTall[ lfstat_wave==1 & seam==T, mean(wavewage,na.rm = T)] - DTall[ lfstat_wave==1 & seam==T & shift(EU_wave,type="lag")==T, mean(wavewage,na.rm = T)]
		sd_uloss   = (0.5888828)**0.5

		pialf(1,nal,d) = 1.- alnorm( ((alfgrid(nal-1,d)+alfgrid(nal,d))/2.-mean_uloss) /sd_uloss,.false.)
		if( nal > 3) then
			do i= 3,(nal-1)
				pialf(1,i,d) = ( &
						&	alnorm( ((alfgrid(i+1,d)+alfgrid(i,d))/2.- mean_uloss ) /sd_uloss,.false.)- &
						&	alnorm( ((alfgrid(i-1,d)+alfgrid(i,d))/2.- mean_uloss ) /sd_uloss,.false.) )
			enddo
		endif
		pialf(1,2,d) = 1.-sum(pialf(1,3:nal ,d))
		pialf(2:nal,1,d) = 0. !exogenous separation is not built into alpha transitions
		!change this!!!!!!!!!!!!!!!
		!pialf(1,3:nal,d) = 0.
		!pialf(1,2,d) = 1.
	enddo

	!~~~~~~~~~~~~~~~~~~~~~~~~~~~
	! Occupation wage component
	do i=1,nj
		do k=2,(Nskill+1)
			occ_onet(i,k-1) = ONET_read(i,k)
		enddo
	enddo

	allocate(tr_knots(Nknots))
	if(Nknots == 5) then
		tr_knots = (/1.,7.,11.,20.,28. /)
	elseif(Nknots == 4) then
		tr_knots = (/1.,9.,16.,28. /)
	else
		tr_knots(1) = 1.
		tr_knots(Nknots) = 28.
		do j=2,Nknots-1
			tr_knots(j) = dble(j)/(28.-1.)+1.
		enddo
	endif

	occwg_datcoef = 0._dp
	occwg_datspline = 0._dp
	occwg_datcoef_sqr   = 0._dp
	if( tr_spline .eqv. .false. ) then
		allocate(wc_guess_nolev(NTpolyT+Nskill+2))
		allocate(wc_guess_lev(NTpolyT+Nskill*2+2))
		if(wglev_0 .eqv. .false.) allocate(wage_coef(NTpolyT+Nskill*2+5))
		if(wglev_0 .eqv. .true.) allocate(wage_coef(NTpolyT+Nskill+5))
		if(NKpolyT >= 2) then
			t=6
			do j=1,(NKpolyT+1)
				do k=1,(Nskill+1)
					if(k > 1 .or. j > 1) then
						if(NKpolyT == 2)	occwg_datcoef_sqr(k,j) = wage_coef_O2_read(t)
						if(NKpolyT == 3)	occwg_datcoef_sqr(k,j) = wage_coef_O3_read(t)
						t = t+1
					else
						if(NKpolyT == 2)	occwg_datcoef_sqr(k,j) = wage_coef_O2_read( (NKpolyT+1)*(Nskill+1) +5 )
						if(NKpolyT == 3)	occwg_datcoef_sqr(k,j) = wage_coef_O3_read( (NKpolyT+1)*(Nskill+1) +5 )
					endif
				enddo
			enddo
		else
			t= 14
			do j=1,NTpolyT !read the NTpolyT in time
				occwg_datcoef(j) = wage_coef_O1_read(t)
				t=t+1
			enddo
			do k=1,(2*Nskill) !level and then trend for each skill
				occwg_datcoef(k+NTpolyT) = wage_coef_O1_read(t)
				t=t+1
			enddo
			occwg_datcoef(1+NTpolyT+2*Nskill) = wage_coef_O1_read(t) !just to get the constant
		endif
	else
		allocate(wc_guess_nolev((Nskill+1)*(Nknots-1)+2))
		allocate(wc_guess_lev( Nknots-1+ Nskill*Nknots +2 ))
		if(wglev_0 .eqv. .true.) allocate(wage_coef((Nskill+1)*(Nknots-1)+5) )
		if(wglev_0 .eqv. .true.) allocate(wage_coef(Nknots-1 + Nskill*Nknots+5) )

		t= 6
		do j=1,Nskill
			occwg_datspline(j) = wage_coef_CS_read(t)
			t = t+1
		enddo
		do j=1,(Nknots-1)
			occwg_datspline(j+Nskill) = wage_coef_CS_read(t)
			t = t+1
		enddo
		do k=1,NSkill
			do j=1,(Nknots-1)
				occwg_datspline(j +(k-1)*(Nknots-1) + Nskill+(Nknots-1)) = wage_coef_CS_read(t)
				t =t+1
			enddo
		enddo
	endif

	!use the coefficients to set the trends (stored in occwg_dattrend):
	if( tr_spline .eqv. .true. ) then
		do j= 1,nj
			do t=1,Tsim
				occwg_dattrend(t,j) =0._dp
				if(t>TossYears*itlen) then
					it = t
				else
					it = TossYears*itlen
				endif
				do k=1,Nskill
					occwg_dattrend(t,j) = occwg_datspline(k)*occ_onet(j,k) +occwg_dattrend(t,j)
				enddo
				tbase = 0._dp
				tbase(1) = (dble(it)/tlen - dble(TossYears))
				do i=1,(Nknots-2)
					if((dble(it)/tlen - dble(TossYears) - tr_knots(i)) > 0.) &
					& 	tbase(i+1) = (dble(it)/tlen - dble(TossYears) - tr_knots(i))**3 + tbase(i+1)
					if( dble(it)/tlen - dble(TossYears) - tr_knots(Nknots-1) >0.) &
					&	tbase(i+1) = -(dble(it)/tlen - dble(TossYears) - tr_knots(Nknots-1))**3 *(tr_knots(Nknots)-tr_knots(i))/(tr_knots(Nknots)-tr_knots(Nknots-1)) &
						&  + tbase(i+1)
					if( dble(it)/tlen - dble(TossYears) - tr_knots(Nknots) >0. ) &
					& 	tbase(i+1) = -(dble(it)/tlen - dble(TossYears) - tr_knots(Nknots) )**3 *(tr_knots(Nknots-1)-tr_knots(i))/(tr_knots(Nknots)-tr_knots(Nknots-1)) &
						&  + tbase(i+1)
					tbase(i+1) = tbase(i+1)*(tr_knots(Nknots)-tr_knots(1))**(-2)
				enddo
				tbase_out(t,:) = tbase
				do i=1,(Nknots-1)
					occwg_dattrend(t,j) = occwg_datspline(i+Nskill)*tbase(i) +occwg_dattrend(t,j)
				enddo
				do k=1,Nskill
					do i=1,(Nknots-1)
						occwg_dattrend(t,j) = occwg_datspline(i+(k-1)*(Nknots-1)+Nskill+Nknots-1)*tbase(i)*occ_onet(j,k) &
								& + occwg_dattrend(t,j)
					enddo
				enddo
			enddo !Tsim
		enddo ! Nj
	else
		if(NKpolyT>=2) then
			do i=1,nj
				do t=1,Tsim
					occwg_dattrend(t,i) = 0._dp
					if(t>TossYears*itlen) then
						do j =1,(NKpolyT+1)
							occwg_dattrend(t,i) =     (dble(t)/tlen-dble(TossYears))**(j-1)*occwg_datcoef_sqr(1,j)                 + occwg_dattrend(t,i)
							do k=1,Nskill
								occwg_dattrend(t,i) = (dble(t)/tlen-dble(TossYears))**(j-1)*occwg_datcoef_sqr(k+1,j)*occ_onet(i,k) + occwg_dattrend(t,i)
							enddo
						enddo
					else
						do j =1,(NKpolyT+1)
							occwg_dattrend(t,i) =     (0.)**(j-1)*occwg_datcoef_sqr(1,j)                 + occwg_dattrend(t,i)
							do k=1,Nskill
								occwg_dattrend(t,i) = (0.)**(j-1)*occwg_datcoef_sqr(k+1,j)*occ_onet(i,k) + occwg_dattrend(t,i)
							enddo
						enddo
					endif
				enddo
			enddo
		else
			do i=1,nj
				do t=1,Tsim
					occwg_dattrend(t,i) = 0._dp
					if(t>TossYears*itlen) then
						do j =1,NTpolyT !unroll time trend
							occwg_dattrend(t,i) = (dble(t)/tlen-dble(TossYears))**j*occwg_datcoef(j) + occwg_dattrend(t,i)
						enddo
						do k=1,Nskill !occupation-specific ONET prices
							occwg_dattrend(t,i) = (dble(t)/tlen-dble(TossYears))*occwg_datcoef(k+NTpolyT+ Nskill)*occ_onet(i,k) &
							& + occwg_datcoef(k+NTpolyT)*occ_onet(i,k) + occwg_dattrend(t,i)
						enddo
					else
						do j =1,NTpolyT !unroll time trend
							occwg_dattrend(t,i) = (0.)**j*occwg_datcoef(j) + occwg_dattrend(t,i)
						enddo
						do k=1,Nskill !occupation-specific ONET prices
							occwg_dattrend(t,i) = (0.)*occwg_datcoef(k+NTpolyT+ Nskill)*occ_onet(i,k) &
							& + occwg_datcoef(k+NTpolyT)*occ_onet(i,k) + occwg_dattrend(t,i)
						enddo
					endif
				enddo
			enddo
		endif !square coefficients or linear
	endif !spline or global polynomial

	!initialize the input to the observed
	do i=1,nj
		occwg_datlev(i) = occwg_dattrend(TossYears*itlen,i)
		occwg_dattrend(:,i) = occwg_dattrend(:,i) - occwg_datlev(i)
	enddo
	if(wglev_0 .eqv. .true.) then
		wage_lev = 0._dp
	else
		wage_lev = occwg_datlev
	endif
	wage_trend = occwg_dattrend

	!Wage-trend grid-setup
	trgrid(1) = (minval(wage_trend(Tsim,:) + wage_lev ))*1.1_dp
	if(ntr>1) then
		trgrid(ntr) = (maxval(wage_trend(Tsim,:)+wage_lev))*1.1_dp
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

	!initialize tr_decls
	tr_decls(1) = minval(trgrid)
	tr_decls(11)= maxval(trgrid)
	do i=2,10
		tr_decls(i) = dble(i-1)/10._dp*(tr_decls(11)-tr_decls(1))
	enddo
	wg_decls(1) = dexp(minval(alfgrid)+minval(trgrid))
	wg_decls(11)= dexp(maxval(alfgrid)+maxval(trgrid))
	do i=2,10
		wg_decls(i) = dble(i-1)/10._dp*(wg_decls(11)-wg_decls(1))
	enddo

	!read these numberrs in already
	seprisk = 0._dp
	fndrate = 0._dp
	do i =1,nz
		do j=1,nj
			k = 3-i !flip, recession is 1 and expansion is 2 in markov chain
			fndrate(i,j) = UE_occ_read(k,j)
			seprisk(i,j) = EU_occ_read(k,j)
		enddo
	enddo

	do i=1,nz
		do j=1,nl

			fndgrid(j,i) = (maxval(fndrate(i,:))-minval(fndrate(i,:))) *dble( j-1 )/dble(nl-1) + minval(fndrate(i,:))
			sepgrid(j,i) = (maxval(seprisk(i,:))-minval(seprisk(i,:))) *dble( j-1 )/dble(nl-1) + minval(seprisk(i,:))

		enddo
	enddo

	! TFP
	zshift = 0._dp
	zscale = 1._dp

	call settfp()


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!Demographics

	! age grid building
	agegrid(1) = youngD/2
	do i = 2,TT-1
		agegrid(i) = youngD + oldD*(i-2) + oldD/2
	enddo
	agegrid(TT) = Longev/2 - (youngD + oldD*oldN) /2 + (youngD + oldD*oldN)

	agegrid = agegrid + 10._dp !started them at 30, but they've been working since 20

	!from Mincer regressions (see Appendix) with Heckman 2-step
	wtau(1) =  0.
	wtau(2) = -0.032/(.5*agegrid(2)+.5*agegrid(3)-agegrid(1))*(agegrid(2)-agegrid(1))
	wtau(3) = -0.032/(.5*agegrid(2)+.5*agegrid(3)-agegrid(1))*(agegrid(3)-agegrid(1))
	wtau(4) = -0.12
	wtau(5) = -0.174

	!Aging Probability (actually, probability of not aging)
	! Mean Duration = (pr(age))^(-1)-1 <--in 1/tlen units
	ptau(1) = 1-(tlen*youngD)**(-1)

	do i=2,TT-1
		ptau(i) = 1-(tlen*oldD)**(-1)
	enddo
	! rate exit retirement (only one way to go.... down)
	ptau(TT) = 1-((Longev-(youngD+oldN*oldD))*tlen)**(-1)

	! prob of death by health, age
	PrDeath(:,1) = Hmat_read(7,7:9)
	do t=1,TT
		k = t+1
		if( t .eq. 2 ) k = t+2
		if( t>1) PrDeath(:,t) = PrDeath(:,1) +  Hmat_read(k,7:9)
		PrDeath(:,t) = 1.d0 - (1.d0 - PrDeath(:,t))**(1.d0/tlen)  !PrDeath_in_read(1+(t-1)*3:3+(t-1)*3)
	enddo


	!Health structure by age
	do t=1,TT
		if(t >2) then
			k = t
		else
			k =t+1
		endif
		if(t==TT) k = t-1
		if( t .ne. 2 .and. t .ne. 3 ) then
			PrDage(:,t) = Hdist_read(k,2:1+nd)
		elseif(t .eq. 2) then
			do i=1,nd
				PrDage(i,t) = (Hdist_read(3,1+i) - Hdist_read(2,1+i))/(.5_dp*(agegrid(3)+agegrid(2)) - agegrid(1))*(agegrid(2)-agegrid(1)) + Hdist_read(2,1+i)
			enddo
		elseif(t .eq. 3) then
			do i=1,nd
				PrDage(i,t) = (Hdist_read(4,1+i) - Hdist_read(3,1+i))/(agegrid(4) - .5_dp*(agegrid(3)+agegrid(2)) )*(agegrid(3)-agegrid(4)) + Hdist_read(4,1+i)
			enddo
		endif
		PrDage(:,t) = PrDage(:,t)/sum(PrDage(:,t))
	enddo
	!health structure extrapolate one period ahead - make the transition rate correct
!~ 	PrDead_tp1 = 0.d0
	do i=1,nd
		do t=1,TT-2
			PrDage_tp1(i,t) = (PrDage(i,t+1)-PrDage(i,t))/(agegrid(t+1)-agegrid(t))+PrDage(i,t)
			!some at each health/age will die. Add those in  too
			PrDage_tp1(i,t) = PrDage_tp1(i,t)*(1.d0 + (1.d0-(1.d0-PrDeath(i,t))**tlen))
		enddo
		PrDage_tp1(i,TT-1) = PrDage(i,TT-1)*(1.d0+ (1.d0-(1.d0-PrDeath(i,TT-1))**tlen))
	enddo

	!age structure extrapolate over periods
	do i=1,TT-1
		call spline( age_read(:,1),age_read(:,i+1),age_read_wkr)
		do t=(TossYears*itlen+1),Tsim
			prob_age(i,t) = splint(age_read(:,1),age_read(:,i+1),age_read_wkr, dble(t-1)/tlen-TossYears+year0 )
		enddo
		do t=1,(TossYears*itlen)
			prob_age(i,t) = prob_age(i,TossYears*itlen+1)
		enddo
	enddo

	!evolution of age-structure!!!!!!!!!!!!!!!!!!!!!!

	if(dieyoung .eqv. .true.) then
		dy = PrDeath(1,1)*PrDage(1,1)+PrDeath(2,1)*PrDage(2,1)+PrDeath(3,1)*PrDage(3,1)
		t=2
		dm = PrDage(1,t)*PrDeath(1,t) + PrDage(2,t)*PrDeath(2,t) + PrDage(nd,t)*PrDeath(nd,t)
		do t=3,TT-1
			dm =  (PrDage(1,t)*PrDeath(1,t) + PrDage(2,t)*PrDeath(2,t) + PrDage(nd,t)*PrDeath(nd,t))*(1._dp-dm) +dm !die
		enddo
		dm = dm/dble(TT-1-2)
		!just using from the simulation
		dm = 0.002250656
		dy = 0.0002027271
	else
		dm = 0.d0
		dy = 0.d0
	endif
	!dm = dm + ((tlen*oldN*oldD)**(-1))*(1.d0-dm)
	!dm = 0.009

	prH = 1.0d0
	prL = 0.d0
	do i =1,maxiter
		!prob of getting born
		t=1
		hazborn_t(t) = 0.5d0*(prH + prL)
		Ny = hazborn_t(t)*dble(Nsim)*prob_age(1,t)
		Nm = hazborn_t(t)*dble(Nsim)*(1.-prob_age(1,t))
		pNy = Ny*ptau(1)*(1.d0-dy)
		pNm	= Nm*(1.d0-dm)*(1.d0-(tlen*oldN*oldD)**(-1)) + Ny*(1.d0-ptau(1))*(1.d0-dy)
		bN(t) = (prob_age(1,1)*(pNy+pNm)-pNy)/(1.-prob_age(1,1))
		totborn = hazborn_t(t)*dble(Nsim)
		do j=1001,(Tsim+999)
			if(j>1000) then
				t=t+1
			endif
			pNy = Ny*ptau(1)*(1.d0-dy)
			pNm	= Nm*(1.d0-dm)*(1.d0-(tlen*oldN*oldD)**(-1)) + Ny*(1.d0-ptau(1))*(1.d0-dy)
			if(t>1) then
				bN(t) = max( (prob_age(1,t)*(pNy+pNm)-pNy)/(1.-prob_age(1,t)), 0.d0)
				hazborn_t(t) = bN(t)/(Nsim - totborn) !hazborn*(remaining unborn) = bN
				totborn = bN(t) + totborn
			else
				bN(t) = max( (prob_age(1,1)*(pNy+pNm)-pNy)/(1.-prob_age(1,1)), 0.d0)
			endif
			Nm = pNm
			Ny = pNy + bN(t)
		enddo
		junk = hazborn_t(1)
		hazborn_t(1) =  (dble(Nsim) - (totborn - hazborn_t(1)*Nsim ))/dble(Nsim) ! need to have some positive mass alive when the survey starts

		if(hazborn_t(1)<0) hazborn_t(1) = 0.d0
		!if(dabs(junk - hazborn_t(1))<1e-8) then
		if(dabs(totborn - dble(Nsim))<1e-6) then
			exit !iterate on the numberr alive in period 1
		elseif( totborn > dble(Nsim) ) then
			prH = junk
		else ! totborn<Nsim
			prL = junk
		endif
		prborn_t(1) = hazborn_t(1)
		prborn_t(2:Tsim) = bN(2:Tsim)/totborn

	enddo
!again for the constant population group
	prH = 1.d0
	prL = 0.d0
	do i =1,maxiter
		!prob of getting born
		t=1
		hazborn_constpop(t) = 0.5d0*(prH + prL)
		Ny = hazborn_constpop(t)*dble(Nsim)*prob_age(1,t)
		Nm = hazborn_constpop(t)*dble(Nsim)*(1.-prob_age(1,t))
		pNy = Ny*ptau(1)*(1.d0-dy)
		pNm	= Nm*(1.d0-dm)*(1.d0-(tlen*oldN*oldD)**(-1)) + Ny*(1.d0-ptau(1))*(1.d0-dy)
		bN(t) = (prob_age(1,1)*(pNy+pNm)-pNy)/(1.-prob_age(1,1))
		totborn = hazborn_constpop(t)*dble(Nsim)
		do j=1001,(Tsim+999)
			pNy = Ny*ptau(1)*(1.d0-dy)
			pNm	= Nm*(1.d0-dm)*(1.d0-(tlen*oldN*oldD)**(-1)) + Ny*(1.d0-ptau(1))*(1.d0-dy)
			if(j>1000) then
				t=t+1
				bN(t) = max( (prob_age(1,1)*(pNy+pNm)-pNy)/(1.-prob_age(1,1)) ,0.d0)
				hazborn_constpop(t) = bN(t)/(Nsim - totborn) !hazborn*(remaining unborn) = bN
				totborn = bN(t) + totborn
			else
				bN(t) = (prob_age(1,1)*(pNy+pNm)-pNy)/(1.-prob_age(1,1))
			endif
			Nm = pNm
			Ny = pNy + bN(t)
		enddo
		junk = hazborn_constpop(1)

		hazborn_constpop(1) =  (dble(Nsim) - (totborn - hazborn_constpop(1)*Nsim ))/dble(Nsim) ! need to have some positive mass alive when the survey starts

		if(hazborn_constpop(1)<0) hazborn_constpop(1) = 0.d0
		if(dabs(totborn - dble(Nsim))<1e-6) then
			exit !iterate on the number alive in period 1
		elseif( totborn > dble(Nsim) ) then
			prH = junk
		else ! totborn<Nsim
			prL = junk
		endif
		prborn_constpop(1) = hazborn_constpop(1)
		prborn_constpop(2:Tsim) = bN(2:Tsim)/totborn
	enddo


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!	occupation structure

	do i=1,nj
		call spline( occbody_trend_read(:,1),occbody_trend_read(:,i+1),occpr_read_wkr)
		do t=(TossYears*itlen+1),Tsim
			occbody_trend_interp(t,i) = splint(occbody_trend_read(:,1),occbody_trend_read(:,i+1),occpr_read_wkr, dble(t-1)/tlen+year0-TossYears)
		enddo
		do t=1,(TossYears*itlen)
			occbody_trend_interp(t,i) = occbody_trend_interp(TossYears*itlen+1,i)
		enddo
	enddo
	summy = 0.
	do j=1,nj
		occpr_trend(1,j) = occbody_trend_interp( 1,j )
		summy = occpr_trend(1,j) + summy
	enddo
	occpr_trend(1,:) = occpr_trend(1,:)/summy
	do t=2,Tsim
		summy = 0.
		if( prborn_t(t) .gt. 0. ) then
			do j=1,nj
				occpr_trend(t,j) = (occbody_trend_interp( t,j ) - occbody_trend_interp(t-1,j )*(1.-ptau(1)))/prborn_t(t)
				summy = occpr_trend(t,j) + summy
			enddo
			occpr_trend(t,:) = occpr_trend(t,:)/summy
		else
			do j=1,nj
				occpr_trend(t,j) = 0.
			enddo
		endif
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
	summy =0._dp
	do j=1,nj
		summy = occsz0(j) + summy
	enddo
	forall(j=1:nj) occsz0(j) = occsz0(j)/summy

! Disability grid
	forall(i=1:nd) dgrid(i) = i

! Disability depreciation by occupation
	!occupation-specific factor
	do j=1,nj
		occdel(j) =  occ_onet(j,1)
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
	if(del_by_occ .eqv. .false.) delgrid = 0.

	!Disability Extent-Specific Things
	!Wage Penalty
	wd(1) = 0.		!Healthy, no penalty
	wd(2) = -0.097	!Partially Disabled, small penalty
	wd(3) = -0.266	!Full Disabled, large penalty
	!Fixed cost of particpation
	Fd(1) = 0.
	Fd(2) = 0.276*5
	Fd(3) = 0.524*5


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

	!Earnings Grid
	!Make linear from lowest possible wage (disabled entrant, lowest types)
	emin = dexp(alfgrid(2,1)+minval(wtau)+minval(wd))
	!... to highest, maximizing over t
	emax = DItest3 !dexp(maxval(alfgrid)+maxval(wtau)+maxval(wd))
	step = (DItest2-emin)/dble(ne-2)
	do i=1,ne-1
		egrid(i) = emin+step*dble(i-1)
	enddo
	egrid(ne) = emax

	!Assets Grid
	do i=1,na
		agrid(i)=dble(i-1)/dble(na-1)
		agrid(i)=agrid(i)**2*(amax-amin)+amin
	enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!Make 2-year Markov transition matrices with all wierd invidual stuff
	!Disability: pid(id,id';i,t) <---indv. type and age specific

	!read in transition matrix

	pid_tmp(1,2:3,1) = Hmat_read(7,1:2)
	pid_tmp(2,1,1)   = Hmat_read(7,3)
	pid_tmp(2,3,1)   = Hmat_read(7,4)
	pid_tmp(3,1:2,1) = Hmat_read(7,5:6)

	do t=2,TT-1
		k = t
		if(t .eq. 2) k = t+1
		pid_tmp(1,2:3,t) = Hmat_read(k,1:2) + pid_tmp(1,2:3,1)
		pid_tmp(2,1,t)   = Hmat_read(k,3)   + pid_tmp(2,1,1)
		pid_tmp(2,3,t)   = Hmat_read(k,4)   + pid_tmp(2,3,1)
		pid_tmp(3,1:2,t) = Hmat_read(k,5:6) + pid_tmp(3,1:2,1)
	enddo

	do t=1,TT-1
		pid_tmp(1,1,t) = 1.d0 - sum(pid_tmp(1,2:3,t))
		pid_tmp(2,2,t) = 1.d0 - pid_tmp(2,1,t) - pid_tmp(2,3,t)
		pid_tmp(3,3,t) = 1.d0 - sum(pid_tmp(3,1:2,t))
	enddo

	pid = 0.
	! multiply by delgrid (was a 2-year transition matrix)
	do i=1,TT-1

		pid1 = pid_tmp(:,:,i)
		if(RAS_pid .eqv. .true.) then
		! implement RAS method to balance matrix and enforce steady-state levels
		! this means the row marginals, u_i = 1 \forall i and column marginals v_j
			do iter=1,maxiter
				do k=1,nd
					r1(k) = ( (1.d0 - PrDeath(k,i))**(tlen) ) /sum(pid1(k,:))
				enddo
				!sdec = diag(r1)*pid1:
				do k=1,nd
				do j=1,nd
					sdec(k,j) = r1(k)*pid1(k,j)
				enddo
				enddo

				do k=1,nd
					s1(k) = PrDage_tp1(k,i)/sum(PrDage(1:nd,i)*sdec(:,k))
				enddo
				!pid1 = sdec*diag(s1);
				do k=1,nd
				do j=1,nd
					pid1(k,j) = s1(j)*sdec(k,j)
				enddo
				enddo
				!will end when r1 and s1 both approach 1
				if (dabs(sum(r1)/dble(nd)+ sum(s1)/dble(nd) -2._dp) < 1e-7 ) then
					exit
				endif

			enddo
			do k=1,nd
				pid1(k,:) = pid1(k,:)/sum(pid1(k,:))
			enddo
		endif !RAS

		do j=1,ndi

			pid1(1,2) = pid1(1,2) + delgrid(j)*Hmat_read(2,1)
			pid1(1,3) = pid1(1,3) + delgrid(j)*Hmat_read(2,2)
			pid1(1,1) = 1._dp-pid1(1,2)-pid1(1,3)

			pid1(2,1) = pid1(2,1) + delgrid(j)*Hmat_read(2,3)
			pid1(2,3) = pid1(2,3) + delgrid(j)*Hmat_read(2,4)
			pid1(2,2) = 1._dp-pid1(2,1)-pid1(2,3)

			pid1(3,1) = pid1(3,1) + delgrid(j)*Hmat_read(2,5)
			pid1(3,2) = pid1(3,2) + delgrid(j)*Hmat_read(2,6)
			pid1(3,3) = 1._dp-pid1(3,1)-pid1(3,2)


		!convert to monthly------------
			sdec = pid1
			lwrk = nd*(nd+6)
			!want to construct right-side eigen vectors into matrix
						!BALANC, JOBVL, JOBVR, SENSE, N , A   , LDA, WR, WI,
			call dgeevx( 'N'   , 'N'  , 'V'  , 'V'  , nd, sdec, nd , wr, wi, &
		&	 vl, nd  , vr, nd  , ilo, ihi, scl, abnrm,rconde, rcondv, wrk , lwrk , iwrk , status )
			!VL, LDVL, VR, LDVR, ILO, IHI, SCALE, ABNRM,RCONDE, RCONDV, WORK, LWORK, IWORK, INFO )

			!replace vl = vr^-1
			call invmat(vr, vl)
			!vr = vr*wr^(1/t)
			do t=1,nd
				vr(:,t) = vr(:,t)*wr(t)**(1._dp/tlen)
			enddo
			call dgemm('N', 'N', nd, nd, nd, 1._dp, vr, nd, vl, nd, 0._dp, pid(:,:,j,i), nd)

			do t=1,nd
				summy = sum(pid(t,:,j,i) )
				pid(t,:,j,i) = pid(t,:,j,i)/summy
			enddo

			!initialize
			!PrDageDel(:,:,j) = PrDage
			sdec = pid1
			call  dgeev( 'V', 'N', nd, sdec, nd, wr, wi, vl, nd, vr, nd, &
		    	& wrk, lwrk, status)
			do t=1,nd
				if( dabs(wr(t)-1._dp)<1.e-4 ) &
				&	PrDageDel(:,i,j) = vl(:,t)
			enddo
			summy = sum(PrDageDel(:,i,j))
			PrDageDel(:,i,j) = PrDageDel(:,i,j)/summy
		enddo !j=1,ndi
	enddo
	PrDageDel(:,TT,:) = PrDageDel(:,TT-1,:)

end subroutine setparams

subroutine settfp()

	integer :: i, k,j,nzblock
	logical, parameter :: lower= .FALSE.
	real(8) ::  summy, zcondsig
	real(8) :: zrhot,zsigt, zcondsigt
	real(8), allocatable :: piblock(:,:),ergpi1(:,:),ergpi2(:,:)
	real(8) :: Zzgrid(nz)

	external dgemm

	nzblock = nz

	allocate(piblock(nzblock,nzblock))
	allocate(ergpi1(nzblock,nzblock))
	allocate(ergpi2(nzblock,nzblock))
	piblock = 0.
	ergpi1  = 0.
	ergpi2  = 0.

	!reset separation/finding probabilities:
	do i=1,nz
		do j=1,nl

			fndgrid(j,i) = (maxval(fndrate(i,:))-minval(fndrate(i,:))) *dble( j-1 )/dble(nl-1) + minval(fndrate(i,:))
			sepgrid(j,i) = (maxval(seprisk(i,:))-minval(seprisk(i,:))) *dble( j-1 )/dble(nl-1) + minval(seprisk(i,:))

		enddo
	enddo

	if( nz>2 ) then !using z process

		zrhot = zrho**(1./tlen)
		zsigt = zsig**(1/tlen)
		zcondsig = ((zsig**2)*(1.-zrho**2))**(0.5)
		zcondsigt = ((zsigt**2)*(1.-zrhot**2))**(0.5)
		!first set transition probabilities at an annual basis
		call rouwenhorst(nzblock,zmu,zrho,zcondsig,zgrid(1:nzblock,1),piz(1:nzblock,1:nzblock))
		piblock= piz(1:nzblock,1:nzblock)
		Zzgrid = zgrid(:,1)
		do j=1,nj
			zgrid(:,j) = Zzgrid*zscale(j)
		enddo
		!DGEMM('N','N',  M,  N,    K, ALPHA,  A,     M,    B,       K,  BETA,     C,       M)
		call dgemm('n','n',nzblock,nzblock,nzblock,1._dp,piblock,nzblock, piblock, nzblock, 0._dp, ergpi1,nzblock)
		do i=1,10000
				call dgemm('n','n',nzblock,nzblock,nzblock,1._dp,piblock,nzblock, ergpi1, nzblock, 0._dp, ergpi2,nzblock)
				ergpi1 = ergpi2
		enddo
	else !nz=2 => taking nber probabilities
		! probability of entering recession  4/(58+12+92+120+73) = 0.011267606
		piz(2,1) = 4./(58.+12.+92.+120.+73.)
		! probability of exit from recession 4/(6+16+8+8+18) = 0.071428571
		piz(1,2) = 4./(6.+16.+8.+8.+18.)
		piz(1,1) = 1.-piz(1,2)
		piz(2,2) = 1.-piz(2,1)
		piblock = piz
		forall(i=1:nz) Zzgrid(i) = dble(i)
		do j=1,nj
			forall(i=1:nz) zgrid(i,j) = dble(i)
		enddo
		ergpi1(1,:) = (/0.136253, 0.863747 /)
		ergpi1(2,:) = (/0.136253, 0.863747 /)
	endif

	ergpiz = ergpi1(1,:)

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


subroutine tauchen(N, mu_eps, rho, sigma_eps, zvect, Pmat,readzvec)

    implicit none

    real(8), intent(in):: rho, mu_eps, sigma_eps
    integer, intent(in):: N
    real(8)            :: zvect(N)
    real(8), intent(out):: Pmat(N,N)
	logical, optional   :: readzvec

    real(8) mu_z, sigma_z, sigma_zcond
    real(8) w, m
    integer i,j

    mu_z = mu_eps/(1-rho)
    sigma_z = sigma_eps/sqrt(1-rho**2)
    sigma_zcond = sigma_eps

	m = 3.
	if(present(readzvec) .eqv. .false.) readzvec = .false.
	if(readzvec .eqv. .false. ) then
		do i=1,N
    		zvect(i) = dble(i-1)/dble(N-1)*(m*sigma_z+m*sigma_z)+ mu_z-m*sigma_z
		enddo
	endif

    do i=1,N
		w = (zvect(2)-zvect(1))/2
		Pmat(i,1) = alnorm((zvect(1) + w - rho*zvect(i) - mu_eps)/sigma_eps, .false. )
        do j=2,N-1
			w = (zvect(j+1)-zvect(j))/2
        	Pmat(i,j) = alnorm((zvect(j) + w- rho*zvect(i) - mu_eps)/sigma_eps, .false.) - &
            	alnorm((zvect(j) - w- rho*zvect(i) - mu_eps)/sigma_eps, .false.)
        end do
		w = (zvect(N)-zvect(N-1))/2
        Pmat(i,N) = 1 - alnorm((zvect(N) - w - rho*zvect(i) - mu_eps)/sigma_eps, .false.)
    end do

end subroutine tauchen

subroutine invmat(A, invA)
	real(dp), dimension(:,:), intent(in) :: A
	real(dp), dimension(size(A,1),size(A,2)), intent(out) :: invA
	real(dp), dimension(size(A,1)*size(A,1)) :: wk
	integer, dimension(size(A,1) + 1) :: ipiv
	integer :: n, info
	external DGETRF
	external DGETRI

	invA = A
	n = size(A,1)
	ipiv = 0
	wk = 0.
	info = 0
	call DGETRF(n,n,invA,n,ipiv,info) ! first computes the LU factorization

	if (info /= 0) then
		print *, 'Matrix is singular.  Make it less singular before proceeding'
	endif
	call DGETRI(n,invA,n,ipiv,wk,n,info)
	if(info /=0) then
		print *, 'Matrix inversion failed, though it is not singular'
	endif
end subroutine invmat


subroutine rand_num_closed(harvest)
	!ensures we're drawing uniform [0,1] on a closed interval
	real(8), intent(inout) :: harvest
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
!	integer, dimension(t+1)        :: ip
	integer, dimension(t)          :: iwork
!	integer, dimension(5*t-6)      :: ia,ja,jlu
	integer                        :: i,info
	real(8), intent(in)               :: phi
	real(8), dimension(t),intent(in)  :: y
	real(8), dimension(t),intent(out) :: yd,yt
	real(8), dimension(t,t)           :: a
!	real(8), dimension(5*t-6)         :: s,alu
!	real(8), dimension(t,nm+1)        :: v
!	real(8), dimension(t-2)           :: x2
!	real(8), dimension(t-3)           :: x3
!	real(8), dimension(t-4)           :: x4

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
	real(8)  a(*), x, t
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


! Local functional approximation by cubic and linear splines
subroutine spline(x,y,y2,yp1,ypn)
! this is the NR version of the spline

	IMPLICIT NONE
	REAL(8), DIMENSION(:), INTENT(IN) :: x,y
	REAL(8), DIMENSION(:), INTENT(OUT) :: y2
	REAL(8), INTENT(IN), OPTIONAL :: yp1,ypn
	REAL(8), DIMENSION(size(x)-1) :: u
	real(8) :: p,qn,si,un
	INTEGER :: n,i,k
	n=size(x)
	IF (size(y)/=n .or. size(y2)/=n) THEN
		PRINT *, 'spline: x,y and y2 must be of the same size'
		STOP 'program terminated by spline'
	END IF

	IF (present(yp1)) THEN
		y2(1)=-0.5
		u(1) = (1.0/(x(2)-x(1) ))*( (y(2)-y(1))/(x(2)-x(1))-yp1)
	ELSE
		y2(1)=0.0
		u(1)=0.0
	END IF
	do i =2,n-1
		si	= (x(i)-x(i-1))/(x(i+1)-x(i-1))
		p	= si*y2(i-1)+2.0
		y2(i)	= (si-1.0)/p
		u(i)	= (6.0*((y(i+1)-y(i))/(x(i+1)-x(i))-(y(i)-y(i-1)) &
			& /(x(i)-x(i-1)))/(x(i+1)-x(i-1))-si*u(i-1))/p
	enddo

	IF (present(ypn)) THEN
		qn = 0.5
		un = (3.0/(x(n)-x(n-1)))*( ypn-(y(n)-y(n-1))/(x(n)-x(n-1)) )
	ELSE
		qn = 0.0
		un = 0.0
!		y2(n)=y2(n-1)
!		a(n-1)=0.0
	END IF
	y2(n) = (un-qn*u(n-1))/(qn*y2(n-1)+1.0)
	do k=n-1,1,-1
		y2(k) = y2(k)*y2(k+1)+u(k)
	enddo



end subroutine spline


FUNCTION splint(x,y,y2,xi)
! cubic interpolation of function y on grid x at interpolation point xi
	IMPLICIT NONE
	REAL(8), DIMENSION(:), INTENT(IN) :: x,y,y2
	REAL(8), INTENT(IN) :: xi
	REAL(8) :: splint
	REAL(8) :: a,b,d,xhr
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
		splint=a*y(i)+b*y(i+1)+((a**3-a)*y2(i)+(b**3-b)*y2(i+1))*(d**2)/6.0
	elseif( xhr .ge. x(n) ) then
		splint = (y(n)-y(n-1))/(x(n)-x(n-1))*(xhr -x(n)) + y(n)
	elseif( xhr .le. x(1) ) then
		splint = (y(2)-y(1))/(x(2)-x(1))*(xhr -x(1)) + y(1)
	endif

END FUNCTION splint

FUNCTION dsplint(x,y,y2,xi)
! derivative implied by cubic interpolation of function y on grid x at interpolation point xi
	IMPLICIT NONE
	REAL(8), DIMENSION(:), INTENT(IN) :: x,y,y2
	REAL(8), INTENT(IN) :: xi
	REAL(8) :: dsplint
	REAL(8) :: a,b,d, xhr
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
	dsplint=(y(i+1)-y(i))/d+((3*b**2-1)*y2(i+1)-(3*a**2-1)*y2(i))*d/6.0
END FUNCTION dsplint

FUNCTION linint(x,y,xi)
! linear interpolation of function y on grid x at interpolation point xi
	IMPLICIT NONE
	REAL(8), DIMENSION(:), INTENT(IN) :: x,y
	REAL(8), INTENT(IN) :: xi
	REAL(8) :: linint
	REAL(8) :: a,b,d,xhr
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
	REAL(8), DIMENSION(:), INTENT(IN) :: x,y
	REAL(8), INTENT(IN) :: xi
	REAL(8) :: dlinint
	REAL(8) :: dx,dy
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
    REAL(8), DIMENSION(:), INTENT(IN)  :: x,y,xi
    REAL(8), DIMENSION(:), INTENT(OUT) :: yi
    REAL(8) :: a,b,d
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
	real(8), dimension(:), intent(in) :: x,y
	real(8), dimension(:,:), intent(in) :: f
	real(8), dimension(:), intent(in) :: xiyi
	real(8) :: fq11,fq21,fq12,fq22, dx,dy

	real(8) :: bilinint
	integer  :: x1,x2,y1,y2

	if(size(x)/=size(f,1)) stop 'x,f grids not the same length in bilinear interpolation'
	if(size(y)/=size(f,2)) stop 'y,f grids not the same length in bilinear interpolation'

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
	real(8), dimension(:), intent(in) :: x,y
	real(8), dimension(:,:), intent(in) :: f
	real(8), dimension(:), intent(in) :: xiyi
	real(8) :: fq11,fq21,fq12,fq22, dx,dy

	real(8), dimension(2) :: dxdy
	integer  :: x1,x2,y1,y2

	if(size(x)/=size(f,1)) stop 'x,f grids not the same length in bilinear interpolation'
	if(size(y)/=size(f,2)) stop 'y,f grids not the same length in bilinear interpolation'

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

	fq11 = f(x1,y1)
	fq21 = f(x2,y1)
	fq12 = f(x1,y2)
	fq22 = f(x2,y2)

	dxdy(1) = ( -fq11*(y(y2)-xiyi(2)) &
		&+ fq21*(y(y2)-xiyi(2))  &
		&- fq12*(xiyi(2)-y(y1))  &
		&+ fq22*(xiyi(2)-y(y1)))/(dx*dy)

	dxdy(2) = (-fq11*(x(x2)-xiyi(1)) &
		&+ fq21*(xiyi(1)-x(x1))  &
		&- fq12*(x(x2)-xiyi(1))  &
		&+ fq22*(xiyi(1)-x(x1)))/(dx*dy)

end subroutine

function bilinint_v(x,y,f,xiyi)
	implicit none
	real(8), dimension(:), intent(in) :: x,y
	real(8), dimension(:), intent(in) :: f
	real(8), dimension(:), intent(in) :: xiyi
	real(8) :: fq11,fq21,fq12,fq22, dx,dy

	real(8) :: bilinint_v
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
	real(8), dimension(:), intent(in) :: x,y
	real(8), dimension(:,:), intent(in) :: f,coefs
	real(8), dimension(:), intent(in) :: xiyi

	real(8) :: bisplint

	bisplint = bilinint(x,y,f,xiyi)

end function bisplint

PURE FUNCTION locate(xx,x)
	IMPLICIT NONE
	REAL(8), DIMENSION(:), INTENT(IN) :: xx
	REAL(8), INTENT(IN) :: x
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

end module V0para
