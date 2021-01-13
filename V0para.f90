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
character(len=18) :: caselabel
character(len=10) :: callog = "callog.log"
character(len=15) :: buffer
integer           :: fcallog = 7

integer, parameter:: dp=kind(0.d0) ! double precision


logical :: dbg_skip = .false. !skip stuff for a minimal debug

!**Environmental Parameters**********************************************************************!
real(8), parameter ::	youngD = 15., &	!Length of initial young period
		oldD = 5., &		!Length of each old period
		tlen =12., &		!Number of periods per year (monthly)
		Longev = 82.- 30.,&	!Median longevity
		UIrr = 0.75, &		!Replacement Rate in UI
		eligY  = 0.834,&	!0.834 : Fraction young who are eligable
		R = dexp(.016/tlen),&	!People can save in the backyard
		upd_zscl = 0.1,&		! rate at which to update zshift
		upd_wgtrnd = 0.01,&		! rate at which update wage_trend
		smth_diaward = 0.05		! number between 0,1 for how much to mix the smoothed awards


!Preferences----------------------------------------------------------------!
! u(c,p,d) = 1/(1-gam)*(c*e^(theta*d)*e^(eta*p))^(1-gam)

real(8) :: 	beta= dexp(-.05/tlen),&	!People are impatient (5% annual discount rate to start)
		gam	= 2. , & !1.5, &	!IES
		eta 	= 0. , & !-0.185, &	!Util cost of participation
		theta 	= 0. !-0.448	!Util cost of disability

integer, parameter :: oldN = 4,&	!4!Number of old periods
		TT = oldN+2,&		!Total number of periods, oldN periods plus young and retired
		itlen = 12,&		! just an integer version of tlen so I don't have to keep casting
		nopt_tgts = 12,&		! number of calibration targets in main program
		nopt_pars = 12
!----------------------------------------------------------------------------!

!**Programming Parameters***********************!
integer, parameter ::	nal = 5,  &!6		!Number of individual alpha types
			ntr = 5,    &!7	        !Number of occupation trend points
			ndi = 2,    &		    !Number of individual disability risk types
			nl	= 1,    &			!Number of finding/separation rates
			nd  = 3,    &		    !Number of disability extents
			ne  = 4,    &!5	        !Points on earnings grid - should be 1 if hearnlw = .true.
			na  = 25,   &!50	    !Points on assets grid
			nz  = 2,    &		    !Number of aggregate shock states
			nj  = 16,   &!16		!Number of occupations
			Nskill  = 2,&			!number of skills that define occupations. First is always physical
			NKpolyT = 1,&			!polynomial order for time trend for occupation or number of spline segments
			NTpolyT = 2,& 			!polynomial order for time trend overall or number of spline segments
			Nknots   = 5,& 			! Number of knots
			max_DFBOLS = 150, &		!Iterations on DFBOLS
			maxiter = 800, &		!Tolerance parameter
			Nsim = 40000,&   !5000*nj	!how many agents to draw
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
					  wglev_0	 = .true. ,&  	!should the initial wage level be 0 for all occupations
					  noDI       = .false. ,&   !cannot get DI, even if they want it
					  nu_byD     = .true.  ,&   ! use different xiz for different health levels??
					  readshocks = .false.	,&	!readshocks from disk?
					  fndsep_seq = .true. , &   ! take finding and separation rates from the data's sequence?
					  allow_wtr_tmean = .false., &
					  nomeantrend = .true. !turn off mean wage trend


! these relate to what's changing over the simulation/across occupation
logical           ::  del_by_occ = .true.,& !delta is fully determined by occupation, right now alternative is fully random
					  j_regimes  = .true.,& !different pref shifts
					  j_rand     = .true.,&! randomly assign j, or let choose.
					  demog_dat	 = .true.,& !do the demographics follow
					  wtr_by_occ = .true.,& ! do we feed in occupation-specific trends for wages
					  occ_dat    = .true.,& ! do we use the occupational sizes?
					  NBER_tseq  = .true.,&	!just feed in NBER recessions?
					  RAS_pid    = .true.,& !balance the health transition matrix
					  buscyc	 = .true.,& !turn on/off business cycles?
					  welfare_cf = .false.  !compute counter-factual w/o DI for welfare comparisons?

logical			  ::  run_experiments = .false., &
					  run_cal = .false., &
					  refine_cal = .false., &
					  elast_xi 	 = .false.   !calculate the elasticity of w.r.t. xi?


real(8), parameter ::  amax 	 = 18.0,   &	!Max on Asset Grid
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
		wtr_tmean_ts(Tsim),&   !average trend among employed
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
		fnd_seq(Tsim+1),&	!sequence of actual finding rates
		sep_seq(Tsim+1),&	!sequence of actual separations rates
!		
		occ_onet(nj,Nskill),&!physical and 3 KSA
		occwg_datcoef_sqr(Nskill+1,NKpolyT+1),& !Only used with NKpolyT>=2. Then these are the data coefficients for wage regression. also includes 0-order and time-only
		occwg_datcoef(Nskill*2+NTpolyT+5),& !!coefficients for wage regression. First is cubic in time, then linear in skill dimension. Then 2 for age profile, 2 for health dummies, 1 const
		occwg_datsplncoef(Nskill+Nknots-1+Nskill*(Nknots-1)+5),& !!coefficients for wage spline regression. First is levels for skills, then cubic spline in time. Then spline for each skill. Then 2 for age profile, 2 for health dummies, 1 const
		occwg_dattrend(Tsim,nj),& !trend in occupation wage
		occwg_datlev(nj),&		!level of occupation wage
		occsz0(nj), &		   !Fraction in each occupation
		occpr_trend(Tsim,nj)   !trend in occupation choice

real(8), allocatable :: wage_coef(:), & !occupation-specific differences in wage-level
						tr_knots(:)  !will be the knot points

integer :: 	dgrid(nd)	! just enumerate the d states
real(8)	::	agegrid(TT)		! the mid points of the ages

!***preferences and technologies that may change
real(8) :: 	nu = 1.e-3, &		!Psychic cost of applying for DI - proportion of potential wage
		nud(nd) = 0. ,&		!psychic cost of applying for DI, if it depends on nu
		nuage   = 0. ,&		!psychic cost of applying for DI, if it depends on age
		util_const = 0.,&	!Give life some value
		Fd(nd,2) = 0.,&			! Fixed cost of participating in labor
!
!	Idiosyncratic income process
		alfrho(nd) = 0.988, &	!Peristence of Alpha_i type
		alfmu(nd) = 0.0,&		!Mean of Alpha_i type
		alfcondsig(nd) = 0.015**0.5,&	!Conditional StdDev of Alpha_i type (Normal)
!
		LTUrr = 0.35,&			!Home production income
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
		proc_time1 = 13.5,&!The average time to decision	(could be 2.5 for 'meets criteria' or 3.64 for initial decision)
		proc_time2 = 13.5,&!The average time to decision	(could be 28.05 for appeal that continues)
		xizcoef    = 0.1, &	!change in acceptance rate with z deterioration
		xid1coef = -0.1, &	!change in acceptance rate with z deterioration if d=1 or 2
		xid2coef = 0.1, &	!change in acceptance rate with z deterioration if d=1 or 2
		xid3coef  = 0.2, &	!change in acceptance rate with z deterioration if d= 3
		xiagezcoef = 0.279,&!OLD/LOWEDU coefficient from Hu, Lahiri, Vaughan & Wixon
		xi_d1shift = -0.,&	!worse hlth stage acceptance for d=1
		xi_d3shift = 0.,&	!better hlth stage acceptance for d=3
		xidscale   = 1.,& 	!scaling term for allowance probability
		wtr_scale  = .01,&	!scaling for wtr in the xi function
		xiVtrend   = 1.0,&  ! trend in vocational awards
!
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



logical  :: cal_on_grad = .false.

!remove this!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
real(8) :: tbase_out(Tsim, Nknots-1)

!**** calibration targets
real(8) :: apprt_target = .01,&	!target for application rates (to be filled below)
		dirt_target = 0.018,&	!target for di rates
		diaward_target = 0.00337,& !target for new award rate (3-yr average)
		d1_diawardfrac_target = 0.16,&
		d2_diawardfrac_target = 0.22,&
		d3_diawardfrac_target = 0.62,&
		voc_acc_target = 0.25,&		!fraction of admissions from vocational criteria, target 1985
		hlth_acc_target = 0.75,&		!fraction taken based on health criteria, target 1985
		hlth_acc_end_target = 0.4, &	!fraction taken based on health criteria, target 2013
		old_target = 0.46,&		!fraction over 55
		avg_unrt = 0.056,&	!average rate of unemployment over the period.
		avg_undur = 3.,&	! average months of unemployment
		avg_frt   = .3242085 ,&	! average rate of job finding from U
		award_age_target  = 0.8/0.3,&	!target for increase in vocation due to age (from Chen & van der Klaauw page 771 ) In levels it's (0.093+0.287)
		p1d1_2545target = 1. ,&	! normalize how much young healthy participate
		p1d2_2545target = 0.7562,&	! how much less d=2 participate: (.927-.226)/.927
		p1d3_2545target = 0.1737,&	! how much less d=3 participate: (.927-.766)/.927
		p1d1_4665target = 0.9364,&	! how much less d=2 participate: (.927-.059)/.927
		p1d2_4665target = 0.5383,&	! how much less d=2 participate: (.927-.428)/.927
		p1d3_4665target = 0.0831,&	! how much less d=3 participate: (.927-.850)/.927
		allowrt_target  = 0.677		! average allowance rate (finally allowed from Autor et al Delay Decay)

real(8) :: p1d1_target, p1d2_target, p1d3_target

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
	real(8) :: wage_coef_O2_read(17),wage_coef_O3_read(21),wage_coef_O1_read(22),wage_coef_CS_read(25), realwagegrowth(36)

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
			close(fread)
		else
			open(unit= fread, file = "OLSWageTrend_CS1.csv")
			do j=1,21 !(7 + Nskill+(Nknots-1)*(Nskill+1))
				read(fread,*) wage_coef_CS_read(j)
			enddo
			close(fread)
		endif

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


	!read the coefficients in:
	occwg_datcoef = 0._dp
	occwg_datsplncoef = 0._dp
	occwg_datcoef_sqr   = 0._dp
	if( tr_spline .eqv. .false. ) then
		if(wglev_0 .eqv. .false.) allocate(wage_coef(NTpolyT+Nskill*2+5))
		if(wglev_0 .eqv. .true.) allocate(wage_coef(NTpolyT+Nskill+5))
		wage_coef = 0._dp
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
			!wage_coef = occwg_datcoef_sqr
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
		wage_coef = occwg_datcoef
	else
	!	if(wglev_0 .eqv. .true.) allocate(wage_coef((Nskill+1)*(Nknots-1)+5) )
	!	if(wglev_0 .eqv. .false.)
		allocate(wage_coef(Nknots-1 + Nskill*Nknots+5) )
		wage_coef = 0._dp

		t= 6
		do j=1,Nskill
			occwg_datsplncoef(j) = wage_coef_CS_read(t)
			t = t+1
		enddo
		! Allow an overall trend or not?
		do j=1,(Nknots-1)
			if(allow_wtr_tmean .eqv. .true.) then
				occwg_datsplncoef(j+Nskill) = wage_coef_CS_read(t)
			else
				occwg_datsplncoef(j+Nskill) = 0._dp
			endif
			t = t+1
		enddo
		do k=1,NSkill
			do j=1,(Nknots-1)
				occwg_datsplncoef(j +(k-1)*(Nknots-1) + Nskill+(Nknots-1)) = wage_coef_CS_read(t)
				t =t+1
			enddo
		enddo


		wage_coef = occwg_datsplncoef
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
					occwg_dattrend(t,j) = occwg_datsplncoef(k)*occ_onet(j,k) +occwg_dattrend(t,j)
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

				if(nomeantrend .eqv. .false. )then
					do i=1,(Nknots-1)
						occwg_dattrend(t,j) = occwg_datsplncoef(i+Nskill)*tbase(i) +occwg_dattrend(t,j)
					enddo
				endif
				do k=1,Nskill
					do i=1,(Nknots-1)
						if (k>1) then
							occwg_dattrend(t,j) = occwg_datsplncoef(i+(k-1)*(Nknots-1)+Nskill+Nknots-1)*tbase(i)*occ_onet(j,k) &
									& + occwg_dattrend(t,j)
						else
							occwg_dattrend(t,j) = occwg_datsplncoef(i+(k-1)*(Nknots-1)+Nskill+Nknots-1)*tbase(i)*occ_onet(j,k) &
									& + occwg_dattrend(t,j)
						endif
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

	!deflate by real wage growth (series LEU0252881600A_PC1, median wage of hourly and salaried workers from the CPS)
	realwagegrowth = (/ 0., 1.91083, 2.1875, 0.61162, -0.91185, -1.22699, -2.17391, -0.63492, 0.31949, 0.95541, -0.63091, -0.31746, -0.31847, 0.31949, 2.2293, 2.49221, 1.51976, 0.5988, 0.59524, -0.29586, 0.29674, -1.47929, 0., 0.6006, 0., 2.98507, -0.86957, -1.75439, -0.29762, -0.59701, 0.3003, 2.09581, 1.75953, 1.15274, 0.5698, 1.69972  /)
	realwagegrowth(1) = exp(realwagegrowth(1))
	do t=2,size(realwagegrowth)
		realwagegrowth(t) = exp(realwagegrowth(t)/100.)*realwagegrowth(t-1)
	enddo
	!do i=1,nj
	!	do t=1,Tsim
	!		if(t>TossYears*itlen) then
	!			occwg_dattrend(t,i) = occwg_dattrend(t,i) - log(realwagegrowth(t/itlen- TossYears+1))
	!		endif
	!	enddo
	!enddo

	!initialize the input to the observed
	if(wglev_0 .eqv. .true.) then
		wage_lev = 0._dp
		do i=1,nj
			occwg_datlev(i) = occwg_dattrend(TossYears*itlen,i)
			occwg_dattrend(:,i) = occwg_dattrend(:,i) - occwg_datlev(i)
		enddo
	else
		occwg_datlev = 0.
		wage_lev = 0._dp
	endif
	wage_trend = occwg_dattrend


	!Wage-trend grid-setup
	trgrid(1) = (minval(wage_trend(Tsim,:) + wage_lev ))*1.1_dp
	do t=1,Tsim
		do i=1,nj
			if(trgrid(1)> wage_trend(t,i) + wage_lev(i) ) trgrid(1)= wage_trend(t,i) + wage_lev(i)
		enddo
	enddo

	if(ntr>1) then
		trgrid(ntr) = (maxval(wage_trend(Tsim,:)+wage_lev))*1.2_dp
		do t=1,Tsim
			do i=1,nj
				if(trgrid(ntr)< wage_trend(t,i) + wage_lev(i) ) trgrid(ntr)= wage_trend(t,i) + wage_lev(i)
			enddo
		enddo

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
			if(i==1) then
				! recession
				fndrate(i,j) = 0.489 ! UE_occ_read(k,j)
				seprisk(i,j) = 0.032 ! 0.0793*fndrate(i,j)/(1.-0.0793) ! EU_occ_read(k,j)
			else
				fndrate(i,j) = 0.578 ! UE_occ_read(k,j)
				seprisk(i,j) = 0.034 !0.05*fndrate(i,j)/(1.-0.05) ! EU_occ_read(k,j)
			endif
		enddo
	enddo

	do i=1,nz
		if(nl>1) then
			do j=1,nl

				fndgrid(j,i) = (maxval(fndrate(i,:))-minval(fndrate(i,:))) *dble( j-1 )/dble(nl-1) + minval(fndrate(i,:))
				sepgrid(j,i) = (maxval(seprisk(i,:))-minval(seprisk(i,:))) *dble( j-1 )/dble(nl-1) + minval(seprisk(i,:))

			enddo
		else
			fndgrid(1,i) = fndrate(i,1)
			sepgrid(1,i) = seprisk(i,1)
		endif
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

	agegrid = agegrid + 10._dp !started them at 30, but they'veeen working since 20

	!from Mincer regressions (see Appendix) with Heckman 2-step
	wtau(1) =  0.0
	!wtau(2) =  0.0757/(.5*agegrid(2)+.5*agegrid(3)-agegrid(1))*(agegrid(2)-agegrid(1))
	!wtau(3) =  0.0757/(.5*agegrid(2)+.5*agegrid(3)-agegrid(1))*(agegrid(3)-agegrid(1))
	!wtau(4) =  0.0157
	!wtau(5) = -0.0661
	do i=1,(TT-1)
		wtau(i)  = 0.0213453024627902*(agegrid(i)-10) - 0.000765398276859*(agegrid(i)-10)**2
	enddo
	do i=(TT-1),1,(-1)
		wtau(i)  = wtau(i) -wtau(1)
	enddo

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
	wd(2) = -0.074520650705695	!Partially Disabled, small penalty
	wd(3) = -0.195639344387911	!Full Disabled, large penalty
	!Fixed cost of particpation
	Fd = 0.
	!Fd(2,1) = 0.276*5
	!Fd(3,1) = 0.524*5
	!Fd(2,2) = 0.276*5
	!Fd(3,2) = 0.524*5


	!DI Acceptance probability for each d,t status
	!xi_d1shift = (1491*.603+2211*0.546)/(1491+2211) - (347*.581+752*.655)/(347+752) !differences from Lahiri, Vaughn, Wixon : denial rate for hlth 1,2 vs. 3,4
	!xi_d3shift = (1491*.603+2211*0.546)/(1491+2211) - .484

	! initialize a few starting values
	xi_d(1) = 0.001!xi_d(2)+xi_d1shift	!
	!xi_d(2) = xi_d(1)-xi_d1shift !xi_d(3)-xi_d3shift
	!xi_d(3) = xi_d(2)+xi_d3shift !

	xi_d(2)= 0.2329423
	xi_d(3)= 0.3083954

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
	! do i=1,nz
	! 	do j=1,nl
	! 		if(nl>1) then
	! 			fndgrid(j,i) = (maxval(fndrate(i,:))-minval(fndrate(i,:))) *dble( j-1 )/dble(nl-1) + minval(fndrate(i,:))
	! 			sepgrid(j,i) = (maxval(seprisk(i,:))-minval(seprisk(i,:))) *dble( j-1 )/dble(nl-1) + minval(seprisk(i,:))
	! 		else
	! 			fndgrid(1,i) = fndrate(i,1)
	! 			sepgrid(1,i) = seprisk(i,1)
	! 		endif
	! 	enddo
	! enddo

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


! fills arrays fndrt_st and seprt_ts
subroutine fndsep_data( fndrt_ts, seprt_ts)

	integer :: it, itr, Tperiods, Tdata
	real(8) :: avg_f, avg_s
	real(8) :: fndrt_ts(:), seprt_ts(:)
	real(8), dimension(361) :: data_frt = &
	!& detrened dimension(370):	(/ 0.434148970867114,0.419868253433264,0.417127880131874,0.440082057505792,0.468106690755282,0.405274919693253,0.480975221303461,0.461209912439115,0.460363421331156,0.493529238919738,0.434839512121959,0.501952681287541,0.490484529960741,0.484862366919243,0.471735294444835,0.498795941689961,0.462705391185101,0.474627313014805,0.522614441948001,0.485653793421429,0.467406402980975,0.517633978567253,0.480287429216672,0.543368003498147,0.409885859089735,0.494069418228195,0.50818915536782,0.488187079018511,0.463837154864789,0.501857860324931,0.516748379110817,0.451194733475091,0.48886521749697,0.483408092119093,0.522457552292026,0.506887603360859,0.499539627731183,0.511106908801297,0.532763144838602,0.51958052960893,0.530767158942067,0.551086809876656,0.529035711252921,0.576752243237703,0.5274991318769,0.57664752724986,0.577635067694792,0.535227184459347,0.536739688795369,0.543501174601163,0.629040871874663,0.535629709866052,0.629670049022239,0.549098808683172,0.540840547471729,0.63810685224787,0.600020142289108,0.601285653236617,0.566562568936927,0.566832547084877,0.702882079269883,0.654324658379769,0.570341009905785,0.622686781733519,0.623416120973217,0.640412891492334,0.598857611977871,0.614916024858344,0.616725945708865,0.603975294576115,0.643596793355234,0.586108362295765,0.637759701328218,0.633602125629276,0.586646738405046,0.589595033287617,0.634500548078047,0.538631637313257,0.544005189311519,0.526751817195876,0.545988151462662,0.508811159037405,0.522533279090492,0.516783025233333,0.502194774740693,0.457770502010148,0.497816777685747,0.465934555840632,0.472166041534478,0.486653085451285,0.463953743435873,0.468997612430869,0.4522384324813,0.447439894310074,0.41243071473571,0.414680901366309,0.375810306151403,0.415525225427736,0.404845310003859,0.367006262851398,0.371355603993997,0.418647529396413,0.409457503729849,0.403082841012528,0.453018333683671,0.379903855969588,0.378850201663036,0.433572156375388,0.461118977187969,0.43328100994551,0.425941812012141,0.427739949463591,0.412349011671066,0.455571858224875,0.443852407616723,0.425331004288582,0.451534578681574,0.448469602749475,0.46669790918234,0.419137100012794,0.426301164490686,0.471837744347915,0.497991910823474,0.52143013931811,0.493516583407679,0.506890230335565,0.476020625386945,0.516772321537307,0.452691066669863,0.545404967804095,0.535052386650851,0.531562441167884,0.545139005783454,0.526137904140962,0.450322742087536,0.532427995464304,0.542754495002531,0.488652504778885,0.539865832751548,0.57011554688734,0.568921173970315,0.549716093850053,0.52311280008088,0.524237934436379,0.572344028185759,0.514892557882938,0.476961605084469,0.551361951773832,0.572178936552446,0.502783822915274,0.597700213584303,0.518268050464634,0.506380717332078,0.591356514994818,0.530124700652908,0.594996777842952,0.53614824541836,0.568888294862495,0.539307577452876,0.603817325625468,0.544552081106339,0.555504598030087,0.584726574284316,0.549284501409509,0.641934770159542,0.624521159419803,0.57208998233013,0.6478396318876,0.644713131997715,0.676992689419396,0.792454035623954,0.684009686143028,0.613818244294894,0.643405925838384,0.67767911950858,0.63724379098181,0.737158193266817,0.675200047766651,0.697075749867245,0.65728263801863,0.643720247408573,0.73922991828563,0.73230962439044,0.733000619739146,0.649560978388891,0.698870072576246,0.737076776802119,0.705349741178399,0.736285031317194,0.746241182783252,0.794261821810841,0.719139094050223,0.701694140387808,0.851299046640755,0.805979387440452,0.681254555820518,0.776630900952006,0.689235031916285,0.690532428576437,0.800372172380119,0.762949951304335,0.714322694417744,0.720495996652179,0.653190099063202,0.779452531790232,0.704914386154815,0.78471042461693,0.714895744363002,0.665674285551345,0.620128404476588,0.612276516916777,0.613014985691863,0.568445443943641,0.564807652869401,0.535364281290891,0.585011888907127,0.554302599325346,0.561929347120111,0.480459626079801,0.534515623259022,0.482906377620008,0.523741730072412,0.549435321673332,0.513985019365197,0.495285526566549,0.502720133115737,0.480535226275818,0.521574291883958,0.47669644049568,0.501280869019331,0.443681172041693,0.493027259640704,0.435290692733977,0.476242506029289,0.485029978183786,0.44848071500299,0.486108407702263,0.481898631565112,0.49016316158153,0.478508912496059,0.466630037914246,0.427956304685269,0.561122827350721,0.501313592822387,0.481906993826952,0.554145234775236,0.528126072800506,0.547064558003159,0.505371500579048,0.526579722279934,0.541822505862964,0.541525526457686,0.509398043939271,0.521032030350452,0.554839226376927,0.568415339082667,0.569808056111195,0.556232301445988,0.555187574353988,0.550919248555709,0.590525390591956,0.579332625708148,0.59809844905575,0.598679764693333,0.548075583203193,0.615058436365429,0.59059542918373,0.599875637755072,0.620505865063057,0.589793881455673,0.592418717287069,0.649796407987504,0.646226299773317,0.558165760567428,0.635689629084961,0.544043837275735,0.623590055416451,0.558671508498143,0.549846186195394,0.591871944362506,0.555533277174797,0.539239373841211,0.596330247198962,0.561481988876863,0.536535980925744,0.580100148562814,0.512130836241452,0.553413438304895,0.58890884030567,0.519058370462914,0.528927870468525,0.540579811627551,0.485630060724597,0.469600489702185,0.485606820211495,0.469316388741029,0.414543462966319,0.446538501990296,0.377588011687704,0.432644415545281,0.343267440223417,0.336600218440503,0.327208981164239,0.29005840204399,0.314278791446328,0.348663690607473,0.275569124466483,0.28736951422546,0.279426786994625,0.297435649817518,0.301291168779797,0.341821135640261,0.283704906860023,0.277343368696956,0.269160608844686,0.322269111095654,0.31461045212038,0.298644644708612,0.275962724184074,0.306706253443921,0.29115283161843,0.264110571249333,0.344438455668298,0.324525509402368,0.295883045142908,0.297041927978628,0.298095776448753,0.310985681279847,0.322807308867574,0.316737246447369,0.299557625077266,0.303031757976876,0.336679295959864,0.326886248511392,0.340026186484514,0.330507950842581,0.322380618766794,0.346065371600795,0.335411706996131,0.320140277493182,0.342066989480064,0.338384041996506,0.370622114592696,0.363926426174905,0.346673555665058,0.359041064801395,0.334562661051898,0.343208794795843,0.406663947370872,0.360925676255309,0.334455721252143,0.377931581367271,0.358760489692887,0.390711315299493,0.359193210362924,0.373645207542555,0.421384715730625,0.393223593068276,0.398572616058768,0.412026958490518,0.357999718428955,0.377298655666592,0.473121566791212,0.410427076751515,0.456872328494582,0.425796611145121,0.446397683907905,0.452117690728835,0.483604794547664,0.434148970866854 /)
	& (/ 0.499931949128846,0.485294552208443,0.4821974994205,0.504794997307864,0.532462951070801,0.469274500522219,0.544618122645874,0.524496134294975,0.523292963700463,0.556102101802492,0.49705569551816,0.563812185197188,0.551987354383835,0.546008511855784,0.532524759894823,0.559228727653396,0.522781497661983,0.534346740005133,0.581977189451776,0.544659861438651,0.526055791511644,0.575926687611369,0.538223458774235,0.600947353569157,0.467108529674191,0.550935409326098,0.56469846697917,0.544339711143308,0.519633107503033,0.557297133476622,0.571830972775955,0.505920647653675,0.543234452189001,0.537420647324571,0.576113428010951,0.560186799593231,0.552482144477002,0.563692746060562,0.584992302611314,0.571453007895089,0.582282957741673,0.602245929189709,0.579838151079421,0.62719800357765,0.577588212730293,0.6263799286167,0.627010789575079,0.584246226853081,0.58540205170255,0.591806858021791,0.676989875808738,0.583222034313573,0.676905693983207,0.595977774157587,0.587362833459591,0.684272458749179,0.645829069303864,0.64673790076482,0.611658136978576,0.611571435639973,0.747264288338426,0.698350187961759,0.614009860001222,0.665998952342403,0.666371612095547,0.683011703128111,0.641099744127095,0.656801477521015,0.658254718884983,0.64514738826568,0.684412207558246,0.626567097012224,0.677861756558123,0.673347501372628,0.626035434661845,0.628627050057863,0.67317588536174,0.576950295110397,0.581967167622105,0.564357116019909,0.583236770800142,0.545703098888332,0.559068539454866,0.552961606111154,0.538016676131961,0.493235723914862,0.532925320103908,0.50068641877224,0.506561224979533,0.520691589409787,0.497635567907822,0.502322757416264,0.485206897980142,0.480051680322363,0.444685821261446,0.446579328405492,0.407352053704033,0.446710293493813,0.435673698583383,0.397477971944368,0.401470633600414,0.448405879516277,0.43885917436316,0.432127832159286,0.481706645343876,0.408235488143239,0.406825154350134,0.461190429575933,0.488380570901961,0.460185924172949,0.452490046753027,0.453931504717924,0.438183887438846,0.481050054506101,0.468973924411396,0.450095841596702,0.475942736503141,0.472521081084489,0.490392708030801,0.442475219374701,0.44928260436604,0.494462504736716,0.520259991725722,0.543341540733805,0.515071305336821,0.528088272778154,0.49686198834298,0.537257005006789,0.472819070652792,0.565176292300471,0.554467031660674,0.550620406691154,0.56384029182017,0.544482510691125,0.468310669151146,0.550059243041361,0.560029063093035,0.505570393382836,0.556427041868946,0.586320076518184,0.584769024114606,0.565207264507791,0.538247291252065,0.539015746121011,0.586765160383838,0.528957010594463,0.490669378309441,0.564713045512251,0.585173350804312,0.515421557680587,0.609981268863063,0.530192426256841,0.517948413637732,0.602567531813918,0.540979037985455,0.605494435688946,0.546289223777801,0.578672593735383,0.548735196839211,0.61288826552525,0.553266341519567,0.563862178956762,0.592727475724438,0.556928723363078,0.649222312626558,0.631452022400266,0.57866416582404,0.654057135894956,0.650573956518518,0.682496834453646,0.797601501171651,0.688800472204172,0.618252350869485,0.647483352926422,0.681399867110064,0.640607859096741,0.740165581895195,0.677850756908476,0.699369779522517,0.659219988187348,0.645300918090738,0.740453909481242,0.733176936099499,0.733511251961652,0.649714931124844,0.698667345825646,0.736517370564966,0.704433655454692,0.735012266106934,0.744611738086439,0.792275697627475,0.716796290380304,0.698994657231336,0.84824288399773,0.802566545310873,0.677485034204386,0.772504699849321,0.684752151327047,0.685692868500646,0.795175932817775,0.757397032255438,0.708413095882293,0.714229718630175,0.646567141554645,0.772472894795122,0.697578069673152,0.777017428648714,0.706846068908232,0.657267930610022,0.611365370048712,0.603156803002348,0.603538592290881,0.558612371056106,0.554617900495313,0.524817849430249,0.574108777559932,0.543042808491598,0.55031287679981,0.468486476272947,0.522185793965615,0.470219868840048,0.510698541805898,0.536035453920265,0.500228472125577,0.481172299840376,0.488250226903011,0.465708640576539,0.506391026698125,0.461156495823294,0.485384244860392,0.427427868396201,0.476417276508659,0.418324030115379,0.458919163924138,0.467349956592081,0.430444013924732,0.467715027137452,0.463148571513748,0.471056422043613,0.459045493471589,0.446809939403223,0.407779526687692,0.540589369866591,0.480423455851704,0.460660177369716,0.532541738831447,0.506165897370164,0.524747703086264,0.482697966175599,0.503549508389932,0.518435612486409,0.517781953594578,0.48529779158961,0.496575098514238,0.530025615054159,0.543245048273346,0.544281085815321,0.530348651663561,0.528947245085008,0.524322239800176,0.56357170234987,0.552022257979509,0.570431401840557,0.570656037991587,0.519695177014894,0.586321350690577,0.561501664022325,0.570425193107114,0.590698740928546,0.559630077834608,0.561898234179451,0.618919245393333,0.614992457692593,0.526575239000151,0.603742428031131,0.511739956735351,0.590929495389514,0.525654268984653,0.516472267195351,0.55814134587591,0.521445999201648,0.504795416381509,0.561529610252706,0.526324672444054,0.501021985006382,0.544229473156899,0.475903481348984,0.516829403925874,0.551968126440095,0.481760977110786,0.491273797629844,0.502569059302317,0.44726262891281,0.430876378403845,0.446526029426602,0.429878918469582,0.374749313208319,0.406387672745743,0.337080502956598,0.391780227327622,0.302046572519205,0.295022671249738,0.285274754486921,0.247767495880118,0.271631205795903,0.305659425470495,0.232208179842952,0.243651890115376,0.235352483397988,0.253004666734327,0.256503506210053,0.296676793583964,0.238961877383637,0.219884125557898,0.238203885317173,0.231485667667553,0.22294622832873,0.275698051093145,0.267682712631317,0.251360225732996,0.228321625721905,0.258708475495199,0.242798374183155,0.215399434327505,0.295370639259917,0.275101013507433,0.24610186976142,0.246904073110587,0.247601242094159,0.2601344674387,0.271599415539874,0.265172673633116,0.247636372776459,0.250753826189516,0.284044684685951,0.273894957750926,0.286678216237495,0.276803301109009,0.268319289546669,0.291647362894116,0.280637018802899,0.265008909813397,0.286578942313726,0.282539315343615,0.314420708453252,0.307368340548907,0.289758790552507,0.301769620202291,0.276934536966241,0.285223991223633,0.348322464312109,0.302227513709993,0.275400879220273,0.318520059848848,0.298992288687911,0.330586434807964,0.298711650384842,0.31280696807792,0.360189796779436,0.331671994630534,0.336664338134473 /)
	real(8), dimension(361) :: data_srt = &
	!& detrended dimension(370) 	(/ 0.032101489603528,0.032545012893885,0.030937450950555,0.029482685514105,0.031343547652577,0.032210979176577,0.035879829005021,0.031465543943862,0.034435618719281,0.033431702815913,0.031864493217095,0.036352936713522,0.034037450404194,0.034408971283663,0.034826298332642,0.034741784074934,0.035568535428875,0.034933333858088,0.034430782424248,0.033950836051939,0.032575359694203,0.034989799971026,0.032975078699234,0.033371663331831,0.032933514716671,0.03534777617956,0.035024396307915,0.035720702896765,0.033058287752522,0.033163882509325,0.034437132257279,0.031769918939131,0.033816032313564,0.032010016705454,0.031405004887245,0.032651492149695,0.032153440866535,0.032991993940142,0.03022962180661,0.03167420510597,0.031016681595715,0.031708339022272,0.029632582594231,0.032014508537329,0.031121218854664,0.030601084050615,0.030729626595657,0.028938065113905,0.029054329102123,0.029487838812159,0.030233110616691,0.029936374021261,0.031046009019306,0.027921225744379,0.030317441901062,0.0316119705072,0.030905988420742,0.029545904962518,0.028365440101326,0.029841256234548,0.033827739215,0.029601089575594,0.029497114401367,0.030896725996244,0.032433160053203,0.031202269662531,0.029670077792518,0.032040798062399,0.03143703794173,0.0322069860692,0.033757695438583,0.030502124671857,0.031996835477474,0.031039083371484,0.032100915444419,0.030805234804876,0.030443558948963,0.030894115587404,0.031770070575607,0.032041992093656,0.031722069291408,0.033624475867288,0.033615748915412,0.033885375654638,0.034939459123752,0.032999094646318,0.032906861587112,0.034206665204313,0.033067774582851,0.032759951523492,0.033311150604636,0.032892775830381,0.033049705087497,0.031889386596238,0.033339662139175,0.03101672524925,0.029307895014445,0.031640200090768,0.030804607704246,0.030458339649007,0.031774474938876,0.03269325586435,0.031424605931212,0.031763225046765,0.031729136318538,0.02979546246248,0.028870239591376,0.03189091871021,0.032139445405185,0.030304685457131,0.031962102041823,0.031283538091962,0.028777048698451,0.03151570844053,0.030121441718528,0.028264404534489,0.032370424630038,0.028613042358948,0.030193554050719,0.029065884003636,0.028756675136286,0.030617537274758,0.031875558035971,0.030082476635777,0.030654596947238,0.031546585584885,0.028193829444853,0.030230525285412,0.025673963937174,0.029404802596184,0.028964716830069,0.030962639849202,0.028202920122363,0.028691951559006,0.029764423061034,0.028772606900555,0.030972738677088,0.029157573756005,0.031476885019766,0.031907027067372,0.031205694143206,0.032285139388672,0.02995886463864,0.030050839734112,0.031508297636744,0.028953945713478,0.028129380186239,0.031760768495442,0.028764954770647,0.029941734401032,0.028171362370497,0.02831552133247,0.026903403019358,0.034445720293175,0.029465906362387,0.031764654023701,0.027863177523867,0.030466147292114,0.027481837263762,0.028897363298794,0.028758400198084,0.027192218779388,0.028084460736776,0.028471657637821,0.029602035052248,0.028803397545486,0.028443003189436,0.029995798721386,0.03060415865687,0.0336582461791,0.032484188817221,0.031789799834171,0.02929678109038,0.029976338774048,0.031620895656629,0.031183343087917,0.033789203681436,0.030131390969444,0.031887110514532,0.028638722107502,0.030223465729051,0.030915301969543,0.033522235446668,0.031393550094066,0.029899624406,0.031415763081413,0.03166569048554,0.03099925608089,0.03089224106978,0.032042973470397,0.03259602029553,0.030167970614543,0.030940064509225,0.035089059553657,0.030810792291875,0.030148321596103,0.032721042007127,0.029111062022567,0.03065515527638,0.031604304768367,0.031456102947717,0.02951457617971,0.029796988115239,0.031394427526126,0.034770618548722,0.033004061215453,0.037437042589917,0.032097157000041,0.033515037418542,0.031362748626609,0.034516524917226,0.033747302302661,0.034860722195649,0.035155531897292,0.034674765961118,0.03611675645843,0.034285319469323,0.034771132255149,0.032565989044461,0.032956866532321,0.030598375471787,0.033137061630645,0.033352239526868,0.032022047797999,0.030916565579934,0.034111260459189,0.032472346472707,0.031545350363421,0.031754716806779,0.032497038181395,0.030232662854962,0.033923119367137,0.032021110478616,0.03175614315885,0.031833813139577,0.030304869722705,0.031421210243287,0.029306144552016,0.030136472184643,0.030274760089106,0.028250919306664,0.029422900763144,0.033164617269826,0.031212778898272,0.030086762880082,0.033012104369597,0.030945676068463,0.032852493388056,0.031848994717202,0.030934247362418,0.032631315937046,0.031268580072131,0.031641144032883,0.028808587430129,0.032296373248906,0.031682942274622,0.031182298189415,0.031275465377783,0.02987619186562,0.031839424981493,0.033163159863571,0.032601689238371,0.032217751047634,0.03025029625915,0.030679454830117,0.031889821677832,0.031471879037251,0.030564731274676,0.032354842699799,0.032289580742354,0.031676939529614,0.031685751244463,0.031617419647812,0.029701358964285,0.031187516503376,0.030455069954638,0.031331710396255,0.027698216379332,0.029461714330563,0.029300637906641,0.03114624125605,0.028815674971594,0.031597667290424,0.031311197595056,0.029323939163345,0.031501212558607,0.032184346520962,0.03195524629619,0.031078750222191,0.033159911127752,0.029945184443859,0.036825666643251,0.031843826154541,0.034069811445745,0.03657775650906,0.034378916489557,0.035126661886848,0.037902961962103,0.036666871066173,0.040599426744264,0.036925312161926,0.036863000634397,0.036119437819408,0.035988789566644,0.035002844463676,0.036373041166357,0.033487585441364,0.033118563312982,0.035083629271682,0.032922906422513,0.03432075968037,0.03492520061046,0.031793383987449,0.032367590461477,0.030502759970134,0.032487465659288,0.031972143930305,0.032224351212095,0.031129806446346,0.033395026494973,0.030512223176673,0.033651292021192,0.031987791118739,0.032077739656289,0.029851943436187,0.030935128873113,0.032306383558238,0.031386599610971,0.034788450433143,0.031980391339391,0.031239451128772,0.031594133696092,0.032377492868636,0.030779840197499,0.032633683889049,0.030163761945107,0.031312467787803,0.032197989321492,0.032186758974536,0.030837151032632,0.032810263618377,0.032495794641126,0.034104260911296,0.030608304612191,0.031910432810538,0.031701911089375,0.033103996277011,0.033265142562745,0.033961776977141,0.030304143079083,0.031256019483569,0.032641758097375,0.032010848061011,0.032050520080693,0.030123015305369,0.032176877283623,0.035899505402815,0.029895445256207,0.030441071303065,0.031925306592362,0.030261403358539,0.030757408538621,0.031057761472228,0.032472959681714,0.032126721428798,0.033131899286874,0.032346076778739,0.031008672887853,0.032097446260227,0.032101489603494 /)
	& (/ 0.04065814614296,0.041055291891098,0.039401352405548,0.03790020942688,0.039714694023132,0.040535748004913,0.044158220291138,0.039697557687759,0.04262125492096,0.041570961475372,0.039957374334335,0.044399440288544,0.042037576436997,0.042362719774246,0.042733669281006,0.042602777481079,0.043383151292801,0.042701572179794,0.042152643203735,0.041626319289208,0.040204465389252,0.042572528123856,0.040511429309845,0.040861636400223,0.040377110242844,0.042744994163513,0.042375236749649,0.04302516579628,0.040316373109818,0.040375590324402,0.041602462530136,0.038888871669769,0.040888607501984,0.039036214351654,0.038384824991226,0.039584934711456,0.039040505886078,0.039832681417465,0.037023931741715,0.038422137498856,0.037718236446381,0.038363516330719,0.036241382360458,0.038576930761337,0.037637263536453,0.037070751190186,0.037152916193009,0.035314977169037,0.035384863615036,0.035771995782852,0.036470890045166,0.036127775907517,0.037191033363342,0.034019872546196,0.03636971116066,0.037617862224579,0.036865502595902,0.035459041595459,0.034232199192047,0.035661637783051,0.039601743221283,0.035328716039658,0.035178363323212,0.03653159737587,0.03802165389061,0.036744385957718,0.035165816545487,0.037490159273148,0.03684002161026,0.037563592195511,0.039067924022675,0.03576597571373,0.037214308977127,0.036210179328919,0.037225633859634,0.035883575677872,0.035475522279739,0.035879701375961,0.036709278821945,0.036934822797775,0.036568522453308,0.038424551486969,0.038369446992874,0.03859269618988,0.039600402116776,0.037613660097122,0.037475049495697,0.038728475570679,0.037543207406998,0.03718900680542,0.037693828344345,0.03722907602787,0.037339627742767,0.03613293170929,0.037536829710007,0.035167515277863,0.033412307500839,0.035698235034943,0.034816265106201,0.034423619508743,0.035693377256394,0.036565780639649,0.035250753164291,0.035542994737625,0.035462528467178,0.033482477068901,0.032510876655579,0.035485178232193,0.035687327384949,0.033806189894676,0.035417228937149,0.034692287445068,0.032139420509338,0.034831702709198,0.033391058444977,0.03148764371872,0.035547286272049,0.03174352645874,0.033277660608292,0.03210361301899,0.031748026609421,0.033562511205673,0.034774154424667,0.032934695482254,0.033460438251495,0.034306049346924,0.030906915664673,0.032897233963013,0.028294295072556,0.031978756189346,0.031492292881012,0.033443838357926,0.030637741088867,0.031080394983292,0.0321064889431,0.031068295240402,0.033222049474716,0.031360507011414,0.033633440732956,0.034017205238342,0.033269494771957,0.034302562475205,0.031929910182953,0.031975507736206,0.033386588096619,0.030785858631134,0.029914915561676,0.033499926328659,0.030457735061646,0.031588137149811,0.029771387577057,0.029869168996811,0.02841067314148,0.035906612873077,0.03088042140007,0.033132791519165,0.029184937477112,0.03174152970314,0.028710842132568,0.030079990625382,0.029894649982452,0.028282091021538,0.029127955436707,0.029468774795532,0.03055277466774,0.029707759618759,0.02930098772049,0.03080740571022,0.031369388103485,0.034377098083496,0.033156663179398,0.032415896654129,0.029876500368118,0.030509680509567,0.03210785984993,0.031623929738999,0.034183412790299,0.030479222536087,0.032188564538956,0.028893798589707,0.030432164669037,0.03107762336731,0.033638179302216,0.031463116407395,0.029922813177109,0.031392574310303,0.031596124172211,0.030883312225342,0.030729919672012,0.031834274530411,0.032340943813324,0.029866516590119,0.030592232942581,0.034694850444794,0.030370205640793,0.029661357402802,0.032187700271607,0.028531342744827,0.030029058456421,0.030931830406189,0.03073725104332,0.028749346733093,0.028985381126404,0.030536442995072,0.033866256475449,0.03205332159996,0.036439925432205,0.03105366230011,0.032425165176392,0.030226498842239,0.033333897590637,0.032518297433853,0.033585339784622,0.033833771944046,0.033306628465653,0.034702241420746,0.03282442688942,0.033263862133026,0.031012341380119,0.03135684132576,0.028951972723007,0.031444281339645,0.031613081693649,0.030236512422562,0.029084652662277,0.032232969999313,0.030547678470612,0.029574304819107,0.029737293720245,0.030433237552643,0.028122484683991,0.031766563653946,0.029818177223206,0.029506832361221,0.029538124799728,0.027962803840637,0.029032766819,0.02687132358551,0.027655273675919,0.027747184038162,0.025676965713501,0.026802569627762,0.030497908592224,0.028499692678452,0.027327299118042,0.030206263065338,0.028093457221985,0.029953896999359,0.028904020786285,0.027942895889282,0.029593586921692,0.028184473514557,0.02851065993309,0.025631725788117,0.029073134064674,0.028413325548172,0.027866303920746,0.027913093566895,0.026467442512512,0.028384298086166,0.029661655426025,0.029053807258606,0.02862349152565,0.026609659194946,0.026992440223694,0.02815642952919,0.02769210934639,0.026738584041596,0.0284823179245,0.028370678424835,0.027711659669876,0.027674093842507,0.027559384703636,0.02559694647789,0.027036726474762,0.026257902383804,0.027088165283203,0.02340829372406,0.025125414133072,0.024917960166931,0.026717185974121,0.024340242147446,0.027075856924057,0.02674300968647,0.02470937371254,0.026840269565582,0.027477025985718,0.027201548218727,0.026278674602509,0.028313457965851,0.025052353739739,0.031886458396912,0.026858240365982,0.029037848114967,0.031499415636063,0.029254198074341,0.029955565929413,0.032685488462448,0.0314030200243,0.035289198160172,0.031568706035614,0.031460016965866,0.030670076608658,0.030493050813675,0.029460728168488,0.030784547328949,0.027852714061737,0.027437314391136,0.029356002807617,0.027148902416229,0.028500378131867,0.029058441519737,0.02566921710968,0.023619890213013,0.026173345889456,0.025435162152151,0.023675705663228,0.028604433557126,0.027772820074331,0.02638588004932,0.02396746347357,0.026841718208111,0.026379424246063,0.022086160300395,0.029569557945712,0.027207792544691,0.02433974536102,0.024717569475317,0.024488034932389,0.026042064397054,0.02686148065779,0.026225868820858,0.024491509395474,0.024195544632322,0.026726305123623,0.02544379388943,0.025947973770678,0.025054170111284,0.023967518238374,0.026051289495989,0.025067794707884,0.023671819830826,0.02559855476005,0.024902812342582,0.026599582710796,0.026002961564875,0.024172726839158,0.025884690549382,0.024081264084021,0.023794417469361,0.028242361971252,0.024858540088701,0.02232980101786,0.025825950798555,0.023545239562263,0.025648947528204,0.023175903909169,0.024269506143977,0.026695054755941,0.02381781740648,0.023789985349974 /)

	Tperiods = size(fndrt_ts)
	Tdata    = size(data_frt)
	avg_s    = 0.0
	avg_f    = 0.0
	
	do it=1,Tdata
		itr = it + (Tperiods-Tdata)
		fndrt_ts(itr) = data_frt(it)
		seprt_ts(itr) = data_srt(it)
		avg_f = data_frt(it) + avg_f
		avg_s = data_srt(it) + avg_s
		
	enddo
	avg_f = avg_f / dble(Tdata)
	avg_s = avg_s / dble(Tdata)

	fndrt_ts(Tperiods) = avg_f
	seprt_ts(Tperiods) = avg_s
	do it = 1,(Tperiods-Tdata)
		fndrt_ts(itr) = avg_f 
		seprt_ts(itr) = avg_s
	enddo

end subroutine fndsep_data



end module V0para



! Desctiption of Hu, Lahiri Vaughan & Wixon paper calculations:
! Average effect of d=1 in stage 3 (conditional prob of getting past 2)
! 0.078*0.61+0.104*0.17+0.149*0.22+.001*140.6
! =
! 0.23864
!
! Average effect of d=2 in stage 3 (conditional prob of getting past 2)
! 0.109*0.24+0.248*0.03+0.161*0.22+0.211*0.06+.195*0.06
! =
! 0.09338
!
! Average effect of d=2 on stage 2
! 0.19*0.03+0.091*0.09+0.087*0.28+0.139*0.13+0.152*0.03+0.049*0.71
! =
! 0.09567
!
! Average effect of d=1 on stage 2
! 0.057*0.18+.107*.07+.127*.03+0.039*.24+.087*.2+.091*.11
! =
! 0.05833
!
!
! Overall probability of going from stage 2->stage 3
! 751/(751+176)
! =0.8101402
!
!
! Prob allowed @stage 3 & d=2:
! (0.8101402*(1+0.09567))*0.09338
! =
! 0.08288841
!
! Prob allowed @stage 3 & d=1:
! (0.8101402*(1-0.05833))*0.23864
! =
! 0.1820548
!
! Overall probability allowed:
! 264/(751+176)
! =
! 0.2847896
!
! D2, avg acceptance rate in T1 periods:
! 0.2847896*(1+0.08288841)
! =
! 0.3083954
!
! D1, avg acceptance rate in T1 periods:
! 0.2847896*(1-0.1820548)
! =
! 0.2329423
