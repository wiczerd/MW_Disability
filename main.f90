! V0main.f90

!************************************************************************************************!


!-----------------------------------------------------
!************************************************************************************************!
! compiler line: gfortran -fopenmp -ffree-line-length-none -g V0para.f90 V0main.f90 -lblas -llapack -lgomp -lnlopt -o V0main.out
!       	     ifort -mkl -qopenmp -parallel -O3 -xhost V0para.f90 V0main.f90 -lnlopt -o V0main.out
!       	     ifort -mkl -init=snan -init=array -g V0para.f90 V0main.f90 -lnlopt -o V0main_dbg.out
!				 ifort -mkl -qopenmp -parallel -O3 -xhost V0para.f90 V0main.f90 -L$HOME/Resources/lib -lnlopt -o V0main.out
! val grind line: valgrind --leak-check=yes --error-limit=no --track-origins=yes --log-file=V0valgrind.log ./V0main_dbg.out &


include './toms659.f90'

module helper_funs

	use V0para
	implicit none

	!**********************************************************!
	!Public Policy Functions
	!	1)UI(e)		Unemployment Insurance
	!	2)SSDI(e,a)	DI (can be asset tested)
	!	3)SSret(e)	Social Sec. Retirement
	! 	4)xifun(d,z,t)	Acceptance probability
	!Utility Functions
	!	5)u(c,d,w)	w=1 if working; 2 if not
	!Earnings Index
	!	6)Hearn(t,e,w)	t=age, e=past index, w=current wage
	!Wage Function
	!	7) wage(bi,ai,d,z,t) indiv(beta,alf),disability,tfp,age
	!Locate Function
	!	8)finder(xx,x)	xx is grid, x is point, returns higher bound
	!Writing Subroutines
	!	9)  mat2csv(A,fname,append)   A=matrix, fname=file, append={0,1}
	!	10)  mati2csv(A,fname,append)  A=matrix, fname=file, append={0,1}: A is integer
	!	11) vec2csv(A,fname,append)   A=matrix, fname=file, append={0,1}
	!	12) veci2csv(A,fname,append)   A=matrix, fname=file, append={0,1}: A is integer
	!User-Defined Types (Structures) for value functiIons and policy functions
	!	a) val_struct: VR, VD, VN, VW, VU, V
	!	b) pol_struct: aR, aD,aU,aN,aW,gapp,gwork, gapp_dif,gwork_dif
	!	c) moments_struct
	!	d) hist_struct
	!**********************************************************!


	!------------------------------------------------------------------
	! a) val_struct: VR, VD, VN, VW, VU, V
	!------------------------------------------------------------------
	type val_struct
		! see alloc_valpol  for dimensions' definitions
		real(dp), allocatable:: 	VR(:,:,:,:), &		!Retirement
					VD(:,:,:,:), &		!Disabled
					VN(:,:,:,:,:,:,:), &	!Long-term Unemployed
					VW(:,:,:,:,:,:,:), &	!Working
					VU(:,:,:,:,:,:,:), &	!Unemployed
					V(:,:,:,:,:,:,:), &		!Participant
					V_CF(:,:,:,:,:,:,:)		!Counter-factual w/o DI
		integer :: alloced 	! status indicator
		integer :: inited 	! intiialized indicator

	end type val_struct

	!------------------------------------------------------------------
	!	b) pol_struct: aR, aD,aU,aN,aW,gapp,gwork, gapp_dif,gwork_dif
	!------------------------------------------------------------------
	type pol_struct

		real(dp), allocatable ::	gapp_dif(:,:,:,:,:,:,:), &
					gwork_dif(:,:,:,:,:,:,:) ! latent value of work/apply
		integer, allocatable ::	aR(:,:,:,:), aD(:,:,:,:), aU(:,:,:,:,:,:,:), &
					aN(:,:,:,:,:,:,:), aW(:,:,:,:,:,:,:)
		integer, allocatable ::	gapp(:,:,:,:,:,:,:), &

					gwork(:,:,:,:,:,:,:) !integer choice of apply/work
		integer :: alloced
	end type pol_struct

	type moments_struct
		real(dp) :: work_coefs(Nk), di_coefs(Nk),ts_emp_coefs(nj+1)
		real(dp) :: di_rate(TT-1), work_rate(TT-1), accept_rate(TT-1) !by age
		real(dp) :: di_rateD(nd), work_rateD(nd), work_rateDage(nd,2), accept_rateD(nd)
		integer :: alloced
		real(dp) :: work_cov_coefs(Nk,Nk),di_cov_coefs(Nk,Nk),ts_emp_cov_coefs(nj+1,nj+1)
		real(dp) :: s2,avg_di,init_di
		real(dp) :: init_hlth_acc,init_diprob,avg_diprob,avg_hlth_acc,d1_diawardfrac, d2_diawardfrac, d3_diawardfrac, &
		&			init_diaward,init_diaward_discr, old_initdiawdfrac, &
		&			diaward_ageeffect,old_diawardfrac, end_hlth_acc
		real(dp) :: hlth_acc_rt(TT-1)

	end type moments_struct

	type hist_struct
		real(dp), allocatable :: work_dif_hist(:,:), app_dif_hist(:,:) !choose work or not, apply or not -- latent value
		real(dp), allocatable :: di_prob_hist(:,:) !choose apply or not * prob of getting it-- latent value
		integer, allocatable :: hlth_voc_hist(:,:)
		real(dp), allocatable :: hlthprob_hist(:,:)
		real(dp), allocatable :: wage_hist(:,:) !realized wages
		real(dp), allocatable :: welfare_hist(:,:)
		integer, allocatable :: z_jt_macroint(:) !endogenous realization of shocks given a sequence
		real(dp), allocatable :: z_jt_panel(:,:)
		! a bunch of explanitory variables to be stacked on each other
		integer, allocatable :: status_hist(:,:) !status: W,U,N,D,R (1:5)
		real(dp), allocatable :: a_hist(:,:)
		real(dp), allocatable :: occgrow_jt(:,:)
		real(dp), allocatable :: occshrink_jt(:,:)
		real(dp), allocatable :: occsize_jt(:,:)
		integer :: alloced
	end type hist_struct

	type shocks_struct
		! arrays related to randomly drawn stuff
		integer, allocatable :: j_i(:)
		real(dp), allocatable :: z_jt_select(:), z_jt_innov(:)
		integer, allocatable :: al_int_hist(:,:)
		integer, allocatable :: d_hist(:,:)
		integer, allocatable :: age_hist(:,:), born_hist(:,:), del_i_int(:), fndsep_i_int(:,:,:)
		real(dp), allocatable:: age_draw(:,:),dead_draw(:,:)
		real(dp), allocatable :: al_hist(:,:), del_i_draw(:),fndsep_i_draw(:,:),fndarrive_draw(:,:)
		real(dp), allocatable :: wtr_it(:,:)
		real(dp), allocatable    :: status_it_innov(:,:),health_it_innov(:,:),al_it_innov(:,:)
		real(dp), allocatable    :: jshock_ij(:,:)
		integer,  allocatable    :: drawi_ititer(:,:)
		integer,  allocatable    :: drawt_ititer(:,:)

		integer :: alloced
		integer :: drawn
	end type shocks_struct


	type val_pol_shocks_struct

		! try to do this with pointers
		type(shocks_struct), pointer :: shk_ptr
		type(val_struct)   , pointer :: vfs_ptr
		type(pol_struct)   , pointer :: pfs_ptr

		integer :: pointed
	end type val_pol_shocks_struct

	contains

	!------------------------------------------------------------------------
	!1)UI(e): Unemployment Insurance
	!----------------
	function UI(ein,invol_un)
	!----------------
		!ein: the AIME level
		!invol_un: indicator = 1 if involunarily separated. This aligns with alpha_int=1

		real(dp), intent(in)	:: ein
		integer , intent(in)    :: invol_un
		real(dp) 		:: UI
		!I'm using a replacement rate of UIrr % for now, can be fancier
		if(invol_un .eq. 1) then
			UI = ein*UIrr
		else
			UI = min(LTUrr*ein,LTUrr) ! same income as for VN
		endif

	end function
	!------------------------------------------------------------------------
	!2) DI (can be asset tested)
	!---------------------
	function SSDI(ein,a0in)
	!---------------------


		real(dp), intent(in)	:: ein ! true average income
		integer, intent(in), optional :: a0in !year of accession to DI
		real(dp) 		:: SSDI
		real(dp)		:: e_eval

		if(present(a0in)) then
			if(a0in< 17) then
				e_eval = ein*dble(a0in)/17._dp
			else
				e_eval = ein
			endif
		else
			e_eval = ein
		endif

		!Follows Pistafferi & Low '15
		if (e_eval<DItest1*wmean) then
			SSDI = 0.9_dp*e_eval
		elseif (e_eval<DItest2*wmean) then
			SSDI = 0.9_dp*DItest1*wmean + 0.32_dp*(e_eval-DItest1*wmean)
		elseif (e_eval<DItest3*wmean) then
			SSDI = 0.9_dp*DItest1*wmean + 0.32_dp*(DItest2-DItest1)*wmean+0.15_dp*(e_eval-DItest2*wmean)
		else
			SSDI = 0.9_dp*DItest1*wmean + 0.32_dp*(DItest2-DItest1)*wmean+0.15_dp*(DItest3*wmean-DItest2*wmean)
		endif

	end function
	!------------------------------------------------------------------------
	!3) Social Sec. Retirement
	!--------------------
	function SSret(ein,a0in)
	!--------------------

		real(dp), intent(in)	:: ein
		real(dp), intent(in), optional :: a0in !year of retirement
		real(dp)			:: SSret
		real(dp)			:: e_eval

		if(present(a0in)) then
			if(a0in .le. agegrid(TT-1)) then
				e_eval = ein*0.8_dp
			else
				e_eval = ein
			endif
		else
			e_eval = ein
		endif

		!Follows Pistafferi & Low '15
		if (e_eval<DItest1*wmean) then
			SSret = 0.9*e_eval
		elseif (e_eval<DItest2*wmean) then
			SSret = 0.9*DItest1*wmean + 0.32*(e_eval-DItest1*wmean)
		elseif (e_eval<DItest3*wmean) then
			SSret = 0.9*DItest1*wmean + 0.32*(DItest2-DItest1)*wmean+0.15*(e_eval-DItest2*wmean)
		else
			SSret = 0.9*DItest1*wmean + 0.32*(DItest2-DItest1)*wmean+0.15*(DItest3*wmean-DItest2*wmean)
		endif

	end function
	!------------------------------------------------------------------------
	! 4) Acceptance probability
	!--------------------
	function xifun(idin,trin,itin,wtr_ref,hlthprob, tcal) result(diprob)

		!idin is the health (integer)
		!trin is the wage trend, as realized by the change in returns to skills (real)
		!itin is the age
		!wtr_ref is the reference trend input. Lower than this spurs vocational grid
		!hlthprob is an output, that is the fraction of the acceptance probability associated with health (not vocational)

		real(dp), intent(in):: trin,wtr_ref
		integer, intent(in):: idin,itin
		real(dp), optional :: hlthprob
		integer,  optional :: tcal
		real(dp) :: xifunH,xifunV, hlthfrac,trqtl,trwt,trhr
		integer :: trqtl_L, trqtl_H
		real(dp) :: diprob, prV=1, tcalhr

		!as a multinomial non-nested logit:

		if( present(tcal) .eqv. .false. ) then 
			tcalhr = 0._dp
		else 
			tcalhr = dble(tcal)/dble(Tsim)
		endif

		prV = 1._dp !max(min( (wtr_ref - trin)/wtr_scale, 1._dp), 0._dp)
		!         weakly rising in age (xiagezcoef) and has an intercept (xizcoef)
		if(trin .le. wtr_ref .and. itin .ge. TT-2) then
			!xifunV =  xiagezcoef*(agegrid(itin)-agegrid(TT-3))/(agegrid(TT-1)-agegrid(1)) + xizcoef
			xifunV =  xiagezcoef + xizcoef
		elseif(trin .le. wtr_ref .and. itin .lt. TT-2 ) then
			xifunV =  xizcoef
		else
			xifunV = -10.
		endif
		if(  idin==2) then
			xifunH = xid2coef
		elseif(idin==3) then
			xifunH = xid3coef
		else
			xifunH =xid1coef
		endif

		!apply xifunV trend:
		xifunV = xifunV * xiVtrend* tcalhr
		xifunH = xifunH - xifunV*(xiVtrend-1)* tcalhr

		diprob = dexp(xifunV)/(1._dp + dexp(xifunV)*prV + dexp(xifunH)) + dexp(xifunH)/(1._dp + dexp(xifunV)*prV + dexp(xifunH))
		if( diprob>0._dp) then
			hlthfrac = dexp(xifunH)/(1._dp + dexp(xifunV)*prV + dexp(xifunH))/diprob
		else
			hlthfrac = 1._dp
		endif
		diprob = max(min(diprob,1._dp),0._dp)

		diprob = 1._dp - max(0.,1.-diprob)**(1._dp / proc_time1)

		if((itin == 1) .and. (ineligNoNu .eqv. .false.)) then
			diprob = diprob*eligY
		endif

		if(present(hlthprob) .eqv. .true.) then
			hlthprob = hlthfrac*diprob
			if( isnan(hlthprob)  .or. hlthprob > HUGE(hlthprob) ) hlthprob = 1._dp
		endif

	end function

	!------------------------------------------------------------------------
	! 4) Utility Function
	!--------------------
	function util(cin,din,wkin)
	!--------------------

		real(dp),intent(in)	:: cin
		integer, intent(in)	:: din, wkin
		real(dp)			:: util


		if(gam> 1+1e-5_dp .or. gam < 1-1e-5_dp) then
			util = ((cin*dexp(theta*dble(din-1)+eta*dble(wkin-1)))**(1._dp-gam) )/(1._dp-gam)
		else
			util = dlog(cin*dexp(theta*dble(din-1)+eta*dble(wkin-1)))
		endif
		util = util + util_const

	end function

	!------------------------------------------------------------------------
	! 5) Earnings Index Function
	!--------------------
	function Hearn(tin,ein,wgin)
	!--------------------

		real(dp), intent(in)	:: wgin
		integer, intent(in)	:: ein, tin
		real(dp)			:: Hearn

		Hearn = (wgin + egrid(ein)*(tlen*agegrid(tin)-1._dp) )/ &
!		&	(tlen*max(dble(agegrid(tin)),35._dp) )
		&	(tlen*agegrid(tin) )
		Hearn = max( egrid(ein),Hearn )

	end function

	!------------------------------------------------------------------------
	! 6) Wage Function
	!--------------------
	function wage(levin,aiin,din,tin)
	!--------------------

		!levin is the level/trend that is occupation specific
		!aiin is the alpha idiosyncratic shock
		!din is integer representing the health state
		!tin is the age integer, which determines where in the age profile (wtau)

		real(dp), intent(in)	:: levin, aiin
		integer, intent(in)	:: din, tin
		real(dp)			:: wage

		wage = dexp( levin + aiin+wd(din)+wtau(tin) )

	end function

	!------------------------------------------------------------------------
	! 7) Locate Function
	!--------------------
	function finder(xx,x)
	!--------------------

		real(dp), dimension(:), intent(IN) :: xx
		real(dp), intent(IN) :: x
		integer :: locate,finder
		integer :: nf,il,im,iu

		nf=size(xx)
		il=0
		iu=nf+1
		do
			if (iu-il <= 1) exit !converged

			im=(iu+il)/2

			if (x >= xx(im)) then
				il=im
			else
				iu=im
			end if
		end do
		if (x <= xx(1)+epsilon(xx(1))) then
			locate=1
		else if (x >= xx(nf)-epsilon(xx(nf))) then
			locate=nf-1
		else
			locate=il
		end if
		finder = locate
	end function

	!------------------------------------------------------------------------
	! 8) Write a Matrix to .csv
	!--------------------
	subroutine mat2csv(A,fname,append)
	!--------------------
	! A: name of matrix
	! fname: name of file, should end in ".csv"
	! append: a 1 or 0.  1 if matrx should add to end of existing file, 0, to overwrite.
	implicit none
	real(dp), dimension(:,:), intent(in) :: A
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
			write(1,'(A1)' , advance='no') ","
		enddo
		write(1,FMT) A(ri,c)
	enddo
!	write(1,*) " "! trailing space
	close(1)

	end subroutine mat2csv

	!------------------------------------------------------------------------
	! 9) Write a Matrix to .csv
	!--------------------
	subroutine mati2csv(A,fname,append)
	!--------------------
	implicit none
	integer, dimension(:,:), intent(in) :: A
	character(LEN=*), intent(in) :: fname
	integer, intent(in), optional :: append
	CHARACTER(LEN=*), PARAMETER  :: FMT = "(I8.1)"
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

	do ri=1,r
		do ci = 1,c-1
			write(1,FMT, advance='no') A(ri,ci)
			write(1,'(A1)' , advance='no') ","
		enddo
		write(1,FMT) A(ri,c)
	enddo
!	write(1,*) " "! trailing space
	close(1)

	end subroutine mati2csv


	!------------------------------------------------------------------------
	! 10) Write a Vector to .csv
	!--------------------
	subroutine vec2csv(A,fname,append)
	!--------------------

	real(dp), dimension(:), intent(in) :: A
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
!	write(1,*) " "! trailing space
	close(1)

	end subroutine vec2csv


	!------------------------------------------------------------------------
	! 11) Write a Vector to .csv
	!--------------------
	subroutine veci2csv(A,fname,append)
	!--------------------

		integer, dimension(:), intent(in) :: A
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
!		write(1,*) " "! trailing space
		close(1)

	end subroutine veci2csv

	subroutine csv2mat(fname,A)
		character(len=*), intent(in) :: fname
		real(dp), dimension(:,:) :: A
		integer :: nrow,ri

		open(unit= fread, file = fname)
		do ri=1,nrow
			read(fread,*) A(ri,:)
		enddo
		close(fread)

	end subroutine csv2mat

	subroutine csv2mati(fname,A)
		character(len=*), intent(in) :: fname
		integer, dimension(:,:) :: A
		integer :: nrow,ri

		open(unit= fread, file = fname)
		do ri=1,nrow
			read(fread,*) A(ri,:)
		enddo
		close(fread)

	end subroutine csv2mati


	subroutine csv2vec(fname,A)
		character(len=*), intent(in) :: fname
		real(dp), dimension(:) :: A
		integer :: nrow,ri

		open(unit= fread, file = fname)
		do ri=1,nrow
			read(fread,*) A(ri)
		enddo
		close(fread)
	end subroutine csv2vec

	subroutine csv2veci(fname,A)
		character(len=*), intent(in) :: fname
		integer, dimension(:) :: A
		integer :: nrow,ri

		open(unit= fread, file = fname)
		do ri=1,nrow
			read(fread,*) A(ri)
		enddo
		close(fread)
	end subroutine csv2veci


	!------------------------------------------------------------------------
	! 12) Run an OLS regression
	!------------------------------------------------------------------------
	subroutine OLS(XX,Y,coefs,cov_coef, hatsig2, status)
		real(dp), dimension(:,:) :: XX
		real(dp), dimension(:) :: Y
		real(dp), dimension(:), intent(out) :: coefs
		real(dp), dimension(:,:), intent(out) :: cov_coef
		real(dp), intent(out) :: hatsig2
		integer, intent(out) :: status

		integer :: nX, nY, nK
		real(dp), dimension(size(XX, dim = 2),size(XX, dim = 2)) :: XpX,XpX_fac,XpX_inv
		real(dp), dimension(:,:), allocatable :: XX1
		real(dp), dimension(:), allocatable :: fitted,resids
		integer :: i, regstatus
		real(dp) :: alpha,beta

		external dgemm,dgemv
		external dpotrs,dpotrf,dpotri

		nK = size(XX, dim = 2)
		nX = size(XX, dim = 1)
		nY = size(Y)

		allocate(fitted(nX))
		allocate(resids(nX))
		allocate(XX1(nX,nK))
		coefs = 0._dp
		cov_coef = 0._dp
		alpha = 1._dp
		beta = 0._dp
		!XX1=XX
		if(nY /= nX ) then
			if(verbose>0 ) print *, 'size of X and Y in regression not compatible'
			status = 1
		else
			XpX = 0._dp
			call dgemm('T', 'N', nK, nK, nX, alpha, XX, nX, XX, nX,beta, XpX, nK)
			XpX_fac = XpX
			call dpotrf('U',Nk,XpX_fac,Nk,regstatus)
			if(regstatus .eq. 0) then
				! 3/ Multiply LHS of regression and solve it
				call dgemv('T', nX, nK, alpha, XX, nX, Y, 1, beta, coefs, 1)
				call dpotrs('U',nK,1,XpX_fac,nK,coefs,Nk,regstatus)
			else
				if(verbose >0 ) print *, "cannot factor XX, error ",regstatus
			endif
		endif

		if(regstatus .eq. 0) then
			XpX_inv = XpX_fac
			call dpotri('U',nK,XpX_inv,nK,regstatus)
			fitted = 0._dp
			call dgemv('N', nX, nK, alpha, XX, nX, coefs, 1, beta, fitted, 1)
			resids = Y - fitted
			hatsig2 = 0._dp
			do i=1,nY
				hatsig2 = (resids(i)**2)/dble(nX-nK) + hatsig2
			enddo
			if(regstatus .eq. 0) cov_coef = hatsig2*XpX_inv
		endif
		status = regstatus

		deallocate(fitted,resids,XX1)

	end subroutine OLS

	subroutine stdev_v(std, dat)

		real(dp), intent(in) :: dat(:)
		real(dp), intent(out) :: std
		real(dp) :: mean
		integer  :: i

		mean = sum(dat)/dble( size(dat) )
		std  = 0._dp
		do i=1,size(dat)
			std = (dat(i) - mean)**2 + std
		enddo
		std =	(std/dble( size(dat) ))**0.5_dp

	end subroutine stdev_v


	subroutine alloc_hist(hst)

		type(hist_struct) :: hst

		allocate(hst%wage_hist(Nsim,Tsim), stat=hst%alloced)
		allocate(hst%welfare_hist(Nsim,Tsim), stat=hst%alloced)
		allocate(hst%work_dif_hist(Nsim,Tsim), stat=hst%alloced)
		allocate(hst%app_dif_hist(Nsim,Tsim), stat=hst%alloced)
		allocate(hst%di_prob_hist(Nsim,Tsim), stat=hst%alloced)
		allocate(hst%hlth_voc_hist(Nsim,Tsim), stat=hst%alloced)
		allocate(hst%hlthprob_hist(Nsim,Tsim), stat=hst%alloced)
		allocate(hst%status_hist(Nsim,Tsim), stat=hst%alloced)
		allocate(hst%a_hist(Nsim,Tsim), stat=hst%alloced)
		allocate(hst%z_jt_macroint(Tsim), stat=hst%alloced)
		allocate(hst%z_jt_panel(nj,Tsim), stat=hst%alloced)
		allocate(hst%occgrow_jt(nj,Tsim), stat=hst%alloced)
		allocate(hst%occshrink_jt(nj,Tsim), stat=hst%alloced)
		allocate(hst%occsize_jt(nj,Tsim), stat=hst%alloced)

	end subroutine alloc_hist

	subroutine clean_hist(hst)
!cleans everything except z_jt_panel,z_jt_macroint
		type(hist_struct) :: hst

		hst%wage_hist = 0.
		hst%welfare_hist = 0.
		hst%work_dif_hist = 0.
		hst%app_dif_hist = 0.
		hst%di_prob_hist  = 0.
		hst%hlth_voc_hist  = 0
		hst%hlthprob_hist  = 0.
		hst%status_hist  = 0
		hst%a_hist  = 0.
		hst%occgrow_jt  = 0.
		hst%occshrink_jt  = 0.
		hst%occsize_jt  = 0.

	end subroutine clean_hist

	subroutine alloc_valpol(vfs, pfs)
		type(val_struct)   :: vfs
		type(pol_struct)   :: pfs

		!************************************************************************************************!
		! Allocate matrices
		!************************************************************************************************!
		! (early retirement, disability extent, earn hist, assets)
		allocate(vfs%VR(2,nd,ne,na), stat=vfs%alloced)
		allocate(pfs%aR(2,nd,ne,na), stat=pfs%alloced)

		! (disability extent, earn hist, assets, age)
		allocate(vfs%VD(nd,ne,na,TT), stat=vfs%alloced)
		allocate(pfs%aD(nd,ne,na,TT-1), stat=pfs%alloced)

		! (occupation X ind exposure, ind disb. risk X ind. wage, disab. extent, earn hist, assets, agg shock, age)
		allocate(vfs%VN(nl*ntr,ndi*nal,nd,ne,na,nz,TT), stat=vfs%alloced)
		allocate(vfs%VU(nl*ntr,ndi*nal,nd,ne,na,nz,TT), stat=vfs%alloced)
		allocate(vfs%VW(nl*ntr,ndi*nal,nd,ne,na,nz,TT), stat=vfs%alloced)
		allocate(vfs%V(nl*ntr,ndi*nal,nd,ne,na,nz,TT), stat=vfs%alloced)
		allocate(vfs%V_CF(nl*ntr,ndi*nal,nd,ne,na,nz,TT), stat=vfs%alloced)
		allocate(pfs%aN(nl*ntr,ndi*nal,nd,ne,na,nz,TT-1), stat=pfs%alloced)
		allocate(pfs%aW(nl*ntr,ndi*nal,nd,ne,na,nz,TT-1), stat=pfs%alloced)
		allocate(pfs%aU(nl*ntr,ndi*nal,nd,ne,na,nz,TT-1), stat=pfs%alloced)
		allocate(pfs%gwork(nl*ntr,ndi*nal,nd,ne,na,nz,TT-1), stat=pfs%alloced)
		allocate(pfs%gapp(nl*ntr,ndi*nal,nd,ne,na,nz,TT-1), stat=pfs%alloced)

		allocate(pfs%gapp_dif(nl*ntr,ndi*nal,nd,ne,na,nz,TT), stat=pfs%alloced)
		allocate(pfs%gwork_dif(nl*ntr,ndi*nal,nd,ne,na,nz,TT), stat=pfs%alloced)

	end subroutine

	subroutine alloc_econ(vfs, pfs,hst)

	! Structure to communicate everything
		type(val_struct)   :: vfs
		type(pol_struct)   :: pfs
		type(hist_struct)  :: hst

		call alloc_valpol(vfs,pfs)
		call alloc_hist(hst)

	end subroutine alloc_econ

	subroutine alloc_shocks(shk)

		type(shocks_struct) :: shk

		allocate(shk%age_hist(Nsim,Tsim), stat=shk%alloced)
		allocate(shk%age_draw(Nsim,Tsim+3), stat=shk%alloced)
		allocate(shk%dead_draw(Nsim,Tsim), stat=shk%alloced)
		allocate(shk%al_hist(Nsim,Tsim), stat=shk%alloced)
		allocate(shk%al_int_hist(Nsim,Tsim), stat=shk%alloced)
		allocate(shk%born_hist(Nsim,Tsim), stat=shk%alloced)
		allocate(shk%del_i_int(Nsim), stat=shk%alloced)
		allocate(shk%d_hist(Nsim,Tsim), stat=shk%alloced)
		allocate(shk%del_i_draw(Nsim), stat=shk%alloced)
		allocate(shk%fndsep_i_draw(Nsim,2), stat=shk%alloced)
		allocate(shk%fndsep_i_int(Nsim,2,nz), stat=shk%alloced)
		allocate(shk%fndarrive_draw(Nsim,Tsim), stat=shk%alloced)
		!this must be big enough that we are sure it's big enough that can always find a worker
		allocate(shk%drawi_ititer(Nsim,Nsim/5), stat=shk%alloced)
		allocate(shk%drawt_ititer(Nsim,Nsim/5), stat=shk%alloced)
		allocate(shk%j_i(Nsim), stat=shk%alloced)
		allocate(shk%jshock_ij(Nsim,nj), stat=shk%alloced)
		allocate(shk%status_it_innov(Nsim,Tsim), stat=shk%alloced)
		allocate(shk%health_it_innov(Nsim,Tsim), stat=shk%alloced)
		allocate(shk%al_it_innov(Nsim,Tsim), stat=shk%alloced)
		allocate(shk%z_jt_innov(Tsim))
		allocate(shk%z_jt_select(Tsim))
		allocate(shk%wtr_it(Nsim,Tsim))




	end subroutine alloc_shocks


	subroutine dealloc_hist(hst)

		type(hist_struct) :: hst

		deallocate(hst%wage_hist     , stat=hst%alloced)
		deallocate(hst%welfare_hist  , stat=hst%alloced)
		deallocate(hst%work_dif_hist , stat=hst%alloced)
		deallocate(hst%app_dif_hist , stat=hst%alloced)
		deallocate(hst%di_prob_hist , stat=hst%alloced)
		deallocate(hst%hlth_voc_hist, stat=hst%alloced)
		deallocate(hst%hlthprob_hist, stat=hst%alloced)
		deallocate(hst%status_hist  , stat=hst%alloced)
		deallocate(hst%a_hist , stat=hst%alloced)
		deallocate(hst%z_jt_macroint, stat=hst%alloced)
		deallocate(hst%z_jt_panel, stat=hst%alloced)
		deallocate(hst%occgrow_jt, stat=hst%alloced)
		deallocate(hst%occshrink_jt, stat=hst%alloced)
		deallocate(hst%occsize_jt, stat=hst%alloced)
		hst%alloced = 0

	end subroutine dealloc_hist

	subroutine dealloc_valpol(vfs,pfs)
		! Structure to communicate everything
			type(val_struct) :: vfs
			type(pol_struct) :: pfs
			deallocate(vfs%VR, stat=vfs%alloced)
			deallocate(pfs%aR, stat=pfs%alloced)

			! (disability extent, earn hist, assets, age)
			deallocate(vfs%VD, stat=vfs%alloced)
			deallocate(pfs%aD, stat=pfs%alloced)

			! (occupation X ind exposure, ind disb. risk X ind. wage, disab. extent, earn hist, assets, agg shock, age)
			deallocate(vfs%VN  , stat=vfs%alloced)
			deallocate(vfs%VU  , stat=vfs%alloced)
			deallocate(vfs%VW  , stat=vfs%alloced)
			deallocate(vfs%V   , stat=vfs%alloced)
			deallocate(vfs%V_CF, stat=vfs%alloced)
			deallocate(pfs%aN, stat=pfs%alloced)
			deallocate(pfs%aW, stat=pfs%alloced)
			deallocate(pfs%aU, stat=pfs%alloced)
			deallocate(pfs%gwork, stat=pfs%alloced)
			deallocate(pfs%gapp, stat=pfs%alloced)

			deallocate(pfs%gapp_dif , stat=pfs%alloced)
			deallocate(pfs%gwork_dif , stat=pfs%alloced)

			vfs%alloced = 0
			pfs%alloced = 0

	end subroutine dealloc_valpol

	subroutine dealloc_econ(vfs,pfs,hst)

	! Structure to communicate everything
		type(val_struct) :: vfs
		type(pol_struct) :: pfs
		type(hist_struct):: hst

		call dealloc_valpol(vfs,pfs)
		call dealloc_hist(hst)


	end subroutine dealloc_econ


	subroutine dealloc_shocks(shk)

		type(shocks_struct) :: shk

		deallocate(shk%age_hist , stat=shk%alloced)
		deallocate(shk%age_draw , stat=shk%alloced)
		deallocate(shk%dead_draw , stat=shk%alloced)
		deallocate(shk%born_hist , stat=shk%alloced)
		deallocate(shk%al_hist , stat=shk%alloced)
		deallocate(shk%al_int_hist , stat=shk%alloced)
		deallocate(shk%j_i , stat=shk%alloced)
		deallocate(shk%jshock_ij, stat=shk%alloced)
		deallocate(shk%status_it_innov , stat=shk%alloced)
		deallocate(shk%health_it_innov , stat=shk%alloced)
		deallocate(shk%al_it_innov, stat=shk%alloced)
		deallocate(shk%d_hist , stat=shk%alloced)
		deallocate(shk%del_i_int, stat=shk%alloced)
		deallocate(shk%del_i_draw, stat=shk%alloced)
		deallocate(shk%fndsep_i_int, stat=shk%alloced)
		deallocate(shk%fndsep_i_draw, stat=shk%alloced)
		deallocate(shk%fndarrive_draw, stat=shk%alloced)
		deallocate(shk%z_jt_innov)
		deallocate(shk%z_jt_select)
		deallocate(shk%wtr_it)
		!this must be big enough that we are sure it's big enough that can always find a worker
		deallocate(shk%drawi_ititer)
		deallocate(shk%drawt_ititer)


	end subroutine dealloc_shocks

end module helper_funs


!**************************************************************************************************************!
!**************************************************************************************************************!
!**************************************************************************************************************!
!**************************************************************************************************************!
! Compute data from the model
!**************************************************************************************************************!
!**************************************************************************************************************!
!**************************************************************************************************************!
!**************************************************************************************************************!


module model_data

	use V0para
	use helper_funs
! should put comp_ustats here too

	implicit none

	contains


	subroutine comp_ustats(hst,shk,urt,Efrt,Esrt)

		! Takes in structures hst and shk, both of which should already have simulation results in them
		! Outputs urt,Efrt,Esrt

		type(shocks_struct) :: shk
		type(hist_struct):: hst
		type(val_struct) :: vfs
		real(dp), intent(out) :: urt,Efrt,Esrt
		real(dp) :: Nunemp,Nlf, Nsep,Nfnd,ltu,Nun,work_dif_hr
		integer :: i, j, it,duri

		ltu  = 0._dp
		urt  = 0._dp
		Nlf  = 0._dp
		Nunemp = 0._dp
		Nsep = 0._dp
		Nfnd = 0._dp
		Nun  = 0._dp
		do i = 1,Nsim
			duri = 0
			do it=1,Tsim
				work_dif_hr = hst%work_dif_hist(i,it)
				if((hst%status_hist(i,it)<=2) .and. &
				&	(hst%status_hist(i,it)>0) .and. (shk%age_hist(i,it)>0)) then
					Nlf = Nlf+1._dp
					if(hst%status_hist(i,it) == 2 .and. (work_dif_hr > 0._dp )) then
						Nunemp = Nunemp + 1._dp
						if(duri == 0 .and. it>1) & !count a new spell
							& Nsep = Nsep+1._dp
						duri = duri+1
						if(duri < 6) ltu = ltu + 1._dp !counted as the fraction of short-term unemployed
					else if( hst%status_hist(i,it) == 2 .and. work_dif_hr .le. 0._dp ) then
						duri = duri+1
					else
						if(duri > 0 ) & !just found (already conditioned on status ==1 or status==2)
							& Nfnd = Nfnd+1._dp
						duri = 0
					endif
				endif
				if(hst%status_hist(i,it)>2 .or. hst%status_hist(i,it)==1) duri = 0
				if(hst%status_hist(i,it)==3) then
					if(it >1) then
						if (hst%status_hist(i,it-1)==2 .and. hst%work_dif_hist(i,it)>0._dp ) Nun = Nun+1._dp
					endif
				endif
			enddo
		enddo
		urt  = Nunemp/Nlf
		Esrt = Nsep/(Nlf-Nunemp)
		if(Nunemp > 0._dp .and. ((Nunemp - Nun) > 0) ) then
			ltu  = 1._dp - ltu/Nunemp
			Efrt = (Nfnd- Nun)/(Nunemp - Nun)
		else
			ltu  = 0._dp
			Efrt = 0._dp
		endif
	end subroutine comp_ustats




	subroutine moments_compute(hst,moments_sim,shk)

		type(moments_struct) 	:: moments_sim
		type(hist_struct)	:: hst
		type(shocks_struct) :: shk

		integer :: i, ij,id,it,ial,st,si,age_hr,status_hr
		integer :: totage(TT),totD(TT),totW(TT),totst(Tsim),total(nal), tot3al(nal), &
				& tot3age(TT-1),totage_st(TT,Tsim),tot_applied(TT-1)

		real(dp) :: dD_age(TT), dD_t(Tsim),a_age(TT),a_t(Tsim),alworkdif(nal),alappdif(nal), &
				& workdif_age(TT-1),  alD(nal), alD_age(nal,TT-1), &
				& status_Nt(5,Tsim),DIatriskpop,napp_t,ninsur_app, dicont_hr=0., &
				& init_diaward_discr, init_diprob_disc, &
				& tot_apppr, vocprob_age(TT),initdiage(TT)
		real(dp) :: end_diprob
		real(dp) :: di_rateD_denom(nd) ,work_rateD_denom(nd),work_rateDage_denom(nd,2),accept_rateD_denom(nd)
		real(dp) :: d1_diawardfrac_disc,d2_diawardfrac_disc,d3_diawardfrac_disc,diawrd_denom, &
		&		old_diawardfrac_disc, avgdiprob, initdenialprob, initdenialprob_denom, old_initdiawdfrac_disc
		integer  :: apl_i, allow_i, apl_periods

		if(hst%alloced /= 0) then
			if(verbose >= 1) print *, "not correctly passing hists_struct to moments"
		endif

		!initialize all the sums
		totage	= 0
		totst	= 0
		totD	= 0
		totW	= 0
		total	= 0
		tot3al	= 0
		tot3age	= 0

		totage_st=0
		tot_applied = 0
		tot_apppr = 0._dp

		dD_age 		= 0._dp
		dD_t 		= 0._dp
		a_age 		= 0._dp
		a_t 		= 0._dp
		alD		= 0._dp
		alD_age		= 0._dp

		workdif_age 	= 0._dp
		vocprob_age = 0._dp
		initdiage 	= 0._dp
		alworkdif	= 0._dp
		alappdif	= 0._dp
		status_Nt 	= 0._dp
		d1_diawardfrac_disc = 0._dp
		d2_diawardfrac_disc = 0._dp
		d3_diawardfrac_disc = 0._dp
		old_diawardfrac_disc = 0._dp
		diawrd_denom = 0._dp
		moments_sim%hlth_acc_rt = 0._dp
		moments_sim%avg_hlth_acc = 0._dp
		moments_sim%d1_diawardfrac = 0._dp
		moments_sim%d2_diawardfrac = 0._dp
		moments_sim%d3_diawardfrac = 0._dp
		moments_sim%old_diawardfrac = 0._dp
		moments_sim%di_rateD = 0._dp
		moments_sim%work_rateD = 0._dp
		moments_sim%work_rateDage = 0._dp
		moments_sim%accept_rateD = 0._dp
		moments_sim%avg_diprob = 0._dp
		di_rateD_denom = 0._dp
		work_rateD_denom = 0._dp
		work_rateDage_denom =0._dp
		accept_rateD_denom = 0._dp
		avgdiprob = 0._dp

		id = 0
		do si=1,Nsim
			do st = (TossYears*itlen),Tsim
				if( (shk%age_hist(si,st) >0) .and. &
				&	 hst%status_hist(si,st) == 3 .and. hst%app_dif_hist(si,st) >0  ) then
					avgdiprob = avgdiprob + hst%di_prob_hist(si,st)
					id = id + 1
				endif
			enddo
		enddo
		avgdiprob = avgdiprob/dble(id)

		moments_sim%avg_diprob = 1._dp - (1._dp - avgdiprob)**tlen
		initdenialprob_denom = 0._dp
		initdenialprob = 0._dp

		do si = 1,Nsim
			apl_i =0
			allow_i =0
			apl_periods =0
			do st = (TossYears*itlen),Tsim
				if((hst%status_hist(si,st)>0) .and. (shk%age_hist(si,st) >0) ) then
					age_hr = shk%age_hist(si,st)
					totage_st(age_hr,st) = totage_st(age_hr,st) + 1
					! savings and disability by time
					a_t(st) = hst%a_hist(si,st) + a_t(st)
					status_hr = hst%status_hist(si,st)
					if(status_hr == 4) dD_t(st) = dD_t(st)+shk%d_hist(si,st)

					status_Nt(status_hr,st) = 1._dp + status_Nt(status_hr,st)

					id=shk%d_hist(si,st)
					if( status_hr<=4 .and. status_hr>0 )then
						it = 1
						if(shk%age_hist(si,st) > 1) it = 2
						work_rateD_denom(id) = work_rateD_denom(id)+1._dp
						work_rateDage_denom(id,it) = work_rateDage_denom(id,it)+1._dp
						di_rateD_denom(id) = di_rateD_denom(id)+1._dp
						if(status_hr ==1) then
							moments_sim%work_rateD(id) = moments_sim%work_rateD(id)+1._dp
							moments_sim%work_rateDage(id,it) = moments_sim%work_rateDage(id,it)+1._dp
						endif
						if(status_hr ==4) moments_sim%di_rateD(id) = moments_sim%di_rateD(id)+1._dp
					endif
					if( status_hr==3 .and. hst%app_dif_hist(si,st) >=0.) then
						accept_rateD_denom(id) =accept_rateD_denom(id) + 1._dp
						moments_sim%accept_rateD(id) = moments_sim%accept_rateD(id) +hst%di_prob_hist(si,st)
					endif

					! disability by age and age X shock
					if(hst%status_hist(si,st) == 1) totW(age_hr) = totW(age_hr) + 1
					if(hst%status_hist(si,st) < 3) &
					&	workdif_age(age_hr) = workdif_age(age_hr) + hst%work_dif_hist(si,st)

					if(hst%hlth_voc_hist(si,st) >0) then
						moments_sim%hlth_acc_rt(age_hr) = 2.- dble(hst%hlth_voc_hist(si,st)) + moments_sim%hlth_acc_rt(age_hr)
						tot_applied(age_hr) = tot_applied(age_hr)+1
						moments_sim%avg_hlth_acc = 2.- dble(hst%hlth_voc_hist(si,st)) + moments_sim%avg_hlth_acc
					endif
					if(hst%status_hist(si,st) == 4) then
						totD(age_hr) = totD(age_hr) + 1
						dD_age(age_hr) = dD_age(age_hr)+shk%d_hist(si,st)
						! associate this with its shock
						do ial = 1,nal
							if(  (shk%al_hist(si,st) <= alfgrid(ial,shk%d_hist(si,st))) &
							&	.and. (shk%al_hist(si,st) >= alfgrid(ial,shk%d_hist(si,st) )) &
							&	.and. (hst%status_hist(si,st) == 4 )) &
							&	alD_age(ial,age_hr) = 1._dp + alD_age(ial,age_hr)
						enddo
					elseif(hst%status_hist(si,st) == 3) then
						if( hst%app_dif_hist(si,st) > 0. ) then
							tot3age(age_hr) = tot3age(age_hr) + 1
							vocprob_age(age_hr) = hst%di_prob_hist(si,st)-hst%hlthprob_hist(si,st) &
							&	+ vocprob_age(age_hr)
						endif
					endif

					! assets by age
					a_age(age_hr) = hst%a_hist(si,st) +a_age(age_hr)

					! Fraction of awards to good health/old
					if(st>1) then
						if( hst%status_hist(si,st) == 4  .and. hst%status_hist(si,st-1) ==3) then
							diawrd_denom = diawrd_denom + 1._dp
							if( shk%d_hist(si,st)==1 ) &
							& 	d1_diawardfrac_disc = d1_diawardfrac_disc+1._dp
							if( shk%d_hist(si,st)==2 ) &
							& 	d2_diawardfrac_disc = d2_diawardfrac_disc+1._dp
							if( shk%d_hist(si,st)==3 ) &
							& 	d3_diawardfrac_disc = d3_diawardfrac_disc+1._dp

							if( shk%age_hist(si,st)>=TT-2 .and. shk%age_hist(si,st)<=TT-1 ) &
							& 	old_diawardfrac_disc = old_diawardfrac_disc + 1._dp
						endif
					endif

					if( hst%status_hist(si,st) == 3 .and. hst%app_dif_hist(si,st) >0 ) then
						tot_apppr = tot_apppr+hst%di_prob_hist(si,st)
						if( shk%d_hist(si,st)==1) then
							moments_sim%d1_diawardfrac = hst%di_prob_hist(si,st) +moments_sim%d1_diawardfrac
						elseif(shk%d_hist(si,st)==2) then
							moments_sim%d2_diawardfrac = hst%di_prob_hist(si,st) +moments_sim%d2_diawardfrac
						elseif(shk%d_hist(si,st)==3) then
							moments_sim%d3_diawardfrac = hst%di_prob_hist(si,st) +moments_sim%d3_diawardfrac
						endif
						if( shk%age_hist(si,st)>=TT-2 .and. shk%age_hist(si,st)<=TT-1 ) then
							moments_sim%old_diawardfrac  = dexp(hst%di_prob_hist(si,st)-avgdiprob)/(1._dp + dexp(hst%di_prob_hist(si,st)-avgdiprob)) &
							&	+ moments_sim%old_diawardfrac
						endif
					endif

					! disability and app choice by shock level
					id = shk%d_hist(si,st)
					do ial = 1,nal
						if( shk%al_hist(si,st) <= alfgrid(ial,id)+2*epsilon(1._dp) &
						&	.and. shk%al_hist(si,st) >= alfgrid(ial,id)-2*epsilon(1._dp)) then
							if(hst%status_hist(si,st) <3) then
								!work choice:
								alworkdif(ial) = alworkdif(ial) + hst%work_dif_hist(si,st)
								total(ial) = total(ial) + 1
							endif
							if(hst%status_hist(si,st) == 3) then
								alappdif(ial) = alappdif(ial) + hst%app_dif_hist(si,st)
								tot3al(ial) = tot3al(ial) + 1
							endif
							if(hst%status_hist(si,st) == 4) &
								& alD(ial) = 1._dp + alD(ial)
						endif
					enddo
					!denial?
					if(hst%status_hist(si,st)==3 .and. hst%app_dif_hist(si,st) >0 .and. st < (init0_yrs+init_yrs)*itlen ) then
						apl_i = 1
						apl_periods = apl_periods + 1
					endif
					if( hst%status_hist(si,st)==4 .and. apl_i ==1 ) then
						apl_i =0
						apl_periods =0
						initdenialprob_denom = initdenialprob_denom + 1._dp
					endif
				endif !status(si,st) >0
			enddo!st = 1,Tsim
			if(apl_i==1) then
				initdenialprob = initdenialprob + dble(apl_periods)/proc_time1
			endif
		enddo ! si=1,Nsim

		initdenialprob = initdenialprob/initdenialprob_denom

		moments_sim%di_rateD   = moments_sim%di_rateD/di_rateD_denom
		moments_sim%work_rateD   = moments_sim%work_rateD/work_rateD_denom
		moments_sim%accept_rateD = moments_sim%accept_rateD/accept_rateD_denom
		moments_sim%work_rateDage(:,1)   = moments_sim%work_rateDage(:,1)/work_rateDage_denom(:,1)
		moments_sim%work_rateDage(:,2)   = moments_sim%work_rateDage(:,2)/work_rateDage_denom(:,2)
		!overall disability rate (the money stat)
		DIatriskpop =0._dp
		moments_sim%avg_di = 0._dp
		do st=1,Tsim
			DIatriskpop = sum(status_Nt(1:3,st)) + DIatriskpop
			moments_sim%avg_di = status_Nt(4,st) + moments_sim%avg_di
		enddo
		moments_sim%avg_di = moments_sim%avg_di/DIatriskpop
		if(sum(tot_applied) >0) then
			moments_sim%avg_hlth_acc = moments_sim%avg_hlth_acc/dble(sum(tot_applied))
		else
			moments_sim%avg_hlth_acc = 0._dp
		endif
		if( diawrd_denom >0 ) then
			d1_diawardfrac_disc = d1_diawardfrac_disc/diawrd_denom
			d2_diawardfrac_disc = d2_diawardfrac_disc/diawrd_denom
			d3_diawardfrac_disc = d3_diawardfrac_disc/diawrd_denom
		endif
		if(tot_apppr>0) then
			moments_sim%d1_diawardfrac = smth_diaward*moments_sim%d1_diawardfrac/tot_apppr &
			&	+ (1._dp - smth_diaward) * d1_diawardfrac_disc
			moments_sim%d2_diawardfrac = smth_diaward*moments_sim%d2_diawardfrac/tot_apppr &
			&	+ (1._dp - smth_diaward) * d2_diawardfrac_disc
			moments_sim%d3_diawardfrac = smth_diaward*moments_sim%d3_diawardfrac/tot_apppr &
			&	+ (1._dp - smth_diaward) * d3_diawardfrac_disc
			moments_sim%old_diawardfrac = old_diawardfrac_disc/diawrd_denom !smth_diaward*moments_sim%old_diawardfrac /tot_apppr &
			!& 	+ (1._dp - smth_diaward)*old_diawardfrac_disc/diawrd_denom

		else
			moments_sim%d1_diawardfrac = 0._dp
			moments_sim%d2_diawardfrac = 0._dp
			moments_sim%d3_diawardfrac = 0._dp
			moments_sim%old_diawardfrac = 0._dp
		endif


		forall(it=1:TT) totage(it) = sum(totage_st(it,:))
		p1d1_target = p1d1_2545target*totage(1)/sum(totage) + p1d1_4665target*(1._dp - totage(1)/sum(totage))
		p1d2_target = p1d2_2545target*totage(1)/sum(totage) + p1d2_4665target*(1._dp - totage(1)/sum(totage))
		p1d3_target = p1d3_2545target*totage(1)/sum(totage) + p1d3_4665target*(1._dp - totage(1)/sum(totage))


		forall(it=1:TT-1) workdif_age(it) = workdif_age(it)/dble(totage(it)-totD(it))
		do ial=1,nal
			if(tot3al(ial) >0 ) then
				alappdif(ial) = alappdif(ial)/dble(tot3al(ial))
			else
				alappdif(ial) =0._dp
			endif
		enddo
		do ial=1,nal
			if(total(ial)>0) then
				alworkdif(ial) = alworkdif(ial)/dble(total(ial))
				!disability distribution by shock and shock,age
				alworkdif(ial) = alworkdif(ial)/dble(total(ial))
				alD(ial) = dble(alD(ial))/dble(total(ial) + alD(ial))
				do it=1,TT-1
					if(totage(it)>0) alD_age(ial,it) = alD_age(ial,it)/dble(total(ial))/dble(totage(it))
				enddo
			endif
		enddo
		! asset distribution by age, time and disability status
		do it=1,TT
			if( it<TT .and. totage(it)>0) then
				moments_sim%di_rate(it) = dble(totD(it))/dble(totage(it))
				moments_sim%work_rate(it) = dble(totW(it))/dble(totage(it))
			endif
			if(totage(it)>0) a_age(it) = a_age(it)/dble(totage(it))
		enddo

		! status distribution by age
		do st=1,Tsim
				if( sum(status_Nt(:,st))>0._dp ) status_Nt(:,st)= status_Nt(:,st)/sum(status_Nt(:,st))
				if( sum(totage_st(:,st)) >0 ) a_t(st) = a_t(st)/dble(sum(totage_st(:,st)))
		enddo

		do it=1,TT-1
			if(tot_applied(it)>0) moments_sim%hlth_acc_rt(it) = moments_sim%hlth_acc_rt(it)/dble(tot_applied(it))
		enddo

		napp_t = 0._dp
		ninsur_app = 0._dp
		moments_sim%init_di = 0._dp
		init_diaward_discr = 0._dp
		old_initdiawdfrac_disc = 0._dp
		moments_sim%init_diaward = 0._dp
		moments_sim%init_hlth_acc = 0._dp
		moments_sim%end_hlth_acc = 0._dp
		end_diprob = 0._dp
		moments_sim%init_diprob = 0._dp
		moments_sim%old_initdiawdfrac = 0._dp

		do i=1,Nsim
			dicont_hr = 0._dp
			do it=((init0_yrs)*itlen+1),((init_yrs+init0_yrs)*itlen)
				if( hst%status_hist(i,it)<5 .and. hst%status_hist(i,it)>0 .and. shk%age_hist(i,it)>0) then
					ninsur_app = 1._dp + ninsur_app
					if(hst%status_hist(i,it) == 3) then
						if(hst%app_dif_hist(i,it)>(-100.) .and. hst%app_dif_hist(i,it)<100.) then
							!latent value of application smoothed by logit
							dicont_hr = dexp(smth_dicont*hst%app_dif_hist(i,it))/(1._dp+dexp(smth_dicont*hst%app_dif_hist(i,it)))
						endif

						if( hst%app_dif_hist(i,it) >=0. .and. hst%di_prob_hist(i,it)>0.) then
							napp_t = napp_t+1._dp
							moments_sim%init_hlth_acc       = hst%hlthprob_hist(i,it) + moments_sim%init_hlth_acc
							moments_sim%init_diprob         = hst%di_prob_hist(i,it)  + moments_sim%init_diprob
							if(shk%age_hist(i,it) .ge. TT-2) moments_sim%old_initdiawdfrac   = hst%di_prob_hist(i,it)  + moments_sim%old_initdiawdfrac
						endif
					endif

					!if get DI then add the latent value when applied
					if(hst%status_hist(i,it) == 4 ) then
						if( it>1 ) then
							if(hst%status_hist(i,it-1) == 3) then
								 init_diaward_discr = init_diaward_discr+1._dp
								 initdiage(shk%age_hist(i,it)) = 1._dp + initdiage(shk%age_hist(i,it))
								 if(shk%age_hist(i,it) .ge. TT-2 .and. shk%age_hist(i,it) .le. TT-1) &
								 &	old_initdiawdfrac_disc = old_initdiawdfrac_disc+1._dp
							endif
							moments_sim%init_di= moments_sim%init_di+1._dp
						endif
						if(hst%status_hist(i,it) == 4 .and. (it==1 .or. dicont_hr==0._dp)) then
							moments_sim%init_di= moments_sim%init_di+1._dp
						endif
					elseif( hst%status_hist(i,it) == 3 ) then
						moments_sim%init_di= moments_sim%init_di+dicont_hr*hst%di_prob_hist(i,it)
						if(hst%app_dif_hist(i,it) >=0 ) moments_sim%init_diaward = hst%di_prob_hist(i,it) + moments_sim%init_diaward
					endif
				endif
			enddo
			do it=(Tsim - init_yrs*itlen+1),Tsim !symetric: if we take 3 yrs for initial stats, take 3 yrs for end stats
				if( hst%status_hist(i,it) == 3 .and. shk%age_hist(i,it)>0) then
						if( hst%app_dif_hist(i,it) >=0. .and. hst%di_prob_hist(i,it)>0.) then
							moments_sim%end_hlth_acc  = hst%hlthprob_hist(i,it) + moments_sim%end_hlth_acc
							end_diprob = hst%di_prob_hist(i,it)  + end_diprob
							
						endif
				endif
			enddo
		enddo
		if(napp_t >0 ) then !napp are the number that applied
			init_diprob_disc = init_diaward_discr / napp_t ! number that got it / number applied
		else
			init_diprob_disc = 0.
		endif
		if(ninsur_app>0._dp) then
			old_initdiawdfrac_disc = old_initdiawdfrac_disc/init_diaward_discr
			moments_sim%init_di= moments_sim%init_di/ninsur_app
			init_diaward_discr = init_diaward_discr/ninsur_app
			moments_sim%init_diaward = moments_sim%init_diaward/ninsur_app
			!make the award rates annual:
			moments_sim%init_diaward = 1._dp-(1._dp-moments_sim%init_diaward)**(tlen)
			init_diaward_discr = 1._dp-(1._dp-init_diaward_discr)**(tlen)
			do i=1,TT
				initdiage(i) = initdiage(i)/sum(initdiage)
			enddo
			if(verbose >2) &
				print *, "init_di, actual and 'smooth' . Ninsured applicants: ", init_diaward_discr,moments_sim%init_diaward, ninsur_app
			! mix the smoothed and discrete:
			moments_sim%init_diaward = smth_diaward*moments_sim%init_diaward + (1._dp -smth_diaward)* init_diaward_discr
		else
			moments_sim%init_diaward = 0._dp
			moments_sim%init_di= 0._dp
			init_diaward_discr = 0._dp
			old_initdiawdfrac_disc = 0._dp
		endif
		moments_sim%init_diaward_discr = init_diaward_discr

		! If we are targeting the age effect in vocational:
		! I am taking the difference in conditional probability of getting on voc, but converting into an annual rate.
		moments_sim%diaward_ageeffect = sum( vocprob_age(TT-2:TT-1)/tot3age(TT-2:TT-1) ) / &
			sum( vocprob_age(1:2)/tot3age(1:2) )
		if(sum( vocprob_age(1:2)/tot3age(1:2) ) <= 0._dp) moments_sim%diaward_ageeffect = 0._dp

		if(napp_t > 0._dp .and. moments_sim%init_diprob>0._dp) then
			!make hlth acc a conditional probability instead of the cumulative probability
			if(verbose>1) then
				print *, moments_sim%init_hlth_acc
				print *, moments_sim%init_diprob
			endif
			moments_sim%init_hlth_acc= moments_sim%init_hlth_acc/moments_sim%init_diprob
			moments_sim%old_initdiawdfrac = moments_sim%old_initdiawdfrac/moments_sim%init_diprob

			!this line needs to be last because init_diprob is the denominator above
			moments_sim%init_diprob  = moments_sim%init_diprob/napp_t
			moments_sim%init_diprob  = smth_diaward* moments_sim%init_diprob  + (1.-smth_diaward)*init_diprob_disc

			moments_sim%init_diprob  = 1._dp - (1._dp - moments_sim%init_diprob)**tlen
		else
			moments_sim%init_hlth_acc= 0._dp
			moments_sim%init_diprob = 0._dp
		endif
		
		if(end_diprob>0._dp) then
			moments_sim%end_hlth_acc = moments_sim%end_hlth_acc/end_diprob
		else 
			moments_sim%end_hlth_acc = 0._dp
		endif
	!	moments_sim%init_hlth_acc= 1._dp - (1._dp - moments_sim%init_hlth_acc)**(moments_sim%init_hlth_acc*proc_time1 + (1._dp - moments_sim%init_hlth_acc)*proc_time2)

		moments_sim%old_initdiawdfrac = smth_diaward*moments_sim%old_initdiawdfrac &
		&	 + (1._dp - smth_diaward)*old_initdiawdfrac_disc
		if(verbose .ge. 2) then
			print *, "Initial denial prob:  ", initdenialprob, " allowance rate: ", moments_sim%init_diprob
		endif
		if(print_lev >= 2) then
			call veci2csv(totst,"pop_st.csv")
			call vec2csv(a_age,"a_age"//trim(caselabel)//".csv")
			call vec2csv(a_t,"a_t"//trim(caselabel)//".csv")
			call vec2csv(moments_sim%di_rate,"di_age"//trim(caselabel)//".csv")

			call vec2csv(moments_sim%work_rate,"work_age"//trim(caselabel)//".csv")
			call vec2csv(workdif_age, "workdif_age"//trim(caselabel)//".csv")
			if(caselabel == "") then
				call vec2csv(alD,"alD"//trim(caselabel)//".csv")
				call mat2csv(alD_age,"alD_age"//trim(caselabel)//".csv")
			endif
			call mat2csv(status_Nt,"status_Nt"//trim(caselabel)//".csv")
			call vec2csv(moments_sim%hlth_acc_rt,"hlth_acc_rt"//trim(caselabel)//".csv")
		endif

	end subroutine moments_compute


end module model_data

!**************************************************************************************************************!
!**************************************************************************************************************!
!**************************************************************************************************************!
!**************************************************************************************************************!
!						Solve the model						       !
!**************************************************************************************************************!
!**************************************************************************************************************!
!**************************************************************************************************************!
!**************************************************************************************************************!


module sol_val

	use V0para
	use helper_funs
! module with subroutines to solve the model and simulate data from it

	implicit none
	real(dp) :: smthV0param
	integer ::  Vevals

	contains

	subroutine maxVR(agei0, id,ie,ia, VR0, iaa0, iaaA, apol,Vout)
		integer, intent(in) :: id,ie,ia
		integer, intent(in) :: iaa0, iaaA
		integer, intent(in) :: agei0 ! the index number of age when first retired
		integer, intent(out) :: apol
		real(dp), intent(out) :: Vout
		real(dp), intent(in) :: VR0(:,:,:,:)
		real(dp) :: Vtest1,Vtest2,chere, Vc1, age0
		integer :: iaa, iw,idd

		age0 = agegrid(agei0+TT-2)
		iw=1
		Vtest1 = -1e6
		apol = iaa0
		do iaa=iaa0,iaaA
			Vevals = Vevals +1
			chere = SSret(egrid(ie),age0)+ R*agrid(ia) - agrid(iaa)
			if( chere .gt. 0.) then !ensure positive consumption
				Vc1 = 0._dp
				do idd =1,nd
					Vc1 = beta*pid(id,idd,1,TT-1)*ptau(TT)*VR0(agei0,id,ie,iaa) + Vc1
				enddo
				Vtest2 = util(chere ,id,iw) + Vc1

				if(Vtest2 > Vtest1  .or. iaa .eq. iaa0 ) then !always replace on the first loop
					Vtest1 = Vtest2
					apol = iaa
				endif
			else   !saved too much and negative consumtion
				exit
			endif
		enddo
		Vout = Vtest1
	end subroutine maxVR

	subroutine maxVD(id,ie,ia,it, VD0, iaa0, iaaA, apol,Vout)
		integer, intent(in) :: id,ie,ia,it
		integer, intent(in) :: iaa0, iaaA
		integer, intent(out) :: apol
		real(dp), intent(out) :: Vout
		real(dp), intent(in) :: VD0(:,:,:,:)
		real(dp) :: pid_here(nd,nd)
		real(dp) :: Vc1,chere,Vtest1,Vtest2,ein
		integer :: iw, iaa,idd

		pid_here = pid(:,:,2,it)

		iw = 1 ! don't work
		Vc1 = beta*((1-ptau(it))*VD0(id,ie,iaa0,it+1)+ptau(it)*VD0(id,ie,iaa0,it))
		ein = egrid(ie)
		chere = SSDI(ein)+R*agrid(ia)-agrid(iaa0)

		Vtest1 = -1e6
		apol = iaa0
		!Find Policy
		do iaa=iaa0,iaaA
			chere = SSDI(egrid(ie))+R*agrid(ia)-agrid(iaa)
			if(chere >0.) then
				Vc1 = 0._dp
				do idd=1,nd
					Vc1 = beta*pid_here(id,idd)*((1-ptau(it))*VD0(idd,ie,iaa,it+1) + ptau(it)*VD0(idd,ie,iaa,it)) + Vc1
				enddo
				Vtest2 = util(chere,id,iw)+ Vc1
				if(Vtest2>Vtest1) then
					Vtest1 = Vtest2
					apol = iaa
				endif
			else
				exit
			endif
		enddo	!iaa
		Vout = Vtest1
		if((Vout <-1e5).and. (verbose >0)) then
			print *,  "SSDI here ",SSDI(egrid(ie))
		endif

	end subroutine maxVD

	subroutine maxVU(il,itr,idi,ial,id,ie,ia,iz,it, VU0,VN0,V0,iaa0,iaaA,apol,Vout)
		integer, intent(in) :: il,itr,idi,ial,id,ie,ia,iz,it
		integer, intent(in) :: iaa0, iaaA
		integer, intent(out) :: apol
		real(dp), intent(out) :: Vout
		real(dp), intent(in) :: VN0(:,:,:,:,:,:,:),V0(:,:,:,:,:,:,:),VU0(:,:,:,:,:,:,:)
		real(dp) :: Vc1,chere,Vtest1,Vtest2, pid_here(nd,nd)
		integer :: iw, iaa,ialal,izz,idd

		pid_here = pid(:,:,idi,it)

		iw = 1 ! don't work
		Vtest1 = -1e6 ! just a very bad number, does not really matter
		apol = iaa0
		do iaa=iaa0,iaaA
			chere = UI(egrid(ie),ial) + R*agrid(ia)-agrid(iaa)
			if(chere>0.) then
				Vtest2 = 0. !Continuation value if don't go on disability
				do izz = 1,nz	 !Loop over z'
				do ialal = ialL,nal !Loop over alpha_i'
				do idd = 1,nd
					if(ial > ialUn) then !unemp by choice
						Vc1 = (1._dp-ptau(it))*(pphi*VN0((il-1)*ntr+itr,(idi-1)*nal+ialal,idd,ie,iaa,izz,it+1) &
							& 	            +(1-pphi)*V0((il-1)*ntr+itr,(idi-1)*nal+ialal,idd,ie,iaa,izz,it+1) )  !Age and might go LTU
						Vc1 = ptau(it)*(pphi*     VN0((il-1)*ntr+itr,(idi-1)*nal+ialal,idd,ie,iaa,izz,it) &
							&	      +(1-pphi)*   V0((il-1)*ntr+itr,(idi-1)*nal+ialal,idd,ie,iaa,izz,it) ) + Vc1    !Don't age, maybe LTU
					else !unemployed exogenously
						Vc1 = (1._dp-ptau(it))*(pphi*     VN0((il-1)*ntr+itr,(idi-1)*nal+ialUn,idd,ie,iaa,izz,it+1) &
							& 	+(1-pphi)*( fndgrid(il,iz)*V0((il-1)*ntr+itr,(idi-1)*nal+ialal,idd,ie,iaa,izz,it+1) +&
									 (1.- fndgrid(il,iz))*VU0((il-1)*ntr+itr,(idi-1)*nal+ialUn,idd,ie,iaa,izz,it+1)  ) )  !Age and might go LTU

						Vc1 = ptau(it)*(pphi*             VN0((il-1)*ntr+itr,(idi-1)*nal+ialUn,idd,ie,iaa,izz,it) &
							  & +(1-pphi)*( fndgrid(il,iz)*V0((il-1)*ntr+itr,(idi-1)*nal+ialal,idd,ie,iaa,izz,it) +&
							  &		  (1.-fndgrid(il,iz))*VU0((il-1)*ntr+itr,(idi-1)*nal+ialUn,idd,ie,iaa,izz,it) ) ) + Vc1    !Don't age, maybe LTU
					endif

					!Vc1 = (1.-fndgrid(il,iz))*Vc1 + fndgrid(il,iz)*V0((ij-1)*ntr+itr,(idi-1)*nal+ialal,id,ie,iaa,izz,it)
					Vtest2 = Vtest2 + beta*piz(iz,izz)*pialf(ial,ialal,id)*pid_here(id,idd) *Vc1  !Probability of alpha_i X z_i draw
				enddo
				enddo
				enddo
				Vtest2 = Vtest2 + util(chere,id,iw)
				if(Vtest2>Vtest1 .or. iaa .eq. iaa0) then !first value or optimal
					apol = iaa! set the policy
					Vtest1 = Vtest2
				elseif( simp_concav .eqv. .true. ) then
					exit
				endif
			else
				exit
			endif
		enddo
		Vout = Vtest1
	end subroutine maxVU

	subroutine maxVN(il,itr,idi,ial,id,ie,ia,iz,it, VN0, VD0,VR0, V0,wagehere, iaa0,iaaA,apol,gapp_pol,gapp_dif,Vout  )
		integer, intent(in) :: il,itr,idi,ial,id,ie,ia,iz,it
		integer, intent(in) :: iaa0, iaaA
		integer, intent(out) :: apol,gapp_pol
		real(dp), intent(out) :: gapp_dif
		real(dp), intent(out) :: Vout
		real(dp), intent(in) :: VN0(:,:,:,:,:,:,:),VR0(:,:,:,:),VD0(:,:,:,:),V0(:,:,:,:,:,:,:),wagehere
		real(dp) :: Vc1,chere,Vtest1,Vtest2,Vapp,VNapp,smthV, VNhr,VRhr, VDhr, maxVNV0, EVD,Ewage,minvalVN, xihr,nuhr
		real(dp) :: pid_here(nd,nd)
		integer :: iw, iaa,ialal,izz,aapp,aNapp, ialalhr,idd

		pid_here = pid(:,:,idi,it)

		iw = 1 ! not working
		!*******************************************
		!**************Value if do not apply for DI

		Vtest1 = -1e6
		apol = iaa0
		do iaa = iaa0,iaaA
			chere = min(LTUrr*egrid(ie),LTUrr)+R*agrid(ia)-agrid(iaa)
			if(chere >0._dp) then
				Vtest2 = 0._dp
				!Continuation if do not apply for DI
				do izz = 1,nz	 !Loop over z'
				do ialal = ialL,nal !Loop over alpha_i'
				do idd = 1,nd
					if(ial == ialUn) ialalhr = ialUn
					if(ial > ialUn)  ialalhr = ialal


					VNhr    =   	VN0((il-1)*ntr+itr,(idi-1)*nal +ialalhr,idd,ie,iaa,izz,it+1)
					!maxVNV0 = max(   V0((il-1)*ntr+itr,(idi-1)*nal +ialalhr,idd,ie,iaa,izz,it+1),VNhr)
					maxVNV0 = max(   V0((il-1)*ntr+itr,(idi-1)*nal +ialal  ,idd,ie,iaa,izz,it+1),VNhr)
					if(it .ge. TT-1) then
						!can choose early retirement
						VRhr = VR0(1,id,ie,iaa)
						VNhr = (max(VNhr,VRhr)+ VNhr)*0.5_dp
						maxVNV0 = (max(maxVNV0,VRhr)+maxVNV0)*0.5_dp
					endif
					Vc1 = (1-ptau(it))*( (1-lrho*fndgrid(il,iz) )*VNhr + lrho*fndgrid(il,iz)*maxVNV0 ) !Age and might work

					VNhr    =   	VN0((il-1)*ntr+itr,(idi-1)*nal +ialalhr,idd,ie,iaa,izz,it)
					!maxVNV0 = max(   V0((il-1)*ntr+itr,(idi-1)*nal +ialalhr,idd,ie,iaa,izz,it),VNhr)
					maxVNV0 = max(   V0((il-1)*ntr+itr,(idi-1)*nal +ialal  ,idd,ie,iaa,izz,it),VNhr)

					Vc1 = Vc1+ptau(it)*((1-lrho*fndgrid(il,iz))*VNhr +lrho*fndgrid(il,iz)*maxVNV0)     !Don't age, might work

					Vtest2 = Vtest2 + beta*piz(iz,izz)*pialf(ial,ialal,id)*pid_here(id,idd) *Vc1
				enddo
				enddo
				enddo
				Vtest2 = Vtest2 + util(chere,id,iw)

				if( Vtest2 > Vtest1 .or. iaa .eq. iaa0) then
					apol = iaa
					Vtest1 = Vtest2
				elseif(simp_concav .eqv. .true.) then
					exit
				endif
			else
				exit
			endif
		enddo !iaa
		Vnapp = Vtest1
		aNapp = apol !agrid(apol)
		if((Vnapp <-1e5) .and. (verbose >0)) then
			print *, "ruh roh!"
			print *, "Vnapp, aNapp: ", Vnapp, aNapp
!~ 			write(*,*) "VD: ",id,ie,iaa,it
!~ 			write(*,*) "VN: ",il,itr,idi,ial,id,ie,ia,iz,it
		endif

		!**********Value if apply for DI


		nuhr = nud(id)
		if(it .ge. TT-2) nuhr = nud(id) + nuage
		!if(it== TT-1) nuhr = nu*(ptau(it)) !only pay nu for non-retired state

		minvalVN = minval(VN0((il-1)*ntr+itr,((idi-1)*nal+1):(idi*nal),:,:,:,:,it))
		Vtest1 = -1e6
		apol = iaa0
		EVD = 0._dp
		Ewage = 0._dp
		do iaa = iaa0,iaaA
			chere = min(LTUrr*egrid(ie),LTUrr)+R*agrid(ia)-agrid(iaa)
			if(chere >0.) then
				Vtest2 = 0.
				!Continuation if apply for DI
				do izz = 1,nz	 !Loop over z'
				do ialal = ialL,nal !Loop over alpha_i'
				do idd = 1,nd
					xihr = xifun(id,trgrid(itr),it,0._dp ) ! <- when it depended on trend
					!xihr = xifun(id, wage(trgrid(itr),alfgrid(ialal,id),idd,it), it, wmean)

					if(ial == ialUn) ialalhr = ialUn
					if(ial > ialUn)  ialalhr = ialal
					!ialalhr = ialal
					VNhr    =   	VN0((il-1)*ntr+itr,(idi-1)*nal +ialalhr,idd,ie,iaa,izz,it+1)
					!maxVNV0 = max(	 V0((il-1)*ntr+itr,(idi-1)*nal +ialalhr,idd,ie,iaa,izz,it+1), VNhr)
					maxVNV0 = max(	 V0((il-1)*ntr+itr,(idi-1)*nal +ialal  ,idd,ie,iaa,izz,it+1), VNhr)

					VDhr    = max(VD0(idd,ie,iaa,it+1),VNhr)
					Vc1 =   (1-ptau(it))*(1-xihr)*( (1-lrho*fndgrid(il,iz) )*VNhr +lrho*fndgrid(il,iz)*maxVNV0 )&
						& + (1-ptau(it))*xihr    *VDhr !Age and might go on DI

					VNhr    =   	VN0((il-1)*ntr+itr,(idi-1)*nal +ialalhr,idd,ie,iaa,izz,it)
					!maxVNV0 = max(	 V0((il-1)*ntr+itr,(idi-1)*nal +ialalhr,idd,ie,iaa,izz,it),VNhr)
					maxVNV0 = max(	 V0((il-1)*ntr+itr,(idi-1)*nal +ialal  ,idd,ie,iaa,izz,it),VNhr)

					VDhr    = max(VD0(idd,ie,iaa,it),VNhr)

					Vc1 = Vc1 +	    ptau(it)*(1-xihr)*( (1-lrho*fndgrid(il,iz))*VNhr +lrho*fndgrid(il,iz)*maxVNV0 ) &
						&     + 	ptau(it)*xihr    * VDhr     !Don't age, might go on DI

					Vtest2 = Vtest2 + beta*piz(iz,izz)*pialf(ial,ialal,id) *pid_here(id,idd) *Vc1
					if(iaa == iaa0) then
						Ewage = Ewage+     piz(iz,izz)*pialf(ial,ialal,id) *pid_here(id,idd) *wage(0._dp,alfgrid(ialal,id),idd,it)
					!	EVD   = EVD + beta*piz(iz,izz)*pialf(ial,ialal,id) *pid_here(id,idd) *VDhr
					endif
				enddo
				enddo
				enddo
				if( ial .eq. ialUn ) Ewage = min(LTUrr*egrid(ie),LTUrr)

				!Vtest2 = util(chere,id,iw) &
				!	&	+ Vtest2 - nuhr*Vtest2
				Vtest2 = util(chere*exp(-nuhr*Ewage),id,iw) &
					&	+ Vtest2
				if (Vtest2>Vtest1  .or. iaa .eq. iaa0) then
					apol = iaa
					Vtest1 = Vtest2
				elseif(simp_concav .eqv. .true.) then
					exit
				endif

			else
				exit
			endif
		enddo !iaa
		Vapp = Vtest1
		aapp = apol !agrid(apol)

		if(Vapp <-1e5) then
			write(*,*) "ruh roh!"
			write(*,*) "Vapp, aapp: ", Vapp, aapp
		endif

		!******************************************************
		!***************** Discrete choice for application
		smthV = dexp(smthV0param*(Vnapp-Vapp))/( dexp(smthV0param*(Vnapp-Vapp)) +1._dp )
		if(smthV .ge. 1 .or. smthV .le. 0 .or. isnan(smthV)) then
			if( Vapp >  Vnapp  ) smthV =0.
			if(Vnapp >=  Vapp  ) smthV =1._dp
		endif
		if (noDI .eqv. .true.) then
			smthV = 1._dp
		endif

		if (Vapp > Vnapp) then
			apol = aapp
			gapp_pol = 1
		else !Don't apply
			apol = aNapp
			gapp_pol = 0
		endif

		gapp_dif = Vapp - Vnapp
		if(verbose .gt. 4) print *, Vapp - Vnapp
		if((it == 1) .and. (ineligNoNu .eqv. .true.)) then
			Vout = smthV*Vnapp + (1._dp - smthV)*(eligY*Vapp + (1._dp- eligY)*Vnapp)
		else
			Vout = smthV*Vnapp + (1._dp - smthV)*Vapp
		endif

	end subroutine maxVN

	subroutine maxVW(il,itr,idi,ial,id,ie,ia,iz,it, VU, V0,wagehere,iee1,iee2,iee1wt,iaa0,iaaA,apol,gwork_pol,gwork_dif,Vout ,VWout )
		integer, intent(in) :: il,itr,idi,ial,id,ie,ia,iz,it
		integer, intent(in) :: iaa0, iaaA,iee1,iee2
		real(dp), intent(in) :: iee1wt,wagehere
		integer, intent(out) :: apol,gwork_pol
		real(dp), intent(out) :: gwork_dif
		real(dp), intent(out) :: Vout, VWout
		real(dp), intent(in) :: VU(:,:,:,:,:,:,:),V0(:,:,:,:,:,:,:)
		real(dp) :: Vc1,utilhere,chere,Vtest1,Vtest2,VWhere,VUhere,smthV, vL, vH, uL, uH, pid_here(nd,nd)
		integer :: iw, iaa,ialal,izz,idd,itt

		pid_here = pid(:,:,idi,it)

		iw = 2 ! working
		Vtest1= -1.e6_dp ! just to initialize, does not matter
		apol = iaa0

		!Find saving Policy
		do iaa=iaa0,iaaA
			!Continuation value if don't go on disability
			chere = wagehere+R*agrid(ia)-agrid(iaa)
			if (chere >0.) then
				Vc1 = 0.
				do izz  = 1,nz	 !Loop over z'
				do ialal = ialL,nal  !Loop over alpha_i'
				do idd  = 1,nd	 !Loop over d'
					!Linearly interpolating on e'
					vL = (1-ptau(it))*V0((il-1)*ntr+itr,(idi-1)*nal+ialal,idd,iee1,iaa,izz,it+1) &
						& +ptau(it)  *V0((il-1)*ntr+itr,(idi-1)*nal+ialal,idd,iee1,iaa,izz,it)
					vH = (1-ptau(it))*V0((il-1)*ntr+itr,(idi-1)*nal+ialal,idd,iee2,iaa,izz,it+1) &
						& +ptau(it)  *V0((il-1)*ntr+itr,(idi-1)*nal+ialal,idd,iee2,iaa,izz,it)
					!if become unemployed here -
					uL = (1-ptau(it))*VU((il-1)*ntr+itr,(idi-1)*nal+ialUn,idd,iee1,iaa,izz,it+1) &
						& +ptau(it)  *VU((il-1)*ntr+itr,(idi-1)*nal+ialUn,idd,iee1,iaa,izz,it)
					uH = (1-ptau(it))*VU((il-1)*ntr+itr,(idi-1)*nal+ialUn,idd,iee2,iaa,izz,it+1) &
						& +ptau(it)  *VU((il-1)*ntr+itr,(idi-1)*nal+ialUn,idd,iee2,iaa,izz,it)

					Vc1 = piz(iz,izz)*pialf(ial,ialal,id)*pid_here(id,idd) &
						& * (  (1.-sepgrid(il,iz))*(vH*(1._dp - iee1wt) + vL*iee1wt) &
							+  sepgrid(il,iz)*     (uH*(1._dp - iee1wt) + uL*iee1wt) ) &
						& + Vc1
				enddo
				enddo
				enddo
				if(it>1) then
					itt = 2
				else
					itt=it
				endif
				utilhere = util(chere,id,iw) - Fd(id,itt)
				Vtest2 = utilhere + beta*Vc1 ! flow utility

				if (Vtest2>Vtest1 .or. iaa .eq. iaa0 ) then  !always replace on the first, or best
					Vtest1 = Vtest2
					apol = iaa
				elseif (simp_concav .eqv. .true.) then
					! if we are imposing concavity
					exit
				endif
			else
				exit
			endif
		enddo	!iaa

		VWhere = Vtest1
		VUhere = VU((il-1)*ntr+itr,(idi-1)*nal+ial,id,ie,ia,iz,it)

		!------------------------------------------------!
		!Calculate V with solved vals of VW and VU -- i.e. can quit into unemployment
		!------------------------------------------------!

		if (VWhere>VUhere) then
			gwork_pol = 1
		else
			gwork_pol = 0
		endif
		smthV = dexp( smthV0param *(VWhere -VUhere) ) &
			& /( dexp(smthV0param *(VWhere -VUhere)) + 1._dp )
		if (smthV <1e-5 .or. smthV>.999999 .or. isnan(smthV) ) then
			if(VWhere >=  VUhere)	smthV = 1._dp
			if(VWhere  <  VUhere)	smthV = 0.
		endif
		Vout = smthV*VWhere + (1._dp-smthV)*VUhere

		gwork_dif = VWhere - VUhere
		VWout = VWhere

	end subroutine maxVW

	subroutine sol(val_funs, pol_funs)

		implicit none

		type(val_struct), intent(inout), target :: val_funs
		type(pol_struct), intent(inout), target :: pol_funs

	!************************************************************************************************!
	! Counters and Indicies
	!************************************************************************************************!

		integer  :: i, j, ia, ie, id, it, ga,gw, anapp,aapp, apol, itr, ial, il , idi,  &
			    iee1, iee2, iz, iw,wo, iter,npara,ipara, iaa_k,iaa0,iaaA,iaN, iter_timeout
		integer, dimension(5) :: maxer_i

		integer :: aa_l(na), aa_u(na), ia_o(na),aa_m(na)

		logical :: ptrsucces
		!************************************************************************************************!
		! Value Functions- Stack z-risk j and indiv. exposure beta_i
		!************************************************************************************************!
		real(dp)  	  	:: Vtest1, maxer_v, smthV, wagehere, iee1wt, gadif,gwdif
		real(dp), allocatable 	:: maxer_hist(:)
		real(dp), allocatable	:: verloc(:,:,:,:,:)
		real(dp), allocatable :: VR0(:,:,:,:), &			!Retirement
					VD0(:,:,:,:), &			!Disabled
					VN0(:,:,:,:,:,:,:), &	!Long-term Unemployed
					VW0(:,:,:,:,:,:,:), &	!Working
					VU0(:,:,:,:,:,:,:), &	!Unemployed
					V0(:,:,:,:,:,:,:)	!Participant

		real(dp), pointer ::	VR(:,:,:,:), &			!Retirement
					VD(:,:,:,:), &			!Disabled
					VN(:,:,:,:,:,:,:), &	!Long-term Unemployed
					VW(:,:,:,:,:,:,:), &	!Working
					VU(:,:,:,:,:,:,:), &	!Unemployed
					V(:,:,:,:,:,:,:)	!Participant

		real(dp), pointer ::	gapp_dif(:,:,:,:,:,:,:), gwork_dif(:,:,:,:,:,:,:) ! latent value of work/apply

		integer, pointer ::	aR(:,:,:,:), aD(:,:,:,:), aU(:,:,:,:,:,:,:), &
					aN(:,:,:,:,:,:,:), aW(:,:,:,:,:,:,:)
		integer, pointer ::	gapp(:,:,:,:,:,:,:), &
					gwork(:,:,:,:,:,:,:)

		!************************************************************************************************!
		! Other
		!************************************************************************************************!
			real(dp)	:: junk,sumer,maxer, eprime, emin, emax, VWhere, V_m(na)
		!************************************************************************************************!

		!************************************************************************************************!
		! Allocate matrices
		!************************************************************************************************!
		! (disability extent, earn hist, assets)
		allocate(VR0(2,nd,ne,na))
		allocate(VD0(nd,ne,na,TT))
		allocate(VN0(nl*ntr,ndi*nal,nd,ne,na,nz,TT))
		allocate(VU0(nl*ntr,ndi*nal,nd,ne,na,nz,TT))
		allocate(VW0(nl*ntr,ndi*nal,nd,ne,na,nz,TT))
		allocate(V0(nl*ntr,ndi*nal,nd,ne,na,nz,TT))
		! there must be a way to use pointers, but it doesn't seem to work
		VR => val_funs%VR
		aR => pol_funs%aR
		VD => val_funs%VD
		aD => pol_funs%aD
		VN => val_funs%VN
		VU => val_funs%VU
		VW => val_funs%VW
		V  => val_funs%V
		aN => pol_funs%aN
		aW => pol_funs%aW
		aU => pol_funs%aU
		gwork => pol_funs%gwork
		gapp => pol_funs%gapp

		gapp_dif => pol_funs%gapp_dif
		gwork_dif => pol_funs%gwork_dif

		allocate(verloc(na,nz,ne,nd,nal))
		allocate(maxer_hist(maxiter))
		maxer_hist = 0._dp

		emin = minval(egrid)
		emax = maxval(egrid)

		ptrsucces = associated(VR,val_funs%VR)

		simp_concav = .false.
		!************************************************************************************************!
		! Caculate things that are independent of occupation/person type
		!	1) Value of Retired:  VR(it0,d,e,a)
		!	2) Value of Disabled: VD(d,e,a)

	!1) Calculate Value of Retired: VR(it0,d,e,a)
		!d in{1,2,3}  : disability extent
		!e inR+       :	earnings index
		!a inR+	      : asset holdings

		Vevals = 0
		!VFI with good guess
		!Initialize
		iw=1
		do it=1,2
		do id=1,nd
		do ie=1,ne
		do ia=1,na
			VR0(it,id,ie,ia) = util(SSret(egrid(ie),agegrid(it+TT-2))+R*agrid(ia),id,iw)* (1._dp/(1._dp-beta*ptau(TT)))
		enddo
		enddo
		enddo
		enddo
		if(print_lev >3) then
			i = 1
			call vec2csv(VR0(i+1,i,i,:),"VR0.csv",0)
		endif
		iter = 1
!		simp_concav = .true. ! use simple concavity here
		do while (iter<=maxiter)
			sumer = 0
			do it =1,2
			do id =1,nd
		  	do ie =1,ne

				iaN =0
				ia = 1
				iaa0 = 1
				iaaA = na
				call maxVR(it,id,ie,ia, VR0,iaa0, iaaA, apol,Vtest1)
				VR(it,id,ie,ia) = Vtest1
				aR(it,id,ie,ia) = apol !agrid(apol)
				iaN = iaN+1
				ia_o(iaN) = ia

				ia = na
				iaa0 = aR(it,id,ie,1)
				iaaA = na
				call maxVR(it,id,ie,ia, VR0,iaa0, iaaA, apol,Vtest1)
				VR(it,id,ie,ia) = Vtest1
				aR(it,id,ie,ia) = apol !agrid(apol)
				iaN = iaN+1
				ia_o(iaN) = ia


				iaa_k = 1
				aa_l(iaa_k) = 1
				aa_u(iaa_k) = na

				!main loop (step 1 of Gordon & Qiu completed)
				outerVR: do
					!Expand list (step 2 of Gordon & Qiu)
					do
						if(aa_u(iaa_k) == aa_l(iaa_k)+1) exit
						iaa_k = iaa_k+1
						aa_l(iaa_k) = aa_l(iaa_k-1)
						aa_u(iaa_k) = (aa_l(iaa_k-1)+aa_u(iaa_k-1))/2
						!search given ia from iaa0 to iaaA
						ia = aa_u(iaa_k)
						iaa0 = aR(it,id,ie, aa_l(iaa_k-1) )
						iaaA = aR(it,id,ie, aa_u(iaa_k-1) )
						call maxVR(it,id,ie,ia, VR0,iaa0, iaaA, apol,Vtest1)
						VR(it,id,ie,ia) = Vtest1
						aR(it,id,ie,ia) = apol
						iaN = iaN+1
						ia_o(iaN) = ia

					enddo
					! Move to a higher interval or stop (step 3 of Gordon & Qiu)
					do
						if(iaa_k==1) exit outerVR
						if( aa_u(iaa_k)/= aa_u(iaa_k - 1) ) exit
						iaa_k = iaa_k -1
					end do
					! more to the right subinterval
					aa_l(iaa_k) = aa_u(iaa_k)
					aa_u(iaa_k) = aa_u(iaa_k-1)
				end do outerVR

				do ia=1,na
					sumer = sumer+ (VR(it,id,ie,ia)-VR0(it,id,ie,ia))**2
				enddo
				VR0(it,id,ie,:) = VR(it,id,ie,:)
			enddo !it
			enddo !ie
			enddo !id
			if (sumer < Vtol) then
				exit	!Converged
			else
				iter=iter+1
				if(print_lev >3) then
					i = 1
					call veci2csv(aR(2,i,i,:),"aR.csv",0)
					call vec2csv(VR(2,i,i,:),"VR.csv",0)
				endif
			endif
		enddo ! iteration iter

		if(sumer >= Vtol) then
			print *, "VR did not converge and distance is", sumer
		endif

		i = 1
		do it =1,2
		do id =2,nd
		do ie =1,ne
		do ia =1,na
			VR(it,id,ie,ia) = VR(it,i,ie,ia)*(dexp(theta*dble(id-1)))**(1-gam)
			aR(it,id,ie,ia) = aR(it,i,ie,ia)
		enddo
		enddo
		enddo
		enddo

		if (print_lev > 2) then
			wo =0
			do id=1,nd
			do ie=1,ne
				call veci2csv(aR(2,i,i,:),"aR.csv",wo)
				call vec2csv(VR(2,i,i,:),"VR.csv",wo)
				if(wo .eq. 0) wo = 1
			enddo
			enddo
		endif

		!----------------------------------------------------------!
		!Set value at t=TT to be VR in all other V-functions
		!----------------------------------------------------------!
		VD(:,:,:,TT) = VR(2,:,:,:)
		VD0(:,:,:,TT) = VR(2,:,:,:)
		VD0(:,:,:,TT-1) = VD0(:,:,:,TT) ! this is necessary for initialization given stochastic aging

	!----------------------------------------------------------------!
	!2) Calculate Value of Disabled: VD(d,e,a,t)
		!d in{1,2,3}  	   :disability extent
		!e inR+       	   :earnings index
		!a inR+	      	   :asset holdings
		!t in[1,2...TT-1]  :age

		!simp_concav = .true.
		npara = nd*ne
		!Work backwards from TT
		do it = TT-1,1,-1
			!Guess will be value at t+1
			VD0(:,:,:,it) = VD(:,:,:,it+1)
			iw = 1 ! not working

			do id =1,nd

				if(print_lev >3) then
					call mat2csv(VD0(id,:,:,it),"VD0.csv",0)
				endif
				!Loop over earnings index
				do ie=1,ne
					!Loop to find V(..,it) as fixed point
					iter=1
					do while (iter<=maxiter)
						sumer = 0

						iaN = 0
						ia = 1
						iaa0 = 1
						iaaA = na
						call maxVD(id,ie,ia,it, VD0, iaa0, iaaA, apol,Vtest1)
						VD(id,ie,ia,it) = Vtest1
						aD(id,ie,ia,it) = apol !agrid(apol)
						iaN = iaN+1
						ia_o(iaN) = ia


						ia = na
						iaa0 = aD(id,ie,1,it)
						iaaA = na
						call maxVD(id,ie,ia,it, VD0, iaa0, iaaA, apol,Vtest1)
						VD(id,ie,ia,it) = Vtest1
						aD(id,ie,ia,it) = apol !agrid(apol)
						iaN = iaN+1
						ia_o(iaN) = ia


						iaa_k = 1
						aa_l(iaa_k) = 1
						aa_u(iaa_k) = na

						!main loop (step 1 of Gordon & Qiu completed)
						outerVD: do
							!Expand list (step 2 of Gordon & Qiu)
							do
								if(aa_u(iaa_k) == aa_l(iaa_k)+1) exit
								iaa_k = iaa_k+1
								aa_l(iaa_k) = aa_l(iaa_k-1)
								aa_u(iaa_k) = (aa_l(iaa_k-1)+aa_u(iaa_k-1))/2
								!search given ia from iaa0 to iaaA
								ia = aa_u(iaa_k)
								iaa0 = aD(id,ie, aa_l(iaa_k-1),it)
								iaaA = aD(id,ie, aa_u(iaa_k-1),it)
								call maxVD(id,ie,ia,it, VD0, iaa0, iaaA, apol,Vtest1)
								VD(id,ie,ia,it) = Vtest1
								aD(id,ie,ia,it) = apol
								iaN = iaN+1
								ia_o(iaN) = ia

							enddo
							! Move to a higher interval or stop (step 3 of Gordon & Qiu)
							do
								if(iaa_k==1) exit outerVD
								if( aa_u(iaa_k)/= aa_u(iaa_k - 1) ) exit
								iaa_k = iaa_k -1
							end do
							! more to the right subinterval
							aa_l(iaa_k) = aa_u(iaa_k)
							aa_u(iaa_k) = aa_u(iaa_k-1)
						end do outerVD

						do ia=1,na
							sumer = sumer+ (VD(id,ie,ia,it)-VD0(id,ie,ia,it))**2
						enddo

						if(print_lev >3) then
							wo =0
							call veci2csv(aD(id,ie,:,it),"aD.csv",wo)
							call vec2csv(VD(id,ie,:,it),"VD.csv",wo)
						endif
						if (sumer < Vtol) then
							exit	!Converged
						endif
						VD0(id,ie,:,it) = VD(id,ie,:,it)	!New guess
						iter=iter+1

					enddo	!iter: V-iter loop
				enddo	!ie

			enddo !id = 1,nd
		enddo	!t loop, going backwards



		VD0 = VD
		if (print_lev >= 2) then
			wo =0
			do id =1,nd
			do ie =1,ne
				call mati2csv(aD(id,ie,:,:),"aD.csv",wo)
				call mat2csv(VD(id,ie,:,:),"VD.csv",wo)
				if(wo == 0) wo =1
			enddo
			enddo
		endif


	!************************************************************************************************!
	!3) Calculate V= max(VW,VN); requires calculating VW and VN

	! initialize

		do id=1,nd
		do ie=1,ne
		do ia=1,na
			VW (:,:,id,ie,ia,:,TT) = VR(2,id,ie,ia)
			VW0(:,:,id,ie,ia,:,TT) = VR(2,id,ie,ia)
			VN (:,:,id,ie,ia,:,TT) = VR(2,id,ie,ia)
			VN0(:,:,id,ie,ia,:,TT) = VR(2,id,ie,ia)
			VU (:,:,id,ie,ia,:,TT) = VR(2,id,ie,ia)
			VU0(:,:,id,ie,ia,:,TT) = VR(2,id,ie,ia)
			V  (:,:,id,ie,ia,:,TT) = VR(2,id,ie,ia)
			V0 (:,:,id,ie,ia,:,TT) = VR(2,id,ie,ia)
		enddo
		enddo
		enddo



		! Begin loop over occupations
		do il = 1,nl
		! And individual disability type
		do idi = 1,ndi

		!************************************************************************************************!
			!Work Backwards TT-1,TT-2...1!
		do it=TT-1,1,-1
			!----Initialize---!

			if( il==1 .and. idi>1) then
				do ial=1,nal
				do itr = 1,ntr
				do id =1,nd
				do ie =1,ne
				do iz =1,nz
				do ia =1,na

				!Guess once, then use last occupation as guess (order occupations intelligently)
				! for it = 1, should be TT-1+1 =TT -> VU,Vw,VN = VR
					VW0((il-1)*ntr+itr,(idi-1)*nal+ial,id,ie,ia,iz,it) = VW0((il-1)*ntr+itr,(idi-2)*nal+ial,id,ie,ia,iz,it)
					VU0((il-1)*ntr+itr,(idi-1)*nal+ial,id,ie,ia,iz,it) = VU0((il-1)*ntr+itr,(idi-2)*nal+ial,id,ie,ia,iz,it)
					VN0((il-1)*ntr+itr,(idi-1)*nal+ial,id,ie,ia,iz,it) = VN0((il-1)*ntr+itr,(idi-2)*nal +ial,id,ie,ia,iz,it)
					V0 ((il-1)*ntr+itr,(idi-1)*nal+ial,id,ie,ia,iz,it) = V0 ((il-1)*ntr+itr,(idi-2)*nal+ial,id,ie,ia,iz,it)

				enddo	!ia
				enddo	!iz
				enddo	!ie
				enddo 	!id
				enddo 	!itr
				enddo	!ial
			elseif( il>1 ) then
				do ial=1,nal
				do itr = 1,ntr
				do id =1,nd
				do ie =1,ne
				do iz =1,nz
				do ia =1,na

				!Guess once, then use last occupation as guess (order occupations intelligently)
				! for it = 1, should be TT-1+1 =TT -> VU,Vw,VN = VR
					VW0((il-1)*ntr+itr,(idi-1)*nal+ial,id,ie,ia,iz,it) = VW0((il-2)*ntr+itr,(idi-1)*nal+ial,id,ie,ia,iz,it)
					VU0((il-1)*ntr+itr,(idi-1)*nal+ial,id,ie,ia,iz,it) = VU0((il-2)*ntr+itr,(idi-1)*nal+ial,id,ie,ia,iz,it)
					VN0((il-1)*ntr+itr,(idi-1)*nal+ial,id,ie,ia,iz,it) = VN0((il-2)*ntr+itr,(idi-1)*nal +ial,id,ie,ia,iz,it)
					V0 ((il-1)*ntr+itr,(idi-1)*nal+ial,id,ie,ia,iz,it) = V0 ((il-2)*ntr+itr,(idi-1)*nal+ial,id,ie,ia,iz,it)

				enddo	!ia
				enddo	!iz
				enddo	!ie
				enddo 	!id
				enddo	!itr
				enddo	!ial
			else  ! should be only (il ==1 .and. idi == 1)  then
				do ial=1,nal
				do itr = 1,ntr
				do id =1,nd
				do ie =1,ne
				do iz =1,nz
				do ia =1,na

				!Guess once, then use next period same occupation/beta as guess
				! for it = 1, should be TT-1+1 =TT -> VU,Vw,VN = VR
					VW0((il-1)*ntr+itr,(idi-1)*nal+ial,id,ie,ia,iz,it)  = VW ((il-1)*ntr+itr,(idi-1)*nal+ial,id,ie,ia,iz,it+1)
					VU0((il-1)*ntr+itr,(idi-1)*nal+ial,id,ie,ia,iz,it)  = VU ((il-1)*ntr+itr,(idi-1)*nal+ial,id,ie,ia,iz,it+1)
					VN0((il-1)*ntr+itr,(idi-1)*nal+ial,id,ie,ia,iz,it)  = VN ((il-1)*ntr+itr,(idi-1)*nal+ial,id,ie,ia,iz,it+1)
					V0 ((il-1)*ntr+itr,(idi-1)*nal+ial,id,ie,ia,iz,it)  = VW0((il-1)*ntr+itr,(idi-1)*nal+ial,id,ie,ia,iz,it)

				enddo	!ia
				enddo	!iz
				enddo	!ie
				enddo 	!id
				enddo	!itr
				enddo	!ial

			endif

			!***********************************************************************************************************
			!Loop over V=max(VU,VW)
			iter=1
			iter_timeout = 0
			smthV0param  =1._dp ! will tighten this down
			do while (iter<=maxiter)
				verloc = 0.
				sumer = 0.	!Use to calc |V-V0|<eps

				simp_concav = .false.

			!------------------------------------------------!
			!Solve VU given guesses on VW, VN, VU and implied V
			!------------------------------------------------!
				sumer = 0.
				wo = 0
				npara = nal*ntr*nd*ne*nz
			!$OMP  parallel do reduction(+:sumer)&
			!$OMP& private(ial,id,ie,iz,iw,itr,apol,iaa0,iaaA,ia,ipara,Vtest1,aa_m,V_m,aa_l,aa_u,iaa_k,iaN,ia_o)
				do ipara = 1,npara
					iz = mod(ipara-1,nz)+1
					ie = mod(ipara-1,nz*ne)/nz + 1
					id = mod(ipara-1,nz*ne*nd)/(nz*ne) +1
					itr= mod(ipara-1,nz*ne*nd*ntr)/(nz*ne*nd) +1
					ial= mod(ipara-1,nz*ne*nd*ntr*nal)/(nz*ne*nd*ntr)+1

					! search over ia
						iaN =0
						ia = 1
						iaa0 = 1
						iaaA = na
						call maxVU(il,itr,idi,ial,id,ie,ia,iz,it, VU0,VN0,V0, iaa0,iaaA,apol,Vtest1)
						V_m(ia) = Vtest1
						aa_m(ia) = apol !agrid(apol)
						iaN = iaN+1
						ia_o(iaN) = ia

						ia = na
						iaa0 = aa_m(1)
						iaaA = na
						call maxVU(il,itr,idi,ial,id,ie,ia,iz,it, VU0,VN0,V0, iaa0,iaaA,apol,Vtest1)
						V_m(ia) = Vtest1
						aa_m(ia) = apol !agrid(apol)
						iaN = iaN+1
						ia_o(iaN) = ia


						iaa_k = 1
						aa_l(iaa_k) = 1
						aa_u(iaa_k) = na

						!main loop (step 1 of Gordon & Qiu completed)
						outerVU: do
							!Expand list (step 2 of Gordon & Qiu)
							do
								if(aa_u(iaa_k) == aa_l(iaa_k)+1) exit
								iaa_k = iaa_k+1
								aa_l(iaa_k) = aa_l(iaa_k-1)
								aa_u(iaa_k) = (aa_l(iaa_k-1)+aa_u(iaa_k-1))/2
								!search given ia from iaa0 to iaaA
								ia = aa_u(iaa_k)
								iaa0 = aa_m( aa_l(iaa_k-1) )
								iaaA = aa_m( aa_u(iaa_k-1) )
								call maxVU(il,itr,idi,ial,id,ie,ia,iz,it, VU0,VN0,V0, iaa0,iaaA,apol,Vtest1)
								V_m(ia) = Vtest1
								aa_m(ia) = apol !agrid(apol)
								iaN = iaN+1
								ia_o(iaN) = ia

							enddo
							! Move to a higher interval or stop (step 3 of Gordon & Qiu)
							do
								if(iaa_k==1) exit outerVU
								if( aa_u(iaa_k)/= aa_u(iaa_k - 1) ) exit
								iaa_k = iaa_k -1
							end do
							! more to the right subinterval
							aa_l(iaa_k) = aa_u(iaa_k)
							aa_u(iaa_k) = aa_u(iaa_k-1)
						end do outerVU

						aU((il-1)*ntr+itr,(idi-1)*nal+ial,id,ie,:,iz,it) = aa_m
						VU((il-1)*ntr+itr,(idi-1)*nal+ial,id,ie,:,iz,it) = V_m

						if((iz>nz .or. ie>ne .or. id>nd .or. ial>nal) .and. verbose > 2) then
							print *, "ipara is not working right", iz,ie,id,ial
						endif
						Vtest1 = sum( (VU((il-1)*ntr+itr,(idi-1)*nal+ial,id,ie,:,iz,it) &
							& - VU0((il-1)*ntr+itr,(idi-1)*nal+ial,id,ie,:,iz,it))**2)
						sumer = Vtest1 + sumer
						if(print_lev >3) then
							call veci2csv(ia_o,"ia_o_VU.csv",wo)
							if(wo == 0) wo =1
						endif

				enddo !ipara
			!$OMP END PARALLEL do
				!update VU0
				do ial=1,nal	!Loop over alpha (ai)
				do itr=1,ntr	!Loop over trend
				do id=1,nd	!Loop over disability index
				do ie=1,ne	!Loop over earnings index
				do iz=1,nz	!Loop over TFP
					do ia =1,na
						VU0((il-1)*ntr+itr,(idi-1)*nal+ial,id,ie,ia,iz,it) = VU((il-1)*ntr+itr,(idi-1)*nal+ial,id,ie,ia,iz,it)
					enddo	!ia
				enddo !ial <- these are not in order
				enddo !itr
				enddo !id
				enddo !ie
				enddo !iz
				if (print_lev > 3) then
					wo = 0
					itr = tri0
					do ial=1,nal	!Loop over alpha (ai)
					do id=1,nd	!Loop over disability index
					do ie=1,ne	!Loop over earnings index
					do iz=1,nz	!Loop over TFP
						call veci2csv(aU((il-1)*ntr+itr,(idi-1)*nal+ial,id,ie,:,iz,it),"aU.csv",wo)
						call vec2csv(VU((il-1)*ntr+itr,(idi-1)*nal+ial,id,ie,:,iz,it),"VU.csv",wo)
						if(wo == 0) wo = 1
					enddo
					enddo
					enddo
					enddo
				endif


			!------------------------------------------------!
			!Solve VN given guesses on VW, VN, and implied V
			!------------------------------------------------!
				sumer = 0._dp
				wo = 0
				aa_u = 0
				npara = nal*ntr*nd*ne*nz
			!$OMP  parallel do reduction(+:sumer)&
			!$OMP& private(ipara,ial,id,ie,iz,itr,apol,ga,gadif,ia,iaa0,iaaA,iaa_k,aa_l,aa_u,aa_m,V_m,Vtest1,wagehere,iaN,ia_o)
				do ipara = 1,npara
					iz = mod(ipara-1,nz)+1
					ie = mod(ipara-1,nz*ne)/nz + 1
					id = mod(ipara-1,nz*ne*nd)/(nz*ne) +1
					itr= mod(ipara-1,nz*ne*nd*ntr)/(nz*ne*nd) +1
					ial= mod(ipara-1,nz*ne*nd*ntr*nal)/(nz*ne*nd*ntr)+1

					wagehere = wage(trgrid(itr),alfgrid(ial,id),id,it)
					!----------------------------------------------------------------
					!Loop over current state: assets
					iaN=0
					ia = 1
					iaa0 = 1
					iaaA = na
					call maxVN(il,itr,idi,ial,id,ie,ia,iz,it, VN0, VD0,VR0,V0,wagehere,iaa0,iaaA,apol,ga,gadif,Vtest1)
					aa_m(ia) = apol
					V_m(ia) = Vtest1
					gapp((il-1)*ntr+itr,(idi-1)*nal+ial,id,ie,ia,iz,it) = ga
					gapp_dif((il-1)*ntr+itr,(idi-1)*nal+ial,id,ie,ia,iz,it) = gadif
					iaN = iaN+1
					ia_o(iaN) = ia


					ia = na
					iaa0 = aa_m(1)
					iaaA = na
					call maxVN(il,itr,idi,ial,id,ie,ia,iz,it, VN0, VD0,VR0,V0,wagehere,iaa0,iaaA,apol,ga,gadif,Vtest1)
					V_m(ia) = Vtest1
					aa_m(ia) = apol !agrid(apol)
					gapp((il-1)*ntr+itr,(idi-1)*nal+ial,id,ie,ia,iz,it) = ga
					gapp_dif((il-1)*ntr+itr,(idi-1)*nal+ial,id,ie,ia,iz,it) = gadif
					iaN = iaN+1
					ia_o(iaN) = ia


					iaa_k = 1
					aa_l(iaa_k) = 1
					aa_u(iaa_k) = na

					outerVN: do
						do
							if(aa_u(iaa_k) == aa_l(iaa_k)+1) exit
							iaa_k = iaa_k+1
							aa_l(iaa_k) = aa_l(iaa_k-1)
							aa_u(iaa_k) = (aa_l(iaa_k-1)+aa_u(iaa_k-1))/2
							!search given ia from iaa0 to iaaA
							ia = aa_u(iaa_k)
							iaa0 = aa_m( aa_l(iaa_k-1) )
							iaaA = aa_m( aa_u(iaa_k-1) )
							call maxVN(il,itr,idi,ial,id,ie,ia,iz,it, VN0,VD0,VR0,V0,wagehere,iaa0,iaaA,apol,ga,gadif,Vtest1)
							V_m(ia) = Vtest1
							aa_m(ia) = apol !agrid(apol)
							gapp((il-1)*ntr+itr,(idi-1)*nal+ial,id,ie,ia,iz,it) = ga
							gapp_dif((il-1)*ntr+itr,(idi-1)*nal+ial,id,ie,ia,iz,it) = gadif
							iaN = iaN+1
							ia_o(iaN) = ia
						enddo
						do
							if(iaa_k==1) exit outerVN
							if( aa_u(iaa_k)/= aa_u(iaa_k - 1) ) exit
							iaa_k = iaa_k -1
						end do
						aa_l(iaa_k) = aa_u(iaa_k)
						aa_u(iaa_k) = aa_u(iaa_k-1)
					end do outerVN

					aN((il-1)*ntr+itr,(idi-1)*nal+ial,id,ie,:,iz,it) = aa_m
					VN((il-1)*ntr+itr,(idi-1)*nal+ial,id,ie,:,iz,it) = V_m

					Vtest1 = sum((VN((il-1)*ntr+itr,(idi-1)*nal+ial,id,ie,:,iz,it) &
						& - VN0((il-1)*ntr+itr,(idi-1)*nal+ial,id,ie,:,iz,it))**2)
					sumer = Vtest1 + sumer
					if(print_lev >3) then
						call veci2csv(ia_o,"ia_o_VN.csv", wo)
						if(wo == 0) wo = 1
					endif
				enddo !ipara
			!$OMP END PARALLEL do

				!------------------------------------------------!
				! Done making VN

				if (print_lev >3) then
					wo = 0
					itr = tri0
					do ial=1,nal	!Loop over alpha (ai)
					do ie=1,ne	!Loop over earnings index
					do iz=1,nz	!Loop over TFP
						! matrix in disability index and assets
						call mat2csv(VN((il-1)*ntr+itr,(idi-1)*nal+ial,:,ie,:,iz,it) ,"VN_it.csv",wo)
						call mati2csv(aN((il-1)*ntr+itr,(idi-1)*nal+ial,:,ie,:,iz,it) ,"aN_it.csv",wo)
						call mat2csv(VN((il-1)*ntr+itr,(idi-1)*nal+ial,:,ie,:,iz,it) ,"VU_it.csv",wo)
						call mati2csv(aN((il-1)*ntr+itr,(idi-1)*nal+ial,:,ie,:,iz,it) ,"aU_it.csv",wo)
						call mati2csv(gapp((il-1)*ntr+itr,(idi-1)*nal+ial,:,ie,:,iz,it) ,"gapp_it.csv",wo)
						call mat2csv(gapp_dif((il-1)*ntr+itr,(idi-1)*nal+ial,:,ie,:,iz,it) ,"gapp_dif_it.csv",wo)

						if(wo == 0 ) wo =1
					enddo !iz
					enddo !id
					enddo !ial
				endif

				!update VN0
				do ial=1,nal	!Loop over alpha (ai)
				do itr=1,ntr	!loop over trend
				do id=1,nd	!Loop over disability index
				do ie=1,ne	!Loop over earnings index
				do iz=1,nz	!Loop over TFP
					do ia =1,na
						VN0((il-1)*ntr+itr,(idi-1)*nal+ial,id,ie,ia,iz,it) = VN((il-1)*ntr+itr,(idi-1)*nal+ial,id,ie,ia,iz,it)
					enddo	!ia
				enddo !ial
				enddo !itr
				enddo !id
				enddo !ie
				enddo !iz

				!------------------------------------------------!
				!Solve VW given guesses on VW, VN, and implied V
				!------------------------------------------------!
				sumer = 0._dp
				maxer = 0._dp
				maxer_v = 0._dp
				wo = 0
				npara = nal*ntr*nd*ne*nz
			!$OMP   parallel do reduction(+:sumer) &
			!$OMP & private(ipara,ial,id,ie,iz,itr,apol,eprime,wagehere,iee1,iee2,iee1wt,ia,iaa0,iaaA,aa_l,aa_u,iaa_k,ia_o,iaN,Vtest1,VWhere,gwdif,gw)
				do ipara = 1,npara
					iz = mod(ipara-1,nz)+1
					ie = mod(ipara-1,nz*ne)/nz + 1
					id = mod(ipara-1,nz*ne*nd)/(nz*ne) +1
					itr= mod(ipara-1,nz*ne*nd*ntr)/(nz*ne*nd) +1
					ial= mod(ipara-1,nz*ne*nd*ntr*nal)/(nz*ne*nd*ntr)+1

					!Earnings evolution independent of choices.
					wagehere = wage(trgrid(itr),alfgrid(ial,id),id,it)
					eprime = Hearn(it,ie,wagehere)
					!linear interpolate for the portion that blocks off bounds on assets
					if((eprime > emin) .and. (eprime < emax)) then  ! this should be the same as if(eprime > minval(egrid) .and. eprime < maxval(egrid))
						iee1 = ne
						do while( (eprime < egrid(iee1)) .and. (iee1>= 1))
							iee1 = iee1 -1
						enddo
						iee2 = min(ne,iee1+1)
						iee1wt = (egrid(iee2)-eprime)/(egrid(iee2)-egrid(iee1))
					elseif( eprime <= emin  ) then
						iee1wt = 1._dp
						iee1 = 1
						iee2 = 1
					else
						iee1wt = 0.
						iee1 = ne
						iee2 = ne
					endif

					!----------------------------------------------------------------
					!Loop over current state: assets
					!----------------------------------------------------------------
					if(ial>1) then
						iaN = 0
						ia = 1
						iaa0 = 1
						iaaA = na
						call maxVW(il,itr,idi,ial,id,ie,ia,iz,it, VU, V0,wagehere,iee1,iee2,iee1wt, &
							& iaa0,iaaA,apol,gw,gwdif,Vtest1 ,VWhere )
						V	((il-1)*ntr+itr,(idi-1)*nal+ial,id,ie,ia,iz,it) = Vtest1
						VW	((il-1)*ntr+itr,(idi-1)*nal+ial,id,ie,ia,iz,it) = VWhere
						gwork	((il-1)*ntr+itr,(idi-1)*nal+ial,id,ie,ia,iz,it) = gw
						gwork_dif((il-1)*ntr+itr,(idi-1)*nal+ial,id,ie,ia,iz,it) = gwdif
						aW	((il-1)*ntr+itr,(idi-1)*nal+ial,id,ie,ia,iz,it) = apol
						iaN = iaN+1
						ia_o(iaN) = ia

						ia = na
						iaa0 = aW((il-1)*ntr+itr,(idi-1)*nal+ial,id,ie,1,iz,it)
						iaaA = na
						call maxVW(il,itr,idi,ial,id,ie,ia,iz,it, VU, V0,wagehere,iee1,iee2,iee1wt, &
							& iaa0,iaaA,apol,gw,gwdif,Vtest1 ,VWhere )
						V	((il-1)*ntr+itr,(idi-1)*nal+ial,id,ie,ia,iz,it) = Vtest1
						VW	((il-1)*ntr+itr,(idi-1)*nal+ial,id,ie,ia,iz,it) = VWhere
						gwork	((il-1)*ntr+itr,(idi-1)*nal+ial,id,ie,ia,iz,it) = gw
						gwork_dif((il-1)*ntr+itr,(idi-1)*nal+ial,id,ie,ia,iz,it) = gwdif
						aW	((il-1)*ntr+itr,(idi-1)*nal+ial,id,ie,ia,iz,it) = apol
						iaN = iaN+1
						ia_o(iaN) = ia

						iaa_k = 1
						aa_l(iaa_k) = 1
						aa_u(iaa_k) = na

						outerVW: do
							do
								if(aa_u(iaa_k) == aa_l(iaa_k)+1) exit
								iaa_k = iaa_k+1
								aa_l(iaa_k) = aa_l(iaa_k-1)
								aa_u(iaa_k) = (aa_l(iaa_k-1)+aa_u(iaa_k-1))/2
								!search given ia from iaa0 to iaaA
								ia = aa_u(iaa_k)
								iaa0 = aW((il-1)*ntr+itr,(idi-1)*nal+ial,id,ie,  aa_l(iaa_k-1)  ,iz,it)
								iaaA = aW((il-1)*ntr+itr,(idi-1)*nal+ial,id,ie,  aa_u(iaa_k-1)  ,iz,it)
								call maxVW(il,itr,idi,ial,id,ie,ia,iz,it, VU, V0,wagehere,iee1,iee2,iee1wt, &
									& iaa0,iaaA,apol,gw,gwdif,Vtest1 ,VWhere )
								V	((il-1)*ntr+itr,(idi-1)*nal+ial,id,ie,ia,iz,it) = Vtest1
								VW	((il-1)*ntr+itr,(idi-1)*nal+ial,id,ie,ia,iz,it) = VWhere
								gwork	((il-1)*ntr+itr,(idi-1)*nal+ial,id,ie,ia,iz,it) = gw
								gwork_dif((il-1)*ntr+itr,(idi-1)*nal+ial,id,ie,ia,iz,it)= gwdif
								aW	((il-1)*ntr+itr,(idi-1)*nal+ial,id,ie,ia,iz,it) = apol
								iaN = iaN+1
								ia_o(iaN) = ia
							enddo
							do
								if(iaa_k==1) exit outerVW
								if( aa_u(iaa_k)/= aa_u(iaa_k - 1) ) exit
								iaa_k = iaa_k -1
							end do
							aa_l(iaa_k) = aa_u(iaa_k)
							aa_u(iaa_k) = aa_u(iaa_k-1)
						end do outerVW
					else
						do ia=1,na
							V	((il-1)*ntr+itr,(idi-1)*nal+ial,id,ie,ia,iz,it) = VU((il-1)*ntr+itr,(idi-1)*nal+ial,id,ie,ia,iz,it)
							VW	((il-1)*ntr+itr,(idi-1)*nal+ial,id,ie,ia,iz,it) = VU((il-1)*ntr+itr,(idi-1)*nal+ial,id,ie,ia,iz,it)
							gwork	((il-1)*ntr+itr,(idi-1)*nal+ial,id,ie,ia,iz,it) = 0
							gwork_dif((il-1)*ntr+itr,(idi-1)*nal+ial,id,ie,ia,iz,it) = 0.
							aW	((il-1)*ntr+itr,(idi-1)*nal+ial,id,ie,ia,iz,it) = aU((il-1)*ntr+itr,(idi-1)*nal+ial,id,ie,ia,iz,it)
						enddo
					endif
					Vtest1 = sum((V((il-1)*ntr+itr,(idi-1)*nal+ial,id,ie,:,iz,it) &
						& - V0((il-1)*ntr+itr,(idi-1)*nal+ial,id,ie,:,iz,it))**2)

					sumer = Vtest1	+ sumer
					maxer = max(Vtest1,maxer)
					verloc(:,iz,ie,id,ial) = (V((il-1)*ntr+itr,(idi-1)*nal+ial,id,ie,:,iz,it)-V0((il-1)*ntr+itr,(idi-1)*nal+ial,id,ie,:,iz,it))**2

					if(print_lev >3) then
						call veci2csv(ia_o,"ia_o_VW.csv", wo)
						if(wo == 0) wo = 1
					endif
	!				enddo !iz  enddo !ie   	enddo !id  	enddo !ial
				enddo
			!$OMP  END PARALLEL do

				if(maxval(verloc)>maxer_v) then
					maxer_i = maxloc(verloc)
					maxer_v = maxval(verloc)
				endif

				!update VW0, V0
				do ial=1,nal	!Loop over alpha (ai)
				do itr=1,ntr	!loop over trend
				do id=1,nd	!Loop over disability index
				do ie=1,ne	!Loop over earnings index
				do iz=1,nz	!Loop over TFP
				do ia =1,na
					VW0((il-1)*ntr+itr,(idi-1)*nal+ial,id,ie,ia,iz,it) = VW((il-1)*ntr+itr,(idi-1)*nal+ial,id,ie,ia,iz,it)
					V0 ((il-1)*ntr+itr,(idi-1)*nal+ial,id,ie,ia,iz,it) =  V((il-1)*ntr+itr,(idi-1)*nal+ial,id,ie,ia,iz,it)
				enddo !ia
				enddo !ial
				enddo !itr
				enddo !id
				enddo !ie
				enddo !iz


				if (print_lev >2) then
					wo = 0
					itr = tri0
					do ial=1,nal	!Loop over alpha (ai)
					do ie=1,ne	!Loop over earnings index
					do iz=1,nz	!Loop over TFP
						! matrix in disability and assets
						call mat2csv(VW((il-1)*ntr+itr,(idi-1)*nal+ial,:,ie,:,iz,it) ,"VW_it.csv",wo)
						call mati2csv(aW((il-1)*ntr+itr,(idi-1)*nal+ial,:,ie,:,iz,it) ,"aW_it.csv",wo)
						call mati2csv(gwork((il-1)*ntr+itr,(idi-1)*nal+ial,:,ie,:,iz,it) ,"gwork_it.csv",wo)
						call mat2csv(gwork_dif((il-1)*ntr+itr,(idi-1)*nal+ial,:,ie,:,iz,it) ,"gwork_idf_it.csv",wo)
						if(wo==0) wo =1
					enddo !iz
					enddo !ie
					enddo !ial
				endif




	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
				! End of iter iteration
				!------------------------------------------------!
				!Check |V-V0|<eps
				!------------------------------------------------!
				if(verbose > 3 .and. mod(iter,100).eq. 0) then
					write(*,*) sumer, iter, it, il
					write(*,*) maxer_v, maxer_i(1)
				endif
				maxer_hist(iter) = sumer
				if (sumer < Vtol .or. maxer_v < Vtol) then
					if(verbose >= 2) then
						write(*,*) sumer, iter, it, il
						write(*,*) maxer_v, maxer_i(1)
					endif
					exit !Converged
				elseif( sumer<1e-4 .and. iter .ge. 300 ) then
					if (sumer + 1e-4 .ge. sum(maxer_hist(iter-10:iter-1))/10._dp) then ! pretty close, but not converging anymore
						if(verbose >= 2) then
							write(*,*) sumer, iter, it, il
							write(*,*) maxer_v, maxer_i(1)
						endif
						exit
					endif
				endif

				iter=iter+1
				if(iter>=maxiter-1) then
					iter_timeout = iter_timeout+1
				endif
!				if(iter<200) smthV0param = smthV0param*1.01_dp !tighten up the discrete choice
			enddo	!iter: V-iter loop
	!WRITE(*,*) il, itr, idi, it
		enddo	!t loop, going backwards
		if(verbose>0 .and. iter_timeout>0) then
			print*, "VW did not converge on il,idi", il, " " , idi, " with error", sumer, " and max ", maxer
			print*, "maxer_v " , maxer_v
			print*, "maxer_i ", maxer_i
		endif
		enddo	!idi

		enddo	!il

		if(verbose>0 .and. iter_timeout>0) print*, "did not converge ", iter_timeout, " times"
		iter_timeout = 0
		! this plots work-rest and di application on the cross product of alphai and deltai and di
		if(print_lev >1 .and. caselabel == "") then
			itr = tri0
			il  = 1
			wo  = 0

			do id  = 1,nd
			do ie  = 1,ne
				call mati2csv(aD(id,ie,:,:),"aD"//trim(caselabel)//".csv",wo)
				call mat2csv (VD(id,ie,:,:),"VD"//trim(caselabel)//".csv",wo)

				call veci2csv(aR(2,id,ie,:),"aR"//trim(caselabel)//".csv",wo)
				call vec2csv (VR(2,id,ie,:),"VR"//trim(caselabel)//".csv",wo)
				if(wo == 0) wo =1
			enddo
			enddo

			wo = 0
			itr=tri0
			do idi=1,ndi	!loop over delta(idi)
			do ial=1,nal	!Loop over alpha (al)
			do ie=1,ne	!Loop over earnings index
			do iz=1,nz	!Loop over TFP
				do it = TT-1,1,-1
					! matrix in disability and assets
					call mat2csv(V((il-1)*ntr+itr,(idi-1)*nal+ial,:,ie,:,iz,it) ,"V"//trim(caselabel)//".csv",wo)

					call mat2csv(VW((il-1)*ntr+itr,(idi-1)*nal+ial,:,ie,:,iz,it) ,"VW"//trim(caselabel)//".csv",wo)
					call mati2csv(aW((il-1)*ntr+itr,(idi-1)*nal+ial,:,ie,:,iz,it) ,"aW"//trim(caselabel)//".csv",wo)

					call mat2csv(VU((il-1)*ntr+itr,(idi-1)*nal+ial,:,ie,:,iz,it) ,"VU"//trim(caselabel)//".csv",wo)
					call mati2csv(aU((il-1)*ntr+itr,(idi-1)*nal+ial,:,ie,:,iz,it) ,"aU"//trim(caselabel)//".csv",wo)

					call mat2csv(VN((il-1)*ntr+itr,(idi-1)*nal+ial,:,ie,:,iz,it) ,"VN"//trim(caselabel)//".csv",wo)
					call mati2csv(aN((il-1)*ntr+itr,(idi-1)*nal+ial,:,ie,:,iz,it) ,"aN"//trim(caselabel)//".csv",wo)

					call mati2csv(gwork((il-1)*ntr+itr,(idi-1)*nal+ial,:,ie,:,iz,it) ,"gwork"//trim(caselabel)//".csv",wo)
					call mat2csv(gwork_dif((il-1)*ntr+itr,(idi-1)*nal+ial,:,ie,:,iz,it) ,"gwork_dif"//trim(caselabel)//".csv",wo)

					call mati2csv(gapp((il-1)*ntr+itr,(idi-1)*nal+ial,:,ie,:,iz,it) ,"gapp"//trim(caselabel)//".csv",wo)
					call mat2csv(gapp_dif((il-1)*ntr+itr,(idi-1)*nal+ial,:,ie,:,iz,it) ,"gapp_dif"//trim(caselabel)//".csv",wo)

					if(wo==0) wo =1
				enddo !it
			enddo !iz
			enddo !ie
			enddo !ial
			enddo !idi

			call vec2csv(maxer_hist, "maxer_hist.csv")
		endif

		deallocate(maxer_hist)
		deallocate(verloc)
		deallocate(VR0,VD0,VN0,VU0,VW0,V0)
	!		deallocate(VR,VD,VN,VU,VW,V)
	!		deallocate(aR,aD,aN,aW,aU,gwork,gapp,gapp_dif,gwork_dif)

	end subroutine sol
end module sol_val


!**************************************************************************************************************!
!**************************************************************************************************************!
!**************************************************************************************************************!
!**************************************************************************************************************!
!						Simulate from solution					       !
!**************************************************************************************************************!
!**************************************************************************************************************!
!**************************************************************************************************************!
!**************************************************************************************************************!


module sim_hists
	use V0para
	use helper_funs
! module with subroutines to solve the model and simulate data from it

	implicit none

	contains


	subroutine draw_fndsepi(fndsep_i_draw, fndsep_i_int, fndarrive_draw, j_i, seed0, success)
	! draws depreciation rates and indices on the delta grid (i.e at the discrete values)
		implicit none

		integer, intent(in) :: seed0
		integer, intent(in), dimension(:) :: j_i
		integer, intent(out) :: success
		real(dp), dimension(:,:) :: fndsep_i_draw,fndarrive_draw
		integer, dimension(:,:,:) :: fndsep_i_int
		integer :: ss=1, Nsim, m,i,it
		real(dp) :: fndgrid_i
		integer, allocatable :: bdayseed(:)

		call random_seed(size = ss)
		allocate(bdayseed(ss))
		forall(m=1:ss) bdayseed(m) = (m-1)*100 + seed0
		call random_seed(put = bdayseed(1:ss) )

		Nsim = size(fndsep_i_draw,1)

		do i=1,Nsim
			call random_number(fndgrid_i) ! draw uniform on 0,1
			fndsep_i_draw(i,1) = fndgrid_i
			call random_number(fndgrid_i) ! draw uniform on 0,1
			fndsep_i_draw(i,2) = fndgrid_i
			do it=1,Tsim
				call random_number(fndgrid_i) ! draw uniform on 0,1
				fndarrive_draw(i,it) = fndgrid_i
			enddo
		enddo

		call set_fndsepi(fndsep_i_int,fndsep_i_draw,j_i)
		success = 0
		deallocate(bdayseed)
	end subroutine draw_fndsepi

	subroutine set_fndsepi( fndsep_i_int,fndsep_i_draw,j_i)

		real(dp), dimension(:,:), intent(in) :: fndsep_i_draw
		integer , dimension(:,:,:), intent(out) :: fndsep_i_int
		integer, intent(in), dimension(:) :: j_i
		integer :: ss=1, si_int,fi_int,m,i,ij,iz
		real(dp) :: fndgridL, fndgridH,fndwtH,fndwtL,fndgrid_i
		real(dp) :: sepgridL, sepgridH,sepwtH,sepwtL,sepgrid_i
		real(dp), dimension(nl+1,nj,nz) :: fndcumwt,sepcumwt

		fndwt	 = 1._dp/dble(nl) ! initialize with equal weight
		sepwt	 = 1._dp/dble(nl) ! initialize with equal weight
		fndcumwt = 0.
		sepcumwt = 0.

		do iz=1,nz
			do i=1,nl
			!		fndgrid(i,iz) = fndgrid(i,iz)*0.5_dp
			!		sepgrid(i,iz) = sepgrid(i,iz)*avg_unrt/(1.-avg_unrt)

			end do
		enddo

		do iz=1,nz
			fndgridL = 0.
			sepgridL = 0.
			do i=1,nl/2
				fndgridL = fndgrid(i,iz)/dble(nl/2) + fndgridL
				sepgridL = sepgrid(i,iz)/dble(nl/2) + sepgridL
			enddo
			fndgridH = 0.
			sepgridH = 0.
			do i= 1+nl/2,nl
				fndgridH = fndgrid(i,iz)/dble(nl-nl/2) + fndgridH
				sepgridH = sepgrid(i,iz)/dble(nl-nl/2) + sepgridH
			enddo

			do ij=1,nj
				! choose the mean to match the target mean by occupation
				fndwtH = (fndrate(iz,ij)-fndgridL)/(fndgridH-fndgridL)
				sepwtH = (seprisk(iz,ij)-sepgridL)/(sepgridH-sepgridL)
				fndwtL = 1._dp - fndwtH
				sepwtL = 1._dp - sepwtH
				do i=1,nl/2
					fndwt(i,ij,iz) = fndwtL/dble(nl/2)
					sepwt(i,ij,iz) = sepwtL/dble(nl/2)
				enddo
				do i=1+nl/2,nl
					fndwt(i,ij,iz) = fndwtH/dble(nl-nl/2)
					sepwt(i,ij,iz) = sepwtH/dble(nl-nl/2)
				enddo
			enddo
			!setup fndcumwt,sepcumwt
			do ij=1,nj
				do i=1,nl
					fndcumwt(i+1,ij,iz) = fndwt(i,ij,iz) + fndcumwt(i,ij,iz)
					sepcumwt(i+1,ij,iz) = sepwt(i,ij,iz) + sepcumwt(i,ij,iz)

				enddo
			enddo
		enddo !iz

		do iz=1,nz
		do i=1,Nsim
			fndgrid_i = fndsep_i_draw(i,1)
			sepgrid_i = fndsep_i_draw(i,2)
			fi_int = finder(fndcumwt(:,j_i(i),iz),fndgrid_i)
			si_int = finder(sepcumwt(:,j_i(i),iz),sepgrid_i)
			fndsep_i_int(i,1,iz) = fi_int
			fndsep_i_int(i,2,iz) = si_int
		enddo
		enddo

		if(print_lev >=2 .and. caselabel == "") &
		&	call mat2csv(fndwt(:,:,2),"fndwt"//trim(caselabel)//".csv")
		if(print_lev >=2 .and. caselabel == "") &
		&	call mat2csv(sepwt(:,:,2),"sepwt"//trim(caselabel)//".csv")

	end subroutine set_fndsepi


	subroutine draw_deli(del_i_draw, del_i_int, j_i, seed0, success)
	! draws depreciation rates and indices on the delta grid (i.e at the discrete values)
		implicit none

		integer, intent(in) :: seed0
		integer, intent(in), dimension(:) :: j_i
		integer, intent(out) :: success
		real(dp), dimension(:) :: del_i_draw
		integer, dimension(:) :: del_i_int
		integer :: ss=1, Ndraw, m,i
		real(dp) :: delgrid_i
		integer, allocatable :: bdayseed(:)

		call random_seed(size = ss)
		allocate(bdayseed(ss))
		forall(m=1:ss) bdayseed(m) = (m-1)*100 + seed0
		call random_seed(put = bdayseed(1:ss) )

		Ndraw = size(del_i_draw)

		do i=1,Ndraw
			call random_number(delgrid_i) ! draw uniform on 0,1
			del_i_draw(i) = delgrid_i
		enddo

		call set_deli(del_i_int,del_i_draw,j_i)

		success = 0
		deallocate(bdayseed)
	end subroutine draw_deli

	subroutine set_deli( del_i_int,del_i_draw,j_i)

		real(dp), dimension(:), intent(in) :: del_i_draw
		integer , dimension(:), intent(out) :: del_i_int
		integer, intent(in), dimension(:) :: j_i
		integer :: ss=1, di_int,m,i,ij,idi
		real(dp) :: delgridL, delgridH,delwtH,delwtL,delgrid_i


		delwt	 = 1._dp/dble(ndi) ! initialize with equal weight
		if(del_by_occ .eqv. .true.) then ! give weight according to mean delta by occ
			delgridL = 0. !average of high cells
			do i=1,ndi/2
				delgridL = delgrid(i)/dble(ndi/2) + delgridL
			enddo
			delgridH = 0. !average of low cells
			do i= 1+ndi/2,ndi
				delgridH = delgrid(i)/dble(ndi-ndi/2) + delgridH
			enddo
			do ij=1,nj
				! choose the mean to match the target mean by occupation
				delwtH = (occdel(ij)-delgridL)/(delgridH-delgridL)
				delwtL = 1._dp - delwtH
				do i=1,ndi/2
					delwt(i,ij) = delwtL/dble(ndi/2)
				enddo
				do i=1+ndi/2,ndi
					delwt(i,ij) = delwtH/dble(ndi-ndi/2)
				enddo
			enddo
		else
			delgridL = delgrid(1)
			delgridH = delgrid(ndi)
			do ij=1,nj
				delwt(:,ij) = 0.
				if(delgridH - delgridL > 1e-4) then
					delwtH = (0. - delgridL)/(delgridH-delgridL)
					delwtL = 1._dp - delwtH
					do i=1,ndi/2
						delwt(i,:) = delwtL/dble(ndi/2)
					enddo
					do i=1+ndi/2,ndi
						delwt(i,:) = delwtH/dble(ndi-ndi/2)
					enddo
				else
					delwt(1,:)=1.
				endif
			enddo
		endif
		!setup delcumwt
		delcumwt = 0.
		do ij=1,nj
			do idi=1,ndi
				delcumwt(idi+1,ij) = delwt(idi,ij) + delcumwt(idi,ij)
			enddo
		enddo

		delgridL = minval(delgrid)
		delgridH = maxval(delgrid)
		do i=1,Nsim
			delgrid_i = del_i_draw(i)
			di_int = finder(delcumwt(:,j_i(i)),delgrid_i)
			di_int = max(1, min(di_int,ndi))
			del_i_int(i) = di_int
		enddo

		if(print_lev >=2) &
		&	call mat2csv(delwt,"delwt"//trim(caselabel)//".csv")

	end subroutine set_deli

	subroutine draw_status_innov(status_it_innov, dead_draw,seed0, success)
	! draws innovations to d, will be used if working and relevant
		implicit none

		integer, intent(in) :: seed0
		integer, intent(out) :: success
		real(dp), dimension(:,:) :: status_it_innov, dead_draw
		integer :: ss=1, m,i,it
		real(dp) :: s_innov
		integer, allocatable :: bdayseed(:)

		call random_seed(size = ss)
		allocate(bdayseed(ss))
		forall(m=1:ss) bdayseed(m) = (m-1)*100 + seed0
		call random_seed(put = bdayseed(1:ss) )

		do i=1,Nsim
			do it=1,Tsim
				call rand_num_closed(s_innov)
				status_it_innov(i,it) = s_innov
				call random_number(s_innov)
				dead_draw(i,it) = s_innov
			enddo
		enddo
		success = 0
		deallocate(bdayseed)
	end subroutine draw_status_innov

	subroutine draw_alit(al_it,al_int_it, al_it_innov, d_it, seed0, success)
		implicit none

		integer, intent(in) :: seed0
		integer, intent(out),optional :: success
		integer, intent(in), dimension(:,:) :: d_it
		real(dp), dimension(:,:) :: al_it
		real(dp), dimension(:,:) :: al_it_innov
		integer, dimension(:,:) :: al_int_it
		integer :: ss=1, Ndraw, t,m,i
		real(dp) :: alf_innov
		integer, allocatable :: bdayseed(:)

		call random_seed(size = ss)
		allocate(bdayseed(ss))
		forall(m=1:ss) bdayseed(m) = (m-1)*100 + seed0
		call random_seed(put = bdayseed(1:ss) )

		Ndraw = size(al_it,1)
		do i=1,Ndraw
			do t=1,Tsim
				if( al_contin .eqv. .true. ) then
					call random_normal(alf_innov)
					al_it_innov(i,t) = alf_innov
				else
					call rand_num_closed(alf_innov)
					al_it_innov(i,t) = alf_innov
				endif
			enddo
		enddo

		call set_alit(al_it,al_int_it,al_it_innov,d_it,success)

		deallocate(bdayseed)
	endsubroutine draw_alit

	subroutine set_alit(al_it,al_int_it,al_it_innov, d_it,  success)
	! draws alpha shocks and idices on the alpha grid (i.e at the discrete values)
		implicit none

		integer, intent(out),optional :: success
		integer, intent(in), dimension(:,:) :: d_it
		real(dp), dimension(:,:) :: al_it
		real(dp), dimension(:,:) :: al_it_innov
		integer, dimension(:,:) :: al_int_it
		integer :: ss=1, Ndraw, alfgrid_int, t,m,i,k,id
		real(dp) :: alfgridL, alfgridH,alf_innov,alfgrid_i,alf_i
		real(dp) :: alfsig(nd),alfrhot(nd),alfsigt(nd),alfcondsigt(nd)
		real(dp) :: alfgrid_minE(nd),alfgrid_maxE(nd),alfgrid_Uval !min max value while employed and val of unemp
		integer, allocatable :: bdayseed(:)
		real(dp), allocatable :: cumpi_al(:,:,:)

		allocate(cumpi_al(nal,nal+1,nd))

		do id=1,nd
			alfgrid_minE(id) = alfgrid(2,id)
			alfgrid_maxE(id) = alfgrid(nal,id)
		enddo
		alfgrid_Uval = alfgrid(1,1)

		cumpi_al =0.
		Ndraw = size(al_it,1)

		do id=1,nd
			alfrhot(id) = alfrho(id)**(1./tlen)
			alfsig(id) = (alfcondsig(id)**2/(1-alfrho(id)**2))**0.5
			alfsigt(id) = (alfsig(id)**2/tlen)**0.5
			alfcondsigt(id) = (alfsigt(id)**2*(1-alfrhot(id)**2))**0.5

			do i=1,nal
				do k=2,nal+1
					cumpi_al(i,k,id) = pialf(i,k-1,id)+cumpi_al(i,k-1,id)
				enddo
			enddo
		enddo

		success =0

		al_it = 0._dp
		al_int_it = nal/2+1
		do i=1,Ndraw

			! draw starting values
			id = 1
			!call random_normal(alf_innov) ! draw normal disturbances on 0,1
			alf_innov = al_it_innov(i,1)
			! transform it by the ergodic distribution for the first period:
			alf_i = alf_innov*alfsigt(id) + alfmu(id)

			if((alf_i >alfgrid_maxE(id)) .or. (alf_i < alfgrid_minE(id)) ) success = 1+success !count how often we truncate
			!impose bounds
			alf_i = max(alf_i,alfgrid_minE(id))
			alf_i = min(alf_i,alfgrid_maxE(id))
			alfgrid_int = finder(alfgrid(:,id),alf_i)
			alfgrid_int = max(2, min(alfgrid_int,nal) )


			! draw sequence (initialize with Tsim values):
			do t=(-2*Tsim),Tsim
				if(t>=1) then
					id = d_it(i,t)
					if(id<1) id =1 !if they're not born, id is initialized to 0. This causes havok.
					alf_innov = al_it_innov(i,t)
				elseif(t<=0) then
					id = 1
					alf_innov = al_it_innov(i,abs(mod(t-1,Tsim-1))+1)
				endif
				if(al_contin .eqv. .true.) then
						alf_i  = alfrhot(id)*alf_i + alfcondsigt(id)*alf_innov + alfmu(id)*(1-alfrhot(id))
						alf_i = max(alf_i,alfgrid_minE(id))
						alf_i = min(alf_i,alfgrid_maxE(id))
						if(t >= 1)  then
							al_it(i,t) = alf_i  ! log of wage shock
							if(alf_i >alfgrid_maxE(id) .or. alf_i < alfgrid_minE(id)) success = 1+success !count how often we truncate
						endif
						alfgrid_int = min(finder(alfgrid(:,id),alf_i),nal-1)
						if( (alf_i - alfgrid(alfgrid_int,id))/(alfgrid(alfgrid_int+1,id)- alfgrid(alfgrid_int,id)) >0.5 ) alfgrid_int = alfgrid_int + 1
						alfgrid_int = max(min(alfgrid_int,nal),2)
				else
					alfgrid_int = finder(cumpi_al(alfgrid_int,:,id), alf_innov )
					alfgrid_int = max(min(alfgrid_int,nal),1)
					if(t>=1) al_it(i,t) = alfgrid(alfgrid_int,id) ! log of wage shock, on grid
				endif
				if(t>=1) al_int_it(i,t) = alfgrid_int
			enddo
		enddo
		if(success > 0.2*Ndraw*Tsim)  success = success
		if(success <= 0.2*Ndraw*Tsim) success = 0

!		call mat2csv(al_it,"drawn_al_it.csv")

		deallocate(cumpi_al)

	end subroutine set_alit

	subroutine draw_dit( d_it,health_it_innov, del_i, age_it, seed0,success )
		integer, allocatable :: bdayseed(:)
		real(dp) :: health_it_innov(:,:)
		integer  :: d_it(:,:)
		integer  :: del_i(:),age_it(:,:)
		integer :: seed0
		integer, optional :: success
		integer :: ss,i,it,m
		real(dp) :: health_draw

		call random_seed(size = ss)
		allocate(bdayseed(ss))
		forall(m=1:ss) bdayseed(m) = (m-1)*100 + seed0
		call random_seed(put = bdayseed(1:ss) )

		do i=1,Nsim
			do it=1,Tsim
				call random_number(health_draw)
				health_it_innov(i,it) = health_draw
			enddo
		enddo

		call set_dit(d_it,health_it_innov,del_i,age_it,success)

		deallocate(bdayseed)

	end subroutine draw_dit

	subroutine set_dit(d_it, health_it_innov, del_i, age_it,success)

		real(dp) :: health_it_innov(:,:)
		integer  :: d_it(:,:)
		integer  :: del_i(:),age_it(:,:)
		integer, optional :: success

		real(dp) :: cumPrDage(nd+1,TT), cumPrDageDel(nd+1,TT,ndi),cumpid(nd,nd+1,ndi,TT-1)
		integer  :: it,i,id,idi,age_hr,del_hr,d_hr


		cumpid = 0.
		cumPrDage = 0.
		cumPrDageDel = 0.
		do idi=1,ndi
		do it =1,TT-1
			do id =1,nd
				do i =1,nd
					cumpid(id,i+1,idi,it) = pid(id,i,idi,it)+cumpid(id,i,idi,it)
				enddo
			enddo
		enddo
		enddo

		do it =1,TT
			do id =1,nd
				cumPrDage(id+1,it) = PrDage(id,it) +cumPrDage(id,it)
				do idi=1,ndi
					cumPrDageDel(id+1,it,idi) = PrDageDel(id,it,idi) + cumPrDageDel(id,it,idi)
				enddo
			enddo
		enddo

		it=1
		do i=1,Nsim
			age_hr = age_it(i,it)
			del_hr = del_i(i)
			if(age_hr>0) then
				d_hr = locate(cumPrDageDel(:,age_hr,del_hr),health_it_innov(i,it) )
				d_it(i,it) = d_hr
			else
				d_it(i,it) = 1
			endif
		enddo

		do i=1,Nsim
			del_hr = del_i(i)
			do it=1,(Tsim-1)
				d_hr = d_it(i,it)
				age_hr = age_it(i,it)
				if(age_hr>0) then
					if(it>1) then
						if( age_it(i,it)>0 .and. (age_it(i,it-1)<= 0)  ) then !born
							d_it(i,it) = locate(cumPrDageDel(:,age_hr,del_hr),health_it_innov(i,it) )
						endif
					endif
					!move forward d
					if(age_hr .lt. TT)  then
						d_it(i,it+1) = locate(   cumpid(d_hr,:,del_hr,age_hr) ,health_it_innov(i,it))
					else
						d_it(i,it+1) = d_hr
					endif
				else
					d_it(i,it+1) = d_hr
				endif
			enddo
		enddo

	end subroutine set_dit

	subroutine draw_ji(j_i,jshock_ij,born_it,seed0, success)
		implicit none
		integer	:: j_i(:),born_it(:,:)
		real(dp) :: jshock_ij(:,:)
		real(dp) :: jwt
		integer	:: i,m,ss=1,ij,it
		integer,intent(in)  :: seed0
		integer,intent(out) :: success
		integer, allocatable :: bdayseed(:)
		real(dp) :: draw_i

		j_i = 0
		if(nj>1) then
			call random_seed(size = ss)
			allocate(bdayseed(ss))
			forall(m=1:ss) bdayseed(m) = (m-1)*100 + seed0
			call random_seed(put = bdayseed(1:ss) )
			if(j_rand .eqv. .false. ) then
				!draw gumbel-distributed shock
				do i = 1,Nsim
					do ij = 1,nj
						call random_gumbel(draw_i)
						jshock_ij(i,ij) = draw_i
					enddo
				enddo
			else
				jshock_ij = 0._dp
				do i = 1,Nsim
					call random_number(draw_i)
					jshock_ij(i,:) = draw_i
				enddo
				call set_ji(j_i,jshock_ij,born_it)
			endif

			deallocate(bdayseed)
		else
			jshock_ij = 0.
			j_i = 1
		endif

		success = 0

	end subroutine draw_ji

	subroutine set_ji(j_i,jshock_ij,born_it)
		implicit none
		integer	:: j_i(:),born_it(:,:)
		real(dp) :: jshock_ij(:,:)
		real(dp) :: jwt
		integer	:: i,m,ss=1,ij,it
		integer, allocatable :: bdayseed(:)
		real(dp) :: Njcumdist(nj+1,Tsim),Njconstcumdist(nj+1)
		real(dp) :: draw_i

		Njcumdist = 0._dp
		Njconstcumdist = 0._dp
		j_i = 0
		it=1
		do ij=1,nj
			Njconstcumdist(ij+1) = occpr_trend(it,ij) + Njconstcumdist(ij)
		enddo
		if(nj>1) then
			!do periods 2:Tsim
			do it=2,Tsim
				if( occ_dat .eqv. .true. ) then
					do ij=1,nj
						Njcumdist(ij+1,it) = occpr_trend(it,ij) + Njcumdist(ij,it)
					enddo
				else
					Njcumdist(:,it) = Njconstcumdist
				endif
				do i=1,Nsim
					if( born_it(i,it)==1) then
						draw_i = jshock_ij(i,1)
						j_i(i) = finder(Njcumdist(:,it),draw_i)
						if(j_i(i) < 1 ) j_i(i) = 1
						if(j_i(i) > nj) j_i(i) = nj
					endif
				enddo
			enddo
			!do period 1
			it=1
			do ij=1,nj
				Njcumdist(ij+1,it) = occpr_trend(it,ij) + Njcumdist(ij,it)
			enddo
			if( occ_dat .eqv. .true. ) then
				do ij=1,nj
					Njcumdist(ij+1,it) = occpr_trend(it,ij) + Njcumdist(ij,it)
				enddo
			else
				Njcumdist(:,it) = Njconstcumdist
			endif
			do i=1,Nsim
				if( j_i(i) ==0) then !not born in 2:Tsim and so doesn't have an occupation yet
					draw_i = jshock_ij(i,1)
					j_i(i) = finder(Njcumdist(:,it),draw_i)
					if(j_i(i) < 1 ) j_i(i) = 1
					if(j_i(i) > nj) j_i(i) = nj
				endif
			enddo

			!make sure everyone has a job:
			do i=1,Nsim
				if( j_i(i)==0 ) then
					draw_i = jshock_ij(i,1)
					j_i(i) = nint( draw_i*(nj-1)+1 )
				endif
			enddo

		else
			jshock_ij = 0.
			j_i = 1
		endif

	end subroutine set_ji

	subroutine draw_zjt(z_jt_select, z_jt_innov, seed0, success)
		implicit none

		real(dp), intent(out) :: z_jt_innov(:) !this will drive the continuous AR process
		real(dp), intent(out) :: z_jt_select(:) !this will select the state if we use a markov chain
		integer	:: it=1,i=1,ij=1,iz=1,izp=1,m=1,ss=1
		integer,intent(in) :: seed0
		integer,intent(out) :: success
		integer, allocatable :: bdayseed(:)
		integer :: NBER_start_stop(5,2)
		real(dp) :: z_innov=1.
		integer  :: nzblock = nz/2

		nzblock = nz

		call random_seed(size = ss)
		allocate(bdayseed(ss))
		forall(m=1:ss) bdayseed(m) = (m-1)*100 + seed0
		call random_seed(put = bdayseed(1:ss) )

		NBER_start_stop = 0
		!compute NBER dates in case of NBER_tseq == 1
		!NO LONGER: 1980 + 0/4 -> 1980 + 2/4
		!NO LONGER: 1981 + 2/4 -> 1982 + 3/4

		! 1990 + 2/4 -> 1991 + 0/4
		! 2001 + 2/4 -> 2001 + 3/4
		! 2007 + 3/4 -> 2009 + 1/4
		NBER_start_stop(1,1) =  7*itlen + 2*( itlen/4 ) +1
		NBER_start_stop(1,2) =  8*itlen + 0*( itlen/4 ) +1
		NBER_start_stop(2,1) = 18*itlen + 2*( itlen/4 ) +1
		NBER_start_stop(2,2) = 18*itlen + 3*( itlen/4 ) +1
		NBER_start_stop(3,1) = 24*itlen + 3*( itlen/4 ) +1
		NBER_start_stop(3,2) = 26*itlen + 1*( itlen/4 ) +1

		NBER_start_stop = NBER_start_stop+itlen*TossYears !shift the whole thing

		!start cycles on 1st

		ss = 1
		do it = 1,Tsim
			if(NBER_tseq .eqv. .true. ) then
				if((it >= NBER_start_stop(ss,1)) .and. (it <= NBER_start_stop(ss,2)) ) then
					z_jt_innov(it) = -1.
					z_jt_select(it) = 0.
				else
					z_jt_innov(it)  = 1.
					z_jt_select(it) = 1.
				endif
				if(it .eq. NBER_start_stop(ss,2) .and. ss < 5) &
					& ss = ss+1
			else
				call random_normal(z_innov)
				z_jt_innov(it) = z_innov
				z_jt_select(it) = alnorm(z_innov,.false.)
			endif
		enddo
		success = 0
		deallocate(bdayseed)
		!call mat2csv(cumpi_z,"cumpi_z.csv")
	end subroutine draw_zjt

	subroutine draw_age_it(age_it, born_it, age_draw,seed0, success)

		integer,intent(out) :: age_it(:,:),born_it(:,:)
		real(dp), intent(out) :: age_draw(:,:)
		integer	:: i,m, Nm,ss=1
		integer,intent(in) :: seed0
		integer,intent(out) :: success
		integer, allocatable :: bdayseed(:)
		real(dp) :: rand_age,rand_born

		call random_seed(size = ss)
		allocate(bdayseed(ss))
		forall(m=1:ss) bdayseed(m) = (m-1)*100 + seed0
		call random_seed(put = bdayseed(1:ss) )

		Nm = size(age_draw,2)

		do i =1,Nsim
			do m = 1,Nm
				call random_number(rand_born)
				age_draw(i,m) = rand_born
			enddo! m=1,50.... if never gets born
		enddo! i=1,Nsim

		call set_age(age_it,born_it,age_draw)

		deallocate(bdayseed)
		success = 0
	end subroutine draw_age_it

	subroutine set_age(age_it, born_it, age_draw)

		integer,intent(out) :: age_it(:,:),born_it(:,:)
		real(dp), intent(in) :: age_draw(:,:)
		integer	:: it,it0,itp,i,m,m1,m2,m3,Nm,bn_i, age_ct, age_jump
		real(dp), dimension(TT) :: cumpi_t0
		real(dp), dimension(TT-1) :: prob_age_nTT
		real(dp) :: rand_age,rand_born
		real(dp), dimension(Tsim) :: hazborn_hr_t,prborn_hr_t,cumprnborn_t


		!set up cumulative probabilities for t0 and conditional draws
		!not begining with anyone from TT
		prob_age_nTT = prob_age(:,1)
		!also include TT
		!prob_age_nTT = prob_age
		if(demog_dat .eqv. .true.) then
			prborn_hr_t  = prborn_t
			hazborn_hr_t = hazborn_t
		else !use a flat age profile
			prborn_hr_t = prborn_constpop
			hazborn_hr_t = hazborn_constpop
		endif
		cumprnborn_t(1) = prborn_t(1)
		do it=2,Tsim
			cumprnborn_t(it) = prborn_hr_t(it) + cumprnborn_t(it-1)
		enddo
		cumprnborn_t = cumprnborn_t/cumprnborn_t(Tsim)

		cumpi_t0 = 0.
		age_it = 0
		do it=1,TT-1
			cumpi_t0(it+1) = prob_age_nTT(it) + cumpi_t0(it)
		enddo

		Nm = size(age_draw,2)
		born_it = 0
		do i =1,Nsim
			!pick the birthdate
			!it = finder(cumprnborn_t, dble(i-1)/dble(Nsim-1))
			it = finder(cumprnborn_t, age_draw(i,1))
			if(it>1) then
				age_it(i,it) =1
				born_it(i,it) = 1
				age_ct = 0
				bn_i = 0 ! initialize, not yet born
			else
				rand_age = age_draw(i,2 )
				age_it(i,it) = finder(cumpi_t0,rand_age)
				born_it(i,it) = 0
				if(age_it(i,it)==1) then
					age_ct = nint(age_draw(i,3)*youngD*tlen)
				else
					age_ct = nint(age_draw(i,3)*oldD*tlen)
				endif
				bn_i = 1 ! initialize, born
			endif
			it0 = max(it,2)
			!count up periods:
			do it=it0,Tsim
				if(age_it(i,it-1)<TT) then
					if( (born_it(i,it)== 1 ).and. ( bn_i == 0) ) then
						age_it(i,it) =1
						born_it(i,it) = 1
						bn_i = 1
						age_ct = 0
					elseif(bn_i == 1) then
						born_it(i,it) = 0
						age_ct = age_ct+1
						if( (age_draw(i,it+3) <1.- ptau( age_it(i,it-1) ) ) .and. (age_it(i,it-1) < TT) ) then
							age_it(i,it) = age_it(i,it-1)+1
							age_ct = 0
						else
							age_it(i,it) = age_it(i,it-1)
						endif
					else
						bn_i = 0
						born_it(i,it) = 0
						age_it(i,it) = 0
					endif
				else
					age_it(i,it) = TT
					born_it(i,it) = 0
				endif
			enddo

		enddo! i=1,Nsim

	end subroutine set_age

	subroutine draw_draw(drawi_ititer, drawt_ititer, age_it,al_int_it,del_i_int, seed0, success)

		integer,intent(in) :: seed0
		integer,intent(out) :: success
		integer,intent(in) :: age_it(:,:),al_int_it(:,:),del_i_int(:)
		integer, allocatable :: bdayseed(:)
		integer,intent(out) :: drawi_ititer(:,:),drawt_ititer(:,:)
		integer :: i,it,id,m,ss=1,drawt,drawi,ndraw,Ncols, seedi,brn_yr(Nsim),iter
		real(dp) :: junk

		call random_seed(size = ss)
		allocate(bdayseed(ss))
		Ndraw = size(drawi_ititer,1)
		Ncols = size(drawi_ititer,2)

		do i =1,Nsim
			it =1
			if(age_it(i,it)>0) then
				brn_yr(i) = 1
			else
				do it=2,Tsim
					if(age_it(i,it)>0) then
						brn_yr(i) = it
						exit
					endif
				enddo
			endif
		enddo

		do m=1,ss
			bdayseed(m) = (m-1)*100 + seed0
		enddo
		call random_seed(put = bdayseed )
		!need to draw these from age- and alpha- specific distributions for iterations > 1
		! \|/ this makes it much slower
		! OMP parallel do private(id, i, it, junk, seedi, drawi, drawt,ss,bdayseed,m)
		do id = 1,Ncols
			do i=1,Ndraw
			!initialize:
				it=brn_yr(i)
				iter = 0
				drawi_ititer(i,id) = i
				drawt_ititer(i,id) = it

				ageloop: do
					call random_number(junk)
					drawi = max(1,idnint(junk*Nsim))
					call random_number(junk)
					drawt = max(1,idnint(junk*Tsim))
					iter = iter+1
					if( (age_it(drawi,drawt) .eq. age_it(i,it)) .and. &
						(del_i_int(drawi)    .eq. del_i_int(i) ) ) then
						drawi_ititer(i,id) = drawi
						drawt_ititer(i,id) = drawt
						exit ageloop
					endif
				end do ageloop
			enddo
		enddo
		! OMP end parallel do

		deallocate(bdayseed)
		 success = 0
	end subroutine draw_draw

	subroutine draw_shocks(shk)

		implicit none
		type(shocks_struct) :: shk
		integer :: seed0,seed1, status
		integer :: time0,timeT,timert

		seed0 = 941987
		seed1 = 12281951

		if( readshocks .eqv. .false. ) then
			call system_clock(count_rate = timert)

			call system_clock(count = time0)
			if(verbose >2) print *, "Drawing types and shocks"
			call draw_age_it(shk%age_hist,shk%born_hist,shk%age_draw,seed0,status)
			seed0 = seed0 + 1
			call draw_ji(shk%j_i,shk%jshock_ij,shk%born_hist,seed1, status)
			seed1 = seed1 + 1
			call draw_deli(shk%del_i_draw, shk%del_i_int, shk%j_i, seed1, status)
			seed0 = seed0 + 1
			call draw_dit(shk%d_hist,shk%health_it_innov,shk%del_i_int,shk%age_hist,seed0,status)
			seed1 = seed1 + 1
			call draw_alit(shk%al_hist,shk%al_int_hist,shk%al_it_innov,shk%d_hist, seed0, status)
			seed0 = seed0 + 1
			call draw_fndsepi(shk%fndsep_i_draw, shk%fndsep_i_int, shk%fndarrive_draw, shk%j_i, seed0, status)
			seed0 = seed0 + 1
			call draw_zjt(shk%z_jt_select,shk%z_jt_innov, seed1, status)
			seed1 = seed1 + 1
			call draw_draw(shk%drawi_ititer, shk%drawt_ititer, shk%age_hist, shk%al_int_hist, shk%del_i_int, seed0, status)
			seed0 = seed0 + 1
			call draw_status_innov(shk%status_it_innov, shk%dead_draw,seed1, status)
			seed1 = seed1 + 1

			shk%drawn = 1
			call system_clock(count = timeT)
			print *, "draws took: ", dble((timeT-time0))/dble(timert)
		else
			call csv2mat("age_draw_hist.csv", shk%age_draw)
			call set_age(shk%age_hist,shk%born_hist,shk%age_draw)
			call csv2mat("al_it_hist.csv",shk%al_hist)
			call csv2mati("al_int_it_hist.csv",shk%al_int_hist)
			call csv2mat("jshock_ij_hist.csv",shk%jshock_ij)
			call set_ji(shk%j_i,shk%jshock_ij,shk%born_hist)
			call csv2vec("del_i_draw_hist.csv",shk%del_i_draw)
			call set_deli(shk%del_i_int,shk%del_i_draw,shk%j_i)
			call csv2mat("fndsep_i_draw.csv",shk%fndsep_i_draw)
			call set_fndsepi(shk%fndsep_i_int,shk%fndsep_i_draw,shk%j_i)
			call csv2vec("z_jt_innov_hist.csv",shk%z_jt_innov)
			call csv2vec("z_jt_select.csv",shk%z_jt_select)
			call csv2mati("drawi_hist.csv",shk%drawi_ititer)
			call csv2mati("drawt_hist.csv",shk%drawt_ititer)
			call csv2mat("status_it_innov.csv",shk%status_it_innov)
			call csv2mat("health_it_innov.csv",shk%health_it_innov)
			call csv2mat("dead_draw.csv",shk%dead_draw)
		endif
		! check the distributions
		if(print_lev > 1 .and. caselabel == "") then
			call mat2csv(shk%jshock_ij,"jshock_ij_hist.csv")
			call mat2csv(shk%age_draw,"age_draw_hist.csv")
			call vec2csv(shk%del_i_draw,"del_i_draw_hist.csv")
			call veci2csv(shk%del_i_int,"del_i_int_hist.csv")
			call veci2csv(shk%j_i,"j_i_hist.csv")
			call mat2csv(shk%al_hist,"al_it_hist.csv")
			call mati2csv(shk%al_int_hist,"al_int_it_hist.csv")
			! call vec2csv(shk%z_jt_innov,"z_jt_innov_hist.csv")
			! call vec2csv(shk%z_jt_innov,"z_jt_select.csv")
			call mati2csv(shk%age_hist,"age_it_hist.csv")
			call mat2csv(shk%fndsep_i_draw,"fndsep_i_draw.csv")
			! call mati2csv(shk%drawi_ititer,"drawi_hist.csv")
			! call mati2csv(shk%drawt_ititer,"drawt_hist.csv")
			call mat2csv(shk%status_it_innov,"status_it_innov.csv")
			call mat2csv(shk%health_it_innov,"health_it_innov.csv")
			call mat2csv(shk%dead_draw,"dead_draw.csv")
		endif

	end subroutine draw_shocks

	subroutine set_zjt(z_jt_macroint, z_jt_panel, shk)

		type(shocks_struct), intent(in) :: shk
		real(8) :: z_jt_panel(:,:)
		integer :: z_jt_macroint(:)
		! for selecting the state of zj
		real(dp) :: cumpi_z(nz,nz+1)
		real(dp) ::  cumergpi(nz+1)

		integer :: it,ij,iz,izp, zi_jt_t , nzblock
		real(8) :: muhere=0._dp, Zz_t=0._dp

		nzblock = nz

		call settfp() ! sets values for piz


		cumpi_z = 0._dp
		cumergpi = 0._dp

		do iz=1,nz
			do izp=1,nz
				cumpi_z(iz,izp+1) = piz(iz,izp) + cumpi_z(iz,izp)
			enddo
		enddo

		it = 1
		zi_jt_t = finder( cumergpi,shk%z_jt_select(it) )
		z_jt_macroint(it) = zi_jt_t
		if( zj_contin .eqv. .true.) then
			Zz_t = zsig*shk%z_jt_innov(it)
			forall (ij = 1:nj)	z_jt_panel(ij,it) = Zz_t*zscale(ij)
		else
			forall (ij = 1:nj) z_jt_panel(ij,it) = zgrid(zi_jt_t,ij)
		endif

		do it= 2,Tsim
			!z_jt_t = finder(cumpi_z(z_jt_t,:),z_innov ) <- random time-block transitions:

			! use conditional probability w/in time block
			zi_jt_t = finder(cumpi_z(zi_jt_t,:),shk%z_jt_select(it) )
			z_jt_macroint(it) = zi_jt_t


			if( zj_contin .eqv. .true.) Zz_t = (1.-zrho) * zmu + zsig*shk%z_jt_innov(it) + zrho*Zz_t
			do ij=1,nj
				muhere = 0.
				if( zj_contin .eqv. .true.) then
					!z_jt_panel(ij,it) = (1.-zrho)*muhere + zsig*shk%z_jt_innov(it) + zrho*z_jt_panel(ij,it-1)
					z_jt_panel(ij,it) = Zz_t*zscale(ij) + (1.-zrho)*muhere
				else
					z_jt_panel(ij,it) = zgrid(z_jt_macroint(it),ij)
				endif
			enddo
		enddo

	end subroutine set_zjt


	subroutine sim(vfs, pfs, hst, shk,occaggs)

		implicit none

		type(val_struct), intent(in), target :: vfs
		type(pol_struct), intent(in), target :: pfs
		type(hist_struct), intent(inout), target :: hst
		type(shocks_struct), intent(inout), target :: shk

		logical, optional :: occaggs

		integer :: i=1, ii=1, iter=1, it=1, it_old=1, ij=1,il=1, idi=1, id=1, Tret=1,wo=1, &
			&  seed0=1, seed1=1, status=1, m=1,ss=1, iter_draws=2,Ncol=100,nomatch=0, &
			&  ndets=5, nregobs


		integer, allocatable :: work_it(:,:), app_it(:,:) !choose work or not, apply or not
		integer, allocatable :: a_it_int(:,:),e_it_int(:,:),invol_it(:,:)
		!integer, allocatable :: hlthvocSSDI(:,:) ! got into ssdi on health or vocational considerations, 0= no ssdi, 1=health, 2=vocation
		real(dp), allocatable :: e_it(:,:)
		integer, allocatable :: brn_drawi_drawt(:,:,:)

		! FOR DEBUGGING
		real(dp), allocatable :: ewt_it(:,:)
        real(dp) ,allocatable :: trindif_it(:,:)

		! because the actual alpha is endgoenous to unemployment and trend

		real(dp), allocatable :: al_it_endog(:,:)
		integer, allocatable  :: al_int_it_endog(:,:)
		real(dp), allocatable :: trX_it(:,:),Xdets(:,:),xidets(:),xjdets(:), &
								& acoef(:),ecoef(:),adep(:),edep(:),&
								& cov_coef(:,:), tr_hist_vec(:)

		! write to hst, get from shk
		real(dp), pointer    :: work_dif_it(:,:), app_dif_it(:,:) !choose work or not, apply or not -- latent value
		real(dp), pointer	 :: di_prob_it(:,:)
		integer, pointer     :: born_it(:,:) ! born status, drawn randomly
		integer, pointer     :: del_i_int(:)  ! integer valued shocks
		integer, pointer     :: fndsep_i_int(:,:,:)  ! integer valued shocks
		integer, pointer     :: status_it(:,:)	!track W,U,N,D,R : 1,2,3,4,5
		integer, pointer     :: age_it(:,:)	! ages, drawn randomly
		integer, pointer     :: j_i(:)		! occupation, maybe random
		integer, pointer     :: d_it(:,:) 	! disability status
		integer, pointer     :: z_jt_macroint(:)! shocks to be drawn
		real(8), pointer     :: z_jt_panel(:,:)! shocks to be drawn

		integer, pointer     :: al_int_it(:,:)! integer valued shocks
		real(dp), pointer    :: occgrow_jt(:,:), occshrink_jt(:,:), occsize_jt(:,:)
		real(dp), pointer    :: a_it(:,:) 	! assets
		real(dp), pointer    ::	al_it(:,:)	! individual shocks
		real(dp), pointer    :: status_it_innov(:,:),health_it_innov(:,:),fndarrive_draw(:,:),dead_draw(:,:)
		real(dp), pointer    :: jshock_ij(:,:)
		integer,  pointer    :: drawi_ititer(:,:)
		integer,  pointer    :: drawt_ititer(:,:)

		real(dp), pointer	 :: z_jt_innov(:)
		real(dp), pointer	 :: z_jt_select(:)


		! read from vals
		real(dp), pointer ::	V(:,:,:,:,:,:,:),V_CF(:,:,:,:,:,:,:)	!Participant

		! read from pols
		real(dp), pointer ::	gapp_dif(:,:,:,:,:,:,:), gwork_dif(:,:,:,:,:,:,:) ! latent value of work/apply

		integer, pointer ::	aR(:,:,:,:), aD(:,:,:,:), aU(:,:,:,:,:,:,:), &
					aN(:,:,:,:,:,:,:), aW(:,:,:,:,:,:,:)
		integer, pointer ::	gapp(:,:,:,:,:,:,:), &
					gwork(:,:,:,:,:,:,:)

		logical :: ptrsuccess=.false., dead = .false.
		real(dp) :: pialf_conddist(nal), cumptau(TT+1),a_mean(TT-1),d_mean(TT-1),a_var(TT-1), &
				& d_var(TT-1),a_mean_liter(TT-1),d_mean_liter(TT-1),a_var_liter(TT-1),d_var_liter(TT-1), &
				& s_mean(TT-1),s_mean_liter(TT-1), occgrow_hr(nj),occsize_hr(nj),occshrink_hr(nj),PrAl1(nz),PrAl1St3(nz),totpopz(nz),&
				& simiter_dist(maxiter),simiter_status_dist(maxiter), a_ap_dist(Na,Na)

		! Other
		real(dp)	:: wage_hr=1.,al_hr=1., junk=1.,a_hr=1., e_hr=1., z_hr=1., alwt=1., ziwt=1., jwt=1., cumval=1., &
					&	work_dif_hr=1., app_dif_hr=1.,js_ij=1., Nworkt=1., ep_hr=1.,apc_hr = 1., sepi=1.,fndi = 1., hlthprob,al_last_invol,&
					&   triwt=1., ewt=0.,hatsig2,a_residj,e_residj, welfare_dif_hr = 1., wtr_avg =0.


		integer :: ali_hr=1,aliH=1,d_hr=1,age_hr=1,del_hr=1, zi_hr=1, ziH=1,il_hr=1 ,j_hr=1,ial=1, ai_hr=1,api_hr=1,ei_hr=1,triH=1,eiH=1, &
			& tri=1, tri_hr=1,fndi_hr(nz),sepi_hr(nz),status_hr=1,status_tmrw=1,drawi=1,drawt=1, invol_un = 0,slice_len=1, brn_yr_hr=1, interp_i, &
			samplestep = 4

		logical :: w_strchng_old = .false., final_iter = .false.,occaggs_hr =.true., converged = .false.

		if(present(occaggs)) then
			occaggs_hr = occaggs
		else
			occaggs_hr = .true.
		endif

		!be sure we're running with a clean history
		call clean_hist(hst)

		!************************************************************************************************!
		! Allocate things
		!************************************************************************************************!

		iter_draws = min(maxiter,8) !globally set variable

		allocate(a_it_int(Nsim,Tsim))
		allocate(e_it(Nsim,Tsim))
		allocate(e_it_int(Nsim,Tsim))
		allocate(work_it(Nsim,Tsim))
		allocate(app_it(Nsim,Tsim))
		allocate(al_int_it_endog(Nsim,Tsim))
		allocate(al_it_endog(Nsim,Tsim))
		allocate(trX_it(Nsim,Tsim))
		allocate(brn_drawi_drawt(Nsim,Tsim,2))
		allocate(invol_it(Nsim,Tsim))
		allocate(tr_hist_vec(Nsim))
		!allocate(hlthvocSSDI(Nsim,Tsim))

		!!!!!!!!!!!!!!! DEBUGGING
		allocate(ewt_it(Nsim,Tsim))
		ewt_it = 0._dp
        allocate(trindif_it(Nsim,Tsim))
        trindif_it = 0._dp
		!************************************************************************************************!
		! Pointers
		!************************************************************************************************!
		! (disability extent, earn hist, assets)

		V => vfs%V
		V_CF => vfs%V_CF
		aR => pfs%aR
		aD => pfs%aD
		aN => pfs%aN
		aW => pfs%aW
		aU => pfs%aU
		gwork => pfs%gwork
		gapp => pfs%gapp

		gapp_dif    => pfs%gapp_dif
		gwork_dif   => pfs%gwork_dif


		z_jt_innov  => shk%z_jt_innov
		z_jt_select => shk%z_jt_select
		del_i_int   => shk%del_i_int
		fndsep_i_int   => shk%fndsep_i_int
		dead_draw   => shk%dead_draw
		j_i         => shk%j_i
		al_it       => shk%al_hist
		al_int_it	=> shk%al_int_hist
		age_it 	    => shk%age_hist
		born_it	    => shk%born_hist
		status_it_innov => shk%status_it_innov
		health_it_innov => shk%health_it_innov
		d_it            => shk%d_hist
		fndarrive_draw  => shk%fndarrive_draw
		jshock_ij  	 	=> shk%jshock_ij
		drawi_ititer    => shk%drawi_ititer
		drawt_ititer    => shk%drawt_ititer


		status_it   => hst%status_hist
		work_dif_it => hst%work_dif_hist
		app_dif_it  => hst%app_dif_hist
		di_prob_it 	=> hst%di_prob_hist
		z_jt_macroint  => hst%z_jt_macroint
		z_jt_panel  => hst%z_jt_panel
		a_it        => hst%a_hist
		occgrow_jt  => hst%occgrow_jt
		occshrink_jt=> hst%occshrink_jt
		occsize_jt  => hst%occsize_jt

		ptrsuccess = associated(a_it,hst%a_hist)
		if(verbose>1) then
			if(ptrsuccess .eqv. .false. ) print *, "failed to associate a_it"
			ptrsuccess = associated(age_it,shk%age_hist)
			if(ptrsuccess .eqv. .false. ) print *, "failed to associate age_it"
			ptrsuccess = associated(V,vfs%V)
			if(ptrsuccess .eqv. .false. ) print *, "failed to associate V"
		endif
		ndets = 2 + nd-1 + oldN+1-1 + ndi-1 + nal-1+1 !status dummies (WUN), d dummies (nd), age dummies (1+oldN), del dummies (ndi), al_int_dummies(nal).
		nregobs = ceiling(dble(Nsim*Tsim)/dble(samplestep))
		allocate(cov_coef(ndets,ndets))
		allocate(acoef(ndets))
		allocate(ecoef(ndets))
		allocate(Xdets(nregobs,ndets))
		allocate(adep(nregobs))
		allocate(edep(nregobs))
		allocate(xidets(ndets))
		allocate(xjdets(ndets))
		nregobs = 0


		a_ap_dist = 0._dp
		hst%wage_hist    = 0.
		Tret = (Longev - youngD - oldD*oldN)*tlen
		work_dif_it      = 0.
		app_dif_it       = 0.
		hst%hlth_voc_hist= 0
		hst%welfare_hist = 0.
		Ncol = size(drawi_ititer,2)
		al_int_it_endog  = al_int_it
		al_it_endog      = al_it
		shk%wtr_it = 0.
		trX_it = 0.
		if(shk%drawn /= 1 )then
			call draw_shocks(shk)
		endif
		do i=1,Nsim
			do it=1,Tsim
				brn_drawi_drawt(i,it,:) = (/ i,it /)
			enddo
		enddo

		!set up cumpid,cumptau
		cumptau = 0.
		it = 1
		cumptau(it+1)=cumptau(it)
		do it =2,TT
			cumptau(it+1) = cumptau(it)+ptau(it)
		enddo

		! will draw these from endogenous distributions the second time around
		a_it = agrid(1)
		a_it_int = 1
		e_it = 0._dp !if not yet alive, then 0
		e_it_int = 1
		status_it = 1 ! just to initialize on the first round
		do i=1,Nsim
			do it=1,Tsim
				if(age_it(i,it)>=TT) status_it(i,it) = 5
				if(age_it(i,it)<= 0) status_it(i,it) = 0
			enddo
		enddo

		a_mean_liter = 0.
		d_mean_liter = 0.
		s_mean_liter = 0.
		a_var_liter  = 0.
		d_var_liter  = 0.

		tri = tri0

		if(verbose >2) print *, "Simulating"
		w_strchng_old = w_strchng
		w_strchng = .false.
		if (iter_draws>1) then
			final_iter = .false.
		else
			final_iter = .true.
		endif
		converged = .false.
		!itertate to get dist of asset/earnings correct at each age from which to draw start conditions
		do iter=1,iter_draws

			di_prob_it = 0.
			app_dif_it = 0.
			work_dif_it = 0.
			hst%hlthprob_hist = 0.
			hst%hlth_voc_hist= 0

			it = 1
			nomatch = 0
			junk =0.
			do i =1,Nsim
				!for the population that is pre-existing in the first period , it=1 and age>0
				!need to draw these from age-specific distributions for iterations > 1
				age_hr = age_it(i,it)
				del_hr = del_i_int(i)
				if(age_hr>0) then
					!use status_it_innov(i,Tsim) to identify the d state, along with age of this guy

					d_hr = d_it(i,it)

					if(iter >1) then

						do ii=1,(Ncol-1)
							drawi = drawi_ititer(i,ii) !drawi = drawi_ititer(i,mod(ii+iter-2,Ncol)+1)
							drawt = drawt_ititer(i,ii) !drawt = drawt_ititer(i,mod(ii+iter-2,Ncol)+1)

							if( (status_it(drawi,drawt)>0) .and. (status_it(drawi,drawt)<4) .and. &
							&	(age_it(drawi,drawt) .eq. age_hr ).and. (d_it(drawi,drawt) .eq. d_hr) ) then

								a_residj = 0._dp
								e_residj = 0._dp

								brn_drawi_drawt(i,it,:) = (/drawi,drawt/)
								status_it(i,it) = status_it(drawi,drawt)


								a_it(i,it)      = a_it(drawi,drawt)
								e_it(i,it)      = e_it(drawi,drawt)
								e_it_int(i,it)  = e_it_int(drawi,drawt)
								a_it_int(i,it)  = a_it_int(drawi,drawt)

								invol_un        = invol_it(drawi,drawt)
								brn_drawi_drawt(i,it,:) = (/drawi,drawt/)
								exit
							elseif(ii>=Ncol-1)then

								status_it(i,it) = 1
								invol_un = 0
								invol_it(drawi,drawt) = invol_un
								d_hr = d_it(i,it)
								brn_drawi_drawt(i,it,:) = (/i,it/)

								if(age_it(i,it)==1) then
									a_it_int(i,it) = 1
									a_it(i,it) = minval(agrid)
									e_it(i,it) = min( max(dexp(al_it(i,it)), minval(egrid)), maxval(egrid))
									e_it_int(i,it) = locate(egrid, e_it(i,it))
								else
									a_it_int(i,it) = na/2
									a_it(i,it) = agrid(na/2)
									e_it(i,it) = min( max(dexp(al_it(i,it)), minval(egrid)), maxval(egrid))
									e_it_int(i,it) = locate(egrid, e_it(i,it))
								endif
								if(iter>1) then
									nomatch = nomatch+1
									if(verbose>0) print *, "t=0 nomatch: d = ", d_hr , "al = ", al_int_it(i,it), "age = ", age_hr
								endif
								exit
							endif
						enddo
					elseif( iter==1 ) then
						status_it(i,it) = 1
						invol_un = 0
						invol_it(drawi,drawt) = invol_un
						d_hr = d_it(i,it)
						brn_drawi_drawt(i,it,:) = (/i,it/)
						if(age_it(i,it)==1) then
							a_it_int(i,it) = 1
							a_it(i,it) = minval(agrid)
							e_it(i,it) = dexp(al_it(i,it))
							e_it_int(i,it) = locate(egrid,e_it(i,it))
						else
							a_it_int(i,it) = na/2
							a_it(i,it) = agrid(na/2)
							e_it(i,it) = al_it(i,it)
							e_it_int(i,it) = locate(egrid,e_it(i,it))
						endif
					endif
				endif !age>0
			enddo !i=1:Nsim

			if( (verbose>0) .and. (nomatch>0) ) print *, "did not find match for t=0 draw ", nomatch, " times"
			nomatch = 0

			!$OMP  parallel do &
			!$OMP& private(i,interp_i,del_hr,j_hr,status_hr,it,it_old,age_hr,al_hr,ali_hr,d_hr,e_hr,a_hr,ei_hr,ai_hr,z_hr,zi_hr,api_hr,tri_hr,apc_hr,ep_hr, &
			!$OMP& ewt, eiH, aliH, alwt, ziwt,ziH,triwt,triH,il,fndi_hr, sepi_hr, il_hr,cumval,jwt,wage_hr,al_last_invol,junk,app_dif_hr,work_dif_hr, welfare_dif_hr, &
			!$OMP& hlthprob,wtr_avg,ii,drawi,drawt,sepi,fndi,invol_un,dead,status_tmrw,brn_yr_hr)
			do i=1,Nsim
				!fixed traits

				!set a j to correspond to the probabilities.  This will get overwritten if born
				j_hr = j_i(i)
				del_hr = del_i_int(i)
				fndi_hr = fndsep_i_int(i,1,:)
				sepi_hr = fndsep_i_int(i,2,:)
				il_hr = fndi_hr(2)
				!initialize stuff
				it_old = 1
				invol_un = 0
				dead = .false.
				do it=1,Tsim
				if(age_it(i,it) > 0 ) then !they've been born

					!set the average wage for the xifunction

					if(it==1) brn_yr_hr=it

					if( (born_it(i,it) .eq. 1) .and. (it> 1) ) then
					! note: no one is ``born'' in the first period

					! draw state from distribution of age 1
						age_hr	= 1
						brn_yr_hr = it

						d_hr = d_it(i,it)
						do ii=1,(Ncol-1)
							! if(ii<(Ncol-1) .and. iter>1) then
							! 	drawi = drawi_ititer(i,ii) !drawi = drawi_ititer(i,mod(ii+iter-2,Ncol)+1)
							! 	drawt = drawt_ititer(i,ii) !drawt = drawt_ititer(i,mod(ii+iter-2,Ncol)+1)
							! 	if((age_it(drawi,drawt) .eq. 1 ) .and. &
							! 	&  (status_it(drawi,drawt) .gt. 0) .and. (status_it(drawi,drawt) .le. 3) .and. &
							! 	&  (d_it(drawi,drawt) .eq. d_hr)  ) then
							! !	&  (al_int_it(drawi,drawt) .eq. al_int_it(i,it) ) ) then
							!
							! 		!ndets = 3-1+nd-1+(oldN+1-1)+ndi-1+nal-1+1 !status dummies (WUN), d dummies, age dummies (1+oldN), del dummies (12), al_int dummies.
							! 		xidets = 0._dp
							! 		xjdets = 0._dp
							! 		m=1
							! 		do il=2,3
							! 			if(status_it(drawi,drawt) == il) xidets(m) = 1._dp
							! 			if(status_it(drawi_ititer(i,ii+1),drawt_ititer(i,ii+1)) == il) xjdets(m) = 1._dp
							! 			m = m+1
							! 		enddo
							! 		do id=2,nd
							! 			if(d_it(drawi,drawt)== id) xidets(m) = 1._dp
							! 			if(d_it(drawi_ititer(i,ii+1),drawt_ititer(i,ii+1)) == id) xjdets(m) = 1._dp
							! 			m = m+1
							! 		enddo
							! 		do il=2,(oldN+1)
							! 			if(age_it(drawi,drawt) == il) xidets(m) = 1._dp
							! 			if(age_it(drawi_ititer(i,ii+1),drawt_ititer(i,ii+1)) == il) xjdets(m) = 1._dp
							! 			m=m+1
							! 		enddo
							! 		do  idi=2,ndi
							! 			if(del_i_int(drawi)==idi) xidets(m)= 1._dp
							! 			if(del_i_int(drawi_ititer(i,ii+1))==idi) xjdets(m)= 1._dp
							! 			m=m+1
							! 		enddo
							! 		do  ial=1,nal
							! 			if(ial .ne. 3)then
							! 				if(al_int_it_endog(drawi,drawt)==ial) xidets(m)= 1._dp
							! 				if(al_int_it_endog(drawi_ititer(i,ii+1),drawt_ititer(i,ii+1))==ial) xjdets(m)= 1._dp
							! 				m=m+1
							! 			endif
							! 		enddo
							! 		xidets(m)= 1._dp
							! 		xjdets(m)= 1._dp
							!
							! 		! a_residj = a_it(drawi_ititer(i,ii+1),drawt_ititer(i,ii+1)) - dot_product(xjdets,acoef)
							! 		! e_residj = e_it(drawi_ititer(i,ii+1),drawt_ititer(i,ii+1)) - dot_product(xjdets,ecoef)
							! 		a_residj = 0._dp
							! 		e_residj = 0._dp
							!
							! 		brn_drawi_drawt(i,it,:) = (/drawi, drawt /)
							! 		status_it(i,it) = status_it(drawi,drawt)
							! 		d_it(i,it)		= d_hr
							! 		! a_it(i,it)		= dot_product(xidets,acoef) + a_residj
							! 		! a_it(i,it)      = min(max(a_it(i,it),minval(agrid)),maxval(agrid))
							! 		! e_it(i,it)      = dot_product(xidets,ecoef) + e_residj
							! 		! e_it(i,it)      = min(max(e_it(i,it),minval(egrid)),maxval(egrid))
							! 		! e_it_int(i,it)  = locate(egrid,e_it(i,it))
							! 		! a_it_int(i,it)  = locate(agrid,a_it(i,it))
							! 		a_it(i,it)      = a_it(drawi,drawt)
							! 		e_it(i,it)      = e_it(drawi,drawt)
							! 		e_it_int(i,it)  = e_it_int(drawi,drawt)
							! 		a_it_int(i,it)  = a_it_int(drawi,drawt)
							! 		status_it(i,it) = status_it(drawi,drawt)
							! 		exit
							! 	endif
							! else
								! brn_drawi_drawt(i,it,1) = i
								! brn_drawi_drawt(i,it,2) = it
								! if(iter >1) then
								! 	a_it(i,it) = dot_product(xidets,acoef) !minval(agrid)
								! 	e_it(i,it) = dot_product(xidets,ecoef)
								! 	e_it_int(i,it) = locate(egrid,e_it(i,it))
								! 	a_it_int(i,it) = locate(agrid,a_it(i,it))
								! else
									a_it(i,it) = minval(agrid)
									e_it(i,it) = min(max(dexp(al_it(i,it)),minval(egrid) ) ,maxval(egrid))
									a_it_int(i,it) = 1
									e_it_int(i,it) = locate(egrid,e_it(i,it))
									exit
								! endif
								status_it(i,it) = 1
								! if(ii>=Ncol-1) then
								! 	nomatch = nomatch+1
								! 	if(verbose >0) print *, "missed at d_hr: ", d_hr, "t= ", it, "al_int ", al_int_it(i,it), "del_i_int ",del_i_int(i), " j_i= ", j_hr
								! endif
								! exit
							! endif
						enddo
					endif !if make decisions when first born?
					! load state - may have been set earlier in the iteration if they're born in the 1st period
					age_hr	= age_it(i,it)
					d_hr	= d_it(i,it)
					a_hr 	= a_it(i,it)
					ei_hr	= e_it_int(i,it)
					e_hr 	= e_it(i,it)
					ai_hr 	= a_it_int(i,it)
					status_hr = status_it(i,it)



					! get set to kill off old (i.e. age_hr ==TT only for Longev - youngD - oldD*oldN )
					if((age_hr .eq. TT) ) then !
						it_old = it_old + 1
						if(it_old >  Tret ) then !
							dead = .true.
						else
							status_it(i,it) = 5
						endif
					endif

					if(dieyoung .eqv. .true.)then
						if(dead_draw(i,it)< PrDeath(d_hr,age_hr)) & !this is using the back end of status: if it's very high, then dead
							dead = .true.
					endif

					if(dead .eqv. .true.) then
							a_it(i,it) = 0.
							a_it_int(i,it) = 1
							app_dif_it(i,it) = 0.
							work_dif_it(i,it) = 0.
							status_it(i,it) = -1
							dead = .true.
							cycle
						!	exit
					endif

					!figure out where to evaluate z
					if(zj_contin .eqv. .false.) then
						if(buscyc .eqv. .true.) then
							zi_hr	= z_jt_macroint(it)
						else
							zi_hr	= 2
						endif
						z_hr	= zgrid(zi_hr,j_hr)
						il_hr   = fndi_hr(zi_hr)
					else
						z_hr	= z_jt_panel(j_hr,it)
						do zi_hr = nz,1,-1
							if(zgrid(zi_hr,j_hr)<z_hr) exit
						enddo
						zi_hr = min( max(zi_hr,1), nz )
						ziH  = min(zi_hr+1,nz)
						ziwt = (zgrid(ziH,j_hr)- z_hr)/( zgrid(ziH,j_hr) -   zgrid(zi_hr,j_hr) )
						if( ziH == zi_hr ) ziwt = 1.
					endif
					
					
					!set the idiosyncratic income state
					al_hr	= al_it(i,it)
					ali_hr	= al_int_it(i,it)
					if( (born_it(i,it) .eq. 1) .or. ((it .eq. 1) .and. (age_it(i,it) .gt. 0)) ) then
						drawi = brn_drawi_drawt(i,it,1)
						drawt = brn_drawi_drawt(i,it,2)
						if(al_int_it_endog(drawi,drawt) .eq. 1) then
							invol_un = 1
							invol_it(i,it) = invol_un
							al_last_invol = al_it(drawi,drawt) ! have to set this to something. Given persistence, it's likely true
						endif
					endif
					if( (w_strchng .eqv. .true.) .and. (it > TossYears*itlen ) )then
						shk%wtr_it(i,it) = wage_trend(it,j_hr) + wage_lev(j_hr)
					else
						shk%wtr_it(i,it) = wage_lev(j_hr)
					endif
					tri_hr = finder(trgrid,shk%wtr_it(i,it))

					tri_hr = min( max(tri_hr, 1), ntr )
					triH = min(tri_hr+1,ntr)
					if(triH>tri_hr) then
						triwt = (trgrid(triH)- shk%wtr_it(i,it))/(trgrid(triH) - trgrid(tri_hr))
					else
						triwt = 1._dp
					endif

					if((invol_un .eq. 1) .and. (status_hr .ne. 1))then !have already been born
						al_hr	= alfgrid(1,d_hr)
						ali_hr	= 1
						aliH    = 2
						alwt   = 0._dp
					else
						invol_it(i,it) = 0
						al_hr	= al_it(i,it)
						ali_hr	= al_int_it(i,it)
					endif

					!figure out where to evaluate alpha
					if(al_contin .eqv. .true.) then
						aliH  = max(1, min(ali_hr+1,nal))
						if((ali_hr>1) .and. (ali_hr .ne. aliH)) then
							alwt = (alfgrid(aliH,d_hr)- al_hr)/( alfgrid(aliH,d_hr) -   alfgrid(ali_hr,d_hr) )
						else
							alwt = 0._dp
						endif
					endif
					!where to evaluate e
					ei_hr = max(min(ei_hr,ne),1)
					eiH  = min(ei_hr+1,ne)
					if(ei_hr .ne. eiH) then
						ewt = (egrid(eiH)- e_hr)/( egrid(eiH) -   egrid(ei_hr) )
					else !at the max
						ewt = 1._dp
					endif
					ewt_it(i,it) = ewt

					junk = 0._dp
					if(w_strchng .eqv. .true.) junk = shk%wtr_it(i,it)
					if(junk == 1000._dp) print *, "big wage trend!"
					if( age_hr <TT ) then
						wage_hr	= wage(junk,al_hr,d_hr,age_hr)
					else
						wage_hr	= wage(junk,al_hr,d_hr,TT-1)
					endif
					if(invol_un ==1 .and. age_hr<TT) then
						!carry forward their old wage if they are now involuntarily unemployed
						wage_hr	= wage(junk,al_last_invol,d_hr,age_hr)
					endif
					hst%wage_hist(i,it) = wage_hr

					!make decisions if not yet retired
                    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                    if( age_hr == TT-1 .and. status_it_innov(i,it) < 0.5_dp ) then
!
!                        if( vfs%VR(1,d_hr,ei_hr,ai_hr)*ewt+ (1._dp - ewt)*vfs%VR(1,d_hr,eiH,ai_hr) > &
!                                & (ewt*(triwt    * alwt  *V( (il_hr-1)*ntr + tri_hr, (del_hr-1)*nal+ali_hr,d_hr,ei_hr,ai_hr,zi_hr,age_hr ) + &
!                                &	       triwt    *(1.-alwt) *V( (il_hr-1)*ntr + tri_hr, (del_hr-1)*nal+aliH   ,d_hr,ei_hr,ai_hr,zi_hr,age_hr ) + &
!                                &	      (1.-triwt)* alwt     *V( (il_hr-1)*ntr + triH  , (del_hr-1)*nal+ali_hr,d_hr,ei_hr,ai_hr,zi_hr,age_hr ) + &
!                                &	      (1.-triwt)*(1.-alwt) *V( (il_hr-1)*ntr + triH  , (del_hr-1)*nal+aliH   ,d_hr,ei_hr,ai_hr,zi_hr,age_hr )) + &
!                                & (1.-ewt)*(triwt    * alwt    *V( (il_hr-1)*ntr + tri_hr,(del_hr-1)*nal+ali_hr,d_hr,eiH  ,ai_hr,zi_hr,age_hr ) + &
!                                &	  		triwt    *(1.-alwt)*V( (il_hr-1)*ntr + tri_hr,(del_hr-1)*nal+aliH   ,d_hr,eiH  ,ai_hr,zi_hr,age_hr ) + &
!                                &	 	   (1.-triwt)* alwt    *V( (il_hr-1)*ntr + triH  ,(del_hr-1)*nal+ali_hr,d_hr,eiH  ,ai_hr,zi_hr,age_hr ) + &
!                                &		   (1.-triwt)*(1.-alwt)*V( (il_hr-1)*ntr + triH  ,(del_hr-1)*nal+aliH   ,d_hr,eiH  ,ai_hr,zi_hr,age_hr )) )) then
!                            age_hr = TT
!                            status_hr = 5
!                        end if
!
!                    end if


					if(age_hr < TT) then
						!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
						! evalutate gwork and gapp to figure out lom of status
						if(status_hr <= 4) then
							!evaluate application choice for diagnostics (would the workers want to apply? even if they can't)
							if((al_contin .eqv. .true.) .and. (zj_contin .eqv. .false.) .and. (w_strchng .eqv. .false.)) then
								app_dif_hr =  ewt*(alwt    *gapp_dif( (il_hr-1)*ntr + tri_hr, (del_hr-1)*nal+ali_hr,d_hr,ei_hr,ai_hr,zi_hr,age_hr ) + &
										&		  (1.-alwt)*gapp_dif( (il_hr-1)*ntr + tri_hr, (del_hr-1)*nal+aliH   ,d_hr,ei_hr,ai_hr,zi_hr,age_hr )) + &
										&(1.-ewt)*(alwt    *gapp_dif( (il_hr-1)*ntr + tri_hr, (del_hr-1)*nal+ali_hr,d_hr,eiH  ,ai_hr,zi_hr,age_hr ) + &
										&		  (1.-alwt)*gapp_dif( (il_hr-1)*ntr + tri_hr, (del_hr-1)*nal+aliH   ,d_hr,eiH  ,ai_hr,zi_hr,age_hr ))
							elseif((al_contin .eqv. .true.) .and. (zj_contin .eqv. .false.) .and. (w_strchng .eqv. .true.)) then
								app_dif_hr = ewt*  (triwt    * alwt    *gapp_dif( (il_hr-1)*ntr + tri_hr, (del_hr-1)*nal+ali_hr,d_hr,ei_hr,ai_hr,zi_hr,age_hr ) + &
										&	 	    triwt    *(1.-alwt)*gapp_dif( (il_hr-1)*ntr + tri_hr, (del_hr-1)*nal+aliH   ,d_hr,ei_hr,ai_hr,zi_hr,age_hr ) + &
										&		   (1.-triwt)*  alwt   *gapp_dif( (il_hr-1)*ntr + triH  , (del_hr-1)*nal+ali_hr,d_hr,ei_hr,ai_hr,zi_hr,age_hr ) + &
										&		   (1.-triwt)*(1.-alwt)*gapp_dif( (il_hr-1)*ntr + triH  , (del_hr-1)*nal+aliH   ,d_hr,ei_hr,ai_hr,zi_hr,age_hr ))+ &
										& (1.-ewt)*(triwt    * alwt    *gapp_dif( (il_hr-1)*ntr + tri_hr, (del_hr-1)*nal+ali_hr,d_hr,eiH  ,ai_hr,zi_hr,age_hr ) + &
										&	 	    triwt    *(1.-alwt)*gapp_dif( (il_hr-1)*ntr + tri_hr, (del_hr-1)*nal+aliH   ,d_hr,eiH  ,ai_hr,zi_hr,age_hr ) + &
										&		   (1.-triwt)*  alwt   *gapp_dif( (il_hr-1)*ntr + triH  , (del_hr-1)*nal+ali_hr,d_hr,eiH  ,ai_hr,zi_hr,age_hr ) + &
										&		   (1.-triwt)*(1.-alwt)*gapp_dif( (il_hr-1)*ntr + triH  , (del_hr-1)*nal+aliH   ,d_hr,eiH  ,ai_hr,zi_hr,age_hr ))
							else ! if((al_contin .eqv. .false.) .and. (zj_contin .eqv. .false.) ) then
								app_dif_hr = gapp_dif( (il_hr-1)*ntr + tri_hr, (del_hr-1)*nal+ali_hr,d_hr,ei_hr,ai_hr,zi_hr,age_hr )
							endif

							app_dif_it(i,it) = app_dif_hr

							!check for rest unemployment
							if((al_contin .eqv. .true.) .and. (zj_contin .eqv. .false.) .and. (w_strchng .eqv. .false.)) then
								work_dif_hr= ewt*   (alwt*gwork_dif( (il_hr-1)*ntr + tri_hr, (del_hr-1)*nal+ali_hr,d_hr,ei_hr,ai_hr,zi_hr,age_hr ) + &
										&	    (1.-alwt)*gwork_dif( (il_hr-1)*ntr + tri_hr, (del_hr-1)*nal+aliH   ,d_hr,ei_hr,ai_hr,zi_hr,age_hr ) )+ &
										& (1.-ewt)* (alwt*gwork_dif( (il_hr-1)*ntr + tri_hr, (del_hr-1)*nal+ali_hr,d_hr,eiH,ai_hr,zi_hr,age_hr ) + &
										&	    (1.-alwt)*gwork_dif( (il_hr-1)*ntr + tri_hr, (del_hr-1)*nal+aliH   ,d_hr,eiH,ai_hr,zi_hr,age_hr ))
							elseif((al_contin .eqv. .true.) .and. (zj_contin .eqv. .false.) .and. (w_strchng .eqv. .true.)) then
								work_dif_hr = ewt*(triwt    * alwt    *gwork_dif( (il_hr-1)*ntr + tri_hr, (del_hr-1)*nal+ali_hr,d_hr,ei_hr,ai_hr,zi_hr,age_hr ) + &
										&	       triwt    *(1.-alwt)*gwork_dif( (il_hr-1)*ntr + tri_hr, (del_hr-1)*nal+aliH   ,d_hr,ei_hr,ai_hr,zi_hr,age_hr ) + &
										&	      (1.-triwt)* alwt    *gwork_dif( (il_hr-1)*ntr + triH  , (del_hr-1)*nal+ali_hr,d_hr,ei_hr,ai_hr,zi_hr,age_hr ) + &
										&	      (1.-triwt)*(1.-alwt)*gwork_dif( (il_hr-1)*ntr + triH  , (del_hr-1)*nal+aliH   ,d_hr,ei_hr,ai_hr,zi_hr,age_hr )) + &
										& (1.-ewt)*(triwt    * alwt    *gwork_dif( (il_hr-1)*ntr + tri_hr,(del_hr-1)*nal+ali_hr,d_hr,eiH  ,ai_hr,zi_hr,age_hr ) + &
										&	  		triwt    *(1.-alwt)*gwork_dif( (il_hr-1)*ntr + tri_hr,(del_hr-1)*nal+aliH   ,d_hr,eiH  ,ai_hr,zi_hr,age_hr ) + &
										&	 	   (1.-triwt)* alwt    *gwork_dif( (il_hr-1)*ntr + triH  ,(del_hr-1)*nal+ali_hr,d_hr,eiH  ,ai_hr,zi_hr,age_hr ) + &
										&		   (1.-triwt)*(1.-alwt)*gwork_dif( (il_hr-1)*ntr + triH  ,(del_hr-1)*nal+aliH   ,d_hr,eiH  ,ai_hr,zi_hr,age_hr ))

							else  !if(al_contin .eqv. .false. ) then
								work_dif_hr = gwork_dif( (il_hr-1)*ntr + tri_hr, (del_hr-1)*nal+ali_hr,d_hr,ei_hr,ai_hr,zi_hr,age_hr )
							endif
							work_dif_it(i,it) = work_dif_hr

							!compute welfare difference if requested
							if(welfare_cf .eqv. .true.) then
								if((al_contin .eqv. .true.) .and. (zj_contin .eqv. .false.) .and. (w_strchng .eqv. .false.)) then
									welfare_dif_hr= (ewt*   (alwt*V( (il_hr-1)*ntr + tri_hr, (del_hr-1)*nal+ali_hr,d_hr,ei_hr,ai_hr,zi_hr,age_hr ) + &
											&	    (1.-alwt)*V( (il_hr-1)*ntr + tri_hr, (del_hr-1)*nal+aliH   ,d_hr,ei_hr,ai_hr,zi_hr,age_hr ) )+ &
											& (1.-ewt)* (alwt*V( (il_hr-1)*ntr + tri_hr, (del_hr-1)*nal+ali_hr,d_hr,eiH,ai_hr,zi_hr,age_hr ) + &
											&	    (1.-alwt)*V( (il_hr-1)*ntr + tri_hr, (del_hr-1)*nal+aliH   ,d_hr,eiH,ai_hr,zi_hr,age_hr )) )&
											&  / &
											& (   ewt*  (alwt*V_CF( (il_hr-1)*ntr + tri_hr, (del_hr-1)*nal+ali_hr,d_hr,ei_hr,ai_hr,zi_hr,age_hr ) + &
											&	    (1.-alwt)*V_CF( (il_hr-1)*ntr + tri_hr, (del_hr-1)*nal+aliH   ,d_hr,ei_hr,ai_hr,zi_hr,age_hr ) )+ &
											& (1.-ewt)* (alwt*V_CF( (il_hr-1)*ntr + tri_hr, (del_hr-1)*nal+ali_hr,d_hr,eiH,ai_hr,zi_hr,age_hr ) + &
											&	    (1.-alwt)*V_CF( (il_hr-1)*ntr + tri_hr, (del_hr-1)*nal+aliH   ,d_hr,eiH,ai_hr,zi_hr,age_hr )) )
								elseif((al_contin .eqv. .true.) .and. (zj_contin .eqv. .false.) .and. (w_strchng .eqv. .true.)) then
									welfare_dif_hr =(ewt*(triwt    * alwt  *V( (il_hr-1)*ntr + tri_hr, (del_hr-1)*nal+ali_hr,d_hr,ei_hr,ai_hr,zi_hr,age_hr ) + &
											&	       triwt    *(1.-alwt) *V( (il_hr-1)*ntr + tri_hr, (del_hr-1)*nal+aliH   ,d_hr,ei_hr,ai_hr,zi_hr,age_hr ) + &
											&	      (1.-triwt)* alwt     *V( (il_hr-1)*ntr + triH  , (del_hr-1)*nal+ali_hr,d_hr,ei_hr,ai_hr,zi_hr,age_hr ) + &
											&	      (1.-triwt)*(1.-alwt) *V( (il_hr-1)*ntr + triH  , (del_hr-1)*nal+aliH   ,d_hr,ei_hr,ai_hr,zi_hr,age_hr )) + &
											& (1.-ewt)*(triwt    * alwt    *V( (il_hr-1)*ntr + tri_hr,(del_hr-1)*nal+ali_hr,d_hr,eiH  ,ai_hr,zi_hr,age_hr ) + &
											&	  		triwt    *(1.-alwt)*V( (il_hr-1)*ntr + tri_hr,(del_hr-1)*nal+aliH   ,d_hr,eiH  ,ai_hr,zi_hr,age_hr ) + &
											&	 	   (1.-triwt)* alwt    *V( (il_hr-1)*ntr + triH  ,(del_hr-1)*nal+ali_hr,d_hr,eiH  ,ai_hr,zi_hr,age_hr ) + &
											&		   (1.-triwt)*(1.-alwt)*V( (il_hr-1)*ntr + triH  ,(del_hr-1)*nal+aliH   ,d_hr,eiH  ,ai_hr,zi_hr,age_hr )) ) &
											&  / &
											&	    (ewt*(triwt    * alwt  *V_CF( (il_hr-1)*ntr + tri_hr, (del_hr-1)*nal+ali_hr,d_hr,ei_hr,ai_hr,zi_hr,age_hr ) + &
											&	       triwt    *(1.-alwt) *V_CF( (il_hr-1)*ntr + tri_hr, (del_hr-1)*nal+aliH   ,d_hr,ei_hr,ai_hr,zi_hr,age_hr ) + &
											&	      (1.-triwt)* alwt     *V_CF( (il_hr-1)*ntr + triH  , (del_hr-1)*nal+ali_hr,d_hr,ei_hr,ai_hr,zi_hr,age_hr ) + &
											&	      (1.-triwt)*(1.-alwt) *V_CF( (il_hr-1)*ntr + triH  , (del_hr-1)*nal+aliH   ,d_hr,ei_hr,ai_hr,zi_hr,age_hr )) + &
											& (1.-ewt)*(triwt    * alwt    *V_CF( (il_hr-1)*ntr + tri_hr,(del_hr-1)*nal+ali_hr,d_hr,eiH  ,ai_hr,zi_hr,age_hr ) + &
											&	  		triwt    *(1.-alwt)*V_CF( (il_hr-1)*ntr + tri_hr,(del_hr-1)*nal+aliH   ,d_hr,eiH  ,ai_hr,zi_hr,age_hr ) + &
											&	 	   (1.-triwt)* alwt    *V_CF( (il_hr-1)*ntr + triH  ,(del_hr-1)*nal+ali_hr,d_hr,eiH  ,ai_hr,zi_hr,age_hr ) + &
											&		   (1.-triwt)*(1.-alwt)*V_CF( (il_hr-1)*ntr + triH  ,(del_hr-1)*nal+aliH   ,d_hr,eiH  ,ai_hr,zi_hr,age_hr )) )

								else  !if(al_contin .eqv. .false. ) then
									welfare_dif_hr = V( (il_hr-1)*ntr + tri_hr, (del_hr-1)*nal+ali_hr,d_hr,ei_hr,ai_hr,zi_hr,age_hr ) &
										& / &
										& V_CF( (il_hr-1)*ntr + tri_hr, (del_hr-1)*nal+ali_hr,d_hr,ei_hr,ai_hr,zi_hr,age_hr )
								endif
								hst%welfare_hist(i,it) = welfare_dif_hr

							endif

							!store the disability probability for the first group
							wtr_avg = shk%wtr_it(i,it)
							!wtr_avg = wage_hr
							di_prob_it(i,it) = xifun(d_hr,wtr_avg,age_hr, wtr_tmean_ts(it),hlthprob, max(it-TossYears*itlen,0))
                            trindif_it(i,it) = wtr_avg - wtr_tmean_ts(it)
							!draws for exog find and sep
							if( zj_contin .eqv. .true.) then
								fndi = ziwt*fndgrid(fndi_hr(zi_hr),zi_hr) + (1.-ziwt)*fndgrid(fndi_hr(ziH),ziH)
								sepi = ziwt*sepgrid(sepi_hr(zi_hr),zi_hr) + (1.-ziwt)*sepgrid(sepi_hr(ziH),ziH)
							else
								fndi = fndgrid(fndi_hr(zi_hr),zi_hr) !fndrate(zi_hr,j_hr)
								sepi = sepgrid(sepi_hr(zi_hr),zi_hr)!seprisk(zi_hr,j_hr)
							endif
							if (fndsep_seq .eqv. .true.) then 
								fndi = fnd_seq(it)
								sepi = sep_seq(it) 
								if (buscyc .eqv. .false.) then 
									fndi = fnd_seq(Tsim+1)
									sepi = sep_seq(Tsim+1) 
								endif
								
							endif
							
						endif ! status<=4

						if(status_hr .le. 3) then !in the labor force
							! figure out status transition and involuntary unemployment
							select case (status_hr)
							case(1) ! working
								if( status_it_innov(i,it) < sepi .or. work_dif_hr<0 ) then !separate? (don't allow endogenous separations in period 1 because their state might be screwy)
									al_last_invol = al_hr
									wage_hr	= wage(0._dp,al_last_invol,d_hr,age_hr)
									if(w_strchng .eqv. .true.) &
										& wage_hr	= wage( shk%wtr_it(i,it),al_last_invol,d_hr,age_hr)
									if(status_it_innov(i,it) < sepi) then !separated involunarily
										ali_hr = 1
										invol_un = 1
										alwt = 0.
										aliH = 2
										al_hr = alfgrid(ali_hr,d_hr)
									endif
									status_tmrw = 2
									status_it(i,it) = 2
								else
									status_tmrw = 1
									status_it(i,it) = 1
								endif
							case(2)
								! unemployed, may stay unemployed or become long-term unemployed
								if(status_it_innov(i,it) <=pphi) then
									status_tmrw = 3
									status_it(i,it) = 2
								elseif( (invol_un == 1)  .and. (fndarrive_draw(i,it) > fndi)) then
								!voluntary or involuntary?
									status_it(i,it) = 2
									status_tmrw =2
									invol_it(i,it) = 1
								elseif((invol_un == 1) .and. (fndarrive_draw(i,it) <= fndi)) then
									invol_un = 0
									invol_it(i,it) = 0
									ali_hr	= al_int_it(i,it)
									al_hr = al_it(i,it)
									! found a job!
									status_it(i,it)= 1
									status_tmrw = 1
								elseif( work_dif_hr < 0 ) then !voluntary unemployment (implies invol_un ==0)
									status_it(i,it) = 2
									status_tmrw = 2
								else ! invol_un != 1 and work_dif>0
									status_tmrw = 1
									status_it(i,it) = 1
									status_hr = 1
								endif
							case(3) ! status_hr eq 3
								!lfstatus updates
								if( (invol_un == 1)  .and. (fndarrive_draw(i,it) > lrho*fndi)) then
								!voluntary or involuntary?
									status_it(i,it) = 3
									status_tmrw =3
									invol_it(i,it) = 1
								elseif((invol_un == 1) .and. (fndarrive_draw(i,it) <= lrho*fndi)) then
									invol_un = 0
									invol_it(i,it) = 0
									ali_hr	= al_int_it(i,it)
									al_hr = al_it(i,it)
									! found a job!
									status_it(i,it)= 3
									status_tmrw = 1
								else ! invol_un != 1 and work_dif>0
									status_tmrw = 1
									status_it(i,it) = 3
									status_hr = 3
								endif
								if( work_dif_hr < 0 ) then !voluntary unemployment (implies invol_un ==0)
									status_it(i,it) = 3
									status_tmrw =3
								endif

								app_dif_it(i,it) = app_dif_hr

								if( app_dif_hr < 0 ) app_it(i,it) = 0
								if(app_dif_hr >= 0) then
									! choose to apply
									app_it(i,it) = 1

									!applying, do you get it?
									!record probabilities
									wtr_avg = shk%wtr_it(i,it)
									!wtr_avg = wage_hr
									junk = xifun(d_hr,wtr_avg, age_hr, wtr_tmean_ts(it),hlthprob, max(it-TossYears*itlen,0))

									if( (age_hr .eq. 1) .and. (ineligNoNu .eqv. .true.) ) then
										di_prob_it(i,it) = junk*eligY
										hst%hlthprob_hist(i,it) = hlthprob*eligY
									else
										di_prob_it(i,it) = junk
										hst%hlthprob_hist(i,it) = hlthprob
									endif

									if(status_it_innov(i,it) < junk) then
										status_tmrw = 4
										if( status_it_innov(i,it) <  hlthprob) then
											hst%hlth_voc_hist(i,it) = 1
										else
											hst%hlth_voc_hist(i,it) = 2
										endif
									else
										status_tmrw = 3
									endif
								endif
							end select

						elseif(status_hr > 3 ) then !absorbing states of D,R
							status_tmrw = status_hr
							! just to fill in values
							!app_dif_it(i,it) = 0.
							!work_dif_it(i,it) = 0.
							if(status_hr==4) di_prob_it(i,it) = 1.
							if( invol_un == 1 .or. al_it_endog(i,it-1) == 1) then
								ali_hr = 1
								al_hr = al_last_invol
								invol_it(i,it) = 1
							endif
						endif
						!evaluate the asset policy
						if(status_hr .eq. 4) then
							api_hr = aD( d_hr,ei_hr,ai_hr,age_hr )
							apc_hr = agrid(api_hr)
						elseif(status_hr .eq. 5) then ! should never be in this condition
							api_hr = aR(2,d_hr,ei_hr,ai_hr)
							apc_hr = agrid(api_hr)
						else
							do interp_i = 1,2
								if( interp_i ==2) then
									ii = ei_hr
									ei_hr = min(ei_hr+1,ne)
								endif
								!INTERPOLATE!!!!!
								if((al_contin .eqv. .true.)  .and. (zj_contin .eqv. .false.) .and. (w_strchng .eqv. .false.)) then
									if(status_hr .eq. 1) then
										apc_hr = alwt    *agrid( aw( (il_hr-1)*ntr + tri_hr, (del_hr-1)*nal+ali_hr,d_hr,ei_hr,ai_hr,zi_hr,age_hr  )) +&
											&	(1.-alwt)*agrid( aw( (il_hr-1)*ntr + tri_hr, (del_hr-1)*nal+aliH   ,d_hr,ei_hr,ai_hr,zi_hr,age_hr  ))
									elseif(status_hr .eq. 2) then
										apc_hr = alwt    *agrid( aU( (il_hr-1)*ntr + tri_hr, (del_hr-1)*nal+ali_hr,d_hr,ei_hr,ai_hr,zi_hr,age_hr )) +&
											&	(1.-alwt)*agrid( aU( (il_hr-1)*ntr + tri_hr, (del_hr-1)*nal+aliH   ,d_hr,ei_hr,ai_hr,zi_hr,age_hr ))
									elseif(status_hr .eq. 3) then
										apc_hr = alwt    *agrid( aN( (il_hr-1)*ntr + tri_hr, (del_hr-1)*nal+ali_hr,d_hr,ei_hr,ai_hr,zi_hr,age_hr )) +&
											&	(1.-alwt)*agrid( aN( (il_hr-1)*ntr + tri_hr, (del_hr-1)*nal+aliH   ,d_hr,ei_hr,ai_hr,zi_hr,age_hr ))
									endif
								elseif((al_contin .eqv. .true.)  .and. (zj_contin .eqv. .false.) .and. (w_strchng .eqv. .true.)) then
									if(status_hr .eq. 1) then
										apc_hr = triwt    * alwt    *agrid( aw( (il_hr-1)*ntr + tri_hr, (del_hr-1)*nal+ali_hr,d_hr,ei_hr,ai_hr,zi_hr,age_hr  )) +&
											&	 triwt    *(1.-alwt)*agrid( aw( (il_hr-1)*ntr + tri_hr, (del_hr-1)*nal+aliH   ,d_hr,ei_hr,ai_hr,zi_hr,age_hr  )) +&
											&	(1.-triwt)* alwt    *agrid( aw( (il_hr-1)*ntr + triH  , (del_hr-1)*nal+ali_hr,d_hr,ei_hr,ai_hr,zi_hr,age_hr  )) +&
											&	(1.-triwt)*(1.-alwt)*agrid( aw( (il_hr-1)*ntr + triH  , (del_hr-1)*nal+aliH   ,d_hr,ei_hr,ai_hr,zi_hr,age_hr  ))
									elseif(status_hr .eq. 2) then
										apc_hr = triwt    * alwt    *agrid( aU( (il_hr-1)*ntr + tri_hr, (del_hr-1)*nal+ali_hr,d_hr,ei_hr,ai_hr,zi_hr,age_hr )) +&
											&	 triwt    *(1.-alwt)*agrid( aU( (il_hr-1)*ntr + tri_hr, (del_hr-1)*nal+aliH   ,d_hr,ei_hr,ai_hr,zi_hr,age_hr )) +&
											&   (1.-triwt)* alwt    *agrid( aU( (il_hr-1)*ntr + triH  , (del_hr-1)*nal+ali_hr,d_hr,ei_hr,ai_hr,zi_hr,age_hr )) +&
											&	(1.-triwt)*(1.-alwt)*agrid( aU( (il_hr-1)*ntr + triH  , (del_hr-1)*nal+aliH   ,d_hr,ei_hr,ai_hr,zi_hr,age_hr ))
									elseif(status_hr .eq. 3) then
										apc_hr = triwt    * alwt    *agrid( aN( (il_hr-1)*ntr + tri_hr, (del_hr-1)*nal+ali_hr,d_hr,ei_hr,ai_hr,zi_hr,age_hr )) +&
											&	 triwt    *(1.-alwt)*agrid( aN( (il_hr-1)*ntr + tri_hr, (del_hr-1)*nal+aliH   ,d_hr,ei_hr,ai_hr,zi_hr,age_hr )) +&
											&   (1.-triwt)* alwt    *agrid( aN( (il_hr-1)*ntr + triH  , (del_hr-1)*nal+ali_hr,d_hr,ei_hr,ai_hr,zi_hr,age_hr )) +&
											&	(1.-triwt)*(1.-alwt)*agrid( aN( (il_hr-1)*ntr + triH  , (del_hr-1)*nal+aliH   ,d_hr,ei_hr,ai_hr,zi_hr,age_hr ))
									endif
								else !if((al_contin .eqv. .false.) .and. (zj_contin .eqv. .false.)) then
									if(status_hr .eq. 1) then
										api_hr = aw( (il_hr-1)*ntr + tri_hr, (del_hr-1)*nal+ali_hr,d_hr,ei_hr,ai_hr,zi_hr,age_hr  )
									elseif(status_hr .eq. 2) then
										api_hr = aU( (il_hr-1)*ntr + tri_hr, (del_hr-1)*nal+ali_hr,d_hr,ei_hr,ai_hr,zi_hr,age_hr )
									elseif(status_hr .eq. 3) then
										api_hr = aN( (il_hr-1)*ntr + tri_hr, (del_hr-1)*nal+ali_hr,d_hr,ei_hr,ai_hr,zi_hr,age_hr )
									endif
									apc_hr = agrid(api_hr)
								endif
								apc_hr = min(max( apc_hr, minval(agrid)),maxval(agrid))
								if(interp_i == 1) then
									junk = apc_hr
								else !interp_i == 2
									apc_hr = ewt*junk+(1._dp-ewt)*junk
									ei_hr = ii
								endif
							enddo !do interp_i
							api_hr=1
							do ii=2,na
								if( dabs(agrid(ii)-apc_hr) <  dabs(agrid(ii)-agrid(api_hr))) api_hr = ii
							enddo
							api_hr = min(max(api_hr,1),na)
							a_ap_dist(ai_hr,api_hr) = a_ap_dist(ai_hr,api_hr)+1._dp
						endif
					! retired
					elseif( (age_hr==TT) .and. (it_old <= Tret)) then
						api_hr      = aR(2, d_hr,ei_hr,ai_hr )
						apc_hr      = agrid(api_hr)
						status_hr   = 5
						status_tmrw = 5
						if(it<Tsim) &
						&	status_it(i,it+1) = status_hr
					endif
					al_int_it_endog(i,it) = ali_hr
					al_it_endog(i,it)     = al_hr
					if(status_hr>1) then
						select case (status_hr)
						case(2)
							trX_it(i,it) = UI(e_hr,ali_hr)
						case(3)
							trX_it(i,it) = min(LTUrr * e_hr,LTUrr)
						case(4)
							trX_it(i,it) = SSDI(e_hr)
						case(5)
							trX_it(i,it) = SSret(e_hr)
						end select
					else
						trX_it(i,it) = 0._dp
					endif

					!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
					!push forward the state:
					if(it<Tsim) then
						! push forward status
						status_it(i,it+1) = status_tmrw
						! push forward asset
						a_it_int(i,it+1) = api_hr
						if( (al_contin .eqv. .true.) .or. (zj_contin .eqv. .true.) ) then
							a_it(i,it+1) = apc_hr
						else
							a_it(i,it+1) = agrid(api_hr)
						endif
						!push forward AIME
						if(status_hr .le. 3) then

							if( it>brn_yr_hr .and. status_hr==1 ) then
								ep_hr = (e_hr*(tlen*agegrid(age_hr)-1._dp) + wage_hr)/ &
								&	(tlen*agegrid(age_hr))
								if( agegrid(age_hr) >17 ) then
									ep_hr = max(e_hr,ep_hr)
								endif
							elseif(it == brn_yr_hr) then
								ep_hr = wage_hr
							else
								ep_hr = e_hr
							endif
							ep_hr = max(min( ep_hr,egrid(ne) ),egrid(1))
							e_it(i,it+1) = ep_hr
							! assign to grid points by floor rule (will later interpolate)
							ei_hr = locate(egrid,ep_hr)
							ei_hr = min(max(ei_hr,1),ne)
							e_it_int(i,it+1) = ei_hr

						else
							e_it(i,it+1) = e_hr
							e_it_int(i,it+1) = ei_hr
						endif
						if(it== brn_yr_hr) then
							e_it_int(i,it) = e_it_int(i,it+1)
							e_it(i,it) = e_it(i,it+1)
						endif
						!push forward d from already drawn d_it
						if(it<Tsim) d_hr = d_it(i,it+1)

					endif !t<Tsim
				else !age_it(i,it) <= 0, they've not been born or they are dead
					a_it(i,it) = 0.
					a_it_int(i,it) = 1
					app_dif_it(i,it) = 0.
					work_dif_it(i,it) = 0.
					status_it(i,it) = -1
				endif ! age_it(i,it)>0
				enddo !1,Tsim
			enddo! 1,Nsim
			!$OMP  end parallel do

			if(verbose >0 .and. nomatch>0)  print *, "did match for new borns in period >1 ", nomatch, " times"

			wtr_scale = 0._dp

			do it=1,Tsim
				wtr_tmean_ts(it) = 0._dp
				junk = 0._dp
				ii = 0
				do i=1,Nsim
					j_hr = j_i(i)
					if( hst%status_hist(i,it) .ge. 1 .and. hst%status_hist(i,it) .le. 3 &
					&  .and. shk%age_hist(i,it)>0) then
						ii = ii+1
						tr_hist_vec(ii) =  wage_trend(it,j_hr) + wage_lev(j_hr)
						wtr_tmean_ts(it) = wtr_tmean_ts(it)+  wage_trend(it,j_hr) + wage_lev(j_hr) !hst%wage_hist(i,it)


					endif
				enddo
				wtr_tmean_ts(it) = wtr_tmean_ts(it)/dble(ii)
				!if(it> (2002-1980)*12 ) wtr_tmean_ts(it) = maxval(trgrid) !wtr_tmean_ts(it-1)
				wtr_avg = 0._dp ! actually the average sum of suqares here
				do i=1,ii
					wtr_avg = (tr_hist_vec(i)-wtr_tmean_ts(it))**2 + wtr_avg
				enddo
				wtr_avg = wtr_avg/dble(ii)
				wtr_scale = wtr_avg + wtr_scale

				call quicksort(tr_hist_vec,1,ii)
				i = floor(  dble(ii-1)/2._dp )
				wtr_tmean_ts(it) = tr_hist_vec(i)
			enddo
			wtr_scale =  wtr_avg/dble(Tsim)
			wtr_scale = wtr_scale**0.5_dp

			if(print_lev>=2) print *, "wtr_scale is : ", wtr_scale

			if(print_lev >=3)then

				write( buffer, '(I3.3)' ) iter
				call vec2csv(wtr_tmean_ts,"wtr_tmean_iter"//trim(buffer)//".csv")
				call mat2csv (e_it,"e_it.csv")
				call mat2csv (a_it,"a_it.csv")
				call mati2csv(a_it_int,"a_it_int.csv")
				call mati2csv(status_it,"status_it.csv")
				call mati2csv(d_it,"d_it.csv")
			endif
			a_mean = 0._dp
			d_mean = 0._dp
			s_mean = 0._dp
			a_var = 0._dp
			d_var = 0._dp
			do age_hr = 1,TT-1
				junk = 1.
				do i=1,Nsim
					do it = ((init0_yrs)*itlen),((init_yrs+init0_yrs)*itlen)
						if( age_hr .eq. age_it(i,it) ) then
							a_mean(age_hr) = a_it(i,it)      + a_mean(age_hr)
							d_mean(age_hr) = d_it(i,it)      + d_mean(age_hr)
							s_mean(age_hr) = status_it(i,it) + s_mean(age_hr)
							junk = junk + 1._dp
						endif
					enddo

				enddo
				a_mean(age_hr) = a_mean(age_hr)/junk
				d_mean(age_hr) = d_mean(age_hr)/junk
				s_mean(age_hr) = s_mean(age_hr)/junk
				do i=1,Nsim
					do it = ((init0_yrs)*itlen),((init_yrs+init0_yrs)*itlen)
						if( (age_hr .eq. age_it(i,it)) .and. (a_it(i,it)>0._dp)) then
							a_var(age_hr) = (dlog(a_it(i,it)) - dlog(a_mean(age_hr)))**2 + a_var(age_hr)
							d_var(age_hr) = (d_it(i,it) - d_mean(age_hr))**2+ d_var(age_hr)
						endif
					enddo
				enddo
				a_var(age_hr) = a_var(age_hr)/junk
				d_var(age_hr) = d_var(age_hr)/junk
			enddo

			if( final_iter .eqv. .true. ) then
				exit !leave the iter loop!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			endif
			slice_len = iter
			simiter_dist(iter) = sum(dabs(a_mean - a_mean_liter) + dabs(s_mean - s_mean_liter)/5 )/TT
			simiter_status_dist(iter) = sum(dabs(s_mean - s_mean_liter))

			if (iter .ge. (iter_draws-1)) then
				if(verbose >=2 ) then
					print *, "done simulating on iteration ", iter, " of ", iter_draws
					print *, "dif a mean, log a var ",  sum(dabs(a_mean - a_mean_liter)), sum(dabs(a_var - a_var_liter))
					! NOTE: this is not actually mean because it does not have demographic weights
					print *, "a mean, a var ",  sum(a_mean), sum(a_var)
					print *, "dif status mean ",  sum(dabs(s_mean - s_mean_liter)/5)
					print *, "status mean ", sum(s_mean)
					print *,  "-------------------------------"
				endif

				final_iter = .true.

				if( w_strchng_old .eqv. .true. ) then
					w_strchng = .true.
				endif
			else
				if(verbose >=4 ) then
					print *, "iter:", iter
					print *, "dif a mean, log a var",  sum(dabs(a_mean - a_mean_liter)), sum((a_var - a_var_liter)**2)
					print *, "a mean, a var",  sum(a_mean), sum(a_var)
					print *, "dif status mean",  sum(dabs(s_mean - s_mean_liter))
					print *,  "-------------------------------"
					slice_len = max(1,slice_len)
					call vec2csv(simiter_dist(1:slice_len), "simiter_dist.csv" )
				endif
			endif
			a_mean_liter = a_mean
			d_mean_liter = d_mean
			s_mean_liter = s_mean
			a_var_liter = a_var
			d_var_liter = d_var
		enddo! iter
		if( print_lev>=2) then
			slice_len = max(1,slice_len)
			call vec2csv(simiter_dist(1:slice_len), "simiter_dist.csv" )
			call vec2csv(simiter_status_dist(1:slice_len), "simiter_status_dist.csv" )
		endif

		! calc occupation growth rates
		if(occaggs_hr) then
			if(verbose >2) print *, "calculating occupation growth rates"
			it = 1
			Nworkt = 0. !labor force in first period are the ones "born" in the first period
			occsize_jt = 0.
			occgrow_jt = 0.
			occshrink_jt = 0.
			do i=1,Nsim
				if( (age_it(i,it) > 0) .and. (status_it(i,it) <= 2) .and. (status_it(i,it)>=1)) Nworkt = 1._dp + Nworkt
				do ij=1,nj
					if((j_i(i) == ij) .and. (age_it(i,it) >= 1) .and. (status_it(i,it) == 1)) &
						& occsize_jt(ij,it) = 1._dp+occsize_jt(ij,it)
				enddo
			enddo
			occsize_jt(:,it) = occsize_jt(:,it) / Nworkt

			!$omp parallel do private(it,Nworkt,i,ij,occgrow_hr,occsize_hr,occshrink_hr)
			do it = 2,Tsim
				Nworkt = 0.
				occgrow_hr   = 0.
				occshrink_hr = 0.
				occsize_hr   = 0.
				do i=1,Nsim
					if(((status_it(i,it) <= 2) .and. (status_it(i,it)>=1)) .and. (age_it(i,it) > 0) ) Nworkt = 1._dp + Nworkt !labor force in this period
					ij =j_i(i)
					if( born_it(i,it) == 1 ) &
						& occgrow_hr(ij) = 1._dp + occgrow_hr(ij)
					if( (status_it(i,it-1) > 1)  .and. (status_it(i,it)== 1)) &
						& occgrow_hr(ij) = 1._dp + occgrow_hr(ij) !wasn't working last period
					if((status_it(i,it-1) == 1) .and. (status_it(i,it) > 1)) &
						& occshrink_hr(ij) = 1._dp + occshrink_hr(ij)
					if(status_it(i,it) == 1) &
						& occsize_hr(ij) = 1._dp + occsize_hr(ij)

				enddo
				do ij =1,nj
					if( occsize_hr(ij) > 0._dp ) then
						occgrow_jt(ij,it) = occgrow_hr(ij)/occsize_hr(ij)
						occshrink_jt(ij,it) = occshrink_hr(ij)/occsize_hr(ij)
					else
						occgrow_jt(ij,it) = 0._dp
						occshrink_jt(ij,it) = 0._dp
					endif
				enddo
				occsize_jt(:,it) = occsize_hr/Nworkt
			enddo
			!$omp end parallel do
		endif

		if(print_lev > 1) then
				if( caselabel == "") then
					if(welfare_cf .eqv. .true.)then
						call mat2csv (hst%welfare_hist,"welfare_hist.csv")
					endif
					call mat2csv (ewt_it,"ewt_it_hist.csv")
                    call mat2csv (trindif_it, "trindif_hist.csv")
					call mati2csv(e_it_int,"e_int_it_hist"//trim(caselabel)//".csv")
					call mati2csv(a_it_int,"a_int_it_hist"//trim(caselabel)//".csv")
					call mati2csv(brn_drawi_drawt(:,:,1),"brn_drawi_drawt.csv",0)
					call mati2csv(brn_drawi_drawt(:,:,2),"brn_drawi_drawt.csv",1)
					call mat2csv (occsize_jt,"occsize_jt_hist"//trim(caselabel)//".csv")
					call mat2csv (occgrow_jt,"occgrow_jt_hist"//trim(caselabel)//".csv")
					call mat2csv (occshrink_jt,"occshrink_jt_hist"//trim(caselabel)//".csv")
					call vec2csv(wtr_tmean_ts, "wtr_tmean_ts"//trim(caselabel)//".csv")
					! call vec2csv(acoef,"acoef.csv")
					! call vec2csv(ecoef,"ecoef.csv")
					! call mat2csv(Xdets(1:nregobs,:),"Xdets.csv")
				endif
				call mat2csv (a_it,"a_it_hist"//trim(caselabel)//".csv")
				call mat2csv (a_ap_dist,"a_ap_dist"//trim(caselabel)//".csv")
				call mat2csv (e_it,"e_it_hist"//trim(caselabel)//".csv")
				call mat2csv (shk%wtr_it,"wtr_it_hist"//trim(caselabel)//".csv")
				call mat2csv (trX_it,"transfer_it_hist"//trim(caselabel)//".csv")
				call mati2csv(al_int_it_endog,"al_int_endog_hist"//trim(caselabel)//".csv")
				call mat2csv (al_it_endog,"al_endog_hist"//trim(caselabel)//".csv")
				call mati2csv(status_it,"status_it_hist"//trim(caselabel)//".csv")
				call mati2csv(d_it,"d_it_hist"//trim(caselabel)//".csv")
				call mati2csv(shk%age_hist,"age_it_hist"//trim(caselabel)//".csv")
				call veci2csv(j_i,"j_i_hist"//trim(caselabel)//".csv")
				call veci2csv(z_jt_macroint,"z_jt_hist"//trim(caselabel)//".csv")
				call mat2csv (hst%wage_hist,"wage_it_hist"//trim(caselabel)//".csv")
				call mat2csv (hst%app_dif_hist,"app_dif_it_hist"//trim(caselabel)//".csv")
				call mat2csv (hst%di_prob_hist,"di_prob_it_hist"//trim(caselabel)//".csv")
				call mat2csv (hst%work_dif_hist,"work_dif_it_hist"//trim(caselabel)//".csv")
				call mati2csv(hst%hlth_voc_hist,"hlth_voc_hist"//trim(caselabel)//".csv")
				call mat2csv (hst%hlthprob_hist,"hlthprob_hist"//trim(caselabel)//".csv")


		endif

		deallocate(al_int_it_endog,al_it_endog)
		deallocate(e_it)
		deallocate(ewt_it)
        deallocate(trindif_it)
		deallocate(a_it_int,e_it_int)
		deallocate(app_it,work_it)
		deallocate(tr_hist_vec)

		deallocate(trX_it)
		deallocate(brn_drawi_drawt)
		deallocate(Xdets)
		deallocate(xidets,xjdets)
		deallocate(adep,edep)
		deallocate(acoef,ecoef)
		deallocate(invol_it)
		deallocate(cov_coef)

	end subroutine sim

end module sim_hists


!**************************************************************************************************************!
!**************************************************************************************************************!
!**************************************************************************************************************!
!**************************************************************************************************************!
!						Searches for parameters						       !
!**************************************************************************************************************!
!**************************************************************************************************************!
!**************************************************************************************************************!
!**************************************************************************************************************!

module find_params

	use V0para
	use helper_funs
	use sol_val
	use sim_hists
	use model_data
	!use mpi

	implicit none

	save

	integer :: mod_solcoefiter !just a counter for iterations on coefs
	type(val_struct), pointer :: mod_vfs
	type(pol_struct), pointer :: mod_pfs
	type(hist_struct), pointer :: mod_hst
	type(shocks_struct), pointer :: mod_shk
	!integer :: dfbols_nuxi_trproc = 1 !indicator for whether dfovec should solve for nuxi or for the zprocess
	real(dp) :: mod_prob_target
	real(dp), allocatable :: mod_xl(:),mod_xu(:)

	contains

	subroutine reg_wgtrend( coefs,vfs, pfs, hst,shk )
	!simulates and pumps out wage coefficients

		type(shocks_struct) :: shk
		type(val_struct) :: vfs
		type(pol_struct) :: pfs
		type(hist_struct):: hst
		real(dp), intent(out) :: coefs(:)

		integer :: Ncoef,status, i, it,ri,ii,ip,ik,ij

		real(dp), allocatable :: XX(:,:), yy(:), coef_est(:), cov_coef(:,:), XX_ii(:,:), yy_ii(:)
		real(dp) :: hatsig2
		real(dp) :: tbase(Tsim,Nknots-1)

		if(tr_spline .eqv. .false.) then
			if( NKpolyT>=2 ) then
				Ncoef = (Nskill+1)*(NKpolyT+1)+2 !NKpolyT*Nskill + Nskill + NKpolyT + 2 + const
			else
				Ncoef = Nskill*2 + NTpolyT+Nnuisance !Nskill-level, Nskill-time trend, cubic time, 2 age, const
			endif
		else
			Ncoef = Nskill*(Nknots-1) + Nskill + (Nknots-1)+Nnuisance
		endif

		allocate(XX(Tsim*Nsim,Ncoef))
		allocate(yy(Tsim*Nsim))
		allocate(coef_est(Ncoef))
		allocate(cov_coef(Ncoef,Ncoef))

		if(Ncoef .ne. size(coefs)) then
			print *, "coefs passed not the same as expected number of coefficients"
		endif

		call sim(vfs, pfs, hst,shk,.false.)

		tbase = 0._dp
		do it=(TossYears*itlen+1),Tsim
			tbase(it,1) = (dble(it)/tlen - dble(TossYears))
			do ip=1,(Nknots-2)
				if((dble(it)/tlen - dble(TossYears) - tr_knots(ip)) > 0._dp) &
				& 	tbase(it,ip+1) = (dble(it)/tlen - dble(TossYears) - tr_knots(ip))**3 + tbase(it,ip+1)
				if( dble(it)/tlen - dble(TossYears) - tr_knots(Nknots-1) >0._dp) &
				&	tbase(it,ip+1) = -(dble(it)/tlen - dble(TossYears) - tr_knots(Nknots-1))**3 *(tr_knots(Nknots)-tr_knots(ip))/(tr_knots(Nknots)-tr_knots(Nknots-1)) &
					&  + tbase(it,ip+1)
				if( dble(it)/tlen - dble(TossYears) - tr_knots(Nknots) >0._dp ) &
				& 	tbase(it,ip+1) = -(dble(it)/tlen - dble(TossYears) - tr_knots(Nknots) )**3 *(tr_knots(Nknots-1)-tr_knots(ip))/(tr_knots(Nknots)-tr_knots(Nknots-1)) &
					&  + tbase(it,ip+1)
				tbase(it,ip+1) = tbase(it,ip+1)*(tr_knots(Nknots)-tr_knots(1))**(-2)
			enddo
		enddo
		ii = 0
		do i=1,Nsim
			do it=(TossYears*itlen+1),Tsim
				if(  (shk%age_hist(i,it) > 0) .and. (hst%status_hist(i,it)==1) ) then
				!& .and. (dble(it)/tlen-dble(TossYears) .ge. tr_knots(1)) .and. (dble(it)/tlen-dble(TossYears) .le. tr_knots(Nknots)) ) then
					ii = 1+ii
					ij = shk%j_i(i) !this guy's occupation
					yy(ii) = log(hst%wage_hist(i,it))
					ri=1

					if(tr_spline .eqv. .false.) then
						if(NKpolyT >= 2 ) then
							do ip=1,(NKpolyT+1)
								do ik=1,(Nskill+1)
									if((ik == 1) .and. (ip == 1)) then
										XX(ii,ri) = 1._dp
									elseif( (ik==1) .and. (ip>1)) then
										XX(ii,ri) = (dble(it)/tlen-dble(TossYears))**(ip-1)
									else
										XX(ii,ri) = occ_onet(ij,ik-1)*(dble(it)/tlen-TossYears)**(ip-1)
									endif
									ri = ri+1
								enddo !ik, skill
							enddo !ip, poly degree
						else
							do ip=1,NTpolyT
								XX(ii,ri) = (dble(it)/tlen-dble(TossYears))**ip
								ri = ri+1
							enddo
							do ik=1,Nskill
								XX(ii,ri) = occ_onet(ij,ik)
								ri = ri+1
							enddo
							do ik=1,Nskill
								XX(ii,ri) = (dble(it)/tlen-dble(TossYears))*occ_onet(ij,ik)
								ri = ri+1
							enddo
							XX(ii,ri) = 1._dp
							ri = ri+1
						endif
					else
						do ik=1,Nskill
							XX(ii,ri) = occ_onet(ij,ik)
							ri = ri+1
						enddo
						do ip=1,(Nknots-1)
							XX(ii,ri) = tbase(it,ip)
							ri = ri+1
						enddo
						do ik=1,Nskill
							do ip=1,(Nknots-1)
								XX(ii,ri) = tbase(it,ip)*occ_onet(ij,ik)
								ri = ri+1
							enddo
						enddo
						XX(ii,ri) = 1._dp
						ri = ri+1
					endif
					XX(ii,ri) = agegrid( shk%age_hist(i,it) )
					ri=ri+1
					XX(ii,ri) = agegrid( shk%age_hist(i,it) )** 2
					!add in health dummies
					do ip=2,nd
						ri = ri+1
						if( shk%d_hist(i,it) == ip ) then
							XX(ii,ri) = 1._dp
						else
							XX(ii,ri) = 0._dp
						endif
					enddo
				endif ! participating
			enddo !it
		enddo ! i

		allocate(XX_ii(ii,ri))
		allocate(yy_ii(ii))
		do i=1,ii
			XX_ii(i,:) = XX(i,:)
			yy_ii(i) = yy(i)
		enddo

		call OLS(XX_ii,yy_ii,coef_est,cov_coef, hatsig2, status)

		do i=1,Ncoef
			coefs(i) = coef_est(i)
		enddo
		if( print_lev .ge. 2) then
			call mat2csv(XX_ii, "XX_ii.csv")
			call vec2csv(yy_ii, "yy_ii.csv")
			call vec2csv(coef_est, "coef_est.csv")
		endif


		deallocate(XX_ii,yy_ii)
		deallocate(XX,yy,coef_est,cov_coef)

	end subroutine reg_wgtrend

	subroutine gen_new_wgtrend(new_wgtrend, wage_coef_in)

		real(dp), intent(in)  :: wage_coef_in(:)
		real(dp), intent(out) :: new_wgtrend(:,:)
		integer :: ij, it, ip,i,ik
		real(dp) :: wage_trend_hr, tbase(Nknots-1), realwagegrowth(36)

		wage_lev = 0._dp
		if(tr_spline .eqv. .true.) then
			do it=1,Tsim
				tbase=0._dp
				if(it>TossYears*itlen) then
					tbase(1) = dble(it)/tlen-dble(TossYears)
					do ip=1,(Nknots-2)
						if((dble(it)/tlen - dble(TossYears) - tr_knots(ip)) > 0.) &
						& 	tbase(ip+1) = (dble(it)/tlen - dble(TossYears) - tr_knots(ip))**3 + tbase(ip+1)
						if( dble(it)/tlen - dble(TossYears) - tr_knots(Nknots-1) >0.) &
						&	tbase(ip+1) = -(dble(it)/tlen - dble(TossYears) - tr_knots(Nknots-1))**3 *(tr_knots(Nknots)-tr_knots(ip))/(tr_knots(Nknots)-tr_knots(Nknots-1)) &
							&  + tbase(ip+1)
						if( dble(it)/tlen - dble(TossYears) - tr_knots(Nknots) >0. ) &
						& 	tbase(ip+1) = -(dble(it)/tlen - dble(TossYears) - tr_knots(Nknots) )**3 *(tr_knots(Nknots-1)-tr_knots(ip))/(tr_knots(Nknots)-tr_knots(Nknots-1)) &
							&  + tbase(ip+1)
						tbase(ip+1) = tbase(ip+1)*(tr_knots(Nknots)-tr_knots(1))**(-2)
					enddo
				endif
				tbase_out(it,:) = tbase
				do ij=1,nj
					wage_trend_hr = 0._dp
					do ik=1,Nskill
						wage_trend_hr = wage_coef_in(ik)*occ_onet(ij,ik) +wage_trend_hr
					enddo
					if(nomeantrend .eqv. .false.) then
						do ip=1,(Nknots-1)
							wage_trend_hr = wage_coef_in(ip+Nskill)*tbase(ip) + wage_trend_hr
						enddo
					endif
					do ik=1,Nskill
						do ip=1,(Nknots-1)
							wage_trend_hr = wage_coef_in(ip + (ik-1)*(Nknots-1)+Nskill+(Nknots-1))*tbase(ip)*occ_onet(ij,ik) +wage_trend_hr
						enddo
					enddo !for each ik
					if(wage_trend_hr .gt. trgrid(ntr)) wage_trend_hr = trgrid(ntr)
					if(wage_trend_hr .lt. trgrid(1)  ) wage_trend_hr = trgrid(1)

					new_wgtrend(it,ij) = wage_trend_hr
				enddo !for ij=1:nj
			enddo !for it=1:Tsim
			do i=1,nj
				if( wglev_0 .eqv. .false.) wage_lev(i) = new_wgtrend(TossYears*itlen,i)
				new_wgtrend(:,i) = new_wgtrend(:,i) - new_wgtrend(TossYears*itlen,i)
			enddo
		else
			if( NKpolyT >= 2 ) then
				do ij=1,nj
					do it=1,Tsim
						wage_trend_hr = 0._dp
						if(it>TossYears*itlen) then
							do ip =1,(NKpolyT+1)
								if((wglev_0 .eqv. .false.) .or. (ip .gt. 1)) then
									if(ip .gt. 1) &
									&	wage_trend_hr  = (dble(it)/tlen-dble(TossYears))**(ip-1)*wage_coef_in( (ip-1)*(Nskill+1)+1 )                   + wage_trend_hr
									do ik=2,(Nskill+1)
										wage_trend_hr  = (dble(it)/tlen-dble(TossYears))**(ip-1)*wage_coef_in( (ip-1)*(Nskill+1)+ik)*occ_onet(ij,ik-1) + wage_trend_hr
									enddo
								endif
							enddo
						else
							do ik=2,(Nskill+1)
								wage_trend_hr  = (0._dp)**(ip-1)*wage_coef_in( (ip-1)*(Nskill+1)+ik)*occ_onet(ij,ik-1) + wage_trend_hr
							enddo
						endif
						if(wage_trend_hr .gt. trgrid(ntr)) wage_trend_hr = trgrid(ntr)
						if(wage_trend_hr .lt. trgrid(1)  ) wage_trend_hr = trgrid(1)

						new_wgtrend(it,ij) = wage_trend_hr !upd_wgtrnd*wage_trend_hr + (1._dp - upd_wgtrnd)*wage_trend(it,ij)
					enddo
				enddo
				do i=1,nj
					if( wglev_0 .eqv. .false.) wage_lev(i) = new_wgtrend(1,i)
					new_wgtrend(:,i) = new_wgtrend(:,i) - new_wgtrend(1,i)
				enddo
			! LINEAR in time for ONET skills
			else
				do ij=1,nj
				do it=1,Tsim
					wage_trend_hr = 0._dp
					if( it>TossYears*itlen) then
						do ip =1,NTpolyT
							wage_trend_hr = (dble(it)/tlen - dble(TossYears))**ip*wage_coef_in(ip) + wage_trend_hr
						enddo
						do ik=1,Nskill
							wage_trend_hr = (dble(it)/tlen - dble(TossYears))*wage_coef_in(ik+NTpolyT+ Nskill)*occ_onet(ij,ik) &
							& + wage_coef_in(ik+NTpolyT)*occ_onet(ij,ik) + wage_trend_hr
						enddo
					else
						do ip =1,NTpolyT
							wage_trend_hr = (0._dp)**ip*wage_coef_in(ip) + wage_trend_hr
						enddo
						do ik=1,Nskill
							wage_trend_hr = (0._dp)*wage_coef_in(ik+NTpolyT+ Nskill)*occ_onet(ij,ik) &
							& + wage_coef_in(ik+NTpolyT)*occ_onet(ij,ik) + wage_trend_hr
						enddo
					endif
					new_wgtrend(it,ij) = wage_trend_hr
				enddo
				enddo
				do ij=1,nj
					if( wglev_0 .eqv. .false.) wage_lev(ij) = new_wgtrend(TossYears*itlen,ij)
					new_wgtrend(:,ij) = new_wgtrend(:,ij) - new_wgtrend(TossYears*itlen,ij)
					do it=1,Tsim
						if(new_wgtrend(it,ij) .gt. trgrid(ntr)) new_wgtrend(it,ij) = trgrid(ntr)
						if(new_wgtrend(it,ij) .lt. trgrid(1)  ) new_wgtrend(it,ij) = trgrid(1)
					enddo
				enddo
			endif !linear in onet skills
		endif !cubic spline
		
		!deflate by real wage growth (series LEU0252881600A_PC1, median wage of hourly and salaried workers from the CPS)
		!realwagegrowth = (/ 0., 1.91083, 2.1875, 0.61162, -0.91185, -1.22699, -2.17391, -0.63492, 0.31949, 0.95541, -0.63091, -0.31746, -0.31847, 0.31949, 2.2293, 2.49221, 1.51976, 0.5988, 0.59524, -0.29586, 0.29674, -1.47929, 0., 0.6006, 0., 2.98507, -0.86957, -1.75439, -0.29762, -0.59701, 0.3003, 2.09581, 1.75953, 1.15274, 0.5698, 1.69972  /)
		!realwagegrowth(1) = exp(realwagegrowth(1))
		!do it=2,size(realwagegrowth)
		!	realwagegrowth(it) = exp(realwagegrowth(it)/100.)*realwagegrowth(it-1)
		!enddo
		
		!! FOR COUNTER FACTUALS !! ++++++++++++++++++
		!! Replace wage trend for everyone with wtr_tsmean.
		if( wtr_by_occ .eqv. .false. ) then
			do ij=1,nj
				do it=1,Tsim
					if(it>TossYears*itlen) then
						new_wgtrend(it,ij) = wtr_tmean_ts(it)
					endif
				enddo
			enddo
		endif

		call mat2csv(tbase_out, "tbase_main.csv")
		if(print_lev .ge. 2) then
			call mat2csv(new_wgtrend,"new_wgtrend.csv")
			call vec2csv(wage_lev,"new_wglev.csv")
		endif

	end subroutine gen_new_wgtrend

	subroutine set_glbparams(paramvec)
		real(dp), intent(in) :: paramvec(:)
		integer :: npar

		npar = size(paramvec)
! Parameter vector
!         nud1,   nud2,     nud3 , nuage  ,  xizcoef,   xiage    , xiVtrend , xid1coef  ,  xid2coef ,  xid3coef ,    Fd(2,2) ,       Fd(3,2)

		if( nu_byD .eqv. .false. ) then
			nu   = paramvec(1)
			xizcoef = paramvec(2)
			xiagezcoef = paramvec(3)
			Fd(2,1) = paramvec(4)
			Fd(3,1) = paramvec(5)

			Fd(2,2) = paramvec(4)
			Fd(3,2) = paramvec(5)
		else
			nud(1)      = paramvec(1)
			nud(2)      = paramvec(2)
			nud(3)      = paramvec(3)
			nuage 		= paramvec(4)
			xizcoef 	= paramvec(5)
			xiagezcoef  = paramvec(6)
			xiVtrend    = paramvec(7)
			xid1coef    = paramvec(8)
			xid2coef    = paramvec(9)
			xid3coef    = paramvec(10)

			Fd(2,1) 	= paramvec(11)
			Fd(3,1) 	= paramvec(12)


			Fd(2,2) = paramvec(11)
			Fd(3,2) = paramvec(12)
		endif

		if(verbose >2) &
		& 	print *, "setting ", npar , " parameters"

	end subroutine set_glbparams

	subroutine readwrite_glbparams(paramvecw, fname, readflag)
		real(dp) :: paramvecw(:) !parameter vector and wmean
		logical, intent(in) :: readflag ! =true if read, =false if write
		character(LEN=*), intent(in) :: fname

		CHARACTER(LEN=*), PARAMETER  :: FMT = "(G20.12)"
		CHARACTER(LEN=20) :: FMT_1

		integer :: i , ui

		ui = fread + 11

		if(readflag .eqv. .true.) then
			open(unit = ui, file=fname, action='read')
			do i=1,nopt_pars
				read(ui,*) paramvecw(i)
			enddo
			read(ui ,*) wmean
			close(ui)

		else
			call vec2csv( paramvecw  , fname)
		endif

	end subroutine readwrite_glbparams


	subroutine cal_dist(paramvec, errvec,shk)
		! the inputs are the values of parameters we're moving in paramvec
		! the outputs in errvec are deviations from targets

		real(dp), intent(in) :: paramvec(:)
		real(dp), intent(out) :: errvec(:)

		type(shocks_struct) :: shk
		type(val_struct) :: vfs,vfs_cf
		type(pol_struct) :: pfs,pfs_cf
		type(hist_struct):: hst
		type(moments_struct):: moments_sim
		real(dp) :: condstd_tsemp,totdi_rt,totapp_dif_hist,ninsur_app,napp_t,nu1,nu0, objval

		real(dp) :: fndgrid0(nl,nz),sepgrid0(nl,nz)
		real(dp) :: errwt(nopt_tgts)
		real(dp) :: dist_urt,dist_frt,Efrt,Esrt,urt !moments fro the unemployment rate calibration
		integer :: ij=1,t0tT(2),it,i,j,npar, nerr
		integer :: rank_hr,ierr,fcal_eval
		character(2) :: rank_str
		npar = size(paramvec)
		nerr = size(errvec)
		cal_niter = cal_niter + 1

		do i=1,nopt_tgts
			errwt(i) =1.
		enddo
		errwt(1) =2._dp**0.5_dp
		errwt(2) =2._dp**0.5_dp

		call set_glbparams(paramvec)

!		call mpi_comm_rank(mpi_comm_world,rank_hr,ierr)
		call alloc_econ(vfs,pfs,hst)

		! set up economy and solve it
		call set_zjt(hst%z_jt_macroint, hst%z_jt_panel, shk) ! includes call settfp() which resets fndgrid and sepgrid

		! these must come after set_zjt because that also resets fndgrid, sepgrid:
		forall(i=1:nl,j=1:nz) fndgrid0(i,j) = fndgrid(i,j)
		forall(i=1:nl,j=1:nz) sepgrid0(i,j) = sepgrid(i,j)
		forall(i=1:nl,j=1:nz) fndgrid(i,j) = fndrt_mul * fndgrid0(i,j)
		forall(i=1:nl,j=1:nz) sepgrid(i,j) = seprt_mul * sepgrid0(i,j)


		if(print_lev .ge. 2) call mat2csv(fndgrid0,"fndgrid0"//trim(caselabel)//".csv")
		if(verbose >2) print *, "In the calibration"


		!solve w/o DI, in case of welfare_cf .eqv. .true.
		if(welfare_cf .eqv. .true. ) then
			call alloc_valpol(vfs_cf,pfs_cf)
			noDI = .true.

			call sol(vfs_cf,pfs_cf)
			vfs%V_CF = vfs_cf%V
			noDI = .false.
			call dealloc_valpol(vfs_cf,pfs_cf)
		endif

		if(verbose >2) print *, "Solving the model"
		call sol(vfs,pfs)

		if(verbose >2) print *, "Simulating the model"
		call sim(vfs, pfs, hst,shk,.false.)
		if(verbose >2) print *, "Computing moments"
		call moments_compute(hst,moments_sim,shk)
		if(verbose >1) print *, "DI rate" , moments_sim%avg_di

		!unemployment stats
		if(verbose>1 .or. print_lev > 1) then
			call comp_ustats(hst,shk,urt,Efrt,Esrt)
			dist_urt = (urt - avg_unrt)/avg_unrt
			dist_frt = (Efrt - avg_frt)/avg_frt
			print * , "missing on fidning rate ", dist_frt, " missing on u rate ", dist_urt

		endif

		totapp_dif_hist = 0.
		ninsur_app = 0.
		do i =1,Nsim
			do it=1,Tsim
				if((hst%status_hist(i,it)<=3) .and. (mod(it,itlen) .eq. 0)) ninsur_app = 1.+ninsur_app ! only count the body once every year, comparable to data
				if(hst%status_hist(i,it)==3) &
				&	totapp_dif_hist = exp(10.*hst%app_dif_hist(i,it))/(1. + exp(10.*hst%app_dif_hist(i,it))) + totapp_dif_hist

			enddo
		enddo


		totapp_dif_hist = totapp_dif_hist/ninsur_app
		if(verbose >1) print *, "App rate (smooth)" , totapp_dif_hist

		if(verbose >1) print *, "age effect here: ", moments_sim%diaward_ageeffect
		errvec(1)  = errwt(1)*(moments_sim%init_diaward - diaward_target)/diaward_target
		errvec(2)  = errwt(2)*(moments_sim%init_hlth_acc - hlth_acc_target)/hlth_acc_target
		errvec(3)  = errwt(3)*(moments_sim%old_initdiawdfrac - old_target)/old_target
		errvec(4)  = errwt(4)*(moments_sim%avg_diprob - allowrt_target)/allowrt_target
		errvec(5)  = errwt(5)*(moments_sim%diaward_ageeffect - award_age_target)/award_age_target

		errvec(6)  = errwt(6)*(moments_sim%work_rateD(1)/moments_sim%work_rateDage(1,1) - p1d1_target  )*2/ ( moments_sim%work_rateD(1)/moments_sim%work_rateDage(1,1) + p1d1_target )
		errvec(7)  = errwt(7)*(moments_sim%work_rateD(2)/moments_sim%work_rateDage(1,1) - p1d2_target  )*2/ ( moments_sim%work_rateD(2)/moments_sim%work_rateDage(1,1) + p1d2_target )
		errvec(8)  = errwt(8)*(moments_sim%work_rateD(3)/moments_sim%work_rateDage(1,1) - p1d3_target  )*2/ ( moments_sim%work_rateD(3)/moments_sim%work_rateDage(1,1) + p1d3_target )

		errvec(9)  = errwt(9)*(moments_sim%d1_diawardfrac - d1_diawardfrac_target)/d1_diawardfrac_target
		errvec(10) =errwt(10)*(moments_sim%d2_diawardfrac - d2_diawardfrac_target)/d2_diawardfrac_target
		errvec(11) =errwt(11)*(moments_sim%d3_diawardfrac - d3_diawardfrac_target)/d3_diawardfrac_target

		errvec(12) =errwt(12)*(moments_sim%end_hlth_acc - hlth_acc_end_target)/hlth_acc_end_target

		!errvec(5) = (moments_sim%diaward_ageeffect - award_age_target)/award_age_target

		!errvec(4) = (moments_sim%work_rateDage(1,1)/moments_sim%work_rateDage(2,1) - 1._dp/p1d2_2545target)*2/(moments_sim%work_rateDage(1,1)/moments_sim%work_rateDage(2,1) + 1._dp/p1d2_2545target)
		!errvec(5) = (moments_sim%work_rateDage(1,1)/moments_sim%work_rateDage(3,1) - 1._dp/p1d3_2545target)*2/(moments_sim%work_rateDage(1,1)/moments_sim%work_rateDage(3,1) + 1._dp/p1d3_2545target)
		!errvec(6) = (moments_sim%work_rateDage(1,1)/moments_sim%work_rateDage(1,2) - 1._dp/p1d1_4665target)*2/(moments_sim%work_rateDage(1,1)/moments_sim%work_rateDage(1,2) + 1._dp/p1d1_4665target)
		!errvec(7) = (moments_sim%work_rateDage(1,1)/moments_sim%work_rateDage(2,2) - 1._dp/p1d2_4665target)*2/(moments_sim%work_rateDage(1,1)/moments_sim%work_rateDage(2,2) + 1._dp/p1d2_4665target)
		!errvec(8) = (moments_sim%work_rateDage(1,1)/moments_sim%work_rateDage(3,2) - 1._dp/p1d3_4665target)*2/(moments_sim%work_rateDage(1,1)/moments_sim%work_rateDage(3,2) + 1._dp/p1d3_4665target)


		objval =0._dp
		do i=1,nerr
			objval =errvec(i)**2 +objval
		enddo

		write(rank_str, '(I2.2)') rank_hr
		fcal_eval = rank_hr+100
		open(unit=fcal_eval, file = "cal_dist_"//trim(rank_str)//".csv" ,ACCESS='APPEND', POSITION='APPEND')
		write(fcal_eval,  "(G20.12)", advance='no') objval
		write(fcal_eval, "(A)", advance='no')  " , "
		do i=1,size(paramvec)
			write(fcal_eval, "(G20.12)", advance='no')  paramvec(i)
			write(fcal_eval, "(A)", advance='no')   " , "
		enddo
		do i=1,size(errvec)
			write(fcal_eval, "(G20.12)", advance='no') errvec(i)
			write(fcal_eval, "(A)", advance='no') " , "
		enddo
		!some diagnostics to print
		if( moments_sim%init_diaward_discr >0._dp ) then
			write(fcal_eval, "(G20.12)", advance='no')  moments_sim%init_diaward/moments_sim%init_diaward_discr
			write(fcal_eval, "(A)", advance='no') " , "
		else
			write(fcal_eval, "(G20.12)", advance='no')  0._dp
			write(fcal_eval, "(A)", advance='no') " , "

		endif
		write(fcal_eval, "(G20.12)", advance='yes')  moments_sim%avg_di

		close(unit=fcal_eval)
		cal_obj = sum( errvec**2 )

		!put things back:
		call dealloc_econ(vfs,pfs,hst)

		fndgrid = fndgrid0
		sepgrid = sepgrid0

	end subroutine cal_dist

	subroutine cal_mlsl( fval, xopt, nopt, xl, xu, shk)

		real(dp), intent(out) :: fval
		real(dp), intent(out) :: xopt(:)
		real(dp), intent(in)  :: xl(:),xu(:)
		integer , intent(in)  :: nopt
		type(shocks_struct), target :: shk

		integer  :: ndraw, nstartpn, ninterppt,ninternalopt,d,dd,i,j,ii
		integer :: ierr, iprint, rank, nnode, nstarts
		integer, allocatable :: seedhr(:)
		real(dp) :: draw(nopt)
		real(dp) :: x0(nopt),x0w(nopt+1), x0hist(nopt,500),xopt_hist(nopt,500),fopt_hist(500),v_err(nopt_tgts)
		real(dp) :: rhobeg, rhoend, EW,W,err0,fdist,xdist,junk,attraction_size
		real(dp), allocatable :: wspace(:),sepgrid0(:,:),fndgrid0(:,:)
		real(dp), allocatable :: internalopt_hist(:,:)
		real(dp), allocatable :: node_fopts(:),node_xopts(:),world_fopts(:), world_xopts(:)
		real(dp), allocatable :: world_internalopt(:),node_internalopt(:)
		character(len=2) :: rank_str
		integer :: fcal_eval,print_lev_old,verbose_old,seed_size
		integer :: c1=1,c2=1,cr=1,cm=1
		real(dp) :: t1=1.,t2=1.
		real(dp) :: zeros(nopt),ones(nopt), initerr



		external dfovec
		call random_seed(size = seed_size)
		allocate( seedhr(seed_size) )
		allocate(sepgrid0(size(sepgrid,1),size(sepgrid,2)))
		allocate(fndgrid0(size(fndgrid,1),size(fndgrid,2)))

		allocate(mod_xl(nopt))
		allocate(mod_xu(nopt))
		mod_xl = xl ! this is so bounds are available in dfovec
		mod_xu = xu
		zeros  = 0._dp
		ones   = 1._dp

		ndraw = size(x0hist,2)
		ninterppt = 2*nopt+1
		ninternalopt = size(wage_coef) + 2 !wage coefficients, fmul smul
		nstartpn = 1 !if it takes lots of starts, probably can speed things up by statically allocating more starts per node

		fopt_hist = -1.e5_dp

		mod_shk => shk

		print_lev_old = print_lev
		verbose_old = verbose
		print_lev = 0
		verbose = 1

!		call mpi_comm_size(mpi_comm_world,nnode,ierr)
!		call mpi_comm_rank(mpi_comm_world,rank,ierr)
		nnode = 1
		rank = 0
		print *, "Starting on node ", rank, "out of ", nnode
		!set ndraw to do it only once:
		ndraw = nnode*nstartpn
		allocate( internalopt_hist(size(wage_coef)+2,500) )

		allocate(wspace((ninterppt+5)*(ninterppt+nopt)+3*nopt*(nopt+5)/2))
		allocate(node_fopts(               nstartpn))
		allocate(node_xopts(               nstartpn*nopt))
		allocate(node_internalopt(         nstartpn*ninternalopt))

		allocate(world_fopts(        nnode*nstartpn))
		allocate(world_xopts(        nnode*nstartpn*nopt))
		allocate(world_internalopt(  nnode*nstartpn*ninternalopt))
		!dfbols_nuxi_trproc = 1

		rhobeg = 0.4_dp/sqrt(dble(nnode)) !trust region starting size is almost the entire space (80% of it)
		rhoend = rhobeg/1.e4_dp

		fval = 0._dp
		xopt = 0._dp
		attraction_size = 0._dp
		call insobl( nopt, ndraw )
		do d = 1,ndraw
			call I4_SOBOL( nopt, draw )
			!x0hist(:,d) = draw
			do dd=1,nopt
				x0hist(dd,d) = xl(dd)+draw(dd)*(xu(dd)-xl(dd))
				if((nnode <= 1) .and. (d==1)) then
					x0hist(dd,d) = xl(dd)+0.5_dp*(xu(dd)-xl(dd))
				endif
			enddo
		enddo

		if( rank .eq. 0 ) then
			call mat2csv(x0hist,"x0hist.csv")
			open(unit=fcallog, file = "cal_mlsl.csv")
			write( fcallog,* ) " "
			close(fcallog)
		endif
		write( rank_str,"(I2.2)" ) rank
		fcal_eval = 100+rank
		open(unit=fcal_eval, file = "cal_dist_"//trim(rank_str)//".csv" )
		write(fcal_eval, "(A)", advance='no') "obj , "
		write(fcal_eval, "(A)", advance='no') "nu1 , "
		write(fcal_eval, "(A)", advance='no') "nu2 , "
		write(fcal_eval, "(A)", advance='no') "nu3 , "
		write(fcal_eval, "(A)", advance='no') "nuage , "
		write(fcal_eval, "(A)", advance='no') "xiz , "
		write(fcal_eval, "(A)", advance='no') "xiage , "
		write(fcal_eval, "(A)", advance='no') "xiVtrend , "
		write(fcal_eval, "(A)", advance='no') "xid1 , "
		write(fcal_eval, "(A)", advance='no') "xid2 , "
		write(fcal_eval, "(A)", advance='no') "xid3 , "
		write(fcal_eval, "(A)", advance='no') "Fd2 , "
		write(fcal_eval, "(A)", advance='no') "Fd3 , "

		write(fcal_eval, "(A)", advance='no') "awd , "
		write(fcal_eval, "(A)", advance='no') "hlth , "
		write(fcal_eval, "(A)", advance='no') "old pct, "
		write(fcal_eval, "(A)", advance='no') "allowrt , "
		write(fcal_eval, "(A)", advance='no') "age effect, "
		write(fcal_eval, "(A)", advance='no') "workrt1 , "
		write(fcal_eval, "(A)", advance='no') "workrt2 , "
		write(fcal_eval, "(A)", advance='no') "workrt3 , "
		write(fcal_eval, "(A)", advance='no') "d1awdfrac , "
		write(fcal_eval, "(A)", advance='no') "d2awdfrac , "
		write(fcal_eval, "(A)", advance='no') "d3awdfrac , "
		write(fcal_eval, "(A)", advance='no') "end_hlth , "

		close( fcal_eval )
		fndgrid0 = fndgrid
		sepgrid0 = sepgrid

		forall(d=1:seed_size) seedhr(d) = rank+671984 + (d-1)*100
		call random_seed(put = seedhr(1:seed_size))
		!loop/distribute stating points for random restarts
		do d =1,(ndraw/nnode/nstartpn)
		!do d=1,1
			do dd= 1,nstartpn
				call system_clock(count_rate=cr)
				call system_clock(count_max=cm)
				call CPU_TIME(t1)
				call SYSTEM_CLOCK(c1)

				x0 = x0hist( :, (d-1)*(nnode*nstartpn) + rank*nstartpn + dd ) !note this indexing looks weird (normally it's rank-1) but rank is base 0

				!if x0 in a basin of attraction, cycle
				do j=1,9

					cal_niter = 0
					xopt = (x0-xl)/(xu-xl)
					call dfovec(nopt,nopt,xopt  ,v_err) !note that dfovec takes the normalized values between 0,1 for x
					initerr = sum( abs(v_err) )

					if( (v_err(1) .ge. 1.99_dp .or. v_err(1) .le. -1.99_dp) .and. nnode .gt. 1 ) then !throw out the ones with crazy DI rates
						! exit with another draw
						call random_number(junk) !burn one
						do i=1,nopt
							call random_number(x0(i))
							x0(i) = x0(i)*( xu(i)-xl(i) ) + xl(i)
						enddo
					else

						! only do BOBYQA if the starting point is good
						if(nopt >=8) then
							print *, "Computing from: ",  x0(1), x0(2), x0(3),x0(4), x0(5), x0(6), x0(7), x0(8)," on node: ", rank, "after ", j, " tries"
						elseif(nopt >=7) then
								print *, "Computing from: ",  x0(1), x0(2), x0(3),x0(4), x0(5), x0(6), x0(7)," on node: ", rank, "after ", j, " tries"
						elseif(nopt >=6) then
							print *, "Computing from: ",  x0(1), x0(2), x0(3),x0(4), x0(5), x0(6)," on node: ", rank, "after ", j, " tries"
						elseif(nopt >=5) then
							print *, "Computing from: ",  x0(1), x0(2), x0(3),x0(4), x0(5)," on node: ", rank, "after ", j, " tries"
						elseif(nopt >=4) then
							print *, "Computing from: ",  x0(1), x0(2), x0(3),x0(4)," on node: ", rank, "after ", j, " tries"
						else
							print *, "Computing from: ",  x0(1), x0(2), x0(3)," on node: ", rank, "after ", j, " tries"
						endif

						cal_niter = 1
						iprint = 1
						xopt = (x0-xl)/(xu-xl) !convert x0 draw into normalized (0,1) units
						call bobyqa_h(nopt,ninterppt,xopt,zeros,ones, &
						&	rhobeg,rhoend,iprint,max_DFBOLS,wspace,nopt)

						call dfovec(nopt,nopt,xopt,v_err)
						exit
					endif
				enddo !j=1,5 to loop over starting points
				attraction_size = sum(dabs(xopt - (x0-xl)/(xu-xl)))/dble(nopt)

				err0 = sum( v_err**2 )
				node_fopts(dd) = err0
				node_xopts(((dd-1)*nopt+1):(dd*nopt)) = xopt*(xu-xl)+xl
				node_internalopt( ((dd - 1)*ninternalopt+1):(dd - 1)*ninternalopt + size(wage_coef) ) = wage_coef
				node_internalopt( ((dd - 1)*ninternalopt+1+ size(wage_coef)):dd*ninternalopt ) = (/ fndrt_mul, seprt_mul /) !are these available globally?
				print *, "-----------------"
				print *, "-----------------"
				print *, "Found min ", err0, "at ", xopt," on node: ", rank, "attraction basin: ", attraction_size

				call CPU_TIME(t2)
				call SYSTEM_CLOCK(c2)
				print *, "System Time on calibration for rank ", rank, ": ", dble(c2-c1)/dble(cr)
				print *, "   CPU Time for rank ", rank_str, ": ", (t2-t1)
				print *, "-----------------"
				print *, "-----------------"

				call vec2csv(wage_coef, "wage_coef_opt"// rank_str //".csv")
				call mat2csv(fndgrid0*fndrt_mul, "fndgrid_opt"// rank_str //".csv")
				call mat2csv(sepgrid0*seprt_mul, "sepgrid_opt"// rank_str //".csv")

				x0 = xopt*(xu-xl)+xl !convert to input space
				call set_glbparams(x0)

				x0w(1:nopt) = x0
				x0w(nopt+1) = wmean
				call  readwrite_glbparams(x0w ,"nuxiw_opt" // rank_str //".csv", .false.)
			enddo ! dd = 1,nstartpn

			open(unit=fcallog, file = "cal_mlsl.csv" ,ACCESS='APPEND', POSITION='APPEND')
			do i=1,(nstartpn)
				write(fcallog, "(I4.2)", advance='no')  (i-1)/nstartpn + rank
				do j=1,nopt
					write(fcallog, "(G20.12)", advance='no')  node_xopts((i-1)*nopt+j)
				enddo
				write(fcallog, "(G20.12)", advance='yes') node_fopts(i)
				print *, node_fopts((i-1)*nopt+1:i*nopt), node_fopts(i)
			enddo
			close(unit=fcallog)

			! COMMENTED TO TAKE OUT MPI STUFF:
			! +++++++++++++
			!pull together all of the optimization arguments

!			call mpi_allgather( node_fopts, nstartpn, MPI_DOUBLE, world_fopts,&
!			&		nstartpn, MPI_DOUBLE,mpi_comm_world, ierr)
!			call mpi_allgather( node_xopts, nopt*nstartpn, MPI_DOUBLE, world_xopts,&
!			&		nopt*nstartpn, MPI_DOUBLE,mpi_comm_world, ierr)
!			call mpi_allgather( node_internalopt, ninternalopt*nstartpn, MPI_DOUBLE, world_internalopt,&
!			&		ninternalopt*nstartpn, MPI_DOUBLE,mpi_comm_world, ierr )


			nstarts = (d-1)*(nnode*nstartpn)
			fopt_hist(nstarts + 1:nstarts + nnode*nstartpn) = node_fopts ! with MPI was : world_fopts
			do j =1,(nnode*nstartpn)
				do i =1,nopt
					xopt_hist(i,nstarts + j) = node_xopts((j-1)*nopt + i) ! with MPI was : node_xopts((j-1)*nopt + i)
				enddo
			!	do i=1,ninternalopt
			!		internalopt_hist(i,nstarts+j) = !world_internalopt( (j-1)*ninternalopt + i )
			!	enddo
			enddo


			if( rank .eq. 0) then
				open(unit=fcallog, file = "cal_mlsl.csv" ,ACCESS='APPEND', POSITION='APPEND')
				do i=1,(nstartpn*nnode)
					write(fcallog, "(I4.2)", advance='no')    (i-1)/nstartpn
					do j=1,nopt
						write(fcallog, "(G20.12)", advance='no')  node_xopts((i-1)*nopt+j)  ! with MPI was: world_xopts((i-1)*nopt+j)
					enddo
					write(fcallog, "(G20.12)", advance='yes') node_fopts(i) ! with MPI was: world_fopts(i)
					print *, node_xopts((i-1)*nopt+1:i*nopt), node_fopts(i) ! with MPI was: world_fopts(i)

				enddo
				close(unit=fcallog)
			endif

			!inspect the optima to see if we should keep searching
			!use Bayesian stopping criteria: EW < W+0.5, where W is # local mins, N is # of searches and E[W] = W(N-1)/(N-W-2)
			W= 0._dp
			nstarts = nstarts+nnode*nstartpn !<-number of starts so far
			fval = fopt_hist(1)
			xopt = xopt_hist(:,1)
			j = 1
			do i=1,nstarts
				fdist = dabs( fopt_hist(i)-fval )/dabs(fval)
				xdist = sum(dabs( xopt_hist(:,i)- xopt ) / dabs(xopt))
				if( (fopt_hist(i) .lt. fval) .or. (i==1)) then
					j = i
					if(xdist .ge. nopt*simtol) W = W + 1._dp
				endif
			enddo
			i =j !j indexed the minimnum value
			fval = fopt_hist(i)
			x0 = xopt_hist(:,i)

			!print optimal parameters
			x0w(1:nopt_pars) = xopt !convert to input space
			x0w(1+nopt_pars) = wmean
			call readwrite_glbparams( x0w , "nuxiw_opt.csv" , .false. )
			call set_glbparams(x0)
			!call vec2csv( (/ nu, xizcoef,xiagezcoef, Fd(2,1)/Fd(3,1),Fd(3,1),Fd(2,2)/Fd(2,2),wmean /)  , "nuxiw_opt.csv")

			call mat2csv(fndgrid, "fndgrid_opt.csv")
			call mat2csv(sepgrid, "sepgrid_opt.csv")


			EW = W*dble(nstarts-1)/max( dble(nstarts) - W-2 ,1._dp)
			if( dble(nstarts) - W-2 .le. 0._dp ) EW = nstarts
			print *, "EW: ", EW, " W: ", W
			if( EW < W + 0.5_dp ) then
				exit
			endif

			call vec2csv(fopt_hist(1:nstarts),"fopt_hist.csv")
			call mat2csv(xopt_hist(:,1:nstarts),"xopt_hist.csv")
		enddo

		print_lev = print_lev_old
		verbose = verbose_old

		call vec2csv(fopt_hist,"fopt_hist.csv")
		call mat2csv(xopt_hist,"xopt_hist.csv")

		deallocate(wspace,node_fopts,node_xopts,world_fopts,world_xopts)
		deallocate(world_internalopt,node_internalopt)
		deallocate(seedhr)
		deallocate(sepgrid0,fndgrid0)
		deallocate(mod_xl,mod_xu)

	end subroutine cal_mlsl


	subroutine cal_dist_nloptwrap(fval, nparam, paramvec, gradvec, need_grad, shk)

		integer :: nparam,need_grad
		real(8) :: fval, paramvec(nparam),gradvec(nparam)
		type(shocks_struct) :: shk
		real(8) :: errvec(nparam),paramwt(nparam),paramvecH(nparam),errvecH(nparam),paramvecL(nparam),errvecL(nparam),gradstep(nparam), fvalH,fvalL
		integer :: i, ii,print_lev_old, verbose_old


		print_lev_old = print_lev
		verbose_old = verbose
		print_lev = 0
		verbose = 1
		!if(smth_dicont .le. 20._dp) smth_dicont = smth_dicont*1.05_dp

		if(verbose_old >=1) print *, "test parameter vector ", paramvec
		if(print_lev_old >=1) then
			open(unit=fcallog, file=callog ,ACCESS='APPEND', POSITION='APPEND')
			write(fcallog,*) "test parameter vector ", paramvec
			close(unit=fcallog)
		endif

		call cal_dist(paramvec, errvec,shk)

		if(verbose_old >=1) print *, "         error vector ", errvec
		if(print_lev_old >=1) then
			open(unit=fcallog, file=callog ,ACCESS='APPEND', POSITION='APPEND')
			write(fcallog,*)  "         error vector ", errvec
			close(unit=fcallog)
		endif


		paramwt = 1./dble(nparam)		! equal weight
		paramwt(1) = 2.0_dp
		paramwt(2) = 1.0_dp
		fval = 0.
		do i = 1,nparam
			fval = errvec(i)**2*paramwt(i) + fval
		enddo
		if( need_grad .ne. 0) then
			cal_on_grad = .true.


			do i=1,nparam
				gradstep (i) = min( dabs( paramvec(i)*(5.e-4_dp) ) ,5.e-4_dp)
				paramvecH(i) = paramvec(i) + gradstep(i)
				call cal_dist(paramvecH, errvecH,shk)
				paramvecL(i) = paramvec(i) - gradstep(i)
				call cal_dist(paramvecL, errvecL,shk)
				fvalH = 0.
				fvalL = 0.
				do ii =1,nparam
					fvalH = paramwt(ii)*errvecH(ii)**2 + fvalH
					fvalL = paramwt(ii)*errvecL(ii)**2 + fvalL
				enddo
				gradvec(i) =  (fvalH - fvalL)/(2._dp * gradstep(i))
			enddo
			if(verbose_old >=1) print *, "             gradient ", gradvec
			if(print_lev_old >=1) then
				open(unit=fcallog, file=callog ,ACCESS='APPEND', POSITION='APPEND')
				write(fcallog,*) "             gradient ", gradvec
				close(unit=fcallog)
			endif

			cal_on_grad = .false.
		endif
		print_lev = print_lev_old
		verbose = verbose_old

		cal_obj = fval

	end subroutine cal_dist_nloptwrap


end module find_params

!**************************************************************************************************************!
!**************************************************************************************************************!
! DFBOLS/HYBRJ OBJECTIVE
!**************************************************************************************************************!
!**************************************************************************************************************!
!objective function(s) for DFBOLS goes here. For multiple objectives, use a flag set at module level in find_params

subroutine dfovec(ntheta, mv, theta0, v_err)
	use V0para
	use find_params

	implicit none

	integer, intent(in) 	:: ntheta,mv
	real(dp), dimension(mv)	:: v_err
	real(dp), dimension(ntheta) :: theta0

	real(dp)  :: coef_loc(ntheta),paramvec(ntheta)
	real(dp) :: fval(mv), errvec(mv)
	integer :: lev_der,ncoef_active
	real(dp), allocatable :: coef_est(:),dif_coef(:),reldist_coef(:),wthr(:)
	real(dp), allocatable :: coef_here(:)
	integer :: Ncoef,i,j,ri,ci,ii
	integer :: print_lev_old,verbose_old
	real(dp):: avwt,avonet(Nskill),urt,Efrt,Esrt
	real(dp):: paramwt(mv),dist_urt,dist_frt
	real(dp):: wage_trend_out(Tsim,nj), junk
	character(len=10) :: char_solcoefiter


	print_lev_old = print_lev
	verbose_old = verbose
	print_lev = 0
	verbose = 1



	paramvec = theta0*(mod_xu-mod_xl)+mod_xl

	call cal_dist(paramvec, errvec,mod_shk)
	!do	i=1,mv
	!	call random_number(junk)
	!	errvec(i) = junk
	!end do

	paramwt = 1._dp
	paramwt(1) = 2.0_dp
	paramwt(2) = 1.0_dp
	do i=1,mv
		v_err(i) = errvec(i)*paramwt(i)
	enddo

	print_lev = print_lev_old
	verbose = verbose_old

end subroutine dfovec


!**************************************************************************************************************!
!**************************************************************************************************************!
!						MAIN DRIVER PROGRAM						       !
!**************************************************************************************************************!
!**************************************************************************************************************!

program V0main

	use V0para
	use helper_funs
	use sol_val
	use sim_hists
	use model_data
	use find_params
!	use mpi

	implicit none


	!************************************************************************************************!
	! Counters and Indicies
	!************************************************************************************************!

		integer  :: id=1, it=1, ij=1, itr=1, ial=1, iz=1, i=1,ii,j=1,narg_in=1, wo=1, idi,iter, status=1,t0tT(2)
		character(len=32) :: arg_in
	!************************************************************************************************!
	! Other
	!************************************************************************************************!
		real(dp)	:: wagehere=1.,utilhere=1., junk=1., totapp_dif_hist,ninsur_app,cumpid(nd,nd+1,ndi,TT-1)
		real(dp), allocatable :: wage_coef_0chng(:),wage_tr_occ(:,:)
		real(dp), allocatable :: tr_hist_vec(:),wg_hist_vec(:)
		real(dp) :: xsec_PrDageDel(nd,TT,ndi)
	!************************************************************************************************!
	! Structure to communicate everything
		type(val_struct), target    :: vfs
		type(pol_struct), target    :: pfs
		type(hist_struct), target   :: hst
		type(shocks_struct), target :: shk
		type(moments_struct):: moments_sim
		type(val_pol_shocks_struct) :: vfs_pfs_shk
	! Timers
		integer :: c1=1,c2=1,cr=1,cm=1
		real(dp) :: t1=1.,t2=1.
		logical :: sol_once = .true.
	! NLopt/BOBYQA stuff
		integer(8) :: calopt=0,ires=0, calopt_loc=0,ires_loc=0
		real(dp) :: lb(nopt_tgts),ub(nopt_tgts), ervec(nopt_tgts), erval, param0(nopt_tgts)=1.,err0(nopt_tgts)=1.,lb_1(1),ub_1(1),parvec_1(1),ervec_1(1)
		real(dp), allocatable :: parvec(:)
		integer  :: nodei,nnode,ierr
		!for refinement
		real(dp) :: rhobeg, rhoend,x0(nopt_tgts),xopt(nopt_tgts),x0w(nopt_tgts+1),zeros(nopt_tgts),ones(nopt_tgts)
		real(dp), allocatable :: wspace(:)
		integer :: nopt,ninterppt
		external dgemm, dfovec
		!for ustats
		real(dp) :: urt,Efrt,Esrt, dist_urt,dist_frt

		character(8)  :: date
		character(10) :: time
		character(5)  :: zone

	moments_sim%alloced = 0

	allocate(parvec(nopt_tgts))
!	occ_dat = .false.
	ierr =0
!	call mpi_init(ierr)
!	call mpi_comm_rank(mpi_comm_world,nodei,ierr)
!	call mpi_comm_size(mpi_comm_world,nnode,ierr)
	nopt = nopt_tgts

	print *, "Running version June 9, 2020"
	print *, "Starting on node ", nodei, "out of ", nnode

    call date_and_time(DATE=date,ZONE=zone)
    call date_and_time(TIME=time)

    print '(a,2x,a,2x,a)', date, time, zone

	call setparams()

	call fndsep_data( fnd_seq, sep_seq)

	print *, fnd_seq(1), fnd_seq(Tsim+1)
	print *, sep_seq(1), sep_seq(Tsim+1)

	call gen_new_wgtrend(wage_trend,wage_coef)
	call mat2csv(wage_trend,"wage_trend.csv")
	call vec2csv(wage_lev,"wage_lev.csv")
	allocate(wage_tr_occ(size(wage_trend,1),size(wage_trend,2)))
	do it=1,Tsim
		wtr_tmean_ts(it) = sum(wage_trend(it,:)+occwg_datlev)/dble(nj)
		wage_tr_occ(it,:) = wage_trend(it,:) + wage_lev(:)
	enddo
	call mat2csv(wage_tr_occ,"wtr_occ.csv")
	deallocate(wage_tr_occ)

	narg_in = iargc()
	if( narg_in > 0 ) then
		call getarg(1, arg_in)
		print *, "initial nu=", arg_in
		read( arg_in, * ) nu
		if(narg_in > 1) then
			call getarg(2, arg_in)
			print *, "initial xiage=", arg_in
			read( arg_in, * ) xiagezcoef

			call getarg(3, arg_in)
			print *, "initial xiz =", arg_in
			read( arg_in, * ) xizcoef

		endif
		if(narg_in > 3) then
			call GETARG(4, arg_in)
			read(arg_in, *) run_experiments
		else
			run_experiments = .false.
		endif
		print *, "run experiments", run_experiments
		if(narg_in > 4) then
			call GETARG(5, arg_in)
			read(arg_in, *) run_cal
		else
			run_cal = .false.
		endif
		print *, "run calibration", run_cal
		if(narg_in > 5) then
			call GETARG(6, arg_in)
			read(arg_in, *) refine_cal
		else
			refine_cal = .false.
		endif
		print *, "refine calibration", refine_cal
		if(narg_in > 6) then
			call GETARG(7, arg_in)
			read(arg_in, *) elast_xi
		else
			elast_xi = .false.
		endif
		print *, "calculate xi elasticity", elast_xi


	endif

	caselabel = ""

	agrid(1) = .05*(agrid(1)+agrid(2))
	if(print_lev >= 2) then
		! plot out a bunch of arrays for analyzing VFs, etc
		wo = 0
		do id=1,nd
			call vec2csv(alfgrid(:,id),'alfgrid.csv',wo)
			call mat2csv(pialf(:,:,id),"pial.csv",wo)
			wo=1
		enddo
		wo=0
		call vec2csv(agrid,"agrid.csv",wo)
		call vec2csv(delgrid,'delgrid.csv',wo)
		call vec2csv(trgrid,'trgrid.csv',wo)
		call vec2csv(egrid,'egrid.csv',wo)
		call mat2csv(zgrid,'zgrid.csv',wo)
		call veci2csv(dgrid,'dgrid.csv',wo)
		call vec2csv(agegrid,'agegrid.csv',wo)
		call mat2csv(piz(:,:),"piz.csv",wo)
		call mat2csv(PrDeath,"PrDeath.csv",wo)
		call mat2csv(PrDage,"PrDage.csv",wo)
		call mat2csv(PrDageDel(:,:,1),"PrDageDelL.csv",wo)
		call mat2csv(PrDageDel(:,:,ndi),"PrDageDelH.csv",wo)
		call vec2csv(hazborn_constpop,"hazborn_constpop.csv",wo)
		call vec2csv(hazborn_t,"hazborn_t.csv",wo)
		call vec2csv(prborn_constpop,"prborn_constpop.csv",wo)
		call vec2csv(prborn_t,"prborn_t.csv",wo)
		call mat2csv(tbase_out,"tbase_out.csv",wo)
		cumpid = 0._dp
		do idi=1,ndi
		do it =1,TT-1
		do id =1,nd
			do i =1,nd
				cumpid(id,i+1,idi,it) = pid(id,i,idi,it)+cumpid(id,i,idi,it)
			enddo
		enddo
		enddo
		enddo

		wo=0
		do it = 1,TT-1
			do ij = 1,ndi
				call mat2csv(pid(:,:,ij,it),"pid.csv",wo)
				call mat2csv(cumpid(:,:,ij,it),"cumpid.csv",wo)

				if(wo==0) wo =1
			enddo
		enddo

		call vec2csv(occsz0, "occsz0.csv")
		call vec2csv(occdel, "occdel.csv")
		call mat2csv(transpose(prob_age), "prob_age.csv")
		call mat2csv(occpr_trend,"occpr_trend.csv")
		call mat2csv(occwg_dattrend,"occwg_dattrend.csv")
		call vec2csv(occwg_datlev,"occwg_datlev.csv")
		if(tr_spline .eqv. .true.) then
			call vec2csv(occwg_datsplncoef,"occwg_datsplncoef.csv")
		else
			if(NKpolyT >=2) then
				call mat2csv(occwg_datcoef_sqr,"occwg_datcoefs.csv")
			else
				call vec2csv(occwg_datcoef,"occwg_datcoefs.csv")
			endif
		endif

		open(1, file="wage_dist.csv")
		itr = tri0
		iz  = 3
		ij = 1
		do it = 1,TT-1
			do ial =1,nal
				do id = 1,nd-1
					wagehere = wage(0._dp,alfgrid(ial,id),id,it)
					write(1, "(G20.12)", advance='no') wagehere
				enddo
				id = nd
				wagehere = wage(0._dp,alfgrid(ial,id),id,it)
				write(1,*) wagehere
			enddo
			write(1,*) " "! trailing space
		enddo
		close(1)

		open(1, file="xi.csv")
		open(2, file="xi_hlth.csv")
		ial = nal/2
		do it=1,TT-1
			do id=1,(nd-1)
				wagehere = wage(minval(trgrid),alfgrid(ial,id),id,it)
				write(1, "(G20.12)", advance='no') xifun(id,minval(trgrid),it,0._dp,junk)!xifun(id,wagehere,it,1._dp,junk) !
				write(2, "(G20.12)", advance='no') junk/xifun(id,minval(trgrid),it,0._dp,junk)!junk/xifun(id,wagehere,it,1._dp,junk) !
			enddo
			id = nd
			wagehere = wage(minval(trgrid),alfgrid(ial,id),id,it)
			write(1,*) xifun(id,minval(trgrid),it,0._dp,junk) !xifun(id,wagehere,it,1._dp,junk)
			write(2,*) junk/xifun(id,minval(trgrid),it,0._dp,junk)!junk/xifun(id,wagehere,it,1._dp,junk)
		enddo
		do it=1,TT-1
			do id=1,(nd-1)
				wagehere = wage(maxval(trgrid),alfgrid(ial,id),id,it)
				write(1, "(G20.12)", advance='no') xifun(id,maxval(trgrid),it,0._dp,junk)!xifun(id,wagehere,it,1._dp,junk)
				write(2, "(G20.12)", advance='no') junk/xifun(id,maxval(trgrid),it,0._dp,junk) !junk/ xifun(id,wagehere,it,1._dp,junk)
			enddo
			id = nd
			wagehere = wage(maxval(trgrid),alfgrid(ial,id),id,it)
			write(1,*) xifun(id,maxval(trgrid),it,0._dp,junk) !xifun(id,wagehere,it,1._dp,junk)
			write(2,*) junk/xifun(id,maxval(trgrid),it,0._dp,junk) !xifun(id,wagehere,it,1._dp,junk)
		enddo
		close(1)
		close(2)

		open(1, file="util_dist.csv")
		junk =0.
		itr =tri0
		iz  =2
		do it = 1,TT-1
			do ial =1,nal
				do id = 1,nd-1
					wagehere = wage(0._dp,alfgrid(ial,id),id,it)
					utilhere = util(wagehere,id,1)
					write(1, "(G20.12)", advance='no') utilhere
					junk = utilhere + junk
					utilhere = util(wagehere,id,2)
					write(1, "(G20.12)", advance='no') utilhere
					junk = utilhere + junk

				enddo
				id = nd
				utilhere = util(wagehere,id,1)
				write(1, "(G20.12)", advance='no') utilhere
				junk = utilhere + junk
				utilhere = util(wagehere,id,2)
				write(1, "(G20.12)", advance='yes') utilhere
				junk = utilhere + junk
			enddo
			write(1,*) " "! trailing space
		enddo
		close(1)
	endif
	junk = junk/dble(2*nd*nal*(TT-1))
	!util_const = - junk - util_const

	allocate(tr_hist_vec(Nsim*Tsim))
	allocate(wg_hist_vec(Nsim*Tsim))
	call alloc_shocks(shk)
	call draw_shocks(shk)
	if( run_experiments .eqv. .true. ) then
		sol_once = .true.
	endif

	!************************************************************************************************!
	!solve it once
	!************************************************************************************************!
	if (sol_once .eqv. .true.) then
		if(verbose >=1) then
			call system_clock(count_rate=cr)
			call system_clock(count_max=cm)
			call CPU_TIME(t1)
			call SYSTEM_CLOCK(c1)
		endif
		! set up economy and solve it
		call alloc_econ(vfs,pfs,hst)

		Vtol = 5e-5
		cal_niter = 0
		! wc_guess_lev = 0._dp
		! wc_guess_nolev = 0._dp
		call vec2csv(shk%z_jt_innov,"z_jt_innov_hist_2.csv")
		call set_zjt(hst%z_jt_macroint, hst%z_jt_panel, shk) ! includes call settfp()

		if(verbose >1) print *, "Solving the model"
		call sol(vfs,pfs)

     ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     ! THIS IS AN ENGINEERING FIX TO GET POPULATIONS STABLE
		call sim(vfs, pfs, hst,shk,.false.) !---- lingering issue that 1st stage seems to have correlated results
		xsec_PrDageDel = 0._dp
		do i=1,Nsim
			do it=100,Tsim !totally ad hoc to start at 100
				if(shk%age_hist(i,it)>0 .and. hst%status_hist(i,it)>0) then
					xsec_PrDageDel( shk%d_hist(i,it),shk%age_hist(i,it), shk%del_i_int(i) ) = 1._dp &
						& + xsec_PrDageDel( shk%d_hist(i,it), shk%age_hist(i,it), shk%del_i_int(i) )
				endif
			enddo
		enddo
		do idi=1,ndi
			do it=1,TT
				junk = sum(xsec_PrDageDel(:,it,idi))
				xsec_PrDageDel(:,it,idi) = xsec_PrDageDel(:,it,idi)/junk
			enddo
		enddo
		call mat2csv( xsec_PrDageDel(:,:,1),"xsec_PrDageDelL.csv")
		call mat2csv( xsec_PrDageDel(:,:,ndi),"xsec_PrDageDelH.csv")
		PrDageDel = xsec_PrDageDel
		call set_dit(shk%d_hist,shk%health_it_innov,shk%del_i_int,shk%age_hist)
		call set_alit(shk%al_hist,shk%al_int_hist, shk%al_it_innov,shk%d_hist, status)

		! if(print_lev>1) then
		! 	call vec2csv(wc_guess_lev,"wc_guess_lev.csv")
		! 	call vec2csv(wc_guess_nolev,"wc_guess_nolev.csv")
		! endif

		if(verbose >=1) print *, "Simulating the model"
		call sim(vfs, pfs, hst,shk)
		if(verbose >=1) print *, "Computing moments"
		call moments_compute(hst,moments_sim,shk)
		if(verbose >=1) print *, "DI rate" , moments_sim%avg_di
		if(verbose >1 .or. print_lev >1 ) then
			call comp_ustats(hst,shk,urt,Efrt,Esrt)
			dist_urt = (urt - avg_unrt)/avg_unrt
			dist_frt = (Efrt - avg_frt)/avg_frt
			print * , "missing on fidning rate ", dist_frt, " missing on u rate ", dist_urt
		endif
		!	set mean wage:
		wmean = 0._dp
		junk = 0._dp
		ii =1
		do i=1,Nsim
			do it=1,Tsim
				wagehere = hst%wage_hist(i,it)
				if( hst%status_hist(i,it)==1) then
					wmean = wagehere + wmean
					junk = junk+1.
					tr_hist_vec(ii) = wage_trend(it,shk%j_i(i))
					wg_hist_vec(ii) = wagehere
					ii = ii+1
				endif
			enddo
		enddo
		wmean = wmean/junk
		print *, "average wage:", wmean
		do i=1,ne
			egrid(i) = egrid(i)*wmean
		enddo
		!fill in tr_decls
		call quicksort(tr_hist_vec,1,ii-1)
		tr_decls(1) = tr_hist_vec(2)
		tr_decls(11)= tr_hist_vec(ii-2)
		do i =2,10
			tr_decls(i) = tr_hist_vec((i*ii-ii)/10 )
		enddo
		call quicksort(wg_hist_vec,1,ii-1)
		wg_decls(1) = wg_hist_vec(2)
		wg_decls(11)= wg_hist_vec(ii-2)
		do i =2,10
			tr_decls(i) = tr_hist_vec((i*ii-ii)/10 )
		enddo
		call vec2csv(tr_decls,"tr_decls.csv")
		call vec2csv(wg_decls,"wg_decls.csv")
		ii =1
		do it=1,Tsim
			wtr_tmean_ts(it) = 0._dp
			junk = 0._dp
			ii = 0
			do i=1,Nsim
				if( hst%status_hist(i,it)>0 .and. shk%age_hist(i,it)>0) then
					tr_hist_vec(ii) = wage_trend(it,shk%j_i(i))
					ii = ii+1
					junk = junk + 1._dp
					wtr_tmean_ts(it) = wtr_tmean_ts(it) + shk%wtr_it(i,it) !wage_trend(it,shk%j_i(i)) + wage_lev(shk%j_i(i))
				endif
			enddo
			!call quicksort(tr_hist_vec,1,ii-1)
			!wtr_tmean_ts(it) = tr_hist_vec((ii-1)/2)
			wtr_tmean_ts(it) = wtr_tmean_ts(it)/dble(ii)
		enddo
		call vec2csv(wtr_tmean_ts, "wtr_tmean_ts0.csv")
		! if(run_experiments .eqv. .false.) then
		! 	call dealloc_econ(vfs,pfs,hst)
		! endif

		if(verbose >=1) then
			call CPU_TIME(t2)
			call SYSTEM_CLOCK(c2)
			if(nodei ==0)  print *, "Sol Once, System Time", dble(c2-c1)/dble(cr)
			if(nodei ==0)  print *, "Sol Once, CPU Time   ", (t2-t1)
		endif
		print *, "error 1:  ",	(moments_sim%init_diaward - diaward_target)/diaward_target
		print *, "error 2:  ",	(moments_sim%init_hlth_acc - hlth_acc_target)/hlth_acc_target
		print *, "error 3:  ",	(moments_sim%old_initdiawdfrac - old_target)/old_target
		print *, "error 4:  ", 	(moments_sim%avg_diprob - allowrt_target)/allowrt_target

		!print *, "error 3: ",	(moments_sim%d1_diawardfrac - d1_diawardfrac_target)/d1_diawardfrac_target
		!print *, "error 3: ",	(moments_sim%work_rateDage(1,1)/moments_sim%work_rateDage(2,1) - 1./p1d2_2545target)*2/(moments_sim%work_rateDage(1,1)/moments_sim%work_rateDage(2,1) + 1./p1d2_2545target)
		!print *, "error 4: ",	(moments_sim%work_rateDage(1,1)/moments_sim%work_rateDage(3,1) - 1./p1d3_2545target)*2/(moments_sim%work_rateDage(1,1)/moments_sim%work_rateDage(3,1) + 1./p1d3_2545target)
		!print *, "error 5: ",	(moments_sim%work_rateDage(1,1)/moments_sim%work_rateDage(1,2) - 1./p1d1_4665target)*2/(moments_sim%work_rateDage(1,1)/moments_sim%work_rateDage(1,2) + 1./p1d1_4665target)
		!print *, "error 6: ",	(moments_sim%work_rateDage(1,1)/moments_sim%work_rateDage(2,2) - 1./p1d2_4665target)*2/(moments_sim%work_rateDage(1,1)/moments_sim%work_rateDage(2,2) + 1./p1d2_4665target)
		!print *, "error 7: ",	(moments_sim%work_rateDage(1,1)/moments_sim%work_rateDage(3,2) - 1./p1d3_4665target)*2/(moments_sim%work_rateDage(1,1)/moments_sim%work_rateDage(3,2) + 1./p1d3_4665target)
		print *, "error 5:  ",	 (moments_sim%diaward_ageeffect - award_age_target)/award_age_target
		print *, "error 6:  ",	 (moments_sim%work_rateD(1)/moments_sim%work_rateDage(1,1) - p1d1_target  )*2/ ( moments_sim%work_rateD(1)/moments_sim%work_rateDage(1,1) + p1d1_target )
		print *, "error 7:  ",	 (moments_sim%work_rateD(2)/moments_sim%work_rateDage(1,1) - p1d2_target  )*2/ ( moments_sim%work_rateD(2)/moments_sim%work_rateDage(1,1)  + p1d2_target )
		print *, "error 8:  ",	 (moments_sim%work_rateD(3)/moments_sim%work_rateDage(1,1) - p1d3_target  )*2/ ( moments_sim%work_rateD(3)/moments_sim%work_rateDage(1,1)  + p1d3_target )
		print *, "error 9:  ",   (moments_sim%d1_diawardfrac - d1_diawardfrac_target)/d1_diawardfrac_target
		print *, "error 10: ",   (moments_sim%d2_diawardfrac - d2_diawardfrac_target)/d2_diawardfrac_target
		print *, "error 11: ",   (moments_sim%d3_diawardfrac - d3_diawardfrac_target)/d3_diawardfrac_target
		print *, "error 12: ",   (moments_sim%end_hlth_acc - hlth_acc_end_target)/hlth_acc_end_target


		call dealloc_econ(vfs,pfs,hst)

	endif !sol_once


	!bounds for paramvec:

	!         nud1,   nud2,     nud3 , nuage  ,  xizcoef,   xiage    , xiVtrend , xid1coef  ,  xid2coef ,  xid3coef ,    Fd(2,2) ,       Fd(3,2)
	lb = (/ 0.00_dp, 0.00_dp,-0.20_dp,-1.00_dp, -2.50_dp,   -2.000_dp,  1.0_dp  , -2.50_dp  , -2.50_dp  , -2.50_dp  ,    0.001_dp,       0.01_dp /)
	ub = (/ 4.00_dp, 2.00_dp, 2.00_dp, 2.00_dp,  2.50_dp,    2.00_dp ,  2.5_dp  , 2.50_dp   ,  2.50_dp  ,  2.50_dp  ,    3.500_dp,       5.01_dp /)

	x0w(1:nopt_pars) = (lb+ub)/2._dp
	x0w(1+nopt_pars) = wmean
	call readwrite_glbparams( x0w, "nuxiw_init.csv", .false. )
	if( (run_cal .eqv. .true.) .and. (dbg_skip .eqv. .false.) ) then
		call system_clock(count_rate=cr)
		call system_clock(count_max=cm)
		call system_clock(c1)
		parvec = 0._dp !will store the optimal values. Will be over-written
		call cal_mlsl( erval, parvec, nopt, lb, ub, shk)
		call set_glbparams(parvec)


		call system_clock(c2)
		if(nodei ==0)  print *, "Calibration, Wall Time in hours ", dble(c2-c1)/dble(cr)/360._dp
	endif

	!****************************************************************************
	!   Now run some experiments:
	print_lev = 2
	verbose = 1
	!read in the optimal inner stuff:
	!wage_trend <- "wage_trend_opt.csv"
	!wage_coef <- "wage_coef_opt.csv")
	!fndgrid <- "fndgrid_opt.csv")
	!sepgrid <- "sepgrid_opt.csv")

	print *, "about to read in parameters "

	!read in xi, nu
	call readwrite_glbparams( x0w, "nuxiw_opt.csv", .true. )
	x0 = x0w(1:nopt_pars)
	wmean = x0w(1+nopt_pars)
	call set_glbparams(x0)
	if(ntr>3) then
		wtr_scale = trgrid(ntr-1) - trgrid(2)
	else
		wtr_scale = trgrid(ntr) - trgrid(1)
	endif

	if (refine_cal .eqv. .true. ) then

		!puts the parameters into x0
		call readwrite_glbparams( x0w, "nuxiw_opt.csv", .true. )
		x0 = x0w(1:nopt_pars)
		wmean = x0w(1+nopt_pars)
		call set_glbparams(x0)
		! if(nopt_tgts==5) &
		! &	x0 = (/ nu,  xizcoef, Fd(2)/Fd(3), Fd(3),xiagecoef /)

		xopt = (x0-lb)/(ub-lb)
		nopt = size(xopt)
		ninterppt = 2*nopt+1
		mod_shk => shk
		allocate(mod_xl(nopt))
		allocate(mod_xu(nopt))
		mod_xl = lb ! this is so bounds are available in dfovec
		mod_xu = ub
		zeros = max(xopt - 0.075_dp, 0._dp)
		ones  = min(xopt + 0.075_dp, 1._dp)

		allocate(wspace((ninterppt+5)*(ninterppt+nopt)+3*nopt*(nopt+5)/2))

		rhobeg = .05
		rhoend = rhobeg*1e-3_dp
		do i=1,nopt
			if( (ones(i)-zeros(i) < 2._dp * rhobeg)  .and. zeros(i) .le. 0._dp ) ones(i)  = zeros(i)+ 0.15_dp
			if( (ones(i)-zeros(i) < 2._dp * rhobeg)  .and. ones(i)  .ge. 1._dp ) zeros(i) = ones(i) - 0.15_dp
		enddo

		call bobyqa_h(nopt,ninterppt,xopt,zeros,ones, &
		&	rhobeg,rhoend,1,50,wspace,nopt)

		x0w(1:nopt_tgts) = xopt*(ub-lb)+ lb
		x0w(1+nopt_tgts) = wmean
		call readwrite_glbparams(x0w,"nuxiw_opt_refine.csv", .false. )


		deallocate(mod_xl,mod_xu)
		deallocate(wspace)

	endif


	if (dbg_skip .eqv. .false.) then
		welfare_cf = .true.
		if(refine_cal .eqv. .false. ) then
			call readwrite_glbparams( x0w, "nuxiw_opt.csv", .true. )
		else
			call readwrite_glbparams( x0w, "nuxiw_opt_refine.csv", .true. )
		endif
		parvec = x0w( 1:nopt_tgts )

		!if(nopt_tgts==6) parvec = (/nu,xizcoef, xiagecoef,Fd(2)/Fd(3),Fd(3),xiagecoef /)
		call gen_new_wgtrend(wage_trend,wage_coef)
		caselabel = ""
	 	print *, caselabel, " ---------------------------------------------------"
		occ_dat = .true.
		wtr_by_occ = .true.
		demog_dat = .true.
		buscyc = .true.
		print_lev = 2
		call cal_dist(parvec,err0,shk)
		print *, "error in initial "
		do i=1,nopt_tgts
			write(*,"(G20.12)", advance = 'no') err0(i)
			write(*,"(G20.12)", advance = 'no') ", "
		enddo
		write(*,*) " "

		print *, "---------------------------------------------------"

		welfare_cf = .false.
	endif

	if( (nodei == 0) .and. (dbg_skip .eqv. .false.) ) then
		!compute counter-factual with 0 correlation between occupation and health risk
		welfare_cf = .true.

		del_by_occ = .false.

		if(refine_cal .eqv. .false. ) then
			call readwrite_glbparams( x0w, "nuxiw_opt.csv", .true. )
		else
			call readwrite_glbparams( x0w, "nuxiw_opt_refine.csv", .true. )
		endif
		parvec = x0w( 1:nopt_pars )

		call gen_new_wgtrend(wage_trend,wage_coef)
		caselabel = "delocc_CF"
	 	print *, caselabel, " ---------------------------------------------------"
		occ_dat = .true.
		wtr_by_occ = .true.
		demog_dat = .true.
		buscyc = .true.

		call set_deli( shk%del_i_int,shk%del_i_draw,shk%j_i)
		call set_dit( shk%d_hist, shk%health_it_innov, shk%del_i_int, shk%age_hist)
		call set_alit(shk%al_hist,shk%al_int_hist, shk%al_it_innov,shk%d_hist, status)

		print_lev = 2
		call cal_dist(parvec,err0,shk)
		print *, "error: "
		do i=1,nopt_tgts
			write(*,"(G20.12)", advance = 'no') err0(i)
			write(*,"(G20.12)", advance = 'no') ", "
		enddo
		write(*,*) " "
		print *, "---------------------------------------------------"

		del_by_occ = .true.
		welfare_cf = .false.

		call set_deli( shk%del_i_int,shk%del_i_draw,shk%j_i)
		call set_dit( shk%d_hist, shk%health_it_innov, shk%del_i_int, shk%age_hist)
		call set_alit(shk%al_hist,shk%al_int_hist, shk%al_it_innov,shk%d_hist, status)

	endif

	if((nodei == 0) .and. (elast_xi .eqv. .true.) .and. (dbg_skip .eqv. .false.) ) then

		allocate(wage_coef_0chng(size(wage_coef)))
   	 	wage_coef_0chng = wage_coef
		if(tr_spline .eqv. .true. )then
			do i=(Nskill+1),size(wage_coef) !turn-off wage trends
				wage_coef_0chng(i)= 0._dp
			enddo
		else
			do i=(1+NTpolyT),size(wage_coef) !only turn off the occupation-specific trends
	   			wage_coef_0chng(i)= 0._dp
	   		enddo
		endif

		!without any drivers
	 	print *, caselabel, " ---------------------------------------------------"

		wtr_by_occ = .false.
	 	occ_dat = .false.
	 	demog_dat  = .false.
		buscyc = .false.

		call gen_new_wgtrend(wage_trend,wage_coef_0chng)

		call set_age(shk%age_hist, shk%born_hist, shk%age_draw)
		call set_ji( shk%j_i,shk%jshock_ij,shk%born_hist)
		call set_fndsepi(shk%fndsep_i_int,shk%fndsep_i_draw,shk%j_i)
	 	call set_deli( shk%del_i_int,shk%del_i_draw,shk%j_i)
		call set_dit( shk%d_hist, shk%health_it_innov, shk%del_i_int, shk%age_hist)
		call set_alit(shk%al_hist,shk%al_int_hist, shk%al_it_innov,shk%d_hist, status)

		caselabel = "xi0L"

		if(refine_cal .eqv. .false. ) then
			call readwrite_glbparams( x0w, "nuxiw_opt.csv", .true. )
		else
			call readwrite_glbparams( x0w, "nuxiw_opt_refine.csv", .true. )
		endif
		parvec = x0w( 1:nopt_pars )

		call cal_dist(parvec,err0,shk)

		print *, "error with low xi ,"
		do i=1,nopt_tgts
			write(*,"(G20.12)", advance = 'no') err0(i)
			write(*,"(G20.12)", advance = 'no') ", "
		enddo
		write(*,*) " "
		print *, "---------------------------------------------------"

		caselabel = "xi0H"


		if(refine_cal .eqv. .false. ) then
			call readwrite_glbparams( x0w, "nuxiw_opt.csv", .true. )
		else
			call readwrite_glbparams( x0w, "nuxiw_opt_refine.csv", .true. )
		endif
		parvec = x0w( 1:nopt_pars )

		call cal_dist(parvec,err0,shk)
		print *, "error with high xi "
		do i=1,nopt_tgts
			write(*,"(G20.12)", advance = 'no') err0(i)
			write(*,"(G20.12)", advance = 'no') ", "
		enddo
		write(*,*) " "
		print *, "---------------------------------------------------"

		call gen_new_wgtrend(wage_trend,wage_coef)

		!now do it with end-period wage trends

		deallocate(wage_coef_0chng)

	endif

	if((nodei == 0) .and. (run_experiments .eqv. .true.) .and. (dbg_skip .eqv. .false.)) then


		if(refine_cal .eqv. .false. ) then
			call readwrite_glbparams( x0w, "nuxiw_opt.csv", .true. )
		else
			call readwrite_glbparams( x0w, "nuxiw_opt_refine.csv", .true. )
		endif
		parvec = x0w( 1:nopt_tgts )
		allocate(wage_coef_0chng(size(wage_coef)))
   	 	wage_coef_0chng = wage_coef
		if(tr_spline .eqv. .true. )then
			do i=(Nskill+1),size(wage_coef) !turn-off wage trends
				wage_coef_0chng(i)= 0._dp
			enddo
		else
			do i=(1+NTpolyT),size(wage_coef) !only turn off the occupation-specific trends
	   			wage_coef_0chng(i)= 0._dp
	   		enddo
		endif

		!without any drivers
		caselabel = "all0"
	 	print *, caselabel, " ---------------------------------------------------"
		if(print_lev .gt. 1) call vec2csv( wage_coef_0chng ,"wage_coef"//trim(caselabel)//".csv")
		if(print_lev .gt. 1) call mat2csv( wage_trend ,"wage_trend"//trim(caselabel)//".csv")
	 	wtr_by_occ = .false.
	 	occ_dat = .false.
	 	demog_dat  = .false.
		buscyc = .false.
		verbose = 1
		print_lev = 2
		print *, "---------------------------------------------------"
		
		call gen_new_wgtrend(wage_trend,wage_coef_0chng)
		
		call set_age(shk%age_hist, shk%born_hist, shk%age_draw)
		call set_ji( shk%j_i,shk%jshock_ij,shk%born_hist)
		call set_fndsepi(shk%fndsep_i_int,shk%fndsep_i_draw,shk%j_i)
	 	call set_deli( shk%del_i_int,shk%del_i_draw,shk%j_i)
		call set_dit( shk%d_hist, shk%health_it_innov, shk%del_i_int, shk%age_hist)
		call set_alit(shk%al_hist,shk%al_int_hist, shk%al_it_innov,shk%d_hist, status)

		call cal_dist(parvec,err0,shk)
	 	print *, "error in targets with ", caselabel,  " ", err0
		call gen_new_wgtrend(wage_trend,wage_coef)
		print *, "---------------------------------------------------"

		! without wage trend
		caselabel = "wchng0"
	 	print *, caselabel, " ---------------------------------------------------"
		if(print_lev .ge. 1) call mat2csv( wage_trend ,"wage_trend"//trim(caselabel)//".csv")
	 	wtr_by_occ = .false.
		buscyc = .true.
	 	occ_dat = .true.
	 	demog_dat  = .true.
		verbose = 1
		print *, "---------------------------------------------------"
		
		call gen_new_wgtrend(wage_trend,wage_coef_0chng)
		
		call set_age(shk%age_hist, shk%born_hist, shk%age_draw)
		call set_ji( shk%j_i,shk%jshock_ij,shk%born_hist)
		call set_fndsepi(shk%fndsep_i_int,shk%fndsep_i_draw,shk%j_i)
	 	call set_deli( shk%del_i_int,shk%del_i_draw,shk%j_i)
		call set_dit( shk%d_hist, shk%health_it_innov, shk%del_i_int, shk%age_hist)
		call set_alit(shk%al_hist,shk%al_int_hist, shk%al_it_innov,shk%d_hist, status)
		call cal_dist(parvec,err0,shk)
	 	print *, "error in targets with ", caselabel,  " ", err0
		call gen_new_wgtrend(wage_trend,wage_coef)
		print *, "---------------------------------------------------"

	 	! without the change in occupation composition
		buscyc = .true.
		occ_dat = .false.
	 	wtr_by_occ = .true.
	 	demog_dat = .true.
	 	caselabel = "deloc0"
	 	print *, caselabel, " ---------------------------------------------------"
		call set_age(shk%age_hist, shk%born_hist, shk%age_draw)
		call set_ji( shk%j_i,shk%jshock_ij,shk%born_hist)
		call set_fndsepi(shk%fndsep_i_int,shk%fndsep_i_draw,shk%j_i)
	 	call set_deli( shk%del_i_int,shk%del_i_draw,shk%j_i)
		call set_dit( shk%d_hist, shk%health_it_innov, shk%del_i_int, shk%age_hist)
		call set_alit(shk%al_hist,shk%al_int_hist, shk%al_it_innov,shk%d_hist, status)
		if(print_lev .ge. 1) call mat2csv( wage_trend ,"wage_trend"//trim(caselabel)//".csv")
		call cal_dist(parvec,err0,shk)
		print *, "error in targets with ", caselabel, " ", err0
		print *, "---------------------------------------------------"

		!without demographic changes
		buscyc = .true.
	 	occ_dat = .true.
	 	wtr_by_occ = .true.
	 	demog_dat = .false.
	 	caselabel = "demog0"
	 	print *, caselabel, " ---------------------------------------------------"
		call set_age(shk%age_hist, shk%born_hist, shk%age_draw)
		call set_ji( shk%j_i,shk%jshock_ij,shk%born_hist)
		call set_fndsepi(shk%fndsep_i_int,shk%fndsep_i_draw,shk%j_i)
	 	call set_deli( shk%del_i_int,shk%del_i_draw,shk%j_i)
		call set_dit( shk%d_hist, shk%health_it_innov, shk%del_i_int, shk%age_hist)
		call set_alit(shk%al_hist,shk%al_int_hist, shk%al_it_innov,shk%d_hist, status)
		call cal_dist(parvec,err0,shk)
		print *, "error in targets with ", caselabel, " ", err0
	 	print *, "---------------------------------------------------"

		!without business cycles
		buscyc = .false.
		occ_dat = .true.
	 	wtr_by_occ = .true.
	 	demog_dat = .true.
	 	caselabel = "bcyc0"
		print *, caselabel, " ---------------------------------------------------"
		call set_age(shk%age_hist, shk%born_hist, shk%age_draw)
		call set_ji( shk%j_i,shk%jshock_ij,shk%born_hist)
		call set_fndsepi(shk%fndsep_i_int,shk%fndsep_i_draw,shk%j_i)
	 	call set_deli( shk%del_i_int,shk%del_i_draw,shk%j_i)
		call set_dit( shk%d_hist, shk%health_it_innov, shk%del_i_int, shk%age_hist)
		call set_alit(shk%al_hist,shk%al_int_hist, shk%al_it_innov,shk%d_hist, status)
		call cal_dist(parvec,err0,shk)
		print *, "error in targets with ", caselabel, " ", err0
	 	print *, "---------------------------------------------------"


	 	! without either the change in occupation or wage trend
	 	occ_dat = .false.
	 	wtr_by_occ = .false.
	 	demog_dat = .true.
		buscyc = .true.
		call gen_new_wgtrend(wage_trend,wage_coef_0chng)
	 	caselabel = "wchng0deloc0"
		print *, caselabel, " ---------------------------------------------------"
		call set_age(shk%age_hist, shk%born_hist, shk%age_draw)
		call set_ji( shk%j_i,shk%jshock_ij,shk%born_hist)
		call set_fndsepi(shk%fndsep_i_int,shk%fndsep_i_draw,shk%j_i)
	 	call set_deli( shk%del_i_int,shk%del_i_draw,shk%j_i)
		call set_dit( shk%d_hist, shk%health_it_innov, shk%del_i_int, shk%age_hist)
		call set_alit(shk%al_hist,shk%al_int_hist, shk%al_it_innov,shk%d_hist, status)
		call cal_dist(parvec,err0,shk)
		call gen_new_wgtrend(wage_trend,wage_coef)
		print *, "error in targets with ", caselabel, " ", err0
	 	print *, "---------------------------------------------------"

		!without either the demograhpic change or wage trend
		occ_dat = .true.
		wtr_by_occ = .false.
		call gen_new_wgtrend(wage_trend,wage_coef_0chng)
		demog_dat = .false.
		buscyc = .true.
		caselabel = "wchng0demog0"
		print *, caselabel, " ---------------------------------------------------"
		call set_age(shk%age_hist, shk%born_hist, shk%age_draw)
		call set_ji( shk%j_i,shk%jshock_ij,shk%born_hist)
		call set_fndsepi(shk%fndsep_i_int,shk%fndsep_i_draw,shk%j_i)
	 	call set_deli( shk%del_i_int,shk%del_i_draw,shk%j_i)
		call set_dit( shk%d_hist, shk%health_it_innov, shk%del_i_int, shk%age_hist)
		call set_alit(shk%al_hist,shk%al_int_hist, shk%al_it_innov,shk%d_hist, status)
		call cal_dist(parvec,err0,shk)
		call gen_new_wgtrend(wage_trend,wage_coef)
		print *, "error in targets with ", caselabel, " ", err0
		print *, "---------------------------------------------------"

	 	occ_dat = .false.
	 	wtr_by_occ = .true.
	 	demog_dat = .false.
		buscyc = .true.
	 	caselabel = "deloc0demog0"
	 	print *, caselabel, " ---------------------------------------------------"
		call set_age(shk%age_hist, shk%born_hist, shk%age_draw)
		call set_ji( shk%j_i,shk%jshock_ij,shk%born_hist)
		call set_fndsepi(shk%fndsep_i_int,shk%fndsep_i_draw,shk%j_i)
	 	call set_deli( shk%del_i_int,shk%del_i_draw,shk%j_i)
		call set_dit( shk%d_hist, shk%health_it_innov, shk%del_i_int, shk%age_hist)
		call set_alit(shk%al_hist,shk%al_int_hist, shk%al_it_innov,shk%d_hist, status)
		call cal_dist(parvec,err0,shk)
		print *, "error in targets with ", caselabel, " ", err0
	 	print *, "---------------------------------------------------"

		! without business cycles or wage trend
	 	occ_dat = .true.
	 	wtr_by_occ = .false.
	 	demog_dat = .true.
		buscyc = .false.
		call gen_new_wgtrend(wage_trend,wage_coef_0chng)
	 	caselabel = "wchng0bcyc0"
		print *, caselabel, " ---------------------------------------------------"
		call set_age(shk%age_hist, shk%born_hist, shk%age_draw)
		call set_ji( shk%j_i,shk%jshock_ij,shk%born_hist)
		call set_fndsepi(shk%fndsep_i_int,shk%fndsep_i_draw,shk%j_i)
	 	call set_deli( shk%del_i_int,shk%del_i_draw,shk%j_i)
		call set_dit( shk%d_hist, shk%health_it_innov, shk%del_i_int, shk%age_hist)
		call set_alit(shk%al_hist,shk%al_int_hist, shk%al_it_innov,shk%d_hist, status)
		call cal_dist(parvec,err0,shk)
		call gen_new_wgtrend(wage_trend,wage_coef)
		print *, "error in targets with ", caselabel, " ", err0
	 	print *, "---------------------------------------------------"

		! without either the change in demographics or business cycles
	 	occ_dat = .true.
	 	wtr_by_occ = .true.
	 	demog_dat = .false.
		buscyc = .false.
		caselabel = "demog0bcyc0"
		print *, caselabel, " ---------------------------------------------------"
		call set_age(shk%age_hist, shk%born_hist, shk%age_draw)
		call set_ji( shk%j_i,shk%jshock_ij,shk%born_hist)
		call set_fndsepi(shk%fndsep_i_int,shk%fndsep_i_draw,shk%j_i)
	 	call set_deli( shk%del_i_int,shk%del_i_draw,shk%j_i)
		call set_dit( shk%d_hist, shk%health_it_innov, shk%del_i_int, shk%age_hist)
		call set_alit(shk%al_hist,shk%al_int_hist, shk%al_it_innov,shk%d_hist, status)
		call cal_dist(parvec,err0,shk)
		print *, "error in targets with ", caselabel, " ", err0
	 	print *, "---------------------------------------------------"

		!without either occupational changes or business cycles
	 	occ_dat = .false.
	 	wtr_by_occ = .true.
	 	demog_dat = .true.
		buscyc = .false.
	 	caselabel = "deloc0bcyc0"
	 	print *, caselabel, " ---------------------------------------------------"
		call set_age(shk%age_hist, shk%born_hist, shk%age_draw)
		call set_ji( shk%j_i,shk%jshock_ij,shk%born_hist)
		call set_fndsepi(shk%fndsep_i_int,shk%fndsep_i_draw,shk%j_i)
	 	call set_deli( shk%del_i_int,shk%del_i_draw,shk%j_i)
		call set_dit( shk%d_hist, shk%health_it_innov, shk%del_i_int, shk%age_hist)
		call set_alit(shk%al_hist,shk%al_int_hist, shk%al_it_innov,shk%d_hist, status)
		call cal_dist(parvec,err0,shk)
		print *, "error in targets with ", caselabel, " ", err0
	 	print *, "---------------------------------------------------"

		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		! Turn on just one at a time
		! with only wage trend
		caselabel = "deloc0demog0bcyc0"
		print *, caselabel, " ---------------------------------------------------"
		wtr_by_occ = .true.
		buscyc = .false.
		occ_dat = .false.
		demog_dat  = .false.
		verbose = 1
		print *, "---------------------------------------------------"
		call set_age(shk%age_hist, shk%born_hist, shk%age_draw)
		call set_ji( shk%j_i,shk%jshock_ij,shk%born_hist)
		call set_fndsepi(shk%fndsep_i_int,shk%fndsep_i_draw,shk%j_i)
		call set_deli( shk%del_i_int,shk%del_i_draw,shk%j_i)
		call set_dit( shk%d_hist, shk%health_it_innov, shk%del_i_int, shk%age_hist)
		call set_alit(shk%al_hist,shk%al_int_hist, shk%al_it_innov,shk%d_hist, status)
		call cal_dist(parvec,err0,shk)
		print *, "error in targets with ", caselabel,  " ", err0
		print *, "---------------------------------------------------"

		! with only the change in occupation composition
		buscyc = .false.
		occ_dat = .true.
		wtr_by_occ = .false.
		demog_dat = .false.
		call gen_new_wgtrend(wage_trend,wage_coef_0chng)
		caselabel = "wchng0demog0bcyc0"
		print *, caselabel, " ---------------------------------------------------"
		call set_age(shk%age_hist, shk%born_hist, shk%age_draw)
		call set_ji( shk%j_i,shk%jshock_ij,shk%born_hist)
		call set_fndsepi(shk%fndsep_i_int,shk%fndsep_i_draw,shk%j_i)
		call set_deli( shk%del_i_int,shk%del_i_draw,shk%j_i)
		call set_dit( shk%d_hist, shk%health_it_innov, shk%del_i_int, shk%age_hist)
		call set_alit(shk%al_hist,shk%al_int_hist, shk%al_it_innov,shk%d_hist, status)
		if(print_lev .ge. 1) call mat2csv( wage_trend ,"wage_trend"//trim(caselabel)//".csv")
		call cal_dist(parvec,err0,shk)
		call gen_new_wgtrend(wage_trend,wage_coef)
		print *, "error in targets with ", caselabel, " ", err0
		print *, "---------------------------------------------------"

		!without demographic changes
		buscyc = .false.
		occ_dat = .false.
		wtr_by_occ = .false.
		demog_dat = .true.
		call gen_new_wgtrend(wage_trend,wage_coef_0chng)
		caselabel = "wchng0deloc0bcyc0"
		call gen_new_wgtrend(wage_trend,wage_coef_0chng)
		print *, caselabel, " ---------------------------------------------------"
		call set_age(shk%age_hist, shk%born_hist, shk%age_draw)
		call set_ji( shk%j_i,shk%jshock_ij,shk%born_hist)
		call set_fndsepi(shk%fndsep_i_int,shk%fndsep_i_draw,shk%j_i)
		call set_deli( shk%del_i_int,shk%del_i_draw,shk%j_i)
		call set_dit( shk%d_hist, shk%health_it_innov, shk%del_i_int, shk%age_hist)
		call set_alit(shk%al_hist,shk%al_int_hist, shk%al_it_innov,shk%d_hist, status)
		call cal_dist(parvec,err0,shk)
		call gen_new_wgtrend(wage_trend,wage_coef)
		print *, "error in targets with ", caselabel, " ", err0
		print *, "---------------------------------------------------"

		!with only business cycles
		buscyc = .true.
		occ_dat = .false.
		wtr_by_occ = .false.
		demog_dat = .false.
		call gen_new_wgtrend(wage_trend,wage_coef_0chng)
		caselabel = "wchng0deloc0demog0"
		print *, caselabel, " ---------------------------------------------------"
		call set_age(shk%age_hist, shk%born_hist, shk%age_draw)
		call set_ji( shk%j_i,shk%jshock_ij,shk%born_hist)
		call set_fndsepi(shk%fndsep_i_int,shk%fndsep_i_draw,shk%j_i)
		call set_deli( shk%del_i_int,shk%del_i_draw,shk%j_i)
		call set_dit( shk%d_hist, shk%health_it_innov, shk%del_i_int, shk%age_hist)
		call set_alit(shk%al_hist,shk%al_int_hist, shk%al_it_innov,shk%d_hist, status)
		call cal_dist(parvec,err0,shk)
		call gen_new_wgtrend(wage_trend,wage_coef)
		print *, "error in targets with ", caselabel, " ", err0
		print *, "---------------------------------------------------"
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

		deallocate(wage_coef_0chng)
	endif

  	!****************************************************************************!
	! IF you love something....
	!****************************************************************************!
	deallocate(tr_hist_vec,wg_hist_vec)
	deallocate(parvec)
	call dealloc_shocks(shk)

!	call mpi_finalize(ierr)



!    .----.   @   @
!   / .-"-.`.  \v/
!   | | '\ \ \_/ )
! ,-\ `-.' /.'  /
!'---`----'----'

End PROGRAM
