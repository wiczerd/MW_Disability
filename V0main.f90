! V0main.f90

!************************************************************************************************!
! @ Amanda Michaud, v1: 10/6/2014; current: 10/31/2014
! @ David Wiczer, v2: 03/31/2015
!-----------------------------------------------------
!This is a sub-program for DIchoices paper 
!	Objective: compute V0(j), the expected value of choosing occupation j in {1,2...J}
!				  before individual types are drawn.
!	
!************************************************************************************************!
! compiler line: gfortran -fopenmp -ffree-line-length-none -g V0para.f90 V0main.f90 -lblas -llapack -lgomp -o V0main.out  
! val grind line: valgrind --leak-check=yes --error-limit=no --log-file=V0valgrind.log ./V0main.out &
module helper_funs
	
	use V0para
	implicit none
	
	!**********************************************************!
	!Public Policy Functions
	!	1)UI(e)		Unemployment Insurance
	!	2)SSDI(e,a)	DI (can be asset tested)
	!	3)SSI(e)	Social Sec. Retirement
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
	!	c) moments_structs
	!	d) hist_struct
	!**********************************************************!
	
	
	!------------------------------------------------------------------
	! a) val_struct: VR, VD, VN, VW, VU, V
	!------------------------------------------------------------------
	type val_struct
		real(dp), allocatable:: 	VR(:,:,:), &		!Retirement
					VD(:,:,:,:), &		!Disabled
					VN(:,:,:,:,:,:,:), &	!Long-term Unemployed
					VW(:,:,:,:,:,:,:), &	!Working
					VU(:,:,:,:,:,:,:), &	!Unemployed
					V(:,:,:,:,:,:,:)	!Participant
		integer :: alloced 	! status indicator
	end type	
	
	!------------------------------------------------------------------
	!	b) pol_struct: aR, aD,aU,aN,aW,gapp,gwork, gapp_dif,gwork_dif
	!------------------------------------------------------------------
	type pol_struct
		
		real(dp), allocatable ::	gapp_dif(:,:,:,:,:,:,:), &
					gwork_dif(:,:,:,:,:,:,:) ! latent value of work/apply
		integer, allocatable ::	aR(:,:,:), aD(:,:,:,:), aU(:,:,:,:,:,:,:), &
					aN(:,:,:,:,:,:,:), aW(:,:,:,:,:,:,:)
		integer, allocatable ::	gapp(:,:,:,:,:,:,:), &

					gwork(:,:,:,:,:,:,:) !integer choice of apply/work
		integer :: alloced
	end type

	type moments_struct
		real(dp) :: work_coefs(Nk), di_coefs(Nk),ts_emp_coefs(nj+2)
		real(dp) :: di_rate(TT-1), work_rate(TT-1), accept_rate(TT-1) !by age
		integer :: alloced
		real(dp) :: work_cov_coefs(Nk,Nk),di_cov_coefs(Nk,Nk),ts_emp_cov_coefs(nj+2,nj+2)
		real(dp) :: s2

	end type 

	type hist_struct
		real(dp), allocatable :: work_dif_hist(:,:), app_dif_hist(:,:) !choose work or not, apply or not -- latent value
		real(dp), allocatable :: wage_hist(:,:) !realized wages
		real(dp), allocatable :: a_hist(:,:), al_hist(:,:), del_i(:)

		! a bunch of explanitory variables to be stacked on each other
		integer, allocatable :: status_hist(:,:) !status: W,U,N,D,R (1:5)

		integer, allocatable :: age_hist(:,:), born_hist(:,:), del_i_int(:)
		real(dp), allocatable :: occgrow_jt(:,:)
		real(dp), allocatable :: occshrink_jt(:,:)
		real(dp), allocatable :: occsize_jt(:,:)
		integer, allocatable :: j_i(:),z_jt_macro(:)
		integer, allocatable :: d_hist(:,:),al_int_hist(:,:)
		integer :: alloced
		integer :: drawn
	end type
	
	contains

	!------------------------------------------------------------------------
	!1)UI(e): Unemployment Insurance
	!----------------
	function UI(ein)
	!----------------
	
		real(dp), intent(in)	:: ein 
		real(dp) 		:: UI
		!I'm using a replacement rate of UIrr % for now, can be fancier
		UI = ein*UIrr

	end function
	!------------------------------------------------------------------------
	!2) DI (can be asset tested)
	!---------------------
	function SSDI(ein)
	!---------------------
	 
	
		real(dp), intent(in)	:: ein
		real(dp) 		:: SSDI
		!Follows Pistafferi & Low '12
		IF (ein<DItest1) then
			SSDI = 0.9*ein
		ELSEIF (ein<DItest2) then
			SSDI = 0.9*DItest1 + 0.32*(ein-DItest1)
		ELSEIF (ein<DItest3) then
			SSDI = 0.9*DItest1 + 0.32*(ein-DItest1)+0.15*(ein-DItest2)
		ELSE
			SSDI = 0.9*DItest1 + 0.32*(ein-DItest1)+0.15*(DItest3-DItest2)
		END IF

	end function
	!------------------------------------------------------------------------
	!3) Social Sec. Retirement	
	!--------------------
	function SSI(ein)
	!--------------------
	 
		real(dp), intent(in)	:: ein
		real(dp)			:: SSI
		
		!Follows Pistafferi & Low '15
		IF (ein<DItest1) then
			SSI = 0.9*ein
		ELSEIF (ein<DItest2) then
			SSI = 0.9*DItest1 + 0.32*(ein-DItest1)
		ELSEIF (ein<DItest3) then
			SSI = 0.9*DItest1 + 0.32*(ein-DItest1)+0.15*(ein-DItest2)
		ELSE
			SSI = 0.9*DItest1 + 0.32*(ein-DItest1)+0.15*(DItest3-DItest2)
		END IF

	end function
	!------------------------------------------------------------------------
	! 4) Acceptance probability
	!--------------------
	function xifun(idin,zin,itin)

		real(dp), intent(in):: zin
		integer, intent(in):: idin,itin
		real(dp) :: xifun
		!if(itin < 4) then
			xifun = 1._dp-(1._dp-xi(idin,itin))**(1._dp/tlen) !+ (1._dp - xi(idin,itin))*(xizcoef*zin)
		!else 
		!	xifun = xi(idin,itin) + (1._dp - xi(idin,itin))*(xizcoef*zin + 0.124)
		!endif
		
	end function

	!------------------------------------------------------------------------
	! 4) Utility Function
	!--------------------
	function util(cin,din,wkin)
	!--------------------
	 
		real(dp),intent(in)	:: cin
		integer, intent(in)	:: din, wkin
		real(dp)			:: util
		
		
		if ((wkin .gt. 1) .or. (din .gt. 1)) then
			if(gam> 1+1e-5_dp .or. gam < 1-1e-5_dp) then
				util = ((cin*dexp(theta*dble(din-1)+eta*dble(wkin-1)))**(1._dp-gam) )/(1._dp-gam)
			else 
				util = dlog(cin*dexp(theta*dble(din-1)+eta*dble(wkin-1)))
			endif
		else 
			if(gam> 1._dp+1e-5_dp .or. gam < 1_dp-1e-5_dp) then
				util = (cin**(1._dp-gam))/(1._dp-gam)
			else
				util = dlog(cin)
			endif
		end if
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

		Hearn = wgin/(tlen*dble(agegrid(tin))) + egrid(ein)*(1._dp-1._dp/(tlen*dble(agegrid(tin))))
	
	end function

	!------------------------------------------------------------------------
	! 6) Wage Function
	!--------------------
	function wage(biin,aiin,din,zin,tin)
	!--------------------
	
		real(dp), intent(in)	:: biin, aiin, zin
		integer, intent(in)	:: din, tin
		real(dp)			:: wage

		wage = dexp( biin*zin + aiin+wd(din)+wtau(tin) ) 

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
		enddo
		write(1,FMT) A(ri,c)
	enddo
	write(1,*) " "! trailing space
	close(1)

	end subroutine mat2csv

	!------------------------------------------------------------------------
	! 9) Write a Matrix to .csv
	!--------------------
	subroutine mati2csv(A,fname,append)
	!--------------------

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
		enddo
		write(1,FMT) A(ri,c)
	enddo
	write(1,*) " "! trailing space
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
	write(1,*) " "! trailing space
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
		write(1,*) " "! trailing space
		close(1)

	end subroutine veci2csv
	
	
	
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
	
		call DGETRF(n,n,invA,n,ipiv,info) ! first computes the LU factorization
	
		if (info /= 0) then
			print *, 'Matrix is singular.  Make it less singular before proceeding'
		endif
		call DGETRI(n,invA,n,ipiv,wk,n,info)
		if(info /=0) then
			print *, 'Matrix inversion failed, though it is not singular'
		endif
	end subroutine invmat
	
	!------------------------------------------------------------------------
	! 12) Run an OLS regression
	!------------------------------------------------------------------------
	subroutine OLS(XX,Y,coefs,cov_coef, hatsig2, status)
		real(dp), dimension(:,:), intent(in) :: XX
		real(dp), dimension(:), intent(in) :: Y
		real(dp), dimension(:), intent(out) :: coefs
		real(dp), dimension(:,:), intent(out) :: cov_coef
		real(dp), intent(out) :: hatsig2
		integer, intent(out) :: status

		integer :: nX, nY, nK, r,c
		real(dp), dimension(:,:), allocatable :: XpX,XpX_fac,XpX_inv
		real(dp), dimension(:), allocatable :: fitted,resids
		integer :: i
		
		external dgemm,dgemv
		external dpotrs,dpotrf,dpotri

		nK = size(XX, dim = 2)
		nX = size(XX, dim = 1)
		nY = size(Y)
		
		allocate(XpX(nK,nK))
		allocate(XpX_fac(nK,nK))
		allocate(XpX_inv(nK,nK))
		allocate(fitted(nX))
		allocate(resids(nX))
		coefs = 0.
		cov_coef = 0.
		
		
		if(nY /= nX ) then
			if(verbose>0 ) print *, 'size of X and Y in regression not compatible'
			status = 1
		else 
			XpX = 0.
			call dgemm('T', 'N', nK, nK, nX, 1._dp, XX, nX, XX, nX, 0., XpX, nK)
			XpX_fac = XpX
			call dpotrf('U',Nk,XpX_fac,Nk,status)
			if(status == 0) then
				! 3/ Multiply LHS of regression and solve it
				call dgemv('T', nX, nK, 1._dp, XX, nX, Y, 1, 0., coefs, 1)
				call dpotrs('U',nK,1,XpX_fac,nK,coefs,Nk,status)
			else 
				if(verbose >0 ) print *, "cannot factor XX"
			endif
		endif
		
		if(status == 0) then 
			XpX_inv = XpX_fac
			call dpotri('U',nK,XpX_inv,nK,status)
			fitted = 0.
			call dgemv('N', nX, nK, 1._dp, XX, nX, coefs, 1, 0., fitted, 1)
			resids = Y - fitted
			hatsig2 = 0.
			do i=1,nY
				hatsig2 = (resids(i)**2)/dble(nX-nK) + hatsig2
			enddo
			if(status == 0) cov_coef = hatsig2*XpX_inv
		endif
		
		deallocate(XpX,XpX_fac,XpX_inv,fitted,resids)
	
	end subroutine OLS



	subroutine alloc_econ(val_sol, pol_sol,hists_sim)

	! Structure to communicate everything
		type(val_struct) :: val_sol
		type(pol_struct) :: pol_sol
		type(hist_struct):: hists_sim

		!************************************************************************************************!
		! Allocate phat matrices
		!************************************************************************************************!
		! (disability extent, earn hist, assets)
		allocate(val_sol%VR(nd,ne,na), stat=val_sol%alloced)
		allocate(pol_sol%aR(nd,ne,na), stat=pol_sol%alloced)

		! (disability extent, earn hist, assets, age)
		allocate(val_sol%VD(nd,ne,na,TT), stat=val_sol%alloced)
		allocate(pol_sol%aD(nd,ne,na,TT-1), stat=pol_sol%alloced)

		! (occupation X ind exposure, ind disb. risk X ind. wage, disab. extent, earn hist, assets, agg shock, age)
		allocate(val_sol%VN(nj*nbi,ndi*nal,nd,ne,na,nz,TT), stat=val_sol%alloced)
		allocate(val_sol%VU(nj*nbi,ndi*nal,nd,ne,na,nz,TT), stat=val_sol%alloced)
		allocate(val_sol%VW(nj*nbi,ndi*nal,nd,ne,na,nz,TT), stat=val_sol%alloced)
		allocate(val_sol%V(nj*nbi,ndi*nal,nd,ne,na,nz,TT), stat=val_sol%alloced)
		allocate(pol_sol%aN(nj*nbi,ndi*nal,nd,ne,na,nz,TT-1), stat=pol_sol%alloced)
		allocate(pol_sol%aW(nj*nbi,ndi*nal,nd,ne,na,nz,TT-1), stat=pol_sol%alloced)
		allocate(pol_sol%aU(nj*nbi,ndi*nal,nd,ne,na,nz,TT-1), stat=pol_sol%alloced)
		allocate(pol_sol%gwork(nj*nbi,ndi*nal,nd,ne,na,nz,TT-1), stat=pol_sol%alloced)
		allocate(pol_sol%gapp(nj*nbi,ndi*nal,nd,ne,na,nz,TT-1), stat=pol_sol%alloced)

		allocate(pol_sol%gapp_dif(nj*nbi,ndi*nal,nd,ne,na,nz,TT), stat=pol_sol%alloced)
		allocate(pol_sol%gwork_dif(nj*nbi,ndi*nal,nd,ne,na,nz,TT), stat=pol_sol%alloced)

		allocate(hists_sim%wage_hist(Nsim,Tsim), stat=hists_sim%alloced)
		allocate(hists_sim%work_dif_hist(Nsim,Tsim), stat=hists_sim%alloced)
		allocate(hists_sim%app_dif_hist(Nsim,Tsim), stat=hists_sim%alloced)
		allocate(hists_sim%age_hist(Nsim,Tsim), stat=hists_sim%alloced)
		allocate(hists_sim%born_hist(Nsim,Tsim), stat=hists_sim%alloced)
		allocate(hists_sim%status_hist(Nsim,Tsim), stat=hists_sim%alloced)
		allocate(hists_sim%al_hist(Nsim,Tsim), stat=hists_sim%alloced)
		allocate(hists_sim%al_int_hist(Nsim,Tsim), stat=hists_sim%alloced)
		allocate(hists_sim%d_hist(Nsim,Tsim), stat=hists_sim%alloced)
		allocate(hists_sim%a_hist(Nsim,Tsim), stat=hists_sim%alloced)
		allocate(hists_sim%j_i(Nsim), stat=hists_sim%alloced)
		allocate(hists_sim%z_jt_macro(Tsim), stat=hists_sim%alloced)
		allocate(hists_sim%del_i_int(Nsim), stat=hists_sim%alloced)
		allocate(hists_sim%del_i(Nsim), stat=hists_sim%alloced)
		allocate(hists_sim%occgrow_jt(nj,Tsim), stat=hists_sim%alloced)
		allocate(hists_sim%occshrink_jt(nj,Tsim), stat=hists_sim%alloced)
		allocate(hists_sim%occsize_jt(nj,Tsim), stat=hists_sim%alloced)

	end subroutine

	subroutine dealloc_econ(val_sol,pol_sol,hists_sim)

	! Structure to communicate everything
		type(val_struct) :: val_sol
		type(pol_struct) :: pol_sol
		type(hist_struct):: hists_sim

		deallocate(pol_sol%aR,pol_sol%aD,pol_sol%aN, pol_sol%aU, pol_sol%aW,pol_sol%gwork, pol_sol%gapp)
		deallocate(pol_sol%gwork_dif,pol_sol%gapp_dif)
		deallocate(val_sol%VR,val_sol%VD,val_sol%VN,val_sol%VU,val_sol%VW,val_sol%V)
		deallocate(hists_sim%wage_hist,hists_sim%work_dif_hist)
		deallocate(hists_sim%app_dif_hist,hists_sim%status_hist,hists_sim%born_hist)
		deallocate(hists_sim%age_hist, hists_sim%al_hist, hists_sim%d_hist)
		deallocate(hists_sim%a_hist,hists_sim%j_i,hists_sim%del_i_int,hists_sim%del_i)
		deallocate(hists_sim%z_jt_macro,hists_sim%al_int_hist)
		deallocate(hists_sim%occgrow_jt, hists_sim%occshrink_jt )
		deallocate(hists_sim%occsize_jt)
		hists_sim%alloced = 0
	end subroutine

end module helper_funs

module model_data

	use V0para
	use helper_funs
	
	implicit none
	
	contains
	
	subroutine ts_employment(hists_sim,moments_sim)

		type(moments_struct)	:: moments_sim
		type(hist_struct)	:: hists_sim

		real(dp), allocatable :: emp_lagtimeconst(:,:), emp_t(:) ! will be a fraction of labor force in each occupation
		integer ij,it,nobs,yr, status

		allocate(emp_t((Tsim-1)*nj))
		allocate(emp_lagtimeconst((Tsim-1)*nj,2+nj))

		do it=2,Tsim
			do ij =1,nj
				emp_t((ij-1)*(Tsim-1) + it-1) = hists_sim%occsize_jt(ij,it)
				emp_lagtimeconst((ij-1)*(Tsim-1) + it-1,1) = hists_sim%occsize_jt(ij,it-1)
				emp_lagtimeconst((ij-1)*(Tsim-1) + it-1,1+ij) = dble(it-1)/tlen
				emp_lagtimeconst((ij-1)*(Tsim-1) + it-1,2+nj) = 1._dp
			enddo
		enddo

		call OLS(emp_lagtimeconst,emp_t,moments_sim%ts_emp_coefs,moments_sim%ts_emp_cov_coefs,moments_sim%s2,status)

		if(print_lev >=2 ) call vec2csv(moments_sim%ts_emp_coefs,"ts_emp_coefs.csv")

		deallocate(emp_lagtimeconst,emp_t)

	end subroutine ts_employment
	
	subroutine LPM_employment(hists_sim,moments_sim)
	
		type(moments_struct)	:: moments_sim
		type(hist_struct)	:: hists_sim

		real(dp), allocatable :: obsX_vars(:,:), work_dif_long(:)
		integer :: i, ij,it,age_hr,id, status, nobs, yr_stride
		real(dp) :: colmean, colvar
		logical :: badcoef(Nk) = .false.

		allocate(obsX_vars(Tsim*Nsim,Nk))
		allocate(work_dif_long(Tsim*Nsim))

		if(hists_sim%alloced /= 0) then
			if(verbose >= 1) print *, "not correctly passing hists_struct to LPM"
		endif

		yr_stride = itlen ! for doing year-length differences
		! build X
		obsX_vars = 0.
		Nobs = 0
		do i=1,Nsim
			do it=1,Tsim
				if( hists_sim%age_hist(i,it) <= TT-1 ) then
					Nobs = Nobs + 1
					work_dif_long( Nobs ) = dexp(hists_sim%work_dif_hist(i,it)*smthELPM )/(1._dp + dexp(hists_sim%work_dif_hist(i,it)*smthELPM) )
					do age_hr=1,TT-1
						if(hists_sim%age_hist(i,it) .eq. age_hr ) obsX_vars( Nobs, age_hr) = 1._dp
					enddo
					do id = 1,nd-1
						if(hists_sim%d_hist(i,it) .eq. id ) obsX_vars(Nobs, (TT-1)+id) = 1._dp
						! lead 1 period health status
						if(it <= Tsim - yr_stride) then
							if(hists_sim%d_hist(i,it+yr_stride) .eq. id ) obsX_vars(Nobs, (TT-1)+nd-1+id) = 1._dp
						endif
					enddo
					do ij = 1,nj
						if(hists_sim%j_i(i) == ij )  then
							if(it<=Tsim - yr_stride .and. it>= yr_stride+1) then
								obsX_vars(Nobs, (TT-1)+(nd-1)*2+1) = (hists_sim%occsize_jt(ij,it+yr_stride) - hists_sim%occsize_jt(ij,it-yr_stride)) &
													& /(hists_sim%occsize_jt(ij,it))
							endif
							if(it<=Tsim - 2*yr_stride .and. it>= 2*yr_stride+1) then
								obsX_vars(Nobs, (TT-1)+(nd-1)*2+2) = hists_sim%occsize_jt(ij,it+2*yr_stride) - hists_sim%occsize_jt(ij,it-2*yr_stride) &
													& /(hists_sim%occsize_jt(ij,it))
							endif
							
						endif
					enddo
				endif
			enddo
		enddo
		! the constant
		obsX_vars(:,Nk) = 1._dp
		moments_sim%work_coefs = 0.

		if(print_lev >=3 ) call mat2csv(obsX_vars,"obsX_vars.csv")
		!check that I don't have constants other than the constant rank < nk
		status = 0
		do ij=1,Nk-1
			colmean = sum(obsX_vars(1:Nobs,ij))/dble(Nobs)
			colvar = 0.
			do i=1,Nobs
				colvar = (obsX_vars(i,ij) - colmean)**2 + colvar
			enddo
			colvar = colvar / dble(Nobs)
			if(colvar < 1e-6) then
				badcoef(ij) = .true.
				status = status+1
				do i=1,Nobs
					call random_number(obsX_vars(i,ij))
					obsX_vars(i,ij) = (obsX_vars(i,ij)-0.5)/10.
				enddo
			endif
		enddo
		if(verbose >=2) print *, "bad cols of obsX", status 
		call OLS(obsX_vars(1:Nobs,:),work_dif_long(1:Nobs),moments_sim%work_coefs,moments_sim%work_cov_coefs,moments_sim%s2,status)

		do ij = 1,Nk
			if(badcoef(ij)) moments_sim%work_coefs(ij) = 0.
		enddo

		if(print_lev >=2 ) call vec2csv(moments_sim%work_coefs,"work_coefs.csv")
		if(print_lev >=2 ) call mat2csv(moments_sim%work_cov_coefs,"work_cov_coefs.csv")

		deallocate(obsX_vars,work_dif_long)
		
	end subroutine LPM_employment

	subroutine moments_compute(hists_sim,moments_sim)
	
		type(moments_struct) 	:: moments_sim
		type(hist_struct)	:: hists_sim
	
		integer :: i, ij,id,it,ial,st,si,age_hr,status_hr
		integer :: totage(TT),totD(TT),totW(TT),totst(Tsim),total(nal), tot3al(nal), &
				& tot3age(TT-1),totage_st(TT,Tsim)

		real(dp) :: dD_age(TT), dD_t(Tsim),a_age(TT),a_t(Tsim),alworkdif(nal),alappdif(nal), &
				& workdif_age(TT-1), appdif_age(TT-1), alD(nal), alD_age(nal,TT-1), &
				& status_Nt(5,Tsim)


		if(hists_sim%alloced /= 0) then
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
		
		dD_age 		= 0.
		dD_t 		= 0.
		a_age 		= 0.
		a_t 		= 0.
		alD		= 0.
		alD_age		= 0.
		appdif_age 	= 0.
		workdif_age 	= 0.
		alworkdif	= 0.
		alappdif	= 0.
		status_Nt 	= 0.
		
		do si = 1,Nsim
			do st = 1,Tsim
				if(hists_sim%status_hist(si,st) >0 .and. hists_sim%age_hist(si,st) >0 ) then
					age_hr = hists_sim%age_hist(si,st)
					totage_st(age_hr,st) = totage_st(age_hr,st) + 1
					! savings and disability by time
					a_t(st) = hists_sim%a_hist(si,st) + a_t(st)
					status_hr = hists_sim%status_hist(si,st)
					if(status_hr == 4) dD_t(st) = dD_t(st)+hists_sim%d_hist(si,st)
					
					status_Nt(status_hr,st) = 1._dp + status_Nt(status_hr,st)

					! disability by age and age X shock
					do it = 1,TT-1
						if(age_hr == it) then
							if(hists_sim%status_hist(si,st) == 1) totW(age_hr) = totW(age_hr) + 1
							if(hists_sim%status_hist(si,st) < 3) &
							&	workdif_age(age_hr) = workdif_age(age_hr) + hists_sim%work_dif_hist(si,st)

							if(hists_sim%status_hist(si,st) == 4) then
								totD(age_hr) = totD(age_hr) + 1
								dD_age(age_hr) = dD_age(age_hr)+hists_sim%d_hist(si,st)
								! associate this with its shock
								do ial = 1,nal
									if(  hists_sim%al_hist(si,st) <= alfgrid(ial)+2*epsilon(1._dp) &
									&	.and. hists_sim%al_hist(si,st) >= alfgrid(ial)-2*epsilon(1._dp) &
									&	.and. hists_sim%status_hist(si,st) == 4 ) &
									&	alD_age(ial,age_hr) = 1._dp + alD_age(ial,it)
								enddo
							elseif(hists_sim%status_hist(si,st) == 3) then
								appdif_age(it) = appdif_age(age_hr) + hists_sim%app_dif_hist(si,st)
								tot3age(age_hr) = tot3age(age_hr) + 1							
							endif
						endif
					enddo !it = 1,TT-1
					! assets by age
					do it=1,TT
						if(age_hr == it) then
							a_age(it) = hists_sim%a_hist(si,st) +a_age(it)
						endif
					enddo
					
					! disability and app choice by shock level
					do ial = 1,nal
						if( hists_sim%al_hist(si,st) <= alfgrid(ial)+2*epsilon(1._dp) &
						&	.and. hists_sim%al_hist(si,st) >= alfgrid(ial)-2*epsilon(1._dp)) then
							if(hists_sim%status_hist(si,st) <3) then
								!work choice:
								alworkdif(ial) = alworkdif(ial) + hists_sim%work_dif_hist(si,st)
								total(ial) = total(ial) + 1
							endif
							if(hists_sim%status_hist(si,st) == 3) then
								alappdif(ial) = alappdif(ial) + hists_sim%app_dif_hist(si,st)
								tot3al(ial) = tot3al(ial) + 1
							endif
							if(hists_sim%status_hist(si,st) == 4) &
								& alD(ial) = 1._dp + alD(ial)
						endif
					enddo
				endif !status(si,st) >0
			enddo!st = 1,Tsim
		enddo ! si=1,Nsim

		!just for convenience
		forall(it=1:TT) totage(it) = sum(totage_st(it,:))
		
		!work-dif app-dif distribution by age, shock
		forall(it=1:TT-1) appdif_age(it) = appdif_age(it)/dble(tot3age(it))
		forall(it=1:TT-1) workdif_age(it) = workdif_age(it)/dble(totage(it)-totD(it))
		forall(ial=1:nal) alappdif(ial) = alappdif(ial)/dble(tot3al(ial))
		forall(ial=1:nal) alworkdif(ial) = alworkdif(ial)/dble(total(ial))
		!disability distribution by shock and shock,age
		forall(ial=1:nal) alworkdif(ial) = alworkdif(ial)/dble(total(ial))
		forall(ial=1:nal) alD(ial) = dble(alD(ial))/dble(total(ial) + alD(ial))
		forall(it=1:TT-1,ial=1:nal) alD_age(ial,it) = alD_age(ial,it)/dble(total(ial))/dble(totage(it))
		! status distribution by age
		forall(it=1:TT-1) moments_sim%di_rate(it) = dble(totD(it))/dble(totage(it))
		forall(it=1:TT-1) moments_sim%work_rate(it) = dble(totW(it))/dble(totage(it))
		! asset distribution by age, time and disability status
		forall(it=1:TT ) a_age(it) = a_age(it)/dble(totage(it))
		forall(st=1:Tsim) a_t(st) = a_t(st)/dble(sum(totage_st(:,st)))

		forall(st=1:Tsim) status_Nt(:,st)= status_Nt(:,st)/sum(status_Nt(:,st))

		if(print_lev >= 2) then
			call veci2csv(totst,"pop_st.csv")
			call vec2csv(a_age,"a_age.csv")
			call vec2csv(a_t,"a_t.csv")
			call vec2csv(moments_sim%di_rate,"di_age.csv")
			call vec2csv(appdif_age, "appdif_age.csv")
			call vec2csv(moments_sim%work_rate,"work_age.csv")
			call vec2csv(workdif_age, "workdif_age.csv")
			call vec2csv(alD,"alD.csv")
			call mat2csv(alD_age,"alD_age.csv")
			call mat2csv(status_Nt,"status_Nt.csv")
		endif

		call LPM_employment(hists_sim,moments_sim)
		call ts_employment(hists_sim, moments_sim)
		
	end subroutine moments_compute


end module model_data

!**************************************************************************************************************!
!**************************************************************************************************************!
!						Solve the model						       !
!**************************************************************************************************************!
!**************************************************************************************************************!

module sol_val

	use V0para
	use helper_funs
! module with subroutines to solve the model and simulate data from it 

	implicit none

	integer ::  Vevals
	
	contains

	subroutine maxVR(id,ie,ia, VR0, iaa0, iaaA, apol,Vout)
		integer, intent(in) :: id,ie,ia
		integer, intent(in) :: iaa0, iaaA 
		integer, intent(out) :: apol
		real(dp), intent(out) :: Vout
		real(dp), intent(in) :: VR0(:,:,:)
		real(dp) :: Vtest1,Vtest2,chere, Vc1
		integer :: iaa, iw
		

		iw=1
		Vtest1 = -1e6
		apol = iaa0
		do iaa=iaa0,iaaA
			Vevals = Vevals +1
			chere = SSI(egrid(ie))+ R*agrid(ia) - agrid(iaa)
			if( chere .gt. 0.) then !ensure positive consumption
				Vc1 = beta*ptau(TT)*VR0(id,ie,iaa)
				Vtest2 = util(chere ,id,iw) + Vc1

				if(Vtest2 > Vtest1  .or. iaa .eq. iaa0 ) then !always replace on the first loop
					Vtest1 = Vtest2
					apol = iaa
				elseif( simp_concav .eqv. .true. ) then ! imposing concavity
					exit
				endif
			else!saved too much and negative consumtion
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
		real(dp) :: Vc1,chere,Vtest1,Vtest2,ein
		integer :: iw, iaa

		iw = 1 ! don't work
		Vc1 = beta*((1-ptau(it))*VD0(id,ie,iaa0,it+1)+ptau(it)*VD0(id,ie,iaa0,it))
		ein = egrid(ie)
		chere = SSDI(ein+R*agrid(ia)-agrid(iaa0))

		Vtest1 = -1e6
		apol = iaa0
		!Find Policy
		do iaa=iaa0,iaaA
			chere = SSDI(egrid(ie))+R*agrid(ia)-agrid(iaa)
			if(chere >0.) then
				Vc1 = beta*((1-ptau(it))*VD0(id,ie,iaa,it+1)+ptau(it)*VD0(id,ie,iaa,it))

				Vtest2 = util(chere,id,iw)+ Vc1
				if(Vtest2>Vtest1) then
					Vtest1 = Vtest2
					apol = iaa
				elseif(simp_concav .eqv. .true.) then
					exit
				endif
			else
				exit
			endif
		enddo	!iaa
		Vout = Vtest1
	end subroutine maxVD

	subroutine maxVU(ij,ibi,idi,ial,id,ie,ia,iz,it, VN0,V0, iaa0,iaaA,apol,Vout)
		integer, intent(in) :: ij,ibi,idi,ial,id,ie,ia,iz,it
		integer, intent(in) :: iaa0, iaaA 
		integer, intent(out) :: apol
		real(dp), intent(out) :: Vout
		real(dp), intent(in) :: VN0(:,:,:,:,:,:,:),V0(:,:,:,:,:,:,:)
		real(dp) :: Vc1,chere,Vtest1,Vtest2
		integer :: iw, iaa,iaai,izz

		iw = 1 ! don't work

		Vtest1 = -1e6 ! just a very bad number, does not really matter
		apol = iaa0
		do iaa=iaa0,iaaA
			chere = UI(egrid(ie))+R*agrid(ia)-agrid(iaa)
			if(chere>0.) then 
				Vtest2 = 0. !Continuation value if don't go on disability
				do izz = 1,nz	 !Loop over z'
				do iaai = 1,nal !Loop over alpha_i'
					Vc1 = (1._dp-ptau(it))*(pphi*VN0((ij-1)*nbi+ibi,(idi-1)*nal+iaai,id,ie,iaa,izz,it+1) &
						& 	+(1-pphi)*   V0((ij-1)*nbi+ibi,(idi-1)*nal+iaai,id,ie,iaa,izz,it+1) )  !Age and might go LTU
					Vc1 = ptau(it)*(pphi*     VN0((ij-1)*nbi+ibi,(idi-1)*nal+iaai,id,ie,iaa,izz,it) & 
						&	+(1-pphi)*   V0((ij-1)*nbi+ibi,(idi-1)*nal+iaai,id,ie,iaa,izz,it) ) + Vc1    !Don't age, maybe LTU
					Vtest2 = Vtest2 + beta*piz(iz,izz)*pialf(ial,iaai)*Vc1  !Probability of alpha_i X z_i draw 
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

	subroutine maxVN(ij,ibi,idi,ial,id,ie,ia,iz,it, VN0, VD0, V0, iaa0,iaaA,apol,gapp_pol,gapp_dif,Vout  )
		integer, intent(in) :: ij,ibi,idi,ial,id,ie,ia,iz,it
		integer, intent(in) :: iaa0, iaaA 
		integer, intent(out) :: apol,gapp_pol
		real(dp), intent(out) :: gapp_dif
		real(dp), intent(out) :: Vout
		real(dp), intent(in) :: VN0(:,:,:,:,:,:,:),VD0(:,:,:,:),V0(:,:,:,:,:,:,:)
		real(dp) :: Vc1,chere,Vtest1,Vtest2,Vapp,VNapp,smthV, maxVNV0
		integer :: iw, iaa,iaai,izz,aapp,aNapp

		iw = 1 ! not working
		!**********Value if apply for DI 
		Vtest1 = -1e6
		apol = iaa0
		do iaa = iaa0,iaaA
			chere = b+R*agrid(ia)-agrid(iaa)
			if(chere >0.) then
				Vtest2 = 0.
				!Continuation if apply for DI
				do izz = 1,nz	 !Loop over z'
				do iaai = 1,nal !Loop over alpha_i'
					Vc1 = (1-ptau(it))*((1-xi(id,it))* &
						& VN0((ij-1)*nbi+ibi,(idi-1)*nal+iaai,id,ie,iaa,izz,it+1)  &
						& +xi(id,it)*VD0(id,ie,iaa,it+1)) !Age and might go on DI
					Vc1 = Vc1 + ptau(it)*((1-xi(id,it))* &
						& VN0((ij-1)*nbi+ibi,(idi-1)*nal+iaai,id,ie,iaa,izz,it)  &
						& +xi(id,it)*VD0(id,ie,iaa,it))     !Don't age, might go on DI		
					Vtest2 = Vtest2 + beta*piz(iz,izz)*pialf(ial,iaai)*Vc1 
				enddo
				enddo
				Vtest2 = util(chere,id,iw) + Vtest2 &
					& - nu*dabs( VN0((ij-1)*nbi+ibi,(idi-1)*nal+ial,id,ie,iaa,iz,it) )
				
				if (Vtest2>Vtest1  .or. iaa .eq. iaa0) then	
					apol = iaa
					Vtest1 = Vtest2
				elseif(simp_concav .eqv. .true.) then
					exit
				endif
			!elseif(iaa<= iaa0 .and. Vtest1 <= -1e5 .and. apol == iaa0) then
			!	iaa = 1 !started too much saving, go back towards zero
			else
				exit
			endif
		enddo !iaa
		Vapp = Vtest1
		aapp = apol !agrid(apol)					

		if(Vapp <-1e5) then
			write(*,*) "ruh roh!"
			write(*,*) "Vapp, aapp: ", Vapp, aapp
			write(*,*) "VD: ",id,ie,iaa,it
			write(*,*) "VN: ",ij,ibi,idi,iaai,id,ie,iaa,izz,it
		endif

		!*******************************************
		!**************Value if do not apply for DI
		Vtest1 = -1e6
		apol = iaa0
		do iaa = iaa0,iaaA
			chere = b+R*agrid(ia)-agrid(iaa)
			if(chere >0.) then
				Vtest2 = 0.
				!Continuation if do not apply for DI
				do izz = 1,nz	 !Loop over z'
				do iaai = 1,nal !Loop over alpha_i'
				
					maxVNV0 = max(		 V0((ij-1)*nbi+ibi,(idi-1)*nal+iaai,id,ie,iaa,izz,it+1), &
							& 	VN0((ij-1)*nbi+ibi,(idi-1)*nal+iaai,id,ie,iaa,izz,it+1))
					Vc1 = (1-ptau(it))*((1-lrho)* &
							&	VN0((ij-1)*nbi+ibi,(idi-1)*nal +iaai,id,ie,iaa,izz,it+1) +lrho*maxVNV0) !Age and might go on DI
					maxVNV0 = max(		 V0((ij-1)*nbi+ibi,(idi-1)*nal+iaai,id,ie,iaa,izz,it), & 
							&	VN0((ij-1)*nbi+ibi,(idi-1)*nal+iaai,id,ie,iaa,izz,it))
					Vc1 = Vc1+ptau(it)*((1-lrho)* & 
							&	VN0((ij-1)*nbi+ibi,(idi-1)*nal +iaai,id,ie,iaa,izz,it) +lrho*maxVNV0)     !Don't age, might go on DI
					Vtest2 = Vtest2 + beta*piz(iz,izz)*pialf(ial,iaai)*Vc1 
				enddo
				enddo
				Vtest2 = Vtest2 + util(chere,id,iw)
				
				if( Vtest2 > Vtest1 .or. iaa .eq. iaa0) then
					apol = iaa 
					Vtest1 = Vtest2
				elseif(simp_concav .eqv. .true.) then
					exit
				endif
			!elseif(iaa .eq. iaa0 .and. Vtest1 <= -1e5 .and. apol .eq. iaa0) then
			!	iaa = 1 !started too much saving, go back towards zero
			else
				exit
			endif
		enddo !iaa
		Vnapp = Vtest1 					
		aNapp = apol !agrid(apol)
		if(Vnapp <-1e5 .and. verbose >0) then
			write(*,*) "ruh roh!"
			write(*,*) "Vnapp, aNapp: ", Vnapp, aNapp
			write(*,*) "VD: ",id,ie,iaa,it
			write(*,*) "VN: ",ij,ibi,idi,ial,id,ie,ia,iz,it
		endif

		!******************************************************
		!***************** Discrete choice for application
		!smthV = dexp(smthV0param*Vnapp)/( dexp(smthV0param*Vnapp) +dexp(smthV0param*Vapp) )
		!if( smthV .lt. 1e-5 .or. smthV .gt. 0.999999 .or. isnan(smthV)) then
			if( Vapp > Vnapp ) smthV =0.
			if(Vnapp > Vapp  ) smthV =1._dp
		!endif
		if (Vapp > Vnapp) then
			apol = aapp
			gapp_pol = 1
		else !Don't apply
			apol = aNapp
			gapp_pol = 0
		endif

		gapp_dif = Vapp - Vnapp
		if(verbose .gt. 4) print *, Vapp - Vnapp
		if(it>1) then
			Vout = smthV*Vnapp + (1._dp - smthV)*Vapp
		else
			Vout = smthV*Vnapp + (1._dp - smthV)*(eligY*Vapp + (1._dp- eligY)*Vnapp)
		endif

	end subroutine maxVN

	subroutine maxVW(ij,ibi,idi,ial,id,ie,ia,iz,it, VU, V0,wagehere,iee1,iee2,iee1wt,iaa0,iaaA,apol,gwork_pol,gwork_dif,Vout ,VWout )
		integer, intent(in) :: ij,ibi,idi,ial,id,ie,ia,iz,it
		integer, intent(in) :: iaa0, iaaA,iee1,iee2
		real(dp), intent(in) :: iee1wt,wagehere
		integer, intent(out) :: apol,gwork_pol
		real(dp), intent(out) :: gwork_dif
		real(dp), intent(out) :: Vout, VWout
		real(dp), intent(in) :: VU(:,:,:,:,:,:,:),V0(:,:,:,:,:,:,:)
		real(dp) :: Vc1,utilhere,chere,Vtest1,Vtest2,VWhere,VUhere,smthV, yL, yH
		integer :: iw, iaa,iaai,izz,idd

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
				do iaai = 1,nal  !Loop over alpha_i'
				do idd  = 1,nd	 !Loop over d'
					!Linearly interpolating on e'
					yL = (1-ptau(it))*V0((ij-1)*nbi+ibi,(idi-1)*nal+iaai,idd,iee1,iaa,izz,it+1) & 
						& +ptau(it)*V0((ij-1)*nbi+ibi,(idi-1)*nal+iaai,idd,iee1,iaa,izz,it)
					yH = (1-ptau(it))*V0((ij-1)*nbi+ibi,(idi-1)*nal+iaai,idd,iee2,iaa,izz,it+1) & 
						& +ptau(it)*V0((ij-1)*nbi+ibi,(idi-1)*nal+iaai,idd,iee2,iaa,izz,it)
					Vc1 = piz(iz,izz)*pialf(ial,iaai)*pid(id,idd,idi,it) &
						& * (yH*(1._dp - iee1wt) + yL*iee1wt )&
						& + Vc1
				enddo
				enddo
				enddo
				utilhere = util(chere,id,iw)
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
		VUhere = VU((ij-1)*nbi+ibi,(idi-1)*nal+ial,id,ie,ia,iz,it)
		
		!------------------------------------------------!
		!Calculate V with solved vals of VW and VU -- i.e. can quit into unemployment
		!------------------------------------------------!
		
		if (VWhere>VUhere) then
			gwork_pol = 1
		else
			gwork_pol = 0
		endif
		!smthV = dexp( smthV0param *VWhere  ) &
		!	& /( dexp(smthV0param * VWhere) + dexp(smthV0param * VUhere ) )
		!if (smthV <1e-5 .or. smthV>.999999 .or. isnan(smthV) ) then
			if(VWhere > VUhere)	smthV = 1._dp 
			if(VWhere < VUhere)	smthV = 0. 							
		!endif
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

		integer  :: i, j, ia, ie, id, it, ga,gw, anapp,aapp, apol, ibi, ial, ij , idi,  &
			    iee1, iee2, iz, iw,wo, iter,npara,ipara, iaa_k,iaa0,iaaA,iaN
		integer, dimension(5) :: maxer_i

		integer :: aa_l(na), aa_u(na), ia_o(na),aa_m(na)

		logical :: ptrsucces
		!************************************************************************************************!
		! Value Functions- Stack z-risk j and indiv. exposure beta_i
		!************************************************************************************************!
		real(dp)  	  	:: Vtest1, maxer_v, smthV,smthV0param, wagehere, iee1wt, gadif,gwdif
		
		real(dp), allocatable	:: maxer(:,:,:,:,:)
		real(dp), allocatable :: VR0(:,:,:), &			!Retirement
					VD0(:,:,:,:), &			!Disabled
					VN0(:,:,:,:,:,:,:), &	!Long-term Unemployed
					VW0(:,:,:,:,:,:,:), &	!Working
					VU0(:,:,:,:,:,:,:), &	!Unemployed
					V0(:,:,:,:,:,:,:)	!Participant
				
		real(dp), pointer ::	VR(:,:,:), &			!Retirement
					VD(:,:,:,:), &			!Disabled
					VN(:,:,:,:,:,:,:), &	!Long-term Unemployed
					VW(:,:,:,:,:,:,:), &	!Working
					VU(:,:,:,:,:,:,:), &	!Unemployed
					V(:,:,:,:,:,:,:)	!Participant
	
		real(dp), pointer ::	gapp_dif(:,:,:,:,:,:,:), gwork_dif(:,:,:,:,:,:,:) ! latent value of work/apply
	
		integer, pointer ::	aR(:,:,:), aD(:,:,:,:), aU(:,:,:,:,:,:,:), &
					aN(:,:,:,:,:,:,:), aW(:,:,:,:,:,:,:)
		integer, pointer ::	gapp(:,:,:,:,:,:,:), &
					gwork(:,:,:,:,:,:,:)
	
		!************************************************************************************************!
		! Other
		!************************************************************************************************!
			real(dp)	:: junk,summer, eprime, emin, emax, VWhere, V_m(na)
		!************************************************************************************************!
		
		!************************************************************************************************!
		! Allocate phat matrices
		!************************************************************************************************!
		! (disability extent, earn hist, assets)
		allocate(VR0(nd,ne,na))
		allocate(VD0(nd,ne,na,TT))
		allocate(VN0(nj*nbi,ndi*nal,nd,ne,na,nz,TT))
		allocate(VU0(nj*nbi,ndi*nal,nd,ne,na,nz,TT))
		allocate(VW0(nj*nbi,ndi*nal,nd,ne,na,nz,TT))
		allocate(V0(nj*nbi,ndi*nal,nd,ne,na,nz,TT))
		! there must be a way to use pointers, but it doesn't seem to work
		VR => val_funs%VR
		aR => pol_funs%aR
		VD => val_funs%VD
		aD => pol_funs%aD
		VN => val_funs%VN
		VU => val_funs%VU
		VW => val_funs%VW
		V =>  val_funs%V
		aN => pol_funs%aN
		aW => pol_funs%aW
		aU => pol_funs%aU
		gwork => pol_funs%gwork
		gapp => pol_funs%gapp

		gapp_dif => pol_funs%gapp_dif
		gwork_dif => pol_funs%gwork_dif

		allocate(maxer(na,nz,ne,nd,nal))
		emin = minval(egrid)
		emax = maxval(egrid)

		ptrsucces = associated(VR,val_funs%VR)
		
		simp_concav = .false.
		!************************************************************************************************!
		! Caculate things that are independent of occupation/person type
		!	1) Value of Retired:  VR(d,e,a)
		!	2) Value of Disabled: VD(d,e,a)
		
	!1) Calculate Value of Retired: VR(d,e,a)
		!d in{1,2,3}  : disability extent
		!e inR+       :	earnings index
		!a inR+	      : asset holdings
		
		Vevals = 0
		!VFI with good guess
		!Initialize
		iw=1
		do id=1,nd
		do ie=1,ne
		do ia=1,na
			VR0(id,ie,ia) = util(SSI(egrid(ie))+R*agrid(ia),id,iw)* (1._dp/(1._dp-beta*ptau(TT)))
		enddo
		enddo
		enddo
		if(print_lev >3) then
			i = 1
			call vec2csv(VR0(i,i,:),"VR0.csv",0)
		endif		
		iter = 1
		simp_concav = .true. ! use simple concavity here
		do while (iter<=maxiter)
			summer = 0
			id =1
		  	do ie=1,ne

				iaN =0
				ia = 1
				iaa0 = 1
				iaaA = na
				call maxVR(id,ie,ia, VR0, iaa0, iaaA, apol,Vtest1)
				VR(id,ie,ia) = Vtest1
				aR(id,ie,ia) = apol !agrid(apol)
				iaN = iaN+1
				ia_o(iaN) = ia
					
				ia = na
				iaa0 = aR(id,ie,1)
				iaaA = na
				call maxVR(id,ie,ia, VR0, iaa0, iaaA, apol,Vtest1)
				VR(id,ie,ia) = Vtest1
				aR(id,ie,ia) = apol !agrid(apol)
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
						iaa0 = aR(id,ie, aa_l(iaa_k-1) )
						iaaA = aR(id,ie, aa_u(iaa_k-1) )
						call maxVR(id,ie,ia, VR0, iaa0, iaaA, apol,Vtest1)
						VR(id,ie,ia) = Vtest1
						aR(id,ie,ia) = apol
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
					summer = summer+ (VR(id,ie,ia)-VR0(id,ie,ia))**2
				enddo
				VR0(id,ie,:) = VR(id,ie,:)
			enddo !ie
			if (summer < Vtol) then
				exit	!Converged				
			else
				iter=iter+1
				if(print_lev >3) then
					i = 1
					call veci2csv(aR(i,i,:),"aR.csv",0)
					call vec2csv(VR(i,i,:),"VR.csv",0)
				endif
			endif
		enddo ! iteration iter

		if(summer >= Vtol) then
			print *, "VR did not converge"
		endif 

		i = 1
		do id =2,nd
		do ie =1,ne
		do ia =1,na
			VR(id,ie,ia) = VR(i,ie,ia)*(dexp(theta*dble(id-1)))**(1-gam)
			aR(id,ie,ia) = aR(i,ie,ia)
		enddo
		enddo
		enddo

		if (print_lev > 2) then
			wo =0
			do id=1,nd
			do ie=1,ne
				call veci2csv(aR(i,i,:),"aR.csv",wo)
				call vec2csv(VR(i,i,:),"VR.csv",wo)
				if(wo .eq. 0) wo = 1
			enddo
			enddo
		endif

		!----------------------------------------------------------!
		!Set value at t=TT to be VR in all other V-functions
		!----------------------------------------------------------!
		VD(:,:,:,TT) = VR
		VD0(:,:,:,TT) = VR
		VD0(:,:,:,TT-1) = VD0(:,:,:,TT) ! this is necessary for initialization given stochastic aging

	!----------------------------------------------------------------!
	!2) Calculate Value of Disabled: VD(d,e,a,t)	 
		!d in{1,2,3}  	   :disability extent
		!e inR+       	   :earnings index
		!a inR+	      	   :asset holdings
		!t in[1,2...TT-1]  :age

		simp_concav = .true.
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
						summer = 0

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
							summer = summer+ (VD(id,ie,ia,it)-VD0(id,ie,ia,it))**2
						enddo

						if(print_lev >3) then
							wo =0
							call veci2csv(aD(id,ie,:,it),"aD.csv",wo)
							call vec2csv(VD(id,ie,:,it),"VD.csv",wo)
						endif
						if (summer < Vtol) then
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
			VW (:,:,id,ie,ia,:,TT) = VR(id,ie,ia)
			VW0(:,:,id,ie,ia,:,TT) = VR(id,ie,ia)
			VN (:,:,id,ie,ia,:,TT) = VR(id,ie,ia)
			VN0(:,:,id,ie,ia,:,TT) = VR(id,ie,ia)
			VU (:,:,id,ie,ia,:,TT) = VR(id,ie,ia)
			VU0(:,:,id,ie,ia,:,TT) = VR(id,ie,ia)
			V  (:,:,id,ie,ia,:,TT) = VR(id,ie,ia)
			V0 (:,:,id,ie,ia,:,TT) = VR(id,ie,ia)	   
		enddo
		enddo
		enddo


	!initial guess for the value functions
		do ij = 1,nj
	! And betas
		do ibi = 1,nbi 
	! And individual disability type
		do idi = 1,ndi
	!************************************************************************************************!
		do it=TT-1,1,-1
			!----Initialize---!
			do ial=1,nal
			do id =1,nd
			do ie =1,ne
			do iz =1,nz
			do ia =1,na
				
			!Guess once, then use next period same occupation/beta as guess
			! for it = 1, should be TT-1+1 =TT -> VU,Vw,VN = VR
				VW0((ij-1)*nbi+ibi,(idi-1)*nal+ial,id,ie,ia,iz,it) = VW((ij-1)*nbi+ibi,(idi-1)*nal+ial,id,ie,ia,iz,it+1)
				VU0((ij-1)*nbi+ibi,(idi-1)*nal+ial,id,ie,ia,iz,it) = VU((ij-1)*nbi+ibi,(idi-1)*nal+ial,id,ie,ia,iz,it+1)
				VN0((ij-1)*nbi+ibi,(idi-1)*nal +ial,id,ie,ia,iz,it) = VN((ij-1)*nbi+ibi,(idi-1)*nal+ial,id,ie,ia,iz,it+1)
				V0((ij-1)*nbi+ibi,(idi-1)*nal+ial,id,ie,ia,iz,it)= VW0((ij-1)*nbi+ibi,(idi-1)*nal+ial,id,ie,ia,iz,it)

			enddo	!ia
			enddo	!iz
			enddo	!ie
			enddo 	!id
			enddo	!ial   
		enddo
		enddo
		enddo
		enddo
		
		! Begin loop over occupations
		do ij = 1,nj
		! And betas
		do ibi = 1,nbi 
		! And individual disability type
		do idi = 1,ndi

		!************************************************************************************************!
			!Work Backwards TT-1,TT-2...1!
		do it=TT-1,1,-1
			!----Initialize---!
			do ial=1,nal
			do id =1,nd
			do ie =1,ne
			do iz =1,nz
			do ia =1,na
				
			!Guess once, then use next period same occupation/beta as guess
			! for it = 1, should be TT-1+1 =TT -> VU,Vw,VN = VR
				VW0((ij-1)*nbi+ibi,(idi-1)*nal+ial,id,ie,ia,iz,it) = VW((ij-1)*nbi+ibi,(idi-1)*nal+ial,id,ie,ia,iz,it+1)
				VU0((ij-1)*nbi+ibi,(idi-1)*nal+ial,id,ie,ia,iz,it) = VU((ij-1)*nbi+ibi,(idi-1)*nal+ial,id,ie,ia,iz,it+1)
				VN0((ij-1)*nbi+ibi,(idi-1)*nal +ial,id,ie,ia,iz,it) = VN((ij-1)*nbi+ibi,(idi-1)*nal+ial,id,ie,ia,iz,it+1)
				V0((ij-1)*nbi+ibi,(idi-1)*nal+ial,id,ie,ia,iz,it)= VW0((ij-1)*nbi+ibi,(idi-1)*nal+ial,id,ie,ia,iz,it)

			enddo	!ia
			enddo	!iz
			enddo	!ie
			enddo 	!id
			enddo	!ial   

			!***********************************************************************************************************
			!Loop over V=max(VU,VW)	
			iter=1
			smthV0param  =1._dp ! will tighten this down
			do while (iter<=maxiter)
				maxer = 0.
				summer = 0.	!Use to calc |V-V0|<eps
				! lots f printing every 100 iterations (mod 1 because Fortran is designed by idiots using base-1 indexing)
				!if(mod(iter,3) .eq. 0) then
					simp_concav = .false.
				!else
				!	simp_concav = .true.
				!endif
			!------------------------------------------------!
			!Solve VU given guesses on VW, VN, VU and implied V
			!------------------------------------------------!
				summer = 0.
				wo = 0
				npara = nal*nd*ne*nz
			!$OMP  parallel do reduction(+:summer)&
			!$OMP& private(ial,id,ie,iz,iw,apol,iaa0,iaaA,ia,ipara,Vtest1,aa_m,V_m,aa_l,aa_u,iaa_k,iaN,ia_o)
				do ipara = 1,npara
					iz = mod(ipara-1,nz)+1
					ie = mod(ipara-1,nz*ne)/nz + 1
					id = mod(ipara-1,nz*ne*nd)/(nz*ne) +1
					ial= mod(ipara-1,nz*ne*nd*nal)/(nz*ne*nd)+1

					! search over ia
						iaN =0
						ia = 1
						iaa0 = 1
						iaaA = na
						call maxVU(ij,ibi,idi,ial,id,ie,ia,iz,it, VN0,V0, iaa0,iaaA,apol,Vtest1)
						V_m(ia) = Vtest1
						aa_m(ia) = apol !agrid(apol)
						iaN = iaN+1
						ia_o(iaN) = ia
					
						ia = na
						iaa0 = aa_m(1)
						iaaA = na
						call maxVU(ij,ibi,idi,ial,id,ie,ia,iz,it, VN0,V0, iaa0,iaaA,apol,Vtest1)
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
								call maxVU(ij,ibi,idi,ial,id,ie,ia,iz,it, VN0,V0, iaa0,iaaA,apol,Vtest1)
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
						
						aU((ij-1)*nbi+ibi,(idi-1)*nal+ial,id,ie,:,iz,it) = aa_m
						VU((ij-1)*nbi+ibi,(idi-1)*nal+ial,id,ie,:,iz,it) = V_m
							
						if((iz>nz .or. ie>ne .or. id>nd .or. ial>nal) .and. verbose > 2) then
							print *, "ipara is not working right", iz,ie,id,ial
						endif
						Vtest1 = sum( (VU((ij-1)*nbi+ibi,(idi-1)*nal+ial,id,ie,:,iz,it) &
							& - VU0((ij-1)*nbi+ibi,(idi-1)*nal+ial,id,ie,:,iz,it))**2)
						summer = Vtest1 + summer
						if(print_lev >3) then
							call veci2csv(ia_o,"ia_o_VU.csv",wo)
							if(wo == 0) wo =1
						endif
						
				enddo !ipara
			!$OMP END PARALLEL do
				!update VU0
				do ial=1,nal	!Loop over alpha (ai)
				do id=1,nd	!Loop over disability index
				do ie=1,ne	!Loop over earnings index
				do iz=1,nz	!Loop over TFP
					do ia =1,na
						VU0((ij-1)*nbi+ibi,(idi-1)*nal+ial,id,ie,ia,iz,it) = VU((ij-1)*nbi+ibi,(idi-1)*nal+ial,id,ie,ia,iz,it)
					enddo	!ia
				enddo !ial
				enddo !id
				enddo !ie
				enddo !iz
				if (print_lev > 3) then
					wo = 0
					do ial=1,nal	!Loop over alpha (ai)
					do id=1,nd	!Loop over disability index
					do ie=1,ne	!Loop over earnings index
					do iz=1,nz	!Loop over TFP
						call veci2csv(aU((ij-1)*nbi+ibi,(idi-1)*nal+ial,id,ie,:,iz,it),"aU.csv",wo)
						call vec2csv(VU((ij-1)*nbi+ibi,(idi-1)*nal+ial,id,ie,:,iz,it),"VU.csv",wo)
						if(wo == 0) wo = 1 
					enddo
					enddo
					enddo
					enddo
				endif


			!------------------------------------------------!
			!Solve VN given guesses on VW, VN, and implied V
			!------------------------------------------------! 
				summer = 0.
				wo = 0
				aa_u = 0
				npara = nal*nd*ne*nz
			!$OMP  parallel do reduction(+:summer)&
			!$OMP& private(ipara,ial,id,ie,iz,apol,ga,gadif,ia,iaa0,iaaA,iaa_k,aa_l,aa_u,aa_m,V_m,Vtest1,iaN,ia_o) 
				do ipara = 1,npara
					iz = mod(ipara-1,nz)+1
					ie = mod(ipara-1,nz*ne)/nz + 1
					id = mod(ipara-1,nz*ne*nd)/(nz*ne) +1
					ial= mod(ipara-1,nz*ne*nd*nal)/(nz*ne*nd)+1

					!----------------------------------------------------------------
					!Loop over current state: assets
					iaN=0
					ia = 1
					iaa0 = 1
					iaaA = na
					call maxVN(ij,ibi,idi,ial,id,ie,ia,iz,it, VN0, VD0,V0,iaa0,iaaA,apol,ga,gadif,Vtest1)
					aa_m(ia) = apol
					V_m(ia) = Vtest1
					gapp((ij-1)*nbi+ibi,(idi-1)*nal+ial,id,ie,ia,iz,it) = ga
					gapp_dif((ij-1)*nbi+ibi,(idi-1)*nal+ial,id,ie,ia,iz,it) = gadif
					iaN = iaN+1
					ia_o(iaN) = ia
					

					ia = na
					iaa0 = aa_m(1)
					iaaA = na
					call maxVN(ij,ibi,idi,ial,id,ie,ia,iz,it, VN0, VD0,V0, iaa0,iaaA,apol,ga,gadif,Vtest1)
					V_m(ia) = Vtest1
					aa_m(ia) = apol !agrid(apol)
					gapp((ij-1)*nbi+ibi,(idi-1)*nal+ial,id,ie,ia,iz,it) = ga
					gapp_dif((ij-1)*nbi+ibi,(idi-1)*nal+ial,id,ie,ia,iz,it) = gadif
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
							call maxVN(ij,ibi,idi,ial,id,ie,ia,iz,it, VN0,VD0,V0, iaa0,iaaA,apol,ga,gadif,Vtest1)
							V_m(ia) = Vtest1
							aa_m(ia) = apol !agrid(apol)
							gapp((ij-1)*nbi+ibi,(idi-1)*nal+ial,id,ie,ia,iz,it) = ga
							gapp_dif((ij-1)*nbi+ibi,(idi-1)*nal+ial,id,ie,ia,iz,it) = gadif
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

					aN((ij-1)*nbi+ibi,(idi-1)*nal+ial,id,ie,:,iz,it) = aa_m
					VN((ij-1)*nbi+ibi,(idi-1)*nal+ial,id,ie,:,iz,it) = V_m
					
					Vtest1 = sum((VN((ij-1)*nbi+ibi,(idi-1)*nal+ial,id,ie,:,iz,it) &
						& - VN0((ij-1)*nbi+ibi,(idi-1)*nal+ial,id,ie,:,iz,it))**2)
					summer = Vtest1 + summer
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
					do ial=1,nal	!Loop over alpha (ai)
					do ie=1,ne	!Loop over earnings index
					do iz=1,nz	!Loop over TFP
						! matrix in disability index and assets
						call mat2csv(VN((ij-1)*nbi+ibi,(idi-1)*nal+ial,:,ie,:,iz,it) ,"VN_it.csv",wo)
						call mati2csv(aN((ij-1)*nbi+ibi,(idi-1)*nal+ial,:,ie,:,iz,it) ,"aN_it.csv",wo)
						call mat2csv(VN((ij-1)*nbi+ibi,(idi-1)*nal+ial,:,ie,:,iz,it) ,"VU_it.csv",wo)
						call mati2csv(aN((ij-1)*nbi+ibi,(idi-1)*nal+ial,:,ie,:,iz,it) ,"aU_it.csv",wo)
						call mati2csv(gapp((ij-1)*nbi+ibi,(idi-1)*nal+ial,:,ie,:,iz,it) ,"gapp_it.csv",wo)
						call mat2csv(gapp_dif((ij-1)*nbi+ibi,(idi-1)*nal+ial,:,ie,:,iz,it) ,"gapp_dif_it.csv",wo)

						if(wo == 0 ) wo =1		  		
					enddo !iz 
					enddo !id 
					enddo !ial 
				endif
			
				!update VN0
				do ial=1,nal	!Loop over alpha (ai)
				do id=1,nd	!Loop over disability index
				do ie=1,ne	!Loop over earnings index
				do iz=1,nz	!Loop over TFP
					do ia =1,na
						VN0((ij-1)*nbi+ibi,(idi-1)*nal+ial,id,ie,ia,iz,it) = VN((ij-1)*nbi+ibi,(idi-1)*nal+ial,id,ie,ia,iz,it)
					enddo	!ia
				enddo !ial
				enddo !id
				enddo !ie
				enddo !iz

				!------------------------------------------------!
				!Solve VW given guesses on VW, VN, and implied V
				!------------------------------------------------!
				summer = 0.
				wo = 0
				npara = nal*nd*ne*nz
			!$OMP   parallel do reduction(+:summer) &
			!$OMP & private(ipara,ial,id,ie,iz,apol,eprime,wagehere,iee1,iee2,iee1wt,ia,iaa0,iaaA,aa_l,aa_u,iaa_k,ia_o,iaN,Vtest1,VWhere,gwdif,gw) 
				do ipara = 1,npara
					iz = mod(ipara-1,nz)+1
					ie = mod(ipara-1,nz*ne)/nz + 1
					id = mod(ipara-1,nz*ne*nd)/(nz*ne) +1
					ial= mod(ipara-1,nz*ne*nd*nal)/(nz*ne*nd)+1

					!Earnings evolution independent of choices.
					wagehere = wage(beti(ibi),alfgrid(ial),id,zgrid(iz,ij),it)
					eprime = Hearn(it,ie,wagehere)
					!linear interpolate for the portion that blocks off bounds on assets
					if(eprime > emin .and. eprime < emax) then  ! this should be the same as if(eprime > minval(egrid) .and. eprime < maxval(egrid))
						iee1 = ne
						do while( eprime < egrid(iee1) .and. iee1>= 1)
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
					iaN = 0
					ia = 1
					iaa0 = 1
					iaaA = na
					call maxVW(ij,ibi,idi,ial,id,ie,ia,iz,it, VU, V0,wagehere,iee1,iee2,iee1wt, &
						& iaa0,iaaA,apol,gw,gwdif,Vtest1 ,VWhere )
					V	((ij-1)*nbi+ibi,(idi-1)*nal+ial,id,ie,ia,iz,it) = Vtest1
					VW	((ij-1)*nbi+ibi,(idi-1)*nal+ial,id,ie,ia,iz,it) = VWhere
					gwork	((ij-1)*nbi+ibi,(idi-1)*nal+ial,id,ie,ia,iz,it) = gw
					gwork_dif((ij-1)*nbi+ibi,(idi-1)*nal+ial,id,ie,ia,iz,it) = gwdif
					aW	((ij-1)*nbi+ibi,(idi-1)*nal+ial,id,ie,ia,iz,it) = apol
					iaN = iaN+1
					ia_o(iaN) = ia
					

					ia = na
					iaa0 = aW((ij-1)*nbi+ibi,(idi-1)*nal+ial,id,ie,1,iz,it)
					iaaA = na
					call maxVW(ij,ibi,idi,ial,id,ie,ia,iz,it, VU, V0,wagehere,iee1,iee2,iee1wt, &
						& iaa0,iaaA,apol,gw,gwdif,Vtest1 ,VWhere )
					V	((ij-1)*nbi+ibi,(idi-1)*nal+ial,id,ie,ia,iz,it) = Vtest1
					VW	((ij-1)*nbi+ibi,(idi-1)*nal+ial,id,ie,ia,iz,it) = VWhere
					gwork	((ij-1)*nbi+ibi,(idi-1)*nal+ial,id,ie,ia,iz,it) = gw
					gwork_dif((ij-1)*nbi+ibi,(idi-1)*nal+ial,id,ie,ia,iz,it) = gwdif
					aW	((ij-1)*nbi+ibi,(idi-1)*nal+ial,id,ie,ia,iz,it) = apol
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
							iaa0 = aW((ij-1)*nbi+ibi,(idi-1)*nal+ial,id,ie,  aa_l(iaa_k-1)  ,iz,it)
							iaaA = aW((ij-1)*nbi+ibi,(idi-1)*nal+ial,id,ie,  aa_u(iaa_k-1)  ,iz,it)
							call maxVW(ij,ibi,idi,ial,id,ie,ia,iz,it, VU, V0,wagehere,iee1,iee2,iee1wt, &
								& iaa0,iaaA,apol,gw,gwdif,Vtest1 ,VWhere )
							V	((ij-1)*nbi+ibi,(idi-1)*nal+ial,id,ie,ia,iz,it) = Vtest1
							VW	((ij-1)*nbi+ibi,(idi-1)*nal+ial,id,ie,ia,iz,it) = VWhere
							gwork	((ij-1)*nbi+ibi,(idi-1)*nal+ial,id,ie,ia,iz,it) = gw
							gwork_dif((ij-1)*nbi+ibi,(idi-1)*nal+ial,id,ie,ia,iz,it)= gwdif
							aW	((ij-1)*nbi+ibi,(idi-1)*nal+ial,id,ie,ia,iz,it) = apol
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
					
					Vtest1 = sum((V((ij-1)*nbi+ibi,(idi-1)*nal+ial,id,ie,:,iz,it) &
						& - V0((ij-1)*nbi+ibi,(idi-1)*nal+ial,id,ie,:,iz,it))**2)

					summer = Vtest1	+ summer

					do ia=1,na
						maxer(ia,iz,ie,id,ial) = (V((ij-1)*nbi+ibi,(idi-1)*nal+ial,id,ie,ia,iz,it)-V0((ij-1)*nbi+ibi,(idi-1)*nal+ial,id,ie,ia,iz,it))**2
					enddo
					if(print_lev >3) then
						call veci2csv(ia_o,"ia_o_VW.csv", wo)
						if(wo == 0) wo = 1
					endif
	!				enddo !iz  enddo !ie   	enddo !id  	enddo !ial
				enddo
			!$OMP  END PARALLEL do
			
				maxer_v = maxval(maxer)
				maxer_i = maxloc(maxer)

				!update VW0, V0
				do ial=1,nal	!Loop over alpha (ai)
				do id=1,nd	!Loop over disability index
				do ie=1,ne	!Loop over earnings index
				do iz=1,nz	!Loop over TFP
				do ia =1,na
					VW0((ij-1)*nbi+ibi,(idi-1)*nal+ial,id,ie,ia,iz,it) = VW((ij-1)*nbi+ibi,(idi-1)*nal+ial,id,ie,ia,iz,it)
					V0 ((ij-1)*nbi+ibi,(idi-1)*nal+ial,id,ie,ia,iz,it) =  V((ij-1)*nbi+ibi,(idi-1)*nal+ial,id,ie,ia,iz,it)
				enddo !ia		  	
				enddo !ial
				enddo !id
				enddo !ie
				enddo !iz

								
				if (print_lev >2) then
					wo = 0
					do ial=1,nal	!Loop over alpha (ai)
					do ie=1,ne	!Loop over earnings index
					do iz=1,nz	!Loop over TFP
						! matrix in disability and assets
						call mat2csv(VW((ij-1)*nbi+ibi,(idi-1)*nal+ial,:,ie,:,iz,it) ,"VW_it.csv",wo)
						call mati2csv(aW((ij-1)*nbi+ibi,(idi-1)*nal+ial,:,ie,:,iz,it) ,"aW_it.csv",wo)
						call mati2csv(gwork((ij-1)*nbi+ibi,(idi-1)*nal+ial,:,ie,:,iz,it) ,"gwork_it.csv",wo)
						call mat2csv(gwork_dif((ij-1)*nbi+ibi,(idi-1)*nal+ial,:,ie,:,iz,it) ,"gwork_idf_it.csv",wo)
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
				if(verbose > 3 .or. (verbose >1 .and. mod(iter,100).eq. 0)) then
					write(*,*) summer, iter, ij, ibi, idi, it
					write(*,*) maxer_v, maxer_i(1), maxer_i(2), maxer_i(3), maxer_i(4), maxer_i(5)
				endif
				if (summer < Vtol ) then
					exit !Converged
				endif

				iter=iter+1
				smthV0param = smthV0param*1.5_dp !tighten up the discrete choice
			enddo	!iter: V-iter loop
	!WRITE(*,*) ij, ibi, idi, it
		enddo	!t loop, going backwards

		enddo	!idi
		enddo	!ibi
		enddo	!ij

			
		! this plots work-rest and di application on the cross product of alphai and deltai and di
		if(print_lev >1) then
			ibi = 1
			ij  = 1
			wo  = 0

			do id  = 1,nd
			do ie  = 1,ne
				call mati2csv(aD(id,ie,:,:),"aD.csv",wo)
				call mat2csv (VD(id,ie,:,:),"VD.csv",wo)

				call veci2csv(aR(id,ie,:),"aR.csv",wo)
				call vec2csv (VR(id,ie,:),"VR.csv",wo)
				if(wo == 0) wo =1
			enddo
			enddo

			wo = 0
			do idi=1,ndi	!loop over delta(idi)
			do ial=1,nal	!Loop over alpha (al)
			do ie=1,ne	!Loop over earnings index
			do iz=1,nz	!Loop over TFP
				do it = TT-1,1,-1
					! matrix in disability and assets
					call mat2csv(V((ij-1)*nbi+ibi,(idi-1)*nal+ial,:,ie,:,iz,it) ,"V.csv",wo)

					call mat2csv(VW((ij-1)*nbi+ibi,(idi-1)*nal+ial,:,ie,:,iz,it) ,"VW.csv",wo)
					call mati2csv(aW((ij-1)*nbi+ibi,(idi-1)*nal+ial,:,ie,:,iz,it) ,"aW.csv",wo)

					call mat2csv(VU((ij-1)*nbi+ibi,(idi-1)*nal+ial,:,ie,:,iz,it) ,"VU.csv",wo)
					call mati2csv(aU((ij-1)*nbi+ibi,(idi-1)*nal+ial,:,ie,:,iz,it) ,"aU.csv",wo)

					call mat2csv(VN((ij-1)*nbi+ibi,(idi-1)*nal+ial,:,ie,:,iz,it) ,"VN.csv",wo)
					call mati2csv(aN((ij-1)*nbi+ibi,(idi-1)*nal+ial,:,ie,:,iz,it) ,"aN.csv",wo)

					call mati2csv(gwork((ij-1)*nbi+ibi,(idi-1)*nal+ial,:,ie,:,iz,it) ,"gwork.csv",wo)
					call mat2csv(gwork_dif((ij-1)*nbi+ibi,(idi-1)*nal+ial,:,ie,:,iz,it) ,"gwork_dif.csv",wo)

					call mati2csv(gapp((ij-1)*nbi+ibi,(idi-1)*nal+ial,:,ie,:,iz,it) ,"gapp.csv",wo)
					call mat2csv(gapp_dif((ij-1)*nbi+ibi,(idi-1)*nal+ial,:,ie,:,iz,it) ,"gapp_dif.csv",wo)

					if(wo==0) wo =1
				enddo !it
			enddo !iz 
			enddo !ie 
			enddo !ial 	
			enddo !idi
		endif

		deallocate(maxer)
		deallocate(VR0,VD0,VN0,VU0,VW0,V0)
	!		deallocate(VR,VD,VN,VU,VW,V)
	!		deallocate(aR,aD,aN,aW,aU,gwork,gapp,gapp_dif,gwork_dif)

	end subroutine sol 
end module sol_val


!**************************************************************************************************************!
!**************************************************************************************************************!
!						Simulate from solution					       !
!**************************************************************************************************************!
!**************************************************************************************************************!
module sim_hists
	use V0para
	use helper_funs
! module with subroutines to solve the model and simulate data from it 

	implicit none
	
	contains


	subroutine draw_deli(del_i, del_i_int, seed0, success)
	! draws depreciation rates and indices on the delta grid (i.e at the discrete values)
		implicit none

		integer, intent(in) :: seed0
		integer, intent(out) :: success
		real(dp), dimension(:) :: del_i
		integer, dimension(:) :: del_i_int
		integer :: ss, Nsim, di_int,m,i,ij,idi
		real(dp) :: delgridL, delgridH,delgrid_i
		real(dp) :: delwtH,delwtL
		integer, dimension(100) :: bdayseed

		call random_seed(size = ss)
		forall(m=1:ss) bdayseed(m) = (m-1)*100 + seed0
		call random_seed(put = bdayseed(1:ss) )

		Nsim = size(del_i)

		delwt	 = 1._dp/dble(ndi) ! initialize with equal weight
		if(del_by_occ .eqv. .true.) then ! give weight according to mean delta by occ
			delgridL = 0.
			do i=1,ndi/2
				delgridL = delgrid(i)/dble(ndi/2) + delgridL
			enddo
			delgridH = 0.
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
			do ij=1,nj
				delwt(:,ij) = 1._dp/dble(idi)
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
			call rand_num_closed(delgrid_i) ! draw uniform on 0,1
			delgrid_i = delgrid_i*(delgridH-delgridL) + delgridL !change domain of uniform
			di_int = finder(delgrid,delgrid_i)
			! round up or down:
			if(di_int < ndi) then
				if( (delgrid_i - delgrid(di_int))/(delgrid(di_int+1) - delgrid(di_int)) > 0.5 ) di_int = di_int + 1
			endif
			di_int = max(min(di_int,nd),1)
			del_i(i) = delgrid_i
			del_i_int(i) = di_int
		enddo
		success = 0
	end subroutine draw_deli
	
	subroutine draw_status_innov(status_it_innov, seed0, success)
	! draws innovations to d, will be used if working and relevant
		implicit none

		integer, intent(in) :: seed0
		integer, intent(out) :: success
		real(dp), dimension(:,:) :: status_it_innov
		integer :: ss, m,i,it
		real(dp) :: s_innov
		integer, dimension(100) :: bdayseed

		call random_seed(size = ss)
		forall(m=1:ss) bdayseed(m) = (m-1)*100 + seed0
		call random_seed(put = bdayseed(1:ss) )

		do i=1,Nsim
			do it=1,Tsim
				call rand_num_closed(s_innov)
				status_it_innov(i,it) = s_innov
			enddo
		enddo
		success = 0
	end subroutine draw_status_innov
	
	subroutine draw_alit(al_it,al_it_int, seed0, success)
	! draws alpha shocks and idices on the alpha grid (i.e at the discrete values)
		implicit none

		integer, intent(in) :: seed0
		integer, intent(out),optional :: success
		real(dp), dimension(:,:) :: al_it
		integer, dimension(:,:) :: al_it_int
		integer :: ss, Nsim, alfgrid_int, t,m,i,k
		real(dp) :: alfgridL, alfgridH,alf_innov,alfgrid_i
		integer, dimension(100) :: bdayseed
		real(dp), allocatable :: cumpi_al(:,:)

		allocate(cumpi_al(nal,nal+1))

		call random_seed(size = ss)
		forall(m=1:ss) bdayseed(m) = (m-1)*100 + seed0
		call random_seed(put = bdayseed(1:ss) )

		cumpi_al =0.
		Nsim = size(al_it,1)

		do i=1,nal
			do k=2,nal+1
				cumpi_al(i,k) = pialf(i,k-1)+cumpi_al(i,k-1)
			enddo
		enddo

		success =0
		do i=1,Nsim

			! draw starting values
			t =1
			
			call random_normal(alf_innov) ! draw normal disturbances on 0,1
			! transform it by the ergodic distribution for the first period:
			alfgrid_i = alf_innov*alfsig + alfmu

			if(alfgrid_i >maxval(alfgrid) .or. alfgrid_i < minval(alfgrid)) success = 1+success !count how often we truncate
			!impose bounds
			alfgrid_i = max(alfgrid_i,minval(alfgrid))
			alfgrid_i = min(alfgrid_i,maxval(alfgrid))
			alfgrid_int = finder(alfgrid,alfgrid_i)
			! round up or down:
			if( (alfgrid_i - alfgrid(alfgrid_int))/(alfgrid(alfgrid_int+1)- alfgrid(alfgrid_int)) >0.5 ) alfgrid_int = alfgrid_int + 1
			if(al_contin .eqv. .true.) then
				al_it(i,t) = alfgrid_i ! log of wage shock
			else
				al_it(i,t) = alfgrid(alfgrid_int) ! log of wage shock, on grid
			endif
			al_it_int(i,t) = alfgrid_int
			
			! draw sequence:

			do t=2,Tsim
				call rand_num_closed(alf_innov)
				alfgrid_int 	= finder(cumpi_al(alfgrid_int,:), alf_innov )
				alfgrid_i	= alfgrid(alfgrid_int)
				if(al_contin .eqv. .true.) then
					al_it(i,t) = alfgrid_i ! log of wage shock
				else
					al_it(i,t) = alfgrid(alfgrid_int) ! log of wage shock, on grid
				endif
				al_it_int(i,t) = alfgrid_int					
			enddo
		enddo
		if(success > 0.2*Nsim*Tsim)  success = success
		if(success <= 0.2*Nsim*Tsim) success = 0
		!call mat2csv(cumpi_al,"cumpi_al.csv")
		deallocate(cumpi_al)
		
	end subroutine draw_alit

	subroutine draw_ji(j_i,seed0, success)
		implicit none
		integer	:: j_i(:)
		real(dp) :: jwt
		integer	:: i,m,ss
		integer,intent(in)  :: seed0
		integer,intent(out) :: success
		integer, dimension(100) :: bdayseed
		real(dp)	:: Njcumdist(nj+1)
		real(dp) :: draw_i

		Njcumdist = 0
		if(nj>1) then
			do i=1,nj
				Njcumdist(i+1) = occsz0(i) + Njcumdist(i)
			enddo
		
			call random_seed(size = ss)
			forall(m=1:ss) bdayseed(m) = (m-1)*100 + seed0
			call random_seed(put = bdayseed(1:ss) )
		
			do i=1,Nsim
				call random_number(draw_i)
				j_i(i) = finder(Njcumdist,draw_i)
				if(j_i(i) < 1 ) j_i(i) = 1
				if(j_i(i) > nj) j_i(i) = nj
			enddo
		else
			j_i = 1
		endif


		success = 0
		
	end subroutine draw_ji

	subroutine draw_zjt(z_jt_macro, seed0, success)
		implicit none

		integer, intent(out) :: z_jt_macro(:) !this will be the panel across occupations -> z_jt by i's j
		integer	:: it,i,ij,iz,izp,m,ss, z_jt_t
		integer,intent(in) :: seed0
		integer,intent(out) :: success
		integer, dimension(100) :: bdayseed
		integer, dimension(5,2) :: NBER_start_stop
		real(dp) :: z_innov
		real(dp) :: cumpi_z(nz,nz+1),cumpi_zblock(nz/2,nz/2+1)

		call random_seed(size = ss)
		forall(m=1:ss) bdayseed(m) = (m-1)*100 + seed0
		call random_seed(put = bdayseed(1:ss) )

		cumpi_z = 0.
		cumpi_zblock = 0.
		! for random transitions of time block
		do iz=1,nz
			do izp=1,nz
				cumpi_z(iz,izp+1) = piz(iz,izp) + cumpi_z(iz,izp)
			enddo
		enddo
		! for deterministic transition of time block
		do iz=1,nz/2
			do izp=1,nz/2
				cumpi_zblock(iz,izp+1) = piz(iz,izp) + cumpi_zblock(iz,izp)
			enddo
			cumpi_zblock(iz,:) = cumpi_zblock(iz,:)/cumpi_zblock(iz,nz/2)
		enddo

		!compute NBER dates in case of NBER_tseq == 1 
		! 1980 + 0/4 -> 1980 + 2/4
		! 1981 + 2/4 -> 1982 + 3/4
		! 1990 + 2/4 -> 1991 + 0/4
		! 2001 + 2/4 -> 2001 + 3/4
		! 2007 + 3/4 -> 2009 + 1/4
		NBER_start_stop(1,1) =  0*itlen + 0*dnint( tlen/4. ) +1
		NBER_start_stop(1,2) =  0*itlen + 2*dnint( tlen/4. ) +1 
		NBER_start_stop(2,1) =  1*itlen + 2*dnint( tlen/4. ) +1
		NBER_start_stop(2,2) =  2*itlen + 3*dnint( tlen/4. ) +1
		NBER_start_stop(3,1) = 10*itlen + 2*dnint( tlen/4. ) +1
		NBER_start_stop(3,2) = 11*itlen + 0*dnint( tlen/4. ) +1
		NBER_start_stop(4,1) = 21*itlen + 2*dnint( tlen/4. ) +1
		NBER_start_stop(4,2) = 21*itlen + 3*dnint( tlen/4. ) +1
		NBER_start_stop(5,1) = 27*itlen + 3*dnint( tlen/4. ) +1
		NBER_start_stop(5,2) = 29*itlen + 1*dnint( tlen/4. ) +1
		!start cycles on 1st
		ss = 1

		
		!draw on zgrid
		! start everyone from middle state, alternatively could start from random draw on ergodic dist
		z_jt_t = (nz/2+1)/2
		do it = 1,Tsim
			if(NBER_tseq .eq. 1 ) then
				if(it >= NBER_start_stop(ss,1) .and. it < NBER_start_stop(ss,2) ) then
					z_innov = 0.
				elseif( it == NBER_start_stop(ss,2) ) then
					z_innov = 1.
					ss = ss+1
				else
					z_innov = 0.5
				endif
			else
				call rand_num_closed(z_innov)
				if(z_innov<0. .or. z_innov >1) then
					call rand_num_closed(z_innov)
				endif
			endif	
			! use conditional probability w/in time block
			z_jt_t = finder(cumpi_zblock(z_jt_t,:),z_innov )
			if( dble(it)<struc_brk*tlen ) then
				z_jt_macro(it) = z_jt_t
			else
				z_jt_macro(it) = z_jt_t + nz/2
			endif
			! allow random time-block transitions
			! z_jt_t = finder(cumpi_z(z_jt_t,:),z_innov )
			
		enddo
		success = 0
		!call mat2csv(cumpi_z,"cumpi_z.csv")
	end subroutine draw_zjt
	
	subroutine draw_age_it(age_it, born_it, seed0, success)

		integer,intent(out) :: age_it(:,:),born_it(:,:)
		integer	:: it,itp,i,m,ss,bn_i
		integer,intent(in) :: seed0
		integer,intent(out) :: success
		integer, dimension(100) :: bdayseed
		real(dp), dimension(TT) :: cumpi_t0
		real(dp), dimension(TT-1) :: prob_age_nTT
		real(dp) :: rand_age,rand_born
		
		call random_seed(size = ss)
		forall(m=1:ss) bdayseed(m) = (m-1)*100 + seed0
		call random_seed(put = bdayseed(1:ss) )

		!set up cumulative probabilities for t0 and conditional draws
		!not begining with anyone from TT
		forall(it=1:TT-1) prob_age_nTT(it) = prob_age(it)/(1-prob_age(TT))
		!also include TT
		!prob_age_nTT = prob_age
		
		cumpi_t0 = 0.

		do it=1,TT-1
			cumpi_t0(it+1) = prob_age_nTT(it) + cumpi_t0(it)
		enddo
		
		born_it = 0
		
		do i =1,Nsim
		do m = 1,50
			bn_i = 0 ! initialize, not yet born
			it = 1
			call rand_num_closed(rand_born)
			if(rand_born < hazborn_t(1)) then 
				born_it(i,it) = 0 ! no one is "born" in the first period
				bn_i = 1
				!draw an age
				call rand_num_closed(rand_age)
				age_it(i,it) = finder(cumpi_t0,rand_age)
			else 
				age_it(i,it) = 0
				born_it(i,it) = 0
			endif
			do it=2,Tsim
				if(age_it(i,it-1)<TT) then
					call rand_num_closed(rand_born)
					if(rand_born < hazborn_t(it) .and. bn_i == 0 ) then
						age_it(i,it) =1
						born_it(i,it) = 1
						bn_i = 1
					elseif(bn_i == 1) then
						born_it(i,it) = 0
						call rand_num_closed(rand_age)
						if(rand_age < 1- ptau(age_it(i,it-1)) .and. age_it(i,it-1) < TT ) then
							age_it(i,it) = age_it(i,it-1)+1
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
			if(bn_i == 1) then !if the born draw never comes, do the whole thing again.
				exit
			endif
		enddo! m=1,50.... if never gets born

		enddo! i=1,Nsim

		success = 0
	end subroutine draw_age_it

	subroutine draw_draw(drawi_ititer, drawt_ititer, age_it, seed0, success)

		integer,intent(in) :: seed0
		integer,intent(out) :: success
		integer,intent(in) :: age_it(:,:)
		integer, dimension(100) :: bdayseed
		integer,intent(out) :: drawi_ititer(:),drawt_ititer(:)
		integer :: i,it,ii,m,ss,drawt,drawi,ndraw
		real(dp) :: junk

		call random_seed(size = ss)
		forall(m=1:ss) bdayseed(m) = (m-1)*100 + seed0
		call random_seed(put = bdayseed(1:ss) )

		Ndraw = size(drawi_ititer)
		!need to draw these from age-specific distributions for iterations > 1
		do i=1,Ndraw
			it=1
			if(age_it(i,it) > 0 )then
				ageloop: do 
					call random_number(junk)
					drawi = max(1,idnint(junk*Nsim))
					call random_number(junk)
					drawt = max(1,idnint(junk*(Tblock_sim*tlen)-1))
					if((age_it(drawi,drawt) .eq. age_it(i,it))) then
						exit ageloop
					endif
				end do ageloop
				drawi_ititer(i) = drawi
				drawt_ititer(i) = drawt
			else
				drawi_ititer(i) = i
				drawt_ititer(i) = it
			endif
		enddo
		success = 0
	end subroutine


	subroutine draw_shocks(hst)

		implicit none
		type(hist_struct) :: hst
		integer :: seed0,seed1, status 

		seed0 = 941987
		seed1 = 12281951

		call draw_age_it(hst%age_hist,hst%born_hist,seed0,status)
		seed0 = seed0 + 1
		call draw_deli(hst%del_i, hst%del_i_int, seed1, status)
		seed1 = seed1 + 1
		call draw_alit(hst%al_hist,hst%al_int_hist, seed0, status)
		seed0 = seed0 + 1
		call draw_ji(hst%j_i,seed1, status)
		seed1 = seed1 + 1
		call draw_zjt(hst%z_jt_macro, seed0, status)
		seed0 = seed0 + 1

		hst%drawn = 1

	end subroutine draw_shocks

	
	subroutine sim(vfs, pfs,hst)
		
		implicit none

		type(val_struct), intent(inout), target :: vfs
		type(pol_struct), intent(inout), target :: pfs	
		type(hist_struct), intent(inout), target :: hst
		


		integer :: i, ii, iter, it, it_old, ij, idi, id, Tret, &
			&  seed0, seed1, status, m,ss, iter_draws=5
		integer :: bdayseed(100)

		real(dp), allocatable ::	del_i(:), jshock_ij(:,:) ! shocks to be drawn
		
		real(dp), allocatable :: status_it_innov(:,:) !innovations to d, drawn randomly

		integer, allocatable :: work_it(:,:), app_it(:,:) !choose work or not, apply or not
		real(dp), allocatable :: e_it(:,:)
		integer, allocatable :: a_it_int(:,:),e_it_int(:,:)
		integer, allocatable :: drawi_ititer(:),drawt_ititer(:)
		
		! write to hst
		real(dp), pointer     :: work_dif_it(:,:), app_dif_it(:,:) !choose work or not, apply or not -- latent value
		integer, pointer     :: born_it(:,:) ! born status, drawn randomly		
		integer, pointer     :: del_i_int(:)  ! integer valued shocks
		integer, pointer     :: status_it(:,:)	!track W,U,N,D,R : 1,2,3,4,5
		integer, pointer     :: age_it(:,:)	! ages, drawn randomly
		integer, pointer     :: j_i(:)		! occupation, maybe random
		integer, pointer     :: d_it(:,:) 	! disability status
		integer, pointer     :: z_jt_macro(:)! shocks to be drawn
		integer, pointer     :: al_it_int(:,:)! integer valued shocks
		real(dp), pointer     :: occgrow_jt(:,:), occshrink_jt(:,:), occsize_jt(:,:)
		real(dp), pointer     :: a_it(:,:) 	! assets
		real(dp), pointer     ::	al_it(:,:)	! individual shocks

		! read from vals
		real(dp), pointer ::	V(:,:,:,:,:,:,:)	!Participant

		! read from pols
		real(dp), pointer ::	gapp_dif(:,:,:,:,:,:,:), gwork_dif(:,:,:,:,:,:,:) ! latent value of work/apply
	
		integer, pointer ::	aR(:,:,:), aD(:,:,:,:), aU(:,:,:,:,:,:,:), &
					aN(:,:,:,:,:,:,:), aW(:,:,:,:,:,:,:)
		integer, pointer ::	gapp(:,:,:,:,:,:,:), &
					gwork(:,:,:,:,:,:,:)

		logical :: ptrsuccess
		real(dp) :: cumpid(nd,nd+1,ndi,TT-1),cumptau(TT+1),a_mean(TT-1),d_mean(TT-1),a_var(TT-1),d_var(TT-1),&
				& a_mean_liter(TT-1),d_mean_liter(TT-1),a_var_liter(TT-1),d_var_liter(TT-1), &
				& s_mean(TT-1),s_mean_liter(TT-1)
	
		! Other
		real(dp)	:: wage_hr,al_hr, junk,a_hr, e_hr, bet_hr,z_hr,j_val,j_val_ij,jwt, vscale,cumval,work_dif_hr, app_dif_hr,js_ij, Nworkt, ep_hr

		integer :: ali_hr,d_hr,age_hr,del_hr, zi_hr, j_hr, ai_hr,api_hr,ei_hr, &
			& beti, status_hr,status_tmrw,drawi,drawt
		!************************************************************************************************!
		! Allocate things
		!************************************************************************************************!

		iter_draws = maxiter !globally set variable
		
		allocate(jshock_ij(Nsim,nj))
		allocate(del_i(Nsim))
		allocate(status_it_innov(Nsim,Tsim))
		
		allocate(a_it_int(Nsim,Tsim))		
		allocate(e_it(Nsim,Tsim))
		allocate(e_it_int(Nsim,Tsim))		
		allocate(work_it(Nsim,Tsim))
		allocate(app_it(Nsim,Tsim))

		!this must be big enough that we are sure it's big enough that can always find a worker
		allocate(drawi_ititer(Nsim*10)) 
		allocate(drawt_ititer(Nsim*10))


		!************************************************************************************************!
		! Pointers
		!************************************************************************************************!
		! (disability extent, earn hist, assets)

		V => vfs%V !need this for the career choice
		aR => pfs%aR
		aD => pfs%aD
		aN => pfs%aN
		aW => pfs%aW
		aU => pfs%aU
		gwork => pfs%gwork
		gapp => pfs%gapp

		gapp_dif    => pfs%gapp_dif
		gwork_dif   => pfs%gwork_dif

		status_it   => hst%status_hist
		age_it 	    => hst%age_hist
		born_it	    => hst%born_hist
		work_dif_it => hst%work_dif_hist
		app_dif_it  => hst%app_dif_hist
		d_it        => hst%d_hist
		del_i_int   => hst%del_i_int
		j_i         => hst%j_i
		z_jt_macro  => hst%z_jt_macro
		a_it        => hst%a_hist
		al_it       => hst%al_hist
		al_it_int	=> hst%al_int_hist
		
		occgrow_jt  => hst%occgrow_jt
		occshrink_jt=> hst%occshrink_jt
		occsize_jt  => hst%occsize_jt

		ptrsuccess = associated(d_it,hst%d_hist)
		if(verbose>1 .and. ptrsuccess .eqv. .false. ) print *, "failed to associate d_it"
		ptrsuccess = associated(age_it,hst%age_hist)
		if(verbose>1 .and. ptrsuccess .eqv. .false. ) print *, "failed to associate age_it"
		ptrsuccess = associated(V,vfs%V)
		if(verbose>1 .and. ptrsuccess .eqv. .false. ) print *, "failed to associate V"
		
		Tret = (Longev - youngD - oldD*oldN)*tlen
		work_dif_it = 0.
		app_dif_it = 0.


		seed0 = 941987
		seed1 = 12281951
		call random_seed(size = ss)
		forall(m=1:ss) bdayseed(m) = (m-1)*100 + seed0
		call random_seed(put = bdayseed(1:ss) )

		if(verbose >2) print *, "Drawing types and shocks"	
		if(j_rand .eqv. .false. ) then
			!draw gumbel-distributed shock
			do i = 1,Nsim
				do ij = 1,nj
					call random_gumbel(js_ij)
					jshock_ij(i,ij)=js_ij
				enddo
			enddo
		endif
		call draw_shocks(hst)
		
		if(hst%drawn /= 1 )then
			call draw_age_it(age_it,born_it,seed0,status)
			seed0 = seed0 + 1
			call draw_deli(del_i, del_i_int, seed1, status)
			seed1 = seed1 + 1
			call draw_alit(al_it,al_it_int, seed0, status)
			seed0 = seed0 + 1
			call draw_ji(j_i,seed1, status)
			seed1 = seed1 + 1
			call draw_zjt(z_jt_macro, seed0, status)
			seed0 = seed0 + 1
		endif
		
		call draw_status_innov(status_it_innov,seed1,status)
		seed1 = seed1 + 1
		call draw_draw(drawi_ititer, drawt_ititer, age_it, seed0, status)
		seed0 = seed0 + 1
		
		! check the distributions
		if(print_lev > 1 ) then
			call vec2csv(del_i,"del_i.csv")
			call veci2csv(j_i,"j_i.csv")
			call mat2csv(al_it,"al_it.csv")
			call veci2csv(z_jt_macro,"z_jt_macro.csv")
			call mati2csv(al_it_int,"al_it_int.csv")
			call mati2csv(age_it,"age_it.csv")
			call veci2csv(drawi_ititer,"drawi.csv")
		endif
		
		!set up cumpid,cumptau
		cumpid = 0.
		cumptau = 0.
		do idi=1,ndi
		do it =1,TT-1
		do id =1,nd
			i=1
			cumpid(id,i+1,idi,it) = pid(id,i,idi,it)
			do i =2,nd
				cumpid(id,i+1,idi,it) = pid(id,i,idi,it)+cumpid(id,i,idi,it)
			enddo
		enddo
		enddo
		enddo
		it = 1
		cumptau(it+1)=cumptau(it)
		do it =2,TT
			cumptau(it+1) = cumptau(it)+ptau(it)
		enddo

		! will draw these from endogenous distributions the second time around
		d_it = 1
		a_it = agrid(1)
		a_it_int = 1
		e_it = egrid(1)
		e_it_int = 1
		status_it = 1 ! just to initialize on the first round - everyone working age starts working
		where(age_it>=TT) status_it = 5
		where(age_it <=0) status_it = 0
		
		a_mean_liter = 0.
		d_mean_liter = 0.
		s_mean_liter = 0.
		a_var_liter  = 0.
		d_var_liter  = 0.


		!use only 1 value of beta
		beti = 1
		bet_hr = 1._dp
		
		if(verbose >2) print *, "Simulating"
		
		!itertate to get dist of asset/earnings correct at each age from which to draw start conditions 
		do iter=1,iter_draws
		if(verbose >3) print *, "iter: ", iter
			it = 1
			ii = 1
			vscale = 0.
			junk =0.
			do i =1,Nsim
				!for the population that is pre-existing in the first period 
				!need to draw these from age-specific distributions for iterations > 1
				if(iter>1 .and. age_it(i,it)  > 0 ) then
					drawi = drawi_ititer(i)!iter-1
					drawt = drawt_ititer(i)!iter-1
					status_it(i,it) = status_it(drawi,drawt)
					d_it(i,it) = d_it(drawi,drawt)
					a_it(i,it) = a_it(drawi,drawt)
					e_it(i,it) = e_it(drawi,drawt)
					e_it_int(i,it) = e_it_int(drawi,drawt)
					a_it_int(i,it) = a_it_int(drawi,drawt)
					ii =ii +1

					if(ii> size(drawi_ititer) ) then
						if(verbose >2) print *, "Make drawi_ititer larger, ran out of workers"
						ii = 1 !should never happen
					endif
				endif
				! get straight the scale of jval shocks (for use with amenityscale)
				if(born_it(i,it).eq. 1 .and. it> 1) then ! no one is ``born'' in the first period, but just to be sure
					junk    = junk+1
					age_hr	= 1
					d_hr	= 1
					a_hr 	= minval(agrid)
					ai_hr	= 1
					ei_hr	= 1
					e_hr 	= 0.
					ai_hr 	= 1
					if(j_rand .eqv. .false. ) then !choose occupation
						do ij = 1,nj
							j_val_ij = 0.
							do idi=1,ndi
								j_val_ij = V((ij-1)*nbi+beti,(idi-1)*nal+ali_hr,d_hr,ei_hr,ai_hr,z_jt_macro(it),age_hr)*delwt(idi,ij) &
												&+ j_val_ij
							enddo
							vscale = j_val_ij + vscale
						enddo
					endif
				endif 
			enddo !i=1:Nsim
			vscale = vscale/junk/dble(nj)
			
			!$OMP  parallel do &
			!$OMP& private(i,del_hr,j_hr,status_hr,it,it_old,age_hr,al_hr,ali_hr,d_hr,e_hr,a_hr,ei_hr,ai_hr,z_hr,zi_hr,api_hr,ep_hr, &
			!$OMP& ij,j_val,j_val_ij,cumval,jwt,wage_hr,junk,app_dif_hr,work_dif_hr,status_tmrw) 
			do i=1,Nsim
				!fixed traits
	
				!set a j to correspond to the probabilities.  This will get overwritten if born
				j_hr = j_i(i)

				del_hr = del_i_int(i) 

				!initialize stuff
				it = 1
				it_old = 1
				
				do it=1,Tsim
				if(age_it(i,it) > 0) then !they've been born 
					!set the state
					al_hr	= al_it(i,it)
					ali_hr	= al_it_int(i,it)
					if(born_it(i,it).eq. 1 .and. it> 1) then ! no one is ``born'' in the first period, but just to be sure
						age_hr	= 1
						d_hr	= 1
						a_hr 	= minval(agrid)
						ai_hr	= 1
						ei_hr	= 1
						e_hr 	= 0.
						ai_hr 	= 1
						if(j_rand .eqv. .false. ) then !choose occupation
							j_val  = -1.e6_dp
							j_hr = 1 !initial value
							do ij = 1,nj
								j_val_ij = 0.
								do idi = 1,ndi
									j_val_ij = V((ij-1)*nbi+beti,(idi-1)*nal+ali_hr,d_hr,ei_hr,ai_hr,z_jt_macro(it),age_hr)*delwt(idi,ij)+j_val_ij
								enddo
								j_val_ij = j_val_ij + jshock_ij(i,ij)*(amenityscale*vscale) + jshift(ij)*(amenityscale*vscale)
								if(j_val< j_val_ij ) then
									j_hr = ij
									j_val = j_val_ij
								endif
							enddo
							j_i(i) = j_hr
							if(del_by_occ .eqv. .true.) then
								do ii=1,ndi
									if( del_i(i)< delcumwt(ii+1,j_hr)  ) then
										del_hr = ii
										exit
									endif
								enddo
							endif

							del_i_int(i) = del_hr
						endif
					else 
						age_hr	= age_it(i,it)
						d_hr	= d_it(i,it)
						a_hr 	= a_it(i,it)
						ei_hr	= e_it_int(i,it)
						e_hr 	= e_it(i,it)
						ai_hr 	= a_it_int(i,it)
					endif
					zi_hr	= z_jt_macro(it)
					z_hr	= zgrid(zi_hr,j_hr)
					
					wage_hr	= wage(bet_hr,al_hr,d_hr,z_hr,age_hr)
					hst%wage_hist(i,it) = wage_hr
					
					status_hr = status_it(i,it)
					! get set to kill off old (i.e. age_hr ==TT only for Longev - youngD - oldD*oldN )
					if((age_hr .eq. TT) ) then
						it_old = it_old + 1
						if(it_old >  Tret ) then
							age_it(i,it:Tsim) = -1
							a_it(i,it:Tsim) = 0.
							a_it_int(i,it:Tsim) = 0
							d_it(i,it:Tsim) = 0
							app_dif_it(i,it:Tsim) = 0.						
							work_dif_it(i,it:Tsim) = 0.
							status_it(i,it:Tsim) = -1
							exit
						else
							status_it(i,it:Tsim) = 5
						endif
						
					endif 

					!make decisions if not yet retired
					if(age_hr < TT) then 
						if(status_hr .eq. 3) then ! choose wait or apply
							app_dif_hr = gapp_dif( (j_hr-1)*nbi + beti, (del_hr-1)*nal+ali_hr,d_hr,ei_hr,ai_hr,zi_hr,age_hr )
							app_dif_it(i,it) = app_dif_hr
							if( app_dif_hr >= 0 ) then
							! choose to apply
								app_it(i,it) = 1
							else
								app_it(i,it) = 0
							endif
						endif
						! evalutate gwork and gapp to figure out lom of status 
						if((status_hr < 3) .or. (status_hr .eq. 3 .and. app_dif_hr < 0 ))then !choose work or rest
							app_dif_it(i,it) = 0. !just to fill in value
							work_dif_hr = gwork_dif( (j_hr-1)*nbi + beti, (del_hr-1)*nal+ali_hr,d_hr,ei_hr,ai_hr,zi_hr,age_hr )
							
							if( work_dif_hr > 0 ) then
							! choose to work
								if(status_hr < 3) then !not LTU
									status_it(i,it) = 1
									status_tmrw = 1
								else !LTU, have to get a good shock 
									if(status_it_innov(i,it) <=lrho) status_tmrw = 1
									if(status_it_innov(i,it) > lrho)  status_tmrw = status_hr
								endif
							elseif(status_hr .le. 2) then
								status_it(i,it) = 2
							! unemployed, may stay unemployed or become long-term unemployed
								if(status_it_innov(i,it) <=pphi) status_tmrw = 3
								if(status_it_innov(i,it) > pphi) status_tmrw = 2
							else
							!NDR no change, though may choose to apply for DI below
								status_tmrw =  status_hr
							endif
							work_dif_it(i,it) = work_dif_hr
						elseif(status_hr .eq. 3 .and. app_dif_hr >=0 ) then
							! eligible to apply?
							if(age_hr > 1 .or. (age_hr ==1 .and. status_it_innov(i,min(it+1,Tsim))<eligY )) then !status_it_innov(i,it+1) is an independent draw
								!applying, do you get it?
								if(status_it_innov(i,it) < xifun(d_hr,z_hr,age_hr)) then 
									status_tmrw = 4
								else	
									status_tmrw = 3
								endif
							else
								status_tmrw = 3
							endif
						elseif(status_hr > 3 ) then !absorbing states of D,R
							status_tmrw = status_hr
							! just to fill in values
							app_dif_it(i,it) = 0.
							work_dif_it(i,it) = 0.

						endif
						!evaluate the asset policy			
						if(status_hr .eq. 1) then
							api_hr = aw( (j_hr-1)*nbi + beti, (del_hr-1)*nal+ali_hr,d_hr,ei_hr,ai_hr,zi_hr,age_hr  )
						elseif(status_hr .eq. 2) then
							api_hr = aU( (j_hr-1)*nbi + beti, (del_hr-1)*nal+ali_hr,d_hr,ei_hr,ai_hr,zi_hr,age_hr )
						elseif(status_hr .eq. 3) then
							api_hr = aN( (j_hr-1)*nbi + beti, (del_hr-1)*nal+ali_hr,d_hr,ei_hr,ai_hr,zi_hr,age_hr )
						elseif(status_hr .eq. 4) then
							api_hr = aD( d_hr,ei_hr,ai_hr,age_hr )
						elseif(status_hr .eq. 5) then ! should never be in this condition
							api_hr = aR(d_hr,ei_hr,ai_hr)
						endif
					! retired
					elseif( (age_hr==TT) .and. (it_old <= Tret)) then
						api_hr = aR( d_hr,ei_hr,ai_hr )
						status_hr = 5
						status_tmrw = 5
						if(it<Tsim) &
						&	status_it(i,it+1) = status_hr
					endif

					!push forward the state:
					if(it<Tsim) then
						! push forward status
						status_it(i,it+1) = status_tmrw
						! push forward asset					
						a_it_int(i,it+1) = api_hr
						a_it(i,it+1) = agrid(api_hr)
						!push forward AIME
						if(status_hr .eq. 1) then
							!here, it is continuous
							ep_hr = min( (e_hr*dble(it-1) + wage_hr)/dble(it),egrid(ne) )
							e_it(i,it+1) = ep_hr
							! assign to grid points by nearest neighbor
							! ei_hr = finder(egrid,e_it(i,it+1)) <- relatively short, but not thread safe
							do ei_hr = ne,1,-1
								if( ep_hr > egrid(ei_hr) ) then
									exit
								endif
							enddo
							ei_hr = max(ei_hr,1) !just to be sure we take base 1
							if(ei_hr < ne) then
								if( (ep_hr - egrid(ei_hr))<  (egrid(ei_hr+1) - ep_hr) ) then
									e_it_int(i,it+1) = ei_hr
								else
									e_it_int(i,it+1) = ei_hr + 1
								endif
							else
								e_it_int(i,it+1) = ne
							endif
						else
							e_it(i,it+1) = e_hr
							e_it_int(i,it+1) = ei_hr
						endif

						!push forward d 
						if(status_hr .eq. 1 .and. d_hr<nd ) then !if working and not already in worst state
							if( status_it_innov(i,it) < pid(d_hr,d_hr+1,del_hr,age_hr) ) then 
								d_it(i,it+1) = d_hr+1 
							else
								d_it(i,it+1) = d_hr
							endif
						else 
							d_it(i,it+1) = d_hr
						endif
					endif
				else !age_it(i,it) <= 0, they've not been born
					a_it(i,it) = 0.
					a_it_int(i,it) = 0
					d_it(i,it) = 0
					app_dif_it(i,it) = 0.						
					work_dif_it(i,it) = 0.
					status_it(i,it) = -1
				endif ! age_it(i,it)>0
				enddo !1,Tsim
			enddo! 1,Nsim
			!$OMP  end parallel do 

			
			if(print_lev >=3)then
				call mat2csv (e_it,"e_it.csv")
				call mat2csv (a_it,"a_it.csv")
				call mati2csv(a_it_int,"a_it_int.csv")
				call mati2csv(status_it,"status_it.csv")
				call mati2csv(d_it,"d_it.csv")
			endif

			a_mean = 0.
			d_mean = 0.
			s_mean = 0.
			a_var = 0.
			d_var = 0.
			do age_hr = 1,TT-1
				junk = 0.
				do i=1,Nsim
					do it = 1,Tsim
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
					do it = 1,Tsim
						if( age_hr .eq. age_it(i,it) ) then
							a_var(age_hr) = (dlog(a_it(i,it)) - dlog(a_mean(age_hr)))**2 + a_var(age_hr)
							d_var(age_hr) = (d_it(i,it) - d_mean(age_hr))**2+ d_var(age_hr)
						endif
					enddo
				enddo
				a_var(age_hr) = a_var(age_hr)/junk
				d_var(age_hr) = d_var(age_hr)/junk

			enddo
			if((sum((a_mean - a_mean_liter)**2) + sum((s_mean - s_mean_liter)**2) + sum((d_mean - d_mean_liter)**2))<1.e-5_dp ) then
				if(verbose >=2 ) then
					print *, "done simulating after convergence in", iter
					print *, "dif a mean, log a var",  sum((a_mean - a_mean_liter)**2), sum((a_var - a_var_liter)**2)
					! NOTE: this is not actually mean because it does not have demographic weights
					print *, "a mean, a var",  sum(a_mean), sum(a_var)
					print *, "dif d mean,     d var",  sum((d_mean - d_mean_liter)**2), sum((d_var - d_var_liter)**2)
					print *, "dif status mean",  sum((s_mean - s_mean_liter)**2)
					print *, "status mean", sum(s_mean)
					print *,  "-------------------------------"
				endif
				exit
			else
				if(verbose >=4 .or. (verbose >=2 .and. mod(iter,100) == 0)) then
					print *, "iter:", iter
					print *, "dif a mean, log a var",  sum((a_mean - a_mean_liter)**2), sum((a_var - a_var_liter)**2)
					print *, "a mean, a var",  sum(a_mean), sum(a_var)
					print *, "dif d mean,     d var",  sum((d_mean - d_mean_liter)**2), sum((d_var - d_var_liter)**2)
					print *, "dif status mean",  sum((s_mean - s_mean_liter)**2)
					print *,  "-------------------------------"
				endif
			endif
			a_mean_liter = a_mean
			d_mean_liter = d_mean
			s_mean_liter = s_mean
			a_var_liter = a_var
			d_var_liter = d_var
		enddo! iter
		
		! calc occupation growth rates
		if(verbose >2) print *, "calculating occupation growth rates"		
		it = 1
		Nworkt = 0. !labor force in first period are the ones "born" in the first period
		occsize_jt = 0.
		occgrow_jt = 0.
		occshrink_jt = 0.
		do i=1,Nsim
			if( age_it(i,it) > 0 .and. status_it(i,it) <= 2 .and. status_it(i,it)>=1) Nworkt = 1._dp + Nworkt
			do ij=1,nj
				if(j_i(i) == ij .and. age_it(i,it) >= 1 .and. status_it(i,it) == 1) &
					& occsize_jt(ij,it) = 1._dp+occsize_jt(ij,it)
			enddo
		enddo
		occsize_jt(:,it) = occsize_jt(:,it) / Nworkt

		!$omp parallel do private(it,Nworkt,i,ij)
		do it = 2,Tsim
			Nworkt = 0.
			occgrow_jt  (:,it) = 0.
			occshrink_jt(:,it) = 0.
			occsize_jt  (:,it) = 0.
			do i=1,Nsim
				if(status_it(i,it) <= 2 .and. status_it(i,it)>=1 .and. age_it(i,it) > 0 ) Nworkt = 1._dp + Nworkt !labor force in this period
				do ij =1,nj
					if(j_i(i) ==ij .and. born_it(i,it) == 1 ) &
						& occgrow_jt(ij,it) = 1._dp + occgrow_jt(ij,it)
					if(j_i(i) ==ij .and. status_it(i,it-1) > 1  .and. status_it(i,it)== 1) &
						& occgrow_jt(ij,it) = 1._dp + occgrow_jt(ij,it) !wasn't working last period
					if(j_i(i) ==ij .and. status_it(i,it-1) == 1 .and. status_it(i,it) > 1) &
						& occshrink_jt(ij,it) = 1._dp + occshrink_jt(ij,it)
					if(j_i(i) ==ij .and. status_it(i,it) == 1) &
						& occsize_jt(ij,it) = 1._dp + occsize_jt(ij,it)
				enddo
				occgrow_jt(ij,it) = occgrow_jt(ij,it)/occsize_jt(ij,it)
				occshrink_jt(ij,it) = occshrink_jt(ij,it)/occsize_jt(ij,it)
			enddo
			forall(ij=1:nj) occsize_jt(ij,it) = occsize_jt(ij,it)/Nworkt
		enddo
		!$omp end parallel do

		
		if(print_lev > 1)then
				call mat2csv (e_it,"e_it.csv")
				call mat2csv (a_it,"a_it.csv")
				call mati2csv(a_it_int,"a_it_int.csv")
				call mati2csv(status_it,"status_it.csv")
				call mati2csv(d_it,"d_it.csv")
				call veci2csv(j_i,"j_i.csv")
				call veci2csv(z_jt_macro,"z_jt.csv")
				call mat2csv (occsize_jt,"occsize_jt.csv")
				call mat2csv (occgrow_jt,"occgrow_jt.csv")
				call mat2csv (occshrink_jt,"occshrink_jt.csv")
				call mat2csv (hst%wage_hist,"wage_hist.csv")
				call mat2csv (hst%app_dif_hist,"app_dif_hist.csv")
				call mat2csv (hst%work_dif_hist,"work_dif_hist.csv")
		endif
	
		deallocate(e_it)
		deallocate(a_it_int,e_it_int)
		deallocate(del_i)
		deallocate(app_it,work_it)
		deallocate(status_it_innov)
		deallocate(drawi_ititer,drawt_ititer)

	end subroutine sim

end module sim_hists



!**************************************************************************************************************!
!**************************************************************************************************************!
!						Searches for parameters						       !
!**************************************************************************************************************!
!**************************************************************************************************************!
module find_params

	use V0para
	use helper_funs
	use sol_val
	use sim_hists
	use model_data

	implicit none

	real(dp), allocatable :: mod_prob_hist_tgt(:,:)

	integer :: mod_ij_obj, mod_it
	type(val_struct), pointer :: mod_val_sol
	type(hist_struct), pointer :: mod_hists_sim
	real(dp) :: mod_prob_target

	contains

	function obj_zj(zj_in,val_sol,hists_sim,prob_target,ij_obj,it0) result(resid)
		
		! for each j in each t, this solves the implied z_j
		! where z_j scales the rest of the markov chain
		
		integer, intent(in) :: ij_obj, it0
		type(val_struct),  intent(in) :: val_sol
		type(hist_struct), intent(in):: hists_sim
		real(dp), intent(in) :: prob_target,zj_in
		real(dp) :: resid, zj_here,nborn,j_val_hi,j_val_lo,j_scale
		real(dp) :: cumval,j_val,zjwt,prob_here,e_hr
		integer :: age_hr,d_hr,ai_hr,ali_hr,ei_hr,i,ii,idi, it, ij,beti
		integer :: zj_hi,zj_lo, simT
		real(dp), allocatable :: j_val_ij(:,:),j_val_here(:)
		real(dp) :: zgrid_min, zgrid_max

		allocate(j_val_ij(Nsim,nj))
		allocate(j_val_here(Nsim))

		j_val_ij = -1.e10_dp
		j_val_here = -1.e10_dp

		
		simT = size(hists_sim%born_hist,2)
		beti = 1
		prob_here = 0._dp
		nborn = 0._dp
		do it = it0,simT
			do i = 1,Nsim

				if(hists_sim%born_hist(i,it).eq. 1 ) then !they've been born this period
					nborn = 1._dp + nborn
					!set the state
					ali_hr 	= hists_sim%al_int_hist(i,it)
					age_hr 	= 1
					d_hr 	= 1
					ai_hr	= 1
					ei_hr	= 1
					e_hr 	= 0._dp
					ai_hr 	= 1

					ij = ij_obj
					zgrid_min = minval(zgrid(nz/2+1:nz,ij))
					zgrid_max = maxval(zgrid(nz/2+1:nz,ij))
					! choose by scaling zgrid using the pre-transition values, 1:nz/2.  The `mod' ensures we're in that range
					zj_here   = zgrid(mod(hists_sim%z_jt_macro(it)-1,nz/2)+1,ij) + zj_in
					! zj_here	  = max(min(zj_here,zgrid_max),zgrid_min) This would prevent extrapolation
					zj_lo	  = finder(zgrid(nz/2+1:nz,ij),zj_here) + nz/2
					zj_lo	  = min(max(zj_lo,nz/2+1) , nz-1)
					zj_hi	  = min(zj_lo +1,nz)
					zjwt 	  = (zgrid(zj_hi,ij) - zj_here)/(zgrid(zj_hi,ij) - zgrid(zj_lo,ij))
					! zjwt 	  = max(min(zjwt,1._dp),0._dp)  This would prevent extrapolation
					j_val_here(i) = 0._dp
					do idi = 1,ndi ! expectation over delta
						j_val_lo  = val_sol%V((ij-1)*nbi+beti,(idi-1)*nal+ali_hr,d_hr,ei_hr,ai_hr,zj_lo,age_hr)
						j_val_hi  = val_sol%V((ij-1)*nbi+beti,(idi-1)*nal+ali_hr,d_hr,ei_hr,ai_hr,zj_hi,age_hr)
						j_val_here(i) = (zjwt*j_val_lo+ (1._dp-zjwt)*j_val_hi)*delwt(idi,ij) + j_val_here(i)
					enddo

					do ij = 1,nj
						if(ij /= ij_obj) then
							j_val_ij(i,ij) = 0.
							do idi = 1,ndi ! expectation over delta
								j_val_ij(i,ij) = val_sol%V((ij-1)*nbi+beti,(idi-1)*nal+ali_hr,d_hr,ei_hr,ai_hr, &
											& hists_sim%z_jt_macro(it),age_hr)*delwt(idi,ij) + j_val_ij(i,ij)
							enddo
						else
								j_val_ij(i,ij) = j_val_here(i)
						endif
					enddo
				endif
			enddo
		enddo

		j_scale = 0.		
		do i = 1,Nsim
		do ij = 1,nj
			if(j_val_ij(i,ij) > -1.e9_dp) then !these guys were actually born
				j_scale = j_val_ij(i,ij) + j_scale
			endif
		enddo
		enddo
		j_scale = j_scale/nborn/dble(nj)

		prob_here = 0.
		do i = 1,Nsim
			if(j_val_here(i) > -1.e9_dp) then !these guys were actually born
				cumval = 0.
				do ij=1,nj
					cumval = dexp((j_val_ij(i,ij)+jshift(ij))/(j_scale*amenityscale))+cumval
				enddo
				prob_here = dexp((j_val_here(i)+jshift(ij_obj))/(j_scale*amenityscale))/cumval + prob_here
			endif
		enddo
		
		prob_here = prob_here/nborn
		resid = prob_target - prob_here


		deallocate(j_val_ij,j_val_here)
	end function obj_zj

	function obj_zj_wrap(zj_here, dummy_params) result(resid)
		real(dp), intent(in) :: zj_here
		real(dp), intent(in) :: dummy_params(:)
		real(dp) :: resid

		resid = obj_zj(zj_here,mod_val_sol,mod_hists_sim,mod_prob_target,mod_ij_obj,mod_it)

	end function obj_zj_wrap


	subroutine newtfpgrid(zgrid_in, zscale_in, zgrid_out,val_sol)
		real(dp), intent(in)	:: zgrid_in(:,:)
		real(dp), intent(out):: zgrid_out(:,:)
		real(dp), intent(in)	:: zscale_in(:)
		type(val_struct)	:: val_sol

		integer :: iz,i,ij,beti,idi,ial,id,ie,ia,it
		real(dp) :: z_min,z_max,zjwt, zj_here
		real(dp) :: j_val_lo,j_val_hi
		integer :: zj_lo, zj_hi

		zgrid_out = zgrid_in
		! apply scaling
		do ij=1,nj
			z_max = maxval(zgrid_in(nz/2+1:nz,ij))
			z_min = minval(zgrid_in(nz/2+1:nz,ij))
			do iz=nz/2+1,nz
				zgrid_out(iz,ij) = zgrid_in(iz - nz/2,ij) + zscale(ij) !scale off the pre-break value
				! this will allow extrapolation
				!zgrid_out(iz,ij) = max(min(zgrid_out(iz,ij),z_max),z_min) ! don't try to extrapolate
			enddo
		enddo

		

		! interpolate V
		beti = 1
		do ij=1,nj
		do idi=1,ndi
		do ial = 1,nal
		do id = 1,nd
		do ie = 1,ne
		do ia=1,na
		do it = 1,TT-1

		!interpolate the value function on z grid
		!and will extrapolate for the top and bottom points (problematic?)
		do iz=nz/2+1,nz
			zj_here   = zgrid_out(iz,ij)
			zj_lo	  = finder(zgrid_in(nz/2+1:nz,ij),zj_here) + nz/2
			zj_lo	  = min( max(zj_lo,nz/2+1), nz-1)
			zj_hi	  = min(zj_lo +1,nz)
			zjwt 	  = (zgrid_in(zj_hi,ij) - zj_here)/(zgrid_in(zj_hi,ij) - zgrid_in(zj_lo,ij))
			! zjwt 	  = max(min(zjwt,1._dp),0.) ! this prevents extrapolation
			j_val_lo  = val_sol%V((ij-1)*nbi+beti,(idi-1)*nal+ial,id,ie,ia,zj_lo,it)
			j_val_hi  = val_sol%V((ij-1)*nbi+beti,(idi-1)*nal+ial,id,ie,ia,zj_hi,it)
			val_sol%V((ij-1)*nbi+beti,(idi-1)*nal+ial,id,ie,ia,iz,it) =  &
				& zjwt*j_val_lo+ (1._dp-zjwt)*j_val_hi
			
		enddo

		enddo
		enddo
		enddo
		enddo
		enddo
		enddo
		enddo

	end subroutine newtfpgrid

	subroutine jshift_sol(val_sol, hists_sim, probj_in, jshift_out)
		! solves for the factors jshift(j)
		type(val_struct),intent(in), target :: val_sol
		type(hist_struct), intent(in):: hists_sim
		
		real(dp), intent(in) :: probj_in(:)
		real(dp), intent(out):: jshift_out(:)
		real(dp) :: vscale,nborn,jshift_prratio,cumval,updjscale,distshift
		real(dp) :: jshift0(nj)
		integer  :: ij, ik, iter
		integer :: age_hr,d_hr,ai_hr,ali_hr,ei_hr,i,ii,idi, it, beti=1
		real(dp), allocatable :: j_val_ij(:,:)
		allocate(j_val_ij(Nsim,nj))

		updjscale  = 0.1_dp
		jshift_out = 0._dp !initialize
		jshift0 = 0._dp
		j_val_ij = -1.e10_dp
		nborn = 0._dp
		vscale= 0._dp
		! first load the value functions for the first period
		it = 1
		do i = 1,Nsim
			if(hists_sim%age_hist(i,it) > 0 ) then !they're alive in the first period
				nborn = 1._dp + nborn
				!set the state
				ali_hr 	= hists_sim%al_int_hist(i,it)
				age_hr 	= 1
				d_hr 	= 1
				ai_hr	= 1
				ei_hr	= 1
				do ij = 1,nj
					j_val_ij(i,ij) = 0.
					do idi = 1,ndi ! expectation over delta
						j_val_ij(i,ij) = val_sol%V((ij-1)*nbi+beti,(idi-1)*nal+ali_hr,d_hr,ei_hr,ai_hr, &
									& hists_sim%z_jt_macro(it),age_hr)*delwt(idi,ij) + j_val_ij(i,ij)
					enddo
					vscale = j_val_ij(i,ij) + vscale
				enddo
			endif
		enddo
		vscale = vscale/nborn/dble(nj)

		!iterate on shift_k
		do iter=1,maxiter*nj
			distshift = 0._dp
			do ij=1,nj
				jshift_prratio = 0.
				do i =1,Nsim
					if(j_val_ij(i,ij) > -1.e9_dp) then
						cumval = 0.
						do ik=1,nj
							cumval = dexp( (j_val_ij(i,ik)+jshift0(ik))/(amenityscale*vscale) ) + cumval
						enddo
						jshift_prratio = dexp(j_val_ij(i,ij)/(amenityscale*vscale))/cumval + jshift_prratio
					endif
				enddo
				jshift_prratio = probj_in(ij)/jshift_prratio*nborn
				jshift_out(ij) = dlog(jshift_prratio)*amenityscale*vscale
				distshift = dabs(jshift_out(ij) - jshift0(ij)) + distshift
				jshift0(ij) = updjscale*jshift_out(ij) + (1._dp - updjscale)*jshift0(ij)
			enddo
			if (distshift<1.e-5_dp) then
				exit
			endif
		enddo

		deallocate(j_val_ij)
		
	end subroutine

	subroutine iter_zproc(val_sol, hists_sim, prob_hist_target)
		integer :: it, ij, ziter
		type(val_struct), target :: val_sol
		type(hist_struct), target:: hists_sim

		real(dp), intent(in) :: prob_hist_target(:,:)
		real(dp) :: zj_here, zscale_mean,zscale_mean1
		real(dp) :: dummy_params(3)
		!integer :: del_i_int(:)
		integer :: flag
		real(dp) :: resid_hi, resid_lo, z_hi,z_lo,zj_opt,zscale_dist, zscale_tol=1.e-4_dp
		real(dp), allocatable :: zscale1(:),zscale0(:),zgrid1(:,:),zgrid0(:,:)
		dummy_params = 0.

		allocate(zscale1(nj))
		allocate(zscale0(nj))
		allocate(zgrid0(nz,nj))
		allocate(zgrid1(nz,nj))
		mod_val_sol => val_sol
		mod_hists_sim => hists_sim
		zscale1 = 0.
		zscale_mean = 0.
		zgrid0 = zgrid
		do ij=1,nj
			zscale0(ij) = zscale(ij) !begin setting it equal to the global
			zscale_mean = zscale(ij)/dble(nj) + zscale_mean
		enddo

		if(print_lev>=2) call mat2csv(zgrid0,"zgrid0.csv")
		do ziter = 1,maxiter
			do ij=1,nj
				mod_ij_obj = ij
				it = struc_brk*itlen
				mod_it = it
				
				mod_prob_target = occprbrk(ij) ! i havd to replace this

				! bounds on the area in which actually solved the problem
				z_hi = zgrid(nz,ij) - zgrid(nz/2+1,ij)
				z_lo = zgrid(nz/2+1,ij) - zgrid(nz,ij)
					
				! this adjusts it to be the bottom point for all grid points
				resid_lo = obj_zj_wrap(z_lo,dummy_params)
				! be at the top point for all grid points
				resid_hi = obj_zj_wrap(z_hi,dummy_params)

				if(resid_hi*resid_lo < 0.) then
					zj_opt = zbrent(obj_zj_wrap,z_lo,z_hi,&
						& dummy_params ,zscale_tol/100.,flag)
				elseif( dabs(resid_hi)<dabs(resid_lo) ) then
					zj_opt = z_hi
				else
					zj_opt = z_lo
				endif
				zscale1(ij) = zj_opt
			enddo
			zscale_mean1 = 0.
			do ij=1,nj
				zscale_mean1 = zscale1(ij)/dble(nj) + zscale_mean1
			enddo
			zscale_dist =0.
			do ij=1,nj
				zscale1(ij) = zscale1(ij) -zscale_mean1 + zscale_mean !re-set the mean of zscale1
				zscale_dist = (zscale1(ij)-zscale(ij))**2 + zscale_dist
			enddo
			zscale_mean1 = 0.
			do ij=1,nj
				zscale(ij) = upd_zscl*zscale1(ij) + (1._dp-upd_zscl)*zscale(ij)
				zscale_mean1 = zscale(ij)/dble(nj) + zscale_mean1 !should not be changing
			enddo
			do ij=1,nj
				zscale(ij) = zscale(ij) -zscale_mean1 + zscale_mean ! hold the mean fixed: prevent drifting (necessary?)
			enddo

			! this will reassign the grid using zscale
			call newtfpgrid(zgrid0,zscale,zgrid1, val_sol)
			if(print_lev>=2) call mat2csv(zgrid1,"zgrid1.csv")
			zgrid = zgrid1

			if (zscale_dist .lt. zscale_tol) then
				exit
			endif
		enddo

		if(print_lev>=1) call mat2csv(zgrid1,"zgrid1.csv")

		deallocate(zscale0,zscale1,zgrid0,zgrid1)
		
	end subroutine iter_zproc
	

	subroutine cal_dist(paramvec, errvec)
		! the inputs are the values of parameters we're moving in paramvec
		! the outputs are deviations from targets
		! 1/ persistence of occupation productivity shock
		! 2/ standard deviation of occupation productivity shock
		! 3/ dispersion of gumbel shock - amenityscale
		! 
		
		real(dp) :: paramvec(:), errvec(:)

		type(val_struct) :: val_sol
		type(pol_struct) :: pol_sol
		type(hist_struct):: hists_sim
		type(moments_struct):: moments_sim
		real(dp) :: condstd_tsemp
		integer :: ij

		zrho = paramvec(1)
		zsig = paramvec(2)
		amenityscale = paramvec(3)

		call settfp()
		call alloc_econ(val_sol,pol_sol,hists_sim)

		! set up economy and solve it
		if(verbose >2) print *, "Solving the model"	
		call sol(val_sol,pol_sol)
		if(verbose >2) print *, "Solving for z process"
		call draw_shocks(hists_sim)
		call jshift_sol(val_sol, hists_sim, occsz0, jshift)
		call iter_zproc(val_sol, hists_sim, mod_prob_hist_tgt)
		
		if(verbose >2) print *, "Simulating the model"	
		call sim(val_sol, pol_sol, hists_sim)
		if(verbose >2) print *, "Computing moments"
		call moments_compute(hists_sim,moments_sim)

		condstd_tsemp = 0.
		do ij = 1,nj
			condstd_tsemp = moments_sim%ts_emp_coefs(ij+1) *occsz0(ij)+ condstd_tsemp
		enddo

		errvec(1) =  moments_sim%ts_emp_coefs(1) - emp_persist! auto-correlation
		!errvec(2) =  


		call dealloc_econ(val_sol,pol_sol,hists_sim)

	end subroutine cal_dist

	subroutine  dfovec(Nin,Nout,paramvec,errvec)
		integer :: Nin, Nout
		real(dp) :: paramvec(:), errvec(:)

		call cal_dist(paramvec,errvec)
		
	end subroutine dfovec

end module find_params



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

	implicit none


	!************************************************************************************************!
	! Counters and Indicies
	!************************************************************************************************!

		integer  :: id, it, ij, ibi, ial, iz, narg_in, wo, status

	!************************************************************************************************!
	! Other
	!************************************************************************************************!
		real(dp)	:: wagehere,utilhere, junk, param0(2),err0(2)
		integer, allocatable :: iz_jt_in(:)
	!************************************************************************************************!
	! Structure to communicate everything
		type(val_struct) :: val_sol
		type(pol_struct) :: pol_sol
		type(hist_struct):: hists_sim
		type(moments_struct):: moments_sim
		
	! Timers
		integer :: c1,c2,cr,cm
		real(dp) :: t1,t2
		
	!************************************************************************************************!
	! Allocate phat matrices
	!************************************************************************************************!
	call alloc_econ(val_sol,pol_sol,hists_sim)
	

	moments_sim%alloced = 0

	narg_in = iargc()


	call setparams()
	agrid(1) = .05*(agrid(1)+agrid(2))
	if(print_lev >= 2) then
		! plot out a bunch of arrays for analyzing VFs, etc
		wo = 0
		call vec2csv(agrid,"agrid.csv",wo)
		call vec2csv(delgrid,'delgrid.csv',wo)
		call vec2csv(alfgrid,'alfgrid.csv',wo)
		call vec2csv(egrid,'egrid.csv',wo)
		call mat2csv(zgrid,'zgrid.csv',wo)
		call veci2csv(dgrid,'dgrid.csv',wo)
		call veci2csv(agegrid,'agegrid.csv',wo)		
		call mat2csv(piz(:,:),"piz.csv",wo)
		call mat2csv(pialf,"pial.csv",wo)
		
		wo=0
		do it = 1,TT-1
			do ij = 1,ndi
				call mat2csv(pid(:,:,ij,it),"pid.csv",wo)
				if(wo==0) wo =1
			enddo
		enddo

		open(1, file="wage_dist.csv")
		ibi = 1
		iz  = 3
		do it = 1,TT-1
			do ial =1,nal
				do id = 1,nd-1
					wagehere = wage(beti(ibi),alfgrid(ial),id,zgrid(iz,ij),it)
					write(1, "(G20.12)", advance='no') wagehere
				enddo
				id = nd
				wagehere = wage(beti(ibi),alfgrid(ial),id,zgrid(iz,ij),it)
				write(1,*) wagehere
			enddo
			write(1,*) " "! trailing space
		enddo	
		close(1)
		open(1, file="util_dist.csv")
		junk =0.
		ibi =1
		iz  =2
		do it = 1,TT-1
			do ial =1,nal
				do id = 1,nd-1
					wagehere = wage(beti(ibi),alfgrid(ial),id,zgrid(iz,ij),it)
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
	util_const = - junk - util_const

	! draw the initial set of targets for z
	allocate(mod_prob_hist_tgt(Tsim,nj))
	allocate(iz_jt_in(Tsim))
	
	mod_prob_hist_tgt = 0.

	call draw_zjt(iz_jt_in, 12092015, status)
	
	do it=1,Tsim
		junk = 0.
		do ij=1,nj
			junk = zgrid( iz_jt_in(it) , ij )/(dble(nj)) + junk
		enddo
		do ij=1,nj!test probabilities are deviations from uniform 
			mod_prob_hist_tgt(it,ij) = (zgrid( iz_jt_in(it) , ij ) -junk) + (1._dp/dble(nj) )
		enddo
	enddo
	call cal_dist(param0,err0)
	

	if(verbose >2) print *, "Solving the model"
	if(verbose >2) then
		call system_clock(count_rate=cr)
		call system_clock(count_max=cm)
		call CPU_TIME(t1)
		call SYSTEM_CLOCK(c1)
	endif
	call sol(val_sol,pol_sol)
	if(verbose > 2) then
		call CPU_TIME(t2)
		call SYSTEM_CLOCK(c2)
		print *, "System Time", dble(c2-c1)/dble(cr)
		print *, "   CPU Time", (t2-t1)
	endif

	
	if(verbose >2) print *, "Simulating the model"	
	call sim(val_sol, pol_sol, hists_sim)
	if(verbose >2) print *, "Computing moments"
	call moments_compute(hists_sim,moments_sim)


!    .----.   @   @
!   / .-"-.`.  \v/
!   | | '\ \ \_/ )
! ,-\ `-.' /.'  /
!'---`----'----'
	!****************************************************************************!
	! IF you love something.... 
	!****************************************************************************!
	call dealloc_econ(val_sol,pol_sol,hists_sim)
	deallocate(mod_prob_hist_tgt)
	deallocate(iz_jt_in)
End PROGRAM






