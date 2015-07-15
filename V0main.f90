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
! compiler line: gfortran -ffree-line-length-none -g V0para.f90 V0main.f90 -o V0main.out 
module helper_funs
	
	use V0para
	implicit none
	
	!**********************************************************!
	!Public Policy Functions
	!	1)UI(e)		Unemployment Insurance
	!	2)SSDI(e,a)	DI (can be asset tested)
	!	3)SSI(e)	Social Sec. Retirement	
	!Utility Functions
	!	4)u(c,d,w)	w=1 if working; 2 if not
	!Earnings Index
	!	5)Hearn(t,e,w)	t=age, e=past index, w=current wage
	!Wage Function
	!	6) wage(bi,ai,d,z,t) indiv(beta,alf),disability,tfp,age
	!Locate Function
	!	7)finder(xx,x)	xx is grid, x is point, returns higher bound
	!Writing Subroutines
	!	8)  mat2csv(A,fname,append)   A=matrix, fname=file, append={0,1}
	!	9)  mati2csv(A,fname,append)  A=matrix, fname=file, append={0,1}: A is integer
	!	10) vec2csv(A,fname,append)   A=matrix, fname=file, append={0,1}
	!	11) veci2csv(A,fname,append)   A=matrix, fname=file, append={0,1}: A is integer
	!User-Defined Types (Structures) for value functiIons and policy functions
	!	a) val_struct: VR, VD, VN, VW, VU, V
	!	b) pol_struct: aR, aD,aU,aN,aW,gapp,gwork, gapp_dif,gwork_dif
	!**********************************************************!
	
	
	!------------------------------------------------------------------
	! a) val_struct: VR, VD, VN, VW, VU, V
	!------------------------------------------------------------------
	type val_struct
		real(8), allocatable:: 	VR(:,:,:), &		!Retirement
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
		
		real(8), allocatable ::	gapp_dif(:,:,:,:,:,:,:), &
					gwork_dif(:,:,:,:,:,:,:) ! latent value of work/apply
		integer, allocatable ::	aR(:,:,:), aD(:,:,:,:), aU(:,:,:,:,:,:,:), &
					aN(:,:,:,:,:,:,:), aW(:,:,:,:,:,:,:)
		integer, allocatable ::	gapp(:,:,:,:,:,:,:), &

					gwork(:,:,:,:,:,:,:) !integer choice of apply/work
		integer :: alloced
	end type

	type moments_struct
		real(8), allocatable :: work_coefs(:), di_coefs(:)
		real(8) :: di_rate, work_rate,accept_rate
		integer :: alloced

	end type 

	type hist_struct
		real(8), allocatable :: work_dif_hist(:,:), app_dif_hist(:,:) !choose work or not, apply or not -- latent value
		real(8), allocatable :: wage_hist(:,:) !realized wages
		real(8), allocatable :: obsX_hist(:,:) ! a bunch of explanitory variables stacked on each other
		integer :: alloced
	end type
	
	contains

	!------------------------------------------------------------------------
	!1)UI(e): Unemployment Insurance
	!----------------
	function UI(ein)	
	!----------------
	
		real(8), intent(in)	:: ein 
		real(8) 		:: UI
		!I'm using a replacement rate of 30% for now, can be fancier
		UI = ein*UIrr

	end function
	!------------------------------------------------------------------------
	!2) DI (can be asset tested)
	!---------------------
	function SSDI(ein)
	!---------------------
	 
	
		real(8), intent(in)	:: ein
		real(8) 		:: SSDI
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
	 
		real(8), intent(in)	:: ein
		real(8)			:: SSI
		
		!Follows Pistafferi & Low '12
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
	! 4) Utility Function
	!--------------------
	function util(cin,din,win)
	!--------------------
	 
		real(8), intent(in)	:: cin
		integer, intent(in)	:: din, win
		real(8)			:: util
		
		
		if ((win .gt. 1) .or. (din .gt. 1)) then
			! win >=2 => disutility from work
			util = ((cin*dexp(theta*dble(din-1)+eta*dble(win-1)))**(1.-gam))/(1.-gam)
		!elseif( din .gt. 1) then
		!	util = ((cin*dexp(theta*dble(din-1)))**(1.-gam))/(1.-gam)
		else 
			util = (cin**(1.-gam))/(1.-gam)
		end if

	end function

	!------------------------------------------------------------------------
	! 5) Earnings Index Function
	!--------------------
	function Hearn(tin,ein,win)
	!--------------------
	
		real(8), intent(in)	:: win
		integer, intent(in)	:: ein, tin
		real(8)			:: Hearn
		
		if (tin .EQ. 1) THEN
		  Hearn = (egrid(ein)*tlength*youngD/2+win)/(tlength*(youngD*oldD*oldN)) 
		else
		  Hearn = (egrid(ein)*tlength*(youngD+(tin-1)*oldD+oldD/2)+win)/(tlength*(youngD*oldD*oldN))
		end if
	end function

	!------------------------------------------------------------------------
	! 6) Wage Function
	!--------------------
	function wage(biin,aiin,din,zin,tin)
	!--------------------
	
		real(8), intent(in)	:: biin, aiin, zin
		integer, intent(in)	:: din, tin
		real(8)			:: wage

		wage = dexp( biin*zin + aiin+wd(din)+wtau(tin) ) 

	end function

	!------------------------------------------------------------------------
	! 7) Locate Function
	!--------------------
	function finder(xx,x)
	!--------------------

		real(8), dimension(:), intent(IN) :: xx
		real(8), intent(IN) :: x
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


end module helper_funs


module sol_sim

	use V0para
	use helper_funs
! module with subroutines to solve the model and simulate data from it 

	implicit none
	
	contains

	subroutine sol(val_funs, pol_funs)

		implicit none
	
		type(val_struct), intent(inout), target :: val_funs
		type(pol_struct), intent(inout), target :: pol_funs	

	!************************************************************************************************!
	! Counters and Indicies
	!************************************************************************************************!

		integer  :: i, j, ia, ie, id, it, iaa,iaa1, iaa1app,iaa1napp, anapp,aapp, apol, ibi, ial, ij , idi, izz, iaai,  &
			    iee1, iee2, iz, iw,wo, iter
		integer, dimension(5) :: maxer_i
	
		!************************************************************************************************!
		! Value Functions- Stack z-risk j and indiv. exposure beta_i
		!************************************************************************************************!
		real(8)  	  	:: Vtest1, Vtest2, utilhere, Vapp, Vc1, Vnapp, maxer_v, smthV,smthV0param, &
					&	iee1wt, maxVNV0
		real(8), allocatable	:: maxer(:,:,:,:,:)
		real(8), allocatable :: VR0(:,:,:), &			!Retirement
					VD0(:,:,:,:), &			!Disabled
					VN0(:,:,:,:,:,:,:), &	!Long-term Unemployed
					VW0(:,:,:,:,:,:,:), &	!Working
					VU0(:,:,:,:,:,:,:), &	!Unemployed
					V0(:,:,:,:,:,:,:)	!Participant
				
		real(8), pointer ::	VR(:,:,:), &			!Retirement
					VD(:,:,:,:), &			!Disabled
					VN(:,:,:,:,:,:,:), &	!Long-term Unemployed
					VW(:,:,:,:,:,:,:), &	!Working
					VU(:,:,:,:,:,:,:), &	!Unemployed
					V(:,:,:,:,:,:,:)	!Participant
	
		real(8), pointer ::	gapp_dif(:,:,:,:,:,:,:), gwork_dif(:,:,:,:,:,:,:) ! latent value of work/apply
	
		integer, pointer ::	aR(:,:,:), aD(:,:,:,:), aU(:,:,:,:,:,:,:), &
					aN(:,:,:,:,:,:,:), aW(:,:,:,:,:,:,:)
		integer, pointer ::	gapp(:,:,:,:,:,:,:), &
					gwork(:,:,:,:,:,:,:)
	
		!************************************************************************************************!
		! Other
		!************************************************************************************************!
			real(8)	:: wagehere,chere, junk,summer, eprime, yL, yH, emin, emax, VUhere, VWhere
		!************************************************************************************************!
		
		!************************************************************************************************!
		! Allocate phat matrices
		!************************************************************************************************!
		! (disability extent, earn hist, assets)
		allocate(VR0(nd,ne,na))
		allocate(VD0(nd,ne,na,TT))
		allocate(VN0(nj*nbi,ndi*nai,nd,ne,na,nz,TT))
		allocate(VU0(nj*nbi,ndi*nai,nd,ne,na,nz,TT))
		allocate(VW0(nj*nbi,ndi*nai,nd,ne,na,nz,TT))
		allocate(V0(nj*nbi,ndi*nai,nd,ne,na,nz,TT))
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

!		! (disability extent, earn hist, assets)
!		allocate(VR(nd,ne,na))
!		allocate(aR(nd,ne,na))
!		! (disability extent, earn hist, assets, age)
!		allocate(VD(nd,ne,na,TT))
!		allocate(aD(nd,ne,na,TT-1))

!		! (occupation X ind exposure, ind disb. risk X ind. wage, disab. extent, earn hist, assets, agg shock, age)
!		allocate(VN(nj*nbi,ndi*nai,nd,ne,na,nz,TT))
!		allocate(VU(nj*nbi,ndi*nai,nd,ne,na,nz,TT))
!		allocate(VW(nj*nbi,ndi*nai,nd,ne,na,nz,TT))
!		allocate(V(nj*nbi,ndi*nai,nd,ne,na,nz,TT))
!		allocate(aN(nj*nbi,ndi*nai,nd,ne,na,nz,TT-1))
!		allocate(aW(nj*nbi,ndi*nai,nd,ne,na,nz,TT-1))
!		allocate(aU(nj*nbi,ndi*nai,nd,ne,na,nz,TT-1))
!		allocate(gwork(nj*nbi,ndi*nai,nd,ne,na,nz,TT-1))
!		allocate(gapp(nj*nbi,ndi*nai,nd,ne,na,nz,TT-1))

!		allocate(gapp_dif(nj*nbi,ndi*nai,nd,ne,na,nz,TT))
!		allocate(gwork_dif(nj*nbi,ndi*nai,nd,ne,na,nz,TT))

		allocate(maxer(na,nz,ne,nd,nai))
		emin = minval(egrid)
		emax = maxval(egrid)

		!************************************************************************************************!
		! Caculate things that are independent of occupation/person type
		!	1) Value of Retired:  VR(d,e,a)
		!	2) Value of Disabled: VD(d,e,a)
		
	!1) Calculate Value of Retired: VR(d,e,a)
		!d in{1,2,3}  : disability extent
		!e inR+       :	earnings index
		!a inR+	      : asset holdings
		
		iw = 1 ! not working
		!VFI with good guess
		!Initialize
		junk = (1.0+(beta*ptau(TT)*R**(1.0-gam))**(-1.0/gam))**(-1.0)
		id=1
		do ie=1,ne
		do ia=1,na
			VR0(id,ie,ia) = util(SSI(egrid(ie))+R*agrid(ia),id,iw)* (1./(1.-beta*ptau(TT)))
		enddo
		enddo
		if(print_lev >3) then
			i = 1
			call vec2csv(VR0(i,i,:),"VR0.csv",0)
		endif		
		iter = 1
		do WHILE (iter<=maxiter)
			summer = 0
			id =1
		  	do ie=1,ne
				iaa1 = 1
			  	do ia=1,na
					Vtest1 = -1e6
					apol = iaa1
					do iaa=iaa1,na
						chere = SSI(egrid(ie))+R*agrid(ia)-agrid(iaa)
						if( chere .gt. 0.) then !ensure positive consumption
							Vc1 = beta*ptau(TT)*VR0(id,ie,iaa)
							Vtest2 = util(chere ,id,iw) + Vc1

							if(Vtest2 > Vtest1  .or. iaa .eq. iaa1 ) then !always replace on the first loop
								Vtest1 = Vtest2
								apol = iaa
							elseif(iaa > apol+iaa_hiwindow) then ! gone some steps w/o new max
								exit
							endif
						else!saved too much and negative consumtion
							exit
						endif
					enddo
					iaa1  = max(apol-iaa_lowindow,1)  	!concave, start next loop here
					VR(id,ie,ia) = Vtest1
					aR(id,ie,ia) = apol !agrid(apol)
					summer = summer+ (VR(id,ie,ia)-VR0(id,ie,ia))**2
					!Polices for disabled are the same, just scale V-function
				enddo !ia
				VR0(id,ie,:) = VR(id,ie,:)	!New guess, recall, only using id =1 
			enddo !ie
			IF (summer < Vtol) THEN
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

			

		do id=1,nd
		do ie=1,ne
		do ia=1,na
			VD(id,ie,ia,TT) = VR(id,ie,ia)
			VD0(id,ie,ia,TT) = VR(id,ie,ia)
		enddo
		enddo
		enddo



	       !----------------------------------------------------------------!
	!2) Calculate Value of Disabled: VD(d,e,a,t)	 
		!d in{1,2,3}  	   :disability extent
		!e inR+       	   :earnings index
		!a inR+	      	   :asset holdings
		!t in[1,2...TT-1]  :age

	
		!Work backwards from TT
		do it = 1,TT-1

			id = 1 ! other values are just multiples thereof
			iw = 1 ! not working
			
			!Guess will be value at t+1
			
			VD0(id,:,:,TT-it) = VD(id,:,:,TT-it+1)

			if(print_lev >3) then
				call mat2csv(VD0(id,:,:,TT-it),"VD0.csv",0)
			endif

			!Loop to find V(..,it) as fixed point
			iter=1
			do WHILE (iter<=maxiter)
				summer = 0
				id =1
				!Loop over earnings index
				do ie=1,ne
					iaa1 = 1
					!Loop over current state: assets
					do ia=1,na
						Vc1 = beta*((1-ptau(TT-it))*VD0(id,ie,apol,TT-it+1)+ptau(TT-it)*VD0(id,ie,iaa1,TT-it))
						chere = SSDI(egrid(ie))+R*agrid(ia)-agrid(iaa1)
						
						Vtest1 = -1e6
						apol = iaa1
						!Find Policy
						do iaa=iaa1,na
							chere = SSDI(egrid(ie))+R*agrid(ia)-agrid(iaa)
							if(chere >0.) then
								Vc1 = beta*((1-ptau(TT-it))*VD0(id,ie,iaa,TT-it+1)+ptau(TT-it)*VD0(id,ie,iaa,TT-it))

								Vtest2 = util(chere,id,iw)+ Vc1
								if(Vtest2>Vtest1) then
									Vtest1 = Vtest2
									apol = iaa
								elseif(iaa>apol+iaa_hiwindow) then
									exit
								endif
							else
								exit
							endif
						enddo	!iaa
						iaa1 = max(apol-iaa_lowindow,1)
						VD(id,ie,ia,TT-it) = Vtest1
						aD(id,ie,ia,TT-it) = apol !agrid(apol)
						summer = summer+ (VD(id,ie,ia,TT-it)-VD0(id,ie,ia,TT-it))**2
					enddo	!ia
				enddo	!ie		
				if(print_lev >3) then
					wo =0
					ie =1
					call veci2csv(aD(id,ie,:,TT-it),"aD.csv",wo)
					call vec2csv(VD(id,ie,:,TT-it),"VD.csv",wo)
				endif
				IF (summer < Vtol) THEN
					exit	!Converged
				EndIF
				VD0(id,:,:,TT-it) = VD(id,:,:,TT-it)	!New guess
				iter=iter+1
			enddo	!iter: V-iter loop
			! set the other value for other disability levels
			i =1 
			do id = 2,nd
				do ie =1,ne
				do ia =1,na
					VD(id,ie,ia,TT-it) = VD(i,ie,ia,TT-it)*( (dexp(theta*dble(id-1)))**(1-gam) )
				enddo
				enddo
				if(TT-it .le. TT-1) aD(id,:,:,TT-it) = aD(i,:,:,TT-it) 
			enddo
		enddo	!t loop, going backwards



		VD0 = VD
		if (print_lev > 2) then
			wo =0 
			id =1
			ie =1
			call mati2csv(aD(id,ie,:,:),"aD.csv",wo)
			call mat2csv(VD(id,ie,:,:),"VD.csv",wo)
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


	
	! Begin loop over occupations
		do ij = 1,nj
	! And betas
		do ibi = 1,nbi 
	! And individual disability type
		do idi = 1,ndi

	!************************************************************************************************!
		!Work Backwards TT-1,TT-2...1!
		do it=1,TT-1
			!----Initialize---!
			do ial=1,nai
			do id =1,nd
			do ie =1,ne
			do iz =1,nz
			do ia =1,na
				
			!Guess once, then use next period same occupation/beta as guess
			! for it = 1, should be TT-1+1 =TT -> VU,Vw,VN = VR
				VW0((ij-1)*nbi+ibi,(idi-1)*nai+ial,id,ie,ia,iz,TT-it) = VW((ij-1)*nbi+ibi,(idi-1)*nai+ial,id,ie,ia,iz,TT-it+1)
				VU0((ij-1)*nbi+ibi,(idi-1)*nai+ial,id,ie,ia,iz,TT-it) = VU((ij-1)*nbi+ibi,(idi-1)*nai+ial,id,ie,ia,iz,TT-it+1)
				VN0((ij-1)*nbi+ibi,(idi-1)*nai +ial,id,ie,ia,iz,TT-it) = VN((ij-1)*nbi+ibi,(idi-1)*nai+ial,id,ie,ia,iz,TT-it+1)
				V0((ij-1)*nbi+ibi,(idi-1)*nai+ial,id,ie,ia,iz,TT-it)= VW0((ij-1)*nbi+ibi,(idi-1)*nai+ial,id,ie,ia,iz,TT-it)

			enddo	!ia
			enddo	!iz
			enddo	!ie
			enddo 	!id
			enddo	!ial   

			!***********************************************************************************************************
			!Loop over V=max(VU,VW)	
			iter=1
			smthV0param  =1. ! will tighten this down
			do WHILE (iter<=maxiter)
				maxer = 0.
				summer = 0.	!Use to calc |V-V0|<eps
				! lots f printing every 100 iterations (mod 1 because Fortran is designed by idiots using base-1 indexing)
				if(mod(iter,100) .eq. 1) then
					print_lev = 4
				else 
					print_lev =1
				endif
			!------------------------------------------------!
			!Solve VU given guesses on VW, VN, VU and implied V
			!------------------------------------------------!  
			!$OMP  parallel do default(shared) &
			!$OMP& private(ial,id,ie,iz,iw,apol,iaa1,ia,iaa,chere,Vtest2,Vtest1,Vc1,iaai,izz,maxVNV0)
			  	do ial=1,nai	!Loop over alpha (ai)
				do id=1,nd	!Loop over disability index
			  	do ie=1,ne	!Loop over earnings index
			  	do iz=1,nz	!Loop over TFP
			  		iw = 1 !not working
					!Restart at bottom of asset grid for each of the above (ai,d,e,z)
					apol = 1
					iaa1 = 1
					!----------------------------------------------------------------
					!Loop over current state: assets
					do ia=1,na
						Vtest1 = -1e5 ! just a very bad number, does not really matter
						do iaa=iaa1,na
							chere = UI(egrid(ie))+R*agrid(ia)-agrid(iaa)
							if(chere>0.) then 
								Vtest2 = 0. !Continuation value if don't go on disability
								do izz = 1,nz	 !Loop over z'
								do iaai = 1,nai !Loop over alpha_i'
									Vc1 = (1.-ptau(TT-it))*(pphi*VN0((ij-1)*nbi+ibi,(idi-1)*nai+iaai,id,ie,iaa,izz,TT-it+1) &
										& 	+(1-pphi)*   VU0((ij-1)*nbi+ibi,(idi-1)*nai+iaai,id,ie,iaa,izz,TT-it+1) )  !Age and might go LTU
									Vc1 = ptau(TT-it)*(pphi*     VN0((ij-1)*nbi+ibi,(idi-1)*nai+iaai,id,ie,iaa,izz,TT-it) & 
										&	+(1-pphi)*   VU0((ij-1)*nbi+ibi,(idi-1)*nai+iaai,id,ie,iaa,izz,TT-it) ) + Vc1    !Don't age, maybe LTU
									Vtest2 = Vtest2 + beta*piz(iz,izz,ij)*pialf(ial,iaai)*Vc1  !Probability of alpha_i X z_i draw 
								enddo
								enddo
								Vtest2 = Vtest2 + util(chere,id,iw)
								if(Vtest2>Vtest1 .or. iaa .eq. iaa1) then !first value or optimal
									apol = iaa! set the policy
									Vtest1 = Vtest2
								elseif(iaa > apol+iaa_hiwindow) then
									exit
								endif
							else
								exit
							endif
						enddo	!iaa
						iaa1 = max(apol - iaa_lowindow,1) !concave, start next loop here
						VU((ij-1)*nbi+ibi,(idi-1)*nai+ial,id,ie,ia,iz,TT-it) = Vtest1
						aU((ij-1)*nbi+ibi,(idi-1)*nai+ial,id,ie,ia,iz,TT-it) = apol !agrid(apol)

					enddo !ia
				
				enddo !id
			  	enddo !ie
			  	enddo !iz
				enddo !ial	
			!$OMP END PARALLEL do

				if (print_lev > 3) then
					wo = 0
					do ial=1,nai	!Loop over alpha (ai)
					do id=1,nd	!Loop over disability index
				  	do ie=1,ne	!Loop over earnings index
				  	do iz=1,nz	!Loop over TFP
						call veci2csv(aU((ij-1)*nbi+ibi,(idi-1)*nai+ial,id,ie,:,iz,TT-it),"aU.csv",wo)
						call vec2csv(VU((ij-1)*nbi+ibi,(idi-1)*nai+ial,id,ie,:,iz,TT-it),"VU.csv",wo)
						if(wo == 0) wo = 1 
					enddo
					enddo
					enddo
					enddo
				endif

	
			!------------------------------------------------!
			!Solve VN given guesses on VW, VN, and implied V
			!------------------------------------------------! 
!			  	do ial=1,nai	!Loop over alpha (ai)
!				do id=1,nd	!Loop over disability index
!			  	do ie=1,ne	!Loop over earnings index
!			  	do iz=1,nz	!Loop over TFP			
!					do ia=1,na
!						VN((ij-1)*nbi+ibi,(idi-1)*nai+ial,id,ie,ia,iz,TT-it) = 	VD(id,ie,ia,TT-it)				
!					enddo
!				enddo
!				enddo
!				enddo
!				enddo



			!$OMP  parallel do default(shared) &
			!$OMP& private(ial,id,ie,iz,iw,apol,iaa1app,iaa1napp,ia,iaa,chere,Vc1,Vtest2,Vtest1,smthV,Vapp,Vnapp,aapp,anapp,iaai,izz,maxVNV0) 
			  	do ial=1,nai	!Loop over alpha (ai)
				do id=1,nd	!Loop over disability index
			  	do ie=1,ne	!Loop over earnings index
			  	do iz=1,nz	!Loop over TFP
					iw = 1 ! not working
					!******************************************************************************
					!---------------------------------------------------------------!
					!First Solve as if do NOT apply (since assets are a joint choice)
					!---------------------------------------------------------------!
					!Restart at bottom of asset grid for each of the above (ai,d,e,z)
					iaa1napp = 1
					iaa1app = 1
					Vapp  = -5e6
					Vnapp = -5e6				
					!----------------------------------------------------------------
					!Loop over current state: assets
					do ia=1,na
						
						!*******************************************
						!**********Value if apply for DI 
						Vtest1 = -1e6
						apol = iaa1app
						iaa = iaa1app
						do while (iaa <= na)
							chere = b+R*agrid(ia)-agrid(iaa)
							if(chere >0.) then
								Vtest2 = 0.
								!Continuation if apply for DI
								do izz = 1,nz	 !Loop over z'
								do iaai = 1,nai !Loop over alpha_i'
									Vc1 = (1-ptau(TT-it))*((1-xi(id))*VN0((ij-1)*nbi+ibi,(idi-1)*nai+iaai,id,ie,iaa,izz,TT-it+1) &
										& +xi(id)*VD0(id,ie,iaa,TT-it+1)) !Age and might go on DI
									Vc1 = Vc1 + ptau(TT-it)*((1-xi(id))*VN0((ij-1)*nbi+ibi,(idi-1)*nai+iaai,id,ie,iaa,izz,TT-it) &
										& +xi(id)*VD0(id,ie,iaa,TT-it))     !Don't age, might go on DI		
									Vtest2 = Vtest2 + beta*piz(iz,izz,ij)*pialf(ial,iaai)*Vc1 
								enddo
								enddo
								Vtest2 = util(chere,id,iw) - nu + Vtest2

								if (Vtest2>Vtest1  .or. iaa .eq. iaa1app) then	
									apol = iaa
									Vtest1 = Vtest2
								elseif(iaa > apol +iaa_hiwindow) then
									exit
								elseif( (Vtest2 < Vtest1) .and. (iaa <= iaa1app+1) .and. (iaa1app > 1)) then 
								!if the second and it's going down, then back up
									iaa = 1
									iaa1app =1
								endif
							elseif(iaa<= iaa1app .and. Vtest1 <= -1e5 .and. apol <= iaa1app) then
								iaa = 1 !started too much saving, go back towards zero
								iaa1app = 1
							else
								exit
							endif
							iaa = iaa + 1
						enddo !iaa
						Vapp = Vtest1
						iaa1app = max(apol -iaa_lowindow,1) !concave, start next loop here
						aapp = apol !agrid(apol)					
						if(Vapp <-1e5) then
							write(*,*) "ruh roh!"
							write(*,*) "VD: ",id,ie,iaa,TT-it
							write(*,*) "VN: ",ij,ibi,idi,iaai,id,ie,iaa,izz,TT-it
						endif

						!*******************************************
						!**************Value if do not apply for DI
						Vtest1 = -1e6
						apol = iaa1napp
						iaa=iaa1napp
						do while (iaa <= na)
							chere = b+R*agrid(ia)-agrid(iaa)
							if(chere >0.) then
								Vtest2 = 0.
								!Continuation if do not apply for DI
								do izz = 1,nz	 !Loop over z'
								do iaai = 1,nai !Loop over alpha_i'
								
								!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
								!  Make this a max s.t. can go to V0 or can go to VN0
									maxVNV0 = max(		 V0((ij-1)*nbi+ibi,(idi-1)*nai+iaai,id,ie,iaa,izz,TT-it+1), &
											& 	VN0((ij-1)*nbi+ibi,(idi-1)*nai+iaai,id,ie,iaa,izz,TT-it+1))
									Vc1 = (1-ptau(TT-it))*((1-rhho)* &
											&	VN0((ij-1)*nbi+ibi,(idi-1)*nai +iaai,id,ie,iaa,izz,TT-it+1) +rhho*maxVNV0) !Age and might go on DI
									maxVNV0 = max(		 V0((ij-1)*nbi+ibi,(idi-1)*nai+iaai,id,ie,iaa,izz,TT-it), & 
											&	VN0((ij-1)*nbi+ibi,(idi-1)*nai+iaai,id,ie,iaa,izz,TT-it))
									Vc1 = Vc1+ptau(TT-it)*((1-rhho)* & 
											&	VN0((ij-1)*nbi+ibi,(idi-1)*nai +iaai,id,ie,iaa,izz,TT-it) +rhho*maxVNV0)     !Don't age, might go on DI
									Vtest2 = Vtest2 + beta*piz(iz,izz,ij)*pialf(ial,iaai)*Vc1 
								enddo
								enddo
								Vtest2 = Vtest2 + util(chere,id,iw)
								
								if( Vtest2 > Vtest1 .or. iaa .eq. iaa1napp) then
									apol = iaa 
									Vtest1 = Vtest2 
								elseif(iaa > apol+iaa_hiwindow) then
									exit
								elseif( (Vtest2 < Vtest1) .and. (iaa <= iaa1napp+1) .and. (iaa1napp > 1)) then 
								!if the second and it's going down, then back up
									iaa = 1
									iaa1napp =1
								endif
							elseif(iaa .eq. iaa1napp .and. Vtest1 <= -1e5 .and. apol .eq. iaa1napp) then
								iaa = 1 !started too much saving, go back towards zero
								iaa1napp = 1
							else
								exit
							endif
							iaa = iaa + 1
						enddo !iaa
						iaa1napp = max(apol-iaa_lowindow,1)
						Vnapp = Vtest1 					
						aNapp = apol !agrid(apol)
						if(Vnapp <-1e5) then
							write(*,*) "ruh roh!"
							write(*,*) "VD: ",id,ie,iaa,TT-it
							write(*,*) "VN: ",ij,ibi,idi,iaai,id,ie,iaa,izz,TT-it
						endif

						!******************************************************
						!***************** Discrete choice for application
						smthV = dexp(smthV0param*Vnapp)/( dexp(smthV0param*Vnapp) +dexp(smthV0param*Vapp) )
						if( smthV .lt. 1e-5 .or. smthV .gt. 0.999999 .or. isnan(smthV)) then
							if( Vapp  > Vnapp ) smthV =0.
							if(Vnapp > Vapp  ) smthV =1.
						endif
						IF (Vapp > Vnapp) THEN
							!VN((ij-1)*nbi+ibi,(idi-1)*nai+ial,id,ie,ia,iz,TT-it) = Vapp
							aN((ij-1)*nbi+ibi,(idi-1)*nai+ial,id,ie,ia,iz,TT-it) = aapp
							gapp((ij-1)*nbi+ibi,(idi-1)*nai+ial,id,ie,ia,iz,TT-it) = 1
						ELSE !Don't apply
							!VN((ij-1)*nbi+ibi,(idi-1)*nai+ial,id,ie,ia,iz,TT-it) = Vnapp
							aN((ij-1)*nbi+ibi,(idi-1)*nai+ial,id,ie,ia,iz,TT-it) = aNapp
							gapp((ij-1)*nbi+ibi,(idi-1)*nai+ial,id,ie,ia,iz,TT-it) = 0
						EndIF

						gapp_dif((ij-1)*nbi+ibi,(idi-1)*nai+ial,id,ie,ia,iz,TT-it) = Vapp - Vnapp
						if(verbose .gt. 4) print *, Vapp - Vnapp
						
						
						VN((ij-1)*nbi+ibi,(idi-1)*nai+ial,id,ie,ia,iz,TT-it) = smthV*Vnapp + (1. - smthV)*Vapp
					enddo !ia
				enddo !iz 
				enddo !ie 
				enddo !id 
				enddo !ial
				
			!$OMP END PARALLEL do
			
				!------------------------------------------------!			
				! Done making VN
				  	
			  	if (print_lev >3) then
			  		wo = 0
				  	do ial=1,nai	!Loop over alpha (ai)
					do ie=1,ne	!Loop over earnings index
				  	do iz=1,nz	!Loop over TFP
						! matrix in disability index and assets
				  		call mat2csv(VN((ij-1)*nbi+ibi,(idi-1)*nai+ial,:,ie,:,iz,TT-it) ,"VN_it.csv",wo)
				  		call mati2csv(aN((ij-1)*nbi+ibi,(idi-1)*nai+ial,:,ie,:,iz,TT-it) ,"aN_it.csv",wo)
				  		call mat2csv(VN((ij-1)*nbi+ibi,(idi-1)*nai+ial,:,ie,:,iz,TT-it) ,"VU_it.csv",wo)
				  		call mati2csv(aN((ij-1)*nbi+ibi,(idi-1)*nai+ial,:,ie,:,iz,TT-it) ,"aU_it.csv",wo)
				  		call mati2csv(gapp((ij-1)*nbi+ibi,(idi-1)*nai+ial,:,ie,:,iz,TT-it) ,"gapp_it.csv",wo)
				  		call mat2csv(gapp_dif((ij-1)*nbi+ibi,(idi-1)*nai+ial,:,ie,:,iz,TT-it) ,"gapp_dif_it.csv",wo)

						if(wo == 0 ) wo =1		  		
					enddo !iz 
					enddo !id 
					enddo !ial 
					
					wo = 0
					do iz=1,nz
					do ie=1,ne
					do ia=1,na
						if(wo == 0 ) wo =1
						call mat2csv(gapp_dif((ij-1)*nbi+ibi,(idi-1)*nai+1:idi*nai,:,ie,ia,iz,TT-it),'dilat0_dalpha_it.csv',wo)
					enddo
					enddo
					enddo
					
						
			  	endif
		  	
			  	!update VN0
			  	do ial=1,nai	!Loop over alpha (ai)
				do id=1,nd	!Loop over disability index
			  	do ie=1,ne	!Loop over earnings index
			  	do iz=1,nz	!Loop over TFP
					do ia =1,na
						VN0((ij-1)*nbi+ibi,(idi-1)*nai+ial,id,ie,ia,iz,TT-it) = VN((ij-1)*nbi+ibi,(idi-1)*nai+ial,id,ie,ia,iz,TT-it)
					enddo	!ia
				enddo !ial
			  	enddo !id
			  	enddo !ie
			  	enddo !iz

				!------------------------------------------------!
				!Solve VW given guesses on VW, VN, and implied V
				!------------------------------------------------!
			!$OMP   parallel do default(shared) reduction(+:summer) &
			!$OMP & private(ial,id,ie,iz,apol,eprime,wagehere,iee1,iee2,egrid,iee1wt,ia,iaa,iaa1,chere,yL,yH,Vc1,utilhere,Vtest2,Vtest1,smthV,VUhere,VWhere,iaai,izz) 
			  	do ial=1,nai	!Loop over alpha (ai)
				do id=1,nd	!Loop over disability index
			  	do ie=1,ne	!Loop over earnings index
			  	do iz=1,nz	!Loop over TFP
					!Earnings evolution independent of choices.
					wagehere = wage(beti(ibi),alfi(ial),id,zgrid(iz),TT-it)
					eprime = Hearn(TT-it,ie,wagehere)
					!linear interpolate for the portion that blocks off bounds on assets
					if(eprime > emin .and. eprime < emax) then  ! this should be the same as if(eprime > minval(egrid) .and. eprime < maxval(egrid))
						iee1 = finder(egrid,eprime) !apparently this causes a race unless egrid is private
						iee2 = min(ne,iee1+1)
						iee1wt = (egrid(iee2)-eprime)/(egrid(iee2)-egrid(iee1))
					elseif( eprime <= emin  ) then 
						iee1wt = 1.
						iee1 = 1
						iee2 = 1
					else 
						iee1wt = 0.
						iee1 = ne
						iee2 = ne
					endif

					!Restart at bottom of asset grid for each of the above (ai,d,e,z)
					iaa1 = 1
					!----------------------------------------------------------------
					!Loop over current state: assets
					do ia=1,na
						Vtest1= -1.e5 ! just to initialize, does not matter
						!Find saving Policy
						iaa = iaa1 
						do while (iaa <= na)
							!Continuation value if don't go on disability
							chere = wagehere+R*agrid(ia)-agrid(iaa)
							if (chere >0.) then
								Vc1 = 0.
								do izz = 1,nz	 !Loop over z'
								do iaai = 1,nai !Loop over alpha_i'
									!Linearly interpolating on e'
									yL = (1-ptau(TT-it))*V0((ij-1)*nbi+ibi,(idi-1)*nai+iaai,id,iee1,iaa,izz,TT-it+1) & 
										& +ptau(TT-it)*V0((ij-1)*nbi+ibi,(idi-1)*nai+iaai,id,iee1,iaa,izz,TT-it)
									yH = (1-ptau(TT-it))*V0((ij-1)*nbi+ibi,(idi-1)*nai+iaai,id,iee2,iaa,izz,TT-it+1) & 
										& +ptau(TT-it)*V0((ij-1)*nbi+ibi,(idi-1)*nai+iaai,id,iee2,iaa,izz,TT-it)
									Vc1 = piz(iz,izz,ij)*pialf(ial,iaai) &
										& * (yH*(1. - iee1wt) + yL*iee1wt )&
										& + Vc1
								enddo
								enddo
								utilhere = util(chere,id,2)
								Vtest2 = utilhere + beta*Vc1 ! flow utility

								if (Vtest2>Vtest1 .or. iaa .eq. iaa1 ) then  !always replace on the first, or best
									Vtest1 = Vtest2
									apol = iaa
								elseif(iaa>apol+iaa_hiwindow) then
									exit
								elseif( (Vtest2 < Vtest1) .and. (iaa <= iaa1+1) .and. (iaa1 > 1)) then 
								!if the second and it's going down, then back up
									iaa = 1
									iaa1 =1
								endif
							elseif(iaa == iaa1 .and. Vtest1<=-5e4 .and. apol == iaa1 ) then
								!back it up because started with unfeasible savings / negative consumption
								iaa = 1
								iaa1 =1
							else 
								exit
							endif
							iaa = iaa+1
						enddo	!iaa
						iaa1 = max(iaa -iaa_lowindow,1)	!concave? start next loop here
						VW((ij-1)*nbi+ibi,(idi-1)*nai+ial,id,ie,ia,iz,TT-it) = Vtest1
						aW((ij-1)*nbi+ibi,(idi-1)*nai+ial,id,ie,ia,iz,TT-it) = agrid(apol)


						!------------------------------------------------!
						!Calculate V with solved vals of VW and VU -- i.e. can quit into unemployment
						!------------------------------------------------!
						VUhere = VU((ij-1)*nbi+ibi,(idi-1)*nai+ial,id,ie,ia,iz,TT-it)
						VWhere = VW((ij-1)*nbi+ibi,(idi-1)*nai+ial,id,ie,ia,iz,TT-it)

						if (VWhere>VUhere) then
							gwork((ij-1)*nbi+ibi,(idi-1)*nai+ial,id,ie,ia,iz,TT-it) = 1
						else
							gwork((ij-1)*nbi+ibi,(idi-1)*nai+ial,id,ie,ia,iz,TT-it) = 0
						endif
						smthV = dexp( smthV0param *VWhere  ) &
							& /( dexp(smthV0param * VWhere) + dexp(smthV0param * VUhere ) )
						if (smthV <1e-5 .or. smthV>.999999 .or. isnan(smthV) ) then
							if(VWhere .gt. VUhere)	smthV = 1. 
							if(VWhere .lt. VUhere)	smthV = 0. 							
						endif
						V((ij-1)*nbi+ibi,(idi-1)*nai+ial,id,ie,ia,iz,TT-it) = smthV*VWhere &
								& + (1.-smthV)*VUhere

						gwork_dif((ij-1)*nbi+ibi,(idi-1)*nai+ial,id,ie,ia,iz,TT-it) = VWhere - VUhere
						
						summer = summer+ & 
							& (V((ij-1)*nbi+ibi,(idi-1)*nai+ial,id,ie,ia,iz,TT-it)-V0((ij-1)*nbi+ibi,(idi-1)*nai+ial,id,ie,ia,iz,TT-it))**2
						
						maxer(ia,iz,ie,id,ial) = (V((ij-1)*nbi+ibi,(idi-1)*nai+ial,id,ie,ia,iz,TT-it)-V0((ij-1)*nbi+ibi,(idi-1)*nai+ial,id,ie,ia,iz,TT-it))**2

					enddo !ia
				enddo !iz
			  	enddo !ie
			  	enddo !id
			  	enddo !ial

			!$OMP  END PARALLEL do
			
				maxer_v = maxval(maxer)
				maxer_i = maxloc(maxer)

			  	!update VW0
			  	do ial=1,nai	!Loop over alpha (ai)
				do id=1,nd	!Loop over disability index
			  	do ie=1,ne	!Loop over earnings index
			  	do iz=1,nz	!Loop over TFP
				do ia =1,na
					VW0((ij-1)*nbi+ibi,(idi-1)*nai+ial,id,ie,ia,iz,TT-it) = VW((ij-1)*nbi+ibi,(idi-1)*nai+ial,id,ie,ia,iz,TT-it)
				enddo !ia		  	
				enddo !ial
			  	enddo !id
			  	enddo !ie
			  	enddo !iz

				if (print_lev >3) then
					wo = 0
					do ial=1,nai	!Loop over alpha (ai)
					do ie=1,ne	!Loop over earnings index
					do iz=1,nz	!Loop over TFP
						! matrix in disability and assets
						call mat2csv(VW((ij-1)*nbi+ibi,(idi-1)*nai+ial,:,ie,:,iz,TT-it) ,"VW_it.csv",wo)
						call mati2csv(aW((ij-1)*nbi+ibi,(idi-1)*nai+ial,:,ie,:,iz,TT-it) ,"aW_it.csv",wo)
						call mati2csv(gwork((ij-1)*nbi+ibi,(idi-1)*nai+ial,:,ie,:,iz,TT-it) ,"gwork_it.csv",wo)
						call mat2csv(gwork_dif((ij-1)*nbi+ibi,(idi-1)*nai+ial,:,ie,:,iz,TT-it) ,"gwork_idf_it.csv",wo)
						if(wo==0) wo =1
					enddo !iz 
					enddo !ie 
					enddo !ial 	
				endif


			  	!update V0
			  	do ial=1,nai	!Loop over alpha (ai)
				do id=1,nd	!Loop over disability index
			  	do ie=1,ne	!Loop over earnings index
			  	do iz=1,nz	!Loop over TFP
				do ia =1,na
					V0((ij-1)*nbi+ibi,(idi-1)*nai+ial,id,ie,ia,iz,TT-it) = V((ij-1)*nbi+ibi,(idi-1)*nai+ial,id,ie,ia,iz,TT-it)
				enddo !ia		  	
				enddo !ial
			  	enddo !id
			  	enddo !ie
			  	enddo !iz


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
				! End of iter iteration
				!------------------------------------------------!
				!Check |V-V0|<eps
				!------------------------------------------------!
				write(*,*) summer, iter, ij, ibi, idi, it
				write(*,*) maxer_v, maxer_i(1), maxer_i(2), maxer_i(3), maxer_i(4), maxer_i(5)
				
				if (summer < Vtol ) then
					exit !Converged
				endif

				iter=iter+1
				smthV0param = smthV0param*1.5 !tighten up the discrete choice
			enddo	!iter: V-iter loop
	!WRITE(*,*) ij, ibi, idi, it
		enddo	!t loop, going backwards

		enddo	!idi
		enddo	!ibi
		enddo	!ij
			! reset V0
		!	do ial=1,nai	!Loop over alpha (ai)
		!	do id=1,nd	!Loop over disability index
		! 	do ie=1,ne	!Loop over earnings index
		!	do iz=1,nz	!Loop over TFP
		!  	do ia=1,na
		!		V0((ij-1)*nbi+ibi,(idi-1)*nai+ial,id,ie,ia,iz,TT-it) = V((ij-1)*nbi+ibi,(idi-1)*nai+ial,id,ie,ia,iz,TT-it)
		!	enddo !ial
		!	enddo !id
		!	enddo !ie
		!	enddo !iz
		! 	enddo !ia
		
		! this plots work-rest and di application on the cross product of alphai and deltai and di
		ibi = 1
		idi = 1
		wo  = 0
		do iz=1,nz
		do ie=1,ne
		do ia=1,na
			do it=1,TT-1
		!				nj*nbi,ndi*nai,nd,ne,na,nz,TT-1
				call mati2csv(gapp(ibi,(idi-1)*nai+1:idi*nai,:,ie,ia,iz,TT-it),'dipol_dalpha.csv',wo)
				call mati2csv(gwork(ibi,(idi-1)*nai+1:idi*nai,:,ie,ia,iz,TT-it),'workpol_dalpha.csv',wo)

				call mat2csv(gapp_dif(ibi,(idi-1)*nai+1:idi*nai,:,ie,ia,iz,TT-it),'dilat_dalpha.csv',wo)
				call mat2csv(gwork_dif(ibi,(idi-1)*nai+1:idi*nai,:,ie,ia,iz,TT-it),'worklat_dalpha.csv',wo)
				if(wo ==0 ) wo =1
			enddo
		enddo
		enddo
		enddo


		call vec2csv(dtype,'DriskGrid.csv',0)
		call vec2csv(alfi(:),'AlfGrid.csv',0)
		call vec2csv(occz(:),'ZriskGrid.csv',0)
		call vec2csv(agrid(:),'Agrid.csv',0)

!		val_funs%VR = VR
!		pol_funs%aR = aR
!		val_funs%VD = VD
!		pol_funs%aD = aD
!		val_funs%VN = VN
!		val_funs%VU = VU
!		val_funs%VW = VW
!		val_funs%V = V 
!		pol_funs%aN = aN
!		pol_funs%aW = aW
!		pol_funs%aU = aU
!		pol_funs%gwork = gwork
!		pol_funs%gapp = gapp

!		pol_funs%gapp_dif = gapp_dif
!		pol_funs%gwork_dif = gwork_dif



		deallocate(maxer)
		deallocate(VR0,VD0,VN0,VU0,VW0,V0)
!		deallocate(VR,VD,VN,VU,VW,V)
!		deallocate(aR,aD,aN,aW,aU,gwork,gapp,gapp_dif,gwork_dif)

	end subroutine sol 


	subroutine draw_deli(del_i,del_i_int, seed0, success)
	! draws depreciation rates and indices on the delta grid (i.e at the discrete values)
		implicit none

		integer, intent(in) :: seed0
		integer, intent(out) :: success
		real(8), dimension(:) :: del_i
		integer, dimension(:) :: del_i_int
		integer :: ss, Nsim, di_int,m,i
		real(8) :: dtypeL, dtypeH,dtype_i
		integer, dimension(100) :: bdayseed

		call random_seed(size = ss)
		forall(m=1:ss) bdayseed(m) = (m-1)*100 + seed0
		call random_seed(put = bdayseed(1:ss) )

		Nsim = size(del_i)
		dtypeL = minval(dtype)
		dtypeH = maxval(dtype)

		do i=1,Nsim
			call random_number(dtype_i) ! draw uniform on 0,1
			dtype_i = dtype_i*(dtypeH-dtypeL) + dtypeL !change domain of uniform
			di_int = finder(dtype,dtype_i)
			! round up or down:
			if( (dtype_i - dtype(di_int))/(dtype(di_int+1)- dtype(di_int+1)) >0.5 ) di_int = di_int + 1
			if(del_contin .eqv. .true.) then
				del_i(i) = dtype_i
			else
				del_i(i) = dtype(di_int)
			endif
			del_i_int(i) = di_int
		enddo
		success = 1
		
	end subroutine draw_deli
	
	subroutine draw_alit(al_it,al_it_int, seed0, success)
	! draws alpha shocks and idices on the alpha grid (i.e at the discrete values)
		implicit none

		integer, intent(in) :: seed0
		integer, intent(out),optional :: success
		real(8), dimension(:,:) :: al_it
		integer, dimension(:,:) :: al_it_int
		integer :: ss, Nsim, alfi_int, t,m,i
		real(8) :: alfiL, alfiH,alfi_innov,alfi_i
		integer, dimension(100) :: bdayseed

		call random_seed(size = ss)
		forall(m=1:ss) bdayseed(m) = (m-1)*100 + seed0
		call random_seed(put = bdayseed(1:ss) )

		Nsim = size(al_it,1)
		alfiL = minval(alfi)
		alfiH = maxval(alfi)

		do i=1,Nsim

			! draw starting values
			t =1
			
			alfi_innov = random_normal() ! draw normal disturbances on 0,1
			! transform it by the ergodic distribution for the first period:
			alfi_i = alfi_innov*alfsig + alfmu

			if(alfi_i >alfiH .or. alfi_i < alfiL) success = 1+success !count how often we truncate
			!impose bounds
			alfi_i = max(alfi_i,alfiL)
			alfi_i = min(alfi_i,alfiH)
			alfi_int = finder(alfi,alfi_i)
			! round up or down:
			if( (alfi_i - alfi(alfi_int))/(alfi(alfi_int+1)- alfi(alfi_int+1)) >0.5 ) alfi_int = alfi_int + 1
			if(al_contin .eqv. .true.) then
				al_it(i,t) = alfi_i ! log of wage shock
			else
				al_it(i,t) = alfi(alfi_int) ! log of wage shock, on grid
			endif
			al_it_int(i,t) = alfi_int
			
			! draw sequence:

			do t=2,Tsim
				alfi_innov = random_normal()
				alfi_innov = alfi_innov*alfsig*(1.-alfrho**2) ! mean 0 with conditional standard deviation implied by (alfsig, alfrho)
				alfi_i = alfrho*alfi_i + (1-alfrho)*alfmu + alfi_innov
				if(alfi_i >alfiH .or. alfi_i < alfiL) success = 1+success !count how often we truncate
				alfi_i = max(alfi_i,alfiL)
				alfi_i = min(alfi_i,alfiH)
				alfi_int = finder(alfi,alfi_i)
				if( (alfi_i - alfi(alfi_int))/(alfi(alfi_int+1)- alfi(alfi_int+1)) >0.5 ) alfi_int = alfi_int + 1
				if(al_contin .eqv. .true.) then
					al_it(i,t) = alfi_i ! log of wage shock
				else
					al_it(i,t) = alfi(alfi_int) ! log of wage shock, on grid
				endif
				al_it_int(i,t) = alfi_int					

			enddo
		enddo
		if(success > 0.2*Nsim*Tsim) success = 0
		if(success <= 0.2*Nsim*Tsim) success = 1
		
	end subroutine draw_alit

	subroutine draw_ji(j_i,seed0, success)
		implicit none
		integer	:: j_i(:)
		integer	:: i,m,ss
		integer,intent(in)  :: seed0
		integer,intent(out) :: success
		integer, dimension(100) :: bdayseed
		real(8)	:: Njcumdist(nj)
		real(8) :: draw_i

		i =1
		Njcumdist = Njdist(i)
		if(nj>1) then
			do i=2,nj
				Njcumdist(i) = Njdist(i) + Njcumdist(i-1)
			enddo
		
			call random_seed(size = ss)
			forall(m=1:ss) bdayseed(m) = (m-1)*100 + seed0
			call random_seed(put = bdayseed(1:ss) )
		
			do i=1,Nsim
				call random_number(draw_i)
				j_i(i) = finder(Njcumdist,draw_i)
			enddo

		else
			j_i = 1
		endif


		success = 1
		
	end subroutine draw_ji

	subroutine draw_zjt(z_jt,z_jt_int, j_i, seed0, success)
		implicit none

		real(8),intent(out) :: z_jt(:,:)
		integer,intent(out) :: z_jt_int(:,:)
		integer,intent(in)  :: j_i(:)
		integer	:: it,i,ij,iz,izp,m,ss, z_jt_t
		integer,intent(in) :: seed0
		integer,intent(out) :: success
		integer, dimension(100) :: bdayseed
		real(8) :: z_innov
		real(8), allocatable :: cumpi_j(:,:)
		integer, allocatable :: z_jt_macro(:,:) !this will be the panel across occupations -> z_jt by i's j

		allocate(z_jt_macro(nj,Tsim))
		allocate(cumpi_j(nz,nz+1))
		
		call random_seed(size = ss)
		forall(m=1:ss) bdayseed(m) = (m-1)*100 + seed0
		call random_seed(put = bdayseed(1:ss) )


		cumpi_j = 0.
		!draw on zgrid
		do ij=1,nj
			do iz=1,nz
				izp = 1
				cumpi_j(iz,izp+1) = piz(iz,izp,ij)
				do izp=2,nz
					cumpi_j(iz,izp+1) = piz(iz,izp,ij) + cumpi_j(iz,izp)
				enddo
			enddo
			! start everyone from good state, alternatively could start from random draw on ergodic dist
			z_jt_t = 3
			do it = 1,Tsim
				call random_number(z_innov)
				! use conditional probability
				z_jt_t = finder(cumpi_j(z_jt_t,:),z_innov ) 
				z_jt_macro(ij,it) = z_jt_t
			enddo
		enddo

		do i=1,Nsim
			do it=1,Tsim
			!fill in shock values from z_jt_macro for each indiv
				z_jt_int(i,it) = z_jt_macro(j_i(i) , it)
				z_jt(i,it) = zgrid(z_jt_int(i,it))
			enddo
		enddo

		success = 1
		
		deallocate(z_jt_macro)
		deallocate(cumpi_j)
	end subroutine draw_zjt
	
	subroutine draw_age_it(age_it, seed0, success)

		integer,intent(out) :: age_it(:,:)
		integer	:: it,itp,i,m,ss
		integer,intent(in) :: seed0
		integer,intent(out) :: success
		integer, dimension(100) :: bdayseed
		real(8), dimension(TT-1) :: cumpi_t0
		real(8), dimension(TT,TT) :: cumpi_t
		real(8), dimension(TT-1) :: prob_t_nTT
		real(8) :: rand_age
		
		call random_seed(size = ss)
		forall(m=1:ss) bdayseed(m) = (m-1)*100 + seed0
		call random_seed(put = bdayseed(1:ss) )

		!set up cumulative probabilities for t0 and conditional draws
		forall(it=1:TT,itp=1:TT) cumpi_t(it,itp) = 0.
		forall(it=1:TT-1) prob_t_nTT(it) = prob_t(it)/(1-prob_t(TT))
		it = 1
		cumpi_t0(it) = prob_t_nTT(it)
		cumpi_t(it,it) = 1.-ptau(it)
		do it=2,TT-1
			cumpi_t0(it) = prob_t_nTT(it) + cumpi_t0(it-1)
		enddo
		
		do i=1,Nsim
			it = 1
			call random_number(rand_age)

			age_it(i,it) = finder(cumpi_t0,rand_age)
			do it=2,Tsim
				call random_number(rand_age)
				if(rand_age < 1- ptau(age_it(i,it-1)) .and. age_it(i,it-1) < TT ) then
					age_it(i,it) = age_it(i,it-1)+1
				else 
					age_it(i,it) = age_it(i,it-1)
				endif
			enddo
		enddo
		success = 1
	end subroutine draw_age_it

	subroutine sim(val_funs, pol_funs,hists)
		
		implicit none
	
		type(val_struct), intent(inout), target :: val_funs
		type(pol_struct), intent(inout), target :: pol_funs	
		type(hist_struct), intent(inout), target :: hists


		integer :: i, ii, iter, it, j, idi, id, &
			&  seed0, seed1, status, m,ss
		integer :: bdayseed(100)
						
		real(8), allocatable ::	del_i(:) ! shocks to be drawn
		real(8), allocatable ::	al_it(:,:),z_jt(:,:)
		integer, allocatable :: z_jt_int(:,:) ! shocks to be drawn
		integer, allocatable :: del_i_int(:), j_i(:) ! integer valued shocks
		integer, allocatable :: al_it_int(:,:)! integer valued shocks
		integer, allocatable :: age_it(:,:) ! ages, drawn randomly

		integer, allocatable :: work_it(:,:), app_it(:,:) !choose work or not, apply or not
		real(8), allocatable :: work_dif_it(:,:), app_dif_it(:,:) !choose work or not, apply or not -- latent value
		integer, allocatable :: status_it(:,:)  !track W,U,N,D,R : 1,2,3,4,5
		real(8), allocatable :: e_it(:,:), a_it(:,:)
		integer, allocatable :: d_it(:,:), a_it_int(:,:),e_it_int(:,:)
		
		real(8), pointer ::	VR(:,:,:), &			!Retirement
					VD(:,:,:,:), &			!Disabled
					VN(:,:,:,:,:,:,:), &	!Long-term Unemployed
					VW(:,:,:,:,:,:,:), &	!Working
					VU(:,:,:,:,:,:,:), &	!Unemployed
					V(:,:,:,:,:,:,:)	!Participant
	
		real(8), pointer ::	gapp_dif(:,:,:,:,:,:,:), gwork_dif(:,:,:,:,:,:,:) ! latent value of work/apply
	
		integer, pointer ::	aR(:,:,:), aD(:,:,:,:), aU(:,:,:,:,:,:,:), &
					aN(:,:,:,:,:,:,:), aW(:,:,:,:,:,:,:)
		integer, pointer ::	gapp(:,:,:,:,:,:,:), &
					gwork(:,:,:,:,:,:,:)
		
		real(8) :: cumpid(nd,nd+1,ndi,TT-1),cumptau(TT+1)
	
		! Other
		real(8)	:: wage_hr,al_hr, junk , a_hr, e_hr, bet_hr,z_hr, work_dif_hr, app_dif_hr

		integer :: work_hr,app_hr
		integer :: ali_hr,d_hr,age_hr,del_hr, zi_hr, j_hr, ai_hr,api_hr,ei_hr, &
			& beti, status_hr,status_tmrw,drawi,drawt
		!************************************************************************************************!

		
		!************************************************************************************************!
		! Pointers
		!************************************************************************************************!
		! (disability extent, earn hist, assets)
		VR => val_funs%VR
		aR => pol_funs%aR
		VD => val_funs%VD
		aD => pol_funs%aD
		VN => val_funs%VN
		VU => val_funs%VU
		VW => val_funs%VW
		V => val_funs%V
		aN => pol_funs%aN
		aW => pol_funs%aW
		aU => pol_funs%aU
		gwork => pol_funs%gwork
		gapp => pol_funs%gapp

		gapp_dif => pol_funs%gapp_dif
		gwork_dif => pol_funs%gwork_dif

		allocate(z_jt(Nsim, Tsim))
		allocate(z_jt_int(Nsim, Tsim))		
		allocate(del_i(Nsim))
		allocate(al_it(Nsim,Tsim))
		allocate(j_i(Nsim))
		allocate(age_it(Nsim,Tsim))

		allocate(del_i_int(Nsim))
		allocate(al_it_int(Nsim,Tsim))


		allocate(a_it(Nsim,Tsim))
		allocate(a_it_int(Nsim,Tsim))		
		allocate(e_it(Nsim,Tsim))
		allocate(e_it_int(Nsim,Tsim))		
		allocate(d_it(Nsim,Tsim))
		allocate(work_it(Nsim,Tsim))
		allocate(work_dif_it(Nsim,Tsim))
		allocate(app_it(Nsim,Tsim))
		allocate(app_dif_it(Nsim,Tsim))
		allocate(status_it(Nsim,Tsim))

		seed0 = 941987
		seed1 = 12281951
		call random_seed(size = ss)
		forall(m=1:ss) bdayseed(m) = (m-1)*100 + seed0
		call random_seed(put = bdayseed(1:ss) )

		if(verbose >2) print *, "Drawing types and shocks"	
		call draw_age_it(age_it,seed0,status)
		call draw_deli(del_i,del_i_int, seed1, status)
		call draw_alit(al_it,al_it_int, seed0, status)
		call draw_ji(j_i,seed1, status)
		call draw_zjt(z_jt,z_jt_int, j_i, seed0, status)

		! check the distributions
		if(print_lev >=3 ) then
			call vec2csv(del_i,"del_i.csv")
			call veci2csv(j_i,"j_i.csv")
			call mat2csv(al_it,"al_it.csv")
			call mat2csv(z_jt,"z_jt.csv")
			call mati2csv(al_it_int,"al_it_int.csv")
			call mati2csv(z_jt_int,"z_jt_int.csv")
			call mati2csv(age_it,"age_it.csv")
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

		!use only 1 value of beta
		beti = 1
		bet_hr = 1.
		
		!itertate to get dist of asset/earnings correct at each age from which to draw start conditions 
		do iter=1,4
			! OMP  parallel do default(shared) &
			! OMP& private(iter, i, del_hr, j_hr, status_hr, it, drawi,drawt ai_hr,a_hr,api_hr,ei_hr,e_hr,wage_hr,junk,z_hr,zi_hr,age_hr,al_hr,ali_hr,d_hr) 
			do i=1,Nsim
				!fixed traits
				del_hr = del_i_int(i)
				j_hr = j_i(i)

				!initialize stuff
				it = 1
				status_hr = 1
				status_it(i,it) = 1

				!need to draw these from age-specific distributions for iterations > 1
				if(iter>1) then
					do ii=1,Nsim*Tsim
						call random_number(junk)
						drawi = max(1,nint(junk*Nsim))
						call random_number(junk)
						drawt = max(1,nint(junk*Tsim))
						if((age_it(drawi,drawt) .eq. age_it(i,it)) & 
							& .and. (status_it(drawi,drawt)  < 4)) then
							exit
						endif
					enddo
					d_it(i,it) = d_it(drawi,drawt)
					a_it(i,it) = a_it(drawi,drawt)
					e_it(i,it) = e_it(drawi,drawt)
					e_it_int(i,it) = e_it_int(drawi,drawt)
					a_it_int(i,it) = a_it_int(drawi,drawt)
					status_it(i,it) = status_it(drawi,drawt)
				endif

				
				do it=1,Tsim
					!set the state
					z_hr	= zgrid(z_jt_int(i,it))
					zi_hr	= z_jt_int(i,it)
					age_hr	= age_it(i,it)
					al_hr	= al_it(i,it)
					ali_hr	= al_it_int(i,it)
					d_hr	= d_it(i,it)
					e_hr 	= e_it(i,it)
					a_hr 	= a_it(i,it)
					ei_hr	= e_it_int(i,it)
					ai_hr 	= a_it_int(i,it)
					
					wage_hr	= wage(bet_hr,al_hr,d_hr,z_hr,age_hr)
					hists%wage_hist(i,it) = wage_hr
					
					status_hr = status_it(i,it)
					
					
					!random number will be used in several potential shocks
					call random_number(junk)
					
					!make decisions if not yet retired
					if(age_hr < TT) then 

						if(status_hr .eq. 3) then ! choose wait or apply
							app_hr = gapp_dif( (j_hr-1)*nbi + beti, (del_hr-1)*nai+ali_hr,d_hr,ei_hr,ai_hr,zi_hr,age_hr )
							
							if( app_dif_hr >= 0 ) then
							! choose to apply
								app_it(i,it) = 1
								if(junk<xi(d_hr)) status_tmrw = 4
							else
								app_it(i,it) = 0
							endif
							app_dif_it(i,it) = app_dif_hr
						endif

						! evalutate gwork and gapp to figure out lom of status 
						if((status_hr < 3) .or. (status_hr .eq. 4 .and. app_dif_hr < 0 ))then !choose work or rest
							work_dif_hr = gwork_dif( (j_hr-1)*nbi + beti, (del_hr-1)*nai+ali_hr,d_hr,ei_hr,ai_hr,zi_hr,age_hr )
							
							if( work_dif_hr > 0 ) then
							! choose to work
								if(status_hr < 3) then !not LTU
									status_it(i,it) = 1
									status_tmrw = 1
								else !LTU, have to get a good shock 
									if(junk <=rhho) status_tmrw = 1
									if(junk > rhho)  status_tmrw = status_hr
								endif
							elseif(status_hr .le. 2) then
								status_it(i,it) = 2
							! unemployed, may stay unemployed or become long-term unemployed
								if(junk <=pphi) status_tmrw = 3
								if(junk > pphi) status_tmrw = 2
							else
							!NDR no change, though may choose to apply for DI below
								status_tmrw =  status_hr
							endif
							work_dif_it(i,it) = work_dif_hr
						endif

						! the status that will go into next period
						if(it<Tsim) &
						&	status_it(i,it+1) = status_tmrw
						!evaluate the asset policy			
						if(status_hr .eq. 1) then
							api_hr = aw( (j_hr-1)*nbi + beti, (del_hr-1)*nai+ali_hr,d_hr,ei_hr,ai_hr,zi_hr,age_hr  )
				
						elseif(status_hr .eq. 2) then
							api_hr = aU( (j_hr-1)*nbi + beti, (del_hr-1)*nai+ali_hr,d_hr,ei_hr,ai_hr,zi_hr,age_hr )
						elseif(status_hr .eq. 3) then
							api_hr = aN( (j_hr-1)*nbi + beti, (del_hr-1)*nai+ali_hr,d_hr,ei_hr,ai_hr,zi_hr,age_hr )
						elseif(status_hr .eq. 4) then
							api_hr = aD( d_hr,ei_hr,ai_hr,age_hr )
						endif
					! retired
					else 
						api_hr = aR( d_hr,ei_hr,ai_hr )
						status_hr = 5
						if(it<Tsim) &
						&	status_it(i,it+1) = status_hr
					endif


					if(it<Tsim) then
						! push forward asset					
						a_it_int(i,it+1) = api_hr
						a_it(i,it+1) = agrid(api_hr)
						!push forward AIME
						if(status_hr .eq. 1) then
							!here, it is continuous
							e_it(i,it+1) = min( (e_hr*dble(it-1) + wage_hr)/dble(it),egrid(ne) )
							
							! assign to grid points by nearest neighbor
							ei_hr = finder(egrid,e_it(i,it+1))
							if(ei_hr < ne) then
							if((e_it(i,it+1) - egrid(ei_hr))<  (egrid(ei_hr+1) - e_it(i,it+1)) ) then
								e_it_int(i,it+1) = ei_hr
							else
								e_it_int(i,it+1) = ei_hr + 1
							endif
							endif
						else
							e_it(i,it+1) = e_it(i,it)
							e_it_int(i,it+1) = e_it_int(i,it)
						endif

						!push forward d 
						if(age_hr<TT .and. d_hr<nd ) then
							call random_number(junk)
							if( junk < pid(d_hr,d_hr+1,del_hr,age_hr) ) d_it(i,it+1) = d_it(i,it)+1 
						else 
							d_it(i,it+1) = d_it(i,it)
						endif
					endif
				enddo !1,Tsim
			enddo! 1,Nsim
			! OMP  end parallel do 

			
			if(print_lev >3)then
				call mat2csv (e_it,"e_it.csv")
				call mat2csv (a_it,"a_it.csv")
				call mati2csv(d_it,"d_it.csv")
			endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! check age - specific distributions of a_it, d_it for convergence
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		enddo! iter
		if(verbose >2 ) print *, "done simulating"
		!fill the histories
		hists%work_dif_hist = work_dif_it
		hists%app_dif_hist  = app_dif_it
		hists%obsX_hist = 0.
		do i=1,Nsim
			do it=1,Tsim
				do j=1,TT-1
					if(age_it(i,it) .eq. j ) hists%obsX_hist(i + (Nsim-1)*j,it) = 1
				enddo
			enddo
		enddo
		
		deallocate(d_it,a_it,e_it)
		deallocate(a_it_int,e_it_int)
		deallocate(z_jt,del_i,al_it,j_i)
		deallocate(app_dif_it,app_it,work_it,work_dif_it)
		deallocate(del_i_int,al_it_int,z_jt_int)

	end subroutine sim

end module sol_sim

!**************************************************************************************************************!
!**************************************************************************************************************!
!						MAIN PROGRAM						       !
!**************************************************************************************************************!
!**************************************************************************************************************!

program V0main

	use V0para
	use helper_funs
	use sol_sim

	implicit none


	!************************************************************************************************!
	! Counters and Indicies
	!************************************************************************************************!

		integer  :: id, it, ij, ibi, ial, iz, narg_in, wo


	!************************************************************************************************!
	! Other
	!************************************************************************************************!
		real(8)	:: wagehere
	!************************************************************************************************!
	! Structure to communicate everything
		type(val_struct) :: val_sol
		type(pol_struct) :: pol_sol
		type(hist_struct):: hists_sim
		type(moments_struct):: moments_sim
		
		
		
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
	allocate(val_sol%VN(nj*nbi,ndi*nai,nd,ne,na,nz,TT), stat=val_sol%alloced)
	allocate(val_sol%VU(nj*nbi,ndi*nai,nd,ne,na,nz,TT), stat=val_sol%alloced)
	allocate(val_sol%VW(nj*nbi,ndi*nai,nd,ne,na,nz,TT), stat=val_sol%alloced)
	allocate(val_sol%V(nj*nbi,ndi*nai,nd,ne,na,nz,TT), stat=val_sol%alloced)
	allocate(pol_sol%aN(nj*nbi,ndi*nai,nd,ne,na,nz,TT-1), stat=pol_sol%alloced)
	allocate(pol_sol%aW(nj*nbi,ndi*nai,nd,ne,na,nz,TT-1), stat=pol_sol%alloced)
	allocate(pol_sol%aU(nj*nbi,ndi*nai,nd,ne,na,nz,TT-1), stat=pol_sol%alloced)
	allocate(pol_sol%gwork(nj*nbi,ndi*nai,nd,ne,na,nz,TT-1), stat=pol_sol%alloced)
	allocate(pol_sol%gapp(nj*nbi,ndi*nai,nd,ne,na,nz,TT-1), stat=pol_sol%alloced)

	allocate(pol_sol%gapp_dif(nj*nbi,ndi*nai,nd,ne,na,nz,TT), stat=pol_sol%alloced)
	allocate(pol_sol%gwork_dif(nj*nbi,ndi*nai,nd,ne,na,nz,TT), stat=pol_sol%alloced)

	allocate(hists_sim%wage_hist(Nsim,Tsim), stat=hists_sim%alloced)
	allocate(hists_sim%work_dif_hist(Nsim,Tsim), stat=hists_sim%alloced)
	allocate(hists_sim%app_dif_hist(Nsim,Tsim), stat=hists_sim%alloced)
	allocate(hists_sim%obsX_hist(Nsim*Nk,Tsim), stat=hists_sim%alloced)

	narg_in = iargc()

	verbose = 4
	print_lev = 4


	call setparams()
	agrid(1) = .05*(agrid(1)+agrid(2))
	if(print_lev >= 3) then
		call vec2csv(agrid,"agrid.csv")
		wo = 0
		do ij = 1,nj
			call mat2csv(piz(:,:,ij),"piz.csv",wo)
			if(wo==0) wo =1
		enddo
		
		wo=0
		do it = 1,TT-1
			do ij = 1,ndi
				call mat2csv(pid(:,:,ij,it),"pid.csv",wo)
				if(wo==0) wo =1
			enddo
		enddo

		open(1, file="wage_dist.csv")
		ibi =1
		iz  =2
		do it = 1,TT-1
			do ial =1,nai
				do id = 1,nd-1
					wagehere = wage(beti(ibi),alfi(ial),id,zgrid(iz),it)
					write(1, "(G20.12)", advance='no') wagehere
				enddo
				id = nd
				wagehere = wage(beti(ibi),alfi(ial),id,zgrid(iz),it)
				write(1,*) wagehere
			enddo
			write(1,*) " "! trailing space
		enddo	
		close(1)
	endif

	if(verbose >2) print *, "Solving the model"
	call sol(val_sol,pol_sol)
	if(verbose >2) print *, "Simulating the model"	
	call sim(val_sol, pol_sol, hists_sim)

!    .----.   @   @
!   / .-"-.`.  \v/
!   | | '\ \ \_/ )
! ,-\ `-.' /.'  /
!'---`----'----'
	!****************************************************************************!
	! IF you love something.... 
	!****************************************************************************!
	deallocate(pol_sol%aR,pol_sol%aD,pol_sol%aN, pol_sol%aU, pol_sol%aW,pol_sol%gwork, pol_sol%gapp)
	deallocate(pol_sol%gwork_dif,pol_sol%gapp_dif)
	deallocate(val_sol%VR,val_sol%VD,val_sol%VN,val_sol%VU,val_sol%VW,val_sol%V)
	deallocate(hists_sim%wage_hist,hists_sim%work_dif_hist,hists_sim%app_dif_hist,hists_sim%obsX_hist)


End PROGRAM






