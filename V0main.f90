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
	!User-Defined Types (Structures) for value functions and policy functions
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
		real(8), allocatable ::	aR(:,:,:), aD(:,:,:,:), aU(:,:,:,:,:,:,:), &
					aN(:,:,:,:,:,:,:), aW(:,:,:,:,:,:,:)
		integer, allocatable ::	gapp(:,:,:,:,:,:,:), &

					gwork(:,:,:,:,:,:,:) !integer choice of apply/work
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

		integer  :: i, j, t, ia, ie, id, it, iaa,iaa1, iaa1app,iaa1napp,apol, ibi, iai, ij , idi, izz, iaai, idd, &
			    iee1, iee2, iz, unitno, print_lev, verbose, narg_in,iw,wo
		integer, dimension(5) :: maxer_i
	
		!************************************************************************************************!
		! Value Functions- Stack z-risk j and indiv. exposure beta_i
		!************************************************************************************************!
		real(8)  	  	:: Vtest1, Vtest2, utilhere, Vapp, VC, app2, Vc1, Vnapp, anapp,aapp, maxer_v, smthV,smthV0param, &
					&	iee1wt
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
	
		real(8), pointer ::	aR(:,:,:), aD(:,:,:,:), aU(:,:,:,:,:,:,:), &
					aN(:,:,:,:,:,:,:), aW(:,:,:,:,:,:,:)
		integer, pointer ::	gapp(:,:,:,:,:,:,:), &
					gwork(:,:,:,:,:,:,:)
	
		!************************************************************************************************!
		! Other
		!************************************************************************************************!
			real(8)	:: wagehere,chere, junk,summer, eprime, yL, yH,VUhere, VWhere
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
		V => val_funs%V
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

		!************************************************************************************************!
		! Caculate things that are independent of occupation/person type
		!	1) Value of Retired:  VR(d,e,a)
		!	2) Value of Disabled: VD(d,e,a)
		!************************************************************************************************!

	!1) Calculate Value of Retired: VR(d,e,a)
		!d in{1,2,3}  : disability extent
		!e inR+       :	earnings index
		!a inR+	      : asset holdings
		
		iw = 1 ! not working
		!VFI with good guess
		!Initialize
		junk = (1.0+(beta*ptau(TT)*R**(1.0-gam))**(-1.0/gam))**(-1.0)
		do ie=1,ne
		do ia=1,na
		do id=1,nd
			VR0(id,ie,ia) = util(SSI(egrid(ie))+R*agrid(ia),id,iw)* (1./(1.-beta*ptau(TT)))
		enddo
		enddo
		enddo
		if(print_lev >3) then
			call vec2csv(VR0(1,1,:),"VR0.csv",0)
		endif		
		j = 1
		do WHILE (j<=maxiter)
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

							if(Vtest2>Vtest1) then
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
					VR(1,ie,ia) = Vtest1
					aR(1,ie,ia) = agrid(apol)
					summer = summer+ (VR(1,ie,ia)-VR0(1,ie,ia))**2
					!Polices for disabled are the same, just scale V-function

				enddo !ia
			enddo !ie
			IF (summer < Vtol) THEN
				do ie =1,ne
				do id =2,ne
					VR(id,ie,:) = VR(1,ie,:)*(dexp(theta*dble(id-1)))**(1-gam)
					aR(id,ie,:) = aR(1,ie,:)
				enddo
				enddo
				exit	!Converged				
			else
				VR0 = VR	!New guess
				j=j+1
				if(print_lev >3) then
					call vec2csv(aR(1,1,:),"aR.csv",0)
					call vec2csv(VR(1,1,:),"VR.csv",0)
				endif
			endif
		enddo ! iteration j
		if (print_lev > 2) then
			call vec2csv(aR(1,1,:),"aR.csv",0)
			call vec2csv(VR(1,1,:),"VR.csv",0)
		endif

		!----------------------------------------------------------!
		!Set value at t=TT to be VR in all other V-functions
		!----------------------------------------------------------!

			
		do ij=1,nj
		do ibi=1,nbi
		do idi=1,ndi
		do iai=1,nai
		do iz=1,nz   
			do id=1,nd
			do ie=1,ne
			do ia=1,na
				VD(id,ie,ia,TT) = VR(id,ie,ia)
				VD0(id,ie,ia,TT) = VR(id,ie,ia)				
				VW ((ij-1)*nbi+ibi,(idi-1)*nai+iai,id,ie,ia,iz,TT) = VR(id,ie,ia)
				VW0((ij-1)*nbi+ibi,(idi-1)*nai+iai,id,ie,ia,iz,TT) = VR(id,ie,ia)
				VN ((ij-1)*nbi+ibi,(idi-1)*nai+iai,id,ie,ia,iz,TT) = VR(id,ie,ia)
				VN0((ij-1)*nbi+ibi,(idi-1)*nai+iai,id,ie,ia,iz,TT) = VR(id,ie,ia)
				VU ((ij-1)*nbi+ibi,(idi-1)*nai+iai,id,ie,ia,iz,TT) = VR(id,ie,ia)
				VU0((ij-1)*nbi+ibi,(idi-1)*nai+iai,id,ie,ia,iz,TT) = VR(id,ie,ia)
				V  ((ij-1)*nbi+ibi,(idi-1)*nai+iai,id,ie,ia,iz,TT) = VR(id,ie,ia)
				V0 ((ij-1)*nbi+ibi,(idi-1)*nai+iai,id,ie,ia,iz,TT) = VR(id,ie,ia)	   
			enddo
			enddo
			enddo
		enddo
		enddo
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
		do it = TT-1,1,-1

			id = 1 ! other values are just multiples thereof
			iw = 1 ! not working
			
			!Guess will be value at t+1
			id = 1
			do ia=1,na
			do ie=1,ne
				VD0(id,ie,ia,TT-it) = VD(id,ie,ia,TT-it+1)
			enddo
			enddo
			if(print_lev >3) then
				call vec2csv(VD0(1,ie,:,TT-it),"VD0.csv",0)
			endif

			!Loop to find V(..,it) as fixed point
			j=1
			do WHILE (j<maxiter)
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
						VD(1,ie,ia,TT-it) = Vtest1
						aD(1,ie,ia,TT-it) = agrid(apol)
						summer = summer+ (VD(1,ie,ia,TT-it)-VD0(1,ie,ia,TT-it))**2
					enddo	!ia
				enddo	!ie		
				if(print_lev >3) then
					wo =0
					call vec2csv(aD(1,1,:,TT-it),"aD.csv",wo)
					call vec2csv(VD(1,1,:,TT-it),"VD.csv",wo)
				endif
				IF (summer < Vtol) THEN
					exit	!Converged
				EndIF
				VD0 = VD	!New guess
				j=j+1
			enddo	!j: V-iter loop
			! set the other value function/asset policies for other disability levels
			do id = 2,nd
				VD(id,ie,ia,TT-it) = VD(id,ie,ia,TT-it)*((dexp(theta*dble(id-1)))**(1-gam))
				aD(id,ie,ia,TT-it) = aD(id,ie,ia,TT-it)
			enddo
		enddo	!t loop, going backwards
		
		if (print_lev > 2) then
			wo =0 
			call mat2csv(aD(1,1,:,:),"aD.csv",wo)
			call mat2csv(VD(1,1,:,:),"VD.csv",wo)
		endif


		smthV0param  =10. ! will tighten this down
	!************************************************************************************************!
	!3) Calculate V= max(VW,VN); requires calculating VW and VN
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
			do iai=1,nai
			do id =1,nd
			do ie =1,ne
			do iz =1,nz
			do ia =1,na

			!Guess once, then use prior occupation/beta as guess
				!IF (j .EQ. 1 .AND. ibi .EQ. 1 .AND. idi .EQ. 1) THEN
				!IF (j .EQ. 1 .AND. ibi .EQ. 1 .AND. idi .EQ. 1) THEN
				 !0) Guess VW0(nj,nbi,nai,nd,ne,na,nz,TT-1)
					VW0((ij-1)*nbi+ibi,(idi-1)*nai+iai,id,ie,ia,iz,TT-it) = VW((ij-1)*nbi+ibi,(idi-1)*nai+iai,id,ie,ia,iz,TT-it+1)
				 !0) Guess VU0(nj,nbi,nai,nd,ne,na,nz,TT-1)
					VU0((ij-1)*nbi+ibi,(idi-1)*nai+iai,id,ie,ia,iz,TT-it) = VU((ij-1)*nbi+ibi,(idi-1)*nai+iai,id,ie,ia,iz,TT-it+1)
				 !0) Guess VN0(nj,nbi,nai,nd,ne,na,nz,TT-1)
					VN0((ij-1)*nbi+ibi,(idi-1)*nai +iai,id,ie,ia,iz,TT-it) = VN((ij-1)*nbi+ibi,(idi-1)*nai+iai,id,ie,ia,iz,TT-it+1)

				IF (it .EQ. 1) THEN
					VU0((ij-1)*nbi+ibi,(idi-1)*nai+iai,1,ie,ia,iz,TT-it) = VW0((ij-1)*nbi+ibi,(idi-1)*nai+iai,1,ie,ia,iz,TT-it)
				EndIF
				!ELSE
				! !0) Guess VW0(nj,nbi,nai,nd,ne,na,nz,TT-1)
				!	VW0((ij-1)*nbi+ibi,idi,iai,id,ie,ia,iz,TT-it) = VW(max(1,ij-1)*nbi+max(1,ibi-1),max(1,idi-1),iai,id,ie,ia,iz,TT-it+1)
				! !0) Guess VN0(nj,nbi,nai,nd,ne,na,nz,TT-1)
				!	VN0((ij-1)*nbi+ibi,idi,iai,id,ie,ia,iz,TT-it) = VN(max(1,ij-1)*nbi+max(1,ibi-1),max(1,idi-1),iai,id,ie,ia,iz,TT-it+1)
				!EndIF !first occupationXbeta loop
				!0) Calculate V0(nj,nbi,nai,nd,ne,na,nz,it) by max:	
				!V0((ij-1)*nbi+ibi,(idi-1)*nai+iai,id,ie,ia,iz,TT-it)= max(VW0((ij-1)*nbi+ibi,(idi-1)*nai+iai,id,ie,ia,iz,TT-it),VU0((ij-1)*nbi+ibi,(idi-1)*nai+iai,id,ie,ia,iz,TT-it))
				!0) Calculate V0(nj,nbi,nai,nd,ne,na,nz,it) so its smooth
				!smthV = dexp( smthV0param * VW0((ij-1)*nbi+ibi,(idi-1)*nai+iai,id,ie,ia,iz,TT-it) ) &
				!	& /(dexp(smthV0param * VW0((ij-1)*nbi+ibi,(idi-1)*nai+iai,id,ie,ia,iz,TT-it)) + dexp(smthV0param * VN0((ij-1)*nbi+ibi,(idi-1)*nai+iai,id,ie,ia,iz,TT-it)) )
				!V0((ij-1)*nbi+ibi,(idi-1)*nai+iai,id,ie,ia,iz,TT-it)= smthV*VW0((ij-1)*nbi+ibi,(idi-1)*nai+iai,id,ie,ia,iz,TT-it) + (1.-smthV)*VN0((ij-1)*nbi+ibi,(idi-1)*nai+iai,id,ie,ia,iz,TT-it)
				V0((ij-1)*nbi+ibi,(idi-1)*nai+iai,id,ie,ia,iz,TT-it)= VW0((ij-1)*nbi+ibi,(idi-1)*nai+iai,id,ie,ia,iz,TT-it)


			enddo	!ia
			enddo	!iz
			enddo	!ie
			enddo 	!id
			enddo	!iai   

			!***********************************************************************************************************
			!Loop over V=max(VU,VW)	
			j=1
			do WHILE (j<maxiter)
				maxer = 0.
				summer = 0.	!Use to calc |V-V0|<eps

				! lots f printing every 100 iterations (mod 1 because Fortran is designed by idiots using base-1 indexing)
				if(mod(j,100) .eq. 1) then
					print_lev = 4
				else 
					print_lev =1
				endif
			!------------------------------------------------!
			!Solve VU given guesses on VW, VN, VU and implied V
			!------------------------------------------------!  
			!$OMP  parallel do default(shared) &
			!$OMP& private(iai,id,ie,iz,iw,apol,iaa1,ia,iaa,chere,Vtest2,Vtest1,Vc1,iaai,izz)
			  	do iai=1,nai	!Loop over alpha (ai)
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
								!Continuation value if don't go on disability
								Vtest2 = 0.
								do izz = 1,nz	 !Loop over z'
								do iaai = 1,nai !Loop over alpha_i'
									Vc1 = (1.-ptau(TT-it))*(pphi*VN0((ij-1)*nbi+ibi,(idi-1)*nai+iaai,id,ie,iaa,izz,TT-it+1)+(1-pphi)*V0((ij-1)*nbi+ibi,(idi-1)*nai+iaai,id,ie,iaa,izz,TT-it+1)) & !Age and might go LTU
										& +ptau(TT-it)*(pphi*VN0((ij-1)*nbi+ibi,(idi-1)*nai+iaai,id,ie,iaa,izz,TT-it)+(1-pphi)*V0((ij-1)*nbi+ibi,(idi-1)*nai+iaai,id,ie,iaa,izz,TT-it))     !Don't age, maybe LTU
									Vtest2 = Vtest2 + beta*piz(iz,izz,ij)*pialf(iai,iaai)*(Vc1)  !Probability of alpha_i X z_i draw 
								enddo
								enddo
								Vtest2 = Vtest2 + util(chere,id,iw)
								if(Vtest2>Vtest1) then
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
						VU((ij-1)*nbi+ibi,(idi-1)*nai+iai,id,ie,ia,iz,TT-it) = Vtest1
						aU((ij-1)*nbi+ibi,(idi-1)*nai+iai,id,ie,ia,iz,TT-it) = agrid(apol)

					enddo !ia
				
				enddo !id
			  	enddo !ie
			  	enddo !iz
			  	enddo !iai	
			!$OMP END PARALLEL do

				if (print_lev > 3) then
					wo = 0
					do iai=1,nai	!Loop over alpha (ai)
					do id=1,nd	!Loop over disability index
				  	do ie=1,ne	!Loop over earnings index
				  	do iz=1,nz	!Loop over TFP
						call vec2csv(aU((ij-1)*nbi+ibi,(idi-1)*nai+iai,id,ie,:,iz,TT-it),"aU.csv",wo)
						call vec2csv(VU((ij-1)*nbi+ibi,(idi-1)*nai+iai,id,ie,:,iz,TT-it),"VU.csv",wo)
						if(wo == 0) wo = 1 
					enddo
					enddo
					enddo
					enddo
				endif

	
			!------------------------------------------------!
			!Solve VN given guesses on VW, VN, and implied V
			!------------------------------------------------! 

			!$OMP  parallel do default(shared) &
			!$OMP& private(iai,id,ie,iz,iw,apol,iaa1app,iaa1napp,ia,iaa,chere,Vc1,Vtest2,Vtest1,smthV,Vapp,Vnapp,aapp,anapp,iaai,izz) 
			  	do iai=1,nai	!Loop over alpha (ai)
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
										& +xi(id)*VD(id,ie,iaa,TT-it+1)) !Age and might go on DI
									Vc1 = Vc1 + ptau(TT-it)*((1-xi(id))*VN0((ij-1)*nbi+ibi,(idi-1)*nai+iaai,id,ie,iaa,izz,TT-it) &
										& +xi(id)*VD(id,ie,iaa,TT-it))     !Don't age, might go on DI		
									Vtest2 = Vtest2 + beta*piz(iz,izz,ij)*pialf(iai,iaai)*Vc1 
								enddo
								enddo
								Vtest2 = util(chere,id,iw) - nu + Vtest2

								if (Vtest2>Vtest1) then	
									apol = iaa
									Vtest1 = Vtest2
								elseif(iaa > apol +iaa_hiwindow) then
									exit
								endif
							elseif(iaa<= iaa1app .and. Vtest1 <= -1e4 .and. apol <= iaa1app) then
								iaa = 1 !started too much saving, go back towards zero
								iaa1app = 1
							else
								exit
							endif
							iaa = iaa + 1
						enddo !iaa
						Vapp = Vtest1
						iaa1app = max(apol -iaa_lowindow,1) !concave, start next loop here
						aapp = agrid(apol)					
						if(Vapp <-1e4) then
							print *, "ruh roh!"
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
									Vc1 = (1-ptau(TT-it))*((1-rhho)*VN0((ij-1)*nbi+ibi,(idi-1)*nai +iaai,id,ie,iaa,izz,TT-it+1) &
										& +rhho*V0((ij-1)*nbi+ibi,(idi-1)*nai+iaai,id,ie,iaa,izz,TT-it+1)) !Age and might go on DI
									Vc1 = Vc1+ptau(TT-it)*((1-rhho)*VN0((ij-1)*nbi+ibi,(idi-1)*nai +iaai,id,ie,iaa,izz,TT-it) &
										& +rhho*V0((ij-1)*nbi+ibi,(idi-1)*nai+iaai,id,ie,iaa,izz,TT-it))     !Don't age, might go on DI
									Vtest2 = Vtest2 + beta*piz(iz,izz,ij)*pialf(iai,iaai)*(Vc1) 
								enddo
								enddo
								Vtest2 = Vtest2 + util(chere,id,iw)
								
								if( Vtest2 > Vtest1) then
									apol = iaa 
									Vtest1 = Vtest2 
								elseif(iaa > apol+iaa_hiwindow) then
									exit
								endif
							elseif(iaa == iaa1napp .and. Vtest1 <= -1e4 .and. apol == iaa1napp) then
								iaa = 1 !started too much saving, go back towards zero
								iaa1napp = 1
							else
								exit
							endif
							iaa = iaa + 1
						enddo !iaa
						iaa1napp = max(apol-iaa_lowindow,1)
						Vnapp = Vtest1 					
						aNapp = agrid(apol)


						!******************************************************
						!***************** Discrete choice for application
						smthV = dexp(smthV0param*Vnapp)/( dexp(smthV0param*Vnapp) +dexp(smthV0param*Vapp) )
						if( smthV .lt. 1e-5 .or. smthV .gt. 0.999999 .or. isnan(smthV)) then
							if( Vapp  > Vnapp ) smthV =0.
							if(Vnapp > Vapp  ) smthV =1.
						endif
						IF (Vapp > Vnapp) THEN
							!VN((ij-1)*nbi+ibi,(idi-1)*nai+iai,id,ie,ia,iz,TT-it) = Vapp
							aN((ij-1)*nbi+ibi,(idi-1)*nai+iai,id,ie,ia,iz,TT-it) = aapp
							gapp((ij-1)*nbi+ibi,(idi-1)*nai+iai,id,ie,ia,iz,TT-it) = 1
						ELSE !Don't apply
							!VN((ij-1)*nbi+ibi,(idi-1)*nai+iai,id,ie,ia,iz,TT-it) = Vnapp
							aN((ij-1)*nbi+ibi,(idi-1)*nai+iai,id,ie,ia,iz,TT-it) = aNapp
							gapp((ij-1)*nbi+ibi,(idi-1)*nai+iai,id,ie,ia,iz,TT-it) = 0
						EndIF

						gapp_dif((ij-1)*nbi+ibi,(idi-1)*nai+iai,id,ie,ia,iz,TT-it) = Vapp - Vnapp
						if(verbose .gt. 4) print *, Vapp - Vnapp
						
						
						VN((ij-1)*nbi+ibi,(idi-1)*nai+iai,id,ie,ia,iz,TT-it) = smthV*Vnapp + (1. - smthV)*Vapp
					enddo !ia
				enddo !iz 
				enddo !ie 
				enddo !id 
				enddo !iai
				
			!$OMP  END PARALLEL do
			
				!------------------------------------------------!			
				! Done making VN
				  	
			  	if (print_lev >3) then
			  		wo = 0
				  	do iai=1,nai	!Loop over alpha (ai)
					do id=1,nd	!Loop over disability index
				  	do ie=1,ne	!Loop over earnings index
				  	do iz=1,nz	!Loop over TFP
				  		call mat2csv(VN((ij-1)*nbi+ibi,(idi-1)*nai+iai,id,ie,:,iz,:) ,"VN.csv",wo)
				  		call mat2csv(aN((ij-1)*nbi+ibi,(idi-1)*nai+iai,id,ie,:,iz,:) ,"aN.csv",wo)
				  		call mat2csv(VN((ij-1)*nbi+ibi,(idi-1)*nai+iai,id,ie,:,iz,:) ,"VU.csv",wo)
				  		call mat2csv(aN((ij-1)*nbi+ibi,(idi-1)*nai+iai,id,ie,:,iz,:) ,"aU.csv",wo)
				  		call mati2csv(gapp((ij-1)*nbi+ibi,(idi-1)*nai+iai,id,ie,:,iz,:) ,"gapp.csv",wo)
				  		call mat2csv(gapp_dif((ij-1)*nbi+ibi,(idi-1)*nai+iai,id,ie,:,iz,:) ,"gapp_dif.csv",wo)

						if(wo == 0 ) wo =1		  		
					enddo !iz 
					enddo !ie 
					enddo !id 
					enddo !iai 
					
					wo = 0
					do iz=1,nz
					do ie=1,ne
					do ia=1,na
						call mat2csv(gapp_dif((ij-1)*nbi+ibi,:,:,ie,ia,iz,1),'dilat0_dalpha.csv',wo)
						if(wo == 0 ) wo =1
						call mat2csv(gapp_dif((ij-1)*nbi+ibi,:,:,ie,ia,iz,2),'dilat0_dalpha.csv',wo)
					enddo
					enddo
					enddo
					
						
			  	endif
		  	
			  	!update VN0
			  	do iai=1,nai	!Loop over alpha (ai)
				do id=1,nd	!Loop over disability index
			  	do ie=1,ne	!Loop over earnings index
			  	do iz=1,nz	!Loop over TFP
					do ia =1,na
						VN0((ij-1)*nbi+ibi,(idi-1)*nai+iai,id,ie,ia,iz,TT-it) = VN((ij-1)*nbi+ibi,(idi-1)*nai+iai,id,ie,ia,iz,TT-it)
					enddo !ia		  	
				enddo !iai
			  	enddo !id
			  	enddo !ie
			  	enddo !iz

				!------------------------------------------------!
				!Solve VW given guesses on VW, VN, and implied V
				!------------------------------------------------!
			!$OMP   parallel do default(shared) reduction(+:summer) &
			!$OMP & private(iai,id,ie,iz,apol,eprime,wagehere,iee1,iee2,egrid,iee1wt,ia,iaa,iaa1,chere,yL,yH,Vc1,utilhere,Vtest2,Vtest1,smthV,VUhere,VWhere,iaai,izz) 
			  	do iai=1,nai	!Loop over alpha (ai)
				do id=1,nd	!Loop over disability index
			  	do ie=1,ne	!Loop over earnings index
			  	do iz=1,nz	!Loop over TFP
					!Earnings evolution independent of choices.
					wagehere = wage(beti(ibi),alfi(iai),id,zgrid(iz),TT-it)
					eprime = Hearn(TT-it,ie,wagehere)
					!linear interpolate for the portion that blocks off bounds on assets
					if(eprime > minval(egrid) .and. eprime < maxval(egrid)) then
						iee1 = finder(egrid,eprime) !apparently this causes a race unless egrid is private
						iee2 = min(ne,iee1+1)
						iee1wt = (egrid(iee2)-eprime)/(egrid(iee2)-egrid(iee1))
					elseif( eprime <= minval(egrid)  ) then 
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
									Vc1 = piz(iz,izz,ij)*pialf(iai,iaai) &
										& * (yH*(1. - iee1wt) + yL*iee1wt )&
										& + Vc1
								enddo
								enddo
								utilhere = util(chere,id,2)
								Vtest2 = utilhere + beta*Vc1 ! flow utility

								IF (Vtest2>Vtest1 ) THEN				     
									Vtest1 = Vtest2
									apol = iaa
								elseif(iaa>apol+iaa_hiwindow) then
									exit
								endif
							elseif(iaa == iaa1 .and. Vtest1<=-5e4 .and. apol == iaa1 ) then
								!back it up
								iaa = 1
								iaa1 =1
							else 
								exit
							endif
							iaa = iaa+1
						enddo	!iaa
						iaa1 = max(iaa -iaa_lowindow,1)	!concave? start next loop here
						VW((ij-1)*nbi+ibi,(idi-1)*nai+iai,id,ie,ia,iz,TT-it) = Vtest1
						aW((ij-1)*nbi+ibi,(idi-1)*nai+iai,id,ie,ia,iz,TT-it) = agrid(apol)


						!------------------------------------------------!
						!Calculate V with solved vals of VW and VU -- i.e. can quit into unemployment
						!------------------------------------------------!
						VUhere = VU((ij-1)*nbi+ibi,(idi-1)*nai+iai,id,ie,ia,iz,TT-it)
						VWhere = VW((ij-1)*nbi+ibi,(idi-1)*nai+iai,id,ie,ia,iz,TT-it)

						IF (VWhere>VUhere) THEN
							gwork((ij-1)*nbi+ibi,(idi-1)*nai+iai,id,ie,ia,iz,TT-it) = 1
						ELSE
							gwork((ij-1)*nbi+ibi,(idi-1)*nai+iai,id,ie,ia,iz,TT-it) = 0
						EndIF
						smthV = dexp( smthV0param *VWhere  ) &
							& /( dexp(smthV0param * VWhere) + dexp(smthV0param * VUhere ) )
						if (smthV <1e-5 .or. smthV>.999999 .or. isnan(smthV) ) then
							if(VWhere .gt. VUhere)	smthV = 1. 
							if(VWhere .lt. VUhere)	smthV = 0. 							
						endif
						V((ij-1)*nbi+ibi,(idi-1)*nai+iai,id,ie,ia,iz,TT-it) = smthV*VWhere &
								& + (1.-smthV)*VUhere

						gwork_dif((ij-1)*nbi+ibi,(idi-1)*nai+iai,id,ie,ia,iz,TT-it) = VWhere - VUhere
						
						summer = summer+ & 
							& (V((ij-1)*nbi+ibi,(idi-1)*nai+iai,id,ie,ia,iz,TT-it)-V0((ij-1)*nbi+ibi,(idi-1)*nai+iai,id,ie,ia,iz,TT-it))**2
						
						maxer(ia,iz,ie,id,iai) = (V((ij-1)*nbi+ibi,(idi-1)*nai+iai,id,ie,ia,iz,TT-it)-V0((ij-1)*nbi+ibi,(idi-1)*nai+iai,id,ie,ia,iz,TT-it))**2

					enddo !ia
				enddo !iz
			  	enddo !ie
			  	enddo !id
			  	enddo !iai

			!$OMP  END PARALLEL do
			
				maxer_v = maxval(maxer)
				maxer_i = maxloc(maxer)

			  	!update VW0
			  	do iai=1,nai	!Loop over alpha (ai)
				do id=1,nd	!Loop over disability index
			  	do ie=1,ne	!Loop over earnings index
			  	do iz=1,nz	!Loop over TFP
				do ia =1,na
					VW0((ij-1)*nbi+ibi,(idi-1)*nai+iai,id,ie,ia,iz,TT-it) = VW((ij-1)*nbi+ibi,(idi-1)*nai+iai,id,ie,ia,iz,TT-it)
				enddo !ia		  	
				enddo !iai
			  	enddo !id
			  	enddo !ie
			  	enddo !iz

				if (print_lev >3) then
					wo = 0
					do iai=1,nai	!Loop over alpha (ai)
					do id=1,nd	!Loop over disability index
					do ie=1,ne	!Loop over earnings index
					do iz=1,nz	!Loop over TFP
						call mat2csv(VW((ij-1)*nbi+ibi,(idi-1)*nai+iai,id,ie,:,iz,:) ,"VW.csv",wo)
						call mat2csv(aW((ij-1)*nbi+ibi,(idi-1)*nai+iai,id,ie,:,iz,:) ,"aW.csv",wo)
						call mati2csv(gwork((ij-1)*nbi+ibi,(idi-1)*nai+iai,id,ie,:,iz,:) ,"gwork.csv",wo)
						if(wo==0) wo =1
					enddo !iz 
					enddo !ie 
					enddo !id 
					enddo !iai 	
				endif


			  	!update V0
			  	do iai=1,nai	!Loop over alpha (ai)
				do id=1,nd	!Loop over disability index
			  	do ie=1,ne	!Loop over earnings index
			  	do iz=1,nz	!Loop over TFP
				do ia =1,na
					V0((ij-1)*nbi+ibi,(idi-1)*nai+iai,id,ie,ia,iz,TT-it) = V((ij-1)*nbi+ibi,(idi-1)*nai+iai,id,ie,ia,iz,TT-it)
				enddo !ia		  	
				enddo !iai
			  	enddo !id
			  	enddo !ie
			  	enddo !iz


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
				! End of j iteration
				!------------------------------------------------!
				!Check |V-V0|<eps
				!------------------------------------------------!
				WRITE(*,*) summer, j, ij, ibi, idi, it
				WRITE(*,*) maxer_v, maxer_i(1), maxer_i(2), maxer_i(3), maxer_i(4), maxer_i(5)
				
				IF (summer < Vtol) THEN
					exit !Converged
				EndIF
				V0 = V	!New guess
				j=j+1
				smthV0param = smthV0param*1.5 !tighten up the discrete choice
			enddo	!j: V-iter loop
	!WRITE(*,*) ij, ibi, idi, it
		enddo	!t loop, going backwards

	enddo	!idi
	enddo	!ibi
	enddo	!ij	
		! reset V0
	!	do iai=1,nai	!Loop over alpha (ai)
	!	do id=1,nd	!Loop over disability index
	! 	do ie=1,ne	!Loop over earnings index
	!	do iz=1,nz	!Loop over TFP
	!  	do ia=1,na
	!		V0((ij-1)*nbi+ibi,(idi-1)*nai+iai,id,ie,ia,iz,TT-it) = V((ij-1)*nbi+ibi,(idi-1)*nai+iai,id,ie,ia,iz,TT-it)
	!	enddo !iai
	!	enddo !id
	!	enddo !ie
	!	enddo !iz
	! 	enddo !ia





	!WRITE(*,*) aD(1,:,:,TT-5)
	!WRITE(*,*) '-----'
	!WRITE(*,*) VR
	!WRITE(*,*) '-----'
	!WRITE(*,*) summer, j

	ie = INT(na/4)
	j = INT(nj/2)
	idi = INT(ndi/2)
	ia = INT(ne/3)

	! this plots work-rest and di application on the cross product of alphai and deltai and di
	wo  = 0
	do iz=1,nz
	do ie=1,ne
	do ia=1,na
		do it=1,TT-1
	!                 nj*nbi,ndi*nai,nd,ne,na,nz,TT-1
			call mati2csv(gapp(ibi,:,:,ie,ia,iz,it),'dipol_dalpha.csv',wo)
			call mati2csv(gwork(ibi,:,:,ie,ia,iz,it),'workpol_dalpha.csv',wo)

			call mat2csv(gapp_dif(ibi,:,:,ie,ia,iz,it),'dilat_dalpha.csv',wo)
			call mat2csv(gwork_dif(ibi,:,:,ie,ia,iz,it),'worklat_dalpha.csv',wo)
			if(wo ==0 ) wo =1
		enddo
	enddo
	enddo
	enddo

	call vec2csv(V(1,1,:,ie,ia,2,2),'Vfun.csv',0)
	call vec2csv(dtype,'DriskGrid.csv',0)
	call vec2csv(alfi(:),'AlfGrid.csv',0)
	call vec2csv(occz(:),'ZriskGrid.csv',0)
	call vec2csv(agrid(:),'Agrid.csv',0)

!	val_funs%VR = VR
!	pol_funs%aR = aR
!	val_funs%VD = VD
!	pol_funs%aD = aD
!	val_funs%VN = VN
!	val_funs%VU = VU
!	val_funs%VW = VW
!	val_funs%V = V 
!	pol_funs%aN = aN
!	pol_funs%aW = aW
!	pol_funs%aU = aU
!	pol_funs%gwork = gwork
!	pol_funs%gapp = gapp

!	pol_funs%gapp_dif = gapp_dif
!	pol_funs%gwork_dif = gwork_dif



	deallocate(maxer)
	deallocate(VR0,VD0,VN0,VU0,VW0,V0)
!	deallocate(VR,VD,VN,VU,VW,V)
!	deallocate(aR,aD,aN,aW,aU,gwork,gapp,gapp_dif,gwork_dif)
	end subroutine sol 


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

		integer  :: i, j, t, ia, ie, id, it, iaa,iaa1, iaa1app,iaa1napp,apol, ibi, iai, ij , idi, izz, iaai, idd, &
			    iee1, iee2, iz, unitno, print_lev, verbose, narg_in,iw,wo
		integer, dimension(5) :: maxer_i
	
	!************************************************************************************************!
	! Value Functions- Stack z-risk j and indiv. exposure beta_i
	!************************************************************************************************!
	real(8)  	  	:: Vtest1, Vtest2, utilhere, Vapp, VC, app2, Vc1, Vnapp, anapp,aapp, maxer_v, smthV,smthV0param, &
				&	iee1wt
        real(8), allocatable	:: maxer(:,:,:,:,:)
	real(8), allocatable	:: VR(:,:,:), &			!Retirement
				VD(:,:,:,:), &			!Disabled
				VN(:,:,:,:,:,:,:), &	!Long-term Unemployed
				VW(:,:,:,:,:,:,:), &	!Working
				VU(:,:,:,:,:,:,:), &	!Unemployed
				V(:,:,:,:,:,:,:)	!Participant
	
	real(8), allocatable	:: gapp_dif(:,:,:,:,:,:,:), gwork_dif(:,:,:,:,:,:,:) ! latent value of work/apply
	

	!************************************************************************************************!
	! Policies objects- Stack z-risk j and indiv. exposure beta_i
	!	+Asset grids correspond to V-funs
	!	+g-grids correspond to discrete choices
	!************************************************************************************************!
	real(8), allocatable	:: aR(:,:,:), aD(:,:,:,:), aU(:,:,:,:,:,:,:), &
				aN(:,:,:,:,:,:,:), aW(:,:,:,:,:,:,:)
	integer, allocatable  	:: aiD(:,:,:,:), gapp(:,:,:,:,:,:,:), &
				gwork(:,:,:,:,:,:,:)
	
	!************************************************************************************************!
	! Other
	!************************************************************************************************!
		real(8)	:: wagehere,chere, junk,summer, eprime, yL, yH,VUhere, VWhere
	!************************************************************************************************!
	! Structure to communicate everything
		type(val_struct) :: val_sol
		type(pol_struct) :: pol_sol
		
		
	!************************************************************************************************!
	! Allocate phat matrices
	!************************************************************************************************!
	! (disability extent, earn hist, assets)
	allocate(VR(nd,ne,na), stat=val_sol%alloced)
	allocate(aR(nd,ne,na), stat=pol_sol%alloced)
	! (disability extent, earn hist, assets, age)
	allocate(VD(nd,ne,na,TT), stat=val_sol%alloced)
	allocate(aD(nd,ne,na,TT-1), stat=pol_sol%alloced)
	allocate(aiD(nd,ne,na,TT-1))	
	! (occupation X ind exposure, ind disb. risk X ind. wage, disab. extent, earn hist, assets, agg shock, age)
	allocate(VN(nj*nbi,ndi*nai,nd,ne,na,nz,TT), stat=val_sol%alloced)
	allocate(VU(nj*nbi,ndi*nai,nd,ne,na,nz,TT), stat=val_sol%alloced)
	allocate(VW(nj*nbi,ndi*nai,nd,ne,na,nz,TT), stat=val_sol%alloced)
	allocate(V(nj*nbi,ndi*nai,nd,ne,na,nz,TT), stat=val_sol%alloced)
	allocate(aN(nj*nbi,ndi*nai,nd,ne,na,nz,TT-1), stat=pol_sol%alloced)
	allocate(aW(nj*nbi,ndi*nai,nd,ne,na,nz,TT-1), stat=pol_sol%alloced)
	allocate(aU(nj*nbi,ndi*nai,nd,ne,na,nz,TT-1), stat=pol_sol%alloced)
	allocate(gwork(nj*nbi,ndi*nai,nd,ne,na,nz,TT-1), stat=pol_sol%alloced)
	allocate(gapp(nj*nbi,ndi*nai,nd,ne,na,nz,TT-1), stat=pol_sol%alloced)

	allocate(gapp_dif(nj*nbi,ndi*nai,nd,ne,na,nz,TT), stat=pol_sol%alloced)
	allocate(gwork_dif(nj*nbi,ndi*nai,nd,ne,na,nz,TT), stat=pol_sol%alloced)


	allocate(maxer(na,nz,ne,nd,nai))

	narg_in = iargc()

	verbose = 4
	print_lev = 4

	! assign structures
	val_sol%VR = VR 
	val_sol%VD = VD 
	val_sol%VW = VW 
	val_sol%VU = VU 
	val_sol%VN = VN
	val_sol%V = V		 

	pol_sol%aR = aR
	pol_sol%aD = aD
	pol_sol%aW = aW
	pol_sol%aN = aN
	pol_sol%aU = aU
	pol_sol%gapp	= gapp
	pol_sol%gwork	= gwork
	pol_sol%gapp_dif= gapp_dif
	pol_sol%gwork_dif= gwork_dif
	

	call setparams()
	agrid(1) = .05*(agrid(1)+agrid(2))
	call vec2csv(agrid,"agrid.csv")

	!**** for diagnostic calculate wages at each wage-determining levels
	open(1, file="wage_dist.csv")
	ibi =1
	iz  =2
	do it = 1,TT-1
		do iai =1,nai
			do id = 1,nd-1
				wagehere = wage(beti(ibi),alfi(iai),id,zgrid(iz),it)
				write(1, "(G20.12)", advance='no') wagehere
			enddo
			id = nd
			wagehere = wage(beti(ibi),alfi(iai),id,zgrid(iz),it)
			write(1,*) wagehere
		enddo
		write(1,*) " "! trailing space
	enddo	
	close(1)


	call sol(val_sol,pol_sol)


!    .----.   @   @
!   / .-"-.`.  \v/
!   | | '\ \ \_/ )
! ,-\ `-.' /.'  /
!'---`----'----'
	!****************************************************************************!
	! IF you love something.... 
	!****************************************************************************!
	deallocate(aR,aD,aN, aU, aW,gwork, gapp)
	deallocate(aiD,maxer)
	deallocate(VR,VD,VN,VU,VW,V)
	deallocate(gwork_dif,gapp_dif)



End PROGRAM






