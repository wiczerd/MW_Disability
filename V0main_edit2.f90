! V0main.f90

!************************************************************************************************!
! @ Amanda Michaud, v1: 10/6/2014; current: 10/31/2014
!-----------------------------------------------------
!This is a sub-program for DIchoices paper w/ David Wiczer
!	Objective: compute V0(j), the expected value of choosing occupation j in {1,2...J}
!				  before individual types are drawn.
!	
!************************************************************************************************!
module helper_funs
	
	use V0para
	implicit none
	contains
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
	!**********************************************************!
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
		
		if (win .EQ. 1) THEN
		  util = ((cin*(exval**(theta*real(din-1)+eta)))**(1-gam))/(1-gam)
		else
		  util = ((cin*(exval**(theta*real(din-1))))**(1-gam))/(1-gam)
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

		wage = aiin+biin*zin+wd(din)+wtau(tin)

	end function

	!------------------------------------------------------------------------
	! 6) Locate Function
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
		else if (x >= xx(nf)-epsilon(xx(nf))) then
			locate=nf-1
		else
			locate=il
		end if
		finder = locate
	end function


end module helper_funs

!**************************************************************************************************************!
!**************************************************************************************************************!
!						MAIN PROGRAM						       !
!**************************************************************************************************************!
!**************************************************************************************************************!

program V0main
	!INCLUDE 'GKSTW_tools.f90'
	!INCLUDE 'nrtype.f90'
	!INCLUDE 'V0para.f90'
	use V0para
	use helper_funs
	!use nrtype
	!use GKSTW_tools 
	implicit none


	!************************************************************************************************!
	! Counters and Indicies
	!************************************************************************************************!

		integer  :: i, j, t, ia, ie, id, it, iaa, apol, ibi, iai, ij , idi, izz, iaai, idd, &
			    iee1, iee2, iz, unitno
	
	!************************************************************************************************!
	! Value Functions- Stack z-risk j and indiv. exposure beta_i
	!************************************************************************************************!
	real(8)  	  ::Vtest1, Vtest2, Vapp, VC, app2
	real(8), allocatable :: VR0(:,:,:), VR(:,:,:), &
				VD0(:,:,:,:), VD(:,:,:,:), &
				VN(:,:,:,:,:,:,:,:), VN0(:,:,:,:,:,:,:,:), &
				VW(:,:,:,:,:,:,:,:), VW0(:,:,:,:,:,:,:,:), &
				V(:,:,:,:,:,:,:,:),  V0(:,:,:,:,:,:,:,:)
	

	!************************************************************************************************!
	! Policies objects- Stack z-risk j and indiv. exposure beta_i
	!************************************************************************************************!
	real(8), allocatable  ::aR(:,:,:), aD(:,:,:,:),  &
				aN(:,:,:,:,:,:,:,:), aW(:,:,:,:,:,:,:,:), &
				gwork(:,:,:,:,:,:,:,:), gapp(:,:,:,:,:,:,:,:)

	!************************************************************************************************!
	! Other
	!************************************************************************************************!
		real(8)	:: junk, summer, eprime, yL, yH 

	!************************************************************************************************!
	! Allocate phat matrices
	!************************************************************************************************!
	
	allocate(VR0(nd,ne,na))
	allocate(VR(nd,ne,na))
	allocate(aR(nd,ne,na))
	allocate(VD0(nd,ne,na,TT))
	allocate(VD(nd,ne,na,TT))
	allocate(aD(nd,ne,na,TT-1))
	allocate(VN(nj*nbi,ndi,nai,nd,ne,na,nz,TT))
	allocate(VN0(nj*nbi,ndi,nai,nd,ne,na,nz,TT))
	allocate(VW(nj*nbi,ndi,nai,nd,ne,na,nz,TT))
	allocate(VW0(nj*nbi,ndi,nai,nd,ne,na,nz,TT))
	allocate(V(nj*nbi,ndi,nai,nd,ne,na,nz,TT))
	allocate(V0(nj*nbi,ndi,nai,nd,ne,na,nz,TT))
	allocate(aN(nj*nbi,ndi,nai,nd,ne,na,nz,TT-1))
	allocate(aW(nj*nbi,ndi,nai,nd,ne,na,nz,TT-1))
	allocate(gwork(nj*nbi,ndi,nai,nd,ne,na,nz,TT-1))
	allocate(gapp(nj*nbi,ndi,nai,nd,ne,na,nz,TT-1))

	call setparams()
	!************************************************************************************************!
	! Caculate things that are independent of occupation/person type
	!	1) Value of Retired:  VR(d,e,a)
	!	2) Value of Disabled: VD(d,e,a)
	!************************************************************************************************!

	!1) Calculate Value of Retired: VR(d,e,a)
		!d in{1,2,3}  : disability extent
		!e inR+       :	earnings index
		!a inR+	      : asset holdings
	
		!VFI with good guess
		!Initialize
		junk = (1.0+(beta*ptau(TT)*R**(1.0-gam))**(-1.0/gam))**-1.0 
		DO ie=1,ne
		DO ia=1,na
		DO id=1,nd

			VR0(id,ie,ia) = (((exval**(theta*real(id)))*(SSI(egrid(ie))+agrid(ia)))**(1-gam))/((1-gam)*((1-junk)**(1-gam))*(1-beta*ptau(TT)*(junk*R)**(1-gam)))
		ENDdo
		ENDdo
		ENDdo
	
		DO WHILE (j<maxiter)
		summer = 0
		  	DO ie=1,ne
			  apol = 1
		  	DO ia=1,na
			  Vtest1 = util(SSI(egrid(ie))+R*agrid(ia)-agrid(apol),1,2)+beta*ptau(TT)*VR0(id,ie,apol)
			  iaa = apol+1
			DO WHILE (iaa<na)
			  Vtest2 = util(SSI(egrid(ie))+R*agrid(ia)-agrid(iaa),1,1)+beta*ptau(TT)*VR0(id,ie,iaa)
	  		  apol = max(iaa-1,1)	!concave, start next loop here
			  IF (Vtest2<Vtest1) THEN
			     iaa = na
			  ELSE
			     Vtest1 = Vtest2	
			  EndIF
			  iaa = iaa+1
			EndDO
			VR(1,ie,ia) = Vtest1
			aR(1,ie,ia) = agrid(min(na,apol+1))
			summer = summer+ abs(VR(1,ie,ia)-VR0(1,ie,ia))	
			!Polices for disabled are the same, just scale V-function
			VR(2,ie,ia) = VR(1,ie,ia)*((exval**(theta*1))**(1-gam))
			VR(3,ie,ia) = VR(1,ie,ia)*((exval**(theta*2))**(1-gam))
			EndDO			
			EndDO
		 IF (summer < Vtol) THEN
		  j = maxiter+100	!Converged
		 EndIF
		  VR0 = VR	!New guess
		  j=j+1
		EndDO

		!----------------------------------------------------------!
		!Set value at t=TT to be VR in all other V-functions
		!----------------------------------------------------------!
		DO id=1,nd
		DO ie=1,ne
		DO ia=1,na
	  	   VD(id,ie,ia,TT) = VR(id,ie,ia)
			DO ij=1,nj
			DO ibi=1,nbi
			DO idi=1,ndi
			DO iai=1,nai
			!DO iz=1,nz   
				VW((ij-1)*nbi+ibi,idi,iai,id,ie,ia,1,TT) = VR(id,ie,ia)
				VN((ij-1)*nbi+ibi,idi,iai,id,ie,ia,1,TT) = VR(id,ie,ia)
				V((ij-1)*nbi+ibi,idi,iai,id,ie,ia,1,TT) = VR(id,ie,ia)	   
			!EndDO
			EndDO
			EndDO
			EndDO
			EndDO
		EndDO
		EndDO
		EndDO

	       !----------------------------------------------------------------!
	!2) Calculate Value of Disabled: VD(d,e,a,t)	 
		!d in{1,2,3}  	   :disability extent
		!e inR+       	   :earnings index
		!a inR+	      	   :asset holdings
		!t in[1,2...TT-1]  :age
	
		!Work backwards from TT
		DO it = 1,TT-1
			!Guess will be value at t+1
			DO ia=1,na
			DO ie=1,ne
			   VD0(1,ie,ia,TT-it) = VD(1,ie,ia,TT-it+1)
			EndDO
			EndDO

		!Loop to find V(..,it) as fixed point
		 j=1
		 DO WHILE (j<maxiter)
		  summer = 0
			!Loop over earnings index
		  	DO ie=1,ne
			   apol = 1
			!Loop over current state: assets
		  	DO ia=1,na
			  Vtest1 = util(SSDI(egrid(ie))+R*agrid(ia)-agrid(apol),1,2)+beta*((1-ptau(TT-it))*VD0(id,ie,apol,TT-it+1)+ptau(TT-it)*VD0(id,ie,apol,TT-it))
			  iaa = apol+1
			!Find Policy
			DO WHILE (iaa<na)
			  Vtest2 = util(SSDI(egrid(ie))+R*agrid(ia)-agrid(iaa),1,2)+beta*((1-ptau(TT-it))*VD0(id,ie,iaa,TT-it+1)+ptau(TT-it)*VD0(id,ie,iaa,TT-it))
			  apol = max(iaa-1,1)		!concave, start next loop here
			  IF (Vtest2<Vtest1) THEN
			     iaa = na
			  ELSE
			     Vtest1 = Vtest2	
			  EndIF
			  iaa = iaa+1
			EndDO	!iaa

			VD(1,ie,ia,TT-it) = Vtest1
			aD(1,ie,ia,TT-it) = agrid(min(apol+1,na))
			summer = summer+ abs(VD(1,ie,ia,TT-it)-VD0(1,ie,ia,TT-it))	
			!Polices for disabled are the same, just scale V-function
			VD(2,ie,ia,TT-it) = VD(1,ie,ia,TT-it)*((exval**(theta*1))**(1-gam))
			VD(3,ie,ia,TT-it) = VD(1,ie,ia,TT-it)*((exval**(theta*2))**(1-gam))
			EndDO	!ie
			EndDO	!ia		

		 IF (summer < Vtol) THEN
		  j = maxiter+100	!Converged
		 EndIF
		  VD0 = VD	!New guess
		  j=j+1
		EndDO	!j: V-iter loop
		EndDO	!t loop, going backwards

	!************************************************************************************************!
	!3) Calculate V= max(VW,VN); requires calculating VW and VN
	! Begin loop over occupations
		DO ij = 1,nj
	! And betas
		DO ibi = 1,nbi 
	! And individual disability type
		DO idi = 1,ndi
	!************************************************************************************************!
		!Work Backwards TT-1,TT-2...1!
		DO it=1,TT-1

		!----Initialize---!
		DO iai=1,nai
		DO id =1,nd
		DO ie =1,ne
		DO ia =1,na
		DO iz =1,nz
		!Guess once, then use prior occupation/beta as guess
			IF (j .EQ. 1 .AND. ibi .EQ. 1 .AND. idi .EQ. 1) THEN
			 !0) Guess VW0(nj,nbi,nai,nd,ne,na,nz,TT-1)
				VW0((ij-1)*nbi+ibi,idi,iai,id,ie,ia,iz,TT-it) = VW((ij-1)*nbi+ibi,idi,iai,id,ie,ia,iz,TT-it+1)
			 !0) Guess VN0(nj,nbi,nai,nd,ne,na,nz,TT-1)
				VN0((ij-1)*nbi+ibi,idi,iai,id,ie,ia,iz,TT-it) = VN0((ij-1)*nbi+ibi,idi,iai,id,ie,ia,iz,TT-it+1)
			IF (it .EQ. 1) THEN
				VN0((ij-1)*nbi+ibi,idi,iai,1,ie,ia,iz,TT-it) = 0.5*VW0((ij-1)*nbi+ibi,idi,iai,1,ie,ia,iz,TT-it+1)
			EndIF
			ELSE
			 !0) Guess VW0(nj,nbi,nai,nd,ne,na,nz,TT-1)
				VW0((ij-1)*nbi+ibi,idi,iai,id,ie,ia,iz,TT-it) = VW0(max(1,ij-1)*nbi+max(1,ibi-1),max(1,idi-1),iai,id,ie,ia,iz,TT-it+1)
			 !0) Guess VN0(nj,nbi,nai,nd,ne,na,nz,TT-1)
				VN0((ij-1)*nbi+ibi,idi,iai,id,ie,ia,iz,TT-it) = VN0(max(1,ij-1)*nbi+max(1,ibi-1),max(1,idi-1),iai,id,ie,ia,iz,TT-it+1)
			EndIF !first occupationXbeta loop
			!0) Calculate V0(nj,nbi,nai,nd,ne,na,nz,it)	
			V0((ij-1)*nbi+ibi,idi,iai,id,ie,ia,iz,TT-it)= max(VW0((ij-1)*nbi+ibi,idi,iai,id,ie,ia,iz,TT-it),VN0((ij-1)*nbi+ibi,idi,iai,id,ie,ia,iz,TT-it))

		EndDO	!iai   
		EndDO 	!id
		EndDO	!ie
		EndDO	!ia
		EndDO	!iz
		!***********************************************************************************************************
		!Loop over V=max(VN,VW)	
		 j=1
		 DO WHILE (j<maxiter)
		 summer = 0	!Use to calc |V-V0|<eps
		!------------------------------------------------!
		!Solve VN given guesses on VW, VN, and implied V
		!------------------------------------------------!  
		  	DO iai=1,nai	!Loop over alpha (ai)
			DO id=1,nd	!Loop over disability index
		  	DO ie=1,ne	!Loop over earnings index
		  	DO iz=1,nz	!Loop over TFP
			!Restart at bottom of asset grid for each of the above (ai,d,e,z)
			apol = 1
			!----------------------------------------------------------------
			!Loop over current state: assets
		  	DO ia=1,na
			!Continuation value if don't go on disability
			 VC = 0	 
			  DO izz = 1,nz	 !Loop over z'
			  DO iaai = 1,nai !Loop over alpha_i'
			  VC = VC + beta*piz(iz,izz,idi)*pialf(iai,iaai)*((1-ptau(TT-it))*V((ij-1)*nbi+ibi,idi,iaai,id,ie,apol,izz,TT-it+1)+ptau(TT-it)*V((ij-1)*nbi+ibi,idi,iaai,id,ie,apol,izz,TT-it)) 
			  EndDO
			  EndDO
			 !Continuation if app is accepted
			  Vapp = beta*((1-ptau(TT-it))*VD(id,ie,apol,TT-it+1)+ptau(TT-it)*VD(id,ie,apol,TT-it))
			 !Apply if xi(id)Vapp-nu>xi(id)VC
			 IF (xi(id)*Vapp > xi(id)*VC) THEN
				Vtest1 = util(UI(egrid(ie))+R*agrid(ia)-agrid(apol),1,2)-nu+(1-xi(id))*VC+xi(id)*Vapp
				gapp((ij-1)*nbi+ibi,idi,iai,id,ie,ia,iz,TT-it) = 1
			 ELSE !Don't apply
				Vtest1 = util(UI(egrid(ie))+R*agrid(ia)-agrid(apol),1,2)+VC
				gapp((ij-1)*nbi+ibi,idi,iai,id,ie,ia,iz,TT-it) = 0
			 EndIF

			  !Guess on higher grid points only
			  iaa = apol+1
			DO WHILE (iaa<na)
			!Continuation value if don't go on disability
			 VC = 0	 
			  DO izz = 1,nz	 !Loop over z'
			  DO iaai = 1,nai !Loop over alpha_i'
			  VC = VC + beta*piz(iz,izz,idi)*pialf(iai,iaai)*((1-ptau(TT-it))*V((ij-1)*nbi+ibi,idi,iaai,id,ie,iaa,izz,TT-it+1)+ptau(TT-it)*V((ij-1)*nbi+ibi,idi,iaai,id,ie,iaa,izz,TT-it)) 
			  EndDO
			  EndDO
			 !Continuation if app is accepted
			  Vapp = beta*((1-ptau(TT-it))*VD(id,ie,iaa,TT-it+1)+ptau(TT-it)*VD(id,ie,iaa,TT-it))
			 !Apply if xi(id)Vapp-nu>xi(id)VC
			 IF (xi(id)*Vapp > xi(id)*VC) THEN
				Vtest2 = util(UI(egrid(ie))+R*agrid(ia)-agrid(iaa),1,2)-nu+(1-xi(id))*VC+xi(id)*Vapp
				app2 = 1
			 ELSE !Don't apply
				Vtest2 = util(UI(egrid(ie))+R*agrid(ia)-agrid(iaa),1,2)+VC
				app2 = 0
			 EndIF
			  apol = max(iaa-1,1)		!concave, start next loop here
			  IF (Vtest2<Vtest1) THEN	
			     iaa = na
			  ELSE
			     Vtest1 = Vtest2
			     gapp((ij-1)*nbi+ibi,idi,iai,id,ie,ia,iz,TT-it) = app2	
			  EndIF
			  iaa = iaa+1
			EndDO	!iaa

			VN((ij-1)*nbi+ibi,idi,iai,id,ie,ia,iz,TT-it) = Vtest1
			aN((ij-1)*nbi+ibi,idi,iai,id,ie,ia,iz,TT-it) = agrid(min(apol+1,na))

		  	EndDO !iai
			EndDO !id
		  	EndDO !ie
		  	EndDO !iz
		  	EndDO !ia		

		!------------------------------------------------!
		!Solve VW given guesses on VW, VN, and implied V
		!------------------------------------------------!
		  	DO iai=1,nai	!Loop over alpha (ai)
			DO id=1,nd	!Loop over disability index
		  	DO ie=1,ne	!Loop over earnings index
		  	DO iz=1,nz	!Loop over TFP
			!Earnings evolution independent of choices.
			junk = wage(beti(ibi),alfi(iai),id,zgrid(iz),TT-it)
			eprime = Hearn(TT-it,ie,junk)
			!I'm going to linear interpolate for the portion that blocks off bounds on assets
			iee2 = finder(egrid,eprime)
			iee1 = max(1,iee2-1)
			!Restart at bottom of asset grid for each of the above (ai,d,e,z)
			apol = 1
			!----------------------------------------------------------------
			!Loop over current state: assets
		  	DO ia=1,na
			Vtest1 = util(junk+R*agrid(ia)-agrid(apol),1,1)	!Flow util
			  DO izz = 1,nz	 !Loop over z'
			  DO iaai = 1,nai !Loop over alpha_i'
				!Linearly interpolating on e'
				yL = beta*piz(iz,izz,idi)*pialf(iai,iaai)*((1-ptau(TT-it))*V((ij-1)*nbi+ibi,idi,iaai,id,iee1,apol,izz,TT-it+1)+ptau(TT-it)*V((ij-1)*nbi+ibi,idi,iaai,id,iee1,apol,izz,TT-it))
				yH = beta*piz(iz,izz,idi)*pialf(iai,iaai)*((1-ptau(TT-it))*V((ij-1)*nbi+ibi,idi,iaai,id,iee2,apol,izz,TT-it+1)+ptau(TT-it)*V((ij-1)*nbi+ibi,idi,iaai,id,iee2,apol,izz,TT-it))
				Vtest1 = Vtest1 + yL+(yH-yL)*(eprime-egrid(iee1))/(egrid(iee2)-egrid(iee1))
			  EndDO
			  EndDO		

	  	       !Guess on higher grid points only
			  iaa = apol+1
			!Find Policy
			DO WHILE (iaa<na)
			!Continuation value if don't go on disability
			Vtest2 = util(junk+R*agrid(ia)-agrid(iaa),1,1)	!Flow util
			  DO izz = 1,nz	 !Loop over z'
			  DO iaai = 1,nai !Loop over alpha_i'
				!Linearly interpolating on e'
				yL = beta*piz(iz,izz,idi)*pialf(iai,iaai)*((1-ptau(TT-it))*V((ij-1)*nbi+ibi,idi,iaai,id,iee1,iaa,izz,TT-it+1)+ptau(TT-it)*V((ij-1)*nbi+ibi,idi,iaai,id,iee1,iaa,izz,TT-it))
				yH = beta*piz(iz,izz,idi)*pialf(iai,iaai)*((1-ptau(TT-it))*V((ij-1)*nbi+ibi,idi,iaai,id,iee2,iaa,izz,TT-it+1)+ptau(TT-it)*V((ij-1)*nbi+ibi,idi,iaai,id,iee2,iaa,izz,TT-it))
				Vtest2 = Vtest2 + yL+(yH-yL)*(eprime-egrid(iee1))/(egrid(iee2)-egrid(iee1))
			  EndDO
			  EndDO	
			  apol = max(iaa-1,1)		!concave, start next loop here
			  IF (Vtest2<Vtest1) THEN	
			     iaa = na
			  ELSE
			     Vtest1 = Vtest2
			  EndIF
			  iaa = iaa+1
			EndDO	!iaa

			VW((ij-1)*nbi+ibi,idi,iai,id,ie,ia,iz,TT-it) = Vtest1
			aW((ij-1)*nbi+ibi,idi,iai,id,ie,ia,iz,TT-it) = agrid(min(apol+1,na))


			!------------------------------------------------!
			!Calculate V with solved vals of VW and VN
			!------------------------------------------------!
			IF (VW((ij-1)*nbi+ibi,idi,iai,id,ie,ia,iz,TT-it)>VN((ij-1)*nbi+ibi,idi,iai,id,ie,ia,iz,TT-it)) THEN
			V((ij-1)*nbi+ibi,idi,iai,id,ie,ia,iz,TT-it) = VW((ij-1)*nbi+ibi,idi,iai,id,ie,ia,iz,TT-it)
			gwork((ij-1)*nbi+ibi,idi,iai,id,ie,ia,iz,TT-it) = 1
			ELSE
			V((ij-1)*nbi+ibi,idi,iai,id,ie,ia,iz,TT-it) = VN((ij-1)*nbi+ibi,idi,iai,id,ie,ia,iz,TT-it)
			gwork((ij-1)*nbi+ibi,idi,iai,id,ie,ia,iz,TT-it) = 0
			EndIF
	
			summer = summer+ abs(V((ij-1)*nbi+ibi,idi,iai,id,ie,ia,iz,TT-it)-V0((ij-1)*nbi+ibi,idi,iai,id,ie,ia,iz,TT-it))
		
		  	EndDO !iai
			EndDO !id
		  	EndDO !ie
		  	EndDO !iz
		  	EndDO !ia

		!------------------------------------------------!
		!Check |V-V0|<eps
		!------------------------------------------------!
!WRITE(*,*) summer, j
		 IF (summer < Vtol) THEN
		  j = maxiter+100	!Converged
		 EndIF
		  V0 = V	!New guess
		  j=j+1
		EndDO	!j: V-iter loop
WRITE(*,*) ij, ibi, idi, it
		EndDO	!t loop, going backwards

	EndDO	!idi
	EndDO	!ibi
	EndDO	!j


!WRITE(*,*) aD(1,:,:,TT-5)
!WRITE(*,*) '-----'
!WRITE(*,*) VR
!WRITE(*,*) '-----'
!WRITE(*,*) summer, j

	open (newunit=unitno,file ='workpol.txt',status ='replace')
	ie = INT(na/4)
	write (unitno,*) gwork(1:nj,:,:,:,floor(egrid/3),ie,2,2)	!((ij-1)*nbi+ibi,idi,iai,id,ie,ia,iz,TT-it)
	close (unitno)
	open (newunit=unitno,file ='driskStat.txt',status ='replace')
	j = INT(nj/2)
	write (unitno,*) gwork(j,:,:,:,floor(egrid/3),ie,:,2)	!((ij-1)*nbi+ibi,idi,iai,id,ie,ia,iz,TT-it)
	close (unitno)
	open (newunit=unitno,file ='zriskStat.txt',status ='replace')
	idi = INT(ndi/2)
	write (unitno,*) gwork(1:nj,idi,:,:,floor(egrid/3),ie,:,2)	!((ij-1)*nbi+ibi,idi,iai,id,ie,ia,iz,TT-it)
	close (unitno)
	open (newunit=unitno,file ='lifeStat.txt',status ='replace')
	write (unitno,*) gwork(j,idi,:,:,floor(egrid/3),ie,:,:)	!((ij-1)*nbi+ibi,idi,iai,id,ie,ia,iz,TT-it)
	close (unitno)
	open (newunit=unitno,file ='Vfun.txt',status ='replace')
	write (unitno,*) V(1:nj,:,:,:,:,:,:,:)	
	close (unitno)


!    .----.   @   @
!   / .-"-.`.  \v/
!   | | '\ \ \_/ )
! ,-\ `-.' /.'  /
!'---`----'----'



End PROGRAM






