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
	! Value Functions
	!************************************************************************************************!
		real(8)  ::	Vtest1, Vtest2, Vapp, VC, VR0(nd,ne,na), VR(nd,ne,na), &
				VD0(nd,ne,na,TT), VD(nd,ne,na,TT), &
				VN(nj,nbi,ndi,nai,nd,ne,na,nz,TT), VN0(nj,nbi,ndi,nai,nd,ne,na,nz,TT), &
				VW(nj,nbi,ndi,nai,nd,ne,na,nz,TT), VW0(nj,nbi,ndi,nai,nd,ne,na,nz,TT), &
				 V(nj,nbi,ndi,nai,nd,ne,na,nz,TT),  V0(nj,nbi,ndi,nai,nd,ne,na,nz,TT)
	

	!************************************************************************************************!
	! Policies objects
	!************************************************************************************************!
		real(8)  ::	aR(nd,ne,na), aD(nd,ne,na,TT-1), app2, &
				aN(nj,nbi,ndi,nai,nd,ne,na,nz,TT-1), aW(nj,nbi,ndi,nai,nd,ne,na,nz,TT-1), &
				gwork(nj,nbi,ndi,nai,nd,ne,na,nz,TT-1), gapp(nj,nbi,ndi,nai,nd,ne,na,nz,TT-1)

	!************************************************************************************************!
	! Other
	!************************************************************************************************!
		real(8)	:: junk, summer, eprime, yL, yH 

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

		VD(:,:,:,TT) = VR
			DO ij=1,nj
			DO ibi=1,nbi
			DO idi=1,ndi
			DO iai=1,nai
			DO iz=1,nz   
				VW(ij,ibi,idi,iai,:,:,:,iz,TT) = VR
				VN(ij,ibi,idi,iai,:,:,:,iz,TT) = VR
				! V(ij,ibi,idi,iai,:,:,:,iz,TT) = VR	   
			EndDO
			EndDO
			EndDO
			EndDO
			EndDO


End PROGRAM






