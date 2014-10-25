! V0main.f90

!************************************************************************************************!
! @ Amanda Michaud, v1: 10/6/2014; current: 10/7/2014
!-----------------------------------------------------
!This is a sub-program for DIchoices paper w/ David Wiczer
!	Objective: compute V0(j), the expected value of choosing occupation j in {1,2...J}
!				  before individual types are drawn.
!	
!************************************************************************************************!

program V0main
use nr_lib, GKSTW_tools
use V0para
implicit none

call setparams()
!************************************************************************************************!
! Counters and Indicies
!************************************************************************************************!

	integer  :: i, j, t, ia, ie, id, it, iaa, apol 

	real(DP), parameter :: pival = 3.14159265

!************************************************************************************************!
! Value Functions
!************************************************************************************************!
	real(DP)  ::	Vtest1, Vtest2, VR0(nd,ne,na), VR(nd,ne,na), &
			VD0(nd,ne,na,TT-1), VD(nd,ne,na,TT-1)

!************************************************************************************************!
! Policies objects
!************************************************************************************************!
	real(DP)	::  

!************************************************************************************************!
! Other
!************************************************************************************************!
	real(DP)	:: junk, summer 


!************************************************************************************************!
! Initialize things
!************************************************************************************************!
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
		VR0(id,ie,ia) = ((exval**(theta*real(id)))*(SSI(egrid(ie))+agrid(ia)))**(1-gam))/((1-gam)*((1-junk)**(1-gam))*(1-beta*ptau(TT)*(junk*R)**(1-gam)))
	ENDdo
	ENDdo
	ENDdo
	
	DO WHILE (j<maxiter)
	summer = 0
	  	DO ie=1,ne
		  apol = 1
	  	DO ia=1,na
		  Vtest1 = util(SSI(egrid(ie))+R*agrid(ia)-agrid(apol))+beta*ptau(TT)*VR0(id,ie,apol)
		  iaa = apol+1
		DO WHILE (iaa<na)
		  Vtest2 = util(SSI(egrid(ie))+R*agrid(ia)-agrid(iaa))+beta*ptau(TT)*VR0(id,ie,iaa)
		  IF Vtest2<Vtest1
		     apol = iaa		!concave, start next loop here
		     iaa = na
		  ELSE
		     Vtest1 = Vtest2	
		  EndIF
		  iaa = iaa+1
		EndDO
		VR(1,ie,ia) = Vtest1
		summer = summer+ abs(VR(1,ie,ia)-VR0(1,ie,ia))	
		!Polices for disabled are the same, just scale V-function
		VR(2,ie,ia) = VR(1,ie,ia)*((exval**(theta*1))**(1-gam))
		VR(3,ie,ia) = VR(1,ie,ia)*((exval**(theta*2))**(1-gam))
		EndDO			
		EndDO
	 IF summer < Vtol
	  j = maxiter+100	!Converged
	 EndIF
	  VR0 = VR	!New guess
	  j=j+1
	EndDO

!2) Calculate Value of Disabled: VD(d,e,a,t)	!<---------------THIS IS WHERE YOU LEFT OFF AMANDA! 
	!d in{1,2,3}  	   :disability extent
	!e inR+       	   :earnings index
	!a inR+	      	   :asset holdings
	!t in[1,2...TT-1]  :age
	
	j=1
	DO WHILE (j<maxiter)
	summer = 0
		DO it = 1,TT-1
	  	DO ie=1,ne
		   apol = 1

		!Guess will be value at t+1
		IF (j==1)
		DO ia=1,na
		IF (it==1)
		   VD0(1,ie,ia,TT-1) = VR(1,ie,ia)	
		ELSE
		   VD0(1,ie,ia,TT-it) = VD0(1,ie,ia,TT-it+1)
		EndIF
		EndDO
		EndIF

		!Loop over current state: assets
	  	DO ia=1,na
		  Vtest1 = util(SSDI(egrid(ie))+R*agrid(ia)-agrid(apol))+beta*((1-ptau(TT-it))*VD0(id,ie,apol,it+1)+ptau(TT-it)*VD0(id,ie,apol,it))
		  iaa = apol+1
		!Find Policy
		DO WHILE (iaa<na)
		  Vtest2 = util(SSDI(egrid(ie))+R*agrid(ia)-agrid(iaa))+beta*((1-ptau(TT-it))*VD0(id,ie,iaa,it+1)+ptau(TT-it)*VD0(id,ie,iaa,it))
		  IF Vtest2<Vtest1
		     apol = iaa		!concave, start next loop here
		     iaa = na
		  ELSE
		     Vtest1 = Vtest2	
		  EndIF
		  iaa = iaa+1
		EndDO

		VD(1,ie,ia,it) = Vtest1
		summer = summer+ abs(VD(1,ie,ia,it)-VD0(1,ie,ia,it))	
		!Polices for disabled are the same, just scale V-function
		VD(2,ie,ia,it) = VD(1,ie,ia)*((exval**(theta*1))**(1-gam))
		VD(3,ie,ia,it) = VD(1,ie,ia)*((exval**(theta*2))**(1-gam))
		EndDO			
		EndDO
	 IF summer < Vtol
	  j = maxiter+100	!Converged
	 EndIF
	  VR0 = VR	!New guess
	  j=j+1
	EndDO



!************************************************************************************************!
! Begin loop over occupations
	!DO j=1,nj
!************************************************************************************************!


!Entire program goes here.

!ENDDO



!************************************************************************************************!
! Simulation objects
!************************************************************************************************!






!************************************************************************************************!
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX!
!************************************************************************************************!
contains
!**********************************************************!
!Public Policy Functions
!	1)UI(e)		Unemployment Insurance
!	2)SSDI(e,a)	DI (can be asset tested)
!	3)SSI(e)	Social Sec. Retirement	
!Utility Functions
!	4)u(c,d,w)	w=1 if working; 2 if not
!**********************************************************!
!------------------------------------------------------------------------
!1)UI(e): Unemployment Insurance
!----------------
function UI(ein)	
!----------------
use V0para
implicit none
	
	real(DP), intent(in)	:: ein 

	!I'm using a replacement rate of 30% for now, can be fancier
	UI = ein*UIrr

end function
!------------------------------------------------------------------------
!2) DI (can be asset tested)
!---------------------
function SSDI(ein,ain)
!---------------------
use V0para
implicit none
	
	real(DP), intent(in)	:: ein, ain 

	!Follows Pistafferi & Low '12
	IF (ein<DItest1) then
	 SSDI = 0.9*ein
	ELSEIF (ein<DItest2) then
	 SSDI = 0.9*DItest1 + 0.32*(ein-DItest1)
	ELSEIF (ein<DItest3) then
	 SSDI = 0.9*DItest1 + 0.32*(ein-DItest1)+0.15*(ein-DItest2)
	ELSE
	 SSDI = 0.9*DItest1 + 0.32*(ein-DItest1)+0.15*(DItest3-DItest2)
	ENDIF

end function
!------------------------------------------------------------------------
!3) Social Sec. Retirement	
!--------------------
function SSI(ein)
!--------------------
use V0para
implicit none
	
	real(DP), intent(in)	:: ein, ain 

	!Follows Pistafferi & Low '12
	IF (ein<DItest1) then
	 SSI = 0.9*ein
	ELSEIF (ein<DItest2) then
	 SSI = 0.9*DItest1 + 0.32*(ein-DItest1)
	ELSEIF (ein<DItest3) then
	 SSI = 0.9*DItest1 + 0.32*(ein-DItest1)+0.15*(ein-DItest2)
	ELSE
	 SSI = 0.9*DItest1 + 0.32*(ein-DItest1)+0.15*(DItest3-DItest2)
	ENDIF

end function

!------------------------------------------------------------------------
! 4) Utility Function
!--------------------
function util(cin,din,win)
!--------------------
use V0para
implicit none
	
	real(DP), intent(in)	:: cin
	integer, intent(in)	:: din, win

	if (win==1)
	  util = ((cin*(exval**(theta*real(din-1)+eta)))**(1-gam))/(1-gam)
	else
	  util = ((cin*(exval**(theta*real(din-1))))**(1-gam))/(1-gam)
	end 

end function

END PROGRAM




