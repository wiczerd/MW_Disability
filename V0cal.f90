	subroutine cal_dist_nloptwrap(fval, nparam, paramvec, gradvec, need_grad, shk)

		use find_params
		use helper_funs

		implicit none 

		integer :: nparam,need_grad
		real(8) :: fval, paramvec(nparam),gradvec(nparam)
		type(shocks_struct) :: shk
		real(8) :: errvec(nparam),paramwt(nparam),paramvecH(nparam),errvecH(nparam),paramvecL(nparam),errvecL(nparam),gradstep(nparam), fvalH,fvalL
		integer :: i, ii,print_lev_old, verbose_old
		
		
		print_lev_old = print_lev 
		verbose_old = verbose
		print_lev = 0
		verbose = 1

		if(verbose_old >=1) print *, "test parameter vector", paramvec
		
		call cal_dist(paramvec, errvec,shk)

		if(verbose_old >=1) print *, "         error vector", errvec

		paramwt = 1./dble(nparam)		! equal weight
		fval = 0.
		do i = 1,nparam
			fval = errvec(i)**2*paramwt(i) + fval
		enddo
		if( need_grad .ne. 0) then

			do i=1,nparam
				gradstep (i) = min( dabs( paramvec(i)*(5.e-5_dp) ) ,5.e-5_dp)
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
			
		endif
		print_lev = print_lev_old
		verbose = verbose_old	

	end subroutine cal_dist_nloptwrap

