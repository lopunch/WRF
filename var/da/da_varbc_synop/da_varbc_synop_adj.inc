   subroutine da_varbc_synop_adj (cv_size, cv, iv, y)

   implicit none

   integer, intent(in)           :: cv_size    
   real, intent(inout)           :: cv(1:cv_size)
   type (iv_type), intent(inout) :: iv        
   type (y_type), intent(in)     :: y             ! y = h (xa)

   integer                       :: iobs,ipred,npred
   real, allocatable             :: varbc_param_adj(:)
   real                          :: yh, cv_local


   if (trace_use) call da_trace_entry("da_varbc_synop_adj")

   !qvapor

   npred = iv%varbc_synop%q%npred

   allocate( varbc_param_adj(npred) )

   if (iv%varbc_synop%q%nobs_sum >= varbc_synop_nobsmin) then

       varbc_param_adj(:) = 0.0

       if (iv%varbc_synop%q%nobs > 0) then
           do iobs = iv%info(synop)%n1, iv%info(synop)%n2
              if (iv%synop(iobs)%q%qc >= 0 .and. iv%info(synop)%proc_domain(1,iobs)) then
                 varbc_param_adj(1:npred) = varbc_param_adj(1:npred) + &
                     y%synop(iobs)%q * iv%varbc_synop%q%pred(1:npred,iobs)
	      end if
           end do
       end if

       do ipred = 1, npred
          varbc_param_adj(ipred) = wrf_dm_sum_real(varbc_param_adj(ipred))
       end do

       if(rootproc) then
           do ipred = 1, npred
              cv_local = dot_product( varbc_param_adj(1:npred), iv%varbc_synop%q%vtox(ipred,1:npred) )
              cv(iv%varbc_synop%q%index(ipred)) = cv(iv%varbc_synop%q%index(ipred)) + cv_local
           end do
       end if

   end if

   deallocate(varbc_param_adj)	     

   !temperature

   npred = iv%varbc_synop%t%npred

   allocate( varbc_param_adj(npred) )

   if (iv%varbc_synop%t%nobs_sum >= varbc_synop_nobsmin) then

       varbc_param_adj(:) = 0.0

       if (iv%varbc_synop%t%nobs > 0) then
           do iobs = iv%info(synop)%n1, iv%info(synop)%n2
              if (iv%synop(iobs)%t%qc >= 0 .and. iv%info(synop)%proc_domain(1,iobs)) then
                 varbc_param_adj(1:npred) = varbc_param_adj(1:npred) + &
                     y%synop(iobs)%t * iv%varbc_synop%t%pred(1:npred,iobs)
	      end if 
           end do
       end if

       do ipred = 1, npred
          varbc_param_adj(ipred) = wrf_dm_sum_real(varbc_param_adj(ipred))
       end do

       if(rootproc) then
           do ipred = 1, npred
              cv_local = dot_product( varbc_param_adj(1:npred), iv%varbc_synop%t%vtox(ipred,1:npred) )
              cv(iv%varbc_synop%t%index(ipred)) = cv(iv%varbc_synop%t%index(ipred)) + cv_local
           end do
       end if

   end if

   deallocate(varbc_param_adj)	     


   !v wind

   npred = iv%varbc_synop%v%npred

   allocate( varbc_param_adj(npred) )

   if (iv%varbc_synop%v%nobs_sum >= varbc_synop_nobsmin) then

       varbc_param_adj(:) = 0.0

       if (iv%varbc_synop%v%nobs > 0) then
           do iobs = iv%info(synop)%n1, iv%info(synop)%n2
              if (iv%synop(iobs)%v%qc >= 0 .and. iv%info(synop)%proc_domain(1,iobs)) then
                 varbc_param_adj(1:npred) = varbc_param_adj(1:npred) + &
                     y%synop(iobs)%v * iv%varbc_synop%v%pred(1:npred,iobs)
	      end if 
           end do
       end if

       do ipred = 1, npred
          varbc_param_adj(ipred) = wrf_dm_sum_real(varbc_param_adj(ipred))
       end do

       if(rootproc) then
           do ipred = 1, npred
              cv_local = dot_product( varbc_param_adj(1:npred), iv%varbc_synop%v%vtox(ipred,1:npred) )
              cv(iv%varbc_synop%v%index(ipred)) = cv(iv%varbc_synop%v%index(ipred)) + cv_local
           end do
       end if

   end if

   deallocate(varbc_param_adj)	     


   !u wind

   npred = iv%varbc_synop%u%npred

   allocate( varbc_param_adj(npred) )

   if (iv%varbc_synop%u%nobs_sum >= varbc_synop_nobsmin) then

       varbc_param_adj(:) = 0.0

       if (iv%varbc_synop%u%nobs > 0) then
           do iobs = iv%info(synop)%n1, iv%info(synop)%n2
              if (iv%synop(iobs)%u%qc >= 0 .and. iv%info(synop)%proc_domain(1,iobs)) then
                 varbc_param_adj(1:npred) = varbc_param_adj(1:npred) + &
                     y%synop(iobs)%u * iv%varbc_synop%u%pred(1:npred,iobs)
	      end if 
           end do
       end if

       do ipred = 1, npred
          varbc_param_adj(ipred) = wrf_dm_sum_real(varbc_param_adj(ipred))
       end do

       if(rootproc) then
           do ipred = 1, npred
              cv_local = dot_product( varbc_param_adj(1:npred), iv%varbc_synop%u%vtox(ipred,1:npred) )
              cv(iv%varbc_synop%u%index(ipred)) = cv(iv%varbc_synop%u%index(ipred)) + cv_local
           end do
       end if

   end if

   deallocate(varbc_param_adj)	     


   if (trace_use) call da_trace_exit("da_varbc_synop_adj")

   end subroutine da_varbc_synop_adj
