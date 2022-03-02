   subroutine da_varbc_synop_direct (iv,ob)

   !-----------------------------------------------!
   !  Apply bias correction to SYNOP innovations  !
   !-----------------------------------------------!

   implicit none

   type (iv_type), intent(inout)   :: iv
   type (y_type),  intent(inout)   :: ob

   real                            :: bc,bias
   real                            :: contri(5)
   integer                         :: i,iobs,ipred,npred
   character(len=30)               :: stringn


   if (trace_use) call da_trace_entry("da_varbc_synop_direct")

   if (rootproc) &
       write(unit=varbc_synop_unit,fmt='(//5X,A/)')'Calculating corrected innovation'

       ! u wind

       bias  = 0.
       npred = iv%varbc_synop%u%npred

       do iobs = iv%info(synop)%n1, iv%info(synop)%n2
          if (abs(iv%synop(iobs)%u%inv-missing_r) > 1.0) then
             bc = dot_product( iv%varbc_synop%u%param(1:npred),  &
                               iv%varbc_synop%u%pred (1:npred,iobs) )
             if(iv%info(synop)%proc_domain(1,iobs)) bias = bias + iv%synop(iobs)%u%inv

             iv%synop(iobs)%u%inv = iv%synop(iobs)%u%inv - bc
	  end if
       end do

       if (iv%varbc_synop%u%nobs_sum > 0) &
           bias = wrf_dm_sum_real(bias)/iv%varbc_synop%u%nobs_sum

       if (rootproc) then
           write(unit=varbc_synop_unit,fmt='(10X,A,2ES15.7)') &
                'bias corrected for synop u: (BC,OMB)', bc, bias
       end if

       ! v wind

       bias  = 0.
       npred = iv%varbc_synop%v%npred

       do iobs = iv%info(synop)%n1, iv%info(synop)%n2
          if (abs(iv%synop(iobs)%v%inv-missing_r) > 1.0) then
             bc = dot_product( iv%varbc_synop%v%param(1:npred),  &
                               iv%varbc_synop%v%pred (1:npred,iobs) )
             if(iv%info(synop)%proc_domain(1,iobs)) bias = bias + iv%synop(iobs)%v%inv

             iv%synop(iobs)%v%inv = iv%synop(iobs)%v%inv - bc
	  end if
       end do

       if (iv%varbc_synop%v%nobs_sum > 0) &
           bias = wrf_dm_sum_real(bias)/iv%varbc_synop%v%nobs_sum

       if (rootproc) then
           write(unit=varbc_synop_unit,fmt='(10X,A,2ES15.7)') &
                'bias corrected for synop v: (BC,OMB)', bc, bias
       end if


       ! temperature

       bias  = 0.
       npred = iv%varbc_synop%t%npred

       do iobs = iv%info(synop)%n1, iv%info(synop)%n2
          if (abs(iv%synop(iobs)%t%inv-missing_r) > 1.0) then
             bc = dot_product( iv%varbc_synop%t%param(1:npred),  &
                               iv%varbc_synop%t%pred (1:npred,iobs) )
             if(iv%info(synop)%proc_domain(1,iobs)) bias = bias + iv%synop(iobs)%t%inv

             iv%synop(iobs)%t%inv = iv%synop(iobs)%t%inv - bc
	  end if
       end do

       if (iv%varbc_synop%t%nobs_sum > 0) &
           bias = wrf_dm_sum_real(bias)/iv%varbc_synop%t%nobs_sum

       if (rootproc) then
           write(unit=varbc_synop_unit,fmt='(10X,A,2ES15.7)') &
                'bias corrected for synop t: (BC,OMB)', bc, bias
       end if


       ! qvapor

       bias  = 0.
       npred = iv%varbc_synop%q%npred

       do iobs = iv%info(synop)%n1, iv%info(synop)%n2
          if (abs(iv%synop(iobs)%q%inv-missing_r) > 1.0) then
             bc = dot_product( iv%varbc_synop%q%param(1:npred),  &
                            iv%varbc_synop%q%pred (1:npred,iobs) )
             if(iv%info(synop)%proc_domain(1,iobs)) bias = bias + iv%synop(iobs)%q%inv

             iv%synop(iobs)%q%inv = iv%synop(iobs)%q%inv - bc
	  end if
       end do

       if (iv%varbc_synop%q%nobs_sum > 0) &
           bias = wrf_dm_sum_real(bias)/iv%varbc_synop%q%nobs_sum

       if (rootproc) then
           write(unit=varbc_synop_unit,fmt='(10X,A,2ES15.7)') &
                'bias corrected for synop q: (BC,OMB)', bc, bias
       end if


       if (rootproc) write(unit=varbc_synop_unit,fmt='(/10X,A)') &
                   "Variable Predictor/Parameter/Bias model (npred)  &  BC"

   if (trace_use) call da_trace_exit("da_varbc_synop_direct")

   end subroutine da_varbc_synop_direct
