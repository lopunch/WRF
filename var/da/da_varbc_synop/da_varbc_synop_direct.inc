   subroutine da_varbc_synop_direct (iv,ob)

   !-----------------------------------------------!
   !  Apply bias correction to TAMDAR innovations  !
   !-----------------------------------------------!

   implicit none

   type (iv_type), intent(inout)   :: iv
   type (y_type),  intent(inout)   :: ob

   real                            :: bc,bias
   real                            :: contri(5)
   integer                         :: i,iobs,ipred,isn,npred
   character(len=30)               :: stringn


   if (trace_use) call da_trace_entry("da_varbc_synop_direct")

   if (rootproc) &
       write(unit=varbc_synop_unit,fmt='(//5X,A/)')'Calculating corrected innovation'

       ! u wind

       bias  = 0.
       npred = iv%varbc_synop%u%npred

       if (iv%varbc_synop%u%nobs > 0) then 
           do iobs = 1, iv%varbc_synop%u%nobs
              bc = dot_product( iv%varbc_synop%u%param(1:npred),  &
                                iv%varbc_synop%u%pred (1:npred,iobs) )
              isn  = iv%varbc_synop%u%obs_sn(iobs)
              bias = bias + iv%synop(isn)%u%inv

              iv%synop(isn)%u%inv = iv%synop(isn)%u%inv - bc
           end do
       end if

       if (iv%varbc_synop%u%nobs_sum > 0) &
           bias = wrf_dm_sum_real(bias)/iv%varbc_synop%u%nobs_sum

       if (rootproc) then
           write(unit=varbc_synop_unit,fmt='(10X,A,2F10.4)') &
                'bias corrected for synop u: (BC,OMB)', bc, bias
       end if

       ! v wind

       bias = 0.
       npred = iv%varbc_synop%v%npred

       if (iv%varbc_synop%v%nobs > 0) then 
           do iobs = 1, iv%varbc_synop%v%nobs
              bc = dot_product( iv%varbc_synop%v%param(1:npred),  &
                                iv%varbc_synop%v%pred (1:npred,iobs) )
              isn  = iv%varbc_synop%v%obs_sn(iobs)
              bias = bias + iv%synop(isn)%v%inv

              iv%synop(isn)%v%inv = iv%synop(isn)%v%inv - bc
           end do
       end if

       if (iv%varbc_synop%v%nobs_sum > 0) &
           bias = wrf_dm_sum_real(bias)/iv%varbc_synop%v%nobs_sum

       if (rootproc) then
           write(unit=varbc_synop_unit,fmt='(10X,A,2F10.4)') &
                'bias corrected for synop v: (BC,OMB)', bc, bias
       end if

       ! temperature

       bias = 0.
       npred = iv%varbc_synop%t%npred

       if (iv%varbc_synop%t%nobs > 0) then 
           do iobs = 1, iv%varbc_synop%t%nobs
              bc = dot_product( iv%varbc_synop%t%param(1:npred),  &
                                iv%varbc_synop%t%pred (1:npred,iobs) )
              isn  = iv%varbc_synop%t%obs_sn(iobs)
              bias = bias + iv%synop(isn)%t%inv

              iv%synop(isn)%t%inv = iv%synop(isn)%t%inv - bc
           end do
       end if

       if (iv%varbc_synop%t%nobs_sum > 0) &
           bias = wrf_dm_sum_real(bias)/iv%varbc_synop%t%nobs_sum

       if (rootproc) then
           write(unit=varbc_synop_unit,fmt='(10X,A,2F10.4)') &
                'bias corrected for synop t: (BC,OMB)', bc, bias
       end if

       ! qvapor

       bias = 0.
       npred = iv%varbc_synop%t%npred

       if (iv%varbc_synop%q%nobs > 0) then 
           do iobs = 1, iv%varbc_synop%q%nobs
              bc = dot_product( iv%varbc_synop%q%param(1:npred),  &
                                iv%varbc_synop%q%pred (1:npred,iobs) )
              isn  = iv%varbc_synop%q%obs_sn(iobs)
              bias = bias + iv%synop(isn)%q%inv

              iv%synop(isn)%q%inv = iv%synop(isn)%q%inv - bc
           end do
       end if

       if (iv%varbc_synop%q%nobs_sum > 0) &
           bias = wrf_dm_sum_real(bias)/iv%varbc_synop%q%nobs_sum

       if (rootproc) then
           write(unit=varbc_synop_unit,fmt='(10X,A,2F10.4)') &
                'bias corrected for synop q: (BC,OMB)', bc, bias
       end if

       if (rootproc) write(unit=varbc_synop_unit,fmt='(/10X,A)') &
                   "Variable Predictor/Parameter/Bias model (npred)  &  BC"

   if (trace_use) call da_trace_exit("da_varbc_synop_direct")

   end subroutine da_varbc_synop_direct
