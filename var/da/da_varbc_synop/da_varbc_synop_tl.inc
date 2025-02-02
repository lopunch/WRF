   subroutine da_varbc_synop_tl (cv_size, cv, iv, y)

   !---------------------------------------------------!
   !  METHOD: delta_d    = y - delta_bc                !
   !          y          = H (delta_x)                 !
   !          delta_bc   = SUM( delta_beta_i Pred_i )  !
   !---------------------------------------------------!

   implicit none

   integer, intent(in)           :: cv_size
   real, intent(in)              :: cv(1:cv_size)
   type (iv_type), intent(in)    :: iv   
   type (y_type), intent(inout)  :: y               ! y = h (xa)

   integer                       :: ipred,iobs
   integer                       :: npred
   real                          :: delta_bc

   real, allocatable             :: varbc_param_tl(:)
   

   if (trace_use) call da_trace_entry("da_varbc_synop_tl")

   !u wind

   npred = iv%varbc_synop%u%npred
   allocate( varbc_param_tl(npred) )

   if (iv%varbc_synop%u%nobs_sum >= varbc_synop_nobsmin) then

       varbc_param_tl(:) = 0.0

       if (rootproc) then
           do ipred = 1, npred
              varbc_param_tl(ipred) = &
                     SUM( cv(iv%varbc_synop%u%index(1:npred)) * &
                             iv%varbc_synop%u%vtox(ipred,1:npred) )
           end do
       end if

       call wrf_dm_bcast_real( varbc_param_tl, npred )

       do iobs = iv%info(synop)%n1, iv%info(synop)%n2
           if (iv%synop(iobs)%u%qc >= 0) then
              delta_bc = dot_product( varbc_param_tl(1:npred), &
                               iv%varbc_synop%u%pred(1:npred,iobs) )
              y%synop(iobs)%u = y%synop(iobs)%u + delta_bc
	   end if
       end do

   end if

   deallocate(varbc_param_tl)

   !v wind

   npred = iv%varbc_synop%v%npred
   allocate( varbc_param_tl(npred) )

   if (iv%varbc_synop%v%nobs_sum >= varbc_synop_nobsmin) then

       varbc_param_tl(:) = 0.0

       if (rootproc) then
           do ipred = 1, npred
              varbc_param_tl(ipred) = &
                     SUM( cv(iv%varbc_synop%v%index(1:npred)) * &
                             iv%varbc_synop%v%vtox(ipred,1:npred) )
           end do
       end if

       call wrf_dm_bcast_real( varbc_param_tl, npred )

       do iobs = iv%info(synop)%n1, iv%info(synop)%n2
           if (iv%synop(iobs)%v%qc >= 0) then
              delta_bc = dot_product( varbc_param_tl(1:npred), &
                               iv%varbc_synop%v%pred(1:npred,iobs) )
              y%synop(iobs)%v = y%synop(iobs)%v + delta_bc
	   end if
       end do

   end if

   deallocate(varbc_param_tl)


   !temperature

   npred = iv%varbc_synop%t%npred
   allocate( varbc_param_tl(npred) )

   if (iv%varbc_synop%t%nobs_sum >= varbc_synop_nobsmin) then

       varbc_param_tl(:) = 0.0

       if (rootproc) then
           do ipred = 1, npred
              varbc_param_tl(ipred) = &
                     SUM( cv(iv%varbc_synop%t%index(1:npred)) * &
                             iv%varbc_synop%t%vtox(ipred,1:npred) )
           end do
       end if

       call wrf_dm_bcast_real( varbc_param_tl, npred )

       do iobs = iv%info(synop)%n1, iv%info(synop)%n2
           if (iv%synop(iobs)%t%qc >= 0) then
              delta_bc = dot_product( varbc_param_tl(1:npred), &
                               iv%varbc_synop%t%pred(1:npred,iobs) )
              y%synop(iobs)%t = y%synop(iobs)%t + delta_bc
	   end if
       end do

   end if

   deallocate(varbc_param_tl)

 
   !qvapor

   npred = iv%varbc_synop%q%npred
   allocate( varbc_param_tl(npred) )

   if (iv%varbc_synop%q%nobs_sum >= varbc_synop_nobsmin) then

       varbc_param_tl(:) = 0.0

       if (rootproc) then
           do ipred = 1, npred
              varbc_param_tl(ipred) = &
                     SUM( cv(iv%varbc_synop%q%index(1:npred)) * &
                             iv%varbc_synop%q%vtox(ipred,1:npred) )
           end do
       end if

       call wrf_dm_bcast_real( varbc_param_tl, npred )

       do iobs = iv%info(synop)%n1, iv%info(synop)%n2
           if (iv%synop(iobs)%q%qc >= 0) then
              delta_bc = dot_product( varbc_param_tl(1:npred), &
                               iv%varbc_synop%q%pred(1:npred,iobs) )
              y%synop(iobs)%q = y%synop(iobs)%q + delta_bc
	   end if
       end do

   end if

   deallocate(varbc_param_tl)



   if (trace_use) call da_trace_exit("da_varbc_synop_tl")

   end subroutine da_varbc_synop_tl
