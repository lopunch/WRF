   subroutine da_varbc_synop_pred(grid, iv, ob)

   !--------------------------------------!
   !  Calculate Predictors in Bias Model  !
   !--------------------------------------!
   
   implicit none

   type (domain),  intent(in)    :: grid    ! model state
   type (iv_type), intent(inout) :: iv
   type (y_type),  intent(inout) :: ob

   integer                       :: iobs,ipred
   real                          :: p1,p2

   real, dimension(iv%info(synop)%n1:iv%info(synop)%n2) :: &
        model_u, model_v, model_znt, model_tsk, model_t, model_hsm

   character(len=40)             :: stringn


   if (trace_use) call da_trace_entry("da_varbc_synop_pred")

   if (rootproc) write(unit=varbc_synop_unit,fmt='(//5X,A/)') &
                      'Calculating predictors'


   call da_interp_lin_2d (grid%xb%u10,  iv%info(synop), 1, model_u(:))
   call da_interp_lin_2d (grid%xb%v10,  iv%info(synop), 1, model_v(:))
   call da_interp_lin_2d (grid%xb%t2,   iv%info(synop), 1, model_t(:))
   call da_interp_lin_2d (grid%xb%rough,iv%info(synop), 1, model_znt(:))
   call da_interp_lin_2d (grid%xb%tsk,  iv%info(synop), 1, model_tsk(:))
   call da_interp_lin_2d (grid%xb%terr, iv%info(synop), 1, model_hsm(:))

   iv%varbc_synop%u%nobs = 0
   iv%varbc_synop%v%nobs = 0
   iv%varbc_synop%t%nobs = 0
   iv%varbc_synop%q%nobs = 0

   do iobs = iv%info(synop)%n1, iv%info(synop)%n2
       !predictors for u,v wind

       p1 =  hypot(model_u(iobs), model_v(iobs))
       p2 =  log(model_znt(iobs))

       if (iv%synop(iobs)%u%qc >= 0) then

          iv%varbc_synop%u%pred(2,iobs) = p1
          iv%varbc_synop%u%pred(3,iobs) = p2
          if(iv%info(synop)%proc_domain(1,iobs)) iv%varbc_synop%u%nobs = iv%varbc_synop%u%nobs + 1

       end if

       if (iv%synop(iobs)%v%qc >= 0) then

          iv%varbc_synop%v%pred(2,iobs) = p1
          iv%varbc_synop%v%pred(3,iobs) = p2
          if(iv%info(synop)%proc_domain(1,iobs)) iv%varbc_synop%v%nobs = iv%varbc_synop%v%nobs + 1

       end if

       !predictors for temperature

       p1 = model_tsk(iobs) - model_t(iobs)
       p2 = iv%synop(iobs)%h  - model_hsm(iobs)

       if (iv%synop(iobs)%t%qc >= 0) then

          iv%varbc_synop%t%pred(2,iobs) = p1
          iv%varbc_synop%t%pred(3,iobs) = p2
          if(iv%info(synop)%proc_domain(1,iobs)) iv%varbc_synop%t%nobs = iv%varbc_synop%t%nobs + 1

       end if

       !predictors for qvapor

       p1 = iv%synop(iobs)%h  - model_hsm(iobs)

       if (iv%synop(iobs)%q%qc >= 0) then

          iv%varbc_synop%q%pred(2,iobs) = p1
          if(iv%info(synop)%proc_domain(1,iobs)) iv%varbc_synop%q%nobs = iv%varbc_synop%q%nobs + 1

       end if
   end do

   !====================================================================

   !u wind

   iv%varbc_synop%u%nobs_sum = wrf_dm_sum_integer(iv%varbc_synop%u%nobs)

   if (iv%varbc_synop%u%nobs_sum < varbc_synop_nobsmin) then

       do iobs = iv%info(synop)%n1, iv%info(synop)%n2
          if (iv%synop(iobs)%u%qc >= 0) iv%synop(iobs)%u%qc = fail_varbc_synop
       end do

   end if

   !v wind

   iv%varbc_synop%v%nobs_sum = wrf_dm_sum_integer(iv%varbc_synop%v%nobs)

   if (iv%varbc_synop%v%nobs_sum < varbc_synop_nobsmin) then

       do iobs = iv%info(synop)%n1, iv%info(synop)%n2
          if (iv%synop(iobs)%v%qc >= 0) iv%synop(iobs)%v%qc = fail_varbc_synop
       end do

   end if

   !temperature

   iv%varbc_synop%t%nobs_sum = wrf_dm_sum_integer(iv%varbc_synop%t%nobs)

   if (iv%varbc_synop%t%nobs_sum < varbc_synop_nobsmin) then

       do iobs = iv%info(synop)%n1, iv%info(synop)%n2
          if (iv%synop(iobs)%t%qc >= 0) iv%synop(iobs)%t%qc = fail_varbc_synop
       end do

   end if

   !qvapor

   iv%varbc_synop%q%nobs_sum = wrf_dm_sum_integer(iv%varbc_synop%q%nobs)

   if (iv%varbc_synop%q%nobs_sum < varbc_synop_nobsmin) then

       do iobs = iv%info(synop)%n1, iv%info(synop)%n2
          if (iv%synop(iobs)%q%qc >= 0) iv%synop(iobs)%q%qc = fail_varbc_synop
       end do

   end if
  

   if (trace_use) call da_trace_exit("da_varbc_synop_pred")

   end subroutine da_varbc_synop_pred
