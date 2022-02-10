   subroutine da_varbc_synop_pred(grid, iv, be, ob)

   !--------------------------------------!
   !  Calculate Predictors in Bias Model  !
   !--------------------------------------!
   
   implicit none

   type (domain),  intent(in)    :: grid    ! model state
   type (iv_type), intent(inout) :: iv
   type (be_type), intent(inout) :: be
   type (y_type),  intent(inout) :: ob

   integer                       :: jsp_start,size_jsp
   integer                       :: iobs,ipred
   real                          :: p1,p2

   real, dimension(iv%info(synop)%n1:iv%info(synop)%n2) :: &
        model_u, model_v, model_znt, model_tsk, model_t, model_hsm

   character(len=40)             :: stringn


   if (trace_use) call da_trace_entry("da_varbc_synop_pred")

   if (rootproc) write(unit=varbc_synop_unit,fmt='(//5X,A/)') &
                      'Calculating predictors'

   size_jsp  = 0
   jsp_start = be%cv%size_jb + be%cv%size_je + be%cv%size_jp + be%cv%size_js + be%cv%size_jl + be%cv%size_jt

   call da_interp_lin_2d (grid%xb%u10,  iv%info(synop), 1, model_u(:))
   call da_interp_lin_2d (grid%xb%v10,  iv%info(synop), 1, model_v(:))
   call da_interp_lin_2d (grid%xb%t2,   iv%info(synop), 1, model_t(:))
   call da_interp_lin_2d (grid%xb%rough,iv%info(synop), 1, model_znt(:))
   call da_interp_lin_2d (grid%xb%tsk,  iv%info(synop), 1, model_tsk(:))
   call da_interp_lin_2d (grid%xb%terr, iv%info(synop), 1, model_hsm(:))

   do iobs = iv%info(synop)%n1, iv%info(synop)%n2
      if( iv%info(synop)%proc_domain(1,iobs) ) then !Don't count HALO data
         !predictors for u,v wind

         p1 =  hypot(model_u(iobs), model_v(iobs))
         p2 =  log(model_znt(iobs))

!        if (ob%synop(iobs)%u > missing_r .and. iv%synop(iobs)%u%qc >= 0) then
         if (iv%synop(iobs)%u%qc >= 0) then

            iv%varbc_synop%u%nobs = iv%varbc_synop%u%nobs + 1
            iv%varbc_synop%u%obs_sn(iv%varbc_synop%u%nobs) = iobs
            iv%varbc_synop%u%pred(2,iv%varbc_synop%u%nobs) = p1
            iv%varbc_synop%u%pred(3,iv%varbc_synop%u%nobs) = p2

         end if

!        if (ob%synop(iobs)%v > missing_r .and. iv%synop(iobs)%v%qc >= 0) then
         if (iv%synop(iobs)%v%qc >= 0) then

            iv%varbc_synop%v%nobs = iv%varbc_synop%v%nobs + 1
            iv%varbc_synop%v%obs_sn(iv%varbc_synop%v%nobs) = iobs
            iv%varbc_synop%v%pred(2,iv%varbc_synop%v%nobs) = p1
            iv%varbc_synop%v%pred(3,iv%varbc_synop%v%nobs) = p2

         end if

         !predictors for temperature

         p1 = model_tsk(iobs) - model_t(iobs)
         p2 = iv%synop(iobs)%h  - model_hsm(iobs)

!        if (ob%synop(iobs)%t > 0.0 .and. iv%synop(iobs)%t%qc >= 0) then
         if (iv%synop(iobs)%t%qc >= 0) then

            iv%varbc_synop%t%nobs = iv%varbc_synop%t%nobs + 1
            iv%varbc_synop%t%obs_sn(iv%varbc_synop%t%nobs) = iobs
            iv%varbc_synop%t%pred(2,iv%varbc_synop%t%nobs) = p1
            iv%varbc_synop%t%pred(3,iv%varbc_synop%t%nobs) = p2

         end if

         !predictors for qvapor

         p1 = iv%synop(iobs)%h  - model_hsm(iobs)

!        if (ob%synop(iobs)%q > 0.0 .and. iv%synop(iobs)%q%qc >= 0) then
         if (iv%synop(iobs)%q%qc >= 0) then

            iv%varbc_synop%q%nobs = iv%varbc_synop%q%nobs + 1
            iv%varbc_synop%q%obs_sn(iv%varbc_synop%q%nobs) = iobs
            iv%varbc_synop%q%pred(2,iv%varbc_synop%q%nobs) = p1

         end if
      end if  !end proc_domain
   end do

   !====================================================================

   !u wind

   iv%varbc_synop%u%nobs_sum = wrf_dm_sum_integer(iv%varbc_synop%u%nobs)

   if (iv%varbc_synop%u%nobs_sum >= varbc_synop_nobsmin) then

       if(rootproc) then
           do ipred = 1, iv%varbc_synop%u%npred 
              size_jsp = size_jsp + 1
              iv%varbc_synop%u%index(ipred) = jsp_start + size_jsp
              iv%varbc_synop%u%vtox(ipred,ipred) = 1.0
           end do
       end if

   else

       if (iv%varbc_synop%u%nobs > 0) then
           do iobs = iv%info(synop)%n1, iv%info(synop)%n2
              if (iv%synop(iobs)%u%qc >= 0) iv%synop(iobs)%u%qc = fail_varbc_synop
           end do
       end if

   end if

   !v wind

   iv%varbc_synop%v%nobs_sum = wrf_dm_sum_integer(iv%varbc_synop%v%nobs)

   if (iv%varbc_synop%v%nobs_sum >= varbc_synop_nobsmin) then

       if( rootproc ) then
           do ipred = 1, iv%varbc_synop%v%npred 
              size_jsp = size_jsp + 1
              iv%varbc_synop%v%index(ipred) = jsp_start + size_jsp
              iv%varbc_synop%v%vtox(ipred,ipred) = 1.0
           end do
       end if

   else

       if (iv%varbc_synop%v%nobs > 0) then
           do iobs = iv%info(synop)%n1, iv%info(synop)%n2
              if (iv%synop(iobs)%v%qc >= 0) iv%synop(iobs)%v%qc = fail_varbc_synop
           end do
       end if

   end if

   !temperature

   iv%varbc_synop%t%nobs_sum = wrf_dm_sum_integer(iv%varbc_synop%t%nobs)

   if (iv%varbc_synop%t%nobs_sum >= varbc_synop_nobsmin) then

       if( rootproc ) then
           do ipred = 1, iv%varbc_synop%t%npred 
              size_jsp = size_jsp + 1
              iv%varbc_synop%t%index(ipred) = jsp_start + size_jsp
              iv%varbc_synop%t%vtox(ipred,ipred) = 1.0
           end do
       end if

   else

       if (iv%varbc_synop%t%nobs > 0) then
           do iobs = iv%info(synop)%n1, iv%info(synop)%n2
              if (iv%synop(iobs)%t%qc >= 0) iv%synop(iobs)%t%qc = fail_varbc_synop
           end do
       end if

   end if

   !qvapor

   iv%varbc_synop%q%nobs_sum = wrf_dm_sum_integer(iv%varbc_synop%q%nobs)

   if (iv%varbc_synop%q%nobs_sum >= varbc_synop_nobsmin) then

       if( rootproc ) then
           do ipred = 1, iv%varbc_synop%q%npred 
              size_jsp = size_jsp + 1
              iv%varbc_synop%q%index(ipred) = jsp_start + size_jsp
              iv%varbc_synop%q%vtox(ipred,ipred) = 1.0
           end do
       end if

   else

       if (iv%varbc_synop%q%nobs > 0) then
           do iobs = iv%info(synop)%n1, iv%info(synop)%n2
              if (iv%synop(iobs)%q%qc >= 0) iv%synop(iobs)%q%qc = fail_varbc_synop
           end do
       end if

   end if
  
   if(rootproc) be % cv % size_jsp = size_jsp

!  if (rootproc .and. iv%varbc_synop%u%nobs_sum >= varbc_synop_nobsmin) &
!      write(unit=varbc_synop_unit,fmt='(/10X,A)') " U    Predictors  &  Obs Count"

!  write(stringn,'(I3)') iv%varbc_synop%u%npred
!  stringn = '(14X,'//trim(ADJUSTL(stringn))//'f8.3,3I5)'
!  stringn = trim(adjustl(stringn))

!  if (rootproc .and. iv%varbc_synop%u%nobs_sum >= varbc_synop_nobsmin) then
!      write(unit=varbc_synop_unit,fmt=stringn) &
!              iv%varbc_synop%u%pred(1:iv%varbc_synop%u%npred), &
!              iv%varbc_synop%u%nobs_sum
!  end if

!  if (rootproc .and. iv%varbc_synop%v%nobs_sum >= varbc_synop_nobsmin) &
!      write(unit=varbc_synop_unit,fmt='(/10X,A)') " V    Predictors  &  Obs Count"

!  write(stringn,'(I3)') iv%varbc_synop%v%npred
!  stringn = '(14X,'//trim(ADJUSTL(stringn))//'f8.3,3I5)'
!  stringn = trim(adjustl(stringn))

!  if (rootproc .and. iv%varbc_synop%v%nobs_sum >= varbc_synop_nobsmin) then
!      write(unit=varbc_synop_unit,fmt=stringn) &
!              iv%varbc_synop%v%pred(1:iv%varbc_synop%v%npred), &
!              iv%varbc_synop%v%nobs_sum
!  end if

!  if (rootproc .and. iv%varbc_synop%t%nobs_sum >= varbc_synop_nobsmin) &
!      write(unit=varbc_synop_unit,fmt='(/10X,A)') " T    Predictors  &  Obs Count"

!  write(stringn,'(I3)') iv%varbc_synop%t%npred
!  stringn = '(14X,'//trim(ADJUSTL(stringn))//'f8.3,3I5)'
!  stringn = trim(adjustl(stringn))

!  if (rootproc .and. iv%varbc_synop%t%nobs_sum >= varbc_synop_nobsmin) then
!      write(unit=varbc_synop_unit,fmt=stringn) &
!              iv%varbc_synop%t%pred(1:iv%varbc_synop%t%npred), &
!              iv%varbc_synop%t%nobs_sum
!  end if

!  if (rootproc .and. iv%varbc_synop%q%nobs_sum >= varbc_synop_nobsmin) &
!      write(unit=varbc_synop_unit,fmt='(/10X,A)') " Q    Predictors  &  Obs Count"

!  write(stringn,'(I3)') iv%varbc_synop%q%npred
!  stringn = '(14X,'//trim(ADJUSTL(stringn))//'f8.3,3I5)'
!  stringn = trim(adjustl(stringn))

!  if (rootproc .and. iv%varbc_synop%q%nobs_sum >= varbc_synop_nobsmin) then
!      write(unit=varbc_synop_unit,fmt=stringn) &
!              iv%varbc_synop%q%pred(1:iv%varbc_synop%q%npred), &
!              iv%varbc_synop%q%nobs_sum
!  end if


   if (trace_use) call da_trace_exit("da_varbc_synop_pred")

   end subroutine da_varbc_synop_pred
