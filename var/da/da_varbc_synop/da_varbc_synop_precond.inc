   subroutine da_varbc_synop_precond (iv)

   !----------------------------------------------------!
   !  Estimate Hessian for Preconditioning              !
   !                                                    !
   !  Preconditioning = inverse square root of Hessian  !
   !                                                    !
   !  Hermite real matrix:                              !
   !                 A = Q*D*Q^T                        !
   !           A^(1/k) = Q*D^(1/K)*Q^T                  !
   !----------------------------------------------------!

   implicit none

   type (iv_type), intent (inout) :: iv

   integer                        :: i,j,iobs,isn,ipred,npred
   real                           :: lhessian,lbgerr,predi,predj,verr

   real, allocatable 		  :: hessian(:,:)
   real*8, allocatable 		  :: eignvec(:,:), eignval(:)
   

   if (trace_use) call da_trace_entry("da_varbc_synop_precond")

   if (rootproc) &
       write(unit=varbc_synop_unit,fmt='(//5X,A)') 'Estimating hessian for preconditioning'

   !u wind

   npred = iv%varbc_synop%u%npred

   allocate ( hessian(npred, npred) )
   allocate ( eignvec(npred, npred) )
   allocate ( eignval(npred)        )

   if (iv%varbc_synop%u%nobs_sum >= varbc_synop_nobsmin) then

       hessian(:,:) = 0.
       eignvec(:,:) = 0.
       eignval(:)   = 0.

       lbgerr = 0.0
       if (iv%varbc_synop%u%nobs > 0) then
           do iobs = 1, iv%varbc_synop%u%nobs
              isn  = iv%varbc_synop%u%obs_sn(iobs)
              verr = iv%synop(isn)%u%error
              lbgerr = lbgerr + verr * verr / varbc_synop_nbgerr
           end do
       end if
       lbgerr = wrf_dm_sum_real(lbgerr)

       do i = 1, npred
          do j = i, npred
             lhessian = 0.0	
   
             if (iv%varbc_synop%u%nobs > 0) then
                 do iobs = 1, iv%varbc_synop%u%nobs
                    predi = iv%varbc_synop%u%pred(i,iobs)
                    predj = iv%varbc_synop%u%pred(j,iobs)

                    isn   = iv%varbc_synop%u%obs_sn(iobs)
                    verr  = iv%synop(isn)%u%error

                    lhessian = lhessian + predi*predj/(verr*verr)
                 end do
             end if

             hessian(i,j) = wrf_dm_sum_real(lhessian)
             hessian(j,i) = hessian(i,j)  	       
          end do
   
          iv%varbc_synop%u%bgerr(i) = lbgerr/iv%varbc_synop%u%nobs_sum
          if (iv%varbc_synop%u%bgerr(i) > 0.) & 
              hessian(i,i)=hessian(i,i)+1.0/iv%varbc_synop%u%bgerr(i)
       end do   

       if (rootproc) then
           call da_eof_decomposition(npred, hessian(1:npred,1:npred), &
                                     eignvec(1:npred,1:npred),eignval(1:npred))

           if (ANY( eignval(1:npred) <= 0 )) then
               write(unit=stdout,fmt='(A,3F12.5)') &
                    'VARBC_TAMDAR: non-positive Hessian for U Eigenvalues =',eignval(1:npred) 

               do i = 1, npred
                  if (hessian(i,i) > 0) iv%varbc_synop%u%vtox(i,i) = 1.0/sqrt(hessian(i,i))
               end do
           else
               do i = 1, npred
               do j = i, npred
                  iv%varbc_synop%u%vtox(i,j) = SUM( eignvec(i,1:npred)* &
                                                   sqrt(1.0/eignval(1:npred)) * &
                                                   eignvec(j,1:npred) )
                  iv%varbc_synop%u%vtox(j,i) = iv%varbc_synop%u%vtox(i,j)
               end do
               end do
           end if
       end if
   end if

   deallocate(hessian, eignvec, eignval)   

   !v wind

   npred = iv%varbc_synop%v%npred

   allocate ( hessian(npred, npred) )
   allocate ( eignvec(npred, npred) )
   allocate ( eignval(npred)        )

   if (iv%varbc_synop%v%nobs_sum >= varbc_synop_nobsmin) then

       hessian(:,:) = 0.
       eignvec(:,:) = 0.
       eignval(:)   = 0.

       lbgerr = 0.0
       if (iv%varbc_synop%v%nobs > 0) then
           do iobs = 1, iv%varbc_synop%v%nobs
              isn  = iv%varbc_synop%v%obs_sn(iobs)
              verr = iv%synop(isn)%v%error
              lbgerr = lbgerr + verr * verr / varbc_synop_nbgerr
           end do
       end if
       lbgerr = wrf_dm_sum_real(lbgerr)

       do i = 1, npred
          do j = i, npred
             lhessian = 0.0	
   
             if (iv%varbc_synop%v%nobs > 0) then
                 do iobs = 1, iv%varbc_synop%v%nobs
                    predi = iv%varbc_synop%v%pred(i,iobs)
                    predj = iv%varbc_synop%v%pred(j,iobs)

                    isn   = iv%varbc_synop%v%obs_sn(iobs)
                    verr  = iv%synop(isn)%v%error

                    lhessian = lhessian + predi*predj/(verr*verr)
                 end do
             end if
             hessian(i,j) = wrf_dm_sum_real(lhessian)
             hessian(j,i) = hessian(i,j)  	       
          end do
   
          iv%varbc_synop%v%bgerr(i) = lbgerr/iv%varbc_synop%v%nobs_sum
          if (iv%varbc_synop%v%bgerr(i) > 0.) & 
              hessian(i,i)=hessian(i,i)+1.0/iv%varbc_synop%v%bgerr(i)
       end do   

       if (rootproc) then
           call da_eof_decomposition(npred, hessian(1:npred,1:npred), &
                                     eignvec(1:npred,1:npred),eignval(1:npred))

           if (ANY( eignval(1:npred) <= 0 )) then
               write(unit=stdout,fmt='(A,3F12.5)') &
                    'VARBC_TAMDAR: non-positive Hessian for U Eigenvalues =',eignval(1:npred) 

               do i = 1, npred
                  if (hessian(i,i) > 0) iv%varbc_synop%v%vtox(i,i) = 1.0/sqrt(hessian(i,i))
               end do
           else
               do i = 1, npred
               do j = i, npred
                  iv%varbc_synop%v%vtox(i,j) = SUM( eignvec(i,1:npred)* &
                                                   sqrt(1.0/eignval(1:npred)) * &
                                                   eignvec(j,1:npred) )
                  iv%varbc_synop%v%vtox(j,i) = iv%varbc_synop%v%vtox(i,j)
               end do
               end do
           end if
       end if
   end if

   deallocate(hessian, eignvec, eignval)   

   !temperature

   npred = iv%varbc_synop%t%npred

   allocate ( hessian(npred, npred) )
   allocate ( eignvec(npred, npred) )
   allocate ( eignval(npred)        )

   if (iv%varbc_synop%t%nobs_sum >= varbc_synop_nobsmin) then

       hessian(:,:) = 0.
       eignvec(:,:) = 0.
       eignval(:)   = 0.

       lbgerr = 0.0
       if (iv%varbc_synop%t%nobs > 0) then
           do iobs = 1, iv%varbc_synop%t%nobs
              isn  = iv%varbc_synop%t%obs_sn(iobs)
              verr = iv%synop(isn)%t%error
              lbgerr = lbgerr + verr * verr / varbc_synop_nbgerr
           end do
       end if
       lbgerr = wrf_dm_sum_real(lbgerr)

       do i = 1, npred
          do j = i, npred
             lhessian = 0.0	
   
             if (iv%varbc_synop%t%nobs > 0) then
                 do iobs = 1, iv%varbc_synop%t%nobs
                    predi = iv%varbc_synop%t%pred(i,iobs)
                    predj = iv%varbc_synop%t%pred(j,iobs)

                    isn   = iv%varbc_synop%t%obs_sn(iobs)
                    verr  = iv%synop(isn)%t%error

                    lhessian = lhessian + predi*predj/(verr*verr)
                 end do
             end if
             hessian(i,j) = wrf_dm_sum_real(lhessian)
             hessian(j,i) = hessian(i,j)  	       
          end do
   
          iv%varbc_synop%t%bgerr(i) = lbgerr/iv%varbc_synop%t%nobs_sum
          if (iv%varbc_synop%t%bgerr(i) > 0.) & 
              hessian(i,i)=hessian(i,i)+1.0/iv%varbc_synop%t%bgerr(i)
       end do   

       if (rootproc) then
           call da_eof_decomposition(npred, hessian(1:npred,1:npred), &
                                     eignvec(1:npred,1:npred),eignval(1:npred))

           if (ANY( eignval(1:npred) <= 0 )) then
               write(unit=stdout,fmt='(A,3F12.5)') &
                    'VARBC_TAMDAR: non-positive Hessian for U Eigenvalues =',eignval(1:npred) 

               do i = 1, npred
                  if (hessian(i,i) > 0) iv%varbc_synop%t%vtox(i,i) = 1.0/sqrt(hessian(i,i))
               end do
           else
               do i = 1, npred
               do j = i, npred
                  iv%varbc_synop%t%vtox(i,j) = SUM( eignvec(i,1:npred)* &
                                                   sqrt(1.0/eignval(1:npred)) * &
                                                   eignvec(j,1:npred) )
                  iv%varbc_synop%t%vtox(j,i) = iv%varbc_synop%t%vtox(i,j)
               end do
               end do
           end if
       end if
   end if

   deallocate(hessian, eignvec, eignval)   

   !qvapor

   npred = iv%varbc_synop%q%npred

   allocate ( hessian(npred, npred) )
   allocate ( eignvec(npred, npred) )
   allocate ( eignval(npred)        )

   if (iv%varbc_synop%q%nobs_sum >= varbc_synop_nobsmin) then

       hessian(:,:) = 0.
       eignvec(:,:) = 0.
       eignval(:)   = 0.

       lbgerr = 0.0
       if (iv%varbc_synop%q%nobs > 0) then
           do iobs = 1, iv%varbc_synop%q%nobs
              isn  = iv%varbc_synop%q%obs_sn(iobs)
              verr = iv%synop(isn)%q%error
              lbgerr = lbgerr + verr * verr / varbc_synop_nbgerr
           end do
       end if
       lbgerr = wrf_dm_sum_real(lbgerr)

       do i = 1, npred
          do j = i, npred
             lhessian = 0.0	
   
             if (iv%varbc_synop%q%nobs > 0) then
                 do iobs = 1, iv%varbc_synop%q%nobs
                    predi = iv%varbc_synop%q%pred(i,iobs)
                    predj = iv%varbc_synop%q%pred(j,iobs)

                    isn   = iv%varbc_synop%q%obs_sn(iobs)
                    verr  = iv%synop(isn)%q%error

                    lhessian = lhessian + predi*predj/(verr*verr)
                 end do
             end if
             hessian(i,j) = wrf_dm_sum_real(lhessian)
             hessian(j,i) = hessian(i,j)  	       
          end do
   
          iv%varbc_synop%q%bgerr(i) = lbgerr/iv%varbc_synop%q%nobs_sum
          if (iv%varbc_synop%q%bgerr(i) > 0.) & 
              hessian(i,i)=hessian(i,i)+1.0/iv%varbc_synop%q%bgerr(i)
       end do   

       if (rootproc) then
           call da_eof_decomposition(npred, hessian(1:npred,1:npred), &
                                     eignvec(1:npred,1:npred),eignval(1:npred))

           if (ANY( eignval(1:npred) <= 0 )) then
               write(unit=stdout,fmt='(A,3F12.5)') &
                    'VARBC_TAMDAR: non-positive Hessian for U Eigenvalues =',eignval(1:npred) 

               do i = 1, npred
                  if (hessian(i,i) > 0) iv%varbc_synop%q%vtox(i,i) = 1.0/sqrt(hessian(i,i))
               end do
           else
               do i = 1, npred
               do j = i, npred
                  iv%varbc_synop%q%vtox(i,j) = SUM( eignvec(i,1:npred)* &
                                                   sqrt(1.0/eignval(1:npred)) * &
                                                   eignvec(j,1:npred) )
                  iv%varbc_synop%q%vtox(j,i) = iv%varbc_synop%q%vtox(i,j)
               end do
               end do
           end if
       end if
   end if

   deallocate(hessian, eignvec, eignval)   


   if (rootproc) &
       write(unit=varbc_synop_unit,fmt='(/5X,A)') 'End estimating hessian for preconditioning'

   if (trace_use) call da_trace_exit("da_varbc_synop_precond")

   end subroutine da_varbc_synop_precond
