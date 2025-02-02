   subroutine da_varbc_synop_update (it, cv_size, cv, iv)

   !-----------------------------------!
   !  Update and write out parameters  !
   !-----------------------------------!

   implicit none

   integer,        intent(in)      :: it        !outer loop counting
   integer,        intent(in)      :: cv_size
   real,           intent(in)      :: cv(1:cv_size)
   type (iv_type), intent(inout)   :: iv   

   integer                         :: i,ipred,npred
   integer                         :: iunit,iost
   character(len=filename_len)     :: filename
   character(len=99)               :: fmt_param
   character(len=30)               :: stringn

   real, allocatable               :: varbc_param_tl(:)


   if (trace_use) call da_trace_entry("da_varbc_synop_update")

   if (rootproc) write(unit=varbc_synop_unit,fmt='(//5X,A)') 'Updating parameters'

       !u wind 
      
       npred = iv%varbc_synop%u%npred 

       write(stringn,'(I3)') iv%varbc_synop%u%npred 
       stringn = '(2X,A,'//trim(ADJUSTL(stringn))//'f9.3)'
       stringn = trim(adjustl(stringn))

       allocate( varbc_param_tl(npred) )

       write(unit=varbc_synop_unit,fmt='(/10X,A)') &
             " ID  Phase  Parameters updated (npred)"

       varbc_param_tl(:) = 0.0

       if (iv%varbc_synop%u%nobs_sum >= varbc_synop_nobsmin) then

           if (rootproc) then
               do ipred = 1, npred
                  varbc_param_tl(ipred) = &
                            SUM( cv(iv%varbc_synop%u%index(1:npred)) * &
                                 iv%varbc_synop%u%vtox(ipred,1:npred) )
               end do

               iv%varbc_synop%u%param(1:npred) = iv%varbc_synop%u%param(1:npred) + varbc_param_tl    

               write(unit=varbc_synop_unit,fmt=stringn) 'U parameters updated : ',varbc_param_tl(:)
           end if

           call wrf_dm_bcast_real( iv%varbc_synop%u%param(1:npred), npred )

       end if

       deallocate(varbc_param_tl)

       !v wind 
      
       npred = iv%varbc_synop%v%npred 

       write(stringn,'(I3)') iv%varbc_synop%v%npred 
       stringn = '(2X,A,'//trim(ADJUSTL(stringn))//'f9.3)'
       stringn = trim(adjustl(stringn))

       allocate( varbc_param_tl(npred) )

       write(unit=varbc_synop_unit,fmt='(/10X,A)') &
             " ID  Phase  Parameters updated (npred)"

       varbc_param_tl(:) = 0.0

       if (iv%varbc_synop%v%nobs_sum >= varbc_synop_nobsmin) then

           if (rootproc) then
               do ipred = 1, npred
                  varbc_param_tl(ipred) = &
                            SUM( cv(iv%varbc_synop%v%index(1:npred)) * &
                                 iv%varbc_synop%v%vtox(ipred,1:npred) )
               end do

               iv%varbc_synop%v%param(1:npred) = iv%varbc_synop%v%param(1:npred) + varbc_param_tl    

               write(unit=varbc_synop_unit,fmt=stringn) 'V parameters updated : ',varbc_param_tl(:)
           end if

           call wrf_dm_bcast_real( iv%varbc_synop%v%param(1:npred), npred )

       end if

       deallocate(varbc_param_tl)

       !temperature
      
       npred = iv%varbc_synop%t%npred 

       write(stringn,'(I3)') iv%varbc_synop%t%npred 
       stringn = '(2X,A,'//trim(ADJUSTL(stringn))//'f9.3)'
       stringn = trim(adjustl(stringn))

       allocate( varbc_param_tl(npred) )

       write(unit=varbc_synop_unit,fmt='(/10X,A)') &
             " ID  Phase  Parameters updated (npred)"

       varbc_param_tl(:) = 0.0

       if (iv%varbc_synop%t%nobs_sum >= varbc_synop_nobsmin) then

           if (rootproc) then
               do ipred = 1, npred
                  varbc_param_tl(ipred) = &
                            SUM( cv(iv%varbc_synop%t%index(1:npred)) * &
                                 iv%varbc_synop%t%vtox(ipred,1:npred) )
               end do

               iv%varbc_synop%t%param(1:npred) = iv%varbc_synop%t%param(1:npred) + varbc_param_tl    

               write(unit=varbc_synop_unit,fmt=stringn) 'T parameters updated : ',varbc_param_tl(:)
           end if

           call wrf_dm_bcast_real( iv%varbc_synop%t%param(1:npred), npred )

       end if

       deallocate(varbc_param_tl)

       !qvapor
      
       npred = iv%varbc_synop%q%npred 

       write(stringn,'(I3)') iv%varbc_synop%q%npred 
       stringn = '(2X,A,'//trim(ADJUSTL(stringn))//'f9.3)'
       stringn = trim(adjustl(stringn))

       allocate( varbc_param_tl(npred) )

       write(unit=varbc_synop_unit,fmt='(/10X,A)') &
             " ID  Phase  Parameters updated (npred)"

       varbc_param_tl(:) = 0.0

       if (iv%varbc_synop%q%nobs_sum >= varbc_synop_nobsmin) then

           if (rootproc) then
               do ipred = 1, npred
                  varbc_param_tl(ipred) = &
                            SUM( cv(iv%varbc_synop%q%index(1:npred)) * &
                                 iv%varbc_synop%q%vtox(ipred,1:npred) )
               end do

               iv%varbc_synop%q%param(1:npred) = iv%varbc_synop%q%param(1:npred) + varbc_param_tl    

               write(unit=varbc_synop_unit,fmt=stringn) 'Q parameters updated : ',varbc_param_tl(:)
           end if

           call wrf_dm_bcast_real( iv%varbc_synop%q%param(1:npred), npred )

       end if

       deallocate(varbc_param_tl)

   if (.not. rootproc) then
       if (trace_use) call da_trace_exit("da_varbc_synop_update")
       return
   end if

   !-----------------------------------!
   !  Write VARBC_SYNOP.tbl.out file  !
   !-----------------------------------!
   
 ! write updated parameters

   call da_get_unit(iunit)

   if( it == max_ext_its ) then
      filename = 'VARBC_SYNOP.tbl.out'
   else
      write(unit=filename, fmt='(a,i2.2)') 'VARBC_SYNOP.tbl.out_',it
   end if

   open(unit=iunit,file=trim(adjustl(filename)),iostat=iost,status='replace')
   if (iost /= 0) then
      message(1)="Cannot open SYNOP bias correction file "//adjustl(filename)
      call da_error(__FILE__,__LINE__,message(1:1))
   end if

   write(iunit, '(A)')'*'
   write(iunit, '(A)')'  PARAMETER TABLE FOR SYNOP VARBC'
   write(iunit, '(A)')' '
   write(iunit, '(A)')'  Table format:'
   write(iunit, '(A)')'- Parameters 3(2X,E15.7)'
   write(iunit, '(A)')' '
   write(iunit, '(A)')'  Preditors(1:3): 1,Hobs-Hter,Tsk-T2 for t'
   write(iunit, '(A)')'  Preditors(1:3): 1,spd,ln(znt)      for u, v'
   write(iunit, '(A)')'  Preditors(1:2): 1,Hobs-Hter        for q'
   write(iunit, '(A)')'*'

   write(iunit,'(A)') '===== U ====='
   write(iunit,fmt=trim(adjustl(iv%varbc_synop%fmt_param))) &
         iv%varbc_synop%u%param(1:iv%varbc_synop%u%nmaxpred)

   write(iunit,'(A)') '===== V ====='
   write(iunit,fmt=trim(adjustl(iv%varbc_synop%fmt_param))) &
         iv%varbc_synop%v%param(1:iv%varbc_synop%v%nmaxpred)

   write(iunit,'(A)') '===== T ====='
   write(iunit,fmt=trim(adjustl(iv%varbc_synop%fmt_param))) &
         iv%varbc_synop%t%param(1:iv%varbc_synop%t%nmaxpred)

   write(iunit,'(A)') '===== Q ====='
   write(iunit,fmt=trim(adjustl(iv%varbc_synop%fmt_param))) &
         iv%varbc_synop%q%param(1:iv%varbc_synop%q%nmaxpred)

   close(iunit)
   call da_free_unit(iunit)

   if (rootproc) write(unit=varbc_synop_unit,fmt='(//5X,A/)') 'VARBC_SYNOP is done'

   if (trace_use) call da_trace_exit("da_varbc_synop_update")

   end subroutine da_varbc_synop_update
