   subroutine da_varbc_synop_init (iv)

   !-----------------------------------------------------!
   !  Initialize Varbc for SYNOP from VARBC_SYNOP.tbl  !
   !-----------------------------------------------------!

   implicit none

   type (iv_type), intent(inout)   :: iv

   integer, parameter              :: nmaxpred = 3  ! 1,Hobs-Hter,Tsk-T2 for t
                                                    ! 1,spd,ln(znt)      for u, v
                                                    ! 1,Hobs-Hter        for q

   integer                         :: npred, num_stn
   integer                         :: iunit,i

   logical                         :: ifexist_table
   character(len=120)              :: filename


   if (trace_use) call da_trace_entry("da_varbc_synop_init")
   
   if (rootproc) then
       write(unit=varbc_synop_unit,fmt='(//5X,A/)')      'VARBC_SYNOP namelist options'
       write(unit=varbc_synop_unit,fmt='(10X,A,7X,L)')   'use_varbc_synop :         ', use_varbc_synop
       write(unit=varbc_synop_unit,fmt='(10X,A,F8.3)')   'varbc_synop_pred0 :       ', varbc_synop_pred0    ! predictor(1) 
       write(unit=varbc_synop_unit,fmt='(10X,A,I8)')     'varbc_synop_nbgerr :      ', varbc_synop_nbgerr
       write(unit=varbc_synop_unit,fmt='(10X,A,I8)')     'varbc_synop_nobsmin :     ', varbc_synop_nobsmin
   end if

   filename = 'VARBC_SYNOP.tbl'
   inquire(file=trim(adjustl(filename)), exist=ifexist_table)

   if (ifexist_table) then

       if (rootproc) write(unit=varbc_synop_unit,fmt='(//5X,A/)') &
          'Reading in VARBC_SYNOP.tbl and initializating'

       iv%varbc_synop%u%nmaxpred = nmaxpred
       iv%varbc_synop%v%nmaxpred = nmaxpred
       iv%varbc_synop%t%nmaxpred = nmaxpred
       iv%varbc_synop%q%nmaxpred = nmaxpred

       call da_get_unit(iunit)
       open(unit=iunit,file=filename, status='old')
       do i = 1, 10; read(iunit,*); end do

       iv%varbc_synop%u%npred = 3
       iv%varbc_synop%v%npred = 3
       iv%varbc_synop%t%npred = 3
       iv%varbc_synop%q%npred = 2
       num_stn = iv%info(synop)%nlocal

       npred = iv%varbc_synop%u%npred
       allocate ( iv%varbc_synop%u%param    (nmaxpred) )
       allocate ( iv%varbc_synop%u%pred     (nmaxpred,num_stn) )
       allocate ( iv%varbc_synop%u%bgerr    (npred) )
       allocate ( iv%varbc_synop%u%index    (npred) )
       allocate ( iv%varbc_synop%u%vtox     (npred,npred) )

       npred = iv%varbc_synop%v%npred
       allocate ( iv%varbc_synop%v%param    (nmaxpred) )
       allocate ( iv%varbc_synop%v%pred     (nmaxpred,num_stn) )
       allocate ( iv%varbc_synop%v%bgerr    (npred) )
       allocate ( iv%varbc_synop%v%index    (npred) )
       allocate ( iv%varbc_synop%v%vtox     (npred,npred) )

       npred = iv%varbc_synop%t%npred
       allocate ( iv%varbc_synop%t%param    (nmaxpred) )
       allocate ( iv%varbc_synop%t%pred     (nmaxpred,num_stn) )
       allocate ( iv%varbc_synop%t%bgerr    (npred) )
       allocate ( iv%varbc_synop%t%index    (npred) )
       allocate ( iv%varbc_synop%t%vtox     (npred,npred) )

       npred = iv%varbc_synop%q%npred
       allocate ( iv%varbc_synop%q%param    (nmaxpred) )
       allocate ( iv%varbc_synop%q%pred     (nmaxpred,num_stn) )
       allocate ( iv%varbc_synop%q%bgerr    (npred) )
       allocate ( iv%varbc_synop%q%index    (npred) )
       allocate ( iv%varbc_synop%q%vtox     (npred,npred) )

       iv%varbc_synop%u%nobs      = 0
       iv%varbc_synop%u%nobs_sum  = 0
       iv%varbc_synop%u%param(:)  = 0.
       iv%varbc_synop%u%pred(:,:) = 0.
       iv%varbc_synop%u%pred(1,:) = varbc_synop_pred0
       iv%varbc_synop%u%bgerr(:)  = 0.
       iv%varbc_synop%u%index(:)  = 0
       iv%varbc_synop%u%vtox(:,:) = 0.

       iv%varbc_synop%v%nobs      = 0
       iv%varbc_synop%v%nobs_sum  = 0
       iv%varbc_synop%v%param(:)  = 0.
       iv%varbc_synop%v%pred(:,:) = 0.
       iv%varbc_synop%v%pred(1,:) = varbc_synop_pred0
       iv%varbc_synop%v%bgerr(:)  = 0.
       iv%varbc_synop%v%index(:)  = 0
       iv%varbc_synop%v%vtox(:,:) = 0.

       iv%varbc_synop%t%nobs      = 0
       iv%varbc_synop%t%nobs_sum  = 0
       iv%varbc_synop%t%param(:)  = 0.
       iv%varbc_synop%t%pred(:,:) = 0.
       iv%varbc_synop%t%pred(1,:) = varbc_synop_pred0
       iv%varbc_synop%t%bgerr(:)  = 0.
       iv%varbc_synop%t%index(:)  = 0
       iv%varbc_synop%t%vtox(:,:) = 0.

       iv%varbc_synop%q%nobs      = 0
       iv%varbc_synop%q%nobs_sum  = 0
       iv%varbc_synop%q%param(:)  = 0.
       iv%varbc_synop%q%pred(:,:) = 0.
       iv%varbc_synop%q%pred(1,:) = varbc_synop_pred0
       iv%varbc_synop%q%bgerr(:)  = 0.
       iv%varbc_synop%q%index(:)  = 0
       iv%varbc_synop%q%vtox(:,:) = 0.

       write(iv%varbc_synop%fmt_param,*) '(3(2X,E15.7))'

       !read u table
       read(iunit, *) !skip variable name
       read(iunit, fmt=trim(adjustl(iv%varbc_synop%fmt_param))) &
           iv%varbc_synop%u%param(1:nmaxpred)

       !read v table
       read(iunit, *) !skip variable name
       read(iunit, fmt=trim(adjustl(iv%varbc_synop%fmt_param))) &
           iv%varbc_synop%v%param(1:nmaxpred)

       !read t table
       read(iunit, *) !skip variable name
       read(iunit, fmt=trim(adjustl(iv%varbc_synop%fmt_param))) &
           iv%varbc_synop%t%param(1:nmaxpred)

       !read q table
       read(iunit, *) !skip variable name
       read(iunit, fmt=trim(adjustl(iv%varbc_synop%fmt_param))) &
           iv%varbc_synop%q%param(1:nmaxpred)

       close(iunit)
       call da_free_unit(iunit)

   else

       use_varbc_synop = .false.
       if (rootproc) write(unit=varbc_synop_unit,fmt='(/5X,A/)') &
                    'Could not find VARBC_SYNOP.tbl file. VARBC_SYNOP is switched off.'

   end if

   if (trace_use) call da_trace_exit("da_varbc_synop_init")

   end subroutine da_varbc_synop_init 
