module da_varbc_synop

   !---------------------------------------------------------------------------
   ! Purpose: module for variational bias correction of SYNOP temperature. 
   !---------------------------------------------------------------------------
   use module_domain, only : domain
   use da_interpolation, only : da_interp_lin_2d

   use module_dm, only : wrf_dm_sum_real, wrf_dm_sum_integer
   use da_wrf_interfaces, only : wrf_dm_bcast_real

   use da_control, only : trace_use,missing_r, qc_varbc_bad, stdout, filename_len, &
      rootproc, varbc_nbgerr, ierr, comm, max_ext_its, &
      obs_qc_pointer, synop, missing_r,print_detail_obs, &
      use_varbc_synop, varbc_synop_nbgerr, varbc_synop_nobsmin, &
      varbc_synop_unit, varbc_synop_pred0, fail_varbc_synop
   use da_define_structures, only : iv_type, y_type, be_type, varbc_synop_type
   use da_reporting, only : da_error, message, da_warning, da_message
   use da_tools, only : da_eof_decomposition, da_diff_seconds
   use da_tools_serial, only : da_free_unit, da_get_unit
   use da_tracing, only : da_trace_entry, da_trace_exit, da_trace, &
      da_trace_int_sort
   
   implicit none

contains

#include "da_varbc_synop_init.inc"
#include "da_varbc_synop_pred.inc"
#include "da_varbc_synop_direct.inc"
#include "da_varbc_synop_precond.inc"
#include "da_varbc_synop_tl.inc"
#include "da_varbc_synop_adj.inc"
#include "da_varbc_synop_update.inc"

end module da_varbc_synop
