subroutine da_sfc_wtq_adj (psfc, tg, ps, ts, qs_spechum, us, vs, regime,& 
   psfc_prime, tg_prime, ps_prime, ts_prime, qs_spechum_prime, &
   us_prime, vs_prime,  hs, roughness, xland, dx,      & 
   u10_prime, v10_prime, t2_prime, q2_prime,           &
   has_lsm, qsfc_wrf, ust_wrf, mol_wrf, hfx, qfx, pblh) 

   !---------------------------------------------------------------------------
   ! Purpose: Calculate the  10m wind, 2m temperature and moisture based on the
   ! similarity theory
   !
   ! Reference:
   ! ---------
   !
   !  input Variables(basic state):
   ! 
   !   psfc, tg               : surface pressure and ground temperature
   !   ps, ts, qs, us, vs, hs : model variable at lowlest half sigma level
   !   regime                 : PBL regime
   !
   !  input Variables(pertubation):
   ! 
   !   psfc_prime, tg_prime   : Surface pressure and ground temperature
   !   ps_prime, ts_prime,    : Model variables at the lowest half sigma
   !   qs_prime, us_prime,    : level 
   !   vs_prime               : 
   !
   !  Constants:
   !
   !   hs                     : height at the lowest half sigma level
   !   roughness              : roughness
   !   xland                  : land-water-mask (=2 water, =1 land)
   !
   !  output Variables(pertubation):
   !
   !   u10_prime, v10_prime   : 10-m high observed wind components
   !   t2_prime , q2_prime    : 2-m high observed temperature and mixing ratio
   !
   !---------------------------------------------------------------------------
   !  
   !                      psim  : mechanical psi at lowlest sigma level
   !                      psim2 : mechanical psi at 2m 
   !                      psimz : mechanical psi at 10m 
   !
   !---------------------------------------------------------------------------

   implicit none

   real, intent (in)          :: regime
   real, intent (in)          :: ps , ts , qs_spechum , us, vs, psfc, tg
   real, intent (inout)       :: ps_prime, ts_prime, qs_spechum_prime, us_prime, vs_prime, psfc_prime, tg_prime
   real, intent (in)          :: hs, roughness, xland
   real, value                :: u10_prime, v10_prime, t2_prime, q2_prime
   logical, intent(in), optional :: has_lsm
   real,    intent(in), optional :: qsfc_wrf, ust_wrf, mol_wrf
   real,    intent(in), optional :: hfx, qfx, pblh
   logical :: use_ust_wrf = .false.
   logical :: vconv_wrf


   ! Maximum number of iterations in computing psim, psih

   integer, parameter :: k_iteration = 10 
   !      integer, parameter :: k_iteration = 1

   ! h10 is the height of 10m where the wind observed
   ! h2  is the height of 2m where the temperature and 
   !        moisture observed.

   real, parameter :: h10 = 10.0, h2 = 2.0
   !
   ! Default roughness over the land

   real, parameter :: zint0 = 0.01 
   !
   ! Von Karman constant

   real, parameter :: k_kar = 0.4
   !
   ! Working variables

   logical :: use_table = .true.
   real :: Vc2, Va2, V2 
   real :: rib, rcp, xx, yy, cc, Pi
   real :: psiw, psiz, mol, ust, hol, holz, hol2
   real :: psim, psimz, psim2, psih, psihz, psih2
   real :: psit, psit2, psiq, psiq2
   real :: gzsoz0, gz10oz0, gz2oz0
   real :: eg, qg, tvg, tvs
   real :: ths, thg, thvs, thvg, vsgd, vsgd2, dx
   real :: zq0, z0
   real :: ust_s, hol_s, holz_s, hol2_s, psim_s, psim2_s, psimz_s, psih_s, psih2_s, psihz_s
   real :: rib_s, fluxc_s
   real :: rr, rz, r2, rhox, fluxc, vc, cqs2, chs2, rho
   integer :: nn, nz, n2

   real :: Vc2_prime, Va2_prime, V2_prime
   real :: rib_prime, xx_prime, yy_prime
   real :: psiw_prime, psiz_prime, mol_prime, ust_prime, &
           hol_prime, holz_prime, hol2_prime
   real :: psim_prime, psimz_prime, psim2_prime, &
           psih_prime, psihz_prime, psih2_prime
   real :: psit_prime, psit2_prime, &
           psiq_prime, psiq2_prime
   real :: qg_prime, qs_prime, tvg_prime, tvs_prime
   real :: ths_prime, thg_prime, thvs_prime, thvg_prime

   !yclo add
   real :: wspd, visc, restar, z0t, z0q
   real :: qg_s, wspd_s, z0t_s, psiq_s, psit_s, psit_ss, psit2_s, psiq2_s
   real :: wspd_prime, visc_prime, restar_prime, z0t_prime, z0q_prime
   real :: rr_prime, rz_prime, r2_prime
   real :: rhox_prime, fluxc_prime, rho_prime, cqs2_prime, chs2_prime
   real :: qs

   real, parameter :: ka = 2.4E-5

   integer :: iregime

   if (trace_use) call da_trace_entry("da_sfc_wtq_adj")

   !-----------------------------------------------------------------------------!
   !  initialize perturbation value
   !------- ----------------------------------------------------------------------!
   !        tg_prime = 0.0
   !        us_prime = 0.0
   !        vs_prime = 0.0
   !        ts_prime = 0.0
   !        ps_prime = 0.0
            qs_prime = 0.0
   !      psfc_prime = 0.0

   psim_prime = 0.0
   psimz_prime = 0.0
   psim2_prime = 0.0
   psih2_prime = 0.0
   psihz_prime = 0.0
   psih_prime = 0.0

   psiw_prime = 0.0
   psiz_prime = 0.0
   psit_prime = 0.0
   psit2_prime = 0.0
   psiq_prime = 0.0
   psiq2_prime = 0.0

   qg_prime = 0.0
   ths_prime = 0.0
   thg_prime = 0.0

   thvs_prime = 0.0
   thvg_prime = 0.0
   V2_prime = 0.0
   rib_prime = 0.0
   ust_prime = 0.0
   tvg_prime = 0.0
   tvs_prime = 0.0
   va2_prime = 0.0
   vc2_prime = 0.0

   wspd_prime = 0.0
   visc_prime = 0.0
   restar_prime = 0.0
   z0t_prime = 0.0
   z0q_prime = 0.0
   rr_prime = 0.0
   rz_prime = 0.0
   r2_prime = 0.0
   rhox_prime = 0.0
   fluxc_prime = 0.0
   rho_prime = 0.0
   cqs2_prime = 0.0
   chs2_prime = 0.0

   !  +++++++++++++++++++++++++++++++++ 
   ! A.0  Calculate Basic state
   !  +++++++++++++++++++++++++++++++++ 
   rcp = gas_constant/cp

   ! 1 Compute the roughness length based upon season and land use 
   ! =====================================

   ! 1.1 Define the rouhness length
   !     -----------------

   z0 = roughness

   if (z0 < 0.0001) z0 = 0.0001

   ! 1.2 Define the rouhgness length for moisture
   !     -----------------

   if (xland .ge. 1.5) then
      zq0 = z0
   else
      zq0 =  zint0
   end if

   ! 1.3 Define the some constant variable for psi
   !     -----------------

   gzsoz0 = log(hs/z0)

   gz10oz0 = log(h10/z0)

   gz2oz0 = log(h2/z0)


   ! 2.0 Calculate the virtual temperature
   ! =====================================

   ! 2.0 Convert specific humidity to mixing ratio
   !     ---------------------------------------------------------

   qs   = qs_spechum / (1.0 - qs_spechum)

   ! 2.1 Compute Virtual temperature on the lowest half sigma level
   !     ---------------------------------------------------------

   tvs  = ts * (1.0 + 0.608 * qs)

   ! 2.2 Convert ground virtual temperature assuming it's saturated
   !     ----------------------------------------------------------
   call da_tp_to_qs(tg, psfc, eg, qg)
   if ( present(qsfc_wrf) ) then
      if ( qsfc_wrf > 0.0 ) then
         qg = qsfc_wrf
      end if
   endif
   qg_s = qg
   qg   = qg/(1.0-qg)                 !hcl convert to mixing ratio

   tvg  = tg * (1.0 + 0.608 * qg)

   ! 3.0  Compute the potential temperature and virtual potential temperature
   ! =======================================================================

   ! 3.1 Potential temperature on the lowest half sigma level
   !     ----------------------------------------------------

   Pi = (100000.0 / ps) ** rcp
   ths  = ts * Pi

   ! 3.2 Virtual potential temperature on the lowest half sigma level
   !     ------------------------------------------------------------

   thvs = tvs * Pi

   ! 3.3 Potential temperature at the ground
   !     -----------------------------------

   Pi = (100000.0 / psfc) ** rcp
   thg  = tg * Pi

   ! 3.4 Virtual potential temperature at ground
   !     ---------------------------------------

   thvg = tvg * Pi

   ! 4.0  BULK RICHARDSON NUMBER AND MONI-OBUKOV LENGTH
   ! =================================================

   ! 4.1 Velocity
   !     --------
   
   !     Wind speed:

   Va2 =   us*us + vs*vs
     
   !     Convective velocity:

   vconv_wrf = .false.
   if ( present(hfx) .and. present(qfx) .and. present(pblh) ) then
      ! calculating vconv over land following wrf method
      if ( pblh > 0.0 ) then
         vconv_wrf = .true.
      end if
   end if


   if (thvg >= thvs) then
      Vc2 = thvg - thvs
   else
      Vc2 = 0.0
   end if
   if ( xland < 1.5 ) then !land
      if ( vconv_wrf ) then
         ! following the calculation as in module_sf_sfclay.F
         rhox = psfc/(gas_constant*tvg)

         fluxc_s = hfx*(thg-ths)/rhox/cp+0.608*tvg*qfx*(qg-qs)/rhox
         if(fluxc_s >= 0.) then
            fluxc = fluxc_s 
         else
            fluxc = 0.
         end if

         vc = (gravity/tg*pblh*fluxc)**0.33
         vc2 = vc*vc
      end if
   end if



   ! Calculate Mahrt and Sun low-res correction         ! Add by Eric Chiang ( AUG 2010 )
   vsgd = 0.32 * (max(dx/5000.-1.,0.))**0.33            ! Add by Eric Chiang ( AUG 2010 )
   vsgd2 = vsgd * vsgd                                  ! Add by Eric Chiang ( AUG 2010 )
   
   V2  = Va2 + Vc2 + vsgd2                              ! Modified by Eric Chiang ( AUG 2010 ) 

   !yclo add
   wspd = sqrt(v2)
   wspd_s = wspd
   wspd = max(wspd,0.1)
   v2 = wspd*wspd
   !yclo end add

   ! 4.2 Bulk richardson number
   !     ----------------------

   rib = (gravity * hs / ths) * (thvs - thvg) / V2
   rib_s = rib
   ! 5.0  CALCULATE PSI BASED UPON REGIME

   ! if previously unstable, do not let into regime 1 and 2
   if ( present(mol_wrf) ) then
      if ( mol_wrf < 0.0 ) then
        !rib = min(rib, 0.0)
        if(rib_s > 0.) then
            rib = 0.0
        end if
      end if
   end if

   ! 5.4.1  Calculate   ust, m/L (mol), h/L (hol)
   !        --------------------------

   psim = 0.0
   psih = 0.0

   ! Friction speed

   if ( present(ust_wrf) ) then
      if ( ust_wrf > 0.0 ) then
         use_ust_wrf = .true.
         ust = ust_wrf
      end if
   end if

   if(.not. use_ust_wrf ) then
      ust = k_kar * sqrt(v2) /(gzsoz0 - psim)
   end if

   ! save ust for adjoint:
   ust_s = ust

   ! Heat flux factor

   if ( present(mol_wrf) ) then
      mol = mol_wrf
   else
      mol = k_kar * (ths - thg)/(gzsoz0 - psih)
   end if

   ! =======================================

   iregime = int(regime)
   select case (iregime) 

   ! 5.1 Stable conditions (REGIME 1)
   !     ---------------------------

   case (1);

      psim = -10.0*gzsoz0
      psim_s  = psim
      psim = max(psim,-10.0)
      psimz = h10/hs * psim
      psimz_s = psimz
      psimz = max(psimz,-10.0)
      psim2 = h2 /hs * psim
      psim2_s = psim2
      psim2 = max(psim2,-10.0)

      psih = psim
      psihz = psimz
      psih2 = psim2

   ! 5.2 Mechanically driven turbulence (REGIME 2)
   !     ------------------------------------------

   case (2);

      psim = (-5.0 * rib) * gzsoz0 / (1.1 - 5.0*rib)
      psim_s  = psim
      if (psim < -10.0) then
         psim = -10.0
      end if

      psimz = h10/hs * psim
      psimz_s = psimz
      if (psimz < -10.0) then
         psimz = -10.0
      end if

      psim2 = h2/hs * psim
      psim2_s = psim2
      if (psim2 < -10.0) then
         psim2 = -10.0
      end if

      psih = psim
      psihz = psimz
      psih2 = psim2

   ! 5.3 Unstable Forced convection (REGIME 3)
   !     -------------------------------------

   case (3);

      psim = 0.0
      psimz = 0.0
      psim2 = 0.0
      psih = psim
      psihz = psimz
      psih2 = psim2


   ! 5.4 Free convection (REGIME 4)
   !     --------------------------

   case (4);

      ! Calculate psi m and pshi h using iteration method

      cc = 2.0 * atan(1.0)
      !
      !        do k = 1 , k_iteration


      ! Ratio of PBL height to Monin-Obukhov length

      if (ust .LT. 0.01) then
         hol = rib * gzsoz0
      else
         hol = k_kar * gravity * hs * mol / (ths * ust * ust)
      end if

      ! 5.4.2  Calculate n, nz, R, Rz
      !        --------------------------

      holz = (h10 / hs) * hol
      hol2 = (h2 / hs) * hol

      hol_s  = hol
      holz_s = holz
      hol2_s = hol2

      if (hol > 0.0) then
         hol = 0.0
      end if
      if (hol < -9.9999) then
         hol = -9.9999
      end if

      if (holz > 0.0) then
         holz = 0.0
      end if
      if (holz < -9.9999) then
         holz = -9.9999
      end if

      if (hol2 > 0.0) then
         hol2 = 0.0
      end if
      if (hol2 < -9.9999) then
         hol2 = -9.9999
      end if

      ! 5.4.3 Calculate Psim & psih
      !        --------------------------

      if ( use_table ) then
         ! Using the look-up table:
         nn = int(-100.0 * hol)
         rr = (-100.0 * hol) - nn
         psim = psimtb(nn) + rr * (psimtb(nn+1) - psimtb(nn))
         psih = psihtb(nn) + rr * (psihtb(nn+1) - psihtb(nn))
      else
         ! Using the continuous function:
         xx = (1.0 - 16.0 * hol) ** 0.25
         yy = log((1.0+xx*xx)/2.0)
         psim = 2.0 * log((1.0+xx)/2.0) + yy - 2.0 * atan(xx) + cc
         psih = 2.0 * yy
      end if

      if ( use_table ) then
         ! Using the look-up table:
         nz = int(-100.0 * holz)
         rz = (-100.0 * holz) - nz
         psimz = psimtb(nz) + rz * (psimtb(nz+1) - psimtb(nz))
         psihz = psihtb(nz) + rz * (psihtb(nz+1) - psihtb(nz))
      else
         ! Using the continuous function:
         xx = (1.0 - 16.0 * holz) ** 0.25
         yy = log((1.0+xx*xx)/2.0)
         psimz = 2.0 * log((1.0+xx)/2.0) + yy - 2.0 * atan(xx) + cc
         psihz = 2.0 * yy
      end if

      if ( use_table ) then
         ! Using the look-up table:
         n2 = int(-100.0 * hol2)
         r2 = (-100.0 * hol2) - n2
         psim2 = psimtb(n2) + r2 * (psimtb(n2+1) - psimtb(n2))
         psih2 = psihtb(n2) + r2 * (psihtb(n2+1) - psihtb(n2))
      else
         ! Using the continuous function:
         xx = (1.0 - 16.0 * hol2) ** 0.25
         yy = log((1.0+xx*xx)/2.0)
         psim2 = 2.0 * log((1.0+xx)/2.0) + yy - 2.0 * atan(xx) + cc
         psih2 = 2.0 * yy
      end if

      ! end do 

      ! 5.4.4 Define the limit value for psim & psih
      !        --------------------------
 
      !  Save the values for adjoint:

      psim_s  = psim
      psimz_s = psimz
      psim2_s = psim2
      psih_s  = psih
      psihz_s = psihz
      psih2_s = psih2

      if (psim > 0.9*gzsoz0) then
         psim = 0.9*gzsoz0
      end if
      if (psimz > 0.9*gz10oz0) then
         psimz = 0.9*gz10oz0
      end if
      if (psim2 > 0.9*gz2oz0) then
         psim2 = 0.9*gz2oz0
      end if
      if (psih > 0.9*gzsoz0) then
         psih = 0.9*gzsoz0
      end if
      if (psihz > 0.9*gz10oz0) then
         psihz = 0.9*gz10oz0
      end if
      if (psih2 > 0.9*gz2oz0) then
         psih2 = 0.9*gz2oz0
      end if

   case default;
      write(unit=message(1),fmt='(A,I2,A)') &
         "Regime=",iregime," is invalid."
      call da_error(__FILE__,__LINE__,message(1:1))

   end select
   
   ! 6.0  CALCULATE PSI FOR WinD, TEMPERATURE AND MOISTURE
   ! =======================================

   psiw = gzsoz0 - psim
   psiz = gz10oz0 - psimz

   if(gzsoz0 - psih < 2.0) then
      psit = 2.0
   else
      psit = gzsoz0 - psih
   end if

   psit2 = gz2oz0 - psih2

   if ( .not. use_ust_wrf ) then
      ! re-calculate ust since psim is now available
      ust = k_kar * sqrt(v2) /(gzsoz0 - psim)
   end if
   psiq  = log(k_kar*ust*hs/ka + hs / zq0) - psih
   psiq2 = log(k_kar*ust*h2/ka + h2 / zq0) - psih2

   !V3.7, as in module_sf_sfclay.F
   if ( xland >= 1.5 ) then !water
      visc = (1.32+0.009*(ts-273.15))*1.e-5
      restar = ust*z0/visc
      z0t = (5.5e-5)*(restar**(-0.60))
      z0t_s = z0t
      z0t = min(z0t,1.0e-4)
      z0t = max(z0t,2.0e-9)
      z0q = z0t
      psiq  = log((hs+z0q)/z0q)-psih
      psit  = log((hs+z0t)/z0t)-psih
      psiq2 = log((2.+z0q)/z0q)-psih2
      psit2 = log((2.+z0t)/z0t)-psih2

      psiq_s  = psiq
      psit_s  = psit
      psiq2_s = psiq2
      psit2_s = psit2

      psiq  = max(psiq,  2.)
      psit  = max(psit,  2.)
      psiq2 = max(psiq2, 2.)
      psit2 = max(psit2, 2.)
   end if


   !  +++++++++++++++++++++++++++++++++ 
   !  B.0  Calculate adjoint solution
   !  +++++++++++++++++++++++++++++++++ 

   ! 7.0  CALCULATE 10M WinD, 2M TEMPERATURE AND MOISTURE
   ! =======================================

   if ( present(has_lsm) ) then
      if ( has_lsm ) then
         !cqs2: 2m surface exchange coefficient for moisture
         !chs2: 2m surface exchange coefficient for heat
         cqs2 = ust*k_kar/psiq2
         if (xland .ge. 1.5) then
            !water
            chs2 = ust*k_kar/psit2
         else
            !land
            chs2 = cqs2 !as in subroutine lsm in phys/module_sf_noahdrv.F
         end if

         !re-calculate T2/Q2 as in module_sf_sfcdiags.F
         rho  = psfc/(gas_constant*tg)

         if ( chs2 < 1.e-5 ) then
            tg_prime   = tg_prime + t2_prime
            t2_prime = 0.0
         else
            if ( present(hfx) ) then
               tg_prime   = tg_prime + t2_prime
               rho_prime  = hfx*(thg-ths)/(rho*cp*chs2)**2 * cp*chs2*t2_prime
               chs2_prime = hfx*(thg-ths)/(rho*cp*chs2)**2 * cp*rho*t2_prime
               thg_prime  = - hfx*t2_prime/(rho*cp*chs2)
               ths_prime  =   hfx*t2_prime/(rho*cp*chs2)
               t2_prime   = 0.0
            end if
         end if

         if ( cqs2 < 1.e-5 ) then
            qg_prime   = q2_prime
            q2_prime = 0.0
         else
            if ( present(qfx) ) then
               qg_prime  = q2_prime
               rho_prime = rho_prime + qfx*(qg-qs)/(rho*cqs2)**2 * cqs2*q2_prime
               cqs2_prime =            qfx*(qg-qs)/(rho*cqs2)**2 * rho*q2_prime
               qg_prime  = qg_prime - qfx*q2_prime/(rho*cqs2)
               qs_prime  =            qfx*q2_prime/(rho*cqs2)
               q2_prime  = 0.0
            end if
         end if

         psfc_prime = psfc_prime + rho_prime/(gas_constant*tg)
         tg_prime   = tg_prime - psfc/(gas_constant*tg**2)*rho_prime

         if( xland .ge. 1.5) then
            !water
            ust_prime = chs2_prime*k_kar/psit2
            psit2_prime = -ust/psit2**2 * chs2_prime*k_kar
         else
            !land
            cqs2_prime  = cqs2_prime + chs2_prime
         end if

         ust_prime   = ust_prime + cqs2_prime*k_kar/psiq2
         psiq2_prime = -cqs2_prime*ust/psiq2**2*k_kar
      end if
   end if

   qg_prime    = qg_prime + (1.0 - psiq2/psiq) * q2_prime
   qs_prime    = qs_prime + psiq2/psiq * q2_prime
   psiq2_prime = psiq2_prime + (qs -qg)/psiq * q2_prime
   psiq_prime  = -(qs -qg)*psiq2/(psiq*psiq) * q2_prime
   ! q2_prime = 0.0
   
   Pi = (psfc/100000.0)**rcp
   thg_prime   = thg_prime + (1.0 - psit2/psit) * Pi * t2_prime
   ths_prime   = ths_prime + (psit2/psit) * Pi * t2_prime
   psit2_prime = psit2_prime + (ths - thg)/psit * Pi * t2_prime
   psit_prime  = - (ths - thg)*psit2/(psit*psit) * Pi * t2_prime
   psfc_prime  = psfc_prime + Pi * rcp*(thg + (ths - thg)*psit2/psit)/psfc * t2_prime 
   ! t2_prime = 0.0
   
   Pi = psiz / psiw
   vs_prime   = vs_prime + Pi * v10_prime
   psiz_prime =   vs / psiz * Pi * v10_prime
   psiw_prime = - vs / psiw * Pi * v10_prime
   ! v10_prime = 0.0 
   
   us_prime   = us_prime + Pi * u10_prime
   psiz_prime =  psiz_prime + us / psiz * Pi * u10_prime
   psiw_prime =  psiw_prime - us / psiw * Pi * u10_prime
   ! u10_prime = 0.0
   
   if ( xland >= 1.5 ) then !water
      psiq  = psiq_s
      psit  = psit_s
      psiq2 = psiq2_s
      psit2 = psit2_s

      if(psiq < 2.) psiq_prime = 0.
      if(psit < 2.) psit_prime = 0.
      if(psiq2 < 2.) psiq2_prime = 0.
      if(psit2 < 2.) psit2_prime = 0.

      z0t_prime = -2. / (z0t * (2.+z0t)) * psit2_prime
      z0q_prime = -2. / (z0q * (2.+z0q)) * psiq2_prime
      psih2_prime = - psiq2_prime - psit2_prime

      z0t_prime = z0t_prime - hs / (z0t * (hs+z0t)) * psit_prime
      z0q_prime = z0q_prime - hs / (z0q * (hs+z0q)) * psiq_prime
      psih_prime = - psiq_prime - psit_prime

      z0t_prime = z0t_prime + z0q_prime
      z0t = z0t_s
      if(z0t > 1.0e-4) z0t_prime = 0.0
      if(z0t < 2.0e-9) z0t_prime = 0.0

      restar_prime = -0.6 * z0t / restar * z0t_prime
      ust_prime    = ust_prime + z0 / visc * restar_prime
      visc_prime   = - z0 * ust / (visc * visc) * restar_prime
      ts_prime     = ts_prime + 0.009*1.e-5*visc_prime

      psiq_prime   = 0.0
      psit_prime   = 0.0
      psiq2_prime  = 0.0
      psit2_prime  = 0.0
   end if




   ! 6.0  CALCULATE PSI FOR WinD, TEMPERATURE AND MOISTURE
   ! =======================================
   
   ! moisture part:
!  psih2_prime = - psiq2_prime
!  psih_prime  = - psiq_prime
!  psiq2_prime = psiq2_prime + psiq_prime
!  ust_prime   = k_kar*hs/(ka*(k_kar*ust*hs/ka + hs / zq0)) * psiq2_prime 
   psih2_prime = psih2_prime - psiq2_prime
   psih_prime  = psih_prime - psiq_prime
!  psiq2_prime = psiq2_prime + psiq_prime
   ust_prime   = ust_prime + k_kar*hs/(ka*(k_kar*ust*hs/ka + hs / zq0)) * psiq_prime &
                           + k_kar*h2/(ka*(k_kar*ust*h2/ka + h2 / zq0)) * psiq2_prime 
   
   if ( .not. use_ust_wrf ) then
      V2_prime   = 0.5/V2 * ust * ust_prime
      psim_prime = ust / (gzsoz0 - psim) * ust_prime
      ust_prime = 0.0
   end if

   ! temperature part:
   psih2_prime = psih2_prime - psit2_prime
   if(gzsoz0 - psih >= 2.0) psih_prime  = psih_prime  - psit_prime

   ! wind part:
   psimz_prime = psimz_prime - psiz_prime
   psim_prime  = psim_prime  - psiw_prime

   ! 5.0  CALCULATE PSI BASED UPON REGIME
   ! =======================================
   select case (iregime) 

   ! 5.1 Stable conditions (REGIME 1)
   !     ---------------------------

   case (1);

      psim2_prime = psim2_prime + psih2_prime
      psimz_prime = psimz_prime + psihz_prime
      psim_prime  = psim_prime  + psih_prime

      psim_prime  = psim_prime + h2/hs  * psim2_prime &
                               + h10/hs * psimz_prime
      psim_prime  = 0.0


   ! 5.2 Mechanically driven turbulence (REGIME 2)
   !      ------------------------------------------

   case (2);

     psim2_prime = psim2_prime + psih2_prime
     psimz_prime = psimz_prime + psihz_prime
     psim_prime  = psim_prime  + psih_prime

     psim  = psim_s
     psimz = psimz_s
     psim2 = psim2_S

     if (psim2 < -10.0) then
        psim2_prime = 0.0
     end if
     psim_prime = psim_prime + h2/hs * psim2_prime

     if (psimz < -10.0) then
        psimz_prime = 0.0
     end if
     psim_prime = psim_prime + h10/hs * psimz_prime

     if (psim < -10.0) then
        psim_prime = 0.0
     end if
     rib_prime = -0.22 * gzsoz0 / (0.22-rib)**2 * psim_prime

     ! 5.3 Unstable Forced convection (REGIME 3)
     !     -------------------------------------

   case (3);

      psim2_prime = psih2_prime
      psimz_prime = psihz_prime
      psim_prime  = psih_prime

      psim2_prime = 0.0
      psimz_prime = 0.0
      psim_prime  = 0.0

   ! 5.4 Free convection (REGIME 4)
   !     --------------------------

   case (4);

      ! 5.4.4 Define the limit value for psim & psih
      !        -------------------------------------

      ! Recover the values:

      psim = psim_s
      psim2 = psim2_s
      psimz = psimz_s
      psihz = psihz_s
      psih  = psih_s
      psih2 = psih2_s

      if (psih2 > 0.9*gz2oz0) then
         psih2_prime = 0.0
      end if
      if (psihz > 0.9*gz10oz0) then
         psihz_prime = 0.0
      end if
      if (psih > 0.9*gzsoz0) then
         psih_prime = 0.0
      end if
      if (psim2 > 0.9*gz2oz0) then
         psim2_prime = 0.0
      end if
      if (psimz > 0.9*gz10oz0) then
         psimz_prime = 0.0
      end if
      if (psim > 0.9*gzsoz0) then
         psim_prime = 0.0
      end if

      ! 5.4.3 Calculate Psim & psih
      !        --------------------------

      ! Recover ust:
      ust = ust_s

      if ( .not. use_table ) then
         xx = (1.0 - 16.0 * hol2) ** 0.25
         yy = log((1.0+xx*xx)/2.0)
         yy_prime = 2.0 * psih2_prime
         psih2_prime = 0.0
         
         xx_prime = 2.0*(1.0/(1.0+xx)-1.0/(1+xx*xx)) * psim2_prime
         yy_prime = yy_prime + psim2_prime
         psim2_prime = 0.0
         
         xx_prime = xx_prime + 2.0* xx/(1.0+xx*xx) * yy_prime
         yy_prime = 0.0
         hol2_prime = -4.0/((1.0 - 16.0 * hol2) ** 0.75) * xx_prime
         xx_prime = 0.0
      else
         n2 = int(-100.0 * hol2)
         r2_prime = psih2_prime * (psihtb(n2+1) - psihtb(n2))
         psih2_prime = 0.0

         r2_prime = r2_prime + psim2_prime * (psimtb(n2+1) - psimtb(n2))
         psim2_prime = 0.0

         hol2_prime = -100.0 * r2_prime
         r2_prime = 0.0
      end if

      if ( .not. use_table ) then
         xx = (1.0 - 16.0 * holz) ** 0.25
         yy = log((1.0+xx*xx)/2.0)
         yy_prime = 2.0 * psihz_prime
         psihz_prime = 0.0
         
         xx_prime = 2.0*(1.0/(1.0+xx)-1.0/(1+xx*xx)) * psimz_prime
         yy_prime = yy_prime + psimz_prime
         psimz_prime = 0.0
         
         xx_prime = xx_prime + 2.0* xx/(1.0+xx*xx) * yy_prime
         yy_prime = 0.0
         holz_prime = -4.0/((1.0 - 16.0 * holz) ** 0.75) * xx_prime
         xx_prime = 0.0
      else
         nz = int(-100.0 * holz)
         rz_prime = psihz_prime * (psihtb(nz+1) - psihtb(nz))
         psihz_prime = 0.0

         rz_prime = rz_prime + psimz_prime * (psimtb(nz+1) - psimtb(nz))
         psimz_prime = 0.0

         holz_prime = -100.0 * rz_prime
         rz_prime = 0.0
      end if

      if ( .not. use_table ) then
         xx = (1.0 - 16.0 * hol) ** 0.25
         yy = log((1.0+xx*xx)/2.0)
         yy_prime = 2.0 * psih_prime
         psih_prime = 0.0
         
         xx_prime = 2.0*(1.0/(1.0+xx)-1.0/(1+xx*xx)) * psim_prime
         yy_prime = yy_prime + psim_prime
         psim_prime = 0.0
         
         xx_prime = xx_prime + 2.0* xx/(1.0+xx*xx)*yy_prime
         yy_prime = 0.0
         hol_prime = -4.0/((1.0 - 16.0 * hol) ** 0.75)*xx_prime
         xx_prime = 0.0
      else
         nn = int(-100.0 * hol)
         rr_prime = psih_prime * (psihtb(nn+1) - psihtb(nn))
         psih_prime = 0.0

         rr_prime = rr_prime + psim_prime * (psimtb(nn+1) - psimtb(nn))
         psim_prime = 0.0

         hol_prime = -100.0 * rr_prime
         rr_prime = 0.0
      end if




      ! 5.4.2  Calculate n, nz, R, Rz
      !        --------------------------

      !    Recover the values:

      hol  = hol_s
      holz = holz_s
      hol2 = hol2_s

      if (hol2 < -9.9999) then
         hol2_prime = 0.0
      end if
      if (hol2 > 0.0) then       
         hol2_prime = 0.0
      end if
      if (holz < -9.9999) then
         holz_prime = 0.0
      end if
      if (holz > 0.0) then
         holz_prime = 0.0
      end if
      if (hol < -9.9999) then
         hol_prime = 0.0
      end if
      if (hol > 0.0) then
         hol_prime = 0.0
      end if

      hol_prime = hol_prime + (h2 / hs)  * hol2_prime &
                            + (h10 / hs) * holz_prime

      ! 5.4.1  Calculate   ust, m/L (mol), h/L (hol)
      !        --------------------------

      !       Ratio of PBL height to Monin-Obukhov length
      if (ust .LT. 0.01) then
         rib_prime = hol_prime * gzsoz0
         hol_prime = 0.0
      else
!        mol_prime =        hol / mol * hol_prime
         mol_prime = k_kar * gravity * hs * hol_prime / (ths * ust * ust)
         ths_prime = ths_prime - hol / ths * hol_prime
         ust_prime = ust_prime - 2.0 * hol / ust * hol_prime
         hol_prime = 0.0
      end if


      ! Calculate psi m and pshi h using iteration method

   case default;
      write(unit=message(1),fmt='(A,I2,A)') &
         "Regime=",iregime," is invalid."
      call da_error(__FILE__,__LINE__,message(1:1))

   end select

   psim = 0.0
   psih = 0.0

   ! Heat flux factor
   if ( .not. present(mol_wrf) ) then
      ths_prime  = ths_prime + mol/(ths - thg) * mol_prime
      thg_prime  = thg_prime - mol/(ths - thg) * mol_prime
      psih_prime = psih_prime + mol/(gzsoz0 - psih) * mol_prime
      mol_prime = 0.0
   end if

      ! Friction speed

   if(.not. use_ust_wrf ) then
      V2_prime   = V2_prime + 0.5 * ust/V2 * ust_prime 
      psim_prime = psim_prime + ust/(gzsoz0 - psim) * ust_prime 
      ust_prime = 0.0
   end if

   psim_prime = 0.0
   psih_prime = 0.0


   ! 4.0  BULK RICHARDSON NUMBER AND MONI-OBUKOV LENGTH
   ! =================================================

   ! 4.2 Bulk richardson number
   !     ----------------------

   ! if previously unstable, do not let into regime 1 and 2
   if ( present(mol_wrf) ) then
      if ( mol_wrf < 0.0 ) then
         if(rib_s > 0.) then
            rib_prime = 0.0
         end if
      end if
   end if

   Pi = gravity * hs / (ths*V2)
   ths_prime = ths_prime - Pi * (thvs-thvg)/ths * rib_prime
   V2_prime  = V2_prime  - Pi * (thvs-thvg)/V2 * rib_prime
   thvs_prime = thvs_prime + Pi * rib_prime
   thvg_prime = thvg_prime - Pi * rib_prime
   rib_prime = 0.0
   
   ! 4.1 Velocity
   !     --------

   !yclo add
   wspd_prime = 2. * wspd * V2_prime
   wspd = wspd_s
   if(wspd < 0.1) wspd_prime = 0.0

   V2_prime = wspd_prime / (2. * wspd)
   !yclo end add


   ! Convective velocity:

   Va2_prime = V2_prime
   Vc2_prime = V2_prime

   if ( xland < 1.5 ) then !land
      if ( vconv_wrf ) then
         ! following the calculation as in module_sf_sfclay.F
         if(vc == 0.0) then
            Vc2_prime = 0.0
         else
            tg_prime  = tg_prime + &
                         2.*vc*0.33/(gravity/tg*pblh*fluxc)**0.67 * &
                    (-(gravity/tg*pblh*fluxc)*Vc2_prime/tg)
            fluxc_prime = 2.*vc*0.33/(gravity/tg*pblh*fluxc)**0.67 * &
                    gravity/tg*pblh*Vc2_prime
            Vc2_prime = 0.0
         end if

         if(fluxc_s >= 0.) then
            rhox_prime = -fluxc*fluxc_prime/rhox
            tvg_prime  = 0.608*qfx*(qg-qs)*fluxc_prime/rhox
            thg_prime  = thg_prime + hfx*fluxc_prime/rhox/cp
            ths_prime  = ths_prime - hfx*fluxc_prime/rhox/cp
            qg_prime   = qg_prime  + 0.608*tvg*qfx*fluxc_prime/rhox
            qs_prime   = qs_prime  - 0.608*tvg*qfx*fluxc_prime/rhox
         end if

         psfc_prime =  psfc_prime + rhox_prime/(gas_constant*tvg)
         tvg_prime  = tvg_prime - psfc/(gas_constant*tvg**2) * rhox_prime
      end if
   end if




   if (thvg >= thvs) then
      thvg_prime = thvg_prime + Vc2_prime
      thvs_prime = thvs_prime - Vc2_prime
      Vc2_prime = 0.0
   else
      Vc2_prime = 0.0
   end if

   ! Wind speed:

   us_prime = us_prime + 2.0 *us * Va2_prime
   vs_prime = vs_prime + 2.0 *vs * Va2_prime
   Va2_prime = 0.0

   ! 3.0 Virtual potential temperature
   ! ==================================

   ! 3.4 Virtual potential temperature at ground
   !     ---------------------------------------

   Pi = (100000.0 / psfc) ** rcp
   tvg_prime = tvg_prime +  Pi * thvg_prime
   psfc_prime = psfc_prime - thvg_prime * rcp * tvg/psfc * Pi
   thvg_prime = 0.0

   ! 3.3 Potential temperature at the ground
   !     -----------------------------------

   tg_prime = tg_prime + Pi * thg_prime 
   psfc_prime = psfc_prime - thg_prime *rcp * tg/psfc * Pi
   thg_prime = 0.0

   ! 3.2 Virtual potential temperature on the lowest half sigma level
   !     ------------------------------------------------------------

   Pi = (100000.0 / ps) ** rcp
   tvs_prime = tvs_prime + PI * thvs_prime
   ps_prime = ps_prime - thvs_prime *rcp * tvs/ps * Pi
   thvs_prime = 0.0

   ! 3.1 Potential temperature on the lowest half sigma level
   !     ----------------------------------------------------
   ts_prime = ts_prime + Pi * ths_prime
   ps_prime = ps_prime - ths_prime *  rcp * ts/ps * Pi
   ths_prime = 0.0

   ! 2.0 Calculate the virtual temperature
   ! =====================================

   ! 2.2 Compute the ground saturated mixing ratio and the ground virtual 
   !     temperature
   !     ----------------------------------------------------------------

   tg_prime = tg_prime + tvg_prime * (1.0 + 0.608 * qg)
   qg_prime = qg_prime + tvg_prime * 0.608 * tg
   tvg_prime = 0.0

   qg       = qg_s
   qg_prime = qg_prime / (1.0 - qg)**2

   if ( present(qsfc_wrf) ) then
      if ( qsfc_wrf > 0.0 ) then
        qg_prime = 0.0
      end if
   endif

   call da_tp_to_qs_adj1(tg, psfc, eg, tg_prime, psfc_prime, qg_prime)

   ! 2.1 Compute Virtual temperature on the lowest half sigma level
   !     ---------------------------------------------------------

   ts_prime = ts_prime + tvs_prime  * (1.0 + 0.608 * qs)
   qs_prime = qs_prime + tvs_prime *  0.608 * ts
   tvs_prime = 0.0

   ! 2.0 Convert specific humidity to mixing ratio
   !     ---------------------------------------------------------

   qs_spechum_prime = qs_spechum_prime + qs_prime / (1.0 - qs_spechum)**2 
   qs_prime = 0.0


   if (trace_use) call da_trace_exit("da_sfc_wtq_adj")

end subroutine da_sfc_wtq_adj


