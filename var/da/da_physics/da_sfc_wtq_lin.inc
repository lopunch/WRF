subroutine da_sfc_wtq_lin(psfc, tg, ps, ts, qs_spechum, us, vs, regime,           &
   psfc_prime,tg_prime,ps_prime,ts_prime,qs_spechum_prime, &
   us_prime,vs_prime, hs,roughness,xland,dx,         & ! Modified by Eric Chiang (JULY 2010)
   u10_prime,v10_prime,t2_prime,q2_prime,          &
   has_lsm, qsfc_wrf, ust_wrf, mol_wrf, hfx, qfx, pblh) 

   !---------------------------------------------------------------------------
   ! Purpose: Calculate the  10m wind, 2m temperature and moisture based on the
   ! similarity theory/
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

   real, intent (in)  :: regime
   real, intent (in)  :: ps , ts , qs_spechum , us, vs, psfc, tg
   real, intent (in)  :: ps_prime, ts_prime, qs_spechum_prime, us_prime, vs_prime, psfc_prime, tg_prime
   real, intent (in)  :: hs, roughness, xland
   real, intent (out) :: u10_prime, v10_prime, t2_prime, q2_prime
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
   real :: rr, rz, r2
   integer :: nn, nz, n2

   real :: Vc2_prime, Va2_prime, V2_prime
   real :: rib_prime, xx_prime, yy_prime
   real :: psiw_prime, psiz_prime, mol_prime, ust_prime, &
           hol_prime, holz_prime, hol2_prime
   real :: psim_prime, psimz_prime, psim2_prime, &
           psih_prime, psihz_prime, psih2_prime
   real :: psit_prime, psit2_prime, &
           psiq_prime, psiq2_prime
   real :: qg_prime, tvg_prime, tvs_prime
   real :: ths_prime, thg_prime, thvs_prime, thvg_prime 
   real :: visc, restar, z0t, z0q
   real :: wspd, wspd_prime, visc_prime, restar_prime, z0t_prime, z0q_prime, &
           rhox, rhox_prime, fluxc, fluxc_s, fluxc_prime, cqs2, cqs2_prime, chs2, chs2_prime, &
           rho, rho_prime, vc
   real :: rr_prime, rz_prime, r2_prime
   real :: qs, qs_prime

   real, parameter :: ka = 2.4E-5

   integer :: iregime

   if (trace_use) call da_trace_entry("da_sfc_wtq_lin")

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


   ! 2. Calculate the virtual temperature
   ! =====================================

   ! 2.0 Convert specific humidity to mixing ratio
   !     ---------------------------------------------------------

   qs = qs_spechum / (1.0 - qs_spechum)
   qs_prime = qs_spechum_prime / (1.0 - qs_spechum)**2

   ! 2.1 Compute Virtual temperature on the lowest half sigma level
   !     ---------------------------------------------------------

   tvs_prime  = ts_prime * (1.0 + 0.608 * qs) + 0.608 * ts * qs_prime
   tvs  = ts * (1.0 + 0.608 * qs)

   ! 2.2 Compute the ground saturated mixing ratio and the ground virtual 
   !     temperature
   !     ----------------------------------------------------------------

   call da_tp_to_qs(tg, psfc, eg, qg)
   call da_tp_to_qs_lin1(tg, psfc, eg, tg_prime, psfc_prime, qg_prime)

   if ( present(qsfc_wrf) ) then
      if ( qsfc_wrf > 0.0 ) then
         qg = qsfc_wrf
         qg_prime = 0.0
      end if
   endif

   qg_prime = qg_prime / (1.0 - qg)**2   !yclo
   qg       = qg/(1.0-qg)                !hcl convert to mixing ratio

   tvg_prime  = tg_prime * (1.0 + 0.608 * qg) + 0.608 * tg * qg_prime
   tvg  = tg * (1.0 + 0.608 * qg)

   ! 3.  Compute the potential temperature and virtual potential temperature
   ! =======================================================================

   ! 3.1 Potential temperature on the lowest half sigma level
   !     ----------------------------------------------------

   Pi = (100000.0 / ps) ** rcp
   ths_prime  = (ts_prime - ps_prime * rcp * ts/ps) * Pi 
   ths  = ts * Pi

   ! 3.2 Virtual potential temperature on the lowest half sigma level
   !     ------------------------------------------------------------

   thvs_prime  = (tvs_prime - ps_prime * rcp * tvs/ps) * Pi 
   thvs = tvs * Pi

   ! 3.3 Potential temperature at the ground
   !     -----------------------------------

   Pi = (100000.0 / psfc) ** rcp
   thg_prime  = (tg_prime - psfc_prime * rcp * tg/psfc) * Pi 
   thg  = tg * Pi

   ! 3.4 Virtual potential temperature at ground
   !     ---------------------------------------

   thvg_prime  = (tvg_prime - psfc_prime * rcp * tvg/psfc) * Pi
   thvg = tvg * Pi

   ! 4.  BULK RICHARDSON NUMBER AND MONI-OBUKOV LENGTH
   ! =================================================

   ! 4.1 Velocity
   !     --------
   
   ! Wind speed:

   Va2_prime =  2.0*us*us_prime + 2.0*vs*vs_prime
   Va2 =   us*us + vs*vs
    
   ! Convective velocity:

   vconv_wrf = .false.
   if ( present(hfx) .and. present(qfx) .and. present(pblh) ) then
      ! calculating vconv over land following wrf method
      if ( pblh > 0.0 ) then
         vconv_wrf = .true.
      end if
   end if


   if (thvg >= thvs) then
      Vc2_prime = thvg_prime - thvs_prime
      Vc2 = thvg - thvs
   else
      Vc2_prime = 0.0
      Vc2 = 0.0
   end if
   if ( xland < 1.5 ) then !land
      if ( vconv_wrf ) then
         ! following the calculation as in module_sf_sfclay.F
         rhox = psfc/(gas_constant*tvg)
         rhox_prime = (psfc_prime - psfc/tvg*tvg_prime)/(gas_constant*tvg)

         fluxc_s = hfx*(thg-ths)/rhox/cp+0.608*tvg*qfx*(qg-qs)/rhox
         if(fluxc_s >= 0.) then
            fluxc = fluxc_s 
            fluxc_prime = (-fluxc*rhox_prime + 0.608*qfx*(qg-qs)*tvg_prime)/rhox + &
                          hfx*(thg_prime-ths_prime)/rhox/cp+0.608*tvg*qfx*(qg_prime-qs_prime)/rhox
         else
            fluxc = 0.
            fluxc_prime = 0.
         end if

         vc = (gravity/tg*pblh*fluxc)**0.33
         vc2 = vc*vc
         if(vc == 0.0) then
            Vc2_prime = 0.0
         else
            Vc2_prime = 2.*vc*( &
                0.33/(gravity/tg*pblh*fluxc)**0.67*gravity*pblh*&
                (-fluxc*tg_prime/tg + fluxc_prime) ) /tg
         end if
      end if
   end if

 
   ! Calculate Mahrt and Sun low-res correction                    ! Add by Eric Chiang ( July 2010 )
   ! dx is a constant, so vsgd is also a constant, 
   ! the perturnations of vsgd and vsgd2 are zero. (YRG, 09/15/2011)
                                     
   vsgd = 0.32 * (max(dx/5000.-1.,0.))**0.33                       ! Add by Eric Chiang ( July 2010 )
   vsgd2 = vsgd * vsgd                                             ! Add by Eric Chiang ( July 2010 )

   ! V2_prime should be computed before used below. (YRG, 09/15/2011)    
   V2_prime = Va2_prime + Vc2_prime
   V2  = Va2 + Vc2 + vsgd2                                         ! Modified by Eric Chiang ( July 2010 )

   !yclo add
   wspd = sqrt(v2)
   wspd_prime = V2_prime / (2. * wspd)

   if(wspd < 0.1) wspd_prime = 0.0
   wspd = max(wspd,0.1)

   v2 = wspd*wspd
   V2_prime = 2. * wspd * wspd_prime
   !yclo end add

   ! 4.2 Bulk richardson number
   !     ----------------------

   Pi = gravity * hs / (ths*V2)
   rib_prime = (thvs_prime - thvg_prime   &
              - (thvs-thvg)/V2  * V2_prime &
              - (thvs-thvg)/ths * ths_prime) * Pi 
   rib = (thvs - thvg) * Pi

   ! if previously unstable, do not let into regime 1 and 2
   if ( present(mol_wrf) ) then
      if ( mol_wrf < 0.0 ) then
        !rib = min(rib, 0.0)
         if(rib > 0.) then
            rib = 0.0
            rib_prime = 0.0
         end if
      end if
   end if

   ! 5.4.1  Calculate   ust, m/L (mol), h/L (hol)
   !        --------------------------

   psim = 0.0
   psih = 0.0
   psim_prime = 0.0
   psih_prime = 0.0

   ! Friction speed

   if ( present(ust_wrf) ) then
      if ( ust_wrf > 0.0 ) then
         use_ust_wrf = .true.
         ust = ust_wrf
      end if
   end if

   if(.not. use_ust_wrf ) then
      ust = k_kar * sqrt(v2) /(gzsoz0 - psim)
      ust_prime = (0.5/V2 * v2_prime + psim_prime /(gzsoz0 - psim)) * ust
   else
      ust_prime = 0.0
   end if

   ! Heat flux factor

   if ( present(mol_wrf) ) then
      mol = mol_wrf
      mol_prime = 0.0
   else 
      mol = k_kar * (ths - thg)/(gzsoz0 - psih)
      mol_prime = k_kar * (ths_prime - thg_prime)/(gzsoz0 - psih)
   end if

   ! 5.  CALCULATE PSI BASED UPON REGIME
   ! =======================================

   iregime = int(regime)

   select case (iregime) 

   ! 5.1 Stable conditions (REGIME 1)
   !     ---------------------------

   case (1);

      psim = -10.0*gzsoz0
      psim_prime = 0.0
      psim = max(psim,-10.0)

      psimz = h10/hs * psim
      psimz_prime = h10/hs * psim_prime
      psimz = max(psimz,-10.0)

      psim2 = h2 /hs * psim
      psim2_prime = h2/hs * psim_prime
      psim2 = max(psim2,-10.0)

      psih = psim
      psihz = psimz
      psih2 = psim2

      psih_prime = psim_prime
      psihz_prime = psimz_prime
      psih2_prime = psim2_prime

   ! 5.2 Mechanically driven turbulence (REGIME 2)
   !     ------------------------------------------

   case (2);

       psim = (-5.0 * rib) * gzsoz0 / (1.1 - 5.0*rib)
!      psim_prime = -5.5 * rib * gzsoz0 / (1.1 - 5.0*rib)**2
       psim_prime = -0.22 * gzsoz0 / (0.22-rib)**2 * rib_prime
       if(psim < -10.0) then
           psim = -10.0
           psim_prime = 0.0
       end if

       psimz = h10/hs * psim
       psimz_prime = h10/hs * psim_prime
       if(psimz < -10.0) then
           psimz = -10.0
           psimz_prime = 0.0
       end if

       psim2 = h2/hs * psim
       psim2_prime = h2/hs * psim_prime
       if(psim2 < -10.0) then
           psim2 = -10.0
           psim2_prime = 0.0
       end if

       psih = psim
       psihz = psimz
       psih2 = psim2
       psih_prime = psim_prime
       psihz_prime = psimz_prime
       psih2_prime = psim2_prime

   ! 5.3 Unstable Forced convection (REGIME 3)
   !     -------------------------------------

   case (3);

      psim_prime  = 0.0
      psimz_prime = 0.0
      psim2_prime = 0.0

      psim  = 0.0
      psimz = 0.0
      psim2 = 0.0

      psih_prime  = psim_prime
      psihz_prime = psimz_prime
      psih2_prime = psim2_prime

      psih  = psim
      psihz = psimz
      psih2 = psim2


      ! 5.4 Free convection (REGIME 4)
      !     --------------------------

   case (4);

      ! Calculate psi m and pshi h using iteration method

      cc = 2.0 * atan(1.0)

      ! do k = 1 , k_iteration


      ! Ratio of PBL height to Monin-Obukhov length

      if (ust .LT. 0.01) then
         hol_prime = rib_prime * gzsoz0
         hol = rib * gzsoz0
      else
         hol = k_kar * gravity * hs * mol / (ths * ust * ust)
         hol_prime = k_kar * gravity * hs * mol_prime / (ths * ust * ust) + &
                    (- ths_prime / ths - 2.0* ust_prime / ust) * hol
      end if

      ! 5.4.2  Calculate n, nz, R, Rz
      !        --------------------------

      holz = (h10 / hs) * hol
      holz_prime = (h10 / hs) * hol_prime
      hol2 = (h2 / hs) * hol
      hol2_prime = (h2 / hs) * hol_prime

      if (hol > 0.0) then
         hol_prime = 0.0
         hol = 0.0
      end if
      if (hol < -9.9999) then
         hol_prime = 0.0
         hol = -9.9999
      end if

      if (holz > 0.0) then
         holz_prime = 0.0
         holz = 0.0
      end if
      if (holz < -9.9999) then
         holz_prime = 0.0
         holz = -9.9999
      end if

      if (hol2 > 0.0) then
         hol2_prime = 0.0
         hol2 = 0.0
      end if
      if (hol2 < -9.9999) then
         hol2_prime = 0.0
         hol2 = -9.9999
      end if

      ! 5.4.3 Calculate Psim & psih
      !        --------------------------

      if ( use_table ) then
         ! Using the look-up table:
         nn = int(-100.0 * hol)
         rr = (-100.0 * hol) - nn
         rr_prime = -100.0 * hol_prime
         psim = psimtb(nn) + rr * (psimtb(nn+1) - psimtb(nn))
         psim_prime = rr_prime * (psimtb(nn+1) - psimtb(nn))
         psih = psihtb(nn) + rr * (psihtb(nn+1) - psihtb(nn))
         psih_prime = rr_prime * (psihtb(nn+1) - psihtb(nn))
      else
         ! Using the continuous function:
         xx_prime = -4.0* hol_prime /((1.0 - 16.0 * hol) ** 0.75)
         xx = (1.0 - 16.0 * hol) ** 0.25
         yy_prime = 2.0* xx * xx_prime /(1.0+xx*xx)
         yy = log((1.0+xx*xx)/2.0)
         psim_prime = 2 * xx_prime *(1.0/(1.0+xx)-1.0/(1+xx*xx)) + yy_prime 
         psim = 2.0 * log((1.0+xx)/2.0) + yy - 2.0 * atan(xx) + cc
         psih_prime = 2.0 * yy_prime
         psih = 2.0 * yy
      end if

      if ( use_table ) then
         ! Using the look-up table:
         nz = int(-100.0 * holz)
         rz = (-100.0 * holz) - nz
         rz_prime = -100.0 * holz_prime
         psimz = psimtb(nz) + rz * (psimtb(nz+1) - psimtb(nz))
         psimz_prime = rz_prime * (psimtb(nz+1) - psimtb(nz))
         psihz = psihtb(nz) + rz * (psihtb(nz+1) - psihtb(nz))
         psihz_prime = rz_prime * (psihtb(nz+1) - psihtb(nz))
      else
         ! Using the continuous function:
         xx_prime = -4.0* holz_prime /((1.0 - 16.0 * holz) ** 0.75)
         xx = (1.0 - 16.0 * holz) ** 0.25
         yy_prime = 2.0* xx * xx_prime /(1.0+xx*xx)
         yy = log((1.0+xx*xx)/2.0)
         psimz_prime = 2.0* xx_prime *(1.0/(1.0+xx)-1.0/(1+xx*xx)) + yy_prime
         psimz = 2.0 * log((1.0+xx)/2.0) + yy - 2.0 * atan(xx) + cc
         psihz_prime = 2.0 * yy_prime
         psihz = 2.0 * yy
      end if

      if ( use_table ) then
         ! Using the look-up table:
         n2 = int(-100.0 * hol2)
         r2 = (-100.0 * hol2) - n2
         r2_prime = -100.0 * hol2_prime
         psim2 = psimtb(n2) + r2 * (psimtb(n2+1) - psimtb(n2))
         psim2_prime = r2_prime * (psimtb(n2+1) - psimtb(n2))
         psih2 = psihtb(n2) + r2 * (psihtb(n2+1) - psihtb(n2))
         psih2_prime = r2_prime * (psihtb(n2+1) - psihtb(n2))
      else
         ! Using the continuous function:
         xx_prime = -4.0* hol2_prime /((1.0 - 16.0 * hol2) ** 0.75)
         xx = (1.0 - 16.0 * hol2) ** 0.25
         yy_prime = 2.0* xx * xx_prime /(1.0+xx*xx)
         yy = log((1.0+xx*xx)/2.0)
         psim2_prime = 2.0* xx_prime *(1.0/(1.0+xx)-1.0/(1+xx*xx)) + yy_prime
         psim2 = 2.0 * log((1.0+xx)/2.0) + yy - 2.0 * atan(xx) + cc
         psih2_prime = 2.0 * yy_prime
         psih2 = 2.0 * yy
      end if

      ! end do 

      ! 5.4.4 Define the limit value for psim & psih
      !        --------------------------

      if (psim > 0.9*gzsoz0) then
         psim = 0.9*gzsoz0
         psim_prime = 0.0
      end if
      if (psimz > 0.9*gz10oz0) then
         psimz_prime = 0.0
         psimz = 0.9*gz10oz0
      end if
      if (psim2 > 0.9*gz2oz0) then
         psim2_prime = 0.0
         psim2 = 0.9*gz2oz0
      end if
      if (psih > 0.9*gzsoz0) then
         psih_prime = 0.0
         psih = 0.9*gzsoz0
      end if
      if (psihz > 0.9*gz10oz0) then
         psihz_prime = 0.0
         psihz = 0.9*gz10oz0
      end if
      if (psih2 > 0.9*gz2oz0) then
         psih2_prime = 0.0
         psih2 = 0.9*gz2oz0
      end if

   case default;
      write(unit=message(1),fmt='(A,I2,A)') &
         "Regime=",iregime," is invalid."
      call da_error(__FILE__,__LINE__,message(1:1))

   end select

   ! 6.  CALCULATE PSI FOR WinD, TEMPERATURE AND MOISTURE
   ! =======================================

   psiw_prime = - psim_prime
   psiw = gzsoz0 - psim
   psiz_prime = - psimz_prime
   psiz = gz10oz0 - psimz

   if(gzsoz0 - psih < 2.0) then
      psit_prime = 0.0
      psit = 2.0
   else
      psit_prime = - psih_prime
      psit = gzsoz0 - psih
   end if

   psit2_prime = - psih2_prime
   psit2 = gz2oz0 - psih2

   if ( .not. use_ust_wrf ) then
      ust = k_kar * sqrt(v2) /(gzsoz0 - psim)
      ust_prime = (0.5/V2 * v2_prime + psim_prime /(gzsoz0 - psim)) * ust
   end if

   psiq_prime  = k_kar*hs/(ka*(k_kar*ust*hs/ka + hs / zq0))*ust_prime - psih_prime
   psiq2_prime = k_kar*h2/(ka*(k_kar*ust*h2/ka + h2 / zq0))*ust_prime - psih2_prime

   psiq  = log(k_kar*ust*hs/ka + hs / zq0) - psih
   psiq2 = log(k_kar*ust*h2/ka + h2 / zq0) - psih2

   !V3.7, as in module_sf_sfclay.F
   if ( xland >= 1.5 ) then !water
      visc = (1.32+0.009*(ts-273.15))*1.e-5
      visc_prime = 0.009*1.e-5*ts_prime

      restar = ust*z0/visc
      restar_prime = z0 / visc * (ust_prime - ust * visc_prime / visc)

      z0t = (5.5e-5)*(restar**(-0.60))
      z0t_prime = -0.6 * z0t / restar * restar_prime

      if(z0t > 1.0e-4) z0t_prime = 0.0
      if(z0t < 2.0e-9) z0t_prime = 0.0
      z0t = min(z0t,1.0e-4)
      z0t = max(z0t,2.0e-9)

      z0q = z0t
      z0q_prime = z0t_prime

      psiq  = log((hs+z0q)/z0q)-psih
      psiq_prime = -hs / (z0q * (hs+z0q)) * z0q_prime - psih_prime

      psit  = log((hs+z0t)/z0t)-psih
      psit_prime = -hs / (z0t * (hs+z0t)) * z0t_prime - psih_prime

      psiq2 = log((2.+z0q)/z0q)-psih2
      psiq2_prime = -2. / (z0q * (2.+z0q)) * z0q_prime - psih2_prime

      psit2 = log((2.+z0t)/z0t)-psih2
      psit2_prime = -2. / (z0t * (2.+z0t)) * z0t_prime - psih2_prime

      if(psiq < 2.) psiq_prime = 0.
      psiq  = max(psiq,  2.)

      if(psit < 2.) psit_prime = 0.
      psit  = max(psit,  2.)

      if(psiq2 < 2.) psiq2_prime = 0.
      psiq2 = max(psiq2, 2.)

      if(psit2 < 2.) psit2_prime = 0.
      psit2 = max(psit2, 2.)
   end if


   ! 7.  CALCULATE THE PERTURBATIONS for 10M WinD, 2M TEMPERATURE AND MOISTURE
   ! =======================================

   Pi = psiz / psiw
   u10_prime= (us_prime + us/psiz * psiz_prime - us/psiw * psiw_prime) * Pi
   v10_prime= (vs_prime + vs/psiz * psiz_prime - vs/psiw * psiw_prime) * Pi

   t2_prime = ((1.0-psit2/psit) * thg_prime + (ths_prime + &
                           (ths - thg)/psit2 * psit2_prime - &
                           (ths - thg)/psit  * psit_prime) * psit2/psit &
             + rcp*(thg + (ths - thg)*psit2/psit)/psfc * psfc_prime) &
             * (psfc/100000.0)**rcp

   q2_prime = (1.0-psiq2/psiq) * qg_prime + psiq2/psiq * qs_prime + &
              (qs -qg)*(psiq2/psiq) * (psiq2_prime/psiq2 - psiq_prime/psiq)

   if ( present(has_lsm) ) then
      if ( has_lsm ) then
         !cqs2: 2m surface exchange coefficient for moisture
         !chs2: 2m surface exchange coefficient for heat
         cqs2 = ust*k_kar/psiq2
         cqs2_prime = (ust_prime - ust/psiq2*psiq2_prime)*k_kar/psiq2
         if (xland .ge. 1.5) then
            !water
            chs2 = ust*k_kar/psit2
            chs2_prime = (ust_prime - ust/psit2*psit2_prime)*k_kar/psit2
         else
            !land
            chs2 = cqs2 !as in subroutine lsm in phys/module_sf_noahdrv.F
            chs2_prime = cqs2_prime
         end if

         !re-calculate T2/Q2 as in module_sf_sfcdiags.F
         rho  = psfc/(gas_constant*tg)
         rho_prime = (psfc_prime - psfc/tg*tg_prime)/(gas_constant*tg)
         if ( cqs2 < 1.e-5 ) then
            q2_prime = qg_prime
         else
            if ( present(qfx) ) then
               q2_prime = qg_prime + qfx*(qg-qs)/(rho*cqs2)**2 * (rho_prime*cqs2+rho*cqs2_prime) &
                                   - qfx*(qg_prime-qs_prime)/(rho*cqs2)
                                      
            end if
         end if
         if ( chs2 < 1.e-5 ) then
            t2_prime = tg_prime
         else
            if ( present(hfx) ) then
               t2_prime = tg_prime + hfx*(thg-ths)/(rho*cp*chs2)**2 * cp*(rho_prime*chs2+rho*chs2_prime) &
                                   - hfx*(thg_prime-ths_prime)/(rho*cp*chs2)
            end if
         end if
      end if
   end if


   if (trace_use) call da_trace_exit("da_sfc_wtq_lin")

end subroutine da_sfc_wtq_lin


