*deck,usermat      USERDISTRIB  parallel                                     gal
      subroutine usermat(
     &     matId, elemId, kDomIntPt, kLayer, kSectPt,
     &     ldstep, isubst, keycut,
     &     nDirect, nShear, ncomp, nStatev, nProp,
     &     Time, dTime, Temp, dTemp,
     &     stress, ustatev, dsdePl, sedEl, sedPl, epseq,
     &     Strain, dStrain, epsPl, prop, coords,
     &     var0, defGrad_t, defGrad,
     &     tsstif, epsZZ,
     &     cutFactor, pVolDer, hrmflg,
     &     var1, var2, var3, var4, var5, var6, var7, var8)
#include "impcom.inc"
      INTEGER matId, elemId, kDomIntPt, kLayer, kSectPt,
     &        ldstep, isubst, keycut, nDirect, nShear, ncomp,
     &        nStatev, nProp
      DOUBLE PRECISION Time, dTime, Temp, dTemp, sedEl, sedPl,
     &                 epseq, epsZZ, cutFactor, hrmflg
      DOUBLE PRECISION stress(ncomp), ustatev(nStatev),
     &                 dsdePl(ncomp,ncomp), pVolDer(3),
     &                 Strain(ncomp), dStrain(ncomp),
     &                 epsPl(ncomp), prop(nProp),
     &                 coords(3), defGrad(3,3), defGrad_t(3,3),
     &                 tsstif(2)
      DOUBLE PRECISION var0,var1,var2,var3,var4,var5,var6,
     &                 var7,var8
      IF (ncomp .GE. 4) THEN
         call usermat3d_epd(
     &        matId, elemId, kDomIntPt, kLayer, kSectPt,
     &        ldstep, isubst, keycut,
     &        nDirect, nShear, ncomp, nStatev, nProp,
     &        Time, dTime, Temp, dTemp,
     &        stress, ustatev, dsdePl, sedEl, sedPl, epseq,
     &        Strain, dStrain, epsPl, prop, coords,
     &        var0, defGrad_t, defGrad,
     &        tsstif, epsZZ, cutFactor,
     &        var1, var2, var3, var4, var5, var6, var7, var8)
      END IF
      return
      end

*deck,usermat3d_epd USERDISTRIB  parallel                                    gal
      subroutine usermat3d_epd(
     &     matId, elemId, kDomIntPt, kLayer, kSectPt,
     &     ldstep, isubst, keycut,
     &     nDirect, nShear, ncomp, nStatev, nProp,
     &     Time, dTime, Temp, dTemp,
     &     stress, ustatev, dsdePl, sedEl, sedPl, epseq,
     &     Strain, dStrain, epsPl, prop, coords,
     &     var0, defGrad_t, defGrad,
     &     tsstif, epsZZ, cutFactor,
     &     var1, var2, var3, var4, var5, var6, var7, var8)
#include "impcom.inc"
      INTEGER matId, elemId, kDomIntPt, kLayer, kSectPt,
     &        ldstep, isubst, keycut, nDirect, nShear, ncomp,
     &        nStatev, nProp
      DOUBLE PRECISION Time, dTime, Temp, dTemp, sedEl, sedPl,
     &                 epseq, epsZZ, cutFactor
      DOUBLE PRECISION stress(ncomp), ustatev(nStatev),
     &                 dsdePl(ncomp,ncomp), Strain(ncomp),
     &                 dStrain(ncomp), epsPl(ncomp), prop(nProp),
     &                 coords(3), defGrad(3,3), defGrad_t(3,3),
     &                 tsstif(2)
      DOUBLE PRECISION var0,var1,var2,var3,var4,var5,var6,
     &                 var7,var8
      DOUBLE PRECISION Ec, nu, st0, et0, stp, etp, Gf, wc
      DOUBLE PRECISION scp, ecp, scu, ecu
      DOUBLE PRECISION Et, sig_floor
      DOUBLE PRECISION Le, etp_obj, etu_obj, ecu_obj, ncf, vtmp
      DOUBLE PRECISION eps_t(6)
      DOUBLE PRECISION stress_old(6), eps_in_old(6)
      DOUBLE PRECISION eps_in(6), eps_e_fin(6), deps_in(6)
      DOUBLE PRECISION kap_t_old, kap_c_old
      DOUBLE PRECISION d_t_new, d_c_new
      DOUBLE PRECISION kap_t_new, kap_c_new
      DOUBLE PRECISION kap_t_pert, kap_c_pert
      DOUBLE PRECISION d_t_pert, d_c_pert
      DOUBLE PRECISION sig_pert(6), eps_pert(6)
      DOUBLE PRECISION eps_in_dummy(6)
      DOUBLE PRECISION De(6,6), sig_trial(6), sig_m(3,3)
      DOUBLE PRECISION eig_val(3), Q(3,3)
      DOUBLE PRECISION eps_m(3,3), eps_eig(3), Qe(3,3)
      DOUBLE PRECISION eps_plus(3), eps_minus(3)
      DOUBLE PRECISION eig_corr(3), eig_plus(3), eig_minus(3)
      DOUBLE PRECISION eig_dmg(3)
      DOUBLE PRECISION ebar_t, ebar_c
      DOUBLE PRECISION d_eff, amax, lam_dom, xval, vden, H_t
      DOUBLE PRECISION delta, epsmax, tau, fac
      INTEGER i, j, k, dom_k

      Ec  = dabs(prop(1))
      nu  = dabs(prop(2))
      st0 = dabs(prop(3))
      et0 = dabs(prop(4))
      stp = dabs(prop(5))
      etp = dabs(prop(6))
      Gf =  dabs(prop(7))   ! NEW: Now reading Fracture Energy (N/mm) instead of etu
      scp = dabs(prop(8))
      ecp = dabs(prop(9))
      scu = dabs(prop(10))
      ecu = dabs(prop(11))
      Et  = st0 / et0
      sig_floor = 1.0d-3 * stp

      do i = 1, 6
         stress_old(i) = stress(i)
         eps_in_old(i) = ustatev(i)
      end do
      kap_t_old = ustatev(7)
      kap_c_old = ustatev(8)

C     ===================================================================
C     INITIALIZE CHARACTERISTIC LENGTH (MESH OBJECTIVITY)
C     ===================================================================
C     Read Le injected by the APDL INISTATE command for this element
      Le = ustatev(12)
      
C     Fail-safe: If INISTATE wasn't run, fallback to prop(12)
      if (Le .lt. 1.0d-10) then
         Le = dabs(prop(12))
         if (Le .lt. 1.0d-10) Le = 1.0d0
         ustatev(12) = Le
      end if

C     --- UPDATED ECC SCALING LOGIC USING FRACTURE ENERGY ---
C     Hardening is a continuum multiple-cracking phenomenon. Do NOT scale by Le.
      etp_obj = etp 

C     Softening represents localized macro-cracking. Scale using Gf and Le.
C     Calculate physical crack opening (wc) assuming linear softening
      wc = (2.0d0 * Gf) / stp
      etu_obj = etp_obj + (wc / Le) 
C     -------------------------------------------------------

      if (etu_obj .le. etp_obj + 1.0d-10) etu_obj = etp_obj + 1.0d-10
    
C     --- COMPRESSION SCALING (Crushing Band Theory) ---
C     --- Scaling the experimental crushing strain by 203 mm specimen height 
      ecu_obj = ecp + ((ecu - ecp)* 203.0d0) / Le
      if (ecu_obj .le. ecp + 1.0d-10) ecu_obj = ecp + 1.0d-10
	  

      vtmp = ecp * ecu_obj
      if (vtmp .gt. 1.0d-15 .and. scp .gt. scu) then
         ncf = (scu/(scp-scu)) * (ecp-ecu_obj)**2 / vtmp
      else
         ncf = 1.5d0
      end if
      if (ncf .lt. 1.0d0) ncf = 1.0d0

      do i = 1, 6
         eps_t(i) = Strain(i) + dStrain(i)
      end do

C ===================================================================
C  STRESS UPDATE (baseline) - with viscous regularization
C ===================================================================
      call get_el_mat_epd(Ec, nu, De)
      do i = 1, 6
         sig_trial(i) = 0.0d0
         do j = 1, 6
            sig_trial(i) = sig_trial(i) + De(i,j) * eps_t(j)
         end do
      end do

      sig_m(1,1) = sig_trial(1)
      sig_m(2,2) = sig_trial(2)
      sig_m(3,3) = sig_trial(3)
      sig_m(1,2) = sig_trial(4)
      sig_m(2,1) = sig_trial(4)
      sig_m(2,3) = sig_trial(5)
      sig_m(3,2) = sig_trial(5)
      sig_m(1,3) = sig_trial(6)
      sig_m(3,1) = sig_trial(6)
      call jacobi_3x3_epd(sig_m, eig_val, Q)

      do k = 1, 3
         if (eig_val(k) .gt. 0.0d0) then
            eig_corr(k) = eig_val(k) * (Et/Ec)
         else
            eig_corr(k) = eig_val(k)
         end if
         eig_plus(k)  = max(eig_corr(k), 0.0d0)
         eig_minus(k) = min(eig_corr(k), 0.0d0)
      end do

      eps_m(1,1) = eps_t(1)
      eps_m(2,2) = eps_t(2)
      eps_m(3,3) = eps_t(3)
      eps_m(1,2) = eps_t(4) * 0.5d0
      eps_m(2,1) = eps_t(4) * 0.5d0
      eps_m(2,3) = eps_t(5) * 0.5d0
      eps_m(3,2) = eps_t(5) * 0.5d0
      eps_m(1,3) = eps_t(6) * 0.5d0
      eps_m(3,1) = eps_t(6) * 0.5d0
      call jacobi_3x3_epd(eps_m, eps_eig, Qe)

      do k = 1, 3
         eps_plus(k)  = max(eps_eig(k), 0.0d0)
         eps_minus(k) = min(eps_eig(k), 0.0d0)
      end do

      call compute_eq_strain_epd(eps_plus, 1.0d0, nu, ebar_t)
      do k = 1, 3
         eps_minus(k) = -eps_minus(k)
      end do
      call compute_eq_strain_epd(eps_minus, 1.0d0, nu, ebar_c)
      do k = 1, 3
         eps_minus(k) = -eps_minus(k)
      end do

C     Viscous regularization for baseline stress update
      tau = 1.0d-5
      if (dTime .gt. 0.0d0) then
         fac = 1.0d0 - exp(-dTime / tau)
      else
         fac = 1.0d0
      end if
      kap_t_new = kap_t_old + (max(kap_t_old, ebar_t) - kap_t_old) * fac
      kap_c_new = kap_c_old + (max(kap_c_old, ebar_c) - kap_c_old) * fac

      call eval_damage_t_epd(kap_t_new, Et, st0, et0, stp,
     &     etp_obj, etu_obj, sig_floor, d_t_new)
      call eval_damage_c_epd(kap_c_new, Ec, scp, ecp, scu,
     &     ecu_obj, ncf, d_c_new)
      d_t_new = max(min(d_t_new, 1.0d0-1.0d-6), 0.0d0)
      d_c_new = max(min(d_c_new, 1.0d0-1.0d-6), 0.0d0)

      do k = 1, 3
         eig_dmg(k) = (1.0d0-d_t_new)*eig_plus(k) +
     &                (1.0d0-d_c_new)*eig_minus(k)
      end do

      do i = 1, 6
         stress(i) = 0.0d0
      end do
      do k = 1, 3
         if (dabs(eig_dmg(k)) .lt. 1.0d-30) goto 30
         stress(1) = stress(1) + eig_dmg(k)*Q(1,k)*Q(1,k)
         stress(2) = stress(2) + eig_dmg(k)*Q(2,k)*Q(2,k)
         stress(3) = stress(3) + eig_dmg(k)*Q(3,k)*Q(3,k)
         stress(4) = stress(4) + eig_dmg(k)*Q(1,k)*Q(2,k)
         stress(5) = stress(5) + eig_dmg(k)*Q(2,k)*Q(3,k)
         stress(6) = stress(6) + eig_dmg(k)*Q(1,k)*Q(3,k)
 30      continue
      end do

      do i = 1, 6
         eps_in(i) = 0.0d0
      end do
      do k = 1, 3
         if (eig_val(k) .ge. 0.0d0) then
            vtmp = d_t_new * eig_val(k) / Ec
         else
            vtmp = d_c_new * eig_val(k) / Ec
         end if
         if (dabs(vtmp) .lt. 1.0d-30) goto 40
         eps_in(1) = eps_in(1) + vtmp*Q(1,k)*Q(1,k)
         eps_in(2) = eps_in(2) + vtmp*Q(2,k)*Q(2,k)
         eps_in(3) = eps_in(3) + vtmp*Q(3,k)*Q(3,k)
         eps_in(4) = eps_in(4) + 2.0d0*vtmp*Q(1,k)*Q(2,k)
         eps_in(5) = eps_in(5) + 2.0d0*vtmp*Q(2,k)*Q(3,k)
         eps_in(6) = eps_in(6) + 2.0d0*vtmp*Q(1,k)*Q(3,k)
 40      continue
      end do

C ===================================================================
C  NUMERICAL CONSISTENT TANGENT (6 perturbations)
C ===================================================================
      epsmax = 0.0d0
      do i = 1, 6
         if (dabs(eps_t(i)) .gt. epsmax) epsmax = dabs(eps_t(i))
      end do
      delta = 1.0d-8 * (1.0d0 + epsmax)

      do j = 1, 6
         do i = 1, 6
            eps_pert(i) = eps_t(i)
         end do
         eps_pert(j) = eps_pert(j) + delta

         kap_t_pert = kap_t_old
         kap_c_pert = kap_c_old

         call stress_only_epd(eps_pert, prop, kap_t_pert, kap_c_pert,
     &        Le, dTime, sig_pert, d_t_pert, d_c_pert,
     &        kap_t_pert, kap_c_pert)

         do i = 1, 6
            dsdePl(i,j) = (sig_pert(i) - stress(i)) / delta
         end do
      end do
      
C     Symmetrize the numerical Jacobian for the standard ANSYS solver
      do i = 1, 6
         do j = i, 6
            vtmp = 0.5d0 * (dsdePl(i,j) + dsdePl(j,i))
            dsdePl(i,j) = vtmp
            dsdePl(j,i) = vtmp
         end do
      end do	  

C ===================================================================
C  OVERRIDE TANGENT TO PREVENT PING-PONG DURING SOFTENING
C ===================================================================
C  Using kap_t_old and kap_c_old which are already in scope.
C  Manually reconstructing the elastic stiffness using prop(1) and prop(2)
C  to bypass the need for undeclared external matrices.

      IF (kap_t_old .gt. 0.0d0 .or. kap_c_old .gt. 0.0d0) THEN
         do i = 1, 6
            do j = 1, 6
               dsdePl(i,j) = 0.0d0
            end do
         end do
         
C        vtmp is already declared. Reusing it to store the E/nu factor.
C        Assuming prop(1) = E, prop(2) = nu
         vtmp = prop(1) / ((1.0d0 + prop(2))*(1.0d0 - 2.0d0 * prop(2)))
         
C        Diagonal axial terms
         dsdePl(1,1) = vtmp * (1.0d0 - prop(2))
         dsdePl(2,2) = dsdePl(1,1)
         dsdePl(3,3) = dsdePl(1,1)
         
C        Off-diagonal axial terms
         dsdePl(1,2) = vtmp * prop(2)
         dsdePl(2,1) = dsdePl(1,2)
         dsdePl(1,3) = dsdePl(1,2)
         dsdePl(3,1) = dsdePl(1,2)
         dsdePl(2,3) = dsdePl(1,2)
         dsdePl(3,2) = dsdePl(1,2)
         
C        Diagonal shear terms
         dsdePl(4,4) = vtmp * (0.5d0 - prop(2))
         dsdePl(5,5) = dsdePl(4,4)
         dsdePl(6,6) = dsdePl(4,4)
      END IF
	  

C ===================================================================
C  STATE UPDATE
C ===================================================================
      do i = 1, 6
         epsPl(i)   = eps_in(i)
         ustatev(i) = eps_in(i)
      end do
      ustatev(7)  = kap_t_new
      ustatev(8)  = kap_c_new
      ustatev(10) = d_t_new
      ustatev(11) = d_c_new
C     NOTE: Do not overwrite ustatev(12) with var1 here! 
C     ustatev(12) now permanently holds your characteristic length.

      sedEl = 0.0d0
      do i = 1, 6
         eps_e_fin(i) = eps_t(i) - eps_in(i)
         sedEl = sedEl + 0.5d0 * stress(i) * eps_e_fin(i)
      end do

      sedPl = 0.0d0
      do i = 1, 6
         deps_in(i) = eps_in(i) - eps_in_old(i)
      end do
      do i = 1, 6
         sedPl = sedPl +
     &        0.5d0*(stress_old(i)+stress(i))*deps_in(i)
      end do

      epseq = kap_t_new + kap_c_new
      tsstif(1) = dsdePl(4,4)
      tsstif(2) = dsdePl(5,5)

      return
      end

C ===================================================================
C  stress_only_epd : stress calculation WITHOUT state update
C ===================================================================
      subroutine stress_only_epd(eps_t, prop, kap_t_old, kap_c_old,
     &     Le, dTime, stress_out, d_t_new, d_c_new,
     &     kap_t_new, kap_c_new)
      DOUBLE PRECISION eps_t(6), prop(12)
      DOUBLE PRECISION kap_t_old, kap_c_old, Le, dTime
      DOUBLE PRECISION stress_out(6)
      DOUBLE PRECISION d_t_new, d_c_new
      DOUBLE PRECISION kap_t_new, kap_c_new

      DOUBLE PRECISION Ec, nu, st0, et0, stp, etp, Gf
      DOUBLE PRECISION scp, ecp, scu, ecu, Et, sig_floor
      DOUBLE PRECISION etp_obj, etu_obj, ecu_obj, ncf, vtmp
      DOUBLE PRECISION De(6,6), sig_trial(6), sig_m(3,3)
      DOUBLE PRECISION eig_val(3), Q(3,3)
      DOUBLE PRECISION eps_m(3,3), eps_eig(3), Qe(3,3)
      DOUBLE PRECISION eps_plus(3), eps_minus(3)
      DOUBLE PRECISION eig_corr(3), eig_plus(3), eig_minus(3)
      DOUBLE PRECISION eig_dmg(3)
      DOUBLE PRECISION ebar_t, ebar_c
      DOUBLE PRECISION tau, fac
      INTEGER i, j, k

      Ec  = dabs(prop(1))
      nu  = dabs(prop(2))
      st0 = dabs(prop(3))
      et0 = dabs(prop(4))
      stp = dabs(prop(5))
      etp = dabs(prop(6))
      Gf  = dabs(prop(7))    ! Read Gf instead of etu
      scp = dabs(prop(8))
      ecp = dabs(prop(9))
      scu = dabs(prop(10))
      ecu = dabs(prop(11))
      Et  = st0 / et0
      sig_floor = 1.0d-3 * stp

C     --- TENSION SCALING ---
      etp_obj = etp 
      etu_obj = etp_obj + ((2.0d0 * Gf) / stp) / Le 
      if (etu_obj .le. etp_obj + 1.0d-10) etu_obj = etp_obj + 1.0d-10

C     --- COMPRESSION SCALING ---
      ecu_obj = ecp + ((ecu - ecp) * 203.0d0) / Le
      if (ecu_obj .le. ecp + 1.0d-10) ecu_obj = ecp + 1.0d-10

      vtmp = ecp * ecu_obj
      if (vtmp .gt. 1.0d-15 .and. scp .gt. scu) then
         ncf = (scu/(scp-scu)) * (ecp-ecu_obj)**2 / vtmp
      else
         ncf = 1.5d0
      end if
      if (ncf .lt. 1.0d0) ncf = 1.0d0

      call get_el_mat_epd(Ec, nu, De)
      do i = 1, 6
         sig_trial(i) = 0.0d0
         do j = 1, 6
            sig_trial(i) = sig_trial(i) + De(i,j) * eps_t(j)
         end do
      end do

      sig_m(1,1) = sig_trial(1)
      sig_m(2,2) = sig_trial(2)
      sig_m(3,3) = sig_trial(3)
      sig_m(1,2) = sig_trial(4)
      sig_m(2,1) = sig_trial(4)
      sig_m(2,3) = sig_trial(5)
      sig_m(3,2) = sig_trial(5)
      sig_m(1,3) = sig_trial(6)
      sig_m(3,1) = sig_trial(6)
      call jacobi_3x3_epd(sig_m, eig_val, Q)

      do k = 1, 3
         if (eig_val(k) .gt. 0.0d0) then
            eig_corr(k) = eig_val(k) * (Et/Ec)
         else
            eig_corr(k) = eig_val(k)
         end if
         eig_plus(k)  = max(eig_corr(k), 0.0d0)
         eig_minus(k) = min(eig_corr(k), 0.0d0)
      end do

      eps_m(1,1) = eps_t(1)
      eps_m(2,2) = eps_t(2)
      eps_m(3,3) = eps_t(3)
      eps_m(1,2) = eps_t(4) * 0.5d0
      eps_m(2,1) = eps_t(4) * 0.5d0
      eps_m(2,3) = eps_t(5) * 0.5d0
      eps_m(3,2) = eps_t(5) * 0.5d0
      eps_m(1,3) = eps_t(6) * 0.5d0
      eps_m(3,1) = eps_t(6) * 0.5d0
      call jacobi_3x3_epd(eps_m, eps_eig, Qe)

      do k = 1, 3
         eps_plus(k)  = max(eps_eig(k), 0.0d0)
         eps_minus(k) = min(eps_eig(k), 0.0d0)
      end do

      call compute_eq_strain_epd(eps_plus, 1.0d0, nu, ebar_t)
      do k = 1, 3
         eps_minus(k) = -eps_minus(k)
      end do
      call compute_eq_strain_epd(eps_minus, 1.0d0, nu, ebar_c)
      do k = 1, 3
         eps_minus(k) = -eps_minus(k)
      end do

C     Small viscous regularization to smooth hardening oscillations
      tau = 1.0d-4
      if (dTime .gt. 0.0d0) then
         fac = 1.0d0 - exp(-dTime / tau)
      else
         fac = 1.0d0
      end if
      kap_t_new = kap_t_old + (max(kap_t_old, ebar_t) - kap_t_old) * fac
      kap_c_new = kap_c_old + (max(kap_c_old, ebar_c) - kap_c_old) * fac

      call eval_damage_t_epd(kap_t_new, Et, st0, et0, stp,
     &     etp_obj, etu_obj, sig_floor, d_t_new)
      call eval_damage_c_epd(kap_c_new, Ec, scp, ecp, scu,
     &     ecu_obj, ncf, d_c_new)
      d_t_new = max(min(d_t_new, 1.0d0-1.0d-6), 0.0d0)
      d_c_new = max(min(d_c_new, 1.0d0-1.0d-6), 0.0d0)

      do k = 1, 3
         eig_dmg(k) = (1.0d0-d_t_new)*eig_plus(k) +
     &                (1.0d0-d_c_new)*eig_minus(k)
      end do

      do i = 1, 6
         stress_out(i) = 0.0d0
      end do
      do k = 1, 3
         if (dabs(eig_dmg(k)) .lt. 1.0d-30) goto 30
         stress_out(1) = stress_out(1) + eig_dmg(k)*Q(1,k)*Q(1,k)
         stress_out(2) = stress_out(2) + eig_dmg(k)*Q(2,k)*Q(2,k)
         stress_out(3) = stress_out(3) + eig_dmg(k)*Q(3,k)*Q(3,k)
         stress_out(4) = stress_out(4) + eig_dmg(k)*Q(1,k)*Q(2,k)
         stress_out(5) = stress_out(5) + eig_dmg(k)*Q(2,k)*Q(3,k)
         stress_out(6) = stress_out(6) + eig_dmg(k)*Q(1,k)*Q(3,k)
 30      continue
      end do

      return
      end

C ===================================================================
C  Auxiliary routines
C ===================================================================
      subroutine get_el_mat_epd(E_val, nu_val, De_mat)
      DOUBLE PRECISION E_val, nu_val, De_mat(6,6)
      DOUBLE PRECISION fval, C11, C12, Gval
      INTEGER i, j
      fval = E_val / ((1.0d0+nu_val)*(1.0d0-2.0d0*nu_val))
      C11  = fval * (1.0d0 - nu_val)
      C12  = fval * nu_val
      Gval = E_val / (2.0d0*(1.0d0+nu_val))
      do i = 1, 6
         do j = 1, 6
            De_mat(i,j) = 0.0d0
         end do
      end do
      De_mat(1,1) = C11
      De_mat(2,2) = C11
      De_mat(3,3) = C11
      De_mat(1,2) = C12
      De_mat(1,3) = C12
      De_mat(2,3) = C12
      De_mat(2,1) = C12
      De_mat(3,1) = C12
      De_mat(3,2) = C12
      De_mat(4,4) = Gval
      De_mat(5,5) = Gval
      De_mat(6,6) = Gval
      return
      end

      subroutine jacobi_3x3_epd(A_in, eigval, Q)
      DOUBLE PRECISION A_in(3,3), eigval(3), Q(3,3)
      DOUBLE PRECISION A(3,3)
      DOUBLE PRECISION tau_r, t_r, c_r, s_r, g_r, h_r
      DOUBLE PRECISION ap_diag, aq_diag, off_sum
      INTEGER sweep, ip, iq, ir, i, j
      do i = 1, 3
         do j = 1, 3
            A(i,j) = A_in(i,j)
         end do
      end do
      do i = 1, 3
         do j = 1, 3
            Q(i,j) = 0.0d0
         end do
         Q(i,i) = 1.0d0
      end do
      do sweep = 1, 50
         off_sum = dabs(A(1,2)) + dabs(A(1,3)) + dabs(A(2,3))
         if (off_sum .lt. 1.0d-14) goto 100
         do ip = 1, 2
            do iq = ip+1, 3
               if (dabs(A(ip,iq)) .lt. 1.0d-15) goto 10
               tau_r = (A(iq,iq)-A(ip,ip)) / (2.0d0*A(ip,iq))
               if (tau_r .ge. 0.0d0) then
                  t_r = 1.0d0/(tau_r + dsqrt(1.0d0+tau_r*tau_r))
               else
                  t_r = -1.0d0/(-tau_r + dsqrt(1.0d0+tau_r*tau_r))
               end if
               c_r = 1.0d0 / dsqrt(1.0d0 + t_r*t_r)
               s_r = t_r * c_r
               do ir = 1, 3
                  if (ir .eq. ip .or. ir .eq. iq) goto 20
                  g_r = A(ir,ip)
                  h_r = A(ir,iq)
                  A(ir,ip) = c_r*g_r - s_r*h_r
                  A(ip,ir) = A(ir,ip)
                  A(ir,iq) = s_r*g_r + c_r*h_r
                  A(iq,ir) = A(ir,iq)
 20               continue
               end do
               ap_diag = A(ip,ip)
               aq_diag = A(iq,iq)
               A(ip,ip) = c_r*c_r*ap_diag - 2.0d0*c_r*s_r*A(ip,iq)
     &                  + s_r*s_r*aq_diag
               A(iq,iq) = s_r*s_r*ap_diag + 2.0d0*c_r*s_r*A(ip,iq)
     &                  + c_r*c_r*aq_diag
               A(ip,iq) = 0.0d0
               A(iq,ip) = 0.0d0
               do ir = 1, 3
                  g_r    = Q(ir,ip)
                  h_r    = Q(ir,iq)
                  Q(ir,ip) = c_r*g_r - s_r*h_r
                  Q(ir,iq) = s_r*g_r + c_r*h_r
               end do
 10            continue
            end do
         end do
      end do
 100  continue
      eigval(1) = A(1,1)
      eigval(2) = A(2,2)
      eigval(3) = A(3,3)
      return
      end

      subroutine compute_eq_strain_epd(eig_v, E_mod, nu_val, ebar)
      DOUBLE PRECISION eig_v(3), E_mod, nu_val, ebar
      DOUBLE PRECISION I1, J2, Y_ep, s1, s2, s3
      I1 = eig_v(1) + eig_v(2) + eig_v(3)
      s1 = eig_v(1) - I1/3.0d0
      s2 = eig_v(2) - I1/3.0d0
      s3 = eig_v(3) - I1/3.0d0
      J2 = 0.5d0*(s1*s1 + s2*s2 + s3*s3)
      Y_ep = (1.0d0/(2.0d0*E_mod)) * (
     &       2.0d0*(1.0d0+nu_val)*J2 +
     &       (1.0d0-2.0d0*nu_val)/3.0d0 * I1*I1 )
      if (Y_ep .lt. 0.0d0) Y_ep = 0.0d0
      ebar = dsqrt(2.0d0 * Y_ep / E_mod)
      return
      end

      subroutine eval_damage_t_epd(kap_t, Et, st0, et0, stp,
     &     etp_obj, etu_obj, sig_floor, d_t)
      DOUBLE PRECISION kap_t, Et, st0, et0, stp
      DOUBLE PRECISION etp_obj, etu_obj, sig_floor, d_t
      DOUBLE PRECISION sigma_t, H_t, denom
      denom = Et * kap_t
      if (kap_t .le. et0) then
         d_t = 0.0d0
      else if (kap_t .lt. etp_obj) then
         H_t = (stp - st0) / (etp_obj - et0)
         sigma_t = st0 + H_t*(kap_t - et0)
         if (denom .gt. 1.0d-12) then
            d_t = 1.0d0 - sigma_t/denom
         else
            d_t = 0.0d0
         end if
      else if (kap_t .lt. etu_obj) then
         sigma_t = stp*(etu_obj-kap_t)/(etu_obj-etp_obj)
         sigma_t = max(sigma_t, sig_floor)
         if (denom .gt. 1.0d-12) then
            d_t = 1.0d0 - sigma_t/denom
         else
            d_t = 0.0d0
         end if
      else
         if (denom .gt. 1.0d-12) then
            d_t = 1.0d0 - sig_floor/denom
         else
            d_t = 1.0d0 - 1.0d-6
         end if
      end if
      d_t = max(d_t, 0.0d0)
      return
      end

      subroutine eval_damage_c_epd(kap_c, Ec, scp, ecp, scu,
     &     ecu_obj, ncf, d_c)
      DOUBLE PRECISION kap_c, Ec, scp, ecp, scu, ecu_obj, ncf, d_c
      DOUBLE PRECISION xval, vden, sigma_c, elim, denom
      elim  = (2.0d0/3.0d0) * ecp
      denom = Ec * kap_c
      if (kap_c .le. elim) then
         d_c = 0.0d0
      else if (kap_c .le. ecu_obj) then
         xval = kap_c / ecp
         vden = 1.0d0 + (ncf-2.0d0)*xval + xval*xval
         if (dabs(vden) .lt. 1.0d-12) vden = 1.0d-12
         sigma_c = scp * ncf * xval / vden
         sigma_c = max(sigma_c, scu)
         if (denom .gt. 1.0d-12) then
            d_c = 1.0d0 - sigma_c/denom
         else
            d_c = 0.0d0
         end if
      else
         if (denom .gt. 1.0d-12) then
            d_c = 1.0d0 - scu/denom
         else
            d_c = 1.0d0 - 1.0d-6
         end if
      end if
      d_c = max(d_c, 0.0d0)
      return
      end