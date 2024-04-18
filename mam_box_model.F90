#include "MAPL_Generic.h"

program mam_box_model

            amc_t    = T              ! temperature at model levels (K)
            amc_pmid = P1              ! pressure at layer center (Pa)
            amc_pdel = P2             ! pressure thickness of layer (Pa)
            amc_zm   = zle           ! altitude (above ground) at layer center (m)
            amc_pblh = zpbl               ! planetary boundary layer depth (m)
            amc_qv   = Q                ! specific humidity (kg/kg)
            amc_cld  = fcld             ! cloud fraction
            amc_rh   = RH               ! relative humidity

            ! current tracer mixing ratios (TMRs)

            !pcnstxx is from gas_pcnst in chem_mods. It is the number of gas-phase species = 37
            amc_qqcw(:pcnstxx)            = tiny(0.0)
            amc_qqcw_precldchem(:pcnstxx) = tiny(0.0)  ! qqcw TMRs before cloud chemistry


            ! units mixing ratios should be 'mol/mol-air' and '#/kmol-air'
            amc_q(i_h2o2)   = tiny(0.0)               ! h2o2
            amc_q(i_h2so4)  = h2so4
            amc_q(i_so2)    = so2
            amc_q(i_dms)    = tiny(0.0)               ! dms
            amc_q(i_nh3)    = nh3
            amc_q(i_soag)   = soa_g
            ! accumulation mode
            amc_q(i_so4_a1) = acc_a_so4 * (mw_air / adv_mass(i_so4_a1))
            amc_q(i_nh4_a1) = acc_a_nh4 * (mw_air / adv_mass(i_nh4_a1))
            amc_q(i_pom_a1) = acc_a_pom * (mw_air / adv_mass(i_pom_a1))
            amc_q(i_soa_a1) = acc_a_soa * (mw_air / adv_mass(i_soa_a1))
            amc_q(i_bc_a1)  = acc_a_bc  * (mw_air / adv_mass(i_bc_a1))
            amc_q(i_ncl_a1) = acc_a_ncl * (mw_air / adv_mass(i_ncl_a1))
            amc_q(i_num_a1) = acc_a_num *  mw_air
            ! aitken mode
            amc_q(i_so4_a2) = ait_a_so4 * (mw_air / adv_mass(i_so4_a2))
            amc_q(i_nh4_a2) = ait_a_nh4 * (mw_air / adv_mass(i_nh4_a2))
            amc_q(i_soa_a2) = ait_a_soa * (mw_air / adv_mass(i_soa_a2))
            amc_q(i_ncl_a2) = ait_a_ncl * (mw_air / adv_mass(i_ncl_a2))
            amc_q(i_num_a2) = ait_a_num *  mw_air
            ! primary carbon mode
            amc_q(i_pom_a3) = pcm_a_pom * (mw_air / adv_mass(i_pom_a3))
            amc_q(i_bc_a3)  = pcm_a_bc  * (mw_air / adv_mass(i_bc_a3))
            amc_q(i_num_a3) = pcm_a_num *  mw_air
            ! fine seasalt mode
            amc_q(i_ncl_a4) = fss_a_ncl * (mw_air / adv_mass(i_ncl_a4))
            amc_q(i_so4_a4) = fss_a_so4 * (mw_air / adv_mass(i_so4_a4))
            amc_q(i_nh4_a4) = fss_a_nh4 * (mw_air / adv_mass(i_nh4_a4))
            amc_q(i_num_a4) = fss_a_num *  mw_air
            ! fine dust mode
            amc_q(i_dst_a5) = fdu_a_dst * (mw_air / adv_mass(i_dst_a5))
            amc_q(i_so4_a5) = fdu_a_so4 * (mw_air / adv_mass(i_so4_a5))
            amc_q(i_nh4_a5) = fdu_a_nh4 * (mw_air / adv_mass(i_nh4_a5))
            amc_q(i_num_a5) = fdu_a_num *  mw_air
            ! coarse seasalt mode
            amc_q(i_ncl_a6) = css_a_ncl * (mw_air / adv_mass(i_ncl_a6))
            amc_q(i_so4_a6) = css_a_so4 * (mw_air / adv_mass(i_so4_a6))
            amc_q(i_nh4_a6) = css_a_nh4 * (mw_air / adv_mass(i_nh4_a6))
            amc_q(i_num_a6) = css_a_num *  mw_air
            ! coarse dust mode
            amc_q(i_dst_a7) = cdu_a_dst * (mw_air / adv_mass(i_dst_a7))
            amc_q(i_so4_a7) = cdu_a_so4 * (mw_air / adv_mass(i_so4_a7))
            amc_q(i_nh4_a7) = cdu_a_nh4 * (mw_air / adv_mass(i_nh4_a7))
            amc_q(i_num_a7) = cdu_a_num *  mw_air


            amc_q_pregaschem(:pcnstxx) = amc_q      ! q TMRs    before gas-phase chemistry

            ! ...or use the pregas exports
            amc_q_pregaschem(i_h2so4) = h2so4_g_

            amc_q_pregaschem(i_so2)   = so2_g_
            amc_q_pregaschem(i_nh3)   = nh3_g_


            amc_q_precldchem(:pcnstxx) = amc_q      ! q TMRs    before cloud chemistry
            
            ! ...or use the preaq exports
            amc_q_precldchem(i_h2so4)  = h2so4_a_
            amc_q_precldchem(i_so2)    = so2_a_
            amc_q_precldchem(i_nh3)    = nh3_a_
            
            amc_q_precldchem(i_so4_a1) = acc_a_so4_ * (mw_air / adv_mass(i_so4_a1))
            amc_q_precldchem(i_nh4_a1) = acc_a_nh4_ * (mw_air / adv_mass(i_nh4_a1))
            amc_q_precldchem(i_so4_a2) = ait_a_so4_ * (mw_air / adv_mass(i_so4_a2))
            amc_q_precldchem(i_nh4_a2) = ait_a_nh4_ * (mw_air / adv_mass(i_nh4_a2))
            amc_q_precldchem(i_so4_a4) = fss_a_so4_ * (mw_air / adv_mass(i_so4_a4))
            amc_q_precldchem(i_nh4_a4) = fss_a_nh4_ * (mw_air / adv_mass(i_nh4_a4))
            amc_q_precldchem(i_so4_a5) = fdu_a_so4_ * (mw_air / adv_mass(i_so4_a5))
            amc_q_precldchem(i_nh4_a5) = fdu_a_nh4_ * (mw_air / adv_mass(i_nh4_a5))
            amc_q_precldchem(i_so4_a6) = css_a_so4_ * (mw_air / adv_mass(i_so4_a6))
            amc_q_precldchem(i_nh4_a6) = css_a_nh4_ * (mw_air / adv_mass(i_nh4_a6))
            amc_q_precldchem(i_so4_a7) = cdu_a_so4_ * (mw_air / adv_mass(i_so4_a7))
            amc_q_precldchem(i_nh4_a7) = cdu_a_nh4_ * (mw_air / adv_mass(i_nh4_a7))


            amc_dgn_a_dry(1) = acc_dgn_dry          ! dry geo. mean dia. (m) of number PSD
            amc_dgn_a_dry(2) = ait_dgn_dry
            amc_dgn_a_dry(3) = pcm_dgn_dry
            amc_dgn_a_dry(4) = fss_dgn_dry
            amc_dgn_a_dry(5) = fdu_dgn_dry
            amc_dgn_a_dry(6) = css_dgn_dry
            amc_dgn_a_dry(7) = cdu_dgn_dry

            amc_dgn_a_wet(1) = acc_dgn_wet          ! wet geo. mean dia. (m) of number PSD
            amc_dgn_a_wet(2) = ait_dgn_wet
            amc_dgn_a_wet(3) = pcm_dgn_wet
            amc_dgn_a_wet(4) = fss_dgn_wet
            amc_dgn_a_wet(5) = fdu_dgn_wet
            amc_dgn_a_wet(6) = css_dgn_wet
            amc_dgn_a_wet(7) = cdu_dgn_wet


            amc_wetdens_host(:ntot_amode) = 1.0e3     ! interstitial aerosol wet density (kg/m3)
            amc_qaerwat(:ntot_amode)      = 0.0       ! optional, aerosol water mixing ratio (kg/kg)
            amc_qaerwat(1) = acc_a_wtr              ! aerosol water
            amc_qaerwat(2) = ait_a_wtr
            amc_qaerwat(3) = pcm_a_wtr
            amc_qaerwat(4) = fss_a_wtr
            amc_qaerwat(5) = fdu_a_wtr
            amc_qaerwat(6) = css_a_wtr
            amc_qaerwat(7) = cdu_a_wtr

            amc_q_coltendaa        = 0.0d0                         ! column integrated tendencies diagnostics
            amc_qqcw_coltendaa     = 0.0d0                         ! --dito-- but for qqcw


            ! the modal_aero_amicphys_intr() subroutine does in the order listed below:
            !
            ! - in clear grid cells
            !    1. condensation / gas-aerosol-exchange of H2SO4, NH3, H2O
            !    2. renaming after "continuous growth"
            !    3. nucleation (new particle formation)
            !    4. coagulation
            !    5. primary carbon aging
            !
            ! - in cloudy grid cells
            !    1. condensation / gas-aerosol-exchange
            !    2. renaming after "continuous growth"
            !    3. primary carbon aging

            call modal_aero_amicphys_intr(amc_do_gasaerexch,   &
                                          amc_do_rename,       &
                                          amc_do_newnuc,       &
                                          amc_do_coag,         &
                                          amc_lchnk,           &
                                          ncol,                &
                                          amc_nstep,           &
                                          amc_loffset,         &
                                          amc_deltat,          &
                                          amc_latndx,          &
                                          amc_lonndx,          &
                                          amc_t,               &
                                          amc_pmid,            &
                                          amc_pdel,            &
                                          amc_zm,              &
                                          amc_pblh,            &
                                          amc_qv,              &
                                          amc_cld,             &
                                          amc_rh,              &
                                          amc_q,               &
                                          amc_qqcw,            &
                                          amc_q_pregaschem,    &
                                          amc_q_precldchem,    &
                                          amc_qqcw_precldchem, &
                                          amc_dgn_a_dry,       &
                                          amc_dgn_a_wet,       &
                                          amc_wetdens_host,    &
                                          amc_q_coltendaa,     &
                                          amc_qqcw_coltendaa)! & amc_qaerwat -- optional)


            ! current tracer mixing ratios (TMRs)
!           h2o2      = amc_q(i_h2o2)
            h2so4     = amc_q(i_h2so4)
            so2       = amc_q(i_so2)
!           dms       = amc_q(i_dms)
            nh3       = amc_q(i_nh3)
            soa_g     = amc_q(i_soag)
            ! accumulation mode
            acc_a_so4 = amc_q(i_so4_a1) * (adv_mass(i_so4_a1) / mw_air)
            acc_a_nh4 = amc_q(i_nh4_a1) * (adv_mass(i_nh4_a1) / mw_air)
            acc_a_pom = amc_q(i_pom_a1) * (adv_mass(i_pom_a1) / mw_air)
            acc_a_soa = amc_q(i_soa_a1) * (adv_mass(i_soa_a1) / mw_air)
            acc_a_bc  = amc_q(i_bc_a1)  * (adv_mass(i_bc_a1)  / mw_air)
            acc_a_ncl = amc_q(i_ncl_a1) * (adv_mass(i_ncl_a1) / mw_air)
            acc_a_num = amc_q(i_num_a1) / mw_air
            ! aitken mode
            ait_a_so4 = amc_q(i_so4_a2) * (adv_mass(i_so4_a2) / mw_air)
            ait_a_nh4 = amc_q(i_nh4_a2) * (adv_mass(i_nh4_a2) / mw_air)
            ait_a_soa = amc_q(i_soa_a2) * (adv_mass(i_soa_a2) / mw_air)
            ait_a_ncl = amc_q(i_ncl_a2) * (adv_mass(i_ncl_a2) / mw_air)
            ait_a_num = amc_q(i_num_a2) / mw_air
            ! primary carbon mode
            pcm_a_pom = amc_q(i_pom_a3) * (adv_mass(i_pom_a3) / mw_air)
            pcm_a_bc  = amc_q(i_bc_a3)  * (adv_mass(i_bc_a3)  / mw_air)
            pcm_a_num = amc_q(i_num_a3) / mw_air
            ! fine seasalt mode
            fss_a_ncl = amc_q(i_ncl_a4) * (adv_mass(i_ncl_a4) / mw_air)
            fss_a_so4 = amc_q(i_so4_a4) * (adv_mass(i_so4_a4) / mw_air)
            fss_a_nh4 = amc_q(i_nh4_a4) * (adv_mass(i_nh4_a4) / mw_air)
            fss_a_num = amc_q(i_num_a4) / mw_air
            ! fine dust mode
            fdu_a_dst = amc_q(i_dst_a5) * (adv_mass(i_dst_a5) / mw_air)
            fdu_a_so4 = amc_q(i_so4_a5) * (adv_mass(i_so4_a5) / mw_air)
            fdu_a_nh4 = amc_q(i_nh4_a5) * (adv_mass(i_nh4_a5) / mw_air)
            fdu_a_num = amc_q(i_num_a5) / mw_air
            ! coarse seasalt mode
            css_a_ncl = amc_q(i_ncl_a6) * (adv_mass(i_ncl_a6) / mw_air)
            css_a_so4 = amc_q(i_so4_a6) * (adv_mass(i_so4_a6) / mw_air)
            css_a_nh4 = amc_q(i_nh4_a6) * (adv_mass(i_nh4_a6) / mw_air)
            css_a_num = amc_q(i_num_a6) / mw_air
            ! coarse dust mode
            cdu_a_dst = amc_q(i_dst_a7) * (adv_mass(i_dst_a7) / mw_air)
            cdu_a_so4 = amc_q(i_so4_a7) * (adv_mass(i_so4_a7) / mw_air)
            cdu_a_nh4 = amc_q(i_nh4_a7) * (adv_mass(i_nh4_a7) / mw_air)
            cdu_a_num = amc_q(i_num_a7) / mw_air

            ! save the colmn-integrated diagnostics
            q_coltend_cond_        = amc_q_coltendaa(iqtend_cond)
            q_coltend_rename_      = amc_q_coltendaa(iqtend_rnam)
            q_coltend_nucl_        = amc_q_coltendaa(iqtend_nnuc)
            q_coltend_coag_        = amc_q_coltendaa(iqtend_coag)
            qqcw_coltend_rename_   = amc_q_coltendaa(iqqcwtend_rnam)
        end do
    end do

    end if AEROSOL_MICROPHYSICS

end program mam_box_model
