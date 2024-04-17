#include "MAPL_Generic.h"

program mam_box_model

    AEROSOL_MICROPHYSICS: if (self%microphysics) then

    do j = 1, jm
        do i = 1, im

            amc_t(ncol, :)    = T(i, j, :)              ! temperature at model levels (K)
            amc_pmid(ncol, :) = 0.5*(ple(i,j,0:lm-1)+ple(i,j,1:lm)) ! pressure at layer center (Pa)
            amc_pdel(ncol, :) = delp(i,j,:)             ! pressure thickness of layer (Pa)
            amc_zm(ncol, :)   = zle(i,j,1:lm)           ! altitude (above ground) at layer center (m)
            amc_pblh(ncol)    = zpbl(i,j)               ! planetary boundary layer depth (m)

            amc_qv(ncol, :)   = Q(i,j,:)                ! specific humidity (kg/kg)
            amc_cld(ncol, :)  = fcld(i,j,:)             ! cloud fraction
            amc_rh(ncol, :)   = RH(i,j,:)               ! relative humidity

            ! current tracer mixing ratios (TMRs)

            amc_qqcw(:ncol,:pver,:pcnstxx)            = tiny(0.0)
            amc_qqcw_precldchem(:ncol,:pver,:pcnstxx) = tiny(0.0)  ! qqcw TMRs before cloud chemistry


            ! units mixing ratios should be 'mol/mol-air' and '#/kmol-air'
            amc_q(ncol,:,i_h2o2)   = tiny(0.0)               ! h2o2
            amc_q(ncol,:,i_h2so4)  = h2so4(i,j,:)
            amc_q(ncol,:,i_so2)    = so2(i,j,:)
            amc_q(ncol,:,i_dms)    = tiny(0.0)               ! dms
            amc_q(ncol,:,i_nh3)    = nh3(i,j,:)
            amc_q(ncol,:,i_soag)   = soa_g(i,j,:)
            ! accumulation mode
            amc_q(ncol,:,i_so4_a1) = acc_a_so4(i,j,:) * (mw_air / adv_mass(i_so4_a1))
            amc_q(ncol,:,i_nh4_a1) = acc_a_nh4(i,j,:) * (mw_air / adv_mass(i_nh4_a1))
            amc_q(ncol,:,i_pom_a1) = acc_a_pom(i,j,:) * (mw_air / adv_mass(i_pom_a1))
            amc_q(ncol,:,i_soa_a1) = acc_a_soa(i,j,:) * (mw_air / adv_mass(i_soa_a1))
            amc_q(ncol,:,i_bc_a1)  = acc_a_bc(i,j,:)  * (mw_air / adv_mass(i_bc_a1))
            amc_q(ncol,:,i_ncl_a1) = acc_a_ncl(i,j,:) * (mw_air / adv_mass(i_ncl_a1))
            amc_q(ncol,:,i_num_a1) = acc_a_num(i,j,:) *  mw_air
            ! aitken mode
            amc_q(ncol,:,i_so4_a2) = ait_a_so4(i,j,:) * (mw_air / adv_mass(i_so4_a2))
            amc_q(ncol,:,i_nh4_a2) = ait_a_nh4(i,j,:) * (mw_air / adv_mass(i_nh4_a2))
            amc_q(ncol,:,i_soa_a2) = ait_a_soa(i,j,:) * (mw_air / adv_mass(i_soa_a2))
            amc_q(ncol,:,i_ncl_a2) = ait_a_ncl(i,j,:) * (mw_air / adv_mass(i_ncl_a2))
            amc_q(ncol,:,i_num_a2) = ait_a_num(i,j,:) *  mw_air
            ! primary carbon mode
            amc_q(ncol,:,i_pom_a3) = pcm_a_pom(i,j,:) * (mw_air / adv_mass(i_pom_a3))
            amc_q(ncol,:,i_bc_a3)  = pcm_a_bc(i,j,:)  * (mw_air / adv_mass(i_bc_a3))
            amc_q(ncol,:,i_num_a3) = pcm_a_num(i,j,:) *  mw_air
            ! fine seasalt mode
            amc_q(ncol,:,i_ncl_a4) = fss_a_ncl(i,j,:) * (mw_air / adv_mass(i_ncl_a4))
            amc_q(ncol,:,i_so4_a4) = fss_a_so4(i,j,:) * (mw_air / adv_mass(i_so4_a4))
            amc_q(ncol,:,i_nh4_a4) = fss_a_nh4(i,j,:) * (mw_air / adv_mass(i_nh4_a4))
            amc_q(ncol,:,i_num_a4) = fss_a_num(i,j,:) *  mw_air
            ! fine dust mode
            amc_q(ncol,:,i_dst_a5) = fdu_a_dst(i,j,:) * (mw_air / adv_mass(i_dst_a5))
            amc_q(ncol,:,i_so4_a5) = fdu_a_so4(i,j,:) * (mw_air / adv_mass(i_so4_a5))
            amc_q(ncol,:,i_nh4_a5) = fdu_a_nh4(i,j,:) * (mw_air / adv_mass(i_nh4_a5))
            amc_q(ncol,:,i_num_a5) = fdu_a_num(i,j,:) *  mw_air
            ! coarse seasalt mode
            amc_q(ncol,:,i_ncl_a6) = css_a_ncl(i,j,:) * (mw_air / adv_mass(i_ncl_a6))
            amc_q(ncol,:,i_so4_a6) = css_a_so4(i,j,:) * (mw_air / adv_mass(i_so4_a6))
            amc_q(ncol,:,i_nh4_a6) = css_a_nh4(i,j,:) * (mw_air / adv_mass(i_nh4_a6))
            amc_q(ncol,:,i_num_a6) = css_a_num(i,j,:) *  mw_air
            ! coarse dust mode
            amc_q(ncol,:,i_dst_a7) = cdu_a_dst(i,j,:) * (mw_air / adv_mass(i_dst_a7))
            amc_q(ncol,:,i_so4_a7) = cdu_a_so4(i,j,:) * (mw_air / adv_mass(i_so4_a7))
            amc_q(ncol,:,i_nh4_a7) = cdu_a_nh4(i,j,:) * (mw_air / adv_mass(i_nh4_a7))
            amc_q(ncol,:,i_num_a7) = cdu_a_num(i,j,:) *  mw_air


            amc_q_pregaschem(:ncol,:pver,:pcnstxx) = amc_q      ! q TMRs    before gas-phase chemistry

! DEBUG-SS
if ((i == 1) .and. (j == 1)) then
print*, "The value of H2SO4 in aerosol microphysics in time step",EmCtr,"is",minval(amc_q(:,:,i_h2so4)) , maxval(amc_q(:,:,i_h2so4))
endif
! DEBUG-SS

#if (0)
            ! compute pregaschem using tendencies
            amc_q_pregaschem(ncol,:,i_h2so4) = h2so4(i,j,:) - (ddt_h2so4_gas(i,j,:) + ddt_h2so4_aq(i,j,:))*self%dt
            amc_q_pregaschem(ncol,:,i_so2)   = so2(i,j,:)   - (ddt_so2_gas(i,j,:)   + ddt_so2_aq(i,j,:)  )*self%dt
            amc_q_pregaschem(ncol,:,i_nh3)   = nh3(i,j,:)   - (ddt_nh3_gas(i,j,:)   + ddt_nh3_aq(i,j,:)  )*self%dt
#else
            ! ...or use the pregas exports
            amc_q_pregaschem(ncol,:,i_h2so4) = h2so4_g_(i,j,:)

! DEBUG-SS
if ((i == 1) .and. (j == 1)) then
print *, "The value of H2SO4 at time step", EmCtr, "after copying the pregaschem value is ", minval(amc_q_pregaschem(:,:,i_h2so4)) , maxval(amc_q_pregaschem(:,:,i_h2so4))
endif
! DEBUG-SS

            amc_q_pregaschem(ncol,:,i_so2)   = so2_g_(i,j,:)
            amc_q_pregaschem(ncol,:,i_nh3)   = nh3_g_(i,j,:)
#endif


            amc_q_precldchem(:ncol,:pver,:pcnstxx) = amc_q      ! q TMRs    before cloud chemistry
#if (0)
            ! compute preaqchem using tendencies
            amc_q_precldchem(ncol,:,i_h2so4)  = h2so4(i,j,:) - (ddt_h2so4_aq(i,j,:))*self%dt
            amc_q_precldchem(ncol,:,i_so2)    = so2(i,j,:)   - (ddt_so2_aq(i,j,:)  )*self%dt
            amc_q_precldchem(ncol,:,i_nh3)    = nh3(i,j,:)   - (ddt_nh3_aq(i,j,:)  )*self%dt
#else
            ! ...or use the preaq exports
            amc_q_precldchem(ncol,:,i_h2so4)  = h2so4_a_(i,j,:)
            amc_q_precldchem(ncol,:,i_so2)    = so2_a_(i,j,:)
            amc_q_precldchem(ncol,:,i_nh3)    = nh3_a_(i,j,:)
#endif
            amc_q_precldchem(ncol,:,i_so4_a1) = acc_a_so4_(i,j,:) * (mw_air / adv_mass(i_so4_a1))
            amc_q_precldchem(ncol,:,i_nh4_a1) = acc_a_nh4_(i,j,:) * (mw_air / adv_mass(i_nh4_a1))
            amc_q_precldchem(ncol,:,i_so4_a2) = ait_a_so4_(i,j,:) * (mw_air / adv_mass(i_so4_a2))
            amc_q_precldchem(ncol,:,i_nh4_a2) = ait_a_nh4_(i,j,:) * (mw_air / adv_mass(i_nh4_a2))
            amc_q_precldchem(ncol,:,i_so4_a4) = fss_a_so4_(i,j,:) * (mw_air / adv_mass(i_so4_a4))
            amc_q_precldchem(ncol,:,i_nh4_a4) = fss_a_nh4_(i,j,:) * (mw_air / adv_mass(i_nh4_a4))
            amc_q_precldchem(ncol,:,i_so4_a5) = fdu_a_so4_(i,j,:) * (mw_air / adv_mass(i_so4_a5))
            amc_q_precldchem(ncol,:,i_nh4_a5) = fdu_a_nh4_(i,j,:) * (mw_air / adv_mass(i_nh4_a5))
            amc_q_precldchem(ncol,:,i_so4_a6) = css_a_so4_(i,j,:) * (mw_air / adv_mass(i_so4_a6))
            amc_q_precldchem(ncol,:,i_nh4_a6) = css_a_nh4_(i,j,:) * (mw_air / adv_mass(i_nh4_a6))
            amc_q_precldchem(ncol,:,i_so4_a7) = cdu_a_so4_(i,j,:) * (mw_air / adv_mass(i_so4_a7))
            amc_q_precldchem(ncol,:,i_nh4_a7) = cdu_a_nh4_(i,j,:) * (mw_air / adv_mass(i_nh4_a7))


            amc_dgn_a_dry(pcols,:,1) = acc_dgn_dry(i,j,:)          ! dry geo. mean dia. (m) of number PSD
            amc_dgn_a_dry(pcols,:,2) = ait_dgn_dry(i,j,:)
            amc_dgn_a_dry(pcols,:,3) = pcm_dgn_dry(i,j,:)
            amc_dgn_a_dry(pcols,:,4) = fss_dgn_dry(i,j,:)
            amc_dgn_a_dry(pcols,:,5) = fdu_dgn_dry(i,j,:)
            amc_dgn_a_dry(pcols,:,6) = css_dgn_dry(i,j,:)
            amc_dgn_a_dry(pcols,:,7) = cdu_dgn_dry(i,j,:)

            amc_dgn_a_wet(pcols,:,1) = acc_dgn_wet(i,j,:)          ! wet geo. mean dia. (m) of number PSD
            amc_dgn_a_wet(pcols,:,2) = ait_dgn_wet(i,j,:)
            amc_dgn_a_wet(pcols,:,3) = pcm_dgn_wet(i,j,:)
            amc_dgn_a_wet(pcols,:,4) = fss_dgn_wet(i,j,:)
            amc_dgn_a_wet(pcols,:,5) = fdu_dgn_wet(i,j,:)
            amc_dgn_a_wet(pcols,:,6) = css_dgn_wet(i,j,:)
            amc_dgn_a_wet(pcols,:,7) = cdu_dgn_wet(i,j,:)


            amc_wetdens_host(:pcols,:pver,:ntot_amode) = 1.0e3     ! interstitial aerosol wet density (kg/m3)
            amc_qaerwat(:pcols,:pver,:ntot_amode)      = 0.0       ! optional, aerosol water mixing ratio (kg/kg)
            amc_qaerwat(pcols,:,1) = acc_a_wtr(i,j,:)              ! aerosol water
            amc_qaerwat(pcols,:,2) = ait_a_wtr(i,j,:)
            amc_qaerwat(pcols,:,3) = pcm_a_wtr(i,j,:)
            amc_qaerwat(pcols,:,4) = fss_a_wtr(i,j,:)
            amc_qaerwat(pcols,:,5) = fdu_a_wtr(i,j,:)
            amc_qaerwat(pcols,:,6) = css_a_wtr(i,j,:)
            amc_qaerwat(pcols,:,7) = cdu_a_wtr(i,j,:)

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

! DEBUG - SS
if ((i == 1) .and. (j == 1)) then
print*, "The value of temperature is", minval(amc_t), maxval(amc_t)
endif
! DEBUG - SS

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
!           h2o2             = amc_q(ncol,:,i_h2o2)
            h2so4(i,j,:)     = amc_q(ncol,:,i_h2so4)

! DEBUG-SS
if ((i == 1) .and. (j == 1)) then
print*, "The H2SO4 after modal_aero_amicphys_intr in time step",EmCtr,"is",minval(h2so4) , maxval(h2so4)
endif
! DEBUG-SS

            so2(i,j,:)       = amc_q(ncol,:,i_so2)
!           dms              = amc_q(ncol,:,i_dms)
            nh3(i,j,:)       = amc_q(ncol,:,i_nh3)
            soa_g(i,j,:)     = amc_q(ncol,:,i_soag)
            ! accumulation mode
            acc_a_so4(i,j,:) = amc_q(ncol,:,i_so4_a1) * (adv_mass(i_so4_a1) / mw_air)
            acc_a_nh4(i,j,:) = amc_q(ncol,:,i_nh4_a1) * (adv_mass(i_nh4_a1) / mw_air)
            acc_a_pom(i,j,:) = amc_q(ncol,:,i_pom_a1) * (adv_mass(i_pom_a1) / mw_air)
            acc_a_soa(i,j,:) = amc_q(ncol,:,i_soa_a1) * (adv_mass(i_soa_a1) / mw_air)
            acc_a_bc(i,j,:)  = amc_q(ncol,:,i_bc_a1)  * (adv_mass(i_bc_a1)  / mw_air)
            acc_a_ncl(i,j,:) = amc_q(ncol,:,i_ncl_a1) * (adv_mass(i_ncl_a1) / mw_air)
            acc_a_num(i,j,:) = amc_q(ncol,:,i_num_a1) / mw_air
            ! aitken mode
            ait_a_so4(i,j,:) = amc_q(ncol,:,i_so4_a2) * (adv_mass(i_so4_a2) / mw_air)
            ait_a_nh4(i,j,:) = amc_q(ncol,:,i_nh4_a2) * (adv_mass(i_nh4_a2) / mw_air)
            ait_a_soa(i,j,:) = amc_q(ncol,:,i_soa_a2) * (adv_mass(i_soa_a2) / mw_air)
            ait_a_ncl(i,j,:) = amc_q(ncol,:,i_ncl_a2) * (adv_mass(i_ncl_a2) / mw_air)
            ait_a_num(i,j,:) = amc_q(ncol,:,i_num_a2) / mw_air
            ! primary carbon mode
            pcm_a_pom(i,j,:) = amc_q(ncol,:,i_pom_a3) * (adv_mass(i_pom_a3) / mw_air)
            pcm_a_bc(i,j,:)  = amc_q(ncol,:,i_bc_a3)  * (adv_mass(i_bc_a3)  / mw_air)
            pcm_a_num(i,j,:) = amc_q(ncol,:,i_num_a3) / mw_air
            ! fine seasalt mode
            fss_a_ncl(i,j,:) = amc_q(ncol,:,i_ncl_a4) * (adv_mass(i_ncl_a4) / mw_air)
            fss_a_so4(i,j,:) = amc_q(ncol,:,i_so4_a4) * (adv_mass(i_so4_a4) / mw_air)
            fss_a_nh4(i,j,:) = amc_q(ncol,:,i_nh4_a4) * (adv_mass(i_nh4_a4) / mw_air)
            fss_a_num(i,j,:) = amc_q(ncol,:,i_num_a4) / mw_air
            ! fine dust mode
            fdu_a_dst(i,j,:) = amc_q(ncol,:,i_dst_a5) * (adv_mass(i_dst_a5) / mw_air)
            fdu_a_so4(i,j,:) = amc_q(ncol,:,i_so4_a5) * (adv_mass(i_so4_a5) / mw_air)
            fdu_a_nh4(i,j,:) = amc_q(ncol,:,i_nh4_a5) * (adv_mass(i_nh4_a5) / mw_air)
            fdu_a_num(i,j,:) = amc_q(ncol,:,i_num_a5) / mw_air
            ! coarse seasalt mode
            css_a_ncl(i,j,:) = amc_q(ncol,:,i_ncl_a6) * (adv_mass(i_ncl_a6) / mw_air)
            css_a_so4(i,j,:) = amc_q(ncol,:,i_so4_a6) * (adv_mass(i_so4_a6) / mw_air)
            css_a_nh4(i,j,:) = amc_q(ncol,:,i_nh4_a6) * (adv_mass(i_nh4_a6) / mw_air)
            css_a_num(i,j,:) = amc_q(ncol,:,i_num_a6) / mw_air
            ! coarse dust mode
            cdu_a_dst(i,j,:) = amc_q(ncol,:,i_dst_a7) * (adv_mass(i_dst_a7) / mw_air)
            cdu_a_so4(i,j,:) = amc_q(ncol,:,i_so4_a7) * (adv_mass(i_so4_a7) / mw_air)
            cdu_a_nh4(i,j,:) = amc_q(ncol,:,i_nh4_a7) * (adv_mass(i_nh4_a7) / mw_air)
            cdu_a_num(i,j,:) = amc_q(ncol,:,i_num_a7) / mw_air

            ! save the colmn-integrated diagnostics
            q_coltend_cond_(i,j,:)   = amc_q_coltendaa(ncol,:,iqtend_cond)
            q_coltend_rename_(i,j,:) = amc_q_coltendaa(ncol,:,iqtend_rnam)
            q_coltend_nucl_(i,j,:)   = amc_q_coltendaa(ncol,:,iqtend_nnuc)
            q_coltend_coag_(i,j,:)   = amc_q_coltendaa(ncol,:,iqtend_coag)
            qqcw_coltend_rename_(i,j,:) = amc_q_coltendaa(ncol,:,iqqcwtend_rnam)
        end do
    end do

    end if AEROSOL_MICROPHYSICS

end program mam_box_model
