% m.file for computing perturbations
%
% for svm_mgen.m, step 3a: network solution
% forms state space model of system

%-----------------------------------------------------------------------------%
% Version history
%
% Version: 2.0
% Date:    2019
% Purpose: Implemented the central difference method. Also made changes
%          required to make svm_mgen consistent with s_simu.
% Author:  Ryan Elliott
%
% Version: 1.0 (initial version)
% Date:    1997
% Author:  Graham Rogers
%-----------------------------------------------------------------------------%

bus_start = bus;  % power flow solution at initial operating point

for kl = 2:4
    bus_sim = bus_start;

    ks = 1;    % required for i_simu call
    flag = 1;  % network solution

    % network-machine interface for generators
    mac_ind(0,kl,bus_sim,flag);
    mac_igen(0,kl,bus_sim,flag);
    mac_sub(0,kl,bus_sim,flag);
    mac_tra(0,kl,bus_sim,flag);
    mac_em(0,kl,bus_sim,flag);
    gfma(0,kl,bus_sim,flag);
    mac_ivm(0,kl,bus_sim,flag);
    mac_ib(0,kl,bus_sim,flag);

    dc_cont(0,kl,10*(kl-1)+1,bus_sim,flag);

    % calling i_simu using the pre-fault system matrices
    % calculates current injections and bus voltages and angles
    line_sim = line;
    g.bus.bus_int = bus_intprf;
    Y1 = Y_gprf;
    Y2 = Y_gncprf;
    Y3 = Y_ncgprf;
    Y4 = Y_ncprf;
    Vr1 = V_rgprf;
    Vr2 = V_rncprf;
    bo = boprf;

    h_sol = i_simu(kl,ks,k_inc,h,bus_sim,Y1,Y2,Y3,Y4,Vr1,Vr2,bo);

    % network interface for generator control models
    dc_cont(0,kl,10*(kl-1)+1,bus_sim,flag);
    dpwf(0,kl,flag);
    pss(0,kl,flag);
    smpexc(0,kl,flag);
    smppi(0,kl,flag);
    exc_st3(0,kl,flag);
    exc_dc12(0,kl,flag);
    tg(0,kl,flag);
    tg_hydro(0,kl,flag);

    % calculate rates of change (dynamics)
    flag = 2;

    if (g.mac.n_mac ~= 0)
        mac_ind(0,kl,bus_sim,flag);
        mac_igen(0,kl,bus_sim,flag);
        mac_sub(0,kl,bus_sim,flag);
        mac_tra(0,kl,bus_sim,flag);
        mac_em(0,kl,bus_sim,flag);
        mac_ivm(0,kl,bus_sim,flag);

        dpwf(0,kl,flag);
        pss(0,kl,flag);
        smpexc(0,kl,flag);
        smppi(0,kl,flag);
        exc_st3(0,kl,flag);
        exc_dc12(0,kl,flag);
        tg(0,kl,flag);
        tg_hydro(0,kl,flag);
    end

    if (g.svc.n_svc ~= 0)
        v_svc = abs(g.bus.bus_v(g.bus.bus_int(g.svc.svc_con(:,2)),kl));
        svc(0,kl,bus_sim,flag,v_svc);
    end

    if (g.tcsc.n_tcsc ~= 0)
        tcsc(0,kl,flag);
    end

    if (g.lmod.n_lmod ~= 0)
        lmod(0,kl,flag);
    end

    if (g.rlmod.n_rlmod ~= 0)
        rlmod(0,kl,flag);
    end

    if (g.pwr.n_pwrmod ~= 0)
        pwrmod_p(0,kl,bus_sim,flag);
        pwrmod_q(0,kl,bus_sim,flag);
    end

    if (g.dc.n_conv ~= 0)
        dc_cont(0,kl,kl,bus_sim,flag);
        dc_line(0,kl,flag);
    end

    if (g.ess.n_ess ~= 0)
        lsc(0,kl,bus_sim,flag);
        reec(0,kl,bus_sim,flag,h_sol);
        ess(0,kl,bus_sim,flag,h_sol);
    end

    if (g.ivm.n_ivm ~= 0)
        gfma(0,kl,bus_sim,flag);
    end

    cur = g.mac.cur_re(1:g.mac.n_mac,kl) + 1j*g.mac.cur_im(1:g.mac.n_mac,kl);
    cur_mag(1:g.mac.n_mac,kl) = abs(cur(1:g.mac.n_mac,1)).*g.mac.mac_pot(:,1);

    g.mac.telect(:,kl) = g.mac.pelect(:,kl).*g.mac.mac_pot(:,1) ...
                         + cur_mag(:,kl).*cur_mag(:,kl).*g.mac.mac_con(:,5);

end

%-----------------------------------------------------------------------------%
% form vector of state derivatives

mac_state = nss.mac_max*g.mac.n_mac;
exc_state = mac_state + nss.exc_max*g.exc.n_exc;
pss_state = exc_state + nss.pss_max*g.pss.n_pss;
dpw_state = pss_state + nss.dpw_max*g.dpw.n_dpw;
tg_state = dpw_state + nss.tg_max*g.tg.n_tg_tot;
mot_state = tg_state + nss.mot_max*g.ind.n_mot;
ig_state = mot_state + nss.ig_max*g.igen.n_ig;
svc_state = ig_state + nss.svc_max*g.svc.n_svc;
tcsc_state = svc_state + nss.tcsc_max*g.tcsc.n_tcsc;
lmod_state = tcsc_state + nss.lmod_max*g.lmod.n_lmod;
rlmod_state = lmod_state + nss.rlmod_max*g.rlmod.n_rlmod;
pwrmod_p_state = rlmod_state + nss.pwrmod_max*g.pwr.n_pwrmod;
pwrmod_q_state = pwrmod_p_state + nss.pwrmod_max*g.pwr.n_pwrmod;
dcl_state = pwrmod_q_state + nss.dcl_max*g.dc.n_dcl;
ess_state = dcl_state + nss.ess_max*g.ess.n_ess;
lsc_state = ess_state + nss.lsc_max*g.lsc.n_lsc;
reec_state = lsc_state + nss.reec_max*g.reec.n_reec;
% gfma_state = reec_state + nss.gfma_max*g.gfma.n_gfma;  % for future expansion

d_vector = zeros(max_state,2);
for kl = 3:4
    if (g.mac.n_mac ~= 0)
        d_vector(0*g.mac.n_mac+1:1*g.mac.n_mac,kl-2) = ...
            g.mac.dmac_ang(:,kl);
        d_vector(1*g.mac.n_mac+1:2*g.mac.n_mac,kl-2) = ...
            g.mac.dmac_spd(:,kl);
        d_vector(2*g.mac.n_mac+1:3*g.mac.n_mac,kl-2) = ...
            g.mac.deqprime(:,kl);
        d_vector(3*g.mac.n_mac+1:4*g.mac.n_mac,kl-2) = ...
            g.mac.dpsikd(:,kl);
        d_vector(4*g.mac.n_mac+1:5*g.mac.n_mac,kl-2) = ...
            g.mac.dedprime(:,kl);
        d_vector(5*g.mac.n_mac+1:6*g.mac.n_mac,kl-2) = ...
            g.mac.dpsikq(:,kl);  % 6
    end

    if (g.exc.n_exc ~= 0)
        d_vector(mac_state+0*g.exc.n_exc+1:mac_state+1*g.exc.n_exc,kl-2) = ...
            g.exc.dV_TR(:,kl);
        d_vector(mac_state+1*g.exc.n_exc+1:mac_state+2*g.exc.n_exc,kl-2) = ...
            g.exc.dV_As(:,kl);
        d_vector(mac_state+2*g.exc.n_exc+1:mac_state+3*g.exc.n_exc,kl-2) = ...
            g.exc.dV_R(:,kl);
        d_vector(mac_state+3*g.exc.n_exc+1:mac_state+4*g.exc.n_exc,kl-2) = ...
            g.exc.dEfd(:,kl);
        d_vector(mac_state+4*g.exc.n_exc+1:mac_state+5*g.exc.n_exc,kl-2) = ...
            g.exc.dR_f(:,kl);  % 11
    end

    if (g.pss.n_pss ~= 0)
        d_vector(exc_state+0*g.pss.n_pss+1:exc_state+1*g.pss.n_pss,kl-2) = ...
            g.pss.dpss1(:,kl);
        d_vector(exc_state+1*g.pss.n_pss+1:exc_state+2*g.pss.n_pss,kl-2) = ...
            g.pss.dpss2(:,kl);
        d_vector(exc_state+2*g.pss.n_pss+1:exc_state+3*g.pss.n_pss,kl-2) = ...
            g.pss.dpss3(:,kl);  % 14
    end

    if (g.dpw.n_dpw ~= 0)
        d_vector(pss_state+0*g.dpw.n_dpw+1:pss_state+1*g.dpw.n_dpw,kl-2) = ...
            g.dpw.dsdpw1(:,kl);
        d_vector(pss_state+1*g.dpw.n_dpw+1:pss_state+2*g.dpw.n_dpw,kl-2) = ...
            g.dpw.dsdpw2(:,kl);
        d_vector(pss_state+2*g.dpw.n_dpw+1:pss_state+3*g.dpw.n_dpw,kl-2) = ...
            g.dpw.dsdpw3(:,kl);
        d_vector(pss_state+3*g.dpw.n_dpw+1:pss_state+4*g.dpw.n_dpw,kl-2) = ...
            g.dpw.dsdpw4(:,kl);
        d_vector(pss_state+4*g.dpw.n_dpw+1:pss_state+5*g.dpw.n_dpw,kl-2) = ...
            g.dpw.dsdpw5(:,kl);
        d_vector(pss_state+5*g.dpw.n_dpw+1:pss_state+6*g.dpw.n_dpw,kl-2) = ...
            g.dpw.dsdpw6(:,kl);  % 20
    end

    if (g.tg.n_tg_tot ~= 0)
        d_vector(dpw_state+0*g.tg.n_tg_tot+1:dpw_state+1*g.tg.n_tg_tot,kl-2) = ...
            g.tg.dtg1(:,kl);
        d_vector(dpw_state+1*g.tg.n_tg_tot+1:dpw_state+2*g.tg.n_tg_tot,kl-2) = ...
            g.tg.dtg2(:,kl);
        d_vector(dpw_state+2*g.tg.n_tg_tot+1:dpw_state+3*g.tg.n_tg_tot,kl-2) = ...
            g.tg.dtg3(:,kl);
        d_vector(dpw_state+3*g.tg.n_tg_tot+1:dpw_state+4*g.tg.n_tg_tot,kl-2) = ...
            g.tg.dtg4(:,kl);
        d_vector(dpw_state+4*g.tg.n_tg_tot+1:dpw_state+5*g.tg.n_tg_tot,kl-2) = ...
            g.tg.dtg5(:,kl);  % 25
    end

    if (g.ind.n_mot ~= 0)
        d_vector(tg_state+0*g.ind.n_mot+1:tg_state+1*g.ind.n_mot,kl-2) = ...
            g.ind.dvdp(:,kl);
        d_vector(tg_state+1*g.ind.n_mot+1:tg_state+2*g.ind.n_mot,kl-2) = ...
            g.ind.dvqp(:,kl);
        d_vector(tg_state+2*g.ind.n_mot+1:tg_state+3*g.ind.n_mot,kl-2) = ...
            g.ind.dslip(:,kl);  % 28
    end

    if (g.igen.n_ig ~= 0)
        d_vector(mot_state+0*g.igen.n_ig+1:mot_state+1*g.igen.n_ig,kl-2) = ...
            g.igen.dvdpig(:,kl);
        d_vector(mot_state+1*g.igen.n_ig+1:mot_state+2*g.igen.n_ig,kl-2) = ...
            g.igen.dvqpig(:,kl);
        d_vector(mot_state+2*g.igen.n_ig+1:mot_state+3*g.igen.n_ig,kl-2) = ...
            g.igen.dslig(:,kl);  % 31
    end

    if (g.svc.n_svc ~= 0)
        d_vector(ig_state+0*g.svc.n_svc+1:ig_state+1*g.svc.n_svc,kl-2) = ...
            g.svc.dB_cv(:,kl);
        d_vector(ig_state+1*g.svc.n_svc+1:ig_state+2*g.svc.n_svc,kl-2) = ...
            g.svc.dB_con(:,kl);  % 33
    end

    if (g.tcsc.n_tcsc ~= 0)
        d_vector(svc_state+1:svc_state+g.tcsc.n_tcsc,kl-2) = ...
            g.tcsc.dB_tcsc(:,kl);  % 34
    end

    if (g.lmod.n_lmod ~= 0)
        d_vector(tcsc_state+1:tcsc_state+g.lmod.n_lmod,kl-2) = ...
            g.lmod.dlmod_st(:,kl);  % 35
    end

    if (g.rlmod.n_rlmod ~= 0)
        d_vector(lmod_state+1:lmod_state+g.rlmod.n_rlmod,kl-2) = ...
            g.rlmod.drlmod_st(:,kl);  % 36
    end

    if (g.pwr.n_pwrmod ~= 0)
        d_vector(rlmod_state+1:rlmod_state+g.pwr.n_pwrmod,kl-2) = ...
            g.pwr.dpwrmod_p_st(:,kl);
        d_vector(pwrmod_p_state+1:pwrmod_p_state+g.pwr.n_pwrmod,kl-2) = ...
            g.pwr.dpwrmod_q_st(:,kl);  % 38
    end

    if (g.dc.n_conv ~= 0)
        d_vector(pwrmod_q_state+0*g.dc.n_dcl+1:pwrmod_q_state+1*g.dc.n_dcl,kl-2) = ...
            g.dc.dv_conr(:,kl);
        d_vector(pwrmod_q_state+1*g.dc.n_dcl+1:pwrmod_q_state+2*g.dc.n_dcl,kl-2) = ...
            g.dc.dv_coni(:,kl);
        d_vector(pwrmod_q_state+2*g.dc.n_dcl+1:pwrmod_q_state+3*g.dc.n_dcl,kl-2) = ...
            g.dc.di_dcr(:,kl);
        d_vector(pwrmod_q_state+3*g.dc.n_dcl+1:pwrmod_q_state+4*g.dc.n_dcl,kl-2) = ...
            g.dc.di_dci(:,kl);
        d_vector(pwrmod_q_state+4*g.dc.n_dcl+1:pwrmod_q_state+5*g.dc.n_dcl,kl-2) = ...
            g.dc.dv_dcc(:,kl);  % 43
    end

    if (g.ess.n_ess ~= 0)
        % note: the ess soc has no effect on linearization, so its
        %       state derivatives are neglected.
        d_vector(dcl_state+0*g.ess.n_ess+1:dcl_state+1*g.ess.n_ess,kl-2) = ...
            g.ess.dess1(:,kl);
        d_vector(dcl_state+1*g.ess.n_ess+1:dcl_state+2*g.ess.n_ess,kl-2) = ...
            g.ess.dess2(:,kl);
        d_vector(dcl_state+2*g.ess.n_ess+1:dcl_state+3*g.ess.n_ess,kl-2) = ...
            g.ess.dess3(:,kl);
        d_vector(dcl_state+3*g.ess.n_ess+1:dcl_state+4*g.ess.n_ess,kl-2) = ...
            g.ess.dess4(:,kl);
        d_vector(dcl_state+4*g.ess.n_ess+1:dcl_state+5*g.ess.n_ess,kl-2) = 0.0;
        %   g.ess.dess5(:,kl);  % 48
    end

    if (g.lsc.n_lsc ~= 0)
        d_vector(ess_state+0*g.lsc.n_lsc+1:ess_state+1*g.lsc.n_lsc,kl-2) = ...
            g.lsc.dlsc1(:,kl);
        d_vector(ess_state+1*g.lsc.n_lsc+1:ess_state+2*g.lsc.n_lsc,kl-2) = ...
            g.lsc.dlsc2(:,kl);
        d_vector(ess_state+2*g.lsc.n_lsc+1:ess_state+3*g.lsc.n_lsc,kl-2) = ...
            g.lsc.dlsc3(:,kl);
        d_vector(ess_state+3*g.lsc.n_lsc+1:ess_state+4*g.lsc.n_lsc,kl-2) = ...
            g.lsc.dlsc4(:,kl);
        d_vector(ess_state+4*g.lsc.n_lsc+1:ess_state+5*g.lsc.n_lsc,kl-2) = ...
            g.lsc.dlsc5(:,kl);
        d_vector(ess_state+5*g.lsc.n_lsc+1:ess_state+6*g.lsc.n_lsc,kl-2) = ...
            g.lsc.dlsc6(:,kl);
        d_vector(ess_state+6*g.lsc.n_lsc+1:ess_state+7*g.lsc.n_lsc,kl-2) = ...
            g.lsc.dlsc7(:,kl);
        d_vector(ess_state+7*g.lsc.n_lsc+1:ess_state+8*g.lsc.n_lsc,kl-2) = ...
            g.lsc.dlsc8(:,kl);
        d_vector(ess_state+8*g.lsc.n_lsc+1:ess_state+9*g.lsc.n_lsc,kl-2) = ...
            g.lsc.dlsc9(:,kl);
        d_vector(ess_state+9*g.lsc.n_lsc+1:ess_state+10*g.lsc.n_lsc,kl-2) = ...
            g.lsc.dlsc10(:,kl);
        d_vector(ess_state+10*g.lsc.n_lsc+1:ess_state+11*g.lsc.n_lsc,kl-2) = ...
            g.lsc.dlsc11(:,kl);
        d_vector(ess_state+11*g.lsc.n_lsc+1:ess_state+12*g.lsc.n_lsc,kl-2) = ...
            g.lsc.dlsc12(:,kl);
        d_vector(ess_state+12*g.lsc.n_lsc+1:ess_state+13*g.lsc.n_lsc,kl-2) = ...
            g.lsc.dlsc13(:,kl);
        d_vector(ess_state+13*g.lsc.n_lsc+1:ess_state+14*g.lsc.n_lsc,kl-2) = ...
            g.lsc.dlsc14(:,kl);
        d_vector(ess_state+14*g.lsc.n_lsc+1:ess_state+15*g.lsc.n_lsc,kl-2) = ...
            g.lsc.dlsc15(:,kl);  % 63
    end

    if (g.reec.n_reec ~= 0)
        d_vector(lsc_state+0*g.reec.n_reec+1:lsc_state+1*g.reec.n_reec,kl-2) = ...
            g.reec.dreec1(:,kl);
        d_vector(lsc_state+1*g.reec.n_reec+1:lsc_state+2*g.reec.n_reec,kl-2) = ...
            g.reec.dreec2(:,kl);
        d_vector(lsc_state+2*g.reec.n_reec+1:lsc_state+3*g.reec.n_reec,kl-2) = ...
            g.reec.dreec3(:,kl);
        d_vector(lsc_state+3*g.reec.n_reec+1:lsc_state+4*g.reec.n_reec,kl-2) = ...
            g.reec.dreec4(:,kl);
        d_vector(lsc_state+4*g.reec.n_reec+1:lsc_state+5*g.reec.n_reec,kl-2) = ...
            g.reec.dreec5(:,kl);
        d_vector(lsc_state+5*g.reec.n_reec+1:lsc_state+6*g.reec.n_reec,kl-2) = ...
            g.reec.dreec6(:,kl);
        d_vector(lsc_state+6*g.reec.n_reec+1:lsc_state+7*g.reec.n_reec,kl-2) = ...
            g.reec.dreec7(:,kl);
        d_vector(lsc_state+7*g.reec.n_reec+1:lsc_state+8*g.reec.n_reec,kl-2) = ...
            g.reec.dreec8(:,kl);
        d_vector(lsc_state+8*g.reec.n_reec+1:lsc_state+9*g.reec.n_reec,kl-2) = ...
            g.reec.dreec9(:,kl);
        d_vector(lsc_state+9*g.reec.n_reec+1:lsc_state+10*g.reec.n_reec,kl-2) = ...
            g.reec.dreec10(:,kl);  % 73
    end

    if (g.gfma.n_gfma ~= 0)
        % note: overload mitigation states not included in linearization;
        %       this includes gfma2, gfma3, gfma6, gfma7
        d_vector(reec_state+0*g.gfma.n_gfma+1:reec_state+1*g.gfma.n_gfma,kl-2) = ...
            g.gfma.dgfma1(:,kl);
        d_vector(reec_state+1*g.gfma.n_gfma+1:reec_state+2*g.gfma.n_gfma,kl-2) = ...
            g.gfma.dgfma4(:,kl);
        d_vector(reec_state+2*g.gfma.n_gfma+1:reec_state+3*g.gfma.n_gfma,kl-2) = ...
            g.gfma.dgfma5(:,kl);
        d_vector(reec_state+3*g.gfma.n_gfma+1:reec_state+4*g.gfma.n_gfma,kl-2) = ...
            g.gfma.dgfma8(:,kl);
        d_vector(reec_state+4*g.gfma.n_gfma+1:reec_state+5*g.gfma.n_gfma,kl-2) = ...
            g.gfma.dgfma9(:,kl);
        d_vector(reec_state+5*g.gfma.n_gfma+1:reec_state+6*g.gfma.n_gfma,kl-2) = ...
            g.gfma.dgfma10(:,kl);  % 77
    end
end

%-----------------------------------------------------------------------------%
% verify that the inputs to the finite difference calculation are sensible

d_verify = [abs(d_vector), zeros(size(d_vector,1),1)];
for kl = 1:size(d_verify,1)
    % checking for pertubation results that move in the same direction
    if (d_verify(kl,1) > 1e-12 && d_verify(kl,2) > 1e-12)
        if (sign(d_vector(kl,1))*sign(d_vector(kl,2)) > 0)
            estr = '\np_file: nonzero finite difference components ';
            estr = [estr, 'have the same sign at state index %0.0f.'];
            warning(sprintf(estr,kl));
        end
    end

    % checking for mismatched perturbation results
    if (d_verify(kl,1) > 1e-15 && d_verify(kl,2) > 1e-15)
        d_row_max = max(d_verify(kl,:));
        d_verify(kl,:) = d_verify(kl,:)/d_row_max;

        if (d_verify(kl,1) > d_verify(kl,2))
            d_verify(kl,3) = 1 - d_verify(kl,2)/d_verify(kl,1);
        else
            d_verify(kl,3) = 1 - d_verify(kl,1)/d_verify(kl,2);
        end

        if (d_verify(kl,3) > 0.95)
            estr = '\np_file: mismatch between finite difference ';
            estr = [estr, 'components is too large at state index %0.0f.'];
            warning(sprintf(estr,kl));
        end
    end
end

d_vector = (d_vector(:,1) - d_vector(:,2))/2;  % center difference method
% d_vector = d_vector(:,1);                    % forward difference method

%-----------------------------------------------------------------------------%
% form state matrices

if (c_state == 0)
    if (k == 1)
        j_state = j;
    else
        j_state = j + sum(state(1:k-1));
    end

    if ~isempty(not_ib_ivm_idx)
        k_nib_idx = find(not_ib_ivm_idx == k);
    else
        k_nib_idx = k;
    end

    if (j == 2)
        if ~isempty(k_nib_idx)
            c_spd(k_nib_idx,j_state) = 1;
        end
    end

    a_mat(:,j_state) = p_mat*d_vector/pert;

    % form output matrices (ivm's have pelect but not telect or pmech)
    c_p(g.mac.not_ib_idx,j_state) = (g.mac.pelect(g.mac.not_ib_idx,3) ...
                                     - g.mac.pelect(g.mac.not_ib_idx,4)) ...
                                    .*g.mac.mac_pot(g.mac.not_ib_idx,1)/(2*pert);
    c_t(not_ib_ivm_idx,j_state) = (g.mac.telect(not_ib_ivm_idx,3) ...
                                   - g.mac.telect(not_ib_ivm_idx,4))/(2*pert);
    c_pm(not_ib_ivm_idx,j_state) = (g.mac.pmech(not_ib_ivm_idx,3) ...
                                    - g.mac.pmech(not_ib_ivm_idx,4))/(2*pert);
    c_v(:,j_state) = (abs(g.bus.bus_v(:,3)) ...
                      - abs(g.bus.bus_v(:,4)))/(2*pert);
    c_ang(:,j_state) = (g.bus.theta(:,3)-g.bus.theta(:,4))/(2*pert);
    c_curd(:,j_state) = (g.mac.curd(:,3)-g.mac.curd(:,4))/(2*pert);
    c_curq(:,j_state) = (g.mac.curq(:,3)-g.mac.curq(:,4))/(2*pert);

    if (g.exc.n_exc ~= 0)
        c_Efd(:,j_state) = (g.exc.Efd(:,3)-g.exc.Efd(:,4))/(2*pert);
    end

    if ~isempty(g.lmon.lmon_con)
        from_idx = g.bus.bus_int(line(g.lmon.lmon_con,1));
        to_idx = g.bus.bus_int(line(g.lmon.lmon_con,2));
        V1 = g.bus.bus_v(from_idx,4);
        V2 = g.bus.bus_v(to_idx,4);
        [s11,s21] = line_pq(V1,V2,R,X,B,tap,phi);
        [l_if1,l_it1] = line_cur(V1,V2,R,X,B,tap,phi);
        V1 = g.bus.bus_v(from_idx,3);
        V2 = g.bus.bus_v(to_idx,3);
        [s12,s22] = line_pq(V1,V2,R,X,B,tap,phi);
        [l_if2,l_it2] = line_cur(V1,V2,R,X,B,tap,phi);
        c_pf1(:,j_state) = real(s12-s11)/(2*pert);
        c_qf1(:,j_state) = imag(s12-s11)/(2*pert);
        c_pf2(:,j_state) = real(s22-s21)/(2*pert);
        c_qf2(:,j_state) = imag(s22-s21)/(2*pert);
        c_ilmf(:,j_state) = (abs(l_if2)-abs(l_if1))/(2*pert);
        c_ilmt(:,j_state) = (abs(l_it2)-abs(l_it1))/(2*pert);
        c_ilrf(:,j_state) = real(l_if2-l_if1)/(2*pert);
        c_ilif(:,j_state) = imag(l_if2-l_if1)/(2*pert);
        c_ilrt(:,j_state) = real(l_it2-l_it1)/(2*pert);
        c_ilit(:,j_state) = imag(l_it2-l_it1)/(2*pert);
    end

    if (g.dc.n_conv ~= 0)
        c_dcir(:,j_state) = (g.dc.i_dcr(:,3)-g.dc.i_dcr(:,4))/(2*pert);
        c_dcii(:,j_state) = (g.dc.i_dci(:,3)-g.dc.i_dci(:,4))/(2*pert);
        c_Vdcr(:,j_state) = (g.dc.Vdc(g.dc.r_idx,3)-g.dc.Vdc(g.dc.r_idx,4))/(2*pert);
        c_Vdci(:,j_state) = (g.dc.Vdc(g.dc.i_idx,3)-g.dc.Vdc(g.dc.i_idx,4))/(2*pert);
    end

    if (g.ess.n_ess ~= 0)
        c_ess_p(:,j_state) = real(g.ess.ess_scmd(:,3)-g.ess.ess_scmd(:,4))/(2*pert);
        c_ess_q(:,j_state) = imag(g.ess.ess_scmd(:,3)-g.ess.ess_scmd(:,4))/(2*pert);
    end

    if (g.reec.n_reec ~= 0)
        c_reec_v(:,j_state) = (g.reec.reec1(:,3) - g.reec.reec1(:,4))/(2*pert);
        c_reec_q(:,j_state) = (g.reec.reec3(:,3) - g.reec.reec3(:,4))/(2*pert);
    end

    if (g.gfma.n_gfma ~= 0)
        c_gfma_p(:,j_state) = (g.gfma.gfma8(:,3) - g.gfma.gfma8(:,4))/(2*pert);
        c_gfma_q(:,j_state) = (g.gfma.gfma9(:,3) - g.gfma.gfma9(:,4))/(2*pert);
    end
else
    % form b and d matrices
    if (c_state == 1)
        b_vr(:,vr_input) = p_mat*d_vector/pert;
        d_pvr(:,vr_input) = (g.mac.pelect(:,3)-g.mac.pelect(:,4)) ...
                            .*g.mac.mac_pot(:,1)/(2*pert);
        d_vvr(:,vr_input) = (abs(g.bus.bus_v(:,3)) ...
                             - abs(g.bus.bus_v(:,4)))/(2*pert);
        d_angvr(:,vr_input) = (g.bus.theta(:,3)-g.bus.theta(:,4))/(2*pert);
    elseif (c_state == 2)
        b_pr(:,pr_input) = p_mat*d_vector/pert;
        d_ppr(:,pr_input) = (g.mac.pelect(:,3)-g.mac.pelect(:,4)) ...
                            .*g.mac.mac_pot(:,1)/(2*pert);
        d_vpr(:,pr_input) = (abs(g.bus.bus_v(:,3)) ...
                             - abs(g.bus.bus_v(:,4)))/(2*pert);
        d_angpr(:,pr_input) = (g.bus.theta(:,3)-g.bus.theta(:,4))/(2*pert);
    elseif (c_state == 3)
        b_svc(:,svc_input) = p_mat*d_vector/pert;
        % note: d_svc is zero because of the time constant
    elseif (c_state == 4)
        b_tcsc(:,tcsc_input) = p_mat*d_vector/pert;
        % note: d_tcsc is zero because of the time constant
    elseif (c_state == 5)
        b_lmod(:,lmod_input) = p_mat*d_vector/pert;
        % note: d_lmod is zero because of the time constant
    elseif (c_state == 6)
        b_rlmod(:,rlmod_input) = p_mat*d_vector/pert;
        % note: d_lmod is zero because of the time constant
    elseif (c_state == 7)
        b_pwrmod_p(:,pwrmod_p_input) = p_mat*d_vector/pert;
        % note: d_pwrmod is zero because of the time constant
    elseif (c_state == 8)
        b_pwrmod_q(:,pwrmod_q_input) = p_mat*d_vector/pert;
        % note: d_pwrmod is zero because of the time constant
    elseif (c_state == 9)
        b_dcr(:,dcmod_input) = p_mat*d_vector/pert;
        d_pdcr(:,dcmod_input) = (g.mac.pelect(:,3)-g.mac.pelect(:,4)) ...
                                .*g.mac.mac_pot(:,1)/(2*pert);
        d_vdcr(:,dcmod_input) = (abs(g.bus.bus_v(:,3)) ...
                                 - abs(g.bus.bus_v(:,4)))/(2*pert);
        d_angdcr(:,dcmod_input) = (g.bus.theta(:,3)-g.bus.theta(:,4))/(2*pert);
        d_pdcr(:,dcmod_input) = (g.mac.pelect(:,3)-g.mac.pelect(:,4)) ...
                                .*g.mac.mac_pot(:,1)/(2*pert);
        d_idcdcr(:,dcmod_input) = (g.dc.i_dcr(:,3)-g.dc.i_dcr(:,4))/(2*pert);
        d_Vdcrdcr(:,dcmod_input) = (g.dc.Vdc(g.dc.r_idx,3) ...
                                    - g.dc.Vdc(g.dc.r_idx,4))/(2*pert);
        d_Vdcidcr(:,dcmod_input) = (g.dc.Vdc(g.dc.i_idx,3) ...
                                    - g.dc.Vdc(g.dc.i_idx,4))/(2*pert);

        if ~isempty(g.lmon.lmon_con)
            from_idx = g.bus.bus_int(line(g.lmon.lmon_con,1));
            to_idx = g.bus.bus_int(line(g.lmon.lmon_con,2));
            V1 = g.bus.bus_v(from_idx,4);
            V2 = g.bus.bus_v(to_idx,4);
            [s11,s21] = line_pq(V1,V2,R,X,B,tap,phi);
            [l_if1,l_it1] = line_cur(V1,V2,R,X,B,tap,phi);
            V1 = g.bus.bus_v(from_idx,3);
            V2 = g.bus.bus_v(to_idx,3);
            [s12,s22] = line_pq(V1,V2,R,X,B,tap,phi);
            [l_if2,l_it2] = line_cur(V1,V2,R,X,B,tap,phi);
            d_pf1cdr(:,dcmod_input) = real(s12-s11)/(2*pert);
            d_qf1dcr(:,dcmod_input) = imag(s12-s11)/(2*pert);
            d_pf2dcr(:,dcmod_input) = real(s22-s21)/(2*pert);
            d_qf2dcr(:,dcmod_input) = imag(s22-s21)/(2*pert);
            d_ilmfdcr(:,dcmod_input) = (abs(l_if2)-abs(l_if1))/(2*pert);
            d_ilmtdcr(:,dcmod_input) = (abs(l_it2)-abs(l_it1))/(2*pert);
            d_ilrfdcr(:,dcmod_input) = real(l_if2-l_if1)/(2*pert);
            d_ilifdcr(:,dcmod_input) = imag(l_if2-l_if1)/(2*pert);
            d_ilrtdcr(:,dcmod_input) = real(l_it2-l_it1)/(2*pert);
            d_ilitdcr(:,dcmod_input) = imag(l_it2-l_it1)/(2*pert);
        end
    elseif (c_state == 10)
        b_dci(:,dcmod_input) = p_mat*d_vector/pert;
        d_pdci(:,dcmod_input) = (g.mac.pelect(:,3)-g.mac.pelect(:,4)) ...
                                .*g.mac.mac_pot(:,1)/(2*pert);
        d_vdci(:,dcmod_input) = (abs(g.bus.bus_v(:,3)) ...
                                 - abs(g.bus.bus_v(:,4)))/(2*pert);
        d_angdci(:,dcmod_input) = (g.bus.theta(:,3)-g.bus.theta(:,4))/(2*pert);
        d_pdci(:,dcmod_input) = (g.mac.pelect(:,3)-g.mac.pelect(:,4)) ...
                                .*g.mac.mac_pot(:,1)/(2*pert);
        d_idcdci(:,dcmod_input) = (g.dc.i_dci(:,3)-g.dc.i_dci(:,4))/(2*pert);
        d_Vdcrdci(:,dcmod_input) = (g.dc.Vdc(g.dc.r_idx,3) ...
                                    - g.dc.Vdc(g.dc.r_idx,4))/(2*pert);
        d_Vdcidci(:,dcmod_input) = (g.dc.Vdc(g.dc.i_idx,3) ...
                                    - g.dc.Vdc(g.dc.i_idx,4))/(2*pert);

        if ~isempty(g.lmon.lmon_con)
            from_idx = g.bus.bus_int(line(g.lmon.lmon_con,1));
            to_idx = g.bus.bus_int(line(g.lmon.lmon_con,2));
            V1 = g.bus.bus_v(from_idx,4);
            V2 = g.bus.bus_v(to_idx,4);
            [s11,s21] = line_pq(V1,V2,R,X,B,tap,phi);
            [l_if1,l_it1] = line_cur(V1,V2,R,X,B,tap,phi);
            V1 = g.bus.bus_v(from_idx,3);
            V2 = g.bus.bus_v(to_idx,3);
            [s12,s22] = line_pq(V1,V2,R,X,B,tap,phi);
            [l_if2,l_it2] = line_cur(V1,V2,R,X,B,tap,phi);
            d_pf1cdi(:,dcmod_input) = real(s12-s11)/(2*pert);
            d_qf1dci(:,dcmod_input) = imag(s12-s11)/(2*pert);
            d_pf2dci(:,dcmod_input) = real(s22-s21)/(2*pert);
            d_qf2dci(:,dcmod_input) = imag(s22-s21)/(2*pert);
            d_ilmfdci(:,dcmod_input) = (abs(l_if2)-abs(l_if1))/(2*pert);
            d_ilmtdci(:,dcmod_input) = (abs(l_it2)-abs(l_it1))/(2*pert);
            d_ilrfdci(:,dcmod_input) = real(l_if2-l_if1)/(2*pert);
            d_ilifdci(:,dcmod_input) = imag(l_if2-l_if1)/(2*pert);
            d_ilrtdci(:,dcmod_input) = real(l_it2-l_it1)/(2*pert);
            d_ilitdci(:,dcmod_input) = imag(l_it2-l_it1)/(2*pert);
        end
    elseif (c_state == 11)
        b_ess_p(:,ess_input) = p_mat*d_vector/pert;
        % note: d_ess is zero because of the time constant
    elseif (c_state == 12)
        b_ess_q(:,ess_input) = p_mat*d_vector/pert;
        % note: d_ess is zero because of the time constant
    elseif (c_state == 13)
        b_reec_v(:,reec_input) = p_mat*d_vector/pert;
        % note: d_reec is zero because of the time constant
    elseif (c_state == 14)
        b_reec_q(:,reec_input) = p_mat*d_vector/pert;
        % note: d_reec is zero because of the time constant
    elseif (c_state == 15)
        b_gfma_p(:,gfma_input) = p_mat*d_vector/pert;
        % note: d_gfma is zero because of the time constant
    elseif (c_state == 16)
        b_gfma_q(:,gfma_input) = p_mat*d_vector/pert;
        % note: d_gfma is zero because of the time constant
    end
end

%-----------------------------------------------------------------------------%
% reset states to initial values

for kl = 2:4
    if (g.bus.n_bus ~= 0)
        g.bus.bus_v(:,kl) = g.bus.bus_v(:,1);
        g.bus.theta(:,kl) = g.bus.theta(:,1);
    end

    if (g.mac.n_mac ~= 0)
        g.mac.pmech(:,kl) = g.mac.pmech(:,1);
        g.mac.telect(:,kl) = g.mac.telect(:,1);
        g.mac.pelect(:,kl) = g.mac.pelect(:,1);
        g.mac.qelect(:,kl) = g.mac.qelect(:,1);
        g.mac.ed(:,kl) = g.mac.ed(:,1);
        g.mac.eq(:,kl) = g.mac.eq(:,1);
        g.mac.curd(:,kl) = g.mac.curd(:,1);
        g.mac.curq(:,kl) = g.mac.curq(:,1);
        g.mac.curdg(:,kl) = g.mac.curdg(:,1);
        g.mac.curqg(:,kl) = g.mac.curqg(:,1);
        g.mac.fldcur(:,kl) = g.mac.fldcur(:,1);
        g.mac.eterm(:,kl) = g.mac.eterm(:,1);
        g.mac.vex(:,kl) = g.mac.vex(:,1);
        g.mac.cur_re(:,kl) = g.mac.cur_re(:,1);
        g.mac.cur_im(:,kl) = g.mac.cur_im(:,1);
        g.mac.psi_re(:,kl) = g.mac.psi_re(:,1);
        g.mac.psi_im(:,kl) = g.mac.psi_im(:,1);
        g.mac.pm_sig(:,kl) = g.mac.pm_sig(:,1);

        g.mac.mac_ang(:,kl) = g.mac.mac_ang(:,1);
        g.mac.mac_spd(:,kl) = g.mac.mac_spd(:,1);
        g.mac.edprime(:,kl) = g.mac.edprime(:,1);
        g.mac.eqprime(:,kl) = g.mac.eqprime(:,1);
        g.mac.psikd(:,kl) = g.mac.psikd(:,1);
        g.mac.psikq(:,kl) = g.mac.psikq(:,1);

        g.mac.dmac_ang(:,kl) = g.mac.dmac_ang(:,1);
        g.mac.dmac_spd(:,kl) = g.mac.dmac_spd(:,1);
        g.mac.dedprime(:,kl) = g.mac.dedprime(:,1);
        g.mac.deqprime(:,kl) = g.mac.deqprime(:,1);
        g.mac.dpsikd(:,kl) = g.mac.dpsikd(:,1);
        g.mac.dpsikq(:,kl) = g.mac.dpsikq(:,1);
    end

    if (g.exc.n_exc ~= 0)
        g.exc.V_A(:,kl) = g.exc.V_A(:,1);  % not a state variable

        g.exc.V_TR(:,kl) = g.exc.V_TR(:,1);
        g.exc.V_As(:,kl) = g.exc.V_As(:,1);
        g.exc.V_R(:,kl) = g.exc.V_R(:,1);
        g.exc.Efd(:,kl) = g.exc.Efd(:,1);
        g.exc.R_f(:,kl) = g.exc.R_f(:,1);

        g.exc.dV_TR(:,kl) = g.exc.dV_TR(:,1);
        g.exc.dV_As(:,kl) = g.exc.dV_As(:,1);
        g.exc.dV_R(:,kl) = g.exc.dV_R(:,1);
        g.exc.dEfd(:,kl) = g.exc.dEfd(:,1);
        g.exc.dR_f(:,kl) = g.exc.dR_f(:,1);
    end

    if (g.pss.n_pss ~= 0)
        g.pss.pss1(:,kl) = g.pss.pss1(:,1);
        g.pss.pss2(:,kl) = g.pss.pss2(:,1);
        g.pss.pss3(:,kl) = g.pss.pss3(:,1);

        g.pss.dpss1(:,kl) = g.pss.dpss1(:,1);
        g.pss.dpss2(:,kl) = g.pss.dpss2(:,1);
        g.pss.dpss3(:,kl) = g.pss.dpss3(:,1);
    end

    if (g.dpw.n_dpw ~= 0)
        g.dpw.sdpw1(:,kl) = g.dpw.sdpw1(:,1);
        g.dpw.sdpw2(:,kl) = g.dpw.sdpw2(:,1);
        g.dpw.sdpw3(:,kl) = g.dpw.sdpw3(:,1);
        g.dpw.sdpw4(:,kl) = g.dpw.sdpw4(:,1);
        g.dpw.sdpw5(:,kl) = g.dpw.sdpw5(:,1);
        g.dpw.sdpw6(:,kl) = g.dpw.sdpw6(:,1);

        g.dpw.dsdpw1(:,kl) = g.dpw.dsdpw1(:,1);
        g.dpw.dsdpw2(:,kl) = g.dpw.dsdpw2(:,1);
        g.dpw.dsdpw3(:,kl) = g.dpw.dsdpw3(:,1);
        g.dpw.dsdpw4(:,kl) = g.dpw.dsdpw4(:,1);
        g.dpw.dsdpw5(:,kl) = g.dpw.dsdpw5(:,1);
        g.dpw.dsdpw6(:,kl) = g.dpw.dsdpw6(:,1);
    end

    if (g.tg.n_tg_tot ~= 0)
        g.tg.tg1(:,kl) = g.tg.tg1(:,1);
        g.tg.tg2(:,kl) = g.tg.tg2(:,1);
        g.tg.tg3(:,kl) = g.tg.tg3(:,1);
        g.tg.tg4(:,kl) = g.tg.tg4(:,1);
        g.tg.tg5(:,kl) = g.tg.tg5(:,1);

        g.tg.dtg1(:,kl) = g.tg.dtg1(:,1);
        g.tg.dtg2(:,kl) = g.tg.dtg2(:,1);
        g.tg.dtg3(:,kl) = g.tg.dtg3(:,1);
        g.tg.dtg4(:,kl) = g.tg.dtg4(:,1);
        g.tg.dtg5(:,kl) = g.tg.dtg5(:,1);
    end

    if (g.ind.n_mot ~= 0)
        g.ind.vdp(:,kl) = g.ind.vdp(:,1);
        g.ind.vqp(:,kl) = g.ind.vqp(:,1);
        g.ind.slip(:,kl) = g.ind.slip(:,1);

        g.ind.dvdp(:,kl) = g.ind.dvdp(:,1);
        g.ind.dvqp(:,kl) = g.ind.dvqp(:,1);
        g.ind.dslip(:,kl) = g.ind.dslip(:,1);
    end

    if (g.igen.n_ig ~= 0)
        g.igen.vdpig(:,kl) = g.igen.vdpig(:,1);
        g.igen.vqpig(:,kl) = g.igen.vqpig(:,1);
        g.igen.slig(:,kl) = g.igen.slig(:,1);

        g.igen.dvdpig(:,kl) = g.igen.dvdpig(:,1);
        g.igen.dvqpig(:,kl) = g.igen.dvqpig(:,1);
        g.igen.dslig(:,kl) = g.igen.dslig(:,1);
    end

    if (g.svc.n_svc ~= 0)
        g.svc.B_cv(:,kl) = g.svc.B_cv(:,1);
        g.svc.B_con(:,kl) = g.svc.B_con(:,1);
        g.svc.dB_cv(:,kl) = g.svc.dB_cv(:,1);
        g.svc.dB_con(:,kl) = g.svc.dB_con(:,1);
        g.svc.svc_sig(:,kl) = g.svc.svc_sig(:,1);
    end

    if (g.tcsc.n_tcsc ~= 0)
        g.tcsc.B_tcsc(:,kl) = g.tcsc.B_tcsc(:,1);
        g.tcsc.dB_tcsc(:,kl) = g.tcsc.dB_tcsc(:,1);
        g.tcsc.tcsc_sig(:,kl) = g.tcsc.tcsc_sig(:,1);
    end

    if (g.lmod.n_lmod ~= 0)
        g.lmod.lmod_st(:,kl) = g.lmod.lmod_st(:,1);
        g.lmod.dlmod_st(:,kl) = g.lmod.dlmod_st(:,1);
        g.lmod.lmod_sig(:,kl) = g.lmod.lmod_sig(:,1);
    end

    if (g.rlmod.n_rlmod ~= 0)
        g.rlmod.rlmod_st(:,kl) = g.rlmod.rlmod_st(:,1);
        g.rlmod.drlmod_st(:,kl) = g.rlmod.drlmod_st(:,1);
        g.rlmod.rlmod_sig(:,kl) = g.rlmod.rlmod_sig(:,1);
    end

    if (g.pwr.n_pwrmod ~= 0)
        g.pwr.pwrmod_p_st(:,kl) = g.pwr.pwrmod_p_st(:,1);
        g.pwr.dpwrmod_p_st(:,kl) = g.pwr.dpwrmod_p_st(:,1);
        g.pwr.pwrmod_p_sig(:,kl) = g.pwr.pwrmod_p_sig(:,1);

        g.pwr.pwrmod_q_st(:,kl) = g.pwr.pwrmod_q_st(:,1);
        g.pwr.dpwrmod_q_st(:,kl) = g.pwr.dpwrmod_q_st(:,1);
        g.pwr.pwrmod_q_sig(:,kl) = g.pwr.pwrmod_q_sig(:,1);
    end

    if (g.dc.n_conv ~= 0)
        g.dc.v_conr(:,kl) = g.dc.v_conr(:,1);
        g.dc.v_coni(:,kl) = g.dc.v_coni(:,1);
        g.dc.i_dcr(:,kl) = g.dc.i_dcr(:,1);
        g.dc.i_dci(:,kl) = g.dc.i_dci(:,1);
        g.dc.v_dcc(:,kl) = g.dc.v_dcc(:,1);

        g.dc.dv_conr(:,kl) = g.dc.dv_conr(:,1);
        g.dc.dv_coni(:,kl) = g.dc.dv_coni(:,1);
        g.dc.di_dcr(:,kl) = g.dc.di_dcr(:,1);
        g.dc.di_dci(:,kl) = g.dc.di_dci(:,1);
        g.dc.dv_dcc(:,kl) = g.dc.dv_dcc(:,1);

        g.dc.Vdc(:,kl) = g.dc.Vdc(:,1);
        g.dc.i_dc(:,kl) = g.dc.i_dc(:,1);
        g.dc.alpha(:,kl) = g.dc.alpha(:,1);
        g.dc.gamma(:,kl) = g.dc.gamma(:,1);
        g.dc.dc_sig(:,kl) = g.dc.dc_sig(:,1);
    end

    if (g.ess.n_ess ~= 0)
        g.ess.ess1(:,kl) = g.ess.ess1(:,1);
        g.ess.ess2(:,kl) = g.ess.ess2(:,1);
        g.ess.ess3(:,kl) = g.ess.ess3(:,1);
        g.ess.ess4(:,kl) = g.ess.ess4(:,1);
        g.ess.ess5(:,kl) = g.ess.ess5(:,1);

        g.ess.dess1(:,kl) = g.ess.dess1(:,1);
        g.ess.dess2(:,kl) = g.ess.dess2(:,1);
        g.ess.dess3(:,kl) = g.ess.dess3(:,1);
        g.ess.dess4(:,kl) = g.ess.dess4(:,1);
        g.ess.dess5(:,kl) = g.ess.dess5(:,1);

        g.ess.ess_sig(:,kl) = g.ess.ess_sig(:,1);
        g.ess.ess_cur(:,kl) = g.ess.ess_cur(:,1);
        g.ess.ess_soc(:,kl) = g.ess.ess_soc(:,1);
        g.ess.ess_vmag(:,kl) = g.ess.ess_vmag(:,1);
        g.ess.ess_scmd(:,kl) = g.ess.ess_scmd(:,1);
        g.ess.ess_sinj(:,kl) = g.ess.ess_sinj(:,1);
        g.ess.ess_iord(:,kl) = g.ess.ess_iord(:,1);
        g.ess.ess_vmag_pade(:,kl) = g.ess.ess_vmag_pade(:,1);
    end

    if (g.lsc.n_lsc ~= 0)
        g.lsc.lsc1(:,kl) = g.lsc.lsc1(:,1);
        g.lsc.lsc2(:,kl) = g.lsc.lsc2(:,1);
        g.lsc.lsc3(:,kl) = g.lsc.lsc3(:,1);
        g.lsc.lsc4(:,kl) = g.lsc.lsc4(:,1);
        g.lsc.lsc5(:,kl) = g.lsc.lsc5(:,1);
        g.lsc.lsc6(:,kl) = g.lsc.lsc6(:,1);
        g.lsc.lsc7(:,kl) = g.lsc.lsc7(:,1);
        g.lsc.lsc8(:,kl) = g.lsc.lsc8(:,1);
        g.lsc.lsc9(:,kl) = g.lsc.lsc9(:,1);
        g.lsc.lsc10(:,kl) = g.lsc.lsc10(:,1);
        g.lsc.lsc11(:,kl) = g.lsc.lsc11(:,1);
        g.lsc.lsc12(:,kl) = g.lsc.lsc12(:,1);
        g.lsc.lsc13(:,kl) = g.lsc.lsc13(:,1);
        g.lsc.lsc14(:,kl) = g.lsc.lsc14(:,1);
        g.lsc.lsc15(:,kl) = g.lsc.lsc15(:,1);

        g.lsc.dlsc1(:,kl) = g.lsc.dlsc1(:,1);
        g.lsc.dlsc2(:,kl) = g.lsc.dlsc2(:,1);
        g.lsc.dlsc3(:,kl) = g.lsc.dlsc3(:,1);
        g.lsc.dlsc4(:,kl) = g.lsc.dlsc4(:,1);
        g.lsc.dlsc5(:,kl) = g.lsc.dlsc5(:,1);
        g.lsc.dlsc6(:,kl) = g.lsc.dlsc6(:,1);
        g.lsc.dlsc7(:,kl) = g.lsc.dlsc7(:,1);
        g.lsc.dlsc8(:,kl) = g.lsc.dlsc8(:,1);
        g.lsc.dlsc9(:,kl) = g.lsc.dlsc9(:,1);
        g.lsc.dlsc10(:,kl) = g.lsc.dlsc10(:,1);
        g.lsc.dlsc11(:,kl) = g.lsc.dlsc11(:,1);
        g.lsc.dlsc12(:,kl) = g.lsc.dlsc12(:,1);
        g.lsc.dlsc13(:,kl) = g.lsc.dlsc13(:,1);
        g.lsc.dlsc14(:,kl) = g.lsc.dlsc14(:,1);
        g.lsc.dlsc15(:,kl) = g.lsc.dlsc15(:,1);

        g.lsc.lsc_scmd(:,kl) = g.lsc.lsc_scmd(:,1);
        g.lsc.theta_sensor(:,kl) = g.lsc.theta_sensor(:,1);
        g.lsc.theta_coi_pade(:,kl) = g.lsc.theta_coi_pade(:,1);
        g.lsc.theta_lsc_pade(:,kl) = g.lsc.theta_lsc_pade(:,1);
      % g.lsc.theta_coi_del(:,kl) = g.lsc.theta_coi_del(:,1);
    end

    if (g.reec.n_reec ~= 0)
        g.reec.reec1(:,kl) = g.reec.reec1(:,1);
        g.reec.reec2(:,kl) = g.reec.reec2(:,1);
        g.reec.reec3(:,kl) = g.reec.reec3(:,1);
        g.reec.reec4(:,kl) = g.reec.reec4(:,1);
        g.reec.reec5(:,kl) = g.reec.reec5(:,1);
        g.reec.reec6(:,kl) = g.reec.reec6(:,1);
        g.reec.reec7(:,kl) = g.reec.reec7(:,1);
        g.reec.reec8(:,kl) = g.reec.reec8(:,1);
        g.reec.reec9(:,kl) = g.reec.reec9(:,1);
        g.reec.reec10(:,kl) = g.reec.reec10(:,1);

        g.reec.dreec1(:,kl) = g.reec.dreec1(:,1);
        g.reec.dreec2(:,kl) = g.reec.dreec2(:,1);
        g.reec.dreec3(:,kl) = g.reec.dreec3(:,1);
        g.reec.dreec4(:,kl) = g.reec.dreec4(:,1);
        g.reec.dreec5(:,kl) = g.reec.dreec5(:,1);
        g.reec.dreec6(:,kl) = g.reec.dreec6(:,1);
        g.reec.dreec7(:,kl) = g.reec.dreec7(:,1);
        g.reec.dreec8(:,kl) = g.reec.dreec8(:,1);
        g.reec.dreec9(:,kl) = g.reec.dreec9(:,1);
        g.reec.dreec10(:,kl) = g.reec.dreec10(:,1);

        g.reec.icmd(:,kl) = g.reec.icmd(:,1);
        g.reec.paux(:,kl) = g.reec.paux(:,1);
        g.reec.reec_sig(:,kl) = g.reec.reec_sig(:,1);

        g.reec.iqmin = -100*ones(g.reec.n_reec,1);
        g.reec.iqmax = 100*ones(g.reec.n_reec,1);

        g.reec.pref(:,kl) = g.reec.pref(:,1);
        g.reec.qref(:,kl) = g.reec.qref(:,1);

      % static references don't get reset here

        g.reec.vdip = false(g.reec.n_reec,1);
        g.reec.vdip_tick = -1*ones(g.reec.n_reec,1);
        g.reec.vdip_time = zeros(g.reec.n_reec,1);
        g.reec.vdip_icmd = zeros(g.reec.n_reec,1);

        g.reec.vblk = false(g.reec.n_reec,1);
        g.reec.vblk_tick = -1*ones(g.reec.n_reec,1);
        g.reec.vblk_time = zeros(g.reec.n_reec,1);
    end

    if (g.gfma.n_gfma ~= 0)
        g.gfma.gfma1(:,kl) = g.gfma.gfma1(:,1);
        g.gfma.gfma2(:,kl) = g.gfma.gfma2(:,1);
        g.gfma.gfma3(:,kl) = g.gfma.gfma3(:,1);
        g.gfma.gfma4(:,kl) = g.gfma.gfma4(:,1);
        g.gfma.gfma5(:,kl) = g.gfma.gfma5(:,1);
        g.gfma.gfma6(:,kl) = g.gfma.gfma6(:,1);
        g.gfma.gfma7(:,kl) = g.gfma.gfma7(:,1);
        g.gfma.gfma8(:,kl) = g.gfma.gfma8(:,1);
        g.gfma.gfma9(:,kl) = g.gfma.gfma9(:,1);
        g.gfma.gfma10(:,kl) = g.gfma.gfma10(:,1);

        g.gfma.dgfma1(:,kl) = g.gfma.dgfma1(:,1);
        g.gfma.dgfma2(:,kl) = g.gfma.dgfma2(:,1);
        g.gfma.dgfma3(:,kl) = g.gfma.dgfma3(:,1);
        g.gfma.dgfma4(:,kl) = g.gfma.dgfma4(:,1);
        g.gfma.dgfma5(:,kl) = g.gfma.dgfma5(:,1);
        g.gfma.dgfma6(:,kl) = g.gfma.dgfma6(:,1);
        g.gfma.dgfma7(:,kl) = g.gfma.dgfma7(:,1);
        g.gfma.dgfma8(:,kl) = g.gfma.dgfma8(:,1);
        g.gfma.dgfma9(:,kl) = g.gfma.dgfma9(:,1);
        g.gfma.dgfma10(:,kl) = g.gfma.dgfma10(:,1);

        g.gfma.gfma_sig(:,kl) = g.gfma.gfma_sig(:,1);

        g.gfma.pset(:,kl) = g.gfma.pset(:,1);
        g.gfma.qset(:,kl) = g.gfma.qset(:,1);
        g.gfma.vset(:,kl) = g.gfma.vset(:,1);
        g.gfma.lim_flag = false(g.gfma.n_gfma,1);
    end
end

% eof
