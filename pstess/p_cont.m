% perturbation control file
%
% Purpose: (step 3) perform perturbation on each state in turn using
%           finite difference methods

%-----------------------------------------------------------------------------%
% Version history
%
% Purpose: modified for center-difference method
% Date:    2019
% Author:  Ryan Elliott
%
% Purpose: modified to include pwrmod
% Date:    2015
% Author:  D. Trudnowski
%
% Purpose: modified to include induction generators
%          modified to include load modulation
%          input disturbance modulation added to svc and lmod
% Author:  Graham Rogers
%
% Version: 1.0 (initial version)
% Date:    December 1996
% Author:  Graham Rogers
%-----------------------------------------------------------------------------%

vr_input = 0;
pr_input = 0;
c_state = 0;
p_ratio = 1e-4;

for k = 1:g.mac.n_mac
    not_ib = 1;
    if ~isempty(g.mac.mac_ib_idx)
        ib_chk = find(g.mac.mac_ib_idx == k);
        if ~isempty(ib_chk)
            not_ib = 0;
        end
    end

    if (not_ib == 1)
        disp('disturbing generator')

        j = 1;
        pert = p_ratio*abs(g.mac.mac_ang(k,1));
        pert = max(pert,p_ratio);
        g.mac.mac_ang(k,3) = g.mac.mac_ang(k,1) + pert;
        g.mac.mac_ang(k,4) = g.mac.mac_ang(k,1) - pert;

        run('p_file');
        st_name(k,j) = 1;

        k_ivm = 0;
        if ~isempty(g.mac.mac_ivm_idx)
            k_idx = find(g.mac.mac_ivm_idx == k);
            if ~isempty(k_idx)
                k_ivm = 1;
            end
        end

        if (k_ivm == 0)
            j = j + 1;
            pert = p_ratio*abs(g.mac.mac_spd(k,1));
            pert = max(pert,p_ratio);
            g.mac.mac_spd(k,3) = g.mac.mac_spd(k,1) + pert;
            g.mac.mac_spd(k,4) = g.mac.mac_spd(k,1) - pert;

            run('p_file');
            st_name(k,j) = 2;
        end

        k_tra = 0;
        k_sub = 0;
        if ~isempty(g.mac.mac_tra_idx)
            k_idx = find(g.mac.mac_tra_idx == k);
            if ~isempty(k_idx)
                k_tra = 1;
            end
        end

        if ~isempty(g.mac.mac_sub_idx)
            k_idx = find(g.mac.mac_sub_idx == k);
            if ~isempty(k_idx)
                k_sub = 1;
            end
        end

        if (k_tra == 1 || k_sub == 1 || k_ivm == 1)
            j = j + 1;
            pert = p_ratio*abs(g.mac.eqprime(k,1));
            pert = max(pert,p_ratio);
            g.mac.eqprime(k,3) = g.mac.eqprime(k,1) + pert;
            g.mac.eqprime(k,4) = g.mac.eqprime(k,1) - pert;

            run('p_file');
            st_name(k,j) = 3;
        end

        if (k_sub == 1)
            j = j + 1;
            pert = p_ratio*abs(g.mac.psikd(k,1));
            pert = max(pert,p_ratio);
            g.mac.psikd(k,3) = g.mac.psikd(k,1) + pert;
            g.mac.psikd(k,4) = g.mac.psikd(k,1) - pert;

            run('p_file');
            st_name(k,j) = 4;
        end

        if (k_tra == 1 || k_sub == 1)
            j = j + 1;
            pert = p_ratio*abs(g.mac.edprime(k,1));
            pert = max(pert,p_ratio);
            g.mac.edprime(k,3) = g.mac.edprime(k,1) + pert;
            g.mac.edprime(k,4) = g.mac.edprime(k,1) - pert;

            run('p_file');
            st_name(k,j) = 5;
        end

        if (k_sub == 1)
            j = j + 1;
            pert = p_ratio*abs(g.mac.psikq(k,1));
            pert = max(pert,p_ratio);
            g.mac.psikq(k,3) = g.mac.psikq(k,1) + pert;
            g.mac.psikq(k,4) = g.mac.psikq(k,1) - pert;

            run('p_file');
            st_name(k,j) = 6;
        end

        % exciters
        if ~isempty(g.exc.exc_con)
            run('p_exc');
        end

        % pss
        if ~isempty(g.pss.pss_con)
            run('p_pss');
        end

        if ~isempty(g.dpw.dpw_con)
            run('p_dpw');
        end

        % turbine/governor
        if ~isempty(g.tg.tg_con)
            run('p_tg');
        end

        % disturb the auxiliary inputs
        if (g.exc.n_exc ~= 0)
            c_state = 1;

            exc_number = find(g.mac.mac_int(g.exc.exc_con(:,2)) == k);
            if ~isempty(exc_number)
                disp('disturbing V_ref')
                vr_input = vr_input + 1;

                pert = p_ratio*abs(g.exc.exc_pot(exc_number,3));
                pert = max(pert,p_ratio);

                g.exc.exc_sig(exc_number,3) = g.exc.exc_sig(exc_number,1) + pert;
                g.exc.exc_sig(exc_number,4) = g.exc.exc_sig(exc_number,1) - pert;

                % resetting the input
                run('p_file');
                g.exc.exc_sig(exc_number,3) = g.exc.exc_sig(exc_number,1);
                g.exc.exc_sig(exc_number,4) = g.exc.exc_sig(exc_number,1);
            end
        end

        if (g.tg.n_tg_tot ~= 0)
            c_state = 2;

            tg_number = find(mac_tg == k);
            if isempty(tg_number)
                tg_number = find(mac_tgh == k);
            end

            if ~isempty(tg_number)
                disp('disturbing P_ref')
                pr_input = pr_input + 1;

                pert = p_ratio*abs(g.tg.tg_pot(tg_number,5));
                pert = max(pert,p_ratio);

                g.tg.tg_sig(tg_number,3) = g.tg.tg_sig(tg_number,1) + pert;
                g.tg.tg_sig(tg_number,4) = g.tg.tg_sig(tg_number,1) - pert;

                % resetting the input
                run('p_file');
                g.tg.tg_sig(tg_number,3) = g.tg.tg_sig(tg_number,1);
                g.tg.tg_sig(tg_number,4) = g.tg.tg_sig(tg_number,1);
            end
        end

        c_state = 0;
    end
end

% disturb induction motor states
if (g.ind.n_mot ~= 0)
    disp('disturbing induction motors')
    for k = g.mac.n_mac+1:n_gm
        k_ind = k - g.mac.n_mac;

        j = 1;
        pert = p_ratio*abs(g.ind.vdp(k_ind,1));
        pert = max(pert,p_ratio);
        g.ind.vdp(k_ind,3) = g.ind.vdp(k_ind,1) + pert;
        g.ind.vdp(k_ind,4) = g.ind.vdp(k_ind,1) - pert;

        run('p_file');
        st_name(k,j) = 26;

        j = j + 1;
        pert = p_ratio*abs(g.ind.vqp(k_ind,1));
        pert = max(pert,p_ratio);
        g.ind.vqp(k_ind,3) = g.ind.vqp(k_ind,1) + pert;
        g.ind.vqp(k_ind,4) = g.ind.vqp(k_ind,1) - pert;

        run('p_file');
        st_name(k,j) = 27;

        j = j + 1;
        pert = p_ratio*abs(g.ind.slip(k_ind,1));
        pert = max(pert,p_ratio);
        g.ind.slip(k_ind,3) = g.ind.slip(k_ind,1) + pert;
        g.ind.slip(k_ind,4) = g.ind.slip(k_ind,1) - pert;

        run('p_file');
        st_name(k,j) = 28;
    end
end

% disturb induction generator states
if (g.igen.n_ig ~= 0)
    disp('disturbing induction generators')
    for k = n_gm+1:n_tot
        k_ig = k - n_gm;

        j = 1;
        pert = p_ratio*abs(g.igen.vdpig(k_ig,1));
        pert = max(pert,p_ratio);
        g.igen.vdpig(k_ig,3) = g.igen.vdpig(k_ig,1) + pert;
        g.igen.vdpig(k_ig,4) = g.igen.vdpig(k_ig,1) - pert;

        run('p_file');
        st_name(k,j) = 29;

        j = j + 1;
        pert = p_ratio*abs(g.igen.vqpig(k_ig,1));
        pert = max(pert,p_ratio);
        g.igen.vqpig(k_ig,3) = g.igen.vqpig(k_ig,1) + pert;
        g.igen.vqpig(k_ig,4) = g.igen.vqpig(k_ig,1) - pert;

        run('p_file');
        st_name(k,j) = 30;

        j = j + 1;
        pert = p_ratio*abs(g.igen.slig(k_ig,1));
        pert = max(pert,p_ratio);
        g.igen.slig(k_ig,3) = g.igen.slig(k_ig,1) + pert;
        g.igen.slig(k_ig,4) = g.igen.slig(k_ig,1) - pert;

        run('p_file');
        st_name(k,j) = 31;
    end
end

% disturb svc states
n_tot_svc = n_tot + g.svc.n_svc;
if (g.svc.n_svc ~= 0)
    disp('disturbing svc')
    for k = n_tot+1:n_tot_svc
        k_svc = k - n_tot;

        j = 1;
        pert = p_ratio*abs(g.svc.B_cv(k_svc,1));
        pert = max(pert,p_ratio);
        g.svc.B_cv(k_svc,3) = g.svc.B_cv(k_svc,1) + pert;
        g.svc.B_cv(k_svc,4) = g.svc.B_cv(k_svc,1) - pert;

        run('p_file');
        st_name(k,j) = 32;

        % disturb B_con
        if ~isempty(g.svc.svcll_idx)
            j = j + 1;
            kcon = find(g.svc.svcll_idx == k_svc);
            if ~isempty(kcon)
                pert = p_ratio*abs(g.svc.B_con(k_svc,1));
                pert = max(pert,p_ratio);
                g.svc.B_con(k_svc,3) = g.svc.B_con(k_svc,1) + pert;
                g.svc.B_con(k_svc,4) = g.svc.B_con(k_svc,1) - pert;

                run('p_file');
                st_name(k,j) = 33;
            end
        end

        % disturb the auxiliary input
        disp('disturbing svc_sig')
        c_state = 3;
        svc_input = k_svc;

        pert = p_ratio;
        g.svc.svc_sig(k_svc,3) = g.svc.svc_sig(k_svc,1) + pert;
        g.svc.svc_sig(k_svc,4) = g.svc.svc_sig(k_svc,1) - pert;

        % resetting the input
        run('p_file');
        g.svc.svc_sig(k_svc,3) = g.svc.svc_sig(k_svc,1);
        g.svc.svc_sig(k_svc,4) = g.svc.svc_sig(k_svc,1);
        c_state = 0;
    end
end

% disturb tcsc states
n_tot_tcsc = n_tot_svc + g.tcsc.n_tcsc;
if (g.tcsc.n_tcsc ~= 0)
    disp('disturbing tcsc')
    for k = n_tot_svc+1:n_tot_tcsc
        k_tcsc = k - n_tot_svc;

        j = 1;
        pert = p_ratio*abs(g.tcsc.B_tcsc(k_tcsc,1));
        pert = max(pert,p_ratio);
        g.tcsc.B_tcsc(k_tcsc,3) = g.tcsc.B_tcsc(k_tcsc,1) + pert;
        g.tcsc.B_tcsc(k_tcsc,4) = g.tcsc.B_tcsc(k_tcsc,1) - pert;

        run('p_file');
        st_name(k,j) = 34;

        % disturb the auxiliary input
        disp('disturbing tcsc_sig')
        c_state = 4;
        tcsc_input = k_tcsc;

        pert = p_ratio;
        g.tcsc.tcsc_sig(k_tcsc,3) = g.tcsc.tcsc_sig(k_tcsc,1) + pert;
        g.tcsc.tcsc_sig(k_tcsc,4) = g.tcsc.tcsc_sig(k_tcsc,1) - pert;

        % resetting the input
        run('p_file');
        g.tcsc.tcsc_sig(k_tcsc,3) = g.tcsc.tcsc_sig(k_tcsc,1);
        g.tcsc.tcsc_sig(k_tcsc,4) = g.tcsc.tcsc_sig(k_tcsc,1);
        c_state = 0;
    end
end

% disturb lmod states
n_tot_lmod = n_tot_tcsc + g.lmod.n_lmod;
if (g.lmod.n_lmod ~= 0)
    disp('disturbing load modulation')
    for k = n_tot_tcsc+1:n_tot_lmod
        k_lmod = k - n_tot_tcsc;

        j = 1;
        pert = p_ratio*abs(g.lmod.lmod_st(k_lmod,1));
        pert = max(pert,p_ratio);
        g.lmod.lmod_st(k_lmod,3) = g.lmod.lmod_st(k_lmod,1) + pert;
        g.lmod.lmod_st(k_lmod,4) = g.lmod.lmod_st(k_lmod,1) - pert;

        run('p_file');
        st_name(k,j) = 35;

        % disturb the auxiliary input
        disp('disturbing lmod_sig')
        c_state = 5;
        lmod_input = k_lmod;

        pert = p_ratio;
        g.lmod.lmod_sig(k_lmod,3) = g.lmod.lmod_sig(k_lmod,1) + pert;
        g.lmod.lmod_sig(k_lmod,4) = g.lmod.lmod_sig(k_lmod,1) - pert;

        % resetting the input
        run('p_file');
        g.lmod.lmod_sig(k_lmod,3) = g.lmod.lmod_sig(k_lmod,1);
        g.lmod.lmod_sig(k_lmod,4) = g.lmod.lmod_sig(k_lmod,1);
        c_state = 0;
    end
end

% disturb rlmod states
n_tot_rlmod = n_tot_lmod + g.rlmod.n_rlmod;
if (g.rlmod.n_rlmod ~= 0)
    disp('disturbing reactive load modulation')
    for k = n_tot_lmod+1:n_tot_rlmod
        k_rlmod = k - n_tot_lmod;

        j = 1;
        pert = p_ratio*abs(g.rlmod.rlmod_st(k_rlmod,1));
        pert = max(pert,p_ratio);
        g.rlmod.rlmod_st(k_rlmod,3) = g.rlmod.rlmod_st(k_rlmod,1) + pert;
        g.rlmod.rlmod_st(k_rlmod,4) = g.rlmod.rlmod_st(k_rlmod,1) - pert;

        run('p_file');
        st_name(k,j) = 36;

        % disturb the auxiliary input
        disp('disturbing rlmod_sig')
        c_state = 6;
        rlmod_input = k_rlmod;

        pert = p_ratio;
        g.rlmod.rlmod_sig(k_rlmod,3) = g.rlmod.rlmod_sig(k_rlmod,1) + pert;
        g.rlmod.rlmod_sig(k_rlmod,4) = g.rlmod.rlmod_sig(k_rlmod,1) - pert;

        % resetting the input
        run('p_file');
        g.rlmod.rlmod_sig(k_rlmod,3) = g.rlmod.rlmod_sig(k_rlmod,1);
        g.rlmod.rlmod_sig(k_rlmod,4) = g.rlmod.rlmod_sig(k_rlmod,1);
        c_state = 0;
    end
end

% disturb pwrmod_p states
n_tot_pwrmod_p = n_tot_rlmod + g.pwr.n_pwrmod;
if (g.pwr.n_pwrmod ~= 0)
    disp('disturbing real power modulation')
    for k = n_tot_rlmod+1:n_tot_pwrmod_p
        k_pwrmod_p = k - n_tot_rlmod;

        j = 1;
        pert = p_ratio*abs(g.pwr.pwrmod_p_st(k_pwrmod_p,1));
        pert = max(pert,p_ratio);
        g.pwr.pwrmod_p_st(k_pwrmod_p,3) = g.pwr.pwrmod_p_st(k_pwrmod_p,1) + pert;
        g.pwr.pwrmod_p_st(k_pwrmod_p,4) = g.pwr.pwrmod_p_st(k_pwrmod_p,1) - pert;

        run('p_file');
        st_name(k,j) = 37;

        % disturb the auxiliary input
        disp('disturbing pwrmod_p_sig')
        c_state = 7;
        pwrmod_p_input = k_pwrmod_p;

        pert = p_ratio;
        g.pwr.pwrmod_p_sig(k_pwrmod_p,3) = g.pwr.pwrmod_p_sig(k_pwrmod_p,1) + pert;
        g.pwr.pwrmod_p_sig(k_pwrmod_p,4) = g.pwr.pwrmod_p_sig(k_pwrmod_p,1) - pert;

        % resetting the input
        run('p_file');
        g.pwr.pwrmod_p_sig(k_pwrmod_p,3) = g.pwr.pwrmod_p_sig(k_pwrmod_p,1);
        g.pwr.pwrmod_p_sig(k_pwrmod_p,4) = g.pwr.pwrmod_p_sig(k_pwrmod_p,1);
        c_state = 0;
    end
end

% disturb pwrmod_q states
n_tot_pwrmod_q = n_tot_pwrmod_p + g.pwr.n_pwrmod;
if (g.pwr.n_pwrmod ~= 0)
    disp('disturbing reactive power modulation')
    for k = n_tot_pwrmod_p+1:n_tot_pwrmod_q
        k_pwrmod_q = k - n_tot_pwrmod_p;

        j = 1;
        pert = p_ratio*abs(g.pwr.pwrmod_q_st(k_pwrmod_q,1));
        pert = max(pert,p_ratio);
        g.pwr.pwrmod_q_st(k_pwrmod_q,3) = g.pwr.pwrmod_q_st(k_pwrmod_q,1) + pert;
        g.pwr.pwrmod_q_st(k_pwrmod_q,4) = g.pwr.pwrmod_q_st(k_pwrmod_q,1) - pert;

        run('p_file');
        st_name(k,j) = 38;

        % disturb the auxiliary input
        disp('disturbing pwrmod_q_sig')
        c_state = 8;
        pwrmod_q_input = k_pwrmod_q;

        pert = p_ratio;
        g.pwr.pwrmod_q_sig(k_pwrmod_q,3) = g.pwr.pwrmod_q_sig(k_pwrmod_q,1) + pert;
        g.pwr.pwrmod_q_sig(k_pwrmod_q,4) = g.pwr.pwrmod_q_sig(k_pwrmod_q,1) - pert;

        % resetting the input
        run('p_file');
        g.pwr.pwrmod_q_sig(k_pwrmod_q,3) = g.pwr.pwrmod_q_sig(k_pwrmod_q,1);
        g.pwr.pwrmod_q_sig(k_pwrmod_q,4) = g.pwr.pwrmod_q_sig(k_pwrmod_q,1);
        c_state = 0;
    end
end

% disturb the hvdc states
n_tot_dcl = n_tot_pwrmod_q + g.dc.n_dcl;
if (g.dc.n_conv ~= 0)
    disp('disturbing hvdc')
    for k = n_tot_pwrmod_q+1:n_tot_dcl
        k_hvdc = k - n_tot_pwrmod_q;

        j = 1;
        pert = p_ratio*abs(g.dc.v_conr(k_hvdc,1));
        pert = max(pert,p_ratio);
        g.dc.v_conr(k_hvdc,3) = g.dc.v_conr(k_hvdc,1) + pert;
        g.dc.v_conr(k_hvdc,4) = g.dc.v_conr(k_hvdc,1) - pert;

        run('p_file');
        st_name(k,j) = 39;

        j = j + 1;
        pert = p_ratio*abs(g.dc.v_coni(k_hvdc,1));
        pert = max(pert,p_ratio);
        g.dc.v_coni(k_hvdc,3) = g.dc.v_coni(k_hvdc,1) + pert;
        g.dc.v_coni(k_hvdc,4) = g.dc.v_coni(k_hvdc,1) - pert;

        run('p_file');
        st_name(k,j) = 40;

        j = j + 1;
        pert = p_ratio*abs(g.dc.i_dcr(k_hvdc,1));
        pert = max(pert,p_ratio);
        g.dc.i_dcr(k_hvdc,3) = g.dc.i_dcr(k_hvdc,1) + pert;
        g.dc.i_dcr(k_hvdc,4) = g.dc.i_dcr(k_hvdc,1) - pert;

        run('p_file');
        st_name(k,j) = 41;

        if ~isempty(g.dc.cap_idx)
            k_cap_idx = find(g.dc.cap_idx == k_hvdc);
            if ~isempty(k_cap_idx)
                pert = p_ratio*abs(g.dc.i_dci(k_hvdc,1));
                pert = max(pert,p_ratio);
                g.dc.i_dci(k_hvdc,3) = g.dc.i_dci(k_hvdc,1) + pert;
                g.dc.i_dci(k_hvdc,4) = g.dc.i_dci(k_hvdc,1) - pert;

                run('p_file');
                st_name(k,j) = 42;

                j = j + 1;
                pert = p_ratio*abs(g.dc.v_dcc(k_hvdc,1));
                pert = max(pert,p_ratio);
                g.dc.v_dcc(k_hvdc,3) = g.dc.v_dcc(k_hvdc,1) + pert;
                g.dc.v_dcc(k_hvdc,4) = g.dc.v_dcc(k_hvdc,1) - pert;

                run('p_file');
                st_name(k,j) = 43;
            end
        end

        % disturb the auxiliary input
        disp('disturbing rectifier dc_sig')
        c_state = 9;
        dcmod_input = k_hvdc;

        pert = p_ratio;
        g.dc.dc_sig(g.dc.r_idx(k_hvdc),3) = g.dc.dc_sig(g.dc.r_idx(k_hvdc),1) + pert;
        g.dc.dc_sig(g.dc.r_idx(k_hvdc),4) = g.dc.dc_sig(g.dc.r_idx(k_hvdc),1) - pert;

        % resetting the input
        run('p_file');
        g.dc.dc_sig(g.dc.r_idx(k_hvdc),3) = g.dc.dc_sig(g.dc.r_idx(k_hvdc),1);
        g.dc.dc_sig(g.dc.r_idx(k_hvdc),4) = g.dc.dc_sig(g.dc.r_idx(k_hvdc),1);
        c_state = 0;

        % disturb the auxiliary input
        disp('disturbing inverter dc_sig')
        c_state = 10;
        dcmod_input = k_hvdc;

        pert = p_ratio;
        g.dc.dc_sig(g.dc.i_idx(k_hvdc),3) = g.dc.dc_sig(g.dc.i_idx(k_hvdc),1) + pert;
        g.dc.dc_sig(g.dc.i_idx(k_hvdc),4) = g.dc.dc_sig(g.dc.i_idx(k_hvdc),1) - pert;

        % resetting the input
        run('p_file');
        g.dc.dc_sig(g.dc.i_idx(k_hvdc),3) = g.dc.dc_sig(g.dc.i_idx(k_hvdc),1);
        g.dc.dc_sig(g.dc.i_idx(k_hvdc),4) = g.dc.dc_sig(g.dc.i_idx(k_hvdc),1);
        c_state = 0;
    end
end

% disturb the ess states
n_tot_ess = n_tot_dcl + g.ess.n_ess;
if (g.ess.n_ess ~= 0)
    disp('disturbing ess')
    for k = n_tot_dcl+1:n_tot_ess
        k_ess = k - n_tot_dcl;

        % transducer for local voltage magnitude (non-bypassable)
        j = 1;
        pert = p_ratio*abs(g.ess.ess1(k_ess,1));
        pert = max(pert,p_ratio);
        g.ess.ess1(k_ess,3) = g.ess.ess1(k_ess,1) + pert;
        g.ess.ess1(k_ess,4) = g.ess.ess1(k_ess,1) - pert;

        run('p_file');
        st_name(k,j) = 44;
        j = j + 1;

        % local voltage magnitude Pade approximation
        if ((g.ess.ess_con(k_ess,4) ~= 0) && (g.ess.ess_con(k_ess,5)/2 >= lbnd))
            pert = p_ratio*abs(g.ess.ess2(k_ess,1));
            pert = max(pert,p_ratio);
            g.ess.ess2(k_ess,3) = g.ess.ess2(k_ess,1) + pert;
            g.ess.ess2(k_ess,4) = g.ess.ess2(k_ess,1) - pert;

            run('p_file');
            st_name(k,j) = 45;
            j = j + 1;
        end

        % active current converter interface (non-bypassable)
        pert = p_ratio*abs(g.ess.ess3(k_ess,1));
        pert = max(pert,p_ratio);
        g.ess.ess3(k_ess,3) = g.ess.ess3(k_ess,1) + pert;
        g.ess.ess3(k_ess,4) = g.ess.ess3(k_ess,1) - pert;

        run('p_file');
        st_name(k,j) = 46;
        j = j + 1;

        % reactive current interface for potential future expansion
        pert = p_ratio*abs(g.ess.ess4(k_ess,1));
        pert = max(pert,p_ratio);
        g.ess.ess4(k_ess,3) = g.ess.ess4(k_ess,1) + pert;
        g.ess.ess4(k_ess,4) = g.ess.ess4(k_ess,1) - pert;

        run('p_file');
        st_name(k,j) = 47;
        j = j + 1;

      % % state of charge integrator (no effect on linearization)
      % pert = p_ratio*abs(g.ess.ess5(k_ess,1));
      % pert = max(pert,p_ratio);
      % g.ess.ess5(k_ess,3) = g.ess.ess5(k_ess,1) + pert;
      % g.ess.ess5(k_ess,4) = g.ess.ess5(k_ess,1) - pert;
      %
      % run('p_file');
      % st_name(k,j) = 48;
      % j = j + 1;

        %--------------------------------------------------%
        % disturbing the converter interface input

        disp('disturbing the real part of ess_sig (active)')
        c_state = 11;

        ess_input = k_ess;
        pert = p_ratio;
        g.ess.ess_sig(k_ess,3) = g.ess.ess_sig(k_ess,1) + pert;
        g.ess.ess_sig(k_ess,4) = g.ess.ess_sig(k_ess,1) - pert;

        % resetting the input
        run('p_file');
        g.ess.ess_sig(k_ess,3) = g.ess.ess_sig(k_ess,1);
        g.ess.ess_sig(k_ess,4) = g.ess.ess_sig(k_ess,1);
        c_state = 0;

        disp('disturbing the imaginary part of ess_sig (reactive)')
        c_state = 12;

        ess_input = k_ess;
        pert = p_ratio;
        g.ess.ess_sig(k_ess,3) = g.ess.ess_sig(k_ess,1) + 1j*pert;
        g.ess.ess_sig(k_ess,4) = g.ess.ess_sig(k_ess,1) - 1j*pert;

        % resetting the input
        run('p_file');
        g.ess.ess_sig(k_ess,3) = g.ess.ess_sig(k_ess,1);
        g.ess.ess_sig(k_ess,4) = g.ess.ess_sig(k_ess,1);
        c_state = 0;
    end
end

% disturb the lsc states
n_tot_lsc = n_tot_ess + g.lsc.n_lsc;
if (g.lsc.n_lsc ~= 0)
    disp('disturbing lsc')
    for k = n_tot_ess+1:n_tot_lsc
        k_lsc = k - n_tot_ess;

        % transducer for remote angle signal (non-bypassable)
        j = 1;
        pert = p_ratio*abs(g.lsc.lsc1(k_lsc,1));
        pert = max(pert,p_ratio);
        g.lsc.lsc1(k_lsc,3) = g.lsc.lsc1(k_lsc,1) + pert;
        g.lsc.lsc1(k_lsc,4) = g.lsc.lsc1(k_lsc,1) - pert;

        run('p_file');
        st_name(k,j) = 49;
        j = j + 1;

        % transducer for local angle signal (non-bypassable)
        pert = p_ratio*abs(g.lsc.lsc2(k_lsc,1));
        pert = max(pert,p_ratio);
        g.lsc.lsc2(k_lsc,3) = g.lsc.lsc2(k_lsc,1) + pert;
        g.lsc.lsc2(k_lsc,4) = g.lsc.lsc2(k_lsc,1) - pert;

        run('p_file');
        st_name(k,j) = 50;
        j = j + 1;

        % remote angle Pade approximation
        if ((g.lsc.lsc_con(k_lsc,4) ~= 0) && (g.lsc.lsc_con(k_lsc,5)/2 >= lbnd))
            pert = p_ratio*abs(g.lsc.lsc3(k_lsc,1));
            pert = max(pert,p_ratio);
            g.lsc.lsc3(k_lsc,3) = g.lsc.lsc3(k_lsc,1) + pert;
            g.lsc.lsc3(k_lsc,4) = g.lsc.lsc3(k_lsc,1) - pert;

            run('p_file');
            st_name(k,j) = 51;
            j = j + 1;
        end

        % local angle Pade approximation
        if ((g.lsc.lsc_con(k_lsc,4) ~= 0) && (g.lsc.lsc_con(k_lsc,6)/2 >= lbnd))
            pert = p_ratio*abs(g.lsc.lsc4(k_lsc,1));
            pert = max(pert,p_ratio);
            g.lsc.lsc4(k_lsc,3) = g.lsc.lsc4(k_lsc,1) + pert;
            g.lsc.lsc4(k_lsc,4) = g.lsc.lsc4(k_lsc,1) - pert;

            run('p_file');
            st_name(k,j) = 52;
            j = j + 1;
        end

        % remote angle setpoint tracking filter
        if (g.lsc.lsc_con(k_lsc,8) >= lbnd)
            pert = p_ratio*abs(g.lsc.lsc5(k_lsc,1));
            pert = max(pert,p_ratio);
            g.lsc.lsc5(k_lsc,3) = g.lsc.lsc5(k_lsc,1) + pert;
            g.lsc.lsc5(k_lsc,4) = g.lsc.lsc5(k_lsc,1) - pert;

            run('p_file');
            st_name(k,j) = 53;
            j = j + 1;
        end

        % local angle setpoint tracking filter
        if (g.lsc.lsc_con(k_lsc,10) >= lbnd)
            pert = p_ratio*abs(g.lsc.lsc6(k_lsc,1));
            pert = max(pert,p_ratio);
            g.lsc.lsc6(k_lsc,3) = g.lsc.lsc6(k_lsc,1) + pert;
            g.lsc.lsc6(k_lsc,4) = g.lsc.lsc6(k_lsc,1) - pert;

            run('p_file');
            st_name(k,j) = 54;
            j = j + 1;
        end

        % local LTV highpass filter state 1 (non-bypassable)
        pert = p_ratio*abs(g.lsc.lsc7(k_lsc,1));
        pert = max(pert,p_ratio);
        g.lsc.lsc7(k_lsc,3) = g.lsc.lsc7(k_lsc,1) + pert;
        g.lsc.lsc7(k_lsc,4) = g.lsc.lsc7(k_lsc,1) - pert;

        run('p_file');
        st_name(k,j) = 55;
        j = j + 1;

        % local LTV highpass filter state 2 (non-bypassable)
        pert = p_ratio*abs(g.lsc.lsc8(k_lsc,1));
        pert = max(pert,p_ratio);
        g.lsc.lsc8(k_lsc,3) = g.lsc.lsc8(k_lsc,1) + pert;
        g.lsc.lsc8(k_lsc,4) = g.lsc.lsc8(k_lsc,1) - pert;

        run('p_file');
        st_name(k,j) = 56;
        j = j + 1;

        % local LTV lead-lag stage 1
        if (g.lsc.lsc_con(k_lsc,17) >= lbnd)
            pert = p_ratio*abs(g.lsc.lsc9(k_lsc,1));
            pert = max(pert,p_ratio);
            g.lsc.lsc9(k_lsc,3) = g.lsc.lsc9(k_lsc,1) + pert;
            g.lsc.lsc9(k_lsc,4) = g.lsc.lsc9(k_lsc,1) - pert;

            run('p_file');
            st_name(k,j) = 57;
            j = j + 1;
        end

        % local LTV lead-lag stage 2
        if (g.lsc.lsc_con(k_lsc,19) >= lbnd)
            pert = p_ratio*abs(g.lsc.lsc10(k_lsc,1));
            pert = max(pert,p_ratio);
            g.lsc.lsc10(k_lsc,3) = g.lsc.lsc10(k_lsc,1) + pert;
            g.lsc.lsc10(k_lsc,4) = g.lsc.lsc10(k_lsc,1) - pert;

            run('p_file');
            st_name(k,j) = 58;
            j = j + 1;
        end

        % center LTI highpass filter state 1 (non-bypassable)
        pert = p_ratio*abs(g.lsc.lsc11(k_lsc,1));
        pert = max(pert,p_ratio);
        g.lsc.lsc11(k_lsc,3) = g.lsc.lsc11(k_lsc,1) + pert;
        g.lsc.lsc11(k_lsc,4) = g.lsc.lsc11(k_lsc,1) - pert;

        run('p_file');
        st_name(k,j) = 59;
        j = j + 1;

        % center LTI highpass filter state 2 (non-bypassable)
        pert = p_ratio*abs(g.lsc.lsc12(k_lsc,1));
        pert = max(pert,p_ratio);
        g.lsc.lsc12(k_lsc,3) = g.lsc.lsc12(k_lsc,1) + pert;
        g.lsc.lsc12(k_lsc,4) = g.lsc.lsc12(k_lsc,1) - pert;

        run('p_file');
        st_name(k,j) = 60;
        j = j + 1;

        % center LTI lead-lag stage 1
        if (g.lsc.lsc_con(k_lsc,28) >= lbnd)
            pert = p_ratio*abs(g.lsc.lsc13(k_lsc,1));
            pert = max(pert,p_ratio);
            g.lsc.lsc13(k_lsc,3) = g.lsc.lsc13(k_lsc,1) + pert;
            g.lsc.lsc13(k_lsc,4) = g.lsc.lsc13(k_lsc,1) - pert;

            run('p_file');
            st_name(k,j) = 61;
            j = j + 1;
        end

        % center LTI lead-lag stage 2
        if (g.lsc.lsc_con(k_lsc,30) >= lbnd)
            pert = p_ratio*abs(g.lsc.lsc14(k_lsc,1));
            pert = max(pert,p_ratio);
            g.lsc.lsc14(k_lsc,3) = g.lsc.lsc14(k_lsc,1) + pert;
            g.lsc.lsc14(k_lsc,4) = g.lsc.lsc14(k_lsc,1) - pert;

            run('p_file');
            st_name(k,j) = 62;
            j = j + 1;
        end

        % lowpass filter (non-bypassable)
        pert = p_ratio*abs(g.lsc.lsc15(k_lsc,1));
        pert = max(pert,p_ratio);
        g.lsc.lsc15(k_lsc,3) = g.lsc.lsc15(k_lsc,1) + pert;
        g.lsc.lsc15(k_lsc,4) = g.lsc.lsc15(k_lsc,1) - pert;

        run('p_file');
        st_name(k,j) = 63;
        j = j + 1;
    end
end

% disturb the reec states
n_tot_reec = n_tot_lsc + g.reec.n_reec;
if (g.reec.n_reec ~= 0)
    disp('disturbing reec')
    for k = n_tot_lsc+1:n_tot_reec
        k_reec = k - n_tot_lsc;

        % reec1 -- voltage magnitude transducer
        j = 1;
        pert = p_ratio*abs(g.reec.reec1(k_reec,1));
        pert = max(pert,p_ratio);
        g.reec.reec1(k_reec,3) = g.reec.reec1(k_reec,1) + pert;
        g.reec.reec1(k_reec,4) = g.reec.reec1(k_reec,1) - pert;

        run('p_file');
        st_name(k,j) = 64;
        j = j + 1;

        % reec2 -- active power transducer is present if pfflag == 1
        if (g.reec.reec_con(k_reec,32) == 1)
            pert = p_ratio*abs(g.reec.reec2(k_reec,1));
            pert = max(pert,p_ratio);
            g.reec.reec2(k_reec,3) = g.reec.reec2(k_reec,1) + pert;
            g.reec.reec2(k_reec,4) = g.reec.reec2(k_reec,1) - pert;

            run('p_file');
            st_name(k,j) = 65;
            j = j + 1;
        end

        % reec3 -- reactive power transducer
        pert = p_ratio*abs(g.reec.reec3(k_reec,1));
        pert = max(pert,p_ratio);
        g.reec.reec3(k_reec,3) = g.reec.reec3(k_reec,1) + pert;
        g.reec.reec3(k_reec,4) = g.reec.reec3(k_reec,1) - pert;

        run('p_file');
        st_name(k,j) = 66;
        j = j + 1;

        % reec4 -- first stage PI loop
        if ((g.reec.reec_con(k_reec,33) == 1) && (g.reec.reec_con(k_reec,34) == 1))
            if (g.reec.reec_con(k_reec,22) > 0)
                pert = p_ratio*abs(g.reec.reec4(k_reec,1));
                pert = max(pert,p_ratio);
                g.reec.reec4(k_reec,3) = g.reec.reec4(k_reec,1) + pert;
                g.reec.reec4(k_reec,4) = g.reec.reec4(k_reec,1) - pert;

                run('p_file');
                st_name(k,j) = 67;
                j = j + 1;
            end
        end

        % reec5 -- second stage PI loop
        if ((g.reec.reec_con(k_reec,34) == 1) && (g.reec.reec_con(k_reec,24) > 0))
            pert = p_ratio*abs(g.reec.reec5(k_reec,1));
            pert = max(pert,p_ratio);
            g.reec.reec5(k_reec,3) = g.reec.reec5(k_reec,1) + pert;
            g.reec.reec5(k_reec,4) = g.reec.reec5(k_reec,1) - pert;

            run('p_file');
            st_name(k,j) = 68;
            j = j + 1;
        end

        % reec6 -- reactive current filter
        if (g.reec.reec_con(k_reec,34) == 0) && (g.reec.reec_con(k_reec,25) >= lbnd)
            pert = p_ratio*abs(g.reec.reec6(k_reec,1));
            pert = max(pert,p_ratio);
            g.reec.reec6(k_reec,3) = g.reec.reec6(k_reec,1) + pert;
            g.reec.reec6(k_reec,4) = g.reec.reec6(k_reec,1) - pert;

            run('p_file');
            st_name(k,j) = 69;
            j = j + 1;
        end

        % reec7 -- active power reference filter
        pert = p_ratio*abs(g.reec.reec7(k_reec,1));
        pert = max(pert,p_ratio);
        g.reec.reec7(k_reec,3) = g.reec.reec7(k_reec,1) + pert;
        g.reec.reec7(k_reec,4) = g.reec.reec7(k_reec,1) - pert;

        run('p_file');
        st_name(k,j) = 70;
        j = j + 1;

        % reec8 -- compensated voltage filter
        if (g.reec.reec_con(k_reec,38) >= lbnd)
            pert = p_ratio*abs(g.reec.reec8(k_reec,1));
            pert = max(pert,p_ratio);
            g.reec.reec8(k_reec,3) = g.reec.reec8(k_reec,1) + pert;
            g.reec.reec8(k_reec,4) = g.reec.reec8(k_reec,1) - pert;

            run('p_file');
            st_name(k,j) = 71;
            j = j + 1;
        end

        % reec9 -- active current command filter
        pert = p_ratio*abs(g.reec.reec9(k_reec,1));
        pert = max(pert,p_ratio);
        g.reec.reec9(k_reec,3) = g.reec.reec9(k_reec,1) + pert;
        g.reec.reec9(k_reec,4) = g.reec.reec9(k_reec,1) - pert;

        run('p_file');
        st_name(k,j) = 72;
        j = j + 1;

        % reec10 -- reactive current command filter
        pert = p_ratio*abs(g.reec.reec10(k_reec,1));
        pert = max(pert,p_ratio);
        g.reec.reec10(k_reec,3) = g.reec.reec10(k_reec,1) + pert;
        g.reec.reec10(k_reec,4) = g.reec.reec10(k_reec,1) - pert;

        run('p_file');
        st_name(k,j) = 73;
        j = j + 1;

        %--------------------------------------------------%
        % disturbing the reec control references

        disp('disturbing the real part of reec_sig (voltage)')
        c_state = 13;

        reec_input = k_reec;
        pert = p_ratio;
        g.reec.reec_sig(k_reec,3) = g.reec.reec_sig(k_reec,1) + pert;
        g.reec.reec_sig(k_reec,4) = g.reec.reec_sig(k_reec,1) - pert;

        % resetting the input
        run('p_file');
        g.reec.reec_sig(k_reec,3) = g.reec.reec_sig(k_reec,1);
        g.reec.reec_sig(k_reec,4) = g.reec.reec_sig(k_reec,1);
        c_state = 0;

        disp('disturbing the imaginary part of reec_sig (reactive)')
        c_state = 14;

        reec_input = k_reec;
        pert = p_ratio;
        g.reec.reec_sig(k_reec,3) = g.reec.reec_sig(k_reec,1) + 1j*pert;
        g.reec.reec_sig(k_reec,4) = g.reec.reec_sig(k_reec,1) - 1j*pert;

        % resetting the input
        run('p_file');
        g.reec.reec_sig(k_reec,3) = g.reec.reec_sig(k_reec,1);
        g.reec.reec_sig(k_reec,4) = g.reec.reec_sig(k_reec,1);
        c_state = 0;
    end
end

% disturb the gfma states
n_tot_gfma = n_tot_reec + g.gfma.n_gfma;
if (g.gfma.n_gfma ~= 0)
    disp('disturbing gfma')
    for k = n_tot_reec+1:n_tot_gfma
        k_gfma = k - n_tot_reec;

        % gfma1 -- commanded voltage angle integrator
        j = 1;
        pert = p_ratio*abs(g.gfma.gfma1(k_gfma,1));
        pert = max(pert,p_ratio);
        g.gfma.gfma1(k_gfma,3) = g.gfma.gfma1(k_gfma,1) + pert;
        g.gfma.gfma1(k_gfma,4) = g.gfma.gfma1(k_gfma,1) - pert;

        run('p_file');
        st_name(k,j) = 74;
        j = j + 1;

        % gfma2 and gfma3 states not included in linearization

        % gfma4 -- commanded voltage magnitude (first-order)
        pert = p_ratio*abs(g.gfma.gfma4(k_gfma,1));
        pert = max(pert,p_ratio);
        g.gfma.gfma4(k_gfma,3) = g.gfma.gfma4(k_gfma,1) + pert;
        g.gfma.gfma4(k_gfma,4) = g.gfma.gfma4(k_gfma,1) - pert;

        run('p_file');
        st_name(k,j) = 75;
        j = j + 1;

        % gfma5 -- voltage regulation integral control state (vflag == 1)
        if ((g.gfma.gfma_con(k_gfma,25) == 1) && (g.gfma.gfma_con(k_gfma,13) > 0))
            pert = p_ratio*abs(g.gfma.gfma5(k_gfma,1));
            pert = max(pert,p_ratio);
            g.gfma.gfma5(k_gfma,3) = g.gfma.gfma5(k_gfma,1) + pert;
            g.gfma.gfma5(k_gfma,4) = g.gfma.gfma5(k_gfma,1) - pert;

            run('p_file');
            st_name(k,j) = 76;
            j = j + 1;
        end

        % gfma6 and gfma7 states not included in linearization

        % gfma8 -- active power transducer
        pert = p_ratio*abs(g.gfma.gfma8(k_gfma,1));
        pert = max(pert,p_ratio);
        g.gfma.gfma8(k_gfma,3) = g.gfma.gfma8(k_gfma,1) + pert;
        g.gfma.gfma8(k_gfma,4) = g.gfma.gfma8(k_gfma,1) - pert;

        run('p_file');
        st_name(k,j) = 77;
        j = j + 1;

        % gfma9 -- reactive power transducer
        pert = p_ratio*abs(g.gfma.gfma9(k_gfma,1));
        pert = max(pert,p_ratio);
        g.gfma.gfma9(k_gfma,3) = g.gfma.gfma9(k_gfma,1) + pert;
        g.gfma.gfma9(k_gfma,4) = g.gfma.gfma9(k_gfma,1) - pert;

        run('p_file');
        st_name(k,j) = 78;
        j = j + 1;

        % gfma10 -- voltage transducer
        pert = p_ratio*abs(g.gfma.gfma10(k_gfma,1));
        pert = max(pert,p_ratio);
        g.gfma.gfma10(k_gfma,3) = g.gfma.gfma10(k_gfma,1) + pert;
        g.gfma.gfma10(k_gfma,4) = g.gfma.gfma10(k_gfma,1) - pert;

        run('p_file');
        st_name(k,j) = 79;
        j = j + 1;

        %--------------------------------------------------%
        % disturbing the gfma control references

        disp('disturbing the real part of gfma_sig (active)')
        c_state = 15;

        gfma_input = k_gfma;
        pert = p_ratio;
        g.gfma.gfma_sig(k_gfma,3) = g.gfma.gfma_sig(k_gfma,1) + pert;
        g.gfma.gfma_sig(k_gfma,4) = g.gfma.gfma_sig(k_gfma,1) - pert;

        % resetting the input
        run('p_file');
        g.gfma.gfma_sig(k_gfma,3) = g.gfma.gfma_sig(k_gfma,1);
        g.gfma.gfma_sig(k_gfma,4) = g.gfma.gfma_sig(k_gfma,1);
        c_state = 0;

        disp('disturbing the imaginary part of gfma_sig (reactive)')
        c_state = 16;

        gfma_input = k_gfma;
        pert = p_ratio;
        g.gfma.gfma_sig(k_gfma,3) = g.gfma.gfma_sig(k_gfma,1) + 1j*pert;
        g.gfma.gfma_sig(k_gfma,4) = g.gfma.gfma_sig(k_gfma,1) - 1j*pert;

        % resetting the input
        run('p_file');
        g.gfma.gfma_sig(k_gfma,3) = g.gfma.gfma_sig(k_gfma,1);
        g.gfma.gfma_sig(k_gfma,4) = g.gfma.gfma_sig(k_gfma,1);
        c_state = 0;
    end
end

% n_tot_foo = n_tot_gfma + g.foo.n_foo;  % template for future expansion

% eof
