% This m-file checks the dynamic data indices and determines the
% total number of states for each device
%
% Called from svm_mgen

%-----------------------------------------------------------------------------%
% Version history
%
% Author:   Ryan Elliott
% Modified: August 2019
% Note:     ess model added
%
% Author:   D. Trudnowski
% Modified: Added pwrmod
% Date:     2015
%
% Author:   Graham Rogers
% Modified: tcsc model added
% Date:     December 1998
%
% Author:   Graham Rogers
% Modified: Induction Generators and load modulation added
% Date:     August 1997
%
% Version:  1.0 (initial version)
% Author:   Graham Rogers
% Date:     September 1996
%-----------------------------------------------------------------------------%

for k = 1:g.mac.n_mac
    no_mac = 0;
    if ~isempty(g.mac.mac_ib_idx)
        k_ib = find(g.mac.mac_ib_idx == k);
        if ~isempty(k_ib)
            state(k) = 0;
            no_mac = 1;
        end
    end

    % generators
    if (no_mac == 0)
        % classical
        if ~isempty(g.mac.mac_em_idx)
            k_em = find(g.mac.mac_em_idx == k);
            if ~isempty(k_em)
                state(k) = nss.mac_em;
            end
        end

        % transient
        if ~isempty(g.mac.mac_tra_idx)
            k_tra = find(g.mac.mac_tra_idx == k);
            if ~isempty(k_tra)
                state(k) = nss.mac_tra;
            end
        end

        % subtransient
        if ~isempty(g.mac.mac_sub_idx)
            k_sub = find(g.mac.mac_sub_idx == k);
            if ~isempty(k_sub)
                state(k) = nss.mac_sub;
            end
        end
    end

    gen_state(k) = state(k);

    % exciters
    if ~isempty(g.exc.exc_con)
        s_TR = 0;
        s_TB = 0;
        s_TA = 0;
        s_TE = 0;

        % simplified
        if (g.exc.n_smp ~= 0)
            k_smp = find(mac_exc(g.exc.smp_idx) == k);
            if ~isempty(k_smp)
                if ~isempty(g.exc.smp_TR_idx)
                    s_TR = sum(g.exc.smp_TR_idx == k_smp);
                end
                if ~isempty(g.exc.smp_TB_idx)
                    s_TB = sum(g.exc.smp_TB_idx == k_smp);
                end
                if ~isempty(g.exc.smp_TA_idx)
                    s_TA = sum(g.exc.smp_TA_idx == k_smp);
                end

                TR_state(k) = s_TR;
                TB_state(k) = s_TB;
                Efd_state(k) = s_TA;
                state(k) = state(k) + s_TR + s_TB + s_TA;
            end
        end

        % simplified proportional/integral
        if (g.exc.n_smppi ~= 0)
            k_smppi = find(mac_exc(g.exc.smppi_idx) == k);
            if ~isempty(k_smppi)
                if ~isempty(g.exc.smppi_TR_idx)
                    s_TR = sum(g.exc.smppi_TR_idx == k_smppi);
                end

                TR_state(k) = s_TR;
                TB_state(k) = 1;
                Efd_state(k) = 1;
                state(k) = state(k) + s_TR + 2;
            end
        end

        % ieee type dc1 and dc2 exciters
        if (g.exc.n_dc ~= 0)
            k_dc = find(mac_exc(g.exc.dc_idx) == k);
            if ~isempty(k_dc)
                if ~isempty(g.exc.dc_TR_idx)
                    s_TR = sum(g.exc.dc_TR_idx == k_dc);
                end
                if ~isempty(g.exc.dc_TB_idx)
                    s_TB = sum(g.exc.dc_TB_idx == k_dc);
                end
                if ~isempty(g.exc.dc_TA_idx)
                    s_TA = sum(g.exc.dc_TA_idx == k_dc);
                end
                if ~isempty(g.exc.dc_TE_idx)
                    s_TE = sum(g.exc.dc_TE_idx == k_dc);
                end

                TR_state(k) = s_TR;
                TB_state(k) = s_TB;
                TA_state(k) = s_TA;
                Efd_state(k) = s_TE;
                R_f_state(k) = 1;

                state(k) = state(k) + s_TR + s_TB + s_TA + s_TE + 1;
            end
        end

        % ieee type st3 compound source rectifier
        if (g.exc.n_st3 ~= 0)
            k_st3 = find(mac_exc(g.exc.st3_idx) == k);
            if ~isempty(k_st3)
                if ~isempty(g.exc.st3_TR_idx)
                    s_TR = sum(g.exc.st3_TR_idx == k_st3);
                end
                if ~isempty(g.exc.st3_TB_idx)
                    s_TB = sum(g.exc.st3_TB_idx == k_st3);
                end
                if ~isempty(g.exc.st3_TA_idx)
                    s_TA = sum(g.exc.st3_TA_idx == k_st3);
                end

                TR_state(k) = s_TR;
                TB_state(k) = s_TB;
                TA_state(k) = s_TA;

                state(k) = state(k) + s_TR + s_TB + s_TA;
            end
        end
    end

    % pss
    if ~isempty(g.pss.pss_con)
        s_T4 = 0;
        k_pss = find(mac_pss == k);
        if ~isempty(k_pss)
            if ~isempty(g.pss.pss_T4_idx)
                s_T4 = sum(g.pss.pss_T4_idx == k_pss);
            else
                s_T4 = 0;
            end

            pss1_state(k) = 1;
            pss2_state(k) = 1;
            pss3_state(k) = s_T4;
            state(k) = state(k) + s_T4 + 2;
        end
    end

    % deltaP/omega filter
    if ~isempty(g.dpw.dpw_con)
        k_dpw = find(mac_dpw == k);
        if ~isempty(k_dpw)
            state(k) = state(k) + nss.dpw;
            dpw_state(k) = 1;
        end
    end

    % governors
    if ~isempty(g.tg.tg_con)
        % conventional
        k_tg = find(mac_tg == k);
        if ~isempty(k_tg)
            state(k) = state(k) + nss.tg;
            tg_state(k) = 1;
        end

        % hydro
        k_tgh = find(mac_tgh == k);
        if ~isempty(k_tgh)
            state(k) = state(k) + nss.tgh;
            tgh_state(k) = 1;
        end
    end
end

% induction motor
n_mot_states = 0;
if (g.ind.n_mot ~= 0)
    state_mot = nss.mot*ones(g.ind.n_mot,1);
    n_mot_states = sum(state_mot);
    state(g.mac.n_mac+1:n_gm) = state_mot;
end

% induction generator
n_ig_states = 0;
if (g.igen.n_ig ~= 0)
    state_ig = nss.ig*ones(g.igen.n_ig,1);
    n_ig_states = sum(state_ig);
    state(n_gm+1:n_tot) = state_ig;
end

% svc
n_svc_states = 0;
n_svc1 = n_tot;
if (g.svc.n_svc ~= 0)
    state_svc = ones(g.svc.n_svc,1);
    state_svccon = 0;
    sc = zeros(g.svc.n_svc,1);

    if ~isempty(g.svc.svcll_idx)
        state_svccon = ones(length(g.svc.svcll_idx),1);
        sc(g.svc.svcll_idx) = state_svccon;
    end

    n_svc_states = sum(state_svc) + sum(state_svccon);
    state(n_svc1+1:n_svc1+g.svc.n_svc) = ones(g.svc.n_svc,1) + sc;
end

% tcsc
n_tcsc_states = 0;
n_tcsc1 = n_svc1 + g.svc.n_svc;
if (g.tcsc.n_tcsc ~= 0)
    state_tcsc = nss.tcsc*ones(g.tcsc.n_tcsc,1);
    n_tcsc_states = sum(state_tcsc);
    state(n_tcsc1+1:n_tcsc1+g.tcsc.n_tcsc) = state_tcsc;
end

% lmod
n_lmod_states = 0;
n_lmod1 = n_tcsc1 + g.tcsc.n_tcsc;
if (g.lmod.n_lmod ~= 0)
    state_lmod = nss.lmod*ones(g.lmod.n_lmod,1);
    n_lmod_states = sum(state_lmod);
    state(n_lmod1+1:n_lmod1+g.lmod.n_lmod) = state_lmod;
end

% rlmod
n_rlmod_states = 0;
n_rlmod1 = n_lmod1 + g.lmod.n_lmod;
if (g.rlmod.n_rlmod ~= 0)
    state_rlmod = nss.rlmod*ones(g.rlmod.n_rlmod,1);
    n_rlmod_states = sum(state_rlmod);
    state(n_rlmod1+1:n_rlmod1+g.rlmod.n_rlmod) = state_rlmod;
end

% pwrmod
n_pwrmod_states = 0;
n_pwrmod_p1 = n_rlmod1 + g.rlmod.n_rlmod;
n_pwrmod_q1 = n_pwrmod_p1 + g.pwr.n_pwrmod;
if (g.pwr.n_pwrmod ~= 0)
    state_p_pwrmod = nss.pwrmod*ones(g.pwr.n_pwrmod,1);
    n_pwrmod_p_states = sum(state_p_pwrmod);
    state(n_pwrmod_p1+1:n_pwrmod_p1+g.pwr.n_pwrmod) = state_p_pwrmod;
    %
    state_q_pwrmod = nss.pwrmod*ones(g.pwr.n_pwrmod,1);
    n_pwrmod_q_states = sum(state_q_pwrmod);
    state(n_pwrmod_q1+1:n_pwrmod_q1+g.pwr.n_pwrmod) = state_q_pwrmod;
end

% hvdc
n_hvdc_states = 0;
n_hvdc1 = n_pwrmod_q1 + g.pwr.n_pwrmod;
if (g.dc.n_conv ~= 0)
    state_hvdc(1:g.dc.n_dcl) = nss.dcl*ones(g.dc.n_dcl,1);
    if (g.dc.l_cap ~= 0)
        state_hvdc(g.dc.cap_idx) = state_hvdc(g.dc.cap_idx) ...
                                   + nss.l_cap*ones(g.dc.l_cap,1);
    end

    n_hvdc_states = sum(state_hvdc);
    state(n_hvdc1+1:n_hvdc1+g.dc.n_dcl) = (nss.dcl + nss.l_cap*g.dc.l_cap) ...
                                          *ones(g.dc.n_dcl,1);
end

% ess
n_ess_states = 0;
n_ess1 = n_hvdc1 + g.dc.n_dcl;
if (g.ess.n_ess ~= 0)
    state_ess = zeros(g.ess.n_ess,1);
    for kk = 1:g.ess.n_ess
        state_ess(kk) = 3;                                  % non-bypassable states

        if (g.ess.ess_con(kk,3) < lbnd)
            estr = '\nsvm_mgen: ess voltage transducer time constant ';
            estr = [estr, 'must be nonzero at ess index %0.0f.'];
            error(sprintf(estr,kk));
        end

        if (g.ess.ess_con(kk,14) < lbnd)
            estr = '\nsvm_mgen: ess converter interface time constant ';
            estr = [estr, 'must be nonzero at ess index %0.0f.'];
            error(sprintf(estr,kk));
        end

        % local voltage magnitude Pade approximation
        if ((g.ess.ess_con(kk,4) ~= 0) && (g.ess.ess_con(kk,5)/2 >= lbnd))
            state_ess(kk) = state_ess(kk) + 1;
        end
    end

    n_ess_states = sum(state_ess);
    state(n_ess1+1:n_ess1+g.ess.n_ess) = state_ess;
end

% lsc
n_lsc_states = 0;
n_lsc1 = n_ess1 + g.ess.n_ess;
if (g.lsc.n_lsc ~= 0)
    state_lsc = zeros(g.lsc.n_lsc,1);
    for kk = 1:g.lsc.n_lsc
        state_lsc(kk) = 7;                                  % non-bypassable states

        if (g.lsc.lsc_con(kk,34) < lbnd)
            estr = '\nsvm_mgen: lsc lowpass filter time constant ';
            estr = [estr, 'must be nonzero at lsc index %0.0f.'];
            error(sprintf(estr,kk));
        end

        % remote angle Pade approximation
        if ((g.lsc.lsc_con(kk,4) ~= 0) && (g.lsc.lsc_con(kk,5)/2 >= lbnd))
            state_lsc(kk) = state_lsc(kk) + 1;
        end

        % local angle Pade approximation
        if ((g.lsc.lsc_con(kk,4) ~= 0) && (g.lsc.lsc_con(kk,6)/2 >= lbnd))
            state_lsc(kk) = state_lsc(kk) + 1;
        end

        for jj = [8,10,17,19,28,30]                         % bypassable states
            if (g.lsc.lsc_con(kk,jj) >= lbnd)
                state_lsc(kk) = state_lsc(kk) + 1;
            end
        end
    end

    n_lsc_states = sum(state_lsc);
    state(n_lsc1+1:n_lsc1+g.lsc.n_lsc) = state_lsc;
end

% eof
