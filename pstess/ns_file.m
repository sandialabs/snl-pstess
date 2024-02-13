% This m-file checks the dynamic data indices and determines the
% total number of states for each device
%
% Called from svm_mgen

%-----------------------------------------------------------------------------%
% Version history
%
% Author:   Ryan Elliott
% Modified: September 2023
% Note:     Added support for ivm, gfma, and reec models
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

        % internal voltage models (ivm)
        if ~isempty(g.mac.mac_ivm_idx)
            k_ivm = find(g.mac.mac_ivm_idx == k);
            if ~isempty(k_ivm)
                state(k) = nss.mac_ivm;

                % preventing mac_ang bypass
                if (g.mac.mac_con(k,9) < lbnd)
                    g.mac.mac_con(k,9) = 5*lbnd;
                    estr = '\nsvm_mgen: increasing ivm commanded angle time ';
                    estr = [estr, 'constant to be nonzero at ivm index %0.0f.'];
                    warning(sprintf(estr,k-length(g.mac.not_ivm_idx)));
                end

                % preventing eqprime bypass
                if (g.mac.mac_con(k,10) < lbnd)
                    g.mac.mac_con(k,10) = 5*lbnd;
                    estr = '\nsvm_mgen: increasing ivm commanded voltage time ';
                    estr = [estr, 'constant to be nonzero at ivm index %0.0f.'];
                    warning(sprintf(estr,k-length(g.mac.not_ivm_idx)));
                end
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
            g.ess.ess_con(kk,3) = 5*lbnd;
            estr = '\nsvm_mgen: increasing ess voltage transducer time ';
            estr = [estr, 'constant to be nonzero at ess index %0.0f.'];
            warning(sprintf(estr,kk));
        end

        if (g.ess.ess_con(kk,14) < lbnd)
            g.ess.ess_con(kk,14) = 5*lbnd;
            estr = '\nsvm_mgen: increasing ess converter interface time ';
            estr = [estr, 'constant to be nonzero at ess index %0.0f.'];
            warning(sprintf(estr,kk));
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

        if (g.lsc.lsc_con(kk,3) < lbnd)
            g.lsc.lsc_con(kk,3) = 5*lbnd;
            estr = '\nsvm_mgen: increasing lsc angle transducer time ';
            estr = [estr, 'constant to be nonzero at lsc index %0.0f.'];
            warning(sprintf(estr,kk));
        end

        if (g.lsc.lsc_con(kk,34) < lbnd)
            g.lsc.lsc_con(kk,34) = 5*lbnd;
            estr = '\nsvm_mgen: increasing lsc lowpass filter time ';
            estr = [estr, 'constant to be nonzero at lsc index %0.0f.'];
            warning(sprintf(estr,kk));
        end

        % remote angle Pade approximation
        if ((g.lsc.lsc_con(kk,4) == 1) && (g.lsc.lsc_con(kk,5)/2 >= lbnd))
            state_lsc(kk) = state_lsc(kk) + 1;
        end

        % local angle Pade approximation
        if ((g.lsc.lsc_con(kk,4) == 1) && (g.lsc.lsc_con(kk,6)/2 >= lbnd))
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

% reec
n_reec_states = 0;
n_reec1 = n_lsc1 + g.lsc.n_lsc;
if (g.reec.n_reec ~= 0)
    state_reec = zeros(g.reec.n_reec,1);
    for kk = 1:g.reec.n_reec
        state_reec(kk) = 5;                                 % non-bypassable states

        % eliminating deadbands from the local voltage control loop
        if ((g.reec.reec_con(kk,6) ~= 0) || (g.reec.reec_con(kk,7) ~= 0))
            g.reec.reec_con(kk,6) = 0;
            g.reec.reec_con(kk,7) = 0;
            estr = '\nsvm_mgen: eliminating local voltage control deadband ';
            estr = [estr, 'at reec index %0.0f.'];
            warning(sprintf(estr,kk));
        end

        % reec 1 -- terminal voltage transducer, reec_con(,5) is trv
        if (g.reec.reec_con(kk,5) < lbnd)
            g.reec.reec_con(kk,5) = 5*lbnd;
            estr = '\nsvm_mgen: increasing reec voltage transducer time ';
            estr = [estr, 'constant trv to be nonzero at reec index %0.0f.'];
            warning(sprintf(estr,kk));
        end

        % reec2 -- active power transducer is present if pfflag == 1
        if (g.reec.reec_con(kk,32) == 1)
            state_reec(kk) = state_reec(kk) + 1;

            if (g.reec.reec_con(kk,14) < lbnd)
                g.reec.reec_con(kk,14) = 5*lbnd;
                estr = '\nsvm_mgen: increasing reec active power filter time ';
                estr = [estr, 'constant tp to be nonzero at reec index %0.0f.'];
                warning(sprintf(estr,kk));
            end
        end

        % reec3 -- reactive power transducer, reec_con(,15) is tq
        if (g.reec.reec_con(kk,15) < lbnd)
            g.reec.reec_con(kk,15) = 5*lbnd;
            estr = '\nsvm_mgen: increasing reec reactive power transducer time ';
            estr = [estr, 'constant tq to be nonzero at reec index %0.0f.'];
            warning(sprintf(estr,kk));
        end

        % reec4 -- first stage PI loop is present if vflag == 1 and qflag == 1
        if ((g.reec.reec_con(kk,33) == 1) && (g.reec.reec_con(kk,34) == 1))
            if (g.reec.reec_con(kk,22) > 0)
                state_reec(kk) = state_reec(kk) + 1;
            end
        end

        % reec5 -- second stage PI loop is present if qflag == 1
        if ((g.reec.reec_con(kk,34) == 1) && (g.reec.reec_con(kk,24) > 0))
            state_reec(kk) = state_reec(kk) + 1;
        end

        % reec6 -- reactive current filter is present if qflag == 0
        if ((g.reec.reec_con(kk,34) == 0) && (g.reec.reec_con(kk,25) >= lbnd))
            state_reec(kk) = state_reec(kk) + 1;
        end

        % reec7 -- active power order filter, reec_con(,30) is tpord
        if (g.reec.reec_con(kk,30) < lbnd)
            g.reec.reec_con(kk,30) = 5*lbnd;
            estr = '\nsvm_mgen: increasing reec active power order time ';
            estr = [estr, 'constant tpord to be nonzero at reec index %0.0f.'];
            warning(sprintf(estr,kk));
        end

        % reec8 -- voltage compensation filter, reec_con(,38) is tr1
        if (g.reec.reec_con(kk,38) >= lbnd)
            state_reec(kk) = state_reec(kk) + 1;
        end

        % reec9 -- active current command filter, reec_con(,43) is tipc
        if (g.reec.reec_con(kk,43) < lbnd)
            g.reec.reec_con(kk,43) = 5*lbnd;
            estr = '\nsvm_mgen: increasing reec active current order time ';
            estr = [estr, 'constant to be nonzero at reec index %0.0f.'];
            warning(sprintf(estr,kk));
        end

        % reec10 -- reactive current command filter, reec_con(,44) is tiqc
        if (g.reec.reec_con(kk,44) < lbnd)
            g.reec.reec_con(kk,44) = 5*lbnd;
            estr = '\nsvm_mgen: increasing reec reactive current order time ';
            estr = [estr, 'constant to be nonzero at reec index %0.0f.'];
            warning(sprintf(estr,kk));
        end
    end

    n_reec_states = sum(state_reec);
    state(n_reec1+1:n_reec1+g.reec.n_reec) = state_reec;
end

% gfma
n_gfma_states = 0;
n_gfma1 = n_reec1 + g.reec.n_reec;
if (g.gfma.n_gfma ~= 0)
    state_gfma = zeros(g.gfma.n_gfma,1);
    for kk = 1:g.gfma.n_gfma
        state_gfma(kk) = 5;                                 % non-bypassable states

        % note: overload mitigation states not included in linearization;
        %       this includes gfma2, gfma3, gfma6, gfma7

        % gfma1 -- commanded voltage angle integrator is not bypassable

        % gfma4 -- commanded voltage magnitude time constant
        if (g.gfma.gfma_con(kk,14) < lbnd)
            g.gfma.gfma_con(kk,14) = 5*lbnd;
            estr = '\nsvm_mgen: increasing gfma voltage magnitude time ';
            estr = [estr, 'constant to be nonzero at gfma index %0.0f.'];
            warning(sprintf(estr,kk));
        end

        % gfma5 -- voltage regulation integral state present if vflag == 1
        if ((g.gfma.gfma_con(kk,25) == 1) && (g.gfma.gfma_con(kk,13) > 0))
            state_gfma(kk) = state_gfma(kk) + 1;
        end

        % gfma8 -- inverter active power transducer
        if (g.gfma.gfma_con(kk,21) < lbnd)
            g.gfma.gfma_con(kk,21) = 5*lbnd;
            estr = '\nsvm_mgen: increasing gfma active power transducer time ';
            estr = [estr, 'constant to be nonzero at gfma index %0.0f.'];
            warning(sprintf(estr,kk));
        end

        % gfma9 -- inverter reactive power transducer
        if (g.gfma.gfma_con(kk,22) < lbnd)
            g.gfma.gfma_con(kk,22) = 5*lbnd;
            estr = '\nsvm_mgen: increasing gfma reactive power transducer time ';
            estr = [estr, 'constant to be nonzero at gfma index %0.0f.'];
            warning(sprintf(estr,kk));
        end

        % gfma10 -- inverter voltage transducer
        if (g.gfma.gfma_con(kk,23) < lbnd)
            g.gfma.gfma_con(kk,23) = 5*lbnd;
            estr = '\nsvm_mgen: increasing gfma voltage mag. transducer time ';
            estr = [estr, 'constant to be nonzero at gfma index %0.0f.'];
            warning(sprintf(estr,kk));
        end
    end

    n_gfma_states = sum(state_gfma);
    state(n_gfma1+1:n_gfma1+g.gfma.n_gfma) = state_gfma;
end

% eof
