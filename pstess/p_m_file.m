% This m-file forms a sparse permutation matrix p_mat
% which converts the vector of rates of change
% of state perturbations (d_vector) into a column of
% the state matrix a_mat
%
% called by svm_mgen

%-----------------------------------------------------------------------------%
% Version history
%
% Modified: added lsc
% Date:     December 2020
% Author:   R. Elliott
%
% Modified: added ess
% Date:     2019
% Author:   R. Elliott
%
% Modified: added pwrmod
% Date:     2015
% Author:   D. Trudnowski
%
% Modified: tcsc model added
% Date:     December 1998
% Author:   Graham Rogers
%
% Modified: hydro turbine added
% Date:     June 1998
% Author:   Graham Rogers
%
% Modified: lead lag added to SVC
% Date:     April 1998
% Author:   Graham Rogers
%
% Modified: Induction Generator added, Load modulation added
% Date:     August 1997
% Author:   Graham Rogers
%
% Modified: HVDC added
% Date:     April 1997
% Author:   Graham Rogers
%
% Modified: November 1996
%           svc added
% Author:   Graham Rogers
% Date:     September 1996
%-----------------------------------------------------------------------------%

k_row = 1;
exc_count = 0;
pss_count = 0;
dpw_count = 0;
tg_count = 0;

for k = 1:g.mac.n_mac
    if (state(k) ~= 0)
        if (k ~= 1)
            k_row = 1 + sum(state(1:k-1));
        end

        k_col = k;

        k_ivm = [];
        if ~isempty(g.mac.mac_ivm_idx)
            k_ivm = find(g.mac.mac_ivm_idx == k);
        end

        if isempty(k_ivm)
            % all synchronous generators (mac_ang, mac_spd)
            for kgs = 1:nss.mac_em
                p_mat(k_row+kgs-1,k_col+(kgs-1)*g.mac.n_mac) = 1;
            end

            % subtransient and transient models (eqprime)
            if (gen_state(k) > nss.mac_em)
                p_mat(k_row+2,k_col+nss.mac_em*g.mac.n_mac) = 1;
            end

            % transient models (edprime)
            if (gen_state(k) == nss.mac_tra)
                p_mat(k_row+3,k_col+nss.mac_tra*g.mac.n_mac) = 1;
            end

            % subtransient models (psikd, edprime, psikq)
            if (gen_state(k) == nss.mac_sub)
                for kgs = nss.mac_tra:nss.mac_sub
                    p_mat(k_row+kgs-1,k_col+(kgs-1)*g.mac.n_mac) = 1;
                end
            end
        else
            % ivm generators (mac_ang)
            p_mat(k_row,k_col) = 1;

            % ivm generators (eqprime)
            p_mat(k_row+1,k_col+nss.mac_em*g.mac.n_mac) = 1;
        end

        % exciters
        k_row = k_row + gen_state(k) - 1;
        k_col = nss.mac_max*g.mac.n_mac;
        if (mac_exc ~= 0)
            k_exc = find(mac_exc == k);
            if ~isempty(k_exc)
                exc_count = exc_count + 1;

                k_cex = k_col + k_exc;
                if (TR_state(k) ~= 0)
                    k_row = k_row + 1;
                    p_mat(k_row,k_cex) = 1;
                end

                k_cex = k_cex + g.exc.n_exc;
                if (TB_state(k) ~= 0)
                    k_row = k_row + 1;
                    p_mat(k_row,k_cex) = 1;
                end

                k_cex = k_cex + g.exc.n_exc;
                if (TA_state(k) ~= 0)
                    k_row = k_row + 1;
                    p_mat(k_row,k_cex) = 1;
                end

                k_cex = k_cex + g.exc.n_exc;
                if (Efd_state(k) ~= 0)
                    k_row = k_row + 1;
                    p_mat(k_row,k_cex) = 1;
                end

                k_cex = k_cex + g.exc.n_exc;
                if (R_f_state(k) ~= 0)
                    k_row = k_row + 1;
                    p_mat(k_row,k_cex) = 1;
                end
            end
        end

        % pss
        k_col = k_col + nss.exc_max*g.exc.n_exc;
        if (mac_pss ~= 0)
            k_pss = find(mac_pss == k);
            if ~isempty(k_pss)
                pss_count = pss_count + 1;

                k_cpss = k_col + k_pss;
                if (pss1_state(k) ~= 0)
                    k_row = k_row + 1;
                    p_mat(k_row,k_cpss) = 1;
                end

                k_cpss = k_cpss + g.pss.n_pss;
                if (pss2_state(k) ~= 0)
                    k_row = k_row + 1;
                    p_mat(k_row,k_cpss) = 1;
                end

                k_cpss = k_cpss + g.pss.n_pss;
                if (pss3_state(k) ~= 0)
                    k_row = k_row + 1;
                    p_mat(k_row,k_cpss) = 1;
                end
            end
        end

        % dpw filter
        k_col = k_col + nss.pss_max*g.pss.n_pss;
        if (mac_dpw ~= 0)
            k_dpw = find(mac_dpw == k);
            if ~isempty(k_dpw)
                dpw_count = dpw_count + 1;

                k_cdpw = k_col + k_dpw;
                k_row = k_row + 1;
                p_mat(k_row,k_cdpw) = 1;

                k_row1 = k_row;
                for kj = 1:nss.dpw_max-1
                    k_cdpw = k_cdpw + g.dpw.n_dpw;
                    k_row = k_row1 + kj;
                    p_mat(k_row,k_cdpw) = 1;
                end
            end
        end

        % thermal governors
        k_col = k_col + nss.dpw_max*g.dpw.n_dpw;
        if (mac_tg ~= 0)
            k_tg = find(mac_tg == k);
            if ~isempty(k_tg)
                tg_count = tg_count + 1;
                k_ctg = k_col + g.tg.tg_idx(k_tg);
                if (tg_state(k) ~= 0)
                    k_row = k_row + 1;
                    p_mat(k_row,k_ctg) = 1;

                    k_row1 = k_row;
                    for kj = 1:nss.tg-1
                        k_ctg = k_ctg + g.tg.n_tg_tot;
                        k_row = k_row1 + kj;
                        p_mat(k_row,k_ctg) = 1;
                    end
                end
            end
        end

        % hydro governors
        if (mac_tgh ~= 0)
            k_tgh = find(mac_tgh == k);
            if ~isempty(k_tgh)
                tg_count = tg_count + 1;
                k_ctg = k_col + g.tg.tgh_idx(k_tgh);
                if (tgh_state(k) ~= 0)
                    k_row = k_row + 1;
                    p_mat(k_row,k_ctg) = 1;

                    k_row1 = k_row;
                    for kj = 1:nss.tgh-1
                        k_ctg = k_ctg + g.tg.n_tg_tot;
                        k_row = k_row1 + kj;
                        p_mat(k_row,k_ctg) = 1;
                    end
                end
            end
        end
    end
end

k_col1 = nss.mac_max*g.mac.n_mac + nss.exc_max*g.exc.n_exc ...
         + nss.pss_max*g.pss.n_pss + nss.dpw_max*g.dpw.n_dpw ...
         + nss.tg_max*g.tg.n_tg_tot;

% induction motors
k_col = k_col1;
if (g.ind.n_mot ~= 0)
    for k = 1:g.ind.n_mot
        k_col = k_col + k;
        k_row = k_row + 1;
        p_mat(k_row,k_col) = 1;

        k_row1 = k_row;
        for kj = 1:nss.mot_max-1
            k_col = k_col + g.ind.n_mot;
            k_row = k_row1 + kj;
            p_mat(k_row,k_col) = 1;
        end

        k_col = k_col1;
    end
end

% induction generators
k_col1 = k_col1 + nss.mot_max*g.ind.n_mot;
k_col = k_col1;
if (g.igen.n_ig ~= 0)
    for k = 1:g.igen.n_ig
        k_col = k_col + k;
        k_row = k_row + 1;
        p_mat(k_row,k_col) = 1;

        k_row1 = k_row;
        for kj = 1:nss.ig_max-1
            k_col = k_col + g.igen.n_ig;
            k_row = k_row1 + kj;
            p_mat(k_row,k_col) = 1;
        end

        k_col = k_col1;
    end
end

% svc
k_col1 = k_col1 + nss.ig_max*g.igen.n_ig;
k_col = k_col1;
if (g.svc.n_svc ~= 0)
    for k = 1:g.svc.n_svc
        k_col = k_col + k;
        k_row = k_row + 1;
        p_mat(k_row,k_col) = 1;

        if ~isempty(g.svc.svcll_idx)
            kcon = find(k == g.svc.svcll_idx);
            if ~isempty(kcon)
                k_col = k_col + g.svc.n_svc;
                k_row = k_row + 1;
                p_mat(k_row,k_col) = 1;
            end
        end

        k_col = k_col1;
    end
end

% tcsc
k_col1 = k_col1 + nss.svc_max*g.svc.n_svc;
k_col = k_col1;
if (g.tcsc.n_tcsc ~= 0)
    for k = 1:g.tcsc.n_tcsc
        k_col = k_col + k;
        k_row = k_row + 1;
        p_mat(k_row,k_col) = 1;

        k_col = k_col1;
    end
end

% lmod
k_col1 = k_col1 + nss.tcsc_max*g.tcsc.n_tcsc;
k_col = k_col1;
if (g.lmod.n_lmod ~= 0)
    for k = 1:g.lmod.n_lmod
        k_col = k_col + k;
        k_row = k_row + 1;
        p_mat(k_row,k_col) = 1;

        k_col = k_col1;
    end
end

% rlmod
k_col1 = k_col1 + nss.lmod_max*g.lmod.n_lmod;
k_col = k_col1;
if (g.rlmod.n_rlmod ~= 0)
    for k = 1:g.rlmod.n_rlmod
        k_col = k_col + k;
        k_row = k_row + 1;
        p_mat(k_row,k_col) = 1;

        k_col = k_col1;
    end
end

% pwrmod_p
k_col1 = k_col1 + nss.rlmod_max*g.rlmod.n_rlmod;
k_col = k_col1;
if (g.pwr.n_pwrmod ~= 0)
    for k = 1:g.pwr.n_pwrmod
        k_col = k_col + k;
        k_row = k_row + 1;
        p_mat(k_row,k_col) = 1;

        k_col = k_col1;
    end
end

% pwrmod_q
k_col1 = k_col1 + nss.pwrmod_max*g.pwr.n_pwrmod;
k_col = k_col1;
if (g.pwr.n_pwrmod ~= 0)
    for k = 1:g.pwr.n_pwrmod
        k_col = k_col + k;
        k_row = k_row + 1;
        p_mat(k_row,k_col) = 1;

        k_col = k_col1;
    end
end

% hvdc
k_col1 = k_col1 + nss.pwrmod_max*g.pwr.n_pwrmod;
k_col = k_col1;
if (g.dc.n_conv ~= 0)
    for k = 1:g.dc.n_dcl
        % converter controls
        k_col = k_col + k;
        k_row = k_row + 1;
        p_mat(k_row,k_col) = 1;

        k_col = k_col + g.dc.n_dcl;
        k_row = k_row + 1;
        p_mat(k_row,k_col) = 1;

        % hvdc line
        k_col = k_col + g.dc.n_dcl;
        k_row = k_row + 1;
        p_mat(k_row,k_col) = 1;

        if (g.dc.l_cap ~= 0)
            if ~isempty(g.dc.cap_idx(k))
                % capacitor in this line model
                k_row1 = k_row;
                for kj = 1:nss.l_cap
                    k_col = k_col + g.dc.n_dcl;
                    k_row = k_row1 + kj;
                    p_mat(k_row,k_col) = 1;
                end
            end
        end

        k_col = k_col1;
    end
end

% ess
k_col1 = k_col1 + nss.dcl_max*g.dc.n_dcl;
k_col = k_col1;
if (g.ess.n_ess ~= 0)
    for k = 1:g.ess.n_ess
        % transducer for local voltage magnitude (non-bypassable)
        k_col = k_col + k;
        k_row = k_row + 1;
        p_mat(k_row,k_col) = 1;

        % local voltage magnitude Pade approximation
        k_col = k_col + g.ess.n_ess;
        if ((g.ess.ess_con(k,4) ~= 0) && (g.ess.ess_con(k,5)/2 >= lbnd))
            k_row = k_row + 1;
            p_mat(k_row,k_col) = 1;
        end

        % active current converter interface (non-bypassable)
        k_col = k_col + g.ess.n_ess;
        k_row = k_row + 1;
        p_mat(k_row,k_col) = 1;

        % reactive current converter interface (non-bypassable)
        k_col = k_col + g.ess.n_ess;
        k_row = k_row + 1;
        p_mat(k_row,k_col) = 1;

      % % state of charge integrator (no effect on linearization)
      % k_col = k_col + g.ess.n_ess;
      % k_row = k_row + 1;
      % p_mat(k_row,k_col) = 1;

        k_col = k_col1;
    end
end

% lsc
k_col1 = k_col1 + nss.ess_max*g.ess.n_ess;
k_col = k_col1;
if (g.lsc.n_lsc ~= 0)
    for k = 1:g.lsc.n_lsc
        % transducer for remote angle signal (non-bypassable)
        k_col = k_col + k;
        k_row = k_row + 1;
        p_mat(k_row,k_col) = 1;

        % transducer for local angle signal (non-bypassable)
        k_col = k_col + g.lsc.n_lsc;
        k_row = k_row + 1;
        p_mat(k_row,k_col) = 1;

        % remote angle Pade approximation
        k_col = k_col + g.lsc.n_lsc;
        if ((g.lsc.lsc_con(k,4) ~= 0) && (g.lsc.lsc_con(k,5)/2 >= lbnd))
            k_row = k_row + 1;
            p_mat(k_row,k_col) = 1;
        end

        % local angle Pade approximation
        k_col = k_col + g.lsc.n_lsc;
        if ((g.lsc.lsc_con(k,4) ~= 0) && (g.lsc.lsc_con(k,6)/2 >= lbnd))
            k_row = k_row + 1;
            p_mat(k_row,k_col) = 1;
        end

        % remote angle setpoint tracking filter
        k_col = k_col + g.lsc.n_lsc;
        if (g.lsc.lsc_con(k,8) >= lbnd)
            k_row = k_row + 1;
            p_mat(k_row,k_col) = 1;
        end

        % local angle setpoint tracking filter
        k_col = k_col + g.lsc.n_lsc;
        if (g.lsc.lsc_con(k,10) >= lbnd)
            k_row = k_row + 1;
            p_mat(k_row,k_col) = 1;
        end

        % local LTV highpass filter state 1 (non-bypassable)
        k_col = k_col + g.lsc.n_lsc;
        k_row = k_row + 1;
        p_mat(k_row,k_col) = 1;

        % local LTV highpass filter state 2 (non-bypassable)
        k_col = k_col + g.lsc.n_lsc;
        k_row = k_row + 1;
        p_mat(k_row,k_col) = 1;

        % local LTV lead-lag stage 1
        k_col = k_col + g.lsc.n_lsc;
        if (g.lsc.lsc_con(k,17) >= lbnd)
            k_row = k_row + 1;
            p_mat(k_row,k_col) = 1;
        end

        % local LTV lead-lag stage 2
        k_col = k_col + g.lsc.n_lsc;
        if (g.lsc.lsc_con(k,19) >= lbnd)
            k_row = k_row + 1;
            p_mat(k_row,k_col) = 1;
        end

        % center LTI highpass filter state 1 (non-bypassable)
        k_col = k_col + g.lsc.n_lsc;
        k_row = k_row + 1;
        p_mat(k_row,k_col) = 1;

        % center LTI highpass filter state 2 (non-bypassable)
        k_col = k_col + g.lsc.n_lsc;
        k_row = k_row + 1;
        p_mat(k_row,k_col) = 1;

        % center LTI lead-lag stage 1
        k_col = k_col + g.lsc.n_lsc;
        if (g.lsc.lsc_con(k,28) >= lbnd)
            k_row = k_row + 1;
            p_mat(k_row,k_col) = 1;
        end

        % center LTI lead-lag stage 2
        k_col = k_col + g.lsc.n_lsc;
        if (g.lsc.lsc_con(k,30) >= lbnd)
            k_row = k_row + 1;
            p_mat(k_row,k_col) = 1;
        end

        % lowpass filter (non-bypassable)
        k_col = k_col + g.lsc.n_lsc;
        k_row = k_row + 1;
        p_mat(k_row,k_col) = 1;

        k_col = k_col1;
    end
end

% reec
k_col1 = k_col1 + nss.lsc_max*g.lsc.n_lsc;
k_col = k_col1;
if (g.reec.n_reec ~= 0)
    for k = 1:g.reec.n_reec
        % reec1 -- terminal voltage transducer
        k_col = k_col + k;
        k_row = k_row + 1;
        p_mat(k_row,k_col) = 1;

        % reec2 -- active power transducer
        k_col = k_col + g.reec.n_reec;
        if (g.reec.reec_con(k,32) == 1)
            k_row = k_row + 1;
            p_mat(k_row,k_col) = 1;
        end

        % reec3 -- reactive power transducer
        k_col = k_col + g.reec.n_reec;
        k_row = k_row + 1;
        p_mat(k_row,k_col) = 1;

        % reec4 -- first stage PI loop
        k_col = k_col + g.reec.n_reec;
        if ((g.reec.reec_con(k,33) == 1) && (g.reec.reec_con(k,34) == 1))
            if (g.reec.reec_con(k,22) > 0)
                k_row = k_row + 1;
                p_mat(k_row,k_col) = 1;
            end
        end

        % reec5 -- second stage PI loop
        k_col = k_col + g.reec.n_reec;
        if ((g.reec.reec_con(k,34) == 1) && (g.reec.reec_con(k,24) > 0))
            k_row = k_row + 1;
            p_mat(k_row,k_col) = 1;
        end

        % reec6 -- reactive current order filter
        k_col = k_col + g.reec.n_reec;
        if ((g.reec.reec_con(k,34) == 0) && (g.reec.reec_con(k,25) >= lbnd))
            k_row = k_row + 1;
            p_mat(k_row,k_col) = 1;
        end

        % reec7 -- active power order filter
        k_col = k_col + g.reec.n_reec;
        k_row = k_row + 1;
        p_mat(k_row,k_col) = 1;

        % reec8 -- voltage magnitude compensation filter
        k_col = k_col + g.reec.n_reec;
        if (g.reec.reec_con(k,38) >= lbnd)
            k_row = k_row + 1;
            p_mat(k_row,k_col) = 1;
        end

        % reec9 -- active current command filter
        k_col = k_col + g.reec.n_reec;
        k_row = k_row + 1;
        p_mat(k_row,k_col) = 1;

        % reec10 -- reactive current command filter
        k_col = k_col + g.reec.n_reec;
        k_row = k_row + 1;
        p_mat(k_row,k_col) = 1;

        k_col = k_col1;
    end
end

% gfma
k_col1 = k_col1 + nss.reec_max*g.reec.n_reec;
k_col = k_col1;
if (g.gfma.n_gfma ~= 0)
    for k = 1:g.gfma.n_gfma
        % gfma1 -- commanded voltage angle integrator
        k_col = k_col + k;
        k_row = k_row + 1;
        p_mat(k_row,k_col) = 1;

        % note: overload mitigation states are not included in linearization;
        %       this includes gfma2, gfma3, gfma6, gfma7

        % gfma4 -- commanded voltage magnitude (first-order)
        k_col = k_col + g.gfma.n_gfma;
        k_row = k_row + 1;
        p_mat(k_row,k_col) = 1;

        % gfma5 -- voltage regulation integral control state
        k_col = k_col + g.gfma.n_gfma;
        if ((g.gfma.gfma_con(k,25) == 1) && (g.gfma.gfma_con(k,13) > 0))
            k_row = k_row + 1;
            p_mat(k_row,k_col) = 1;
        end

        % gfma8 -- active power transducer
        k_col = k_col + g.gfma.n_gfma;
        k_row = k_row + 1;
        p_mat(k_row,k_col) = 1;

        % gfma9 -- reactive power transducer
        k_col = k_col + g.gfma.n_gfma;
        k_row = k_row + 1;
        p_mat(k_row,k_col) = 1;

        % gfma10 -- voltage transducer
        k_col = k_col + g.gfma.n_gfma;
        k_row = k_row + 1;
        p_mat(k_row,k_col) = 1;

        k_col = k_col1;
    end
end

% k_col1 = k_col1 + nss.gfma_max*g.gfma.n_gfma;  % for future expansion

% eof
