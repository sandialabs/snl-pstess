function [ip_inj,iq_inj] = ess_sat(k,p_inj,q_inj,v_term,flag,varargin)
% Syntax: [ip_inj,iq_inj] = ess_sat(k,p_inj,q_inj,v_term,flag) OR
%                           ess_sat(k,p_inj,q_inj,v_term,flag,h_sol)
%
% Purpose: saturation routine for inverter-based resources
%          and energy storage systems,
%          NOTE - this is an auxiliary function, not a dynamic model
%
% Input:   p_inj - active power command (pu on syst. base)
%          q_inj - reactive power command (pu on syst. base)
%          v_term - terminal voltage magnitude (pu)
%          flag - 2 - generator dynamics computation
%                 3 - generator dynamics computation for linearization
%                 4 - iterative network interface computation (for nc_load)
%
% Output:  ip_inj - post-saturation active current command (pu on syst. base)
%          iq_inj - post-saturation reactive current command (pu on syst. base)
%
% Note:    complex power modulation priority mode is specified as either
%          1 = active; 2 = reactive; 3 = proportional

%-----------------------------------------------------------------------------%
% Version history
%
% Version:  1.1
% Date:     February 2024
% Author:   Ryan Elliott
% Purpose:  Updated to handle REEC_D protection logic
%
% Version:  1.0
% Date:     November 2020
% Author:   Ryan Elliott
% Purpose:  Inverter-based resource saturation with priority modes
%-----------------------------------------------------------------------------%

global g;  % declaring struct of global variables

lbnd = 1e-3;  % lower bound to prevent division by zero

% pre-saturation complex current
ip_inj = p_inj./max(v_term,lbnd);
iq_inj = q_inj./max(v_term,lbnd);
i_inj = ip_inj + 1j*iq_inj;

%-----------------------------------------------------------------%
% interpolating the vdl curves

vdlp = zeros(g.ess.n_ess,1);
vdlq = zeros(g.ess.n_ess,1);
for ii = 1:g.ess.n_ess
    vdl = g.ess.vdl{ii};

    % interpolating the active current curves
    mp = diff(vdl.ip)./diff(vdl.vp);

    if (v_term(ii) <= vdl.vp(1))
        vdlp(ii) = vdl.ip(1);
    elseif (v_term(ii) >= vdl.vp(end))
        vdlp(ii) = vdl.ip(end);
    else
        % binary search, leftmost element procedure
        jj = 0;
        ub = length(vdl.vp);
        while (jj < ub)
            mb = floor((jj + ub)/2);
            if (v_term(ii) > vdl.vp(mb+1))
                jj = mb + 1;
            else
                ub = mb;
            end
        end

        vdlp(ii) = mp(jj)*(v_term(ii) - vdl.vp(jj)) + vdl.ip(jj);
    end

    %-------------------------------------------------------------%
    % interpolating the reactive current curves

    mq = diff(vdl.iq)./diff(vdl.vq);

    if (v_term(ii) <= vdl.vq(1))
        vdlq(ii) = vdl.iq(1);
    elseif (v_term(ii) >= vdl.vq(end))
        vdlq(ii) = vdl.iq(end);
    else
        % binary search, leftmost element procedure
        jj = 0;
        ub = length(vdl.vq);
        while (jj < ub)
            mb = floor((jj + ub)/2);
            if (v_term(ii) > vdl.vq(mb+1))
                jj = mb + 1;
            else
                ub = mb;
            end
        end

        vdlq(ii) = mq(jj)*(v_term(ii) - vdl.vq(jj)) + vdl.iq(jj);
    end
end

%-----------------------------------------------------------------%
% current limiting

ip_max = min(vdlp,g.ess.ess_pot(:,3));
iq_max = min(vdlq,g.ess.ess_pot(:,3));

mask = (vdlq < 0);
if any(mask)
    % binding the iq command when vdlq is negative
    iq_max(mask) = max(vdlq(mask),-g.ess.ess_pot(mask,3));
    iq_inj(mask) = iq_max(mask);
end

p_mask = (g.ess.ess_con(:,9) == 1) & (iq_max >= 0);   % active priority
q_mask = (g.ess.ess_con(:,9) == 2) | (iq_max < 0);    % reactive priority
pq_mask = (g.ess.ess_con(:,9) == 3) & (iq_max >= 0);  % proportional scaling

% scaling the current command in proportional mode
mask = pq_mask & (abs(i_inj) > g.ess.ess_pot(:,3));
if any(mask)
    pq_scale = g.ess.ess_pot(mask,3)./abs(i_inj(mask));
    ip_inj(mask) = pq_scale.*ip_inj(mask);
    iq_inj(mask) = pq_scale.*iq_inj(mask);
end

ip_cap = sqrt(g.ess.ess_pot(:,3).^2 - iq_inj.^2);
iq_cap = sqrt(g.ess.ess_pot(:,3).^2 - ip_inj.^2);

ip_max(q_mask) = min(ip_max(q_mask),real(ip_cap(q_mask)));
iq_max(p_mask) = min(iq_max(p_mask),real(iq_cap(p_mask)));

ip_min = -g.ess.ess_con(:,22).*ip_max;
iq_min = -iq_max;
iq_min(iq_max < 0) = iq_max(iq_max < 0);

ip_inj = min(ip_inj,ip_max);
ip_inj = max(ip_inj,ip_min);

iq_inj = min(iq_inj,iq_max);
iq_inj = max(iq_inj,iq_min);

%-----------------------------------------------------------------%
% rate limiting and reec protection logic

if (flag == 2 || flag == 3)  % dynamics calculation
    % determining the appropriate time constant
    if (nargin ~= 6)
        error('ess_sat: incorrect number of arguments passed to function.');
    else
        h_sol = varargin{1};
    end

    %-------------------------------------------------------------%
    % reec protection logic

    if (g.reec.n_reec ~= 0)
        % passing iq limits back to the reec controls
        g.reec.iqmin = iq_min;
        g.reec.iqmax = iq_max;

        %---------------------------------------------------------%
        % voltage dip logic

        % vt_filt > vdiph or vt_filt < vdipl
        mask = (g.reec.reec1(:,k) > g.reec.reec_con(:,3)) ...
               | (g.reec.reec1(:,k) < g.reec.reec_con(:,4));

        if any(mask)
            g.reec.vdip(mask) = true;
            g.reec.vdip_tick(mask) = k;
            g.reec.vdip_time(mask) = 0;
        end

        % recovery stage
        mask = g.reec.vdip ...
               & (g.reec.reec1(:,k) <= g.reec.reec_con(:,3)) ...
               & (g.reec.reec1(:,k) >= g.reec.reec_con(:,4));

        if any(mask)
            % r_mask -- check for recovery transition, store vdip_icmd
            r_mask = mask & (g.reec.vdip_time == 0);
            if any(r_mask)
                iq_hld = iq_inj;
                qfrz_mask = (g.reec.reec_con(:,12) < 0);
                iq_hld(qfrz_mask) = g.reec.reec_con(qfrz_mask, 11);
                g.reec.vdip_icmd(r_mask) = ip_inj(r_mask) + 1j*iq_hld(r_mask);
            end

            % q_mask -- hold the reactive current until thldq
            q_mask = mask & (g.reec.vdip_time <= abs(g.reec.reec_con(:,12)));
            iq_inj(q_mask) = imag(g.reec.vdip_icmd(q_mask));

            % p_mask -- hold the active current until thldp
            p_mask = mask & (g.reec.vdip_time <= g.reec.reec_con(:,13));
            ip_inj(p_mask) = real(g.reec.vdip_icmd(p_mask));

            % t_mask -- track the time since recovering
            t_mask = mask & (g.reec.vdip_tick < k);
            if any(t_mask)
                g.reec.vdip_tick(t_mask) = k;
                g.reec.vdip_time(t_mask) = g.reec.vdip_time(t_mask) + h_sol;
            end

            % s_mask -- check the reset condition
            s_mask = mask & (g.reec.vdip_time > max(abs(g.reec.reec_con(:,12)), ...
                                                    g.reec.reec_con(:,13)));

            if any(s_mask)
                g.reec.vdip(s_mask) = false;
                g.reec.vdip_tick(s_mask) = -1;
                g.reec.vdip_time(s_mask) = 0;
                g.reec.vdip_icmd(s_mask) = 0;
            end
        end

        %---------------------------------------------------------%
        % voltage blocking logic

        mask = (g.reec.reec1(:,k) > g.reec.reec_con(:,40)) ...
               | (g.reec.reec1(:,k) < g.reec.reec_con(:,41));

        if any(mask)
            ip_inj(mask) = 0;
            iq_inj(mask) = 0;

            g.reec.vblk(mask) = true;
            g.reec.vblk_tick(mask) = k;
            g.reec.vblk_time(mask) = 0;
        end

        % blocking timer
        mask = g.reec.vblk ...
               & (g.reec.reec1(:,k) <= g.reec.reec_con(:,40)) ...
               & (g.reec.reec1(:,k) >= g.reec.reec_con(:,41));

        if any(mask)
            % b_mask -- block the current command until tblkdelay
            b_mask = mask & (g.reec.vblk_time <= g.reec.reec_con(:,42));
            ip_inj(b_mask) = 0;
            iq_inj(b_mask) = 0;

            % t_mask -- track the time since recovering
            t_mask = mask & (g.reec.vblk_tick < k);
            if any(t_mask)
                g.reec.vblk_tick(t_mask) = k;
                g.reec.vblk_time(t_mask) = g.reec.vblk_time(t_mask) + h_sol;
            end

            % s_mask -- check the reset condition
            s_mask = mask & (g.reec.vblk_time > g.reec.reec_con(:,42));
            if any(s_mask)
                g.reec.vblk(s_mask) = false;
                g.reec.vblk_tick(s_mask) = -1;
                g.reec.vblk_time(s_mask) = 0;
            end
        end
    end

    %-------------------------------------------------------------%
    % rate limiting

    rate_time_con = g.ess.ess_con(:,14);
    rate_time_con(rate_time_con < lbnd) = h_sol;

    % active current
    ip_diff = ip_inj - g.ess.ess3(:,k);

    ip_diff_ub = g.ess.ess_con(:,15).*g.ess.ess_pot(:,1).*rate_time_con;
    ip_diff_lb = -g.ess.ess_con(:,15).*g.ess.ess_pot(:,1).*rate_time_con;

    ip_diff = max(ip_diff,ip_diff_lb);
    ip_diff = min(ip_diff,ip_diff_ub);

    ip_inj = ip_diff + g.ess.ess3(:,k);

    % reactive current
    iq_diff = iq_inj - g.ess.ess4(:,k);

    iq_diff_ub = g.ess.ess_con(:,16).*g.ess.ess_pot(:,1).*rate_time_con;
    iq_diff_lb = -g.ess.ess_con(:,16).*g.ess.ess_pot(:,1).*rate_time_con;

    iq_diff = max(iq_diff,iq_diff_lb);
    iq_diff = min(iq_diff,iq_diff_ub);

    iq_inj = iq_diff + g.ess.ess4(:,k);
end

%-----------------------------------------------------------------%
% LVPL logic (limits active current only)

% i_lvpl1 -- breakpoint current (pu on syst. base)
% v_zerox -- zero crossing volt (pu)
% v_break -- breakpoint voltage (pu)
i_lvpl1 = g.ess.ess_con(:,17).*g.ess.ess_pot(:,1);
v_zerox = g.ess.ess_con(:,18);
v_break = g.ess.ess_con(:,19);

a_slope = i_lvpl1./(v_break - v_zerox);
b_inter = -i_lvpl1.*v_zerox./(v_break - v_zerox);

ip_lvpl_ub = i_lvpl1;

mask = (v_term < v_break) & (i_lvpl1 > 0);        % linear roll-off
if any(mask)
    ip_lvpl_ub(mask) = a_slope(mask).*max(v_term(mask),0) + b_inter(mask);

    zero_mask = (v_term < v_zerox);               % zero after this point
    if any(zero_mask)
        ip_lvpl_ub(zero_mask) = 0.0;
    end

    ip_lvpl_lb = -ip_lvpl_ub;                     % symmetric limits

    ip_inj(mask) = max(ip_inj(mask),ip_lvpl_lb(mask));
    ip_inj(mask) = min(ip_inj(mask),ip_lvpl_ub(mask));
end

%-----------------------------------------------------------------%
% respecting converter power factor limits

iq_pf_ub = g.ess.ess_pot(:,3);                    % initializing the bound

mask = (g.ess.ess_con(:,10) ~= 0);
if any(mask)
    iq_pf_ub(mask) = sqrt((ip_inj(mask)./g.ess.ess_con(mask,10)).^2 ...
                          - ip_inj(mask).^2);
end

iq_pf_lb = -iq_pf_ub;                             % symmetric limits

iq_inj = max(iq_inj,iq_pf_lb);
iq_inj = min(iq_inj,iq_pf_ub);

end  % function end

% eof
