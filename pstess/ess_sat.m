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

% checking for current limit violations
mask = (abs(i_inj) > g.ess.ess_pot(:,3));
if any(mask)
    ip_sign = sign(ip_inj);                       % storing injection signs
    iq_sign = sign(iq_inj);

    ip_inj(mask) = abs(ip_inj(mask));             % one-quadrant perspective
    iq_inj(mask) = abs(iq_inj(mask));

    ip_lim = g.ess.ess_pot(:,3);                  % initialize the limits
    iq_lim = g.ess.ess_pot(:,3);

    p_mask = mask & (g.ess.ess_con(:,9) == 1);    % active priority
    q_mask = mask & (g.ess.ess_con(:,9) == 2);    % reactive priority
    pq_mask = mask & (g.ess.ess_con(:,9) == 3);   % proportional scaling

    % active priority
    ip_inj(p_mask) = min(ip_inj(p_mask),g.ess.ess_pot(p_mask,3));
    iq_lim(p_mask) = sqrt(g.ess.ess_pot(p_mask,3).^2 - ip_inj(p_mask).^2);
    iq_inj(p_mask) = min(iq_inj(p_mask),iq_lim(p_mask));

    % reactive priority
    iq_inj(q_mask) = min(iq_inj(q_mask),g.ess.ess_pot(q_mask,3));
    ip_lim(q_mask) = sqrt(g.ess.ess_pot(q_mask,3).^2 - iq_inj(q_mask).^2);
    ip_inj(q_mask) = min(ip_inj(q_mask),ip_lim(q_mask));

    % proportional scaling
    pq_scale = g.ess.ess_pot(pq_mask,3)./abs(i_inj(pq_mask));
    ip_inj(pq_mask) = pq_scale.*ip_inj(pq_mask);
    iq_inj(pq_mask) = pq_scale.*iq_inj(pq_mask);

    % accounting for positive and negative injections
    ip_inj(mask) = ip_sign(mask).*abs(ip_inj(mask));
    iq_inj(mask) = iq_sign(mask).*abs(iq_inj(mask));
end

%-----------------------------------------------------------------%
% rate limiting

if (flag == 2 || flag == 3)  % dynamics calculation
    % determining the appropriate time constant
    if (nargin ~= 6)
        error('ess_sat: incorrect number of arguments passed to function.');
    else
        h_sol = varargin{1};
    end

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

mask = (v_term < v_break);                        % linear roll-off
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
    iq_pf_ub(mask) = ...
        sqrt((ip_inj(mask)./g.ess.ess_con(mask,10)).^2 - ip_inj(mask).^2);
end

iq_pf_lb = -iq_pf_ub;                             % symmetric limits

iq_inj = max(iq_inj,iq_pf_lb);
iq_inj = min(iq_inj,iq_pf_ub);

end  % function end

% eof
