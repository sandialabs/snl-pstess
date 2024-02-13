function [x] = ivm_sud(i,k,bus,flag,x)
% Syntax: x = ivm_sud(i,k,bus,flag,x)
%
% Purpose: ivm user-defined control model
%
% Input:   i - ivm_sud controller number
%          k - integer time
%          flag - 0 - initialization
%                 1 - network interface computation
%                 2 - system dynamics computation
%          x - state object (struct)
%              x.s  - dynamic states
%              x.ds - state derivatives
%              x.v  - algebraic storage variables
%
% Output:  x - updated state object
%
% See Also: ivm

%-----------------------------------------------------------------------------%
% Note:    The ivm power base is inherited from the ivm to which the
%          controller is connected.
%
% ivmud_con matrix format
%    col    data                                              units
%     1     ivm_sud controller number                         integer
%     2     bus number                                        integer
%     3     P-f droop gain, mp                                pu
%     4     active power upper limit, Pmax                    pu
%     5     active power lower limit, Pmin                    pu
%     6     active overload proportional gain, kppmax         pu
%     7     active overload integral gain, kipmax             pu
%     8     active overload windup flag (1 = anti-windup)     binary
%     9     Q-V droop gain, mq                                pu
%    10     active power upper limit, Qmax                    pu
%    11     active power lower limit, Qmin                    pu
%    12     voltage controller proportional gain, kpv         pu
%    13     voltage controller integral gain, kiv             pu
%    14     commanded voltage magnitude time constant         s
%    15     reactive overload proportional gain, kpqmax       pu
%    16     reactive overload integral gain, kiqmax           pu
%    17     reactive overload windup flag (1 = anti-windup)   binary
%    18     voltage control loop upper limit, Emax            pu
%    19     voltage control loop upper limit, Emin            pu
%    20     (not used)                                        -
%    21     active power transducer time constant             s
%    22     reactive power transducer time constant           s
%    23     voltage magnitude transducer time constant        s
%    24     Q-V initialization flag, qvflag (note 1)          binary
%    25     voltage control flag, vflag (note 2)              binary
%
% note 1: when ivmud_con(,24) = 0, the reactive power setpoint is
%         also set to 0, as in the PNNL specification; however,
%         when ivmud_con(,24) = 1, the reactive power setpoint is
%         initialized based on the power flow solution and Vset
%         is calculated accordingly.
%
% note 2: when ivmud_con(,25) = 0, the PI section of the voltage
%         control loop is bypassed; however, when
%         ivmud_con(,25) = 1, the commanded voltage magnitude
%         is taken from the output of the PI section (after
%         limiting).
%
%-----------------------------------------------------------------------------%
% Version history
%
% Version: Example
% Author:  Ryan T. Elliott
% Date:    November 2023
%
% Note:    In this example, ivm_sud.m is a simplified version of the
%          native PSTess model gfma.m without current limiting.
%          See data/d2asbegp_ivm.m for information about the data and
%          parameters stored in g.ivm.ivmud_con.
%
%          NOTE: A given ivm generator should only have a single
%          controller providing modulation commands to it. For example,
%          an ivm generator may receive commands from ivm_sud
%          or a native PSTess model, such as gfma.m, but not both.
%
% Ref:     For additional information about the control strategy
%          implemented here, see the specification below:
%
%          "Model Specification of Droop-Controlled, Grid-Forming
%           Inverters (REGFM_A1)," PNNL-32278, September 2023.
%
%-----------------------------------------------------------------------------%
% ivm_sud states
%   var      description
%   x.s{1}   commanded voltage angle (integrator)
%   x.s{2}   max P overload mitigation control state
%   x.s{3}   min P overload mitigation control state
%   x.s{4}   commanded voltage magnitude (first-order)
%   x.s{5}   voltage regulation integral control state
%   x.s{6}   max Q overload mitigation control state
%   x.s{7}   min Q overload mitigation control state
%   x.s{8}   inverter active power transducer
%   x.s{9}   inverter reactive power transducer
%   x.s{10}  inverter voltage transducer
%
%-----------------------------------------------------------------------------%

global g;  % declaring struct of global variables

lbnd = 1e-3;  % lower bound to prevent division by zero

if (~isempty(g.ivm.ivmud_con) && i ~= 0)
    % non-vectorized computation
    error('\nivm_sud: all user-defined ivm calculations must be vectorized.');
elseif (~isempty(g.ivm.ivmud_con) && i == 0)
    if (flag == 0)  % initialization
        % busnum -- vector of bus numbers with ivm_sud controllers
        busnum = g.bus.bus_int(g.ivm.ivmud_con(:,2));

        % x.s{1} is assigned in s_simu as z_ivmud where
        % z_ivmud = zeros(g.ivm.n_ivmud,n_ktime);
        z_ivmud = x.s{1};

        for idx = 1:10
            x.s{idx} = z_ivmud;
            x.ds{idx} = z_ivmud;
        end

        % s{1}  -- commanded voltage angle (integrator)
        % s{2}  -- max P overload mitigation control state
        % s{3}  -- min P overload mitigation control state
        % s{4}  -- commanded voltage magnitude (first-order)
        % s{5}  -- voltage regulation integral control state
        % s{6}  -- max Q overload mitigation control state
        % s{7}  -- min Q overload mitigation control state
        % s{8}  -- inverter active power transducer
        % s{9}  -- inverter reactive power transducer
        % s{10} -- inverter voltage transducer

        x.s{1}(:,1) = g.mac.mac_ang(g.ivm.divmud_mac_idx,1);
        x.s{2}(:,1) = 0.0;
        x.s{3}(:,1) = 0.0;
        x.s{4}(:,1) = g.mac.eqprime(g.ivm.divmud_mac_idx,1);
        x.s{5}(:,1) = g.mac.eqprime(g.ivm.divmud_mac_idx,1);
        x.s{6}(:,1) = 0.0;
        x.s{7}(:,1) = 0.0;
        x.s{8}(:,1) = g.mac.pelect(g.ivm.divmud_mac_idx,1) ...  % P_inv
                      .*g.mac.mac_pot(g.ivm.divmud_mac_idx,1);
        x.s{9}(:,1) = g.mac.qelect(g.ivm.divmud_mac_idx,1) ...  % Q_inv
                      .*g.mac.mac_pot(g.ivm.divmud_mac_idx,1);
        x.s{10}(:,1) = g.mac.eterm(g.ivm.divmud_mac_idx,1);     % V_inv

        x.pset = diag(x.s{8}(:,1))*ones(g.ivm.n_ivmud,size(x.s{1},2));
        x.qset = zeros(g.ivm.n_ivmud,size(x.s{1},2));

        % x.vref is a single column, not a matrix
        x.vref = x.s{4}(:,1);                        % pi bypass
        mask = (g.ivm.ivmud_con(:,25) ~= 0);
        if any(mask)
            x.vref(mask) = x.s{10}(mask,1);          % pi control
        end

        x.vset = diag(x.vref + g.ivm.ivmud_con(:,9).*x.s{9}(:,1)) ...
                 *ones(g.ivm.n_ivmud,size(x.s{1},2));

        mask = (g.ivm.ivmud_con(:,24) == 0);
        if any(mask)
            x.qset(mask,:) = diag(x.s{9}(mask,1)) ...
                             *ones(size(x.qset(mask,:)));
            x.vset(mask,:) = diag(x.vref(mask)) ...
                             *ones(size(x.vset(mask,:)));
        end

        %-----------------------------------------------------------------%
        % initialization checks

        mask = (x.s{8}(:,1) > g.ivm.ivmud_con(:,4));
        if any(mask)
            estr = '\nivm_sud: controller initialized above Pmax at ';
            estr = [estr, 'bus %0.0f '];
            error(sprintf(estr,bus(busnum(mask),1)));
        end

        mask = (x.s{8}(:,1) < g.ivm.ivmud_con(:,5));
        if any(mask)
            estr = '\nivm_sud: controller initialized below Pmin at ';
            estr = [estr, 'bus %0.0f '];
            error(sprintf(estr,bus(busnum(mask),1)));
        end

        mask = (x.s{9}(:,1) > g.ivm.ivmud_con(:,10));
        if any(mask)
            estr = '\nivm_sud: controller initialized above Qmax at ';
            estr = [estr, 'bus %0.0f '];
            error(sprintf(estr,bus(busnum(mask),1)));
        end

        mask = (x.s{9}(:,1) < g.ivm.ivmud_con(:,11));
        if any(mask)
            estr = '\nivm_sud: controller initialized below Qmin at ';
            estr = [estr, 'bus %0.0f '];
            error(sprintf(estr,bus(busnum(mask),1)));
        end

        %-----------------------------------------------------------------%
        % parameter checks

        mask = (g.ivm.ivmud_con(:,4) < g.ivm.ivmud_con(:,5));
        if any(mask)
            estr = '\nivm_sud: invalid active power limits selected at ';
            estr = [estr, 'bus %0.0f '];
            error(sprintf(estr,bus(busnum(mask),1)));
        end

        mask = (g.ivm.ivmud_con(:,10) < g.ivm.ivmud_con(:,11));
        if any(mask)
            estr = '\nivm_sud: invalid reactive power limits selected at ';
            estr = [estr, 'bus %0.0f '];
            error(sprintf(estr,bus(busnum(mask),1)));
        end

        mask = (g.ivm.ivmud_con(:,18) < g.ivm.ivmud_con(:,19));
        if any(mask)
            estr = '\nivm_sud: invalid commanded voltage magnitude limits ';
            estr = [estr, 'selected at bus %0.0f '];
            error(sprintf(estr,bus(busnum(mask),1)));
        end

        %-----------------------------------------------------------------%
        % flag checks

        mask = ~ismember(g.ivm.ivmud_con(:,8),[0,1]);
        if any(mask)
            estr = '\nivm_sud: invalid active overload windup flag selected at ';
            estr = [estr, 'bus %0.0f '];
            error(sprintf(estr,bus(busnum(mask),1)));
        end

        mask = ~ismember(g.ivm.ivmud_con(:,17),[0,1]);
        if any(mask)
            estr = '\nivm_sud: invalid reactive overload windup flag selected at ';
            estr = [estr, 'bus %0.0f '];
            error(sprintf(estr,bus(busnum(mask),1)));
        end

        mask = ~ismember(g.ivm.ivmud_con(:,24),[0,1]);
        if any(mask)
            estr = '\nivm_sud: invalid Q-V initialization flag (qvflag) selected at';
            estr = [estr, 'bus %0.0f '];
            error(sprintf(estr,bus(busnum(mask),1)));
        end

        mask = ~ismember(g.ivm.ivmud_con(:,25),[0,1]);
        if any(mask)
            estr = '\nivm_sud: invalid voltage control flag (vflag) selected at ';
            estr = [estr, 'bus %0.0f '];
            error(sprintf(estr,bus(busnum(mask),1)));
        end
    end

    if (flag == 1)  % network interface
        g.mac.vex(g.ivm.divmud_mac_idx,k) = x.s{4}(:,k);
        g.mac.fldcur(g.ivm.divmud_mac_idx,k) = x.s{1}(:,k);
    end

    if (flag == 2 || flag == 3)  % dynamics calculation
        %---------------------------------------------------------------------%
        % transducers

        % active power
        x.ds{8}(:,k) = (g.mac.pelect(g.ivm.divmud_mac_idx,k) ...
                        .*g.mac.mac_pot(g.ivm.divmud_mac_idx,1) ...
                        - x.s{8}(:,k))./max(g.ivm.ivmud_con(:,21),lbnd);

        mask = (g.ivm.ivmud_con(:,21) < lbnd);
        if any(mask)
            x.ds{8}(mask,k) = 0.0;
            x.s{8}(mask,k) = g.mac.pelect(g.ivm.divmud_mac_idx(mask),k) ...
                             .*g.mac.mac_pot(g.ivm.divmud_mac_idx(mask),1);
        end

        % reactive power
        x.ds{9}(:,k) = (g.mac.qelect(g.ivm.divmud_mac_idx,k) ...
                        .*g.mac.mac_pot(g.ivm.divmud_mac_idx,1) ...
                        - x.s{9}(:,k))./max(g.ivm.ivmud_con(:,22),lbnd);

        mask = (g.ivm.ivmud_con(:,22) < lbnd);
        if any(mask)
            x.ds{9}(mask,k) = 0.0;
            x.s{9}(mask,k) = g.mac.qelect(g.ivm.divmud_mac_idx(mask),k) ...
                             .*g.mac.mac_pot(g.ivm.divmud_mac_idx(mask),1);
        end

        % voltage magnitude
        x.ds{10}(:,k) = ...
            (g.mac.eterm(g.ivm.divmud_mac_idx,k) - x.s{10}(:,k)) ...
            ./max(g.ivm.ivmud_con(:,23),lbnd);

        mask = (g.ivm.ivmud_con(:,23) < lbnd);
        if any(mask)
            x.ds{10}(mask,k) = 0.0;
            x.s{10}(mask,k) = g.mac.eterm(g.ivm.divmud_mac_idx(mask),k);
        end

        %---------------------------------------------------------------------%
        % P-f droop control loop

        P1 = g.ivm.ivmud_con(:,3).*(x.pset(:,k) - x.s{8}(:,k));

        x.ds{2}(:,k) = g.ivm.ivmud_con(:,7).*(g.ivm.ivmud_con(:,4) - x.s{8}(:,k));

        % anti-windup limits (upper bound)
        mask = (x.s{2}(:,k) > 0);
        if any(mask)
            x.s{2}(mask,k) = 0;

            anti_windup = mask & (g.ivm.ivmud_con(:,8) == 1) & (x.ds{2}(:,k) > 0);

            if any(anti_windup)
                x.ds{2}(anti_windup,k) = 0;
            end
        end

        P2 = g.ivm.ivmud_con(:,6) ...
             .*(g.ivm.ivmud_con(:,4) - x.s{8}(:,k)) + x.s{2}(:,k);

        P2 = min(P2,0);

        x.ds{3}(:,k) = g.ivm.ivmud_con(:,7).*(g.ivm.ivmud_con(:,5) - x.s{8}(:,k));

        % anti-windup limits (lower bound)
        mask = (x.s{3}(:,k) < 0);
        if any(mask)
            x.s{3}(mask,k) = 0;

            anti_windup = mask & (g.ivm.ivmud_con(:,8) == 1) & (x.ds{3}(:,k) < 0);
            if any(anti_windup)
                x.ds{3}(anti_windup,k) = 0;
            end
        end

        P3 = g.ivm.ivmud_con(:,6) ...
             .*(g.ivm.ivmud_con(:,5) - x.s{8}(:,k)) + x.s{3}(:,k);

        P3 = max(P3,0);

        % rate of change of the commanded angle
        x.ds{1}(:,k) = g.sys.basrad.*(P1 + P2 + P3);

        %---------------------------------------------------------------------%
        % Q-V droop control loop

        Q1 = g.ivm.ivmud_con(:,9).*(x.qset(:,k) - x.s{9}(:,k));

        x.ds{6}(:,k) = g.ivm.ivmud_con(:,16).*(g.ivm.ivmud_con(:,10) - x.s{9}(:,k));

        % anti-windup limits
        mask = (x.s{6}(:,k) > 0);
        if any(mask)
            x.s{6}(mask,k) = 0;

            anti_windup = mask & (g.ivm.ivmud_con(:,17) == 1) ...
                          & (x.ds{6}(:,k) > 0);

            if any(anti_windup)
                x.ds{6}(anti_windup,k) = 0;
            end
        end

        Q2 = g.ivm.ivmud_con(:,15).*(g.ivm.ivmud_con(:,10) - x.s{9}(:,k)) ...
             + x.s{6}(:,k);

        Q2 = min(Q2,0);

        x.ds{7}(:,k) = g.ivm.ivmud_con(:,16).*(g.ivm.ivmud_con(:,11) - x.s{9}(:,k));

        mask = (x.s{7}(:,k) < 0);
        if any(mask)
            x.s{7}(mask,k) = 0;

            anti_windup = mask & (g.ivm.ivmud_con(:,17) == 1) ...
                          & (x.ds{7}(:,k) < 0);

            if any(anti_windup)
                x.ds{7}(anti_windup,k) = 0;
            end
        end

        Q3 = g.ivm.ivmud_con(:,15).*(g.ivm.ivmud_con(:,11) - x.s{9}(:,k)) ...
             + x.s{7}(:,k);

        Q3 = max(Q3,0);

        x.vref = x.vset(:,k) + Q1 + Q2 + Q3;

        x.ds{5}(:,k) = g.ivm.ivmud_con(:,13).*(x.vref - x.s{10}(:,k));

        % voltage integral state upper bound
        mask = (x.s{5}(:,k) > g.ivm.ivmud_con(:,18));
        if any(mask)
            x.s{5}(mask,k) = g.ivm.ivmud_con(mask,18);

            anti_windup = mask & (g.ivm.ivmud_con(:,17) == 1) ...
                          & (x.ds{5}(:,k) > 0);

            if any(anti_windup)
                x.ds{5}(anti_windup,k) = 0;
            end
        end

        % voltage integral state lower bound
        mask = (x.s{5}(:,k) < g.ivm.ivmud_con(:,19));
        if any(mask)
            x.s{5}(mask,k) = g.ivm.ivmud_con(mask,19);

            anti_windup = mask & (g.ivm.ivmud_con(:,17) == 1) ...
                          & (x.ds{5}(:,k) < 0);

            if any(anti_windup)
                x.ds{5}(anti_windup,k) = 0;
            end
        end

        %---------------------------------------------------------------------%
        % commanded voltage magnitude time constant

        x.ds{4}(:,k) = (x.vref - x.s{4}(:,k))./max(g.ivm.ivmud_con(:,14),lbnd);

        mask = (g.ivm.ivmud_con(:,25) == 1);
        if any(mask)
            x.ds{4}(mask,k) = (g.ivm.ivmud_con(mask,12) ...
                               .*(x.vref(mask) - x.s{10}(mask,k)) ...
                               + x.s{5}(mask,k) - x.s{4}(mask,k)) ...
                              ./max(g.ivm.ivmud_con(mask,14),lbnd);
        end

        mask = (g.ivm.ivmud_con(:,14) < lbnd);
        if any(mask)
            x.ds{4}(mask,k) = 0.0;
            x.s{4}(mask,k) = x.vref(mask);

            pi_cascade = mask & (g.ivm.ivmud_con(:,25) == 1);
            if any(pi_cascade)
                x.s{4}(pi_cascade,k) = g.ivm.ivmud_con(pi_cascade,12) ...
                                       .*(x.vref(pi_cascade) ...
                                          - x.s{10}(pi_cascade,k)) ...
                                       + x.s{5}(pi_cascade,k);
            end
        end

        x.s{4}(:,k) = min(x.s{4}(:,k),g.ivm.ivmud_con(:,18));
        x.s{4}(:,k) = max(x.s{4}(:,k),g.ivm.ivmud_con(:,19));
    end
end

end  % function end

% eof
