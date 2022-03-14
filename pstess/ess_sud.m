function [x] = ess_sud(i,k,bus,flag,x)
% Syntax: x = ess_sud(i,k,bus,flag,x)
%
% Purpose: ess user-defined control model
%
% Input:   i - ess number
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
% See Also: ess

%-----------------------------------------------------------------------------%
% Version history
%
% Version: Example
% Author:  Ryan T. Elliott
% Date:    June 2021
%
% Note:    In this example, ess_sud.m is a simplified version of the
%          native PSTess model lsc.m. See data/d2asbegp_ess.m for
%          information about the data and parameters stored in
%          g.ess.essud_con.
%
%          NOTE: A given energy storage system should only have a single
%          controller providing modulation commands to it. For example,
%          an energy storage device may receive commands from ess_sud
%          or a native PSTess model, such as lsc.m, but not both.
%
% Ref:     For additional information about the control strategy
%          implemented here, see the paper below (available free on
%          arXiv at https://arxiv.org/abs/2012.03161).
%
%          R. T. Elliott, P. Arabshahi, and D. S. Kirschen,
%          "Stabilizing transient disturbances with utility-scale
%          inverter-based resources," IET Generation, Transmission &
%          Distribution, vol. 14, pp. 6534–6544, December 2020.
%
%-----------------------------------------------------------------------------%
% ess_sud states
%   var      description
%   x.s{1}   transducer for remote angle signal
%   x.s{2}   transducer for local angle signal
%   x.s{3}   local LTV highpass filter state 1
%   x.s{4}   local LTV highpass filter state 2
%   x.s{5}   local LTV lead-lag stage 1
%   x.s{6}   local LTV lead-lag stage 2
%   x.s{7}   center LTI highpass filter state 1
%   x.s{8}   center LTI highpass filter state 2
%   x.s{9}   center LTI lead-lag stage 1
%   x.s{10}  center LTI lead-lag stage 2
%   x.s{11}  lowpass filter
%
% storage variables
%   var      description
%   x.v{1}   recorded angle measurements
%   x.v{2}   apparent power command (complex)
%-----------------------------------------------------------------------------%

global g;  % declaring struct of global variables

if ~isempty(g.ess.essud_con)
    lbnd = 1e-3;  % lower bound to prevent division by zero

    % essud_pot(:,1) -- conversion factor into system MVA base
    essud_pot(:,1) = g.ess.ess_con(g.ess.dessud_idx,6)/g.sys.basmva;

    % sensor_set -- buses with angle sensors (internal)
    % sud_sensor_idx -- indexes in sensor_set with ess_sud controllers
    sensor_set = [3, 5, 9, 10, 12, 13];  % external buses 3, 10, 14, 20, 110, 120
    sud_sensor_idx = [2, 4, 5, 6];       % external buses 10, 20, 110, 120
    n_sensor = length(sensor_set);

    if (flag == 0)  % initialization
        if (i ~= 0)
            % caution: non-vectorized computation is not supported!
            error('ess_sud: initialization must be vectorized!');
        else
            % vectorized computation (initialization)

            % busnum -- ess_sud bus number vector
            busnum = g.bus.bus_int(g.ess.essud_con(:,2));

            % x.s{1} is assigned in s_simu as z_essud where
            % z_essud = zeros(max(g.ess.n_essud,1),n_ktime)
            z_essud = x.s{1};                           % initialization vector
            n_ktime = size(z_essud,2);                  % number of time indexes

            % initialize the algebraic storage variables
            x.v{1} = zeros(n_sensor,n_ktime);           % theta_sensor
            x.v{2} = zeros(g.ess.n_essud,n_ktime);      % apparent power command

            % v1, theta_sensor(:,1) -- initial sensor voltage angles (rad)
            x.v{1}(:,1) = bus(sensor_set,3)*pi/180;

            % theta_coi_ini -- initial center-of-inertia angle
            % theta_sud_ini -- initial ess_sud voltage angle
            theta_coi_ini = sum(x.v{1}(:,1))/n_sensor;
            theta_sud_ini = x.v{1}(sud_sensor_idx,1);

            % initialize the dynamic states
            for idx = 1:11                              % 11 dynamic states
                x.s{idx} = z_essud;
                x.ds{idx} = z_essud;
            end

            x.s{1}(:,1) = theta_coi_ini;                % center-of-inertia angle
            x.s{2}(:,1) = theta_sud_ini;                % ess_sud angles

            % parameter checks
            mask = (g.ess.essud_con(:,13) >= g.ess.essud_con(:,14));
            if any(mask)
                estr = '\ness_sud: impermissible local ltv modulation limits at ';
                estr = [estr, 'bus %0.0f, min >= max.'];
                error(sprintf(estr,bus(busnum(mask),1)));
            end

            mask = (g.ess.essud_con(:,24) >= g.ess.essud_con(:,25));
            if any(mask)
                estr = '\ness_sud: impermissible center lti modulation limits at ';
                estr = [estr, 'bus %0.0f, min >= max.'];
                error(sprintf(estr,bus(busnum(mask),1)));
            end

            mask = (g.ess.essud_con(:,28) >= g.ess.essud_con(:,29));
            if any(mask)
                estr = '\ness_sud: impermissible total modulation limits at ';
                estr = [estr, 'bus %0.0f, min >= max.'];
                error(sprintf(estr,bus(busnum(mask),1)));
            end
        end
    end

    if (flag == 1)  % network interface computation
                    % interface calculation handled by ess
    end

    if (flag == 2)  % dynamics calculation
        if (i ~= 0)
            % caution: non-vectorized computation is not supported!
            error('ess_sud: dynamics calculation must be vectorized!');
        else
            % vectorized computation (dynamics calculation)

            % unwrapping the phase angles (phase tracking)
            %
            % Algorithm adapted from:
            %   "Analysis of the phase unwrapping algorithm," K. Itoh,
            %   Journal of Applied Optics, 1982.

            k_i = max(1,k-1);  % initial sample index for unwrapping

            theta_sensor_raw = angle(g.bus.bus_v(sensor_set,k));

            % step 1 -- take the diff between the last 2 samples
            % step 2 -- wrap the diff between +pi and -pi
            % step 3 -- perform discrete-time integration to unwrap the signal
            angle_jump = theta_sensor_raw - x.v{1}(:,k_i);
            angle_jump = angle_jump - 2*pi*floor((angle_jump+pi)/(2*pi));
            x.v{1}(:,k) = x.v{1}(:,k_i) + angle_jump;

            %-----------------------------------------------------------------%
            % remote and local angle transducers

            % removing buses with low voltage from the center-of-inertia calculation
            safe_volt_mask = (abs(g.bus.bus_v(sensor_set,k)) >= 0.70);

            u_theta_coi = sum(x.v{1}(safe_volt_mask,k))/sum(safe_volt_mask);
            u_theta_sud = x.v{1}(sud_sensor_idx,k);

            x.ds{1}(:,k) = (u_theta_coi - x.s{1}(:,k))./max(g.ess.essud_con(:,3),lbnd);
            x.ds{2}(:,k) = (u_theta_sud - x.s{2}(:,k))./max(g.ess.essud_con(:,3),lbnd);

            mask = (g.ess.essud_con(:,3) < lbnd);
            if any(mask)  % integrator bypass
                x.ds{1}(mask,k) = 0.0;
                x.ds{2}(mask,k) = 0.0;

                % x.s{1} -- center-of-inertia angle
                % x.s{2} -- local voltage angle
                x.s{1}(mask,k) = u_theta_coi(mask);
                x.s{2}(mask,k) = u_theta_sud(mask);
            end

            u_coi_set = x.s{1}(:,1);                    % center of angle setpoint
            u_local_set = x.s{2}(:,1);                  % local angle setpoint

            %-----------------------------------------------------------------%
            % local LTV compensation path (alpha 1)

            u_a1hp = u_local_set - u_coi_set + x.s{1}(:,k) - x.s{2}(:,k);

            % washout filter
            x.ds{3}(:,k) = g.ess.essud_con(:,5).*u_a1hp ...
                           - g.ess.essud_con(:,6).*x.s{3}(:,k) + x.s{4}(:,k);

            x.ds{4}(:,k) = g.ess.essud_con(:,7).*u_a1hp - g.ess.essud_con(:,8).*x.s{3}(:,k);

            u_a1c1 = g.ess.essud_con(:,4).*(u_a1hp + x.s{3}(:,k));

            % lead-lag compensator
            x.ds{5}(:,k) = (u_a1c1 - x.s{5}(:,k))./max(g.ess.essud_con(:,10),lbnd);

            mask = (g.ess.essud_con(:,10) < lbnd);
            if any(mask)                                % integrator bypass
                x.ds{5}(mask,k) = 0.0;
                x.s{5}(mask,k) = u_a1c1(mask);          % lead-lag state
            end

            tmp_a1c1 = x.s{5}(:,k) ...
                       .*(1 - g.ess.essud_con(:,9)./max(g.ess.essud_con(:,10),lbnd));

            u_a1c2 = tmp_a1c1 ...
                     + u_a1c1.*(g.ess.essud_con(:,9)./max(g.ess.essud_con(:,10),lbnd));

            x.ds{6}(:,k) = (u_a1c2 - x.s{6}(:,k))./max(g.ess.essud_con(:,12),lbnd);

            mask = (g.ess.essud_con(:,12) < lbnd);
            if any(mask)                                % integrator bypass
                x.ds{6}(mask,k) = 0.0;
                x.s{6}(mask,k) = u_a1c2(mask);          % lead-lag state
            end

            tmp_a1c2 = x.s{6}(:,k) ...
                       .*(1 - g.ess.essud_con(:,11)./max(g.ess.essud_con(:,12),lbnd));

            y_a1c2 = tmp_a1c2 ...
                     + u_a1c2.*(g.ess.essud_con(:,11)./max(g.ess.essud_con(:,12),lbnd));

            y_a1 = max(y_a1c2,g.ess.essud_con(:,13).*essud_pot(:,1));
            y_a1 = min(y_a1,g.ess.essud_con(:,14).*essud_pot(:,1));

            %-----------------------------------------------------------------%
            % center LTI compensation path (alpha 2)

            u_a2hp = u_coi_set - x.s{1}(:,k);

            % washout filter
            x.ds{7}(:,k) = g.ess.essud_con(:,16).*u_a2hp ...
                           - g.ess.essud_con(:,17).*x.s{7}(:,k) + x.s{8}(:,k);

            x.ds{8}(:,k) = g.ess.essud_con(:,18).*u_a2hp - g.ess.essud_con(:,19).*x.s{7}(:,k);

            u_a2c1 = g.ess.essud_con(:,15).*(u_a2hp + x.s{7}(:,k));

            % lead-lag compensator
            x.ds{9}(:,k) = (u_a2c1 - x.s{9}(:,k))./max(g.ess.essud_con(:,21),lbnd);

            mask = (g.ess.essud_con(:,21) < lbnd);
            if any(mask)                                % integrator bypass
                x.ds{9}(mask,k) = 0.0;
                x.s{9}(mask,k) = u_a2c1(mask);          % lead-lag state
            end

            tmp_a2c1 = x.s{9}(:,k) ...
                       .*(1 - g.ess.essud_con(:,20)./max(g.ess.essud_con(:,21),lbnd));

            u_a2c2 = tmp_a2c1 ...
                     + u_a2c1.*(g.ess.essud_con(:,20)./max(g.ess.essud_con(:,21),lbnd));

            x.ds{10}(:,k) = (u_a2c2 - x.s{10}(:,k))./max(g.ess.essud_con(:,23),lbnd);

            mask = (g.ess.essud_con(:,23) < lbnd);
            if any(mask)                                % integrator bypass
                x.ds{10}(mask,k) = 0.0;
                x.s{10}(mask,k) = u_a2c2(mask);         % lead-lag state
            end

            tmp_a2c2 = x.s{10}(:,k) ...
                       .*(1 - g.ess.essud_con(:,22)./max(g.ess.essud_con(:,23),lbnd));

            y_a2c2 = tmp_a2c2 ...
                     + u_a2c2.*(g.ess.essud_con(:,22)./max(g.ess.essud_con(:,23),lbnd));

            y_a2 = max(y_a2c2,g.ess.essud_con(:,24).*essud_pot(:,1));
            y_a2 = min(y_a2,g.ess.essud_con(:,25).*essud_pot(:,1));

            %-----------------------------------------------------------------%
            % computing the ess_sud control command

            u_lp = g.ess.essud_con(:,26).*(y_a1 + y_a2);

            x.ds{11}(:,k) = (u_lp - x.s{11}(:,k))./max(g.ess.essud_con(:,27),lbnd);

            mask = (g.ess.essud_con(:,27) < lbnd);
            if any(mask)                                % integrator bypass
                x.ds{11}(mask,k) = 0.0;
                x.s{11}(mask,k) = u_lp(mask);           % lowpass state
            end

            pcmd_lb = g.ess.essud_con(:,28).*essud_pot(:,1);
            pcmd_ub = g.ess.essud_con(:,29).*essud_pot(:,1);

            % active power command
            ess_sud_pcmd = max(x.s{11}(:,k),pcmd_lb);
            ess_sud_pcmd = min(ess_sud_pcmd,pcmd_ub);

            % reactive power command
            ess_sud_qcmd = zeros(size(ess_sud_pcmd));

            % sending the modulation command to the inverter-based resource
            x.v{2}(:,k) = ess_sud_pcmd + 1j*ess_sud_qcmd;
            g.ess.ess_dsig(g.ess.dessud_idx,k) = x.v{2}(:,k);

        end
    end
end

end  % function end

% eof
