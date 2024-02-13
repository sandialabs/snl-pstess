function lsc(i,k,bus,flag)
% Syntax: lsc(i,k,bus,flag)
%
% Purpose: LTV synchronizing torque controller,
%          with vectorized computation option
%          NOTE - there must be a matching inverter-based resource,
%                 such as an energy storage system
%
% Input:   i - LTV synchronizing torque control number (index)
%              if i = 0, vectorized computation
%          k - integer time
%          bus - solved loadflow bus data
%          flag - 0 - initialization
%                 1 - network interface computation
%                 2 - generator dynamics computation
%                 3 - generator dynamics computation for linearization
%
% Note:    The lsc power base is inherited from the ess to which the
%          controller is connected.
%
% lsc_con matrix format
%    col   data                                        units
%     1    LTV synchronizing controller index number   integer
%     2    bus number                                  integer
%     3    angle transducer time constant              sec
%     4    use Pade approximation flag (0 = bypass)    binary
%     5    remote signal time delay (Pade)             sec
%     6    local signal time delay (Pade)              sec
%     7    remote angle setpoint rate limit            rad/sec
%     8    remote angle setpoint time constant         sec
%     9    local angle setpoint rate limit             rad/sec
%    10    local angle setpoint time constant          sec
%    11    local LTV tuning weight (alpha 1)           pu
%    12    local LTV highpass numerator 1              sec
%    13    local LTV highpass denominator 1            sec
%    14    local LTV highpass numerator 2              sec
%    15    local LTV highpass denominator 2            sec
%    16    local LTV lead-lag numerator 1              sec
%    17    local LTV lead-lag denominator 1            sec
%    18    local LTV lead-lag numerator 2              sec
%    19    local LTV lead-lag denominator 2            sec
%    20    local LTV modulation command lower bound    pu on lsc base
%    21    local LTV modulation command upper bound    pu on lsc base
%    22    center LTI tuning weight (alpha 2)          pu
%    23    center LTI highpass numerator 1             sec
%    24    center LTI highpass denominator 1           sec
%    25    center LTI highpass numerator 2             sec
%    26    center LTI highpass denominator 2           sec
%    27    center LTI lead-lag numerator 1             sec
%    28    center LTI lead-lag denominator 1           sec
%    29    center LTI lead-lag numerator 2             sec
%    30    center LTI lead-lag denominator 2           sec
%    31    center LTI modulation command lower bound   pu on lsc base
%    32    center LTI modulation command upper bound   pu on lsc base
%    33    transient stability control gain            pu
%    34    lowpass filter time constant                sec
%    35    total modulation command lower bound        pu on lsc base
%    36    total modulation command upper bound        pu on lsc base
%
% lsc states
%   var    description
%   lsc1   transducer for remote angle signal
%   lsc2   transducer for local angle signal
%   lsc3   Pade approx. for remote angle signal
%   lsc4   Pade approx. for local angle signal
%   lsc5   remote angle signal setpoint tracking filter
%   lsc6   local angle signal setpoint tracking filter
%   lsc7   local LTV highpass filter state 1
%   lsc8   local LTV highpass filter state 2
%   lsc9   local LTV lead-lag stage 1
%   lsc10  local LTV lead-lag stage 2
%   lsc11  center LTI highpass filter state 1
%   lsc12  center LTI highpass filter state 2
%   lsc13  center LTI lead-lag stage 1
%   lsc14  center LTI lead-lag stage 2
%   lsc15  lowpass filter

%-----------------------------------------------------------------------------%
% Version history
%
% Version: 2.1
% Author:  Ryan T. Elliott
% Date:    November 2020
% Note:    Split the control functionality into a dedicated model
%
% Version: 2.0
% Author:  Ryan T. Elliott
% Date:    September 2019
% Note:    Modified compensation paths and implemented setpoint tracking
%
% Version: 1.0
% Author:  Ryan T. Elliott
% Date:    August 2019
% Note:    Initial version
%-----------------------------------------------------------------------------%

global g;  % declaring struct of global variables

lbnd = 1e-3;  % lower bound to prevent division by zero

if ~isempty(g.lsc.lsc_con)
    if (flag == 0)  % initialization
        if (i ~= 0)
            % caution: non-vectorized computation is not supported!
            error('lsc: initialization must be vectorized!');
        else
            % vectorized computation (initialization)

            % busnum -- lsc bus number vector
            busnum = g.bus.bus_int(g.lsc.lsc_con(:,2));

            % theta_sensor(:,1) -- initial sensor voltage angles (rad)
            g.lsc.theta_sensor(:,1) = bus(g.lsc.sensor_set,3)*pi/180;

            % theta_coi_ini -- initial center-of-inertia angle
            % theta_lsc_ini -- initial lsc voltage angle
            theta_coi_ini = sum(g.lsc.theta_sensor(:,1))/g.lsc.n_sensor;
            theta_lsc_ini = g.lsc.theta_sensor(g.lsc.lsc_sensor_idx,1);

            % theta_coi_pade(:,1) -- initial delayed center-of-inertia angle (Pade)
            % theta_lsc_pade(:,1) -- initial delayed lsc voltage angle (Pade)
            g.lsc.theta_coi_pade(:,1) = theta_coi_ini;
            g.lsc.theta_lsc_pade(:,1) = theta_lsc_ini;

            % lsc_pot(:,1) -- lsc power capacity (pu on syst. base)
            g.lsc.lsc_pot(:,1) = g.ess.ess_con(g.lsc.lsc_ess_idx,6)/g.sys.basmva;

            % state variables
            g.lsc.lsc1(:,1) = theta_coi_ini;
            g.lsc.lsc2(:,1) = theta_lsc_ini;
            g.lsc.lsc3(:,1) = theta_coi_ini;
            g.lsc.lsc4(:,1) = theta_lsc_ini;
            g.lsc.lsc5(:,1) = theta_coi_ini;
            g.lsc.lsc6(:,1) = theta_lsc_ini;
            g.lsc.lsc7(:,1) = 0.0;
            g.lsc.lsc8(:,1) = 0.0;
            g.lsc.lsc9(:,1) = 0.0;
            g.lsc.lsc10(:,1) = 0.0;
            g.lsc.lsc11(:,1) = 0.0;
            g.lsc.lsc12(:,1) = 0.0;
            g.lsc.lsc13(:,1) = 0.0;
            g.lsc.lsc14(:,1) = 0.0;
            g.lsc.lsc15(:,1) = 0.0;

            % state derivatives
            g.lsc.dlsc1(:,1) = 0.0;
            g.lsc.dlsc2(:,1) = 0.0;
            g.lsc.dlsc3(:,1) = 0.0;
            g.lsc.dlsc4(:,1) = 0.0;
            g.lsc.dlsc5(:,1) = 0.0;
            g.lsc.dlsc6(:,1) = 0.0;
            g.lsc.dlsc7(:,1) = 0.0;
            g.lsc.dlsc8(:,1) = 0.0;
            g.lsc.dlsc9(:,1) = 0.0;
            g.lsc.dlsc10(:,1) = 0.0;
            g.lsc.dlsc11(:,1) = 0.0;
            g.lsc.dlsc12(:,1) = 0.0;
            g.lsc.dlsc13(:,1) = 0.0;
            g.lsc.dlsc14(:,1) = 0.0;
            g.lsc.dlsc15(:,1) = 0.0;

            % parameter checks
            mask = (g.lsc.lsc_con(:,20) >= g.lsc.lsc_con(:,21));
            if any(mask)
                estr = '\nlsc: impermissible local ltv modulation limits at ';
                estr = [estr, 'bus %0.0f, min >= max.'];
                error(sprintf(estr,bus(busnum(mask),1)));
            end

            mask = (g.lsc.lsc_con(:,31) >= g.lsc.lsc_con(:,32));
            if any(mask)
                estr = '\nlsc: impermissible center lti modulation limits at ';
                estr = [estr, 'bus %0.0f, min >= max.'];
                error(sprintf(estr,bus(busnum(mask),1)));
            end

            mask = (g.lsc.lsc_con(:,35) >= g.lsc.lsc_con(:,36));
            if any(mask)
                estr = '\nlsc: impermissible total modulation limits at ';
                estr = [estr, 'bus %0.0f, min >= max.'];
                error(sprintf(estr,bus(busnum(mask),1)));
            end
        end
    end

    if (flag == 1)  % network interface computation
        % interface calculation handled by ess
    end

    if (flag == 2 || flag == 3)  % dynamics calculation
        if (i ~= 0)
            % caution: non-vectorized computation is not supported!
            error('lsc: dynamics calculation must be vectorized!');
        else
            % vectorized computation (dynamics calculation)

            % unwrapping the phase angles (phase tracking)
            %
            % Algorithm adapted from:
            %   "Analysis of the phase unwrapping algorithm," K. Itoh,
            %   Journal of Applied Optics, 1982.

            k_i = max(1,k-1);  % initial sample index for unwrapping

            theta_sensor_raw = angle(g.bus.bus_v(g.lsc.sensor_set,k));

            % step 1 -- take the diff between the last 2 samples
            % step 2 -- wrap the diff between +pi and -pi
            % step 3 -- perform discrete-time integration to unwrap the signal
            angle_jump = theta_sensor_raw - g.lsc.theta_sensor(:,k_i);
            angle_jump = angle_jump - 2*pi*floor((angle_jump+pi)/(2*pi));
            g.lsc.theta_sensor(:,k) = g.lsc.theta_sensor(:,k_i) + angle_jump;

            %-----------------------------------------------------------------%
            % remote and local angle transducers

            % removing buses with low voltage from the center-of-inertia calculation
            safe_volt_mask = (abs(g.bus.bus_v(g.lsc.sensor_set,k)) >= 0.70);

            u_theta_coi = sum(g.lsc.theta_sensor(safe_volt_mask,k)) ...
                          /sum(safe_volt_mask);

            u_theta_lsc = g.lsc.theta_sensor(g.lsc.lsc_sensor_idx,k);

            g.lsc.dlsc1(:,k) = (u_theta_coi - g.lsc.lsc1(:,k)) ...
                               ./max(g.lsc.lsc_con(:,3),lbnd);

            g.lsc.dlsc2(:,k) = (u_theta_lsc - g.lsc.lsc2(:,k)) ...
                               ./max(g.lsc.lsc_con(:,3),lbnd);

            mask = (g.lsc.lsc_con(:,3) < lbnd);
            if any(mask)  % integrator bypass
                g.lsc.dlsc1(mask,k) = 0.0;
                g.lsc.dlsc2(mask,k) = 0.0;

                % lsc1 -- center-of-inertia angle
                % lsc2 -- local voltage angle
                g.lsc.lsc1(mask,k) = u_theta_coi(mask);
                g.lsc.lsc2(mask,k) = u_theta_lsc(mask);
            end

            %-----------------------------------------------------------------%
            % time delay via Pade approximation

            % input is the unwrapped center-of-inertia angle
            g.lsc.dlsc3(:,k) = (g.lsc.lsc1(:,k) - g.lsc.lsc3(:,k)) ...
                               ./max(g.lsc.lsc_con(:,5)/2,lbnd);

            mask = (g.lsc.lsc_con(:,5)/2 < lbnd);
            if any(mask)                                  % integrator bypass
                g.lsc.dlsc3(mask,k) = 0.0;
                g.lsc.lsc3(mask,k) = g.lsc.lsc1(mask,k);  % average voltage angle
            end

            g.lsc.theta_coi_pade(:,k) = 2*g.lsc.lsc3(:,k) - g.lsc.lsc1(:,k);

            % input is the unwrapped lsc voltage angle
            g.lsc.dlsc4(:,k) = (g.lsc.lsc2(:,k) - g.lsc.lsc4(:,k)) ...
                               ./max(g.lsc.lsc_con(:,6)/2,lbnd);

            mask = (g.lsc.lsc_con(:,6)/2 < lbnd);
            if any(mask)                                  % integrator bypass
                g.lsc.dlsc4(mask,k) = 0.0;
                g.lsc.lsc4(mask,k) = g.lsc.lsc2(mask,k);  % local voltage angle
            end

            g.lsc.theta_lsc_pade(:,k) = 2*g.lsc.lsc4(:,k) - g.lsc.lsc2(:,k);

            %-----------------------------------------------------------------%
            % setpoint tracking

            % remote setpoint tracking filter input
            u_coi_set = g.lsc.lsc1(:,k);                  % center of angle input

            mask = (g.lsc.lsc_con(:,4) ~= 0);             % check Pade flag
            if any(mask)
                u_coi_set(mask) = g.lsc.theta_coi_pade(mask,k);
            end

            % rate limiting the remote setpoint
            u_coi_err = u_coi_set - g.lsc.lsc5(:,k);

            u_coi_err_ub = g.lsc.lsc_con(:,7);
            u_coi_err_lb = -g.lsc.lsc_con(:,7);

            u_coi_err = max(u_coi_err,u_coi_err_lb);
            u_coi_err = min(u_coi_err,u_coi_err_ub);

            u_coi_set = u_coi_err + g.lsc.lsc5(:,k);

            g.lsc.dlsc5(:,k) = (u_coi_set - g.lsc.lsc5(:,k)) ...
                               ./max(g.lsc.lsc_con(:,8),lbnd);

            mask = (g.lsc.lsc_con(:,8) < lbnd);
            if any(mask)                                  % integrator bypass
                g.lsc.dlsc5(mask,k) = 0.0;
                g.lsc.lsc5(mask,k) = u_coi_set(mask);     % center of angle setpoint
            end

            % local setpoint tracking filter input
            u_local_set = g.lsc.lsc2(:,k);                % local angle input

            mask = (g.lsc.lsc_con(:,4) ~= 0);             % check Pade flag
            if any(mask)
                u_local_set(mask) = g.lsc.theta_lsc_pade(mask,k);
            end

            % rate limiting the local setpoint
            u_local_err = u_local_set - g.lsc.lsc6(:,k);

            u_local_err_ub = g.lsc.lsc_con(:,9);
            u_local_err_lb = -g.lsc.lsc_con(:,9);

            u_local_err = max(u_local_err,u_local_err_lb);
            u_local_err = min(u_local_err,u_local_err_ub);

            u_local_set = u_local_err + g.lsc.lsc6(:,k);

            g.lsc.dlsc6(:,k) = (u_local_set - g.lsc.lsc6(:,k)) ...
                               ./max(g.lsc.lsc_con(:,10),lbnd);

            mask = (g.lsc.lsc_con(:,10) < lbnd);
            if any(mask)                                  % integrator bypass
                g.lsc.dlsc6(mask,k) = 0.0;
                g.lsc.lsc6(mask,k) = u_local_set(mask);   % local angle setpoint
            end

            %-----------------------------------------------------------------%
            % local LTV compensation path (alpha 1)

            u_a1hp = g.lsc.lsc6(:,k) - g.lsc.lsc5(:,k) ...
                     - g.lsc.lsc2(:,k) + g.lsc.lsc1(:,k);

            mask = (g.lsc.lsc_con(:,4) ~= 0);             % check Pade flag
            if any(mask)
                u_a1hp(mask) = g.lsc.lsc6(mask,k) - g.lsc.lsc5(mask,k) ...
                               - g.lsc.theta_lsc_pade(mask,k) ...
                               + g.lsc.theta_coi_pade(mask,k);
            end

            % washout filter
            g.lsc.dlsc7(:,k) = g.lsc.lsc_con(:,12).*u_a1hp ...
                               - g.lsc.lsc_con(:,13).*g.lsc.lsc7(:,k) ...
                               + g.lsc.lsc8(:,k);

            g.lsc.dlsc8(:,k) = g.lsc.lsc_con(:,14).*u_a1hp ...
                               - g.lsc.lsc_con(:,15).*g.lsc.lsc7(:,k);

            u_a1c1 = g.lsc.lsc_con(:,11).*(u_a1hp + g.lsc.lsc7(:,k));

            % lead-lag compensator
            g.lsc.dlsc9(:,k) = (u_a1c1 - g.lsc.lsc9(:,k)) ...
                               ./max(g.lsc.lsc_con(:,17),lbnd);

            mask = (g.lsc.lsc_con(:,17) < lbnd);
            if any(mask)                                  % integrator bypass
                g.lsc.dlsc9(mask,k) = 0.0;
                g.lsc.lsc9(mask,k) = u_a1c1(mask);        % lead-lag state
            end

            tmp_a1c1 = g.lsc.lsc9(:,k) ...
                       .*(1 - g.lsc.lsc_con(:,16)./max(g.lsc.lsc_con(:,17),lbnd));

            u_a1c2 = tmp_a1c1 ...
                     + u_a1c1.*(g.lsc.lsc_con(:,16)./max(g.lsc.lsc_con(:,17),lbnd));

            g.lsc.dlsc10(:,k) = (u_a1c2 - g.lsc.lsc10(:,k)) ...
                                ./max(g.lsc.lsc_con(:,19),lbnd);

            mask = (g.lsc.lsc_con(:,19) < lbnd);
            if any(mask)                                  % integrator bypass
                g.lsc.dlsc10(mask,k) = 0.0;
                g.lsc.lsc10(mask,k) = u_a1c2(mask);       % lead-lag state
            end

            tmp_a1c2 = g.lsc.lsc10(:,k) ...
                       .*(1 - g.lsc.lsc_con(:,18)./max(g.lsc.lsc_con(:,19),lbnd));

            y_a1c2 = tmp_a1c2 ...
                     + u_a1c2.*(g.lsc.lsc_con(:,18)./max(g.lsc.lsc_con(:,19),lbnd));

            y_a1 = max(y_a1c2,g.lsc.lsc_con(:,20).*g.lsc.lsc_pot(:,1));
            y_a1 = min(y_a1,g.lsc.lsc_con(:,21).*g.lsc.lsc_pot(:,1));

            %-----------------------------------------------------------------%
            % center LTI compensation path (alpha 2)

            u_a2hp = g.lsc.lsc5(:,k) - g.lsc.lsc1(:,k);

            mask = (g.lsc.lsc_con(:,4) ~= 0);             % check Pade flag
            if any(mask)
                u_a2hp(mask) = g.lsc.lsc5(mask,k) - g.lsc.theta_coi_pade(mask,k);
            end

            % washout filter
            g.lsc.dlsc11(:,k) = g.lsc.lsc_con(:,23).*u_a2hp ...
                                - g.lsc.lsc_con(:,24).*g.lsc.lsc11(:,k) ...
                                + g.lsc.lsc12(:,k);

            g.lsc.dlsc12(:,k) = g.lsc.lsc_con(:,25).*u_a2hp ...
                                - g.lsc.lsc_con(:,26).*g.lsc.lsc11(:,k);

            u_a2c1 = g.lsc.lsc_con(:,22).*(u_a2hp + g.lsc.lsc11(:,k));

            % lead-lag compensator
            g.lsc.dlsc13(:,k) = (u_a2c1 - g.lsc.lsc13(:,k)) ...
                                ./max(g.lsc.lsc_con(:,28),lbnd);

            mask = (g.lsc.lsc_con(:,28) < lbnd);
            if any(mask)                                  % integrator bypass
                g.lsc.dlsc13(mask,k) = 0.0;
                g.lsc.lsc13(mask,k) = u_a2c1(mask);       % lead-lag state
            end

            tmp_a2c1 = g.lsc.lsc13(:,k) ...
                       .*(1 - g.lsc.lsc_con(:,27)./max(g.lsc.lsc_con(:,28),lbnd));

            u_a2c2 = tmp_a2c1 ...
                     + u_a2c1.*(g.lsc.lsc_con(:,27)./max(g.lsc.lsc_con(:,28),lbnd));

            g.lsc.dlsc14(:,k) = (u_a2c2 - g.lsc.lsc14(:,k)) ...
                                ./max(g.lsc.lsc_con(:,30),lbnd);

            mask = (g.lsc.lsc_con(:,30) < lbnd);
            if any(mask)                                  % integrator bypass
                g.lsc.dlsc14(mask,k) = 0.0;
                g.lsc.lsc14(mask,k) = u_a2c2(mask);       % lead-lag state
            end

            tmp_a2c2 = g.lsc.lsc14(:,k) ...
                       .*(1 - g.lsc.lsc_con(:,29)./max(g.lsc.lsc_con(:,30),lbnd));

            y_a2c2 = tmp_a2c2 ...
                     + u_a2c2.*(g.lsc.lsc_con(:,29)./max(g.lsc.lsc_con(:,30),lbnd));

            y_a2 = max(y_a2c2,g.lsc.lsc_con(:,31).*g.lsc.lsc_pot(:,1));
            y_a2 = min(y_a2,g.lsc.lsc_con(:,32).*g.lsc.lsc_pot(:,1));

            %-----------------------------------------------------------------%
            % computing the lsc control command

            u_lp = g.lsc.lsc_con(:,33).*(y_a1 + y_a2);

            g.lsc.dlsc15(:,k) = (u_lp - g.lsc.lsc15(:,k)) ...
                                ./max(g.lsc.lsc_con(:,34),lbnd);

            mask = (g.lsc.lsc_con(:,34) < lbnd);
            if any(mask)                                  % integrator bypass
                g.lsc.dlsc15(mask,k) = 0.0;
                g.lsc.lsc15(mask,k) = u_lp(mask);         % lowpass state
            end

            pcmd_lb = g.lsc.lsc_con(:,35).*g.lsc.lsc_pot(:,1);
            pcmd_ub = g.lsc.lsc_con(:,36).*g.lsc.lsc_pot(:,1);

            % active power command
            lsc_pcmd = max(g.lsc.lsc15(:,k),pcmd_lb);
            lsc_pcmd = min(lsc_pcmd,pcmd_ub);

            % reactive power command
            lsc_qcmd = zeros(size(lsc_pcmd));

            % sending the modulation command to the inverter-based resource
            g.lsc.lsc_scmd(:,k) = lsc_pcmd + 1j*lsc_qcmd;
            g.ess.ess_scmd(g.lsc.lsc_ess_idx,k) = g.lsc.lsc_scmd(g.lsc.lsc_idx,k);

        end
    end
end

end  % function end

% eof
