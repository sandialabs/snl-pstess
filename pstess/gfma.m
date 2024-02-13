function gfma(i,k,bus,flag)
% Syntax: gfma(i,k,bus,flag)
%
% Purpose: Droop-based grid-forming inverter control model
%          WECC-approved as REGFM_A1, an augmented version of the
%          CERTS droop model
%
%          NOTE - there must be a matching ivm model that accompanies gfma
%
% Input:   i - gfma controller number (index)
%              if i = 0, vectorized computation
%          k - integer time
%          flag - 0 - initialization
%                 1 - network interface computation
%                 2 - generator dynamics computation
%                 3 - generator dynamics computation for linearization
%
% Note:    The gfma power base is inherited from the ivm to which the
%          controller is connected.
%
% gfma_con matrix format
%    col    data                                              units
%     1     gfma controller number                            integer
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
%    20     inverter maximum output current                   pu
%    21     active power transducer time constant             s
%    22     reactive power transducer time constant           s
%    23     voltage magnitude transducer time constant        s
%    24     Q-V initialization flag, qvflag (note 1)          binary
%    25     voltage control flag, vflag (note 2)              binary
%    26     current limiting state behavior flag (note 3)     binary
%
% note 1: when gfma_con(,24) = 0, the reactive power setpoint is
%         initialized based on the power flow solution and Vset
%         is calculated accordingly; however, when
%         gfma_con(,24) = 1, the reactive power setpoint is set
%         to 0.
%
% note 2: when gfma_con(,25) = 0, the PI section of the voltage
%         control loop is bypassed; however, when
%         gfma_con(,25) = 1, the commanded voltage magnitude
%         is taken from the output of the PI section (after
%         limiting).
%
% note 3: when gfma_con(,26) = 0, the control states are not
%         updated when the current limit is reached; however,
%         when gfma_con(,26) = 1, the control states
%         (gfma1, gfma4, gfma5) are updated to reflect the
%         modified commanded voltage magnitude and angle.
%
% gfma states
%   var     description
%   gfma1   commanded voltage angle (integrator)
%   gfma2   max P overload mitigation control state
%   gfma3   min P overload mitigation control state
%   gfma4   commanded voltage magnitude (first-order)
%   gfma5   voltage regulation integral control state
%   gfma6   max Q overload mitigation control state
%   gfma7   min Q overload mitigation control state
%   gfma8   inverter active power transducer
%   gfma9   inverter reactive power transducer
%   gfma10  inverter voltage transducer

%-----------------------------------------------------------------------------%
% Version history
%
% Version: 1.0
% Author:  Ryan T. Elliott
% Date:    October 2022
% Note:    Initial version
%-----------------------------------------------------------------------------%

global g;  % declaring struct of global variables

lbnd = 1e-3;  % lower bound to prevent division by zero

if (~isempty(g.gfma.gfma_con) && i ~= 0)
    % non-vectorized computation
    error('gfma: all gfma calculations must be vectorized.');
elseif (~isempty(g.gfma.gfma_con) && i == 0)
    if (flag == 0)  % initialization
        % busnum -- vector of bus numbers with gfma controllers
        busnum = g.bus.bus_int(g.gfma.gfma_con(:,2));

        g.gfma.gfma1(:,1) = g.mac.mac_ang(g.gfma.mac_idx,1);
        g.gfma.gfma2(:,1) = 0.0;
        g.gfma.gfma3(:,1) = 0.0;
        g.gfma.gfma4(:,1) = g.mac.eqprime(g.gfma.mac_idx,1);
        g.gfma.gfma5(:,1) = g.mac.eqprime(g.gfma.mac_idx,1);
        g.gfma.gfma6(:,1) = 0.0;
        g.gfma.gfma7(:,1) = 0.0;
        g.gfma.gfma8(:,1) = g.mac.pelect(g.gfma.mac_idx,1) ...  % P_inv
                            .*g.mac.mac_pot(g.gfma.mac_idx,1);
        g.gfma.gfma9(:,1) = g.mac.qelect(g.gfma.mac_idx,1) ...  % Q_inv
                            .*g.mac.mac_pot(g.gfma.mac_idx,1);
        g.gfma.gfma10(:,1) = g.mac.eterm(g.gfma.mac_idx,1);     % V_inv

        g.gfma.dgfma1(:,1) = 0.0;   % s1  -- commanded voltage angle (integrator)
        g.gfma.dgfma2(:,1) = 0.0;   % s2  -- max P overload mitigation control state
        g.gfma.dgfma3(:,1) = 0.0;   % s3  -- min P overload mitigation control state
        g.gfma.dgfma4(:,1) = 0.0;   % s4  -- commanded voltage magnitude (first-order)
        g.gfma.dgfma5(:,1) = 0.0;   % s5  -- voltage regulation integral control state
        g.gfma.dgfma6(:,1) = 0.0;   % s6  -- max Q overload mitigation control state
        g.gfma.dgfma7(:,1) = 0.0;   % s7  -- min Q overload mitigation control state
        g.gfma.dgfma8(:,1) = 0.0;   % s8  -- inverter active power transducer
        g.gfma.dgfma9(:,1) = 0.0;   % s9  -- inverter reactive power transducer
        g.gfma.dgfma10(:,1) = 0.0;  % s10 -- inverter voltage transducer

        g.gfma.pset = diag(g.gfma.gfma8(:,1))*ones(size(g.gfma.pset));
        g.gfma.qset = zeros(size(g.gfma.qset));

        % g.gfma.vref is a single column, not a matrix
        g.gfma.vref = g.gfma.gfma4(:,1);

        mask = (g.gfma.gfma_con(:,25) == 1);
        if any(mask)
            g.gfma.vref(mask) = g.gfma.gfma10(mask,1);          % pi control
        end

        g.gfma.vset = ...
            diag(g.gfma.vref + g.gfma.gfma_con(:,9).*g.gfma.gfma9(:,1)) ...
            *ones(size(g.gfma.vset));

        mask = (g.gfma.gfma_con(:,24) == 0);
        if any(mask)
            g.gfma.qset(mask,:) = diag(g.gfma.gfma9(mask,1)) ...
                                  *ones(size(g.gfma.qset(mask,:)));
            g.gfma.vset(mask,:) = diag(g.gfma.vref(mask)) ...
                                  *ones(size(g.gfma.vset(mask,:)));
        end

        %-----------------------------------------------------------------%
        % initialization checks

        mask = (g.gfma.gfma8(:,1) > g.gfma.gfma_con(:,4));
        if any(mask)
            estr = '\ngfma: controller initialized above Pmax at ';
            estr = [estr, 'bus %0.0f '];
            error(sprintf(estr,bus(busnum(mask),1)));
        end

        mask = (g.gfma.gfma8(:,1) < g.gfma.gfma_con(:,5));
        if any(mask)
            estr = '\ngfma: controller initialized below Pmin at ';
            estr = [estr, 'bus %0.0f '];
            error(sprintf(estr,bus(busnum(mask),1)));
        end

        mask = (g.gfma.gfma9(:,1) > g.gfma.gfma_con(:,10));
        if any(mask)
            estr = '\ngfma: controller initialized above Qmax at ';
            estr = [estr, 'bus %0.0f '];
            error(sprintf(estr,bus(busnum(mask),1)));
        end

        mask = (g.gfma.gfma9(:,1) < g.gfma.gfma_con(:,11));
        if any(mask)
            estr = '\ngfma: controller initialized below Qmin at ';
            estr = [estr, 'bus %0.0f '];
            error(sprintf(estr,bus(busnum(mask),1)));
        end

        %-----------------------------------------------------------------%
        % parameter checks

        mask = (g.gfma.gfma_con(:,4) < g.gfma.gfma_con(:,5));
        if any(mask)
            estr = '\ngfma: invalid active power limits selected at ';
            estr = [estr, 'bus %0.0f '];
            error(sprintf(estr,bus(busnum(mask),1)));
        end

        mask = (g.gfma.gfma_con(:,10) < g.gfma.gfma_con(:,11));
        if any(mask)
            estr = '\ngfma: invalid reactive power limits selected at ';
            estr = [estr, 'bus %0.0f '];
            error(sprintf(estr,bus(busnum(mask),1)));
        end

        mask = (g.gfma.gfma_con(:,18) < g.gfma.gfma_con(:,19));
        if any(mask)
            estr = '\ngfma: invalid commanded voltage magnitude limits ';
            estr = [estr, 'selected at bus %0.0f '];
            error(sprintf(estr,bus(busnum(mask),1)));
        end

        %-----------------------------------------------------------------%
        % flag checks

        mask = ~ismember(g.gfma.gfma_con(:,8),[0,1]);
        if any(mask)
            estr = '\ngfma: invalid active overload windup flag selected at ';
            estr = [estr, 'bus %0.0f '];
            error(sprintf(estr,bus(busnum(mask),1)));
        end

        mask = ~ismember(g.gfma.gfma_con(:,17),[0,1]);
        if any(mask)
            estr = '\ngfma: invalid reactive overload windup flag selected at ';
            estr = [estr, 'bus %0.0f '];
            error(sprintf(estr,bus(busnum(mask),1)));
        end

        mask = ~ismember(g.gfma.gfma_con(:,24),[0,1]);
        if any(mask)
            estr = '\ngfma: invalid Q-V initialization flag (qvflag) selected at';
            estr = [estr, 'bus %0.0f '];
            error(sprintf(estr,bus(busnum(mask),1)));
        end

        mask = ~ismember(g.gfma.gfma_con(:,25),[0,1]);
        if any(mask)
            estr = '\ngfma: invalid voltage control flag (vflag) selected at ';
            estr = [estr, 'bus %0.0f '];
            error(sprintf(estr,bus(busnum(mask),1)));
        end

        % setting kiv to zero if the PI section is bypassed
        mask = ((g.gfma.gfma_con(:,25) == 0) & (g.gfma.gfma_con(:,13) > 0));
        if any(mask)
            g.gfma.gfma_con(mask,13) = 0;
            wstr = '\ngfma: setting gfma kiv integral gain to zero at ';
            wstr = [wstr, 'bus %0.0f '];
            warning(sprintf(wstr,bus(busnum(mask),1)));
        end

        mask = ~ismember(g.gfma.gfma_con(:,26),[0,1]);
        if any(mask)
            estr = '\ngfma: invalid limiting state behavior flag selected at ';
            estr = [estr, 'bus %0.0f '];
            error(sprintf(estr,bus(busnum(mask),1)));
        end
    end

    if (flag == 1)  % network interface
        %---------------------------------------------------------------------%
        % current limiting

        mask = (g.gfma.gfma_con(:,14) < lbnd);
        if any(mask)
            g.gfma.gfma4(mask,k) = g.gfma.vref(mask);

            pi_cascade = mask & (g.gfma.gfma_con(:,25) == 1);
            if any(pi_cascade)
                g.gfma.gfma4(pi_cascade,k) = g.gfma.gfma_con(pi_cascade,12) ...
                                             .*(g.gfma.vref(pi_cascade) ...
                                                - g.gfma.gfma10(pi_cascade,k)) ...
                                             + g.gfma.gfma5(pi_cascade,k);
            end

            g.gfma.gfma4(mask,k) = min(g.gfma.gfma4(mask,k), ...
                                       g.gfma.gfma_con(mask,18));

            g.gfma.gfma4(mask,k) = max(g.gfma.gfma4(mask,k), ...
                                       g.gfma.gfma_con(mask,19));
        end

        g.mac.vex(g.gfma.mac_idx,k) = g.gfma.gfma4(:,k);
        g.mac.fldcur(g.gfma.mac_idx,k) = g.gfma.gfma1(:,k);

        % busnum -- vector of bus numbers where gfma models are connected
        busnum = g.bus.bus_int(g.gfma.gfma_con(:,2));

        Vt = g.bus.bus_v(busnum,k);

        E = g.gfma.gfma4(:,k).*exp(1j*g.gfma.gfma1(:,k));

        It = (E - Vt)./(g.mac.mac_con(g.gfma.mac_idx,5) ...
                        + 1j*g.mac.mac_con(g.gfma.mac_idx,7));

        mask = (abs(It) > g.gfma.gfma_con(:,20));
        if any(mask)
            It(mask) = (It(mask)./abs(It(mask))).*g.gfma.gfma_con(mask,20);

            E(mask) = Vt(mask) + It(mask) ...
                                 .*(g.mac.mac_con(g.gfma.mac_idx(mask),5) ...
                                    + 1j*g.mac.mac_con(g.gfma.mac_idx(mask),7));

            lim_state = mask & (g.gfma.gfma_con(:,26) == 1);
            if any(lim_state)
                g.gfma.gfma1(lim_state,k) = angle(E(lim_state));
                g.gfma.gfma4(lim_state,k) = abs(E(lim_state));
                g.gfma.gfma5(lim_state,k) = abs(E(lim_state)) ...
                                            - g.gfma.gfma_con(lim_state,12) ...
                                              .*(g.gfma.vref(lim_state) ...
                                                 - g.gfma.gfma10(lim_state,k));

                % zero rate conditions would go here
                %
                % if foo
                %     g.gfma.dgfma1(lim_state,k) = 0;
                %     g.gfma.dgfma4(lim_state,k) = 0;
                %     g.gfma.dgfma5(lim_state,k) = 0;
                % end
            end

            g.mac.vex(g.gfma.mac_idx(mask),k) = abs(E(mask));
            g.mac.fldcur(g.gfma.mac_idx(mask),k) = angle(E(mask));
        end

        %---------------------------------------------------------------------%
        % checking the limit condition

        curg_mag = abs(g.mac.cur_re(g.gfma.mac_idx,k) ...
                       + 1j*g.mac.cur_im(g.gfma.mac_idx,k)) ...
                   .*g.mac.mac_pot(g.gfma.mac_idx,1);

        if (sum(curg_mag) ~= 0)  % 0.5 pct tolerance
            g.gfma.lim_flag = (curg_mag > 1.005*g.gfma.gfma_con(:,20));
        end
    end

    if (flag == 2 || flag == 3)  % dynamics calculation
        %---------------------------------------------------------------------%
        % transducers

        % active power
        g.gfma.dgfma8(:,k) = (g.mac.pelect(g.gfma.mac_idx,k) ...
                              .*g.mac.mac_pot(g.gfma.mac_idx,1) ...
                              - g.gfma.gfma8(:,k))./max(g.gfma.gfma_con(:,21),lbnd);

        mask = (g.gfma.gfma_con(:,21) < lbnd);
        if any(mask)
            g.gfma.dgfma8(mask,k) = 0.0;
            g.gfma.gfma8(mask,k) = g.mac.pelect(g.gfma.mac_idx(mask),k) ...
                                   .*g.mac.mac_pot(g.gfma.mac_idx(mask),1);
        end

        % reactive power
        g.gfma.dgfma9(:,k) = (g.mac.qelect(g.gfma.mac_idx,k) ...
                              .*g.mac.mac_pot(g.gfma.mac_idx,1) ...
                              - g.gfma.gfma9(:,k))./max(g.gfma.gfma_con(:,22),lbnd);

        mask = (g.gfma.gfma_con(:,22) < lbnd);
        if any(mask)
            g.gfma.dgfma9(mask,k) = 0.0;
            g.gfma.gfma9(mask,k) = g.mac.qelect(g.gfma.mac_idx(mask),k) ...
                                   .*g.mac.mac_pot(g.gfma.mac_idx(mask),1);
        end

        % voltage magnitude
        g.gfma.dgfma10(:,k) = ...
            (g.mac.eterm(g.gfma.mac_idx,k) - g.gfma.gfma10(:,k)) ...
            ./max(g.gfma.gfma_con(:,23),lbnd);

        mask = (g.gfma.gfma_con(:,23) < lbnd);
        if any(mask)
            g.gfma.dgfma10(mask,k) = 0.0;
            g.gfma.gfma10(mask,k) = g.mac.eterm(g.gfma.mac_idx(mask),k);
        end

        %---------------------------------------------------------------------%
        % P-f droop control loop

        P1 = g.gfma.gfma_con(:,3).*(g.gfma.pset(:,k) - g.gfma.gfma8(:,k) ...
                                    + real(g.gfma.gfma_sig(:,k)));

        g.gfma.dgfma2(:,k) = g.gfma.gfma_con(:,7) ...
                             .*(g.gfma.gfma_con(:,4) - g.gfma.gfma8(:,k));

        % anti-windup limits (upper bound)
        mask = (g.gfma.gfma2(:,k) > 0);
        if any(mask)
            g.gfma.gfma2(mask,k) = 0;

            anti_windup = mask & (g.gfma.gfma_con(:,8) == 1) ...
                          & (g.gfma.dgfma2(:,k) > 0);

            if any(anti_windup)
                g.gfma.dgfma2(anti_windup,k) = 0;
            end
        end

        P2 = g.gfma.gfma_con(:,6).*(g.gfma.gfma_con(:,4) - g.gfma.gfma8(:,k)) ...
             + g.gfma.gfma2(:,k);

        P2 = min(P2,0);

        g.gfma.dgfma3(:,k) = g.gfma.gfma_con(:,7) ...
                             .*(g.gfma.gfma_con(:,5) - g.gfma.gfma8(:,k));

        % anti-windup limits (lower bound)
        mask = (g.gfma.gfma3(:,k) < 0);
        if any(mask)
            g.gfma.gfma3(mask,k) = 0;

            anti_windup = mask & (g.gfma.gfma_con(:,8) == 1) ...
                          & (g.gfma.dgfma3(:,k) < 0);

            if any(anti_windup)
                g.gfma.dgfma3(anti_windup,k) = 0;
            end
        end

        P3 = g.gfma.gfma_con(:,6) ...
             .*(g.gfma.gfma_con(:,5) - g.gfma.gfma8(:,k)) + g.gfma.gfma3(:,k);

        P3 = max(P3,0);

        % rate of change of the commanded angle
        g.gfma.dgfma1(:,k) = g.sys.basrad.*(P1 + P2 + P3);

        %---------------------------------------------------------------------%
        % Q-V droop control loop

        Q1 = g.gfma.gfma_con(:,9).*(g.gfma.qset(:,k) - g.gfma.gfma9(:,k) ...
                                    + imag(g.gfma.gfma_sig(:,k)));

        g.gfma.dgfma6(:,k) = g.gfma.gfma_con(:,16) ...
                             .*(g.gfma.gfma_con(:,10) - g.gfma.gfma9(:,k));

        % anti-windup limits
        mask = (g.gfma.gfma6(:,k) > 0);
        if any(mask)
            g.gfma.gfma6(mask,k) = 0;

            anti_windup = mask & (g.gfma.gfma_con(:,17) == 1) ...
                          & (g.gfma.dgfma6(:,k) > 0);

            if any(anti_windup)
                g.gfma.dgfma6(anti_windup,k) = 0;
            end
        end

        Q2 = g.gfma.gfma_con(:,15).*(g.gfma.gfma_con(:,10) - g.gfma.gfma9(:,k)) ...
             + g.gfma.gfma6(:,k);

        Q2 = min(Q2,0);

        g.gfma.dgfma7(:,k) = g.gfma.gfma_con(:,16) ...
                             .*(g.gfma.gfma_con(:,11) - g.gfma.gfma9(:,k));

        mask = (g.gfma.gfma7(:,k) < 0);
        if any(mask)
            g.gfma.gfma7(mask,k) = 0;

            anti_windup = mask & (g.gfma.gfma_con(:,17) == 1) ...
                          & (g.gfma.dgfma7(:,k) < 0);

            if any(anti_windup)
                g.gfma.dgfma7(anti_windup,k) = 0;
            end
        end

        Q3 = g.gfma.gfma_con(:,15).*(g.gfma.gfma_con(:,11) - g.gfma.gfma9(:,k)) ...
             + g.gfma.gfma7(:,k);

        Q3 = max(Q3,0);

        g.gfma.vref = g.gfma.vset(:,k) + Q1 + Q2 + Q3;

        g.gfma.dgfma5(:,k) = g.gfma.gfma_con(:,13) ...
                             .*(g.gfma.vref - g.gfma.gfma10(:,k));

        % voltage integral state upper bound
        mask = (g.gfma.gfma5(:,k) > g.gfma.gfma_con(:,18));
        if any(mask)
            g.gfma.gfma5(mask,k) = g.gfma.gfma_con(mask,18);

            anti_windup = mask & (g.gfma.gfma_con(:,17) == 1) ...
                          & (g.gfma.dgfma5(:,k) > 0);

            if any(anti_windup)
                g.gfma.dgfma5(anti_windup,k) = 0;
            end
        end

        % voltage integral state lower bound
        mask = (g.gfma.gfma5(:,k) < g.gfma.gfma_con(:,19));
        if any(mask)
            g.gfma.gfma5(mask,k) = g.gfma.gfma_con(mask,19);

            anti_windup = mask & (g.gfma.gfma_con(:,17) == 1) ...
                          & (g.gfma.dgfma5(:,k) < 0);

            if any(anti_windup)
                g.gfma.dgfma5(anti_windup,k) = 0;
            end
        end

        %---------------------------------------------------------------------%
        % commanded voltage magnitude time constant

        g.gfma.dgfma4(:,k) = (g.gfma.vref - g.gfma.gfma4(:,k)) ...
                             ./max(g.gfma.gfma_con(:,14),lbnd);

        mask = (g.gfma.gfma_con(:,25) == 1);
        if any(mask)
            g.gfma.dgfma4(mask,k) = (g.gfma.gfma_con(mask,12) ...
                                     .*(g.gfma.vref(mask) ...
                                        - g.gfma.gfma10(mask,k)) ...
                                     + g.gfma.gfma5(mask,k) ...
                                     - g.gfma.gfma4(mask,k)) ...
                                    ./max(g.gfma.gfma_con(mask,14),lbnd);
        end

        mask = (g.gfma.gfma_con(:,14) < lbnd);
        if any(mask)
            g.gfma.dgfma4(mask,k) = 0.0;
            % gfma4 is updated in the network interface when bypassed
        end

        g.gfma.gfma4(:,k) = min(g.gfma.gfma4(:,k),g.gfma.gfma_con(:,18));
        g.gfma.gfma4(:,k) = max(g.gfma.gfma4(:,k),g.gfma.gfma_con(:,19));
    end
end

end  % function end

% eof
