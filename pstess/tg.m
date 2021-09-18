function tg(i,k,flag)
% Syntax: tg(i,k,flag)
%
% Purpose: Simple turbine governor model
%
% Input:   i - generator number
%          k - integer time
%          flag - 0 - initialization
%                 1 - network interface computation
%                 2 - system dynamics computation
%
% Data format for tg_con
%          col    data                                unit
%            1    turbine model number (=1)
%            2    machine number
%            3    speed set point, wf                 pu
%            4    steady state gain, 1/R              pu
%            5    maximum power order, Tmax           pu on generator base
%            6    servo time constant, Ts             sec
%            7    governor time constant, Tc          sec
%            8    transient gain time constant, T3    sec
%            9    HP section time constant, T4        sec
%           10    reheater time constant T5           sec

%-----------------------------------------------------------------------------%
% Version history
%
% Version: 1.0
% Author:  Joe H. Chow
% Date:    August 1993
%-----------------------------------------------------------------------------%

global g;  % declaring struct of global variables

if (flag == 0)  % initialization
    if (i ~= 0)
        if (g.tg.tg_con(i,1) ~= 1)
            error('tg: it is required that tg_con(i,1) = 1.');
        end
    end

    if (i ~= 0)  % scalar computation
        n = g.mac.mac_int(g.tg.tg_con(i,2));  % machine number

        if (g.mac.pmech(n,k) > g.tg.tg_con(i,5))
            estr = '\ntg: pmech exceeds maximum limit at ';
            estr = [estr, 'generator %0.0f.'];
            error(sprintf(estr,n));
        end
        %
        if (g.mac.pmech(n,k) < 0)
            estr = '\ntg: pmech less than zero at ';
            estr = [estr, 'generator %0.0f.'];
            error(sprintf(estr,n));
        end

        g.tg.tg1(i,1) = g.mac.pmech(n,k);
        %
        g.tg.tg_pot(i,1) = g.tg.tg_con(i,8)/g.tg.tg_con(i,7);
        a1 = 1 - g.tg.tg_pot(i,1);
        g.tg.tg_pot(i,2) = a1;
        g.tg.tg2(i,1) = a1*g.mac.pmech(n,k);
        %
        g.tg.tg_pot(i,3) = g.tg.tg_con(i,9)/g.tg.tg_con(i,10);
        a2 = 1 - g.tg.tg_pot(i,3);
        g.tg.tg_pot(i,4) = a2;
        g.tg.tg3(i,1) = a2*g.mac.pmech(n,k);
        %
        g.tg.tg_pot(i,5) = g.mac.pmech(n,k);
        %
        g.tg.tg_sig(i,1) = 0;
    else
        % vectorized computation
        if (g.tg.n_tg ~= 0)
            n = g.mac.mac_int(g.tg.tg_con(g.tg.tg_idx,2));  % machine number

            maxlmt = find(g.mac.pmech(n,1)>g.tg.tg_con(g.tg.tg_idx,5));
            if ~isempty(maxlmt)
                estr = '\ntg: pmech exceeds maximum limit at ';
                estr = [estr, 'generator %0.0f.'];
                error(sprintf(estr,n(maxlmt)));
            end

            minlmt = find(g.mac.pmech(n,1)<zeros(g.tg.n_tg,1));
            if ~isempty(minlmt)
                estr = '\ntg: pmech less than zero at ';
                estr = [estr, 'generator %0.0f.'];
                error(sprintf(estr,n(minlmt)));
            end

            g.tg.tg1(g.tg.tg_idx,1) = g.mac.pmech(n,1);
            %
            g.tg.tg_pot(g.tg.tg_idx,1) = g.tg.tg_con(g.tg.tg_idx,8) ...
                                         ./g.tg.tg_con(g.tg.tg_idx,7);
            a1 = ones(g.tg.n_tg,1) - g.tg.tg_pot(g.tg.tg_idx,1);
            g.tg.tg_pot(g.tg.tg_idx,2) = a1;
            g.tg.tg2(g.tg.tg_idx,1) = a1.*g.mac.pmech(n,k);
            %
            g.tg.tg_pot(g.tg.tg_idx,3) = g.tg.tg_con(g.tg.tg_idx,9) ...
                                         ./g.tg.tg_con(g.tg.tg_idx,10);
            a2 = ones(g.tg.n_tg,1) - g.tg.tg_pot(g.tg.tg_idx,3);
            g.tg.tg_pot(g.tg.tg_idx,4) = a2;
            g.tg.tg3(g.tg.tg_idx,1) = a2.*g.mac.pmech(n,k);
            %
            g.tg.tg_pot(g.tg.tg_idx,5) = g.mac.pmech(n,k);  % set reference value
            g.tg.tg_sig(g.tg.tg_idx,1) = zeros(g.tg.n_tg,1);
        end
    end
end

if (flag == 1)   % network interface computation
    if (i ~= 0)  % scalar computation
        n = g.mac.mac_int(g.tg.tg_con(i,2));  % machine number
        % needed because pmech depends on tg1, tg2, and tg3
        g.mac.pmech(n,k) = ...
            g.tg.tg3(i,k) + g.tg.tg_pot(i,3) ...
                            *(g.tg.tg2(i,k) + g.tg.tg_pot(i,1)*g.tg.tg1(i,k));
    else
        if (g.tg.n_tg ~= 0)
            n = g.mac.mac_int(g.tg.tg_con(g.tg.tg_idx,2));  % machine number
            g.mac.pmech(n,k) = ...
                g.tg.tg3(g.tg.tg_idx,k) ...
                + g.tg.tg_pot(g.tg.tg_idx,3) ...
                  .*(g.tg.tg2(g.tg.tg_idx,k) ...
                     + g.tg.tg_pot(g.tg.tg_idx,1).*g.tg.tg1(g.tg.tg_idx,k));
        end
    end
end

if (flag == 2)   % turbine governor dynamics calculation
    if (i ~= 0)  % scalar computation
        n = g.mac.mac_int(g.tg.tg_con(i,2));  % machine number
        spd_err = g.tg.tg_con(i,3) - g.mac.mac_spd(n,k);
        demand = g.tg.tg_pot(i,5) + spd_err*g.tg.tg_con(i,4) + g.tg.tg_sig(i,k);
        demand = min(max(demand,0),g.tg.tg_con(i,5));

        g.tg.dtg1(i,k) = (demand - g.tg.tg1(i,k))/g.tg.tg_con(i,6);
        %
        g.tg.dtg2(i,k) = (g.tg.tg_pot(i,2)*g.tg.tg1(i,k) - g.tg.tg2(i,k)) ...
                         /g.tg.tg_con(i,7);
        %
        g.tg.dtg3(i,k) = ((g.tg.tg2(i,k) ...
                           + g.tg.tg_pot(i,1)*g.tg.tg1(i,k))*g.tg.tg_pot(i,4) ...
                          - g.tg.tg3(i,k))/g.tg.tg_con(i,10);

        g.mac.pmech(n,k) = g.tg.tg3(i,k) ...
                           + g.tg.tg_pot(i,3) ...
                             *(g.tg.tg2(i,k) + g.tg.tg_pot(:,1)*g.tg.tg1(i,k));
    else
        % vectorized computation
        if (g.tg.n_tg ~= 0)
            n = g.mac.mac_int(g.tg.tg_con(g.tg.tg_idx,2));  % machine number
            spd_err = g.tg.tg_con(g.tg.tg_idx,3) - g.mac.mac_spd(n,k);
            demand = g.tg.tg_pot(g.tg.tg_idx,5) ...
                     + spd_err.*g.tg.tg_con(g.tg.tg_idx,4) ...
                     + g.tg.tg_sig(g.tg.tg_idx,k);
            demand = min(max(demand,zeros(g.tg.n_tg,1)),g.tg.tg_con(g.tg.tg_idx,5));

            g.tg.dtg1(g.tg.tg_idx,k) = (demand - g.tg.tg1(g.tg.tg_idx,k)) ...
                                       ./g.tg.tg_con(g.tg.tg_idx,6);
            %
            g.tg.dtg2(g.tg.tg_idx,k) = ...
                (g.tg.tg1(g.tg.tg_idx,k).*g.tg.tg_pot(g.tg.tg_idx,2) ...
                 - g.tg.tg2(g.tg.tg_idx,k))./g.tg.tg_con(g.tg.tg_idx,7);
            %
            g.tg.dtg3(g.tg.tg_idx,k) = ...
                ((g.tg.tg2(g.tg.tg_idx,k) ...
                  + g.tg.tg_pot(g.tg.tg_idx,1).*g.tg.tg1(g.tg.tg_idx,k)) ...
                 .*g.tg.tg_pot(g.tg.tg_idx,4) - g.tg.tg3(g.tg.tg_idx,k)) ...
                ./g.tg.tg_con(g.tg.tg_idx,10);

            g.mac.pmech(n,k) = ...
                g.tg.tg3(g.tg.tg_idx,k) ...
                + g.tg.tg_pot(g.tg.tg_idx,3) ...
                  .*(g.tg.tg2(g.tg.tg_idx,k) ...
                     + g.tg.tg_pot(g.tg.tg_idx,1).*g.tg.tg1(g.tg.tg_idx,k));
        end
    end
end

end  % function end

% eof
