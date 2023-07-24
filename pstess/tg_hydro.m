function tg_hydro(i,k,flag)
% Syntax: tg_hydro(i,k,flag)
%
% Purpose: hydraulic turbine governor model
%
% Input:   i - generator number
%          k - integer time
%          flag - 0 - initialization
%                 1 - network interface computation
%                 2 - system dynamics computation
%
% Data format for tg_con
%          col    data                          unit
%            1    turbine model number (=2)
%            2    machine number
%            3    speed set point, wf           pu
%            4    permanent droop, Rp           pu
%            5    transient droop, Rt           pu
%            6    maximum power order, Tmax     pu on generator base
%            7    maximum rate limit            pu on gen base/sec
%            8    minimum rate limit            pu on gen base per sec
%            9    servo time constant, Ts       sec
%           10    servo gain, Ks
%           11    governor time constant, Tg    sec
%           12    reset time constant, Tr       sec
%           13    water starting time, Tw       sec

%-----------------------------------------------------------------------------%
% Version history
%
% Version: 1.0
% Date:    June 1998
% Author:  Graham Rogers
% Purpose: Model of hydraulic turbine governor
%-----------------------------------------------------------------------------%

global g;  % declaring struct of global variables

if (g.tg.n_tgh ~= 0)
    if (flag == 0)  % initialization
        if (i ~= 0)
            error('tg_hydro: initialization must be vectorized.');
        else
            % vectorized computation
            if (g.tg.n_tgh ~= 0)
                n = g.mac.mac_int(g.tg.tg_con(g.tg.tgh_idx,2));  % machine number

                maxlmt = find(g.mac.pmech(n,1)>g.tg.tg_con(g.tg.tgh_idx,6));
                if ~isempty(maxlmt)
                    estr = '\ntg_hydro: pmech exceeds maximum limit at ';
                    estr = [estr, 'generator %0.0f.'];
                    error(sprintf(estr,n(maxlmt)));
                end

                minlmt = find(g.mac.pmech(n,1)<zeros(g.tg.n_tgh,1));
                if ~isempty(minlmt)
                    estr = '\ntg_hydro: pmech less than zero at ';
                    estr = [estr, 'generator %0.0f.'];
                    error(sprintf(estr,n(minlmt)));
                end

                g.tg.tg1(g.tg.tgh_idx,1) = zeros(g.tg.n_tgh,1);
                g.tg.tg2(g.tg.tgh_idx,1) = g.mac.pmech(n,1);
                g.tg.tg3(g.tg.tgh_idx,1) = g.mac.pmech(n,1);
                g.tg.tg4(g.tg.tgh_idx,1) = g.mac.pmech(n,1);
                g.tg.tg5(g.tg.tgh_idx,1) = 3*g.mac.pmech(n,1);

                % reference value
                g.tg.tg_pot(g.tg.tgh_idx,5) = g.mac.pmech(n,1) ...
                                              .*g.tg.tg_con(g.tg.tgh_idx,4);

                g.tg.tg_sig(g.tg.tgh_idx,1) = zeros(g.tg.n_tgh,1);
            end
        end
    end

    if (flag == 1)   % network interface computation
        if (i ~= 0)  % scalar computation
            % vector computation only for tg_hydro
            error('tg_hydro: network interface calculation must be vectorized.');
        else
            if (g.tg.n_tgh ~= 0)
                n = g.mac.mac_int(g.tg.tg_con(g.tg.tgh_idx,2));  % machine number
                g.mac.pmech(n,k) = g.tg.tg5(g.tg.tgh_idx,k) ...
                                   - 2*g.tg.tg4(g.tg.tgh_idx,k);
            end
        end
    end

    if (flag == 2)   % turbine governor dynamics calculation
        if (i ~= 0)  % scalar computation
            % vector computation only
            error('tg_hydro: dynamics calculation must be vectorized.');
        else
            % vectorized computation
            if (g.tg.n_tgh ~= 0)
                n = g.mac.mac_int(g.tg.tg_con(g.tg.tgh_idx,2));  % machine number

                spd_err = g.tg.tg_con(g.tg.tgh_idx,3) - g.mac.mac_spd(n,k);
                demand = spd_err + g.tg.tg_sig(g.tg.tgh_idx,k) ...
                         + g.tg.tg_pot(g.tg.tgh_idx,5);

                g.tg.dtg1(g.tg.tgh_idx,k) = ...
                    g.tg.tg_con(g.tg.tgh_idx,10) ...
                    .*(demand-g.tg.tg1(g.tg.tgh_idx,k) ...
                       - (g.tg.tg_con(g.tg.tgh_idx,4) ...
                          + g.tg.tg_con(g.tg.tgh_idx,5)) ...
                         .*g.tg.tg2(g.tg.tgh_idx,k) ...
                       + g.tg.tg_con(g.tg.tgh_idx,5).*g.tg.tg3(g.tg.tgh_idx,k)) ...
                    ./g.tg.tg_con(g.tg.tgh_idx,9);

                % apply rate limit
                rmax = find(g.tg.dtg1(g.tg.tgh_idx,k) > g.tg.tg_con(g.tg.tgh_idx,7));
                if ~isempty(rmax)
                    g.tg.dtg1(g.tg.tgh_idx(rmax),k) = ...
                        g.tg.tg_con(g.tg.tgh_idx(rmax),7);
                end

                rmin = find(g.tg.dtg1(g.tg.tgh_idx,k) < g.tg.tg_con(g.tg.tgh_idx,8));
                if ~isempty(rmin)
                    g.tg.dtg1(g.tg.tgh_idx(rmin),k) = ...
                        g.tg.tg_con(g.tg.tgh_idx(rmin),8);
                end

                g.tg.dtg2(g.tg.tgh_idx,k) = g.tg.tg1(g.tg.tgh_idx,k);

                % check non-wind-up limit
                smax = find(g.tg.tg2(g.tg.tgh_idx,k) >= g.tg.tg_con(g.tg.tgh_idx,6));
                if ~isempty(smax)
                    g.tg.tg2(g.tg.tgh_idx(smax),k) = ...
                        g.tg.tg_con(g.tg.tgh_idx(smax),6);

                    if (g.tg.dtg2(g.tg.tgh_idx(smax),k) > 0)
                        g.tg.dtg2(g.tg.tgh_idx(smax),k) = zeros(length(smax),1);
                    end
                end

                smin = find(g.tg.tg2(g.tg.tgh_idx,k) <= 0);
                if ~isempty(smin)
                    g.tg.tg2(g.tg.tgh_idx(smax),k) = zeros(length(smin),1);
                    if (g.tg.dtg2(g.tg.tgh_idx(smin),k) < 0)
                        g.tg.dtg2(g.tg.tgh_idx(smin),k) = zeros(length(smin),1);
                    end
                end

                % transient droop
                g.tg.dtg3(g.tg.tgh_idx,k) = ...
                    (g.tg.tg2(g.tg.tgh_idx,k) - g.tg.tg3(g.tg.tgh_idx,k)) ...
                    ./g.tg.tg_con(g.tg.tgh_idx,12);

                % gate servo
                g.tg.dtg4(g.tg.tgh_idx,k) = ...
                    (g.tg.tg2(g.tg.tgh_idx,k) - g.tg.tg4(g.tg.tgh_idx,k)) ...
                    ./g.tg.tg_con(g.tg.tgh_idx,11);

                % hydraulic turbine
                g.tg.dtg5(g.tg.tgh_idx,k) = ...
                    2*(3*g.tg.tg4(g.tg.tgh_idx,k) - g.tg.tg5(g.tg.tgh_idx,k)) ...
                    ./g.tg.tg_con(g.tg.tgh_idx,13);

                g.mac.pmech(n,k) = g.tg.tg5(g.tg.tgh_idx,k) ...
                                   - 2*g.tg.tg4(g.tg.tgh_idx,k);
            end
        end
    end
end

end  % function end

% eof
