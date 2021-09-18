function [bus_new] = svc(i,k,bus,flag,v_sbus)
% Syntax: [bus_new] = svc(i,k,bus,flag,v_sbus)
%
% Purpose: static var system,
%          with vectorized computation option
%
%          NOTE - static var bus must be declared as a
%                 non-conforming load bus
%
% Input: i - static var number
%            if i= 0, vectorized computation
%        k - integer time
%        bus - solved loadflow bus data
%        flag - 0 - initialization
%               1 - network interface computation
%               2 - generator dynamics computation
%        v_sbus - svc bus voltage
%
% Output: bus_new - on initialization bus_new is bus matrix with the
%                   reactive generation at the svc buses set to zero;
%                   otherwise, bus_new = bus

%-----------------------------------------------------------------------------%
% Version history
%
% Version:  1.2
% Date:     May 1998
% Author:   Graham Rogers
% Purpose:  Add lead lag element
%
% Version:  1.1
% Date:     July 1995
% Author:   Graham Rogers
% Purpose:  Add vectorization and correct bugs
%
% Version:  1.0
% Author:   Joe H. Chow
% Date:     June 1991
%-----------------------------------------------------------------------------%

global g;  % declaring struct of global variables

bus_new = bus;

if ~isempty(g.svc.svc_con)
    if (flag == 0)  % initialization
        if (i ~= 0)
            % svc_pot(i,1) -- B_cv max on system base
            % svc_pot(i,2) -- B_cv min on system base
            g.svc.svc_pot(i,1) = g.svc.svc_con(i,4)*g.svc.svc_con(i,3) ...
                                 /g.sys.basmva;
            g.svc.svc_pot(i,2) = g.svc.svc_con(i,5)*g.svc.svc_con(i,3) ...
                                 /g.sys.basmva;

            % initial B_cv based on generation
            j = g.bus.bus_int(g.svc.svc_con(i,2));            % bus number
            g.svc.B_cv(i,1) = bus(j,5)/(bus(j,2).*bus(j,2));
            bus_new(j,5) = 0;

            if (g.svc.B_cv(i,1) > g.svc.svc_pot(i,1))
                estr = '\nsvc: B_cv exceeds maximum at initialization at ';
                estr = [estr, 'svc index %0.0f.'];
                error(sprintf(estr,i));
            end
            %
            if (g.svc.B_cv(i,1) < g.svc.svc_pot(i,2))
                estr = '\nsvc: B_cv below minimum at initialization at ';
                estr = [estr, 'svc index %0.0f.'];
                error(sprintf(estr,i));
            end

            % store initial value of B_cv
            g.svc.svc_pot(i,3) = g.svc.B_cv(i,1);
            g.svc.svc_pot(i,4) = bus(j,2) + g.svc.B_cv(i,1)/g.svc.svc_con(i,6);

            % reference voltage
            if (g.svc.svc_con(1,9) ~= 0)
                g.svc.svc_pot(i,5) = g.svc.svc_con(i,8)/g.svc.svc_con(i,9);
            else
                g.svc.svc_pot(i,5) = 1;
            end

            g.svc.B_con(i,1) = g.svc.B_cv(i,1)*(1-g.svc.svc_pot(i,5)) ...
                               /g.svc.svc_con(i,6);
        else  % vectorized calculation
            % svc_pot(:,1) -- B_cv max on system base
            % svc_pot(:,2) -- B_cv min on system base
            g.svc.svc_pot(:,1) = g.svc.svc_con(:,4).*g.svc.svc_con(:,3) ...
                                 /g.sys.basmva;
            g.svc.svc_pot(:,2) = g.svc.svc_con(:,5).*g.svc.svc_con(:,3) ...
                                 /g.sys.basmva;

            % initial B_cv based on generation
            jsvc = g.bus.bus_int(g.svc.svc_con(:,2));                   % bus number
            g.svc.B_cv(:,1) = bus(jsvc,5)./(bus(jsvc,2).*bus(jsvc,2));
            bus_new(jsvc,5) = zeros(g.svc.n_svc,1);

            mask = (g.svc.B_cv(:,1) > g.svc.svc_pot(:,1));
            if any(mask)
                estr = '\nsvc: B_cv exceeds maximum at initialization at ';
                estr = [estr, 'svc index %0.0f.'];
                error(sprintf(estr,find(mask)));
            end

            mask = (g.svc.B_cv(:,1) < g.svc.svc_pot(:,2));
            if any(mask)
                estr = '\nsvc: B_cv below minimum at initialization at ';
                estr = [estr, 'svc index %0.0f.'];
                error(sprintf(estr,find(mask)));
            end

            % store initial value of B_cv
            g.svc.svc_pot(:,3) = g.svc.B_cv(:,1);

            % reference voltage
            g.svc.svc_pot(:,4) = bus(jsvc,2) + g.svc.B_cv(:,1)./g.svc.svc_con(:,6);
            g.svc.svc_pot(:,5) = ones(g.svc.n_svc,1);
            if ~isempty(g.svc.svcll_idx)
                g.svc.svc_pot(g.svc.svcll_idx,5) = ...
                    g.svc.svc_con(g.svc.svcll_idx,8) ...
                    ./g.svc.svc_con(g.svc.svcll_idx,9);
            end

            g.svc.B_con(:,1) = ...
                g.svc.B_cv(:,1).*(ones(g.svc.n_svc,1) ...
                                  - g.svc.svc_pot(:,5))./g.svc.svc_con(:,6);
        end
    end

    if (flag == 1) % network interface computation
                   % no interface calculation required - done in nc_load
    end

    if (flag == 2) % exciter dynamics calculation
                   % for linearization with operating condition at limits,
                   % additional code will be needed
        if (i ~= 0)
            err =  g.svc.svc_sig(i,k) + g.svc.svc_pot(i,4) ...
                   + g.svc.svc_dsig(i,k) - v_sbus;

            if (g.svc.svc_con(i,9) ~= 0)
                g.svc.dB_con(i,k) = ...
                    (-g.svc.B_con(i,k) + err*(1 - g.svc.svc_pot(i,5))) ...
                    /g.svc.svc_con(i,9);
            else
                g.svc.dB_con(i,k) = 0;
            end

            g.svc.dB_cv(i,k) = ...
                (-g.svc.B_cv(i,k) ...
                 + g.svc.svc_con(i,6)*(err*g.svc.svc_pot(i,5) ...
                                       + g.svc.B_con(i,k)))/g.svc.svc_con(i,7);

            % anti-windup reset
            if (g.svc.B_cv(i,k) > g.svc.svc_pot(i,1))
                if (g.svc.dB_cv(i,k) > 0)
                    g.svc.dB_cv(i,k) = 0;
                end
            end
            %
            if (g.svc.B_cv(i,k) < g.svc.svc_pot(i,2))
                if (g.svc.dB_cv(i,k) < 0)
                    g.svc.dB_cv(i,k) = 0;
                end
            end
        else  % vectorized computation
            lv_sbus = find(v_sbus < 0.9 & g.svc.svc_dsig(:,k) < 0);
            d_sigin = g.svc.svc_dsig(:,k);
            if ~isempty(lv_sbus)
                d_sigin(lv_sbus) = zeros(length(lv_sbus),1);
            end

            err = g.svc.svc_sig(:,k) + g.svc.svc_pot(:,4) + d_sigin - v_sbus;

            g.svc.dB_con(:,k) = zeros(g.svc.n_svc,1);
            if ~isempty(g.svc.svcll_idx)
                nll = length(g.svc.svcll_idx);
                g.svc.dB_con(g.svc.svcll_idx,k) = ...
                    (-g.svc.B_con(g.svc.svcll_idx,k) ...
                     + (ones(nll,1) - g.svc.svc_pot(g.svc.svcll_idx,5)).*err) ...
                    ./g.svc.svc_con(g.svc.svcll_idx,9);
            end

            g.svc.dB_cv(:,k) = ...
                (-g.svc.B_cv(:,k) ...
                 + g.svc.svc_con(:,6).*(err.*g.svc.svc_pot(:,5) ...
                                        + g.svc.B_con(:,k)))./g.svc.svc_con(:,7);

            % anti-windup reset
            indmx = find(g.svc.B_cv(:,k) > g.svc.svc_pot(:,1));
            if ~isempty(indmx)
                g.svc.B_cv(indmx,k) = g.svc.svc_pot(indmx,1);
                indrate = find(g.svc.dB_cv(indmx,k) > 0);
                if ~isempty(indrate)
                    % set rate to zero
                    g.svc.dB_cv(indmx(indrate),k) = zeros(length(indrate),1);
                end
            end

            indmn = find(g.svc.B_cv(:,k) < g.svc.svc_pot(:,2));
            if ~isempty(indmn)
                g.svc.B_cv(indmn,k) = g.svc.svc_pot(indmn,2);
                indrate = find(g.svc.dB_cv(indmn)<0);
                if ~isempty(indrate)
                    % set rate to zero
                    g.svc.dB_cv(indmn(indrate),k) = zeros(length(indrate),1);
                end
            end
        end
    end
end

end  % function end

% eof
