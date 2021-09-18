function dpwf(i,k,flag)
% Syntax: dpwf(i,k,flag)
%
% Purpose: filter model for dpw stabilizer
%
% Input: i - generator number
%        k - integer time
%        flag - 0 - initialization
%               1 - network interface computation
%               2 - system dynamics computation
%
% dpw_con - input data format (5-stage filter assumed)
% col1: generator number
% col2: number of lead elements (maximum 4)
% col3: lead time constant
% col4: lag time constant

%-----------------------------------------------------------------------------%
% Version history
%
% Version: 1.0
% Date:    July 1998
% Author:  Graham Rogers
%-----------------------------------------------------------------------------%

global g;  % declaring struct of global variables

if (g.dpw.n_dpw ~= 0)
    if (flag == 0)   % initialization
        if (i ~= 0)  % scalar computation
            error('dpwf: initialization must be vectorized.');
        else
            % vector computation
            g.dpw.dpw_pot = ones(g.dpw.n_dpw,2);
            n = g.dpw.dpw_mb_idx(g.dpw.dpw_pss_idx);  % generator numbers

            if isempty(n)
                estr = 'dpwf: there are no generators associated with the ';
                estr = [estr, 'power system stabilizers.'];
                error(estr);
            end

            g.dpw.dpw_out(g.dpw.dpw_pss_idx,1) = zeros(g.dpw.n_dpw,1);
            g.dpw.dpw_pot(:,1) = g.dpw.dpw_con(:,3)./g.dpw.dpw_con(:,4);
            g.dpw.dpw_pot(:,2) = 1/2./g.mac.mac_con(n,16);          % 1/2/H
            g.dpw.dpw_pot(:,3) = g.sys.basmva./g.mac.mac_con(n,3);  % pu base ratio
            g.dpw.dpw_pot(:,4) = g.dpw.dpw_pot(:,1).*g.dpw.dpw_Td_idx(:,1);
            g.dpw.dpw_pot(:,5) = g.dpw.dpw_pot(:,1).*g.dpw.dpw_Td_idx(:,2);
            g.dpw.dpw_pot(:,6) = g.dpw.dpw_pot(:,1).*g.dpw.dpw_Td_idx(:,3);
            g.dpw.dpw_pot(:,7) = g.dpw.dpw_pot(:,1).*g.dpw.dpw_Td_idx(:,4);

            % initialize first state
            g.dpw.sdpw1(:,1) = 10*g.mac.pelect(n,1) ...
                               .*g.dpw.dpw_pot(:,3).*g.dpw.dpw_pot(:,2);

            % initialize all filter states
            g.dpw.sdpw2(:,1) = g.dpw.sdpw1(:,1) ...
                               .*(ones(g.dpw.n_dpw,1) - g.dpw.dpw_pot(:,4));

            var1 = g.dpw.sdpw1(:,1).*g.dpw.dpw_pot(:,4);
            g.dpw.sdpw3(:,1) = (var1 + g.dpw.sdpw2(:,1)) ...
                               .*(ones(g.dpw.n_dpw,1) - g.dpw.dpw_pot(:,5));

            var2 = (var1 + g.dpw.sdpw2(:,1)).*g.dpw.dpw_pot(:,5);
            g.dpw.sdpw4(:,1) = (var2 + g.dpw.sdpw3(:,1)) ...
                               .*(ones(g.dpw.n_dpw,1) - g.dpw.dpw_pot(:,6));

            var3 = (var2 + g.dpw.sdpw3(:,1)).*g.dpw.dpw_pot(:,6);
            g.dpw.sdpw5(:,1) = (var3 + g.dpw.sdpw4(:,1)) ...
                               .*(ones(g.dpw.n_dpw,1) - g.dpw.dpw_pot(:,7));

            var4 = (var3 + g.dpw.sdpw4(:,1)).*g.dpw.dpw_pot(:,6);
            g.dpw.sdpw6(:,1) = var4 + g.dpw.sdpw5(:,1);

            g.dpw.dpw_out(g.dpw.dpw_pss_idx,1) = g.dpw.sdpw6(:,1) - g.dpw.sdpw1(:,1);
        end
    end

    if (flag == 1)   % network interface computation
                     % no interface calculation is required
    end

    if (flag == 2)   % dpw filter dynamics calculation
        if (i ~= 0)  % scalar computation
            error('dpwf: dynamics calculation must be vectorized.');
        else
            % vector computation
            if (g.dpw.n_dpw ~= 0)
                n = g.dpw.dpw_mb_idx;  % machine number vector
                var1 = g.mac.mac_spd(n,k) - ones(g.dpw.n_dpw,1) + g.dpw.sdpw1(:,k);
                var2 = var1.*g.dpw.dpw_pot(:,4) + g.dpw.sdpw2(:,k);
                var3 = var2*g.dpw.dpw_pot(:,5) + g.dpw.sdpw3(:,k);
                var4 = var3*g.dpw.dpw_pot(:,6) + g.dpw.sdpw4(:,k);
                var5 = var4*g.dpw.dpw_pot(:,7) + g.dpw.sdpw5(:,k);

                % integrator for power input 10/(1+10s)
                g.dpw.dsdpw1(:,k) = -g.dpw.sdpw1(:,k)/10 ...
                                    + g.mac.pelect(n,k).*g.dpw.dpw_pot(:,3) ...
                                      .*g.dpw.dpw_pot(:,2);

                % filter state rate of change
                g.dpw.dsdpw2(:,k) = ...
                    (var1.*(ones(g.dpw.n_dpw,1) - g.dpw.dpw_pot(:,4)) ...
                     - g.dpw.sdpw2(:,k))./g.dpw.dpw_con(:,4);

                g.dpw.dsdpw3(:,k) = ...
                    (var2.*(ones(g.dpw.n_dpw,1) - g.dpw.dpw_pot(:,5)) ...
                     - g.dpw.sdpw3(:,k))./g.dpw.dpw_con(:,4);

                g.dpw.dsdpw4(:,k) = ...
                    (var3.*(ones(g.dpw.n_dpw,1) - g.dpw.dpw_pot(:,6)) ...
                     - g.dpw.sdpw4(:,k))./g.dpw.dpw_con(:,4);

                g.dpw.dsdpw5(:,k) = ...
                    (var4.*(ones(g.dpw.n_dpw,1) - g.dpw.dpw_pot(:,7)) ...
                     - g.dpw.sdpw5(:,k))./g.dpw.dpw_con(:,4);

                g.dpw.dsdpw6(:,k) = (var5 - g.dpw.sdpw6(:,k))./g.dpw.dpw_con(:,4);

                g.dpw.dpw_out(g.dpw.dpw_pss_idx,k) = g.dpw.sdpw6(:,k) ...
                                                     - g.dpw.sdpw1(:,k);
            end
        end
    end
end

end  % function end

% eof
