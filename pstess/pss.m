function pss(i,k,flag)
% Syntax: pss(i,k,flag)
%
% Purpose: power system stabilization model
%
% Input: i - generator number
%        k - integer time
%        flag - 0 - initialization
%               1 - network interface computation
%               2 - system dynamics computation
%
% See Also: exc_dc12, exc_st3

%-----------------------------------------------------------------------------%
% Version history
%
% Version:  2.2
% Date:     July 1998
% Purpose:  added interface to deltaP/omega filter
% Author:   Graham Rogers
%
% Version:  2.1
% Date:     October 1997
% Purpose:  changed length to isempty in index checks
% Author:   Graham Rogers
%
% Version:  2.0
% Date:     June 1996
% Author:   Graham Rogers
% Purpose:  To allow vector calculation with pss on only some units
% Modified: modified vector code, added negative output limit
%
% Version:  1.0
% Author:   Joe H. Chow
% Date:     April 1992
%-----------------------------------------------------------------------------%

global g;  % declaring struct of global variables

if (g.pss.n_pss ~= 0)
    if (flag == 0)   % initialization
        if (i ~= 0)  % scalar computation
            n = g.pss.pss_mb_idx(i);  % machine number
            if (g.pss.pss_con(i,1) == 1)
                g.pss.pss1(i,1) = g.mac.mac_spd(n,1);
            else
                g.pss.pss1(i,1) = g.mac.pelect(n,1)*g.sys.basmva/g.mac.mac_con(n,3);
            end

            if (g.dpw.n_dpw ~= 0)
                i_dpw = find(g.dpw.dpw_pss_idx == i);
                if ~isempty(i_dpw)
                    g.pss.pss1(i,1) = g.dpw.dpw_out(i_dpw,1);
                end
            end

            g.pss.pss2(i,1) = 0.0;
            g.pss.pss3(i,1) = 0.0;
            g.pss.pss_pot(i,1) = g.pss.pss_con(i,5)/g.pss.pss_con(i,6);

            % pss output uses the exciter table index
            g.exc.pss_out(g.pss.pss_exc_idx(i),1) = 0.0;

            g.pss.pss_pot(i,2) = 1.0;
            if (g.pss.pss_con(i,8) ~= 0)
                g.pss.pss_pot(i,2) = g.pss.pss_con(i,7)/g.pss.pss_con(i,8);
            end

            if (g.pss.pss_con(i,1) ~= 1 & g.pss.pss_con(i,1) ~= 2)
                estr = '\npss: specified type unrecognized for the ';
                estr = [estr, 'stabilizer at pss index %0.0f.'];
                error(sprintf(estr,i));
            end
        else
            % vectorized computation

            g.pss.pss_pot = ones(g.pss.n_pss,2);
            n = g.pss.pss_mb_idx;
            if ~isempty(g.pss.pss_sp_idx)
                n_sp = g.mac.mac_int(g.pss.pss_con(g.pss.pss_sp_idx,2));
                g.pss.pss1(g.pss.pss_sp_idx,1) = g.mac.mac_spd(n_sp,1);
            end

            if ~isempty(g.pss.pss_p_idx)
                n_p = g.mac.mac_int(g.pss.pss_con(g.pss.pss_p_idx,2));
                g.pss.pss1(g.pss.pss_p_idx,1) = g.mac.pelect(n_p,1)*g.sys.basmva ...
                                                ./g.mac.mac_con(n_p,3);
            end

            if (g.dpw.n_dpw ~= 0)
                g.pss.pss1(g.dpw.dpw_pss_idx,1) = g.dpw.dpw_out(:,1);
            end

            g.pss.pss2(g.pss.pss_idx,1) = zeros(g.pss.n_pss,1);
            g.pss.pss3(g.pss.pss_idx,1) = zeros(g.pss.n_pss,1);
            g.pss.pss_pot(:,1) = g.pss.pss_con(:,5)./g.pss.pss_con(:,6);

            % pss output uses the exciter table index
            g.exc.pss_out(g.pss.pss_exc_idx,1) = zeros(g.pss.n_pss,1);

            if ~isempty(g.pss.pss_T4_idx)
                g.pss.pss_pot(g.pss.pss_T4_idx,2) = ...
                    g.pss.pss_con(g.pss.pss_T4_idx,7) ...
                    ./g.pss.pss_T4(g.pss.pss_T4_idx);
            end
        end
    end

    if (flag == 1)   % network interface computation
        if (i ~= 0)  % scalar computation
            n = g.pss.pss_mb_idx(i);  % machine number
            if (g.pss.pss_con(i,1) == 1)
                var1 = (g.mac.mac_spd(i,k) - g.pss.pss1(i,k))/g.pss.pss_con(i,4);
            else
                n = g.mac.mac_int(g.pss.pss_con(i,2));  % machine number
                var1 = (g.mac.pelect(i,k)*g.sys.basmva/g.mac.mac_con(n,3) ...
                        - g.pss.pss1(i,k))/g.pss.pss_con(i,4);
            end

            if (g.dpw.n_dpw ~= 0)
                if (g.dpw.n_dpw ~= 0)
                    i_dpw = find(g.dpw.dpw_pss_idx == i);
                    if ~isempty(i_dpw)
                        var1 = (g.dpw.dpw_out(i_dpw,k) - g.pss.pss1(i,k)) ...
                               /g.pss.pss_con(i,4);
                    end
                end
            end

            var2 = g.pss.pss_pot(i,1)*g.pss.pss_con(i,3)*var1 + g.pss.pss2(i,k);

            if (g.pss.pss_con(i,8) == 0)
                var3 = var2;
            else
                var3 = g.pss.pss_pot(i,2)*var2 + g.pss.pss3(i,k);
            end

            g.exc.pss_out(g.pss.pss_exc_idx(i),k) = ...
                min(g.pss.pss_con(i,9),max(var3,-g.pss.pss_con(i,9)));
        else
            % vectorized computation

            if (g.pss.n_pss ~= 0)
                n = g.pss.pss_mb_idx;  % machine number vector

                var1 = zeros(g.pss.n_pss,1);
                var2 = var1;
                var3 = var1;

                if ~isempty(g.pss.pss_sp_idx)
                    n_sp = g.mac.mac_int(g.pss.pss_con(g.pss.pss_sp_idx,2));
                    var1(g.pss.pss_sp_idx) = ...
                        (g.mac.mac_spd(n_sp,k) - g.pss.pss1(g.pss.pss_sp_idx,k)) ...
                        ./g.pss.pss_con(g.pss.pss_sp_idx,4);
                end

                if ~isempty(g.pss.pss_p_idx)
                    n_p = g.mac.mac_int(g.pss.pss_con(g.pss.pss_p_idx,2));
                    var1(g.pss.pss_p_idx) = ...
                        (g.mac.pelect(n_p,k)*g.sys.basmva./g.mac.mac_con(n_p,3) ...
                         - g.pss.pss1(g.pss.pss_p_idx,k)) ...
                        ./g.pss.pss_con(g.pss.pss_p_idx,4);
                end

                if (g.dpw.n_dpw ~= 0)
                    var1 = (g.dpw.dpw_out(:,k) - g.pss.pss1(g.dpw.dpw_pss_idx,k)) ...
                           ./g.pss.pss_con(g.dpw.dpw_pss_idx,4);
                end
            end

            var2(g.pss.pss_idx) = g.pss.pss_pot(g.pss.pss_idx,1) ...
                                  .*(g.pss.pss_con(g.pss.pss_idx,3).*var1) ...
                                  + g.pss.pss2(g.pss.pss_idx,k);
            var3 = var2;

            if ~isempty(g.pss.pss_T4_idx)
                var3(g.pss.pss_T4_idx,1) = g.pss.pss_pot(g.pss.pss_T4_idx,2) ...
                                           .*var2(g.pss.pss_T4_idx,1) ...
                                           + g.pss.pss3(g.pss.pss_T4_idx,k);
            end

            g.exc.pss_out(g.pss.pss_exc_idx,k) = ...
                min(g.pss.pss_con(g.pss.pss_idx,9), ...
                    max(var3,g.pss.pss_con(g.pss.pss_idx,10)));
        end
    end

    if (flag == 2)   % pss dynamics calculation
        if (i ~= 0)  % scalar computation
            n = g.pss.pss_mb_idx(i);  % machine number
            if (g.pss.pss_con(i,1) == 1)
                var1 = (g.mac.mac_spd(i,k) - g.pss.pss1(i,k))/g.pss.pss_con(i,4);
            else
                n = g.mac.mac_int(g.pss.pss_con(i,2));  % machine number
                var1 = (g.mac.pelect(i,k)*g.sys.basmva./g.mac.mac_con(n,3) ...
                        - g.pss.pss1(i,k))/g.pss.pss_con(i,4);
            end

            if (g.dpw.n_dpw ~= 0)
                if (g.dpw.n_dpw ~= 0)
                    i_dpw = find(g.dpw.dpw_pss_idx == i);
                    if ~isempty(i_dpw)
                        var1 = (g.dpw.dpw_out(i_dpw,k) - g.pss.pss1(i,k)) ...
                               /g.pss.pss_con(i,4);
                    end
                end
            end

            g.pss.dpss1(i,k) = var1;

            var2 = g.pss.pss_pot(i,1)*g.pss.pss_con(i,3)*var1 + g.pss.pss2(i,k);
            g.pss.dpss2(i,k) = ((1 - g.pss.pss_pot(i,1))*g.pss.pss_con(i,3)*var1 ...
                                - g.pss.pss2(i,k))/g.pss.pss_con(i,6);

            if (g.pss.pss_con(i,8) == 0)
                var3 = var2;
                g.pss.dpss3(i,k) = g.pss.dpss2(i,k);
            else
                var3 = g.pss.pss_pot(i,2)*var2 + g.pss.pss3(i,k);
                g.pss.dpss3(i,k) = ((1 - g.pss.pss_pot(i,2))*var2 ...
                                    - g.pss.pss3(i,k))/g.pss.pss_con(i,8);
            end

            g.exc.pss_out(g.pss.pss_exc_idx(i),k) = ...
                min(g.pss.pss_con(i,9),max(var3,-g.pss.pss_con(i,9)));
        else
            % vectorized computation
            if (g.pss.n_pss ~= 0)
                n = g.pss.pss_mb_idx;  % machine number vector

                var1 = zeros(g.pss.n_pss,1);
                var2 = var1;
                var3 = var1;

                if ~isempty(g.pss.pss_sp_idx)
                    n_sp = g.mac.mac_int(g.pss.pss_con(g.pss.pss_sp_idx,2));
                    var1(g.pss.pss_sp_idx) = ...
                        (g.mac.mac_spd(n_sp,k) - g.pss.pss1(g.pss.pss_sp_idx,k)) ...
                        ./g.pss.pss_con(g.pss.pss_sp_idx,4);
                end

                if ~isempty(g.pss.pss_p_idx)
                    n_p = g.mac.mac_int(g.pss.pss_con(g.pss.pss_p_idx,2));
                    var1(g.pss.pss_p_idx) = ...
                        (g.mac.pelect(n_p,k)*g.sys.basmva./g.mac.mac_con(n_p,3) ...
                         - g.pss.pss1(g.pss.pss_p_idx,k)) ...
                        ./g.pss.pss_con(g.pss.pss_p_idx,4);
                end

                if (g.dpw.n_dpw ~= 0)
                    var1 = (g.dpw.dpw_out(:,k) - g.pss.pss1(g.dpw.dpw_pss_idx,k)) ...
                           ./g.pss.pss_con(g.dpw.dpw_pss_idx,4);
                end
            end

            g.pss.dpss1(g.pss.pss_idx,k) = var1;

            var2 = g.pss.pss_pot(g.pss.pss_idx,1) ...
                   .*(g.pss.pss_con(g.pss.pss_idx,3).*var1) ...
                   + g.pss.pss2(g.pss.pss_idx,k);

            g.pss.dpss2(g.pss.pss_idx,k) = ...
                ((ones(g.pss.n_pss,1) - g.pss.pss_pot(g.pss.pss_idx,1)) ...
                 .*(g.pss.pss_con(g.pss.pss_idx,3).*var1) ...
                 - g.pss.pss2(g.pss.pss_idx,k))./g.pss.pss_con(g.pss.pss_idx,6);

            var3 = var2;

            g.pss.dpss3(:,k) = g.pss.dpss2(:,k);
            if ~isempty(g.pss.pss_T4_idx)
                var3(g.pss.pss_T4_idx) = ...
                    g.pss.pss_pot(g.pss.pss_T4_idx,2).*var2(g.pss.pss_T4_idx) ...
                    + g.pss.pss3(g.pss.pss_T4_idx,k);

                g.pss.dpss3(g.pss.pss_T4_idx,k) = ...
                    ((ones(length(g.pss.pss_T4_idx),1) ...
                      - g.pss.pss_pot(g.pss.pss_T4_idx,2)) ...
                     .*var2(g.pss.pss_T4_idx) - g.pss.pss3(g.pss.pss_T4_idx,k)) ...
                    ./g.pss.pss_T4(g.pss.pss_T4_idx);
            end

            g.exc.pss_out(g.pss.pss_exc_idx,k) = ...
                min(g.pss.pss_con(g.pss.pss_idx,9), ...
                    max(var3,g.pss.pss_con(g.pss.pss_idx,10)));
        end
    end
end

end  % function end

% eof
