function smppi(i,k,flag)
% Syntax: smppi(i,k,flag)
%
% Purpose: simple excitation system with pi avr, (exc_con(i,1) = 4)
%          with vectorized computation option
%
% Input:   i - generator number
%              if 0 - vectorized computation
%          k - integer time
%          flag - 0 - initialization
%                 1 - network interface computation
%                 2 - generator dynamics computation
%
% See Also: smpexc, exc_dc12, exc_st3

%-----------------------------------------------------------------------------%
% Version history
%
% Version: 1.0
% Date:    September 1999
% Author:  Graham Rogers
%-----------------------------------------------------------------------------%

global g;  % declaring struct of global variables

if (i ~= 0)
    if (g.exc.exc_con(i,1) ~= 4)
        estr = '\nsmppi: the specified exciter type does not match the ';
        estr = [estr, 'dynamic model at exciter index %0.0f.'];
        error(sprintf(estr,i));
    end
end

[nexc,~] = size(g.exc.exc_con);

if (flag == 0)   % initialization
    % n           -- index of machine(s) with simple PI exciters
    % err         -- summing junction error
    % Efd         -- exciter output voltage (same as generator field voltage)
    % V_As        -- integrator state
    % V_TR        -- input filter state
    % exc_pot(,3) -- reference voltage
    if (i ~= 0)  % scalar computation
        n = g.mac.mac_int(g.exc.exc_con(i,2));
        g.exc.Efd(i,1) = g.mac.vex(n,1);
        g.exc.V_As(i,1) = g.exc.Efd(i,1);
        g.exc.V_TR(i,1) = g.mac.eterm(n,1);
        err = 0;
        g.exc.exc_pot(i,3) = g.mac.eterm(n,1);
    else  % vectorized computation
        if (g.exc.n_smppi ~= 0)
            n = g.mac.mac_int(g.exc.exc_con(g.exc.smppi_idx,2));
            g.exc.Efd(g.exc.smppi_idx,1) = g.mac.vex(n,1);
            g.exc.V_As(g.exc.smppi_idx,1) = g.exc.Efd(g.exc.smppi_idx,1);
            g.exc.V_TR(g.exc.smppi_idx,1) = g.mac.eterm(n,1);
            g.exc.exc_pot(g.exc.smppi_idx,3) = g.mac.eterm(n,1);
        end
    end
end

if (flag == 1)   % network interface computation
    if (i ~= 0)  % scalar computation
        n = g.mac.mac_int(g.exc.exc_con(i,2));
        g.mac.vex(n,k) = g.exc.Efd(i,k);        % field voltage on machine base
    else  % vectorized computation
        if (g.exc.n_smppi ~= 0)                 % check for simple PI exciters
            n = g.mac.mac_int(g.exc.exc_con(g.exc.smppi_idx,2));
            g.mac.vex(n,k) = g.exc.Efd(g.exc.smppi_idx,k);
        end
    end
end

if (flag == 2)   % exciter dynamics calculation
    if (i ~= 0)  % scalar computation
        n = g.mac.mac_int(g.exc.exc_con(i,2));  % machine number
        if (g.exc.exc_con(i,3) == 0)            % no input filter
            g.exc.dV_TR(i,k) = 0;
            g.exc.V_TR(i,k) = g.mac.eterm(n,k);
        else
            g.exc.dV_TR(i,k) = (g.mac.eterm(n,k) - g.exc.V_TR(i,k)) ...
                               /g.exc.exc_con(i,3);
        end

        err = g.exc.exc_sig(i,k) + g.exc.exc_pot(i,3) - g.exc.V_TR(i,k) ...
              + g.exc.pss_out(i,k);

        g.exc.dV_As(i,k) = err*g.exc.exc_con(i,4);

        g.exc.dEfd(i,k) = (-g.exc.Efd(i,k) + g.exc.V_A(i,k) ...
                           + g.exc.exc_con(i,6)*err)/g.exc.exc_con(i,5);

        % anti-windup reset
        if (g.exc.Efd(i,k) > g.exc.exc_con(i,8))
            g.exc.Efd(i,k) = g.exc.exc_con(i,8);
            if (g.exc.dEfd(i,k) > 0)
                g.exc.dEfd(i,k) = 0;
            end
        end

        if (g.exc.Efd(i,k) < g.exc.exc_con(i,9))
            g.exc.Efd(i,k) = g.exc.exc_con(i,9);
            if (g.exc.dEfd(i,k) < 0)
                g.exc.dEfd(i,k) = 0;
            end
        end

        g.exc.R_f(i,k) = 0;
        g.exc.dR_f(i,k) = 0;
        %
        g.exc.V_R(i,k) = 0;
        g.exc.dV_R(i,k) = 0;

    else % vectorized computation

        if (g.exc.n_smppi ~= 0)
            % machine numbers
            n = g.mac.mac_int(g.exc.exc_con(g.exc.smppi_idx,2));
            TR = g.exc.smppi_TR_idx;
            no_TR = g.exc.smppi_noTR_idx;

            % for exciters with zero TR
            if ~isempty(no_TR)
                n_nTR = n(no_TR);
                g.exc.dV_TR(g.exc.smppi_idx(no_TR),k) = zeros(length(no_TR),1);
                g.exc.V_TR(g.exc.smppi_idx(no_TR),k) = g.mac.eterm(n_nTR,k);
            end

            % for exciters with nonzero TR
            if ~isempty(TR)
                n_TR = n(TR);
                g.exc.dV_TR(g.exc.smppi_idx(TR),k) = ...
                    (g.mac.eterm(n_TR,k) - g.exc.V_TR(g.exc.smppi_idx(TR),k)) ...
                    ./g.exc.exc_con(g.exc.smppi_idx(TR),3);
            end

            % error defined for all simple exciters
            err = g.exc.exc_sig(g.exc.smppi_idx,k) ...
                  + g.exc.exc_pot(g.exc.smppi_idx,3) ...
                  - g.exc.V_TR(g.exc.smppi_idx,k) ...
                  + g.exc.pss_out(g.exc.smppi_idx,k);

            g.exc.dV_As(g.exc.smppi_idx,k) = err*g.exc.exc_con(g.exc.smppi_idx,4);

            g.exc.dEfd(g.exc.smppi_idx,k) = ...
                (-g.exc.Efd(g.exc.smppi_idx,k) + g.exc.V_As(g.exc.smppi_idx,k) ...
                 + g.exc.exc_con(g.exc.smppi_idx,6)*err) ...
                ./g.exc.exc_con(g.exc.smppi_idx,5);

            TA_max = find(g.exc.Efd(g.exc.smppi_idx,k) ...
                          > g.exc.exc_con(g.exc.smppi_idx,8));
            TA_min = find(g.exc.Efd(g.exc.smppi_idx,k) ...
                          < g.exc.exc_con(g.exc.smppi_idx,9));

            if ~isempty(TA_max)
                g.exc.Efd(g.exc.smppi_idx(TA_max),k) = ...
                    g.exc.exc_con(g.exc.smppi_idx(TA_max),8);

                dTA_max = find(g.exc.dEfd(g.exc.smppi_idx(TA_max),k) > 0);
                if ~isempty(dTA_max)
                    n_dTA = length(dTA_max);
                    dEFD(g.exc.smppi_idx(TA_max(dTA_max)),k) = zeros(n_dTA,1);
                end
            end

            if ~isempty(TA_min)
                g.exc.Efd(g.exc.smppi_idx(TA_min),k) = ...
                    g.exc.exc_con(g.exc.smppi_idx(TA_min),9);

                dTA_min = find(g.exc.dEfd(g.exc.smppi_idx(TA_min),k) < 0);
                if ~isempty(dTA_min)
                    n_dTA = length(dTA_min);
                    dEFD(g.exc.smppi_idx(TA_min(dTA_min)),k) = zeros(n_dTA,1);
                end
            end

            g.exc.R_f(g.exc.smppi_idx,k) = zeros(g.exc.n_smppi,1);
            g.exc.dR_f(g.exc.smppi_idx,k) = zeros(g.exc.n_smppi,1);
            %
            g.exc.V_R(g.exc.smppi_idx,k) = zeros(g.exc.n_smppi,1);
            g.exc.dV_R(g.exc.smppi_idx,k) = zeros(g.exc.n_smppi,1);
        end
    end
end

end  % function end

% eof
