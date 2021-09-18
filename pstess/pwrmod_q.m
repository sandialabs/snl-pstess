function pwrmod_q(i,k,bus,flag)
% Syntax: pwrmod_q(i,k,bus,flag)
%
% Purpose: reac-power modulation model
%          with vectorized computation option
%          NOTE - load modulation bus must be declared as a
%                 non-conforming constant-power load bus
%
% Input:   i - power modulation number
%              if i= 0, vectorized computation
%          k - integer time
%          bus - solved loadflow bus data
%          flag - 0 - initialization
%                 1 - network interface computation
%                 2 - generator dynamics computation

%-----------------------------------------------------------------------------%
% Version history
%
% Version: 1.0 (initial version)
% Date:    Feb 2015
% Author:  Dan Trudnowski
%-----------------------------------------------------------------------------%

global g;  % declaring struct of global variables

if ~isempty(g.pwr.pwrmod_con)
    if (flag == 0)  % initialization
        if (i ~= 0)
            jj = g.bus.bus_int(g.pwr.pwrmod_con(i,1));
            if ((g.ncl.load_con(g.pwr.pwrmod_idx(i),2) == 1) ...
                && (g.ncl.load_con(g.pwr.pwrmod_idx(i),3) == 1))
                g.pwr.pwrmod_q_st(i,1) = bus(jj,5);
            else
                g.pwr.pwrmod_q_st(i,1) = bus(jj,5)/abs(bus(jj,2));
            end
            clear jj
        else  % vectorized calculation
            for ii = 1:g.pwr.n_pwrmod
                jj = g.bus.bus_int(g.pwr.pwrmod_con(ii,1));
                if ((g.ncl.load_con(g.pwr.pwrmod_idx(ii),2) == 1) ...
                    && (g.ncl.load_con(g.pwr.pwrmod_idx(ii),3) == 1))
                    g.pwr.pwrmod_q_st(ii,1) = bus(jj,5);
                else
                    g.pwr.pwrmod_q_st(ii,1) = bus(jj,5)/abs(bus(jj,2));
                end
            end
            clear jj
        end
    end

    if (flag == 1)  % network interface computation
        % no interface calculation required - done in nc_load
    end

    if (flag == 2)  % dynamics calculation
        % for linearization with operating condition at limits,
        % additional code will be needed
        if (i ~= 0)
            g.pwr.dpwrmod_q_st(i,k) = ...
                (-g.pwr.pwrmod_q_st(i,k) + g.pwr.pwrmod_q_sig(i,k)) ...
                /g.pwr.pwrmod_con(i,5);

            % anti-windup reset
            if (g.pwr.pwrmod_q_st(i,k) > g.pwr.pwrmod_con(i,6))
                if (g.pwr.dpwrmod_q_st(i,k) > 0)
                    g.pwr.dpwrmod_q_st(i,k) = 0;
                end
            end

            if (g.pwr.pwrmod_q_st(i,k) < g.pwr.pwrmod_con(i,7))
                if (g.pwr.dpwrmod_q_st(i,k) < 0)
                    g.pwr.dpwrmod_q_st(i,k) = 0;
                end
            end
        else  % vectorized computation
            g.pwr.dpwrmod_q_st(:,k) = ...
                (-g.pwr.pwrmod_q_st(:,k) + g.pwr.pwrmod_q_sig(:,k)) ...
                ./g.pwr.pwrmod_con(:,5);

            % anti-windup reset
            indmx = find(g.pwr.pwrmod_q_st(:,k) > g.pwr.pwrmod_con(:,6));
            if ~isempty(indmx)
                indrate = find(g.pwr.dpwrmod_q_st(indmx,k) > 0);
                if ~isempty(indrate)
                    % set rate to zero
                    g.pwr.dpwrmod_q_st(indmx(indrate),k) = zeros(length(indrate),1);
                end
            end

            indmn = find(g.pwr.pwrmod_q_st(:,k) < g.pwr.pwrmod_con(:,7));
            if ~isempty(indmn)
                indrate = find(g.pwr.dpwrmod_q_st(indmn) < 0);
                if ~isempty(indrate)
                    % set rate to zero
                    g.pwr.dpwrmod_q_st(indmn(indrate),k) = zeros(length(indrate),1);
                end
            end
        end
    end
end

end  % function end

% eof
