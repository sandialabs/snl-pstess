function lmod(i,k,flag)
% Syntax: lmod(i,k,flag)
%
% Purpose: load modulation model
%          with vectorized computation option
%          NOTE - load modulation bus must be declared as a
%                 non-conforming load bus
%
% Input: i - load modulation number
%            if i= 0, vectorized computation
%        k - integer time
%        flag - 0 - initialization
%               1 - network interface computation
%               2 - generator dynamics computation

%-----------------------------------------------------------------------------%
% Version history
%
% Version: 1.0 (initial version)
% Date:    August 1997
% Author:  Graham Rogers
%-----------------------------------------------------------------------------%

global g;  % declaring struct of global variables

if ~isempty(g.lmod.lmod_con)
    if (flag == 0)  % initialization
        if (i ~= 0)
            % lmod_pot(i,1) -- max modulation on system base
            % lmod_pot(i,2) -- min modulation on system base
            g.lmod.lmod_pot(i,1) = g.lmod.lmod_con(i,4) ...
                                   *g.lmod.lmod_con(i,3)/g.sys.basmva;
            g.lmod.lmod_pot(i,2) = g.lmod.lmod_con(i,5) ...
                                   *g.lmod.lmod_con(i,3)/g.sys.basmva;
            g.lmod.lmod_st(i,1) = 0;
        else  % vectorized calculation
            % lmod_pot(:,1) -- max modulation on system base
            % lmod_pot(:,2) -- min modulation on system base
            g.lmod.lmod_pot(:,1) = g.lmod.lmod_con(:,4) ...
                                   .*g.lmod.lmod_con(:,3)/g.sys.basmva;
            g.lmod.lmod_pot(:,2) = g.lmod.lmod_con(:,5) ...
                                   .*g.lmod.lmod_con(:,3)/g.sys.basmva;
            g.lmod.lmod_st(:,1) = zeros(g.lmod.n_lmod,1);
        end
    end

    if (flag == 1)  % network interface computation
        % no interface calculation required - done in nc_load
    end

    if (flag == 2)  % dynamics calculation
        % for linearization with operating condition at limits,
        % additional code will be needed
        if (i ~= 0)
            u_lmod = g.lmod.lmod_con(i,6)*g.lmod.lmod_sig(i,k);
            g.lmod.dlmod_st(i,k) = (u_lmod - g.lmod.lmod_st(i,k)) ...
                                   /g.lmod.lmod_con(i,7);

            if (g.lmod.lmod_st(i,k) > g.lmod.lmod_pot(i,1))     % anti-windup reset
                g.lmod.lmod_st(i,k) = g.lmod.lmod_pot(i,1);
                if (g.lmod.dlmod_st(i,k) > 0)
                    g.lmod.dlmod_st(i,k) = 0;
                end
            end

            if (g.lmod.lmod_st(i,k) < g.lmod.lmod_pot(i,2))     % anti-windup reset
                g.lmod.lmod_st(i,k) = g.lmod.lmod_pot(i,2);
                if (g.lmod.dlmod_st(i,k) < 0)
                    g.lmod.dlmod_st(i,k) = 0;
                end
            end
        else  % vectorized computation
            u_lmod = g.lmod.lmod_con(:,6).*g.lmod.lmod_sig(:,k);
            g.lmod.dlmod_st(:,k) = (u_lmod - g.lmod.lmod_st(:,k)) ...
                                   ./g.lmod.lmod_con(:,7);

            indmx = find(g.lmod.lmod_st(:,k) > g.lmod.lmod_pot(:,1));
            if ~isempty(indmx)                                  % anti-windup reset
                g.lmod.lmod_st(indmx,k) = g.lmod.lmod_pot(indmx,1);
                indrate = find(g.lmod.dlmod_st(indmx,k) > 0);
                if ~isempty(indrate)
                    g.lmod.dlmod_st(indmx(indrate),k) = zeros(length(indrate),1);
                end
            end

            indmn = find(g.lmod.lmod_st(:,k) < g.lmod.lmod_pot(:,2));
            if ~isempty(indmn)                                  % anti-windup reset
                g.lmod.lmod_st(indmn,k) = g.lmod.lmod_pot(indmn,2);
                indrate = find(g.lmod.dlmod_st(indmn) < 0);
                if ~isempty(indrate)
                    g.lmod.dlmod_st(indmn(indrate),k) = zeros(length(indrate),1);
                end
            end
        end
    end
end

end  % function end

% eof
