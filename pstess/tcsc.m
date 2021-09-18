function tcsc(i,k,flag)
% Syntax: tcsc(i,k,flag)
%
% Purpose: thyristor controlled series capacitor,
%          with vectorized computation option
%
%          NOTE - TCSC buses must be declared as non-conforming load buses
%                 The initial capacitance must be inserted as a loss free
%                 negative reactance between the two tcsc buses.
%                 These buses are inserted in the compensated line
%
% Input:   i - tcsc number: if i = 0, vectorized computation
%          k - integer time
%          flag - 0 - initialization
%                 1 - network interface computation
%                 2 - generator dynamics computation

%-----------------------------------------------------------------------------%
% Version history
%
% Version: 1.0
% Date:    December 1998
% Author:  Graham Rogers
%-----------------------------------------------------------------------------%

global g;  % declaring struct of global variables

if ~isempty(g.tcsc.tcsc_con)
    if (flag == 0)  % initialization
        if (i ~= 0)
            g.tcsc.B_tcsc(i,1) = 0;
            g.tcsc.dB_tcsc(i,1) = 0;
        else  % vectorized calculation
            g.tcsc.B_tcsc(:,1) = zeros(g.tcsc.n_tcsc,1);
            g.tcsc.dB_tcsc(:,1) = zeros(g.tcsc.n_tcsc,1);
        end
    end
    if (flag == 1)  % network interface computation
        % no interface calculation required - done in nc_load
    end
    if (flag == 2)  % tcsc dynamics calculation
        if (i ~= 0)
            err = g.tcsc.tcsc_sig(i,k) + g.tcsc.tcsc_dsig(i,k);
            g.tcsc.dB_tcsc(i,k) = (-g.tcsc.B_tcsc(i,k) ...
                                   + g.tcsc.tcsc_con(i,4)*err) ...
                                  /g.tcsc.tcsc_con(i,5);

            % anti-windup reset
            if (g.tcsc.B_tcsc(i,k) > g.tcsc.tcsc_con(i,6))
                g.tcsc.B_tcsc(i,k) = g.tcsc.tcsc_con(i,6);
                if (g.tcsc.dB_tcsc(i,k) > 0)
                    g.tcsc.dB_tcsc(i,k) = 0;
                end
            end
            if (g.tcsc.B_tcsc(i,k) < g.tcsc.tcsc_con(i,7))
                g.tcsc.B_tcsc(i,k) = g.tcsc.tcsc_con(i,7);
                if (g.tcsc.dB_tcsc(i,k) < 0)
                    g.tcsc.dB_tcsc(i,k) = 0;
                end
            end
        else  % vectorized computation
            err = g.tcsc.tcsc_sig(:,k) + g.tcsc.tcsc_dsig(:,k);
            g.tcsc.dB_tcsc(:,k) = (-g.tcsc.B_tcsc(:,k) ...
                                   + g.tcsc.tcsc_con(:,4).*err) ...
                                  ./g.tcsc.tcsc_con(:,5);

            % anti-windup reset
            indmx = find(g.tcsc.B_tcsc(:,k) > g.tcsc.tcsc_con(:,6));
            if ~isempty(indmx)
                g.tcsc.B_tcsc(indmx,k) = g.tcsc.tcsc_con(indmx,6);
                indrate = find(g.tcsc.dB_tcsc(indmx,k)>0);
                if ~isempty(indrate)
                    % set rate to zero
                    g.tcsc.dB_tcsc(indmx(indrate),k) = zeros(length(indrate),1);
                end
            end
            indmn = find(g.tcsc.B_tcsc(:,k) < g.tcsc.tcsc_con(:,7));
            if ~isempty(indmn)
                g.tcsc.B_tcsc(indmn,k) = g.tcsc.tcsc_con(indmn,7);
                indrate = find(g.tcsc.dB_tcsc(indmn,k) < 0);
                if ~isempty(indrate)
                    % set rate to zero
                    g.tcsc.dB_tcsc(indmn(indrate),k) = zeros(length(indrate),1);
                end
            end
        end
    end
end

end  % function end

% eof
