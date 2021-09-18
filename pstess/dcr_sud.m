function [y,xn,dx] = dcr_sud(i,k,flag,s,d_sig,ymax,ymin,x)
% Syntax: [y,xn,dx] = dcr_sud(i,k,flag,s,d_sig,x)
%
% Purpose: rectifier user-defined damping control
%
% Input:   i - converter number
%          k - integer time
%          flag - 0 - initialization
%                 1 - network interface computation
%                 2 - system dynamics computation
%          s - control state space object
%          d_sig - damping control input
%          x - state
%
% Output:  y - controller output
%         xn - state with limits applied
%         dx - rate of change of state
%
% See Also: dc_cont, dci_sud

%-----------------------------------------------------------------------------%
% Version history
%
% Author: Graham Rogers
% Date:   June 1999
%-----------------------------------------------------------------------------%

y = 0;

if (flag == 0)   % initialization
    if (i ~= 0)  % scalar computation
        % check user-defined control
        if (s.NumInputs ~= 1)
            error('dcr_sud: ud rectifier stabilizer must be single-input.');
        end

        if (s.NumOutputs ~= 1)
            error('dcr_sud: ud rectifier stabilizer must be single-output.');
        end

        [y,xn] = eval(s,0,d_sig);
        if (y > ymax)
            warning('dcr_sud: y outside maximum limit initially.');
            y = ymax;
        end

        if (y < ymin)
            warning('dcr_sud: y outside minimum limit initially.');
            y = ymin;
        end

        dx = zeros(size(xn));
    else
        error('dcr_sud: initialization must be vectorized.');
    end
end

if (flag == 1)   % network interface computation
    if (i ~= 0)  % scalar computation
        ns = s.NumStates;
        xmax = 1e5*ones(ns,1);
        xmin = -xmax;
        y = dstate(s,x,d_sig,xmax,xmin,ymax,ymin);
    else
        error('dcr_sud: network interface calculation must be vectorized.');
    end
end

if (flag == 2)   % rectifier damping control dynamics calculation
    if (i ~= 0)  % scalar computation
        ns = s.NumStates;
        xmax = 1e5*ones(ns,1);
        xmin = -xmax;
        [y,xn,dx] = dstate(s,x,d_sig,xmax,xmin,ymax,ymin);
    else
        error('dcr_sud: dynamics calculation must be vectorized.');
    end
end

end  % function end

% eof
