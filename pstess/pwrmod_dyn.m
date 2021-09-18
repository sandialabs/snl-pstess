function [P,Q,dP_states,dQ_states,P_statesIni,Q_statesIni] = pwrmod_dyn(P_states,Q_states,bus,Time,kSim,Flag)
% Implement state or output variables to model power injection
%
% Inputs:
%   n_pwrmod = number of injection buses (number of rows of pwrmod_con).
%   Time = vector of simulation time
%   kSim = current simulation index.  Current time = Time(kSim).
%   Flag:
%       If Flag==0, Initialize P_statesIni, Q_statesIni at t = 0.
%       If Flag==1, Calculate P, Q at Time(kSim)
%       If Flag==2, Calculate dP_states, dQ_states at Time(kSim)
%   P_state = cell array of P injection states
%       P_state{k} = column vector of states corresponding to row k of
%       pwrmod_con.
%   Q_state = cell array of Q injection states.  Same format as P_states.
%   bus = initial bus matrix from the solved power flow.
%
% Outputs:
%   P_statesIni == cell array of initial of P_states
%   Q_statesIni == cell array of initial of Q_states
%   dP_state = cell array of d/dt of P_states.
%   dQ_state = cell array of d/dt of Q_states.
%   P = n_pwrmod by 1 column vector of P commands at t = Time(kSim).
%       P(k) corresponds to row k of pwrmod_con.
%   Q = n_pwrmod by 1 column vector of Q commands at t = Time(kSim).
%       Q(k) corresponds to row k of pwrmod_con.
%   Note that injections are either power or current depending on load_con.
%
% Global:
%   pwrmod_data = general variable for storing data when necessary.
%   bus_v = solved bus voltages.  bus_v(n,k) is the solved voltage at bus n
%       at time t(k).  Note that for k>kSim, bus_v is not solved.
%   load_con = see system data file.
%   pwrmod_con = see system data file.
% D. Trudnowski, 2015

global g;  % declaring struct of global variables

%% Parameters
nOrderP = [1;1];  % Order of state equations for P modulation
nOrderQ = [1;1];  % Order of state equations for Q modulation

%% Initialize output variables
P = zeros(g.pwr.n_pwrmod,1);
Q = zeros(g.pwr.n_pwrmod,1);
dP_states = cell(g.pwr.n_pwrmod,1);
dQ_states = cell(g.pwr.n_pwrmod,1);
P_statesIni = cell(g.pwr.n_pwrmod,1);
Q_statesIni = cell(g.pwr.n_pwrmod,1);

%% Define and initialize state derivatives at t = 0 (set to zero).
if (Flag == 0)
    for k = 1:g.pwr.n_pwrmod
        Q_statesIni{k} = zeros(nOrderQ(k),1);
        P_statesIni{k} = zeros(nOrderP(k),1);
    end

    g.pwr.pwrmod_data = zeros(length(Time),2);  % Store Pref in pwrmod_data

%% Calculate P and Q
elseif (Flag == 1)
    for k = 1:length(P)
        n = find(g.pwr.pwrmod_con(k,1) == bus(:,1));
        m = find(g.pwr.pwrmod_con(k,1) == g.ncl.load_con(:,1));
        if (g.ncl.load_con(m,2) == 1)
            % Initial power injection
            P(k) = bus(n,4);  % Real power
            Q(k) = bus(n,5);  % Reactive power
        else
            % Initial current injection
            P(k) = bus(n,4)/abs(bus(n,2));  % Real current
            Q(k) = bus(n,5)/abs(bus(n,2));  % Reactive current
        end
    end

    if ((Time(kSim) >= 1) && (Time(kSim) < 1.5))
        P(1) = P(1) - 0.0001;
    end

    if ((Time(kSim) >= 4) && (Time(kSim) < 4.5))
        P(2) = P(2) + 0.0002;
    end

%% Calculate derivatives
elseif (Flag == 2)
    for k = 1:g.pwr.n_pwrmod
        dP_states{k} = zeros(nOrderP(k),1);
        dQ_states{k} = zeros(nOrderQ(k),1);
    end
end

end  % function end

% eof
