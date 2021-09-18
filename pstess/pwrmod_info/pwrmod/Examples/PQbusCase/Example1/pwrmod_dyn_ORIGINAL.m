function [P,Q,dP_states,dQ_states,P_statesIni,Q_statesIni] = pwrmod_dyn(P_states,Q_states,bus,Time,kSim,Flag,n_pwrmod)
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
%   bus = initially solved bus matrix.
%   
% Outputs:
%   P_statesIni = = cell array of initial of P_states
%   Q_statesIni = = cell array of initial of Q_states
%   dP_state = cell array of d/dt of P_states.
%   dQ_state = cell array of d/dt of Q_states.
%   P = n_pwrmod by 1 column vector of P injections at t = Time(kSim).
%   Q = n_pwrmod by 1 column vector of Q injections at t = Time(kSim).
%   Note that injections are either power or current depending on load_con.
%
% Global:
%   pwrmod_data = general variable for storing data when necessary.
%   bus_v = solved bus voltages.  bus_v(n,k) is the solved voltage at bus n
%       at time t(k).  Note that for k>kSim, bus_v is not solved. 
% D. Trudnowski, 2015

global pwrmod_data bus_v pwrmod_con

%% Parameters
nOrderP = 1; %order of state equations for P modulation
nOrderQ = 1; %order of state equations for Q modulation

%% Initialize output variables
%% Initialize output variables
P = zeros(n_pwrmod,1);
Q = zeros(n_pwrmod,1);
dP_states = cell(n_pwrmod,1);
dQ_states = cell(n_pwrmod,1);
P_statesIni = cell(n_pwrmod,1);
Q_statesIni = cell(n_pwrmod,1);

%% Define and initialize state derivatives at t = 0 (set to zero).
if Flag==0
    for k=1:n_pwrmod
        Q_statesIni{k} = zeros(nOrderQ,1);
        P_statesIni{k} = zeros(nOrderP,1);
    end
    pwrmod_data = zeros(length(Time),2); %Store Pref in pwrmod_data
    clear k

%% Calculate P and Q
elseif Flag==1
    %Not used

%% Calculate derivatives
elseif Flag==2
    %Not used
    
end


end

