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
nOrderP = [2;2]; %order of state equations for P modulation
nOrderQ = [1;1]; %order of state equations for Q modulation

%P modulation parameters
Tpord = [0.25; 0.35]; %Power order filter time constant
Tv = [0.05; 0.15]; %Voltage measurment filter
Ipmax = [0.1;0.1]; %Max modulated P current
Ipmin = [-0.1;-0.1]; %Min modulated P current
dIpmax = [100;100]; %Max dIp
dIpmin = [-100;-100]; %Min dIp

%% Initialize output variables
P = zeros(n_pwrmod,1);
Q = zeros(n_pwrmod,1);
dP_states = cell(n_pwrmod,1);
dQ_states = cell(n_pwrmod,1);
P_statesIni = cell(n_pwrmod,1);
Q_statesIni = cell(n_pwrmod,1);

%% Define and initialize states at t = 0.
if Flag==0
    for k=1:n_pwrmod
        Q_statesIni{k} = zeros(nOrderQ(k),1);
        P_statesIni{k} = zeros(nOrderP(k),1);
        n = find(bus(:,1)==pwrmod_con(k,1));
        P_statesIni{k}(2) = bus(n,2);
    end
    pwrmod_data = zeros(length(Time),2); %Store Pref in pwrmod_data
    clear k n

%% Calculate P and Q
elseif Flag==1
    Q = zeros(n_pwrmod,1);
    for k=1:n_pwrmod
        if abs(P_states{k}(2))>1e-8
            P(k) = P_states{k}(1)/P_states{k}(2);
        else
            P(k) = P_states{k}(1)/0.01;
        end
    end

%% Calculate derivatives
elseif Flag==2
    for k=1:n_pwrmod
        %No dynamics for Q, initialize dP
        dQ_states{k} = zeros(nOrderQ(k),1);
        dP_states{k} = zeros(nOrderP(k),1);
        
        %Set Pref and store result in pwrmod_data
        if k==1
            Pref = 0;
            if Time(kSim)>=1 && Time(kSim)<1.5; 
                Pref = Pref-0.0001; %Pulse Pref
            end 
        elseif k==2
            Pref = 0;
            if Time(kSim)>=4 && Time(kSim)<4.5; 
                Pref = Pref + 0.0002; %Pulse Pref
            end 
        end
        pwrmod_data(kSim,k) = Pref;
        
        %Calculate dP_states
        if P_states{k}(1)>Ipmax(k) || P_states{k}(1)<Ipmin(k) 
            dP_states{k}(1) = 0; %Limit set using anti-windup
            dP_states{k}(1) = 0; %Limit set using anti-windup
        else
            dP_states{k}(1) = (-P_states{k}(1) + Pref)/Tpord(k);
            if dP_states{k}(1) > dIpmax(k); dP_states{k}(1) = dIpmax(k); end %Rate limit
            if dP_states{k}(1) < dIpmin(k); dP_states{k}(1) = dIpmin(k); end %Rate limit
        end
        n = find(bus(:,1)==pwrmod_con(k,1));
        VT = abs(bus_v(n,kSim)); %Voltage magnitude at the injection bus
        dP_states{k}(2) = (-P_states{k}(2) + VT)/Tv(k);
        if P_states{k}(2)<=0.01; dP_states{k}(2) = 0; end %Limit set using anti-windup
    end
    clear k
    
end


end

