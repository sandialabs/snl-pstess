function [tripOut,mac_trip_states] = mac_trip_logic(tripStatus,mac_trip_states,t,kT)
% Purpose: trip generators.
%
% Inputs:
%   tripStatus = n_mac x 1 bool vector of current trip status.  If
%       tripStatus(n) is true, then the generator corresponding to the nth
%       row of mac_con is already tripped.  Else, it is false.
%   mac_trip_states = storage matrix defined by user.
%   t = vector of simulation time (sec.).
%   kT = current integer time (sample).  Corresponds to t(kT)
%
% Output:
%   tripOut = n_mac x 1 bool vector of desired trips.  If
%       tripOut(n)==1, then the generator corresponding to the nth
%       row of mac_con is will be tripped.  Note that each element of
%       tripOut must be either 0 or 1.

%-----------------------------------------------------------------------------%
% Version: 1.0
% Author:  Dan Trudnowski
% Date:    Jan 2017
%-----------------------------------------------------------------------------%

global g;  % declaring struct of global variables

% index of global variables frequently used for protection logic
% g.mac.n_mac        % number of synchronous generators in the base case
% g.mac.mac_con      % matrix of generator dynamic model parameters
% g.bus.bus_v        % bus_v(n,kT)   -- bus n pu voltage at time t(kT)
% g.mac.pelect       % pelect(n,kT)  -- gen n pu real power at time t(kT)
% g.mac.qelect       % qelect(n,kT)  -- gen n pu reactive power at time t(kT)
% g.mac.mac_spd      % mac_spd(n,kT) -- gen n pu speed at time t(kT)
% g.mac.cur_re       % cur_re        -- generator real current (pu)
% g.mac.cur_im       % cur_im        -- generator imag current (pu)

% Initialize
tripOut = false(g.mac.n_mac,1);
mac_trip_states = 0;

end  % function end

% eof
