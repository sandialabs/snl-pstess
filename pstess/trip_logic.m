function [gen_trip_out,load_trip_out,line_trip_out] = trip_logic(t,k,bus,line)
% Purpose: trip generators.
%
% Inputs:
%   bus = the bus matrix at time t(k)
%   line = the line matrix at time t(k)
%   t = simulation time (sec.)
%   k = current integer time (sample)
%
% Output:
%   gen_trip_out = n_mac x 1 bool vector of desired trips at time t(k). If
%       gen_trip_out(n) == true, then the generator corresponding to the nth
%       row of mac_con will be tripped. Note that each element of
%       gen_trip_out must be either 0 or 1. Once a generator is tripped,
%       it remains tripped for the rest of the simulation.
%   load_trip_out = N x 1 bool vector of desired trips at time kT where N
%       is the number of buses. If load_trip_out(n) == true, then the load
%       corresponding to the nth row of bus_con  will be tripped.
%       Note that each element of load_trip_out must be either false or true.
%       If load_trip_out(n) == false, the load is re-connected.
%   line_trip_out = N x 1 bool vector of desired trips at time kT where N
%       is the number of rows of the line matrix. If load_trip_out(n)==true,
%       then line(n,:) will be tripped. Each element of line_trip_out must
%       be either false or true. If line_trip_out(n)==false, the load is re-connected.

%-----------------------------------------------------------------------------%
% Version: 1.0
% Author:  Dan Trudnowski
% Date:    Jan 2017
%
% Version: 2.0 - updated to use load_trip_out too
% Author:  D. Trudnowski
% Date:    Feb. 2023
%-----------------------------------------------------------------------------%

global g;  % declaring struct of global variables

% index of global variables frequently used for protection logic
% g.mac.n_mac             number of synchronous generators in the base case
% g.mac.mac_con           matrix of generator dynamic model parameters
% g.bus.bus_v             bus_v(n,k)   -- bus n pu voltage at time t(k)
% g.mac.pelect            pelect(n,k)  -- gen n pu real power at time t(k)
% g.mac.qelect            qelect(n,k)  -- gen n pu reactive power at time t(k)
% g.mac.mac_spd           mac_spd(n,k) -- gen n pu speed at time t(k)
% g.mac.cur_re            cur_re       -- generator real current (pu)
% g.mac.cur_im            cur_im       -- generator imag current (pu)
% g.trip.mac_trip_flags   g.trip.mac_trip_flags(n,k) is the gen_trip_out(n) at
%                         time t(k). Note that it is latched.
%                         DO NOT CHANGE VALUES. THIS IS DONE IN s_simu.
% g.trip.load_trip_flags  g.trip.load_trip_flags(n,k) is the load_trip_out(n) at
%                         time t(k). DO NOT CHANGE VALUES. THIS IS DONE IN
%                         s_simu.
% g.trip.line_trip_flags  g.trip.line_trip_flags(n,k) is the line_trip_out(n) at
%                         time t(k). DO NOT CHANGE VALUES. THIS IS DONE IN
%                         s_simu.
% g.trip.user_variables   This is a general storage area used by the user
%                         of this function.

% Initialize key variables
nB = length(find(g.bus.bus_int>0));   % number of buses
nL = size(line,1);                    % number of lines
gen_trip_out = false(g.mac.n_mac,1);
load_trip_out = false(nB,1);
line_trip_out = false(nL,1);

% if (k == 1)
%     g.trip.user_variables = 0;
% end

end  % function end

% eof
