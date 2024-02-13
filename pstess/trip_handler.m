function [bus,line] = trip_handler(t,k,flag,bus,line)
% Syntax: [bus,line] = trip_handler(t,k,flag,bus,line)
%
% Purpose: Supports protection and RAS tripping of generators, loads, and lines
%
% Input:   t - simulation time at index k
%          k - integer time
%          flag - true - predictor step (time k within s_simu)
%                 false - corrector step (time j within s_simu)
%
% See Also: trip_logic, trip_load
%-----------------------------------------------------------------------------%

global g;  % declaring struct of global variables

if flag    % flag is true in the predictor step, false otherwise
    % trip/re-close gens, loads and lines
    [f1,f2,f3] = trip_logic(t,k,bus,line);

    % there is an extra bus in pf1 and pf2 inserted at the fault location
    g.trip.mac_trip_flags(:,k) = (g.trip.mac_trip_flags(:,max(k-1,1)) | logical(f1));
    g.trip.load_trip_flags(:,k) = logical(f2(1:size(g.trip.load_trip_flags,1)));
    g.trip.line_trip_flags(:,k) = logical(f3(1:size(g.trip.line_trip_flags,1)));
else
    % apply gen, load, and line trips (re-closure is automatic in bus and line)
    kk = max(k-1,1);
    g.trip.mac_trip_flags(:,k) = g.trip.mac_trip_flags(:,kk);
    g.trip.load_trip_flags(:,k) = g.trip.load_trip_flags(:,kk);
    g.trip.line_trip_flags(:,k) = g.trip.line_trip_flags(:,kk);
end

trip_flag = any(g.trip.mac_trip_flags(:,k)) || any(g.trip.load_trip_flags(:,k));
trip_flag = trip_flag || any(g.trip.line_trip_flags(:,k));

if trip_flag
    if any(g.trip.mac_trip_flags(:,k))   % gen trips
        nB = g.mac.mac_con(g.trip.mac_trip_flags(:,k),2);  % tripped gen buses
        for kB = 1:length(nB)
            % make reactance infinite for all lines connected to tripped gen
            nL = find(nB(kB) == line(:,1) | nB(kB) == line(:,2));
            line(nL,4) = 1e7;
        end
    end

    if any(g.trip.load_trip_flags(:,k))  % load trips
        n_ncl = find(g.trip.load_trip_flags(:,k) & g.trip.bus_lshed_type);
        n_zc = find(g.trip.load_trip_flags(:,k) & ~g.trip.bus_lshed_type);

        % updating non-conforming loads (passed to nc_load)
        if ~isempty(n_ncl)
            jj = g.trip.bus_lshed_map(n_ncl);
            g.trip.load_trip_ncl(jj,k) = (bus(n_ncl,6) + 1j*bus(n_ncl,7)) ...
                                         .*(-1*g.trip.load_trip_frac(n_ncl,k));
        end

        % updating constant-impedance loads
        if ~isempty(n_zc)
            bus(n_zc,[6,7]) = bus(n_zc,[6,7]).*(1 - g.trip.load_trip_frac(n_zc,k));
        end
    end

    if any(g.trip.line_trip_flags(:,k))  % line trips
        n = find(g.trip.line_trip_flags(:,k));
        line(n,4) = 1e7*ones(length(n),1);
    end
end

% eof
