function trip_indx(bus)
% Syntax: trip_indx(bus)
%
% Purpose: Determines the relationship between tripped buses and nc loads

%-----------------------------------------------------------------------------%

global g;  % declaring struct of global variables

if g.trip.enable
    nbus = size(bus,1);
    n_trip_ncl = 0;

    g.trip.bus_lshed_idx = zeros(nbus,1);   % index relative to load_con
    g.trip.bus_lshed_type = false(nbus,1);  % non-conforming load flag
    g.trip.bus_lshed_map = zeros(nbus,1);   % index relative to trip signal

    for j = 1:nbus
        index = find(bus(j,1) == g.ncl.load_con(:,1));
        if ~isempty(index)
            % can also check any(g.ncl.load_con(index,[2:4]) > 0)
            n_trip_ncl = n_trip_ncl + 1;

            g.trip.bus_lshed_idx(j) = index;
            g.trip.bus_lshed_type(j) = true;
            g.trip.bus_lshed_map(j) = n_trip_ncl;
        end
    end

    % number of trippable non-conforming loads
    g.trip.n_trip_ncl = n_trip_ncl;

end  % function end

% eof
