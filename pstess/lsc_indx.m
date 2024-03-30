function lsc_indx(bus)
% Syntax: lsc_indx(bus)
%
% Purpose: Forms indexes for the lsc and determines indexes for the
%          inverter-based resources to which the lsc are connected

%-----------------------------------------------------------------------------%
% Version history
%
% Version: 1.0 (initial version)
% Author:  Ryan Elliott
% Date:    November 2020
%-----------------------------------------------------------------------------%

global g;  % declaring struct of global variables

g.lsc.n_sensor = 0;
g.lsc.sensor_set = [];

g.lsc.lsc_idx = [];
g.lsc.lsc_ess_idx = [];
g.lsc.lsc_sensor_idx = [];

if ~isempty(g.lsc.lsc_con)
    g.lsc.n_lsc = size(g.lsc.lsc_con,1);
    g.lsc.lsc_idx = (1:g.lsc.n_lsc).';
    g.lsc.lsc_ess_idx = zeros(g.lsc.n_lsc,1);
    g.lsc.lsc_sensor_idx = zeros(g.lsc.n_lsc,1);

    % finding list of buses without generators or lsc
    bus_nomac = bus(1:end-1,1);
    for ii = 1:length(bus_nomac)
        if ismember(bus_nomac(ii),g.mac.mac_con(:,2))
            bus_nomac(ii) = 0;
        end
        %
        if ismember(bus_nomac(ii),g.lsc.lsc_con(:,2))
            bus_nomac(ii) = 0;
        end
    end

    bus_nomac = bus_nomac(bus_nomac ~= 0);

    % the sensor set includes all buses with lsc
    sensor_set = [g.bus.bus_int(bus_nomac(1:3:end));
                  g.bus.bus_int(g.lsc.lsc_con(:,2))];

    % g.lsc.sensor_set -- index vector for sensor buses
    % g.lsc.n_sensor -- number of sensors
    g.lsc.sensor_set = sort(sensor_set);
    g.lsc.n_sensor = length(sensor_set);

    % building index vectors corresponding to ess_con and sensor_set
    for jj = 1:g.lsc.n_lsc
        ess_idx = find(g.lsc.lsc_con(jj,2) == g.ess.ess_con(:,2));
        if ~isempty(ess_idx)
            g.lsc.lsc_ess_idx(jj) = ess_idx;        % index relative to ess_con
        else
            estr = '\nlsc_indx: the lsc at bus %0.0f must be connected to ';
            estr = [estr, 'an ess model.'];
            error(sprintf(estr,g.lsc.lsc_con(jj,2)));
        end
        %
        sensor_idx = find(g.bus.bus_int(g.lsc.lsc_con(jj,2)) == g.lsc.sensor_set);
        if ~isempty(sensor_idx)
            g.lsc.lsc_sensor_idx(jj) = sensor_idx;  % index relative to sensor_set
        else
            estr = '\nlsc_indx: the lsc at bus %0.0f must also have an ';
            estr = [estr, 'angle sensor.'];
            error(sprintf(estr,g.lsc.lsc_con(jj,2)));
        end
    end
end

end  % function end

% eof
