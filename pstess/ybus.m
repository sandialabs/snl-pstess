function [Y,nSW,nPV,nPQ,SB] = ybus(bus,line)
% Syntax: [Y,nSW,nPV,nPQ,SB] = ybus(bus,line)
%
% Purpose: build admittance matrix Y from the line data
%
% Input:   bus  - bus data
%          line - line data
%
% Output:  Y    - admittance matrix
%          nSW  - total number of swing buses
%          nPV  - total number generator buses
%          nPQ  - total number of load buses
%          SB   - bus number of swing bus
%
% Called by: loadflow

%-----------------------------------------------------------------------------%
% Version history
%
% Version: 1.0
% Author:  Kwok W. Cheung, Joe H. Chow
% Date:    March 1991
%-----------------------------------------------------------------------------%

global g;  % declaring struct of global variables

% note: bus_type is a file local variable that is not the same as g.lfac.bus_type

swing_bus = 1;
gen_bus = 2;
load_bus = 3;

n_line = length(line(:,1));  % number of lines
n_bus = length(bus(:,1));    % number of buses

r = zeros(n_line,1);
rx = zeros(n_line,1);
chrg = zeros(n_line,1);
z = zeros(n_line,1);
y = zeros(n_line,1);

Y = sparse(1,1,0,n_bus,n_bus);

% set up internal bus numbers for second indexing of buses
busmax = max(bus(:,1));
g.bus.bus_int = zeros(busmax,1);
ibus = [1:1:n_bus]';
for ii = 1:n_bus
    g.bus.bus_int(round(bus(ii,1))) = ii;
end

% process line data and build admittance matrix Y
r = line(:,3);
rx = line(:,4);
chrg = line(:,5);
z = r + 1j*rx;               % line impedance
y = ones(n_line,1)./z;

for ii = 1:n_line
    from_bus = line(ii,1);
    from_int = g.bus.bus_int(round(from_bus));

    to_bus = line(ii,2);
    to_int = g.bus.bus_int(round(to_bus));

    tap_ratio = line(ii,6);
    if (tap_ratio == 0)      % this line has no transformer
        tap_ratio = 1;
    end

    phase_shift = line(ii,7);
    tps = tap_ratio*exp(1j*phase_shift*pi/180);
    tps2 = tap_ratio*tap_ratio;

    % sparse matrix formulation
    jy1(1,1) = from_int;
    jy2(1,1) = to_int;
    w(1,1) = -y(ii)/conj(tps);

    jy1(2,1) = to_int;
    jy2(2,1) = from_int;
    w(2,1) = -y(ii)/tps;

    jy1(3,1) = from_int;
    jy2(3,1) = from_int;
    w(3,1) = (y(ii) + 1j*chrg(ii)/2)/tps2;

    jy1(4,1) = to_int;
    jy2(4,1) = to_int;
    w(4,1) = y(ii) + 1j*chrg(ii)/2;

    Y = Y + sparse(jy1,jy2,w,n_bus,n_bus);
end

Gb = bus(:,8);  % bus conductance
Bb = bus(:,9);  % bus susceptance

% add diagonal shunt admittances
Y = Y + sparse(ibus,ibus,Gb+1j*Bb,n_bus,n_bus);

if (nargout > 1)
    % count buses of different types
    nSW = 0;
    nPV = 0;
    nPQ = 0;
    for ii = 1:n_bus,
        bus_type(ii) = bus(ii,10);
        if (bus_type(ii) == swing_bus)
            SB = g.bus.bus_int(round(bus(ii,1)));  % swing bus number
            nSW = nSW + 1;                         % increment swing bus counter
        elseif (bus_type(ii) == gen_bus)
            nPV = nPV + 1;                         % increment generator bus counter
        else
            nPQ = nPQ + 1;                         % increment load bus counter
        end
    end
end

end  % function end

% eof
