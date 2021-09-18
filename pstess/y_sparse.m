function [Y,nSW,nPV,nPQ,SB] = y_sparse(bus,line)
% Syntax: [Y,nSW,nPV,nPQ,SB] = y_sparse(bus,line)
%
% Purpose: build sparse admittance matrix Y from the line data
%
% Input:   bus  - bus data
%          line - line data
%
% Output:  Y    - admittance matrix
%          nSW  - total number of swing buses
%          nPV  - total number generator buses
%          nPQ  - total number of load buses
%          SB   - internal bus numbers of swing bus
%
% Called By: loadflow, form_j, calc

%-----------------------------------------------------------------------------%
% Version history
%
% Version:   2.0
% Author:    Graham Rogers
% Date:      April 1994
% Version:   1.0
% Author:    Kwok W. Cheung, Joe H. Chow
% Date:      March 1991
%-----------------------------------------------------------------------------%

global g;  % declaring struct of global variables

% note: bus_type is a file local variable that is not the same as g.lfac.bus_type

swing_bus = 1;
gen_bus = 2;
load_bus = 3;

n_line = length(line(:,1));             % number of lines
n_bus = length(bus(:,1));               % number of buses
r = zeros(n_line,1);
rx = zeros(n_line,1);
chrg = zeros(n_line,1);
z = zeros(n_line,1);
y = zeros(n_line,1);

Y = sparse(1,1,0,n_bus,n_bus);

% set up internal bus numbers for second indexing of buses
busmax = max(bus(:,1));
g.bus.bus_int = zeros(busmax,1);
ibus = (1:n_bus)';
g.bus.bus_int(round(bus(:,1))) = ibus;

% process line data and build admittance matrix Y
r = line(:,3);
rx = line(:,4);
chrg = 1j*sparse(diag(0.5*line(:,5)));
z = r + 1j*rx;                          % line impedance
y = sparse(diag(ones(n_line,1)./z));

% determine connection matrices including tap changers and phase shifters
from_bus = round(line(:,1));
from_int = g.bus.bus_int(from_bus);
to_bus = round(line(:,2));
to_int = g.bus.bus_int(to_bus);
tap_index = find(abs(line(:,6))>0);
tap = ones(n_line,1);
tap(tap_index) = 1./line(tap_index,6);  % syntax correction
phase_shift = line(:,7);
tap = tap.*exp(-1j*phase_shift*pi/180);

% sparse matrix formulation
iline = [1:1:n_line]';
C_from = sparse(from_int,iline,tap,n_bus,n_line,n_line);
C_to = sparse(to_int,iline,ones(n_line,1),n_bus,n_line,n_line);
C_line = C_from - C_to;

% form Y matrix from primative line ys and connection matrices
Y = C_from*chrg*C_from' + C_to*chrg*C_to';
Y = Y + C_line*y*C_line';
Gb = bus(:,8);                          % bus conductance
Bb = bus(:,9);                          % bus susceptance

% add diagonal shunt admittances
Y = Y + sparse(ibus,ibus,Gb+1j*Bb,n_bus,n_bus);

if (nargout > 1)
    % count buses of different types
    nSW = 0;
    nPV = 0;
    nPQ = 0;
    bus_type = round(bus(:,10));
    load_index = find(bus_type == 3);
    gen_index = find(bus_type == 2);
    SB = find(bus_type == 1);
    nSW = length(SB);
    nPV = length(gen_index);
    nPQ = length(load_index);
end

end  % function end

% eof
