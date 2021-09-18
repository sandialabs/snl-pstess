function     [Y,nSW,nPV,nPQ,SB] = ybus(bus,line)
% Syntax:    [Y,nSW,nPV,nPQ,SB] = ybus(bus,line) 
%
% Purpose:   build admittance matrix Y from the line data
%
% Input:     bus  - bus data
%            line - line data
%
% Output:    Y    - admittance matrix
%            nSW  - total number of swing buses
%            nPV  - total number generator buses
%            nPQ  - total number of load buses
%            SB   - bus number of swing bus
%
% See also:  
%
% Calls:     calc, form_jac
%
% Call By:   loadflow

% (c) Copyright 1991 Joe H. Chow - All Rights Reserved
%
% History (in reverse chronological order)
%
% Version:   1.0
% Author:    Kwok W. Cheung, Joe H. Chow
% Date:      March 1991
%
% ************************************************************
global bus_int

jay = sqrt(-1);
swing_bus = 1;
gen_bus = 2;
load_bus = 3;

nline = length(line(:,1));     % number of lines
nbus = length(bus(:,1));     % number of buses
r=zeros(nline,1);
rx=zeros(nline,1);
chrg=zeros(nline,1);
z=zeros(nline,1);
y=zeros(nline,1);

Y = sparse(1,1,0,nbus,nbus);

% set up internal bus numbers for second indexing of buses
busmax = max(bus(:,1));
bus_int = zeros(busmax,1);
ibus = [1:1:nbus]';
for i = 1:nbus
  bus_int(round(bus(i,1))) = i;
end

% process line data and build admittance matrix Y
  r = line(:,3);
  rx = line(:,4);
  chrg = line(:,5);
  z = r + jay*rx;     % line impedance
  y = ones(nline,1)./z;

for i = 1:nline
  from_bus = line(i,1);
  from_int = bus_int(round(from_bus));
  to_bus = line(i,2);
  to_int = bus_int(round(to_bus));
  tap_ratio = line(i,6);
  if tap_ratio == 0,     % this line has no transformer
    tap_ratio = 1;
  end
  phase_shift = line(i,7);
  tps = tap_ratio*exp(jay*phase_shift*pi/180);
  tps2=tap_ratio*tap_ratio;
  % sparse matrix formulation
  j1(1,1) = from_int; j2(1,1) = to_int;
  w(1,1) = - y(i)/conj(tps);
  j1(2,1) = to_int; j2(2,1) = from_int;
  w(2,1) = - y(i)/tps;
  j1(3,1) = from_int; j2(3,1) = from_int;
  w(3,1) = (y(i) + jay*chrg(i)/2)/tps2;
  j1(4,1) = to_int; j2(4,1) = to_int;
  w(4,1) = y(i) + jay*chrg(i)/2;
  Y = Y + sparse(j1,j2,w,nbus,nbus);
end; %
Gb = bus(:,8);     % bus conductance
Bb = bus(:,9);     % bus susceptance
% add diagonal shunt admittances
Y = Y + sparse(ibus,ibus,Gb+jay*Bb,nbus,nbus);


if nargout > 1
  % count buses of different types
  nSW = 0;
  nPV = 0;
  nPQ = 0;
  for i = 1:nbus,
    bus_type(i) = bus(i,10);
    if bus_type(i) == swing_bus,
	SB = bus_int(round(bus(i,1))); % swing bus number 
	nSW = nSW + 1;     % increment swing bus counter
      elseif bus_type(i) == gen_bus,
	nPV = nPV +1;      % increment generator bus counter
      else
	nPQ = nPQ + 1;     % increment load bus counter
    end
  end; %
end

return

