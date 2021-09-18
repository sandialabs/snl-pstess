function     [Y,nSW,nPV,nPQ,SB] = y_sparse(bus,line)
% Syntax:    [Y,nSW,nPV,nPQ,SB] = y_sparse(bus,line) 
%
% Purpose:   build sparse admittance matrix Y from the line data
%
% Input:     bus  - bus data
%            line - line data
%
% Output:    Y    - admittance matrix
%            nSW  - total number of swing buses
%            nPV  - total number generator buses
%            nPQ  - total number of load buses
%            SB   - internal bus numbers of swing bus
%
% See also:  
%
% Calls:    
%
% Called By:   loadflow, form_j, calc

% (c) Copyright 1994-1996 Joe Chow - All Rights Reserved
%
% History (in reverse chronological order)
%
% Version:   2.0
% Author:    Graham Rogers
% Date:      April 1994
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
ibus = (1:nbus)';
bus_int(round(bus(:,1))) = ibus;


% process line data and build admittance matrix Y
  r = line(:,3);
  rx = line(:,4);
  chrg =jay*sparse(diag( 0.5*line(:,5)));
  z = r + jay*rx;     % line impedance
  y = sparse(diag(ones(nline,1)./z));

% determine connection matrices including tap changers and phase shifters
  from_bus = round(line(:,1));
  from_int = bus_int(from_bus);
  to_bus = round(line(:,2));
  to_int = bus_int(to_bus);
  tap_index = find(abs(line(:,6))>0);
  tap=ones(nline,1);
  tap(tap_index)=1. ./line(tap_index,6);
  phase_shift = line(:,7);
  tap = tap.*exp(-jay*phase_shift*pi/180); % Graham's code with - sign
  % tap = tap.*exp(jay*phase_shift*pi/180);   % Joe's original code 

  % sparse matrix formulation
  iline = [1:1:nline]';
  C_from = sparse(from_int,iline,tap,nbus,nline,nline);
  C_to = sparse(to_int,iline,ones(nline,1),nbus,nline,nline);
  C_line = C_from - C_to;

% form Y matrix from primative line ys and connection matrices
  Y=C_from*chrg*C_from' + C_to*chrg*C_to' ;
  Y = Y + C_line*y*C_line';
  Gb = bus(:,8);     % bus conductance
  Bb = bus(:,9);     % bus susceptance

% add diagonal shunt admittances
  Y = Y + sparse(ibus,ibus,Gb+jay*Bb,nbus,nbus);

if nargout > 1
  % count buses of different types
  nSW = 0;
  nPV = 0;
  nPQ = 0;
  bus_type=round(bus(:,10));
  load_index=find(bus_type==3);
  gen_index=find(bus_type==2);
  SB=find(bus_type==1);
  nSW=length(SB);
  nPV=length(gen_index);
  nPQ=length(load_index);
end

return

