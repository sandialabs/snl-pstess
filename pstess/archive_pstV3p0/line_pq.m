function [S1,S2] = line_pq(V1,V2,R,X,B,tap,phi)
% Syntax:   [S1,S2] = line_pq(V1,V2,R,X,B,tap,phi) 
%
% Purpose:  Compute line flows. Inputs can be vectors and matrices.
%	
% Input:    V1        - from bus complex voltage matrix
%           V2        - to bus complex voltage matrix
%           R         - line resistance vector
%           X         - line reactance vector
%           B         - line charging vector
%           tap       - tap ratio vector
%           phi       - phase shifter angle vector in degrees
% Output:   S1        - complex power injection matrix at from bus
%           S2        - complex power injection matrix at to bus
%
% See also:  
%
% Algorithm: Assumes that V1 and V2 are matrices of bus voltages
%            in the form v(:,1:k) where each column is the voltage at
%            a time step j. V1 and V2 must have the same size.
% 	     The tap is at the from bus and represents the step down
%            ratio i.e. V1' = V1/t*exp(jphi*pi/180); 
%            i1' = i1*t*exp(jphi*pi/180)
% Application: To calculate the complex power flow from transient simulation 
%              records
%	       Set V1 = bus_v(bus_int(line(:,1)),:) the from bus voltages
%              Set V2 = bus_v(bus_int(line(:,2)),:) the to   bus voltages
%              Set R = line(;,3); X = line(:,4); B = line(;,5)
%              Set tap = line(:,6); phi = line(:,7)
%              The flow on any line may then be plotted using
%              plot(t,real(S1(line_nu,:)) for example    
%
% (c) Copyright 1991-2 Joe H. Chow - All Rights Reserved
%
% History (in reverse chronological order)
% Version:   2.0
% Author:    Graham Rogers
% date:      October 1996
% Version:   1.0
% Authors:   Joe H. Chow
% Date:      March 1992
%
% ***********************************************************
jay = sqrt(-1);
[nline,dummy] = size(V1);
for i = 1:nline
  if tap(i) == 0
    tap(i) = 1;
  end
end
tps = tap.*exp(jay*phi*pi/180);
tpsi = diag(ones(nline,1)./tps);
tps = diag(tps);
z = R + jay*X;
y = diag(ones(nline,1)./z);
chg = diag(jay*B/2);
cur1 = tps*(y*(tpsi*V1-V2) + chg*V1);
cur2 = y*(V2 - tpsi*V1) + chg*V2;
S1 = V1.*conj(cur1);
S2 = V2.*conj(cur2);
