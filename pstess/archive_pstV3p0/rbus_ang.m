function [caf,cat] = rbus_ang(V1,V2,R,X,B,tap,phi,c_vf,c_vt,c_if,c_it)
% Syntax:   [caf,cat] = rbus_ang(V1,V2,R,X,B,tap,phi,c_vf,c_vt,c_if,c_it) 
%
% Purpose:  Compute bus angle changes. Inputs can be vectors and matrices.
%	
% Input:    V1        - from bus complex voltage matrix
%           V2        - to bus complex voltage matrix
%           R         - line resistance vector
%           X         - line reactance vector
%           B         - line charging vector
%           tap       - tap ratio vector
%           phi       - phase shifter angle vector in degrees
%           c_vf       - complex vector for small changes in voltage 
%                       at V1
%           c_vt       - complex vector for small changes in voltage 
%                       at V2
%           c_if       - complex vector for small changes in current 
%                       at from bus
%           c_it       - complex vector for small changes in current 
%                       at to bus
% Output:   caf        - c matrix for change in angle at from bus due to
%                        a change in V2
%           cat        - c matix for change in angle at to bus due to
%                        a change in V1
%               
%
% See also:  
%
% Algorithm: Assumes that V1 and V2 are matrices of bus voltages
%            in the form v(:,1:k) where each column is the voltage at
%            a time step j. V1 and V2 must have the same size.
% 	     The tap is at the from bus and represents the step down
%            ratio i.e. V1' = V1/t*exp(jphi*pi/180); 
%            i1' = i1*t*exp(jphi*pi/180)
% Application: To calculate the line current from transient simulation 
%              records
%	       Set V1 = bus_v(bus_int(line(:,1)),:) the from bus voltages
%              Set V2 = bus_v(bus_int(line(:,2)),:) the to   bus voltages
%              Set R = line(;,3); X = line(:,4); B = line(;,5)
%              Set tap = line(:,6); phi = line(:,7)
%              The flow on any line may then be plotted using
%              plot(t,abs(l_i(line_nu,:)) for example    
%
% (c) Copyright 1991-8 Cherry Tree Scientific Software/Joe H. Chow - All Rights Reserved
%
% History (in reverse chronological order)
% Version:   1.0
% Author:    Graham Rogers
% date:      March 1998
%
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
z = diag(z);
chg = diag(jay*B/2);
il = y*(tpsi*V1-V2);
cf = tps*(c_vt-z*chg*c_vt-z*c_it);
caf = imag(diag(1./V1)*cf);
ct = tpsi*(c_vf + z*chg*c_vf) - tps*z*c_if;
cat = imag(diag(1./V2)*ct);

