function [delP,delQ,P,Q,conv_flag] = ...
                 calc(V,ang,Y,Pg,Qg,Pl,Ql,sw_bno,g_bno,tol)
% Syntax:  [delP,delQ,P,Q,conv_flag] = 
%                calc(V,ang,Y,Pg,Qg,Pl,Ql,sw_bno,g_bno,tol)
%
% Purpose: calculates power mismatch and checks convergence
%          also determines the values of P and Q based on the 
%          supplied values of voltage magnitude and angle
% Version: 2.0 eliminates do loop
% Input:   nbus      - total number of buses
%          bus_type  - load_bus(3), gen_bus(2), swing_bus(1)
%          V         - magnitude of bus voltage
%          ang       - angle(rad) of bus voltage
%          Y         - admittance matrix
%          Pg        - real power of generation
%          Qg        - reactive power  of generation
%          Pl        - real power of load
%          Ql        - reactive power of load
%	 sw_bno - a vector having zeros at all  swing_bus locations ones otherwise
%	 g_bno  - a vector having zeros at all  generator bus locations ones otherwise
%          tol       - a tolerance of computational error
%
% Output:  delP      - real power mismatch
%          delQ      - reactive power mismatch
%          P         - calculated real power
%          Q         - calculated reactive power
%          conv_flag - 0, converged
%                      1, not yet converged
%
% See also:  
%
% Calls:
%
% Called By:   loadflow

% (c) Copyright 1991 Joe H. Chow - All Rights Reserved
%
% History (in reverse chronological order)
% Version:   2.0
% Author:    Graham Rogers
% Date:      July 1994
%
% Version:   1.0
% Author:    Kwok W. Cheung, Joe H. Chow
% Date:      March 1991
%
% ************************************************************
jay = sqrt(-1);
swing_bus = 1;
gen_bus = 2;
load_bus = 3;
% voltage in rectangular coordinate
V_rect = V.*exp(jay*ang);  
% bus current injection
cur_inj = Y*V_rect;
% power output based on voltages 
S = V_rect.*conj(cur_inj);
P = real(S); Q = imag(S);
delP = Pg - Pl - P;
delQ = Qg - Ql - Q;
% zero out mismatches on swing bus and generation bus
delP=delP.*sw_bno;
delQ=delQ.*sw_bno;
delQ=delQ.*g_bno;
%  total mismatch
[pmis,ip]=max(abs(delP));
[qmis,iq]=max(abs(delQ));
mism = pmis+qmis;
if mism > tol,
    conv_flag = 1;
  else
    conv_flag = 0;
end
return

