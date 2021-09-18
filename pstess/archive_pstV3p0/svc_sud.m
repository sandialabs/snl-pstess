function [y,xn,dx] = svc_sud(i,k,flag,s,d_sig,ymax,ymin,x)
% Syntax: [y,xn,dx] = svc_sud(i,k,flag,s,d_sig,ymax,ymin,x)
% 8:13 am 11/12/98
%
% Purpose: svc user defined damping control           
% Input: i - svc number
%        k - integer time
%        flag - 0 - initialization
%               1 - network interface computation
%               2 - system dynamics computation
%        s - control state space object
%        d_sig - damping control input
%        x - state
%
% Output: 
%        y -  controller output
%        xn - state with limits applied
%        dx - rate of change of state
%
% Files:
%
% See Also: svc
% Algorithm:
%
% Calls:
%
% Called By:

% (c) Copyright Graham Rogers 1998-2008 - All Rights Reserved

% History (in reverse chronological order)
% Author:   Graham Rogers
% Date:     November 1998

y=0;

if flag == 0; % initialization
   if i ~= 0  % scalar computation
      % check user defined control
      sv = get(s);
      if sv.NumInputs~=1;error('ud svc stabilizer must be single input');end
      if sv.NumOutputs~=1;error('ud svc stabilizer must be single output');end 
      [y,xn]=eval(s,0,d_sig);
      if y>ymax;warning('y outside max limit initialy');y=ymax;end
      if y<ymin;warning('y outside minimum limit initially');y=ymin;end   
      dx=zeros(size(xn));         
   else
      error('no vector computation in user defined svc stabilizer')
   end
end

if flag == 1 % network interface computation
   %no interface computation in user defined svc stabilizer
end

if flag == 2 % svc damping control dynamics calculation
   if i ~= 0 % scalar computation
      xmax = 1e5*ones(s.NumStates,1);xmin = -xmax;
      [y,xn,dx]=dstate(s,x,d_sig,xmax,xmin,ymax,ymin);
   else
      error('no vector computation for user defined svc stabilizer')
   end
end

