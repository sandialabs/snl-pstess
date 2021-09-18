function f = tg(i,k,bus,flag)
% Syntax: f = tg(i,k,bus,flag)  
% 1:19 PM 15/08/97
% Purpose: simple turbine governor model
%           
% Input: i - generator number
%        k - integer time
%        bus - solved loadflow bus data
%        flag - 0 - initialization
%               1 - network interface computation
%               2 - system dynamics computation
%
% Output: f - dummy variable 
%
% tg_con matrix format
%column	       data			unit
%  1	turbine model number (=1)	
%  2	machine number	
%  3	speed set point   wf		pu
%  4	steady state gain 1/R		pu
%  5	maximum power order  Tmax	pu on generator base
%  6	servo time constant   Ts	sec
%  7	governor time constant  Tc	sec
%  8	transient gain time constant T3	sec
%  9	HP section time constant   T4	sec
% 10	reheater time constant    T5	sec

%
% Files:
%
% See Also: pst_var

% Algorithm:
%
% Calls:
%
% Call By:

% (c) Copyright 1991-3 Joe H. Chow - All Rights Reserved

% History (in reverse chronological order)
%
% Version:
% Date:
% Author:
% Purpose:
% Modification:

% Version:  1.0
% Author:   Joe H. Chow
% Date:     August 1993

global  basmva mac_int mac_con pelect pmech mac_spd 
global  tg_con tg_pot 
global  tg1 tg2 tg3 dtg1 dtg2 dtg3 
global  tg_idx n_tg tg_sig


jay = sqrt(-1);
f = 0;
if flag == 0; % initialization
   if i ~= 0
      if tg_con(i,1) ~= 1
         error('TG: requires tg_con(i,1) = 1')
      end
   end
   if i ~= 0  % scalar computation
      n = mac_int(tg_con(i,2)); % machine number
      
      if pmech(n,k) > tg_con(i,5) 
         error('TG init: pmech > upper limit, check machine base')
      end
      if pmech(n,k) < 0
         error('TG init: pmech < 0, check data')
      end
      tg1(i,1) = pmech(n,k);
      %
      tg_pot(i,1) = tg_con(i,8)/tg_con(i,7);
      a1 = 1 - tg_pot(i,1);
      tg_pot(i,2) = a1;
      tg2(i,1) = a1*pmech(n,k);
      %    
      tg_pot(i,3) = tg_con(i,9)/tg_con(i,10);
      a2 = 1 - tg_pot(i,3);
      tg_pot(i,4) = a2;
      tg3(i,1) = a2*pmech(n,k);
      %
      tg_pot(i,5) = pmech(n,k);
      %
      tg_sig(i,1)=0;
   else
      %  vectorized computation
      if n_tg~=0
         n = mac_int(tg_con(tg_idx,2)); % machine number
         maxlmt = find(pmech(n,1)>tg_con(tg_idx,5));
         if ~isempty(maxlmt)
            n(maxlmt)
            error(' pmech excedes maximum limit')
         end
         minlmt = find(pmech(n,1)<zeros(n_tg,1));
         if ~isempty(minlmt)
            n(minlmt)
            error('pmech less than zero')
         end
         tg1(tg_idx,1) = pmech(n,1);
         %
         tg_pot(tg_idx,1) = tg_con(tg_idx,8)./tg_con(tg_idx,7);
         a1 = ones(n_tg,1) - tg_pot(tg_idx,1);
         tg_pot(tg_idx,2) = a1;
         tg2(tg_idx,1) = a1.*pmech(n,k);
         %    
         tg_pot(tg_idx,3) = tg_con(tg_idx,9)./tg_con(tg_idx,10);
         a2 = ones(n_tg,1) - tg_pot(tg_idx,3);
         tg_pot(tg_idx,4) = a2;
         tg3(tg_idx,1) = a2.*pmech(n,k);
         %
         tg_pot(tg_idx,5) = pmech(n,k);% set reference value
         tg_sig(tg_idx,1) = zeros(n_tg,1);
      end
   end
end

if flag == 1 % network interface computation
   if i ~= 0 % scalar computation
      n = mac_int(tg_con(i,2)); % machine number
      % the following update is needed because pmech depends on 
      %   the output of the states tg1, tg2 and tg3
      pmech(n,k) = tg3(i,k) + ...
         tg_pot(i,3)*( tg2(i,k) + tg_pot(i,1)*tg1(i,k) );
   else
      if n_tg~=0
         n = mac_int(tg_con(tg_idx,2)); % machine number
         pmech(n,k) = tg3(tg_idx,k) + ...
            tg_pot(tg_idx,3).*( tg2(tg_idx,k) + tg_pot(tg_idx,1).*tg1(tg_idx,k) );
      end
   end
end

if flag == 2 % turbine governor dynamics calculation
   if i ~= 0 % scalar computation
      n = mac_int(tg_con(i,2)); % machine number
      spd_err = tg_con(i,3) - mac_spd(n,k);
      demand = tg_pot(i,5) + spd_err*tg_con(i,4) + tg_sig(i,k);
      demand = min( max(demand,0),tg_con(i,5) );
      dtg1(i,k) = (demand - tg1(i,k))/tg_con(i,6);
      %
      dtg2(i,k) = (tg_pot(i,2)* tg1(i,k)-tg2(i,k))/tg_con(i,7);
      % 
      dtg3(i,k) = ( (tg2(i,k)+tg_pot(i,1)*tg1(i,k))*tg_pot(i,4) ...
         -tg3(i,k) )/tg_con(i,10);
      pmech(n,k) = tg3(i,k) + ...
         tg_pot(i,3)*(tg2(i,k) + tg_pot(:,1)*tg1(i,k));
   else
      % vectorized computation
      if n_tg ~=0
         n = mac_int(tg_con(tg_idx,2)); % machine number
         spd_err = tg_con(tg_idx,3) - mac_spd(n,k);
         demand = tg_pot(tg_idx,5) + spd_err.*tg_con(tg_idx,4) + tg_sig(tg_idx,k);
         demand = min( max(demand,zeros(n_tg,1)),tg_con(tg_idx,5) );
         dtg1(tg_idx,k) = (demand - tg1(tg_idx,k))./tg_con(tg_idx,6);
         %
         dtg2(tg_idx,k) = ( tg1(tg_idx,k).*tg_pot(tg_idx,2)-tg2(tg_idx,k))...
            ./tg_con(tg_idx,7);
         % 
         dtg3(tg_idx,k) = ((tg2(tg_idx,k)+tg_pot(tg_idx,1).*tg1(tg_idx,k))...
            .*tg_pot(tg_idx,4)...
            -tg3(tg_idx,k))./tg_con(tg_idx,10);
         pmech(n,k) = tg3(tg_idx,k) + ...
            tg_pot(tg_idx,3).*(tg2(tg_idx,k)...
            + tg_pot(tg_idx,1).*tg1(tg_idx,k));
      end
   end
end

