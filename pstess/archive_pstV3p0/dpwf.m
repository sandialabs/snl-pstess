function [f] = dpwf(i,k,bus,flag)
% Syntax: [f] = dpwf(i,k,bus,flag)
% 10:51 am July 8, 1998
%
% Purpose: filter model for deltaP/w stabilizer
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
% dpw_con - input data format
% col1: generator number
% col2: number of lead elements (maximum 4)
% col3: lead time constant
% col4: lag time constant
% 5 stage filter assumed
% History (in reverse chronological order)
% Version: 1.0
% Date: July 1998
% Author:   Graham Rogers
% (c) copyright Cherry Tree Scientific Software 1998 - All rights reserved

% deltaP/w variables 
global  dpw_con dpw_pot dpw_out dpw_pss_idx dpw_mb_idx dpw_idx n_dpw dpw_Td_idx
global  sdpw1 sdpw2 sdpw3 sdpw4 sdpw5 sdpw6
global  dsdpw1 dsdpw2 dsdpw3 dsdpw4 dsdpw5 dsdpw6 
global  mac_con mac_int mac_spd pelect basmva 
f=0;
if n_dpw~=0
   jay = sqrt(-1);   
   if flag == 0; % initialization
      if i ~= 0  % scalar computation
         error('vector computation only for dpwf')
      else
         %vector computation
         ezeros = find(dpw_con(:,2)>4);
         if ~isempty(ezeros)
            error('number of dpw filter zeros must be 4 or less')
         end   
         dpw_pot=ones(n_dpw,2);
         n=dpw_mb_idx(dpw_pss_idx);%generator numbers
         if isempty(n)
            error('there are no generators associated with the power system stabilizers')
         end  
                  
         dpw_out(dpw_pss_idx,1)=zeros(n_dpw,1);
         dpw_pot(:,1)=dpw_con(:,3)./dpw_con(:,4);
         dpw_pot(:,2) = 1/2./mac_con(n,16);% 1/2/H
         dpw_pot(:,3) = basmva./mac_con(n,3); % base power ratio
         dpw_pot(:,4) = dpw_pot(:,1).*dpw_Td_idx(:,1);
         dpw_pot(:,5) = dpw_pot(:,1).*dpw_Td_idx(:,2);
         dpw_pot(:,6) = dpw_pot(:,1).*dpw_Td_idx(:,3);
         dpw_pot(:,7) = dpw_pot(:,1).*dpw_Td_idx(:,4);
         % initialize first state
         sdpw1(:,1)=10*pelect(n,1).*dpw_pot(:,3).*dpw_pot(:,2);
         %initialize all filter states 
         sdpw2(:,1)=sdpw1(:,1).*(ones(n_dpw,1)-dpw_pot(:,4));
         var1 = sdpw1(:,1).*dpw_pot(:,4);
         sdpw3(:,1)=(var1+sdpw2(:,1)).*(ones(n_dpw,1)-dpw_pot(:,5));
         var2 = (var1 + sdpw2(:,1)).*dpw_pot(:,5);
         sdpw4(:,1) = (var2 +sdpw3(:,1)).*(ones(n_dpw,1)-dpw_pot(:,6));
         var3 = (var2 + sdpw3(:,1)).*dpw_pot(:,6);
         sdpw5(:,1)=(var3 +sdpw4(:,1)).*(ones(n_dpw,1)-dpw_pot(:,7));
         var4 = (var3 + sdpw4(:,1)).*dpw_pot(:,6);
         sdpw6(:,1)= var4 + sdpw5(:,1);
         dpw_out(dpw_pss_idx,1)= sdpw6(:,1)-sdpw1(:,1);
      end
   end
   
   if flag == 1 % network interface computation
      % no interface calculation is required
   end
   
   if flag == 2 % deltP/w filter dynamics calculation
      if i ~= 0 % scalar computation
         error('vector computation only for dpwf')
      else 
         % vector computation
         if n_dpw~=0
            n = dpw_mb_idx; % machine number vector       
            var1 = mac_spd(n,k)-ones(n_dpw,1) +sdpw1(:,k);
            var2 = var1.*dpw_pot(:,4) + sdpw2(:,k);
            var3 = var2*dpw_pot(:,5) + sdpw3(:,k);
            var4 = var3*dpw_pot(:,6) + sdpw4(:,k);
            var5 = var4*dpw_pot(:,7) + sdpw5(:,k);
            % integrator for power input 10/(1+10s)
            dsdpw1(:,k) = (-sdpw1(:,k)/10 + pelect(n,k).*dpw_pot(:,3).*dpw_pot(:,2));
            % filter state rate of change
            dsdpw2(:,k) = (var1.*(ones(n_dpw,1)-dpw_pot(:,4))-sdpw2(:,k))./dpw_con(:,4);
            dsdpw3(:,k) = (var2.*(ones(n_dpw,1)-dpw_pot(:,5))-sdpw3(:,k))./dpw_con(:,4);
            dsdpw4(:,k) = (var3.*(ones(n_dpw,1)-dpw_pot(:,6))-sdpw4(:,k))./dpw_con(:,4);
            dsdpw5(:,k) = (var4.*(ones(n_dpw,1)-dpw_pot(:,7))- sdpw5(:,k))./dpw_con(:,4);
            dsdpw6(:,k) = (var5-sdpw6(:,k))./dpw_con(:,4);
            dpw_out(dpw_pss_idx,k) = sdpw6(:,k)-sdpw1(:,k);   
         end
      end
   end
end
