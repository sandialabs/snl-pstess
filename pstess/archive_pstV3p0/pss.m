function [f] = pss(i,k,bus,flag)
% Syntax: [f] = pss(i,k,bus,flag)
% 8:05 am 23/10/97
%
% Purpose: power system stabilization model
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
% Files:
%
% See Also: exc_dc12, exc_st3

% Algorithm:
%
% Calls:
%
% Called By:

% (c) Copyright 1991-2 Joe H. Chow - All Rights Reserved

% History (in reverse chronological order)
% Version 3.1
% April 2011
% corrected washout filter calculation 
%
% Version 2.2
% July 1998
% added interface to deltaP/omega filter
%
% Version 2.1
% October 1997
% changed length to isempty in index checks
%
% Version:  2.0
% Date:     June 1996
% Author:   Graham Rogers
% Purpose:  To allow vector calculation with pss on only some units
% Modification: modified vector code
%               added negative output limit

% Version:  1.0
% Author:   Joe H. Chow
% Date:     April 1992

% pss variables
global  pss_con pss_pot pss_mb_idx pss_exc_idx
global  pss1 pss2 pss3 dpss1 dpss2 dpss3 pss_out
global  pss_idx n_pss pss_sp_idx pss_p_idx;
global  pss_T  pss_T2 pss_T4 pss_T4_idx  pss_noT4_idx;
global  dpw_pss_idx n_dpw dpw_out
global  mac_con mac_int mac_spd pelect basmva 
f=0;
if n_pss~=0
   if i ~= 0
      if pss_con(i,1) ~= 1 & pss_con(i,1) ~= 2 
         error('PSS: inappropriate power system stablizer model')
      end
   end
   
   
   
   jay = sqrt(-1);
   num_mac=length(mac_con(:,1));
   
   if flag == 0; % initialization
      if i ~= 0  % scalar computation
         n = pss_mb_idx(i); % machine number
         if pss_con(i,1) == 1
            pss1(i,1) = mac_spd(n,1);
         else
            pss1(i,1) = pelect(n,1)*basmva/mac_con(n,3);
         end
         if n_dpw ~= 0
            i_dpw = find(dpw_pss_idx==i);
            if ~isempty(i_dpw)
               pss1(i,1)= dpw_out(i_dpw,1);
            end
         end
         pss2(i,1) = 0.;
         pss3(i,1) = 0.;
         pss_out(pss_exc_idx(i),1) = 0.;
         pss_pot(i,1) = pss_con(i,5)/pss_con(i,6);
         pss_pot(i,2) = 1.0;
         if pss_con(i,8) ~= 0
            pss_pot(i,2) = pss_con(i,7)/pss_con(i,8);
         end
      else
         
         %vector computation
         pss_pot=ones(n_pss,2);
         n=pss_mb_idx;
         if ~isempty(pss_sp_idx)
            n_sp = mac_int(pss_con(pss_sp_idx,2));
            pss1(pss_sp_idx,1)=mac_spd(n_sp,1);
         end
         if ~isempty(pss_p_idx)
            n_p = mac_int(pss_con(pss_p_idx,2));
            pss1(pss_p_idx,1)=pelect(n_p,1)*basmva./mac_con(n_p,3);
         end
         if n_dpw ~=0
            pss1(dpw_pss_idx,1) = dpw_out(:,1);
         end
         pss2(pss_idx,1)=zeros(n_pss,1);
         pss3(pss_idx,1)=zeros(n_pss,1);
         pss_out(pss_exc_idx,1)=zeros(n_pss,1);
         pss_pot(:,1)=pss_con(:,5)./pss_con(:,6);
         if ~isempty(pss_T4_idx)
            pss_pot(pss_T4_idx,2)=pss_con(pss_T4_idx,7)./pss_T4(pss_T4_idx);
         end
      end
   end
   
   if flag == 1 % network interface computation
      if i ~= 0 % scalar computation
         n = pss_mb_idx(i); % machine number
         if pss_con(i,1) == 1
            var1 = mac_spd(i,k)-pss1(i,k); % do not divide by pss_con(i,4)  JHC April 17, 2011
         else
            n = mac_int(pss_con(i,2)); % machine number 
            var1 = pelect(i,k)*basmva/mac_con(n,3)-pss1(i,k); % do not divide by pss_con(i,4)  JHC April 17, 2011
         end
         if n_dpw~=0
            if n_dpw ~= 0
               i_dpw = find(dpw_pss_idx==i);
               if ~isempty(i_dpw)
                  var1 = dpw_out(i_dpw,k)-pss1(i,k); % do not divide by pss_con(i,4)  JHC April 17, 2011
               end
            end
         end  
         var2 = pss_pot(i,1)*pss_con(i,3)*var1 + pss2(i,k);
         
         if pss_con(i,8) == 0
            var3 = var2;
         else
            var3 = pss_pot(i,2)*var2 + pss3(i,k);
         end
         pss_out(pss_exc_idx(i),k) = min(pss_con(i,9),max(var3,-pss_con(i,9)));
      else
         % vector computation
         if n_pss~=0
            n = pss_mb_idx; % machine number vector
            
            var1 = zeros(n_pss,1);
            var2 = var1; var3 = var1;
            if length(pss_sp_idx)~=0
               n_sp = mac_int(pss_con(pss_sp_idx,2));
               var1(pss_sp_idx) = mac_spd(n_sp,k)-pss1(pss_sp_idx,k); % do not divide by pss_con(i,4)  JHC April 17, 2011
            end
            if ~isempty(pss_p_idx)
               n_p = mac_int(pss_con(pss_p_idx,2));
               var1(pss_p_idx) = pelect(n_p,k)*basmva./mac_con(n_p,3)-pss1(pss_p_idx,k); % do not divide by pss_con(i,4)  JHC April 17, 2011
            end
            if n_dpw ~= 0
               var1 = dpw_out(:,k)-pss1(dpw_pss_idx,k); % do not divide by pss_con(i,4)  JHC April 17, 2011
            end
         end  
         
         var2(pss_idx) = pss_pot(pss_idx,1).*(pss_con(pss_idx,3).*var1) + pss2(pss_idx,k);
         var3 = var2;
         
         if ~isempty(pss_T4_idx)
            var3(pss_T4_idx,1) = pss_pot(pss_T4_idx,2).*var2(pss_T4_idx,1)...
               + pss3(pss_T4_idx,k);
            
         end
         
         pss_out(pss_exc_idx,k) = min(pss_con(pss_idx,9),max(var3,pss_con(pss_idx,10)));    
      end
   end
   
   if flag == 2 % pss dynamics calculation
      if i ~= 0 % scalar computation
         n = pss_mb_idx(i); % machine number
         if pss_con(i,1) == 1
            var1 = mac_spd(i,k)-pss1(i,k);  % do not divide by pss_con(i,4)  JHC April 17, 2011
         else
            n = mac_int(pss_con(i,2)); % machine number 
            var1 = pelect(i,k)*basmva./mac_con(n,3)-pss1(i,k); % do not divide by pss_con(i,4)  JHC April 17, 2011
         end
         if n_dpw~=0
            if n_dpw ~= 0
               i_dpw = find(dpw_pss_idx==i);
               if ~isempty(i_dpw)
                  var1 = (dpw_out(i_dpw,k)-pss1(i,k))/pss_con(i,4);
               end
            end
         end  
         
         dpss1(i,k) = var1/pss_con(i,4);   % divide by pss_con(i,4)  JHC April 17, 2011
         
         var2 = pss_pot(i,1)*pss_con(i,3)*var1 + pss2(i,k);
         dpss2(i,k) = ((1-pss_pot(i,1))*pss_con(i,3)*var1 - pss2(i,k))/pss_con(i,6);
         
         if pss_con(i,8) == 0
            var3 = var2;
            dpss3(i,k) = dpss2(i,k);
         else
            var3 = pss_pot(i,2)*var2 + pss3(i,k);
            dpss3(i,k) = ((1-pss_pot(i,2))*var2 - pss3(i,k))/pss_con(i,8);
         end
         pss_out(pss_exc_idx(i),k) = min(pss_con(i,9),max(var3,-pss_con(i,9)));
         
      else
         
         % vector computation
         if n_pss~=0
            n = pss_mb_idx; % machine number vector
            var1 = zeros(n_pss,1);
            var2 = var1; var3 = var1;
            if ~isempty(pss_sp_idx)
               n_sp = mac_int(pss_con(pss_sp_idx,2));
               var1(pss_sp_idx) = mac_spd(n_sp,k)-pss1(pss_sp_idx,k); % do not divide by pss_con(pss_sp_idx,4)  JHC April 17, 2011
            end
            if ~isempty(pss_p_idx)
               n_p = mac_int(pss_con(pss_p_idx,2));
               var1(pss_p_idx) = pelect(n_p,k)*basmva./mac_con(n_p,3)...
                  -pss1(pss_p_idx,k); % do not divide by pss_con(pss_sp_idx,4)  JHC April 17, 2011
            end
            if n_dpw ~= 0
               var1 = (dpw_out(:,k)-pss1(dpw_pss_idx,k))./pss_con(dpw_pss_idx,4);
            end
         end  
         
         dpss1(pss_idx,k) = var1./pss_con(pss_sp_idx,4); % divide by pss_con(pss_sp_idx,4)  JHC April 17, 2011
         
         var2 = pss_pot(pss_idx,1).*(pss_con(pss_idx,3).*var1) + pss2(pss_idx,k);
         dpss2(pss_idx,k) = ((ones(n_pss,1)-pss_pot(pss_idx,1))...
            .*(pss_con(pss_idx,3).*var1 )...
            - pss2(pss_idx,k))./pss_con(pss_idx,6);
         
         var3 = var2;
         dpss3(:,k) = dpss2(:,k);
         if ~isempty(pss_T4_idx)
            var3(pss_T4_idx) = pss_pot(pss_T4_idx,2).*var2(pss_T4_idx)...
               + pss3(pss_T4_idx,k);
            dpss3(pss_T4_idx,k) = ((ones(length(pss_T4_idx),1)...
               -pss_pot(pss_T4_idx,2)).*var2(pss_T4_idx)...
               - pss3(pss_T4_idx,k))./pss_T4(pss_T4_idx);
         end
         pss_out(pss_exc_idx,k) = min(pss_con(pss_idx,9),max(var3,pss_con(pss_idx,10)));    
      end
   end
end   
   