function f = tcsc(i,k,bus,flag)
% Syntax: f = tcsc(i,k,bus,flag)
% 29/12/98
% Purpose: thyristor controlled series capacitor, 
%          with vectorized computation option
%          NOTE - TCSC buses must be declared as 
%                 non-conforming load buses
%                 The initial capacitance must be inserted as a loss free negative reactance 
%                 between the two tcsc buses.
%                 These buses are inserted in the compensated line 
% Input: i - tcsc number: if i= 0, vectorized computation
%        	  k - integer time
%        	  bus - solved loadflow bus data
%             flag - 0 - initialization
%                    1 - network interface computation
%                    2 - generator dynamics computation
%       
%
% Output: f - a dummy variable
%                    



% (c) Copyright 1998 Joe H. Chow/Cherry Tree Scientific Software
%     All Rights Reserved

% History (in reverse chronological order)
%
% Version 1
% Date: December 1998
% Author: Graham Rogers
% system variables
global  basmva bus_int

% tcsc variables
global  tcsc_con n_tcsc tcsvf_idx tcsct_idx 
global  B_tcsc dB_tcsc 
global  tcsc_sig tcsc_dsig


jay = sqrt(-1);
f=0;
if ~isempty(tcsc_con)
   if flag == 0; % initialization
      if i~=0
         B_tcsc(i,1)=0;dB_tcsc(i,1)=0;
      else % vectorized calculation
         B_tcsc(:,1)=zeros(n_tcsc,1);
         dB_tcsc(:,1)=zeros(n_tcsc,1);
      end
   end
   if flag == 1 % network interface computation
      % no interface calculation required - done in nc_load
   end  
   if flag == 2 % tcsc dynamics calculation
      if i ~= 0
         err =  tcsc_sig(i,k) + tcsc_dsig(i,k);
         dB_cv(i,k) = (-B_tcsc(i,k)+tcsc_con(i,4)*err)/tcsc_con(i,5);
         % anti-windup reset
         if B_tcsc(i,k) > tcsc_con(i,6)
            B_tcsc(i,k) = tcsc_con(i,6);
            if dB_tcsc(i,k)>0
               dB_tcsc(i,k) = 0;
            end
         end
         if B_tcsc(i,k) < tcsc_con(i,7)
            B_tcsc(i,k)= tcsc_con(i,7);
            if dB_tcsc(i,k)<0
               dB_tcsc(i,k) = 0;
            end
         end
      else %vectorized computation
         err =  tcsc_sig(:,k) + tcsc_dsig(:,k);
         dB_tcsc(:,k) = (-B_tcsc(:,k)+tcsc_con(:,4).*err)./tcsc_con(:,5);
         % anti-windup reset
         indmx =find( B_tcsc(:,k) > tcsc_con(:,6));
         if ~isempty(indmx)
            B_tcsc(indmx,k) = tcsc_con(indmx,6);
            indrate = find(dB_tcsc(indmx,k)>0);
            if ~isempty(indrate)
               % set rate to zero
               dB_tcsc(indmx(indrate),k) = zeros(length(indrate),1);
            end
         end
         indmn = find(B_tcsc(:,k) < tcsc_con(:,7));
         if ~isempty(indmn)
            B_tcsc(indmn,k) = tcsc_con(indmn,7);
            indrate = find(dB_tcsc(indmn,k)<0);
            if ~isempty(indrate)
               % set rate to zero
               dB_tcsc(indmn(indrate),k) = zeros(length(indrate),1);
            end
         end
      end
   end
end
