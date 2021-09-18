function [f] = smppi(i,k,bus,flag)
% Syntax: [f] = smpexc(i,k,bus,flag) 
% 8:51 am January 29, 1999
%
% Purpose: simple excitation system with pi avr, (exc_con(i,1)=4)
%            with vectorized computation option
%           
% Input: i - generator number
%            0, vectorized computation
%        k - integer time
%        bus - solved loadflow bus data
%        flag - 0 - initialization
%               1 - network interface computation
%               2 - generator dynamics computation
%        
%
% Output: f - dummy variable 
%

% Files:
%
% See Also: smpexc,exc_dc12, exc_st3

% Algorithm:
%
% Calls:
%
% Call By:

% (c) Copyright 1999 Joe H. Chow/Cherry Tree Scientific Software - All Rights Reserved

% Version:1.0
% Date:September 1999
% Author:Graham Rogers

% system variables
global  psi_re psi_im cur_re cur_im

% synchronous machine variables
global  vex eterm mac_int

% excitation system variables
global  exc_con exc_pot Efd V_R V_A V_As R_f V_FB V_TR V_B
global  dEfd dV_R dV_As dR_f dV_TR
global  exc_sig
global  smppi_idx n_smppi ;
global  smppi_TR smppi_TR_idx  smppi_noTR_idx;

% pss variables
global pss_out 
if i ~= 0
   if exc_con(i,1) ~= 4
      error('SMPPI: inappropriate exciter model')
   end
end
f = 0;

[nexc dum] =size(exc_con);
jay = sqrt(-1);
if flag == 0; % initialization
   if i ~= 0  % scalar computation
      n = mac_int(exc_con(i,2)); % machine number
      Efd(i,1) = vex(n,1);
      V_As(i,1) = Efd(i,1); % integrator state
      V_TR(i,1)=eterm(n,1); % input filter state
      err = 0; % summing junction error
      exc_pot(i,3) = eterm(n,1); % reference voltage
   else  % vectorized computation
      if n_smppi ~= 0
         n = mac_int(exc_con(smppi_idx,2)); % machine number with simple exciters
         Efd(smppi_idx,1) = vex(n,1);
         V_As(smppi_idx,1) = Efd(smppi_idx,1); % integrator
         V_TR(smppi_idx,1) = eterm(n,1); % input filter state
         exc_pot(smppi_idx,3) = eterm(n,1); % reference voltage
      end
   end
end

if flag == 1 % network interface computation
   if i ~= 0 % scalar computation
      n = mac_int(exc_con(i,2)); % machine number
      vex(n,k) = Efd(i,k); % field voltage for machines
   else      % vectorized computation
      if n_smppi ~=0 % check for any simple pi exciters
         n = mac_int(exc_con(smppi_idx,2)); % machine numbers for simple exciters
         vex(n,k) = Efd(smppi_idx,k);
      end
   end
end

if flag == 2 % exciter dynamics calculation
   if i ~= 0 % scalar computation
      n = mac_int(exc_con(i,2)); % machine number
      if exc_con(i,3)== 0 %no input filter
         dV_TR(i,k) = 0;
         V_TR(i,k)=eterm(n,k);
      else
         dV_TR(i,k)=(eterm(n,k)-V_TR(i,k))/exc_con(i,3);
      end
      err = exc_sig(i,k)+exc_pot(i,3)-V_TR(i,k)...
         + pss_out(i,k);
      dV_As(i,k) = err*exc_con(i,4);
      dEfd(i,k) = (-Efd(i,k)+V_A(i,k)+exc_con(i,6)*err)...
         /exc_con(i,5);
      % anti-windup reset
      if Efd(i,k) > exc_con(i,8)
         Efd(i,k) = exc_con(i,8);
         if dEfd(i,k)>0
            dEfd(i,k) = 0;
         end
      end
      if Efd(i,k) < exc_con(i,9)
         Efd(i,k) = exc_con(i,9);
         if dEfd(i,k) < 0
            dEfd(i,k) =0;
         end
      end
      R_f(i,k) = 0; dR_f(i,k) = 0;
      V_R(i,k) = 0; dV_R(i,k) = 0;
  
      
   else % vectorized computation
      
      if n_smppi~=0
         % machine numbers
         n = mac_int(exc_con(smppi_idx,2));
         TR = smppi_TR_idx;
         no_TR = smppi_noTR_idx;
         if ~isempty(no_TR) % some exciters have zero TR
            n_nTR = n(no_TR);
            dV_TR(smppi_idx(no_TR),k)=zeros(length(no_TR),1);
            V_TR(smppi_idx(no_TR),k)=eterm(n_nTR,k);
         end
         if ~isempty(TR) %some exciters have nonzero TR
            n_TR = n(TR);
            dV_TR(smppi_idx(TR),k)=(eterm(n_TR,k)-V_TR(smppi_idx(TR),k))...
               ./exc_con(smppi_idx(TR),3);
         end
         % error defined for all simple exciters
         err = exc_sig(smppi_idx,k)+exc_pot(smppi_idx,3)-V_TR(smppi_idx,k)...
            + pss_out(smppi_idx,k);
         dV_As(smppi_idx,k) = err*exc_con(smppi_idx,4);
         dEfd(smppi_idx,k) = (-Efd(smppi_idx,k)+V_As(smppi_idx,k)+exc_con(smppi_idx,6)*err)...
            ./exc_con(smppi_idx,5);
         TA_max = find(Efd(smppi_idx,k)>exc_con(smppi_idx,8));
         TA_min = find(Efd(smppi_idx,k)<exc_con(smppi_idx,9));
         if ~isempty(TA_max)
            Efd(smppi_idx(TA_max),k) = exc_con(smppi_idx(TA_max),8);
            dTA_max = find(dEfd(smppi_idx(TA_max),k)>0);
            if ~isempty(dTA_max)
               n_dTA = length(dTA_max);
               dEFD(smppi_idx(TA_max(dTA_max)),k) = zeros(n_dTA,1);
            end
         end
         if ~isempty(TA_min)
            Efd(smppi_idx(TA_min),k) = exc_con(smppi_idx(TA_min),9);
            dTA_min = find(dEfd(smppi_idx(TA_min),k)<0);
            if ~isempty(dTA_min)
               n_dTA = length(dTA_min);
               dEFD(smppi_idx(TA_min(dTA_min)),k) = zeros(n_dTA,1);
            end
         end
         R_f(smppi_idx,k) = zeros(n_smppi,1); dR_f(smppi_idx,k) = zeros(n_smppi,1);
         V_R(smppi_idx,k) =zeros(n_smppi,1); dV_R(smppi_idx,k) = zeros(n_smppi,1);
      end
   end
end
