function [bus_new] = svc(i,k,bus,flag,v_sbus)
% Syntax: [bus_new] = svc(i,k,bus,flag,v_sbus)
% 06/05/98
% Purpose: static var system, 
%          with vectorized computation option
%          NOTE - static var bus must be declared as a
%                 non-conforming load bus
% Input: i - static var number
%            if i= 0, vectorized computation
%        k - integer time
%        bus - solved loadflow bus data
%        flag - 0 - initialization
%               1 - network interface computation
%               2 - generator dynamics computation
%        v_sbus - svc bus voltage
%       
%
% Output: bus_new - on initialization bus_new is bus matrix 
%                   with the reactive generation at the 
%                   svc buses set to zero
%                 - otherwise bus_new = bus
%                    



% (c) Copyright 1991-1998 Joe H. Chow/Cherry Tree Scientific Software
%     All Rights Reserved

% History (in reverse chronological order)
%
% Version 1.2
% Date: May 1998
% Author: Graham Rogers
% Purpose: Add lead lag element

% Version:1.1
% Date:July 1995
% Author:Graham Rogers
% Purpose:Add vectorization and correct bugs
% Modification:

% Version:  1.0
% Author:   Joe H. Chow
% Date:     June 1991

% system variables
global  basmva bus_int

% svc variables
global  svc_con n_svc svc_idx svcll_idx
global svc_pot B_cv dB_cv B_con dB_con
global  svc_sig svc_dsig

bus_new  = bus;

jay = sqrt(-1);
if ~isempty(svc_con)
   if flag == 0; % initialization
      if i~=0
         svc_pot(i,1) = svc_con(i,4)*svc_con(i,3)/basmva;
         % B_cv max on system base
         svc_pot(i,2) = svc_con(i,5)*svc_con(i,3)/basmva;
         % B_cv min on system base
         j = bus_int(svc_con(i,2)); % bus number
         B_cv(i,1) = bus(j,5)/(bus(j,2).*bus(j,2)); % initial B_cv based on generation
         bus_new(j,5) = 0;
         if B_cv(i,1) > svc_pot(i,1)
            error('SVC: BCV exceeds maximum at initialization')
         end
         if B_cv(i,1) < svc_pot(i,2)
            error('SVC: BCV below minimum at initialization')
         end
         svc_pot(i,3) = B_cv(i,1); % store initial value of B_cv
         svc_pot(i,4) = bus(j,2) + B_cv(i,1)/svc_con(i,6);
         % reference voltage
         if  svc_con(1,9)~=0
            svc_pot(i,5) = svc_con(i,8)/svc_con(i,9);
         else
            svc_pot(i,5)=1;
         end
         B_con(i,1) = B_cv(i,1)*(1-svc_pot(i,5))/svc_con(i,6);
      else % vectorized calculation
         svc_pot(:,1) = svc_con(:,4).*svc_con(:,3)/basmva;
         % B_cv max on system base
         svc_pot(:,2) = svc_con(:,5).*svc_con(:,3)/basmva;
         % B_cv min on system base
         jsvc = bus_int(svc_con(:,2)); % bus number
         B_cv(:,1) = bus(jsvc,5)./(bus(jsvc,2).*bus(jsvc,2)); % initial B_cv
         bus_new(jsvc,5) = zeros(n_svc,1);
         testmxlmt=max( B_cv(:,1) > svc_pot(:,1));
         if testmxlmt
            error('SVC: BCV exceeds maximum at initialization')
         end
         testmnlmt=max(B_cv(:,1) < svc_pot(:,2));
         if testmnlmt
            error('SVC: BCV below minimum at initialization')
         end
         svc_pot(:,3) = B_cv(:,1); % store initial value of B_cv
         svc_pot(:,4) = bus(jsvc,2) + B_cv(:,1)./svc_con(:,6);% reference voltage  
         svc_pot(:,5) = ones(n_svc,1);
         if ~isempty(svcll_idx)
            svc_pot(svcll_idx,5)= svc_con(svcll_idx,8)./svc_con(svcll_idx,9);
         end
         B_con(:,1) = B_cv(:,1).*(ones(n_svc,1)-svc_pot(:,5))./svc_con(:,6);
      end
   end
   if flag == 1 % network interface computation
      % no interface calculation required - done in nc_load
   end
   
   if flag == 2 % exciter dynamics calculation
      % for linearization with operating condition at limits,
      % additional code will be needed
      if i ~= 0
         err =  svc_sig(i,k) + svc_pot(i,4) + svc_dsig(i,k) - v_sbus;
         if svc_con(i,9)~=0
            dB_con(i,k) = (-B_con(i,k)+err*(1-svc_pot(i,5)))/svc_con(i,9);
         else
            dB_con(i,k)=0;
         end 
         dB_cv(i,k) = (-B_cv(i,k)+svc_con(i,6)*(err*svc_pot(i,5)+B_con(i,k)))/svc_con(i,7);
            % anti-windup reset
            if B_cv(i,k) > svc_pot(i,1)
               if dB_cv(i,k)>0
                  dB_cv(i,k) = 0;
               end
            end
            if B_cv(i,k) < svc_pot(i,2)
               if dB_cv(i,k)<0
                  dB_cv(i,k) = 0;
               end
            end
         else %vectorized computation
            lv_sbus=find(v_sbus<.9&svc_dsig(:,k)<0);
            d_sigin = svc_dsig(:,k);
            if ~isempty(lv_sbus)
               d_sigin(lv_sbus)=zeros(length(lv_sbus),1);
            end
            err =  svc_sig(:,k) + svc_pot(:,4) + d_sigin - v_sbus;
            
            dB_con(:,k)= zeros(n_svc,1);
            if ~isempty(svcll_idx)
               nll = length(svcll_idx);
               dB_con(svcll_idx,k) = (-B_con(svcll_idx,k)+(ones(nll,1)-svc_pot(svcll_idx,5)).*err)./svc_con(svcll_idx,9);
            end
            dB_cv(:,k) = (-B_cv(:,k)+svc_con(:,6).*(err.*svc_pot(:,5)+B_con(:,k)))./svc_con(:,7);
            % anti-windup reset
            indmx =find( B_cv(:,k) > svc_pot(:,1));
            if ~isempty(indmx)
               B_cv(indmx,k) = svc_pot(indmx,1);
               indrate = find(dB_cv(indmx,k)>0);
               if ~isempty(indrate)
                  % set rate to zero
                  dB_cv(indmx(indrate),k) = zeros(length(indrate),1);
               end
            end
            indmn = find(B_cv(:,k) < svc_pot(:,2));
            if ~isempty(indmn)
               B_cv(indmn,k) = svc_pot(indmn,2);
               indrate = find(dB_cv(indmn)<0);
               if ~isempty(indrate)
                  % set rate to zero
                  dB_cv(indmn(indrate),k) = zeros(length(indrate),1);
               end
            end
         end
      end
   end
   