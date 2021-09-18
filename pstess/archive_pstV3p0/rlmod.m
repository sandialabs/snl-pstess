function f = rlmod(i,k,bus,flag)
% Syntax: f = rlmod(i,k,bus,flag)
% 5:20 PM  27/8/97
% Purpose: reactive load modulation model 
%          with vectorized computation option
%          NOTE - the reactive load modulation bus must be declared as a
%                 non-conforming load bus
% Input: i - reactive load modulation number
%            if i= 0, vectorized computation
%        k - integer time
%        bus - solved loadflow bus data
%        flag - 0 - initialization
%               1 - network interface computation
%               2 - generator dynamics computation
%
% Output: f is a dummy variable
%                    


% (c) Copyright 1991-1997 Joe H. Chow/ Cherry Tree Scientific Software - All Rights Reserved

% History (in reverse chronological order)
%
% Version:1
% Date:August 1997
% Author:Graham Rogers

% system variables
global  basmva bus_int

% load modulation variables
global  rlmod_con n_rlmod rlmod_idx
global  rlmod_pot rlmod_st drlmod_st
global  rlmod_sig

f = 0;% dummy variable


jay = sqrt(-1);
if ~isempty(rlmod_con)
   if flag == 0; % initialization
      if i~=0
         rlmod_pot(i,1) = rlmod_con(i,4)*rlmod_con(i,3)/basmva;
         % max modulation on system base
         rlmod_pot(i,2) = rlmod_con(i,5)*rlmod_con(i,3)/basmva;
         % min modulation on system base
         rlmod_st(i,1) = 0;
      else % vectorized calculation
         rlmod_pot(:,1) = rlmod_con(:,4).*rlmod_con(:,3)/basmva;
         % max modulation on system base
         rlmod_pot(:,2) = rlmod_con(:,5).*rlmod_con(:,3)/basmva;
         % min modulation on system base
         rlmod_st(:,1) = zeros(n_rlmod,1);
      end
   end
   if flag == 1 % network interface computation
      % no interface calculation required - done in nc_load
   end
   
   if flag == 2 %  dynamics calculation
      % for linearization with operating condition at limits,
      % additional code will be needed
      if i ~= 0
         drlmod_st(i,k) = (-rlmod_st(i,k)+rlmod_con(i,6)*rlmod_sig(i,k))/rlmod_con(i,7);
         % anti-windup reset
         if rlmod_st(i,k) > rlmod_pot(i,1)
            if drlmod_st(i,k)>0
               drlmod_st(i,k) = 0;
            end
         end
         if rlmod_st(i,k) < rlmod_pot(i,2)
            if drlmod_st(i,k)<0
               drlmod_st(i,k) = 0;
            end
         end
      else %vectorized computation
         drlmod_st(:,k) = (-rlmod_st(:,k)+rlmod_con(:,6).*rlmod_sig(:,k))./rlmod_con(:,7);
         % anti-windup reset
         indmx =find( rlmod_st(:,k) > rlmod_pot(:,1));
         if ~isempty(indmx)
            rlmod_st(indmx,k) = rlmod_pot(indmx,1);
            indrate = find(drlmod_st(indmx,k)>0);
            if ~isempty(indrate)
               % set rate to zero
               drlmod_st(indmx(indrate),k) = zeros(length(indrate),1);
            end
         end
         indmn = find(rlmod_st(:,k) < rlmod_pot(:,2));
         if ~isempty(indmn)
            rl_mod_st(indmn,k) = rlmod_pot(indmn,2);
            indrate = find(drlmod_st(indmn)<0);
            if ~isempty(indrate)
               % set rate to zero
               drlmod_st(indmn(indrate),k) = zeros(length(indrate),1);
            end
         end
      end
   end
end
