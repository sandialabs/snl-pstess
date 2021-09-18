function f = lmod(i,k,bus,flag)
% Syntax: f = lmod(i,k,bus,flag)
% 4:58 PM 15/08/97
% Purpose: load modulation model 
%          with vectorized computation option
%          NOTE - load modulation bus must be declared as a
%                 non-conforming load bus
% Input: i - load modulation number
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
global  lmod_con n_lmod lmod_idx
global lmod_pot lmod_st dlmod_st
global  lmod_sig

f = 0;% dummy variable


jay = sqrt(-1);
if ~isempty(lmod_con)
   if flag == 0; % initialization
      if i~=0
         lmod_pot(i,1) = lmod_con(i,4)*lmod_con(i,3)/basmva;
         % max modulation on system base
         lmod_pot(i,2) = lmod_con(i,5)*lmod_con(i,3)/basmva;
         % min modulation on system base
         lmod_st(i,1) = 0;
      else % vectorized calculation
         lmod_pot(:,1) = lmod_con(:,4).*lmod_con(:,3)/basmva;
         % max modulation on system base
         lmod_pot(:,2) = lmod_con(:,5).*lmod_con(:,3)/basmva;
         % min modulation on system base
         lmod_st(:,1) = zeros(n_lmod,1);
      end
   end
   if flag == 1 % network interface computation
      % no interface calculation required - done in nc_load
   end
   
   if flag == 2 %  dynamics calculation
      % for linearization with operating condition at limits,
      % additional code will be needed
      if i ~= 0
         dlmod_st(i,k) = (-lmod_st(i,k)+lmod_con(i,6)*lmod_sig(i,k))/lmod_con(i,7);
         % anti-windup reset
         if lmod_st(i,k) > lmod_pot(i,1)
            if dlmod_st(i,k)>0
               dlmod_st(i,k) = 0;
            end
         end
         if lmod_st(i,k) < lmod_pot(i,2)
            if dlmod_st(i,k)<0
               dmod_st(i,k) = 0;
            end
         end
      else %vectorized computation
         dlmod_st(:,k) = (-lmod_st(:,k)+lmod_con(:,6).*lmod_sig(:,k))./lmod_con(:,7);
         % anti-windup reset
         indmx =find( lmod_st(:,k) > lmod_pot(:,1));
         if ~isempty(indmx)
            lmod_st(indmx,k) = lmod_pot(indmx,1);
            indrate = find(dlmod_st(indmx,k)>0);
            if ~isempty(indrate)
               % set rate to zero
               dlmod_st(indmx(indrate),k) = zeros(length(indrate),1);
            end
         end
         indmn = find(lmod_st(:,k) < lmod_pot(:,2));
         if ~isempty(indmn)
            l_mod_st(indmn,k) = lmod_pot(indmn,2);
            indrate = find(dlmod_st(indmn)<0);
            if ~isempty(indrate)
               % set rate to zero
               dlmod_st(indmn(indrate),k) = zeros(length(indrate),1);
            end
         end
      end
   end
end
