function [f] = tg_hydro(i,k,bus,flag)
% Syntax: [f] = tg_hydro(i,k,bus,flag)  
% 1:19 PM 15/08/97
% Purpose: hydraulic turbine governor model
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
%  1	turbine model number (=2)	
%  2	machine number	
%  3	speed set point   wf		pu
%  4	permanent droop	Rp	pu
%  5  transient droop   Rt pu
%  6	maximum power order  Tmax	pu on generator base
%  7  maximum rate limit pu on gen base/sec
%  8  minimum rate limit pu on gen base per sec
%  9	servo time constant   Ts	sec
%  10  servo gain  Ks
%  11	governor time constant  Tg	sec
%  12	reset time constant Tr	sec
%  13	water starting time    Tw	sec


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
% Version:   1.0
% Date:  June 1998
% Author: Graham Rogers
% Purpose: Model of hydraulic turbine
% Modification:
% (c)copyright Cherry Tree Scientific Software 1998 All Rights Reserved
global  basmva mac_int mac_con pelect pmech mac_spd 
global  tg_con tg_pot 
global  tg1 tg2 tg3 tg4 tg5 dtg1 dtg2 dtg3 dtg4 dtg5
global  tgh_idx n_tgh tg_sig


jay = sqrt(-1);
f = 0;
if n_tgh~=0
   if flag == 0; % initialization
      if i~=0
         error('vector calculation ony in TG_Hydro')
      else
         %  vectorized computation
         if n_tgh~=0
            n = mac_int(tg_con(tgh_idx,2)); % machine number
            maxlmt = find(pmech(n,1)>tg_con(tgh_idx,6));
            if ~isempty(maxlmt)
               disp(' pmech excedes maximum limit at generators')
               disp(n(maxlmt))
               error('stop')   
            end
            minlmt = find(pmech(n,1)<zeros(n_tgh,1));
            if ~isempty(minlmt)
               disp(n(minlmt))
               error('pmech less than zero at some generators ')
            end
            tg1(tgh_idx,1)=zeros(n_tgh,1);
            tg2(tgh_idx,1) = pmech(n,1);
            tg3(tgh_idx,1) = pmech(n,1);
            tg4(tgh_idx,1) = pmech(n,1);
            tg5(tgh_idx,1) = 3*pmech(n,1);
            tg_pot(tgh_idx,5) = pmech(n,1).*tg_con(tgh_idx,4);% reference value
            tg_sig(tgh_idx,1) = zeros(n_tgh,1);
         end
      end
   end
   
   if flag == 1 % network interface computation
      if i ~= 0 % scalar computation
         %vector computation only for tg_hydro
         error('vector computation only for tg_hydro')
      else
         if n_tgh~=0
            n = mac_int(tg_con(tgh_idx,2)); % machine number
            pmech(n,k) = tg5(tgh_idx,k) -2*tg4(tgh_idx,k);
         end
      end
   end
   
   if flag == 2 % turbine governor dynamics calculation
      if i ~= 0 % scalar computation
         % vector computaion only
         error('vector computation only for tg_hydro')
      else
         % vectorized computation
         if n_tgh ~=0
            n = mac_int(tg_con(tgh_idx,2)); % machine number
            spd_err = tg_con(tgh_idx,3) - mac_spd(n,k);
            demand = spd_err + tg_sig(tgh_idx,k)+tg_pot(tgh_idx,5);
            dtg1(tgh_idx,k) = tg_con(tgh_idx,10).*(demand-tg1(tgh_idx,k)-...
               (tg_con(tgh_idx,4)+tg_con(tgh_idx,5)).*tg2(tgh_idx,k)+...
               tg_con(tgh_idx,5).*tg3(tgh_idx,k))./tg_con(tgh_idx,9) ;
            %apply rate limit
            rmax = find(dtg1(tgh_idx,k)>tg_con(tgh_idx,7));
            if ~isempty(rmax)
               dtg1(tgh_idx(rmax),k)=tg_con(tgh_idx(rmax),7);
            end
            rmin = find(dtg1(tgh_idx,k)<tg_con(tgh_idx,8));
            if ~isempty(rmin)
               dtg1(tgh_idx(rmin),k)=tg_con(tgh_idx(rmin),8);
            end
            dtg2(tgh_idx,k) = tg1(tgh_idx,k);
            %check non-wind-up limit
            smax= find(tg2(tgh_idx,k)>=tg_con(tgh_idx,6));
            if ~isempty(smax)
               tg2(tgh_idx(smax),k)=tg_con(tgh_idx(smax),6);
               if dtg2(tgh_idx(smax),k)>0
                  dtg2(tgh_idx(smax),k)=zeros(length(smax),1);
               end
            end
            smin= find(tg2(tgh_idx,k)<=0);
            if ~isempty(smin)
               tg2(tgh_idx(smax),k)=zeros(length(smin),1);
               if dtg2(tgh_idx(smin),k)<0
                  dtg2(tgh_idx(smin),k)=zeros(length(smin),1);
               end
            end
            % transient droop
            dtg3(tgh_idx,k)=(tg2(tgh_idx,k)-tg3(tgh_idx,k))./tg_con(tgh_idx,12);
            % gate servo
            dtg4(tgh_idx,k) = (tg2(tgh_idx,k)-tg4(tgh_idx,k))./tg_con(tgh_idx,11);
            % hydraulic turbine
            dtg5(tgh_idx,k) = 2*(3*tg4(tgh_idx,k)-tg5(tgh_idx,k))./tg_con(tgh_idx,13);
            
            pmech(n,k) = tg5(tgh_idx,k)-2*tg4(tgh_idx,k);
         end
      end
   end
end
