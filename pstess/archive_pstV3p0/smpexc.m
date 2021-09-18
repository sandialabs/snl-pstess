function [f] = smpexc(i,k,bus,flag)
% Syntax: [f] = smpexc(i,k,bus,flag) 
% 8:51 am January 29, 1999
%
% Purpose: simple excitation system, (exc_con(i,1)=0)
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
% See Also: exc_dc12, exc_st3

% Algorithm:
%
% Calls:
%
% Call By:

% (c) Copyright 1991 Joe H. Chow - All Rights Reserved

% History (in reverse chronological order)
% Error correction January 1999
% in non-vector dynamics calculation
% Efd(i,k) = exc_cin(i,9); changed to Efd(i,k) = exc_con(i,9);
% Version:2.1
% Date:August 1997
% Author:Graham Rogers
% Purpose:Reversed change in exc_sig referencing
%         Added pss_out to free exc_sig for other control functions
%
% Version:2.0
% Date:June 1996
% Author:Graham Rogers
% Purpose:Modified vectorization to allow different exciters
% Modification:	note change in way exc_sig is referenced
%               i.e., in terms of the generator number and not the
%               exciter number

% Version:  1.0
% Author:   Joe H. Chow
% Date:     March 1991


% system variables
global  psi_re psi_im cur_re cur_im

% synchronous machine variables
global  vex eterm mac_int

% excitation system variables
global  exc_con exc_pot Efd V_R V_A V_As R_f V_FB V_TR V_B
global  dEfd dV_R dV_As dR_f dV_TR
global  exc_sig
global  smp_idx n_smp ;
global  smp_TA smp_TA_idx smp_noTA_idx s smp_TB smp_TB_idx smp_noTB_idx;
global  smp_TR smp_TR_idx  smp_noTR_idx;

% pss variables
global pss_out 
if i ~= 0
  if exc_con(i,1) ~= 0
    error('SMPEXC: inappropriate exciter model')
  end
end
f = 0;

[nexc dum] =size(exc_con);
jay = sqrt(-1);
if flag == 0; % initialization
   if i ~= 0  % scalar computation
     n = mac_int(exc_con(i,2)); % machine number
     Efd(i,1) = vex(n,1);
     V_A(i,1) = Efd(i,1)/exc_con(i,4); % laglead
     V_As(i,1) = V_A(i,1); % leadlag state variable
     V_TR(i,1)=eterm(n,1); % input filter state
     err = V_A(i,1); % summing junction error
     if exc_con(i,6)~=0
       exc_pot(i,5) = exc_con(i,7)/exc_con(i,6);
     else
       exc_pot(i,5)=1;
     end
     exc_pot(i,3) = eterm(n,1)+err; % reference voltage
    
   else  % vectorized computation

     if n_smp ~= 0
       n = mac_int(exc_con(smp_idx,2)); % machine number with simple exciters
       Efd(smp_idx,1) = vex(n,1);
       V_A(smp_idx,1) = Efd(smp_idx,1)./exc_con(smp_idx,4); % laglead
       V_As(smp_idx,1) = V_A(smp_idx,1); % leadlag state variable
       V_TR(smp_idx,1) = eterm(n,1); % input filter state
       err = V_A(smp_idx,1); % summing junction error
       exc_pot(smp_idx,5) = ones(n_smp,1);
       TB = smp_TB_idx;
       if length(TB)~=0
         exc_pot(smp_idx(TB),5) = exc_con(smp_idx(TB),7)./exc_con(smp_idx(TB),6);
       end
       exc_pot(smp_idx,3) = eterm(n,1)+err; % reference voltage
     end
   end
end

if flag == 1 % network interface computation
   if i ~= 0 % scalar computation
     n = mac_int(exc_con(i,2)); % machine number
     if exc_con(i,5)~=0
       vex(n,k) = Efd(i,k); % field voltage for machines
     else
       if exc_con(i,6)==0
         vex(n,k)=(exc_sig(i,k)+exc_pot(i,3)-V_TR(i,k)...
                   + pss_out(i,k))*exc_con(i,4);
       else
         V_A(i,k)=(exc_sig(i,k)+exc_pot(i,3)-V_TR(i,k)+pss_out(i,k))...
                  *exc_pot(i,5) + (1.-exc_pot(i,5))*V_As(i,k);
         vex(n,k) = V_A(i,k)*exc_con(i,4);
       end
     end 
   else      % vectorized computation

     if n_smp ~=0 % check for any simple exciters
        n = mac_int(exc_con(smp_idx,2)); % machine numbers for simple exciters
        % field voltage for machines having TA and TB  zero
        not_TATB = find((smp_TA + smp_TB)<0.001);
        if length(not_TATB) ~=0
          n_nTATB = n(not_TATB);% machine numbers for exciters with TA & TB zero
          vex(n_nTATB,k) = (exc_sig(smp_idx(not_TATB),k)...
                           +exc_pot(smp_idx(not_TATB),3)-...
                           V_TR(smp_idx(not_TATB),k) + ...
                           pss_out(smp_idx(not_TATB),k))...
                           .*exc_con(smp_idx(not_TATB),4);
        end
       
 % field voltage for machines with TB non zero and TA zero      
        TB = find((smp_TA<0.001)&(smp_TB>0.001));
        if ~isempty(TB)   
           n_TB = n(TB);
           V_A(smp_idx(TB),k) = (exc_sig(smp_idx(TB),k)...
                                +exc_pot(smp_idx(TB),3)...
                                -V_TR(smp_idx(TB),k)...
                                +pss_out(smp_idx(TB),k))...
                                .*exc_pot(smp_idx(TB),5) + ...
                                (ones(length(TB),1)-...
                                exc_pot(smp_idx(TB),5)).*V_As(smp_idx(TB),k);
           vex(n_TB,k) = V_A(smp_idx(TB),k).*exc_con(smp_idx(TB),4);
        end 

        % field voltage for non zero TA
        TA=smp_TA_idx;
        if ~isempty(TA)
           n_TA = n(TA); % machine number 
           vex(n_TA,k) = Efd(smp_idx(TA),k); % field voltage for machines with TA
        end 
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
     if exc_con(i,6) == 0 % no leadlag
        dV_As(i,k) = 0;
        V_As(i,k) = err;
        V_A(i,k) = err;
     else
        dV_As(i,k) = (-V_As(i,k)+err)/exc_con(i,6);
        V_A(i,k) = exc_pot(i,5)*err + ...
                   (1-exc_pot(i,5))*V_As(i,k);
     end
     if exc_con(i,5) == 0 % Ta zero
        dEfd(i,k) = 0;
        Efd(i,k) = exc_con(i,4)*V_A(i,k);
        Efd(i,k) = max(exc_con(i,9),...
                   min(Efd(i,k),exc_con(i,8))); % voltage limit
     else
        dEfd(i,k) = (-Efd(i,k)+exc_con(i,4)*V_A(i,k))...
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
     end
     R_f(i,k) = 0; dR_f(i,k) = 0;
     V_R(i,k) = 0; dV_R(i,k) = 0;

  else % vectorized computation
  
    if n_smp~=0
         % machine numbers
         n = mac_int(exc_con(smp_idx,2));
         TR = smp_TR_idx;
         no_TR = smp_noTR_idx;
         if ~isempty(no_TR) % some exciters have zero TR
           n_nTR = n(no_TR);
           dV_TR(smp_idx(no_TR),k)=zeros(length(no_TR),1);
           V_TR(smp_idx(no_TR),k)=eterm(n_nTR,k);
         end
         if ~isempty(TR) %some exciters have nonzero TR
           n_TR = n(TR);
           dV_TR(smp_idx(TR),k)=(eterm(n_TR,k)-V_TR(smp_idx(TR),k))...
                              ./exc_con(smp_idx(TR),3);
         end
         % error defined for all simple exciters
         err = exc_sig(smp_idx,k)+exc_pot(smp_idx,3)-V_TR(smp_idx,k)...
               + pss_out(smp_idx,k);
         no_TB = smp_noTB_idx;
         if ~isempty(no_TB)
           dV_As(smp_idx(no_TB),k) = zeros(length(no_TB),1);
           V_As(smp_idx(no_TB),k) = err(no_TB);
           V_A(smp_idx(no_TB),k) = err(no_TB);
         end
         TB = smp_TB_idx;
         if ~isempty(TB)
           dV_As(smp_idx(TB),k) = (-V_As(smp_idx(TB),k)+err(TB))...
                                 ./exc_con(smp_idx(TB),6);
           V_A(smp_idx(TB),k) = exc_pot(smp_idx(TB),5).*err(TB)+ ...
                     (ones(length(TB),1)-exc_pot(smp_idx(TB),5)).*V_As(smp_idx(TB),k);
         end
         no_TA = smp_noTA_idx;
         if ~isempty(no_TA)
           dEfd(smp_idx(no_TA),k) = zeros(length(no_TA),1);
           Efd(smp_idx(no_TA),k) = exc_con(smp_idx(no_TA),4)...
                                  .*V_A(smp_idx(no_TA),k);
           % apply output limits
           Efd(smp_idx(no_TA),k) = min(exc_con(smp_idx(no_TA),8)...
                                  ,max(Efd(smp_idx(no_TA),k)...
                                  ,exc_con(smp_idx(no_TA),9)));
         end
         TA = smp_TA_idx;
         if ~isempty(TA)
           dEfd(smp_idx(TA),k) = (-Efd(smp_idx(TA),k) + exc_con(smp_idx(TA),4)...
                               .*V_A(smp_idx(TA),k))...
                               ./exc_con(smp_idx(TA),5);
           %apply non-windup limits
           TA_max = find(Efd(smp_idx(TA),k)>exc_con(smp_idx(TA),8));
           TA_min = find(Efd(smp_idx(TA),k)<exc_con(smp_idx(TA),9));
           if ~isempty(TA_max)
             Efd(smp_idx(TA(TA_max)),k) = exc_con(smp_idx(TA(TA_max)),8);
             dTA_max = find(dEfd(smp_idx(TA(TA_max)),k)>0);
             if ~isempty(dTA_max)
               n_dTA = length(dTA_max);
               dEFD(smp_idx(TA(TA_max(dTA_max))),k) = zeros(n_dTA,1);
             end
           end
           if ~isempty(TA_min)
             Efd(smp_idx(TA(TA_min)),k) = exc_con(smp_idx(TA(TA_min)),9);
             dTA_min = find(dEfd(smp_idx(TA(TA_min)),k)<0);
             if ~isempty(dTA_min)
               n_dTA = length(dTA_min);
               dEFD(smp_idx(TA(TA_min(dTA_min))),k) = zeros(n_dTA,1);
             end
           end
         end
         R_f(smp_idx,k) = zeros(n_smp,1); dR_f(smp_idx,k) = zeros(n_smp,1);
         V_R(smp_idx,k) =zeros(n_smp,1); dV_R(smp_idx,k) = zeros(n_smp,1);
    end
  end
end
