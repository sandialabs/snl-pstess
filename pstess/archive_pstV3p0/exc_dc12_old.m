function [f] = exc_dc12(i,k,bus,flag)
% Syntax: [f] = exc_dc12(i,k,bus,flag)  
% 9:07 am 2/7/98
% Purpose: excitation system, model DC1 (exc_con(i,1)=1)
%            and model DC2 (exc_con(i,1)=2)
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
% See Also: exc_st3 smp_exc
%
% (c) Copyright 1991-1996 Joe H. Chow - All Rights Reserved
% History (in reverse chronological order)
%
% Corrected error with no_TE - GJR
% Date:         2/7/98
% Version:      2.2
% Date:         30/6/98
% Purpose:      Corrected satuaration Asat modified
% Version:	2.1
% Date:		14/8/97
% Author:	Graham Rogers
% Purpose:	Reverted to exciter number for exc_sig
%               Added pss_out so that exc_sig can be used for other
%               control inputs if desired
%               Corrected sign of exc_sig in error signal
% Version:	2.0
% Date:		22/6/96
% Author:	Graham Rogers
% Purpose:	To enable vector calculations with different
%               exciter models on generators
%               This means that the vector option (i=0) 
%               is the mode normally used
%               Indexes formed in exc_indx are required
% Modification:	changed way exc_sig is referenced
%               i.e., in terms of the generator number and not the
%               exciter number
% Version:	1.1
% Date:		22/5/95
% Author:	GJR
% Purpose:	Correct errors	
% Modification:	Saturation model and Lead/Lag initialization

% Version:  1.0
% Author:   Joe H. Chow
% Date:     March 1991

% excitation system variables
global  exc_con exc_pot Efd V_R V_A V_As R_f V_FB V_TR V_B
global  dEfd dV_R dV_As dR_f dV_TR

% indexing variables
global  dc_idx n_dc dc2_idx n_dc2 ;
global dc_TA dc_TA_idx dc_noTA_idx dc_TB dc_TB_idx dc_noTB_idx;
global dc_TE dc_TE_idx dc_noTE_idx;
global dc_TF dc_TF_idx  dc_TR dc_TR_idx dc_noTR_idx;

% pss variables

global pss_out 
%interface variable
global  exc_sig 

% machine variables
global  mac_int vex eterm

if i ~= 0
  if exc_con(i,1) ~= 1 | exc_con(i,1) ~= 2 
    disp('EXC_DC12: inappropriate exciter model')
    stop
  end
end
f = 0;

[nexc dum] =size(exc_con);
jay = sqrt(-1);
if flag == 0; % initialization
  if i ~= 0  % scalar computation
    n = mac_int(exc_con(i,2)); % machine number
    Efd(i,1) = vex(n,1);
    if exc_con(i,11) == 0,  % T_E = 0
        V_R(i,1) = vex(n,1); 
        exc_pot(i,1) = 0; exc_pot(i,2) = 0;
    else
        if exc_con(i,14) == 0  % Asat and Bsat specified
            exc_pot(i,1) = exc_con(i,13); exc_pot(i,2) = exc_con(i,15);
        else
            exc_pot(i,2) =log(exc_con(i,15)*exc_con(i,14)/exc_con(i,12)/exc_con(i,13))...
                          /(exc_con(i,14)-exc_con(i,12)); %B
            exc_pot(i,1) = exc_con(i,13)/exp(exc_pot(i,2)...
                          *exc_con(i,12));                %A
        end
        if exc_con(i,8) == 0 % Vrmax = 0 and compute Vrmax 
                             % assuming E2 is max Efd
          exc_con(i,8) = exc_con(i,10)*exc_con(i,14) + ...
                         exc_pot(i,1)*exp(exc_pot(i,2)...
                         *abs(exc_con(i,14)))*sign(exc_con(i,14));
          if exc_con(i,1) == 2  % bus fed voltage limit
            exc_con(i,8) = exc_con(i,8)/eterm(n,1);
          end
          exc_con(i,9) = -exc_con(i,8);
        end
    end
    if exc_con(i,10) == 0 % KE = 0, calculate KE to make V_R zero
        exc_con(i,10) = - exc_pot(i,1)*exp(exc_pot(i,2)...
                          *abs(Efd(i,1)))/abs(Efd(i,1));
    end
    V_R(i,1) = exc_con(i,10)*Efd(i,1) + exc_pot(i,1)*...
                   exp(exc_pot(i,2)*abs(Efd(i,1)))*sign(Efd(i,1));
    if exc_con(i,1) == 1
      mult = 1;
    else
      mult = eterm(n,1);
    end
    if V_R(i,1) > exc_con(i,8)*mult
         error('EXC_DC12: V_R exceeds upper limit')
    elseif V_R(i,1) < exc_con(i,9)*mult
            error('EXC_DC12: V_R under lower limit')
    end
    if exc_con(i,4) == 0 % KA = 0
      exc_con(i,4) == 1; % reset to 1
    end
    V_A(i,1) = V_R(i,1)/exc_con(i,4); % laglead
    V_As(i,1) = V_A(i,1); % leadlag state variable
    if exc_con(i,6) ~= 0
      exc_pot(i,4) = exc_con(i,7)/exc_con(i,6);
    end
    err = V_A(i,1); % summing junction error
    exc_pot(i,5) = exc_con(i,16)/exc_con(i,17);
    R_f(i,1) = Efd(i,1); % rate feedback state variable
    V_FB(i,1) = 0;
    exc_pot(i,3) = eterm(n,1)+err; % reference voltage
    V_TR(i,1) = eterm(n,1); % transducer state var
    
  else  % vectorized computation
    if n_dc~=0 
       n = mac_int(exc_con(dc_idx,2)); % dc machine numbers
       if n_dc2~=0
        n2 = mac_int(exc_con(dc2_idx,2)); %dc type 2 exciters gen numb
       end
       Efd(dc_idx,1) = vex(n,1);
       TE = dc_TE(dc_TE_idx);
              
       expsat= find(exc_con(dc_idx,14)==0);
       if ~isempty(expsat)
         exc_pot(dc_idx(expsat),1) = exc_con(dc_idx(expsat),13); 
         exc_pot(dc_idx(expsat),2) = exc_con(dc_idx(expsat),15);         
       end
       sesat = find(exc_con(dc_idx,14)~=0);
       if ~isempty(sesat)
         exc_pot(dc_idx(sesat),2)=log(exc_con(dc_idx(sesat),15).*exc_con(dc_idx(sesat),14)...
                                  ./exc_con(dc_idx(sesat),12)...
                                  ./exc_con(dc_idx(sesat),13))...
                                  ./(exc_con(dc_idx(sesat),14)...
                                  -exc_con(dc_idx(sesat),12)); %Bsat
         exc_pot(dc_idx(sesat),1) = exc_con(dc_idx(sesat),13)...
                                    ./exp(exc_pot(dc_idx(sesat),2)...
                                    .*exc_con(dc_idx(sesat),12)); %Asat               
       end
       if ~isempty(dc_noTE_idx)
          no_TE = dc_noTE_idx;
           V_R(dc_idx(no_TE),1) = vex(n(no_TE),1); 
           exc_pot(dc_idx(no_TE),1) = zeros(length(no_TE),1); 
           exc_pot(dc_idx(no_TE),2) = exc_pot(dc_idx(no_TE),1);
       end

       no_Vrmax=find( exc_con(dc_idx,8) == 0); % Vrmax = 0 and compute Vrmax 
                                         % assuming E2 is max Efd
       if ~isempty(no_Vrmax)    
         exc_con(dc_idx(no_Vrmax),8) = exc_con(dc_idx(no_Vrmax),10)...
                                    .*exc_con(dc_idx(no_Vrmax),14) + ...
                                    exc_pot(dc_idx(noVrmax),1)...
                                    .*exp(exc_pot(dc_idx(no_Vrmax),2)...
                                    .*abs(exc_con(dc_idx(no_Vrmax),14)))...
                                    .*sign(exc_con(dc_idx(no_Vrmax),14));
         exc_con(dc_idx(no_Vrmax),9) =  - exc_con(dc_idx(no_Vrmax),8);
         bus_fed = find(exc_con(dc_idx(no_Vrmax),1)==2);
         if ~isempty(bus_fed)
              exc_con(dc_idx(no_Vrmax(bus_fed)),8) = exc_con(dc_idx(no_Vrmax(bus_fed)),8)...
                                                     ./eterm(n(no_Vrmax(bus_fed)),1);
              exc_con(dc_idx(no_Vrmax(bus_fed)),9) = - exc_con(dc_idx(no_Vrmax(bus_fed)),8);
         end
       end
      
       no_KE =find( exc_con(dc_idx,10) == 0); % KE = 0
       if ~isempty(no_KE)  %calculate KE to make initial V_R zero
          asat =exc_pot(dc_idx(no_KE),1);
          bsat = exc_pot(dc_idx(no_KE),2);
          vf = abs(Efd(dc_idx(no_KE),k)); 
          exc_con(dc_idx(no_KE),10) = -asat.*exp(bsat.*vf)./vf;
       end
       if ~isempty(dc_TE_idx)
          R_idx = dc_idx(dc_TE_idx);
          V_R(R_idx,1) = exc_con(R_idx,10).*Efd(R_idx,1)...
                      +exc_pot(R_idx,1).*...
                      exp(exc_pot(R_idx,2).*abs(Efd(R_idx,1)))...
                      .*sign(Efd(R_idx,1));
       end
       %check limits
       mult = ones(n_dc,1);
       if n_dc2~=0
          mult(find(exc_con(dc_idx,1)==2)) = eterm(n2,1);
       end
       over_lmt = find(V_R(dc_idx,1) > exc_con(dc_idx,8).*mult);
       if ~isempty(over_lmt)
         mac_int(exc_con(dc_idx(over_lmt),2))
         error('EXC_DC12: V_R exceeds upper limit')
       end
       under_lmt = find(V_R(dc_idx,1)<exc_con(dc_idx,9).*mult);
       if ~isempty(under_lmt)
           mac_int(exc_con(dc_idx(under_lmt),2))
           error('EXC_DC12: V_R below lower limit')
       end
    
       no_KA = find(exc_con(dc_idx,4) == 0); % KA = 0
       if ~isempty(no_KA)
         exc_con(dc_idx(no_KA),4) == ones(length(no_KA),1); % reset to 1
       end
  
       V_A(dc_idx,1) = V_R(dc_idx,1)./exc_con(dc_idx,4); % laglead
       V_As(dc_idx,1) = V_A(dc_idx,1); % leadlag state variable
       TB = dc_TB_idx;
       if ~isempty(TB)
         exc_pot(dc_idx(TB),4) = exc_con(dc_idx(TB),7)./exc_con(dc_idx(TB),6);
       end
       err = V_A(dc_idx,1); % summing junction error
       exc_pot(dc_idx,5) = exc_con(dc_idx,16)./exc_con(dc_idx,17);
       R_f(dc_idx,1) = Efd(dc_idx,1); % rate feedback state variable
       V_FB(dc_idx,1) = zeros(n_dc,1);
       exc_pot(dc_idx,3) = eterm(n,1)+err; % reference voltage
       V_TR(dc_idx,1) = eterm(n,1); % transducer state var
    end
  end
end

if flag == 1 % network interface computation
   if i ~= 0 % scalar computation
     n = mac_int(exc_con(i,2)); % machine number
     vex(n,k) = Efd(i,k); % set field voltage for machines
   else      % vectorized computation
     if n_dc~=0
       n = mac_int(exc_con(dc_idx,2)); % machine number
       vex(n,k) = Efd(dc_idx,k); % set field voltage for machines
     end
   end
end

if flag == 2 % exciter dynamics calculation
 if i ~= 0 % scalar computation
    n = mac_int(exc_con(i,2)); % machine number
    if exc_con(i,3) == 0  % transducer time constant = 0
      dV_TR(i,k) = 0; 
      V_TR(i,k) = eterm(n,k);
    else 
      dV_TR(i,k) = (-V_TR(i,k)+eterm(n,k))/exc_con(i,3);
    end
    V_FB(i,k) = exc_pot(i,5)*(Efd(i,k)-R_f(i,k));  
    err = exc_sig(i,k) - V_FB(i,k) + exc_pot(i,3) - V_TR(i,k);
    err = err + pss_out(i,k);
    if exc_con(i,6) == 0 % no leadlag
      dV_As(i,k) = 0;
      V_As(i,k) = err;
      V_A(i,k) = err;
    else
      dV_As(i,k) = (-V_As(i,k)+err)/exc_con(i,6);
      V_A(i,k) = exc_pot(i,4)*err + ...
                         (1-exc_pot(i,4))*V_As(i,k);
    end
    if exc_con(i,1) == 1
      mult = 1;  % solid fed
    else
      mult = eterm(n,k);  % bus fed
    end
    if exc_con(i,5) == 0 % no Ta
      dV_R(i,k) = 0.0;
      V_R(i,k) = exc_con(i,4)*V_A(i,k);
      V_R(i,k) = max(exc_con(i,9)*mult,...
           min(V_R(i,k),exc_con(i,8)*mult)); % voltage limit
    else
      dV_R(i,k) = (-V_R(i,k)+exc_con(i,4)*V_A(i,k))...
                  /exc_con(i,5);
      % anti-windup reset
      if V_R(i,k)>exc_con(i,8)*mult
         V_R(i,k) = exc_con(i,8)*mult;
         if dV_R(i,k)>0.0
            dV_R(i,k) = 0.0;
         end
      end
      if V_R(i,k)<exc_con(i,9)*mult
         V_R(i,k) = exc_con(i,9)*mult;
         if dV_R(i,k)<0.0
            dV_R(i,k) = 0.0;
         end
      end
    end

    if exc_con(i,11) == 0 % no exciter dynamics
      dEfd(i,k) = 0.0;
      Efd(i,k) = V_R(i,k);
    else
      SE = exc_pot(i,1)*exp(exc_pot(i,2)*abs(Efd(i,k)))...
               *sign(Efd(i,k));
      dEfd(i,k) = (V_R(i,k) - exc_con(i,10)*Efd(i,k) - SE)...
                  /exc_con(i,11);
    end
    dR_f(i,k) = (-R_f(i,k)+Efd(i,k))/exc_con(i,17); 
                              % rate feedback state variable

 else % vectorized computation
   if n_dc~=0
     n = mac_int(exc_con(dc_idx,2)); % machine number
     if n_dc2~=0
       n2 = mac_int(exc_con(dc2_idx,2));
     end
     TR = dc_TR_idx;
     no_TR = dc_noTR_idx;
     if ~isempty(no_TR)
       dV_TR(dc_idx(no_TR),k) = zeros(length(no_TR),1); 
       V_TR(dc_idx(no_TR),k) = eterm(n(no_TR),k);
     end
     if ~isempty(TR)
       dV_TR(dc_idx(TR),k) = (-V_TR(dc_idx(TR),k)+eterm(n(TR),k))./exc_con(dc_idx(TR),3);
     end
     V_FB(dc_idx,k) = exc_pot(dc_idx,5).*(Efd(dc_idx,k)-R_f(dc_idx,k));  
     err = exc_sig(dc_idx,k) - V_FB(dc_idx,k) + exc_pot(dc_idx,3) - V_TR(dc_idx,k);
     err = err + pss_out(dc_idx,k);
     no_TB = dc_noTB_idx;
     if ~isempty(no_TB)
       dV_As(dc_idx(no_TB),k) = zeros(length(no_TB),1);
       V_As(dc_idx(no_TB),k) = err(no_TB);
       V_A(dc_idx(no_TB),k) = err(no_TB);
     end
     TB= dc_TB_idx;
     if ~isempty(TB)
       dV_As(dc_idx(TB),k) = (-V_As(dc_idx(TB),k)+err(TB))./exc_con(dc_idx(TB),6);
       V_A(dc_idx(TB),k) = exc_pot(dc_idx(TB),4).*err(TB) + ...
                    (ones(length(TB),1)-exc_pot(dc_idx(TB),4))...
                    .*V_As(dc_idx(TB),k);
     end
     mult = ones(n_dc,1); 
     if n_dc2 ~=0
       mult(find(exc_con(dc_idx,1)==2)) = eterm(n2,k);
     end
     no_TA = dc_noTA_idx;
     if ~isempty(no_TA)
       dV_R(dc_idx(no_TA),k) = zeros(length(no_TA),1);
       V_R(dc_idx(no_TA),k) = exc_con(dc_idx(no_TA),4).*V_A(dc_idx(no_TA),k);
       V_R(dc_idx(no_TA),k) = max(exc_con(dc_idx(no_TA),9).*mult(no_TA),...
          min(V_R(dc_idx(no_TA),k),exc_con(dc_idx(no_TA),8).*mult(no_TA))); % voltage limit
     end
     TA = dc_TA_idx;
     if ~isempty(TA)
       dV_R(dc_idx(TA),k) = (-V_R(dc_idx(TA),k)+exc_con(dc_idx(TA),4).*V_A(dc_idx(TA),k))...
                  ./exc_con(dc_idx(TA),5);
       % anti-windup reset
       maxlmt =find( V_R(dc_idx(TA),k)>exc_con(dc_idx(TA),8).*mult(TA));
       if ~isempty(maxlmt)
          V_R(dc_idx(TA(maxlmt)),k) = exc_con(dc_idx(TA(maxlmt)),8).*mult(TA(maxlmt));
         pos_rate = find(dV_R(dc_idx(TA(maxlmt)),k)>0.0);
         prl = length(pos_rate);
         if prl ~=0
           dV_R(dc_idx(TA(maxlmt(pos_rate))),k) = zeros(prl,1);
         end
       end
       minlmt = find(V_R(dc_idx(TA),k)<exc_con(dc_idx(TA),9).*mult(TA));
       if ~isempty(minlmt)
         V_R(dc_idx(TA(minlmt)),k) = exc_con(dc_idx(TA(minlmt)),9).*mult(TA(minlmt));
         neg_rate = find(dV_R(dc_idx(TA(minlmt)),k)<0.0);
         nrl = length(neg_rate);
         if nrl ~=0
           dV_R(dc_idx(TA(minlmt(neg_rate))),k) = zeros(nrl,1);
         end
       end
     end
     no_TE = dc_noTE_idx;
     if ~isempty(no_TE)  
         dEfd(dc_idx(no_TE),k) = zeros(length(no_TE),1);
         Efd(dc_idx(no_TE),k) = V_R(dc_idx(no_TE),k);
     end
     TE = dc_TE_idx;
     if ~isempty(TE)
       SE = exc_pot(dc_idx(TE),1).*exp(exc_pot(dc_idx(TE),2)...
          .*abs(Efd(dc_idx(TE),k)))...
          .*sign(Efd(dc_idx(TE),k));
       dEfd(dc_idx(TE),k) = (V_R(dc_idx(TE),k) - exc_con(dc_idx(TE),10)...
                            .*Efd(dc_idx(TE),k)- SE)./exc_con(dc_idx(TE),11);
     end
     dR_f(dc_idx,k) = (-R_f(dc_idx,k)+Efd(dc_idx,k))./exc_con(dc_idx,17); 
                              % rate feedback state variable
   end
 end
end
