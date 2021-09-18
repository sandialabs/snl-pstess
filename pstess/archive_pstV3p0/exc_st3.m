function f = exc_st3(i,k,bus,flag)
% Syntax: f = exc_st3(i,k,bus,flag) 
% 9:42 AM 14/08/97
%
% Purpose: excitation system, model ST3 (exc_con(i,1)=3)
%            Compound source controlled rectifier exciter,
%            with vectorized computation option
%            state variables are: V_TR, V_As, V_R
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
% See Also: exc_dc12, smp_exc, exc_indx
% History (in reverse chronological order)
%
% Version:  3.1
% Date:     August 1997
% Author:   Graham Rogers
% Purpose:  revert to exciter number for exc_sig
%           add pss_out so that exc_sig is available for other control
%           functions
% Version:  3.0
% Date:     June 1996
% Author:   Graham Rogers
% Purpose:  Changes to allow vectorization with different exciters
%           vector option (i=0) is the normal use.
%           Uses indexes formed by exc_indx.
% Modification:	Change in the way exc_sig is referenced
%               i.e., in terms of the generator number and not the
%               exciter number
% (c) copyright Joe Chow 1991-1996. All rights reserved

% Version:  2.0
% Author:   Graham Rogers
% Date:     October 1995
% This version has vector capability
% The svc regulation has been altered to comply
% with IEEE excitation system standard

% Version:  1.0
% Author:   Joe H. Chow
% Date:     April 1991

% system variables
global  psi_re psi_im cur_re cur_im bus_int

% synchronous machine variables
global  mac_con mac_pot mac_int
global  mac_ang mac_spd eqprime edprime psikd psikq
global  curd curq curdg curqg fldcur
global  psidpp psiqpp vex eterm theta ed eq
global  pmech pelect qelect
global  mac_tra_idx  mac_sub_idx n_tra n_sub ;

% excitation system variables
global  exc_con exc_pot Efd V_R V_A V_As R_f V_FB V_TR V_B
global  Efd dEfd dV_R dV_As dR_f dV_TR
global  exc_sig
global  st3_idx n_st3 ;
global  st3_TA st3_TA_idx st3_noTA_idx st3_TB st3_TB_idx st3_noTB_idx;
global  st3_TR st3_TR_idx  st3_noTR_idx;

% pss variables
global pss_out 

f = 0;

jay = sqrt(-1);
if n_st3~=0
  if flag == 0; % initialization
  %check for either mac_sub or mac_tra models
    if (n_sub==0)&&(n_tra==0)
      disp('error in exc_st3')
      error(' you must have either subtransient or transient generator models')
    end
    if i ~= 0; % not vectorized
      if exc_con(i,1) ~= 3
        error('EXC_ST3: inappropriate exciter model')
      end
      n = mac_int(exc_con(i,2)); % machine number
      n_bus = bus_int(mac_con(n,2)); % bus number
      Efd(i,1) = vex(n,1);
      if Efd(i,1) > exc_con(i,18)
        error('EXC_ST3: Efd exceeds maximum in initialization')
      end
      exc_pot(i,1) = exc_con(i,13)*cos(exc_con(i,14)*pi/180);
      exc_pot(i,2) = exc_con(i,13)*sin(exc_con(i,14)*pi/180);
      iterm = (pelect(n,1)-jay*qelect(n,1))/(eterm(n,1)...
               *exp(-jay*theta(n_bus,1)));
      vep = eterm(n,1)*exp(jay*theta(n_bus,1))*(exc_pot(i,1) + jay*exc_pot(i,2));
      ve = vep+jay*(exc_con(i,15)...
           + (exc_pot(i,1)+jay*exc_pot(i,2))*exc_con(i,16))*iterm;
      V_E = abs(ve);
      % calculate regulation of potential source
      I_N = exc_con(i,17)*fldcur(n,1)/V_E;
      if I_N < 0.433
         F_EX = 1-0.5771*I_N;
      elseif I_N < .75
         F_EX = (.75-I_N^2)^.5;
      else
         F_EX = 1.708*(1-I_N);
      end
      if F_EX <=0
	error('I_N zero or negative: check data input')
      end
      V_B(i,1) = V_E*F_EX;
      V_R(i,1) = Efd(i,1)/V_B(i,1);
      if exc_con(i,4) == 0 % KA = 0
        exc_con(i,4) == 1; % reset to 1
      end
      V_A(i,1) = V_R(i,1)/exc_con(i,4)+min(...
                 exc_con(i,20),exc_con(i,19)*Efd(i,1)); % laglead
      V_As(i,1) = V_A(i,1); % leadlag state variable
      if exc_con(i,6) ~= 0 %check for non-zero TB
        exc_pot(i,5) = exc_con(i,7)/exc_con(i,6); %TC/TB
      end
      V_I = V_A(i,1)/exc_con(i,12);
      if V_I > exc_con(i,10)
        error('EXC_ST3: V_I above maximum in initialization')
      elseif V_I < exc_con(i,11)
        error('EXC_ST3: V_I below minimum in initialization')
      end
      exc_pot(i,3) = eterm(n,1)+V_I; % reference voltage
      V_TR(i,1) = eterm(n,1); % transducer state var
      R_f(i,1) = 0;    % zero out unnecessary state variables
  
    else % vector computation
      V_E = zeros(n_st3,1); I_N = V_E; iterm = V_E; vep = V_E;
      ve = V_E;F_EX = V_E;
      n = mac_int(exc_con(st3_idx,2)); % machine number vector
      n_bus = bus_int(mac_con(n,2));%gen bus number vector
      Efd(st3_idx,1) = vex(n,1);
      max_lim = find(Efd(st3_idx,1) > exc_con(st3_idx,18)); %check maximum limit
      if ~isempty(max_lim)
        disp('EXC_ST3: Efd exceeds maximum in initialization at ')
        n_error = mac_int(exc_con(max_lim,2))  % machine number
        error('stop')
      end
      exc_pot(st3_idx,1) = exc_con(st3_idx,13).*cos(exc_con(st3_idx,14)*pi/180);
      exc_pot(st3_idx,2) = exc_con(st3_idx,13).*sin(exc_con(st3_idx,14)*pi/180);
      iterm =(pelect(n,1)-jay*qelect(n,1))./...
             (eterm(n,1).*exp(-jay*theta(n_bus,1))).*mac_pot(n,1);
      vep = eterm(n,1).*exp(jay*theta(n,1)).*(exc_pot(st3_idx,1) + jay*exc_pot(st3_idx,2));
      ve = vep+jay*(exc_con(st3_idx,15)...
           + (exc_pot(st3_idx,1)+jay*exc_pot(st3_idx,2)).*exc_con(st3_idx,16)).*iterm;
      V_E = abs(ve);
 
      % this is the equivalent terminal voltage at the rectifier
      %terminals
      ve_low=find(V_E < 1e-6);
      if ~isempty(ve_low)
        disp('excitation system error:no supply voltage at the following')
        n_excerr = n(ve_low)
        error('stop')
      end

      I_N = exc_con(st3_idx,17).*fldcur(n,1)./V_E;
      % select operating point on the rectifier regulation characteristic
      low_IN=find(I_N < 0.433);
      if ~isempty(low_IN)
        F_EX(low_IN) = ones(length(low_IN),1)-0.5771*I_N(low_IN);
      end

      big_IN=find(I_N > .75);
      if ~isempty(big_IN)
        F_EX(big_IN) = 1.708*(ones(length(big_IN),1)-I_N(big_IN));
      end
      mid_IN =find((I_N > 0.433) & (I_N < .75));
      if ~isempty(mid_IN)
        F_EX(mid_IN) = (0.75*ones(length(mid_IN),1)-(I_N(mid_IN))^.2).^.5;
      end
      fex_error=find(F_EX<=0);
      if ~isempty(fex_error)
	     disp('error: F_EX zero or negative at the following generators')
	     n_error = n(fex_error)
        error('stop')
      end

      V_B(st3_idx,1) = V_E.*F_EX;
      V_R(st3_idx,1) = Efd(st3_idx,1)./V_B(st3_idx,1);

      nKA_idx =find(exc_con(st3_idx,4) == 0); % KA = 0
      if ~isempty(nKA_idx)
        exc_con(st3_idx(nKA_idx),4) = ones(length(nKA_idx),1); % reset to 1
      end

      V_A(st3_idx,1) = V_R(st3_idx,1)./exc_con(st3_idx,4)+min(...
           exc_con(st3_idx,20),exc_con(st3_idx,19).*Efd(st3_idx,1)); % laglead
      V_As(st3_idx,1) = V_A(st3_idx,1); % leadlag state variable
      exc_pot(st3_idx,5) = ones(n_st3,1);

      TB = st3_TB_idx;
      if ~isempty(TB)
        exc_pot(st3_idx(TB),5) = exc_con(st3_idx(TB),7)./exc_con(st3_idx(TB),6);
      end
      V_I = V_A(st3_idx,1)./exc_con(st3_idx,12);
      max_VI=find(V_I > exc_con(st3_idx,10));
      if ~isempty(max_VI)
        disp('EXC_ST3: V_I above maximum in initialization at')
        n_error=mac_int(exc_con(max_VI,2))
        error('stop')
      end

      min_VI = find(V_I < exc_con(st3_idx,11));
      if ~isempty(min_VI)
        disp('EXC_ST3: V_I below minimum in initialization at')
        n_error=mac_int(exc_con(min_VI,2))
        error('stop')
      end

      exc_pot(st3_idx,3) = eterm(n,1)+V_I; % reference voltage
      V_TR(st3_idx,1) = eterm(n,1); % transducer state var
      R_f(st3_idx,1) = zeros(n_st3,1);    % zero out unnecessary state variables
   end
  %end initialization
  end

  if flag == 1 % network interface computation
    if i~=0 %exciter-by-exciter calculation
      n = mac_int(exc_con(i,2)); % machine number
      n_bus = bus_int(mac_con(n,2));
      curd(n,k) = sin(mac_ang(n,k))*cur_re(n,k) - ...
            cos(mac_ang(n,k))*cur_im(n,k); % d-axis current
      curq(n,k) = cos(mac_ang(n,k))*cur_re(n,k) + ...
            sin(mac_ang(n,k))*cur_im(n,k); % q-axis current
      curdg(n,k) = curd(n,k)*mac_pot(n,1);
      curqg(n,k) = curq(n,k)*mac_pot(n,1);
      E_Isat = mac_pot(n,3)*eqprime(n,1)^2 ...
            + mac_pot(n,4)*eqprime(n,1) + mac_pot(n,5);
      E_Isat = max(eqprime(n,1),E_Isat); % select higher
                                       % voltage
      if n_sub~=0 %check for any sub synchronous machine model
          n_mac_number = find(mac_sub_idx == n);
          if length(n_mac_number ~=0)
            psiqpp = mac_pot(n,14)*edprime(n,k) + ...
                     mac_pot(n,15)*psikq(n,k);
            psidpp = mac_pot(n,9)*eqprime(n,k) + ...
                     mac_pot(n,10)*psikd(n,k);
            fldcur(n,k) = E_Isat + mac_pot(n,6)...
                          *(eqprime(n,k)-psikd(n,k)) + mac_pot(n,7)...
                        *curdg(n,k);
            ed(n,k) = -mac_con(n,5)*curdg(n,k) - (psiqpp...
                      -mac_con(n,13)*curqg(n,k));
            eq(n,k) = -mac_con(n,5)*curqg(n,k) + (psidpp...
                      -mac_con(n,8)*curdg(n,k));
          end
      end
      if n_tra~=0 %check for any transient generator models
          n_mac_number = find(mac_tra_idx == n);
          if length(n_mac_number)~=0
             fldcur(n,k) = E_Isat + (mac_con(n,6) - mac_con(n,7))*curdg(n,k);
             ed(n,k) = edprime(n,k) + mac_con(n,7)*curqg(n,k);
             eq(n,k) = eqprime(n,k) - mac_con(n,7)*curdg(n,k);
          end
      end
      eterm(n,k) = sqrt(ed(n,k)^2+eq(n,k)^2);
      pelect(n,k) = eq(n,k)*curq(n,k) + ed(n,k)*curd(n,k);
      qelect(n,k) = eq(n,k)*curd(n,k) - ed(n,k)*curq(n,k);
      iterm = (pelect(n,1)-jay*qelect(n,1))/...
              (eterm(n,1)*exp(-jay*theta(n_bus,1)))*mac_pot(n,1);
      vep = eterm(n,1)*exp(jay*theta(n_bus,1))*(exc_pot(i,1) + jay*exc_pot(i,2));
      ve = vep+jay*(exc_con(i,15)...
           + (exc_pot(i,1)+jay*exc_pot(i,2))*exc_con(i,16))*iterm;
      V_E = abs(ve);
      if V_E < 1e-6
        disp('excitation system error:no supply voltage')
      I_N = 2;
      else
        I_N = exc_con(i,17)*fldcur(n,k)/V_E;
      end
      if I_N < 0.433
        F_EX = 1-0.5771*I_N;
      elseif I_N < .75
        F_EX = (.75-I_N^2)^.5;
      else
        F_EX = 1.732*(1-I_N);
        F_EX=max(F_EX,0);
      end
      V_B(i,k) = V_E*F_EX;
      % set V_R limit
      if V_R(i,k) > exc_con(i,8)
        V_R(i,k) = exc_con(i,8);
      elseif V_R(i,k) < exc_con(i,9)
        V_R(i,k) = exc_con(i,9);
      end
      Efd(i,k) = min(V_R(i,k)*V_B(i,k),exc_con(i,18));
      vex(n,k) = Efd(i,k); %set field voltage for machines
   
    else
    %vector calculation
      n = mac_int(exc_con(st3_idx,2)); % machine number vector
      nst3_tra = zeros(n_st3,1);
      nst3_sub = zeros(n_st3,1);
      for j = 1:n_st3
         if ~isempty(mac_tra_idx)
         	test = find(mac_tra_idx==n(j));
            if ~isempty(test); nst3_tra(j) = test; end
         end
         if ~isempty(mac_sub_idx)
         	test = find(mac_sub_idx==n(j));
            if ~isempty(test); nst3_sub(j) = test; end
         end
      end
      nst3_tra = find(nst3_tra~=0);
      nst3_sub = find(nst3_sub~=0);
      n_bus = bus_int(mac_con(n,2));
      V_E = zeros(n_st3,1);iterm = V_E;F_EX = V_E;ve = V_E;
      I_N = V_E; E_Isat=V_E;vep=V_E;
      curd(n,k) = sin(mac_ang(n,k)).*cur_re(n,k) - ...
                  cos(mac_ang(n,k)).*cur_im(n,k); % d-axis current
      curq(n,k) = cos(mac_ang(n,k)).*cur_re(n,k) + ...
                  sin(mac_ang(n,k)).*cur_im(n,k); % q-axis current
      curdg(n,k) = curd(n,k).*mac_pot(n,1);
      curqg(n,k) = curq(n,k).*mac_pot(n,1);
      E_Isat(n) = mac_pot(n,3).*eqprime(n,k).^2 ...
               + mac_pot(n,4).*eqprime(n,k) + mac_pot(n,5);
      E_Isat(n) = max(eqprime(n,k),E_Isat(n)); % select higher
                                         % voltage
      if ~isempty(nst3_sub) %check for any sub synchronous machine model
        n_mac_number = n(nst3_sub); 
          psiqpp = mac_pot(n_mac_number,14)...
                   *edprime(n_mac_number,k) + ...
                   mac_pot(n_mac_number,15)...
                   *psikq(n_mac_number,k);
          psidpp = mac_pot(n_mac_number,9)...
                   *eqprime(n_mac_number,k) + ...
                   mac_pot(n_mac_number,10)...
                   *psikd(n_mac_number,k);
          fldcur(n_mac_number,k) = E_Isat(n_mac_number)...
                                   + mac_pot(n_mac_number,6)...
                                   *(eqprime(n_mac_number,k)...
                                   -psikd(n_mac_number,k))...
                                   + mac_pot(n_mac_number,7)...
                                   *curdg(n_mac_number,k);
          ed(n_mac_number,k) = -mac_con(n_mac_number,5)...
                               *curdg(n_mac_number,k)...
                               - (psiqpp - mac_con(n_mac_number,13)...
                               *curqg(n_mac_number,k));
          eq(n_mac_number,k) = -mac_con(n_mac_number,5)...
                               *curqg(n_mac_number,k) +...
                               (psidpp - mac_con(n_mac_number,8)...
                               *curdg(n_mac_number,k));
      end
      if ~isempty(nst3_tra) %check for any transient generator models
        n_mac_number = n(nst3_tra);
        fldcur(n_mac_number,k) = E_Isat(n_mac_number) +...
                                 (mac_con(n_mac_number,6)...
                                 - mac_con(n_mac_number,7))...
                                 *curdg(n_mac_number,k);
        ed(n_mac_number,k) = edprime(n_mac_number,k)...
                             + mac_con(n_mac_number,7)...
                             *curqg(n_mac_number,k);
        eq(n_mac_number,k) = eqprime(n_mac_number,k)...
                             - mac_con(n_mac_number,7)...
                             *curdg(n_mac_number,k);
      end
      eterm(n,k) = sqrt(ed(n,k).^2+eq(n,k).^2);
      pelect(n,k) = eq(n,k).*curq(n,k) + ed(n,k).*curd(n,k);
      qelect(n,k) = eq(n,k).*curd(n,k) - ed(n,k).*curq(n,k);
      iterm =(pelect(n,k)-jay*qelect(n,k))./...
             (eterm(n,k).*exp(-jay*theta(n_bus,k))).*mac_pot(n,1);
      vep = eterm(n,k).*exp(jay*theta(n_bus,k)).*(exc_pot(st3_idx,1) + jay*exc_pot(st3_idx,2));
      ve = vep+jay*(exc_con(st3_idx,15)...
           + (exc_pot(st3_idx,1)+jay*exc_pot(st3_idx,2)).*exc_con(st3_idx,16)).*iterm;
      V_E = abs(ve); % this is the equivalent terminal voltage at the rectifier
                     %terminals
      ve_low=find(V_E < 1e-6);
      ve_norm= find(V_E > 1e-6);
      if ~isempty(ve_norm)
        n_ven = length(ve_norm);
        n_norm = n(ve_norm);
        I_N(ve_norm) = exc_con(st3_idx(ve_norm),17).*fldcur(n_norm,k)./V_E(ve_norm);
      end
      if ~isempty(ve_low)
        disp('excitation system error:no supply voltage at')
        n_excerr = n(ve_low)
        I_N(ve_low) = 2*ones(length(ve_low),1);
        error
      end
      
      % select operating point on the inverter
      low_IN=find(I_N < 0.433);
      if ~isempty(low_IN)
        F_EX(low_IN) = ones(length(low_IN),1)-0.5771*I_N(low_IN);
      end
      big_IN=find(I_N > .75);
      if ~isempty(big_IN)
        bigl=length(big_IN);
        F_EX(big_IN) = 1.732*(ones(bigl,1)-I_N(big_IN));
        F_EX(big_IN)=max(F_EX(big_IN),zeros(bigl,1));
      end
      mid_IN =find((I_N > 0.433) & (I_N < .75));
      if ~isempty(mid_IN)
        F_EX(mid_IN) = (0.75*ones(length(mid_IN),1)-(I_N(mid_IN)).^2).^.5;
      end
      V_B(st3_idx,k) = V_E.*F_EX;
      Efd(st3_idx,k) = min(V_R(st3_idx,k).*V_B(st3_idx,k),exc_con(st3_idx,18));
      vex(n,k) = Efd(st3_idx,k); %set field voltage for machines
    end
  %end interface
  end
  if flag == 2 % exciter dynamics calculation
    if i~=0
      n = mac_int(exc_con(i,2)); % machine number
      if exc_con(i,3) == 0  % transducer time constant = 0
        dV_TR(i,k) = 0;
        V_TR(i,k) = eterm(n,k);
      else
        dV_TR(i,k) = (-V_TR(i,k)+eterm(n,k))/exc_con(i,3);
      end
      V_I = exc_sig(i,k) + exc_pot(i,3) - V_TR(i,k);
      V_I = V_I + pss_out(i,k);
      V_I = min(exc_con(i,10),max(V_I,exc_con(i,11)));
      if exc_con(i,6) == 0 % no leadlag
        dV_As(i,k) = 0;
        V_As(i,k) = exc_con(i,12)*V_I;
        V_A(i,k) = V_As(i,k);
      else
        dV_As(i,k) = (-V_As(i,k)+exc_con(i,12)*V_I)...
                     /exc_con(i,6);
        V_A(i,k) = exc_pot(i,5)*exc_con(i,12)*V_I + ...
                   (1-exc_pot(i,5))*V_As(i,k);
      end
      dV_R(i,k) = (-V_R(i,k)+exc_con(i,4)*(V_A(i,k)...
                  -min(exc_con(i,20),exc_con(i,19)...
                *Efd(i,k))))/exc_con(i,5);
      % anti-windup reset
      if V_R > exc_con(i,8)
        V_R(i,k) = exc_con(i,8);
        if dV_R(i,k)>0.0
          dV_R(i,k)=0.0;
        end
      end
      if V_R < exc_con(i,9)
        V_R(i,k) = exc_con(i,9);
        if dV_R(i,k)<0.0
           dV_R(i,k)= 0.0;
        end
      end
      dEfd(i,k) = 0;  % zero out unnecessary
      dR_f(i,k) = 0;  % state variable derivatives
 
    else
    % vector calculation
     
      n = mac_int(exc_con(st3_idx,2)); % machine number vector
      no_TR = st3_noTR_idx;
      if ~isempty(no_TR)
        n_nTR = n(no_TR);
        dV_TR(st3_idx(no_TR),k)=zeros(length(no_TR),1);
        V_TR(st3_idx(no_TR),k)=eterm(n_nTR,k);
      end
      TR = st3_TR_idx;
      if ~isempty(TR)
        n_TR = mac_int(exc_con(st3_idx(TR),2));
        dV_TR(st3_idx(TR),k)=(eterm(n_TR,k)-V_TR(st3_idx(TR),k))...
                             ./exc_con(st3_idx(TR),3);
      end  
      V_I = exc_sig(st3_idx,k) + exc_pot(st3_idx,3) - V_TR(st3_idx,k);
      V_I = V_I + pss_out(st3_idx,k);
      V_I = min(exc_con(st3_idx,10),max(V_I,exc_con(st3_idx,11)));

 
      no_TB = st3_noTB_idx;
      if ~isempty(no_TB)
        dV_As(st3_idx(no_TB),k) = zeros(length(no_TB),1);
        V_As(st3_idx(no_TB),k) = exc_con(st3_idx(no_TB),12).*V_I(no_TB);
        V_A(st3_idx(no_TB),k) = V_As(st3_idx(no_TB),k);
      end
      TB = st3_TB_idx;
      if ~isempty(TB)
        dV_As(st3_idx(TB),k) = (-V_As(st3_idx(TB),k)+exc_con(st3_idx(TB),12)...
                               .*V_I(TB))./exc_con(st3_idx(TB),6);
        V_A(st3_idx(TB),k) = exc_pot(st3_idx(TB),5).*exc_con(st3_idx(TB),12)...
                             .*V_I(TB) + (ones(length(TB),1)-...
                             exc_pot(st3_idx(TB),5)).*V_As(st3_idx(TB),k);
      end
  
      dV_R(st3_idx,k) = (-V_R(st3_idx,k)+exc_con(st3_idx,4).*(V_A(st3_idx,k)...
                        -min(exc_con(st3_idx,20),exc_con(st3_idx,19)...
                        .*Efd(st3_idx,k))))./exc_con(st3_idx,5);
   
      % anti-windup reset
      max_lim=find(V_R(st3_idx,k) > exc_con(st3_idx,8));
      if ~isempty(max_lim)
        V_R(st3_idx(max_lim),k) = exc_con(st3_idx(max_lim),8);
        pos_rate = find(dV_R(st3_idx(max_lim),k)>0);
        n_pos = length(pos_rate);
        if n_pos~=0
           dV_R(st3_idx(max_lim(pos_rate)),k) = zeros(n_pos,1);
        end
      end
      min_lim=find(V_R(st3_idx,k) < exc_con(st3_idx,9));
      if ~isempty(min_lim)
        V_R(st3_idx(min_lim),k) = exc_con(st3_idx(min_lim),9);
        neg_rate = find(dV_R(st3_idx(min_lim),k)<0);
        n_neg = length(neg_rate);
        if n_neg~=0
           dV_R(st3_idx(min_lim(neg_rate)),k) = zeros(n_neg,1);
        end
      end
      dEfd(st3_idx,k) = zeros(n_st3,1);  % zero out unnecessary
      dR_f(st3_idx,k) = zeros(n_st3,1);  % state variable derivatives
    end
  %end rate calculation
  end
end

