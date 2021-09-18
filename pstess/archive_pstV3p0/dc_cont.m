function f = dc_cont(i,k,kdc,bus,flag)
% Syntax  f = dc_cont(i,k,kdc,bus,flag)
% 5:10 PM 15/08/97
% Purpose - models hvdc pole controls
% Input: i - 0 vector computaion only for HVDC controls
%        k   - integer time
%        kdc - integer time for dc
%        bus - solved loadflow bus data
%        flag - 0 - initialization
%               1 - network interface computation
%               2 - hvdc dynamics computation 
%               3 - state matrix building 
%
% Output: f - dummy variable 
%
% Calls:
%
% Called By:

% (c) Copyright 1991 Joe H. Chow - All Rights Reserved

% History (in reverse chronological order)
%
% Version:  1.0
% Date:     April 1997
% Author:   Graham Rogers

%define global variables
global  basmva  bus_v
global  dcsp_con  dcl_con dcc_con 
global  r_idx  i_idx n_dcl  n_conv  ac_bus rec_ac_bus  inv_ac_bus
global  Vdc  i_dc  dc_pot  alpha  gamma Vdc_ref
global  cur_ord dc_sig dcc_pot dcr_dsig dci_dsig
% States
%line
global i_dcr i_dci  v_dcc
global di_dcr  di_dci  dv_dcc  
%rectifier
global v_conr  
global dv_conr  
%inverter
global v_coni
global dv_coni 



jay = sqrt(-1);
f = 0;
% check that dc controls are defined
if ~isempty(dcc_con)
   if flag == 0
      % initialize controls
      if i == 0 
         % vector computation
         % rectifier controls
         % calculate current order
         cur_ord(r_idx,1) = i_dc(r_idx,1);
         cur_ord(i_idx,1) = i_dc(i_idx,1).*(ones(n_dcl,1)-dcl_con(:,9)/100.0);
         Vdc_ref = Vdc(i_idx,1);
         v_conr(:,1) = alpha(:,1)./dcc_con(r_idx,4);
         vrmax=find(v_conr(:,1)>dcc_con(r_idx,5));
         vrmin=find(v_conr(:,1)<dcc_con(r_idx,6));
         % check rectifier integrator limit violation
         if ~isempty(vrmax)
            svrm=int2str(vrmax);
            disp('v_conr greater than maximum limit at the following rectifiers')
            disp(svrm)
            error
         end
         if ~isempty(vrmin)
            svrm=int2str(vrmin);
            disp('v_conr less than minimum limit at the following rectifiers')
            disp(svrm)
            error
         end
         
         v_coni(:,1) = gamma(:,1)./dcc_con(i_idx,4);
         vimax=find(v_coni(:,1)>dcc_con(i_idx,5));
         vimin=find(v_coni(:,1)<dcc_con(i_idx,6));
         % check rectifier integrator limit violation
         if ~isempty(vimax)
            svim=int2str(vimax);
            disp('v_coni greater than maximum limit at the following inverters')
            disp(svim)
            error
         end
         if ~isempty(vimin)
            svim=int2str(vimin);
            disp('v_coni less than minimum limit at the following inverters')
            disp(svim)
            error
         end
         
         dcc_pot(:,1) = gamma(:,1); %gamma reference value
         dcc_pot(:,2) = dcsp_con(r_idx,5)./dcsp_con(r_idx,6); 
         dcc_pot(:,2) = dcc_pot(:,2)*basmva./bus(rec_ac_bus,13)./bus(rec_ac_bus,13); % xeqpu rec
         dcc_pot(:,3) = dcsp_con(r_idx,6).*dcsp_con(r_idx,5)*3/pi;  % Rc rectifiers   
         dcc_pot(:,4) = dcsp_con(i_idx,5)./dcsp_con(i_idx,6); 
         dcc_pot(:,4) = dcc_pot(:,4)*basmva./bus(inv_ac_bus,13)./bus(inv_ac_bus,13); %xeqpu inv
         dcc_pot(:,5) = dcsp_con(i_idx,6).*dcsp_con(i_idx,5)*3/pi;  % Rc inverters  
         dcc_pot(:,6) = dcc_pot(:,5); 
         % multiplier ideal rec dc voltage 
         dcc_pot(:,7) = 3*sqrt(2).*dcsp_con(r_idx,6).*bus(rec_ac_bus,13)/pi;
         % multiplier ideal inv dc voltage
         dcc_pot(:,8) = 3*sqrt(2).*dcsp_con(i_idx,6).*bus(inv_ac_bus,13)/pi;
         % multiplier acpu to dc amps
         dcc_pot(:,9) = pi*basmva/3/sqrt(2)./bus(rec_ac_bus,13)./dcsp_con(r_idx,6);%rectifier
         dcc_pot(:,10) = pi*basmva/3/sqrt(2)./bus(inv_ac_bus,13)./dcsp_con(i_idx,6);%inverter
         dc_dsig(:,1)=zeros(n_conv,1); % zero damping control signals
         if dc_sig(:,1)~=zeros(n_conv,1)
            % reset initial values of alpha and gamma
            alpha(:,1) = ((-cur_ord(r_idx,1) + dc_sig(r_idx,1) + dcr_dsig(:,1)...
               + i_dcr(:,1)).*dcc_con(r_idx,2) + ...    
               v_conr(:,1)).*dcc_con(r_idx,4); 
            gamma(:,1) = ((Vdc(i_idx,1)-Vdc_ref)./Vdc_ref.*dcc_con(i_idx,2)...
               + v_coni(:,1)).*dcc_con(i_idx,4);
         end
         
      else
         error('vector computation only in dc controls')
      end
      % end of initialization
   end
   if flag== 1
      % network interface
      if i ~= 0
         error('vector computation only with dc')
      else
         % vector computation 
         % i_dc vcdc and the control states are fixed
         
         % determine firing and extinction angles 
         alpha(:,kdc) = ((-cur_ord(r_idx,k) + dc_sig(r_idx,k) + dcr_dsig(:,k)...
            + i_dcr(:,kdc)).*dcc_con(r_idx,2) + ...    
            v_conr(:,kdc)).*dcc_con(r_idx,4); 
         % check for alpha limits
         alpha(:,kdc) = max(alpha(:,kdc), dcc_con(r_idx,8)*pi/180);
         alpha(:,kdc) = min(alpha(:,kdc), dcc_con(r_idx,7)*pi/180);
         
         gamma(:,kdc) = ((Vdc(i_idx,kdc)-Vdc_ref)./Vdc_ref.*dcc_con(i_idx,2)...
            + v_coni(:,kdc)).*dcc_con(i_idx,4);
         cur_error = i_dci(:,kdc) - cur_ord(i_idx,k);
         ce_idx = find(cur_error < 0);
         if ~isempty(ce_idx)
            gamma(ce_idx,kdc) = gamma(ce_idx,kdc) + cur_error(ce_idx).*...
               dcc_con(i_idx(ce_idx),2).*dcc_con(i_idx(ce_idx),4);
         end
         % check gamma limits
         gamma(:,kdc) = max(gamma(:,kdc), dcc_con(i_idx,8)*pi/180);
         gamma(:,kdc) = min( gamma(:,kdc), dcc_con(i_idx,7)*pi/180);
      end 
   end
   if flag == 2
      % calculate rates of change of states
      if i==0
         % vector computation 
         % rectifier 
         dv_conr(:,kdc) = (-cur_ord(r_idx,k) + i_dcr(:,kdc) + dc_sig(r_idx,k) + dcr_dsig(:,k))...
            .*dcc_con(r_idx,3);
         %check for state limits
         recmx  = find(v_conr(:,kdc)>dcc_con(r_idx,5));
         if ~isempty(recmx)
            v_conr(recmx,kdc) = dcc_con(r_idx(recmx),5);
            recdmx = find(dv_conr(recmx,k)>0);
            if ~isempty(recdmx)
               dv_conr(recmx(recdmx),kdc) = zeros(length(recdmx),1);
            end
         end
         recmn  = find(v_conr(:,kdc)<dcc_con(r_idx,6));
         if ~isempty(recmn)
            v_conr(recmn,kdc) = dcc_con(r_idx(recmn),6);
            recdmn = find(dv_conr(recmn,kdc)<0);
            if ~isempty(recdmn)
               dv_conr(recmn(recdmn),kdc) = zeros(length(recdmn),1);
            end
         end  
         
         %inverter
         cur_err = cur_ord(i_idx,k) - i_dci(:,kdc);
         n_ce = find(cur_err<0);
         if ~isempty(n_ce)
            cur_err(n_ce) = zeros(length(n_ce),1);
         end
         inv_err = (Vdc(i_idx,kdc) - Vdc_ref)./Vdc_ref - cur_err;
         dv_coni(:,kdc) = (inv_err + dc_sig(i_idx,k)+dci_dsig(:,k))...
            .*dcc_con(i_idx,3);
         % check state limits
         invmx  = find(v_coni(:,kdc)>dcc_con(i_idx,5));
         if ~isempty(invmx)
            v_coni(invmx,kdc) = dcc_con(i_idx(invmx),5);
            invdmx = find(dv_coni(invmx,k)>0);
            if ~isempty(invdmx)
               dv_coni(recmx(invdmx),kdc) = zeros(length(invdmx),1);
            end
         end
         invmn  = find(v_coni(:,kdc)<dcc_con(i_idx,6));
         if ~isempty(invmn)
            v_coni(invmn,kdc) = dcc_con(i_idx(invmn),6);
            invdmn = find(dv_coni(invmn,kdc)<0);
            if ~isempty(invdmn)
               dv_conr(invmn(invdmn),kdc) = zeros(length(invdmn),1);
            end
         end  
      end
   end
end
















