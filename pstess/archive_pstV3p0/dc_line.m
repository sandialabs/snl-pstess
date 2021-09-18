function  f = dc_line(i,k,kdc,bus,flag)
%Syntax:  f = dc_line(i,kdc,bus,flag)
% 5:14 PM 15/08/97
%Purpose: Models HVDC line dynamics
% Input: i - 0 for vectorized computation only option
%        k  - integer time
%        kc - integer time for dc
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

% (c) Copyright 1991-1997 Joe H. Chow - All Rights Reserved

% History (in reverse chronological order)
%
% Version:  1.0
% Date:     February 1997
% Author:   Graham Rogers

%define global variables
global dcsp_con  dcl_con  dcc_con  dcc_pot  dc_pot
global  r_idx  i_idx n_dcl  n_conv  
global  Vdc  i_dc  
global  no_cap_idx  cap_idx  no_ind_idx  l_no_cap  l_cap

% States
%line
global i_dcr i_dci  v_dcc
global di_dcr  di_dci  dv_dcc  


f=0;
jay = sqrt(-1);
k=fix(1+(kdc+1)/10);
% check for dcline data
if ~isempty(dcsp_con)
  if flag==0
  % initialization
    % check that inductances are specified
    if ~isempty(no_ind_idx)
        error('you must specify inductances for dc lines and smoothing reactors')
    end

    if i==0
    % vector computation
      i_dcr(:,1) = i_dc(r_idx,1);
      i_dci(:,1) = i_dc(i_idx,1);
      v_dcc(:,1) = Vdc(r_idx,1) - i_dcr(:,1).*dcl_con(:,3)/2.0;
      if l_no_cap~=0
        dc_pot(no_cap_idx,3) = zeros(l_no_cap,1);
        dc_pot(no_cap_idx,1) = 1000*ones(l_no_cap,1)./(dcl_con(no_cap_idx,4) ...
                               + dcl_con(no_cap_idx,6) + dcl_con(no_cap_idx,7));
        dc_pot(no_cap_idx,2) = dcl_con(no_cap_idx,3);
        dc_pot(no_cap_idx,4) = dc_pot(no_cap_idx,1);
        dc_pot(no_cap_idx,5) = dc_pot(no_cap_idx,2);
      end
      if l_cap ~= 0
        dc_pot(cap_idx,1) = 1000*ones(l_cap,1)./(dcl_con(cap_idx,4)*.5 + dcl_con(cap_idx,6));
        dc_pot(cap_idx,2) = dcl_con(cap_idx,3)/2;
        dc_pot(cap_idx,3) = 1e6*ones(l_cap,1)./dcl_con(cap_idx,5);
        dc_pot(cap_idx,4) = 1000*ones(l_cap,1)./(dcl_con(cap_idx,4)*.5 + dcl_con(cap_idx,7));
        dc_pot(cap_idx,5) = dcl_con(cap_idx,3)/2;
      end
    else
      error(' no non-vector computation in HVDC')
    end
  end
  if flag == 1
  % network inter face - no calculation required
  end
  if flag == 2
  % rate of change of states
    if i== 0
    %vector compuation
      if l_cap~=0 
        di_dcr(cap_idx,kdc) = - dc_pot(cap_idx,1)...
                           .*(dc_pot(cap_idx,2).*i_dcr(cap_idx,kdc) + v_dcc(cap_idx,kdc) - Vdc(r_idx(cap_idx),kdc));
        di_dci(cap_idx,kdc) = -dc_pot(cap_idx,4)...
                           .*(dc_pot(cap_idx,5).*i_dci(cap_idx,kdc) - v_dcc(cap_idx,kdc) + Vdc(i_idx(cap_idx),kdc));
        dv_dcc(cap_idx,kdc) = dc_pot(cap_idx,3).*(i_dcr(cap_idx,kdc) - i_dci(cap_idx,kdc));
      end
      if l_no_cap~=0
        di_dcr(no_cap_idx,kdc) = (-dc_pot(no_cap_idx,2).*i_dcr(no_cap_idx,kdc) + ...
                               Vdc(r_idx(no_cap_idx),kdc)-Vdc(i_idx(no_cap_idx),kdc))... 
                               .*dc_pot(no_cap_idx,1);
        di_dci(no_cap_idx,kdc) = di_dcr(no_cap_idx,kdc);
        dv_dcc(no_cap_idx) = 0;
      end
    else
      error(' no non-vector calculation in HVDC')
    end
  end
end
