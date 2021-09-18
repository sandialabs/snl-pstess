% s_simu.m  
% 9:59 AM 14 June 1999
% An m.file to simulate power system transients  
% using the Matlab Power System Toolbox
% This m-file takes the dynamic and load flow data and
% calculates the response of the power system to a fault
% which is specified in a switching file
% see one of the supplied data files (data2a.m) for the 
% switching file format

% version 1.8
% Author: Rui
% July 19, 2017
% add code for multiple HVDC systems 

% version 1.7
% Author: Graham Rogers
% Date September 1999
% Add smple exciter with pi avr - smppi
% version 1.6
% Author: Graham Rogers
% Date: June 1999
% Modification: add user defined hvdc control
%               modify dc so that it is integrated with a sub-multiple time constant

% version 1.5
% Author: Graham Rogers
% Date: December 1998/ January 1999
% Modification: Add svc and tcsc user defined damping controls, and add tcsc model

% version 1.4
% Author: Graham Rogers
% Date:  July 1998
% Modification: Add deltaP/omega filter

% version 1.3
% Author: Graham Rogers
% Date:  June 1998
% Modification: Add hydraulic turbine/governor

% Version: 1.2
% Author:  Graham Rogers
% Date:    August 1997
% Modification : add induction generator
% Version: 1.1
% Author:  Graham Rogers
% Date:    August 1997
% Modification : add modulation for load, exciter and governor
% Version  1.0
% Author:  Graham Rogers
% Date: February 1997
%
% (c) Copyright: Joe Chow/ Cherry Tree Scientific Software 1991 to 2020 - All rights reserved
% 
% Permission is hereby granted, free of charge, to any person obtaining a copy
% of this software and associated documentation files (the "Software"), to deal
% in the Software without restriction, including without limitation the rights
% to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
% copies of the Software, and to permit persons to whom the Software is
% furnished to do so.
% 
% The above copyright notice and this permission notice shall be included in all
% copies or substantial portions of the Software.
% 
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
% IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
% FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
% AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
% LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
% OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
% SOFTWARE.
%
clear all
clear global
close % close graphics windows
tic % set timer
plot_now=0;
jay = sqrt(-1);
pst_var % set up global variables 
svc_dc=[];tcsc_dc=[];dcr_dc=[];dci_dc=[];
% load input data from m.file
disp('non-linear simulation')
% input data file
[dfile,pathname]=uigetfile('d*.m','Select Data File');
if pathname == 0
   error(' you must select a valid data file')
else
   lfile =length(dfile);
   % strip off .m and convert to lower case
   dfile = lower(dfile(1:lfile-2));
   eval(dfile);
end
% check for valid dynamic data file
if isempty(mac_con)
   error(' the selected file is not a valid data file')
end
if isempty(sw_con)
   error(' the selected file has no switching data')
end
basdat = inputdlg({'Base MVA:','Base Frequency Hz:'},'Input Base Data',1,{'100','60'}); 
sys_freq = str2double(basdat{2});
basrad = 2*pi*sys_freq; % default system frequency is 60 Hz
basmva = str2double(basdat{1});
syn_ref = 0 ;     % synchronous reference frame
ibus_con = []; % ignore infinite buses in transient simulation

% solve for loadflow - loadflow parameter
lfs = inputdlg('Do you want to solve loadflow > (y/n)[y] ','s');
if isempty(lfs{1});lfs{1}='y';end
if lfs{1}=='Y';lfs{1}='y';end
if lfs{1} == 'y'
   if isempty(dcsp_con)
      n_conv = 0;
      n_dcl = 0;
      ndcr_ud=0;
      ndci_ud=0;
      tol = 1e-8;   % tolerance for convergence
      iter_max = 30; % maximum number of iterations
      acc = 1.0;   % acceleration factor
      [bus_sol,line,line_flw] = ...
         loadflow(bus,line,tol,iter_max,acc,'n',2);
      bus = bus_sol;  % solved loadflow solution needed for
      % initialization
      save sim_fle bus line 
   else
      [bus_sol,line,line_flw,rec_par,inv_par, line_par] = lfdcs(bus,line,dci_dc,dcr_dc);
      bus = bus_sol;
      save sim_fle bus line rec_par  inv_par line_par
   end 
else
   load sim_fle 
end
%set indexes
% note: dc index set in dc load flow
mac_indx;
exc_indx;
tg_indx;
dpwf_indx;
pss_indx;
svc_indx(svc_dc);
tcsc_indx(tcsc_dc);
lm_indx;
rlm_indx;
n_mot = size(ind_con,1);
n_ig = size(igen_con,1);
if isempty(n_mot); n_mot = 0;end
if isempty(n_ig); n_ig = 0; end
ntot = n_mac+n_mot+n_ig;
ngm = n_mac + n_mot;
n_pm = n_mac;
disp(' ')
disp('Performing simulation.')
%
% construct simulation switching sequence as defined in sw_con
tswitch(1) = sw_con(1,1);
k = 1;kdc=1;
n_switch = length(sw_con(:,1));
k_inc = zeros(n_switch-1,1);k_incdc=k_inc;
t_switch = zeros(n_switch,1);
h=t_switch;h_dc=h;
for sw_count = 1:n_switch-1
   h(sw_count) = sw_con(sw_count,7);%specified time step
   if h(sw_count)==0, h(sw_count) = 0.01;end % default time step
   k_inc(sw_count) = fix((sw_con(sw_count+1,1)-sw_con(sw_count,1))/h(sw_count));%nearest lower integer
   if k_inc(sw_count)==0;k_inc(sw_count)=1;end% minimum 1
   h(sw_count) = (sw_con(sw_count+1,1)-sw_con(sw_count,1))/k_inc(sw_count);%step length
   h_dc(sw_count) = h(sw_count)/10;
   k_incdc(sw_count) = 10*k_inc(sw_count);
   t_switch(sw_count+1) =t_switch(sw_count) +  k_inc(sw_count)*h(sw_count);
   t(k:k-1+k_inc(sw_count)) = t_switch(sw_count):h(sw_count):t_switch(sw_count+1)-h(sw_count);
   t_dc(kdc:kdc-1+k_incdc(sw_count)) = t_switch(sw_count):h_dc(sw_count):t_switch(sw_count+1)-h_dc(sw_count);
   k=k+k_inc(sw_count);kdc=kdc+k_incdc(sw_count);
end
t_dc(kdc)=t_dc(kdc-1)+h_dc(sw_count);
for kk=1:10;kdc=kdc+1;t_dc(kdc)=t_dc(kdc-1)+h_dc(sw_count);end

k = sum(k_inc)+1; % k is the total number of time steps in the simulation

t(k) = sw_con(n_switch,1);
[n dummy]=size(mac_con) ;
n_bus = length(bus(:,1));
%
% create zero matrices for variables to make algorithm more efficient?
z = zeros(n,k);
z1 = zeros(1,k);
zm = zeros(1,k);if n_mot>1;zm = zeros(n_mot,k);end
zig = zeros(1,k);if n_ig>1;zig = zeros(n_ig,k);end
zdc = zeros(2,kdc);if n_conv>2; zdc = zeros(n_conv,kdc);end
zdcl = zeros(1,kdc);if n_dcl>1;zdcl=zeros(n_dcl,kdc);end
% set dc parameters  
Vdc = zeros(n_conv,kdc);
i_dc = zdc;  
P_dc = z; cur_ord = z;
alpha = zdcl; 
gamma = zdcl;  
dc_sig = zeros(n_conv,k);dcr_dsig = zeros(n_dcl,k);dci_dsig=zeros(n_dcl,k);
i_dcr = zdcl; i_dci = zdcl; v_dcc = zdcl;
di_dcr = zdcl; di_dci = zdcl; dv_dcc = zdcl;
v_conr = zdcl; v_coni = zdcl; dv_conr = zdcl; dv_coni = zdcl;
if n_conv~=0
    % Modified by Rui on Oct. 5, 2016 
    for ihvdc_count=1:n_dcl
        ire=r_idx(ihvdc_count);
        Vdc(ire,:)=rec_par(ihvdc_count,2);
        i_dc(ire,:)=line_par(ihvdc_count);    % for PDCI
        i_dcr(ihvdc_count,:) = i_dc(ire,:);
        alpha(ihvdc_count,:) = rec_par(ihvdc_count,1)*pi/180;
    end
            
    for ihvdc_count=1:n_dcl
        iin=i_idx(ihvdc_count);
        Vdc(iin,:)=inv_par(ihvdc_count,2);
        i_dc(iin,:) = line_par(ihvdc_count);
        i_dci(ihvdc_count,:) = i_dc(iin,:);
        gamma(ihvdc_count,:) = inv_par(ihvdc_count,1)*pi/180;
    end
    % end modification by Rui

   Vdc(r_idx,:) = rec_par(:,2); Vdc(i_idx,:) = inv_par(:,2);
   i_dc(r_idx,:) = line_par; i_dc(i_idx,:) = line_par;
   i_dcr(:,:) = i_dc(r_idx,:); i_dci(:,:) = i_dc(i_idx,:);
   alpha(:,:) = rec_par(:,1)*pi/180;
   gamma(:,:) = inv_par(:,1)*pi/180;
   if ndcr_ud~=0
      for j = 1:ndcr_ud
         sv = get(dcr_dc{j,1});
         if j==1
            xdcr_dc =zeros(sv.NumStates,kdc); 
            dxdcr_dc = zeros(sv.NumStates,kdc);
         else
            xdcr_dc = [xdcr_dc;zeros(sv.NumStates,kdc)];
            dxdcr_dc = [dxdcr_dc;zeros(sv.NumStates,kdc)];
         end
      end
      dcrd_sig=zeros(ndcr_ud,k);
      angdcr = zeros(ndcr_ud,k);
   else
      xdcr_dc = zeros(1,kdc);
      dxdcr_dc = zeros(1,kdc);
      dcrd_sig = zeros(1,k);
   end
   if ndci_ud~=0
      for j = 1:ndci_ud
         sv = get(dci_dc{j,1});
         if j==1
            xdci_dc =zeros(sv.NumStates,kdc); 
            dxdci_dc = zeros(sv.NumStates,kdc);
         else
            xdci_dc = [xsvc_dc;zeros(sv.NumStates,kdc)];
            dxdci_dc = [dxsvc_dc;zeros(sv.NumStates,kdc)];
         end
      end
      dcid_sig=zeros(ndcr_ud,k);
      angdci = zeros(ndci_dc,k);
   else
      xdci_dc = zeros(1,kdc);
      dxdci_dc = zeros(1,kdc);
      dcid_sig = zeros(1,k);
   end   
else
   xdcr_dc = zeros(1,kdc);
   dxdcr_dc = zeros(1,kdc);
   xdci_dc = zeros(1,kdc);
   dxdci_dc = zeros(1,kdc);
end
v_p = z1;
theta = zeros(n_bus+1,k);bus_v = zeros(n_bus+1,k);
mac_ang = z; mac_spd = z; dmac_ang = z; dmac_spd = z;
pmech = z; pelect = z; mac_ref = z1;  sys_ref = z1; 
edprime = z; eqprime = z; dedprime = z; deqprime = z;
psikd = z; psikq = z; dpsikd = z; dpsikq = z;
pm_sig = z;
z_tg = zeros(1,k);if n_tg+n_tgh~=0;z_tg = zeros(n_tg+n_tgh,k);end
tg1 = z_tg; tg2 = z_tg; tg3 = z_tg; tg4 = z_tg; tg5 = z_tg;
dtg1 = z_tg; dtg2 = z_tg; dtg3 = z_tg;dtg4 = z_tg; dtg5 = z_tg;
tg_sig = z_tg;
z_pss = zeros(1,k);if n_pss~=0;z_pss = zeros(n_pss,k);end
pss1 = z_pss; pss2 = z_pss; pss3 = z_pss;
dpss1 = z_pss; dpss2 = z_pss; dpss3 = z_pss;
z_dpw = zeros(1,k);if n_dpw~=0; z_dpw = zeros(n_dpw,k);end
sdpw1 = z_dpw; sdpw2 = z_dpw; sdpw3 = z_dpw; sdpw4 = z_dpw; sdpw5 = z_dpw; sdpw6 = z_dpw; dpw_out = z_dpw;
dsdpw1 = z_dpw; dsdpw2 = z_dpw; dsdpw3 = z_dpw; dsdpw4 = z_dpw; dsdpw5 = z_dpw; dsdpw6 = z_dpw;
curd = z; curq = z; curdg = z; curqg = z; fldcur = z;
ed = z; eq = z; eterm = z; qelect = z;
vex = z; cur_re = z; cur_im = z; psi_re = z; psi_im = z;
ze = zeros(1,k);if n_exc~=0; ze = zeros(n_exc,k);end
V_B = ze;exc_sig = ze;
V_TR = ze; V_R = ze; V_A = ze; V_As = ze; Efd = ze; R_f = ze;
dV_TR = ze; dV_R = ze; dV_As = ze; dEfd = ze; dR_f = ze;
pss_out = ze;
vdp = zm; vqp = zm; slip = zm; 
dvdp = zm; dvqp = zm; dslip = zm;
s_mot = zm; p_mot = zm; q_mot = zm;
vdpig = zig; vqpig = zig; slig = zig; 
dvdpig = zig; dvqpig = zig; dslig = zig;
s_igen = zig; pig = zig; qig = zig; tmig = zig;
if n_svc~=0
   B_cv = zeros(n_svc,k); dB_cv = zeros(n_svc,k);svc_sig = zeros(n_svc,k);svc_dsig=zeros(n_svc,k);
   B_con = zeros(n_svc,k); dB_con=zeros(n_svc,k);
   if n_dcud~=0
      d_sig = zeros(n_dcud,k);
      for j = 1:n_dcud
         sv = get(svc_dc{j,1});
         if j==1
            xsvc_dc =zeros(sv.NumStates,k); 
            dxsvc_dc = zeros(sv.NumStates,k);
         else
            xsvc_dc = [xsvc_dc;zeros(sv.NumStates,k)];
            dxsvc_dc = [dxsvc_dc;zeros(sv.NumStates,k)];
         end
      end
   else
      xsvc_dc = zeros(1,k);
      dxsvc_dc = zeros(1,k);
   end
else
   B_cv = zeros(1,k);dB_cv = zeros(1,k); svc_sig = zeros(1,k);svc_dsig = zeros(1,k);
   B_con = zeros(1,k);dB_con=zeros(1,k);
   xsvc_dc = zeros(1,k);dxsvc_dc = zeros(1,k);
   d_sig = zeros(1,k);
end
if n_tcsc~=0
   B_tcsc = zeros(n_tcsc,k); dB_tcsc = zeros(n_tcsc,k);tcsc_sig = zeros(n_tcsc,k);tcsc_dsig=zeros(n_tcsc,k);
   if n_tcscud~=0
      td_sig = zeros(n_tcscud,k);%input to tcsc damping control
      for j = 1:n_tcscud
         sv = get(tcsc_dc{j,1});% damping control state space object
         if j==1
            xtcsc_dc =zeros(sv.NumStates,k); % tcsc damping control states
            dxtcsc_dc = zeros(sv.NumStates,k);% tcsc dc rates of chage of states
         else
            xtcsc_dc = [xtcsc_dc;zeros(sv.NumStates,k)];% in order of damping controls
            dxtcsc_dc = [dxtcsc_dc;zeros(sv.NumStates,k)];
         end
      end
   else
      xtcsc_dc = zeros(1,k);
      dxtcsc_dc = zeros(1,k);
   end
else
   B_tcsc = zeros(1,k);dB_tcsc = zeros(1,k); tcsc_sig = zeros(1,k);tcsc_dsig = zeros(1,k);
   xtcsc_dc = zeros(1,k);dxtcsc_dc = zeros(1,k);
   td_sig = zeros(1,k);
end

if n_lmod ~= 0
   lmod_st = zeros(n_lmod,k); dlmod_st = lmod_st; lmod_sig = lmod_st;
else
   lmod_st = zeros(1,k); dlmod_st = lmod_st; lmod_sig = lmod_st;
end
if n_rlmod ~= 0
   rlmod_st = zeros(n_rlmod,k); drlmod_st = rlmod_st; rlmod_sig = rlmod_st;
else
   rlmod_st = zeros(1,k); drlmod_st = rlmod_st; rlmod_sig = rlmod_st;
end

sys_freq = ones(1,k);
disp('constructing reduced y matrices')
% step 1: construct reduced Y matrices 

disp('initializing motor,induction generator, svc and dc control models')       
bus = mac_ind(0,1,bus,0);% initialize induction motor
bus = mac_igen(0,1,bus,0); % initialize induction generator
bus = svc(0,1,bus,0);%initialize svc
f = dc_cont(0,1,1,bus,0);% initialize dc controls
% this has to be done before red_ybus is used since the motor and svc 
% initialization alters the bus matrix and dc parameters are required

y_switch % calculates the reduced y matrices for the different switching conditions
disp('initializing other models')
% step 2: initialization
theta(1:n_bus,1) = bus(:,3)*pi/180;
bus_v(1:n_bus,1) = bus(:,2).*exp(jay*theta(1:n_bus,1));
if n_dcud ~=0
   % calculate the initial magnitude of line current for svc damping controls
   for j=1:n_dcud
      l_num = svc_dc{j,3};svc_num = svc_dc{j,2};
      from_bus = bus_int(line(l_num,1)); to_bus=bus_int(line(l_num,2));
      svc_bn = bus_int(svc_con(svc_num,2));
      if svc_bn~=from_bus&&svc_bn~=to_bus;error('the svc is not at the end of the specified line');end
      V1 = bus_v(from_bus,1);
      V2 = bus_v(to_bus,1);
      R = line(l_num,3);X=line(l_num,4);B=line(l_num,5);tap = line(l_num,6);phi = line(l_num,7);
      [l_if,l_it] = line_cur(V1,V2,R,X,B,tap,phi);l_if0(j)=l_if;l_it0(j)=l_it;
      if svc_bn == from_bus; d_sig(j,1)=abs(l_if);elseif svc_bn==to_bus;d_sig(j,1)=abs(l_it);end
   end
end
if n_tcscud ~=0
   % calculate the initial magnitude of bus voltage magnitude for tcsc damping controls
   for j=1:n_tcscud
      b_num = tcsc_dc{j,3};tcsc_num = tcsc_dc{j,2};
      td_sig(j,1)=abs(bus_v(bus_int(b_num),1));
   end
end

if n_conv~=0
   % change dc buses from LT to HT
   Pr = bus(rec_ac_bus,6);
   Pi = bus(inv_ac_bus,6);
   Qr = bus(rec_ac_bus,7);
   Qi = bus(inv_ac_bus,7);
   VLT= bus_v(ac_bus,1);
   i_acr = (Pr-jay*Qr)./conj(VLT(r_idx));
   i_aci = (Pi - jay*Qi)./conj(VLT(i_idx));
   IHT(r_idx,1)=i_acr;
   IHT(i_idx,1)=i_aci;
   VHT(r_idx,1) = (VLT(r_idx) + jay*dcc_pot(:,2).*i_acr);
   VHT(i_idx,1) = (VLT(i_idx) + jay*dcc_pot(:,4).*i_aci);
   bus_v(ac_bus,1) = VHT;
   theta(ac_bus,1) = angle(bus_v(ac_bus,1));
   % modify the bus matrix to the HT buses
   bus(ac_bus,2) = abs(bus_v(ac_bus,1));
   bus(ac_bus,3) = theta(ac_bus,1)*180/pi;
   SHT = VHT.*conj(IHT);
   bus(ac_bus,6) = real(SHT);
   bus(ac_bus,7) = imag(SHT);
   if ndcr_ud~=0
      % calculate the initial value of bus angles rectifier user defined control
      for j = 1:ndcr_ud
         b_num1 = dcr_dc{j,3};b_num2 = dcr_dc{j,4};conv_num = dcr_dc{j,2};
         angdcr(j,:)=theta(bus_int(b_num1),1)-theta(bus_int(b_num2),1);
         dcrd_sig(j,:)=angdcr(j,:);
      end
   end
   if ndci_ud~=0
      % calculate the initial value of bus angles inverter user defined control
      for j = 1:ndci_ud
         b_num1 = dci_dc{j,3};b_num2 = dci_dc{j,4};conv_num = dci_dc{j,2};
         angdci(j,:)=theta(bus_int(b_num1),1)-theta(bus_int(b_num2),1);
         dcid_sig(j,:)=angdci(j,:);
      end
   end      
end
flag = 0;
bus_int = bus_intprf;% pre-fault system
disp('generators')
mac_sub(0,1,bus,flag);   
mac_tra(0,1,bus,flag);
mac_em(0,1,bus,flag);
disp('generator controls')
dpwf(0,1,bus,flag);
pss(0,1,bus,flag);
smpexc(0,1,bus,flag);
smppi(0,1,bus,flag);
exc_st3(0,1,bus,flag);
exc_dc12(0,1,bus,flag);
tg(0,1,bus,flag); 
tg_hydro(0,1,bus,flag);

% initialize svc damping controls
if n_dcud~=0
   tot_states=0;
   for i = 1:n_dcud
      ysvcmx = svc_dc{i,4};ysvcmn = svc_dc{i,5};
      svc_num = svc_dc{i,2};
      st_state = tot_states+1; svc_states = svc_dc{i,6}; tot_states = tot_states+svc_states; 
      [svc_dsig(svc_num,1),xsvc_dc(st_state:tot_states,1),dxsvc_dc(st_state:tot_states,1)] =...
         svc_sud(i,1,flag,svc_dc{i,1},d_sig(i,1),ysvcmx,ysvcmn);
   end
end
% initialize tcsc damping controls
if n_tcscud~=0
   tot_states=0;
   for i = 1:n_tcscud
      ytcscmx = tcsc_dc{i,4};ytcscmn = tcsc_dc{i,5};
      tcsc_num = tcsc_dc{i,2};
      st_state = tot_states+1; tcsc_states = tcsc_dc{i,6}; tot_states = tot_states+tcsc_states; 
      [tcsc_dsig(tcsc_num,1),xtcsc_dc(st_state:tot_states,1),dxtcsc_dc(st_state:tot_states,1)] =...
         tcsc_sud(i,1,flag,tcsc_dc{i,1},td_sig(i,1),ytcscmx,ytcscmn);
   end
end
if ndcr_ud~=0
   % initialize rectifier damping controls
   tot_states=0;
   for i = 1:ndcr_ud
      ydcrmx = dcr_dc{i,5};ydcrmn = dcr_dc{i,6};
      rec_num = dcr_dc{i,2};
      st_state = tot_states+1; dcr_states = dcr_dc{i,7}; tot_states = tot_states+dcr_states; 
      [dcr_dsig(rec_num,1),xdcr_dc(st_state:tot_states,1),dxdcr_dc(st_state:tot_states,1)] = ...
         dcr_sud(i,1,flag,dcr_dc{i,1},dcrd_sig(i,1),ydcrmx,ydcrmn);
   end
end
if ndci_ud~=0
   % initialize inverter damping controls
   tot_states=0;
   for i = 1:ndci_ud
      ydcimx = dci_dc{i,5};ydcrmn = dci_dc{i,6};
      inv_num = dci_dc{i,2};
      st_state = tot_states+1; dci_states = dci_dc{i,7}; tot_states = tot_states+dci_states; 
      [dci_dsig(inv_num,1),xdci_dc(st_state:tot_states,1),dxdci_dc(st_state:tot_states,1)] =...
         dci_sud(i,1,flag,dci_dc{i,1},dcid_sig(i,1),ydcimx,ydcimn);
   end
end


% initialize load modulation control
if ~isempty(lmod_con)
   disp('load modulation')
   lmod(0,1,bus,flag);
end
if ~isempty(rlmod_con)
   disp('reactive load modulation')
   rlmod(0,1,bus,flag);
end

% initialize non-linear loads
if ~isempty(load_con)
   disp('non-linear loads')
   vnc = nc_load(bus,flag,Y_ncprf,Y_ncgprf);
else
   nload = 0;
end
if ~isempty(dcsp_con)
   bus_sim = bus;
   bus_int = bus_intprf;
   Y1 = Y_gprf;
   Y2 = Y_gncprf;
   Y3 = Y_ncgprf;
   Y4 = Y_ncprf;
   Vr1 = V_rgprf; 
   Vr2 = V_rncprf;
   bo = boprf;
   h_sol = i_simu(1,1,k_inc,h,bus_sim,Y1,Y2,Y3,Y4,Vr1,Vr2,bo);
   % reinitialize dc controls
   mdc_sig(0,1);
   dc_cont(0,1,1,bus,flag);
   % initialize dc line
   dc_line(0,1,1,bus,flag);  
end
H_sum = sum(mac_con(:,16)./mac_pot(:,1));
% step 3: perform a predictor-corrector integration 
kt = 0;
ks = 1;
k_tot = sum(k_inc);
lswitch = length(k_inc);
ktmax = k_tot-k_inc(lswitch);
bus_sim = bus;
plot_now = 0;
while (kt<=ktmax)
   k_start = kt+1;
   if kt==ktmax
      k_end = kt + k_inc(ks);
   else
      k_end = kt + k_inc(ks) + 1;
   end
   for k = k_start:k_end
      % step 3a: network solution
      % mach_ref(k) = mac_ang(syn_ref,k);
      mach_ref(k) = 0;
      pmech(:,k+1) = pmech(:,k);
      tmig(:,k+1) = tmig(:,k);
      
      if n_conv~=0;cur_ord(:,k+1) = cur_ord(:,k);end
      
      flag = 1;
      timestep = int2str(k);
      % network-machine interface
      mac_ind(0,k,bus_sim,flag);
      mac_igen(0,k,bus_sim,flag);
      mac_sub(0,k,bus_sim,flag); 
      mac_tra(0,k,bus_sim,flag);
      mac_em(0,k,bus_sim,flag);
      mdc_sig(t(k),k);
      dc_cont(0,k,10*(k-1)+1,bus_sim,flag);
      % Calculate current injections and bus voltages and angles
      if k >= sum(k_inc(1:3))+1
         % fault cleared
         line_sim = line_pf2;
         bus_sim = bus_pf2;
         bus_int = bus_intpf2;
         Y1 = Y_gpf2;
         Y2 = Y_gncpf2;
         Y3 = Y_ncgpf2;
         Y4 = Y_ncpf2;
         Vr1 = V_rgpf2; 
         Vr2 = V_rncpf2;
         bo = bopf2;
         h_sol = i_simu(k,ks,k_inc,h,bus_sim,Y1,Y2,Y3,Y4,Vr1,Vr2,bo);     
      elseif k >=sum(k_inc(1:2))+1
         % near bus cleared
         line_sim = line_pf1;
         bus_sim = bus_pf1;
         bus_int = bus_intpf1;
         Y1 = Y_gpf1;
         Y2 = Y_gncpf1;
         Y3 = Y_ncgpf1;
         Y4 = Y_ncpf1;
         Vr1 = V_rgpf1; 
         Vr2 = V_rncpf1;
         bo = bopf1;
         h_sol = i_simu(k,ks,k_inc,h,bus_sim,Y1,Y2,Y3,Y4,Vr1,Vr2,bo);   
      elseif k>=k_inc(1)+1
         % fault applied
         line_sim = line_f;
         bus_sim = bus_f;
         bus_int = bus_intf;
         Y1 = Y_gf;
         Y2 = Y_gncf;
         Y3 = Y_ncgf;
         Y4 = Y_ncf;
         Vr1 = V_rgf; 
         Vr2 = V_rncf;
         bo = bof;
         h_sol = i_simu(k,ks,k_inc,h,bus_sim,Y1,Y2,Y3,Y4,Vr1,Vr2,bo);     
      elseif k<k_inc(1)+1
         % pre fault
         line_sim = line;
         bus_sim = bus;
         bus_int = bus_intprf;
         Y1 = Y_gprf;
         Y2 = Y_gncprf;
         Y3 = Y_ncgprf;
         Y4 = Y_ncprf;
         Vr1 = V_rgprf; 
         Vr2 = V_rncprf;
         bo = boprf;
         h_sol = i_simu(k,ks,k_inc,h,bus_sim,Y1,Y2,Y3,Y4,Vr1,Vr2,bo);  
      end
      % HVDC
      if ndcr_ud~=0
         % calculate the new value of bus angles rectifier user defined control
         tot_states=0;
         for jj = 1:ndcr_ud
            b_num1 = dcr_dc{jj,3};b_num2 = dcr_dc{jj,4};conv_num = dcr_dc{jj,2};
            angdcr(jj,k)=(theta(bus_int(b_num1),k)-theta(bus_int(b_num2),k));
            dcrd_sig(jj,k)=angdcr(jj,k);
            st_state = tot_states+1; dcr_states = dcr_dc{jj,7}; tot_states = tot_states+dcr_states; 
            ydcrmx=dcr_dc{jj,5};ydcrmn = dcr_dc{jj,6};
            dcr_dsig(jj,k) = ...
               dcr_sud(jj,k,flag,dcr_dc{jj,1},dcrd_sig(jj,k),ydcrmx,ydcrmn,xdcr_dc(st_state:tot_states,10*(k-1)+1));
         end
      end
      if ndci_ud~=0
         % calculate the new value of bus angles inverter user defined control
         for jj = 1:ndci_ud
            tot_states=0;
            b_num1 = dci_dc{jj,3};b_num2 = dci_dc{jj,4};conv_num = dci_dc{jj,2};
            angdci(jj,k)=theta(bus_int(b_num1),k)-theta(bus_int(b_num2),k);
            dcid_sig(jj,k)=(angdci(jj,k)-angdci(jj,k-1))/(t(k)-t(k-1));
            st_state = tot_states+1; dci_states = dci_dc{jj,7}; tot_states = tot_states+dci_states; 
            ydcimx=dci_dc{jj,5};ydcimn = dci_dc{jj,6};
            dci_dsig(jj,k) = ...
               dci_sud(jj,k,flag,dci_dc{jj,1},dcid_sig(jj,k),ydcimx,ydcimn,xdci_dc(st_state:tot_states,10*(k-1)+1));
         end
      end 
      dc_cont(0,k,10*(k-1)+1,bus_sim,flag);
      % network interface for control models
      dpwf(0,k,bus_sim,flag);
      mexc_sig(t(k),k);
      smpexc(0,k,bus_sim,flag);
      smppi(0,k,bus_sim,flag);
      exc_st3(0,k,bus_sim,flag);
      exc_dc12(0,k,bus_sim,flag);
      mtg_sig(t(k),k);
      tg(0,k,bus_sim,flag);
      tg_hydro(0,k,bus_sim,flag);
      if n_dcud~=0
         % set the new line currents
         for jj=1:n_dcud
            l_num = svc_dc{jj,3};svc_num = svc_dc{jj,2};
            from_bus = bus_int(line_sim(l_num,1)); to_bus=bus_int(line_sim(l_num,2));
            svc_bn = bus_int(svc_con(svc_num,2));
            V1 = bus_v(from_bus,k);
            V2 = bus_v(to_bus,k);
            R = line_sim(l_num,3);X=line_sim(l_num,4);B=line_sim(l_num,5);tap = line_sim(l_num,6);phi = line_sim(l_num,7);
            [l_if,l_it] = line_cur(V1,V2,R,X,B,tap,phi);
            if svc_bn == from_bus; d_sig(jj,k)=abs(l_if);elseif svc_bn==to_bus;d_sig(jj,k)=abs(l_it); end
         end
      end
      if n_tcscud~=0
         % set the new bus voltages
         for jj=1:n_tcscud
            b_num = tcsc_dc{jj,3};tcsc_num = tcsc_dc{jj,2};
            td_sig(jj,k)=abs(bus_v(bus_int(b_num),k));
         end
      end
      
      i_plot=k-plot_now;
      if i_plot == 10
         plot_now = k;
	     v_p(1:k)=abs(bus_v(bus_idx(1),1:k));
	     % plot the voltage of the faulted bus
	     plot(t(1:k),v_p(1:k),'r')
	     title(['Voltage Magnitude at ' num2str(bus(bus_idx(1),1)) ' ' dfile]);
	     xlabel('time (s)');
         drawnow
      end
      % step 3b: compute dynamics and integrate
      flag = 2;
      sys_freq(k) = 1.0;
      mpm_sig(t(k),k);
      mac_ind(0,k,bus_sim,flag);
      mac_igen(0,k,bus_sim,flag);
      mac_sub(0,k,bus_sim,flag); 
      mac_tra(0,k,bus_sim,flag);
      mac_em(0,k,bus_sim,flag);
      dpwf(0,k,bus_sim,flag);
      pss(0,k,bus_sim,flag);
      mexc_sig(t(k),k);
      smpexc(0,k,bus_sim,flag);
      smppi(0,k,bus_sim,flag);
      exc_st3(0,k,bus_sim,flag);
      exc_dc12(0,k,bus_sim,flag);
      mtg_sig(t(k),k);
      tg(0,k,bus_sim,flag);
      tg_hydro(0,k,bus_sim,flag);
      if n_svc~=0
         v_svc = abs(bus_v(bus_int(svc_con(:,2)),k));
         if n_dcud~=0
            tot_states=0;
            for jj = 1:n_dcud
               ysvcmx = svc_dc{jj,4};ysvcmn = svc_dc{jj,5};
               svc_num = svc_dc{jj,2};
               st_state = tot_states+1; svc_states = svc_dc{jj,6}; tot_states = tot_states+svc_states; 
               [svc_dsig(svc_num,k),xsvc_dc(st_state:tot_states,k),dxsvc_dc(st_state:tot_states,k)] =...
                  svc_sud(jj,k,flag,svc_dc{jj,1},d_sig(jj,k),ysvcmx,ysvcmn,xsvc_dc(st_state:tot_states,k));
            end
         end
         msvc_sig(t(k),k);
         svc(0,k,bus_sim,flag,v_svc);
      end
      if n_tcsc~=0
         if n_tcscud~=0
            tot_states=0;
            for jj = 1:n_tcscud
               ytcscmx = tcsc_dc{jj,4};ytcscmn = tcsc_dc{jj,5};
               tcsc_num = tcsc_dc{jj,2};
               st_state = tot_states+1; tcsc_states = tcsc_dc{jj,6}; tot_states = tot_states+tcsc_states; 
               [tcsc_dsig(tcsc_num,k),xtcsc_dc(st_state:tot_states,k),dxtcsc_dc(st_state:tot_states,k)] =...
                  tcsc_sud(jj,k,flag,tcsc_dc{jj,1},td_sig(jj,k),ytcscmx,ytcscmn,xtcsc_dc(st_state:tot_states,k));
            end
         end
         mtcsc_sig(t(k),k);
         tcsc(0,k,bus_sim,flag);
      end
      
      if n_lmod~=0
         ml_sig(t(k),k);
         lmod(0,k,bus_sim,flag);
      end
      if n_rlmod~=0
         rml_sig(t(k),k);
         rlmod(0,k,bus_sim,flag);
      end
      % integrate dc at ten times rate
      mdc_sig(t(k),k);
      if n_conv~=0
         hdc_sol = h_sol/10;
         for kk = 1:10
            kdc=10*(k-1)+kk;
            [xdcr_dc(:,kdc:kdc+1),dxdcr_dc(:,kdc:kdc+1),xdci_dc(:,kdc:kdc+1),dxdci_dc(:,kdc:kdc+1)] = ...
               dc_sim(k,kk,dcr_dc,dci_dc,xdcr_dc(:,kdc),xdci_dc(:,kdc),bus_sim,hdc_sol); 
         end
      else
         dc_cont(0,k,k,bus_sim,2);
         dc_line(0,k,k,bus_sim,2);
      end
      
      j = k+1;
     
      % following statements are predictor steps
      mac_ang(:,j) = mac_ang(:,k) + h_sol*dmac_ang(:,k); 
      mac_spd(:,j) = mac_spd(:,k) + h_sol*dmac_spd(:,k);
      edprime(:,j) = edprime(:,k) + h_sol*dedprime(:,k);
      eqprime(:,j) = eqprime(:,k) + h_sol*deqprime(:,k);
      psikd(:,j) = psikd(:,k) + h_sol*dpsikd(:,k);
      psikq(:,j) = psikq(:,k) + h_sol*dpsikq(:,k);
      Efd(:,j) = Efd(:,k) + h_sol*dEfd(:,k);
      V_R(:,j) = V_R(:,k) + h_sol*dV_R(:,k);
      V_As(:,j) = V_As(:,k) + h_sol*dV_As(:,k);
      R_f(:,j) = R_f(:,k) + h_sol*dR_f(:,k);
      V_TR(:,j) = V_TR(:,k) + h_sol*dV_TR(:,k);
      sdpw1(:,j) = sdpw1(:,k) + h_sol*dsdpw1(:,k);
      sdpw2(:,j) = sdpw2(:,k) + h_sol*dsdpw2(:,k);
      sdpw3(:,j) = sdpw3(:,k) + h_sol*dsdpw3(:,k);
      sdpw4(:,j) = sdpw4(:,k) + h_sol*dsdpw4(:,k);
      sdpw5(:,j) = sdpw5(:,k) + h_sol*dsdpw5(:,k);
      sdpw6(:,j) = sdpw6(:,k) + h_sol*dsdpw6(:,k);
      pss1(:,j) = pss1(:,k) + h_sol*dpss1(:,k);
      pss2(:,j) = pss2(:,k) + h_sol*dpss2(:,k);
      pss3(:,j) = pss3(:,k) + h_sol*dpss3(:,k);
      tg1(:,j) = tg1(:,k) + h_sol*dtg1(:,k);
      tg2(:,j) = tg2(:,k) + h_sol*dtg2(:,k);
      tg3(:,j) = tg3(:,k) + h_sol*dtg3(:,k);
      tg4(:,j) = tg4(:,k) + h_sol*dtg4(:,k);
      tg5(:,j) = tg5(:,k) + h_sol*dtg5(:,k);
      vdp(:,j) = vdp(:,k) + h_sol*dvdp(:,k);
      vqp(:,j) = vqp(:,k) + h_sol*dvqp(:,k);
      slip(:,j) = slip(:,k) + h_sol*dslip(:,k);
      vdpig(:,j) = vdpig(:,k) + h_sol*dvdpig(:,k);
      vqpig(:,j) = vqpig(:,k) + h_sol*dvqpig(:,k);
      slig(:,j) = slig(:,k) + h_sol*dslig(:,k);
      B_cv(:,j) = B_cv(:,k) + h_sol*dB_cv(:,k);
      B_con(:,j) = B_con(:,k) + h_sol*dB_con(:,k);
      xsvc_dc(:,j) = xsvc_dc(:,k) + h_sol* dxsvc_dc(:,k);
      B_tcsc(:,j) = B_tcsc(:,k) + h_sol*dB_tcsc(:,k);
      xtcsc_dc(:,j) = xtcsc_dc(:,k) + h_sol* dxtcsc_dc(:,k);
      lmod_st(:,j) = lmod_st(:,k) + h_sol*dlmod_st(:,k);
      rlmod_st(:,j) = rlmod_st(:,k)+h_sol*drlmod_st(:,k);
      flag = 1;
      % mach_ref(j) = mac_ang(syn_ref,j);
      mach_ref(j) = 0;
      % perform network interface calculations again with predicted states
      mpm_sig(t(j),j);
      mac_ind(0,j,bus_sim,flag);
      mac_igen(0,j,bus_sim,flag);
      mac_sub(0,j,bus_sim,flag);
      mac_tra(0,j,bus_sim,flag);
      mac_em(0,j,bus_sim,flag);
      % assume Vdc remains unchanged for first pass through dc controls interface
      mdc_sig(t(j),j);
      dc_cont(0,j,10*(j-1)+1,bus_sim,flag);
      
      % Calculate current injections and bus voltages and angles
      if j >= sum(k_inc(1:3))+1
         % fault cleared
         bus_sim = bus_pf2;
         bus_int = bus_intpf2;
         Y1 = Y_gpf2;
         Y2 = Y_gncpf2;
         Y3 = Y_ncgpf2;
         Y4 = Y_ncpf2;
         Vr1 = V_rgpf2; 
         Vr2 = V_rncpf2;
         bo = bopf2;
         h_sol = i_simu(j,ks,k_inc,h,bus_sim,Y1,Y2,Y3,Y4,Vr1,Vr2,bo);     
      elseif j >=sum(k_inc(1:2))+1
         % near bus cleared
         bus_sim = bus_pf1;
         bus_int = bus_intpf1;
         Y1 = Y_gpf1;
         Y2 = Y_gncpf1;
         Y3 = Y_ncgpf1;
         Y4 = Y_ncpf1;
         Vr1 = V_rgpf1; 
         Vr2 = V_rncpf1;
         bo = bopf1;
         h_sol = i_simu(j,ks,k_inc,h,bus_sim,Y1,Y2,Y3,Y4,Vr1,Vr2,bo);   
      elseif j>=k_inc(1)+1
         % fault applied
         bus_sim = bus_f;
         bus_int = bus_intf;
         Y1 = Y_gf;
         Y2 = Y_gncf;
         Y3 = Y_ncgf;
         Y4 = Y_ncf;
         Vr1 = V_rgf; 
         Vr2 = V_rncf;
         bo = bof;
         h_sol = i_simu(j,ks,k_inc,h,bus_sim,Y1,Y2,Y3,Y4,Vr1,Vr2,bo);     
      elseif k<k_inc(1)+1  % JHC - DKF thinks k should be j
         % pre fault
         bus_sim = bus;
         bus_int = bus_intprf;
         Y1 = Y_gprf;
         Y2 = Y_gncprf;
         Y3 = Y_ncgprf;
         Y4 = Y_ncprf;
         Vr1 = V_rgprf; 
         Vr2 = V_rncprf;
         bo = boprf;
         h_sol = i_simu(j,ks,k_inc,h,bus_sim,Y1,Y2,Y3,Y4,Vr1,Vr2,bo);     
      end
      vex(:,j)=vex(:,k);
      cur_ord(:,j) = cur_ord(:,k);
      % calculate the new value of bus angles rectifier user defined control
      if ndcr_ud~=0
         tot_states=0;
         for jj = 1:ndcr_ud
            b_num1 = dcr_dc{jj,3};b_num2 = dcr_dc{jj,4};conv_num = dcr_dc{jj,2};
            angdcr(jj,j)=theta(bus_int(b_num1),j)-theta(bus_int(b_num2),j);
            dcrd_sig(jj,j)=angdcr(jj,j);
            st_state = tot_states+1; dcr_states = dcr_dc{jj,7}; tot_states = tot_states+dcr_states; 
            ydcrmx=dcr_dc{jj,5};ydcrmn = dcr_dc{jj,6};
            dcr_dsig(jj,j) = ...
               dcr_sud(jj,j,flag,dcr_dc{jj,1},dcrd_sig(jj,j),ydcrmx,ydcrmn,xdcr_dc(st_state:tot_states,10*(j-1)+1));
         end
      end
      if ndci_ud~=0
         % calculate the new value of bus angles inverter user defined control
         for jj = 1:ndci_ud
            tot_states=0;
            b_num1 = dci_dc{jj,3};b_num2 = dci_dc{jj,4};conv_num = dci_dc{jj,2};
            angdci(jj,j)=theta(bus_int(b_num1),j)-theta(bus_int(b_num2),j);
            dcid_sig(jj,j)=(angdci(jj,j)-angdci(jj,k))/(t(j)-t(k));
            st_state = tot_states+1; dci_states = dci_dc{jj,7}; tot_states = tot_states+dci_states; 
            ydcimx=dci_dc{jj,5};ydcimn = dci_dc{jj,6};
            dci_dsig(jj,j) = ...
               dci_sud(jj,j,flag,dci_dc{jj,1},dcid_sig(jj,j),ydcimx,ydcimn,xdci_dc(st_state:tot_states,10*(j-1)+1));
         end
      end 
      dc_cont(0,j,10*(j-1)+1,bus_sim,flag);
      dpwf(0,j,bus_sim,flag);
      pss(0,j,bus_sim,flag);
      mexc_sig(t(j),j);
      smpexc(0,j,bus_sim,flag);
      smppi(0,j,bus_sim,flag);
      exc_st3(0,j,bus_sim,flag);
      exc_dc12(0,j,bus_sim,flag);
      tg(0,j,bus_sim,flag); 
      tg_hydro(0,j,bus_sim,flag);
      if n_dcud~=0
         % set the new line currents
         for jj=1:n_dcud
            l_num = svc_dc{jj,3};svc_num = svc_dc{jj,2};
            from_bus = bus_int(line_sim(l_num,1)); to_bus=bus_int(line_sim(l_num,2));
            svc_bn = bus_int(svc_con(svc_num,2));
            V1 = bus_v(from_bus,j);
            V2 = bus_v(to_bus,j);
            R = line_sim(l_num,3);X=line_sim(l_num,4);B=line_sim(l_num,5);tap = line_sim(l_num,6);phi = line_sim(l_num,7);
            [l_if,l_it] = line_cur(V1,V2,R,X,B,tap,phi);
            if svc_bn == from_bus; 
               d_sig(jj,j)=abs(l_if);            
            elseif svc_bn==to_bus;
               d_sig(jj,j)=abs(l_it);
            end
         end
      end
      if n_tcscud~=0
         % set the new line currents
         for jj=1:n_tcscud
            b_num = tcsc_dc{jj,3};tcsc_num = tcsc_dc{jj,2};
            td_sig(jj,j)=abs(bus_v(bus_int(b_num),j));
         end
      end
      
      
      flag = 2;
      mac_ind(0,j,bus_sim,flag);
      mac_igen(0,j,bus_sim,flag);
      mac_sub(0,j,bus_sim,flag);
      mac_tra(0,j,bus_sim,flag);
      mac_em(0,j,bus_sim,flag);
      dpwf(0,j,bus_sim,flag);
      pss(0,j,bus_sim,flag);
      mexc_sig(t(j),j);
      smpexc(0,j,bus_sim,flag);
      smppi(0,j,bus_sim,flag);
      exc_st3(0,j,bus_sim,flag);
      exc_dc12(0,j,bus_sim,flag);
      mtg_sig(t(j),j);
      tg(0,j,bus_sim,flag);
      tg_hydro(0,j,bus_sim,flag);
      if n_svc~=0
         msvc_sig(t(j),j);
         if n_dcud~=0
            tot_states=0;
            for jj = 1:n_dcud
               ysvcmx = svc_dc{jj,4};ysvcmn = svc_dc{jj,5};
               svc_num = svc_dc{jj,2};
               st_state = tot_states+1; svc_states = svc_dc{jj,6}; tot_states = tot_states+svc_states; 
               [svc_dsig(svc_num,j),xsvc_dc(st_state:tot_states,j),dxsvc_dc(st_state:tot_states,j)] =...
                  svc_sud(jj,j,flag,svc_dc{jj,1},d_sig(jj,j),ysvcmx,ysvcmn,xsvc_dc(st_state:tot_states,j));
            end
         end
         v_svc = abs(bus_v(bus_int(svc_con(:,2)),j));
         bus_sim = svc(0,j,bus_sim,flag,v_svc);
      end
      if n_tcsc~=0
         mtcsc_sig(t(j),j);
         if n_tcscud~=0
            tot_states=0;
            for jj = 1:n_tcscud
               ytcscmx = tcsc_dc{jj,4};ytcscmn = tcsc_dc{jj,5};
               tcsc_num = tcsc_dc{jj,2};
               st_state = tot_states+1; tcsc_states = tcsc_dc{jj,6}; tot_states = tot_states+tcsc_states; 
               [tcsc_dsig(tcsc_num,j),xtcsc_dc(st_state:tot_states,j),dxtcsc_dc(st_state:tot_states,j)] =...
                  tcsc_sud(jj,j,flag,tcsc_dc{jj,1},td_sig(jj,j),ytcscmx,ytcscmn,xtcsc_dc(st_state:tot_states,j));
            end
         end
         tcsc(0,j,bus_sim,flag);
      end
      
      if n_lmod~=0
         ml_sig(t(j),j);
         lmod(0,j,bus_sim,flag);
      end
      if n_rlmod~=0
         rml_sig(t(j),j);
         rlmod(0,j,bus_sim,flag);
      end
      if n_conv~=0
         hdc_sol = h_sol/10;
         for kk = 1:10
            jdc=10*(j-1)+kk;
            [xdcr_dc(:,jdc:jdc+1),dxdcr_dc(:,jdc:jdc+1),xdci_dc(:,jdc:jdc+1),dxdci_dc(:,jdc:jdc+1)] = ...
               dc_sim(j,kk,dcr_dc,dci_dc,xdcr_dc(:,jdc),xdci_dc(:,jdc),bus_sim,hdc_sol); 
         end
      else
         dc_cont(0,j,j,bus_sim,2);
         dc_line(0,j,j,bus_sim,2);
      end
      % following statements are corrector steps
      mac_ang(:,j) = mac_ang(:,k) +...
         h_sol*(dmac_ang(:,k)+dmac_ang(:,j))/2.;
      mac_spd(:,j) = mac_spd(:,k) +...
         h_sol*(dmac_spd(:,k)+dmac_spd(:,j))/2.;
      edprime(:,j) = edprime(:,k) +...
         h_sol*(dedprime(:,k)+dedprime(:,j))/2.;
      eqprime(:,j) = eqprime(:,k) +...
         h_sol*(deqprime(:,k)+deqprime(:,j))/2.;
      psikd(:,j) = psikd(:,k) +...
         h_sol*(dpsikd(:,k)+dpsikd(:,j))/2.;
      psikq(:,j) = psikq(:,k) +...
         h_sol*(dpsikq(:,k)+dpsikq(:,j))/2.;
      Efd(:,j) = Efd(:,k) +...
         h_sol*(dEfd(:,k)+dEfd(:,j))/2.;
      V_R(:,j) = V_R(:,k) +...
         h_sol*(dV_R(:,k)+dV_R(:,j))/2.;
      V_As(:,j) = V_As(:,k) +...
         h_sol*(dV_As(:,k)+dV_As(:,j))/2.;
      R_f(:,j) = R_f(:,k) +...
         h_sol*(dR_f(:,k)+dR_f(:,j))/2.;
      V_TR(:,j) = V_TR(:,k) +...
         h_sol*(dV_TR(:,k)+dV_TR(:,j))/2.;
      sdpw11(:,j) = sdpw1(:,k) +h_sol*(dsdpw1(:,k)+dsdpw1(:,j))/2.;
      sdpw12(:,j) = sdpw2(:,k) +h_sol*(dsdpw2(:,k)+dsdpw2(:,j))/2.;
      sdpw13(:,j) = sdpw3(:,k) +h_sol*(dsdpw3(:,k)+dsdpw3(:,j))/2.;
      sdpw14(:,j) = sdpw4(:,k) +h_sol*(dsdpw4(:,k)+dsdpw4(:,j))/2.;
      sdpw15(:,j) = sdpw5(:,k) +h_sol*(dsdpw5(:,k)+dsdpw5(:,j))/2.;
      sdpw16(:,j) = sdpw6(:,k) +h_sol*(dsdpw6(:,k)+dsdpw6(:,j))/2.;
      pss1(:,j) = pss1(:,k) +h_sol*(dpss1(:,k)+dpss1(:,j))/2.;
      pss2(:,j) = pss2(:,k) +h_sol*(dpss2(:,k)+dpss2(:,j))/2.;
      pss3(:,j) = pss3(:,k) +h_sol*(dpss3(:,k)+dpss3(:,j))/2.;
      tg1(:,j) = tg1(:,k) + h_sol*(dtg1(:,k) + dtg1(:,j))/2.;
      tg2(:,j) = tg2(:,k) + h_sol*(dtg2(:,k) + dtg2(:,j))/2.;
      tg3(:,j) = tg3(:,k) + h_sol*(dtg3(:,k) + dtg3(:,j))/2.;
      tg4(:,j) = tg4(:,k) + h_sol*(dtg4(:,k) + dtg4(:,j))/2.;
      tg5(:,j) = tg5(:,k) + h_sol*(dtg5(:,k) + dtg5(:,j))/2.;
      vdp(:,j) = vdp(:,k) + h_sol*(dvdp(:,j) + dvdp(:,k))/2.;
      vqp(:,j) = vqp(:,k) + h_sol*(dvqp(:,j) + dvqp(:,k))/2.;
      slip(:,j) = slip(:,k) + h_sol*(dslip(:,j) + dslip(:,k))/2.;
      vdpig(:,j) = vdpig(:,k) + h_sol*(dvdpig(:,j) + dvdpig(:,k))/2.;
      vqpig(:,j) = vqpig(:,k) + h_sol*(dvqpig(:,j) + dvqpig(:,k))/2.;
      slig(:,j) = slig(:,k) + h_sol*(dslig(:,j) + dslig(:,k))/2.;
      B_cv(:,j) = B_cv(:,k) + h_sol*(dB_cv(:,j) + dB_cv(:,k))/2.;
      B_con(:,j) = B_con(:,k) + h_sol*(dB_con(:,j) + dB_con(:,k))/2.;
      xsvc_dc(:,j) = xsvc_dc(:,k) + h_sol*(dxsvc_dc(:,j) + dxsvc_dc(:,k))/2.;
      B_tcsc(:,j) = B_tcsc(:,k) + h_sol*(dB_tcsc(:,j) + dB_tcsc(:,k))/2.;
      xtcsc_dc(:,j) = xtcsc_dc(:,k) + h_sol*(dxtcsc_dc(:,j) + dxtcsc_dc(:,k))/2.;
      lmod_st(:,j) = lmod_st(:,k) + h_sol*(dlmod_st(:,j) + dlmod_st(:,k))/2.;
      rlmod_st(:,j) = rlmod_st(:,k) + h_sol*(drlmod_st(:,j) + drlmod_st(:,k))/2.;
   end
   kt = kt + k_inc(ks);
   ks = ks+1; 
end
V1 = bus_v(bus_int(line(:,1)),:);
V2 = bus_v(bus_int(line(:,2)),:);
R = line(:,3); X = line(:,4); B = line(:,5);
tap = line(:,6); phi = line(:,7);
[ilf,ilt]=line_cur(V1,V2,R,X,B,tap,phi);%line currents
et = toc;
ets = num2str(et);
disp(['elapsed time = ' ets 's'])

flag = 0;
while(flag == 0)
   disp('You can examine the system response')
   disp('Type 1 to see all machine angles in 3D')
   disp('     2 to see all machine speed deviation in 3D')
   disp('     3 to see all machine turbine powers')
   disp('     4 to see all machine electrical powers')
   disp('     5 to see all field voltages')
   disp('     6 to see all bus voltage magnitude in 3D')
   disp('     7 to see the line power flows')
   disp('     0 to quit and plot your own curves')
   sel = input('enter selection >> ');
   if isempty(sel);sel = 0;end
   if sel == 1
      mesh(t,[1:1:n_mac],mac_ang*180/pi)
      title('machine angles')
      xlabel('time in seconds')
      ylabel('internal generator number')
      zlabel('angle in degrees')
      disp('paused: press any key to continue')
      pause
   elseif sel == 2
      lt = length(t);
      mesh(t,[1:1:n_mac],mac_spd- ones(n_mac,lt))
      title('machine speed deviations')
      xlabel('time in seconds')
      ylabel('internal generator number')
      zlabel('speed in pu')
      disp('paused: press any key to continue')
      pause
   elseif sel == 3
      plot(t,pmech)
      title('turbine power')
      xlabel('time in seconds')
      ylabel('power in pu on generator base')
      disp('paused: press any key to continue')
      pause
   elseif sel == 4
      plot(t,pelect)
      title('generator electrical power')
      xlabel('time in seconds')
      ylabel('power in pu on system base')
      disp('paused: press any key to continue')
      pause
   elseif sel == 5
      plot(t,Efd)
      title('Generator field voltages')
      xlabel('time in seconds')
      ylabel('voltage in pu')
      disp('paused: press any key to continue')
      pause
   elseif sel == 6
      [nbus,dum] = size(bus_v);
      mesh(t,(1:1:nbus),abs(bus_v))
      xlabel('time in seconds')
      ylabel('bus number')
      zlabel('voltage in pu')
      disp('paused: press any key to continue')
      pause
   elseif sel == 7
      nline = length(line(:,1));
      if nline<50
         [S1,S2] = line_pq(V1,V2,R,X,B,tap,phi);
         plot(t,real(S1));
         xlabel('time in seconds')
         ylabel(' power flow per unit')
         disp('paused: press any key to continue')
         pause
      else
         % ask for lines to be plotted
         line_range = input(' enter a rangle of lines less than 50, ie 1:49   >');
         if isempty(line_range);line_range = 1:49;end
         V1 = bus_v(bus_int(line(line_range,1)),:);
         V2 = bus_v(bus_int(line(line_range,2)),:);
         R = line(line_range,3); X = line(line_range,4); 
         B = line(line_range,5);
         tap = line(line_range,6); phi = line(line_range,7);
         [S1,S2] = line_pq(V1,V2,R,X,B,tap,phi);
         plot(t,real(S1));
         xlabel('time in seconds')
         ylabel(' power flow per unit')
         disp('paused: press any key to continue')
         pause
      end
   elseif sel == 0
      flag = 1;
   else
      error('invalid selection')
   end
end
t_dc=t_dc(1:length(t_dc)-10);
i_dc=i_dc(:,1:length(t_dc));
i_dcr=i_dcr(:,1:length(t_dc));
i_dci=i_dci(:,1:length(t_dc));
alpha=alpha(:,1:length(t_dc));
gamma=gamma(:,1:length(t_dc));
Vdc=Vdc(:,1:length(t_dc));
v_conr=v_conr(:,1:length(t_dc));
v_coni= v_coni(:,1:length(t_dc));
dv_conr= dv_conr(:,1:length(t_dc));
dv_coni= dv_coni(:,1:length(t_dc));
xdcr_dc= xdcr_dc(:,1:length(t_dc));
xdci_dc= xdci_dc(:,1:length(t_dc));
dxdcr_dc= dxdcr_dc(:,1:length(t_dc));
dxdci_dc= dxdci_dc(:,1:length(t_dc));
dv_dcc= dv_dcc(:,1:length(t_dc));
v_dcc= v_dcc(:,1:length(t_dc));
di_dci= di_dci(:,1:length(t_dc));
di_dcr=di_dcr(:,1:length(t_dc));

% tidy workspace
clear B H_sum IHT  R  SHT   VLT          
clear V_rgf V_rgpf1 V_rgpf2 V_rgprf V_rncf V_rncpf1 V_rncpf2 V_rncprf Vdc_ref      
clear Vr1 Vr2 X Y1 Y2 Y3 Y4           
clear Y_gf Y_gncf Y_gncpf1 Y_gncpf2 Y_gncprf Y_gpf1       
clear Y_gpf2 Y_gprf Y_ncf Y_ncgf Y_ncgpf1 Y_ncgpf2 Y_ncgprf Y_ncpf1 Y_ncpf2 Y_ncprf      
clear ac_bus ac_line ans b_num1 b_num2  bo bof bopf1 bopf2 boprf                  
clear bus_ang bus_f bus_idx bus_int bus_intf bus_intpf1 bus_intpf2 bus_intprf   
clear bus_pf1 bus_pf2 bus_sim  cap_idx conv_num  dc2_idx dc_TA  dc_TA_idx    
clear dc_TB dc_TB_idx dc_TE dc_TE_idx dc_TF dc_TF_idx dc_TR dc_TR_idx    
clear dc_idx dc_noTB_idx dc_noTE_idx dc_noTR_idx dc_pot dcc_con      
clear dcc_pot dci_dc dciud_idx dcl_con  dcli_idx dcr_dc dcr_states dcrud_idx    
clear dcsp_con dcud_idx dfile dpw_Td_idx dpw_con dpw_idx dpw_mb_idx   
%clear dpw_pot dpw_pss_idx dtcscud_idx  dummy et ets exc_con exc_pot      
clear dpw_pot dpw_pss_idx dtcscud_idx  dummy et ets 
      % JHC 1/19/2017 keep exc_pot for Vref and exc_con for K_E     
clear f f_nearbus f_type flag h h_sol i i_aci i_acr i_dcinj      
clear i_idx i_plot ibus_con idig idmot igbus igen_con igen_int igen_pot ind_con      
clear ind_int ind_pot inv_ac_bus inv_ac_line inv_par iqig iqmot        
clear j jay jj k k_end k_inc k_start k_tot ks kt           
clear ktmax l_cap l_no_cap lfile line_f line_flw line_par line_pf1 line_pf2 line_sim     
clear lmod_con lmod_idx lmod_pot lmon_con  load_con load_pot lswitch lt mac_con mac_em_idx   
clear mac_ib_em mac_ib_idx mac_ib_sub mac_ib_tra mac_int  mac_pot mac_sub_idx mac_tra_idx  
clear mld_con motbus n n_bus n_conv n_dc n_dc2 n_dcl n_dcud n_dpw n_em n_exc n_ib n_ib_em      
clear n_ib_sub n_ib_tra n_ig n_lmod n_mac n_mot n_pm n_pss n_rlmod n_smp n_st3 n_sub n_svc        
clear n_switch n_tcsc n_tcscud n_tg n_tgh n_tra nbus         
clear ndci_ud ndcr_ud netg_con ngm nload no_cap_idx no_ind_idx not_ib_idx   
clear ntot pathname  phi pig plot_now psidpp psiqpp pss_T pss_T2 pss_T4  pss_T4_idx pss_con      
clear pss_exc_idx pss_idx pss_mb_idx pss_noT4_idx pss_p_idx pss_pot pss_sp_idx r_idx        
clear rec_ac_bus rec_ac_line rec_num rec_par ric_idx rlmod_con rlmod_idx rlmod_pot rpc_idx scr_con sel          
clear smp_TA smp_TA_idx smp_TB smp_TB_idx smp_TR smp_TR_idx smp_idx smp_noTA_idx smp_noTB_idx 
clear smp_no_TR_idx st3_TA st3_TA_idx st3_TB st3_TB_idx st3_TR st3_TR_idx st3_idx st3_noTA_idx 
clear st3_noTB_idx st3_noTR_idx st_state stab_con sv svc_con svc_dc       
clear svc_idx svc_pot svcll_idx sw_con sw_count syn_ref t_init t_switch tap          
clear tapi tapr tcsc_con tcsc_dc tcsct_idx  tcsvf_idx tg_con tg_idx tg_pot       
clear tg_sig tgh_idx timestep tload tmax tmaxi tmaxr tmin tmini tminr        
clear tot_states tstep tstepi tstepr tswitch vdig vdmot vnc vqig vqmot ydcrmn ydcrmx  
clear z z1 z_dpw z_pss z_tg zdc zdcl ze zig zm 
