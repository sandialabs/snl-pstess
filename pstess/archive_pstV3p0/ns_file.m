% This m-file checks the dynamic data indices and determines the 
% total number of states for each device
% Called from svm_mgen
% 5:39 pm 29/12/98
% Author: Graham Rogers
% Modified: December 1998
% tcsc model added
% Modified: August 1997
% Induction Generators added
% Modified: August 1997
% Load modulation added
% Date: September 1996
% Copyright: Joe Chow/ Cherry Tree Scientific Software 1991 to 1997

for k = 1:n_mac
   % generators
   no_mac = 0;
   if ~isempty(mac_ib_idx);
      k_ib = find(mac_ib_idx==k);
      if ~isempty(k_ib); 
         state(k) = 0;
         no_mac = 1;
      end;
   end
   if no_mac ==0
      if ~isempty(mac_em_idx);
         k_em = find(mac_em_idx == k);
         if ~isempty(k_em); 
            state(k) = 2; 
         end;
      end
      if ~isempty(mac_tra_idx);
         k_tra = find(mac_tra_idx==k);
         if ~isempty(k_tra);
            state(k) = 4; 
         end;
      end;
      if ~isempty(mac_sub_idx);
         k_sub = find(mac_sub_idx == k);
         if ~isempty(k_sub);
            state(k)=6;
         end;
      end;
   end
   gen_state(k) = state(k);
   % exciters
   if ~isempty(exc_con)
      s_TR = 0; s_TB = 0; s_TA=0; s_TE=0;
      if n_smp~=0
         k_smp = find(mac_int(mac_exc(smp_idx))==k);
         if ~isempty(k_smp)
            k_exc = find(mac_exc(smp_idx)==k);
            if ~isempty(smp_TR_idx);s_TR = sum(smp_TR_idx == k_exc);end
            if ~isempty(smp_TB_idx);s_TB = sum(smp_TB_idx == k_exc);end
            if ~isempty(smp_TA_idx);s_TA = sum(smp_TA_idx == k_exc);end
            TR_state(k) = s_TR;
            TB_state(k) = s_TB;
            Efd_state(k) = s_TA;
            state(k) = state(k) + s_TR + s_TB +s_TA;
         end
      end
      if n_smppi~=0
         k_smppi = find(mac_int(mac_exc(smppi_idx))==k);
         if ~isempty(k_smppi)
            k_exc = find(mac_exc(smppi_idx)==k);
            if ~isempty(smppi_TR_idx);s_TR = sum(smppi_TR_idx == k_exc);end
            TR_state(k) = s_TR;
            TB_state(k) = 1;
            Efd_state(k) = 1;
            state(k) = state(k) + s_TR + 2;
         end
      end
      if n_dc~=0
         k_dc = find(mac_int(mac_exc(dc_idx))==k);
         if ~isempty(k_dc)
            k_exc = find(mac_exc(dc_idx)==k);
            if ~isempty(dc_TR_idx);s_TR = sum(dc_TR_idx == k_exc);end
            if ~isempty(dc_TB_idx);s_TB = sum(dc_TB_idx == k_exc);end
            if ~isempty(dc_TA_idx);s_TA = sum(dc_TA_idx == k_exc);end
            if ~isempty(dc_TE_idx);s_TE = sum(dc_TE_idx == k_exc);end
            TR_state(k) = s_TR;
            TB_state(k) = s_TB;
            TA_state(k) = s_TA;
            Efd_state(k) = s_TE;
            R_f_state(k) = 1;
            state(k) = state(k) + s_TR + s_TB +s_TA + s_TE + 1;
         end
      end
      if n_st3 ~=0 
         k_st3 = find(mac_int(mac_exc(st3_idx))==k);
         if ~isempty(k_st3)~=0
            k_exc = find(mac_exc(st3_idx)==k);
            if ~isempty(st3_TR_idx);s_TR = sum(st3_TR_idx == k_exc);end
            if ~isempty(st3_TB_idx);s_TB = sum(st3_TB_idx == k_exc);end
            if ~isempty(st3_TA_idx);s_TA = sum(st3_TA_idx == k_exc);end
            TR_state(k) = s_TR;
            TB_state(k) = s_TB;
            TA_state(k) = s_TA;
            state(k) = state(k) + s_TR + s_TB +s_TA;
         end
      end
   end
   %pss
   if ~isempty(pss_con)
      k_pss = find(mac_int(mac_pss)==k);
      s_T4 = 0;
      if ~isempty(k_pss)
         k_p = find(mac_pss==k);
         if ~isempty(pss_T4_idx); s_T4 =sum( pss_T4_idx == k_p);else;s_T4=0;end
         pss1_state(k) = 1;
         pss2_state(k) = 1;
         pss3_state(k) = s_T4;
         state(k) = state(k) + s_T4 + 2;
      end
   end
   %deltaP/omega filter
   if ~isempty(dpw_con)
      k_dpw = find(mac_int(mac_dpw)==k);
      if ~isempty(k_dpw)
         state(k) = state(k) + 6;
         dpw_state(k)=1;
      end
   end

   % turbine governors
   if ~isempty(tg_con)
      k_tg = find(mac_int(mac_tg)==k);
      if ~isempty(k_tg)
         state(k) = state(k) + 3;
         tg_state(k) = 1;
      end
      k_tgh = find(mac_int(mac_tgh)==k);
      if ~isempty(k_tgh)
         state(k) = state(k) + 5;
         tgh_state(k) = 1;
      end
   end
end 
% induction motor
n_mot_states = 0;
if n_mot~=0
   state_mot(1:n_mot) = 3*ones(n_mot,1);
   n_mot_states = sum(state_mot);
   state(n_mac+1:ngm) = 3*ones(n_mot,1);
end 

% induction generator
n_ig_states = 0;
if n_ig~=0
   state_ig(1:n_ig) = 3*ones(n_ig,1);
   n_ig_states = sum(state_ig);
   state(ngm+1:ntot) = 3*ones(n_ig,1);
end 

% svc
n_svc_states = 0;
if n_svc~=0
   state_svc(1:n_svc) = ones(n_svc,1);
   state_svccon = 0;
   sc = zeros(n_svc,1);
   if ~isempty(svcll_idx)
      state_svccon = ones(length(svcll_idx),1);
      sc(svcll_idx)=state_svccon;
   end
   n_svc_states = sum(state_svc)+sum(state_svccon);
   n_svc1 = ntot;
   state(n_svc1+1:n_svc1+n_svc) = ones(n_svc,1)+sc;
end;
% tcsc
n_tcsc_states = 0;
if n_tcsc~=0
   state_tcsc(1:n_tcsc) = ones(n_tcsc,1);
   n_svc_states = sum(state_tcsc);
   n_tcsc1 = ntot+n_svc;
   state(n_tcsc1+1:n_tcsc1+n_tcsc) = ones(n_tcsc,1);
end;

% lmod
n_lmod_states = 0;
if n_lmod~=0
   state_lmod(1:n_lmod) = ones(n_lmod,1);
   n_lmod_states = sum(state_lmod);
   n_lmod1 = ntot+n_svc+n_tcsc;
   state(n_lmod1+1:n_lmod1+n_lmod) = ones(n_lmod,1);
end;
% rlmod
n_rlmod_states = 0;
if n_rlmod~=0
   state_rlmod(1:n_rlmod) = ones(n_rlmod,1);
   n_rlmod_states = sum(state_rlmod);
   n_rlmod1 = ntot+n_svc+n_tcsc+n_lmod;
   state(n_rlmod1+1:n_rlmod1+n_rlmod) = ones(n_rlmod,1);
end;

% HVDC
n_hvdc_states = 0;
if n_conv~= 0
   state_hvdc(1:n_dcl) = 3*ones(n_dcl,1);
   if l_cap ~=0
      state_hvdc(cap_idx) =state_hvdc(cap_idx) + 2*ones(l_cap,1);
   end
   n_hvdc_states = sum(state_hvdc);
   n_hvdc1 = ntot +n_svc + n_tcsc+ n_lmod+n_rlmod;
   state(n_hvdc1+1:n_hvdc1+n_dcl) = (3 + 2*l_cap)*ones(n_dcl,1);
end
