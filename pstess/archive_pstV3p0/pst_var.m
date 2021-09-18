% file: pst_var.m  
% 9:58 am 10/6/99
% Syntax: pst_var
%
% Purpose: Define global variables for power system 
%          simulation.

%


% (c) Copyright 1991-1999 Joe H. Chow/Cherry Tree Scientific Software
%     All Rights Reserved

% History (in reverse chronological order)
% 
% version 1.6
% date: December 1998
% author: Graham Rogers
% modification: added TCSC

% version 1.5
% date:   July 1998
% author: Graham Rogers
% modification: add deltaP/omega filter

% version 1.4
% date:   June 1998
% author: Graham Rogers
% modification: add hydraulic turbine/governor model

% Version: 1.3
% Date:         August 1997         
% Author:       Graham Rogers	
% Modification: Added load modulation variables

% Version: 1.2
% Date:         November 1996
% Author:       Graham Rogers	
% Modification: Added dc variables

% Version: 1.1
% Date:    November 1995
% Author:  Graham Rogers	
% Purpose: Added Induction Motor Variables
% Modification:

% Version:  1.0
% Author:   Joe H. Chow
% Date:     January 1991

% system variables
global  basmva basrad syn_ref mach_ref sys_freq
global  bus_v bus_ang psi_re psi_im cur_re cur_im bus_int
global  lmon_con

% synchronous machine variables
global  mac_con mac_pot mac_int ibus_con
global  mac_ang mac_spd eqprime edprime psikd psikq
global  curd curq curdg curqg fldcur
global  psidpp psiqpp vex eterm theta ed eq 
global  pmech pelect qelect
global  dmac_ang dmac_spd deqprime dedprime dpsikd dpsikq
global  n_mac n_em n_tra n_sub n_ib
global  mac_em_idx mac_tra_idx mac_sub_idx mac_ib_idx not_ib_idx
global  mac_ib_em mac_ib_tra mac_ib_sub n_ib_em n_ib_tra n_ib_sub
  

% excitation system variables
global  exc_con exc_pot n_exc
global  Efd V_R V_A V_As R_f V_FB V_TR V_B
global  dEfd dV_R dV_As dR_f dV_TR
global  exc_sig pm_sig n_pm
global smp_idx n_smp dc_idx n_dc  dc2_idx n_dc2 st3_idx n_st3;
global smppi_idx n_smppi smppi_TR smppi_TR_idx smppi_no_TR_idx ;
global smp_TA smp_TA_idx smp_noTA_idx smp_TB smp_TB_idx smp_noTB_idx;
global smp_TR smp_TR_idx smp_no_TR_idx ;
global dc_TA dc_TA_idx dc_noTR_idx dc_TB dc_TB_idx dc_noTB_idx;
global dc_TE  dc_TE_idx dc_noTE_idx;
global dc_TF dc_TF_idx dc_TR dc_TR_idx 
global st3_TA st3_TA_idx st3_noTA_idx st3_TB st3_TB_idx st3_noTB_idx;
global st3_TR st3_TR_idx st3_noTR_idx;

% non-conforming load variables
global  load_con load_pot nload

% induction motor variables
global  tload t_init p_mot q_mot vdmot vqmot  idmot iqmot ind_con ind_pot
global  motbus ind_int mld_con n_mot t_mot
% states
global  vdp vqp slip 
% dstates
global dvdp dvqp dslip 

% induction genertaor variables
global  tmig  pig qig vdig vqig  idig iqig igen_con igen_pot
global  igen_int igbus n_ig
%states
global  vdpig vqpig slig 
%dstates
global dvdpig dvqpig dslig

% svc variables
global  svc_con n_svc svc_idx svc_pot svcll_idx
global  svc_sig
% svc user defined damping controls
global n_dcud dcud_idx svc_dsig
%states
global B_cv B_con
%dstates
global dB_cv dB_con

% tcsc variables
global  tcsc_con n_tcsc tcsvf_idx tcsct_idx 
global  B_tcsc dB_tcsc 
global  tcsc_sig tcsc_dsig
global  n_tcscud dtcscud_idx  %user defined damping controls

% load modulation variables
global  lmod_con n_lmod lmod_idx
global  lmod_pot lmod_st dlmod_st
global  lmod_sig
% reactive load modulation variables
global  rlmod_con n_rlmod rlmod_idx
global  rlmod_pot rlmod_st drlmod_st
global  rlmod_sig


% pss variables
global  pss_con pss_pot pss_mb_idx pss_exc_idx
global  pss1 pss2 pss3 dpss1 dpss2 dpss3 pss_out
global  pss_idx n_pss pss_sp_idx pss_p_idx;
global  pss_T  pss_T2 pss_T4 pss_T4_idx  pss_noT4_idx;

% DeltaP/omega filter variables
global  dpw_con dpw_out dpw_pot dpw_pss_idx dpw_mb_idx dpw_idx n_dpw dpw_Td_idx dpw_Tz_idx
global  sdpw1 sdpw2 sdpw3 sdpw4 sdpw5 sdpw6
global  dsdpw1 dsdpw2 dsdpw3 dsdpw4 dsdpw5 dsdpw6 

% turbine-governor variables
global  tg_con tg_pot 
global  tg1 tg2 tg3 tg4 tg5 dtg1 dtg2 dtg3 dtg4 dtg5
global  tg_idx  n_tg tg_sig tgh_idx n_tgh

%HVDC link variables
global  dcsp_con  dcl_con  dcc_con
global  r_idx  i_idx n_dcl  n_conv  ac_bus rec_ac_bus  inv_ac_bus
global  inv_ac_line  rec_ac_line ac_line dcli_idx
global  tap tapr tapi tmax tmin tstep tmaxr tmaxi tminr tmini tstepr tstepi
global  Vdc  i_dc P_dc i_dcinj dc_pot alpha gamma VHT dc_sig  cur_ord dcr_dsig dci_dsig
global  ric_idx  rpc_idx Vdc_ref dcc_pot 
global  no_cap_idx  cap_idx  no_ind_idx  l_no_cap  l_cap
global  ndcr_ud ndci_ud dcrud_idx dciud_idx dcrd_sig dcid_sig


% States
%line
global i_dcr i_dci  v_dcc 
global di_dcr  di_dci  dv_dcc 
%rectifier
global v_conr dv_conr  
%inverter
global v_coni dv_coni



% simulation control
global sw_con  scr_con

% pss design
global ibus_con  netg_con  stab_con