%-----------------------------------------------------------------------------%
% List of global variables for power system simulation
%
% Power and Energy Storage Systems Toolbox, PSTess, v1.1
%-----------------------------------------------------------------------------%

% bus // algebraic network variables
% g.bus.n_bus
% g.bus.bus_v
% g.bus.theta
% g.bus.bus_int

% dc // hvdc variables
% g.dc.dcsp_con
% g.dc.dcl_con
% g.dc.dcc_con
% g.dc.dcc_pot
% g.dc.dc_pot
% g.dc.ac_bus
% g.dc.rec_ac_bus
% g.dc.inv_ac_bus
% g.dc.ac_line
% g.dc.rec_ac_line
% g.dc.inv_ac_line
% g.dc.n_conv
% g.dc.n_dcl
% g.dc.r_idx
% g.dc.i_idx
% g.dc.dcli_idx
% g.dc.ldc_idx
% g.dc.ric_idx
% g.dc.rpc_idx
% g.dc.l_cap
% g.dc.l_no_cap
% g.dc.cap_idx
% g.dc.no_cap_idx
% g.dc.no_ind_idx
% g.dc.alpha
% g.dc.gamma
% g.dc.Vdc
% g.dc.Vdc_ref
% g.dc.i_dc
% g.dc.cur_ord
% g.dc.i_dci
% g.dc.i_dcr
% g.dc.v_dcc
% g.dc.v_conr
% g.dc.v_coni
% g.dc.dv_dcc
% g.dc.di_dcr
% g.dc.di_dci
% g.dc.dv_conr
% g.dc.dv_coni
% g.dc.tapr
% g.dc.tapi
% g.dc.tminr
% g.dc.tmini
% g.dc.tmaxr
% g.dc.tmaxi
% g.dc.tstepr
% g.dc.tstepi
% g.dc.ndcr_ud
% g.dc.ndci_ud
% g.dc.dcrud_idx
% g.dc.dciud_idx
% g.dc.dc_sig
% g.dc.dcr_dsig
% g.dc.dci_dsig

% dpw // delta-p-omega pss filters
% g.dpw.dpw_pss_idx
% g.dpw.dpw_con
% g.dpw.dpw_pot
% g.dpw.n_dpw
% g.dpw.dpw_Td_idx
% g.dpw.dpw_Tz_idx
% g.dpw.dpw_mb_idx
% g.dpw.dpw_out
% g.dpw.sdpw1
% g.dpw.sdpw2
% g.dpw.sdpw3
% g.dpw.sdpw4
% g.dpw.sdpw5
% g.dpw.sdpw6
% g.dpw.dsdpw1
% g.dpw.dsdpw2
% g.dpw.dsdpw3
% g.dpw.dsdpw4
% g.dpw.dsdpw5
% g.dpw.dsdpw6

% ess // converter-interfaced energy storage systems
% g.ess.n_ess
% g.ess.ess_idx
% g.ess.ess_con
% g.ess.ess_pot
% g.ess.ess_soc
% g.ess.ess_vmag
% g.ess.ess_vmag_pade
% g.ess.ess_sinj
% g.ess.ess_scmd
% g.ess.ess_iord
% g.ess.ess_cur
% g.ess.ess1--5
% g.ess.dess1--5
% g.ess.n_essud
% g.ess.dessud_idx
% g.ess.essud_con
% g.ess.ess_sig
% g.ess.ess_dsig

% exc // exciters
% g.exc.exc_con
% g.exc.exc_pot
% g.exc.n_exc
% g.exc.V_A
% g.exc.V_B
% g.exc.V_FB
% g.exc.Efd
% g.exc.R_f
% g.exc.V_R
% g.exc.V_As
% g.exc.V_TR
% g.exc.dEfd
% g.exc.dR_f
% g.exc.dV_R
% g.exc.dV_As
% g.exc.dV_TR
% g.exc.exc_sig
% g.exc.pss_out
% g.exc.n_smp
% g.exc.n_smppi
% g.exc.smp_idx
% g.exc.smppi_idx
% g.exc.smp_TA
% g.exc.smp_TA_idx
% g.exc.smp_noTA_idx
% g.exc.smp_TB
% g.exc.smp_TB_idx
% g.exc.smp_noTB_idx
% g.exc.smp_TR
% g.exc.smp_TR_idx
% g.exc.smp_noTR_idx
% g.exc.smppi_TR
% g.exc.smppi_TR_idx
% g.exc.smppi_noTR_idx
% g.exc.n_dc
% g.exc.n_dc1
% g.exc.n_dc2
% g.exc.dc_idx
% g.exc.dc1_idx
% g.exc.dc2_idx
% g.exc.dc_TA
% g.exc.dc_TA_idx
% g.exc.dc_noTA_idx
% g.exc.dc_TB
% g.exc.dc_TB_idx
% g.exc.dc_noTB_idx
% g.exc.dc_TE
% g.exc.dc_TE_idx
% g.exc.dc_noTE_idx
% g.exc.dc_TF
% g.exc.dc_TF_idx
% g.exc.dc_TR
% g.exc.dc_TR_idx
% g.exc.dc_noTR_idx
% g.exc.n_st3
% g.exc.st3_idx
% g.exc.st3_TA
% g.exc.st3_TA_idx
% g.exc.st3_noTA_idx
% g.exc.st3_TB
% g.exc.st3_TB_idx
% g.exc.st3_noTB_idx
% g.exc.st3_TR
% g.exc.st3_TR_idx
% g.exc.st3_noTR_idx

% freq // frequency measurement filters
% g.freq.bus_freq
% g.freq.bf_hpf
% g.freq.dbf_hpf
% g.freq.kpx
% g.freq.o2p
% g.freq.dx1_snlf
% g.freq.dx2_snlf
% g.freq.x2_snlf
% g.freq.x1_snlf
% g.freq.bus_freqf
% g.freq.bus_freqsnl

% gfma // droop-based grid-forming inverter controls
% g.gfma.dgfma1
% g.gfma.dgfma2
% g.gfma.dgfma3
% g.gfma.dgfma4
% g.gfma.dgfma5
% g.gfma.dgfma6
% g.gfma.dgfma7
% g.gfma.dgfma8
% g.gfma.dgfma9
% g.gfma.dgfma10
% g.gfma.gfma1
% g.gfma.gfma2
% g.gfma.gfma3
% g.gfma.gfma4
% g.gfma.gfma5
% g.gfma.gfma6
% g.gfma.gfma7
% g.gfma.gfma8
% g.gfma.gfma9
% g.gfma.gfma10
% g.gfma.gfma_con
% g.gfma.gfma_sig
% g.gfma.ivm_idx
% g.gfma.lim_flag
% g.gfma.mac_idx
% g.gfma.pset
% g.gfma.qset
% g.gfma.vset
% g.gfma.vref
% g.gfma.n_gfma

% igen // induction generators
% g.igen.tmig
% g.igen.igen_con
% g.igen.igen_pot
% g.igen.n_ig
% g.igen.igbus
% g.igen.igen_int
% g.igen.vdig
% g.igen.vqig
% g.igen.idig
% g.igen.iqig
% g.igen.pig
% g.igen.qig
% g.igen.s_igen
% g.igen.vdpig
% g.igen.vqpig
% g.igen.slig
% g.igen.dvdpig
% g.igen.dvqpig
% g.igen.dslig

% ind // induction motors
% g.ind.tload
% g.ind.t_init
% g.ind.ind_con
% g.ind.ind_pot
% g.ind.mld_con
% g.ind.n_mot
% g.ind.sat_idx
% g.ind.db_idx
% g.ind.dbc_idx
% g.ind.motbus
% g.ind.ind_int
% g.ind.t_mot
% g.ind.p_mot
% g.ind.q_mot
% g.ind.s_mot
% g.ind.vdmot
% g.ind.vqmot
% g.ind.idmot
% g.ind.iqmot
% g.ind.vdp
% g.ind.vqp
% g.ind.slip
% g.ind.dvdp
% g.ind.dvqp
% g.ind.dslip

% lfac // ac load flow
% g.lfac.bus_type
% g.lfac.gen_chg_idx
% g.lfac.PQV_no
% g.lfac.PQ_no
% g.lfac.volt_red
% g.lfac.ang_red

% lmod // active load modulation
% g.lmod.n_lmod
% g.lmod.lmod_idx
% g.lmod.lmod_con
% g.lmod.lmod_pot
% g.lmod.lmod_st
% g.lmod.dlmod_st
% g.lmod.lmod_sig
% g.lmod.lmod_data

% lmon // line monitoring
% g.lmon.lmon_con

% lsc // linear time-varying synchronizing controllers
% g.lsc.n_lsc
% g.lsc.lsc_idx
% g.lsc.lsc_con
% g.lsc.lsc_pot
% g.lsc.n_sensor
% g.lsc.sensor_set
% g.lsc.lsc_ess_idx
% g.lsc.lsc_sensor_idx
% g.lsc.theta_coi_pade
% g.lsc.theta_lsc_pade
% g.lsc.theta_sensor
% g.lsc.lsc1--15
% g.lsc.dlsc1--15
% g.lsc.lsc_scmd

% ncl // non-conforming loads
% g.ncl.n_load
% g.ncl.load_con
% g.ncl.load_pot

% mac // machines
% g.mac.mac_con
% g.mac.mac_pot
% g.mac.n_mac
% g.mac.n_em
% g.mac.n_tra
% g.mac.n_sub
% g.mac.n_ib
% g.mac.n_ivm
% g.mac.mac_int
% g.mac.mac_em_idx
% g.mac.mac_tra_idx
% g.mac.mac_sub_idx
% g.mac.mac_ib_idx
% g.mac.mac_ivm_idx
% g.mac.ibus_con
% g.mac.not_ib_idx
% g.mac.n_ib_em
% g.mac.n_ib_tra
% g.mac.n_ib_sub
% g.mac.n_ib_ivm
% g.mac.mac_ib_em
% g.mac.mac_ib_tra
% g.mac.mac_ib_sub
% g.mac.mac_ib_ivm
% g.mac.n_pm
% g.mac.pm_sig
% g.mac.pmech
% g.mac.telect
% g.mac.pelect
% g.mac.qelect
% g.mac.vex
% g.mac.eterm
% g.mac.ed
% g.mac.eq
% g.mac.psidpp
% g.mac.psiqpp
% g.mac.curd
% g.mac.curdg
% g.mac.curq
% g.mac.curqg
% g.mac.fldcur
% g.mac.psi_re
% g.mac.psi_im
% g.mac.cur_re
% g.mac.cur_im
% g.mac.mac_ang
% g.mac.mac_spd
% g.mac.edprime
% g.mac.eqprime
% g.mac.psikd
% g.mac.psikq
% g.mac.dmac_ang
% g.mac.dmac_spd
% g.mac.dedprime
% g.mac.deqprime
% g.mac.dpsikd
% g.mac.dpsikq
% g.mac.ivmmod_data
% g.mac.ivmmod_d_sig
% g.mac.ivmmod_e_sig

% pss // power system stabilizers
% g.pss.pss_con
% g.pss.pss_pot
% g.pss.n_pss
% g.pss.pss_idx
% g.pss.pss_sp_idx
% g.pss.pss_p_idx
% g.pss.pss_mb_idx
% g.pss.pss_exc_idx
% g.pss.pss_T
% g.pss.pss_T2
% g.pss.pss_T4
% g.pss.pss_T4_idx
% g.pss.pss_noT4_idx
% g.pss.pss1
% g.pss.pss2
% g.pss.pss3
% g.pss.dpss1
% g.pss.dpss2
% g.pss.dpss3

% pwr // pwrmod injections
% g.pwr.pwrmod_con
% g.pwr.n_pwrmod
% g.pwr.pwrmod_idx
% g.pwr.pwrmod_data
% g.pwr.pwrmod_p_sig
% g.pwr.pwrmod_p_st
% g.pwr.dpwrmod_p_st
% g.pwr.pwrmod_q_sig
% g.pwr.pwrmod_q_st
% g.pwr.dpwrmod_q_st

% reec // renewable energy electrical control
% g.reec.icmd
% g.reec.pref
% g.reec.qref
% g.reec.vblk
% g.reec.vdip
% g.reec.iqmax
% g.reec.iqmin
% g.reec.reec1
% g.reec.reec2
% g.reec.reec3
% g.reec.reec4
% g.reec.reec5
% g.reec.reec6
% g.reec.reec7
% g.reec.reec8
% g.reec.reec9
% g.reec.reec10
% g.reec.vref0
% g.reec.vref1
% g.reec.dreec1
% g.reec.dreec2
% g.reec.dreec3
% g.reec.dreec4
% g.reec.dreec5
% g.reec.dreec6
% g.reec.dreec7
% g.reec.dreec8
% g.reec.dreec9
% g.reec.dreec10
% g.reec.n_reec
% g.reec.pfaref
% g.reec.ess_idx
% g.reec.reec_con
% g.reec.reec_pot
% g.reec.reec_sig
% g.reec.vblk_tick
% g.reec.vblk_time
% g.reec.vdip_icmd
% g.reec.vdip_tick
% g.reec.vdip_time

% rlmod // reactive load modulation
% g.rlmod.n_rlmod
% g.rlmod.rlmod_idx
% g.rlmod.rlmod_con
% g.rlmod.rlmod_pot
% g.rlmod.rlmod_st
% g.rlmod.drlmod_st
% g.rlmod.rlmod_sig

% svc // static var compensators
% g.svc.svc_con
% g.svc.n_svc
% g.svc.svc_idx
% g.svc.svc_pot
% g.svc.svcll_idx
% g.svc.n_svcud
% g.svc.svcud_idx
% g.svc.svc_sig
% g.svc.svc_dsig
% g.svc.B_cv
% g.svc.dB_cv
% g.svc.B_con
% g.svc.dB_con

% sys // system variables
% g.sys.basmva
% g.sys.basrad
% g.sys.sys_freq
% g.sys.syn_ref
% g.sys.mach_ref
% g.sys.sw_con

% tcsc // thyristor-controlled series capacitors
% g.tcsc.tcsc_con
% g.tcsc.n_tcsc
% g.tcsc.tcscf_idx
% g.tcsc.tcsct_idx
% g.tcsc.B_tcsc
% g.tcsc.dB_tcsc
% g.tcsc.tcsc_sig
% g.tcsc.tcsc_dsig
% g.tcsc.n_tcscud
% g.tcsc.dtcscud_idx

% tg // turbine governors
% g.tg.tg_con
% g.tg.tg_pot
% g.tg.n_tg
% g.tg.n_tgh
% g.tg.n_tg_tot
% g.tg.tg_idx
% g.tg.tgh_idx
% g.tg.tg_sig
% g.tg.tg1
% g.tg.tg2
% g.tg.tg3
% g.tg.tg4
% g.tg.tg5
% g.tg.dtg1
% g.tg.dtg2
% g.tg.dtg3
% g.tg.dtg4
% g.tg.dtg5

% tripping status // variables for tracking the trip and reclosure status
% g.trip.mac_trip_flags
% g.trip.load_trip_flags
% g.trip.line_trip_flags

% eof
