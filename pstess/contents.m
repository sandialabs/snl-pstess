% Power and Energy Storage Systems Toolbox, Version 1.1
% February 2024
%
% Contents
%
% Folders
% /archive_pstV3p0 : archive of PST Version 3.0, on which PSTess is based
% /pwrmod_info     : information about pwrmod model developed by Montana Tech
%
% Utility functions
% cdps       : changes the current working directory
% get_dypar  : stores information about PSTess application settings
% get_nstate : stores data about the number of states associated with each model
% get_path   : stores information about the PSTess data path and working file
% set_path   : companion function that changes the info returned by get_path()
%
% Load flow functions
% calc     : calculates load flow mismatch
% chq_lim  : checks for generator var limits
% dc_lf    : performs hvdc load flow calculations
% form_jac : forms the load flow Jacobian
% inv_lf   : load flow computations for inverter
% lfdc     : ac and hvdc load flow driver (standalone)
% lfdcs    : solves load flow when dc lines are present in simulation
% lfdemo   : ac load flow driver (standalone)
% lftap    : modifies transformer tap settings
% loadflow : performs ac load flow calculations
% rec_lf   : load flow computations for the rectifiers
% vsdemo   : voltage stability driver
% y_sparse : forms the sparse admittance matrix of the network
%
% Dynamic models and helper functions
% dbcage      : finds the equiv. rotor resistance and reactance (for mac_ind)
% dc_cont     : dc converter
% dc_cur      : calculates dc line currents - used in nc_load
% dc_indx     : checks dc data for consistency and sets numbers and options
% dc_line     : hvdc line
% dc_load     : used in nc_load
% dc_sim      : dc only simulation routine (not well documented)
% dc_vidc     : updates Vdc and i_dc assuming constant ac voltages (for dc_sim)
% dci_sud     : hvdc inverter user-defined damping control
% dcr_sud     : hvdc rectifier user-defined damping control
% deepbar     : finds the equiv. rotor leakage reactance and resistance (for mac_ind)
% dessat      : vectorized describing funtion for saturation (for mac_ind)
% dpwf        : deltaP/omega filter for power system stabilizer units
% dpwf_indx   : index for deltaP/omega filters
% ess         : energy storage-based transient stability control model
% ess_indx    : index for energy storage systems
% ess_sat     : auxiliary function for saturation of inverter-based resources
% ess_sud     : user-defined energy storage damping control
% exc_dc12    : dc exciter and AVR
% exc_indx    : checks exciter data and presets exciter numbers and options
% exc_st3     : IEEE ST3 static exciter and AVR
% freqcalc    : dynamic model for simulating frequency measurements
% gfma        : droop-based grid-forming inverter control model
% gfma_indx   : index for droop-based grid-forming inverters
% i_simu      : forms the network interface variables (boundary current injections)
% import_var  : initializes all global model parameters required for simulation
% ind_ldto    : calculates induction motor load torque (for mac_ind)
% ivm_sud     : user-defined internal voltage control model
% ivm_indx    : index for internal voltage model generators
% line_cur    : calculates currents in lines from transient voltage records
% line_pq     : calculates real and reactive powers in lines
% lmod        : modulates specified real loads
% lmod_indx   : index for real load modulation
% lsc         : LTV synchronizing torque controller model (for use with ess)
% lsc_indx    : index for LTV synchronizing torque controllers
% mac_em      : two-state 'classical' generator model (swing dynamics only)
% mac_ib      : infinite bus generator model (no dynamics)
% mac_igen    : induction generator model
% mac_ind     : double-cage induction motor model (alternatively mac_ind_single_cage)
% mac_indx    : index for machines, also checks data and options
% mac_ivm     : voltage behind reactance converter model
% mac_sub     : subtransient generator
% mac_tra     : transient generator
% mdc_sig     : modulation function for hvdc
% mess_sig    : modulation function for ess reference input(s)
% mexc_sig    : modulation function for exciter Vref
% mgfma_sig   : modulation function for gfma reference input(s)
% mlmod_sig   : modulation function for active loads
% mpm_sig     : modulation function for generator mechanical power
% mreec_sig   : modulation function for reec reference input(s)
% mrlmod_sig  : modulation function for reactive loads
% msvc_sig    : modulation function for svc reference input
% mtcsc_sig   : modulation function for tcsc reference input
% mtg_sig     : modulation function for turbine governor reference input
% nc_load     : non-conforming loads, hvdc, svc, load modulation network interface
% ns_file     : determines total number of states in small signal stability
% p_cont      : controls perturbations for linear model development
% p_dpw       : perturbs deltaP/omega filter models
% p_exc       : perturbs exciter models
% p_file      : forms columns of state matrix, b, c, and d matrices
% p_m_file    : forms permutation matrix for state matrix calculation
% p_pss       : perturbs pss models
% p_tg        : perturbs turbine/governor models
% pss         : power system stabilizer
% pss_des     : power system stabilizer design
% pss_indx    : checks pss data and sets numbers and options
% pss_phse    : calculates phase shift through pss
% pwrmod_dyn  : specifies dynamics for pwrmod_p and pwrmod_q (also pwrmod_dyn_alt)
% pwrmod_indx : index for constant power/current injection
% pwrmod_p    : real power injection at spectified buses
% pwrmod_q    : reactive power injection at spectified buses
% rbus_ang    : computes bus angle changes
% red_ybus    : calculates reduced y matrices
% reec        : renewable energy electrical control model (for use with ess)
% reec_indx   : index for renewable energy electrical control
% rlmod       : modulates selected reactive loads
% rlmod_indx  : index for reactive load modulation
% rltf        : calculates root locus for transfer function feedback
% rootloc     : calculates rootlocus for scalar feedback around state space system
% s_simu      : transient simulation driver
% sd_torque   : calculates generator synchronizing and damping torques
% smpexc      : simple exciter model
% smppi       : simple excitation system with proportional-integral avr
% stab_d      : interactive pss design
% stab_f      : pss frequency response
% state_f     : frequency response from state space
% step_res    : step response from state space
% svc         : static var compensator
% svc_indx    : index for static var compensators (svc)
% svc_sud     : user-defined svc damping control
% svm_mgen    : small signal stability driver
% tcsc        : thryristor-controlled series capacitor model
% tcsc_indx   : index for thryristor-controlled series capacitors (tcsc)
% tcsc_sud    : user-defined tcsc damping control
% tg          : turbine/governor
% tg_hydro    : hydraulic turbine/governor
% tg_indx     : turbine/governor index
% y_switch    : organizes reduced y matrices
%
% Protection and RAS logic
% trip_handler : auxiliary function for managing protection and RAS actions
% trip_indx    : creates mappings necessary to implement load shedding
% trip_logic   : specifies protection and RAS logic (e.g., UFLS)
%
% Other
% insimit    : simultaneous iteration on inverse of a matrix (uncalled)
% nm_if      : functionalized network machine interface (uncalled)
% swcap      : specialized simulation switching routine (uncalled)
% time_stamp : generates a time stamp using the cpu clock (uncalled)
% ybus       : build admittance matrix Y from the line data (uncalled)

% eof
