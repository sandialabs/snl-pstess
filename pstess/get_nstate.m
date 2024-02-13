function [ nss ] = get_nstate( )
% Syntax: nss = get_nstate()
%
% Purpose: Returns the number of states associated with each
%          dynamic model implemented in PSTess. This information
%          is required in the linearization routine.
%
% Output:  nss - a struct containing information about the number of states
%
% Called by: svm_mgen

%-----------------------------------------------------------------------------%
% Version history
%
% Version: 1.0
% Author:  Ryan T. Elliott
% Date:    August 2019
% Note:    Initial version
%-----------------------------------------------------------------------------%

% states per model
ns_mac_ib = 0;                  % generator
ns_mac_em = 2;
ns_mac_tra = 4;
ns_mac_sub = 6;
ns_mac_ivm = 2;                 % internal voltage model
ns_mac_igen = 3;                % induction generator
ns_mac_ind = 3;                 % induction motor
ns_exc_smp = 3;                 % exciter
ns_exc_smppi = 3;
ns_exc_dc = 5;
ns_exc_st3 = 3;
ns_tg = 3;                      % governor
ns_tgh = 5;
ns_pss = 3;                     % pss
ns_dpw = 6;
ns_gfma = 6;                    % gfma
ns_svc = 2;                     % svc
ns_tcsc = 1;                    % tcsc
ns_lsc = 15;                    % lsc
ns_reec = 10;                   % reec
ns_ess = 5;                     % ess
ns_lmod = 1;                    % load modulation
ns_rlmod = 1;
ns_pwrmod = 1;                  % pwrmod
ns_dcl = 3;                     % hvdc link
ns_l_cap = 2;                   % hvdc capacitance

% output struct
nss.mac_ib = ns_mac_ib;
nss.mac_em = ns_mac_em;
nss.mac_tra = ns_mac_tra;
nss.mac_sub = ns_mac_sub;
nss.mac_ivm = ns_mac_ivm;
nss.ig = ns_mac_igen;
nss.mot = ns_mac_ind;
nss.exc_smp = ns_exc_smp;
nss.exc_smppi = ns_exc_smppi;
nss.exc_dc = ns_exc_dc;
nss.exc_st3 = ns_exc_st3;
nss.tg = ns_tg;
nss.tgh = ns_tgh;
nss.pss = ns_pss;
nss.dpw = ns_dpw;
nss.gfma = ns_gfma;
nss.svc = ns_svc;
nss.tcsc = ns_tcsc;
nss.lsc = ns_lsc;
nss.reec = ns_reec;
nss.ess = ns_ess;
nss.lmod = ns_lmod;
nss.rlmod = ns_rlmod;
nss.pwrmod = ns_pwrmod;
nss.dcl = ns_dcl;
nss.l_cap = ns_l_cap;

% maximum states per model type
nss.mac_max = max([ns_mac_ib, ns_mac_em, ns_mac_tra, ns_mac_sub, ns_mac_ivm]);
nss.mot_max = ns_mac_ind;
nss.ig_max = ns_mac_igen;
nss.tg_max = max([ns_tg, ns_tgh]);
nss.exc_max = max([ns_exc_smp, ns_exc_smppi, ns_exc_dc, ns_exc_st3]);
nss.pss_max = ns_pss;
nss.dpw_max = ns_dpw;
nss.gfma_max = ns_gfma;
nss.svc_max = ns_svc;
nss.tcsc_max = ns_tcsc;
nss.lsc_max = ns_lsc;
nss.reec_max = ns_reec;
nss.ess_max = ns_ess;
nss.lmod_max = ns_lmod;
nss.rlmod_max = ns_rlmod;
nss.pwrmod_max = ns_pwrmod;
nss.dcl_max = ns_dcl + ns_l_cap;

end  % function end

% eof
