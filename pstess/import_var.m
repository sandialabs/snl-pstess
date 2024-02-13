%
% Purpose: After the data file has been read, this script imports variables
%          from the workspace into the global struct g.

%-----------------------------------------------------------------------------%
% Version history
%
% Version:  1.0
% Purpose:  Initial version
% Date:     June 2021
% Author:   Ryan T. Elliott
%-----------------------------------------------------------------------------%

global g;  % declaring struct of global variables

% active load modulation
if exist('lmod_con','var')
    g.lmod.lmod_con = lmod_con;
    clear('lmod_con');
else
    g.lmod.lmod_con = [];
    g.lmod.n_lmod = 0;
end

% reactive load modulation
if exist('rlmod_con','var')
    g.rlmod.rlmod_con = rlmod_con;
    clear('rlmod_con');
else
    g.rlmod.rlmod_con = [];
    g.lmod.n_rlmod = 0;
end

% turbine governors
if exist('tg_con','var')
    g.tg.tg_con = tg_con;
    clear('tg_con');
else
    g.tg.tg_con = [];
end

% exciters
if exist('exc_con','var')
    g.exc.exc_con = exc_con;
    clear('exc_con');
else
    g.exc.exc_con = [];
end

% synchronous machines
if exist('mac_con','var')
    g.mac.mac_con = mac_con;
    clear('mac_con');

    g.mac.not_ivm_idx = 1:size(g.mac.mac_con,1);
    g.mac.n_not_ivm = length(g.mac.not_ivm_idx);
else
    g.mac.mac_con = [];
    g.mac.not_ivm_idx = [];
    g.mac.n_not_ivm = 0;
end

% infinite buses
if exist('ibus_con','var')
    g.mac.ibus_con = ibus_con;
    clear('ibus_con');
else
    g.mac.ibus_con = [];
    g.mac.n_ib = 0;
end

% power modulation
if exist('pwrmod_con','var')
    g.pwr.pwrmod_con = pwrmod_con;
    g.pwr.n_pwrmod = length(g.pwr.pwrmod_con(:,1));
    clear('pwrmod_con');
else
    g.pwr.pwrmod_con = [];
    g.pwr.n_pwrmod = 0;
end

% non-conforming loads
if exist('load_con','var')
    g.ncl.load_con = load_con;

    if ~isempty(g.ncl.load_con)
        g.ncl.n_load = size(g.ncl.load_con,1);
    else
        g.ncl.n_load = 0;
    end

    clear('load_con');
else
    g.ncl.load_con = [];
    g.ncl.n_load = 0;
end

% switching sequence
if exist('sw_con','var')
    g.sys.sw_con = sw_con;
    clear('sw_con');
else
    g.sys.sw_con = [];
end

% power system stabilizers
if exist('pss_con','var')
    g.pss.pss_con = pss_con;
    clear('pss_con');
else
    g.pss.pss_con = [];
end

% filters for dual-input stabilizers
if exist('dpw_con','var')
    g.dpw.dpw_con = dpw_con;
    clear('dpw_con');
else
    g.dpw.dpw_con = [];
end

% static var compensators
if exist('svc_con','var')
    g.svc.svc_con = svc_con;
    clear('svc_con');
else
    g.svc.svc_con = [];
end

% thyristor-controlled series capacitors
if exist('tcsc_con','var')
    g.tcsc.tcsc_con = tcsc_con;
    clear('tcsc_con');
else
    g.tcsc.tcsc_con = [];
    g.tcsc.n_tcsc = 0;
end

% induction generators
if exist('igen_con','var')
    g.igen.igen_con = igen_con;

    if ~isempty(g.igen.igen_con)
        g.igen.n_ig = size(g.igen.igen_con,1);
    else
        g.igen.n_ig = 0;
    end

    clear('igen_con');
else
    g.igen.igen_con = [];
    g.igen.n_ig = 0;
end

% induction motors
if exist('ind_con','var')
    g.ind.ind_con = ind_con;

    if ~isempty(g.ind.ind_con)
        g.ind.n_mot = size(g.ind.ind_con,1);
    else
        g.ind.n_mot = 0;
    end

    clear('ind_con');
else
    g.ind.ind_con = [];
    g.ind.n_mot = 0;
end

% induction motor loads
if exist('mld_con','var')
    g.ind.mld_con = mld_con;
    clear('mld_con');
else
    g.ind.mld_con = [];
end

% hvdc converters
if exist('dcsp_con','var')
    g.dc.dcsp_con = dcsp_con;
    clear('dcsp_con');
else
    g.dc.dcsp_con = [];
end

% hvdc converter controls
if exist('dcc_con','var')
    g.dc.dcc_con = dcc_con;
    clear('dcc_con');
else
    g.dc.dcc_con = [];
end

% hvdc lines
if exist('dcl_con','var')
    g.dc.dcl_con = dcl_con;
    clear('dcl_con');
else
    g.dc.dcl_con = [];
end

% line monitoring
if exist('lmon_con','var')
    g.lmon.lmon_con = lmon_con;
    clear('lmon_con');
else
    g.lmon.lmon_con = [];
end

% protection and RAS tripping
if exist('trip_enable','var')
    g.trip.enable = true;
    clear('trip_enable');
else
    g.trip.enable = false;
end

% synchronizing torque controllers
if exist('lsc_con','var')
    g.lsc.lsc_con = lsc_con;
    clear('lsc_con');
else
    g.lsc.lsc_con = [];
    g.lsc.n_lsc = 0;
end

% renewable energy electrical control models
if exist('reec_con','var')
    g.reec.reec_con = reec_con;
    clear('reec_con');
else
    g.reec.reec_con = [];
    g.reec.n_reec = 0;
end

% energy storage/converter interfaces
if exist('ess_con','var')
    g.ess.ess_con = ess_con;
    clear('ess_con');

    if exist('ess_vdl','var')
        g.ess.vdl = ess_vdl;
        clear('ess_vdl');
    end
else
    g.ess.ess_con = [];
    g.ess.n_ess = 0;
end

% user-defined ess controls
if exist('essud_con','var')
    g.ess.essud_con = essud_con;
    clear('essud_con');
else
    g.ess.essud_con = [];
    g.ess.n_essud = 0;
end

% internal voltage model generators
if exist('ivm_con','var')
    g.ivm.ivm_con = ivm_con;
    clear('ivm_con');
    ivm_indx(true);           % pre-indexing
else
    g.ivm.ivm_con = [];
    g.ivm.n_ivm = 0;
    g.mac.mac_ivm_idx = [];
    g.mac.n_ivm = 0;
end

% user-defined ivm controls
if exist('ivmud_con','var')
    g.ivm.ivmud_con = ivmud_con;
    clear('ivmud_con');
else
    g.ivm.ivmud_con = [];
    g.ivm.n_ivmud = 0;
end

% droop-based grid-forming inverter controls
if exist('gfma_con','var')
    g.gfma.gfma_con = gfma_con;
    clear('gfma_con');
else
    g.gfma.gfma_con = [];
    g.gfma.n_gfma = 0;
end

% eof
