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
g.lmod.lmod_con = [];
g.lmod.n_lmod = 0;

if exist('lmod_con','var')
    if ~isempty(lmod_con)
        g.lmod.lmod_con = lmod_con;
        g.lmod.n_lmod = size(lmod_con,1);
    end
    clear('lmod_con');
end

% reactive load modulation
g.rlmod.rlmod_con = [];
g.rlmod.n_rlmod = 0;

if exist('rlmod_con','var')
    if ~isempty(rlmod_con)
        g.rlmod.rlmod_con = rlmod_con;
        g.rlmod.n_rlmod = size(rlmod_con,1);
    end
    clear('rlmod_con');
end

% turbine governors
g.tg.tg_con = [];
g.tg.n_tg_tot = 0;

if exist('tg_con','var')
    if ~isempty(tg_con)
        g.tg.tg_con = tg_con;
        g.tg.n_tg_tot = size(tg_con,1);
    end
    clear('tg_con');
end

% exciters
g.exc.exc_con = [];
g.exc.n_exc = 0;

if exist('exc_con','var')
    if ~isempty(exc_con)
        g.exc.exc_con = exc_con;
        g.exc.n_exc = size(exc_con,1);
    end
    clear('exc_con');
end

% synchronous machines
g.mac.mac_con = [];
g.mac.n_mac = 0;
g.mac.not_ivm_idx = [];
g.mac.n_not_ivm = 0;

if exist('mac_con','var')
    if ~isempty(mac_con)
        g.mac.mac_con = mac_con;
        g.mac.n_mac = size(mac_con,1);
        g.mac.not_ivm_idx = 1:size(g.mac.mac_con,1);
        g.mac.n_not_ivm = length(g.mac.not_ivm_idx);
    end
    clear('mac_con');
end

% infinite buses
g.mac.ibus_con = [];
g.mac.n_ib = 0;

if exist('ibus_con','var')
    if ~isempty(ibus_con)
        g.mac.ibus_con = ibus_con;
        g.mac.n_ib = sum(ibus_con);
    end
    clear('ibus_con');
end

% power modulation
g.pwr.pwrmod_con = [];
g.pwr.n_pwrmod = 0;

if exist('pwrmod_con','var')
    if ~isempty(pwrmod_con)
        g.pwr.pwrmod_con = pwrmod_con;
        g.pwr.n_pwrmod = size(pwrmod_con,1);
    end
    clear('pwrmod_con');
end

% non-conforming loads
g.ncl.load_con = [];
g.ncl.n_load = 0;

if exist('load_con','var')
    if ~isempty(load_con)
        g.ncl.load_con = load_con;
        g.ncl.n_load = size(load_con,1);
    end
    clear('load_con');
end

% switching sequence
g.sys.sw_con = [];

if exist('sw_con','var')
    if ~isempty(sw_con)
        g.sys.sw_con = sw_con;
    end
    clear('sw_con');
end

% power system stabilizers
g.pss.pss_con = [];
g.pss.n_pss = 0;

if exist('pss_con','var')
    if ~isempty(pss_con)
        g.pss.pss_con = pss_con;
        g.pss.n_pss = size(pss_con,1);
    end
    clear('pss_con');
end

% filters for dual-input stabilizers
g.dpw.dpw_con = [];
g.dpw.n_dpw = 0;

if exist('dpw_con','var')
    if ~isempty(dpw_con)
        g.dpw.dpw_con = dpw_con;
        g.dpw.n_dpw = size(dpw_con,1);
    end
    clear('dpw_con');
end

% static var compensators
g.svc.svc_con = [];
g.svc.n_svc = 0;

if exist('svc_con','var')
    if ~isempty(svc_con)
        g.svc.svc_con = svc_con;
        g.svc.n_svc = size(svc_con,1);
    end
    clear('svc_con');
end

% thyristor-controlled series capacitors
g.tcsc.tcsc_con = [];
g.tcsc.n_tcsc = 0;

if exist('tcsc_con','var')
    if ~isempty(tcsc_con)
        g.tcsc.tcsc_con = tcsc_con;
        g.tcsc.n_tcsc = size(tcsc_con,1);
    end
    clear('tcsc_con');
end

% induction generators
g.igen.igen_con = [];
g.igen.n_ig = 0;

if exist('igen_con','var')
    if ~isempty(igen_con)
        g.igen.igen_con = igen_con;
        g.igen.n_ig = size(igen_con,1);
    end
    clear('igen_con');
end

% induction motors
g.ind.ind_con = [];
g.ind.n_mot = 0;

if exist('ind_con','var')
    if ~isempty(ind_con)
        g.ind.ind_con = ind_con;
        g.ind.n_mot = size(ind_con,1);
    end
    clear('ind_con');
end

% induction motor loads
g.ind.mld_con = [];

if exist('mld_con','var')
    if ~isempty(mld_con)
        g.ind.mld_con = mld_con;
    end
    clear('mld_con');
end

% hvdc converters
g.dc.dcsp_con = [];
g.dc.n_conv = 0;

if exist('dcsp_con','var')
    if ~isempty(dcsp_con)
        g.dc.dcsp_con = dcsp_con;
        g.dc.n_conv = size(dcsp_con,1);
    end
    clear('dcsp_con');
end

% hvdc converter controls
g.dc.dcc_con = [];

if exist('dcc_con','var')
    if ~isempty(dcc_con)
        g.dc.dcc_con = dcc_con;
    end
    clear('dcc_con');
end

% hvdc lines
g.dc.dcl_con = [];
g.dc.n_dcl = 0;

if exist('dcl_con','var')
    if ~isempty(dcl_con)
        g.dc.dcl_con = dcl_con;
        g.dc.n_dcl = size(dcl_con,1);
    end
    clear('dcl_con');
end

% line monitoring
g.lmon.lmon_con = [];

if exist('lmon_con','var')
    if ~isempty(lmon_con)
        g.lmon.lmon_con = lmon_con;
    end
    clear('lmon_con');
end

% protection and RAS tripping
g.trip.enable = false;

if exist('trip_enable','var')
    if ~isempty(trip_enable)
        g.trip.enable = logical(trip_enable);
    end
    clear('trip_enable');
end

% synchronizing torque controllers
g.lsc.lsc_con = [];
g.lsc.n_lsc = 0;

if exist('lsc_con','var')
    if ~isempty(lsc_con)
        g.lsc.lsc_con = lsc_con;
        g.lsc.n_lsc = size(lsc_con,1);
    end
    clear('lsc_con');
end

% renewable energy electrical control models
g.reec.reec_con = [];
g.reec.n_reec = 0;

if exist('reec_con','var')
    if ~isempty(reec_con)
        g.reec.reec_con = reec_con;
        g.reec.n_reec = size(reec_con,1);
    end
    clear('reec_con');
end

% energy storage/converter interfaces
g.ess.ess_con = [];
g.ess.n_ess = 0;

if exist('ess_con','var')
    if ~isempty(ess_con)
        g.ess.ess_con = ess_con;
        g.ess.n_ess = size(ess_con,1);
    end
    clear('ess_con');
end

if exist('ess_vdl','var')
    if (~isempty(ess_vdl) && ~isempty(g.ess.ess_con))
        g.ess.vdl = ess_vdl;
    end
    clear('ess_vdl');
end

% user-defined ess controls
g.ess.essud_con = [];
g.ess.n_essud = 0;

if exist('essud_con','var')
    if ~isempty(essud_con)
        g.ess.essud_con = essud_con;
        g.ess.n_essud = size(essud_con,1);
    end
    clear('essud_con');
end

% internal voltage model generators
g.ivm.ivm_con = [];
g.ivm.n_ivm = 0;
g.mac.mac_ivm_idx = [];
g.mac.n_ivm = 0;

if exist('ivm_con','var')
    if ~isempty(ivm_con)
        g.ivm.ivm_con = ivm_con;
        g.ivm.n_ivm = size(ivm_con,1);
        g.mac.n_ivm = g.ivm.n_ivm;
        ivm_indx(true);                 % pre-indexing
    end
    clear('ivm_con');
end

% user-defined ivm controls
g.ivm.ivmud_con = [];
g.ivm.n_ivmud = 0;

if exist('ivmud_con','var')
    if ~isempty(ivmud_con)
        g.ivm.ivmud_con = ivmud_con;
        g.ivm.n_ivmud = size(ivmud_con,1);
    end
    clear('ivmud_con');
end

% droop-based grid-forming inverter controls
g.gfma.gfma_con = [];
g.gfma.n_gfma = 0;

if exist('gfma_con','var')
    if ~isempty(gfma_con)
        g.gfma.gfma_con = gfma_con;
        g.gfma.n_gfma = size(gfma_con,1);
    end
    clear('gfma_con');
end

% eof
