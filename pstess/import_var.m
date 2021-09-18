%
% Purpose: After the data file has been read, this script imports variables
%          from the workspace into the global struct `g`.

%-----------------------------------------------------------------------------%
% Version history
%
% Version:  1.0
% Purpose:  Initial version
% Date:     June 2021
% Author:   Ryan T. Elliott
%-----------------------------------------------------------------------------%

global g;  % declaring struct of global variables

% lmod
if exist('lmod_con','var')
    g.lmod.lmod_con = lmod_con;
    clear('lmod_con');
else
    g.lmod.lmod_con = [];
    g.lmod.n_lmod = 0;
end

% rlmod
if exist('rlmod_con','var')
    g.rlmod.rlmod_con = rlmod_con;
    clear('rlmod_con');
else
    g.rlmod.rlmod_con = [];
    g.lmod.n_rlmod = 0;
end

% tg
if exist('tg_con','var')
    g.tg.tg_con = tg_con;
    clear('tg_con');
else
    g.tg.tg_con = [];
end

% exc
if exist('exc_con','var')
    g.exc.exc_con = exc_con;
    clear('exc_con');
else
    g.exc.exc_con = [];
end

% mac
if exist('mac_con','var')
    g.mac.mac_con = mac_con;
    clear('mac_con');
else
    g.mac.mac_con = [];
end

% ibus
if exist('ibus_con','var')
    g.mac.ibus_con = ibus_con;
    clear('ibus_con');
else
    g.mac.ibus_con = [];
end

% pwrmod
if exist('pwrmod_con','var')
    g.pwr.pwrmod_con = pwrmod_con;
    g.pwr.n_pwrmod = length(g.pwr.pwrmod_con(:,1));
    clear('pwrmod_con');
else
    g.pwr.pwrmod_con = [];
    g.pwr.n_pwrmod = 0;
end

% load_con
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

% sw_con
if exist('sw_con','var')
    g.sys.sw_con = sw_con;
    clear('sw_con');
else
    g.sys.sw_con = [];
end

% pss
if exist('pss_con','var')
    % account for washout gain alteration between version 2 and 3
    if exist('pss_scale_gain','var')
        if pss_scale_gain
            pss_con(:,3) = pss_con(:,3)*pss_con(:,4);
        end
    else
        pss_scale_gain = 0;
    end

    g.pss.pss_con = pss_con;
    g.pss.pss_scale_gain = pss_scale_gain;
    clear('pss_con','pss_scale_gain');
else
    g.pss.pss_con = [];
    g.pss.pss_scale_gain = nan;
end

% dpw filter
if exist('dpw_con','var')
    g.dpw.dpw_con = dpw_con;
    clear('dpw_con');
else
    g.dpw.dpw_con = [];
end

% svc_con
if exist('svc_con','var')
    g.svc.svc_con = svc_con;
    clear('svc_con');
else
    g.svc.svc_con = [];
end

% tcsc_con
if exist('tcsc_con','var')
    g.tcsc.tcsc_con = tcsc_con;
    clear('tcsc_con');
else
    g.tcsc.tcsc_con = [];
    g.tcsc.n_tcsc = 0;
end

% igen_con
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

% ind_con
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

% mld_con // for mac_ind
if exist('mld_con','var')
    g.ind.mld_con = mld_con;
    clear('mld_con');
else
    g.ind.mld_con = [];
end

% hvdc
if exist('dcsp_con','var')  % dc converter specifications
    g.dc.dcsp_con = dcsp_con;
    clear('dcsp_con');
else
    g.dc.dcsp_con = [];
end

if exist('dcl_con','var')   % dc lines
    g.dc.dcl_con = dcl_con;
    clear('dcl_con');
else
    g.dc.dcl_con = [];
end

if exist('dcc_con','var')   % dc converter controls
    g.dc.dcc_con = dcc_con;
    clear('dcc_con');
else
    g.dc.dcc_con = [];
end

% line monitoring
if exist('lmon_con','var')
    g.lmon.lmon_con = lmon_con;
    clear('lmon_con');
else
    g.lmon.lmon_con = [];
end

% lsc
if exist('lsc_con','var')
    g.lsc.lsc_con = lsc_con;
    clear('lsc_con');
else
    g.lsc.lsc_con = [];
    g.lsc.n_lsc = 0;
end

% ess
if exist('ess_con','var')
    g.ess.ess_con = ess_con;
    clear('ess_con');
else
    g.ess.ess_con = [];
    g.ess.n_ess = 0;
end

% ess_sud, user-defined ess controls
if exist('essud_con','var')
    g.ess.essud_con = essud_con;
    clear('essud_con');
else
    g.ess.essud_con = [];
    g.ess.n_essud = 0;
end

% eof
