% s_simu.m
%
% An m.file to simulate power system transients
% using the Matlab Power System Toolbox
% This m-file takes the dynamic and load flow data and
% calculates the response of the power system to a fault
% which is specified in a switching file
% see one of the supplied data files (data2a.m) for the
% switching file format

%-----------------------------------------------------------------------------%
% Version history
%
% Version: 1.92
% Author:  D. Trudnowski
% Date:    2023
% Purpose: Updated to allow for generator, load, and line tripping
%
% Version: 1.91
% Author:  Ryan Elliott
% Date:    2020
% Purpose: Updated to accommodate batch processing mode
%
% Version: 1.9
% Author:  Thad Haines
% Date:    June 2020
% Purpose: Modified from PSTv3 to accomodate batch style running
%          pwrmod, ivmmod, and gen trip sections from 2.3 copied
%
% Version: 1.86
% Author:  Ryan Elliott
% Date:    2019
% Purpose: Added code for energy storage/inverter-based resources
%
% Version: 1.85
% Author:  D. Trudnowski
% Date:    2019
% Purpose: Added ivmmod code
%
% Version: 1.8
% Author:  Rui
% Date:    July 2017
% Purpose: add code for multiple hvdc systems
%
% Version: 1.75
% Author:  D. Trudnowski
% Date:    2015
% Purpose: Added pwrmod code
%
% Version: 1.7
% Author:  Graham Rogers
% Date:    September 1999
% Purpose: Add simple exciter with pi avr - smppi
%
% Version: 1.6
% Author:  Graham Rogers
% Date:    June 1999
% Purpose: add user defined hvdc control, modify dc so that it is
%          integrated with a sub-multiple time constant
%
% Version: 1.5
% Author:  Graham Rogers
% Date:    December 1998/January 1999
% Purpose: Add svc and tcsc user defined damping controls, and add tcsc model
%
% Version: 1.4
% Author:  Graham Rogers
% Date:    July 1998
% Purpose: Add deltaP/omega filter
%
% Version: 1.3
% Author:  Graham Rogers
% Date:    June 1998
% Purpose: Add hydraulic turbine/governor
%
% Version: 1.2
% Author:  Graham Rogers
% Date:    August 1997
% Purpose: add induction generator
%
% Version: 1.1
% Author:  Graham Rogers
% Date:    August 1997
% Purpose: add modulation for load, exciter and governor
%
% Version  1.0 (initial version)
% Author:  Graham Rogers
% Date:    February 1997
%-----------------------------------------------------------------------------%

clear all; close all; clc;

%-----------------------------------------------------------------------------%

%% Starting s_simu
warning('*** Starting s_simu')

tic                                      % set timer
plot_now = 0;                            % plot flag

dypar = get_dypar();                     % get parameter settings
batch_mode = dypar.batch_mode;           % batch processing mode

svc_dc = {};                             % svc damping control
tcsc_dc = {};                            % tcsc damping control
dcr_dc = {};                             % hvdc rectifier damping control
dci_dc = {};                             % hvdc inverter damping control

% note: user-defined energy storage control data is entered in the
%       working case data file

% load input data from m.file
disp('non-linear simulation')
% input data file
if batch_mode
    [dfile,pathname] = get_path();
else
    [dfile,pathname] = uigetfile('d*.m','Select Data File');
end

% import base case data
if (pathname == 0)
    error('s_simu: you must select a valid data file.');
else
    full_dfile = fullfile(pathname, dfile);
    run(full_dfile);
    run('import_var');                   % set up global variables
end

% check for valid dynamic data file
if (isempty(g.mac.mac_con) && isempty(g.ivm.ivm_con))
    error('s_simu: the selected file is not a valid data file (mac_con).');
end

% check for valid switching sequence
if isempty(g.sys.sw_con)
    error('s_simu: the selected file has no switching data (sw_con).');
else
    if any(diff(g.sys.sw_con(:,1)) < 0)
         error('s_simu: time points in sw_con must be non-decreasing.');
    end
end

% specify the system mva base and nominal frequency
if batch_mode
    basdat = {num2str(dypar.basmva),num2str(dypar.sys_freq)}.';
else
    basdat = inputdlg({'Base MVA:','Base Frequency Hz:'}, ...
                      'Input Base Data',1,{'100','60'});
end

g.sys.sys_freq = str2double(basdat{2});  % nominal system freq, default 60 Hz
g.sys.basrad = 2*pi*g.sys.sys_freq;      % nominal system freq, default ~377 rad/s
g.sys.basmva = str2double(basdat{1});    % system complex power base, default 100 MVA
g.sys.syn_ref = 0;                       % angle reference machine index

% ignore infinite buses in time-domain simulation
g.mac.ibus_con = [];

%% Determine initial power flow solution
warning('*** Determining initial power flow solution')

if batch_mode
    lfs = 'y';
else
    lfs = inputdlg('Do you want to solve the power flow, (y/n)[y]: ','s');
    lfs = lower(lfs);
    if isempty(lfs)
        lfs = 'y';  % default
    end
end

if strcmp(lfs,'y')
    if isempty(g.dc.dcsp_con)
        % ac power flow
        g.dc.n_conv = 0;
        g.dc.n_dcl = 0;
        g.dc.ndcr_ud = 0;
        g.dc.ndci_ud = 0;
        g.dc.ldc_idx = [];
        %
        tol = 1e-8;                      % tolerance for convergence
        itermax = 50;                    % maximum number of iterations
        acc = 1.0;                       % acceleration factor

        [bus,line,line_flw] = loadflow(bus,line,tol,itermax,acc,'n',2);

        save('sim_fle','bus','line');
    else
        % has hvdc, use dc load flow
        [bus,line,line_flw,rec_par,inv_par,line_par] = lfdcs(bus,line,dci_dc,dcr_dc);

        save('sim_fle','bus','line','rec_par','inv_par','line_par');
    end
else
    load('sim_fle');
end

% adjusting bus matrix for ess instances specified as generation
if ~isempty(g.ess.ess_con)
    j_ess = g.bus.bus_int(g.ess.ess_con(:,2));

    mask = (bus(j_ess,10) ~= 3);
    if any(mask)
        bus(j_ess(mask),6) = -bus(j_ess(mask),4);
        bus(j_ess(mask),7) = -bus(j_ess(mask),5);
        bus(j_ess(mask),4) = 0.0;
        bus(j_ess(mask),5) = 0.0;
    end
end

% note: dc_indx() called in dc load flow (lfdcs)
mac_indx();
exc_indx();
dpwf_indx();
pss_indx();
tg_indx();
ivm_indx();
gfma_indx();
svc_indx(svc_dc);
tcsc_indx(tcsc_dc);
lsc_indx(bus);
reec_indx();
ess_indx();
lmod_indx();
rlmod_indx();
pwrmod_indx(bus);
trip_indx(bus);

% make sure bus max/min Q is the same as the pwrmod_con max/min Q
if ~isempty(g.pwr.n_pwrmod)
    for kk = 1:g.pwr.n_pwrmod
        n = find(g.pwr.pwrmod_con(kk,1) == bus(:,1));
        bus(n,11:12) = g.pwr.pwrmod_con(kk,6:7);
    end
    clear('kk','n');
end

n_tot = g.mac.n_mac + g.ind.n_mot + g.igen.n_ig;
n_gm = g.mac.n_mac + g.ind.n_mot;

% construct simulation switching sequence as defined in sw_con
tswitch(1) = g.sys.sw_con(1,1);
k = 1;
kdc = 1;
n_switch = length(g.sys.sw_con(:,1));
k_inc = zeros(n_switch-1,1);
k_incdc = k_inc;
t_switch = zeros(n_switch,1);
h = t_switch;
h_dc = h;

for sw_count = 1:n_switch-1
    h(sw_count) = g.sys.sw_con(sw_count,7);  % specified time step

    if (h(sw_count) == 0)
        h(sw_count) = 0.01;                  % default
    end

    % nearest lower integer
    k_inc(sw_count) = fix((g.sys.sw_con(sw_count+1,1) ...
                           - g.sys.sw_con(sw_count,1))/h(sw_count));

    if (k_inc(sw_count) == 0)
        k_inc(sw_count) = 1;                 % minimum 1
    end

    h(sw_count) = (g.sys.sw_con(sw_count+1,1) - g.sys.sw_con(sw_count,1)) ...
                  ./k_inc(sw_count);         % step length
    h_dc(sw_count) = h(sw_count)/10;
    k_incdc(sw_count) = 10*k_inc(sw_count);
    t_switch(sw_count+1) = t_switch(sw_count) + k_inc(sw_count)*h(sw_count);
    t(k:k-1+k_inc(sw_count)) = ...
        t_switch(sw_count):h(sw_count):t_switch(sw_count+1)-h(sw_count);
    t_dc(kdc:kdc-1+k_incdc(sw_count)) = ...
        t_switch(sw_count):h_dc(sw_count):t_switch(sw_count+1)-h_dc(sw_count);
    k = k + k_inc(sw_count);
    kdc = kdc + k_incdc(sw_count);
end

t_dc(kdc) = t_dc(kdc-1) + h_dc(sw_count);
for kk = 1:10
    kdc = kdc + 1;
    t_dc(kdc) = t_dc(kdc-1) + h_dc(sw_count);
end

k = sum(k_inc) + 1;                          % total number of time steps in the sim.
t(k) = g.sys.sw_con(n_switch,1);

[n,~] = size(g.mac.mac_con);
n_bus = length(bus(:,1));
g.bus.n_bus = n_bus;

% create zero matrices to initialize variables
z = zeros(n,k);
z1 = zeros(1,k);

zm = zeros(1,k);
if (g.ind.n_mot > 1)
    zm = zeros(g.ind.n_mot,k);
end

zig = zeros(1,k);
if (g.igen.n_ig > 1)
    zig = zeros(g.igen.n_ig,k);
end

zdc = zeros(2,kdc);
if (g.dc.n_conv > 2)
    zdc = zeros(g.dc.n_conv,kdc);
end

zdcl = zeros(1,kdc);
if (g.dc.n_dcl > 1)
    zdcl = zeros(g.dc.n_dcl,kdc);
end

% set dc parameters
g.dc.Vdc = zdc;
g.dc.i_dc = zdc;
g.dc.cur_ord = zdc;
g.dc.alpha = zdcl;
g.dc.gamma = zdcl;
g.dc.dc_sig = zeros(g.dc.n_conv,k);
g.dc.dcr_dsig = zeros(g.dc.n_dcl,k);
g.dc.dci_dsig = zeros(g.dc.n_dcl,k);
g.dc.i_dcr = zdcl;
g.dc.i_dci = zdcl;
g.dc.v_dcc = zdcl;
g.dc.v_conr = zdcl;
g.dc.v_coni = zdcl;
g.dc.di_dcr = zdcl;
g.dc.di_dci = zdcl;
g.dc.dv_dcc = zdcl;
g.dc.dv_conr = zdcl;
g.dc.dv_coni = zdcl;

if (g.dc.n_conv ~= 0)
    % modified by Rui on Oct. 5, 2016
    for ihvdc_count = 1:g.dc.n_dcl
        ire = g.dc.r_idx(ihvdc_count);
        g.dc.Vdc(ire,:) = rec_par(ihvdc_count,2);
        g.dc.i_dc(ire,:) = line_par(ihvdc_count);    % for pdci
        g.dc.i_dcr(ihvdc_count,:) = g.dc.i_dc(ire,:);
        g.dc.alpha(ihvdc_count,:) = rec_par(ihvdc_count,1)*pi/180;
    end

    for ihvdc_count = 1:g.dc.n_dcl
        iin = g.dc.i_idx(ihvdc_count);
        g.dc.Vdc(iin,:) = inv_par(ihvdc_count,2);
        g.dc.i_dc(iin,:) = line_par(ihvdc_count);
        g.dc.i_dci(ihvdc_count,:) = g.dc.i_dc(iin,:);
        g.dc.gamma(ihvdc_count,:) = inv_par(ihvdc_count,1)*pi/180;
    end
    % end modification by Rui

    g.dc.Vdc(g.dc.r_idx,:) = rec_par(:,2);
    g.dc.Vdc(g.dc.i_idx,:) = inv_par(:,2);

    g.dc.i_dc(g.dc.r_idx,:) = line_par;
    g.dc.i_dc(g.dc.i_idx,:) = line_par;

    g.dc.i_dcr(:,:) = g.dc.i_dc(g.dc.r_idx,:);
    g.dc.i_dci(:,:) = g.dc.i_dc(g.dc.i_idx,:);

    g.dc.alpha(:,:) = rec_par(:,1)*pi/180;
    g.dc.gamma(:,:) = inv_par(:,1)*pi/180;

    if (g.dc.ndcr_ud ~= 0)
        for j = 1:g.dc.ndcr_ud
            sv = get(dcr_dc{j,1});
            if (j == 1)
                xdcr_dc = zeros(sv.NumStates,kdc);
                dxdcr_dc = zeros(sv.NumStates,kdc);
            else
                xdcr_dc = [xdcr_dc;zeros(sv.NumStates,kdc)];
                dxdcr_dc = [dxdcr_dc;zeros(sv.NumStates,kdc)];
            end
        end
        dcrd_sig = zeros(g.dc.ndcr_ud,k);
        angdcr = zeros(g.dc.ndcr_ud,k);
    else
        xdcr_dc = zeros(1,kdc);
        dxdcr_dc = zeros(1,kdc);
        dcrd_sig = zeros(1,k);
    end
    %
    if (g.dc.ndci_ud ~= 0)
        for j = 1:g.dc.ndci_ud
            sv = get(dci_dc{j,1});
            if (j == 1)
                xdci_dc = zeros(sv.NumStates,kdc);
                dxdci_dc = zeros(sv.NumStates,kdc);
            else
                xdci_dc = [xsvc_dc;zeros(sv.NumStates,kdc)];
                dxdci_dc = [dxsvc_dc;zeros(sv.NumStates,kdc)];
            end
        end
        dcid_sig = zeros(g.dc.ndci_ud,k);
        angdci = zeros(g.dc.ndci_ud,k);
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
g.bus.theta = zeros(n_bus+1,k);
g.bus.bus_v = zeros(n_bus+1,k);

g.mac.pmech = z;
g.mac.telect = z;
g.mac.pelect = z;
g.mac.qelect = z;
g.mac.pm_sig = z;
g.mac.ed = z;
g.mac.eq = z;
g.mac.curd = z;
g.mac.curq = z;
g.mac.curdg = z;
g.mac.curqg = z;
g.mac.fldcur = z;
g.mac.vex = z;
g.mac.eterm = z;
g.mac.cur_re = z;
g.mac.cur_im = z;
g.mac.psi_re = z;
g.mac.psi_im = z;
g.mac.mac_ang = z;   % states
g.mac.mac_spd = z;
g.mac.edprime = z;
g.mac.eqprime = z;
g.mac.psikd = z;
g.mac.psikq = z;
g.mac.dmac_ang = z;  % state derivatives
g.mac.dmac_spd = z;
g.mac.dedprime = z;
g.mac.deqprime = z;
g.mac.dpsikd = z;
g.mac.dpsikq = z;

z_tg = zeros(1,k);
if (g.tg.n_tg_tot ~= 0)
    z_tg = zeros(g.tg.n_tg_tot,k);
end

g.tg.tg1 = z_tg;
g.tg.tg2 = z_tg;
g.tg.tg3 = z_tg;
g.tg.tg4 = z_tg;
g.tg.tg5 = z_tg;
g.tg.dtg1 = z_tg;
g.tg.dtg2 = z_tg;
g.tg.dtg3 = z_tg;
g.tg.dtg4 = z_tg;
g.tg.dtg5 = z_tg;
g.tg.tg_sig = z_tg;

z_pss = zeros(1,k);
if (g.pss.n_pss ~= 0)
    z_pss = zeros(g.pss.n_pss,k);
end

g.pss.pss1 = z_pss;
g.pss.pss2 = z_pss;
g.pss.pss3 = z_pss;
g.pss.dpss1 = z_pss;
g.pss.dpss2 = z_pss;
g.pss.dpss3 = z_pss;

z_dpw = zeros(1,k);
if (g.dpw.n_dpw ~= 0)
    z_dpw = zeros(g.dpw.n_dpw,k);
end

g.dpw.sdpw1 = z_dpw;
g.dpw.sdpw2 = z_dpw;
g.dpw.sdpw3 = z_dpw;
g.dpw.sdpw4 = z_dpw;
g.dpw.sdpw5 = z_dpw;
g.dpw.sdpw6 = z_dpw;
g.dpw.dsdpw1 = z_dpw;
g.dpw.dsdpw2 = z_dpw;
g.dpw.dsdpw3 = z_dpw;
g.dpw.dsdpw4 = z_dpw;
g.dpw.dsdpw5 = z_dpw;
g.dpw.dsdpw6 = z_dpw;
g.dpw.dpw_out = z_dpw;

ze = zeros(1,k);
if (g.exc.n_exc ~= 0)
    ze = zeros(g.exc.n_exc,k);
end

g.exc.V_A = ze;
g.exc.V_B = ze;
g.exc.V_TR = ze;
g.exc.V_R = ze;
g.exc.V_As = ze;
g.exc.Efd = ze;
g.exc.R_f = ze;
g.exc.dV_TR = ze;
g.exc.dV_R = ze;
g.exc.dV_As = ze;
g.exc.dEfd = ze;
g.exc.dR_f = ze;
g.exc.exc_sig = ze;
g.exc.pss_out = ze;

% ivm user-defined control variable declarations
if ((g.ivm.n_ivm ~= 0) && (g.ivm.n_ivmud ~= 0))
    xivm_dc.s{1} = zeros(g.ivm.n_ivmud,k);
end

% gfma variable declarations
if (g.gfma.n_gfma ~= 0)
    z_gfma = zeros(g.gfma.n_gfma,k);
    o_gfma = ones(g.gfma.n_gfma,k);

    g.gfma.gfma1 = z_gfma;
    g.gfma.gfma2 = z_gfma;
    g.gfma.gfma3 = z_gfma;
    g.gfma.gfma4 = z_gfma;
    g.gfma.gfma5 = z_gfma;
    g.gfma.gfma6 = z_gfma;
    g.gfma.gfma7 = z_gfma;
    g.gfma.gfma8 = z_gfma;
    g.gfma.gfma9 = z_gfma;
    g.gfma.gfma10 = z_gfma;

    g.gfma.dgfma1 = z_gfma;
    g.gfma.dgfma2 = z_gfma;
    g.gfma.dgfma3 = z_gfma;
    g.gfma.dgfma4 = z_gfma;
    g.gfma.dgfma5 = z_gfma;
    g.gfma.dgfma6 = z_gfma;
    g.gfma.dgfma7 = z_gfma;
    g.gfma.dgfma8 = z_gfma;
    g.gfma.dgfma9 = z_gfma;
    g.gfma.dgfma10 = z_gfma;

    g.gfma.gfma_sig = z_gfma;

    g.gfma.pset = o_gfma;
    g.gfma.qset = o_gfma;
    g.gfma.vset = o_gfma;
    g.gfma.lim_flag = false(g.gfma.n_gfma,1);
end

g.ind.t_mot = zm;
g.ind.p_mot = zm;
g.ind.q_mot = zm;
g.ind.s_mot = zm;
g.ind.vdp = zm;
g.ind.vqp = zm;
g.ind.slip = zm;
g.ind.dvdp = zm;
g.ind.dvqp = zm;
g.ind.dslip = zm;

g.igen.tmig = zig;
g.igen.pig = zig;
g.igen.qig = zig;
g.igen.s_igen = zig;
g.igen.vdpig = zig;
g.igen.vqpig = zig;
g.igen.slig = zig;
g.igen.dvdpig = zig;
g.igen.dvqpig = zig;
g.igen.dslig = zig;

if (g.svc.n_svc ~= 0)
    g.svc.B_cv = zeros(g.svc.n_svc,k);
    g.svc.dB_cv = zeros(g.svc.n_svc,k);
    g.svc.svc_sig = zeros(g.svc.n_svc,k);
    g.svc.svc_dsig = zeros(g.svc.n_svc,k);
    g.svc.B_con = zeros(g.svc.n_svc,k);
    g.svc.dB_con = zeros(g.svc.n_svc,k);
    %
    if (g.svc.n_svcud ~= 0)
        d_sig = zeros(g.svc.n_svcud,k);
        for j = 1:g.svc.n_svcud
            sv = get(svc_dc{j,1});
            if (j == 1)
                xsvc_dc = zeros(sv.NumStates,k);
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
    g.svc.B_cv = zeros(1,k);
    g.svc.dB_cv = zeros(1,k);
    g.svc.svc_sig = zeros(1,k);
    g.svc.svc_dsig = zeros(1,k);
    g.svc.B_con = zeros(1,k);
    g.svc.dB_con = zeros(1,k);
    %
    xsvc_dc = zeros(1,k);
    dxsvc_dc = zeros(1,k);
    d_sig = zeros(1,k);
end

if (g.tcsc.n_tcsc ~= 0)
    g.tcsc.B_tcsc = zeros(g.tcsc.n_tcsc,k);
    g.tcsc.dB_tcsc = zeros(g.tcsc.n_tcsc,k);
    g.tcsc.tcsc_sig = zeros(g.tcsc.n_tcsc,k);
    g.tcsc.tcsc_dsig = zeros(g.tcsc.n_tcsc,k);
    %
    if (g.tcsc.n_tcscud ~= 0)
        td_sig = zeros(g.tcsc.n_tcscud,k);          % input to tcsc damping control
        for j = 1:g.tcsc.n_tcscud
            sv = get(tcsc_dc{j,1});                 % damp. ctrl. state space object
            if (j == 1)
                xtcsc_dc = zeros(sv.NumStates,k);   % tcsc damping control states
                dxtcsc_dc = zeros(sv.NumStates,k);  % tcsc dc state derivatives
            else
                % in order of damping controls
                xtcsc_dc = [xtcsc_dc;zeros(sv.NumStates,k)];
                dxtcsc_dc = [dxtcsc_dc;zeros(sv.NumStates,k)];
            end
        end
    else
        xtcsc_dc = zeros(1,k);
        dxtcsc_dc = zeros(1,k);
    end
else
    g.tcsc.B_tcsc = zeros(1,k);
    g.tcsc.dB_tcsc = zeros(1,k);
    g.tcsc.tcsc_sig = zeros(1,k);
    g.tcsc.tcsc_dsig = zeros(1,k);
    %
    xtcsc_dc = zeros(1,k);
    dxtcsc_dc = zeros(1,k);
    td_sig = zeros(1,k);
end

% lsc variable declarations
if (g.lsc.n_lsc ~= 0)
    z_lsc = zeros(g.lsc.n_lsc,k);

    g.lsc.lsc1 = z_lsc;
    g.lsc.lsc2 = z_lsc;
    g.lsc.lsc3 = z_lsc;
    g.lsc.lsc4 = z_lsc;
    g.lsc.lsc5 = z_lsc;
    g.lsc.lsc6 = z_lsc;
    g.lsc.lsc7 = z_lsc;
    g.lsc.lsc8 = z_lsc;
    g.lsc.lsc9 = z_lsc;
    g.lsc.lsc10 = z_lsc;
    g.lsc.lsc11 = z_lsc;
    g.lsc.lsc12 = z_lsc;
    g.lsc.lsc13 = z_lsc;
    g.lsc.lsc14 = z_lsc;
    g.lsc.lsc15 = z_lsc;

    g.lsc.dlsc1 = z_lsc;
    g.lsc.dlsc2 = z_lsc;
    g.lsc.dlsc3 = z_lsc;
    g.lsc.dlsc4 = z_lsc;
    g.lsc.dlsc5 = z_lsc;
    g.lsc.dlsc6 = z_lsc;
    g.lsc.dlsc7 = z_lsc;
    g.lsc.dlsc8 = z_lsc;
    g.lsc.dlsc9 = z_lsc;
    g.lsc.dlsc10 = z_lsc;
    g.lsc.dlsc11 = z_lsc;
    g.lsc.dlsc12 = z_lsc;
    g.lsc.dlsc13 = z_lsc;
    g.lsc.dlsc14 = z_lsc;
    g.lsc.dlsc15 = z_lsc;

    g.lsc.theta_sensor = zeros(max(g.lsc.n_sensor,1),k);
    g.lsc.theta_coi_pade = z_lsc;
    g.lsc.theta_lsc_pade = z_lsc;
    % g.lsc.theta_coi_del = z_lsc;  % for shift-register delay model

    g.lsc.lsc_scmd = z_lsc;
end

% reec variable declarations
if (g.reec.n_reec ~= 0)
    % reec variable declarations
    z_reec = zeros(g.reec.n_reec,k);
    o_reec = ones(g.reec.n_reec,k);

    g.reec.reec1 = z_reec;                        % states
    g.reec.reec2 = z_reec;
    g.reec.reec3 = z_reec;
    g.reec.reec4 = z_reec;
    g.reec.reec5 = z_reec;
    g.reec.reec6 = z_reec;
    g.reec.reec7 = z_reec;
    g.reec.reec8 = z_reec;
    g.reec.reec9 = z_reec;
    g.reec.reec10 = z_reec;

    g.reec.dreec1 = z_reec;                       % derivatives
    g.reec.dreec2 = z_reec;
    g.reec.dreec3 = z_reec;
    g.reec.dreec4 = z_reec;
    g.reec.dreec5 = z_reec;
    g.reec.dreec6 = z_reec;
    g.reec.dreec7 = z_reec;
    g.reec.dreec8 = z_reec;
    g.reec.dreec9 = z_reec;
    g.reec.dreec10 = z_reec;

    g.reec.icmd = z_reec;                         % commands
    g.reec.paux = z_reec;
    g.reec.reec_sig = z_reec;                     % modulation

    g.reec.iqmin = ones(g.reec.n_reec,1);
    g.reec.iqmax = ones(g.reec.n_reec,1);

    g.reec.pref = o_reec;                         % references
    g.reec.qref = o_reec;
    g.reec.pfaref = zeros(g.reec.n_reec,1);

    g.reec.vref0 = ones(g.reec.n_reec,1);
    g.reec.vref1 = ones(g.reec.n_reec,1);

    g.reec.vdip = false(g.reec.n_reec,1);
    g.reec.vdip_tick = -1*ones(g.reec.n_reec,1);
    g.reec.vdip_time = zeros(g.reec.n_reec,1);
    g.reec.vdip_icmd = zeros(g.reec.n_reec,1);

    g.reec.vblk = false(g.reec.n_reec,1);
    g.reec.vblk_tick = -1*ones(g.reec.n_reec,1);
    g.reec.vblk_time = zeros(g.reec.n_reec,1);
end

% ess variable declarations
if (g.ess.n_ess ~= 0)
    z_ess = zeros(g.ess.n_ess,k);

    g.ess.ess1 = z_ess;
    g.ess.ess2 = z_ess;
    g.ess.ess3 = z_ess;
    g.ess.ess4 = z_ess;
    g.ess.ess5 = z_ess;

    g.ess.dess1 = z_ess;
    g.ess.dess2 = z_ess;
    g.ess.dess3 = z_ess;
    g.ess.dess4 = z_ess;
    g.ess.dess5 = z_ess;

    g.ess.ess_cur = z_ess;
    g.ess.ess_soc = z_ess;
    g.ess.ess_vmag = z_ess;
    g.ess.ess_scmd = z_ess;
    g.ess.ess_sinj = z_ess;
    g.ess.ess_iord = z_ess;
    g.ess.ess_vmag_pade = z_ess;

    g.ess.ess_sig = z_ess;
    g.ess.ess_dsig = z_ess;
end

% ess user-defined control variable declarations
if ((g.ess.n_ess ~= 0) && (g.ess.n_essud ~= 0))
    z_essud = zeros(g.ess.n_essud,k);
    xess_dc.s{1} = z_essud;  % dynamic states
else
    xess_dc.s{1} = [];       % empty array
end

% initialize lmod and rlmod
if (g.lmod.n_lmod ~= 0)
    g.lmod.lmod_st = zeros(g.lmod.n_lmod,k);
    g.lmod.dlmod_st = g.lmod.lmod_st;
    g.lmod.lmod_sig = g.lmod.lmod_st;
else
    g.lmod.lmod_st = zeros(1,k);
    g.lmod.dlmod_st = g.lmod.lmod_st;
    g.lmod.lmod_sig = g.lmod.lmod_st;
end

if (g.rlmod.n_rlmod ~= 0)
    g.rlmod.rlmod_st = zeros(g.rlmod.n_rlmod,k);
    g.rlmod.drlmod_st = g.rlmod.rlmod_st;
    g.rlmod.rlmod_sig = g.rlmod.rlmod_st;
else
    g.rlmod.rlmod_st = zeros(1,k);
    g.rlmod.drlmod_st = g.rlmod.rlmod_st;
    g.rlmod.rlmod_sig = g.rlmod.rlmod_st;
end

% initialize pwrmod
if (g.pwr.n_pwrmod ~= 0)
    g.pwr.pwrmod_p_st = zeros(g.pwr.n_pwrmod,k);
    g.pwr.dpwrmod_p_st = g.pwr.pwrmod_p_st;
    g.pwr.pwrmod_p_sig = g.pwr.pwrmod_p_st;

    g.pwr.pwrmod_q_st = zeros(g.pwr.n_pwrmod,k);
    g.pwr.dpwrmod_q_st = g.pwr.pwrmod_q_st;
    g.pwr.pwrmod_q_sig = g.pwr.pwrmod_q_st;
else
    g.pwr.pwrmod_p_st = zeros(1,k);
    g.pwr.dpwrmod_p_st = g.pwr.pwrmod_p_st;
    g.pwr.pwrmod_p_sig = g.pwr.pwrmod_p_st;

    g.pwr.pwrmod_q_st = zeros(1,k);
    g.pwr.dpwrmod_q_st = g.pwr.pwrmod_q_st;
    g.pwr.pwrmod_q_sig = g.pwr.pwrmod_q_st;
end

% initialize trip flags
if g.trip.enable
    g.trip.mac_trip_flags = false(g.mac.n_mac,k);
    g.trip.load_trip_flags = false(n_bus,k);
    g.trip.load_trip_frac = zeros(n_bus,k);
    g.trip.line_trip_flags = false(size(line,1),k);
    g.trip.load_trip_ncl = zeros(g.trip.n_trip_ncl,k);
    g.trip.user_variables = [];
end

g.sys.sys_freq = ones(1,k);

%% Initializing Y matrices and dynamic models
warning('*** Initializing Y matrices and dynamic models')

disp('constructing reduced y matrices')
disp('initializing motor and induction generator models')
bus = mac_ind(0,1,bus,0);   % initialize induction motor models
bus = mac_igen(0,1,bus,0);  % initialize induction generator models

% this has to be done before red_ybus is called since the motor and svc
% initialization alters the bus matrix and dc parameters are required
disp('initializing svc, tcsc, lsc, ess, and hvdc models')
bus = svc(0,1,bus,0);       % initialize svc models
tcsc(0,1,0);                % initialize tcsc models
lsc(0,1,bus,0);             % initialize lsc models
reec(0,1,bus,0);            % initialize reec models
ess(0,1,bus,0);             % initialize ess models
dc_cont(0,1,1,bus,0);       % initialize hvdc control models

% initialize user-defined energy storage control models
xess_dc = ess_sud(0,1,bus,0,xess_dc);

% calculate the reduced y matrices for the different switching conditions
run('y_switch');

% step 2: initialization
disp('initializing other models')

g.bus.theta(1:n_bus,1) = bus(:,3)*pi/180;
g.bus.bus_v(1:n_bus,1) = bus(:,2).*exp(1j*g.bus.theta(1:n_bus,1));
freqcalc(k,t,0);            % Felipe's frequency measurement model

if (g.svc.n_svcud ~= 0)
    % calculate the initial magnitude of line current for svc damping controls
    for j = 1:g.svc.n_svcud
        l_num = svc_dc{j,3};
        svc_num = svc_dc{j,2};
        from_bus = g.bus.bus_int(line(l_num,1));
        to_bus = g.bus.bus_int(line(l_num,2));

        svc_bn = g.bus.bus_int(g.svc.svc_con(svc_num,2));
        if (svc_bn ~= from_bus && svc_bn ~= to_bus)
            error('s_simu: the svc is not at the end of the specified line.');
        end

        V1 = g.bus.bus_v(from_bus,1);
        V2 = g.bus.bus_v(to_bus,1);
        R = line(l_num,3);
        X = line(l_num,4);
        B = line(l_num,5);
        tap = line(l_num,6);
        phi = line(l_num,7);

        [l_if,l_it] = line_cur(V1,V2,R,X,B,tap,phi);
        l_if0(j) = l_if;
        l_it0(j) = l_it;

        if (svc_bn == from_bus)
            d_sig(j,1) = abs(l_if);
        elseif (svc_bn == to_bus)
            d_sig(j,1) = abs(l_it);
        end
    end
end

if (g.tcsc.n_tcscud ~= 0)
    % calculate the initial bus voltage magnitude for tcsc damping controls
    for j = 1:g.tcsc.n_tcscud
        b_num = tcsc_dc{j,3};
        tcsc_num = tcsc_dc{j,2};
        td_sig(j,1) = abs(g.bus.bus_v(g.bus.bus_int(b_num),1));
    end
end

if (g.dc.n_conv ~= 0)
    % change dc buses from LT to HT
    Pr = bus(g.dc.rec_ac_bus,6);
    Qr = bus(g.dc.rec_ac_bus,7);

    Pi = bus(g.dc.inv_ac_bus,6);
    Qi = bus(g.dc.inv_ac_bus,7);

    VLT = g.bus.bus_v(g.dc.ac_bus,1);
    i_acr = (Pr-1j*Qr)./conj(VLT(g.dc.r_idx));
    i_aci = (Pi-1j*Qi)./conj(VLT(g.dc.i_idx));

    IHT(g.dc.r_idx,1) = i_acr;
    IHT(g.dc.i_idx,1) = i_aci;

    VHT(g.dc.r_idx,1) = VLT(g.dc.r_idx) + 1j*g.dc.dcc_pot(:,2).*i_acr;
    VHT(g.dc.i_idx,1) = VLT(g.dc.i_idx) + 1j*g.dc.dcc_pot(:,4).*i_aci;

    g.bus.bus_v(g.dc.ac_bus,1) = VHT;
    g.bus.theta(g.dc.ac_bus,1) = angle(g.bus.bus_v(g.dc.ac_bus,1));

    % modify the bus matrix to the HT buses
    bus(g.dc.ac_bus,2) = abs(g.bus.bus_v(g.dc.ac_bus,1));
    bus(g.dc.ac_bus,3) = g.bus.theta(g.dc.ac_bus,1)*180/pi;

    SHT = VHT.*conj(IHT);
    bus(g.dc.ac_bus,6) = real(SHT);
    bus(g.dc.ac_bus,7) = imag(SHT);

    if (g.dc.ndcr_ud ~= 0)
        % calculate the initial value of bus angles rectifier user defined control
        for j = 1:g.dc.ndcr_ud
            b_num1 = dcr_dc{j,3};
            b_num2 = dcr_dc{j,4};
            conv_num = dcr_dc{j,2};
            angdcr(j,:) = g.bus.theta(g.bus.bus_int(b_num1),1) ...
                          - g.bus.theta(g.bus.bus_int(b_num2),1);
            dcrd_sig(j,:) = angdcr(j,:);
        end
    end
    %
    if (g.dc.ndci_ud ~= 0)
        % calculate the initial value of bus angles inverter user defined control
        for j = 1:g.dc.ndci_ud
            b_num1 = dci_dc{j,3};
            b_num2 = dci_dc{j,4};
            conv_num = dci_dc{j,2};
            angdci(j,:) = g.bus.theta(g.bus.bus_int(b_num1),1) ...
                          - g.bus.theta(g.bus.bus_int(b_num2),1);
            dcid_sig(j,:) = angdci(j,:);
        end
    end
end

%% Initializing dynamic models
warning('*** Initializing dynamic models')

flag = 0;
g.bus.bus_int = bus_intprf;     % pre-fault system

disp('generators')
mac_sub(0,1,bus,flag);
mac_tra(0,1,bus,flag);
mac_em(0,1,bus,flag);
mac_ivm(0,1,bus,flag);

disp('generator controls')
dpwf(0,1,flag);
pss(0,1,flag);
smpexc(0,1,flag);
smppi(0,1,flag);
exc_st3(0,1,flag);
exc_dc12(0,1,flag);
tg(0,1,flag);
tg_hydro(0,1,flag);

disp('ivm generator controls')
gfma(0,1,bus,flag);

% user-defined ivm control models
if ((g.ivm.n_ivm ~= 0) && (g.ivm.n_ivmud ~= 0))
    xivm_dc = ivm_sud(0,1,bus,flag,xivm_dc);
end

% initialize svc damping controls
if (g.svc.n_svcud ~= 0)
    tot_states = 0;
    for i = 1:g.svc.n_svcud
        ysvcmx = svc_dc{i,4};
        ysvcmn = svc_dc{i,5};
        svc_num = svc_dc{i,2};
        st_state = tot_states + 1;
        svc_states = svc_dc{i,6};
        tot_states = tot_states + svc_states;
        [g.svc.svc_dsig(svc_num,1),xsvc_dc(st_state:tot_states,1),dxsvc_dc(st_state:tot_states,1)] = ...
            svc_sud(i,1,flag,svc_dc{i,1},d_sig(i,1),ysvcmx,ysvcmn);
    end
end

% initialize tcsc damping controls
if (g.tcsc.n_tcscud ~= 0)
    tot_states = 0;
    for i = 1:g.tcsc.n_tcscud
        ytcscmx = tcsc_dc{i,4};
        ytcscmn = tcsc_dc{i,5};
        tcsc_num = tcsc_dc{i,2};
        st_state = tot_states + 1;
        tcsc_states = tcsc_dc{i,6};
        tot_states = tot_states + tcsc_states;
        [g.tcsc.tcsc_dsig(tcsc_num,1),xtcsc_dc(st_state:tot_states,1),dxtcsc_dc(st_state:tot_states,1)] = ...
            tcsc_sud(i,1,flag,tcsc_dc{i,1},td_sig(i,1),ytcscmx,ytcscmn);
    end
end

% initialize hvdc rectifier damping controls
if (g.dc.ndcr_ud ~= 0)
    tot_states = 0;
    for i = 1:g.dc.ndcr_ud
        ydcrmx = dcr_dc{i,5};
        ydcrmn = dcr_dc{i,6};
        rec_num = dcr_dc{i,2};
        st_state = tot_states + 1;
        dcr_states = dcr_dc{i,7};
        tot_states = tot_states + dcr_states;
        [g.dc.dcr_dsig(rec_num,1),xdcr_dc(st_state:tot_states,1),dxdcr_dc(st_state:tot_states,1)] = ...
            dcr_sud(i,1,flag,dcr_dc{i,1},dcrd_sig(i,1),ydcrmx,ydcrmn);
    end
end

% initialize hvdc inverter damping controls
if (g.dc.ndci_ud ~= 0)
    tot_states = 0;
    for i = 1:g.dc.ndci_ud
        ydcimx = dci_dc{i,5};
        ydcrmn = dci_dc{i,6};
        inv_num = dci_dc{i,2};
        st_state = tot_states + 1;
        dci_states = dci_dc{i,7};
        tot_states = tot_states + dci_states;
        [g.dc.dci_dsig(inv_num,1),xdci_dc(st_state:tot_states,1),dxdci_dc(st_state:tot_states,1)] = ...
            dci_sud(i,1,flag,dci_dc{i,1},dcid_sig(i,1),ydcimx,ydcimn);
    end
end

% initialize load modulation control
if ~isempty(g.lmod.lmod_con)
    disp('load modulation')
    lmod(0,1,flag);
end
%
if ~isempty(g.rlmod.rlmod_con)
    disp('reactive load modulation')
    rlmod(0,1,flag);
end

% initialize power modulation control - copied from v2.3 06/01/20 - thad
if ~isempty(g.pwr.pwrmod_con)
    disp('power modulation')
    pwrmod_p(0,1,bus,flag);
    pwrmod_q(0,1,bus,flag);
    [~,~,~,~,Pini,Qini] = pwrmod_dyn([],[],bus,t,0,0);
    if (~iscell(Pini) || ~iscell(Qini))
        estr = 's_simu: error in pwrmod_dyn, P_statesIni and P_statesIni ';
        estr = [estr, 'must be cells.'];
        error(estr);
    end
    if (any(size(Pini)-[g.pwr.n_pwrmod 1]) || any(size(Qini)-[g.pwr.n_pwrmod 1]))
        error('s_simu: dimension error in pwrmod_dyn.');
    end
    pwrmod_p_sigst = cell(g.pwr.n_pwrmod,1);
    pwrmod_q_sigst = pwrmod_p_sigst;
    dpwrmod_p_sigst = pwrmod_p_sigst;
    dpwrmod_q_sigst = pwrmod_p_sigst;
    for index = 1:g.pwr.n_pwrmod
        if ((size(Pini{index},2) ~= 1) || (size(Qini{index},2) ~= 1))
            error('s_simu: dimension error in pwrmod_dyn.');
        end
        dpwrmod_p_sigst{index} = zeros(size(Pini{index},1),k);
        pwrmod_p_sigst{index} = Pini{index}*ones(1,k);
        dpwrmod_q_sigst{index} = zeros(size(Qini{index},1),k);
        pwrmod_q_sigst{index} = Qini{index}*ones(1,k);
    end
    clear('index','dp','dq','Pini','Qini');
end

% initialize non-conforming loads
if ~isempty(g.ncl.load_con)
    disp('non-linear loads')
    vnc = nc_load(bus,flag,Y_ncprf,Y_ncgprf);
end

% reinitialize hvdc controls
if ~isempty(g.dc.dcsp_con)
    disp('dc converter specification')
    bus_sim = bus;
    g.bus.bus_int = bus_intprf;
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
    dc_line(0,1,flag);
end

H_sum = sum(g.mac.mac_con(:,16)./g.mac.mac_pot(:,1));

% step 3: perform a predictor-corrector integration
kt = 0;
ks = 1;
k_tot = sum(k_inc);
lswitch = length(k_inc);
ktmax = k_tot - k_inc(lswitch);
bus_sim = bus;

%% Starting simulation loop
warning('*** Starting simulation loop')
while (kt <= ktmax)
    k_start = kt + 1;

    if (kt == ktmax)
        k_end = kt + k_inc(ks);
    else
        k_end = kt + k_inc(ks) + 1;
    end

    for k = k_start:k_end
        % step 3a: network solution

        % g.sys.mach_ref(k) = g.mac.mac_ang(g.sys.syn_ref,k);
        g.sys.mach_ref(k) = 0;

        g.mac.pmech(:,k+1) = g.mac.pmech(:,k);
        g.igen.tmig(:,k+1) = g.igen.tmig(:,k);

        if (g.dc.n_conv ~= 0)
            g.dc.cur_ord(:,k+1) = g.dc.cur_ord(:,k);
        end

        flag = 1;

        % required to track actual ivm speed
        mac_ang_ivm_old = g.mac.mac_ang(g.mac.mac_ivm_idx,k);

        % network-machine interface for generators
        mac_ind(0,k,bus_sim,flag);
        mac_igen(0,k,bus_sim,flag);
        mac_sub(0,k,bus_sim,flag);
        mac_tra(0,k,bus_sim,flag);
        mac_em(0,k,bus_sim,flag);

        % network interface for internal voltage models
        gfma(0,k,bus_sim,flag);

        if (g.ivm.n_ivmud ~= 0)
            xivm_dc = ivm_sud(0,k,bus_sim,flag,xivm_dc);
        end

        mac_ivm(0,k,bus_sim,flag);

        % network interface for hvdc
        mdc_sig(t(k),k);
        dc_cont(0,k,10*(k-1)+1,bus_sim,flag);

        % Calculate current injections and bus voltages and angles
        if (k >= sum(k_inc(1:3))+1)
            % fault cleared
            line_sim = line_pf2;
            bus_sim = bus_pf2;
            g.bus.bus_int = bus_intpf2;
            Y1 = Y_gpf2;
            Y2 = Y_gncpf2;
            Y3 = Y_ncgpf2;
            Y4 = Y_ncpf2;
            Vr1 = V_rgpf2;
            Vr2 = V_rncpf2;
            bo = bopf2;
        elseif (k >= sum(k_inc(1:2))+1)
            % near bus cleared
            line_sim = line_pf1;
            bus_sim = bus_pf1;
            g.bus.bus_int = bus_intpf1;
            Y1 = Y_gpf1;
            Y2 = Y_gncpf1;
            Y3 = Y_ncgpf1;
            Y4 = Y_ncpf1;
            Vr1 = V_rgpf1;
            Vr2 = V_rncpf1;
            bo = bopf1;
        elseif (k >= k_inc(1)+1)
            % fault applied
            line_sim = line_f;
            bus_sim = bus_f;
            g.bus.bus_int = bus_intf;
            Y1 = Y_gf;
            Y2 = Y_gncf;
            Y3 = Y_ncgf;
            Y4 = Y_ncf;
            Vr1 = V_rgf;
            Vr2 = V_rncf;
            bo = bof;
        elseif (k < k_inc(1)+1)
            % pre fault
            line_sim = line;
            bus_sim = bus;
            g.bus.bus_int = bus_intprf;
            Y1 = Y_gprf;
            Y2 = Y_gncprf;
            Y3 = Y_ncgprf;
            Y4 = Y_ncprf;
            Vr1 = V_rgprf;
            Vr2 = V_rncprf;
            bo = boprf;
        end

        if g.trip.enable
            [bus_sim,line_sim] = trip_handler(t(k),k,true,bus_sim,line_sim);
            [Y1,Y2,Y3,Y4,Vr1,Vr2,bo] = red_ybus(bus_sim,line_sim);
        end

        % form the network interface variables
        h_sol = i_simu(k,ks,k_inc,h,bus_sim,Y1,Y2,Y3,Y4,Vr1,Vr2,bo);

        % checking for gfma current limit violations
        if (g.gfma.n_gfma ~= 0)
            if any(g.gfma.lim_flag)
                fprintf('entered gfma current limiting at time index %0.0f\n', k);
                for jj = 1:50
                    gfma(0,k,bus_sim,flag);
                    mac_ivm(0,k,bus_sim,flag);
                    h_sol = i_simu(k,ks,k_inc,h,bus_sim,Y1,Y2,Y3,Y4,Vr1,Vr2,bo);

                    if ~any(g.gfma.lim_flag)
                        break
                    end
                end
                fprintf('exited after %0.0f iterations\n', jj);
            end
        end

        freqcalc(k,t,1);               % measuring frequency (Felipe Wilches)

        % hvdc
        if (g.dc.ndcr_ud ~= 0)
            % calculate the new value of bus angles rectifier user defined control
            tot_states = 0;
            for jj = 1:g.dc.ndcr_ud
                b_num1 = dcr_dc{jj,3};
                b_num2 = dcr_dc{jj,4};
                conv_num = dcr_dc{jj,2};
                angdcr(jj,k) = g.bus.theta(g.bus.bus_int(b_num1),k) ...
                               - g.bus.theta(g.bus.bus_int(b_num2),k);
                dcrd_sig(jj,k) = angdcr(jj,k);
                st_state = tot_states + 1;
                dcr_states = dcr_dc{jj,7};
                tot_states = tot_states + dcr_states;
                ydcrmx = dcr_dc{jj,5};
                ydcrmn = dcr_dc{jj,6};
                g.dc.dcr_dsig(jj,k) = ...
                    dcr_sud(jj,k,flag,dcr_dc{jj,1},dcrd_sig(jj,k), ...
                            ydcrmx,ydcrmn,xdcr_dc(st_state:tot_states,10*(k-1)+1));
            end
        end
        %
        if (g.dc.ndci_ud ~= 0)
            % calculate the new value of bus angles inverter user defined control
            for jj = 1:g.dc.ndci_ud
                tot_states = 0;
                b_num1 = dci_dc{jj,3};
                b_num2 = dci_dc{jj,4};
                conv_num = dci_dc{jj,2};
                angdci(jj,k) = g.bus.theta(g.bus.bus_int(b_num1),k) ...
                               - g.bus.theta(g.bus.bus_int(b_num2),k);
                dcid_sig(jj,k) = (angdci(jj,k)-angdci(jj,k-1))/(t(k)-t(k-1));
                st_state = tot_states + 1;
                dci_states = dci_dc{jj,7};
                tot_states = tot_states + dci_states;
                ydcimx = dci_dc{jj,5};
                ydcimn = dci_dc{jj,6};
                g.dc.dci_dsig(jj,k) = ...
                    dci_sud(jj,k,flag,dci_dc{jj,1},dcid_sig(jj,k), ...
                            ydcimx,ydcimn,xdci_dc(st_state:tot_states,10*(k-1)+1));
            end
        end

        % network interface for control models
        mdc_sig(t(k),k);
        dc_cont(0,k,10*(k-1)+1,bus_sim,flag);

        dpwf(0,k,flag);
        pss(0,k,flag);

        mexc_sig(t(k),k);
        smpexc(0,k,flag);
        smppi(0,k,flag);
        exc_st3(0,k,flag);
        exc_dc12(0,k,flag);

        mtg_sig(t(k),k);
        tg(0,k,flag);
        tg_hydro(0,k,flag);

        if (g.svc.n_svcud ~= 0)
            % set the new line currents
            for jj = 1:g.svc.n_svcud
                l_num = svc_dc{jj,3};
                svc_num = svc_dc{jj,2};
                from_bus = g.bus.bus_int(line_sim(l_num,1));
                to_bus = g.bus.bus_int(line_sim(l_num,2));
                svc_bn = g.bus.bus_int(g.svc.svc_con(svc_num,2));
                V1 = g.bus.bus_v(from_bus,k);
                V2 = g.bus.bus_v(to_bus,k);
                R = line_sim(l_num,3);
                X = line_sim(l_num,4);
                B = line_sim(l_num,5);
                tap = line_sim(l_num,6);
                phi = line_sim(l_num,7);

                [l_if,l_it] = line_cur(V1,V2,R,X,B,tap,phi);
                if (svc_bn == from_bus)
                    d_sig(jj,k) = abs(l_if);
                elseif (svc_bn == to_bus)
                    d_sig(jj,k) = abs(l_it);
                end
            end
        end

        if (g.tcsc.n_tcscud ~= 0)
            % set the new bus voltages
            for jj = 1:g.tcsc.n_tcscud
                b_num = tcsc_dc{jj,3};
                tcsc_num = tcsc_dc{jj,2};
                td_sig(jj,k) = abs(g.bus.bus_v(g.bus.bus_int(b_num),k));
            end
        end

        i_plot = k - plot_now;
        if (i_plot == 10)
            plot_now = k;
            v_p(1:k) = abs(g.bus.bus_v(bus_idx(1),1:k));
            % plot the voltage of the faulted bus
            plot(t(1:k),v_p(1:k),'r')
            tdfile = strrep(dfile,'_','\_');
            title(['Voltage Magnitude at ', num2str(bus(bus_idx(1),1)), ...
                   ' ', tdfile]);
            xlabel('time (s)');
            drawnow
        end

        g.sys.sys_freq(k) = 1.0;

        % step 3b: compute dynamics and integrate
        flag = 2;

        mpm_sig(t(k),k);
        mac_ind(0,k,bus_sim,flag);
        mac_igen(0,k,bus_sim,flag);
        mac_sub(0,k,bus_sim,flag);
        mac_tra(0,k,bus_sim,flag);
        mac_em(0,k,bus_sim,flag);
        mac_ivm(0,k,bus_sim,flag);

        dpwf(0,k,flag);
        pss(0,k,flag);

        mexc_sig(t(k),k);
        smpexc(0,k,flag);
        smppi(0,k,flag);
        exc_st3(0,k,flag);
        exc_dc12(0,k,flag);

        mtg_sig(t(k),k);
        tg(0,k,flag);
        tg_hydro(0,k,flag);

        if (g.svc.n_svc ~= 0)
            v_svc = abs(g.bus.bus_v(g.bus.bus_int(g.svc.svc_con(:,2)),k));
            if (g.svc.n_svcud ~= 0)
                tot_states = 0;
                for jj = 1:g.svc.n_svcud
                    ysvcmx = svc_dc{jj,4};
                    ysvcmn = svc_dc{jj,5};
                    svc_num = svc_dc{jj,2};
                    st_state = tot_states + 1;
                    svc_states = svc_dc{jj,6};
                    tot_states = tot_states + svc_states;
                    [g.svc.svc_dsig(svc_num,k),xsvc_dc(st_state:tot_states,k),dxsvc_dc(st_state:tot_states,k)] = ...
                        svc_sud(jj,k,flag,svc_dc{jj,1},d_sig(jj,k),ysvcmx,ysvcmn,xsvc_dc(st_state:tot_states,k));
                end
            end
            msvc_sig(t(k),k);
            svc(0,k,bus_sim,flag,v_svc);
        end

        if (g.tcsc.n_tcsc ~= 0)
            if (g.tcsc.n_tcscud ~= 0)
                tot_states = 0;
                for jj = 1:g.tcsc.n_tcscud
                    ytcscmx = tcsc_dc{jj,4};
                    ytcscmn = tcsc_dc{jj,5};
                    tcsc_num = tcsc_dc{jj,2};
                    st_state = tot_states + 1;
                    tcsc_states = tcsc_dc{jj,6};
                    tot_states = tot_states + tcsc_states;
                    [g.tcsc.tcsc_dsig(tcsc_num,k),xtcsc_dc(st_state:tot_states,k),dxtcsc_dc(st_state:tot_states,k)] = ...
                        tcsc_sud(jj,k,flag,tcsc_dc{jj,1},td_sig(jj,k),ytcscmx,ytcscmn,xtcsc_dc(st_state:tot_states,k));
                end
            end
            mtcsc_sig(t(k),k);
            tcsc(0,k,flag);
        end

        if (g.ess.n_ess ~= 0)
            mess_sig(t(k),k);
            mreec_sig(t(k),k);
            lsc(0,k,bus_sim,flag);
            reec(0,k,bus_sim,flag,h_sol);
            xess_dc = ess_sud(0,k,bus_sim,flag,xess_dc);
            ess(0,k,bus_sim,flag,h_sol);
        end

        if (g.lmod.n_lmod ~= 0)
            mlmod_sig(t(k),k);
            lmod(0,k,flag);
        end

        if (g.rlmod.n_rlmod ~= 0)
            mrlmod_sig(t(k),k);
            rlmod(0,k,flag);
        end

        % pwrmod - copied from v2.3 - 06/01/20 - thad
        if (g.pwr.n_pwrmod ~= 0)
            Pst = cell(g.pwr.n_pwrmod,1);
            Qst = Pst;
            for index = 1:g.pwr.n_pwrmod
                Pst{index} = pwrmod_p_sigst{index}(:,k);
                Qst{index} = pwrmod_q_sigst{index}(:,k);
            end

            [~,~,dp,dq,~,~] = pwrmod_dyn(Pst,Qst,bus,t,k,flag);
            if (~iscell(dp) || ~iscell(dq))
                error('s_simu: error in pwrmod_dyn, dp and dq must be cells.');
            end

            if (any(size(dp)-[g.pwr.n_pwrmod 1]) || any(size(dq)-[g.pwr.n_pwrmod 1]))
                error('s_simu: dimension error in pwrmod_dyn.');
            end

            for index = 1:g.pwr.n_pwrmod
                if (size(dp{index},2) ~= 1 || size(dq{index},2) ~= 1)
                    error('s_simu: dimension error in pwrmod_dyn.');
                end
                if (size(dp{index},1) ~= size(dpwrmod_p_sigst{index},1))
                    error('s_simu: dimension error in pwrmod_dyn.');
                end
                if (size(dq{index},1) ~= size(dpwrmod_q_sigst{index},1))
                    error('s_simu: dimension error in pwrmod_dyn.');
                end
                dpwrmod_p_sigst{index}(:,k) = dp{index};
                dpwrmod_q_sigst{index}(:,k) = dq{index};
            end

            % update pwrmod_p_sig and pwrmod_q_sig
            [P,Q,~,~] = pwrmod_dyn(Pst,Qst,bus,t,k,1);
            if (length(P) ~= g.pwr.n_pwrmod || length(Q) ~= g.pwr.n_pwrmod)
                error('s_simu: dimension error in pwrmod_dyn.');
            end

            g.pwr.pwrmod_p_sig(:,k) = P;
            g.pwr.pwrmod_q_sig(:,k) = Q;
            pwrmod_p(0,k,bus_sim,flag);
            pwrmod_q(0,k,bus_sim,flag);
            clear('P','Q','Pst','Qst','dp','dq','index');
        end

        % ivm control models
        if (g.ivm.n_ivm ~= 0)
            mgfma_sig(t(k),k);
            gfma(0,k,bus_sim,flag);
            if (g.ivm.n_ivmud ~= 0)
                xivm_dc = ivm_sud(0,k,bus_sim,flag,xivm_dc);
            end
        end

        % integrate hvdc at ten times rate
        if (g.dc.n_conv ~= 0)
            hdc_sol = h_sol/10;
            for kk = 1:10
                kdc = 10*(k-1) + kk;
                [xdcr_dc(:,kdc:kdc+1),dxdcr_dc(:,kdc:kdc+1),xdci_dc(:,kdc:kdc+1),dxdci_dc(:,kdc:kdc+1)] = ...
                    dc_sim(k,kk,dcr_dc,dci_dc,dcrd_sig,dcid_sig,xdcr_dc(:,kdc),xdci_dc(:,kdc),bus_sim,hdc_sol);
            end
        else
            dc_cont(0,k,k,bus_sim,2);
            dc_line(0,k,2);
        end

        j = k + 1;

        % following statements are predictor steps
        %
        % mac
        %
        g.mac.mac_ang(:,j) = g.mac.mac_ang(:,k) + h_sol*g.mac.dmac_ang(:,k);
        g.mac.mac_spd(:,j) = g.mac.mac_spd(:,k) + h_sol*g.mac.dmac_spd(:,k);
        g.mac.edprime(:,j) = g.mac.edprime(:,k) + h_sol*g.mac.dedprime(:,k);
        g.mac.eqprime(:,j) = g.mac.eqprime(:,k) + h_sol*g.mac.deqprime(:,k);
        g.mac.psikd(:,j) = g.mac.psikd(:,k) + h_sol*g.mac.dpsikd(:,k);
        g.mac.psikq(:,j) = g.mac.psikq(:,k) + h_sol*g.mac.dpsikq(:,k);
        %
        % exc
        %
        g.exc.Efd(:,j) = g.exc.Efd(:,k) + h_sol*g.exc.dEfd(:,k);
        g.exc.R_f(:,j) = g.exc.R_f(:,k) + h_sol*g.exc.dR_f(:,k);
        g.exc.V_R(:,j) = g.exc.V_R(:,k) + h_sol*g.exc.dV_R(:,k);
        g.exc.V_As(:,j) = g.exc.V_As(:,k) + h_sol*g.exc.dV_As(:,k);
        g.exc.V_TR(:,j) = g.exc.V_TR(:,k) + h_sol*g.exc.dV_TR(:,k);
        %
        % dpw
        %
        g.dpw.sdpw1(:,j) = g.dpw.sdpw1(:,k) + h_sol*g.dpw.dsdpw1(:,k);
        g.dpw.sdpw2(:,j) = g.dpw.sdpw2(:,k) + h_sol*g.dpw.dsdpw2(:,k);
        g.dpw.sdpw3(:,j) = g.dpw.sdpw3(:,k) + h_sol*g.dpw.dsdpw3(:,k);
        g.dpw.sdpw4(:,j) = g.dpw.sdpw4(:,k) + h_sol*g.dpw.dsdpw4(:,k);
        g.dpw.sdpw5(:,j) = g.dpw.sdpw5(:,k) + h_sol*g.dpw.dsdpw5(:,k);
        g.dpw.sdpw6(:,j) = g.dpw.sdpw6(:,k) + h_sol*g.dpw.dsdpw6(:,k);
        %
        % pss
        %
        g.pss.pss1(:,j) = g.pss.pss1(:,k) + h_sol*g.pss.dpss1(:,k);
        g.pss.pss2(:,j) = g.pss.pss2(:,k) + h_sol*g.pss.dpss2(:,k);
        g.pss.pss3(:,j) = g.pss.pss3(:,k) + h_sol*g.pss.dpss3(:,k);
        %
        % tg
        %
        g.tg.tg1(:,j) = g.tg.tg1(:,k) + h_sol*g.tg.dtg1(:,k);
        g.tg.tg2(:,j) = g.tg.tg2(:,k) + h_sol*g.tg.dtg2(:,k);
        g.tg.tg3(:,j) = g.tg.tg3(:,k) + h_sol*g.tg.dtg3(:,k);
        g.tg.tg4(:,j) = g.tg.tg4(:,k) + h_sol*g.tg.dtg4(:,k);
        g.tg.tg5(:,j) = g.tg.tg5(:,k) + h_sol*g.tg.dtg5(:,k);
        %
        % ind
        %
        g.ind.vdp(:,j) = g.ind.vdp(:,k) + h_sol*g.ind.dvdp(:,k);
        g.ind.vqp(:,j) = g.ind.vqp(:,k) + h_sol*g.ind.dvqp(:,k);
        g.ind.slip(:,j) = g.ind.slip(:,k) + h_sol*g.ind.dslip(:,k);
        %
        % igen
        %
        g.igen.vdpig(:,j) = g.igen.vdpig(:,k) + h_sol*g.igen.dvdpig(:,k);
        g.igen.vqpig(:,j) = g.igen.vqpig(:,k) + h_sol*g.igen.dvqpig(:,k);
        g.igen.slig(:,j) = g.igen.slig(:,k) + h_sol*g.igen.dslig(:,k);
        %
        % svc
        %
        g.svc.B_cv(:,j) = g.svc.B_cv(:,k) + h_sol*g.svc.dB_cv(:,k);
        g.svc.B_con(:,j) = g.svc.B_con(:,k) + h_sol*g.svc.dB_con(:,k);
        xsvc_dc(:,j) = xsvc_dc(:,k) + h_sol*dxsvc_dc(:,k);
        %
        % tcsc
        %
        g.tcsc.B_tcsc(:,j) = g.tcsc.B_tcsc(:,k) + h_sol*g.tcsc.dB_tcsc(:,k);
        xtcsc_dc(:,j) = xtcsc_dc(:,k) + h_sol*dxtcsc_dc(:,k);
        %
        % lsc
        %
        if (g.lsc.n_lsc ~= 0)
            g.lsc.lsc1(:,j) = g.lsc.lsc1(:,k) + h_sol*g.lsc.dlsc1(:,k);
            g.lsc.lsc2(:,j) = g.lsc.lsc2(:,k) + h_sol*g.lsc.dlsc2(:,k);
            g.lsc.lsc3(:,j) = g.lsc.lsc3(:,k) + h_sol*g.lsc.dlsc3(:,k);
            g.lsc.lsc4(:,j) = g.lsc.lsc4(:,k) + h_sol*g.lsc.dlsc4(:,k);
            g.lsc.lsc5(:,j) = g.lsc.lsc5(:,k) + h_sol*g.lsc.dlsc5(:,k);
            g.lsc.lsc6(:,j) = g.lsc.lsc6(:,k) + h_sol*g.lsc.dlsc6(:,k);
            g.lsc.lsc7(:,j) = g.lsc.lsc7(:,k) + h_sol*g.lsc.dlsc7(:,k);
            g.lsc.lsc8(:,j) = g.lsc.lsc8(:,k) + h_sol*g.lsc.dlsc8(:,k);
            g.lsc.lsc9(:,j) = g.lsc.lsc9(:,k) + h_sol*g.lsc.dlsc9(:,k);
            g.lsc.lsc10(:,j) = g.lsc.lsc10(:,k) + h_sol*g.lsc.dlsc10(:,k);
            g.lsc.lsc11(:,j) = g.lsc.lsc11(:,k) + h_sol*g.lsc.dlsc11(:,k);
            g.lsc.lsc12(:,j) = g.lsc.lsc12(:,k) + h_sol*g.lsc.dlsc12(:,k);
            g.lsc.lsc13(:,j) = g.lsc.lsc13(:,k) + h_sol*g.lsc.dlsc13(:,k);
            g.lsc.lsc14(:,j) = g.lsc.lsc14(:,k) + h_sol*g.lsc.dlsc14(:,k);
            g.lsc.lsc15(:,j) = g.lsc.lsc15(:,k) + h_sol*g.lsc.dlsc15(:,k);
        end
        %
        % reec
        %
        if (g.reec.n_reec ~= 0)
            g.reec.reec1(:,j) = g.reec.reec1(:,k) + h_sol*g.reec.dreec1(:,k);
            g.reec.reec2(:,j) = g.reec.reec2(:,k) + h_sol*g.reec.dreec2(:,k);
            g.reec.reec3(:,j) = g.reec.reec3(:,k) + h_sol*g.reec.dreec3(:,k);
            g.reec.reec4(:,j) = g.reec.reec4(:,k) + h_sol*g.reec.dreec4(:,k);
            g.reec.reec5(:,j) = g.reec.reec5(:,k) + h_sol*g.reec.dreec5(:,k);
            g.reec.reec6(:,j) = g.reec.reec6(:,k) + h_sol*g.reec.dreec6(:,k);
            g.reec.reec7(:,j) = g.reec.reec7(:,k) + h_sol*g.reec.dreec7(:,k);
            g.reec.reec8(:,j) = g.reec.reec8(:,k) + h_sol*g.reec.dreec8(:,k);
            g.reec.reec9(:,j) = g.reec.reec9(:,k) + h_sol*g.reec.dreec9(:,k);
            g.reec.reec10(:,j) = g.reec.reec10(:,k) + h_sol*g.reec.dreec10(:,k);
        end
        %
        % ess
        %
        if (g.ess.n_ess ~= 0)
            g.ess.ess1(:,j) = g.ess.ess1(:,k) + h_sol*g.ess.dess1(:,k);
            g.ess.ess2(:,j) = g.ess.ess2(:,k) + h_sol*g.ess.dess2(:,k);
            g.ess.ess3(:,j) = g.ess.ess3(:,k) + h_sol*g.ess.dess3(:,k);
            g.ess.ess4(:,j) = g.ess.ess4(:,k) + h_sol*g.ess.dess4(:,k);
            g.ess.ess5(:,j) = g.ess.ess5(:,k) + h_sol*g.ess.dess5(:,k);
        end

        if ((g.ess.n_ess ~= 0) && (g.ess.n_essud ~= 0))
            for idx = 1:numel(xess_dc.s)
                xess_dc.s{idx}(:,j) = xess_dc.s{idx}(:,k) ...
                                      + h_sol*xess_dc.ds{idx}(:,k);
            end
        end
        %
        % lmod, rlmod
        %
        g.lmod.lmod_st(:,j) = g.lmod.lmod_st(:,k) + h_sol*g.lmod.dlmod_st(:,k);
        g.rlmod.rlmod_st(:,j) = g.rlmod.rlmod_st(:,k) + h_sol*g.rlmod.drlmod_st(:,k);
        %
        % freqcalc
        %
        g.freq.bf_hpf(:,j) = g.freq.bf_hpf(:,k) + h_sol*g.freq.dbf_hpf(:,k);
        g.freq.x1_snlf(:,j) = g.freq.x1_snlf(:,k) + h_sol*g.freq.dx1_snlf(:,k);
        g.freq.x2_snlf(:,j) = g.freq.x2_snlf(:,k) + h_sol*g.freq.dx2_snlf(:,k);
        %
        % pwrmod_p, pwrmod_q
        %
        g.pwr.pwrmod_p_st(:,j) = g.pwr.pwrmod_p_st(:,k) ...
                                 + h_sol*g.pwr.dpwrmod_p_st(:,k);
        g.pwr.pwrmod_q_st(:,j) = g.pwr.pwrmod_q_st(:,k) ...
                                 + h_sol*g.pwr.dpwrmod_q_st(:,k);
        %
        if (g.pwr.n_pwrmod ~= 0)
            for index = 1:g.pwr.n_pwrmod
                pwrmod_p_sigst{index}(:,j) = pwrmod_p_sigst{index}(:,k) ...
                                             + h_sol*dpwrmod_p_sigst{index}(:,k);
                pwrmod_q_sigst{index}(:,j) = pwrmod_q_sigst{index}(:,k) ...
                                             + h_sol*dpwrmod_q_sigst{index}(:,k);
            end
        end
        %
        % gfma
        %
        if (g.gfma.n_gfma ~= 0)
            g.gfma.gfma1(:,j) = g.gfma.gfma1(:,k) + h_sol*g.gfma.dgfma1(:,k);
            g.gfma.gfma2(:,j) = g.gfma.gfma2(:,k) + h_sol*g.gfma.dgfma2(:,k);
            g.gfma.gfma3(:,j) = g.gfma.gfma3(:,k) + h_sol*g.gfma.dgfma3(:,k);
            g.gfma.gfma4(:,j) = g.gfma.gfma4(:,k) + h_sol*g.gfma.dgfma4(:,k);
            g.gfma.gfma5(:,j) = g.gfma.gfma5(:,k) + h_sol*g.gfma.dgfma5(:,k);
            g.gfma.gfma6(:,j) = g.gfma.gfma6(:,k) + h_sol*g.gfma.dgfma6(:,k);
            g.gfma.gfma7(:,j) = g.gfma.gfma7(:,k) + h_sol*g.gfma.dgfma7(:,k);
            g.gfma.gfma8(:,j) = g.gfma.gfma8(:,k) + h_sol*g.gfma.dgfma8(:,k);
            g.gfma.gfma9(:,j) = g.gfma.gfma9(:,k) + h_sol*g.gfma.dgfma9(:,k);
            g.gfma.gfma10(:,j) = g.gfma.gfma10(:,k) + h_sol*g.gfma.dgfma10(:,k);
        end
        %
        % ivm
        %
        if (g.mac.mac_ang(g.mac.mac_ivm_idx,j) ~= mac_ang_ivm_old)
            g.mac.mac_spd(g.mac.mac_ivm_idx,k) = ...
               (g.mac.mac_ang(g.mac.mac_ivm_idx,j) ...
                - mac_ang_ivm_old)./(h_sol*g.sys.basrad) ...
               + g.mac.mac_spd(g.mac.mac_ivm_idx,1);
        else
            g.mac.mac_spd(g.mac.mac_ivm_idx,k) = ...
                g.mac.mac_spd(g.mac.mac_ivm_idx,max(k-1,1));
        end
        %
        % ivm_sud
        %
        if (g.ivm.n_ivmud ~= 0)
            for idx = 1:numel(xivm_dc.s)
                xivm_dc.s{idx}(:,j) = xivm_dc.s{idx}(:,k) ...
                                      + h_sol*xivm_dc.ds{idx}(:,k);
            end
        end

        % end of predictor steps

        %---------------------------------------------------------------------%
        flag = 1;

        % g.sys.mach_ref(j) = g.mac.mac_ang(g.sys.syn_ref,j);
        g.sys.mach_ref(j) = 0;

        % updating ivm commands to allow time constants to be bypassed
        if (g.gfma.n_gfma ~= 0)
            g.mac.vex(g.mac.mac_ivm_idx,j) = g.mac.vex(g.mac.mac_ivm_idx,k);
            g.mac.fldcur(g.mac.mac_ivm_idx,j) = g.mac.fldcur(g.mac.mac_ivm_idx,k);

            busnum = g.bus.bus_int(g.gfma.gfma_con(:,2));
            g.bus.bus_v(busnum,j) = g.bus.bus_v(busnum,k);
        end

        % perform network interface calculations again with predicted states
        mac_ind(0,j,bus_sim,flag);
        mac_igen(0,j,bus_sim,flag);
        mac_sub(0,j,bus_sim,flag);
        mac_tra(0,j,bus_sim,flag);
        mac_em(0,j,bus_sim,flag);

        gfma(0,j,bus_sim,flag);

        if (g.ivm.n_ivmud ~= 0)
            xivm_dc = ivm_sud(0,j,bus_sim,flag,xivm_dc);
        end

        mac_ivm(0,j,bus_sim,flag);

        % assume Vdc remains unchanged for first pass through dc controls interface
        mdc_sig(t(j),j);
        dc_cont(0,j,10*(j-1)+1,bus_sim,flag);

        % Calculate current injections and bus voltages and angles
        if (j >= sum(k_inc(1:3))+1)
            % fault cleared
            bus_sim = bus_pf2;
            g.bus.bus_int = bus_intpf2;
            Y1 = Y_gpf2;
            Y2 = Y_gncpf2;
            Y3 = Y_ncgpf2;
            Y4 = Y_ncpf2;
            Vr1 = V_rgpf2;
            Vr2 = V_rncpf2;
            bo = bopf2;
        elseif (j >= sum(k_inc(1:2))+1)
            % near bus cleared
            bus_sim = bus_pf1;
            g.bus.bus_int = bus_intpf1;
            Y1 = Y_gpf1;
            Y2 = Y_gncpf1;
            Y3 = Y_ncgpf1;
            Y4 = Y_ncpf1;
            Vr1 = V_rgpf1;
            Vr2 = V_rncpf1;
            bo = bopf1;
        elseif (j >= k_inc(1)+1)
            % fault applied
            bus_sim = bus_f;
            g.bus.bus_int = bus_intf;
            Y1 = Y_gf;
            Y2 = Y_gncf;
            Y3 = Y_ncgf;
            Y4 = Y_ncf;
            Vr1 = V_rgf;
            Vr2 = V_rncf;
            bo = bof;
        elseif (j < k_inc(1)+1)  % JHC changed k to j per DKF
            % pre fault
            bus_sim = bus;
            g.bus.bus_int = bus_intprf;
            Y1 = Y_gprf;
            Y2 = Y_gncprf;
            Y3 = Y_ncgprf;
            Y4 = Y_ncprf;
            Vr1 = V_rgprf;
            Vr2 = V_rncprf;
            bo = boprf;
        end

        if g.trip.enable
            [bus_sim,line_sim] = trip_handler(t(j),j,false,bus_sim,line_sim);
            [Y1,Y2,Y3,Y4,Vr1,Vr2,bo] = red_ybus(bus_sim,line_sim);
        end

        % form the network interface variables
        h_sol = i_simu(j,ks,k_inc,h,bus_sim,Y1,Y2,Y3,Y4,Vr1,Vr2,bo);

        freqcalc(j,t,1);               % measuring frequency (Felipe Wilches)

        g.mac.vex(:,j) = g.mac.vex(:,k);
        g.dc.cur_ord(:,j) = g.dc.cur_ord(:,k);

        % hvdc
        if (g.dc.ndcr_ud ~= 0)
            % calculate the new value of bus angles rectifier user defined control
            tot_states = 0;
            for jj = 1:g.dc.ndcr_ud
                b_num1 = dcr_dc{jj,3};
                b_num2 = dcr_dc{jj,4};
                conv_num = dcr_dc{jj,2};
                angdcr(jj,j) = g.bus.theta(g.bus.bus_int(b_num1),j) ...
                               - g.bus.theta(g.bus.bus_int(b_num2),j);
                dcrd_sig(jj,j) = angdcr(jj,j);
                st_state = tot_states + 1;
                dcr_states = dcr_dc{jj,7};
                tot_states = tot_states + dcr_states;
                ydcrmx = dcr_dc{jj,5};
                ydcrmn = dcr_dc{jj,6};
                g.dc.dcr_dsig(jj,j) = ...
                    dcr_sud(jj,j,flag,dcr_dc{jj,1},dcrd_sig(jj,j), ...
                            ydcrmx,ydcrmn,xdcr_dc(st_state:tot_states,10*(j-1)+1));
            end
        end

        if (g.dc.ndci_ud ~= 0)
            % calculate the new value of bus angles inverter user defined control
            for jj = 1:g.dc.ndci_ud
                tot_states = 0;
                b_num1 = dci_dc{jj,3};
                b_num2 = dci_dc{jj,4};
                conv_num = dci_dc{jj,2};
                angdci(jj,j) = g.bus.theta(g.bus.bus_int(b_num1),j) ...
                               - g.bus.theta(g.bus.bus_int(b_num2),j);
                dcid_sig(jj,j) = (angdci(jj,j)-angdci(jj,k))/(t(j)-t(k));
                st_state = tot_states + 1;
                dci_states = dci_dc{jj,7};
                tot_states = tot_states + dci_states;
                ydcimx = dci_dc{jj,5};
                ydcimn = dci_dc{jj,6};
                g.dc.dci_dsig(jj,j) = ...
                    dci_sud(jj,j,flag,dci_dc{jj,1},dcid_sig(jj,j), ...
                            ydcimx,ydcimn,xdci_dc(st_state:tot_states,10*(j-1)+1));
            end
        end

        mdc_sig(t(j),j);
        dc_cont(0,j,10*(j-1)+1,bus_sim,flag);

        dpwf(0,j,flag);
        pss(0,j,flag);

        mexc_sig(t(j),j);
        smpexc(0,j,flag);
        smppi(0,j,flag);
        exc_st3(0,j,flag);
        exc_dc12(0,j,flag);

        mtg_sig(t(j),j);
        tg(0,j,flag);
        tg_hydro(0,j,flag);

        if (g.svc.n_svcud ~= 0)
            % set the new line currents
            for jj = 1:g.svc.n_svcud
                l_num = svc_dc{jj,3};
                svc_num = svc_dc{jj,2};
                from_bus = g.bus.bus_int(line_sim(l_num,1));
                to_bus = g.bus.bus_int(line_sim(l_num,2));
                svc_bn = g.bus.bus_int(g.svc.svc_con(svc_num,2));
                V1 = g.bus.bus_v(from_bus,j);
                V2 = g.bus.bus_v(to_bus,j);
                R = line_sim(l_num,3);
                X = line_sim(l_num,4);
                B = line_sim(l_num,5);
                tap = line_sim(l_num,6);
                phi = line_sim(l_num,7);

                [l_if,l_it] = line_cur(V1,V2,R,X,B,tap,phi);
                if (svc_bn == from_bus)
                    d_sig(jj,j) = abs(l_if);
                elseif (svc_bn == to_bus)
                    d_sig(jj,j) = abs(l_it);
                end
            end
        end

        if (g.tcsc.n_tcscud ~= 0)
            % set the new bus voltages
            for jj = 1:g.tcsc.n_tcscud
                b_num = tcsc_dc{jj,3};
                tcsc_num = tcsc_dc{jj,2};
                td_sig(jj,j) = abs(g.bus.bus_v(g.bus.bus_int(b_num),j));
            end
        end

        % corrector step
        flag = 2;

        mpm_sig(t(j),j);
        mac_ind(0,j,bus_sim,flag);
        mac_igen(0,j,bus_sim,flag);
        mac_sub(0,j,bus_sim,flag);
        mac_tra(0,j,bus_sim,flag);
        mac_em(0,j,bus_sim,flag);
        mac_ivm(0,j,bus_sim,flag);

        dpwf(0,j,flag);
        pss(0,j,flag);

        mexc_sig(t(j),j);
        smpexc(0,j,flag);
        smppi(0,j,flag);
        exc_st3(0,j,flag);
        exc_dc12(0,j,flag);

        mtg_sig(t(j),j);
        tg(0,j,flag);
        tg_hydro(0,j,flag);

        if (g.svc.n_svc ~= 0)
            msvc_sig(t(j),j);
            if (g.svc.n_svcud ~= 0)
                tot_states = 0;
                for jj = 1:g.svc.n_svcud
                    ysvcmx = svc_dc{jj,4};
                    ysvcmn = svc_dc{jj,5};
                    svc_num = svc_dc{jj,2};
                    st_state = tot_states + 1;
                    svc_states = svc_dc{jj,6};
                    tot_states = tot_states + svc_states;
                    [g.svc.svc_dsig(svc_num,j),xsvc_dc(st_state:tot_states,j),dxsvc_dc(st_state:tot_states,j)] = ...
                        svc_sud(jj,j,flag,svc_dc{jj,1},d_sig(jj,j),ysvcmx,ysvcmn,xsvc_dc(st_state:tot_states,j));
                end
            end
            v_svc = abs(g.bus.bus_v(g.bus.bus_int(g.svc.svc_con(:,2)),j));
            bus_sim = svc(0,j,bus_sim,flag,v_svc);
        end

        if (g.tcsc.n_tcsc ~= 0)
            mtcsc_sig(t(j),j);
            if (g.tcsc.n_tcscud ~= 0)
                tot_states = 0;
                for jj = 1:g.tcsc.n_tcscud
                    ytcscmx = tcsc_dc{jj,4};
                    ytcscmn = tcsc_dc{jj,5};
                    tcsc_num = tcsc_dc{jj,2};
                    st_state = tot_states + 1;
                    tcsc_states = tcsc_dc{jj,6};
                    tot_states = tot_states + tcsc_states;
                    [g.tcsc.tcsc_dsig(tcsc_num,j),xtcsc_dc(st_state:tot_states,j),dxtcsc_dc(st_state:tot_states,j)] = ...
                        tcsc_sud(jj,j,flag,tcsc_dc{jj,1},td_sig(jj,j),ytcscmx,ytcscmn,xtcsc_dc(st_state:tot_states,j));
                end
            end
            tcsc(0,j,flag);
        end

        if (g.ess.n_ess ~= 0)
            mess_sig(t(j),j);
            mreec_sig(t(j),j);
            lsc(0,j,bus_sim,flag);
            reec(0,j,bus_sim,flag,h_sol);
            xess_dc = ess_sud(0,j,bus_sim,flag,xess_dc);
            ess(0,j,bus_sim,flag,h_sol);
        end

        if (g.lmod.n_lmod ~= 0)
            mlmod_sig(t(j),j);
            lmod(0,j,flag);
        end

        if (g.rlmod.n_rlmod ~= 0)
            mrlmod_sig(t(j),j);
            rlmod(0,j,flag);
        end

        % pwrmod - copied from v2.3 - 06/01/20 - thad
        if (g.pwr.n_pwrmod ~= 0)
            Pst = cell(g.pwr.n_pwrmod,1);
            Qst = Pst;
            for index = 1:g.pwr.n_pwrmod
                Pst{index} = pwrmod_p_sigst{index}(:,j);
                Qst{index} = pwrmod_q_sigst{index}(:,j);
            end

            [~,~,dp,dq,~,~] = pwrmod_dyn(Pst,Qst,bus,t,j,flag);
            if (~iscell(dp) || ~iscell(dq))
                error('s_simu: error in pwrmod_dyn, dp and dq must be cells.');
            end

            if (any(size(dp)-[g.pwr.n_pwrmod 1]) || any(size(dq)-[g.pwr.n_pwrmod 1]))
                error('s_simu: dimension error in pwrmod_dyn.');
            end

            for index = 1:g.pwr.n_pwrmod
                if (size(dp{index},2) ~= 1 || size(dq{index},2) ~= 1)
                    error('s_simu: dimension error in pwrmod_dyn.');
                end
                if (size(dp{index},1) ~= size(dpwrmod_p_sigst{index},1))
                    error('s_simu: dimension error in pwrmod_dyn.');
                end
                if (size(dq{index},1) ~= size(dpwrmod_q_sigst{index},1))
                    error('s_simu: dimension error in pwrmod_dyn.');
                end
                dpwrmod_p_sigst{index}(:,j) = dp{index};
                dpwrmod_q_sigst{index}(:,j) = dq{index};
            end

            % update pwrmod_p_sig and pwrmod_q_sig
            [P,Q,~,~,~,~] = pwrmod_dyn(Pst,Qst,bus,t,j,1);
            if (length(P) ~= g.pwr.n_pwrmod || length(Q) ~= g.pwr.n_pwrmod)
                error('s_simu: dimension error in pwrmod_dyn.');
            end

            g.pwr.pwrmod_p_sig(:,j) = P;
            g.pwr.pwrmod_q_sig(:,j) = Q;
            pwrmod_p(0,j,bus_sim,flag);
            pwrmod_q(0,j,bus_sim,flag);
            clear('P','Q','Pst','Qst','dp','dq','index');
        end

        % ivm control models
        if (g.ivm.n_ivm ~= 0)
            mgfma_sig(t(j),j);
            gfma(0,j,bus_sim,flag);
            if (g.ivm.n_ivmud ~= 0)
                xivm_dc = ivm_sud(0,j,bus_sim,flag,xivm_dc);
            end
        end

        if (g.dc.n_conv ~= 0)
            hdc_sol = h_sol/10;
            for kk = 1:10
                jdc = 10*(j-1) + kk;
                [xdcr_dc(:,jdc:jdc+1),dxdcr_dc(:,jdc:jdc+1),xdci_dc(:,jdc:jdc+1),dxdci_dc(:,jdc:jdc+1)] = ...
                    dc_sim(j,kk,dcr_dc,dci_dc,dcrd_sig,dcid_sig,xdcr_dc(:,jdc),xdci_dc(:,jdc),bus_sim,hdc_sol);
            end
        else
            dc_cont(0,j,j,bus_sim,2);
            dc_line(0,j,2);
        end

        % following statements are corrector steps
        %
        % mac
        %
        g.mac.mac_ang(:,j) = g.mac.mac_ang(:,k) ...
                             + h_sol*(g.mac.dmac_ang(:,k) + g.mac.dmac_ang(:,j))/2;
        g.mac.mac_spd(:,j) = g.mac.mac_spd(:,k) ...
                             + h_sol*(g.mac.dmac_spd(:,k) + g.mac.dmac_spd(:,j))/2;
        g.mac.edprime(:,j) = g.mac.edprime(:,k) ...
                             + h_sol*(g.mac.dedprime(:,k) + g.mac.dedprime(:,j))/2;
        g.mac.eqprime(:,j) = g.mac.eqprime(:,k) ...
                             + h_sol*(g.mac.deqprime(:,k) + g.mac.deqprime(:,j))/2;
        g.mac.psikd(:,j) = g.mac.psikd(:,k) ...
                           + h_sol*(g.mac.dpsikd(:,k) + g.mac.dpsikd(:,j))/2;
        g.mac.psikq(:,j) = g.mac.psikq(:,k) ...
                           + h_sol*(g.mac.dpsikq(:,k) + g.mac.dpsikq(:,j))/2;
        %
        % exc
        %
        g.exc.Efd(:,j) = g.exc.Efd(:,k) ...
                         + h_sol*(g.exc.dEfd(:,k) + g.exc.dEfd(:,j))/2;
        g.exc.R_f(:,j) = g.exc.R_f(:,k) ...
                         + h_sol*(g.exc.dR_f(:,k) + g.exc.dR_f(:,j))/2;
        g.exc.V_R(:,j) = g.exc.V_R(:,k) ...
                         + h_sol*(g.exc.dV_R(:,k) + g.exc.dV_R(:,j))/2;
        g.exc.V_As(:,j) = g.exc.V_As(:,k) ...
                          + h_sol*(g.exc.dV_As(:,k) + g.exc.dV_As(:,j))/2;
        g.exc.V_TR(:,j) = g.exc.V_TR(:,k) ...
                          + h_sol*(g.exc.dV_TR(:,k) + g.exc.dV_TR(:,j))/2;
        %
        % dpw
        %
        g.dpw.sdpw1(:,j) = g.dpw.sdpw1(:,k) ...
                           + h_sol*(g.dpw.dsdpw1(:,k) + g.dpw.dsdpw1(:,j))/2;
        g.dpw.sdpw2(:,j) = g.dpw.sdpw2(:,k) ...
                           + h_sol*(g.dpw.dsdpw2(:,k) + g.dpw.dsdpw2(:,j))/2;
        g.dpw.sdpw3(:,j) = g.dpw.sdpw3(:,k) ...
                           + h_sol*(g.dpw.dsdpw3(:,k) + g.dpw.dsdpw3(:,j))/2;
        g.dpw.sdpw4(:,j) = g.dpw.sdpw4(:,k) ...
                           + h_sol*(g.dpw.dsdpw4(:,k) + g.dpw.dsdpw4(:,j))/2;
        g.dpw.sdpw5(:,j) = g.dpw.sdpw5(:,k) ...
                           + h_sol*(g.dpw.dsdpw5(:,k) + g.dpw.dsdpw5(:,j))/2;
        g.dpw.sdpw6(:,j) = g.dpw.sdpw6(:,k) ...
                           + h_sol*(g.dpw.dsdpw6(:,k) + g.dpw.dsdpw6(:,j))/2;
        %
        % pss
        %
        g.pss.pss1(:,j) = g.pss.pss1(:,k) ...
                          + h_sol*(g.pss.dpss1(:,k) + g.pss.dpss1(:,j))/2;
        g.pss.pss2(:,j) = g.pss.pss2(:,k) ...
                          + h_sol*(g.pss.dpss2(:,k) + g.pss.dpss2(:,j))/2;
        g.pss.pss3(:,j) = g.pss.pss3(:,k) ...
                          + h_sol*(g.pss.dpss3(:,k) + g.pss.dpss3(:,j))/2;
        %
        % tg
        %
        g.tg.tg1(:,j) = g.tg.tg1(:,k) + h_sol*(g.tg.dtg1(:,k) + g.tg.dtg1(:,j))/2;
        g.tg.tg2(:,j) = g.tg.tg2(:,k) + h_sol*(g.tg.dtg2(:,k) + g.tg.dtg2(:,j))/2;
        g.tg.tg3(:,j) = g.tg.tg3(:,k) + h_sol*(g.tg.dtg3(:,k) + g.tg.dtg3(:,j))/2;
        g.tg.tg4(:,j) = g.tg.tg4(:,k) + h_sol*(g.tg.dtg4(:,k) + g.tg.dtg4(:,j))/2;
        g.tg.tg5(:,j) = g.tg.tg5(:,k) + h_sol*(g.tg.dtg5(:,k) + g.tg.dtg5(:,j))/2;
        %
        % ind
        %
        g.ind.vdp(:,j) = g.ind.vdp(:,k) ...
                         + h_sol*(g.ind.dvdp(:,j) + g.ind.dvdp(:,k))/2;
        g.ind.vqp(:,j) = g.ind.vqp(:,k) ...
                         + h_sol*(g.ind.dvqp(:,j) + g.ind.dvqp(:,k))/2;
        g.ind.slip(:,j) = g.ind.slip(:,k) ...
                          + h_sol*(g.ind.dslip(:,j) + g.ind.dslip(:,k))/2;
        %
        % igen
        %
        g.igen.vdpig(:,j) = g.igen.vdpig(:,k) ...
                            + h_sol*(g.igen.dvdpig(:,j) + g.igen.dvdpig(:,k))/2;
        g.igen.vqpig(:,j) = g.igen.vqpig(:,k) ...
                            + h_sol*(g.igen.dvqpig(:,j) + g.igen.dvqpig(:,k))/2;
        g.igen.slig(:,j) = g.igen.slig(:,k) ...
                           + h_sol*(g.igen.dslig(:,j) + g.igen.dslig(:,k))/2;
        %
        % svc
        %
        g.svc.B_cv(:,j) = g.svc.B_cv(:,k) ...
                          + h_sol*(g.svc.dB_cv(:,j) + g.svc.dB_cv(:,k))/2;
        g.svc.B_con(:,j) = g.svc.B_con(:,k) ...
                           + h_sol*(g.svc.dB_con(:,j) + g.svc.dB_con(:,k))/2;
        xsvc_dc(:,j) = xsvc_dc(:,k) + h_sol*(dxsvc_dc(:,j) + dxsvc_dc(:,k))/2;
        %
        % tcsc
        %
        g.tcsc.B_tcsc(:,j) = g.tcsc.B_tcsc(:,k) ...
                             + h_sol*(g.tcsc.dB_tcsc(:,j) + g.tcsc.dB_tcsc(:,k))/2;
        xtcsc_dc(:,j) = xtcsc_dc(:,k) + h_sol*(dxtcsc_dc(:,j) + dxtcsc_dc(:,k))/2;
        %
        % lsc
        %
        if (g.lsc.n_lsc ~= 0)
            g.lsc.lsc1(:,j) = g.lsc.lsc1(:,k) ...
                              + h_sol*(g.lsc.dlsc1(:,k) + g.lsc.dlsc1(:,j))/2;
            g.lsc.lsc2(:,j) = g.lsc.lsc2(:,k) ...
                              + h_sol*(g.lsc.dlsc2(:,k) + g.lsc.dlsc2(:,j))/2;
            g.lsc.lsc3(:,j) = g.lsc.lsc3(:,k) ...
                              + h_sol*(g.lsc.dlsc3(:,k) + g.lsc.dlsc3(:,j))/2;
            g.lsc.lsc4(:,j) = g.lsc.lsc4(:,k) ...
                              + h_sol*(g.lsc.dlsc4(:,k) + g.lsc.dlsc4(:,j))/2;
            g.lsc.lsc5(:,j) = g.lsc.lsc5(:,k) ...
                              + h_sol*(g.lsc.dlsc5(:,k) + g.lsc.dlsc5(:,j))/2;
            g.lsc.lsc6(:,j) = g.lsc.lsc6(:,k) ...
                              + h_sol*(g.lsc.dlsc6(:,k) + g.lsc.dlsc6(:,j))/2;
            g.lsc.lsc7(:,j) = g.lsc.lsc7(:,k) ...
                              + h_sol*(g.lsc.dlsc7(:,k) + g.lsc.dlsc7(:,j))/2;
            g.lsc.lsc8(:,j) = g.lsc.lsc8(:,k) ...
                              + h_sol*(g.lsc.dlsc8(:,k) + g.lsc.dlsc8(:,j))/2;
            g.lsc.lsc9(:,j) = g.lsc.lsc9(:,k) ...
                              + h_sol*(g.lsc.dlsc9(:,k) + g.lsc.dlsc9(:,j))/2;
            g.lsc.lsc10(:,j) = g.lsc.lsc10(:,k) ...
                               + h_sol*(g.lsc.dlsc10(:,k) + g.lsc.dlsc10(:,j))/2;
            g.lsc.lsc11(:,j) = g.lsc.lsc11(:,k) ...
                               + h_sol*(g.lsc.dlsc11(:,k) + g.lsc.dlsc11(:,j))/2;
            g.lsc.lsc12(:,j) = g.lsc.lsc12(:,k) ...
                               + h_sol*(g.lsc.dlsc12(:,k) + g.lsc.dlsc12(:,j))/2;
            g.lsc.lsc13(:,j) = g.lsc.lsc13(:,k) ...
                               + h_sol*(g.lsc.dlsc13(:,k) + g.lsc.dlsc13(:,j))/2;
            g.lsc.lsc14(:,j) = g.lsc.lsc14(:,k) ...
                               + h_sol*(g.lsc.dlsc14(:,k) + g.lsc.dlsc14(:,j))/2;
            g.lsc.lsc15(:,j) = g.lsc.lsc15(:,k) ...
                               + h_sol*(g.lsc.dlsc15(:,k) + g.lsc.dlsc15(:,j))/2;
        end
        %
        % reec
        %
        if (g.reec.n_reec ~= 0)
            g.reec.reec1(:,j) = g.reec.reec1(:,k) ...
                                + h_sol*(g.reec.dreec1(:,k) + g.reec.dreec1(:,j))/2;
            g.reec.reec2(:,j) = g.reec.reec2(:,k) ...
                                + h_sol*(g.reec.dreec2(:,k) + g.reec.dreec2(:,j))/2;
            g.reec.reec3(:,j) = g.reec.reec3(:,k) ...
                                + h_sol*(g.reec.dreec3(:,k) + g.reec.dreec3(:,j))/2;
            g.reec.reec4(:,j) = g.reec.reec4(:,k) ...
                                + h_sol*(g.reec.dreec4(:,k) + g.reec.dreec4(:,j))/2;
            g.reec.reec5(:,j) = g.reec.reec5(:,k) ...
                                + h_sol*(g.reec.dreec5(:,k) + g.reec.dreec5(:,j))/2;
            g.reec.reec6(:,j) = g.reec.reec6(:,k) ...
                                + h_sol*(g.reec.dreec6(:,k) + g.reec.dreec6(:,j))/2;
            g.reec.reec7(:,j) = g.reec.reec7(:,k) ...
                                + h_sol*(g.reec.dreec7(:,k) + g.reec.dreec7(:,j))/2;
            g.reec.reec8(:,j) = g.reec.reec8(:,k) ...
                                + h_sol*(g.reec.dreec8(:,k) + g.reec.dreec8(:,j))/2;
            g.reec.reec9(:,j) = g.reec.reec9(:,k) ...
                                + h_sol*(g.reec.dreec9(:,k) + g.reec.dreec9(:,j))/2;
            g.reec.reec10(:,j) = g.reec.reec10(:,k) ...
                                + h_sol*(g.reec.dreec10(:,k) + g.reec.dreec10(:,j))/2;
        end
        %
        % ess
        %
        if (g.ess.n_ess ~= 0)
            g.ess.ess1(:,j) = g.ess.ess1(:,k) ...
                              + h_sol*(g.ess.dess1(:,k) + g.ess.dess1(:,j))/2;
            g.ess.ess2(:,j) = g.ess.ess2(:,k) ...
                              + h_sol*(g.ess.dess2(:,k) + g.ess.dess2(:,j))/2;
            g.ess.ess3(:,j) = g.ess.ess3(:,k) ...
                              + h_sol*(g.ess.dess3(:,k) + g.ess.dess3(:,j))/2;
            g.ess.ess4(:,j) = g.ess.ess4(:,k) ...
                              + h_sol*(g.ess.dess4(:,k) + g.ess.dess4(:,j))/2;
            g.ess.ess5(:,j) = g.ess.ess5(:,k) ...
                              + h_sol*(g.ess.dess5(:,k) + g.ess.dess5(:,j))/2;
        end
        %
        if (g.ess.n_essud ~= 0)
            for idx = 1:numel(xess_dc.s)
                xess_dc.s{idx}(:,j) = ...
                    xess_dc.s{idx}(:,k) ...
                    + h_sol*(xess_dc.ds{idx}(:,k) + xess_dc.ds{idx}(:,j))/2;
            end
        end
        %
        % lmod, rlmod
        %
        g.lmod.lmod_st(:,j) = ...
            g.lmod.lmod_st(:,k) ...
            + h_sol*(g.lmod.dlmod_st(:,j) + g.lmod.dlmod_st(:,k))/2;
        g.rlmod.rlmod_st(:,j) = ...
            g.rlmod.rlmod_st(:,k) ...
            + h_sol*(g.rlmod.drlmod_st(:,j) + g.rlmod.drlmod_st(:,k))/2;
        %
        % freqcalc
        %
        g.freq.bf_hpf(:,j) = ...
            g.freq.bf_hpf(:,k) ...
            + h_sol*(g.freq.dbf_hpf(:,j) + g.freq.dbf_hpf(:,k))/2;
        g.freq.x1_snlf(:,j) = ...
            g.freq.x1_snlf(:,k) ...
            + h_sol*(g.freq.dx1_snlf(:,j) + g.freq.dx1_snlf(:,k))/2;
        g.freq.x2_snlf(:,j) = ...
            g.freq.x2_snlf(:,k) ...
            + h_sol*(g.freq.dx2_snlf(:,j) + g.freq.dx2_snlf(:,k))/2;
        %
        % pwrmod_p, pwrmod_q
        %
        g.pwr.pwrmod_p_st(:,j) = ...
            g.pwr.pwrmod_p_st(:,k) ...
            + h_sol*(g.pwr.dpwrmod_p_st(:,j) + g.pwr.dpwrmod_p_st(:,k))/2;

        g.pwr.pwrmod_q_st(:,j) = ...
            g.pwr.pwrmod_q_st(:,k) ...
            + h_sol*(g.pwr.dpwrmod_q_st(:,j) + g.pwr.dpwrmod_q_st(:,k))/2;

        if (g.pwr.n_pwrmod ~= 0)
            for index = 1:g.pwr.n_pwrmod
                pwrmod_p_sigst{index}(:,j) = ...
                    pwrmod_p_sigst{index}(:,k) ...
                    + h_sol*(dpwrmod_p_sigst{index}(:,j) ...
                             + dpwrmod_p_sigst{index}(:,k))/2;

                pwrmod_q_sigst{index}(:,j) = ...
                    pwrmod_q_sigst{index}(:,k) ...
                    + h_sol*(dpwrmod_q_sigst{index}(:,j) ...
                             + dpwrmod_q_sigst{index}(:,k))/2;
            end
        end
        %
        % gfma
        %
        if (g.gfma.n_gfma ~= 0)
            g.gfma.gfma1(:,j) = g.gfma.gfma1(:,k) ...
                              + h_sol*(g.gfma.dgfma1(:,k) + g.gfma.dgfma1(:,j))/2;
            g.gfma.gfma2(:,j) = g.gfma.gfma2(:,k) ...
                              + h_sol*(g.gfma.dgfma2(:,k) + g.gfma.dgfma2(:,j))/2;
            g.gfma.gfma3(:,j) = g.gfma.gfma3(:,k) ...
                              + h_sol*(g.gfma.dgfma3(:,k) + g.gfma.dgfma3(:,j))/2;
            g.gfma.gfma4(:,j) = g.gfma.gfma4(:,k) ...
                              + h_sol*(g.gfma.dgfma4(:,k) + g.gfma.dgfma4(:,j))/2;
            g.gfma.gfma5(:,j) = g.gfma.gfma5(:,k) ...
                              + h_sol*(g.gfma.dgfma5(:,k) + g.gfma.dgfma5(:,j))/2;
            g.gfma.gfma6(:,j) = g.gfma.gfma6(:,k) ...
                              + h_sol*(g.gfma.dgfma6(:,k) + g.gfma.dgfma6(:,j))/2;
            g.gfma.gfma7(:,j) = g.gfma.gfma7(:,k) ...
                              + h_sol*(g.gfma.dgfma7(:,k) + g.gfma.dgfma7(:,j))/2;
            g.gfma.gfma8(:,j) = g.gfma.gfma8(:,k) ...
                              + h_sol*(g.gfma.dgfma8(:,k) + g.gfma.dgfma8(:,j))/2;
            g.gfma.gfma9(:,j) = g.gfma.gfma9(:,k) ...
                              + h_sol*(g.gfma.dgfma9(:,k) + g.gfma.dgfma9(:,j))/2;
            g.gfma.gfma10(:,j) = g.gfma.gfma10(:,k) ...
                               + h_sol*(g.gfma.dgfma10(:,k) + g.gfma.dgfma10(:,j))/2;
        end
        %
        % ivm
        %
        if (g.mac.mac_ang(g.mac.mac_ivm_idx,j) ~= mac_ang_ivm_old)
            g.mac.mac_spd(g.mac.mac_ivm_idx,j) = ...
               2*((g.mac.mac_ang(g.mac.mac_ivm_idx,j) ...
                   - mac_ang_ivm_old)./(h_sol*g.sys.basrad) ...
                  + g.mac.mac_spd(g.mac.mac_ivm_idx,1)) ...
               - g.mac.mac_spd(g.mac.mac_ivm_idx,k);
        else
            g.mac.mac_spd(g.mac.mac_ivm_idx,j) = ...
                g.mac.mac_spd(g.mac.mac_ivm_idx,max(j-1,1));
        end
        %
        % ivm_sud
        %
        if (g.ivm.n_ivmud ~= 0)
            for idx = 1:numel(xivm_dc.s)
                xivm_dc.s{idx}(:,j) = ...
                    xivm_dc.s{idx}(:,k) ...
                    + h_sol*(xivm_dc.ds{idx}(:,k) + xivm_dc.ds{idx}(:,j))/2;
            end
        end
    end
    % counter increment
    kt = kt + k_inc(ks);
    ks = ks + 1;
end  % end simulation loop

% calculation of line currents post sim
V1 = g.bus.bus_v(g.bus.bus_int(line(:,1)),:);
V2 = g.bus.bus_v(g.bus.bus_int(line(:,2)),:);
R = line(:,3);
X = line(:,4);
B = line(:,5);
tap = line(:,6);
phi = line(:,7);

[ilf,ilt] = line_cur(V1,V2,R,X,B,tap,phi);     % line currents
[sInjF,sInjT] = line_pq(V1,V2,R,X,B,tap,phi);  % line flows (compl. power injections)

% simulation timing
et = toc;
disp(sprintf('elapsed time = %0.3fs',et));

%% Ending simulation routine
warning('*** Ending simulation routine')

if batch_mode
    flag = 1;  % skip interactive plotting
else
    flag = 0;
end

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

    sel = input('Enter your selection, (0--7)[0]: ');
    if isempty(sel)
        sel = 0;  % default
    end

    if (sel == 1)
        mesh(t,[1:1:g.mac.n_mac],g.mac.mac_ang*180/pi)
        title('machine angles')
        xlabel('time in seconds')
        ylabel('internal generator number')
        zlabel('angle in degrees')
        disp('paused: press any key to continue')
        pause
    elseif (sel == 2)
        lt = length(t);
        mesh(t,[1:1:g.mac.n_mac],g.mac.mac_spd-ones(g.mac.n_mac,lt))
        title('machine speed deviations')
        xlabel('time in seconds')
        ylabel('internal generator number')
        zlabel('speed in pu')
        disp('paused: press any key to continue')
        pause
    elseif (sel == 3)
        plot(t,g.mac.pmech)
        title('turbine power')
        xlabel('time in seconds')
        ylabel('power in pu on generator base')
        disp('paused: press any key to continue')
        pause
    elseif (sel == 4)
        plot(t,g.mac.pelect)
        title('generator electrical power')
        xlabel('time in seconds')
        ylabel('power in pu on system base')
        disp('paused: press any key to continue')
        pause
    elseif (sel == 5)
        plot(t,g.exc.Efd)
        title('Generator field voltages')
        xlabel('time in seconds')
        ylabel('voltage in pu')
        disp('paused: press any key to continue')
        pause
    elseif (sel == 6)
        [n_bus,~] = size(g.bus.bus_v);
        mesh(t,(1:1:n_bus),abs(g.bus.bus_v))
        xlabel('time in seconds')
        ylabel('bus number')
        zlabel('voltage in pu')
        disp('paused: press any key to continue')
        pause
    elseif (sel == 7)
        n_line = length(line(:,1));
        if (n_line < 50)
            [S1,S2] = line_pq(V1,V2,R,X,B,tap,phi);
            plot(t,real(S1));
            xlabel('time in seconds')
            ylabel('power flow per unit')
            disp('paused: press any key to continue')
            pause
        else
            % ask for lines to be plotted
            line_range = input('Enter a range of fewer than 50 lines, 1:10 for example [1:49]: ');
            if isempty(line_range)
                line_range = 1:49;  % default
            end

            V1 = g.bus.bus_v(g.bus.bus_int(line(line_range,1)),:);
            V2 = g.bus.bus_v(g.bus.bus_int(line(line_range,2)),:);
            R = line(line_range,3);
            X = line(line_range,4);
            B = line(line_range,5);
            tap = line(line_range,6);
            phi = line(line_range,7);
            [S1,S2] = line_pq(V1,V2,R,X,B,tap,phi);
            plot(t,real(S1));
            xlabel('time in seconds')
            ylabel('power flow per unit')
            disp('paused: press any key to continue')
            pause
        end
    elseif (sel == 0)
        flag = 1;
    else
        error('s_simu: invalid plotting selection.');
    end
end

% hvdc information
t_dc = t_dc(1:length(t_dc)-10);
g.dc.Vdc = g.dc.Vdc(:,1:length(t_dc));
g.dc.i_dc = g.dc.i_dc(:,1:length(t_dc));
g.dc.alpha = g.dc.alpha(:,1:length(t_dc));
g.dc.gamma = g.dc.gamma(:,1:length(t_dc));
g.dc.i_dcr = g.dc.i_dcr(:,1:length(t_dc));
g.dc.i_dci = g.dc.i_dci(:,1:length(t_dc));
g.dc.v_dcc = g.dc.v_dcc(:,1:length(t_dc));
g.dc.v_conr = g.dc.v_conr(:,1:length(t_dc));
g.dc.v_coni = g.dc.v_coni(:,1:length(t_dc));
g.dc.di_dci = g.dc.di_dci(:,1:length(t_dc));
g.dc.di_dcr = g.dc.di_dcr(:,1:length(t_dc));
g.dc.dv_dcc = g.dc.dv_dcc(:,1:length(t_dc));
g.dc.dv_conr = g.dc.dv_conr(:,1:length(t_dc));
g.dc.dv_coni = g.dc.dv_coni(:,1:length(t_dc));
xdcr_dc = xdcr_dc(:,1:length(t_dc));
xdci_dc = xdci_dc(:,1:length(t_dc));
dxdcr_dc = dxdcr_dc(:,1:length(t_dc));
dxdci_dc = dxdci_dc(:,1:length(t_dc));

%-----------------------------------------------------------------------------%
% tidy workspace

% clear global
clear ans B b_num1 b_num2 bo bof bopf1 bopf2 boprf bus_f bus_idx bus_intf
clear bus_intpf1 bus_intpf2 bus_intprf bus_pf1 bus_pf2 bus_sim conv_num dci_dc
clear dcr_dc dcr_states dfile et ets f f1 f2 f3 f_nearbus f_type flag h h_sol
clear H_sum i i_aci i_acr i_plot inv_par IHT j jj k k_end k_inc k_start k_tot
clear ks kt ktmax lfile line_f line_flw line_par line_pf1 line_pf2 line_sim
clear lswitch lt mac_ang_ivm_old n n_switch o_gfma o_reec pathname phi
clear plot_now rec_num rec_par R sel st_state sv svc_dc sw_count SHT t_switch
clear tcsc_dc tdfile tot_states tswitch vnc VLT V_rgf V_rgpf1 V_rgpf2 V_rgprf
clear V_rncf V_rncpf1 V_rncpf2 V_rncprf Vr1 Vr2 X ydcrmn ydcrmx Y1 Y2 Y3 Y4
clear Y_gf Y_gncf Y_gncpf1 Y_gncpf2 Y_gncprf Y_gpf1 Y_gpf2 Y_gprf Y_ncf Y_ncgf
clear Y_ncgpf1 Y_ncgpf2 Y_ncgprf Y_ncpf1 Y_ncpf2 Y_ncprf z z1 z_dpw z_gfma
clear z_pss z_reec z_tg zdc zdcl ze zig zm

% eof
