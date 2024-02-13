% svm_mgen.m
%
% m.file to generate state variable models
% This m file takes the load flow and dynamic data from a data m file
% and calulates the state space matrices:
%     A matrix in a_mat
%
%     B matrices
%       for a change in Exciter Vref in b_vr
%       for a change in Turbine governor Pref b_pr
%       for a change in generator mechanical torque b_pm
%       for a change in svc reference voltage b_svc
%       for a change in tcsc input b_tcsc
%       for a change in active power modulation b_lmod
%       for a change in reactive load modulation b_rlmod
%       for a change in active power modulation b_pwrmod_p
%       for a change in reactive power modulation b_pwrmod_q
%       for a change in ess active power reference b_ess_p
%       for a change in ess reactive power reference b_ess_q
%       for a change in reec voltage reference b_reec_v
%       for a change in reec reactive power reference b_reec_q
%       for a change in gfma active power reference b_gfma_p
%       for a change in gfma reactive power reference b_gfma_q
%
%     C matrices
%       change in generator speed c_spd
%       change in generator electrical torque c_t on generator base
%       change in generator electrical power c_p on generator base
%       change in bus voltage angles c_ang
%       change in bus voltage magnitude c_v
%       change in from_bus line active power c_pf1
%       change in from_bus line reactive power c_qf1
%       change in to_bus line active power c_pf2
%       change in to_bus line reactive power c_qf2
%       change in from bus current magnitude c_ilmf
%       change in to bus current magnitude c_ilmt
%       change in from bus real current c_ilrf
%       change in to bus real current c_ilrt
%       change in from bus imaginary current c_ilif
%       change in to bus imaginary current c_ilit
%       change in ess active power command c_ess_p
%       change in ess reactive power command c_ess_q
%       change in reec measured terminal voltage c_reec_v
%       change in reec measured reactive power c_reec_q
%       change in gfma measured real power c_gfma_p
%       change in gfma measured reactive power c_gfma_q
%
%     D matrices
%       combination of output and input, e.g.,
%       for power out and p_ref in --- d_ppr
%
% l is the eigenvalue vector of a_mat
% u is the right eigenvector matrix of a_mat
% v is the left eigenvector matrix of a_mat (vu = I)
% p is the unscaled participation matrix
% p_norm is the scaled participation matrix
% the maximum value of each column in p_norm is unity
% all scaled participations less than 0.1 are set to zero
% to find the states associated with the jth eigenvalue
% use sparse(abs(p_norm(:,j)))
% pr gives the participation factors of the eigenvalues
% associated with the rotor angle states of each generator
% use sparse(abs(pr(k,:))) to find the modes associated
% with the rotor angle of the kth generator

%-----------------------------------------------------------------------------%
% Version history
%
% Author:  Ryan Elliott
% Date:    September 2023
% Purpose: added support for ivm, reec, and gfma models
%
% Author:  Ryan Elliott
% Date:    July 2020
% Purpose: revised for center-difference method and consistency with s_simu,
%          and added batch processing mode to bypass user queries
%
% Author:  Dan Trudnowski
% Date:    2015
% Purpose: pwrmod_p and pwrmod_q models added
%
% Author:  Graham Rogers
% Date:    December 1998
% Purpose: tcsc model added
%
% Author:  Graham Rogers
% Date:    July 1998
% Purpose: dpw filter added
%
% Author:  Graham Rogers
% Date:    June 1998
% Purpose: hydraulic turbine/governor added
%
% Author:  Graham Rogers
% Date:    August 1997
% Purpose: induction Generator added
%
% Author:  Graham Rogers
% Date:    August 1997
% Purpose: load modulation and output matrices for line flow added
%
% Author:  Graham Rogers
% Date:    April 1997
% Purpose: HVDC added
%
% Version: 1.0
% Date:    September 1996
% Author:  Joe Chow
%-----------------------------------------------------------------------------%

dypar = get_dypar();                     % get parameter settings
batch_mode = dypar.batch_mode;           % batch processing mode

if batch_mode
    % preserve some data for root locus analysis
    clear global; close all;
else
    clear all; close all; clc;
    dypar = get_dypar();
    batch_mode = dypar.batch_mode;
end

nss = get_nstate();                      % get information about model states
lbnd = 1e-3;                             % time constant lower bound

% set user-defined damping control models to empty
svc_dc = {};                             % svc damping control
tcsc_dc = {};                            % tcsc damping control
dcr_dc = {};                             % hvdc rectifier damping control
dci_dc = {};                             % hvdc inverter damping control

% load input data from m.file
disp('linearized model development by perturbation of the non-linear model')

% input data file
if batch_mode
    [dfile,pathname] = get_path();
else
    [dfile,pathname] = uigetfile('d*.m','Select Data File');
end

% import base case data
if (pathname == 0)
    error('svm_mgen: you must select a valid data file.');
else
    full_dfile = fullfile(pathname, dfile);
    run(full_dfile);
    run('import_var');                   % set up global variables
end

% check for valid dynamic data file
if (isempty(g.mac.mac_con) && isempty(g.ivm.ivm_con))
    error('svm_mgen: the selected file is not a valid data file (mac_con).');
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

% determine the initial power flow solution
if batch_mode
    lfs = 'y';                           % solve load flow?
else
    disp(' ')
    lfpf = inputdlg('Do you want to perform a post-fault load flow, (y/n)[n]: ','s');
    lfpf = lower(lfpf);
    if isempty(lfpf)
        lfpf = 'n';
    end

    if strcmp(lfpf,'y')
        disp('Enter the changes to bus and line required to give the post fault condition.')
        disp('When you have finished, type ``return'' and press enter.')
        keyboard
    end

    % solve for loadflow - loadflow parameter
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

g.exc.n_exc = 0;
g.exc.n_smp = 0;
g.exc.n_smppi = 0;
g.exc.n_dc = 0;
g.exc.n_dc1 = 0;
g.exc.n_dc2 = 0;
g.exc.n_st3 = 0;
g.pss.n_pss = 0;
g.dpw.n_dpw = 0;
g.tg.n_tg = 0;
g.tg.n_tgh = 0;
g.tg.n_tg_tot = 0;
g.gfma.n_gfma = 0;
g.svc.n_svc = 0;
g.tcsc.n_tcsc = 0;
g.lsc.n_lsc = 0;
g.reec.n_reec = 0;
g.ess.n_ess = 0;
g.lmod.n_lmod = 0;
g.rlmod.n_rlmod = 0;
g.pwr.n_pwrmod = 0;

% note: dc_indx() called in dc load flow (lfdcs)
mac_indx();
ivm_indx();

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

%-----------------------------------------------------------------------------%
% remove controls associated with infinite bus generators

if (g.mac.n_ib ~= 0)
    % remove exciters
    if ~isempty(g.exc.exc_con)
        net_idx = zeros(length(g.exc.exc_con(:,1)),1);
        for j = 1:g.mac.n_ib
            net_idx = (net_idx | g.exc.exc_con(:,2) ...
                                 == g.mac.mac_con(g.mac.mac_ib_idx(j),1));
        end

        if (length(net_idx) == 1)
            if (net_idx == 1)
                g.exc.exc_con = [];
            end
        else
            perm = diag(~net_idx);
            perm = perm(~net_idx,:);
            g.exc.exc_con = perm*g.exc.exc_con;
        end
    end

    % remove dpw filters
    if ~isempty(g.dpw.dpw_con)
        net_idx = zeros(length(g.dpw.dpw_con(:,1)),1);
        for j = 1:g.mac.n_ib
            net_idx = (net_idx | g.dpw.dpw_con(:,1) ...
                                 == g.mac.mac_con(g.mac.mac_ib_idx(j),1));
        end

        if (length(net_idx) == 1)
            if (net_idx == 1)
                g.dpw.dpw_con = [];
            end
        else
            perm = diag(~net_idx);
            perm = perm(~net_idx,:);
            g.dpw.dpw_con = perm*g.dpw.dpw_con;
        end
    end

    % remove pss
    if ~isempty(g.pss.pss_con)
        net_idx = zeros(length(g.pss.pss_con(:,1)),1);
        for j = 1:g.mac.n_ib
            net_idx = (net_idx | g.pss.pss_con(:,2) ...
                                 == g.mac.mac_con(g.mac.mac_ib_idx(j),1));
        end

        if (length(net_idx) == 1)
            if (net_idx == 1)
                g.pss.pss_con = [];
            end
        else
            perm = diag(~net_idx);
            perm = perm(~net_idx,:);
            g.pss.pss_con = perm*g.pss.pss_con;
        end
    end

    % remove turbine governors
    if ~isempty(g.tg.tg_con)
        net_idx = zeros(length(g.tg.tg_con(:,1)),1);
        for j = 1:g.mac.n_ib
            net_idx = (net_idx | g.tg.tg_con(:,2) ...
                                 == g.mac.mac_con(g.mac.mac_ib_idx(j),1));
        end

        if (length(net_idx) == 1)
            if (net_idx == 1)
                g.tg.tg_con = [];
            end
        else
            perm = diag(~net_idx);
            perm = perm(~net_idx,:);
            g.tg.tg_con = perm*g.tg.tg_con;
        end
    end

    % remove gfma controllers
    if ~isempty(g.gfma.gfma_con)
        net_idx = zeros(length(g.gfma.gfma_con(:,1)),1);
        for j = 1:g.mac.n_ib
            net_idx = (net_idx | g.gfma.gfma_con(:,2) ...
                                 == g.mac.mac_con(g.mac.mac_ib_idx(j),2));
        end

        if (length(net_idx) == 1)
            if (net_idx == 1)
                g.gfma.gfma_con = [];
            end
        else
            perm = diag(~net_idx);
            perm = perm(~net_idx,:);
            g.gfma.gfma_con = perm*g.gfma.gfma_con;
        end
    end
end

% identify exciters
mac_exc = 0;
if ~isempty(g.exc.exc_con)
    exc_indx();
    mac_exc = g.mac.mac_int(g.exc.exc_con(:,2));
else
    g.exc.n_exc = 0;
end

% identify dpw filters
mac_dpw = 0;
if ~isempty(g.dpw.dpw_con)
    dpwf_indx();
    mac_dpw = g.mac.mac_int(g.dpw.dpw_con(:,1));
else
    g.dpw.n_dpw = 0;
end

% identify pss
mac_pss = 0;
if ~isempty(g.pss.pss_con)
    pss_indx();
    mac_pss = g.mac.mac_int(g.pss.pss_con(:,2));
else
    g.pss.n_pss = 0;
end

% identify turbine governors
mac_tg = 0;
mac_tgh = 0;
if ~isempty(g.tg.tg_con)
    tg_indx();
    mac_tg = g.mac.mac_int(g.tg.tg_con(g.tg.tg_idx,2));
    mac_tgh = g.mac.mac_int(g.tg.tg_con(g.tg.tgh_idx,2));
else
    g.tg.n_tg = 0;
    g.tg.n_tgh = 0;
    g.tg.n_tg_tot = 0;
end

% identify gfma controllers
if ~isempty(g.gfma.gfma_con)
    gfma_indx();
else
    g.gfma.n_gfma = 0;
    g.gfma.ivm_idx = [];
    g.gfma.mac_idx = [];
end

% identify svc
if ~isempty(g.svc.svc_con)
    svc_indx(svc_dc);
else
    g.svc.n_svc = 0;
    g.svc.n_svcud = 0;
end

% identify tcsc
if ~isempty(g.tcsc.tcsc_con)
    tcsc_indx(tcsc_dc);
else
    g.tcsc.n_tcsc = 0;
    g.tcsc.n_tcscud = 0;
end

% identify lsc
if ~isempty(g.lsc.lsc_con)
    lsc_indx(bus);
else
    g.lsc.n_lsc = 0;
    g.lsc.n_sensor = 0;
    g.lsc.lsc_idx = [];
end

% identify reec
if ~isempty(g.reec.reec_con)
    reec_indx();
else
    g.reec.n_reec = 0;
    g.reec.ess_idx = [];
end

% identify ess
if ~isempty(g.ess.ess_con)
    ess_indx();
else
    g.ess.n_ess = 0;
    g.ess.n_essud = 0;
    g.ess.ess_idx = [];
end

% identify lmod
if ~isempty(g.lmod.lmod_con)
    lmod_indx();
else
    g.lmod.n_lmod = 0;
end

% identify rlmod
if ~isempty(g.rlmod.rlmod_con)
    rlmod_indx();
else
    g.rlmod.n_rlmod = 0;
end

% identify pwrmod
if ~isempty(g.pwr.pwrmod_con)
    pwrmod_indx(bus);
else
    g.pwr.n_pwrmod = 0;
end

n_bus = length(bus(:,1));
g.bus.n_bus = n_bus;

%-----------------------------------------------------------------------------%
% variable declarations

k = 4;                     %
n_kfd = k;                 % total number of finite-difference data points
kdc = 10*k + 1;            % configuring indexing for hvdc
k_inc = 1;                 %
k_incdc = 10*k_inc;        %
t_switch = 1/(60*8);       % step size
h = t_switch;              %
h_dc = h;                  %

% create zero matrices to initialize variables
z = zeros(size(g.mac.mac_con,1),k);
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
        error('svm_mgen: user-defined hvdc rectifier damping control not allowed.');
    else
        xdcr_dc = zeros(1,kdc);
        dxdcr_dc = zeros(1,kdc);
        dcrd_sig = zeros(1,k);
    end

    if (g.dc.ndci_ud ~= 0)
        error('svm_mgen: user-defined hvdc inverter damping control not allowed.');
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

% ivm models
xivm_dc.s{1} = zeros(max(g.ivm.n_ivmud,1),k);  % user-defined states

if (g.ivm.n_ivmud ~= 0)
    error('svm_mgen: user-defined ivm controllers not allowed.');
end

% gfma controllers
z_gfma = zeros(max(g.gfma.n_gfma,1),k);
o_gfma = ones(max(g.gfma.n_gfma,1),k);

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
g.gfma.lim_flag = false(max(g.gfma.n_gfma,1),1);

% induction motors
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

% induction generators
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

    if (g.svc.n_svcud ~= 0)
        error('svm_mgen: user-defined svc damping control not allowed.');
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

    xsvc_dc = zeros(1,k);
    dxsvc_dc = zeros(1,k);
    d_sig = zeros(1,k);
end

if (g.tcsc.n_tcsc ~= 0)
    g.tcsc.B_tcsc = zeros(g.tcsc.n_tcsc,k);
    g.tcsc.dB_tcsc = zeros(g.tcsc.n_tcsc,k);
    g.tcsc.tcsc_sig = zeros(g.tcsc.n_tcsc,k);
    g.tcsc.tcsc_dsig = zeros(g.tcsc.n_tcsc,k);

    if (g.tcsc.n_tcscud ~= 0)
        error('svm_mgen: user-defined tcsc damping control not allowed.');
    else
        xtcsc_dc = zeros(1,k);
        dxtcsc_dc = zeros(1,k);
    end
else
    g.tcsc.B_tcsc = zeros(1,k);
    g.tcsc.dB_tcsc = zeros(1,k);
    g.tcsc.tcsc_sig = zeros(1,k);
    g.tcsc.tcsc_dsig = zeros(1,k);

    xtcsc_dc = zeros(1,k);
    dxtcsc_dc = zeros(1,k);
    td_sig = zeros(1,k);
end

% lsc variable declarations
z_lsc = zeros(max(g.lsc.n_lsc,1),k);

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

% reec variable declarations
z_reec = zeros(max(g.reec.n_reec,1),k);
o_reec = ones(max(g.reec.n_reec,1),k);

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

g.reec.iqmin = ones(max(g.reec.n_reec,1),1);
g.reec.iqmax = ones(max(g.reec.n_reec,1),1);

g.reec.pref = o_reec;                         % references
g.reec.qref = o_reec;
g.reec.pfaref = zeros(max(g.reec.n_reec,1),1);

g.reec.vref0 = ones(max(g.reec.n_reec,1),1);
g.reec.vref1 = ones(max(g.reec.n_reec,1),1);

g.reec.vdip = false(max(g.reec.n_reec,1),1);
g.reec.vdip_tick = -1*ones(max(g.reec.n_reec,1),1);
g.reec.vdip_time = zeros(max(g.reec.n_reec,1),1);
g.reec.vdip_icmd = zeros(max(g.reec.n_reec,1),1);

g.reec.vblk = false(max(g.reec.n_reec,1),1);
g.reec.vblk_tick = -1*ones(max(g.reec.n_reec,1),1);
g.reec.vblk_time = zeros(max(g.reec.n_reec,1),1);

% ess variable declarations
z_ess = zeros(max(g.ess.n_ess,1),k);

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

% ess user-defined control variable declarations
xess_dc.s{1} = zeros(max(g.ess.n_essud,1),k);  % user-defined states

if (g.ess.n_essud ~= 0)
    error('svm_mgen: user-defined ess controllers not allowed.');
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

g.sys.mach_ref = zeros(1,k);
g.sys.sys_freq = ones(1,k);

%-----------------------------------------------------------------------------%
% initialization

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

% set line parameters
if ~isempty(g.lmon.lmon_con)
    R = line(g.lmon.lmon_con,3);
    X = line(g.lmon.lmon_con,4);
    B = line(g.lmon.lmon_con,5);
    tap = line(g.lmon.lmon_con,6);
    phi = line(g.lmon.lmon_con,7);
end

% construct reduced Y matrices
[Y_gprf,Y_gncprf,Y_ncgprf,Y_ncprf,V_rgprf,V_rncprf,boprf] = red_ybus(bus,line);
bus_intprf = g.bus.bus_int;  % store the internal bus numbers for the pre_fault system

disp('initializing other models')

g.bus.theta(1:n_bus,1) = bus(:,3)*pi/180;
g.bus.bus_v(1:n_bus,1) = bus(:,2).*exp(1j*g.bus.theta(1:n_bus,1));

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
end

% setting initial power flow across all finite difference points
for kl = 2:n_kfd
    g.bus.bus_v(:,kl) = g.bus.bus_v(:,1);
    g.bus.theta(:,kl) = g.bus.theta(:,1);
end

flag = 0;
g.bus.bus_int = bus_intprf;     % pre-fault system

disp('generators')
mac_sub(0,1,bus,flag);
mac_tra(0,1,bus,flag);
mac_em(0,1,bus,flag);
mac_ivm(0,1,bus,flag);
mac_ib(0,1,bus,flag);

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

% user-defined controls not permitted

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

%-----------------------------------------------------------------------------%
% linearization routine

disp(' ')
disp('performing linearization')

state = zeros(g.mac.n_mac,1);
gen_state = state;
TR_state = state;
TB_state = state;
TA_state = state;
Efd_state = state;
R_f_state = state;
pss1_state = state;
pss2_state = state;
pss3_state = state;
dpw_state = state;
tg_state = state;

% generator controls (e.g., exc, pss, dpw, tg) are intentionally
% not counted in n_device_tot because the ``state'' vector tracks
% the number total number of states per machine/device

n_device_tot = g.mac.n_mac + g.ind.n_mot + g.igen.n_ig ...
               + g.svc.n_svc + g.tcsc.n_tcsc ...
               + g.lmod.n_lmod + g.rlmod.n_rlmod + 2*g.pwr.n_pwrmod ...
               + g.dc.n_dcl + g.ess.n_ess + g.lsc.n_lsc + g.reec.n_reec ...
               + g.gfma.n_gfma;

state = zeros(n_device_tot,1);
max_state = nss.mac_max*g.mac.n_mac + nss.exc_max*g.exc.n_exc ...
            + nss.pss_max*g.pss.n_pss + nss.dpw_max*g.dpw.n_dpw ...
            + nss.tg_max*g.tg.n_tg_tot + nss.mot_max*g.ind.n_mot ...
            + nss.ig_max*g.igen.n_ig + nss.svc_max*g.svc.n_svc ...
            + nss.tcsc_max*g.tcsc.n_tcsc + nss.lmod_max*g.lmod.n_lmod ...
            + nss.rlmod_max*g.rlmod.n_rlmod + 2*nss.pwrmod_max*g.pwr.n_pwrmod ...
            + nss.dcl_max*g.dc.n_dcl + nss.ess_max*g.ess.n_ess ...
            + nss.lsc_max*g.lsc.n_lsc + nss.reec_max*g.reec.n_reec ...
            + nss.gfma_max*g.gfma.n_gfma;

% find total number of states
run('ns_file');
NumStates = sum(state);

% set dimensions for A matrix and permutation matrix
a_mat = zeros(NumStates);
p_mat = sparse(zeros(NumStates,max_state));

% set dimensions for generator C matrices
c_curd = zeros(g.mac.n_mac,NumStates);
c_curq = zeros(g.mac.n_mac,NumStates);
c_spd = zeros(g.mac.n_mac,NumStates);
c_pm = zeros(g.mac.n_mac,NumStates);
c_p = zeros(g.mac.n_mac,NumStates);
c_t = zeros(g.mac.n_mac,NumStates);

% determine p_mat
% converts the vector of length max_states to a column of a_mat or b

run('p_m_file');

%-----------------------------------------------------------------------------%
% set initial states

cur = g.mac.cur_re(1:g.mac.n_mac,1) + 1j*g.mac.cur_im(1:g.mac.n_mac,1);
cur_mag(1:g.mac.n_mac,1) = abs(cur(1:g.mac.n_mac,1)).*g.mac.mac_pot(:,1);

g.mac.telect(:,1) = g.mac.pelect(:,1).*g.mac.mac_pot(:,1) ...
                    + cur_mag(:,1).*cur_mag(:,1).*g.mac.mac_con(:,5);

for kl = 2:n_kfd
    if (g.mac.n_mac ~= 0)
        g.mac.mac_ang(:,kl) = g.mac.mac_ang(:,1);
        g.mac.mac_spd(:,kl) = g.mac.mac_spd(:,1);
        g.mac.edprime(:,kl) = g.mac.edprime(:,1);
        g.mac.eqprime(:,kl) = g.mac.eqprime(:,1);
        g.mac.psikd(:,kl) = g.mac.psikd(:,1);
        g.mac.psikq(:,kl) = g.mac.psikq(:,1);

        g.mac.pmech(:,kl) = g.mac.pmech(:,1);
        g.mac.telect(:,kl) = g.mac.telect(:,1);
        g.mac.pelect(:,kl) = g.mac.pelect(:,1);
        g.mac.qelect(:,kl) = g.mac.qelect(:,1);
        g.mac.pm_sig(:,kl) = g.mac.pm_sig(:,1);
        g.mac.ed(:,kl) = g.mac.ed(:,1);
        g.mac.eq(:,kl) = g.mac.eq(:,1);
        g.mac.curd(:,kl) = g.mac.curd(:,1);
        g.mac.curq(:,kl) = g.mac.curq(:,1);
        g.mac.curdg(:,kl)  = g.mac.curdg(:,1);
        g.mac.curqg(:,kl)  = g.mac.curqg(:,1);
        g.mac.fldcur(:,kl) = g.mac.fldcur(:,1);
        g.mac.vex(:,kl) = g.mac.vex(:,1);
        g.mac.eterm(:,kl) = g.mac.eterm(:,1);
        g.mac.cur_re(:,kl) = g.mac.cur_re(:,1);
        g.mac.cur_im(:,kl) = g.mac.cur_im(:,1);
        g.mac.psi_re(:,kl) = g.mac.psi_re(:,1);
        g.mac.psi_im(:,kl) = g.mac.psi_im(:,1);
    end

    if (g.exc.n_exc ~= 0)
        g.exc.V_TR(:,kl) = g.exc.V_TR(:,1);
        g.exc.V_As(:,kl) = g.exc.V_As(:,1);
        g.exc.V_A(:,kl) = g.exc.V_A(:,1);
        g.exc.V_R(:,kl) = g.exc.V_R(:,1);
        g.exc.Efd(:,kl) = g.exc.Efd(:,1);
        g.exc.R_f(:,kl) = g.exc.R_f(:,1);
    end

    if (g.pss.n_pss ~= 0)
        g.pss.pss1(:,kl) = g.pss.pss1(:,1);
        g.pss.pss2(:,kl) = g.pss.pss2(:,1);
        g.pss.pss3(:,kl) = g.pss.pss3(:,1);
    end

    if (g.dpw.n_dpw ~= 0)
        g.dpw.sdpw1(:,kl) = g.dpw.sdpw1(:,1);
        g.dpw.sdpw2(:,kl) = g.dpw.sdpw2(:,1);
        g.dpw.sdpw3(:,kl) = g.dpw.sdpw3(:,1);
        g.dpw.sdpw4(:,kl) = g.dpw.sdpw4(:,1);
        g.dpw.sdpw5(:,kl) = g.dpw.sdpw5(:,1);
        g.dpw.sdpw6(:,kl) = g.dpw.sdpw6(:,1);
    end

    if (g.tg.n_tg_tot ~= 0)
        g.tg.tg1(:,kl) = g.tg.tg1(:,1);
        g.tg.tg2(:,kl) = g.tg.tg2(:,1);
        g.tg.tg3(:,kl) = g.tg.tg3(:,1);
        g.tg.tg4(:,kl) = g.tg.tg4(:,1);
        g.tg.tg5(:,kl) = g.tg.tg5(:,1);
    end

    if (g.gfma.n_gfma ~= 0)
        g.gfma.gfma1(:,kl) = g.gfma.gfma1(:,1);
        g.gfma.gfma2(:,kl) = g.gfma.gfma2(:,1);
        g.gfma.gfma3(:,kl) = g.gfma.gfma3(:,1);
        g.gfma.gfma4(:,kl) = g.gfma.gfma4(:,1);
        g.gfma.gfma5(:,kl) = g.gfma.gfma5(:,1);
        g.gfma.gfma6(:,kl) = g.gfma.gfma6(:,1);
        g.gfma.gfma7(:,kl) = g.gfma.gfma7(:,1);
        g.gfma.gfma8(:,kl) = g.gfma.gfma8(:,1);
        g.gfma.gfma9(:,kl) = g.gfma.gfma9(:,1);
        g.gfma.gfma10(:,kl) = g.gfma.gfma10(:,1);

        g.gfma.dgfma1(:,kl) = g.gfma.dgfma1(:,1);
        g.gfma.dgfma2(:,kl) = g.gfma.dgfma2(:,1);
        g.gfma.dgfma3(:,kl) = g.gfma.dgfma3(:,1);
        g.gfma.dgfma4(:,kl) = g.gfma.dgfma4(:,1);
        g.gfma.dgfma5(:,kl) = g.gfma.dgfma5(:,1);
        g.gfma.dgfma6(:,kl) = g.gfma.dgfma6(:,1);
        g.gfma.dgfma7(:,kl) = g.gfma.dgfma7(:,1);
        g.gfma.dgfma8(:,kl) = g.gfma.dgfma8(:,1);
        g.gfma.dgfma9(:,kl) = g.gfma.dgfma9(:,1);
        g.gfma.dgfma10(:,kl) = g.gfma.dgfma10(:,1);
    end

    if (g.ind.n_mot ~= 0)
        g.ind.vdp(:,kl) = g.ind.vdp(:,1);
        g.ind.vqp(:,kl) = g.ind.vqp(:,1);
        g.ind.slip(:,kl) = g.ind.slip(:,1);
    end

    if (g.igen.n_ig ~= 0)
        g.igen.vdpig(:,kl) = g.igen.vdpig(:,1);
        g.igen.vqpig(:,kl) = g.igen.vqpig(:,1);
        g.igen.slig(:,kl) = g.igen.slig(:,1);
        g.igen.tmig(:,kl) = g.igen.tmig(:,1);
    end

    if (g.svc.n_svc ~= 0)
        g.svc.B_cv(:,kl) = g.svc.B_cv(:,1);
        g.svc.B_con(:,kl) = g.svc.B_con(:,1);
    end

    if (g.tcsc.n_tcsc ~= 0)
        g.tcsc.B_tcsc(:,kl) = g.tcsc.B_tcsc(:,1);
    end

    if (g.ess.n_ess ~= 0)
        g.ess.ess1(:,kl) = g.ess.ess1(:,1);
        g.ess.ess2(:,kl) = g.ess.ess2(:,1);
        g.ess.ess3(:,kl) = g.ess.ess3(:,1);
        g.ess.ess4(:,kl) = g.ess.ess4(:,1);
        g.ess.ess5(:,kl) = g.ess.ess5(:,1);

        g.ess.dess1(:,kl) = g.ess.dess1(:,1);
        g.ess.dess2(:,kl) = g.ess.dess2(:,1);
        g.ess.dess3(:,kl) = g.ess.dess3(:,1);
        g.ess.dess4(:,kl) = g.ess.dess4(:,1);
        g.ess.dess5(:,kl) = g.ess.dess5(:,1);
    end

    if (g.lsc.n_lsc ~= 0)
        g.lsc.lsc1(:,kl) = g.lsc.lsc1(:,1);
        g.lsc.lsc2(:,kl) = g.lsc.lsc2(:,1);
        g.lsc.lsc3(:,kl) = g.lsc.lsc3(:,1);
        g.lsc.lsc4(:,kl) = g.lsc.lsc4(:,1);
        g.lsc.lsc5(:,kl) = g.lsc.lsc5(:,1);
        g.lsc.lsc6(:,kl) = g.lsc.lsc6(:,1);
        g.lsc.lsc7(:,kl) = g.lsc.lsc7(:,1);
        g.lsc.lsc8(:,kl) = g.lsc.lsc8(:,1);
        g.lsc.lsc9(:,kl) = g.lsc.lsc9(:,1);
        g.lsc.lsc10(:,kl) = g.lsc.lsc10(:,1);
        g.lsc.lsc11(:,kl) = g.lsc.lsc11(:,1);
        g.lsc.lsc12(:,kl) = g.lsc.lsc12(:,1);
        g.lsc.lsc13(:,kl) = g.lsc.lsc13(:,1);
        g.lsc.lsc14(:,kl) = g.lsc.lsc14(:,1);
        g.lsc.lsc15(:,kl) = g.lsc.lsc15(:,1);

        g.lsc.dlsc1(:,kl) = g.lsc.dlsc1(:,1);
        g.lsc.dlsc2(:,kl) = g.lsc.dlsc2(:,1);
        g.lsc.dlsc3(:,kl) = g.lsc.dlsc3(:,1);
        g.lsc.dlsc4(:,kl) = g.lsc.dlsc4(:,1);
        g.lsc.dlsc5(:,kl) = g.lsc.dlsc5(:,1);
        g.lsc.dlsc6(:,kl) = g.lsc.dlsc6(:,1);
        g.lsc.dlsc7(:,kl) = g.lsc.dlsc7(:,1);
        g.lsc.dlsc8(:,kl) = g.lsc.dlsc8(:,1);
        g.lsc.dlsc9(:,kl) = g.lsc.dlsc9(:,1);
        g.lsc.dlsc10(:,kl) = g.lsc.dlsc10(:,1);
        g.lsc.dlsc11(:,kl) = g.lsc.dlsc11(:,1);
        g.lsc.dlsc12(:,kl) = g.lsc.dlsc12(:,1);
        g.lsc.dlsc13(:,kl) = g.lsc.dlsc13(:,1);
        g.lsc.dlsc14(:,kl) = g.lsc.dlsc14(:,1);
        g.lsc.dlsc15(:,kl) = g.lsc.dlsc15(:,1);
    end

    if (g.reec.n_reec ~= 0)
        g.reec.reec1(:,kl) = g.reec.reec1(:,1);
        g.reec.reec2(:,kl) = g.reec.reec2(:,1);
        g.reec.reec3(:,kl) = g.reec.reec3(:,1);
        g.reec.reec4(:,kl) = g.reec.reec4(:,1);
        g.reec.reec5(:,kl) = g.reec.reec5(:,1);
        g.reec.reec6(:,kl) = g.reec.reec6(:,1);
        g.reec.reec7(:,kl) = g.reec.reec7(:,1);
        g.reec.reec8(:,kl) = g.reec.reec8(:,1);
        g.reec.reec9(:,kl) = g.reec.reec9(:,1);
        g.reec.reec10(:,kl) = g.reec.reec10(:,1);

        g.reec.dreec1(:,kl) = g.reec.dreec1(:,1);
        g.reec.dreec2(:,kl) = g.reec.dreec2(:,1);
        g.reec.dreec3(:,kl) = g.reec.dreec3(:,1);
        g.reec.dreec4(:,kl) = g.reec.dreec4(:,1);
        g.reec.dreec5(:,kl) = g.reec.dreec5(:,1);
        g.reec.dreec6(:,kl) = g.reec.dreec6(:,1);
        g.reec.dreec7(:,kl) = g.reec.dreec7(:,1);
        g.reec.dreec8(:,kl) = g.reec.dreec8(:,1);
        g.reec.dreec9(:,kl) = g.reec.dreec9(:,1);
        g.reec.dreec10(:,kl) = g.reec.dreec10(:,1);
    end

    if (g.lmod.n_lmod ~= 0)
        g.lmod.lmod_st(:,kl) = g.lmod.lmod_st(:,1);
    end

    if (g.rlmod.n_rlmod ~= 0)
        g.rlmod.rlmod_st(:,kl) = g.rlmod.rlmod_st(:,1);
    end

    if (g.pwr.n_pwrmod ~= 0)
        g.pwr.pwrmod_p_st(:,kl) = g.pwr.pwrmod_p_st(:,1);
        g.pwr.pwrmod_q_st(:,kl) = g.pwr.pwrmod_q_st(:,1);

        g.pwr.pwrmod_p_sig(:,1) = g.pwr.pwrmod_p_st(:,1);
        g.pwr.pwrmod_q_sig(:,1) = g.pwr.pwrmod_q_st(:,1);
    end

    if (g.dc.n_conv ~= 0)
        g.dc.v_conr(:,kl) = g.dc.v_conr(:,1);
        g.dc.v_coni(:,kl) = g.dc.v_coni(:,1);
        g.dc.i_dcr(:,kl) = g.dc.i_dcr(:,1);
        g.dc.i_dci(:,kl) = g.dc.i_dci(:,1);
        g.dc.v_dcc(:,kl) = g.dc.v_dcc(:,1);

        g.dc.Vdc(:,kl) = g.dc.Vdc(:,1);
        g.dc.i_dc(:,kl) = g.dc.i_dc(:,1);
        g.dc.dc_sig(:,kl) = g.dc.dc_sig(:,1);
        g.dc.cur_ord(:,kl) = g.dc.cur_ord(:,1);
        g.dc.alpha(:,kl) = g.dc.alpha(:,1);
        g.dc.gamma(:,kl) = g.dc.gamma(:,1);
    end

    % specify the auxiliary inputs
    g.exc.exc_sig(:,kl) = g.exc.exc_sig(:,1);
    g.tg.tg_sig(:,kl) = g.tg.tg_sig(:,1);
    g.gfma.gfma_sig(:,kl) = g.gfma.gfma_sig(:,1);
    g.svc.svc_sig(:,kl) = g.svc.svc_sig(:,1);
    g.tcsc.tcsc_sig(:,kl) = g.tcsc.tcsc_sig(:,1);
    g.ess.ess_sig(:,kl) = g.ess.ess_sig(:,1);
    g.reec.reec_sig(:,kl) = g.reec.reec_sig(:,1);
    g.lmod.lmod_sig(:,kl) = g.lmod.lmod_sig(:,1);
    g.rlmod.rlmod_sig(:,kl) = g.rlmod.rlmod_sig(:,1);
    g.pwr.pwrmod_p_sig(:,kl) = g.pwr.pwrmod_p_sig(:,1);
    g.pwr.pwrmod_q_sig(:,kl) = g.pwr.pwrmod_q_sig(:,1);
end

%-----------------------------------------------------------------------------%
% perform perturbation of state variables

not_ib_ivm_idx = g.mac.not_ib_idx(ismember(g.mac.not_ib_idx,g.mac.not_ivm_idx));

run('p_cont');

% setup matrix giving state numbers for generators
mac_state = zeros(sum(state(1:g.mac.n_mac)),3);
for k = 1:g.mac.n_mac
    if (state(k) ~= 0)
        if (k == 1)
            j = 1;
        else
            j = 1 + sum(state(1:k-1));
        end

        jj = sum(state(1:k));
        mac_state(j:jj,1) = (j:jj)';
        mac_state(j:jj,2) = st_name(k,1:state(k))';
        mac_state(j:jj,3) = k*ones(state(k),1);
    end
end

ang_idx = find(mac_state(:,2)==1);
spd_idx = find(mac_state(:,2)==2);

b_pm = zeros(NumStates,g.mac.n_mac);
b_pm(spd_idx,not_ib_ivm_idx) = diag(0.5./g.mac.mac_con(not_ib_ivm_idx,16));

st_vec = reshape(st_name.',numel(st_name),1);
st_vec = st_vec(st_vec ~= 0);

% Form transformation matrix to get rid of zero eigenvalue
% Use generator 1 as reference and check for infinite buses
if ~any(g.mac.ibus_con)
    if batch_mode
        ref_gen = 'n';
    else
        ref_gen = inputdlg('Set gen 1 as reference, (y/n)[n]: ');
        ref_gen = lower(ref_gen)
        if isempty(ref_gen)
            ref_gen = 'n';  % default
        end
    end

    if strcmp(ref_gen,'y')
        p_ang = eye(NumStates);
        p_ang(ang_idx,1) = -ones(length(ang_idx),1);
        p_ang(1,1) = 1;
        p_angi = inv(p_ang);

        % transform state matrix
        a_mat = p_ang*a_mat*p_angi;

        % transform the c matrices
        c_v = c_v*p_angi;
        c_ang = c_ang*p_angi;
        c_spd = c_spd*p_angi;
        c_pm = c_pm*p_angi;
        c_t = c_t*p_angi;
        c_p = c_p*p_angi;

        if ~isempty(g.lmon.lmon_con)
            c_pf1 = c_pf1*p_angi;
            c_qf1 = c_qf1*p_angi;
            c_pf2 = c_pf2*p_angi;
            c_qf2 = c_qf2*p_angi;
            c_ilmf = c_ilmf*p_angi;
            c_ilmt = c_ilmt*p_angi;
            c_ilrf = c_ilrf*p_angi;
            c_ilrt = c_ilrt*p_angi;
            c_ilif = c_ilif*p_angi;
            c_ilit = c_ilit*p_angi;
        end

        if (g.dc.n_conv ~= 0)
            c_dcir = c_dcir*p_angi;
            c_dcii = c_dcii*p_angi;
            c_Vdcr = c_dcVr*p_angi;
            c_Vdci = c_Vdci*p_angi;
        end

        % transform the b matrices
        if (g.mac.n_mac ~= 0)
            b_pm = p_ang*b_pm;
        end

        if (g.tg.n_tg ~= 0)
            b_pr = p_ang*b_pr;
        end

        if (g.exc.n_exc ~= 0)
            b_vr = p_ang*b_vr;
        end

        if (g.gfma.n_gfma ~= 0)
            b_gfma_p = p_ang*b_gfma_p;
            b_gfma_q = p_ang*b_gfma_q;
        end

        if (g.svc.n_svc ~= 0)
            b_svc = p_ang*b_svc;
        end

        if (g.tcsc.n_tcsc ~= 0)
            b_tcsc = p_ang*b_tcsc;
        end

        if (g.ess.n_ess ~= 0)
            b_ess = p_ang*b_ess;
        end

        if (g.reec.n_reec ~= 0)
            b_reec_v = p_ang*b_reec_v;
            b_reec_q = p_ang*b_reec_q;
        end

        if (g.lmod.n_lmod ~= 0)
            b_lmod = p_ang*b_lmod;
        end

        if (g.rlmod.n_rlmod ~= 0)
            b_rlmod = p_ang*b_rlmod;
        end

        if (g.pwr.n_pwrmod ~= 0)
            b_pwrmod_p = p_ang*b_pwrmod_p;
            b_pwrmod_q = p_ang*b_pwrmod_q;
        end

        if (g.dc.n_dcl ~= 0)
            b_dcr = p_ang*b_dcr;
            b_dci = p_ang*b_dci;
        end
    end
end

disp('calculating eigenvalues and eigenvectors')
% eigenvectors and eigenvalues of a_mat
[u,l] = eig(a_mat);  % u is the right eigenvector

% sort the eigenvalues
[l,l_idx] = sort(diag(l));

% reorder the eigenvector matrix
u = u(:,l_idx);

for j = 1:NumStates
    if (imag(l(j)) ~= 0)
        % scale the complex eigenvectors so that the maximum element is 1+j0
        [maxu,mu_idx] = max(abs(u(:,j)));
        u(:,j) = u(:,j)/u(mu_idx,j);
    end
end

v = inv(u);  % left eigenvectors

% find the participation factors
p = zeros(NumStates);
p_norm = p;
for j = 1:NumStates
    p(:,j) = (conj(v(j,:)))'.*u(:,j);     % unnormalized participation vectors
    [p_max,p_max_idx] = max((p(:,j)));    %
    p_norm(:,j) = p(:,j)/p(p_max_idx,j);  % p_norm has biggest element = 1
    p_big = abs(p_norm(:,j)) > 0.1;       % big sorts out normalized participation > 1
    p_norm(:,j) = p_big.*p_norm(:,j);     % p_norm now contains only values of p_norm > 0.1
end

% find states associated with the generator angles
pr = p_norm(ang_idx,:);

% frequency and damping ratio
freq = abs(imag(l))/2/pi;

zero_eig = find(abs(l)<=1e-4);
if ~isempty(zero_eig)
    damp(zero_eig,1) = ones(length(zero_eig),1);
end

nz_eig = find(abs(l)>1e-4);
damp(nz_eig,1) = -real(l(nz_eig))./abs(l(nz_eig));

if ~batch_mode
    figure;
    hold;
    box on;
    stab_idx = find(damp>=0.05);
    plot(damp(stab_idx),freq(stab_idx),'k+')
    fmax = ceil(max(freq));
    plot([0.05 0.05],[0 fmax],'r')
    title(['Calculated Modes ' dfile])
    xlabel('damping ratio')
    ylabel('frequency Hz')
    us_idx = find(damp<0);
    plot(damp(us_idx),freq(us_idx),'r+')
    ud_idx = find(damp>0&damp<0.05);
    plot(damp(ud_idx),freq(ud_idx),'g+')
end

%-----------------------------------------------------------------------------%
% tidy work space

% clear global
clear ans B b_v boprf bus_intprf bvnc c_state chk_dc chk_smp cur cur_mag
clear d_vector dc_start dci_dc dcmod_input dfile dpw_count dpw_state
clear Efd_state exc_count exc_number exc_state f flag from_idx gen_state
clear gh IHT i_aci i_acr int_volt inv_par j j_state jj k kj kl k_cex
clear k_col k_colg k_ctg k_dc k_exc k_exc_idx k_hvdc k_idx k_lmod
clear k_nib_idx k_rlmod k_row k_smp k_sub k_tg k_tgh k_tra kgs l_idx
clear l_if1 l_if2 l_it1 l_it2 lf lfile line_flw line_par lmod_input
clear mac_dpw mac_exc mac_pss mac_tg mac_tgh max_state maxu mot_state
clear mu_idx n_hvdc1 n_hvdc_states n_ig_states n_mot ngt no_mac nominal
clear not_TA not_TB not_TE not_TF not_TR not_ib not_ib_ivm_idx ntdc ntf
clear ntl ntpwr_p ntpwr_q ntrl nts nz_eig Pi Pr p_ang p_angi p_big p_mat
clear p_max p_max_idx p_ratio pathname pert phi pr pr_input psi
clear pss1_state pss2_state pss3_state pss_count pss_state Qi Qr R
clear R_f_state rec_par ref_gen rlmod_input SHT s11 s12 s21 s22 s_TA
clear s_TB s_TE s_TR s_idx sel st_name state state_hvdc state_lmod
clear state_rlmod TA_state TB_state TR_state tg_count tg_number tg_state
clear to_idx V0 V1 V2 VLT V_rgprf V_rncprf vnc vr_input X Y_gncprf
clear Y_gprf Y_ncgprf Y_ncprf zero_eig

% eof
