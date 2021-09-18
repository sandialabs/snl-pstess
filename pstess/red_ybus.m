function [Y11,Y12,Y21,Y22,rec_V1,rec_V2,bus_order] = red_ybus(bus_sol,line)
% Syntax: [red_Y,rec_V] = red_ybus(bus_sol,line)
%         [Y11,Y12,Y21,Y22,rec_V1,rec_V2,bus_order] = red_ybus(bus_sol,line)
%
% Purpose: form the reduced admittance matrix
%
% Input:   bus_sol   - bus solution (generated by loadflow)
%          line      - line data
%
% Output:  red_Y     - reduced admittance matrix
%          rec_V     - voltage reconstruction matrix
%          Y11,Y12,  - reduced admittance matrix for
%          Y21,Y22     systems with non-conforming loads
%          rec_V1,   - voltage reconstruction bus
%          rec_V2
%          bus_order - vector of bus number for recovering bus voltages
%
% See also: loadflow, ybus
%
% Calls: y_sparse
%
% Called By: s_simu, svm_mgen

%-----------------------------------------------------------------------------%
% Version history
%
% Version:  2.35
% Date:     August 2019
% Author:   R. Elliott
% Purpose:  Add energy storage
%
% Version:  2.34
% Date:     Feb 2015
% Author:   D. Trudnowski
% Purpose:  Add power modulation
%
% Version:  2.3
% Modified: Add induction generators, correction to dc
% Date:     August 1997
% Author:   Graham Rogers
%
% Version:  2.2
% Modified: add dc model
% Date:     March 1997
% Author:   Graham Rogers
%
% Version:  2.1
% Modified: add capability to have more than one generator or induction
%           motor on a single bus modification of non-conforming load
%           section to allow call with a single driver
% Date:     November 1996
% Author:   Graham Rogers
%
% Version   2.0
% Modified: remove loops, add induction motors
% Date:     July 1995
% Author:   Graham Rogers
%
% Version:  1.0
% Author:   Kwok W. Cheung, Joe H. Chow
% Date:     March 1991
%-----------------------------------------------------------------------------%

global g;  % declaring struct of global variables

swing_bus = 1;
gen_bus = 2;
load_bus = 3;

n_line = length(line(:,1));                       % number of lines
n_bus = length(bus_sol(:,1));                     % number of buses

n_tot = g.mac.n_mac + g.ind.n_mot + g.igen.n_ig;  % total number of machines
n_gm = g.mac.n_mac + g.ind.n_mot;                 % generators + induction motors
n = g.mac.n_mac;                                  % number of generators

xd = zeros(n,1);

% build sparse admittance matrix Y
Y_d = y_sparse(bus_sol,line);                     % bus admittance matrix
V = bus_sol(:,2);                                 % magnitude of terminal voltage

% compute constant impedance component of non-conforming loads
if (nargout > 2)  % checking number of output arguments
    % non-conforming load ajustments
    % subtract non-conforming loads from bus P and Q loads

    if (g.ncl.n_load ~= 0)
        j = g.bus.bus_int(g.ncl.load_con(:,1));

        bus_sol(j,6) = (ones(g.ncl.n_load,1) - g.ncl.load_con(:,2) ...
                        - g.ncl.load_con(:,4)).*bus_sol(j,6);
        bus_sol(j,7) = (ones(g.ncl.n_load,1) - g.ncl.load_con(:,3) ...
                        - g.ncl.load_con(:,5)).*bus_sol(j,7);

        % uncomment if ess injections are specified as generation
        % if (g.ess.n_ess ~= 0)
        %     j_ess = g.bus.bus_int(g.ess.ess_con(:,2));
        %     bus_sol(j_ess,6) = bus_sol(j_ess,4);
        %     bus_sol(j_ess,7) = bus_sol(j_ess,5);
        % end

        if (g.dc.n_conv ~= 0)
            % remove dc loads from LT bus
            bus_sol(g.dc.ac_bus,6) = zeros(g.dc.n_conv,1);
            bus_sol(g.dc.ac_bus,7) = zeros(g.dc.n_conv,1);
        end
    end
end

% adjust load for pwrmod buses as these are PV buses
if (g.pwr.n_pwrmod ~= 0)
    j_pwrmod = g.bus.bus_int(g.pwr.pwrmod_con(:,1));
    bus_sol(j_pwrmod,6) = bus_sol(j_pwrmod,4);
    bus_sol(j_pwrmod,7) = bus_sol(j_pwrmod,5);
end

% add load components to Y matrix
Pl = bus_sol(:,6);  % real power of loads
Ql = bus_sol(:,7);  % reactive power of loads

% modify load component to take into account generation
gen_exist = zeros(max(bus_sol(:,1)),1);  % buses with no generator data
gen_exist(round(g.mac.mac_con(1:n,2))) = 1:n;

% convert generation with no associated dynamic data to negative load
netgen = find(gen_exist(round(bus_sol(:,1))) < 1);
Pl(netgen) = Pl(netgen) - bus_sol(netgen,4);
Ql(netgen) = Ql(netgen) - bus_sol(netgen,5);

% form constant impedance load admittance for all buses
yl = (Pl - 1j*Ql)./V.^2;
ii = [1:1:n_bus]';
y1 = sparse(ii,ii,yl,n_bus,n_bus);
Y_d = Y_d + y1;  % add to system y matrix

% initialize matrix for machine internal admittances
iin = [1:1:n_tot]';
Y_b = sparse(1,1,0,n_tot,n_bus);

% extract appropriate xdprime and xdpprime from machine data
ra = g.mac.mac_con(:,5)*g.sys.basmva./g.mac.mac_con(:,3);
testxpp = (g.mac.mac_con(:,8) ~= zeros(n,1));
testxp = ~testxpp;

txpp = find(testxpp);
if ~isempty(txpp)
    xd(txpp,1) = g.mac.mac_con(txpp,8)*g.sys.basmva./g.mac.mac_con(txpp,3);  % xppd
end

txp = find(testxp);
if ~isempty(txp)
    xd(txp,1) = g.mac.mac_con(txp,7)*g.sys.basmva./g.mac.mac_con(txp,3);     % xpd
end

y(1:n,1) = ones(n,1)./(ra+1j*xd);

% check for multiple generators at a bus
perm = eye(n);
jg = g.bus.bus_int(round(g.mac.mac_con(:,2)));  % generator buses
for k = 1:n
    mg_idx = find(jg == jg(k));
    lmg = length(mg_idx);
    if (lmg > 1)
        % set 2nd or higher occurences to zero
        jg(mg_idx(2:lmg)) = zeros(lmg-1,1);
        perm(k,mg_idx) = ones(1,lmg);
    end
end

% remove zero elements from jg
jgz_idx = find(jg == 0);
jg(jgz_idx) = [];

% remove zero rows from permutaion matrix
perm(jgz_idx,:) = [];

Ymod = (diag(y))*perm';
Y_b(1:n,jg) = -Ymod;
Y_d = full(Y_d);
Y_d(jg,jg) = Y_d(jg,jg) + perm*Ymod;
Y_d = sparse(Y_d);

% extract appropriate xsp from induction motor data
motmax = 0;
if (length(g.ind.ind_con) ~= 0)
    xsp = g.ind.ind_pot(:,5).*g.ind.ind_pot(:,1);
    rs = g.ind.ind_con(:,4).*g.ind.ind_pot(:,1);
    y(n+1:n_gm,1) = ones(g.ind.n_mot,1)./(rs+1j*xsp);

    % check for multiple induction motors at a bus
    perm = eye(g.ind.n_mot);
    jm = g.bus.bus_int(round(g.ind.ind_con(:,2)));  % motor buses
    for k = 1:g.ind.n_mot
        mm_idx = find(jm == jm(k));
        lmm = length(mm_idx);
        if (lmm > 1)
            % set 2nd or higher occurences to zero
            jm(mm_idx(2:lmm)) = zeros(lmm-1,1);
            perm(k,mm_idx) = ones(1,lmm);
        end
    end

    % remove zero elements from jm
    jmz_idx = find(jm == 0);
    jm(jmz_idx) = [];

    % remove zero rows from permutation matrix
    perm(jmz_idx,:) = [];
    Ymmod = diag(y(n+1:n_gm,1))*perm';
    Y_b(n+1:n_gm,jm) = -Ymmod;
    Y_d(jm,jm) = Y_d(jm,jm) + perm*Ymmod;
    motmax = max(g.ind.ind_con(:,1));
end

% extract appropriate xsp from induction generator data
igmax = 0;
if (g.igen.n_ig ~= 0)
    xsp = g.igen.igen_pot(:,5).*g.igen.igen_pot(:,1);
    rs = g.igen.igen_con(:,4).*g.igen.igen_pot(:,1);
    y(n_gm+1:n_tot,1) = ones(g.igen.n_ig,1)./(rs+1j*xsp);

    % check for multiple induction generators at a bus
    perm = eye(g.igen.n_ig);
    jm = g.bus.bus_int(round(g.igen.igen_con(:,2)));  % induction generator buses
    for k = 1:g.igen.n_ig
        mm_idx = find(jm == jm(k));
        lmm = length(mm_idx);
        if (lmm > 1)
            % set 2nd or higher occurences to zero
            jm(mm_idx(2:lmm)) = zeros(lmm-1,1);
            perm(k,mm_idx) = ones(1,lmm);
        end
    end

    % remove zero elements from jm
    jmz_idx = find(jm == 0);
    jm(jmz_idx) = [];

    % remove zero rows from permutaion matrix
    perm(jmz_idx,:) = [];
    Ymmod = diag(y(n_gm+1:n_tot,1))*perm';
    Y_b(n_gm+1:n_tot,jm) = -Ymmod;
    Y_d(jm,jm) = Y_d(jm,jm) + perm*Ymmod;
    igmax = max(g.igen.igen_con(:,1));
end

Y_a = sparse(iin,iin,y,n_tot,n_tot);
Y_c = Y_b.';  % ordinary transpose (i.e., non-conjugate)

% form the reduced admittance matrix
if (nargout <= 2)
    Y12 = -Y_d\Y_c;
    Y11 = full(Y_a + Y_b*Y12);
    Y12 = full(Y12);  % rec_V
else
    if (g.ncl.n_load ~= 0)
        % non-conforming load Y matrix reduction
        % make vector with non-conforming load buses first
        % note: dc buses must be the last entries in load_con
        bus_order = zeros(n_bus,1);
        bus_conf = zeros(n_bus,1);
        bus_order(1:g.ncl.n_load) = g.bus.bus_int(g.ncl.load_con(:,1));

        % constant impedance buses
        bus_conf(bus_order(1:g.ncl.n_load)) = ones(g.ncl.n_load,1);
        bus_order(g.ncl.n_load+1:n_bus) = find(~bus_conf);

        % make permutation matrix
        P = sparse(1,1,0,n_bus,n_bus);
        P(1:n_bus,bus_order) = eye(n_bus);

        % apply permutation matrix to Y matrix
        % this puts the nonconforming buses in the first
        % n_load by n_load block of Y
        Y_b = Y_b*P';
        Y_c = P*Y_c;
        Y_d = P*Y_d*P';
        % Im = Y_a E_m + Y_b Vb
        % 0  = Y_c E_m + Y_d Vb

        % partition Y matrices
        Y_b1 = Y_b(:,1:g.ncl.n_load);
        Y_b2 = Y_b(:,g.ncl.n_load+1:n_bus);
        Y_c1 = Y_c(1:g.ncl.n_load,:);
        Y_c2 = Y_c(g.ncl.n_load+1:n_bus,:);
        Y_d1 = Y_d(1:g.ncl.n_load,:);
        Y_d2 = Y_d(g.ncl.n_load+1:n_bus,:);
        Y_d11 = Y_d1(:,1:g.ncl.n_load);
        Y_d12 = Y_d1(:,g.ncl.n_load+1:n_bus);
        Y_d21 = Y_d2(:,1:g.ncl.n_load);
        Y_d22 = Y_d2(:,g.ncl.n_load+1:n_bus);

        % yinv = inv(Y_d22);
        % rec_V1 = -yinv*Y_c2;
        % rec_V2 = -yinv*Y_d21;
        rec_V2 = -Y_d22\[Y_c2 Y_d21];
        rec_V1 = rec_V2(:,1:n_tot);
        rec_V2 = rec_V2(:,n_tot+1:n_tot+g.ncl.n_load);
        Y11 = full(Y_a + Y_b2*rec_V1);
        Y12 = full(Y_b1 + Y_b2*rec_V2);
        Y21 = full(Y_c1 + Y_d12*rec_V1);
        Y22 = full(Y_d11 + Y_d12*rec_V2);

        if (g.dc.n_conv ~= 0)
            % modify so that the HV dc voltage replaces the LT dc voltage
            x_dc(g.dc.r_idx) = g.dc.dcc_pot(:,2);
            x_dc(g.dc.i_idx) = g.dc.dcc_pot(:,4);

            % the dc LT buses are after the non-conforming load buses
            n_start = g.ncl.n_load - g.dc.n_conv + 1;
            if (n_start > 1)
                n_fin = n_start - 1;  % end of non-dc non-conforming loads
                y33 = Y22(n_start:g.ncl.n_load,n_start:g.ncl.n_load);
                y31 = Y21(n_start:g.ncl.n_load,:);
                y32 = Y22(n_start:g.ncl.n_load,1:n_fin);
                y21 = Y21(1:n_fin,:);
                y22 = Y22(1:n_fin,1:n_fin);
                y23 = Y22(1:n_fin,n_start:g.ncl.n_load);
                y11 = Y11;
                y12 = Y12(:,1:n_fin);
                y13 = Y12(:,n_start:g.ncl.n_load);
                vr1 = rec_V1;
                vr2 = rec_V2(:,1:n_fin);
                vr3 = rec_V2(:,n_start:g.ncl.n_load);

                % make modifications
                kdc = eye(g.dc.n_conv) - 1j*y33*diag(x_dc);
                kdc = inv(kdc);
                y31 = kdc*y31;
                y32 = kdc*y32;
                y33 = kdc*y33;

                y11 = y11 + 1j*y13*diag(x_dc)*y31;
                y12 = y12 + 1j*y13*diag(x_dc)*y32;
                y13 = y13 + 1j*y13*diag(x_dc)*y33;
                y21 = y21 + 1j*y23*diag(x_dc)*y31;
                y22 = y22 + 1j*y23*diag(x_dc)*y32;
                y23 = y23 + 1j*y23*diag(x_dc)*y33;
                vr1 = vr1 + 1j*vr3*diag(x_dc)*y31;
                vr2 = vr2 + 1j*vr3*diag(x_dc)*y32;
                vr3 = vr3 + 1j*vr3*diag(x_dc)*y33;

                Y11 = y11;
                Y12(:,1:n_fin) = y12;
                Y12(:,n_start:g.ncl.n_load) = y13;
                Y21(1:n_fin,:) = y21;
                Y22(1:n_fin,1:n_fin) = y22;
                Y22(1:n_fin,n_start:g.ncl.n_load)= y23;
                Y21(n_start:g.ncl.n_load,:) = y31;
                Y22(n_start:g.ncl.n_load,1:n_fin) = y32;
                Y22(n_start:g.ncl.n_load,n_start:g.ncl.n_load) = y33;
                rec_V1 = vr1;
                rec_V2(:,1:n_fin) = vr2;
                rec_V2(:,n_start:g.ncl.n_load) = vr3;
            else
                % dc buses only
                y22 = Y22;
                y12 = Y12;
                y11 = Y11;
                y21 = Y21;

                vr1 = rec_V1;
                vr2 = rec_V2;

                kdc = eye(g.dc.n_conv) - 1j*y22*diag(x_dc);
                kdc = inv(kdc);

                y21 = kdc*y21;
                y22 = kdc*y22;
                y11 = y11 + 1j*y12*diag(x_dc)*y21;
                y12 = y12 + 1j*y12*diag(x_dc)*y22;

                vr1 = vr1 + 1j*vr2*diag(x_dc)*y21;
                vr2 = vr2 + 1j*vr2*diag(x_dc)*y22;

                Y11 = y11;
                Y12 = y12;
                Y21 = y21;
                Y22 = y22;

                rec_V1 = vr1;
                rec_V2 = vr2;
            end
        end
    else
        Y12 = -Y_d\Y_c;
        Y11 = full(Y_a + Y_b*Y12);
        Y12 = full(Y12);  % rec_V
        rec_V1 = Y12;
        Y12 = [];
        Y21 = [];
        Y22 = [];
        rec_V2 = [];
        bus_order = g.bus.bus_int(bus_sol(:,1));
    end
end

end  % function end

% eof
