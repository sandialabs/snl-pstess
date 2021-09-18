function h_sol = i_simu(k,ks,k_inc,h,bus_sim,Y_g,Y_gnc,Y_ncg,Y_nc,rec_V1,rec_V2,bo)
% Syntax: h_sol = i_simu(k,ks,k_inc,h,bus_sim,...
%                 Y_g,Y_gnc,Y_ncg,Y_nc,rec_V1,rec_V2,bo)
%
% Purpose: forms the network interface variables
%
% Input:   k - the current time step
%          ks - indicates the switching times
%          k_inc - the number of time seps between switching points
%          h vector of time steps
%          bus_sim value of bus matrix at this switching time
%          Y_g - reduced Y matrix for generators
%          Y_gnc - mutual reduced Y generators-nc loads
%          Y_ncg - mutual reduced Y nc loads generators
%          Y_nc - reduced Y matrix nc loads
%          rec_V1 - voltage recovery matrix generators
%          rec_V2 - voltage recovery matrix nc loads
%          bo bus order for this switching time
%
% Output: h_sol - the time step at this value of ks
%
% Called by: s_simu

%-----------------------------------------------------------------------------%
% Version history
%
% Version: 1.1
% Date:    August 1997
% Purpose: add induction generator
%
% Version: 1.0
% Date:    March 1997
% Author:  Graham Rogers
%-----------------------------------------------------------------------------%

global g;  % declaring struct of global variables

flag = 1;
kdc = 10*(k-1) + 1;

if isempty(g.ind.n_mot)
    g.ind.n_mot = 0;
end

if isempty(g.igen.n_ig)
    g.igen.n_ig = 0;
end

psi = g.mac.psi_re(:,k) + 1j*g.mac.psi_im(:,k);
vmp = g.ind.vdp(:,k) + 1j*g.ind.vqp(:,k);
vmpig = g.igen.vdpig(:,k) + 1j*g.igen.vqpig(:,k);

% specifying the vector of internal voltages
if (g.ind.n_mot ~= 0 & g.igen.n_ig == 0)
    n_tot = g.mac.n_mac + g.ind.n_mot;
    n_gm = g.mac.n_mac + g.ind.n_mot;
    int_volt = [psi; vmp];                  % generators and motors
elseif (g.ind.n_mot == 0 & g.igen.n_ig ~= 0)
    n_tot = g.mac.n_mac + g.igen.n_ig;
    n_gm = g.mac.n_mac;
    int_volt = [psi; vmpig];                % generators and induction generators
elseif  (g.ind.n_mot ~= 0 & g.igen.n_ig ~= 0)
    n_tot = g.mac.n_mac + g.ind.n_mot + g.igen.n_ig;
    n_gm = g.mac.n_mac + g.ind.n_mot;
    int_volt = [psi; vmp; vmpig];           % generators, motors, and ind. gens
else
    int_volt = psi;                         % synchronous generators only
end

h_sol = h(ks);
n_bus = length(bus_sim(:,1));
cur = Y_g*int_volt;                         % network soln currents into machines

% bus voltage reconstruction
b_v(bo(g.ncl.n_load+1:n_bus),1) = rec_V1*int_volt;

if (g.ncl.n_load ~= 0)
    if (k ~= 1)
        kk = k - 1;
    else
        kk = k;
    end

    % initializing nc load voltages
    vnc = g.bus.bus_v(g.bus.bus_int(g.ncl.load_con(:,1)),kk);
    vnc = nc_load(bus_sim,flag,Y_nc,Y_ncg,int_volt,vnc,1e-6,k,kdc);

    % set nc load voltages
    b_v(bo(1:g.ncl.n_load),1) = vnc;
    b_v(bo(g.ncl.n_load+1:n_bus),1) = b_v(bo(g.ncl.n_load+1:n_bus),1) + rec_V2*vnc;

    % modify generator currents for nc loads
    cur = cur + Y_gnc*vnc;
end

% note: the dc bus voltages are the equivalent HT bus voltages
%       and not the LT bus voltages
g.bus.bus_v(g.bus.bus_int(bus_sim(:,1)),k) = b_v;
g.bus.theta(g.bus.bus_int(bus_sim(:,1)),k) = angle(b_v);

g.mac.cur_re(:,k) = real(cur(1:g.mac.n_mac));
g.mac.cur_im(:,k) = imag(cur(1:g.mac.n_mac));             % generator currents

if (g.ind.n_mot ~= 0)
    g.ind.idmot(:,k) = -real(cur(g.mac.n_mac+1:n_gm));    % induction motor currents
    g.ind.iqmot(:,k) = -imag(cur(g.mac.n_mac+1:n_gm));    % current out of network

    g.ind.s_mot(:,k) = g.bus.bus_v(g.bus.bus_int(g.ind.ind_con(:,2)),k) ...
                       .*(g.ind.idmot(:,k)-1j*g.ind.iqmot(:,k));

    g.ind.p_mot(:,k) = real(g.ind.s_mot(:,k));
    g.ind.q_mot(:,k) = imag(g.ind.s_mot(:,k));
end

if (g.igen.n_ig ~= 0)
    g.igen.idig(:,k) = -real(cur(n_gm+1:n_tot));          % induction gen currents
    g.igen.iqig(:,k) = -imag(cur(n_gm+1:n_tot));          % current out of network

    g.igen.s_igen(:,k) = g.bus.bus_v(g.bus.bus_int(g.igen.igen_con(:,2)),k) ...
                         .*(g.igen.idig(:,k)-1j*g.igen.iqig(:,k));

    g.igen.pig(:,k) = real(g.igen.s_igen(:,k));
    g.igen.qig(:,k) = imag(g.igen.s_igen(:,k));
end

if (g.dc.n_conv ~= 0)
    % calculate dc voltage and current
    V0(g.dc.r_idx,1) = abs(g.bus.bus_v(g.dc.rec_ac_bus,k)).*g.dc.dcc_pot(:,7);
    V0(g.dc.i_idx,1) = abs(g.bus.bus_v(g.dc.inv_ac_bus,k)).*g.dc.dcc_pot(:,8);

    g.dc.Vdc(g.dc.r_idx,kdc) = V0(g.dc.r_idx,1).*cos(g.dc.alpha(:,kdc)) ...
                               - g.dc.i_dcr(:,kdc).*g.dc.dcc_pot(:,3);
    g.dc.Vdc(g.dc.i_idx,kdc) = V0(g.dc.i_idx,1).*cos(g.dc.gamma(:,kdc)) ...
                               - g.dc.i_dci(:,kdc).*g.dc.dcc_pot(:,5);

    g.dc.i_dc(g.dc.r_idx,kdc) = g.dc.i_dcr(:,kdc);
    g.dc.i_dc(g.dc.i_idx,kdc) = g.dc.i_dci(:,kdc);
end

end  % function end

% eof
