function i_ac = dc_cur(V,kdc)
% Syntax: i_ac = dc_cur(V,kdc)
%
% Purpose: calculates the ac current injected into the network
%          as a function of the equivalent HT voltage
% Input:   V the complex equivalent HT voltage
%          kdc is the dc time step indicator
% Output:  i_ac the pu current injection at the HT bus
%
% Called by nc_load

%-----------------------------------------------------------------------------%
% Version history
%
% Version: 1.0
% Date:    March 1997
%-----------------------------------------------------------------------------%

global g;  % declaring struct of global variables

% note: Vdc is a file local variable that is not the same as g.dc.Vdc

V0(g.dc.r_idx,1) = g.dc.dcc_pot(:,7).*abs(V(g.dc.r_idx));
V0(g.dc.i_idx,1) = g.dc.dcc_pot(:,8).*abs(V(g.dc.i_idx));

dc_ang(g.dc.r_idx,1) = g.dc.alpha(:,kdc);
dc_ang(g.dc.i_idx,1) = g.dc.gamma(:,kdc);

Rc(g.dc.r_idx,1) = g.dc.dcc_pot(:,3);
Rc(g.dc.i_idx,1) = g.dc.dcc_pot(:,5);

idc(g.dc.r_idx,1) = g.dc.i_dcr(:,kdc);
idc(g.dc.i_idx,1) = g.dc.i_dci(:,kdc);

Vdc = V0.*cos(dc_ang) - idc.*Rc;
cphi = Vdc./V0;
sphi = sqrt(ones(g.dc.n_conv,1) - cphi.*cphi);

P = Vdc.*idc/g.sys.basmva;
Q = P.*sphi./cphi;
P(g.dc.i_idx) = -P(g.dc.i_idx);

i_ac = (P - 1j*Q)./conj(V);

end  % function end

% eof
